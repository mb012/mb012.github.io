/* i.i.d. vs joint (non-i.i.d.) guidance on an ensemble-averaged measurement. 1D toy.
 *
 * A loop can sit in conformation A (x ~ -1.5) or B (x ~ +1.5). The crystal contains both,
 * so the observed density F_o is bimodal. The prior (our stand-in for AlphaFold) strongly
 * prefers B. Each ensemble member x_k produces a blurred density; the *measured* density is
 * the average over the ensemble. We sample with annealed Langevin dynamics on
 *     prior score  -  eta * grad(data loss)
 * where the data loss is either per-member (i.i.d.) or on the ensemble average (joint).
 */
(function () {
  "use strict";
  var el = document.getElementById("ensemble-demo");
  if (!el) return;
  var canvas = el.querySelector("canvas");
  var modeBtns = el.querySelectorAll("[data-mode]");
  var nSlider = el.querySelector("#ens-n");
  var nOut = el.querySelector("#ens-n-out");
  var runBtn = el.querySelector("#ens-run");
  var outCos = el.querySelector("#ens-cos");
  var outA = el.querySelector("#ens-a");

  var XA = -1.15, XB = 1.15, W = 0.85;                  // conformations and density blur
  var OBS_A = 0.5;                                     // fraction of A in the crystal
  var PRIOR = [{ w: 0.93, m: XB, s: 0.3 }, { w: 0.07, m: XA, s: 0.3 }];
  var G = 160, LO = -3.6, HI = 3.6;
  var grid = new Float64Array(G);
  for (var i = 0; i < G; i++) grid[i] = LO + (HI - LO) * i / (G - 1);
  function gauss(x, m, s) { var d = (x - m) / s; return Math.exp(-0.5 * d * d) / (s * 2.5066); }
  var Fo = new Float64Array(G);
  for (i = 0; i < G; i++) Fo[i] = OBS_A * gauss(grid[i], XA, W) + (1 - OBS_A) * gauss(grid[i], XB, W);

  var mode = "joint";
  var SEED = 11;
  var xs = [];                 // current ensemble
  var history = [];            // frames for animation
  var frame = 0, raf = null;

  function priorScore(x, sig) {
    var num = 0, den = 0;
    for (var j = 0; j < PRIOR.length; j++) {
      var v = PRIOR[j].s * PRIOR[j].s + sig * sig;
      var p = PRIOR[j].w * Math.exp(-0.5 * (x - PRIOR[j].m) * (x - PRIOR[j].m) / v) / Math.sqrt(v);
      num += p * (PRIOR[j].m - x) / v; den += p;
    }
    return num / (den + 1e-300);
  }

  function density(arr) {
    var F = new Float64Array(G), n = arr.length;
    for (var k = 0; k < n; k++) for (var i = 0; i < G; i++) F[i] += gauss(grid[i], arr[k], W) / n;
    return F;
  }

  // gradient of 0.5*||Fo - Fc||^2 wrt each member, Fc = average density of `group`
  function dataGrad(arr, joint) {
    var n = arr.length, g = new Float64Array(n), dx = (HI - LO) / (G - 1);
    if (joint) {
      var F = density(arr);
      for (var k = 0; k < n; k++) {
        var s = 0;
        for (var i = 0; i < G; i++) {
          var r = Fo[i] - F[i];
          s += r * gauss(grid[i], arr[k], W) * (grid[i] - arr[k]) / (W * W);
        }
        g[k] = s * dx;          // = -dL/dx (ascent direction on log-likelihood)
      }
    } else {
      for (k = 0; k < n; k++) {
        var s2 = 0;
        for (i = 0; i < G; i++) {
          var Fk = gauss(grid[i], arr[k], W);
          s2 += (Fo[i] - Fk) * Fk * (grid[i] - arr[k]) / (W * W);
        }
        g[k] = s2 * dx;
      }
    }
    return g;
  }

  function cosine(F) {
    var a = 0, b = 0, c = 0;
    for (var i = 0; i < G; i++) { a += F[i] * Fo[i]; b += F[i] * F[i]; c += Fo[i] * Fo[i]; }
    return a / Math.sqrt(b * c);
  }

  function simulate() {
    var n = parseInt(nSlider.value, 10);
    var r = Demo.rng(SEED + n * 7);
    var x = [];
    for (var k = 0; k < n; k++) x.push(2.5 * r.normal());
    history = [x.slice()];
    var L = 60, STEPS = 4, eta = mode === "none" ? 0 : 150;
    for (var l = 0; l < L; l++) {
      var sig = 2.5 * Math.pow(0.03 / 2.5, l / (L - 1));
      for (var s = 0; s < STEPS; s++) {
        var eps = 0.12 * sig * sig + 0.002;
        var g = eta ? dataGrad(x, mode === "joint") : null;
        for (k = 0; k < n; k++) {
          // guidance strength grows as noise shrinks, like guidance on the denoised estimate
          var drift = priorScore(x[k], sig) + (g ? eta * g[k] / (1 + sig * sig) : 0);
          var stepx = Math.max(-0.35, Math.min(0.35, eps * drift));
          x[k] += stepx + Math.sqrt(2 * eps) * r.normal() * 0.9;
        }
      }
      history.push(x.slice());
    }
    xs = x;
  }

  function readouts(arr) {
    var F = density(arr);
    outCos.textContent = cosine(F).toFixed(2);
    var a = 0, b = 0, mid = 0;
    arr.forEach(function (v) { if (v < XA / 2) a++; else if (v > XB / 2) b++; else mid++; });
    outA.textContent = a + " · " + mid + " · " + b;
  }

  function draw(arr) {
    var P = Demo.palette();
    var g = Demo.fit(canvas, canvas.clientWidth < 520 ? 1.25 : 2.1);
    var ctx = g.ctx, Wd = g.w, H = g.h;
    ctx.fillStyle = P.bg; ctx.fillRect(0, 0, Wd, H);
    var padL = 14, padR = 14, top = 34, base = H - 64;
    function X(v) { return padL + (v - LO) / (HI - LO) * (Wd - padL - padR); }
    var maxF = 0; for (var i = 0; i < G; i++) maxF = Math.max(maxF, Fo[i]);
    var F = density(arr);
    for (i = 0; i < G; i++) maxF = Math.max(maxF, F[i] * 0.85);
    function Y(v) { return base - v / maxF * (base - top) * 0.92; }

    // labels for conformations
    ctx.font = "600 12px Inter, sans-serif"; ctx.textAlign = "center";
    ctx.fillStyle = P.muted;
    ctx.fillText("conformation A", X(XA), H - 12);
    ctx.fillText("conformation B (the prior's favourite)", X(XB), H - 12);
    // prior curve
    ctx.beginPath();
    for (i = 0; i < G; i++) {
      var p = 0; PRIOR.forEach(function (c) { p += c.w * gauss(grid[i], c.m, c.s); });
      var yy = base - p / 1.0 * (base - top) * 0.32;
      if (i === 0) ctx.moveTo(X(grid[i]), yy); else ctx.lineTo(X(grid[i]), yy);
    }
    ctx.strokeStyle = Demo.alpha(P.muted, 0.55); ctx.setLineDash([2, 4]); ctx.lineWidth = 1.2; ctx.stroke(); ctx.setLineDash([]);

    // observed density (filled)
    ctx.beginPath(); ctx.moveTo(X(LO), base);
    for (i = 0; i < G; i++) ctx.lineTo(X(grid[i]), Y(Fo[i]));
    ctx.lineTo(X(HI), base); ctx.closePath();
    ctx.fillStyle = Demo.alpha(P.olive, 0.18); ctx.fill();
    ctx.strokeStyle = Demo.alpha(P.olive, 0.7); ctx.lineWidth = 1.5; ctx.stroke();

    // calculated density from the ensemble
    ctx.beginPath();
    for (i = 0; i < G; i++) { if (i === 0) ctx.moveTo(X(grid[i]), Y(F[i])); else ctx.lineTo(X(grid[i]), Y(F[i])); }
    ctx.strokeStyle = P.ochre; ctx.lineWidth = 2.6; ctx.stroke();

    // baseline + ensemble members
    ctx.strokeStyle = P.line; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(padL, base); ctx.lineTo(Wd - padR, base); ctx.stroke();
    var rows = {};
    arr.forEach(function (v) {
      var key = Math.round(v * 6);
      rows[key] = (rows[key] || 0) + 1;
      var cx = X(Math.max(LO, Math.min(HI, v)));
      ctx.beginPath(); ctx.arc(cx, base + 14, 6, 0, 2 * Math.PI);
      ctx.fillStyle = Demo.alpha(P.oliveDeep, 0.55); ctx.fill();
      ctx.strokeStyle = P.paper; ctx.lineWidth = 1.5; ctx.stroke();
    });

    ctx.textAlign = "left"; ctx.font = "500 12px Inter, sans-serif"; ctx.fillStyle = P.muted;
    var title = { none: "unguided: just the prior", iid: "guided one member at a time (i.i.d.)", joint: "guided as an ensemble (non-i.i.d.)" }[mode];
    ctx.fillText(title, 12, 20);
  }

  function stop() { if (raf) cancelAnimationFrame(raf); raf = null; }
  function animate() {
    stop(); frame = 0;
    if (Demo.reducedMotion()) { frame = history.length - 1; draw(history[frame]); readouts(history[frame]); return; }
    var last = 0;
    (function step(now) {
      if (now - last > 30) { last = now; draw(history[frame]); readouts(history[frame]); frame++; }
      if (frame < history.length) raf = requestAnimationFrame(step);
    })(0);
  }
  function rerun() { simulate(); animate(); }

  modeBtns.forEach(function (b) {
    b.addEventListener("click", function () {
      mode = b.getAttribute("data-mode");
      modeBtns.forEach(function (x) { x.setAttribute("aria-pressed", x === b ? "true" : "false"); });
      rerun();
    });
  });
  nSlider.addEventListener("input", function () { nOut.textContent = nSlider.value; rerun(); });
  runBtn.addEventListener("click", function () { SEED += 1; rerun(); });

  simulate();
  frame = history.length - 1;
  draw(xs); readouts(xs);
  var started = false;
  Demo.whenVisible(el, function () { if (!started) { started = true; animate(); } });
  Demo.onResize(function () { draw(history[Math.min(frame, history.length - 1)] || xs); });
})();
