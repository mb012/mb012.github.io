/* Correcting flows with marginal matching: a 2D toy.
 *
 * Target pi_1: six Gaussians on a ring. Source pi_0: standard Gaussian.
 * The *ideal* flow for linear interpolation x_t = (1-t) x0 + t x1 has a closed-form
 * velocity u_t(x) = E[x1 - x0 | x_t = x]. We make an *imperfect* model v_theta by
 * computing the posterior over mixture components with an inflated variance (as if the
 * network had blurred the modes). Then:
 *   forward  = integrate v_theta from Gaussian noise          (what the model generates)
 *   backward = integrate -v_theta from real data samples     (a path that ends exactly on data)
 *   corrected = forward, but at chosen times take one score step toward the backward
 *               marginal (score estimated with a kernel density on backward samples).
 */
(function () {
  "use strict";
  var el = document.getElementById("flows-demo");
  if (!el) return;

  var canvas = el.querySelector("canvas");
  var viewBtns = el.querySelectorAll("[data-view]");
  var tSlider = el.querySelector("#flows-t");
  var tOut = el.querySelector("#flows-t-out");
  var impSlider = el.querySelector("#flows-imp");
  var impOut = el.querySelector("#flows-imp-out");
  var playBtn = el.querySelector("#flows-play");
  var corrBoxes = el.querySelectorAll("[data-corr]");
  var outFwd = el.querySelector("#flows-off-fwd");
  var outCor = el.querySelector("#flows-off-cor");
  var outBwd = el.querySelector("#flows-off-bwd");

  var K = 6, R = 2.2, S = 0.17;            // target modes
  var MU = [];
  for (var k = 0; k < K; k++) MU.push([R * Math.cos(2 * Math.PI * k / K + 0.3), R * Math.sin(2 * Math.PI * k / K + 0.3)]);
  var N = 900;                              // particles
  var F = 100;                              // frames (t = f / F)
  var SUB = 3;                              // RK4 substeps per frame
  var view = "forward";
  var frame = F;
  var traj = {};                            // forward / backward / corrected : Float32Array[F+1] of 2N

  // ---------- imperfect velocity field ----------
  function velocity(x, y, t, temp, out) {
    var v = t * t * S * S + (1 - t) * (1 - t);
    var vb = v * temp;                      // blurred posterior variance (the "imperfection")
    var wsum = 0, ux = 0, uy = 0, maxl = -1e9, l = [];
    for (var k = 0; k < K; k++) {
      var dx = x - t * MU[k][0], dy = y - t * MU[k][1];
      l[k] = -(dx * dx + dy * dy) / (2 * vb);
      if (l[k] > maxl) maxl = l[k];
    }
    for (k = 0; k < K; k++) {
      var w = Math.exp(l[k] - maxl);
      var mx = MU[k][0], my = MU[k][1];
      var ex = x - t * mx, ey = y - t * my;
      // E[x1|x,k] - E[x0|x,k]
      var c1 = t * S * S / v, c0 = (1 - t) / v;
      ux += w * (mx + c1 * ex - c0 * ex);
      uy += w * (my + c1 * ey - c0 * ey);
      wsum += w;
    }
    out[0] = ux / wsum; out[1] = uy / wsum;
  }

  var k1 = [0, 0], k2 = [0, 0], k3 = [0, 0], k4 = [0, 0];
  function rk4(x, y, t, h, temp, out) {
    velocity(x, y, t, temp, k1);
    velocity(x + h / 2 * k1[0], y + h / 2 * k1[1], t + h / 2, temp, k2);
    velocity(x + h / 2 * k2[0], y + h / 2 * k2[1], t + h / 2, temp, k3);
    velocity(x + h * k3[0], y + h * k3[1], t + h, temp, k4);
    out[0] = x + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    out[1] = y + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
  }

  function integrate(start, fromF, toF, temp, frames, correctAt, bwdFrames) {
    var dir = toF > fromF ? 1 : -1;
    var h = dir / (F * SUB);
    var cur = Float32Array.from(start);
    frames[fromF] = Float32Array.from(cur);
    var o = [0, 0];
    for (var f = fromF; f !== toF; f += dir) {
      for (var s = 0; s < SUB; s++) {
        var t = (f + s / SUB * dir) / F;
        for (var i = 0; i < N; i++) {
          rk4(cur[2 * i], cur[2 * i + 1], t, h, temp, o);
          cur[2 * i] = o[0]; cur[2 * i + 1] = o[1];
        }
      }
      var nf = f + dir;
      if (correctAt && correctAt.indexOf(nf) !== -1) correct(cur, bwdFrames[nf], nf / F);
      frames[nf] = Float32Array.from(cur);
    }
    return frames;
  }

  // One score step toward the backward marginal: x <- x + eta * sigma^2 * grad log p_b(x),
  // with p_b a Gaussian KDE over backward samples. sigma^2 * score = (weighted mean) - x.
  function correct(cur, ref, t) {
    var sig = 0.12 + 0.35 * (1 - t);
    var eta = 0.85;
    var inv = 1 / (2 * sig * sig);
    var nx = new Float32Array(cur.length);
    for (var i = 0; i < N; i++) {
      var x = cur[2 * i], y = cur[2 * i + 1];
      var ws = 0, mx = 0, my = 0, best = 1e9;
      for (var j = 0; j < N; j++) {
        var dx = ref[2 * j] - x, dy = ref[2 * j + 1] - y;
        var d = dx * dx + dy * dy; if (d < best) best = d;
      }
      for (j = 0; j < N; j++) {
        var ex = ref[2 * j] - x, ey = ref[2 * j + 1] - y;
        var w = Math.exp(-(ex * ex + ey * ey - best) * inv);
        ws += w; mx += w * ref[2 * j]; my += w * ref[2 * j + 1];
      }
      nx[2 * i] = x + eta * (mx / ws - x);
      nx[2 * i + 1] = y + eta * (my / ws - y);
    }
    cur.set(nx);
  }

  function offTarget(arr) {
    var off = 0, lim = (3.2 * S) * (3.2 * S);
    for (var i = 0; i < N; i++) {
      var ok = false;
      for (var k = 0; k < K; k++) {
        var dx = arr[2 * i] - MU[k][0], dy = arr[2 * i + 1] - MU[k][1];
        if (dx * dx + dy * dy < lim) { ok = true; break; }
      }
      if (!ok) off++;
    }
    return off / N;
  }

  // ---------- simulate ----------
  var noise, data;
  (function sample() {
    var r = Demo.rng(7);
    noise = new Float32Array(2 * N); data = new Float32Array(2 * N);
    for (var i = 0; i < N; i++) {
      noise[2 * i] = r.normal(); noise[2 * i + 1] = r.normal();
      var k = Math.floor(r() * K);
      data[2 * i] = MU[k][0] + S * r.normal(); data[2 * i + 1] = MU[k][1] + S * r.normal();
    }
  })();

  function corrFrames() {
    var out = [];
    corrBoxes.forEach(function (b) { if (b.checked) out.push(Math.round(parseFloat(b.getAttribute("data-corr")) * F)); });
    return out;
  }

  function simulate() {
    var temp = parseFloat(impSlider.value);
    traj.forward = integrate(noise, 0, F, temp, new Array(F + 1));
    traj.backward = integrate(data, F, 0, temp, new Array(F + 1));
    traj.corrected = integrate(noise, 0, F, temp, new Array(F + 1), corrFrames(), traj.backward);
    // Start the *same* model from where the backward path begins (p_0^b): it lands on the data.
    var fromB = integrate(traj.backward[0], 0, F, temp, new Array(F + 1));
    var fwdOff = offTarget(traj.forward[F]);
    var corOff = offTarget(traj.corrected[F]);
    var bwdOff = offTarget(fromB[F]);
    outFwd.textContent = Math.round(fwdOff * 100) + "%";
    outCor.textContent = Math.round(corOff * 100) + "%";
    outBwd.textContent = Math.round(bwdOff * 100) + "%";
  }

  // ---------- draw ----------
  var g;
  function draw() {
    var P = Demo.palette();
    g = Demo.fit(canvas, canvas.clientWidth < 520 ? 1 : 1.45);
    var ctx = g.ctx, W = g.w, H = g.h;
    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = P.bg; ctx.fillRect(0, 0, W, H);
    var sc = Math.min(W, H) / 7.2, cx = W / 2, cy = H / 2;
    function X(x) { return cx + x * sc; } function Y(y) { return cy - y * sc; }

    var t = frame / F;
    // target modes, fading in with t
    ctx.lineWidth = 1;
    for (var k = 0; k < K; k++) {
      ctx.beginPath(); ctx.arc(X(MU[k][0]), Y(MU[k][1]), 3.2 * S * sc, 0, 2 * Math.PI);
      ctx.strokeStyle = Demo.alpha(P.olive, 0.18 + 0.3 * t); ctx.setLineDash([3, 4]); ctx.stroke();
    }
    ctx.setLineDash([]);
    // source 2-sigma circle
    ctx.beginPath(); ctx.arc(X(0), Y(0), 2 * sc, 0, 2 * Math.PI);
    ctx.strokeStyle = Demo.alpha(P.muted, 0.25 * (1 - t) + 0.05); ctx.setLineDash([2, 5]); ctx.stroke(); ctx.setLineDash([]);

    function dots(arr, color, r, a) {
      ctx.fillStyle = Demo.alpha(color, a);
      for (var i = 0; i < N; i++) {
        ctx.beginPath(); ctx.arc(X(arr[2 * i]), Y(arr[2 * i + 1]), r, 0, 2 * Math.PI); ctx.fill();
      }
    }
    var rr = W < 520 ? 1.5 : 1.9;
    if (view === "forward") dots(traj.forward[frame], P.ochre, rr, 0.62);
    if (view === "backward") dots(traj.backward[frame], P.olive, rr, 0.62);
    if (view === "both") { dots(traj.forward[frame], P.ochre, rr, 0.5); dots(traj.backward[frame], P.olive, rr, 0.5); }
    if (view === "corrected") { dots(traj.backward[frame], P.olive, rr, 0.12); dots(traj.corrected[frame], P.oliveDeep, rr, 0.66); }

    // correction markers on a time bar
    ctx.font = "500 12px Inter, sans-serif"; ctx.fillStyle = P.muted;
    ctx.fillText("t = " + t.toFixed(2), 12, 20);
    var label = { forward: "model's path, from noise", backward: "path from the data, run backwards", both: "both paths at the same time t", corrected: "model's path + corrections" }[view];
    ctx.fillText(label, 12, 38);
    if (view === "corrected") {
      corrFrames().forEach(function (cf) {
        if (frame >= cf) { ctx.fillStyle = P.oliveDeep; ctx.fillText("✓ corrected at t = " + (cf / F).toFixed(1), 12, 56 + 16 * corrFrames().indexOf(cf)); }
      });
    }
  }

  // ---------- controls ----------
  viewBtns.forEach(function (b) {
    b.addEventListener("click", function () {
      view = b.getAttribute("data-view");
      viewBtns.forEach(function (x) { x.setAttribute("aria-pressed", x === b ? "true" : "false"); });
      draw();
    });
  });
  tSlider.addEventListener("input", function () {
    frame = Math.round(parseFloat(tSlider.value) * F); tOut.textContent = (frame / F).toFixed(2); stop(); draw();
  });
  var simTimer;
  function resim() { clearTimeout(simTimer); simTimer = setTimeout(function () { simulate(); draw(); }, 60); }
  impSlider.addEventListener("input", function () { impOut.textContent = parseFloat(impSlider.value).toFixed(1) + "×"; resim(); });
  corrBoxes.forEach(function (b) { b.addEventListener("change", resim); });

  var playing = false, raf;
  function stop() { playing = false; playBtn.textContent = "▶ Play"; cancelAnimationFrame(raf); }
  function play() {
    playing = true; playBtn.textContent = "❚❚ Pause";
    if (frame >= F) frame = 0;
    var last = performance.now();
    (function step(now) {
      if (!playing) return;
      if (now - last > 28) { frame = Math.min(F, frame + 1); last = now; tSlider.value = frame / F; tOut.textContent = (frame / F).toFixed(2); draw(); }
      if (frame >= F) { stop(); return; }
      raf = requestAnimationFrame(step);
    })(last);
  }
  playBtn.addEventListener("click", function () { playing ? stop() : play(); });

  simulate();
  draw();
  Demo.onResize(draw);
})();
