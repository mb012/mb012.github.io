/* Seek and You Shall Fold: gradient-free guidance with a genetic algorithm. 2D toy.
 *
 * A "protein" is a chain of beads. Each residue has one latent number z_i. A fixed decoder (the stand-in
 * for the generator) mixes each latent with its neighbours' and maps the result to a bend angle, so every
 * latent gives a valid chain and changing one latent bends the chain locally around that residue, like a
 * generator decoding a perturbed residue together with its surroundings.
 * The score is NMR-style and non-differentiable:  (fraction of violated distance restraints) x
 * (mean violation).  The search starts from a "prior prediction" whose second half is misfolded.
 *
 * Genetic algorithm in the spirit of the paper (smaller, so it runs in the browser): no crossover, keep
 * elites, tournament selection, and perturb a single residue per child by re-noising its latent to level t:
 *   z_a <- sqrt(1 - t^2) z_a + t * eps.
 * t is either annealed (held high to explore, then lowered every generation to refine), or kept fixed.
 */
(function () {
  "use strict";
  var el = document.getElementById("ga-demo");
  if (!el) return;
  var cChain = el.querySelector("#ga-chain");
  var cCurve = el.querySelector("#ga-curve");
  var schedBtns = el.querySelectorAll("[data-sched]");
  var runBtn = el.querySelector("#ga-run");
  var outGen = el.querySelector("#ga-gen");
  var outViol = el.querySelector("#ga-viol");
  var outScore = el.querySelector("#ga-score");

  var N = 16, TMAX = 1.55;
  var TARGET_ANG = [0.25, 0.3, 0.2, 1.35, 1.3, 0.15, 0.1, 0.25, -0.2, -1.25, -1.2, -0.3, -0.2, -0.15];
  function atanh(x) { return 0.5 * Math.log((1 + x) / (1 - x)); }
  // decoder coupling: each bend angle depends on its own latent and, more weakly, on its neighbours'
  var KER = [1, 0.6, 0.25];
  function kw(i, j) { var d = Math.abs(i - j); return d < KER.length ? KER[d] : 0; }
  function mix(z) {
    return z.map(function (_, i) { var u = 0; for (var j = 0; j < z.length; j++) u += kw(i, j) * z[j]; return u; });
  }
  // latents that decode to the given bend angles (solve the small linear system mix(z) = atanh(angle / TMAX))
  function latentsFor(angles) {
    var n = angles.length, M = [];
    for (var i = 0; i < n; i++) {
      var row = [];
      for (var j = 0; j < n; j++) row.push(kw(i, j));
      row.push(atanh(angles[i] / TMAX));
      M.push(row);
    }
    for (var c = 0; c < n; c++) {
      var piv = c;
      for (var r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r;
      var tmp = M[c]; M[c] = M[piv]; M[piv] = tmp;
      var pv = M[c][c];
      for (var k = c; k <= n; k++) M[c][k] /= pv;
      for (var r2 = 0; r2 < n; r2++) {
        if (r2 === c) continue;
        var f = M[r2][c];
        for (var k2 = c; k2 <= n; k2++) M[r2][k2] -= f * M[c][k2];
      }
    }
    return M.map(function (row) { return row[n]; });
  }
  var zStar = latentsFor(TARGET_ANG);
  // prior prediction: right first half, misfolded second half (the tail curls the wrong way)
  var PRIOR_ANG = TARGET_ANG.map(function (a, i) { return i < 8 ? a : TMAX * Math.tanh(0.9); });
  var zPrior = latentsFor(PRIOR_ANG);

  function decode(z) {
    var pts = [[0, 0], [1, 0]], ang = 0, u = mix(z);
    for (var i = 0; i < u.length; i++) {
      ang += TMAX * Math.tanh(u[i]);
      var p = pts[pts.length - 1];
      pts.push([p[0] + Math.cos(ang), p[1] + Math.sin(ang)]);
    }
    return pts;
  }
  function dist(a, b) { var dx = a[0] - b[0], dy = a[1] - b[1]; return Math.sqrt(dx * dx + dy * dy); }

  // restraints: close pairs in the target (|i-j| >= 3); each wants its target distance within +-TOL
  var TOL = 0.2, CLASH = 0.75;
  var target = decode(zStar), R = [];
  for (var i = 0; i < N; i++) for (var j = i + 3; j < N; j++) {
    var d = dist(target[i], target[j]);
    if (d < 2.6) R.push([i, j, d]);
  }
  function viol(dd, ref) { return Math.abs(dd - ref) - TOL; }

  function evaluate(pts) {
    var nv = 0, sv = 0;
    for (var k = 0; k < R.length; k++) {
      var v = viol(dist(pts[R[k][0]], pts[R[k][1]]), R[k][2]);
      if (v > 0) { nv++; sv += v; }
    }
    var restraintViol = nv;
    // chains can't pass through themselves: clashing beads count as violations too
    for (var a = 0; a < N; a++) for (var b = a + 2; b < N; b++) {
      var dd = dist(pts[a], pts[b]);
      if (dd < CLASH) { nv++; sv += 3 * (CLASH - dd); }
    }
    return { nv: restraintViol, score: nv ? (nv / R.length) * (sv / nv) : 0 };
  }

  var sched = "anneal";
  var GENS = 300, POP = 50, ELITE = 3;
  var T_HI = 0.95, T_LO = 0.02, HOLD = 0.2;
  function tAt(k) {
    // annealed: explore at high t for the first 20% of generations, then lower t every generation (geometrically)
    if (sched === "anneal") {
      var f = k / GENS;
      return f < HOLD ? T_HI : T_HI * Math.pow(T_LO / T_HI, (f - HOLD) / (1 - HOLD));
    }
    if (sched === "big") return 0.95;
    return 0.03;
  }

  var runs = { anneal: null, big: null, small: null }; // stored curves for comparison
  var cur = null, anim = null, seed = 8;

  function runGA(s, seedv) {
    var r = Demo.rng(seedv);
    var pop = [];
    for (var p = 0; p < POP; p++) pop.push(zPrior.slice());
    var curve = [], bestHist = [];
    var saved = sched; sched = s;
    for (var k = 0; k <= GENS; k++) {
      var scored = pop.map(function (z) { var e = evaluate(decode(z)); return { z: z, e: e }; });
      scored.sort(function (a, b) { return a.e.score - b.e.score; });
      curve.push(scored[0].e.score);
      bestHist.push({ z: scored[0].z.slice(), e: scored[0].e });
      if (k === GENS) break;
      var next = scored.slice(0, ELITE).map(function (x) { return x.z.slice(); });
      var t = tAt(k);
      while (next.length < POP) {
        // tournament of 3
        var best = null;
        for (var q = 0; q < 3; q++) { var c = scored[Math.floor(r() * POP)]; if (!best || c.e.score < best.e.score) best = c; }
        var child = best.z.slice();
        var a = Math.floor(r() * child.length);
        child[a] = Math.sqrt(1 - t * t) * child[a] + t * r.normal() * 1.3;
        next.push(child);
      }
      pop = next;
    }
    sched = saved;
    return { curve: curve, best: bestHist };
  }

  // ---------- drawing ----------
  function drawChain(z, e, k) {
    var P = Demo.palette();
    var g = Demo.fit(cChain, 1.15);
    var ctx = g.ctx, W = g.w, H = g.h;
    ctx.fillStyle = P.bg; ctx.fillRect(0, 0, W, H);
    var pts = decode(z);
    // frame on target + current
    var all = target.concat(pts), minx = 1e9, maxx = -1e9, miny = 1e9, maxy = -1e9;
    all.forEach(function (p) { minx = Math.min(minx, p[0]); maxx = Math.max(maxx, p[0]); miny = Math.min(miny, p[1]); maxy = Math.max(maxy, p[1]); });
    var sc = Math.min((W - 40) / (maxx - minx + 1), (H - 60) / (maxy - miny + 1));
    var ox = W / 2 - sc * (minx + maxx) / 2, oy = H / 2 + 10 + sc * (miny + maxy) / 2;
    function X(p) { return ox + p[0] * sc; } function Y(p) { return oy - p[1] * sc; }

    // target ghost
    ctx.lineWidth = Math.max(6, sc * 0.32); ctx.lineCap = "round"; ctx.lineJoin = "round";
    ctx.strokeStyle = Demo.alpha(P.muted, 0.14);
    ctx.beginPath(); target.forEach(function (p, i) { if (i) ctx.lineTo(X(p), Y(p)); else ctx.moveTo(X(p), Y(p)); }); ctx.stroke();

    // restraints
    ctx.lineWidth = 1.5;
    R.forEach(function (rr) {
      var a = pts[rr[0]], b = pts[rr[1]], v = viol(dist(a, b), rr[2]);
      ctx.strokeStyle = v > 0 ? Demo.alpha(P.rust, 0.75) : Demo.alpha(P.olive, 0.25);
      ctx.setLineDash(v > 0 ? [] : [2, 3]);
      ctx.beginPath(); ctx.moveTo(X(a), Y(a)); ctx.lineTo(X(b), Y(b)); ctx.stroke();
    });
    ctx.setLineDash([]);

    // chain
    ctx.lineWidth = 3.2; ctx.strokeStyle = P.oliveDeep;
    ctx.beginPath(); pts.forEach(function (p, i) { if (i) ctx.lineTo(X(p), Y(p)); else ctx.moveTo(X(p), Y(p)); }); ctx.stroke();
    pts.forEach(function (p, i) {
      ctx.beginPath(); ctx.arc(X(p), Y(p), i === 0 ? 5.5 : 4.2, 0, 2 * Math.PI);
      ctx.fillStyle = i === 0 ? P.ochre : P.olive; ctx.fill();
      ctx.strokeStyle = P.bg; ctx.lineWidth = 1.5; ctx.stroke();
    });
    ctx.font = "500 12px Inter, sans-serif"; ctx.fillStyle = P.muted; ctx.textAlign = "left";
    ctx.fillText(k === 0 ? "start: the prior's prediction" : "best structure, generation " + k, 12, 20);
  }

  function drawCurve(upto) {
    var P = Demo.palette();
    var g = Demo.fit(cCurve, 1.15);
    var ctx = g.ctx, W = g.w, H = g.h;
    ctx.fillStyle = P.bg; ctx.fillRect(0, 0, W, H);
    var L = 40, Rr = 14, T = 30, B = H - 34;
    var y0 = runs.anneal ? runs.anneal.curve[0] : (cur ? cur.curve[0] : 1);
    var ymax = y0 * 1.05;
    function X(k) { return L + k / GENS * (W - L - Rr); }
    function Y(v) { return B - (v / ymax) * (B - T); }
    // axes
    ctx.strokeStyle = P.line; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(L, T - 6); ctx.lineTo(L, B); ctx.lineTo(W - Rr, B); ctx.stroke();
    ctx.fillStyle = P.muted; ctx.font = "500 11px Inter, sans-serif"; ctx.textAlign = "center";
    ctx.fillText("generation", (L + W - Rr) / 2, H - 10);
    ctx.save(); ctx.translate(13, (T + B) / 2); ctx.rotate(-Math.PI / 2); ctx.fillText("violation score (lower is better)", 0, 0); ctx.restore();
    ctx.textAlign = "left"; ctx.font = "500 12px Inter, sans-serif";
    ctx.fillText("best score so far", 12, 18);

    var colors = { anneal: P.olive, big: P.ochre, small: P.rust };
    var names = { anneal: "annealed t", big: "fixed, large t", small: "fixed, small t" };
    var ly = T + 4;
    Object.keys(runs).forEach(function (s) {
      var run = runs[s]; if (!run) return;
      var n = (cur && s === cur.s) ? upto : GENS;
      ctx.strokeStyle = colors[s]; ctx.lineWidth = s === (cur && cur.s) ? 2.8 : 1.6;
      ctx.globalAlpha = s === (cur && cur.s) ? 1 : 0.55;
      ctx.beginPath();
      for (var k = 0; k <= n; k++) { if (k) ctx.lineTo(X(k), Y(run.curve[k])); else ctx.moveTo(X(k), Y(run.curve[k])); }
      ctx.stroke(); ctx.globalAlpha = 1;
      ctx.fillStyle = colors[s]; ctx.fillRect(W - Rr - 118, ly, 10, 3);
      ctx.fillStyle = P.inkSoft; ctx.fillText(names[s], W - Rr - 102, ly + 5); ly += 16;
    });
  }

  function show(k) {
    var b = cur.best[k];
    drawChain(b.z, b.e, k);
    drawCurve(k);
    outGen.textContent = k;
    outViol.textContent = b.e.nv + " / " + R.length;
    outScore.textContent = b.e.score.toFixed(3);
  }

  function start() {
    if (anim) cancelAnimationFrame(anim);
    cur = runGA(sched, seed); cur.s = sched;
    runs[sched] = cur;
    if (Demo.reducedMotion()) { show(GENS); return; }
    var k = 0, last = 0;
    (function step(now) {
      if (now - last > 22) { last = now; show(Math.min(k, GENS)); k += 2; }
      if (k <= GENS + 1) anim = requestAnimationFrame(step);
    })(0);
  }

  schedBtns.forEach(function (b) {
    b.addEventListener("click", function () {
      sched = b.getAttribute("data-sched");
      schedBtns.forEach(function (x) { x.setAttribute("aria-pressed", x === b ? "true" : "false"); });
      start();
    });
  });
  runBtn.addEventListener("click", function () { seed += 1; runs = { anneal: null, big: null, small: null }; start(); });

  // initial still frame, start animating when visible
  cur = runGA(sched, seed); cur.s = sched; runs[sched] = cur;
  show(0);
  var started = false;
  Demo.whenVisible(el, function () { if (!started) { started = true; start(); } });
  Demo.onResize(function () { if (cur) show(Math.min(GENS, cur.best.length - 1)); });

})();
