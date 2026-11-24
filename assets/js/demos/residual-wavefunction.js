/* Why factor a wave function into "analytic envelope x learned residual"? A 1D radial toy.
 *
 * Radial Schrodinger equation (hbar = m = 1, l = 0):  -1/2 u'' + V(r) u = E u,  u(0) = u(R) = 0,
 * with a hard short-range core and an attractive long-range tail: V(r) = C / r^6 - D / r^3 (smoothly cut off).
 * We solve it exactly on a grid (finite differences, Sturm bisection + inverse iteration).
 *
 * Envelope (what physics tells us without solving anything):
 *   short range  (WKB under the C/r^6 wall): u ~ exp(-sqrt(2C) / (2 r^2))
 *   long range   (bound state, energy E):    u ~ exp(-kappa r),  kappa = sqrt(-2E)
 * The residual u / envelope is what a network would still need to learn.
 */
(function () {
  "use strict";
  var el = document.getElementById("wf-demo");
  if (!el) return;
  var cTop = el.querySelector("#wf-top");
  var cBot = el.querySelector("#wf-bot");
  var slider = el.querySelector("#wf-d");
  var out = el.querySelector("#wf-d-out");
  var scaleBtns = el.querySelectorAll("[data-scale]");
  var outE = el.querySelector("#wf-E");
  var outRangeU = el.querySelector("#wf-range-u");
  var outRangeH = el.querySelector("#wf-range-h");

  var C = 1, R0 = 0.2, RMAX = 16, N = 1400;
  var scale = "linear";
  var sol = null;

  function solve(D) {
    var h = (RMAX - R0) / (N + 1), e = -0.5 / (h * h);
    var r = new Float64Array(N), d = new Float64Array(N);
    for (var i = 0; i < N; i++) {
      var x = R0 + (i + 1) * h; r[i] = x;
      d[i] = 1 / (h * h) + C / Math.pow(x, 6) - D / Math.pow(x, 3) * (1 - Math.exp(-Math.pow(x / 0.6, 6)));
    }
    function count(E) {
      var c = 0, q = d[0] - E; if (q < 0) c++;
      for (var i = 1; i < N; i++) { q = d[i] - E - e * e / q; if (q < 0) c++; }
      return c;
    }
    var lo = Infinity; for (i = 0; i < N; i++) lo = Math.min(lo, d[i]);
    lo -= 2 * Math.abs(e);
    var hi = 0;
    if (count(hi) < 1) return null;
    for (var k = 0; k < 70; k++) { var m = (lo + hi) / 2; if (count(m) >= 1) hi = m; else lo = m; }
    var E = (lo + hi) / 2, s = E + 1e-9;
    var u = new Float64Array(N).fill(1);
    var a = new Float64Array(N), b = new Float64Array(N), x2 = new Float64Array(N);
    for (var it = 0; it < 4; it++) {
      a[0] = d[0] - s; b[0] = u[0];
      for (i = 1; i < N; i++) { var w = e / a[i - 1]; a[i] = d[i] - s - w * e; b[i] = u[i] - w * b[i - 1]; }
      x2[N - 1] = b[N - 1] / a[N - 1];
      for (i = N - 2; i >= 0; i--) x2[i] = (b[i] - e * x2[i + 1]) / a[i];
      var nrm = 0; for (i = 0; i < N; i++) nrm += x2[i] * x2[i]; nrm = Math.sqrt(nrm);
      for (i = 0; i < N; i++) u[i] = x2[i] / nrm;
    }
    var mx = 0, imx = 0; for (i = 0; i < N; i++) if (Math.abs(u[i]) > mx) { mx = Math.abs(u[i]); imx = i; }
    var sign = u[imx] < 0 ? -1 : 1;
    for (i = 0; i < N; i++) u[i] = sign * u[i] / mx;
    var kappa = Math.sqrt(-2 * E), aS = Math.sqrt(2 * C) / 2;
    var g = new Float64Array(N), res = new Float64Array(N), gmax = 0;
    for (i = 0; i < N; i++) { g[i] = Math.exp(-aS / (r[i] * r[i]) - kappa * r[i]); if (g[i] > gmax) gmax = g[i]; }
    for (i = 0; i < N; i++) { g[i] /= gmax; res[i] = u[i] / g[i]; }
    var rShow = Math.min(RMAX - 2.5, Math.max(4, r[imx] + 7 / kappa));
    // ranges over the shown window, excluding the last bit where the box wall forces u -> 0
    var umin = Infinity, umax = 0, hmin = Infinity, hmax = 0;
    for (i = 0; i < N; i++) {
      if (r[i] > rShow) break;
      var au = Math.abs(u[i]); if (au > 0) { umin = Math.min(umin, au); umax = Math.max(umax, au); }
      hmin = Math.min(hmin, Math.abs(res[i])); hmax = Math.max(hmax, Math.abs(res[i]));
    }
    var hm = 0; for (i = 0; i < N && r[i] <= rShow; i++) hm = Math.max(hm, Math.abs(res[i]));
    for (i = 0; i < N; i++) res[i] /= hm;
    return { E: E, r: r, u: u, g: g, res: res, kappa: kappa, rShow: rShow,
             uOrders: Math.log10(umax / umin), hOrders: Math.log10(hmax / hmin) };
  }

  function axes(ctx, W, H, L, Rr, T, B, xmax, P, xlabel) {
    ctx.strokeStyle = P.line; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(L, T - 4); ctx.lineTo(L, B); ctx.lineTo(W - Rr, B); ctx.stroke();
    ctx.fillStyle = P.muted; ctx.font = "500 10.5px Inter, sans-serif"; ctx.textAlign = "center";
    var step = xmax > 10 ? 2 : 1;
    for (var x = 0; x <= xmax + 1e-9; x += step) {
      var px = L + x / xmax * (W - L - Rr);
      ctx.fillText(x.toFixed(0), px, B + 14);
    }
    if (xlabel) ctx.fillText(xlabel, (L + W - Rr) / 2, H - 4);
  }

  function drawTop() {
    var P = Demo.palette();
    var g = Demo.fit(cTop, cTop.clientWidth < 520 ? 1.5 : 2.6);
    var ctx = g.ctx, W = g.w, H = g.h;
    ctx.fillStyle = P.bg; ctx.fillRect(0, 0, W, H);
    if (!sol) return;
    var L = 44, Rr = 14, T = 30, B = H - 30, xmax = sol.rShow;
    function X(r) { return L + r / xmax * (W - L - Rr); }
    var lg = scale === "log";
    var ymin = lg ? -30 : -0.05, ymax = lg ? 0.5 : 1.12;
    function Y(v) { return B - (v - ymin) / (ymax - ymin) * (B - T); }
    function val(v) { return lg ? Math.log10(Math.max(1e-30, Math.abs(v))) : v; }
    axes(ctx, W, H, L, Rr, T, B, xmax, P, "distance between the molecules, r");
    // y ticks
    ctx.textAlign = "right"; ctx.fillStyle = P.muted;
    (lg ? [0, -10, -20, -30] : [0, 0.5, 1]).forEach(function (v) {
      ctx.fillText(lg ? (v === 0 ? "1" : "10^" + v) : v.toFixed(1), L - 6, Y(v) + 3);
      ctx.strokeStyle = P.lineSoft; ctx.beginPath(); ctx.moveTo(L, Y(v)); ctx.lineTo(W - Rr, Y(v)); ctx.stroke();
    });
    // wave function (filled in linear)
    ctx.beginPath();
    var started = false;
    for (var i = 0; i < sol.r.length && sol.r[i] <= xmax; i++) {
      var y = Y(Math.max(ymin, val(sol.u[i])));
      if (!started) { ctx.moveTo(X(sol.r[i]), y); started = true; } else ctx.lineTo(X(sol.r[i]), y);
    }
    if (!lg) { ctx.lineTo(X(xmax), Y(0)); ctx.lineTo(X(sol.r[0]), Y(0)); ctx.closePath(); ctx.fillStyle = Demo.alpha(P.olive, 0.16); ctx.fill(); }
    ctx.beginPath(); started = false;
    for (i = 0; i < sol.r.length && sol.r[i] <= xmax; i++) {
      var yy = Y(Math.max(ymin, val(sol.u[i])));
      if (!started) { ctx.moveTo(X(sol.r[i]), yy); started = true; } else ctx.lineTo(X(sol.r[i]), yy);
    }
    ctx.strokeStyle = P.olive; ctx.lineWidth = 2.6; ctx.stroke();
    // envelope
    ctx.beginPath(); started = false;
    for (i = 0; i < sol.r.length && sol.r[i] <= xmax; i++) {
      var ye = Y(Math.max(ymin, val(sol.g[i])));
      if (!started) { ctx.moveTo(X(sol.r[i]), ye); started = true; } else ctx.lineTo(X(sol.r[i]), ye);
    }
    ctx.strokeStyle = P.ochre; ctx.lineWidth = 2; ctx.setLineDash([6, 5]); ctx.stroke(); ctx.setLineDash([]);
    ctx.textAlign = "left"; ctx.font = "500 12px Inter, sans-serif"; ctx.fillStyle = P.muted;
    ctx.fillText("the exact wave function and the analytic envelope" + (lg ? " (log scale)" : ""), 12, 18);
  }

  function drawBot() {
    var P = Demo.palette();
    var g = Demo.fit(cBot, cBot.clientWidth < 520 ? 2.2 : 4.2);
    var ctx = g.ctx, W = g.w, H = g.h;
    ctx.fillStyle = P.bg; ctx.fillRect(0, 0, W, H);
    if (!sol) return;
    var L = 44, Rr = 14, T = 28, B = H - 22, xmax = sol.rShow;
    function X(r) { return L + r / xmax * (W - L - Rr); }
    function Y(v) { return B - v / 1.1 * (B - T); }
    axes(ctx, W, H, L, Rr, T, B, xmax, P, null);
    ctx.textAlign = "right"; ctx.fillStyle = P.muted; ctx.font = "500 10.5px Inter, sans-serif";
    [0, 1].forEach(function (v) { ctx.fillText(v.toFixed(0), L - 6, Y(v) + 3); });
    ctx.beginPath();
    for (var i = 0; i < sol.r.length && sol.r[i] <= xmax; i++) {
      if (i === 0) ctx.moveTo(X(sol.r[i]), Y(sol.res[i])); else ctx.lineTo(X(sol.r[i]), Y(sol.res[i]));
    }
    ctx.strokeStyle = P.oliveDeep; ctx.lineWidth = 2.6; ctx.stroke();
    ctx.textAlign = "left"; ctx.font = "500 12px Inter, sans-serif"; ctx.fillStyle = P.muted;
    ctx.fillText("what's left for the network to learn: wave function ÷ envelope", 12, 16);
  }

  function update() {
    var D = parseFloat(slider.value);
    out.textContent = D.toFixed(1);
    sol = solve(D);
    if (sol) {
      outE.textContent = sol.E.toFixed(2);
      outRangeU.innerHTML = "10<sup>" + Math.round(sol.uOrders) + "</sup>";
      outRangeH.innerHTML = "10<sup>" + sol.hOrders.toFixed(1) + "</sup>";
    }
    drawTop(); drawBot();
  }

  slider.addEventListener("input", update);
  scaleBtns.forEach(function (b) {
    b.addEventListener("click", function () {
      scale = b.getAttribute("data-scale");
      scaleBtns.forEach(function (x) { x.setAttribute("aria-pressed", x === b ? "true" : "false"); });
      drawTop();
    });
  });
  update();
  Demo.onResize(function () { drawTop(); drawBot(); });
})();
