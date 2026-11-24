/* Ring currents around a rotating phenylalanine ring: the physics a good shift predictor should respect.
 *
 * Point-dipole (Pople) approximation:  delta(ppm) = B * (1 - 3 cos^2 theta) / r^3
 * r: distance from the ring centre, theta: angle between r and the ring normal.
 * The ring rotates about its CB-CG axis (x axis here) by chi2, so the normal is n = (0, cos chi, sin chi).
 * Because n and -n give the same cos^2, the effect is 180-degree periodic; because of 1/r^3 it is local.
 */
(function () {
  "use strict";
  var el = document.getElementById("ring-demo");
  if (!el) return;
  var cScene = el.querySelector("#ring-scene");
  var cPlot = el.querySelector("#ring-plot");
  var slider = el.querySelector("#ring-chi");
  var out = el.querySelector("#ring-chi-out");
  var spinBtn = el.querySelector("#ring-spin");

  var B = 27.6; // ppm * A^3, order of magnitude for a benzene ring current
  var ATOMS = [
    { name: "CA 50", p: [-2.6, 1.2, 1.0] },
    { name: "CA 51", p: [1.0, 3.2, -2.4] },
    { name: "CA 49", p: [-3.5, -3.0, 2.6] },
    { name: "CA 40", p: [6.2, -5.6, 4.6] }
  ];
  ATOMS.forEach(function (a) { a.r = Math.hypot(a.p[0], a.p[1], a.p[2]); });

  function shift(a, chi) {
    var ny = Math.cos(chi), nz = Math.sin(chi);
    var c = (a.p[1] * ny + a.p[2] * nz) / a.r;
    return B * (1 - 3 * c * c) / (a.r * a.r * a.r);
  }

  var chiDeg = parseFloat(slider.value);

  // simple fixed camera: rotate world a bit for a 3/4 view, orthographic projection
  var yaw = -0.55, pitch = 0.42;
  function proj(p) {
    var x = p[0] * Math.cos(yaw) - p[2] * Math.sin(yaw);
    var z = p[0] * Math.sin(yaw) + p[2] * Math.cos(yaw);
    var y = p[1] * Math.cos(pitch) - z * Math.sin(pitch);
    var d = p[1] * Math.sin(pitch) + z * Math.cos(pitch);
    return [x, y, d];
  }

  function colors(P) { return [P.olive, P.ochre, P.rust, P.muted]; }

  function drawScene() {
    var P = Demo.palette();
    var g = Demo.fit(cScene, 1.15);
    var ctx = g.ctx, W = g.w, H = g.h;
    ctx.fillStyle = P.bg; ctx.fillRect(0, 0, W, H);
    var sc = Math.min(W, H) / 19, cx = W * 0.47, cy = H * 0.52;
    function S(p) { var q = proj(p); return [cx + q[0] * sc, cy - q[1] * sc, q[2]]; }
    var chi = chiDeg * Math.PI / 180;

    // 7 A sphere outline (projected circle)
    ctx.beginPath(); ctx.arc(cx, cy, 7 * sc, 0, 2 * Math.PI);
    ctx.strokeStyle = Demo.alpha(P.muted, 0.35); ctx.setLineDash([3, 5]); ctx.lineWidth = 1; ctx.stroke(); ctx.setLineDash([]);
    ctx.font = "500 11px Inter, sans-serif"; ctx.fillStyle = P.muted; ctx.textAlign = "left";
    ctx.fillText("7 Å", cx + 7 * sc * 0.72 + 4, cy - 7 * sc * 0.72);

    // CB-CG axis
    var a0 = S([-2.6, 0, 0]), a1 = S([0, 0, 0]);
    ctx.strokeStyle = Demo.alpha(P.inkSoft, 0.6); ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(a0[0], a0[1]); ctx.lineTo(a1[0], a1[1]); ctx.stroke();

    // ring: hexagon in plane spanned by x-axis and u = (0, -sin chi, cos chi)
    var u = [0, -Math.sin(chi), Math.cos(chi)], rr = 1.39;
    // ring centred at x = 1.39 so that one vertex (CG) sits on the axis at the origin
    var hex = [];
    for (var k = 0; k < 6; k++) {
      var tt = Math.PI + k * Math.PI / 3;
      hex.push(S([1.39 + rr * Math.cos(tt), rr * Math.sin(tt) * u[1], rr * Math.sin(tt) * u[2]]));
    }
    ctx.beginPath(); hex.forEach(function (h, i) { if (i) ctx.lineTo(h[0], h[1]); else ctx.moveTo(h[0], h[1]); }); ctx.closePath();
    ctx.fillStyle = Demo.alpha(P.ochre, 0.22); ctx.fill();
    ctx.strokeStyle = P.ochre; ctx.lineWidth = 2.4; ctx.stroke();
    // ring normal
    var c0 = S([1.39, 0, 0]), n1 = S([1.39, 1.8 * Math.cos(chi), 1.8 * Math.sin(chi)]), n2 = S([1.39, -1.8 * Math.cos(chi), -1.8 * Math.sin(chi)]);
    ctx.strokeStyle = Demo.alpha(P.ochre, 0.7); ctx.lineWidth = 1.3; ctx.setLineDash([3, 3]);
    ctx.beginPath(); ctx.moveTo(n2[0], n2[1]); ctx.lineTo(n1[0], n1[1]); ctx.stroke(); ctx.setLineDash([]);

    // atoms (positions are relative to the ring centre; offset by ring centre x=1.39)
    var cols = colors(P);
    ATOMS.forEach(function (a, ai) {
      var p = S([a.p[0] + 1.39, a.p[1], a.p[2]]);
      var v = shift(a, chi);
      ctx.strokeStyle = Demo.alpha(P.muted, 0.35); ctx.lineWidth = 1; ctx.setLineDash([2, 4]);
      ctx.beginPath(); ctx.moveTo(c0[0], c0[1]); ctx.lineTo(p[0], p[1]); ctx.stroke(); ctx.setLineDash([]);
      var rad = 5 + Math.min(10, Math.abs(v) * 7);
      ctx.beginPath(); ctx.arc(p[0], p[1], rad, 0, 2 * Math.PI);
      ctx.fillStyle = Demo.alpha(cols[ai], 0.9); ctx.fill();
      ctx.strokeStyle = P.bg; ctx.lineWidth = 1.5; ctx.stroke();
      var left = p[0] > W * 0.62 || a.name === "CA 50";
      var tx = left ? p[0] - rad - 4 : p[0] + rad + 4;
      ctx.textAlign = left ? "right" : "left";
      ctx.fillStyle = P.ink; ctx.font = "600 11px Inter, sans-serif";
      ctx.fillText(a.name, tx, p[1] - 2);
      ctx.fillStyle = P.muted; ctx.font = "500 10.5px Inter, sans-serif";
      ctx.fillText(a.r.toFixed(1) + " Å · " + (v >= 0 ? "+" : "") + v.toFixed(2) + " ppm", tx, p[1] + 11);
    });

    ctx.fillStyle = P.muted; ctx.font = "500 12px Inter, sans-serif"; ctx.textAlign = "left";
    ctx.fillText("phenylalanine ring, χ2 = " + Math.round(chiDeg) + "°", 12, 20);
  }

  function drawPlot() {
    var P = Demo.palette();
    var g = Demo.fit(cPlot, 1.15);
    var ctx = g.ctx, W = g.w, H = g.h;
    ctx.fillStyle = P.bg; ctx.fillRect(0, 0, W, H);
    var L = 46, R = 12, T = 30, Bt = H - 36;
    var ymin = -0.9, ymax = 1.15;
    function X(d) { return L + (d + 180) / 360 * (W - L - R); }
    function Y(v) { return T + (ymax - v) / (ymax - ymin) * (Bt - T); }
    // grid
    ctx.strokeStyle = P.lineSoft; ctx.lineWidth = 1; ctx.font = "500 10.5px Inter, sans-serif"; ctx.fillStyle = P.muted;
    ctx.textAlign = "right";
    [-0.5, 0, 0.5, 1].forEach(function (v) { ctx.beginPath(); ctx.moveTo(L, Y(v)); ctx.lineTo(W - R, Y(v)); ctx.stroke(); ctx.fillText(v.toFixed(1), L - 6, Y(v) + 3); });
    ctx.textAlign = "center";
    [-180, -90, 0, 90, 180].forEach(function (d) { ctx.fillText(d + "°", X(d), Bt + 16); });
    ctx.fillText("ring angle χ2", (L + W - R) / 2, H - 4);
    ctx.save(); ctx.translate(12, (T + Bt) / 2); ctx.rotate(-Math.PI / 2); ctx.fillText("ring-current shift (ppm)", 0, 0); ctx.restore();
    ctx.strokeStyle = P.line; ctx.beginPath(); ctx.moveTo(L, Y(0)); ctx.lineTo(W - R, Y(0)); ctx.stroke();

    var cols = colors(P);
    ATOMS.forEach(function (a, i) {
      ctx.strokeStyle = cols[i]; ctx.lineWidth = i === 3 ? 1.4 : 2.2;
      ctx.beginPath();
      for (var d = -180; d <= 180; d += 3) {
        var v = shift(a, d * Math.PI / 180);
        if (d === -180) ctx.moveTo(X(d), Y(v)); else ctx.lineTo(X(d), Y(v));
      }
      ctx.stroke();
    });
    // marker
    ctx.strokeStyle = P.ochre; ctx.lineWidth = 1.5; ctx.setLineDash([4, 3]);
    ctx.beginPath(); ctx.moveTo(X(chiDeg), T); ctx.lineTo(X(chiDeg), Bt); ctx.stroke(); ctx.setLineDash([]);
    ATOMS.forEach(function (a, i) {
      var v = shift(a, chiDeg * Math.PI / 180);
      ctx.beginPath(); ctx.arc(X(chiDeg), Y(v), 4.5, 0, 2 * Math.PI); ctx.fillStyle = cols[i]; ctx.fill();
    });
    // legend
    ctx.textAlign = "left"; ctx.font = "500 11px Inter, sans-serif";
    ATOMS.forEach(function (a, i) {
      var lx = L + 8 + i * Math.min(92, (W - L - R) / 4), ly = T - 14;
      ctx.fillStyle = cols[i]; ctx.fillRect(lx, ly - 4, 12, 3);
      ctx.fillStyle = P.inkSoft; ctx.fillText(a.name, lx + 16, ly);
    });
  }

  function draw() { drawScene(); drawPlot(); out.textContent = Math.round(chiDeg) + "°"; }

  slider.addEventListener("input", function () { chiDeg = parseFloat(slider.value); stopSpin(); draw(); });
  var spinning = false, raf;
  function stopSpin() { spinning = false; spinBtn.textContent = "▶ Rotate"; cancelAnimationFrame(raf); }
  spinBtn.addEventListener("click", function () {
    if (spinning) { stopSpin(); return; }
    spinning = true; spinBtn.textContent = "❚❚ Pause";
    (function step() {
      if (!spinning) return;
      chiDeg += 1.5; if (chiDeg > 180) chiDeg -= 360;
      slider.value = chiDeg; draw();
      raf = requestAnimationFrame(step);
    })();
  });

  draw();
  Demo.onResize(draw);
})();
