/* Shared site behaviour: theme, header, reveal-on-scroll, lightbox, filters, CV nav. */
(function () {
  "use strict";

  var root = document.documentElement;

  function store(key, val) {
    try {
      if (val === undefined) return localStorage.getItem(key);
      localStorage.setItem(key, val);
    } catch (e) { return null; }
  }

  /* ---- theme toggle ---- */
  function currentTheme() {
    var t = root.getAttribute("data-theme");
    if (t) return t;
    return window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  }
  document.querySelectorAll(".theme-toggle").forEach(function (btn) {
    btn.addEventListener("click", function () {
      var next = currentTheme() === "dark" ? "light" : "dark";
      root.setAttribute("data-theme", next);
      store("theme", next);
      document.dispatchEvent(new CustomEvent("themechange"));
    });
  });
  if (window.matchMedia) {
    var mq = window.matchMedia("(prefers-color-scheme: dark)");
    var onChange = function () { document.dispatchEvent(new CustomEvent("themechange")); };
    if (mq.addEventListener) mq.addEventListener("change", onChange);
  }

  /* ---- header border on scroll ---- */
  var header = document.querySelector(".site-header");
  if (header) {
    var onScroll = function () { header.classList.toggle("scrolled", window.scrollY > 8); };
    window.addEventListener("scroll", onScroll, { passive: true });
    onScroll();
  }

  /* ---- reveal on scroll ---- */
  var reveals = document.querySelectorAll(".reveal");
  if ("IntersectionObserver" in window) {
    var io = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (e.isIntersecting) { e.target.classList.add("in"); io.unobserve(e.target); }
      });
    }, { rootMargin: "0px 0px -8% 0px", threshold: 0.05 });
    reveals.forEach(function (el) { io.observe(el); });
  } else {
    reveals.forEach(function (el) { el.classList.add("in"); });
  }

  /* ---- figure lightbox ---- */
  var zoomables = document.querySelectorAll("img.zoomable");
  if (zoomables.length) {
    var lb = document.createElement("div");
    lb.className = "lightbox";
    lb.setAttribute("role", "dialog");
    lb.setAttribute("aria-label", "Enlarged figure");
    var lbImg = document.createElement("img");
    lbImg.alt = "";
    lb.appendChild(lbImg);
    document.body.appendChild(lb);
    var close = function () { lb.classList.remove("open"); };
    lb.addEventListener("click", close);
    document.addEventListener("keydown", function (e) { if (e.key === "Escape") close(); });
    zoomables.forEach(function (img) {
      img.addEventListener("click", function () {
        lbImg.src = img.currentSrc || img.src;
        lbImg.alt = img.alt;
        lb.classList.add("open");
      });
    });
  }

  /* ---- publication filters ---- */
  var bar = document.querySelector(".filter-bar");
  if (bar) {
    var buttons = bar.querySelectorAll("button");
    var cards = document.querySelectorAll(".pub-grid .pub-card");
    buttons.forEach(function (b) {
      b.addEventListener("click", function () {
        var f = b.getAttribute("data-filter");
        buttons.forEach(function (x) { x.setAttribute("aria-pressed", x === b ? "true" : "false"); });
        cards.forEach(function (c) {
          var tags = (c.getAttribute("data-tags") || "").split(" ");
          c.hidden = !(f === "all" || tags.indexOf(f) !== -1);
        });
      });
    });
  }

  /* ---- CV side-nav highlight ---- */
  var cvLinks = document.querySelectorAll(".cv-nav a");
  if (cvLinks.length && "IntersectionObserver" in window) {
    var map = {};
    cvLinks.forEach(function (a) { map[a.getAttribute("href").slice(1)] = a; });
    var secIO = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (e.isIntersecting) {
          cvLinks.forEach(function (a) { a.classList.remove("active"); });
          var a = map[e.target.id];
          if (a) a.classList.add("active");
        }
      });
    }, { rootMargin: "-20% 0px -70% 0px" });
    document.querySelectorAll(".cv-section[id]").forEach(function (s) { secIO.observe(s); });
  }

  /* ---- show the CV download button only if the PDF is actually there ---- */
  document.querySelectorAll("[data-requires-file]").forEach(function (el) {
    var url = el.getAttribute("data-requires-file");
    if (location.protocol === "file:") { el.hidden = false; return; }
    fetch(url, { method: "HEAD" }).then(function (r) { if (r.ok) el.hidden = false; }).catch(function () {});
  });
})();

/* Small helper shared by the demos: read theme colours from CSS variables
   and set up crisp, resizable canvases. */
window.Demo = (function () {
  "use strict";
  function css(name) {
    return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  }
  function palette() {
    return {
      bg: css("--bg"), bgSoft: css("--bg-soft"), paper: css("--paper"),
      ink: css("--ink"), inkSoft: css("--ink-soft"), muted: css("--muted"),
      line: css("--line"), lineSoft: css("--line-soft"),
      olive: css("--olive"), oliveDeep: css("--olive-deep"), oliveMid: css("--olive-mid"),
      oliveLight: css("--olive-light"), oliveWash: css("--olive-wash"),
      ochre: css("--ochre"), rust: css("--rust"), sage: css("--sage")
    };
  }
  /* Size a canvas to its CSS width with a fixed aspect ratio; returns ctx + logical size. */
  function fit(canvas, aspect) {
    var w = canvas.clientWidth || canvas.parentNode.clientWidth;
    var h = Math.round(w / aspect);
    var dpr = Math.min(window.devicePixelRatio || 1, 2);
    canvas.style.height = h + "px";
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
    var ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx: ctx, w: w, h: h };
  }
  function onResize(fn) {
    var t;
    window.addEventListener("resize", function () { clearTimeout(t); t = setTimeout(fn, 120); });
    document.addEventListener("themechange", fn);
  }
  /* Seeded RNG so demos are reproducible between visits. */
  function rng(seed) {
    var s = seed >>> 0 || 1;
    var r = function () {
      s ^= s << 13; s >>>= 0; s ^= s >>> 17; s ^= s << 5; s >>>= 0;
      return s / 4294967296;
    };
    r.normal = function () {
      var u = 1 - r(), v = r();
      return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
    };
    return r;
  }
  function alpha(color, a) {
    // works for #rrggbb
    if (color[0] === "#" && color.length === 7) {
      var n = parseInt(color.slice(1), 16);
      return "rgba(" + (n >> 16) + "," + ((n >> 8) & 255) + "," + (n & 255) + "," + a + ")";
    }
    return color;
  }
  function reducedMotion() {
    return window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  }
  /* Run fn only while the element is on screen. */
  function whenVisible(el, onShow, onHide) {
    if (!("IntersectionObserver" in window)) { onShow(); return; }
    new IntersectionObserver(function (es) {
      es.forEach(function (e) { if (e.isIntersecting) onShow(); else if (onHide) onHide(); });
    }, { threshold: 0.1 }).observe(el);
  }
  return { css: css, palette: palette, fit: fit, onResize: onResize, rng: rng, alpha: alpha, reducedMotion: reducedMotion, whenVisible: whenVisible };
})();
