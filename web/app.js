"use strict";
// Copyright (c) 2026 Sharif Mutasim
// SPDX-License-Identifier: MIT
/* Cycler Catalog web UI. Data comes from catalog.js (window.CATALOG),
   written by scripts/export_web.py. No dependencies. */

const D = window.CATALOG;
const $ = (s) => document.querySelector(s);
const TAU = 2 * Math.PI;
const SYS = ["EM", "EV", "MV", "EMV"];
const SYS_NAME = { EM: "Earth–Mars", EV: "Earth–Venus", MV: "Mars–Venus", EMV: "Venus–Earth–Mars" };
const PNAME = { V: "Venus", E: "Earth", M: "Mars" };
const FAM_NAME = { single_home: "Single-assist", sequence: "Multi-assist", home_only: "Home-only" };
const DV_FLOOR = 1e-3; // km/s; ballistic floor on the log axis

const cssVar = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const sysOf = (s) => [...new Set(s.replace(/[^VEM]/g, ""))].sort().join("");
const fmt = (v, d = 2) => (v == null || !isFinite(v) ? "—" : Number(v).toFixed(d));
const fmtInt = (v) => Number(v).toLocaleString("en-US");
const esc = (s) => String(s ?? "").replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
const groupBy = (xs, f) => xs.reduce((m, x) => ((m[f(x)] ||= []).push(x), m), {});
const j2000Date = (days) => new Date(Date.UTC(2000, 0, 1, 12) + days * 86400000);
const isoDay = (days) => j2000Date(days).toISOString().slice(0, 10);
const decYear = (iso) => { const d = new Date(iso + "T00:00:00Z"); const y = d.getUTCFullYear(); return y + (d - Date.UTC(y, 0, 1)) / (Date.UTC(y + 1, 0, 1) - Date.UTC(y, 0, 1)); };
const swatch = (sys, sq = false) => `<i class="sw${sq ? " sq" : ""}" style="background:var(--sys-${sys})"></i>`;

// ---------------------------------------------------------------------------
// Data model
// ---------------------------------------------------------------------------
// windows arrive column-wise; trajectory files are web/ephem/<label>.js (see export_web.py)
const WINDOWS = D.ephemeris.rows.map((row) => Object.fromEntries(D.ephemeris.cols.map((c, i) => [c, row[i]])));
const ephemFile = (label) => "ephem/" + label.replace(/[^A-Za-z0-9._-]/g, "_") + ".js";
const EPH_BY_LABEL = groupBy(WINDOWS, (r) => r.label);
const entry = (r) => ({ ev: r.earth_vinf_kms ?? null, c3: r.earth_c3_km2s2 ?? null, leo: r.leo_departure_dv_kms ?? null });
const IDEAL = [
  ...D.single_home.map((r) => ({
    id: r.cycler_id, fam: r.reaches_target ? "single_home" : "home_only", sys: sysOf(r.home + r.target),
    seq: r.sequence, notation: r.notation, home: r.home, period: r.period_years, dv: r.dv_per_cycle_kms,
    ballistic: !!r.ballistic, vmax: r.reaches_target ? Math.max(r.vinf_home_kms ?? 0, r.vinf_target_kms ?? 0) : r.vinf_home_kms,
    aph: r.aphelion_au, per: r.perihelion_au, nenc: r.reaches_target ? 2 : 1, raw: r, traj: r.traj,
    real: r.reaches_target ? r.real : null, ...entry(r),
  })),
  ...D.sequence.map((r) => ({
    id: r.cycler_id, fam: "sequence", sys: r.planets, seq: r.sequence, notation: "", home: r.sequence[0],
    period: r.period_years, dv: r.dv_per_cycle_kms, ballistic: !!r.ballistic, vmax: r.max_vinf_kms,
    aph: r.aphelion_max_au, per: r.perihelion_min_au, nenc: r.n_encounters, raw: r,
    traj: r.traj, real: r.real, ...entry(r),
  })),
];
const fmtC3 = (v) => (v == null ? "—" : fmt(v, 1));
const IDEAL_BY_ID = Object.fromEntries(IDEAL.map((r) => [r.id, r]));
// "M:7.42|E:10.17" -> 10.17 (highest encounter speed over the window's planets)
function maxVinf(s) {
  const v = String(s || "").split("|").map((kv) => parseFloat(kv.split(":")[1])).filter(isFinite);
  return v.length ? Math.max(...v) : null;
}
const EPH = WINDOWS.map((r, i) => ({
  id: r.window_id || `${r.label}#${i}`, label: r.label, sys: sysOf(r.sequence), seq: r.sequence,
  start: r.start_epoch, year: r.start_epoch ? decYear(r.start_epoch) : null, cycles: r.n_cycles,
  dv: r.dv_per_cycle_kms, worst: r.worst_cycle_dv_kms, total: r.total_dv_kms, dev: r.max_vinf_deviation_kms,
  conv: !!r.converged, vmax: r.max_vinf_kms, period: r.period_years, raw: r, hasTraj: !!r.has_traj, file: ephemFile(r.label),
  fam: r.label.startsWith("SQ-") ? "sequence" : "single_home",
  vmaxNum: maxVinf(r.max_vinf_kms),
  ...entry(r),
}));
const EPH_BY_ID = Object.fromEntries(EPH.map((r) => [r.id, r]));

// ---------------------------------------------------------------------------
// Two-body propagation (AU, days). Safeguarded Newton on the anomaly change.
// ---------------------------------------------------------------------------
const MU = D.meta.mu_sun_au3_day2;
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

function kepler(r0, v0, dt) {
  const r = Math.hypot(r0[0], r0[1], r0[2]);
  const alpha = 2 / r - dot(v0, v0) / MU; // 1/a
  const sig = dot(r0, v0) / Math.sqrt(MU);
  let f, g;
  if (alpha > 1e-12) {
    const a = 1 / alpha, n = Math.sqrt(MU * alpha ** 3), sa = Math.sqrt(a);
    const t = ((dt % (TAU / n)) + TAU / n) % (TAU / n);
    const M = n * t, c1 = sig / sa, c2 = 1 - r / a;
    const F = (E) => E + c1 * (1 - Math.cos(E)) - c2 * Math.sin(E) - M;
    const E = solveMono(F, (E) => 1 + c1 * Math.sin(E) - c2 * Math.cos(E), 0, TAU, M);
    f = 1 - (a / r) * (1 - Math.cos(E));
    g = t - (E - Math.sin(E)) / n;
  } else {
    const a = 1 / Math.min(alpha, -1e-12), n = Math.sqrt(MU / (-a) ** 3), sa = Math.sqrt(-a);
    const M = n * dt, c1 = sig / sa, c2 = 1 - r / a;
    const F = (H) => -H + c1 * (Math.cosh(H) - 1) + c2 * Math.sinh(H) - M;
    let hi = 1;
    while (F(hi) < 0 && hi < 100) hi *= 2;
    const H = solveMono(F, (H) => -1 + c1 * Math.sinh(H) + c2 * Math.cosh(H), 0, hi, hi / 2);
    f = 1 - (a / r) * (1 - Math.cosh(H));
    g = dt - (Math.sinh(H) - H) / n;
  }
  return [f * r0[0] + g * v0[0], f * r0[1] + g * v0[1], f * r0[2] + g * v0[2]];
}

function solveMono(F, dF, lo, hi, x) {
  for (let i = 0; i < 80; i++) {
    const fx = F(x);
    if (fx > 0) hi = x; else lo = x;
    let xn = x - fx / dF(x);
    if (!(xn > lo && xn < hi)) xn = 0.5 * (lo + hi);
    if (Math.abs(xn - x) < 1e-13) return xn;
    x = xn;
  }
  return x;
}

// ---------------------------------------------------------------------------
// Trajectory models for the orbit view
// ---------------------------------------------------------------------------
function sampleLegs(legs, rotRad = 0, tOffset = 0, out = []) {
  const c = Math.cos(rotRad), s = Math.sin(rotRad);
  for (const L of legs) {
    const [t0, x, y, z, vx, vy, vz, dt] = L;
    const r0 = [x, y, z], v0 = [vx, vy, vz];
    const steps = Math.max(40, Math.ceil(dt / 2));
    for (let k = 0; k <= steps; k++) {
      const p = kepler(r0, v0, (dt * k) / steps);
      out.push({ t: tOffset + t0 + (dt * k) / steps, x: c * p[0] - s * p[1], y: s * p[0] + c * p[1] });
    }
  }
  return out;
}

function idealModel(row, nCycles) {
  const tr = row.traj, P = tr.period_days, rot = (tr.rotate_deg * Math.PI) / 180;
  const BOD = tr.bodies || D.meta.bodies; // validation rows use the paper's idealised planets
  const ph = Object.fromEntries(Object.entries(tr.phases).map(([k, v]) => [k, (v * Math.PI) / 180]));
  const planets = Object.keys(ph);
  const samples = [];
  for (let c = 0; c < nCycles; c++) sampleLegs(tr.legs, c * rot, c * P, samples);
  const planetPos = (code, t) => {
    const b = BOD[code], th = ph[code] + (TAU * t) / b.period_days;
    return [b.a_au * Math.cos(th), b.a_au * Math.sin(th)];
  };
  const enc = [];
  for (let c = 0; c < nCycles; c++) {
    for (const [code, te] of tr.enc) enc.push({ code, t: c * P + te });
    for (const [code, te] of tr.marks || []) enc.push({ code, t: c * P + te, onPath: true });
  }
  enc.push({ code: row.home, t: nCycles * P });
  return {
    kind: "ideal", T0: 0, T1: nCycles * P, samples, planets, planetPos, home: row.home,
    homeAngle: (t) => ph[row.home] + (TAU * t) / BOD[row.home].period_days,
    enc, orbits: planets.map((c) => ({ code: c, circle: BOD[c].a_au })),
    clock: (t) => `t = ${(t / 365.25).toFixed(2)} yr · cycle ${Math.min(nCycles, Math.floor(t / P) + 1)} of ${nCycles}`,
  };
}

const PT = D.planets;
function tablePos(code, t) {
  const arr = PT[code];
  let u = (t - PT.t0_days) / PT.step_days;
  u = Math.max(0, Math.min(PT.n - 1.000001, u));
  const i = Math.floor(u), w = u - i;
  return [arr[3 * i] * (1 - w) + arr[3 * i + 3] * w, arr[3 * i + 1] * (1 - w) + arr[3 * i + 4] * w];
}

function ephemModel(tr, seq) {
  const samples = sampleLegs(tr.legs);
  const planets = [...new Set(seq)];
  const T0 = tr.enc[0][1], T1 = tr.enc[tr.enc.length - 1][1];
  const home = seq[0];
  const orbits = planets.map((c) => {
    const P = D.meta.bodies[c].period_days, pts = [];
    for (let k = 0; k <= 180; k++) pts.push(tablePos(c, T0 + (P * k) / 180));
    return { code: c, track: pts, circle: D.meta.bodies[c].a_au };
  });
  return {
    kind: "ephem", T0, T1, samples, planets, planetPos: tablePos, home,
    homeAngle: (t) => { const p = tablePos(home, t); return Math.atan2(p[1], p[0]); },
    enc: tr.enc.map(([code, t]) => ({ code, t, onPath: true })), orbits,
    clock: (t) => `${isoDay(t)} · day ${Math.round(t - T0)}`,
  };
}

// ---------------------------------------------------------------------------
// Orbit view (canvas)
// ---------------------------------------------------------------------------
class OrbitView {
  constructor(ids) {
    this.cv = $(ids.canvas); this.empty = $(ids.empty); this.btn = $(ids.play);
    this.slider = $(ids.slider); this.clockEl = $(ids.clock); this.speedSel = $(ids.speed);
    this.frame = "inertial"; this.model = null; this.t = 0; this.playing = false; this.last = 0;
    segmented($(ids.frame), (v) => { this.frame = v; this.prepare(); this.draw(); });
    this.btn.onclick = () => this.toggle();
    this.slider.oninput = () => {
      if (!this.model) return;
      this.t = this.model.T0 + (this.slider.value / 1000) * (this.model.T1 - this.model.T0);
      this.draw();
    };
    new ResizeObserver(() => this.draw()).observe(this.cv);
  }
  set(model) {
    this.model = model;
    this.empty.hidden = !!model;
    if (!model) { this.stop(); this.draw(); return; }
    this.t = model.T0;
    this.prepare();
    this.draw();
  }
  prepare() {
    const m = this.model;
    if (!m) return;
    const rotating = this.frame === "rotating";
    const tf = (x, y, t) => {
      if (!rotating) return [x, y];
      const a = -m.homeAngle(t), c = Math.cos(a), s = Math.sin(a);
      return [c * x - s * y, s * x + c * y];
    };
    this.tf = tf;
    this.path = m.samples.map((p) => tf(p.x, p.y, p.t));
    this.encPts = m.enc.map((e) => {
      if (!e.onPath) return { code: e.code, p: tf(...m.planetPos(e.code, e.t), e.t) };
      const S = m.samples, i = this.scAt(e.t), a = S[i], b = S[Math.min(i + 1, S.length - 1)];
      const w = b.t > a.t ? (e.t - a.t) / (b.t - a.t) : 0;
      return { code: e.code, p: tf(a.x + (b.x - a.x) * w, a.y + (b.y - a.y) * w, e.t) };
    });
    let R = Math.max(...m.orbits.map((o) => o.circle));
    for (const p of this.path) R = Math.max(R, Math.hypot(p[0], p[1]));
    this.R = Math.min(R, 8) * 1.06;
  }
  toggle() { this.playing ? this.stop() : this.play(); }
  play() {
    if (!this.model) return;
    if (this.t >= this.model.T1) this.t = this.model.T0;
    this.playing = true; this.btn.textContent = "❚❚"; this.last = performance.now();
    const step = (now) => {
      if (!this.playing) return;
      const m = this.model, span = m.T1 - m.T0;
      this.t += ((now - this.last) / 1000) * (span / 24) * Number(this.speedSel.value);
      this.last = now;
      if (this.t >= m.T1) { this.t = m.T1; this.stop(); }
      this.draw();
      if (this.playing) requestAnimationFrame(step);
    };
    requestAnimationFrame(step);
  }
  stop() { this.playing = false; this.btn.textContent = "▶"; }
  scAt(t) {
    const S = this.model.samples;
    let lo = 0, hi = S.length - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (S[mid].t <= t) lo = mid; else hi = mid; }
    return lo;
  }
  draw() {
    const cv = this.cv, dpr = window.devicePixelRatio || 1;
    const w = cv.clientWidth, h = cv.clientHeight;
    if (!w || !h) return;
    cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
    const g = cv.getContext("2d");
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, w, h);
    const m = this.model;
    if (!m) { this.clockEl.textContent = ""; return; }
    const col = {
      surface: cssVar("--surface-1"), grid: cssVar("--grid"), muted: cssVar("--text-muted"),
      text: cssVar("--text-primary"), text2: cssVar("--text-secondary"), traj: cssVar("--traj"), sun: cssVar("--sun"),
      V: cssVar("--pl-V"), E: cssVar("--pl-E"), M: cssVar("--pl-M"),
    };
    const cx = w / 2, cy = h / 2, s = (Math.min(w, h) / 2 - 18) / this.R;
    const X = (p) => cx + p[0] * s, Y = (p) => cy - p[1] * s;
    const rotating = this.frame === "rotating";

    // planet orbits
    g.lineWidth = 1;
    for (const o of m.orbits) {
      g.strokeStyle = col[o.code]; g.globalAlpha = 0.35;
      g.beginPath();
      if (o.track && !rotating) {
        o.track.forEach((p, i) => (i ? g.lineTo(X(p), Y(p)) : g.moveTo(X(p), Y(p))));
      } else {
        g.arc(cx, cy, o.circle * s, 0, TAU);
      }
      g.stroke();
    }
    g.globalAlpha = 1;
    // sun
    g.fillStyle = col.sun; g.beginPath(); g.arc(cx, cy, 4, 0, TAU); g.fill();

    // full path (faint) and travelled path
    const idx = this.scAt(this.t);
    g.lineJoin = g.lineCap = "round";
    g.strokeStyle = col.traj; g.globalAlpha = 0.28; g.lineWidth = 1.25;
    g.beginPath();
    this.path.forEach((p, i) => (i ? g.lineTo(X(p), Y(p)) : g.moveTo(X(p), Y(p))));
    g.stroke();
    g.globalAlpha = 1; g.lineWidth = 2;
    g.beginPath();
    for (let i = 0; i <= idx; i++) { const p = this.path[i]; i ? g.lineTo(X(p), Y(p)) : g.moveTo(X(p), Y(p)); }
    const S = m.samples, a = S[idx], b = S[Math.min(idx + 1, S.length - 1)];
    const w8 = b.t > a.t ? (this.t - a.t) / (b.t - a.t) : 0;
    const sc = this.tf(a.x + (b.x - a.x) * w8, a.y + (b.y - a.y) * w8, this.t);
    g.lineTo(X(sc), Y(sc));
    g.stroke();

    // encounter markers
    for (const e of this.encPts) {
      g.strokeStyle = col[e.code]; g.lineWidth = 1.5;
      g.beginPath(); g.arc(X(e.p), Y(e.p), 3.5, 0, TAU); g.stroke();
    }
    // planets now
    g.font = "12px system-ui, -apple-system, 'Segoe UI', sans-serif";
    for (const c of m.planets) {
      const p = this.tf(...m.planetPos(c, this.t), this.t);
      g.fillStyle = col.surface; g.beginPath(); g.arc(X(p), Y(p), 7.5, 0, TAU); g.fill();
      g.fillStyle = col[c]; g.beginPath(); g.arc(X(p), Y(p), 5.5, 0, TAU); g.fill();
      g.fillStyle = col.text2; g.fillText(PNAME[c], X(p) + 9, Y(p) - 7);
    }
    // spacecraft
    g.fillStyle = col.surface; g.beginPath(); g.arc(X(sc), Y(sc), 6.5, 0, TAU); g.fill();
    g.fillStyle = col.text; g.beginPath(); g.arc(X(sc), Y(sc), 4.5, 0, TAU); g.fill();

    // scale bar + frame label
    const L = s * 1;
    g.strokeStyle = col.muted; g.lineWidth = 1;
    g.beginPath(); g.moveTo(14, h - 16); g.lineTo(14 + L, h - 16); g.stroke();
    g.fillStyle = col.muted; g.fillText("1 AU", 14, h - 22);
    g.fillText(rotating ? `Rotating with ${PNAME[m.home]}` : "Inertial (J2000 ecliptic)", 14, 18);

    this.clockEl.textContent = m.clock(this.t);
    this.slider.value = Math.round(((this.t - m.T0) / (m.T1 - m.T0)) * 1000);
  }
}

// ---------------------------------------------------------------------------
// Scatter plot (SVG)
// ---------------------------------------------------------------------------
function niceTicks(lo, hi, n = 6) {
  const span = hi - lo || 1, step0 = span / n, mag = 10 ** Math.floor(Math.log10(step0));
  const step = [1, 2, 2.5, 5, 10].map((k) => k * mag).find((s) => span / s <= n) || 10 * mag;
  const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(+v.toFixed(10));
  return out;
}

// One axis: linear (padded data range, or fixed min/max) or log10 with values at or
// below DV_FLOOR drawn on the floor (used for ΔV, where ≤ 1 m/s means ballistic).
function axisScale(vals, { log = false, min, max, fmt: f } = {}) {
  if (log) {
    const top = Math.max(1, ...vals);
    const l0 = Math.log10(DV_FLOOR), l1 = Math.ceil(Math.log10(top));
    const ticks = [];
    for (let k = l0; k <= l1; k++) ticks.push(10 ** k);
    return { T: (v) => (Math.log10(Math.max(v, DV_FLOOR)) - l0) / (l1 - l0 || 1), ticks,
      label: (v) => (v <= DV_FLOOR ? "≤0.001" : v >= 1 ? String(v) : String(+v.toPrecision(1))), clipped: 0 };
  }
  let lo = min ?? (vals.length ? Math.min(...vals) : 0), hi = max ?? (vals.length ? Math.max(...vals) : 1);
  // A handful of extreme points (e.g. C3 of 10,000 km²/s² against a bulk under 2,500) would
  // squash everything else into one strip: stop the axis just past the 99th percentile and
  // pin the few beyond it to the edge.
  let clipped = 0;
  if (max == null && vals.length >= 50) {
    const p99 = [...vals].sort((a, b) => a - b)[Math.floor(0.99 * (vals.length - 1))];
    if (hi - lo > 2.5 * (p99 - lo) && p99 > lo) { hi = p99 + (p99 - lo) * 0.08; clipped = vals.filter((v) => v > hi).length; }
  }
  if (min == null || max == null) { const pad = (hi - lo) * 0.04 || 0.5; lo -= pad; if (!clipped) hi += pad; }
  return { T: (v) => Math.min((v - lo) / (hi - lo), 1), ticks: (n) => niceTicks(lo, hi, n), label: f || ((v) => String(v)), clipped };
}

function scatter(el, pts, o) {
  const w = el.clientWidth, h = el.clientHeight;
  if (!w || !h) return;
  const m = { l: 56, r: 14, t: 10, b: 42 };
  const ax = axisScale(pts.map((p) => p.x).filter(isFinite), { log: o.logX, min: o.xMin, max: o.xMax, fmt: o.xFmt });
  const ay = axisScale(pts.map((p) => p.y).filter(isFinite), { log: o.logY });
  const PX = (v) => m.l + ax.T(v) * (w - m.l - m.r);
  const PY = (v) => h - m.b - ay.T(v) * (h - m.t - m.b);
  const xticks = typeof ax.ticks === "function" ? ax.ticks(Math.max(3, Math.floor(w / 90))) : ax.ticks;
  const yticks = typeof ay.ticks === "function" ? ay.ticks(6) : ay.ticks;
  let svg = `<svg viewBox="0 0 ${w} ${h}" role="img" aria-label="${esc(o.aria)}">`;
  for (const v of yticks) {
    const y = PY(v);
    svg += `<line class="gridline" x1="${m.l}" x2="${w - m.r}" y1="${y}" y2="${y}"/><text class="tick" x="${m.l - 8}" y="${y + 4}" text-anchor="end">${ay.label(v)}</text>`;
  }
  for (const v of xticks) {
    const x = PX(v);
    if (o.logX) svg += `<line class="gridline" x1="${x}" x2="${x}" y1="${m.t}" y2="${h - m.b}"/>`;
    svg += `<text class="tick" x="${x}" y="${h - m.b + 16}" text-anchor="middle">${ax.label(v)}</text>`;
  }
  svg += `<line class="axisline" x1="${m.l}" x2="${w - m.r}" y1="${h - m.b}" y2="${h - m.b}"/>`;
  svg += `<text class="axis-title" x="${(m.l + w - m.r) / 2}" y="${h - 6}" text-anchor="middle">${esc(o.xLabel)}</text>`;
  svg += `<text class="axis-title" transform="translate(13 ${(m.t + h - m.b) / 2}) rotate(-90)" text-anchor="middle">${esc(o.yLabel)}</text>`;
  const off = [ax.clipped && `${ax.clipped} beyond the right edge`, ay.clipped && `${ay.clipped} above the top`].filter(Boolean);
  if (off.length) svg += `<text class="tick" x="${w - m.r}" y="${m.t + 10}" text-anchor="end">${fmtInt(ax.clipped + ay.clipped)} extreme point${ax.clipped + ay.clipped > 1 ? "s" : ""} pinned to the edge (${off.join(", ")})</text>`;
  svg += `</svg><canvas class="pts"></canvas>`;
  el.innerHTML = svg;

  // Points go on a canvas: tens of thousands of SVG nodes are too slow to rebuild on
  // every filter change. Hover uses a coarse grid so lookups stay cheap.
  const cv = el.querySelector("canvas"), dpr = window.devicePixelRatio || 1;
  cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
  const g = cv.getContext("2d");
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  const n = pts.length, r = n > 6000 ? 2.6 : n > 1500 ? 3.4 : 4.5;
  const surface = cssVar("--surface-1"), ink = cssVar("--text-primary");
  const color = Object.fromEntries(SYS.map((k) => [k, cssVar(`--sys-${k}`)]));
  const CELL = 16, grid = new Map(), placed = [];
  let sel = null;
  const mark = (x, y, rr, p) => {
    g.beginPath();
    if (p.shape === "square") g.rect(x - rr, y - rr, 2 * rr, 2 * rr); else g.arc(x, y, rr, 0, TAU);
  };
  g.lineWidth = n > 6000 ? 0.75 : 1.25;
  g.strokeStyle = surface;
  for (const p of pts) {
    if (!isFinite(p.x) || !isFinite(p.y)) continue;
    const x = PX(p.x), y = PY(p.y);
    const q = { x, y, p };
    placed.push(q);
    const key = `${Math.floor(x / CELL)},${Math.floor(y / CELL)}`;
    (grid.get(key) || grid.set(key, []).get(key)).push(q);
    if (p.id === o.selected) { sel = q; continue; }
    g.fillStyle = color[p.sys];
    mark(x, y, r, p); g.fill(); g.stroke();
  }
  if (sel) {   // selected point on top, larger, with an ink ring
    g.fillStyle = color[sel.p.sys]; g.strokeStyle = surface; g.lineWidth = 3;
    mark(sel.x, sel.y, 6.5, sel.p); g.fill(); g.stroke();
    g.strokeStyle = ink; g.lineWidth = 2; mark(sel.x, sel.y, 7.5, sel.p); g.stroke();
  }

  const s = el.querySelector("svg");
  const nearest = (ev) => {
    const b = s.getBoundingClientRect(), mx = ev.clientX - b.left, my = ev.clientY - b.top;
    const cx = Math.floor(mx / CELL), cy = Math.floor(my / CELL);
    let best = null, bd = 144; // 12px hit radius
    for (let i = -1; i <= 1; i++) for (let j = -1; j <= 1; j++) {
      for (const q of grid.get(`${cx + i},${cy + j}`) || []) {
        const d = (q.x - mx) ** 2 + (q.y - my) ** 2;
        if (d < bd) { bd = d; best = q; }
      }
    }
    return best;
  };
  s.onmousemove = (ev) => { const q = nearest(ev); q ? showTip(ev, o.tip(q.p.row)) : hideTip(); s.style.cursor = q ? "pointer" : "default"; };
  s.onmouseleave = hideTip;
  s.onclick = (ev) => { const q = nearest(ev); if (q) o.onPick(q.p.row); };
  return placed.length;
}

// Coalesce bursts of calls (slider drags, resize events) into one per animation frame.
function perFrame(fn) {
  let queued = false;
  return () => { if (queued) return; queued = true; requestAnimationFrame(() => { queued = false; fn(); }); };
}

const tipEl = $("#tip");
function showTip(ev, html) {
  tipEl.innerHTML = html; tipEl.hidden = false;
  const r = tipEl.getBoundingClientRect();
  let x = ev.clientX + 14, y = ev.clientY + 14;
  if (x + r.width > innerWidth - 8) x = ev.clientX - r.width - 14;
  if (y + r.height > innerHeight - 8) y = ev.clientY - r.height - 14;
  tipEl.style.left = x + "px"; tipEl.style.top = y + "px";
}
function hideTip() { tipEl.hidden = true; }
const tipRow = (k, v) => `<div class="t-row"><span>${k}</span><b>${v}</b></div>`;

// ---------------------------------------------------------------------------
// Tables
// ---------------------------------------------------------------------------
function table(el, cols, rows, st, selected, onPick) {
  const sorted = [...rows];
  const c = cols.find((c) => c.key === st.sort);
  if (c) {
    const get = c.sortGet || c.get;
    sorted.sort((a, b) => {
      const va = get(a), vb = get(b);
      if (va == null || (typeof va === "number" && !isFinite(va))) return 1;
      if (vb == null || (typeof vb === "number" && !isFinite(vb))) return -1;
      return (va < vb ? -1 : va > vb ? 1 : 0) * st.dir;
    });
  }
  const shown = sorted.slice(0, st.limit);
  let html = "<thead><tr>" + cols.map((c) =>
    `<th data-k="${c.key}" class="${c.num ? "num" : ""}" ${c.key === st.sort ? `aria-sort="${st.dir > 0 ? "ascending" : "descending"}"` : ""}>${c.label}</th>`).join("") + "</tr></thead><tbody>";
  for (const r of shown) {
    html += `<tr data-id="${esc(r.id)}" class="${r.id === selected ? "sel" : ""}">` +
      cols.map((c) => `<td class="${c.num ? "num" : ""}">${c.html ? c.html(r) : esc(c.fmt ? c.fmt(c.get(r)) : c.get(r))}</td>`).join("") + "</tr>";
  }
  el.innerHTML = html + "</tbody>";
  el.querySelectorAll("th").forEach((th) => (th.onclick = () => {
    const k = th.dataset.k;
    st.dir = st.sort === k ? -st.dir : 1; st.sort = k; st.rerender();
  }));
  el.querySelectorAll("tbody tr").forEach((tr) => (tr.onclick = () => onPick(tr.dataset.id)));
  return { shown: shown.length, total: rows.length };
}

function downloadCSV(name, rows) {
  if (!rows.length) return;
  const keys = Object.keys(rows[0]).filter((k) => k !== "traj");
  const q = (v) => { const s = v == null ? "" : String(v); return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s; };
  const csv = [keys.join(","), ...rows.map((r) => keys.map((k) => q(r[k])).join(","))].join("\n");
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([csv], { type: "text/csv" }));
  a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

// ---------------------------------------------------------------------------
// Small UI helpers
// ---------------------------------------------------------------------------
function segmented(el, onChange) {
  el.querySelectorAll("button").forEach((b) => (b.onclick = () => {
    el.querySelectorAll("button").forEach((x) => x.setAttribute("aria-pressed", String(x === b)));
    onChange(b.dataset.v);
  }));
}
function setSegmented(el, v) {
  el.querySelectorAll("button").forEach((x) => x.setAttribute("aria-pressed", String(x.dataset.v === v)));
}
function sysChips(el, set, onChange, systems = SYS) {
  el.innerHTML = systems.map((s) => `<button class="chip" data-s="${s}" aria-pressed="${set.has(s)}">${swatch(s)}${SYS_NAME[s]}</button>`).join("");
  el.querySelectorAll(".chip").forEach((b) => (b.onclick = () => {
    const s = b.dataset.s;
    set.has(s) ? set.delete(s) : set.add(s);
    b.setAttribute("aria-pressed", String(set.has(s)));
    onChange();
  }));
}
function legend(el, systems, shapes) {
  el.innerHTML = systems.map((s) => `<span>${swatch(s)}${SYS_NAME[s]}</span>`).join("") +
    (shapes ? `<span style="margin-left:8px">${swatch("EMV")} circle = single-assist</span><span>${swatch("EMV", true)} square = multi-assist</span>` : "");
}
const kv = (pairs) => `<div class="kv">${pairs.map(([k, v]) => `<div><span>${k}</span><b>${v}</b></div>`).join("")}</div>`;
const badge = (ok, yes, no) => `<span class="badge ${ok ? "good" : "bad"}">${ok ? yes : no}</span>`;

// ---------------------------------------------------------------------------
// Real-planet trajectories: one small script per cycler, loaded on demand
// ---------------------------------------------------------------------------
const EPH_TRAJ = {};
const EPH_WAIT = {};
window.__ephemLoaded = (label, trajs) => {
  EPH_TRAJ[label] = trajs;
  (EPH_WAIT[label] || []).forEach((f) => f(trajs));
  delete EPH_WAIT[label];
};
function loadEphem(label, file) {
  return new Promise((resolve) => {
    if (EPH_TRAJ[label]) return resolve(EPH_TRAJ[label]);
    const first = !EPH_WAIT[label];
    (EPH_WAIT[label] ||= []).push(resolve);
    if (first) {
      const s = document.createElement("script");
      s.src = file;
      s.onerror = () => { (EPH_WAIT[label] || []).forEach((f) => f(null)); delete EPH_WAIT[label]; };
      document.head.appendChild(s);
    }
  });
}
async function windowTraj(w) {
  if (!w || !w.hasTraj) return null;
  const all = await loadEphem(w.label, w.file);
  return all ? all[w.id] : null;
}
const bestWindow = (label) => (EPH_BY_LABEL[label] || []).filter((w) => w.converged && w.has_traj)
  .sort((a, b) => a.dv_per_cycle_kms - b.dv_per_cycle_kms)[0];

// ---------------------------------------------------------------------------
// Catalog tab (ideal / real-planet toggle)
// ---------------------------------------------------------------------------
let MODE = "ideal";
const num = (v) => (v == null ? NaN : v);
const M_ = {
  dv: (r) => (MODE === "ideal" ? r.dv : num(r.real?.real_dv_median_kms)),
  vmax: (r) => (MODE === "ideal" ? r.vmax : num(r.real?.real_vinf_max_median_kms)),
  c3: (r) => (MODE === "ideal" ? num(r.c3) : num(r.real?.real_c3_median_km2s2)),
  leo: (r) => (MODE === "ideal" ? num(r.leo) : num(r.real?.real_leo_dv_median_kms)),
  dev: (r) => (MODE === "ideal" ? NaN : num(r.real?.real_vinf_deviation_median_kms)),
  ball: (r) => (MODE === "ideal" ? r.ballistic : (r.real?.real_ballistic_windows || 0) > 0),
};
const hasReal = (r) => (r.real?.n_converged || 0) > 0;
const fmtDv = (v) => (v == null || !isFinite(v) ? "—" : v <= DV_FLOOR ? "≤0.001" : fmt(v, 3));
const MODE_NOTE = {
  ideal: "<b>Ideal model:</b> circular, coplanar planet orbits, in which a cycler repeats exactly. This is how cyclers are defined and named, and the model of the published tables. <button class=\"link\" data-goto=\"method#m-models\">Why two models?</button>",
  real: "<b>Real planets:</b> each cycler re-solved on the real, eccentric and inclined orbits (VSOP2013) for every launch window 2030–2060. Values are the <b>median over windows</b>; the lowest and highest window are in the tooltip and details, and can fall either side of the ideal value. <button class=\"link\" data-goto=\"method#m-real\">How?</button>",
};

const IF_DEFAULT = { q: "", fam: "all", dvMax: 2.3, vinfMax: 40, perMax: 16, c3Max: 150, ball: false, earth: false };
const IF = { ...IF_DEFAULT, sys: new Set(SYS) };
const ideal = { sel: null, st: { sort: "dv", dir: 1, limit: 200 } };
const AX = {
  period: { key: "period", label: "Period (years)", get: (r) => r.period },
  vmax: { key: "vmax", label: "Max v∞ (km/s)", get: (r) => M_.vmax(r) },
  aph: { key: "aph", label: "Aphelion (AU)", get: (r) => r.aph },
  per: { key: "per", label: "Perihelion (AU)", get: (r) => r.per },
  dv: { key: "dv", label: "ΔV per cycle (km/s)", get: (r) => M_.dv(r), log: true },
  c3: { key: "c3", label: "Launch C3 at Earth (km²/s²)", get: (r) => M_.c3(r) },
  leo: { key: "leo", label: "ΔV from 200 km LEO (km/s)", get: (r) => M_.leo(r) },
  dev: { key: "dev", label: "v∞ drift from ideal (km/s)", get: (r) => M_.dev(r) },
};
const orbit = new OrbitView({ canvas: "#orbit", empty: "#orbit-empty", play: "#o-play", slider: "#o-time", clock: "#o-clock", frame: "#o-frame", speed: "#o-speed" });

function idealFiltered() {
  const q = IF.q.trim().toUpperCase();
  const within = (v, max, any) => max >= any || (isFinite(v) && v <= max);
  return IDEAL.filter((r) =>
    (IF.fam === "all" ? r.fam !== "home_only" : r.fam === IF.fam) && IF.sys.has(r.sys) &&
    within(M_.dv(r), 10 ** IF.dvMax, 10 ** 2.3) && within(M_.vmax(r), IF.vinfMax, 40) &&
    (IF.perMax >= 16 || r.period <= IF.perMax) && within(M_.c3(r), IF.c3Max, 150) &&
    (!IF.ball || M_.ball(r)) && (!IF.earth || r.ev != null) &&
    (!q || r.id.toUpperCase().includes(q) || r.seq.toUpperCase().includes(q) || r.notation.toUpperCase() === q));
}

function idealCols() {
  const base = [
    { key: "id", label: "ID", get: (r) => r.id },
    { key: "notation", label: "Notation", get: (r) => r.notation },
    { key: "fam", label: "Family", get: (r) => FAM_NAME[r.fam] },
    { key: "sys", label: "Planets", get: (r) => r.sys, html: (r) => `${swatch(r.sys, r.fam === "sequence")} ${SYS_NAME[r.sys]}` },
    { key: "seq", label: "Sequence", get: (r) => r.seq },
    { key: "period", label: "Period (yr)", num: true, get: (r) => r.period, fmt: (v) => fmt(v, 2) },
  ];
  if (MODE === "ideal") return base.concat([
    { key: "dv", label: "ΔV / cycle (km/s)", num: true, get: (r) => r.dv, fmt: fmtDv },
    { key: "ballistic", label: "Ballistic", get: (r) => (r.ballistic ? 1 : 0), html: (r) => (r.ballistic ? badge(true, "yes") : "") },
    { key: "vmax", label: "Max v∞ (km/s)", num: true, get: (r) => r.vmax, fmt: (v) => fmt(v, 2) },
    { key: "c3", label: "Launch C3 (km²/s²)", num: true, get: (r) => r.c3, fmt: fmtC3 },
    { key: "leo", label: "ΔV from LEO (km/s)", num: true, get: (r) => r.leo, fmt: (v) => fmt(v, 2) },
    { key: "aph", label: "Aphelion (AU)", num: true, get: (r) => r.aph, fmt: (v) => fmt(v, 3) },
    { key: "per", label: "Perihelion (AU)", num: true, get: (r) => r.per, fmt: (v) => fmt(v, 3) },
    { key: "rdv", label: "Real ΔV, median", num: true, get: (r) => r.real?.real_dv_median_kms, fmt: fmtDv },
  ]);
  return base.concat([
    { key: "dv", label: "ΔV / cycle, median", num: true, get: (r) => r.real?.real_dv_median_kms, fmt: fmtDv },
    { key: "best", label: "Best window", num: true, get: (r) => r.real?.real_dv_best_kms, fmt: fmtDv },
    { key: "ballistic", label: "Ballistic windows", num: true, get: (r) => r.real?.real_ballistic_windows ?? null,
      html: (r) => (r.real ? `${r.real.real_ballistic_windows} / ${r.real.n_converged}` : "—") },
    { key: "vmax", label: "Max v∞, median", num: true, get: (r) => r.real?.real_vinf_max_median_kms, fmt: (v) => fmt(v, 2) },
    { key: "dev", label: "v∞ drift, median", num: true, get: (r) => r.real?.real_vinf_deviation_median_kms, fmt: (v) => fmt(v, 2) },
    { key: "c3", label: "Launch C3, median", num: true, get: (r) => r.real?.real_c3_median_km2s2, fmt: fmtC3 },
    { key: "idv", label: "Ideal ΔV", num: true, get: (r) => r.dv, fmt: fmtDv },
    { key: "win", label: "Windows", num: true, get: (r) => r.real?.n_windows ?? 0, html: (r) => (r.real ? `${r.real.n_converged} / ${r.real.n_windows}` : "—") },
  ]);
}

const idealName = (r) => (r.notation ? `${r.notation} ${r.seq}` : r.seq);

function idealTip(r) {
  let h = `<div class="t-head">${swatch(r.sys, r.fam === "sequence")}${esc(r.id)}</div>` +
    tipRow(r.notation ? "Notation" : "Sequence", esc(idealName(r))) + tipRow("Period", fmt(r.period) + " yr");
  if (MODE === "ideal") {
    h += tipRow("ΔV / cycle", fmtDv(r.dv) + " km/s") + tipRow("Max v∞", fmt(r.vmax) + " km/s") +
      tipRow("Launch C3", r.c3 == null ? "no Earth encounter" : fmtC3(r.c3) + " km²/s²");
    if (r.real) h += tipRow("Real planets, median", fmtDv(r.real.real_dv_median_kms) + " km/s");
  } else if (r.real) {
    const x = r.real;
    h += tipRow("ΔV / cycle, median", fmtDv(x.real_dv_median_kms) + " km/s") +
      tipRow("best / worst window", `${fmtDv(x.real_dv_best_kms)} / ${fmtDv(x.real_dv_worst_kms)}`) +
      tipRow("Ballistic windows", `${x.real_ballistic_windows} of ${x.n_converged}`) +
      tipRow("v∞ drift, median", fmt(x.real_vinf_deviation_median_kms) + " km/s") +
      tipRow("Ideal ΔV / cycle", fmtDv(r.dv) + " km/s");
  } else {
    h += tipRow("Real planets", "no converged window");
  }
  return h;
}

// Geometry (period, aphelion, perihelion) is the ideal value in both modes; the rest
// are medians over launch windows in Real-planets mode.
const GEOMETRY = new Set(["period", "aph", "per"]);
const shortName = (a) => a.label.replace(/ \(.*\)$/, "").replace(/^v∞/, "V∞");
const axisTitle = (a) => a.label + (MODE === "real" && !GEOMETRY.has(a.key) ? ", median" : "") + (a.log ? " — log" : "");

function renderIdeal() {
  const rows = idealFiltered();
  const xa = AX[$("#ax-x").value], ya = AX[$("#ax-y").value];
  $("#ideal-plot-title").textContent = `${shortName(ya)} vs ${shortName(xa)}`;
  legend($("#ideal-legend"), SYS.filter((s) => IF.sys.has(s)), true);
  const pts = rows.map((r) => ({ id: r.id, x: xa.get(r), y: ya.get(r), sys: r.sys, shape: r.fam === "sequence" ? "square" : "circle", row: r }));
  scatter($("#ideal-scatter"), pts, {
    xLabel: axisTitle(xa), yLabel: axisTitle(ya),
    logX: !!xa.log, logY: !!ya.log, selected: ideal.sel?.id, aria: "Scatter of cyclers", onPick: selectIdeal, tip: idealTip,
  });
  const missing = pts.filter((p) => !isFinite(p.x) || !isFinite(p.y)).length;
  $("#ideal-plot-note").innerHTML = (MODE === "ideal"
    ? "Each point is one cycler in the ideal (circular, coplanar) model. "
    : "Each point is one cycler: the median of its real-planet launch windows. Period, aphelion and perihelion are the ideal values. ") +
    "ΔV ≤ 1 m/s (ballistic) is drawn on the floor of the log axis. Click a point to see its trajectory." +
    (missing ? ` ${missing} cycler${missing > 1 ? "s" : ""} not plotted (no value on these axes${MODE === "real" ? ", e.g. no converged real window or no Earth encounter" : ""}).` : "");
  ideal.st.rerender = renderIdeal;
  const n = table($("#ideal-table"), idealCols(), rows, ideal.st, ideal.sel?.id, (id) => selectIdeal(IDEAL_BY_ID[id]));
  $("#ideal-count").textContent = `${fmtInt(n.total)} match · showing ${fmtInt(n.shown)}`;
  $("#ideal-more").hidden = n.shown >= n.total;
  $("#f-dv-v").textContent = IF.dvMax >= 2.3 ? "any" : `≤ ${+(10 ** IF.dvMax).toPrecision(2)} km/s`;
  $("#f-vinf-v").textContent = IF.vinfMax >= 40 ? "any" : `≤ ${IF.vinfMax} km/s`;
  $("#f-per-v").textContent = IF.perMax >= 16 ? "any" : `≤ ${IF.perMax} yr`;
  $("#f-c3-v").textContent = IF.c3Max >= 150 ? "any" : `≤ ${IF.c3Max} km²/s²`;
}

let selToken = 0;
async function selectIdeal(r) {
  ideal.sel = r;
  const token = ++selToken;
  renderIdeal();
  $("#detail").innerHTML = r ? idealDetail(r) : "";
  const b = $("#detail .to-ephem");
  if (b) b.onclick = () => openEphemFor(r.id);
  $("#o-cycles-l").hidden = MODE === "real";
  if (!r) { orbit.set(null); $("#orbit-title").textContent = "Trajectory"; return; }
  if (MODE === "real") {
    const w = bestWindow(r.id);
    $("#orbit-title").textContent = w ? `${idealName(r)} — best real window, launch ${w.start_epoch}` : `${idealName(r)} — no real-planet window`;
    const tr = await windowTraj(w && EPH_BY_ID[w.window_id]);
    if (token !== selToken) return;
    // single-assist sequences read "E(M)"; the model wants planet letters only
    orbit.set(tr ? ephemModel(tr, r.seq.replace(/[^VEM]/g, "")) : null);
    if (!tr) $("#orbit-empty").textContent = "No converged real-planet window for this cycler";
  } else {
    $("#orbit-title").textContent = `${r.id} — ${idealName(r)} (ideal model)`;
    orbit.set(r.traj ? idealModel(r, Number($("#o-cycles").value)) : null);
  }
}

function boardingHtml(ev, c3, leo, where) {
  if (ev == null) return `<h3>Boarding from Earth</h3><p class="note">This cycler never meets Earth, so it can't be boarded by a direct launch; a crew would first need a transfer to one of its planets.</p>`;
  return `<h3>Boarding from Earth</h3>` + kv([
    ["Earth v∞", fmt(ev) + " km/s"], ["Launch C3", fmtC3(c3) + " km²/s²"], ["ΔV from 200 km LEO", fmt(leo) + " km/s"],
  ]) + `<p class="note">One-time cost to join the cycler: a launch straight onto its outbound hyperbola at ${where}. C3 = v∞²; the LEO burn is √(v∞² + 2μ/r) − √(μ/r) from a 200 km circular orbit.</p>`;
}

function compareHtml(r) {
  const x = r.real;
  const row = (k, i, med, best, worst) => `<tr><td>${k}</td><td class="num">${i}</td><td class="num">${med}</td><td class="num">${best}</td><td class="num">${worst}</td></tr>`;
  let h = `<h3>Ideal model vs real planets</h3>`;
  if (!x) return h + `<p class="note">${r.fam === "home_only" ? "Home-only orbits are not cyclers, so they were not re-solved on the real planets."
    : r.sys === "EMV" ? "No launch date between 2030 and 2060 lines up all three real planets within 15° of this cycler's ideal geometry, so it has no real-planet windows."
    : "No real-planet results for this cycler."}</p>`;
  const n = x.n_converged;
  // Lowest/highest are extremes over n launch windows, so they bracket the ideal value
  // more often than not; say so whenever a real value comes in under the ideal one.
  const under = (real, ideal) => real != null && ideal != null && real < ideal - 1e-9;
  const beats = under(x.real_dv_best_kms, r.dv) || under(x.real_c3_best_km2s2, r.c3);
  h += `<div class="table-wrap"><table class="vt"><thead><tr><th></th><th class="num">Ideal</th><th class="num">Real, median</th><th class="num">Real, lowest</th><th class="num">Real, highest</th></tr></thead><tbody>` +
    row("ΔV / cycle (km/s)", fmtDv(r.dv), fmtDv(x.real_dv_median_kms), fmtDv(x.real_dv_best_kms), fmtDv(x.real_dv_worst_kms)) +
    row("Max v∞ (km/s)", fmt(r.vmax), fmt(x.real_vinf_max_median_kms), "", "") +
    row("Launch C3 (km²/s²)", fmtC3(r.c3), fmtC3(x.real_c3_median_km2s2), fmtC3(x.real_c3_best_km2s2), fmtC3(x.real_c3_highest_km2s2)) +
    row("ΔV from LEO (km/s)", fmt(r.leo), fmt(x.real_leo_dv_median_kms), "", "") +
    `</tbody></table></div>` +
    `<p class="note">Real values are taken over ${n} launch window${n === 1 ? "" : "s"} (2030–2060). On the real, eccentric and tilted planet orbits the geometry changes from one launch date to the next, and the optimiser only minimises ΔV per cycle, so encounter speeds (and with them C3) come out above or below the ideal value window by window.` +
    (beats && n > 1 ? ` That is why the lowest real value can be under the ideal one; the median is the fair comparison.` : "") +
    ` Ideal C3 uses ${r.fam === "sequence" ? "the cycler's lowest-v∞ Earth flyby" : "its Earth flyby"}; real C3 uses the first Earth encounter of each window.</p>` + kv([
      ["Launch windows (converged)", `${x.n_converged} of ${x.n_windows}`],
      ["Ballistic windows (≤ 1 m/s per cycle)", x.real_ballistic_windows],
      ["Windows under 0.1 km/s per cycle", x.real_cheap_windows],
      ["v∞ drift from ideal, median", fmt(x.real_vinf_deviation_median_kms) + " km/s"],
    ]) + `<p class="note">${(x.real_vinf_deviation_median_kms ?? 0) > 2 ? "The real-planet encounter speeds differ from the ideal ones by more than 2 km/s, so the cheap real windows may be a different cycler that happens to fit those dates. " : ""}Real-planet costs come from a local search (encounter dates within ±45 days, powered flybys only), so they are upper bounds. <button class="link to-ephem">Open all windows →</button></p>`;
  return h;
}

function singleHomeGeometry(x) {
  const branch = x.branch === "hi" ? "long" : x.branch === "lo" ? "short" : "unique";
  return kv([
    ["Notation", `${x.notation} (${x.n_synodic} synodic, ${branch}-period, ${x.n_revs} rev)`],
    [`v∞ at ${PNAME[x.home]}`, fmt(x.vinf_home_kms) + " km/s"], [`v∞ at ${PNAME[x.target]}`, fmt(x.vinf_target_kms) + " km/s"],
    [`${x.home}→${x.target} transfer`, x.transfer_days == null ? "—" : fmt(x.transfer_days, 1) + " d"],
    ["Time to aphelion", fmt(x.aphelion_days, 1) + " d"],
    ["Turn required / max", `${fmt(x.turn_required_deg, 1)}° / ${fmt(x.turn_max_deg, 1)}°`],
    ["Semi-major axis / e", `${fmt(x.a_au, 3)} AU / ${fmt(x.e, 3)}`],
  ]);
}

function idealDetail(r) {
  const x = r.raw;
  let html = kv([
    ["Family", FAM_NAME[r.fam]], ["Planets", SYS_NAME[r.sys]], ["Period", fmt(r.period, 3) + " yr"],
    ["ΔV / cycle (ideal)", fmtDv(r.dv) + " km/s"], ["Max v∞ (ideal)", fmt(r.vmax) + " km/s"],
    ["Aphelion", fmt(r.aph, 3) + " AU"], ["Perihelion", fmt(r.per, 3) + " AU"], ["Ballistic (ideal)", badge(r.ballistic, "yes", "no — powered flybys")],
  ]);
  if (r.fam === "home_only") {
    const B = D.meta.bodies, out = B[x.target].a_au > B[x.home].a_au;
    html += `<p class="note"><b>Home-only orbit.</b> The arc never reaches ${PNAME[x.target]}'s orbit: its ${out ? "aphelion" : "perihelion"} is ${fmt(out ? r.aph : r.per, 3)} AU, against ${fmt(B[x.target].a_au, 3)} AU for ${PNAME[x.target]}. ` +
      `It is a valid periodic orbit of ${PNAME[x.home]} that the search turns up, but it only ever revisits ${PNAME[x.home]}, so it is kept apart from the cyclers. ` +
      (x.reaches_target_eccentric ? `It does overlap the range of ${PNAME[x.target]}'s real, eccentric orbit, so a real-ephemeris version might still meet it. ` : "") +
      `${PNAME[x.target]} is drawn where it comes closest.</p>`;
  }
  html += compareHtml(r);
  if (r.fam !== "sequence") {
    html += `<h3>Single-assist geometry (ideal)</h3>` + singleHomeGeometry(x);
  } else {
    const sp = (s) => String(s).split("|");
    const v = sp(x.vinf_kms), f = sp(x.flyby_dv_kms), L = sp(x.leg_days), R = sp(x.leg_revs);
    let t = 0;
    html += `<h3>Encounters, one ideal cycle</h3><div class="table-wrap"><table><thead><tr><th>#</th><th>Planet</th><th class="num">t (days)</th><th class="num">v∞ (km/s)</th><th class="num">Flyby ΔV (km/s)</th><th class="num">Next leg (days)</th><th class="num">Revs</th></tr></thead><tbody>` +
      [...x.sequence].map((c, i) => { const row = `<tr><td>${i + 1}</td><td>${PNAME[c]}</td><td class="num">${fmt(t, 1)}</td><td class="num">${v[i]}</td><td class="num">${fmt(+f[i], 3)}</td><td class="num">${L[i]}</td><td class="num">${R[i]}</td></tr>`; t += +L[i]; return row; }).join("") +
      `</tbody></table></div>` + kv([["Phase drift / cycle", x.drift_deg.split("|").map((d) => d.replace(":", " ") + "°").join(", ")], ["Search seed", x.seed]]);
    const rdv = r.traj?.recomputed_dv_kms;
    if (rdv != null && Math.abs(rdv - r.dv) > 0.01) {
      html += `<p class="note">The drawn ideal trajectory was re-solved from the catalog's rounded leg times and costs ${fmt(rdv, 3)} km/s per cycle rather than ${fmt(r.dv, 3)}.</p>`;
    }
  }
  if (r.fam !== "home_only") html += boardingHtml(r.ev, r.c3, r.leo, r.fam === "sequence" ? "its lowest-v∞ Earth flyby (ideal model)" : "the Earth flyby (ideal model)");
  return html;
}

function setMode(m) {
  MODE = m;
  setSegmented($("#model"), m);
  $("#model-note").innerHTML = MODE_NOTE[m];
  $("#f-ball-l").textContent = m === "ideal" ? "Ballistic only" : "Has a ballistic window";
  // v∞ drift only exists on real planets
  for (const [sel, fallback] of [["#ax-x", "period"], ["#ax-y", "dv"]]) {
    $(`${sel} option[value=dev]`).disabled = m === "ideal";
    if (m === "ideal" && $(sel).value === "dev") $(sel).value = fallback;
  }
  ideal.st.sort = "dv";
  wireGotos($("#model-note"));
  selectIdeal(ideal.sel);
}

function initIdeal() {
  const dv = $("#f-dv"), vi = $("#f-vinf"), pe = $("#f-per"), c3 = $("#f-c3");
  const sync = () => { dv.value = IF.dvMax; vi.value = IF.vinfMax; pe.value = IF.perMax; c3.value = IF.c3Max; $("#f-q").value = IF.q;
    $("#f-ball").checked = IF.ball; $("#f-earth").checked = IF.earth; setSegmented($("#f-fam"), IF.fam); };
  const chips = () => sysChips($("#f-sys"), IF.sys, renderIdeal);
  sync(); chips();
  segmented($("#model"), setMode);
  segmented($("#f-fam"), (v) => { IF.fam = v; renderIdeal(); });
  const soon = perFrame(renderIdeal);
  dv.oninput = () => { IF.dvMax = +dv.value; soon(); };
  vi.oninput = () => { IF.vinfMax = +vi.value; soon(); };
  pe.oninput = () => { IF.perMax = +pe.value; soon(); };
  c3.oninput = () => { IF.c3Max = +c3.value; soon(); };
  $("#f-q").oninput = (e) => { IF.q = e.target.value; soon(); };
  $("#f-ball").onchange = (e) => { IF.ball = e.target.checked; renderIdeal(); };
  $("#f-earth").onchange = (e) => { IF.earth = e.target.checked; renderIdeal(); };
  $("#f-reset").onclick = () => { Object.assign(IF, IF_DEFAULT, { sys: new Set(SYS) }); sync(); chips(); renderIdeal(); };
  $("#ax-x").onchange = $("#ax-y").onchange = renderIdeal;
  $("#ideal-more").onclick = () => { ideal.st.limit += 300; renderIdeal(); };
  $("#ideal-csv").onclick = () => downloadCSV("cyclers_filtered.csv", idealFiltered().map((r) => ({ ...r.raw, ...(r.real || {}), traj: undefined, real: undefined })));
  $("#o-cycles").onchange = () => selectIdeal(ideal.sel);
  new ResizeObserver(perFrame(() => { if (!$("#tab-ideal").hidden) renderIdeal(); })).observe($("#ideal-scatter"));
}

// ---------------------------------------------------------------------------
// Real-planet windows tab
// ---------------------------------------------------------------------------
const EF_DEFAULT = { label: "", fam: "all", dvMax: 2, vinfMax: 40, c3Max: 150, devMax: 15, y0: 2030, y1: 2060, conv: true, ball: false, earth: false };
const EF = { ...EF_DEFAULT, sys: new Set(SYS) };
const eph = { sel: null, st: { sort: "dv", dir: 1, limit: 200 } };
const eorbit = new OrbitView({ canvas: "#eorbit", empty: "#eorbit-empty", play: "#e-play", slider: "#e-time", clock: "#e-clock", frame: "#e-frame", speed: "#e-speed" });
const EAX = {
  // optimised launch dates can fall a few weeks either side of 2030–2060, so the axis fits the data
  date: { key: "date", label: "Launch date", get: (r) => num(r.year), fmt: (v) => String(Math.round(v)) },
  dv: { key: "dv", label: "ΔV per cycle (km/s)", get: (r) => num(r.dv), log: true },
  worst: { key: "worst", label: "Worst-cycle ΔV (km/s)", get: (r) => num(r.worst), log: true },
  vmax: { key: "vmax", label: "Max v∞ (km/s)", get: (r) => num(r.vmaxNum) },
  c3: { key: "c3", label: "Launch C3 at Earth (km²/s²)", get: (r) => num(r.c3) },
  leo: { key: "leo", label: "ΔV from LEO (km/s)", get: (r) => num(r.leo) },
  dev: { key: "dev", label: "v∞ drift from ideal (km/s)", get: (r) => num(r.dev) },
};
const labelName = (l) => (IDEAL_BY_ID[l]?.notation ? `${l} (${IDEAL_BY_ID[l].notation})` : l);

function ephemFiltered() {
  const within = (v, max, any) => max >= any || (v != null && isFinite(v) && v <= max);
  return EPH.filter((r) => (!EF.label || r.label === EF.label) && (EF.fam === "all" || r.fam === EF.fam) && EF.sys.has(r.sys) &&
    (!EF.conv || (r.conv && r.hasTraj)) && within(r.dv, 10 ** EF.dvMax, 100) && within(r.vmaxNum, EF.vinfMax, 40) &&
    within(r.c3, EF.c3Max, 150) && within(r.dev, EF.devMax, 15) &&
    // the slider ends are open: 2030 means "from the start", 2060 "to the end"
    (r.year == null || ((EF.y0 <= 2030 || r.year >= EF.y0) && (EF.y1 >= 2060 || r.year < EF.y1 + 1))) &&
    (!EF.ball || (r.dv != null && r.dv <= DV_FLOOR)) && (!EF.earth || r.ev != null));
}

const EPH_COLS = [
  { key: "label", label: "Cycler", get: (r) => r.label, html: (r) => `${swatch(r.sys)} ${esc(labelName(r.label))}` },
  { key: "seq", label: "Sequence", get: (r) => r.seq },
  { key: "start", label: "Launch", get: (r) => r.start || "—" },
  { key: "cycles", label: "Cycles", num: true, get: (r) => r.cycles },
  { key: "dv", label: "ΔV / cycle (km/s)", num: true, get: (r) => r.dv, fmt: fmtDv },
  { key: "worst", label: "Worst cycle (km/s)", num: true, get: (r) => r.worst, fmt: (v) => fmt(v, 3) },
  { key: "vmax", label: "Max v∞ (km/s)", get: (r) => (r.vmax || "—").replace(/\|/g, "  ").replace(/:/g, " ") },
  { key: "dev", label: "v∞ drift (km/s)", num: true, get: (r) => r.dev, fmt: (v) => fmt(v, 2) },
  { key: "c3", label: "Launch C3 (km²/s²)", num: true, get: (r) => r.c3, fmt: fmtC3 },
  { key: "leo", label: "ΔV from LEO (km/s)", num: true, get: (r) => r.leo, fmt: (v) => fmt(v, 2) },
  { key: "conv", label: "Converged", get: (r) => (r.conv ? 1 : 0), html: (r) => badge(r.conv, "yes", "no") },
];

function renderEphem() {
  const rows = ephemFiltered();
  const xa = EAX[$("#e-ax-x").value], ya = EAX[$("#e-ax-y").value];
  const title = (a) => a.label + (a.log ? " — log" : "");
  $("#ephem-plot-title").textContent = xa.key === "date" && ya.key === "dv" ? "Launch windows 2030–2060"
    : `${ya.label.replace(/ \(.*\)$/, "")} vs ${xa.label.replace(/ \(.*\)$/, "")}`;
  legend($("#ephem-legend"), SYS.filter((s) => EF.sys.has(s)), false);
  scatter($("#ephem-scatter"), rows.map((r) => ({ id: r.id, x: xa.get(r), y: ya.get(r), sys: r.sys, shape: r.fam === "sequence" ? "square" : "circle", row: r })), {
    xLabel: title(xa), yLabel: title(ya), logX: !!xa.log, logY: !!ya.log, xMin: xa.min, xMax: xa.max, xFmt: xa.fmt,
    selected: eph.sel?.id, aria: "Scatter of real-planet launch windows", onPick: selectEphem,
    tip: (r) => `<div class="t-head">${swatch(r.sys, r.fam === "sequence")}${esc(labelName(r.label))}</div>` + tipRow("Launch", r.start) +
      tipRow("ΔV / cycle", fmtDv(r.dv) + " km/s") + tipRow("Worst cycle", fmt(r.worst, 3) + " km/s") +
      tipRow("Max v∞", esc(r.vmax)) + tipRow("v∞ drift", fmt(r.dev) + " km/s") +
      tipRow("Launch C3", r.c3 == null ? "no Earth encounter" : fmtC3(r.c3) + " km²/s²"),
  });
  eph.st.rerender = renderEphem;
  const n = table($("#ephem-table"), EPH_COLS, rows, eph.st, eph.sel?.id, (id) => selectEphem(EPH_BY_ID[id]));
  $("#ephem-count").textContent = `${fmtInt(n.total)} match · showing ${fmtInt(n.shown)}`;
  $("#ephem-more").hidden = n.shown >= n.total;
  $("#e-dv-v").textContent = EF.dvMax >= 2 ? "any" : `≤ ${+(10 ** EF.dvMax).toPrecision(2)} km/s`;
  $("#e-vinf-v").textContent = EF.vinfMax >= 40 ? "any" : `≤ ${EF.vinfMax} km/s`;
  $("#e-c3-v").textContent = EF.c3Max >= 150 ? "any" : `≤ ${EF.c3Max} km²/s²`;
  $("#e-dev-v").textContent = EF.devMax >= 15 ? "any" : `≤ ${EF.devMax} km/s`;
  $("#e-yr-v").textContent = EF.y0 === 2030 && EF.y1 === 2060 ? "all" : `${EF.y0}–${EF.y1}`;
}

let eToken = 0;
async function selectEphem(r) {
  eph.sel = r;
  const token = ++eToken;
  if (!$("#tab-ephem").hidden) renderEphem();
  $("#eorbit-title").textContent = r ? `${labelName(r.label)} — launch ${r.start}` : "Real-planet trajectory";
  const show = (tr) => {
    $("#edetail").innerHTML = r ? ephemDetail(r, tr) : "";
    const b = $("#edetail .to-ideal");
    if (b) b.onclick = () => { switchTab("ideal"); selectIdeal(IDEAL_BY_ID[r.label]); };
  };
  show(null);
  const tr = await windowTraj(r);
  if (token !== eToken) return;
  show(tr);
  eorbit.set(tr ? ephemModel(tr, r.seq) : null);
}

function ephemDetail(r, tr) {
  const x = r.raw;
  let html = kv([
    ["Cycler", esc(labelName(r.label))], ["Sequence", esc(r.seq)], ["Launch", r.start || "—"], ["Cycles flown", r.cycles],
    ["ΔV / cycle", fmtDv(r.dv) + " km/s"], ["Worst cycle", fmt(r.worst, 3) + " km/s"], ["Total ΔV", fmt(r.total, 3) + " km/s"],
    ["Max v∞", esc((r.vmax || "—").replace(/\|/g, ", ").replace(/:/g, " ")) + " km/s"],
    ["v∞ drift from ideal", fmt(r.dev) + " km/s"], ["SEP accel. (lower bound)", fmt(x.sep_min_accel_mms2, 4) + " mm/s²"],
    ["Converged", badge(r.conv, "yes", "no")],
  ]);
  const dates = String(tr?.dates || "").split("|").filter(Boolean).map((d) => d.split("@"));
  if (!tr && r.hasTraj) html += `<p class="note">Loading encounters…</p>`;
  if (dates.length) {
    html += `<h3>Encounters</h3><div class="table-wrap"><table><thead><tr><th>#</th><th>Planet</th><th>Date</th><th class="num">Days from launch</th><th class="num">v∞ in / out (km/s)</th></tr></thead><tbody>` +
      dates.map(([c, d], i) => {
        const e = tr?.enc?.[i];
        const v = e ? `${fmt(e[2])} / ${fmt(e[3])}` : "";
        return `<tr><td>${i + 1}</td><td>${PNAME[c]}</td><td>${d}</td><td class="num">${Math.round((new Date(d) - new Date(dates[0][1])) / 864e5)}</td><td class="num">${v}</td></tr>`;
      }).join("") + `</tbody></table></div>`;
  }
  html += boardingHtml(r.ev, r.c3, r.leo, "the first Earth encounter of this window");
  if (IDEAL_BY_ID[r.label]) html += `<p><button class="link to-ideal">View the cycler in the catalog →</button></p>`;
  return html;
}

function openEphemFor(label) {
  $("#e-reset").click();   // clear any filters that could hide this cycler's windows
  EF.label = label; $("#e-label").value = label;
  switchTab("ephem");
  const rows = ephemFiltered().sort((a, b) => (a.dv ?? 1e9) - (b.dv ?? 1e9));
  selectEphem(rows[0] || null);
}

function initEphem() {
  const labels = [...new Set(EPH.map((r) => r.label))].sort();
  $("#e-label").innerHTML += labels.map((l) => `<option value="${esc(l)}">${esc(labelName(l))} · ${EPH_BY_LABEL[l].length} windows</option>`).join("");
  const soon = perFrame(renderEphem);
  const sliders = { dv: "dvMax", vinf: "vinfMax", c3: "c3Max", dev: "devMax", y0: "y0", y1: "y1" };
  const sync = () => {
    $("#e-label").value = EF.label;
    for (const [id, k] of Object.entries(sliders)) $(`#e-${id}`).value = EF[k];
    $("#e-conv").checked = EF.conv; $("#e-ball").checked = EF.ball; $("#e-earth").checked = EF.earth;
    setSegmented($("#e-fam"), EF.fam);
    sysChips($("#e-sys"), EF.sys, renderEphem);
  };
  sync();
  $("#e-label").onchange = (e) => { EF.label = e.target.value; renderEphem(); };
  segmented($("#e-fam"), (v) => { EF.fam = v; renderEphem(); });
  for (const [id, k] of Object.entries(sliders)) {
    $(`#e-${id}`).oninput = (e) => {
      EF[k] = +e.target.value;
      // keep the launch-year range ordered
      if (k === "y0" && EF.y0 > EF.y1) { EF.y1 = EF.y0; $("#e-y1").value = EF.y1; }
      if (k === "y1" && EF.y1 < EF.y0) { EF.y0 = EF.y1; $("#e-y0").value = EF.y0; }
      soon();
    };
  }
  for (const [id, k] of [["conv", "conv"], ["ball", "ball"], ["earth", "earth"]]) $(`#e-${id}`).onchange = (e) => { EF[k] = e.target.checked; renderEphem(); };
  $("#e-reset").onclick = () => { Object.assign(EF, EF_DEFAULT, { sys: new Set(SYS) }); sync(); renderEphem(); };
  $("#e-ax-x").onchange = $("#e-ax-y").onchange = renderEphem;
  $("#ephem-more").onclick = () => { eph.st.limit += 300; renderEphem(); };
  $("#ephem-csv").onclick = () => downloadCSV("real_planet_windows_filtered.csv", ephemFiltered().map((r) => ({ ...r.raw })));
  new ResizeObserver(perFrame(() => { if (!$("#tab-ephem").hidden) renderEphem(); })).observe($("#ephem-scatter"));
}

// ---------------------------------------------------------------------------
// Validation tab
// ---------------------------------------------------------------------------
const VD = D.validation;
const vorbit = new OrbitView({ canvas: "#vorbit", empty: "#vorbit-empty", play: "#v-play", slider: "#v-time", clock: "#v-clock", frame: "#v-frame", speed: "#v-speed" });
const valid = { kind: null, key: null };
const FOOT = {
  b: "aphelion below Mars' semimajor axis (1.52 AU)",
  c: "v∞ at Mars = Mars' orbital speed − spacecraft aphelion speed (the arc never reaches Mars' circular orbit)",
  d: "transfer time = \"time to transfer from Earth to aphelion\" (see note)",
  e: "ballistic: required turn angle below the maximum",
};
const po = (p, o, ok, d = 2) => `<span class="po"><span class="p">${p ?? "—"}</span><span>${o == null ? "—" : typeof o === "number" ? fmt(o, d) : o}</span>${ok === false ? '<span class="x" title="outside tolerance">✗</span>' : ""}</span>`;
const passCount = (rows) => `${rows.filter((r) => r.pass).length}/${rows.length}`;

function renderValidation() {
  const t4 = VD.table4, t6 = VD.table6_2002, s1 = VD.s1l1_de405, al = VD.aldrin_real, ch = VD.checks;
  $("#valid-intro").innerHTML = `<h2 style="margin:0 0 6px">Validation against James Longuski's group (Purdue / JPL)</h2>
    <p>Four published results are reproduced, each in the same model the authors used. The two McConaghy, Longuski &amp; Byrnes papers (AIAA 2002-4420 and J. Spacecraft &amp; Rockets 41(4), 2004) state their tables' assumptions explicitly: <b>Earth and Mars on circular, coplanar orbits; synodic period 2 1/7 years (Mars period 1 7/8 years); gravity assists at Earth only, at least 200 km up</b>. Those tables are therefore checked in that idealised model, not with eccentric orbits. The papers' one real-ephemeris result (DE405) is checked on real planets.</p>
    <div class="sumrow">
      <span>${badge(t4.every((r) => r.pass), "pass", "fail")} Table 4 (19 cyclers, every column): <b>${passCount(t4)}</b></span>
      <span>${badge(t6.every((r) => r.pass), "pass", "fail")} 2002 Table 6 (extra Earth flyby): <b>${passCount(t6)}</b></span>
      <span>${badge(s1.rows.every((r) => r.pass), "pass", "fail")} 2004 Table 6 on real ephemerides: <b>${passCount(s1.rows)}</b> encounters</span>
      <span>${badge(ch.every((c) => c.passed), "pass", "fail")} Aldrin cycler (Byrnes, Longuski &amp; Aldrin 1993): <b>${ch.filter((c) => c.passed).length}/${ch.length}</b></span>
    </div>
    <p style="margin-top:8px"><button class="link" data-goto="method#m-valid">How the comparison works →</button> Click any row below to plot it.</p>`;

  // Table 4
  const T = VD.tolerances;
  $("#v-t4").innerHTML = `<div class="card-head"><h2>Table 4 — the most promising cyclers repeating every 1–6 synodic periods</h2><span class="count">paper value, then ours · ✗ = outside tolerance</span></div>
    <div class="table-wrap"><table class="vt" id="t4"><thead><tr><th>Cycler</th><th class="num">Aphelion (AU)</th><th class="num">v∞ Earth (km/s)</th><th class="num">v∞ Mars (km/s)</th><th class="num">Shortest transfer (d)</th><th class="num">Required turn (°)</th><th class="num">Max turn (°)</th><th class="c">Ballistic</th><th>Result</th></tr></thead><tbody>` +
    t4.map((r) => `<tr data-k="${r.notation}" class="${valid.kind === "t4" && valid.key === r.notation ? "sel" : ""}"><td><b>${r.notation}</b>${r.footnotes ? `<sup class="dag">${r.footnotes}</sup>` : ""} <span class="count">${esc(r.note)}</span></td>
      <td class="num">${po(r.paper_aphelion_au.toFixed(2), r.ours_aphelion_au, r.pass_aphelion, 3)}</td>
      <td class="num">${po(r.paper_vinf_earth_kms.toFixed(2), r.ours_vinf_earth_kms, r.pass_vinf_earth)}</td>
      <td class="num">${po(r.paper_vinf_mars_kms.toFixed(2), r.ours_vinf_mars_kms, r.pass_vinf_mars)}</td>
      <td class="num" title="${r.footnotes.includes("d") ? `exact time to aphelion: ${fmt(r.ours_time_to_aphelion_days, 1)} d` : ""}">${po(r.paper_transfer_days, r.ours_transfer_days, r.pass_transfer, 1)}</td>
      <td class="num">${po(r.paper_turn_required_deg, r.ours_turn_required_deg, null, 1)}</td>
      <td class="num">${po(r.paper_turn_max_deg, r.ours_turn_max_deg, r.pass_turn, 1)}</td>
      <td class="c">${po(r.paper_ballistic ? "yes" : "no", r.ours_ballistic ? "yes" : "no", r.pass_ballistic)}</td>
      <td>${badge(r.pass, "pass", "fail")}</td></tr>`).join("") + `</tbody></table></div>
    <ul class="foot-notes">${Object.entries(FOOT).map(([k, v]) => `<li><sup>${k}</sup> ${v}</li>`).join("")}
      <li>Tolerances follow the paper's rounding: aphelion ±${T.aphelion_au} AU, speeds ±${T.vinf_kms} km/s, times ±${T.transfer_days} d, angles ±${T.turn_deg}°.</li>
      <li><b>Note on footnote d.</b> The exact time from Earth to aphelion is 8–10 days longer than each printed value. The time to reach ${(VD.aphelion_fraction * 100).toFixed(1)}% of the aphelion radius reproduces all five within 0.4 days, which is what a radius-based root-finder with a 0.1% tolerance would return (near aphelion the radius hardly changes, so the time moves a lot). That is the value compared here; the exact time to aphelion is in each cell's tooltip. The paper doesn't state this; it is inferred from the numbers.</li></ul>
    <h3 style="margin:14px 0 6px">Scalar checks</h3>
    <div class="table-wrap"><table class="vt"><thead><tr><th>Check</th><th>Cycler</th><th>Source</th><th class="num">Paper</th><th class="num">Ours</th><th class="num">Tolerance</th><th>Result</th></tr></thead><tbody>` +
    ch.map((c) => `<tr data-k="${c.notation}"><td>${esc(c.label)}</td><td>${esc(c.notation)}</td><td>${esc(c.source)}</td><td class="num">${c.paper}</td><td class="num">${c.ours}</td><td class="num">±${c.tolerance}</td><td>${badge(c.passed, "pass", "fail")}</td></tr>`).join("") + `</tbody></table></div>`;

  // 2002 Table 6
  const cr = (v) => (v == null ? "no crossing" : v.map((x) => x.toFixed(3)).join(", "));
  $("#v-t6").innerHTML = `<div class="card-head"><h2>2002 Table 6 — two-synodic-period cyclers with an intermediate Earth flyby</h2><span class="count">idealised model</span></div>
    <p class="note" style="margin:0 0 8px">Each cycler is two Earth-to-Earth legs (P<sub>1</sub>r<sub>1</sub> then P<sub>2</sub>r<sub>2</sub>) that together last 2 synodic periods, with a gravity assist at time τ in between. Because both ends of each leg are at 1 AU, each leg has a single v∞; τ is found where the two legs' v∞ are equal and both flybys can be unpowered. Times are years from the start of the repeat interval.</p>
    <div class="table-wrap"><table class="vt"><thead><tr><th>Cycler</th><th class="num">τ (yr)</th><th class="num">v∞ Earth</th><th>Leg 1 Mars-orbit crossings (yr)</th><th class="num">Leg 1 v∞ Mars</th><th>Leg 2 crossings (yr)</th><th class="num">Leg 2 v∞ Mars</th><th class="num">Turn req./max (°)</th><th>Result</th></tr></thead><tbody>` +
    t6.map((r) => `<tr data-k="${r.cycler}" class="${valid.kind === "t6" && valid.key === r.cycler ? "sel" : ""}"><td><b>${r.cycler}</b></td>
      <td class="num">${po(r.paper_tau_yr.toFixed(3), r.ours_tau_yr, null, 3)}</td>
      <td class="num">${po(r.paper_vinf_earth_kms.toFixed(1), r.ours_vinf_earth_kms)}</td>
      <td><span class="po"><span class="p">${cr(r.paper_leg1_crossings_yr)}</span><span>${cr(r.ours_leg1_crossings_yr)}</span></span></td>
      <td class="num">${po(r.paper_leg1_vinf_mars_kms?.toFixed(1), r.ours_leg1_vinf_mars_kms)}</td>
      <td><span class="po"><span class="p">${cr(r.paper_leg2_crossings_yr)}</span><span>${cr(r.ours_leg2_crossings_yr)}</span></span></td>
      <td class="num">${po(r.paper_leg2_vinf_mars_kms?.toFixed(1), r.ours_leg2_vinf_mars_kms)}</td>
      <td class="num">${fmt(r.turn_required_deg, 1)} / ${fmt(r.turn_max_deg, 1)}</td>
      <td>${badge(r.pass, "pass", "fail")}</td></tr>`).join("") + `</tbody></table></div>`;

  // 2004 Table 6 on real ephemerides
  const S = s1.summary;
  $("#v-s1l1").innerHTML = `<div class="card-head"><h2>2004 Table 6 — 30-year ballistic S1L1 itinerary on real ephemerides</h2><button class="ghost" id="s1l1-plot">Plot the itinerary</button></div>
    <p class="note" style="margin:0 0 8px">The paper flew this cycler on JPL's DE405 ephemeris; this check uses VSOP2013. Lambert arcs connect the real planet positions on the paper's 22 encounter dates. Because the paper gives dates without times, and the ephemerides differ slightly, each date may move by up to 1.5 days to make the chain as ballistic as possible. At the paper's dates exactly, the chain needs ${fmt(S.total_flyby_dv_at_paper_dates_kms, 2)} km/s of flyby ΔV in total; after moving dates by at most ${fmt(S.max_date_shift_days, 2)} days it needs ${fmt(S.total_flyby_dv_kms * 1000, 0)} m/s over 30 years, i.e. essentially ballistic. Largest v∞ difference ${fmt(S.max_vinf_error_kms, 3)} km/s; flyby altitudes agree to a median ${fmt(S.median_altitude_error_pct, 1)}% (largest ${fmt(S.max_altitude_error_pct, 1)}%). Tolerances: v∞ ±0.1 km/s, altitude ±10%.</p>
    <div class="table-wrap"><table class="vt"><thead><tr><th>#</th><th>Body</th><th>Paper date</th><th class="num">Date shift (d)</th><th class="num">Approach v∞ (km/s)</th><th class="num">Closest approach altitude (km)</th><th class="num">Flyby ΔV (km/s)</th><th>Result</th></tr></thead><tbody>` +
    s1.rows.map((r) => `<tr data-k="s1l1"><td>${r.encounter}</td><td>${PNAME[r.body]}</td><td>${r.paper_date}</td><td class="num">${fmt(r.date_shift_days, 2)}</td>
      <td class="num">${po(r.paper_vinf_kms.toFixed(2), r.ours_vinf_kms, r.pass_vinf)}</td>
      <td class="num">${r.paper_altitude_km == null ? (r.encounter === 1 ? "launch" : "—") : po(fmtInt(r.paper_altitude_km), r.ours_altitude_km == null ? null : fmtInt(r.ours_altitude_km), r.pass_altitude)}</td>
      <td class="num">${r.flyby_dv_kms == null ? "" : fmt(r.flyby_dv_kms, 4)}</td><td>${badge(r.pass, "pass", "fail")}</td></tr>`).join("") + `</tbody></table></div>`;

  // Aldrin on real planets (comparison)
  const P = al.paper, rng = (a) => (a ? `${a[0]}–${a[1]}` : "—");
  $("#v-aldrin").innerHTML = `<div class="card-head"><h2>Aldrin cycler on real planets — comparison, not pass/fail</h2></div>
    <p class="note" style="margin:0 0 8px">Byrnes, Longuski &amp; Aldrin (1993), as summarised in McConaghy et al. (2004), optimised the Aldrin cycler over a full 15-year period in an accurate solar-system model. Their dates are not given, and their optimiser was more thorough than this catalog's (whole-trajectory optimisation, not a ±45-day local search with powered flybys only, over 4 cycles per window). So the numbers below are a sanity check of scale, not a reproduction.</p>
    <div class="table-wrap"><table class="vt"><thead><tr><th></th><th class="num">Byrnes et al. (1993)</th><th class="num">This catalog, 1L1, ${al.n_windows} windows 2030–2060</th></tr></thead><tbody>
      <tr><td>ΔV per synodic period (m/s)</td><td class="num">${P.outbound_dv_per_synodic_ms} outbound · ${P.inbound_dv_per_synodic_ms} inbound (15-yr averages)</td><td class="num">median ${al.ours_dv_per_synodic_ms_median ?? "—"} · best ${al.ours_dv_per_synodic_ms_best ?? "—"} · worst ${al.ours_dv_per_synodic_ms_worst ?? "—"}</td></tr>
      <tr><td>Earth → Mars transfer (days)</td><td class="num">${rng(P.transfer_days)}</td><td class="num">${rng(al.ours_transfer_days)}</td></tr>
      <tr><td>v∞ at Earth (km/s)</td><td class="num">${rng(P.vinf_earth_kms)}</td><td class="num">${rng(al.ours_vinf_earth_kms)}</td></tr>
      <tr><td>v∞ at Mars (km/s)</td><td class="num">${rng(P.vinf_mars_kms)}</td><td class="num">${rng(al.ours_vinf_mars_kms)}</td></tr>
    </tbody></table></div>`;

  for (const id of ["#v-t4", "#v-t6", "#v-s1l1"]) {
    $(id).querySelectorAll("tbody tr[data-k]").forEach((tr) => (tr.onclick = () => {
      const k = tr.dataset.k;
      if (k === "s1l1") selectValidation("s1l1", "s1l1");
      else if (id === "#v-t6") selectValidation("t6", k);
      else selectValidation("t4", k);
    }));
  }
  $("#s1l1-plot").onclick = () => selectValidation("s1l1", "s1l1");
  wireGotos($("#valid-intro"));
}

function selectValidation(kind, key) {
  valid.kind = kind; valid.key = key;
  $("#v-cycles-l").hidden = kind === "s1l1";
  let model = null, title = "", html = "";
  if (kind === "t4") {
    const r = VD.table4.find((x) => x.notation === key);
    if (!r) return;
    model = idealModel({ home: "E", traj: r.traj }, Number($("#v-cycles").value));
    title = `${r.notation} — Earth–Mars, the paper's idealised model`;
    const branch = r.branch === "hi" ? "long" : r.branch === "lo" ? "short" : "unique";
    html = kv([
      ["Notation", `${r.notation} (${r.n_synodic} synodic, ${branch}-period, ${r.n_revs} rev)`], ["Period", fmt(r.period_years, 3) + " yr"], ["Result", badge(r.pass, "pass", "fail")],
      ["Aphelion paper / ours", `${r.paper_aphelion_au} / ${fmt(r.ours_aphelion_au, 3)} AU`],
      ["v∞ Earth paper / ours", `${r.paper_vinf_earth_kms} / ${fmt(r.ours_vinf_earth_kms, 3)} km/s`],
      ["v∞ Mars paper / ours", `${r.paper_vinf_mars_kms} / ${fmt(r.ours_vinf_mars_kms, 3)} km/s`],
      ["Transfer paper / ours", `${r.paper_transfer_days} / ${fmt(r.ours_transfer_days, 1)} d`],
      ["Time to aphelion (exact)", fmt(r.ours_time_to_aphelion_days, 1) + " d"],
      ["Turn required paper / ours", `${r.paper_turn_required_deg}° / ${fmt(r.ours_turn_required_deg, 1)}°`],
      ["Max turn paper / ours", `${r.paper_turn_max_deg}° / ${fmt(r.ours_turn_max_deg, 1)}°`],
      ["ΔV per cycle (ours)", fmtDv(r.ours_dv_per_cycle_kms) + " km/s"],
      ["Perihelion / e", `${fmt(r.perihelion_au, 3)} AU / ${fmt(r.e, 3)}`],
    ]) + `<p class="note">Mars v∞: ${esc(r.vinf_mars_definition)}. Transfer: ${esc(r.transfer_definition)}.</p>` +
      (IDEAL_BY_ID[r.catalogue_id] ? `<p><button class="link to-cat">Open the catalog version (mean planetary elements, real planets) →</button></p>` : "");
  } else if (kind === "t6") {
    const r = VD.table6_2002.find((x) => x.cycler === key);
    if (!r || !r.traj) return;
    model = idealModel({ home: "E", traj: r.traj }, Number($("#v-cycles").value));
    title = `${r.cycler} — two-synodic cycler with an Earth flyby at τ = ${fmt(r.ours_tau_yr, 3)} yr`;
    html = kv([["τ paper / ours", `${r.paper_tau_yr} / ${fmt(r.ours_tau_yr, 4)} yr`], ["v∞ Earth paper / ours", `${r.paper_vinf_earth_kms} / ${fmt(r.ours_vinf_earth_kms, 3)} km/s`],
      ["Turn required / max", `${fmt(r.turn_required_deg, 1)}° / ${fmt(r.turn_max_deg, 1)}°`], ["Result", badge(r.pass, "pass", "fail")]]) +
      `<p class="note">Earth flybys are the encounter rings; Mars-orbit crossings are marked with small rings on the path. Mars is placed at the first crossing; the paper does not fix Mars' phase ("we choose the launch date so that the spacecraft encounters Mars").</p>`;
  } else {
    const s = VD.s1l1_de405;
    model = ephemModel(s.traj, "EEM");
    title = "Ballistic S1L1, 2008–2038 on real planets (VSOP2013)";
    html = kv([["Encounters", s.rows.length], ["Total flyby ΔV after date adjustment", fmt(s.summary.total_flyby_dv_kms * 1000, 0) + " m/s over 30 yr"],
      ["Largest date shift", fmt(s.summary.max_date_shift_days, 2) + " d"], ["Largest v∞ difference", fmt(s.summary.max_vinf_error_kms, 3) + " km/s"],
      ["Altitude difference, median / max", `${fmt(s.summary.median_altitude_error_pct, 1)}% / ${fmt(s.summary.max_altitude_error_pct, 1)}%`],
      ["Result", badge(s.rows.every((r) => r.pass), `${s.rows.length}/${s.rows.length} pass`, "fail")]]) +
      `<p class="note">Three legs per 2-synodic-period cycle: Earth → Earth (≈ 1.5 yr), Earth → Mars, Mars → Earth. Encounter dates from McConaghy, Longuski &amp; Byrnes (2004), Table 6.</p>`;
  }
  $("#vorbit-title").textContent = title;
  $("#vdetail").innerHTML = html;
  vorbit.set(model);
  const b = $("#vdetail .to-cat");
  if (b) b.onclick = () => { const r = VD.table4.find((x) => x.notation === key); switchTab("ideal"); selectIdeal(IDEAL_BY_ID[r.catalogue_id]); };
  document.querySelectorAll("#tab-valid tbody tr").forEach((tr) => tr.classList.toggle("sel",
    (kind === "s1l1" && tr.dataset.k === "s1l1") || (kind !== "s1l1" && tr.dataset.k === key && tr.closest(kind === "t6" ? "#v-t6" : "#v-t4"))));
  const v = $("#vviewer").getBoundingClientRect();
  if (v.top < 0 || v.top > innerHeight * 0.6) $("#vviewer").scrollIntoView({ behavior: "smooth", block: "start" });
}

function initValidation() {
  renderValidation();
  $("#v-cycles").onchange = () => valid.kind && selectValidation(valid.kind, valid.key);
}

// ---------------------------------------------------------------------------
// Method tab, tiles, tabs, theme
// ---------------------------------------------------------------------------
function initMethod() {
  const heads = [...document.querySelectorAll("#doc h2, #doc h3")];
  $("#toc").innerHTML = `<b style="display:block;margin:2px 6px 6px">Contents</b>` +
    heads.map((h) => `<a href="#${h.id}" data-id="${h.id}" style="${h.tagName === "H3" ? "padding-left:18px" : ""}">${esc(h.textContent)}</a>`).join("");
  $("#toc").querySelectorAll("a").forEach((a) => (a.onclick = (e) => { e.preventDefault(); document.getElementById(a.dataset.id).scrollIntoView({ behavior: "smooth" }); }));
  wireGotos($("#doc"));
}

// links like data-goto="method#m-real" switch tab and scroll
function wireGotos(root) {
  root.querySelectorAll("[data-goto]").forEach((a) => (a.onclick = (e) => {
    e.preventDefault();
    const [tab, anchor] = a.dataset.goto.split("#");
    switchTab(tab);
    if (anchor) setTimeout(() => document.getElementById(anchor)?.scrollIntoView({ behavior: "smooth" }), 30);
    else window.scrollTo({ top: 0, behavior: "smooth" });
  }));
}

function initTiles() {
  // Each tile counts one kind of thing, and its sub-line adds up to the big number.
  const cyc = IDEAL.filter((r) => r.fam !== "home_only");
  const sh = cyc.filter((r) => r.fam === "single_home"), sq = cyc.filter((r) => r.fam === "sequence");
  const homeOnly = IDEAL.length - cyc.length;
  const conv = EPH.filter((r) => r.conv);
  const withWindows = new Set(EPH.map((r) => r.label)).size;
  const ballWin = conv.filter((r) => r.dv <= DV_FLOOR).length;
  const realBall = cyc.filter((r) => (r.real?.real_ballistic_windows || 0) > 0).length;
  const t4 = VD.table4, t6 = VD.table6_2002, s1 = VD.s1l1_de405.rows, ch = VD.checks;
  const results = [t4.every((r) => r.pass), t6.every((r) => r.pass), s1.every((r) => r.pass), ch.every((c) => c.passed)];
  const tile = (label, value, hint, title) =>
    `<div class="tile" title="${esc(title)}"><div class="label">${label}</div><div class="value">${value}</div><div class="hint">${hint}</div></div>`;
  $("#tiles").innerHTML =
    tile("Cyclers", fmtInt(cyc.length), `${fmtInt(sh.length)} single-assist + ${sq.length} multi-assist`,
      `${homeOnly} more single-assist orbits never reach the target planet; they are listed under Family → Home-only and not counted as cyclers.`) +
    tile("Ballistic cyclers, ideal model", `${cyc.filter((r) => r.ballistic).length} <small>of ${fmtInt(cyc.length)}</small>`,
      `${sh.filter((r) => r.ballistic).length} single-assist + ${sq.filter((r) => r.ballistic).length} multi-assist`,
      "Need ≤ 1 m/s of ΔV per cycle when the planets are on circular, coplanar orbits.") +
    tile("Ballistic cyclers, real planets", `${realBall} <small>of ${fmtInt(withWindows)}</small>`,
      `cyclers with at least one launch window at ≤ 1 m/s per cycle`,
      `Out of the ${withWindows} cyclers that have real-planet launch windows. ${cyc.length - withWindows} three-planet cyclers never line up with the real planets in 2030–2060.`) +
    tile("Real-planet launch windows", fmtInt(EPH.length), `${fmtInt(conv.length)} converged · ${fmtInt(ballWin)} of them ballistic`,
      `Every launch opportunity 2030–2060 for ${withWindows} of the ${cyc.length} cyclers. A window "converges" when the optimiser finds a continuous trajectory.`) +
    tile("Literature validation", `${results.filter(Boolean).length} <small>of 4 reproduced</small>`,
      `Table 4 ${passCount(t4)} · Table 6 ${passCount(t6)} · DE405 ${passCount(s1)} · Aldrin ${ch.filter((c) => c.passed).length}/${ch.length}`,
      "Published results from McConaghy, Longuski & Byrnes (2002, 2004) and Byrnes, Longuski & Aldrin (1993). See the Validation tab.");
  $("#meta-line").textContent = `exported ${D.meta.generated} from ${D.meta.source}`;
}

function switchTab(name) {
  document.querySelectorAll(".tabs button").forEach((b) => b.setAttribute("aria-selected", String(b.dataset.tab === name)));
  for (const t of ["ideal", "ephem", "valid", "method"]) $(`#tab-${t}`).hidden = t !== name;
  $("#tiles").hidden = name === "method";
  history.replaceState(null, "", "#" + name);
  if (name === "ideal") { renderIdeal(); orbit.draw(); }
  if (name === "ephem") { renderEphem(); eorbit.draw(); }
  if (name === "valid") vorbit.draw();
}

function initTheme() {
  const root = document.documentElement;
  try { const t = localStorage.getItem("cycler-theme"); if (t) root.dataset.theme = t; } catch (e) { /* storage unavailable */ }
  const redraw = () => { orbit.draw(); eorbit.draw(); vorbit.draw(); };
  $("#theme").onclick = () => {
    const dark = root.dataset.theme ? root.dataset.theme === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
    root.dataset.theme = dark ? "light" : "dark";
    try { localStorage.setItem("cycler-theme", root.dataset.theme); } catch (e) { /* ignore */ }
    redraw();
  };
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change", redraw);
}

function main() {
  if (!D) { document.body.innerHTML = "<p style='padding:24px'>catalog.js not found. Run <code>python scripts/export_web.py</code>.</p>"; return; }
  initTheme(); initTiles(); initIdeal(); initEphem(); initValidation(); initMethod();
  document.querySelectorAll(".tabs button").forEach((b) => (b.onclick = () => switchTab(b.dataset.tab)));
  wireGotos(document.querySelector("#tab-ephem .intro"));
  const [tab, sub] = location.hash.slice(1).split("/");
  setMode(sub === "real" ? "real" : "ideal");
  selectIdeal(IDEAL_BY_ID["SH-EM-n1r1hi"] || IDEAL[0]);
  const best = EPH.filter((r) => r.conv && r.hasTraj).sort((a, b) => a.dv - b.dv)[0];
  if (best) selectEphem(best);
  selectValidation("t4", "1L1");
  if (["ideal", "ephem", "valid", "method"].includes(tab)) switchTab(tab);
}
main();
