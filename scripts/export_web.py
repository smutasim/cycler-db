# Copyright (c) 2026 Sharif Mutasim
# SPDX-License-Identifier: MIT
"""Export the catalogue (data/cyclers.sqlite) to web/catalog.js for the web UI.

Each cycler gets a compact trajectory: the departure state and flight time of
every Lambert leg (the browser propagates the conics itself), the encounter
times, and the planet phases. Real-ephemeris windows are listed in the main
file; their trajectories go to web/ephem/<cycler>.js, loaded on demand, and
reference a shared table of VSOP2013 planet positions sampled every 4 days.

Units in the JSON: AU, AU/day, days (circular-model trajectories start at t=0;
ephemeris times are days past J2000).

The file is a script (``window.CATALOG = {...}``) rather than bare JSON so that
web/index.html also works when opened straight from disk.

Usage:  python scripts/export_web.py [--db data/cyclers.sqlite] [--out web/catalog.js]
"""
import argparse
import json
import math
import sqlite3
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import cyclerdb  # noqa: E402,F401  (import before pandas)
import numpy as np  # noqa: E402
import pykep as pk  # noqa: E402
from scipy.optimize import minimize  # noqa: E402
from cyclerdb.astro import lambert_all  # noqa: E402
from cyclerdb.bodies import (BODIES, AU, DAY, YEAR, MU_SUN, EphemerisModel, epoch_seconds,  # noqa: E402
                             idealized_bodies, synodic_period)
from cyclerdb.continuation import IdealCycler  # noqa: E402
from cyclerdb.sequence import SequenceProblem, describe  # noqa: E402
from cyclerdb.single_home import enumerate_single_home  # noqa: E402

VEL = DAY / AU  # m/s -> AU/day


def _clean(v):
    """JSON-safe scalar: numpy -> python, NaN/inf -> None."""
    if isinstance(v, (np.floating, float)):
        v = float(v)
        return float(f"{v:.7g}") if math.isfinite(v) else None
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.bool_):
        return bool(v)
    return v


def _rows(con, table):
    try:
        cur = con.execute(f"SELECT * FROM {table}")
    except sqlite3.OperationalError:
        return []
    cols = [c[0] for c in cur.description]
    return [{k: _clean(v) for k, v in zip(cols, r)} for r in cur.fetchall()]


def _leg(t0_days, r, v, dt_days):
    r = np.asarray(r) / AU
    v = np.asarray(v) * VEL
    return [round(t0_days, 4), *[round(float(x), 7) for x in r], *[round(float(x), 9) for x in v], round(dt_days, 4)]


# ---------------------------------------------------------------------------
# Single-assist cyclers
# ---------------------------------------------------------------------------
def sh_traj(c, bodies=BODIES, with_bodies=False):
    """Trajectory of one single-assist cycler (circular model, home at longitude 0 at t=0)."""
    h, t = c.home, c.target
    H = bodies[h]
    Tc = c.period_years * YEAR
    r0 = np.array([H.a, 0.0, 0.0])
    vH = np.array([0.0, H.v_circ, 0.0])
    th = 2 * np.pi * Tc / H.period
    r1 = H.a * np.array([np.cos(th), np.sin(th), 0.0])
    best = None
    for n_rev, _b, v0, _v1 in lambert_all(r0, r1, Tc, c.n_revs):
        if n_rev != c.n_revs:
            continue
        d = abs(np.linalg.norm(v0 - vH) / 1e3 - c.vinf_home_kms)
        if best is None or d < best[0]:
            best = (d, v0)
    if best is None:
        return None
    v0 = best[1]
    enc = [[h, 0.0]]
    reaches = c.transfer_days is not None and c.reaches_target
    if reaches:
        tf = c.transfer_days * DAY
        enc.append([t, round(c.transfer_days, 3)])
    else:
        # never meets the target: place it at the arc's closest approach to its orbit
        P = 2 * np.pi * np.sqrt((c.a_au * AU) ** 3 / MU_SUN)
        tf = c.aphelion_days * DAY if bodies[t].a > H.a else (c.aphelion_days * DAY + P / 2) % P
    rT, _ = pk.propagate_lagrangian([list(r0), list(v0)], tf, MU_SUN)
    phase_t = (math.atan2(rT[1], rT[0]) - 2 * np.pi * tf / bodies[t].period) % (2 * np.pi)
    out = dict(period_days=round(Tc / DAY, 4), rotate_deg=round(math.degrees(th) % 360, 4),
               legs=[_leg(0.0, r0, v0, Tc / DAY)], enc=enc,
               phases={h: 0.0, t: round(math.degrees(phase_t) % 360, 4)})
    if with_bodies:
        out["bodies"] = {k: dict(a_au=b.a_au, period_days=b.period / DAY) for k, b in bodies.items()}
    return out


def single_home_traj(rows):
    cache = {}
    out = {}
    for r in rows:
        h, t = r["home"], r["target"]
        if (h, t) not in cache:
            S = synodic_period(BODIES[h], BODIES[t])
            nmax = max(1, int(16 * YEAR // S))
            cache[(h, t)] = {(c.n_synodic, c.n_revs, c.branch): c
                             for c in enumerate_single_home(h, t, max_synodic=nmax, max_revs=12)}
        c = cache[(h, t)].get((r["n_synodic"], r["n_revs"], r["branch"]))
        if c is not None:
            out[r["cycler_id"]] = sh_traj(c)
    return out


# ---------------------------------------------------------------------------
# Multi-assist sequence cyclers
# ---------------------------------------------------------------------------
def _seq_problem(row):
    seq, Tc = row["sequence"], row["period_years"] * YEAR
    revs = [int(v) for v in str(row["leg_revs"]).split("|")]
    max_revs = max(max(revs), int(min(5, np.ceil(Tc / YEAR / len(seq)) + 1)))
    prob = SequenceProblem(seq, Tc, max_revs=max_revs)
    legs = [float(v) for v in str(row["leg_days"]).split("|")]
    ph = {kv.split(":")[0]: math.radians(float(kv.split(":")[1])) for kv in str(row["phases_deg"]).split("|")}
    x = [d * DAY / Tc for d in legs] + [ph[c] for c in prob.others]
    return prob, x, revs


def sequence_traj(rows):
    out = {}
    for r in rows:
        prob, x, revs = _seq_problem(r)
        # leg_days/phases are stored rounded; a short local polish recovers the catalogued cost
        # (describe() gives the pure flyby ΔV; evaluate() adds the soft v∞-cap penalty)
        c0 = describe(prob, x)
        if c0 is None:
            continue
        dv = c0.dv_per_cycle_kms
        if dv > r["dv_per_cycle_kms"] + 1e-3:
            loc = minimize(prob.evaluate, x, method="Nelder-Mead",
                           options=dict(maxiter=600, xatol=1e-10, fatol=1e-6))
            c = describe(prob, loc.x)
            if c is not None and c.dv_per_cycle_kms < dv:
                x, dv = list(loc.x), c.dv_per_cycle_kms
        _, det = prob.evaluate(x, detail=True)
        legs, path, t, phases = det
        out_legs, enc = [], []
        for i in range(prob.m):
            r0, _ = prob._state(prob.seq[i], t[i], phases)
            v0 = legs[i][path[i]][3]
            out_legs.append(_leg(t[i] / DAY, r0, v0, (t[i + 1] - t[i]) / DAY))
            enc.append([prob.seq[i], round(t[i] / DAY, 3)])
        out[r["cycler_id"]] = dict(
            period_days=round(prob.Tc / DAY, 4), rotate_deg=round(math.degrees(prob.dth) % 360, 4),
            legs=out_legs, enc=enc, recomputed_dv_kms=round(dv, 4),
            phases={k: round(math.degrees(v) % 360, 4) for k, v in phases.items()})
    return out


# ---------------------------------------------------------------------------
# Ephemeris windows
# ---------------------------------------------------------------------------
def ephemeris_traj(rows, sh_rows, sq_rows, model):
    """Trajectories of the real-ephemeris windows, rebuilt from the encounter dates.

    Each leg is the Lambert arc between the real planet states on consecutive
    encounter dates (dates are stored to the day, so the drawn path can differ
    slightly from the optimised one). Encounters carry the V∞ in and out.
    """
    sh_by = {r["cycler_id"]: r for r in sh_rows}
    sq_by = {r["cycler_id"]: r for r in sq_rows}
    ideals, out = {}, []
    for i, r in enumerate(rows):
        r["window_id"] = f"{r['label']}@{r['start_epoch']}" if r.get("start_epoch") else f"{r['label']}#{i}"
        dates = [d.split("@") for d in str(r.get("encounter_dates") or "").split("|") if "@" in d]
        if len(dates) < 2:
            out.append(None)
            continue
        if r["label"] not in ideals:
            ideals[r["label"]] = _ideal_for_label(r["label"], sh_by, sq_by)
        ic = ideals[r["label"]]
        if ic is None:
            out.append(None)
            continue
        m = len(ic.seq)
        codes = [c for c, _ in dates]
        ts = [epoch_seconds(d) for _, d in dates]
        legs, vin, vout = [], [None] * len(ts), [None] * len(ts)
        for k in range(len(ts) - 1):
            r0, vp = model.state(codes[k], ts[k])
            r1, vq = model.state(codes[k + 1], ts[k + 1])
            want = ic.leg_revs[k % m]
            vi = ic.vinf_kms[k % m] * 1e3 if ic.vinf_kms else None
            sols = [(abs(np.linalg.norm(v0 - vp) - vi) if vi else 0.0, v0, v1)
                    for n, _b, v0, v1 in lambert_all(r0, r1, ts[k + 1] - ts[k], want)
                    if n == want and np.isfinite(v0).all()]
            if not sols:
                legs = None
                break
            _, v0, v1 = min(sols, key=lambda s: s[0])
            vout[k] = round(float(np.linalg.norm(v0 - vp)) / 1e3, 3)
            vin[k + 1] = round(float(np.linalg.norm(v1 - vq)) / 1e3, 3)
            legs.append(_leg(ts[k] / DAY, r0, v0, (ts[k + 1] - ts[k]) / DAY))
        out.append(None if legs is None else dict(
            legs=legs, enc=[[c, round(t / DAY, 3), a, b] for c, t, a, b in zip(codes, ts, vin, vout)]))
    return out


def _ideal_for_label(label, sh_by, sq_by):
    if label.startswith("SQ-"):
        row = sq_by.get(label)
        return None if row is None else IdealCycler.from_sequence_row(row)
    row = sh_by.get(label)
    if row is None:
        return None
    c = next((c for c in enumerate_single_home(row["home"], row["target"], row["n_synodic"], 12)
              if c.n_synodic == row["n_synodic"] and c.n_revs == row["n_revs"] and c.branch == row["branch"]), None)
    return None if c is None else IdealCycler.from_single_home(c)


def planet_table(model, start="2005-01-01", end="2100-01-01", step_days=4.0):
    t0, t1 = epoch_seconds(start), epoch_seconds(end)
    ts = np.arange(t0, t1, step_days * DAY)
    tab = dict(t0_days=round(t0 / DAY, 4), step_days=step_days, n=len(ts))
    for c in "VEM":
        xs = []
        for t in ts:
            rr, _ = model.state(c, t)
            xs += [round(float(v) / AU, 5) for v in rr]
        tab[c] = xs
    return tab


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
def _bodies_json(bodies):
    return {k: dict(a_au=b.a_au, period_days=b.period / DAY) for k, b in bodies.items()}


def validation_export(real_windows):
    """All literature reproductions, with trajectories for the UI."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from validate_mcconaghy import (run, run_table6_2002, run_s1l1_de405, CHECK_INFO, TOL,
                                    APHELION_FRACTION)
    from cyclerdb.validation import ALDRIN_REAL_B93
    ib = idealized_bodies()
    E, M = ib["E"], ib["M"]

    # Table 4 (circular coplanar, idealised planets)
    t4, checks = run()
    sols = {c.notation: c for c in enumerate_single_home("E", "M", 6, 10, bodies=ib)}
    table4 = []
    for r in t4.to_dict("records"):
        r = {k: _clean(v) for k, v in r.items()}
        c = sols[r["notation"]]
        r["traj"] = sh_traj(c, ib, with_bodies=True)
        r["catalogue_id"] = f"SH-EM-n{c.n_synodic}r{c.n_revs}{c.branch}"
        table4.append(r)

    # [M02] Table 6: two Earth-to-Earth legs, Mars phased at the first crossing
    t6 = run_table6_2002(with_traj=True)
    table6 = []
    for r in t6.to_dict("records"):
        legs = r.pop("_legs", None)
        r.pop("_bodies", None)
        T = r.pop("_T", None)
        r = {k: (list(v) if isinstance(v, tuple) else _clean(v)) for k, v in r.items()}
        if legs:
            crossings = (r.get("ours_leg1_crossings_yr") or []) + (r.get("ours_leg2_crossings_yr") or [])
            phases = {"E": 0.0, "M": 0.0}
            if crossings:
                tc = crossings[0] * YEAR
                L = legs[0] if tc < legs[1][0] else legs[1]
                rr, _ = pk.propagate_lagrangian([list(L[1]), list(L[2])], tc - L[0], MU_SUN)
                phases["M"] = math.degrees(math.atan2(rr[1], rr[0]) - 2 * np.pi * tc / M.period) % 360
            r["traj"] = dict(period_days=T / DAY, rotate_deg=math.degrees(2 * np.pi * T / E.period) % 360,
                             legs=[_leg(t0 / DAY, r0, v0, dt / DAY) for t0, r0, v0, dt in legs],
                             enc=[["E", 0.0], ["E", round(legs[1][0] / DAY, 3)]],
                             phases=phases, bodies=_bodies_json(ib),
                             marks=[["M", round(x * YEAR / DAY, 3)] for x in crossings])
        table6.append(r)

    # [M04] Table 6: the 30-year ballistic S1L1 itinerary on real ephemerides
    s1, summ = run_s1l1_de405(with_traj=True)
    legs = summ.pop("_legs")
    enc = summ.pop("_enc")
    s1l1 = dict(rows=[{k: _clean(v) for k, v in r.items()} for r in s1.to_dict("records")],
                summary={k: _clean(v) for k, v in summ.items()},
                traj=dict(legs=[_leg(t0 / DAY, r0, v0, dt / DAY) for t0, r0, v0, dt in legs],
                          enc=[[c, round(t / DAY, 3)] for c, t in enc]))

    # [B93] Aldrin cycler on real planets vs our continuation of 1L1 (comparison, not pass/fail)
    al = [w for w in real_windows if w["label"] == "SH-EM-n1r1hi" and w.get("converged")]
    aldrin = dict(paper=ALDRIN_REAL_B93, n_windows=len(al))
    if al:
        dv = np.array([w["dv_per_cycle_kms"] for w in al])
        tofs, vE, vM = [], [], []
        for w in al:
            tr = w.get("traj")
            if not tr:
                continue
            e = tr["enc"]
            for k in range(len(e) - 1):
                if e[k][0] == "E" and e[k + 1][0] == "M":
                    tofs.append(e[k + 1][1] - e[k][1])
            vE += [x for c, _, a, b in e if c == "E" for x in (a, b) if x is not None]
            vM += [x for c, _, a, b in e if c == "M" for x in (a, b) if x is not None]
        aldrin.update(ours_dv_per_synodic_ms_median=round(float(np.median(dv)) * 1e3),
                      ours_dv_per_synodic_ms_best=round(float(dv.min()) * 1e3),
                      ours_dv_per_synodic_ms_worst=round(float(dv.max()) * 1e3),
                      ours_transfer_days=[round(min(tofs)), round(max(tofs))] if tofs else None,
                      ours_vinf_earth_kms=[round(min(vE), 2), round(max(vE), 2)] if vE else None,
                      ours_vinf_mars_kms=[round(min(vM), 2), round(max(vM), 2)] if vM else None)
    scalar = [dict(key=k, label=CHECK_INFO[k][0], source=CHECK_INFO[k][1], notation=CHECK_INFO[k][2],
                   ours=round(float(ours), 3), paper=ref, tolerance=tol, passed=bool(abs(ours - ref) <= tol))
              for k, (ours, ref, tol) in checks.items()]
    return dict(table4=table4, table6_2002=table6, s1l1_de405=s1l1, aldrin_real=aldrin, checks=scalar,
                tolerances=TOL, aphelion_fraction=APHELION_FRACTION)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def _safe(label):
    return "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in label)


# Window table, stored column-wise to keep catalog.js small. Encounter dates live in
# the per-cycler trajectory files (web/ephem/<label>.js), loaded on demand.
WINDOW_COLS = ["window_id", "label", "sequence", "start_epoch", "n_cycles", "dv_per_cycle_kms",
               "worst_cycle_dv_kms", "total_dv_kms", "max_vinf_kms", "converged", "max_vinf_deviation_kms",
               "sep_min_accel_mms2", "earth_vinf_kms", "earth_c3_km2s2", "leo_departure_dv_kms", "has_traj"]


def _short(v):
    return float(f"{v:.5g}") if isinstance(v, float) else v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=str(ROOT / "data" / "cyclers.sqlite"))
    ap.add_argument("--out", default=str(ROOT / "web" / "catalog.js"))
    a = ap.parse_args()
    t_start = time.time()
    with sqlite3.connect(a.db) as con:
        sh = _rows(con, "single_home")
        sq = _rows(con, "sequence")
        ep = _rows(con, "ephemeris")
        real = {r["label"]: r for r in _rows(con, "real_summary")}
    model = EphemerisModel()
    sh_traj_by = single_home_traj(sh)
    for r in sh:
        r["traj"] = sh_traj_by.get(r["cycler_id"])
        r["real"] = real.get(r["cycler_id"])
    sq_traj_by = sequence_traj(sq)
    for r in sq:
        r["traj"] = sq_traj_by.get(r["cycler_id"])
        r["real"] = real.get(r["cycler_id"])
    print(f"ideal trajectories done ({time.time() - t_start:.0f}s)", flush=True)

    # real-ephemeris windows: compact rows in the main file, trajectories per cycler on demand
    for r, tr in zip(ep, ephemeris_traj(ep, sh, sq, model)):
        r["traj"] = tr
    print(f"ephemeris trajectories done ({time.time() - t_start:.0f}s)", flush=True)
    out = Path(a.out)
    edir = out.parent / "ephem"
    edir.mkdir(parents=True, exist_ok=True)
    by_label = {}
    for r in ep:
        tr = r["traj"]
        if tr is not None:
            tr["dates"] = r.get("encounter_dates")
        by_label.setdefault(r["label"], {})[r["window_id"]] = tr
    for label, trajs in by_label.items():
        (edir / f"{_safe(label)}.js").write_text(
            f"window.__ephemLoaded({json.dumps(label)}, " + json.dumps(trajs, separators=(",", ":"), allow_nan=False) + ");\n",
            encoding="utf-8")
    for r in ep:
        r["has_traj"] = r["traj"] is not None
    windows = dict(cols=WINDOW_COLS, rows=[[_short(r.get(k)) for k in WINDOW_COLS] for r in ep])

    validation = validation_export(ep)
    print(f"validation done ({time.time() - t_start:.0f}s)", flush=True)
    bodies = {c: dict(name=b.name, a_au=b.a_au, period_days=b.period / DAY, min_alt_km=b.min_alt / 1e3)
              for c, b in BODIES.items()}
    doc = dict(
        meta=dict(generated=time.strftime("%Y-%m-%d %H:%M"), source=str(Path(a.db).name),
                  mu_sun_au3_day2=MU_SUN / AU**3 * DAY**2, j2000_iso="2000-01-01T12:00:00",
                  bodies=bodies, ephemeris="VSOP2013 (J2000 ecliptic)"),
        single_home=sh, sequence=sq, ephemeris=windows, validation=validation,
        planets=planet_table(model))
    out.write_text("window.CATALOG = " + json.dumps(doc, separators=(",", ":"), allow_nan=False) + ";\n",
                   encoding="utf-8")
    n_missing = sum(r["traj"] is None for r in sh + sq)
    print(f"wrote {out} ({out.stat().st_size / 1e6:.2f} MB) and {len(by_label)} files in {edir}: "
          f"{len(sh)} single-assist, {len(sq)} sequence, {len(ep)} real-planet windows; "
          f"{n_missing} ideal cyclers without trajectory ({time.time() - t_start:.0f}s)")
    sys.stdout.flush()
    import os; os._exit(0)  # skip pykep teardown crash (see README known issues)


if __name__ == "__main__":
    main()
