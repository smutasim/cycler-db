# Copyright (c) 2026 Sharif Mutasim
# SPDX-License-Identifier: MIT
"""Reproduce the published cycler tables of McConaghy, Longuski & Byrnes.

  run()                -> [M02]/[M04] Table 4 (19 cyclers, circular coplanar model)
  run_table6_2002()    -> [M02] Table 6 (two-synodic cyclers with an intermediate
                          Earth flyby, circular coplanar model)
  run_s1l1_de405()     -> [M04] Table 6 (30-year ballistic S1L1 itinerary on real
                          ephemerides; the paper used DE405, we use VSOP2013)

Sources and the paper's model are described in cyclerdb/validation.py.

Usage:  python scripts/validate_mcconaghy.py  [--out data/validation_mcconaghy.csv]
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cyclerdb  # noqa: E402,F401  (must precede pandas, see cyclerdb/__init__.py)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import brentq, minimize  # noqa: E402
from cyclerdb.astro import (angle_between, conic, flyby_dv, lambert_all, max_turn,  # noqa: E402
                            radius_crossings, rot_z)
from cyclerdb.bodies import BODIES, DAY, MU_SUN, YEAR, EphemerisModel, epoch_seconds, idealized_bodies  # noqa: E402
from cyclerdb.single_home import departure_state, enumerate_single_home  # noqa: E402
from cyclerdb.validation import (ALDRIN_1993, NOTABLE_2SYN, S1L1_DE405, TABLE4, TABLE4_NOTES,  # noqa: E402
                                 TABLE6_2002)

# Tolerances follow the paper's rounding: aphelion to 0.01 AU, speeds to 0.01 km/s,
# times to whole days, angles to whole degrees (each with a little margin).
TOL = dict(aphelion_au=0.01, vinf_kms=0.02, transfer_days=2.0, turn_deg=1.0)

# Footnote d rows: the paper says "time to transfer from Earth to aphelion", but
# its numbers are 8-10 days shorter than the exact time to aphelion. All five
# are reproduced (within 0.4 d) by the time to reach 99.9% of the aphelion
# radius, which is what a radius-based root find with a 0.1% tolerance returns.
# This is a reconstruction from the numbers, not something the paper states.
APHELION_FRACTION = 0.999

# label, source and cycler (notation) of each scalar check in run()
CHECK_INFO = {
    "aldrin_vinf_earth": ("Aldrin cycler: v∞ at Earth (km/s)", "Byrnes, Longuski & Aldrin (1993)", "1L1"),
    "aldrin_vinf_mars": ("Aldrin cycler: v∞ at Mars (km/s)", "Byrnes, Longuski & Aldrin (1993)", "1L1"),
    "aldrin_transfer_days": ("Aldrin cycler: Earth→Mars transfer (days)", "Byrnes, Longuski & Aldrin (1993)", "1L1"),
    "notable2syn_vinf_earth": ("Two-synodic cycler: v∞ at Earth (km/s)", "McConaghy, Longuski & Byrnes (2002)", "2L3"),
}


# ---------------------------------------------------------------------------
# Table 4
# ---------------------------------------------------------------------------
def run(out_csv=None):
    ib = idealized_bodies()
    M = ib["M"]
    sols = {c.notation: c for c in enumerate_single_home("E", "M", max_synodic=6, max_revs=10, bodies=ib)}
    rows = []
    for nota, aph, vE, vM, tof, req, mx, flags in TABLE4:
        c = sols[nota]
        r0, v0 = departure_state(c, ib)
        a, e, rp, ra, P = conic(r0, v0)
        if "c" in flags:   # footnote c: Mars' speed minus spacecraft aphelion speed
            v_aph = np.sqrt(MU_SUN * (2 / ra - 1 / a))
            vM_ours = (M.v_circ - v_aph) / 1e3
        else:
            vM_ours = c.vinf_target_kms
        t_aph = c.aphelion_days
        if "d" in flags:   # footnote d: "time to transfer from Earth to aphelion"
            tof_ours = radius_crossings(r0, v0, APHELION_FRACTION * ra)["t_out"] / DAY
        else:
            tof_ours = c.transfer_days
        ballistic = c.dv_per_cycle_kms <= 1e-3
        d = dict(
            notation=nota, n_synodic=c.n_synodic, n_revs=c.n_revs, branch=c.branch, footnotes=flags,
            paper_aphelion_au=aph, ours_aphelion_au=round(c.aphelion_au, 3),
            paper_vinf_earth_kms=vE, ours_vinf_earth_kms=round(c.vinf_home_kms, 3),
            paper_vinf_mars_kms=vM, ours_vinf_mars_kms=round(vM_ours, 3),
            vinf_mars_definition="Mars speed − aphelion speed (footnote c)" if "c" in flags else "v∞ at the Mars-orbit crossing",
            paper_transfer_days=tof, ours_transfer_days=round(tof_ours, 1),
            transfer_definition=(f"time to {APHELION_FRACTION:.1%} of aphelion radius (footnote d, reconstructed)"
                                 if "d" in flags else "Earth → first Mars-orbit crossing"),
            ours_time_to_aphelion_days=round(t_aph, 1),
            paper_turn_required_deg=req, ours_turn_required_deg=round(c.turn_required_deg, 1),
            paper_turn_max_deg=mx, ours_turn_max_deg=round(c.turn_max_deg, 1),
            paper_ballistic="e" in flags, ours_ballistic=bool(ballistic),
            ours_dv_per_cycle_kms=round(c.dv_per_cycle_kms, 3),
            period_years=c.period_years, a_au=c.a_au, e=c.e, perihelion_au=c.perihelion_au,
            note=TABLE4_NOTES.get(nota, ""),
        )
        d["pass_aphelion"] = abs(d["ours_aphelion_au"] - aph) <= TOL["aphelion_au"]
        d["pass_vinf_earth"] = abs(d["ours_vinf_earth_kms"] - vE) <= TOL["vinf_kms"]
        d["pass_vinf_mars"] = abs(d["ours_vinf_mars_kms"] - vM) <= TOL["vinf_kms"]
        d["pass_transfer"] = abs(d["ours_transfer_days"] - tof) <= TOL["transfer_days"]
        d["pass_turn"] = (abs(d["ours_turn_required_deg"] - req) <= TOL["turn_deg"]
                          and abs(d["ours_turn_max_deg"] - mx) <= TOL["turn_deg"])
        d["pass_ballistic"] = d["paper_ballistic"] == d["ours_ballistic"]
        d["pass"] = all(d[k] for k in d if k.startswith("pass_"))
        rows.append(d)
    df = pd.DataFrame(rows)

    aldrin, notable = sols["1L1"], sols["2L3"]
    checks = {
        "aldrin_vinf_earth": (aldrin.vinf_home_kms, ALDRIN_1993["vinf_earth_kms"], 0.02),
        "aldrin_vinf_mars": (aldrin.vinf_target_kms, ALDRIN_1993["vinf_mars_kms"], 0.02),
        "aldrin_transfer_days": (aldrin.transfer_days, ALDRIN_1993["transfer_days"], 1.0),
        "notable2syn_vinf_earth": (notable.vinf_home_kms, NOTABLE_2SYN["vinf_earth_kms"], 0.02),
    }
    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)
    return df, checks


# ---------------------------------------------------------------------------
# [M02] Table 6: two-synodic cyclers with an intermediate Earth flyby
# ---------------------------------------------------------------------------
def _earth_leg(E, t0, t1, P, r):
    th0, th1 = 2 * np.pi * t0 / E.period, 2 * np.pi * t1 / E.period
    r0 = E.a * np.array([np.cos(th0), np.sin(th0), 0.0])
    r1 = E.a * np.array([np.cos(th1), np.sin(th1), 0.0])
    vE0 = E.v_circ * np.array([-np.sin(th0), np.cos(th0), 0.0])
    vE1 = E.v_circ * np.array([-np.sin(th1), np.cos(th1), 0.0])
    sols = [(conic(r0, v0)[0], v0, v1) for n, _, v0, v1 in lambert_all(r0, r1, t1 - t0, r)
            if n == r and np.isfinite(v0).all()]
    sols = sorted((s for s in sols if s[0] > 0), key=lambda s: -s[0])   # L (larger a) first
    if P == "U" and sols:
        s = sols[0]
    elif P in "LS" and len(sols) == 2:
        s = sols[0] if P == "L" else sols[1]
    else:
        return None
    return dict(r0=r0, v0=s[1], vi0=s[1] - vE0, vi1=s[2] - vE1, t0=t0, t1=t1)


def _mars_crossings(L, M):
    x = radius_crossings(L["r0"], L["v0"], M.a)
    if x is None:
        return None, None
    dur, out = L["t1"] - L["t0"], []
    for base in (x["t_out"], x["t_in"]):
        k = 0
        while base + k * x["period"] < dur:
            out.append((L["t0"] + base + k * x["period"]) / YEAR)
            k += 1
    s, g = x["speed"], x["gamma"]
    vM = np.sqrt(s**2 + M.v_circ**2 - 2 * s * M.v_circ * np.cos(g)) / 1e3
    return sorted(out), vM


def run_table6_2002(with_traj=False):
    """Find the intermediate flyby time tau at which both Earth-to-Earth legs have the
    same V∞ (each leg has equal V∞ at both ends, since both ends are at 1 AU), keep the
    roots where both flybys can be made unpowered, and compare with the paper."""
    ib = idealized_bodies()
    E, M = ib["E"], ib["M"]
    T = 2 / (1 / E.period - 1 / M.period)
    rot = rot_z(2 * np.pi * T / E.period)
    rows = []
    for P1r1, P2r2, tau_p, vE_p, c1_p, vm1_p, c2_p, vm2_p in TABLE6_2002:
        P1, r1, P2, r2 = P1r1[0], int(P1r1[1:]), P2r2[0], int(P2r2[1:])

        def f(tau):
            a, b = _earth_leg(E, 0, tau, P1, r1), _earth_leg(E, tau, T, P2, r2)
            return np.nan if a is None or b is None else np.linalg.norm(a["vi0"]) - np.linalg.norm(b["vi0"])
        grid = np.linspace(T / 2, T, 1500)[1:-1]
        vals = np.array([f(t) for t in grid])
        best = None
        for i in range(len(grid) - 1):
            if not (np.isfinite(vals[i]) and np.isfinite(vals[i + 1])) or np.sign(vals[i]) == np.sign(vals[i + 1]):
                continue
            if abs(vals[i] - vals[i + 1]) > 2000:
                continue
            tau = brentq(f, grid[i], grid[i + 1], xtol=1.0)
            L1, L2 = _earth_leg(E, 0, tau, P1, r1), _earth_leg(E, tau, T, P2, r2)
            v = np.linalg.norm(L1["vi0"])
            turns = (angle_between(L1["vi1"], L2["vi0"]), angle_between(L2["vi1"], rot @ L1["vi0"]))
            if max(turns) > max_turn(v, E):
                continue  # not ballistic: the paper lists ballistic cyclers only
            if best is None or abs(tau / YEAR - tau_p) < abs(best[0] / YEAR - tau_p):
                best = (tau, L1, L2, v, turns)
        d = dict(cycler=P1r1 + P2r2, paper_tau_yr=tau_p, paper_vinf_earth_kms=vE_p,
                 paper_leg1_crossings_yr=c1_p, paper_leg1_vinf_mars_kms=vm1_p,
                 paper_leg2_crossings_yr=c2_p, paper_leg2_vinf_mars_kms=vm2_p)
        if best is None:
            d.update(found=False, **{"pass": False})
            rows.append(d)
            continue
        tau, L1, L2, v, turns = best
        c1, m1 = _mars_crossings(L1, M)
        c2, m2 = _mars_crossings(L2, M)
        d.update(found=True, ours_tau_yr=round(tau / YEAR, 4), ours_vinf_earth_kms=round(v / 1e3, 3),
                 ours_leg1_crossings_yr=None if c1 is None else tuple(round(x, 4) for x in c1),
                 ours_leg1_vinf_mars_kms=None if m1 is None else round(m1, 3),
                 ours_leg2_crossings_yr=None if c2 is None else tuple(round(x, 4) for x in c2),
                 ours_leg2_vinf_mars_kms=None if m2 is None else round(m2, 3),
                 turn_required_deg=round(np.degrees(max(turns)), 1), turn_max_deg=round(np.degrees(max_turn(v, E)), 1))

        def same(p, o, tol):
            if p is None or o is None:
                return p is None and o is None
            if isinstance(p, tuple):
                return len(p) == len(o) and all(abs(a - b) <= tol for a, b in zip(p, o))
            return abs(p - o) <= tol
        # paper rounds times to 0.001 yr and speeds to 0.1 km/s
        d["pass"] = bool(same(tau_p, d["ours_tau_yr"], 0.002) and same(vE_p, d["ours_vinf_earth_kms"], 0.05)
                         and same(c1_p, d["ours_leg1_crossings_yr"], 0.002) and same(vm1_p, d["ours_leg1_vinf_mars_kms"], 0.05)
                         and same(c2_p, d["ours_leg2_crossings_yr"], 0.002) and same(vm2_p, d["ours_leg2_vinf_mars_kms"], 0.05))
        if with_traj:
            d["_legs"] = [(L1["t0"], L1["r0"], L1["v0"], L1["t1"] - L1["t0"]), (L2["t0"], L2["r0"], L2["v0"], L2["t1"] - L2["t0"])]
            d["_bodies"] = ib
            d["_T"] = T
        rows.append(d)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# [M04] Table 6: ballistic S1L1 itinerary on real ephemerides
# ---------------------------------------------------------------------------
def run_s1l1_de405(shift_days=1.5, with_traj=False):
    """Fly Lambert arcs between real planet states (VSOP2013) on the paper's dates.

    The paper gives dates only (no times) and used DE405, so each encounter epoch
    may move by up to ``shift_days`` to minimise the total flyby ΔV. Revolutions
    per leg and the Lambert branches are chosen by a dynamic program, exactly as
    in the catalogue's ephemeris continuation. Returns (per-encounter table, summary).
    """
    model = EphemerisModel()
    codes = [b for b, *_ in S1L1_DE405]
    t_paper = np.array([epoch_seconds(d) for _, d, *_ in S1L1_DE405])
    cache = {}

    def state(c, t):
        k = (c, round(t, 0))
        if k not in cache:
            cache[k] = model.state(c, t)
        return cache[k]

    def chain(t, revs=None):
        L = []
        for k in range(len(t) - 1):
            r0, vp = state(codes[k], t[k])
            r1, vq = state(codes[k + 1], t[k + 1])
            sols = [(n, v0 - vp, v1 - vq, v0) for n, _, v0, v1 in lambert_all(r0, r1, t[k + 1] - t[k], 3)
                    if np.isfinite(v0).all() and (revs is None or n == revs[k])]
            if not sols:
                return 1e9, None
            L.append(sols)
        cost = {b: (0.0, [b]) for b in range(len(L[0]))}
        for k in range(1, len(L)):
            body = BODIES[codes[k]]
            cost = {bj: min(((cp + flyby_dv(L[k - 1][bp][2], sj[1], body), path + [bj])
                             for bp, (cp, path) in cost.items()), key=lambda x: x[0])
                    for bj, sj in enumerate(L[k])}
        tot, path = min(cost.values(), key=lambda x: x[0])
        return tot, [L[k][b] for k, b in enumerate(path)]

    dv0, legs0 = chain(t_paper)
    revs = [l[0] for l in legs0]
    res = minimize(lambda x: chain(t_paper + np.asarray(x) * DAY, revs)[0], np.zeros(len(t_paper)),
                   method="Powell", bounds=[(-shift_days, shift_days)] * len(t_paper),
                   options=dict(maxiter=4000, xtol=1e-3, ftol=1e-6))
    t = t_paper + res.x * DAY
    dv, L = chain(t, revs)
    rows = []
    for k, (b, date, v_p, alt_p, _leg) in enumerate(S1L1_DE405):
        v_in = L[0][1] if k == 0 else L[k - 1][2]          # launch: departure V∞
        v = np.linalg.norm(v_in) / 1e3
        alt = None
        if 0 < k < len(L):
            body = BODIES[b]
            delta = angle_between(L[k - 1][2], L[k][1])
            vv = np.linalg.norm(L[k - 1][2])
            alt = (body.mu / vv**2 * (1 / np.sin(delta / 2) - 1) - body.radius) / 1e3
        d = dict(encounter=k + 1, body=b, paper_date=date, date_shift_days=round(float(res.x[k]), 2),
                 paper_vinf_kms=v_p, ours_vinf_kms=round(v, 3),
                 paper_altitude_km=alt_p, ours_altitude_km=None if alt is None else round(alt),
                 flyby_dv_kms=None if not 0 < k < len(L) else round(flyby_dv(L[k - 1][2], L[k][1], BODIES[b]) / 1e3, 4),
                 leg_revs=revs[k] if k < len(revs) else None)
        d["pass_vinf"] = abs(d["ours_vinf_kms"] - v_p) <= 0.1
        d["pass_altitude"] = alt_p is None or (alt is not None and abs(alt - alt_p) <= 0.10 * alt_p)
        d["pass"] = d["pass_vinf"] and d["pass_altitude"]
        rows.append(d)
    df = pd.DataFrame(rows)
    summary = dict(total_flyby_dv_at_paper_dates_kms=round(dv0 / 1e3, 3), total_flyby_dv_kms=round(dv / 1e3, 4),
                   max_date_shift_days=round(float(np.abs(res.x).max()), 2),
                   max_vinf_error_kms=round(float((df.ours_vinf_kms - df.paper_vinf_kms).abs().max()), 3),
                   median_altitude_error_pct=round(float(((df.ours_altitude_km - df.paper_altitude_km).abs()
                                                          / df.paper_altitude_km).dropna().median() * 100), 1),
                   max_altitude_error_pct=round(float(((df.ours_altitude_km - df.paper_altitude_km).abs()
                                                       / df.paper_altitude_km).dropna().max() * 100), 1),
                   rows_passing=int(df["pass"].sum()), rows=len(df))
    if with_traj:
        summary["_legs"] = [(t[k], state(codes[k], t[k])[0], L[k][3], t[k + 1] - t[k]) for k in range(len(L))]
        summary["_enc"] = list(zip(codes, t))
    return df, summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/validation_mcconaghy.csv")
    args = ap.parse_args()
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    df, checks = run(args.out)
    cols = ["notation", "paper_aphelion_au", "ours_aphelion_au", "paper_vinf_earth_kms", "ours_vinf_earth_kms",
            "paper_vinf_mars_kms", "ours_vinf_mars_kms", "paper_transfer_days", "ours_transfer_days",
            "paper_turn_required_deg", "ours_turn_required_deg", "paper_turn_max_deg", "ours_turn_max_deg", "pass"]
    print("[M02]/[M04] Table 4\n" + df[cols].to_string(index=False))
    print(f"Table 4 rows passing (all columns): {int(df['pass'].sum())}/{len(df)}\n")
    for k, (ours, ref, tol) in checks.items():
        print(f"{k:26s} ours={ours:8.3f} paper={ref:8.3f} {'PASS' if abs(ours - ref) <= tol else 'FAIL'}")
    t6 = run_table6_2002()
    print("\n[M02] Table 6\n" + t6.drop(columns=[c for c in t6.columns if c.startswith("_")]).to_string(index=False))
    s1, summ = run_s1l1_de405()
    print("\n[M04] Table 6 (DE405 in the paper, VSOP2013 here)\n" + s1.to_string(index=False))
    print(summ)
    sys.stdout.flush()
    import os; os._exit(0)  # skip pykep teardown crash (see README known issues)
