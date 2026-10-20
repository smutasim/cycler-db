# Copyright (c) 2026 Sharif Mutasim
# SPDX-License-Identifier: MIT
"""Build the V/E/M cycler catalogue.

Stages (each resumable; results are cached under data/):
  1. single-assist cyclers for every ordered planet pair (fast, exhaustive)
  2. multi-assist sequence cyclers for EM, EV, MV and EMV (global search)
  3. ephemeris continuation (VSOP2013) of the best candidates, 2030-2060

Usage:
  python scripts/build_catalog.py [--workers 2] [--seeds 2] [--stage all|1|2|3]
"""
import argparse
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import cyclerdb  # noqa: E402,F401  (import before pandas)
import pandas as pd  # noqa: E402
from cyclerdb.bodies import BODIES, YEAR  # noqa: E402
from cyclerdb.single_home import enumerate_single_home  # noqa: E402
from cyclerdb.sequence import solve  # noqa: E402
from cyclerdb.catalog import search_space, dedupe, write_sqlite, add_entry_costs, real_summary  # noqa: E402
import numpy as np  # noqa: E402
from cyclerdb.continuation import IdealCycler, continue_cycler  # noqa: E402
from cyclerdb.lowthrust import min_accel_screen  # noqa: E402

DATA = ROOT / "data"


def stage1():
    rows = []
    for h in "VEM":
        for t in "VEM":
            if h == t:
                continue
            S = 1 / abs(1 / BODIES[h].period - 1 / BODIES[t].period)
            nmax = max(1, int(16 * YEAR // S))
            rows += [c.as_row() for c in enumerate_single_home(h, t, max_synodic=nmax, max_revs=12)]
    df = pd.DataFrame(rows)
    df.insert(0, "cycler_id", [f"SH-{r.home}{r.target}-n{r.n_synodic}r{r.n_revs}{r.branch}" for r in df.itertuples()])
    df.to_csv(DATA / "single_home.csv", index=False)
    print(f"stage 1: {len(df)} single-assist cyclers, {int(df.ballistic.sum())} ballistic")
    return df


def _solve_one(args):
    seq, Tc, max_revs, seeds, maxiter = args
    t0 = time.time()
    try:
        res = solve(seq, Tc, seeds=seeds, max_revs=max_revs, maxiter=maxiter, popsize=15)
        rows = [r.as_row() for r in res]
    except Exception as e:  # keep the batch going
        rows = [dict(sequence=seq, period_years=Tc / YEAR, error=repr(e))]
    return dict(seq=seq, Tc=Tc, rows=rows, secs=time.time() - t0)


def stage2(workers, seeds, maxiter):
    cache = DATA / "sequence_raw.jsonl"
    done = set()
    if cache.exists():
        for line in cache.read_text().splitlines():
            d = json.loads(line)
            done.add((d["seq"], round(d["Tc"], 0)))
    probs = [p for p in search_space() if (p[0], round(p[1], 0)) not in done]
    print(f"stage 2: {len(probs)} problems to solve ({len(done)} cached)")
    jobs = [(s, Tc, mr, tuple(range(seeds)), maxiter) for s, Tc, mr in probs]
    t0 = time.time()
    with Pool(workers) as pool, cache.open("a") as fh:
        for i, out in enumerate(pool.imap_unordered(_solve_one, jobs), 1):
            fh.write(json.dumps(out) + "\n"); fh.flush()
            best = min((r.get("dv_per_cycle_kms", 9e9) for r in out["rows"]), default=float("nan"))
            print(f"  [{i}/{len(jobs)}] {out['seq']:5s} {out['Tc']/YEAR:5.2f} yr  best ΔV {best:.3f} km/s  "
                  f"({out['secs']:.0f}s, elapsed {time.time()-t0:.0f}s)", flush=True)
    rows = []
    for line in cache.read_text().splitlines():
        rows += [r for r in json.loads(line)["rows"] if "error" not in r]
    rows = dedupe(rows)
    df = pd.DataFrame(rows)
    df = df.sort_values(["planets", "dv_per_cycle_kms", "max_vinf_kms"]).reset_index(drop=True)
    df.insert(0, "cycler_id", [f"SQ-{r.sequence}-{r.period_years:.2f}y-{i}" for i, r in enumerate(df.itertuples())])
    df.to_csv(DATA / "sequence.csv", index=False)
    print(f"stage 2: {len(df)} distinct sequence cyclers, {int(df.ballistic.sum())} ballistic")
    return df


def _continue_one(args):
    ic, n_cycles = args
    t0 = time.time()
    out = []
    for r in continue_cycler(ic, n_cycles=n_cycles):
        d = r.as_row()
        d["period_years"] = ic.Tc / YEAR
        # SEP screen: worst-cycle ΔV spread over one cycle
        d["sep_min_accel_mms2"] = min_accel_screen(r.worst_cycle_dv_kms, ic.Tc / 86400.0)
        out.append(d)
    return ic.label, out, time.time() - t0


def _candidates(sh, sq, scope, top_per_system=4, vinf_max=9.0):
    """Ideal cyclers to continue into the real ephemeris.

    scope="all": every single-assist cycler that reaches its target plus every
    sequence cycler (home-only orbits are not cyclers and are skipped).
    scope="top": the original quick screen, the few cheapest per planet set.
    """
    if scope == "all":
        sh_sel = sh[sh.reaches_target.astype(bool)]
        sq_sel = sq
    else:
        ok = sh[(sh.reaches_target) & (sh.dv_per_cycle_kms <= 0.5) & (sh.vinf_home_kms <= vinf_max)
                & (sh.vinf_target_kms <= vinf_max)]
        sh_sel = pd.concat([g.sort_values(["dv_per_cycle_kms", "period_years"]).head(2)
                            for _, g in ok.groupby(["home", "target"])])
        good = sq[(sq.dv_per_cycle_kms <= 0.5) & (sq.max_vinf_kms <= vinf_max)]
        sq_sel = pd.concat([g.head(top_per_system) for _, g in good.groupby("planets")])
    cands = []
    sols = {}
    for _, r in sh_sel.iterrows():
        key = (r.home, r.target)
        if key not in sols:
            S = 1 / abs(1 / BODIES[r.home].period - 1 / BODIES[r.target].period)
            sols[key] = {(c.n_synodic, c.n_revs, c.branch): c
                         for c in enumerate_single_home(r.home, r.target, max(1, int(16 * YEAR // S)), 12)}
        c = sols[key].get((r.n_synodic, r.n_revs, r.branch))
        try:
            ic = IdealCycler.from_single_home(c)
        except Exception as e:  # noqa: BLE001  (record and skip; never fatal)
            print(f"  skip {r.cycler_id}: {e!r}")
            continue
        ic.label = r.cycler_id
        cands.append(ic)
    for _, r in sq_sel.iterrows():
        cands.append(IdealCycler.from_sequence_row(r))
    return cands


def _keep_awake():
    """Stop Windows from sleeping while the (long) continuation runs."""
    if sys.platform == "win32":
        import ctypes
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000000 | 0x00000001)  # CONTINUOUS | SYSTEM_REQUIRED


def stage3(sh, sq, workers=2, n_cycles=4, scope="all"):
    """Continue ideal cyclers into the real ephemeris (resumable, cached per cycler)."""
    _keep_awake()
    cache = DATA / "ephemeris_raw.jsonl"
    done = {}
    if cache.exists():
        for line in cache.read_text(encoding="utf-8").splitlines():
            d = json.loads(line)
            done[d["label"]] = d["rows"]
    cands = _candidates(sh, sq, scope)
    todo = [ic for ic in cands if ic.label not in done]
    print(f"stage 3: {len(cands)} cyclers into VSOP2013 ({len(done)} cached, {len(todo)} to run), "
          f"{n_cycles} cycles per window (2 if period > 8 yr), 2030-2060", flush=True)
    # longest jobs first (more encounters per cycle) so the pool drains evenly
    todo.sort(key=lambda ic: -len(ic.seq) * (n_cycles if ic.Tc < 8 * YEAR else 2))
    jobs = [(ic, n_cycles if ic.Tc < 8 * YEAR else 2) for ic in todo]
    t0 = time.time()
    with Pool(workers) as pool, cache.open("a", encoding="utf-8") as fh:
        for i, (label, res_rows, secs) in enumerate(pool.imap_unordered(_continue_one, jobs), 1):
            res_rows = [{k: (None if isinstance(v, float) and not np.isfinite(v) else v) for k, v in r.items()}
                        for r in res_rows]
            fh.write(json.dumps(dict(label=label, rows=res_rows)) + "\n"); fh.flush()
            done[label] = res_rows
            best = min((r["dv_per_cycle_kms"] for r in res_rows if r["dv_per_cycle_kms"] is not None), default=float("nan"))
            print(f"  [{i}/{len(jobs)}] {label:28s} {len(res_rows):2d} windows  best {best:.3f} km/s/cycle "
                  f"({secs:.0f}s, elapsed {time.time() - t0:.0f}s)", flush=True)
    labels = {ic.label for ic in cands}
    rows = [r for lab, rs in done.items() if lab in labels for r in rs]
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "ephemeris.csv", index=False)
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--maxiter", type=int, default=150)
    ap.add_argument("--stage", default="all")
    ap.add_argument("--real-scope", choices=["all", "top"], default="all",
                    help="stage 3: continue every cycler (all, ~1 h on 15 cores) or only a quick top few")
    a = ap.parse_args()
    DATA.mkdir(exist_ok=True)
    sh = stage1() if a.stage in ("all", "1") else pd.read_csv(DATA / "single_home.csv")
    if a.stage in ("all", "2"):
        sq = stage2(a.workers, a.seeds, a.maxiter)
    else:
        sq = pd.read_csv(DATA / "sequence.csv") if (DATA / "sequence.csv").exists() else None
    if a.stage in ("all", "3") and sq is not None:
        ep = stage3(sh, sq, a.workers, scope=a.real_scope)
    else:
        ep = pd.read_csv(DATA / "ephemeris.csv") if (DATA / "ephemeris.csv").exists() else pd.DataFrame()
    tables = {"single_home": add_entry_costs(sh, "single_home")}
    if sq is not None:
        tables["sequence"] = add_entry_costs(sq, "sequence")
    if len(ep):
        tables["ephemeris"] = add_entry_costs(ep, "ephemeris")
        tables["real_summary"] = real_summary(tables["ephemeris"])
    sys.path.insert(0, str(ROOT / "scripts"))
    from validate_mcconaghy import run as validate, run_table6_2002, run_s1l1_de405  # noqa: E402
    tables["validation_table4"], _ = validate(DATA / "validation_mcconaghy.csv")
    t6 = run_table6_2002()
    tables["validation_table6_2002"] = t6.astype({c: str for c in t6.columns if c.endswith("_yr") and "crossings" in c})
    tables["validation_s1l1_de405"], _ = run_s1l1_de405()
    write_sqlite(DATA / "cyclers.sqlite", tables)
    print("wrote", DATA / "cyclers.sqlite", {k: len(v) for k, v in tables.items()})
    sys.stdout.flush()
    import os; os._exit(0)  # skip pykep teardown crash (see README known issues)
