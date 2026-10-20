# Copyright (c) 2026 Sharif Mutasim
# SPDX-License-Identifier: MIT
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import cyclerdb  # noqa: E402,F401
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from cyclerdb.bodies import BODIES, YEAR, DAY, idealized_bodies, synodic_period, EphemerisModel, epoch_seconds  # noqa: E402
from cyclerdb.astro import flyby_dv, max_turn, rot_z  # noqa: E402
from cyclerdb.single_home import enumerate_single_home  # noqa: E402
from cyclerdb.sequence import SequenceProblem, solve, candidate_periods, drift_per_cycle  # noqa: E402
from cyclerdb.continuation import IdealCycler, launch_epochs, EphemerisChain  # noqa: E402
from cyclerdb.catalog import sequences, canonical  # noqa: E402


# --- literature validation -------------------------------------------------
def test_mcconaghy_table1_all_rows():
    from validate_mcconaghy import run
    df, checks = run()
    assert df["pass"].all(), df[~df["pass"]]
    for name, (ours, ref, tol) in checks.items():
        assert abs(ours - ref) <= tol, name


def test_mcconaghy_notation():
    """Aldrin is "1L1"; the 2.81 AU six-synodic cycler is "6S4" (L = higher energy branch)."""
    from validate_mcconaghy import run
    df, _ = run()
    lab = {(r.n_synodic, r.n_revs, r.paper_aphelion_au): r.notation for r in df.itertuples()}
    assert lab[(1, 1, 2.23)] == "1L1"
    assert lab[(6, 4, 2.81)] == "6S4"
    assert lab[(3, 5, 1.45)] == "3L5" and lab[(3, 5, 1.52)] == "3S5"


def test_mcconaghy_2002_table6_intermediate_flyby():
    from validate_mcconaghy import run_table6_2002
    df = run_table6_2002()
    assert df["pass"].all(), df


@pytest.mark.slow
def test_mcconaghy_2004_s1l1_real_ephemeris():
    """The 30-year DE405 S1L1 itinerary is reproduced on VSOP2013: every V∞ within
    0.1 km/s and nearly ballistic after moving each date by at most 1.5 days."""
    from validate_mcconaghy import run_s1l1_de405
    df, s = run_s1l1_de405()
    assert df["pass"].all(), df[~df["pass"]]
    assert s["total_flyby_dv_kms"] < 0.1


def test_single_home_leg_revs_follow_the_true_arc():
    """Regression: the per-leg revolution counts used for ephemeris continuation must
    reproduce the single-assist arc itself (the sequence-solver route sometimes picked
    a cheaper, different Lambert branch for the second leg)."""
    import pykep as pk
    from cyclerdb.single_home import departure_state
    from cyclerdb.astro import lambert_all
    for H, T in [("E", "M"), ("M", "E"), ("V", "E"), ("E", "V")]:
        for c in enumerate_single_home(H, T, 4, 12):
            if not c.reaches_target:
                continue
            ic = IdealCycler.from_single_home(c)
            r0, v0 = departure_state(c)
            tf = ic.t_enc[1]
            rT, vT = pk.propagate_lagrangian([list(r0), list(v0)], tf, cyclerdb.bodies.MU_SUN)
            th = 2 * np.pi * ic.Tc / BODIES[H].period
            r1 = BODIES[H].a * np.array([np.cos(th), np.sin(th), 0.0])
            sols = [v for n, _, v, _ in lambert_all(np.array(rT), r1, ic.Tc - tf, ic.leg_revs[1]) if n == ic.leg_revs[1]]
            assert min(np.linalg.norm(v - np.array(vT)) for v in sols) < 1.0, c.notation


def test_synodic_periods():
    assert synodic_period(BODIES["E"], BODIES["M"]) / YEAR == pytest.approx(2.135, abs=2e-3)
    assert synodic_period(BODIES["E"], BODIES["V"]) / YEAR == pytest.approx(1.599, abs=2e-3)
    ib = idealized_bodies()
    assert 7 * synodic_period(ib["E"], ib["M"]) / YEAR == pytest.approx(15.0, abs=1e-9)


# --- physics helpers --------------------------------------------------------
def test_flyby_dv_zero_inside_turn_limit():
    E = BODIES["E"]
    v = np.array([5000.0, 0, 0])
    d = 0.5 * max_turn(5000.0, E)
    assert flyby_dv(v, rot_z(d) @ v, E) < 1e-6
    assert flyby_dv(v, rot_z(2 * max_turn(5000.0, E)) @ v, E) > 100.0


def test_ephemeris_is_ecliptic():
    m = EphemerisModel()
    r, v = m.state("E", epoch_seconds("2035-06-01"))
    h = np.cross(r, v)
    assert np.degrees(np.arccos(h[2] / np.linalg.norm(h))) < 0.05


# --- sequence optimiser ------------------------------------------------------
def test_sequence_matches_single_home_aldrin():
    """The EM sequence objective evaluated at the Aldrin geometry equals the
    single-assist result, and the optimiser can only do as well or better."""
    al = [c for c in enumerate_single_home("E", "M", 1, 1) if c.n_revs == 1 and c.branch == "hi"][0]
    ic = IdealCycler.from_single_home(al)
    p = SequenceProblem("EM", ic.Tc, max_revs=2)
    x = [ic.t_enc[1] / ic.Tc, 1 - ic.t_enc[1] / ic.Tc, ic.phases["M"] % (2 * np.pi)]
    assert p.evaluate(x) / 1e3 == pytest.approx(al.dv_per_cycle_kms, abs=1e-6)
    best = min(c.dv_per_cycle_kms for c in solve("EM", ic.Tc, seeds=(0,), maxiter=80))
    assert best <= al.dv_per_cycle_kms + 1e-6


def test_triple_period_drift_small():
    Ts = candidate_periods("EMV")
    assert any(abs(T / YEAR - 6.41) < 0.02 for T in Ts)
    d = drift_per_cycle("EMV", Ts[0])
    assert abs(d["M"]) < 1e-9 and abs(d["V"]) < 10


def test_sequence_enumeration_canonical():
    assert canonical("MEV") == "EVM"
    assert sequences("EM", 3) == ["EM", "EEM", "EMM"]


# --- ephemeris continuation --------------------------------------------------
def test_launch_windows_one_per_synodic_period():
    al = [c for c in enumerate_single_home("E", "M", 1, 1) if c.n_revs == 1 and c.branch == "hi"][0]
    ic = IdealCycler.from_single_home(al)
    eps = np.array(launch_epochs(ic, "2030-01-01", "2045-01-01"))
    gaps = np.diff(eps) / YEAR
    assert len(eps) == 7
    assert np.all(np.abs(gaps - 2.135) < 0.2)  # eccentric orbits make the spacing uneven


def test_continuation_reduces_cost():
    al = [c for c in enumerate_single_home("E", "M", 1, 1) if c.n_revs == 1 and c.branch == "hi"][0]
    ic = IdealCycler.from_single_home(al)
    tau = launch_epochs(ic)[0]
    ch = EphemerisChain(ic, tau, n_cycles=2)
    x0 = np.zeros(len(ch.t_guess))
    res = ch.optimise()
    assert ch.cost(res.x) <= ch.cost(x0)
    r = ch.result(res.x)
    assert np.isfinite(r.total_dv_kms) and set(r.max_vinf_kms) == {"E", "M"}
