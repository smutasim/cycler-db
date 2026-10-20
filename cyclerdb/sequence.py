# Copyright (c) 2026 Sharif Mutasim
# SPDX-License-Identifier: MIT
"""Multi-assist periodic cyclers in the circular-coplanar model.

A cycler is described by an encounter sequence, e.g. ``"EM"``, ``"EVM"``,
``"EMEV"``. The spacecraft flies Lambert arcs P0 -> P1 -> ... -> P_{m-1} -> P0
and the whole picture repeats, rotated by P0's advance angle, every period
``Tc``. Every encounter is a gravity assist; the cost of a cycler is the sum
of powered-flyby ΔVs needed to patch the arcs (pykep's ``fb_dv``), so a
ballistic cycler has cost 0.

Two-planet sequences are exactly periodic when ``Tc`` is a multiple of the
synodic period. Three-planet (V-E-M) sequences can only be quasi-periodic:
``Tc`` is picked at a near-commensurability (e.g. ~12.8 yr ≈ 6 E-M ≈ 8 E-V ≈
14 V-M synodic periods) and the per-cycle phase drift of the planets that
are not P0 is reported as ``drift_deg``.

Free variables: leg durations (summing to Tc) and the phases of the planets
other than P0. The discrete choice of Lambert solution on each leg (revolution
count and branch) is solved exactly inside the objective with a cyclic
dynamic program, so the continuous optimiser only sees a scalar cost.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict, field
import itertools

import numpy as np
from scipy.optimize import differential_evolution, minimize

from .astro import lambert_all, flyby_dv, rot_z, conic
from .bodies import BODIES, AU, DAY, YEAR, synodic_period

MIN_LEG_DAYS = 20.0


# ---------------------------------------------------------------------------
# Periods
# ---------------------------------------------------------------------------
def drift_per_cycle(seq: str, Tc: float, bodies=BODIES) -> dict[str, float]:
    """Phase drift (deg/cycle) of each non-home planet relative to the rotated frame."""
    home = bodies[seq[0]]
    dth = 2 * np.pi * Tc / home.period
    out = {}
    for c in sorted(set(seq) - {seq[0]}):
        d = 2 * np.pi * Tc / bodies[c].period - dth
        out[c] = float(np.degrees((d + np.pi) % (2 * np.pi) - np.pi))
    return out


def candidate_periods(planets: str, max_years: float = 16.0, max_drift_deg: float = 10.0,
                      bodies=BODIES) -> list[float]:
    """Cycle periods (s) at which all planets in ``planets`` nearly realign.

    For two planets this is just multiples of the synodic period (drift = 0).
    For three it scans multiples of the largest-period pair's synodic period
    and keeps those where every other planet drifts < ``max_drift_deg``.
    """
    ps = sorted(set(planets))
    if len(ps) == 1:
        raise ValueError("need at least two planets")
    pairs = list(itertools.combinations(ps, 2))
    base_pair = max(pairs, key=lambda p: synodic_period(bodies[p[0]], bodies[p[1]]))
    S = synodic_period(bodies[base_pair[0]], bodies[base_pair[1]])
    out = []
    for k in range(1, int(max_years * YEAR / S) + 1):
        Tc = k * S
        if len(ps) == 2:
            out.append(Tc)
            continue
        d = drift_per_cycle(base_pair[0] + "".join(ps), Tc, bodies)
        if max(abs(v) for v in d.values()) <= max_drift_deg:
            out.append(Tc)
    return out


# ---------------------------------------------------------------------------
# Objective
# ---------------------------------------------------------------------------
@dataclass
class SequenceProblem:
    seq: str
    Tc: float
    max_revs: int = 3
    vinf_cap_kms: float = 12.0      # soft cap: penalise encounters faster than this
    bodies: dict = field(default_factory=lambda: BODIES)

    def __post_init__(self):
        self.m = len(self.seq)
        self.others = sorted(set(self.seq) - {self.seq[0]})
        self.dim = self.m + len(self.others)  # m duration weights + phases
        home = self.bodies[self.seq[0]]
        self.dth = 2 * np.pi * self.Tc / home.period
        self.R = rot_z(self.dth)

    # decision vector -> encounter times and phases
    def decode(self, x):
        w = np.asarray(x[: self.m])
        d = w / w.sum() * self.Tc
        t = np.concatenate([[0.0], np.cumsum(d)])
        phases = {self.seq[0]: 0.0}
        for i, c in enumerate(self.others):
            phases[c] = float(x[self.m + i])
        return t, d, phases

    def _state(self, code, t, phases):
        b = self.bodies[code]
        th = phases[code] + 2 * np.pi * t / b.period
        return (b.a * np.array([np.cos(th), np.sin(th), 0.0]),
                b.v_circ * np.array([-np.sin(th), np.cos(th), 0.0]))

    def legs(self, x):
        """For each leg: list of (n_rev, v_dep_inf, v_arr_inf). None if infeasible."""
        t, d, phases = self.decode(x)
        if d.min() < MIN_LEG_DAYS * DAY:
            return None, t, phases
        legs = []
        for i in range(self.m):
            p, q = self.seq[i], self.seq[(i + 1) % self.m]
            r0, vp = self._state(p, t[i], phases)
            r1, vq = self._state(q, t[i + 1], phases)
            if i == self.m - 1:  # closing leg arrives at P0 at Tc
                r1, vq = self._state(self.seq[0], self.Tc, phases)
            sols = []
            for n_rev, _br, v0, v1 in lambert_all(r0, r1, d[i], self.max_revs):
                if not (np.all(np.isfinite(v0)) and np.all(np.isfinite(v1))):
                    continue
                vi0, vi1 = v0 - vp, v1 - vq
                if np.linalg.norm(vi0) < 1.0 or np.linalg.norm(vi1) < 1.0:
                    continue  # trivial: riding the planet's own orbit
                sols.append((n_rev, vi0, vi1, v0))
            if not sols:
                return None, t, phases
            legs.append(sols)
        return legs, t, phases

    def _encounter_cost(self, body_code, vin, vout):
        b = self.bodies[body_code]
        dv = flyby_dv(vin, vout, b)
        cap = self.vinf_cap_kms * 1e3
        pen = max(0.0, np.linalg.norm(vin) - cap) + max(0.0, np.linalg.norm(vout) - cap)
        return dv + 0.5 * pen

    def evaluate(self, x, detail=False):
        legs, t, phases = self.legs(x)
        if legs is None:
            return (1e6, None) if detail else 1e6
        m = self.m
        best_total, best_path = np.inf, None
        # cyclic DP over branch choices; encounter i joins leg i-1 (arrival) and leg i (departure)
        for b0, s0 in enumerate(legs[0]):
            cost = {b0: (0.0, [b0])}
            for i in range(1, m):
                new = {}
                for bj, sj in enumerate(legs[i]):
                    best = None
                    for bp, (cp, path) in cost.items():
                        c = cp + self._encounter_cost(self.seq[i], legs[i - 1][bp][2], sj[1])
                        if best is None or c < best[0]:
                            best = (c, path + [bj])
                    new[bj] = best
                cost = new
            for bl, (cl, path) in cost.items():
                vin = legs[m - 1][bl][2]
                vout = self.R @ s0[1]
                c = cl + self._encounter_cost(self.seq[0], vin, vout)
                if c < best_total:
                    best_total, best_path = c, path
        if not detail:
            return best_total
        return best_total, (legs, best_path, t, phases)

    def bounds(self):
        return [(0.05, 1.0)] * self.m + [(0.0, 2 * np.pi)] * len(self.others)


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------
@dataclass
class SequenceCycler:
    family: str
    sequence: str
    planets: str
    n_encounters: int
    period_years: float
    dv_per_cycle_kms: float
    ballistic: bool
    max_vinf_kms: float
    vinf_kms: str            # per encounter, '|' separated
    flyby_dv_kms: str        # per encounter
    leg_days: str            # per leg
    leg_revs: str            # full revolutions per leg
    aphelion_max_au: float
    perihelion_min_au: float
    drift_deg: str           # per-cycle phase drift of non-home planets (0 for pairs)
    phases_deg: str
    seed: int

    def as_row(self):
        return asdict(self)


def describe(prob: SequenceProblem, x, seed=0, dv_tol_kms=1e-3) -> SequenceCycler | None:
    total, det = prob.evaluate(x, detail=True)
    if det is None:
        return None
    legs, path, t, phases = det
    m = prob.m
    vinf, fdv = [], []
    for i in range(m):
        vin = legs[i - 1][path[i - 1]][2] if i > 0 else legs[m - 1][path[m - 1]][2]
        vout = legs[i][path[i]][1] if i > 0 else prob.R @ legs[0][path[0]][1]
        b = prob.bodies[prob.seq[i]]
        vinf.append(max(np.linalg.norm(vin), np.linalg.norm(vout)) / 1e3)
        fdv.append(flyby_dv(vin, vout, b) / 1e3)
    aph, per = [], []
    t_enc, _, _ = prob.decode(x)
    for i in range(m):
        r0, _ = prob._state(prob.seq[i], t_enc[i], phases)
        a, e, rp, ra, _ = conic(r0, legs[i][path[i]][3])
        aph.append(ra / AU); per.append(rp / AU)
    dv = float(sum(fdv))
    drift = drift_per_cycle(prob.seq, prob.Tc, prob.bodies)
    return SequenceCycler(
        family="sequence", sequence=prob.seq, planets="".join(sorted(set(prob.seq))),
        n_encounters=m, period_years=prob.Tc / YEAR, dv_per_cycle_kms=dv,
        ballistic=dv <= dv_tol_kms, max_vinf_kms=float(max(vinf)),
        vinf_kms="|".join(f"{v:.3f}" for v in vinf),
        flyby_dv_kms="|".join(f"{v:.4f}" for v in fdv),
        leg_days="|".join(f"{d / DAY:.1f}" for d in np.diff(t)),
        leg_revs="|".join(str(legs[i][path[i]][0]) for i in range(m)),
        aphelion_max_au=float(max(aph)), perihelion_min_au=float(min(per)),
        drift_deg="|".join(f"{k}:{v:.2f}" for k, v in drift.items()),
        phases_deg="|".join(f"{k}:{np.degrees(v) % 360:.3f}" for k, v in phases.items()),
        seed=seed)


def solve(seq: str, Tc: float, seeds=(0, 1, 2), max_revs=3, vinf_cap_kms=12.0,
          maxiter=200, popsize=20, bodies=BODIES, polish=True):
    """Global search (differential evolution, several seeds) + local polish.

    Returns a list of SequenceCycler, one per seed (duplicates are removed later
    when the catalogue is built).
    """
    prob = SequenceProblem(seq, Tc, max_revs=max_revs, vinf_cap_kms=vinf_cap_kms, bodies=bodies)
    out = []
    for s in seeds:
        res = differential_evolution(prob.evaluate, prob.bounds(), seed=s, maxiter=maxiter,
                                     popsize=popsize, tol=1e-8, polish=False,
                                     mutation=(0.5, 1.0), recombination=0.9, init="sobol")
        x = res.x
        if polish:
            loc = minimize(prob.evaluate, x, method="Nelder-Mead",
                           options=dict(maxiter=4000, xatol=1e-10, fatol=1e-6))
            if loc.fun < res.fun:
                x = np.clip(loc.x, [b[0] for b in prob.bounds()], [b[1] for b in prob.bounds()])
        c = describe(prob, x, seed=s)
        if c is not None:
            out.append(c)
    return out
