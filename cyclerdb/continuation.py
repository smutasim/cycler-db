# Copyright (c) 2026 Sharif Mutasim
# SPDX-License-Identifier: MIT
"""Continue an idealised cycler into the real (eccentric, inclined) solar system.

Steps
1. Express the cycler as an encounter sequence with ideal encounter times
   ``t_i`` within one period ``Tc`` and ideal planet phases (longitude of each
   planet relative to the home planet at t = 0).
2. Find real epochs in a date window where the actual planets match those
   relative phases (one opportunity per synodic period for a pair cycler).
3. Starting from each such epoch, chain ``n_cycles`` cycles of Lambert arcs
   between real planet positions (VSOP2013 or DE440) and optimise all
   encounter epochs to minimise the total powered-flyby ΔV. This is the
   standard patched-conic "ephemeris version" of a cycler, as in Byrnes et
   al. (1993) and Russell & Ocampo (2004/2005).

Result: per launch window, ΔV per cycle, worst-case v-infinity at each planet
and the encounter dates — the data a mission designer actually needs.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
import pykep as pk
from scipy.optimize import brentq, minimize

from .astro import lambert_all, flyby_dv, conic
from .bodies import BODIES, DAY, YEAR, MU_SUN, EphemerisModel, epoch_seconds


@dataclass
class IdealCycler:
    """Minimal description needed to continue a cycler."""
    seq: str                     # encounter sequence, home first
    Tc: float                    # period, s
    t_enc: list[float]           # encounter times in [0, Tc), first is 0, s
    phases: dict[str, float]     # ecliptic longitude of each planet minus home's at t=0, rad
    leg_revs: list[int]          # full revolutions on each leg
    label: str = ""
    vinf_kms: list[float] | None = None  # ideal v∞ at each encounter (keeps the same family)

    @staticmethod
    def from_sequence_row(row) -> "IdealCycler":
        legs = [float(v) * DAY for v in str(row["leg_days"]).split("|")]
        t = list(np.concatenate([[0.0], np.cumsum(legs)[:-1]]))
        ph = {kv.split(":")[0]: np.radians(float(kv.split(":")[1]))
              for kv in str(row["phases_deg"]).split("|")}
        revs = [int(v) for v in str(row["leg_revs"]).split("|")]
        vinf = [float(v) for v in str(row["vinf_kms"]).split("|")]
        return IdealCycler(row["sequence"], row["period_years"] * YEAR, t, ph, revs,
                           label=str(row.get("cycler_id", "")), vinf_kms=vinf)

    @staticmethod
    def from_single_home(c, bodies=BODIES) -> "IdealCycler":
        """Convert a single-assist cycler (one Lambert arc) into sequence form."""
        from .astro import lambert_all as _lam
        H, T = bodies[c.home], bodies[c.target]
        Tc = c.period_years * YEAR
        r0 = np.array([H.a, 0.0, 0.0]); vH = np.array([0.0, H.v_circ, 0.0])
        th = 2 * np.pi * Tc / H.period
        r1 = H.a * np.array([np.cos(th), np.sin(th), 0.0])
        # pick the matching Lambert solution
        best = None
        for n_rev, _b, v0, _v1 in _lam(r0, r1, Tc, c.n_revs):
            if n_rev != c.n_revs:
                continue
            d = abs(np.linalg.norm(v0 - vH) / 1e3 - c.vinf_home_kms)
            if best is None or d < best[0]:
                best = (d, v0)
        v0 = best[1]
        tf = c.transfer_days * DAY
        rT, _ = pk.propagate_lagrangian([list(r0), list(v0)], tf, MU_SUN)
        lam_T = np.arctan2(rT[1], rT[0])
        phase_T = lam_T - 2 * np.pi * tf / T.period
        # full revolutions on each sub-leg: a conic arc of duration dt completes
        # floor(dt / P) whole revolutions plus a partial one (< 2π), which is the
        # Lambert revolution count. (Some V–E arcs reach the target after only
        # ~13 days, too short for the sequence solver's 20-day minimum leg.)
        P = conic(r0, v0)[4]
        revs = [int(tf // P), int((Tc - tf) // P)]
        return IdealCycler(c.home + c.target, Tc, [0.0, tf], {c.home: 0.0, c.target: phase_T},
                           revs, label=f"{c.sequence} n{c.n_synodic} r{c.n_revs} {c.branch}",
                           vinf_kms=[c.vinf_home_kms, c.vinf_target_kms])


def _wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def launch_epochs(ic: IdealCycler, start: str = "2030-01-01", end: str = "2060-01-01",
                  model: EphemerisModel | None = None, step_days: float = 2.0,
                  tol_deg: float = 15.0) -> list[float]:
    """Real epochs (s past J2000) where planet phases match the ideal cycler."""
    model = model or EphemerisModel()
    home = ic.seq[0]
    others = [c for c in sorted(set(ic.seq)) if c != home]
    t0, t1 = epoch_seconds(start), epoch_seconds(end)
    ts = np.arange(t0, t1, step_days * DAY)

    def err(t, c):
        return _wrap(model.longitude(c, t) - model.longitude(home, t) - ic.phases[c])

    if len(others) == 1:
        c = others[0]
        e = np.array([err(t, c) for t in ts])
        out = []
        for i in range(len(ts) - 1):
            if np.sign(e[i]) != np.sign(e[i + 1]) and abs(e[i] - e[i + 1]) < np.pi:
                out.append(brentq(lambda t: err(t, c), ts[i], ts[i + 1], xtol=60.0))
        return out
    # three planets: local minima of the worst phase error
    e = np.array([max(abs(err(t, c)) for c in others) for t in ts])
    out = []
    for i in range(1, len(ts) - 1):
        if e[i] <= e[i - 1] and e[i] <= e[i + 1] and np.degrees(e[i]) <= tol_deg:
            out.append(float(ts[i]))
    return out


@dataclass
class EphemerisResult:
    label: str
    sequence: str
    start_epoch: str
    n_cycles: int
    total_dv_kms: float
    dv_per_cycle_kms: float
    worst_cycle_dv_kms: float
    max_vinf_kms: dict
    encounter_dates: str
    converged: bool
    max_vinf_deviation_kms: float = float("nan")  # vs. the ideal cycler (same family check)
    earth_vinf_kms: float = float("nan")  # departure v∞ at the first Earth encounter (boarding)

    def as_row(self):
        d = asdict(self)
        d["max_vinf_kms"] = "|".join(f"{k}:{v:.2f}" for k, v in self.max_vinf_kms.items())
        return d


class EphemerisChain:
    def __init__(self, ic: IdealCycler, tau0: float, n_cycles: int = 4,
                 model: EphemerisModel | None = None, window_days: float = 45.0,
                 rev_slack: int = 0):
        self.ic, self.N = ic, n_cycles
        self.model = model or EphemerisModel()
        m = len(ic.seq)
        self.codes, self.t_guess, self.revs, self.vinf_ideal = [], [], [], []
        for c in range(n_cycles):
            for i in range(m):
                self.codes.append(ic.seq[i])
                self.vinf_ideal.append(ic.vinf_kms[i] * 1e3 if ic.vinf_kms else None)
                self.t_guess.append(tau0 + c * ic.Tc + ic.t_enc[i])
                self.revs.append(ic.leg_revs[i])
        self.codes.append(ic.seq[0])
        self.t_guess.append(tau0 + n_cycles * ic.Tc)
        self.t_guess = np.array(self.t_guess)
        self.m = m
        self.vinf_ideal.append(ic.vinf_kms[0] * 1e3 if ic.vinf_kms else None)
        self.window = window_days * DAY
        self.family_tol = 1500.0  # m/s: keep Lambert branches whose v∞ stay near the ideal
        self.rev_slack = rev_slack
        self._cache = {}

    def _state(self, code, t):
        key = (code, round(t, 1))
        if key not in self._cache:
            self._cache[key] = self.model.state(code, t)
        return self._cache[key]

    def _legs(self, t):
        legs = []
        for k in range(len(t) - 1):
            if t[k + 1] - t[k] < 10 * DAY:
                return None
            r0, vp = self._state(self.codes[k], t[k])
            r1, vq = self._state(self.codes[k + 1], t[k + 1])
            want = self.revs[k]
            sols = [(vi0, vi1) for n, _b, v0, v1 in lambert_all(r0, r1, t[k + 1] - t[k], want + self.rev_slack)
                    if abs(n - want) <= self.rev_slack
                    for vi0, vi1 in [(v0 - vp, v1 - vq)]
                    if np.isfinite(vi0).all() and np.linalg.norm(vi0) > 1 and np.linalg.norm(vi1) > 1]
            vi = self.vinf_ideal[k]
            if sols and vi is not None:
                dev = [abs(np.linalg.norm(a) - vi) for a, _ in sols]
                near = [sl for sl, d in zip(sols, dev) if d <= self.family_tol]
                sols = near or [sols[int(np.argmin(dev))]]
            if not sols:
                return None
            legs.append(sols)
        return legs

    def cost(self, dt_days, detail=False):
        t = self.t_guess + np.asarray(dt_days) * DAY
        legs = self._legs(t)
        if legs is None:
            return (1e6, None) if detail else 1e6
        # linear DP over branch choices; interior nodes are flybys
        cost = {b: (0.0, [b]) for b in range(len(legs[0]))}
        for k in range(1, len(legs)):
            body = BODIES[self.codes[k]]
            new = {}
            for bj, sj in enumerate(legs[k]):
                best = None
                for bp, (cp, path) in cost.items():
                    c = cp + flyby_dv(legs[k - 1][bp][1], sj[0], body)
                    if best is None or c < best[0]:
                        best = (c, path + [bj])
                new[bj] = best
            cost = new
        tot, path = min(cost.values(), key=lambda x: x[0])
        if not detail:
            return tot
        return tot, (t, legs, path)

    def optimise(self, maxiter=400):
        n = len(self.t_guess)
        w = self.window / DAY
        bounds = [(-w, w)] * n
        x0 = np.zeros(n)
        res = minimize(self.cost, x0, method="Powell", bounds=bounds,
                       options=dict(maxiter=maxiter, xtol=1e-3, ftol=1e-4))
        return res

    def result(self, x, converged=True) -> EphemerisResult:
        tot, det = self.cost(x, detail=True)
        if det is None:
            return EphemerisResult(self.ic.label, self.ic.seq, "", self.N, np.inf, np.inf, np.inf, {}, "", False)
        t, legs, path = det
        per_node = [0.0] * len(t)
        vmax = {}
        for k in range(1, len(legs)):
            vin, vout = legs[k - 1][path[k - 1]][1], legs[k][path[k]][0]
            per_node[k] = flyby_dv(vin, vout, BODIES[self.codes[k]])
            c = self.codes[k]
            vmax[c] = max(vmax.get(c, 0), np.linalg.norm(vin) / 1e3, np.linalg.norm(vout) / 1e3)
        dev = 0.0
        for k in range(len(legs)):
            if self.vinf_ideal[k] is not None:
                dev = max(dev, abs(np.linalg.norm(legs[k][path[k]][0]) - self.vinf_ideal[k]))
            if self.vinf_ideal[k + 1] is not None:
                dev = max(dev, abs(np.linalg.norm(legs[k][path[k]][1]) - self.vinf_ideal[k + 1]))
        c0 = self.codes[0]
        vmax[c0] = max(vmax.get(c0, 0), np.linalg.norm(legs[0][path[0]][0]) / 1e3)
        earth = next((np.linalg.norm(legs[k][path[k]][0]) / 1e3
                      for k in range(len(legs)) if self.codes[k] == "E"), float("nan"))
        per_cycle = [sum(per_node[c * self.m + 1:(c + 1) * self.m + 1]) for c in range(self.N)]
        dates = "|".join(f"{self.codes[k]}@{pk.epoch(t[k] / DAY).__repr__()[:10]}" for k in range(len(t)))
        return EphemerisResult(
            label=self.ic.label, sequence=self.ic.seq,
            start_epoch=pk.epoch(t[0] / DAY).__repr__()[:10], n_cycles=self.N,
            total_dv_kms=tot / 1e3, dv_per_cycle_kms=tot / 1e3 / self.N,
            worst_cycle_dv_kms=max(per_cycle) / 1e3, max_vinf_kms=vmax,
            encounter_dates=dates, converged=bool(converged and np.isfinite(tot) and tot < 1e5),
            max_vinf_deviation_kms=dev / 1e3, earth_vinf_kms=earth)


def continue_cycler(ic: IdealCycler, n_cycles: int = 4, start="2030-01-01", end="2060-01-01",
                    max_windows: int | None = None, model=None, window_days=45.0):
    """Ephemeris versions of ``ic`` for every launch opportunity in [start, end)."""
    model = model or EphemerisModel()
    out = []
    epochs = launch_epochs(ic, start, end, model)
    if max_windows:
        epochs = epochs[:max_windows]
    for tau in epochs:
        ch = EphemerisChain(ic, tau, n_cycles, model, window_days)
        res = ch.optimise()
        out.append(ch.result(res.x, converged=bool(res.success)))
    return out
