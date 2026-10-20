# Copyright (c) 2026 Sharif Mutasim
# SPDX-License-Identifier: MIT
"""Single-assist cyclers (Hollister / Byrnes / McConaghy construction).

A cycler is one Lambert arc that leaves the *home* planet H and returns to it
after exactly n synodic periods of the (H, target) pair, making r full solar
revolutions on the way. Because the H–target geometry repeats every synodic
period, the same arc rotated by H's advance angle is flown on every cycle; the
gravity assist at H must rotate the v-infinity vector by that angle. If the
planet can't turn it far enough, a powered-flyby ΔV is needed each cycle.

The target is met wherever the arc crosses its orbit; the epoch is chosen so
the target is there. With circular orbits one encounter per cycle is
guaranteed (outbound or inbound, depending on phasing).

With ``bodies=idealized_bodies()`` and H = Earth, target = Mars this
reproduces Table 1 of McConaghy, Longuski & Byrnes, AIAA 2002-4420.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

from .astro import (lambert_all, conic, radius_crossings, max_turn, angle_between,
                    flyby_dv, rot_z)
from .bodies import BODIES, AU, DAY, YEAR, ECCENTRICITY, CircularModel, synodic_period


@dataclass
class SingleHomeCycler:
    family: str
    sequence: str            # e.g. "E(M)": home planet, target in brackets
    home: str
    target: str
    n_synodic: int
    n_revs: int              # full solar revolutions per cycle
    branch: str              # 'hi' / 'lo' energy solution for this n_revs
    notation: str            # McConaghy et al. (2002) "nXr": n synodic periods, L/S branch, r revs
    period_years: float
    a_au: float
    e: float
    perihelion_au: float
    aphelion_au: float
    vinf_home_kms: float
    vinf_target_kms: float | None
    transfer_days: float | None      # home -> first crossing of target orbit
    aphelion_days: float             # home -> first aphelion (closest approach when the arc
                                     # stops short of the target's circular orbit)
    turn_required_deg: float
    turn_max_deg: float
    dv_per_cycle_kms: float          # powered-flyby ΔV at the home flyby
    reaches_target: bool             # crosses the circular target orbit
    reaches_target_eccentric: bool   # overlaps the real target's [q, Q] range
    ballistic: bool

    def as_row(self):
        return asdict(self)


def departure_state(c: "SingleHomeCycler", bodies: dict | None = None):
    """Heliocentric (r0, v0) of the cycler leaving its home planet at longitude 0, t = 0."""
    bodies = bodies or BODIES
    H = bodies[c.home]
    Tc = c.period_years * YEAR
    th = 2 * np.pi * Tc / H.period
    r0 = np.array([H.a, 0.0, 0.0])
    r1 = H.a * np.array([np.cos(th), np.sin(th), 0.0])
    vH = np.array([0.0, H.v_circ, 0.0])
    sols = [v0 for n, _b, v0, _v1 in lambert_all(r0, r1, Tc, c.n_revs) if n == c.n_revs]
    v0 = min(sols, key=lambda v: abs(np.linalg.norm(v - vH) / 1e3 - c.vinf_home_kms))
    return r0, v0


def notation(n_synodic: int, branch: str, n_revs: int) -> str:
    """McConaghy, Longuski & Byrnes (2002) cycler label, e.g. "1L1" (Aldrin) or "6S4".

    n synodic periods per cycle, then L (long period: the higher-energy, larger
    semi-major axis Lambert solution) or S (short period), then r full solar
    revolutions. The single solution with fewer than one revolution is "U"
    (unique period), e.g. "1U0", as in the paper.
    """
    return f"{n_synodic}{dict(hi='L', lo='S').get(branch, 'U')}{n_revs}"


def enumerate_single_home(home: str, target: str, max_synodic: int = 7,
                          max_revs: int = 12, dv_tol_kms: float = 1e-3,
                          bodies: dict | None = None) -> list[SingleHomeCycler]:
    bodies = bodies or BODIES
    H, T = bodies[home], bodies[target]
    S = synodic_period(H, T)
    model = CircularModel(bodies=bodies)
    q_T = T.a * (1 - ECCENTRICITY[target])
    Q_T = T.a * (1 + ECCENTRICITY[target])
    out = []
    for n in range(1, max_synodic + 1):
        Tc = n * S
        dtheta = 2 * np.pi * Tc / H.period
        frac = (Tc / H.period) % 1.0
        if min(frac, 1 - frac) < 1e-4 or abs(frac - 0.5) < 1e-4:
            continue  # 0 or 180 deg: a resonant return, not a Lambert arc
        r0, vH0 = model.state(home, 0.0)
        r1, vH1 = model.state(home, Tc)
        sols = lambert_all(r0, r1, Tc, max_revs)
        # label the two solutions of each revolution count by energy
        by_rev: dict[int, list] = {}
        for n_rev, _, v0, v1 in sols:
            a = conic(r0, v0)[0]
            if a > 0:
                by_rev.setdefault(n_rev, []).append((a, v0, v1))
        for n_rev, lst in by_rev.items():
            lst.sort(key=lambda x: -x[0])
            for i, (a, v0, v1) in enumerate(lst):
                branch = "hi" if (len(lst) == 2 and i == 0) else ("lo" if len(lst) == 2 else "-")
                _, e, rp, ra, _ = conic(r0, v0)
                vinf_dep = v0 - vH0
                vinf_arr = v1 - vH1
                vinf_next = rot_z(dtheta) @ vinf_dep
                vh = float(np.linalg.norm(vinf_dep))
                if vh < 1.0:
                    continue  # trivial solution: the planet's own orbit
                dv = flyby_dv(vinf_arr, vinf_next, H)
                cross = radius_crossings(r0, v0, T.a)
                if cross is not None:
                    t_first = min(cross["t_out"], cross["t_in"])
                    s, g, vT = cross["speed"], cross["gamma"], T.v_circ
                    vinf_T = np.sqrt(s**2 + vT**2 - 2 * s * vT * np.cos(g))
                    reaches = t_first < Tc
                else:
                    t_first, vinf_T, reaches = None, None, False
                reaches_ecc = bool(rp <= Q_T and ra >= q_T)
                aph = radius_crossings(r0, v0, ra * (1 - 1e-9))
                out.append(SingleHomeCycler(
                    family="single_home", sequence=f"{home}({target})", home=home, target=target,
                    n_synodic=n, n_revs=n_rev, branch=branch, notation=notation(n, branch, n_rev),
                    period_years=Tc / YEAR,
                    a_au=a / AU, e=e, perihelion_au=rp / AU, aphelion_au=ra / AU,
                    vinf_home_kms=vh / 1e3,
                    vinf_target_kms=None if vinf_T is None else vinf_T / 1e3,
                    transfer_days=None if t_first is None else t_first / DAY,
                    aphelion_days=None if aph is None else aph["t_out"] / DAY,
                    turn_required_deg=np.degrees(angle_between(vinf_arr, vinf_next)),
                    turn_max_deg=np.degrees(max_turn(vh, H)),
                    dv_per_cycle_kms=dv / 1e3, reaches_target=bool(reaches),
                    reaches_target_eccentric=reaches_ecc,
                    ballistic=bool((reaches or reaches_ecc) and dv / 1e3 <= dv_tol_kms)))
    return out
