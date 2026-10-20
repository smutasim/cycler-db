# Copyright (c) 2026 Sharif Mutasim
# SPDX-License-Identifier: MIT
"""Low-thrust cycler maintenance (experimental, first cut).

Status: the leg solver converges (homotopy + SLSQP), but this *local*
replacement — keep the impulsive solution's encounter dates and fix both
leg endpoints — is usually MORE expensive than the impulsive flyby it
replaces (Aldrin 2031 window: 2.7 km/s impulsive vs ~6.6 km/s low-thrust at
10 N / 20 t). A fair low-thrust cycler needs the whole chain re-optimised
with thrust arcs (an MGA-LT problem); that is the natural next step.
``min_accel_screen`` gives a quick lower bound on the SEP acceleration a
cycler's maintenance ΔV would need and is what the catalogue reports.


Idea: a powered flyby is replaced by solar-electric thrusting on the leg that
leads into the flyby. The flyby itself is kept unpowered: the spacecraft
arrives with a v-infinity the planet *can* turn into the next leg's departure
v-infinity (same magnitude, turn ≤ the maximum at the safe radius). The
low-thrust leg then has to bridge from the previous departure state to that
arrival state in the same flight time. Each leg is a Sims–Flanagan
point-to-point problem (pykep ``sf_point2point``) solved with IPOPT (pygmo).

ΔV_LT = v_eff · ln(m0 / mf) for that leg. Comparing ΔV_LT with the impulsive
ΔV shows whether deep-space low thrust is cheaper (it often is, because the
correction can be spread along the arc and placed where it is most useful),
and the thrust-to-mass needed says whether a given SEP system can do it.
"""
from __future__ import annotations

import numpy as np
import pykep as pk
import pygmo as pg

from .astro import max_turn
from .bodies import BODIES, DAY

G0 = 9.80665


def feasible_arrival_vinf(vinf_in, vinf_out, body):
    """Closest incoming v∞ that an unpowered flyby can turn into ``vinf_out``.

    Same magnitude as ``vinf_out``; direction is ``vinf_in`` rotated towards
    ``vinf_out`` until the angle between them equals the maximum turn.
    """
    vo = np.asarray(vinf_out, float)
    vi = np.asarray(vinf_in, float)
    mag = np.linalg.norm(vo)
    uo, ui = vo / mag, vi / np.linalg.norm(vi)
    dmax = max_turn(mag, body)
    ang = np.arccos(np.clip(ui @ uo, -1, 1))
    if ang <= dmax:
        return ui * mag
    # rotate uo towards ui by dmax in the plane they span
    k = np.cross(uo, ui)
    if np.linalg.norm(k) < 1e-12:
        k = np.cross(uo, [0, 0, 1.0])
    k /= np.linalg.norm(k)
    u = uo * np.cos(dmax) + np.cross(k, uo) * np.sin(dmax)
    return u * mag


def solve_leg(r0, v0, rf, vf, tof_s, m0=20000.0, max_thrust=10.0, isp=3000.0,
              nseg=None, steps=10, v_ballistic=None, max_eval=3000):
    """Minimum-propellant Sims–Flanagan leg with (nearly) fixed time of flight.

    Solved by homotopy: start from the ballistic arrival velocity
    ``v_ballistic`` (zero thrust is an exact solution there) and walk the
    target arrival velocity to ``vf`` in ``steps`` increments, warm-starting
    SLSQP (pygmo/NLopt) each time. IPOPT from random guesses was found not to
    converge reliably on these legs.

    Returns dict(dv_lt_kms, mf_kg, feasible, accel_mms2, lambda_reached).
    """
    tof_d = tof_s / DAY
    nseg = nseg or int(np.clip(tof_d / 8, 20, 120))
    veff = isp * G0
    r0, v0, rf, vf = map(np.asarray, (r0, v0, rf, vf))
    if v_ballistic is None:
        _, vb = pk.propagate_lagrangian([list(r0), list(v0)], tof_s, pk.MU_SUN)
        v_ballistic = np.asarray(vb)
    x, f, ok, lam_ok = None, None, False, 0.0
    for lam in np.linspace(0, 1, steps + 1):
        target = v_ballistic + lam * (vf - v_ballistic)
        udp = pk.trajopt.sf_point2point(
            rvs=[r0, v0], rvf=[rf, target], ms=m0, mu=pk.MU_SUN, max_thrust=max_thrust,
            veff=veff, tof_bounds=[tof_d - 0.5, tof_d + 0.5], mf_bounds=[0.05 * m0, m0],
            nseg=nseg, mass_scaling=m0)
        prob = pg.problem(udp)
        prob.c_tol = [1e-6] * (prob.get_nec() + prob.get_nic())
        lb, ub = prob.get_bounds()
        pop = pg.population(prob, 0)
        if x is None:
            x0 = np.full(prob.get_nx(), 1e-4)
            x0[0], x0[-1] = 0.999 * m0, tof_d
            pop.push_back(np.clip(x0, lb, ub))
        else:
            pop.push_back(np.clip(x, lb, ub))
        a = pg.nlopt("slsqp"); a.maxeval = max_eval; a.ftol_rel = 1e-12; a.xtol_rel = 1e-12
        pop = pg.algorithm(a).evolve(pop)
        x, f = pop.champion_x, pop.champion_f
        ok = bool(prob.feasibility_f(f))
        if not ok:
            break
        lam_ok = float(lam)
    mf = -f[0] * m0
    dv = veff * np.log(m0 / mf) if mf > 0 else np.inf
    return dict(dv_lt_kms=dv / 1e3 if ok else None, mf_kg=mf if ok else None,
                feasible=ok, accel_mms2=max_thrust / m0 * 1e3, lambda_reached=lam_ok)


def min_accel_screen(dv_kms: float, leg_days: float) -> float:
    """Lower bound (mm/s^2) on the constant acceleration that could deliver
    ``dv_kms`` over the preceding leg. A quick feasibility screen for SEP."""
    return dv_kms * 1e3 / (leg_days * DAY) * 1e3


def low_thrust_chain(chain, x, **kw):
    """Low-thrust replacement for every powered flyby in an EphemerisChain solution.

    Returns (rows, summary). Each row is one leg that needed thrusting.
    """
    tot, det = chain.cost(x, detail=True)
    t, legs, path = det
    rows = []
    for k in range(1, len(legs)):
        body = BODIES[chain.codes[k]]
        vin = legs[k - 1][path[k - 1]][1]
        vout = legs[k][path[k]][0]
        dv_imp = pk.fb_dv(list(vin), list(vout), body.mu, body.safe_radius)
        if dv_imp < 1.0:
            continue
        r0, vp0 = chain._state(chain.codes[k - 1], t[k - 1])
        rf, vpf = chain._state(chain.codes[k], t[k])
        v_dep = vp0 + legs[k - 1][path[k - 1]][0]
        v_arr = vpf + feasible_arrival_vinf(vin, vout, body)
        r = solve_leg(r0, v_dep, rf, v_arr, t[k] - t[k - 1], v_ballistic=vpf + vin, **kw)
        r.update(node=k, planet=chain.codes[k], dv_impulsive_kms=dv_imp / 1e3,
                 leg_days=(t[k] - t[k - 1]) / DAY)
        rows.append(r)
    summary = dict(
        dv_impulsive_kms=sum(r["dv_impulsive_kms"] for r in rows),
        dv_low_thrust_kms=(sum(r["dv_lt_kms"] for r in rows)
                           if all(r["feasible"] for r in rows) else None),
        all_feasible=all(r["feasible"] for r in rows),
        legs_thrusting=len(rows))
    return rows, summary
