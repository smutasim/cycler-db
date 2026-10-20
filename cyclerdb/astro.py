# Copyright (c) 2026 Sharif Mutasim
# SPDX-License-Identifier: MIT
"""Two-body helpers built on pykep: Lambert solutions, flybys, conic geometry."""
from __future__ import annotations

import numpy as np
import pykep as pk

from .bodies import MU_SUN, Body


def lambert_all(r0, r1, tof, max_revs, cw=False, mu=MU_SUN):
    """All Lambert solutions up to ``max_revs`` full revolutions.

    Returns a list of (n_revs, branch, v0, v1) where branch is 0 for the
    single zero-rev solution, otherwise 'L' (long period / low energy) or 'R'
    (short period) following pykep's ordering. Empty list on failure.
    """
    try:
        lp = pk.lambert_problem(r0=list(r0), r1=list(r1), tof=float(tof), mu=mu,
                                cw=cw, multi_revs=int(max_revs))
    except (ValueError, RuntimeError):
        return []
    v0s, v1s = lp.v0, lp.v1
    out = [(0, "0", np.asarray(v0s[0]), np.asarray(v1s[0]))]
    for k in range(1, len(v0s)):
        n = (k + 1) // 2
        branch = "L" if k % 2 == 1 else "R"
        out.append((n, branch, np.asarray(v0s[k]), np.asarray(v1s[k])))
    return out


def max_turn(vinf: float, body: Body) -> float:
    """Largest unpowered turn angle (rad) for a flyby at the safe radius."""
    return 2.0 * np.arcsin(1.0 / (1.0 + body.safe_radius * vinf**2 / body.mu))


def angle_between(a, b) -> float:
    c = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.arccos(np.clip(c, -1.0, 1.0)))


def flyby_dv(vinf_in, vinf_out, body: Body) -> float:
    """Impulsive ΔV (m/s) to patch in/out v-infinity at a flyby (pykep, Mivovitch model)."""
    return float(pk.fb_dv(list(vinf_in), list(vinf_out), body.mu, body.safe_radius))


def conic(r, v, mu=MU_SUN):
    """Return (a, e, rp, ra, period) of the heliocentric conic."""
    r = np.asarray(r); v = np.asarray(v)
    rn, vn = np.linalg.norm(r), np.linalg.norm(v)
    energy = vn**2 / 2 - mu / rn
    a = -mu / (2 * energy)
    h = np.cross(r, v)
    evec = np.cross(v, h) / mu - r / rn
    e = float(np.linalg.norm(evec))
    rp = a * (1 - e)
    ra = a * (1 + e) if e < 1 else np.inf
    period = 2 * np.pi * np.sqrt(a**3 / mu) if a > 0 else np.inf
    return a, e, rp, ra, period


def true_anomaly(r, v, mu=MU_SUN) -> float:
    r = np.asarray(r); v = np.asarray(v)
    rn = np.linalg.norm(r)
    h = np.cross(r, v)
    evec = np.cross(v, h) / mu - r / rn
    e = np.linalg.norm(evec)
    if e < 1e-12:
        return 0.0
    nu = np.arccos(np.clip(np.dot(evec, r) / (e * rn), -1, 1))
    if np.dot(r, v) < 0:
        nu = 2 * np.pi - nu
    return float(nu)


def time_since_periapsis(nu: float, a: float, e: float, mu=MU_SUN) -> float:
    """Time from periapsis to true anomaly nu in [0, 2π) on an ellipse."""
    E = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu / 2), np.sqrt(1 + e) * np.cos(nu / 2))
    M = E - e * np.sin(E)
    M = M % (2 * np.pi)
    return M / np.sqrt(mu / a**3)


def radius_crossings(r0, v0, radius: float, mu=MU_SUN):
    """Times after the state (r0, v0) at which the conic first reaches ``radius``.

    Returns (t_outbound, t_inbound) for the first crossing with radial velocity
    > 0 and < 0 respectively, or (None, None) if the orbit never reaches it.
    Also returns the heliocentric speed and flight path angle there.
    """
    a, e, rp, ra, P = conic(r0, v0, mu)
    if not (rp <= radius <= ra) or a <= 0:
        return None
    p = a * (1 - e**2)
    cosnu = (p / radius - 1) / e
    nu_out = float(np.arccos(np.clip(cosnu, -1, 1)))  # radial velocity > 0
    nu_in = 2 * np.pi - nu_out
    nu0 = true_anomaly(r0, v0, mu)
    t0 = time_since_periapsis(nu0, a, e, mu)
    t_out = (time_since_periapsis(nu_out, a, e, mu) - t0) % P
    t_in = (time_since_periapsis(nu_in, a, e, mu) - t0) % P
    speed = np.sqrt(mu * (2 / radius - 1 / a))
    h = np.sqrt(mu * p)
    gamma = np.arccos(np.clip(h / (radius * speed), -1, 1))  # flight path angle magnitude
    return dict(t_out=t_out, t_in=t_in, speed=speed, gamma=gamma, period=P, a=a, e=e, rp=rp, ra=ra)


def rot_z(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
