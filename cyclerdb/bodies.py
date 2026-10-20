# Copyright (c) 2026 Sharif Mutasim
# SPDX-License-Identifier: MIT
"""Planet constants and ephemeris back-ends.

Two models are provided:

* ``CircularModel`` – circular, coplanar planetary orbits. This is the model in
  which the classic cycler literature (Hollister, Byrnes, McConaghy, Russell &
  Ocampo) defines and catalogues cyclers, and in which a cycler is exactly
  periodic.
* ``EphemerisModel`` – real planetary positions. Uses pykep's VSOP2013
  analytic theory by default (works offline, arcsecond-level accuracy), or a
  JPL DE440 SPICE kernel when one is available on disk.

All quantities are SI (m, s, m^3/s^2) unless a name says otherwise.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import os

import numpy as np
import pykep as pk

MU_SUN = pk.MU_SUN
AU = pk.AU
DAY = pk.DAY2SEC
YEAR = 2.0 * np.pi * np.sqrt(AU**3 / MU_SUN)  # one Earth sidereal year, s


@dataclass(frozen=True)
class Body:
    code: str          # one-letter code used in sequences: V, E, M
    name: str
    mu: float          # m^3/s^2
    radius: float      # m
    min_alt: float     # m, minimum flyby altitude
    a_au: float        # mean semi-major axis, AU (circular model radius)
    vsop: str          # pykep vsop2013 body name
    spice: str         # SPICE body name

    @property
    def safe_radius(self) -> float:
        return self.radius + self.min_alt

    @property
    def a(self) -> float:
        return self.a_au * AU

    @property
    def period(self) -> float:
        return 2.0 * np.pi * np.sqrt(self.a**3 / MU_SUN)

    @property
    def v_circ(self) -> float:
        return np.sqrt(MU_SUN / self.a)


# Mean semi-major axes from Standish (JPL "approximate positions", J2000).
# Minimum flyby altitudes: 300 km Venus (atmosphere), 200 km Earth and Mars,
# matching common practice in the cycler literature. All are overridable.
BODIES: dict[str, Body] = {
    "V": Body("V", "Venus", 3.24858592e14, 6051.8e3, 300e3, 0.72333566, "venus", "VENUS BARYCENTER"),
    "E": Body("E", "Earth", 3.98600436e14, 6378.137e3, 200e3, 1.00000261, "earth_moon", "EARTH BARYCENTER"),
    "M": Body("M", "Mars", 4.28283752e13, 3396.19e3, 200e3, 1.52371034, "mars", "MARS BARYCENTER"),
}


def idealized_bodies(mars_period_years: float = 15 / 8) -> dict[str, Body]:
    """Bodies for the classic idealised Earth–Mars model used by Byrnes,
    McConaghy et al.: Earth on a 1 AU circle and a Mars period of exactly 15/8
    years, so 7 synodic periods = 15 years. Venus is left at its mean value."""
    out = dict(BODIES)
    e, m = BODIES["E"], BODIES["M"]
    out["E"] = Body("E", e.name, e.mu, e.radius, e.min_alt, 1.0, e.vsop, e.spice)
    out["M"] = Body("M", m.name, m.mu, m.radius, m.min_alt,
                    mars_period_years ** (2 / 3), m.vsop, m.spice)
    return out


# Real orbital eccentricity of each planet: used only to flag cyclers whose
# arc misses the circular orbit but would still reach the real (eccentric) one.
ECCENTRICITY = {"V": 0.00677672, "E": 0.01671123, "M": 0.09339410}


def synodic_period(b1: Body, b2: Body) -> float:
    return 1.0 / abs(1.0 / b1.period - 1.0 / b2.period)


def with_min_alt(code: str, min_alt_m: float) -> Body:
    """Return a copy of a body with a different minimum flyby altitude."""
    b = BODIES[code]
    return Body(b.code, b.name, b.mu, b.radius, min_alt_m, b.a_au, b.vsop, b.spice)


# ---------------------------------------------------------------------------
# Circular coplanar model
# ---------------------------------------------------------------------------
class CircularModel:
    """Planets on circular coplanar orbits, longitude = phase + n * t."""

    def __init__(self, phases: dict[str, float] | None = None, bodies: dict | None = None):
        self.phases = phases or {}
        self.bodies = bodies or BODIES

    def state(self, code: str, t: float, phase: float | None = None):
        b = self.bodies[code]
        th = (self.phases.get(code, 0.0) if phase is None else phase) + 2 * np.pi * t / b.period
        r = b.a * np.array([np.cos(th), np.sin(th), 0.0])
        v = b.v_circ * np.array([-np.sin(th), np.cos(th), 0.0])
        return r, v


# ---------------------------------------------------------------------------
# Ephemeris model
# ---------------------------------------------------------------------------
@lru_cache(maxsize=None)
def _vsop(code: str):
    return pk.planet(pk.udpla.vsop2013(BODIES[code].vsop))


@lru_cache(maxsize=None)
def _spice(code: str, kernel: str):
    return pk.planet(pk.udpla.spice(kernel, BODIES[code].spice, "ECLIPJ2000", "SSB"))


_EPS = np.radians(23.4392911)  # J2000 obliquity
_ICRF_TO_ECLIP = np.array([[1, 0, 0],
                           [0, np.cos(_EPS), np.sin(_EPS)],
                           [0, -np.sin(_EPS), np.cos(_EPS)]])


class EphemerisModel:
    """Real planet positions in the J2000 ecliptic frame, heliocentric.

    ``t`` is seconds past J2000 (TDB). VSOP2013 states come out of pykep in the
    ICRF (equatorial) frame and are rotated to the ecliptic so that z is out of
    the ecliptic plane and longitudes are ecliptic longitudes.
    """

    def __init__(self, backend: str = "vsop2013", kernel: str | None = None):
        self.backend = backend
        self.kernel = kernel or os.environ.get("CYCLERDB_DE440")
        if backend == "spice" and not (self.kernel and os.path.exists(self.kernel)):
            raise FileNotFoundError(
                "SPICE backend needs a DE440 kernel; set CYCLERDB_DE440 or pass kernel=")

    def state(self, code: str, t: float):
        ep = pk.epoch(t / DAY)  # pykep epochs are MJD2000 days
        if self.backend == "spice":
            r, v = _spice(code, self.kernel).eph(ep)
            return np.asarray(r), np.asarray(v)
        r, v = _vsop(code).eph(ep)
        return _ICRF_TO_ECLIP @ np.asarray(r), _ICRF_TO_ECLIP @ np.asarray(v)

    def longitude(self, code: str, t: float) -> float:
        r, _ = self.state(code, t)
        return float(np.arctan2(r[1], r[0]))


def epoch_seconds(iso: str) -> float:
    """ISO date -> seconds past J2000."""
    return pk.epoch(iso).mjd2000 * DAY
