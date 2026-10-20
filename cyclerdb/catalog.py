# Copyright (c) 2026 Sharif Mutasim
# SPDX-License-Identifier: MIT
"""Catalogue schema, search-space enumeration and storage (SQLite + CSV)."""
from __future__ import annotations

import itertools
import sqlite3
from pathlib import Path

import numpy as np

from .bodies import BODIES, YEAR
from .sequence import candidate_periods

SYSTEMS = ["EM", "EV", "MV", "EMV"]   # planet sets, sorted letters


def canonical(seq: str) -> str:
    """Canonical rotation: start at Earth if present, then lexicographically smallest."""
    rots = [seq[i:] + seq[:i] for i in range(len(seq))]
    e_rots = [r for r in rots if r[0] == "E"]
    return min(e_rots or rots)


def sequences(planets: str, max_len: int = 4) -> list[str]:
    """All distinct cyclic encounter sequences that visit every planet in ``planets``."""
    out = set()
    for L in range(len(planets), max_len + 1):
        for tup in itertools.product(planets, repeat=L):
            s = "".join(tup)
            if set(s) != set(planets):
                continue
            out.add(canonical(s))
    return sorted(out, key=lambda s: (len(s), s))


def search_space(max_len: int = 4, pair_max_years: float = 6.5, triple_max_years: float = 13.0,
                 triple_max_drift_deg: float = 10.0):
    """(sequence, Tc_seconds, max_revs) problems for the sequence search."""
    probs = []
    for sysp in SYSTEMS:
        if len(sysp) == 2:
            periods = candidate_periods(sysp, max_years=pair_max_years)
        else:
            periods = candidate_periods(sysp, max_years=triple_max_years,
                                        max_drift_deg=triple_max_drift_deg)
        for seq in sequences(sysp, max_len):
            for Tc in periods:
                max_revs = int(min(5, np.ceil(Tc / YEAR / len(seq)) + 1))
                probs.append((seq, float(Tc), max_revs))
    return probs


def dedupe(rows: list[dict], dv_tol=0.002, vinf_tol=0.05) -> list[dict]:
    """Remove duplicate sequence solutions (same seq/period, ~same ΔV and v∞)."""
    kept: list[dict] = []
    for r in sorted(rows, key=lambda r: r["dv_per_cycle_kms"]):
        v = np.array([float(x) for x in r["vinf_kms"].split("|")])
        dup = False
        for k in kept:
            if k["sequence"] != r["sequence"] or abs(k["period_years"] - r["period_years"]) > 1e-6:
                continue
            vk = np.array([float(x) for x in k["vinf_kms"].split("|")])
            if abs(k["dv_per_cycle_kms"] - r["dv_per_cycle_kms"]) < dv_tol and np.allclose(np.sort(v), np.sort(vk), atol=vinf_tol):
                dup = True
                break
        if not dup:
            kept.append(r)
    return kept


LEO_ALT_KM = 200.0   # circular parking orbit used for the departure ΔV


def leo_departure_dv(vinf_kms, alt_km: float = LEO_ALT_KM):
    """Impulsive ΔV (km/s) from a circular Earth parking orbit onto a hyperbola with v∞."""
    E = BODIES["E"]
    r = E.radius + alt_km * 1e3
    v = np.asarray(vinf_kms, dtype=float) * 1e3
    return (np.sqrt(v**2 + 2 * E.mu / r) - np.sqrt(E.mu / r)) / 1e3


def add_entry_costs(df, kind: str):
    """Add the cost of boarding the cycler from Earth.

    ``earth_vinf_kms`` is the cycler's v∞ at an Earth encounter (for sequence
    cyclers the lowest one, i.e. the cheapest boarding point; for ephemeris
    windows the departure v∞ at the first Earth encounter). From it:
    launch energy ``earth_c3_km2s2`` = v∞² and ``leo_departure_dv_kms``, the
    burn from a 200 km circular orbit. NaN when the cycler never meets Earth.
    """
    df = df.copy()
    if kind == "single_home":
        v = np.where(df["home"] == "E", df["vinf_home_kms"],
                     np.where((df["target"] == "E") & df["reaches_target"].astype(bool),
                              df["vinf_target_kms"], np.nan))
    elif kind == "sequence":
        v = [min((float(x) for c, x in zip(s, str(vi).split("|")) if c == "E"), default=np.nan)
             for s, vi in zip(df["sequence"], df["vinf_kms"])]
    elif kind == "ephemeris":
        v = df["earth_vinf_kms"] if "earth_vinf_kms" in df else np.full(len(df), np.nan)
    else:
        raise ValueError(kind)
    v = np.asarray(v, dtype=float)
    df["earth_vinf_kms"] = v
    df["earth_c3_km2s2"] = v**2
    df["leo_departure_dv_kms"] = leo_departure_dv(v)
    return df


def _max_vinf(s) -> float:
    vals = [float(kv.split(":")[1]) for kv in str(s or "").split("|") if ":" in kv]
    return max(vals) if vals else np.nan


def real_summary(ep):
    """One row per cycler summarising its real-ephemeris launch windows.

    Statistics are over converged windows only. ``real_ballistic_windows``
    counts windows at ≤ 1 m/s per cycle (the same threshold as the ideal
    ``ballistic`` flag).
    """
    import pandas as pd
    ep = ep.copy()
    ep["max_vinf_any_kms"] = ep["max_vinf_kms"].map(_max_vinf)
    rows = []
    for label, g in ep.groupby("label"):
        ok = g[g["converged"].astype(bool) & np.isfinite(g["dv_per_cycle_kms"].astype(float))]
        dv = ok["dv_per_cycle_kms"].astype(float)
        c3 = ok["earth_c3_km2s2"].astype(float)
        rows.append(dict(
            label=label, sequence=g["sequence"].iloc[0], n_windows=len(g), n_converged=len(ok),
            real_dv_best_kms=dv.min() if len(ok) else np.nan,
            real_dv_median_kms=dv.median() if len(ok) else np.nan,
            real_dv_worst_kms=dv.max() if len(ok) else np.nan,
            real_ballistic_windows=int((dv <= 1e-3).sum()),
            real_cheap_windows=int((dv <= 0.1).sum()),
            real_vinf_max_median_kms=ok["max_vinf_any_kms"].median() if len(ok) else np.nan,
            real_vinf_deviation_median_kms=ok["max_vinf_deviation_kms"].astype(float).median() if len(ok) else np.nan,
            real_c3_best_km2s2=c3.min() if c3.notna().any() else np.nan,
            real_c3_median_km2s2=c3.median() if c3.notna().any() else np.nan,
            real_c3_highest_km2s2=c3.max() if c3.notna().any() else np.nan,
            real_leo_dv_median_kms=ok["leo_departure_dv_kms"].astype(float).median() if c3.notna().any() else np.nan,
            first_window=ok["start_epoch"].min() if len(ok) else None,
        ))
    return pd.DataFrame(rows)


def write_sqlite(path: str | Path, tables: dict[str, "pandas.DataFrame"]):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(path) as con:
        for name, df in tables.items():
            df.to_sql(name, con, if_exists="replace", index=False)
