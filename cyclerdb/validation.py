# Copyright (c) 2026 Sharif Mutasim
# SPDX-License-Identifier: MIT
"""Reference values from the literature used to validate the code.

Sources (both in the repository root as PDFs):
  [M02] McConaghy, Longuski & Byrnes, "Analysis of a Broad Class of Earth-Mars
        Cycler Trajectories", AIAA 2002-4420 (2002).
  [M04] McConaghy, Longuski & Byrnes, "Analysis of a Class of Earth-Mars Cycler
        Trajectories", J. Spacecraft & Rockets 41(4), 622-628 (2004).
  [B93] Byrnes, Longuski & Aldrin, "Cycler Orbit Between Earth and Mars",
        J. Spacecraft & Rockets 30(3), 334-336 (1993), as summarised in [M04].

Model used by [M02]/[M04] for their tables (stated in both papers): Earth and
Mars on circular, coplanar orbits; S = 2 1/7 yr (Mars period 1 7/8 yr);
conic, prograde trajectory; gravity assists at Earth only, instantaneous;
minimum Earth flyby altitude 200 km. ``bodies.idealized_bodies()``.
"""

# [M02] Table 4 = [M04] Table 4, "The most promising cyclers that repeat every
# one to six synodic periods". Columns: notation, aphelion (AU), V∞ Earth,
# V∞ Mars (km/s), shortest transfer time (days), required turn angle (deg),
# maximum possible turn angle (deg), footnote flags.
# Footnotes: b = aphelion below Mars' semimajor axis (1.52 AU);
#            c = V∞ at Mars is "difference between Mars' speed and spacecraft
#                aphelion speed" (no actual crossing);
#            d = transfer time is "time to transfer from Earth to aphelion";
#            e = ballistic (required turn < maximum turn).
TABLE4 = [
    ("1L1", 2.23, 6.54, 9.75, 146, 84, 72, ""),
    ("2L2", 2.33, 10.06, 11.27, 158, 134, 44, ""),
    ("2L3", 1.51, 5.65, 3.05, 280, 135, 82, "bcd"),
    ("3L4", 1.89, 11.78, 9.68, 189, 167, 35, ""),
    ("3L5", 1.45, 7.61, 2.97, 274, 167, 62, "bcd"),
    ("3S5", 1.52, 12.27, 5.45, 134, 167, 33, "bcd"),
    ("4S5", 1.82, 11.23, 8.89, 88, 167, 38, ""),
    ("4S6", 1.53, 8.51, 4.07, 157, 167, 54, ""),
    ("5S4", 2.49, 10.62, 12.05, 75, 134, 41, ""),
    ("5S5", 2.09, 9.08, 9.87, 89, 134, 50, ""),
    ("5S6", 1.79, 7.51, 7.32, 111, 135, 62, ""),
    ("5S7", 1.54, 5.86, 3.67, 170, 135, 79, ""),
    ("5S8", 1.34, 4.11, 0.71, 167, 136, 103, "bcd"),
    ("6S4", 2.81, 7.93, 12.05, 87, 83, 59, ""),
    ("6S5", 2.37, 6.94, 10.44, 97, 84, 68, ""),
    ("6S6", 2.04, 5.96, 8.69, 111, 84, 78, ""),
    ("6S7", 1.78, 4.99, 6.66, 133, 85, 90, "e"),
    ("6S8", 1.57, 4.02, 3.90, 179, 85, 104, "e"),
    ("6S9", 1.40, 3.04, 1.21, 203, 86, 120, "bcde"),
]
TABLE4_NOTES = {
    "1L1": "Aldrin cycler",
    "2L3": "Case 1 two-synodic cycler of Byrnes et al. (2002)",
}

# [M02] Table 6: ballistic two-synodic-period cyclers with one intermediate
# Earth flyby at time tau (years from the start of the 2S repeat interval).
# (legs P1r1, P2r2, tau, V∞ Earth, leg-1 Mars-orbit crossing times (yr),
#  leg-1 V∞ at crossing, leg-2 crossing times, leg-2 V∞); None = no crossing.
TABLE6_2002 = [
    ("U0", "L1", 2.754, 11.3, (0.188, 2.567), 14.0, (3.449, 3.592), 5.4),
    ("L2", "U0", 2.541, 8.8, None, None, (2.781, 4.046), 10.3),
    ("S1", "L1", 2.828, 4.7, (0.419, 0.920, 1.908, 2.409), 5.0, None, None),
]

# [M04] Table 6: outbound ballistic S1L1 cycler itinerary using JPL DE405
# ephemerides (a real, eccentric and inclined solar system).
# (body, date, approach V∞ km/s, closest-approach altitude km, next-leg days)
S1L1_DE405 = [
    ("E", "2008-06-09", 6.89, None, 541),   # launch
    ("E", "2009-12-03", 6.90, 31114, 186),
    ("M", "2010-06-06", 4.31, 17704, 809),
    ("E", "2012-08-24", 6.42, 26490, 540),
    ("E", "2014-02-14", 6.43, 41524, 139),
    ("M", "2014-07-03", 7.14, 12179, 890),
    ("E", "2016-12-09", 4.01, 27726, 530),
    ("E", "2018-05-22", 4.03, 19923, 115),
    ("M", "2018-09-15", 6.47, 11570, 934),
    ("E", "2021-04-06", 4.61, 22992, 532),
    ("E", "2022-09-20", 4.59, 14780, 223),
    ("M", "2023-05-01", 2.77, 7593, 793),
    ("E", "2025-07-02", 7.08, 23858, 542),
    ("E", "2026-12-26", 7.09, 35164, 170),
    ("M", "2027-06-14", 5.26, 13751, 830),
    ("E", "2029-09-21", 5.78, 26818, 537),
    ("E", "2031-03-12", 5.78, 39044, 125),
    ("M", "2031-07-15", 7.70, 10566, 915),
    ("E", "2034-01-15", 3.78, 22988, 529),
    ("E", "2035-06-28", 3.76, 9586, 138),
    ("M", "2035-11-13", 4.68, 15525, 907),
    ("E", "2038-05-07", 5.55, None, None),
]

# [B93] via [M04]: Aldrin cycler optimised over 15 years in an accurate solar
# system model (dates not given in [M04]).
ALDRIN_REAL_B93 = dict(outbound_total_dv_kms=1.73, outbound_dv_per_synodic_ms=247,
                       inbound_total_dv_kms=2.04, inbound_dv_per_synodic_ms=291,
                       transfer_days=(147, 170), vinf_earth_kms=(5.39, 6.19),
                       vinf_mars_kms=(6.05, 11.74))

# Kept for the original scalar checks.
ALDRIN_1993 = dict(vinf_earth_kms=6.54, vinf_mars_kms=9.75, aphelion_au=2.23, transfer_days=146)
NOTABLE_2SYN = dict(n_synodic=2, n_revs=3, vinf_earth_kms=5.65, vinf_mars_kms=3.05)
