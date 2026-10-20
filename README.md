# cyclerdb — a Venus / Earth / Mars cycler catalogue

Search, validate and catalogue **cycler trajectories**: orbits that shuttle a spacecraft
indefinitely between planets using gravity assists, ideally with no propulsive maneuvers.
Built on ESA's [pykep](https://github.com/esa/pykep) (Lambert solver, flyby model, VSOP2013
ephemeris, Sims–Flanagan low-thrust legs) and [pygmo](https://github.com/esa/pygmo2), plus SciPy.
Comes with a browser UI for filtering, plotting and animating every cycler.

MIT licensed: free to use, modify and redistribute, keeping the copyright notice. If you use
the code or the catalogue, please credit the project (see [Citing](#citing)).

## What it does

| Stage | Module | Model | Output table |
|---|---|---|---|
| 1. Single-assist cyclers (Hollister/Byrnes/McConaghy construction), all 6 ordered pairs | `single_home.py` | circular coplanar | `single_home` |
| 2. Multi-assist periodic cyclers for EM, EV, MV and the V-E-M triple | `sequence.py` | circular coplanar | `sequence` |
| 3. Real-planet continuation of **every** cycler, every launch window 2030–2060 | `continuation.py` | VSOP2013 (or DE440) | `ephemeris`, `real_summary` |
| Low-thrust maintenance (experimental) | `lowthrust.py` | Sims–Flanagan | screen column in `ephemeris` |

Everything lands in `data/cyclers.sqlite` (and CSVs). The web UI's Method tab explains the
whole process in plain language; this README is the short version.

## Two models

Stages 1–2 use the **ideal model**: circular, coplanar planet orbits, in which a cycler repeats
exactly. That is how cyclers are defined and named in the literature, and it is the model of
the published tables used for validation. It ignores Mars' eccentricity (0.093) and the planets'
inclinations, so stage 3 re-solves every cycler on the **real planets** (VSOP2013) for every
launch window 2030–2060. The `real_summary` table and the web UI's *Ideal / Real planets*
switch give both numbers side by side; real-planet values are medians over launch windows.

## Validation

`python scripts/validate_mcconaghy.py` reproduces four published results from James Longuski's
group (both McConaghy papers are cited in the references; the PDFs are not redistributed here). Each is checked in the model
the authors used; both McConaghy papers state that their tables assume circular, coplanar
orbits with a synodic period of exactly 2 1/7 years (`bodies.idealized_bodies()`).

| Source | What is compared | Result |
|---|---|---|
| McConaghy, Longuski & Byrnes 2002/2004, **Table 4** | 19 Earth–Mars cyclers: aphelion, V∞ at Earth and Mars, shortest transfer time, required and maximum turn angle, ballistic flag | 19/19, every column within the paper's rounding |
| McConaghy et al. 2002, **Table 6** | 3 two-synodic cyclers with an intermediate Earth flyby: τ, V∞, every Mars-orbit crossing time, V∞ at Mars | 3/3 (times within 0.002 yr) |
| McConaghy et al. 2004, **Table 6** | 30-year ballistic S1L1 itinerary on real ephemerides (paper: DE405; here VSOP2013), 22 encounters | 22/22: V∞ within 0.085 km/s, flyby altitudes median 0.2% (max 7.3%), 45 m/s total residual after moving dates ≤ 1.5 d |
| Byrnes, Longuski & Aldrin 1993 | Aldrin cycler V∞ at Earth/Mars, transfer time; real-ephemeris cost as a comparison | 4/4 |

Two details matter:
- **Footnote c.** For the five Table 4 cyclers whose aphelion is inside Mars' orbit, the paper's
  "V∞ at Mars" is Mars' orbital speed minus the spacecraft's aphelion speed. That is what is
  compared.
- **Footnote d (reconstructed).** Their transfer time is labelled "time to transfer from Earth to
  aphelion", but every printed value is 8–10 days shorter than the exact time to aphelion. The
  time to reach 99.9% of the aphelion radius matches all five within 0.4 days, consistent with a
  radius-based root-finder at 0.1% tolerance. The paper doesn't say so; it is inferred from the
  numbers, and both values are reported.

Other checks (`pytest tests`, add `-m "not slow"` to skip the 20-second real-ephemeris test):
the multi-assist optimiser reproduces the single-assist Aldrin cost exactly at the Aldrin geometry;
launch windows recur once per synodic period; ephemeris states are in the J2000 ecliptic frame;
per-leg revolution counts used for real-planet continuation reproduce the cycler arc itself.

## Catalogue (default settings)

- **Single-assist:** 1,089 solutions across the 6 ordered pairs. 769 cross the target planet's
  orbit and are cyclers; the other 320 are *home-only* orbits (see [Definitions](#definitions)).
  41 of the cyclers are ballistic in the ideal model.
- **Multi-assist:** 177 distinct cyclers, 56 ballistic in the ideal model: EM 6, EV 19, MV 16,
  V-E-M 15.
- **Real planets (2030–2060):** all 946 cyclers re-solved; 21,917 launch windows, 17,680
  converged. 9 three-planet cyclers never line up with the real planets within 15° in that
  period, so they have no windows.
  - **Ideal ballistic ≠ real ballistic.** Of the 97 cyclers that are ballistic in the ideal
    model, 30 keep at least one ballistic window on real planets. Overall 56 cyclers have a
    ballistic window and 94 have a window under 0.1 km/s per cycle (929 such windows).
  - Median real-planet ΔV per cycle (median over each cycler's windows, then over cyclers): EV
    7.7 km/s, MV 11.7, EM 12.1, V-E-M 4.2. Typical v∞ drift from the ideal is 4–8 km/s, because
    the planets' eccentricity and inclination break the ideal geometry.
  - The 15-year Earth–Mars family (7 synodic periods) stays ballistic in nearly every window
    (e.g. 7S8: 13 of 14), though at v∞ about 2 km/s off the ideal.
  - The Aldrin cycler (1L1): median 1.69 km/s per cycle across 14 windows, best window
    ballistic, worst 2.78 km/s. Byrnes et al. (1993) report 247–291 m/s per synodic period with a
    full 15-year optimisation; this catalogue's local ±45-day, powered-flyby-only search is an
    upper bound.

`max_vinf_deviation_kms` in the `ephemeris` table (the "v∞ drift" in the UI) says how far each
window strayed from its ideal family. Treat large values as "a different cycler that happens to
be cheap".

## Run it

Linux / macOS (pykep 3.0.1 wheels are on PyPI for Linux only):

```bash
pip install -r requirements.txt
python scripts/validate_mcconaghy.py
python scripts/build_catalog.py --workers 2 --seeds 2    # stages 1-2 ~25 min on 2 cores;
                                                         # stage 3 (all cyclers) ~80 min on 15
                                                         # (--real-scope top for a quick subset)
python scripts/export_web.py                             # -> web/catalog.js
pytest tests
```

Windows (conda-forge pykep 3.0.0; a self-contained micromamba lives in `.tools/`):

```powershell
$env:MAMBA_ROOT_PREFIX = "$PWD\.tools\mamba-root"
.\.tools\micromamba.exe create -p .\.env -f environment.yml
.\.tools\micromamba.exe run -p .\.env python scripts/validate_mcconaghy.py
.\.tools\micromamba.exe run -p .\.env python scripts/build_catalog.py --workers 15 --seeds 2
.\.tools\micromamba.exe run -p .\.env python scripts/export_web.py
.\.tools\micromamba.exe run -p .\.env python -m pytest tests
```

Use `micromamba run` (or activate the env) rather than calling `.env\python.exe` directly, or
pykep's DLLs won't be found. If micromamba is missing, download it from
`https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-win-64`.

Stages cache to `data/` and can be rerun individually with `--stage 1|2|3`. Stages 2 and 3 are
resumable: they skip anything already in `data/sequence_raw.jsonl` / `data/ephemeris_raw.jsonl`
(move a file away to re-solve from scratch). On Windows, stage 3 stops the machine from sleeping
while it runs.

## Web UI

`web/index.html` is a static, dependency-free browser app. Open it directly from disk (no
server needed) after running `scripts/export_web.py`, which turns `data/cyclers.sqlite` into
`web/catalog.js` plus one small file per cycler in `web/ephem/` (real-planet trajectories,
loaded when you click a window).

- **Catalogue**, with an *Ideal model / Real planets* switch: filter cyclers by family, planet
  pair, ΔV per cycle, v∞, period and launch C3; plot any two of period, ΔV, v∞, aphelion,
  perihelion, C3, ΔV from LEO and v∞ drift; sortable table with CSV export. The detail panel
  compares ideal and real values (median / best / worst window). The trajectory view animates
  the ideal cycler, or in Real-planets mode its best real launch window, in the inertial frame or
  rotating with the home planet. Home-only orbits have their own family button. Single-assist
  cyclers carry their McConaghy notation, which the search box accepts (`1L1`).
- **Real-planet windows**: every launch window 2030–2060, plotted against launch date, with
  encounter dates, v∞ in/out at each encounter and the real trajectory.
- **Validation**: all four reproductions above, paper value beside ours for every number; click
  any row to plot it (Table 4 and 2002 Table 6 in the paper's idealised model, the 2008–2038
  itinerary on real planets).
- **Method**: documentation of the whole process: models, pipeline, notation, optimisation,
  costs, validation, limitations, references.

`#ideal`, `#ideal/real`, `#ephem`, `#valid` and `#method` in the URL open a view directly.

## Definitions

- **ΔV per cycle**: impulsive ΔV to patch arcs at the flybys (pykep `fb_dv`, applied to the
  v∞ vector, which upper-bounds a periapsis burn). Ballistic = ≤ 1 m/s.
- **Safe flyby altitudes**: Venus 300 km, Earth 200 km, Mars 200 km (`bodies.py`).
- **Triple (V-E-M) cyclers** cannot be exactly periodic. They use near-commensurate periods
  (6.41 yr ≈ 3 E-M ≈ 4 E-V ≈ 7 V-M synodic periods; 12.81 yr) and report `drift_deg`, the
  per-cycle phase slip of the non-home planet. Ephemeris continuation absorbs the drift.
- **Single-assist cyclers** (table `single_home`) get exactly one gravity assist per cycle,
  always at the home planet. They visit the target where the arc crosses its orbit; that
  encounter needs no turn, so all cost is at the home flyby. **Multi-assist cyclers** (table
  `sequence`) use a gravity assist at every encounter of a longer pattern such as EEM or EMEV.
  (Earlier versions called these "single-home" and "multi-flyby".)
- **Notation** (McConaghy, Longuski & Byrnes 2002): `nPr`, i.e. *n* synodic periods per cycle,
  **L** (long period: the higher-energy Lambert solution with the larger semi-major axis),
  **S** (short period) or **U** (unique: the single solution with fewer than one revolution),
  then *r* full solar revolutions. The Aldrin cycler is `1L1`; the n = 1 family is 1U0, 1L1,
  1S1, 1L2, 1S2, 1L3, 1S3. Column `notation` in `single_home`; used here for every planet
  pair, not just E–M. (Cycler IDs keep their older form, e.g. `SH-EM-n1r1hi`, `SH-VE-n2r0-`.)
- **Home-only orbits**: single-assist solutions whose arc never crosses the target's circular
  orbit (`reaches_target = 0`). They are genuine periodic orbits of the home planet that the
  enumeration turns up, but they never visit the target, so they are not cyclers. They stay in
  the database for completeness; the UI lists them separately. 40 of them still overlap the
  target's real, eccentric orbit (`reaches_target_eccentric`).
- **Boarding cost**: `earth_vinf_kms` is the cycler's v∞ at an Earth encounter (for multi-assist
  cyclers the lowest one; for ephemeris windows the departure v∞ at the first Earth encounter).
  `earth_c3_km2s2` = v∞² is the launch energy to board it directly, and `leo_departure_dv_kms`
  = √(v∞² + 2μ/r) − √(μ/r) is the burn from a 200 km circular parking orbit. Empty for cyclers
  that never meet Earth. This is a one-time cost, separate from the ΔV per cycle.

## Known limits / next steps

1. **Coplanar ideal model.** Stage 1–2 are 2-D; inclination enters only in stage 3.
2. **Local continuation.** Stage 3 optimises encounter dates within ±45 days of the ideal
   guess (Powell). Global re-search per window would find more ballistic windows.
3. **No deep-space maneuvers.** Only powered flybys. Adding a DSM per leg (MGA-1DSM, pykep has
   the encoding) would lower ephemeris costs, as in Russell & Ocampo.
4. **Low thrust is experimental.** The leg solver converges (homotopy + SLSQP), but replacing a
   powered flyby with thrust on the preceding leg, holding dates fixed, is usually *more*
   expensive than the impulse. A proper low-thrust cycler needs the whole chain re-optimised
   with thrust arcs (MGA-LT). The catalogue reports only a lower-bound SEP acceleration screen.
5. **Resonant / half-revolution "generic return" legs** (planet → same planet after an integer
   or half-integer number of years) are not yet modelled; they matter for VISIT-style cyclers
   (the n = 7 family of McConaghy et al. Table 5, which is degenerate in the Lambert sense).
6. **Rounded dates in drawn real trajectories.** Encounter dates are stored to the day; the web
   UI rebuilds real-planet arcs from them, so drawings can differ slightly from the optimum.
7. **Fixed in this version:** the conversion of single-assist cyclers for real-planet
   continuation used to take per-leg revolution counts from the multi-assist solver, which could
   pick a cheaper but different Lambert arc for the second leg (24 of 74 test cases; up to
   13 km/s off). Revolution counts now come from the arc itself, floor(Δt / P), with a
   regression test. Older stage-3 results for single-assist cyclers should be regenerated.
8. pykep 3.0.1's wheel is missing some data files: four TOPS benchmark JSONs (stub them with
   `{}` in `pykep/trajopt/gym/tops/`) and `de440s.bsp` (hence VSOP2013 by default; set
   `CYCLERDB_DE440=/path/de440s.bsp` and `EphemerisModel("spice")` for DE440).
9. On Windows the only pykep build is conda-forge 3.0.0. It passes the full validation and test
   suite, but its `kep3.dll` needs `openblas.dll`, so `libopenblas` must be installed alongside
   it (included in `environment.yml`).
10. pykep 3.0.1 sometimes aborts at interpreter exit ("corrupted double-linked list"), after
   all work is done. Importing pykep before pandas makes it rarer. The scripts end with
   `os._exit(0)` to avoid it. Results are never affected.

## Prior work this builds on

Hollister (1969); Byrnes, Longuski & Aldrin (1993); McConaghy, Longuski & Byrnes (2002, 2004);
Russell & Ocampo (2004, 2005); Chen et al. (2002, powered cyclers); Barrow & Holzinger
(Mars-Earth-Venus tours); Hernandez et al. (2017, V-E-M triple cyclers).

## Citing

If you use the code or the catalogue, please credit the project, for example:

> Mutasim, S. (2026). *cyclerdb: a Venus / Earth / Mars cycler trajectory catalogue*
> (version 0.1.0).

`CITATION.cff` carries the same information in machine-readable form; GitHub shows it as
"Cite this repository". Please also cite the original cycler papers listed above where
relevant.

## License

Copyright (c) 2026 Sharif Mutasim. Released under the [MIT License](LICENSE): you may use,
copy, modify and distribute this software and the catalogue it produces, provided the
copyright notice and license text are kept in all copies or substantial portions.

Third-party dependencies keep their own licenses and are not bundled: pykep and pygmo
(MPL-2.0), pagmo (LGPL-3.0-or-later / GPL-3.0-or-later), NumPy, SciPy and pandas
(BSD-3-Clause).
