# B Current wetland optimization and dynamic comparisons

[Project overview](../README.md) · [A Initial model and baseline](../A/README.md)

B is the current project: wetland design optimization under fixed gate operations, together with the no-wetland dynamic comparison. This page presents the final accepted model, method, results, figures and files for B.

The final accepted analysis asks how wetland design affects the downstream river peak when gate operation is fixed. It includes an eight-state wetland design search and a six-state no-wetland dynamic comparison. Both use the no-dam-spillway scope, the same prescribed gate openings and the approved regional event-weather reconstruction.

## Model and study conditions

The two wetland branches follow:

```text
Upland -> Dam pool 1 -> Gate 1 -> Wetland 1 -> River
Upland -> Dam pool 2 -> Gate 2 -> Wetland 2 -> River
```

The model retains soil and groundwater pathways and direct upland-to-river runoff. Dam pools have no independent spillway path; wetland normal drainage and wetland overflow remain represented. Water storage is not clipped to enforce a safety limit.

The installation footprint is an assumed 0.66 km². The event covers 2022-09-09 00:00 through 2022-09-15 00:00 UTC. Regional station observations, solar grid estimates, spatial transfers, trace-zero lower bounds and assumed boundary/initial conditions retain their separate evidence labels. USGS 04087088 is contextual evidence only. The [source and assumption document](data_sources_and_assumptions.md) defines the accepted scope and its limitations.

## Data, parameter origins and initial conditions

The [complete source and assumption record](data_sources_and_assumptions.md#active-parameter-provenance) lists active values, units, derivations, executable locations and applicability. An implemented number is not automatically an observed or literature-derived value. The evidence labels distinguish observations, values derived from observations, external physical/hydraulic support, design assumptions, numerical settings and unresolved empirical support. Simulation outputs have their own result-table sources.

| Input or quantity | Value or series used | Source and transformation | Evidence / limitation |
|---|---|---|---|
| Precipitation `P(t)` | 145.3 mm quantified event total | KMKE / NOAA ISD 72640014839; accepted `AA1–AA4` accumulations converted from mm per interval to m/s | Regional station observations spatially transferred; five event trace intervals use a declared zero lower bound |
| Air temperature `T_a(t)` | KMWC event history, degrees C | NOAA ISD 72640594869 `TMP`, accepted reports interpolated in UTC | Observed regional proxy; not site temperature |
| Wind `v(t)` | KMWC event history, m/s | Same station, `WND` speed and quality fields | Observed regional proxy; no unsupported height correction |
| Shortwave `R_s(t)` | 21,166.81 Wh/m² event energy | NASA POWER `ALLSKY_SFC_SW_DWN`, SYN1DEG; hourly energy conservatively reconstructed as W/m² | Gridded estimate; label-as-hour-start is the approved processing assumption |
| Humidity `r_H` | 0.7084916611553665 | Event-time mean of RH derived from KMWC `TMP` and `DEW` using the retained saturation-pressure relation | Derived from observations, then fixed for this study |
| Pressure `p_0` | 98,692.51171875 Pa | Event-time mean of KMWC `MA1` station pressure | Derived from observations; not sea-level pressure |
| Context discharge/level | USGS 04087088 `00060`, `00065`, `63160` | Converted from agency units, with native datum and approval flags retained | Context only; not forcing or calibration data for this installation |
| External surface inflow | 0 m³/s | Approved study boundary: no additional surface contributing area | Explicit assumption; USGS discharge is not added |
| External channel/groundwater levels | Both 97.15 m | Assumed model datum | Explicit assumptions, not observed NAVD88 elevations |
| Outlet blockage | 0 | Assumed maintained outlet | No local blockage observation |

Weather rows are documented individually in the [weather record](docs/05_weather_observations.md), including exact queries, dates, fields, QC, trace handling and hashes. Stream records are in the [USGS record](docs/04_usgs_observations.md). The prepared weather has 4,321 knots at 120 s. The same knots are used for every run and refinement; precipitation volume and solar energy are preserved.

| State | Initial storage (m³) | Derivation and origin |
|---|---:|---|
| Upland `S_u` | 600 | Assumed 0.001 m depth × 600,000 m² |
| Soil `S_m` | 90,000 | Assumed 0.75 × 120,000 m³ capacity; not reconstructed September 2022 soil wetness |
| Wetland 1 `S_w1` | 3,000 | Assumed 0.25 m × 12,000 m² |
| Wetland 2 `S_w2` | 2,500 | Assumed 0.25 m × 10,000 m² |
| Pool 1 `S_d1` | 3,000 | Assumed 0.30 m × 10,000 m² |
| Pool 2 `S_d2` | 2,400 | Assumed 0.30 m × 8,000 m² |
| Groundwater `S_g` | 409,600 | Assumed `(97.20 − 94.00) m × 0.20 × 640,000 m²` |
| River `S_r` | 7,000 | Assumed 0.35 m × 20,000 m² |

**Source:** inherited initial-condition definitions in [A class c](../A/wetland_flood_model.py), retained explicitly in [no_spillway_config.json](no_spillway_config.json), `initial_conditions.x0_m3`. The common eight-state initial total is 518,100 m³. Arithmetic derived from an assumed depth remains an assumption, not a local measurement.

The source record also lists every soil, groundwater, outlet, wetland-drainage, overflow, ET and availability coefficient. Their local empirical support is missing unless explicitly stated; no SSURGO, terrain, wetland inventory or plant database was used to manufacture coefficient values. The archived source catalog/evidence CSV retains historical proposal statuses and obsolete spillway/linear-gate entries, so it must not replace the current active parameter record.

## Optimization problem and numerical method

The submitted B.1 formulation defines the peak-stage objective and includes gate histories together with wetland designs. The accepted analysis here fixes those gate histories and optimizes only wetland designs, as explicitly authorized.

The wetland designs are `Di=(tau_i,rho_i)`. The synthetic categories A/A, A/B, B/A and B/B are enumerated, and the continuous decision vector contains only the two cover fractions `[rho1,rho2]`, each in [0,1]. Categories and covers determine the paired hydraulic-resistance and ET response. Wetland geometry remains fixed, and the selected design is constant during each simulation.

```text
Fixed inputs: g1(t), g2(t), weather, shared geometry, boundaries and initial conditions
Decisions: D1 and D2
Objective: minimize the maximum simulated river stage Hr(t) over the event
```

The official multivarious `ode4u` solver integrates each candidate, and the project-selected official `nms` routine searches the covers. Four category pairs and three initial cover guesses give 12 searches. Seven met the convergence conditions and five reached the evaluation budget. The reported best feasible design is the best found in these bounded searches, not a proven global optimum.

Nonnegative storage, soil capacity, fixed upstream safety levels and shared applicability limits remain enforced. The four maximum allowable levels are 101.20/101.00 m for the pools and 98.90/98.80 m for the wetlands, all in the assumed model datum. G5 is checked during RHS/RK evaluations as hydraulic applicability for the prescribed gate schedule; it is excluded from the wetland-only OPT residual vector. The remaining 18 residuals and native cover bounds retain their applicable roles. Physical river peaks are reported separately from numerical penalty values.

## Accepted wetland designs

| Design | Synthetic categories | Cover 1 | Cover 2 | River peak at 30 s |
|---|---|---:|---:|---:|
| Reference wetlands | A/A | 0.65 | 0.60 | 97.558443 m |
| Best-found wetland design | B/B | 0.2714117624187232 | 0.2141279193188096 | 97.553431 m |

The best-found design lowers the reference peak by **5.012 mm** under identical prescribed gate operation. Its lower effective drainage resistance accompanies more normal drainage and less concentrated wetland overflow. The result does not imply that more vegetation cover or greater stored volume always reduces the downstream peak. A/B are synthetic response classes, not identified plant species or local vegetation recommendations.

**Result provenance:** covers and search termination come from [optimization_summary.csv](tables/wetland_only/optimization_summary.csv); 30 s peak values come from [fixed_design_timesteps.csv](tables/wetland_only/fixed_design_timesteps.csv). They are computed model outputs, not externally observed water levels.

The [full optimization report](docs/12_wetland_only_optimization.md) provides the search results, margins, mechanisms, solver termination and numerical-resolution checks.

## Dynamic no-wetland comparison

The six-state counterfactual removes the wetland compartments and routes each gate directly into the river. Rainfall over the former wetland footprints enters the river immediately; the river's storage-stage and evaporation areas remain unchanged. Gate discharge is calculated from internal river tailwater rather than replayed from a saved wetland case. This is a dynamic bypass assumption, not a calibrated replacement land-cover scenario.

| Fixed case | Dynamic states | River peak at 30 s | Bypass peak minus case peak |
|---|---:|---:|---:|
| No-wetland dynamic bypass | 6 | 97.637815 m | 0 mm |
| Reference wetlands | 8 | 97.558443 m | 79.373 mm |
| Best-found wetland design | 8 | 97.553431 m | 84.385 mm |

**Result provenance:** the table above uses [comparison.csv](tables/wetland_counterfactual/comparison.csv), `timestep_s=30`, primary cases `bypass`, `reference`, `best`; the [readable table](tables/wetland_counterfactual/comparison.md) retains additional water-balance and safety values.

![Simulated river-stage comparison](plots/wetland_counterfactual/01_simulated_river_stage.png)

This fixed comparison does not run OPT. The 79.373/84.385 mm differences include the represented wetland storage, routing, seepage and ET mechanisms. Only the 5.012 mm reference-to-best difference is the additional wetland-design optimization improvement.

The primary wetland cases contain 5,500 m³ more initial water because their wetland compartments retain their approved initial storages. A separate diagnostic sets those two storages to zero without changing the shared initial states. It changes the reference/best peaks by approximately 0.2875/0.2217 mm and preserves the qualitative comparison.

## Storage, drainage and overflow mechanism

| Computed quantity at 30 s | Reference | Best-found wetland design |
|---|---:|---:|
| Gate 1 released volume (m³) | 33,255.472 | 33,255.472 |
| Gate 2 released volume (m³) | 27,224.010 | 27,224.010 |
| Peak wetland storage, branch 1 / 2 (m³) | 8,226.768 / 6,984.349 | 7,558.433 / 6,358.098 |
| Total normal wetland drainage (m³) | 57,041.508 | 59,587.321 |
| Total wetland overflow (m³) | 2,194.918 | 235.449 |
| Total wetland-to-river release centroid (h) | 79.267520 | 78.203764 |
| Physical river peak time (h) | 80.041667 | 80.050000 |

**Source:** [mechanism_metrics.csv](tables/wetland_only/mechanism_metrics.csv), rows `fixed_reference` and `wetland_best`, both `dt_s=30`. Volumes integrate the saved flow histories; a release centroid is a flow-weighted mean time, not water-particle residence time. Internal transfers are not added again to the whole-system external water balance.

Normal drainage increases and concentrated overflow decreases under the same gates. The best design has smaller peak wetland storage and an earlier wetland-release centroid, yet a slightly later, lower river peak. These are different diagnostics; the result is not a claim that maximizing storage or delay always minimizes the downstream peak.

![Wetland storage and routing mechanisms](plots/wetland_only/01_wetland_design_mechanisms.png)

## Verification and interpretation

The fixed configurations were checked at 120, 60 and 30 seconds using the same underlying forcing and gate functions. All reported comparison cases pass the applicable feasibility and hydraulic checks. The largest 30-second integrated water-balance residual is approximately 0.00631 m³. Actual gate releases agree within floating-point precision. The software suite passes 100 tests, including four focused contracts for the standalone analysis interface. The fixed reference, best-found wetland design and six-state bypass have been reproduced at 120 seconds with exact agreement in their saved time, state, input, output and balance arrays. Calling the public analysis function also exactly reproduces the reference/best objective and all 18 optimizer residuals. Saved numerical files are not overwritten.

The event maxima occur inside the selected window, but water remains stored and the river rises slightly at the endpoint. These conditional, uncalibrated simulations do not establish complete recession, measured Milwaukee flood reductions, engineering certification or universal wetland benefits. The [full dynamic comparison report](docs/13_no_wetland_counterfactual.md) provides the initialization sensitivity, conservation accounting and limits.

## Direct analysis and design bounds

[The analysis guide](docs/analysis_and_boundary_checks.md) provides the standalone `analysis(v,c)` interface. For one fixed category pair, `v=[rho1,rho2]`; the constants contain the existing gate schedule, forcing, initial conditions and fixed parameters. The function returns the unchanged dimensionless conditioned objective and an 18-element constraint vector. Its diagnostic record separately exposes the physical river peak, validity and feasibility. Both direct calls and the optimizer use this same implementation.

```python
from B.wetland_only_study import make_analysis_constants, analysis
c = make_analysis_constants(("type_a", "type_b"))
f, g = analysis([0.0, 1.0], c)
print(c.last_evaluation["physical_peak_m"])
```

The [design-bound table](tables/analysis_checks/boundary_checks.md) contains 32 actual dynamic evaluations: four coverage corners for all four category pairs, each at 120 and 60 seconds. All 32 complete and satisfy the unchanged applicability and feasibility checks. The maximum coarse/fine peak difference is 0.010340 mm. These checks use 34 integrations including two saved-case reproductions; they do not run OPT or replace the existing best design. Endpoint tests do not prove stability or feasibility throughout the interior of the design box.

Boundary verification figures are in [plots/analysis_checks/](plots/analysis_checks/):

- [River-stage time histories](plots/analysis_checks/01_boundary_river_peaks.png): four coverage-corner curves per class pair, with peak markers, from saved 60 s trajectories.
- [Safety and gate applicability](plots/analysis_checks/02_boundary_safety_margins.png): minimum headroom to the four upstream limits through time for each corner, plus the two shared G5 histories.
- [Numerical checks](plots/analysis_checks/03_boundary_numerical_checks.png): envelopes of matched-time 120/60 s stage differences and cumulative water-balance residuals; sampled peak differences remain a separate statistic.

![River-stage histories for the coverage-bound cases](plots/analysis_checks/01_boundary_river_peaks.png)

All three figures are regenerated from saved checks by `python3 -B -m B.wetland_analysis_checks --report`; this command runs no simulation or optimization.

The [assumption table](data_sources_and_assumptions.md) states the actual resistance/ET endpoints, their units and interpolation, and the basis for full coverage bounds [0,1]. Those coefficients remain explicit synthetic assumptions, not locally measured vegetation properties.

## Final files and reproduction

| File or directory | Purpose |
|---|---|
| [no_spillway_model.py](no_spillway_model.py) and [configuration](no_spillway_config.json) | Accepted eight-state equations and fixed conditions |
| [wetland_only_study.py](wetland_only_study.py) and [configuration](wetland_only_config.json) | Fixed-gate wetland optimization |
| [wetland_counterfactual.py](wetland_counterfactual.py) | Six-state dynamic bypass and fixed comparisons |
| [wetland_analysis_checks.py](wetland_analysis_checks.py) and [analysis guide](docs/analysis_and_boundary_checks.md) | Standalone analysis verification and full coverage-bound dynamic checks |
| [Optimization report](docs/12_wetland_only_optimization.md) | Full wetland design analysis |
| [Comparison report](docs/13_no_wetland_counterfactual.md) | Full three-case analysis |
| [Optimization table](tables/wetland_only/optimization_summary.csv) | Search designs, physical peaks and termination |
| [Comparison table](tables/wetland_counterfactual/comparison.md) and [CSV](tables/wetland_counterfactual/comparison.csv) | Primary and sensitivity results |
| [Mechanism table](tables/wetland_counterfactual/mechanisms.csv) | Flow volumes, timing, ET and remaining water |
| [Wetland figure](plots/wetland_only/01_wetland_design_mechanisms.png) and [comparison figure](plots/wetland_counterfactual/01_simulated_river_stage.png) | Simulated trajectories and mechanisms |

Run commands from the project root in the existing environment with the verified official solver. The detailed reports provide the full numerical execution instructions. Existing completed runs are preserved by their execution guards. To regenerate only the final reports, figures and tables from saved local results:

```sh
MPLCONFIGDIR=B/work/mpl-cache python3 -B -m B.wetland_only_study --report
MPLCONFIGDIR=B/work/mpl-cache python3 -B -m B.wetland_counterfactual --report
python3 -B -m B.wetland_analysis_checks --report
```

Numerical arrays, source locks and callback histories are local under `B/results/`; source observations and processing files are local under `B/work/`. Both stages share the official solver checkout at `A/work/multivarious/`. Saved configuration and provenance JSON files retain their original path labels and recorded hashes. [project_paths.py](project_paths.py) resolves those labels to the present A/B locations and checks recorded import/path-only source changes against the local relocation manifest. Numerical outputs and scientific values are unchanged.

The [A environment guide](../A/README.md#install-and-run) documents the verified dependency revision and installation status. The detailed final reports above are the entry points for the accepted analysis; this guide does not list individual edits or intermediate study attempts.

## Public files and local reproduction prerequisites

The GitHub repository contains the source code, configurations, English explanations, PNG figures and readable/CSV result tables. Raw observations and processed weather remain local in `B/work/`; completed trajectories and source locks remain local in `B/results/`. Public source links, retrieval hashes and per-value derivations document origin without pretending those excluded artifacts are in a fresh clone.

Running the present B workflow elsewhere additionally requires the exact prepared forcing and manifest, frozen reference schedule NPZ, recorded source/relocation checks and the pinned official solver checkout. Report-only commands require the corresponding completed result files. The current repository has not been tested as a self-contained fresh-clone B reproduction bundle. The A instructions separately document the installable synthetic example.
