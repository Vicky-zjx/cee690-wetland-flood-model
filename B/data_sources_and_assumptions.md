# B data sources and assumptions: fixed-gate wetland design

The accepted B analysis optimizes the two wetland designs under one fixed gate schedule. Water passes from each upstream pool to its associated wetland **only through the modeled gate**. The project is a hypothetical engineered 0.66 km² installation driven by regional event weather, not a calibrated model of the Underwood Creek gauge catchment or a certified local design. The [B overview](README.md) presents the final analysis and results.

The active scientific model is [no_spillway_model.py](no_spillway_model.py), and the optimization driver is [wetland_only_study.py](wetland_only_study.py). The fixed physical inputs are recorded in [no_spillway_config.json](no_spillway_config.json); [wetland_only_config.json](wetland_only_config.json) freezes that configuration and the saved reference gate schedule by checksum. The study loader verifies both hashes before execution. Saved provenance paths retain their original labels and are resolved by [project_paths.py](project_paths.py). These source records remain immutable.

Current optimization results, tables and figures are under `B/results/wetland_only/`, `B/tables/wetland_only/` and `B/plots/wetland_only/`. The [optimization report](docs/12_wetland_only_optimization.md) records actual trajectories, feasibility, search termination and timestep checks. The separate fixed no-wetland comparison is documented below and in the [comparison report](docs/13_no_wetland_counterfactual.md).

## How to read the evidence

This page is the current value-by-value source record. **Origin in code and empirical evidence are different:** a code link identifies the implemented value; it does not establish that the value was measured or calibrated. All entries retain their actual evidence category:

| Category | Meaning in this study |
|---|---|
| Observed | Agency station observations, with station, variable, date range and quality treatment identified. NOAA airport observations require a separate spatial transfer to the hypothetical installation. |
| Derived from observations | A stated calculation on those records, such as unit conversion, rainfall totals or time-weighted RH/pressure means. |
| External physical support / engineering analogue | A specific constant, formula or limited engineering concept supported by an identified primary source; not a measurement of the proposed facility. NASA radiation is separately labeled a gridded estimate. |
| Explicit design or boundary assumption | A selected hypothetical geometry, coefficient, initial condition or operating condition. Approval allows the exploratory case; it does not supply empirical support. |
| Numerical setting | Grid size, optimization start, tolerance, scaling, budget or display choice. |
| Unsupported / unresolved | A local measurement, calibration or site-specific justification has not been obtained. An approved assumed number may still lack local empirical support; unavailable observations are not invented. |
| Computed model result | A value from a specific saved simulation or optimizer record. These are not source observations and do not become field evidence through numerical verification. |

[Weather provenance](docs/05_weather_observations.md) provides exact agency queries, fields, dates, reconstruction assumptions and hashes. [USGS provenance](docs/04_usgs_observations.md) separately records the contextual gauge data. The [active parameter inventory](#active-parameter-provenance) below covers every inherited constructor attribute, active B override, fixed gate knot, initial state, boundary and numerical setting, with inactive fields marked explicitly.

The [source catalog](source_catalog.json), [original evidence register](tables/evidence_register.csv) and [first decision table](tables/decision_values.csv) retain their historical preparation-stage text and mixed old/new values. Their unapproved/not-simulated flags and old gate/spillway entries are not the current state of the accepted study. Frozen configuration and source records have not been rewritten to erase that history. The implemented A source is now [A/wetland_flood_model.py](../A/wetland_flood_model.py), replacing the catalog's historical relative path.

## Accepted optimization problem

The continuous vector is `v = [rho_1, rho_2]`, with each cover fraction in the full closed interval [0,1]. The category pairs A/A, A/B, B/A and B/B are enumerated separately; category labels are not continuous variables. All 24 gate-opening knots and their times, geometry, initial storages, weather, boundary conditions and physical limits are fixed. The schedule comes from the saved hydraulically admissible reference identified in the frozen study specification, rather than the original gate values retained in the physical configuration's provenance.

The physical objective is `min_v max_t H_r(t)`. It is computed from the completed river-stage history outside the ODE; it is not an accumulated cost, discharge or volume objective. The existing OPT callback conditions it as `(H_peak - 97 m)/(1 m)` without changing the minimizer. Physical peaks and numerical penalty values are reported separately.

Exactly 18 residuals are supplied to OPT: eight nonnegative-storage conditions, one soil-capacity limit, four upstream water-level limits, four model-domain conditions and one integration-status residual. Each physical residual uses a fixed positive scale and is satisfied at `g_j <= 0`, subject to the declared numerical tolerance of 1e-6. Cover bounds are enforced separately as the native two-variable box. G5 is an RHS/Runge–Kutta hydraulic-applicability guard and a separate diagnostic, **not an OPT residual**. A failed integration has no physical river peak; its finite numerical failure response rejects the evaluation without completing an invalid trajectory.

## Topology and model scope

The eight states remain `S_u, S_m, S_w1, S_w2, S_d1, S_d2, S_g, S_r`. The current surface-water paths are:

```text
                              gate 1 -> wetland 1 -> normal drainage / wetland overflow --+
upland surface runoff -> pool 1                                                          |
                      -> pool 2 -> gate 2 -> wetland 2 -> normal drainage / wetland overflow -> river -> outlet
                      -> direct upland-to-river bypass ----------------------------------+
```

Soil infiltration, soil-to-groundwater transfer, wetland-to-groundwater seepage, groundwater-to-river flow, external groundwater exchange, direct precipitation and evapotranspiration remain in the eight-state system.

**There is no modeled dam spillway, emergency pool overflow, automatic gate closure or water-level clipping in the accepted model.** The former pool-to-wetland dam-spillway fluxes `Q_dsp1` and `Q_dsp2` are zero and do not appear as active transfers in the balances. Their legacy output columns may remain for array compatibility and verification. Historical dam-spillway coefficients and onset elevations are not an active discharge law. Their retention in archived parameter descriptions must not be read as an emergency-release mechanism.

The two pool and wetland balances therefore include, for branch `i`:

```text
dS_di/dt = P A_di + Q_udi - Q_Gi - E_di

dS_wi/dt = P A_wi + Q_Gi - Q_wri - Q_wspi - Q_wgi - E_wi
```

**Wetland overflow `Q_wspi` is retained.** It remains an internal transfer from its wetland to the receiving river; it is distinct from the removed dam spillway. Normal wetland drainage and wetland overflow must be reported separately when interpreting retention and downstream timing. The direct upland-to-river bypass is also retained.

The unchanged water-level constraints remain mandatory. A numerically valid trajectory above a design limit is an infeasible design, not a successful emergency-spill simulation. Removing the dam spillway does not justify relaxing limits or continuing unsupported gate hydraulics. The model cannot describe real embankment overtopping, breach or other failure once the assumed operating design is exceeded.

## Geography, fixed assumptions and provenance

The existing [source catalog](source_catalog.json), [data-source audit](docs/03_data_sources.md), [USGS observations](docs/04_usgs_observations.md) and [weather observations](docs/05_weather_observations.md) retain the detailed source locators, retrieved records, quality checks and known gaps. The fixed-gate study uses those existing records and assumptions; it introduces no calibration.

| Quantity or input | Value or treatment retained | Evidence category and applicability |
|---|---|---|
| Hypothetical surface footprint | 0.660 km²: upland 0.600; pools 0.010/0.008; wetlands 0.012/0.010; river 0.020 km² | Approved engineering area assumptions; not a delineated parcel or the 46.879 km² station catchment |
| Event window | 2022-09-09 00:00 through 2022-09-15 00:00 UTC | Approved six-day study window |
| USGS 04087088 discharge and stage observations | Contextual event evidence only | Local observations; not inserted as external inflow and not a simulated-outlet calibration target |
| KMKE rainfall | 145.3 mm quantified prepared event depth | Regional station observations spatially transferred to the hypothetical footprint |
| Five event rainfall trace intervals | Zero quantified contribution in prepared forcing | Explicit lower-bound assumption; trace precipitation is not established to be physically zero |
| KMWC temperature and wind | Prepared histories from accepted observations | Regional station observations with documented interpolation and spatial transfer; no unsupported wind-height correction |
| Relative humidity | `r_H = 0.7084916611553665` | Event mean derived from accepted KMWC temperature/dewpoint using the existing model relation; fixed for this study |
| Station pressure | `p_0 = 98692.51171875 Pa` | Event mean derived from accepted KMWC station-pressure observations; not a sea-level pressure substitution |
| NASA POWER shortwave radiation | 21166.81 Wh/m² prepared event integral | Gridded estimate spatially transferred to the installation; accepted hourly timestamp convention and accumulation reconstruction retained |
| External surface inflow | `Q_ext = 0 m³/s` | Approved boundary assumption: no outside surface contributing area is modeled |
| Downstream and external groundwater heads | `h_d = h_g = 97.15 m` | Fixed assumptions in the common model datum; not observed Underwood elevations |
| Outlet blockage | `B = 0` | Assumed maintained outlet, not an observed blockage record |
| Eight initial storages | `[600,90000,3000,2500,3000,2400,409600,7000] m³` | Approved assumed conditions; initial soil storage is not an observed reconstruction of September 2022 |
| Other effective geometry, soil and groundwater coefficients | Existing approved configuration and A.2 parameter values, except removal of active dam-spillway discharge | Fixed uncalibrated assumptions; no claim of local measurements or parameter fitting |

Prepared forcing remains `B/work/event_data/proposed_forcing_120s.csv`, with provenance and checksum in its companion manifest. The source rainfall accumulation intervals and solar-energy integrals are preserved by the existing reconstruction. Finer integration steps interpolate the same prepared 120-second knots; they do not reconstruct new weather observations or change accumulated input energy/water. UTC alignment, rainfall and solar integrals should be reported from actual new runs.

## Fixed design limits and gate applicability

The fixed maximum allowable water levels remain:

```text
pool 1:    101.20 m      pool 2:    101.00 m
wetland 1:  98.90 m      wetland 2:  98.80 m
```

Wetlands are **proposed embanked/contained units**. Settled tops of 101.55/101.35 m for the pools and 99.25/99.15 m for the wetland units are separate engineering assumptions. The 0.35 m reserved freeboard is distinct from those absolute elevations. The previously reviewed [Wisconsin NRCS Pond 378 document](https://www.nrcs.usda.gov/sites/default/files/2024-12/450-11%20FOTG%20Notice%20WI-108%20Combined.pdf) informs a freeboard concept only; it does not supply these elevations, establish a site-specific hazard classification or certify the installation. The former relationship between a dam-spillway onset and a crest is historical context, not an active outlet design in the accepted model.

The unchanged assumed gate widths are 1.0/0.8 m; maximum openings are 0.30/0.30 m; sills are 100.00/99.80 m. The free-flow coefficient is 0.60, the retained submerged-flow coefficient is 0.80, and approach velocity head is neglected. These are documented engineering assumptions, not local calibration. The retained [USACE sluice-gate relationships](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.6/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways/sluice-gate) and [transition convention](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/7.0/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways) support the implemented hydraulic form and conventions, not surveyed geometry or locally calibrated coefficients.

The G5 hydraulic-applicability condition is:

```text
1.25 * opening_height - upstream_depth_above_sill <= 0
```

It applies to open, wet, forward-flow gates and is checked during integration, including the existing Runge–Kutta stage evaluations. The prescribed gates are input constants; these two G5 checks are excluded from the 18-row OPT residual vector. A violating candidate stops at the first detected invalid hydraulic evaluation and is recorded as infeasible; no remainder of that trajectory is manufactured. The gate setting, geometry and restriction are not silently modified to complete a run.

Under the retained wetland-tailwater assumption, every safety-feasible wetland level is below its associated gate sill: `98.90 < 100.00` and `98.80 < 99.80 m`. Thus wetland-induced gate submergence is excluded for safety-feasible trajectories. The remaining gate-applicability issue is insufficient upstream head for the prescribed opening.

## Synthetic vegetation parameters and cover bounds

The A/B labels are synthetic hydraulic/evapotranspiration response classes. They are not named species, field-estimated vegetation properties or recommendations specific to Milwaukee. The following values come from the existing A `class c` and the frozen B configuration, rather than local measurements or an empirical fit:

| Endpoint | Resistance, branch 1 / 2 (s/m²) | ET multiplier, branch 1 / 2 (dimensionless) | Origin and rationale |
|---|---|---|---|
| Bare contained wetland | 0.8 / 1.0 | 1.00 / 1.00 | `R_bare` and `K_bare` in [A class c](../A/wetland_flood_model.py); inherited illustrative baseline |
| A (`type_a`) | 2.8 / 3.0 | 1.10 / 1.10 | `categories.endpoints.type_a` in [no_spillway_config.json](no_spillway_config.json), matching A; assumed lower resistance and smaller ET increase |
| B (`type_b`) | 4.0 / 4.5 | 1.25 / 1.25 | `categories.endpoints.type_b` in the same configuration, matching A; assumed larger resistance and ET increase |

The paired endpoints provide a simple, internally consistent contrast for a course design exercise. Their numerical separation is an explicit synthetic design assumption; no literature-derived typical range or numerical uncertainty distribution is claimed. Positive resistances prevent division by zero in the drainage law, and finite nonnegative ET multipliers retain the model's demand convention. These properties support evaluating the assumed design box, not a claim that either endpoint is realistic for a particular plant community.

The inherited `class c.design()` relation is reused without a new response law:

```text
R_i(rho_i, tau_i) = (1 - rho_i) R_bare,i + rho_i R_veg[tau_i],i
K_i(rho_i, tau_i) = (1 - rho_i) K_bare,i + rho_i K_veg[tau_i],i
```

Here `rho_i` is dimensionless, so interpolation preserves resistance units of s/m² and dimensionless `K_i`. In the drainage rating, a head difference in metres divided by resistance in s/m² gives m³/s. Multiplying the existing evaporation demand in m/s by wetland area in m² and `K_i` gives an ET rate in m³/s before the model's dimensionless availability factor.

For class A, the resulting resistances are `R_1 = 0.8 + 2.0 rho_1` and `R_2 = 1.0 + 2.0 rho_2` s/m²; for B they are `R_1 = 0.8 + 3.2 rho_1` and `R_2 = 1.0 + 3.5 rho_2` s/m². The corresponding ET multipliers are `1 + 0.10 rho_i` for A and `1 + 0.25 rho_i` for B. The reference A/A covers of 0.65/0.60 therefore give `R = [2.10, 2.20] s/m²` and `K = [1.065, 1.060]`.

The full bounds `0 <= rho_i <= 1` follow the definition of an area fraction: zero is the bare endpoint and one is the full synthetic-class endpoint. No observed planting constraint or independently justified narrower interval is available, so the inherited full interval is retained. Both coefficients stay within their stated endpoints throughout it. A zero cover keeps the wetland compartment, geometry, storage, drainage, seepage and overflow; it does not create the no-wetland counterfactual. Physical feasibility still requires simulating and checking each design. Bounds alone do not establish numerical stability or acceptable water levels.

The existing [FAO56 relationships](https://www.fao.org/4/x0490e/x0490e07.htm) support selected vapour-pressure relationships and physical constants, not these synthetic `K` endpoints or the full reduced evaporation model. FAO crop coefficients defined relative to FAO reference evapotranspiration cannot be inserted as this model's `K` without a compatible derivation. Similarly, Manning roughness is not numerically interchangeable with the lumped resistance `R`. Neither source concept supplies local evidence for the endpoint numbers above.

## Attribution under fixed gates

Each category-specific search may change both covers, while the gate schedule, forcing, initial states, fixed parameters and limits stay identical. The reference is A/A at covers 0.65/0.60. Reference-to-design differences therefore quantify wetland-design effects under those fixed operations; a comparison of optimized category candidates also includes their different covers and does not isolate category at equal cover.

Class and cover changes jointly alter normal drainage resistance and ET, with consequent changes in wetland storage, overflow and groundwater exchange. They identify effects of the synthetic response package in this model, not an isolated resistance effect or a measured ecosystem benefit. Feasibility is assessed for every simulated design. Results from earlier joint gate-and-vegetation searches are separate historical studies and do not supply the current optimization benefit.

## Mechanism and conservation reporting

Use complete valid trajectories to report each branch's gate-release volume, release-weighted centroid time, peak-release time, wetland peak storage/level and time, normal-drainage volume, wetland-overflow volume, wetland ET and seepage, together with the physical river peak and its time. Release centroids are undefined for zero integrated release and must not be invented. Integrated volumes and centroids use the declared reporting-grid quadrature; grid sensitivity remains a numerical check.

The fixed-gate study checks gate openings and actual gate releases across designs so its attribution can be assessed. Changes in wetland storage, drainage, overflow and ET explain downstream timing and retention within the modeled system. Peak time and release centroid answer different questions and should be labeled separately. Delaying a peak is not the same as reducing its maximum value.

For total system storage, internal transfers cancel. The whole-system balance is:

```text
change in total storage = integral(P_total + Q_ext + Q_gx - Q_out - E_tot) dt
```

`Q_gx` is signed external groundwater exchange. Gate releases, normal wetland drainage, retained wetland overflow, upland bypass, infiltration and internal seepage must not be added again as external inputs or outputs. Report both numerical balance residuals and changes in final storage; a finite event can redistribute water beyond its end without proving a reduction over an indefinitely extended horizon.

## Historical documents and reproducibility

The [first exploratory configuration](docs/01_study_configuration.md), [first implementation](docs/02_model_and_optimization.md), [first results](docs/09_exploratory_results.md) and [initial-guess check](docs/10_multistart_check.md) describe the preserved joint-search model with active dam-to-wetland spillways. The later [no-spillway joint study](docs/11_no_spillway_revision.md), executed by [no_spillway_study.py](no_spillway_study.py), retained joint gate/cover decisions. Their numerical results and method descriptions are historical records, not the final two-cover optimization or its reference comparison.

The accepted fixed gate schedule is preserved in the source trajectory identified by `wetland_only_config.json`; the wetland-only study plan records its 24 values. The original B configuration, official solver checkout, downloaded observations and completed numerical outputs retain their values and hashes. The current reproduction instructions and completed fixed-design resolution checks are in the [wetland-only report](docs/12_wetland_only_optimization.md) and [B overview](README.md).

## Fixed six-state no-wetland counterfactual

The separately authorized [dynamic counterfactual](docs/13_no_wetland_counterfactual.md) adds one reduced model in [wetland_counterfactual.py](wetland_counterfactual.py). This is a comparison assumption, not an externally sourced local engineering design. It does not revise the eight-state scientific model above or its archived results, and it does not run OPT.

In this case the six states are `S_u, S_m, S_d1, S_d2, S_g, S_r`; each gate discharges directly into the river. Wetland storage, normal drainage, overflow, seepage and wetland ET are absent. No independent dam spillway is introduced. Rain over the former 12,000/10,000 m² wetland footprints goes immediately to the river as `P*(A_w1+A_w2)`. This adopted instantaneous bypass preserves the 660,000 m² rainfall footprint, but does not change the 20,000 m² river storage/evaporation area or the 640,000 m² aquifer footprint. It is not a calibrated replacement land-cover scenario or an isolation of storage alone.

The same saved gate-opening histories drive all cases; gate discharge is recalculated using internal river stage as bypass tailwater and internal wetland stage in the eight-state cases. The existing gate law, G5 validation at RHS/RK stages, shared limits and initial states remain unchanged. Wetland-specific constraints become N/A. The independent inherited groundwater ceiling `H_g_valid_max=97.9 m` remains despite removal of the duplicate wetland-bed separation check; no shared applicability limit is loosened.

Primary wetland cases retain 3,000/2,500 m³ of initial wetland water. The six-state case has no such water, so total initial water is 512,600 m³ versus 518,100 m³. None is transferred into the river. A separately labeled sensitivity sets only the two wetland initial storages to zero, equalizing initial total water without replacing approved initial conditions. All comparisons reuse the existing UTC regional weather reconstruction and 120-second source knots at 120/60/30-second integration steps. No new observations, citations, coefficients or safety thresholds are introduced. Results are conditional simulations; see the report for actual runs, conservation, feasibility, initial-water sensitivity and finite-window limitations.

## Active parameter provenance

This inventory describes the **accepted fixed-gate wetland-only B model**, not the historical joint-search input table. Arrays are ordered branch 1 / branch 2 unless a state order is stated. All elevations use the assumed model datum. “Derived” is separated into derivation from observations and derivation from assumptions; converting an assumed depth into a storage does not make it observed.

Implementation precedence is: inherited A constants → `no_spillway_model.base_parameters()` plus `no_spillway_config.json` overrides → category endpoints and frozen initial state → the **frozen 24 gate values** loaded by `wetland_only_study.load_study()` → current wetland-only numerical settings. A source-code link establishes where an assumption is implemented; it is not empirical evidence for that value. An approved assumption can remain unsupported by site measurements.

### Geometry

| Quantity / code field | Active value and units | Evidence status | Exact origin, derivation and limits |
|---|---|---|---|
| `A_u` | 600000 m² | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Hypothetical upland footprint. |
| `A_d` | [10000,8000] m² | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Pool surface areas, branch 1 / 2; constant-area reservoirs. |
| `A_w` | [12000,10000] m² | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Contained wetland surface areas, branch 1 / 2; not mapped natural wetlands. |
| `A_r` | 20000 m² | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). River storage and evaporation footprint; not a surveyed cross-section. |
| `A_g` | 640000 m² | Derived from approved engineering assumptions; not derived from local measurements | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). A_u + sum(A_d) + sum(A_w); effective aquifer area. Not an additional rainfall footprint. |
| `A_rain` | 660000 m² | Derived from approved engineering assumptions; not derived from local measurements | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`) (`A_rain` property). A_u + sum(A_d) + sum(A_w) + A_r = 0.66 km²; soil/aquifer do not receive duplicate rainfall. |
| `S_y` | 0.2 1 | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Assumed effective groundwater specific yield; no well or soil-data fit. |
| `z_d` | [100.0,99.8] m | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Pool beds in an assumed common datum. |
| `z_w` | [98.0,97.9] m | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Wetland beds in the same assumed datum. |
| `z_r`, `z_g` | [97.0,94.0] m | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). River and groundwater zero-storage reference elevations; not NAVD88. |

### Soil and runoff

| Quantity / code field | Active value and units | Evidence status | Exact origin, derivation and limits |
|---|---|---|---|
| `C_m` | 120000 m³ | Derived from approved engineering assumptions; not derived from local measurements | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). A_u × assumed 0.20 m water capacity; not an SSURGO estimate. |
| `S_fc` | 60000 m³ | Derived from approved engineering assumptions; not derived from local measurements | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). A_u × assumed 0.10 m field-capacity storage; not an observed initial condition. |
| `K_u` | 1e-06 m/s | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Uncalibrated lumped infiltration coefficient; no local measured conductivity. |
| `k_mg` | 2.3148148148148148e-06 s⁻¹ | Derived from approved engineering assumptions; not derived from local measurements | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). 1/(5 days), from an assumed drainage timescale. |
| `k_ud` | [3.472222222222222e-05,2.777777777777778e-05] s⁻¹ | Derived from approved engineering assumptions; not derived from local measurements | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Pool routing from assumed 8 h / 10 h time constants. |
| `k_ur` | 1.388888888888889e-05 s⁻¹ | Derived from approved engineering assumptions; not derived from local measurements | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Direct river bypass routing from an assumed 20 h time constant. |
| `Q_udcap`, `Q_urcap` | {"Q_udcap":[1.2,1.0],"Q_urcap":0.7} m³/s | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Caps on upland-to-pool and direct upland-to-river transfers; not measured capacities. |

### Wetland and outlet hydraulics

| Quantity / code field | Active value and units | Evidence status | Exact origin, derivation and limits |
|---|---|---|---|
| `z_wo` | [98.05,97.95] m | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Normal-drainage sill levels; assumed 0.05 m above each wetland bed. |
| `z_wsp` | [98.6,98.5] m | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Retained wetland-overflow onset levels; assumed 0.60 m above each wetland bed. Not removed dam-spillway elevations. |
| `c_wsp` | [3.0,2.5] m^(3/2)/s | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`); [Revised B model](no_spillway_model.py) (`make_b2_rhs`). Retained Q_wsp coefficient in c_wsp × positive head^(3/2); no measured rating. |
| `Q_wrcap` | [0.65,0.55] m³/s | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Normal-wetland-drainage capacity; not the overflow capacity. |
| `C_wg` | [0.015,0.012] m²/s | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Lumped wetland-to-groundwater conductance; no observed seepage fit. |
| `C_gr`, `C_gx` | {"C_gr":0.08,"C_gx":0.035} m²/s | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Groundwater-to-river and external groundwater conductances; no local calibration. |
| `z_ro` | 97.1 m | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). River-outlet sill, assumed 0.10 m above river zero-storage bed. |
| `R_r`, `Q_rcap` | {"R_r":0.45,"Q_rcap":3.0} R_r: s/m²; Q_rcap: m³/s | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Lumped final-outlet resistance/capacity. A head divided by R gives m³/s. |

### Availability volumes

| Quantity / code field | Active value and units | Evidence status | Exact origin, derivation and limits |
|---|---|---|---|
| `V_u`, `V_m` | {"V_u":600,"V_m":6000} m³ | Derived from approved engineering assumptions; not derived from local measurements | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). A_u × assumed 1 mm / 10 mm water depths. These are constitutive half-availability volumes, not numerical clipping thresholds. |
| `V_w`, `V_d`, `V_r` | {"V_w":[120,100],"V_d":[100,80],"V_r":200} m³ | Derived from approved engineering assumptions; not derived from local measurements | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). 0.01 m × each wetland/pool/river area; constitutive assumptions, not calibrated capacities. |
| `V_g` | 6400 m³ | Derived from approved engineering assumptions; not derived from local measurements | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Assumed 0.05 m head × S_y × A_g; not an observed groundwater volume. |

### Evaporation physics and closure

| Quantity / code field | Active value and units | Evidence status | Exact origin, derivation and limits |
|---|---|---|---|
| `lambda_v` | 2450000.0 J/kg | Externally supported physical approximation | [FAO 56, Chapter 3](https://www.fao.org/4/x0490e/x0490e07.htm) (latent heat; 2.45 MJ/kg); [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Unit conversion × 10^6. Does not validate the complete reduced evaporation closure. |
| `c_p` | 1013 J/(kg K) | Externally supported physical approximation | [FAO 56, Chapter 3](https://www.fao.org/4/x0490e/x0490e07.htm) (psychrometric constant, Eq. 8); [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). 1.013×10^-3 MJ/(kg °C) converted to J/(kg K) for temperature differences. |
| `epsilon` | 0.622 1 | Externally supported physical approximation | [FAO 56, Chapter 3](https://www.fao.org/4/x0490e/x0490e07.htm) (Eq. 8); [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Water-vapour/dry-air molecular-weight ratio; not a tuned basin parameter. |
| `rho_w` | 1000 kg/m³ | Externally supported physical approximation | [FAO 56, Chapter 3](https://www.fao.org/4/x0490e/x0490e07.htm) (Box 8); [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Constant water-density approximation. |
| `R_v` | 461.5 J/(kg K) | Externally supported physical approximation | [JRA-25 physical-constant table, NCAR/UCAR](https://osdata.gdex.ucar.edu/web/datasets/d625000/docs/jra-25_model_constants.html); [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Water-vapour gas constant only; other constants in that source are not adopted wholesale. |
| `saturation_pressure_constants` | [610.8,17.27,237.3] Pa, 1, °C | Externally supported relation | [FAO 56, Chapter 3](https://www.fao.org/4/x0490e/x0490e07.htm) (Eq. 11; Eq. 13 derivative); [evaporative_demand](../A/wetland_flood_model.py). e_s = 610.8 exp(17.27 T/(T+237.3)); 0.6108 kPa converted to Pa. T_K=T+273.15 is the unit conversion. |
| `p_0` | 98692.51171875 Pa | Derived from regional station observations; spatial transfer | [Approved B physical configuration](no_spillway_config.json) (`parameter_overrides.p_0`); [weather derivation](docs/05_weather_observations.md); [NOAA ISD](https://www.ncei.noaa.gov/products/land-based-station/integrated-surface-database). KMWC MA1 event-time mean of station pressure. Overrides A default 101325 Pa. Raw/processed records and integration are retained locally in B/work/event_data/weather/. |
| `r_H` | 0.7084916611553665 1 | Derived from regional station observations; spatial transfer | [Approved B physical configuration](no_spillway_config.json) (`parameter_overrides.r_H`); [weather derivation](docs/05_weather_observations.md); [NOAA ISD](https://www.ncei.noaa.gov/products/land-based-station/integrated-surface-database). KMWC temperature/dewpoint → relative humidity using retained saturation relation → time mean. Overrides A default 0.65; not observed installation humidity. |
| `eta`, `k_0`, `k_v` | {"eta":0.6,"k_0":0.001,"k_v":0.002} eta, k_v: 1; k_0: m/s | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Radiation multiplier and wind-transfer closure k_a=k_0+k_v v; external physical references do not supply these numbers. |
| `K_ETm`, `K_ETd`, `K_ETr` | {"K_ETm":0.85,"K_ETd":[1.0,1.0],"K_ETr":1.0} 1 | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Soil, pool and river multipliers for this reduced demand; not FAO crop coefficients or local ET observations. |

### Wetland design

| Quantity / code field | Active value and units | Evidence status | Exact origin, derivation and limits |
|---|---|---|---|
| `R_bare`, `K_bare` | {"R_bare":[0.8,1.0],"K_bare":[1.0,1.0]} R: s/m²; K: 1 | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Bare contained-wetland endpoint; zero cover still retains the wetland compartment. |
| `R_veg`, `K_veg` | A: R=[2.8,3.0], K=[1.10,1.10]; B: R=[4.0,4.5], K=[1.25,1.25] R: s/m²; K: 1 | Approved engineering assumption; no local numeric measurement or calibration | [Approved B physical configuration](no_spillway_config.json) (`categories.endpoints`), applied in [Revised B model](no_spillway_model.py) (`simulate`). A: R=[2.8,3.0], K=[1.10,1.10]; B: R=[4.0,4.5], K=[1.25,1.25]. Synthetic response classes, not species. Active endpoint values are explicitly loaded from B config even though numerically identical to A defaults. |
| `tau`, `rho` | {"pairs":["A/A","A/B","B/A","B/B"],"rho_bounds":[0,1],"reference_pair":"A/A","reference_rho":[0.65,0.6]} category; cover fraction 1 | Approved design definition / assumed reference | [Frozen wetland-only specification](wetland_only_config.json); [Wetland-only driver](wetland_only_study.py) (`assemble_design`, `run_study`); [A interpolation](../A/wetland_flood_model.py) (`design`). R=(1-rho)R_bare+rho R_veg and K=(1-rho)K_bare+rho K_veg. No empirically justified narrower cover range is asserted. Optimized covers are model outputs, not source observations. |

### Gates

| Quantity / code field | Active value and units | Evidence status | Exact origin, derivation and limits |
|---|---|---|---|
| `gravity` | 9.80665 m/s² | Externally supported conventional constant | [NIST SP 1038](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication1038.pdf) (section 4.3.3, footnote 6; section 5.1.8); [Approved B physical configuration](no_spillway_config.json) (`gates.gravity`). Standard gravity, not measured Milwaukee acceleration. |
| `gate_width`, `gate_max_opening`, `gate_sill_z` | {"width":[1.0,0.8],"maximum_opening":[0.3,0.3],"sill":[100.0,99.8]} m | Approved engineering assumption; no local numeric measurement or calibration | [Approved B physical configuration](no_spillway_config.json) (`gates`); [gate law](gate_hydraulics.py). Hypothetical gate geometry; sill equals pool bed in this case. Physical opening is gate_max_opening × normalized fixed gate control. |
| `Cd_gate_free`, `Cd_gate_submerged` | {"free":[0.6,0.6],"submerged":[0.8,0.8]} 1 | Approved engineering coefficients informed by external hydraulic analogues; not calibrated | [USACE HEC-RAS 6.6 sluice gate](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.6/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways/sluice-gate); [submerged orifice reference](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.6/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways/radial-gates); [Approved B physical configuration](no_spillway_config.json). Sources inform the law/coefficient convention; these specific hypothetical coefficients are not local measurements. Submerged coefficient is implemented but cannot affect safety-feasible trajectories here because wetland levels stay below the sills. |
| `neglect_approach_velocity_head`, `one_way_gate` | [true,true] Boolean | Approved engineering assumption; no local numeric measurement or calibration | [Approved B physical configuration](no_spillway_config.json); [sluice_gate_flow](gate_hydraulics.py). Upstream water surface is the energy elevation. Closed/dry/adverse-head cases return zero; no reverse gate flow. |
| `G5_head_opening_ratio`, `submergence_lower`, `submergence_upper` | [1.25,0.67,0.8] 1 | Hydraulic-method applicability and external regime conventions | [USACE HEC-RAS 6.0 manual](https://www.hec.usace.army.mil/software/hec-ras/documentation/HEC-RAS_6.0_Reference_Manual.pdf) pp. 8-7–8-12; [HEC-RAS 7.0 transition convention](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/7.0/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways); [implementation](gate_hydraulics.py). Require H/a>=1.25 for open/wet/forward gates. Submergence cutoffs 0.67/0.80 and documented weighted transition are retained; in the safety-feasible case submergence is excluded. G5 is checked at RHS/RK evaluations and is NOT an active wetland-only OPT residual. |
| `gate1_hours`, `gate2_hours` | {"gate1_hours":[0,17,18,26,28,50,65,66,76,80,110,144],"gate2_hours":[0,18,20,27,30,52,66,68,78,82,112,144]} h | Prescribed schedule parameterization / numerical knot times | [Approved B physical configuration](no_spillway_config.json) (`controls`); [Wetland-only driver](wetland_only_study.py) (`load_study`). Times retained from A, linearly interpolated. The opening values in no_spillway_config.controls are provenance, not the current frozen openings. |
| `fixed_gate1_values`, `fixed_gate2_values` | {"g1":[0.005,0.005,0.005,0.005,0.005,0.005,0.336,0.384,0.696,0.36,0.005,0.005],"g2":[0.005,0.005,0.005,0.005,0.005,0.005,0.36,0.42,0.78,0.324,0.005,0.005]} 1 | Prescribed hydraulically admissible engineered schedule; not observations | [Frozen wetland-only specification](wetland_only_config.json) (`fixed_gate_source`, SHA256); [Wetland-only driver](wetland_only_study.py) (`load_study`); [reference-selection report](docs/11_no_spillway_revision.md). Read from local `B/results/no_spillway_revision/reference.npz` and retained verbatim in the public [optimization table](tables/wetland_only/optimization_summary.csv), column `fixed_gate_knots_json`, as well as the local study plan. Selected label `retention_storm_seed_x1.2`: first fully feasible candidate in a declared ordered finite pool, selected after 8 actual attempts. Dry-period knots are 0.005 and storm knots were a prescribed retention/release proposal informed by event timing; this is a replacement feasibility reference, not observed operations or a new optimized gate history. |

### Safety and model domain

| Quantity / code field | Active value and units | Evidence status | Exact origin, derivation and limits |
|---|---|---|---|
| `dam_m`, `wetland_m` | {"dam":[101.2,101.0],"wetland":[98.9,98.8]} m | Approved fixed hypothetical design limits | [Approved B physical configuration](no_spillway_config.json) (`safety`). These four maximum levels remain hard feasibility limits; high pool storage is not clipped or relieved by a dam spillway. |
| `settled_crest_m`, `physical_freeboard_m` | {"dam_crests":[101.55,101.35],"wetland_crests":[99.25,99.15],"freeboard":0.35} m | Explicit engineering design assumptions, not local certification | [Approved B physical configuration](no_spillway_config.json) (`safety`); [Wisconsin NRCS Pond 378 concept](https://www.nrcs.usda.gov/sites/default/files/2024-12/450-11%20FOTG%20Notice%20WI-108%20Combined.pdf). H_safe = assumed settled crest − assumed 0.35 m reserve. NRCS informs the freeboard concept only; it does not provide these elevations or certify this installation. Crest metadata is a derivation of limits, not an active overflow law. |
| `H_r_valid_max`, `H_g_valid_max`, `soil_bottom` | {"river_ceiling":99.0,"groundwater_ceiling":97.9,"soil_bottom":100.0} m | Approved engineering assumption; no local numeric measurement or calibration | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`); [Revised B model](no_spillway_model.py) (`make_b2_rhs`); [constraint audit](optimize_wetland.py). Separate model-domain constraints: H_r<99; H_g<97.9, H_g<100, H_g<min(z_w)=97.9. Duplicate groundwater bounds remain explicitly reported; these are not observed safety levels. |
| `temperature_domain`, `storage_domain` | {"temperature_C":">0","states_m3":">=0","soil_storage_m3":"<=120000"} °C; m³ | Model scope / physically admissible storage domain | [Revised B model](no_spillway_model.py) (`validate`, `make_b2_rhs`). Liquid-water case only; no snow/ice model. Invalid RK storage or model-domain excursions are rejected, not clipped. |

### Initial conditions

| Quantity / code field | Active value and units | Evidence status | Exact origin, derivation and limits |
|---|---|---|---|
| `x0_m3` | [600.0,90000.0,3000.0,2500.0,3000.0,2400.0,409600.00000000035,7000.0] m³ | Approved assumed physical initial state | [Approved B physical configuration](no_spillway_config.json) (`initial_conditions.x0_m3`); [Revised B model](no_spillway_model.py) (`simulate`). Order [S_u,S_m,S_w1,S_w2,S_d1,S_d2,S_g,S_r]. The frozen vector overrides constructor initial-condition defaults; not reconstructed event observations. |
| `u_depth_init`, `soil_fraction_init`, `w_depth_init`, `d_depth_init`, `H_g_init`, `r_depth_init` | {"upland_depth":0.001,"soil_fraction":0.75,"wetland_depths":[0.25,0.25],"pool_depths":[0.3,0.3],"groundwater_head":97.2,"river_depth":0.35} m; soil fraction 1 | Derived from approved engineering assumptions; not derived from local measurements | [Approved B physical configuration](no_spillway_config.json); [Revised B model](no_spillway_model.py) (`simulate`); [A state-to-volume relations](../A/wetland_flood_model.py). 600=A_u×0.001; 90000=C_m×0.75; [3000,2500]=A_w×0.25; [3000,2400]=A_d×0.30; 409600=(97.20−94)×0.20×640000; 7000=A_r×0.35. Runtime depth attributes are reconstructed from frozen x0; their equality with A values is intentional. |

### Forcing and boundaries

| Quantity / code field | Active value and units | Evidence status | Exact origin, derivation and limits |
|---|---|---|---|
| `event_start_utc`, `event_end_utc`, `duration` | ["2022-09-09T00:00:00Z","2022-09-15T00:00:00Z",518400] UTC; s | Approved event selection / numerical horizon | [Approved B physical configuration](no_spillway_config.json) (`geography`, `numerics.duration_s`). Six days / 144 hours. Date selection is a study choice, not a measured parameter. |
| `P(t)` | 145.3 mm event accumulation | Regional station observations + spatial-transfer and trace-zero assumptions | [NOAA ISD](https://www.ncei.noaa.gov/products/land-based-station/integrated-surface-database); [weather processing record](docs/05_weather_observations.md); [Approved B physical configuration](no_spillway_config.json) (`forcing`). KMKE accepted non-overlapping accumulation intervals. Quantified depth is preserved in prepared knots; five event trace intervals contribute zero as an explicit lower-bound assumption. Not measured areal rainfall at the hypothetical facility. |
| `T_a(t)`, `v(t)` | Prepared 4,321-node event histories °C; m/s | Regional station observations + interpolation / spatial transfer | [NOAA ISD](https://www.ncei.noaa.gov/products/land-based-station/integrated-surface-database); [weather processing record](docs/05_weather_observations.md); [Approved B physical configuration](no_spillway_config.json) (`forcing`). KMWC accepted temperature and wind reports. Interpolation adds no observations; no unsupported wind-height adjustment. |
| `R_s(t)` | 21166.81 Wh/m² event energy | Gridded radiation estimate + spatial transfer / reconstruction assumption | [NASA POWER actual event request](https://power.larc.nasa.gov/api/temporal/hourly/point?parameters=ALLSKY_SFC_SW_DWN&community=RE&longitude=-88.0461111&latitude=43.0547222&start=20220902&end=20220914&format=JSON&time-standard=UTC); [weather record](docs/05_weather_observations.md). ALLSKY_SFC_SW_DWN; returned SYN1DEG product, not local pyranometer measurements. Prepared forcing retains approved hourly-label interpretation and energy-preserving 120 s reconstruction. |
| `Q_ext`, `h_d`, `h_g`, `B` | {"Q_ext":0,"h_d":97.15,"h_g":97.15,"B":0} Q_ext: m³/s; heads: m; B: 1 | Approved engineering assumption; no local numeric measurement or calibration | [Approved B physical configuration](no_spillway_config.json) (`forcing.proposed_fixed_boundaries`); [fixed forcing loader](optimize_wetland.py). Closed surface-drainage area, fixed downstream/aquifer exterior heads, maintained unblocked outlet. No USGS hydrograph or gauge elevation is inserted into the model. |

### Numerics

| Quantity / code field | Active value and units | Evidence status | Exact origin, derivation and limits |
|---|---|---|---|
| `hour`, `day` | [3600,86400] s | Exact unit conversions | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`). Not tuned constants or observed data. |
| `time_step`, `refine_factor`, `comparison_grids` | {"search_dt":120,"default_refine_factor":2,"reported_dt":[120,60,30]} s; refine factor 1 | Numerical setting; not a physical observation | [Approved B physical configuration](no_spillway_config.json) (`numerics`); [Wetland-only driver](wetland_only_study.py) (`compare`). 30 s checks use the same 120 s forcing knots. The boundary check uses 120/60 s and only triggers 30 s if its documented diagnostic fails. |
| `balance_rtol` | 0.0001 1 | Numerical setting; not a physical observation | [A parameter constructor](../A/wetland_flood_model.py) (`class c`); [B inheritance](no_spillway_model.py) (`base_parameters`); [Revised B model](no_spillway_model.py) (`simulate`). Max absolute cumulative balance residual divided by integrated external positive inflow (or initial total storage if no inflow); must be <1e-4. Not physical freeboard or acceptable real water loss. |
| `height_scale_m` | 1.0 m | Numerical setting; not a physical observation | [Approved B physical configuration](no_spillway_config.json) (`numerics`); [Wetland-only driver](wetland_only_study.py). Positive scale for level residuals and objective. |
| `volume_scales_m3` | [600.0,90000.0,3000.0,2500.0,3000.0,2400.0,409600.00000000035,7000.0] m³ | Numerical setting; not a physical observation | [Approved B physical configuration](no_spillway_config.json) (`numerics`); [Wetland-only driver](wetland_only_study.py). Eight fixed scales, equal to the assumed initial storages; do not alter physical storage constraints. |
| `strict_margin_m` | 0.001 m | Numerical setting; not a physical observation | [Approved B physical configuration](no_spillway_config.json) (`numerics`); [Wetland-only driver](wetland_only_study.py). Numerical inward margin for strict model-domain checks only; distinct from 0.35 m physical freeboard. |
| `feasibility_tol` | 1e-06 1 | Numerical setting; not a physical observation | [Approved B physical configuration](no_spillway_config.json) (`numerics`); [Wetland-only driver](wetland_only_study.py). Scaled-residual acceptance tolerance; metre equivalent is 1e-6 × the 1 m level scale. |
| `objective_origin_m` | 97.0 m | Numerical setting; not a physical observation | [Approved B physical configuration](no_spillway_config.json) (`numerics`); [Wetland-only driver](wetland_only_study.py). Numerical conditioning f=(H_peak−97)/1; physical H_peak always reported separately. |
| `failure_objective`, `failure_residual` | 1000.0 1 | Numerical setting; not a physical observation | [Approved B physical configuration](no_spillway_config.json) (`numerics`); [Wetland-only driver](wetland_only_study.py). Invalid integration returns f=1000; only integration-status residual=1000 in the current 18-row callback. No physical river peak is assigned. |
| `peak_grid_tolerance_m` | 0.001 m | Numerical setting; not a physical observation | [Approved B physical configuration](no_spillway_config.json) (`numerics`); [Wetland-only driver](wetland_only_study.py). 0.001 m boundary-comparison diagnostic threshold; not a relaxed safety limit. |
| `OPT_continuous_dimension`, `OPT_residual_count`, `native_bounds`, `cover_starts` | {"dimension":2,"residual_count":18,"bounds":[[0,1],[0,1]],"starts":[[0.65,0.6],[0.3,0.8],[0.8,0.35]]} 1 | Numerical setting; not a physical observation | [Frozen wetland-only specification](wetland_only_config.json); [Wetland-only driver](wetland_only_study.py). Four categorical pairs × three cover starts = 12 bounded runs. Bounds follow cover definition; starts are numerical choices, not observations. Fixed gate controls never enter the continuous optimization vector. |
| `max_evals`, `max_seconds`, `tol_x`, `tol_f`, `optimizer_feasibility_tol`, `penalty`, `penalty_exponent`, `display` | {"max_evals":100,"max_seconds":180,"tol_x":0.0001,"tol_f":1e-06,"feasibility_tol":1e-06,"penalty":100.0,"penalty_exponent":2.0,"display":0} callback count; s; otherwise 1 | Numerical setting; not a physical observation | [Frozen wetland-only specification](wetland_only_config.json) (`settings`); [official OPT adapter](official_opt.py) (`run_opt`). Actual values: 100 callbacks/run, 180 s/run, tol_x=1e-4, tol_f=1e-6, feasibility_tol=1e-6, penalty=100, exponent=2, display=0. These override generic adapter defaults and historical 180-callback settings. |
| `maximum_batch_seconds`, `maximum_study_ode_attempts`, `reserved_comparison_attempts`, `cache_entries` | {"seconds":1500,"ode_attempts":1250,"reserve":20,"cache_entries":8} s; counts | Numerical setting; not a physical observation | [Frozen wetland-only specification](wetland_only_config.json); [Wetland-only driver](wetland_only_study.py) (`run_study`). Wetland-only batch controls. no_spillway_config total_ode_budget=1300 and total_seconds=1500 belong to the earlier joint-search configuration; the current driver uses its own specification. |
| `OPT_hyp7_8_9_18`, `OPT_coordinate_map` | {"hyp7":1,"hyp8":0.1,"hyp9":0,"hyp18":0,"cover_to_solver":"2*rho-1"} 1 | Numerical setting; not a physical observation | [official OPT adapter](official_opt.py) (`run_opt`); [verified solver audit](docs/08_official_opt_audit.md). Adapter assigns first 10 official hyperparameters and inequality-only hyp[18]=0; other entries come from pinned official opt_hyp(). Cover [0,1] is mapped to solver [-1,1]. No joint-study gate odds-cube map is used for covers. |
| `resolution_comparison_rule` | {"relative_effect_threshold":0.1,"absolute_floor_m":1e-10} 1; m | Numerical setting; not a physical observation | [Wetland-only driver](wetland_only_study.py) (`compare`). Combined peak grid changes must be <max(0.1×effect magnitude,1e-10 m), with stable sign and feasibility. This is a resolution diagnostic, not proof of global optimum. |
| `official_revision` | 06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02 Git revision | Verified external solver implementation | [solver setup and audit](docs/08_official_opt_audit.md); [loader](official_opt.py). Official multivarious 1.2.1 ode4u/nms, local checkout excluded from upload. Solver coefficients belong to this version, not independent hydrologic parameters. |

### Historical values that are not active B inputs

| Historical fields / records | Status in the accepted B model |
|---|---|
| `z_dsp`, `c_dsp`, `C_d`; `Q_dsp1`, `Q_dsp2` | Dam-spillway attributes are removed from active parameters. `C_d` is not evaluated by the B model; the legacy output positions contain literal zeros only. No dam-spillway relief or storage cap exists. |
| A linear-gate `z_G=[100.05,99.85]`, `R_G=[1.6,1.8]`, `Q_Gcap=[1.0,0.8]` | Replaced by explicit B gate geometry and sluice law. These inherited object attributes are not used by the revised B RHS. |
| `C_w=A_w*(z_wsp-z_w)` | A derived wetland-overflow onset storage, approximately `[7200,6000] m³`; not a hard storage cap or active B constraint. B retains the actual `z_wsp`, `c_wsp` overflow law. |
| `H_flood=97.9 m` | Historical illustrative reporting threshold. Current `reporting.flood_threshold_m=null`; no locally supported flood-duration threshold is inserted. |
| A `p_0=101325 Pa`, `r_H=0.65` | Overridden by the KMWC-derived event means listed above. |
| A `rain_pulses`, `ext_pulses`, `Q_ext_base`, `weather_period_h`, `temperature_phase_h`, `sunrise_phase_h`, `wind_period_h`, `groundwater_period_h`, `T_mean`, `T_amplitude`, `R_s_peak`, `wind_mean`, `wind_amplitude`, `h_d_base`, `h_d_pulse`, `h_g_base`, `h_g_amplitude`, `B_base` | Synthetic A forcing parameters remain historical and do not drive the prepared B event. B uses the saved event forcing and four explicit fixed boundaries above. |
| A `gate1_values`, `gate2_values`; corresponding `no_spillway_config.controls` values | Historical original openings. Knot times are retained, but the actual fixed B openings are the two arrays listed above. Do not use the old arrays to recreate the accepted study. |
| `tables/evidence_register.csv`, `tables/decision_values.csv` | First exploratory study records; mixed old/new columns and original gate decisions. They are not a complete current-input manifest. Their A values, deleted spillway entries and original optimizer budgets must not be read as active B settings. |

No dataset currently supplies a surveyed installation polygon, absolute facility elevations, local gate calibration, measured soil/groundwater coefficients, observed initial storages, or species-specific resistance/ET endpoints. Those empirical gaps are explicit; the exploratory case proceeds with the approved assumptions above, without attributing them to unprocessed soil, terrain, vegetation or engineering databases. The context USGS station hydrograph is not used to fill those gaps.

## Provenance of reported output values

All reported performance numbers are generated model outputs. An external paper or gauge record is not their numerical source. Units appear in the output column names; `dt_s` or `timestep_s` identifies the reporting grid. Rounded narrative values do not replace the full-precision stored values.

| Reported value | Public numerical source and row/field | Calculation / interpretation |
|---|---|---|
| Optimized covers and selected category | [optimization_summary.csv](tables/wetland_only/optimization_summary.csv), `B/B`, `cover_start_3`, `final_cover1/2` | Best feasible design found by the bounded 120 s search; not a global-optimality proof or species recommendation |
| Physical river peaks, time-to-peak and bypass differences | [comparison.csv](tables/wetland_counterfactual/comparison.csv), primary cases and `timestep_s=30`; `river_peak_m`, `peak_time_h`, `difference_vs_bypass_mm` | Maximum of each saved stage history; peak time from its sampled maximum; differences are subtraction × 1,000 to mm |
| Reference-to-best improvement | Same table, primary `reference` minus `best` at 30 s | `(97.55844251385638 − 97.55343050504051) × 1000 = 5.012008815868 mm`; the arbitrary absolute datum is not used as a percentage denominator |
| Gate volumes/timing, normal drainage, wetland overflow, peak wetland storage | [mechanism_metrics.csv](tables/wetland_only/mechanism_metrics.csv), `fixed_reference`, `wetland_best`, `dt_s=30` | Trapezoidal event integrals; maxima of saved storage/flow; release-weighted centroids. Internal flows are not double-counted as external water |
| Physical/scaled constraint residuals and margins | [constraint_audit.csv](tables/wetland_only/constraint_audit.csv); [dynamic comparison](tables/wetland_counterfactual/comparison.md) | Saved node extrema versus unchanged limits/scales; RHS/RK applicability is separately checked |
| Solver stopping reasons, callbacks and convergence | [optimization_summary.csv](tables/wetland_only/optimization_summary.csv), `termination`, `solver_converged`, `callbacks`, `actual_ode_attempts` | Actual recorded search outcomes: 7 converged / 5 evaluation-budget stops across 12 searches |
| Time-step sensitivity | [fixed_design_timesteps.csv](tables/wetland_only/fixed_design_timesteps.csv) and [comparison.csv](tables/wetland_counterfactual/comparison.csv) | Frozen designs and forcing checked at 120/60/30 s; no refinement optimization |
| Cover-bound peaks, feasibility and balance residuals | [boundary_checks.md](tables/analysis_checks/boundary_checks.md) | Four corners × four pairs × two time steps; 32 actual dynamic evaluations, all valid and feasible; maximum sampled peak change 0.010340 mm |
| Time-history lines in the three boundary figures | [boundary analysis guide](docs/analysis_and_boundary_checks.md); plotting implementation [wetland_analysis_checks.py](wetland_analysis_checks.py), `plot_checks` | Existing 60 s trajectories for response/safety; exact shared 120 s nodes for time-aligned differences; local NPZ arrays remain excluded |

Raw source responses, processed forcing and simulation NPZ files remain in local `B/work/` and `B/results/`. The public tables and recorded source queries establish traceability; they do not claim that the current fresh clone includes every artifact required to execute B. No new observations, parameter fits, simulations or optimization results were introduced by this documentation update.
