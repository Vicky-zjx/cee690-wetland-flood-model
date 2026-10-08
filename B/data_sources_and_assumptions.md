# B data sources and assumptions: fixed-gate wetland design

The accepted B analysis optimizes the two wetland designs under one fixed gate schedule. Water passes from each upstream pool to its associated wetland **only through the modeled gate**. The project is a hypothetical engineered 0.66 km² installation driven by regional event weather, not a calibrated model of the Underwood Creek gauge catchment or a certified local design. The [B overview](README.md) presents the final analysis and results.

The active scientific model is [no_spillway_model.py](no_spillway_model.py), and the optimization driver is [wetland_only_study.py](wetland_only_study.py). The fixed physical inputs are recorded in [no_spillway_config.json](no_spillway_config.json); [wetland_only_config.json](wetland_only_config.json) freezes that configuration and the saved reference gate schedule by checksum. The study loader verifies both hashes before execution. Saved provenance paths retain their original labels and are resolved by [project_paths.py](project_paths.py). These source records remain immutable.

Current optimization results, tables and figures are under `B/results/wetland_only/`, `B/tables/wetland_only/` and `B/plots/wetland_only/`. The [optimization report](docs/12_wetland_only_optimization.md) records actual trajectories, feasibility, search termination and timestep checks. The separate fixed no-wetland comparison is documented below and in the [comparison report](docs/13_no_wetland_counterfactual.md).

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
