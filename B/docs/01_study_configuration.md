# First exploratory B.2 study configuration

**Historical scope.** This document records the first joint gate/cover study, including its active dam-to-wetland spillways and original execution conditions. The final B analysis fixes all 24 gate knots, uses the no-dam-spillway model, and optimizes only two wetland covers while enumerating four synthetic category pairs. Its current configuration, assumptions and 18-residual/G5 roles are described in the [B overview](../README.md), [data sources and assumptions](../data_sources_and_assumptions.md) and [wetland-only report](12_wetland_only_optimization.md). The historical account below does not define the final study.

**Approved for the first exploratory case, 2026-10-01 (G1–G6).** The user authorizes a hypothetical engineered system driven by regional event-weather data. This approval does not establish a calibrated Underwood Creek model or local engineering safety certification. Execution results are reported separately; approval alone is not a successful baseline or optimization. The objective remains min max_t H_r(t), with the submitted B.1 bounds, eight balances, fixed initial conditions, storage/soil constraints and four fixed upstream maximum levels. No safety level was selected from a simulated peak or adjusted to make a run feasible.

## One mapping and its claim boundary

Use the approved **regional-event-weather-forced hypothetical 0.66 km² drainage installation in the Underwood Creek geographic context** over 2022-09-09 00:00 UTC through 2022-09-15 00:00 UTC. This is an unsited analytical installation, not a delineated real parcel, a reconstruction of County Grounds, or a prediction at USGS04087088. Its effective surface areas are explicit design assumptions. The actual gauge catchment remains the 18.1 mi² (46.879 km²) observation context; it is not scaled into these eight states.

The eight states map as follows:

- S_u: 0.600 km² upland surface runoff storage.
- S_m: soil water beneath that same upland; it receives infiltration, not a second rainfall footprint.
- S_d1/S_d2: hypothetical upstream pools of 0.010/0.008 km².
- S_w1/S_w2: associated downstream **proposed embanked/contained wetland units** of 0.012/0.010 km²; these are designed storage cells, not unsurveyed natural wetlands.
- S_g: shared shallow aquifer below 0.640 km², an overlapping subsurface footprint.
- S_r: a hypothetical local receiving channel with 0.020 km² effective water surface, not the entire Underwood mainstem.

The disjoint surface footprints sum to 0.660 km². Upland runoff splits between both dam pools and the retained direct river bypass. Each gate AND its independent dam spillway transfer to its associated wetland. Both wetlands drain/overflow to the local channel; soil and groundwater paths stay intact. Original constant-area storage relationships and topology are preserved.

Use event weather over the hypothetical footprints. Set **Q_ext=0 m³/s** because no outside surface-drainage area is included. The observed USGS hydrograph is contextual evidence only: it is neither Q_ext nor a calibration/validation target for S_r. Use assumed model-datum boundaries **h_d=97.15 m**, **h_g=97.15 m**, and **B=0** for a maintained, unobstructed local outlet. These are explicit design/boundary assumptions, not observed Underwood levels or blockage. The original synthetic h_d/h_g pulses and Q_ext pulses are not used. All modeled elevations remain in one assumed common datum; do not offset or mix them with the observed NAVD88 series.

Downloaded official context includes the WBD HUC12 040400030404 polygon (agency area50.72km²), NLDI station point43.0547222,-88.0461111 and NLDI COMID12163882. The split-at-gauge basin endpoint failed with HTTP502; a whole-catchments simplified polygon was retrieved but is approximately50.56km² by a spherical diagnostic, so it is explicitly not the exact46.879km² station catchment. The structured station coordinate agrees with NLDI, whereas legacy narrative latitude differs by17arcseconds. A field-sited or full-catchment prediction would still require the precise boundary, facility/reach survey and a different evidence-supported geometry/rating configuration. No arbitrary rectangle has been labeled a real catchment.

## Four approved hypothetical maximum levels, fixed before any run

Specify **assumed hypothetical settled crest elevations** for the dam embankments and wetland containment tops at dam1=101.55m, dam2=101.35m, wetland1=99.25m and wetland2=99.15m in the model datum. Reserve **0.35m physical freeboard**. This gives:

```text
H_d1_safe = 101.20 m       H_d2_safe = 101.00 m
H_w1_safe =  98.90 m       H_w2_safe =  98.80 m
H_safe = specified settled top - 0.35 m
```

Thus the maximum modeled depths are1.20m in both pools and0.90m in both wetland cells. Each specified top is0.65m above its existing spill-onset level; each maximum water level is0.30m above that onset. Spill-onset elevations remain flow thresholds, not safety ceilings. The0.001m numerical applicability margin is unrelated to the0.35m physical freeboard.

**Evidence category: approved engineering assumptions for this exploratory design, with a limited Wisconsin pond analogue for the freeboard concept.** The absolute crest elevations above are separately assumed; they are not supplied by NRCS. Wisconsin NRCS Pond378, November2024, printed378-CPS-5 (PDFp87/zero-index86), specifies at least1ft high-water-to-settled-top freeboard for certain pond classes and2ft auxiliary-crest-to-top separation when drainage exceeds20acres. The proposal's0.35m and0.65m exceed those two dimensions (0.3048m and0.6096m). The document includes flow detention but requires site/hazard classification and other design checks. It does not supply these absolute elevations or certify these hypothetical facilities. Applying the same reserve to constructed wetland containment is an explicit analogy, not a natural-wetland safety standard. [WI NRCS378](https://www.nrcs.usda.gov/sites/default/files/2024-12/450-11%20FOTG%20Notice%20WI-108%20Combined.pdf).

These are fixed hypothetical containment limits, not demonstrated ecological depth tolerances. There is no surveyed crest, exposed-asset elevation, plant flood-tolerance study, or real facility safety certification here. Any actual-site safety claim would require those records and appropriate design checks. The user has accepted these limits for this first exploratory case only; they must not be raised to obtain feasibility.

## Fixed gates and the method decision

Use the approved assumed gate widths[1.0,0.8]m; maximum opening heights[0.30,0.30]m; bottom-mounted sills[100.00,99.80]m; gravity9.80665m/s²; Cd_free=[0.60,0.60] and Cd_submerged=[0.80,0.80]. Width scales with existing pool-area ratio; travel is one third of the0.90m bed-to-spill depth. These are design dimensions, not observations. Free Cd0.60 is the midpoint of the USACE typical0.5–0.7 range; submerged0.80 is the separately stated typical value, with no verified range supplied. The pool-depth approximation neglects approach velocity head.

The user has approved the following **additional hydraulic-applicability constraint**: for each open, wet, forward-flow gate, require `H_pool - z_sill >= 1.25 * gate_max_opening * g`. At full opening this requires0.375m depth; at initial0.30m depth it permits g≤0.80. It must hold during the entire integration, not only initially. No opening is clipped and no automatic feedback is added. Closed/dry/adverse routing remains the existing one-way assumption. This explicit physical modeling constraint narrows the admissible control/state domain and is additional to the unchanged B.1 objective and safety limits; it is not merely a numerical check. Test the intended original baseline gate schedule first. If it violates this condition, report the first violation with branch, time, upstream water depth, opening height and residual, without clipping the opening, changing geometry or loosening the restriction. Any admissible replacement must be identified as a replacement reference schedule and its search method reported; a feasibility search is not the original baseline.

Under the model's gate-tailwater assumption, the associated wetland surface supplies the tailwater. Every safety-feasible trajectory satisfies `H_w1 <= 98.90 < 100.00 = sill_1` and `H_w2 <= 98.80 < 99.80 = sill_2`. Wetland-induced gate submergence is therefore excluded from all safety-feasible trajectories under this assumption. The remaining relevant coverage issue is low upstream head. The existing free/submerged/transition helper is retained without expanding irrelevant regimes. A complete low-head/weir law remains unresolved: crest type/coefficient, tailwater correction and the precise iterative weir–gate transition are not determined by width/opening/orifice coefficients. [gate_method_proposal.md](06_gate_proposal.md) states the verified source limitations and the exact alternative information needed. The approved restriction avoids evaluating an unsupported low-head law; no new closure is invented.

## Vegetation and remaining parameters

Retain the approved original two **synthetic hydraulic/ET response classes**, admissible for both branches, with their paired endpoints kept together:

```text
                    R_veg [s/m²], branch1/2       K_veg [-], branch1/2
class type_a        [2.8, 3.0]                   [1.10, 1.10]
class type_b        [4.0, 4.5]                   [1.25, 1.25]
bare endpoint       [0.8, 1.0]                   [1.00, 1.00]
```

They represent lower/higher resistance and ET response assumptions, **not locally identified species or measured plant properties**. Covers remain the B.1 decision variables in[0,1], with reference[0.65,0.60]. These classes are not a Milwaukee-specific vegetation recommendation. A quantitative species mapping is not defensible from the downloaded data. It would require site vegetation/stem-density and depth-dependent hydraulic measurements plus an ET relation compatible with this reduced demand model. FAO56's wetland Kc values apply to FAO Penman–Monteith ETo, not automatically to this E_star; Manning n is not numerically equal to the lumped R. Neither is used to disguise the existing endpoint choices as measurements.

Retain the approved other A.2 physical constants, effective areas, beds, spill/drainage ratings, soil/groundwater coefficients and initial storages as **fixed, uncalibrated assumptions for this hypothetical installation**. Their exact values are enumerated in the existing evidence register's B2_new_value column, with unused old gate/synthetic-weather quantities marked accordingly. Important examples are soil capacity120000m³ (0.20m water-equivalent over the upland), S_fc60000m³ (0.10m), K_u1e-6m/s (3.6mm/h), groundwater S_y0.20, and initial soil storage90000m³ (75%capacity). These are not SSURGO estimates or measured antecedent moisture.

Fixed x0 remains `[600,90000,3000,2500,3000,2400,409600,7000] m³` in the original state order. Antecedent rain is reported but does not uniquely determine eight storages without calibration or observed states. No warm-up-generated condition is quietly substituted. The original120s integration grid, six-day horizon, gate knot times,26continuous decisions and official optimization budget remain unchanged.

## Actual data, accepted assumptions and execution conditions

See [usgs_event_findings.md](04_usgs_observations.md) and [weather_event_findings.md](05_weather_observations.md) for actual downloaded records, dates, QC, missingness and proposed transformations. The prepared forcing and explicitly declared boundary/initial assumptions are accepted for this exploratory case only. No observation is both a prescribed boundary and an independent validation target. The original six-day horizon is retained; observed recession beyond September15 is reported rather than silently changing T or adding a terminal constraint.

The user has accepted G1–G6 subject to the execution checks below; [approval_gaps.csv](../tables/approval_gaps.csv) is the machine-generated status table. First test the original gate schedule with the revised event forcing and approved fixed configuration. If hydraulically inadmissible, report the failure and seek a clearly labeled admissible reference schedule under the unchanged limits. Run that reference before OPT and assess its safety feasibility. A valid trajectory that exceeds a safety limit is distinct from a hydraulically invalid trajectory. Use exactly the same hydraulic model, forcing, x0 and fixed parameters for the reference and official OPT.

During OPT, reject an invalid hydraulic candidate as infeasible through the verified OPT interface, record its diagnostic, and continue the search. Stop that candidate before evaluating unsupported hydraulics; do not invent the rest of its trajectory. Report physical river-peak objectives separately from any optimizer penalty, all relevant constraint residuals, and solver termination information. If no feasible solution is found, report that observed search outcome without claiming proven infeasibility. Never raise limits or change fixed parameters to manufacture success.

## Accepted exploratory forcing and required verification

`B/work/event_data/proposed_forcing_120s.csv` contains4,321 finite rows over the fixed six-day window. Its companion manifest records the case approval separately from source provenance. KMKE rainfall contributes145.3mm quantified depth; five event trace intervals retain an unquantified positive amount, treated as zero only for the accepted exploratory lower-bound input. KMWC supplies T/wind; accepted constant `r_H=0.7084916611553665` and `p_0=98692.51171875Pa` are event-time means derived from its accepted observations. POWER contributes21,166.81Wh/m² of shortwave over the event, with its UTC label provisionally interpreted as the start of each hourly bin.

Both accumulation inputs use one120s ramp at each interval boundary and an adjusted interior plateau, preserving every source interval's quantified integral while retaining zero intervals. Boundary rates are the smaller adjoining means. Clipped first/last rainfall intervals assume uniform source-interval intensity. This accepted subinterval shape is a numerical assumption, not measured variability; a finer ODE grid must reuse the same prepared knots. Temperature/wind are linearly sampled from bracketing observations, with no extrapolation and at most70minutes between brackets. No wind-height correction is claimed because the returned records do not establish sensor height.

```sh
python3 -m B.work.fetch_usgs_event --offline
python3 -m B.work.fetch_weather_event --process
python3 -m B.prepare_b2_data --event-forcing
python3 -m B.prepare_b2_data --evidence-register
python3 -m B.optimize_wetland --preflight
```

The preparation command does not run the ODE. Before execution, verify UTC alignment, source accumulation-window handling, and preservation of quantified rainfall depth and solar energy. Keep KMKE/KMWC station observations, POWER grid estimates, transfer to the hypothetical footprint, trace-as-zero rainfall lower bounds, and assumed boundary/initial conditions separately labeled. In particular, the retained initial soil storage is an assumed condition, not an observed reconstruction of September 2022. Physical approval flags, parameter overrides and prepared-manifest approval must consistently represent this first exploratory case. The run entry point is `python3 -m B.optimize_wetland --run`; its reference and optimizer results must be reported before B.2 can be described as completed.
