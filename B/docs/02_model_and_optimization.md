# HW B.2: peak-stage optimization — implementation and evidence status

**Historical scope.** The 26-variable, 20-residual implementation and numerical results below belong to the first joint gate/cover study with active dam-to-wetland spillways. The final B analysis uses fixed gates and the no-dam-spillway model, with `v = [rho_1, rho_2]`, four enumerated synthetic category pairs and 18 OPT residuals; G5 remains an RHS/Runge–Kutta applicability guard only. See the [B overview](../README.md), [data sources and assumptions](../data_sources_and_assumptions.md) and [wetland-only report](12_wetland_only_optimization.md) for the accepted formulation and results. This historical account retains its original method and execution evidence.

**CEE690, Jiaxuan Zhang. Reviewed 2026-10-01. Status: G1–G6 approved for the first exploratory hypothetical case; a valid and safety-feasible replacement reference has run, and the bounded official OPT search and refined verification are complete.** Final physical results, verification and termination belong in the [execution report](09_exploratory_results.md). This is not local calibration or engineering certification. A.2 remains reproducible and its historical outputs are unchanged.

## B.1 formulation and implementation

Both pages of the submitted *SMCR-Jiaxuan Zhang-HW B.1-2026.pdf* were read visually. Page 1 defines gate histories g1(t), g2(t), wetland choices D_i=(tau_i,rho_i), and **min J(theta), J=max over [0,T] H_r(t)**. Page 2 defines gate/cover/category bounds, the dynamics, nonnegative storage, soil capacity and four upstream safety levels. It does not supply numerical safe levels. The current task explicitly adds fixed initial conditions, evidence handling, numerical applicability checks and an independent finer-grid check. These are implementation requirements from the current task, not invented quotations from B.1. No separate B.2 rubric or professor-specific OPT handout was found in the targeted supplied-file inventory.

`B/optimize_wetland.py` (run from the project root) uses direct shooting: the decisions are **12 opening values for each gate plus two covers (26 continuous variables)**. The original A.2 knot times and piecewise-linear interpolation are retained. The original schedule is tested first and retained as the documented failed baseline check; the replacement reference uses the declared uniform scaling described below. The knot grid is a numerical approximation, not an observed operating rule. All four approved synthetic category pairs are enumerated with their paired resistance/ET endpoints intact. `type_a/type_b` remain synthetic labels; no Milwaukee species mapping is inferred. The same explicitly labeled replacement-reference vector initializes each pair, preserving the original covers and using 0.10 times the original gate values. No optional gate-rate limit, budget constraint, new state, terminal-storage condition or downstream hard flood threshold is added.

Each evaluation uses fresh parameter arrays, the same fixed initial state, geometry, forcing definition and horizon. The reference case and any candidate use the **same revised B.2 gate model**. A.2 results are never the denominator of a claimed B.2 optimization benefit. Numerical conditioning is `(peak stage − 97 m)/(1 m)` in the approved assumed-datum configuration; this fixed shift/positive scale does not change the objective. A future local-datum configuration must document its own consistent origin. All reported peaks remain metres in the declared datum. No percentage reduction of absolute elevation is reported.

## Official modules actually inspected and executed

The existing official checkout is **multivarious 1.2.1**, revision **06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02**. It was neither updated nor edited. The inspected routines are:

```python
nms(func, v_init, v_lb=None, v_ub=None, hyp_in=None, consts=1.0)
ode4u(odefun, time, x0, u=None, c=None)
```

The objective callback returns `(finite_scalar, fixed_constraint_vector)`, with feasibility `g <= tol_g`. The implemented derivative-free course routine is official `nms`, suitable for exploring nonsmooth grid extrema without differentiating category labels. It is not a guarantee of global optimality. `B/official_opt.py` wraps the actual routine; it does not copy or replace its algorithm. The full [source audit](08_official_opt_audit.md) records all settings, return order, history and stopping behavior.

Important verified details:

- Official NMS rewrites zero bounds; the adapter maps each declared physical search box to [-1,1] to preserve its exact bounds (the full box is [0,1]). It also moves a start to the inner part of its box if necessary. The driver evaluates the exact reference-valued starting point separately and reports that distinction.
- The actual return order is `(v, penalized_f, g, history, n_evaluations, iteration)`, differing from the source docstring's return order.
- The module penalty is `f + penalty * sum(g[g > tol_g])**exponent`. Both f and residuals are dimensionless before that calculation. Raw peak stage, module penalty score and independently assessed feasibility remain separate.
- Each of the four category pairs receives up to 40 official NMS callbacks over the full physical [0,1] box, followed by up to 140 in an explicitly narrower numerical neighborhood of the reference. The reference seed is evaluated independently once per category pair before the two phases because official NMS can move it to the interior of its box. Local gate bounds are `0 <= g_k <= min(1, max(0.02, 2*g_reference,k))`; cover bounds remain [0,1]. For this nonzero replacement schedule the local upper gates equal twice the reference values. This is numerical local exploration, not a new physical constraint; the full-box phase retains the complete physical bounds. There is no convergence or global-optimum claim. The hard shared budget remains 800 ODE calls and 900 s, including reference trials and verification; three calls are reserved for fresh candidate, refined candidate and refined reference. Wall-clock checks cannot interrupt an in-flight ODE call.
- Formal NMS convergence is checked only after more than n² iterations: over 676 for 26 variables. A modest budget may end much earlier by stall/budget. Such a result may be a best feasible point found, but must not be called converged automatically.

The actual official Example A from `examples/opt_example.py` was run first with seed 20261001: **66 callbacks, 21 iterations**, raw objective about **0.018410663**, with the direct-example maximum constraint **0.0009287613** within that example's **0.001** tolerance. This is an optimizer smoke example, not a wetland result. The smoke suite passed 12 checks. Its logs are local at `B/results/official_opt_smoke_direct.json` and `official_opt_smoke_adapter.json`.

## Constraints, scaling and failures

The callback has **20 fixed, labeled residuals**. The original 18 remain: eight nonnegative storage, one soil upper bound, four upstream safe levels, four applicability limits and one integration-status residual. Two additional approved physical gate-head/opening residuals enforce hydraulic applicability, one per branch. Native gate/cover bounds are separately audited on the reconstructed histories. No state is an optimization variable.

- Storage: `-min S_j / Vref_j`; soil: `(max S_m − C_m)/Vref_m`.
- Upstream safety: `(max H_di − H_di_safe)/Href` and `(max H_wi − H_wi_safe)/Href`, separately for both branches.
- Applicability: groundwater below `H_g_valid_max`, `soil_bottom` and both wetland beds, and river stage below `H_r_valid_max`. Original geometry/order checks remain in place. The approved explicit gate sill replaces the legacy gate-threshold coordinate for geometry checks only.
- Gate hydraulic applicability: `max(1.25*a_i-H_i)/Href <= 0` on open, wet, forward-flow portions, with inactive closed/dry/adverse cases contributing zero. Here `a_i=gate_max_opening_i*g_i` and `H_i=H_di-z_sill_i`. The RHS rejects a violating intermediate RK stage before unsupported flow is evaluated; this physical condition is not loosened by the numerical strict-bound margin.
- Volume scales are fixed to the positive A.2 initial storage vector, solely as numerical scales; **they are not storage limits**. Href=1 m. Feasibility tolerance is 1e-6, equal in the optimizer and independent audit.
- A 0.001 m **numerical strict-bound margin**, greater than the 1e-6 m tolerated residual, keeps accepted sampled applicability levels strictly below their physical domain bound. It is not freeboard or a safety recommendation. Initial values are included; the RHS also checks intermediate RK states against the un-margined applicability domain.
- Spill crests remain flow-onset levels. No dam or wetland spill-onset storage is imposed as a hard ceiling. `H_flood` is not a constraint; B.2 flood-duration reporting stays disabled until its optional reporting threshold is explicitly defined.

Extrema use the full integration grid and carry the worst time, signed margin, unit, fixed scale, threshold/source and tolerance into the constraint audit. Grid equivalence is not a continuous-time guarantee. The candidate is re-evaluated without cache and both reference/candidate are checked on a 60 s grid against the 120 s search grid, using the same forcing knots and gate interpolation. Both peak differences must pass the separately declared 0.001 m numerical check. Failed verification prevents a feasible/optimized label.

Expected storage/domain errors, nonfinite ODE results, failed mass balance and explicitly unsupported gate regimes are rejected with finite **numerical sentinel** objective=1000 and 20 residuals=1000. They are not simulated stages or safety thresholds; logs leave the physical peak blank and record the reason. Programming errors escape instead of being swallowed. The user has explicitly accepted this restricted hydraulic domain for the first exploratory case. A rejected candidate carries its failure diagnostic and blank physical peak; the search continues. It is not evaluated through an unsupported law, and rejection does not prove physical infeasibility under a different, unimplemented hydraulic method.

A small cache is scoped to the complete configuration hash, category pair and exact decision vector. Cache hits are logged separately from ODE executions. The best independently feasible raw peak is tracked across all evaluations, not selected solely from the optimizer's final iterate. Budget exhaustion and no feasible point found are distinct from mathematical infeasibility.

## Targeted gate revision

Only Q_G1 and Q_G2 change in the B.2 wrapper. The old A.2 rating is `g * min(Q_Gcap, positive_head/R_G)`. The new helper uses **a=gate_max_opening*g**, then:

- Free sluice: `Cd_gate_free * gate_width * a * sqrt(2*gravity*H_energy)`.
- Fully submerged: `Cd_gate_submerged * gate_width * a * sqrt(2*gravity*(H_pool − H_wetland))`.
- Intermediate submergence: the explicit USACE 7.0 weighted treatment, documented in [gate_hydraulics.md](07_gate_hydraulics.md), following the 6.6 head definitions and original manual checks.

H_energy is energy above the sill, approximated by pool depth only under an explicitly approved negligible-approach-velocity-head assumption. Old z_G is not treated as surveyed sill geometry. No extra g multiplier, R_G, unsupported discharge cap or global `head_flow()` change remains in the B.2 gate flux. Each gate flow is computed once, subtracted from its dam and added to its wetland. Independent spillways remain active with gates closed.

Closed, dry and adverse-head cases explicitly return zero under the retained one-way routing assumption. The documented submergence blend is implemented, including its approximately 0.501% lower-end jump from the source's rounded 0.67 criterion; it is not silently smoothed. **Low upstream depth, partially wetted/weir openings and the weir–gate transition are not implemented.** The user has approved the explicit additional condition `H_energy >= 1.25*a` on open, wet, forward-flow gates. Widths, openings, sills and coefficients in [01_study_configuration.md](01_study_configuration.md) and [b2_config.json](../b2_config.json) are accepted fixed engineering assumptions, not local measurements. Feasible wetlands satisfy `H_w1<=98.90<100.00=sill_1` and `H_w2<=98.80<99.80=sill_2`; wetland-induced submergence is excluded from every safety-feasible trajectory under the associated-wetland tailwater assumption. Low upstream head is the relevant remaining domain issue. The missing copied “svg” equation has not been reconstructed.

The helper/wrapper passed **22 fixture tests**, including dimensions, two distinct coefficients, datum invariance, closure/dry/adverse conditions, sourced transition behavior, exact conservative transfer, unchanged other outputs, preserved spillways, and independence from old R_G/Q_Gcap. Fixture values are arithmetic tests only. The original schedule was evaluated first under the revised B.2 model and found hydraulically inadmissible. Its first encountered failure is saved with branch, time, upstream depth, opening height and signed physical residual. The declared feasibility search then tests uniform factors `[0, 0.25, 0.10, 0.05]` on all original gate knots, preserving covers, classes, knot times and fixed parameters; it selects the largest hydraulically admissible tested factor. Factor 0.10 was selected and the resulting replacement reference is numerically valid and safety-feasible. This is a newly prescribed open-loop reference, not state-dependent clipping and not the original baseline. The actual-helper gate-rating diagnostic and replacement trajectory are generated before official OPT. A hydraulically valid but safety-infeasible result would remain a distinct outcome from an invalid trajectory.

## Milwaukee context and evidence

The single [evidence register](../tables/evidence_register.csv) covers every `class c` attribute plus derived geometry/initial state, structural assumptions, safety levels, gate settings and numerical settings. Its source IDs resolve to [source_catalog.json](../source_catalog.json), containing full agency/title/date/version/URL/locator/access/status metadata. A2_SOURCE records where an illustrative value came from; it is **not** numeric local evidence. Potential data sources are separated from sources actually supporting a value.

The approved first case is a regional-event-weather-forced hypothetical0.66km²installation in Underwood context, **not** the full46.879km²gauge catchment. The actual September2–15observations and separateSeptember15–18USGSrecession data have now been downloaded/processed. See [study_configuration.md](01_study_configuration.md) for the concrete state/area/boundary mapping, fixed hypothetical containment levels, assumed gate dimensions, synthetic response classes and inherited fixed assumptions; [data_sources.md](03_data_sources.md) records actual availability and suitability.

The accepted exploratory forcing has4,321120snodes,145.3mmquantifiedKMKErain and21166.81Wh/m²POWERsolar, with KMWCtemperature/wind and fixedRH/pressure derived from its records. Five eventtrace intervals are represented by a flagged zero lower bound, accepted for this first exploratory case only. The selected3-hourtrace record resolves a strict hourly data gap at coarser resolution. The prepared-input manifest retains QC, time conventions, checksums, rainfall/energy conservation, datum andQ_extaccounting. It records the exploratory approval separately from data provenance. The replacement reference uses these processed regional event-weather data; the retained initial soil storage is an assumed condition, not an observed reconstruction of September 2022.

FAO56 supports specific vapour-pressure relationships and several physical constants, not the full reduced evaporation model or arbitrary transfer/vegetation coefficients. SSURGO conductivity is not automatically lumped K_u; available water capacity is not saturated soil storage; Manning n is not R_w. Wetland classification does not identify resistance/ET endpoints. DNR Standard1001's reviewed 10/07 document explicitly excludes flood-control application, so it supplies none of the four safety limits.

## Approval and preflight status

The user has accepted G1–G6 for this first exploratory case subject to baseline, applicability, forcing and result checks. Preflight regenerates [approval_gaps.csv](../tables/approval_gaps.csv) and reports readiness when the approved configuration, declared fixed parameters and forcing manifest agree. The gate restriction is an explicitly approved physical modeling condition, not a hidden numerical filter. This approval does not establish a calibrated Underwood model, real species recommendation or engineering safety certification.

## Commands, outputs and execution record

Use the existing A.2 environment and official checkout described in [README.md](../../README.md); no package installation was performed for B.2.

```sh
python3 -B -m B.official_opt_smoke
python3 -B -m unittest B.test_gate_hydraulics B.test_workflow B.test_data_preparation -v
python3 -B -m B.prepare_b2_data --evidence-register
python3 -B -m B.optimize_wetland --preflight
python3 -B -m B.optimize_wetland --plumbing-check
# Approved first exploratory case; writes separate B.2 results:
python3 -B -m B.optimize_wetland --run
```

Actually executed: official OPT smoke and its 12 checks; **67 current gate/workflow/data tests**; evidence-register generation; preflight; one unchanged original A.2 six-day ODE evaluation. The latter took approximately 0.66 s in this environment and reproduced peak stage **97.87940653911606 m**, peak outflow **1.3578956737526984 m³/s**, and maximum integrated balance residual **0.07673464775143657 m³**. These are A.2 plumbing/regression results only. They do not measure revised-gate cost or B.2 feasibility. Its output is in `B/results/a2_plumbing_check.json`; original A.2 files were not overwritten.

Preflight now returns readiness for the approved configuration. The original-schedule hydraulic test and declared reference feasibility search have run, producing a valid, safety-feasible 0.10-scale replacement reference. Official optimization and refined verification are complete; actual physical metrics, residuals and termination information are recorded in [09_exploratory_results.md](09_exploratory_results.md). No penalized score is reported as a river peak, and no optimization completion or convergence is inferred from approval or code checks.

Actual outputs include `B/plots/B2_0_gate_rating.png` and five river/control/safety/history/gate-applicability figures, with decision/performance/constraint/category tables under `B/tables/`. Callback histories label evaluations, including cache hits, rather than iterations. Raw arrays, configuration/data manifests and logs go to ignored `B/results/`; raw source material remains in ignored `B/work/`. Evidence/approval/data-availability tables are kept distinct from simulated performance, decisions, residuals and callback diagnostics. The current execution check is local `B/results/exploratory_verification.json`. Earlier `verification.json` and `unit_tests.log` are historical records from the pre-approval stage. Historical protection checks recorded unchanged scientific files. The current A.2 scientific files and outputs remain protected; README navigation and B.2 documentation are updated separately. The Git index, HEAD and remote are unchanged, and the official dependency has no tracked source edits.

Every future report must distinguish numerical/model feasibility from real-world safety validation. It will report end-of-horizon storage, late peaks and any rising tail, without silently imposing a terminal constraint. Event-specific optimized gates have perfect foresight and are not a real-time operating policy. A finite-grid, budget-limited search can at most report the **best feasible solution found under the stated model, data and parameterization**.

## References and traceability

All external references accessed **2026-10-01**. Full records (including actual dataset versions and field/section locators) are in the source catalog; data-specific scope/access limitations are in the data-sources note.

- **B1:** Jiaxuan Zhang, *SMCR–HW B.1*, due2026-10-01, pp1–2; user-supplied submission, not redistributed.
- **OFFICIAL_OPT:** Henri P. Gavin, *multivarious*, version1.2.1, [pinned source](https://github.com/hpgavin/multivarious/tree/06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02); NMS, ode4u, settings/history and ExampleA inspected locally.
- **G1/G2/G4/G9:** USACE HEC, *HEC-RAS Hydraulic Reference Manual*, [6.6 sluice section](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.6/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways/sluice-gate), original CPD-69v6.0 (May2021), and [7.0 explicit transition](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/7.0/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways). Complete conventions and source discrepancy are documented in the gate note.
- **G5–G8/G10:** USBR, *Water Measurement Manual*, 3rd edition, revised reprint2001, [Chapter9](https://www.usbr.gov/tsc/techreferences/mands/wmm/chap09_01.html); measuring-orifice geometry and coefficient conventions compared, not indiscriminately adopted.
- **FAO56:** Allen, Pereira, Raes and Smith, *Crop Evapotranspiration*, FAO Irrigation and Drainage Paper56,1998, [Chapter3](https://www.fao.org/4/x0490e/x0490e07.htm), Eq8/Eq11/Box8 and unit conversions.
- **USGS_04087088 / USGS_WBD_UNDERWOOD:** USGS [Underwood Creek at Wauwatosa](https://waterdata.usgs.gov/monitoring-location/USGS-04087088/) and WBD HUC12 record; metadata/context only, not modeled observations.
- **WIDNR_1001:** Wisconsin DNR, *Wet Detention Pond Technical Standard1001*,10/07 version, [document](https://dnr.wisconsin.gov/sites/default/files/topic/Stormwater/1001WetDetentionPond.pdf); its scope does not establish flood-control safety levels.
