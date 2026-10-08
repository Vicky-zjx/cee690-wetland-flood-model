# Approved exploratory B.2 gate restriction and fixed geometry

**Status: approved by the user for the first exploratory case, 2026-10-01 (G4–G5).** These are hypothetical gates for the approved 0.66 km² installation, not surveyed Milwaukee structures or locally calibrated gates. Dimensions and coefficients are fixed engineering assumptions. This document records the approved method and required checks; actual reference/search outcomes are reported separately. A.2 remains unchanged.

## Approved additional hydraulic-applicability constraint

Retain the documented sluice/orifice helper in [gate_hydraulics.py](../gate_hydraulics.py), with the following explicit additional physical hydraulic-applicability constraint. For each gate and every evaluated time/state, define

```text
a_i = a_max_i * g_i
H_i = H_dam_i - z_sill_i
DeltaH_i = H_dam_i - H_wetland_i
```

Whenever `g_i > 0`, `H_i > 0`, and `DeltaH_i > 0`, require

```text
H_i >= 1.25 * a_max_i * g_i.
```

Equivalently, an open gate under forward wet conditions must have `g_i <= min(1, H_i/(1.25*a_max_i))`. This expression defines an admissibility condition; it is **not** an instruction to clip the gate, add automatic feedback, or change its prescribed schedule. Closed gates, a dry opening, and equal/adverse head retain the explicitly assumed zero-flow routes. The independent dam spillways remain active.

This approved constraint excludes low-head/weir and weir–gate transition states and explicitly narrows the original admissible control/state domain. It is additional to the unchanged B.1 objective, bounds and safety limits, and must not be hidden as a numerical check. Enforce and report it throughout integration; satisfying it initially is insufficient. A rejected hydraulic trajectory has no complete physical objective or downstream safety assessment. A hydraulically valid simulation that exceeds a safety limit is a different outcome. Neither outcome establishes proven infeasibility of the design problem.

For the approved geometry and model tailwater assumptions, the gate tailwater is its associated wetland surface. On every safety-feasible trajectory, `H_w1 <= 98.90 m < sill_1 = 100.00 m` and `H_w2 <= 98.80 m < sill_2 = 99.80 m`. Wetland-induced gate submergence is therefore excluded from all safety-feasible trajectories. This conclusion depends on those modeled tailwaters and does not cover an unmodeled external backwater connection. Prioritize low upstream head; no additional submergence regime is needed for the feasible set.

Within this domain the existing law is complete for the three documented tailwater regimes. With `s = (H_wetland-z_sill)/H`, `A = width*a`, and `gamma` denoting gravity:

| Condition | Q [m³/s] |
|---|---|
| `s <= .67` | `Cd_free*A*sqrt(2*gamma*H)` |
| `.67 < s < .80` | `(1-f)*Cd_free*A*sqrt(2*gamma*3*DeltaH) + f*Cd_sub*A*sqrt(2*gamma*DeltaH)`, where `f=(s-.67)/(.80-.67)` |
| `s >= .80` | `Cd_sub*A*sqrt(2*gamma*DeltaH)` |

The free and submerged formulas follow the [USACE 6.6 sluice section](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.6/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways/sluice-gate) and its linked [submerged-orifice equation](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.6/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways/radial-gates). The explicit transition is from the [USACE 7.0 technical overview](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/7.0/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways). Its rounded .67 boundary produces the already documented approximately 0.501% downward jump; it is not silently smoothed.

## Why an all-regime implementation is not supplied

The complete [official CPD-69 manual](https://www.hec.usace.army.mil/software/hec-ras/documentation/HEC-RAS_6.0_Reference_Manual.pdf), printed pp.8-7 and 8-11–8-12, specifies weir flow for `H/a<=1`, using `Q=Cw*L*H^(3/2)`. Its coefficient and tailwater correction depend on crest geometry: broad-crested and ogee structures use different curves; sharp-crested structures use the Villemonte correction. Width, maximum opening, and two orifice coefficients do not supply those missing choices.

For `1<H/a<1.25`, the manual describes an iterative weighted average of upstream heads computed from both laws. The 7.0 overview repeats this head-averaging description but also calls the process a linear flow transition. Those statements do not uniquely specify the direct-discharge algorithm, its weighting evaluation, and the interaction with crest-specific tailwater correction. Choosing a flow blend, a head blend, or a minimum would therefore introduce an additional modeling decision. No such choice is presented as a fully verified USACE implementation here. The alternative to the restricted proposal is to resolve those details and approve the necessary crest geometry before implementing and testing an all-regime law.

## Approved fixed engineering assumptions

Both gates are hypothesized to be vertical-lift gates over outlet-channel inverts, with linear opening maps and negligible upstream approach-velocity head. The [USACE gate-input documentation](https://www.hec.usace.army.mil/confluence/cwmsdocs/rasum/7.0/inline-structures-dams-weirs-and-gated-spillways-125763625.html) distinguishes sill-referenced channel discharge from a freely jetting orifice, which uses an opening-center reference. This exploratory case adopts the channel-discharge convention. Pool plan area does not establish approach-flow cross-sectional area or prove negligible velocity head; this remains an explicit engineering assumption accepted for this case, not a verified local property.

| Parameter / units | Gate 1 | Gate 2 | Basis, range, and uncertainty |
|---|---:|---:|---|
| `gate_width` [m] | 1.00 | 0.80 | Hypothetical design dimensions; positive mathematical domain, no surveyed range. The 0.8 ratio matches the existing pool-area ratio, giving equal width per pool plan area. |
| `gate_max_opening` [m] | 0.30 | 0.30 | Hypothetical shared actuator travel, one third of each existing 0.90 m bed-to-spill depth; not fitted to a hydrograph. Positive domain; no evidence-based local interval. |
| `gate_sill_z` [m, existing relative datum] | 100.00 | 99.80 | Hypothetical bottom-mounted gates, exactly at their respective zero-storage beds. Each remains 0.90 m below its existing dam spill crest. Old empirical `z_G` is not treated as a surveyed sill. |
| `Cd_gate_free` [dimensionless] | 0.60 | 0.60 | Midpoint of the USACE sluice section's typical 0.5–0.7 range. This is transferable engineering context, not local calibration or an optimization bound. |
| `Cd_gate_submerged` [dimensionless] | 0.80 | 0.80 | USACE's separately stated typical submerged-orifice value. The reviewed gate section does **not** give a numerical interval; no empirical range is invented. Mathematical domain is positive. |
| `gravity` [m/s²] | 9.80665 | 9.80665 | Conventional standard gravity in the [NIST constants table](https://physics.nist.gov/cuu/Constants/Table/allascii.txt), not a local gravity measurement; distinct from normalized gate controls. |
| `neglect_approach_velocity_head` | `True` | `True` | Accepted lumped-pool approximation; not established by the dimensions alone. |

All dimensions and coefficients remain fixed during reference evaluation and optimization. Their selection uses existing storage geometry and general coefficient guidance, not the desired peak, constraint feasibility, or an optimizer's answer. No range for either coefficient is a calibrated confidence interval. The separate `.6` and `.8` choices preserve the manuals' different head/regime conventions; the USBR `.61` measuring-orifice coefficient is not substituted.

The approved restriction requires at least 0.375 m upstream depth above the sill at full opening. At the existing 0.30 m initial pool depths, the allowed opening under forward-flow conditions is at most 0.80. The existing initial controls `.35` and `.45` meet that initial condition. These arithmetic checks establish neither later admissibility nor event feasibility. A declining pool can enter an unsupported regime even when the initial point is covered.

## Required reference and optimization checks

The approved package includes the hypothetical dimensions and independent coefficients above, the channel/sill and negligible-approach-head assumptions, existing one-way routing, and the explicit state-dependent hydraulic-applicability constraint. The original intended gate schedule must be tested first. If it violates the constraint, report the first violation, branch, elapsed time, upstream depth `H`, opening height `a`, and physical residual `1.25*a-H` (positive means a violation). Do not clip controls, change geometry, or loosen the restriction.

If the original schedule is inadmissible, seek a hydraulically admissible reference under the unchanged limits. Clearly label any replacement reference and explain how it was obtained; do not describe a feasibility search as the original baseline. A replacement reference must be simulated and its numerical, hydraulic and safety outcomes reported before optimization proceeds. Reference and optimization must share the same hydraulic model, initial conditions, forcing and fixed parameters.

Within official OPT, an invalid candidate must be returned as infeasible through its verified interface with a recorded diagnostic. Reject the candidate before evaluating unsupported hydraulics, without terminating the entire search or completing an invented trajectory. Report the physical river peak independently of optimizer penalties, constraint residuals and solver termination status. Failure to find a feasible solution must be reported without changing limits or claiming proven infeasibility.

The existing 22 helper/wrapper fixture tests cover the implemented branches and conservative routing; they do not validate these assumed project values or substitute for the required trajectory checks. See [gate_hydraulics.md](07_gate_hydraulics.md) for exact API, source discrepancies, limits, and previously executed checks. All primary sources linked here were accessed on 2026-10-01.
