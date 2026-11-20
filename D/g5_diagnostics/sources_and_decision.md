# Stage D gate sources and hydraulic decision

## Decision and scope

No low-head discharge closure is implemented. The existing hypothetical gate specification does not identify the crest geometry, coefficient or transition algorithm needed to extend its discharge law defensibly. Stage D therefore evaluates a documented, common **time-only operating schedule** with the original Stage B hydraulic function. This is an operating-domain evaluation, not an all-regime hydraulic correction or a claim of feasibility under arbitrary stochastic weather. Pilot findings and any acceptance decision belong in the diagnostic results, not in this source assessment.

A/B/C equations, configurations and saved results remain historical evidence. No independent dam spillway is introduced. Low-head flow through a gate would be another regime of the same passage, not a separate overflow path.

## Existing specification

| Property | Gate 1 | Gate 2 | Status |
|---|---:|---:|---|
| Rectangular gate width | 1.00 m | 0.80 m | Accepted hypothetical design |
| Maximum opening height | 0.30 m | 0.30 m | Actual opening is `a = 0.30 g(t)` |
| Sill elevation | 100.00 m | 99.80 m | Equal to the corresponding zero-storage pool bed |
| Free-flow coefficient | 0.60 | 0.60 | Dimensionless orifice coefficient |
| Submerged coefficient | 0.80 | 0.80 | Separate dimensionless orifice coefficient |
| Gravity | 9.80665 m/s² | 9.80665 m/s² | Existing setting |
| Upstream energy approximation | Pool water surface | Pool water surface | Approach velocity head neglected |

The gates are interpreted as vertical-lift gates discharging into channels. Upstream head is measured above the sill. These are engineering assumptions, not measured Milwaukee gate characteristics. Pool plan area does not establish the approach channel cross-section.

Under the accepted eight-state tailwater assumption, the wetland safety limits (98.90 and 98.80 m) lie below their respective gate sills. Wetland-induced submergence is consequently excluded from safety-feasible eight-state trajectories. The six-state model uses river tailwater and requires its own check.

## What remains supported

The original function returns zero for a closed gate, a dry donor above-sill section, or nonpositive forward head. For positive forward flow, `H/a >= 1.25` is supported: free flow at submergence ratio `s <= 0.67`, the existing free/submerged blend for `0.67 < s < 0.80`, and submerged orifice flow at `s >= 0.80`. Positive-flow states with `H/a <= 1` or `1 < H/a < 1.25` remain explicitly unsupported.

The original rounded submergence boundary is preserved. At `s = 0.67`, the right-hand transition limit is `sqrt(0.99)` times the free-flow value, approximately 0.501% lower. This existing small jump is reported rather than silently smoothed. It is distinct from the missing low-head regime. A test documenting it must not be presented as evidence of exact continuity.

## Source findings

The [USACE HEC-RAS 7.0 gated-flow technical reference](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/7.0/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways) identifies the weir, transition and orifice regimes. Its overview describes iterative head blending and subsequently describes linear discharge transition. These descriptions alone do not specify a unique direct `Q(H, a)` implementation for this lumped model.

The [USACE low-flow guidance](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/7.0/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways/low-flow-through-gates) uses `Q = Cw L H^(3/2)` with crest-dependent coefficients and submergence treatment. Broad-crested, ogee and sharp-crested options are physically different. From dimensional analysis, SI `Cw` has units `m^(1/2)/s`; the dimensionless orifice coefficient cannot replace it. Matching a guessed coefficient to another branch at one point would establish continuity, not physical applicability.

The [USACE inline-structure user manual](https://www.hec.usace.army.mil/confluence/rasdocs/rasum/7.0/entering-and-editing-geometric-data/inline-structures-dams-weirs-and-gated-spillways) requires separate weir choices. Ogee treatment needs design-head information; sharp-crested coefficient methods can require approach height and width. This page also gives a 1.0–1.1 transition statement, differing from the technical reference's 1.0–1.25. That discrepancy is not grounds to relax the project's existing 1.25 restriction.

## Evidence and decisions needed for an extension

| Missing item | Required evidence or explicit engineering decision |
|---|---|
| Crest class and physical profile | Specify whether the gate sill is broad-crested, sharp-crested or ogee, with the relevant streamwise geometry and edge/contraction conditions. Gate width and sill elevation alone are insufficient. |
| Applicable weir coefficient | Supply a traceable coefficient or validated coefficient method in a verified unit convention, supported by the selected geometry and approach conditions. |
| Approach geometry/losses | Supply the dimensions and assumptions required by the selected method; do not infer them from pool storage area or assume zero loss implicitly. |
| Weir-gate transition | Select a documented algorithm consistent with that geometry and the intended head convention, then verify limiting behavior and continuity. Resolve the source ambiguity before implementation. |
| Low-head tailwater effects | Establish the applicable crest-specific correction wherever tailwater can submerge the passage; check the river-tailwater case separately. |

Until those decisions are supported, the fallback must use one explicitly recorded schedule shared across both physical models and all selected scenarios. No state-dependent opening clipping, altered geometry, raised safety limit or dropped scenario is permitted. A schedule can pass G5 yet exceed dam safety limits; both checks are required. Failure of tested schedules establishes only those failures, not proven infeasibility of all schedules.
