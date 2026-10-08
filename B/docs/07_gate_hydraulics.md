# B.2 gate-law revision: implemented branches and approved restricted domain

Access/review date: **2026-10-01**. This targeted revision concerns Q_G1 and Q_G2 only. A.2 sources/results, the eight storage balances, other flow laws, and separate dam spillways are unchanged. The copied “svg” equation was unavailable; no claim is made to have recovered it.

The old rating was `g_i * min(Q_Gcap[i], positive_head/R_G[i])`. The new helper uses the independently verified USACE relationships below. It contains neither R_G nor Q_Gcap, no extra opening multiplier, and no unexplained minimum of free/submerged formulas. **It is not a complete all-regime gate model.** The user has approved fixed hypothetical dimensions/coefficients and an explicit additional head/opening applicability constraint for the first exploratory case. Unsupported low-head regimes remain outside that approved domain; passing fixtures alone did not provide this authorization.

## Interface and heads

`GateParameters(gravity, gate_width, gate_max_opening, gate_sill_z, Cd_gate_free, Cd_gate_submerged, neglect_approach_velocity_head)` describes one branch. No numeric field has a project default: omitted/`None`, nonfinite, or invalid values raise `GateConfigurationError`. The head-approximation flag must explicitly be `True`; the user accepted that engineering assumption for this exploratory case; it is not a local measurement.

`sluice_gate_flow(upstream_level, downstream_level, normalized_opening, params)` returns `GateFlow`, with `discharge` [m³/s], `regime`, `opening_height` [m], `energy_head` [m], `differential_head` [m], and dimensionless `submergence`. `GateInputError` identifies invalid numerical inputs. The implementation is scalar, pure, and independent of the ODE solver.

The explicit modeling assumptions and definitions are:

- `a = gate_max_opening * normalized_opening`, with normalized opening in [0,1]. This linear mapping uses opening once; width and maximum opening remain fixed configuration parameters, not optimization variables.
- `E_up ≈ H_pool`: neglect upstream approach-velocity head. The HEC-RAS energy head is then `H_energy = H_pool − gate_sill_z`, not an absolute elevation and not the head above the opening centroid [G1,G2].
- `dH = H_pool − H_wetland`; `s = (H_wetland − gate_sill_z)/H_energy` when H_energy>0. A negative downstream depth simply belongs to the free-discharge side; it is not substituted into a square root.
- `gate_sill_z` must be independently documented in the common datum. The old empirical z_G is not automatically surveyed sill geometry. The integrating configuration must check sill≥pool zero-storage bed so an empty donor cannot release water through a lower fictitious opening.

## Implemented and unsupported regimes

For an open gate with positive upstream depth and forward head:

| Regime | Applicability used by the helper | Discharge / behavior |
|---|---|---|
| Free sluice | H_energy/a≥1.25 and s≤0.67 | `Cd_gate_free * gate_width * a * sqrt(2*gravity*H_energy)` [G1; G2 for definitions] |
| Fully submerged orifice | H_energy/a≥1.25 and s≥0.80 | `Cd_gate_submerged * gate_width * a * sqrt(2*gravity*dH)` [G2, Eq. 3] |
| Partially wetted opening / weir | H_energy/a≤1 | `UnsupportedGateRegime('weir_or_partially_wetted_opening', ...)` [G3,G4,G5] |
| Weir–gate transition | 1<H_energy/a<1.25 | `UnsupportedGateRegime('weir_gate_transition', ...)` [G3,G4] |
| Free–submerged transition | H_energy/a≥1.25 and 0.67<s<0.80 | Explicit weighted law below [G9], consistent with the transition head in G1/G2. |

At H_energy/a=1.25, the gate branch is used. At s=.67 and .80, the corresponding free and fully submerged endpoint is used. Because .80×1.25=1, the submerged branch's tailwater reaches at least the opening top under these assumptions.

The 6.6 sluice page supplies a transition expression involving `sqrt(2*gravity*3*dH)`; the radial page describes a weighted combination with its Eq.3 but does not specify the weights. The **complete official 7.0 overview [G9] explicitly supplies them**, so the implemented transition uses that documented clarification:

```text
fraction = (s - 0.67)/(0.80 - 0.67)
Q_transition = (1-fraction)*Cd_gate_free*width*a*sqrt(2*gravity*3*dH)
             + fraction*Cd_gate_submerged*width*a*sqrt(2*gravity*dH)
```

This is not an invented minimum or an interpolation of differently defined heads. It is the source's stated combination. The rounded boundary 0.67 gives a right-hand transition limit sqrt(0.99)≈0.994987 times the pure free-flow value: approximately a **0.501% downward jump**. The code preserves and tests that source approximation instead of replacing 0.67 with 2/3 or silently smoothing it. The upper boundary tends continuously to submerged discharge. Numerical trajectory/refinement checks must account for this small lower-boundary jump. The low-depth weir and weir–gate transition still need additional geometry and a complete applicable treatment.

**Source discrepancy:** the 6.6 parent overview [G3] calls H/a≥1.25 “weir” flow, while its neighboring text places weir flow at H/a≤1. The original official CPD-69 v6.0 PDF, May 2021, p.8-7 [G4], explicitly calls the ≥1.25 branch gate flow; its pp.8-8–8-9 reproduce the gate equations. The implemented cutoff follows that original manual, rather than silently adopting the contradictory word in the web overview. The current official 7.0 overview also corroborates gate/orifice flow above that cutoff [G9]; this does not claim the helper reproduces the full 7.0 algorithm.

Closure returns zero with regime `closed`. Upstream level at/below the sill returns zero with `dry`. Equal/adverse upstream-to-wetland head returns zero with `adverse_head`, enforcing the project's existing **one-way routing assumption**, not asserting that physical gates cannot reverse flow. No storage/head/control is clipped. Low upstream depth relative to opening is explicitly unsupported as above, rather than treated as a full orifice.

The approved first case explicitly adds `H_energy >= 1.25*a` for every open, wet, forward-flow gate as a physical hydraulic-applicability constraint. It supplements the original B.1 formulation rather than silently narrowing it. A violating RHS/RK stage stops that candidate before unsupported flow is evaluated. The driver records branch, time, upstream depth, opening height and `1.25*a-H_energy`; a rejected candidate has no complete physical peak. Official OPT receives finite numerical rejection sentinels and continues. A valid trajectory that violates a safety limit is reported separately.

The targeted recheck and concrete hypothetical dimensions are recorded in [gate_method_proposal.md](06_gate_proposal.md). Its approved condition is `H_energy >= 1.25*gate_max_opening*g` whenever the opening is positive, the opening is wet, and the head is forward. This is enforced throughout integration and reported among the 20 labeled callback residuals (the original 18 plus one hydraulic constraint per gate). The low-head review found crest-dependent weir/submergence inputs and an unresolved head-versus-flow transition description; it does not justify inventing a universal closure.

## Approved parameters and evidence limits

| Item | Units / mathematical domain | First exploratory case | Evidence and limits |
|---|---|---|---|
| gravity | m/s²; >0 | 9.80665 | Stated standard-gravity convention, not a local measurement. |
| gate_width[i] | m; >0 | [1.00, 0.80] | Approved hypothetical design dimensions; not surveyed gates. |
| gate_max_opening[i] | m; >0 | [0.30, 0.30] | Approved hypothetical travel and linear actuator mapping. |
| gate_sill_z[i] | m; finite | [100.00, 99.80] | Approved hypothetical bottom-mounted sills in the assumed common datum. |
| Cd_gate_free[i] | dimensionless; >0 | [0.60, 0.60] | Midpoint of G1's typical 0.5–0.7 engineering range; not local calibration. |
| Cd_gate_submerged[i] | dimensionless; >0 | [0.80, 0.80] | G2's separately stated typical orifice value; no verified numerical range invented. |
| Approach-head approximation | explicit Boolean assumption | True, approved | Pool stage approximates upstream energy grade; non-negligible approach velocity remains unimplemented. |
| Sill/head-ratio applicability | H/a≥1.25 for positive flow | Additional approved physical constraint | G4 resolves G3 discrepancy; user accepted the restriction for this hypothetical case. |
| Submergence criteria | .67 and .80 | Fixed reference criteria | G1/G2; G9 supplies transition weights/head expression. Not tuned. |
| Remaining low-head regimes | documented laws/geometry | Excluded by approved applicability condition | No fabricated closure or continuation through rejected states. |

The numerical domains are mathematical checks, not local evidence. No Milwaukee gate dimensions or calibrated coefficients are claimed. With modeled gate tailwater equal to the associated wetland stage, the fixed safety-feasible bounds `98.90 < 100.00` and `98.80 < 99.80` place both wetlands below their respective sills. Wetland-induced submergence is excluded from every safety-feasible trajectory under this assumption; low upstream head is the relevant method issue.

USBR Chapter 9 was checked independently [G5–G8], including its full PDF equations [G10, pp.9-3–9-6]. It distinguishes contraction, velocity/friction, and approach-velocity factors. Its 0.61 coefficient concerns a specified fully contracted rectangular measuring orifice, with detailed edge/clearance/approach conditions and differential head at least 0.2 ft (0.06096 m). Its free-orifice head is referenced to the opening center, unlike the HEC sluice's sill reference. Section 8 discusses suppressed contractions and requires an approach area at least eight opening areas for its approximate adjustment. Those geometries/coefficient conventions are not established for the hypothetical project gates. Consequently, neither 0.61, its suppression correction, nor its instrument-specific minimum head is transferred into the HEC helper. The HEC submerged formula mathematically tends to zero as differential head tends to zero; that limiting test does not establish coefficient accuracy at very small head.

## Conservation and checks actually executed

Compute each new gate flow once per RHS evaluation. Subtract that identical value from its dam pool and add it to its associated wetland; replace, rather than add to, the old gate contribution. Keep the independent spillway path active even when its gate is closed. These integration requirements belong in the B.2 adapter; the helper alone cannot certify the full system's water balance or trajectory applicability.

Run `python3 -m unittest B.test_gate_hydraulics -v` from the project root. **16 helper tests and 6 B.2 wrapper tests passed (22 total)** on 2026-10-01: independent free/submerged arithmetic fixtures, distinct coefficients, opening-height mapping and absence of duplicate opening/cap, datum invariance, dimensional scaling, closure/dry/adverse routing, parameter/input checks, transition midpoint/endpoints and the documented lower-boundary jump, and explicit failure for unsupported weir/partial-opening states. Test constants—including gravity=10 m/s², width=2 m, and coefficients .6/.8—are isolated arithmetic fixtures, not approved project values.

Dimensions: width×opening has m²; sqrt(gravity×head) has m/s; discharge is m³/s. The original gate schedule was tested under the revised B.2 hydraulics and failed the head/opening constraint. A declared search of uniform open-loop factors `[0, 0.25, 0.10, 0.05]` selected the largest hydraulically admissible tested scale, 0.10, as an explicitly labeled replacement reference. It has been simulated and is numerically valid and safety-feasible; fixed limits and geometry were unchanged. Official OPT and refined verification are complete. See [09_exploratory_results.md](09_exploratory_results.md) for actual failure details, final metrics and solver termination. A.2 regression values are not calibration targets for this revised model.

## Primary references

All accessed 2026-10-01; complete relevant sections/equations were read, not just search snippets.

- **G1:** USACE HEC, *HEC-RAS Hydraulic Reference Manual*, v6.6, “Sluice Gate,” Eq.1, Eq.2 and submergence paragraph. [Official section](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.6/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways/sluice-gate).
- **G2:** USACE HEC, same manual/version, “Radial Gates,” head definitions, transition paragraph, and Eq.3 (fully submerged orifice). Its radial free-flow exponent law is not adopted. [Official section](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.6/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways/radial-gates).
- **G3:** USACE HEC, same manual/version, “Hydraulic Computations through Gated Spillways,” overview; contradictory wording noted above. [Official section](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.6/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways).
- **G4:** Gary W. Brunner / USACE HEC, *HEC-RAS River Analysis System Hydraulic Reference Manual*, CPD-69, v6.0, May 2021, pp.8-7–8-9 and 8-11–8-12; Eqs.8-3, 8-4 and low-flow discussion. [Official PDF](https://www.hec.usace.army.mil/software/hec-ras/documentation/HEC-RAS_6.0_Reference_Manual.pdf).
- **G5:** U.S. Bureau of Reclamation, *Water Measurement Manual*, third edition 1997, revised reprint 2001, Chapter 9 §1, definition, suppression, and partial-opening caveats. [Official section](https://www.usbr.gov/tsc/techreferences/mands/wmm/chap09_01.html).
- **G6:** USBR, same manual, Chapter 9 §4, conditions (6), (12)–(14): contraction clearances, differential head, free-orifice head, approach area. [Official section](https://www.usbr.gov/tsc/techreferences/mands/wmm/chap09_04.html).
- **G7:** USBR, same manual, Chapter 9 §5, Eq.9-1b and coefficient definitions. [Official section](https://www.usbr.gov/tsc/techreferences/mands/wmm/chap09_05.html).
- **G8:** USBR, same manual, Chapter 9 §8, Eq.9-2 suppression adjustment and applicability. [Official section](https://www.usbr.gov/tsc/techreferences/mands/wmm/chap09_08.html).
- **G9:** USACE HEC, *HEC-RAS Hydraulic Reference Manual*, v7.0, “Hydraulic Computations through Gated Spillways,” corroborating overview, regime definitions, and explicit transition equation/fraction used here. [Official section](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/7.0/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways).
- **G10:** USBR, *Water Measurement Manual*, third edition, revised reprint 2001, Chapter 9 (especially printed pp.9-3–9-6; PDF pp.162–165). Used to verify equations whose HTML versions are GIF images. [Official PDF](https://www.usbr.gov/tsc/techreferences/mands/wmm/WMM_3rd_2001.pdf).
