# Stage D broad-crested extension: source and compatibility review

Date: 2026-11-17. This is a read-only engineering review. No model, configuration, solver result or hydraulic law was changed.

## Finding

The newly authorized hypothetical broad-crested geometry cannot, by itself, produce a physically consistent extension while **exactly retaining the existing same-width, sill-head free-sluice law for every H/a >= 1.25**. Under a genuine rectangular broad-crested control, negligible approach velocity head and a hydrostatic critical section, the retained law exceeds the ideal available-energy discharge near that endpoint. A fitted weir coefficient or continuous interpolating curve would conceal, not resolve, that inconsistency.

This is a conditional incompatibility for the specified broad-crested, same-width, same-datum interpretation. It is not a universal theorem about short crests, ogee crests, nonhydrostatic jets, wider upstream control sections, or different conduit geometries. Those are materially different physical specifications and cannot be substituted without evidence and authorization. It also does not invalidate every historical Stage B result: the historical model and its declared assumptions remain unchanged.

## Verified primary sources

1. [USACE HEC-RAS 7.0, Hydraulic Computations through Gated Spillways](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/7.0/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways). The technical reference identifies weir flow, transition over H/a = 1–1.25, and the sluice formula using upstream energy head. It allows several crest shapes, so its generic regime labels do not establish that every independent coefficient choice is compatible with every crest. Its earlier transition description blends heads iteratively; the later overview describes a linear flow transition. An interpolation statement is not an energy-consistency proof for this chosen coefficient/geometry combination.

2. [USBR Water Measurement Manual, Chapter 2, section 13](https://www.usbr.gov/tsc/techreferences/mands/wmm/chap02_13.html). The rectangular broad-crested derivation uses the critical-flow relation and specific energy. It distinguishes total upstream energy head from measured water depth: replacing one by the other requires a velocity correction. Friction, curvature and velocity-distribution corrections also require explicit treatment. The project currently equates pool water level with energy level; an unreported velocity correction would change that assumption.

3. [USBR Water Measurement Manual, Chapter 7, section 3](https://www.usbr.gov/tsc/techreferences/mands/wmm/chap07_03.html). True broad-crested behavior is associated approximately with head/streamwise-crest-length ratios between 1/20 and 1/2. Shape and regime therefore depend on head as well as physical crest dimensions. The same wall can behave differently outside that range. This is an applicability criterion, not a calibrated coefficient for this hypothetical installation.

4. [FHWA HDS 7, 2012, page 3.27, equations 3.60–3.61](https://www.fhwa.dot.gov/engineering/hydraulics/pubs/hif12018.pdf#page=81). This archived edition explicitly derives the ideal rectangular broad-crested coefficient as approximately 1.705 in SI, with hydrostatic pressure at the critical section. It lists typical practical SI coefficients around 1.44–1.68, dependent on geometry, roughness and head. The derivation is used here as basic hydraulics, not current site-design approval; the document is marked superseded by the 2024 edition.

5. [USACE HEC-RAS 7.0, High Flow Computations, Weir Flow Computations](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/7.0/modeling-bridges/hydraulic-computations-through-the-bridge/high-flow-computations). This provides the broad-crested SI coefficient range 1.38–1.71 and explicitly uses energy elevation to include approach velocity. The bridge/roadway examples are coefficient-range analogues, not a specification or calibration of these gates. Their range does not support Cw approximately 2.126 for the selected broad-crested interpretation.

## Independent dimensional and energy check

Use b for transverse clear width, a for actual opening, H for upstream energy head above the unchanged sill, and g = 9.80665 m/s². All following numbers are derived from the existing parameters, not fitted.

For a hydrostatic, rectangular free-surface section of unit discharge q = Q/b,

    E(y) = y + q²/(2 g y²).

The minimum occurs at critical depth. With available specific energy no greater than H and no energy source,

    Q <= b sqrt(g) (2H/3)^(3/2)
      = Ccrit b H^(3/2),
    Ccrit = 1.704603788 m^(1/2)/s.

The Stage B free-sluice branch is

    Q_B = 0.60 b a sqrt(2 g H).

Writing r = H/a gives its equivalent H^(3/2) coefficient:

    C_B(r) = (0.60/r) sqrt(2g).

| Derived check | Value |
|---|---:|
| C_B at r = 1.25 | 2.125771465 m^(1/2)/s |
| C_B / Ccrit at r = 1.25 | 1.247076581 |
| Excess over ideal broad-crested capacity | 24.707658% |
| r where the retained law first meets the ideal bound | 1.558845727 |
| Largest Cd compatible with the ideal bound at r = 1.25 | 0.481125224 |
| Energy head multiplier needed to support that endpoint at the ideal bound | 1.158587262 |

Thus the conflict covers **1.25 <= H/a < 1.558845727**, not only one floating-point boundary. A real crest with positive losses would generally be more restrictive than the ideal calculation. For b = 1 m, a = 0.10 m and H = 0.125 m, the retained branch gives 0.093946714 m³/s while the ideal broad-crested capacity is 0.075333556 m³/s. Both scale with transverse width, so gate 2 has the same relative conflict.

A continuous transition ending at the retained value necessarily exceeds the ideal broad-crested bound near H/a = 1.25. This holds regardless of whether a linear or smooth interpolation is chosen. At fixed H, an ordinary same-width broad crest with Cw near 1.7 would also show an anomalous increase in modeled discharge when a wide-open gate is partly closed to a = 0.8H if its endpoint were forced to the old value. Continuity alone is therefore an insufficient test.

## Geometry and small-head qualifications

The accepted widths are 1.00 and 0.80 m; maximum openings are 0.30 m; sills are 100.00 and 99.80 m; the free-sluice coefficient is 0.60. The existing approximation neglects approach velocity head. A longer, rougher or more contracted crest cannot increase the above ideal capacity. Conversely, invoking substantial approach energy, a wider control section, an ogee/narrow crest, a lower control elevation, or nonhydrostatic suction changes the physical interpretation and requires its own coupled derivation. Those changes are not established merely by choosing a broad-crest length.

A candidate streamwise length Lc must also have an explicit head range. The USBR criterion is approximately 0.05 <= H/Lc <= 0.5. No positive fixed length satisfies that criterion for every H approaching zero. An ideal formula's limit Q -> 0 is mathematically sensible but does not validate its quantitative accuracy for arbitrarily shallow, friction-dominated flow. A defensible variant would need a declared operating range or an additional supported small-head treatment. No length is selected here because the higher-head compatibility conflict must first be resolved.

## Required decision

The exact-preservation requirement prevents implementing a physically justified broad-crested all-regime variant from these assumptions alone. Report this specific conflict before changing the model. Possible future decisions are to authorize a Stage D-only revision of the overlapping gate-law regime or to select and substantiate a materially different physical gate/crest configuration. Neither change should be silently made. Historical B files can remain byte-identical under either future option.

Do not use Cw = 2.125771465 merely to match the old endpoint, silently cap the discharge, raise the 1.25 boundary, or claim that the USACE interpolation text validates this combination. No universal physical extension, successful pilot or G5 resolution is claimed by this source review.
