# Stage D gate geometry and hydraulic method

## Scope and status

This opt-in Stage D rating changes the gate equations and their applicability regions. It is a hypothetical engineering model, not a calibrated Milwaukee structure, an exact HEC-RAS implementation, or an experimentally validated transition law. A/B/C equations and results, gate schedules, physical initial states, forcing, other model parameters and safety limits are preserved. The gate remains the only dam-pool release route; no independent dam spillway is added.

This document specifies the method and algebraic checks. Actual simulation, safety and numerical-resolution outcomes are reported separately; none is implied by selecting the rating. The rating remains a candidate, not a validated G5 correction or evidence of readiness for full HW D.2 uncertainty propagation.

## Geometry and head convention

| Property | Gate 1 | Gate 2 | Basis |
|---|---:|---:|---|
| Rectangular clear width, b | 1.00 m | 0.80 m | Original geometry |
| Gate sill elevation | 100.00 m | 99.80 m | Original hypothetical datum |
| Maximum opening | 0.30 m | 0.30 m | Original geometry; a = 0.30 g(t) |
| Streamwise flat sill/apron length, Lc | 1.00 m | 1.00 m | Explicit hypothetical dimension, retained from the prior audit |
| Effective transverse width | 1.00 m | 0.80 m | Parallel sidewalls; no end-contraction deduction |

The proposed passage has a vertical-lift gate at the upstream edge of a flat rectangular sill/apron, with a square upstream edge and an open downstream channel. The apron discharges toward the existing wetland or river. There is no ogee suction, hidden lower control elevation, wider parallel opening, added energy source or new storage compartment. This positioning is a Stage D assumption; it does not retrospectively specify the historical B structure.

Use H = upstream pool water-surface elevation minus sill elevation as total upstream energy head, retaining the negligible pool approach-velocity approximation. The [USACE inline-structure manual](https://www.hec.usace.army.mil/confluence/rasdocs/rasum/7.0/entering-and-editing-geometric-data/inline-structures-dams-weirs-and-gated-spillways) distinguishes sill reference for discharge into a channel from opening-center reference for a free atmospheric jet. Changing the reference alone would not justify a different physical outfall. The revised gate equation instead retains the downstream contracted depth explicitly.

## Why the original gate branch changes

For a uniform, hydrostatic rectangular section, specific energy is

    E = y + q²/(2 g y²),  q = Q/b.

The ideal free critical-flow coefficient is

    Ccrit = sqrt(g) (2/3)^(3/2) = 1.704603788 m^(1/2)/s,
    Qcritical = Ccrit b H^(3/2),  g = 9.80665 m/s².

This follows the [USBR rectangular broad-crest derivation](https://www.usbr.gov/tsc/techreferences/mands/wmm/chap02_13.html). The old Stage B equation, Qold = 0.60 b a sqrt(2gH), implies Ceff = 2.125771465 m^(1/2)/s at H/a = 1.25, **24.707658% above** this ideal capacity. The conflict extends over 1.25 <= H/a < 1.558845727. It is conditional on the proposed same-width broad-crest geometry and energy convention, not a universal statement about all sluice geometries. The [previous feasibility audit](../low_head_feasibility/README.md) remains intact.

## Weir coefficient and explicit loss hypothesis

For the clear-gate branch, use Qw = Cw b H^(3/2). The [USACE crest-coefficient table](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/7.0/modeling-gated-spillways-weirs-and-drop-structures/uncontrolled-overflow-weirs) supplies the broad-crested customary range 2.6–3.1. Values 2.6, 2.8 and 3.0 are declared before simulation, with 2.8 nominal. They are engineering analogues, not local measurements or fitted optimal values. The dimensional conversion is Cw_SI = Cw_US sqrt(0.3048).

Assume uniform velocity and a nonnegative head loss xi V²/(2g). At an actual critical section V² = g yc, so

    H = yc + (1 + xi) yc/2,
    yc = 2H/(3 + xi),
    Cw = sqrt(g) [2/(3 + xi)]^(3/2),
    xi = 3[(Ccrit/Cw)^(2/3) - 1].

| Cw customary | Cw SI, m^(1/2)/s | Derived xi | Role |
|---:|---:|---:|---|
| 2.6 | 1.435426069 | 0.364209256 | Lower coefficient sensitivity |
| 2.8 | 1.545843459 | 0.202038473 | Nominal assumption |
| 3.0 | 1.656260849 | 0.058095639 | Upper coefficient sensitivity |

Using this same xi in the gate branch is an explicit shared-loss hypothesis for the hypothetical passage. It is not a measured gate loss or a USACE/USBR instruction to equate the losses in different regimes. The alternative mapping xi = (Ccrit/Cw)² - 1 does not follow the actual critical-flow derivation above and is not used. Source coefficient ranges do not independently validate this loss hypothesis.

## Supported gate energy relation

[Clemmens, Wahl and Belaud (2012), USBR PAP-1094](https://www.usbr.gov/tsc/techreferences/hydraulics_lab/pubs/PAP/PAP-1094.pdf), equations 1–2, provides the energy-based gate relation and relative-opening contraction polynomial. For s = a/H,

    delta(s) = 0.617 - 0.06s + 0.026s² + 0.026s³,
    yj = delta(s) a,
    Qg = b yj sqrt[2g(H - yj)/(1 + xi)].

Here delta is dimensionless contraction, not the old discharge coefficient 0.60; yj is contracted jet depth. The denominator represents the declared loss/velocity convention. The publication restricts free-orifice flow to s < 2/3; this implementation uses the continuous limiting value at s = 2/3. It does not evaluate the polynomial at larger relative openings. The polynomial and its 0.617 limiting value were checked in the original PDF, rather than inferred from text extraction. The paper's separate Reynolds-dependent loss calibration is not claimed to have been implemented.

The frozen B parameter object still records the historical free/submerged discharge coefficients 0.60/0.80. They remain provenance values and are not active coefficients in this D rating. D uses the contraction polynomial and declared loss hypothesis for free gate flow; positive submerged flow is explicitly unsupported. Thus preserving the B configuration does not mean retaining its discharge law in D.

[Wahl (2004), USBR PAP-0937](https://www.usbr.gov/tsc/techreferences/hydraulics_lab/pubs/PAP/PAP-0937.pdf), page 4, also gives the contracted-depth free-gate energy form and treats contraction and losses as separate inputs. Keeping H - yj explicitly avoids the old near-transition omission of jet-depth energy.

## Declared transition approximation

Define r = H/a and r_end = 1.5. Evaluate the supported gate formula only at its endpoint to obtain a dimensional coefficient:

    s_end = 2/3,
    delta_end = delta(s_end),
    C_end = delta_end s_end
            sqrt[2g(1 - delta_end s_end)/(1 + xi)].

The complete positive-forward-flow rating is

    r <= 1:       Q = Cw b H^(3/2)
    1 < r < 1.5:  w = (r - 1)/0.5
                   F(r) = (1 - w) Cw + w C_end
                   Q = b H^(3/2) F(r)
    r >= 1.5:     Q = Qg.

This interpolates the discharge coefficient after factoring out b H^(3/2); it does not extrapolate the gate law into its unsupported relative-opening domain. The [USACE gated-flow overview](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/7.0/modeling-gated-spillways-weirs-and-drop-structures/hydraulic-computations-through-gated-spillways) describes linear transition between weir and gate behavior. The present normalized endpoint interpolation and 1–1.5 band are **new engineering assumptions**, not that manual's exact algorithm or original 1–1.25 limits.

The upstream-edge interpretation allows the plate to affect entry from the pool when H exceeds a, while the downstream crest controls unimpeded flow. Local drawdown, gate-lip pressure and the detailed transition free surface are unresolved. This is a physically constrained, quasi-steady reduced-order rating; its continuity and energy bounds do not establish experimental validity or uniqueness.

## Physical checks and rejected alternative

C_end/Cw is approximately 0.8150, 0.8063 and 0.8021 for the declared sensitivity set. Therefore 0 < C_end < Cw < Ccrit: transition discharge is nonnegative, increases with opening at fixed H, and stays below the ideal crest capacity. At fixed opening,

    dQ/dH = b sqrt(H) [1.5 F(r) + r F'(r)].

For the linear band above this is positive when C_end/Cw >= 2/3; all three declared cases satisfy that bound. The source gate branch must also pass head/opening monotonicity checks with the variable contraction polynomial, positive energy loss and supercritical contracted flow. Branch joins are continuous in Q; derivative continuity is not claimed.

An earlier algebraic candidate blended the current weir and gate discharges over the narrow band from r = (3 + xi)/2 to r = 2. It was rejected because discharge could decrease as upstream head increased at fixed opening, despite passing opening-monotonicity and capacity checks. It is not used in either ODE. No coefficient, crest length or operating schedule was changed to obtain simulation completion.

## Applicability and remaining limitations

- The [USBR crest classification](https://www.usbr.gov/tsc/techreferences/mands/wmm/chap07_03.html) associates broad-crested flow approximately with 0.05 <= H/Lc <= 0.5. For the fixed 1 m crest, the weir and transition branches require H = 0.05–0.50 m. The approximate source classification is adopted as a conservative, explicit modeling-domain bound. An out-of-range positive-flow state raises an applicability error; this does not mean that real flow ceases or becomes physically impossible at the boundary. A mathematical dry limit does not validate shallow, friction-dominated ratings.
- Positive forward flow requires tailwater below the sill. Tailwater at or above the sill is unsupported, not assigned an invented submergence correction. Existing closed, dry and equal/adverse-head one-way behavior returns zero. Reverse-flow physics is not added.
- Wetland safety bounds below the sills and the six-state river geometry limit exclude drowned gates within those admissible states. This does not remove the runtime tailwater check or relax any safety limit.
- The rating is instantaneous and quasi-steady. It does not resolve gate jets, hydraulic jumps, hysteresis, local momentum transients, sediment/debris, seal leakage or a Reynolds-dependent gate calibration. The original schedule includes openings as small as 1.5 mm; coefficient transfer at that scale remains an unvalidated approximation.
- Identical revised hydraulics must be used in the wetland and no-wetland cases. Full-horizon safety, mass balance, scenario completion and numerical consistency require separate actual runs. No discharge cap, opening clipping, storage clipping or increased safety limit is part of this method.
