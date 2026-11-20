# Basis for exploratory parameter uncertainty and sensitivity

The available evidence does not establish calibrated parameter probability distributions. The user has explicitly approved an **exploratory 64-pair Monte Carlo experiment**, using independent log-uniform multipliers on `[0.5, 2]` for `K_u`, one common multiplier on all three upland routing rates, and `R_r`. The marginal laws, support endpoints and independence are declared study assumptions, not estimated confidence limits, measured Milwaukee ranges, or hard limits on physically possible parameters.

The user also approved a separate deterministic factorial study crossing `0.5, 1, 2` for all three coordinates: 27 configurations. A small numerical pilot must precede the full Monte Carlo execution. Each configuration/draw uses a matched wetland and no-wetland pair under the same saved September 2022 event. Weather uncertainty is separate. This document records authorization and statistical meaning; it does not claim that any planned simulations have completed. Actual outcomes must come from saved execution results.

## What the existing project establishes

The [parameter mapping](parameter_mapping.md) identifies the exact laws and baseline values. The [B source documentation](../../../B/data_sources_and_assumptions.md) classifies these values as uncalibrated engineering assumptions. The [existing source catalog](../../../B/source_catalog.json) identifies SSURGO and soil-data access resources, but it does not contain an extracted installation-specific soil profile or an upscaled estimate of `K_u`. The [data-source review](../../../B/docs/03_data_sources.md) records the absence of a surveyed installation, local parameter measurements and compatible initial-state observations. USGS 04087088 remains geographic/event context for a larger catchment, not a calibration target for this hypothetical 0.66 km^2 system.

The [submitted D.1](/Users/zhangjiaxuan/Downloads/SMCR-JiaxuanZhang-HWD.1-2026.pdf) proposes sensitivity analysis where physical-parameter data are insufficient; it supplies no numerical sampling law. The current three-category selection is the user's later refinement.

## Primary-source review

Sources below were checked on November 18, 2026. They support interpretation and identify missing evidence; none establishes a local probability law for the selected model coefficients.

| Quantity | Verified source information | Applicability to this study |
|---|---|---|
| `K_u` | USACE's Green-Ampt tutorial, Table 4, lists generic saturated-conductivity estimates of 0.1 in/h for loam and 0.3 in/h for silt loam, equivalent to 2.54 and 7.62 mm/h. Its parameter-estimation note requires calibration and validation of initial estimates. [USACE tutorial](https://www.hec.usace.army.mil/confluence/hmsdocs/hmsguides/applying-loss-methods-in-hec-hms/applying-the-green-and-ampt-loss-method) | The accepted 3.6 mm/h coefficient is of a comparable numerical order, but this does not identify the installation's texture or validate its value. Our law suppresses infiltration as finite soil storage fills and is not the Green-Ampt law. Texture-table conductivity cannot be assigned directly as its calibrated distribution. |
| `K_u` | EPA SWMM Hydrology, Section 4.4.4, describes texture-based estimates from approximately 5000 US soils as preliminary and not site specific. It explicitly notes that the cited Table 4-7 does not supply conductivity ranges, despite ranges for other listed properties. [EPA manual, pp. 113-114](https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=P100NYRA.TXT) | A texture mean or another property's spread cannot furnish a variance or confidence interval for this project's `K_u`. No statistical fit is performed. |
| `K_u` | NRCS NEH Part 630, Chapter 7, bases hydrologic soil groups on the least-transmissive layer and restrictive-layer/water-table depth. One deep-profile Group C criterion uses conductivity greater than 0.40 and less than 4.0 micrometers/s, with restriction and water table deeper than 100 cm. [NRCS Chapter 7, pp. 7-1 to 7-4](https://directives.nrcs.usda.gov/sites/default/files2/1712930597/11905.pdf) | Baseline `1e-6 m/s` equals 1 micrometer/s, but no required profile/depth evidence identifies a soil group here. Classification thresholds are not probability bounds for a lumped infiltration coefficient. |
| Routing rates | The USACE Clark technical reference uses linear storage proportional to discharge and estimates response times using basin characteristics, flow lengths/velocities, regional relations and observed responses. [USACE Clark model](https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/transform/clark-unit-hydrograph-model) | This supports the time-constant interpretation of an uncapped `Q=kS` link. It does not identify this project's 8/10/20 h branch times or a distribution of their common multiplier. No branch travel-time observations or installation flow-path geometry are available. |
| Sensitivity levels | A USACE tutorial explicitly compares half, baseline and double time-of-concentration/storage-coefficient values. [USACE hydrograph sensitivity tutorial](https://www.hec.usace.army.mil/confluence/hmsdocs/hmsguides/applying-transform-methods-in-hec-hms/estimating-clark-unit-hydrograph-parameters/time-of-concentration-and-storage-coefficient-impacts-to-the-hydrograph) | This is precedent for a transparent multiplicative sensitivity experiment. It does not turn the user's levels into a local uncertainty interval. In this project the routing coordinate multiplies rates, so it divides the corresponding timescales. |
| `R_r` | USACE's outflow-curve tutorial constructs a relation using observed stage/discharge and checks simulated against observed elevation. [USACE outflow-curve tutorial](https://www.hec.usace.army.mil/confluence/hmsdocs/hmsguides/modeling-reservoirs-in-hec-hms/modeling-reservoir-operations-with-hec-hms/outflow-curve-option-for-an-hec-hms-reservoir-element) | Our `R_r` is the inverse slope of an assumed uncapped head-discharge relation, not Manning roughness or the Clark coefficient. No matching local outlet rating, geometry-based derivation, or residual fit supplies its range or probability law. No numerical coefficients from that tutorial's example dam are imported. |

## Concrete experimental settings

| Quantity | Factor 0.5 | Factor 1 | Factor 2 | Status |
|---|---:|---:|---:|---|
| `K_u` (m/s) | `5e-7` | `1e-6` | `2e-6` | User-approved sensitivity levels |
| `K_u` (mm/h) | 1.8 | 3.6 | 7.2 | Unit conversion of the same levels |
| `alpha_route` (dimensionless) | 0.5 | 1 | 2 | One confirmed common factor for `k_ud[0]`, `k_ud[1]`, `k_ur` |
| Corresponding branch timescales (h) | `[16, 20, 40]` | `[8, 10, 20]` | `[4, 5, 10]` | Derived reciprocal-rate scales; not observed parcel travel times |
| `R_r` (s/m^2) | 0.225 | 0.45 | 0.9 | User-approved sensitivity levels |
| Isolated uncapped outlet timescale `A_r*R_r` (h) | 1.25 | 2.5 | 5 | Algebraic interpretation with `A_r=20000 m^2`; not full-system spin-up time |

Keep routing/outlet caps, geometry, boundary paths, gate hydraulic law, controller policy, initial storages and wetland design fixed. All three routing rates move together; their fixed ratios are an explicit approved coupling assumption, not measured correlation. In the approved Monte Carlo model this produces correlation +1 among those rate parameters, without requiring correlation +1 among their realized fluxes. Independence between the three groups is a separately approved exploration assumption; it has not been estimated from observations.

A full factorial crosses every level without assuming that physical parameters are statistically independent. Treat its 27 points as experimental configurations, not 27 random environmental observations. Minima, maxima and contrasts over this grid are configuration summaries; the fraction of safe or beneficial configurations is not a population probability. Do not infer a median or 95th-percentile uncertainty outcome by assigning equal probability to grid points without a separately approved probability model. Keep incomplete, invalid, unsafe and negative-effect cases in the accounting and do not shrink the design after seeing their outcomes.

## Approved exploratory Monte Carlo specification

Generate 64 independent, identically distributed parameter vectors. Within each vector the three log-multipliers are mutually independent under the approved hypothetical law:

```text
For j in {infiltration, routing, outlet resistance}:
    Z_j ~ Uniform(-ln(2), +ln(2))
    M_j = exp(Z_j)

K_u  = M_infiltration * K_u,baseline
k_ud = M_routing      * k_ud,baseline
k_ur = M_routing      * k_ur,baseline
R_r  = M_resistance   * R_r,baseline
```

This specification gives positive log-uniform factors on `[0.5, 2]`, with median 1 and equal probability in reciprocal multiplicative intervals. It is not uniform on the original parameter scale, and its arithmetic mean is `(2-0.5)/ln(4)`, approximately 1.082, rather than 1. It encodes the user-approved symmetric multiplicative exploration rule, not a distribution estimated from hydrologic observations or D.1.

Independence of the three `Z_j` is an explicit simplifying scenario assumption. It is not established by the source review or by the deterministic factorial design. The hydrologic outcomes remain coupled through the ODE, and the wetland/no-wetland outputs for each vector are paired rather than independent observations. Use the same vector for both models. Retain the reproducible random seed and all realized vectors, and keep the 64 Monte Carlo vectors separate from the 27 factorial configurations and numerical-refinement repeats.

**Current status: approved for exploratory execution after the numerical pilot; no further approval is needed for this declared design.** Monte Carlo summaries describe propagation of the specified hypothetical parameter law conditional on one saved storm, model, policy, design and initial state. They are not calibrated confidence intervals for local physical parameters, flood return periods or proof of universal wetland benefit. Preserve all failed, invalid, unsafe and negative-effect outcomes; any summary limited to complete valid pairs must state that conditional denominator and exclusions.
