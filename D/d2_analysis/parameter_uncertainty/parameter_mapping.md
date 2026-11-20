# Parameter uncertainty: scope and mapping

The user selected three categories for the next Stage D analysis: **soil infiltration capacity, runoff routing speed, and river drainage conditions**, and confirmed their mapping to `K_u`, one shared routing-rate multiplier, and `R_r`. The study holds the September 9-15, 2022 UTC event forcing fixed. Weather uncertainty remains separate. The user has approved 64 paired iid Monte Carlo vectors with independent log-uniform group multipliers on `[0.5, 2]`, following a small numerical pilot, and a separate 27-configuration factorial sensitivity study at factors `0.5, 1, 2`. These probability laws and design levels are exploratory assumptions, not empirically calibrated bounds. See [uncertainty_basis.md](uncertainty_basis.md).

## Submitted D.1 and the current decision

The submitted source is [SMCR-JiaxuanZhang-HWD.1-2026.pdf](/Users/zhangjiaxuan/Downloads/SMCR-JiaxuanZhang-HWD.1-2026.pdf), a one-page HW D.1 submission dated November 12, 2026. Its physical-parameter discussion includes:

> I can address parameters for which data are insufficient through sensitivity analysis.

D.1 identifies weather as the principal uncertainty and mentions wetland drainage resistance and seepage as examples of uncertain physical properties. It does **not** specify parameter names, numerical ranges, probability distributions, covariance/correlation values, or a sampling law. The user's later selection of the three categories above governs the present scope; the earlier wetland examples do not replace that selection. A numerical distribution cannot be attributed to the submitted document.

## Existing quantities and confirmed coordinates

All baseline values below are retained engineering assumptions or values derived from assumed timescales, not locally measured/calibrated parameter estimates. See the [B evidence register in the model documentation](../../../B/data_sources_and_assumptions.md), particularly its soil/runoff and outlet sections, and the [A parameter constructor](../../../A/wetland_flood_model.py).

| Selected category | Existing quantity and active baseline | Confirmed mapping | Approved exploratory probability specification |
|---|---|---|---|
| Soil infiltration capacity | `K_u = 1e-6 m/s = 3.6 mm/h` | Multiply the effective upland-to-soil infiltration-rate coefficient by the selected factor. This is not a measured saturated hydraulic conductivity. | Log-uniform factor on `[0.5, 2]`, independent of the routing and resistance factors by assumption. |
| Runoff routing speed | `k_ud = [1/(8 h), 1/(10 h)] = [3.472222222222222e-5, 2.777777777777778e-5] s^-1`; `k_ur = 1/(20 h) = 1.388888888888889e-5 s^-1` | Multiply all three rates by the same positive `alpha_route`. This confirmed coupling preserves their ratios. | `alpha_route` is a new uncertainty coordinate: log-uniform on `[0.5, 2]`, independent of the other two group factors by assumption. It is not a D.1-specified distribution. |
| River drainage conditions | `R_r = 0.45 s/m^2` | Multiply final river-outlet resistance by the selected factor. Increasing it reduces uncapped outlet discharge at a fixed effective head. Other drainage quantities remain fixed. | Log-uniform factor on `[0.5, 2]`, independent of the infiltration and routing factors by assumption. |

The baseline is `alpha_route = 1`. The approved three-level experiment crosses `0.5, 1, 2` for each of the three coordinates, giving 27 configurations. This deterministic grid is separate from the 64 iid Monte Carlo vectors. Probability support, laws and independence come from the user's explicit Monte Carlo approval, not from the grid or empirical evidence.

## Equations and distinctions

Write `[z]_+ = max(z, 0)` and `phi_u = S_u/(S_u + V_u)`. The retained infiltration law is:

```text
Q_inf = A_u * K_u * phi_u * [1 - S_m/C_m]_+                 [m^3/s]
```

`K_u` is a rate scale, whereas `C_m = 120000 m^3` is soil water-storage capacity. `V_u = 600 m^3` controls surface-water availability; `S_fc = 60000 m^3` and `k_mg = 1/(5 days)` govern soil-to-groundwater recharge. Those are separate parameters and are not automatically varied by the selected infiltration category. Actual infiltration depends on surface supply and soil saturation, so it does not equal `K_u` alone.

The routing laws are:

```text
Q_ud,i = min(Q_udcap,i, k_ud,i * S_u)                       [m^3/s]
Q_ur   = min(Q_urcap,   k_ur   * S_u)                       [m^3/s]

Confirmed common coordinate:
k_ud,i = alpha_route * k_ud,i,baseline
k_ur   = alpha_route * k_ur,baseline
```

Here `Q_udcap = [1.2, 1.0] m^3/s` and `Q_urcap = 0.7 m^3/s` are distinct capacity assumptions and remain fixed. A larger rate multiplier means shorter uncapped routing timescales, `tau = tau_baseline/alpha_route`; it is not the same as multiplying travel times by that factor. Caps can flatten the discharge response to increasing rates. The shared factor preserves all three parameter ratios and makes those three rates perfectly positively correlated under the approved Monte Carlo law. That is an explicit user-approved coupling assumption, not measured correlation. It does not imply perfectly correlated realized flows after caps and dynamics act. Cross-group independence is also an approved exploration assumption; empirical dependence remains unknown.

The river outlet law is:

```text
Q_out = (1 - B) * min(Q_rcap,
                     [H_r - max(h_d, z_ro)]_+ / R_r)       [m^3/s]
```

Other possible interpretations of drainage conditions are `Q_rcap = 3 m^3/s` (unblocked outlet cap), `z_ro = 97.1 m` (outlet sill), `h_d = 97.15 m` (prescribed downstream head), and `B = 0` (prescribed blockage fraction). The confirmed resistance-only study leaves these unchanged. `h_d` and `B` are boundary inputs, not interchangeable with `R_r`; their values come from the [approved boundary configuration](../../../B/no_spillway_config.json). No downstream flood threshold has been established.

## Application to the paired models

The [eight-state equations](../../../B/no_spillway_model.py) and [six-state equations](../../../B/wetland_counterfactual.py) use the same infiltration, upland routing and final river-outlet laws. The [Stage D hydraulic wrapper](../../gate_hydraulics_revision/model.py) replaces gate fluxes only; it does not replace these laws. The accepted B configuration currently overrides only humidity and atmospheric pressure, so the baseline values above remain active in D.

Any approved parameter realization must be shared by the wetland and no-wetland cases, with the same saved event weather, shared initial storages, boundary paths, wetland design and controller policy. The wetland case retains its documented additional 5500 m^3 initial wetland water. A common feedback policy need not imply identical realized gate histories after parameters change; actual histories must still be compared. No parameter draw should trigger reoptimization, relaxed safety limits, state clipping or exclusion merely because a wetland effect is negative.

If `R_r` varies, initial `Q_out` also varies at the unchanged initial river stage. Metrics relative to initial discharge must therefore use each trajectory's actual `Q_out(0)`, rather than the baseline value of approximately `0.444444 m^3/s`. The complete event must retain its initial sample and report safety and model-validity outcomes separately.

**Current status:** mappings, shared routing coupling, the 27 factorial settings, and the 64-pair independent log-uniform Monte Carlo design are explicitly approved. A small numerical pilot precedes full execution. The assumptions remain uncalibrated, but no further approval is needed for this design. Existing A/B/C files, D hydraulic equations and saved results remain unchanged by this mapping document.
