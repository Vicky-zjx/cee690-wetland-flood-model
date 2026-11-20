# Antecedent conditions: technical findings

The predeclared **12 paired physical experiments** are complete: August 2019/August 2024 × 24 h/72 h prefixes × `mc_00016`, `mc_00024`, `mc_00018`. All 24 primary trajectories completed, were hydraulically applicable and satisfied the unchanged dam/wetland safety checks. One of these pairs was repeated at 7.5 s for numerical refinement; both additional trajectories also passed. Six original no-prefix pairs were read from saved results, not rerun. There were no new parameter draws, weather scenarios, fitted models or optimization runs. This is a bounded interpretation check, **not the final HW D.2 report**.

## Main interpretation: conditioning matters, but does not erase the full-horizon initial peak

All 12 conditioned **original 144 h event windows** show positive wetland peak attenuation, with no event-window river-stage maximum at its initial point. The river has evolved before the event: its initial level is approximately **97.151–97.161 m**, rather than the original assumed 97.35 m. This evolution was produced by the unchanged equations, real cached antecedent forcing and the unchanged feedback policy; no state was removed or reset at the event boundary.

However, **four of the 12 complete prefix-plus-event horizons still have zero peak attenuation**: `mc_00016` in both events under both prefixes. Both models retain their maximum at the beginning of the entire continuous integration. The earlier 87 zero outcomes were not rerun or reclassified. This small experiment shows how their initial-peak mechanism can persist in a full-horizon statistic while a physically evolved event window shows a different response.

Peak reductions below are **no wetland minus wetland**, in millimetres. Every maximum includes its window's initial endpoint.

| Event | Saved configuration | Original 144 h reference | Event window after 24 h | Event window after 72 h | Full 168 h horizon | Full 216 h horizon |
|---|---|---:|---:|---:|---:|---:|
| August 2019 | mc_00016 | 0.000 | 31.889 | 31.905 | 0.000 | 0.000 |
| August 2019 | mc_00024 | 30.471 | 81.305 | 79.391 | 24.287 | 14.932 |
| August 2019 | mc_00018 | 104.340 | 102.715 | 98.183 | 92.319 | 73.289 |
| August 2024 | mc_00016 | 0.000 | 40.968 | 40.404 | 0.000 | 0.000 |
| August 2024 | mc_00024 | 64.330 | 104.405 | 101.899 | 58.069 | 50.945 |
| August 2024 | mc_00018 | 135.998 | 132.917 | 128.551 | 132.917 | 119.816 |

Conditioning therefore does **not** uniformly increase estimated mitigation. For the largest predeclared mean-log configuration, event-window attenuation falls relative to the original reference. Full-horizon reductions can fall while event-window reductions rise in another configuration. Longer conditioning also does not systematically increase benefit or establish equilibrium.

The antecedent interval alone has four small negative peak contrasts, all for `mc_00018`: approximately **−0.552 to −0.872 mm**. These are retained normal response outcomes, not invalid simulations. No claim is made that wetlands reduce stage at every instant or in every window. See [all window metrics](tables/paired_window_metrics.csv), [compact peak table](tables/peak_reduction_comparison.csv), [changes from the saved reference](tables/change_from_saved_baseline.csv), and the [hydrograph captions](figure_captions.md).

## Weather and continuous state/control checks

The 2019 prefixes contain 0/0 mm of quantified station rain over 24/72 h; the 2024 prefixes contain 0/18.5 mm. All four joint inputs are complete and above the non-freezing temperature limit. The 2024 72 h input retains six trace-affected UTC hours under the existing trace-zero lower-bound convention. The original event's rainfall, temperature, solar and wind records and solver inputs remain unchanged.

The hourly mapper is applied separately to the prefix and original event, preserving the event's existing half-hour endpoint holds. Temperature/wind representation jumps at the event seam are explicit (2019: −1.5 °C/−1.65 m/s; 2024: +1.27631579 °C/0 m/s), not invented instantaneous observations. They occur on the unchanged hourly/120 s control boundary. Every state is passed continuously through that boundary; no artificial reset accompanies the input seam. Rainfall and solar integrals pass both whole-run and separate-phase conservation checks. [Native intervals, source provenance and input tests](forcing/README.md) retain station versus grid-estimate distinctions.

The independent audit replayed **302,400 controller-branch updates** against the saved accepted states. Event-boundary state differences are exactly zero; actual opening and pre-update latch memories carry forward. Maximum gate replay discrepancy is **2.78×10⁻¹⁷ m**. Within every wetland/no-wetland pair, gates agree within **5.70×10⁻¹⁵ m**, with no latch mismatches. Thus the paired contrast is not created by unequal gate policies or unequal realized gate operations. Across prefix durations, the unchanged state-dependent policy can produce different histories; these are part of the conditioning experiment, not a new control design.

These runs encounter only closed and supported energy-gate operation; the minimum positive-opening H/a is **1.9893**. They do not newly validate unvisited hydraulic regimes. The minimum recorded dam/wetland safety margin is **0.5921 m** across primary and refinement outputs/RHS checks. [Execution accounting](tables/execution_accounting.csv), [successful-RHS safety](tables/safety_at_successful_rhs.csv), [node safety](tables/safety_at_output_nodes.csv), [controller summary](tables/controller_summary.csv), and full per-case diagnostics preserve every planned outcome.

## Initial water, outlet volume and storage

At each prefix start, the wetland system contains 518,100 m³ and the no-wetland system 512,600 m³. All shared initial states are identical; the additional 5,500 m³ remains in the accepted wetland compartments. That difference is not deleted, redistributed, or subtracted from an output as a correction.

By event start, combined wetland-compartment storage is **567.0–604.2 m³ after 24 h**, and **0.0037–0.6441 m³ after 72 h** across these runs. This is the current wetland stock, not tagged original water. The total system still contains more water in the wetland case, much of it in groundwater. After 72 h, mean extra groundwater storage is about **1,943 m³ in 2019** and **2,975 m³ in 2024** across the three configurations.

The following are simple means over the **three selected configurations**, not probability-weighted population estimates. Positive outlet differences mean less outlet water in the wetland case.

| Event / prefix | Extra whole-system water at event start, wet−no (m³) | Extra final water, wet−no (m³) | Antecedent outlet difference, no−wet (m³) | Event outlet difference, no−wet (m³) | Complete-horizon outlet difference, no−wet (m³) |
|---|---:|---:|---:|---:|---:|
| August 2019 / 24 h | 2,289.44 | 3,832.64 | −3,107.09 | +2,202.05 | −905.03 |
| August 2019 / 72 h | 1,955.72 | 3,545.71 | −3,331.74 | +2,189.68 | −1,142.06 |
| August 2024 / 24 h | 2,325.92 | 3,872.63 | −3,043.51 | +2,261.30 | −782.21 |
| August 2024 / 72 h | 2,995.08 | 4,252.74 | −2,145.03 | +1,993.41 | −151.62 |

All 12 event windows have lower cumulative outlet volume in the wetland model, but **all 12 entire continuous runs have higher cumulative outlet volume in that model**. The extra antecedent drainage and subsequent event retention explain why the signs can differ. Lower event peaks or volumes are not evidence that the wetland system released less water over every accounting horizon.

For each phase, with Δ = no wetland minus wetland,

`ΔV_out = ΔS_start − ΔS_end + ΔN_nonoutlet + Δr`.

The event phase uses its actual evolved initial-storage contrast, not a fixed −5,500 m³. N includes precipitation, external inflow, groundwater-boundary exchange and evapotranspiration; internal transfers are not additional external water. The [paired accounting table](tables/paired_water_accounting.csv) and [compartment stocks](tables/compartment_storages.csv) provide all terms. The reconstructed identity is arithmetic consistency; independent numerical evidence comes from the direct block/phase flux integration and derivative identity checks. Water is not tagged by source, so this experiment still cannot isolate the causal contribution of exactly the original 5,500 m³ from all other antecedent-state changes.

## Numerical evidence and limits

The predeclared longest 72 h August 2024 / `mc_00018` pair was checked at 15 s and 7.5 s under the same 120 s controller clock. Both models pass the inherited tolerances: maximum common-grid stage difference **1.27×10⁻⁹ m**, maximum native peak difference across all reported stages **4.19×10⁻⁸ m**, event river-peak difference **1.47×10⁻¹⁰ m**, and no latch mismatch. No safety limit or numerical tolerance was relaxed. This supports the checked pair; it does not establish convergence for every parameter value or an equilibrium spin-up.

Across all 26 executions, maximum separate-phase relative mass residual is **3.39×10⁻⁷**, below the unchanged **10⁻⁴** criterion. The saved derivative mass identity error is at most **1.78×10⁻¹⁵ m³/s**. Phase-net accounting agrees within **1.13×10⁻¹⁰ m³**. The observer introduced no additional physical RHS evaluations; it only recorded the existing official solver calls. See [numerical pilot](results/pilot_summary.json) and [independent audit](independent_audit.json).

The observed centimetre-scale event-response differences are much larger than the checked timestep effects. Nevertheless, 24/72 h are chosen conditioning durations, not proven equilibration times. Every prefix begins from the same assumed states at a different date. The 2024 longer prefix also introduces actual additional rain, so its effect cannot be assigned to elapsed time alone. The three configurations do not estimate a probability distribution or replace the original 64-vector ensemble. Peak discharge is reported alongside stage; the existing outlet rating still ties their differences and is not independently field-validated. Own-window rise/duration uses each case's evolved starting level and is not a common flood-threshold comparison.

## Importance-reweighting verification

Using all 64 saved Monte Carlo vectors and no new simulations, the independent log-uniform-to-uniform reweighting gives **ESS = 43.9935/64**. The largest normalized weight is **4.98%** (`mc_00018`); the largest five/ten weights hold **20.89%/35.65%** of total weight. No single weight dominates the total probability mass.

Outcome concentration is stronger. Ranking by contribution to the uniform-weighted peak-reduction mean, the largest ten contributors supply **73.9% for August 2019** and **65.9% for August 2024**; the largest five supply **49.7%/43.5%**. Thus moderate weight ESS does not make the weak-event mean estimates insensitive to a limited subset. These outcome-ranked contributors are not necessarily the top-weight-ranked samples. The diagnostic is not a replacement sample size for confidence intervals, a calibration, or an effective number of independent historical storms. Self-normalized estimates retain finite-sample error and dependence on the exploratory input law. [Weight audit and per-sample contributions](statistics/README.md) retain all five events and signed outcomes.

## What is supported, and what remains conditional

For the selected configurations, wetland routing reduces the original-event peak after hydrological conditioning, with matched actual gates and preserved forcing. Initial-condition dominance and the choice of accounting window materially affect apparent benefit, including the sign of cumulative outlet-volume differences. The experiment does not establish universal wetland benefit, observed antecedent storage, field safety certification, flood probabilities or full Monte Carlo convergence. No further simulations are needed to state these bounded findings; expanding them to a conditioned 64-vector distribution would be a separate decision.

All **5,103** protected prior project/source files remain unchanged. A/B/C, earlier Stage D diagnostics, the original 320 pairs and the earlier robustness checks were preserved. No commit, push or final submission report was produced.
