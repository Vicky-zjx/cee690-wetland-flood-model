# Historical-storm robustness: technical findings

This is a focused analysis record, not the final HW D.2 submission. The saved 64 parameter vectors were applied to five separate, fixed 144-hour weather windows. The original September 2022 results were reused once; four additional event ensembles were completed. All **320 pairs** are complete, hydraulically valid and safe under the existing model checks. There are only **64 distinct parameter vectors** and five deliberately selected weather windows, not 320 independent weather observations.

The observed result is **nonnegative full-horizon peak attenuation across the tested samples, with substantial event dependence and many zero effects**. It does not establish strictly positive benefits under every storm or improvement in every response metric.

## Weather selection and scope

The [selection rule and source audit](event_selection.md) rank observed 24-hour rainfall cores, require complete joint precipitation/temperature/solar/wind coverage and positive temperature, and select nonoverlapping windows without using simulation outcomes. The benchmark is reserved once. The unchanged new-event rows retain NOAA station observations, their original interpolation/trace flags, and NASA POWER solar grid estimates. No rainfall is scaled, no missing hour is made dry, and no new parameter samples are drawn.

| Window start (UTC) | Total rain, 144 h (mm) | Maximum UTC-hour depth (mm) | Maximum 24 h (mm) |
|---|---:|---:|---:|
| 2022-09-09 00:00, benchmark | 145.300 | 15.906 | 141.400 |
| 2023-08-13 08:00 | 58.900 | 16.693 | 53.300 |
| 2022-07-04 21:00 | 66.427 | 26.147 | 47.200 |
| 2019-08-25 02:00 | 42.800 | 11.367 | 42.300 |
| 2024-08-04 15:00 | 41.200 | 20.020 | 41.200 |

These are heavy-rainfall windows, potentially containing multiple wet spells. All additional windows have less total and 24-hour rain than the benchmark; some have larger hourly bursts. This is not evidence of robustness to storms with greater 24-hour depth than September 2022. Heavier incomplete and freezing candidates remain in the [inventory](events/candidate_inventory.csv.gz), with explicit exclusion reasons. Completeness and nonfreezing screens introduce selection bias.

## Peak effects and initial storage

All differences below are **no-wetland minus wetland**. The full horizon includes time zero. Counts are sample descriptions, not empirical flood probabilities.

| Window | Mean stage reduction (mm) | Median (mm) | Range (mm) | Positive / zero / negative | Mean peak-discharge reduction (m³/s) | Initial peak: wet / no wet |
|---|---:|---:|---:|---|---:|---|
| September 2022 | 91.573 | 94.055 | 34.284–159.218 | 64 / 0 / 0 | 0.176506 | 0 / 0 |
| August 2023 | 91.045 | 86.859 | 0–213.837 | 55 / 9 / 0 | 0.158558 | 38 / 9 |
| July 2022 | 73.814 | 69.598 | 0–176.461 | 53 / 11 / 0 | 0.123844 | 39 / 11 |
| August 2019 | 26.178 | 0 | 0–121.525 | 27 / 37 / 0 | 0.035132 | 54 / 37 |
| August 2024 | 41.003 | 7.341 | 0–164.744 | 34 / 30 / 0 | 0.058315 | 52 / 30 |

Every zero peak difference in the added events has both maxima at the initial river stage, **97.35 m**. Their later hydrographs can differ even though their full-horizon maxima coincide. No initial peak was removed and no antecedent storage was adjusted. Maximum rise above the initial stage is saved separately; because both cases share the initial stage and include time zero, its paired difference equals the peak-stage difference. It does not eliminate the initial-stage floor.

Rain intensity alone does not order attenuation: July 2022 has the highest hourly burst here but a smaller mean reduction than August 2023. August 2023's mean is close to the September benchmark, but its matched parameter effects are larger in 29 draws and smaller in 35. The corresponding larger-effect counts are 19 for July 2022, one for August 2019 and eight for August 2024. Rainfall timing/duration, other joint weather variables and initial-storage effects differ together; five selected windows do not isolate a causal rainfall-intensity relationship.

The [illustrative hydrographs](plots/05_fixed_draw_hydrographs.png) use the first saved draw, `mc_00000`, selected in advance. They show attenuation and later recession, including visible event responses below the initial peak. The [storage/timing check](mechanism_check.md) keeps full-window peaks primary and labels any separate storm-response timing explicitly. This one draw does not represent all parameter settings. Original initial states remain hypothetical assumptions for every event, not observed antecedent reconstructions or spin-up results.

## Other outcomes and mechanism limits

Peak attenuation does not imply shorter duration or less total outflow. The wetland case spends **longer above the initial stage** in 27, 8, 4, 1 and 1 draws for September 2022, August 2023, July 2022, August 2019 and August 2024, respectively. This reference level is not a flood threshold.

Mean cumulative outlet-volume differences over 144 hours are **+1,490.67, −523.71, −668.12, −908.95 and −741.40 m³**, in the same order. Every new-event draw has more cumulative outflow from the wetland model. The wetland case starts with an accepted extra **5,500 m³** of wetland water and has different final storage and evaporation pathways. Consequently, these volume differences are not pure storm-water retention or permanent flood-volume removal. [Storage accounting](tables/water_storage_balance.csv) retains initial/final total storage, integrated external net flux and residuals without adding internal transfers twice. Final-storage fields unavailable in the original benchmark summaries are explicitly left unavailable, not reconstructed.

Both models use the same controller, and their realized gate openings agree to approximately **5.73×10⁻¹⁵ m** across the retained pairs, with no latch disagreements. Paired effects therefore do not arise from a changed gate policy in this experiment. The comparison includes the existing wetland compartments, routing and accepted initial water; it does not separately identify an effect of adding initially empty wetlands.

## Parameter associations and the outlet equation

The assumed outlet law remains `Q = min(3, max(H_r − 97.15, 0) / R_r)` in SI units. All retained peak outlets are on its uncapped branch. The identity **ΔH_peak = R_r × ΔQ_peak** holds to within **2.86×10⁻¹⁴ m**. Stage and discharge are therefore not independent confirmations of benefit, and a strong stage association with resistance includes the rating conversion itself. Resistance also changes drainage dynamics, so its effect cannot be dismissed entirely as a plotting artifact.

The [rank-association table](tables/parameter_rank_associations.csv) shows an event-dependent pattern. For the benchmark, Spearman associations of infiltration, routing and resistance with peak-discharge reduction are **0.743, −0.542 and 0.038**. Across the four added events they span **−0.186 to −0.077**, **0.190 to 0.631**, and **0.470 to 0.820**, respectively. Resistance's stage association is 0.970 for the benchmark and 0.828–0.865 for new events. These are descriptive marginal associations within 64 matched exploratory samples, affected by ties and interactions; they are not causal effects, variance-based sensitivity indices or fitted physical laws.

## Numerical, hydraulic and statistical checks

The September 2022 regression at saved draw `mc_00000` reproduces all eight saved trajectory arrays exactly for both models. Before ensemble execution, three existing low/central/high mean-log-multiplier vectors were tested for each new event at **15 s and 7.5 s**, keeping the controller's **120 s decision clock**. All 24 model-specific resolution comparisons passed: maximum common-time stage difference **6.65×10⁻⁹ m**, maximum native sampled stage-peak difference **5.30×10⁻⁸ m**, and zero latch mismatches. This verifies the pilot points, not numerical convergence at every possible parameter setting.

Across the five retained ensembles there are no hydraulic failures, numerical failures or safety violations. The smallest dam/wetland safety margin is **0.153490 m** in the benchmark and **0.401331 m** in the new events. The largest recorded relative mass residual is **1.601×10⁻⁷**, below the unchanged 10⁻⁴ tolerance. New-event gate opening/head checks remain inside the documented formulation, with minimum open-gate `H/a` **1.98904**. Only closed and energy-gate regimes occur; these runs do not provide additional validation of weir/transition flow. Safety limits and applicability criteria were not relaxed. Full per-branch and per-trajectory evidence is retained in [safety checks](tables/safety_checks.csv) and [diagnostics](tables/trajectory_diagnostics.csv).

[ECDF figures](plots/02_per_event_peak_distributions.png) are separate for each fixed event. The professor's official rank/(N+1) plotting positions are retained. Because zero effects create probability atoms under the exploratory parameter law, exact count/N curves and fixed-threshold **95% Clopper–Pearson intervals** are used for the displayed uncertainty bands; see [method and verification](cdf_interval_method.md). These are pointwise finite-sample intervals, not simultaneous bands, calibrated physical uncertainty or weather-occurrence probabilities. No event-by-parameter pooling or outcome filtering is used.

## Remaining limits

The experiment supports conditional, event-dependent peak attenuation under the unchanged hypothetical geometry, controls and exploratory parameter law. It does not establish universal wetland benefit, an observed Milwaukee infrastructure response, a calibrated distribution, or a flood-exceedance probability. No justified downstream flood threshold is available.

The benchmark retains its original trace-zero lower bounds and 120 s rainfall/solar ramps, while new windows retain the existing C hourly representation; this processing difference is documented rather than silently normalized. Regional station/grid inputs are spatial transfers to the hypothetical installation. Fixed initial conditions, extra wetland initial water, incomplete-data screening, rainfall temporal resolution, an assumed outlet rating and only five weather windows limit generalization. No colder, missing-data or more extreme synthetic storm was substituted to force an outcome.

All A/B/C and prior D code/results are protected by the [preservation check](preservation_check.json). No final HW D.2 report, Git commit or push is part of this work. [README](README.md) provides reproduction commands; [figure captions](figure_captions.md) and [complete outcome tables](tables/outcome_summaries.csv) provide the retained evidence.
