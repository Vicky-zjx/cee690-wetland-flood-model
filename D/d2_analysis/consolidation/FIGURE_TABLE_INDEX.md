# D.2 figure and table index

This index links the existing results; it introduces no simulations or new statistical population. All 30 PNGs were visually reviewed. The [review record](figure_review.md) and [complete inventory with SHA-256 hashes](figure_inventory.csv) document that check.

| Study | Count and interpretation |
|---|---|
| Original stochastic weather | 20 original IDs: 16 warm complete pairs (10 Copula/VAR, 6 bootstrap), 4 cold exclusions. Deterministic controls and numerical repeats are separate. |
| Fixed-weather parameters | 27 factorial combinations and a separate 64-draw iid exploratory Monte Carlo study, using the September 2022 event. |
| Historical storms | Five fixed weather windows × the same 64 parameter vectors = 320 paired outcomes. The original benchmark 64 are reused once; only 256 additional event pairs were run. |
| Saved-result robustness | No new ODE runs: all 320 paired summaries/endpoints, 258 retained paired waveforms, and the same 64 parameter vectors for subset/reweighting diagnostics. |
| Antecedent conditions | 12 new pairs from two events × two prefixes × three predeclared vectors; six saved no-prefix pairs are reused. Refinement repeats are not additional draws. |

A pair contains an eight-state wetland simulation and a six-state no-wetland simulation. These study counts must not be added into one iid ensemble. Peak differences use no wetland minus wetland. Unless labeled otherwise, maxima cover the entire 144 h event including time zero. No accepted downstream flood threshold has been supplied; initial-reference exceedances are not flood probabilities.

## Original stochastic weather

| Figure | Related tables |
|---|---|
| [Initial-peak diagnosis](../plots/01_initial_peak_diagnosis.png) | [metrics](../results/metrics.csv); [weather_summary](../results/weather_summary.csv) |
| [Full-window peak distributions](../plots/02_primary_peak_ecdfs.png) | [ecdf_values](../results/ecdf_values.csv); [paired_metrics](../results/paired_metrics.csv) |
| [Initial-reference response distributions](../plots/03_initial_reference_ecdfs.png) | [ecdf_values](../results/ecdf_values.csv); [paired_metrics](../results/paired_metrics.csv) |
| [Numerical resolution and controller checks](../plots/04_resolution_and_controller_checks.png) | [convergence](../results/convergence.csv); [controller_audit](../results/controller_audit.csv); [safety_audit](../results/safety_audit.csv) |

[Captions / study explanation](../figure_captions.md).

## Fixed-weather parameter study

| Figure | Related tables |
|---|---|
| [27-point factorial peak comparison](../parameter_uncertainty/plots/factorial_peak_reduction.png) | [sensitivity_paired_metrics](../parameter_uncertainty/tables/sensitivity_paired_metrics.csv) |
| [Parameter/output scatter](../parameter_uncertainty/plots/mc_parameter_effects.png) | [monte_carlo_paired_metrics](../parameter_uncertainty/tables/monte_carlo_paired_metrics.csv) |
| [Official paired peak-stage distribution](../parameter_uncertainty/plots/official_mc/official_ecdf_delta_peak_stage_m.png) | [official_ecdf_delta_peak_stage_m](../parameter_uncertainty/plots/official_mc/official_ecdf_delta_peak_stage_m.csv) |
| [Official paired peak-discharge distribution](../parameter_uncertainty/plots/official_mc/official_ecdf_delta_peak_discharge_m3_s.png) | [official_ecdf_delta_peak_discharge_m3_s](../parameter_uncertainty/plots/official_mc/official_ecdf_delta_peak_discharge_m3_s.csv) |
| [Official wetland peak-stage distribution](../parameter_uncertainty/plots/official_mc/official_ecdf_wetland_peak_stage_m.png) | [official_ecdf_wetland_peak_stage_m](../parameter_uncertainty/plots/official_mc/official_ecdf_wetland_peak_stage_m.csv) |
| [Official no-wetland peak-stage distribution](../parameter_uncertainty/plots/official_mc/official_ecdf_no_wetland_peak_stage_m.png) | [official_ecdf_no_wetland_peak_stage_m](../parameter_uncertainty/plots/official_mc/official_ecdf_no_wetland_peak_stage_m.csv) |
| [Official wetland peak-discharge distribution](../parameter_uncertainty/plots/official_mc/official_ecdf_wetland_peak_discharge_m3_s.png) | [official_ecdf_wetland_peak_discharge_m3_s](../parameter_uncertainty/plots/official_mc/official_ecdf_wetland_peak_discharge_m3_s.csv) |
| [Official no-wetland peak-discharge distribution](../parameter_uncertainty/plots/official_mc/official_ecdf_no_wetland_peak_discharge_m3_s.png) | [official_ecdf_no_wetland_peak_discharge_m3_s](../parameter_uncertainty/plots/official_mc/official_ecdf_no_wetland_peak_discharge_m3_s.csv) |
| [Official parameter/output scatter matrix](../parameter_uncertainty/plots/official_mc/official_parameter_output_scatter.png) | [monte_carlo_paired_metrics](../parameter_uncertainty/tables/monte_carlo_paired_metrics.csv) |

[Captions / study explanation](../parameter_uncertainty/README.md).
The six official marginal plots use rank/(N+1), with pointwise order-statistic limits. See [official-figure captions](../parameter_uncertainty/plots/official_mc/figure_captions.md). Factorial points have no probability weights.

## Historical storms

| Figure | Related tables |
|---|---|
| [Event rainfall summaries](../historical_storm_robustness/plots/01_event_rainfall_summary.png) | [selected_events](../historical_storm_robustness/events/selected_events.csv); [event_weather_summary](../historical_storm_robustness/tables/event_weather_summary.csv) |
| [Per-event peak distributions](../historical_storm_robustness/plots/02_per_event_peak_distributions.png) | [all_event_paired_metrics](../historical_storm_robustness/tables/all_event_paired_metrics.csv); [fixed_grid_pointwise_cdf_intervals.csv](../historical_storm_robustness/tables/fixed_grid_pointwise_cdf_intervals.csv.gz) |
| [Matched event contrasts](../historical_storm_robustness/plots/03_matched_event_peak_effects.png) | [matched_event_changes](../historical_storm_robustness/tables/matched_event_changes.csv) |
| [Per-event parameter rank associations](../historical_storm_robustness/plots/04_parameter_rank_associations.png) | [parameter_rank_associations](../historical_storm_robustness/tables/parameter_rank_associations.csv); [outlet_rating_identity](../historical_storm_robustness/tables/outlet_rating_identity.csv) |
| [One fixed draw across all events](../historical_storm_robustness/plots/05_fixed_draw_hydrographs.png) | [illustrative_hydrograph_sources](../historical_storm_robustness/tables/illustrative_hydrograph_sources.csv); [illustrative_mechanism](../historical_storm_robustness/tables/illustrative_mechanism.csv) |

[Captions / study explanation](../historical_storm_robustness/figure_captions.md).
Use [event selection and provenance](../historical_storm_robustness/event_selection.md) and [CDF interval method](../historical_storm_robustness/cdf_interval_method.md) when interpreting cross-event distributions. The benchmark retains its original processing; the other four use the documented Stage C hourly mapping.

## Saved-result robustness

| Figure | Related tables |
|---|---|
| [Zero peaks versus nonidentical trajectories](../robustness_checks/plots/01_zero_peak_diagnostics.png) | [trajectory_sample_audit](../robustness_checks/tables/trajectory_sample_audit.csv); [trajectory_event_summary](../robustness_checks/tables/trajectory_event_summary.csv) |
| [Predeclared zero and positive examples](../robustness_checks/plots/02_event_hydrographs.png) | [trajectory_example_selection](../robustness_checks/tables/trajectory_example_selection.json); [trajectory_peak_summary](../robustness_checks/tables/trajectory_peak_summary.csv) |
| [Initial/final storage and outlet accounting](../robustness_checks/plots/03_storage_balance.png) | [storage_paired_decomposition](../robustness_checks/tables/storage_paired_decomposition.csv); [storage_event_summary](../robustness_checks/tables/storage_event_summary.csv); [storage_wetland_timing](../robustness_checks/tables/storage_wetland_timing.csv) |
| [Exact zero atoms and pointwise intervals](../robustness_checks/plots/04_ecdf_zero_atoms.png) | [statistical_zero_atoms](../robustness_checks/tables/statistical_zero_atoms.csv); [statistical_saved_grid_checks](../robustness_checks/tables/statistical_saved_grid_checks.csv) |
| [Existing-sample subset diagnostics](../robustness_checks/plots/05_sample_size_stability.png) | [statistical_nested_prefixes](../robustness_checks/tables/statistical_nested_prefixes.csv); [statistical_subset_summaries](../robustness_checks/tables/statistical_subset_summaries.csv) |
| [Parameter-law reweighting diagnostics](../robustness_checks/plots/06_parameter_law_checks.png) | [statistical_law_outcomes](../robustness_checks/tables/statistical_law_outcomes.csv); [statistical_law_diagnostics](../robustness_checks/tables/statistical_law_diagnostics.csv) |

[Captions / study explanation](../robustness_checks/figure_captions.md).
The [event/model peak and timing summary](../robustness_checks/tables/trajectory_peak_summary.csv) contains all 64 summaries per event and model. There are 62 benchmark pairs without retained full waveforms; all benchmark endpoint summaries remain available.

## Antecedent conditioning

| Figure | Related tables |
|---|---|
| [August 2019 continuous hydrographs](../antecedent_conditions/plots/01_historical_20190825T0200Z_hydrographs.png) | [model_window_metrics](../antecedent_conditions/tables/model_window_metrics.csv); [paired_window_metrics](../antecedent_conditions/tables/paired_window_metrics.csv) |
| [August 2024 continuous hydrographs](../antecedent_conditions/plots/01_historical_20240804T1500Z_hydrographs.png) | [model_window_metrics](../antecedent_conditions/tables/model_window_metrics.csv); [paired_window_metrics](../antecedent_conditions/tables/paired_window_metrics.csv) |
| [Original-event-window outcomes](../antecedent_conditions/plots/02_paired_event_outcomes.png) | [paired_window_metrics](../antecedent_conditions/tables/paired_window_metrics.csv); [change_from_saved_baseline](../antecedent_conditions/tables/change_from_saved_baseline.csv) |
| [Actual feedback gate openings](../antecedent_conditions/plots/03_gate_histories.png) | [controller_summary](../antecedent_conditions/tables/controller_summary.csv) |
| [Whole-system storage contrasts](../antecedent_conditions/plots/04_storage_contrasts.png) | [compartment_storages](../antecedent_conditions/tables/compartment_storages.csv); [paired_water_accounting](../antecedent_conditions/tables/paired_water_accounting.csv) |
| [Joint saved meteorological inputs](../antecedent_conditions/plots/05_weather_inputs.png) | [event_conditioning_summary](../antecedent_conditions/tables/event_conditioning_summary.csv) |

[Captions / study explanation](../antecedent_conditions/figure_captions.md).
Read [execution accounting](../antecedent_conditions/tables/execution_accounting.csv), [event conditioning summary](../antecedent_conditions/tables/event_conditioning_summary.csv) and [forcing provenance](../antecedent_conditions/forcing/README.md). Event-window statistics use evolved states at time zero; full 168 h / 216 h statistics retain the earlier assumed initial state.
