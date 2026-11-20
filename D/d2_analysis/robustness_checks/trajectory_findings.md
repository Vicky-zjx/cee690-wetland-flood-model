# Saved-trajectory robustness checks

This analysis reads existing saved results only. It does not integrate an ODE, redraw parameters, alter weather, remove time zero, or modify any historical model, result or controller. The 320 paired outcomes use the same 64 exploratory parameter vectors under five selected fixed weather windows; they are not 320 independent weather observations.

## Coverage and primary result

All 320 primary peak values, peak times, peak discharges and classifications were checked against their saved per-model JSON records. Full arrays are retained for all 256 added-event pairs and benchmark draws `mc_00000` and `mc_00001`: 258 paired waveform audits. The other 62 benchmark pairs remain included in primary statistics but have no retained full arrays for this check. No unavailable waveform is reconstructed.

| Event | Positive / zero / negative ΔH | Both initial maxima | Zero pairs with different waveforms | Array pairs |
|---|---:|---:|---:|---:|
| Sep 2022 benchmark | 64 / 0 / 0 | 0 | 0 / 0 | 2 |
| Aug 2023 | 55 / 9 / 0 | 9 | 9 / 9 | 64 |
| Jul 2022 | 53 / 11 / 0 | 11 | 11 / 11 | 64 |
| Aug 2019 | 27 / 37 / 0 | 37 | 37 / 37 | 64 |
| Aug 2024 | 34 / 30 / 0 | 30 | 30 / 30 | 64 |

Every exact zero in the saved primary metric is equivalent, within these results, to both first full-window maxima occurring at the shared initial stage 97.35 m. All 87 zero pairs have full arrays and later stage differences larger than 1e−9 m; none has an identical hydrograph. Their maximum absolute waveform differences range from 49.722 to 102.834 mm. These remain valid zero full-window peak-reduction outcomes. Their later dynamical differences are masked by the initial maximum in that metric; secondary differences are not a recovered or replacement primary benefit. This does not mean the initial stored wetland water has been causally isolated or that a different primary objective should replace the accepted one.

[Primary peak summaries](tables/trajectory_peak_summary.csv) report mean, median, minimum and maximum stage peak, discharge peak and their first-occurrence times for each event and model, using all 64 saved samples. The ten rows retain time-zero maxima and give their counts; mean or median peak time across such samples is a descriptive full-window statistic, not an isolated storm-delay estimate.

Peak-discharge difference is also exactly zero in these same pairs. With the unchanged monotone outlet rating, equal maximum stages imply equal peak discharges; these are not independent confirmations that the hydrographs coincide. For all array-available pairs, the saved stage and discharge maxima and their first-occurrence times match the source records.

Realized gate openings match between models to 5.73e-15 m in all saved paired history summaries, with 0 latch disagreements. Direct array comparisons cover the 258 retained pairs and agree with those summaries within 1e−12 m. Every full-array pair has distinct river waveforms; the difference is not a changed gate-operation policy. The wetland model's additional accepted initial water and different storage/routing remain bundled in the comparison, not separately identified causal effects.

## Separate descriptive event-response window

For timing diagnostics only, a secondary window begins at the existing manifest's largest quantified wet-spell start and ends at 144 h. The spell was selected from rainfall, without model outcomes; positive quantified hours are joined unless separated by at least six quantified-zero hours, and trace-zero is not verified dry. Secondary peaks are recorded separately and do not replace primary full-window peaks. The primary peak delay is left unavailable when either full-window peak occurs at t=0.

| Event | Zero-primary secondary difference range (mm) | Secondary peak-delay range (h), all available pairs | Later wet > no pairs |
|---|---:|---:|---:|
| Sep 2022 benchmark | No zero draws | 1.983–2.196 | 2 / 2 |
| Aug 2023 | 52.570–79.747 | 1.767–4.013 | 64 / 64 |
| Jul 2022 | 44.154–72.582 | 2.875–3.842 | 64 / 64 |
| Aug 2019 | 31.663–71.877 | 2.325–5.000 | 64 / 64 |
| Aug 2024 | 40.850–83.450 | 0.433–4.008 | 64 / 64 |

The later comparison begins after both secondary peaks. A later positive difference means wetland stage exceeds no-wetland stage by more than 1e−6 m at a saved node; its first time is resolved to 15 s, not an interpolated crossing. This is a descriptive reporting tolerance, not a safety or applicability relaxation. Smaller primary peaks can coexist with a higher later recession. Some later rain may still be present, so the check does not imply that rainfall has ceased.

Secondary maxima occur at the secondary left boundary in 0 wetland and 0 no-wetland array trajectories. Any such time is retained and flagged, not interpreted as an interior storm-response peak.

## Figures and deterministic example selection

[Figure 1](plots/01_zero_peak_diagnostics.png) retains all 320 full-window outcomes and displays exact zero counts; its waveform panels use all 258 retained paired arrays. No zero-effect samples are removed. The secondary comparison is visibly labeled as a different descriptive statistic.

[Figure 2](plots/02_event_hydrographs.png) uses the lowest saved sample ID with an exact zero primary difference and the lowest ID with a positive difference in each event. These IDs were written to `tables/trajectory_example_selection.json` before reading arrays for the audit or generating the plots. This rule does not select the largest benefit. The benchmark has no zero pair and its panel remains empty. Both columns retain the entire 144-hour window and can use different parameter vectors, so their contrast is not an intervention on a single parameter. Within each event row the two plotted panels share a stage scale; different event rows may use different stage scales.

| Event | First zero ID | First positive ID |
|---|---|---|
| Sep 2022 benchmark | None | mc_00000 |
| Aug 2023 | mc_00002 | mc_00000 |
| Jul 2022 | mc_00002 | mc_00000 |
| Aug 2019 | mc_00000 | mc_00001 |
| Aug 2024 | mc_00000 | mc_00001 |

## Rainfall and limits

Rainfall totals, derived UTC-hour and 24-hour maxima, wet-hour counts and wet-spell duration are copied from the unchanged event manifest into `tables/trajectory_event_summary.csv`. All added events have lower 24-hour depth than the September benchmark. Four weather variables vary jointly; the benchmark uses its original 120 s reconstruction whereas the added events use preserved Stage C hourly inputs. These selected, complete nonfreezing windows do not identify a causal rainfall-intensity relationship or represent the distribution of all regional storms.

This audit preserves the initial-condition floor instead of claiming it has been removed. It does not determine how an initially empty wetland, an observed antecedent state, a different event start, or a separate warm-up would behave. Such claims require controlled simulations and are outside this saved-result check. No downstream flood threshold is supplied or invented.

## Reproduction and output evidence

Run `/opt/miniconda3/bin/python3 -B D/d2_analysis/robustness_checks/trajectory_checks.py` from the project root (the tested existing environment). The system `python3` lacks Matplotlib; no new dependency was installed. Only `trajectory_*` tables, this note and the two assigned figures are written under the new robustness-check directory. `tables/trajectory_sample_audit.csv` contains all 320 rows, source paths and hashes, complete coverage flags, primary/secondary timing, waveform differences and gate checks. `tables/trajectory_validation.json` records the aggregate audit and input hashes.
