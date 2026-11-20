# Continuous antecedent-condition experiment

This separate Stage D.2 experiment follows the predeclared design in [the earlier robustness findings](../robustness_checks/FINDINGS.md). It is an interpretation and numerical check, not the final HW D.2 submission report. Earlier A/B/C and validated D code and results are protected by a new file-hash manifest.

## Fixed design

- Two selected historical events: August 2019 and August 2024; each original event lasts 144 hours.
- Prefixes: the actual preceding 24 or 72 hours of cached joint meteorological data.
- Three existing parameter vectors: `mc_00016`, `mc_00024`, `mc_00018`, selected previously as low, central and high mean log multiplier, not by benefit or event severity.
- Twelve physical pairs, 24 primary model executions at 15 s. The longest 72 h August 2024 case with `mc_00018` is also checked at 7.5 s (two additional model executions of the same physical pair).
- Six original no-prefix pairs are read from saved results. None of the previous 320 pairs is rerun. No new parameter vectors or meteorological scenarios are generated.

Every pair uses the same saved soil infiltration capacity, common runoff-routing multiplier and river outlet resistance. The original feedback policy, 120 s controller period, actuator limits, hydraulic laws, geometry, wetland design, boundary inputs and safety limits remain unchanged. Completion, hydraulic validity and safety are separate outcomes.

## Continuous initialization and forcing

At the **beginning of each prefix**, the accepted eight-state initial vector is `[600, 90000, 3000, 2500, 3000, 2400, 409600.00000000035, 7000]` m³, ordered as upland, soil, Wetland 1, Wetland 2, Dam 1, Dam 2, groundwater and river. The six-state model retains the corresponding shared entries. The wetland system therefore begins with 5,500 m³ more water; none is deleted, redistributed or subtracted as an output correction.

Both gates initially have actual opening 0.0015 m and the original enabled latches. One uninterrupted frozen feedback runner advances through all 168 or 216 hours. There is no controller or water-state reinitialization at hour 24 or 72. A private, opt-in adapter extends the numerical duration; it does not edit the historical 144-hour loader or physical equations. All water states, gate positions, hysteresis memory and the continuous controller clock cross the event boundary.

The prefix and original event are mapped separately with the unchanged hourly input mapper so the original event's solver weather values remain bit-for-bit unchanged. This preserves its original first-half-hour temperature/wind endpoint holds. The explicit weather representation can jump at the seam; the physical states and actuators remain continuous. Precipitation depth and solar energy are integrated within each block before concatenation. See [weather provenance, actual native intervals and seam checks](forcing/README.md). Station-derived rainfall/temperature/wind and gridded solar estimates are not all direct observations at the hypothetical installation.

Different prefix lengths start the same assumed initial state at different dates and expose it to different weather and drainage durations. The experiment is not a reconstruction of observed antecedent storage, an equilibrium proof, or an isolated intervention on the extra 5,500 m³.

## Metrics and accounting

[Model-window metrics](tables/model_window_metrics.csv) retain three windows explicitly: the full continuous prefix-plus-event horizon, the antecedent interval, and the unchanged original 144-hour event. The saved no-prefix reference has identical full/event windows. All windows include their initial endpoint and retain any peak there. Absolute stage, rise above each window's own initial stage, discharge, first peak time, volume and storage are reported separately. No flood threshold or flood occurrence probability is inferred.

[Paired metrics](tables/paired_window_metrics.csv) use **no wetland minus wetland**. Positive stage/discharge differences indicate lower wetland-case peaks. Positive volume differences mean lower wetland-case cumulative outlet volume. After conditioning, the event-start storage difference is measured from the carried states; it is not assumed to remain −5,500 m³.

For each window, `S_end − S_start = N − V + r`, where `N` is net non-outlet external input, `V` is river outlet volume, and `r` is the recorded integration residual. Thus `ΔV = ΔS_start − ΔS_end + ΔN + Δr`. Internal gate/wetland transfers are never counted as new external water. The aggregate `N` reconstructed from the saved external-net residual is an accounting check, not a separate independent flux estimate or source-water tracer. New runs additionally retain direct per-block/phase integrals in JSON diagnostics.

[Importance-weight verification](statistics/README.md) reuses all 64 original samples. Weight concentration and outcome-contribution concentration are distinguished; this does not turn the three conditioning configurations into a probability ensemble.

## Reproduction

Use the project's existing Python environment and official solver checkout, as described in [the Stage D dependency guide](../../reproducibility_inputs/README.md). Commands below run from the project root. The observed weather cache and original saved pilot arrays are required; no data are downloaded. Completed experiment destinations are protected from accidental overwrite. To reproduce from scratch, use a separate copy with an empty **new experiment** directory while retaining its source scripts and all prior inputs; do not erase the accepted results.

```sh
/opt/miniconda3/bin/python3 -B D/d2_analysis/antecedent_conditions/preservation.py snapshot
/opt/miniconda3/bin/python3 -B D/d2_analysis/antecedent_conditions/forcing/prepare.py
/opt/miniconda3/bin/python3 -B D/d2_analysis/antecedent_conditions/run_study.py --phase prepare
/opt/miniconda3/bin/python3 -B D/d2_analysis/antecedent_conditions/run_study.py --phase pilot
/opt/miniconda3/bin/python3 -B D/d2_analysis/antecedent_conditions/run_study.py --phase remaining --workers 3
/opt/miniconda3/bin/python3 -B D/d2_analysis/antecedent_conditions/analyze_results.py
/opt/miniconda3/bin/python3 -B D/d2_analysis/antecedent_conditions/forcing/plot_weather.py
/opt/miniconda3/bin/python3 -B D/d2_analysis/antecedent_conditions/statistics/verify_importance_weights.py
/opt/miniconda3/bin/python3 -B D/d2_analysis/antecedent_conditions/verify_results.py
/opt/miniconda3/bin/python3 -B D/d2_analysis/antecedent_conditions/preservation.py verify
```

For an existing completed directory, run only the saved-data analysis/audits. See [technical findings](FINDINGS.md), [figure captions](figure_captions.md), the frozen [study plan](study_plan.json), and [numerical pilot results](results/pilot_summary.json). Every original unsafe or failed attempt must remain in the accounting; no feasibility filter is applied to manufacture positive benefits.
