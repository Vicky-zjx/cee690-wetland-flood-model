# Historical-storm robustness analysis

This directory extends the completed Stage D.2 parameter experiment to separately selected historical Milwaukee-area storms. It is an analysis workspace and technical findings record, **not the final HW D.2 submission report**.

The unchanged eight-state wetland and six-state no-wetland equations, Stage D gate hydraulics, 120 s feedback controller, actuator limits, physical initial states, wetland design and safety thresholds are reused. No physical parameter is tuned to a storm or to its outcomes. The original 64 parameter vectors are read from the previous immutable study plan, without new sampling. The same vector is used in both models and matched across events.

The original initial storages and other boundary conditions are assumed for every event. They are not reconstructions of each historical storm's antecedent soil, groundwater or river conditions. The wetland model retains its documented extra 5,500 m³ initial wetland water. Humidity, pressure, downstream/groundwater head and external-flow assumptions remain unchanged. These are historical-weather forcing experiments for a hypothetical installation, not simulations of measured Milwaukee infrastructure.

In inherited JSON metadata, `fixed_configuration.schedule_description` retains the old loader's historical Stage B fixed-schedule wording. The operative schedule is identified by `schedule_kind`, `controller_policy`, the 120 s clock and the saved actual controller histories: the validated Stage D feedback policy. The historical description does not select the executed control law.

## Files

- `event_selection.py`, `event_selection.md`, `events/`: rainfall-only selection rule, complete candidate/exclusion inventory, source provenance and saved joint weather rows.
- `adapter.py`: a private event provider around the existing frozen parameter/controller runner. No original function, cache or file is replaced.
- `study_plan.json`: selected events, original parameter-plan hash, actual unchanged vectors and predeclared numerical pilot.
- `run_study.py`: original-reference regression, event pilots, and complete paired ensembles.
- `results/`: per-event outcomes, failures/safety checks, controller histories and saved trajectories.
- `analyze_results.py`, `tables/`, `plots/`: saved-result analysis; no ODE execution or new fitting.
- `mechanism_check.py`, `mechanism_check.md`: one preselected draw's full-window peaks, wetland storage and later recession, calculated from saved trajectories.
- `FINDINGS.md`: completed scientific findings and limitations; all five event ensembles have 64 valid, safety-feasible pairs.
- `preservation_before.json`, `preservation_check.json`: hash checks protecting all prior work and supplied sources.
- `independent_review.py`, `independent_review.json`: read-only audit of retained outputs, matching, safety/applicability, mass balance and benchmark preservation.

## Execution

Use the existing Python environment and pinned local `multivarious` checkout described in the [previous parameter study](../parameter_uncertainty/README.md). The official solver revision remains `06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02`.

```bash
cd /Users/zhangjiaxuan/Documents/Codex/Modeling
/opt/miniconda3/bin/python3 -B D/d2_analysis/historical_storm_robustness/event_selection.py
/opt/miniconda3/bin/python3 -B D/d2_analysis/historical_storm_robustness/run_study.py --phase prepare
/opt/miniconda3/bin/python3 -B D/d2_analysis/historical_storm_robustness/run_study.py --phase regression
/opt/miniconda3/bin/python3 -B D/d2_analysis/historical_storm_robustness/run_study.py --phase pilot --workers 4
/opt/miniconda3/bin/python3 -B D/d2_analysis/historical_storm_robustness/run_study.py --phase ensemble --workers 4
/opt/miniconda3/bin/python3 -B D/d2_analysis/historical_storm_robustness/analyze_results.py
/opt/miniconda3/bin/python3 -B D/d2_analysis/historical_storm_robustness/mechanism_check.py
/opt/miniconda3/bin/python3 -B D/d2_analysis/historical_storm_robustness/preservation.py verify --output D/d2_analysis/historical_storm_robustness/preservation_check.json
```

Simulation phases refuse to overwrite existing result directories. Preserve an existing phase by archiving it within this analysis directory before repeating it. Keep the immutable plan and original sources. The analysis command can regenerate derived tables/figures without rerunning models. No package installation, data download, Git commit or push is required by this workflow.

## Numerical and data conventions

The original September 2022 benchmark retains its original prepared forcing and previously saved 64-pair results; it is counted once. A fresh paired regression checks complete arrays against a saved original trajectory. New storms use the unchanged Stage C hourly rows: rain and radiation are held over each actual hour; temperature and wind use the existing interpolation between hourly centers. The prior conservative hourly adapter is used without modification. Rain and solar integrals are checked within control blocks before joining output arrays; integrating a concatenated discontinuous rainfall curve would introduce a different numerical convention.

For each new event, three **existing** vectors represent low, central and high mean log multiplier—not low/median/high flood severity. Both models are checked at 15 s and 7.5 s with the same 120 s controller clock. All full-horizon safety and applicability results remain distinct from numerical agreement. The complete ensemble runs only after the pilot passes. Valid unsafe trajectories and negative or zero wetland effects are retained. An invalid trajectory has a recorded cause and unavailable physical full-horizon metrics, not a fabricated peak or penalty.

Full arrays are retained for the numerical pilot and original-regression check. New ensemble trajectories retain time, all water states, water levels, downstream discharge, actual openings and balance histories; their unchanged event forcing is saved separately. Reused benchmark metrics and diagnostics point back to the original results; missing benchmark final-storage arrays are not invented.

The saved external-net storage histories recover the unchanged runner's per-control-block balance integrals. They expose that accounting for inspection; they are not an independent second integration-based conservation test.

## Interpretation

Positive paired differences mean **no-wetland minus wetland**. Full-horizon peaks retain time zero. Maximum rise and duration above the initial river stage are separate reference-response metrics. There is no accepted downstream flood threshold, so no flood-exceedance probability, flood duration or flood-excess volume is inferred.

The 64 vectors retain their approved independent exploratory log-uniform multipliers on [0.5, 2], with one shared multiplier for all three routing rates. They are not calibrated local distributions. Each event has its own distribution conditional on this assumed parameter law. The same vectors induce dependence across events; comparisons use matched sample IDs. Rainfall-based event selection supplies no occurrence probabilities, and event-by-parameter rows are never pooled as independent weather draws.

Official `multivarious` utilities provide the professor's rank/(N+1) positions and Beta order-statistic limits, retained in the saved diagnostic tables. Initial-peak plateaus produce exact zero-effect atoms, so the figure uses a documented statistical adaptation: the solid curve is the exact ties-aware count/N empirical CDF, the dotted curve is the official plotting-position output, and shading gives 95% Clopper–Pearson binomial intervals at predeclared fixed thresholds. This handles atoms without pretending they have a continuous distribution. The limits concern finite parameter sampling at each fixed threshold, not simultaneous coverage or empirical weather risk. See [cdf_interval_method.md](cdf_interval_method.md). Sample outcome percentiles are not confidence limits. The assumed outlet rating links stage and discharge reductions, so parameter associations are reported for both and are not interpreted as causal sensitivity indices.
