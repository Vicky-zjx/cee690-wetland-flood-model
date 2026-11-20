# Saved-result interpretation and robustness checks

This directory strengthens the interpretation of the completed Stage D.2 studies. It is not the final HW D.2 report. It reads the original 27-case sensitivity study, 64 saved parameter vectors and five fixed historical-event ensembles. No new physical parameter vector, weather series, initial condition, ODE solution or control policy is generated here. Existing initial peaks and all outcome classifications are retained.

The primary contrast remains **no-wetland minus wetland over the entire 144-hour horizon, including time zero**. Any separate rainfall-response timing or waveform comparison is a secondary diagnostic and never replaces that contrast. Event windows are selected historical cases, not random observations from a known storm-frequency distribution.

## Analyses

| Script | Purpose | Outputs |
|---|---|---|
| `trajectory_checks.py` | Distinguish initial-peak floors from identical river responses; compare rainfall, stage/discharge peaks and timing | `tables/trajectory_*`, figures 01–02, `trajectory_findings.md` |
| `storage_checks.py` | Recover recorded final states and audit paired storage/outlet/external-flux balances; distinguish accounting from causal initial-water attribution | `tables/storage_*`, figure 03, `storage_findings.md` |
| `statistical_checks.py` | Verify zero atoms and CDF intervals; inspect finite-sample stability, exploratory law reweighting and existing factorial contrasts | `tables/statistical_*`, figures 04–06, `statistical_findings.md` |
| `antecedent_check.py` | Check existing weather coverage before proposing a conditioning experiment; no forcing or simulation is generated | `tables/antecedent_*` |
| `preservation.py` | Check that all prior project files remain byte-identical to this task's starting snapshot | `preservation_before.json`, `preservation_check.json` |

`FINDINGS.md` combines the scientific conclusions and remaining gaps. `analysis_scope.json` records the starting Git revision, key source hashes and the no-simulation scope. Each analysis records its source files and diagnostics. Scripts write only within this new directory; rerunning them replaces only their own derived outputs.

## Data availability and interpretation

All 320 pairs have saved metrics and per-case JSON diagnostics. All four additional events have complete saved compact trajectories. Only selected benchmark Monte Carlo cases have full saved trajectories; a benchmark summary is not treated as a missing-data-free waveform. Recorded final states and completed-block balance diagnostics can support endpoint accounting even when a full NPZ history was not retained. Availability is reported explicitly in the tables.

Wetland initial storage exceeds the no-wetland initial total by 5,500 m³. This exact initial-stock difference enters the mass-balance identity. It does not identify how much of downstream outflow consists of that initial water, nor the counterfactual response with different initial storage. Those questions require a justified additional experiment or a separately designed water-origin analysis. No such experiment is executed here.

The original distributions are independent exploratory log-uniform multipliers on [0.5, 2] for infiltration capacity, common runoff routing speed and river resistance. Reweighting evaluates explicitly hypothetical alternative laws on this same support using existing outputs only. It is a finite-sample approximation, not a new calibrated distribution or a new ensemble. Effective sample size diagnoses weight concentration; it is not a replacement binomial sample count. Subset diagnostics describe sensitivity to the existing finite sample and do not prove convergence.

## Reproduction

Use the existing Stage D Python environment. The required packages and external solver setup are documented in [the publication guide](../../reproducibility_inputs/README.md). These checks use the full locally retained results; the reduced GitHub snapshot intentionally excludes most dense NPZ trajectories and cannot perform all waveform checks without those arrays. No data download is performed.

From the project root:

```sh
/opt/miniconda3/bin/python3 -B D/d2_analysis/robustness_checks/trajectory_checks.py
/opt/miniconda3/bin/python3 -B D/d2_analysis/robustness_checks/storage_checks.py
/opt/miniconda3/bin/python3 -B D/d2_analysis/robustness_checks/statistical_checks.py
/opt/miniconda3/bin/python3 -B D/d2_analysis/robustness_checks/antecedent_check.py
/opt/miniconda3/bin/python3 -B D/d2_analysis/robustness_checks/preservation.py verify --output D/d2_analysis/robustness_checks/preservation_check.json
```

These commands use the verified local Miniconda interpreter. The system `python3` lacks Matplotlib; no package was installed for this check. On another workstation, use an environment with the documented dependencies and resolve historical source paths before running.

The starting preservation snapshot covers 5,048 pre-existing files, including prior publication documentation, local results and the external professor/D.1 source documents. It excludes this new directory and generated caches. Preserve the recorded snapshot rather than replacing it to conceal a change. The original study preservation records remain unchanged.

No A/B/C or earlier D files, scientific constraints, feedback settings or initial states are edited. No final submission report, commit or push is part of this task.
