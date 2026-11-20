# Stage D publication and reproduction guide

This guide describes the repository snapshot of the completed analyses. It does not revise the physical model, saved numerical results, parameter vectors, uncertainty assumptions or original validation records. Publication preparation ran component and saved-result checks, not new ODE simulations.

## Included and local files

The repository includes all Stage D Python source, Markdown documentation, study/configuration JSON files, per-run JSON diagnostics, numerical-resolution summaries, CSV tables, figures and source/hash records. The 27-case sensitivity study, the original 64-vector Monte Carlo study and the five separate historical-event ensembles remain distinct. The four selected historical weather CSVs and complete compressed selection inventory are included.

Ten unchanged trajectory fixtures are included for the final historical hydrographs and mechanism analysis: both models for `mc_00000` in the four added events, and the original September 2022 `pilot_mc_00000` 15 s pair. These are about 45.9 MB in total. They are selected by the existing figure readers, not by favorable outcomes. All other dense Stage D NPZ trajectories and `*_controller.csv.gz` histories remain local under `.gitignore`; they were not deleted. Their summary outcomes, failed/unsafe classifications, existing numerical checks and provenance remain published. The older Stage D development evidence is retained, including unsuccessful hydraulic alternatives.

This bundle adds exact copies of five previously ignored inputs required by the existing loaders:

| Original restore path | Purpose |
|---|---|
| `B/results/no_spillway_revision/reference.npz` | Frozen original gate schedule and its source hash |
| `B/results/wetland_only/wetland_best_120s.npz` | Accepted B-selected wetland design |
| `B/results/wetland_only/fixed_reference_120s.npz` | Original A/A reference used by historical comparison paths |
| `B/work/event_data/proposed_forcing_120s.csv` | Original September 2022 forcing, unchanged |
| `B/work/event_data/proposed_forcing_manifest.json` | Approval, units, provenance and conservation metadata |

[manifest.json](manifest.json) records byte counts and SHA-256 hashes. Original A/B/C directories remain unchanged and retain their existing ignore rules. This is an input bundle, not a second model or a replacement for source observations. The professor's example, submitted D.1 PDF, raw meteorological archives, entire official solver checkout, credentials, caches and virtual environments are not bundled.

## Environment and restoring inputs

The verified local environment used Python 3.13.5 and the direct dependencies in [D/requirements.txt](../requirements.txt). A fresh installation has not been tested as part of this publication task. From the repository root, a separate environment can be prepared with:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r D/requirements.txt
```

The official [multivarious source](https://github.com/hpgavin/multivarious) remains an external dependency at `A/work/multivarious`, revision `06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02`. Existing loaders verify its revision, source hashes and lack of tracked modifications. No similarly named unverified PyPI package is substituted. If that checkout is absent:

```sh
git clone https://github.com/hpgavin/multivarious.git A/work/multivarious
git -C A/work/multivarious checkout 06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02
```

The project imports this pinned checkout directly. The clone/installation commands above were documented, not executed during publication. Existing local dependencies were used for the checks. Rebuilding historical weather selection additionally requires the existing Stage C environment and cached observations described in [event_selection.md](../d2_analysis/historical_storm_robustness/event_selection.md); those complete raw archives are not needed merely to inspect the included results.

Check the frozen input bundle without writing anything:

```sh
python -B D/reproducibility_inputs/restore_inputs.py
```

On a new checkout, restore only missing inputs to their original ignored locations:

```sh
python -B D/reproducibility_inputs/restore_inputs.py --restore
```

The helper verifies every source and existing destination before copying. It refuses to overwrite an existing different file, never rewrites a configuration, and runs no simulation. Both check-only behavior on this workstation and restoration to an isolated temporary directory were tested.

## Reproduction levels and limits

**Inspect the completed work:** all final figures, captions, outcome tables, plans, safety diagnostics and findings are directly readable. No simulation is required. Start at [Stage D](../README.md).

**Regenerate the final parameter-study summaries and plots:** the existing commands read included JSON/CSV records, with the pinned official plotting utilities for the second command:

```sh
python -B D/d2_analysis/parameter_uncertainty/analyze_results.py
python -B D/d2_analysis/parameter_uncertainty/run_official_analysis.py
```

**Regenerate final historical-storm figures/mechanism checks:** the required ten trajectory fixtures and numerical inputs are included. The existing `analyze_results.py` and `mechanism_check.py` commands are documented in that study's README. They write derived tables/figures and were deliberately not rerun for this push. Saved manifests contain original absolute paths under `/Users/zhangjiaxuan/Documents/Codex/Modeling`; on a different workstation those paths must be resolved in a separate relocation/configuration copy before these commands can use the supplied forcing. Preserve the original manifests and hashes as provenance. The unmodified runners are not claimed to be immediately portable to an arbitrary directory.

**Rerun simulations or full trajectory audits:** use the original frozen plans and existing drivers, not new sampling. Existing phase directories are protected against overwrite. Full-local-data audits such as `historical_storm_robustness/independent_review.py` require all locally retained dense arrays; the published passed audit JSON records the completed check, but the reduced Git checkout cannot rerun that full audit without restoring or regenerating the omitted arrays. Older weather-pilot plotting and nominal adapter-regression paths likewise need their omitted local trajectories/controller records. They remain reproducible from the declared inputs and model after those local artifacts and path dependencies are supplied; this publication did not regenerate them or verify a fresh-machine rerun.

**Repeat event selection:** requires the complete original Stage C standardized calendar and cached NOAA/NASA raw responses. The selected weather rows, inventory, public source URLs, interval audit and hashes are published; the large raw archive remains local. Downloaded source bytes and retrieval metadata must be checked before claiming exact source reproduction.

## Preservation records and publication checks

Historical `preservation_before.json` / `preservation_check.json` files describe checks made when each experiment completed. Some refer to workstation-wide files, earlier dense trajectories and two private Downloads documents. They are unchanged historical evidence, not fresh-clone integrity checks. This publication updates only navigation/reproduction documentation and ignore rules outside the frozen scientific artifacts, so rerunning an old whole-project preservation guard after publication can correctly flag those administrative changes or absent local files. Do not rewrite old manifests to conceal that distinction.

Before publication, 178 existing component checks passed (46 controller, 123 hydraulic/RHS and nine response-metric checks), saved controller/hydraulic reports matched, and the 27/64/five-by-64 summary counts and matched vectors were verified. The complete historical trajectory audit and the pre-publication 3,130-file preservation check passed. No fitting, optimization, meteorological generation or ODE integration was performed for publication.
