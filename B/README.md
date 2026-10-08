# B Current wetland optimization and dynamic comparisons

[Project overview](../README.md) · [A Initial model and baseline](../A/README.md)

B is the current project: wetland design optimization under fixed gate operations, together with the no-wetland dynamic comparison. This page presents the final accepted model, method, results, figures and files for B.

The final accepted analysis asks how wetland design affects the downstream river peak when gate operation is fixed. It includes an eight-state wetland design search and a six-state no-wetland dynamic comparison. Both use the no-dam-spillway scope, the same prescribed gate openings and the approved regional event-weather reconstruction.

## Model and study conditions

The two wetland branches follow:

```text
Upland -> Dam pool 1 -> Gate 1 -> Wetland 1 -> River
Upland -> Dam pool 2 -> Gate 2 -> Wetland 2 -> River
```

The model retains soil and groundwater pathways and direct upland-to-river runoff. Dam pools have no independent spillway path; wetland normal drainage and wetland overflow remain represented. Water storage is not clipped to enforce a safety limit.

The installation footprint is an assumed 0.66 km². The event covers 2022-09-09 00:00 through 2022-09-15 00:00 UTC. Regional station observations, solar grid estimates, spatial transfers, trace-zero lower bounds and assumed boundary/initial conditions retain their separate evidence labels. USGS 04087088 is contextual evidence only. The [source and assumption document](data_sources_and_assumptions.md) defines the accepted scope and its limitations.

## Optimization problem and numerical method

The submitted B.1 formulation defines the peak-stage objective and includes gate histories together with wetland designs. The accepted analysis here fixes those gate histories and optimizes only wetland designs, as explicitly authorized.

The wetland designs are `Di=(tau_i,rho_i)`. The synthetic categories A/A, A/B, B/A and B/B are enumerated, and the continuous decision vector contains only the two cover fractions `[rho1,rho2]`, each in [0,1]. Categories and covers determine the paired hydraulic-resistance and ET response. Wetland geometry remains fixed, and the selected design is constant during each simulation.

```text
Fixed inputs: g1(t), g2(t), weather, shared geometry, boundaries and initial conditions
Decisions: D1 and D2
Objective: minimize the maximum simulated river stage Hr(t) over the event
```

The official multivarious `ode4u` solver integrates each candidate, and the project-selected official `nms` routine searches the covers. Four category pairs and three initial cover guesses give 12 searches. Seven met the convergence conditions and five reached the evaluation budget. The reported best feasible design is the best found in these bounded searches, not a proven global optimum.

Nonnegative storage, soil capacity, fixed upstream safety levels and shared applicability limits remain enforced. The four maximum allowable levels are 101.20/101.00 m for the pools and 98.90/98.80 m for the wetlands, all in the assumed model datum. G5 is checked during RHS/RK evaluations as hydraulic applicability for the prescribed gate schedule; it is excluded from the wetland-only OPT residual vector. The remaining 18 residuals and native cover bounds retain their applicable roles. Physical river peaks are reported separately from numerical penalty values.

## Accepted wetland designs

| Design | Synthetic categories | Cover 1 | Cover 2 | River peak at 30 s |
|---|---|---:|---:|---:|
| Reference wetlands | A/A | 0.65 | 0.60 | 97.558443 m |
| Best-found wetland design | B/B | 0.2714117624187232 | 0.2141279193188096 | 97.553431 m |

The best-found design lowers the reference peak by **5.012 mm** under identical prescribed gate operation. Its lower effective drainage resistance accompanies more normal drainage and less concentrated wetland overflow. The result does not imply that more vegetation cover or greater stored volume always reduces the downstream peak. A/B are synthetic response classes, not identified plant species or local vegetation recommendations.

The [full optimization report](docs/12_wetland_only_optimization.md) provides the search results, margins, mechanisms, solver termination and numerical-resolution checks.

## Dynamic no-wetland comparison

The six-state counterfactual removes the wetland compartments and routes each gate directly into the river. Rainfall over the former wetland footprints enters the river immediately; the river's storage-stage and evaporation areas remain unchanged. Gate discharge is calculated from internal river tailwater rather than replayed from a saved wetland case. This is a dynamic bypass assumption, not a calibrated replacement land-cover scenario.

| Fixed case | Dynamic states | River peak at 30 s | Bypass peak minus case peak |
|---|---:|---:|---:|
| No-wetland dynamic bypass | 6 | 97.637815 m | 0 mm |
| Reference wetlands | 8 | 97.558443 m | 79.373 mm |
| Best-found wetland design | 8 | 97.553431 m | 84.385 mm |

This fixed comparison does not run OPT. The 79.373/84.385 mm differences include the represented wetland storage, routing, seepage and ET mechanisms. Only the 5.012 mm reference-to-best difference is the additional wetland-design optimization improvement.

The primary wetland cases contain 5,500 m³ more initial water because their wetland compartments retain their approved initial storages. A separate diagnostic sets those two storages to zero without changing the shared initial states. It changes the reference/best peaks by approximately 0.2875/0.2217 mm and preserves the qualitative comparison.

## Verification and interpretation

The fixed configurations were checked at 120, 60 and 30 seconds using the same underlying forcing and gate functions. All reported comparison cases pass the applicable feasibility and hydraulic checks. The largest 30-second integrated water-balance residual is approximately 0.00631 m³. Actual gate releases agree within floating-point precision. The software suite passes 100 tests, including four focused contracts for the standalone analysis interface. The fixed reference, best-found wetland design and six-state bypass have been reproduced at 120 seconds with exact agreement in their saved time, state, input, output and balance arrays. Calling the public analysis function also exactly reproduces the reference/best objective and all 18 optimizer residuals. Saved numerical files are not overwritten.

The event maxima occur inside the selected window, but water remains stored and the river rises slightly at the endpoint. These conditional, uncalibrated simulations do not establish complete recession, measured Milwaukee flood reductions, engineering certification or universal wetland benefits. The [full dynamic comparison report](docs/13_no_wetland_counterfactual.md) provides the initialization sensitivity, conservation accounting and limits.

## Direct analysis and design bounds

[The analysis guide](docs/analysis_and_boundary_checks.md) provides the standalone `analysis(v,c)` interface. For one fixed category pair, `v=[rho1,rho2]`; the constants contain the existing gate schedule, forcing, initial conditions and fixed parameters. The function returns the unchanged dimensionless conditioned objective and an 18-element constraint vector. Its diagnostic record separately exposes the physical river peak, validity and feasibility. Both direct calls and the optimizer use this same implementation.

```python
from B.wetland_only_study import make_analysis_constants, analysis
c = make_analysis_constants(("type_a", "type_b"))
f, g = analysis([0.0, 1.0], c)
print(c.last_evaluation["physical_peak_m"])
```

The [design-bound table](tables/analysis_checks/boundary_checks.md) contains 32 actual dynamic evaluations: four coverage corners for all four category pairs, each at 120 and 60 seconds. All 32 complete and satisfy the unchanged applicability and feasibility checks. The maximum coarse/fine peak difference is 0.010340 mm. These checks use 34 integrations including two saved-case reproductions; they do not run OPT or replace the existing best design. Endpoint tests do not prove stability or feasibility throughout the interior of the design box.

Boundary verification figures are in [plots/analysis_checks/](plots/analysis_checks/):

- [River-stage time histories](plots/analysis_checks/01_boundary_river_peaks.png): four coverage-corner curves per class pair, with peak markers, from saved 60 s trajectories.
- [Safety and gate applicability](plots/analysis_checks/02_boundary_safety_margins.png): minimum headroom to the four upstream limits through time for each corner, plus the two shared G5 histories.
- [Numerical checks](plots/analysis_checks/03_boundary_numerical_checks.png): envelopes of matched-time 120/60 s stage differences and cumulative water-balance residuals; sampled peak differences remain a separate statistic.

All three figures are regenerated from saved checks by `python3 -B -m B.wetland_analysis_checks --report`; this command runs no simulation or optimization.

The [assumption table](data_sources_and_assumptions.md) states the actual resistance/ET endpoints, their units and interpolation, and the basis for full coverage bounds [0,1]. Those coefficients remain explicit synthetic assumptions, not locally measured vegetation properties.

## Final files and reproduction

| File or directory | Purpose |
|---|---|
| [no_spillway_model.py](no_spillway_model.py) and [configuration](no_spillway_config.json) | Accepted eight-state equations and fixed conditions |
| [wetland_only_study.py](wetland_only_study.py) and [configuration](wetland_only_config.json) | Fixed-gate wetland optimization |
| [wetland_counterfactual.py](wetland_counterfactual.py) | Six-state dynamic bypass and fixed comparisons |
| [wetland_analysis_checks.py](wetland_analysis_checks.py) and [analysis guide](docs/analysis_and_boundary_checks.md) | Standalone analysis verification and full coverage-bound dynamic checks |
| [Optimization report](docs/12_wetland_only_optimization.md) | Full wetland design analysis |
| [Comparison report](docs/13_no_wetland_counterfactual.md) | Full three-case analysis |
| [Optimization table](tables/wetland_only/optimization_summary.csv) | Search designs, physical peaks and termination |
| [Comparison table](tables/wetland_counterfactual/comparison.md) and [CSV](tables/wetland_counterfactual/comparison.csv) | Primary and sensitivity results |
| [Mechanism table](tables/wetland_counterfactual/mechanisms.csv) | Flow volumes, timing, ET and remaining water |
| [Wetland figure](plots/wetland_only/01_wetland_design_mechanisms.png) and [comparison figure](plots/wetland_counterfactual/01_simulated_river_stage.png) | Simulated trajectories and mechanisms |

Run commands from the project root in the existing environment with the verified official solver. The detailed reports provide the full numerical execution instructions. Existing completed runs are preserved by their execution guards. To regenerate only the final reports, figures and tables from saved local results:

```sh
MPLCONFIGDIR=B/work/mpl-cache python3 -B -m B.wetland_only_study --report
MPLCONFIGDIR=B/work/mpl-cache python3 -B -m B.wetland_counterfactual --report
python3 -B -m B.wetland_analysis_checks --report
```

Numerical arrays, source locks and callback histories are local under `B/results/`; source observations and processing files are local under `B/work/`. Both stages share the official solver checkout at `A/work/multivarious/`. Saved configuration and provenance JSON files retain their original path labels and recorded hashes. [project_paths.py](project_paths.py) resolves those labels to the present A/B locations and checks recorded import/path-only source changes against the local relocation manifest. Numerical outputs and scientific values are unchanged.

The [A environment guide](../A/README.md#install-and-run) documents the verified dependency revision and installation status. The detailed final reports above are the entry points for the accepted analysis; this guide does not list individual edits or intermediate study attempts.
