# D.2 consolidated evidence and readiness audit

This is a navigation and consistency audit of completed work, **not the final HW D.2 report**. No new ODE integration, physical experiment, parameter sampling, fitting or optimization was performed. Scientific outputs remain preserved. Read the [study ledger](study_ledger.csv), [scientific consistency checks](scientific_consistency.md), [professor-method and requirement review](methodology_audit.md), and [master figure/table index](FIGURE_TABLE_INDEX.md) together.

## Keep the studies separate

| Study | Primary evidence | What must not be added to its statistical sample count |
|---|---|---|
| Original stochastic weather pilot | 16 valid non-freezing pairs from 20 original weather IDs; four cold IDs retained as exclusions | The second integration resolution, deterministic controls, or a presumed random heavy-storm population |
| Fixed-storm sensitivity | 27 deterministic combinations of three multiplier levels | Probability weights or a claim that these are 27 iid draws |
| Parameter Monte Carlo | 64 matched vectors under the reference September 2022 storm | The same vectors' pilot repeats or numerical refinements |
| Historical storms | Five fixed 144 h windows × the same 64 vectors = 320 pairs | Adding the original 64 again: the benchmark is already included. Events are not independent random flood-frequency draws |
| Saved-result robustness checks | Reanalysis of the above tables/trajectories; 87 exact zero primary outcomes diagnosed | Any new physical run, parameter draw or causal experiment on the initial wetland water |
| Antecedent conditions | 12 physical pairs: two events × two prefixes × three predeclared configurations | Its one refinement pair, six reused no-prefix references, or treating three selected configurations as a Monte Carlo distribution |

No cross-study pooled probability estimate is formed. The same parameter draw identifier in different events means matched physical parameter values, not independent replications. Earlier README statements such as “next experiment” or “not yet implemented” describe their original study's completion time; the present ledger and linked completed studies supply the current status without rewriting historical evidence.

## Interpretation rules retained

- Every primary full-horizon maximum includes its initial endpoint. The robustness study's later rainfall-response window and the antecedent study's original 144 h event window are separate secondary/conditional definitions, not replacements for the full metric.
- Lower peaks do not imply shorter duration above a reference level, lower recession stages, or lower cumulative outlet volume. No justified downstream flood threshold has been supplied; reference-duration and excess-volume metrics are not flood probabilities or avoided flood volumes.
- The wetland model starts with 5,500 m³ more water. Antecedent event windows use the actual carried storage contrast. Mass accounting does not identify the origin of water or isolate a causal effect of that initial stock.
- The uncapped outlet rating makes paired stage and discharge peak differences mathematically dependent through river resistance. They are not independent validation evidence.
- Independent log-uniform parameter multipliers on [0.5,2] are approved exploratory assumptions, not calibrated D.1 distributions. Uniform reweighting has ESS about 44/64, but outcome contributions can be concentrated in fewer draws. Ties are retained in count/N ECDFs and fixed-threshold binomial intervals.
- Antecedent numerical refinement covers one declared physical pair. It supports that pair and does not prove convergence for every configuration or equilibrium after 24/72 h.

## Checks and method status

The [saved-result audit](cross_study_audit.json) records comparisons, overlaps and source hashes. A separate [per-run diagnostic cross-check](diagnostics_crosscheck.json) passes 9,652 comparisons between reported tables and 944 saved run JSON records, including copied benchmark evidence and numerical/control records. This is an evidence-file count, not 944 independent simulations. The existing component checks were rerun with only their output destination privately redirected: **46 controller  + 123 hydraulic/RHS  + 9 response-metric checks = 178**, all passing and matching the original saved reports. This ran component function evaluations, not ODE trajectories. [Current component evidence](component_checks.json).

The professor's physical-analysis, sampling, Monte Carlo, ECDF and scatter methods are mapped to actual project code and saved executed outputs. A polynomial/surrogate fit to D.2 Monte Carlo responses is not present. The [method review](methodology_audit.md) distinguishes the example's optional code paths, the separate Stage C weather fit, and requirements stated in the locally supplied final-project draft. It does not invent a dedicated D.2 rubric or claim that all final-project requirements are complete. No missing experiment or surrogate was added during this audit.

All 30 existing D.2 PNGs were visually reviewed; [review notes](figure_review.md) preserve units, denominators, signs, scopes and caption limitations. Existing numerical curves and figures were not regenerated merely for publication.

## Reproduce this audit

Run from the repository root with the documented Stage D environment:

```sh
python -B D/d2_analysis/consolidation/cross_study_audit.py
python -B D/d2_analysis/consolidation/diagnostics_crosscheck.py
python -B D/d2_analysis/consolidation/component_checks.py
python -B D/d2_analysis/consolidation/preservation.py verify
```

The first two commands check saved evidence, not simulations. Both passed in an [isolated checkout containing only indexed files](reduced_checkout_check.json): 3,355 cross-study checks and 9,652 per-run diagnostic checks. This used the existing installed Python dependencies and did not test a fresh dependency installation or the omitted full-array audits. Component checks additionally need the existing official solver/input environment. The preservation check is local and compares the current workstation with its task-start snapshot, allowing only three navigation files and the existing ignore file to change. It is not a fresh-clone integrity test. Historical preservation records are untouched.

The [publication audit](publication_audit.json) records included/excluded file counts, credential screening, preservation and reduced-checkout checks. [Stage D reproduction guidance](../../reproducibility_inputs/README.md) distinguishes included plans/forcing/tables/figures from local dense arrays and source caches. The reduced repository is not claimed to reproduce every full-trajectory audit without those local artifacts. No source observation, discarded failure, initial peak or unsafe outcome is silently replaced by a fabricated value.
