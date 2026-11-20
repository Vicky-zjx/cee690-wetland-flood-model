# Stage D: upstream-level gate feedback pilot

The tested feedback policy is a viable computational candidate for a **bounded, non-freezing, hypothetical-system study**. At a uniform 15 s ODE step, the September 2022 reference and all 16 original non-freezing weather pairs complete, remain within the checked hydraulic domain, satisfy the unchanged dam/wetland limits, and pass the unchanged water-balance check. This is conditional model evidence, not field validation, engineering certification, or proof of safety under arbitrary weather.

The policy avoids the problematic low-head regimes rather than validating them. The original Stage D hydraulic candidate remains unchanged and available as an alternative. A/B/C and all earlier D files are protected by the [preservation manifest](preservation_before.json); the [final verification](preservation_check.json) records the byte-level check. There is no new spillway, storage clipping, safety-limit adjustment, policy optimization, weather generation, or full Monte Carlo analysis.

## Operating policy

Both gates use the same policy, relative to their existing sill elevations. It reads only the accepted upstream pool level at a fixed 120 s controller clock. Hysteresis prevents immediate reversal at a single threshold; actual opening moves continuously at a bounded rate.

| Setting | Declared value |
|---|---:|
| Close when pool depth above sill is at or below | 0.20 m |
| Reopen when depth is at or above | 0.30 m |
| Between these levels | Retain ON/OFF mode |
| ON target at depth at or below 0.30 m | 0.10 m |
| ON target from 0.30 to 0.80 m depth | Linear increase from 0.10 to 0.30 m |
| Target at depth at or above 0.80 m | 0.30 m |
| Maximum opening/closing speed | 0.0005 m/s |
| Initial actual opening | Original 0.0015 m for each gate |
| Unchanged dam safety depth above sill | 1.20 m |

These are explicit hypothetical operational assumptions, selected before the simulations. The [policy basis](policy_basis.md) supplies the equations, USACE references, finite actuator timing, assumptions and low-head buffer calculation. The machine-readable settings are in [policy.json](policy.json). Position limits are implemented by the controller; the unchanged discharge function receives the actual opening and retains its applicability exceptions.

## Actual findings

The canonical nominal run uses the unchanged Stage D hydraulic formulation with customary weir coefficient 2.8, the official `ode4u` solver, a 15 s integration step and a 120 s controller clock. The event is the original 144-hour window beginning 2022-09-09 00:00 UTC.

| Evaluation | Hydraulic completion | Numerical acceptance | Dam/wetland safety |
|---|---:|---:|---|
| September 2022: A/A reference wetland, historical B/B design, no wetland | 3/3 | 3/3 | All checked limits satisfied |
| Original non-freezing scenarios: eight-state wetland | 16/16 | 16/16 | All checked limits satisfied |
| Same scenarios: six-state no-wetland | 16/16 | 16/16 | Both dam limits satisfied; no wetland compartments |
| Four original freezing scenarios | Not run | Not applicable | Retained as liquid-water-domain exclusions |

The 16 pairs use the same saved weather, original seeds, shared physical initial states and policy. Their actual gate histories also agree within the reported numerical tolerance. The wetland system retains its accepted extra 5,500 m³ initial wetland storage; total initial water was not artificially equalized.

The canonical minimum open-gate H/a is **1.98969**, above the Stage D energy-gate boundary of 1.5. Only `energy_gate` and `closed` regimes occur. Therefore no low-head weir equation is extrapolated and the assumed transition is not exercised. Using the original B hydraulics with this same feedback policy also avoids the original interruptions; that supports operating-domain avoidance as the mechanism, not validation of the low-head extension.

In the deterministic reference, Dam 1 now peaks at **100.791729 m**, below its unchanged **101.20 m** limit. The prior D fixed-schedule model peaked near **101.225533 m**, above that limit. Across the canonical deterministic cases, the smallest dam/wetland margin is **0.143216 m**; across the 16 weather pairs it is **0.404256 m**. Safety is checked at accepted output and successful Runge–Kutta evaluations, not certified at every instant between evaluations.

### Upstream safety has a downstream cost

| Deterministic case, new D feedback at 15 s | River peak, m | Peak outlet discharge, m³/s |
|---|---:|---:|
| A/A reference wetland | 97.621748 | 1.048329 |
| Historical B/B wetland design | 97.615517 | 1.034483 |
| No wetland | 97.695105 | 1.211345 |

The B/B label identifies the retained synthetic response-class design selected previously in B; it is not an optimum recomputed for D. Its river peak is about **79.588 mm lower** than the no-wetland peak under the new common policy. However, on the common 120 s grid, switching from D fixed gates to D feedback **raises** its river peak by **86.997 mm**. Relative to the original B fixed case, the net increase is **62.088 mm**. This policy protects upstream storage through more timely release; it is not a downstream-peak optimization.

The [attribution table](results/deterministic_attribution.csv) separates the signed hydraulic change at fixed gates from the subsequent control change under D hydraulics. The additional B-feedback comparison can have different realized gate histories and is labeled accordingly. Failed fixed-schedule stochastic runs do not have invented full-event peaks.

For **15 of the 16 stochastic pairs**, both river maxima equal the initial 97.35 m at time zero. Only `block_bootstrap_00000` has a later higher peak: 97.393277 m with wetlands versus 97.503606 m without, a 110.329 mm difference. The zero contrasts are dominated by the initial condition; they do not establish that wetlands have no routing effect. These 16 cases do not support a general flood-benefit probability or a rainfall-extreme validation claim.

## Numerical checks and coefficient sensitivity

All [46 controller checks](results/policy_checks.json) pass. Full-ODE instantaneous mass identities are consistent to a maximum recorded 1.33e-15 m³/s. The largest canonical nominal relative integrated balance residual is **3.46243e-5**, below the unchanged **1e-4** tolerance.

The nearly dry `block_bootstrap_00008` contains only 0.04 mm rainfall, or 26.4 m³ external rain volume. Its original coarse runs completed hydraulically and safely but failed the integrated balance criterion: the relative tolerance permits only 0.00264 m³ absolute residual. Residuals decreased approximately fourfold when the step halved, consistent with the trapezoidal flux quadrature. At 30 s the no-wetland run still failed. A uniform 15 s nominal rerun passed for every warm case. All earlier failed checks remain in the tables; no tolerance, physical state or weather value was changed.

The delivered study contains **307 evaluation records**, including cold preflight exclusions and repeated numerical/hydraulic arms; it still uses only **20 unique original weather scenarios**. The 120/60/30/15 s checks hold the controller clock fixed. For accepted comparable trajectories, the largest common-grid river difference is 2.81e-8 m for 60→30 s and 7.49e-9 m for 30→15 s. These are numerical consistency comparisons, not a formal convergence proof or empirical validation.

The declared coefficient cases 2.6 and 3.0 also complete safely: 120 s evaluations are supplemented by explicit 15 s refinements for the same nearly dry pair. The deterministic B/B river peak changes by approximately -5.738/+3.840 mm relative to nominal. Because the coefficient also affects the adjoining energy-gate law, this sensitivity remains relevant even though no weir regime is reached. Most stochastic peak comparisons are insensitive only because their maximum is the prescribed initial value. No coefficient or controller parameter was tuned to obtain completion.

## Files and reproducibility

- [README_results.md](README_results.md): four figures, captions and detailed interpretation.
- [results/case_results.csv](results/case_results.csv): every completion, exclusion and failure, with solver termination and mass balance.
- [results/safety_results.csv](results/safety_results.csv): separate dam and wetland peaks, limits and margins.
- [results/physical_outcomes.csv](results/physical_outcomes.csv): river/outlet peaks and timing, gate releases, wetland drainage and overflow. Internal volumes must not be added as independent inputs.
- [results/controller_results.csv](results/controller_results.csv): actual gate movement, mode changes, minimum H/a and executed regimes.
- [results/paired_results.csv](results/paired_results.csv): input/state/grid checks and paired outcomes.
- `results/<suite>/*_controller.csv.gz`: complete command and actual travel histories; valid full trajectories are saved as NPZ alongside per-run JSON diagnostics.
- [results/resolution_comparison.csv](results/resolution_comparison.csv), [results/coefficient_sensitivity.csv](results/coefficient_sensitivity.csv), and [results/source_checks.csv](results/source_checks.csv): numerical, assumption and provenance checks.

The tested environment is the existing `/opt/miniconda3/bin/python3`, with the project's existing NumPy, pandas and Matplotlib dependencies and external official `multivarious` checkout. No replacement solver was installed. The verified solver revision is `06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02`; each suite records its source hashes. The early runner allowed only 120/60/30 s; adding the 15 s option required no model or policy change. Its exact earlier source is retained in [source_snapshots](source_snapshots/README.md), so the earlier recorded hashes remain reproducible.

```bash
cd /Users/zhangjiaxuan/Documents/Codex/Modeling
export MPLCONFIGDIR=/private/tmp/cee690-feedback-matplotlib

# Focused controller checks; no ODE ensemble.
/opt/miniconda3/bin/python3 -B -m D.feedback_control.verify

# Reproduce nominal cases in a NEW output tag; existing tags are protected.
/opt/miniconda3/bin/python3 -B -m D.feedback_control.run_policy --suite all --dt 15 --tag replay_nominal

# Rebuild tables and figures from the delivered saved suites; no simulations.
/opt/miniconda3/bin/python3 -B D/feedback_control/report.py

# Confirm A/B/C and prior D files remain unchanged.
/opt/miniconda3/bin/python3 -B -m D.feedback_control.preservation verify
```

The delivered execution plans also record the lower/upper coefficient and original-B-hydraulics settings. Those arms use `--suite all --dt 120` with `--Cw-US 2.6`, `--Cw-US 3.0`, or `--hydraulics B`, respectively, and separate unused tags. The targeted dry-case refinements are implemented in `refine_dry_case.py`. Replay tags must be unused; the report reads the named delivered suites. Original weather files and solver source must remain available locally.

## Bounded HW D.2 scope

A defensible next scope is **Uncertainty Propagation under Non-Freezing Meteorological Conditions with a Fixed Upstream-Level Gate Policy** for this hypothetical installation. Preserve the policy, safety limits, initial conditions and forcing definitions; retain any future unsupported, numerically invalid or unsafe outcome in the accounting. Gate control is not an optimization variable, and success here does not justify dropping unfavorable future trajectories.

The principal limitations remain hypothetical geometry and actuator behavior, ideal level sensing, uncalibrated discharge coefficients (especially tiny openings), the unvalidated low-head transition, omitted freezing, initial-condition-dominated stochastic peaks, and the limited weather sample. The pilot supports proceeding toward this bounded computational study; it does not establish universal hydraulic resolution or real-infrastructure safety. Full HW D.2 Monte Carlo propagation has not been started. No Git commit or push was made.
