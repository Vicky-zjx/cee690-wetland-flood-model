# Statistical robustness checks using saved results

This analysis reads the existing five fixed-event, 64-pair result matrices and the existing September 2022 27-point parameter grid. It does not run an ODE, generate a new physical parameter vector, fit a surrogate, revise a parameter or alter an earlier result. All 320 saved pairs remain in the analysis: 320 valid, zero invalid and zero upstream safety violations. Zero or unfavorable response metrics are not excluded.

## ECDF audit and zero effects

All 42,030 saved fixed-grid CDF rows match fresh count/N and Clopper–Pearson calculations. All 65 possible success counts for N=64 agree with SciPy's exact-binomial interface to a maximum absolute endpoint difference of 4.35e-13. The original official rank/(N+1) arrays remain correct as plotting positions, but they are not the ties-aware empirical CDF. The statistical audit checks all saved outcome-summary metrics as well as the peak-stage differences.

| Fixed event | Positive / zero / negative peak differences | Empirical median (mm) | Exact-zero fraction; pointwise 95% CP interval |
|---|---:|---:|---|
| Sep 2022 benchmark | 64/0/0 | 94.0549 | 0.0000; [0.0000, 0.0560] |
| Aug 2023 | 55/9/0 | 86.8593 | 0.1406; [0.0664, 0.2502] |
| Jul 2022 | 53/11/0 | 69.5981 | 0.1719; [0.0890, 0.2868] |
| Aug 2019 | 27/37/0 | 0.0000 | 0.5781; [0.4482, 0.7006] |
| Aug 2024 | 34/30/0 | 7.3411 | 0.4688; [0.3428, 0.5977] |

Every exact-zero peak difference coincides with both models attaining the same full-window maximum at their initial river level. This does not mean the hydrographs or all other response metrics coincide. No near-zero tolerance was used to manufacture zero counts. Because no negative differences were observed, the empirical zero mass and empirical F(0) have identical counts and numerical intervals here; their population meanings remain different. Absence of an observed negative value does not establish zero population probability below zero. The pointwise F(0) intervals for 2019 [0.4482, 0.7006] and 2024 [0.3428, 0.5977] both include one-half. Their empirical medians are valid descriptions of these 64 draws, but the data do not precisely distinguish a nonpositive population median from a positive one under the assumed law.

The CDF intervals describe Monte Carlo sampling uncertainty at predeclared fixed thresholds under an assumed parameter distribution. They are not the spread of physical outcomes, not simultaneous CDF bands, and not probabilities of historical weather occurrence. The five selected weather windows are fixed comparisons, not an iid weather sample. No 320-pair pooled risk bound is computed because events reuse the same 64 vectors and are not a random weather sample. See the [SciPy exact interval API](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats._result_classes.BinomTestResult.proportion_ci.html).

`statistical_sampling_uncertainty.csv` separates the sample standard deviation and 5th–95th outcome quantiles from the estimated Monte Carlo standard error of each mean, s/√64. The latter is a precision diagnostic for iid draws under the assumed law, not a physical variability interval or a verified confidence interval for a tail quantile. Estimates across events are dependent because their vectors are matched.

## Sample-size diagnostic

Original-order prefixes use 16, 32, 48 and 64 saved vectors. A separate diagnostic uses 1,000 fixed-seed random permutations (seed 69020261119), taking the same nested prefixes for every event. No replacement is used within a subset; the identical subset indices are shared across events. All permutations and replicate statistics are saved. The middle 95% ranges describe possible subsets of this specific finite set, not confidence intervals for an unknown population. At N=64 the range collapses mechanically; it does not establish convergence.

| Fixed event | Mean ΔH, original prefixes 16 / 32 / 48 / 64 (mm) | 48-draw subset range for mean ΔH (mm) |
|---|---:|---:|
| Sep 2022 benchmark | 96.19 / 91.61 / 92.58 / 91.57 | [86.62, 96.90] |
| Aug 2023 | 104.21 / 95.93 / 95.29 / 91.04 | [82.39, 99.41] |
| Jul 2022 | 86.45 / 78.49 / 77.62 / 73.81 | [66.16, 81.27] |
| Aug 2019 | 42.16 / 32.01 / 29.00 / 26.18 | [20.50, 31.22] |
| Aug 2024 | 61.01 / 48.94 / 44.97 / 41.00 | [33.63, 47.84] |

The subset tables also retain median and 95th-percentile peak differences, zero fractions, discharge differences, duration differences and total outlet-volume differences. These diagnostics can reveal dependence on which saved vectors are included, but cannot discover unvisited parameter-space behavior or verify a rare-event tail.

## Parameter-law sensitivity without new simulations

The original law is independent log-uniform on [0.5, 2] for three categories: upland infiltration, a single common routing-rate multiplier, and river resistance. Independent uniform multipliers on the same support are compared by self-normalized importance weights proportional to the product of the three multipliers. This changes the assumed input measure, not any physical parameter vector or response.

Two modest descriptive dependence stresses retain the original log-uniform marginals and use latent-normal correlations −0.25 and +0.25 between infiltration and common routing; river resistance stays independent. The upland connection motivates examining dependence, but neither its sign nor magnitude is measured. Symmetric stress values were specified as hypothetical settings before reweighted outcome summaries, not estimated from the response. They do not represent new calibrated priors. No dependence involving river resistance is tested, so small changes in these two stresses cannot establish that independence assumptions are generally unimportant.

For u_j=ln(m_j/0.5)/ln(4), z_j=Φ⁻¹(u_j), the copula density ratio is det(R)⁻¹/² exp[−zᵀ(R⁻¹−I)z/2]. The implementation independently checks this formula against the [SciPy multivariate-normal log-density](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.multivariate_normal.html). Full formulas, per-draw weights and marginal checks are saved in the statistical law tables.

| Law | Weight ESS / 64 | Largest normalized weight | Mean raw density ratio |
|---|---:|---:|---:|
| Independent log-uniform | 64.00 | 0.0156 | 1.0000 |
| Independent uniform | 43.99 | 0.0498 | 0.9643 |
| Log-uniform; latent rho = -0.25 | 57.72 | 0.0448 | 1.0610 |
| Log-uniform; latent rho = +0.25 | 60.41 | 0.0257 | 0.9520 |

These estimates reuse only 64 proposal draws. Self-normalization is finite-sample biased; a high ESS cannot establish adequate coverage of every output tail. ESS=1/Σw² measures weight concentration only; it is not substituted for N in binomial intervals. No weighted confidence bands are computed. The reweighting and its limitations follow [Owen, Monte Carlo theory, importance-sampling chapter, Sections 9.2–9.3](https://artowen.su.domains/mc/Ch-var-is.pdf).

| Fixed event | Baseline mean ΔH (mm) | Uniform mean ΔH (mm) | Latent rho −0.25 / +0.25 mean ΔH (mm) | Uniform zero mass |
|---|---:|---:|---:|---:|
| Sep 2022 benchmark | 91.573 | 107.881 | 91.119 / 90.785 | 0.0000 |
| Aug 2023 | 91.045 | 117.158 | 92.804 / 88.927 | 0.0599 |
| Jul 2022 | 73.814 | 97.922 | 74.903 / 72.187 | 0.0728 |
| Aug 2019 | 26.178 | 38.840 | 26.652 / 25.246 | 0.4127 |
| Aug 2024 | 41.003 | 59.616 | 41.473 / 39.958 | 0.2991 |

The 2019 exact-zero mass changes from 0.5781 under equal weights to 0.4127 under the uniform-multiplier stress; its weighted inverse-CDF median becomes 20.651 mm. Thus the zero/positive median distinction is sensitive both to finite-sample uncertainty and to the uncalibrated parameter-law assumption. This is not evidence that one of the laws is more realistic.

All reweighted peak differences stay nonnegative because every original saved peak difference is nonnegative and every weight is positive. That sign cannot be an independent validation of another parameter law or evidence about parameter combinations absent from the saved set. Output quantile ranges describe variability across hypothetical parameter values. Their estimation uncertainty is a different quantity. A weighted inverse-CDF median uses the first value reaching cumulative weight 0.5, which can differ from the interpolated even-N empirical median; both definitions are recorded explicitly.

## Parameter patterns and limits of attribution

Saved Spearman rank associations are independently reproduced, without significance testing or a new fit. They are marginal descriptions within each event, not causal effects or Sobol indices. The saved 27-point benchmark grid additionally supplies nine matched endpoint contrasts per factor: change the factor from 0.5 to 2 while holding the other two exactly fixed. No extra grid was evaluated for another event.

| Factor | Outcome | Positive / zero / negative endpoint changes over nine contexts | Endpoint range |
|---|---|---:|---:|
| m_K_u | ΔH | 8/0/1 | -30.17694 to 25.80374 mm |
| m_K_u | ΔQ | 8/0/1 | -0.13412 to 0.08483 m³/s |
| m_routing | ΔH | 4/0/5 | -20.71476 to 28.07183 mm |
| m_routing | ΔQ | 4/0/5 | -0.09207 to 0.12476 m³/s |
| m_R_r | ΔH | 9/0/0 | 96.18539 to 135.33685 mm |
| m_R_r | ΔQ | 4/0/5 | -0.04526 to 0.09505 m³/s |

The grid shows why one rank correlation is insufficient: increasing infiltration raises the stage/discharge peak difference in eight contexts but lowers it in one; increasing common routing raises it in four and lowers it in five. River resistance raises ΔH in all nine contexts, while its ΔQ effect changes sign. These are controlled finite differences within the saved hypothetical benchmark model, not measured field sensitivities or general monotonicity results.

The outlet rating is uncapped at these peaks, so ΔH_peak=R_r ΔQ_peak for a paired realization. A strong association of resistance with ΔH is partly the conversion from discharge to stage; it must not be read as an independent increase in hydraulic mitigation. Common routing remains a single category rather than three separately perturbed rates.

All duration and volume signs remain available in the full outcome tables. Negative differences are normal modeled outcomes, not invalidity. In particular, all four additional storms have higher total outlet volume in the wetland arm. The wetland initial state contains 5,500 m³ more water, so the 144-hour volume contrast is not pure avoided floodwater. H0 is an initial-state reference, not a flood threshold. Initial states, event windows and the benchmark/non-benchmark weather-processing differences are unchanged.

## What remains unresolved

The saved data support a conditional robustness description, not calibrated distributions, a proven sample-size convergence result, or an environmental risk probability. A useful next approved experiment would append an independent batch of parameter draws under the original law, using exactly the same five fixed storms and paired model settings, and prospectively compare the mean, zero mass and tail quantiles. This would test whether the current estimates persist beyond reuse of the same finite 64-point set. It is a proposal only: no new draws or ODE evaluations were made here.

## Figures and reproducibility

- `plots/04_ecdf_zero_atoms.png`: exact jumps at zero and pointwise binomial intervals; whole-window maxima retained.
- `plots/05_sample_size_stability.png`: matched finite-set subset diagnostics; ranges are not population confidence intervals.
- `plots/06_parameter_law_checks.png`: reweighted means, zero mass and weight concentration; no weighted intervals.

Run `python D/d2_analysis/robustness_checks/statistical_checks.py` from the repository root. Inputs are existing CSV/JSON files; outputs are confined to this directory. `tables/statistical_audit.json` records source hashes, library versions and statistical checks. The original results are never rewritten.
