# Pointwise CDF intervals when outcomes contain atoms

The historical-storm study can produce an exact initial-state maximum in both models, yielding a point mass at zero paired peak reduction. These are valid outcomes. They are retained without jittering, removal, or replacement. The official continuous-order-statistic interpretation of `plot_ECDF_ci` does not automatically apply to such atoms.

The distribution figure therefore makes an explicit statistical adaptation while retaining the professor's unchanged utility:

- The solid curve is the exact ties-aware empirical CDF, `count(Y <= x) / n`; its first jump starts at zero.
- The dotted curve uses the unchanged official `rank/(n+1)` output. The utility is still called separately for every event and metric.
- The shaded limits are **95% two-sided Clopper–Pearson intervals for the CDF at fixed numerical grid locations**, not the official order-statistic limits.
- The official sorted values, plotting positions and Beta order-statistic limits remain in `tables/official_pointwise_ecdf_arrays.csv`, with `displayed_as_shading=false`. No exact continuous-coverage claim is made for tied/atomic outputs.

For a fixed location x, each indicator `1{Y_i <= x}` is Bernoulli with success probability `F(x)`, whether or not Y has atoms. With `k = count(Y_i <= x)` among n valid iid parameter realizations, the displayed endpoints are

```
lower = 0                                      if k = 0
        Beta.quantile(0.025; k, n-k+1)          otherwise
upper = 1                                      if k = n
        Beta.quantile(0.975; k+1, n-k)          otherwise
```

These are exact-binomial, generally conservative **pointwise** intervals. They are not simultaneous bands over all x, metrics or events. The shaded connection between fixed grid points is a visual guide; exact pointwise coverage is asserted only at the specified grid locations. The two-sided exact method is documented by [SciPy's `BinomTestResult.proportion_ci`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats._result_classes.BinomTestResult.proportion_ci.html). The implementation compares its vectorized Beta-quantile endpoints against `binomtest(k, n).proportion_ci(0.95, method='exact')` for every possible k at each sample size used. No hypothesis-test p-value is reported or used to select a result.

The following grids are fixed in the new analysis code before ensemble completion. They are not observed order statistics, sample quantiles, flood thresholds, or a procedure for changing physical limits:

| Output | Fixed range | Increment | Locations |
|---|---:|---:|---:|
| Wetland/no-wetland peak stage | 97–101 m | 0.002 m | 2,001 |
| Wetland/no-wetland peak discharge | 0–3 m³/s | 0.005 m³/s | 601 |
| Paired stage difference, no-wetland minus wetland | −1–1 m | 0.001 m | 2,001 |
| Paired discharge difference, no-wetland minus wetland | −3–3 m³/s | 0.005 m³/s | 1,201 |

The plotting window can zoom to the observed values for readability, without changing those computation grids. Values are never clipped to a grid; coverage is not extrapolated beyond its range. Counts, endpoints, denominators and fixed locations are saved in `tables/fixed_grid_pointwise_cdf_intervals.csv.gz`. Runtime comparison with the installed SciPy implementation is saved in `tables/clopper_pearson_implementation_checks.csv`.

Each event remains a separate fixed forcing case. The original 64 iid parameter vectors are reused across events, so comparisons between events are dependent and matched by sample ID. There are no weather-occurrence probabilities or empirical flood-risk estimates. The intervals address finite parameter-sampling uncertainty conditional on the explicit exploratory laws, physical model, event forcing and initial conditions. They do not account for calibration error, model discrepancy, uncertain event frequency or uncertainty in the assumed laws.

If any pair is invalid, its cause and planned denominator remain reported. A plot based on the retained complete pairs is explicitly conditional on a valid pair; it does not describe the unconditional output law with missing invalid trajectories. Physically valid unsafe cases, negative effects, initial peaks and zero effects are not filtered. No physical equation, controller, parameter vector, forcing file or saved simulation is changed by this statistical adaptation.
