# Meteorological disturbance model and validation

This module fits and diagnoses historical precipitation, air temperature, solar radiation and wind speed. The retained workflow is independent of the A/B ODE models. It contains custom NumPy/SciPy statistical methods; official FIT/mimo_srs has not yet been applied. These methods are retained as shared fitting work pending the detailed C.1/C.2 specification.

## Data and clocks

The standardized input is an hourly, timezone-aware UTC table with columns `timestamp_utc`, `precipitation_mm_per_h`, `temperature_C`, `solar_W_per_m2`, and `wind_m_per_s`. A timestamp starts `[t,t+1 h)`. Source-specific meanings, overlap handling, station identifiers, quality flags and missingness are documented in the data-source and quality reports. Training uses 2019–2023 only; chronological held-out validation uses 2024. The generator does not fill missing records. Nonconsecutive hours are never VAR neighbors or contiguous bootstrap blocks.

Seasonal and diurnal features use `America/Chicago`, including daylight saving time. NOAA station observations and NASA POWER gridded radiation remain distinct evidence types. A zero derived from a trace report is a **quantified lower bound**, not a verified absence of precipitation. Quantified-zero probabilities and spell diagnostics inherit this limitation. Unquantified trace is not given an invented positive amount.

## Conditional marginal distributions

Four meteorological seasons (DJF, MAM, JJA, SON) provide enough samples without fitting a separate distribution to every hour of every month. An inadequately sampled season is pooled across all training seasons and explicitly flagged. A variable with zero or undefined training variance is rejected rather than standardized.

| Variable | Modeled quantity | Candidate distributions |
| --- | --- | --- |
| Precipitation | Point mass at quantified zero, then positive hourly amount | Gamma, lognormal, empirical |
| Temperature | Anomaly from a smooth training-only seasonal and local-hour harmonic mean | Normal, empirical |
| Radiation | Daytime radiation divided by a geometric solar-elevation factor; daytime zeros have a mass | Gamma, lognormal, empirical |
| Wind | Calm-zero mass, then positive speed | Weibull, gamma, empirical |

Temperature has an intercept, annual and semiannual harmonics, diurnal and semidiurnal harmonics, and two annual–diurnal interaction terms. It remains in degrees Celsius, including negative temperatures. Negative and zero temperatures remain valid meteorological observations; no hydrologic applicability filter is applied to fitting.

Each candidate is fitted on the earliest 70% of the relevant chronological training sample and scored on the last 30%, using mean pinball loss at probabilities 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95 and 0.99. Positive-amount comparisons exclude zeros identically for every candidate. The lowest-loss candidate is then refitted on the full training subset. Empirical and parametric candidates use the same physical-unit loss; there is no invented empirical AIC. Saved candidates include inner-fit parameters or empirical sample count, validation counts, score and failure reason, with full final parameters or ordered empirical samples for the selected marginal. This small inner split is a pragmatic comparison, not exhaustive hyperparameter validation.

## Radiation, daytime and latent missing values

Approximate solar elevation follows [NOAA's published solar-position equations](https://gml.noaa.gov/grad/solcalc/solareqns.PDF), evaluated at twelve five-minute midpoints in each UTC hourly interval. The mean positive sine of elevation is the geometric factor. Its divisor is floored at 0.05 in twilight to avoid dividing by near-zero values; this floor is a numerical conditioning choice, not a calibrated clear-sky coefficient. Solar geometry is a deterministic astronomical covariate, not an alternative observed radiation product. The approximation uses the source's 365-day fractional-year convention and omits refraction; any discrepancy between geometry and observed twilight energy is recorded. The NASA representative coordinates are inherited from Stage B.

An hourly interval is called fully night if all twelve evaluation points are below the geometric horizon. Generated radiation is exactly zero in these bins. A partly sunlit hourly bin may have a positive hourly mean even when its midpoint is dark. Hourly means do not establish minute-level sunrise/sunset timing; the five-minute geometry is approximate. The retained hourly weather representation preserves the source interval-energy convention.

Nighttime radiation is **not** fitted as ordinary cloudiness. Its latent value is missing during dependence estimation. Source-missing daytime observations also remain missing. No random latent nighttime draws are inserted into the observations. The latent solar process is simulated continuously through the night to provide a state for the next day, but its nighttime values are unobserved and never converted to positive nighttime radiation.

## Gaussian copula and stable temporal model

For each observed training value, the selected seasonal marginal maps it to a CDF interval. The probability integral transform draws reproducibly inside point masses: precipitation zeros, calm winds and empirical ties. Applying the inverse standard-Normal CDF yields latent Gaussian values. `pit_seed` makes these draws reproducible. A randomized zero interval can weaken inferred wet/dry persistence: this is a known model limitation to be tested, not evidence that hourly precipitation can be treated independently.

The saved contemporaneous latent correlation is estimated after marginal transformation, not from raw physical weather. The code centers and scales the training latent values once. Pairwise counts and overlap limitations matter because nighttime solar is absent. Original and adjusted matrices, eigenvalues and maximum/Frobenius adjustments are saved. A symmetric eigenvalue floor, followed by diagonal normalization for a correlation matrix, repairs an invalid matrix if necessary; no identity substitution is used.

A first-order covariance-moment VAR is fitted as

`A = C1 @ inv(C0)`

where each entry of `C1` uses only genuinely adjacent, valid observed variable pairs. This pairwise approach allows nighttime precipitation/temperature/wind observations to contribute without pretending nighttime cloudiness was observed. It is a pragmatic missing-data estimator, not a full maximum-likelihood state-space fit. Different matrix entries may use different samples; full pair counts are saved. Residual lag-one diagnostics are separately calculated on complete observed daylight transitions, with their smaller sample size reported.

If the spectral radius of the fitted `A` exceeds the preregistered 0.98 maximum, it is reduced by uniform scaling and both matrices/radii are retained. The innovation covariance starts from

`Q = C0 - A @ C0 @ A.T`.

Any required positive-definite repair of `Q` is separately quantified. Generation is `z[t] = A @ z[t-1] + chol(Q) @ epsilon[t]`. Cholesky is applied to **innovation covariance Q once**. The original physical correlations are not applied again. The stationary covariance solves `Sigma = A Sigma A.T + Q`, using [SciPy's discrete Lyapunov solver](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_discrete_lyapunov.html).

After 336 burn-in hours, latent values are divided by their **fitted stationary** standard deviations before the standard-Normal CDF and conditional marginal inverse CDF. No generated scenario is centered or scaled by its own realized mean/variance. Consequently, scenario-to-scenario differences in weather intensity remain. The fitted zero process mean is explicit. Burn-in, PIT seed, matrices, stationary moments, training hash and fitted-model hash are saved.

## Aligned joint block-bootstrap benchmark

The benchmark selects complete contiguous 48-hour training blocks. All four variables come from the same source rows; partially missing blocks are ineligible. Block starts match the target local hour and are within one calendar month of the target month, cyclically. Adjacent-month pooling is a declared seasonal sample-size tradeoff. Every selected source start, block length and candidate-pool count is recorded in generation metadata.

The preregistered 48-hour length covers two diurnal cycles. Its justification is checked against the actual training latent autocorrelations at 1, 6, 12, 24, 48 and 72 hours. The fit records the maximum absolute correlation at the block length and whether it is below 0.2. If this screen fails, remaining long persistence is a documented limitation; the block is not silently retuned after seeing validation. The held-out spell and storm checks determine whether this simple benchmark is adequate for the proposed use.

Blocks preserve within-block joint weather variation. They can introduce jumps, sever long storms/dry spells and join unrelated weather at boundaries. Full-block eligibility can also overrepresent periods with good station coverage. The source calendar and target calendar may have different daylight, so the only physical-variable modification is setting target fully-night solar to zero; altered hours and removed `Wh/m²` are recorded. No source precipitation values are adjusted. Generated component extremes stay within training samples, although joining blocks can create storm totals not present as complete historical events. This is a benchmark, not a model of unprecedented weather.

## Held-out validation and acceptance

A full 2024 generated calendar is compared with the unshuffled held-out year. For every diagnostic, generated values receive the **same variable-specific missingness mask** as observations. Thus storms break at the same missing-rain hours and correlations use the same calendar support. Full generated arrays are retained separately for physical constraints, and full-sample descriptive statistics are also saved. Counts and censoring rules are explicit.

The reports include physical-unit means, standard deviations, quantiles and marginal CDFs; quantified-zero frequency; positive rainfall tails; raw Pearson/Spearman correlation matrices; month/hour-demeaned physical anomaly correlation matrices; all four-by-four lagged correlation matrices; wet/zero-amount spells; storm totals, durations and peaks; and local monthly and hourly mean patterns. Descriptive anomaly means in validation are computed separately for observed and generated samples; they do not refit any generator parameter.

Storms are joined across fewer than six quantified-zero hours and separated by at least six. Duration runs from the first through last positive-amount hour. Missing hours terminate events, never count as dry. At data gaps and calendar boundaries some spells/storms are censored, and both samples use the same gap mask. Trace lower bounds can understate occurrence and storm totals.

Configuration fixes the following exploratory screens **before** generated results are inspected:

- Wet-hour frequency absolute error ≤ 0.03.
- Positive hourly precipitation 99th percentile relative error ≤ 0.50.
- Storm total, duration and peak 95th percentiles relative error ≤ 0.50.
- Wet- and quantified-zero spell 95th percentiles relative error ≤ 0.50.
- Physical lag-one correlation absolute error ≤ 0.25 for each variable.
- No nonfinite generated inputs, negative precipitation/radiation/wind, or nonzero fully-night radiation.

Missing required comparison statistics fail their checks. The diagnostic matrices and seasonal/diurnal plots are reported without inventing post hoc pass thresholds. An overall pass requires every declared screen. A single generated held-out year also contains Monte Carlo/interannual variation; these screens are an exploratory adequacy check, not calibration, validation of rare floods, or a guarantee of correct joint tails.

Failed storm/persistence checks remain failed in the weather-validation report. The Gaussian copula lacks asymptotic upper-tail dependence; fitting a heavy marginal does not fix joint extremes. Empirical marginals cannot extrapolate beyond training support. A successful statistical diagnostic does not establish rare-event accuracy. No weather criterion is loosened to manufacture success.

## Interface

```python
fit = fit_weather(train_df, config)
weather = generate_weather(
    fit, target_hourly_utc_index,
    rng=np.random.default_rng(independent_seed_sequence),
    method="copula_var",       # or "block_bootstrap"
    training_df=train_df,
    config=config,
)
validation = validate_weather(heldout_df, generated_heldout_df, config, outdir)
```

The fitting runner assigns reproducible validation streams and records fitted-model/data hashes. One generated held-out calendar per method is retained for weather-only diagnostic comparisons; no batch is sent to a hydrologic model. Generation diagnostics preserve bootstrap source blocks and nighttime radiation adjustments. Synthetic test fixtures are never used as fallback training observations.
