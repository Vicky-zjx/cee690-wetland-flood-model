# Retained Stage C fitting and diagnostic results

Scope: meteorological data preparation, existing statistical fitting and weather-only figures/tables. The original results below were produced on 28 October 2026 and are retained for the current HW C.2 boundary. This retained report describes the earlier weather-model work. The subsequent official FIT/mimo_srs application is recorded separately in [C.2 fit summary](c2_fit_summary.md); the existing results below remain unchanged.

## Actual inputs and fitted weather models

Actual downloads cover 2019–2024: NOAA ISD/global-hourly KMKE `72640014839` supplies precipitation, KMWC `72640594869` supplies temperature/wind, and NASA POWER `ALLSKY_SFC_SW_DWN` supplies solar radiation at 43.0547222° N, −88.0461111° E. NOAA station observations and POWER grid estimates have different spatial support. The [data-source report](data_sources.md) contains official links, QC rules, timestamp meanings and yearly completeness. Exact URLs, retrieval times and raw-response hashes are retained locally under `C/data/raw/`.

The standardized calendar has 52,608 hourly UTC bins: 43,824 training hours in 2019–2023 and 8,784 held-out hours in 2024. Missing observations remain missing. In 2024, precipitation is missing in 9.35% of hours and wind in 12.92%; December precipitation has only 27 valid hours. Observed/generated validation uses identical variable-specific availability masks. Complete calendar windows do not imply complete observations.

The fitted models are seasonal mixed-marginal Gaussian-copula VAR(1) and aligned multivariate 48-hour block bootstrap. Training-only chronological quantile-loss selection chose lognormal precipitation in DJF, empirical in MAM, and gamma in JJA/SON. Temperature anomalies use empirical distributions except normal in JJA. Wind and daytime normalized radiation selected empirical distributions in every season. Full parameters, sample counts, selection scores, empirical support and dependence matrices are in `C/data/processed/weather_model.json`; [weather methods](weather_method.md) explain the choices.

VAR spectral radius is **0.9724845824**, below one. Latent-correlation and innovation-covariance Frobenius repairs are approximately `1.17e-15` and `4.24e-16`. Stability does not imply adequate temporal reproduction: complete-daylight residual lag-one correlations are −0.0502, 0.4418, 0.5228 and 0.0157 for precipitation, temperature, radiation and wind. Temperature and radiation retain substantial residual structure.

Global seed is `6902026`; method and scenario IDs define independent `SeedSequence` streams. Randomized PIT seed is `1729`; VAR burn-in is 336 hours. Every scenario has saved metadata, training/fitted-model versions and forcing-file hash. Both ODEs read the same saved CSV. All 20 scenario pairs passed forcing-file and canonical-array equality checks, including those later excluded.

## Held-out weather validation

One independently generated 2024 calendar per method was compared with the chronological holdout. Original acceptance tolerances were retained.

| Diagnostic | Observed | Copula/VAR | Block bootstrap |
|---|---:|---:|---:|
| Quantified wet-hour frequency | 0.08678 | 0.07535; pass | 0.07648; pass |
| Positive hourly precipitation p99, mm/h | 11.1207 | 8.6089; pass | 11.1021; pass |
| Storm-total p95, mm | 25.5880 | 7.7262; **fail** | 30.3770; pass |
| Wet-spell p95, h | 12.0 | 2.0; **fail** | 11.4; pass |
| Quantified-dry-spell p95, h | 149.2 | 45.55; **fail** | 163.35; pass |
| Precipitation lag-one correlation | 0.49067 | 0.05055; **fail** | 0.60976; pass |
| Physical-constraint violations | — | 0; pass | 0; pass |

Copula/VAR failed four of 12 configured held-out screens. Storm-duration/peak-intensity and temperature/radiation/wind lag-one screens passed their declared tolerances. The randomized precipitation point-mass transform and simple latent VAR do not retain adequate rain/spell persistence. Descriptive plots also show a poorly reproduced wind diurnal cycle and understated mean daytime radiation; no post hoc pass threshold was invented for these plots.

Bootstrap passed the 12 held-out screens. However, maximum absolute training latent autocorrelation at the 48-hour block length is **0.49434**, above the prespecified 0.2 persistence screen; temperature alone is 0.30631. Material long persistence is truncated at block boundaries. The benchmark is useful diagnostically but is not qualified for reliable flood-probability estimation. A single validation realization does not quantify validation Monte Carlo uncertainty; neither method establishes annual or extreme-tail flood probabilities.

Exact checks and diagnostic figures are in [`C/tables/weather_validation/`](../tables/weather_validation/) and [`C/plots/weather_validation/`](../plots/weather_validation/). The `validate` command regenerates these copies and the detailed local reports.

![Copula VAR storm and spell validation](../plots/weather_validation/copula_var/storms_and_spells.png)

![Bootstrap storm and spell validation](../plots/weather_validation/block_bootstrap/storms_and_spells.png)

## Reproduction and assignment status

The retained entry point is `python3 -B -m C.run_fitting all`. It fits the same training data and reruns the same held-out weather diagnostics with the original seeds, methods and thresholds. It does not execute A/B models or generate downstream response statistics. The two generated held-out calendars are retained solely as weather-validation examples.

The weather-only entry point was executed successfully after scope reduction. All **21 retained tests passed**. The fitted model, both generated validation calendars, both check tables and all **12 PNG figures** reproduced the pre-cleanup files byte-for-byte. Original source responses, acquisition metadata, standardized observations, QC and chronological split records remained unchanged. Existing failed weather checks remain failed; they were not retuned.

The existing environment is `/opt/miniconda3/bin/python3` (Python 3.13.5), NumPy 2.2.6, SciPy 1.16.2, pandas 2.3.2, Matplotlib 3.10.6 and requests 2.32.4. No packages were installed or removed. The execution log is `C/results/logs/fitting.log`; retained-test output is `C/results/fitting_tests.log`.

The four weather/data implementation modules and their statistical methods are unchanged. Marginal/correlation/Cholesky/VAR/bootstrap components remain because they support the existing fitting diagnostics and the detailed C.1/C.2 specification has not been located. They are not described as the professor's FIT implementation.

The subsequent [official FIT application](c2_fit_summary.md) adds the project-specific predictor/response definition, actual official call, coefficient/scaling tables and external fit/residual figures. Its approved causal data alignment and results are separate from the weather-model outputs retained here. The original official implementation and example remain unchanged in `A/work/multivarious`.
