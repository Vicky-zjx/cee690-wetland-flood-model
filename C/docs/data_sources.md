# Stage C meteorological data and quality control

Stage C uses actual multi-year regional data. These records are distinct from the approved September 2022 deterministic forcing in Stage B, which is unchanged. The initial fitting split is **2019–2023 training and 2024 chronological validation**. Six complete calendar-year windows were requested; this does not mean every station hour was observed. Missing data remain in the standardized table.

## Sources and geographic support

| Quantity | Product and verified identifier | Native information | Use and limitation |
|---|---|---|---|
| Precipitation | NOAA ISD / `global-hourly`, USAF `726400`, WBAN `14839`; concatenated identifier `72640014839`, call sign KMKE; General Mitchell International Airport | Routine METAR (`FM-15`) fields `AA1`–`AA4`, accumulation duration, depth in tenths of mm, condition and quality codes | Use only accepted independent one-hour accumulations. This regional airport proxy is not observed areal rainfall at the hypothetical installation. |
| Air temperature | NOAA ISD `72640594869`, USAF `726405`, WBAN `94869`; call sign KMWC, Lawrence J. Timmerman Airport | `TMP` in tenths of degrees Celsius, quality flag, UTC observation time | Station observations interpolated to each UTC hour's midpoint; no temperature clipping at freezing. |
| Wind | Same KMWC product and identifier | Speed component of `WND` in tenths of m/s with a separate quality flag | Regional reported wind, including accepted calm zeros. No verified sensor-height/exposure adjustment or conversion to a locally measured wetland wind is asserted. |
| Solar radiation | NASA POWER hourly point API, `ALLSKY_SFC_SW_DWN`, requested **43.0547222° N, −88.0461111° E**, `time-standard=UTC` | Returned hourly energy `Wh/m^2`; response identifies `SYN1DEG` | Gridded satellite-derived estimate, not station observation. Energy per hour becomes the numerically identical hourly mean W/m². The location is the existing Underwood context point, not a surveyed location of the hypothetical facility. |

NOAA's [ISD station history](https://www.ncei.noaa.gov/pub/data/noaa/isd-history.csv) was downloaded and checked by the USAF/WBAN identifiers, rather than assuming airport codes identify every NOAA product. It lists KMKE at 42.955, −87.905 and KMWC at 43.109, −88.031; more precise locations appear in each response. The history record supports the requested years. Airport stations represent regional conditions at different sites, while the solar grid has much broader support. The combined table is a regional modeling input, not four colocated instruments.

The exact successful NOAA requests use the [official access service](https://www.ncei.noaa.gov/support/access-data-service-api-user-documentation), one station/calendar year at a time. Field definitions and flags follow the [NOAA ISD format document](https://www.ncei.noaa.gov/pub/data/noaa/isd-format-document.pdf), especially pages 4–5, 7–10 and 12–13. POWER's [hourly documentation](https://power.larc.nasa.gov/docs/services/api/temporal/hourly/) identifies hourly averages and explicitly distinguishes UTC from Local Solar Time. The [SYN1deg product documentation](https://power.larc.nasa.gov/docs/methodology/energy-fluxes/syn1deg/) describes the solar product; no unretrieved edition number is asserted.

## Reproducible acquisition and saved artifacts

Run the data-preparation command from the Stage C README. The implementation is [`uncertainty/data.py`](../uncertainty/data.py), with a public `prepare_data(config, root)` function; `root` is the `C` directory. It uses `data` settings from [`configs/fitting.json`](../configs/fitting.json).

- `data/raw/isd-history.csv`: actual station metadata.
- `data/raw/NOAA_<station>_<year>.csv`: twelve actual station responses for the six years.
- `data/raw/POWER_solar_<year>.json`: six actual solar responses.
- Each raw file has a `.request.json` sidecar: exact request URL, returned URL, UTC retrieval timestamp, status, byte count and SHA-256.
- `data/processed/weather_hourly.csv`: 52,608 hourly timestamps, four weather variables, precipitation coverage and trace flags, interpolation/gap-fill flags, and chronological split.
- `data/processed/data_quality_report.json`: source manifest, output hash, per-year counts, QC exclusions, imputation fractions, remaining gap lengths and all conversion rules.
- `data/raw/retrieval_issues.json`: initial retrieval problems and their resolved alternatives. Failed future downloads generate `data/processed/download_failures.json`; there is no synthetic fallback.

The standardized CSV SHA-256 for the initial acquisition is `7edf24675b891d680863079323efb5da46e3d0c89a0640ed5852378a4be58fc4`. A subsequent intentional data refresh must retain its own hashes rather than treating changed bytes as the same dataset. Successful caches are reused only when their URL and content hash match the sidecar. Raw observations are local data, not copied into historical Stage A/B folders.

The annual static NOAA URL first returned HTML rather than data. Initial service requests spanning more than one calendar year unexpectedly returned only the trailing day for most years. Both were detected; the accepted acquisition uses separate January 1–December 31 requests and verifies all twelve months are represented. No fitting used those incomplete initial responses. Temperature/wind and precipitation edge hours without a same-year supporting observation remain missing rather than being extrapolated across the split.

## Time, accumulation and quality rules

All internal timestamps are timezone-aware UTC; Milwaukee seasonal and diurnal summaries use `America/Chicago`, including daylight-saving transitions. A standardized `timestamp_utc` labels **the beginning of `[t,t+1 hour)`**. Solar hour labels are interpreted as starts of the corresponding hourly energy intervals, preserving the explicit Stage B processing assumption: the JSON provides no independent interval bounds verifying that anchor. POWER Local Solar Time is never substituted for civil local time.

For precipitation, the actual NOAA report timestamp ends its accumulation. Accepted quality codes are `0`, `1`, `4`, `5`, `9`; missing sentinels and unsupported/suspect codes are rejected. Only independent one-hour routine reports are used. Identical duplicate intervals are deduplicated, conflicting duplicates invalidate the interval, and overlapping intervals are excluded with counts. Special reports and longer accumulations are not added to the selected one-hour series. This avoids both double counting and treating a six-hour total as an hourly intensity.

The selected one-hour intervals are conservatively overlap-weighted onto UTC hours, assuming a constant rate within each source accumulation. For example, a report ending at 00:52 covers 23:52–00:52; its contribution to either UTC bin is proportional to overlap. This preserves represented depth but can spread short wet/dry transitions between two bins, because unobserved subhour rainfall timing is unknown. An hour lacking full selected coverage is `NaN`, even if some depth was observed. No daily totals are disaggregated and no absent rain field becomes zero.

Trace condition `2` with zero encoded depth is retained as a **flagged quantified lower bound**, not a known exact dry hour. The generator models quantified precipitation; trace amounts remain unresolved. Zero-frequency and storm diagnostics therefore concern measured/quantified amounts and must not be described as exact meteorological occurrence probabilities. Trace-affected bin flags may coexist with positive precipitation contributed by an adjacent interval. Condition `9` alone does not invalidate an otherwise accepted numeric depth. Other conditions are excluded.

Temperature and wind use accepted routine and special observations, interpolated at the standardized hour's midpoint. Maximum permitted bracketing separation is two hours, with no extrapolation. Brackets above 75 minutes are counted separately as short-gap filling; ordinary interpolation to the midpoint is also explicitly counted. Conflicting duplicate meteorological observations are excluded. Temperature quality `C` is accepted for whole-degree AWOS reporting. Numeric sanity checks reject temperatures outside the ISD format range, negative or >90 m/s wind, and negative or >1500 W/m² solar. These are input-QC bounds, not hydrologic outcome selection. Longer gaps remain missing; no constant-value replacement is used.

POWER energy is divided by one hour to produce mean irradiance, so `sum(solar_W_per_m2)*1 hour` retains the reported energy. Nighttime physical constraints in generation are separate from raw source QC. Missing or deterministic-zero solar observations are not treated as random cloudiness samples by the fitted generator.

## Actual completeness of the initial acquisition

| Year | Expected hours | Valid precipitation | Valid temperature | Valid solar | Valid wind | All four valid |
|---|---:|---:|---:|---:|---:|---:|
| 2019 | 8,760 | 8,484 | 8,749 | 8,760 | 8,746 | 8,470 |
| 2020 | 8,784 | 8,549 | 8,740 | 8,784 | 8,737 | 8,502 |
| 2021 | 8,760 | 8,632 | 8,700 | 8,760 | 8,714 | 8,574 |
| 2022 | 8,760 | 8,580 | 8,597 | 8,760 | 8,591 | 8,414 |
| 2023 | 8,760 | 8,685 | 8,316 | 8,760 | 8,313 | 8,233 |
| 2024, holdout | 8,784 | 7,963 | 8,602 | 8,784 | 7,649 | 7,076 |

The 2024 holdout has **9.35% missing precipitation and 12.92% missing wind**. December precipitation is particularly incomplete (27/744 valid UTC hours); wind is reduced in September–December. September 2024 has 701/720 valid precipitation hours, 720/720 temperature and solar hours, and 665/720 wind hours. These limitations must remain visible in validation sample counts. Annual complete-case validation does not represent the unobserved hours, and nonconsecutive records must never be joined as adjacent VAR observations or bootstrap blocks.

All twelve monthly windows are retained, including winter measurements; the ODE analysis is a separate nonfreezing/rainfall applicability domain. The observed winter precipitation record must not be interpreted as modeled snowmelt. The weather analysis does not establish annual flood probability, local calibration, or complete extreme-event coverage.
