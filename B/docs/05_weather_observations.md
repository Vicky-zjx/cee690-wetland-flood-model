# September 2022 weather: sources, processing and accepted use

**Current status: records were downloaded and processed, and the documented reconstruction was approved for the hypothetical exploratory installation. The completed fixed-gate wetland optimization and dynamic comparison use this same forcing.** See the [current B overview](../README.md), [accepted assumptions](../data_sources_and_assumptions.md), [optimization results](12_wetland_only_optimization.md) and [comparison results](13_no_wetland_counterfactual.md). The retrieval covers **2022-09-02 00:00 through 2022-09-15 00:00 UTC**, including antecedent conditions. The comparison event is **September 9–15 UTC**. NOAA requests include one additional hour on each side to obtain interpolation/accumulation brackets. No annual or regional bulk files were needed. A.2 code, configuration, equations and figures were not changed.

## Concrete available inputs

| Input | Actual source and coverage | Accepted use and limitation |
|---|---|---|
| Air temperature and wind | NOAA ISD **KMWC / 72640594869**, Timmerman Airport: 585 records in the requested window, including 567 FM-15 meteorological observations; 569 with the two boundary observations. The event contains 262 observations. Accepted temperature, dew point, wind and station pressure bracket the complete event and antecedent windows, with a maximum 70-minute gap. | Piecewise linear interpolation of the separate KMWC observations, with a 90-minute maximum permitted bracket. The airport is approximately 6.1 km from the USGS context point. Wind sensor height is not established by these responses; no height correction or site calibration was applied. |
| Rainfall | NOAA ISD **KMKE / 72640014839**, Mitchell Airport: 481 records in the requested window, including 312 hourly FM-15 reports at minute 52. Two reports lack an independent one-hour rain amount. The selected, nonoverlapping accumulation series contains **310 one-hour intervals and one three-hour trace interval**, including the end bracket. | Explicit spatial proxy approximately 16.0 km from the USGS context point, not observed Underwood areal rainfall. The quantified event depth is **145.3 mm**; the full 13-day depth is **145.6 mm**, excluding unquantified traces. The largest reported one-hour amount is **16.0 mm**, September 12 05:52–06:52 UTC. |
| Shortwave radiation | NASA POWER `ALLSKY_SFC_SW_DWN`, point **43.0547222, −88.0461111**, `time-standard=UTC`: **312/312 valid hourly values**, API v2.10.2, returned source `SYN1DEG`. | Approximately 1° satellite-derived grid estimate rather than a local pyranometer observation. Returned units are Wh/m² per hourly record; dividing by one hour gives numerically identical mean W/m². The event's sum is **21,166.81 Wh/m²**. |
| Fixed relative humidity | Derived at each accepted KMWC observation as `exp(17.27 Td/(Td+237.3) − 17.27 T/(T+237.3))`, matching the existing saturation-pressure relation. | Retain **`r_H = 0.7084916611553665`**: the event-time mean of linearly interpolated derived RH, not the observation-count median. It remains a fixed approximation. |
| Fixed pressure | KMWC `MA1` station pressure, separately retained from altimeter setting and sea-level pressure. | Retain **`p_0 = 98692.51171875 Pa`**, using the event-time mean. This is station pressure at the airport, not pressure inferred from the model's assumed elevation datum. |

The station metadata and observations come from the bounded [NOAA access service](https://www.ncei.noaa.gov/support/access-data-service-api-user-documentation). POWER's [hourly service](https://power.larc.nasa.gov/docs/services/api/temporal/hourly/) and [SYN1deg documentation](https://power.larc.nasa.gov/docs/methodology/energy-fluxes/syn1deg/) identify the solar product and its spatial scale. The exact SYN1deg edition is not included in the response and is not asserted here. RH uses the [FAO saturation-pressure relationship](https://www.fao.org/4/x0490e/x0490e07.htm); its fixed-event treatment is an approved modeling approximation for this case.

## Rainfall timing, traces and missingness

KMWC rainfall is unsuitable as a complete event series: only 43 FM-15 reports carry a one-hour accumulation, one has an unsupported accumulated-period condition, and many intervals overlap. Absent `AA1` fields remain missing. No sum of those overlapping reports is treated as event rainfall.

For KMKE, the strict one-hour series has a **two-hour gap, September 13 12:52–14:52 UTC**. The 14:52 report supplies an actual **three-hour trace for 11:52–14:52 UTC**, quality code 1. The accepted reconstruction retains this three-hour record and removes the overlapping one-hour trace ending at 12:52. Its hourly distribution is unknown. Across the selected 13-day series, **nine trace intervals** remain flagged; **five intersect the six-day model event**.

**Treating trace amounts as numerical zero is an approved quantified-lower-bound approximation for this exploratory case, not evidence of dry weather.** No trace upper bound is asserted. Source amounts and durations remain in the interval CSV. The optional hourly CSV allocates each source accumulation uniformly over its own interval, crops boundary intervals to the requested window and preserves quantified mass to floating-point precision. It contains 312 bins; 16 bins carry trace flags. The actual rain report at **September 15 00:52 UTC** supplies the final eight minutes of the requested window. Two missing independent hourly measurements remain disclosed even though a coarser trace accumulation covers their period.

The [NOAA format document, January 12, 2018](https://www.ncei.noaa.gov/pub/data/noaa/isd-format-document.pdf) establishes UTC report times, scaling and quality fields. Parsing rejects numeric missing sentinels and suspect/error quality codes. It accepts ordinary valid quality codes and AWOS whole-degree `C` temperatures/dew points. Rain condition codes are separate from quality codes: condition 9 means unspecified condition, not necessarily missing depth; condition 2 denotes trace. Original fields and all rejected accumulation rows are preserved. `MA1` distinguishes altimeter setting from station pressure.

## MMSD and solar access findings

The [MMSD rain page](https://www.mmsd.com/about-us/milwaukee-rain-facility-information) exposes a public historical date-range query. A read-only POST for September 2–14, 2022 returned HTTP 200 and the requested dates, but **all 17 gauges returned 0.00 with October 2026 `LastTimeStamp` fields**. WS1229 appears as `NewGaugeWS1229`, Wauwatosa, 43.037460, −88.065840. The response gives neither interval records nor usable historical QC/completeness evidence. It was **not accepted as 2022 rainfall**. The linked export page describes the last four days. This establishes a failed usable-event retrieval through the tested public interfaces, not proof that MMSD lacks a historical archive.

POWER supplied the bounded solar request without login. Its JSON identifies UTC hourly labels but contains no explicit time bounds. **Interpreting each label as the start of its corresponding hour is an approved, disclosed processing assumption**, not a verified timestamp anchor. Preserve the hourly energy when constructing a model input; interpolating rates can otherwise alter it. The [GES DISC NLDAS access guidance](https://hydro1.gesdisc.eosdis.nasa.gov/data/NLDAS/) reviewed during data acquisition required Earthdata registration and noted a cloud transition. No NLDAS granules, credentials or account configuration were used. A nearby SURFRAD series was not pursued after the bounded POWER response succeeded.

## Source of each weather variable

The links below are the **exact primary bounded requests recorded at retrieval**, not generic references attached to assumed numbers. All NOAA `DATE` report timestamps are UTC. Requests cover September 2–15, 2022; separate start/end bracket responses are retained locally. Only September 9–15 drives the model. A fresh API response may differ from the retained bytes, so the recorded hashes below identify what was actually used.

| Model quantity | Source, returned field and unit | Derivation and current role | Limitation |
|---|---|---|---|
| `P(t)` precipitation | [NOAA KMKE 72640014839 exact request](https://www.ncei.noaa.gov/access/services/data/v1?dataset=global-hourly&stations=72640014839&startDate=2022-09-02T00%3A00%3A00Z&endDate=2022-09-15T00%3A00%3A00Z&format=json&includeAttributes=true&includeStationName=true&includeStationLocation=true); `AA1`–`AA4` provide accumulation duration, depth, condition and quality; depth scaling is 0.1 mm | Select accepted nonoverlapping intervals, divide each depth by its duration, then preserve each represented depth in the 120 s piecewise-linear input; 145.3 mm quantified event total | Regional station transfer; two missing independent one-hour amounts, five event trace lower bounds and assumed within-interval timing |
| `T_a(t)` temperature | [NOAA KMWC 72640594869 exact request](https://www.ncei.noaa.gov/access/services/data/v1?dataset=global-hourly&stations=72640594869&startDate=2022-09-02T00%3A00%3A00Z&endDate=2022-09-15T00%3A00%3A00Z&format=json&includeAttributes=true&includeStationName=true&includeStationLocation=true); `TMP`, scaling 0.1 °C, with quality field | Interpolate accepted report temperatures in °C; the evaporation function separately forms kelvin temperature as `T_K=T_a+273.15`; 262 event observations, 12–28 °C observed range | Maximum accepted interpolation bracket 70 min; hypothetical site is not at the airport |
| `v(t)` wind speed | [Same KMWC exact request](https://www.ncei.noaa.gov/access/services/data/v1?dataset=global-hourly&stations=72640594869&startDate=2022-09-02T00%3A00%3A00Z&endDate=2022-09-15T00%3A00%3A00Z&format=json&includeAttributes=true&includeStationName=true&includeStationLocation=true); speed component of `WND`, scaling 0.1 m/s, with quality field | Interpolate accepted wind reports; 262 event observations, 0–12.4 m/s observed range | Instrument height/exposure adjustment is not supported by the retrieved response; no height correction applied |
| `r_H` fixed relative humidity | [Same KMWC exact request](https://www.ncei.noaa.gov/access/services/data/v1?dataset=global-hourly&stations=72640594869&startDate=2022-09-02T00%3A00%3A00Z&endDate=2022-09-15T00%3A00%3A00Z&format=json&includeAttributes=true&includeStationName=true&includeStationLocation=true); `TMP` and `DEW`, both scaled to °C | Compute `exp(17.27 Td/(Td+237.3) - 17.27 T/(T+237.3))` at accepted observations, linearly interpolate the derived ratios, and time-average over the six-day event: 0.7084916611553665 | Derived from regional observations; holding that event mean fixed is an assumption, not measured constant humidity |
| `p_0` fixed pressure | [Same KMWC exact request](https://www.ncei.noaa.gov/access/services/data/v1?dataset=global-hourly&stations=72640594869&startDate=2022-09-02T00%3A00%3A00Z&endDate=2022-09-15T00%3A00%3A00Z&format=json&includeAttributes=true&includeStationName=true&includeStationLocation=true); station-pressure component of `MA1`, scaled to Pa | Time-average accepted, linearly interpolated station pressures: 98692.51171875 Pa | Fixed airport station-pressure approximation; neither sea-level pressure nor altimeter setting is substituted |
| `R_s(t)` shortwave radiation | [NASA POWER exact point/hourly request](https://power.larc.nasa.gov/api/temporal/hourly/point?parameters=ALLSKY_SFC_SW_DWN&community=RE&longitude=-88.0461111&latitude=43.0547222&start=20220902&end=20220914&format=JSON&time-standard=UTC); `ALLSKY_SFC_SW_DWN`, Wh/m² for each hourly record | Divide hourly energy by one hour to obtain the mean W/m²; reconstruct a nonnegative input that preserves each hourly energy; event integral 21166.81 Wh/m² | Approximately 1° `SYN1DEG` estimate; hour-start labeling is an accepted assumption, not established by explicit response time bounds |

Temperature and wind retain time variation. Their event-time means, 17.794950810185185 °C and 3.5009403935185186 m/s, are descriptive statistics only, not the fixed inputs used in the ODE. The quality filtering and field scaling follow the [NOAA ISD format document](https://www.ncei.noaa.gov/pub/data/noaa/isd-format-document.pdf); the humidity calculation uses the [FAO saturation-pressure relationship](https://www.fao.org/4/x0490e/x0490e07.htm).

## Accepted reconstruction and integral checks

The prepared input contains **4,321 finite knots at 120 s**, from September 9 00:00 to September 15 00:00 UTC inclusive. For rainfall and shortwave radiation, endpoints use the minimum of adjacent interval-mean rates and interior plateaus are adjusted to preserve represented interval depth/energy with 120 s linear ramps. This is assumed subinterval timing, not additional observations. Temperature/wind use accepted bracketing observations, with no extrapolation. The numerical reconstruction is implemented in [prepare_b2_data.py](../prepare_b2_data.py).

The completed 120/60/30 s comparisons use the same prepared knots. Every comparison preserves **145.300000 mm precipitation** and **21166.810000 Wh/m² shortwave energy**; refinement interpolates the existing input instead of rebuilding a different storm. See the [actual timestep results](13_no_wetland_counterfactual.md#time-resolution-provenance-and-remaining-storage). Fixed `Q_ext=0 m³/s`, `h_d=h_g=97.15 m` and `B=0` are separately approved boundary assumptions; NOAA, POWER and USGS do not supply these four values.

## Retrieval and retained-byte provenance

| Primary response | Recorded retrieval UTC | SHA-256 |
|---|---|---|
| NOAA KMKE main bounded JSON | `2026-10-02T01:41:04.705190+00:00` | `5126d339c0c47dd832adb65f73a713b9937f66aec29d9a1a96128581fc9de730` |
| NOAA KMWC main bounded JSON | `2026-10-02T01:41:04.704898+00:00` | `b70ea2dd40c6123d9832adf6265aba3c6e74dbda74ad4565d1942b316b7092a6` |
| NASA POWER hourly JSON | `2026-10-02T01:41:05.270771+00:00` | `cc8e3fd43f59727e11c126d24e2724c1cab6e3ebf85a6f5b81e8d942a2d7be00` |

The corresponding response filenames and exact URLs are recorded in the local `.request.json` sidecars. No new observation requests were made to prepare this documentation. The processed inputs used by the accepted reconstruction are identified below; hashes refer to the complete retained files, including the antecedent window where applicable.

| Processed input | SHA-256 |
|---|---|
| `KMKE_rain_selected_intervals_utc.csv` | `db9c2c3315e301e22a466b8cfebb811cdc7534ae839f200322e162513f1d2215` |
| `KMWC_meteorology_observations_utc.csv` | `fadc5f2c256f9dc3256d5c92088d3151458d1c9a4c7170ee565707c947550ec1` |
| `NASA_POWER_solar_hourly_utc.csv` | `702a80a86bde204118dc6e7a44e412efc0c5298ec4317a8d1e5edb8fff5930fe` |
| `proposed_forcing_120s.csv` (six-day model input) | `883135c877699ec0354f7a4abdcaffe326953293f88fe1dc7d5d6d7824336394` |

These files remain local under `B/work/event_data/`; the public source links, processing description and hashes provide traceability without implying that ignored files can be downloaded from this repository. The prepared manifest records approval and processing choices separately from the unchanged historical acquisition reports.

## Local reproducibility and files

Everything downloaded is retained under `B/work/event_data/weather/`. Each raw response has a `.request.json` sidecar recording its exact URL, POST payload where applicable, retrieval time, HTTP status, byte count and SHA-256. `weather_file_manifest.json` hashes the retained raw and processed files; `weather_coverage_gap_report.json` and the two station reports record counts, gaps, quality treatment, statistics and the proposals as they stood at preparation time. Their historical `proposal_approved=False` and `not_model_used` fields have not been rewritten. Current approval is recorded separately in `B/work/event_data/proposed_forcing_manifest.json` and the frozen [B model configuration](../no_spillway_config.json). These old flags do not describe the present completed study.

Main processed files:

- `KMWC_meteorology_observations_utc.csv` and `KMKE_meteorology_observations_utc.csv`: separate station observations, units, raw fields and quality status; no station splicing.
- `KMWC_rain_accumulations_audit_utc.csv` and `KMKE_rain_accumulations_audit_utc.csv`: every returned AA1–AA4 accumulation, including rejected or overlapping records.
- `KMKE_rain_selected_intervals_utc.csv`: `interval_start_utc`, `interval_end_utc`, `duration_h`, `depth_mm_quantified`, `trace`, condition/quality and `proposal_approved=False`.
- `KMKE_rain_hourly_utc_proposal.csv`: optional conservative UTC-hour aggregation with trace and coverage flags; no missing-to-zero filling.
- `NASA_POWER_solar_hourly_utc.csv`: original UTC hour label, hourly energy, equivalent hourly mean irradiance and source-status flag.

The commands below were used with the retained local processing scripts in `B/work/`. Those scripts and cached data are intentionally absent from the GitHub upload, so these commands are **local reproduction instructions**, not a complete fresh-clone data setup. Run from the project root; the first two commands need network access and reuse successful cached requests:

```sh
python3 -m B.work.fetch_weather_event
python3 -m B.work.fetch_weather_event --supplement
python3 -m B.work.fetch_weather_event --process
```

The processing uses Python's standard library. Checks verified parseable JSON, unique meteorological report keys, nonoverlapping selected rain intervals, full target coverage at the disclosed mixed resolution, rain mass preservation and all expected solar hour labels. These are data-processing checks, not model calibration or physical validation. Spatial representativeness, the trace treatment, solar-hour interpretation and fixed RH/pressure approximations remain part of the interpretation of every completed run.
