# September 2022 weather: downloaded evidence and forcing proposals

**Status: actual records retrieved and processed; forcing choices remain unapproved and have not been used in a model run.** The retrieval covers **2022-09-02 00:00 through 2022-09-15 00:00 UTC**, including antecedent conditions. The comparison event is **September 9–15 UTC**. NOAA requests include one additional hour on each side to obtain interpolation/accumulation brackets. No annual or regional bulk files were needed. A.2 code, configuration, equations and figures were not changed.

## Concrete available inputs

| Input | Actual source and coverage | Proposed use and limitation |
|---|---|---|
| Air temperature and wind | NOAA ISD **KMWC / 72640594869**, Timmerman Airport: 585 records in the requested window, including 567 FM-15 meteorological observations; 569 with the two boundary observations. The event contains 262 observations. Accepted temperature, dew point, wind and station pressure bracket the complete event and antecedent windows, with a maximum 70-minute gap. | Piecewise linear interpolation of the separate KMWC observations, with a 90-minute maximum permitted bracket. The airport is approximately 6.1 km from the proposed outlet point. Wind sensor height is not established by these responses; no height correction or site calibration was applied. |
| Rainfall | NOAA ISD **KMKE / 72640014839**, Mitchell Airport: 481 records in the requested window, including 312 hourly FM-15 reports at minute 52. Two reports lack an independent one-hour rain amount. The selected, nonoverlapping accumulation series contains **310 one-hour intervals and one three-hour trace interval**, including the end bracket. | Explicit spatial proxy approximately 16.0 km from the proposed outlet, not observed Underwood areal rainfall. The quantified event depth is **145.3 mm**; the full 13-day depth is **145.6 mm**, excluding unquantified traces. The largest reported one-hour amount is **16.0 mm**, September 12 05:52–06:52 UTC. |
| Shortwave radiation | NASA POWER `ALLSKY_SFC_SW_DWN`, point **43.0547222, −88.0461111**, `time-standard=UTC`: **312/312 valid hourly values**, API v2.10.2, returned source `SYN1DEG`. | Approximately 1° satellite-derived grid estimate rather than a local pyranometer observation. Returned units are Wh/m² per hourly record; dividing by one hour gives numerically identical mean W/m². The event's sum is **21,166.81 Wh/m²**. |
| Fixed relative humidity | Derived at each accepted KMWC observation as `exp(17.27 Td/(Td+237.3) − 17.27 T/(T+237.3))`, matching the existing saturation-pressure relation. | Propose **`r_H = 0.7084916612`**: the event-time mean of linearly interpolated derived RH, not the observation-count median. It remains a fixed approximation. |
| Fixed pressure | KMWC `MA1` station pressure, separately retained from altimeter setting and sea-level pressure. | Propose **`p_0 = 98692.51172 Pa`**, using the event-time mean. This is station pressure at the airport, not pressure inferred from the model's assumed elevation datum. |

The station metadata and observations come from the bounded [NOAA access service](https://www.ncei.noaa.gov/support/access-data-service-api-user-documentation). POWER's [hourly service](https://power.larc.nasa.gov/docs/services/api/temporal/hourly/) and [SYN1deg documentation](https://power.larc.nasa.gov/docs/methodology/energy-fluxes/syn1deg/) identify the solar product and its spatial scale. The exact SYN1deg edition is not included in the response and is not asserted here. RH uses the [FAO saturation-pressure relationship](https://www.fao.org/4/x0490e/x0490e07.htm); its fixed-event treatment is a modeling proposal.

## Rainfall timing, traces and missingness

KMWC rainfall is unsuitable as a complete event series: only 43 FM-15 reports carry a one-hour accumulation, one has an unsupported accumulated-period condition, and many intervals overlap. Absent `AA1` fields remain missing. No sum of those overlapping reports is treated as event rainfall.

For KMKE, the strict one-hour series has a **two-hour gap, September 13 12:52–14:52 UTC**. The 14:52 report supplies an actual **three-hour trace for 11:52–14:52 UTC**, quality code 1. The selected proposal retains this three-hour record and removes the overlapping one-hour trace ending at 12:52. Its hourly distribution is unknown. Across the selected series, **nine trace intervals** remain flagged.

**Treating trace amounts as numerical zero is an unapproved quantified-lower-bound choice, not evidence of dry weather.** No trace upper bound is asserted. Source amounts and durations remain in the interval CSV. The optional hourly CSV allocates each source accumulation uniformly over its own interval, crops boundary intervals to the requested window and preserves quantified mass to floating-point precision. It contains 312 bins; 16 bins carry trace flags. The actual rain report at **September 15 00:52 UTC** supplies the final eight minutes of the requested window. Two missing independent hourly measurements remain disclosed even though a coarser trace accumulation covers their period.

The [NOAA format document, January 12, 2018](https://www.ncei.noaa.gov/pub/data/noaa/isd-format-document.pdf) establishes UTC report times, scaling and quality fields. Parsing rejects numeric missing sentinels and suspect/error quality codes. It accepts ordinary valid quality codes and AWOS whole-degree `C` temperatures/dew points. Rain condition codes are separate from quality codes: condition 9 means unspecified condition, not necessarily missing depth; condition 2 denotes trace. Original fields and all rejected accumulation rows are preserved. `MA1` distinguishes altimeter setting from station pressure.

## MMSD and solar access findings

The [MMSD rain page](https://www.mmsd.com/about-us/milwaukee-rain-facility-information) exposes a public historical date-range query. A read-only POST for September 2–14, 2022 returned HTTP 200 and the requested dates, but **all 17 gauges returned 0.00 with October 2026 `LastTimeStamp` fields**. WS1229 appears as `NewGaugeWS1229`, Wauwatosa, 43.037460, −88.065840. The response gives neither interval records nor usable historical QC/completeness evidence. It was **not accepted as 2022 rainfall**. The linked export page describes the last four days. This establishes a failed usable-event retrieval through the tested public interfaces, not proof that MMSD lacks a historical archive.

POWER supplied the bounded solar request without login. Its JSON identifies UTC hourly labels but contains no explicit time bounds. **Interpreting each label as the start of its corresponding hour is a disclosed processing proposal**, not a verified timestamp anchor. Preserve the hourly energy when constructing a model input; interpolating rates can otherwise alter it. Current [GES DISC NLDAS access guidance](https://hydro1.gesdisc.eosdis.nasa.gov/data/NLDAS/) requires Earthdata registration and notes the ongoing cloud transition. No NLDAS granules, credentials or account configuration were used. A nearby SURFRAD series was not pursued after the bounded POWER response succeeded.

## Local reproducibility and files

Everything downloaded is retained under `B/work/event_data/weather/`. Each raw response has a `.request.json` sidecar recording its exact URL, POST payload where applicable, retrieval time, HTTP status, byte count and SHA-256. `weather_file_manifest.json` hashes the retained raw and processed files; `weather_coverage_gap_report.json` and the two station reports record counts, gaps, quality treatment, statistics and unapproved proposals.

Main processed files:

- `KMWC_meteorology_observations_utc.csv` and `KMKE_meteorology_observations_utc.csv`: separate station observations, units, raw fields and quality status; no station splicing.
- `KMWC_rain_accumulations_audit_utc.csv` and `KMKE_rain_accumulations_audit_utc.csv`: every returned AA1–AA4 accumulation, including rejected or overlapping records.
- `KMKE_rain_selected_intervals_utc.csv`: `interval_start_utc`, `interval_end_utc`, `duration_h`, `depth_mm_quantified`, `trace`, condition/quality and `proposal_approved=False`.
- `KMKE_rain_hourly_utc_proposal.csv`: optional conservative UTC-hour aggregation with trace and coverage flags; no missing-to-zero filling.
- `NASA_POWER_solar_hourly_utc.csv`: original UTC hour label, hourly energy, equivalent hourly mean irradiance and source-status flag.

Run from the project root; the first two commands need network access and reuse successful cached requests:

```sh
python3 -m B.work.fetch_weather_event
python3 -m B.work.fetch_weather_event --supplement
python3 -m B.work.fetch_weather_event --process
```

The processing uses Python's standard library. Checks verified parseable JSON, unique meteorological report keys, nonoverlapping selected rain intervals, full target coverage at the disclosed mixed resolution, rain mass preservation and all expected solar hour labels. These are data-processing checks, not model calibration or physical validation. Spatial representativeness, the trace treatment, solar-hour interpretation and fixed RH/pressure approximations must remain visible in any subsequent proposed run.
