# Underwood Creek: observed September 2022 hydrograph

Retrieved and processed **USGS 04087088** observations for **2022-09-02 00:00 UTC through 2022-09-15 00:00 UTC**, including both endpoints. These records provide geographic/event context only for the completed hypothetical B study. They are not used as an inflow, boundary level, calibration target or validation hydrograph. The [B overview](../README.md) and [accepted assumptions](../data_sources_and_assumptions.md) explain that distinction.

## Availability and observed peaks

The [official instantaneous-values query](https://waterservices.usgs.gov/nwis/iv/?format=json&sites=04087088&startDT=2022-09-02T00%3A00%3A00Z&endDT=2022-09-15T00%3A00%3A00Z&parameterCd=00060%2C00065%2C63160&siteStatus=all) returned **1,250 values for each of three parameters**. All carry agency qualifier **A**, meaning processing/review is complete and the values are approved for publication. The agency labels these as derived values; in particular, approved discharge is not a direct measurement at every timestamp.

| Parameter | Maximum in the six-day event window, September 9–15 UTC | UTC time on September 12 |
|---|---|---|
| 00060 discharge | 996 ft³/s = **28.2036 m³/s** | 08:45 **and** 08:50 |
| 00065 gage height | 15.80 ft = **4.81584 m above the gage datum** | 08:50 |
| 63160 water-surface elevation | 695.94 ft = **212.122512 m NAVD88** | 08:50 |

These are maxima of the returned observations, without interpolation. Each series contains every expected 15-minute timestamp plus an additional observation at September 12 08:50 UTC. There are **no missing values, duplicate timestamps or gaps exceeding 15 minutes**. The six-day window contains 578 values per parameter. Original timestamps retain the source's **CDT / −05:00** offset; processed timestamps are UTC.

A separate bounded [00045 precipitation query](https://waterservices.usgs.gov/nwis/iv/?format=json&sites=04087088&startDT=2022-09-02T00%3A00%3A00Z&endDT=2022-09-15T00%3A00%3A00Z&parameterCd=00045&siteStatus=all) returned HTTP 200 with **zero time series**. This is unavailable rainfall for the requested window, not zero rainfall.

## Datum, location and suitable use

The [station history](https://waterdata.usgs.gov/nwis/uv/?referred_module=sw&site_no=04087088) dates the 680.14 ft NAVD88 gage datum to October 3, 2016 onward. All 1,250 paired event values independently satisfy `63160 − 00065 = 680.14 ft` to floating-point precision. The SI file converts the supplied 63160 series directly; it does not mistake gage height for elevation.

The [expanded site metadata](https://waterservices.usgs.gov/nwis/site/?format=rdb&sites=04087088&siteOutput=expanded&siteStatus=all) supplies **43.0547222, −88.0461111 in NAD83**, corresponding to 43°03′17″ N, 88°02′46″ W. The legacy narrative instead says 43°03′00″ N, approximately 525 m south. Current structured metadata is a usable provisional point; the narrative discrepancy remains unresolved, so this is not a surveyed/snapped outlet. Preserve native NAD83 information even though the IV JSON labels its geographic coordinates EPSG:4326.

The station drains **18.1 mi² ≈ 46.879 km²**, about 71 times the A.2 0.66 km² footprint. These observations can characterize the surrounding event. They cannot serve as the small installation's measured output, as `Q_ext`, or as a receiving-channel boundary without a defensible mapping. Using observed outlet discharge as inflow while separately adding catchment rainfall would double count runoff. Neither these records nor the observed peak defines any dam/wetland safety limit, shallow-groundwater head or blockage factor.

Over September 9–15, trapezoidal integration on the native timestamps gives **2,300,664.88 m³** of passing discharge. This is a derived total-flow estimate, with no baseflow subtraction, not a simulated or excess-flood volume. Discharge remains 41.7 ft³/s at the final endpoint versus 2.11 ft³/s at the initial endpoint; this window captures the observed peak but not a complete return to the starting flow.

A separately retrieved [September 15–18 UTC recession extension](https://waterservices.usgs.gov/nwis/iv/?format=json&sites=04087088&startDT=2022-09-15T00%3A00%3A00Z&endDT=2022-09-18T00%3A00%3A00Z&parameterCd=00060%2C00065%2C63160&siteStatus=all) contains 289 approved, complete 15-minute values per parameter. At September 18 00:00 UTC, discharge is **11.5 ft³/s = 0.325644 m³/s**, still **5.45 times** the September 9 starting discharge of 2.11 ft³/s = 0.0597486 m³/s. The extension remains an observed-tail check in a separate file; the six-day B horizon is unchanged.

## Retrieval provenance and field-level derivations

Each hydrologic variable is supplied by the same exact main query above, with its own agency parameter code: **00060** discharge, **00065** gage height and **63160** water-surface elevation. The reported peak is a maximum of that parameter's unchanged timestamped values. Conversion factors are `0.028316846592 m³/s per ft³/s` for 00060 and `0.3048 m per ft` for 00065/63160. The native 63160 series supplies NAVD88 elevation directly; no assumed model elevation is added. The 680.14 ft datum offset is also independently recovered by paired subtraction. The reported event-flow volume uses trapezoidal integration of 00060 on its actual UTC timestamps; it is not a direct volume measurement.

The 18.1 mi² catchment area and NAD83 coordinates come from the exact expanded site-metadata query linked above. `18.1 × 2.589988110336 = 46.8787847970816 km²` is the international-unit area conversion. None of these geographic values locates or surveys the hypothetical 0.66 km² facilities.

| Retained response | Recorded retrieval UTC or explicitly labeled server time | SHA-256 |
|---|---|---|
| `usgs_04087088_iv_raw.json` | Client retrieval time not recorded; server response Date: Fri, 02 Oct 2026 01:39:55 GMT | `018be3469d87a86225e02cf8d33a4ea68035f68de3cb45d4a699562f0ace68d0` |
| `usgs_04087088_site_raw.rdb` | 2026-10-02T01:40:55.569535+00:00 | `73820f766ab15d4d3e3664991b757d3a20f7db8255db32945943d3316205f7f5` |
| `usgs_04087088_history_raw.html` | 2026-10-02T01:41:00.880475+00:00 | `cf3af69ce203c80fc07e5f8015345a5f91e6d8574e6192f86e4a6c173c0b031e` |
| `usgs_04087088_recession_raw.json` | 2026-10-02T01:46:12.458153+00:00 | `02f7e0821eef165d8643196aa4e6eeba6287fb1b8bcaf9b893061a7f95819b75` |
| `usgs_04087088_precip00045_raw.json` | 2026-10-02T01:43:46.027999+00:00 | `8eede38cf07fc81e20c73010e18618ad12e3ffdf98af2f05c840a767b9509107` |

The initial main-series sidecar did **not** record a client retrieval timestamp; its server response date is reported above rather than invented. The manifest records processing at `2026-10-02T01:47:07Z`. Retained processed main-series CSV SHA-256 is `3254488acc5a22e80a83192f5be4ebfc29aaed2234ed9a1b67c86a7f3ff6d7e5`; the separate recession CSV SHA-256 is `381b2beca025888a128a6db19c703e19c29921edba8831bf5ce11b80cc45f77d`. The queries, source values, quality flags and checksums identify the actual records consulted even if a future agency response changes. Raw and processed files remain local; the public URLs are agency queries, not links to ignored repository files.

## Reproduction and retained files

This command requires the retained local processing script and downloaded bytes under ignored `B/work/`; those files are not included in a fresh GitHub clone. Run from the project root:

```sh
/opt/miniconda3/bin/python3 -m B.work.fetch_usgs_event --offline
```

This tested command processes retained raw files without network access. Omitting `--offline` downloads any absent files through the exact recorded official queries; existing raw bytes are reused. The individual download requests and cached processing were executed successfully. The raw/processed files remain local in ignored `B/work/event_data/`:

- `usgs_04087088_20220902_20220915_utc_si.csv`: 1,250 unchanged observation times, SI values, original numeric values and qualifiers.
- `usgs_04087088_manifest.json`: exact queries, response-metadata locations, SHA-256 checksums, processing script hash and conversions.
- `usgs_04087088_recession_20220915_20220918_utc_si.csv` and `usgs_04087088_recession_availability.json`: separate 289-row observed-tail records and audit.
- `usgs_04087088_event_availability.json`: per-series/window counts, extrema, QC, gap checks and paired-datum verification.
- `usgs_04087088_agency_metadata.json`, raw JSON/RDB/HTML responses and HTTP metadata.

Conversions are `ft × 0.3048` and `ft³/s × 0.3048³`, using [conventional USGS hydrologic units](https://pubs.usgs.gov/circ/2003/circ1262/htdocs/conversions.html) and the [NIST international-foot definition](https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors). Source labels say ft, not US survey ft; raw feet are retained. The alternative survey-foot conversion would change the peak elevation by approximately 0.000424 m, below its 0.01 ft reporting increment. This does not authorize mixing terrain datums or foot conventions.

No filling, resampling, smoothing or rating-curve reconstruction was applied. The used legacy USGS service is scheduled for retirement in early 2027; cached files/checksums preserve this retrieval, while a future fresh download may require the replacement API.
