# Underwood Creek: observed September 2022 hydrograph

Retrieved and processed **USGS 04087088** observations for **2022-09-02 00:00 UTC through 2022-09-15 00:00 UTC**, including both endpoints. This updates the earlier data-source review's “not downloaded” status for these USGS series only. No model simulation, calibration, parameter assignment or optimization was performed.

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

A separately retrieved **September 15–18 UTC recession extension** contains 289 approved, complete 15-minute values per parameter. At September 18 00:00 UTC, discharge is **11.5 ft³/s = 0.325644 m³/s**, still **5.45 times** the September 9 starting discharge of 2.11 ft³/s = 0.0597486 m³/s. The extension remains an observed-tail check in a separate file; the six-day B.2 horizon is unchanged.

## Reproduction and retained files

Run from the project root:

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
