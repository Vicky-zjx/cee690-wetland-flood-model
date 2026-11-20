# Historical rainfall-window selection and source audit

Five fixed 144-hour weather windows are prepared: the original September 2022 benchmark **exactly once**, plus four additional historical windows selected from the existing 2019–2024 Stage C records. Selection uses rainfall and declared data/model applicability screens only. It uses no ODE output, wetland benefit, parameter realization, safety margin or optimization result. This file describes event preparation, not completed hydraulic validation or a new flood-frequency estimate.

## Fixed selection rule

1. Reserve the original September 9 00:00–September 15 00:00, 2022 UTC benchmark and its already accepted 120-second forcing file.
2. Enumerate every possible 144-hour window of the contiguous hourly 2019–2024 source calendar. Its ranking rainfall core is fixed at hours 24–47: a 24-hour accumulation preceded by 24 hours of window context.
3. Require all 144 standardized P, T, solar and wind values to be finite, every rain hour to have complete native-interval coverage, and every standardized temperature to exceed 0°C. Existing documented interpolation and trace flags remain visible. These are completeness and liquid-water applicability screens, not selections based on modeled outcomes.
4. Rank candidates by the core's quantified 24-hour rain total, descending; resolve exact stored-numeric ties by the earliest UTC core start. All arithmetic uses the retained standardized values, which were originally saved to ten significant digits.
5. Greedily retain four eligible windows that do not overlap the reserved benchmark or an already selected window. No event-specific parameter tuning, weather scaling or change in initial water storage enters this rule.

The fixed window may contain more than one wet spell. Calling these **rainfall windows** avoids asserting that a six-day interval is one independently verified meteorological storm. The July 2022 window is a separate historical event; it is not a second execution or replacement of the September benchmark.

[The manifest](events/manifest.json) records five event IDs, timestamps, hashes, sources, original row locations and QC. [The complete compressed candidate inventory](events/candidate_inventory.csv.gz) records every candidate, including missingness, freezing, rain coverage, ranking core, full-window statistics, reserved-benchmark overlap and selected status. No screened candidate is erased from that inventory.

There are **52,465 candidate windows**, including 45,684 with a complete ranking core, 20,516 with complete standardized four-variable coverage, and 9,440 passing complete-input/nonfreezing eligibility. These overlapping-window counts are not counts of independent storms or probabilities.

## Selected windows

All times below are UTC; end times are exclusive. Rain totals are quantified amounts and may understate unresolved traces. “UTC-hour maximum” is a derived clock-bin depth, not the maximum original station-report amount. The benchmark UTC-hour numbers integrate its unchanged 120-second forcing; additional windows sum their unchanged hourly rows.

| Event ID | Start UTC | End UTC | 144h rain (mm) | Maximum UTC hour (mm) | Maximum 24h (mm) | Maximum native 1h observation (mm) |
|---|---|---|---:|---:|---:|---:|
| benchmark_september_2022 | 2022-09-09 00:00 | 2022-09-15 00:00 | 145.300 | 15.906 | 141.400 | 16.0 |
| historical_20230813T0800Z | 2023-08-13 08:00 | 2023-08-19 08:00 | 58.900 | 16.693 | 53.300 | 18.8 |
| historical_20220704T2100Z | 2022-07-04 21:00 | 2022-07-10 21:00 | 66.427 | 26.147 | 47.200 | 29.0 |
| historical_20190825T0200Z | 2019-08-25 02:00 | 2019-08-31 02:00 | 42.800 | 11.367 | 42.300 | 11.9 |
| historical_20240804T1500Z | 2024-08-04 15:00 | 2024-08-10 15:00 | 41.200 | 20.020 | 41.200 | 23.1 |

All four additional windows have smaller 24-hour totals than the 141.4 mm benchmark. This test therefore checks robustness across distinct observed heavy-rain windows; it does not claim to stress the model with events more severe than the benchmark in 24-hour depth. Several have larger **native one-hour** bursts than the benchmark's 16.0 mm; duration-specific rankings differ. The rule does not guarantee inclusion of each year's largest hourly burst.

[Selected-event summary](events/selected_events.csv) and individual event files are linked through the manifest. New weather files retain the exact original CSV header and 144 row byte strings, including the original precision and every quality flag. Reading the subsets back verifies exact numeric values; pandas may infer an all-ones coverage column as integer rather than the source table's float type, which changes neither values nor bytes.

## Native observations versus reconstructed hourly weather

- **Rainfall:** NOAA ISD/global-hourly KMKE, station 72640014839, independent FM-15 one-hour AA accumulations. The report time ends the observed interval. For these four additional windows all selected native intervals end at minute 52, e.g. 05:52–06:52, rather than on a clock hour. Accepted conditions, quality codes, duplicate and overlap treatment are the unchanged Stage C parser rules. The audit exports each intersecting native interval, observed depth, trace flag, actual duration, boundary overlap and attributed window depth.
- **UTC-bin rainfall:** Stage C conservatively overlaps the native intervals onto [t,t+1h), assuming a constant rate within each observed interval. This preserves represented depth but can spread a short burst or wet/dry transition across two clock bins. The unobserved subhour timing is not recovered. No long accumulation is disaggregated for these new windows, and no missing rainfall becomes zero.
- **Temperature and wind:** NOAA KMWC, station 72640594869, quality-screened reports interpolated to hourly midpoints. The existing maximum bracketing separation is two hours; brackets exceeding 75 minutes are separately flagged as short-gap interpolation. These are regional interpolated inputs, not 144 direct on-the-hour measurements. No new interpolation is added here.
- **Solar:** NASA POWER ALLSKY_SFC_SW_DWN at the retained Underwood context point, returned source SYN1DEG. This is a gridded estimate, not a station pyranometer. Hourly energy in Wh/m² is numerically the hourly mean W/m²; the inherited interpretation of each UTC label as an hour start remains an explicit processing assumption because the response provides no independent interval bounds.

The four variables stay in the same original rows; there is no independent resampling or shuffling. This is retrospective historical forcing, not an operational forecast and not the separately lagged FIT predictor dataset. Airport rainfall and meteorology and the solar grid have different spatial support; none is asserted to be measured areal weather at the hypothetical 0.66 km² installation.

The native rainfall reconstruction agrees with saved C precipitation to at most **3.34×10⁻⁹ mm per hour**, consistent with the existing ten-significant-digit serialization. This is a processing consistency check, not measurement accuracy. Each new window intersects 145 native one-hour intervals: 143 lie wholly inside it and two are cropped at its edges. All original hourly source rows and QC files remain unchanged.

Exact raw-file hashes, retrieval timestamps and official request URLs for all 18 cached station-year/solar-year responses are retained in the manifest; every raw response was checked against its request sidecar. The [existing source/QC report](../../../C/docs/data_sources.md) and [processing code](../../../C/uncertainty/data.py) explain the unchanged NOAA field and POWER unit rules. No data was downloaded for this extension.

### Retained flags in the additional windows

“Interpolated” includes ordinary interpolation to the hour midpoint. Short-gap counts are subsets of those interpolation counts. Trace-affected clock bins can contain positive amounts contributed by adjacent native intervals.

| Event ID | Trace-affected hours | T midpoint interpolation | Wind midpoint interpolation | T short-gap interpolation | Wind short-gap interpolation |
|---|---:|---:|---:|---:|---:|
| historical_20230813T0800Z | 15 | 144 | 144 | 0 | 0 |
| historical_20220704T2100Z | 8 | 143 | 143 | 0 | 0 |
| historical_20190825T0200Z | 8 | 143 | 143 | 2 | 2 |
| historical_20240804T1500Z | 9 | 144 | 144 | 0 | 0 |

The 2019 window includes two existing short-gap interpolations for each of temperature and wind. Those values are retained with their flags rather than presented as direct observations or silently re-filled.

## Coverage and excluded heavier candidates

| Year | Calendar hours | Valid rain | Valid T | Valid solar | Valid wind | Four-variable complete hours | Eligible 144h windows |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2019 | 8760 | 8484 | 8749 | 8760 | 8746 | 8470 | 1777 |
| 2020 | 8784 | 8549 | 8740 | 8784 | 8737 | 8502 | 1434 |
| 2021 | 8760 | 8632 | 8700 | 8760 | 8714 | 8574 | 2110 |
| 2022 | 8760 | 8580 | 8597 | 8760 | 8591 | 8414 | 1431 |
| 2023 | 8760 | 8685 | 8316 | 8760 | 8313 | 8233 | 1224 |
| 2024 | 8784 | 7963 | 8602 | 8784 | 7649 | 7076 | 1464 |

[Annual coverage](events/annual_coverage.csv) also records nonpositive-temperature hours and known annual hourly/24-hour rainfall maxima. Candidate-year counts refer to window start years; the final 143 hours cannot start a complete window inside the cached period. High completeness at individual hours does not imply a complete six-day window. Selecting only fully supported, warm windows introduces coverage and applicability bias, especially where 2024 wind or rainfall is missing. No absent hour is treated as dry and no winter precipitation is silently converted to modeled snowmelt.

Examples from the [highest distinct screening windows](events/highest_nonoverlap_screening_windows.csv):

| Ranking core date | Known core depth (mm/24h) | Why the corresponding fixed window cannot be used as-is |
|---|---:|---|
| 2020-05-17 | 76.967 | Two missing standardized rain hours |
| 2023-02-26 | 52.400 | Three missing T and wind hours, and temperature reaches −5°C |
| 2024-05-28 | 47.867 | Two missing standardized rain hours |
| 2019-05-22 | 44.120 | Two missing standardized rain hours |

A highest-ranked July 2022 candidate beginning July 4 08:00 has two missing rain hours; a later eligible window beginning July 4 21:00 retains the same 47.2 mm ranking-core depth. The declared rule ranks eligible windows, including fixed-window alternatives at tied rainfall totals; it does not invent the missing early rain. Both candidates remain in the complete inventory. The compact nonoverlapping screening display is only a readable subset of that inventory and is not a second selection algorithm.

## Rain-duration diagnostics

A quantified wet spell joins positive-rain hours unless at least six quantified-zero hours separate them. Trace zeros are unresolved lower bounds, so this is an exploratory operational description, not an official storm definition or evidence of meteorological dryness. The first-to-last positive span can include long dry gaps; it is not the wet-hour count.

| Event ID | Positive-rain hours | First-to-last positive span (h) | Quantified wet spells | Largest-depth spell span (h) | Its rain (mm) |
|---|---:|---:|---:|---:|---:|
| benchmark_september_2022 | 27 | 49 | 2 | 21 | 141.400 |
| historical_20230813T0800Z | 19 | 71 | 2 | 19 | 53.300 |
| historical_20220704T2100Z | 15 | 93 | 3 | 9 | 47.200 |
| historical_20190825T0200Z | 17 | 41 | 2 | 15 | 42.300 |
| historical_20240804T1500Z | 9 | 9 | 1 | 9 | 41.200 |

## Preserve the original September 2022 benchmark

The benchmark points directly to `B/work/event_data/proposed_forcing_120s.csv` (SHA-256 `883135c877699ec0354f7a4abdcaffe326953293f88fe1dc7d5d6d7824336394`). Its manifest flag tells the runner to use the existing deterministic path rather than the new hourly adapter. The event ID occurs once in the selection. No C-derived replacement event or second September ensemble is created.

The [benchmark comparison JSON](events/benchmark_C_comparison.json) and [hour-by-hour comparison](events/benchmark_C_comparison.csv) document the difference at the same September 9–15 interval:

- Original input integrates to **145.300 mm** rain and **21,166.810 Wh/m²** solar energy.
- C has three missing rain bins at September 13 12:00, 13:00 and 14:00 UTC. The sum **over its observed bins** is 145.300 mm; that is not a complete observed total and does not establish that the absent bins were dry.
- Original B processing uses an actual three-hour trace record spanning the two missing native one-hour rain reports. Its numerical trace-zero treatment was explicitly approved as a quantified lower bound. The strict C one-hour selection rejects that longer accumulation, so the three partially unsupported clock bins remain missing.
- Original B rain and solar are reconstructed with 120-second mass-preserving linear ramps. New events retain C's hourly interval amounts. On the 141 common rain hours, the maximum difference is **0.206351 mm** and the mean absolute difference is **0.008586 mm**. The original benchmark maximum clock-hour amount is **15.906351 mm**, versus **15.700000 mm** in C.
- Hourly solar energy agrees within **4.7×10⁻¹⁰ Wh/m²**. T and wind at hourly midpoints agree within **4.6×10⁻⁹ °C** and **4.4×10⁻¹⁰ m/s**. Agreement at these points does not make the entire between-point forcing identical.

Thus the material cross-event processing caveat is rainfall timing/trace coverage, with a different continuous-time representation for the benchmark. No rainfall correction, scaling or missing-value fill is applied to make the representations match. [The accepted B source report](../../../B/docs/05_weather_observations.md) records the original evidence and reconstruction separately.

## Scope and reproducibility

The same model initial storages are assumed for each event; these are not reconstructed historical antecedent states or a spin-up result. Other physical settings, including the original fixed relative humidity and pressure, are not re-estimated by this event-preparation script. The cached Stage C requests do not supply the dewpoint/pressure fields needed for event-specific reconstructions. Hydraulic feasibility, safety and numerical reliability must be assessed by the separate model runs, with all failures retained; they did not influence this event selection.

Run from the project root with the existing environment:

```sh
/opt/miniconda3/bin/python3 -B D/d2_analysis/historical_storm_robustness/event_selection.py
```

The selector reads only existing local sources and writes only its new `events/` artifacts. It refuses to overwrite an already completed manifest. The current delivered manifest and CSVs are the frozen selection; no new event can be substituted after viewing a model result. [Preservation verification](events/preservation_check.json) confirms all **3,130** protected pre-existing files remain unchanged. No ODE was run by event preparation, and no Git operation was performed.
