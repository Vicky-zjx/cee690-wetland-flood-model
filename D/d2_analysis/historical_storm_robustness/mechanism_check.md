# Illustrative paired hydrograph and storage check

This read-only analysis uses the already completed `mc_00000` parameter draw in all five historical event windows. It does not run an ODE, refit a model, change forcing, or choose a new parameter vector. It is one illustrative draw, not an estimate of event probabilities or a causal explanation of ensemble averages.

Run from the project root: `python3 -B D/d2_analysis/historical_storm_robustness/mechanism_check.py`. The script reads the original parameter-study pilot NPZ/JSON files for the September 2022 benchmark and the new compact ensemble NPZ/JSON files for the other four events. It writes only this note and [the detailed CSV](tables/illustrative_mechanism.csv). That CSV records source paths and SHA-256 hashes. No saved simulation result is modified.

Event manifest SHA-256: `66371c784b4a2f9f7ea1c23f142bd7626e505c3558b5b86b4d4b1009631d2115`. Analysis script SHA-256: `2a2dc5c0dcd7d715ad18a7f0f2ad77b0159fcfce286723d40fbcf583a5807606`.

## Primary full-window quantities

Peaks and maximum rises retain the complete 0–144 h window and the first maximum on the native 15 s output grid. A peak at the initial condition remains the primary peak. The physical peak delay is left unavailable if either maximum occurs at time zero; a negative time difference involving the initial maximum must not be interpreted as an advance of the storm response.

| Event | Wet peak (m) / time (h) | No-wetland peak (m) / time (h) | Peak attenuation, no minus wet (mm) | Maximum rise wet / no (mm) |
|---|---:|---:|---:|---:|
| benchmark_september_2022 | 97.580137 / 82.983 | 97.693139 / 80.787 | 113.002 | 230.137 / 343.139 |
| historical_20230813T0800Z | 97.350000 / 0.000 | 97.417090 / 40.742 | 67.090 | 0.000 / 67.090 |
| historical_20220704T2100Z | 97.350000 / 0.000 | 97.406624 / 34.000 | 56.624 | 0.000 / 56.624 |
| historical_20190825T0200Z | 97.350000 / 0.000 | 97.350000 / 0.000 | 0.000 | 0.000 / 0.000 |
| historical_20240804T1500Z | 97.350000 / 0.000 | 97.350000 / 0.000 | 0.000 | 0.000 / 0.000 |

## Secondary timing and recession description

For descriptive storm timing only, a secondary maximum is taken from the start of each event's largest quantified wet spell through 144 h. Spell definitions and dates are inherited unchanged from the weather-only event manifest: positive quantified rain hours are joined unless at least six quantified-zero hours separate them; trace-zero is not verified dry weather. These secondary quantities do not replace the full-window outcome, redefine a flood threshold, or revise any simulation.

The later comparison starts after both secondary river maxima. Its crossover time is the first saved node where wetland stage minus no-wetland stage exceeds 0.000001 m (a reporting tolerance, not a hydraulic or safety limit). It is resolved to 15 s and is not an interpolated crossing. A positive later difference is retained as a normal outcome. The largest later difference can occur while small rainfall inputs continue, so this is a hydrograph/recession comparison rather than an assertion that all rainfall has ceased.

| Event | Secondary peak attenuation (mm) | Wet peak delay (h) | First later wet > no (h) | Largest later wet − no (mm) | Full sum-wetland storage peak (m³) / time (h) |
|---|---:|---:|---:|---:|---:|
| benchmark_september_2022 | 113.002 | 2.196 | 85.533 | 115.476 | 13656.995 / 82.142 |
| historical_20230813T0800Z | 100.346 | 3.229 | 47.117 | 57.282 | 6610.341 / 42.829 |
| historical_20220704T2100Z | 90.036 | 3.004 | 40.167 | 56.496 | 6510.825 / 35.879 |
| historical_20190825T0200Z | 61.522 | 5.000 | 49.046 | 35.668 | 5500.000 / 0.000 |
| historical_20240804T1500Z | 77.031 | 3.112 | 55.029 | 36.941 | 5500.000 / 0.000 |

## Interpretation and limits

All ten source simulations are complete and classified safety feasible. Each wet/no pair has identical saved time arrays. The four added-event actual gate-opening arrays are bit-for-bit identical; the benchmark differs by at most 5.69e-15 m (floating-point roundoff). The script verifies agreement within 1e-12 m, purely as an output-comparison tolerance. Shared upland, soil and dam-pool states differ by at most 3.64e-12 m³. Thus differences in actual gate operations do not explain the within-pair hydrograph differences in this draw. Both models use the same physical parameter vector. The wetland case contains an additional 5,500 m³ of accepted initial wetland water, so initial recession and the contribution of that stored water must not be presented as a rain-only or initially empty-wetland effect.

The saved wetland storage peaks and later river-stage crossovers are consistent with storage and delayed routing: a lower full-window peak can coexist with a higher later stage. The added-event compact arrays do not contain the complete normal-drainage, overflow, groundwater-exchange and ET flux histories, so this check does not quantify how much of the difference is caused by each individual flux.

For August 2023, the storm-period secondary attenuation is 100.346 mm versus 113.002 mm in September 2022, while their primary full-window attenuations are 67.090 mm and 113.002 mm. The August wetland storm response remains below the initial river stage. Its derived maximum clock-hour rainfall (16.693 mm) is comparable with the benchmark (15.906 mm) despite much lower total rainfall (58.9 versus 145.3 mm). These observations illustrate why cumulative rainfall alone cannot determine a peak-stage difference. They do not establish why ensemble averages may be similar, isolate an individual weather driver, or generalize beyond this sample.

All four meteorological variables vary jointly across the selected events. The benchmark retains its original mass-preserving 120 s reconstruction, whereas added events use the preserved Stage C hourly reconstruction. The selected warm complete windows are not a random sample of storms. Initial-state dominance, incomplete earlier-event coverage, and differences in forcing reconstruction remain material limits. No flood threshold is introduced by this analysis.
