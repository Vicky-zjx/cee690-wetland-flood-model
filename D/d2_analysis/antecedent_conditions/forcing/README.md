# Antecedent weather preparation

Four joint-weather inputs prepend the immediately preceding 24 or 72 hours to the unchanged August 2019 and August 2024 events. Only existing cached Stage C observations and grid estimates are read. No weather is generated, scaled, downloaded or filled. Original event CSV rows and all QC flags remain byte-for-byte unchanged. Each additional row is copied exactly from the original standardized CSV.

## Inputs and time semantics

| Event | Prefix (h) | Prefix rain (mm) | Solar energy (Wh/m²) | Minimum T (°C) |
|---|---:|---:|---:|---:|
| historical_20190825T0200Z | 24 | 0.000 | 5494.37 | 11.75 |
| historical_20190825T0200Z | 72 | 0.000 | 16390.58 | 11.75 |
| historical_20240804T1500Z | 24 | 0.000 | 6741.29 | 20.00 |
| historical_20240804T1500Z | 72 | 18.500 | 17522.47 | 18.00 |

All requested prefixes contain complete finite joint weather, complete native rainfall coverage and positive temperature. The 2019 prefixes and 2024 24 h prefix contain zero quantified rain, with no trace flags in those prefixes; the 2024 72 h prefix contains 18.5 mm and six trace-affected UTC hours. These statements describe the cached station record transferred to the hypothetical installation, not spatially exact local precipitation.

Native one-hour precipitation intervals, their actual starts/ends, boundary-overlap allocation and trace flags are saved in each `_native_rain.csv`. Native reports end at minute 52, so the UTC-hour rates are derived overlap allocations rather than native clock-hour observations. Each `_native_met_brackets.csv` records actual accepted temperature/wind observation times bracketing the midpoint target and the interpolation span. The raw file and request hashes, station sources, POWER units and all reconstruction errors are retained in each metadata file. No new interpolation changes the standardized rows.

## Preserve the original event mapping

The provider maps the prefix and 144 h event separately with the unchanged original hourly-input function, then shifts the event time grid by the prefix length. Precipitation and solar are constant within each UTC interval; temperature and wind are linear between their hourly midpoints with a half-hour endpoint hold. A naive single concatenated mapping would alter the original event's first half-hour temperature/wind values and is therefore not used.

The deliberate mapping seam retains the independent endpoint holds. At the event start the 2019 temperature/wind changes are −1.5 °C/−1.65 m/s; the 2024 changes are +1.27631579 °C/0 m/s. These are explicit discontinuities of the retained representation, not measured instantaneous jumps. They occur at an exact hourly and 120 s controller boundary. The caller must continue the physical state, actual actuator opening and controller latch state across this seam without a reset.

All eight preparation checks (four inputs at 15 s and 7.5 s) conserve rainfall depth and solar energy and reproduce every original event weather/fixed-boundary input value bit-for-bit after the time shift. These are input-mapping checks; no ODE, safety, model convergence or antecedent-state validation is performed here. The first two rows of the original input mapper are historical gate placeholders; the continuous feedback runner replaces them with its unchanged policy, so event equality is checked on rows 2 onward.

## Provider contract

`record_for(event_id, hours)` returns one JSON-serializable manifest record. `read_forcing(path)` returns the checked combined frame and diagnostics, with seam metadata attached to `frame.attrs`. `hourly_input_segments(frame,fixed)` and `forcing_conservation(frame,fixed)` match the original adapter signatures. The caller privately extends the numerical duration to 168 h or 216 h and delegates all other model functions unchanged.

Run `/opt/miniconda3/bin/python3 -B D/d2_analysis/antecedent_conditions/forcing/prepare.py` from the project root to prepare a new, empty forcing destination. Existing completed manifests are not overwritten.

The 24/72 h conditioning periods are controlled experimental choices, not evidence of equilibration or observed initial soil/groundwater/wetland storage. All existing fixed boundaries, including humidity/pressure assumptions, remain unchanged. These prefixes do not alter the original 144 h event metrics window.
