# Underwood event: actual data status and scope

**Updated2026-10-01: bounded observations have been downloaded and processed; G1–G6 are approved for the first exploratory case. A valid and safety-feasible replacement reference has run; the bounded official OPT search and refined verification are complete.** This replaces the earlier metadata-only assessment. The concrete study mapping and approved fixed assumptions are in [study_configuration.md](01_study_configuration.md); the one approval table is [approval_gaps.csv](../tables/approval_gaps.csv).

The core retrieval window is2022-09-02T00:00Z–2022-09-15T00:00Z. The six-day model event remains September9–15UTC. NOAA includes one bracketing report on either side. A separate USGS-only September15–18recession extract checks the observed tail without extending the model horizon.

## Downloaded and processed

- **USGS04087088:** discharge00060, gageheight00065 and NAVD88 elevation63160, each1,250 records in the main window, all agency flagA, no missing nominal15-minute timestamps. PeakQ28.203579m³/s occurs September12 at08:45 and08:50UTC; peak elevation212.122512mNAVD88 at08:50. Paired elevation minus gageheight verifies the2022datum offset680.14ft. The separate recession extract has289records per series; September18flow is still5.45times September9startingflow. These observations are contextual and excluded from model forcing/validation. [Detailed USGS findings](04_usgs_observations.md).
- **KMWC meteorology:** NOAA ISD72640594869, approximately6.1km from the station context point. Temperature, dewpoint, wind and station pressure cover the event;262event observations, maximum70-minute interpolation bracket. Accepted event-time means set fixedRH0.7084916611553665 and stationpressure98692.51171875Pa. Neither sea-level pressure nor altimeter setting is substituted for station pressure.
- **KMKE precipitation:** NOAA ISD72640014839, approximately16km away. The selected raw accumulation sequence has310one-hour intervals plus one3-hour trace interval. The six-day quantified depth is145.3mm, excluding trace amounts. The13-day quantified depth is145.6mm. Five trace intervals intersect the event. The approved exploratory mixed-resolution/trace-lower-bound representation resolves the input representation without calling missing rain dry.
- **NASA POWER shortwave:**312valid UTC hourly values, requested at43.0547222,-88.0461111; APIv2.10.2, sourceSYN1DEG, approximately1degree grid. The event sum is21166.81Wh/m². This is a satellite-derived gridded estimate, not a local radiometer. Label-as-start-of-hour remains a disclosed timing assumption because the response lacks explicit bounds. [Detailed weather findings](05_weather_observations.md).

Raw bytes, exact queries/payloads, response times, checksums, QC, rejected/overlapping rain records and UTC/SI CSVs remain under ignored`B/work/event_data/`. Retrieval and processing scripts are`B/work/fetch_usgs_event.py` and`B/work/fetch_weather_event.py`. No annual NOAA files or basin-wide bulk data were downloaded.

## Missing or unsuitable data

MMSD's tested2022historical query returned all17gauges as0.00 with2026last-timestamp fields and no interval completeness/QC. It is unsuitable as verified2022rainfall and was excluded; this does not establish that MMSD has no historical archive. An authoritative WS1229historical interval export would improve the spatial rainfall evidence.

KMWC rain is too sparse and overlapping for a complete event series. USGS04087088returned no00045precipitation series for the window. KMKE's strict one-hour rain has a September13 12:52–14:52UTC gap; the actual3-hour trace ending14:52 covers that interval at coarser resolution, with the overlapping one-hour trace removed. Numerical zero for trace is an explicit lower-bound approximation accepted for this exploratory case, not an observation of no rain.

No local radiation instrument, hydraulically connected shallow groundwater well, initial soil-moisture/storage record, outlet blockage inspection, surveyed facility levels or quantitative vegetation response record was obtained. NLDAS was not downloaded through an authenticated service; the accepted POWER alternative avoids inventing inaccessible radiation. LiDAR/SSURGO/WWI source metadata remain identified but are not used to claim a surveyed hypothetical facility, bathymetry, species coefficient or lumped soil parameter.

## Geographic and boundary interpretation

Station metadata specify18.1mi²=46.879km². The downloaded WBD HUC12 Underwood Creek polygon040400030404 has an agency attribute area50.72km² and downstreamHUC040400030405. The NLDI point agrees with structured USGS coordinates43.0547222,-88.0461111; legacy narrative latitude43°03′00″ differs from the structured43°03′17″. NLDI split-at-station basin requests failedHTTP502. A retrieved whole-catchments simplified polygon is approximately50.56km² by a spherical diagnostic and is not substituted for the station catchment. Exact request records are in`geography_retrieval.json`.

The approved exploratory model remains an **unsited hypothetical0.66km²installation**, with its own local receiving channel. The real County Grounds diversion/tunnel arrangement is not the student's two dam–gate–wetland branches. No polygon is fabricated for those facilities and no total gauge discharge is used asQ_ext. Model elevations stay in one assumed datum; observed212mNAVD88 values are never mixed with model97–101mlevels. A real-site or whole-basin application still needs a suitable boundary, topographic/hydraulic survey and parameterization.

## Prepared forcing and actual checks

`B/prepare_b2_data.py --event-forcing` now builds4,321finite120srows from the processed actual data and approved hypothetical constant boundaries. Its manifest records the first-case approval separately from provenance. Each rain accumulation and solar hourly energy is conserved under a declared nonnegative piecewise-linear reconstruction with120sboundary ramps. A dry interval stays zero; trace rows remain separately flagged in source/manifest metadata. The resulting rainfall total is145.3mm and solar total21166.81Wh/m². T/wind use QC-accepted bracketing observations without extrapolation. Finer ODE checks must reuse these same knots.

The NOAA and USGS offline processing/checksum checks passed; forcing integral/coverage checks passed. These are data preparation checks, not physical calibration or model validation. The approved data drive the actual replacement reference and completed bounded official OPT search; run metrics and termination are reported separately in [09_exploratory_results.md](09_exploratory_results.md). Data availability and suitability are summarized in[B/tables/data_availability.csv](../tables/data_availability.csv).

Station observations, POWER grid estimates, spatial transfers to the unsited hypothetical footprint, trace-zero rainfall lower bounds, and assumed conditions retain separate evidence labels. Initial soil storage is a fixed inherited assumption, not an observed September 2022 reconstruction. USGS observations remain contextual and are neither external inflow nor an outlet calibration target.

## Narrow physical support retained

FAO56Chapter3 supports the existing vapor-pressure relation and selected constants: lambda_v2.45e6J/kg, cp1013J/(kgK),epsilon0.622 and water-density1000kg/m³. The JMA/UCAR constant table supportsRv461.5J/(kgK); NIST supports standardgravity9.80665m/s². These references do not validate the reduced ET closure. FAO wetlandKc values cannot be imported directly as itsK_veg, and Manningn cannot be equated toR_w. DNR1001does not provide the proposed flood-containment limits; their limited WisconsinPond378freeboard analogy and explicit design assumptions are documented separately.

All IDs resolve through[source_catalog.json](../source_catalog.json). Only three new identifiers were needed for this step: the concrete proposal itself, the Wisconsin freeboard analogue and the actual POWEReventextract. Existing NOAA/USGS/MMSD/WBD records were updated to their actual retrieval/processing status.
