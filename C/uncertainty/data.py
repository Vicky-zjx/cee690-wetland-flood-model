"""Acquire real Milwaukee observations and conservatively align hourly weather.

The public API is prepare_data(config, root), where root is the Stage C folder.
No missing precipitation is converted to zero and no synthetic fallback exists.
"""
from __future__ import annotations

from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
import pandas as pd
import requests

VARIABLES = ["precipitation_mm_per_h", "temperature_C", "solar_W_per_m2", "wind_m_per_s"]
QC_OK = {"0", "1", "4", "5", "9"}
NOAA_URL = "https://www.ncei.noaa.gov/access/services/data/v1"
HISTORY_URL = "https://www.ncei.noaa.gov/pub/data/noaa/isd-history.csv"
POWER_URL = "https://power.larc.nasa.gov/api/temporal/hourly/point"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fetch_cached(path, url, download=True, timeout=120):
    """Cache successful, type-checked responses with exact URL and hash sidecars."""
    path = Path(path)
    meta_path = path.with_suffix(path.suffix + ".request.json")
    if path.exists() and meta_path.exists():
        meta = json.loads(meta_path.read_text())
        if meta.get("url") == url and meta.get("sha256") == sha256(path) and meta.get("ok"):
            return path
    if not download:
        raise FileNotFoundError(f"Required real data cache unavailable: {path}; enable data.download.")
    path.parent.mkdir(parents=True, exist_ok=True)
    meta = {"url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat(), "ok": False}
    try:
        response = requests.get(url, timeout=timeout, headers={"User-Agent": "CEE690-research-weather/1.0"})
        meta.update(http_status=response.status_code, final_url=response.url, bytes=len(response.content))
        response.raise_for_status()
        if response.text.lstrip().startswith("<"):
            raise ValueError("Expected observational CSV/JSON, received HTML; not accepted as data.")
        if path.suffix == ".json":
            obj = response.json()
            if "properties" not in obj:
                raise ValueError("POWER response has no properties; not an accepted data response.")
        path.write_bytes(response.content)
        meta.update(ok=True, sha256=sha256(path))
        print(f"Downloaded {path.name}: {len(response.content)} bytes", flush=True)
    except Exception as exc:
        meta["error"] = f"{type(exc).__name__}: {exc}"
        raise RuntimeError(f"Real observation download failed for {url}: {exc}") from exc
    finally:
        meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    return path


def parse_value(raw, value_index, quality_index, temperature=False):
    """Decode ISD tenths and reject missing/suspect/unsupported quality codes."""
    fields = str(raw).split(",")
    if len(fields) <= max(value_index, quality_index):
        return np.nan
    try:
        value = int(fields[value_index])
    except ValueError:
        return np.nan
    good = fields[quality_index] in QC_OK or (temperature and fields[quality_index] == "C")
    if value == 9999 or not good:
        return np.nan
    result = value / 10.0
    if temperature and not -93.2 <= result <= 61.8:
        return np.nan
    if not temperature and not 0 <= result <= 90.0:
        return np.nan
    return result


def precipitation_intervals(frame):
    """Retain QC-approved one-hour routine reports; never sum longer overlaps.

    Different report sources can repeat the same interval. Identical copies are
    deduplicated; conflicting amounts invalidate that interval. Greedy overlap
    removal is deterministic and reported. No hour-length reconstruction from
    daily or six-hour totals is made.
    """
    counts = Counter()
    candidates = {}
    routine = frame[frame["REPORT_TYPE"].str.strip() == "FM-15"]
    for row in routine.to_dict("records"):
        end = pd.to_datetime(row["DATE"], utc=True)
        for name in ("AA1", "AA2", "AA3", "AA4"):
            raw = str(row.get(name, ""))
            fields = raw.split(",")
            if len(fields) != 4:
                counts["absent_fields"] += 1
                continue
            try:
                hours, depth = int(fields[0]), int(fields[1])
            except ValueError:
                counts["invalid_numeric"] += 1
                continue
            condition, quality = fields[2:]
            if hours != 1:
                counts["not_one_hour"] += 1
                continue
            if depth == 9999 or depth < 0 or quality not in QC_OK or condition not in {"2", "9"}:
                counts["missing_bad_quality_or_unsupported_condition"] += 1
                continue
            trace = condition == "2"
            if trace and depth != 0:
                counts["inconsistent_trace_depth"] += 1
                continue
            value = (depth / 10., trace)
            key = (end - pd.Timedelta(hours=1), end)
            if key in candidates and value in candidates[key]:
                counts["identical_duplicate_intervals_removed"] += 1
            candidates.setdefault(key, set()).add(value)
    selected = []
    previous_end = None
    for (start, end), values in sorted(candidates.items()):
        if len(values) != 1:
            counts["conflicting_duplicate_intervals"] += 1
            continue
        if previous_end is not None and start < previous_end:
            counts["overlapping_intervals_excluded"] += 1
            continue
        depth, trace = next(iter(values))
        selected.append((start, end, depth, trace))
        previous_end = end
    counts["selected_intervals"] = len(selected)
    counts["trace_intervals"] = sum(int(r[3]) for r in selected)
    return selected, dict(counts)


def align_precipitation(intervals, index):
    """Overlap-weight actual nonoverlapping intervals, retaining incomplete NaN.

    Within each observed one-hour interval the rate is assumed constant. The
    target UTC-hour timestamp is its beginning. Trace depth is a flagged zero
    quantified lower bound, not a dry observation.
    """
    index = pd.DatetimeIndex(index)
    if index.tz is None:
        raise ValueError("Timezone-aware UTC index required")
    starts = index.as_unit("s").asi8.astype(float)
    depth = np.zeros(len(index)); coverage = np.zeros(len(index)); trace = np.zeros(len(index), bool)
    previous_end = None
    for lo, hi, amount, is_trace in intervals:
        lo, hi = pd.Timestamp(lo), pd.Timestamp(hi)
        if lo.tz is None or hi.tz is None or hi <= lo or amount < 0 or not np.isfinite(amount):
            raise ValueError("Invalid precipitation interval")
        if previous_end is not None and lo < previous_end:
            raise ValueError("Overlapping precipitation intervals are not allowed")
        previous_end = hi
        a, b = lo.timestamp(), hi.timestamp()
        left = max(0, np.searchsorted(starts, a, side="right") - 1)
        right = min(len(index), np.searchsorted(starts, b, side="left") + 1)
        overlap = np.maximum(0, np.minimum(b, starts[left:right] + 3600) - np.maximum(a, starts[left:right]))
        depth[left:right] += amount * overlap / (b-a)
        coverage[left:right] += overlap
        trace[left:right] |= (overlap > 0) & is_trace
    complete = np.isclose(coverage, 3600., atol=1e-6, rtol=0)
    depth[~complete] = np.nan
    return depth, coverage / 3600., trace


def interpolate_met(times, values, targets, max_gap_hours=2):
    """Interpolate only accepted bracketing observations, never extrapolate."""
    series = pd.Series(values, index=pd.DatetimeIndex(times), dtype=float).dropna().sort_index()
    # Conflicting simultaneous observations are missing rather than averaged.
    series = series.groupby(level=0).agg(lambda a: a.iloc[0] if a.nunique() == 1 else np.nan).dropna()
    x = series.index.as_unit("s").asi8.astype(float)
    y = series.to_numpy()
    target = pd.DatetimeIndex(targets).as_unit("s").asi8.astype(float)
    out = np.full(len(target), np.nan); interpolated = np.zeros(len(target), bool)
    short_gap = np.zeros(len(target), bool)
    if len(x) < 2:
        return out, interpolated, short_gap
    right = np.searchsorted(x, target, side="left")
    exact = (right < len(x)) & (x[np.minimum(right, len(x)-1)] == target)
    out[exact] = y[right[exact]]
    valid = (~exact) & (right > 0) & (right < len(x))
    ids = np.flatnonzero(valid)
    spans = x[right[ids]] - x[right[ids]-1]
    ids = ids[spans <= float(max_gap_hours)*3600]
    span = x[right[ids]] - x[right[ids]-1]
    fraction = (target[ids]-x[right[ids]-1])/span
    out[ids] = y[right[ids]-1] + fraction*(y[right[ids]]-y[right[ids]-1])
    interpolated[ids] = True
    short_gap[ids] = span > 75*60
    return out, interpolated, short_gap


def load_prepared(config, root):
    cfg = config.get("data", config)
    path = Path(root) / cfg.get("processed_path", "data/processed/weather_hourly.csv")
    frame = pd.read_csv(path)
    frame["timestamp_utc"] = pd.to_datetime(frame["timestamp_utc"], utc=True)
    return frame


def gap_statistics(values):
    """Count missing runs without bridging nonobserved hours."""
    missing = np.asarray(pd.isna(values), dtype=bool)
    edges = np.diff(np.r_[False, missing, False].astype(int))
    lengths = np.flatnonzero(edges == -1) - np.flatnonzero(edges == 1)
    return {"runs":len(lengths),"longest_hours":int(lengths.max()) if len(lengths) else 0,
            "hours_in_gaps_longer_than_2h":int(lengths[lengths>2].sum())}


def prepare_data(config, root):
    """Download/cache multiple complete calendar years and save hourly QC data."""
    cfg = config.get("data", config)
    root = Path(root)
    raw = root / cfg.get("raw_dir", "data/raw")
    output = root / cfg.get("processed_path", "data/processed/weather_hourly.csv")
    output.parent.mkdir(parents=True, exist_ok=True)
    train = [int(y) for y in cfg.get("train_years", [2019,2020,2021,2022,2023])]
    hold = [int(y) for y in cfg.get("validation_years", [2024])]
    if not train or not hold or max(train) >= min(hold) or set(train) & set(hold):
        raise ValueError("Require chronological disjoint complete training and held-out years")
    years = sorted(set(train + hold))
    rain_station = str(cfg.get("precipitation_station", "72640014839"))
    met_station = str(cfg.get("met_station", "72640594869"))
    lat, lon = float(cfg.get("latitude", 43.0547222)), float(cfg.get("longitude", -88.0461111))
    download = bool(cfg.get("download", True))
    metadata = fetch_cached(raw / "isd-history.csv", HISTORY_URL, download)
    history = pd.read_csv(metadata, dtype=str, keep_default_na=False)
    station_rows = []
    for station in (rain_station, met_station):
        match = history[(history["USAF"] + history["WBAN"]) == station]
        if len(match) != 1:
            raise ValueError(f"Cannot unambiguously verify NOAA product station {station}")
        station_rows.append(match.iloc[0].to_dict())
    tasks=[]
    for year in years:
        for station in sorted(set((rain_station, met_station))):
            # Keep each request within its calendar year. The service can return
            # only the trailing day for a query spanning more than one year.
            params = dict(dataset="global-hourly", stations=station,
                          startDate=f"{year}-01-01", endDate=f"{year}-12-31",
                          format="csv", dataTypes="TMP,WND,AA1,AA2,AA3,AA4", includeAttributes="true",
                          includeStationName="true", includeStationLocation="true")
            tasks.append((raw / f"NOAA_{station}_{year}.csv", NOAA_URL + "?" + urlencode(params)))
        params = dict(parameters="ALLSKY_SFC_SW_DWN", community="RE", longitude=lon, latitude=lat,
                      start=f"{year}0101", end=f"{year}1231", format="JSON", **{"time-standard":"UTC"})
        tasks.append((raw / f"POWER_solar_{year}.json", POWER_URL + "?" + urlencode(params)))
    # Bounded parallel read-only requests, with all failed requests reported.
    failures=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = [(path, pool.submit(fetch_cached, path, url, download)) for path,url in tasks]
        for path,future in futures:
            try:
                future.result()
            except Exception as exc:
                failures.append({"file":str(path.relative_to(root)), "error":str(exc)})
    if failures:
        (output.parent/"download_failures.json").write_text(json.dumps(failures, indent=2)+"\n")
        raise RuntimeError(f"{len(failures)} real datasets unavailable; see data/processed/download_failures.json. No synthetic substitute.")
    reports=[]; frames=[]
    for year in years:
        start = pd.Timestamp(f"{year}-01-01", tz="UTC")
        stop = pd.Timestamp(f"{year+1}-01-01", tz="UTC")
        index = pd.date_range(start,stop, freq="h", inclusive="left")
        target = index + pd.Timedelta(minutes=30)
        rain = pd.read_csv(raw/f"NOAA_{rain_station}_{year}.csv", dtype=str,keep_default_na=False)
        met = pd.read_csv(raw/f"NOAA_{met_station}_{year}.csv", dtype=str,keep_default_na=False)
        for table, station in ((rain,rain_station),(met,met_station)):
            if "STATION" not in table or set(table.STATION) != {station}:
                raise ValueError(f"Unexpected station in NOAA data: {station}")
            available_months = set(pd.to_datetime(table.DATE).dt.month)
            if len(available_months) != 12:
                raise ValueError(f"NOAA {station} {year}: response covers only {sorted(available_months)} months; incomplete retrieval, not a usable annual response")
        intervals, rain_counts = precipitation_intervals(rain)
        p, coverage, trace = align_precipitation(intervals,index)
        obs = met[met.REPORT_TYPE.str.strip().isin(["FM-15","FM-16"])]
        times = pd.to_datetime(obs.DATE, utc=True)
        tv = [parse_value(v,0,1,True) for v in obs.TMP]
        wv = [parse_value(v,3,4) for v in obs.WND]
        temperature, t_interp, t_fill = interpolate_met(times,tv,target,cfg.get("max_met_gap_hours",2))
        wind, w_interp, w_fill = interpolate_met(times,wv,target,cfg.get("max_met_gap_hours",2))
        solar_data = json.loads((raw/f"POWER_solar_{year}.json").read_text())
        header = solar_data.get("header",{})
        if header.get("time_standard") != "UTC":
            raise ValueError("POWER did not return explicitly UTC data")
        units = solar_data["parameters"]["ALLSKY_SFC_SW_DWN"]["units"]
        if units not in {"Wh/m^2", "Wh/m2", "W-hr/m^2", "W/m^2"}:
            raise ValueError(f"Unverified POWER solar units: {units}")
        vals = solar_data["properties"]["parameter"]["ALLSKY_SFC_SW_DWN"]
        solar = pd.Series({pd.to_datetime(k, format="%Y%m%d%H",utc=True):float(v) for k,v in vals.items()}).reindex(index)
        solar[(solar < 0) | (solar > 1500)] = np.nan
        f = pd.DataFrame({"timestamp_utc":index,"precipitation_mm_per_h":p,"temperature_C":temperature,
                          "solar_W_per_m2":solar.to_numpy(),"wind_m_per_s":wind,
                          "precipitation_coverage_fraction":coverage,"precipitation_trace_lower_bound":trace,
                          "temperature_interpolated":t_interp,"wind_interpolated":w_interp,
                          "temperature_short_gap_filled":t_fill,"wind_short_gap_filled":w_fill,
                          "split":"train" if year in train else "validation"})
        reports.append({"year":year,"expected_hours":len(index),"valid_counts":{v:int(f[v].notna().sum()) for v in VARIABLES},
                        "missing_fraction":{v:float(f[v].isna().mean()) for v in VARIABLES},
                        "remaining_missing_runs":{v:gap_statistics(f[v]) for v in VARIABLES},
                        "complete_joint_hours":int(f[VARIABLES].notna().all(axis=1).sum()),"rainfall_qc":rain_counts,
                        "noaa_raw_rows":{"precipitation_station":len(rain),"met_station":len(met)},
                        "accepted_met_report_count":len(obs),
                        "temperature_rejected_or_missing_reports":int(np.isnan(tv).sum()),
                        "wind_rejected_or_missing_reports":int(np.isnan(wv).sum()),
                        "trace_affected_hours":int(trace.sum()),"temperature_interpolation_fraction":float(t_interp.mean()),
                        "wind_interpolation_fraction":float(w_interp.mean()),"temperature_short_gap_imputed":int(t_fill.sum()),
                        "wind_short_gap_imputed":int(w_fill.sum()),"solar_source":header.get("sources"),"solar_api":header.get("api"),
                        "solar_native_units":units,"quantified_precipitation_on_complete_hours_mm":float(np.nansum(p)),
                        "solar_energy_Wh_per_m2":float(solar.sum())})
        frames.append(f)
    result=pd.concat(frames,ignore_index=True)
    result.to_csv(output,index=False,float_format="%.10g")
    # Read saved precision back: subsequent steps consume exactly this artifact.
    result=load_prepared(cfg,root)
    report={"status":"real_data_prepared","train_years":train,"validation_years":hold,
            "timezone":"UTC","analysis_timezone":"America/Chicago","hour_label":"interval_start",
            "station_metadata":station_rows,"power_location":{"latitude":lat,"longitude":lon},"years":reports,
            "processed_file":str(output.relative_to(root)),"processed_sha256":sha256(output),
            "source_manifest":[json.loads(metadata.with_suffix(metadata.suffix+".request.json").read_text())] +
                              [json.loads(path.with_suffix(path.suffix+".request.json").read_text()) for path,_ in tasks],
            "rules":["Only independently reported 1 h FM-15 accumulations; exclude multihour/special reports and overlapping intervals.",
                     "Actual NOAA accumulation end times retained; UTC-hour overlap weighting assumes constant rate within each observed hour.",
                     "A target hour without complete observation coverage remains NaN; never missing-to-zero precipitation.",
                     "Trace is flagged zero quantified lower bound; it is not observed dry weather or a known exact amount.",
                     f"Temperature and wind are linear midpoint estimates from accepted observations, maximum {cfg.get('max_met_gap_hours',2)} h bracket; no extrapolation.",
                     "Bracket >75 min counted as short-gap filling separately from ordinary grid interpolation.",
                     "POWER hourly energy divided by 1 h; hour-start anchor is an inherited explicit processing assumption, not independently established time bounds.",
                     "Station observations and gridded solar estimates have different spatial support; spatial transfers are not local calibration."]}
    (output.parent/"data_quality_report.json").write_text(json.dumps(report,indent=2)+"\n")
    return result
