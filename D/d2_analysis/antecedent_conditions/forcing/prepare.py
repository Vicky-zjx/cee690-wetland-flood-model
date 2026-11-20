"""Prepare four exact cached-weather prefixes; no downloads or model execution."""
from __future__ import annotations
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
from C.uncertainty.data import VARIABLES, precipitation_intervals, align_precipitation, parse_value, interpolate_met
from D.g5_diagnostics import runner as historical
from D.d2_analysis.antecedent_conditions.forcing import provider
SOURCE = ROOT/'C/data/processed/weather_hourly.csv'
HIST = ROOT/'D/d2_analysis/historical_storm_robustness'
EVENTS = ('historical_20190825T0200Z', 'historical_20240804T1500Z')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean(value):
    if isinstance(value, dict): return {k:clean(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)): return [clean(v) for v in value]
    if isinstance(value, np.generic): return clean(value.item())
    if isinstance(value, (pd.Timestamp, Path)): return str(value)
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def write(path, obj):
    path.write_text(json.dumps(clean(obj), indent=2, allow_nan=False)+'\n')


@lru_cache(maxsize=2)
def raw_year(year):
    files = {}
    for name in (f'NOAA_72640014839_{year}.csv', f'NOAA_72640594869_{year}.csv', f'POWER_solar_{year}.json'):
        path = ROOT/'C/data/raw'/name
        request = json.loads(Path(str(path)+'.request.json').read_text())
        if not request.get('ok') or sha(path) != request['sha256']:
            raise ValueError('Original raw cached source checksum failed: '+name)
        files[name] = dict(file=str(path.relative_to(ROOT)), sha256=sha(path), request=request)
    rain = pd.read_csv(ROOT/f'C/data/raw/NOAA_72640014839_{year}.csv', dtype=str, keep_default_na=False)
    intervals, qc = precipitation_intervals(rain)
    met = pd.read_csv(ROOT/f'C/data/raw/NOAA_72640594869_{year}.csv', dtype=str, keep_default_na=False)
    obs = met[met.REPORT_TYPE.str.strip().isin(['FM-15', 'FM-16'])]
    times = pd.to_datetime(obs.DATE, utc=True)
    series = {}
    for name, values in [('temperature_C', [parse_value(v,0,1,True) for v in obs.TMP]),
                         ('wind_m_per_s', [parse_value(v,3,4) for v in obs.WND])]:
        values = pd.Series(values, index=pd.DatetimeIndex(times)).dropna().sort_index()
        series[name] = values.groupby(level=0).agg(lambda a: a.iloc[0] if a.nunique()==1 else np.nan).dropna()
    solar = json.loads((ROOT/f'C/data/raw/POWER_solar_{year}.json').read_text())
    assert solar['header']['time_standard'] == 'UTC'
    return files, intervals, qc, series, solar


def native_audit(frame, cid):
    start, end = frame.timestamp_utc.iloc[0], frame.timestamp_utc.iloc[-1]+pd.Timedelta(hours=1)
    files, intervals, qc, series, solar = raw_year(int(start.year))
    selected = [(a,b,p,tr) for a,b,p,tr in intervals if a<end and b>start]
    records = []
    for a,b,amount,trace in selected:
        overlap = (min(b,end)-max(a,start)).total_seconds()/3600
        records.append(dict(interval_start_utc=a, interval_end_utc=b, duration_h=(b-a).total_seconds()/3600,
            observed_depth_mm=amount, trace_quantified_lower_bound=trace,
            overlap_with_antecedent_h=overlap, attributed_antecedent_depth_mm=amount*overlap,
            rule='Uniform rate within actual native1h interval; overlap cropping only'))
    pd.DataFrame(records).to_csv(HERE/f'{cid}_native_rain.csv', index=False)
    p, coverage, trace = align_precipitation(intervals, pd.DatetimeIndex(frame.timestamp_utc))
    rain_error = float(np.max(abs(p-frame.precipitation_mm_per_h.to_numpy())))
    assert rain_error<=5e-8 and np.all(np.isclose(coverage,1,atol=1e-6,rtol=0))
    assert np.array_equal(trace, frame.precipitation_trace_lower_bound.to_numpy())
    targets = pd.DatetimeIndex(frame.timestamp_utc)+pd.Timedelta(minutes=30)
    metrows, errors = [], {}
    for name, values in series.items():
        estimate, interpolated, short_gap = interpolate_met(values.index,values.to_numpy(),targets,2)
        errors[name] = float(np.max(abs(estimate-frame[name].to_numpy())))
        assert errors[name]<=5e-8
        for target, result in zip(targets, estimate):
            i=int(values.index.searchsorted(target)); exact=i<len(values) and values.index[i]==target
            lo=hi=i if exact else i-1
            if not exact: hi=i
            assert lo>=0 and hi<len(values)
            left,right=values.index[lo],values.index[hi]
            metrows.append(dict(variable=name,target_midpoint_utc=target,left_observation_utc=left,right_observation_utc=right,
                left_observed_value=float(values.iloc[lo]),right_observed_value=float(values.iloc[hi]),
                bracketing_span_minutes=(right-left).total_seconds()/60,interpolated_value=float(result),
                exact_observation=bool(exact),method='Existing C accepted bracketing interpolation, maximum2h, no extrapolation'))
    pd.DataFrame(metrows).to_csv(HERE/f'{cid}_native_met_brackets.csv',index=False)
    vals=solar['properties']['parameter']['ALLSKY_SFC_SW_DWN']
    solar_values=np.array([vals[t.strftime('%Y%m%d%H')] for t in frame.timestamp_utc],dtype=float)
    solar_error=float(np.max(abs(solar_values-frame.solar_W_per_m2.to_numpy())))
    assert solar_error<=5e-8
    return dict(raw_sources=files,native_rain_selected_intervals=len(selected),
        native_rain_end_minutes=sorted({b.minute for _,b,_,_ in selected}),
        native_trace_intervals=sum(bool(tr) for _,_,_,tr in selected),
        native_rain_reconstruction_error_mm=rain_error, native_met_reconstruction_error=errors,
        solar_cached_grid_reconstruction_error_W_m2=solar_error,
        solar_native_units=solar['parameters']['ALLSKY_SFC_SW_DWN']['units'],
        solar_source=solar['header'].get('sources'), rain_parser_annual_qc=qc,
        actual_interval_files=[f'{cid}_native_rain.csv',f'{cid}_native_met_brackets.csv'])


def main():
    if (HERE/'manifest.json').exists():
        raise FileExistsError('Prepared antecedent manifest exists; do not overwrite frozen inputs')
    source_manifest=HIST/'events/manifest.json'
    old=json.loads(source_manifest.read_text())
    if sha(SOURCE)!=old['standardized_source']['sha256']:
        raise ValueError('Stage C hourly source differs from frozen historical-event provenance')
    frame=pd.read_csv(SOURCE,float_precision='round_trip');frame.timestamp_utc=pd.to_datetime(frame.timestamp_utc,utc=True)
    lines=SOURCE.read_bytes().splitlines(keepends=True)
    assert len(lines)==len(frame)+1
    rows=[];validation=[]
    for event in old['events']:
        if event['event_id'] not in EVENTS:continue
        original=Path(event['forcing_file']);event_bytes=original.read_bytes()
        assert sha(original)==event['forcing_sha256']
        start=pd.Timestamp(event['start_utc']); first=int(np.flatnonzero((frame.timestamp_utc==start).to_numpy())[0])
        for hours in (24,72):
            window=frame.iloc[first-hours:first].copy().reset_index(drop=True)
            assert len(window)==hours and window.timestamp_utc.iloc[0]==start-pd.Timedelta(hours=hours)
            assert np.isfinite(window[VARIABLES].to_numpy()).all() and (window.temperature_C>0).all()
            assert np.all(np.isclose(window.precipitation_coverage_fraction,1,rtol=0,atol=1e-6))
            cid=f"{event['event_id']}__antecedent_{hours}h"
            prefix=lines[0]+b''.join(lines[first-hours+1:first+1])
            path=HERE/f'{cid}.csv'
            original_lines=event_bytes.splitlines(keepends=True)
            assert lines[0]==original_lines[0]
            path.write_bytes(prefix+b''.join(original_lines[1:]))
            saved=pd.read_csv(path,float_precision='round_trip');saved.timestamp_utc=pd.to_datetime(saved.timestamp_utc,utc=True)
            pd.testing.assert_frame_equal(saved.iloc[:hours].reset_index(drop=True),window,check_exact=True,check_dtype=False)
            audit=native_audit(window,cid)
            record=dict(conditioning_id=cid, scenario_id=cid,event_id=event['event_id'],antecedent_hours=hours,
                event_offset_s=hours*3600,total_duration_s=(hours+144)*3600,event_duration_s=144*3600,
                forcing_file=str(path),forcing_sha256=sha(path),start_utc=window.timestamp_utc.iloc[0],
                event_start_utc=start,end_utc=event['end_utc'],seed_sequence_entropy=None,
                original_event_forcing_file=str(original),original_event_forcing_sha256=event['forcing_sha256'],
                original_event_rows_byte_exact=True,antecedent_rows_source_byte_exact=True,
                antecedent_rain_mm=float(window.precipitation_mm_per_h.sum()),
                antecedent_solar_energy_Wh_m2=float(window.solar_W_per_m2.sum()),
                antecedent_min_temperature_C=float(window.temperature_C.min()),
                antecedent_missing_variables=window[VARIABLES].isna().sum().to_dict(),
                antecedent_qc={k:int(window[k].sum()) for k in ['precipitation_trace_lower_bound','temperature_interpolated','wind_interpolated','temperature_short_gap_filled','wind_short_gap_filled']},
                native_audit=audit,
                mapping='Prefix and original event separately mapped by unchanged hourly adapter; event shifted only; no T/wind interpolation across seam',
                source_standardized_csv=str(SOURCE.relative_to(ROOT)),source_standardized_sha256=sha(SOURCE))
            write(path.with_suffix('.meta.json'),record)
            for dt in (15.,7.5):
                fixed=historical.load_fixed(dt=dt,design='best')
                fixed['model_config']['numerics']['duration_s']=record['total_duration_s']
                fixed['manifest']['duration_s']=record['total_duration_s']
                combined,info=provider.read_forcing(path)
                check=provider.forcing_conservation(combined,fixed)
                assert check['pass_conservation'] and check['event_weather_and_boundary_inputs_bit_exact']
                validation.append(dict(conditioning_id=cid,dt_s=dt,**check))
            record['mapping_validation']=validation[-2:]
            write(path.with_suffix('.meta.json'),record)
            rows.append(record)
    manifest=dict(cases=rows,case_count=len(rows),model_runs_executed=0,downloads=0,
        source_event_manifest=str(source_manifest.relative_to(ROOT)),source_event_manifest_sha256=sha(source_manifest),
        source_code_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'provider.py',ROOT/'C/uncertainty/data.py',ROOT/'D/g5_diagnostics/recovered/model_adapters.py']},
        event_mapping_preserved=True,original_event_duration_h=144,
        interpretation='Continuous dynamical conditioning is implemented by caller; these are observed regional weather prefixes, not observed initial storages or a spin-up convergence claim',
        limitations=['Original C rainfall clock-hour redistribution assumes constant rate within actual one-hour accumulations ending at minute52.',
        'Trace rain remains a zero quantified lower bound; missing rain is never zero-filled.',
        'Temperature/wind use accepted regional station interpolation at hourly midpoints; no operational availability claim.',
        'Solar is a POWER grid estimate with inherited UTC start-label assumption; hourly energy is preserved.',
        'Separate prefix/event endpoint holds cause an explicit seam jump in T/wind; this preserves every original event input value.',
        'Fixed humidity/pressure and other boundaries remain unchanged rather than reconstructed for each antecedent date.',
        'Starting storages and controller state are assumptions; conditioning duration alone does not establish observed antecedent conditions.'])
    write(HERE/'manifest.json',manifest)
    write(HERE/'mapping_validation.json',dict(passed=True,checks=validation,checks_count=len(validation),model_runs_executed=0))
    pd.DataFrame([{k:r[k] for k in ['conditioning_id','event_id','antecedent_hours','start_utc','event_start_utc','end_utc','antecedent_rain_mm','antecedent_solar_energy_Wh_m2','antecedent_min_temperature_C','forcing_sha256']} for r in rows]).to_csv(HERE/'forcing_summary.csv',index=False)
    note=['# Antecedent weather preparation','',
        'Four joint-weather inputs prepend the immediately preceding 24 or72 hours to the unchanged August2019 and August2024 events. Only existing cached StageC observations and grid estimates are read. No weather is generated, scaled, downloaded or filled. Original event CSV rows and all QC flags remain byte-for-byte unchanged. Each additional row is copied exactly from the original standardized CSV.','',
        '## Inputs and time semantics','',
        '| Event | Prefix (h) | Prefix rain (mm) | Solar energy (Wh/m²) | Minimum T (°C) |','|---|---:|---:|---:|---:|']
    for r in rows:note.append(f"| {r['event_id']} | {r['antecedent_hours']} | {r['antecedent_rain_mm']:.3f} | {r['antecedent_solar_energy_Wh_m2']:.2f} | {r['antecedent_min_temperature_C']:.2f} |")
    note += ['',
        'All requested prefixes contain complete finite joint weather, complete native rainfall coverage and positive temperature. The2019 prefixes and2024 24h prefix contain zero quantified rain, with no trace flags in those prefixes; the2024 72h prefix contains18.5mm and six trace-affected UTC hours. These statements describe the cached station record transferred to the hypothetical installation, not spatially exact local precipitation.','',
        'Native one-hour precipitation intervals, their actual starts/ends, boundary-overlap allocation and trace flags are saved in each `_native_rain.csv`. Native reports end at minute52, so the UTC-hour rates are derived overlap allocations rather than native clock-hour observations. Each `_native_met_brackets.csv` records actual accepted temperature/wind observation times bracketing the midpoint target and the interpolation span. The raw file and request hashes, station sources, POWER units and all reconstruction errors are retained in each metadata file. No new interpolation changes the standardized rows.','',
        '## Preserve the original event mapping','',
        'The provider maps the prefix and144h event separately with the unchanged original hourly-input function, then shifts the event time grid by the prefix length. Precipitation and solar are constant within each UTC interval; temperature and wind are linear between their hourly midpoints with a half-hour endpoint hold. A naive single concatenated mapping would alter the original event\'s first half-hour temperature/wind values and is therefore not used.','',
        'The deliberate mapping seam retains the independent endpoint holds. At the event start the2019 temperature/wind changes are−1.5°C/−1.65m/s; the2024 changes are+1.27631579°C/0m/s. These are explicit discontinuities of the retained representation, not measured instantaneous jumps. They occur at an exact hourly and120s controller boundary. The caller must continue the physical state, actual actuator opening and controller latch state across this seam without a reset.','',
        'All eight preparation checks (four inputs at15s and7.5s) conserve rainfall depth and solar energy and reproduce every original event weather/fixed-boundary input value bit-for-bit after the time shift. These are input-mapping checks; no ODE, safety, model convergence or antecedent-state validation is performed here. The first two rows of the original input mapper are historical gate placeholders; the continuous feedback runner replaces them with its unchanged policy, so event equality is checked on rows2 onward.','',
        '## Provider contract','',
        '`record_for(event_id, hours)` returns one JSON-serializable manifest record. `read_forcing(path)` returns the checked combined frame and diagnostics, with seam metadata attached to `frame.attrs`. `hourly_input_segments(frame,fixed)` and `forcing_conservation(frame,fixed)` match the original adapter signatures. The caller privately extends the numerical duration to168h or216h and delegates all other model functions unchanged.','',
        'Run `/opt/miniconda3/bin/python3 -B D/d2_analysis/antecedent_conditions/forcing/prepare.py` from the project root to prepare a new, empty forcing destination. Existing completed manifests are not overwritten.','',
        'The24/72h conditioning periods are controlled experimental choices, not evidence of equilibration or observed initial soil/groundwater/wetland storage. All existing fixed boundaries, including humidity/pressure assumptions, remain unchanged. These prefixes do not alter the original144h event metrics window.','']
    (HERE/'README.md').write_text('\n'.join(note))
    print(json.dumps(clean(dict(cases=len(rows),mapping_checks=len(validation),event_input_exact=True,model_runs=0,prefix_rain_mm=[r['antecedent_rain_mm'] for r in rows]))))

if __name__=='__main__':main()
