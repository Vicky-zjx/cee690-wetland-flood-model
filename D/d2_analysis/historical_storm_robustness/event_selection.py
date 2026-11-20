"""Select fixed historical rain windows from frozen Stage C observations only.

No model output, parameter outcome or forecast enters selection. No network
access, missing-value replacement, rainfall scaling or weather generation occurs.
The original September 2022 benchmark remains its distinct approved forcing.
"""
from __future__ import annotations
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
from C.uncertainty.data import VARIABLES,precipitation_intervals,align_precipitation
from B.project_paths import resolve_path
from D.d2_analysis.historical_storm_robustness import preservation
SOURCE=ROOT/'C/data/processed/weather_hourly.csv'
BENCHMARK_START=pd.Timestamp('2022-09-09T00:00:00Z')
BENCHMARK_END=BENCHMARK_START+pd.Timedelta(hours=144)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def clean(v):
    if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple,np.ndarray)):return [clean(x) for x in v]
    if isinstance(v,np.generic):return clean(v.item())
    if isinstance(v,(pd.Timestamp,Path)):return str(v)
    if isinstance(v,float) and not np.isfinite(v):return None
    return v

def write(path,value):Path(path).write_text(json.dumps(clean(value),indent=2,allow_nan=False)+'\n')


def candidates(frame):
    """Every144h window; rank its fixed hours24..47 rainfall core."""
    n=len(frame)-143;data={}
    rain=frame[VARIABLES[0]];rolling24=rain.rolling(24,min_periods=24).sum()
    for variable in VARIABLES:
        data['missing_'+variable]=frame[variable].isna().astype(int).rolling(144).sum().shift(-143).iloc[:n].astype(int).to_numpy()
    data.update(window_start_index=np.arange(n),start_utc=frame.timestamp_utc.iloc[:n].to_numpy(),
        end_utc=(frame.timestamp_utc.iloc[:n]+pd.Timedelta(hours=144)).to_numpy(),
        ranking_core_start_utc=(frame.timestamp_utc.iloc[:n]+pd.Timedelta(hours=24)).to_numpy(),
        ranking_core_end_utc=(frame.timestamp_utc.iloc[:n]+pd.Timedelta(hours=48)).to_numpy(),
        ranking_24h_mm=rolling24.shift(-47).iloc[:n].to_numpy(),
        max_24h_mm=rolling24.rolling(121,min_periods=1).max().shift(-143).iloc[:n].to_numpy(),
        max_hourly_mm=rain.rolling(144,min_periods=1).max().shift(-143).iloc[:n].to_numpy(),
        rainfall_total_mm=rain.rolling(144,min_periods=144).sum().shift(-143).iloc[:n].to_numpy(),
        partial_observed_rainfall_mm=rain.rolling(144,min_periods=1).sum().shift(-143).iloc[:n].to_numpy(),
        minimum_temperature_C=frame.temperature_C.rolling(144,min_periods=1).min().shift(-143).iloc[:n].to_numpy(),
        nonpositive_temperature_hours=(frame.temperature_C<=0).astype(int).rolling(144).sum().shift(-143).iloc[:n].astype(int).to_numpy(),
        incomplete_rain_coverage_hours=pd.Series((~np.isclose(frame.precipitation_coverage_fraction,1.,atol=1e-6,rtol=0)).astype(int)).rolling(144).sum().shift(-143).iloc[:n].astype(int).to_numpy())
    table=pd.DataFrame(data)
    table['complete_joint_window']=table[['missing_'+v for v in VARIABLES]].sum(axis=1)==0
    table['eligible']=table.complete_joint_window&(table.nonpositive_temperature_hours==0)&(table.incomplete_rain_coverage_hours==0)&table.ranking_24h_mm.notna()
    table['overlaps_reserved_benchmark']=(table.start_utc<BENCHMARK_END)&(table.end_utc>BENCHMARK_START)
    table['qc_exclusion_reasons']=['; '.join([*(v+' missing' for v in VARIABLES if r['missing_'+v]>0),
        *(['nonpositive temperature'] if r.nonpositive_temperature_hours else []),
        *(['incomplete native-rain coverage'] if r.incomplete_rain_coverage_hours else []),
        *(['incomplete24h ranking core'] if not np.isfinite(r.ranking_24h_mm) else [])]) for _,r in table.iterrows()]
    ranked=table.sort_values(['ranking_24h_mm','ranking_core_start_utc'],ascending=[False,True],kind='stable',na_position='last').copy()
    ranked['rank']=np.arange(1,len(ranked)+1)
    chosen=[]
    for _,row in ranked[ranked.eligible&~ranked.overlaps_reserved_benchmark].iterrows():
        if any(row.start_utc<old.end_utc and row.end_utc>old.start_utc for old in chosen):continue
        chosen.append(row)
        if len(chosen)==4:break
    if len(chosen)!=4:raise RuntimeError('Cached data do not support four distinct complete nonfreezing additional windows')
    ids={int(r.window_start_index):'historical_'+r.start_utc.strftime('%Y%m%dT%H%MZ') for r in chosen}
    ranked['selected_event_id']=ranked.window_start_index.map(ids)
    return ranked,chosen


def wet_spells(frame):
    positions=np.flatnonzero(frame.precipitation_mm_per_h.to_numpy()>0)
    if not len(positions):return dict(positive_rain_hours=0,first_to_last_positive_span_h=0,largest_quantified_wet_spell=None)
    groups=np.split(positions,np.flatnonzero(np.diff(positions)>=7)+1)
    spells=[]
    for group in groups:
        a,b=int(group[0]),int(group[-1]);subset=frame.iloc[a:b+1]
        spells.append(dict(start_utc=frame.timestamp_utc.iloc[a],end_utc=frame.timestamp_utc.iloc[b]+pd.Timedelta(hours=1),
            span_h=b-a+1,positive_hours=len(group),depth_mm=float(subset.precipitation_mm_per_h.sum())))
    largest=sorted(spells,key=lambda s:(-s['depth_mm'],s['start_utc']))[0]
    return dict(positive_rain_hours=len(positions),first_to_last_positive_span_h=int(positions[-1]-positions[0]+1),
        quantified_wet_spell_count=len(spells),largest_quantified_wet_spell=largest,
        spell_definition='Positive quantified hours joined unless separated by at least6 quantified-zero hours. Trace-zero is not verified dry; exploratory operational definition, not certified storm delineation.')


@lru_cache(maxsize=6)
def native_year(year):
    path=ROOT/f'C/data/raw/NOAA_72640014839_{year}.csv'
    raw=pd.read_csv(path,dtype=str,keep_default_na=False)
    intervals,qc=precipitation_intervals(raw)
    return path,intervals,qc


def native_audit(window,event_id,out):
    year=int(window.timestamp_utc.iloc[0].year);path,intervals,qc=native_year(year)
    start=window.timestamp_utc.iloc[0];stop=start+pd.Timedelta(hours=144)
    selected=[r for r in intervals if r[0]<stop and r[1]>start]
    records=[]
    for a,b,depth,trace in selected:
        overlap=(min(b,stop)-max(a,start)).total_seconds()/3600
        records.append(dict(interval_start_utc=a,interval_end_utc=b,duration_h=(b-a).total_seconds()/3600,
            observed_depth_mm=depth,trace_quantified_lower_bound=bool(trace),overlap_with_model_window_h=overlap,
            attributed_window_depth_mm=depth*overlap,attribution_rule='Constant rate within the observed native one-hour interval; only boundary cropping'))
    table=pd.DataFrame(records);table.to_csv(out/(event_id+'_native_rain_intervals.csv'),index=False)
    reconstructed,coverage,trace=align_precipitation(intervals,pd.DatetimeIndex(window.timestamp_utc))
    error=float(np.max(np.abs(reconstructed-window.precipitation_mm_per_h.to_numpy())))
    if not np.isfinite(error) or error>5e-8 or not np.all(np.isclose(coverage,1,atol=1e-6,rtol=0)):
        raise RuntimeError('Stored C rain does not match audited native observations within saved ten-significant-digit precision')
    complete=[r for r in selected if r[0]>=start and r[1]<=stop]
    ends=pd.DatetimeIndex([r[1] for r in complete]);depth=np.array([r[2] for r in complete])
    rolls=pd.Series(depth).rolling(24).sum().to_numpy();span=np.full(len(complete),np.nan)
    if len(complete)>=24:span[23:]=(ends[23:].asi8-pd.DatetimeIndex([r[0] for r in complete[:-23]]).asi8)/3.6e12
    native24=rolls[np.isclose(span,24)]
    return dict(raw_file=str(path.relative_to(ROOT)),raw_sha256=sha(path),native_selected_intervals=len(selected),
        complete_native_one_hour_intervals=len(complete),native_end_minutes=sorted(set(r[1].minute for r in selected)),
        maximum_native_one_hour_depth_mm=max(r[2] for r in complete),
        maximum_24_consecutive_native_hours_mm=float(np.max(native24)) if len(native24) else None,
        native_trace_intervals_intersecting_window=sum(bool(r[3]) for r in selected),
        overlap_reconstruction_max_error_mm=error,represented_native_window_depth_mm=sum(r['attributed_window_depth_mm'] for r in records),
        interpretation='Native interval depth is an observation; UTC-hour amount is a conservative redistribution with assumed within-hour timing.',
        annual_parser_qc=qc)


def benchmark(frame,out):
    config=json.loads((ROOT/'B/no_spillway_config.json').read_text())
    path=resolve_path(config['forcing']['prepared_csv']);meta_path=resolve_path(config['forcing']['manifest_json'])
    original=pd.read_csv(path,float_precision='round_trip');t=original.time_s.to_numpy();P=original.P.to_numpy()*1000
    native_path=ROOT/'B/work/event_data/weather/KMKE_rain_selected_intervals_utc.csv'
    native=pd.read_csv(native_path,float_precision='round_trip')
    native_start=pd.to_datetime(native.interval_start_utc,utc=True);native_end=pd.to_datetime(native.interval_end_utc,utc=True)
    native_event=native[(native_start<BENCHMARK_END)&(native_end>BENCHMARK_START)]
    native_peak=float(native_event.loc[native_event.duration_h==1,'depth_mm_quantified'].max())
    expected=json.loads(meta_path.read_text())['sha256']
    if sha(path)!=expected:raise RuntimeError('Original benchmark forcing checksum failed')
    hourly=[]
    for k in range(144):
        mask=(t>=3600*k)&(t<=3600*(k+1));local=t[mask]
        hourly.append(dict(timestamp_utc=BENCHMARK_START+pd.Timedelta(hours=k),
            original_precipitation_mm=float(np.trapezoid(P[mask],local)),
            original_solar_Wh_m2=float(np.trapezoid(original.R_s.to_numpy()[mask],local)/3600),
            original_temperature_at_midpoint_C=float(np.interp(3600*k+1800,t,original.T_a)),
            original_wind_at_midpoint_m_s=float(np.interp(3600*k+1800,t,original.v))))
    hourly=pd.DataFrame(hourly);c=frame[(frame.timestamp_utc>=BENCHMARK_START)&(frame.timestamp_utc<BENCHMARK_END)].reset_index(drop=True)
    for variable in VARIABLES:hourly['C_'+variable]=c[variable]
    hourly.to_csv(out/'benchmark_C_comparison.csv',index=False)
    pairs=[('rainfall_hour_mm','original_precipitation_mm','C_precipitation_mm_per_h'),
        ('solar_hour_Wh_m2','original_solar_Wh_m2','C_solar_W_per_m2'),
        ('temperature_midpoint_C','original_temperature_at_midpoint_C','C_temperature_C'),
        ('wind_midpoint_m_s','original_wind_at_midpoint_m_s','C_wind_m_per_s')]
    differences={}
    for name,a,b in pairs:
        valid=hourly[b].notna();difference=hourly.loc[valid,a]-hourly.loc[valid,b]
        differences[name]=dict(common_hours=int(valid.sum()),C_missing_hours=int((~valid).sum()),
            max_absolute_difference=float(np.abs(difference).max()),mean_absolute_difference=float(np.abs(difference).mean()))
    comparison=dict(window_start_utc=BENCHMARK_START,window_end_utc=BENCHMARK_END,
        original_forcing_file=str(path.relative_to(ROOT)),original_forcing_sha256=sha(path),
        original_forcing_manifest_sha256=sha(meta_path),original_quantified_rain_mm=float(np.trapezoid(P,t)),
        C_quantified_sum_over_observed_hours_mm=float(c.precipitation_mm_per_h.sum()),
        C_rainfall_total_complete=False,C_missing_rain_hours=c.loc[c.precipitation_mm_per_h.isna(),'timestamp_utc'].tolist(),
        original_solar_energy_Wh_m2=float(np.trapezoid(original.R_s,t)/3600),C_solar_energy_Wh_m2=float(c.solar_W_per_m2.sum()),
        common_hour_comparison=differences,
        known_processing_difference='Original approved benchmark includes an actual3h trace accumulation spanning the two missing native one-hour reports, uses120s mass-preserving rain/solar ramps and direct met interpolation. C uses strict independent1h rain intervals, overlap-weighted clock hours, midpoint T/wind estimates, and preserves three incomplete clock-hour rain bins as missing. No missing C rain is filled or used to replace the benchmark.')
    write(out/'benchmark_C_comparison.json',comparison)
    p=hourly.original_precipitation_mm
    record=dict(event_id='benchmark_september_2022',scenario_id='benchmark_september_2022',benchmark=True,
        start_utc=BENCHMARK_START,end_utc=BENCHMARK_END,forcing_file=str(path),forcing_sha256=sha(path),seed_sequence_entropy=None,
        forcing_schema='Original prepared120s model inputs; benchmark flag selects original deterministic execution, not the hourly event adapter',
        rainfall_total_mm=float(np.trapezoid(P,t)),max_hourly_mm=float(p.max()),max_24h_mm=float(p.rolling(24).sum().max()),
        solar_energy_Wh_m2=float(np.trapezoid(original.R_s,t)/3600),minimum_temperature_C=float(original.T_a.min()),
        complete_forcing=True,hours=144,source='Original approved September2022 KMKE/KMWC/POWER reconstruction, preserved byte-for-byte',
        native_maximum_one_hour_mm=native_peak,
        provenance=dict(prepared_manifest=str(meta_path.relative_to(ROOT)),prepared_manifest_sha256=sha(meta_path),comparison_report='benchmark_C_comparison.json',
            native_peak_source=str(native_path.relative_to(ROOT)),native_intervals_sha256=sha(native_path)))
    temp=pd.DataFrame(dict(timestamp_utc=hourly.timestamp_utc,precipitation_mm_per_h=p));record.update(wet_spells(temp))
    return record,comparison


def main():
    guard=preservation.verify()
    if not guard['passed']:raise RuntimeError('Protected files changed before event preparation')
    out=HERE/'events'
    if (out/'manifest.json').exists():raise FileExistsError('Completed event manifest already exists; never overwrite it')
    frame=pd.read_csv(SOURCE,float_precision='round_trip');frame.timestamp_utc=pd.to_datetime(frame.timestamp_utc,utc=True)
    quality=json.loads((ROOT/'C/data/processed/data_quality_report.json').read_text())
    if sha(SOURCE)!=quality['processed_sha256']:raise RuntimeError('Standardized source hash differs from retained QC report')
    if len(frame)!=52608 or not np.all(np.diff(frame.timestamp_utc.astype('int64'))==3600*10**9):raise RuntimeError('Unexpected source coverage/calendar gaps')
    sources=[]
    for year in range(2019,2025):
        for name in (f'NOAA_72640014839_{year}.csv',f'NOAA_72640594869_{year}.csv',f'POWER_solar_{year}.json'):
            path=ROOT/'C/data/raw'/name;meta=json.loads(Path(str(path)+'.request.json').read_text())
            if not meta.get('ok') or sha(path)!=meta['sha256']:raise RuntimeError('Raw cached source checksum failed: '+name)
            sources.append(dict(file=str(path.relative_to(ROOT)),sha256=sha(path),request=meta))
    inventory,chosen=candidates(frame);out.mkdir(parents=True,exist_ok=True)
    inventory.to_csv(out/'candidate_inventory.csv.gz',index=False,compression='gzip')
    annual=[]
    for year,g in frame.groupby(frame.timestamp_utc.dt.year):
        pool=inventory[pd.to_datetime(inventory.start_utc).dt.year==year]
        annual.append(dict(year=int(year),expected_hours=len(g),**{'valid_'+v:int(g[v].notna().sum()) for v in VARIABLES},
            complete_joint_hours=int(g[VARIABLES].notna().all(axis=1).sum()),nonpositive_temperature_hours=int((g.temperature_C<=0).sum()),
            known_max_hourly_mm=float(g.precipitation_mm_per_h.max()),known_max_24h_mm=float(g.precipitation_mm_per_h.rolling(24).sum().max()),
            candidate_windows=len(pool),complete_joint_windows=int(pool.complete_joint_window.sum()),eligible_complete_warm_windows=int(pool.eligible.sum())))
    pd.DataFrame(annual).to_csv(out/'annual_coverage.csv',index=False)
    events=[];baseline,comparison=benchmark(frame,out);events.append(baseline)
    lines=SOURCE.read_bytes().splitlines(keepends=True)
    if len(lines)!=len(frame)+1:raise RuntimeError('Cannot preserve source rows as exactCSV lines')
    for rank,row in enumerate(chosen,1):
        index=int(row.window_start_index);window=frame.iloc[index:index+144].copy();eid='historical_'+row.start_utc.strftime('%Y%m%dT%H%MZ')
        path=out/(eid+'.csv');path.write_bytes(lines[0]+b''.join(lines[index+1:index+145]))
        saved=pd.read_csv(path,float_precision='round_trip');saved.timestamp_utc=pd.to_datetime(saved.timestamp_utc,utc=True)
        pd.testing.assert_frame_equal(saved.reset_index(drop=True),window.reset_index(drop=True),check_exact=True,check_dtype=False)
        native=native_audit(window,eid,out)
        counts={name:int(window[name].sum()) for name in ('precipitation_trace_lower_bound','temperature_interpolated','wind_interpolated','temperature_short_gap_filled','wind_short_gap_filled')}
        record=dict(event_id=eid,scenario_id=eid,benchmark=False,start_utc=row.start_utc,end_utc=row.end_utc,
            forcing_file=str(path),forcing_sha256=sha(path),seed_sequence_entropy=None,hours=144,complete_forcing=True,
            forcing_schema='Original StageC joint hourly rows, UTC interval starts, all QC flags retained exactly',
            rainfall_total_mm=float(window.precipitation_mm_per_h.sum()),max_hourly_mm=float(window.precipitation_mm_per_h.max()),
            max_24h_mm=float(window.precipitation_mm_per_h.rolling(24).sum().max()),ranking_24h_mm=float(row.ranking_24h_mm),
            ranking_core_start_utc=row.ranking_core_start_utc,ranking_core_end_utc=row.ranking_core_end_utc,
            rainfall_rank=int(row['rank']),selected_additional_rank=rank,minimum_temperature_C=float(window.temperature_C.min()),
            solar_energy_Wh_m2=float(window.solar_W_per_m2.sum()),source='Observed KMKE rainfall and KMWC meteorology with POWER grid solar; existing StageC conservative hourly alignment',
            qc_counts=counts,native_rain_audit=native,provenance=dict(standardized_file=str(SOURCE.relative_to(ROOT)),standardized_sha256=sha(SOURCE),
                first_source_row=index,last_source_row=index+143,row_preservation='Original CSV header and144 row byte strings copied without rewriting values'))
        record.update(wet_spells(window));events.append(record)
    # A compact display of the highest distinct screening windows supplements the
    # full inventory. It is not a second, different selection algorithm.
    screen=[]
    for _,row in inventory[inventory.ranking_24h_mm.notna()].iterrows():
        if any(row.start_utc<old.end_utc and row.end_utc>old.start_utc for old in screen):continue
        screen.append(row)
        if len(screen)==25:break
    pd.DataFrame(screen).to_csv(out/'highest_nonoverlap_screening_windows.csv',index=False)
    selected=[]
    for event in events:
        selected.append({k:event.get(k) for k in ('event_id','benchmark','start_utc','end_utc','rainfall_total_mm','max_hourly_mm','max_24h_mm','minimum_temperature_C','positive_rain_hours','first_to_last_positive_span_h','solar_energy_Wh_m2','forcing_file','forcing_sha256')})
    pd.DataFrame(selected).to_csv(out/'selected_events.csv',index=False)
    manifest=dict(events=events,total_events=5,benchmark_count=1,additional_events=4,
        selection_rule='Reserve original Sep9–15,2022. Enumerate all144h UTC windows from2019–2024; rank the fixed24h core at hours24..47 by quantified precipitation descending, earliest core UTC start on exact ties. Require all144 original joint weather rows finite, full rain coverage and T>0. Greedily select4 eligible windows not overlapping the benchmark or each other. No ODE result enters selection.',
        ranking_is_rainfall_only=True,complete_joint_and_liquid_domain_are_eligibility_filters=True,
        candidate_counts=dict(total=len(inventory),valid24h_core=int(inventory.ranking_24h_mm.notna().sum()),
            complete_joint=int(inventory.complete_joint_window.sum()),eligible_complete_warm=int(inventory.eligible.sum()),
            overlapping_benchmark=int(inventory.overlaps_reserved_benchmark.sum())),
        standardized_source=dict(file=str(SOURCE.relative_to(ROOT)),sha256=sha(SOURCE),rows=len(frame),
            quality_report_sha256=sha(ROOT/'C/data/processed/data_quality_report.json'),processing_code_sha256=sha(ROOT/'C/uncertainty/data.py')),
        selection_code_sha256=sha(__file__),raw_source_manifest=sources,
        important_limitations=['Selection favors well-observed, nonfreezing windows; it does not rank unobserved or freezing storms.',
            'Rainfall is a regional airport proxy and hourly clock-bin amounts assume constant rates within observed native1h accumulations.',
            'Trace zeros are quantified lower bounds; six-zero-hour spell separation is not proof of meteorological dry periods.',
            'Temperature/wind are midpoint interpolation from accepted observations, including flagged gaps up to2h; solar is grid-derived with inherited hour-start convention.',
            'No joint weather vector is independently shuffled or scaled. New events retain source rows; benchmark retains its different approved120s reconstruction.',
            'All four selected new windows have smaller24h totals than the September2022 benchmark; this is not a more-extreme-than-benchmark or return-period study.',
            'The same model initial conditions are assumptions, not reconstructed antecedent storages at each historical date.'],
        benchmark_comparison=comparison,preservation_before=guard)
    write(out/'manifest.json',manifest)
    after=preservation.verify();write(out/'preservation_check.json',after)
    if not after['passed']:raise RuntimeError('Protected files changed during event preparation')
    print(json.dumps(clean(dict(total_events=5,selected=selected,candidate_counts=manifest['candidate_counts'],preservation=after)),indent=2))


if __name__=='__main__':main()
