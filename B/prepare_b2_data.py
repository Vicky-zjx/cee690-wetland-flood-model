#!/usr/bin/env python3
"""Small B.2 evidence builder and conservative rainfall utility; no downloads.

UTC offset-aware times and complete QC-approved intervals are mandatory.
Resampled interval means are derived estimates, never additional observations.
This does NOT silently turn interval means into ode4u's linear node values.
"""
from datetime import datetime, timezone
import argparse
import ast
import csv
import json
from pathlib import Path
import numpy as np
from A import wetland_flood_model as a2
ROOT=Path(__file__).resolve().parent


def utc_seconds(value):
    dt=datetime.fromisoformat(value.replace('Z','+00:00'))
    if dt.tzinfo is None or dt.utcoffset() is None:
        raise ValueError('Explicit timezone/UTC offset required; ambiguous local/DST time is not inferred.')
    return dt.astimezone(timezone.utc).timestamp()


def resample_rainfall(intervals,edges):
    """Overlap-weight depths assuming uniform rate WITHIN each source interval.

    intervals: [(start_seconds,end_seconds,depth_mm,quality_ok), ...].
    Returns target interval depths [mm]. Gaps, overlap, invalid QC, missing and
    negative values fail. Partial source intervals are proportional assumptions.
    """
    edges=np.asarray(edges,dtype=float)
    if len(edges)<2 or not np.all(np.isfinite(edges)) or np.any(np.diff(edges)<=0): raise ValueError('Invalid target interval edges.')
    records=sorted(intervals,key=lambda r:r[0]);out=np.zeros(len(edges)-1);coverage=np.zeros_like(out)
    previous_end=None
    for start,end,depth,quality_ok in records:
        if not all(np.isfinite(v) for v in (start,end,depth)) or end<=start or depth<0 or quality_ok is not True:
            raise ValueError('Rainfall interval is missing, negative, bad-QC or invalid; never replace with zero.')
        if previous_end is not None and start<previous_end: raise ValueError('Overlapping rainfall intervals double count accumulations.')
        previous_end=end
        overlap=np.maximum(0,np.minimum(end,edges[1:])-np.maximum(start,edges[:-1]))
        out+=depth*overlap/(end-start);coverage+=overlap
    if not np.allclose(coverage,np.diff(edges),rtol=0,atol=1e-7): raise ValueError('Unobserved rainfall interval; no zero filling or implicit interpolation.')
    return out


def reconcile_head(gage_height_m,datum_offset_m):
    if datum_offset_m is None: raise ValueError('Gage height is not elevation; surveyed datum offset required.')
    values=np.asarray(gage_height_m,dtype=float)
    if not np.all(np.isfinite(values)) or not np.isfinite(datum_offset_m): raise ValueError('Invalid head or datum.')
    return values+datum_offset_m


def stringify(value):
    if isinstance(value,np.ndarray):value=value.tolist()
    return json.dumps(value,default=lambda v:v.tolist() if isinstance(v,np.ndarray) else v,ensure_ascii=False)


def evidence_register():
    params=a2.c();catalog=json.loads((ROOT/'source_catalog.json').read_text());source={row['id']:row for row in catalog}
    cfg=json.loads((ROOT/'b2_config.json').read_text())
    numeric={'hour','day','duration','time_step','balance_rtol','gate1_hours','gate2_hours'}
    physical={'lambda_v':('FAO56','Ch3 latent heat; 2.45 MJ/kg times 1e6'),
              'c_p':('FAO56','Ch3 Eq8: 1.013e-3 MJ/(kg degC) times 1e6; K and degC intervals equal'),
              'epsilon':('FAO56','Ch3 Eq8 molecular mass ratio .622'),
              'rho_w':('FAO56','Ch3 Box8 rho_w=1000 kg/m3 approximation'),
              'R_v':('UCAR_JRA25_CONSTANTS','Gas constant for water vapor row: 461.5 J/(kg K)')}
    unit={}
    for units,names in {
        's':'hour day duration time_step', 'm':'H_flood H_r_valid_max H_g_valid_max soil_bottom z_d z_w z_r z_g z_G z_dsp z_wo z_wsp z_ro u_depth_init w_depth_init d_depth_init r_depth_init H_g_init h_d_base h_g_base h_g_amplitude',
        'm2':'A_u A_d A_w A_r A_g', 'm3':'C_m S_fc V_u V_m V_w V_d V_g V_r',
        'm/s':'K_u k_0 wind_mean wind_amplitude', '1/s':'k_mg k_ud k_ur',
        'm3/s':'Q_udcap Q_urcap Q_Gcap Q_wrcap Q_rcap Q_ext_base',
        's/m2':'R_G R_bare R_veg R_r', 'm^(3/2)/s':'c_dsp c_wsp', 'm2/s':'C_wg C_gr C_gx',
        'J/kg':'lambda_v', 'J/(kg K)':'c_p R_v', 'kg/m3':'rho_w', 'Pa':'p_0',
        'h':'weather_period_h temperature_phase_h sunrise_phase_h wind_period_h groundwater_period_h gate1_hours gate2_hours',
        'degC':'T_mean T_amplitude', 'W/m2':'R_s_peak', 'category':'tau',
        '(h,h,h,mm/h)':'rain_pulses','(h,h,h,m3/s)':'ext_pulses','(h,h,h,m)':'h_d_pulse',
        '1':'S_y balance_rtol K_bare K_veg rho epsilon r_H eta k_v K_ETm K_ETd K_ETr soil_fraction_init B_base gate1_values gate2_values'
    }.items():
        for name in names.split(): unit[name]=units
    tree=ast.parse((a2.ROOT/'wetland_flood_model.py').read_text());locations={}
    for node in ast.walk(tree):
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Attribute) and isinstance(target.value,ast.Name) and target.value.id=='self':locations.setdefault(target.attr,node.lineno)
    rows=[]
    def add(identifier,name,value,units,category,sid,locator,derivation,rationale,limits,status,old=None,new=None,station='not selected',version=None,potential=''):
        meta=source.get(sid,{})
        rows.append(dict(item_id=identifier,symbol_code_name=name,value_range=stringify(value),units=units,evidence_category=category,
            source_id=sid,original_field_table_page_equation=locator,dataset_version=version or meta.get('date_or_version','local 82e9a78'),
            station_polygon_dates=station,transformation_derivation=derivation,rationale=rationale,applicability_limitations=limits,status=status,
            A2_old_value=stringify(old),B2_new_value=stringify(new),potential_source_ids=potential))
    for name,value in vars(params).items():
        if name not in unit: raise ValueError('Missing unit in evidence register: '+name)
        category='unsupported placeholder';sid='A2_SOURCE';locator=f'class c, line {locations[name]}'
        derivation='Existing illustrative value; no local measurement or calibrated derivation.'
        status='implemented_A2; NOT approved/used in project B2'
        limitation='A local source portal or a mechanism citation does not support this numeric value.'
        if name in numeric:
            category='numerical setting';status='A2 setting retained as explicit B2 numerical proposal';limitation='Not observed operating history, physical freeboard, or evidence of accuracy.'
            derivation='Exact unit conversion or documented numerical horizon, grid, knot times or tolerance.'
        if name in physical:
            category='physical relation';sid,derivation=physical[name];locator=derivation
            status='source numeric relation verified; implemented A2; inherited B2 approximation pending configuration'
            limitation='Supports only this constant; not calibration or validation of reduced ET.'
        potential=''
        if name in {'A_u','A_d','A_w','A_r','z_d','z_w','z_r','z_g','z_dsp','z_wsp','z_wo','z_G','soil_bottom'}:potential='USGS_WBD;WIDNR_WWI;USGS_SEWRPC_LIDAR2015'
        if name in {'C_m','S_fc','K_u','S_y'}:potential='USDA_SSURGO;USDA_SSURGO_FIELDS'
        if name in {'R_veg','K_veg','tau','rho'}:potential='WIDNR_WWI_CLASS'
        if name in {'rain_pulses'}:potential='MMSD_RAIN'
        if name in {'T_mean','T_amplitude','wind_mean','wind_amplitude','p_0','r_H'}:potential='NOAA_ISD'
        if name=='R_s_peak':potential='NASA_NLDAS2'
        if name in {'Q_ext_base','ext_pulses','h_d_base','h_d_pulse'}:potential='USGS_04087088'
        if name in {'h_g_base','h_g_amplitude'}:
            limitation='Requires a hydraulically connected shallow-well record with surveyed common datum; the candidate streamgage cannot establish groundwater head.'
        if name=='B_base':
            limitation='Requires outlet inspection or a justified event-specific engineering blockage assumption; streamgage discharge cannot directly identify blockage.'
        new=value
        if category=='unsupported placeholder':
            category='proposed engineering assumption';sid='B2_PROPOSAL'
            derivation='Retain this explicit uncalibrated A.2 value for the proposed 0.66km2 hypothetical installation; no local derivation.'
            status='concrete proposal; approval pending; not used in B2 simulation'
            limitation+=' This is a selected exploratory assumption, not a field-supported property.'
        unused={'R_G','Q_Gcap','rain_pulses','ext_pulses','Q_ext_base','weather_period_h','temperature_phase_h','sunrise_phase_h','wind_period_h','groundwater_period_h','T_mean','T_amplitude','R_s_peak','wind_mean','wind_amplitude','h_d_base','h_d_pulse','h_g_base','h_g_amplitude','B_base'}
        if name in unused:
            new=None;status='A2 reproducibility only; NOT used in proposed B2 revised gates / prepared forcing'
            derivation='Replaced by explicit gate parameters or prepared event forcing and separately declared constant boundaries.'
        if name=='H_flood':new=None;status='Optional flood reporting disabled; not a hard constraint'
        if name=='z_G':
            new=cfg['gates']['gate_sill_z'];sid='B2_PROPOSAL';derivation='New hypothetical bed-mounted sill used for geometry checks; old threshold not interpreted as surveyed sill.'
        if name in cfg['parameter_overrides']:
            rec=cfg['parameter_overrides'][name];new=rec['value'];sid=rec['source_id'];category=rec.get('evidence_category','proposed engineering assumption');derivation=rec['reason'];status='derived/proposed fixed event setting; approval pending';locator=rec.get('locator',locator);limitation='Event-time airport observation-derived constant; spatial transfer and a fixed temporal approximation remain unapproved, not calibration.'
        add('A2-'+name,name,value,unit[name],category,sid,locator,derivation,'Fixed before comparison; same parameter in revised reference and all candidates.',limitation,status,old=value,new=new,potential=potential)
    for name,value,units,formula in [('A_rain',params.A_rain,'m2','A_u+sum(A_d)+sum(A_w)+A_r'),('C_d',params.C_d,'m3','A_d*(z_dsp-z_d); spill onset, NOT safety cap'),('C_w',params.C_w,'m3','A_w*(z_wsp-z_w); spill onset, NOT safety cap'),('x_init',params.x_init,'m3','Initial depths, soil fraction and groundwater level converted using fixed geometry'),('R_w',params.design()[0],'s/m2','(1-rho)*R_bare+rho*R_veg[tau]'),('K_ETw',params.design()[1],'1','(1-rho)*K_bare+rho*K_veg[tau]')]:
        add('DER-'+name,name,value,units,'proposed engineering assumption','A2_SOURCE','class c properties',formula,'Derived from synthetic A.2 quantities, NOT locally derived data.','Inherits all unsupported geometry/design inputs.','explicit hypothetical design/initial-state proposal; approval pending',old=value,new=value)
    for identifier,name,value,category,sid,locator,limits in [
        ('REL-es','saturation pressure constants',[610.8,17.27,237.3],'physical relation','FAO56','Ch3 Eq11; 0.6108kPa to610.8Pa','Only vapour-pressure subrelationship; empirical valid-temperature scope.'),
        ('STR-balance','eight ODE balances',None,'physical relation','B1','p2 constraint C and A2 system','Conservation does not validate closures.'),
        ('STR-branches','two hypothetical dam-gate-wetland branches',2,'proposed engineering assumption','B1','p1 design decisions; existing model topology','Not the as-built Underwood/County Grounds layout.'),
        ('STR-areas','fixed effective surface areas',None,'proposed engineering assumption','A2_SOURCE','stage()','Wetland polygons are not stage-storage curves; DEM is not bathymetry.'),
        ('STR-oneway','one-way internal gate and surface exchange',None,'proposed engineering assumption','A2_SOURCE','head_flow and system','Ignores possible actual backflow; h_d may suppress forward flow.'),
        ('STR-ET','reduced Penman-type closure',None,'proposed engineering assumption','A2_SOURCE','evaporative_demand; compare FAO56 only for subterms','Not FAO reference ET0; fixed humidity and transfer coefficients unsupported.'),
        ('STR-design','linear bare/vegetation endpoint interpolation',None,'proposed engineering assumption','A2_SOURCE','class c.design','WWI category cannot establish R or ET endpoint; Manning n is not R_w.'),
        ('STR-rain','six disjoint rain footprints; no duplicate soil/aquifer rainfall',None,'proposed engineering assumption','A2_SOURCE','A_rain and system','Must spatially delineate footprints and external inflow origin.'),
        ('STR-groundwater','lumped shallow groundwater with fixed specific yield',None,'proposed engineering assumption','A2_SOURCE','stage and Qgx','Deep confined well observations cannot substitute for shallow head without evidence.')]:
        units='[Pa, 1, degC]' if identifier=='REL-es' else 'see definition'
        add(identifier,name,value,units,category,sid,locator,'No newly fitted relation.','Model structure stated transparently.',limits,'current structural assumption; local applicability unverified')
    cfg=json.loads((ROOT/'b2_config.json').read_text())
    for group in ['dam','wetland']:
        for i in range(2):
            name=f'H_{group}{i+1}_safe'
            value=cfg['safety'][group+'_m'][i]
            add('SAFE-'+name,name,value,'m assumed model datum','proposed engineering assumption','B2_PROPOSAL','Four proposed maximum levels; WI NRCS378 printed p5 freeboard analogue',cfg['safety']['derivation'],'Fixed containment constraint independent of simulated response.','Hypothetical settled top, not surveyed; no ecological or real facility safety certification.','specific proposal; approval pending',new=value,potential='WI_NRCS378_2024')
    for name,units in [('gravity','m/s2'),('gate_width','m'),('gate_max_opening','m'),('gate_sill_z','m assumed model datum'),('Cd_gate_free','1'),('Cd_gate_submerged','1'),('neglect_approach_velocity_head','Boolean')]:
        sid=cfg['gates']['source_ids'][name]
        category='physical relation' if name=='gravity' else 'external analogue' if name.startswith('Cd_') else 'proposed engineering assumption'
        value=cfg['gates'][name]
        locator='NIST standard gravity' if name=='gravity' else 'USACE free sluice typical0.5–0.7; proposed midpoint0.60' if name=='Cd_gate_free' else 'USACE fully submerged Eq3; typical0.80, no verified interval' if name=='Cd_gate_submerged' else 'b2/docs/01_study_configuration.md and b2/docs/06_gate_proposal.md fixed gates'
        add('GATE-'+name,name,value,units,category,sid,locator,'Dimensions based on hypothetical pool-area ratio and bed-to-spill depth; no fitting or local survey.','Explicit fixed gate parameter; never optimized.','Reference coefficients are not local calibration; approach-head and limited-domain assumptions require approval.','specific proposal; approval pending',new=value)
    for name,value,sid,locator in [('gate_flow_H_over_a',1.25,'G4','p8-7; gate domain'),('submergence_lower',.67,'G1','Submergence paragraph'),('submergence_upper',.80,'G2','Transition to full submergence'),('submergence_blend',None,'G9','Explicit fraction and weighted law; reviewed complete page')]:
        add('GATE-'+name,name,value,'1','external analogue',sid,locator,'Use stated criteria with head conventions.','Implement documented gate regimes.','Rounded .67 yields source-level ~0.501% lower endpoint jump; low-head/weir regimes not supplied.','implemented; fixture-tested; local applicability unapproved')
    for name,value in cfg['numerics'].items():
        units='m3' if name=='volume_scales_m3' else 'm' if name.endswith('_m') else 's' if name.endswith('_s') or name=='total_seconds' else '1'
        add('NUM-'+name,name,value,units,'numerical setting','B2_CONFIG','numerics.'+name,'Fixed before optimization; not a physical observation.','Conditioning, finite failure response, grid checking or run budget.','Strict margin is numerical, not physical freeboard; grid checks do not prove intersample feasibility.','proposed numeric setting; optimization not executed',new=value)
    hyp_index={'display':0,'tol_x':1,'tol_f':2,'feasibility_tol':3,'max_evals':4,'penalty':5,'penalty_exponent':6}
    for name,value in cfg['optimizer'].items():
        mapping=('Adapter wall-clock check before each callback; NOT an official NMS hyperparameter.'
                 if name=='max_seconds' else f'Local choice mapped by b2/official_opt.py run_opt to official hyp[{hyp_index[name]}].')
        if name=='max_evals':mapping+=' Adapter also enforces a hard callback cap because native simplex operations may overshoot.'
        add('OPT-'+name,name,value,'s' if name=='max_seconds' else '1','numerical setting','B2_CONFIG','optimizer.'+name+'; b2/official_opt.py run_opt',mapping,'Bounded derivative-free search of nonsmooth grid extrema.','Penalty does not establish feasibility; n^2 iteration requirement may preclude convergence under budget.','adapter tested on official example only',new=value,potential='OFFICIAL_OPT')
    add('CTX-gage_area','USGS04087088 drainage area',46.879,'km2','locally derived','USGS_04087088','Drainage area18.1 square miles','18.1*2.589988110336 = 46.8787847970816; rounded to 46.879','Check scale compatibility, NOT change A2 model area.','Whole HUC12 != gage catchment; 0.66km2 A2 footprint is much smaller.','metadata reviewed; NOT model input',station='USGS04087088; metadata review2026-10-01')
    add('CTX-huc_area','Underwood HUC12 polygon area',50.72,'km2','locally derived','USGS_WBD_UNDERWOOD','huc12=040400030404; areasqkm','Official WBD attribute, no outlet catchment delineation.','Provisional geographic context only.','Do not substitute whole HUC12 for gage catchment.','small attribute downloaded; NOT model input',station='040400030404; attribute query2026-10-01')
    # Actual preparation evidence; these arrays have not been used in a model run.
    prepared_path=ROOT/'work/event_data/proposed_forcing_manifest.json'
    if prepared_path.exists():
        prepared=json.loads(prepared_path.read_text())
        for name,value,units,category,sid,locator,derivation,limits in [
            ('P(t)',prepared['rain_depth_m'],'m event accumulation','local observation','NOAA_ISD','KMKE AA1 selected nonoverlapping1h/3h intervals','Quantified145.3mm converted to SI; each interval integral preserved; trace0 lower-bound proposal.','Airport16km from gauge; not measured arealrain; trace amount unquantified.'),
            ('T_a(t)','prepared4321-node series','degC','local observation','NOAA_ISD','KMWC TMP; observationUTC CSV','QC-accepted reports linearly sampled to120s; brackets<=70min, no extrapolation.','Interpolation is numerical, not additional observations; spatial proxy.'),
            ('v(t)','prepared4321-node series','m/s','local observation','NOAA_ISD','KMWC WND speed','QC-accepted reported wind linearly sampled to120s.','Sensorheight not established; no height correction or local exposure calibration.'),
            ('R_s(t)',prepared['solar_energy_Wh_m2'],'Wh/m2 event energy','local gridded/reanalysis estimate','NASA_POWER_EVENT','ALLSKY_SFC_SW_DWN; SYN1DEG UTC hourlylabels','Wh/m2 perhour to meanW/m2, energy-preserving120s nodal reconstruction.','Coarse~1degree satellitegrid, not local measurement; hour-start convention unverified and unapproved.'),
            ('Q_ext(t)',0.,'m3/s','proposed engineering assumption','B2_PROPOSAL','One mapping: closed surface drainage unit','No external surface contributing area; gaugeflow excluded.','Applies only to hypotheticalinstallation, not actualUnderwoodmainstem.'),
            ('h_d(t)',97.15,'m assumed model datum','proposed engineering assumption','B2_PROPOSAL','One mapping: hypothetical downstream boundary','Fixed modeldatum head, no syntheticpulse.','Not NAVD88 or gaugestage; no surveyed downstream boundary.'),
            ('h_g(t)',97.15,'m assumed model datum','proposed engineering assumption','B2_PROPOSAL','One mapping: hypothetical shallow-aquifer boundary','Fixed exteriorhead; no syntheticoscillation.','No connected shallowwell observation; actualgroundwater unknown.'),
            ('B(t)',0.,'1','proposed engineering assumption','B2_PROPOSAL','One mapping: maintained unobstructed outlet','No blocked fraction in this specified clean-outlet design condition.','No eventinspection supports actual blockage; unapproved design assumption.')]:
            add('INPUT-'+name,name,value,units,category,sid,locator,derivation,'Concrete proposed fixed input, same for reference/candidates.',limits,'prepared/proposed; NOT approved or simulated',new=value,station='Event2022-09-09T00Z through2022-09-15T00Z; antecedentSep2 onward')
    # Approval is a study decision, not evidence of observation or calibration.
    case_approved=cfg.get('mode')=='approved_illustrative' and all(cfg[k].get('approved') for k in
        ['geography','safety','gates','categories','fixed_model','forcing','initial_conditions'])
    if case_approved:
        run_status_path=ROOT/'results/run_status.json'
        execution=json.loads(run_status_path.read_text()).get('status','not yet run') if run_status_path.exists() else 'not yet run'
        for row in rows:
            row['approval_scope']='First exploratory hypothetical installation only'
            if row['item_id'].startswith('CTX-'):
                continue
            if row['B2_new_value'] in ('null','None','') and 'NOT used' in row['status']:
                continue
            if row['symbol_code_name']=='H_flood':
                continue
            row['status']='Approved/declared for first exploratory case; execution: '+execution
            row['applicability_limitations']=row['applicability_limitations'].replace('unapproved','approved for this case only').replace('remain unapproved','are approved for this case only')
            if row['item_id'].startswith('SAFE-'):
                row['original_field_table_page_equation']='Assumed crests in b2_config.json; NRCS reference supports freeboard concept only'
                row['applicability_limitations']='Proposed embanked/contained units. All elevations are assumed, not supplied or certified by NRCS; no local engineering certification.'
            if row['symbol_code_name'] in ['gate1_values','gate2_values','rho','tau']:
                row['rationale']='Control/design decisions vary between reference and optimizer candidates; see decision_values.csv.'
                row['status']='Original A.2 decision values retained for provenance; executed replacement and optimized decisions are recorded separately.'
            if row['item_id']=='OPT-max_evals':
                row['transformation_derivation']+=' First exploratory case allocates this180-callback total per category as40 full-box plus140 local-box callbacks; exact numerical boxes in run_config/search and each optimizer JSON.'
    path=ROOT/'tables/evidence_register.csv';path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    # Every A2 attribute must be represented; all IDs must resolve to catalog entries.
    assert {r['symbol_code_name'] for r in rows if r['item_id'].startswith('A2-')}==set(vars(params))
    unknown={r['source_id'] for r in rows}-set(source)
    if unknown:raise ValueError('Missing source IDs: '+str(unknown))
    print(f'Wrote {len(rows)} evidence rows; covers every class c attribute plus derived, structural, gate, safety and numerical items.')


def interval_rates_to_nodes(intervals, times):
    """Conserve each interval integral with nonnegative linear node values.

    Each tuple is start,end,mean_rate in UTC epoch seconds. Interior plateau
    and one-grid-step ramps meet the smaller adjoining rate at a boundary.
    This is a declared within-interval numerical reconstruction, not observed
    subinterval variability. Unlike interpolation of means, dry bins stay dry.
    """
    dt=float(times[1]-times[0]);out=np.full(len(times),np.nan);target=0.0
    parts=[(max(a,times[0]),min(b,times[-1]),r) for a,b,r in intervals if b>times[0] and a<times[-1]]
    if not parts or parts[0][0]!=times[0] or parts[-1][1]!=times[-1]:raise ValueError('Source intervals do not bracket the requested event.')
    for i,(a,b,r) in enumerate(parts):
        if not np.isfinite(r) or r<0 or b-a<2*dt:raise ValueError('Invalid interval rate or interval too short for declared reconstruction.')
        if i and a!=parts[i-1][1]:raise ValueError('Source interval gap/overlap; no implicit zero fill.')
        start,end=(np.array([a,b])-times[0])/dt
        if not np.allclose([start,end],np.round([start,end]),atol=1e-7,rtol=0):raise ValueError('Source interval endpoints do not align with model grid.')
        j,k=int(round(start)),int(round(end))
        left=r if i==0 else min(parts[i-1][2],r)
        right=r if i==len(parts)-1 else min(r,parts[i+1][2])
        plateau=(r*(b-a)-dt*(left+right)/2)/((b-a)-dt)
        out[j:k+1]=plateau;out[j]=left;out[k]=right
        if not np.isclose(np.trapezoid(out[j:k+1],times[j:k+1]),r*(b-a),rtol=1e-12,atol=1e-10):raise ValueError('Interval mass/energy reconstruction failed.')
        target+=r*(b-a)
    if not np.all(np.isfinite(out)) or np.any(out<0):raise ValueError('Missing or negative reconstructed input.')
    assert np.isclose(np.trapezoid(out,times),target,rtol=1e-12,atol=1e-10)
    return out,target


def prepare_event_forcing():
    """Build this event's proposed forcing from downloaded, QC-processed CSVs.

    No network, simulation or implicit approval. Raw retrieval and QC provenance
    remain in work/event_data and the accompanying findings documents.
    """
    import hashlib
    cfg=json.loads((ROOT/'b2_config.json').read_text())
    base=ROOT/'work/event_data/weather'
    def read(name):
        with (base/name).open() as f:return list(csv.DictReader(f))
    start=utc_seconds(cfg['geography']['event_start_utc']);end=utc_seconds(cfg['geography']['event_end_utc'])
    dt=cfg['numerics']['dt_s'];times=np.arange(start,end+dt/2,dt)
    if end-start!=cfg['numerics']['duration_s']:raise ValueError('Event and fixed six-day model horizon disagree.')
    rain_rows=read('KMKE_rain_selected_intervals_utc.csv')
    intervals=[];trace_count=0
    for row in rain_rows:
        a=utc_seconds(row['interval_start_utc']);b=utc_seconds(row['interval_end_utc'])
        if row['status'] not in {'accepted','accepted_trace_lower_bound_proposal'} and not row['status'].startswith('accepted'):raise ValueError('Rainfall QC not accepted: '+row['status'])
        depth=float(row['depth_mm_quantified'])
        intervals.append((a,b,depth*.001/(b-a)))
        if b>start and a<end and row['trace']=='True':trace_count+=1
    rain,rain_depth=interval_rates_to_nodes(intervals,times)
    met=read('KMWC_meteorology_observations_utc.csv');met_t=np.array([utc_seconds(r['timestamp_utc']) for r in met])
    forcing={}
    for field,name in [('temperature_C','T_a'),('wind_speed_m_s','v')]:
        keep=np.array([bool(row[field]) for row in met]);tt=met_t[keep];vv=np.array([float(row[field]) for row in met if row[field]])
        j=max(0,int(np.searchsorted(tt,start))-1);k=min(len(tt),int(np.searchsorted(tt,end))+1)
        if tt[j]>start or tt[k-1]<end or np.max(np.diff(tt[j:k]))>70*60:raise ValueError('Weather lacks accepted bracketing observations or exceeds70minute interpolation proposal.')
        forcing[name]=np.interp(times,tt,vv)
    solar=read('NASA_POWER_solar_hourly_utc.csv')
    radiation,solar_energy=interval_rates_to_nodes([(utc_seconds(r['source_hour_label_utc']),utc_seconds(r['source_hour_label_utc'])+3600,float(r['hourly_mean_W_m2'])) for r in solar],times)
    bounds=cfg['forcing']['proposed_fixed_boundaries']
    values=np.column_stack([times-start,rain,forcing['T_a'],radiation,forcing['v'],np.full(len(times),bounds['Q_ext_m3_s']),np.full(len(times),bounds['h_d_m']),np.full(len(times),bounds['h_g_m']),np.full(len(times),bounds['B'])])
    if not np.all(np.isfinite(values)):raise ValueError('Proposed forcing contains gaps/nonfinite entries.')
    out=ROOT/'work/event_data/proposed_forcing_120s.csv'
    np.savetxt(out,values,delimiter=',',header='time_s,P,T_a,R_s,v,Q_ext,h_d,h_g,B',comments='',fmt='%.12g')
    loaded=np.genfromtxt(out,delimiter=',',names=True)
    if not np.isclose(np.trapezoid(loaded['P'],loaded['time_s']),rain_depth,rtol=0,atol=1e-10):raise ValueError('Written rainfall CSV lost accumulated depth.')
    previous_path=ROOT/'work/event_data/proposed_forcing_manifest.json'
    previous=json.loads(previous_path.read_text()) if previous_path.exists() else {}
    same_approved=(cfg['forcing'].get('approved') is True and previous.get('approved') is True
        and previous.get('sha256')==hashlib.sha256(out.read_bytes()).hexdigest()
        and previous.get('event_start_utc')==cfg['geography']['event_start_utc']
        and previous.get('event_end_utc')==cfg['geography']['event_end_utc'])
    meta=dict(approved=same_approved,status=('Identical previously approved forcing regenerated for first exploratory case only' if same_approved else 'Proposed forcing prepared; changed or new inputs require explicit approval'),
      sha256=hashlib.sha256(out.read_bytes()).hexdigest(),source_ids=['NOAA_ISD','NASA_POWER_EVENT','B2_PROPOSAL'],
      event_start_utc=cfg['geography']['event_start_utc'],event_end_utc=cfg['geography']['event_end_utc'],timezone='UTC',
      original_timezone='NOAA ISD DATE is UTC; POWER requested time-standard=UTC; USGS localoffsets retained separately and not used as forcing',
      timestamp_convention='Rain AA periods end at actual report times; radiation UTC labels provisionally denote start of hour; T and wind instantaneous-report linear interpolation',
      qc_summary='NOAA processed accepted records only; KMKE two unavailable hourly amounts represented by retained nonoverlapping observed3hour TRACE interval; subinterval amounts remain unresolved; no bad/missing value becomes an observed zero',
      gap_policy='No missing-data zero filling. Trace depth uses quantified lower bound0 as an explicit processing assumption governed by the manifest approval flag. Meteorology interpolation gaps <=70min. MMSD historicalresponse unsuitable and excluded.',
      derivations=['Rain mm/interval converted to m/s; intervals clipped by proportional within-interval assumption; exactly preserved each represented interval depth.',
       'T and reported airport wind interpolated to120s model knots; no unsupported wind-height conversion; wind exposure remains an approximation in the empirical ET closure.',
       'POWER Wh/m2 per hour converted to meanW/m2; coarse SYN1DEG estimate, not a local radiometer.',
       'Rain/radiation node reconstruction: endpoints=min(adjacent mean rates), interior plateau adjusted to preserve integral,120s linear ramps. Dry/no-energy intervals remain zero. Not measured subhour variability.',
       'Refined ODE evaluations must reuse these same120s forcing knots, not reconstruct a new storm shape.'],
      vertical_datum=cfg['geography']['vertical_datum'],Q_ext_accounting='Q_ext=0: no external surface catchment in proposed hypothetical installation; USGS totalflow excluded',
      rain_depth_m=rain_depth,rain_mass_tolerance_m=1e-10,solar_energy_Wh_m2=solar_energy/3600,trace_intervals_in_event=trace_count,
      rows=len(times),dt_s=dt,local_observations_used_for_preparation=True,model_simulation_performed=False,
      source_files=[dict(path=str((base/n).relative_to(ROOT)),sha256=hashlib.sha256((base/n).read_bytes()).hexdigest()) for n in ['KMKE_rain_selected_intervals_utc.csv','KMWC_meteorology_observations_utc.csv','NASA_POWER_solar_hourly_utc.csv']])
    manifest=ROOT/'work/event_data/proposed_forcing_manifest.json';manifest.write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps(dict(rows=len(times),rain_mm=rain_depth*1000,solar_Wh_m2=solar_energy/3600,trace_intervals=trace_count,approved=False),indent=2))

def main():
    parser=argparse.ArgumentParser(description=__doc__);action=parser.add_mutually_exclusive_group(required=True);action.add_argument('--evidence-register',action='store_true');action.add_argument('--event-forcing',action='store_true')
    args=parser.parse_args()
    if args.evidence_register:evidence_register()
    elif args.event_forcing:prepare_event_forcing()
if __name__=='__main__':main()
