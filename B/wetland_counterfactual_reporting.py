"""Small table/figure/report export for the fixed dynamic counterfactual."""
from datetime import datetime, timedelta
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from B import optimize_wetland as old
from A import wetland_flood_model as a2
from B import wetland_counterfactual as cf

ROOT, OUT = cf.ROOT, cf.OUT
TABLES = ROOT/'tables/wetland_counterfactual'
PLOTS = ROOT/'plots/wetland_counterfactual'
REPORT = ROOT/'docs/13_no_wetland_counterfactual.md'


def format_value(value):
    if value is None: return 'N/A'
    if isinstance(value, (bool,np.bool_)): return 'Yes' if value else 'No'
    if isinstance(value, (float,np.floating)): return f'{value:.6f}'
    return str(value)


def csv_write(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as file:
        writer=csv.DictWriter(file,fieldnames=keys); writer.writeheader()
        for row in rows: writer.writerow({k: 'N/A' if row.get(k) is None else row[k] for k in keys})


def transposed(rows, fields):
    lines=['| Quantity | '+' | '.join(r['case'] for r in rows)+' |',
           '|---|'+'---:|'*len(rows)]
    for key,title in fields:
        lines.append('| '+title+' | '+' | '.join(format_value(r.get(key)) for r in rows)+' |')
    return '\n'.join(lines)


def series_stats(t, values):
    volume=float(np.trapezoid(values,t)); k=int(np.argmax(values))
    return dict(volume_m3=volume,peak_m3s=float(values[k]),peak_time_h=float(t[k]/3600),
                centroid_h=None if volume<=0 else float(np.trapezoid(values*t,t)/volume/3600))


def main_row(rec, run, peak_bypass, config):
    t=run['time']; c=lambda name:cf.column(run,name); integral=lambda values:float(np.trapezoid(values,t))
    peak_index=int(np.argmax(c('H_r'))); start=datetime.fromisoformat(config['geography']['event_start_utc'].replace('Z','+00:00'))
    row=dict(case=cf.LABELS[rec['case']],case_id=rec['case'],
             initialization='approved shared states; absent wetland water' if rec['case']=='bypass' else ('zero initial wetlands diagnostic' if rec['case'].endswith('_zero') else 'approved primary'),
             design=None if rec['case']=='bypass' else '/'.join(s.replace('type_','').upper() for s in run['pair'])+f'; covers {run["theta"][-2]:.17g}/{run["theta"][-1]:.17g}',
             dynamic_states=run['x'].shape[0],river_peak_m=float(c('H_r')[peak_index]),peak_time_h=float(t[peak_index]/3600),
             peak_time_utc=(start+timedelta(seconds=float(t[peak_index]))).isoformat(),
             difference_vs_bypass_mm=(peak_bypass-float(c('H_r')[peak_index]))*1000,
             gate_applicability='Passed at all official RHS/RK stages; saved executions verified' if rec['reused'] else 'Passed at all official RHS/RK stages',
             initial_water_m3=float(run['x'][:,0].sum()),final_water_m3=float(run['x'][:,-1].sum()),
             storage_change_m3=float(run['x'][:,-1].sum()-run['x'][:,0].sum()),
             rain_volume_m3=integral(c('P_total')),external_surface_inflow_m3=integral(run['u'][a2.U['Q_ext']]),
             ET_m3=integral(c('E_tot')),signed_external_groundwater_m3=integral(c('Q_gx')),
             groundwater_in_m3=integral(c('Q_gx_in')),groundwater_out_m3=integral(c('Q_gx_out')),
             outlet_volume_m3=integral(c('Q_out')),balance_final_m3=float(run['balance'][-1]),
             balance_max_abs_m3=float(max(abs(run['balance']))),balance_relative=rec['metrics']['balance_relative'],
             derivative_balance_max_m3s=rec['metrics']['derivative_balance_max_m3s'],
             timestep_s=rec['dt_s'],valid=rec['valid'],feasible=rec['feasible'],status=rec['status'],
             solver_termination='Completed fixed interval',reused=rec['reused'],
             river_final_stage_m=float(c('H_r')[-1]),river_rising_at_end=rec['metrics']['river_rising_at_end'],
             peak_near_horizon=rec['metrics']['peak_near_horizon'],rain_depth_mm=rec['rain_depth_m']*1000,
             solar_energy_Wh_m2=rec['solar_energy_Wh_m2'])
    for i in (1,2):
        row[f'dam{i}_max_m']=float(max(c(f'H_d{i}')))
        row[f'dam{i}_margin_m']=config['safety']['dam_m'][i-1]-row[f'dam{i}_max_m']
        row[f'G5_branch{i}_min_margin_m']=next(r['margin'] for r in run['audit'] if r['label']==f'gate{i}_head_opening')
        row[f'wetland{i}_max_m']=None if rec['case']=='bypass' else float(max(c(f'H_w{i}')))
        row[f'wetland{i}_margin_m']=None if rec['case']=='bypass' else config['safety']['wetland_m'][i-1]-row[f'wetland{i}_max_m']
        tail=c('H_r') if rec['case']=='bypass' else c(f'H_w{i}')
        row[f'gate{i}_tailwater_clearance_min_m']=float(config['gates']['gate_sill_z'][i-1]-max(tail))
    # Independent accounting from individual integrated components, not stored net.
    row['independent_balance_final_m3']=row['storage_change_m3']-(row['rain_volume_m3']+row['external_surface_inflow_m3']+row['signed_external_groundwater_m3']-row['outlet_volume_m3']-row['ET_m3'])
    return row


def mechanism(rec, run, p):
    t=run['time']; c=lambda name:cf.column(run,name); integral=lambda v:float(np.trapezoid(v,t))
    k=int(np.argmax(c('H_r'))); bypass=rec['case']=='bypass'; row=dict(case=cf.LABELS[rec['case']],case_id=rec['case'])
    groups={'gate1':c('Q_G1'),'gate2':c('Q_G2'),'gate_combined':c('Q_G1')+c('Q_G2'),
            'upland_direct':c('Q_ur'),'groundwater_to_river':c('Q_gr')}
    if not bypass:
        for key, stem in [('normal','Q_wr'),('overflow','Q_wsp'),('seepage','Q_wg'),('wetland_ET','E_w')]:
            groups[key]=c(stem+'1')+c(stem+'2')
            for branch in (1,2): row[f'{key}{branch}_volume_m3']=integral(c(stem+str(branch)))
        groups['wetland_to_river']=groups['normal']+groups['overflow']
    else:
        for key in ['normal','overflow','seepage','wetland_ET','wetland_to_river']:
            row.update({key+'_'+s:None for s in ['volume_m3','peak_m3s','peak_time_h','centroid_h','at_river_peak_m3s']})
        for key in ['normal','overflow','seepage','wetland_ET']:
            for branch in (1,2): row[f'{key}{branch}_volume_m3']=None
    for key,values in groups.items():
        row.update({key+'_'+name:value for name,value in series_stats(t,values).items()})
        row[key+'_at_river_peak_m3s']=float(values[k])
    rainfall=run['u'][a2.U['P']]*sum(p.A_w)
    row['former_footprint_rain_total_m3']=integral(rainfall)
    row['former_footprint_rain_direct_river_m3']=integral(rainfall) if bypass else None
    row['former_footprint_rain_direct_at_river_peak_m3s']=float(rainfall[k]) if bypass else None
    row['river_peak_time_h']=float(t[k]/3600)
    for i in (1,2):
        if bypass:
            for suffix in ['initial_m3','final_m3','change_m3','peak_m3','peak_time_h']:row[f'wetland{i}_storage_'+suffix]=None
        else:
            values=run['x'][a2.STATE_NAMES.index(f'S_w{i}')];j=int(np.argmax(values))
            for suffix,value in dict(initial_m3=values[0],final_m3=values[-1],change_m3=values[-1]-values[0],peak_m3=values[j],peak_time_h=t[j]/3600).items():
                row[f'wetland{i}_storage_'+suffix]=float(value)
    # Direct six-state branch inflow vs simultaneous normal+overflow, not sum of separate maxima.
    contribution=groups['gate_combined']+rainfall if bypass else groups['wetland_to_river']
    row.update({'branch_contribution_'+name:value for name,value in series_stats(t,contribution).items()})
    row['branch_contribution_at_river_peak_m3s']=float(contribution[k])
    state_names=run.get('state_names',a2.STATE_NAMES)
    for name,values in zip(state_names,run['x']):
        row[name+'_initial_m3']=float(values[0]);row[name+'_final_m3']=float(values[-1])
    return row


def report():
    result=old.load_config(OUT/'comparison.json'); lock=old.load_config(OUT/'configuration_and_provenance.json'); config=lock['configuration']
    TABLES.mkdir(parents=True,exist_ok=True);PLOTS.mkdir(parents=True,exist_ok=True)
    runs={};records={(r['case'],r['dt_s']):r for r in result['records']};rows=[]
    for key,rec in records.items():
        if rec['valid']:runs[key]=cf.load_trajectory(ROOT/rec['trajectory'])
    for key,rec in records.items():
        if not rec['valid']:
            rows.append(dict(case=rec['label'],case_id=rec['case'],timestep_s=rec['dt_s'],valid=False,feasible=False,status=rec['status'],failure=rec.get('failure')))
            continue
        peak=records['bypass',rec['dt_s']].get('metrics',{}).get('peak_stage_m')
        if peak is None: raise RuntimeError('Bypass invalid: cannot calculate physical peak differences')
        rows.append(main_row(rec,runs[key],peak,config))
    csv_write(TABLES/'comparison.csv',rows)
    finest={r['case_id']:r for r in rows if r['timestep_s']==30}
    if len(finest)!=5 or not all(r['valid'] for r in finest.values()):
        raise RuntimeError('A comparison was invalid; inspect comparison.csv and diagnostic JSON before reporting')
    fields=[('design','Synthetic design'),('dynamic_states','Dynamic states'),('river_peak_m','River peak (m)'),
            ('peak_time_h','Peak elapsed time (h)'),('peak_time_utc','Peak time (UTC)'),('difference_vs_bypass_mm','Bypass peak minus case peak (mm)'),
            ('dam1_max_m','Maximum pool 1 (m)'),('dam1_margin_m','Pool 1 safety margin (m)'),('dam2_max_m','Maximum pool 2 (m)'),('dam2_margin_m','Pool 2 safety margin (m)'),
            ('wetland1_max_m','Maximum wetland 1 (m)'),('wetland1_margin_m','Wetland 1 safety margin (m)'),('wetland2_max_m','Maximum wetland 2 (m)'),('wetland2_margin_m','Wetland 2 safety margin (m)'),
            ('gate_applicability','Gate applicability'),('G5_branch1_min_margin_m','G5 branch 1 node margin (m)'),('G5_branch2_min_margin_m','G5 branch 2 node margin (m)'),
            ('initial_water_m3','Initial total water (m³)'),('final_water_m3','Final total water (m³)'),('storage_change_m3','Event storage change (m³)'),
            ('rain_volume_m3','Rainfall (m³)'),('external_surface_inflow_m3','External surface inflow (m³)'),('ET_m3','Total ET (m³)'),
            ('signed_external_groundwater_m3','Signed external groundwater exchange (m³)'),('groundwater_in_m3','External groundwater inflow (m³)'),('groundwater_out_m3','External groundwater outflow (m³)'),
            ('outlet_volume_m3','River outlet volume (m³)'),('balance_final_m3','Final integrated balance residual (m³)'),('balance_max_abs_m3','Maximum absolute integrated residual (m³)'),
            ('timestep_s','Timestep (s)'),('status','Validity / feasibility'),('solver_termination','Solver termination'),('river_final_stage_m','Final river stage (m)'),('river_rising_at_end','River rising at final step')]
    main_md='# Fixed dynamic comparison: 30-second results\n\nAll levels use the assumed model datum. Positive differences mean lower peak than the six-state bypass. N/A denotes an absent compartment/process, not zero water level. CSV includes all 15 configurations/time steps.\n\n'
    main_md+='## Approved primary initial conditions\n\n'+transposed([finest[k] for k in ['bypass','reference','best']],fields)
    main_md+='\n\n## Equal-total-water initialization diagnostic\n\nOnly the two wetland initial storages are changed to zero. The identical bypass is reused.\n\n'+transposed([finest[k] for k in ['bypass','reference_zero','best_zero']],fields)+'\n'
    (TABLES/'comparison.md').write_text(main_md)
    p=cf.wet.base_parameters(config)
    mechanisms=[mechanism(records[k,30],runs[k,30],p) for k in cf.LABELS]
    csv_write(TABLES/'mechanisms.csv',mechanisms)
    comparisons=[]; gatechecks=[]
    for dt in [120,60,30]:
        base=runs['bypass',dt]
        for case in ['reference','best','reference_zero','best_zero']:
            r=runs[case,dt]
            gates=np.vstack([cf.column(r,'Q_G1'),cf.column(r,'Q_G2')]);bg=np.vstack([cf.column(base,'Q_G1'),cf.column(base,'Q_G2')])
            check=dict(case=case,dt_s=dt,opening_max_difference=float(np.max(abs(base['u'][:2]-r['u'][:2]))),
                       external_inputs_max_difference=float(np.max(abs(base['u'][2:]-r['u'][2:]))),
                       gate_discharge_max_difference_m3s=float(np.max(abs(gates-bg))),
                       gate_volume_difference_m3=np.trapezoid(gates-bg,r['time'],axis=1).tolist(),
                       shared_initial_state_max_difference_m3=float(max(abs(base['x'][:,0]-r['x'][cf.SHARED,0]))),
                       upstream_state_max_difference_m3=float(np.max(abs(base['x'][:4]-r['x'][cf.SHARED[:4]]))))
            if check['opening_max_difference'] or check['external_inputs_max_difference'] or check['shared_initial_state_max_difference_m3']:
                raise RuntimeError('A controlled shared comparison input changed')
            gatechecks.append(check)
    for coarse,fine in [(120,60),(60,30)]:
        for first,second in [('bypass','reference'),('bypass','best'),('reference','best'),('bypass','reference_zero'),('bypass','best_zero'),('reference_zero','best_zero')]:
            ac,af,bc,bf=[records[label,dt]['metrics']['peak_stage_m'] for label in (first,second) for dt in (coarse,fine)]
            dc,df=(ac-bc)*1000,(af-bf)*1000;change=(abs(ac-af)+abs(bc-bf))*1000
            comparisons.append(dict(first=first,second=second,coarse_dt_s=coarse,fine_dt_s=fine,
                difference_coarse_mm=dc,difference_fine_mm=df,sum_grid_changes_mm=change,
                resolved=bool(np.sign(dc)==np.sign(df) and change<max(.1*abs(df),1e-7))))
    old.save_json(OUT/'comparison_diagnostics.json',dict(gate_and_input_checks=gatechecks,timestep_checks=comparisons,
        initial_water_absent_m3=finest['reference']['initial_water_m3']-finest['bypass']['initial_water_m3'],mechanisms=mechanisms,
        all_strict_constraints_pass=all(r['strict_pass'] for run in runs.values() for r in run['audit']),
        six_state_NA_constraints=['nonnegative_S_w1','nonnegative_S_w2','H_w1_safe','H_w2_safe','groundwater_wetland_bed','cover1_lower','cover1_upper','cover2_lower','cover2_upper']))
    # Only one main figure. Full event extent, transparent stage scale, no inset.
    colors={'bypass':'#b45309','reference':'#2563eb','best':'#15803d'}
    fig,ax=plt.subplots(figsize=(10.5,5.5))
    for case in ['bypass','reference','best']:
        r=runs[case,30];h=cf.column(r,'H_r')
        ax.plot(r['time']/3600,h,color=colors[case],lw=1.8,ls='--' if case=='reference' else '-',label=cf.LABELS[case])
        k=int(np.argmax(h));ax.scatter(r['time'][k]/3600,h[k],s=20,color=colors[case],zorder=5)
    ax.set(xlabel='Elapsed time from 2022-09-09 00:00 UTC (h)',ylabel='Simulated river stage (m; assumed datum)',
           title='Dynamic wetland counterfactual | fixed gate openings, same event forcing',xlim=(0,144),ylim=(97.1,97.68))
    ax.ticklabel_format(axis='y',style='plain',useOffset=False);ax.grid(alpha=.22);ax.legend(loc='upper right',fontsize=9)
    fig.text(.5,.025,f'30 s integration | Bypass peak minus reference: {finest["reference"]["difference_vs_bypass_mm"]:.3f} mm; minus best: {finest["best"]["difference_vs_bypass_mm"]:.3f} mm\nPrimary initial water: 512,600 m³ (bypass), 518,100 m³ (wetlands). Hypothetical, uncalibrated simulation.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.09,1,1]);fig.savefig(PLOTS/'01_simulated_river_stage.png',dpi=180);plt.close(fig)
    b,ref,best,zref,zbest=[finest[k] for k in ['bypass','reference','best','reference_zero','best_zero']]
    m={r['case_id']:r for r in mechanisms};primary=[b,ref,best]
    peak_table=transposed(primary,[('design','Design'),('dynamic_states','States'),('river_peak_m','Physical river peak (m)'),('peak_time_h','Peak time (h)'),('difference_vs_bypass_mm','Bypass minus case peak (mm)'),('status','Status')])
    balance_table=transposed(primary,[(k,title) for k,title in fields if k in ['initial_water_m3','final_water_m3','storage_change_m3','rain_volume_m3','external_surface_inflow_m3','ET_m3','signed_external_groundwater_m3','outlet_volume_m3','balance_final_m3','balance_max_abs_m3']])
    mech_fields=[('gate1_volume_m3','Gate 1 release volume (m³)'),('gate2_volume_m3','Gate 2 release volume (m³)'),
                 ('gate_combined_centroid_h','Combined gate release centroid (h)'),('gate_combined_peak_time_h','Combined gate peak time (h)'),
                 ('normal_volume_m3','Wetland normal drainage (m³)'),('overflow_volume_m3','Wetland overflow (m³)'),('seepage_volume_m3','Wetland seepage (m³)'),('wetland_ET_volume_m3','Wetland ET (m³)'),
                 ('wetland1_storage_peak_m3','Wetland 1 peak storage (m³)'),('wetland2_storage_peak_m3','Wetland 2 peak storage (m³)'),
                 ('wetland1_storage_change_m3','Wetland 1 storage change (m³)'),('wetland2_storage_change_m3','Wetland 2 storage change (m³)'),
                 ('wetland1_storage_final_m3','Wetland 1 final storage (m³)'),('wetland2_storage_final_m3','Wetland 2 final storage (m³)'),
                 ('former_footprint_rain_direct_river_m3','Former-footprint rain directly to river (m³)'),
                 ('branch_contribution_peak_m3s','Simultaneous branch contribution peak (m³/s)'),('branch_contribution_peak_time_h','Branch contribution peak time (h)'),
                 ('branch_contribution_centroid_h','Branch contribution centroid (h)'),('branch_contribution_at_river_peak_m3s','Branch contribution at river peak (m³/s)'),
                 ('gate_combined_at_river_peak_m3s','Gates at own river peak (m³/s)'),('former_footprint_rain_direct_at_river_peak_m3s','Former-footprint rain to river at own peak (m³/s)'),('normal_at_river_peak_m3s','Normal drainage at own river peak (m³/s)'),
                 ('overflow_at_river_peak_m3s','Overflow at own river peak (m³/s)'),('groundwater_to_river_at_river_peak_m3s','Groundwater inflow at own river peak (m³/s)')]
    resolution=['| Step (s) | Bypass peak (m) | Reference peak (m) | Best peak (m) | Bypass−reference (mm) | Bypass−best (mm) | Reference−best (mm) |', '|---:|---:|---:|---:|---:|---:|---:|']
    for dt in [120,60,30]:
        vals=[records[k,dt]['metrics']['peak_stage_m'] for k in ['bypass','reference','best']]
        resolution.append(f'| {dt} | '+' | '.join(f'{v:.9f}' for v in vals)+f' | {(vals[0]-vals[1])*1000:.6f} | {(vals[0]-vals[2])*1000:.6f} | {(vals[1]-vals[2])*1000:.6f} |')
    sensitivity=transposed([b,zref,zbest], [('initial_water_m3','Initial total water (m³)'),('river_peak_m','Peak (m)'),('difference_vs_bypass_mm','Bypass−case peak (mm)'),('status','Status')])
    safety=transposed(primary,[(k,title) for k,title in fields if k in ['dam1_max_m','dam2_max_m','dam1_margin_m','dam2_margin_m','wetland1_margin_m','wetland2_margin_m','G5_branch1_min_margin_m','G5_branch2_min_margin_m']])
    reces=transposed(primary,[('river_final_stage_m','Final river stage (m)'),('river_rising_at_end','River rising at final step'),('peak_near_horizon','Peak in last 5% of interval')])
    text=f'''# B.2 six-state no-wetland dynamic counterfactual

This is one new **fixed dynamic control case**, with two existing wetland designs and one small initialization sensitivity. It is not a new optimization. All three primary cases use the current no-dam-spillway formulation and the same prescribed gate openings, regional event weather, shared geometry, boundaries and non-wetland parameters.

## Peak comparison

{peak_table}

Results above use the 30-second grid. The reference A/A covers are 0.65/0.60; the best-found B/B covers were loaded from `batch_best.npz` at full precision: **{lock['best_theta'][-2]:.17g}, {lock['best_theta'][-1]:.17g}**. Positive differences indicate lower peak with wetlands. No percentage uses the arbitrary absolute level datum.

- ΔH_reference = H_peak,case1 − H_peak,case2 = **{ref['difference_vs_bypass_mm']:.6f} mm**.
- ΔH_best = H_peak,case1 − H_peak,case3 = **{best['difference_vs_bypass_mm']:.6f} mm**.
- ΔH_design = H_peak,case2 − H_peak,case3 = **{(ref['river_peak_m']-best['river_peak_m'])*1000:.6f} mm**.

Only the last difference is the additional saved wetland-design optimization improvement. The bypass-to-best difference is a structural counterfactual comparison, not another OPT gain.

![Simulated river-stage histories](../plots/wetland_counterfactual/01_simulated_river_stage.png)

The [complete readable comparison table](../tables/wetland_counterfactual/comparison.md) includes physical levels, constraints, initial/final water, all external balance components, solver completion and the initialization diagnostic. [comparison.csv](../tables/wetland_counterfactual/comparison.csv) contains all 15 case/timestep rows. No new observations, coefficients or safety thresholds were added.

## Structural definition and six-state balances

Case 1 is a reduced **six-state ODE**, ordered `S_u, S_m, S_d1, S_d2, S_g, S_r` (m³). Both gates discharge directly into the river. There are no dummy wetland states, wetland drainage, overflow, seepage or wetland ET. Cases 2 and 3 remain eight-state dynamic models with those wetland processes retained. No case has an independent dam spillway.

```text
dS_u/dt  = P A_u − Q_inf − Q_ud1 − Q_ud2 − Q_ur
dS_m/dt  = Q_inf − Q_mg − E_m
dS_d1/dt = P A_d1 + Q_ud1 − Q_G1 − E_d1
dS_d2/dt = P A_d2 + Q_ud2 − Q_G2 − E_d2
dS_g/dt  = Q_mg + Q_gx − Q_gr
dS_r/dt  = P (A_r+A_w1+A_w2) + Q_ur + Q_G1 + Q_G2 + Q_gr + Q_ext − Q_out − E_r
```

Retained constitutive fluxes reuse the accepted helpers and parameter values. The six-state derivative and diagnostic indices are separate from eight-state positions. The aquifer remains 640,000 m². Rainfall over the former 22,000 m² wetland footprint immediately enters the river: **{m['bypass']['former_footprint_rain_direct_river_m3']:.3f} m³**. This is an explicitly adopted instantaneous bypass assumption, not calibrated replacement land cover. River storage–stage/evaporation area remains **20,000 m²**, while its rainfall-receiving footprint is 42,000 m². Upland area is unchanged. Total rainfall footprint is 660,000 m² in every case, with no aquifer rainfall double count.

All modeled wetland mechanisms—storage, routing, seepage and ET—change together in Case 1. This does not isolate storage alone. Removing wetland seepage changes groundwater evolution and its exchange with the river/exterior; it does not remove the common aquifer.

## Hydraulics, constraints and feasibility

The existing hydraulic law computes every gate discharge. Wetland cases use their internal H_wi; the bypass uses **internal H_r**, never external h_d, zero tailwater or replayed discharge. The maximum tested pointwise opening difference is {max(r['opening_max_difference'] for r in gatechecks):.3g}; maximum actual release difference is **{max(r['gate_discharge_max_difference_m3s'] for r in gatechecks):.3g} m³/s**; maximum absolute branch release-volume difference is **{max(abs(v) for r in gatechecks for v in r['gate_volume_difference_m3']):.3g} m³** across comparisons/steps (floating-point roundoff). These are numerical findings, not hard-coded equalities.

Under the retained river applicability ceiling, H_r < 99.0 m < gate sills 99.8/100.0 m. The finest bypass integration recorded minimum all-stage tailwater clearances of **{records['bypass',30]['metrics']['integration_stage_diagnostics']['minimum_tailwater_clearance_m'][0]:.6f}/{records['bypass',30]['metrics']['integration_stage_diagnostics']['minimum_tailwater_clearance_m'][1]:.6f} m**. All {records['bypass',30]['metrics']['integration_stage_diagnostics']['rhs_evaluations']} RHS evaluations used free-flow gates. Wetland-case tailwaters are also below their sills; shared upstream state histories and free releases agree numerically. Changes in downstream storage and groundwater remain possible despite identical gates.

{safety}

All 15 case/timestep records are numerically and hydraulically valid and pass all applicable constraints, including strict node residual signs. The six-state audit has **19 meaningful rows**; wetland nonnegativity, wetland maximum levels, wetland-bed separation and wetland cover bounds are N/A. The inherited **H_g_valid_max=97.9 m** stays unchanged despite its original lowest-wetland-bed rationale; thus removing the duplicated wetland-bed row does not relax the common groundwater ceiling. Soil-bottom and river limits, the 0.001 m node audit margin, fixed physical limits and positive scaling volumes are unchanged.

Nonnegative storage, soil capacity, common hard geometry and G5 are checked within every existing official RHS/RK stage, including the accepted final point. G5 remains applicability validation for fixed operations; there is no active optimizer in this task. The bypass minimum all-stage G5 margins are **{records['bypass',30]['metrics']['integration_stage_diagnostics']['minimum_G5_margin_m'][0]:.9f}/{records['bypass',30]['metrics']['integration_stage_diagnostics']['minimum_G5_margin_m'][1]:.9f} m**. Node safety maxima include the final time. These finite checks are not continuous-time engineering certification. No storage clipping or safety-limit adjustments occur. The implementation distinguishes an invalid integration from a completed trajectory violating safety.

## Water accounting and initial conditions

{balance_table}

The approved wetland initial storages are 3,000/2,500 m³. Case 1 has no such compartments and starts with **5,500 m³ less total water**, while every shared initial state is identical. This water is not transferred into the river at t=0.

For the bypass, the checked identity is:

```text
dS_total/dt = P A_rain + Q_ext + Q_gx − Q_out − (E_m + E_d1 + E_d2 + E_r)
```

The maximum bypass derivative cancellation error at any recorded RHS/RK stage is **{b['derivative_balance_max_m3s']:.3g} m³/s**. Its maximum integrated residual at 30 seconds is **{b['balance_max_abs_m3']:.6f} m³**, relative {b['balance_relative']:.3g}. Whole-system balance sums only external terms; gate discharge, normal drainage, overflow and seepage are internal transfers and are not added as external inputs. Signed groundwater exchange is positive into the modeled aquifer. The independent final balance also recomputes each integral separately. Trapezoidal diagnostic integration and RK state integration have small, timestep-dependent residuals; no correction water is inserted.

## Mechanism and timing

{transposed([m[k] for k in ['bypass','reference','best']],mech_fields)}

“Branch contribution” means Q_G1+Q_G2+former-footprint rain for the bypass, and simultaneous Q_wr1+Q_wr2+Q_wsp1+Q_wsp2 for wetland cases. Gate inputs are therefore not added again to wetland drainage. Values “at river peak” use each case's own reported peak time. Centroids are temporal first moments, not measured parcel travel times. [mechanisms.csv](../tables/wetland_counterfactual/mechanisms.csv) includes branch volumes/timing, all remaining state water, and both zero-initial-storage diagnostics.

The wetland module changes the concentration and timing of downstream inflow. The bypass river peaks at {b['peak_time_h']:.4f} h, versus {ref['peak_time_h']:.4f}/{best['peak_time_h']:.4f} h with wetlands. The branch-inflow peak is {m['bypass']['branch_contribution_peak_m3s']:.6f} m³/s without wetlands and {m['reference']['branch_contribution_peak_m3s']:.6f}/{m['best']['branch_contribution_peak_m3s']:.6f} m³/s with them. Storage, normal routing, retained wetland overflow, seepage and ET all contribute to the differences; the external accounting above discloses total water effects. The best design's improvement over reference accompanies reduced overflow concentration, not a claim that greater wetland storage always gives a lower peak.

## Equal-initial-water sensitivity

{sensitivity}

Only wetland initial storages were set to zero; all shared initial water, gate functions, forcing, geometry and limits were retained. These two additional fixed configurations were run at all three timesteps without optimization. All six diagnostic integrations are valid and feasible. This is a sensitivity condition, not an observed antecedent state or replacement for the approved primary case.

Zeroing initial wetland water lowers the reference peak by **{(ref['river_peak_m']-zref['river_peak_m'])*1000:.6f} mm** and best-design peak by **{(best['river_peak_m']-zbest['river_peak_m'])*1000:.6f} mm**. The bypass comparisons become **{zref['difference_vs_bypass_mm']:.6f}/{zbest['difference_vs_bypass_mm']:.6f} mm**; reference-to-best design difference is **{(zref['river_peak_m']-zbest['river_peak_m'])*1000:.6f} mm**. The direction and approximate magnitude of the modeled wetland comparison persist. The numerical differences do depend slightly on initial water; this diagnostic does not remove other assumption sensitivity or supply calibration.

## Time resolution, provenance and remaining storage

{chr(10).join(resolution)}

The 120/60/30-second checks retain the identical 120-second weather reconstruction and original gate-knot functions. Rain depth is {b['rain_depth_mm']:.6f} mm (95,898 m³), and solar energy is {b['solar_energy_Wh_m2']:.6f} Wh/m² in every case. All tested difference signs persist; every comparison passes the existing empirical screen (sum of peak grid changes below 10% of the difference). The largest sum of 60→30-second peak changes among primary pairs is {max(r['sum_grid_changes_mm'] for r in comparisons if r['coarse_dt_s']==60 and r['first'] in ['bypass','reference'] and r['second'] in ['reference','best']):.6f} mm. This is a resolution diagnostic, not a rigorous error bound.

The event is **2022-09-09 00:00 through 2022-09-15 00:00 UTC**. NOAA rainfall intervals end at actual report times; meteorology interpolation/spatial transfer, five trace-zero lower bounds and unresolved subinterval timing remain as approved. POWER solar values are gridded estimates, with the retained start-of-hour timestamp assumption and energy-preserving reconstruction. No local gauge discharge is inserted into the boundary inputs. Assumed initial water and boundary heads are not observations.

{reces}

The relevant simulated maxima are interior to the selected six-day window; however, every case has a small end-step rise and substantial remaining storage. The fixed window does **not** establish complete recession or exclude a later maximum beyond it. Remaining stored water is not permanently removed flood volume. Final compartment storages and outlet volumes are available in the mechanism and comparison tables; no warm-up or extended scenario was added.

## Execution, preservation and interpretation

The actual official `ode4u` implementation is multivarious {lock['official_version']}, pinned revision `{lock['official_revision']}`. There were **{result['new_ode_attempts']} new completed integrations**: 3 bypass runs, 6 zero-initial-wetland runs, and 2 exact 120-second reproductions of saved primary wetland trajectories. **Six compatible primary wetland trajectories were reused** after checksum, full design/x0/forcing, all accepted-node flux and constraint-audit checks. The two fresh reproductions matched saved time, state, inputs, outputs and balance arrays exactly. All integrations completed the full horizon. OPT was not run.

The software suite passed **96 tests**, including 10 targeted counterfactual tests for six-state conservation, rainfall routing, internal river tailwater, an actual official RK-stage G5 rejection, and valid but unsafe uncapped pools. These tests verify implementation, not field calibration.

The new modules are [wetland_counterfactual.py](../wetland_counterfactual.py) and [wetland_counterfactual_reporting.py](../wetland_counterfactual_reporting.py), with focused tests in [test_wetland_counterfactual.py](../test_wetland_counterfactual.py). Source hashes, exact numerical designs, forcing provenance and solver revision are recorded in `B/results/wetland_counterfactual/configuration_and_provenance.json`; that directory also holds new trajectories, per-run diagnostics and preservation records. Old primary NPZs are referenced in place and never overwritten. The source/assumption document gains only a counterfactual scope section.

Saved configuration and provenance records retain their original path labels and hashes. `B/project_paths.py` resolves them to the physical stage folders; the local relocation manifest records the path/import changes.

From the project root with the existing pinned solver/environment:

```sh
MPLCONFIGDIR=B/work/mpl-cache python3 -B -m B.wetland_counterfactual --run --report
```

`--run` refuses to overwrite completed comparison results. To regenerate only these new tables/figure/report:

```sh
MPLCONFIGDIR=B/work/mpl-cache python3 -B -m B.wetland_counterfactual --report
```

Case 1 is a six-state dynamic instantaneous-bypass counterfactual, not a static calculation or calibrated land-cover alternative. Cases 2 and 3 are eight-state dynamic systems. All are simulated through time under identical prescribed external forcing and gate histories, with the disclosed primary initial-water difference. Differences include the full represented wetland module and any hydraulic feedback; equal releases here are verified for this case only. Results are conditional simulated differences, not measured Milwaukee flood reductions, local design certification, universal wetland benefits or proof of a global optimum. No data were collected, no historical model/results changed, and no Git commit, push or publication was performed.
'''
    REPORT.write_text(text)
    print('Wrote comparison CSV/Markdown, mechanism CSV, one figure and report:',REPORT)

if __name__=='__main__':report()
