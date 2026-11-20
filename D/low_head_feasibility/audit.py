"""Stage D broad-crest compatibility audit, NOT a new hydraulic/ODE model.

Evaluate the user's frozen high-head law against a proposed broad rectangular
crest. Reproduce the original deterministic cases with the unchanged solver.
No candidate blend is connected to either ODE; no weather is regenerated.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import os
import sys
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/cee690-low-head-feasibility-mpl')
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from D.g5_diagnostics import runner
from B.gate_hydraulics import sluice_gate_flow, UnsupportedGateRegime

G=9.80665
CD=0.6
CREST_LENGTH_M=1.0
C_US=np.array([2.6,2.8,3.0,3.1])
C_SI=C_US*math.sqrt(.3048)
C_IDEAL=math.sqrt(G)*(2/3)**1.5
K=CD*math.sqrt(2*G)
SWITCH=1.25
R_MIN=K/C_IDEAL

def write_json(path,value):
    path.write_text(json.dumps(runner.json_clean(value),indent=2,allow_nan=False)+'\n')

def checks():
    tests=[]
    def add(name,ok,evidence):
        tests.append(dict(check=name,pass_regression=bool(ok),evidence=evidence))
    fixed=runner.load_fixed();gates=runner.adapter().wet.gate_parameters(fixed['model_config'])
    add('coefficient_dimensions_and_SI_conversion',math.isclose(3*.3048**3/(.3048*.3048**1.5),C_SI[2]),
        'Q[L^3/T]=C[L^.5/T]*b[L]*H[L]^1.5; C_SI=C_US*sqrt(0.3048).')
    H=.125;b=1.;yc=2*H/3;qcrit=b*math.sqrt(G)*(2*H/3)**1.5
    energy=lambda y,q: y+q*q/(2*G*b*b*y*y)
    add('critical_capacity_derivation',math.isclose(energy(yc,qcrit),H,rel_tol=1e-14),dict(H_m=H,y_critical_m=yc,Q_critical_m3_s=qcrit))
    ys=np.linspace(.05*H,4*H,10000)
    add('specific_energy_minimum',bool(np.all([energy(y,qcrit)>=H-1e-14 for y in ys])), 'Specific energy has its minimum at critical depth.')
    add('length_unit_covariance',math.isclose(100*b*math.sqrt(100*G)*(2*100*H/3)**1.5/qcrit,100**3), 'All lengths and gravity length unit converted by100 -> Q by100^3.')
    add('incompatible_endpoint_detected',K/SWITCH>C_IDEAL,dict(old_effective_coefficient=K/SWITCH,ideal_coefficient=C_IDEAL))
    add('necessary_domain_threshold',math.isclose(K/R_MIN,C_IDEAL,rel_tol=1e-14),dict(necessary_not_sufficient_H_over_a=R_MIN))
    coverage=[]
    for i,gate in enumerate(gates,1):
        a=.125;opening=a/gate.gate_max_opening
        for ratio in (.5,1,1.125,1.25,1.5,R_MIN,2):
            head=ratio*a;up=gate.gate_sill_z+head;tail=gate.gate_sill_z-1
            try:
                f=sluice_gate_flow(up,tail,opening,gate)
                status=f.regime;q=f.discharge
            except UnsupportedGateRegime as e:status=e.regime;q=None
            cap=gate.gate_width*math.sqrt(G)*(2*head/3)**1.5
            coverage.append(dict(branch=i,H_over_a=ratio,H_m=head,a_m=a,original_status=status,Q_original_m3_s=q,
                ideal_broad_crest_capacity_m3_s=cap,original_over_ideal=None if q is None else q/cap,
                H_over_crest_length=head/CREST_LENGTH_M,
                proposed_geometry_nominal_head_range=(.05<=head/CREST_LENGTH_M<=.5)))
            add(f'gate{i}_original_coverage_r{ratio:g}',(q is None)==(ratio<1.25),dict(regime=status,Q_m3_s=q))
        for label,up,down,op in [('closed',gate.gate_sill_z+1,gate.gate_sill_z-1,0),
                                 ('dry',gate.gate_sill_z,gate.gate_sill_z-1,.4),
                                 ('negative_head',gate.gate_sill_z-.1,gate.gate_sill_z-1,.4),
                                 ('equal_head',gate.gate_sill_z+1,gate.gate_sill_z+1,.4),
                                 ('adverse_head',gate.gate_sill_z+1,gate.gate_sill_z+2,.4)]:
            f=sluice_gate_flow(up,down,op,gate);add(f'gate{i}_{label}',f.discharge==0,dict(Q=f.discharge,regime=f.regime))
    pd.DataFrame(coverage).to_csv(HERE/'results/regime_compatibility.csv',index=False)
    sensitivity=[]
    for us,si in zip(C_US,C_SI):
        # A diagnostic direct-Q interpolation discussed by USACE, not a physical
        # implementation. Its failed energy check prevents its use in any ODE.
        r=np.linspace(1,1.25,1001);w=(r-1)/.25
        blended=(1-w)*si+w*K/r
        eps=1e-8
        lower=((1-eps/.25)*si+(eps/.25)*K/(1+eps))
        upper=((eps/.25)*si+(1-eps/.25)*K/(1.25-eps))
        cont=abs(lower-si)<1e-6 and abs(upper-K/1.25)<1e-6
        sensitivity.append(dict(Cw_US_ft_half_per_s=float(us),Cw_SI_m_half_per_s=float(si),
            coefficient_status='Published-range algebraic screen only; not selected/calibrated for an ODE',
            nominal_transition_continuity=cont,endpoint_Q_over_ideal=(K/1.25)/C_IDEAL,
            endpoint_Q_over_free_weir=(K/1.25)/si,
            maximum_blended_effective_C=float(blended.max()),
            transition_sample_count=len(r),samples_exceeding_ideal=int(np.sum(blended>C_IDEAL*(1+1e-12))),
            necessary_gate_onset_H_over_a_for_this_Cw=K/si,
            fixed_head_opening_from_0p8H_to_H_relative_Q_change=si/(K/1.25)-1,
            physical_extension_accepted=False))
        add(f'diagnostic_blend_continuity_Cus{us:g}',cont,'Algebraic continuity only; the energy compatibility test fails.')
        add(f'diagnostic_blend_incompatibility_Cus{us:g}',bool(blended[-1]>C_IDEAL),'Correctly detects incompatible endpoint. This is not a passing physical model.')
    pd.DataFrame(sensitivity).to_csv(HERE/'results/coefficient_sensitivity.csv',index=False)
    result=dict(checks=len(tests),passed=sum(v['pass_regression'] for v in tests),records=tests,
        physical_compatibility_pass=False,new_hydraulic_model_implemented=False,
        note='Regression checks validate the diagnostic calculation and historical behavior; they do not validate a low-head model.')
    write_json(HERE/'results/audit_checks.json',result)
    assert all(v['pass_regression'] for v in tests)

def saved_capacity_audit():
    m=runner.adapter();cfg=runner.load_fixed()['model_config'];gates=m.wet.gate_parameters(cfg)
    sources=[('wetland_reference','wetland','B/results/wetland_only/fixed_reference_120s.npz'),
             ('wetland_best','wetland','B/results/wetland_only/wetland_best_120s.npz'),
             ('no_wetland','no_wetland','B/results/wetland_counterfactual/bypass_120s.npz')]
    rows=[];points=[]
    for label,case,rel in sources:
        saved=m.load_trajectory(ROOT/rel);names=list(m.a2.OUTPUT_NAMES if case=='wetland' else m.dry.OUTPUT_NAMES)
        for i,gate in enumerate(gates):
            H=saved['y'][names.index(f'H_d{i+1}')]-gate.gate_sill_z
            a=gate.gate_max_opening*saved['u'][i]
            q=saved['y'][names.index(f'Q_G{i+1}')]
            cap=gate.gate_width*np.sqrt(G)*(2*H/3)**1.5
            ratio=q/cap;mask=ratio>1+1e-12;k=int(ratio.argmax())
            rows.append(dict(case=label,branch=i+1,source=rel,source_sha256=hashlib.sha256((ROOT/rel).read_bytes()).hexdigest(),
                nodes=len(H),nodes_exceeding_ideal=int(mask.sum()),max_Q_ratio=float(ratio[k]),max_ratio_time_s=float(saved['time'][k]),
                H_m=float(H[k]),a_m=float(a[k]),H_over_a=float(H[k]/a[k]),Q_original_m3_s=float(q[k]),Q_ideal_m3_s=float(cap[k]),
                first_exceedance_s=float(saved['time'][mask][0]) if mask.any() else None,last_exceedance_s=float(saved['time'][mask][-1]) if mask.any() else None))
            if label=='wetland_reference' and i==0:
                points=[dict(time_s=float(t),H_m=float(h),a_m=float(aa),H_over_a=float(h/aa),Q_original_m3_s=float(qq),Q_ideal_m3_s=float(cc),ratio=float(rr)) for t,h,aa,qq,cc,rr in zip(saved['time'],H,a,q,cap,ratio)]
    pd.DataFrame(rows).to_csv(HERE/'results/saved_reference_capacity.csv',index=False)
    pd.DataFrame(points).to_csv(HERE/'results/reference_gate1_capacity_timeseries.csv',index=False)
    return rows,points

def original_regression():
    rows=[]
    for case,design in [('wetland','reference'),('wetland','best'),('no_wetland','best')]:
        r=runner.run_deterministic(case,dt=120,design=design);d=r['diagnostics']
        write_json(HERE/'results'/f'original_{case}_{design}_120s.json',d)
        row=dict(case=case,design=design,valid=r['valid'],physical_peak_m=d.get('physical_river_peak_m'),
                 physical_peak_discharge_m3_s=float(max(r['downstream_discharge'])) if r['valid'] else None,
                 safety_violation=d['safety_violation_any_successful_rhs'],balance_relative=d.get('balance_relative'),
                 derivative_mass_identity_m3_s=d['all_rhs_observations']['derivative_physical_mass_identity_max_abs_m3_s'],
                 **d['regression'])
        rows.append(row);print('Fresh original regression',case,design,row['valid'],row['pass_regression'],flush=True)
    write_json(HERE/'results/original_regression.json',rows)
    assert all(r['valid'] and r['pass_regression'] and not r['safety_violation'] for r in rows)

def scenario_accounting():
    rows=[];base=ROOT/'D/g5_diagnostics/recovered/scenarios'
    for p in sorted(base.glob('*.csv')):
        weather=pd.read_csv(p);meta=json.loads(p.with_suffix('.json').read_text());digest=hashlib.sha256(p.read_bytes()).hexdigest()
        assert digest==meta['forcing_sha256'];cold=bool((weather.temperature_C<=0).any())
        for case in ['wetland','no_wetland']:
            old=json.loads((ROOT/'D/g5_diagnostics/results/original_20_120s'/f'{p.stem}__{case}__best.json').read_text())
            rows.append(dict(scenario_id=p.stem,case=case,seed_sequence_entropy=json.dumps(meta['seed_sequence_entropy']),
                weather_sha256=digest,minimum_temperature_C=float(weather.temperature_C.min()),
                nonpositive_temperature_hours=int((weather.temperature_C<=0).sum()),
                original_status=old['classification'],original_failure_s=old['first_violation']['time_s'],
                proposed_variant_run_status='not_run_cold_outside_liquid_model' if cold else 'not_run_physical_compatibility_blocker',
                new_simulation_started=False,new_peak_river_m=None,new_peak_Q_m3_s=None,new_safety_feasible=None,
                note='No replacement weather, no new seed, no fabricated trajectory or constraint result'))
    pd.DataFrame(rows).to_csv(HERE/'results/original_scenario_accounting.csv',index=False)
    assert len(rows)==40 and sum(r['minimum_temperature_C']>0 for r in rows)==32

def plot(points):
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
    r=np.linspace(.5,2.2,1100);colors=['#6699aa','#407296','#d19b3e','#975848']
    for us,si,col in zip(C_US,C_SI,colors):
        low=r<=1;mid=(r>1)&(r<1.25)
        axes[0].plot(r[low],np.full(low.sum(),si),color=col,lw=1.3)
        axes[0].plot(r[mid],(1-(r[mid]-1)/.25)*si+(r[mid]-1)/.25*K/r[mid],ls='--',color=col,label=f'C US={us:g}: diagnostic blend')
    high=r>=1.25;axes[0].plot(r[high],K/r[high],color='#253852',lw=2,label='Unchanged gate law')
    axes[0].axhline(C_IDEAL,color='black',ls=':',label='Ideal broad-crest capacity')
    axes[0].axvspan(1.25,R_MIN,color='#b94242',alpha=.12)
    axes[0].set(xlabel='H/a',ylabel='Q / (b H$^{3/2}$)  [m$^{1/2}$/s]',title='Continuity does not ensure energy compatibility',xlim=(.5,2.2))
    axes[0].legend(fontsize=7,loc='upper right');axes[0].grid(alpha=.2)
    data=pd.DataFrame(points);part=data[(data.time_s>=57*3600)&(data.time_s<=65*3600)]
    axes[1].plot(part.time_s/3600,part.ratio,color='#253852',lw=1.7)
    bad=part.ratio>1+1e-12;axes[1].scatter(part.time_s[bad]/3600,part.ratio[bad],s=14,color='#b94242',label='24 saved nodes exceed capacity')
    axes[1].axhline(1,color='black',ls=':');axes[1].set(xlabel='Hours since event start',ylabel='Original gate Q / ideal broad-crest capacity',title='Original reference, Gate 1 (saved 120 s grid)')
    axes[1].legend(fontsize=8);axes[1].grid(alpha=.2)
    fig.suptitle('Proposed same-width broad crest: incompatible with the frozen gate branch\nDiagnostic only; no candidate blend was used in an ODE',fontsize=12)
    fig.savefig(HERE/'plots/compatibility.png',dpi=180);plt.close(fig)

def main():
    checks();rows,points=saved_capacity_audit();original_regression();scenario_accounting();plot(points)
    manifest=json.loads((HERE/'preservation_before.json').read_text())
    changed=[p for p,h in manifest.items() if not (ROOT/p).is_file() or hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
    write_json(HERE/'preservation_check.json',dict(files_checked=len(manifest),changed_or_missing=changed,passed=not changed))
    summary=dict(status='blocked_physical_compatibility',C_ideal_m_half_s=C_IDEAL,old_endpoint_C_m_half_s=K/SWITCH,
        endpoint_excess_percent=100*((K/SWITCH)/C_IDEAL-1),necessary_not_sufficient_r_min=R_MIN,
        necessary_not_sufficient_Cd_max=SWITCH*C_IDEAL/math.sqrt(2*G),extra_total_energy_fraction_if_otherwise_fixed=((K/SWITCH)/C_IDEAL)**(2/3)-1,
        warm_scenarios_preserved=16,cold_scenarios_preserved=4,new_warm_simulations_run=0,
        fresh_original_deterministic_runs=3,original_fixed_schedules_preserved=True,
        new_hydraulic_model_implemented=False,hydraulic_pilot_passed=False,ready_for_bounded_D2=False,
        preservation_passed=not changed)
    write_json(HERE/'results/summary.json',summary);print(json.dumps(summary,indent=2))
    assert not changed

if __name__=='__main__':main()
