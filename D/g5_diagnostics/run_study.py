"""Reproduce archived G5 failures and evaluate declared fixed schedules only.

No weather generation, parameter optimization or hydraulic closure changes.
Run from the project root with Python -B. All outputs remain beside this file.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time
import traceback
import numpy as np
import pandas as pd
import runner

HERE = Path(__file__).resolve().parent

def save_result(result, directory, name):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / (name+'.json')).write_text(json.dumps(runner.json_clean(result['diagnostics']), indent=2, allow_nan=False)+'\n')
    if result['valid']:
        np.savez_compressed(directory/(name+'.npz'), **{k:result[k] for k in ('time','x','dx','y','u','balance','downstream_stage','downstream_discharge','state_names','output_names')})
    row=runner.failure_record(result)
    d=result['diagnostics']; ob=d['all_rhs_observations']; fc=d.get('forcing_conservation',{})
    row.update(ode_started=d.get('ode_started',False), safety_violation=d.get('safety_violation'),
               rhs_safety_violation=any(v['exceeded'] for v in ob['safety_successful_rhs']) if ob['rhs_successes'] else None,
               balance_relative=d.get('balance_relative'), balance_max_m3=d.get('balance_max_m3'),
               derivative_balance_max_m3s=ob.get('derivative_physical_mass_identity_max_abs_m3_s'),
               rain_depth_m=fc.get('source_precipitation_m'),rain_mapping_error_m=fc.get('precipitation_error_m'),
               solar_mapping_error_Wh_m2=fc.get('solar_error_Wh_m2'),forcing_conservation_pass=fc.get('pass_conservation'),
               forcing_sha256=d.get('forcing_file_sha256'),
               first_exception=(ob['first_rhs_exception'] or {}).get('exception_type'),
               full_output_samples=len(result['time']))
    for v in ob['safety_successful_rhs']:
        row[v['variable']+'_peak_m']=v['peak_m'];row[v['variable']+'_margin_m']=v['margin_m']
    return row

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['deterministic','original','pilot','selected'])
    parser.add_argument('--dt',type=int,default=120)
    parser.add_argument('--scale',type=float,default=1)
    parser.add_argument('--constant',type=float,choices=[.005],default=None)
    parser.add_argument('--tag',default=None)
    parser.add_argument('--ids',nargs='*')
    args=parser.parse_args()
    tag=args.tag or f'{args.mode}_scale_{args.scale:g}_constant_{args.constant}_{args.dt}s'
    output=HERE/'results'/tag;output.mkdir(parents=True,exist_ok=True)
    plans=[]
    if args.mode=='deterministic':
        plans=[('deterministic_september_2022','wetland','reference'),('deterministic_september_2022','wetland','best'),('deterministic_september_2022','no_wetland','best')]
    else:
        ids=[p.stem for p in sorted((HERE/'recovered/scenarios').glob('*.csv'))]
        if args.mode=='pilot': ids=[i for i in ids if int(i.rsplit('_',1)[1])<5]
        if args.mode=='selected': ids=args.ids or []
        plans=[(i,c,'best') for i in ids for c in ('wetland','no_wetland')]
    config=runner.load_fixed(args.dt,args.scale,constant_opening=args.constant)
    (output/'execution_plan.json').write_text(json.dumps(runner.json_clean(dict(mode=args.mode,dt_s=args.dt,scale=args.scale,constant_opening=args.constant,plans=plans,fixed_configuration=config['manifest'])),indent=2)+'\n')
    rows=[]
    for sid,case,design in plans:
        name=f'{sid}__{case}__{design}'
        start=time.monotonic()
        print('START',tag,name,flush=True)
        try:
            if args.mode=='deterministic': result=runner.run_deterministic(case,args.dt,args.scale,design,constant_opening=args.constant)
            else:
                source=HERE/'recovered/scenarios'/(sid+'.csv')
                meta=json.loads(source.with_suffix('.json').read_text())
                if hashlib.sha256(source.read_bytes()).hexdigest()!=meta['forcing_sha256']:raise ValueError('Saved forcing hash mismatch')
                result=runner.run_saved(case,source,args.dt,args.scale,design,scenario_id=sid,seed=meta['seed_sequence_entropy'],constant_opening=args.constant)
            row=save_result(result,output,name)
        except Exception as exc:
            error=dict(scenario_id=sid,case_name=case,design_source_label=design,dt_s=args.dt,gate_scale=args.scale,
                       classification='unexpected_execution_error',valid=False,exception_type=type(exc).__name__,exception_message=str(exc),traceback=traceback.format_exc())
            (output/(name+'.error.json')).write_text(json.dumps(error,indent=2)+'\n')
            row={k:v for k,v in error.items() if k!='traceback'}
        row['elapsed_s']=time.monotonic()-start
        rows.append(row)
        pd.DataFrame(rows).to_csv(output/'scenario_diagnostics.csv',index=False)
        print('DONE',name,row.get('classification'),'peak',row.get('physical_river_peak_m'),'failure_t',row.get('time_s'),'elapsed',round(row['elapsed_s'],2),flush=True)
    counts=pd.DataFrame(rows).groupby(['case_name','classification'],dropna=False).size().to_dict() if rows else {}
    (output/'counts.json').write_text(json.dumps({f'{k[0]} / {k[1]}':int(v) for k,v in counts.items()},indent=2)+'\n')

if __name__=='__main__':main()
