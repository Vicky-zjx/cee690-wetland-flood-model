"""Immutable paired parameter plans using the inspected official uniform sampler."""
from __future__ import annotations
import importlib
import inspect
import itertools
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
from D.d2_analysis.parameter_uncertainty import adapter
from D.feedback_control import run_policy as frozen
CONFIG=HERE/'config.json'


def physical(multiplier,nominal):
    a,b,c=map(float,multiplier)
    return dict(K_u=nominal['K_u']*a,k_ud=[v*b for v in nominal['k_ud']],
                k_ur=nominal['k_ur']*b,R_r=nominal['R_r']*c)


def point(sample_id,kind,multiplier,nominal,**extra):
    return dict(sample_id=sample_id,study_kind=kind,
                **dict(zip(('m_K_u','m_routing','m_R_r'),map(float,multiplier))),
                physical_values=physical(multiplier,nominal),**extra)


def prepare():
    cfg=json.loads(CONFIG.read_text());target=HERE/'results/study_plan.json'
    if target.exists():
        plan=json.loads(target.read_text())
        if plan['config_sha256']!=frozen.sha(CONFIG):raise RuntimeError('Configuration changed after immutable sampling plan was saved')
        if plan['sampling_source_sha256']!=frozen.sha(__file__):raise RuntimeError('Sampler changed after the plan was saved')
        return plan
    nominal=adapter.nominal_values()
    for key,value in nominal.items():
        if not np.array_equal(value,cfg['nominal_physical_values'][key]):raise ValueError('Declared nominal differs from actual model: '+key)
    official_root=ROOT/'A/work/multivarious'
    if str(official_root) not in sys.path:sys.path.insert(0,str(official_root))
    uniform=importlib.import_module('multivarious.rvs.uniform')
    if not Path(uniform.__file__).resolve().is_relative_to(official_root.resolve()):raise RuntimeError('Wrong multivarious installation')
    settings=cfg['monte_carlo'];children=np.random.SeedSequence(settings['seed']).spawn(3)
    logs=[];seeds=[]
    for name,child in zip(cfg['multiplier_order'],children):
        logs.append(np.asarray(uniform.rnd(np.log(settings['lower_multiplier']),np.log(settings['upper_multiplier']),settings['sample_count'],seed=child),float))
        seeds.append(dict(parameter=name,entropy=child.entropy,spawn_key=list(child.spawn_key),pool_size=child.pool_size,
            state_words=child.generate_state(4).tolist(),argument_type='numpy.random.SeedSequence',
            generator='numpy.random.default_rng inside unchanged official correlated_rvs'))
    multipliers=np.exp(np.column_stack(logs))
    if multipliers.shape!=(64,3) or np.any(multipliers<.5) or np.any(multipliers>2):raise RuntimeError('Unexpected official sample shape/support')
    mc=[point(f'mc_{i:05d}',settings['study_kind'],row,nominal,draw_index=i) for i,row in enumerate(multipliers)]
    grid=[point(f'grid_{i:03d}',cfg['grid']['study_kind'],row,nominal) for i,row in enumerate(itertools.product(cfg['grid']['levels'],repeat=3))]
    pilot=[dict(mc[i],sample_id=f'pilot_mc_{i:05d}',study_kind='numerical_pilot',source_sample_id=mc[i]['sample_id']) for i in range(2)]
    pilot += [point('pilot_low_corner','numerical_pilot',[.5]*3,nominal,source_sample_id=grid[0]['sample_id']),
              point('pilot_high_corner','numerical_pilot',[2.]*3,nominal,source_sample_id=grid[-1]['sample_id'])]
    official=frozen.historical.load_fixed(dt=15,design='best')['official']
    plan=dict(config=cfg,config_sha256=frozen.sha(CONFIG),sampling_source_sha256=frozen.sha(__file__),
        adapter_sha256=frozen.sha(adapter.__file__),official_revision=official['revision'],official_version=official['version'],
        sampler=dict(module=str(Path(uniform.__file__).relative_to(ROOT)),sha256=frozen.sha(uniform.__file__),signature=str(inspect.signature(uniform.rnd)),
            helper_sha256=frozen.sha(official_root/'multivarious/utl/correlated_rvs.py'),calls=3,draws_per_call=64,
            method='Independent official uniform draws in log-multiplier coordinates, exponentiated once'),
        seed=settings['seed'],seed_children=seeds,grid=grid,mc=mc,pilot=pilot,
        nominal_physical_values=nominal,empirical_multiplier_correlation=np.corrcoef(multipliers,rowvar=False).tolist(),
        sampling_interpretation='Independence is a declared exploratory assumption; finite-sample correlations are not forcibly zeroed.',
        stored_before_ODE_execution=True)
    target.parent.mkdir(parents=True,exist_ok=True)
    # Exclusive creation prevents accidental redraws or silent plan replacement.
    with target.open('x') as stream:json.dump(frozen.historical.json_clean(plan),stream,indent=2,allow_nan=False);stream.write('\n')
    rows=[]
    for item in grid+mc+pilot:
        p=item['physical_values'];rows.append({k:v for k,v in item.items() if k!='physical_values'}|dict(K_u=p['K_u'],k_ud_1=p['k_ud'][0],k_ud_2=p['k_ud'][1],k_ur=p['k_ur'],R_r=p['R_r']))
    pd.DataFrame(rows).to_csv(HERE/'results/planned_samples.csv',index=False)
    return plan


if __name__=='__main__':
    p=prepare();print(json.dumps(dict(grid=len(p['grid']),mc=len(p['mc']),pilot=len(p['pilot']),seed=p['seed'],ODE_calls=0)))
