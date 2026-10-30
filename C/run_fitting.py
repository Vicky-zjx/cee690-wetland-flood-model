"""Existing meteorological fitting and diagnostics only; no ODE propagation.

Run from the repository root: python3 -B -m C.run_fitting all
Use prepare separately to acquire/process source data when required.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import logging
import os
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'results/mpl-cache'))
os.environ.setdefault('MPLBACKEND','Agg')
import numpy as np
import pandas as pd

METHOD_CODES={'copula_var':11,'block_bootstrap':22}
LOG=logging.getLogger('stage_c_fitting')

def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.integer,)): return int(value)
    if isinstance(value, (np.bool_,)): return bool(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, Path): return str(value)
    if isinstance(value, pd.Timestamp): return value.isoformat()
    return value

def save_json(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_safe(value), indent=2, allow_nan=False)+"\n")

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(data): return hashlib.sha256(json.dumps(json_safe(data),sort_keys=True).encode()).hexdigest()


def weather_source_hashes():
    return {f'C/uncertainty/{name}':sha(ROOT/'uncertainty'/name)
            for name in ('marginals.py','weather_generator.py')}


def training_hash(frame):
    """Exactly the row/column hashing convention used by fit_weather."""
    from C.uncertainty.weather_generator import COLUMNS
    values=pd.util.hash_pandas_object(frame[['timestamp_utc']+COLUMNS],index=False).values
    return hashlib.sha256(values.tobytes()).hexdigest()


def verify_fit(fitted,train,config):
    if fitted.get('training_data_hash')!=training_hash(train):
        raise ValueError('Current training observations differ from the fitted model; rerun fit before validation/generation')
    if fitted.get('implementation_sha256')!=weather_source_hashes():
        raise ValueError('Weather implementation changed or its version is absent; rerun fit and validate')
    for key,value in config.get('weather',{}).items():
        if fitted.get('settings',{}).get(key)!=value:
            raise ValueError(f'Weather setting {key} differs from the fitted model; rerun fit')


def training_temporal_report(fitted):
    """Expose existing training diagnostics; do not tune or invent a residual test."""
    var=fitted.get('var',{})
    block=fitted.get('bootstrap_block_diagnostics',{})
    block_ok=block.get('acf_below_prespecified_point2_screen') is True
    return dict(fitted_model_version=fitted.get('model_hash'),training_data_hash=fitted.get('training_data_hash'),
        var=dict(estimator=var.get('estimator'),spectral_radius=var.get('spectral_radius'),
                 stable=bool(var.get('spectral_radius',1)<1),
                 residual_lag1=var.get('residual_lag1'),
                 residual_complete_daylight_count=var.get('residual_complete_daylight_count'),
                 residual_interpretation='Descriptive lag-one diagnostics on genuinely adjacent complete daylight residuals. No unconfigured residual-white-noise acceptance threshold is invented; missing nighttime cloudiness is not filled as observed data.',
                 innovation_covariance_repair=var.get('Q_repair')),
        bootstrap=dict(**block,passed=block_ok,criterion='Maximum absolute training latent ACF at the configured block length must be < 0.2.',
                       criterion_source='Existing prespecified fit_weather block-adequacy screen; unchanged.',
                       interpretation='A failed screen means this block length truncates material remaining persistence. The benchmark has not passed the retained weather-persistence screen.' if not block_ok else 'This persistence screen passed; held-out weather validation is still required.'),
        latent_acf=fitted.get('latent_acf'),bootstrap_block_basis=fitted.get('bootstrap_block_basis'))


def read_config(path):
    config=json.loads(Path(path).read_text())
    if set(config['data']['train_years']) & set(config['data']['validation_years']):
        raise ValueError('Training and validation years overlap')
    if max(config['data']['train_years'])>=min(config['data']['validation_years']):
        raise ValueError('Require chronological whole-year holdout after training')
    for m in config['methods']:
        if m not in METHOD_CODES:raise ValueError('Unknown weather method '+m)
    return config


def data_split(config):
    path=ROOT/config['data']['processed_path']
    data=pd.read_csv(path,float_precision='round_trip')
    data['timestamp_utc']=pd.to_datetime(data.timestamp_utc,utc=True)
    if data.timestamp_utc.duplicated().any():raise ValueError('Duplicated standardized timestamps')
    data=data.sort_values('timestamp_utc').reset_index(drop=True)
    year=data.timestamp_utc.dt.year
    train=data[year.isin(config['data']['train_years'])].copy().reset_index(drop=True)
    valid=data[year.isin(config['data']['validation_years'])].copy().reset_index(drop=True)
    if not len(train) or not len(valid):raise ValueError('Training or chronological holdout has no data')
    return train,valid


def prepare(config):
    from C.uncertainty.data import prepare_data
    result=prepare_data(config,ROOT)
    LOG.info('Prepared %s hourly rows',len(result))


def fit(config):
    from C.uncertainty.weather_generator import fit_weather
    train,valid=data_split(config)
    fitted=fit_weather(train,config)
    fitted['implementation_sha256']=weather_source_hashes()
    fitted.pop('model_hash',None)
    fitted['model_hash']=hashlib.sha256(json.dumps(fitted,sort_keys=True).encode()).hexdigest()
    save_json(ROOT/config['paths']['fitted_model'],fitted)
    save_json(ROOT/'data/processed/training_temporal_diagnostics.json',training_temporal_report(fitted))
    split=dict(train_years=config['data']['train_years'],validation_years=config['data']['validation_years'],
        train_rows=len(train),heldout_rows=len(valid),training_start=train.timestamp_utc.min(),training_end=train.timestamp_utc.max(),
        validation_start=valid.timestamp_utc.min(),validation_end=valid.timestamp_utc.max(),split_rule='Whole UTC calendar years; no shuffled hours. Chicago civil calendar used for conditional weather distributions.',
        standardized_data_sha256=sha(ROOT/config['data']['processed_path']),configuration_sha256=digest(config))
    split_path=ROOT/'data/processed/chronological_split.json'
    if split_path.exists():
        previous=json.loads(split_path.read_text())
        if any(previous.get(k)!=v for k,v in json_safe(split).items() if k!='configuration_sha256'):
            raise ValueError('Existing chronological source/split record differs; review it before an intentional refresh')
        # Preserve original data/split provenance, including its historical config hash.
    else:
        save_json(split_path,split)
    LOG.info('Fitted weather model using %s training hours only',len(train))
    return fitted


def read_fit(config,train=None):
    fitted=json.loads((ROOT/config['paths']['fitted_model']).read_text())
    if train is None:train,_=data_split(config)
    verify_fit(fitted,train,config)
    return fitted


def seed_stream(config, method, scenario_number, phase=1):
    entropy=[int(config['seed']),phase,METHOD_CODES[method],int(scenario_number)]
    return np.random.default_rng(np.random.SeedSequence(entropy)),dict(seed_sequence_entropy=entropy,bit_generator='PCG64')


def validate(config):
    from C.uncertainty.weather_generator import generate_weather
    from C.uncertainty.validation import validate_weather
    train,observed=data_split(config);fitted=read_fit(config,train)
    reports={}
    for method in config['methods']:
        rng,seed=seed_stream(config,method,0,phase=2)
        generated=generate_weather(fitted,pd.DatetimeIndex(observed.timestamp_utc),rng=rng,method=method,training_df=train,config=config)
        out=ROOT/config['paths']['validation_dir']/method
        out.mkdir(parents=True,exist_ok=True)
        reports[method]=validate_weather(observed,generated,config,outdir=out)
        reports[method]['random_stream']=seed
        reports[method]['generation_diagnostics']=generated.attrs.get('generation_diagnostics',{})
        reports[method]['fitted_model_sha256']=sha(ROOT/config['paths']['fitted_model'])
        reports[method]['configuration_sha256']=digest(config)
        reports[method]['training_data_hash']=fitted['training_data_hash']
        reports[method]['weather_implementation_sha256']=fitted['implementation_sha256']
        reports[method]['training_temporal_diagnostics']=training_temporal_report(fitted)
        block_ok=method!='block_bootstrap' or reports[method]['training_temporal_diagnostics']['bootstrap']['passed']
        reports[method]['weather_diagnostics_passed']=bool(reports[method]['passed'] and block_ok)
        reports[method]['diagnostic_scope']='Weather-only held-out checks plus the existing training persistence screen; not proof of completed official FIT/mimo_srs work.'
        save_json(out/'validation.json',reports[method])
        generated.to_csv(out/'generated_heldout_calendar.csv',index=False,float_format='%.17g')
        LOG.info('%s held-out weather validation: passed=%s, failed=%s',method,reports[method]['passed'],reports[method].get('failed_checks'))
    save_json(ROOT/config['paths']['validation_dir']/'summary.json',reports)
    export_weather_validation(config)
    return reports


def export_weather_validation(config):
    """Copy existing validation artifacts to the readable stage output folders."""
    source=ROOT/config['paths']['validation_dir']
    for method in config['methods']:
        plot_dir=ROOT/config['paths']['plot_root']/'weather_validation'/method
        table_dir=ROOT/config['paths']['table_root']/'weather_validation'/method
        plot_dir.mkdir(parents=True,exist_ok=True)
        table_dir.mkdir(parents=True,exist_ok=True)
        for path in (source/method).glob('*.png'):
            shutil.copyfile(path,plot_dir/path.name)
        for name in ('validation_checks.csv','validation.json'):
            shutil.copyfile(source/method/name,table_dir/name)
    shutil.copyfile(source/'summary.json',ROOT/config['paths']['table_root']/'weather_validation/summary.json')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['prepare','fit','validate','all','mimo','mimo-history','mimo-square'])
    parser.add_argument('--config',type=Path,default=ROOT/'configs/fitting.json')
    args=parser.parse_args()
    config=read_config(args.config)
    if args.stage == 'mimo':
        from C.mimo_fit import run
        raise SystemExit(run(config))
    if args.stage == 'mimo-history':
        from C.mimo_history import run
        raise SystemExit(run(config))
    if args.stage == 'mimo-square':
        from C.mimo_square import run
        raise SystemExit(run(config))
    logs=ROOT/'results/logs';logs.mkdir(parents=True,exist_ok=True)
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s',
        handlers=[logging.StreamHandler(),logging.FileHandler(logs/'fitting.log')])
    LOG.info('Starting weather-only %s, configuration=%s',args.stage,digest(config))
    operations=dict(prepare=lambda:prepare(config),fit=lambda:fit(config),validate=lambda:validate(config))
    for name in (['fit','validate'] if args.stage=='all' else [args.stage]):
        operations[name]()
    LOG.info('Finished weather-only %s',args.stage)


if __name__=='__main__':main()
