"""Official FIT conditional-response application; no hydrologic simulation.

The source-time preflight is deliberate: hourly row adjacency alone does not
prove that an existing aggregated lag is free of target-hour observations.
"""
from __future__ import annotations

from contextlib import chdir, redirect_stdout, redirect_stderr
import importlib
import inspect
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import tomllib

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
COLUMNS = ['precipitation_mm_per_h', 'temperature_C', 'solar_W_per_m2', 'wind_m_per_s']
OUTPUTS = ['P', 'T', 'Rs', 'W']
UNITS = ['mm/h', 'degC', 'W/m2', 'm/s']
PREDICTORS = ['P_completed_1h_at_cutoff', 'T_latest_at_cutoff', 'Rs_completed_clock_hour', 'W_latest_at_cutoff',
              'sin_annual', 'cos_annual', 'sin_diurnal', 'cos_diurnal']
DEFAULT_SETTINGS = dict(max_order=1, pTrain=70, L1_pnlty=1.0, basis_fctn='P',
                        seed=6902026, external_bound=0.9, smoke_samples=512,
                        timezone='America/Chicago', data_alignment='approved_fit_only_causal_hour_end')
MODEL_KEYS = ['ordr', 'coeff', 'meanX', 'meanY', 'invTX', 'TY']


def sha(path):
    import hashlib
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, value):
    from C.run_fitting import save_json as save
    save(path, value)


def official_module():
    """Import the verified shared checkout, refusing a different installed copy."""
    checkout = REPO / 'A/work/multivarious'
    expected = checkout / 'multivarious/fit/mimo_srs.py'
    if not expected.is_file():
        raise FileNotFoundError(f'Required official checkout unavailable: {expected}')
    sys.path.insert(0, str(checkout)) if str(checkout) not in sys.path else None
    module = importlib.import_module('multivarious.fit.mimo_srs')
    actual = Path(inspect.getsourcefile(module.mimo_srs)).resolve()
    if actual != expected.resolve():
        raise RuntimeError(f'Wrong official FIT import: {actual}; expected {expected}')
    l1_path = Path(inspect.getsourcefile(module.l1)).resolve()
    if not l1_path.is_relative_to(checkout.resolve()):
        raise RuntimeError(f'Wrong L1 dependency import: {l1_path}')
    version = tomllib.loads((checkout / 'pyproject.toml').read_text())['project']['version']
    revision = subprocess.run(['git', '-C', str(checkout), 'rev-parse', 'HEAD'],
                              capture_output=True, text=True, check=True).stdout.strip()
    provenance = dict(source_path=str(actual), source_sha256=sha(actual),
                      package_version=version, checkout_revision=revision,
                      l1_source_path=str(l1_path), l1_sha256=sha(l1_path),
                      example_sha256=sha(checkout / 'examples/mimo_srs_example.py'))
    return module, provenance


def calendar_features(timestamps, timezone='America/Chicago'):
    """Local civil calendar; actual 365/366-day local year, UTC adjacency elsewhere."""
    local = pd.DatetimeIndex(timestamps).tz_convert(timezone)
    hour = local.hour.to_numpy() + local.minute.to_numpy()/60 + local.second.to_numpy()/3600
    year_days = np.where(local.is_leap_year, 366., 365.)
    annual = 2*np.pi*(local.dayofyear.to_numpy()-1+hour/24)/year_days
    diurnal = 2*np.pi*hour/24
    return np.column_stack([np.sin(annual), np.cos(annual), np.sin(diurnal), np.cos(diurnal)])


def make_pairs(weather, timezone='America/Chicago'):
    """Exact adjacent rows only; no interpolation, zero filling or extreme filter."""
    frame = weather.copy()
    frame['timestamp_utc'] = pd.to_datetime(frame['timestamp_utc'], utc=True)
    frame = frame.sort_values('timestamp_utc').reset_index(drop=True)
    if frame['timestamp_utc'].duplicated().any():
        raise ValueError('Duplicate hourly timestamps are not accepted')
    times = frame['timestamp_utc']
    previous = frame[COLUMNS].shift(1)
    exact = times.diff().eq(pd.Timedelta(hours=1))
    current_ok = np.isfinite(frame[COLUMNS]).all(axis=1)
    lag_ok = np.isfinite(previous).all(axis=1)
    reasons = np.select([~exact, ~current_ok, ~lag_ok],
                        ['no_exact_previous_hour', 'missing_or_nonfinite_response',
                         'missing_or_nonfinite_lag'], default='eligible_hourly_pair')
    audit = pd.DataFrame(dict(timestamp_utc=times, predictor_timestamp_utc=times.shift(1),
                              reason=reasons, exact_hour=exact,
                              response_complete=current_ok, lag_complete=lag_ok))
    for column in COLUMNS:
        audit[f'response_valid_{column}'] = np.isfinite(frame[column])
        audit[f'lag_valid_{column}'] = np.isfinite(previous[column])
    valid = exact & current_ok & lag_ok
    pairs = audit.loc[valid, ['timestamp_utc', 'predictor_timestamp_utc']].copy()
    pairs[PREDICTORS[:4]] = previous.loc[valid].to_numpy()
    pairs[PREDICTORS[4:]] = calendar_features(times[valid], timezone)
    pairs[OUTPUTS] = frame.loc[valid, COLUMNS].to_numpy()
    # Retain original standardized source row positions for traceability.
    pairs['source_response_row'] = np.flatnonzero(valid)
    pairs['source_predictor_row'] = pairs['source_response_row']-1
    return pairs.reset_index(drop=True), audit


def attach_pair_support(pairs, weather):
    """Keep cutoff, actual station support and interpreted solar bounds per pair.

    No numeric predictor/response values or sample membership are changed.
    Report-time bounds for T/W and assumed solar bounds are explicitly labelled.
    """
    result=pairs.copy()
    result['prediction_cutoff_utc']=result['predictor_timestamp_utc']
    result['target_reference_time_utc']=result['timestamp_utc']
    result['nominal_target_clock_start_utc']=result['prediction_cutoff_utc']
    result['nominal_target_clock_end_utc']=result['target_reference_time_utc']
    source=weather.copy()
    source['timestamp_utc']=pd.to_datetime(source['timestamp_utc'],utc=True)
    source=source.set_index('timestamp_utc')
    support_types={
        'P':'observed native accumulation interval reconstructed from report end and one-hour duration',
        'T':'reported observation timestamp; instrument averaging-window bounds unavailable',
        'Rs':'interpreted hourly interval; hour-start anchor assumed, not supplied by cached metadata',
        'W':'reported observation timestamp; instrument averaging-window bounds unavailable',
    }
    summary={'scope':'Station report/native accumulation timing is checked; solar interval timing is conditional on the declared hour-start assumption.',
             'all_reported_or_interpreted_predictor_ends_by_cutoff':True,'variables':{}}
    for name,raw in zip(OUTPUTS,['precipitation','temperature','solar','wind']):
        for role,lookup in [('predictor','prediction_cutoff_utc'),('target','target_reference_time_utc')]:
            for bound in ['start','end']:
                values=pd.to_datetime(source[f'{raw}_source_{bound}_utc'],utc=True)
                result[f'{name}_{role}_start_utc' if bound=='start' else f'{name}_{role}_end_utc']=result[lookup].map(values)
            result[f'{name}_{role}_age_minutes']=result[lookup].map(source[f'{raw}_source_age_minutes'])
        result[f'{name}_support_type']=support_types[name]
        end=result[f'{name}_predictor_end_utc']
        if end.isna().any() or not end.le(result.prediction_cutoff_utc).all():
            raise ValueError(f'{name} predictor has missing support or ends after its cutoff')
        if not result[f'{name}_target_end_utc'].le(result.target_reference_time_utc).all():
            raise ValueError(f'{name} target support ends after its selection reference time')
        result[f'{name}_target_end_after_cutoff_minutes']=(result[f'{name}_target_end_utc']-result.prediction_cutoff_utc).dt.total_seconds()/60
        result[f'{name}_target_start_relative_to_cutoff_minutes']=(result[f'{name}_target_start_utc']-result.prediction_cutoff_utc).dt.total_seconds()/60
        summary['variables'][name]={
            'support_type':support_types[name],
            'predictor_age_min_minutes':result[f'{name}_predictor_age_minutes'].min(),
            'predictor_age_max_minutes':result[f'{name}_predictor_age_minutes'].max(),
            'target_end_after_cutoff_min_minutes':result[f'{name}_target_end_after_cutoff_minutes'].min(),
            'target_end_after_cutoff_max_minutes':result[f'{name}_target_end_after_cutoff_minutes'].max(),
            'targets_starting_before_cutoff':int(result[f'{name}_target_start_relative_to_cutoff_minutes'].lt(0).sum()),
            'target_observation_or_interval_ends_after_cutoff':int(result[f'{name}_target_end_after_cutoff_minutes'].gt(0).sum()),
        }
    result['P_predictor_trace_lower_bound']=result.prediction_cutoff_utc.map(source.precipitation_trace_lower_bound)
    result['P_target_trace_lower_bound']=result.target_reference_time_utc.map(source.precipitation_trace_lower_bound)
    selections=[('first_development',result[result.timestamp_utc.dt.year.between(2019,2023)]),
                ('first_2024_cutoff',result[result.prediction_cutoff_utc.dt.year.eq(2024)]),
                ('first_at_or_after_DST_transition',result[result.prediction_cutoff_utc.ge(pd.Timestamp('2024-03-10T08:00Z'))]),
                ('first_nonmodal_rain_predictor_age',result[result.P_predictor_age_minutes.ne(result.P_predictor_age_minutes.mode().iloc[0])])]
    examples=[]
    for label,eligible in selections:
        if eligible.empty: continue
        row=eligible.iloc[0]
        for i,name in enumerate(OUTPUTS):
            examples.append(dict(example=label,variable=name,units=UNITS[i],
                prediction_cutoff_utc=row.prediction_cutoff_utc,
                target_reference_time_utc=row.target_reference_time_utc,
                predictor_name=PREDICTORS[i],predictor_value=row[PREDICTORS[i]],
                predictor_start_utc=row[f'{name}_predictor_start_utc'],
                predictor_end_utc=row[f'{name}_predictor_end_utc'],
                predictor_age_minutes=row[f'{name}_predictor_age_minutes'],
                target_value=row[name],target_start_utc=row[f'{name}_target_start_utc'],
                target_end_utc=row[f'{name}_target_end_utc'],
                target_end_after_cutoff_minutes=row[f'{name}_target_end_after_cutoff_minutes'],
                support_type=support_types[name]))
    summary['example_rule']='First eligible development row; first eligible cutoff within 2024; first eligible cutoff at/after 2024-03-10 08:00 UTC (DST transition); first nonmodal rain-predictor age. No fit-quality selection.'
    summary['exact_future_clock_hour_rainfall_available']=False
    summary['rainfall_limitation']='Native one-hour targets include eight minutes before cutoff in this dataset. Exact clock-hour totals cannot be recovered without adding a redistribution assumption.'
    return result,pd.DataFrame(examples),summary


def split_periods(pairs):
    years = pairs['timestamp_utc'].dt.year
    development = pairs[years.between(2019, 2023)].copy()
    evaluation = pairs[years.eq(2024)].copy()
    if not development.empty and (development['predictor_timestamp_utc'].dt.year > 2023).any():
        raise ValueError('External-period predictor found in development')
    return development, evaluation


def fit_affine(values, bound=0.9):
    values = np.asarray(values, float)
    if values.ndim != 2 or not np.isfinite(values).all() or not 0 < bound < 1:
        raise ValueError('Finite feature-by-sample input and bound inside (0,1) required')
    low, high = values.min(axis=1), values.max(axis=1)
    midpoint = (low+high)/2
    scale = (high-low)/(2*bound)
    scale[scale == 0] = 1.
    return dict(midpoint=midpoint, scale=scale, minimum=low, maximum=high, bound=bound)


def transform(values, parameters, inverse=False):
    midpoint = np.asarray(parameters['midpoint'])[:, None]
    scale = np.asarray(parameters['scale'])[:, None]
    return np.asarray(values)*scale+midpoint if inverse else (np.asarray(values)-midpoint)/scale


def run_official(X, Y, settings, output_dir, purpose):
    """Unmodified official call; capture real counts and verify its random split."""
    module, provenance = official_module()
    output_dir = Path(output_dir).resolve(); output_dir.mkdir(parents=True, exist_ok=True)
    status = dict(status='started', purpose=purpose, official=provenance,
                  arguments={key:settings[key] for key in ['max_order','pTrain','L1_pnlty','basis_fctn']},
                  scaling=[0,0], random_seed=settings['seed'], X_shape=list(X.shape), Y_shape=list(Y.shape))
    save_json(output_dir/'execution.json', status)
    old_rng = np.random.get_state()
    started = time.monotonic()
    try:
        np.random.seed(settings['seed'])
        with (output_dir/'official_stdout_stderr.log').open('w') as log, chdir(output_dir), redirect_stdout(log), redirect_stderr(log):
            returned = module.mimo_srs(X, Y, max_order=settings['max_order'], pTrain=settings['pTrain'],
                scaling=np.array([0,0]), L1_pnlty=settings['L1_pnlty'], basis_fctn=settings['basis_fctn'],
                var_names={'X':PREDICTORS, 'Y':OUTPUTS})
        counts = [int(v) for v in re.findall(r'^\s*(\d+) data values\s*$',
                    (output_dir/'official_stdout_stderr.log').read_text(), re.MULTILINE)]
        if counts != [X.shape[1], X.shape[1]]:
            raise RuntimeError(f'Official count check failed: supplied/retained log counts {counts}')
        train_count = int(np.floor(settings['pTrain']/100*X.shape[1]))
        permutation = np.random.RandomState(settings['seed']).permutation(X.shape[1])
        test_ids = permutation[train_count:]
        np.testing.assert_array_equal(returned[6], X[:,test_ids])
        np.testing.assert_array_equal(returned[7], Y[:,test_ids])
        if not all(np.isfinite(a).all() for a in returned):
            raise FloatingPointError('Official function returned nonfinite arrays')
        model = dict(zip(MODEL_KEYS, returned[:6]))
        reproduced, _ = module.compute_model(*returned[:6], returned[6], np.array([0,0]), settings['basis_fctn'])
        np.testing.assert_allclose(reproduced, returned[8], rtol=0, atol=1e-13)
        status.update(status='completed', returned_shapes=[list(a.shape) for a in returned],
                      supplied_samples=counts[0], retained_samples=counts[1], removed_samples=0,
                      internal_training_samples=train_count, internal_testing_samples=len(test_ids),
                      exact_internal_split_verified=True,
                      internal_filter_threshold=0.8+0.4*np.log(X.shape[1]),
                      termination='Function returned finite arrays; no convergence flag exposed by mimo_srs.',
                      regularization='1.0 is the initial positive penalty; official l1 adapts it internally.')
        np.savez(output_dir/'official_model.npz', **model)
        return model, status, permutation
    except BaseException as exc:
        status.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        status['elapsed_seconds'] = time.monotonic()-started
        save_json(output_dir/'execution.json', status)
        np.random.set_state(old_rng)
        import matplotlib.pyplot as plt
        plt.close('all'); plt.ioff()


def predict(model, X, scaling, settings):
    module, _ = official_module()
    values, _ = module.compute_model(*(model[key] for key in MODEL_KEYS),
                   transform(X, scaling['X']), np.array([0,0]), settings['basis_fctn'])
    return transform(values, scaling['Y'], inverse=True)


def metrics(observed, predicted):
    observed, predicted = np.asarray(observed), np.asarray(predicted)
    if observed.shape != predicted.shape or not np.isfinite([observed,predicted]).all():
        raise ValueError('Metrics require identical finite observations and predictions')
    error = predicted-observed
    denominator = np.sum((observed-observed.mean())**2)
    return dict(N=len(observed), MAE=np.mean(abs(error)), RMSE=np.sqrt(np.mean(error**2)),
                mean_prediction_error=np.mean(error),
                R_squared=1-np.sum(error**2)/denominator if denominator > 0 else np.nan)


def evaluation_products(development, evaluation, settings, result_dir):
    """Fit only the development pool; external responses enter metrics after fit."""
    result_dir = Path(result_dir)
    X = development[PREDICTORS].to_numpy().T
    Y = development[OUTPUTS].to_numpy().T
    scaling = dict(X=fit_affine(X, settings['external_bound']),
                   Y=fit_affine(Y, settings['external_bound']))
    save_json(result_dir/'scaling.json', scaling)
    model, status, permutation = run_official(transform(X,scaling['X']), transform(Y,scaling['Y']),
                                            settings, result_dir/'fit', 'development_pool_fit')
    membership = development.copy()
    membership['internal_split'] = 'internal_diagnostic'
    membership.loc[membership.index[permutation[:status['internal_training_samples']]], 'internal_split'] = 'coefficient_training'
    membership.to_csv(result_dir/'development_pairs.csv', index=False, float_format='%.17g')
    external_X = evaluation[PREDICTORS].to_numpy().T
    predicted = predict(model, external_X, scaling, settings).T
    predictions = evaluation.drop(columns=OUTPUTS+['source_response_row','source_predictor_row'], errors='ignore').copy()
    rows = []
    for i, name in enumerate(OUTPUTS):
        observed = evaluation[name].to_numpy()
        baseline = evaluation[PREDICTORS[i]].to_numpy()
        predictions[f'{name}_observed'] = observed
        predictions[f'{name}_predicted'] = predicted[:,i]
        predictions[f'{name}_persistence'] = baseline
        predictions[f'{name}_residual'] = observed-predicted[:,i]
        for method, values in [('official_mimo_srs',predicted[:,i]), ('persistence',baseline)]:
            rows.append(dict(variable=name, units=UNITS[i], method=method, **metrics(observed,values)))
    from C.uncertainty.weather_generator import solar_geometry
    _, daylight = solar_geometry(pd.DatetimeIndex(predictions['timestamp_utc'])-pd.Timedelta(hours=1),43.0547222,-88.0461111)
    predictions['night'] = ~daylight
    predictions.to_csv(result_dir/'external_predictions.csv', index=False, float_format='%.17g')
    saved = pd.read_csv(result_dir/'external_predictions.csv', float_precision='round_trip')
    with np.load(result_dir/'fit/official_model.npz') as artifact:
        replay_model = {key:artifact[key] for key in MODEL_KEYS}
    replay = predict(replay_model, saved[PREDICTORS].to_numpy().T,
                     json.loads((result_dir/'scaling.json').read_text()), settings).T
    np.testing.assert_array_equal(replay, saved[[f'{name}_predicted' for name in OUTPUTS]].to_numpy())
    return model, scaling, status, predictions, pd.DataFrame(rows)


def run(config):
    """Run the requested task only if the saved lagged data are causal."""
    settings = dict(DEFAULT_SETTINGS, **config.get('mimo', {}))
    if (settings['max_order'], settings['basis_fctn'], settings['L1_pnlty']) != (1,'P',1.0):
        raise ValueError('This focused application uses order 1, power basis, initial L1=1.0; no search.')
    if settings['data_alignment'] != 'approved_fit_only_causal_hour_end':
        raise ValueError('This application requires the user-approved causal FIT-only alignment')
    result_dir = ROOT/'results/mimo_srs'
    table_dir = ROOT/'tables/mimo_srs'
    for directory in (result_dir, table_dir):
        directory.mkdir(parents=True, exist_ok=True)
    path = ROOT/config['data']['processed_path']
    weather = pd.read_csv(path, float_precision='round_trip')
    from C.mimo_source_audit import audit_source_times
    source_audit = audit_source_times(ROOT, weather, config)
    source_dir = result_dir/'source_time_audit'; source_dir.mkdir(exist_ok=True)
    source_audit.to_csv(source_dir/'source_time_audit.csv', index=False)
    save_json(source_dir/'source_time_audit.json', source_audit.attrs['summary'])
    from C.mimo_causal_data import build_causal_weather, ALIGNMENT_DEFINITION
    save_json(result_dir/'alignment_definition.json', ALIGNMENT_DEFINITION)
    causal, causal_manifest = build_causal_weather(ROOT, config)
    causal.to_csv(result_dir/'causal_weather_hourly.csv', index=False, float_format='%.17g')
    save_json(result_dir/'causal_alignment.json', causal_manifest)
    # Read saved precision back: fit the exact artifact offered for reproduction.
    causal = pd.read_csv(result_dir/'causal_weather_hourly.csv', float_precision='round_trip', low_memory=False)
    pairs, eligibility = make_pairs(causal, settings['timezone'])
    pairs, alignment_examples, support_summary = attach_pair_support(pairs, causal)
    alignment_examples.to_csv(table_dir/'alignment_examples.csv',index=False)
    save_json(result_dir/'alignment_support_summary.json',support_summary)
    eligibility.to_csv(result_dir/'hourly_pair_eligibility.csv', index=False)
    source_times = pd.to_datetime(causal['source_latest_utc'], utc=True)
    stamps = pd.to_datetime(causal['timestamp_utc'], utc=True)
    available = pd.Series(source_times.le(stamps).to_numpy(), index=stamps)
    pairs['sources_available_at_forecast_origin'] = pairs['predictor_timestamp_utc'].map(available)
    pairs['response_sources_end_by_response_timestamp'] = pairs['timestamp_utc'].map(available)
    pairs['sources_available_at_forecast_origin'] &= pairs['response_sources_end_by_response_timestamp']
    pairs.to_csv(result_dir/'candidate_hourly_pairs.csv', index=False, float_format='%.17g')
    development, evaluation = split_periods(pairs)
    _, provenance = official_module()
    counts = []
    for label, years, subset in [('development',range(2019,2024),development),
                                 ('external_2024',[2024],evaluation)]:
        rows = eligibility[eligibility['timestamp_utc'].dt.year.isin(years)]
        reasons = rows['reason'].value_counts()
        counts.append(dict(period=label, hourly_rows=len(rows), eligible_adjacent_finite_pairs=len(subset),
                           no_exact_previous_hour=int(reasons.get('no_exact_previous_hour',0)),
                           missing_or_nonfinite_response=int(reasons.get('missing_or_nonfinite_response',0)),
                           missing_or_nonfinite_lag=int(reasons.get('missing_or_nonfinite_lag',0)),
                           future_source_exclusions=int((~subset['sources_available_at_forecast_origin']).sum()),
                           eligible_aligned_pairs=int(subset['sources_available_at_forecast_origin'].sum())))
    pd.DataFrame(counts).to_csv(table_dir/'sample_counts.csv', index=False)
    development = development[development['sources_available_at_forecast_origin']].reset_index(drop=True)
    evaluation = evaluation[evaluation['sources_available_at_forecast_origin']].reset_index(drop=True)
    status = dict(status='source_time_preflight_complete', official=provenance, settings=settings,
                  input_file=str(path), input_sha256=sha(path), sample_counts=counts,
                  real_project_fit_executed=False, external_evaluation_executed=False,
                  original_data_modified=False, causal_alignment=causal_manifest, alignment_support=support_summary,
                  original_hourly_source_time_audit=source_audit.attrs['summary']['counts'])
    if development.empty or evaluation.empty:
        status.update(status='blocked_by_inherited_future_source_support',
            reason='No causal pairs remain in a required period after the approved FIT-only realignment.',
            no_2024_response_fitted_or_scaled=True)
        save_json(result_dir/'run_status.json', status)
        print('BLOCKED: no causal pairs in a required period; see C/results/mimo_srs/run_status.json.', flush=True)
        return 3
    # Executed only with eligible causal observations. No hidden relaxed mode.
    status['status']='running_official_interface_check'
    save_json(result_dir/'run_status.json',status)
    print(f'Causal pairs: development={len(development)}, external={len(evaluation)}. Starting official interface check.', flush=True)
    indices = np.linspace(0,len(development)-1,min(settings['smoke_samples'],len(development)),dtype=int)
    smoke = development.iloc[indices]
    sx,sy = smoke[PREDICTORS].to_numpy().T, smoke[OUTPUTS].to_numpy().T
    run_official(transform(sx,fit_affine(sx)),transform(sy,fit_affine(sy)),settings,
                 result_dir/'smoke','real_development_interface_check')
    status['status']='running_official_fit_and_external_evaluation'
    save_json(result_dir/'run_status.json',status)
    print('Official interface check passed. Starting development fit and external evaluation.', flush=True)
    model, scaling, fit_status, predictions, external = evaluation_products(development,evaluation,settings,result_dir)
    external.to_csv(table_dir/'external_evaluation.csv',index=False,float_format='%.17g')
    from C.mimo_figures import make_figures
    figures = make_figures(predictions, ROOT/'plots/mimo_srs', table_dir)
    pd.DataFrame(model['coeff'],index=['intercept']+PREDICTORS,columns=OUTPUTS).rename_axis('basis_term').to_csv(table_dir/'coefficients_scaled.csv',float_format='%.17g')
    save_json(table_dir/'model_specification.json',dict(settings=settings, predictors=PREDICTORS, outputs=OUTPUTS,
        units=UNITS, alignment=causal_manifest['alignment_definition'], polynomial_orders=model['ordr'], scaling=scaling, official=provenance,
        coefficient_equation='y_scaled = intercept + sum(coefficient_j * x_scaled_j)',
        internal_split='Seeded random 70/30 development split; internal diagnostic only, not an independent chronological test.'))
    negatives = {name:int(predictions[f'{name}_predicted'].lt(0).sum()) for name in ['P','Rs','W']}
    night = predictions.loc[predictions['night'],'Rs_predicted']
    status.update(status='completed',real_project_fit_executed=True,external_evaluation_executed=True,
                  fit=fit_status, figures=figures, negative_predictions=negatives,
                  night_radiation=dict(N=len(night),mean=night.mean(),minimum=night.min(),maximum=night.max(),
                                       positive_count=int(night.gt(0).sum())),
                  saved_predictions_reproduced_exactly=True)
    status['implementation_sha256']={str(p.relative_to(REPO)):sha(p) for p in [ROOT/'run_fitting.py',ROOT/'mimo_fit.py',ROOT/'mimo_causal_data.py',ROOT/'mimo_source_audit.py',ROOT/'mimo_figures.py']}
    status['causal_data_sha256']=sha(result_dir/'causal_weather_hourly.csv')
    write_tables_and_summary(model, scaling, status, predictions, external, settings, table_dir, result_dir)
    save_json(result_dir/'run_status.json',status)
    print(external.to_string(index=False),flush=True)
    return 0


def markdown_table(frame):
    """Small local CSV-to-Markdown rendering without optional table dependencies."""
    def value(item):
        return f'{item:.6g}' if isinstance(item, (float,np.floating)) else str(item)
    return '\n'.join(['| '+' | '.join(map(str,frame.columns))+' |',
                      '| '+' | '.join(['---']*len(frame.columns))+' |']+
                     ['| '+' | '.join(value(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None)])


def write_tables_and_summary(model, scaling, status, predictions, evaluation, settings, table_dir, result_dir):
    scaling_rows=[]
    for role,names,units in [('X',PREDICTORS,UNITS+['dimensionless']*4),('Y',OUTPUTS,UNITS)]:
        p=scaling[role]
        for i,name in enumerate(names):
            scaling_rows.append(dict(role=role,variable=name,units=units[i],development_minimum=p['minimum'][i],
                development_maximum=p['maximum'][i],midpoint=p['midpoint'][i],divisor=p['scale'][i],bound=p['bound']))
    pd.DataFrame(scaling_rows).to_csv(table_dir/'affine_scaling.csv',index=False,float_format='%.17g')
    coefficients=model['coeff']
    xs,xm=scaling['X']['scale'],scaling['X']['midpoint']
    ys,ym=scaling['Y']['scale'],scaling['Y']['midpoint']
    physical=coefficients.copy()
    physical[1:,:]=coefficients[1:,:]*ys[None,:]/xs[:,None]
    physical[0,:]=ym+ys*coefficients[0,:]-xm@physical[1:,:]
    rows=[]
    for term,name in enumerate(['intercept']+PREDICTORS):
        for out,output in enumerate(OUTPUTS):
            input_unit='1' if term==0 else (UNITS+['dimensionless']*4)[term-1]
            rows.append(dict(basis_term=name,output=output,polynomial_order=settings['max_order'],
                basis=settings['basis_fctn'],initial_L1_penalty=settings['L1_pnlty'],
                coefficient_scaled=coefficients[term,out],coefficient_physical=physical[term,out],
                physical_coefficient_units=UNITS[out] if term==0 else f'{UNITS[out]} per ({input_unit})'))
    pd.DataFrame(rows).to_csv(table_dir/'coefficients.csv',index=False,float_format='%.17g')
    unscaled_check=np.column_stack([np.ones(len(predictions)),predictions[PREDICTORS].to_numpy()])@physical
    np.testing.assert_allclose(unscaled_check,predictions[[f'{n}_predicted' for n in OUTPUTS]].to_numpy(),rtol=1e-12,atol=1e-11)
    # Out-of-development-range values remain present and unmodified.
    X=predictions[PREDICTORS].to_numpy().T
    outside=(X<scaling['X']['minimum'][:,None])|(X>scaling['X']['maximum'][:,None])
    status['external_input_values_outside_development_range']={n:int(outside[i].sum()) for i,n in enumerate(PREDICTORS)}
    status['physical_coefficient_equivalence_verified']=True
    diagnostics=[]
    for name,unit in zip(OUTPUTS,UNITS):
        values=predictions[f'{name}_predicted']
        diagnostics.append(dict(variable=name,units=unit,N=len(values),negative_predictions=int(values.lt(0).sum()),
                                minimum=values.min(),maximum=values.max()))
    pd.DataFrame(diagnostics).to_csv(table_dir/'physical_prediction_diagnostics.csv',index=False,float_format='%.17g')
    night=status['night_radiation']
    pd.DataFrame([night]).to_csv(table_dir/'night_radiation.csv',index=False,float_format='%.17g')
    acf=pd.read_csv(table_dir/'residual_acf.csv')
    acf1=acf.loc[acf.lag_hours.eq(1),['variable','pair_count','autocorrelation']]
    counts=pd.DataFrame(status['sample_counts'])
    fit=status['fit']; official=status['official']; window=status['figures']['example_window']
    examples=pd.read_csv(table_dir/'alignment_examples.csv')
    alignment_display=examples.loc[examples.variable.eq('P'),['example','prediction_cutoff_utc','predictor_start_utc','predictor_end_utc','predictor_age_minutes','target_start_utc','target_end_utc']]
    summary=f'''# HW C.2: Official FIT / mimo_srs conditional meteorological response

## Question and scope

How well does a first-order conditional response estimate the observations selected at the next hourly reference time, compared with persistence? This is retrospective, observation-time-aligned evaluation; it is not an exact forecast of four uniform future clock-hour quantities. This is a project-specific application of the professor's tool, not a predictor specification claimed to be instructor-approved. The statistical X/Y labels do not change weather's role as uncontrollable disturbances in A/B.

Let the prediction cutoff be `c` and the next target reference time be `t=c+1 hour`. `Y(t)` contains the four observations selected at t, with their actual supports retained. `X(c) = [P_completed_1h_at_cutoff, T_latest_at_cutoff, Rs_completed_clock_hour, W_latest_at_cutoff, sin(annual_t), cos(annual_t), sin(diurnal_t), cos(diurnal_t)]`. The completed native rainfall accumulation is **not named P(t−1)** because it does not cover the preceding clock hour. Units are mm/h, °C, W/m² and m/s. There is one fitted response surface with nine power-basis terms per output. There is no hydrologic propagation or recursive weather generation.

## Causal alignment approved for this application

The original standardized weather remains unchanged. Auditing its source intervals revealed that every otherwise eligible lagged precipitation value used a report ending 52 minutes inside its target UTC hour: 41,801 development pairs and 6,949 external pairs would violate a strict no-future-source rule. The concrete example at 2019-01-01 02:00 UTC has lag rainfall 0.04 = (8/60)×0.3 and response rainfall 0.26 = (52/60)×0.3 from the same 02:52 report. This is shared target-period source support, not just publication delay. See the [original source-time audit](../results/mimo_srs/source_time_audit/source_time_audit.json).

The user approved a separate FIT-only causal alignment from existing cached observations. The [causal data](../results/mimo_srs/causal_weather_hourly.csv) and [alignment manifest](../results/mimo_srs/causal_alignment.json) preserve each source time, support and age. At hourly UTC **selection/reference time t**, P is the latest accepted native one-hour accumulation ending within the preceding hour; its interval is not split between artificial UTC bins. T and W are the latest valid station reports at or before t, with age no greater than 60 minutes. The 60-minute freshness limit is an explicit eligibility choice; values are neither interpolated nor extrapolated. Rs is the previously standardized solar energy interpreted as covering the completed UTC interval (t−1 h,t], numerically equal to its one-hour mean flux. Each used solar value is preserved exactly; no energy is redistributed.

Thus the components have different physical supports: each native P window usually ends eight minutes before its own selection/reference time, T/W have recorded observation ages, and Rs represents the just-completed UTC hour. The evaluation estimates the next **observation-selection vector**, at a reference time one hour after the predictor cutoff. The reference-time separation does not establish a one-hour physical forecast horizon. It must not be described as four perfectly simultaneous point measurements or uniformly aligned hour averages. Station report times and rainfall accumulation ends are checked against the prediction cutoff c=t−1 h. The solar check remains conditional on the declared hour-start interval interpretation. T/W source start=end records a nominal report timestamp; the cache does not give instrument averaging-window bounds. Operational data-publication latency is unknown, particularly for retrospective NASA POWER grid estimates; this is an offline one-step evaluation, not a demonstrated real-time forecast service. The inherited POWER hour-start support assumption remains explicit and unresolved by cached interval-bound metadata. The [official hourly API description](https://power.larc.nasa.gov/docs/services/api/temporal/hourly/) identifies hourly averages; it does not replace missing bounds in the cached responses.

Missing values stay missing, trace precipitation remains a flagged quantified-zero lower bound, and negative temperatures are valid. No storms or negative predictions are removed. Exact UTC adjacency is required before row pairing. All four outputs and persistence use the same external timestamps. The [pair audit](../results/mimo_srs/hourly_pair_eligibility.csv) gives mutually exclusive primary exclusion reasons and separate variable-validity flags.

## Explicit cutoff and actual source intervals

The [alignment definition](../results/mimo_srs/alignment_definition.json) is saved before rebuilding pairs. Predictors use only records ending by **c**; observations after c may enter targets selected at **t=c+1 hour**. Native rain is not redistributed. The nominal next clock window (c,t] is retained for reference, separately from each actual target support.

In 48,952 retained pairs, the rain predictor covers **(c−68 min,c−8 min]**, with an 8-minute age. One development pair has a 9-minute age and covers **(c−69 min,c−9 min]**; its original timestamp is retained. The rain target covers **(c−8 min,c+52 min]**. Its first 8 minutes precede c. The cached native accumulations therefore cannot provide an exact total over (c,c+1 hour] without an additional redistribution assumption, which this application does not introduce. Results must not be described as an exact one-hour rainfall forecast.

Every pair and prediction CSV now includes `prediction_cutoff_utc`, `target_reference_time_utc`, nominal clock-window bounds, and actual source-specific predictor/target bounds, ages and target end offsets. Four deterministic examples are shown below; [alignment_examples.csv](../tables/mimo_srs/alignment_examples.csv) includes T/W report times and interpreted solar intervals for the same rows. T/W zero-duration entries mean report timestamps, not independently known instrument averaging windows. Solar bounds in that table are **assumed**, not observed bounds.

{markdown_table(alignment_display)}

For retrospective selection, station observations are assumed usable at their reported observation/accumulation end and POWER values at the assumed hour end. No publication-delay or operational availability model is added. The original raw and standardized data remain unchanged. No exact clock-hour rainfall or fully verified solar-interval forecast is claimed.

## Calendar features and separation

Target reference timestamps are converted to `America/Chicago`. Let h be fractional local civil hour, d be local day of year and L be 365 or 366 for that local calendar year. Annual phase = `2π(d−1+h/24)/L`; diurnal phase = `2πh/24`. Leap days are included. Repeated fall-back civil hours share calendar features, and spring-forward features jump; UTC adjacency still requires exactly one elapsed hour.

Target reference years 2019–2023 form the development pool; target reference year 2024 is external. An external row may use its observed lag from late 2023; no 2024 response is supplied to fitting, scaling or the internal split. The first 2019 row has no preceding cached solar hour and remains missing. The last solar hour beginning in 2024 would end in 2025, so it is outside the 2024 target-reference-year evaluation. Each external prediction uses observed lagged inputs, never previous model predictions.

## Official implementation and numerical settings

Actual source: `{official['source_path']}`. Package version `{official['package_version']}`, checkout revision `{official['checkout_revision']}`, source SHA256 `{official['source_sha256']}`. The shared implementation and professor's example were read and remain unchanged.

The exact call is `mimo_srs(X_scaled, Y_scaled, max_order=1, pTrain=70, scaling=np.array([0,0]), L1_pnlty=1.0, basis_fctn='P', var_names=...)`, seed `{settings['seed']}`. Power basis P gives an intercept and eight linear terms. The initial positive L1 penalty 1.0 is a single declared numerical choice, not selected using 2024 or a search. The official L1 routine adapts its penalty internally; it is not described as a certified fixed-penalty solution.

External affine scaling uses development minima/maxima: `z=(value−midpoint)/divisor`, with development range mapped to [−0.9,0.9]. Constant columns would use divisor 1. Scales are reversible and saved. All external values use the same scales without clipping. Official internal scaling is disabled. The official outlier criterion is `0.8+0.4 ln(N)`, rather than simply ±1; bounded development inputs avoid triggering removal. The function printed **{fit['supplied_samples']:,} supplied and {fit['retained_samples']:,} retained**. Actual matrix shapes were X `{fit['X_shape']}` and Y `{fit['Y_shape']}`; coefficient shape `{fit['returned_shapes'][1]}`. Its seeded random split used **{fit['internal_training_samples']:,}** samples for coefficient fitting and **{fit['internal_testing_samples']:,}** for internal diagnostics. This is not an independent time-series holdout. The returned test arrays were checked against the exact seeded permutation.

A small real-development interface call preceded the full call; [smoke execution](../results/mimo_srs/smoke/execution.json) and [full execution](../results/mimo_srs/fit/execution.json) record actual status, arguments, shapes and hashes. Both called the official implementation. It returned finite arrays but exposes no convergence flag; no certified convergence is asserted. Built-in training/testing plots use scaled units and are retained only inside the run folders. Their internal R² is not used here because its denominator differs from the requested observed-data definition.

## Samples and external results

{markdown_table(counts)}

The exclusion columns are sequential, mutually exclusive primary reasons. Source-time checks follow finite, exact-adjacent pairing. Sample selection is unrelated to response magnitude or model errors. `eligible_aligned_pairs` counts samples passing the station report-time checks and the explicitly conditional solar interval check; it does not certify unknown instrument or upstream satellite averaging supports.

{markdown_table(evaluation)}

Metrics are independently calculated in physical units. Mean prediction error is predicted minus observed; residuals in the CSV and figures are observed minus predicted. `R² = 1 − Σ(predicted−observed)² / Σ(observed−mean(observed))²`. Persistence predicts each variable's observed preceding-cutoff value. Neither missing samples nor failures are silently removed differently between methods.

All four external RMSE values are lower than persistence, but improvement is not uniform across metrics: precipitation and wind MAEs are higher. Precipitation R² is only {evaluation.loc[(evaluation.variable=='P') & (evaluation.method=='official_mimo_srs'),'R_squared'].iloc[0]:.6g}. Solar residuals retain substantial diurnal structure; high R² does not establish independent residuals or physical validity.

## Residual and physical limitations

{markdown_table(pd.DataFrame(diagnostics))}

For {night['N']:,} astronomical night target intervals, predicted Rs has mean **{night['mean']:.6g} W/m²**, minimum **{night['minimum']:.6g}**, maximum **{night['maximum']:.6g}**, and **{night['positive_count']:,} positive predictions**. Night classification reuses the existing approximate solar-geometry helper for the completed interval; it is a diagnostic, not a constraint. Predictions remain unchanged.

Lag-one residual correlations, using exact UTC time separation:

{markdown_table(acf1)}

These descriptive residuals are not assumed independent. Negative rainfall/radiation/wind and nonzero night radiation reveal the limitations of an unconstrained affine conditional response. Weak rainfall skill is reported without adding a hurdle model, log transform, tuning campaign or weather generator. The fit does not establish dry-hour probabilities, precipitation tails or a complete joint weather distribution. Station observations and gridded radiation have different spatial support. Source-time alignment does not establish local calibration or operational publication availability.

The retained marginal-distribution, Copula/VAR and bootstrap work is a separate component. Its parameters, results and existing failed diagnostic checks were not altered by this fit. A conditional response cannot replace marginal or distributional evaluation.

## Figures and reproducible artifacts

The example rule is: {window['rule']} The selected window is **{window['start_utc']} through {window['end_utc']}**, {window['n_hours']} hours. Selection never uses prediction error.

![Observed versus predicted](../plots/mimo_srs/observed_vs_predicted.png)

![Continuous example](../plots/mimo_srs/continuous_example.png)

![Residual diagnostics](../plots/mimo_srs/residual_diagnostics.png)

- [Evaluation table](../tables/mimo_srs/external_evaluation.csv), [coefficient table](../tables/mimo_srs/coefficients.csv), [affine scales](../tables/mimo_srs/affine_scaling.csv), [model specification](../tables/mimo_srs/model_specification.json).
- [Timestamped predictions and residuals](../results/mimo_srs/external_predictions.csv), [development rows and actual split](../results/mimo_srs/development_pairs.csv), [saved official arrays](../results/mimo_srs/fit/official_model.npz).
- [Complete execution status](../results/mimo_srs/run_status.json), [official log](../results/mimo_srs/fit/official_stdout_stderr.log), [residual ACF and pair counts](../tables/mimo_srs/residual_acf.csv).

From the project root, the command executed was:

```sh
/opt/miniconda3/bin/python3 -B -m C.run_fitting mimo
```

The existing Python 3.13 environment and existing shared multivarious checkout were used; no package installation or data download was performed. In that environment `python3 -B -m C.run_fitting mimo` is equivalent. `all` continues to run only the retained weather fitting/validation. The new command reads cached sources, creates its own aligned dataset, runs the interface check and official fit, evaluates 2024, and regenerates these figures/tables. Saved predictions are reproduced from saved coefficients, saved scaling and saved predictors; physical-unit coefficient equivalence is also checked.

Relevant retained and focused tests are run with `python3 -B -m unittest discover -s C/tests -v`; the actual latest output is in [tests.log](../results/mimo_srs/tests.log). The focused checks cover exact hourly adjacency, source-time causality, reused solar values, no 2024 development/scaling leakage, official import provenance and replay of saved predictions.

This completes the requested focused statistical fit only. Detailed original assignment/grading instructions were unavailable; no claim is made that every HW C.2 grading requirement or the entire course project is satisfied. A/B, original C weather data, the official checkout and existing weather-generator results remain unchanged. No Git commit or push was made.
'''
    (ROOT/'docs/c2_fit_summary.md').write_text(summary)
