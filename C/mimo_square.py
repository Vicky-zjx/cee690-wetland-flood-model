"""One added squared-rain feature; official FIT and frozen chronological splits."""
from __future__ import annotations

from contextlib import chdir, redirect_stderr, redirect_stdout
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd

from C import mimo_fit as base
from C.mimo_history import rainfall_diagnostics

ROOT = base.ROOT
SOURCE = ROOT/'results/mimo_rain_history'
RESULT = ROOT/'results/mimo_rain_square'
TABLE = ROOT/'tables/mimo_rain_square'
SQUARE = 'P_completed_1h_squared'
FEATURES = base.PREDICTORS + [SQUARE]
PERIOD_COUNTS = {'fit_2019_2022': 33412, 'validation_2023': 8078, 'external_2024': 6951}
PERIOD_YEARS = {'fit_2019_2022': (2019, 2022), 'validation_2023': (2023, 2023), 'external_2024': (2024, 2024)}


def add_square(pairs):
    """Square the observed completed-hour predictor in physical units, not Y."""
    result = pairs.copy()
    rain = result[base.PREDICTORS[0]].to_numpy(dtype=float)
    if not np.isfinite(rain).all() or (rain < 0).any():
        raise ValueError('Squared rainfall requires finite nonnegative observed predictors')
    result[SQUARE] = rain**2
    if not np.isfinite(result[SQUARE]).all():
        raise ValueError('Squared rainfall overflow; no clipping is applied')
    return result


def load_comparison_pairs():
    """Use exactly the saved comparison cohorts; never construct a new mask."""
    periods = {}
    for name, count in PERIOD_COUNTS.items():
        path = SOURCE/f'{name}_pairs.csv'
        frame = pd.read_csv(path, float_precision='round_trip', low_memory=False)
        for column in frame:
            if column.endswith('_utc'):
                frame[column] = pd.to_datetime(frame[column], utc=True)
        lo, hi = PERIOD_YEARS[name]
        if len(frame) != count or frame.timestamp_utc.duplicated().any() or not frame.timestamp_utc.is_monotonic_increasing:
            raise ValueError(f'{name}: original common cohort differs')
        if not frame.timestamp_utc.dt.year.between(lo, hi).all():
            raise ValueError('Chronological period boundary mismatch')
        if not frame.P_predictor_end_utc.le(frame.prediction_cutoff_utc).all():
            raise ValueError('Rain predictor source ends after its cutoff')
        if not np.isfinite(frame[base.PREDICTORS+base.OUTPUTS]).all().all():
            raise ValueError('Saved common cohort contains nonfinite values')
        periods[name] = add_square(frame)
    return periods


def choose_on_validation(metrics):
    selected = metrics.loc[metrics.period.eq('validation_2023') & metrics.variable.eq('P')]
    score = selected.set_index('model').RMSE
    if not np.isfinite(score[['one_hour', 'rain_square']]).all():
        raise ValueError('Both validation rainfall RMSEs are required')
    return 'rain_square' if score['rain_square'] < score['one_hour'] else 'one_hour'


def read_model(directory):
    with np.load(directory/'official_model.npz') as saved:
        model = {key: saved[key] for key in base.MODEL_KEYS}
    return model, json.loads((directory/'scaling.json').read_text())


def fit_square(train, settings):
    """Fit one candidate with the unchanged official nine-input call."""
    out = RESULT/'rain_square'; out.mkdir(parents=True, exist_ok=True)
    X, Y = train[FEATURES].to_numpy().T, train[base.OUTPUTS].to_numpy().T
    scales = {'X': base.fit_affine(X, settings['external_bound']),
              'Y': base.fit_affine(Y, settings['external_bound'])}
    base.save_json(out/'scaling.json', scales)
    Xs, Ys = base.transform(X, scales['X']), base.transform(Y, scales['Y'])
    module, provenance = base.official_module()
    status = dict(status='started', official=provenance, predictors=FEATURES,
        arguments={k:settings[k] for k in ['max_order','pTrain','L1_pnlty','basis_fctn']},
        seed=settings['seed'], scaling=[0,0], X_shape=list(X.shape), Y_shape=list(Y.shape),
        interpretation='Affine in nine expanded inputs; restricted quadratic in original rainfall.')
    base.save_json(out/'execution.json', status)
    ntrain = int(np.floor(len(train)*settings['pTrain']/100))
    order = np.random.RandomState(settings['seed']).permutation(len(train))
    membership = train[['timestamp_utc']].copy()
    membership['internal_split'] = 'internal_diagnostic'
    membership.loc[order[:ntrain], 'internal_split'] = 'coefficient_training'
    old_membership = pd.read_csv(SOURCE/'one_hour/training_membership.csv')
    old_membership['timestamp_utc'] = pd.to_datetime(old_membership.timestamp_utc, utc=True)
    pd.testing.assert_frame_equal(membership, old_membership)
    membership.to_csv(out/'training_membership.csv', index=False)
    B = np.column_stack([np.ones(ntrain), Xs[:,order[:ntrain]].T])
    singular = np.linalg.svd(B, compute_uv=False)
    rank = int((singular > singular[0]*max(B.shape)*np.finfo(float).eps).sum())
    status['coefficient_training_basis'] = dict(terms=B.shape[1], rank=rank,
        condition_number=float(singular[0]/singular[-1]))
    if rank != B.shape[1]:
        raise ValueError('Restricted quadratic design is rank deficient')
    started = time.monotonic(); rng = np.random.get_state()
    try:
        np.random.seed(settings['seed'])
        with (out/'official_stdout_stderr.log').open('w') as log, chdir(out), redirect_stdout(log), redirect_stderr(log):
            returned = module.mimo_srs(Xs, Ys, max_order=1, pTrain=70, scaling=np.array([0,0]),
                L1_pnlty=1.0, basis_fctn='P', var_names={'X':FEATURES, 'Y':base.OUTPUTS})
        counts = [int(v) for v in re.findall(r'^\s*(\d+) data values\s*$',
            (out/'official_stdout_stderr.log').read_text(), re.MULTILINE)]
        if counts != [len(train), len(train)]:
            raise RuntimeError(f'Official filtering changed the common fit pool: {counts}')
        np.testing.assert_array_equal(returned[6], Xs[:,order[ntrain:]])
        np.testing.assert_array_equal(returned[7], Ys[:,order[ntrain:]])
        if not all(np.isfinite(a).all() for a in returned):
            raise FloatingPointError('Official FIT returned nonfinite arrays')
        model = dict(zip(base.MODEL_KEYS, returned[:6]))
        replay, _ = module.compute_model(*returned[:6], returned[6], np.array([0,0]), 'P')
        np.testing.assert_allclose(replay, returned[8], atol=1e-13, rtol=0)
        np.savez(out/'official_model.npz', **model)
        status.update(status='completed', supplied_samples=counts[0], retained_samples=counts[1],
            removed_samples=0, internal_training_samples=ntrain, internal_testing_samples=len(train)-ntrain,
            identical_internal_membership_to_reference=True, returned_shapes=[list(a.shape) for a in returned],
            squared_feature_scaled_coefficients=dict(zip(base.OUTPUTS,model['coeff'][-1])),
            termination='Function returned finite arrays; official mimo_srs exposes no convergence flag.',
            regularization='Initial L1 penalty 1.0; official adaptive penalty behavior unchanged.')
        return model, scales, status
    except BaseException as exc:
        status.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        status['elapsed_seconds'] = time.monotonic()-started
        base.save_json(out/'execution.json', status)
        np.random.set_state(rng)
        import matplotlib.pyplot as plt
        plt.close('all'); plt.ioff()


def evaluate(pairs, models, settings, period):
    from C.uncertainty.weather_generator import solar_geometry
    _, daylight = solar_geometry(pd.DatetimeIndex(pairs.timestamp_utc)-pd.Timedelta(hours=1),43.0547222,-88.0461111)
    frames, rows = [], []
    for name in ['one_hour','rain_square','persistence']:
        features = FEATURES if name=='rain_square' else base.PREDICTORS
        if name == 'persistence':
            values = pairs[base.PREDICTORS[:4]].to_numpy()
        else:
            model, scales = models[name]
            values = base.predict(model, pairs[features].to_numpy().T, scales, settings).T
        f = pairs.copy(); f['period'] = period; f['model'] = name; f['night'] = ~daylight
        for i, variable in enumerate(base.OUTPUTS):
            observed = pairs[variable].to_numpy()
            f[f'{variable}_observed'] = observed
            f[f'{variable}_predicted'] = values[:,i]
            f[f'{variable}_residual'] = observed-values[:,i]
            rows.append(dict(period=period, model=name, variable=variable, units=base.UNITS[i],
                             **base.metrics(observed,values[:,i])))
        frames.append(f)
    return pd.concat(frames,ignore_index=True), pd.DataFrame(rows)


def protect_files():
    path = RESULT/'protected_files_before.json'
    if path.exists():
        return
    records = {}
    skip = {'C/run_fitting.py','C/README.md','C/docs/c2_rain_square_comparison.md'}
    for folder in ['A','B','C']:
        for p in (base.REPO/folder).rglob('*'):
            if not p.is_file() or any(s in {'.git','.venv','venv','__pycache__','mpl-cache','mimo_rain_square'} for s in p.parts):
                continue
            name = p.relative_to(base.REPO).as_posix()
            if name not in skip:
                records[name] = base.sha(p)
    base.save_json(path, records)


def verify_protection():
    previous = json.loads((RESULT/'protected_files_before.json').read_text())
    changed = [name for name,digest in previous.items()
               if not (base.REPO/name).is_file() or base.sha(base.REPO/name)!=digest]
    report = dict(checked_files=len(previous), changed_or_missing=changed)
    base.save_json(RESULT/'preservation_check.json', report)
    if changed:
        raise RuntimeError(f'Protected artifacts changed: {changed}')
    return report


def write_summary(metrics, diagnostics, peak, events, selection, threshold, status):
    rain = metrics.loc[metrics.variable.eq('P')]
    event_summary = events.groupby(['period','model']).agg(complete_events=('total_error_mm','size'),
        total_MAE_mm=('total_error_mm',lambda x:x.abs().mean()),
        peak_MAE_mm_per_h=('peak_error_mm_per_h',lambda x:x.abs().mean()),
        within_window_peak_time_MAE_hours=('peak_time_error_hours',lambda x:x.abs().mean())).reset_index()
    event_summary.to_csv(TABLE/'rain_event_summary.csv', index=False, float_format='%.17g')
    score = rain.set_index(['period','model'])
    changes=[]
    for period in ['validation_2023','external_2024']:
        a,b = score.loc[(period,'one_hour')], score.loc[(period,'rain_square')]
        changes.append(dict(period=period,RMSE_change_percent=(b.RMSE/a.RMSE-1)*100,
                            MAE_change_percent=(b.MAE/a.MAE-1)*100,
                            mean_error_change_mm_per_h=b.mean_prediction_error-a.mean_prediction_error))
    changes = pd.DataFrame(changes); changes.to_csv(TABLE/'rain_metric_changes.csv',index=False,float_format='%.17g')
    coeff_table = pd.read_csv(TABLE/'coefficients.csv',float_precision='round_trip')
    square_terms = coeff_table.loc[coeff_table.term.eq(SQUARE),['output','coefficient_scaled','coefficient_physical','coefficient_units']]
    all_zero = square_terms.coefficient_scaled.eq(0).all()
    term_finding = ('The added squared-rain coefficient is **exactly zero for all four outputs**. '
        'The official fit did not retain this feature. Consequently the actual fitted responses remain affine in the original eight inputs; '
        'changes in scores come from refitted intercepts and other original coefficients, not active quadratic rainfall curvature.'
        if all_zero else 'The actual added-feature coefficients are recorded below; only nonzero values represent active quadratic terms.')
    er = score.loc[('external_2024','one_hour')]; en = score.loc[('external_2024','rain_square')]
    negative = diagnostics.loc[diagnostics.period.eq('external_2024') & diagnostics.variable.eq('P')].set_index('model').negative_predictions
    text = f'''# Restricted quadratic rainfall response: one added squared-rain input

## Controlled change

This exploratory case adds **only `P_completed_1h_squared = P_completed_1h_at_cutoff ** 2`** to the original eight predictors. Squaring is applied to the physical completed native-hour rain predictor before reversible training-only affine scaling. No target value enters the feature. No older rainfall lags or interaction terms are fitted. The official call remains `mimo_srs(..., max_order=1, pTrain=70, scaling=[0,0], L1_pnlty=1.0, basis_fctn='P')`, seed 6902026. The candidate family is affine in the nine expanded inputs and **permits a restricted quadratic response in original rainfall**. There are ten candidate basis terms per output; the coefficient finding below distinguishes the permitted family from the fitted functions.

**Key finding:** {term_finding}

{base.markdown_table(square_terms)}

The original one-hour comparator is reused exactly from the previous controlled fit (2019-2022), including its coefficients and scales; its predictions are replayed and checked against the saved results. The older historical 2019-2023 fit and three-hour trial remain unchanged. The new response has four outputs P, T, Rs and W, with only the one added input differing. No observation, raw file, A/B model or official package is modified.

## Same data, chronology and numerical settings

Exactly the saved common cohorts are used: **33,412** rows in the 2019-2022 fit pool, **8,078** validation rows in 2023 and **6,951** external rows in 2024. The inherited common mask is kept for a controlled comparison even though this candidate does not use the two older rain lags. There is no new exclusion by magnitude, error or feature value. Training scales come only from 2019-2022. The official internal diagnostic split has **23,388** coefficient-training and **10,024** internal diagnostic rows, identical to the comparator; it is not the chronological validation.

The squared term's source interval is the same as its original rain predictor, usually (c-68,c-8] minutes relative to cutoff c. Squaring does not shift time or add future information. The rain target remains (c-8,c+52] rather than the exact next clock hour. Solar intervals retain the documented hour-start assumption; station averaging windows and operational publication latency are not verified. This is retrospective observation-time-aligned evaluation, not an operational or recursive forecast.

The fit uses the unchanged official checkout revision `{status['official']['checkout_revision']}` (source SHA256 `{status['official']['source_sha256']}`). Supplied and retained counts are both **{status['retained_samples']:,}**. On the actual coefficient-training subset, the expanded basis has {status['coefficient_training_basis']['terms']} terms, rank {status['coefficient_training_basis']['rank']}, and condition number {status['coefficient_training_basis']['condition_number']:.6g}. The full two-degree expansion of all calendar predictors is deliberately not used because sin-squared plus cos-squared creates redundant terms. The official function returned finite arrays but exposes no convergence flag; certified convergence is not claimed.

## Validation decision and actual results

Before predicting 2024, the previously declared primary criterion selected **{selection}**: lower 2023 rainfall RMSE, with ties retaining the simpler one-hour comparator. All coefficients and scales were then frozen. No refit with 2023 or tuning using 2024 was performed. Preference under that one criterion is not proof of overall or extreme-rain improvement; the secondary diagnostics below remain part of the interpretation. Both 2023 and 2024 have already been inspected in earlier exploratory work, so this is not a new untouched blind test.

{base.markdown_table(rain)}

Changes are squared-rain candidate minus one-hour comparator; a negative percentage means lower error:

{base.markdown_table(changes)}

The result is a tradeoff, not a demonstrated overall improvement: external rainfall RMSE decreases by {(1-en.RMSE/er.RMSE)*100:.3f}%, but MAE increases by {(en.MAE/er.MAE-1)*100:.2f}%. Mean prediction error moves from {er.mean_prediction_error:.6g} to {en.mean_prediction_error:.6g} mm/h, and negative rain predictions fall from {int(negative['one_hour']):,} to {int(negative['rain_square']):,}. Positive predictions on quantified-zero/trace observations increase, as the strata table shows; fewer negative values do not prove correct rain occurrence. Temperature performance worsens, while the maximum rain remains severely underestimated. Selection on the prespecified rainfall-RMSE criterion is preserved, but the case is not promoted as a generally improved weather model or proof that the quadratic hypothesis works. No significance claim is made.

Mean error is predicted minus observed; plotted residual is observed minus predicted. The [full evaluation](../tables/mimo_rain_square/evaluation.csv) includes all four outputs and persistence on identical rows. Raw negative predictions and out-of-training-range values are not clipped. The [rain strata](../tables/mimo_rain_square/rain_strata.csv) show quantified-zero/trace, positive-rain and heavier-rain performance. The heavier-rain diagnostic threshold remains **{threshold:g} mm/h**, the previous 95th percentile of positive rain in the common fit pool; it is not a hazard threshold.

At the largest observed sample in each comparison year:

{base.markdown_table(peak)}

## Physical diagnostics and event limitations

{base.markdown_table(diagnostics[diagnostics.variable.isin(['P','Rs'])])}

Negative temperatures are valid and are not treated as physical violations. Rainfall and radiation negatives remain explicit. The squared feature is zero when the preceding rain is zero, so it supplies no new information on rain onset. It neither enforces nonnegative precipitation nor provides approaching-storm observations.

{base.markdown_table(event_summary)}

Event diagnostics reuse the unchanged definition: complete consecutive positive native-rain intervals bounded by quantified-zero/trace observations. They measure totals and peak timing **inside the observed wet window**; false rain outside it and delayed predicted peaks beyond it are not captured by event metrics, although their eligible hours remain in overall evaluation. Incomplete boundary/gap windows are counted separately. These are simple wet-spell diagnostics, not certified storm separation or extreme-risk estimates.

## Figures and reproducibility

![Matched comparison](../plots/mimo_rain_square/observed_vs_predicted.png)

![Earliest complete week](../plots/mimo_rain_square/continuous_example.png)

![Rain residuals and observed maximum](../plots/mimo_rain_square/rainfall_diagnostics.png)

The week is the earliest complete 168-hour common 2024 interval, chosen without prediction errors. The maximum-rain view is deliberately selected by observed rainfall, not favorable fit. Missing hours remain gaps.

```sh
cd /Users/zhangjiaxuan/Documents/Codex/Modeling
/opt/miniconda3/bin/python3 -B -m C.run_fitting mimo-square
```

The existing saved comparison cohorts and one-hour model are required; no download is attempted. [Study plan](../results/mimo_rain_square/study_plan.json), [frozen validation decision](../results/mimo_rain_square/frozen_validation_decision.json), [execution status](../results/mimo_rain_square/run_status.json), [coefficients](../tables/mimo_rain_square/coefficients.csv), and [preservation check](../results/mimo_rain_square/preservation_check.json) record reproducible provenance. No commit or push is performed.
'''
    (ROOT/'docs/c2_rain_square_comparison.md').write_text(text)


def run(config):
    settings = dict(base.DEFAULT_SETTINGS, **config.get('mimo',{}))
    if any(settings[k] != v for k,v in {'max_order':1,'pTrain':70,'L1_pnlty':1.0,'basis_fctn':'P','seed':6902026,'external_bound':0.9}.items()):
        raise ValueError('This one-feature comparison preserves the exact previous numerical settings')
    for folder in [RESULT,TABLE,ROOT/'plots/mimo_rain_square']:
        folder.mkdir(parents=True,exist_ok=True)
    protect_files()
    files = [SOURCE/f'{name}_pairs.csv' for name in PERIOD_COUNTS] + [SOURCE/'one_hour'/n for n in
        ['official_model.npz','scaling.json','execution.json','training_membership.csv']]
    source_hashes = {str(p.relative_to(base.REPO)):base.sha(p) for p in files}
    base.save_json(RESULT/'study_plan.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),
        candidate='Only squared latest completed native-hour rain added to original eight inputs.',
        predictors=FEATURES, source_sha256=source_hashes, settings=settings, sample_counts=PERIOD_COUNTS,
        primary_choice='Lower validation_2023 rainfall RMSE; ties keep one_hour. Secondary diagnostics do not disappear.',
        no_final_refit=True, exploratory_reuse='2023 and 2024 previously inspected; not a new blind test.'))
    base.save_json(RESULT/'run_status.json',dict(status='started'))
    periods = load_comparison_pairs()
    train = periods['fit_2019_2022']
    reference, reference_scales = read_model(SOURCE/'one_hour')
    reference_status = json.loads((SOURCE/'one_hour/execution.json').read_text())
    for role,columns in [('X',base.PREDICTORS),('Y',base.OUTPUTS)]:
        expected = base.fit_affine(train[columns].to_numpy().T,settings['external_bound'])
        for key in expected:
            np.testing.assert_array_equal(expected[key],reference_scales[role][key])
    print('Verified identical saved cohorts and reference training-only scales.',flush=True)
    candidate, candidate_scales, status = fit_square(train,settings)
    if status['official'] != reference_status['official']:
        raise RuntimeError('Official source differs from comparator provenance')
    models={'one_hour':(reference,reference_scales),'rain_square':(candidate,candidate_scales)}
    validation, validation_metrics = evaluate(periods['validation_2023'],models,settings,'validation_2023')
    selection = choose_on_validation(validation_metrics)
    frozen_paths = [SOURCE/'one_hour'/n for n in ['official_model.npz','scaling.json']] + [RESULT/'rain_square'/n for n in ['official_model.npz','scaling.json']]
    frozen = {str(p.relative_to(base.REPO)):base.sha(p) for p in frozen_paths}
    base.save_json(RESULT/'frozen_validation_decision.json',dict(frozen_utc=datetime.now(timezone.utc).isoformat(),
        selected=selection, criterion='2023 rainfall RMSE only; ties keep one_hour.',
        validation_metrics=validation_metrics.to_dict('records'), frozen_artifacts=frozen))
    validation.to_csv(RESULT/'validation_predictions.csv',index=False,float_format='%.17g')
    validation_metrics.to_csv(TABLE/'validation_evaluation.csv',index=False,float_format='%.17g')
    print(f'Frozen validation selection: {selection}. Evaluating 2024 without refitting.',flush=True)
    external, external_metrics = evaluate(periods['external_2024'],models,settings,'external_2024')
    external.to_csv(RESULT/'external_predictions.csv',index=False,float_format='%.17g')
    for period in ['validation','external']:
        saved = pd.read_csv(RESULT/f'{period}_predictions.csv',float_precision='round_trip')
        old = pd.read_csv(SOURCE/f'{period}_predictions.csv',float_precision='round_trip')
        columns = [f'{v}_{kind}' for v in base.OUTPUTS for kind in ['observed','predicted']]
        for name in ['one_hour','persistence']:
            a,b = saved.loc[saved.model.eq(name)],old.loc[old.model.eq(name)]
            np.testing.assert_array_equal(a.timestamp_utc.to_numpy(),b.timestamp_utc.to_numpy())
            np.testing.assert_array_equal(a[columns].to_numpy(),b[columns].to_numpy())
        candidate_saved = saved.loc[saved.model.eq('rain_square')]
        model,scales=read_model(RESULT/'rain_square')
        replay=base.predict(model,candidate_saved[FEATURES].to_numpy().T,scales,settings).T
        np.testing.assert_array_equal(replay,candidate_saved[[f'{v}_predicted' for v in base.OUTPUTS]].to_numpy())
    metrics = pd.concat([validation_metrics,external_metrics],ignore_index=True)
    predictions = pd.concat([validation,external],ignore_index=True)
    metrics.to_csv(TABLE/'evaluation.csv',index=False,float_format='%.17g')
    threshold = json.loads((SOURCE/'diagnostic_definition.json').read_text())['heavy_rain_threshold_mm_per_h']
    np.testing.assert_allclose(threshold,train.loc[train.P.gt(0),'P'].quantile(0.95),rtol=0,atol=1e-14)
    strata,events,coverage,physical,peak=rainfall_diagnostics(predictions,threshold)
    for name,frame in [('rain_strata',strata),('rain_events',events),('rain_event_coverage',coverage),
                       ('physical_diagnostics',physical),('largest_observed_rain',peak)]:
        frame.to_csv(TABLE/f'{name}.csv',index=False,float_format='%.17g')
    coefficients=[]; scaling_rows=[]
    for name,(model,scales) in models.items():
        features = FEATURES if name=='rain_square' else base.PREDICTORS
        coeff=model['coeff'].copy()
        coeff[1:] *= np.asarray(scales['Y']['scale'])[None,:]/np.asarray(scales['X']['scale'])[:,None]
        coeff[0]=np.asarray(scales['Y']['midpoint'])+np.asarray(scales['Y']['scale'])*model['coeff'][0]-np.asarray(scales['X']['midpoint'])@coeff[1:]
        check=np.column_stack([np.ones(len(train)),train[features].to_numpy()])@coeff
        np.testing.assert_allclose(check,base.predict(model,train[features].to_numpy().T,scales,settings).T,rtol=1e-12,atol=1e-11)
        feature_units = base.UNITS+['dimensionless']*4+(['(mm/h)^2'] if name=='rain_square' else [])
        for i,term in enumerate(['intercept']+features):
            for j,var in enumerate(base.OUTPUTS):
                coefficients.append(dict(model=name,term=term,output=var,coefficient_scaled=model['coeff'][i,j],
                    coefficient_physical=coeff[i,j],coefficient_units=base.UNITS[j] if i==0 else f'{base.UNITS[j]} per ({feature_units[i-1]})'))
        for role,names,units in [('X',features,feature_units),('Y',base.OUTPUTS,base.UNITS)]:
            for i,var in enumerate(names):
                scaling_rows.append(dict(model=name,role=role,variable=var,units=units[i],
                    training_min=scales[role]['minimum'][i],training_max=scales[role]['maximum'][i],
                    midpoint=scales[role]['midpoint'][i],divisor=scales[role]['scale'][i]))
    pd.DataFrame(coefficients).to_csv(TABLE/'coefficients.csv',index=False,float_format='%.17g')
    pd.DataFrame(scaling_rows).to_csv(TABLE/'affine_scaling.csv',index=False,float_format='%.17g')
    from C.mimo_square_figures import make_figures
    figures=make_figures(predictions,ROOT/'plots/mimo_rain_square',square_coefficients=status['squared_feature_scaled_coefficients'])
    base.save_json(RESULT/'figures.json',figures)
    write_summary(metrics,physical,peak,events,selection,threshold,status)
    if any(base.sha(base.REPO/path)!=digest for path,digest in {**source_hashes,**frozen}.items()):
        raise RuntimeError('Source or frozen model artifact changed during the comparison')
    protection=verify_protection()
    base.save_json(RESULT/'run_status.json',dict(status='completed',official_fit=status,
        reused_reference=str(SOURCE/'one_hour'),reference_and_persistence_replayed_identically=True,
        saved_candidate_predictions_reproduced_exactly=True,source_and_frozen_hashes_unchanged=True,
        selected_on_validation=selection,sample_counts=PERIOD_COUNTS,preservation=protection,
        implementation_sha256={str(p.relative_to(base.REPO)):base.sha(p) for p in
            [Path(__file__),ROOT/'mimo_square_figures.py',ROOT/'run_fitting.py']}))
    print(metrics.loc[metrics.variable.eq('P')].to_string(index=False),flush=True)
    return 0
