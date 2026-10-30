"""Controlled one- versus three-native-hour rainfall history comparison.

Only the official first-order FIT is used. Existing data and fits are read-only.
"""
from __future__ import annotations

from contextlib import chdir, redirect_stdout, redirect_stderr
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd

from C import mimo_fit as base

ROOT = base.ROOT
RESULT = ROOT / 'results/mimo_rain_history'
TABLE = ROOT / 'tables/mimo_rain_history'
OLDER = ['P_completed_1h_at_cutoff_minus_1h', 'P_completed_1h_at_cutoff_minus_2h']
FEATURES = {'one_hour': base.PREDICTORS, 'three_hour': base.PREDICTORS + OLDER}
PERIODS = {'fit_2019_2022': (2019, 2022), 'validation_2023': (2023, 2023),
           'external_2024': (2024, 2024)}


def build_history_pairs(weather):
    """Match exact UTC selections and contiguous actual native rain supports."""
    weather = weather.copy()
    weather['timestamp_utc'] = pd.to_datetime(weather.timestamp_utc, utc=True)
    pairs, audit = base.make_pairs(weather)
    source = weather.set_index('timestamp_utc')
    cutoff = pairs.predictor_timestamp_utc
    starts, ends = [], []
    complete = pd.Series(True, index=pairs.index)
    valid_support = pd.Series(True, index=pairs.index)
    for offset, name in enumerate([base.PREDICTORS[0]] + OLDER):
        selected_at = cutoff - pd.Timedelta(hours=offset)
        value = selected_at.map(source.precipitation_mm_per_h)
        start = selected_at.map(pd.to_datetime(source.precipitation_source_start_utc, utc=True))
        end = selected_at.map(pd.to_datetime(source.precipitation_source_end_utc, utc=True))
        age = (selected_at-end).dt.total_seconds()/60
        complete &= np.isfinite(value)
        valid_support &= ((end-start).eq(pd.Timedelta(hours=1)) & age.ge(0) & age.lt(60))
        starts.append(start); ends.append(end)
        if offset:
            pairs[name] = value
            pairs[name+'_start_utc'] = start
            pairs[name+'_end_utc'] = end
            pairs[name+'_age_at_cutoff_minutes'] = (cutoff-end).dt.total_seconds()/60
            pairs[name+'_trace_lower_bound'] = selected_at.map(source.precipitation_trace_lower_bound)
    continuous = starts[0].eq(ends[1]) & starts[1].eq(ends[2])
    target_start = pairs.timestamp_utc.map(pd.to_datetime(source.precipitation_source_start_utc, utc=True))
    target_end = pairs.timestamp_utc.map(pd.to_datetime(source.precipitation_source_end_utc, utc=True))
    target_age = (pairs.timestamp_utc-target_end).dt.total_seconds()/60
    target_valid = (target_end-target_start).eq(pd.Timedelta(hours=1)) & target_age.ge(0) & target_age.lt(60)
    reasons = pd.Series(np.select([~complete, ~valid_support, ~continuous, ~target_valid],
        ['missing_older_rain_history', 'invalid_native_rain_support',
         'nonconsecutive_native_rain_intervals', 'invalid_target_rain_support'], default='eligible'), index=pairs.timestamp_utc)
    audit['one_hour_eligible'] = audit.reason.eq('eligible_hourly_pair')
    audit.loc[audit.one_hour_eligible, 'reason'] = audit.loc[audit.one_hour_eligible, 'timestamp_utc'].map(reasons)
    audit['common_eligible'] = audit.reason.eq('eligible')
    selected = complete & valid_support & continuous & target_valid
    pairs, _, _ = base.attach_pair_support(pairs.loc[selected].reset_index(drop=True), weather)
    return pairs, audit


def split_periods(pairs):
    years = pairs.timestamp_utc.dt.year
    return {name: pairs.loc[years.between(lo, hi)].reset_index(drop=True)
            for name, (lo, hi) in PERIODS.items()}


def choose_on_validation(evaluation):
    """Predeclared primary criterion; ties retain the simpler comparator."""
    rows = evaluation.loc[evaluation.period.eq('validation_2023') & evaluation.variable.eq('P')]
    scores = rows.set_index('model').RMSE
    if not np.isfinite(scores[['one_hour', 'three_hour']]).all():
        raise ValueError('Both validation rainfall RMSEs must be finite')
    return 'three_hour' if scores['three_hour'] < scores['one_hour'] else 'one_hour'


def fit_official(train, name, settings):
    """Call the unchanged official function with the correct feature labels."""
    features = FEATURES[name]
    X, Y = train[features].to_numpy().T, train[base.OUTPUTS].to_numpy().T
    scales = {'X': base.fit_affine(X, settings['external_bound']),
              'Y': base.fit_affine(Y, settings['external_bound'])}
    out = RESULT/name; out.mkdir(parents=True, exist_ok=True)
    base.save_json(out/'scaling.json', scales)
    Xs, Ys = base.transform(X, scales['X']), base.transform(Y, scales['Y'])
    module, provenance = base.official_module()
    status = dict(status='started', model=name, official=provenance, predictors=features,
        arguments={key: settings[key] for key in ['max_order','pTrain','L1_pnlty','basis_fctn']},
        scaling=[0, 0], X_shape=list(X.shape), Y_shape=list(Y.shape), seed=settings['seed'])
    base.save_json(out/'execution.json', status)
    rng_state = np.random.get_state(); started = time.monotonic()
    try:
        np.random.seed(settings['seed'])
        with (out/'official_stdout_stderr.log').open('w') as log, chdir(out), redirect_stdout(log), redirect_stderr(log):
            returned = module.mimo_srs(Xs, Ys, max_order=1, pTrain=settings['pTrain'],
                scaling=np.array([0, 0]), L1_pnlty=1.0, basis_fctn='P',
                var_names={'X': features, 'Y': base.OUTPUTS})
        counts = [int(v) for v in re.findall(r'^\s*(\d+) data values\s*$',
                  (out/'official_stdout_stderr.log').read_text(), re.MULTILINE)]
        if counts != [len(train), len(train)]:
            raise RuntimeError(f'Official sample count mismatch: {counts}')
        ntrain = int(np.floor(len(train)*settings['pTrain']/100))
        permutation = np.random.RandomState(settings['seed']).permutation(len(train))
        np.testing.assert_array_equal(returned[6], Xs[:, permutation[ntrain:]])
        np.testing.assert_array_equal(returned[7], Ys[:, permutation[ntrain:]])
        if not all(np.isfinite(a).all() for a in returned):
            raise FloatingPointError('Nonfinite official output')
        model = dict(zip(base.MODEL_KEYS, returned[:6]))
        check, _ = module.compute_model(*returned[:6], returned[6], np.array([0, 0]), 'P')
        np.testing.assert_allclose(check, returned[8], rtol=0, atol=1e-13)
        np.savez(out/'official_model.npz', **model)
        membership = train[['timestamp_utc']].copy()
        membership['internal_split'] = 'internal_diagnostic'
        membership.loc[permutation[:ntrain], 'internal_split'] = 'coefficient_training'
        membership.to_csv(out/'training_membership.csv', index=False)
        status.update(status='completed', supplied_samples=counts[0], retained_samples=counts[1],
            removed_samples=0, internal_training_samples=ntrain, internal_testing_samples=len(train)-ntrain,
            returned_shapes=[list(a.shape) for a in returned], exact_internal_split_verified=True,
            termination='Official function returned finite arrays; no convergence flag is exposed.',
            regularization='Initial positive L1=1.0; official routine adapts the penalty internally.')
        print(f'{name}: official FIT completed, {counts[1]} samples retained.', flush=True)
        return model, scales, status
    except BaseException as exc:
        status.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        status['elapsed_seconds'] = time.monotonic()-started
        base.save_json(out/'execution.json', status)
        np.random.set_state(rng_state)
        import matplotlib.pyplot as plt
        plt.close('all'); plt.ioff()


def evaluate(pairs, models, settings, period):
    frames, rows = [], []
    from C.uncertainty.weather_generator import solar_geometry
    _, daylight = solar_geometry(pd.DatetimeIndex(pairs.timestamp_utc)-pd.Timedelta(hours=1),
                                 43.0547222, -88.0461111)
    for name in ['one_hour', 'three_hour', 'persistence']:
        if name == 'persistence':
            predicted = pairs[base.PREDICTORS[:4]].to_numpy()
        else:
            model, scales, _ = models[name]
            predicted = base.predict(model, pairs[FEATURES[name]].to_numpy().T, scales, settings).T
        frame = pairs.copy(); frame['period'] = period; frame['model'] = name; frame['night'] = ~daylight
        for i, variable in enumerate(base.OUTPUTS):
            observed = pairs[variable].to_numpy()
            frame[f'{variable}_observed'] = observed
            frame[f'{variable}_predicted'] = predicted[:, i]
            frame[f'{variable}_residual'] = observed-predicted[:, i]
            rows.append(dict(period=period, model=name, variable=variable, units=base.UNITS[i],
                             **base.metrics(observed, predicted[:, i])))
        frames.append(frame)
    return pd.concat(frames, ignore_index=True), pd.DataFrame(rows)


def rainfall_diagnostics(predictions, threshold):
    """Descriptive strata and complete observed wet spells; never filter fits."""
    strata, events, event_audit, physical, maxima = [], [], [], [], []
    for (period, name), f in predictions.groupby(['period', 'model'], sort=False):
        f = f.sort_values('timestamp_utc').reset_index(drop=True)
        obs = f.P_observed
        for label, mask in [('zero_or_trace_lower_bound', obs.eq(0)),
                            ('positive_rain', obs.gt(0)),
                            ('upper_5pct_training_positive_rain', obs.ge(threshold))]:
            part = f.loc[mask]
            if len(part):
                strata.append(dict(period=period, model=name, stratum=label,
                    threshold_mm_per_h=threshold if 'upper' in label else 0,
                    positive_predictions=int(part.P_predicted.gt(0).sum()),
                    negative_predictions=int(part.P_predicted.lt(0).sum()),
                    **base.metrics(part.P_observed, part.P_predicted)))
        for variable in base.OUTPUTS:
            p = f[f'{variable}_predicted']
            physical.append(dict(period=period, model=name, variable=variable, N=len(p),
                negative_predictions=int(p.lt(0).sum()), minimum=p.min(), maximum=p.max(),
                night_positive_predictions=int((p.gt(0) & f.night).sum()) if variable=='Rs' else None))
        largest = f.loc[obs.idxmax()]
        maxima.append(dict(period=period, model=name, timestamp_utc=largest.timestamp_utc,
            target_start_utc=largest.P_target_start_utc, target_end_utc=largest.P_target_end_utc,
            observed_mm_per_h=largest.P_observed, predicted_mm_per_h=largest.P_predicted,
            error_mm_per_h=largest.P_predicted-largest.P_observed))
        # Each spell is a run of positive native-hour observations, bounded by
        # known quantified-zero hours. Trace flags do not imply true dry weather.
        wet = obs.gt(0).to_numpy()
        adjacent = (f.timestamp_utc.diff().eq(pd.Timedelta(hours=1)) &
                    f.P_target_start_utc.eq(f.P_target_end_utc.shift(1))).to_numpy()
        starts = np.flatnonzero(wet & ~np.r_[False, wet[:-1] & adjacent[1:]])
        complete_count = 0
        for begin in starts:
            end = begin
            while end+1 < len(f) and wet[end+1] and adjacent[end+1]:
                end += 1
            bounded = (begin>0 and end+1<len(f) and adjacent[begin] and adjacent[end+1]
                       and not wet[begin-1] and not wet[end+1])
            if not bounded:
                continue
            complete_count += 1
            event = f.iloc[begin:end+1]
            oi, pi = event.P_observed.idxmax(), event.P_predicted.idxmax()
            events.append(dict(period=period, model=name, event_start_utc=event.P_target_start_utc.iloc[0],
                event_end_utc=event.P_target_end_utc.iloc[-1], duration_hours=len(event),
                observed_total_mm=event.P_observed.sum(), predicted_total_mm=event.P_predicted.sum(),
                total_error_mm=event.P_predicted.sum()-event.P_observed.sum(),
                observed_peak_mm_per_h=event.P_observed.max(), predicted_peak_mm_per_h=event.P_predicted.max(),
                peak_error_mm_per_h=event.P_predicted.max()-event.P_observed.max(),
                peak_time_error_hours=(f.loc[pi,'timestamp_utc']-f.loc[oi,'timestamp_utc']).total_seconds()/3600))
        event_audit.append(dict(period=period, model=name, observed_positive_runs=len(starts),
            complete_bounded_runs=complete_count, incomplete_boundary_or_gap_runs=len(starts)-complete_count))
    return [pd.DataFrame(x) for x in (strata, events, event_audit, physical, maxima)]


def verify_protected():
    before = json.loads((RESULT/'protected_files_before.json').read_text())
    changed = [name for name, digest in before.items()
               if not (base.REPO/name).is_file() or base.sha(base.REPO/name) != digest]
    report = dict(checked_files=len(before), changed_or_missing=changed,
                  scope='A/B, cached data, historical FIT and pre-existing C scientific artifacts; caches excluded.')
    base.save_json(RESULT/'preservation_check.json', report)
    if changed:
        raise RuntimeError(f'Protected files changed: {changed}')
    return report


def prepare_protection_record():
    """Allow reproduction with an existing aligned cache but no prior manifest."""
    path = RESULT/'protected_files_before.json'
    if path.exists():
        return
    protected = {}
    for folder in ['A', 'B', 'C']:
        for p in (base.REPO/folder).rglob('*'):
            if not p.is_file() or any(part in {'.git','.venv','venv','__pycache__','mpl-cache','mimo_rain_history'} for part in p.parts):
                continue
            if p == ROOT/'docs/c2_rain_history_comparison.md':
                continue
            protected[str(p.relative_to(base.REPO))] = base.sha(p)
    base.save_json(path, protected)


def write_report(evaluation, counts, choice, threshold, physical, maxima, events, figures, models):
    rain = evaluation[evaluation.variable.eq('P')]
    compared = rain.set_index(['period','model'])
    old = compared.loc[('external_2024','one_hour')]
    new = compared.loc[('external_2024','three_hour')]
    negative = physical.loc[physical.period.eq('external_2024') & physical.variable.eq('P')].set_index('model').negative_predictions
    event_summary = events.groupby(['period','model']).agg(
        complete_events=('total_error_mm','size'), total_MAE_mm=('total_error_mm',lambda x: x.abs().mean()),
        peak_MAE_mm_per_h=('peak_error_mm_per_h',lambda x: x.abs().mean()),
        peak_time_MAE_hours=('peak_time_error_hours',lambda x: x.abs().mean())).reset_index()
    event_summary.to_csv(TABLE/'rain_event_summary.csv', index=False, float_format='%.17g')
    text = f'''# Controlled rainfall-history comparison: official first-order FIT

## Design and chronology

The preserved historical fit used 2019-2023 development data. Its code, data, figures and results remain unchanged. Here **one_hour** means that same eight-input specification refitted on the common 2019-2022 pool; **three_hour** adds only two completed native-hour rainfall inputs. Both use the same observations, initial L1 penalty 1.0, power basis, order 1, seed 6902026, and official random 70/30 internal diagnostic split. This is a controlled predictor comparison, not a comparison against historical coefficients trained through 2023. No high-order, penalty or random-seed search was run.

The two extra inputs are selected at cutoff c-1 hour and c-2 hours. Retained actual native supports are (c-68,c-8], (c-128,c-68], (c-188,c-128] minutes. All three must be finite, each exactly 60 minutes, completed by its own selection cutoff, and touch without gaps/overlaps. No missing rain is filled. A real one-minute source gap is excluded explicitly. Shared eligibility is applied to both models and persistence; its reason is unrelated to target magnitude or error. See [pair audit](../results/mimo_rain_history/pair_eligibility.csv).

{base.markdown_table(counts)}

2019-2022 alone supplies the fit pool and reversible affine scales. Both models retain all {models['one_hour'][2]['retained_samples']:,} supplied rows; the official seeded split uses {models['one_hour'][2]['internal_training_samples']:,} coefficient-training and {models['one_hour'][2]['internal_testing_samples']:,} internal diagnostic rows. Internal random splitting is not a chronological validation. 2023 is the separate chronological validation year. Before evaluating 2024, the predeclared rainfall-RMSE rule selected **{choice}** (ties retain one_hour). This is the preferred candidate on that criterion, not a claim of general superiority. Both specifications and coefficients were frozen before evaluation; no refit using 2023 was performed. The same frozen models were evaluated on 2024 for the requested comparison. The previously inspected 2024 historical results are disclosed: this is not a newly untouched blind test.

## Results in physical units

{base.markdown_table(rain)}

**The three-hour candidate does not demonstrate an overall improvement in this controlled case.** It loses the predeclared 2023 rainfall-RMSE comparison, so the one-hour specification remains preferred. In 2024 its RMSE decreases by only {(old.RMSE-new.RMSE)/old.RMSE*100:.4f}%, while MAE increases by {(new.MAE-old.MAE)/old.MAE*100:.2f}%. Negative rain predictions increase from {int(negative['one_hour']):,} to {int(negative['three_hour']):,}. The small external RMSE difference is not used to reverse the frozen validation decision. Neither a statistical-significance claim nor a general conclusion that extra rainfall history can never help is made. This finding applies to these fixed first-order settings, data and years.

Full four-output metrics are in [evaluation.csv](../tables/mimo_rain_history/evaluation.csv). Error is predicted minus observed; plotted residual is observed minus predicted. No predictions are clipped, and no storms are discarded as outliers.

The heavier-rain diagnostic threshold is the 95th percentile of positive rain in the common 2019-2022 fit pool: **{threshold:g} mm/h**. It is a descriptive sample threshold, not an engineering or meteorological warning level. It was fixed before validation/external comparisons. [Rain strata](../tables/mimo_rain_history/rain_strata.csv) distinguish quantified zero/trace lower bounds, positive rain, and values at/above that threshold; positive predictions at trace-zero targets are not proven false alarms.

Complete observed positive-rain spells require adjacent native supports and known quantified-zero/trace boundaries on both sides. This simple diagnostic definition is not a certified storm separation rule. Totals and peak timing are evaluated inside those same observed windows, without including rain from other events. Ties use the earliest peak. Incomplete boundary/gap spells are counted in [event coverage](../tables/mimo_rain_history/rain_event_coverage.csv), not silently completed; all their eligible hourly samples remain in overall metrics. Original negative predictions remain in totals.

These are within-observed-window diagnostics: they do not measure false rain outside the window or delayed predicted peaks occurring beyond it. Those hours remain in overall hourly evaluation and the displayed event detail where eligible.

{base.markdown_table(event_summary)}

At the largest observed rainfall sample in each comparison period:

{base.markdown_table(maxima)}

## Physical and temporal limitations

{base.markdown_table(physical[physical.variable.isin(['P','Rs'])])}

Negative precipitation and radiation and nonzero nighttime radiation remain explicit. Negative temperatures are valid. Extra lagged rain may describe persistence but does not introduce an occurrence model, nonlinear transitions, a nonnegative response, or information about approaching storms. Small changes in aggregate metrics do not establish reliable extreme-rain prediction or flood probabilities.

This remains retrospective observation-time-aligned evaluation. The rain target spans c-8 to c+52 minutes rather than the exact future clock hour. Solar support remains conditional on the existing hour-start interpretation; T/W averaging-window bounds and operational publication latency are unknown. Observed histories are supplied anew at every evaluation time, not recursively simulated. Station and grid quantities retain different spatial supports. A/B, original weather data, weather-model results and the historical FIT are untouched.

## Figures and reproduction

![Matched observed versus predicted](../plots/mimo_rain_history/observed_vs_predicted.png)

![First complete 168-hour example](../plots/mimo_rain_history/continuous_example.png)

![Rainfall residuals and largest observed sample](../plots/mimo_rain_history/rainfall_diagnostics.png)

The week is selected using timestamp coverage only. The largest-rain detail is selected by observed rainfall, not fit quality, and is explicitly an extreme-case diagnostic. Plot rules and timestamps are in [figure metadata](../results/mimo_rain_history/figures.json).

```sh
cd /Users/zhangjiaxuan/Documents/Codex/Modeling
/opt/miniconda3/bin/python3 -B -m C.run_fitting mimo-history
```

The command reads the existing FIT-only aligned cache, calls the unchanged official implementation at revision `{models['one_hour'][2]['official']['checkout_revision']}`, evaluates validation first, saves the decision, then evaluates 2024. The official function returned finite arrays and exposes no convergence flag; certified convergence is not claimed. Saved coefficients and scales reproduce the saved predictions. [Study plan](../results/mimo_rain_history/study_plan.json), [frozen decision](../results/mimo_rain_history/frozen_validation_decision.json), [execution status](../results/mimo_rain_history/run_status.json), and [preservation check](../results/mimo_rain_history/preservation_check.json) provide provenance. No data were downloaded and no Git commit/push was made.
'''
    (ROOT/'docs/c2_rain_history_comparison.md').write_text(text)


def run(config):
    settings = dict(base.DEFAULT_SETTINGS, **config.get('mimo', {}))
    if (settings['max_order'],settings['L1_pnlty'],settings['basis_fctn'],settings['pTrain']) != (1,1.0,'P',70):
        raise ValueError('This comparison fixes order 1, initial L1=1.0, power basis and internal 70/30.')
    for folder in [RESULT, TABLE, ROOT/'plots/mimo_rain_history']:
        folder.mkdir(parents=True, exist_ok=True)
    prepare_protection_record()
    source = ROOT/'results/mimo_srs/causal_weather_hourly.csv'
    if not source.is_file():
        raise FileNotFoundError('Existing FIT-only aligned cache is required; no download/rebuild is attempted.')
    plan = dict(created_utc=datetime.now(timezone.utc).isoformat(), source=str(source), source_sha256=base.sha(source),
        candidates=FEATURES, settings=settings, periods=PERIODS,
        primary_selection='Smaller 2023 rainfall RMSE; ties choose one_hour. Other metrics are reported, not used to tune.',
        heavy_rain_threshold='95th percentile of positive target rainfall in common 2019-2022 fit pool.',
        final_refit=False, shared_training_and_evaluation_timestamps=True,
        external_caveat='Historical 2024 results already inspected; no repeated model selection using 2024.',
        scope='Exactly two first-order fits; no seed, penalty, order or additional feature search.')
    base.save_json(RESULT/'study_plan.json', plan)
    state = {'status':'started', 'source_sha256':plan['source_sha256']}
    base.save_json(RESULT/'run_status.json', state)
    weather = pd.read_csv(source, float_precision='round_trip', low_memory=False)
    pairs, audit = build_history_pairs(weather)
    splits = split_periods(pairs)
    audit.to_csv(RESULT/'pair_eligibility.csv', index=False)
    counts = []
    for period, (lo, hi) in PERIODS.items():
        a = audit.loc[audit.timestamp_utc.dt.year.between(lo, hi)]
        row = dict(period=period, hourly_rows=len(a), one_hour_eligible=int(a.one_hour_eligible.sum()),
                   common_eligible=len(splits[period]))
        for reason in ['missing_older_rain_history','invalid_native_rain_support','nonconsecutive_native_rain_intervals','invalid_target_rain_support']:
            row[reason] = int(a.reason.eq(reason).sum())
        counts.append(row)
        splits[period].to_csv(RESULT/f'{period}_pairs.csv', index=False, float_format='%.17g')
    counts = pd.DataFrame(counts); counts.to_csv(TABLE/'sample_counts.csv', index=False)
    print(counts.to_string(index=False), flush=True)
    train = splits['fit_2019_2022']
    threshold = float(train.loc[train.P.gt(0), 'P'].quantile(0.95))
    base.save_json(RESULT/'diagnostic_definition.json', {'heavy_rain_threshold_mm_per_h': threshold,
        'basis':'95th percentile of positive common 2019-2022 target rain, before held-out evaluation'})
    models = {name: fit_official(train, name, settings) for name in FEATURES}
    # Freeze the decision using 2023 only, before even predicting the external rows.
    valid, valid_metrics = evaluate(splits['validation_2023'], models, settings, 'validation_2023')
    choice = choose_on_validation(valid_metrics)
    frozen = dict(frozen_utc=datetime.now(timezone.utc).isoformat(), selected=choice,
        criterion=plan['primary_selection'], validation_metrics=valid_metrics.to_dict('records'),
        artifacts={name: {p.name: base.sha(p) for p in (RESULT/name).iterdir()
                          if p.name in ['official_model.npz','scaling.json']} for name in models})
    base.save_json(RESULT/'frozen_validation_decision.json', frozen)
    valid.to_csv(RESULT/'validation_predictions.csv', index=False, float_format='%.17g')
    valid_metrics.to_csv(TABLE/'validation_evaluation.csv', index=False, float_format='%.17g')
    print(f'Frozen validation selection: {choice}. Now evaluating 2024 without refitting.', flush=True)
    external, external_metrics = evaluate(splits['external_2024'], models, settings, 'external_2024')
    external.to_csv(RESULT/'external_predictions.csv', index=False, float_format='%.17g')
    for name, saved in frozen['artifacts'].items():
        if any(base.sha(RESULT/name/file) != digest for file, digest in saved.items()):
            raise RuntimeError('Frozen fit/scaling artifact changed during external evaluation')
    for period in ['validation', 'external']:
        replay_data = pd.read_csv(RESULT/f'{period}_predictions.csv', float_precision='round_trip')
        for name in FEATURES:
            saved = replay_data.loc[replay_data.model.eq(name)]
            with np.load(RESULT/name/'official_model.npz') as archive:
                model = {key: archive[key] for key in base.MODEL_KEYS}
            scales = json.loads((RESULT/name/'scaling.json').read_text())
            replay = base.predict(model, saved[FEATURES[name]].to_numpy().T, scales, settings).T
            np.testing.assert_array_equal(replay, saved[[f'{v}_predicted' for v in base.OUTPUTS]].to_numpy())
    predictions = pd.concat([valid, external], ignore_index=True)
    evaluation = pd.concat([valid_metrics, external_metrics], ignore_index=True)
    evaluation.to_csv(TABLE/'evaluation.csv', index=False, float_format='%.17g')
    strata, events, event_audit, physical, maxima = rainfall_diagnostics(predictions, threshold)
    for name, frame in [('rain_strata',strata),('rain_events',events),('rain_event_coverage',event_audit),
                        ('physical_diagnostics',physical),('largest_observed_rain',maxima)]:
        frame.to_csv(TABLE/f'{name}.csv', index=False, float_format='%.17g')
    coefficients=[]
    for name, (model, scales, _) in models.items():
        physical_coef = model['coeff'].copy()
        physical_coef[1:] *= scales['Y']['scale'][None,:]/scales['X']['scale'][:,None]
        physical_coef[0] = scales['Y']['midpoint']+scales['Y']['scale']*model['coeff'][0]-scales['X']['midpoint']@physical_coef[1:]
        values = np.column_stack([np.ones(len(train)),train[FEATURES[name]].to_numpy()])@physical_coef
        np.testing.assert_allclose(values, base.predict(model,train[FEATURES[name]].to_numpy().T,scales,settings).T,
                                   rtol=1e-12,atol=1e-11)
        for i, term in enumerate(['intercept']+FEATURES[name]):
            for j, variable in enumerate(base.OUTPUTS):
                coefficients.append(dict(model=name, term=term, output=variable,
                    coefficient_scaled=model['coeff'][i,j], coefficient_physical=physical_coef[i,j]))
    pd.DataFrame(coefficients).to_csv(TABLE/'coefficients.csv', index=False, float_format='%.17g')
    # Historical scores are context only, never the controlled improvement estimate.
    historical = pd.read_csv(ROOT/'results/mimo_srs/external_predictions.csv', float_precision='round_trip')
    historical['timestamp_utc'] = pd.to_datetime(historical.timestamp_utc, utc=True)
    historical = historical.set_index('timestamp_utc').loc[splits['external_2024'].timestamp_utc]
    np.testing.assert_array_equal(historical[[f'{v}_observed' for v in base.OUTPUTS]].to_numpy(),
                                  splits['external_2024'][base.OUTPUTS].to_numpy())
    history_rows = [dict(model='historical_2019_2023_fit_context_only', variable=v,
        **base.metrics(historical[f'{v}_observed'], historical[f'{v}_predicted'])) for v in base.OUTPUTS]
    pd.DataFrame(history_rows).to_csv(TABLE/'historical_context.csv', index=False, float_format='%.17g')
    from C.mimo_history_figures import make_figures
    figures = make_figures(predictions, ROOT/'plots/mimo_rain_history')
    base.save_json(RESULT/'figures.json', figures)
    write_report(evaluation,counts,choice,threshold,physical,maxima,events,figures,models)
    preservation = verify_protected()
    state.update(status='completed', selected_on_validation=choice, sample_counts=counts.to_dict('records'),
        official_fits={k:v[2] for k,v in models.items()}, frozen_artifacts_unchanged=True,
        saved_predictions_reproduced_exactly=True, preservation=preservation,
        source_hash_unchanged=base.sha(source)==plan['source_sha256'],
        implementation_sha256={str(p.relative_to(base.REPO)):base.sha(p) for p in
            [Path(__file__), ROOT/'mimo_history_figures.py', ROOT/'run_fitting.py']})
    base.save_json(RESULT/'run_status.json', state)
    print(evaluation.loc[evaluation.variable.eq('P')].to_string(index=False), flush=True)
    return 0
