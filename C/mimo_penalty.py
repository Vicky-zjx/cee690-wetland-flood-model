"""Three initial positive penalties, unchanged official FIT and saved cohorts.

Run: python3 -B -m C.mimo_penalty
Tracing observes official Python frames; it never replaces or mutates a solver.
"""
from __future__ import annotations

from contextlib import chdir, redirect_stderr, redirect_stdout
from datetime import datetime, timezone
import inspect
import json
import os
from pathlib import Path
import re
import sys
import time

os.environ.setdefault('MPLBACKEND', 'Agg')
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parent/'results/mpl-cache'))
import numpy as np
import pandas as pd

from C import mimo_fit as base
from C import mimo_square as square
from C.mimo_history import rainfall_diagnostics

ROOT = base.ROOT
RESULT = ROOT/'results/mimo_rain_penalty'
TABLE = ROOT/'tables/mimo_rain_penalty'
PLOT = ROOT/'plots/mimo_rain_penalty'
CASES = {'l1_1': 1.0, 'l1_0p1': 0.1, 'l1_0p01': 0.01}
KKT_TOLERANCE = 1e-6  # Numerical diagnostic in scaled summed-loss gradient units.


def kkt_check(B, y, c, alpha):
    """Stationarity for SSE + alpha*sum(abs(c)), including the intercept."""
    c = np.asarray(c).reshape(-1)
    gradient = 2*B.T@(B@c-np.asarray(y).reshape(-1))
    active = c != 0
    residual = np.where(active, np.abs(gradient+alpha*np.sign(c)),
                        np.maximum(np.abs(gradient)-alpha, 0))
    implied = -gradient[active]*np.sign(c[active])
    return dict(kkt_max_abs_residual=float(residual.max()),
        kkt_pass=bool(residual.max() <= KKT_TOLERANCE),
        square_gradient=float(gradient[-1]),
        implied_alpha_min=float(implied.min()) if len(implied) else None,
        implied_alpha_max=float(implied.max()) if len(implied) else None)


class OfficialTrace:
    """Read-only line/return observation of the verified l1 function only."""
    def __init__(self, function):
        self.code = function.__code__
        lines, start = inspect.getsourcelines(function)
        def locate(fragment):
            matches = [start+i for i, line in enumerate(lines) if fragment in line]
            if len(matches) != 1:
                raise RuntimeError(f'Official trace location is ambiguous: {fragment}')
            return matches[0]
        self.branch = locate('if err_norm < err_norm_old:')
        self.zeroing = locate('c[np.abs(c) < 1e-6] = 0.0')
        self.calls = []
        self.current = None

    def __call__(self, frame, event, arg):
        if frame.f_code is not self.code:
            return None
        loc = frame.f_locals
        if event == 'call':
            self.current = {'steps': []}
        elif event == 'line' and frame.f_lineno == self.branch:
            accepted = bool(loc['err_norm'] < loc['err_norm_old'])
            self.current['steps'].append(dict(iteration=int(loc['iter_idx']),
                alpha_used=float(loc['alfa']), accepted=accepted,
                fit_error_candidate=float(loc['err_norm']),
                candidate_square=float(loc['c_new'][-1]),
                c_square_before_zeroing=None, c_square_returned=None))
        elif event == 'line' and frame.f_lineno == self.zeroing:
            step = self.current['steps'][-1]
            value = float(loc['c'][-1])
            step.update(alpha_next=float(loc['alfa']), c_square_before_zeroing=value,
                c_square_returned=0.0 if abs(value) < 1e-6 else value)
        elif event == 'return':
            if arg is None:
                # Do not replace an original solver exception with a tracing error.
                self.current['capture_failed'] = True
                return None
            c, mu, nu, history = arg
            update_stop = bool(np.linalg.norm(loc['u']) <= np.linalg.norm(loc['p'])/100
                and np.min(loc['p']) > -1e-4
                and np.linalg.norm(loc['v']) <= np.linalg.norm(loc['q'])/100
                and np.min(loc['q']) > -1e-4)
            reason = ('penalty_below_1e-4' if loc['alfa'] < 1e-4 else
                      'official_update_size_and_split_feasibility' if update_stop else
                      'iteration_limit' if loc['iter_idx'] == loc['max_iter']-1 else
                      'unclassified')
            # A floor exit occurs before zeroing/history storage; preserve its actual state.
            step = self.current['steps'][-1]
            if step['c_square_returned'] is None:
                step.update(alpha_next=float(loc['alfa']), c_square_returned=float(c[-1]))
            self.current.update(c=c.copy(), history=history.copy(),
                termination=reason, next_alpha=float(loc['alfa']),
                update_u_norm=float(np.linalg.norm(loc['u'])),
                update_v_norm=float(np.linalg.norm(loc['v'])),
                split_p_min=float(np.min(loc['p'])), split_q_min=float(np.min(loc['q'])))
            self.calls.append(self.current)
        return self


def run_case(name, penalty, train, scales, settings, module, provenance):
    out = RESULT/name
    out.mkdir(parents=True, exist_ok=True)
    X = base.transform(train[square.FEATURES].to_numpy().T, scales['X'])
    Y = base.transform(train[base.OUTPUTS].to_numpy().T, scales['Y'])
    order = np.random.RandomState(settings['seed']).permutation(len(train))
    n = int(np.floor(len(train)*settings['pTrain']/100))
    B = np.column_stack([np.ones(n), X[:, order[:n]].T])
    expected_membership = pd.read_csv(square.RESULT/'rain_square/training_membership.csv')
    training_times = pd.to_datetime(expected_membership.loc[
        expected_membership.internal_split.eq('coefficient_training'), 'timestamp_utc'], utc=True)
    np.testing.assert_array_equal(np.sort(train.timestamp_utc.iloc[order[:n]].to_numpy()),
                                  np.sort(training_times.to_numpy()))
    base.save_json(out/'scaling.json', scales)
    trace = OfficialTrace(module.l1)
    status = dict(status='started', official=provenance, initial_penalty=penalty,
        arguments={**settings, 'L1_pnlty':penalty}, supplied_samples=len(train),
        training_samples=n, diagnostic_samples=len(train)-n,
        capture='sys.settrace reads the unchanged official l1 frame; no solver replacement or mutation.')
    base.save_json(out/'execution.json', status)
    started = time.monotonic()
    old_rng, old_trace = np.random.get_state(), sys.gettrace()
    try:
        np.random.seed(settings['seed'])
        with (out/'official_stdout_stderr.log').open('w') as log, chdir(out), redirect_stdout(log), redirect_stderr(log):
            sys.settrace(trace)
            try:
                returned = module.mimo_srs(X, Y, max_order=1, pTrain=70,
                    scaling=np.array([0,0]), L1_pnlty=penalty, basis_fctn='P',
                    var_names={'X':square.FEATURES, 'Y':base.OUTPUTS})
            finally:
                sys.settrace(old_trace)
        counts = [int(v) for v in re.findall(r'^\s*(\d+) data values\s*$',
            (out/'official_stdout_stderr.log').read_text(), re.MULTILINE)]
        if counts != [len(train), len(train)] or len(trace.calls) != 4:
            raise RuntimeError(f'Official filtering/capture mismatch: {counts}; {len(trace.calls)} calls')
        np.testing.assert_array_equal(returned[6], X[:,order[n:]])
        np.testing.assert_array_equal(returned[7], Y[:,order[n:]])
        if not all(np.isfinite(a).all() for a in returned):
            raise FloatingPointError('Nonfinite official result')
        model = dict(zip(base.MODEL_KEYS, returned[:6]))
        replay, _ = module.compute_model(*returned[:6], returned[6], np.array([0,0]), 'P')
        np.testing.assert_allclose(replay, returned[8], rtol=0, atol=1e-13)
        np.savez(out/'official_model.npz', **model)
        rows, iterations = [], []
        for j, (output, captured) in enumerate(zip(base.OUTPUTS, trace.calls)):
            np.testing.assert_array_equal(captured['c'], model['coeff'][:,j])
            accepted = [s for s in captured['steps'] if s['accepted']]
            if not accepted:
                raise RuntimeError('No accepted official step; needs explicit treatment')
            last_alpha = accepted[-1]['alpha_used']
            check = kkt_check(B, Y[j,order[:n]], captured['c'], last_alpha)
            zero_steps = [s for s in captured['steps'] if s['c_square_returned'] == 0]
            row = dict(model=name, output=output, initial_penalty=penalty,
                square_coefficient_scaled=float(captured['c'][-1]),
                initial_ols_square=float(captured['history'][len(captured['c'])-1,0]),
                last_accepted_alpha=last_alpha, next_alpha=captured['next_alpha'],
                termination=captured['termination'], iterations=len(captured['steps']),
                accepted_steps=len(accepted), first_square_zero_iteration=zero_steps[0]['iteration'] if zero_steps else None,
                first_zero_prethreshold_value=zero_steps[0]['c_square_before_zeroing'] if zero_steps else None,
                active_coefficients=int(np.count_nonzero(captured['c'])), **check)
            rows.append(row)
            iterations.extend(dict(model=name, output=output, **s) for s in captured['steps'])
            np.savez(out/f'{output}_official_history.npz', history=captured['history'])
        status.update(status='completed', retained_samples=counts[1],
            termination_by_output={r['output']:r['termination'] for r in rows},
            numerically_screened_rain_fit=rows[0]['kkt_pass'])
        return model, rows, iterations
    except BaseException as exc:
        status.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        sys.settrace(old_trace)
        np.random.set_state(old_rng)
        status['elapsed_seconds'] = time.monotonic()-started
        base.save_json(out/'execution.json', status)
        import matplotlib.pyplot as plt
        plt.close('all'); plt.ioff()


def evaluate(pairs, models, scales, settings, period):
    from C.uncertainty.weather_generator import solar_geometry
    _, daylight = solar_geometry(pd.DatetimeIndex(pairs.timestamp_utc)-pd.Timedelta(hours=1),43.0547222,-88.0461111)
    frames, rows = [], []
    original, original_scales = square.read_model(square.SOURCE/'one_hour')
    for name in [*CASES, 'one_hour', 'persistence']:
        if name == 'persistence':
            values = pairs[base.PREDICTORS[:4]].to_numpy()
        elif name == 'one_hour':
            values = base.predict(original, pairs[base.PREDICTORS].to_numpy().T, original_scales, settings).T
        else:
            values = base.predict(models[name], pairs[square.FEATURES].to_numpy().T, scales, settings).T
        f = pairs.copy(); f['period'] = period; f['model'] = name; f['night'] = ~daylight
        for j, var in enumerate(base.OUTPUTS):
            f[f'{var}_observed'] = pairs[var].to_numpy()
            f[f'{var}_predicted'] = values[:,j]
            f[f'{var}_residual'] = pairs[var].to_numpy()-values[:,j]
            rows.append(dict(period=period, model=name, variable=var, units=base.UNITS[j],
                **base.metrics(pairs[var].to_numpy(), values[:,j])))
        frames.append(f)
    return pd.concat(frames,ignore_index=True), pd.DataFrame(rows)


def validation_choice(metrics, diagnostics):
    scores = metrics.loc[metrics.period.eq('validation_2023') & metrics.variable.eq('P')
        & metrics.model.isin(CASES)].set_index('model').RMSE
    raw = min(CASES, key=lambda name:scores[name])
    passed = diagnostics.loc[diagnostics.output.eq('P') & diagnostics.kkt_pass, 'model'].tolist()
    selected = min(passed, key=lambda name:(scores[name], -CASES[name])) if passed else None
    return dict(lowest_validation_RMSE_descriptive=raw, selected_after_numerical_screen=selected,
        recommendation='No numerically screened replacement' if selected is None else
            'Provisional rainfall-only selection; inspect all-output and physical diagnostics',
        criterion='2023 rainfall RMSE among cases passing the predeclared rain KKT diagnostic; ties keep larger initial penalty.',
        numerical_rule=f'Maximum stationarity residual <= {KKT_TOLERANCE} for SSE + alpha*L1 at the last accepted-step alpha, all scaled coefficients including intercept.',
        caveat='Diagnostic added here, not a convergence certificate from official mimo_srs; all raw outputs remain reported.')


def verify_protection():
    old = json.loads((RESULT/'protected_files_before.json').read_text())
    changed = [p for p,h in old.items() if not (base.REPO/p).is_file() or base.sha(base.REPO/p) != h]
    report = dict(checked_files=len(old), changed_or_missing=changed)
    base.save_json(RESULT/'preservation_check.json', report)
    if changed:
        raise RuntimeError(f'Protected files changed: {changed}')
    return report


def run():
    for p in [RESULT,TABLE,PLOT]:
        p.mkdir(parents=True, exist_ok=True)
    if not (RESULT/'protected_files_before.json').exists():
        # Reproduction checks all previous source/model artifacts used by this run.
        paths = [*square.SOURCE.glob('*_pairs.csv'), * (square.RESULT/'rain_square').glob('*')]
        base.save_json(RESULT/'protected_files_before.json', {
            p.relative_to(base.REPO).as_posix():base.sha(p) for p in paths if p.is_file()})
    settings = base.DEFAULT_SETTINGS.copy()
    module, provenance = base.official_module()
    old_status = json.loads((square.RESULT/'rain_square/execution.json').read_text())
    if provenance != old_status['official']:
        raise RuntimeError('Official source provenance changed')
    original, scales = square.read_model(square.RESULT/'rain_square')
    periods = square.load_comparison_pairs()
    train = periods['fit_2019_2022']
    for role, columns in [('X',square.FEATURES),('Y',base.OUTPUTS)]:
        expected = base.fit_affine(train[columns].to_numpy().T, settings['external_bound'])
        for key in expected:
            np.testing.assert_array_equal(expected[key], scales[role][key])
    base.save_json(RESULT/'study_plan.json', dict(created_utc=datetime.now(timezone.utc).isoformat(),
        cases=CASES, settings=settings, predictors=square.FEATURES, sample_counts=square.PERIOD_COUNTS,
        official=provenance, kkt_tolerance=KKT_TOLERANCE,
        scope='Only initial positive L1 penalty varies. Initial 1.0 replay captures absent histories and must reproduce saved model exactly.',
        selection='Lowest 2023 P RMSE among rain fits passing KKT at last accepted-step alpha; no replacement if none passes.',
        external='Freeze coefficients and decision before evaluating 2024; already inspected years remain exploratory, not a new blind test.',
        zero_penalty='Not run: official zero branch uses incompatible B and untransposed Zy dimensions. No source fix is authorized here.'))
    base.save_json(RESULT/'run_status.json', dict(status='started'))
    models, all_diag, all_steps = {}, [], []
    for name, penalty in CASES.items():
        print(f'Running unchanged official mimo_srs: {name}',flush=True)
        model, diagnostic, iterations = run_case(name, penalty, train, scales, settings, module, provenance)
        if name == 'l1_1':
            for key in base.MODEL_KEYS:
                np.testing.assert_array_equal(model[key],original[key])
        models[name]=model; all_diag.extend(diagnostic); all_steps.extend(iterations)
    diagnostics, iterations = pd.DataFrame(all_diag), pd.DataFrame(all_steps)
    diagnostics.to_csv(TABLE/'solver_diagnostics.csv',index=False,float_format='%.17g')
    iterations.to_csv(TABLE/'solver_iterations.csv',index=False,float_format='%.17g')
    validation, vm = evaluate(periods['validation_2023'],models,scales,settings,'validation_2023')
    decision = validation_choice(vm,diagnostics)
    frozen = {str(p.relative_to(base.REPO)):base.sha(p) for name in CASES for p in
        [RESULT/name/'official_model.npz',RESULT/name/'scaling.json']}
    base.save_json(RESULT/'frozen_validation_decision.json',dict(frozen_utc=datetime.now(timezone.utc).isoformat(),
        **decision, validation_metrics=vm.to_dict('records'), frozen_artifacts=frozen))
    print(f'Frozen validation decision: {decision}',flush=True)
    external, em = evaluate(periods['external_2024'],models,scales,settings,'external_2024')
    for name, frame in [('validation',validation),('external',external)]:
        frame.to_csv(RESULT/f'{name}_predictions.csv',index=False,float_format='%.17g')
        saved = pd.read_csv(RESULT/f'{name}_predictions.csv',float_precision='round_trip')
        old = pd.read_csv(square.RESULT/f'{name}_predictions.csv',float_precision='round_trip')
        fields = ['timestamp_utc']+[f'{v}_{k}' for v in base.OUTPUTS for k in ['observed','predicted']]
        for current, previous in [('l1_1','rain_square'),('one_hour','one_hour'),('persistence','persistence')]:
            pd.testing.assert_frame_equal(saved.loc[saved.model.eq(current),fields].reset_index(drop=True),
                old.loc[old.model.eq(previous),fields].reset_index(drop=True),check_exact=True)
        for case in CASES:
            f = saved.loc[saved.model.eq(case)]
            model, saved_scales = square.read_model(RESULT/case)
            pred = base.predict(model,f[square.FEATURES].to_numpy().T,saved_scales,settings).T
            np.testing.assert_array_equal(pred,f[[f'{v}_predicted' for v in base.OUTPUTS]].to_numpy())
    predictions = pd.concat([validation,external],ignore_index=True)
    metrics = pd.concat([vm,em],ignore_index=True)
    metrics.to_csv(TABLE/'evaluation.csv',index=False,float_format='%.17g')
    threshold = json.loads((square.SOURCE/'diagnostic_definition.json').read_text())['heavy_rain_threshold_mm_per_h']
    products = rainfall_diagnostics(predictions,threshold)
    for name, frame in zip(['rain_strata','rain_events','rain_event_coverage','physical_diagnostics','largest_observed_rain'],products):
        frame.to_csv(TABLE/f'{name}.csv',index=False,float_format='%.17g')
    coefficients=[]
    for name, model in models.items():
        c = model['coeff'].copy()
        c[1:] *= np.asarray(scales['Y']['scale'])[None,:]/np.asarray(scales['X']['scale'])[:,None]
        c[0] = np.asarray(scales['Y']['midpoint'])+np.asarray(scales['Y']['scale'])*model['coeff'][0]-np.asarray(scales['X']['midpoint'])@c[1:]
        for i, term in enumerate(['intercept']+square.FEATURES):
            for j, var in enumerate(base.OUTPUTS):
                coefficients.append(dict(model=name,term=term,output=var,coefficient_scaled=model['coeff'][i,j],coefficient_physical=c[i,j]))
    pd.DataFrame(coefficients).to_csv(TABLE/'coefficients.csv',index=False,float_format='%.17g')
    if any(base.sha(base.REPO/p)!=h for p,h in frozen.items()):
        raise RuntimeError('Frozen model changed during external evaluation')
    from C.mimo_penalty_figures import make_figures
    figures=make_figures(predictions,metrics,iterations,diagnostics,PLOT)
    base.save_json(RESULT/'figures.json',figures)
    protection=verify_protection()
    base.save_json(RESULT/'run_status.json',dict(status='completed',official=provenance,
        original_square_model_replayed_exactly=True,saved_predictions_replayed_exactly=True,
        original_reference_predictions_unchanged=True,sample_counts=square.PERIOD_COUNTS,
        validation_decision=decision,preservation=protection,
        implementation_sha256={p.name:base.sha(p) for p in [Path(__file__),ROOT/'mimo_penalty_figures.py']}))
    print(metrics.loc[metrics.variable.eq('P')].to_string(index=False),flush=True)
    print(diagnostics.loc[diagnostics.output.eq('P')].to_string(index=False),flush=True)


if __name__ == '__main__':
    run()
