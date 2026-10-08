"""Direct analysis(v,c) checks at wetland design bounds; never runs OPT.

Run from the project root: python3 -B -m B.wetland_analysis_checks --run --report
Existing completed checks and all earlier numerical studies are preserved.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import time
import numpy as np
from B import wetland_only_study as study
from B import optimize_wetland as old
from B.multistart_check import save_trajectory
from B.no_spillway_study import load_trajectory
from B.official_opt import load_official

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results' / 'analysis_checks'
TABLE = ROOT / 'tables' / 'analysis_checks' / 'boundary_checks.md'
REPORT = ROOT / 'docs' / 'analysis_and_boundary_checks.md'
PLOTS = ROOT / 'plots' / 'analysis_checks'
CORNERS = ((0., 0.), (0., 1.), (1., 0.), (1., 1.))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def protected_files(config):
    """Scientific inputs and earlier outputs; no historical hashes are rewritten."""
    paths = list((ROOT.parent / 'A').glob('*.py'))
    paths += list((ROOT.parent / 'A/plots').glob('*.png'))
    paths += [p for p in (ROOT.parent / 'A/results').rglob('*') if p.is_file()]
    paths += [p for p in (ROOT / 'results').rglob('*')
              if p.is_file() and OUT not in p.parents and p.name != '.DS_Store']
    paths += [ROOT / name for name in ('no_spillway_config.json', 'wetland_only_config.json',
              'no_spillway_model.py', 'gate_hydraulics.py', 'official_opt.py',
              'optimize_wetland.py', 'wetland_only_study.py', 'wetland_analysis_checks.py')]
    paths += [old.resolve_path(config['forcing'][key]) for key in ('prepared_csv', 'manifest_json')]
    manifest = old.load_config(config['forcing']['manifest_json'])
    paths += [old.resolve_path(row['path']) for row in manifest['source_files']]
    return {str(p.relative_to(ROOT.parent)): digest(p) for p in paths}


def run_checks():
    if OUT.exists() and any(OUT.iterdir()):
        raise RuntimeError('Analysis checks already exist; use --report to read them, without overwriting.')
    config, spec, gates = study.load_study()
    official = load_official()
    forcing = study.model.fixed_forcing(config, study.model.base_parameters(config))
    before = protected_files(config)
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    tolerance = config['numerics']['peak_grid_tolerance_m']
    plan = dict(scope='Direct public analysis function verification; no optimization',
        corners=CORNERS, category_pairs=list(itertools.product(*config['categories']['admissible'])),
        primary_timesteps_s=[120, 60], optional_timestep_s=30,
        optional_refinement_rule='Run 30 s if validity/feasibility differs between 120 and 60 s, or their physical peak difference exceeds the existing peak_grid_tolerance_m.',
        peak_grid_tolerance_m=tolerance, tolerance_role='Numerical diagnostic only; no physical limit is changed.',
        configuration=config, specification=spec, fixed_gate_knots=gates,
        official_revision=official['revision'], official_version=official['version'],
        objective='Physical max H_r; callback returns (H_peak-objective_origin_m)/height_scale_m',
        optimizer_constraint_count=18, G5='Separate unchanged RHS/RK applicability check',
        protected_sha256=before, optimizer_runs=0)
    old.save_json(OUT / 'plan.json', plan)
    reproduction = []
    # Test the public function against the saved full reference and best arrays.
    for label, name in [('reference', 'fixed_reference_120s.npz'), ('best', 'wetland_best_120s.npz')]:
        saved = load_trajectory(ROOT / 'results/wetland_only' / name)
        c = study.AnalysisConstants(config, gates, tuple(saved['pair']), forcing)
        f, g = study.analysis(saved['theta'][-2:], c)
        run = c.last_evaluation['run']
        if run is None:
            raise RuntimeError('Saved valid case became invalid: ' + label)
        errors = {k: float(np.max(np.abs(run[k] - saved[k]))) for k in ['time', 'x', 'u', 'y', 'balance']}
        expected_g = study.optimization_view(saved)['residuals']
        expected_f = (saved['metrics']['peak_stage_m'] - config['numerics']['objective_origin_m']) / config['numerics']['height_scale_m']
        if any(value != 0. for value in errors.values()) or f != expected_f or not np.array_equal(g, expected_g):
            raise RuntimeError('Extracted analysis changed the saved evaluation: ' + str(errors))
        reproduction.append(dict(case=label, max_array_errors=errors, same_f_and_g=True,
            physical_peak_m=run['metrics']['peak_stage_m']))
        print('Exact public-API reproduction:', label, flush=True)
    records = []

    def evaluate(pair, covers, dt):
        c = study.AnalysisConstants(config, gates, tuple(pair), forcing, refine=120 // dt)
        tick = time.monotonic()
        f, g = study.analysis(covers, c)
        detail = c.last_evaluation
        row = dict(pair=pair, covers=covers, timestep_s=dt, finite_callback=True,
            conditioned_objective=f, optimizer_residuals=g, valid=detail['valid'],
            feasible=detail['feasible'], physical_peak_m=detail['physical_peak_m'],
            failure=detail['failure'], diagnostic=detail['diagnostic'], seconds=time.monotonic()-tick,
            status='invalid_model_evaluation' if not detail['valid'] else
                   ('valid_feasible' if detail['feasible'] else 'valid_constraint_violation'))
        if detail['valid']:
            run, view = detail['run'], detail['view']
            assert run['time'][-1] == config['numerics']['duration_s']
            assert np.array_equal(run['theta'][:-2], gates)
            assert np.array_equal(run['x'][:, 0], config['initial_conditions']['x0_m3'])
            np.testing.assert_array_equal(run['u'][2:], np.vstack([
                np.interp(run['time'], forcing[0], channel) for channel in forcing[1]]))
            violated = [a['label'] for a in view['active_audit'] + view['native_cover_audit']
                        if not a['pass_at_tolerance']]
            name = '_'.join(pair) + f'_rho{int(covers[0])}{int(covers[1])}_{dt}s.npz'
            save_trajectory(OUT / name, run)
            row.update(metrics=run['metrics'], active_audit=view['active_audit'],
                applicability_audit=view['applicability_audit'], native_cover_audit=view['native_cover_audit'],
                violated_constraints=violated, max_scaled_residual=float(max(g)),
                min_G5_margin_m=float(min(a['margin'] for a in view['applicability_audit'])),
                trajectory='B/results/analysis_checks/' + name,
                residual_kind='computed_physical_residuals')
        else:
            row.update(violated_constraints=None, max_scaled_residual=None,
                residual_kind='finite_invalid_integration_sentinel_not_physical_residuals')
        records.append(row)
        old.save_json(OUT / 'progress.json', dict(records=records))
        print('/'.join(pair), covers, dt, row['status'], row['physical_peak_m'], flush=True)
        return row

    comparisons = []
    for pair in plan['category_pairs']:
        for covers in CORNERS:
            coarse, fine = (evaluate(pair, covers, dt) for dt in (120, 60))
            same_status = coarse['status'] == fine['status']
            delta = abs(coarse['physical_peak_m'] - fine['physical_peak_m']) if coarse['valid'] and fine['valid'] else None
            needs30 = not same_status or (delta is not None and delta > tolerance)
            comparison = dict(pair=pair, covers=covers, same_status_120_60=same_status,
                peak_change_120_60_m=delta, thirty_second_required=needs30)
            if needs30:
                finer = evaluate(pair, covers, 30)
                comparison.update(status_30s=finer['status'], same_status_60_30=fine['status']==finer['status'],
                    peak_change_60_30_m=abs(fine['physical_peak_m']-finer['physical_peak_m']) if fine['valid'] and finer['valid'] else None)
            comparisons.append(comparison)
    after = protected_files(config)
    if before != after:
        raise RuntimeError('A protected input, source or earlier result changed during checks.')
    result = dict(records=records, timestep_comparisons=comparisons, reproduction_checks=reproduction,
        actual_ode_attempts=2+len(records), optimizer_runs=0, protected_files_unchanged=True,
        seconds=time.monotonic()-started, completion='All planned direct boundary evaluations completed; individual validity and feasibility are reported separately.')
    old.save_json(OUT / 'checks.json', result)
    print(result['completion'], flush=True)



def plot_checks(result, plan):
    """Plot saved time histories only; no new analysis, ODE or OPT calls."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    pairs = [tuple(pair) for pair in plan['category_pairs']]
    labels = ['/'.join(x.removeprefix('type_').upper() for x in pair) for pair in pairs]
    corners = [tuple(x) for x in plan['corners']]
    cover_labels = [f'ρ = ({a:g}, {b:g})' for a, b in corners]
    records = {(tuple(row['pair']), tuple(row['covers']), row['timestep_s']): row
               for row in result['records']}
    expected = {(pair, corner, dt) for pair in pairs for corner in corners for dt in (120, 60)}
    if not expected <= records.keys():
        raise ValueError('A planned boundary time series is missing.')
    runs = {}
    for key in expected:
        record = records[key]
        if not record['valid']:
            raise ValueError('Invalid trajectories cannot be drawn as complete physical histories.')
        run = load_trajectory(old.resolve_path(record['trajectory']))
        assert tuple(run['pair']) == key[0] and tuple(run['theta'][-2:]) == key[1]
        assert np.max(run['y'][old.a2.Y['H_r']]) == record['physical_peak_m']
        runs[key] = run
    fine_time = runs[pairs[0], corners[0], 60]['time']
    coarse_time = runs[pairs[0], corners[0], 120]['time']
    assert np.array_equal(coarse_time, fine_time[::2])
    for (pair, corner, dt), run in runs.items():
        assert np.array_equal(run['time'], fine_time if dt == 60 else coarse_time)
    fine_hours, coarse_hours = fine_time/3600, coarse_time/3600
    PLOTS.mkdir(parents=True, exist_ok=True)
    paths = []
    colors = ['#2166ac', '#db7c25', '#238b74', '#9659a1']
    linestyles = ['-', '--', '-.', ':']
    styles = list(zip(colors, linestyles))
    handles = [Line2D([0], [0], color=color, ls=ls, lw=2, label=label)
               for (color, ls), label in zip(styles, cover_labels)]

    def decorate(ax, ylabel=None):
        ax.set_xlim(0, 144)
        ax.set_xticks(np.arange(0, 145, 24))
        ax.grid(alpha=.19, lw=.65)
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines[['bottom', 'left']].set_color('#b8c5ca')
        ax.tick_params(colors='#294553')
        ax.ticklabel_format(axis='y', style='plain', useOffset=False)
        if ylabel:
            ax.set_ylabel(ylabel)

    def save(fig, name):
        path = PLOTS/name
        fig.savefig(path, dpi=180, facecolor='white', bbox_inches='tight', pad_inches=.22)
        plt.close(fig)
        paths.append(path)

    style = {'font.family': 'DejaVu Sans', 'font.size': 10.5,
             'axes.titleweight': 'bold', 'text.color': '#213b4b',
             'axes.labelcolor': '#213b4b', 'axes.edgecolor': '#bac8ce',
             'legend.frameon': False}
    with plt.rc_context(style):
        fig, axes = plt.subplots(2, 2, figsize=(12.2, 9), sharex=True, sharey=True)
        fig.subplots_adjust(left=.09, right=.975, top=.78, bottom=.16, hspace=.27, wspace=.18)
        fig.text(.09, .955, 'Stage B | Boundary river-stage histories', fontsize=19, weight='bold')
        fig.text(.09, .910, '60 s saved trajectories · fixed gates and event forcing · each panel shows one class pair', fontsize=11)
        fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.52, .88), ncol=4, columnspacing=2.5)
        for ax, pair, label in zip(axes.flat, pairs, labels):
            for corner, (color, ls) in zip(corners, styles):
                run = runs[pair, corner, 60]
                stages = run['y'][old.a2.Y['H_r']]
                ax.plot(fine_hours, stages, color=color, ls=ls, lw=1.7)
                k = int(np.argmax(stages))
                ax.plot(fine_hours[k], stages[k], marker='o', color=color, ms=4, linestyle='none')
            ax.set_title('Response classes '+label, loc='left', pad=10)
            decorate(ax)
        for ax in axes[:, 0]:
            ax.set_ylabel('River level Hᵣ (m, assumed datum)')
        for ax in axes[1]:
            ax.set_xlabel('Elapsed time from 2022-09-09 00:00 UTC (h)')
        fig.text(.09, .08,
            'Dots mark each simulated peak. Lines join saved time samples; exact peaks remain in the boundary table.\n'
            'A/B are synthetic classes. At zero cover the wetland remains and its class has no effect; 16 labeled cases represent 9 physical designs.',
            fontsize=9.5, va='top', linespacing=1.6)
        save(fig, '01_boundary_river_peaks.png')

        config = plan['configuration']
        limits = [*config['safety']['dam_m'], *config['safety']['wetland_m']]
        stage_labels = ['H_d1', 'H_d2', 'H_w1', 'H_w2']
        safety = {}
        gate_histories = []
        for pair in pairs:
            for corner in corners:
                run = runs[pair, corner, 60]
                margins = np.array([limit-run['y'][old.a2.Y[name]]
                                    for name, limit in zip(stage_labels, limits)])
                saved = {a['label']: a for a in records[pair, corner, 60]['active_audit']}
                for k, name in enumerate(stage_labels):
                    assert np.isclose(margins[k].min(), saved[name+'_safe']['margin'], rtol=0, atol=1e-12)
                safety[pair, corner] = margins.min(axis=0)
                heads = np.array([run['y'][old.a2.Y[f'H_d{i+1}']]-config['gates']['gate_sill_z'][i] for i in range(2)])
                openings = np.asarray(config['gates']['gate_max_opening'])[:, None]*run['u'][:2]
                tailwater = np.array([run['y'][old.a2.Y[f'H_w{i+1}']] for i in range(2)])
                donor = np.array([run['y'][old.a2.Y[f'H_d{i+1}']] for i in range(2)])
                if not np.all((heads > 0) & (openings > 0) & (donor > tailwater)):
                    raise ValueError('G5 plot assumes the recorded gates are active at all displayed nodes.')
                margin = heads-1.25*openings
                saved_g5 = {a['label']: a for a in records[pair, corner, 60]['applicability_audit']}
                for i in range(2):
                    assert np.isclose(margin[i].min(), saved_g5[f'gate{i+1}_head_opening']['margin'], rtol=0, atol=1e-12)
                gate_histories.append(margin)
        if not all(np.array_equal(gate_histories[0], g) for g in gate_histories[1:]):
            raise ValueError('G5 trajectories differ across cases; do not label them as shared.')
        fig = plt.figure(figsize=(12.2, 11.2))
        grid = fig.add_gridspec(3, 2, left=.09, right=.975, top=.815, bottom=.15,
                               height_ratios=[1, 1, .9], hspace=.40, wspace=.18)
        axes = [fig.add_subplot(grid[i//2, i % 2]) for i in range(4)]
        gax = fig.add_subplot(grid[2, :])
        fig.text(.09, .964, 'Stage B | Safety headroom through the event', fontsize=19, weight='bold')
        fig.text(.09, .925, 'Each upper curve is the smallest headroom to the four upstream level limits at that time.', fontsize=11)
        fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.52, .90), ncol=4, columnspacing=2.5)
        upper = max(np.max(v) for v in safety.values())
        for ax, pair, label in zip(axes, pairs, labels):
            for corner, (color, ls) in zip(corners, styles):
                ax.plot(fine_hours, safety[pair, corner], color=color, ls=ls, lw=1.6)
            ax.axhline(0, color='#ab3b42', ls='--', lw=1.1)
            ax.set_ylim(-.035, upper*1.10)
            ax.set_title('Response classes '+label, loc='left', pad=8)
            decorate(ax, 'Minimum safety headroom (m)')
        for i, color in enumerate(('#2166ac', '#db7c25')):
            gax.plot(fine_hours, gate_histories[0][i], color=color, ls=('-', '--')[i],
                     lw=1.65, label=f'Gate {i+1}')
        gax.axhline(0, color='#ab3b42', ls='--', lw=1.1, label='Zero headroom')
        gax.set_title('G5 applicability: the two gate curves are identical across all 16 cases', loc='left', pad=8)
        gax.set_ylim(bottom=-.035)
        decorate(gax, 'G5 head margin (m)')
        gax.set_xlabel('Elapsed time from 2022-09-09 00:00 UTC (h)')
        gax.legend(loc='upper right', ncol=3, fontsize=9)
        fig.text(.09, .081,
            'Safety headroom = min(H_limit − H) over Pool 1, Pool 2, Wetland 1 and Wetland 2 at each saved 60 s node.\n'
            'Positive headroom satisfies all four sampled level limits. G5 margin = depth above sill − 1.25 × gate opening height.\n'
            'G5 is an RHS/RK applicability guard, separate from the 18 OPT residuals; other constraints remain in the full saved audit.',
            fontsize=9.3, va='top', linespacing=1.6)
        save(fig, '02_boundary_safety_margins.png')

        differences, balance120, balance60 = [], [], []
        for pair in pairs:
            for corner in corners:
                coarse, fine = runs[pair, corner, 120], runs[pair, corner, 60]
                differences.append(np.abs(coarse['y'][old.a2.Y['H_r']] - fine['y'][old.a2.Y['H_r'], ::2])*1e6)
                balance120.append(np.abs(coarse['balance']))
                balance60.append(np.abs(fine['balance']))
                for dt, run in ((120, coarse), (60, fine)):
                    assert np.max(np.abs(run['balance'])) == records[pair, corner, dt]['metrics']['balance_max_m3']
        discrepancy = np.max(differences, axis=0)
        b120, b60 = np.max(balance120, axis=0), np.max(balance60, axis=0)
        peak_change = max(row['peak_change_120_60_m'] for row in result['timestep_comparisons'])*1000
        fig, axes = plt.subplots(2, 1, figsize=(12.2, 8.6), sharex=True)
        fig.subplots_adjust(left=.11, right=.97, top=.81, bottom=.19, hspace=.38)
        fig.text(.11, .957, 'Stage B | Numerical differences through time', fontsize=19, weight='bold')
        fig.text(.11, .905, 'Curves show the largest discrepancy across the 16 saved cases at each time.', fontsize=11)
        axes[0].plot(coarse_hours, discrepancy, color='#2166ac', lw=1.8)
        axes[0].set_title('River-stage difference at matching time nodes', loc='left', pad=10)
        decorate(axes[0], 'max |Hᵣ,120(t) − Hᵣ,60(t)| (µm)')
        axes[0].set_ylim(bottom=0)
        axes[0].text(.02, .91, f'Maximum at common nodes: {discrepancy.max():.6f} µm',
                     transform=axes[0].transAxes, va='top', fontsize=10)
        axes[1].plot(coarse_hours, b120, color='#2166ac', lw=1.8, label='120 s grid')
        axes[1].plot(fine_hours, b60, color='#db7c25', ls='--', lw=1.8, label='60 s grid')
        axes[1].set_title('Cumulative water-balance residual: maximum absolute value across cases', loc='left', pad=10)
        decorate(axes[1], 'max |balance residual(t)| (m³)')
        axes[1].set_ylim(bottom=0)
        axes[1].set_xlabel('Elapsed time from 2022-09-09 00:00 UTC (h)')
        axes[1].legend(loc='upper left', ncol=2)
        fig.text(.11, .107,
            'Top: shared 120 s time nodes only, with no interpolation. Extra 60 s samples can capture a peak between those nodes.\n'
            f'Maximum change between separately sampled peaks = {peak_change:.6f} mm; this differs from the matched-time curve above.\n'
            'These envelopes summarize saved simulations, not new trajectories. Finite-grid agreement is not a field-accuracy certificate.',
            fontsize=9.4, va='top', linespacing=1.6)
        save(fig, '03_boundary_numerical_checks.png')
    return paths


def report():
    result, plan = (old.load_config(OUT / name) for name in ('checks.json', 'plan.json'))
    rows = result['records']
    valid = [r for r in rows if r['valid']]
    feasible = [r for r in valid if r['feasible']]
    comparable = [r['peak_change_120_60_m'] for r in result['timestep_comparisons'] if r['peak_change_120_60_m'] is not None]
    table = ['# Direct wetland coverage-bound checks', '',
        'A/B are synthetic response classes. All values below come from the public `analysis(v,c)` function and official `ode4u`; no optimizer is run. Negative physical/scaled inequality residuals satisfy their corresponding limit. Feasibility uses the existing declared numerical tolerance. The integration-status residual is -1 for a valid completed run. The maximum residual also includes the storage constraints, so a value close to zero can indicate a nearly empty storage rather than an upstream safety limit.', '',
        '| Pair | Covers | Step (s) | Status | Physical river peak (m) | Max scaled residual | Max balance error (m³) | Violated constraints |',
        '|---|---|---:|---|---:|---:|---:|---|']
    for r in rows:
        pair='/'.join(x.replace('type_', '').upper() for x in r['pair'])
        peak='N/A' if r['physical_peak_m'] is None else f"{r['physical_peak_m']:.9f}"
        residual='N/A' if not r['valid'] else f"{r['max_scaled_residual']:.6g}"
        balance='N/A' if not r['valid'] else f"{r['metrics']['balance_max_m3']:.6g}"
        violations=', '.join(r['violated_constraints'] or []) if r['valid'] else r['failure']
        table.append(f"| {pair} | {r['covers'][0]:g}, {r['covers'][1]:g} | {r['timestep_s']} | {r['status']} | {peak} | {residual} | {balance} | {violations or 'None'} |")
    table += ['', '## Timestep comparison', '',
        '| Pair | Covers | Absolute peak change, 120 to 60 s (mm) | Same validity/feasibility classification | 30 s needed |',
        '|---|---|---:|---|---|']
    for r in result['timestep_comparisons']:
        delta='N/A' if r['peak_change_120_60_m'] is None else f"{1000*r['peak_change_120_60_m']:.6f}"
        table.append(f"| {'/'.join(r['pair'])} | {r['covers']} | {delta} | {r['same_status_120_60']} | {r['thirty_second_required']} |")
    TABLE.parent.mkdir(parents=True, exist_ok=True)
    TABLE.write_text('\n'.join(table)+'\n')
    lines = ['# Stage B analysis function and design-bound verification', '',
        'The final study fixes gate operations and varies only wetland designs. Each category pair is held in the constants; the continuous vector is `v=[rho1,rho2]`. The four synthetic pairs A/A, A/B, B/A and B/B are enumerated outside the continuous search. This validation follows the professor\'s recommendation to call the analysis on distinct and extreme designs and check the design bounds. Testing all four two-dimensional corners is this project\'s explicit coverage choice.', '',
        '## Callable analysis', '',
        'The implementation is [wetland_only_study.py](../wetland_only_study.py). Both standalone calls and the optimizer callback use `analysis(v,c)` with the same simulation and constraint implementation. No optimizer algorithm was added.', '',
        '```python',
        'from B.wetland_only_study import make_analysis_constants, analysis',
        'c = make_analysis_constants(("type_a", "type_b"))',
        'f, g = analysis([0.0, 1.0], c)',
        'print(f, g)  # dimensionless scalar and 1-D array of 18 residuals',
        'print(c.last_evaluation["physical_peak_m"])',
        'print(c.last_evaluation["valid"], c.last_evaluation["feasible"])', '```', '',
        '`c` contains the fixed model configuration, readonly copy of the gate schedule, category pair, forcing and timestep refinement. `last_evaluation` is diagnostic output replaced on every call; it is not a simulation input. The call does not run OPT or write files.', '',
        'The physical objective is `max_t H_r(t)` in metres, evaluated after ODE integration. The returned objective is the unchanged numerical conditioning `f=(H_peak-97 m)/(1 m)`. Thus `f` and all 18 entries of `g` are dimensionless. Water flows remain m³/s, storage m³ and time s; volume constraints use their fixed m³ scales and level constraints their fixed 1 m scale. These numerical scales do not change physical safety levels. Native cover bounds remain [0,1]. G5 stays an RHS/RK applicability check outside `g`.', '',
        'A valid but unsafe simulation returns its actual physical peak and positive violated-constraint residuals. An invalid hydraulic/numerical evaluation returns the pre-existing finite failure sentinel with positive integration-status residual; it has no physical peak or complete trajectory. The sentinel is never reported as a physical river level. Programming/configuration errors propagate rather than being converted into design outcomes.', '',
        '## Actual dynamic checks', '',
        f"The four corners `[0,0]`, `[0,1]`, `[1,0]`, `[1,1]` were evaluated for all four category pairs. There were **{len(rows)} boundary integrations** plus **2 exact reproductions** of the saved reference and best-design trajectories, for **{result['actual_ode_attempts']} official ODE attempts**. No OPT search was performed.", '',
        f"Boundary outcomes: **{len(valid)} valid complete trajectories**, **{len(feasible)} feasible**, **{len(valid)-len(feasible)} valid but constraint-infeasible**, and **{len(rows)-len(valid)} invalid evaluations**.", '',
        f"The prescribed 120/60 s comparisons use the same event, forcing knots, initial states, geometry, coefficients and gates. The existing numerical peak-grid tolerance is {plan['peak_grid_tolerance_m']:g} m; it is a diagnostic threshold, not a safety allowance. A 30 s check is triggered only by changed validity/feasibility or a larger peak change. Maximum comparable 120-to-60 s peak change: {max(comparable)*1000:.6f} mm." if comparable else 'No two complete coarse/fine trajectories were available for peak comparison.', '',
        'The full [boundary table](../tables/analysis_checks/boundary_checks.md) reports every case, constraint failure and balance residual. Passing these finite checks does not establish stability for every interior point, continuous-time constraint satisfaction, calibration or a global optimum.', '',
        'Both reference and best-design calls exactly reproduce the saved time, state, input, output and balance arrays, as well as the prior conditioned objective and all 18 residuals. The checks verify reuse of the existing hydraulic model; they do not introduce a new baseline or new optimized result.', '',
        '## Reproduction and evidence', '',
        'Run from the project root in the existing environment:', '', '```sh',
        'python3 -B -m B.wetland_analysis_checks --run --report',
        '# Once results exist, regenerate figures and documentation only:',
        'python3 -B -m B.wetland_analysis_checks --report', '```', '',
        'The driver refuses to overwrite existing checks. Numerical results, complete valid trajectories, frozen configuration, actual source hashes, diagnostics and preservation checks are local under `B/results/analysis_checks/`. The retained `executed_check_driver.py` matches the driver hash recorded in the run plan; subsequent report wording is separate from numerical execution. Earlier results and scientific inputs are checked byte-for-byte and left intact. A and all previous optimized results remain unchanged.', '',
        'Coverage 0 means an unvegetated wetland compartment; 1 means complete coverage by the selected synthetic class. Neither endpoint removes wetland storage or routing. The [source and assumption document](../data_sources_and_assumptions.md) distinguishes the assumed resistance/ET endpoints and their interpolation from observed weather, sourced hydraulic concepts and numerical settings. These tests add no observational support for the assumed coefficients.', '',
        '[Final Stage B guide](../README.md) · [Optimization results](12_wetland_only_optimization.md)', '']
    figures = plot_checks(result, plan)
    lines += ['## Boundary verification figures', '',
        'All three figures use elapsed simulation time on the horizontal axis and read only saved trajectories. River-stage panels show the four coverage corners for each synthetic class pair. Upper safety panels show, for each case, the smallest headroom to the four upstream level limits at each time; the lower panel shows the shared G5 histories. Numerical panels show maximum discrepancies across the cases: river-stage differences at common time nodes and cumulative water-balance residuals. Matched-time differences are distinct from differences between independently sampled peak values. No lines interpolate across coverage or class choices. The complete audit and exact peaks remain in the table and saved records.', '']
    titles = ['River-stage time histories and peak markers', 'Safety headroom and G5 through time', 'Matched-time numerical differences and water-balance histories']
    for title, path in zip(titles, figures):
        lines += [f'![{title}](../plots/analysis_checks/{path.name})', '']
    REPORT.write_text('\n'.join(lines))
    print('Wrote', TABLE, REPORT, 'and', len(figures), 'boundary verification figures.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--report', action='store_true')
    args = parser.parse_args()
    if args.run:
        run_checks()
    if args.report:
        report()
    if not (args.run or args.report):
        parser.print_help()


if __name__ == '__main__':
    main()
