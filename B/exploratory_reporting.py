"""Report actual exploratory B.2 trajectories without substituting OPT penalties."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import re

import numpy as np
from A import wetland_flood_model as a2

ROOT = Path(__file__).resolve().parent
CASE_NOTE = 'Hypothetical engineered installation; assumed elevation datum; synthetic A/B response classes.'


def _serial(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def _table(name, rows, fields=None):
    path = ROOT / 'tables' / name
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return str(path.relative_to(ROOT))


def save_run(run, label):
    """Save an existing complete simulation; never construct an invalid trajectory."""
    slug = re.sub(r'[^A-Za-z0-9_-]+', '_', label).strip('_')
    if not slug:
        raise ValueError('A nonempty run label is required.')
    path = ROOT / 'results' / (slug + '.npz')
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {key: run[key] for key in ('time', 'x', 'dx', 'u', 'y', 'theta', 'residuals', 'balance') if key in run}
    metadata = dict(label=label, metrics=run['metrics'], audit=run['audit'], pair=run['pair'], interpretation=CASE_NOTE)
    arrays['metadata_json'] = np.asarray(json.dumps(metadata, default=_serial, allow_nan=False))
    np.savez_compressed(path, **arrays)
    return path


def export_case(reference, candidate, config, records, reference_label, candidate_label):
    """Export computed simulation results; candidate=None explicitly means none found."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plotdir = ROOT / 'plots'
    plotdir.mkdir(parents=True, exist_ok=True)
    runs = [(reference, reference_label, '#2166ac')]
    if candidate is not None:
        runs.append((candidate, candidate_label, '#b35806'))
    figures = []
    tables = []

    def finish(fig, name):
        fig.text(.5, .01, CASE_NOTE, ha='center', fontsize=8, color='#444444')
        fig.tight_layout(rect=(0, .04, 1, 1))
        path = plotdir / name
        fig.savefig(path, dpi=160)
        plt.close(fig)
        figures.append(str(path.relative_to(ROOT)))

    fig, ax = plt.subplots(figsize=(10, 4.8))
    for run, label, color in runs:
        hours = np.asarray(run['time']) / 3600
        stage = run['y'][a2.Y['H_r']]
        metric = run['metrics']
        ax.plot(hours, stage, color=color, label=f"{label}; peak {metric['peak_stage_m']:.6f} m")
        ax.scatter([metric['peak_time_s'] / 3600], [metric['peak_stage_m']], color=color, s=24, zorder=3)
    if candidate is None:
        title = 'Actual reference trajectory; no feasible optimized solution found'
    else:
        reduction = reference['metrics']['peak_stage_m'] - candidate['metrics']['peak_stage_m']
        title = f'Actual physical river peak comparison; reduction {reduction:.6f} m'
    ax.set(xlabel='Elapsed time from event start (h)', ylabel='River level (m, assumed datum)', title=title)
    ax.legend(fontsize=8)
    ax.grid(alpha=.2)
    finish(fig, 'B2_1_river_stage.png')

    fig, axes = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True)
    for i, ax in enumerate(axes):
        offset = 0 if i == 0 else len(config['controls']['gate1_hours'])
        knots = config['controls'][f'gate{i+1}_hours']
        for run, label, color in runs:
            ax.plot(np.asarray(run['time']) / 3600, run['u'][i], color=color, label=label)
            ax.scatter(knots, run['theta'][offset:offset+len(knots)], color=color, s=15)
        shown_max=max(float(np.max(run['u'][i])) for run,_,_ in runs)
        ax.set(ylabel=f'Gate {i+1} opening fraction', ylim=(0, max(.08,1.15*shown_max)))
        ax.axhline(0, color='gray', lw=.6)
        ax.legend(fontsize=8)
        ax.grid(alpha=.2)
    axes[0].set_title('Actual prescribed schedules; zoomed y axes (physical opening range 0–1)')
    axes[-1].set_xlabel('Elapsed time (h)')
    finish(fig, 'B2_2_controls.png')

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    limits = config['safety']['dam_m'] + config['safety']['wetland_m']
    names = ['H_d1', 'H_d2', 'H_w1', 'H_w2']
    titles = ['Dam pool 1', 'Dam pool 2', 'Contained wetland 1', 'Contained wetland 2']
    for name, title, limit, ax in zip(names, titles, limits, axes.flat):
        for run, label, color in runs:
            ax.plot(np.asarray(run['time']) / 3600, run['y'][a2.Y[name]], color=color, label=label)
        ax.axhline(limit, color='#a50026', ls='--', label=f'Fixed hypothetical limit {limit:.2f} m')
        ax.set(title=title, xlabel='Elapsed time (h)', ylabel='Water level (m, assumed datum)')
        ax.legend(fontsize=7)
        ax.grid(alpha=.2)
    finish(fig, 'B2_3_upstream_safety.png')

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    if records:
        pairs = list(dict.fromkeys(row.get('category_pair', 'unspecified') for row in records))
        for pair in pairs:
            indexed = [(j + 1, row) for j, row in enumerate(records) if row.get('category_pair', 'unspecified') == pair]
            indices = [row.get('global_evaluation', j) for j, row in indexed]
            physical = [np.nan if row.get('raw_peak_stage_m') is None else row['raw_peak_stage_m'] for _, row in indexed]
            line, = axes[0].plot(indices, physical, '.', ms=3, label=pair)
            best = [np.nan if row.get('best_feasible_peak_m') is None else row['best_feasible_peak_m'] for _, row in indexed]
            axes[0].plot(indices, best, color=line.get_color(), lw=1)
            violations = [np.nan if row.get('raw_peak_stage_m') is None else row.get('max_scaled_violation', np.nan) for _, row in indexed]
            axes[1].plot(indices, violations, '.', ms=3, color=line.get_color())
            statuses = [2 if row.get('raw_peak_stage_m') is None else (0 if row.get('feasible') else 1) for _, row in indexed]
            axes[2].plot(indices, statuses, '.', ms=3, color=line.get_color())
        axes[0].legend(title='Synthetic class pair; line = best feasible', fontsize=8, title_fontsize=8)
    else:
        axes[0].text(.5, .5, 'No optimizer callbacks recorded', transform=axes[0].transAxes, ha='center')
    axes[0].set(title='Actual OPT callbacks; invalid candidates have no physical peak', ylabel='Physical peak (m)')
    axes[0].ticklabel_format(axis='y',style='plain',useOffset=False)
    axes[1].axhline(config['numerics']['feasibility_tol'], color='#a50026', ls='--', label='Reported feasibility tolerance')
    axes[1].set(ylabel='Max scaled constraint violation\n(valid simulations only)')
    axes[1].set_yscale('symlog', linthresh=1e-6)
    axes[1].legend(fontsize=8)
    axes[2].set(yticks=[0, 1, 2], yticklabels=['Valid, feasible', 'Valid, infeasible', 'Invalid / rejected'], ylim=(-.3, 2.3),
                xlabel='Global callback count (including cached evaluations)')
    for ax in axes:
        ax.grid(alpha=.2)
    finish(fig, 'B2_4_optimization_history.png')

    gate_rows = []
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for i, ax in enumerate(axes):
        for run, label, color in runs:
            time = np.asarray(run['time'])
            depth = run['y'][a2.Y[f'H_d{i+1}']] - config['gates']['gate_sill_z'][i]
            opening = config['gates']['gate_max_opening'][i] * run['u'][i]
            forward = run['y'][a2.Y[f'H_d{i+1}']] > run['y'][a2.Y[f'H_w{i+1}']]
            active = (run['u'][i] > 0) & (depth > 0) & forward
            physical = np.where(active, 1.25 * opening - depth, 0.)
            ax.plot(time / 3600, depth, color=color, label=f'{label}: upstream depth H')
            ax.plot(time / 3600, 1.25 * opening, color=color, ls='--', label=f'{label}: 1.25 a')
            for j in range(len(time)):
                gate_rows.append(dict(run=label, branch=i+1, time_s=float(time[j]), upstream_depth_m=float(depth[j]),
                    opening_height_m=float(opening[j]), required_depth_m=float(1.25*opening[j]), active=bool(active[j]),
                    physical_residual_m=float(physical[j]), strict_pass=bool(physical[j] <= 0)))
        ax.set(ylabel=f'Gate {i+1}: depth / height (m)')
        ax.legend(fontsize=7)
        ax.grid(alpha=.2)
    axes[0].set_title('Additional hydraulic-applicability constraint: H ≥ 1.25 a when open, wet and forward')
    axes[-1].set_xlabel('Elapsed time (h); plotted samples supplement the RHS-stage checks')
    finish(fig, 'B2_5_gate_applicability.png')
    tables.append(_table('gate_applicability_timeseries.csv', gate_rows))

    metrics = [
        ('peak_stage_m', 'm (assumed datum)', 'Physical objective max H_r; no optimizer penalty'),
        ('peak_time_s', 's', 'Sampled physical peak time'),
        ('peak_outlet_m3s', 'm3/s', 'Maximum physical outlet discharge'),
        ('outlet_volume_m3', 'm3', 'Time-integrated outlet discharge'),
        ('balance_max_m3', 'm3', 'Maximum absolute numerical water-balance residual'),
        ('balance_relative', '1', 'Relative numerical water-balance residual'),
        ('peak_near_horizon', 'boolean', 'Peak occurs within final 5% of horizon'),
        ('river_rising_at_end', 'boolean', 'Positive final sampled river-level change'),
        ('flood_duration_h', 'h', 'Blank when no flood threshold was supplied'),
    ]
    performance = []
    for key, units, explanation in metrics:
        performance.append(dict(metric=key, units=units, reference_label=reference_label, reference_value=reference['metrics'].get(key),
            candidate_label=candidate_label if candidate is not None else 'No feasible optimized solution found',
            candidate_value=None if candidate is None else candidate['metrics'].get(key), interpretation=explanation))
    for j, state in enumerate(a2.STATE_NAMES):
        performance.append(dict(metric='end_storage_' + state, units='m3', reference_label=reference_label,
            reference_value=reference['metrics']['end_storage_m3'][j], candidate_label=candidate_label if candidate is not None else 'No feasible optimized solution found',
            candidate_value=None if candidate is None else candidate['metrics']['end_storage_m3'][j], interpretation='Physical end-of-horizon storage'))
    performance.append(dict(metric='absolute_peak_reduction', units='m', reference_label=reference_label, reference_value=0.,
        candidate_label=candidate_label if candidate is not None else 'No feasible optimized solution found',
        candidate_value=None if candidate is None else reference['metrics']['peak_stage_m']-candidate['metrics']['peak_stage_m'],
        interpretation='Reference physical peak minus candidate physical peak; blank if no candidate'))
    tables.append(_table('performance.csv', performance))

    audit_rows = []
    for run, label, _ in runs:
        for row in run['audit']:
            audit_rows.append(dict(run=label, physical_residual=-float(row['margin']), **row))
    tables.append(_table('constraint_audit.csv', audit_rows))

    decisions = []
    for i in range(2):
        decisions.append(dict(decision=f'tau_{i+1}', time_h='', reference_label=reference_label, reference_value=reference['pair'][i],
            candidate_label=candidate_label if candidate is not None else 'None found', candidate_value=None if candidate is None else candidate['pair'][i],
            interpretation='Synthetic hydraulic/ET response class; not a species'))
        decisions.append(dict(decision=f'rho_{i+1}', time_h='', reference_label=reference_label, reference_value=float(reference['theta'][-2+i]),
            candidate_label=candidate_label if candidate is not None else 'None found', candidate_value=None if candidate is None else float(candidate['theta'][-2+i]),
            interpretation='Fixed-design vegetation cover fraction'))
        offset = 0 if i == 0 else len(config['controls']['gate1_hours'])
        for j, hour in enumerate(config['controls'][f'gate{i+1}_hours']):
            decisions.append(dict(decision=f'g{i+1}', time_h=hour, reference_label=reference_label, reference_value=float(reference['theta'][offset+j]),
                candidate_label=candidate_label if candidate is not None else 'None found', candidate_value=None if candidate is None else float(candidate['theta'][offset+j]),
                interpretation='Prescribed gate opening fraction; linear interpolation between knots'))
    tables.append(_table('decision_values.csv', decisions))
    tables.append(_table('optimization_evaluations.csv', records,
        None if records else ['global_evaluation', 'category_pair', 'raw_peak_stage_m', 'conditioned_objective', 'max_scaled_violation', 'feasible', 'failure']))
    rejected=[]
    for row in records:
        if row.get('raw_peak_stage_m') is not None:
            continue
        error=row.get('failure','');diagnostic={}
        if '{' in error:
            try:diagnostic=json.loads(error[error.index('{'):])
            except json.JSONDecodeError:pass
        rejected.append(dict(global_evaluation=row.get('global_evaluation'),category_pair=row.get('category_pair'),
            branch=diagnostic.get('branch'),time_s=diagnostic.get('time_s'),
            water_depth_m=diagnostic.get('water_depth_m'),opening_height_m=diagnostic.get('opening_height_m'),
            residual_m=diagnostic.get('residual_m'),regime=diagnostic.get('regime'),
            physical_peak_m=None,interpretation='Rejected incomplete trajectory; no physical peak exists',failure=error))
    tables.append(_table('rejected_candidates.csv',rejected,None if rejected else
        ['global_evaluation','category_pair','branch','time_s','water_depth_m','opening_height_m','residual_m','regime','physical_peak_m','interpretation','failure']))
    return dict(figures=figures, tables=tables, candidate_available=candidate is not None)
