"""Audit saved full-window peak differences without running either ODE model."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
HIST = ROOT / 'D/d2_analysis/historical_storm_robustness'
PARAM = ROOT / 'D/d2_analysis/parameter_uncertainty'
TABLES = HERE / 'tables'
PLOTS = HERE / 'plots'
LABELS = ['Sep 2022 benchmark', 'Aug 2023', 'Jul 2022', 'Aug 2019', 'Aug 2024']
COLORS = ['#6f4c9b', '#008c95', '#d98324', '#3569a4', '#b94164']
COMPARE_TOL_M = 1e-12
DIFFERENT_TOL_M = 1e-9


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def npz_path(event, sample, case):
    if event['benchmark']:
        path = PARAM / 'results/pilot' / f'pilot_{sample}__{case}__15s.npz'
    else:
        path = HIST / 'results/ensemble' / event['event_id'] / f'{sample}__{case}__15s.npz'
    return path if path.is_file() else None


def read_array(path, benchmark):
    with np.load(path, allow_pickle=False) as z:
        result = {k: z[k].reshape(-1) for k in ['time', 'downstream_stage', 'downstream_discharge']}
        result['opening'] = z['u'][:2] * 0.3 if benchmark else z['actual_opening_m']
    assert result['time'][0] == 0 and result['time'][-1] == 518400
    assert np.all(np.diff(result['time']) == 15)
    assert all(np.isfinite(v).all() for v in result.values())
    return result


def secondary_start(event):
    start = datetime.fromisoformat(event['largest_quantified_wet_spell']['start_utc'])
    origin = datetime.fromisoformat(event['start_utc'])
    return (start - origin).total_seconds()


def audit(summary, events):
    rows, examples = [], {}
    event_map = {e['event_id']: e for e in events}
    for source in summary.to_dict('records'):
        event = event_map[source['event_id']]
        row = {key: source[key] for key in ['event_id', 'sample_id', 'draw_index', 'm_K_u', 'm_routing', 'm_R_r', 'R_r', 'paired_valid', 'wetland_classification', 'no_wetland_classification', 'wetland_peak_stage_m', 'no_wetland_peak_stage_m', 'wetland_initial_stage_m', 'no_wetland_initial_stage_m', 'wetland_peak_stage_time_s', 'no_wetland_peak_stage_time_s', 'wetland_peak_discharge_m3_s', 'no_wetland_peak_discharge_m3_s', 'delta_peak_stage_m', 'delta_peak_discharge_m3_s', 'maximum_paired_gate_opening_difference_m', 'paired_latch_mismatch_count']}
        assert row['paired_valid']
        row['primary_exact_zero'] = row['delta_peak_stage_m'] == 0
        row['both_primary_peaks_initial'] = row['wetland_peak_stage_time_s'] == 0 and row['no_wetland_peak_stage_time_s'] == 0
        row['zero_iff_both_initial_from_saved_summary'] = row['primary_exact_zero'] == row['both_primary_peaks_initial']
        assert row['zero_iff_both_initial_from_saved_summary']
        row['both_primary_peak_stages_equal_initial'] = row['wetland_peak_stage_m'] == row['wetland_initial_stage_m'] and row['no_wetland_peak_stage_m'] == row['no_wetland_initial_stage_m']
        row['peak_discharge_exact_zero'] = row['delta_peak_discharge_m3_s'] == 0
        assert row['primary_exact_zero'] == row['peak_discharge_exact_zero']
        row['primary_peak_delay_h_when_both_after_initial'] = (row['wetland_peak_stage_time_s'] - row['no_wetland_peak_stage_time_s']) / 3600 if row['wetland_peak_stage_time_s'] > 0 and row['no_wetland_peak_stage_time_s'] > 0 else np.nan
        paths = {case: npz_path(event, row['sample_id'], case) for case in ('wetland', 'no_wetland')}
        row['full_trajectory_pair_available'] = all(p is not None for p in paths.values())
        for case in paths:
            jp = HIST / 'results/ensemble' / event['event_id'] / f"{row['sample_id']}__{case}__15s.json"
            record = json.loads(jp.read_text())
            assert record['full_valid_trajectory_available']
            assert record['classification'] == row[f'{case}_classification']
            assert abs(record['physical_river_peak_m'] - row[f'{case}_peak_stage_m']) <= COMPARE_TOL_M
            assert abs(record['physical_river_peak_discharge_m3_s'] - row[f'{case}_peak_discharge_m3_s']) <= COMPARE_TOL_M
            assert record['physical_river_peak_time_s'] == row[f'{case}_peak_stage_time_s']
            row[f'{case}_peak_discharge_time_s'] = record['physical_river_peak_discharge_time_s']
            row[f'{case}_diagnostic_file'] = str(jp.relative_to(ROOT))
            row[f'{case}_diagnostic_sha256'] = digest(jp)
            row[f'{case}_array_file'] = str(paths[case].relative_to(ROOT)) if paths[case] else None
            row[f'{case}_array_sha256'] = digest(paths[case]) if paths[case] else None
        if row['full_trajectory_pair_available']:
            w, n = (read_array(paths[c], event['benchmark']) for c in ('wetland', 'no_wetland'))
            t, wh, nh, wq, nq = w['time'], w['downstream_stage'], n['downstream_stage'], w['downstream_discharge'], n['downstream_discharge']
            assert np.array_equal(t, n['time'])
            for case, stage, discharge in [('wetland', wh, wq), ('no_wetland', nh, nq)]:
                assert abs(float(stage.max()) - row[f'{case}_peak_stage_m']) <= COMPARE_TOL_M
                assert abs(float(discharge.max()) - row[f'{case}_peak_discharge_m3_s']) <= COMPARE_TOL_M
                assert float(t[stage.argmax()]) == row[f'{case}_peak_stage_time_s']
                assert float(t[discharge.argmax()]) == row[f'{case}_peak_discharge_time_s']
            row['zero_iff_both_initial_from_arrays'] = bool((nh.max() - wh.max() == 0) == (wh.argmax() == 0 and nh.argmax() == 0))
            assert row['zero_iff_both_initial_from_arrays']
            row['waveform_stage_arrays_exactly_identical'] = bool(np.array_equal(wh, nh))
            row['waveform_max_abs_stage_difference_m'] = float(np.max(abs(wh - nh)))
            row['waveform_max_abs_discharge_difference_m3_s'] = float(np.max(abs(wq - nq)))
            row['waveform_stage_difference_exceeds_1e_minus_9_m'] = row['waveform_max_abs_stage_difference_m'] > DIFFERENT_TOL_M
            row['actual_gate_arrays_exactly_identical'] = bool(np.array_equal(w['opening'], n['opening']))
            row['actual_gate_max_abs_difference_from_arrays_m'] = float(np.max(abs(w['opening'] - n['opening'])))
            assert abs(row['actual_gate_max_abs_difference_from_arrays_m'] - row['maximum_paired_gate_opening_difference_m']) <= COMPARE_TOL_M
            row['maximum_wet_minus_no_stage_m'] = float(np.max(wh - nh))
            row['maximum_no_minus_wet_stage_m'] = float(np.max(nh - wh))
            i = np.flatnonzero(t >= secondary_start(event))
            iw, ino = int(i[wh[i].argmax()]), int(i[nh[i].argmax()])
            row['secondary_start_s'] = secondary_start(event)
            row['secondary_wetland_peak_stage_m'] = float(wh[iw])
            row['secondary_no_wetland_peak_stage_m'] = float(nh[ino])
            row['secondary_delta_peak_stage_m'] = float(nh[ino] - wh[iw])
            row['secondary_wetland_peak_time_h'] = float(t[iw] / 3600)
            row['secondary_no_wetland_peak_time_h'] = float(t[ino] / 3600)
            row['secondary_peak_delay_h'] = float((t[iw] - t[ino]) / 3600)
            row['secondary_wetland_peak_is_left_boundary'] = bool(iw == i[0])
            row['secondary_no_wetland_peak_is_left_boundary'] = bool(ino == i[0])
            after = np.flatnonzero(t >= max(t[iw], t[ino]))
            later = after[wh[after] - nh[after] > 1e-6]
            row['later_wet_above_no_by_1e_minus_6_m'] = bool(later.size)
            row['later_first_wet_above_no_time_h'] = float(t[later[0]] / 3600) if later.size else np.nan
            row['later_maximum_wet_minus_no_m'] = float(np.max((wh - nh)[after]))
            key = (event['event_id'], row['sample_id'])
            if key in EXAMPLE_KEYS:
                examples[key] = (w, n)
        rows.append(row)
    return pd.DataFrame(rows), examples


def summarize_samples(data, events):
    rows = []
    for e in events:
        group = data[data.event_id == e['event_id']]
        a = group[group.full_trajectory_pair_available]
        z = a[a.primary_exact_zero]
        row = dict(event_id=e['event_id'], samples=len(group), array_pairs=len(a),
                   primary_positive=int((group.delta_peak_stage_m > 0).sum()), primary_zero=int(group.primary_exact_zero.sum()), primary_negative=int((group.delta_peak_stage_m < 0).sum()),
                   wetland_initial_peaks=int((group.wetland_peak_stage_time_s == 0).sum()), no_wetland_initial_peaks=int((group.no_wetland_peak_stage_time_s == 0).sum()),
                   zero_iff_both_initial=bool(group.zero_iff_both_initial_from_saved_summary.all()),
                   zero_with_waveform_arrays=len(z), zero_with_nonidentical_hydrographs=int((z.waveform_stage_difference_exceeds_1e_minus_9_m == True).sum()),
                   exactly_identical_array_pairs=int((a.waveform_stage_arrays_exactly_identical == True).sum()),
                   zero_waveform_max_abs_difference_min_m=z.waveform_max_abs_stage_difference_m.min(), zero_waveform_max_abs_difference_max_m=z.waveform_max_abs_stage_difference_m.max(),
                   zero_secondary_peak_difference_min_m=z.secondary_delta_peak_stage_m.min(), zero_secondary_peak_difference_max_m=z.secondary_delta_peak_stage_m.max(),
                   secondary_peak_delay_min_h=a.secondary_peak_delay_h.min(), secondary_peak_delay_median_h=a.secondary_peak_delay_h.median(), secondary_peak_delay_max_h=a.secondary_peak_delay_h.max(),
                   secondary_peak_wet_at_left_boundary=int((a.secondary_wetland_peak_is_left_boundary == True).sum()), secondary_peak_no_at_left_boundary=int((a.secondary_no_wetland_peak_is_left_boundary == True).sum()),
                   waveform_max_gate_difference_m=a.actual_gate_max_abs_difference_from_arrays_m.max(), saved_history_max_gate_difference_m=group.maximum_paired_gate_opening_difference_m.max(),
                   saved_latch_mismatches=int(group.paired_latch_mismatch_count.sum()), later_wet_above_no_pairs=int((a.later_wet_above_no_by_1e_minus_6_m == True).sum()),
                   rainfall_total_mm=e['rainfall_total_mm'], max_hourly_mm=e['max_hourly_mm'], max_24h_mm=e['max_24h_mm'], positive_rain_hours=e['positive_rain_hours'], largest_wet_spell_duration_h=e['largest_quantified_wet_spell']['span_h'])
        rows.append(row)
    return pd.DataFrame(rows)


def primary_peak_summary(data, events):
    """Summarize all 64 saved samples, retaining peaks at the initial time."""
    rows = []
    for event in events:
        group = data[data.event_id == event['event_id']]
        for case in ('wetland', 'no_wetland'):
            row = {
                'event_id': event['event_id'], 'case': case, 'samples': len(group),
                'initial_stage_peak_count': int((group[f'{case}_peak_stage_time_s'] == 0).sum()),
                'initial_discharge_peak_count': int((group[f'{case}_peak_discharge_time_s'] == 0).sum()),
                'timing_includes_initial_peaks': True,
            }
            quantities = {
                'peak_stage_m': group[f'{case}_peak_stage_m'],
                'peak_discharge_m3_s': group[f'{case}_peak_discharge_m3_s'],
                'peak_stage_time_h': group[f'{case}_peak_stage_time_s'] / 3600,
                'peak_discharge_time_h': group[f'{case}_peak_discharge_time_s'] / 3600,
            }
            for quantity, values in quantities.items():
                for statistic in ('mean', 'median', 'min', 'max'):
                    row[f'{quantity}_{statistic}'] = float(getattr(values, statistic)())
            rows.append(row)
    return pd.DataFrame(rows)


def figures(data, stats, events, choices, examples):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 180})
    fig, ax = plt.subplots(2, 2, figsize=(13, 9.7), layout='constrained')
    for e, label, color in zip(events, LABELS, COLORS):
        g = data[data.event_id == e['event_id']]
        ax[0, 0].scatter(g.draw_index, 1000*g.delta_peak_stage_m, s=15, alpha=.7, color=color, label=label)
    ax[0, 0].set(title='A. All 320 full-window outcomes', xlabel='Saved parameter draw index (same 64 per event)', ylabel='Peak stage: no wetland − wetland (mm)')
    ax[0, 0].axhline(0, color='0.4', lw=.8)
    ax[0, 0].legend(fontsize=8, ncols=2)
    bars = ax[0, 1].bar(np.arange(5), stats.primary_zero, color=COLORS, alpha=.85)
    ax[0, 1].set_xticks(np.arange(5), ['Sep22', 'Aug23', 'Jul22', 'Aug19', 'Aug24'])
    ax[0, 1].set(title='B. Exact zero peak difference', ylabel='Count out of 64 draws', ylim=(0, 64))
    for b, count in zip(bars, stats.primary_zero):
        ax[0, 1].text(b.get_x()+b.get_width()/2, count+1, str(count), ha='center')
    ax[0, 1].text(.03, .95, 'Every zero has both maxima at t = 0.\nAll 87 zero-pair hydrographs differ later.', transform=ax[0, 1].transAxes, va='top', fontsize=10)
    a = data[data.full_trajectory_pair_available]
    for iszero, color, label in [(False, '#3569a4', 'Positive primary difference'), (True, '#c24b46', 'Zero primary difference')]:
        g = a[a.primary_exact_zero == iszero]
        ax[1, 0].scatter(1000*g.delta_peak_stage_m, 1000*g.waveform_max_abs_stage_difference_m, s=24, alpha=.65, color=color, label=label)
        ax[1, 1].scatter(1000*g.delta_peak_stage_m, 1000*g.secondary_delta_peak_stage_m, s=24, alpha=.65, color=color)
    ax[1, 0].set(title='C. Waveforms can differ despite equal maxima', xlabel='Primary full-window peak difference (mm)', ylabel='Largest absolute paired stage difference (mm)')
    ax[1, 0].legend(fontsize=9)
    ax[1, 1].set(title='D. Separate descriptive storm-response metric', xlabel='Primary full-window peak difference (mm)', ylabel='Secondary peak difference after wet-spell start (mm)')
    limit = max(a.delta_peak_stage_m.max(), a.secondary_delta_peak_stage_m.max())*1000
    ax[1, 1].plot([0, limit], [0, limit], '--', color='0.4', lw=1, label='Equal primary / secondary')
    ax[1, 1].legend(fontsize=9)
    for axis in ax.flat:
        axis.grid(alpha=.18)
    fig.suptitle('Zero full-window peak attenuation is not an identical hydrograph\nPanels A–B: 320 pairs; C–D: 258 pairs with retained full arrays (256 added-event + 2 benchmark)', fontsize=13)
    fig.savefig(PLOTS/'01_zero_peak_diagnostics.png')
    plt.close(fig)

    fig, axes = plt.subplots(5, 2, figsize=(14, 17), sharex=True, layout='constrained')
    for index, (event, label) in enumerate(zip(events, LABELS)):
        for col, kind in enumerate(['first_zero', 'first_positive']):
            axis = axes[index, col]
            sample = choices[event['event_id']][kind]
            if sample is None:
                axis.text(.5, .5, 'No zero-effect draw exists\n0 / 64 in the benchmark', transform=axis.transAxes, ha='center', va='center', fontsize=12)
                axis.set_title(f'{label} | first zero: unavailable', loc='left', fontsize=11)
                axis.set_axis_off()
                continue
            w, n = examples[(event['event_id'], sample)]
            t = w['time']/3600
            axis.plot(t, n['downstream_stage'], color='#ce7527', lw=1.5, label='No wetland')
            axis.plot(t, w['downstream_stage'], color='#187a8c', lw=1.5, label='Wetland')
            axis.axhline(w['downstream_stage'][0], color='0.35', ls='--', lw=.9, label='Initial stage 97.35 m')
            spell = event['largest_quantified_wet_spell']
            start = secondary_start(event)/3600
            axis.axvspan(start, start+spell['span_h'], color='0.7', alpha=.17, label='Largest quantified wet spell')
            row = data[(data.event_id == event['event_id']) & (data.sample_id == sample)].iloc[0]
            axis.set_title(f"{label} | {kind.replace('_', ' ')}: {sample}\nFull-window ΔH = {1000*row.delta_peak_stage_m:.3f} mm", loc='left', fontsize=11)
            for series, color in [(w['downstream_stage'], '#187a8c'), (n['downstream_stage'], '#ce7527')]:
                i = int(series.argmax()); axis.scatter([t[i]], [series[i]], color=color, s=18, zorder=5)
            axis.set_ylabel('River stage (m)')
            axis.grid(alpha=.16)
            axis.set_xlim(0, 144)
        active = [axes[index, k] for k, kind in enumerate(['first_zero', 'first_positive']) if choices[event['event_id']][kind] is not None]
        limits = [axis.get_ylim() for axis in active]
        lower, upper = min(pair[0] for pair in limits), max(pair[1] for pair in limits)
        for axis in active:
            axis.set_ylim(lower, upper)
            axis.ticklabel_format(axis='y', style='plain', useOffset=False)
    handles, labels = axes[0, 1].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncols=4, frameon=False, fontsize=11)
    for axis in axes[-1]:
        axis.set_xlabel('Time from saved window start (h)')
    fig.suptitle('Predeclared examples: lowest sample ID with zero / positive primary attenuation\nAll 144 hours and time-zero peaks retained; columns may use different parameter vectors', fontsize=14)
    fig.savefig(PLOTS/'02_event_hydrographs.png')
    plt.close(fig)


def notes(data, stats, events, choices):
    zero = data[data.primary_exact_zero]
    available = data[data.full_trajectory_pair_available]
    lines = ['# Saved-trajectory robustness checks', '',
        'This analysis reads existing saved results only. It does not integrate an ODE, redraw parameters, alter weather, remove time zero, or modify any historical model, result or controller. The 320 paired outcomes use the same 64 exploratory parameter vectors under five selected fixed weather windows; they are not 320 independent weather observations.', '',
        '## Coverage and primary result', '',
        'All 320 primary peak values, peak times, peak discharges and classifications were checked against their saved per-model JSON records. Full arrays are retained for all 256 added-event pairs and benchmark draws `mc_00000` and `mc_00001`: 258 paired waveform audits. The other 62 benchmark pairs remain included in primary statistics but have no retained full arrays for this check. No unavailable waveform is reconstructed.', '',
        '| Event | Positive / zero / negative ΔH | Both initial maxima | Zero pairs with different waveforms | Array pairs |',
        '|---|---:|---:|---:|---:|']
    for label, row in zip(LABELS, stats.to_dict('records')):
        lines.append(f"| {label} | {row['primary_positive']} / {row['primary_zero']} / {row['primary_negative']} | {row['no_wetland_initial_peaks']} | {row['zero_with_nonidentical_hydrographs']} / {row['zero_with_waveform_arrays']} | {row['array_pairs']} |")
    lines += ['',
        f"Every exact zero in the saved primary metric is equivalent, within these results, to both first full-window maxima occurring at the shared initial stage 97.35 m. All {len(zero)} zero pairs have full arrays and later stage differences larger than 1e−9 m; none has an identical hydrograph. Their maximum absolute waveform differences range from {1000*zero.waveform_max_abs_stage_difference_m.min():.3f} to {1000*zero.waveform_max_abs_stage_difference_m.max():.3f} mm. These remain valid zero full-window peak-reduction outcomes. Their later dynamical differences are masked by the initial maximum in that metric; secondary differences are not a recovered or replacement primary benefit. This does not mean the initial stored wetland water has been causally isolated or that a different primary objective should replace the accepted one.", '',
        '[Primary peak summaries](tables/trajectory_peak_summary.csv) report mean, median, minimum and maximum stage peak, discharge peak and their first-occurrence times for each event and model, using all 64 saved samples. The ten rows retain time-zero maxima and give their counts; mean or median peak time across such samples is a descriptive full-window statistic, not an isolated storm-delay estimate.', '',
        'Peak-discharge difference is also exactly zero in these same pairs. With the unchanged monotone outlet rating, equal maximum stages imply equal peak discharges; these are not independent confirmations that the hydrographs coincide. For all array-available pairs, the saved stage and discharge maxima and their first-occurrence times match the source records.', '',
        f"Realized gate openings match between models to {data.maximum_paired_gate_opening_difference_m.max():.3g} m in all saved paired history summaries, with {int(data.paired_latch_mismatch_count.sum())} latch disagreements. Direct array comparisons cover the 258 retained pairs and agree with those summaries within 1e−12 m. Every full-array pair has distinct river waveforms; the difference is not a changed gate-operation policy. The wetland model's additional accepted initial water and different storage/routing remain bundled in the comparison, not separately identified causal effects.", '',
        '## Separate descriptive event-response window', '',
        'For timing diagnostics only, a secondary window begins at the existing manifest\'s largest quantified wet-spell start and ends at 144 h. The spell was selected from rainfall, without model outcomes; positive quantified hours are joined unless separated by at least six quantified-zero hours, and trace-zero is not verified dry. Secondary peaks are recorded separately and do not replace primary full-window peaks. The primary peak delay is left unavailable when either full-window peak occurs at t=0.', '',
        '| Event | Zero-primary secondary difference range (mm) | Secondary peak-delay range (h), all available pairs | Later wet > no pairs |',
        '|---|---:|---:|---:|']
    for label, row in zip(LABELS, stats.to_dict('records')):
        low, high = row['zero_secondary_peak_difference_min_m'], row['zero_secondary_peak_difference_max_m']
        interval = 'No zero draws' if pd.isna(low) else f'{1000*low:.3f}–{1000*high:.3f}'
        lines.append(f"| {label} | {interval} | {row['secondary_peak_delay_min_h']:.3f}–{row['secondary_peak_delay_max_h']:.3f} | {row['later_wet_above_no_pairs']} / {row['array_pairs']} |")
    lines += ['',
        'The later comparison begins after both secondary peaks. A later positive difference means wetland stage exceeds no-wetland stage by more than 1e−6 m at a saved node; its first time is resolved to 15 s, not an interpolated crossing. This is a descriptive reporting tolerance, not a safety or applicability relaxation. Smaller primary peaks can coexist with a higher later recession. Some later rain may still be present, so the check does not imply that rainfall has ceased.', '',
        f"Secondary maxima occur at the secondary left boundary in {int(stats.secondary_peak_wet_at_left_boundary.sum())} wetland and {int(stats.secondary_peak_no_at_left_boundary.sum())} no-wetland array trajectories. Any such time is retained and flagged, not interpreted as an interior storm-response peak.", '',
        '## Figures and deterministic example selection', '',
        '[Figure 1](plots/01_zero_peak_diagnostics.png) retains all 320 full-window outcomes and displays exact zero counts; its waveform panels use all 258 retained paired arrays. No zero-effect samples are removed. The secondary comparison is visibly labeled as a different descriptive statistic.', '',
        '[Figure 2](plots/02_event_hydrographs.png) uses the lowest saved sample ID with an exact zero primary difference and the lowest ID with a positive difference in each event. These IDs were written to `tables/trajectory_example_selection.json` before reading arrays for the audit or generating the plots. This rule does not select the largest benefit. The benchmark has no zero pair and its panel remains empty. Both columns retain the entire 144-hour window and can use different parameter vectors, so their contrast is not an intervention on a single parameter. Within each event row the two plotted panels share a stage scale; different event rows may use different stage scales.', '',
        '| Event | First zero ID | First positive ID |', '|---|---|---|']
    for event, label in zip(events, LABELS):
        c = choices[event['event_id']]
        lines.append(f"| {label} | {c['first_zero'] or 'None'} | {c['first_positive'] or 'None'} |")
    lines += ['',
        '## Rainfall and limits', '',
        'Rainfall totals, derived UTC-hour and 24-hour maxima, wet-hour counts and wet-spell duration are copied from the unchanged event manifest into `tables/trajectory_event_summary.csv`. All added events have lower 24-hour depth than the September benchmark. Four weather variables vary jointly; the benchmark uses its original 120 s reconstruction whereas the added events use preserved Stage C hourly inputs. These selected, complete nonfreezing windows do not identify a causal rainfall-intensity relationship or represent the distribution of all regional storms.', '',
        'This audit preserves the initial-condition floor instead of claiming it has been removed. It does not determine how an initially empty wetland, an observed antecedent state, a different event start, or a separate warm-up would behave. Such claims require controlled simulations and are outside this saved-result check. No downstream flood threshold is supplied or invented.', '',
        '## Reproduction and output evidence', '',
        'Run `/opt/miniconda3/bin/python3 -B D/d2_analysis/robustness_checks/trajectory_checks.py` from the project root (the tested existing environment). The system `python3` lacks Matplotlib; no new dependency was installed. Only `trajectory_*` tables, this note and the two assigned figures are written under the new robustness-check directory. `tables/trajectory_sample_audit.csv` contains all 320 rows, source paths and hashes, complete coverage flags, primary/secondary timing, waveform differences and gate checks. `tables/trajectory_validation.json` records the aggregate audit and input hashes.', '']
    (HERE/'trajectory_findings.md').write_text('\n'.join(lines))


def main():
    global EXAMPLE_KEYS
    TABLES.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)
    mp = HIST/'events/manifest.json'
    sp = HIST/'tables/all_event_paired_metrics.csv'
    events = json.loads(mp.read_text())['events']
    summary = pd.read_csv(sp).sort_values(['event_id', 'sample_id'])
    assert len(summary) == 320 and summary.groupby('event_id').size().eq(64).all()
    choices = {}
    for event in events:
        group = summary[summary.event_id == event['event_id']].sort_values('sample_id')
        choices[event['event_id']] = {}
        for key, selection in [('first_zero', group.delta_peak_stage_m.eq(0)), ('first_positive', group.delta_peak_stage_m.gt(0))]:
            ids = group.loc[selection, 'sample_id'].tolist()
            choices[event['event_id']][key] = ids[0] if ids else None
    (TABLES/'trajectory_example_selection.json').write_text(json.dumps({'rule': 'Lowest saved sample ID in each primary exact-zero / positive class; fixed before trajectory audit and plotting; no largest-benefit selection', 'source_summary_sha256': digest(sp), 'choices': choices}, indent=2)+'\n')
    EXAMPLE_KEYS = {(event, sample) for event, groups in choices.items() for sample in groups.values() if sample is not None}
    data, examples = audit(summary, events)
    stats = summarize_samples(data, events)
    data.to_csv(TABLES/'trajectory_sample_audit.csv', index=False)
    stats.to_csv(TABLES/'trajectory_event_summary.csv', index=False)
    peaks = primary_peak_summary(data, events)
    peaks.to_csv(TABLES/'trajectory_peak_summary.csv', index=False)
    figures(data, stats, events, choices, examples)
    notes(data, stats, events, choices)
    result = {
        'status': 'saved_result_checks_passed', 'no_ODE_executed': True,
        'summary_pairs_checked': len(data), 'primary_peak_summary_rows': len(peaks), 'full_array_pairs_checked': int(data.full_trajectory_pair_available.sum()),
        'exact_zero_primary_pairs': int(data.primary_exact_zero.sum()),
        'all_zero_iff_both_initial': bool(data.zero_iff_both_initial_from_saved_summary.all()),
        'zero_pair_array_coverage_complete': bool(data.loc[data.primary_exact_zero, 'full_trajectory_pair_available'].all()),
        'zero_pairs_with_nonidentical_hydrographs': int((data.loc[data.primary_exact_zero, 'waveform_stage_difference_exceeds_1e_minus_9_m'] == True).sum()),
        'summary_source': str(sp.relative_to(ROOT)), 'summary_sha256': digest(sp), 'manifest_sha256': digest(mp), 'analysis_script_sha256': digest(Path(__file__)),
        'primary_definition': 'max over saved0..144h including initial; difference no_wetland minus wetland; zero exact saved floating value0',
        'secondary_definition': 'max after largest quantified wet-spell start through144h; separate descriptive metric, never replaces primary',
        'missing_waveform_pairs': data.loc[~data.full_trajectory_pair_available, ['event_id', 'sample_id']].to_dict('records'),
        'source_peak_comparison_tolerance_m': COMPARE_TOL_M, 'waveform_different_reporting_tolerance_m': DIFFERENT_TOL_M,
        'actual_gate_max_difference_from_arrays_m': float(data.actual_gate_max_abs_difference_from_arrays_m.max()),
        'saved_gate_history_max_difference_m': float(data.maximum_paired_gate_opening_difference_m.max()),
        'saved_latch_mismatch_count': int(data.paired_latch_mismatch_count.sum()),
    }
    (TABLES/'trajectory_validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'missing_waveform_pairs'}, indent=2))
    print(stats.to_string(index=False))


if __name__ == '__main__':
    main()
