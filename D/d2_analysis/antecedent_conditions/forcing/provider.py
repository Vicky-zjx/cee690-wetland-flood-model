"""Separate hourly mapping across antecedent/event seam; no simulation here."""
from __future__ import annotations
import copy
import hashlib
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import numpy as np
from D.g5_diagnostics import runner as historical


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def resolve(path):
    path = Path(path)
    return path if path.is_absolute() else ROOT / path


def record_for(event_id, hours):
    matches = [r for r in json.loads((HERE/'manifest.json').read_text())['cases']
               if r['event_id'] == event_id and r['antecedent_hours'] == hours]
    if len(matches) != 1:
        raise ValueError('Expected one prepared antecedent forcing record')
    return copy.deepcopy(matches[0])


def read_forcing(path):
    """Original preflight plus exact event-tail identity and attached seam record."""
    path = resolve(path)
    meta = json.loads(path.with_suffix('.meta.json').read_text())
    if sha(path) != meta['forcing_sha256']:
        raise ValueError('Combined weather bytes changed after preparation')
    original = resolve(meta['original_event_forcing_file'])
    if sha(original) != meta['original_event_forcing_sha256']:
        raise ValueError('Original historical event weather changed')
    combined_lines = path.read_bytes().splitlines(keepends=True)
    hours = int(meta['antecedent_hours'])
    event_bytes = combined_lines[0] + b''.join(combined_lines[hours+1:])
    if event_bytes != original.read_bytes():
        raise ValueError('Combined input does not retain original event rows byte-for-byte')
    frame, diagnostics = historical.adapter().read_forcing(path)
    if len(frame) != hours + 144:
        raise ValueError('Require requested prefix and original complete 144-hour event')
    frame.attrs['antecedent_forcing'] = meta
    diagnostics.update(antecedent_hours=hours, event_offset_s=hours*3600,
                       original_event_rows_byte_exact=True,
                       mapping_seam='Separate unchanged hourly mappings; original event endpoint holds preserved')
    return frame, diagnostics


def _partition_fixed(fixed, hours):
    local = dict(fixed)
    local['model_config'] = copy.deepcopy(fixed['model_config'])
    local['model_config']['numerics']['duration_s'] = hours * 3600.
    local['manifest'] = copy.deepcopy(fixed['manifest'])
    local['manifest']['duration_s'] = hours * 3600.
    return local


def hourly_input_segments(frame, fixed):
    """Prefix and original event map independently; shift event clock, not values.

    The frozen feedback driver replaces rows0/1 with continuous controller
    outputs. Rows2 onward are weather and fixed boundaries. Their event values
    retain the original event's first/last half-hour temperature/wind holds.
    """
    meta = frame.attrs.get('antecedent_forcing')
    if not meta:
        raise ValueError('Use this provider read_forcing so the seam is explicit')
    hours = int(meta['antecedent_hours'])
    duration = float(fixed['model_config']['numerics']['duration_s'])
    if duration != (hours+144)*3600 or len(frame) != hours+144:
        raise ValueError('Extended duration and saved hourly weather must agree')
    original = historical.adapter()
    for offset, section, length in [(0., frame.iloc[:hours], hours),
                                     (hours*3600., frame.iloc[hours:], 144)]:
        section = section.reset_index(drop=True)
        for time, inputs in original.hourly_input_segments(section, _partition_fixed(fixed, length)):
            yield time + offset, inputs


def forcing_conservation(frame, fixed):
    """Conserve segment rain/solar integrals; assert original event inputs exact."""
    meta = frame.attrs['antecedent_forcing']; hours = int(meta['antecedent_hours'])
    segments = list(hourly_input_segments(frame, fixed))
    original_frame, _ = historical.adapter().read_forcing(resolve(meta['original_event_forcing_file']))
    reference = list(historical.adapter().hourly_input_segments(original_frame, _partition_fixed(fixed, 144)))
    for (time, inputs), (old_time, old_inputs) in zip(segments[hours:], reference):
        if not np.array_equal(time-hours*3600, old_time) or not np.array_equal(inputs[2:], old_inputs[2:]):
            raise ValueError('Event weather/boundary inputs differ from unchanged original event mapping')
    rain = sum(float(np.trapezoid(u[2], t)) for t, u in segments)
    solar = sum(float(np.trapezoid(u[4], t)/3600) for t, u in segments)
    source_rain = float(frame.precipitation_mm_per_h.sum()/1000)
    source_solar = float(frame.solar_W_per_m2.sum())
    prefix, event = segments[hours-1][1][:, -1], segments[hours][1][:, 0]
    checks = dict(
        source_precipitation_m=source_rain, integrated_precipitation_m=rain,
        precipitation_error_m=rain-source_rain,
        source_solar_Wh_m2=source_solar, integrated_solar_Wh_m2=solar,
        solar_error_Wh_m2=solar-source_solar,
        event_weather_and_boundary_inputs_bit_exact=True, event_output_knots_shift_exact=True,
        event_offset_s=hours*3600, antecedent_hours=hours,
        seam_temperature_jump_C=float(event[3]-prefix[3]),
        seam_wind_jump_m_s=float(event[5]-prefix[5]),
        seam_precipitation_rate_jump_m_s=float(event[2]-prefix[2]),
        seam_solar_mean_jump_W_m2=float(event[4]-prefix[4]),
        mapping='Unchanged mapper applied independently to antecedent and original144h event; time shifted only; controller/state remain continuous in caller',
        pass_conservation=bool(abs(rain-source_rain)<=1e-12 and abs(solar-source_solar)<=1e-7))
    return checks
