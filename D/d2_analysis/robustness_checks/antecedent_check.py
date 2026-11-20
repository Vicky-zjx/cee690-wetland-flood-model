"""Inspect cached antecedent-weather coverage for a proposed experiment only.

Does not generate forcing, change initial conditions, fill gaps or run an ODE.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / 'C/data/processed/weather_hourly.csv'
MANIFEST = HERE.parent / 'historical_storm_robustness/events/manifest.json'
VARIABLES = ['precipitation_mm_per_h', 'temperature_C', 'solar_W_per_m2', 'wind_m_per_s']
HORIZONS = [24, 48, 72, 168, 336]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads(MANIFEST.read_text())
    if sha(SOURCE) != manifest['standardized_source']['sha256']:
        raise ValueError('Original standardized weather hash changed')
    weather = pd.read_csv(SOURCE, float_precision='round_trip')
    times = pd.to_datetime(weather.timestamp_utc, utc=True)
    if not np.all(np.diff(times.astype('int64')) == 3600 * 10**9):
        raise ValueError('Unexpected source time grid')
    finite = np.isfinite(weather[VARIABLES].to_numpy()).all(axis=1)
    covered = np.isclose(weather.precipitation_coverage_fraction, 1., rtol=0, atol=1e-6)
    valid = finite & covered & (weather.temperature_C.to_numpy() > 0)
    rows, continuous = [], []
    for event in manifest['events']:
        start = pd.Timestamp(event['start_utc'])
        index = np.flatnonzero(times == start)
        if len(index) != 1:
            raise ValueError('Event start missing from cached source')
        stop = int(index[0]); first = stop
        while first > 0 and valid[first - 1]:
            first -= 1
        reasons = []
        if first > 0:
            bad = weather.iloc[first - 1]
            reasons.extend(v + ' missing/nonfinite' for v in VARIABLES if not np.isfinite(bad[v]))
            if not covered[first - 1]:
                reasons.append('incomplete native rainfall coverage')
            if np.isfinite(bad.temperature_C) and bad.temperature_C <= 0:
                reasons.append('outside nonfreezing domain')
        else:
            reasons.append('start of cached calendar')
        continuous.append(dict(event_id=event['event_id'], event_start_utc=str(start),
            contiguous_supported_antecedent_hours=stop-first,
            contiguous_start_utc=str(times.iloc[first]) if first < stop else str(start),
            previous_unusable_hour_utc=str(times.iloc[first-1]) if first > 0 else None,
            limiting_reason='; '.join(reasons)))
        for hours in HORIZONS:
            wanted = start - pd.Timedelta(hours=hours)
            mask = (times >= wanted) & (times < start)
            part = weather.loc[mask]
            missing = {v:int((~np.isfinite(part[v])).sum()) for v in VARIABLES}
            okay = bool(len(part) == hours and valid[mask].all())
            rows.append(dict(event_id=event['event_id'], antecedent_hours=hours,
                start_utc=str(wanted), end_utc=str(start), actual_rows=len(part),
                **{'missing_'+k:v for k,v in missing.items()},
                incomplete_rain_coverage_hours=int((~covered[mask]).sum()),
                nonpositive_temperature_hours=int((part.temperature_C <= 0).sum()),
                minimum_temperature_C=float(part.temperature_C.min()),
                quantified_known_rain_mm=float(part.precipitation_mm_per_h.sum()),
                complete_nonfreezing_supported_window=okay,
                scope='Coverage check only; not a selected or simulated warm-up. Missing rainfall is not zero.'))
    (HERE/'tables').mkdir(exist_ok=True)
    pd.DataFrame(rows).to_csv(HERE/'tables/antecedent_windows.csv', index=False)
    pd.DataFrame(continuous).to_csv(HERE/'tables/antecedent_continuous_coverage.csv', index=False)
    metadata = dict(source=str(SOURCE.relative_to(ROOT)), source_sha256=sha(SOURCE),
        manifest_sha256=sha(MANIFEST), script_sha256=sha(Path(__file__)),
        window_hours=HORIZONS, model_simulations=0, initial_states_changed=False,
        rule='Same finite joint weather, full native-rain coverage and T>0 checks as the original selection; no interpolation/filling added.',
        limitation='Coverage does not demonstrate dynamic equilibration or independence from starting storage. No antecedent experiment has been approved or executed.')
    (HERE/'tables/antecedent_audit.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(pd.DataFrame(continuous).to_string(index=False))


if __name__ == '__main__':
    main()
