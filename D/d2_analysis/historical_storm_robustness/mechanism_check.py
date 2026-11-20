"""Read saved mc_00000 arrays; do not integrate or modify existing results."""
from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SAMPLE_ID = "mc_00000"
CROSSOVER_TOLERANCE_M = 1e-6


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(event: dict, case: str) -> dict:
    if event["benchmark"]:
        stem = ROOT / "D/d2_analysis/parameter_uncertainty/results/pilot" / f"pilot_{SAMPLE_ID}__{case}__15s"
    else:
        stem = HERE / "results/ensemble" / event["event_id"] / f"{SAMPLE_ID}__{case}__15s"
    array_path, diagnostic_path = stem.with_suffix(".npz"), stem.with_suffix(".json")
    with np.load(array_path, allow_pickle=False) as source:
        data = {key: source[key].copy() for key in source.files}
    data["diagnostics"] = json.loads(diagnostic_path.read_text())
    if event["benchmark"]:
        names = list(data["output_names"])
        selected = [i for i, name in enumerate(names) if str(name).startswith("H_")]
        data["stages"] = data["y"][selected]
        data["stage_names"] = [names[i] for i in selected]
        # Original stored u contains normalized actual opening, with unchanged 0.3 m geometry.
        data["actual_opening_m"] = data["u"][:2] * 0.3
    data["time"] = data["time"].reshape(-1)
    data["downstream_stage"] = data["downstream_stage"].reshape(-1)
    data["downstream_discharge"] = data["downstream_discharge"].reshape(-1)
    assert data["diagnostics"]["full_valid_trajectory_available"], array_path
    assert data["time"][0] == 0 and data["time"][-1] == 144 * 3600, array_path
    assert np.all(np.diff(data["time"]) == 15), array_path
    for field in ("x", "stages", "downstream_stage", "downstream_discharge", "actual_opening_m"):
        assert np.isfinite(data[field]).all(), (array_path, field)
    data["array_file"] = str(array_path.relative_to(ROOT))
    data["array_sha256"] = sha256(array_path)
    data["diagnostic_sha256"] = sha256(diagnostic_path)
    return data


def summarize(event: dict) -> dict:
    wet, dry = load(event, "wetland"), load(event, "no_wetland")
    t, hw, hn = wet["time"], wet["downstream_stage"], dry["downstream_stage"]
    assert np.array_equal(t, dry["time"]), event["event_id"]
    ws, ns = list(wet["state_names"]), list(dry["state_names"])
    iw, ino = int(np.argmax(hw)), int(np.argmax(hn))
    params_w = wet["diagnostics"]["parameter_evaluation"]["effective_physical_values"]
    params_n = dry["diagnostics"]["parameter_evaluation"]["effective_physical_values"]
    assert params_w == params_n, event["event_id"]
    row = {
        "event_id": event["event_id"], "sample_id": SAMPLE_ID, "dt_s": 15,
        "rainfall_total_mm": event["rainfall_total_mm"],
        "maximum_derived_hourly_rain_mm": event["max_hourly_mm"],
        "maximum_24h_rain_mm": event["max_24h_mm"],
        "K_u": params_w["K_u"], "k_ud_1": params_w["k_ud"][0],
        "k_ud_2": params_w["k_ud"][1], "k_ur": params_w["k_ur"], "R_r": params_w["R_r"],
        "wetland_classification": wet["diagnostics"]["classification"],
        "no_wetland_classification": dry["diagnostics"]["classification"],
        "wetland_full_peak_stage_m": float(hw[iw]),
        "wetland_full_peak_time_h": float(t[iw] / 3600),
        "no_wetland_full_peak_stage_m": float(hn[ino]),
        "no_wetland_full_peak_time_h": float(t[ino] / 3600),
        "full_peak_attenuation_no_minus_wet_m": float(hn[ino] - hw[iw]),
        "wetland_maximum_rise_from_initial_m": float(hw[iw] - hw[0]),
        "no_wetland_maximum_rise_from_initial_m": float(hn[ino] - hn[0]),
        "wetland_full_peak_is_initial": bool(iw == 0),
        "no_wetland_full_peak_is_initial": bool(ino == 0),
        "full_peak_delay_h_if_both_after_initial": float((t[iw] - t[ino]) / 3600) if iw > 0 and ino > 0 else None,
        "wetland_full_peak_discharge_m3_s": float(np.max(wet["downstream_discharge"])),
        "no_wetland_full_peak_discharge_m3_s": float(np.max(dry["downstream_discharge"])),
        "actual_gate_opening_arrays_identical": bool(np.array_equal(wet["actual_opening_m"], dry["actual_opening_m"])),
        "actual_gate_opening_max_difference_m": float(np.max(np.abs(wet["actual_opening_m"] - dry["actual_opening_m"]))),
        "shared_upstream_state_max_difference_m3": max(float(np.max(np.abs(wet["x"][ws.index(name)] - dry["x"][ns.index(name)]))) for name in ("S_u", "S_m", "S_d1", "S_d2")),
        "extra_initial_wetland_case_water_m3": float(wet["x"][:, 0].sum() - dry["x"][:, 0].sum()),
    }
    sw = wet["x"][[ws.index("S_w1"), ws.index("S_w2")]]
    summed = sw.sum(axis=0)
    isum = int(np.argmax(summed))
    row.update({"sum_wetland_storage_full_peak_m3": float(summed[isum]),
                "sum_wetland_storage_full_peak_time_h": float(t[isum] / 3600)})
    for branch in (1, 2):
        values = sw[branch - 1]
        peak = int(np.argmax(values))
        row.update({f"wetland_{branch}_storage_full_peak_m3": float(values[peak]),
                    f"wetland_{branch}_storage_full_peak_time_h": float(t[peak] / 3600),
                    f"wetland_{branch}_level_full_peak_m": float(np.max(wet["stages"][list(wet["stage_names"]).index(f"H_w{branch}")])),
                    f"wetland_{branch}_storage_at_no_wetland_full_peak_m3": float(values[ino]),
                    f"wetland_{branch}_final_storage_m3": float(values[-1])})
    spell = event["largest_quantified_wet_spell"]
    origin = datetime.fromisoformat(event["start_utc"])
    start_h = (datetime.fromisoformat(spell["start_utc"]) - origin).total_seconds() / 3600
    end_h = (datetime.fromisoformat(spell["end_utc"]) - origin).total_seconds() / 3600
    eligible = np.flatnonzero(t >= start_h * 3600)
    jw, jn = int(eligible[np.argmax(hw[eligible])]), int(eligible[np.argmax(hn[eligible])])
    js = int(eligible[np.argmax(summed[eligible])])
    later = np.flatnonzero(t >= max(t[jw], t[jn]))
    delta = hw - hn
    above = later[delta[later] > CROSSOVER_TOLERANCE_M]
    jm = int(later[np.argmax(delta[later])])
    row.update({
        "secondary_largest_wet_spell_start_h": start_h,
        "secondary_largest_wet_spell_end_h": end_h,
        "secondary_wetland_peak_after_spell_start_m": float(hw[jw]),
        "secondary_wetland_peak_after_spell_start_time_h": float(t[jw] / 3600),
        "secondary_no_wetland_peak_after_spell_start_m": float(hn[jn]),
        "secondary_no_wetland_peak_after_spell_start_time_h": float(t[jn] / 3600),
        "secondary_peak_attenuation_no_minus_wet_m": float(hn[jn] - hw[jw]),
        "secondary_peak_delay_h": float((t[jw] - t[jn]) / 3600),
        "secondary_wet_minus_no_at_no_wetland_peak_m": float(delta[jn]),
        "secondary_sum_wetland_storage_peak_after_spell_start_m3": float(summed[js]),
        "secondary_sum_wetland_storage_peak_after_spell_start_time_h": float(t[js] / 3600),
        "later_first_saved_wet_above_no_time_h": float(t[above[0]] / 3600) if above.size else None,
        "later_maximum_wet_minus_no_stage_m": float(delta[jm]),
        "later_maximum_wet_minus_no_stage_time_h": float(t[jm] / 3600),
        "final_wet_minus_no_stage_m": float(delta[-1]),
        "wetland_array_file": wet["array_file"], "wetland_array_sha256": wet["array_sha256"],
        "wetland_diagnostic_sha256": wet["diagnostic_sha256"],
        "no_wetland_array_file": dry["array_file"], "no_wetland_array_sha256": dry["array_sha256"],
        "no_wetland_diagnostic_sha256": dry["diagnostic_sha256"],
        "event_forcing_sha256": event["forcing_sha256"],
    })
    return row


def main() -> None:
    manifest_path = HERE / "events/manifest.json"
    events = json.loads(manifest_path.read_text())["events"]
    rows = [summarize(event) for event in events]
    parameters = ("K_u", "k_ud_1", "k_ud_2", "k_ur", "R_r")
    assert all(all(row[key] == rows[0][key] for key in parameters) for row in rows)
    assert all(row["actual_gate_opening_max_difference_m"] <= 1e-12 for row in rows)
    assert all(row[key] == "completed_valid_safety_feasible" for row in rows for key in ("wetland_classification", "no_wetland_classification"))
    output = HERE / "tables/illustrative_mechanism.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    lines = [
        "# Illustrative paired hydrograph and storage check", "",
        "This read-only analysis uses the already completed `mc_00000` parameter draw in all five historical event windows. It does not run an ODE, refit a model, change forcing, or choose a new parameter vector. It is one illustrative draw, not an estimate of event probabilities or a causal explanation of ensemble averages.", "",
        "Run from the project root: `python3 -B D/d2_analysis/historical_storm_robustness/mechanism_check.py`. The script reads the original parameter-study pilot NPZ/JSON files for the September 2022 benchmark and the new compact ensemble NPZ/JSON files for the other four events. It writes only this note and [the detailed CSV](tables/illustrative_mechanism.csv). That CSV records source paths and SHA-256 hashes. No saved simulation result is modified.", "",
        f"Event manifest SHA-256: `{sha256(manifest_path)}`. Analysis script SHA-256: `{sha256(Path(__file__))}`.", "",
        "## Primary full-window quantities", "",
        "Peaks and maximum rises retain the complete 0–144 h window and the first maximum on the native 15 s output grid. A peak at the initial condition remains the primary peak. The physical peak delay is left unavailable if either maximum occurs at time zero; a negative time difference involving the initial maximum must not be interpreted as an advance of the storm response.", "",
        "| Event | Wet peak (m) / time (h) | No-wetland peak (m) / time (h) | Peak attenuation, no minus wet (mm) | Maximum rise wet / no (mm) |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(f"| {row['event_id']} | {row['wetland_full_peak_stage_m']:.6f} / {row['wetland_full_peak_time_h']:.3f} | {row['no_wetland_full_peak_stage_m']:.6f} / {row['no_wetland_full_peak_time_h']:.3f} | {row['full_peak_attenuation_no_minus_wet_m']*1000:.3f} | {row['wetland_maximum_rise_from_initial_m']*1000:.3f} / {row['no_wetland_maximum_rise_from_initial_m']*1000:.3f} |")
    lines += ["", "## Secondary timing and recession description", "",
        "For descriptive storm timing only, a secondary maximum is taken from the start of each event's largest quantified wet spell through 144 h. Spell definitions and dates are inherited unchanged from the weather-only event manifest: positive quantified rain hours are joined unless at least six quantified-zero hours separate them; trace-zero is not verified dry weather. These secondary quantities do not replace the full-window outcome, redefine a flood threshold, or revise any simulation.", "",
        "The later comparison starts after both secondary river maxima. Its crossover time is the first saved node where wetland stage minus no-wetland stage exceeds 0.000001 m (a reporting tolerance, not a hydraulic or safety limit). It is resolved to 15 s and is not an interpolated crossing. A positive later difference is retained as a normal outcome. The largest later difference can occur while small rainfall inputs continue, so this is a hydrograph/recession comparison rather than an assertion that all rainfall has ceased.", "",
        "| Event | Secondary peak attenuation (mm) | Wet peak delay (h) | First later wet > no (h) | Largest later wet − no (mm) | Full sum-wetland storage peak (m³) / time (h) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        cross = row['later_first_saved_wet_above_no_time_h']
        cross_text = "unavailable" if cross is None else f"{cross:.3f}"
        lines.append(f"| {row['event_id']} | {row['secondary_peak_attenuation_no_minus_wet_m']*1000:.3f} | {row['secondary_peak_delay_h']:.3f} | {cross_text} | {row['later_maximum_wet_minus_no_stage_m']*1000:.3f} | {row['sum_wetland_storage_full_peak_m3']:.3f} / {row['sum_wetland_storage_full_peak_time_h']:.3f} |")
    lines += ["", "## Interpretation and limits", "",
        "All ten source simulations are complete and classified safety feasible. Each wet/no pair has identical saved time arrays. The four added-event actual gate-opening arrays are bit-for-bit identical; the benchmark differs by at most 5.69e-15 m (floating-point roundoff). The script verifies agreement within 1e-12 m, purely as an output-comparison tolerance. Shared upland, soil and dam-pool states differ by at most " + f"{max(row['shared_upstream_state_max_difference_m3'] for row in rows):.3g}" + " m³. Thus differences in actual gate operations do not explain the within-pair hydrograph differences in this draw. Both models use the same physical parameter vector. The wetland case contains an additional 5,500 m³ of accepted initial wetland water, so initial recession and the contribution of that stored water must not be presented as a rain-only or initially empty-wetland effect.", "",
        "The saved wetland storage peaks and later river-stage crossovers are consistent with storage and delayed routing: a lower full-window peak can coexist with a higher later stage. The added-event compact arrays do not contain the complete normal-drainage, overflow, groundwater-exchange and ET flux histories, so this check does not quantify how much of the difference is caused by each individual flux.", "",
        "For August 2023, the storm-period secondary attenuation is 100.346 mm versus 113.002 mm in September 2022, while their primary full-window attenuations are 67.090 mm and 113.002 mm. The August wetland storm response remains below the initial river stage. Its derived maximum clock-hour rainfall (16.693 mm) is comparable with the benchmark (15.906 mm) despite much lower total rainfall (58.9 versus 145.3 mm). These observations illustrate why cumulative rainfall alone cannot determine a peak-stage difference. They do not establish why ensemble averages may be similar, isolate an individual weather driver, or generalize beyond this sample.", "",
        "All four meteorological variables vary jointly across the selected events. The benchmark retains its original mass-preserving 120 s reconstruction, whereas added events use the preserved Stage C hourly reconstruction. The selected warm complete windows are not a random sample of storms. Initial-state dominance, incomplete earlier-event coverage, and differences in forcing reconstruction remain material limits. No flood threshold is introduced by this analysis.", "",
    ]
    (HERE / "mechanism_check.md").write_text("\n".join(lines))
    print(json.dumps({"events": len(rows), "paired_complete": len(rows), "exact_gate_array_matches": sum(row["actual_gate_opening_arrays_identical"] for row in rows), "maximum_gate_difference_m": max(row["actual_gate_opening_max_difference_m"] for row in rows), "csv": str(output), "note": str(HERE / 'mechanism_check.md')}))


if __name__ == "__main__":
    main()
