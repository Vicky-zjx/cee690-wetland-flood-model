"""Paired meteorological forcing around the accepted Stage B ODEs.

No hydrologic equation, decision, geometry, boundary, or initial state is changed.
Hourly rain/radiation are interval constants; official ode4u is restarted at
hour boundaries so its midpoint input averaging cannot smear those quantities.
"""
from __future__ import annotations

import copy
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import tempfile
import traceback

# Importing the historical helper sets up Matplotlib; keep its cache outside A.
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "cee690-stage-c-matplotlib"))
import numpy as np
import pandas as pd
from A import wetland_flood_model as a2
from B import no_spillway_model as wet
from B import wetland_counterfactual as dry
from B.no_spillway_study import load_trajectory
from B.official_opt import load_official
from B.wetland_only_study import load_study

ROOT = Path(__file__).resolve().parents[2]
WEATHER_COLUMNS = ["precipitation_mm_per_h", "temperature_C", "solar_W_per_m2", "wind_m_per_s"]
CASE_ALIASES = {"no_wetland": "no_wetland", "bypass": "no_wetland", "case1": "no_wetland",
                "wetland": "wetland", "restoration": "wetland", "case2": "wetland"}


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _array_hash(*arrays):
    h = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(array)
        h.update(str(value.shape).encode())
        h.update(value.dtype.str.encode())
        h.update(value.tobytes())
    return h.hexdigest()


@lru_cache(maxsize=1)
def _accepted_inputs():
    config, spec, gates = load_study()
    official = load_official()
    return config, spec, gates, official


def load_fixed_config(options=None):
    """Load frozen decisions from saved B; options are numerical scope only.

    Recognized options: dt_s, duration_hours, design_source, min_temperature_C.
    An already loaded dict may be passed directly to run_case. Nothing is optimized.
    """
    options = {} if options is None else dict(options)
    allowed = {"dt_s", "duration_hours", "design_source", "min_temperature_C"}
    unknown = set(options) - allowed
    if unknown:
        raise ValueError("Unknown adapter settings: " + ", ".join(sorted(unknown)))
    original, spec, gates, official = _accepted_inputs()
    config = copy.deepcopy(original)
    dt = float(options.get("dt_s", config["numerics"]["dt_s"]))
    duration = float(options.get("duration_hours", config["numerics"]["duration_s"] / 3600)) * 3600
    if not np.isfinite(dt) or dt <= 0 or not np.isclose(3600 / dt, round(3600 / dt), rtol=0, atol=1e-10):
        raise ValueError("dt_s must divide one hour exactly; no partial rain intervals.")
    if not np.isfinite(duration) or duration <= 0 or duration > original["numerics"]["duration_s"] or duration % 3600:
        raise ValueError("Duration must be whole hours within the accepted 144-hour gate schedule.")
    minimum = float(options.get("min_temperature_C", 0))
    if minimum != 0:
        raise ValueError("The accepted liquid-water applicability limit is T > 0 C; not redefined here.")
    path = Path(options.get("design_source", "B/results/wetland_only/wetland_best_120s.npz"))
    if not path.is_absolute():
        path = ROOT / path
    design = load_trajectory(path)
    if not np.array_equal(design["theta"][:-2], gates):
        raise ValueError("Saved design does not use the accepted fixed gate schedule.")
    pair = list(design["pair"])
    covers = np.asarray(design["theta"][-2:], dtype=float).copy()
    if len(pair) != 2 or any(label not in config["categories"]["endpoints"] for label in pair) or np.any((covers < 0) | (covers > 1)):
        raise ValueError("Saved wetland design outside accepted domains.")
    config["numerics"]["dt_s"] = dt
    config["numerics"]["duration_s"] = duration
    initial = np.array(config["initial_conditions"]["x0_m3"], dtype=float)
    sources = [ROOT / "A/wetland_flood_model.py", ROOT / "B/no_spillway_model.py",
               ROOT / "B/wetland_counterfactual.py", ROOT / "B/gate_hydraulics.py",
               ROOT / "A/work/multivarious/multivarious/ode/ode4u.py", ROOT / "B/no_spillway_config.json",
               ROOT / "B/wetland_only_config.json",
               ROOT / "B/results/no_spillway_revision/reference.npz", path]
    manifest = dict(
        source_sha256={str(p.relative_to(ROOT)): sha256_file(p) for p in sources},
        solver={k: official[k] for k in ("version", "revision", "ode4u_signature")},
        model_scope="Accepted Stage B no-independent-dam-spillway models; wetland overflow retained",
        design_categories=pair, design_covers=covers.tolist(), fixed_gate_values=gates.tolist(),
        gate1_hours=config["controls"]["gate1_hours"], gate2_hours=config["controls"]["gate2_hours"],
        boundaries=copy.deepcopy(config["forcing"]["proposed_fixed_boundaries"]),
        initial_wetland_m3=initial.tolist(), initial_no_wetland_m3=initial[dry.SHARED].tolist(),
        initial_extra_wetland_water_m3=float(initial[2:4].sum()),
        initial_condition_note="Matched shared storages. Wetland case retains 5500 m3 additional accepted wetland water; no spin-up or redistribution.",
        datum=config["geography"]["vertical_datum"], dt_s=dt, duration_s=duration,
        conditional_on="Fixed models, parameters, boundary paths, gate schedule, wetland design and initial assumptions; meteorology only varies.")
    return dict(_stage_c_fixed=True, model_config=config, fixed_gates=gates.copy(), pair=pair,
                covers=covers, theta=np.r_[gates, covers], manifest=manifest, official=official)


def read_forcing(forcing):
    """Read one already generated weather history, preserving its saved identity."""
    path = Path(forcing) if isinstance(forcing, (str, Path)) else None
    frame = pd.read_csv(path, float_precision="round_trip") if path is not None else forcing.copy(deep=True)
    required = ["timestamp_utc"] + WEATHER_COLUMNS
    if not set(required).issubset(frame.columns):
        raise ValueError("Hourly forcing requires columns: " + ", ".join(required))
    raw_times = pd.to_datetime(frame["timestamp_utc"])
    if getattr(raw_times.dt, "tz", None) is None:
        raise ValueError("Forcing timestamps must declare their timezone; naive timestamps are rejected.")
    times = raw_times.dt.tz_convert("UTC")
    if len(times) == 0 or not times.is_monotonic_increasing or times.duplicated().any():
        raise ValueError("Require nonempty, unique increasing hourly timestamps.")
    if len(times) > 1 and not np.all(np.diff(times.astype("int64").to_numpy()) == 3600 * 10**9):
        raise ValueError("Missing or non-hourly intervals in saved forcing; no gap filling in the adapter.")
    if any((times.dt.minute != 0) | (times.dt.second != 0) | (times.dt.microsecond != 0)):
        raise ValueError("Each forcing row must start exactly on a UTC hour.")
    values = frame[WEATHER_COLUMNS].to_numpy(dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("Forcing contains missing/nonfinite values; missing precipitation is never zero-filled.")
    if np.any(values[:, [0, 2, 3]] < 0):
        raise ValueError("Precipitation, radiation and wind must be nonnegative; no clipping.")
    frame = frame.copy()
    frame["timestamp_utc"] = times
    digest = _array_hash(times.astype("int64").to_numpy(), values)
    return frame, dict(forcing_file=None if path is None else str(path.resolve()),
                       forcing_file_sha256=None if path is None else sha256_file(path),
                       forcing_array_sha256=digest,
                       timestamp_meaning="UTC beginning of a one-hour accumulation/mean interval",
                       forcing_rows=len(frame))


def hourly_input_segments(frame, fixed):
    """Yield official-solver inputs; rain depth and solar energy are conservative.

    P is mm per hour / 3.6e6 -> m/s, constant on each interval. R_s is
    interval-mean W/m2, also constant. T and wind are linear between hourly
    centers with endpoint hold over the first/last half hour. Gates keep their
    accepted linear knot schedule; all other boundaries remain fixed.
    """
    cfg = fixed["model_config"]
    duration = cfg["numerics"]["duration_s"]
    dt = cfg["numerics"]["dt_s"]
    if len(frame) != round(duration / 3600):
        raise ValueError(f"Expected {round(duration / 3600)} full hourly intervals, found {len(frame)}.")
    weather = frame[WEATHER_COLUMNS].to_numpy(dtype=float)
    centers = np.arange(len(frame)) * 3600 + 1800
    bounds = cfg["forcing"]["proposed_fixed_boundaries"]
    n1 = len(cfg["controls"]["gate1_hours"])
    for index, row in enumerate(weather):
        time = index * 3600 + np.arange(round(3600 / dt) + 1) * dt
        u = np.vstack([
            np.interp(time, np.asarray(cfg["controls"]["gate1_hours"]) * 3600, fixed["fixed_gates"][:n1]),
            np.interp(time, np.asarray(cfg["controls"]["gate2_hours"]) * 3600, fixed["fixed_gates"][n1:]),
            np.full(len(time), row[0] / 3.6e6),
            np.interp(time, centers, weather[:, 1]),
            np.full(len(time), row[2]),
            np.interp(time, centers, weather[:, 3]),
            *[np.full(len(time), bounds[k]) for k in ("Q_ext_m3_s", "h_d_m", "h_g_m", "B")]])
        yield time, u


def forcing_conservation(frame, fixed):
    segments = list(hourly_input_segments(frame, fixed))
    depth = sum(np.trapezoid(u[a2.U["P"]], t) for t, u in segments)
    energy = sum(np.trapezoid(u[a2.U["R_s"]], t) / 3600 for t, u in segments)
    source_depth = float(frame[WEATHER_COLUMNS[0]].sum() / 1000)
    source_energy = float(frame[WEATHER_COLUMNS[2]].sum())
    return dict(source_precipitation_m=source_depth, integrated_precipitation_m=float(depth),
                precipitation_error_m=float(depth - source_depth), source_solar_Wh_m2=source_energy,
                integrated_solar_Wh_m2=float(energy), solar_error_Wh_m2=float(energy - source_energy),
                mapping="Official ode4u hourly restarts; constant P/R_s per segment, linear center-time T/wind",
                pass_conservation=bool(abs(depth-source_depth) <= 1e-12 and abs(energy-source_energy) <= 1e-7))


def _parameters(fixed):
    config = fixed["model_config"]
    params = wet.base_parameters(config)
    params.rho = fixed["covers"].copy()
    params.tau = list(fixed["pair"])
    for label in fixed["pair"]:
        params.R_veg[label] = np.asarray(config["categories"]["endpoints"][label]["R_veg"], dtype=float)
        params.K_veg[label] = np.asarray(config["categories"]["endpoints"][label]["K_veg"], dtype=float)
    initial = np.asarray(config["initial_conditions"]["x0_m3"], dtype=float)
    params.u_depth_init = initial[0] / params.A_u
    params.soil_fraction_init = initial[1] / params.C_m
    params.w_depth_init = initial[2:4] / params.A_w
    params.d_depth_init = initial[4:6] / params.A_d
    params.H_g_init = a2.stage(initial[6], params.S_y * params.A_g, params.z_g)
    params.r_depth_init = initial[7] / params.A_r
    return params, initial, wet.gate_parameters(config)


def _safety_summary(case, y, names, config, time):
    rows = []
    for kind in ("dam", "wetland") if case == "wetland" else ("dam",):
        prefix = "d" if kind == "dam" else "w"
        for i, bound in enumerate(config["safety"][kind + "_m"]):
            values = y[names.index(f"H_{prefix}{i+1}")]
            k = int(np.argmax(values))
            rows.append(dict(variable=f"H_{prefix}{i+1}", threshold_m=bound,
                             peak_m=float(values[k]), peak_time_s=float(time[k]),
                             margin_m=float(bound-values[k]), exceeded=bool(values[k] > bound)))
    return rows


def _integrate(case, segments, fixed, diagnostics):
    params, initial8, gates = _parameters(fixed)
    initial = initial8.copy() if case == "wetland" else initial8[dry.SHARED].copy()
    names = list(a2.OUTPUT_NAMES if case == "wetland" else dry.OUTPUT_NAMES)
    states = list(a2.STATE_NAMES if case == "wetland" else dry.STATE_NAMES)
    raw_rhs = wet.make_b2_rhs(gates) if case == "wetland" else dry.make_rhs(gates)
    last_rhs_time = 0.0
    last_rhs_state = None
    last_rhs_input = None
    derivative_balance = 0.0
    minimum_g5 = [float("inf"), float("inf")]
    def checked_rhs(t, x, u, p):
        nonlocal last_rhs_time, last_rhs_state, last_rhs_input, derivative_balance
        last_rhs_time = float(t)
        last_rhs_state, last_rhs_input = x, u
        diagnostics["rhs_evaluations"] += 1
        dx, y = raw_rhs(t, x, u, p)
        diagnostics["last_successful_rhs_time_s"] = float(t)
        derivative_balance = max(derivative_balance, abs(float(sum(dx)-y[names.index("external_net")])))
        for i in range(2):
            Hd = y[names.index(f"H_d{i+1}")]
            tail = y[names.index(f"H_w{i+1}")] if case == "wetland" else y[names.index("H_r")]
            H = Hd - gates[i].gate_sill_z
            if u[i] > 0 and H > 0 and Hd > tail:
                minimum_g5[i] = min(minimum_g5[i], float(H-1.25*gates[i].gate_max_opening*u[i]))
        return dx, y
    chunks = []; total_net = 0.0; total_inflow = 0.0
    max_balance = 0.0; previous_end = 0.0
    diagnostics.update(state_names=states, initial_storage_m3=initial.tolist(),
                       rainfall_area_m2=float(params.A_rain), former_wetland_rainfall_area_m2=float(sum(params.A_w)),
                       former_wetland_rainfall_routing="direct to river" if case == "no_wetland" else "wetland storages",
                       safety_is_a_reported_outcome=True, completed_until_s=0.0,
                       ode_started=False, ode_segment_calls=0, rhs_evaluations=0,
                       expected_output_samples=int(round(fixed['manifest']['duration_s']/fixed['manifest']['dt_s']))+1,
                       downstream_output_columns={name:names.index(name) for name in ('H_r','Q_out')})
    def record_exception(error, stage):
        """Record the actual failure, without changing the RHS or its handling."""
        cause=error
        while cause.__cause__ is not None:
            cause=cause.__cause__
        stack=traceback.extract_tb(cause.__traceback__)
        location=stack[-1] if stack else None
        diagnostics.update(failure_stage=stage, exception_type=type(error).__name__,
                           root_exception_type=type(cause).__name__, root_exception_message=str(cause),
                           exception_traceback=''.join(traceback.format_exception(error)),
                           failure_location=None if location is None else dict(file=location.filename,line=location.lineno,function=location.name),
                           last_completed_state_m3=initial.tolist(),
                           failure_rhs_state_m3=None if last_rhs_state is None else np.asarray(last_rhs_state).tolist(),
                           failure_rhs_inputs=None if last_rhs_input is None else np.asarray(last_rhs_input).tolist())
    try:
        for index, (time, u) in enumerate(segments):
            # A full-shape validator checks every input once before integration;
            # its p.time length is temporarily the hourly segment length.
            if index == 0:
                original_duration = params.duration
                params.duration = time[-1] - time[0]
                if case == "wetland":
                    wet.validate(params, u, gates)
                else:
                    dry.validate_six(params, u, initial, gates)
                params.duration = original_duration
            diagnostics['ode_started']=True
            diagnostics['ode_segment_calls']+=1
            t, x, dx, y = fixed["official"]["ode4u"](checked_rhs, time, initial, u, params)
            if not all(np.all(np.isfinite(v)) for v in (x, dx, y)) or x.shape != (len(states), len(t)):
                raise FloatingPointError("Official ode4u returned incomplete/nonfinite arrays")
            if np.any(x < 0) or np.any(x[1] > params.C_m):
                raise FloatingPointError("Computed storage is outside its unchanged physical domain")
            net = a2.cumulative_trapezoid(y[names.index("external_net")], t)
            balance = x.sum(axis=0) - np.sum(initial8 if case == "wetland" else initial8[dry.SHARED]) - (total_net + net)
            max_balance = max(max_balance, float(np.max(np.abs(balance))))
            total_net += float(net[-1])
            total_inflow += float(np.trapezoid(y[names.index("P_total")]+u[a2.U["Q_ext"]]+y[names.index("Q_gx_in")], t))
            chunks.append(dict(time=t, x=x, dx=dx, y=y, u=u, balance=balance))
            initial = x[:, -1].copy()
            previous_end = float(t[-1])
            diagnostics["completed_until_s"] = previous_end
    except wet.HydraulicApplicabilityFailure as error:
        record_exception(error, 'ode_rhs_model_applicability')
        diagnostics.update(failure_reason="unsupported_gate_regime", first_violation=error.diagnostic,
                           failure_message=str(error), attempted_until_s=last_rhs_time)
        return _empty_result(case, "applicability_excluded", diagnostics, states)
    except wet.TrajectoryFailure as error:
        if "Applicability domain exceeded" not in str(error):
            raise
        record_exception(error, 'ode_rhs_model_applicability')
        diagnostics.update(failure_reason="hydrologic_geometry_applicability", failure_message=str(error),
                           attempted_until_s=last_rhs_time)
        return _empty_result(case, "applicability_excluded", diagnostics, states)
    except FloatingPointError as error:
        record_exception(error, 'ode_or_state_validation')
        diagnostics.update(failure_reason="nonfinite_or_inadmissible_storage", failure_message=str(error),
                           attempted_until_s=last_rhs_time)
        return _empty_result(case, "numerical_failure", diagnostics, states)
    # Duplicate boundaries carry the SAME states. Keep the next interval's
    # right-hand input/output; integrate net flux per segment before merging.
    arrays = {}
    for key in ("time", "x", "dx", "y", "u", "balance"):
        arrays[key] = np.concatenate([chunk[key][..., :-1] for chunk in chunks[:-1]] + [chunks[-1][key]], axis=-1)
    scale = total_inflow if total_inflow > 0 else np.sum(initial8 if case == "wetland" else initial8[dry.SHARED])
    relative = max_balance / scale if scale > 0 else (0. if max_balance == 0 else float("inf"))
    safety = _safety_summary(case, arrays["y"], names, fixed["model_config"], arrays["time"])
    diagnostics.update(balance_max_m3=max_balance, balance_relative=float(relative),
                       balance_relative_tolerance=float(params.balance_rtol), derivative_balance_max_m3s=derivative_balance,
                       balance_method="Cumulative trapezoidal external net flux within each official-solver segment; jumps not smeared",
                       safety=safety, safety_violation=any(r["exceeded"] for r in safety),
                       minimum_G5_margin_m=[None if np.isinf(v) else v for v in minimum_g5],
                       mapped_input_array_sha256=_array_hash(arrays["time"], arrays["u"]),
                       completed_until_s=previous_end, solver_termination="Completed fixed interval",
                       min_storage_m3=float(arrays["x"].min()))
    if relative >= params.balance_rtol:
        diagnostics.update(failure_reason="excessive_mass_balance_residual", failure_stage="post_integration_balance_validation",
                           failure_message="Unchanged model balance_rtol exceeded")
        return _empty_result(case, "numerical_failure", diagnostics, states)
    return dict(case_name=case, **arrays, state_names=states, output_names=names,
                downstream_stage=arrays["y"][names.index("H_r")],
                downstream_discharge=arrays["y"][names.index("Q_out")],
                solver_status="completed", valid=True, diagnostics=diagnostics)


def _empty_result(case, status, diagnostics, states=None):
    states = list(a2.STATE_NAMES if case == "wetland" else dry.STATE_NAMES) if states is None else states
    return dict(case_name=case, time=np.empty(0), timestamp_utc=[], downstream_stage=np.empty(0),
                downstream_discharge=np.empty(0), x=np.empty((len(states or []), 0)),
                state_names=states or [], solver_status=status, valid=False, diagnostics=diagnostics)


def run_case(case_name, forcing, fixed_config=None):
    """Run one existing dynamic model on the SAME saved paired weather history.

    Safety exceedance is a completed valid simulation, distinct from a G5/model
    applicability exclusion or numerical failure. Invalid runs have no full peak.
    Unexpected errors propagate; no random draws occur anywhere in this adapter.
    """
    if case_name not in CASE_ALIASES:
        raise ValueError("case_name must be no_wetland or wetland")
    case = CASE_ALIASES[case_name]
    fixed = fixed_config if isinstance(fixed_config, dict) and fixed_config.get("_stage_c_fixed") else load_fixed_config(fixed_config)
    frame, diagnostics = read_forcing(forcing)
    diagnostics.update(ode_started=False, ode_segment_calls=0, rhs_evaluations=0)
    diagnostics["fixed_configuration"] = fixed["manifest"]
    diagnostics["forcing_conservation"] = forcing_conservation(frame, fixed)
    if not diagnostics["forcing_conservation"]["pass_conservation"]:
        raise RuntimeError("Conservative input mapping failed before integration")
    cold = frame["temperature_C"] <= 0
    if cold.any():
        first = int(np.flatnonzero(cold.to_numpy())[0])
        diagnostics.update(failure_reason="nonfreezing_liquid_water_domain", completed_until_s=0.,
                           failure_stage="preflight_model_applicability",
                           first_violation=dict(timestamp_utc=frame["timestamp_utc"].iloc[first].isoformat(),
                                                time_s=first*3600, temperature_C=float(frame["temperature_C"].iloc[first])),
                           nonpositive_temperature_hours=int(cold.sum()),
                           failure_message="Accepted hydrologic model requires T > 0 C; no snowmelt or clipping")
        return _empty_result(case, "applicability_excluded", diagnostics)
    result = _integrate(case, list(hourly_input_segments(frame, fixed)), fixed, diagnostics)
    if result["valid"]:
        result["timestamp_utc"] = [value.isoformat() for value in frame["timestamp_utc"].iloc[0] + pd.to_timedelta(result["time"], unit="s")]
    return result


def deterministic_regression(output_dir=None):
    """Fresh 120 s official-solver reproduction of unchanged saved B cases.

    The historical forcing already has a documented conservative 120 s linear
    representation; it is intentionally NOT aggregated/remapped to hourly rain.
    One common full-grid segment verifies adapter math against original wrappers.
    """
    fixed = load_fixed_config()
    config = fixed["model_config"]
    forcing = wet.fixed_forcing(config, wet.base_parameters(config))
    time = wet.base_parameters(config).time
    u = dry.inputs(config, fixed["fixed_gates"], time, forcing)
    records = []
    for case, saved_path in [("no_wetland", ROOT / "B/results/wetland_counterfactual/bypass_120s.npz"),
                             ("wetland", ROOT / "B/results/wetland_only/wetland_best_120s.npz")]:
        saved = load_trajectory(saved_path)
        before = sha256_file(saved_path)
        actual = _integrate(case, [(time, u)], fixed, {"forcing_mapping": "Unchanged historical 120-second linear knots"})
        if not actual["valid"]:
            raise RuntimeError("Deterministic regression did not complete: " + str(actual["diagnostics"]))
        differences = {key: float(np.max(np.abs(actual[key]-saved[key]))) for key in ("x", "dx", "y", "u") if key in saved}
        row = dict(case_name=case, saved_trajectory=str(saved_path.relative_to(ROOT)), saved_sha256=before,
                   numerical_differences=differences, not_saved_columns=[key for key in ("dx",) if key not in saved],
                   storage_atol_m3=1e-7, output_atol=1e-9,
                   peak_stage_m=float(max(actual["downstream_stage"])),
                   pass_regression=bool(differences["x"] <= 1e-7 and max(differences[k] for k in ("dx", "y", "u") if k in differences) <= 1e-9),
                   saved_file_unchanged=before == sha256_file(saved_path))
        records.append(row)
    source_unchanged = all(sha256_file(ROOT / path) == digest for path, digest in fixed["manifest"]["source_sha256"].items())
    report = dict(records=records, source_files_unchanged=source_unchanged,
                  passed=bool(source_unchanged and all(r["pass_regression"] and r["saved_file_unchanged"] for r in records)),
                  note="Actual fresh official ode4u runs; no original code/config/trajectory files edited.")
    if output_dir is not None:
        directory = Path(output_dir); directory.mkdir(parents=True, exist_ok=True)
        (directory / "deterministic_regression.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
