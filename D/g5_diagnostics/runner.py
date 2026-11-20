"""Read-only reuse of the recovered Stage C adapter for Stage D diagnostics.

No fitting, optimization, random draws, historical writes, hydraulic clipping,
or changes to the Stage B equations occur here. The observer wraps ode4u's RHS
calls, including intermediate RK stages; it does not alter their arguments or
return values. JSON writers are intentionally left to the calling D workflow.
"""
from __future__ import annotations

import importlib.util
import math
import os
from pathlib import Path
import sys
import tempfile

sys.dont_write_bytecode = True
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "cee690-stage-d-matplotlib"))
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ADAPTER_PATH = Path(__file__).resolve().parent / "recovered/model_adapters.py"
_ADAPTER = None


def configure(project_root, adapter_path=None):
    """Set paths for an imported temporary draft; installed D defaults suffice."""
    global ROOT, ADAPTER_PATH, _ADAPTER
    ROOT = Path(project_root).resolve()
    ADAPTER_PATH = Path(adapter_path).resolve() if adapter_path else ROOT / "D/g5_diagnostics/recovered/model_adapters.py"
    _ADAPTER = None


def adapter():
    """Load the recovered file without editing it; point its data root at repo."""
    global _ADAPTER
    if _ADAPTER is None:
        if not ADAPTER_PATH.is_file():
            raise FileNotFoundError(f"Recovered original adapter is required: {ADAPTER_PATH}")
        if str(ROOT) not in sys.path:
            sys.path.insert(0, str(ROOT))
        spec = importlib.util.spec_from_file_location("_stage_d_recovered_adapter", ADAPTER_PATH)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.ROOT = ROOT
        _ADAPTER = module
    return _ADAPTER


def json_clean(value):
    """JSON-safe values; missing/nonfinite quantities are null, never zero."""
    if isinstance(value, dict):
        return {str(k): json_clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_clean(value.tolist())
    if isinstance(value, np.generic):
        return json_clean(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    return value


def _design_path(design):
    names = {"best": "wetland_best_120s.npz", "reference": "fixed_reference_120s.npz"}
    if design not in names:
        raise ValueError("design must be 'best' or 'reference'.")
    path = ROOT / "B/results/wetland_only" / names[design]
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def load_fixed(dt=120, gate_scale=1, design="best", constant_opening=None):
    """Load accepted initial states/design; optionally declare one scaled schedule.

    Scaling is an explicit common open-loop operating alternative, never an
    adaptive response or silently clipped opening. Values outside [0,1] fail.
    The only constant alternative permitted in this bounded diagnostic study
    is the existing 0.005 base opening at every original knot on both gates.
    """
    module = adapter()
    fixed = module.load_fixed_config(dict(dt_s=float(dt), design_source=str(_design_path(design))))
    fixed["official"] = dict(fixed["official"])  # Do not mutate its cached dictionary.
    scale = float(gate_scale)
    original = fixed["fixed_gates"].copy()
    changed = original * scale
    if not np.isfinite(scale) or scale < 0 or not np.all(np.isfinite(changed)) or np.any((changed < 0) | (changed > 1)):
        raise ValueError("Declared schedule scale must yield finite gate fractions in [0,1]; no clipping.")
    if constant_opening is not None:
        constant_opening = float(constant_opening)
        if constant_opening != 0.005 or scale != 1:
            raise ValueError("This bounded study permits only constant_opening=0.005 with gate_scale=1; alternatives are not combined.")
        if float(original.min()) != constant_opening:
            raise ValueError("The proposed constant opening is not the original schedule's saved base opening.")
        changed = np.full_like(original, constant_opening)
    schedule_kind = ("constant_base_opening" if constant_opening is not None else
                     "original" if scale == 1 else "uniform_scale")
    fixed["fixed_gates"] = changed
    fixed["theta"] = np.r_[changed, fixed["covers"]]
    fixed["manifest"].update(
        original_fixed_gate_values=original.tolist(), fixed_gate_values=changed.tolist(),
        schedule_scale=scale if constant_opening is None else None,
        schedule_is_original=bool(scale == 1 and constant_opening is None), design_source_label=design,
        schedule_kind=schedule_kind, constant_opening=constant_opening,
        schedule_description=("Explicit weather-independent constant base-opening alternative: both gate fractions remain 0.005 at every original knot; no search or state feedback"
                              if constant_opening is not None else
                              "Original accepted Stage B fixed schedule" if scale == 1 else
                              "Explicit common fixed operating alternative: every original knot multiplied by declared scale; knot times unchanged"),
        recovered_adapter_sha256=module.sha256_file(ADAPTER_PATH))
    return fixed


class RHSObserver:
    """Observe every actual official-solver RHS call, including failed calls."""
    def __init__(self, case, fixed):
        self.case = case
        self.fixed = fixed
        self.names = list(adapter().a2.OUTPUT_NAMES if case == "wetland" else adapter().dry.OUTPUT_NAMES)
        self.states = list(adapter().a2.STATE_NAMES if case == "wetland" else adapter().dry.STATE_NAMES)
        self.gates = adapter().wet.gate_parameters(fixed["model_config"])
        self.calls = self.successes = self.segment_calls = 0
        self.safety_attempted = {}
        self.safety_successful = {}
        self.geometry = {}
        self.g5 = [None, None]
        self.tailwater_clearance = [None, None]
        self.minimum_state = {name: None for name in self.states}
        self.soil_excess_max = None
        self.derivative_error = None
        self.declared_net_error = None
        self.first_error = None
        self.last_time = None

    @staticmethod
    def _maximum(target, name, level, bound, t):
        if not np.isfinite(level):
            return
        if name not in target or level > target[name]["peak_m"]:
            target[name] = dict(variable=name, threshold_m=float(bound), peak_m=float(level),
                                peak_time_s=float(t), margin_m=float(bound-level), exceeded=bool(level > bound))

    def _levels(self, x, p):
        levels = {}
        indices = {name: i for i, name in enumerate(self.states)}
        for i in range(2):
            levels[f"H_d{i+1}"] = float(x[indices[f"S_d{i+1}"]] / p.A_d[i] + p.z_d[i])
            if self.case == "wetland":
                levels[f"H_w{i+1}"] = float(x[indices[f"S_w{i+1}"]] / p.A_w[i] + p.z_w[i])
        levels["H_g"] = float(x[indices["S_g"]] / (p.S_y*p.A_g) + p.z_g)
        levels["H_r"] = float(x[indices["S_r"]] / p.A_r + p.z_r)
        return levels

    def _safety(self, target, levels, t):
        for kind, prefix in (("dam", "d"), ("wetland", "w")):
            if kind == "wetland" and self.case != "wetland":
                continue
            for i, bound in enumerate(self.fixed["model_config"]["safety"][kind+"_m"]):
                name = f"H_{prefix}{i+1}"
                self._maximum(target, name, levels[name], bound, t)

    def wrap_solver(self, original_solver):
        def observed_solver(rhs, time, initial, inputs, params):
            self.segment_calls += 1
            segment_rhs_calls = 0
            def observed_rhs(t, x, u, p):
                nonlocal segment_rhs_calls
                segment_rhs_calls += 1
                self.calls += 1
                self.last_time = float(t)
                levels = self._levels(x, p)
                self._safety(self.safety_attempted, levels, t)
                for name, val in zip(self.states, x):
                    if np.isfinite(val) and (self.minimum_state[name] is None or val < self.minimum_state[name]):
                        self.minimum_state[name] = float(val)
                excess = float(x[self.states.index("S_m")]-p.C_m)
                self.soil_excess_max = excess if self.soil_excess_max is None else max(excess, self.soil_excess_max)
                bounds = [("groundwater_ceiling", "H_g", p.H_g_valid_max),
                          ("groundwater_soil_bottom", "H_g", p.soil_bottom),
                          ("river_geometry", "H_r", p.H_r_valid_max)]
                if self.case == "wetland":
                    bounds.append(("groundwater_wetland_bed", "H_g", min(p.z_w)))
                for label, name, bound in bounds:
                    self._maximum(self.geometry, label, levels[name], bound, t)
                    if label in self.geometry:
                        self.geometry[label]["variable"] = name
                        self.geometry[label]["invalid_at_equality"] = True
                        self.geometry[label]["strict_margin_m"] = self.fixed["model_config"]["numerics"]["strict_margin_m"]
                        self.geometry[label]["applicability_violated"] = self.geometry[label]["peak_m"] >= bound
                        self.geometry[label]["margin_after_configured_strict_margin_m"] = (self.geometry[label]["margin_m"]-
                            self.fixed["model_config"]["numerics"]["strict_margin_m"])
                for i, gate in enumerate(self.gates):
                    upstream = levels[f"H_d{i+1}"]
                    tail = levels[f"H_w{i+1}"] if self.case == "wetland" else levels["H_r"]
                    head = upstream-gate.gate_sill_z
                    clearance = gate.gate_sill_z-tail
                    self.tailwater_clearance[i] = clearance if self.tailwater_clearance[i] is None else min(clearance, self.tailwater_clearance[i])
                    if u[i] > 0 and head > 0 and upstream > tail:
                        margin = float(head-1.25*gate.gate_max_opening*u[i])
                        if self.g5[i] is None or margin < self.g5[i]["margin_m"]:
                            opening = float(gate.gate_max_opening*u[i])
                            self.g5[i] = dict(branch=i+1, margin_m=margin, residual_m=-margin,
                                             time_s=float(t), upstream_stage_m=upstream, tailwater_stage_m=tail,
                                             water_depth_m=head, opening_height_m=opening, H_over_a=head/opening)
                try:
                    dx, y = rhs(t, x, u, p)
                except Exception as exc:
                    if self.first_error is None:
                        self.first_error = dict(time_s=float(t), call_index=self.calls,
                                                segment_index=self.segment_calls-1,
                                                segment_rhs_index=segment_rhs_calls,
                                                rk_evaluation=("initial" if segment_rhs_calls == 1 else
                                                               ["k2", "k3", "k4", "accepted_endpoint"][(segment_rhs_calls-2) % 4]),
                                                exception_type=type(exc).__name__, message=str(exc),
                                                stage_levels_m=levels,
                                                diagnostic=getattr(exc, "diagnostic", None))
                    raise
                self.successes += 1
                self._safety(self.safety_successful, levels, t)
                output = dict(zip(self.names, y))
                external = float(u[adapter().a2.U["P"]]*p.A_rain + u[adapter().a2.U["Q_ext"]]
                                 + output["Q_gx"]-output["Q_out"]-output["E_tot"])
                error = abs(float(np.sum(dx))-external)
                net_error = abs(float(output["external_net"])-external)
                self.derivative_error = error if self.derivative_error is None else max(error, self.derivative_error)
                self.declared_net_error = net_error if self.declared_net_error is None else max(net_error, self.declared_net_error)
                return dx, y
            return original_solver(observed_rhs, time, initial, inputs, params)
        return observed_solver

    def summary(self):
        return json_clean(dict(
            ode_segment_calls=self.segment_calls, rhs_attempts=self.calls, rhs_successes=self.successes,
            rhs_failures=self.calls-self.successes, last_rhs_time_s=self.last_time,
            observation_scope="Every official ode4u RHS/RK stage; separate preflight validator RHS is not counted. Attempted maxima may include the first rejected stage.",
            safety_successful_rhs=list(self.safety_successful.values()),
            safety_attempted_rhs=list(self.safety_attempted.values()),
            wetland_safety_applicable=self.case == "wetland",
            geometry_attempted_rhs=list(self.geometry.values()),
            min_storage_by_state_m3=self.minimum_state, max_soil_capacity_residual_m3=self.soil_excess_max,
            minimum_G5_by_branch=self.g5, minimum_gate_sill_above_tailwater_m=self.tailwater_clearance,
            derivative_physical_mass_identity_max_abs_m3_s=self.derivative_error,
            declared_external_net_identity_max_abs_m3_s=self.declared_net_error,
            physical_mass_identity="sum(dx) = P*A_rain + Q_ext + Q_gx - Q_out - E_tot; internal transfers cancel once",
            first_rhs_exception=self.first_error))


def _attach_observer(case, fixed):
    observer = RHSObserver(case, fixed)
    fixed["official"]["ode4u"] = observer.wrap_solver(fixed["official"]["ode4u"])
    return observer


def _finish(result, observer, scenario_id, start, dt, gate_scale, design, seed=None, constant_opening=None):
    diagnostics = result["diagnostics"]
    diagnostics.update(scenario_id=scenario_id, original_seed=seed, dt_s=float(dt),
                       gate_scale=float(gate_scale) if constant_opening is None else None,
                       schedule_kind=observer.fixed["manifest"]["schedule_kind"],
                       constant_opening=observer.fixed["manifest"]["constant_opening"],
                       design_source_label=design, all_rhs_observations=observer.summary())
    diagnostics["safety_violation_at_output_nodes"] = diagnostics.get("safety_violation")
    diagnostics["safety_violation_any_successful_rhs"] = (any(row["exceeded"] for row in observer.safety_successful.values())
                                                         if observer.successes else None)
    violation = diagnostics.get("first_violation")
    if isinstance(violation, dict):
        violation = dict(violation)
        head, opening = violation.get("water_depth_m"), violation.get("opening_height_m")
        if head is not None and opening is not None:
            violation["H_over_a"] = float(head/opening) if opening != 0 else None
        if "time_s" in violation and start is not None:
            violation["timestamp_utc"] = (pd.Timestamp(start)+pd.to_timedelta(violation["time_s"], unit="s")).isoformat()
        violation.update(scenario_id=scenario_id, case_name=result["case_name"])
        diagnostics["first_violation"] = violation
    reason = diagnostics.get("failure_reason")
    if result["valid"]:
        unsafe = diagnostics.get("safety_violation") or diagnostics["safety_violation_any_successful_rhs"]
        diagnostics["classification"] = "completed_valid_but_safety_exceeded" if unsafe else "completed_valid_safety_feasible"
        peak_index = int(np.argmax(result["downstream_stage"]))
        diagnostics.update(physical_river_peak_m=float(result["downstream_stage"][peak_index]),
                           physical_river_peak_time_s=float(result["time"][peak_index]))
    elif reason == "unsupported_gate_regime":
        diagnostics["classification"] = "unsupported_gate_hydraulics"
    elif reason == "hydrologic_geometry_applicability":
        diagnostics["classification"] = "other_physical_geometry_applicability"
    elif reason == "nonfreezing_liquid_water_domain":
        diagnostics["classification"] = "meteorological_liquid_water_applicability"
    else:
        diagnostics["classification"] = "numerical_or_state_failure"
    diagnostics["full_valid_trajectory_available"] = bool(result["valid"])
    diagnostics["full_integration_horizon_completed"] = bool(
        diagnostics.get("completed_until_s", 0) == observer.fixed["manifest"]["duration_s"])
    return result


def run_saved(case, weatherpath, dt=120, gate_scale=1, design="best", *, scenario_id=None, seed=None, constant_opening=None):
    """Run one saved hourly scenario, preserving the original hourly restarts.

    Input errors intentionally raise before simulation; the orchestrator must
    report them as input-validation failures rather than hydraulic exclusions.
    """
    module = adapter()
    canonical = module.CASE_ALIASES.get(case)
    if canonical is None:
        raise ValueError("Unknown case alias: " + str(case))
    fixed = load_fixed(dt, gate_scale, design, constant_opening)
    observer = _attach_observer(canonical, fixed)
    frame, _ = module.read_forcing(weatherpath)
    result = module.run_case(canonical, weatherpath, fixed)
    return _finish(result, observer, scenario_id or Path(weatherpath).stem,
                   frame["timestamp_utc"].iloc[0], dt, gate_scale, design, seed, constant_opening)


def run_deterministic(case="wetland", dt=120, gate_scale=1, design="best", *, check_regression=True, constant_opening=None):
    """Run the original September 2022 knots, optionally compare exact 120s output.

    A/A reference and B/B best are separate eight-state runs. No-wetland uses
    the same shared physical initial states; its wetland design is inapplicable.
    """
    module = adapter()
    canonical = module.CASE_ALIASES.get(case)
    if canonical is None:
        raise ValueError("Unknown case alias: " + str(case))
    fixed = load_fixed(dt, gate_scale, design, constant_opening)
    config = fixed["model_config"]
    params = module.wet.base_parameters(config)
    forcing = module.wet.fixed_forcing(config, params)
    time = params.time
    inputs = module.dry.inputs(config, fixed["fixed_gates"], time, forcing)
    observer = _attach_observer(canonical, fixed)
    diagnostics = dict(forcing_mapping="Original documented deterministic linear weather knots; no hourly remapping",
                       fixed_configuration=fixed["manifest"], event_start_utc=config["geography"]["event_start_utc"],
                       forcing_conservation=dict(
                           source_precipitation_m=float(np.trapezoid(forcing[1][0], forcing[0])),
                           mapped_precipitation_m=float(np.trapezoid(inputs[module.a2.U["P"]], time)),
                           source_solar_Wh_m2=float(np.trapezoid(forcing[1][2], forcing[0])/3600),
                           mapped_solar_Wh_m2=float(np.trapezoid(inputs[module.a2.U["R_s"]], time)/3600)))
    result = module._integrate(canonical, [(time, inputs)], fixed, diagnostics)
    result = _finish(result, observer, "deterministic_september_2022", config["geography"]["event_start_utc"],
                     dt, gate_scale, design, constant_opening=constant_opening)
    if result["valid"]:
        start = pd.Timestamp(config["geography"]["event_start_utc"])
        result["timestamp_utc"] = [(start+pd.to_timedelta(t, unit="s")).isoformat() for t in result["time"]]
    if check_regression:
        path = ROOT / "B/results/wetland_counterfactual/bypass_120s.npz" if canonical == "no_wetland" else _design_path(design)
        original_sha = module.sha256_file(path)
        original_schedule = gate_scale == 1 and constant_opening is None
        eligible = dt == 120 and original_schedule and result["valid"]
        record = dict(saved_trajectory=str(path.relative_to(ROOT)), saved_sha256=original_sha,
                      exact_original_grid_and_schedule=bool(dt == 120 and original_schedule),
                      pass_regression=None, storage_atol_m3=1e-7, output_atol=1e-9)
        if eligible:
            saved = module.load_trajectory(path)
            differences = {key: float(np.max(np.abs(result[key]-saved[key])))
                           for key in ("x", "dx", "y", "u") if key in saved}
            record.update(numerical_differences=differences,
                          pass_regression=bool(differences["x"] <= 1e-7 and
                                               all(v <= 1e-9 for key, v in differences.items() if key != "x")),
                          missing_saved_columns=[key for key in ("dx",) if key not in saved])
        else:
            record["reason_not_tested"] = "Only complete original-schedule 120 s runs are exact regression comparisons."
        record["saved_file_unchanged"] = original_sha == module.sha256_file(path)
        result["diagnostics"]["regression"] = record
    return result


def failure_record(result):
    """Compact CSV-ready first-failure row; no invalid partial physical peak."""
    diagnostics = result["diagnostics"]
    row = {name: diagnostics.get(name) for name in
           ("scenario_id", "original_seed", "dt_s", "gate_scale", "schedule_kind", "constant_opening", "design_source_label", "classification",
            "failure_stage", "failure_reason", "failure_message", "completed_until_s",
            "attempted_until_s", "physical_river_peak_m", "physical_river_peak_time_s")}
    row.update(case_name=result["case_name"], valid=bool(result["valid"]), solver_status=result["solver_status"])
    violation = diagnostics.get("first_violation") or {}
    row.update({name: violation.get(name) for name in
                ("branch", "time_s", "timestamp_utc", "upstream_stage_m", "tailwater_stage_m",
                 "water_depth_m", "opening_height_m", "normalized_opening", "H_over_a", "residual_m", "regime")})
    return json_clean(row)
