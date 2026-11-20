"""Matched 7.5 s refinement of the frozen nominal Stage D feedback pilot.

The only in-memory adaptation of the existing driver is one timestep whitelist
entry. The original file, hydraulic equations, policy, 120 s controller clock,
weather observations/scenarios, states and official ode4u implementation stay
unchanged. This is a numerical refinement, not another Monte Carlo sample.
"""
from __future__ import annotations

import argparse
import ast
from functools import lru_cache
import hashlib
import inspect
import json
from pathlib import Path
import sys
import textwrap
import time
import traceback

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from D.feedback_control import run_policy as frozen

COMPARISON = ROOT / 'D/feedback_control/results/fine_nominal'
DT_S = 7.5


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@lru_cache(maxsize=1)
def refined_runner():
    """Compile a private function; never mutate imported modules or their code."""
    source = textwrap.dedent(inspect.getsource(frozen.run_case))
    tree = ast.parse(source)
    guards = [node for node in ast.walk(tree)
              if isinstance(node, ast.Compare)
              and isinstance(node.left, ast.Name) and node.left.id == 'dt'
              and len(node.ops) == 1 and isinstance(node.ops[0], ast.NotIn)
              and len(node.comparators) == 1
              and isinstance(node.comparators[0], ast.Tuple)]
    if len(guards) != 1:
        raise RuntimeError('Expected exactly one timestep-membership guard; review the frozen driver before adapting')
    allowed = guards[0].comparators[0]
    if [ast.literal_eval(item) for item in allowed.elts] != [120, 60, 30, 15]:
        raise RuntimeError('Frozen timestep whitelist differs from the audited 120/60/30/15 s version')
    allowed.elts.append(ast.Constant(value=DT_S))
    ast.fix_missing_locations(tree)
    namespace = dict(frozen.run_case.__globals__)
    code = compile(tree, str(HERE / 'run_refinement.py') + ' [frozen run_case; whitelist-only adaptation]', 'exec')
    exec(code, namespace)
    adapted = namespace['run_case']
    metadata = dict(
        kind='In-memory AST adaptation: add 7.5 only to the timestep-membership whitelist',
        frozen_driver=str(Path(inspect.getsourcefile(frozen.run_case)).relative_to(ROOT)),
        frozen_driver_sha256=sha(inspect.getsourcefile(frozen.run_case)),
        frozen_run_case_source_sha256=hashlib.sha256(source.encode()).hexdigest(),
        adapted_ast_sha256=hashlib.sha256(ast.dump(tree, include_attributes=False).encode()).hexdigest(),
        refinement_driver_sha256=sha(__file__), timestep_s=DT_S,
        controller_period_s=frozen.CONTROL_PERIOD_S,
        no_integer_timestep_cast=True, no_original_module_mutation=True,
        changed_physics_or_policy=False)
    return adapted, metadata


def run_case(case, scenario_id=None, design='best'):
    """Run only the existing nominal hydraulic arm at the actual 7.5 s step."""
    adapted, metadata = refined_runner()
    result = adapted(case=case, scenario_id=scenario_id, design=design,
                     dt=DT_S, hydraulics='D', Cw_US=2.8)
    diagnostics = result['diagnostics']
    if diagnostics['dt_s'] != DT_S:
        raise RuntimeError('Refinement timestep was altered')
    if result['valid'] and not np.all(np.diff(result['time']) == DT_S):
        raise RuntimeError('Saved integration grid does not have exact 7.5 s intervals')
    diagnostics['refinement_adapter'] = metadata
    diagnostics['comparison_15s_directory'] = str(COMPARISON.relative_to(ROOT))
    return result


def check_inputs():
    """Validate saved provenance and floating clocks; perform no ODE calls."""
    _, metadata = refined_runner()
    comparison = json.loads((COMPARISON / 'execution_plan.json').read_text())
    weather = frozen.original_weather()
    if (comparison['hydraulics'] != 'D' or comparison['Cw_US'] != 2.8
            or comparison['control_period_s'] != 120
            or comparison['policy_sha256'] != sha(frozen.POLICY_PATH)):
        raise RuntimeError('The saved 15 s comparison is not the same nominal policy/physics')
    if len(weather) != 20 or sum(w['original_liquid_domain_admissible'] for w in weather) != 16:
        raise RuntimeError('Expected the original 16 warm and 4 cold scenario identifiers')
    if weather != comparison['weather']:
        raise RuntimeError('Saved scenario identities, seeds, forcing hashes or intervals differ from the 15 s comparison')
    module = frozen.historical.adapter()
    fixed = frozen.historical.load_fixed(dt=DT_S, gate_scale=1, design='best', constant_opening=None)
    if fixed['manifest']['dt_s'] != DT_S or 120 / DT_S != 16 or 3600 / DT_S != 480:
        raise RuntimeError('The original input adapter did not preserve the fractional timestep')
    for item in weather:
        if sha(item['forcing_file']) != item['forcing_sha256']:
            raise RuntimeError('Saved weather bytes changed: ' + item['scenario_id'])
        frame, _ = module.read_forcing(item['forcing_file'])
        conservation = module.forcing_conservation(frame, fixed)
        if not conservation['pass_conservation']:
            raise RuntimeError('Original precipitation/solar mapping failed at 7.5 s')
        first_time, _ = next(module.hourly_input_segments(frame, fixed))
        if len(first_time) != 481 or not np.all(np.diff(first_time) == DT_S):
            raise RuntimeError('Hourly input clock is not the exact requested floating grid')
    return weather, comparison, metadata


def verify_15s_regression(destination):
    """Run the whitelist-adapted function at 15 s against its saved control."""
    adapted, metadata = refined_runner()
    reference_path = COMPARISON / 'deterministic_september_2022__wetland__best__15s.npz'
    start = time.monotonic()
    result = adapted(case='wetland', scenario_id=None, design='best',
                     dt=15, hydraulics='D', Cw_US=2.8)
    checks = {}
    if result['valid']:
        with np.load(reference_path, allow_pickle=False) as reference:
            for key in ('time', 'x', 'dx', 'y', 'u', 'balance',
                        'downstream_stage', 'downstream_discharge'):
                previous = reference[key]
                current = np.asarray(result[key])
                same_shape = previous.shape == current.shape
                checks[key] = dict(same_shape=same_shape,
                    exactly_equal=bool(same_shape and np.array_equal(previous, current)),
                    maximum_absolute_difference=None if not same_shape else float(np.max(np.abs(previous-current))))
    report = dict(passed=bool(result['valid'] and checks and all(v['exactly_equal'] for v in checks.values())),
        scope='Actual official-ode4u deterministic wetland historical-B-selected design at unchanged 15 s, using the in-memory whitelist-adapted driver',
        reference_file=str(reference_path.relative_to(ROOT)), reference_sha256=sha(reference_path),
        actual_run_classification=result['diagnostics']['classification'],
        compared_arrays=checks, elapsed_s=time.monotonic()-start,
        refinement_adapter=metadata, run_diagnostics=result['diagnostics'])
    frozen.write_json(destination, report)
    if not report['passed']:
        raise RuntimeError('The adapted 15 s control did not reproduce the saved arrays; refinement requires review')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only', action='store_true', help='Read-only provenance and clock checks; no ODE calls')
    args = parser.parse_args()
    weather, comparison, metadata = check_inputs()
    if args.check_only:
        print(json.dumps(dict(passed=True, timestep_s=DT_S, controller_period_s=120,
            original_weather_ids=len(weather), warm_ids=16, cold_ids=4,
            planned_ode_cases=35, no_ODE_calls=True, adaptation=metadata), indent=2))
        return
    from D.d2_analysis import preservation
    guard = preservation.verify()
    if not guard['passed']:
        raise RuntimeError('Protected files changed: ' + json.dumps(guard))
    output = HERE / 'results/fine_7p5'
    if output.exists():
        raise FileExistsError(output)
    plans = []
    for plan in comparison['plans']:
        if float(plan['dt']) != 15:
            raise RuntimeError('Comparison contains a non-15 s run')
        plans.append(dict(plan, dt=DT_S))
    output.mkdir(parents=True)
    print('START exact 15 s deterministic regression control', flush=True)
    regression = verify_15s_regression(output / 'regression_15s.json')
    print('DONE exact 15 s regression', regression['passed'], round(regression['elapsed_s'], 2), flush=True)
    frozen.write_json(output / 'execution_plan.json', dict(plans=plans, weather=weather,
        policy=comparison['policy'], policy_sha256=comparison['policy_sha256'],
        hydraulics='D', Cw_US=2.8, control_period_s=120, refinement_adapter=metadata,
        comparison_15s_directory=str(COMPARISON.relative_to(ROOT)),
        scope='Same original 20 scenarios; 16 warm pairs and 4 cold IDs accounted without ODE, plus three deterministic controls; no new weather draws'))
    rows = []
    for plan in plans:
        sid = plan['scenario_id'] or 'deterministic_september_2022'
        name = f'{sid}__{plan["case"]}__{plan["design"]}__7p5s'
        print('START', name, flush=True)
        start = time.monotonic()
        try:
            result = run_case(case=plan['case'], scenario_id=plan['scenario_id'], design=plan['design'])
            row = frozen.save_run(result, output, name)
        except Exception as exc:
            error = dict(scenario_id=sid, case_name=plan['case'], design=plan['design'], dt_s=DT_S,
                hydraulic_arm='D', classification='unexpected_execution_error', valid=False,
                exception_type=type(exc).__name__, failure_message=str(exc), exception_traceback=traceback.format_exc())
            frozen.write_json(output / (name + '.error.json'), error)
            row = {key:value for key,value in error.items() if key != 'exception_traceback'}
        row['elapsed_s'] = time.monotonic() - start
        rows.append(row)
        pd.DataFrame(rows).to_csv(output / 'scenario_diagnostics.csv', index=False)
        print('DONE', name, row.get('classification'), round(row['elapsed_s'], 2), flush=True)
    counts = pd.DataFrame(rows).groupby(['case_name', 'classification'], dropna=False).size().to_dict()
    frozen.write_json(output / 'counts.json', {f'{key[0]} / {key[1]}':int(value) for key,value in counts.items()})
    guard = preservation.verify()
    frozen.write_json(output / 'preservation_check.json', guard)
    if not guard['passed']:
        raise RuntimeError('Protected files changed during execution')


if __name__ == '__main__':
    main()
