"""Wetland design search: fixed gates, two continuous covers, four class pairs.

G5 remains an unchanged model/RK-stage validation guard, never an OPT residual.
Run --run, --compare, --report; completed numerical stages are not overwritten.
"""
from pathlib import Path
from collections import Counter,OrderedDict
from dataclasses import dataclass, field
import argparse,hashlib,itertools,json,time
import numpy as np
from B import optimize_wetland as old
from B import no_spillway_model as model
from B.official_opt import run_opt,load_official
from B.multistart_check import save_trajectory
from B.no_spillway_study import load_trajectory

ROOT=old.ROOT
OUT=ROOT/'results/wetland_only'
SPEC=ROOT/'wetland_only_config.json'
G5_LABELS=frozenset({'gate1_head_opening','gate2_head_opening'})
OPT_CONSTRAINT_COUNT=18


def save(name,data):old.save_json(OUT/name,data)


def load_study():
    spec=old.load_config(SPEC)
    for path_key,hash_key in [('model_config','model_config_sha256'),('fixed_gate_source','fixed_gate_source_sha256')]:
        if hashlib.sha256(old.resolve_path(spec[path_key]).read_bytes()).hexdigest()!=spec[hash_key]:
            raise RuntimeError('Frozen wetland-only input changed: '+path_key)
    config=old.load_config(ROOT/spec['model_config'])
    if model.preflight(config):raise RuntimeError('Unapproved model configuration')
    saved=load_trajectory(ROOT/spec['fixed_gate_source'])
    fixed_gates=saved['theta'][:-2].copy();fixed_gates.setflags(write=False)
    return config,spec,fixed_gates


def assemble_design(covers,fixed_gates):
    covers=np.asarray(covers,dtype=float)
    if covers.shape!=(2,) or not np.all(np.isfinite(covers)):
        raise ValueError('Only two finite continuous wetland cover decisions are accepted.')
    if np.any(covers<0) or np.any(covers>1):raise old.TrajectoryFailure('Cover outside full native [0,1] bounds; not clipped.')
    gates=np.asarray(fixed_gates,dtype=float)
    if gates.shape!=(24,):raise ValueError('Expected the unchanged 24 fixed gate knots.')
    return np.r_[gates,covers]


def optimization_view(run):
    """Separate active residuals, fixed-schedule applicability and native bounds."""
    original=run['audit'][:model.CONSTRAINT_COUNT]
    active=[r for r in original if r['label'] not in G5_LABELS]
    applicability=[r for r in original if r['label'] in G5_LABELS]
    if len(active)!=OPT_CONSTRAINT_COUNT or len(applicability)!=2 or active[-1]['label']!='integration_status':raise RuntimeError('Unexpected residual schema')
    native_covers=[r for r in run['audit'][model.CONSTRAINT_COUNT:] if r['label'].startswith('cover')]
    fixed_bounds=[r for r in run['audit'][model.CONSTRAINT_COUNT:] if r['label'].startswith('gate')]
    if len(native_covers)!=4 or len(fixed_bounds)!=4:raise RuntimeError('Unexpected bound audit')
    # Stage validation has already guarded hydraulics. These are separate
    # model-applicability records and are deliberately absent from OPT's vector.
    if not all(r['strict_pass'] for r in applicability) or not all(r['pass_at_tolerance'] for r in fixed_bounds):
        raise old.TrajectoryFailure('Fixed gate schedule failed model applicability validation.')
    feasible=all(r['pass_at_tolerance'] for r in active+native_covers)
    return dict(residuals=np.array([r['scaled_residual'] for r in active]),active_audit=active,
                applicability_audit=applicability,fixed_gate_bound_audit=fixed_bounds,native_cover_audit=native_covers,
                feasible=bool(feasible))


def invalid_response(config):
    """An invalid integration has no river peak; no G5 residual is fabricated."""
    residuals=np.zeros(OPT_CONSTRAINT_COUNT)
    residuals[-1]=config['numerics']['failure_residual']  # retained integration_status
    return float(config['numerics']['failure_objective']),residuals


def simulate_covers(covers,pair,config,fixed_gates,*,forcing=None,refine=1):
    theta=assemble_design(covers,fixed_gates)
    run=model.simulate(theta,pair,config,forcing=forcing,refine=refine)
    view=optimization_view(run)
    assert np.array_equal(run['theta'][:-2],fixed_gates)
    return run,view



@dataclass
class AnalysisConstants:
    """Fixed inputs for one category pair; covers are the only decisions.

    last_evaluation is replaced on every call and exposes physical results or an
    invalid-integration diagnostic. It is never used as an input to simulation.
    """
    config: dict
    fixed_gates: np.ndarray
    pair: tuple
    forcing: object = None
    refine: int = 1
    last_evaluation: dict | None = field(default=None, init=False)

    def __post_init__(self):
        self.fixed_gates = np.array(self.fixed_gates, dtype=float, copy=True)
        self.fixed_gates.setflags(write=False)
        self.pair = tuple(self.pair)
        if len(self.pair) != 2 or any(
                label not in self.config['categories']['admissible'][i]
                for i, label in enumerate(self.pair)):
            raise ValueError('Select one approved response class for each wetland.')
        if self.refine not in (1, 2, 4):
            raise ValueError('Use the existing 120, 60 or 30 second study grid.')


def make_analysis_constants(pair=None, *, refine=1):
    """Load the frozen study without starting OPT or writing any output."""
    config, spec, gates = load_study()
    return AnalysisConstants(config, gates,
        tuple(spec['reference_categories'] if pair is None else pair),
        model.fixed_forcing(config, model.base_parameters(config)), refine)


def analysis(v, c):
    """Return (f, g) for v=[rho1,rho2], using the official dynamic model.

    f=(physical river peak - fixed origin)/fixed positive height scale; g is the
    18-element dimensionless inequality vector, with g<=0 physically satisfied
    and the existing declared numerical tolerance used for feasibility. G5 is
    an RHS/RK model-applicability check, not an entry in this optimizer vector.

    Expected invalid trajectories return the unchanged finite failure response;
    c.last_evaluation then has no physical peak or complete trajectory. A valid
    unsafe trajectory retains its physical objective and positive constraint
    residuals. Shape/configuration/programming errors are not hidden.
    """
    c.last_evaluation = None
    try:
        run, view = simulate_covers(v, c.pair, c.config, c.fixed_gates,
                                   forcing=c.forcing, refine=c.refine)
        f = (run['metrics']['peak_stage_m'] - c.config['numerics']['objective_origin_m']) / c.config['numerics']['height_scale_m']
        g = view['residuals']
        error = None
        diagnostic = None
    except (old.TrajectoryFailure, FloatingPointError) as exc:
        f, g = invalid_response(c.config)
        run = view = None
        error = type(exc).__name__ + ': ' + str(exc)
        diagnostic = getattr(exc, 'diagnostic', None)
    g = np.asarray(g, dtype=float)
    if not np.isfinite(f) or g.shape != (OPT_CONSTRAINT_COUNT,) or not np.all(np.isfinite(g)):
        raise ValueError('Analysis must return one finite objective and 18 finite residuals.')
    c.last_evaluation = dict(run=run, view=view, failure=error, diagnostic=diagnostic,
        physical_peak_m=None if run is None else run['metrics']['peak_stage_m'],
        valid=run is not None, feasible=bool(view is not None and view['feasible']))
    return float(f), g.copy()


def run_study():
    if (OUT/'batch_status.json').exists():raise RuntimeError('Completed wetland-only batch exists; preserving it.')
    config,spec,gates=load_study();forcing=model.fixed_forcing(config,model.base_parameters(config))
    official=load_official()
    save('study_plan.json',dict(specification=spec,physical_model_configuration=config,
        fixed_gate_knots=gates,continuous_design=['rho1','rho2'],continuous_dimension=2,
        optimizer_constraints=18,model_applicability_only=sorted(G5_LABELS),native_continuous_bounds=[[0,1],[0,1]],
        official_source={k:v for k,v in official.items() if k not in ['nms','ode4u','opt_hyp']}))
    reference,view=simulate_covers(spec['reference_covers'],spec['reference_categories'],config,gates,forcing=forcing)
    if not view['feasible']:raise RuntimeError('The fixed-schedule reference is not feasible; no optimization started.')
    save_trajectory(OUT/'reference.npz',reference)
    save('fixed_gate_validation.json',dict(status='passed_before_search',gate_knots=gates,metrics=reference['metrics'],
        applicability_audit=view['applicability_audit'],fixed_gate_bound_audit=view['fixed_gate_bound_audit'],
        gate_tailwater_separation_m=np.array(config['gates']['gate_sill_z'])-np.array(config['safety']['wetland_m']),
        scope='G5 guarded at every RHS/RK stage of every simulation; reporting margins are supplementary. No G5 residual sent to OPT.'))
    records=[];accepted=[]
    for pair in itertools.product(*config['categories']['admissible']):
        for i,covers in enumerate(spec['cover_starts'],1):
            row=dict(category_pair=list(pair),start_id=f'cover_start_{i}',requested_covers=covers)
            try:
                result,v=simulate_covers(covers,pair,config,gates,forcing=forcing)
                row.update(feasible=v['feasible'],status='feasible' if v['feasible'] else 'valid_constraint_violation',
                           metrics=result['metrics'],**{k:v[k] for k in ['active_audit','applicability_audit','native_cover_audit']})
                if v['feasible']:accepted.append(dict(category_pair=list(pair),start_id=row['start_id'],covers=covers))
            except (old.TrajectoryFailure,FloatingPointError) as exc:
                row.update(feasible=False,status='invalid_model_evaluation',failure=str(exc),diagnostic=getattr(exc,'diagnostic',None))
            records.append(row);print('Screen',pair,covers,row['status'],flush=True)
    save('start_screening.json',dict(records=records,accepted=accepted))
    started=time.monotonic();total_odes=1+len(records);summaries=[];global_best=None;aborted=False
    for start in accepted:
        pair=start['category_pair'];initial=np.array(start['covers']);label=start['start_id'];slug='_'.join(pair)+'_'+label
        log=[];cache=OrderedDict();best=None;first=[];odes=0;before=time.monotonic();solver={}
        constants = AnalysisConstants(config, gates, tuple(pair), forcing)
        def callback(covers,forwarded_constants=None):
            nonlocal best,odes
            if time.monotonic()-started>spec['maximum_batch_seconds'] or total_odes+odes>=spec['maximum_study_ode_attempts']-20:
                raise old.RunBudgetExceeded('Declared wetland-only batch budget exhausted.')
            if not first:
                first.extend(np.asarray(covers).tolist())
                if not np.allclose(covers,initial,rtol=0,atol=1e-15):raise RuntimeError('Official initialization moved cover start.')
            key=np.asarray(covers,dtype=float).tobytes();cached=key in cache
            if cached:f,g,result,v,error,diagnostic=cache[key]
            else:
                odes+=1
                f,g=analysis(covers,constants)
                evaluation=constants.last_evaluation
                result=evaluation['run'];v=evaluation['view']
                error=evaluation['failure'];diagnostic=evaluation['diagnostic']
                cache[key]=(f,g,result,v,error,diagnostic)
                if len(cache)>8:cache.popitem(last=False)
            feasible=result is not None and v['feasible']
            if feasible and (best is None or result['metrics']['peak_stage_m']<best['metrics']['peak_stage_m']):best=result
            log.append(dict(evaluation=len(log)+1,covers=np.asarray(covers),full_design=assemble_design(covers,gates),cached=cached,
                raw_peak_stage_m=None if result is None else result['metrics']['peak_stage_m'],conditioned_objective=float(f),
                optimizer_residuals=g,optimizer_residual_kind='computed_physical_residuals' if result is not None else 'synthetic_invalid_integration_response_not_physical_residuals',feasible=feasible,max_scaled_violation=float(max(0,max(g))),
                best_feasible_peak_m=None if best is None else best['metrics']['peak_stage_m'],failure=error or '',diagnostic=diagnostic,
                model_applicability_passed=result is not None,
                G5_validation=None if v is None else v['applicability_audit']))
            return float(f),g.copy()
        try:
            solver=run_opt(callback,initial,np.zeros(2),np.ones(2),spec['settings']);termination=solver['status']
        except old.RunBudgetExceeded as exc:termination='batch_budget: '+str(exc);aborted=True
        total_odes+=odes
        if best is not None:
            save_trajectory(OUT/(slug+'_best.npz'),best)
            if global_best is None or best['metrics']['peak_stage_m']<global_best['metrics']['peak_stage_m']:global_best=best
        bv=None if best is None else optimization_view(best)
        record=dict(category_pair=pair,start_id=label,requested_initial_covers=initial,first_optimizer_covers=first,
            fixed_gate_knots=gates,first_design_feasible=bool(log and log[0]['feasible']),optimizer_result=solver,
            physical_callback_log=log,termination=termination,callbacks=len(log),cached_callbacks=sum(r['cached'] for r in log),
            actual_ode_attempts=odes,seconds=time.monotonic()-before,rejection_counts=dict(Counter(r['failure'].split(':',1)[0] for r in log if r['failure'])),
            best_covers=None if best is None else best['theta'][-2:],best_design=None if best is None else best['theta'],
            best_metrics=None if best is None else best['metrics'],best_audit=None if best is None else best['audit'],
            best_active_audit=None if bv is None else bv['active_audit'],best_applicability_audit=None if bv is None else bv['applicability_audit'])
        save(slug+'_run.json',record)
        summaries.append({k:v for k,v in record.items() if k not in ['optimizer_result','physical_callback_log','best_audit','best_active_audit','best_applicability_audit']})
        print('OPT',pair,label,termination,'peak',None if best is None else best['metrics']['peak_stage_m'],
              'covers',None if best is None else best['theta'][-2:],'ODE',odes,flush=True)
        if aborted:break
    if global_best is not None:save_trajectory(OUT/'batch_best.npz',global_best)
    save('batch_status.json',dict(status='bounded_batch_stopped' if aborted else 'bounded_batch_finished',runs=summaries,
        seconds=time.monotonic()-started,actual_ode_attempts_including_reference_and_screening=total_odes,
        global_best_pair=None if global_best is None else global_best['pair'],global_best_design=None if global_best is None else global_best['theta']))


def compare():
    if (OUT/'fixed_design_timesteps.json').exists():raise RuntimeError('Fixed-design checks exist; preserving results.')
    config,spec,gates=load_study();plan=old.load_config(OUT/'study_plan.json')
    if plan['specification']!=spec or plan['physical_model_configuration']!=config:raise RuntimeError('Case changed since search.')
    forcing=model.fixed_forcing(config,model.base_parameters(config));designs=[('fixed_reference',load_trajectory(OUT/'reference.npz'))]
    if (OUT/'batch_best.npz').exists():designs.append(('wetland_best',load_trajectory(OUT/'batch_best.npz')))
    # Resolve the very small category gaps without reoptimizing or changing inputs.
    by_pair={}
    for path in OUT.glob('*_run.json'):
        rec=old.load_config(path)
        if rec['best_metrics'] is None:continue
        pair=tuple(rec['category_pair'])
        if pair not in by_pair or rec['best_metrics']['peak_stage_m']<by_pair[pair][1]['best_metrics']['peak_stage_m']:by_pair[pair]=(path,rec)
    category_labels=[]
    for pair,(path,rec) in sorted(by_pair.items()):
        existing=next((label for label,r in designs if r['pair']==list(pair) and np.array_equal(r['theta'],rec['best_design'])),None)
        if existing is None:
            label='category_best_'+'_'.join(pair)
            designs.append((label,load_trajectory(path.with_name(path.name.replace('_run.json','_best.npz')))))
        else:label=existing
        category_labels.append(label)
    rows=[]
    def evaluate(label,r,dt):
        try:
            run,view=simulate_covers(r['theta'][-2:],r['pair'],config,gates,forcing=forcing,refine=int(120/dt))
            filename=f'{label}_{dt}s.npz';save_trajectory(OUT/filename,run)
            row=dict(label=label,dt_s=dt,design=run['theta'],pair=run['pair'],covers=run['theta'][-2:],status='valid',feasible=view['feasible'],
                metrics=run['metrics'],audit=run['audit'],active_audit=view['active_audit'],applicability_audit=view['applicability_audit'],
                trajectory_file=filename,rain_depth_m=float(np.trapezoid(run['u'][old.a2.U['P']],run['time'])),
                solar_energy_Wh_m2=float(np.trapezoid(run['u'][old.a2.U['R_s']],run['time'])/3600))
        except (old.TrajectoryFailure,FloatingPointError) as exc:
            row=dict(label=label,dt_s=dt,design=r['theta'],pair=r['pair'],status='invalid',feasible=False,failure=str(exc),diagnostic=getattr(exc,'diagnostic',None))
        rows.append(row);print('Resolution',label,dt,row['status'],row['feasible'],flush=True)
    for label,r in designs:
        for dt in [120,60]:evaluate(label,r,dt)
    comparisons=[];category_comparisons=[];need30=False
    def comparison(coarse,fine):
        by={(r['label'],r['dt_s']):r for r in rows}
        group=[by.get((label,dt)) for label in ['fixed_reference','wetland_best'] for dt in [coarse,fine]]
        if any(r is None or r['status']!='valid' for r in group):return dict(coarse_dt_s=coarse,fine_dt_s=fine,resolved=False)
        refc,reff,bestc,bestf=group
        reduction_coarse=refc['metrics']['peak_stage_m']-bestc['metrics']['peak_stage_m']
        reduction_fine=reff['metrics']['peak_stage_m']-bestf['metrics']['peak_stage_m']
        change=abs(refc['metrics']['peak_stage_m']-reff['metrics']['peak_stage_m'])+abs(bestc['metrics']['peak_stage_m']-bestf['metrics']['peak_stage_m'])
        return dict(coarse_dt_s=coarse,fine_dt_s=fine,reduction_coarse_m=reduction_coarse,reduction_fine_m=reduction_fine,
            sum_grid_changes_m=change,resolved=bool(all(r['feasible'] for r in group) and np.sign(reduction_coarse)==np.sign(reduction_fine) and change<max(.1*abs(reduction_fine),1e-10)))
    def category_check(coarse,fine):
        by={(r['label'],r['dt_s']):r for r in rows};records=[]
        for first_label,second_label in itertools.combinations(category_labels,2):
            g=[by.get((label,dt)) for label in [first_label,second_label] for dt in [coarse,fine]]
            row=dict(first_label=first_label,second_label=second_label,coarse_dt_s=coarse,fine_dt_s=fine,resolved=False)
            if all(r is not None and r['status']=='valid' for r in g):
                ac,af,bc,bf=g
                dc=bc['metrics']['peak_stage_m']-ac['metrics']['peak_stage_m'];df=bf['metrics']['peak_stage_m']-af['metrics']['peak_stage_m']
                change=abs(ac['metrics']['peak_stage_m']-af['metrics']['peak_stage_m'])+abs(bc['metrics']['peak_stage_m']-bf['metrics']['peak_stage_m'])
                row.update(peak_difference_coarse_m=dc,peak_difference_fine_m=df,sum_grid_changes_m=change,
                    resolved=bool(all(r['feasible'] for r in g) and np.sign(dc)==np.sign(df) and change<max(.1*abs(df),1e-10)))
            records.append(row)
        return records
    if len(designs)>1:
        first=comparison(120,60);comparisons.append(first)
        category_comparisons.extend(category_check(120,60))
        need30=not first['resolved'] or any(not r['resolved'] for r in category_comparisons)
        if need30:
            for label,r in designs:evaluate(label,r,30)
            comparisons.append(comparison(60,30))
            category_comparisons.extend(category_check(60,30))
    save('fixed_design_timesteps.json',dict(records=rows,comparisons=comparisons,category_comparisons=category_comparisons,thirty_second_required=need30))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',action='store_true');p.add_argument('--compare',action='store_true');p.add_argument('--report',action='store_true');args=p.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    if args.run:run_study()
    if args.compare:compare()
    if args.report:
        from B.wetland_only_reporting import export_study
        export_study(OUT)

if __name__=='__main__':main()
