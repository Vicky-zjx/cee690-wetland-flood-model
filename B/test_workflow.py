"""Unit checks of B.2 workflow contracts; no project ODE or optimization run.

All numerical safety/trajectory fixtures below are artificial TEST-ONLY arrays.
They are not proposed thresholds, approved parameters, or project evidence.
Run: python3 -B -m unittest B.test_workflow -v
"""
from __future__ import annotations

import copy
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from B import optimize_wetland as workflow
from B.gate_hydraulics import UnsupportedGateRegime
from B.official_opt import DEFAULT_SETTINGS


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.config = workflow.load_config()
        # Approval tests use an explicit unapproved fixture, independent of the
        # user's approved live exploratory configuration.
        self.config['mode'] = None
        for group in ('geography', 'safety', 'categories', 'gates', 'fixed_model',
                      'forcing', 'initial_conditions'):
            self.config[group]['approved'] = False
        self.config['gates']['limited_regimes_approved'] = False
        for override in self.config['parameter_overrides'].values():
            override['approved'] = False
        self.params = workflow.a2.c()
        self.theta = workflow.reference_decision(self.config)

    def artificial_audit_arrays(self):
        """Pure array fixture; intentionally includes violations, never simulated."""
        cfg = copy.deepcopy(self.config)
        # This artificial trajectory fixture is independent of unapproved event proposals.
        cfg['parameter_overrides'] = {}
        cfg['forcing']['kind'] = 'a2_synthetic_fixed_samples'
        cfg['safety'].update(dam_m=[10.,20.], wetland_m=[30.,40.],
                             source_ids=['TEST_ONLY_d1','TEST_ONLY_d2','TEST_ONLY_w1','TEST_ONLY_w2'])
        cfg['numerics'].update(height_scale_m=2., strict_margin_m=.02,
                               volume_scales_m3=[2.,4.,10.,20.,50.,60.,70.,80.])
        t = np.array([0.,60.,120.])
        x = np.repeat(self.params.x_init[:,None],3,axis=1)
        x[0] = [10.,-2.,5.]
        x[1] = [self.params.C_m-1.,self.params.C_m+8.,self.params.C_m+2.]
        y = np.zeros((len(workflow.a2.OUTPUT_NAMES),3))
        for name,values in {'H_d1':[9.,11.,10.], 'H_d2':[19.,18.,21.],
                            'H_w1':[29.,32.,30.], 'H_w2':[39.,39.5,39.8]}.items():
            y[workflow.a2.Y[name]] = values
        y[workflow.a2.Y['H_g']] = self.params.H_g_valid_max-.02+np.array([-.2,0.,.4])
        y[workflow.a2.Y['H_r']] = self.params.H_r_valid_max-.02+np.array([-.1,0.,.6])
        u = np.zeros((len(workflow.a2.INPUT_NAMES),3))
        u[:2] = .5
        return cfg,t,x,y,u

    @staticmethod
    def mock_result(peak=98., residual=-.1):
        """Minimal deterministic simulator contract; not physical model results."""
        return {'metrics':{'peak_stage_m':peak}, 'residuals':np.full(workflow.CONSTRAINT_COUNT,residual),
                'audit':[{'pass_at_tolerance':residual <= 1e-6}]}

    def evaluator(self):
        budget = {'evaluations':0,'start':workflow.clock.monotonic()}
        return workflow.Evaluator(self.config,('type_a','type_a'),None,budget)

    def test_project_preflight_exposes_exactly_six_unapproved_groups(self):
        untouched = copy.deepcopy(self.config)
        gaps = workflow.preflight(self.config)
        self.assertEqual([g['id'] for g in gaps], ['G1','G2','G3','G4','G5','G6'])
        self.assertTrue(all(g['status']=='blocked' for g in gaps))
        self.assertEqual(self.config, untouched)
        self.assertTrue(all(v is not None for v in self.config['safety']['dam_m']+self.config['safety']['wetland_m']))
        self.assertFalse(self.config['safety']['approved'])
        self.assertTrue(all(not record['approved'] for record in self.config['parameter_overrides'].values()))
        self.assertIsNone(self.config['controls']['rate_limits'])
        self.assertIsNone(self.config['reporting']['flood_threshold_m'])

    def test_approval_flags_alone_do_not_fill_missing_safety(self):
        cfg = copy.deepcopy(self.config)
        cfg['safety']['dam_m'] = [None,None]
        cfg['safety']['wetland_m'] = [None,None]
        cfg['safety']['approved'] = True
        cfg['safety']['source_ids'] = ['TEST_ONLY']*4
        gaps = workflow.preflight(cfg)
        self.assertIn('G2',[g['id'] for g in gaps])
        self.assertEqual(cfg['safety']['dam_m'], [None,None])

    def test_native_twenty_six_decisions_represent_a2_schedule(self):
        self.assertEqual(self.theta.shape,(26,))
        g1,g2,covers = workflow.decode(self.theta,self.config)
        np.testing.assert_array_equal(g1,self.params.gate1_values)
        np.testing.assert_array_equal(g2,self.params.gate2_values)
        np.testing.assert_array_equal(covers,self.params.rho)
        original = workflow.a2.make_inputs(self.params)
        for i,values in enumerate((g1,g2)):
            hours = self.config['controls'][f'gate{i+1}_hours']
            self.assertEqual(hours,getattr(self.params,f'gate{i+1}_hours'))
            reconstructed = np.interp(self.params.time,np.array(hours)*3600.,values)
            # Seconds-vs-hours interpolation differs only by roundoff, not schedule.
            np.testing.assert_allclose(reconstructed,original[i],rtol=0,atol=1e-15)
            np.testing.assert_array_equal(np.interp(np.array(hours)*3600.,np.array(hours)*3600.,values),values)

    def test_decode_rejects_instead_of_clipping_decisions(self):
        for bad in (-.001,1.001):
            theta = self.theta.copy();theta[0]=bad
            with self.assertRaises(workflow.TrajectoryFailure):
                workflow.decode(theta,self.config)
            self.assertEqual(theta[0],bad)
        with self.assertRaises(ValueError):
            workflow.decode(self.theta[:-1],self.config)

    def test_all_zero_and_all_one_knots_bound_interpolated_histories(self):
        t = np.linspace(0,self.params.duration,20001)
        for value in (0.,1.):
            theta = np.full(26,value)
            g1,g2,covers = workflow.decode(theta,self.config)
            for i,values in enumerate((g1,g2)):
                actual = np.interp(t,np.array(self.config['controls'][f'gate{i+1}_hours'])*3600.,values)
                np.testing.assert_array_equal(actual,np.full_like(t,value))
            np.testing.assert_array_equal(covers,[value,value])

    def test_constraints_include_two_explicit_gate_applicability_residuals(self):
        cfg,t,x,y,u = self.artificial_audit_arrays()
        residuals,rows = workflow.constraint_audit(t,x,y,u,self.theta,self.params,cfg)
        labels = [r['label'] for r in rows[:18]]
        self.assertEqual(residuals.shape,(workflow.CONSTRAINT_COUNT,))
        self.assertEqual(len(rows),workflow.CONSTRAINT_COUNT+8)
        self.assertEqual(labels[:8],['nonnegative_'+s for s in workflow.a2.STATE_NAMES])
        self.assertEqual(labels[8:13],['soil_capacity','H_d1_safe','H_d2_safe','H_w1_safe','H_w2_safe'])
        self.assertEqual(labels[13:],['groundwater_geometry','groundwater_soil','groundwater_wetland_bed','river_geometry','integration_status'])
        self.assertEqual(len(set(r['label'] for r in rows)),workflow.CONSTRAINT_COUNT+8)
        self.assertEqual(workflow.CONSTRAINT_COUNT,20)
        self.assertTrue(all('applicability' in r['kind'] for r in rows[18:20]))
        self.assertNotIn('H_flood',labels)
        self.assertNotIn('terminal_storage',labels)
        self.assertNotIn('gate_rate',labels)

    def test_storage_soil_and_safety_signs_scaling_and_worst_time(self):
        cfg,t,x,y,u = self.artificial_audit_arrays()
        residuals,rows = workflow.constraint_audit(t,x,y,u,self.theta,self.params,cfg)
        by_name = {r['label']:r for r in rows}
        self.assertEqual(residuals[0],1.)  # -min(-2)/2
        self.assertEqual(residuals[8],2.)  # soil excess8 / volume scale4
        self.assertEqual(by_name['nonnegative_S_u']['worst_time_s'],60.)
        self.assertEqual(by_name['nonnegative_S_u']['margin'],-2.)
        self.assertEqual(by_name['H_d1_safe']['scaled_residual'],.5)
        self.assertEqual(by_name['H_d2_safe']['worst_time_s'],120.)
        self.assertEqual(by_name['H_w1_safe']['scaled_residual'],1.)
        self.assertAlmostEqual(by_name['H_w2_safe']['margin'],.2)
        self.assertTrue(by_name['H_w2_safe']['strict_pass'])
        self.assertFalse(by_name['H_d1_safe']['pass_at_tolerance'])
        self.assertEqual(by_name['H_d1_safe']['source_id'],'TEST_ONLY_d1')

    def test_applicability_margin_is_explicit_and_separate_from_safety(self):
        cfg,t,x,y,u = self.artificial_audit_arrays()
        _,rows = workflow.constraint_audit(t,x,y,u,self.theta,self.params,cfg)
        by_name={r['label']:r for r in rows}
        self.assertAlmostEqual(by_name['groundwater_geometry']['scaled_residual'],.2)
        self.assertAlmostEqual(by_name['river_geometry']['scaled_residual'],.3)
        self.assertIn('numerical margin only',by_name['river_geometry']['kind'])
        self.assertEqual(by_name['H_d1_safe']['threshold'],10.)
        self.assertNotEqual(by_name['H_d1_safe']['threshold'],self.params.z_dsp[0])
        cfg['numerics']['strict_margin_m']=cfg['numerics']['feasibility_tol']*cfg['numerics']['height_scale_m']
        with self.assertRaisesRegex(workflow.PreflightError,'margin'):
            workflow.preflight(cfg)

    def test_low_head_guard_reports_first_invalid_hydraulic_state(self):
        # Local synthetic state chosen to violate only the approved gate law's
        # applicability domain. No completed trajectory or physical peak exists.
        cfg=copy.deepcopy(self.config)
        cfg['parameter_overrides']={}
        params=workflow.base_parameters(cfg)
        gates=workflow.gate_parameters(cfg)
        state=params.x_init.copy()
        state[4]=params.A_d[0]*(gates[0].gate_sill_z+.1-params.z_d[0])
        inputs=workflow.a2.make_inputs(params)[:,0].copy()
        inputs[0]=.5
        with self.assertRaises(workflow.HydraulicApplicabilityFailure) as raised:
            workflow.make_b2_rhs(gates)(123.,state,inputs,params)
        diagnostic=raised.exception.diagnostic
        self.assertEqual(diagnostic['branch'],1)
        self.assertEqual(diagnostic['time_s'],123.)
        self.assertAlmostEqual(diagnostic['water_depth_m'],.1)
        self.assertAlmostEqual(diagnostic['opening_height_m'],.15)
        self.assertAlmostEqual(diagnostic['residual_m'],1.25*.15-.1)
        self.assertGreater(diagnostic['residual_m'],0.)

    def test_gate_applicability_has_physical_residual_and_inactive_closed_case(self):
        cfg,t,x,y,u=self.artificial_audit_arrays()
        y[workflow.a2.Y['H_d1']]=[100.1,100.2,100.3]
        y[workflow.a2.Y['H_w1']]=98.25
        u[0]=.5
        u[1]=0.
        residuals,rows=workflow.constraint_audit(t,x,y,u,self.theta,self.params,cfg)
        by_name={r['label']:r for r in rows}
        row=by_name['gate1_head_opening']
        self.assertAlmostEqual(row['scaled_residual'],(.1875-.1)/2)
        self.assertAlmostEqual(row['margin'],-(.1875-.1))
        self.assertEqual(row['worst_time_s'],0.)
        self.assertFalse(row['pass_at_tolerance'])
        self.assertEqual(by_name['gate2_head_opening']['scaled_residual'],0.)
        self.assertTrue(by_name['gate2_head_opening']['pass_at_tolerance'])

    def test_audit_catches_reconstructed_bounds_outside_residual_vector(self):
        cfg,t,x,y,u=self.artificial_audit_arrays()
        u[0]=[-.1,.2,1.2]
        theta=self.theta.copy();theta[-2:]=[-.2,1.3]
        residuals,rows=workflow.constraint_audit(t,x,y,u,theta,self.params,cfg)
        self.assertEqual(residuals.shape,(workflow.CONSTRAINT_COUNT,))
        by_name={r['label']:r for r in rows}
        self.assertAlmostEqual(by_name['gate1_lower']['scaled_residual'],.1)
        self.assertAlmostEqual(by_name['gate1_upper']['scaled_residual'],.2)
        self.assertAlmostEqual(by_name['cover1_lower']['scaled_residual'],.2)
        self.assertAlmostEqual(by_name['cover2_upper']['scaled_residual'],.3)
        self.assertFalse(workflow.is_feasible({'audit':rows},cfg))

    def test_missing_safety_stays_an_error_at_constraint_assembly(self):
        cfg,t,x,y,u=self.artificial_audit_arrays()
        cfg['safety']['dam_m'] = [None,None]
        cfg['safety']['wetland_m'] = [None,None]
        with self.assertRaisesRegex(workflow.PreflightError,'Missing safety'):
            workflow.constraint_audit(t,x,y,u,self.theta,self.params,cfg)

    def test_numerical_scales_must_be_positive(self):
        cfg=copy.deepcopy(self.config);cfg['numerics']['volume_scales_m3'][2]=0.
        with self.assertRaises(workflow.PreflightError):workflow.preflight(cfg)
        cfg=copy.deepcopy(self.config);cfg['numerics']['height_scale_m']=-1.
        with self.assertRaises(workflow.PreflightError):workflow.preflight(cfg)

    def test_optimizer_settings_match_adapter_and_feasibility_units(self):
        self.assertFalse(set(self.config['optimizer'])-set(DEFAULT_SETTINGS))
        self.assertEqual(self.config['optimizer']['feasibility_tol'],self.config['numerics']['feasibility_tol'])
        self.assertGreaterEqual(self.config['optimizer']['max_evals'],self.theta.size+2)

    def test_base_parameters_are_fresh_and_overrides_require_provenance(self):
        cfg=copy.deepcopy(self.config)
        cfg['parameter_overrides']={}  # TEST-ONLY baseline object construction, no project approval.
        first=workflow.base_parameters(cfg);second=workflow.base_parameters(cfg)
        first.R_veg['type_a'][0]=999.
        first.A_w[0]=1.
        self.assertEqual(second.R_veg['type_a'][0],2.8)
        self.assertEqual(second.A_w[0],12000.)
        cfg['parameter_overrides']['A_u']={'value':1.,'source_id':'TEST_ONLY','reason':'unit test','approved':False}
        with self.assertRaises(workflow.PreflightError):workflow.base_parameters(cfg)
        cfg['parameter_overrides']={'rho':{'value':[0.,0.],'source_id':'TEST_ONLY','reason':'unit test','approved':True}}
        with self.assertRaises(workflow.PreflightError):workflow.base_parameters(cfg)

    def test_simulation_passes_fixed_initial_state_and_fresh_parameters_to_solver(self):
        # Mock the ODE and gates. No artificial gate/safety values enter a saved config.
        cfg,t_stub,x_stub,y_stub,u_stub=self.artificial_audit_arrays()
        original=copy.deepcopy(cfg)
        observed=[]
        def fake_solver(rhs,t,initial,u,params):
            observed.append((params,initial.copy(),u.copy()))
            n=len(t);x=np.repeat(initial[:,None],n,axis=1);dx=np.zeros_like(x)
            y=np.zeros((len(workflow.a2.OUTPUT_NAMES),n))
            for label,value in [('H_w1',98.25),('H_w2',98.15),('H_d1',100.3),('H_d2',100.1),('H_g',97.2),('H_r',97.35)]:
                y[workflow.a2.Y[label]]=value
            # Deliberately poison one local object: a later candidate must start fresh.
            params.R_veg['type_a'][0]=999.
            return t,x,dx,y
        fake_gates=[SimpleNamespace(gate_sill_z=100.05),SimpleNamespace(gate_sill_z=99.85)]
        with patch.object(workflow,'gate_parameters',return_value=fake_gates), \
             patch.object(workflow.a2,'validate'), \
             patch.object(workflow.a2,'load_ode4u',return_value=fake_solver):
            workflow.simulate(self.theta,('type_a','type_a'),cfg)
            changed=self.theta.copy();changed[0]=.1;changed[-2]=.2
            workflow.simulate(changed,('type_a','type_a'),cfg)
        self.assertIsNot(observed[0][0],observed[1][0])
        np.testing.assert_array_equal(observed[0][1],cfg['initial_conditions']['x0_m3'])
        np.testing.assert_array_equal(observed[1][1],observed[0][1])
        np.testing.assert_array_equal(observed[0][2][2:],observed[1][2][2:])
        self.assertEqual(observed[0][0].rho[0],.65)
        self.assertEqual(observed[1][0].rho[0],.2)
        self.assertEqual(cfg,original)

    def test_evaluator_logs_cache_hit_and_returns_residual_copies(self):
        evaluator=self.evaluator()
        with patch.object(workflow,'simulate',return_value=self.mock_result()) as simulate:
            f,first=evaluator(self.theta)
            first[0]=99.
            f2,second=evaluator(self.theta.copy())
        self.assertEqual(simulate.call_count,1)
        self.assertEqual(evaluator.budget['evaluations'],1)
        self.assertEqual([r['cached'] for r in evaluator.records],[False,True])
        self.assertEqual(f,f2)
        self.assertEqual(second[0],-.1)
        self.assertEqual(evaluator.records[-1]['raw_peak_stage_m'],98.)
        self.assertEqual(evaluator.records[-1]['max_scaled_violation'],0.)

    def test_evaluator_uses_raw_peak_for_best_feasible_selection(self):
        evaluator=self.evaluator()
        outcomes=[self.mock_result(98.),self.mock_result(97.,.1),self.mock_result(97.5)]
        with patch.object(workflow,'simulate',side_effect=outcomes):
            for i in range(3):
                theta=self.theta.copy();theta[0]+=.01*i;evaluator(theta)
        self.assertEqual(evaluator.best['metrics']['peak_stage_m'],97.5)
        self.assertEqual([r['best_feasible_peak_m'] for r in evaluator.records],[98.,98.,97.5])
        self.assertEqual([r['feasible'] for r in evaluator.records],[True,False,True])

    def test_expected_failures_have_finite_constant_vector_and_no_fake_peak(self):
        failures=[FloatingPointError('TEST_ONLY inadmissible stage'),
                  workflow.TrajectoryFailure('TEST_ONLY rejected trajectory'),
                  UnsupportedGateRegime('TEST_ONLY unsupported',1.1,.2)]
        for error in failures:
            with self.subTest(error=type(error).__name__):
                evaluator=self.evaluator()
                with patch.object(workflow,'simulate',side_effect=error):f,g=evaluator(self.theta)
                self.assertTrue(np.isfinite(f));self.assertEqual(g.shape,(workflow.CONSTRAINT_COUNT,))
                self.assertTrue(np.all(np.isfinite(g)));self.assertTrue(np.all(g>0))
                self.assertIsNone(evaluator.best)
                self.assertIsNone(evaluator.records[0]['raw_peak_stage_m'])
                self.assertFalse(evaluator.records[0]['feasible'])
                self.assertIn(type(error).__name__,evaluator.records[0]['failure'])

    def test_hydraulically_invalid_candidate_does_not_abort_next_evaluation(self):
        evaluator=self.evaluator()
        diagnostic=dict(branch=1,time_s=123.,water_depth_m=.1,
                        opening_height_m=.15,residual_m=.0875)
        error=workflow.HydraulicApplicabilityFailure(diagnostic)
        valid_theta=self.theta.copy();valid_theta[0]*=.1
        with patch.object(workflow,'simulate',side_effect=[error,self.mock_result(98.)]):
            bad_f,bad_c=evaluator(self.theta)
            good_f,good_c=evaluator(valid_theta)
        self.assertEqual(bad_c.shape,good_c.shape)
        self.assertEqual(bad_c.shape,(workflow.CONSTRAINT_COUNT,))
        self.assertTrue(np.all(bad_c>0.))
        self.assertTrue(np.all(good_c<=0.))
        self.assertIsNone(evaluator.records[0]['raw_peak_stage_m'])
        self.assertTrue(evaluator.records[0]['hydraulically_invalid'])
        self.assertIn('123.0',evaluator.records[0]['failure'])
        self.assertEqual(evaluator.records[1]['raw_peak_stage_m'],98.)
        self.assertFalse(evaluator.records[1]['hydraulically_invalid'])
        self.assertTrue(evaluator.records[1]['feasible'])
        self.assertEqual(evaluator.budget['evaluations'],2)

    def test_unexpected_programming_errors_escape_unchanged(self):
        for error in (KeyError('TEST_ONLY bug'),ValueError('TEST_ONLY bad shape')):
            with self.subTest(error=type(error).__name__):
                evaluator=self.evaluator()
                with patch.object(workflow,'simulate',side_effect=error):
                    with self.assertRaises(type(error)) as caught:evaluator(self.theta)
                self.assertIs(caught.exception,error)
                self.assertEqual(evaluator.records,[])

    def test_evaluator_copies_config_and_includes_pair_in_identity(self):
        evaluator=self.evaluator()
        self.config['numerics']['objective_origin_m']=-1000.
        self.assertEqual(evaluator.config['numerics']['objective_origin_m'],97.)
        with patch.object(workflow,'simulate',return_value=self.mock_result()) as simulate:
            evaluator(self.theta)
            evaluator.pair=('type_b','type_a')
            evaluator(self.theta)
        self.assertEqual(simulate.call_count,2)
        self.assertEqual([r['category_pair'] for r in evaluator.records],['type_a/type_a','type_b/type_a'])

    def test_shared_evaluation_and_wall_budget_stop_before_simulation(self):
        for kind in ('evaluations','time'):
            with self.subTest(kind=kind):
                evaluator=self.evaluator()
                if kind=='evaluations':evaluator.budget['evaluations']=self.config['numerics']['total_ode_budget']
                else:evaluator.budget['start']-=self.config['numerics']['total_seconds']+1
                with patch.object(workflow,'simulate') as simulate:
                    with self.assertRaises(workflow.RunBudgetExceeded):evaluator(self.theta)
                simulate.assert_not_called()

    def test_refinement_reuses_fixed_forcing_knots(self):
        cfg=copy.deepcopy(self.config)
        cfg['parameter_overrides']={}
        cfg['forcing']['kind']='a2_synthetic_fixed_samples'
        base=workflow.base_parameters(cfg)
        t,u=workflow.fixed_forcing(cfg,base)
        fine=workflow.base_parameters(cfg);fine.time_step/=2
        finer=np.vstack([np.interp(fine.time,t,row) for row in u])
        np.testing.assert_array_equal(finer[:,::2],u)
        self.assertAlmostEqual(np.trapezoid(finer[0],fine.time),np.trapezoid(u[0],t),places=14)


    def test_optimizer_and_independent_tolerance_mismatch_is_blocked(self):
        cfg=copy.deepcopy(self.config)
        cfg['optimizer']['feasibility_tol']=2*cfg['numerics']['feasibility_tol']
        with self.assertRaisesRegex(workflow.PreflightError,'tolerances must match'):
            workflow.preflight(cfg)

    def test_gap_text_matches_implemented_submergence_blending(self):
        gaps={g['id']:g for g in workflow.preflight(self.config)}
        self.assertIn('submergence blending is implemented',gaps['G5']['missing'])
        self.assertIn('weir',gaps['G5']['missing'])

    def test_search_reserves_three_fresh_verification_integrations(self):
        evaluator=self.evaluator()
        evaluator.budget['evaluations']=self.config['numerics']['total_ode_budget']-3
        with patch.object(workflow,'simulate') as simulate:
            with self.assertRaises(workflow.RunBudgetExceeded):evaluator(self.theta)
        simulate.assert_not_called()
        self.assertEqual(evaluator.budget['evaluations'],self.config['numerics']['total_ode_budget']-3)

    def test_unapproved_project_cannot_start_optimizer_or_write_run_files(self):
        with self.assertRaisesRegex(workflow.PreflightError,'Unapproved/unsourced'):
            workflow.base_parameters(self.config)
        with patch('B.official_opt.run_opt') as optimizer, \
             patch.object(workflow,'simulate') as simulate, \
             patch.object(workflow,'save_json') as save:
            with self.assertRaises(workflow.PreflightError):workflow.run_project(self.config)
        optimizer.assert_not_called();simulate.assert_not_called();save.assert_not_called()

    def mocked_project_run(self, reference_fine_peak=98.0001, expire_after_fresh=False,
                           reference_feasible=True, events=None, original_hydraulic_failure=False):
        """Exercise run orchestration entirely in memory, with every simulation mocked."""
        cfg=copy.deepcopy(self.config)
        cfg['categories']['admissible']=[['type_a'],['type_a']]
        now=[0.]
        calls=[]
        saved={}
        outcomes=[98.,97.5,97.5,97.5001,reference_fine_peak]
        expected_start=self.theta.copy()
        if original_hydraulic_failure:
            failure=workflow.HydraulicApplicabilityFailure(dict(branch=2,time_s=123.,
                water_depth_m=.1,opening_height_m=.15,residual_m=.0875))
            outcomes=[failure,99.,failure,98.,98.5,97.5,97.5,97.5001,reference_fine_peak]
            expected_start[:-2]*=.1
            self.assertEqual(cfg['reference_search']['uniform_scales'],[0.,.25,.1,.05])
        def fake_simulate(theta,pair,config,**kwargs):
            index=len(calls)
            if events is not None:events.append('reference' if index==0 else f'simulation_{index}')
            calls.append({'theta':np.array(theta),'pair':tuple(pair),'refine':kwargs.get('refine',1)})
            if isinstance(outcomes[index],Exception):raise outcomes[index]
            result=self.mock_result(outcomes[index],.1 if index==0 and not reference_feasible else -.1)
            result.update(theta=np.array(theta),pair=list(pair),time=np.array([0.,1.]),
                          x=np.ones((8,2)),u=np.zeros((10,2)),y=np.zeros((37,2)))
            if expire_after_fresh and index==2:now[0]=cfg['numerics']['total_seconds']+1
            return result
        def fake_opt(callback,x0,lower,upper,settings,constants=None):
            if events is not None:events.append('optimizer')
            self.assertIsInstance(callback,workflow.Evaluator)
            np.testing.assert_array_equal(x0,expected_start)
            np.testing.assert_array_equal(lower,np.zeros(26))
            expected_upper=np.ones(26)
            phase='full_box' if optimizer.call_count==1 else 'reference_neighborhood'
            if phase=='reference_neighborhood':
                expected_upper[:-2]=np.minimum(1.,np.maximum(cfg['search']['minimum_local_upper'],
                    cfg['search']['local_gate_upper_multiplier']*expected_start[:-2]))
            np.testing.assert_array_equal(upper,expected_upper)
            expected_settings=dict(cfg['optimizer'])
            expected_settings['max_evals']=cfg['search'][phase+'_evaluations']
            self.assertEqual(settings,expected_settings)
            # Deliberately supply only guaranteed adapter status: driver must use
            # independently feasible callback records rather than a fake last x.
            return {'status':'evaluation_budget'}
        def save(path,data):saved[path.name]=copy.deepcopy(data)
        def gate_plot(reference,config):
            self.assertEqual(reference['metrics']['peak_stage_m'],98.)
            self.assertEqual(config,cfg)
            if events is not None:events.append('gate_rating')
        with patch.object(workflow,'preflight',return_value=[]), \
             patch.object(workflow,'base_parameters',return_value=self.params), \
             patch.object(workflow,'fixed_forcing',return_value=(np.array([0.,1.]),np.zeros((8,2)))), \
             patch.object(workflow,'simulate',side_effect=fake_simulate), \
             patch('B.official_opt.run_opt',side_effect=fake_opt) as optimizer, \
             patch('B.official_opt.load_official',return_value={'revision':'TEST_ONLY','version':'TEST_ONLY'}), \
             patch.object(workflow,'save_json',side_effect=save), \
             patch.object(workflow,'write_table'), \
             patch('B.exploratory_reporting.export_case') as export, \
             patch('B.exploratory_reporting.save_run'), \
             patch.object(workflow,'export_gate_rating',side_effect=gate_plot) as gate_export, \
             patch.object(workflow.np,'savez_compressed'), \
             patch.object(workflow.clock,'monotonic',side_effect=lambda:now[0]), \
             patch.object(workflow.Path,'mkdir'):
            workflow.run_project(cfg)
        self.assertEqual(optimizer.call_count,2)
        self.assertEqual(gate_export.call_count,1)
        if saved['run_status.json']['status']=='refined_verification_failed':
            self.assertIn('NOT verified',export.call_args.args[-1])
        else:
            self.assertIn('Best feasible found',export.call_args.args[-1])
        return saved,calls,export.call_count

    def test_invalid_original_baseline_uses_labeled_reference_search_without_clipping(self):
        saved,calls,exports=self.mocked_project_run(original_hydraulic_failure=True)
        self.assertEqual(len(calls),9)
        self.assertEqual(saved['original_baseline_check.json']['status'],'hydraulically_invalid')
        self.assertIsNone(saved['original_baseline_check.json']['physical_peak_m'])
        search=saved['reference_search.json']
        self.assertEqual(search['selected_scale'],.1)
        self.assertIn('Replacement reference',search['label'])
        self.assertIn('feasibility search',search['label'])
        self.assertEqual([row['status']=='rejected' for row in search['trials']],[False,True,False,False])
        for call,scale in zip(calls[1:5],[0.,.25,.1,.05]):
            np.testing.assert_array_equal(call['theta'][:-2],self.theta[:-2]*scale)
            np.testing.assert_array_equal(call['theta'][-2:],self.theta[-2:])
            self.assertEqual(call['pair'],tuple(self.config['reference']['categories']))
        self.assertTrue(saved['run_status.json']['reference_replaced'])
        self.assertEqual(saved['run_status.json']['ode_evaluations'],9)
        self.assertEqual(exports,2)

    def test_reference_grid_sensitivity_failure_blocks_optimized_outputs(self):
        saved,calls,exports=self.mocked_project_run(reference_fine_peak=98.1)
        self.assertEqual(len(calls),5)
        self.assertEqual(saved['run_status.json']['status'],'refined_verification_failed')
        self.assertIn('Reference or candidate',saved['run_status.json']['reason'])
        self.assertFalse(saved['run_status.json']['reported_feasible'])
        self.assertEqual(exports,2)  # reference plus explicitly unverified coarse report

    def test_verification_wall_budget_is_checked_before_each_integration(self):
        saved,calls,exports=self.mocked_project_run(expire_after_fresh=True)
        self.assertEqual(len(calls),3)  # reference, exact start, fresh candidate
        self.assertEqual(saved['run_status.json']['status'],'refined_verification_failed')
        self.assertIn('budget exhausted',saved['run_status.json']['reason'])
        self.assertEqual(exports,2)  # reference plus explicitly unverified coarse report

    def test_driver_uses_adapter_status_and_counts_all_verification_calls(self):
        saved,calls,exports=self.mocked_project_run()
        status=saved['run_status.json']
        self.assertEqual(status['status'],'best_feasible_solution_found_under_stated_configuration')
        self.assertEqual(status['ode_evaluations'],5)
        self.assertEqual([r['refine'] for r in calls],[1,1,1,2,2])
        self.assertAlmostEqual(status['peak_grid_difference_m'],.0001)
        self.assertAlmostEqual(status['reference_peak_grid_difference_m'],.0001)
        self.assertEqual(exports,2)  # baseline report precedes final comparison
        self.assertFalse(status['continuous_time_guarantee'])
        self.assertTrue(status['perfect_foresight'])


    def test_constraint_sources_resolve_to_a2_and_config_without_overrides(self):
        cfg,t,x,y,u=self.artificial_audit_arrays()
        _,rows=workflow.constraint_audit(t,x,y,u,self.theta,self.params,cfg)
        sources={r['label']:r['source_id'] for r in rows}
        for label in ['soil_capacity','groundwater_geometry','groundwater_soil',
                      'groundwater_wetland_bed','river_geometry']:
            self.assertEqual(sources[label],'A2_SOURCE')
        self.assertEqual(sources['integration_status'],'B2_CONFIG')

    def test_constraint_sources_follow_each_parameter_override_independently(self):
        mapping={'C_m':'soil_capacity','H_g_valid_max':'groundwater_geometry',
                 'soil_bottom':'groundwater_soil','z_w':'groundwater_wetland_bed',
                 'H_r_valid_max':'river_geometry'}
        for parameter,label in mapping.items():
            with self.subTest(parameter=parameter):
                cfg,t,x,y,u=self.artificial_audit_arrays()
                cfg['parameter_overrides'][parameter]={
                    'value':copy.deepcopy(getattr(self.params,parameter)),
                    'source_id':'TEST_ONLY_'+parameter,
                    'reason':'Unit-test provenance dispatch only; no project parameter change',
                    'approved':True}
                _,rows=workflow.constraint_audit(t,x,y,u,self.theta,self.params,cfg)
                sources={r['label']:r['source_id'] for r in rows}
                self.assertEqual(sources[label],'TEST_ONLY_'+parameter)
                for other in set(mapping.values())-{label}:
                    self.assertEqual(sources[other],'A2_SOURCE')
                self.assertEqual(sources['integration_status'],'B2_CONFIG')

    def test_gate_rating_is_after_revised_reference_and_before_optimizer(self):
        events=[]
        self.mocked_project_run(events=events)
        self.assertLess(events.index('reference'),events.index('gate_rating'))
        self.assertLess(events.index('gate_rating'),events.index('optimizer'))
        self.assertEqual(events.count('gate_rating'),1)

    def test_constraint_infeasible_reference_still_allows_search(self):
        events=[]
        saved,calls,exports=self.mocked_project_run(reference_feasible=False,events=events)
        self.assertFalse(saved['pilot.json']['reference_feasible'])
        self.assertIn('optimizer',events)
        self.assertIn('gate_rating',events)
        self.assertEqual(len(calls),5)
        self.assertEqual(exports,2)  # baseline report precedes final comparison

    def test_failed_reference_stops_before_gate_rating_and_optimizer(self):
        failures=[workflow.TrajectoryFailure('TEST_ONLY reference failure'),
                  FloatingPointError('TEST_ONLY invalid RK state'),
                  UnsupportedGateRegime('TEST_ONLY unimplemented reference',1.1,.2)]
        for error in failures:
            with self.subTest(error=type(error).__name__):
                cfg=copy.deepcopy(self.config)
                cfg['categories']['admissible']=[['type_a'],['type_a']]
                saved={}
                def save(path,data):saved[path.name]=copy.deepcopy(data)
                with patch.object(workflow,'preflight',return_value=[]), \
                     patch.object(workflow,'base_parameters',return_value=self.params), \
                     patch.object(workflow,'fixed_forcing',return_value=(np.array([0.,1.]),np.zeros((8,2)))), \
                     patch.object(workflow,'simulate',side_effect=error) as simulate, \
                     patch('B.official_opt.run_opt') as optimizer, \
                     patch('B.official_opt.load_official',return_value={'revision':'TEST_ONLY','version':'TEST_ONLY'}), \
                     patch.object(workflow,'export_gate_rating') as gate_export, \
                     patch('B.exploratory_reporting.export_case') as comparison, \
                     patch('B.exploratory_reporting.save_run'), \
                     patch.object(workflow,'save_json',side_effect=save), \
                     patch.object(workflow,'write_table'), \
                     patch.object(workflow.Path,'mkdir'):
                    workflow.run_project(cfg)
                self.assertEqual(simulate.call_count,1)
                optimizer.assert_not_called();gate_export.assert_not_called();comparison.assert_not_called()
                self.assertEqual(saved['run_status.json']['status'],'reference_simulation_failed')
                self.assertEqual(saved['run_status.json']['ode_evaluations'],1)
                self.assertFalse(saved['pilot.json']['reference_feasible'])


if __name__=='__main__':
    unittest.main()
