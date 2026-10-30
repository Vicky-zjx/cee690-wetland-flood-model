"""Retained fitting/data-provenance tests; synthetic fixtures are not observations."""
import copy
import hashlib
import json
import unittest
from unittest.mock import patch
import pandas as pd
from C import run_fitting as workflow
from C.tests.test_weather import synthetic_training
from C.uncertainty.weather_generator import fit_weather, generate_weather


class FittingWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config={'seed':6902026,'weather':{'pit_seed':123,'block_hours':48,'burn_in_hours':336},
                    'paths':{'validation_dir':'validation'},'data':{}}
        cls.training=synthetic_training()
        cls.fitted=fit_weather(cls.training,cls.config)
        cls.fitted['implementation_sha256']=workflow.weather_source_hashes()
        cls.timestamps=pd.date_range('2022-09-09',periods=144,freq='h',tz='UTC')

    def test_weather_streams_are_reproducible_and_order_invariant(self):
        def build(task):
            method,number=task
            rng,stream=workflow.seed_stream(self.config,method,number)
            frame=generate_weather(self.fitted,self.timestamps,rng=rng,method=method,
                                   training_df=self.training,config=self.config)
            checksum=hashlib.sha256(frame.to_csv(index=False,float_format='%.17g').encode()).hexdigest()
            return task,(checksum,stream)
        tasks=[(method,i) for method in workflow.METHOD_CODES for i in range(4)]
        serial=dict(build(task) for task in tasks)
        reordered=dict(build(task) for task in reversed(tasks))
        self.assertEqual(serial,reordered)
        self.assertEqual(len({result[0] for result in serial.values()}),len(tasks))
        self.assertEqual(len({tuple(result[1]['seed_sequence_entropy']) for result in serial.values()}),len(tasks))

    def test_fit_requires_exact_training_data_and_weather_settings(self):
        workflow.verify_fit(self.fitted,self.training,self.config)
        changed=self.training.copy();changed.loc[5,'temperature_C']+=1
        with self.assertRaisesRegex(ValueError,'training observations differ'):
            workflow.verify_fit(self.fitted,changed,self.config)
        config=copy.deepcopy(self.config);config['weather']['block_hours']=72
        with self.assertRaisesRegex(ValueError,'Weather setting'):
            workflow.verify_fit(self.fitted,self.training,config)

    def test_fit_requires_implementation_version(self):
        old=copy.deepcopy(self.fitted);old.pop('implementation_sha256')
        with self.assertRaisesRegex(ValueError,'implementation changed'):
            workflow.verify_fit(old,self.training,self.config)
        with patch.object(workflow,'weather_source_hashes',return_value={'changed':'hash'}):
            with self.assertRaisesRegex(ValueError,'implementation changed'):
                workflow.verify_fit(self.fitted,self.training,self.config)

    def test_training_reports_preserve_failed_block_screen(self):
        fitted=copy.deepcopy(self.fitted)
        fitted['bootstrap_block_diagnostics']['acf_below_prespecified_point2_screen']=False
        report=workflow.training_temporal_report(fitted)
        self.assertFalse(report['bootstrap']['passed'])
        self.assertEqual(report['var']['residual_lag1'],fitted['var']['residual_lag1'])
        self.assertEqual(report['bootstrap']['criterion'],
                         'Maximum absolute training latent ACF at the configured block length must be < 0.2.')
        json.dumps(report,allow_nan=False)

if __name__=='__main__':unittest.main()
