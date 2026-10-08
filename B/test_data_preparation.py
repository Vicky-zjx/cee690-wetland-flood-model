"""Artificial QC/unit-conversion fixtures; no local observations fabricated."""
import unittest
import numpy as np
from B.prepare_b2_data import utc_seconds,resample_rainfall,reconcile_head

class DataPreparationTests(unittest.TestCase):
    def test_depth_conserved_when_split(self):
        result=resample_rainfall([(0,3600,12,True)],[0,900,1800,2700,3600])
        np.testing.assert_allclose(result,[3,3,3,3]);self.assertAlmostEqual(sum(result),12)
    def test_depth_conserved_when_aggregated(self):
        np.testing.assert_allclose(resample_rainfall([(0,60,2,True),(60,120,5,True)],[0,120]),[7])
    def test_zero_is_observed_not_missing(self):
        np.testing.assert_array_equal(resample_rainfall([(0,60,0,True)],[0,60]),[0])
        with self.assertRaises(ValueError):resample_rainfall([(0,60,np.nan,True)],[0,60])
    def test_gap_qc_and_overlap_fail(self):
        for rows in [[(0,30,1,True)],[(0,60,1,False)],[(0,60,1,True),(30,90,1,True)]]:
            with self.subTest(rows=rows),self.assertRaises(ValueError):resample_rainfall(rows,[0,60])
    def test_partial_source_interval_is_explicit_uniform_assumption(self):
        np.testing.assert_allclose(resample_rainfall([(0,100,10,True)],[20,60]),[4])
    def test_dst_offsets_not_inferred(self):
        self.assertEqual(utc_seconds('2022-11-06T01:30:00-06:00')-utc_seconds('2022-11-06T01:30:00-05:00'),3600)
        with self.assertRaises(ValueError):utc_seconds('2022-11-06T01:30:00')
    def test_utc_equivalence(self):
        self.assertEqual(utc_seconds('2022-09-10T00:00:00Z'),utc_seconds('2022-09-09T19:00:00-05:00'))
    def test_datum_requires_documented_offset(self):
        with self.assertRaises(ValueError):reconcile_head([1,2],None)
        np.testing.assert_allclose(reconcile_head([1,2],100),[101,102])
if __name__=='__main__':unittest.main()
