"""Synthetic parser fixtures only; no synthetic data enter the study."""
import unittest

import numpy as np
import pandas as pd

from C.uncertainty.data import align_precipitation, interpolate_met, parse_value, precipitation_intervals


class DataPreparationTests(unittest.TestCase):
    def test_qc_missing_negative_temperature_and_calm(self):
        self.assertEqual(parse_value("-0150,1", 0, 1, True), -15)
        self.assertEqual(parse_value("+0020,C", 0, 1, True), 2)
        self.assertTrue(np.isnan(parse_value("+9999,1", 0, 1, True)))
        self.assertTrue(np.isnan(parse_value("+0200,3", 0, 1, True)))
        self.assertEqual(parse_value("999,9,C,0000,5", 3, 4), 0)

    def test_no_multihour_or_overlapping_rain(self):
        records = [dict(DATE="2022-09-01T01:00:00",REPORT_TYPE="FM-15",AA1="01,0010,9,1",AA2="06,0200,9,1"),
                   dict(DATE="2022-09-01T01:30:00",REPORT_TYPE="FM-15",AA1="01,0020,9,1"),
                   dict(DATE="2022-09-01T02:00:00",REPORT_TYPE="FM-16",AA1="01,0030,9,1"),
                   dict(DATE="2022-09-01T03:00:00",REPORT_TYPE="FM-15",AA1="01,0000,2,1")]
        intervals,counts=precipitation_intervals(pd.DataFrame(records).fillna(""))
        self.assertEqual(len(intervals),2)
        self.assertEqual(counts["not_one_hour"],1)
        self.assertEqual(counts["overlapping_intervals_excluded"],1)
        self.assertTrue(intervals[-1][3])

    def test_conservative_interval_alignment_and_missing(self):
        idx=pd.date_range("2022-09-01",periods=3,freq="h",tz="UTC")
        intervals=[(idx[0]-pd.Timedelta(minutes=8),idx[1]-pd.Timedelta(minutes=8),3.,False),
                   (idx[1]-pd.Timedelta(minutes=8),idx[2]-pd.Timedelta(minutes=8),6.,False),
                   (idx[2]-pd.Timedelta(minutes=8),idx[2]+pd.Timedelta(minutes=52),0.,True)]
        values,coverage,trace=align_precipitation(intervals,idx)
        self.assertAlmostEqual(values[0],3.*52/60+6.*8/60)
        self.assertAlmostEqual(values[1],6.*52/60)
        self.assertTrue(np.isnan(values[2]))
        self.assertTrue(trace[1])
        self.assertLess(coverage[-1],1)
        self.assertAlmostEqual(values[:2].sum(),3.*52/60+6.)

    def test_overlap_rejected_and_naive_time_rejected(self):
        idx=pd.date_range("2022-09-01",periods=2,freq="h",tz="UTC")
        with self.assertRaises(ValueError):
            align_precipitation([(idx[0],idx[1],1,False),(idx[0],idx[1],1,False)],idx)
        with self.assertRaises(ValueError):
            align_precipitation([],idx.tz_localize(None))

    def test_long_gap_not_interpolated(self):
        obs=pd.to_datetime(["2022-01-01T00:00Z","2022-01-01T01:00Z","2022-01-01T05:00Z"])
        target=pd.to_datetime(["2022-01-01T00:30Z","2022-01-01T03:00Z","2021-12-31T23:00Z"])
        values,interp,gap=interpolate_met(obs,[-10,-8,2],target,2)
        self.assertEqual(values[0],-9)
        self.assertTrue(np.isnan(values[1:]).all())
        self.assertTrue(interp[0])
        self.assertFalse(gap[0])

    def test_conflicting_duplicate_is_missing(self):
        records=[dict(DATE="2022-09-01T01:00:00",REPORT_TYPE="FM-15",AA1="01,0010,9,1"),
                 dict(DATE="2022-09-01T01:00:00",REPORT_TYPE="FM-15",AA1="01,0020,9,1")]
        values,counts=precipitation_intervals(pd.DataFrame(records))
        self.assertEqual(values,[])
        self.assertEqual(counts["conflicting_duplicate_intervals"],1)


if __name__ == "__main__":
    unittest.main()
