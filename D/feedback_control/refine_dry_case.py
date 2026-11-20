"""Retest the identified near-dry mass-audit case at 15 s; no other changes.

This is a numerical follow-up to the stored 120 s low/high/B feedback runs.
The canonical nominal study separately repeats all original inputs at 15 s.
"""
import argparse
from collections import Counter
import pandas as pd
from .run_policy import HERE,run_case,save_run,write_json,original_weather,sha
from .policy import POLICY_PATH
from . import preservation


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('member',choices=['low','high','B_feedback'])
    args=parser.parse_args()
    arm='B' if args.member=='B_feedback' else 'D'
    cw={'low':2.6,'high':3.0,'B_feedback':2.8}[args.member]
    out=HERE/'results'/('fine_'+args.member)
    if out.exists():raise FileExistsError(out)
    guard=preservation.verify()
    if not guard['passed']:raise RuntimeError(str(guard))
    out.mkdir()
    plans=[dict(case=c,scenario_id='block_bootstrap_00008',design='best',dt=15) for c in ('wetland','no_wetland')]
    write_json(out/'execution_plan.json',dict(plans=plans,weather=original_weather(),hydraulics=arm,Cw_US=cw,
        control_period_s=120,policy_sha256=sha(POLICY_PATH),entrypoint_sha256=sha(__file__),
        scope='Targeted paired refinement of the near-dry balance-audit failure; no scenario exclusion or coefficient adjustment'))
    rows=[]
    for plan in plans:
        name=f'{plan["scenario_id"]}__{plan["case"]}__best__15s'
        print('START',name,flush=True)
        result=run_case(**plan,hydraulics=arm,Cw_US=cw)
        rows.append(save_run(result,out,name));pd.DataFrame(rows).to_csv(out/'scenario_diagnostics.csv',index=False)
        print('DONE',name,result['diagnostics']['classification'],result['diagnostics'].get('balance_relative'),flush=True)
    write_json(out/'counts.json',dict(Counter(f'{r["case_name"]} / {r["classification"]}' for r in rows)))
    write_json(out/'preservation_check.json',preservation.verify())


if __name__=='__main__':main()
