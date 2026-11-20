"""Run existing component checks into this audit directory; no ODE integration."""
from pathlib import Path
import hashlib
import json
import sys
from types import FunctionType
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from D.feedback_control import verify as policy
from D.gate_hydraulics_revision import verify as hydraulic
from D.d2_analysis import analysis as metrics


def run():
    out=HERE/'component_results';out.mkdir(exist_ok=True)
    records=[]
    for label,module,key,original in (
        ('policy',policy,'HERE',ROOT/'D/feedback_control/results/policy_checks.json'),
        ('hydraulics',hydraulic,'HERE',ROOT/'D/gate_hydraulics_revision/results/hydraulic_checks.json'),
        ('metrics',metrics,'OUT',ROOT/'D/d2_analysis/results/metric_checks.json')):
        folder=out/label;folder.mkdir(exist_ok=True)
        base=module.metric_checks if label=='metrics' else module.run
        private=FunctionType(base.__code__,dict(base.__globals__,**{key:folder}),base.__name__,base.__defaults__,base.__closure__)
        private()
        destination=folder/original.name if label=='metrics' else folder/'results'/original.name
        result=json.loads(destination.read_text());prior=json.loads(original.read_text())
        records.append(dict(component=label,passed=result['passed'],count=len(result['checks']),
            saved_report_equal=result==prior,original_source=str(Path(module.__file__).relative_to(ROOT)),
            original_source_sha256=hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest(),
            saved_report=str(original.relative_to(ROOT)),saved_report_sha256=hashlib.sha256(original.read_bytes()).hexdigest(),
            current_report=str(destination.relative_to(HERE))))
    summary=dict(passed=all(r['passed'] and r['saved_report_equal'] for r in records),
        component_checks=sum(r['count'] for r in records),checks=records,ODE_integrations=0,
        method='Existing component function body, private output destination only; original modules/files unchanged')
    (HERE/'component_checks.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    if not summary['passed']:raise AssertionError('Component checks or saved-report equality failed')


if __name__=='__main__':run()
