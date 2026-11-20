"""Check the exact Git staging area for this publication; no network or mutation."""
from pathlib import Path
from collections import Counter
import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import subprocess
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from D.d2_analysis.consolidation import preservation
PREFIXES=('D/d2_analysis/robustness_checks/','D/d2_analysis/antecedent_conditions/','D/d2_analysis/consolidation/')
ADMIN={'README.md','D/README.md','D/reproducibility_inputs/README.md','.gitignore'}
PATTERNS={
 'private_key':r'-----BEGIN (?:RSA |EC |DSA |OPENSSH )?PRIVATE KEY-----',
 'github_token':r'(?<![A-Za-z0-9])(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,})',
 'aws_access_key':r'(?<![A-Z0-9])(?:AKIA|ASIA)[A-Z0-9]{16}(?![A-Z0-9])',
 'slack_token':r'xox[baprs]-[A-Za-z0-9-]{20,}',
 'openai_key':r'\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{30,}',
 'credential_url':r'https?://[^\s/:@]+:[^\s/@]{4,}@',
 'assigned_secret':r'(?i)(?:api[_-]?key|access[_-]?token|password|client[_-]?secret)\s*[=:]\s*[\"\x27][^\"\x27\s]{10,}[\"\x27]'
}


def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT)


def check():
    status=git('diff','--cached','--name-status','-z').decode().split('\0');entries=[]
    for i in range(0,len(status)-1,2):entries.append((status[i],status[i+1]))
    assert entries,'No staged files to audit'
    assert all(s in ('A','M') for s,_ in entries),'Deletion/rename requires separate review'
    assert all(p in ADMIN or p.startswith(PREFIXES) for _,p in entries),'Unrelated staged file'
    findings=[];inventory=[];text_count=0;compressed_count=0
    for state,path in entries:
        payload=git('show',':'+path);p=Path(path)
        assert payload==(ROOT/path).read_bytes(),'Staged and working bytes differ: '+path
        assert len(payload)<10*1024*1024,'Unexpected large individual file: '+path
        assert p.suffix not in ('.npz','.pyc','.pyo','.log','.pem','.key','.pptx')
        assert not path.endswith('_controller.csv.gz')
        assert not any(part in {'.git','.venv','venv','__pycache__','.cache','secrets','private'} for part in p.parts)
        if p.suffix=='.gz':payload_text=gzip.decompress(payload).decode();compressed_count+=1
        else:
            try:payload_text=payload.decode()
            except UnicodeDecodeError:payload_text=None
        if payload_text is not None:
            text_count+=1
            for name,pattern in PATTERNS.items():
                for match in re.finditer(pattern,payload_text):
                    findings.append(dict(path=path,line=payload_text[:match.start()].count('\n')+1,pattern=name))
        inventory.append(dict(status=state,path=path,bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest()))
    guard=preservation.verify()
    assert guard['passed'],guard
    # Report-only artifacts are omitted from their own digest manifest to avoid
    # recursive self-hashes. The eventual Git commit records their exact bytes.
    omit={str((HERE/'publication_audit.json').relative_to(ROOT)),str((HERE/'staged_files.csv').relative_to(ROOT))}
    rows=[r for r in inventory if r['path'] not in omit]
    with (HERE/'staged_files.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=['status','path','bytes','sha256'],lineterminator='\n');writer.writeheader();writer.writerows(rows)
    ignored=[]
    for folder in (ROOT/PREFIXES[0],ROOT/PREFIXES[1]):
        for p in folder.rglob('*'):
            if p.is_file() and subprocess.run(['git','check-ignore','-q',str(p.relative_to(ROOT))],cwd=ROOT).returncode==0:
                ignored.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size))
    report=dict(passed=not findings,base_commit=git('rev-parse','HEAD').decode().strip(),branch=git('branch','--show-current').decode().strip(),
        upstream=git('rev-parse','--abbrev-ref','@{upstream}').decode().strip(),remote=git('remote','get-url','origin').decode().strip(),
        staged_review_scope='Exact index bytes; report/manifest self-digests omitted; Git commit is final authority',
        scientific_studies_preserved=guard,nonreport_staged_files=len(rows),nonreport_staged_bytes=sum(r['bytes'] for r in rows),
        stages_by_status=dict(Counter(r['status'] for r in rows)),file_types=dict(Counter(Path(r['path']).suffix for r in rows)),
        secret_scan=dict(scope='Known credential/private-key patterns and assigned-secret literals; plain and gzip text',
            text_files_checked=text_count,gzip_text_files_checked=compressed_count,findings=findings,not_a_security_certification=True),
        ignored_new_study_files=len(ignored),ignored_new_study_bytes=sum(r['bytes'] for r in ignored),ignored_new_study_paths=ignored,
        intended_untracked_preserved=['B/presentations/Current_Wetland_Research_Comparison.png','B/presentations/Current_Wetland_Research_Comparison.pptx'],
        ODE_integrations=0,final_submission_report_written=False)
    (HERE/'publication_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('ignored_new_study_paths','scientific_studies_preserved')},indent=2))
    if findings:raise RuntimeError('Potential credentials found; stop before commit')


if __name__=='__main__':check()
