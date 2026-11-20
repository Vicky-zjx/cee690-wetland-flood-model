"""Preserve earlier project studies and the original professor/D.1 sources."""
from pathlib import Path
import argparse
import hashlib
import json
import os

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
MANIFEST=HERE/'preservation_before.json'
EXTERNAL=[Path('/Users/zhangjiaxuan/Downloads/fair_lite_demo.py'),
          Path('/Users/zhangjiaxuan/Downloads/SMCR-JiaxuanZhang-HWD.1-2026.pdf')]
SKIP={'.git','.venv','venv','__pycache__','.pytest_cache','.mypy_cache','.ruff_cache','.cache','node_modules'}


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def files():
    for directory,dirs,names in os.walk(ROOT):
        current=Path(directory)
        dirs[:]=sorted(d for d in dirs if d not in SKIP and current/d != HERE)
        for name in sorted(names):
            p=current/name
            if p.is_file() and p.suffix not in {'.pyc','.pyo'} and name!='.DS_Store':yield p
    yield from EXTERNAL


def snapshot():
    values={str(p):digest(p) for p in files()}
    with MANIFEST.open('x') as f:json.dump(values,f,indent=2,sort_keys=True);f.write('\n')
    return dict(files_recorded=len(values))


def verify():
    before=json.loads(MANIFEST.read_text());now={str(p):p for p in files()}
    changed=sorted(k for k in before.keys() & now.keys() if digest(now[k])!=before[k])
    missing=sorted(before.keys()-now.keys());added=sorted(now.keys()-before.keys())
    return dict(files_checked=len(before),changed=changed,missing=missing,new_protected_files=added,
        passed=not(changed or missing or added),scope='All prior project files plus original professor source and submitted D.1; current antecedent_conditions directory excluded')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['snapshot','verify']);p.add_argument('--output',type=Path)
    args=p.parse_args();result=snapshot() if args.action=='snapshot' else verify();text=json.dumps(result,indent=2)+'\n'
    if args.output:
        if HERE not in args.output.resolve().parents:raise ValueError('Save only in current study directory')
        args.output.write_text(text)
    print(text,end='')
    raise SystemExit(0 if result.get('passed',True) else 1)
