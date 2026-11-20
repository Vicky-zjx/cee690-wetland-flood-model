"""Hash protection for A/B/C, prior D studies and root files.

The current D/gate_hydraulics_revision directory is deliberately excluded.
Bytecode, virtual environments, Git internals and runtime cache directories are
excluded; scientific data, sources, figures, results and external solver source
files are included. No protected file is written by this module.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MANIFEST = HERE / 'preservation_before.json'
SKIP_DIRS = frozenset({'.git', '.venv', 'venv', '__pycache__', '.pytest_cache',
                       '.mypy_cache', '.ruff_cache', '.cache', 'node_modules'})


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def protected_files(root=ROOT):
    root = Path(root).resolve()
    own = root / 'D/gate_hydraulics_revision'
    for stage in ('A', 'B', 'C', 'D'):
        for directory, dirs, files in os.walk(root/stage):
            current = Path(directory)
            dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS and current/d != own)
            for name in sorted(files):
                path = current/name
                if path.suffix not in {'.pyc', '.pyo'} and path.is_file():
                    yield path
    for path in sorted(root.iterdir()):
        if path.is_file() and path.name != '.DS_Store':
            yield path


def snapshot(root=ROOT, manifest=MANIFEST):
    """Create the baseline once; refuse to replace any existing baseline."""
    root, manifest = Path(root).resolve(), Path(manifest)
    values = {str(path.relative_to(root)): digest(path) for path in protected_files(root)}
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open('x') as handle:
        json.dump(values, handle, indent=2, sort_keys=True)
        handle.write('\n')
    return {'files_recorded': len(values), 'manifest': str(manifest)}


def verify(root=ROOT, manifest=MANIFEST):
    """Return exact modified/missing/new protected paths; never bless changes."""
    root, manifest = Path(root).resolve(), Path(manifest)
    before = json.loads(manifest.read_text())
    now = {str(path.relative_to(root)): path for path in protected_files(root)}
    missing = sorted(set(before)-set(now))
    added = sorted(set(now)-set(before))
    changed = sorted(rel for rel in set(before)&set(now) if digest(now[rel]) != before[rel])
    return {'files_checked':len(before), 'changed':changed, 'missing':missing,
            'new_protected_files':added, 'passed':not (changed or missing or added),
            'scope':'A/B/C, earlier D studies and root files; current revision directory excluded'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['snapshot', 'verify'])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = snapshot() if args.action=='snapshot' else verify()
    text = json.dumps(report, indent=2)+'\n'
    if args.output:
        output = args.output.resolve()
        if HERE not in output.parents:
            raise ValueError('Verification output must remain inside the new revision directory')
        output.write_text(text)
    print(text, end='')
    return 0 if report.get('passed', True) else 1


if __name__=='__main__':
    raise SystemExit(main())
