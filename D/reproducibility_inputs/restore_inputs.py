"""Check, or restore missing byte-identical inputs to existing ignored B paths.

No simulation, resampling, overwrite of different files or configuration editing.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore', action='store_true', help='Copy missing inputs; existing different files are never overwritten.')
    args = parser.parse_args()
    records = json.loads((HERE / 'manifest.json').read_text())['files']
    checked = []
    for row in records:
        source = (HERE / row['bundle_path']).resolve()
        target = (ROOT / row['restore_path']).resolve()
        if not source.is_relative_to(HERE) or not target.is_relative_to(ROOT / 'B'):
            raise ValueError('Input path leaves its declared bundle or Stage B directory')
        if hashlib.sha256(source.read_bytes()).hexdigest() != row['sha256']:
            raise ValueError('Bundled input hash mismatch: ' + row['bundle_path'])
        if target.exists() and hashlib.sha256(target.read_bytes()).hexdigest() != row['sha256']:
            raise ValueError('Existing input differs; preserved without overwrite: ' + row['restore_path'])
        checked.append((source, target, row))
    missing = [row['restore_path'] for _, target, row in checked if not target.exists()]
    if args.restore:
        for source, target, row in checked:
            if not target.exists():
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
                if hashlib.sha256(target.read_bytes()).hexdigest() != row['sha256']:
                    raise RuntimeError('Restored input failed hash check: ' + row['restore_path'])
    print(json.dumps(dict(bundle_files_checked=len(checked), missing_before=missing,
                         restored=len(missing) if args.restore else 0,
                         note='Original model/configuration files are not edited.'), indent=2))


if __name__ == '__main__':
    main()
