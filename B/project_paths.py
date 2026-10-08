"""Resolve preserved provenance paths after the physical A/B stage organization.

Saved configurations and numerical records retain their original bytes and hashes.
New execution paths are stage-local. No scientific values are translated here.
"""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = ROOT.parent
A_ROOT = PROJECT_ROOT / 'A'


def resolve_path(value):
    """Accept a current path or a historical project-relative path."""
    path = Path(value)
    if path.is_absolute():
        try:
            name = path.relative_to(ROOT).as_posix()
        except ValueError:
            try:
                name = path.relative_to(PROJECT_ROOT).as_posix()
            except ValueError:
                return path
    else:
        name = path.as_posix()
    if name == '.gitignore' or name.startswith(('A/', 'B/', 'C/', 'D/')):
        return PROJECT_ROOT / name
    for old, new in (
        ('results/b2/', 'results/'), ('plots/b2/', 'plots/'),
        ('tables/b2/', 'tables/'), ('work/b2/', 'work/'), ('b2/', ''),
    ):
        if name.startswith(old):
            return ROOT / (new + name[len(old):])
    if name.startswith('work/multivarious/'):
        return A_ROOT / name
    if name.startswith('work/reference/'):
        return A_ROOT / 'reference' / name[len('work/reference/'):]
    if name in {'wetland_flood_model.py', 'algebraic_plots.py',
                'model_notes.md', 'requirements.txt',
                'results/baseline.npz', 'results/baseline_summary.json'}:
        return A_ROOT / name
    if name.startswith('plots/') and Path(name).name[:2] in {
            letter + '_' for letter in 'ABCDEFGHIJ'}:
        return A_ROOT / name
    return ROOT / name


def matches_historical_source(path, expected_sha256):
    """Require original bytes or an explicitly recorded path/import-only edit.

    Historical hashes are never replaced. The local relocation record must match
    both the expected original digest and the actual relocated source digest.
    """
    actual_path = resolve_path(path)
    actual = hashlib.sha256(actual_path.read_bytes()).hexdigest()
    if actual == expected_sha256:
        return True
    manifest = ROOT / 'work' / 'stage_relocation_manifest.json'
    if not manifest.exists():
        return False
    record = json.loads(manifest.read_text())['files'].get(str(path), {})
    return (record.get('original_sha256') == expected_sha256
            and record.get('current_sha256') == actual
            and PROJECT_ROOT / record.get('path', '') == actual_path)
