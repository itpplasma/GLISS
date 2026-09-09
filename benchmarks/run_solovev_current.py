#!/usr/bin/env python3
"""Run current production FEEC on independently qualified Solov'ev exports.

The input CSV must contain q0,export. DCON Newcomb counts and GLISS negative
inertia have a common stability interpretation; raw eigenvalues do not share
a normalization. External references must therefore remain separate.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--mmax', type=int, default=8)
    parser.add_argument('--degree', type=int, default=2)
    parser.add_argument('--eigenpair', action='store_true')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    env = dict(os.environ, FO_JOBS='2', OMP_NUM_THREADS='1',
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    provenance = {'base_commit': subprocess.check_output(
        ['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip(),
        'runner_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'patch_sha256': hashlib.sha256(subprocess.check_output(
            ['git', 'diff', 'HEAD'], cwd=repo)).hexdigest(), 'runs': []}
    with args.inputs.open() as handle:
        inputs = list(csv.DictReader(handle))
    for index, item in enumerate(inputs):
        export = Path(item['export']).resolve()
        command = ['fo', 'exec', 'gliss_axisymmetric', str(export), '1',
                   str(args.mmax), '--degree', str(args.degree)]
        if not args.eigenpair:
            command.append('--count-only')
        result = subprocess.run(command, cwd=repo, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (args.output / f'run_{index:03d}.log').write_text(result.stdout)
        provenance['runs'].append({'q0': item['q0'], 'command': command,
            'export': str(export), 'returncode': result.returncode,
            'export_sha256': hashlib.sha256(export.read_bytes()).hexdigest()})
        (args.output / 'provenance.json').write_text(
            json.dumps(provenance, indent=2) + '\n')
        result.check_returncode()
        lines = result.stdout.splitlines()
        start = next(i for i, line in enumerate(lines)
                     if line.startswith('chart_metric,'))
        row = next(csv.DictReader(lines[start:start + 2]))
        row = {'q0': item['q0'], **row}
        rows.append(row)
        with (args.output / 'current_gliss.csv').open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        (args.output / 'provenance.json').write_text(
            json.dumps(provenance, indent=2) + '\n')


if __name__ == '__main__':
    main()
