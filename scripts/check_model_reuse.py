"""Explicit 2 -> 1 -> 1 raw requests: verify reuse on flexible ethanol.

This diagnostic bypasses acceptance/quota selection only to exercise additional
batches. It does not alter the production search policy or optimize molecules.
"""
import argparse
import json
from pathlib import Path
from rdkit import Chem
from qcforever.conformer_search.config import SearchConfig
from qcforever.conformer_search.model_session import ModelSession


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--model', choices=['ditmc', 'torsional_diffusion'], required=True)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    cfg = SearchConfig.resolve('high', args.config)
    session = ModelSession(args.output, cfg.generators[args.model], cfg.threads)
    ref = Chem.MolFromSmiles('CCO')
    rows = []
    try:
        for i, count in enumerate((2, 1, 1)):
            folder = args.output/f'batch_{i:03d}'
            folder.mkdir()
            records = session.generate(ref, count, cfg.seed+i*1009, folder, 2)
            assert len(records) == count, 'Not all requested raw candidates returned'
            for path in sorted(folder.glob('worker_*/model_execution.json')):
                row = json.loads(path.read_text())
                row['worker'] = path.parent.name
                row['batch'] = i
                rows.append(row)
        initial = {r['worker']: r['worker_pid'] for r in rows if r['batch'] == 0}
        assert all(r['worker_pid'] == initial[r['worker']] for r in rows)
        assert all(r['model_reused'] and r['initialization_seconds'] == 0 for r in rows if r['batch'] > 0)
        (args.output/'summary.json').write_text(json.dumps({'state': 'passed', 'smiles': 'CCO',
            'requests': [2, 1, 1], 'rows': rows}, indent=2))
    finally:
        session.close()


if __name__ == '__main__':
    main()
