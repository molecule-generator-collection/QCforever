"""Replay saved generation batches; never regenerate or run xTB/PM6.

Compare the full-filter algorithm with incremental filtering on identical loaded
coordinates, using the current validator for both. This is not a re-evaluation
of historical energies or a claim that older RDKit/validator versions agree.
"""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import time
from unittest.mock import patch

import numpy as np
from rdkit import Chem, rdBase

from qcforever.conformer_search.config import ValidationSettings
from qcforever.conformer_search import validation


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ReplayRMSDCache:
    """Test-only memoization of an ordered pair; never used by production.

    Include the atom-ordered graph, exact coordinates and RMSD settings. This
    only speeds differential control-flow validation, not timing benchmarks.
    """
    def __init__(self):
        self.values = {}
        self.actual_calls = 0
        self.rmsd = validation.duplicate_rmsd

    def __call__(self, probe, reference, settings):
        flags = int(Chem.PropertyPickleOptions.NoProps) | int(Chem.PropertyPickleOptions.NoConformers)
        def key(mol):
            return (mol.ToBinary(flags), mol.GetConformer().GetPositions().tobytes())
        pair = (key(probe), key(reference), settings)
        if pair not in self.values:
            self.values[pair] = self.rmsd(probe, reference, settings)
            self.actual_calls += 1
        return self.values[pair]


def replay(root):
    config_path, status_path = root/'resolved_config.json', root/'status.json'
    reference_path = root/'reference.sdf'
    reference = next(iter(Chem.SDMolSupplier(str(reference_path), removeHs=False)))
    config = json.loads(config_path.read_text())
    settings = ValidationSettings(**config['validation'])
    maximum = json.loads(status_path.read_text())['budget']['maximum_candidates']
    incremental = validation.IncrementalCandidateFilter(reference, maximum, settings)
    full = []
    rows = []
    for path in sorted(root.glob('*/batch_*/raw_candidates.sdf')):
        raw = list(Chem.SDMolSupplier(str(path), removeHs=False))
        generator = path.parent.parent.name.split('_', 1)[1]
        batch_index = int(path.parent.name.split('_')[1])
        batch_status_path = path.parent/'status.json'
        seed = json.loads(batch_status_path.read_text())['seed']
        for i, mol in enumerate(raw):
            if mol is not None:
                mol.SetProp('candidate_id', f'{generator}:b{batch_index:03d}:r{i:05d}')
                mol.SetProp('generator', generator)
                mol.SetIntProp('generator_seed', seed)
        row = {'batch': str(path.relative_to(root)), 'raw_sha256': digest(path),
               'status_sha256': digest(batch_status_path), 'prefix_count': len(full), 'new_count': len(raw)}
        calls = [('full', lambda: validation.filter_candidates(full + raw, reference, maximum, settings)),
                 ('incremental', lambda: incremental.extend(raw))]
        # Alternate order to avoid always favouring one path with warm caches.
        outputs = {}
        for name, call in calls if len(rows) % 2 == 0 else reversed(calls):
            with patch.object(validation, 'duplicate_rmsd', wraps=validation.duplicate_rmsd) as counter:
                start = time.perf_counter()
                outputs[name] = call()
                row[name+'_seconds'] = time.perf_counter() - start
                row[name+'_rmsd_comparisons'] = counter.call_count
        full, old_audit = outputs['full']
        new, new_audit = outputs['incremental']
        equal = old_audit == new_audit and len(full) == len(new)
        equal = equal and all(a.GetPropsAsDict() == b.GetPropsAsDict() and
                             Chem.MolToSmiles(a) == Chem.MolToSmiles(b) and
                             np.array_equal(a.GetConformer().GetPositions(), b.GetConformer().GetPositions())
                             for a, b in zip(full, new))
        row.update(identical=equal, accepted_ids=[m.GetProp('candidate_id') for m in new],
                   rejected=new_audit['rejected'])
        rows.append(row)
        if not equal:
            row.update(full_audit=old_audit, incremental_audit=new_audit)
            break
    return {'root': str(root.resolve()), 'config_sha256': digest(config_path),
            'reference_sha256': digest(reference_path), 'maximum': maximum,
            'identical': bool(rows) and all(r['identical'] for r in rows), 'batches': rows}


def benchmark_last_batch(root):
    """Time steady-state append; validate the trusted prefix with the old path.

    Seeding private state is confined to this measurement: the full-filter run
    first proves that every prefix record survives, and supplies its audit.
    Production always builds its pool via extend().
    """
    paths = sorted(root.glob('*/batch_*/raw_candidates.sdf'))
    previous, last = paths[-2], paths[-1]
    prefix = list(Chem.SDMolSupplier(str(previous.parent/'accepted_pool.sdf'), removeHs=False))
    raw = list(Chem.SDMolSupplier(str(last), removeHs=False))
    reference = next(iter(Chem.SDMolSupplier(str(root/'reference.sdf'), removeHs=False)))
    settings = ValidationSettings(**json.loads((root/'resolved_config.json').read_text())['validation'])
    maximum = json.loads((root/'status.json').read_text())['budget']['maximum_candidates']
    row = {'scope': 'uncached local last-batch filter only; excludes prior initialization, generation and relaxation',
           'root': str(root.resolve()), 'prefix_count': len(prefix), 'new_count': len(raw),
           'raw_sha256': digest(last), 'prefix_sha256': digest(previous.parent/'accepted_pool.sdf')}
    with patch.object(validation, 'duplicate_rmsd', wraps=validation.duplicate_rmsd) as counter:
        start = time.perf_counter()
        full, audit = validation.filter_candidates(prefix+raw, reference, maximum, settings)
        row.update(full_seconds=time.perf_counter()-start, full_rmsd_comparisons=counter.call_count)
    assert all(r['accepted'] for r in audit['candidates'][:len(prefix)])
    incremental = validation.IncrementalCandidateFilter(reference, maximum, settings)
    incremental._accepted = [Chem.Mol(mol) for mol in full[:len(prefix)]]
    incremental._accepted_audit = deepcopy(audit['candidates'][:len(prefix)])
    with patch.object(validation, 'duplicate_rmsd', wraps=validation.duplicate_rmsd) as counter:
        start = time.perf_counter()
        new, new_audit = incremental.extend(raw)
        row.update(incremental_seconds=time.perf_counter()-start,
                   incremental_rmsd_comparisons=counter.call_count)
    row['identical'] = audit == new_audit and len(full) == len(new) and all(
        a.GetPropsAsDict() == b.GetPropsAsDict() and
        np.array_equal(a.GetConformer().GetPositions(), b.GetConformer().GetPositions())
        for a, b in zip(full, new))
    assert row['identical']
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cache-rmsd', action='store_true',
                        help='Reuse exact pair values for equivalence only; reported times are not speed benchmarks')
    parser.add_argument('--benchmark-last', action='store_true',
                        help='Interpret results as one conformer_search directory; time only its last append')
    args = parser.parse_args()
    if args.benchmark_last:
        report = benchmark_last_batch(args.results)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps(report, indent=2), flush=True)
        return
    report = {'rdkit_version': rdBase.rdkitVersion,
              'validation_source_sha256': digest(Path(validation.__file__)), 'cases': []}
    report['timing_scope'] = ('cached equivalence check, not a speed benchmark' if args.cache_rmsd
                              else 'uncached local filter replay, not GENKAI or full pipeline timing')
    cache = ReplayRMSDCache() if args.cache_rmsd else None
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for config in sorted(args.results.glob('*/*/*/conformer_search/resolved_config.json')):
        if cache is None:
            case = replay(config.parent)
        else:
            with patch.object(validation, 'duplicate_rmsd', cache):
                case = replay(config.parent)
            report['actual_rmsd_evaluations'] = cache.actual_calls
        report['cases'].append(case)
        report['all_identical'] = all(c['identical'] for c in report['cases'])
        args.output.write_text(json.dumps(report, indent=2)+'\n')
        print(f"{config.parent.relative_to(args.results)}: {len(case['batches'])} batches, "
              f"identical={case['identical']}", flush=True)
    if not report['cases'] or not report['all_identical']:
        raise SystemExit('Missing batches or equivalence failure; inspect report')


if __name__ == '__main__':
    main()
