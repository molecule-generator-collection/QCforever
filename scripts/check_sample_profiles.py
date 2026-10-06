"""Four Samples, middle/high: generation/MM and optional real QC backends.

Generate once per profile, then feed the exact same MM candidates to each
backend. Never replace a failed learned generator with a hidden substitute.
"""
import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter

from rdkit import Chem
from qcforever.conformer_search.config import SearchConfig
from qcforever.conformer_search.pipeline import prepare_candidates, save, write_sdf, PreparationResult
from qcforever.conformer_search.bridge import working_directory
from qcforever.conformer_search.relaxation import relax_candidates
from qcforever.util import job_timeout


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cores', type=int, default=8)
    parser.add_argument('--backends', nargs='*', choices=['xtb', 'pm6'], default=[])
    parser.add_argument('--qc-cores', type=int, default=4,
                        help='Cores per xTB/PM6 calculation; --cores is the total allocation')
    parser.add_argument('--samples', nargs='+', default=['ch2o', 'Chlorobenzene', 'ethanol', '200_11'])
    parser.add_argument('--profiles', nargs='+', choices=['middle', 'high'], default=['middle', 'high'])
    parser.add_argument('--smoke-candidates', type=int, choices=[1, 2],
                        help='Test-only override: never edits the normal candidate budget')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    args.output.mkdir(parents=True, exist_ok=False)
    results = []
    for name in args.samples:
        source = repo/'Samples'/f'{name}.sdf'
        reference = Chem.SDMolSupplier(str(source), removeHs=False)[0]
        for profile in args.profiles:
            config = SearchConfig.resolve(profile, args.config)
            config = SearchConfig.resolve(profile, {**config.to_mapping(), 'relaxation': {
                **config.relaxation, 'cores_per_calculation': args.qc_cores}})
            if args.smoke_candidates is not None:
                config = SearchConfig.resolve(profile, {
                    **config.to_mapping(), 'budget': {
                        'formula': 'fixed', 'fixed': args.smoke_candidates,
                        'maximum': args.smoke_candidates, 'raw_attempt_multiplier': 2}})
            started = perf_counter()
            row = {'sample': name, 'profile': profile, 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                   'electronic_relaxation': 'not_requested',
                   'test_only_candidate_override': args.smoke_candidates,
                   'normal_budget': SearchConfig.resolve(profile, args.config).budget.resolve(Chem.RemoveHs(reference))}
            try:
                out = prepare_candidates(reference, args.output/name/profile, config, allocated_cores=args.cores)
                status = json.loads(out.status_path.read_text())
                row.update(state=out.state, candidates=len(out.candidates), maximum_candidates=status['budget']['maximum_candidates'],
                           stages=status['stages'], mm_state=status.get('mm_state'),
                           post_mm_audit=status.get('post_mm_audit'), status_file=str(out.status_path))
                model_stages = [s for s in status['stages'] if s['generator'] != 'etkdgv3']
                row['real_model_raw_generated'] = sum(s['returned_raw'] for s in model_stages)
                row['has_model_stage_failure'] = any(s['state'] in ('unavailable_dependency', 'generation_failed') for s in model_stages)
                row['electronic'] = {}
                charge = Chem.GetFormalCharge(reference)
                multiplicity = 1 + sum(a.GetNumRadicalElectrons() for a in reference.GetAtoms())
                if out.candidates:
                    for backend in args.backends:
                        qc_root = out.status_path.parent/backend
                        qc_root.mkdir()
                        save(qc_root/'preparation_status.json', status)
                        write_sdf(qc_root/'initial_structures.sdf', out.candidates)
                        qc_prepared = PreparationResult(out.state, out.candidates,
                            qc_root/'initial_structures.sdf', qc_root/'preparation_status.json')
                        tick = perf_counter()
                        try:
                            with working_directory(qc_root), job_timeout.overall_timeout(1200):
                                audit = relax_candidates(qc_prepared, reference, config, charge,
                                    multiplicity, backend, args.cores, '1GB')
                            primary = audit['selected_index']
                            row['electronic'][backend] = {
                                'state': 'succeeded_with_structure_warning' if audit['selected_structure_warning'] else 'succeeded',
                                'structure_check_warning': audit['selected_structure_warning'],
                                'allocated_cores': audit['allocated_cores'],
                                'cores_per_calculation': audit['cores_per_calculation'],
                                'parallel_workers': audit['parallel_workers'],
                                'converged_candidates': audit['converged_candidates'],
                                'failed_candidates': audit['failed_candidates'],
                                'primary_valid_candidates': audit['primary_valid_candidates'],
                                'energy_hartree': audit['candidates'][primary]['energy_hartree'],
                                'audit_file': str(qc_root/'electronic/audit.json')}
                        except Exception as exc:
                            row['electronic'][backend] = {'state': 'failed', 'error': f'{type(exc).__name__}: {exc}'}
                        row['electronic'][backend]['wall_seconds'] = perf_counter()-tick
                        save(args.output/'current_case.json', row)
                if args.backends:
                    row['electronic_relaxation'] = 'requested'
                    states = [v['state'] for v in row['electronic'].values()]
                    row['state'] = ('failed' if not states or 'failed' in states else
                        'succeeded_with_structure_warning' if 'succeeded_with_structure_warning' in states else 'succeeded')
            except Exception as exc:
                row.update(state='failed', error=f'{type(exc).__name__}: {exc}')
            row['wall_seconds'] = perf_counter()-started
            results.append(row)
            save(args.output/'results.json', results)
            print(name, profile, row['state'], row.get('candidates'), f"{row['wall_seconds']:.2f}s", flush=True)


if __name__ == '__main__':
    main()
