"""Regression tests for shared geometry bounds and readable native failures."""
import json

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from qcforever.conformer_search.config import SearchConfig
from qcforever.conformer_search.validation import geometry_failure
from qcforever.conformer_search.pm6_errors import failure_diagnostic, PM6CalculationError
from qcforever.conformer_search.relaxation import pm6_relax


def formaldehyde():
    mol = Chem.AddHs(Chem.MolFromSmiles('C=O'))
    assert AllChem.EmbedMolecule(mol, randomSeed=4) == 0
    return mol


def test_compressed_bond_is_rejected_without_altering_normal_structure():
    mol = formaldehyde()
    cfg = SearchConfig.resolve().validation
    assert geometry_failure(mol, mol, cfg) is None
    h = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1)
    conf = mol.GetConformer()
    carbon = conf.GetAtomPosition(0)
    direction = conf.GetAtomPosition(h) - carbon
    conf.SetAtomPosition(h, carbon + direction * (0.40567 / direction.Length()))
    assert geometry_failure(mol, mol, cfg) == 'bond_compressed'
    with pytest.raises(ValueError, match='Minimum bonded'):
        SearchConfig.resolve(override={'validation': {'minimum_bonded_covalent_ratio': 2}})


@pytest.mark.parametrize('log,reason', [
    ('Small interatomic distances encountered\nError termination via Lnk1e in l202.exe', 'invalid_interatomic_distances'),
    ('Convergence failure -- run terminated.', 'scf_not_converged'),
    ('Number of steps exceeded', 'optimization_step_limit'),
    ('', 'missing_native_log'),
    ('Unrecognized output', 'missing_scf_energy'),
])
def test_pm6_failure_classification(log, reason):
    assert failure_diagnostic(log)['reason'] == reason


def test_pm6_parser_error_preserves_native_reason_and_trace(tmp_path, monkeypatch):
    from qcforever.laqa_fafoom.pyg16 import g16Object
    from pathlib import Path
    monkeypatch.setattr('shutil.which', lambda _: '/fake/g16')
    def failed_native(self):
        Path('Gau_molecule.log').write_text('Small interatomic distances encountered\nError termination\n')
        raise UnboundLocalError('energy was not assigned')
    monkeypatch.setattr(g16Object, 'run_g16', failed_native)
    with pytest.raises(PM6CalculationError, match='invalid_interatomic_distances'):
        pm6_relax(formaldehyde(), tmp_path, 0, 1, 4, '1GB', {})
    diagnostic = json.loads((tmp_path/'failure.json').read_text())
    assert 'UnboundLocalError' in diagnostic['adapter_error']
    assert (tmp_path/'Gau_molecule.log').exists() and (tmp_path/'trace.json').exists()


def test_optional_worker_cli_needs_no_core_or_ml_imports():
    import subprocess
    import sys
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    code = ('import sys; import qcforever_model_workers.worker; '
            'assert all(x not in sys.modules for x in ("qcforever", "torch", "jax", "rdkit"))')
    # -S disables site-packages, so accidental eager model/core imports fail.
    subprocess.run([sys.executable, '-S', '-c', code], cwd=root, check=True)
    result = subprocess.run([sys.executable, '-S', '-m', 'qcforever_model_workers.worker', '--help'],
                            cwd=root, check=True, capture_output=True, text=True)
    assert '--checkpoint' in result.stdout and '--legacy-workflow' not in result.stdout
