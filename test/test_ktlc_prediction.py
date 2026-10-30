"""Keep saved-model predictions and the Gaussian handoff reproducible."""
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from qcforever import gaussian_run
from qcforever.gaussian_run.GaussianRunPack import GaussianDFTRun
from qcforever.gaussian_run.Gen_IOPline_functional import rounding_num
from qcforever.gaussian_run.predict_ktlc_mu import (
    calculate_prediction_features, load_prediction_models, predict_ktlc_mu,
)
from qcforever.util import check_resource


def test_original_five_heldout_predictions():
    # First held-out prediction in each original fold, before integration.
    cases = [('N/N=C/CNO', 0.3516414595689523),
             ('Cn1c[c-][n+](C)c1', 0.33254173723160385),
             ('OC1=CCN=CC1', 0.33081475175025554),
             ('CC(=[N-])[CH+]C(C)=O', 0.35019398094713167),
             ('[N-]=C1C=C[O+]=C1', 0.3439222560679436)]
    metadata, models = load_prediction_models()
    for model, (smiles, expected) in zip(models, cases):
        features = calculate_prediction_features(Chem.MolFromSmiles(smiles), metadata)
        assert features.shape == (1, 2265)
        assert model.predict(features, num_threads=1)[0] == pytest.approx(expected, abs=1e-12, rel=0)


def test_prediction_mean_sd_and_gaussian_handoff(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(check_resource, 'respec_cores', lambda n: n)
    mol = Chem.AddHs(Chem.MolFromSmiles('CO'))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    with Chem.SDWriter('methanol.sdf') as writer:
        writer.write(mol)
    metadata, models = load_prediction_models()
    features = calculate_prediction_features(Chem.RemoveHs(mol), metadata)
    predictions = [model.predict(features, num_threads=1)[0] for model in models]
    report = predict_ktlc_mu('methanol.sdf')
    assert set(report) == {'model_version', 'mu_mean', 'mu_std'}
    assert report['mu_mean'] == pytest.approx(np.mean(predictions), abs=1e-12)
    assert report['mu_std'] == pytest.approx(np.std(predictions, ddof=0), abs=1e-12)
    Chem.MolToXYZFile(mol, 'methanol.xyz')
    assert predict_ktlc_mu('methanol.xyz') == report

    runner = GaussianDFTRun('KTLC-BLYP-Pred', '6-31G*', 1, 'energy', 'methanol.sdf')
    runner.mem = '1GB'
    inputs = []

    def execute(name, *args, **kwargs):
        inputs.append(Path(name + '.com').read_text())
        return 'normal'

    monkeypatch.setattr(gaussian_run.Exe_Gaussian, 'exe_Gaussian', execute)
    monkeypatch.setattr(runner, 'Extract_values', lambda *args, **kwargs: {'Energy': [-1.0, 0.0]})
    monkeypatch.setattr(runner, 'LC_para_BOopt', lambda *args: pytest.fail('Prediction must not run BO'))
    result = runner.run_gaussian()
    assert result['functional_parameter_prediction'] == report
    assert result['functional_param'] == [rounding_num(report['mu_mean'])]
    assert len(inputs) == 1
    assert 'lc-blyp' in inputs[0].lower() and 'ktlc-blyp-pred' not in inputs[0].lower()
    expected = gaussian_run.Gen_IOPline_functional.functional_para('lc-blyp', result['functional_param'])
    assert expected.strip() in inputs[0]
    with pytest.warns(UserWarning, match='outside its training conditions'):
        predict_ktlc_mu('methanol.sdf', basis='sto-3g')
    with pytest.raises(ValueError, match='requires SDF or XYZ'):
        GaussianDFTRun('KTLC-BLYP-Pred', '6-31G*', 1, 'energy', 'unknown.chk').run_gaussian()
