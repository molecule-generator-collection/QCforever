"""Predict LC-BLYP mu (bohr^-1) with the five existing RDKit/Morgan models.

The models are the original 60/20/20 train/validation/test fits, not retrained
models. Their population standard deviation describes ensemble disagreement,
not a calibrated prediction uncertainty. No quantum calculation is performed.
"""
from functools import lru_cache
import json
from pathlib import Path
import warnings

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdFingerprintGenerator

MODEL_DIRECTORY = Path(__file__).with_name('ktlc_mu_models')


def predict_ktlc_mu(filename, charge=0, multiplicity=1, basis='6-31g*', solvent=None):
    """Read SDF/XYZ and return model_version, mu_mean and mu_std only."""
    from qcforever.conformer_search.structure_file_io import read_input_structure

    if Path(filename).suffix.lower() not in {'.sdf', '.xyz'}:
        raise ValueError('KTLC-BLYP-Pred requires SDF or XYZ molecular input')
    mol = Chem.RemoveHs(read_input_structure(filename, charge))
    if Chem.GetFormalCharge(mol) != charge:
        raise ValueError('KTLC-BLYP-Pred: specified charge differs from the molecular graph')
    outside_training = (
        charge != 0 or multiplicity != 1 or Descriptors.NumRadicalElectrons(mol) != 0
        or len(Chem.GetMolFrags(mol)) != 1 or mol.GetNumAtoms() == 0
        or not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
        or any(atom.GetSymbol() not in {'H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Br'}
               for atom in mol.GetAtoms())
        or Descriptors.MolWt(mol) > 1000 or basis.lower() not in {'6-31g*', '6-31g(d)'}
        or solvent not in (None, '0')
    )
    if outside_training:
        warnings.warn('KTLC-BLYP-Pred is outside its training conditions: neutral closed-shell '
                      'single organic molecules up to MW 1000, H/B/C/N/O/F/Si/P/S/Cl/Br, '
                      'gas-phase LC-BLYP/6-31G*. Prediction accuracy is unvalidated here.',
                      stacklevel=2)
    metadata, models = load_prediction_models()
    features = calculate_prediction_features(mol, metadata)
    predictions = np.array([model.predict(features, num_threads=1)[0] for model in models])
    if not np.isfinite(predictions).all():
        raise ValueError('KTLC-BLYP-Pred returned a non-finite mu')
    return {'model_version': metadata['model_version'],
            'mu_mean': float(predictions.mean()), 'mu_std': float(predictions.std(ddof=0))}


def calculate_prediction_features(mol, metadata):
    """Match training descriptor order, float32 precision and NaN handling."""
    values = []
    for name in metadata['descriptor_names']:
        descriptor = getattr(Descriptors, name, None)
        if descriptor is None:
            raise RuntimeError(f'RDKit descriptor {name} is unavailable; '
                               f'install RDKit >= {metadata["rdkit_version"]}')
        try:
            values.append(descriptor(mol))
        except Exception:
            # Training CalcMolDescriptors used missingVal=NaN as well.
            values.append(float('nan'))
    descriptors = np.asarray(values, dtype=np.float32)
    descriptors[~np.isfinite(descriptors)] = np.nan
    generator = rdFingerprintGenerator.GetMorganGenerator(**metadata['morgan'])
    fingerprint = generator.GetFingerprintAsNumPy(mol)
    return np.concatenate([descriptors, fingerprint]).reshape(1, -1)


@lru_cache(maxsize=1)
def load_prediction_models():
    """Load the bundled text models once; ordinary DFT/BO needs no model load."""
    import lightgbm as lgb

    metadata = json.loads((MODEL_DIRECTORY / 'metadata.json').read_text())
    models = [lgb.Booster(model_file=str(MODEL_DIRECTORY / f'fold_{fold}.txt'))
              for fold in range(5)]
    return metadata, models
