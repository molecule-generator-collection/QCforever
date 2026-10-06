"""MM parameters must cover the whole molecule; never silently switch MM."""
from rdkit import Chem
from rdkit.Chem import AllChem


class UnsupportedParametersError(RuntimeError):
    pass


class MissingMMDependency(RuntimeError):
    pass


def preoptimize(records, method, options):
    if method == 'none':
        return [Chem.Mol(m) for m in records], [None]*len(records)
    limit = options.get('maximum_iterations', 5000)
    if not isinstance(limit, int) or limit < 1:
        raise ValueError('MM maximum_iterations must be a positive integer')
    if method in ('uff', 'mmff94s'):
        has = AllChem.UFFHasAllMoleculeParams if method == 'uff' else AllChem.MMFFHasAllMoleculeParams
        if not all(has(m) for m in records):
            raise UnsupportedParametersError(f'{method} does not cover all molecular parameters')
        result, statuses = [], []
        for source in records:
            mol = Chem.Mol(source)
            if method == 'uff':
                status = AllChem.UFFOptimizeMolecule(mol, maxIters=limit)
            else:
                props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
                ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=0)
                status = ff.Minimize(maxIts=limit)
            result.append(mol)
            statuses.append(int(status))
        return result, statuses
    raise ValueError(f'Unknown MM method: {method}')
