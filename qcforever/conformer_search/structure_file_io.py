"""Read input structures and write SDF records without changing chemical identity."""
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds


def read_input_structure(filename, charge):
    """Read the first SDF record or determine XYZ bonds using the supplied charge."""
    path = Path(filename).resolve()
    if path.suffix.lower() == '.sdf':
        records = [mol for mol in Chem.SDMolSupplier(str(path), removeHs=False) if mol is not None]
        if not records:
            raise ValueError('No readable input molecule')
        return records[0]
    if path.suffix.lower() == '.xyz':
        reference = Chem.MolFromXYZFile(str(path))
        if reference is None:
            raise ValueError('Unreadable XYZ input')
        rdDetermineBonds.DetermineBonds(reference, charge=charge)
        return reference
    raise ValueError('Configured optconf accepts SDF or XYZ')


def write_sdf(path, records):
    """Write available records in input order, keeping hydrogens and properties."""
    with Chem.SDWriter(str(path)) as writer:
        for mol in records:
            if mol is not None:
                writer.write(mol)
