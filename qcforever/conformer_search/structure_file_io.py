"""Read molecular inputs and save structures or trajectories without changing them."""
import io
import os
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


def write_sdf(path, records, *, verify=False):
    """Keep record order and properties; verify final pools before scratch removal."""
    output = io.StringIO() if verify else str(path)
    with Chem.SDWriter(output) as writer:
        for mol in records:
            if mol is not None:
                writer.write(mol)
    if verify:
        write_verified_text(path, output.getvalue())


def write_verified_text(path, text):
    """Flush and read back a temporary file before replacing its destination."""
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('x') as stream:
        stream.write(text)
        stream.flush()
        os.fsync(stream.fileno())
    if temporary.read_text() != text:
        raise OSError(f'Saved content differs from the intended output: {path}')
    temporary.replace(path)
