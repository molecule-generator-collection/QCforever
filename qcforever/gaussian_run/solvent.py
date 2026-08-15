"""Resolve Gaussian PCM solvent inputs against the bundled solvent table."""

import csv
import math
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, rdBase


SOLVENT_MAP_PATH = Path(__file__).with_name("solvent_map_gaussian_184.csv")


def _canonical_smiles(smiles):
    # An unsupported solvent name is also attempted as SMILES. Keep RDKit's
    # parser diagnostics from obscuring the actionable ValueError below.
    with rdBase.BlockLogs():
        molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def _load_solvent_maps():
    names = {}
    smiles = defaultdict(list)
    with SOLVENT_MAP_PATH.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            gaussian_name = row["Gaussian_Name"].strip()
            names[gaussian_name.casefold()] = gaussian_name
            canonical = _canonical_smiles(row["SMILES"].strip())
            if canonical is not None:
                smiles[canonical].append(gaussian_name)
    return names, dict(smiles)


_GAUSSIAN_NAMES, _SMILES_TO_GAUSSIAN_NAMES = _load_solvent_maps()


def resolve_gaussian_solvent(solvent):
    """Return a dielectric value or a Gaussian solvent name.

    Numeric inputs retain the existing Generic/Read behavior. Text names must
    occur in the ``Gaussian_Name`` column. SMILES are canonicalized before
    lookup so equivalent valid representations can be used.
    """
    if solvent is None:
        return None
    if isinstance(solvent, bool):
        raise ValueError("Solvent must be a dielectric constant, Gaussian name, or SMILES.")

    value = str(solvent).strip()
    if not value:
        raise ValueError("Solvent must not be empty.")

    try:
        numeric_value = float(value)
    except ValueError:
        numeric_value = None
    if numeric_value is not None:
        if not math.isfinite(numeric_value):
            raise ValueError("The solvent dielectric constant must be finite.")
        return value

    gaussian_name = _GAUSSIAN_NAMES.get(value.casefold())
    if gaussian_name is not None:
        return gaussian_name

    canonical = _canonical_smiles(value)
    matches = _SMILES_TO_GAUSSIAN_NAMES.get(canonical, []) if canonical else []
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        choices = ", ".join(matches)
        raise ValueError(
            f"Solvent SMILES {value!r} is ambiguous in the Gaussian solvent map "
            f"({choices}). Specify one of those Gaussian_Name values instead."
        )

    raise ValueError(
        f"Unsupported Gaussian solvent {value!r}. Specify a dielectric constant, "
        "or a name/SMILES present in solvent_map_gaussian_184.csv."
    )
