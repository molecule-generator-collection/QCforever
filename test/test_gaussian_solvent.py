import pytest

from qcforever.gaussian_run.GaussianRunPack import GaussianDFTRun
from qcforever.gaussian_run.solvent import resolve_gaussian_solvent


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("water", "Water"),
        ("  DiMethylSulfoxide  ", "DiMethylSulfoxide"),
        ("O", "Water"),
        ("C(C)O", "Ethanol"),
        (78.3553, "78.3553"),
        ("0", "0"),
        (None, None),
    ],
)
def test_resolve_gaussian_solvent(value, expected):
    assert resolve_gaussian_solvent(value) == expected


def test_unknown_name_stops_before_gaussian():
    with pytest.raises(ValueError, match="Unsupported Gaussian solvent"):
        resolve_gaussian_solvent("NotAGaussianSolvent")


def test_ambiguous_smiles_requires_gaussian_name():
    with pytest.raises(ValueError, match="ambiguous"):
        resolve_gaussian_solvent("C(=CCl)Cl")


@pytest.mark.parametrize("value", ["", float("nan"), float("inf")])
def test_invalid_solvent_input(value):
    with pytest.raises(ValueError):
        resolve_gaussian_solvent(value)


def test_gaussian_run_defaults_to_vacuum():
    calculation = GaussianDFTRun("B3LYP", "STO-3G", 1, "energy", "dummy.xyz")
    assert calculation.solvent is None
    assert calculation.MakeSolventLine() == ("", "")
