"""Preserve the chosen electronic state across the optconf SDF handoff."""

from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from qcforever import gamess_run, gaussian_run
from qcforever.conformer_search import conformer_search as workflow
from qcforever.gamess_run.GamessRunPack import GamessDFTRun
from qcforever.gaussian_run.GaussianRunPack import GaussianDFTRun
from qcforever.util import check_resource, read_mol_file


@pytest.fixture(params=["gaussian", "gamess"])
def backend(request, monkeypatch):
    monkeypatch.setattr(check_resource, "respec_cores", lambda count: count)
    monkeypatch.setattr(check_resource, "respec_memory", lambda **kwargs: "1GB")
    if request.param == "gaussian":
        monkeypatch.setattr(gaussian_run.Exe_Gaussian, "exe_Gaussian", lambda *args: "normal")
        monkeypatch.setattr(gaussian_run.chk2fchk, "Get_chklist", lambda *args: None)
        run_class, method = GaussianDFTRun, "run_gaussian"
    else:
        monkeypatch.setattr(gamess_run.Exe_Gamess, "exe_Gamess", lambda *args: "")
        # Cleanup of the user's GAMESS scratch directories is unrelated to this test.
        monkeypatch.setattr(gamess_run.GamessRunPack.job_cleanup, "cleanup_gamess", lambda *args, **kwargs: None)
        run_class, method = GamessDFTRun, "run_gamess"
    monkeypatch.setattr(run_class, "Extract_values", lambda *args: {})
    return request.param, run_class, method


def run_handoff(tmp_path, monkeypatch, backend, *, charge=np.nan, multiplicity=np.nan,
                smiles="C", xyz_state=None, tune=False, scan_spin=None, search_fails=False):
    """Run real file readers and input writers; replace only native computation."""
    name, run_class, method = backend
    monkeypatch.chdir(tmp_path)
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=22) == 0
    selected = Chem.Mol(mol)
    for index, position in enumerate(selected.GetConformer().GetPositions()):
        selected.GetConformer().SetAtomPosition(index, position + [4, 2, -3])
    if xyz_state is None:
        input_file = "molecule.sdf"
        with Chem.SDWriter(input_file) as writer:
            writer.write(mol)
    else:
        input_file = "molecule.xyz"
        lines = Chem.MolToXYZBlock(mol).splitlines()
        lines[1] = "%s %s" % xyz_state
        Path(input_file).write_text("\n".join(lines) + "\n")

    searched, tuned = [], []

    def search(filename, total_charge, spin_multi, *args, **kwargs):
        searched.append((total_charge, spin_multi))
        if search_fails:
            return {"state": "failed", "error": "fixture search failure"}
        with Chem.SDWriter("optimized_structures.sdf") as writer:
            writer.write(selected)
        return {"state": "succeeded"}

    monkeypatch.setattr(workflow, "run_conformer_search", search)
    option = "energy optconf=xtb" + (" optspin" if scan_spin is not None else "")
    functional = ("ktlc-blyp-bo" if name == "gaussian" else "KTLC-BLYP") if tune else "B3LYP"
    runner = run_class(functional, "STO-3G", 1, option, input_file)
    runner.SpecTotalCharge = charge
    runner.SpecSpinMulti = multiplicity
    if scan_spin is not None:
        monkeypatch.setattr(runner, "SpinMulti_scan", lambda **kwargs: scan_spin)
    if tune:
        def optimize(*args):
            tuned.append(tuple(args[-2:]))
            return []
        monkeypatch.setattr(runner, "LC_para_BOopt", optimize)

    result = getattr(runner, method)()
    if name == "gaussian":
        inputs = list((tmp_path / "molecule").glob("*.com"))
        assert len(inputs) == 1
    else:
        inputs = [tmp_path / "molecule/molecule.inp"]
    return result, searched, tuned, inputs[0].read_text()


def assert_input_state(backend, content, expected):
    if backend[0] == "gaussian":
        assert list(map(str, expected)) in [line.split() for line in content.splitlines()]
    else:
        assert "ICHARG=%s MULT=%s" % expected in content


@pytest.mark.parametrize("tune", [False, True])
@pytest.mark.parametrize("charge,multiplicity,smiles,expected", [
    (-1, 2, "C", (-1, 2)),
    (-2, np.nan, "C", (-2, 1)),
    (np.nan, 3, "C", (0, 3)),
    (0, 2, "[NH4+]", (0, 2)),
    (np.nan, np.nan, "C", (0, 1)),
])
def test_optconf_preserves_explicit_or_input_electronic_state(
        tmp_path, monkeypatch, backend, charge, multiplicity, smiles, expected, tune):
    result, searched, tuned, content = run_handoff(
        tmp_path, monkeypatch, backend, charge=charge, multiplicity=multiplicity,
        smiles=smiles, tune=tune)
    assert searched == [expected]
    assert tuned == ([expected] if tune else [])
    assert_input_state(backend, content, expected)
    assert result["optconf"] is True
    # The SDF's own formal charge/radicals remain unchanged: it supplies geometry.
    selected_state = read_mol_file.read_sdf(str(tmp_path / "molecule/optimized_structures.sdf"))[4:6]
    assert selected_state == ((1, 1) if smiles == "[NH4+]" else (0, 1))


@pytest.mark.parametrize("tune", [False, True])
def test_optconf_preserves_xyz_comment_state(tmp_path, monkeypatch, backend, tune):
    result, searched, tuned, content = run_handoff(
        tmp_path, monkeypatch, backend, xyz_state=(-1, 2), tune=tune)
    assert searched == [(-1, 2)]
    assert tuned == ([(-1, 2)] if tune else [])
    assert_input_state(backend, content, (-1, 2))
    assert result["optconf"] is True


@pytest.mark.parametrize("backend", ["gaussian"], indirect=True)
def test_gaussian_optconf_preserves_spin_scan_selection(tmp_path, monkeypatch, backend):
    result, searched, tuned, content = run_handoff(
        tmp_path, monkeypatch, backend, multiplicity=1, scan_spin=3, tune=True)
    assert searched == tuned == [(0, 3)]
    assert_input_state(backend, content, (0, 3))
    assert result["spinmulti"] == 3


def test_failed_optconf_keeps_electronic_state(tmp_path, monkeypatch, backend):
    result, searched, tuned, content = run_handoff(
        tmp_path, monkeypatch, backend, charge=-1, multiplicity=2,
        tune=True, search_fails=True)
    assert searched == tuned == [(-1, 2)]
    assert_input_state(backend, content, (-1, 2))
    assert result["optconf"] is False
