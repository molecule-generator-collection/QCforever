import io
import unittest
from contextlib import redirect_stdout
from unittest import mock

import numpy as np
from rdkit import Chem

from qcforever.laqa_fafoom import (
    get_parameters,
    laqa_optgeom,
    structure,
    utilities,
)
from qcforever.util import job_timeout


def _molecule(smiles="CC"):
    sdf = get_parameters.template_sdf(
        smiles, 1.3, 2.15, max_attempts=2, opt_steps=20, random_seed=7)
    return Chem.MolFromMolBlock(sdf, removeHs=False)


def _input_params(macro_opt_cycle=2):
    return {
        "inp_base": "initial_structures",
        "sdf_inp": "initial_structures.sdf",
        "struct_conv": 10,
        "macro_opt_cycle": macro_opt_cycle,
        "micro_opt_cycle": 5,
        "thr_econv": 5.0e-4,
        "thr_fconv": 1.0e-6,
        "printlevel": 0,
        "energy_function": "xtb",
        "gauss_exedir": "g16",
        "gauss_scrdir": ".",
        "qcmethod": "pm6",
        "xtb_call": "xtb",
        "gfn": "2",
        "charge": 0,
        "mult": 1,
        "solvmethod": None,
        "solvent": "water",
        "nprocs": 1,
        "memory": "1GB",
    }


class LAQAErrorRecoveryTest(unittest.TestCase):
    def test_invalid_optimized_geometry_is_not_written(self):
        mol_info = structure.MoleculeDescription(smiles="CC")
        mol_info.get_parameters()
        mol_info.create_template_sdf(random_seed=7)
        conformer = structure.Structure(mol_info)
        conformer.generate_structure()

        xtb_object = mock.Mock()
        xtb_object.get_energy.return_value = -1.0
        xtb_object.get_sdf_string_opt.return_value = conformer.sdf_string
        with mock.patch.object(
                structure.laqa_fafoom.pyxtb, "xTBObject",
                return_value=xtb_object), mock.patch.object(
                    conformer, "is_geometry_valid", return_value=False):
            with self.assertRaisesRegex(ValueError, "invalid interatomic"):
                conformer.perform_xtb("xtb")

        xtb_object.save_to_file.assert_not_called()
        self.assertEqual(xtb_object.clean.call_count, 2)

    def test_initial_failure_skips_only_the_bad_conformer(self):
        mol = _molecule()
        force = np.zeros((mol.GetNumAtoms(), 3))

        with mock.patch.object(
                laqa_optgeom.laqa_fafoom.pyxtb, "xtb_exec",
                side_effect=[FileNotFoundError("gradient"), (-1.0, force)]) as run:
            output = io.StringIO()
            with redirect_stdout(output):
                _, optimized = laqa_optgeom.LAQA_do_opt(
                    _input_params(macro_opt_cycle=0),
                    [Chem.Mol(mol), Chem.Mol(mol)],
                )

        self.assertEqual(run.call_count, 2)
        self.assertEqual(optimized, {})
        self.assertIn("Skipping structure ID 0", output.getvalue())

    def test_local_failure_drops_bad_conformer_and_continues(self):
        mol = _molecule()
        force = np.zeros((mol.GetNumAtoms(), 3))
        xyz = utilities.sdf2xyz(Chem.MolToMolBlock(mol))

        responses = [
            (-2.0, force),              # initial gradient, conformer 0
            (-1.0, force),              # initial gradient, conformer 1
            RuntimeError("xTB failed"), # local opt, conformer 0
            (-1.0, xyz),                # local opt, conformer 1
            (-1.0, force),              # new gradient, conformer 1
        ]
        with mock.patch.object(
                laqa_optgeom.laqa_fafoom.pyxtb, "xtb_exec",
                side_effect=responses) as run:
            output = io.StringIO()
            with redirect_stdout(output):
                _, optimized = laqa_optgeom.LAQA_do_opt(
                    _input_params(), [Chem.Mol(mol), Chem.Mol(mol)])

        self.assertEqual(run.call_count, 5)
        self.assertEqual(optimized, {2: -1.0})
        self.assertIn("Skipping structure ID 0", output.getvalue())

    def test_overall_timeout_is_not_treated_as_a_conformer_failure(self):
        mol = _molecule()
        with mock.patch.object(
                laqa_optgeom.laqa_fafoom.pyxtb, "xtb_exec",
                side_effect=job_timeout.QCforeverTimeoutError("deadline")):
            with self.assertRaises(job_timeout.QCforeverTimeoutError):
                laqa_optgeom.LAQA_do_opt(_input_params(), [mol])

    def test_rejects_partial_gradient(self):
        with self.assertRaisesRegex(ValueError, "expected 2"):
            laqa_optgeom.calc_mae_rms_max_force([[0.0, 0.0, 0.0]], 2)


if __name__ == "__main__":
    unittest.main()
