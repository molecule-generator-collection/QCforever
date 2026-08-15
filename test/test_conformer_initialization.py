import os
import tempfile
import unittest
from unittest import mock

from rdkit import Chem

from qcforever.laqa_fafoom import get_parameters
from qcforever.laqa_fafoom.laqa_confopt_QCforever import make_laqa_input


class ConformerInitializationTest(unittest.TestCase):
    def test_template_sdf_generates_valid_3d_geometry_without_writing_file(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            old_cwd = os.getcwd()
            os.chdir(tmp_dir)
            try:
                sdf = get_parameters.template_sdf(
                    "CCCC", 1.3, 2.15, max_attempts=2, opt_steps=50,
                    random_seed=7,
                )
            finally:
                os.chdir(old_cwd)

            mol = Chem.MolFromMolBlock(sdf, removeHs=False)
            self.assertIsNotNone(mol)
            self.assertEqual(mol.GetNumConformers(), 1)
            self.assertFalse(os.path.exists(os.path.join(tmp_dir, "mol.sdf")))

    def test_template_sdf_stops_after_configured_embedding_attempts(self):
        with mock.patch.object(
                get_parameters.AllChem, "EmbedMolecule", return_value=-1) as embed:
            with self.assertRaisesRegex(RuntimeError, "after 3 attempts"):
                get_parameters.template_sdf("CC", 1.3, 2.15, max_attempts=3)

        self.assertEqual(embed.call_count, 3)

    def test_generated_population_size_is_bounded(self):
        for rotatable_bonds, expected in [(0, 1), (2, 6), (20, 30)]:
            with self.subTest(rotatable_bonds=rotatable_bonds):
                with tempfile.TemporaryDirectory() as tmp_dir:
                    old_cwd = os.getcwd()
                    os.chdir(tmp_dir)
                    try:
                        make_laqa_input(
                            "CC", 1, 0, rotatable_bonds, "xtb", 1, "")
                        with open("laqa_setting.inp") as settings_file:
                            settings = settings_file.read()
                    finally:
                        os.chdir(old_cwd)

                self.assertIn(f"popsize = {expected}", settings)


if __name__ == "__main__":
    unittest.main()
