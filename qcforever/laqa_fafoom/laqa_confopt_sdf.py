# -*- coding: utf-8 -*-

import sys
import time
import datetime

from rdkit import rdBase, Chem
from rdkit.Chem import AllChem

from qcforever import laqa_fafoom


def LAQA_confopt_main(param_file, insdffile):

    t_laqaopt_bgn = time.time()
    print(f"\nStart LAQA conforation searh job at {datetime.datetime.now()}")

    sdfmol = Chem.SDMolSupplier(insdffile, removeHs=False)
    mols = [x for x in sdfmol if x is not None]
    mol = mols[0]

    SMILES = Chem.MolToSmiles(mol)

    laqa_fafoom.initgeom.LAQA_initgeom(param_file, SMILES)
    laqa_fafoom.laqa_optgeom.LAQA_optgeom(param_file)
	
    t_laqaopt_end = time.time()
    print(f"\nFinish LAQA conforation searh job at {datetime.datetime.now()}")
    print(f"Wall time of LAQA conformation optimization job: {t_laqaopt_end - t_laqaopt_bgn:20.2f} sec.")


