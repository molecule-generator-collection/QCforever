# -*- coding: utf-8 -*-

import sys
import time
import datetime
import shutil

from rdkit import rdBase, Chem
from rdkit.Chem import AllChem, rdDetermineBonds, rdMolDescriptors

from qcforever import laqa_fafoom


def LAQA_confopt_main(infilename, TotalCharge, SpinMulti, method, nproc, mem):

    PreInput = infilename.split('.')

    t_laqaopt_bgn = time.time()
    print(f"\nStart LAQA conforation searh job at {datetime.datetime.now()}")
    if PreInput[-1] == "sdf":
        sdfmol = Chem.SDMolSupplier(infilename, removeHs=False)
        mols = [x for x in sdfmol if x is not None]
        mol = mols[0]
    elif PreInput[-1] == "xyz":
        mol = Chem.MolFromXYZFile(infilename)
        rdDetermineBonds.DetermineBonds(mol)
    else:
        exit()

    #get the number of rotatable bonds
    RotBond = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)

    SMILES = Chem.MolToSmiles(mol)

    make_laqa_input(SMILES, SpinMulti, TotalCharge, RotBond, method, nproc, mem)

    laqa_fafoom.initgeom.LAQA_initgeom('laqa_setting.inp', SMILES)
    laqa_fafoom.laqa_optgeom.LAQA_optgeom('laqa_setting.inp')
	
    t_laqaopt_end = time.time()
    print(f"\nFinish LAQA conforation searh job at {datetime.datetime.now()}")
    print(f"Wall time of LAQA conformation optimization job: {t_laqaopt_end - t_laqaopt_bgn:20.2f} sec.")


def make_laqa_input(SMILES, SpinMulti, TotalCharge, RotBond, method, nproc, mem):

    Num_popsize = 3*RotBond

    input_s = ''

    input_s += '[Molecule]\n'
    input_s += f'\nsmiles = "{SMILES}"\n'

    input_s += f'\n[Initial geometry]\n'
    input_s += f'popsize = {Num_popsize}\n'

    input_s += '\n[LAQA settings]\n'

    input_s += f'\ncharge = {TotalCharge}\n'
    input_s += f'mult = {SpinMulti}\n'


    if method == 'xtb':
        input_s += f'\nenergy_function = "{method}"\n'
        input_s += f'gfn = "2"\n'
    if method == 'pm6':
        input_s += '\nenergy_function = "g16"\n'
        #Get the directry that include 'g16'
        Gaussian_exedir = shutil.which('g16')
        #Delete the final binary file name of 'g16'
        input_s += f'gauss_exedir  = "{Gaussian_exedir[:-4]}"\n'
        input_s += f'qcmethod = "{method}" \n'
        if nproc > 1:
            input_s += f'nprocs = "{nproc}" \n'
        if mem != '':
            input_s += f'memory = "{mem}" \n'

    with open('laqa_setting.inp', 'w') as laqa_infile:
        
            laqa_infile.write(input_s)


