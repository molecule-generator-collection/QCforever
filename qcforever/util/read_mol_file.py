import sys

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdmolops


def read_sdf(infile):
    with open(infile, 'r') as ifile:
        sdfinfo = Chem.SDMolSupplier(infile, removeHs=False)
        for mol in sdfinfo:
            N = mol.GetNumAtoms()
            N_Bond = mol.GetNumBonds()
            N_radical = Descriptors.NumRadicalElectrons(mol)
            Charge = rdmolops.GetFormalCharge(mol) 
        
        count = 0
        X = []
        Y = []
        Z = []
        element_symbol = []
        Bond_pair1 = []
        Bond_pair2 = []
        Bond_type = []
        sdfCharge = 0 
        CHG_atom = []
        CHG = []
        sdfSpinMulti = 0
        RAD_atom =[]
        RAD = []
        
        for line in ifile:
            if count == 0:
                Header1 = line  # not used
                count += 1
                continue
            if count == 1:
                Header2 = line  # not used
                count += 1
                continue
            if count == 2:
                Header3 = line  # not used
                count += 1
                continue
            if count == 3:
                a = line.split()  # not used
                count += 1
                continue
            if 3 < count <= N+4:
                i_atom = line.split()
                if len(i_atom) != 0:
                    X.append(float(i_atom[0]))
                    Y.append(float(i_atom[1]))
                    Z.append(float(i_atom[2]))
                    element_symbol.append(i_atom[3])
                count += 1
                continue
            if N+4 < count <= N+N_Bond+3:
                bond_info = line.split()
                bond_info = line.split()
                Bond_pair1.append(int(bond_info[0]))
                Bond_pair2.append(int(bond_info[1]))
                Bond_type.append(int(bond_info[2]))
                count += 1
                continue
            if count > N+N_Bond+3:
                mol_info = line.split()
                if (mol_info[0] == "M"):
                    if (mol_info[1] == "END"):
                        break
                    if (mol_info[1] == "CHG"):
                        Num_CHGInfo = int(mol_info[2])
                        for k in range(Num_CHGInfo):
                            CHG_atom.append(int(mol_info[3+2*k]))
                            CHG.append(int(mol_info[4+2*k]))
                            sdfCharge += int(mol_info[4+2*k])
                    if (mol_info[1] == "RAD"):
                        Num_RADInfo = int(mol_info[2])
                        for l in range(Num_RADInfo):
                            RAD_atom.append(int(mol_info[3+2*l]))
                            RAD.append(int(mol_info[4+2*l]))
                            sdfSpinMulti += int(mol_info[4+2*l])
                else:
                    print("The sdf file is invalid!")
                    sys.exit()
                count +=1
        
    # Copy to array of numpy
    Mol_atom = []
    Mol_CartX = np.zeros(N)
    Mol_CartY = np.zeros(N)
    Mol_CartZ = np.zeros(N)
    CHG_atom = np.array(CHG_atom)
    CHG = np.array(CHG)
    for j in range(N):
        Mol_CartX[j] = X[j] 
        Mol_CartY[j] = Y[j] 
        Mol_CartZ[j] = Z[j] 
        Mol_atom.append(element_symbol[j])

    print('Reading the sdf file has finished')
    SpinMulti = N_radical + 1
    return Mol_atom, Mol_CartX, Mol_CartY, Mol_CartZ, Charge, SpinMulti, Bond_pair1, Bond_pair2


def read_xyz(infile):
    with open(infile, 'r') as ifile:
        count = 0
        X = []
        Y = []
        Z = []
        element_symbol = []
        Charge = 0 
        SpinMulti = 0
        for line in ifile:
        #    print(line)
            if count == 0:
                Atom_num = line.split()
                N = int(Atom_num[0])
                count += 1
                continue
            if count == 1:
                Comment = line.split()
                Charge = int(Comment[0])
                SpinMulti = int (Comment[1])
                count += 1
                continue
            if 1 < count:
                i_atom = line.split()
                if len(i_atom) > 3:
                    element_symbol.append(i_atom[0])
                    X.append(float(i_atom[1]))
                    Y.append(float(i_atom[2]))
                    Z.append(float(i_atom[3]))
                    count += 1
                    continue
                else:
                    count = 0
                    break
    # Copy to array of numpy
    Mol_atom = []
    Mol_CartX = np.zeros(N)
    Mol_CartY = np.zeros(N)
    Mol_CartZ = np.zeros(N)
    for j in range(N):
        Mol_CartX[j] = X[j] 
        Mol_CartY[j] = Y[j] 
        Mol_CartZ[j] = Z[j] 
        Mol_atom.append(element_symbol[j])

    print('Reading the xyz file has finished')
    return Mol_atom, Mol_CartX, Mol_CartY, Mol_CartZ, Charge, SpinMulti
