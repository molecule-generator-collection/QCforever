import sys
import gc
import numpy as np
from numpy import *

from qcforever import gaussian_run


def Make_xyzfile(atom, X, Y, Z, ofile, charge, spin):
    s = ""
    s += f"{len(atom)}\n"
    s += f"{charge} {spin}\n"
    with open(ofile,'w') as f:
        f.write(s)
        for i in range(len(atom)):
            f.write('%-4s % 10.5f % 10.5f % 10.5f \n' % (atom[i], X[i], Y[i], Z[i]))


def MaxBondLength(lines, bondpair1, bondpair2):
    _, Mol_X, Mol_Y, Mol_Z = Extract_Coordinate(lines, None)
    BondLength = zeros(len(bondpair1))
    BondYN = zeros(len(bondpair1))  # not used
    for i in range(len(bondpair1)):
        BondLength[i] = sqrt((Mol_X[int(bondpair2[i])-1] - Mol_X[int(bondpair1[i])-1])**2
                             + (Mol_Y[int(bondpair2[i])-1] - Mol_Y[int(bondpair1[i])-1])**2
                             + (Mol_Z[int(bondpair2[i])-1] - Mol_Z[int(bondpair1[i])-1])**2)
    return  amax(BondLength)


def Extract_Coordinate(lines, outfile=None):
    print ("Start finding coordinates...")
    Atom_index = []
    NumElement = []
    AtomicType = []
    X = []
    Y = []
    Z = []
    count = 0
    count_Station = 0
    Total_NumStation =0
    Index_Station = 0

    # Counting Stationary points
    for line in lines:
        if line.find("Charge =") >= 0:
            line_charge_spin = line.split()
            given_Charge = int(line_charge_spin[2]) 
            given_SpinMulti = int(line_charge_spin[-1])
        if line.find("-- Stationary point found.") >= 0:
            count_Station += 1

    Total_NumStation = count_Station
    if Total_NumStation > 0:
        count_Station = 0
        print ("Total number of stationary points: ", Total_NumStation)
        for line in lines:
            if line.find("-- Stationary point found.") >= 0:
                count_Station += 1
                Index_Station += 1
                print ("A stationary point is found! #", Index_Station)
                continue
            if line.find("Standard orientation:") >=0 and count_Station > 0:
                count += 1
                continue
            if count == 1:
                Border_1 = line  # not used
                count += 1
                continue
            if count == 2:
                Index_1 = line  # not used
                count += 1
                continue
            if count == 3:
                Index_2 = line  # not used
                count += 1
                continue
            if count == 4:
                Border_2 = line  # not used
                count += 1
                continue
            if count >= 5:
                i_atom = line.split()
                if len(i_atom) == 6:
                    Atom_index.append(int(i_atom[0]))
                    NumElement.append(int(i_atom[1]))
                    AtomicType.append(int(i_atom[2]))
                    X.append(float(i_atom[3]))
                    Y.append(float(i_atom[4]))
                    Z.append(float(i_atom[5]))
                    count += 1
                    continue
                else :
                    print ("Reading atom coordinates is finished...")
                    N = count-5
                    print ("Number of atoms: ", N)
                    count = 0
                    count_Station = 0
                continue
    else:
        for line in lines:
            if line.find("Standard orientation:") >=0:
                print ("Standard orientaion was found")
                print ("Start reading coordinate")
                count += 1
                continue
            if count == 1:
                Border_1 = line  # not used
                count += 1
                continue
            if count == 2:
                Index_1 = line  # not used
                count += 1
                continue
            if count == 3:
                Index_2 = line  # not used
                count += 1
                continue
            if count == 4:
                Border_2 = line  # not used
                count += 1
                continue
            if count >= 5:
                i_atom = line.split()
                if len(i_atom) == 6:
                    Atom_index.append(int(i_atom[0]))
                    NumElement.append(int(i_atom[1]))
                    AtomicType.append(int(i_atom[2]))
                    X.append(float(i_atom[3]))
                    Y.append(float(i_atom[4]))
                    Z.append(float(i_atom[5]))
                    count += 1
                    continue
                else :
                    print ("Reading atom coordinates is finished...")
                    N = count-5
                    print ("Number of atoms: ", N)
                    count = 0
                continue

    # Translating  atomic number to element symbol
    Mol_atom = []
    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    for i in range(N):
        Mol_atom.append(gaussian_run.AtomInfo.AtomicNumElec(NumElement[i]))

    if outfile != None:
        Make_xyzfile(Mol_atom, X, Y, Z, outfile , given_Charge, given_SpinMulti)

    return Mol_atom, X, Y, Z


if __name__ == '__main__':
    usage = 'Usage; %s infilename' % sys.argv[0]
    try:
        infilename = sys. argv[1]
    except:
        print (usage); sys.exit()
    with open(infilename, 'r') as ifile:
        lines = ifile.readlines()
    print(Extract_Coordinate(lines))
