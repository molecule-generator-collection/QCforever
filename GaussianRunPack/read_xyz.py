import gc
import sys
import numpy as np


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
            print(line)
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
    # del element_symbol[N:TotalStep]
    del element_symbol[:]
    del X[:]
    del Y[:]
    del Z[:]
    gc.collect()

    print('Reading the xyz file has finished')
    return Mol_atom, Mol_CartX, Mol_CartY, Mol_CartZ, Charge, SpinMulti


if __name__ == '__main__':
    usage = 'Usage; %s sdf_file' % sys.argv[0]
    try:
        infilename = sys. argv[1]
    except:
        print (usage); sys.exit()
    print (read_xyz(infilename))
