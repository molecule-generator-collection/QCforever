import glob
import os
import re

import GaussianRunPack.AtomInfo


Bohr2Ang = 0.529177


def Make_xyzfile(atom, X, Y, Z, ofile, charge, spin):
    s = ""
    s += str(len(atom)) + "\n"
    s += str(charge)+" "+ str(spin)+"\n"

    with open(ofile, 'w') as f:
        f.write(s)
        for i in range(len(atom)):
            f.write('%-4s % 10.5f % 10.5f % 10.5f \n' % (atom[i], X[i], Y[i], Z[i]))


def Get_fchklist2xyz():
    for f in glob.glob('./*.fchk'):
        print(f)
        filename = re.split('[/.]', f)
        xyzfile = f"{filename[-2]}.xyz"
        print()
        with open(f, 'r') as ifile:
            lines = ifile.readlines()
        Extract_MolCoord(lines, xyzfile)
        os.remove(os.path.join('.', f))


def Extract_MolCoord(lines, out = None):
        count = 0
        Charge = 0
        SpinMulti = 0
        NuclearCharge = []
        MolCoord = []
        for line in lines:
            if line.find("Charge ") >= 0:
                line_ChargeInfo = line.split()
                Charge = int(line_ChargeInfo[-1])
                continue
            if line.find("Multiplicity ") >= 0:
                line_SpinMultiInfo = line.split()
                SpinMulti = int(line_SpinMultiInfo[-1])
                continue
            if line.find("Nuclear charges ") >= 0:
                line_NuclearCharges = line.split()
                NumAtom = int(line_NuclearCharges[-1])
                count += 1
                continue
            if count == 1:
                NuclearCharge_info = line.split()
                for i in range(len(NuclearCharge_info)):
                    NuclearCharge.append(float(NuclearCharge_info[i]))
                if len(NuclearCharge) == NumAtom:
                    count = 0
            if line.find("Current cartesian coordinates ") >= 0:
                line_Numdata = line.split()
                Numdata = int(line_Numdata[-1])
                count += -1
                continue
            if count == -1:
                Coordinate_info = line.split()
                for i in range(len(Coordinate_info)):
                    MolCoord.append(float(Coordinate_info[i]))
                if len(MolCoord) == Numdata: 
                    count = 0

        Mol_atom = []
        Mol_X = []
        Mol_Y = []
        Mol_Z = []
        for i in range(NumAtom):
            Mol_atom.append(GaussianRunPack.AtomInfo.AtomicNumElec(int(NuclearCharge[i])))
            # Mol_atom.append(AtomInfo.AtomicNumElec(int(NuclearCharge[i])))
        for i in range(0, len(MolCoord)-1,3):
            Mol_X.append(MolCoord[i]*Bohr2Ang)
            Mol_Y.append(MolCoord[i+1]*Bohr2Ang)
            Mol_Z.append(MolCoord[i+2]*Bohr2Ang)
        if out !=None:
            Make_xyzfile(Mol_atom, Mol_X, Mol_Y, Mol_Z, out, Charge, SpinMulti)
        return Charge, SpinMulti, Mol_atom, Mol_X, Mol_Y, Mol_Z

if __name__ == '__main__':
    Get_fchklist2xyz()      
