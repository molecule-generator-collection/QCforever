import os, sys, math
import re, gc
from numpy import *


def Extract_ChargeSpin(lines):

    print("Start finding charge & spin...")

    count = 0

    for line in lines:
        if line.find("Mulliken charges and spin densities:") >= 0:
            Atom_index = []
            Charge = []
            Spin = []
            count += 1
            continue

        if count == 1:
            count += 1
            continue

        if count == 2:
            if line.find("Sum of Mulliken charges = ") >= 0:
                count = 0 

            else:
#                print(line)
                Charge_spin_line = line.split()
                Atom_index.append(Charge_spin_line[1])
                Charge.append(Charge_spin_line[2])
                Spin.append(Charge_spin_line[3])
                continue


    return Atom_index, Charge, Spin


def find_MaxProtic(Atom, Charge):

    HCharge = []
    H_index = []

    for i in range(len(Atom)):
        if Atom[i] == 'H':
            HCharge.append(Charge[i])
            H_index.append(i)
 
    if H_index != []:
        Max_HCharge = max(HCharge)
        Max_HCharge_index = HCharge.index(Max_HCharge)

        print (Max_HCharge)
    
        return H_index[Max_HCharge_index]

    else:

        return None

         

if __name__ == '__main__':
        
    usage = 'Usage; %s jobname' % sys.argv[0]

    try:
        infilename = sys. argv[1]
    except:
        print (usage); sys.exit()


    with open(infilename, 'r') as ifile:
        lines = ifile.readlines()

    Atom, Charge, Spin = Extract_ChargeSpin(lines)
    
    print(Atom, Charge)
       

    print(find_MaxProtic(Atom, Charge))
