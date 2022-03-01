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
         

if __name__ == '__main__':
        
    usage = 'Usage; %s jobname' % sys.argv[0]

    try:
        infilename = sys. argv[1]
    except:
        print (usage); sys.exit()


    with open(infilename, 'r') as ifile:
        lines = ifile.readlines()

       
    print (Extract_ChargeSpin(lines))

