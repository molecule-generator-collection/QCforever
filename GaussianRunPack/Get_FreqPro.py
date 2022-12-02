import sys
from numpy import *


def Extract_Freq(lines):
    count = 0
    Freq = []
    IR_int = []
    Raman_int = []

    print ("Get information about  Frequencies")
    for line in lines:
        if line.find("Frequencies --") >= 0:
            line_Freq = line.split()
            for i in range(2, len(line_Freq)):
                Freq.append(float(line_Freq[i]))
        if line.find("IR Inten    --") >= 0:
            line_ir = line.split()
            for i in range(3, len(line_ir)):
                IR_int.append(float(line_ir[i]))
        if line.find("Raman Activ --") >= 0:
            line_raman = line.split()
            for i in range(3, len(line_raman)):
                Raman_int.append(float(line_raman[i]))
        if line.find("Sum of electronic and zero-point Energies=") >= 0:
            line_zp = line.split()
            E_ZP = float(line_zp[-1])
        if line.find("Sum of electronic and thermal Energies=") >= 0:
            line_t = line.split()
            E_t = float(line_t[-1])
        if line.find("Sum of electronic and thermal Enthalpies=") >= 0:
            line_t = line.split()
            E_enth = float(line_t[-1])
        if line.find("Sum of electronic and thermal Free Energies=") >= 0:
            line_free = line.split()
            E_free = float(line_free[-1])
        if line.find("KCal/Mol") >= 0:
            count += 1
            continue
        if count == 1:
            line_ECvS = line.split()
            Ei = line_ECvS[1]
            Cv = line_ECvS[2]
            St= line_ECvS[3]
            count = 0
    return  Freq, IR_int, Raman_int, E_ZP,  E_t, E_enth, E_free, Ei, Cv, St


if __name__ == '__main__':
    usage = 'Usage; %s jobname' % sys.argv[0]
    try:
        infilename = sys. argv[1]
    except:
        print (usage); sys.exit()
    with open(infilename, 'r') as ifile:
        lines = ifile.readlines()
    print (Extract_Freq(lines))
