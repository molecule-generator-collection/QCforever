import re
import sys
import numpy as np

def Extract_Freq(lines):
    count = 0
    Freq = []
    IR_int = []
    Raman_int = []

    print ("Get information about Frequencies")
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


def Extract_vibvec(lines):
    modes = []
    for i, line in enumerate(lines):
        fields = line.split()
        if fields[:2] != ["Atom", "AN"] or not fields[2:]:
            continue
        width, remainder = divmod(len(fields) - 2, 3)
        if remainder or fields[2:] != ["X", "Y", "Z"] * width:
            continue
        block = []
        for row in lines[i + 1:]:
            values = row.split()
            if len(values) != 2 + 3 * width or not all(v.isdigit() for v in values[:2]):
                break
            block.append([float(v) for v in values[2:]])
        if block:
            modes.extend(np.asarray(block).reshape(-1, width, 3).transpose(1, 0, 2))
    if not modes:
        raise ValueError("Missing Gaussian Cartesian normal modes")
    return tuple(np.asarray(modes).transpose(2, 0, 1))


def Extract_polar(lines):

    activate_extraction = False
    iso_values = []
    aniso_values = []

    for line in lines:
       # if line.strip() == "Dipole polarizability, Alpha (dipole orientation).":
       #     activate_extraction = True
        if line.strip() == "Dipole polarizability, Alpha (input orientation).":
            activate_extraction = True
        if activate_extraction:
            if line.strip().startswith("iso"):
                values = re.findall(r"\d+\.\d+[Dd][\+\-]?\d+", line)
                #iso_values.extend(values)
                iso_values.extend([float(val.replace('D', 'E')) for val in values])
            elif line.strip().startswith("aniso"):
                values = re.findall(r"\d+\.\d+[Dd][\+\-]?\d+", line)
                #aniso_values.extend(values)
                aniso_values.extend([float(val.replace('D', 'E')) for val in values])

    return iso_values, aniso_values


if __name__ == '__main__':
    usage = 'Usage; %s jobname' % sys.argv[0]
    try:
        infilename = sys. argv[1]
    except:
        print (usage); sys.exit()
    with open(infilename, 'r') as ifile:
        lines = ifile.readlines()
    print (Extract_Freq(lines))

    print (Extract_polar(lines))

    print (Extract_vibvec(lines))
