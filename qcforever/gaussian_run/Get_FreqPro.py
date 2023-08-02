import re
import sys
import numpy as np

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


def Extract_vibvec(lines):
    result_lines = []
    vibX = []
    vibY = []
    vibZ = []

    count = 0

    for line in lines:
        pattern_match = re.match(r'\s*\bAtom\s+AN\s+X\s+Y\s+Z\s+X\s+Y\s+Z\s+X\s+Y\s+Z\b', line)
        if pattern_match:
            count += 1
        match = re.match(r'\s*(\d+)\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+', line)
        if match:
            vibvec_info = line.strip().split()
            #print(vibvec_info[2:])
            for i in range(2, len(vibvec_info), 3): 
                vibX.append(float(vibvec_info[i]))
                vibY.append(float(vibvec_info[i+1]))
                vibZ.append(float(vibvec_info[i+2]))

            result_lines.append(line)

    vibX = np.array(vibX)
    vibY = np.array(vibY)
    vibZ = np.array(vibZ)

    vibX_tmp = vibX.reshape(count, -1)
    vibY_tmp = vibY.reshape(count, -1)
    vibZ_tmp = vibZ.reshape(count, -1)

    vibX_reshape = vibX_tmp[0].reshape(-1, 3).T
    vibY_reshape = vibY_tmp[0].reshape(-1, 3).T
    vibZ_reshape = vibZ_tmp[0].reshape(-1, 3).T

    for i in range(1, count):
        vibX_reshape = np.vstack([vibX_reshape, vibX_tmp[i].reshape(-1, 3).T])
        vibY_reshape = np.vstack([vibY_reshape, vibY_tmp[i].reshape(-1, 3).T])
        vibZ_reshape = np.vstack([vibZ_reshape, vibZ_tmp[i].reshape(-1, 3).T])

#    print(vibX_reshape)
#    print(vibY_reshape)
#    print(vibZ_reshape)

    return vibX_reshape, vibY_reshape, vibZ_reshape

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
