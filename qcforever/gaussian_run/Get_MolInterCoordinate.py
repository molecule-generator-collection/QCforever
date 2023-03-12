import sys
import re
import numpy as np


def Extract_InterMol(lines, outfile=None):
    print ("Start finding internal coordinates...")
    count_Station = 0
    ifrag = 0
    icount = 0
    fcount = 0
    # Counting Stationary points
    for line in lines:
        if line.find("Charge =") >= 0:
            line_charge_spin = line.split()
            given_Charge = int(line_charge_spin[2])  # not used
            given_SpinMulti = int(line_charge_spin[-1])  # not used
        if line.find("-- Stationary point found.") >= 0:
            count_Station += 1

    Total_NumStation = count_Station
    if Total_NumStation > 0:
        print ("Total number of stationary points: ", Total_NumStation)
        for line in lines:
            if line.find("!    Initial Parameters    !") >= 0 and ifrag == 0 :
                Init_BondLength    = []
                Init_BondAngle     = []
                Init_DihedralAngle = []
                BL = []
                BA = []
                DA = []
                ifrag += 1
                icount += 1
                print ("Initial internal coordinate parameters were found.")
                continue
            if icount >= 1 and icount <= 4:
                icount += 1
                continue
            if icount >= 5:
                i_intermol = line.split()
                icount += 1
                if len(i_intermol) > 2:
                    if re.search('R',i_intermol[1]):
                        Init_BondLength.append(float(i_intermol[3]))
                        BL.append(re.sub(r'[A-Z()]', '', i_intermol[2]))
                    if re.search('A',i_intermol[1]):
                        Init_BondAngle.append(float(i_intermol[3]))
                        BA.append(re.sub(r'[A-Z()]', '', i_intermol[2]))
                    if re.search('D',i_intermol[1]):
                        tmp_D = float(i_intermol[3])
                        DA.append(re.sub(r'[A-Z()]', '', i_intermol[2]))
                        if abs(tmp_D)/180 < 1:
                            Init_DihedralAngle.append(tmp_D)
                        else:
                            if tmp_D > 0:
                                Init_DihedralAngle.append(180-(abs(tmp_D)-180)*(-1))
                            else:
                                Init_DihedralAngle.append(180-(abs(tmp_D)-180))
                elif len(i_intermol) < 2  :
                    print ("Reading the initial internal coordinates is finished...")
                    N = icount-6
                    print ("Number of internal coordinate: ", N)
                    icount = 0
                continue
            if line.find("!   Optimized Parameters   !") >= 0:
                fcount += 1
                print ("Final internal coordinate parameters were found.")
                Final_BondLength    = []
                Final_BondAngle     = []
                Final_DihedralAngle = []
                continue
            if fcount >= 1 and fcount <= 4:
                fcount += 1
                continue
            if fcount >= 5:
                i_intermol = line.split()
                fcount += 1
                if len(i_intermol) > 2:
                    if re.search('R',i_intermol[1]):
                        Final_BondLength.append(float(i_intermol[3]))
                    if re.search('A',i_intermol[1]):
                        Final_BondAngle.append(float(i_intermol[3]))
                    if re.search('D',i_intermol[1]):
                        tmp_D = float(i_intermol[3])
                        if abs(tmp_D)/180 < 1:
                            Final_DihedralAngle.append(tmp_D)
                        else:
                            if tmp_D > 0:
                                Final_DihedralAngle.append(180-(abs(tmp_D)-180)*(-1))
                            else:
                                Final_DihedralAngle.append(180-(abs(tmp_D)-180))
                elif len(i_intermol) < 2:
                    print ("Reading the final internal coordinates is finished...")
                    N = fcount-6
                    print ("Number of internal coordinate: ", N)
                    fcount = 0
                continue
    else:
        print("Stational points are not found!")

    MaxDisplacementB = 0.0
    MaxDisplacementA = 0.0
    MaxDisplacementD = 0.0

    if BL != []:
        Init_BondLength     = np.array(Init_BondLength)
        Final_BondLength    = np.array(Final_BondLength)
        DisplacementB = Final_BondLength -  Init_BondLength
        MaxDisplacementB = DisplacementB.max()
        iMaxB = BL[DisplacementB.argmax()]

    if BA != []:
        Init_BondAngle      = np.array(Init_BondAngle)
        Final_BondAngle     = np.array(Final_BondAngle)
        DisplacementA = Final_BondAngle -  Init_BondAngle
        MaxDisplacementA = DisplacementA.max()
        iMaxA = BA[DisplacementA.argmax()]

    if DA != []:
        Init_DihedralAngle  = np.array(Init_DihedralAngle)
        Final_DihedralAngle = np.array(Final_DihedralAngle)
        DisplacementD = Final_DihedralAngle -  Init_DihedralAngle

        for i in range(len(DisplacementD)):
            if abs(DisplacementD[i]) > 180:
                DisplacementD[i] = (360-abs(DisplacementD[i]))

        MaxDisplacementD = DisplacementD.max()
        iMaxD = DA[DisplacementD.argmax()]

    if outfile != None:
        out_line =''
        for i in range (len(BL)):
            out_line += f'{BL[i]:10s} {Init_BondLength[i]: 10.5f} {Final_BondLength[i]: 10.5f}\n'
        for i in range (len(BA)):
            out_line += f'{BA[i]:10s} {Init_BondAngle[i]: 10.5f} {Final_BondAngle[i]: 10.5f}\n'
        for i in range (len(DA)):
            out_line += f'{DA[i]:10s} {Init_DihedralAngle[i]: 10.5f} {Final_DihedralAngle[i]: 10.5f}\n'
        with open(outfile, 'w') as f:
            f.write(out_line)

    return [MaxDisplacementB, MaxDisplacementA, MaxDisplacementD]


if __name__ == '__main__':
    usage = 'Usage; %s infilename' % sys.argv[0]
    try:
        infilename = sys. argv[1]
    except:
        print (usage); sys.exit()
    with open(infilename, 'r') as ifile:
        lines = ifile.readlines()
    print(Extract_InterMol(lines,'internalcoord.txt'))
