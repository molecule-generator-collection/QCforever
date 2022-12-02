import re
from numpy import *

import GaussianRunPack.Estimate_SpinContami


def classify_SpinStates(Ideal_SS, SS, WaveLength, OS, CD_L, CD_OS):
    State_allowed = []
    State_forbidden = []
    WL_allowed = []
    WL_forbidden = []
    OS_allowed = []
    OS_forbidden = []
    CD_L_allowed = []
    CD_L_forbidden = []
    CD_OS_allowed = []
    CD_OS_forbidden = []
    for i in range(len(SS)):
        if abs(SS[i]-Ideal_SS) <= 0.1:
            State_allowed.append(i+1)
            WL_allowed.append(WaveLength[i])                   
            OS_allowed.append(OS[i])
            CD_L_allowed.append(CD_L[i])
            CD_OS_allowed.append(CD_OS[i])
        else:
            State_forbidden.append(i+1)
            WL_forbidden.append(WaveLength[i])                   
            OS_forbidden.append(OS[i])
            CD_L_forbidden.append(CD_L[i])
            CD_OS_forbidden.append(CD_OS[i])
    return State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
            CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden 


def extract_CD(lines):
    count = 0
    CD_length = []
    CD_Osc = []
    for line in lines:
        if line.find("1/2[<0|r|b>*<b|rxdel|0> + (<0|rxdel|b>*<b|r|0>)*]") >= 0: 
            CD_length = []
            CD_Osc = []
            count += 1
            continue
        if count == 1:
            count += 1
            continue
        if count == 2:
            count += 1
            continue
        if count == 3:
            if line.find(" 1/2[<0|del|b>*<b|r|0> + (<0|r|b>*<b|del|0>)*] (Au)") >= 0: 
                count += 1
                continue
            else:
                info_R =  line.split()
                cd_value = float(info_R[-1])
                CD_length.append(cd_value)    
                continue
        if count == 4:
            count += 1 
            continue
        if count == 5:
            if line == "": 
                count = 0
                continue
            else:
                info_Osc =  line.split()
                osc_value = float(info_Osc[-1])
                CD_Osc.append(osc_value)    
                continue
    return CD_length, CD_Osc


def Extract_ExcitedState(lines):
    Egrd = 0.0
    Eext = 0.0
    Found = False
    WaveLength = []
    V_OS = []
    SS = []
    _, Ideal_SS = GaussianRunPack.Estimate_SpinContami.Estimate_SpinDiff(lines)

    print ("Get information about excited state")
    for line in lines:
        if line.find("SCF Done:  ") >=0:
            line_SCFEnergy = re.split("\s+", line)
            Egrd = float(line_SCFEnergy[5])
        if line.find("Total Energy, E(TD-HF/TD-DFT)") >=0:
            line_totalenergy = line.split('=')
            Eext = float(line_totalenergy[1])
        if line.find("Excitation energies and oscillator strengths:") >=0:
             WaveLength = []
             V_OS = []
             SS = []
        if line.find("Excited State  ") >=0:
             line_StateInfo = line.split()
             WaveLength.append(float(line_StateInfo[6]))
             OS_info = line_StateInfo[8].split('=')
             V_OS.append(float(OS_info[1]))
             SS_info = line_StateInfo[9].split('=')
             SS.append(float(SS_info[1]))
        if line.find("-- Stationary point found.") >=0:
            Found = True

    CD_length, CD_OS = extract_CD(lines)
    State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
        CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden = classify_SpinStates(Ideal_SS, SS, WaveLength, V_OS, CD_length, CD_OS)
    return Found, Egrd, Eext, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
        CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden 
