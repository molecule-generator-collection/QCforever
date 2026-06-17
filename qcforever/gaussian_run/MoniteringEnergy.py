import re
import numpy as np

eV2Eh = 0.03674932

class Monitor_log:

    def __init__(self, outfile):
        with open (outfile) as output_line:
            self.lines = output_line.readlines()

    def Extract_ConversionProcess(self):
        Egrd = 0.0
        Etarget = 0.0
        EeV = 0.0
        Found = False
        Es = []
        Eprofile = []
        OSprofile = []
        Index_TState = 0 #Index of the target state for optimization, default is the ground state(0)

        print ("Get information of convergence energy profile")
        for line in self.lines:
            if line.find("SCF Done:  ") >=0:
                Es = []
                V_OS = []
                SS = []
                line_SCFEnergy = re.split(r"\s+", line)
                Egrd = float(line_SCFEnergy[5])
                Es.append(Egrd)
            if line.find("Total Energy, E(TD-HF/TD-DFT)") >=0:
                line_totalenergy = line.split('=')
                Etarget = float(line_totalenergy[1])
                Index_TState = len(Es)-1
#                print('Target state: ', Index_TState)
            if line.find("Excited State  ") >= 0:
                line_StateInfo = line.split()
                EeV = float(line_StateInfo[4])
                Es.append(EeV*eV2Eh + Egrd)
                OS_info = line_StateInfo[8].split('=')
                V_OS.append(float(OS_info[1]))
                SS_info = line_StateInfo[9].split('=')
                SS.append(float(SS_info[1]))
            if line.find("Cartesian Forces:") >= 0:
                Eprofile.append(Es)
                OSprofile.append(V_OS)

        return Index_TState, Eprofile, OSprofile


    def Extract_ConversionCriterion(self):

        MaxForce = []
        RMSForce = []
        MaxDisp = []
        RMSDisp = []

        MaxForce_threshold = []
        RMSForce_threshold = []
        MaxDisp_threshold = []
        RMSDisp_threshold = []
        
        print ("Get information of convergence force & displacement profile")
        for line in self.lines:
            if line.find("Maximum Force") >=0:
                line_MaxForce = re.split(r"\s+", line) 
                if line_MaxForce[3] == '********':
                    line_MaxForce[3] = 10.000000
                    MaxForce_threshold.append(float(line_MaxForce[4]))
                else:
                    MaxForce.append(float(line_MaxForce[3]))
                    MaxForce_threshold.append(float(line_MaxForce[4]))
            if line.find("RMS     Force") >=0:
                line_RMSForce = re.split(r"\s+", line) 
                if line_RMSForce[3] == '********':
                    line_RMSForce[3] = 10.000000
                    RMSForce_threshold.append(float(line_RMSForce[4]))
                else:
                    RMSForce.append(float(line_RMSForce[3]))
                    RMSForce_threshold.append(float(line_RMSForce[4]))
            if line.find("Maximum Displacement") >=0:
                line_MaxDisp = re.split(r"\s+", line) 
                if line_MaxDisp[3] == '********':
                    line_MaxDisp[3] = 10.000000
                    MaxDisp_threshold.append(float(line_MaxDisp[4]))
                else:
                    MaxDisp.append(float(line_MaxDisp[3]))
                    MaxDisp_threshold.append(float(line_MaxDisp[4]))
            if line.find("RMS     Displacement") >=0:
                line_RMSDisp = re.split(r"\s+", line) 
                if line_RMSDisp ==  '********':
                    line_RMSDisp = 10.000000
                    RMSDisp_threshold.append(float(line_RMSDisp[4]))
                else:
                    RMSDisp.append(float(line_RMSDisp[3]))
                    RMSDisp_threshold.append(float(line_RMSDisp[4]))
    
        return (MaxForce, 
                RMSForce, 
                MaxDisp, 
                RMSDisp, 
                MaxForce_threshold[0], 
                RMSForce_threshold[0], 
                MaxDisp_threshold[0], 
                RMSDisp_threshold[0]
        )

if __name__ == '__main__':

    import sys
    import pandas as pd
    import matplotlib.pyplot as plt

    import ConvergenceJudge

    usage =f'Usage; {sys.argv} gaussian_log' 

    try:
        infilename = sys.argv[1]
        del sys.argv[1]
    except:
        print (usage); sys.exit()


    parselog = Monitor_log(infilename)
    Index_TState, Eprofile, OSprofile = parselog.Extract_ConversionProcess()
    Values_ForceDisp = parselog.Extract_ConversionCriterion()
    
    print(f'Number of step: {len(Eprofile)}')

    Eprofile = np.array(Eprofile)
    OSprofile = np.array(OSprofile)


    #For drowing#################
    fig = plt.figure()

    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)


    for i in range(len(Eprofile[0,:])):
        if i == Index_TState:
            ax1.plot((Eprofile[:,i]), marker='o')
        else:
            ax1.plot((Eprofile[:,i]), marker='x')


    for i in range(len(Values_ForceDisp[0])):
        for j in range(len(Values_ForceDisp)):
            ax2.plot((Values_ForceDisp[j]), marker='v')

    #For drowing################

    #For Judge convergence possibility
    
    Judge = ConvergenceJudge.ConvergenceJudge(Index_TState, Eprofile, Values_ForceDisp)
    print(Judge.judge())
    ################################


    plt.savefig("Energy_profile.pdf")
    plt.show()

