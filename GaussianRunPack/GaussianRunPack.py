import os, sys, math
import re, gc
import shutil
import subprocess
from numpy import *
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdmolops

import GaussianRunPack.AtomInfo
import GaussianRunPack.read_sdf
import GaussianRunPack.Exe_Gaussian
import GaussianRunPack.Get_ExcitedState
import GaussianRunPack.Get_MolCoordinate
import GaussianRunPack.Estimate_SpinContami
import GaussianRunPack.Get_ChargeSpin
import GaussianRunPack.Get_MOEnergy
import GaussianRunPack.chk2fchk

class GaussianDFTRun:

    def __init__(self, functional, basis, nproc, value, solvent, in_file, error):

        self.in_file = in_file
        self.functional = functional.lower()
        self.basis = basis.lower()
        self.nproc = nproc
        self.value = value.lower()
        self.solvent = solvent.lower()
        self.error = error

        self.mem = ''

    def Extract_SCFEnergy(self, lines):

        Energy=[]

        for line in lines:
            if line.find("SCF Done:  ") >=0:
                line_StateInfo = line.split()
        #print (line_StateInfo[4])
                Energy.append(float(line_StateInfo[4]))

        Comp_SS, Ideal_SS = GaussianRunPack.Estimate_SpinContami.Estimate_SpinDiff(lines)

        return Energy[-1]


    def Extract_values(self, infilename, option_array, Bondpair1, Bondpair2):

        opt       = int(option_array[0]) 
        nmr       = int(option_array[1]) 
        uv        = int(option_array[2])
        energy    = int(option_array[3])
        gap       = int(option_array[4])
        dipole    = int(option_array[5])
        deen      = int(option_array[6])
        stable2o2 = int(option_array[7])
        fluor     = int(option_array[8])
        tadf      = int(option_array[9])
        ipe       = int(option_array[10])
        eae       = int(option_array[11])
        pne       = int(option_array[12])
        nne       = int(option_array[13])
        cden      = int(option_array[14])

        with open(infilename, 'r') as ifile:
            lines = ifile.readlines()

        output = {}

        for line in lines:
            if line.find("Error termination") >=0:
                print("Gausssian is stopped due to some errors")
                return output

        print("Spliting links....")
        Links = self.SplitLinks(infilename)
        n = len(Links)
        print("The number of linkes = ", n)
  
        GS_lines = Links[1].splitlines()
        Links[1] = ""

        lines = ""

        if opt == 1:
            if Bondpair1 != []:
                print ("Optimization was performed...Check geometry...")
                MaxBondLength = GaussianRunPack.Get_MolCoordinate.MaxBondLength(GS_lines, Bondpair1, Bondpair2)
                output["GS_MaxBondLength"] = MaxBondLength
            else:
                output["GS_MaxBondLength"] = 0

        if gap == 1:
  
            NumAlphaElec, NumBetaElec, AlphaEigenVal, BetaEigenVal = GaussianRunPack.Get_MOEnergy.Extract_MO(GS_lines)

            if BetaEigenVal == []:
                Alpha_gap = 27.211*(AlphaEigenVal[NumAlphaElec]-AlphaEigenVal[NumAlphaElec-1])
                output["gap"] = Alpha_gap
            #    return Alpha_gap
            else:
                Alpha_gap = 27.211*(AlphaEigenVal[NumAlphaElec]-AlphaEigenVal[NumAlphaElec-1])
                Beta_gap = 27.211*(BetaEigenVal[NumBetaElec]-BetaEigenVal[NumBetaElec-1])
  
                output["gap"] = [Alpha_gap, Beta_gap]
            #    return Alpha_gap, Beta_gap
  
  
        if dipole == 1:
            Dipole_X = []
            Dipole_Y = []
            Dipole_Z = []
            Dipole_Total = []

            for line in GS_lines:
                if line.find(" X= ") >=0:
                    line_StateInfo = line.split()
                    #print (line_StateInfo[1])
                    Dipole_X.append(float(line_StateInfo[1]))
                    #print (line_StateInfo[3])
                    Dipole_Y.append(float(line_StateInfo[3]))
                    #print (line_StateInfo[5])
                    Dipole_Z.append(float(line_StateInfo[5]))
                    #print (line_StateInfo[7])
                    Dipole_Total.append(float(line_StateInfo[7]))
                
            output["dipole"] = [Dipole_X[-1], Dipole_Y[-1], Dipole_Z[-1], Dipole_Total[-1]]
        #    return Dipole_X[-1], Dipole_Y[-1], Dipole_Z[-1], Dipole_Total[-1]
  
        if energy == 1:

            output["Energy"] = self.Extract_SCFEnergy(GS_lines)


        if deen == 1:
            try:
                GS_Energy = output["Energy"]
            except KeyError:
                GS_Energy = self.Extract_SCFEnergy(GS_lines)
        
            Mol_atom, Mol_X, Mol_Y, Mol_Z = GaussianRunPack.Get_MolCoordinate.Extract_Coordinate(GS_lines)

#######Calculating Decomposed atoms total energy#######################
            decomposed_Energy = 0
#
            for i in range(len(Mol_atom)):    
                #print (Mol_atom[i], GaussianRunPack.AtomInfo.One_Atom_Energy(Mol_atom[i], self.functional, self.basis))
                decomposed_Energy += GaussianRunPack.AtomInfo.One_Atom_Energy(Mol_atom[i], self.functional, self.basis)
            print("Decomposed energy: ", decomposed_Energy)
#
#############################################################
  
        #    return deen 
            output["deen"] = GS_Energy - (decomposed_Energy)

        if stable2o2 == 1:
            NumAlphaElec, NumBetaElec, AlphaEigenVal, BetaEigenVal = GaussianRunPack.Get_MOEnergy.Extract_MO(GS_lines)
            O2_SOMO, O2_LUMO = GaussianRunPack.AtomInfo.O2_MO_refer(self.functional, self.basis)

            if BetaEigenVal == []:
                
               # OxidizedbyO2  >  0 oxidation by O2 is hard to occure.
               # OxidizedbyO2 <=  0 oxidation by O2 is easy to occure.

               # ReducedbyO2  >  0 reduction by O2 is hard to occure.
               # ReducedbyO2 <=  0 reduction by O2 is easy to occure.

                OxidizedbyO2 = O2_LUMO -AlphaEigenVal[NumAlphaElec-1]

                ReducedbyO2 =  AlphaEigenVal[NumAlphaElec] - O2_SOMO

                output["stable2o2"] = [OxidizedbyO2,ReducedbyO2]

            else:

                OxidizedbyO2 = O2_LUMO - BetaEigenVal[NumBetaElec-1] 

                ReducedbyO2 = AlphaEigenVal[NumAlphaElec] - O2_SOMO

                output["stable2o2"] = [OxidizedbyO2,ReducedbyO2]


        if cden == 1:
            output["cden"] = GaussianRunPack.Get_ChargeSpin.Extract_ChargeSpin(GS_lines)
  
        ##########################################################################
        # Index of Links
        ##########################################################################
        # Index = 0           : (always blank)
        # Index = 1           : Ground state             [if opt==1]
        # Index = 1+nmr       : NMR chemical shift of S0       [if nmr==1]
        # Index = 1+nmr+ipe   : Ionization potential      [if ipe==1]
        # Index = 1+nmr+ipe+epe   : Electronic affinity      [if eae==1]
        # Index = 1+nmr+ipe+epe+pne   : neutrization energy from cation  [if pne==1]
        # Index = 1+nmr+ipe+epe+pne+nne  : neutrization energy from anion  [if nne==1]
        # Index = 1+nmr+ipe+epe+pne+nne+1 : Virtical excitation (S0 -> S1) [uv]
        #########################################################################
        # Index = 2+nmr+ipe+epe+pne+nne+1 : Optimization of S1             [fluor or tadf] 
        # Index = 3+nmr+ipe+epe+pne+nne+1 : Optimization of T1             [tadf]
        ##########################################################################
  

        if nmr == 1:
            Element = []
            ppm = []
  
            Index = 1+nmr
            nmr_lines = Links[Index].splitlines()
            Links[Index]=""

            for line in nmr_lines:
                if line.find("Isotropic =  ") >=0:
                    line_Info = line.split()
                    #print (line_Info[1])
                    Element.append(line_Info[1])
                    #print(line_Info[4])
                    ppm.append(float(line_Info[4]))

        # calculating chemical shift for H, C, or Si

            for i in range(len(Element)):
                if Element[i] =="H" or  Element[i] =="C" or  Element[i] =="Si":
                    ppm[i] = GaussianRunPack.AtomInfo.One_TMS_refer(Element[i],  self.functional, self.basis)-ppm[i]
  
        #    return Element, ppm 
            output["nmr"] = [Element, ppm]
  
        #print (output)

        if ipe == 1 or eae == 1:

            try:
                GS_Energy = output["Energy"]
            except KeyError:
                GS_Energy = self.Extract_SCFEnergy(GS_lines)

            print (GS_Energy)

            if ipe == 1:

                Index = 1+nmr+ipe
                IP_lines = Links[Index].splitlines()
                Links[Index]=""

                IP_Energy =  self.Extract_SCFEnergy(IP_lines)
                IP_Comp_SS, IP_Ideal_SS = GaussianRunPack.Estimate_SpinContami.Estimate_SpinDiff(IP_lines)

                print (IP_Energy)

#                output["ipe"] = 27.211*(GS_Energy - IP_Energy)

######Normal ionization potential calculation####################
                output["ipe"] = [27.211*(IP_Energy - GS_Energy), IP_Comp_SS-IP_Ideal_SS]
#################################################################

            if eae == 1:

                Index = 1+nmr+ipe+eae
                EA_lines = Links[Index].splitlines()
                Links[Index]=""

                EA_Energy =  self.Extract_SCFEnergy(EA_lines)
                EA_Comp_SS, EA_Ideal_SS = GaussianRunPack.Estimate_SpinContami.Estimate_SpinDiff(EA_lines)

                print (EA_Energy)

#                output["eae"] = 27.211*(EA_Energy - GS_Energy)

######Normal electronic affinity calculation####################
                output["eae"] = [27.211*(GS_Energy - EA_Energy), EA_Comp_SS-EA_Ideal_SS]
#################################################################


        if pne >= 1 or nne >= 1:

            if pne == 1:

                Index = 1+nmr+ipe+eae+pne
                PC_lines = Links[Index].splitlines()
                Links[Index]=""

                VNP_lines = Links[Index+1].splitlines()
                Links[Index+1]=""
                pne += 1

                PC_Energy =  self.Extract_SCFEnergy(PC_lines)
                PC_Comp_SS, PC_Ideal_SS = GaussianRunPack.Estimate_SpinContami.Estimate_SpinDiff(PC_lines)
                VNP_Energy =  self.Extract_SCFEnergy(VNP_lines)
                VNP_Comp_SS, VNP_Ideal_SS = GaussianRunPack.Estimate_SpinContami.Estimate_SpinDiff(VNP_lines)

                print (PC_Energy)
                print (VNP_Energy)

                if Bondpair1 != []:
                    MaxBondLength = GaussianRunPack.Get_MolCoordinate.MaxBondLength(VNP_lines, Bondpair1, Bondpair2)
                else:
                    MaxBondLength = 0

                output["pne_MaxBondLength"] = MaxBondLength
                output["pne"] = [27.211*(PC_Energy - VNP_Energy), PC_Comp_SS-PC_Ideal_SS, VNP_Comp_SS-VNP_Ideal_SS]

            if nne == 1:

                Index = 1+nmr+ipe+eae+pne+nne 
                NC_lines = Links[Index].splitlines()
                Links[Index]=""

                VNN_lines = Links[Index+1].splitlines()
                Links[Index+1]=""
                nne +=1

                NC_Energy =  self.Extract_SCFEnergy(NC_lines)
                NC_Comp_SS, NC_Ideal_SS = GaussianRunPack.Estimate_SpinContami.Estimate_SpinDiff(NC_lines)
                VNN_Energy =  self.Extract_SCFEnergy(VNN_lines)
                VNN_Comp_SS, NC_Ideal_SS = GaussianRunPack.Estimate_SpinContami.Estimate_SpinDiff(VNN_lines)

                print (NC_Energy)
                print (VNN_Energy)

#                output["nne"] = 27.211*(NC_Energy - VNN_Energy)

                if Bondpair1 != []:
                    MaxBondLength = GaussianRunPack.Get_MolCoordinate.MaxBondLength(VNN_lines, Bondpair1, Bondpair2)
                else:
                    MaxBondLength = 0


                output["nne_MaxBondLength"] = MaxBondLength

######Normal electronic affinity calculation####################
                output["nne"] = [27.211*(VNN_Energy - NC_Energy), NC_Comp_SS-NC_Ideal_SS, VNN_Comp_SS-NC_Ideal_SS]
#################################################################

  
        if uv == 1:
            Index = 1+nmr+ipe+eae+pne+nne+1 
            lines = "" if Index >= n else Links[Index].splitlines()

            Found, Egrd, Eext, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
            CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden  \
            = GaussianRunPack.Get_ExcitedState.Extract_ExcitedState(lines)
#            print (Found)
#            print (Egrd)
#            print (Eext)

            output["uv"] = [WL_allowed, OS_allowed, CD_L_allowed, CD_OS_allowed]
            output["state_index"] = [State_allowed, State_forbidden]
            Links[Index]=""
  
        if fluor == 1 or tadf == 1:
            Index = 1+uv+nmr+ipe+eae+pne+nne+fluor 
            lines = "" if Index >= n else Links[Index].splitlines()

            S_Found, S_Egrd, S_Eext, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
            CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden  \
            = GaussianRunPack.Get_ExcitedState.Extract_ExcitedState(lines)

            if Bondpair1 != []:
                MaxBondLength = GaussianRunPack.Get_MolCoordinate.MaxBondLength(lines, Bondpair1, Bondpair2)
            else:
                MaxBondLength = 0

            output["MinEtarget"] = S_Eext
            output["Min_MaxBondLength"] = MaxBondLength
            output["fluor"] = [WL_allowed, OS_allowed, CD_L_allowed, CD_OS_allowed]
            Links[Index]=""
  
        if tadf == 1:
            Index = 1+uv+nmr+ipe+eae+pne+nne+fluor+tadf
            lines = "" if Index >= n else Links[Index].splitlines()

            T_Found, T_Egrd, T_Eext, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
            CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden  \
            = GaussianRunPack.Get_ExcitedState.Extract_ExcitedState(lines)

            if Bondpair1 != []:
                MaxBondLength = GaussianRunPack.Get_MolCoordinate.MaxBondLength(lines, Bondpair1, Bondpair2)
            else:
                MaxBondLength = 0

            output["T_Min"] = T_Eext
            output["T_Min_MaxBondLength"] = MaxBondLength
            output["T_Phos"] = [WL_forbidden, OS_forbidden, CD_L_forbidden, CD_OS_forbidden]
            TADF_Eng = 0.0
            if S_Found and T_Found:
               TADF_Eng = S_Eext - T_Eext
            output["Delta(S-T)"] = TADF_Eng
            Links[Index]=""
        
        del lines
        del Links
  
        return output

    def MakeSolventLine(self):

        s = ''
        s_solvent = ''

        try:
            float(self.solvent)
        except ValueError:
            print ('Solvent effect is included by PCM')
            s = s +'SCRF=(PCM, solvent=' + self.solvent + ')'
            s = s +'\n'

        else:
            if self.solvent == '0':
                return s, s_solvent

            else:
                print ('Solvent effect is included by PCM')
                s = s +'SCRF=(PCM, solvent=Generic, Read)'
                s = s +'\n'
                s_solvent = s_solvent+'EPS='+self.solvent
                s_solvent = s_solvent+'\n'
                s_solvent = s_solvent+'Radii=UA0'
                s_solvent = s_solvent+'\n'
                s_solvent = s_solvent+'\n'

        return s, s_solvent

    def MakeLinkTD(self, line_chk, line_method, State, SCRF, SCRF_read, Opt, targetstate=1):

        self.MakeSolventLine()

        Jobname_line =line_chk.split('=')
        Jobname = Jobname_line[-1]

### Setting NState######
        if Opt == True:
            NState = targetstate + 4
        else:
            NState = 20
#######################

        line_oldchk = '%Oldchk='+Jobname
        if State == None:
            line_newchk = '%chk='+ Jobname + '_ExOptState' + str(targetstate)
        else:
            line_newchk = '%chk='+ Jobname + '_ExOptState' + State

        if State == None:
            line_method_TD = 'TD(Nstate='+ str(NState) + ', root='+ str(targetstate) +')'
        else:
            line_method_TD = 'TD(Nstate='+ str(NState) +', ' + State + ', root='+ str(targetstate) +')'
        line_readMOGeom = 'Geom=AllCheck Guess=Read'
        s = ''
        s = s + '--Link1--\n'
        if self.mem != '':
            line_mem = '%mem='+str(self.mem)
            s = s + line_mem      
            s = s + '\n'
        if self.nproc > 1 :
            line_proc = '%nproc='+str(self.nproc)
            s = s + line_proc
            s = s + '\n'
        if Opt == True:
            s = s + line_oldchk + '\n'
            s = s + line_newchk + '\n'
        else:
            s = s + line_chk +'\n'
        s = s + line_method
        s = s + '\n'
        s = s + line_method_TD
        s = s + '\n'
        s = s + line_readMOGeom
        s = s + '\n'
        s = s + SCRF
        if Opt == True:
#            s = s + "Opt=(MaxCycles=50)"
            s = s + "Opt"
            s = s + '\n'
        s = s + '\n' 

        s = s + SCRF_read

        return s


    def run_gaussian(self):
        
        infilename = self.in_file

        option_line = self.value    

        options = option_line.split()

        option_array = zeros(15)
        option_array_Ex = zeros(15)

        targetstate = 1

        PreGauInput = infilename.split('.')
    
        GauInputName = PreGauInput[0]+'.com'    

####File type of input?######################
        ReadFromchk = 0 
        ReadFromsdf = 0 

        if PreGauInput[1] == "sdf":
            ReadFromsdf = 1 
            Mol_atom, X, Y, Z, TotalCharge, SpinMulti, Bondpair1, Bondpair2 = GaussianRunPack.read_sdf.read_sdf(infilename)
        elif PreGauInput[1] == "chk":
            ReadFromchk = 1 
            Bondpair1 = []
            Bondpair2 = []
        else:
            print ("Invalid input file")
############################################


        for i in range(len(options)):
            option = options[i]
            if option == 'opt':
#                opt = 1
                option_array[0] = 1
            elif option == 'nmr':
#                nmr = 1
                option_array[1] = 1
                #print ('nmr')
            elif option == 'uv':
#                uv = 1
                option_array[2] = 1
                #print ('uv')
            elif option == 'energy':
#                energy = 1
                option_array[3] = 1
                #print ('energy')
            elif option == 'homolumo':
#                gap = 1
                option_array[4] = 1
                #print ('HOMO/LUMO')
            elif option == 'dipole':
#                dipole = 1
                option_array[5] = 1
                #print ('dipole')
            elif option == 'deen':
#                deen = 1
                option_array[6] = 1
                #print ('Decomposition energy')
            elif option == 'stable2o2':
#                stable2o2 = 1
                option_array[7] = 1
                #print ('Stability to O2')
#            elif option == 'fluor':
            elif 'fluor' in option:
#                fluor = 1
                option_array[2] = 1
                if SpinMulti == 1:
                    option_array[8] = 1
                else:
                    option_array_Ex[8] = 1
                if '=' in option:
                    in_target = option.split("=")
                    targetstate = int(in_target[-1])
                #print ('Fluorescence')
            elif option == 'tadf':
#                tadf = 1
                option_array[2] = 1
                if SpinMulti == 1:
                    option_array[8] = 1
                    option_array[9] = 1
                else:
                    option_array_Ex[8] = 1
                    option_array_Ex[9] = 1
                #print ('Thermally Activated Delayed Fluorescence')
            elif option == 'ipe':
#                ipe = 1
                option_array[10] = 1
                #print ('Ionization potential')
            elif option == 'eae':
#                eae = 1
                option_array[11] = 1
                #print ('Electronic affinity')
            elif option == 'pne':
#                pne = 1
                option_array[12] = 1
                #print ('Neutraization energy from cation')
            elif option == 'nne':
#                nne = 1
                option_array[13] = 1
                #print ('Neutraization energy from anion')
            elif option == 'cden':
#                nne = 1
                option_array[14] = 1
                #print ('Neutraization energy from anion')
            else:
                print('invalid option: ', option)


        line_chk = '%chk='+PreGauInput[0]
        line_oldchk = '%Oldchk='+PreGauInput[0]

        line_method = '#'+self.functional+'/'+self.basis+' test'
        line_o_method = '#u'+self.functional+'/'+self.basis+' test'
        line_c_method = '#r'+self.functional+'/'+self.basis+' test'

        line_comment = infilename

        ofile = open(GauInputName ,'w')

#####For reading geomerty and MO from checkpoint file###
        line_readMOGeom = 'Geom=AllCheck Guess=Read'
        line_readOnlyMOGeom = 'Geom=CheckPoint Guess=Read'
        line_readGeom = 'Geom=Checkpoint'
########################################################

        if self.mem != '':
            line_mem = '%mem='+str(self.mem)
            ofile.write(line_mem+'\n')

        if self.nproc > 1 :
            line_proc = '%nproc='+str(self.nproc)
            ofile.write(line_proc+'\n')
    
        ofile.write(line_chk+'\n')
        ofile.write(line_o_method+'\n')
        if option_array[0] == 1: # opt == 1
            ofile.write('Opt=(MaxCycles=100)\n')


#######Solvent effect############################
        SCRF, SCRF_read = self.MakeSolventLine()
        ofile.write(SCRF)

#####Reading Geometry and MO from Checkpoint file
        if ReadFromchk == 1:
            ofile.write(line_readMOGeom+'\n')
#################################################

        ofile.write('\n')

#####Reading Geometry from sdf file#####

        if ReadFromsdf == 1:

            ofile.write(line_comment+'\n')
            ofile.write('\n')

            ofile.write('%5d %5d \n' % (TotalCharge, SpinMulti))

            for j in range(len(Mol_atom)):
                ofile.write('%-4s % 10.5f  % 10.5f  % 10.5f \n'
                 % (Mol_atom[j],X[j], Y[j], Z[j]))

            ofile.write('\n')
            ofile.write(SCRF_read)
#######################################

        if option_array[1] == 1: # nmr == 1

            ofile.write('--Link1--\n')

            if self.mem != '':
                line_mem = '%mem='+str(self.mem)
                ofile.write(line_mem+'\n')


            if self.nproc > 1 :
                line_proc = '%nproc='+str(self.nproc)
                ofile.write(line_proc+'\n')

            ofile.write(line_chk+'\n')
            ofile.write(line_method+'\n')

            line_method_nmr = 'NMR'
            ofile.write(line_method_nmr+'\n')
            ofile.write(SCRF)

            ofile.write(line_readMOGeom+'\n')
            ofile.write('\n')
            ofile.write(SCRF_read)

        if option_array[10] == 1: # ipe == 1

            ofile.write('--Link1--\n')

            if self.mem != '':
                line_mem = '%mem='+str(self.mem)
                ofile.write(line_mem+'\n')


            if self.nproc > 1 :
                line_proc = '%nproc='+str(self.nproc)
                ofile.write(line_proc+'\n')

            ofile.write(line_chk+'_IP\n')
            ofile.write(line_oldchk+'\n')
            ofile.write(line_o_method+'\n')
            ofile.write(SCRF)

            ofile.write(line_readOnlyMOGeom+'\n')
            ofile.write('\n')

            ofile.write('ionization potential calculation\n')
            ofile.write('\n')

            ofile.write('%5d %5d \n' % (TotalCharge+1, SpinMulti+1))
            ofile.write('\n')
            ofile.write(SCRF_read)

        if option_array[11] == 1: # eae == 1

            ofile.write('--Link1--\n')

            if self.mem != '':
                line_mem = '%mem='+str(self.mem)
                ofile.write(line_mem+'\n')


            if self.nproc > 1 :
                line_proc = '%nproc='+str(self.nproc)
                ofile.write(line_proc+'\n')

            ofile.write(line_oldchk+'\n')
            ofile.write(line_chk+'_EA\n')
            ofile.write(line_o_method+'\n')
            ofile.write(SCRF)

            ofile.write(line_readOnlyMOGeom+'\n')
            ofile.write('\n')

            ofile.write('electronic affinity calculation\n')
            ofile.write('\n')

            ofile.write('%5d %5d \n' % (TotalCharge-1, SpinMulti+1))
            ofile.write('\n')
            ofile.write(SCRF_read)

        if option_array[12] == 1: # pne == 1

            ofile.write('--Link1--\n')

            if self.mem != '':
                line_mem = '%mem='+str(self.mem)
                ofile.write(line_mem+'\n')


            if self.nproc > 1 :
                line_proc = '%nproc='+str(self.nproc)
                ofile.write(line_proc+'\n')

            ofile.write(line_chk+'_PC\n')
            ofile.write(line_oldchk+'\n')
            ofile.write(line_o_method+' Opt'+'\n')
            ofile.write(SCRF)

            ofile.write(line_readOnlyMOGeom+'\n')
            ofile.write('\n')

            ofile.write('neutrization energy calculation from a cation\n')
            ofile.write('\n')

            ofile.write('%5d %5d \n' % (TotalCharge+1, SpinMulti+1))
            ofile.write('\n')
            ofile.write(SCRF_read)

            ofile.write('--Link1--\n')

            if self.mem != '':
                line_mem = '%mem='+str(self.mem)
                ofile.write(line_mem+'\n')

            if self.nproc > 1 :
                line_proc = '%nproc='+str(self.nproc)
                ofile.write(line_proc+'\n')

            ofile.write(line_chk+'_PC\n')
            ofile.write(line_oldchk+'\n')
            ofile.write(line_o_method+'\n')
            ofile.write(SCRF)

            ofile.write(line_readGeom+'\n')
            ofile.write('\n')

            ofile.write('neutrization energy calculation from a cation\n')
            ofile.write('\n')

            ofile.write('%5d %5d \n' % (TotalCharge, SpinMulti))
            ofile.write('\n')
            ofile.write(SCRF_read)

        if option_array[13] == 1: #nne == 1

            ofile.write('--Link1--\n')

            if self.mem != '':
                line_mem = '%mem='+str(self.mem)
                ofile.write(line_mem+'\n')


            if self.nproc > 1 :
                line_proc = '%nproc='+str(self.nproc)
                ofile.write(line_proc+'\n')

            ofile.write(line_oldchk+'\n')
            ofile.write(line_chk+'_NC\n')
            ofile.write(line_o_method+' Opt'+'\n')
            ofile.write(SCRF)

            ofile.write(line_readOnlyMOGeom+'\n')
            ofile.write('\n')

            ofile.write('neutrization energy calculation from an anion\n')
            ofile.write('\n')

            ofile.write('%5d %5d \n' % (TotalCharge-1, SpinMulti+1))
            ofile.write('\n')
            ofile.write(SCRF_read)

            ofile.write('--Link1--\n')

            if self.mem != '':
                line_mem = '%mem='+str(self.mem)
                ofile.write(line_mem+'\n')

            if self.nproc > 1 :
                line_proc = '%nproc='+str(self.nproc)
                ofile.write(line_proc+'\n')

            ofile.write(line_oldchk+'\n')
            ofile.write(line_chk+'_NC\n')
            ofile.write(line_o_method+'\n')
            ofile.write(SCRF)

            ofile.write(line_readGeom+'\n')
            ofile.write('\n')

            ofile.write('neutrization energy calculation from an anion\n')
            ofile.write('\n')

            ofile.write('%5d %5d \n' % (TotalCharge, SpinMulti))
            ofile.write('\n')
            ofile.write(SCRF_read)


        if option_array[2] == 1 : #uv == 1 or fluor==1 or tadf == 1
            sTD = self.MakeLinkTD(line_chk, line_method, None, SCRF, SCRF_read, False, targetstate )
            ofile.write(sTD) 
        if option_array[8] == 1 :
            sTD = self.MakeLinkTD(line_chk, line_c_method, 'Singlet',  SCRF, SCRF_read, True, targetstate )
            ofile.write(sTD) 
        if option_array[9] == 1 :
            sTD = self.MakeLinkTD(line_chk, line_c_method, 'Triplet',  SCRF, SCRF_read, True, targetstate )
            ofile.write(sTD) 

        ofile.write('\n') 
        ofile.close()


############Run Gaussian##############################

        if os.path.isdir(PreGauInput[0]):
            shutil.rmtree(PreGauInput[0])
        os.mkdir(PreGauInput[0])
        shutil.move(GauInputName, PreGauInput[0])

        if ReadFromchk == 1 :
            shutil.move(infilename, PreGauInput[0]) 

        os.chdir(PreGauInput[0])
        
        GaussianRunPack.Exe_Gaussian.exe_Gaussian(PreGauInput[0])

        logfile = PreGauInput[0]+'.log'

        output_dic = self.Extract_values(logfile, option_array,Bondpair1, Bondpair2)

        if option_array_Ex[8] == 1 or option_array_Ex[9] == 1: #fluor == 1 or tadf == 1 for open shel
    
            compute_state = output_dic["state_index"][0][int(targetstate)-1] 
            JobName_ExOpt = PreGauInput[0] + "_ExOpt"
            GauInputName_ExOpt = JobName_ExOpt + ".com"
            logfile_ExOpt = JobName_ExOpt + ".log"
            ofile_ExOpt = open(GauInputName_ExOpt ,'w')
    
            if self.mem != '':
                line_mem = '%mem='+str(self.mem)
                ofile_ExOpt.write(line_mem+'\n')
    
            if self.nproc > 1 :
                line_proc = '%nproc='+str(self.nproc)
                ofile_ExOpt.write(line_proc+'\n')
       
                ofile_ExOpt.write(line_chk+'\n')
                ofile_ExOpt.write(line_o_method+'\n')
    
            ofile_ExOpt.write(line_readMOGeom+'\n')
            ofile_ExOpt.write('\n')
    
            sTD = self.MakeLinkTD(line_chk, line_o_method, None, SCRF, SCRF_read, True, compute_state)
            ofile_ExOpt.write(sTD)
    
            if option_array_Ex[9] == 1: #tadf == 1
                compute_state = output_dic["state_index"][1][int(targetstate)-1] 
                sTD = self.MakeLinkTD(line_chk, line_o_method, None, SCRF, SCRF_read, True, compute_state)
                ofile_ExOpt.write(sTD) 
    
            ofile_ExOpt.close()
    
            GaussianRunPack.Exe_Gaussian.exe_Gaussian(JobName_ExOpt)
    
            output_dic_Ex = self.Extract_values(logfile_ExOpt, option_array_Ex, Bondpair1, Bondpair2)
    
            output_dic.update(output_dic_Ex)

        GaussianRunPack.chk2fchk.Get_chklist()

        os.chdir("..")

        return(output_dic)

    def SplitLinks(self, logfile):
        with open(logfile, 'r') as f:
            lines = f.readlines()

        Links = []
        Link = ""

        for line in lines:
            if line.find("Initial command:")>0:
                Links.append(Link)
                Link = ""

            Link = Link + line

        Links.append(Link)
        return Links


        
