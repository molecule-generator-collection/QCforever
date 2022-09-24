import os, sys, math
import re, gc
import shutil
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdmolops

import GaussianRunPack.AtomInfo
import GaussianRunPack.read_sdf
import GaussianRunPack.read_xyz
import GaussianRunPack.Exe_Gaussian
import GaussianRunPack.Get_ExcitedState
import GaussianRunPack.Get_MolCoordinate
import GaussianRunPack.Get_MolCoordinate_fchk
import GaussianRunPack.Estimate_SpinContami
import GaussianRunPack.Get_ChargeSpin
import GaussianRunPack.Get_MOEnergy
import GaussianRunPack.Get_MOEnergy_fchk
import GaussianRunPack.chk2fchk
import GaussianRunPack.fchk2chk
import GaussianRunPack.UV_similarity
import GaussianRunPack.Get_FreqPro

Eh2kJmol = 2625.5

class GaussianDFTRun:

    def __init__(self, functional, basis, nproc, value, in_file, solvent="0", error=0,  restart=True):

        self.in_file = in_file
        self.functional = functional.lower()
        self.basis = basis.lower()
        self.nproc = nproc
#        self.value = value.lower()
        self.value = value
        self.solvent = solvent.lower()
        self.restart = restart
        self.error = error

        self.mem = ''
        self.timexe = 60*60*80
        self.SpecTotalCharge = np.nan
        self.SpecSpinMulti = np.nan

        self.ref_uv_path = ''

    def Extract_SCFEnergy(self, lines):

        Energy=[]

        for line in lines:
            if line.find("SCF Done:  ") >=0:
                line_StateInfo = line.split()
        #print (line_StateInfo[4])
                Energy.append(float(line_StateInfo[4]))

        Comp_SS, Ideal_SS = GaussianRunPack.Estimate_SpinContami.Estimate_SpinDiff(lines)

        Energy_Spin = [Energy[-1], Comp_SS-Ideal_SS]

        return Energy_Spin


    def Extract_values(self, jobname, option_array, Bondpair1, Bondpair2):

        opt       = int(option_array[0]) 
        freq      = int(option_array[1]) 
        nmr       = int(option_array[2]) 
        uv        = int(option_array[3])
        energy    = int(option_array[4])
        homolumo  = int(option_array[5])
        dipole    = int(option_array[6])
        deen      = int(option_array[7])
        stable2o2 = int(option_array[8])
        fluor     = int(option_array[9])
        tadf      = int(option_array[10])
        vip       = int(option_array[11])
        vea       = int(option_array[12])
        aip       = int(option_array[13])
        aea       = int(option_array[14])
        cden      = int(option_array[15])
        pka       = int(option_array[16])
        symm      = int(option_array[17])

        infilename = jobname+'.log'
        fchkname = jobname+'.fchk'

        with open(infilename, 'r') as ifile:
            lines = ifile.readlines()

        output = {}

#        for line in lines:
#            if line.find("Error termination") >=0:
#                print("Gausssian is stopped due to some errors")
#                return output

        print("Spliting links....")
        Links = self.SplitLinks(infilename)
        n = len(Links)
        print("The number of linkes = ", n)

#####Clean lines
        lines = ""

#####For extracting symmetry of molecule without SCF
        if symm ==1:
            Symm_lines = Links[symm].splitlines()
            pGroup = "C1"
            for line in Symm_lines:
                if line.find("Full point group  ") >=0:
                    line_symmInfo = line.split()
                    pGroup = line_symmInfo[3]

            output["symm"] = pGroup
####################################################


#####For extracting properties of molecule after SCF        
        if 1 in option_array[0:17]: 
            GS_lines = Links[1+symm].splitlines()
            Links[1+symm] = ""

        if opt == 1:
            if Bondpair1 != []:
                print ("Optimization was performed...Check geometry...")
                MaxBondLength = GaussianRunPack.Get_MolCoordinate.MaxBondLength(GS_lines, Bondpair1, Bondpair2)
                output["GS_MaxBondLength"] = MaxBondLength
            else:
                output["GS_MaxBondLength"] = 0

        if freq == 1:
            print("For getting frequency")
            Freq, IR, Raman, E_zp,  E_t, E_enth, E_free, Ei, Cv, St = GaussianRunPack.Get_FreqPro.Extract_Freq(GS_lines) 
            output["freq"] = Freq 
            output["IR"] = IR
            output["Raman"] = Raman
            output["Ezp"] = E_zp
            output["Et"] = E_t
            output["E_enth"] = E_enth
            output["E_free"] = E_free
            output["Ei"] = Ei
            output["Cv"] = Cv
            output["Si"] = St

        if homolumo == 1:

            with open(fchkname, 'r') as ifile:
                fchk_lines = ifile.readlines()
  
            NumAlphaElec, NumBetaElec, AlphaEigenVal, BetaEigenVal = GaussianRunPack.Get_MOEnergy_fchk.Extract_MO(fchk_lines)

          #  NumAlphaElec, NumBetaElec, AlphaEigenVal, BetaEigenVal = GaussianRunPack.Get_MOEnergy.Extract_MO(GS_lines)
            if BetaEigenVal == []:
                Alpha_gap = 27.211*(AlphaEigenVal[NumAlphaElec]-AlphaEigenVal[NumAlphaElec-1])
                output["homolumo"] = Alpha_gap
            #    return Alpha_gap
            else:
                Alpha_gap = 27.211*(AlphaEigenVal[NumAlphaElec]-AlphaEigenVal[NumAlphaElec-1])
                Beta_gap = 27.211*(BetaEigenVal[NumBetaElec]-BetaEigenVal[NumBetaElec-1])
  
                output["homolumo"] = [Alpha_gap, Beta_gap]
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
                GS_Energy = output["Energy"][0]
            except KeyError:
                output["Energy"] = self.Extract_SCFEnergy(GS_lines)
                GS_Energy = output["Energy"][0] 
        
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

            try:
                NumAlphaElec
        
            except:
                with open(fchkname, 'r') as fchkfile:
                    fchk_lines = fchkfile.readlines()
  
                NumAlphaElec, NumBetaElec, AlphaEigenVal, BetaEigenVal = GaussianRunPack.Get_MOEnergy_fchk.Extract_MO(fchk_lines)
             #   NumAlphaElec, NumBetaElec, AlphaEigenVal, BetaEigenVal = GaussianRunPack.Get_MOEnergy.Extract_MO(GS_lines)
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
        # Index = 1+symm           : Ground state             [if opt==1]
        # Index = 1+symm+nmr       : NMR chemical shift of S0       [if nmr==1]
        # Index = 1+symm+nmr+vip   : Ionization potential      [if vip==1]
        # Index = 1+symm+nmr+vip+vea   : Electronic affinity      [if vea==1]
        # Index = 1+symm+nmr+vip+vea+aip   : adiabatic ionization potential  [if aip==1]
        # Index = 1+symm+nmr+vip+vea+aip+aea  : adiabatic electronic affinity  [if aea==1]
        # Index = 1+symm+nmr+vip+vea+aip+aea+1 : Virtical excitation (S0 -> S1) [uv]
        #########################################################################
        # Index = 2+symm+nmr+vip+vea+aip+aea+1 : Optimization of S1             [fluor or tadf] 
        # Index = 3+symm+nmr+vip+vea+aip+aea+1 : Optimization of T1             [tadf]
        ##########################################################################
  

        if nmr == 1:
            Element = []
            ppm = []
  
            Index = 1+symm+nmr
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

        if vip == 1 or vea == 1:

            try:
                GS_Energy = output["Energy"][0]
            except KeyError:
                output["Energy"] = self.Extract_SCFEnergy(GS_lines)
                GS_Energy = output["Energy"][0] 

            print (GS_Energy)

            if vip == 1:

                Index = 1+symm+nmr+vip
                IP_lines = Links[Index].splitlines()
                Links[Index]=""

                IP_Energy_SS =  self.Extract_SCFEnergy(IP_lines)
#                IP_Comp_SS, IP_Ideal_SS = GaussianRunPack.Estimate_SpinContami.Estimate_SpinDiff(IP_lines)

######Normal ionization potential calculation####################
                output["vip"] = [27.211*(IP_Energy_SS[0] - GS_Energy), IP_Energy_SS[1]]
#################################################################

            if vea == 1:

                Index = 1+symm+nmr+vip+vea
                EA_lines = Links[Index].splitlines()
                Links[Index]=""

                EA_Energy_SS =  self.Extract_SCFEnergy(EA_lines)
#                EA_Comp_SS, EA_Ideal_SS = GaussianRunPack.Estimate_SpinContami.Estimate_SpinDiff(EA_lines)

######Normal electronic affinity calculation####################
                output["vea"] = [27.211*(GS_Energy - EA_Energy_SS[0]), EA_Energy_SS[1]]
#################################################################


        if aip >= 1 or aea >= 1:

            try:
                GS_Energy = output["Energy"][0]
            except KeyError:
                output["Energy"] = self.Extract_SCFEnergy(GS_lines)
                GS_Energy = output["Energy"][0] 

            if aip == 1:

                Index = 1+symm+nmr+vip+vea+aip
                PC_lines = Links[Index].splitlines()
                Links[Index]=""

                VNP_lines = Links[Index+1].splitlines()
                Links[Index+1]=""
                aip += 1

                PC_Energy_SS =  self.Extract_SCFEnergy(PC_lines)
#                PC_Comp_SS, PC_Ideal_SS = GaussianRunPack.Estimate_SpinContami.Estimate_SpinDiff(PC_lines)
                VNP_Energy_SS =  self.Extract_SCFEnergy(VNP_lines)
#                VNP_Comp_SS, VNP_Ideal_SS = GaussianRunPack.Estimate_SpinContami.Estimate_SpinDiff(VNP_lines)

                if Bondpair1 != []:
                    MaxBondLength = GaussianRunPack.Get_MolCoordinate.MaxBondLength(VNP_lines, Bondpair1, Bondpair2)
                else:
                    MaxBondLength = 0

                output["relaxedIP_MaxBondLength"] = MaxBondLength
                output["aip"] = [27.211*(PC_Energy_SS[0] - GS_Energy), 27.211*(PC_Energy_SS[0] - VNP_Energy_SS[0]), PC_Energy_SS[1], VNP_Energy_SS[1]]

            if aea == 1:

                Index = 1+symm+nmr+vip+vea+aip+aea 
                NC_lines = Links[Index].splitlines()
                Links[Index]=""

                VNN_lines = Links[Index+1].splitlines()
                Links[Index+1]=""
                aea +=1

                NC_Energy_SS =  self.Extract_SCFEnergy(NC_lines)
#                NC_Comp_SS, NC_Ideal_SS = GaussianRunPack.Estimate_SpinContami.Estimate_SpinDiff(NC_lines)
                VNN_Energy_SS =  self.Extract_SCFEnergy(VNN_lines)
#                VNN_Comp_SS, NC_Ideal_SS = GaussianRunPack.Estimate_SpinContami.Estimate_SpinDiff(VNN_lines)

                if Bondpair1 != []:
                    MaxBondLength = GaussianRunPack.Get_MolCoordinate.MaxBondLength(VNN_lines, Bondpair1, Bondpair2)
                else:
                    MaxBondLength = 0


                output["relaxedEA_MaxBondLength"] = MaxBondLength

######Normal electronic affinity calculation####################
                output["aea"] = [27.211*(GS_Energy - NC_Energy_SS[0]), 27.211*(VNN_Energy_SS[0] - NC_Energy_SS[0]), NC_Energy_SS[1], VNN_Energy_SS[1]]
#################################################################

  
        if uv == 1:
            Index = 1+symm+nmr+vip+vea+aip+aea+1 
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

            if self.ref_uv_path != '':
                ref_uv = {}

                print ("Read the reference spectrum...", self.ref_uv_path)

                ref_uv["uv"] = GaussianRunPack.UV_similarity.read_data(self.ref_uv_path)
    
                S, D = GaussianRunPack.UV_similarity.smililarity_dissimilarity(ref_uv["uv"][0], ref_uv["uv"][1], \
                output["uv"][0], output["uv"][1])

                output["Similality/Disdimilarity"] =[S, D]

  
        if fluor == 1 or tadf == 1:
            Index = 1+symm+uv+nmr+vip+vea+aip+aea+fluor 
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
            Index = 1+symm+uv+nmr+vip+vea+aip+aea+fluor+tadf
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

######Setting NState######
        if Opt == True:
            NState = targetstate + 4
        else:
            NState = 20
##########################


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

        option_array = np.zeros(18)
        option_array_Ex = np.zeros(18)
        option_array_pka = np.zeros(18)

        targetstate = 1

        PreGauInput = infilename.split('.')
    
        GauInputName = PreGauInput[0]+'.com'    

####File type of input?######################
        ReadFromchk = 0 
        ReadFromsdf = 0 
        ReadFromxyz = 0 

        if PreGauInput[1] == "sdf":
            ReadFromsdf = 1 
            Mol_atom, X, Y, Z, TotalCharge, SpinMulti, Bondpair1, Bondpair2 = GaussianRunPack.read_sdf.read_sdf(infilename)
        elif PreGauInput[1] == "xyz":
            ReadFromxyz = 1
            Mol_atom, X, Y, Z, TotalCharge, SpinMulti = GaussianRunPack.read_xyz.read_xyz(infilename)
            Bondpair1 = []
            Bondpair2 = []
        elif PreGauInput[1] == "chk":
            ReadFromchk = 1 
            Bondpair1 = []
            Bondpair2 = []
        elif PreGauInput[1] == "fchk":
            TotalCharge, SpinMulti = GaussianRunPack.fchk2chk.Get_fchk(PreGauInput[0])
            ReadFromchk = 1 
            Bondpair1 = []
            Bondpair2 = []
        else:
            print ("Invalid input file")

##############################################


#######Setting electronic structure###########

        if np.isnan(self.SpecTotalCharge) !=  True:
            TotalCharge = self.SpecTotalCharge

        if np.isnan(self.SpecSpinMulti) != True:
            SpinMulti = self.SpecSpinMulti

##############################################

        for i in range(len(options)):
            option = options[i]
            if option.lower() == 'opt':
#                opt = 1
                option_array[0] = 1
            elif option.lower() == 'freq':
#                nmr = 1
                option_array[1] = 1
                #print ('freq')
            elif option.lower() == 'nmr':
#                nmr = 1
                option_array[2] = 1
                #print ('nmr')
            elif 'uv' in option.lower():
#                uv = 1
                option_array[3] = 1
                if '=' in option:
                    ref_spectrum = option.split("=")
                    self.ref_uv_path = ref_spectrum[-1]
                #print ('uv')
            elif option.lower() == 'energy':
#                energy = 1
                option_array[4] = 1
                #print ('energy')
            elif option.lower() == 'homolumo':
#                homolumo = 1
                option_array[5] = 1
                #print ('HOMO/LUMO')
            elif option.lower() == 'dipole':
#                dipole = 1
                option_array[6] = 1
                #print ('dipole')
            elif option.lower() == 'deen':
#                deen = 1
                option_array[7] = 1
                #print ('Decomposition energy')
            elif option.lower() == 'stable2o2':
#                stable2o2 = 1
                option_array[8] = 1
                #print ('Stability to O2')
#            elif option == 'fluor':
            elif 'fluor' in option.lower():
#                fluor = 1
                option_array[3] = 1
                if SpinMulti == 1:
                    option_array[9] = 1
                else:
                    option_array_Ex[9] = 1
                if '=' in option:
                    in_target = option.split("=")
                    targetstate = int(in_target[-1])
                #print ('Fluorescence')
            elif option.lower() == 'tadf':
#                tadf = 1
                option_array[2] = 1
                if SpinMulti == 1:
                    option_array[9] = 1
                    option_array[10] = 1
                else:
                    option_array_Ex[9] = 1
                    option_array_Ex[10] = 1
                #print ('Thermally Activated Delayed Fluorescence')
            elif option.lower() == 'vip':
#                vip = 1
                option_array[11] = 1
                #print ('Ionization potential')
            elif option.lower() == 'vea':
#                vea = 1
                option_array[12] = 1
                #print ('Electronic affinity')
            elif option.lower() == 'aip':
#                aip = 1
                option_array[11] = 1
                option_array[13] = 1
                #print ('Neutraization energy from cation')
            elif option.lower() == 'aea':
#                aea = 1
                option_array[12] = 1
                option_array[14] = 1
                #print ('Neutraization energy from anion')
            elif option.lower() == 'cden':
#                nne = 1
                option_array[15] = 1
                #print ('Neutraization energy from anion')
            elif option.lower() == 'pka':
                option_array[4] = 1
                option_array[15] = 1
                option_array[16] = 1
                option_array_pka[4] = 1
            elif option.lower() == 'symm':
                option_array[17] = 1
                #print ('Calculate symmetry')
            else:
                print('invalid option: ', option)

        line_system = ''
        if self.mem !='' or self.nproc > 1:

            if self.mem != '':
                line_mem = '%mem='+str(self.mem)+'\n'
                line_system += line_mem

            if self.nproc > 1 :
                line_proc = '%nproc='+str(self.nproc)+'\n'
                line_system += line_proc

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

        if option_array[17] == 1:
            ofile.write(line_system)
            ofile.write(line_chk+'\n')
            ofile.write(line_o_method+'\n')
            ofile.write('Guess=Only Symmetry=loose'+'\n')
           # ofile.write('Guess=Only Symmetry=on'+'\n')
#####Reading Geometry and MO from Checkpoint file
            if ReadFromchk == 1:
                ofile.write(line_readMOGeom+'\n')
#################################################

            ofile.write('\n')
            #ofile.write(line_comment+'\n')
            #ofile.write('\n')

            #ofile.write('%5d %5d \n' % (TotalCharge, SpinMulti))

            if ReadFromsdf == 1 or ReadFromxyz == 1:

                ofile.write(line_comment+'\n')
                ofile.write('\n')

                ofile.write('%5d %5d \n' % (TotalCharge, SpinMulti))

                for j in range(len(Mol_atom)):
                    ofile.write('%-4s % 10.5f  % 10.5f  % 10.5f \n'
                    % (Mol_atom[j],X[j], Y[j], Z[j]))

            ofile.write('\n')
#####For adding other jobs when other options are required...

            if 1 in option_array[0:17]: 

                ofile.write('--Link1--\n')
##########################################################

        if 1 in option_array[0:17]: 

            ofile.write(line_system)
            ofile.write(line_chk+'\n')
            ofile.write(line_o_method+'\n')

            if option_array[0] == 1: # opt == 1
                ofile.write('Opt=(MaxCycles=100)\n')

           
            if option_array[1] == 1: # opt == 1
                ofile.write('Freq=(Raman)\n')

#######Solvent effect############################
            SCRF, SCRF_read = self.MakeSolventLine()
            ofile.write(SCRF)

#####Reading Geometry and MO from Checkpoint file
            if ReadFromchk == 1:
                ofile.write(line_readMOGeom+'\n')
#################################################

            ofile.write('\n')

#####Reading Geometry from sdf file#####

            if ReadFromsdf == 1 or ReadFromxyz == 1:

                ofile.write(line_comment+'\n')
                ofile.write('\n')

                ofile.write('%5d %5d \n' % (TotalCharge, SpinMulti))

                for j in range(len(Mol_atom)):
                    ofile.write('%-4s % 10.5f  % 10.5f  % 10.5f \n'
                     % (Mol_atom[j],X[j], Y[j], Z[j]))

                ofile.write('\n')
                ofile.write(SCRF_read)
#######################################

            if option_array[2] == 1: # nmr == 1

                ofile.write('--Link1--\n')

                ofile.write(line_system)
                ofile.write(line_chk+'\n')
                ofile.write(line_method+'\n')

                line_method_nmr = 'NMR'
                ofile.write(line_method_nmr+'\n')
                ofile.write(SCRF)

                ofile.write(line_readMOGeom+'\n')
                ofile.write('\n')
                ofile.write(SCRF_read)

            if option_array[11] == 1: # vip == 1

                ofile.write('--Link1--\n')
                ofile.write(line_system)

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

            if option_array[12] == 1: # vea == 1

                ofile.write('--Link1--\n')
                ofile.write(line_system)
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

            if option_array[13] == 1: # aip == 1

                ofile.write('--Link1--\n')
                ofile.write(line_system)
                ofile.write(line_chk+'_PC\n')
                ofile.write(line_oldchk+'_IP\n')
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
                ofile.write(line_system)
                ofile.write(line_chk+'_PC\n')
                #ofile.write(line_oldchk+'\n')
                ofile.write(line_o_method+'\n')
                ofile.write(SCRF)

                ofile.write(line_readGeom+'\n')
                ofile.write('\n')

                ofile.write('neutrization energy calculation from a cation\n')
                ofile.write('\n')

                ofile.write('%5d %5d \n' % (TotalCharge, SpinMulti))
                ofile.write('\n')
                ofile.write(SCRF_read)

            if option_array[14] == 1: #aea == 1

                ofile.write('--Link1--\n')
                ofile.write(line_system)
                ofile.write(line_oldchk+'_EA\n')
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
                ofile.write(line_system)
                #ofile.write(line_oldchk+'\n')
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

            if option_array[3] == 1 : #uv == 1 or fluor==1 or tadf == 1
                sTD = self.MakeLinkTD(line_chk, line_method, None, SCRF, SCRF_read, False, targetstate )
                ofile.write(sTD) 
            if option_array[9] == 1 :
                sTD = self.MakeLinkTD(line_chk, line_c_method, 'Singlet',  SCRF, SCRF_read, True, targetstate )
                ofile.write(sTD) 
            if option_array[10] == 1 :
                sTD = self.MakeLinkTD(line_chk, line_c_method, 'Triplet',  SCRF, SCRF_read, True, targetstate )
                ofile.write(sTD) 

            ofile.write('\n') 

        ofile.close()

############Run Gaussian##############################

        job_state = "normal"
        output_dic = {}

        if os.path.isdir(PreGauInput[0]):
            shutil.rmtree(PreGauInput[0])

        os.mkdir(PreGauInput[0])
        shutil.move(GauInputName, PreGauInput[0])

        if ReadFromchk == 1 :
            inchkfile = PreGauInput[0]+".chk"
            shutil.move(inchkfile, PreGauInput[0]) 

        os.chdir(PreGauInput[0])
        
        job_state = GaussianRunPack.Exe_Gaussian.exe_Gaussian(PreGauInput[0], self.timexe)
        GaussianRunPack.chk2fchk.Get_chklist()
      #  output_dic = self.Extract_values(PreGauInput[0], option_array, Bondpair1, Bondpair2)
        try:
            output_dic = self.Extract_values(PreGauInput[0], option_array, Bondpair1, Bondpair2)
        except:
            job_state = "error"
            pass

#####for pka computation#######################

        if option_array[16] == 1:

            output_dic_pka = {}
            
            E_pH = output_dic["Energy"][0]
            Atom = output_dic["cden"][0]
            MullCharge = output_dic["cden"][1]

            Index_MaxProtic = GaussianRunPack.Get_ChargeSpin.find_MaxProtic(Atom, MullCharge)

            print(Index_MaxProtic)

            GS_fchk = PreGauInput[0]+".fchk"

            with open(GS_fchk,'r') as ifile:
                GS_lines = ifile.readlines()
            

            TotalCharge, SpinMulti, Mol_atom, Mol_X, Mol_Y, Mol_Z = GaussianRunPack.Get_MolCoordinate_fchk.Extract_MolCoord(GS_lines)

            DeHMol_atom = np.delete(Mol_atom,Index_MaxProtic) 
            DeHMol_X = np.delete(Mol_X,Index_MaxProtic)
            DeHMol_Y = np.delete(Mol_Y,Index_MaxProtic)
            DeHMol_Z = np.delete(Mol_Z,Index_MaxProtic)

            JobName_DeHMol = PreGauInput[0] + "_DeH"
            GauInputName_DeHOpt = JobName_DeHMol + ".com"
            GauchkName_DeHOpt = JobName_DeHMol + ".chk"

            ofile_DeHMol = open(GauInputName_DeHOpt, 'w')

            line_DeHchk = '%chk='+ JobName_DeHMol

            ####making input for a deprotonated molecule######
            ofile_DeHMol.write(line_system)
            ofile_DeHMol.write(line_DeHchk + '\n')
            ofile_DeHMol.write(line_o_method+'\n')
            if option_array[0] == 1: # opt == 1
                ofile_DeHMol.write('Opt=(MaxCycles=100)\n')
            ofile_DeHMol.write(SCRF)

            ofile_DeHMol.write("\n")
            ofile_DeHMol.write("Deprotonated molecule \n")
            ofile_DeHMol.write("\n")

            ofile_DeHMol.write('%5d %5d \n' % (TotalCharge-1, SpinMulti))

            for j in range(len(DeHMol_atom)):
                ofile_DeHMol.write('%-4s % 10.5f  % 10.5f  % 10.5f \n'
                 % (DeHMol_atom[j], DeHMol_X[j], DeHMol_Y[j], DeHMol_Z[j]))

            ofile_DeHMol.write(SCRF_read)
            ofile_DeHMol.write('\n')

            ofile_DeHMol.close()

            job_state = GaussianRunPack.Exe_Gaussian.exe_Gaussian(JobName_DeHMol, self.timexe)
            GaussianRunPack.chk2fchk.Get_chklist()
#            output_dic_pka = self.Extract_values(JobName_DeHMol, option_array_pka, Bondpair1, Bondpair2)
            try:
                output_dic_pka = self.Extract_values(JobName_DeHMol, option_array_pka, Bondpair1, Bondpair2)
            except:
                job_state = "error"
                pass

            #print("pka: ", output_dic_pka)
            E_dH = output_dic_pka["Energy"][0]

            output_dic["pka"] = (E_dH - E_pH)*Eh2kJmol

#####for fluor == 1 or tadf == 1 for open shell#############
        if option_array_Ex[9] == 1 or option_array_Ex[10] == 1: 
    
            output_dic_Ex = {}
        
            compute_state = output_dic["state_index"][0][int(targetstate)-1] 
            JobName_ExOpt = PreGauInput[0] + "_ExOpt"
            GauInputName_ExOpt = JobName_ExOpt + ".com"
            ofile_ExOpt = open(GauInputName_ExOpt ,'w')

            ofile.write(line_system)
    
            ofile_ExOpt.write(line_readMOGeom+'\n')
            ofile_ExOpt.write('\n')
    
            sTD = self.MakeLinkTD(line_chk, line_o_method, None, SCRF, SCRF_read, True, compute_state)
            ofile_ExOpt.write(sTD)
    
            if option_array_Ex[10] == 1: #tadf == 1
                compute_state = output_dic["state_index"][1][int(targetstate)-1] 
                sTD = self.MakeLinkTD(line_chk, line_o_method, None, SCRF, SCRF_read, True, compute_state)
                ofile_ExOpt.write(sTD) 
    
            ofile_ExOpt.close()

            job_state = GaussianRunPack.Exe_Gaussian.exe_Gaussian(JobName_ExOpt, self.timexe)
            GaussianRunPack.chk2fchk.Get_chklist()
            try:
                output_dic_Ex = self.Extract_values(JobName_ExOpt, option_array_Ex, Bondpair1, Bondpair2)
            except:
                job_state = "error"
                pass


            output_dic.update(output_dic_Ex)

        output_dic["log"] = job_state


###Convert fchk to xyz ###########################################
        if self.restart == False:
            GaussianRunPack.Get_MolCoordinate_fchk.Get_fchklist2xyz()
##################################################################

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


        
