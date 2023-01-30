import os
import shutil

import numpy as np

from qcforever import gaussian_run
from qcforever.util import read_mol_file


Eh2kJmol = 2625.5
Eh2eV = 27.211


class GaussianDFTRun:
    def __init__(self, functional, basis, nproc, value, in_file, solvent="0", error=0, restart=True):
        self.in_file = in_file
        self.functional = functional.lower()
        self.basis = basis.lower()
        self.nproc = nproc
        self.value = value
        self.solvent = solvent.lower()
        self.restart = restart
        self.error = error
        self.mem = ''
        self.timexe = 60 * 60 * 80
        self.SpecTotalCharge = np.nan
        self.SpecSpinMulti = np.nan
        self.ref_uv_path = ''
        self.para_functional = []

    def Extract_SCFEnergy(self, lines):
        Energy = []
        for line in lines:
            if line.find("SCF Done:  ") >= 0:
                line_StateInfo = line.split()
                Energy.append(float(line_StateInfo[4]))
        Comp_SS, Ideal_SS = gaussian_run.Estimate_SpinContami.Estimate_SpinDiff(lines)
        Energy_Spin = [Energy[-1], Comp_SS-Ideal_SS]
        return Energy_Spin

    def Extract_values(self, jobname, option_dict, Bondpair1, Bondpair2):
        """
        MEMO: Bondpair1 and Bondpair2 are not used
        """
        is_opt = option_dict['opt'] if 'opt' in option_dict else False
        is_freq = option_dict['freq'] if 'freq' in option_dict else False
        is_nmr = option_dict['nmr'] if  'nmr' in option_dict else False
        is_uv = option_dict['uv'] if 'uv' in option_dict else False
        is_energy = option_dict['energy'] if 'energy' in option_dict else False
        is_homolumo = option_dict['homolumo'] if 'homolumo' in option_dict else False
        is_dipole = option_dict['dipole'] if 'dipole' in option_dict else False
        is_deen = option_dict['deen'] if 'deen' in option_dict else False
        is_stable2o2 = option_dict['stable2o2'] if 'stable2o2' in option_dict else False
        is_fluor = option_dict['fluor'] if 'fluor' in option_dict else False
        is_tadf = option_dict['tadf'] if 'tadf' in option_dict else False
        is_vip = option_dict['vip'] if 'vip' in option_dict else False
        is_vea = option_dict['vea'] if 'vea' in option_dict else False
        is_aip = option_dict['aip'] if 'aip' in option_dict else False
        is_aea = option_dict['aea'] if 'aea' in option_dict else False
        is_cden = option_dict['cden'] if 'cden' in option_dict else False
        is_pka = option_dict['pka']  if 'pka' in option_dict else False # not used
        is_satkoopmans = option_dict['satkoopmans'] if 'satkoopmans' in option_dict else False
        is_symm = option_dict['symm'] if 'symm' in option_dict else False
        infilename = f"{jobname}.log"
        fchkname = f"{jobname}.fchk"

        output = {}
        with open(infilename, 'r') as ifile:
            lines = ifile.readlines()
        print("Spliting links....")
        Links = self.SplitLinks(infilename)
        n = len(Links)
        print(f"The number of linkes = {n}")

        # Clean lines
        lines = ""

        if len(option_dict.keys()) == 1 and option_dict.keys() == 'symm': 
            scf_done = False 
        else:
            scf_done = True

        # For extracting symmetry of molecule without SCF
        if is_symm:
            Symm_lines = Links[is_symm].splitlines()
            pGroup = "C1"
            for line in Symm_lines:
                if line.find("Full point group  ") >= 0:
                    line_symmInfo = line.split()
                    pGroup = line_symmInfo[3]
            output["symm"] = pGroup

        # For extracting properties of molecule after SCF        
        if scf_done: 
            GS_lines = Links[1+is_symm].splitlines()
            Links[1+is_symm] = ""

        if is_opt:
            print ("Optimization was performed...Check geometry...")
            MaxDisplace = gaussian_run.Get_MolInterCoordinate.Extract_InterMol(GS_lines)
            output["GS_MaxDisplace"] = MaxDisplace

        if is_freq:
            print("For getting frequency")
            Freq, IR, Raman, E_zp,  E_t, E_enth, E_free, Ei, Cv, St = gaussian_run.Get_FreqPro.Extract_Freq(GS_lines) 
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

        if is_homolumo:
            with open(fchkname, 'r') as ifile:
                fchk_lines = ifile.readlines()
            NumAlphaElec, NumBetaElec, AlphaEigenVal, BetaEigenVal = gaussian_run.Get_MOEnergy_fchk.Extract_MO(fchk_lines)
            # NumAlphaElec, NumBetaElec, AlphaEigenVal, BetaEigenVal = gaussian_run.Get_MOEnergy.Extract_MO(GS_lines)
            if BetaEigenVal == []:
                Alpha_gap = Eh2eV * (AlphaEigenVal[NumAlphaElec]-AlphaEigenVal[NumAlphaElec-1])
                output["homolumo"] = Alpha_gap
                # return Alpha_gap
            else:
                Alpha_gap = Eh2eV * (AlphaEigenVal[NumAlphaElec]-AlphaEigenVal[NumAlphaElec-1])
                Beta_gap = Eh2eV * (BetaEigenVal[NumBetaElec]-BetaEigenVal[NumBetaElec-1])
                output["homolumo"] = [Alpha_gap, Beta_gap]
                # return Alpha_gap, Beta_gap

        if is_dipole:
            Dipole_X = []
            Dipole_Y = []
            Dipole_Z = []
            Dipole_Total = []
            for line in GS_lines:
                if line.find(" X= ") >= 0:
                    line_StateInfo = line.split()
                    Dipole_X.append(float(line_StateInfo[1]))
                    Dipole_Y.append(float(line_StateInfo[3]))
                    Dipole_Z.append(float(line_StateInfo[5]))
                    Dipole_Total.append(float(line_StateInfo[7]))
            output["dipole"] = [Dipole_X[-1], Dipole_Y[-1], Dipole_Z[-1], Dipole_Total[-1]]
            # return Dipole_X[-1], Dipole_Y[-1], Dipole_Z[-1], Dipole_Total[-1]

        if is_energy:
            output["Energy"] = self.Extract_SCFEnergy(GS_lines)

        if is_deen:
            try:
                GS_Energy = output["Energy"][0]
            except KeyError:
                output["Energy"] = self.Extract_SCFEnergy(GS_lines)
                GS_Energy = output["Energy"][0] 
            Mol_atom, _, _, _ = gaussian_run.Get_MolCoordinate.Extract_Coordinate(GS_lines)
            # Calculating Decomposed atoms total energy
            decomposed_Energy = 0
            for i in range(len(Mol_atom)):    
                decomposed_Energy += gaussian_run.AtomInfo.One_Atom_Energy(Mol_atom[i], self.functional, self.basis)
            print("Decomposed energy: ", decomposed_Energy)
            # return deen 
            output["deen"] = GS_Energy - (decomposed_Energy)

        if is_stable2o2:
            try:
                NumAlphaElec
            except:
                with open(fchkname, 'r') as fchkfile:
                    fchk_lines = fchkfile.readlines()
                NumAlphaElec, NumBetaElec, AlphaEigenVal, BetaEigenVal = gaussian_run.Get_MOEnergy_fchk.Extract_MO(fchk_lines)
                # NumAlphaElec, NumBetaElec, AlphaEigenVal, BetaEigenVal = gaussian_run.Get_MOEnergy.Extract_MO(GS_lines)
            O2_SOMO, O2_LUMO = gaussian_run.AtomInfo.O2_MO_refer(self.functional, self.basis)
            if BetaEigenVal == []:
                """
                OxidizedbyO2  >  0 oxidation by O2 is hard to occure.
                OxidizedbyO2 <=  0 oxidation by O2 is easy to occure.
                ReducedbyO2  >  0 reduction by O2 is hard to occure.
                ReducedbyO2 <=  0 reduction by O2 is easy to occure.
                """
                OxidizedbyO2 = O2_LUMO -AlphaEigenVal[NumAlphaElec-1]
                ReducedbyO2 =  AlphaEigenVal[NumAlphaElec] - O2_SOMO
                output["stable2o2"] = [OxidizedbyO2,ReducedbyO2]
            else:
                OxidizedbyO2 = O2_LUMO - BetaEigenVal[NumBetaElec-1] 
                ReducedbyO2 = AlphaEigenVal[NumAlphaElec] - O2_SOMO
                output["stable2o2"] = [OxidizedbyO2,ReducedbyO2]

        if is_cden:
            output["cden"] = gaussian_run.Get_ChargeSpin.Extract_ChargeSpin(GS_lines)

        """ Index of Links
        Index = 0 : (always blank)
        Index = 1+symm : Ground state [if opt==1]
        Index = 1+symm+nmr : NMR chemical shift of S0 [if nmr==1]
        Index = 1+symm+nmr+vip : Ionization potential [if vip==1]
        Index = 1+symm+nmr+vip+vea : Electronic affinity [if vea==1]
        Index = 1+symm+nmr+vip+vea+aip : adiabatic ionization potential [if aip==1]
        Index = 1+symm+nmr+vip+vea+aip+aea : adiabatic electronic affinity [if aea==1]
        Index = 1+symm+nmr+vip+vea+aip+aea+1 : Virtical excitation (S0 -> S1) [uv]
        Index = 2+symm+nmr+vip+vea+aip+aea+1 : Optimization of S1 [fluor or tadf] 
        Index = 3+symm+nmr+vip+vea+aip+aea+1 : Optimization of T1 [tadf]
        """
        if is_nmr:
            Element = []
            ppm = []
            Index = 1 + is_symm + is_nmr
            nmr_lines = Links[Index].splitlines()
            Links[Index] = ""
            for line in nmr_lines:
                if line.find("Isotropic =  ") >= 0:
                    line_Info = line.split()
                    Element.append(line_Info[1])
                    ppm.append(float(line_Info[4]))
            # calculating chemical shift for H, C, or Si
            for i in range(len(Element)):
                if Element[i]=="H" or Element[i]=="C" or Element[i]=="Si":
                    ppm[i] = gaussian_run.AtomInfo.One_TMS_refer(Element[i], self.functional, self.basis) - ppm[i]
            # return Element, ppm 
            output["nmr"] = [Element, ppm]

        if is_vip or is_vea:
            try:
                GS_Energy = output["Energy"][0]
            except KeyError:
                output["Energy"] = self.Extract_SCFEnergy(GS_lines)
                GS_Energy = output["Energy"][0] 
           # print (GS_Energy)
            if is_vip:
                Index = 1 + is_symm + is_nmr + is_vip
                IP_lines = Links[Index].splitlines()
                Links[Index] = ""
                IP_Energy_SS = self.Extract_SCFEnergy(IP_lines)
                # IP_Comp_SS, IP_Ideal_SS = gaussian_run.Estimate_SpinContami.Estimate_SpinDiff(IP_lines)
                # Normal ionization potential calculation
                output["vip"] = [Eh2eV*(IP_Energy_SS[0]-GS_Energy), IP_Energy_SS[1]]
            if is_vea:
                Index = 1 + is_symm + is_nmr + is_vip + is_vea
                EA_lines = Links[Index].splitlines()
                Links[Index] = ""
                EA_Energy_SS =  self.Extract_SCFEnergy(EA_lines)
                # EA_Comp_SS, EA_Ideal_SS = gaussian_run.Estimate_SpinContami.Estimate_SpinDiff(EA_lines)
                # Normal electronic affinity calculation
                output["vea"] = [Eh2eV*(GS_Energy-EA_Energy_SS[0]), EA_Energy_SS[1]]

        if is_aip or is_aea:
            try:
                GS_Energy = output["Energy"][0]
            except KeyError:
                output["Energy"] = self.Extract_SCFEnergy(GS_lines)
                GS_Energy = output["Energy"][0] 
            if is_aip:
                Index = 1 + is_symm + is_nmr + is_vip + is_vea + is_aip
                PC_lines = Links[Index].splitlines()
                Links[Index] = ""
                VNP_lines = Links[Index+1].splitlines()
                Links[Index+1] = ""
                is_aip += 1
                PC_Energy_SS = self.Extract_SCFEnergy(PC_lines)
                # PC_Comp_SS, PC_Ideal_SS = gaussian_run.Estimate_SpinContami.Estimate_SpinDiff(PC_lines)
                VNP_Energy_SS = self.Extract_SCFEnergy(VNP_lines)
                # VNP_Comp_SS, VNP_Ideal_SS = gaussian_run.Estimate_SpinContami.Estimate_SpinDiff(VNP_lines)
                # For Check internal coordinate
                MaxDisplace = gaussian_run.Get_MolInterCoordinate.Extract_InterMol(PC_lines)
                output["relaxedIP_MaxDisplace"] = MaxDisplace
                output["aip"] = [Eh2eV*(PC_Energy_SS[0]-GS_Energy), Eh2eV*(PC_Energy_SS[0]-VNP_Energy_SS[0]), PC_Energy_SS[1], VNP_Energy_SS[1]]
            if is_aea:
                Index = 1 + is_symm + is_nmr + is_vip + is_vea + is_aip + is_aea 
                NC_lines = Links[Index].splitlines()
                Links[Index] = ""
                VNN_lines = Links[Index+1].splitlines()
                Links[Index+1] = ""
                is_aea += 1
                NC_Energy_SS = self.Extract_SCFEnergy(NC_lines)
                # NC_Comp_SS, NC_Ideal_SS = gaussian_run.Estimate_SpinContami.Estimate_SpinDiff(NC_lines)
                VNN_Energy_SS = self.Extract_SCFEnergy(VNN_lines)
                # VNN_Comp_SS, NC_Ideal_SS = gaussian_run.Estimate_SpinContami.Estimate_SpinDiff(VNN_lines)
                # For Check internal coordinate
                MaxDisplace = gaussian_run.Get_MolInterCoordinate.Extract_InterMol(NC_lines)
                output["relaxedEA_MaxDisplace"] = MaxDisplace
                # Normal electronic affinity calculation
                output["aea"] = [Eh2eV*(GS_Energy-NC_Energy_SS[0]), Eh2eV*(VNN_Energy_SS[0]-NC_Energy_SS[0]), NC_Energy_SS[1], VNN_Energy_SS[1]]

        if is_satkoopmans:
            dSCF_vip = output["vip"][0]
            dSCF_vea = output["vea"][0]
            Alpha_eHOMO = AlphaEigenVal[NumAlphaElec-1] 
            Beta_eHOMO = BetaEigenVal[NumBetaElec-1]
            Alpha_eLUMO = AlphaEigenVal[NumAlphaElec] 
            Beta_eLUMO = BetaEigenVal[NumBetaElec]
            Delta_vip = dSCF_vip + Eh2eV*max(Alpha_eHOMO, Beta_eHOMO)
            Delta_vea = dSCF_vea + Eh2eV*min(Alpha_eLUMO, Beta_eLUMO)
            output["satkoopmans"] = [Delta_vip, Delta_vea]
  
        if is_uv:
            Index = 1 + is_symm + is_nmr + is_vip + is_vea + is_aip + is_aea + 1 
            lines = "" if Index >= n else Links[Index].splitlines()
            _, _, _, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
                CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden = gaussian_run.Get_ExcitedState.Extract_ExcitedState(lines)
            output["uv"] = [WL_allowed, OS_allowed, CD_L_allowed, CD_OS_allowed]
            output["state_index"] = [State_allowed, State_forbidden]
            Links[Index] = ""
            if self.ref_uv_path != '':
                ref_uv = {}
                print(f"Read the reference spectrum...{self.ref_uv_path}")
                ref_uv["uv"] = gaussian_run.UV_similarity.read_data(self.ref_uv_path)
                S, D = gaussian_run.UV_similarity.smililarity_dissimilarity(ref_uv["uv"][0], ref_uv["uv"][1], output["uv"][0], output["uv"][1])
                output["Similality/Disdimilarity"] = [S, D]
  
        if is_fluor or is_tadf:
            Index = 1 + is_symm + is_uv + is_nmr + is_vip + is_vea + is_aip + is_aea + is_fluor 
            lines = "" if Index >= n else Links[Index].splitlines()
            S_Found, S_Egrd, S_Eext, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
                CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden = gaussian_run.Get_ExcitedState.Extract_ExcitedState(lines)
            # For Check internal coordinate
            MaxDisplace = gaussian_run.Get_MolInterCoordinate.Extract_InterMol(lines)
            output["MinEtarget"] = S_Eext
            output["Min_MaxDisplace"] = MaxDisplace
            output["fluor"] = [WL_allowed, OS_allowed, CD_L_allowed, CD_OS_allowed]
            Links[Index] = ""
  
        if is_tadf:
            Index = 1 + is_symm + is_uv + is_nmr + is_vip + is_vea + is_aip + is_aea + is_fluor + is_tadf
            lines = "" if Index >= n else Links[Index].splitlines()
            T_Found, _, T_Eext, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
                CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden  = gaussian_run.Get_ExcitedState.Extract_ExcitedState(lines)
            # For Check internal coordinate
            MaxDisplace = gaussian_run.Get_MolInterCoordinate.Extract_InterMol(lines)
            output["T_Min"] = T_Eext
            output["T_Min_MaxDisplace"] = MaxDisplace
            output["T_Phos"] = [WL_forbidden, OS_forbidden, CD_L_forbidden, CD_OS_forbidden]
            TADF_Eng = 0.0
            if S_Found and T_Found:
               TADF_Eng = S_Eext - T_Eext
            output["Delta(S-T)"] = TADF_Eng
            Links[Index] = ""
        del lines
        del Links
        return output

    def MakeSolventLine(self):
        s = ''
        s_solvent = ''
        try:
            float(self.solvent)
        except ValueError:
            print('Solvent effect is included by PCM')
            s += f'SCRF=(PCM, solvent={self.solvent})\n'
        else:
            if self.solvent == '0':
                return s, s_solvent
            else:
                print('Solvent effect is included by PCM')
                s += 'SCRF=(PCM, solvent=Generic, Read)\n'
                s_solvent += f'EPS={self.solvent}\n'
                s_solvent += f'Radii=UA0\n\n'
        return s, s_solvent

    def MakeLinkTD(self, line_chk, line_method, State, SCRF, SCRF_read, Opt, targetstate=1):
        self.MakeSolventLine()
        Jobname_line =line_chk.split('=')
        Jobname = Jobname_line[-1]
        # Setting NState######
        NState = targetstate+4 if Opt == True else 20
        line_oldchk = f'%Oldchk={Jobname}'
        # Potaintially cause a bug because '==' is used instead of 'is' to check 'None'. Need to check
        line_newchk = f"%chk={Jobname}_ExOptState{targetstate}" if State == None else f"%chk={Jobname}_ExOptState{State}"
        # Potaintially cause a bug because '==' is used instead of 'is' to check 'None'. Need to check
        line_method_TD = f"TD(Nstate={NState}, root={targetstate})" if State == None else f"TD(Nstate={NState}, {State}, root={targetstate})"
        line_readMOGeom = 'Geom=AllCheck Guess=Read'
        s = ''
        s += '--Link1--\n'
        if self.mem != '':
            line_mem = f'%mem={self.mem}'
            s += f"{line_mem}\n"
        if self.nproc > 1:
            line_proc = f'%nproc={self.nproc}'
            s += f"{line_proc}\n"
        if Opt == True:
            s += f"{line_oldchk}\n"
            s += f"{line_newchk}\n"
        else:
            s += f"{line_chk}\n"
        s += f"{line_method}\n"
        s += f"{line_method_TD}\n"
        s += f"{line_readMOGeom}\n"
        s += SCRF
        if Opt == True:
            s += "Opt=(MaxCycles=50)\n"
            # s = s + "Opt"
        s += '\n' 
        s += SCRF_read
        return s

    def run_gaussian(self):
        infilename = self.in_file
        option_line = self.value    
        options = option_line.split()
        option_dict = {}
        option_dict_Ex = {}
        option_dict_pka = {}
        targetstate = 1
        PreGauInput = infilename.split('.')
        GauInputName = PreGauInput[0]+'.com'    
        # File type of input?
        ReadFromchk = 0 
        ReadFromsdf = 0 
        ReadFromxyz = 0 
        if PreGauInput[1] == "sdf":
            ReadFromsdf = 1 
            Mol_atom, X, Y, Z, TotalCharge, SpinMulti, Bondpair1, Bondpair2 = read_mol_file.read_sdf(infilename)
        elif PreGauInput[1] == "xyz":
            ReadFromxyz = 1
            Mol_atom, X, Y, Z, TotalCharge, SpinMulti = read_mol_file.read_xyz(infilename)
            Bondpair1 = []
            Bondpair2 = []
        elif PreGauInput[1] == "chk":
            ReadFromchk = 1 
            Bondpair1 = []
            Bondpair2 = []
        elif PreGauInput[1] == "fchk":
            TotalCharge, SpinMulti = gaussian_run.fchk2chk.Get_fchk(PreGauInput[0])
            ReadFromchk = 1 
            Bondpair1 = []
            Bondpair2 = []
        else:
            print("Invalid input file")

        # Setting electronic structure
        if np.isnan(self.SpecTotalCharge) !=  True:
            TotalCharge = self.SpecTotalCharge
        if np.isnan(self.SpecSpinMulti) != True:
            SpinMulti = self.SpecSpinMulti

        # Setting options
        for i in range(len(options)):
            option = options[i]
            if option.lower() == 'opt':
                # opt = 1
                option_dict['opt'] = True
            elif option.lower() == 'freq':
                # nmr = 1
                option_dict['freq'] = True
            elif option.lower() == 'nmr':
                # nmr = 1
                option_dict['nmr'] = True
            elif 'uv' in option.lower():
                # uv = 1
                option_dict['uv'] = True
                if '=' in option:
                    ref_spectrum = option.split("=")
                    self.ref_uv_path = ref_spectrum[-1]
            elif option.lower() == 'energy':
                # energy = 1
                option_dict['energy'] = True
            elif option.lower() == 'homolumo':
                # homolumo = 1
                option_dict['homolumo'] = True
            elif option.lower() == 'dipole':
                # dipole = 1
                option_dict['dipole'] = True
            elif option.lower() == 'deen':
                # deen = 1
                option_dict['deen'] = True
            elif option.lower() == 'stable2o2':
                # stable2o2 = 1
                option_dict['stable2o2'] = True
            elif 'fluor' in option.lower():
                # fluor = 1
                option_dict['uv'] = True
                if SpinMulti == 1:
                    option_dict['fluor'] = True
                else:
                    option_dict_Ex['fluor'] = True
                if '=' in option:
                    in_target = option.split("=")
                    targetstate = int(in_target[-1])
            elif option.lower() == 'tadf':
                # tadf = 1
                option_dict['uv'] = True
                if SpinMulti == 1:
                    option_dict['fluor'] = True
                    option_dict['tadf'] = True
                else:
                    option_dict_Ex['fluor'] = True
                    option_dict_Ex['tadf'] = True
            elif option.lower() == 'vip':
                # vip = 1
                option_dict['vip'] = True
            elif option.lower() == 'vea':
                # vea = 1
                option_dict['vea'] = True
            elif option.lower() == 'aip':
                # aip = 1
                option_dict['vip'] = True
                option_dict['aip'] = True
            elif option.lower() == 'aea':
                # aea = 1
                option_dict['vea'] = True
                option_dict['aea'] = True
            elif option.lower() == 'cden':
                # nne = 1
                #print ('Neutraization energy from anion')
                option_dict['cden'] = True
            elif option.lower() == 'pka':
                option_dict['energy'] = True
                option_dict['cden'] = True
                option_dict['pka'] = True
                option_dict_pka['energy'] = True
            elif option.lower() == 'satkoopmans':
                option_dict['homolumo'] = True
                option_dict['vip'] = True
                option_dict['vea'] = True
                option_dict['satkoopmans'] = True
            elif option.lower() == 'symm':
                option_dict['symm'] = True
            else:
                print('invalid option: ', option)


        if len(option_dict.keys()) == 1 and option_dict.keys() == 'symm': 
            scf_needed = False 
        else:
            scf_needed = True

        line_system = ''
        if self.mem !='' or self.nproc > 1:
            if self.mem != '':
                line_mem = f'%mem={self.mem}\n'
                line_system += line_mem
            if self.nproc > 1 :
                line_proc = f'%nproc={self.nproc}\n'
                line_system += line_proc
        line_chk = f'%chk={PreGauInput[0]}'
        line_oldchk = f'%Oldchk={PreGauInput[0]}'
        line_iop_functional = ''
        if self.para_functional == []:
            pass
        else:
            line_iop_functional += gaussian_run.Gen_IOPline_functional.functional_para(self.functional, self.para_functional)
        line_method = f'#{self.functional}/{self.basis} test {line_iop_functional} '
        line_o_method = f'#u{self.functional}/{self.basis} test {line_iop_functional} '
        line_c_method = f'#r{self.functional}/{self.basis} test {line_iop_functional} '
        line_comment = infilename

        ofile = open(GauInputName ,'w')
        # For reading geomerty and MO from checkpoint file
        line_readMOGeom = 'Geom=AllCheck Guess=Read'
        line_readOnlyMOGeom = 'Geom=CheckPoint Guess=Read'
        line_readGeom = 'Geom=Checkpoint'

        if 'symm' in option_dict:
            ofile.write(line_system)
            ofile.write(line_chk+'\n')
            ofile.write(line_o_method+'\n')
            ofile.write('Guess=Only Symmetry=loose'+'\n')
           # ofile.write('Guess=Only Symmetry=on'+'\n')
        # Reading Geometry and MO from Checkpoint file
            if ReadFromchk == 1:
                ofile.write(line_readMOGeom+'\n')
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
            # For adding other jobs when other options are required...
            if scf_needed: 
                ofile.write('--Link1--\n')

        if scf_needed:
            ofile.write(line_system)
            ofile.write(line_chk+'\n')
            ofile.write(line_o_method+'\n')

            if 'opt' in option_dict: # opt == 1
                ofile.write('Opt=(MaxCycles=100)\n')

            if 'freq' in option_dict: # freq == 1
                ofile.write('Freq=(Raman)\n')
            # Solvent effect
            SCRF, SCRF_read = self.MakeSolventLine()
            ofile.write(SCRF)

            # Reading Geometry and MO from Checkpoint file
            if ReadFromchk == 1:
                ofile.write(line_readMOGeom+'\n')
            ofile.write('\n')

            # Reading Geometry from sdf or xyz file
            if ReadFromsdf == 1 or ReadFromxyz == 1:
                ofile.write(line_comment+'\n')
                ofile.write('\n')
                ofile.write('%5d %5d \n' % (TotalCharge, SpinMulti))
                for j in range(len(Mol_atom)):
                    ofile.write('%-4s % 10.5f  % 10.5f  % 10.5f \n'
                     % (Mol_atom[j],X[j], Y[j], Z[j]))
                ofile.write('\n')
                ofile.write(SCRF_read)
            #
            if 'nmr' in option_dict: # nmr == 1
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

            if 'vip' in option_dict: # vip == 1
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
                if SpinMulti == 1:
                    ofile.write('%5d %5d \n' % (TotalCharge+1, SpinMulti+1))
                else:
                    ofile.write('%5d %5d \n' % (TotalCharge+1, SpinMulti-1))
                ofile.write('\n')
                ofile.write(SCRF_read)

            if 'vea' in option_dict: # vea == 1
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
                if SpinMulti == 1:
                    ofile.write('%5d %5d \n' % (TotalCharge-1, SpinMulti+1))
                else:
                    ofile.write('%5d %5d \n' % (TotalCharge-1, SpinMulti-1))
                ofile.write('\n')
                ofile.write(SCRF_read)

            if 'aip' in option_dict: # aip == 1
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
                if SpinMulti == 1:
                    ofile.write('%5d %5d \n' % (TotalCharge+1, SpinMulti+1))
                else:
                    ofile.write('%5d %5d \n' % (TotalCharge+1, SpinMulti-1))
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

            if 'aea' in  option_dict: #aea == 1
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
                if SpinMulti == 1:
                    ofile.write('%5d %5d \n' % (TotalCharge-1, SpinMulti+1))
                else:
                    ofile.write('%5d %5d \n' % (TotalCharge-1, SpinMulti-1))
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

            if 'uv' in option_dict: #uv == 1 or fluor==1 or tadf == 1
                sTD = self.MakeLinkTD(line_chk, line_method, None, SCRF, SCRF_read, False, targetstate )
                ofile.write(sTD) 

            if 'fluor' in option_dict:
                sTD = self.MakeLinkTD(line_chk, line_c_method, 'Singlet',  SCRF, SCRF_read, True, targetstate )
                ofile.write(sTD) 

            if 'tadf' in option_dict:
                sTD = self.MakeLinkTD(line_chk, line_c_method, 'Triplet',  SCRF, SCRF_read, True, targetstate )
                ofile.write(sTD) 
            ofile.write('\n') 
        ofile.close()

        # Run Gaussian
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
        job_state = gaussian_run.Exe_Gaussian.exe_Gaussian(PreGauInput[0], self.timexe)
        gaussian_run.chk2fchk.Get_chklist()
        #output_dic = self.Extract_values(PreGauInput[0], option_dict, Bondpair1, Bondpair2)
        try:
            output_dic = self.Extract_values(PreGauInput[0], option_dict, Bondpair1, Bondpair2)
        except Exception as e:
            job_state = "error"
            print(e)
            pass
        
        # for pka computation
        if 'pka' in option_dict:
            output_dic_pka = {}
            E_pH = output_dic["Energy"][0]
            Atom = output_dic["cden"][0]
            MullCharge = output_dic["cden"][1]
            Index_MaxProtic = gaussian_run.Get_ChargeSpin.find_MaxProtic(Atom, MullCharge)
            print(Index_MaxProtic)
            GS_fchk = PreGauInput[0]+".fchk"
            with open(GS_fchk,'r') as ifile:
                GS_lines = ifile.readlines()
            TotalCharge, SpinMulti, Mol_atom, Mol_X, Mol_Y, Mol_Z = gaussian_run.Get_MolCoordinate_fchk.Extract_MolCoord(GS_lines)
            DeHMol_atom = np.delete(Mol_atom,Index_MaxProtic) 
            DeHMol_X = np.delete(Mol_X,Index_MaxProtic)
            DeHMol_Y = np.delete(Mol_Y,Index_MaxProtic)
            DeHMol_Z = np.delete(Mol_Z,Index_MaxProtic)
            JobName_DeHMol = PreGauInput[0] + "_DeH"
            GauInputName_DeHOpt = JobName_DeHMol + ".com"
            GauchkName_DeHOpt = JobName_DeHMol + ".chk"

            ofile_DeHMol = open(GauInputName_DeHOpt, 'w')
            line_DeHchk = '%chk='+ JobName_DeHMol
            # making input for a deprotonated molecule
            ofile_DeHMol.write(line_system)
            ofile_DeHMol.write(line_DeHchk + '\n')
            ofile_DeHMol.write(line_o_method+'\n')
            if 'opt' in option_dict: # opt == 1
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

            job_state = gaussian_run.Exe_Gaussian.exe_Gaussian(JobName_DeHMol, self.timexe)
            gaussian_run.chk2fchk.Get_chklist()
            # output_dic_pka = self.Extract_values(JobName_DeHMol, option_dict_pka, Bondpair1, Bondpair2)
            try:
                output_dic_pka = self.Extract_values(JobName_DeHMol, option_dict_pka, Bondpair1, Bondpair2)
            except Exception as e:
                job_state = "error"
                print (e)
                pass
            #print("pka: ", output_dic_pka)
            E_dH = output_dic_pka["Energy"][0]
            output_dic["pka"] = (E_dH - E_pH)*Eh2kJmol

        # for fluor == 1 or tadf == 1 for open shell
        if 'fluor' in option_dict_Ex or 'tadf' in option_dict_Ex: 
            TotalCharge, SpinMulti = gaussian_run.fchk2chk.Get_fchk(PreGauInput[0])
            output_dic_Ex = {}
            compute_state = output_dic["state_index"][0][int(targetstate)-1] 
            JobName_ExOpt = PreGauInput[0] + "_ExOpt"
            GauInputName_ExOpt = JobName_ExOpt + ".com"

            ofile_ExOpt = open(GauInputName_ExOpt ,'w')
            ofile_ExOpt.write(line_system)
            ofile_ExOpt.write(line_chk+'\n')
            ofile_ExOpt.write(line_o_method+'\n')
            ofile_ExOpt.write(line_readMOGeom+'\n')
            ofile_ExOpt.write('\n')
            sTD = self.MakeLinkTD(line_chk, line_o_method, None, SCRF, SCRF_read, True, compute_state)
            ofile_ExOpt.write(sTD)
            if 'tadf' in option_dict_Ex: #tadf == 1
                compute_state = output_dic["state_index"][1][int(targetstate)-1] 
                sTD = self.MakeLinkTD(line_chk, line_o_method, None, SCRF, SCRF_read, True, compute_state)
                ofile_ExOpt.write(sTD) 
            ofile_ExOpt.close()

            job_state = gaussian_run.Exe_Gaussian.exe_Gaussian(JobName_ExOpt, self.timexe)
            gaussian_run.chk2fchk.Get_chklist()
            # output_dic_Ex = self.Extract_values(JobName_ExOpt, option_dict_Ex, Bondpair1, Bondpair2)
            try:
                output_dic_Ex = self.Extract_values(JobName_ExOpt, option_dict_Ex, Bondpair1, Bondpair2)
            except Exception as e:
                job_state = "error"
                print (e)
                pass
            output_dic.update(output_dic_Ex)
        output_dic["log"] = job_state

        # Convert fchk to xyz 
        if self.restart == False:
            gaussian_run.Get_MolCoordinate_fchk.Get_fchklist2xyz()
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
