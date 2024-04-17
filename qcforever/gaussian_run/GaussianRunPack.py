import os
import glob
import sys
import copy
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
        self.geom_spec = {}
        self.para_functional = []

    def check_scf_need(self, option):
        all_keys = set(option.keys())
        nonSCF_keys = ['symm', 'volume']
        if len(all_keys - set(nonSCF_keys)) >= 1:
            scf_need = True
        else:
            scf_need = False

        return scf_need

    def Extract_values(self, jobname, option_dict, Bondpair1, Bondpair2):
        """
        MEMO: Bondpair1 and Bondpair2 are not used
        """
        is_opt = option_dict['opt'] if 'opt' in option_dict else False
        is_polar = option_dict['polar'] if 'polar' in option_dict else False
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
        is_volume = option_dict['volume'] if 'volume' in option_dict else False
        infilename = f"{jobname}.log"
        fchkname = f"{jobname}.fchk"

        print (option_dict)
        
        output = {}

        parse_log = gaussian_run.parse_log.parse_log(infilename)
        job_index, Links_split = parse_log.Check_task()

        # For extracting properties of molecule after SCF        
        if 'ts' in job_index:
            GS_lines = Links_split[job_index['ts']]
        
        if is_opt:
            print ("Optimization was performed...Check geometry...")
            MaxDisplace = gaussian_run.Get_MolInterCoordinate.Extract_InterMol(GS_lines)
            output["GS_MaxDisplace"] = MaxDisplace
        
        if is_homolumo:
            with open(fchkname, 'r') as ifile:
                fchk_lines = ifile.readlines()
            NumAlphaElec, NumBetaElec, AlphaEigenVal, BetaEigenVal = gaussian_run.Get_MOEnergy_fchk.Extract_MO(fchk_lines)
            # NumAlphaElec, NumBetaElec, AlphaEigenVal, BetaEigenVal = gaussian_run.Get_MOEnergy.Extract_MO(GS_lines)
            if BetaEigenVal == []:
                Alpha_gap = Eh2eV * (AlphaEigenVal[NumAlphaElec]-AlphaEigenVal[NumAlphaElec-1])
                output["homolumo"] = Alpha_gap
            else:
                Alpha_gap = Eh2eV * (AlphaEigenVal[NumAlphaElec]-AlphaEigenVal[NumAlphaElec-1])
                Beta_gap = Eh2eV * (BetaEigenVal[NumBetaElec]-BetaEigenVal[NumBetaElec-1])
                output["homolumo"] = [Alpha_gap, Beta_gap]
        
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
        
        if is_energy:
            output["Energy"] = parse_log.Extract_SCFEnergy(GS_lines)
        
        if is_deen:
            try:
                GS_Energy = output["Energy"][0]
            except KeyError:
                output["Energy"] = parse_log.Extract_SCFEnergy(GS_lines)
                GS_Energy = output["Energy"][0] 
            Mol_atom, _, _, _ = gaussian_run.Get_MolCoordinate.Extract_Coordinate(GS_lines)
            # Calculating Decomposed atoms total energy
            decomposed_Energy = 0
            for i in range(len(Mol_atom)):    
                decomposed_Energy += gaussian_run.AtomInfo.One_Atom_Energy(Mol_atom[i], self.functional, self.basis)
            #print("Decomposed energy: ", decomposed_Energy)
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
        
        if is_symm:
            Symm_lines = Links_split[job_index['symm']]
            output["symm"] = parse_log.Extract_symm(Symm_lines)
        
        if is_volume:
            Volume_lines = Links_split[job_index['volume']]
            output["volume"] = parse_log.Extract_volume(Volume_lines)
        
        if is_polar:
            polar_lines = Links_split[job_index['polar_line']]
            #Polar_iso, Polar_aniso = gaussian_run.Get_FreqPro.Extract_polar(polar_lines) 
            Polar_tens, Polar_iso, Polar_aniso = parse_log.extract_polar(polar_lines) 
            output["polar_tens"] = Polar_tens
            output["polar_iso"] = Polar_iso[1]
            output["polar_aniso"] = Polar_aniso[1]
        
        if is_freq:
            freq_lines = Links_split[job_index['freq_line']]
            Freq, IR, Raman, E_zp,  E_t, E_enth, E_free, Ei, Cv, St = gaussian_run.Get_FreqPro.Extract_Freq(freq_lines) 
            VibX, VibY, VibZ = gaussian_run.Get_FreqPro.Extract_vibvec(freq_lines) 
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
            output["freqmode"] = [VibX, VibY, VibZ]
        
        if is_nmr:
            Element = []
            ppm = []
            nmr_lines = Links_split[job_index['nmr_line']]
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
                output["Energy"] = parse_log.Extract_SCFEnergy(GS_lines)
                GS_Energy = output["Energy"][0] 
           # print (GS_Energy)
            if is_vip:
                IP_lines = Links_split[job_index['IP_line']]
                IP_Energy_SS = parse_log.Extract_SCFEnergy(IP_lines)
                # Normal ionization potential calculation
                output["vip"] = [Eh2eV*(IP_Energy_SS[0]-GS_Energy), IP_Energy_SS[1]]
            if is_vea:
                EA_lines = Links_split[job_index['EA_line']]
                EA_Energy_SS =  parse_log.Extract_SCFEnergy(EA_lines)
                # Normal electronic affinity calculation
                output["vea"] = [Eh2eV*(GS_Energy-EA_Energy_SS[0]), EA_Energy_SS[1]]
        
        if is_aip or is_aea:
            try:
                GS_Energy = output["Energy"][0]
            except KeyError:
                output["Energy"] = parse_log.Extract_SCFEnergy(GS_lines)
                GS_Energy = output["Energy"][0] 
            if is_aip:
                PC_lines = Links_split[job_index['PC_line']]
                VNP_lines = Links_split[job_index['VNP_line']]
                PC_Energy_SS = parse_log.Extract_SCFEnergy(PC_lines)
                VNP_Energy_SS = parse_log.Extract_SCFEnergy(VNP_lines)
                # For Check internal coordinate
                MaxDisplace = gaussian_run.Get_MolInterCoordinate.Extract_InterMol(PC_lines)
                output["relaxedIP_MaxDisplace"] = MaxDisplace
                output["aip"] = [Eh2eV*(PC_Energy_SS[0]-GS_Energy), Eh2eV*(PC_Energy_SS[0]-VNP_Energy_SS[0]), PC_Energy_SS[1], VNP_Energy_SS[1]]
            if is_aea:
                NC_lines = Links_split[job_index['NC_line']]
                VNN_lines = Links_split[job_index['VNN_line']]
                NC_Energy_SS = parse_log.Extract_SCFEnergy(NC_lines)
                VNN_Energy_SS = parse_log.Extract_SCFEnergy(VNN_lines)
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
            lines = Links_split[job_index['uv_line']] 
            _, _, _, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
                CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden, \
                mu_allowed, mu_forbidden, theta_allowed, theta_forbidden, g_allowed, g_forbidden = parse_log.Extract_ExcitedState(lines)
            output["uv"] = [WL_allowed, OS_allowed, CD_L_allowed, CD_OS_allowed]
            output["state_index"] = [State_allowed, State_forbidden]
            output["CD_mu"] = [mu_allowed, mu_forbidden]
            output["CD_theta"] = [theta_allowed, theta_forbidden]
            output["CD_g"] = [g_allowed, g_forbidden]
            if self.ref_uv_path != '':
                ref_uv = {}
                print(f"Read the reference spectrum...{self.ref_uv_path}")
                ref_uv["uv"] = gaussian_run.UV_similarity.read_data(self.ref_uv_path)
                S, D = gaussian_run.UV_similarity.smililarity_dissimilarity(ref_uv["uv"][0], ref_uv["uv"][1], output["uv"][0], output["uv"][1])
                output["Similality/Disdimilarity"] = [S, D]
        
        if is_fluor:
            lines = Links_split[job_index['relaxAEstate']]
            S_Found, S_Egrd, S_Eext, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
                CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden, \
                mu_allowed, mu_forbidden, theta_allowed, theta_forbidden, g_allowed, g_forbidden = parse_log.Extract_ExcitedState(lines)
            # For Check internal coordinate
            MaxDisplace = gaussian_run.Get_MolInterCoordinate.Extract_InterMol(lines)
            output["MinEtarget"] = S_Eext
            output["Min_MaxDisplace"] = MaxDisplace
            output["fluor"] = [WL_allowed, OS_allowed, CD_L_allowed, CD_OS_allowed]
            #mu, theta, g_cal = parse_log.extract_transtion_EMmoment(lines)
            output["CPL_mu"] = mu_allowed 
            output["CPL_theta"] = theta_allowed
            output["CPL_g"] = g_allowed
        
        if is_tadf:
            lines = Links_split[job_index['relaxFEstate']] 
            T_Found, _, T_Eext, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
                CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden, \
                mu_allowed, mu_forbidden, theta_allowed, theta_forbidden, g_allowed, g_forbidden = parse_log.Extract_ExcitedState(lines)
            # For Check internal coordinate
            MaxDisplace = gaussian_run.Get_MolInterCoordinate.Extract_InterMol(lines)
            output["T_Min"] = T_Eext
            output["T_Min_MaxDisplace"] = MaxDisplace
            output["T_Phos"] = [WL_forbidden, OS_forbidden, CD_L_forbidden, CD_OS_forbidden]
            output["CPL_T_mu"] = mu_forbidden 
            output["CPL_T_theta"] = theta_forbidden
            output["CPL_T_g"] = g_forbidden
            #TADF_Eng = 0.0
            #if S_Found and T_Found:
            #   TADF_Eng = S_Eext - T_Eext
            #output["Delta(S-T)"] = TADF_Eng

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

    def make_input(
        self, JobName, TotalCharge, SpinMulti, 
        scf='open', run_type=None, optoption='', Newinput=False, 
        Mol_atom=[], X=[], Y=[], Z=[], geom_spec=False,
        TDDFT=False, TDstate=None, target=1,
        readchk=None, oldchk=None, newchk=None, solvent='0'):

        #Section for system control 
        line_system = ''
        if self.mem !='' or self.nproc > 1:
            if self.mem != '':
                line_mem = f'%mem={self.mem}\n'
                line_system += line_mem
            if self.nproc > 1 :
                line_proc = f'%nproc={self.nproc}\n'
                line_system += line_proc

        #Section for checkpoint files
        if oldchk == None and newchk == None:
            line_chk = f'%chk={JobName}\n'
        elif oldchk == None and newchk != None:
            line_chk = f'%chk={newchk}\n'
            line_chk += f'%Oldchk={JobName}\n'
        elif oldchk != None and newchk == None:
            line_chk = f'%chk={oldchk}\n'
        else:
            line_chk = f'%chk={newchk}\n'
            line_chk += f'%Oldchk={oldchk}\n'

        #Section for DFT functional parameter
        line_iop_functional = ''
        if self.para_functional == []:
            pass
        else:
            line_iop_functional += gaussian_run.Gen_IOPline_functional.functional_para(self.functional, self.para_functional)

        #Section for method
        if scf == 'open' or SpinMulti != 1:
            line_method   = f'#u{self.functional}/{self.basis} test {line_iop_functional}\n'
        elif scf == 'close' and SpinMulti == 1:
            line_method = f'#r{self.functional}/{self.basis} test {line_iop_functional}\n' 
        else:
            line_method = f'#{self.functional}/{self.basis} test {line_iop_functional}\n' 

        if TDDFT:
            NState = target+4 if (run_type == 'opt' or run_type == 'freq') else 20
            line_method += f'TD(Nstate={NState}, root={target})\n'if TDstate == None else f'TD(Nstate={NState}, {TDstate}, root={target})\n'
        
        #Section for solvent
        SCRF = ''
        SCRF_read = ''
        if solvent != '0':
            self.solvent = solvent
            SCRF, SCRF_read = self.MakeSolventLine()

        #Geometry information
        if geom_spec: 
            line_GeomConstrain = 'Geom=ModRedundant\n'
            line_readMOConstrain = 'Geom=ModRedundant Guess=Read\n'
            line_readGeomConstrain = 'Geom=(Checkpoint,ModRedundant)\n'
            line_readAllMOGeomConstrain = 'Geom=(AllCheckpoint,ModRedundant) Guess=Read\n'
            line_readMOGeomConstrain = 'Geom=(Checkpoint,ModRedundant) Guess=Read\n'
            line_GeomConstrainSpec = gaussian_run.ConstrainInfo_line.get_IntCoorddata(self.geom_spec)
        else:
            line_readAllMOGeom = 'Geom=AllCheck Guess=Read\n'
            line_readMOGeom = 'Geom=CheckPoint Guess=Read\n'
            line_readMO = 'Guess=Read\n'
            line_readGeom = 'Geom=Checkpoint\n'

        #Comment line
        line_comment = f'\n {JobName} {run_type} \n \n'

        GauInputName = JobName+'.com' 

        # A new file will be made if the specified file name does not exist.
        # If the same name input file exist, that file will be open and new lines will be added.
        if Newinput == True:
            Input_mode = 'w'
        else:
            if not os.path.isfile(GauInputName):
                Input_mode = 'w'
            else:
                Input_mode = 'a'

        with open(GauInputName, Input_mode) as ofile:
            if Input_mode == 'a':
                input_s = '--Link1--\n'
            else:
                input_s = ''

            input_s += line_system
            input_s += line_chk
            input_s += line_method

            #For non SCF run_type
            if run_type == 'symm' or run_type == 'volume':
                input_s += 'Guess=Only'
                if run_type == 'volume':
                    input_s += ' volume\n'
                if run_type == 'symm':
                    input_s += ' Symmetry=loose\n'
            
            if run_type == 'opt':
                if TDDFT:
                    input_s += 'Opt=(Maxcycles=60)\n'
                else:
                    if optoption == '':
                        input_s += 'Opt=(MaxCycles=100)\n'
                    else:
                        input_s += f'Opt=({optoption}, MaxCycles=100)\n'

            if run_type == 'polar':
                input_s += 'polar CPHF=Static\n'
            
            if run_type == 'freq':
                input_s += 'Freq=(Raman)\n' 
            
            if run_type == 'nmr':
                input_s += 'NMR\n'

            if solvent != '0':
                input_s += SCRF
            
            #Get geometry and guess information
            if  len(Mol_atom) == 0 and readchk != None:
                if geom_spec:
                    if readchk == 'all':
                        input_s += line_readAllMOGeomConstrain
                    elif readchk == 'geomguess':
                        input_s += line_readMOGeomConstrain  
                    elif readchk == 'geom':
                        input_s += line_readGeomConstrain
                else:
                    if readchk == 'all':
                        input_s += line_readAllMOGeom
                    elif readchk == 'geomguess':
                        input_s += line_readMOGeom
                    elif readchk == 'geom':
                        input_s += line_readGeom
            elif len(Mol_atom) != 0:
                if geom_spec:
                    if readchk == None:
                        input_s += line_GeomConstrain
                    if readchk == 'guess':
                        input_s += line_readMOConstrain
                else:
                    if readchk == 'guess':
                        input_s +=  line_readMO
            else:
                print ('No information about a molecules!')
                exit

            if readchk != 'all':
                input_s += line_comment
            else:
                input_s += '\n' 

            if len(Mol_atom) != 0 or readchk != 'all':
                input_s += f'{TotalCharge: 5d} {SpinMulti: 5d} \n'
                for j in range(len(Mol_atom)):
                    input_s += f'{Mol_atom[j]:4s} {X[j]: 10.5f}  {Y[j]: 10.5f} {Z[j]: 10.5f}\n'
                input_s += '\n'
            if solvent != '0':
                input_s += SCRF_read
            if geom_spec:
                input_s += line_GeomConstrainSpec
                input_s += '\n'

            ofile.write(input_s)

    def chain_job(self, JobName, scf_need, job_dict, TotalCharge, SpinMulti, targetstate,
                        ReadFrom, element=[], atomX=[], atomY=[], atomZ=[], optoption='', TDstate_info=[]):

        option_dict = copy.deepcopy(job_dict)

        print(f'Charge: {TotalCharge}, Spin Multiplicity: {SpinMulti}')

        # Setting molecular charge and multiplicity
        Is_ChargeSpec = False
        Is_SpinMultiSpec = False
        if np.isnan(self.SpecTotalCharge) !=  True:
            Is_ChargeSpec = True
            TotalCharge = self.SpecTotalCharge
        if np.isnan(self.SpecSpinMulti) != True:
            Is_SpinMultiSpec = True
            SpinMulti = self.SpecSpinMulti

        # Setting  geometric constrain
        if self.geom_spec != {}: 
            Is_geom_spec = True
        else:
            Is_geom_spec = False

        #For setting TD-DFT
        scf_tag='open'
        td_tag = False
        TDstate_tag = None 
        if targetstate != 0:
            td_tag = True
            if SpinMulti == 1:
                TDstate_tag = 'Singlet'
                scf_tag='close'
            elif SpinMulti == 3:
                SpinMulti = 1
                TDstate_tag = 'Triplet'
                scf_tag='close'
            else:
                if TDstate_info != []:
                    targetstate =  TDstate_info[int(targetstate)-1]
                    TDstate_tag = None 
                else:
                    print('No information about the target state, we set the first excited state the target.')
                    targetstate = 1

        if 'fluor' in option_dict:
               del option_dict['fluor'] 
               option_dict['opt'] = True
        if 'tadf' in option_dict:
               del option_dict['tadf'] 
               option_dict['opt'] = True

        if scf_need == False:
        #Only no-scf jobs
            if 'symm' in option_dict:
                if ReadFrom == 'chk':
                    self.make_input(JobName, TotalCharge, SpinMulti, 
                                    run_type='symm', Newinput=True, 
                                    geom_spec=Is_geom_spec, readchk='all') 
                if ReadFrom == 'sdf' or ReadFrom == 'xyz':
                    self.make_input(JobName, TotalCharge, SpinMulti, 
                                    run_type='symm', Newinput=True, 
                                    Mol_atom=element, X=atomX, Y=atomY, Z=atomZ, 
                                    geom_spec=Is_geom_spec) 
            if 'volume' in option_dict:
                if ReadFrom == 'chk':
                    self.make_input(JobName, TotalCharge, SpinMulti, 
                                    run_type='volume', Newinput=False, 
                                    geom_spec=Is_geom_spec, readchk='all') 
                if ReadFrom == 'sdf' or ReadFrom == 'xyz':
                    self.make_input(JobName, TotalCharge, SpinMulti, 
                                    run_type='volume', Newinput=False, 
                                    Mol_atom=element, X=atomX, Y=atomY, Z=atomZ, 
                                    geom_spec=Is_geom_spec) 
        else:
        #The first scf job
            run_opt = 'opt' if 'opt' in option_dict else ''
            if ReadFrom == 'chk' and element == []:
                if Is_ChargeSpec == False and Is_SpinMultiSpec == False:
                    if optoption == '':
                        self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag,
                                        run_type=run_opt, TDDFT=td_tag, TDstate=TDstate_tag, Newinput=True, 
                                        geom_spec=Is_geom_spec, readchk='all', solvent=self.solvent) 
                    else:
                        self.make_input(JobName, TotalCharge, SpinMulti, 
                                        run_type=run_opt, TDDFT=td_tag, TDstate=TDstate_tag, optoption=optoption, Newinput=True, 
                                        geom_spec=Is_geom_spec, readchk='all', solvent=self.solvent) 
                else:
                    if optoption == '':
                        self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag, 
                                        run_type=run_opt, TDDFT=td_tag, TDstate=TDstate_tag, Newinput=True, 
                                        geom_spec=Is_geom_spec, readchk='geomguess', solvent=self.solvent) 
                    else:
                        self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag, 
                                        run_type=run_opt, TDDFT=td_tag, TDstate=TDstate_tag, optoption=optoption, Newinput=True, 
                                        geom_spec=Is_geom_spec, readchk='geomguess', solvent=self.solvent) 
            elif ReadFrom == 'chk' and element != []:
                if optoption == '':
                    self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag, 
                                    run_type=run_opt, TDDFT=td_tag, TDstate=TDstate_tag, Newinput=True, 
                                    Mol_atom=element, X=atomX, Y=atomY, Z=atomZ, 
                                    geom_spec=Is_geom_spec, readchk='guess', solvent=self.solvent) 
                else:
                    self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag, 
                                    run_type=run_opt, TDDFT=td_tag, TDstate=TDstate_tag, optoption=optoption, Newinput=True, 
                                    Mol_atom=element, X=atomX, Y=atomY, Z=atomZ, 
                                    geom_spec=Is_geom_spec, readchk='guess', solvent=self.solvent) 
            elif ReadFrom == 'sdf' or ReadFrom == 'xyz':
                if optoption == '':
                    self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag, 
                                    run_type=run_opt, TDDFT=td_tag, TDstate=TDstate_tag, Newinput=True, 
                                    Mol_atom=element, X=atomX, Y=atomY, Z=atomZ, 
                                    geom_spec=Is_geom_spec, solvent=self.solvent) 
                else:
                    self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag, 
                                    run_type=run_opt, TDDFT=td_tag, TDstate=TDstate_tag, optoption=optpotion, Newinput=True, 
                                    Mol_atom=element, X=atomX, Y=atomY, Z=atomZ, 
                                    geom_spec=Is_geom_spec, solvent=self.solvent) 

        #The post job after getting wavefunction
            if 'symm' in option_dict:
                self.make_input(JobName, TotalCharge, SpinMulti, 
                                run_type='symm', Newinput=False, 
                                readchk='all', solvent=self.solvent) 

            if 'volume' in option_dict:
                self.make_input(JobName, TotalCharge, SpinMulti, 
                                run_type='volume', Newinput=False, 
                                readchk='all', solvent=self.solvent) 

            if 'freq' in option_dict: # freq == 1
                self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag, 
                                run_type='freq', TDDFT=td_tag, TDstate=TDstate_tag, Newinput=False, 
                                readchk='all', solvent=self.solvent) 

            if 'polar' in option_dict: # polar == 1
                self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag, 
                                run_type='polar', TDDFT=td_tag, TDstate=TDstate_tag, Newinput=False, 
                                readchk='all', solvent=self.solvent) 

            if 'nmr' in option_dict: # nmr == 1
                self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag, 
                                run_type='nmr', TDDFT=td_tag, TDstate=TDstate_tag, Newinput=False, 
                                readchk='all', solvent=self.solvent) 

            if 'vip' in option_dict: # vip == 1
                IP_chk = JobName+'_IP'

                if SpinMulti == 1:
                    IP_TotalCharge = TotalCharge + 1
                    IP_SpinMulti = SpinMulti + 1
                else:
                    IP_TotalCharge = TotalCharge + 1
                    IP_SpinMulti = SpinMulti - 1

                self.make_input(JobName, IP_TotalCharge, IP_SpinMulti, 
                                run_type='', Newinput=False, 
                                readchk='geomguess', newchk=IP_chk, solvent=self.solvent) 

            if 'vea' in option_dict: # vea == 1
                EA_chk = JobName+'_EA'

                if SpinMulti == 1:
                    EA_TotalCharge = TotalCharge - 1
                    EA_SpinMulti = SpinMulti + 1
                else:
                    EA_TotalCharge = TotalCharge - 1
                    EA_SpinMulti = SpinMulti - 1

                self.make_input(JobName, EA_TotalCharge, EA_SpinMulti, 
                                run_type='', Newinput=False, 
                                readchk='geomguess', newchk=EA_chk, solvent=self.solvent) 

            if 'aip' in option_dict: # aip == 1
                AIP_chk = JobName+'_AIP'
                #AIP energy
                self.make_input(JobName, IP_TotalCharge, IP_SpinMulti, 
                                run_type='opt', 
                                readchk='geomguess', oldchk=IP_chk, newchk=AIP_chk, solvent=self.solvent) 

                OE_AIP_chk = JobName+'_AIP_OE'
                #Neutrization energy
                self.make_input(JobName, TotalCharge, SpinMulti, 
                                run_type='', 
                                readchk='geomguess', oldchk=AIP_chk, newchk=OE_AIP_chk, solvent=self.solvent) 

            if 'aea' in  option_dict: #aea == 1
                AEA_chk = JobName+'_AEA'
                #AEA energy
                self.make_input(JobName, EA_TotalCharge, EA_SpinMulti, 
                                run_type='opt', 
                                readchk='geomguess', oldchk=EA_chk, newchk=AEA_chk, solvent=self.solvent) 

                OE_AEA_chk = JobName+'_AEA_OE'
                #Neutrization energy
                self.make_input(JobName, TotalCharge, SpinMulti, 
                                run_type='', readchk='geomguess', oldchk=AEA_chk, newchk=OE_AEA_chk, solvent=self.solvent) 

            if 'uv' in option_dict: #uv == 1 of fluor ==1 or tadf ==1
                self.make_input(JobName, TotalCharge, SpinMulti, 
                                run_type='', TDDFT=True, 
                                readchk='all', solvent=self.solvent) 

           # if 'fluor' in option_dict:
           #     TDOpt_chk = f'{JobName}_ExOptAState{targetstate}' if targetstate != 1 else f'{JobName}_ExOptStateSinglet'
           #     self.make_input(JobName, TotalCharge, SpinMulti, 
           #                     scf='close', run_type='opt', TDDFT=True, TDstate='Singlet', target=targetstate, 
           #                     readchk='all', newchk=TDOpt_chk, solvent=self.solvent) 

           # if 'tadf' in option_dict:
           #     TTDOpt_chk = f'{JobName}_ExOptFState{targetstate}' if targetstate != 1 else f'{JobName}_ExOptStateTriplet'
           #     self.make_input(JobName, TotalCharge, SpinMulti, 
           #                     scf='close', run_type='opt', TDDFT=True, TDstate='Triplet', target=targetstate, 
           #                     readchk='all', newchk=TTDOpt_chk, solvent=self.solvent) 

    def SpinMulti_scan(self, JobName,  targetstate, ReadFrom, GivenSpinMulti, TotalCharge, atm, X, Y, Z, TDstate_info):
        print('Try to optimize the spin state of the ground state!')
        
        SpinMulti_list = []
        SpinDiff_list = []
        Energy_list = []
        
        if GivenSpinMulti%2 != 0:
            print("The number of electron is even.")
            InitSpinMulti = 1
        else:
            print("The number of electron is odd.")
            InitSpinMulti = 2
        
        for i in range(InitSpinMulti, InitSpinMulti+2*3, 2):
        
            JobNameState = JobName + f'_State{targetstate}_{i}_{TotalCharge}'
        
            option_dict_spincheck = {'energy': True}
            scf_need=True
            self.chain_job(JobNameState, scf_need, option_dict_spincheck, TotalCharge, i, 
                         targetstate, ReadFrom, 
                         element=atm, atomX=X, atomY=Y, atomZ=Z, optoption='', TDstate_info=TDstate_info)
            job_state = "normal"
            job_state = gaussian_run.Exe_Gaussian.exe_Gaussian(JobNameState, self.timexe)
        
            try:
                output_prop = self.Extract_values(JobNameState, option_dict_spincheck, Bondpair1=[], Bondpair2=[])
            except Exception as e:
                job_state = "error"
                print(e)
                pass

            if job_state == "normal":
                SpinMulti_list.append(i)
                Energy_list.append(output_prop["Energy"][0])
                SpinDiff_list.append(output_prop["Energy"][1])
            else:
                pass

        StdSpinDiff = np.std(SpinDiff_list)

        print(SpinDiff_list)
        print(StdSpinDiff)
        print(Energy_list)
        print(np.std(Energy_list))
        Index_LowSpinDiff = SpinDiff_list.index(min(SpinDiff_list))
        Index_LowEnergy = Energy_list.index(min(Energy_list))
        
        if StdSpinDiff < 0.1:
            SpinMulti = SpinMulti_list[Index_LowEnergy]
        else:
            SpinMulti = SpinMulti_list[Index_LowSpinDiff]

        return SpinMulti

    def run_gaussian(self):
        infilename = self.in_file
        option_line = self.value    
        options = option_line.split()
        job_eachState = []
        job_eachState.append({})
        option_dict = {}
        TargetStates = [0]
        TargetSpinMulti = []
        TargetTotalCharge = []
        PreGauInput = infilename.split('.')
        JobName = PreGauInput[0]
        #GauInputName = PreGauInput[0]+'.com'    
        # File type of input?
        ReadFrom = None
        # Initialization of molecular coordinate
        atm, X, Y, Z = [], [], [], []
        Bondpair1 = []
        Bondpair2 = []
        # Classify input
        if PreGauInput[1] == "sdf":
            ReadFrom = 'sdf' 
            atm, X, Y, Z, TotalCharge, SpinMulti, Bondpair1, Bondpair2 = read_mol_file.read_sdf(infilename)
        elif PreGauInput[1] == "xyz":
            ReadFrom = 'xyz' 
            atm, X, Y, Z, TotalCharge, SpinMulti = read_mol_file.read_xyz(infilename)
        elif PreGauInput[1] == "chk":
            ReadFrom = 'chk' 
        elif PreGauInput[1] == "fchk":
            TotalCharge, SpinMulti = gaussian_run.fchk2chk.Get_fchk(PreGauInput[0])
            ReadFrom = 'chk' 
        else:
            print("Invalid input file")

        # Specify the electronic structure (spin multiplicity and charge) of the target state
        if np.isnan(self.SpecTotalCharge) !=  True:
            Is_ChargeSpec = True
            TotalCharge = self.SpecTotalCharge
        if np.isnan(self.SpecSpinMulti) != True:
            Is_SpinMultiSpec = True
            SpinMulti = self.SpecSpinMulti

        TargetSpinMulti.append(SpinMulti)
        TargetTotalCharge.append(TotalCharge)

        # Setting options
        for i in range(len(options)):
            option = options[i]
            if option.lower() == 'opt':
                option_dict['opt'] = True
                job_eachState[0]['opt'] = True
                if '=' in option:
                    in_target = option.split('=')
                    optoption = in_target[-1]
            elif 'optspin' in option.lower():
                option_dict['optspin'] = True
                option_dict['energy'] = True
                job_eachState[0]['optspin'] = True
                job_eachState[0]['energy'] = True
            elif option.lower() == 'freq':
                option_dict['freq'] = True
                job_eachState[0]['freq'] = True
            elif option.lower() == 'polar':
                option_dict['polar'] =  True
                job_eachState[0]['polar'] = True
            elif option.lower() == 'nmr':
                option_dict['nmr'] = True
                job_eachState[0]['nmr'] = True
            elif 'uv' in option.lower():
                option_dict['uv'] = True
                job_eachState[0]['uv'] = True
                if '=' in option:
                    ref_spectrum = option.split('=')
                    self.ref_uv_path = ref_spectrum[-1]
            elif option.lower() == 'energy':
                option_dict['energy'] = True
                job_eachState[0]['energy'] = True
            elif option.lower() == 'homolumo':
                option_dict['homolumo'] = True
                job_eachState[0]['homolumo'] = True
            elif option.lower() == 'dipole':
                option_dict['dipole'] = True
                job_eachState[0]['dipole'] = True
            elif option.lower() == 'deen':
                option_dict['deen'] = True
                job_eachState[0]['deen'] = True
            elif option.lower() == 'stable2o2':
                option_dict['stable2o2'] = True
                job_eachState[0]['stable2o2'] = True
            elif 'fluor' in option.lower():
                option_dict['uv'] = True
                option_dict['fluor'] = True
                job_eachState.append({})
                job_eachState[0]['uv'] = True
                job_eachState[len(job_eachState)-1]['fluor'] = True
                #setting the target state and its spin
                if '=' in option:
                    in_target = option.split('=')
                    TargetStates.append(int(in_target[-1]))
                    TargetSpinMulti.append(SpinMulti)
                    TargetTotalCharge.append(TotalCharge)
                else:
                    TargetStates.append(1)
                    TargetSpinMulti.append(SpinMulti)
                    TargetTotalCharge.append(TotalCharge)
            elif option.lower() == 'tadf':
                option_dict['uv'] = True
                job_eachState[0]['uv'] = True
                job_eachState.append({})
                if 'fluor' in option_dict:
                    option_dict['tadf'] = True
                    job_eachState[len(job_eachState)-1]['tadf'] = True
                    #setting the target state and its spin
                    if '=' in option:
                        in_target = option.split('=')
                        TargetStates.append(int(in_target[-1]))
                        TargetSpinMulti.append(SpinMulti+2)
                        TargetTotalCharge.append(TotalCharge)
                    else:
                        TargetStates.append(1)
                        TargetSpinMulti.append(SpinMulti+2)
                        TargetTotalCharge.append(TotalCharge)
                else:
                    option_dict['fluor'] = True
                    option_dict['tadf'] = True
                    job_eachState.append({})
                    job_eachState[len(job_eachState)-2]['fluor'] = True
                    job_eachState[len(job_eachState)-1]['tadf'] = True
                    #setting the target state and its spin
                    if '=' in option:
                        in_target = option.split('=')
                        TargetStates.append(int(in_target[-1]))
                        TargetSpinMulti.append(SpinMulti)
                        TargetTotalCharge.append(TotalCharge)
                        TargetStates.append(int(in_target[-1]))
                        TargetSpinMulti.append(SpinMulti+2)
                        TargetTotalCharge.append(TotalCharge)
                    else:
                        TargetStates.append(1)
                        TargetStates.append(1)
                        TargetSpinMulti.append(SpinMulti)
                        TargetSpinMulti.append(SpinMulti+2)
                        TargetTotalCharge.append(TotalCharge)
                        TargetTotalCharge.append(TotalCharge)
            elif option.lower() == 'vip':
                job_eachState[0]['vip'] = True
                option_dict['vip'] = True
            elif option.lower() == 'vea':
                job_eachState[0]['vea'] = True
                option_dict['vea'] = True
            elif option.lower() == 'aip':
                job_eachState[0]['vip'] = True
                job_eachState[0]['aip'] = True
                option_dict['vip'] = True
                option_dict['aip'] = True
            elif option.lower() == 'aea':
                job_eachState[0]['vea'] = True
                job_eachState[0]['aea'] = True
                option_dict['vea'] = True
                option_dict['aea'] = True
            elif option.lower() == 'cden':
                job_eachState[0]['cden'] = True
                option_dict['cden'] = True
            elif option.lower() == 'pka':
                TargetStates.append(0)
                job_eachState[0]['energy'] = True
                job_eachState[0]['cden'] = True
                job_eachState.append({})
                job_eachState[len(job_eachState)-1]['energy'] = True
                job_eachState[len(job_eachState)-1]['opt'] = True
                job_eachState[len(job_eachState)-1]['pka'] = True
                TargetSpinMulti.append(SpinMulti)
                TargetTotalCharge.append(TotalCharge-1)
                option_dict['energy'] = True
                option_dict['cden'] = True
                option_dict['pka'] = True
            elif option.lower() == 'satkoopmans':
                job_eachState[0]['homolumo'] = True
                job_eachState[0]['vip'] = True
                job_eachState[0]['vea'] = True
                job_eachState[0]['satkoopmans'] = True
                option_dict['homolumo'] = True
                option_dict['vip'] = True
                option_dict['vea'] = True
                option_dict['satkoopmans'] = True
            elif option.lower() == 'symm':
                job_eachState[0]['symm'] = True
                option_dict['symm'] = True
            elif option.lower() == 'volume':
                job_eachState[0]['volume'] = True
                option_dict['volume'] = True
            elif 'stable' in option.lower():
                option_dict['stable'] = True
            else:
                print('Invalid option: ', option)

        # Make work directory and move to the directory
        if os.path.isdir(JobName):
            shutil.rmtree(JobName)
        os.mkdir(JobName)
        if ReadFrom == 'chk':
            inchkfile = JobName + '.chk'
            shutil.move(inchkfile, JobName) 
        os.chdir(JobName)

        # Initialization of otuput dictionary.
        output_dic = []

        TDstate_info = []
        #print(job_eachState)
        #print(TargetTotalCharge)
        #print(TargetSpinMulti)

        job_state = ""

        for i in range(len(TargetStates)):

            print (f'~~~~~~~~~~~~~Computing the state/species {i}~~~~~~~~~~~~~~~~~~~')

            output_dic.append({})

            job_thisstate = job_eachState[i]
            targetstate = TargetStates[i]
            SpinMulti = TargetSpinMulti[i]
            TotalCharge = TargetTotalCharge[i]
            if 'pka' in job_thisstate:
                JobNameState = JobName + f'_DeH'
            else:
            #For scanning the spin state of the target 
                if 'optspin' in job_thisstate:
                    SpinMulti  = self.SpinMulti_scan(targetstate=targetstate, JobName=JobName, ReadFrom=ReadFrom, GivenSpinMulti=SpinMulti, TotalCharge=TotalCharge, 
                                atm=atm, X=X, Y=Y, Z=Z, TDstate_info=TDstate_info)
                    JobNameState = JobName + f'_State{targetstate}_{SpinMulti}_{TotalCharge}'
                    ReadFrom == 'chk'
                    output_dic[i]['spinmulti'] = SpinMulti
                else:
                    JobNameState = JobName + f'_State{targetstate}_{SpinMulti}_{TotalCharge}'

            Jobchk = JobNameState  + '.chk' 
            AttmptStable = True

            if ReadFrom == 'chk' and i == 0:
                inchkfile = JobName + '.chk'
                for f in glob.glob('./*.chk'):
                    File_name = f.replace('./','')
                    if File_name == inchkfile:
                        shutil.move(File_name, Jobchk) 
                    else:
                        print ('Suitable chk file was not found!')
            
            elif i > 0 and ('pka' not in job_thisstate):
                for f in glob.glob('./*.fchk'):
                    print(f)
                    File_name = f.replace('./','')
                    if File_name == JobName + f'_State0_{SpinMulti}_{TotalCharge}.fchk' and ('fluor' in job_thisstate):
                        TDstate_info = output_dic[0]["state_index"][0] 
                        #print(TDstate_info)
                        print ('The ground state fchk file was found!')
                        GS_TotalCharge, GS_SpinMulti = gaussian_run.fchk2chk.Get_fchk(File_name.replace('.fchk',''))
                        shutil.copy(File_name.replace('.fchk','.chk'), Jobchk) 
                        ReadFrom = 'chk'
                        # Initialization of molecular coordinate
                        atm, X, Y, Z = [], [], [], []
                    #if File_name == JobName + f'_State1_{SpinMulti-2}_{TotalCharge}.fchk' and  ('tadf' in job_thisstate):
                    if File_name == JobName + f'_State0_{SpinMulti-2}_{TotalCharge}.fchk' and  ('tadf' in job_thisstate):
                        TDstate_info = output_dic[0]["state_index"][1] 
                        #print(TDstate_info)
                        print ('The ground state fchk file was found!')
                        GS_TotalCharge, GS_SpinMulti = gaussian_run.fchk2chk.Get_fchk(File_name.replace('.fchk',''))
                        shutil.copy(File_name.replace('.fchk','.chk'), Jobchk) 
                        ReadFrom = 'chk'
                        # Initialization of molecular coordinate
                        atm, X, Y, Z = [], [], [], []
           #         else:
           #             print ('Suitable chk file was not found!')

            elif 'pka' in job_thisstate:
                E_pH = output_dic[0]["Energy"][0]
                Atom = output_dic[0]["cden"][0]
                MullCharge = output_dic[0]["cden"][1]
                Index_MaxProtic = gaussian_run.Get_ChargeSpin.find_MaxProtic(Atom, MullCharge)
                print(Index_MaxProtic)
                GS_fchk = JobName+f'_State0_{SpinMulti}_{TotalCharge+1}.fchk'
                with open(GS_fchk,'r') as ifile:
                    GS_fchk_lines = ifile.readlines()
                _, _, Mol_atom, Mol_X, Mol_Y, Mol_Z = gaussian_run.Get_MolCoordinate_fchk.Extract_MolCoord(GS_fchk_lines)
                atm = np.delete(Mol_atom,Index_MaxProtic) 
                X = np.delete(Mol_X,Index_MaxProtic)
                Y = np.delete(Mol_Y,Index_MaxProtic)
                Z = np.delete(Mol_Z,Index_MaxProtic)
                ReadFrom = 'xyz'

            if 'stable' in option_dict:

                if self.geom_spec != {}: 
                    for k in self.geom_spec.keys():
                        if self.geom_spec[k][1] == 'F':
                            print('Specified frozen molecular freedoms are found.')
                            print('We will not try to find a stable structure.')
                            #setting the interation limit, relese from the cycle to find stable structure.
                            Stable_Iteration = 9
                        else:
                            Stable_Iteration = 0
                else:
                    Stable_Iteration = 0

                AttmptStable = False
                option_dict_optcheck = {'opt': True, 'freq': True, 'energy': True}
                scf_need=True
                self.chain_job(JobNameState, scf_need, option_dict_optcheck, TotalCharge, SpinMulti, 
                             targetstate, ReadFrom, 
                             element=atm, atomX=X, atomY=Y, atomZ=Z, optoption='', TDstate_info=TDstate_info)
                job_state = gaussian_run.Exe_Gaussian.exe_Gaussian(JobNameState, self.timexe)
            
                while AttmptStable == False:

                    try:
                        output_prop = self.Extract_values(JobNameState, option_dict_optcheck, Bondpair1, Bondpair2)
                    except Exception as e:
                        job_state = "error"
                        print(e)
                        pass
            
                    #print(output_prop)
                    Freq_CK = np.array(output_prop['freq'])
                    Freq_min = Freq_CK.min()
                    Freq_min_index = np.argmin(Freq_CK)
                    
                    #For the next job, change  RaedFrom to 'chk' 
                    ReadFrom = 'chk'
            
                    if Freq_min > 0:
                        print (f'Mnimum frequency is positive...{Freq_min} cm^-1')
                        #remove the stable option
                        output_dic[i][f'stable_{i}'] = True
                        AttmptStable = True
                        Stable_Iteration += 1
            
                        #Initializing molecular coordinate since the next job will use chk file
                        atm, X, Y, Z = [], [], [], [] 
            
                        #Checking the remain task (option)
                        all_options = set(job_thisstate.keys())
                        remain_keys = all_options - set(option_dict_optcheck.keys())
            
                        if len(remain_keys) >= 1:
                            output_dic[i].update(output_prop)
                            print(remain_keys)
                            if 'freq' in job_thisstate:
                                del job_thisstate['freq']
                            if 'opt' in job_thisstate:
                                del job_thisstate['opt']
                            #for the next computation
                        else:
                            output_dic[i].update(output_prop)

                        print(f'Remaining job....{job_thisstate}')
            
                    else:
                        print (f'Negative frequency is detected...{Freq_min} cm^-1')
                
                        Stable_Iteration += 1
                        if  Stable_Iteration == 10: # quit after over 10 iterations
                            output_dic[i][f'stable_{i}'] = False
                            output_dic[i].update(output_prop)
                            AttmptStable = True

                            #Checking the remain task (option)
                            all_options = set(job_thisstate.keys())
                            remain_keys = all_options - set(option_dict_optcheck.keys())
            
                            if len(remain_keys) >= 1:
                                output_dic[i].update(output_prop)
                                print(remain_keys)
                                if 'freq' in job_thisstate:
                                    del job_thisstate['freq']
                                if 'opt' in job_thisstate:
                                    del job_thisstate['opt']
                            #for the next computation
                            else:
                                output_dic[i].update(output_prop)
                            print(f'Remaining job....{job_thisstate}')
                            break

                        #Get molecular coordinate
                        check_logfile = JobNameState + '.log'
                        with open(check_logfile,'r') as ifile:
                            CKF_lines = ifile.readlines()
                        Mol_atom, Mol_X, Mol_Y, Mol_Z = gaussian_run.Get_MolCoordinate.Extract_Coordinate(CKF_lines)
            
                        #try to erase the imaginary frequency mode...
                        ImFreq_modeX = output_prop['freqmode'][0][Freq_min_index]
                        ImFreq_modeY = output_prop['freqmode'][1][Freq_min_index]
                        ImFreq_modeZ = output_prop['freqmode'][2][Freq_min_index]
            
                        fX = Mol_X + 0.05 * ImFreq_modeX
                        fY = Mol_Y + 0.05 * ImFreq_modeY
                        fZ = Mol_Z + 0.05 * ImFreq_modeZ
            
                        rX = Mol_X - 0.05 * ImFreq_modeX
                        rY = Mol_Y - 0.05 * ImFreq_modeY
                        rZ = Mol_Z - 0.05 * ImFreq_modeZ
            
                        JobName_ChStableF = JobNameState+'_ChStableF'+f'_State{targetstate}'
                        JobName_ChStableR = JobNameState+'_ChStableR'+f'_State{targetstate}'
                        check_logfileF = JobName_ChStableF+'.log'
                        check_logfileR = JobName_ChStableR+'.log'
                        JobName_ChStableF_chk = JobName_ChStableF + '.chk'
                        JobName_ChStableR_chk = JobName_ChStableR + '.chk'
            
                        shutil.copy(Jobchk, JobName_ChStableF_chk)
                        shutil.copy(Jobchk, JobName_ChStableR_chk)
            
                        #explore the direction of forward mode 
                        self.chain_job(JobName_ChStableF, scf_need, option_dict_optcheck, TotalCharge, SpinMulti, 
                                            targetstate, ReadFrom, 
                                            element=atm, atomX=fX, atomY=fY, atomZ=fZ, optoption='ReadFC', TDstate_info=TDstate_info)
                        job_state = gaussian_run.Exe_Gaussian.exe_Gaussian(JobName_ChStableF, self.timexe)
            
                        try:
                            output_prop_f = self.Extract_values(JobName_ChStableF, option_dict_optcheck, Bondpair1, Bondpair2)
                        except Exception as e:
                            job_state = "error"
                            print(e)
                            output_prop_f['Energy'] = [0, 0]
                            pass
            
                        #explore the direction of reverse mode 
                        self.chain_job(JobName_ChStableR, scf_need, option_dict_optcheck, TotalCharge, SpinMulti, 
                                            targetstate, ReadFrom, 
                                            element=atm, atomX=rX, atomY=rY, atomZ=rZ, optoption='ReadFC', TDstate_info=TDstate_info)
                        job_state = gaussian_run.Exe_Gaussian.exe_Gaussian(JobName_ChStableR, self.timexe)
            
                        try:
                            output_prop_r = self.Extract_values(JobName_ChStableR, option_dict_optcheck, Bondpair1, Bondpair2)
                        except Exception as e:
                            job_state = "error"
                            print(e)
                            output_prop_r['Energy'] = [0, 0]
                            pass
            
                        if output_prop_f["Energy"][0] <=  output_prop_r["Energy"][0]:
                            output_prop = output_prop_f.copy()
                            shutil.copy(JobName_ChStableF_chk, Jobchk)
                            shutil.copy(check_logfileF, check_logfile)
                            output_prop_f = {}
                        else:
                            output_prop = output_prop_r.copy()
                            shutil.copy(JobName_ChStableR_chk, Jobchk)
                            shutil.copy(check_logfileR, check_logfile)
                            output_prop_r = {}

            if ('stable' not in option_dict) or AttmptStable:
                check_logfile = JobName+'.log'
                scf_need = self.check_scf_need(job_thisstate)
                self.chain_job(JobNameState, scf_need, job_thisstate, TotalCharge, SpinMulti, 
                                    targetstate, ReadFrom, 
                                    element=atm, atomX=X, atomY=Y, atomZ=Z, optoption='', TDstate_info=TDstate_info)
            
                job_state = gaussian_run.Exe_Gaussian.exe_Gaussian(JobNameState, self.timexe)
                print (job_state)
                # When the scf is performed, the obtained wavefunction is saved to chk file.
                # But for 'symm' and 'volume' that information is emply except for when the input is chk or fchk files.
                if scf_need or ReadFromchk:
                    gaussian_run.chk2fchk.Get_chklist()
                elif scf_need != True and ReadFromchk != True:
                    for f in glob.glob('./*.chk'):
                        os.remove(os.path.join('.', f))

                #output_prop = self.Extract_values(JobName, option_dict, Bondpair1, Bondpair2)
                #output_dic[i].update(output_prop)
                try:
                    output_prop = self.Extract_values(JobNameState, job_thisstate, Bondpair1, Bondpair2)
                    output_dic[i].update(output_prop)
                    output_dic[i][f'log_{i}'] = job_state
                except Exception as e:
                    job_state = "error"
                    print(e)
                    pass
            
                # for pka computation
                if 'pka' in job_thisstate:
                    #print("pka: ", output_dic_pka)
                    E_dH = output_dic[i]["Energy"][0]
                    output_dic[i]["pka"] = (E_dH - E_pH)*Eh2kJmol

                if 'tadf' in job_thisstate:
                    TADF_Eng = 0.0
                    S_Eext = output_dic[i-1]["MinEtarget"]
                    T_Eext = output_dic[i]["T_Min"]
                    TADF_Eng = S_Eext - T_Eext
                    output_dic[i]["Delta(S-T)"] = TADF_Eng

        output_sum = {}
        for i, j in enumerate(output_dic):
            #print(i, j)
            for k in output_dic[i].keys():
                output_sum.setdefault(k, output_dic[i][k])
            #for k in range(len(option_dict)):
    
                
        #output_dic["log"] = job_state

        # Convert fchk to xyz 
        if self.restart == False and scf_need == True:
            gaussian_run.Get_MolCoordinate_fchk.Get_fchklist2xyz()

        os.chdir("..")

        return(output_sum)


if __name__ == '__main__':
    usage ='Usage; %s infile' % sys.argv[0]

    try:
        infilename = sys.argv[1]
    except:
        print (usage); sys.exit()

    test_sdf = GaussianDFTRun('B3LYP', '3-21g*',8, 'opt',infilename)

    Mol = ['C', 'H', 'O', 'H']
    Molx = [0.54867, 1.16419, -0.70255, 1.16421]
    Moly = [0.00000, -0.94126, 0.00000, 0.94124]
    Molz = [0.00008, 0.00004, -0.00007, 0.00004]

    JobName = infilename
    test_sdf.make_input(JobName, 0, 1,  run_type='opt', Newinput=True,  Mol_atom=Mol, X=Molx, Y=Moly, Z=Molz,  newchk=None)

    test_sdf.make_input(JobName, 0, 1, run_type='freq', Newinput=False, readchk='all', newchk=None, solvent='0.52')

    newchk = JobName+'_IP'
    test_sdf.make_input(JobName, 1, 2, run_type='opt', Newinput=False, readchk='gomeguess', newchk=newchk, solvent='0')

    test_sdf.make_input(JobName, 0, 1, run_type='symm', Newinput=False, readchk='gome', newchk=None, solvent='0')

    newchk = JobName+'_exopt'
    test_sdf.make_input(JobName, 0, 1, run_type='opt', Newinput=False, TDDFT=True, TDstate=None, target=3, readchk='all', newchk=newchk)


