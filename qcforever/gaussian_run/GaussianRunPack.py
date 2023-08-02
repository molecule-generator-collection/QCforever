import os
import glob
import sys
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

        output = {}

        parse_log = gaussian_run.parse_log.parse_log(infilename)
        job_index, Links_split = parse_log.Check_task()

        # For extracting properties of molecule after SCF        
        if 'gs' in job_index:
            GS_lines = Links_split[job_index['gs']]
        
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
            Polar_iso, Polar_aniso = gaussian_run.Get_FreqPro.Extract_polar(polar_lines) 
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
                CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden = parse_log.Extract_ExcitedState(lines)
            output["uv"] = [WL_allowed, OS_allowed, CD_L_allowed, CD_OS_allowed]
            output["state_index"] = [State_allowed, State_forbidden]
            if self.ref_uv_path != '':
                ref_uv = {}
                print(f"Read the reference spectrum...{self.ref_uv_path}")
                ref_uv["uv"] = gaussian_run.UV_similarity.read_data(self.ref_uv_path)
                S, D = gaussian_run.UV_similarity.smililarity_dissimilarity(ref_uv["uv"][0], ref_uv["uv"][1], output["uv"][0], output["uv"][1])
                output["Similality/Disdimilarity"] = [S, D]
        
        if is_fluor or is_tadf:
            lines = Links_split[job_index['relaxAEstate']]
            S_Found, S_Egrd, S_Eext, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
                CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden = parse_log.Extract_ExcitedState(lines)
            # For Check internal coordinate
            MaxDisplace = gaussian_run.Get_MolInterCoordinate.Extract_InterMol(lines)
            output["MinEtarget"] = S_Eext
            output["Min_MaxDisplace"] = MaxDisplace
            output["fluor"] = [WL_allowed, OS_allowed, CD_L_allowed, CD_OS_allowed]
        
        if is_tadf:
            lines = Links_split[job_index['relaxFEstate']] 
            T_Found, _, T_Eext, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
                CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden  = parse_log.Extract_ExcitedState(lines)
            # For Check internal coordinate
            MaxDisplace = gaussian_run.Get_MolInterCoordinate.Extract_InterMol(lines)
            output["T_Min"] = T_Eext
            output["T_Min_MaxDisplace"] = MaxDisplace
            output["T_Phos"] = [WL_forbidden, OS_forbidden, CD_L_forbidden, CD_OS_forbidden]
            TADF_Eng = 0.0
            if S_Found and T_Found:
               TADF_Eng = S_Eext - T_Eext
            output["Delta(S-T)"] = TADF_Eng

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
            NState = target+4 if run_type == 'opt' else 20
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

    def job_construction(self, JobName, scf_need, option_dict, TotalCharge, SpinMulti, targetstate,
                        ReadFromchk, ReadFromsdf, ReadFromxyz,  
                        element=[], atomX=[], atomY=[], atomZ=[], optoption=''):

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

        if scf_need == False:
        #Only no-scf jobs
            if 'symm' in option_dict:
                if ReadFromchk:
                    self.make_input(JobName, TotalCharge, SpinMulti, 
                                    run_type='symm', Newinput=True, 
                                    geom_spec=Is_geom_spec, readchk='all') 
                if ReadFromsdf or ReadFromxyz:
                    self.make_input(JobName, TotalCharge, SpinMulti, 
                                    run_type='symm', Newinput=True, 
                                    Mol_atom=element, X=atomX, Y=atomY, Z=atomZ, 
                                    geom_spec=Is_geom_spec) 
            if 'volume' in option_dict:
                if ReadFromchk:
                    self.make_input(JobName, TotalCharge, SpinMulti, 
                                    run_type='volume', Newinput=False, 
                                    geom_spec=Is_geom_spec, readchk='all') 
                if ReadFromsdf or ReadFromxyz:
                    self.make_input(JobName, TotalCharge, SpinMulti, 
                                    run_type='volume', Newinput=False, 
                                    Mol_atom=element, X=atomX, Y=atomY, Z=atomZ, 
                                    geom_spec=Is_geom_spec) 
        else:
        #The first scf job
            run_opt = 'opt' if 'opt' in option_dict else ''
            if ReadFromchk and element == []:
                if Is_ChargeSpec == False and Is_SpinMultiSpec == False:
                    if optoption == '':
                        self.make_input(JobName, TotalCharge, SpinMulti, 
                                        run_type=run_opt, Newinput=True, 
                                        geom_spec=Is_geom_spec, readchk='all', solvent=self.solvent) 
                    else:
                        self.make_input(JobName, TotalCharge, SpinMulti, 
                                        run_type=run_opt, optoption=optoption, Newinput=True, 
                                        geom_spec=Is_geom_spec, readchk='all', solvent=self.solvent) 
                else:
                    if optoption == '':
                        self.make_input(JobName, TotalCharge, SpinMulti, 
                                        run_type=run_opt, Newinput=True, 
                                        geom_spec=Is_geom_spec, readchk='geomguess', solvent=self.solvent) 
                    else:
                        self.make_input(JobName, TotalCharge, SpinMulti, 
                                        run_type=run_opt, optoption=optoption, Newinput=True, 
                                        geom_spec=Is_geom_spec, readchk='geomguess', solvent=self.solvent) 
            elif ReadFromchk and element != []:
                if optoption == '':
                    self.make_input(JobName, TotalCharge, SpinMulti, 
                                    run_type=run_opt, Newinput=True, 
                                    Mol_atom=element, X=atomX, Y=atomY, Z=atomZ, 
                                    geom_spec=Is_geom_spec, readchk='guess', solvent=self.solvent) 
                else:
                    self.make_input(JobName, TotalCharge, SpinMulti, 
                                    run_type=run_opt, optoption=optoption, Newinput=True, 
                                    Mol_atom=element, X=atomX, Y=atomY, Z=atomZ, 
                                    geom_spec=Is_geom_spec, readchk='guess', solvent=self.solvent) 
            elif ReadFromsdf or ReadFromxyz:
                if optoption == '':
                    self.make_input(JobName, TotalCharge, SpinMulti, 
                                    run_type=run_opt, Newinput=True, 
                                    Mol_atom=element, X=atomX, Y=atomY, Z=atomZ, 
                                    geom_spec=Is_geom_spec, solvent=self.solvent) 
                else:
                    self.make_input(JobName, TotalCharge, SpinMulti, 
                                    run_type=run_opt, optoption=optpotion, Newinput=True, 
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
                self.make_input(JobName, TotalCharge, SpinMulti, 
                                run_type='freq', Newinput=False, 
                                readchk='all', solvent=self.solvent) 

            if 'polar' in option_dict: # polar == 1
                self.make_input(JobName, TotalCharge, SpinMulti, 
                                run_type='polar', Newinput=False, 
                                readchk='all', solvent=self.solvent) 

            if 'nmr' in option_dict: # nmr == 1
                self.make_input(JobName, TotalCharge, SpinMulti, 
                                run_type='nmr', Newinput=False, 
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

            if 'fluor' in option_dict:
                TDOpt_chk = f'{JobName}_ExOptAState{targetstate}' if targetstate != 1 else f'{JobName}_ExOptStateSinglet'
                self.make_input(JobName, TotalCharge, SpinMulti, 
                                scf='close', run_type='opt', TDDFT=True, TDstate='Singlet', target=targetstate, 
                                readchk='all', newchk=TDOpt_chk, solvent=self.solvent) 

            if 'tadf' in option_dict:
                TTDOpt_chk = f'{JobName}_ExOptFState{targetstate}' if targetstate != 1 else f'{JobName}_ExOptStateTriplet'
                #self.make_input(JobName, TotalCharge, SpinMulti, scf='close', run_type='opt', TDDFT=True, TDstate='Triplet', target=targetstate, readchk='all', oldchk=TDOpt_chk, newchk=TTDOpt_chk, solvent=self.solvent) 
                self.make_input(JobName, TotalCharge, SpinMulti, 
                                scf='close', run_type='opt', TDDFT=True, TDstate='Triplet', target=targetstate, 
                                readchk='all', newchk=TTDOpt_chk, solvent=self.solvent) 

    def run_gaussian(self):
        infilename = self.in_file
        option_line = self.value    
        options = option_line.split()
        option_dict = {}
        option_dict_Ex = {}
        option_dict_pka = {}
        targetstate = 1
        PreGauInput = infilename.split('.')
        JobName = PreGauInput[0]
        #GauInputName = PreGauInput[0]+'.com'    
        # File type of input?
        ReadFromchk = False 
        ReadFromsdf = False
        ReadFromxyz = False
        # Initialization of molecular coordinate
        atm = [] 
        X = [] 
        Y = [] 
        Z = []
        Bondpair1 = []
        Bondpair2 = []
        # Classify input
        if PreGauInput[1] == "sdf":
            ReadFromsdf = True 
            atm, X, Y, Z, TotalCharge, SpinMulti, Bondpair1, Bondpair2 = read_mol_file.read_sdf(infilename)
        elif PreGauInput[1] == "xyz":
            ReadFromxyz = True
            atm, X, Y, Z, TotalCharge, SpinMulti = read_mol_file.read_xyz(infilename)
        elif PreGauInput[1] == "chk":
            ReadFromchk = True 
        elif PreGauInput[1] == "fchk":
            TotalCharge, SpinMulti = gaussian_run.fchk2chk.Get_fchk(PreGauInput[0])
            ReadFromchk = True 
        else:
            print("Invalid input file")

        # Setting options
        for i in range(len(options)):
            option = options[i]
            if option.lower() == 'opt':
                option_dict['opt'] = True
                if '=' in option:
                    in_target = option.split("=")
                    optoption = in_target[-1]
            elif option.lower() == 'freq':
                option_dict['freq'] = True
            elif option.lower() == 'polar':
                option_dict['polar'] =  True
            elif option.lower() == 'nmr':
                option_dict['nmr'] = True
            elif 'uv' in option.lower():
                option_dict['uv'] = True
                if '=' in option:
                    ref_spectrum = option.split("=")
                    self.ref_uv_path = ref_spectrum[-1]
            elif option.lower() == 'energy':
                option_dict['energy'] = True
            elif option.lower() == 'homolumo':
                option_dict['homolumo'] = True
            elif option.lower() == 'dipole':
                option_dict['dipole'] = True
            elif option.lower() == 'deen':
                option_dict['deen'] = True
            elif option.lower() == 'stable2o2':
                option_dict['stable2o2'] = True
            elif 'fluor' in option.lower():
                option_dict['uv'] = True
                if SpinMulti == 1:
                    option_dict['fluor'] = True
                else:
                    option_dict_Ex['fluor'] = True
                if '=' in option:
                    in_target = option.split("=")
                    targetstate = int(in_target[-1])
            elif option.lower() == 'tadf':
                option_dict['uv'] = True
                if SpinMulti == 1:
                    option_dict['fluor'] = True
                    option_dict['tadf'] = True
                else:
                    option_dict_Ex['fluor'] = True
                    option_dict_Ex['tadf'] = True
            elif option.lower() == 'vip':
                option_dict['vip'] = True
            elif option.lower() == 'vea':
                option_dict['vea'] = True
            elif option.lower() == 'aip':
                option_dict['vip'] = True
                option_dict['aip'] = True
            elif option.lower() == 'aea':
                option_dict['vea'] = True
                option_dict['aea'] = True
            elif option.lower() == 'cden':
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
            elif option.lower() == 'volume':
                option_dict['volume'] = True
            else:
                print('invalid option: ', option)

        # Make work directory and change to the directory
        if os.path.isdir(JobName):
            shutil.rmtree(JobName)
        os.mkdir(JobName)
        if ReadFromchk:
            inchkfile = JobName + ".chk"
            shutil.move(inchkfile, JobName) 
        os.chdir(JobName)

        # Initialization of otuput dictionary.
        output_dic = {}

        if self.error == 1:
            check_logfile = JobName+'.log'
            option_dict_optcheck = {'opt': True, 'freq': True}
            scf_need=True
            self.job_construction(JobName, scf_need, option_dict_optcheck, TotalCharge, SpinMulti, 
                                targetstate, ReadFromchk, ReadFromsdf, ReadFromxyz, 
                                element=atm, atomX=X, atomY=Y, atomZ=Z, optoption='')
            job_state = "normal"
            job_state = gaussian_run.Exe_Gaussian.exe_Gaussian(JobName, self.timexe)

            try:
                output_prop = self.Extract_values(JobName, option_dict_optcheck, Bondpair1, Bondpair2)
            except Exception as e:
                job_state = "error"
                print(e)
                pass

            with open(check_logfile,'r') as ifile:
                CKF_lines = ifile.readlines()
            Mol_atom, Mol_X, Mol_Y, Mol_Z = gaussian_run.Get_MolCoordinate.Extract_Coordinate(CKF_lines)

            print(output_prop)
            Freq_CK = np.array(output_prop['freq'])
            Freq_min = Freq_CK.min()
            Freq_min_index = np.argmin(Freq_CK)

            ReadFromchk=True 
            ReadFromsdf=False 
            ReadFromxyz=False
        
            if Freq_min > 0:
                print (f'Mnimum frequency is positive...{Freq_min} cm^-1')
                #Checking the remain task (option)
                all_options = set(option_dict.keys())
                remain_keys = all_options - set(option_dict_optcheck.keys())
        
                if len(remain_keys) >= 1:
                    self.error -= 1
                    output_dic.update(output_prop)
                    if 'freq' in option_dict:
                        del option_dict['freq']
                    if 'opt' in option_dict:
                        del option_dict['opt']
                    #for the next computation
                else:
                    output_dic.update(output_prop)
                    return(output_dic)
            else:
                print (f'Negative frequency is detected...{Freq_min} cm^-1')

                ImFreq_modeX = output_prop['freqmode'][0][Freq_min_index]
                ImFreq_modeY = output_prop['freqmode'][1][Freq_min_index]
                ImFreq_modeZ = output_prop['freqmode'][2][Freq_min_index]

                X = Mol_X+ 0.05* ImFreq_modeX
                Y = Mol_Y+ 0.05* ImFreq_modeY
                Z = Mol_Z+ 0.05* ImFreq_modeZ

                self.job_construction(JobName, scf_need, option_dict_optcheck, TotalCharge, SpinMulti, 
                                    targetstate, ReadFromchk, ReadFromsdf, ReadFromxyz, 
                                    element=atm, atomX=X, atomY=Y, atomZ=Z, optoption='')
                job_state = "normal"
                job_state = gaussian_run.Exe_Gaussian.exe_Gaussian(JobName, self.timexe)

                try:
                    output_prop = self.Extract_values(JobName, option_dict_optcheck, Bondpair1, Bondpair2)
                except Exception as e:
                    job_state = "error"
                    print(e)
                    pass

                output_dic.update(output_prop)

                if 'freq' in option_dict:
                    del option_dict['freq']
                if 'opt' in option_dict:
                    del option_dict['opt']

                self.error -= 1

        # Initialization of molecular coordinate
            atm = [] 
            X = [] 
            Y = [] 
            Z = []

        if self.error == 0:
            scf_need = self.check_scf_need(option_dict)
            self.job_construction(JobName, scf_need, option_dict, TotalCharge, SpinMulti, 
                                targetstate, ReadFromchk, ReadFromsdf, ReadFromxyz, 
                                element=atm, atomX=X, atomY=Y, atomZ=Z, optoption='')

            job_state = "normal"
            job_state = gaussian_run.Exe_Gaussian.exe_Gaussian(JobName, self.timexe)
            # When the scf is performed, the obtained wavefunction is saved to chk file.
            # But for 'symm' and 'volume' that information is emply except for when the input is chk or fchk files.
            if scf_need or ReadFromchk:
                gaussian_run.chk2fchk.Get_chklist()
            elif scf_need != True and ReadFromchk != True:
                for f in glob.glob('./*.chk'):
                    os.remove(os.path.join('.', f))

            #output_prop = self.Extract_values(JobName, option_dict, Bondpair1, Bondpair2)
            #output_dic.update(output_prop)
            try:
                output_prop = self.Extract_values(JobName, option_dict, Bondpair1, Bondpair2)
                output_dic.update(output_prop)
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
                GS_fchk = JobName+".fchk"
                with open(GS_fchk,'r') as ifile:
                    GS_lines = ifile.readlines()
                TotalCharge, SpinMulti, Mol_atom, Mol_X, Mol_Y, Mol_Z = gaussian_run.Get_MolCoordinate_fchk.Extract_MolCoord(GS_lines)
                DeHMol_atom = np.delete(Mol_atom,Index_MaxProtic) 
                DeHMol_X = np.delete(Mol_X,Index_MaxProtic)
                DeHMol_Y = np.delete(Mol_Y,Index_MaxProtic)
                DeHMol_Z = np.delete(Mol_Z,Index_MaxProtic)
                JobName_DeHMol = JobName + "_DeH"
            
                self.make_input(JobName_DeHMol, TotalCharge-1, SpinMulti, 
                            run_type='opt', Newinput=True, 
                            Mol_atom=DeHMol_atom, X=DeHMol_X, Y=DeHMol_Y, Z=DeHMol_Z, readchk=False, solvent=self.solvent) 
            
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
            
             #for fluor == 1 or tadf == 1 for open shell
            if 'fluor' in option_dict_Ex or 'tadf' in option_dict_Ex: 
                TotalCharge, SpinMulti = gaussian_run.fchk2chk.Get_fchk(PreGauInput[0])
                output_dic_Ex = {}
                compute_state = output_dic["state_index"][0][int(targetstate)-1] 
                JobName_ExOpt = JobName + '_ExOptAState'+f'{targetstate}'
            
                self.make_input(JobName_ExOpt, TotalCharge, SpinMulti, 
                                run_type='opt', Newinput=True, TDDFT=True, target=compute_state, 
                                readchk='all', oldchk=JobName, newchk=JobName_ExOpt, solvent=self.solvent) 
            
                if 'tadf' in option_dict_Ex: #tadf == 1
                    ExOptFState_chk=  JobName + '_ExOptFState'+f'{targetstate}'
                    compute_state = output_dic["state_index"][1][int(targetstate)-1] 
                    self.make_input(JobName_ExOpt, TotalCharge, SpinMulti, 
                                    run_type='opt', TDDFT=True, target=compute_state, 
                                    readchk='all', oldchk=JobName_ExOpt, newchk=ExOptFState_chk, solvent=self.solvent) 
            
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
        if self.restart == False and scf_need == True:
            gaussian_run.Get_MolCoordinate_fchk.Get_fchklist2xyz()

        os.chdir("..")

        return(output_dic)


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


