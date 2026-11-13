import os
import re
import glob
import sys
import copy
import shutil
import pickle
from subprocess import CalledProcessError

import numpy as np

from pathlib import Path

from qcforever import gaussian_run
from qcforever.util import calculation_cache
from qcforever.util.calculation_options import cden_fukui, excited_options, requires_hessian
from qcforever.util.fukui import condensed_fukui
from qcforever.util.imaginary_modes import NO_STABLE_MINIMUM, StationaryPoint, follow_imaginary_modes
from qcforever.util import (
    read_mol_file,
    check_resource,
    UV_similarity,
    Spectrum_similarity,
    job_cleanup,
    job_timeout,
)
from qcforever.gaussian_run import frequency_seed
from qcforever.gaussian_run.solvent import resolve_gaussian_solvent
from qcforever.util.quadrupole import quadrupole_at_center_of_mass, quadrupole_properties


Eh2kJmol = 2625.5
Eh2eV = 27.211


class GaussianDFTRun:
    def __init__(self, 
                functional, 
                basis, 
                nproc, 
                value, 
                in_file, 
                solvent=None,
                error=0, 
                restart=True, 
                pklsave=False):
        self.in_file = str(Path(in_file).expanduser().absolute())
        self.functional = functional.lower()
        self.basis = basis.lower()
        self.nproc = check_resource.respec_cores(nproc)
        self.value = value
        # Validate before any Gaussian input or process can be created. Named
        # solvents are normalized to the exact spelling expected by Gaussian;
        # mapped SMILES are converted to that same name.
        self.solvent = resolve_gaussian_solvent(solvent)
        self.restart = restart
        self.error = error
        self.pklsave = pklsave
        self.conformer_config = None
        self.mem = ''
        self.timexe = 60 * 60 * 80
        # Overall QCforever wall-clock limit in seconds.  None means unlimited.
        # ``timexe`` remains the per-Gaussian-calculation limit.
        self.timejob = None
        self.SpecTotalCharge = np.nan
        self.SpecSpinMulti = np.nan
        self.ref_uv_path = ''
        self.ref_spectrum_paths = {}
        self.geom_spec = {}
        self.para_functional = []
        self.functional_parameter_optimization = None
        self.functional_parameter_prediction = None
        self.reuse = True
        self.cleanup_cache = None
        self.stable_from = None

    def check_scf_need(self, option):
        all_keys = set(option.keys())
        nonSCF_keys = ['symm', 'volume']
        if len(all_keys - set(nonSCF_keys)) >= 1:
            scf_need = True
        else:
            scf_need = False

        return scf_need

    #@staticmethod
    #def _normal_termination(jobname):
    #    """Return whether the final Gaussian Link1 job terminated normally."""
    #    try:
    #        with open(f"{jobname}.log", encoding="utf-8", errors="replace") as logfile:
    #            logtext = logfile.read()
    #        return (
    #            logtext.rfind("Normal termination of Gaussian")
    #            > logtext.rfind("Error termination")
    #        )
    #    except OSError:
    #        return False

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
        is_nac = option_dict['nac'] if 'nac' in option_dict else False
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
            output["GS_MaxDisplace"] = gaussian_run.Get_MolInterCoordinate.Extract_InterMol(GS_lines) 
        
        if is_homolumo:
            with open(fchkname, 'r') as ifile:
                fchk_lines = ifile.readlines()
            NumAlphaElec, NumBetaElec, AlphaEigenVal, BetaEigenVal = gaussian_run.Get_MOEnergy_fchk.Extract_MO(fchk_lines)
            # NumAlphaElec, NumBetaElec, AlphaEigenVal, BetaEigenVal = gaussian_run.Get_MOEnergy.Extract_MO(GS_lines)
            if BetaEigenVal == []:
                Alpha_gap = Eh2eV * (AlphaEigenVal[NumAlphaElec]-AlphaEigenVal[NumAlphaElec-1])
                output["homolumo"] = Alpha_gap
                output["MO"] = [AlphaEigenVal[:NumAlphaElec], AlphaEigenVal[NumAlphaElec:]]
            else:
                Alpha_gap = Eh2eV * (AlphaEigenVal[NumAlphaElec]-AlphaEigenVal[NumAlphaElec-1])
                Beta_gap = Eh2eV * (BetaEigenVal[NumBetaElec]-BetaEigenVal[NumBetaElec-1])
                output["homolumo"] = [Alpha_gap, Beta_gap]
                output["Alpha_MO"] = [AlphaEigenVal[:NumAlphaElec], AlphaEigenVal[NumAlphaElec:]]
                output["Beta_MO"] = [BetaEigenVal[:NumBetaElec], BetaEigenVal[NumBetaElec:]]
        
        if is_dipole:
            output["dipole"] = parse_log.Extract_DipoleMoment(GS_lines)
            try:
                moments = parse_log.Extract_QuadrupoleMoment(GS_lines, fchkname)
                output.update(quadrupole_at_center_of_mass(**moments, source="Gaussian"))
            except ValueError as error:
                output.update(quadrupole_properties(
                    source="Gaussian", origin="center_of_mass",
                    status="unavailable", reason=str(error)))
        
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
            if getattr(self, 'pklsave', False) and any(' orientation:' in line for line in freq_lines):
                output['frequency_context'] = frequency_seed.frequency_context(parse_log, freq_lines)
            for spectrum_type, intensities in (("IR", IR), ("Raman", Raman)):
                reference_path = self.ref_spectrum_paths.get(spectrum_type)
                if reference_path:
                    output.setdefault("Spectrum_similarity", {})[spectrum_type] = (
                        Spectrum_similarity.compare_with_file(
                            reference_path, Freq, intensities, spectrum_type
                        )
                    )
        
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
            reference_path = self.ref_spectrum_paths.get("NMR")
            if reference_path:
                nmr_shifts = [
                    shift
                    for element, shift in zip(Element, ppm)
                    if element in ("H", "C", "Si")
                ]
                output.setdefault("Spectrum_similarity", {})["NMR"] = (
                    Spectrum_similarity.compare_with_file(
                        reference_path,
                        nmr_shifts,
                        [1.0] * len(nmr_shifts),
                        "NMR",
                    )
                )
        
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
            if option_dict.get('fukui'):
                extract_charges = gaussian_run.Get_ChargeSpin.Extract_ChargeSpin
                output['fukui'] = condensed_fukui(
                    output['cden'], extract_charges(IP_lines), extract_charges(EA_lines))
        
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
                output["relaxedIP_MaxDisplace"] = gaussian_run.Get_MolInterCoordinate.Extract_InterMol(PC_lines)
                output["aip"] = [Eh2eV*(PC_Energy_SS[0]-GS_Energy), Eh2eV*(PC_Energy_SS[0]-VNP_Energy_SS[0]), PC_Energy_SS[1], VNP_Energy_SS[1]]
            if is_aea:
                NC_lines = Links_split[job_index['NC_line']]
                VNN_lines = Links_split[job_index['VNN_line']]
                NC_Energy_SS = parse_log.Extract_SCFEnergy(NC_lines)
                VNN_Energy_SS = parse_log.Extract_SCFEnergy(VNN_lines)
                # For Check internal coordinate
                output["relaxedEA_MaxDisplace"] = gaussian_run.Get_MolInterCoordinate.Extract_InterMol(NC_lines)
                # Normal electronic affinity calculation
                output["aea"] = [Eh2eV*(GS_Energy-NC_Energy_SS[0]), Eh2eV*(VNN_Energy_SS[0]-NC_Energy_SS[0]), NC_Energy_SS[1], VNN_Energy_SS[1]]
        
        if is_satkoopmans:
            dSCF_vip = output["vip"][0]
            dSCF_vea = output["vea"][0]
            Alpha_eHOMO = AlphaEigenVal[NumAlphaElec-1] 
            Beta_eHOMO = BetaEigenVal[NumBetaElec-1]
            Alpha_eLUMO = AlphaEigenVal[NumAlphaElec] 
            Beta_eLUMO = BetaEigenVal[NumBetaElec]
            # max(finite, NaN) can hide an invalid spin-channel HOMO.
            if not np.all(np.isfinite([dSCF_vip, Alpha_eHOMO, Beta_eHOMO])):
                raise ValueError("The IP and both HOMO energies must be finite")
            Delta_vip = dSCF_vip + Eh2eV*max(Alpha_eHOMO, Beta_eHOMO)
            Delta_vea = dSCF_vea + Eh2eV*min(Alpha_eLUMO, Beta_eLUMO)
            output["satkoopmans"] = [Delta_vip, Delta_vea]
        
        if is_uv:
            lines = Links_split[job_index['uv_line']] 
            _, _, _, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
                CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden, \
                mu_allowed, mu_forbidden, mm_allowed, mm_forbidden, theta_allowed, theta_forbidden, g_allowed, g_forbidden \
                = parse_log.Extract_ExcitedState(lines)
            output["uv"] = [WL_allowed, OS_allowed, CD_L_allowed, CD_OS_allowed]
            output["state_index"] = [State_allowed, State_forbidden]
            output["CD_mu"] = [mu_allowed, mu_forbidden]
            output["CD_mm"] = [mm_allowed, mm_forbidden]
            output["CD_theta"] = [theta_allowed, theta_forbidden]
            output["CD_g"] = [g_allowed, g_forbidden]
            if self.ref_uv_path != '':
                ref_uv = {}
                print(f"Read the reference spectrum...{self.ref_uv_path}")
                ref_uv["uv"] = UV_similarity.read_data(self.ref_uv_path)
                S, D, Wasser= UV_similarity.smililarity_dissimilarity(ref_uv["uv"][0], ref_uv["uv"][1], output["uv"][0], output["uv"][1])
                output["Similarity/Dissimilarity"] = [S, D]
                output["Wasserstein_distance"] = [Wasser]
        
        if is_fluor:
            lines = Links_split[job_index['relaxAEstate']]
            S_Found, S_Egrd, S_Eext, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
                CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden, \
                mu_allowed, mu_forbidden, mm_allowed, mm_forbidden, theta_allowed, theta_forbidden, g_allowed, g_forbidden \
                = parse_log.Extract_ExcitedState(lines)
            # For Check internal coordinate
            output["MinEtarget"] = S_Eext
            output["Min_MaxDisplace"] = gaussian_run.Get_MolInterCoordinate.Extract_InterMol(lines)
            output["fluor"] = [WL_allowed, OS_allowed, CD_L_allowed, CD_OS_allowed]
            #mu, mm, theta, g_cal = parse_log.extract_transtion_EMmoment(lines)
            output["CPL_mu"] = mu_allowed 
            output["CPL_mm"] = mm_allowed 
            output["CPL_theta"] = theta_allowed
            output["CPL_g"] = g_allowed
        
        if is_tadf:
            lines = Links_split[job_index['relaxFEstate']] 
            T_Found, _, T_Eext, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
                CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden, \
                mu_allowed, mu_forbidden, mm_allowed, mm_forbidden, theta_allowed, theta_forbidden, g_allowed, g_forbidden \
                = parse_log.Extract_ExcitedState(lines)
            # For Check internal coordinate
            output["T_Min"] = T_Eext
            output["T_Min_MaxDisplace"] =gaussian_run.Get_MolInterCoordinate.Extract_InterMol(lines)
            output["T_Phos"] = [WL_forbidden, OS_forbidden, CD_L_forbidden, CD_OS_forbidden]
            output["CPL_T_mu"] = mu_forbidden 
            output["CPL_T_mm"] = mm_forbidden 
            output["CPL_T_theta"] = theta_forbidden
            output["CPL_T_g"] = g_forbidden
            #TADF_Eng = 0.0
            #if S_Found and T_Found:
            #   TADF_Eng = S_Eext - T_Eext
            #output["Delta(S-T)"] = TADF_Eng

        if is_nac:
            # parse_log.classify_task stores the final NAC Link1 under the
            # uppercase ``NAC`` key and Extract_NAC returns (max, rms).
            # Keep errors local so the fluorescence data remains serializable.
            try:
                nac_lines = Links_split[job_index['NAC']]
                nac_max, nac_rms = parse_log.Extract_NAC(nac_lines)
                if nac_max is None or nac_rms is None:
                    raise ValueError('NAC values were not found in the NAC Link1 output.')
                output['NAC'] = [nac_max, nac_rms]
            except (AttributeError, IndexError, KeyError, ValueError) as error:
                output['NAC_error'] = str(error)

        return output

    def MakeSolventLine(self):
        s = ''
        s_solvent = ''
        if self.solvent is None or self.solvent == '0':
            return s, s_solvent
        try:
            float(self.solvent)
        except ValueError:
            print('Solvent effect is included by PCM')
            s += f'SCRF=(PCM, solvent={self.solvent})\n'
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
        TDDFT=False, TDstate=None, target=1, nac=False,
        readchk=None, oldchk=None, newchk=None, solvent=None, mode_scan=False):

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
        method_basis = ' ChkBasis' if self.basis == 'fchk' else f'/{self.basis}'
        if scf == 'open' or SpinMulti != 1:
            line_method   = f'#u{self.functional}{method_basis} test {line_iop_functional}\n'
        elif scf == 'close' and SpinMulti == 1:
            line_method = f'#r{self.functional}{method_basis} test {line_iop_functional}\n'
        else:
            line_method = f'#{self.functional}{method_basis} test {line_iop_functional}\n'

        if TDDFT:
            NState = target+4 if (run_type == 'opt' or run_type == 'freq') else 20
            td_options = f'NStates={NState},root={target}'
            if TDstate is not None:
                td_options = f'NStates={NState},{TDstate},root={target}'
            if nac:
                td_options += ',NAC'
            line_method += f'TD({td_options})\n'
#            line_method += f'TD(Nstate={NState}, root={target})\n'if TDstate == None else f'TD(Nstate={NState}, {TDstate}, root={target})\n'
        
        #Section for solvent
        SCRF = ''
        SCRF_read = ''
        solvent = resolve_gaussian_solvent(solvent)
        if solvent is not None and solvent != '0':
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
            if not Path(GauInputName).is_file():
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
            if mode_scan:
                input_s += 'NoSymm SCF=VeryTight\n'

            #For non SCF run_type
            if run_type == 'symm' or run_type == 'volume':
                input_s += 'Guess=Only'
                if run_type == 'volume':
                    input_s += ' volume\n'
                if run_type == 'symm':
                    input_s += ' Symmetry=loose\n'
            
            if run_type == 'opt':
                if TDDFT:
                    if optoption == '':
                        input_s += 'Opt=(Maxcycles=100) NoSymm\n'
                    else:
                        input_s += f'Opt=({optoption}, MaxCycles=100) NoSymm\n'
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

            if solvent is not None and solvent != '0':
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
            if solvent is not None and solvent != '0':
                input_s += SCRF_read
            if geom_spec:
                input_s += line_GeomConstrainSpec
                input_s += '\n'

            ofile.write(input_s)

    def chain_job(self, JobName, scf_need, job_dict, TotalCharge, SpinMulti, targetstate,
                        ReadFrom, element=[], atomX=[], atomY=[], atomZ=[], optoption='', TDstate_info=[],
                        guess_chk=None):

        option_dict = copy.deepcopy(job_dict)
        # opt=... belongs to the ground-state optimization. An excited state
        # takes force-constant controls exclusively from its own option group.
        if 'fluor' in option_dict or 'tadf' in option_dict:
            optoption = option_dict.pop('excited_optoption', '')

        print (option_dict)

        print(f'Target state: {targetstate}, Charge: {TotalCharge}, Spin Multiplicity: {SpinMulti}')

        # Setting molecular charge and multiplicity
        Is_ChargeSpec = False
        Is_SpinMultiSpec = False
        if np.isnan(self.SpecTotalCharge) !=  True:
            Is_ChargeSpec = True
        #    TotalCharge = self.SpecTotalCharge
        if np.isnan(self.SpecSpinMulti) != True:
            Is_SpinMultiSpec = True
        #    SpinMulti = self.SpecSpinMulti

        # Setting  geometric constrain
        if self.geom_spec != {}: 
            Is_geom_spec = True
        else:
            Is_geom_spec = False

        #For setting TD-DFT
        scf_tag='open'
        td_tag = False
        TDstate_tag = None 
        run_nac = False
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
        if 'nac' in option_dict:
                del option_dict['nac']
                run_nac = True
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
                                        run_type=run_opt, TDDFT=td_tag, TDstate=TDstate_tag, target=targetstate, 
                                        Newinput=True, geom_spec=Is_geom_spec, readchk='all', solvent=self.solvent) 
                    else:
                        self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag, 
                                        run_type=run_opt, TDDFT=td_tag, TDstate=TDstate_tag, target=targetstate, 
                                        optoption=optoption, 
                                        Newinput=True, geom_spec=Is_geom_spec, readchk='all', solvent=self.solvent) 
                else:
                    if optoption == '':
                        self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag, 
                                        run_type=run_opt, TDDFT=td_tag, TDstate=TDstate_tag, target=targetstate, 
                                        Newinput=True, geom_spec=Is_geom_spec, readchk='geomguess', solvent=self.solvent) 
                    else:
                        self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag, 
                                        run_type=run_opt, TDDFT=td_tag, TDstate=TDstate_tag, target=targetstate,  
                                        optoption=optoption, Newinput=True, 
                                        geom_spec=Is_geom_spec, readchk='geomguess', solvent=self.solvent) 
            elif ReadFrom == 'chk' and element != []:
                if optoption == '':
                    self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag, 
                                    run_type=run_opt, TDDFT=td_tag, TDstate=TDstate_tag, target=targetstate,  
                                    Newinput=True, 
                                    Mol_atom=element, X=atomX, Y=atomY, Z=atomZ, 
                                    geom_spec=Is_geom_spec, readchk='guess', solvent=self.solvent) 
                else:
                    self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag, 
                                    run_type=run_opt, TDDFT=td_tag, TDstate=TDstate_tag, target=targetstate, 
                                    optoption=optoption, Newinput=True, 
                                    Mol_atom=element, X=atomX, Y=atomY, Z=atomZ, 
                                    geom_spec=Is_geom_spec, readchk='guess', solvent=self.solvent) 
            elif ReadFrom == 'sdf' or ReadFrom == 'xyz':
                self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag,
                                run_type=run_opt, TDDFT=td_tag, TDstate=TDstate_tag, target=targetstate,
                                optoption=optoption, Newinput=True,
                                Mol_atom=element, X=atomX, Y=atomY, Z=atomZ,
                                geom_spec=Is_geom_spec, solvent=self.solvent,
                                readchk='guess' if guess_chk else None, oldchk=guess_chk,
                                newchk=JobName if guess_chk else None)

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
                                run_type='freq', TDDFT=td_tag, TDstate=TDstate_tag, target=targetstate,  
                                Newinput=False, 
                                readchk='all', solvent=self.solvent) 

            if 'polar' in option_dict: # polar == 1
                self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag, 
                                run_type='polar', TDDFT=td_tag, TDstate=TDstate_tag, target=targetstate, 
                                Newinput=False, 
                                readchk='all', solvent=self.solvent) 

            if 'nmr' in option_dict: # nmr == 1
                self.make_input(JobName, TotalCharge, SpinMulti, scf=scf_tag, 
                                run_type='nmr', TDDFT=td_tag, TDstate=TDstate_tag, target=targetstate, 
                                Newinput=False, 
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

            if run_nac:
                self.make_input(
                JobName, TotalCharge, SpinMulti, scf=scf_tag,
                run_type='nac',
                TDDFT=True, TDstate=TDstate_tag, target=targetstate,
                nac=True,
                Newinput=False, readchk='all',
                solvent=self.solvent,
                )


    def _stable_ground(self, name, charge, mult, read_from, atoms, x, y, z, optoption, stable_from=None):
        options = {'opt': True, 'freq': True, 'energy': True}

        def execute(job, error=0):
            status = gaussian_run.Exe_Gaussian.exe_Gaussian(job, self.timexe, error)
            if status == 'timeout':
                raise job_timeout.QCforeverTimeoutError(f'Gaussian timed out: {job}.log')
            if status != 'normal':
                raise RuntimeError(f'Gaussian ended with {status!r}: {job}.log')

        def calculate(job, source, elements, coordinates, optimizer, guess_chk=None):
            self.chain_job(job, True, options, charge, mult, 0, source, element=elements,
                           atomX=coordinates[0], atomY=coordinates[1], atomZ=coordinates[2],
                           optoption=optimizer, guess_chk=guess_chk)
            execute(job, self.error)
            parsed = gaussian_run.parse_log.parse_log(job + '.log')
            index, links = parsed.Check_task()
            if 'freq_line' not in index or 'ts' not in index:
                raise ValueError(f'Missing optimization or frequency step: {job}.log')
            if not any('Stationary point found' in line for line in links[index['ts']]):
                raise ValueError(f'Ground-state optimization did not converge: {job}.log')
            props = self.Extract_values(job, options, [], [])
            context = frequency_seed.frequency_context(parsed, links[index['freq_line']])
            props['frequency_context'] = context
            return StationaryPoint(job, context['atoms'], np.array(context['coordinates_angstrom']),
                                   props['Energy'][0], props['freq'], props['freqmode'], props)

        point = None
        if stable_from or (self.reuse and len(atoms) and not self.geom_spec):
            if not len(atoms) or self.geom_spec or self.basis == 'fchk':
                raise ValueError('stable_from requires explicit, unconstrained input coordinates and basis')
            previous_log = calculation_cache.previous_output(name + '.log')
            sources = [calculation_cache.previous_output(stable_from)] if stable_from else [
                calculation_cache.previous_output(Path(self.in_file).stem + '.pkl'), previous_log,
            ]
            solvent_route, solvent_read = self.MakeSolventLine()
            parameters = (gaussian_run.Gen_IOPline_functional.functional_para(self.functional, self.para_functional)
                          if self.para_functional else '')
            for source in sources:
                if not stable_from and not source.is_file():
                    continue
                try:
                    point = frequency_seed.load(
                        source, name, atoms=atoms, coordinates=np.column_stack((x, y, z)),
                        charge=charge, mult=mult, method=f'u{self.functional}/{self.basis}',
                        solvent_route=solvent_route, solvent_read=solvent_read,
                        parameter_route=parameters, companion=previous_log,
                        allow_completed_stable=not stable_from)
                except (OSError, KeyError, IndexError, TypeError, ValueError, EOFError, pickle.UnpicklingError) as exc:
                    if stable_from:
                        raise ValueError(f'Cannot reuse {source}: {exc}') from exc
                    continue
                print(f'Reusing ground-state vibrations from {source}')
                break
        if point is None:
            point = calculate(name, read_from, atoms, (x, y, z), optoption)

        probe_checkpoints = set()

        def reference_guess(center, job):
            checkpoint = Path(center.name + '.chk')
            if checkpoint.is_file():
                return checkpoint.name
            source = center.properties.get('stable_source')
            if source:
                native = Path(source).with_suffix('.chk')
                formatted = native.with_suffix('.fchk')
                has_native = native.is_file() and native.stat().st_size > 0
                if has_native or formatted.is_file():
                    local = Path(job + '_guess.chk')
                    probe_checkpoints.add(local)
                    try:
                        if has_native:
                            shutil.copyfile(native, local)
                        else:
                            gaussian_run.fchk2chk.convert(formatted, local)
                        return local.name
                    except job_timeout.QCforeverTimeoutError:
                        raise
                    except (OSError, CalledProcessError, RuntimeError) as exc:
                        print(f'Could not prepare saved orbitals from {source}: {exc}')
            print(f'{job}: starting the reference single point with a fresh orbital guess')
            return None

        def single_point(label, center, coordinates, reference_label):
            job = f'{name}_stable_{label}'
            guess = (f'{name}_stable_{reference_label}.chk' if reference_label else reference_guess(center, job))
            probe_checkpoints.add(Path(job + '.chk'))
            self.make_input(job, charge, mult, run_type='', Newinput=True, Mol_atom=center.atoms,
                            X=coordinates[:, 0], Y=coordinates[:, 1], Z=coordinates[:, 2],
                            solvent=self.solvent, mode_scan=True, readchk='guess' if guess else None,
                            oldchk=guess, newchk=job if guess else None)
            execute(job)
            parsed = gaussian_run.parse_log.parse_log(job + '.log')
            try:
                return parsed.Extract_SCFEnergy(parsed.lines)[0]
            except IndexError as exc:
                raise ValueError(f'Missing single-point energy: {job}.log') from exc

        def evaluate(label, elements, coordinates, guess_label):
            # A new input uses the displaced Cartesian geometry directly.
            return calculate(f'{name}_stable_{label}', 'xyz', elements, coordinates.T, '',
                             guess_chk=f'{name}_stable_{guess_label}.chk')

        reused_source = point.properties.get('stable_source')
        try:
            if not any(spec[1] == 'F' for spec in self.geom_spec.values()):
                point = follow_imaginary_modes(point, single_point, evaluate)
        finally:
            # The calculation cache retains successful probes for a restart.
            for checkpoint in probe_checkpoints:
                checkpoint.unlink(missing_ok=True)
        if reused_source:
            point.properties['stable_source'] = reused_source
        if len(atoms) and not self.geom_spec:
            point.properties['stable_input_geometry'] = np.column_stack((x, y, z)).tolist()
        if point.name != name:
            for suffix in ('.chk', '.log', '.com'):
                shutil.copyfile(point.name + suffix, name + suffix)
        return point


    def SpinMulti_scan(self, JobName,  ReadFrom, GivenSpinMulti, TotalCharge, atm, X, Y, Z, targetstate=0, TDstate_info=[]):
        print('Try to optimize the spin state of the ground state!')
        
        SpinMulti_list = []
        SpinDiff_list = []
        Energy_list = []
        
        if GivenSpinMulti%2 != 0:
            print("The number of electrons is even.")
            InitSpinMulti = 1
        else:
            print("The number of electrons is odd.")
            InitSpinMulti = 2
        
        for i in range(InitSpinMulti, InitSpinMulti+2*3, 2):
        
            JobNameState = JobName + f'_State{targetstate}_{i}_{TotalCharge}'

            if ReadFrom == 'chk':
                inchkfile = JobName + '.chk'
                Jobchk = JobNameState  + '.chk' 
                shutil.copy(inchkfile, Jobchk) 
        
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
        
        if Index_LowSpinDiff == Index_LowEnergy:
            SpinMulti = SpinMulti_list[Index_LowEnergy]
        elif max(SpinDiff_list) < 0.1:
            SpinMulti = SpinMulti_list[Index_LowEnergy]
        else:
            SpinMulti = SpinMulti_list[Index_LowSpinDiff]

        #overwrite/save the chk file as input chk file
        ReNameState =  JobName + f'_State{targetstate}_{SpinMulti}_{TotalCharge}'
        ReNameStateChk = ReNameState + '.chk'
        inchkfile = JobName + '.chk'
        shutil.copy(ReNameStateChk, inchkfile) 

        return SpinMulti

    def LC_para_BOopt(self, JobName, ReadFrom, atm, X, Y, Z, TotalCharge, SpinMulti):
        """Tune mu using the IP residual (eV), retaining the best valid trial."""
        from bayes_opt import BayesianOptimization
        from bayes_opt import UtilityFunction

        tolerance_ev = 0.01
        best_valid_parameters = None
        best_valid_error_ev = float('inf')
        valid_evaluations = 0

        def KTLC_BB(mu):
            nonlocal best_valid_parameters, best_valid_error_ev, valid_evaluations
            self.para_functional = [mu]

            option_dict_lcparacheck = {
                'homolumo': True, 'vip': True, 'vea': True, 'satkoopmans': True}
            self.chain_job(JobName, True, option_dict_lcparacheck, TotalCharge, SpinMulti,
                           0, ReadFrom, element=atm, atomX=X, atomY=Y, atomZ=Z,
                           optoption='', TDstate_info=[])
            job_state = gaussian_run.Exe_Gaussian.exe_Gaussian(JobName, self.timexe)
            if job_state != "normal":
                return -10

            try:
                # Only the neutral checkpoint supplies the HOMO. A failed
                # conversion must not let a previous trial's fchk be reused.
                job_timeout.run(['formchk', f'{JobName}.chk'], check=True)
                output_prop = self.Extract_values(
                    JobName, option_dict_lcparacheck, Bondpair1=[], Bondpair2=[])
                error_ev = abs(float(output_prop['satkoopmans'][0]))
                target = -error_ev**2
                if not np.isfinite(error_ev) or not np.isfinite(target):
                    raise ValueError("Non-finite Koopmans residual")
            except job_timeout.QCforeverTimeoutError:
                raise
            except Exception as error:
                print(f"KTLC evaluation failed at mu={mu}: {error}")
                return -10

            # The BO failure penalty is not a measured error and must never
            # participate in selection of the parameter used for production.
            valid_evaluations += 1
            if error_ev < best_valid_error_ev:
                best_valid_error_ev = error_ev
                best_valid_parameters = [mu]
            return target

        optimizer = BayesianOptimization(
            f=KTLC_BB, pbounds={'mu': (0, 1)}, verbose=1, random_state=1,
        )

        # Preserve the existing initial points and maximum of 52 evaluations.
        init_mu = [0.15, 0.35, 0.65]
        for mu in init_mu:
            optimizer.set_bounds(new_bounds={"mu": (mu, mu)})
            optimizer.maximize(init_points=1, n_iter=0)
        optimizer.set_bounds(new_bounds={"mu": (0, 1)})

        if best_valid_error_ev > tolerance_ev:
            i = 3
            acquisition_function = UtilityFunction(kind="ei", xi=1e-4)
            while best_valid_error_ev > tolerance_ev:
                optimizer.maximize(
                    init_points=0, n_iter=1,
                    acquisition_function=acquisition_function,
                )
                i += 1
                if i > 51:
                    break

        if best_valid_parameters is None:
            selected_parameters = []
            state = 'default_no_valid'
            selection_reason = 'no_valid_evaluations'
            print("WARNING: KTLC has no valid evaluation; the default parameter will be used.")
        else:
            selected_parameters = best_valid_parameters
            if best_valid_error_ev <= tolerance_ev:
                state = 'converged'
                selection_reason = 'tolerance_reached'
                print("Successful parameter optimization!")
                for i, res in enumerate(optimizer.res):
                    print(f"Iteration {i}: {res}")
            else:
                state = 'best_evaluated'
                selection_reason = 'evaluation_limit'
                print("WARNING: KTLC did not converge; using the best valid evaluated mu.")
            print(f"Selected mu is {selected_parameters[0]} "
                  f"with Koopmans error {best_valid_error_ev} eV.")

        self.para_functional = selected_parameters.copy()
        self.functional_parameter_optimization = {
            'state': state,
            'selected_mu': selected_parameters[0] if selected_parameters else None,
            'error_ev': best_valid_error_ev if best_valid_parameters is not None else None,
            'tolerance_ev': tolerance_ev,
            'valid_evaluations': valid_evaluations,
            'selection_reason': selection_reason,
        }
        return selected_parameters

    @calculation_cache.workflow('gaussian')
    def run_gaussian(self):
        """Run the complete workflow, enforcing ``timejob`` when configured."""
        from qcforever.conformer_search.settings import parse_conformer_options
        self._conformer_settings = parse_conformer_options(getattr(self, 'value', ''), getattr(self, 'conformer_config', None))
        self._conformer_summary = None
        self.functional_parameter_optimization = None
        original_directory = os.getcwd()
        self._active_job_directory = None
        try:
            with job_timeout.overall_timeout(self.timejob):
                return self._run_gaussian()
        except job_timeout.QCforeverTimeoutError as exc:
            print(exc)
            return {"log": "timeout", "error": str(exc)}
        finally:
            os.chdir(original_directory)
            if self._active_job_directory is not None:
                job_cleanup.cleanup_gaussian(
                    self._active_job_directory,
                    preserve_pickle=self.pklsave,
                    preserve_conformers=self._conformer_settings is not None,
                    preserve_xyz=not self.restart,
                )
            self._active_job_directory = None

    def _run_gaussian(self):
        infilename = self.in_file
        stable_from = Path(self.stable_from).expanduser().resolve() if self.stable_from else None
        option_line = self.value    
        from qcforever.conformer_search.settings import calculation_tokens
        options = calculation_tokens(option_line)
        job_eachState = []
        job_eachState.append({})
        option_dict = {}
        TargetStates = [0]
        TargetSpinMulti = []
        TargetTotalCharge = []
        p = Path(infilename)
        JobName, JobNameSuffix = p.stem, p.suffix
        JobNameSuffix = JobNameSuffix.lower()
        if self.basis == 'fchk' and JobNameSuffix not in ('.fchk', '.fch'):
            raise ValueError('basis=FCHK requires a Gaussian formatted checkpoint')
        if self.functional == 'ktlc-blyp-pred' and JobNameSuffix not in {'.sdf', '.xyz'}:
            raise ValueError('KTLC-BLYP-Pred requires SDF or XYZ molecular input')
        ReadFrom = None
        # Initialization of molecular coordinate
        atm, X, Y, Z = [], [], [], []
        Bondpair1 = []
        Bondpair2 = []
        # Classify input
        if JobNameSuffix == ".sdf":
            ReadFrom = 'sdf' 
            atm, X, Y, Z, TotalCharge, SpinMulti, Bondpair1, Bondpair2 = read_mol_file.read_sdf(infilename)
        elif JobNameSuffix == ".xyz":
            ReadFrom = 'xyz' 
            atm, X, Y, Z, TotalCharge, SpinMulti = read_mol_file.read_xyz(infilename)
        elif JobNameSuffix == ".chk":
            ReadFrom = 'chk' 
        elif JobNameSuffix in ('.fchk', '.fch'):
            TotalCharge, SpinMulti = gaussian_run.fchk2chk.ChargeSpin_fchk(infilename)
            ReadFrom = 'chk' 
        else:
            print("Invalid input file")
            exit()

        # Specify the electronic structure (spin multiplicity and charge) of the target state
        Is_ChargeSpec = False
        Is_SpinMultiSpec = False
        if np.isnan(self.SpecTotalCharge) !=  True:
            Is_ChargeSpec = True
            TotalCharge = self.SpecTotalCharge
        if np.isnan(self.SpecSpinMulti) != True:
            Is_SpinMultiSpec = True
            SpinMulti = self.SpecSpinMulti

        if self.functional == 'ktlc-blyp-pred':
            from qcforever.gaussian_run.predict_ktlc_mu import predict_ktlc_mu
            self.functional_parameter_prediction = predict_ktlc_mu(
                infilename, TotalCharge, SpinMulti, self.basis, self.solvent)
            self.para_functional = [gaussian_run.Gen_IOPline_functional.rounding_num(
                self.functional_parameter_prediction['mu_mean'])]
            self.functional = 'lc-blyp'
            print('KTLC-BLYP-Pred:', self.functional_parameter_prediction)

        TargetTotalCharge.append(TotalCharge)
        TargetSpinMulti.append(SpinMulti)

        # Setting options
        optoption = ''
        options, excited_specs = excited_options(options, 'gaussian')
        for i in range(len(options)):
            option = options[i]
            if re.match(r'^opt(?:=|$)',  option.lower()):
            #if 'optgeom' in option.lower():
                option_dict['opt'] = True
                job_eachState[0]['opt'] = True
                if '=' in option:
                    optoption = option.split('=', 1)[1].strip()
                    if optoption.startswith('(') and optoption.endswith(')'):
                        optoption = optoption[1:-1]
                    print(f'Use {optoption} for geometry optimization in the ground state.')
                else:
                    optoption = ''
            elif 'optconf' in option.lower():
                option_dict['optconf'] = True
                if '=' in option:
                    in_target = option.split('=')
                    optconfoption = in_target[-1]
                else:
                    optconfoption = 'pm6'
            elif 'optspin' in option.lower():
                option_dict['optspin'] = True
                option_dict['energy'] = True
                job_eachState[0]['optspin'] = True
                job_eachState[0]['energy'] = True
            elif re.match(r'^freq(?:=|$)', option.lower()):
                option_dict['freq'] = True
                job_eachState[0]['freq'] = True
                if '=' in option:
                    reference_paths = option.split('=', 1)[1].split(',')
                    if not 1 <= len(reference_paths) <= 3 or any(
                        not path for path in reference_paths
                    ):
                        raise ValueError(
                            "freq accepts IR[,Raman[,NMR]] reference files"
                        )
                    for spectrum_type, path in zip(
                        ("IR", "Raman", "NMR"), reference_paths
                    ):
                        self.ref_spectrum_paths[spectrum_type] = str(
                            Path(path).expanduser().resolve()
                        )
                    if len(reference_paths) == 3:
                        option_dict['nmr'] = True
                        job_eachState[0]['nmr'] = True
            elif option.lower() == 'polar':
                option_dict['polar'] =  True
                job_eachState[0]['polar'] = True
            elif re.match(r'^nmr(?:=|$)', option.lower()):
                option_dict['nmr'] = True
                job_eachState[0]['nmr'] = True
                if '=' in option:
                    reference_path = option.split('=', 1)[1]
                    if not reference_path:
                        raise ValueError("nmr requires a reference file after '='")
                    self.ref_spectrum_paths["NMR"] = str(
                        Path(reference_path).expanduser().resolve()
                    )
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
            elif re.match(r'^fluor(?:=|$)', option.lower()):
                option_dict['uv'] = True
                option_dict['fluor'] = True
                job_eachState.append({})
                job_eachState[0]['uv'] = True
                job_eachState[len(job_eachState)-1]['fluor'] = True
                job_eachState[-1]['excited_optoption'] = excited_specs['fluor']['fc']
                TargetStates.append(excited_specs['fluor']['state'])
                TargetSpinMulti.append(SpinMulti)
                TargetTotalCharge.append(TotalCharge)
            elif option.lower() == 'nac':
                option_dict['nac'] = True
                fluor_job = next((job for job in job_eachState if 'fluor' in job), None)
                if fluor_job is None:
                    raise ValueError("'nac' requires 'fluor' (e.g. fluor=1 nac).")
                fluor_job['nac'] = True
            elif re.match(r'^tadf(?:=|$)', option.lower()):
                option_dict['uv'] = True
                job_eachState[0]['uv'] = True
                if 'fluor' not in option_dict:
                    option_dict['fluor'] = True
                    job_eachState.append({'fluor': True, 'excited_optoption': ''})
                    TargetStates.append(1)
                    TargetSpinMulti.append(SpinMulti)
                    TargetTotalCharge.append(TotalCharge)
                option_dict['tadf'] = True
                job_eachState.append({'tadf': True, 'excited_optoption': excited_specs['tadf']['fc']})
                TargetStates.append(excited_specs['tadf']['state'])
                TargetSpinMulti.append(SpinMulti+2)
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
            elif option.partition('=')[0].lower() == 'cden':
                job_eachState[0]['cden'] = True
                option_dict['cden'] = True
                if cden_fukui(option):
                    job_eachState[0].update(fukui=True, vip=True, vea=True)
                    option_dict.update(fukui=True, vip=True, vea=True)
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
            elif option.lower() == 'stable':
                option_dict['stable'] = True
            else:
                print('Invalid option: ', option)

        self.mem = check_resource.respec_memory(
            spec_memory=self.mem,
            atoms=atm,
            basis=self.basis,
            total_charge=TotalCharge,
            nproc=self.nproc,
            options=option_dict.keys(),
            calculate_hessian='stable' in option_dict or requires_hessian(
                optoption, *(spec['fc'] for spec in excited_specs.values())),
        )

        # Make work directory and move to the directory
        pwd = os.getcwd()
        print(f'Current dir: {pwd}')
        calculation_cache.prepare_directory(JobName, infilename)
        self._active_job_directory = Path(pwd) / JobName
        if ReadFrom == 'chk':
            inchkfile = JobName + '.chk'
            if JobNameSuffix in ('.fchk', '.fch'):
                gaussian_run.fchk2chk.convert(infilename, Path(JobName) / inchkfile)
            else:
                shutil.copy2(infilename, Path(JobName) / inchkfile)
        os.chdir(JobName)
        print(f'Calculation dir: {os.getcwd()}')

        # Initialization of otuput dictionary.
        output_dic = []

        TDstate_info = []
        #print(job_eachState)
        #print(TargetTotalCharge)
        #print(TargetSpinMulti)

        job_state = ""

        # Setting the KTLC functionals
        functional_para_opt = False
        if self.functional == 'ktlc-blyp-bo' or  self.functional == 'ktlc-wpbe-bo' or self.functional == 'ktlc-whpbe-bo':
            functional_para_opt = True
            if self.functional == 'ktlc-blyp-bo':
                self.functional = 'lc-blyp' 
            elif self.functional == 'ktlc-wpbe-bo':
                self.functional = 'lc-wpbe' 
            elif self.functional == 'ktlc-whpbe-bo':
                self.functional = 'lc-whpbe' 
            else:
                print('Your specified functional is not avaiable!')
                print('Instead, we set KTLC-BLYP-BO.')
                self.functional = 'lc-blyp' 

        #Scaning spin state
        if 'optspin' in option_dict:
            ReSpinMulti  = self.SpinMulti_scan( 
                            JobName=JobName, 
                            ReadFrom=ReadFrom, 
                            GivenSpinMulti=SpinMulti, 
                            TotalCharge=TotalCharge, 
                            atm=atm, X=X, Y=Y, Z=Z)
            #output_dic.append({})
            #output_dic['spinmulti'] = ReSpinMulti
            if ReSpinMulti != SpinMulti:
                for i_sm in range(len(TargetSpinMulti)): 
                    TargetSpinMulti[i_sm] += ReSpinMulti - SpinMulti
                SpinMulti = ReSpinMulti
                print(f"Target States' Spin multiplicity: {TargetSpinMulti}")
            else:
                pass 

            if 'optconf' not in option_dict:
                ReadFrom = 'chk'
                    
        #Conformation search
        if 'optconf' in option_dict:
            print('Try to find stable conformation...')
            if ReadFrom == 'sdf' or ReadFrom == 'xyz':
                print('The input is ready for the conformation search...')
                original_file = str(Path("..") / infilename)
                #original_file = '../' + infilename
                print (type(original_file))
                try:
                    from qcforever.conformer_search.conformer_search import run_conformer_search
                    self._conformer_summary = run_conformer_search(
                        original_file, TotalCharge, SpinMulti, optconfoption, self.nproc, self.mem,
                        config=self._conformer_settings)
                except job_timeout.QCforeverTimeoutError:
                    raise
                except Exception as exc:
                    self._conformer_summary = {'state': 'failed', 'error': f'{type(exc).__name__}: {exc}',
                                               'details_directory': 'conformer_search'}
                    reconf = False
                    pass
            else:
                print('Conformation search is only possible when the input file is sdf or xyz.')
                reconf = False
                pass

            try:
                from qcforever.conformer_search.conformer_search import read_selected_structure
                # The SDF supplies geometry/bonds; keep the electronic state
                # passed to optconf (including overrides or an optspin result).
                atm, X, Y, Z, _, _, Bondpair1, Bondpair2 = read_selected_structure(self._conformer_summary)
                reconf = True
            except Exception as e:
                print('Conformation search is failed!')
                print(e)
                reconf = False
                pass

        #When ktlc-blyp-bo is specified as a functionl, try to optimize mu prameter with BO
        if functional_para_opt:
            self.para_functional = self.LC_para_BOopt(JobName, ReadFrom, atm, X, Y, Z, TotalCharge, SpinMulti)
            ReadFrom = 'chk'
        else:
            pass

        for i in range(len(TargetStates)):

            print (f'~~~~~~~~~~~~~Computing the state/species {i}~~~~~~~~~~~~~~~~~~~')

            output_dic.append({})

            job_thisstate = job_eachState[i]
            if TargetStates[i] == -1:
                #Get the index of the state whose OS is max. 
                targetstate = np.argmax(output_dic[0]["uv"][1]) + 1
            else:
                targetstate = TargetStates[i]
            SpinMulti = TargetSpinMulti[i]
            TotalCharge = TargetTotalCharge[i]
            if 'pka' in job_thisstate:
                JobNameState = JobName + f'_DeH'
            else:
                JobNameState = JobName + f'_State{targetstate}_{SpinMulti}_{TotalCharge}'

            '''#For scanning the spin state of the target 
                if 'optspin' in job_thisstate:
                    SpinMulti  = self.SpinMulti_scan(JobName=JobName, 
                                ReadFrom=ReadFrom, 
                                GivenSpinMulti=SpinMulti, 
                                TotalCharge=TotalCharge, 
                                atm=atm, X=X, Y=Y, Z=Z, 
                                targetstate=targetstate, 
                                TDstate_info=TDstate_info)
                    JobNameState = JobName + f'_State{targetstate}_{SpinMulti}_{TotalCharge}'
                    ReadFrom = 'chk'
                    output_dic[i]['spinmulti'] = SpinMulti
                else:
                    JobNameState = JobName + f'_State{targetstate}_{SpinMulti}_{TotalCharge}'
            '''

            Jobchk = JobNameState  + '.chk' 

            if ReadFrom == 'chk' and i == 0:
                inchkfile = JobName + '.chk'
                for f in glob.glob('./*.chk'):
                    File_name = f.replace('./','')
                    if File_name == inchkfile:
                        shutil.move(File_name, Jobchk) 
                        print (f'Get the information from {File_name}')
                    else:
                        print (f'{File_name} is not suitable for this job!')
            
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

            if i == 0 and 'stable' in option_dict:
                try:
                    point = self._stable_ground(JobNameState, TotalCharge, SpinMulti, ReadFrom,
                                                atm, X, Y, Z, optoption, stable_from)
                except job_timeout.QCforeverTimeoutError:
                    raise
                except (OSError, ValueError, RuntimeError) as exc:
                    output_dic[i].update(log_0='error', error_0=str(exc), failed_job=JobNameState,
                                         logfile=str(Path(JobNameState + '.log').resolve()))
                    break
                output_dic[i].update(point.properties, stable_0=point.minimum_frequency >= 0)
                if not output_dic[i]['stable_0']:
                    output_dic[i].update(log_0='error', error_0=point.properties.get('stable_error', NO_STABLE_MINIMUM))
                    break
                for key in ('opt', 'freq', 'energy'):
                    job_thisstate.pop(key, None)
                has_checkpoint = Path(Jobchk).is_file()
                if not job_thisstate:
                    scf_need = has_checkpoint
                    if has_checkpoint:
                        gaussian_run.chk2fchk.Get_chklist(1)
                    output_dic[i]['log_0'] = 'normal'
                    continue
                if has_checkpoint:
                    ReadFrom = 'chk'
                    atm, X, Y, Z = [], [], [], []
                else:
                    ReadFrom = 'xyz'
                    atm = point.atoms
                    X, Y, Z = point.coordinates.T

            check_logfile = JobName+'.log'
            scf_need = self.check_scf_need(job_thisstate)

            self.chain_job(JobNameState, scf_need, job_thisstate, TotalCharge, SpinMulti,
                                targetstate, ReadFrom,
                                element=atm, atomX=X, atomY=Y, atomZ=Z,
                                optoption=optoption if i == 0 else '', TDstate_info=TDstate_info)

            job_state = gaussian_run.Exe_Gaussian.exe_Gaussian(JobNameState, self.timexe, self.error)
            if job_state != 'normal':
                logfile = str(Path(JobNameState + '.log').resolve())
                output_dic[i].update({
                    f'log_{i}': job_state or 'abnormal',
                    f'error_{i}': f'Gaussian calculation ended with {job_state or "abnormal"}; see {logfile}',
                    'failed_job': JobNameState,
                    'logfile': logfile,
                })
                # Preserve incomplete checkpoints for diagnosis. Neither
                # formchk nor property extraction can repair a failed job.
                break
            # When the scf is performed, the obtained wavefunction is saved to chk file.
            # But for 'symm' and 'volume' that information is emply except for when the input is chk or fchk files.
            if scf_need or ReadFrom=='chk':
                gaussian_run.chk2fchk.Get_chklist(1)
            elif scf_need != True and ReadFrom != 'chk':
                for f in glob.glob('./*.chk'):
                    #os.remove(os.path.join('.', f))
                    Path(f).unlink()

            #output_prop = self.Extract_values(JobName, option_dict, Bondpair1, Bondpair2)
            #output_dic[i].update(output_prop)
            try:
                output_prop = self.Extract_values(JobNameState, job_thisstate, Bondpair1, Bondpair2)
                output_dic[i].update(output_prop)
                output_dic[i][f'log_{i}'] = job_state
            #except Exception as e:
            #    print(e)
            #    job_state = f"{e}"
            #    output_dic[i][f'log_{i}'] = job_state
            #    break
            except Exception as e:
                print(e)
                output_dic[i][f"log_{i}"] = "error"
                output_dic[i][f"error_{i}"] = (
                    f"{type(e).__name__}: {e}"
                )
                break

            # for pka computation
            if 'pka' in job_thisstate:
                #print("pka: ", output_dic_pka)
                E_dH = output_dic[i]["Energy"][0]
                output_dic[i]["pka"] = (E_dH - E_pH)*Eh2kJmol

            if 'fluor' in job_thisstate :
                if job_state == 'normal':
                    output_dic[i]["fluorescence"] = 'Likely'
                else:
                    output_dic[i]["fluorescence"] = 'Unlikely'


            if 'tadf' in job_thisstate:
                if job_state == 'normal':
                    output_dic[i]["phosphorescence"] = 'Likely'
                    TADF_Eng = 0.0
                    S_Eext = output_dic[i-1]["MinEtarget"]
                    T_Eext = output_dic[i]["T_Min"]
                    TADF_Eng = S_Eext - T_Eext
                    output_dic[i]["Delta(S-T)"] = TADF_Eng
                else:
                    output_dic[i]["phosphorescence"] = 'Unlikely'

        output_sum = {}
        for i, j in enumerate(output_dic):
            #print(i, j)
            for k in output_dic[i].keys():
                output_sum.setdefault(k, output_dic[i][k])
            #for k in range(len(option_dict)):

        if 'optspin' in option_dict:
            output_sum['spinmulti'] = ReSpinMulti

        if 'optconf' in option_dict:
            output_sum['optconf'] = reconf
            output_sum['conformer_search'] = self._conformer_summary or {'state': 'failed', 'error': 'Unsupported input'}

        #If functional parameter is not default ones, add them to the output dictionary.
        if self.para_functional != []:
            output_sum['functional_param'] = self.para_functional
        if self.functional_parameter_optimization is not None:
            output_sum['functional_parameter_optimization'] = self.functional_parameter_optimization
        if self.functional_parameter_prediction is not None:
            output_sum['functional_parameter_prediction'] = self.functional_parameter_prediction
    
        #Decide the log state        
        keylist = list(output_sum.keys())

        messages = {
            'error': 'Value extraction failed! ',
            'abnormal': 'Gaussian computation failed! ',
            'timeout': 'Not enough time for Gaussian calculations! ',
            'transition': 'Optimization may not succeed! ',
            'stagnation': 'Optimization may not succeed! ',
            'divergence': 'Optimization may not succeed! ',
            'oscillation': 'Optimization may not succeed! ',
            '': 'Wrong point is not identified! ',
        }

        logs = []

        #for key in keylist:
        #    if not re.match(r'log', key):
        #        continue
        #
        #    status = output_sum[key]
        #
        #    if status in messages:
        #        logs.append(messages[status])

        for key in keylist:
            if not re.fullmatch(r"log_\d+", key):
                continue

            status = output_sum[key]

            if status == "normal":
                continue

            logs.append(
                messages.get(
                    status,
                    f"Unexpected job state: {status}! "
                )
            )

        output_sum['log'] = ''.join(logs) if logs else 'normal'


        #if output_sum['log'] == 'normal' and 'fluor' in option_dict and 'opt' in option_dict and 'energy' in option_dict:
        #    print(output_sum['MinEtarget'])
        #    print(output_sum['Energy'][0])
        #    output_sum['E0-0'] = Eh2eV*(output_sum['MinEtarget']-output_sum['Energy'][0])
        #else:
        #    pass

        needs_e00 = (
            "fluor" in option_dict
            and "opt" in option_dict
            and "energy" in option_dict
        )

        if output_sum["log"] == "normal" and needs_e00:
            missing = [
                key for key in ("MinEtarget", "Energy")
                if key not in output_sum
            ]

            if missing:
                output_sum["log"] = (
                    "Value extraction failed! "
                    f"Missing result: {', '.join(missing)}"
                )
            elif not output_sum["Energy"]:
                output_sum["log"] = (
                    "Value extraction failed! Energy is empty."
                )
            else:
                output_sum["E0-0"] = Eh2eV * (
                    output_sum["MinEtarget"]
                    - output_sum["Energy"][0]
                )

        #Save as pickle
        if self.pklsave:
            result_dic = JobName+".pkl"
            with open(result_dic, 'wb') as f:
                pickle.dump(output_sum, f)

        # Convert fchk to xyz 
        if self.restart == False and scf_need == True and output_sum['log'] == 'normal':
            try:
                gaussian_run.Get_MolCoordinate_fchk.Get_fchklist2xyz()
            except:
                print('There is error in fchk')
                pass

        os.chdir(pwd)

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
