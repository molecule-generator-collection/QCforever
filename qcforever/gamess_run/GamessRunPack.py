import os
import re
import sys
import math
import shutil
import pickle

import numpy as np

from pathlib import Path

from qcforever import gamess_run
from qcforever.util.quadrupole import quadrupole_from_buckingham, quadrupole_properties
from qcforever.util import (
    read_mol_file,
    check_resource,
    job_cleanup,
    job_timeout,
    Spectrum_similarity,
)


byte2words = 1/8
Eh2kJmol = 2625.5
Eh2eV = 27.211
kcalmol2Eh = 1.593601E-3


class GamessDFTRun:

    def __init__(self, 
                functional, 
                basis, 
                nproc, 
                value, 
                in_file, 
                solvent="0",
                error=0, 
                pklsave=False):
        #self.in_file = os.path.basename(in_file)
        self.in_file = Path(in_file).name
        self.functional = functional
        self.basis = basis
        self.nproc = check_resource.respec_cores(nproc)
        self.value = value
        self.solvent = solvent.upper()
        self.error = error
        self.pklsave = pklsave
        self.conformer_config = None
        self.gamessversion = '00'
        self.mem = ''
        self.timeexe = 60 * 60 * 80
        # Overall QCforever wall-clock limit in seconds.  None means unlimited.
        # ``timeexe`` remains the value written to GAMESS input files.
        self.timejob = None
        self.SpecTotalCharge = np.nan
        self.SpecSpinMulti = np.nan
#        self.ref_uv_path = ''
        self.ref_spectrum_paths = {}
        self.para_functional = []
        self.functional_parameter_optimization = None

    def conversion_memory(self, mem_byte):
        r = re.compile("([0-9]+)([a-zA-Z]+)")
        m = r.match(mem_byte)
        memory = int(m.group(1))
        unit = (m.group(2)).upper()

        if unit == 'GB':
            memory_MB = memory*(10**3)

        if unit == 'MB':
            memory_MB = memory

        return memory_MB/8

    def Extract_values(self, jobname, option_dict):

        opt = option_dict['opt'] if 'opt' in option_dict else False
        is_energy_specified = option_dict['energy'] if 'energy' in option_dict else False
        is_homolumo_specified = option_dict['homolumo'] if 'homolumo' in option_dict else False
        is_dipole_specified = option_dict['dipole'] if 'dipole' in option_dict else False 
        is_cden_specified = option_dict['cden'] if 'cden' in option_dict else False 
        is_uv_specified = option_dict['uv'] if 'uv' in option_dict else False 
        is_fluor_specified = option_dict['fluor'] if 'fluor' in option_dict else False 
        is_vip_specified = option_dict['vip'] if 'vip' in option_dict else False 
        is_vea_specified = option_dict['vea'] if 'vea' in option_dict else False 
        is_aip_specified = option_dict['aip'] if 'aip' in option_dict else False 
        is_aea_specified = option_dict['aea'] if 'aea' in option_dict else False 
        is_freq_specified = option_dict['freq'] if 'freq' in option_dict else False 

        infilename = f"{jobname}.log"

        output = {}
        parselog = gamess_run.parse_log.parse_log(infilename)
        job_scfstate = parselog.Check_SCF()

        if is_energy_specified or is_freq_specified:
            if job_scfstate == True:
                output['energy'] = parselog.getEnergy()
            else:
                output['energy'] = ''
                output['log'] = job_scfstate

        if is_homolumo_specified:
            #To get the number of electrons.
            num_occu_alpha, num_occu_beta  = parselog.getNumberElectron()
            #To get the block for molecular orbitals
            bb = parselog.getBlock("MOLECULAR ORBITALS")
            if bb == []:
                bb = parselog.getBlock("EIGENVECTORS")
                alpha_values = parselog.getMO_single(bb[-2])
                beta_values = parselog.getMO_single(bb[-1])
            else:
                alpha_values, beta_values = parselog.getMO_set(bb[-1])

            alpha_gap, beta_gap = parselog.gethomolumogap(alpha_values, beta_values, num_occu_alpha, num_occu_beta)
            output['homolumo'] = [alpha_gap, beta_gap]
            output['Alpha_MO'] = [alpha_values[:num_occu_alpha],alpha_values[num_occu_alpha:]]
            output['Beta_MO'] = [beta_values[:num_occu_beta],alpha_values[num_occu_beta:]]

        if is_dipole_specified:
            if self.solvent == '0':
                dd = parselog.getBlock("ELECTROSTATIC MOMENTS")    
            else:
                dd = parselog.getBlock("SYSTEM ELECTROSTATIC MOMENTS")    
            dval = parselog.getDipoleMoment(dd[-1])
            output['dipole'] = dval
            try:
                moments = parselog.getQuadrupoleMoment(dd[-1])
                output.update(quadrupole_from_buckingham(
                    **moments, source="GAMESS", origin="center_of_mass"))
            except ValueError as error:
                output.update(quadrupole_properties(
                    source="GAMESS", origin="center_of_mass", source_convention="Buckingham",
                    status="unavailable", reason=str(error)))

        if is_cden_specified:
            output['cden'] = parselog.getChargeSpin()

        if is_uv_specified:
            exjobname = jobname + '_TD'
            infilename_ex = f"{exjobname}.log"
            parse_exlog = gamess_run.parse_log.parse_log(infilename_ex)
            wavel, os = parse_exlog.getTDDFT()
            output['uv'] = [wavel, os]

        if is_fluor_specified:
            fluorjobname = jobname + '_TDopt'
            infilename_fluor = f"{fluorjobname}.log"
            parse_fluorlog = gamess_run.parse_log.parse_log(infilename_fluor)
            wavel, os = parse_fluorlog.getTDDFT()
            output['fluor'] = [wavel, os]

        if is_vip_specified:
            vipjobname = jobname + '_VIP'
            infilename_vip = f"{vipjobname}.log"
            parse_viplog = gamess_run.parse_log.parse_log(infilename_vip)
            vip_scfstate = parse_viplog.Check_SCF()
            if vip_scfstate == True:
                vip_E = parse_viplog.getEnergy()
                output['vip'] = [Eh2eV * (vip_E[0] - output['energy'][0]), vip_E[1]]
            else:
                output['vip'] = ''
                output['log'] = vip_scfstate

        if is_vea_specified:
            veajobname = jobname + '_VEA'
            infilename_vea = f"{veajobname}.log"
            parse_vealog = gamess_run.parse_log.parse_log(infilename_vea)
            vea_scfstate = parse_vealog.Check_SCF()
            if vea_scfstate == True:
                vea_E = parse_vealog.getEnergy()
                output['vea'] = [Eh2eV * (output['energy'][0] - vea_E[0]), vea_E[1]]
            else:
                output['vea'] = ''
                output['log'] = vea_scfstate

        if is_aip_specified:
            aipjobname = jobname + '_AIP'
            infilename_aip = f"{aipjobname}.log"
            parse_aiplog = gamess_run.parse_log.parse_log(infilename_aip)
            aip_scfstate = parse_aiplog.Check_SCF()
            if aip_scfstate == True:
                aip_E = parse_aiplog.getEnergy()
                output['aip'] = [Eh2eV * (aip_E[0] - output['energy'][0]), aip_E[1]]
            else:
                output['vip'] = ''
                output['log'] = aip_scfstate

        if is_aea_specified:
            aeajobname = jobname + '_AEA'
            infilename_aea = f"{aeajobname}.log"
            parse_aealog = gamess_run.parse_log.parse_log(infilename_aea)
            aea_scfstate = parse_aealog.Check_SCF()
            if aea_scfstate == True:
                aea_E = parse_aealog.getEnergy()
                output['aea'] = [Eh2eV * (output['energy'][0] - aea_E[0]), aea_E[1]]
            else:
                output['aea'] = ''
                output['log'] = vea_scfstate

        if is_freq_specified:
            freqjobname = jobname + '_raman'
            infilename_freq = f"{freqjobname}.log"
            parse_freqlog = gamess_run.parse_log.parse_log(infilename_freq)
            freq_scfstate = parse_freqlog.Check_VIB()
            if freq_scfstate == True:
                ff = parse_freqlog.getBlock("NORMAL COORDINATE ANALYSIS IN THE HARMONIC APPROXIMATION")
                output['freq'], output['IR'], output['Raman'] = parse_freqlog.getFreq(ff[0])
                for spectrum_type in ('IR', 'Raman'):
                    reference_path = self.ref_spectrum_paths.get(spectrum_type)
                    if reference_path:
                        output.setdefault('Spectrum_similarity', {})[spectrum_type] = (
                            Spectrum_similarity.compare_with_file(
                                reference_path,
                                output['freq'],
                                output[spectrum_type],
                                spectrum_type,
                            )
                        )
                tt = parse_freqlog.getBlock("THERMOCHEMISTRY AT ")
                E_0, U, H, G, Cv, Cp, S =  parse_freqlog.getThermo(tt[0])
                output['Ezp'] = output['energy'][0] + E_0[0]
                output['Et'] = output['energy'][0] + U[1]*kcalmol2Eh
                output['E_enth'] = output['energy'][0] + H[1]*kcalmol2Eh
                output['E_free'] = output['energy'][0] + G[1]*kcalmol2Eh
                output['Ei'] = U[1]
                output['Cv'] = Cv[1]
                output['Cp'] = Cp[1]
                output['Si'] = S[1]
            else:
                output['freq'] = []
                output['IR'] = []
                output['Raman'] = []
                output['log'] = freq_scfstate

        lines = [] 
        exlines = []

        return output

    def MakeSolventLine(self):
        s = ''
        try:
            float(self.solvent)
        except ValueError:
            print('Solvent effect is included by PCM')
            s += f' $PCM SOLVNT={self.solvent} $END\n'
        else:
            if self.solvent == '0':
                return s, s_solvent
            else:
                print('Solvent effect is included by PCM')
                s += f' $PCM RSOLV=1.38 EPS={self.solvent} $END\n'
        return s

    def make_input(
        self, run_type, 
        TotalCharge, SpinMulti, 
        GamInputName, Mol_atom=[], X=[], Y=[], Z=[], 
        TDDFT=False, target=[1, 1], datfile=None, solvent='0'):
        #target[0]=target state index, target[1]=spin multiplicity of the target state
        from qcforever.conformer_search.settings import calculation_tokens
        is_dipole = any(token.lower() == 'dipole' for token in calculation_tokens(self.value))
#setting for memory
        if self.mem == '':
            mem_words = 125000000/10**6 #1 GB on 64 bit PC as the unit of 10**6 words 
        else:
            mem_words = self.conversion_memory(self.mem)

#setting for functional
        lc_option = False
        functional = ''
        if re.match('LC-', self.functional):
            kind_functional = self.functional.split('-')
            functional = kind_functional[-1]
            lc_option = True
            if self.para_functional != []:
                mu_param = self.para_functional[0]
            else:
                mu_param = None
        else:
            functional = self.functional

#setting for basis set
        GBASIS, NGAUSS, NDFUNC = gamess_run.make_basis.basis_dissection(self.basis)

        if datfile != None and Mol_atom == []:
            datlines = gamess_run.read_dat.read_dat(datfile)
            try:
                cc = gamess_run.read_dat.getMolBlock(datlines,"COORDINATES OF SYMMETRY UNIQUE ATOMS")
            except:
                cc = gamess_run.read_dat.get_dataBlock(datlines, "DATA")
            orb_bb = gamess_run.read_dat.get_dataBlock(datlines, "VEC")
            Norb = gamess_run.read_dat.count_orbital(orb_bb)

        if run_type == 'RAMAN':
            hess_bb = gamess_run.read_dat.get_dataBlock(datlines, "HESS")

        if TDDFT and run_type=='OPTIMIZE' and target[1]==1:
            scftype = "RHF"
        elif TDDFT and solvent != '0':
            scftype = "RHF"
        else:
            scftype = "UHF"

        CPCM = ''
        if solvent != '0':
            self.solvent = solvent
            CPCM = self.MakeSolventLine()

#make input
        with open(GamInputName ,'w') as ofile:
            line_input = f' $CONTRL SCFTYP={scftype} RUNTYP={run_type} DFTTYP={functional} MAXIT=100'
            if TDDFT == True:
                line_input += ' TDDFT=EXCITE \n' 
            elif TDDFT == 'SPNFLP':
                line_input += ' TDDFT=SPNFLP \n' 
            else: 
                line_input += '\n'
            line_input += f'  COORD=UNIQUE NZVAR=0 NOSYM=1 ICHARG={TotalCharge} MULT={SpinMulti} $END\n'
            if TDDFT:
                line_input += f' $TDDFT NSTATE=10 IROOT={target[0]}'
                if scftype == "RHF":
                    line_input += f' MULT={target[1]} $END\n' 
                else:
                    line_input += f' $END\n' 
            if lc_option == True:
                line_input += f' $DFT LC=.TRUE.'
                if mu_param != None:
                    line_input += f' MU={mu_param}'
                else:
                    pass
                line_input += f' $END\n'
            # Solvent effect cannot be combined with the analytical vibrational analsysis
            line_input += CPCM # for solvent effect line
            # Request Cartesian quadrupoles at a molecular origin for all jobs,
            # including inputs regenerated from a .dat file after optimization.
            if is_dipole:
                line_input += ' $ELMOM IEMOM=2 WHERE=COMASS $END\n'
            line_input += ' $SCF damp=.TRUE.'
            if run_type == 'HESSIAN' or run_type == 'RAMAN':
                line_input += ' DIRSCF=.TRUE. $END\n'
            else:
                line_input += ' $END\n'
            time_limit = self.timeexe/60
            line_input += f' $SYSTEM TIMLIM={time_limit} MWORDS={mem_words} $END\n'
            if run_type == 'OPTIMIZE':
                line_input += ' $STATPT OPTTOL=1.0E-5 NSTEP=100 $END\n'
            if run_type == 'HESSIAN' or run_type == 'RAMAN':
                line_input +=' $FORCE METHOD=ANALYTIC $END\n' 
            line_input += f' $BASIS GBASIS={GBASIS} NGAUSS={NGAUSS} NDFUNC={NDFUNC} $END\n'
            if datfile == None:
                line_input += ' $GUESS GUESS=HUCKEL $END\n' 
            else:
                line_input += f' $GUESS GUESS=MOREAD NORB={Norb} $END\n'
            line_input += ' $DATA\n'
            line_input += f'{GamInputName} {self.functional}/{self.basis}\n'
            line_input += 'C1\n'
            if datfile == None and Mol_atom != []:
                for i in range(len(Mol_atom)):
                    AtomicNum = float(gamess_run.AtomInfo.AtomicNumElec(Mol_atom[i]))
                    line_input += f'{Mol_atom[i]}   {AtomicNum:.1f}   {X[i]:.10f}   {Y[i]:.10f}   {Z[i]:.10f}\n' 
                line_input += ' $END\n'
            else:
                for i in range(len(cc)):
                    if len(cc[i].split()) == 5:
                        line_input += cc[i] 
                line_input += ' $END\n'
                line_input += ' $VEC\n'
                for i in range(len(orb_bb)):
                    line_input += orb_bb[i]
                line_input += ' $END\n'
                if run_type == 'RAMAN':
                    line_input += ' $HESS\n'
                    for i in range(len(hess_bb)):
                        line_input += hess_bb[i]
                    line_input += ' $END\n'

            ofile.write(line_input)

    def LC_para_BOopt(self, jobname, Mol_atom, X, Y, Z, TotalCharge, SpinMulti):
        """Tune mu and retain the best successfully evaluated LC parameters."""
        tolerance_ev = 0.01
        best_valid_parameters = None
        best_valid_error_ev = math.inf
        valid_evaluations = 0

        def KTLC_BB(mu):
            nonlocal best_valid_parameters, best_valid_error_ev, valid_evaluations

            print("mu parameter:", mu)
            # Zero retains the existing meaning: evaluate the GAMESS default.
            if mu == 0.:
                self.para_functional = []
            else:
                self.para_functional = [float(mu)]
            evaluated_parameters = self.para_functional.copy()

            if not math.isfinite(mu):
                raise ValueError('The trial mu must be finite')

            # Calculate the neutral, cation, and anion at the same geometry.
            GamInputName = jobname + '.inp'
            run_type = 'ENERGY'
            self.make_input(run_type, TotalCharge, SpinMulti, GamInputName, Mol_atom, X, Y, Z, TDDFT=False, datfile=None, solvent=self.solvent)
            gamess_run.Exe_Gamess.exe_Gamess(jobname, self.gamessversion, self.nproc)

            GSdatfile = jobname + '.dat'
            vipjobname = jobname + '_VIP'
            GamInputName = vipjobname + '.inp'
            IPTotalCharge = TotalCharge + 1
            if SpinMulti == 1:
                IPSpinMulti = SpinMulti + 1
            else:
                IPSpinMulti = SpinMulti - 1
            self.make_input(run_type, IPTotalCharge, IPSpinMulti, GamInputName, datfile=GSdatfile, solvent=self.solvent)
            gamess_run.Exe_Gamess.exe_Gamess(vipjobname, self.gamessversion, self.nproc)

            veajobname = jobname + '_VEA'
            GamInputName = veajobname + '.inp'
            EATotalCharge = TotalCharge - 1
            if SpinMulti == 1:
                EASpinMulti = SpinMulti + 1
            else:
                EASpinMulti = SpinMulti - 1
            self.make_input(run_type, EATotalCharge, EASpinMulti, GamInputName, datfile=GSdatfile, solvent=self.solvent)
            gamess_run.Exe_Gamess.exe_Gamess(veajobname, self.gamessversion, self.nproc)

            try:
                option_dict_lcparacheck = {'energy': True, 'homolumo': True, 'vip': True, 'vea': True}
                outdic = self.Extract_values(jobname, option_dict_lcparacheck)
                vip_ev = float(outdic['vip'][0])
                alpha_homo = float(outdic['Alpha_MO'][0][-1])
                beta_homo = float(outdic['Beta_MO'][0][-1])
                if not all(math.isfinite(value) for value in (vip_ev, alpha_homo, beta_homo)):
                    raise ValueError('VIP and both HOMO energies must be finite')
                Delta_IP = vip_ev + max(alpha_homo, beta_homo)*Eh2eV
                Delta_EA = outdic['vea'][0] + min(outdic['Alpha_MO'][1][0],outdic['Beta_MO'][1][0])*Eh2eV
                square_diff_satkoopmans = -Delta_IP**2
                error_ev = abs(float(Delta_IP))
                if not math.isfinite(error_ev) or not math.isfinite(square_diff_satkoopmans):
                    raise ValueError('The Koopmans error and objective must be finite')
            except job_timeout.QCforeverTimeoutError:
                raise
            except Exception as exc:
                print(f'Warning: mu={mu} evaluation failed: {exc}')
                # Preserve the BO failure score, but never select it as a result.
                return -10

            valid_evaluations += 1
            if error_ev < best_valid_error_ev:
                best_valid_parameters = evaluated_parameters
                best_valid_error_ev = error_ev
            return square_diff_satkoopmans

        from bayes_opt import BayesianOptimization
        from bayes_opt import UtilityFunction

        optimizer = BayesianOptimization(
            f=KTLC_BB,
            pbounds={'mu': (0, 1)},
            verbose=1,
            random_state=1,
        )

        init_mu = [0.15, 0.35, 0.65]
        for mu in init_mu:
            optimizer.set_bounds(new_bounds={'mu': (mu, mu)})
            optimizer.maximize(
                init_points=1,
                n_iter=0,
            )
        optimizer.set_bounds(new_bounds={'mu': (0, 1)})

        if best_valid_error_ev > tolerance_ev:
            i = 3
            acquisition_function = UtilityFunction(kind='ei', xi=1e-4)
            while best_valid_error_ev > tolerance_ev:
                optimizer.maximize(
                    init_points=0,
                    n_iter=1,
                    acquisition_function=acquisition_function,
                )
                i += 1
                if i > 51:
                    break

        for i, res in enumerate(optimizer.res):
            print(f'Iteration {i}: {res}')

        if best_valid_parameters is None:
            state = 'default_no_valid'
            selected_parameters = []
            selected_error_ev = None
            selection_reason = 'no_valid_evaluations'
            print('Warning: no valid mu evaluation was obtained; the default parameter will be used.')
        else:
            selected_parameters = best_valid_parameters
            selected_error_ev = best_valid_error_ev
            if best_valid_error_ev <= tolerance_ev:
                state = 'converged'
                selection_reason = 'tolerance_reached'
                print('Successful parameter optimization!')
            else:
                state = 'best_evaluated'
                selection_reason = 'evaluation_limit'
                print('Warning: parameter optimization did not converge; the best valid evaluated parameter will be used.')
            selected_label = selected_parameters[0] if selected_parameters else 'default (evaluated at mu=0)'
            print(f'Selected mu is {selected_label} with error {selected_error_ev} eV.')
        self.para_functional = selected_parameters.copy()
        self.functional_parameter_optimization = {
            'state': state,
            'selected_mu': selected_parameters[0] if selected_parameters else None,
            'error_ev': selected_error_ev,
            'tolerance_ev': tolerance_ev,
            'valid_evaluations': valid_evaluations,
            'selection_reason': selection_reason,
        }
        return selected_parameters


    def run_gamess(self):
        """Run the complete workflow, enforcing ``timejob`` when configured."""
        from qcforever.conformer_search.settings import parse_conformer_options
        self._conformer_settings = parse_conformer_options(getattr(self, 'value', ''), getattr(self, 'conformer_config', None))
        self._conformer_summary = None
        self.functional_parameter_optimization = None
        original_directory = os.getcwd()
        self._active_job_directory = None
        self._active_job_name = None
        try:
            with job_timeout.overall_timeout(self.timejob):
                return self._run_gamess()
        except job_timeout.QCforeverTimeoutError as exc:
            print(exc)
            return {"log": "timeout", "error": str(exc)}
        finally:
            os.chdir(original_directory)
            if self._active_job_directory is not None:
                job_cleanup.cleanup_gamess(
                    self._active_job_directory,
                    self._active_job_name,
                    preserve_pickle=self.pklsave,
                    preserve_conformers=self._conformer_settings is not None,
                )
            self._active_job_directory = None
            self._active_job_name = None

    def _run_gamess(self):
        infilename = self.in_file
        option_line = self.value    
        from qcforever.conformer_search.settings import calculation_tokens
        options = calculation_tokens(option_line)
        option_dict = {}
        #option_dict_Ex = np.zeros(19)  # not used
        #option_dict_pka = np.zeros(19) # not used
        targetstate = 1
        p = Path(infilename)
        jobname, jobnameSuffix  = p.stem, p.suffix 
        GamInputName = jobname +'.inp'    
        # File type of input?
        ReadFrom = '' 
        if jobnameSuffix == ".sdf":
            ReadFrom = "sdf" 
            Mol_atom, X, Y, Z, TotalCharge, SpinMulti, Bondpair1, Bondpair2 = read_mol_file.read_sdf(infilename)
        elif jobnameSuffix == ".xyz":
            ReadFrom = "xyz" 
            Mol_atom, X, Y, Z, TotalCharge, SpinMulti = read_mol_file.read_xyz(infilename)
            Bondpair1 = []
            Bondpair2 = []
        else:
            print("Invalid input file")
            exit()

        # Setting electronic structure
        if np.isnan(self.SpecTotalCharge) !=  True:
            TotalCharge = self.SpecTotalCharge
        if np.isnan(self.SpecSpinMulti) != True:
            SpinMulti = self.SpecSpinMulti

        # Setting options
        for i in range(len(options)):
            option = options[i]
            if option.lower() == 'opt':
                option_dict['opt'] = True
            elif 'optconf' in option.lower():
                option_dict['optconf'] = True
                if '=' in option:
                    in_target = option.split('=')
                    optconfoption = in_target[-1]
                else:
                    optconfoption = 'pm6'
            elif option.lower() == 'energy':
                option_dict['energy'] = True
            elif option.lower() == 'homolumo':
                option_dict['homolumo'] = True
            elif option.lower() == 'dipole':
                option_dict['dipole'] = True
            elif option.lower() == 'cden':
                option_dict['cden'] = True
            elif option.lower() == 'uv':
                option_dict['uv'] = True
            elif option.lower() == 'fluor':
                option_dict['uv'] = True
                option_dict['fluor'] = True
                if '=' in option:
                    in_target = option.split('=')
                    targetstate = int(in_target[-1])
            elif option.lower() == 'vip':
                option_dict['energy'] = True
                option_dict['vip'] = True
            elif option.lower() == 'vea':
                option_dict['energy'] = True
                option_dict['vea'] = True
            elif option.lower() == 'aip':
                option_dict['energy'] = True
                option_dict['vip'] = True
                option_dict['aip'] = True
            elif option.lower() == 'aea':
                option_dict['energy'] = True
                option_dict['vea'] = True
                option_dict['aea'] = True
            elif re.match(r'^freq(?:=|$)', option.lower()):
                option_dict['freq'] = True
                if '=' in option:
                    reference_paths = option.split('=', 1)[1].split(',')
                    if not 1 <= len(reference_paths) <= 2 or any(
                        not path for path in reference_paths
                    ):
                        raise ValueError("freq accepts IR[,Raman] reference files")
                    for spectrum_type, path in zip(('IR', 'Raman'), reference_paths):
                        self.ref_spectrum_paths[spectrum_type] = str(
                            Path(path).expanduser().resolve()
                        )
            else:
                print('invalid option: ', option)

        self.mem = check_resource.respec_memory(
            spec_memory=self.mem,
            atoms=Mol_atom,
            basis=self.basis,
            total_charge=TotalCharge,
            nproc=self.nproc,
            options=option_dict.keys(),
        )

        output_dic = {}
        # Make work directory and move to the directory
        pwd = os.getcwd()
        print(pwd)

        job_dir = Path(jobname)
        if job_dir.is_dir():
            shutil.rmtree(job_dir)

        job_dir.mkdir()
        self._active_job_directory = Path(pwd) / job_dir
        self._active_job_name = jobname
        os.chdir(job_dir)


        #Conformation search
        if 'optconf' in option_dict:
            print('Try to find stable conformation...')
            if ReadFrom == 'sdf' or ReadFrom == 'xyz':
                print('The input is ready for the conformation search...')
                original_file = str(Path('..') / infilename)
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
                # passed to optconf rather than its formal charge/radicals.
                Mol_atom, X, Y, Z, _, _, Bondpair1, Bondpair2 = read_selected_structure(self._conformer_summary)
                reconf = True
            except Exception as e:
                print('Conformation search is failed!')
                print(e)
                reconf = False
                pass

#setting for KTLC-functional
        if re.match('KTLC-', self.functional):
            base_functional = self.functional.split('-')
            self.functional = "LC-" + base_functional[1]
            self.para_functional = self.LC_para_BOopt(jobname, Mol_atom, X, Y, Z, TotalCharge, SpinMulti)

#setting of run type for the ground state
        if 'opt' in  option_dict:
            run_type = 'OPTIMIZE'
        else:
            run_type = 'ENERGY'

        self.make_input(run_type, TotalCharge, SpinMulti, GamInputName, Mol_atom, X, Y, Z, TDDFT=False, datfile=None, solvent=self.solvent)
        job_state = gamess_run.Exe_Gamess.exe_Gamess(jobname, self.gamessversion, self.nproc)

#for uv computation
        if 'uv' in  option_dict:
            run_type = 'ENERGY'
        
            GSdatfile = jobname + '.dat'
            exjobname = jobname + '_TD'
            GamInputName = exjobname +'.inp'    

            self.make_input(run_type, TotalCharge, SpinMulti, GamInputName, TDDFT=True, datfile=GSdatfile, solvent=self.solvent)
            job_state = gamess_run.Exe_Gamess.exe_Gamess(exjobname, self.gamessversion, self.nproc)

        if 'fluor' in  option_dict:
            run_type = 'OPTIMIZE'

            GSdatfile = jobname + '.dat'
            exoptjobname = jobname + '_TDopt'
            GamInputName = exoptjobname +'.inp'    

            self.make_input(run_type, TotalCharge, SpinMulti, GamInputName, TDDFT=True, target=[targetstate, SpinMulti], datfile=GSdatfile, solvent=self.solvent)
            job_state = gamess_run.Exe_Gamess.exe_Gamess(exoptjobname, self.gamessversion, self.nproc)

        if 'vip' in option_dict:
            run_type = 'ENERGY'

            GSdatfile = jobname + '.dat'
            vipjobname = jobname + '_VIP'
            GamInputName = vipjobname + '.inp'

            IPTotalCharge = TotalCharge + 1
        
            if SpinMulti == 1:
                IPSpinMulti = SpinMulti + 1
            else:
                IPSpinMulti = SpinMulti - 1

            self.make_input(run_type, IPTotalCharge, IPSpinMulti, GamInputName, datfile=GSdatfile, solvent=self.solvent)
            job_state = gamess_run.Exe_Gamess.exe_Gamess(vipjobname, self.gamessversion, self.nproc)

        if 'vea' in option_dict:
            run_type = 'ENERGY'

            GSdatfile = jobname + '.dat'
            veajobname = jobname + '_VEA'
            GamInputName = veajobname + '.inp'

            EATotalCharge = TotalCharge - 1
        
            if SpinMulti == 1:
                EASpinMulti = SpinMulti + 1
            else:
                EASpinMulti = SpinMulti - 1

            self.make_input(run_type, EATotalCharge, EASpinMulti, GamInputName, datfile=GSdatfile, solvent=self.solvent)
            job_state = gamess_run.Exe_Gamess.exe_Gamess(veajobname, self.gamessversion, self.nproc)

        if 'aip' in option_dict:
            run_type = 'OPTIMIZE'

            VIPdatfile = vipjobname + '.dat'
            aipjobname = jobname + '_AIP'
            GamInputName = aipjobname + '.inp'

            self.make_input(run_type, IPTotalCharge, IPSpinMulti, GamInputName, datfile=VIPdatfile, solvent=self.solvent)
            job_state = gamess_run.Exe_Gamess.exe_Gamess(aipjobname, self.gamessversion, self.nproc)

        if 'aea' in option_dict:
            run_type = 'OPTIMIZE'

            VEAdatfile = veajobname + '.dat'
            aeajobname = jobname + '_AEA'
            GamInputName = aeajobname + '.inp'

            self.make_input(run_type, EATotalCharge, EASpinMulti, GamInputName, datfile=VEAdatfile, solvent=self.solvent)
            job_state = gamess_run.Exe_Gamess.exe_Gamess(aeajobname, self.gamessversion, self.nproc)

        if 'freq' in  option_dict:
            run_type = 'HESSIAN'
        
            GSdatfile = jobname + '.dat'
            freqjobname = jobname + '_freq'
            GamInputName = freqjobname +'.inp'    

            #self.make_input(run_type, TotalCharge, SpinMulti, GamInputName, datfile=GSdatfile, solvent=self.solvent)
            self.make_input(run_type, TotalCharge, SpinMulti, GamInputName, datfile=GSdatfile, solvent='0')
            job_state = gamess_run.Exe_Gamess.exe_Gamess(freqjobname, self.gamessversion, self.nproc)

            run_type= 'RAMAN'

            freqdatfile = jobname + '_freq.dat'
            ramanjobname = jobname + '_raman'
            GamInputName = ramanjobname + '.inp'

            #self.make_input(run_type, TotalCharge, SpinMulti, GamInputName, datfile=freqdatfile, solvent=self.solvent)
            self.make_input(run_type, TotalCharge, SpinMulti, GamInputName, datfile=freqdatfile, solvent='0')
            job_state = gamess_run.Exe_Gamess.exe_Gamess(ramanjobname, self.gamessversion, self.nproc)

#        output_dic = self.Extract_values(jobname, option_dict)
        try:
            output_dic = self.Extract_values(jobname, option_dict)
        except Exception as e: 
            job_state = "error"
            print(e)
            pass

        if self.para_functional != []:
            output_dic['functional_param'] = self.para_functional
        if self.functional_parameter_optimization is not None:
            output_dic['functional_parameter_optimization'] = self.functional_parameter_optimization.copy()

        if 'optconf' in option_dict:
            output_dic['optconf'] = reconf
            output_dic['conformer_search'] = self._conformer_summary or {'state': 'failed', 'error': 'Unsupported input'}

        if 'log' not in output_dic:
            if job_state == '':
                output_dic["log"] = 'normal' 
            else:
                output_dic["log"] = job_state

        #Save as pickle
        if self.pklsave:
            result_dic = jobname+".pkl"
            with open(result_dic, 'wb') as f:
                pickle.dump(output_dic, f)

        os.chdir(pwd)

        return(output_dic)
            
				
if __name__ == '__main__':
    usage ='Usage; %s infile' % sys.argv[0]

    try:
        infilename = sys.argv[1]
    except:
        print (usage); sys.exit()

    test_sdf = GamessDFTRun('B3LYP', '3-21g*',8, 'opt',infilename,0)

    test_sdf.run_gamess()
