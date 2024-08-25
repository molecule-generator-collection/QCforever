import os
import re
import sys
import math
import shutil
import pickle

import numpy as np

from qcforever import gamess_run
from qcforever.util import read_mol_file, check_resource
from qcforever.laqa_fafoom import laqa_confopt_sdf


byte2words = 1/8
Eh2kJmol = 2625.5
Eh2eV = 27.211
kcalmol2Eh = 1.593601E-3


class GamessDFTRun:

    def __init__(self, functional, basis, nproc, value, in_file, error=0, pklsave=False):
        self.in_file = in_file
        self.functional = functional
        self.basis = basis
        self.nproc = check_resource.respec_cores(nproc)
        self.value = value
        self.error = error
        self.pklsave = pklsave
        self.gamessversion = '00'
        self.mem = ''
        self.timeexe = 60 * 60 * 80
        self.SpecTotalCharge = np.nan
        self.SpecSpinMulti = np.nan
#        self.ref_uv_path = ''
        self.para_functional = []

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
            dd = parselog.getBlock("ELECTROSTATIC MOMENTS")    
            dval = parselog.getDipoleMoment(dd[-1])
            output['dipole'] = dval

        if is_dipole_specified:
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
                output['vip'] = Eh2eV * (vip_E - output['energy'])
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
                output['vea'] = Eh2eV * (output['energy'] - vea_E)
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
                output['aip'] = Eh2eV * (aip_E - output['energy'])
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
                output['aea'] = Eh2eV * (output['energy'] - aea_E)
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
                tt = parse_freqlog.getBlock("THERMOCHEMISTRY AT ")
                E_0, U, H, G, Cv, Cp, S =  parse_freqlog.getThermo(tt[0])
                output['Ezp'] = output['energy'] + E_0[0]
                output['Et'] = output['energy'] + U[1]*kcalmol2Eh
                output['E_enth'] = output['energy'] + H[1]*kcalmol2Eh
                output['E_free'] = output['energy'] + G[1]*kcalmol2Eh
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

    def make_input(
        self, run_type, 
        TotalCharge, SpinMulti, 
        GamInputName, Mol_atom=[], X=[], Y=[], Z=[], 
        TDDFT=False, target=[1, 1], datfile=None):
        #target[0]=target state index, target[1]=spin multiplicity of the target state
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
        else:
            scftype = "UHF"

#make input
        with open(GamInputName ,'w') as ofile:
            line_input = f' $CONTRL SCFTYP={scftype} RUNTYP={run_type} DFTTYP={functional}'
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
        
        from bayes_opt import BayesianOptimization
        from bayes_opt import UtilityFunction


        def KTLC_BB(mu):
            self.para_functional = [mu]

            #Calculating the ground state
            GamInputName = jobname +'.inp'
            run_type = 'ENERGY'

            self.make_input(run_type, TotalCharge, SpinMulti, GamInputName, Mol_atom, X, Y, Z, TDDFT=False, datfile=None)
            job_state = gamess_run.Exe_Gamess.exe_Gamess(jobname, self.gamessversion, self.nproc)
            ###########################

            #vip computation
            run_type = 'ENERGY'

            GSdatfile = jobname + '.dat'
            vipjobname = jobname + '_VIP'
            GamInputName = vipjobname + '.inp'

            IPTotalCharge = TotalCharge + 1
        
            if SpinMulti == 1:
                IPSpinMulti = SpinMulti + 1
            else:
                IPSpinMulti = SpinMulti - 1

            self.make_input(run_type, IPTotalCharge, IPSpinMulti, GamInputName, datfile=GSdatfile)
            job_state = gamess_run.Exe_Gamess.exe_Gamess(vipjobname, self.gamessversion, self.nproc)
            ###########################

            #'vea' computation
            run_type = 'ENERGY'

            GSdatfile = jobname + '.dat'
            veajobname = jobname + '_VEA'
            GamInputName = veajobname + '.inp'

            EATotalCharge = TotalCharge - 1
        
            if SpinMulti == 1:
                EASpinMulti = SpinMulti + 1
            else:
                EASpinMulti = SpinMulti - 1

            self.make_input(run_type, EATotalCharge, EASpinMulti, GamInputName, datfile=GSdatfile)
            job_state = gamess_run.Exe_Gamess.exe_Gamess(veajobname, self.gamessversion, self.nproc)
            ###########################

            job_state = "normal"
            option_dict_lcparacheck= {'energy': True, 'homolumo': True, 'vip': True, 'vea': True}
            try:
                outdic = self.Extract_values(jobname, option_dict_lcparacheck)
            except Exception as e: 
                job_state = "error"
                print(e)
                pass
            ##For debug##########################
            #outdic = self.Extract_values(jobname, option_dict_lcparacheck)
            #print(outdic)
            #####################################

            square_diff_satkoopmans = -10
            if job_state == "normal":
                Delta_IP = outdic['vip'] + max(outdic['Alpha_MO'][0][-1],outdic['Beta_MO'][0][-1])*Eh2eV
                Delta_EA = outdic['vea'] + min(outdic['Alpha_MO'][1][0],outdic['Beta_MO'][1][0])*Eh2eV

                square_diff_satkoopmans = -Delta_IP**2 
            else:
                square_diff_satkoopmans = -10
                pass
            
            return square_diff_satkoopmans

        pbounds = {'mu': (0,1)}

        optimizer = BayesianOptimization( 
        f=KTLC_BB,
        pbounds=pbounds,
        verbose=1, 
        random_state=1,
        )

        #initial mu values 
        init_mu = [0.15, 0.35, 0.65]
   
        #initial try...
        for i in range(len(init_mu)):
            optimizer.set_bounds(new_bounds={"mu": (init_mu[i],init_mu[i])})
            optimizer.maximize(
                init_points=1,
                n_iter=0,
            )

        diff_koopmans = np.sqrt(abs(optimizer.max['target']))
        optimizer.set_bounds(new_bounds={"mu": (0,1)})

        #Extra try...
        if diff_koopmans <= 0.01:
            print("Successful parameter optimization!")
            pass
        else:
            i = 3 
            acquisition_function = UtilityFunction(kind="ei", xi=1e-4)
            while diff_koopmans > 0.01:
                optimizer.maximize(
                    init_points=0,
                    n_iter=1,
                    acquisition_function=acquisition_function,
                )
                i += 1
                diff_koopmans = np.sqrt(abs(optimizer.max['target']))
                if i > 51:
                    break

        if diff_koopmans > 0.01:
            print("The prameter optimization failed!")
            print("The default parameter will be used!")
            return []
        else:
            print("Successful parameter optimization!")

            for i, res in enumerate(optimizer.res):
                print(f"Iteration {i}: {res}")

            diff_koopmans = np.sqrt(abs(optimizer.max['target']))
            print(f"Optimized mu is {optimizer.max['params']['mu']} with the difference {diff_koopmans}.")

            return [optimizer.max['params']['mu']]


    def run_gamess(self):
        infilename = self.in_file
        option_line = self.value    
        options = option_line.split()
        option_dict = {}
        #option_dict_Ex = np.zeros(19)  # not used
        #option_dict_pka = np.zeros(19) # not used
        targetstate = 1
        PreGamInput = infilename.split('.')
        jobname = PreGamInput[0]
        GamInputName = jobname +'.inp'    
        # File type of input?
        ReadFrom = '' 
        if PreGamInput[1] == "sdf":
            ReadFrom = "sdf" 
            Mol_atom, X, Y, Z, TotalCharge, SpinMulti, Bondpair1, Bondpair2 = read_mol_file.read_sdf(infilename)
        elif PreGamInput[1] == "xyz":
            ReadFrom = "xyz" 
            Mol_atom, X, Y, Z, TotalCharge, SpinMulti = read_mol_file.read_xyz(infilename)
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
            elif option.lower() == 'freq':
                option_dict['freq'] = True
            else:
                print('invalid option: ', option)

        output_dic = {}
        # Make work directory and move to the directory
        pwd = os.getcwd()
        print(pwd)
        if os.path.isdir(PreGamInput[0]):
            shutil.rmtree(PreGamInput[0])
        os.mkdir(PreGamInput[0])
        #shutil.move(GamInputName, PreGamInput[0])
        os.chdir(PreGamInput[0])

        #Conformation search
        if 'optconf' in option_dict:
            print('Try to find stable conformation...')
            if ReadFrom == 'sdf':
                print('The input is a sdf file, OK...')
                original_sdf = '../' + infilename
                try:
                    laqa_confopt_sdf.LAQA_confopt_main(original_sdf, TotalCharge, SpinMulti, 
                                                        optconfoption, self.nproc, self.mem)
                except:
                    reconf = False
                    pass
                    
            else:
                print('Conformation search is only possible when the input file is sdf.')
                reconf = False
                pass

            try:
                Mol_atm, X, Y, Z, TotalCharge, SpinMulti, Bondpair1, Bondpair2 = read_mol_file.read_sdf("./optimized_structures.sdf")
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
            #try:
            #    self.para_functional = [self.LC_para_BOopt(jobname, Mol_atom, X, Y, Z, TotalCharge, SpinMulti)]
            #except: 
            #    self.para_functional = []

#setting of run type for the ground state
        if 'opt' in  option_dict:
            run_type = 'OPTIMIZE'
        else:
            run_type = 'ENERGY'

        self.make_input(run_type, TotalCharge, SpinMulti, GamInputName, Mol_atom, X, Y, Z, TDDFT=False, datfile=None)
        job_state = gamess_run.Exe_Gamess.exe_Gamess(jobname, self.gamessversion, self.nproc)

#for uv computation
        if 'uv' in  option_dict:
            run_type = 'ENERGY'
        
            GSdatfile = jobname + '.dat'
            exjobname = jobname + '_TD'
            GamInputName = exjobname +'.inp'    

            self.make_input(run_type, TotalCharge, SpinMulti, GamInputName, TDDFT=True, datfile=GSdatfile)
            job_state = gamess_run.Exe_Gamess.exe_Gamess(exjobname, self.gamessversion, self.nproc)

        if 'fluor' in  option_dict:
            run_type = 'OPTIMIZE'

            GSdatfile = jobname + '.dat'
            exoptjobname = jobname + '_TDopt'
            GamInputName = exoptjobname +'.inp'    

            self.make_input(run_type, TotalCharge, SpinMulti, GamInputName, TDDFT=True, target=[targetstate, SpinMulti], datfile=GSdatfile)
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

            self.make_input(run_type, IPTotalCharge, IPSpinMulti, GamInputName, datfile=GSdatfile)
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

            self.make_input(run_type, EATotalCharge, EASpinMulti, GamInputName, datfile=GSdatfile)
            job_state = gamess_run.Exe_Gamess.exe_Gamess(veajobname, self.gamessversion, self.nproc)

        if 'aip' in option_dict:
            run_type = 'OPTIMIZE'

            VIPdatfile = vipjobname + '.dat'
            aipjobname = jobname + '_AIP'
            GamInputName = aipjobname + '.inp'

            self.make_input(run_type, IPTotalCharge, IPSpinMulti, GamInputName, datfile=VIPdatfile)
            job_state = gamess_run.Exe_Gamess.exe_Gamess(aipjobname, self.gamessversion, self.nproc)

        if 'aea' in option_dict:
            run_type = 'OPTIMIZE'

            VEAdatfile = veajobname + '.dat'
            aeajobname = jobname + '_AEA'
            GamInputName = aeajobname + '.inp'

            self.make_input(run_type, EATotalCharge, EASpinMulti, GamInputName, datfile=VEAdatfile)
            job_state = gamess_run.Exe_Gamess.exe_Gamess(aeajobname, self.gamessversion, self.nproc)

        if 'freq' in  option_dict:
            run_type = 'HESSIAN'
        
            GSdatfile = jobname + '.dat'
            freqjobname = jobname + '_freq'
            GamInputName = freqjobname +'.inp'    

            self.make_input(run_type, TotalCharge, SpinMulti, GamInputName, datfile=GSdatfile)
            job_state = gamess_run.Exe_Gamess.exe_Gamess(freqjobname, self.gamessversion, self.nproc)

            run_type= 'RAMAN'

            freqdatfile = jobname + '_freq.dat'
            ramanjobname = jobname + '_raman'
            GamInputName = ramanjobname + '.inp'

            self.make_input(run_type, TotalCharge, SpinMulti, GamInputName, datfile=freqdatfile)
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

        if 'optconf' in option_dict:
            output_dic['optconf'] = reconf

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


