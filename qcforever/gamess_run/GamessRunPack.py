import os
import re
import sys
import math
import shutil
import numpy as np


from qcforever import gamess_run
from qcforever.util import read_mol_file


byte2words = 1/8
Eh2kJmol = 2625.5
Eh2eV = 27.211


class GamessDFTRun:

    def __init__(self, functional, basis, nproc, value, in_file, error=0):
        self.in_file = in_file
        self.functional = functional
        self.basis = basis
        self.nproc = nproc
        self.value = value
        self.error = error
        self.gamessversion = '00'
        self.mem = ''
        self.timeexe = 60 * 60 * 80
        self.SpecTotalCharge = np.nan
        self.SpecSpinMulti = np.nan
#        self.ref_uv_path = ''
#        self.para_functional = []

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

        if is_energy_specified:
            if job_scfstate == True:
                output['energy'] = parselog.getEnergy()
            else:
                output['energy'] = ''
                output['log'] = job_scfstate

        if is_homolumo_specified:
            num_occu_alpha, num_occu_beta  = parselog.getNumberElectron()
            bb = parselog.getBlock("MOLECULAR ORBITALS")
            if bb == []:
                bb = parselog.getBlock("EIGENVECTORS")
                alpha_values = parselog.getMO_single(bb[-2])
                beta_values = parselog.getMO_single(bb[-1])
            else:
                alpha_values, beta_values = parselog.getMO_set(bb[-1])
            alpha_gap, beta_gap = parselog.gethomolumogap(alpha_values, beta_values, num_occu_alpha, num_occu_beta)
            output['homolumo'] = [alpha_gap, beta_gap]

        if is_dipole_specified:
            dd = parselog.getBlock("ELECTROSTATIC MOMENTS")    
            dval = parselog.getDipoleMoment(dd[-1])
            output['dipole'] = dval

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
            line_input = f' $CONTRL SCFTYP={scftype} RUNTYP={run_type} DFTTYP={self.functional}'
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
        ReadFromsdf = 0 
        ReadFromxyz = 0 
        if PreGamInput[1] == "sdf":
            ReadFromsdf = 1 
            Mol_atom, X, Y, Z, TotalCharge, SpinMulti, Bondpair1, Bondpair2 = read_mol_file.read_sdf(infilename)
        elif PreGamInput[1] == "xyz":
            ReadFromxyz = 1
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
            elif option.lower() == 'energy':
                option_dict['energy'] = True
            elif option.lower() == 'homolumo':
                option_dict['homolumo'] = True
            elif option.lower() == 'dipole':
                option_dict['dipole'] = True
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

#setting for run type
        if 'opt' in  option_dict:
            run_type = 'OPTIMIZE'
        else:
            run_type = 'ENERGY'

        self.make_input(run_type, TotalCharge, SpinMulti, GamInputName, Mol_atom, X, Y, Z, TDDFT=False, datfile=None)

        output_dic = {}
        if os.path.isdir(PreGamInput[0]):
            shutil.rmtree(PreGamInput[0])
        os.mkdir(PreGamInput[0])
        shutil.move(GamInputName, PreGamInput[0])
        os.chdir(PreGamInput[0])
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

        try:
            output_dic = self.Extract_values(jobname, option_dict)
        except Exception as e: 
            job_state = "error"
            print(e)
            pass

        if 'log' not in output_dic:
            output_dic["log"] = job_state

        os.chdir("..")

        return(output_dic)
            
				
if __name__ == '__main__':
    usage ='Usage; %s infile' % sys.argv[0]

    try:
        infilename = sys.argv[1]
    except:
        print (usage); sys.exit()

    test_sdf = GamessDFTRun('B3LYP', '3-21g*',8, 'opt',infilename,0)

    test_sdf.run_gamess()


