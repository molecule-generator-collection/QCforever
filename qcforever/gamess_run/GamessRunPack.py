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

        infilename = f"{jobname}.log"

        output = {}
        lines = gamess_run.read_log.read_log(infilename)

        if is_energy_specified:
            output['energy'] = gamess_run.read_log.getEnergy(lines)

        if is_homolumo_specified:
            num_occu_alpha, num_occu_beta  = gamess_run.read_log.getNumberElectron(lines)
            bb = gamess_run.read_log.getBlock(lines,"MOLECULAR ORBITALS")
            if bb == []:
                bb = gamess_run.read_log.getBlock(lines,"EIGENVECTORS")
                alpha_values = gamess_run.read_log.getMO_single(bb[-2])
                beta_values = gamess_run.read_log.getMO_single(bb[-1])
            else:
                alpha_values, beta_values = gamess_run.read_log.getMO_set(bb[-1])
            alpha_gap, beta_gap = gamess_run.read_log.gethomolumogap(alpha_values, beta_values, num_occu_alpha, num_occu_beta)
            output['homolumo'] = [alpha_gap, beta_gap]

        if is_dipole_specified:
            dd = gamess_run.read_log.getBlock(lines,"ELECTROSTATIC MOMENTS")    
            dval = gamess_run.read_log.getDipoleMoment(dd[-1])
            output['dipole'] = dval

        if is_uv_specified:
            exjobname = jobname+'_TD'
            infilename_ex = f"{exjobname}.log"
            exlines = gamess_run.read_log.read_log(infilename_ex)
            wavel, os = gamess_run.read_log.getTDDFT(exlines)

            output['uv'] = [wavel, os]
        
        lines = [] 
        exlines = []

        return output

    def make_input(
        self, run_type, 
        TotalCharge, SpinMulti, 
        GamInputName, Mol_atom=[], X=[], Y=[], Z=[], 
        TDDFT=False, datfile=None):

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

#make input
        with open(GamInputName ,'w') as ofile:
            line_input = f' $CONTRL SCFTYP=UHF RUNTYP={run_type} DFTTYP={self.functional}'
            if TDDFT:
                line_input += ' TDDFT=EXCITE \n' 
            else: 
                line_input += '\n'
            line_input += f'  COORD=UNIQUE NZVAR=0 ICHARG={TotalCharge} MULT={SpinMulti} $END\n'
            if TDDFT:
                line_input += ' $TDDFT NSTATE=10 $END\n' 
            line_input += ' $SCF damp=.TRUE. $END\n'
            time_limit = self.timeexe/60
            line_input += f' $SYSTEM TIMLIM={time_limit} MWORDS={mem_words} $END\n'
            if run_type == 'OPTIMIZE':
                line_input += ' $STATPT OPTTOL=1.0E-5 NSTEP=100 $END\n'
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
        
            GSdatfile = jobname+'.dat'
            exjobname = jobname+'_TD'
            GamInputName = exjobname +'.inp'    

            self.make_input(run_type, TotalCharge, SpinMulti, GamInputName, TDDFT=True, datfile=GSdatfile)

            job_state = gamess_run.Exe_Gamess.exe_Gamess(exjobname, self.gamessversion, self.nproc)

        try:
            output_dic = self.Extract_values(jobname, option_dict)
        except Exception as e: 
            job_state = "error"
            print(e)
            pass

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


