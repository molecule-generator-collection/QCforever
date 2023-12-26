import re
import sys
import os
import numpy as np

bohr2cm = 5.2918e-9
Ac = 6.0221e+23

class parse_log:

    def __init__(self, outfile):
        with open (outfile) as output_line:
            self.lines = output_line.readlines()

    def SplitLinks(self):
        Links = []
        Link = ""
        for line in self.lines:
            if line.find("Initial command:")>0:
                Links.append(Link)
                Link = ""
            Link = Link + line
        Links.append(Link)
        return Links

    def extract_method_text(self, lines):
        text = '\n'.join(lines)
        #pattern = r'\s#.*?(?=\s-{35,70})'  
        pattern = r'\s#.*?(?=\s-{4,71})'  
        enclosed_text = re.findall(pattern, text, re.DOTALL)
        method_text = ''
        for line in enclosed_text:
            method_text += re.sub(r'\n\s', '', line.strip())
        return method_text

    def extract_functionalbasis(self, line):
        output = []
        # Extract the part between '#' and the first space
        part = line.split(" ")[0]
        # Remove '#r' or '#u' if present
        part = part.replace("#r", "").replace("#u", "")
        split_method = part.split('/')
        for i in range(len(split_method)):
            output.append(split_method[i])
        return output

    def extract_MolCoordlog(self, lines):
        text = '\n'.join(lines)
        pattern_block = r"Standard orientation:(?:\s+(-+\n\s+Center\s+Atomic\s+Atomic\s+Coordinates \(Angstroms\)\n\s+Number\s+Number\s+Type\s+X\s+Y\s+Z\n\s+(-+\n\s+\d+\s+\d+\s+\d+\s+[-.\d\s]+\n)+)(?=\s+-+\n))"
        Coord_blocks = re.findall(pattern_block, text)
        #print(Coord_blocks)

        for block in Coord_blocks:
            print ('Atomic coordinates are found.')
            #CentNum = []
            #atomicNum = []
            pattern = r'\s*\d+\s+\d+\s+\d+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+'
            coordinates_inf = re.findall(pattern, str(block[1]))
            coordinates_tmp = []
            for mcoord in coordinates_inf:
                coordinates_tmp.append(mcoord.strip())
            CentNum_AtomicN = [list(map(float, line.split()[:2])) for line in coordinates_tmp]
            coordinates = [list(map(float, line.split()[3:])) for line in coordinates_tmp]
            #CentNum = np.array(CentNum).reshape(1, -1)
            #atomicNum = np.array(atomicNum).reshape(1, -1)
            CentNum_AtomicN = np.array(CentNum_AtomicN)
            coordinates = np.array(coordinates)
            #print(CentNum_AtomicN)
            #print(atomicNum)
            #print(coordinates)

        #print(CentNum_AtomicN[:,0])
        #print(CentNum_AtomicN[:,1])
        NumAtom = int(np.amax(CentNum_AtomicN[:,0]))
        print(f'Number of atoms: {NumAtom}')
        Coordinfo = coordinates.shape
        #print(f'Coord: {Coordinfo}')

        atomicNum = CentNum_AtomicN[:,1]
        X = coordinates[:,0]
        Y = coordinates[:,1]
        Z = coordinates[:,2]

        #print(coordinates)
        #print (X)
        #print (Y)
        #print (Z)

        Mol_atom =[]
        for i in range(NumAtom):
            #print(atomicNum[0][i])
            if __name__ == '__main__':
                Mol_atom.append(AtomInfo.AtomicNumElec(int(atomicNum[i])))
            else:
                Mol_atom.append(gaussian_run.AtomInfo.AtomicNumElec(int(atomicNum[i])))
        #print(Mol_atom)

        return Mol_atom, X, Y, Z 

    def Estimate_SpinDiff(self, lines):
        print ("Start evaluating spin contamination")    
        Ideal_SS = 0.0
        Computed_SS = 0.0
        InputCharge, InputSpinMulti = self.extract_inputChargeSpin(lines)
        TotalS = (InputSpinMulti-1) / 2
        Ideal_SS = TotalS * (TotalS+1)
        for line in lines:
            if re.match("\s+<Sx>=", line):
                line_computed_Spin = line.split()
                Computed_SS = float(line_computed_Spin[-3])
                continue

        return Computed_SS, Ideal_SS

    def Extract_SCFEnergy(self, lines):
        Energy = []
        Eext = []
        for line in lines:
            if line.find("SCF Done:  ") >= 0:
                line_StateInfo = line.split()
                Energy.append(float(line_StateInfo[4]))
            if line.find("Total Energy, E(TD-HF/TD-DFT)") >=0:
                line_totalenergy = line.split('=')
                Eext.append(float(line_totalenergy[1]))

        Comp_SS, Ideal_SS = self.Estimate_SpinDiff(lines)
        if Eext == []:
            Energy_Spin = [Energy[-1], Comp_SS-Ideal_SS]
        else:
            Energy_Spin = [Eext[-1], Comp_SS-Ideal_SS]

        return Energy_Spin

    def Extract_symm(self, lines):
        pGroup = "C1"
        for line in lines:
            if line.find("Full point group  ") >= 0:
                line_symmInfo = line.split()
                pGroup = line_symmInfo[3]

        return pGroup

    def Extract_volume(self, lines):
        Mvolume =  0.0
        for line in lines:
            if line.find("Molar volume ") >= 0:
                line = line.replace('(', '').replace(')', '').replace('=','')
                line_volumeInfo = line.split()
                print(line_volumeInfo)
        try:    
            Mvolume = float(line_volumeInfo[2])*(bohr2cm**3)*Ac
        except:
            Mvolume = 0.0

        return  Mvolume

    def extract_inputChargeSpin(self, lines):
        charge = 0
        spinmulti  = 1
        text = '\n'.join(lines)
        pattern = r'\sCharge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)'
        match = re.search(pattern, text)

        if match:
            charge = int(match.group(1))
            spinmulti = int(match.group(2))

        return charge, spinmulti

    def classify_SpinStates(self, Ideal_SS, SS, WaveLength, OS, CD_L, CD_OS, mu, theta, g):
        State_allowed = []
        State_forbidden = []
        WL_allowed = []
        WL_forbidden = []
        OS_allowed = []
        OS_forbidden = []
        CD_L_allowed = []
        CD_L_forbidden = []
        CD_OS_allowed = []
        CD_OS_forbidden = []
        mu_allowed = []
        mu_forbidden = []
        theta_allowed = []
        theta_forbidden = []
        g_allowed = []
        g_forbidden = []
        for i in range(len(SS)):
            if abs(SS[i]-Ideal_SS) <= 0.1:
                State_allowed.append(i+1)
                WL_allowed.append(WaveLength[i])                   
                OS_allowed.append(OS[i])
                CD_L_allowed.append(CD_L[i])
                CD_OS_allowed.append(CD_OS[i])
                mu_allowed.append(mu[i])
                theta_allowed.append(theta[i])
                g_allowed.append(g[i])
            else:
                State_forbidden.append(i+1)
                WL_forbidden.append(WaveLength[i])                   
                OS_forbidden.append(OS[i])
                CD_L_forbidden.append(CD_L[i])
                CD_OS_forbidden.append(CD_OS[i])
                mu_forbidden.append(mu[i])
                theta_forbidden.append(theta[i])
                g_forbidden.append(g[i])

        return State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
                CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden, \
                mu_allowed, mu_forbidden, theta_allowed, theta_forbidden, g_allowed, g_forbidden

    def extract_transtion_EMmoment(self, lines):
        text = '\n'.join(lines)
        # Define start and end patterns for transition electric dipole moments
        e_start_pattern = "Ground to excited state transition electric dipole moments (Au):"
        e_end_pattern = "Ground to excited state transition velocity dipole moments (Au):"

        # Define start and end patterns for transition magnetic dipole moments
        m_start_pattern = "Ground to excited state transition magnetic dipole moments (Au):"
        m_end_pattern = "Ground to excited state transition velocity quadrupole moments (Au):"

        # Extract the relevant portion of the text
        e_relevant_text = re.search(f"{re.escape(e_start_pattern)}(.*?){re.escape(e_end_pattern)}", text, re.DOTALL).group(1)
        m_relevant_text = re.search(f"{re.escape(m_start_pattern)}(.*?){re.escape(m_end_pattern)}", text, re.DOTALL).group(1)

        # Split the relevant text into lines and extract the values
        lines = e_relevant_text.strip().split('\n')
        header = lines[0].split()
        te_dipole_values = [line.split() for line in lines[1:]]

        lines = m_relevant_text.strip().split('\n')
        header = lines[0].split()
        tm_dipole_values = [line.split() for line in lines[1:]]

        #print(te_dipole_values)
        #print(tm_dipole_values)
        mu = []
        mm = []
        
        for i in range(len(te_dipole_values)):
            mu.append(list(map(float, te_dipole_values[i][1:4])))
            mm.append(list(map(float, tm_dipole_values[i][1:4])))

        #print(mu)
        #print(mm)

        T_const = 2.54174E-18
        Y_const = 274.0717293
        V_const = 9.27401E-21 

        state_mu = []
        state_g = []
        state_theta = []

        for i in range(len(mu)):
            mu_state = np.array(mu[i])
            mm_state = np.array(mm[i])
            
            norm_mu_state = np.linalg.norm(mu_state)
            norm_mm_state = np.linalg.norm(mm_state)

            in_mu_mm = np.inner(mu_state,mm_state)

            mu_state_converted = norm_mu_state*T_const*Y_const
            state_mu.append(mu_state_converted)
             
            if norm_mu_state != 0 and norm_mm_state != 0:
                cos_theta = in_mu_mm/(norm_mu_state*norm_mm_state)
                valence = 4*(norm_mu_state*Y_const*norm_mm_state/((norm_mu_state**2)*Y_const**2+norm_mm_state**2))
            else:
                cos_theta = 0.0
                valence = 0.0

            theta = (np.arccos(cos_theta)/np.pi)*180
            g_calc = cos_theta*valence

            state_theta.append(theta)
            state_g.append(g_calc)
    
            #print(theta)

        return  state_mu, state_theta, state_g

    def extract_CD(self, lines):
        count = 0
        CD_length = []
        CD_Osc = []
        for line in lines:
            if line.find("1/2[<0|r|b>*<b|rxdel|0> + (<0|rxdel|b>*<b|r|0>)*]") >= 0: 
                CD_length = []
                CD_Osc = []
                count += 1
                continue
            if count == 1:
                count += 1
                continue
            if count == 2:
                count += 1
                continue
            if count == 3:
                if line.find(" 1/2[<0|del|b>*<b|r|0> + (<0|r|b>*<b|del|0>)*] (Au)") >= 0: 
                    count += 1
                    continue
                else:
                    info_R =  line.split()
                    cd_value = float(info_R[-1])
                    CD_length.append(cd_value)    
                    continue
            if count == 4:
                count += 1 
                continue
            if count == 5:
                if line == "": 
                    count = 0
                    continue
                else:
                    info_Osc =  line.split()
                    osc_value = float(info_Osc[-1])
                    CD_Osc.append(osc_value)    
                    continue
        return CD_length, CD_Osc

    def Extract_ExcitedState(self, lines):
        Egrd = 0.0
        Eext = 0.0
        Found = False
        WaveLength = []
        V_OS = []
        SS = []
        _, Ideal_SS = self.Estimate_SpinDiff(lines)
        
        print ("Get information about excited state")
        for line in lines:
            if line.find("SCF Done:  ") >=0:
                line_SCFEnergy = re.split("\s+", line)
                Egrd = float(line_SCFEnergy[5])
            if line.find("Total Energy, E(TD-HF/TD-DFT)") >=0:
                line_totalenergy = line.split('=')
                Eext = float(line_totalenergy[1])
            if line.find("Excitation energies and oscillator strengths:") >=0:
                 WaveLength = []
                 V_OS = []
                 SS = []
            if line.find("Excited State  ") >=0:
                 line_StateInfo = line.split()
                 WaveLength.append(float(line_StateInfo[6]))
                 OS_info = line_StateInfo[8].split('=')
                 V_OS.append(float(OS_info[1]))
                 SS_info = line_StateInfo[9].split('=')
                 SS.append(float(SS_info[1]))
            if line.find("-- Stationary point found.") >=0:
                Found = True
        
        CD_length, CD_OS = self.extract_CD(lines)
        mu, theta, g = self.extract_transtion_EMmoment(lines)
        State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
                CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden, \
                mu_allowed, mu_forbidden, theta_allowed, theta_forbidden, g_allowed, g_forbidden \
                = self.classify_SpinStates(Ideal_SS, SS, WaveLength, V_OS, CD_length, CD_OS, mu, theta, g)


        return Found, Egrd, Eext, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
            CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden, \
            mu_allowed, mu_forbidden, theta_allowed, theta_forbidden, g_allowed, g_forbidden  

    def extract_Freq(self, lines):
        count = 0
        Freq = []
        IR_int = []
        Raman_int = []

        print ("Get information about  Frequencies")
        for line in lines:
            if line.find("Frequencies --") >= 0:
                line_Freq = line.split()
                for i in range(2, len(line_Freq)):
                    Freq.append(float(line_Freq[i]))
            if line.find("IR Inten    --") >= 0:
                line_ir = line.split()
                for i in range(3, len(line_ir)):
                    IR_int.append(float(line_ir[i]))
            if line.find("Raman Activ --") >= 0:
                line_raman = line.split()
                for i in range(3, len(line_raman)):
                    Raman_int.append(float(line_raman[i]))
            if line.find("Sum of electronic and zero-point Energies=") >= 0:
                line_zp = line.split()
                E_ZP = float(line_zp[-1])
            if line.find("Sum of electronic and thermal Energies=") >= 0:
                line_t = line.split()
                E_t = float(line_t[-1])
            if line.find("Sum of electronic and thermal Enthalpies=") >= 0:
                line_t = line.split()
                E_enth = float(line_t[-1])
            if line.find("Sum of electronic and thermal Free Energies=") >= 0:
                line_free = line.split()
                E_free = float(line_free[-1])
            if line.find("KCal/Mol") >= 0:
                count += 1
                continue
            if count == 1:
                line_ECvS = line.split()
                Ei = line_ECvS[1]
                Cv = line_ECvS[2]
                St= line_ECvS[3]
                count = 0
        
        return  Freq, IR_int, Raman_int, E_ZP,  E_t, E_enth, E_free, Ei, Cv, St

    def extract_vibvec(self, lines):
        result_lines = []
        vibX = []
        vibY = []
        vibZ = []
        
        count = 0
        
        for line in lines:
            pattern_match = re.match(r'\s*\bAtom\s+AN\s+X\s+Y\s+Z\s+X\s+Y\s+Z\s+X\s+Y\s+Z\b', line)
            if pattern_match:
                count += 1
            match = re.match(r'\s*(\d+)\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+', line)
            if match:
                vibvec_info = line.strip().split()
                #print(vibvec_info[2:])
                for i in range(2, len(vibvec_info), 3): 
                    vibX.append(float(vibvec_info[i]))
                    vibY.append(float(vibvec_info[i+1]))
                    vibZ.append(float(vibvec_info[i+2]))
        
                result_lines.append(line)
        
        vibX = np.array(vibX)
        vibY = np.array(vibY)
        vibZ = np.array(vibZ)
        
        vibX_tmp = vibX.reshape(count, -1)
        vibY_tmp = vibY.reshape(count, -1)
        vibZ_tmp = vibZ.reshape(count, -1)
        
        vibX_reshape = vibX_tmp[0].reshape(-1, 3).T
        vibY_reshape = vibY_tmp[0].reshape(-1, 3).T
        vibZ_reshape = vibZ_tmp[0].reshape(-1, 3).T
        
        for i in range(1, count):
            vibX_reshape = np.vstack([vibX_reshape, vibX_tmp[i].reshape(-1, 3).T])
            vibY_reshape = np.vstack([vibY_reshape, vibY_tmp[i].reshape(-1, 3).T])
            vibZ_reshape = np.vstack([vibZ_reshape, vibZ_tmp[i].reshape(-1, 3).T])
        
        #print(vibX_reshape)
        #print(vibY_reshape)
        #print(vibZ_reshape)
        
        return vibX_reshape, vibY_reshape, vibZ_reshape

    def extract_polar(self, lines):

        activate_extraction = False
        iso_values = []
        aniso_values = []
        
        for line in lines:
           # if line.strip() == "Dipole polarizability, Alpha (dipole orientation).":
           #     activate_extraction = True
            if line.strip() == "Dipole polarizability, Alpha (input orientation).":
                activate_extraction = True
            if activate_extraction:
                if line.strip().startswith("iso"):
                    values = re.findall(r"\d+\.\d+[Dd][\+\-]?\d+", line)
                    #iso_values.extend(values)
                    iso_values.extend([float(val.replace('D', 'E')) for val in values])
                elif line.strip().startswith("aniso"):
                    values = re.findall(r"\d+\.\d+[Dd][\+\-]?\d+", line)
                    #aniso_values.extend(values)
                    aniso_values.extend([float(val.replace('D', 'E')) for val in values])
        
        return iso_values, aniso_values
        
    def classify_task(self, lines, charge, spinmulti):
        count = 0
        GScharge = 0
        GSspin = 1
    
        job_index = {}
    
        for i in range(len(lines)):
            print(lines[i])
            if re.search('\sGuess=Only', lines[i]):
                #print('SCF energy is not available')
                if re.search('\sSymmetry', lines[i]):
                    symm = True
                    job_index['symm'] = i
                if re.search('\svolume', lines[i]):
                    is_volume = True
                    job_index['volume'] = i
            else:
                count += 1
                if count == 1:
                    #print (f'GS electronic structure is charge {charge[i]} and spin multiplicity {spinmulti[i]}.')
                    GScharge = charge[i]
                    GSspin = spinmulti[i]
                    #if re.search('\sOpt', lines[i]) and re.search('\sTD', lines[i]) == None:
                    if re.search('\sTD', lines[i]) == None:
                        is_opt = True
                        job_index['ts'] = i
                        job_index['state_0'] = i
                    if re.search('\sOpt', lines[i]) and re.search('\sTD', lines[i]) and re.search('Singlet', lines[i]):
                        is_fluor = True
                        match_root = re.search('\s+root=\d+', lines[i])
                        if match_root:    
                            root_line = match_root.group().split('=')
                            target = root_line[1]
                        #job_index[f'relaxAEstate{target}_line'] = i
                        job_index['ts'] = i
                        job_index['relaxAEstate'] = i
                        job_index[f'state_{target}'] = i
                    if re.search('\sOpt', lines[i]) and re.search('\sTD', lines[i]) and re.search('Triplet', lines[i]):
                        is_tadf = True
                        match_root = re.search('\s+root=\d+', lines[i])
                        if match_root:    
                            root_line = match_root.group().split('=')
                            target = root_line[1]
                        #job_index[f'relaxFEstate{target}_line'] = i
                        job_index['ts'] = i
                        job_index['relaxFEstate'] = i
                        job_index[f'state_{target}'] = i
                    if re.search('\sOpt', lines[i]) and re.search('\sTD', lines[i]) and re.search('Singlet', lines[i]) == None and re.search('Triplet', lines[i]) == None: 
                        is_fluor = True
                        match_root = re.search('\s+root=', lines[i])
                        if match_root:    
                            root_line = match_root.group().split('=')
                            target = root_line[1]
                        job_index['ts'] = i
                        job_index['relaxAEstate'] = i
                        job_index[f'state_{target}'] = i
                        #if i+1 < len(lines):
                        #    job_index[f'relaxFEstate'] = i+1
                if count > 1:
                    if re.search('\sOpt', lines[i]) and re.search('\sTD', lines[i]) and re.search('Singlet', lines[i]):
                        is_fluor = True
                        match_root = re.search('\s+root=\d+', lines[i])
                        if match_root:    
                            root_line = match_root.group().split('=')
                            target = root_line[1]
                        #job_index[f'relaxAEstate{target}_line'] = i
                        job_index['relaxAEstate'] = i
                        job_index[f'state_{target}'] = i
                    if charge[i] == GScharge + 1:
                        #print ('Ionization computation')
                        if re.search('\sOpt', lines[i]):
                            is_aip = True
                            job_index['PC_line'] = i
                            if i+1 <len(lines):
                                job_index['VNP_line'] = i+1
                        else:
                            is_vip = True
                            job_index['IP_line'] = i
                    if charge[i] == GScharge - 1:
                        #print ('Electronic affinity computation')
                        if re.search('\sOpt', lines[i]):
                            is_aea = True
                            job_index['NC_line'] = i
                            if i+1 <len(lines):
                                job_index['VNN_line'] = i+1
                        else:
                            is_vea = True
                            job_index['EA_line'] = i
                if re.search('\spolar', lines[i]):
                    is_polar = True
                    job_index['polar_line'] = i
                if re.search('\sFreq', lines[i]):
                    is_freq = True
                    is_polar = True
                    job_index['freq_line'] = i
                if re.search('\sNMR', lines[i]):
                    is_nmr = True
                    job_index['nmr_line'] = i
                if re.search('\sTD', lines[i]) and re.search('\sOpt', lines[i]) == None:
                    is_uv = True
                    job_index['uv_line'] = i
    
        count = 0
    
        return job_index 

    def Check_task(self):
        Links = self.SplitLinks()
        Links.pop(0)
        method = [] 
        functional = []
        basis = []
        charge = []
        spinmulti = []
        Links_split = []
        AtomNum =  []
        MolX =  []
        MolY =  []
        MolZ =  []
        for i in range(len(Links)):
            lines = Links[i].splitlines()
            #extract method line
            method_text = self.extract_method_text(lines)
            method.append(method_text) 
            #extract functional and basis set from method lines
            afunctional, abasis = self.extract_functionalbasis(method_text)
            functional.append(afunctional)
            basis.append(abasis)
            #extract charge and spin multiplicity
            charge_link, spinmulti_link = self.extract_inputChargeSpin(lines)
            charge.append(charge_link)           
            spinmulti.append(spinmulti_link)
            Links_split.append(lines)

            AtomNum_tmp, Mol_tmpX, Mol_tmpY, Mol_tmpZ = self.extract_MolCoordlog(lines)
            AtomNum.append(AtomNum_tmp)
            MolX.append(Mol_tmpX)
            MolY.append(Mol_tmpY)
            MolZ.append(Mol_tmpZ)


        print (AtomNum)
        print (MolX)
        print (MolY)
        print (MolZ)

        #print (functional, basis)
        job_index = self.classify_task(method, charge, spinmulti)

###test...
#        print(self.Extract_SCFEnergy(Links_split[job_index['gs']]))
#        print(self.Extract_SCFEnergy(Links_split[job_index['PC_line']]))
#        print(self.Extract_symm(Links_split[job_index['symm']]))
######            
        return job_index, Links_split

if __name__ == '__main__':
    usage ='Usage; %s infile' % sys.argv[0]
    import AtomInfo

    try:
        infilename = sys.argv[1]
    except:
        print (usage); sys.exit()

    parselog = parse_log(infilename)
    job_index, Links_split = parselog.Check_task()
    print(job_index)

    output = {}

    if 'uv_line' in job_index:
        lines = Links_split[job_index['uv_line']] 
        _, _, _, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
            CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden, \
            mu_allowed, mu_forbidden, theta_allowed, theta_forbidden, g_allowed, g_forbidden = parselog.Extract_ExcitedState(lines)
        output["uv"] = [WL_allowed, OS_allowed, CD_L_allowed, CD_OS_allowed]
        output["state_index"] = [State_allowed, State_forbidden]
        output["CD_mu"] = [mu_allowed, mu_forbidden]
        output["CD_theta"] = [theta_allowed, theta_forbidden]
        output["CD_g"] = [g_allowed, g_forbidden]

    if 'relaxAEstate' in job_index:
        lines = Links_split[job_index['relaxAEstate']] 
        S_Found, S_Egrd, S_Eext, State_allowed, State_forbidden, WL_allowed, WL_forbidden, OS_allowed, OS_forbidden, \
            CD_L_allowed, CD_L_forbidden, CD_OS_allowed, CD_OS_forbidden, \
            mu_allowed, mu_forbidden, theta_allowed, theta_forbidden, g_allowed, g_forbidden = parselog.Extract_ExcitedState(lines)
        # For Check internal coordinate
        output["MinEtarget"] = S_Eext
        output["fluor"] = [WL_allowed, OS_allowed, CD_L_allowed, CD_OS_allowed]
        #mu, theta, g_cal = parse_log.extract_transtion_EMmoment(lines)
        output["CDL_mu"] = mu_allowed 
        output["CDL_theta"] = theta_allowed
        output["CDL_g"] = g_allowed

    if 'freq_line' in job_index:
        lines = Links_split[job_index['freq_line']] 
        Polar_iso, Polar_aniso = parselog.extract_polar(lines)
        output["polar_iso"] = Polar_iso[1]
        output["polar_aniso"] = Polar_aniso[1]

    print (output)

else:
    from qcforever import gaussian_run

