import re
import sys
import os
import numpy as np

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

    def extract_enclosed_text(self, lines):
        text = '\n'.join(lines)
        pattern = r'\s#.*?(?=\s-{35,70})'  
        enclosed_text = re.findall(pattern, text, re.DOTALL)
        return enclosed_text

    def extract_MolCoordlog(self, lines):
        text = '\n'.join(lines)
        pattern = r"Standard orientation:(?:\s+(-+\n\s+Center\s+Atomic\s+Atomic\s+Coordinates \(Angstroms\)\n\s+Number\s+Number\s+Type\s+X\s+Y\s+Z\n\s+(-+\n\s+\d+\s+\d+\s+\d+\s+[-.\d\s]+\n)+)(?=\s+-+\n))"
        Coord_blocks = re.findall(pattern, text)

        for block in Coord_blocks:
            print ('Atomic coordinates are found.')
            pattern = r'\s*\d+\s+\d+\s+\d+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+'
            coordinates_inf = re.findall(pattern, str(block[1]))
            coordinates_tmp = []
            for mcoord in coordinates_inf:
                coordinates_tmp.append(mcoord.strip())
                #print(coordinates_tmp)
            coordinates = [list(map(float, line.split()[3:])) for line in coordinates_tmp]
            coordinates = np.array(coordinates)
            #print(coordinates)

        return coordinates

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
        for line in lines:
            if line.find("SCF Done:  ") >= 0:
                line_StateInfo = line.split()
                Energy.append(float(line_StateInfo[4]))
        Comp_SS, Ideal_SS = self.Estimate_SpinDiff(lines)
        Energy_Spin = [Energy[-1], Comp_SS-Ideal_SS]

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
                line = line.replace('(', '').replace(')', '')
                line_volumeInfo = line.split()
        Mvolume = float(line_volumeInfo[-2])

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
                    if re.search('\sOpt', lines[i]):
                        is_opt = True
                    job_index['gs'] = i
                if count > 1:
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
                if re.search('\sOpt', lines[i]) and re.search('\sTD', lines[i]) and re.search('Singlet', lines[i]):
                    is_fluor = True
                    match_root = re.search('\s+root=\d+', lines[i])
                    if match_root:    
                        root_line = match_root.group().split('=')
                        #print(root_line)
                        target = root_line[1]
                    #job_index[f'relaxAEstate{target}_line'] = i
                    job_index[f'relaxAEstate'] = i
                if re.search('\sOpt', lines[i]) and re.search('\sTD', lines[i]) and re.search('Triplet', lines[i]):
                    is_tadf = True
                    match_root = re.search('\s+root=\d+', lines[i])
                    if match_root:    
                        root_line = match_root.group().split('=')
                        target = root_line[1]
                    #job_index[f'relaxFEstate{target}_line'] = i
                    job_index[f'relaxFEstate'] = i
                if re.search('\sOpt', lines[i]) and re.search('\sTD', lines[i]) and re.search('Singlet', lines[i]) == None and re.search('Triplet', lines[i]) == None: 
                    is_fluor = True
                    match_root = re.search('root=', lines[i])
                    if match_root:    
                        root_line = match_root.group().split('=')
                        target = root_line[1]
                    job_index[f'relaxAEstate{target}_line'] = i

        count = 0

        return job_index 

    def Check_task(self):
        Links = self.SplitLinks()
        Links.pop(0)
        method = [] 
        charge = []
        spinmulti = []
        Links_split = []
        #MolCoord =  []
        for i in range(len(Links)):
            lines = Links[i].splitlines()
            method_line = self.extract_enclosed_text(lines)
            method_text = ''
            for text in method_line:
                method_text += re.sub(r'\n\s', '', text.strip())
            method.append(method_text) 

            charge_link, spinmulti_link = self.extract_inputChargeSpin(lines)
            charge.append(charge_link)           
            spinmulti.append(spinmulti_link)
            Links_split.append(lines)

            #MolCoord.append(self.extract_MolCoordlog(lines))

        job_index = self.classify_task(method, charge, spinmulti)

###test...
#        print(self.Extract_SCFEnergy(Links_split[job_index['gs']]))
#        print(self.Extract_SCFEnergy(Links_split[job_index['PC_line']]))
#        print(self.Extract_symm(Links_split[job_index['symm']]))
######            
        return job_index, Links_split

if __name__ == '__main__':
    usage ='Usage; %s infile' % sys.argv[0]

    try:
        infilename = sys.argv[1]
    except:
        print (usage); sys.exit()

    parselog = parse_log(infilename)
    job_index, Links_split = parselog.Check_task()
    print(job_index)

