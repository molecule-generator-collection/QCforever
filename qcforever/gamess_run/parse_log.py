import re
import sys
import os


class parse_log:

    def __init__(self, outfile):

        with open (outfile) as output_line:
            self.lines = output_line.readlines()

    def Check_SCF(self):

        scf_state = True

        for ii in range(len(self.lines)):
            ll = self.lines[ii]
            if (re.match("\s*SCF IS UNCONVERGED",ll)):
                scf_state = "scf error"

        return scf_state
        
    def getNumberElectron(self):

        for ii in range(len(self.lines)):
            ll = self.lines[ii]

            mmat = re.search("^[\s]*NUMBER OF OCCUPIED ORBITALS[^=]+=[\s]*([^\s]+)",ll)
            if(mmat):
                if(re.search("ALPHA",ll)):
                    num_occu_alpha = int(mmat.group(1))
                if(re.search("BETA",ll)):
                    num_occu_beta = int(mmat.group(1))

        return num_occu_alpha, num_occu_beta

    def getEnergy(self):

        energy = []

        for ii in range(len(self.lines)):
            ll = self.lines[ii]
       
            if(re.search("^[\s]*FINAL", ll)):
                #print(ll)
                sline = ll.split()
                energy.append(sline[-4])

        return float(energy[-1])

    def getTDDFT(self):

        Wavel = []
        OS = []
        flag = 0
        ii = 0
        while(ii < len(self.lines)):
            ll = self.lines[ii]
        
            if (flag == 0 and re.search("SUMMARY OF TDDFT RESULTS", ll)):
                flag = 1
                Wavel = []
                OS = []
                ii += 2
            if (flag == 1 and re.match("\s+[0-9]",ll) and not re.search('HARTREE', ll)):
                sline = ll.split()
                if len(sline) > 4: 
                    Wavel.append(1240/float(sline[-5]))
                    OS.append(float(sline[-1]))
            if (flag == 1 and  re.match("\s+\n", ll)):
                flag = 0
            ii += 1
        
        return Wavel, OS

    def getBlock(self, label):
        flag = 0
        ret = []
        currentlist = []
        ii = 0
        while(ii < len(self.lines)):
            ll = self.lines[ii]
            
            if(flag == 1 and (re.search("^[\s]*----",ll) or re.search("\.\.\.\.\.\.",ll))):
                flag = 0
                ret.append(currentlist)
                currentlist = []
                
            if(re.search("^[\s]*"+label,ll)):
                if(re.search("^[\s]*----",self.lines[ii-1]) and re.search("^[\s]*----",self.lines[ii+1])):
                    flag = 1
                    ii += 2
                    continue
                    
            if(flag == 1):
                currentlist.append(ll)
            
            ii+=1
        if(len(currentlist)):
            ret.append(currentlist)
        return ret
    
    def getMO_set(self, block):
        flag = 0
        ret = []
        currentlist = []
        ii = 0
        
        hitflag = 0
        
        elec_flag = 0
        alpha_indices = []
        alpha_values = []
        beta_indices = []
        beta_values = []
        while(ii < len(block)-1):
            ll = block[ii]
            nex = block[ii+1]
            if(re.search(r"ALPHA SET",ll)):
                #print (ll)
                elec_flag = 1
            if(re.search(r"BETA SET",ll)):
                #print (ll)
                elec_flag = -1
            if(re.search("^          +[0-9]+ ",ll) and re.search("^          ",nex) and not re.search("[A-DF-Za-df-z]",nex)):
                ipt = re.split("[\s]+",re.sub("[\s]+$","",re.sub("^[\s]+","",ll)))
                vpt = re.split("[\s]+",re.sub("[\s]+$","",re.sub("^[\s]+","",nex)))
                for pp in range(len(ipt)):
                    if elec_flag == 1:
                        alpha_indices.append(int(ipt[pp]))
                        alpha_values.append(float(vpt[pp]))
                    if elec_flag == -1:
                        beta_indices.append(int(ipt[pp]))
                        beta_values.append(float(vpt[pp]))
            ii += 1
        
        #print(alpha_indices)
        #print(alpha_values)
        #print(beta_indices)
        #print(beta_values)
        
        if(len(alpha_indices) != len(alpha_values) or len(beta_indices) != len(beta_values)):
            raise Exception("???different length bet index list and value list. parsing error???")
        
        return alpha_values, beta_values

    def getMO_single(self, block):
        ret = []
        currentlist = []
        ii = 0
        
        hitflag = 0
        
        elec_flag = 0
        indices = []
        values = []
        while(ii < len(block)-1):
            ll = block[ii]
            nex = block[ii+1]
            if(re.search("^          +[0-9]+ ",ll) and re.search("^          ",nex) and not re.search("[A-DF-Za-df-z]",nex)):
                ipt = re.split("[\s]+",re.sub("[\s]+$","",re.sub("^[\s]+","",ll)))
                vpt = re.split("[\s]+",re.sub("[\s]+$","",re.sub("^[\s]+","",nex)))
                for pp in range(len(ipt)):
                        indices.append(int(ipt[pp]))
                        values.append(float(vpt[pp]))
            ii += 1
        
        #print(indices)
        #print(values)
        
        if(len(indices) != len(values)) :
            raise Exception("???different length bet index list and value list. parsing error???")
        
        return values

    def gethomolumogap(self, alpha_values, beta_values, num_alpha_elec, num_beta_elec):
    
        ret1 = None
        ret2 = None
        for ii in range(len(alpha_values)):
            if(ii == num_alpha_elec-1):
                ret1 = alpha_values[ii]
                ret2 = alpha_values[ii+1]
        alpha_gap  =(float(ret2)-float(ret1))*27.211
        
        for ii in range(len(beta_values)):
            if(ii == num_beta_elec-1):
                ret1 = beta_values[ii]
                ret2 = beta_values[ii+1]
        beta_gap  =(float(ret2)-float(ret1))*27.211
        
        return  alpha_gap, beta_gap

    def getDipoleMoment(self, block):
        ii = 0
        ret = None
        while(ii < len(block)-1):
            ll = block[ii]
            nex = block[ii+1]
            #/D/
            if(re.search("/D/",ll)):
                ipt = re.split("[\s]+",re.sub("[\s]+$","",re.sub("^[\s]+","",ll)))
                vpt = re.split("[\s]+",re.sub("[\s]+$","",re.sub("^[\s]+","",nex)))
                
                for pp in range(len(ipt)):
                    if(ipt[pp] == "/D/"):#/D/
                        ret =float(vpt[pp])
            ii += 1
#    if(ret == None):
#        raise Exception("??? data pasing error?"+"\n".join(block));
    
        return ret


if __name__ == '__main__':
    usage ='Usage; %s infile' % sys.argv[0]

    try:
        infilename = sys.argv[1]
    except:
        print (usage); sys.exit()

    parselog = parse_log(infilename)
    print(parselog.Check_SCF())

    flag_homolumo = False
    flag_dipole = False

    num_occu_alpha, num_occu_beta  = parselog.getNumberElectron()

    print(parselog.getEnergy())

    Wavel, OS = parselog.getTDDFT()

    print ("Wave length:", Wavel)
    print ("OS:", OS)

    bb = parselog.getBlock("EIGENVECTORS")
    #print(bb[-2])
    #print(bb[-1])
    alpha_values = parselog.getMO_single(bb[-2])
    beta_values = parselog.getMO_single(bb[-1])
    #bb = getBlock(llines,"MOLECULAR ORBITALS")
    #alpha_values, beta_values = getMO_set(bb[-1])
    alpha_gap, beta_gap = parselog.gethomolumogap(alpha_values, beta_values, num_occu_alpha, num_occu_beta)
    dd = parselog.getBlock("ELECTROSTATIC MOMENTS")    
    dval = parselog.getDipoleMoment(dd[-1])
    #print("HOMO:\t"+str(hval[0]))
    #print("LUMO:\t"+str(hval[1]))
    print("HOMO - LUMO gap:\t"+str(alpha_gap))
    print("HOMO - LUMO gap:\t"+str(beta_gap))
    print("DIPOLEMOMENT:\t"+str(dval))

