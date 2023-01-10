import re
import sys
import os


def read_log(outfile):

    with open (outfile) as output_line:
        lines = output_line.readlines()

    return lines


def getNumberElectron(lines):

    for ii in range(len(lines)):
        ll = lines[ii]

        mmat = re.search("^[\s]*NUMBER OF OCCUPIED ORBITALS[^=]+=[\s]*([^\s]+)",ll)
        if(mmat):
            if(re.search("ALPHA",ll)):
                num_occu_alpha = int(mmat.group(1))
            if(re.search("BETA",ll)):
                num_occu_beta = int(mmat.group(1))

    return num_occu_alpha, num_occu_beta


def getEnergy(lines):

    energy = []

    for ii in range(len(lines)):
        ll = lines[ii]
       
        if(re.search("^[\s]*FINAL", ll)):
            #print(ll)
            sline = ll.split()
            energy.append(sline[-4])

    return energy[-1]


def getBlock(lines, label):
    flag = 0
    ret = []
    currentlist = []
    ii = 0
    while(ii < len(lines)):
        ll = lines[ii]
        
        if(flag == 1 and (re.search("^[\s]*----",ll) or re.search("\.\.\.\.\.\.",ll))):
            flag = 0
            ret.append(currentlist)
            currentlist = []
            
        if(re.search("^[\s]*"+label,ll)):
            if(re.search("^[\s]*----",lines[ii-1]) and re.search("^[\s]*----",lines[ii+1])):
                flag = 1
                ii += 2
                continue
                
        if(flag == 1):
            currentlist.append(ll)
        
        ii+=1
    if(len(currentlist)):
        ret.append(currentlist)
    return ret
    
def getMO_set(block):
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

def getMO_single(block):
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


def gethomolumogap(alpha_values, beta_values, num_alpha_elec, num_beta_elec):
    
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
#    if(ret1 == None):
#        raise Exception("???? "+str(targetindex)+" not found");
#    return [ret1,ret2];

    return  alpha_gap, beta_gap


def getDipoleMoment(block):
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


    flag_homolumo = False
    flag_dipole = False

    llines = read_log(sys.argv[1])

    num_occu_alpha, num_occu_beta  = getNumberElectron(llines)

    print(getEnergy(llines))

    bb = getBlock(llines,"EIGENVECTORS")
    #print(bb[-2])
    #print(bb[-1])
    alpha_values = getMO_single(bb[-2])
    beta_values = getMO_single(bb[-1])
    #bb = getBlock(llines,"MOLECULAR ORBITALS")
    #alpha_values, beta_values = getMO_set(bb[-1])
    alpha_gap, beta_gap = gethomolumogap(alpha_values, beta_values, num_occu_alpha, num_occu_beta)
    dd = getBlock(llines,"ELECTROSTATIC MOMENTS")    
    dval = getDipoleMoment(dd[-1])
    #print("HOMO:\t"+str(hval[0]))
    #print("LUMO:\t"+str(hval[1]))
    print("HOMO - LUMO gap:\t"+str(alpha_gap))
    print("HOMO - LUMO gap:\t"+str(beta_gap))
    print("DIPOLEMOMENT:\t"+str(dval))

