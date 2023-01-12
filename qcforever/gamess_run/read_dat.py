import re
import sys
import os


def read_dat(datfile):

    with open(datfile) as dat_line:
        lines = dat_line.readlines()

    return lines


def getMolBlock(lines, label):
    flag = 0
    ret = []
    currentlist = []
    ii = 0
    while(ii < len(lines)):
        ll = lines[ii]
        
        if(flag == 1 and (re.search("^[\s]*---",ll))):
            flag = 0
            ret.append(currentlist)
            currentlist = []
            
       # if(re.search("^[\s]*"+label,ll)):
        if(re.search("^[\s]"+label,ll)):
            flag = 1
            ii += 3
            continue
                
        if(flag == 1):
            currentlist.append(ll)
        
        ii+=1
    if(len(currentlist)):
        ret.append(currentlist)
    return ret[-1]


def get_dataBlock(lines, label):
    flag = 0
    ret = []
    currentlist = []
    ii = 0
    while(ii < len(lines)):
        ll = lines[ii]

        if(flag == 1 and (re.search("\$END", ll))):
            flag = 0
            ret.append(currentlist)
            currentlist = []

        if(re.search("^[\s]*"+"\$"+label, ll)):
            flag = 1    
            ii += 1
            continue

        if(flag == 1):
            currentlist.append(ll)

        ii += 1

    if(len(currentlist)):
        ret.append(currentlist)

    return ret[-1]


if __name__ == '__main__':
    usage ='Usage; %s infile' % sys.argv[0]

    try:
        infilename = sys.argv[1]
    except:
        print (usage); sys.exit()


    llines = read_dat(sys.argv[1])

#Molecular coordinate
    cc = getMolBlock(llines,"COORDINATES OF SYMMETRY UNIQUE ATOMS")
    print (cc)

    dat_line = ' $DATA\n'
    dat_line += 'Comment\n'
    dat_line += 'C1\n'
    for i in range(len(cc)):
        dat_line += cc[i]
    dat_line += ' $END\n'

#Molecular orbitals
    bb = get_dataBlock(llines,"VEC")

    MOvec_line = ' $VEC\n'
    for i in range(len(bb)):
        MOvec_line += bb[i]
    MOvec_line += ' $END\n'

    with open('mol_test.txt', 'w') as ofile:
        ofile.write(dat_line)
        ofile.write(MOvec_line)




