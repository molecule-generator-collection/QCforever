import re
import sys
import os


def read_dat(datfile):

    with open(datfile) as dat_line:
        lines = dat_line.readlines()

    return lines


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


    bb = get_dataBlock(llines,"VEC")

    MOvec_line = ' $VEC\n'
    with open('vec_test.txt', 'w') as ofile:
        for i in range(len(bb)):
            MOvec_line += bb[i]
        MOvec_line += ' $END\n'

        ofile.write(MOvec_line)

