import os, sys
from pprint import pprint
import GaussianRunPack
				
usage ='Usage; %s infile' % sys.argv[0]

try:
    infile = sys.argv[1]
except:
    print (usage); sys.exit()

ifile = open(infile, 'r')
lines = ifile.readlines()
ifile.close()

option = "deen UV "
test = GaussianRunPack.GaussianDFTRun('LC-BLYP', 'STO-3G', 20, option, infile ,0)

#test.Extract_Coordinate(lines)

test.Extract

print (test)

