import sys
from qcforever.gamess_run import GamessRunPack
				
usage ='Usage; %s infile' % sys.argv[0]

try:
    infilename = sys.argv[1]
except:
    print (usage); sys.exit()


#for starting from gamess calculation
option = "opt homolumo energy dipole uv fluor"
#option = "opt uv vip vea"

test = GamessRunPack.GamessDFTRun('LC-BLYP', '3-21G', 8, option, infilename)
#test = GamessRunPack.GamessDFTRun('BHHLYP', '3-21G', 8, option, infilename)

test.mem = '5GB'
test.timexe = 60*60 #second unit
test.para_functional = [0.65]
#test.SpecSpinMulti = 3
#test.SpecTotalCharge = 3
outdic = test.run_gamess()

print (outdic)
