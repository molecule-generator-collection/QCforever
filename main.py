import os, sys, gc
#from pprint import pprint
import GaussianRunPack
				
usage ='Usage; %s infile' % sys.argv[0]

try:
    infilename = sys.argv[1]
except:
    print (usage); sys.exit()


#if len(sys.argv)==3: option = sys.argv[2]

#######For reading log file ####
#option = "eae homolumo"
#option_array = [1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#bond1 = []
#bond2 = []
#test_sdf = GaussianRunPack.GaussianDFTRun('B3LYP', 'STO-3G', 12, option, infilename)
#outdic = test_sdf.Extract_values(infilename,option_array, bond1 , bond2 )
#print (outdic)
###################################


#for stargint from gaussian calculation
#option = "symm uv"
option = "opt Freq"
#option = "opt deen vip vea homolumo dipole uv symm"
#option = "symm opt freq nmr uv=/home/sumita/GaussianRun_2.2beta/Samples/UV_peak.dat" 
#option = "opt homolumo energy dipole deen stable2o2 uv=/home/sumita/GaussianRun_2.2beta/Samples/UV_peak.dat vip vea aip aea" 
#option = "opt homolumo energy dipole deen stable2o2 uv vip vea aip aea" 
#option = "opt homolumo energy dipole deen stable2o2 fluor=3" 
#option = "opt cden homolumo energy dipole deen stable2o2 tadf" 
#option = "opt nmr uv energy homolumo dipole deen stable2o2 vip vea cden symm"
#solvent = "water"

test_sdf = GaussianRunPack.GaussianDFTRun('B3LYP', 'STO-3G', 8, option, infilename)

test_sdf.mem = '5GB'
test_sdf.timexe = 60*60
outdic = test_sdf.run_gaussian()

#del test_sdf
#gc.collect()

#print(os.getcwd())

print (outdic)
#######################################
