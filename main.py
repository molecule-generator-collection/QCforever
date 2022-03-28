import os, sys, gc
#from pprint import pprint
import GaussianRunPack
				
usage ='Usage; %s infile' % sys.argv[0]

try:
    infilename = sys.argv[1]
except:
    print (usage); sys.exit()


#if len(sys.argv)==3: option = sys.argv[2]

#For reading log file ####
#option = "eae homolumo"
#test_sdf = GaussianRunPack.GaussianDFTRun('B3LYP', 'STO-3G', 12, option,infilename,0)
#outdic = test_sdf.Extract_values(infilename,0,0,0,0,0,1,1,0,0,0,0)
#print (outdic)
#########################


#for stargint from gaussian calculation
#option = "symm uv"
#option = "opt deen nmr stable2o2"
#option = "opt deen vip vea homolumo dipole uv symm"
#option = "uv=/home/sumita/GaussianRun_2.1beta/Samples/UV_peak.dat" 
#option = "opt homolumo energy dipole deen stable2o2 uv=/home/sumita/GaussianRun_2.1beta/Samples/UV_peak.dat vip vea aip aea" 
#option = "opt homolumo energy dipole deen stable2o2 uv vip vea aip aea" 
option = "opt homolumo energy dipole deen stable2o2 fluor=3" 
#option = "opt cden homolumo energy dipole deen stable2o2 tadf" 
#solvent = "water"

test_sdf = GaussianRunPack.GaussianDFTRun('B3LYP', 'STO-3G', 8, option, infilename)

test_sdf.mem = '5GB'
test_sdf.timexe = 60*60
outdic = test_sdf.run_gaussian()

del test_sdf
gc.collect()

print(os.getcwd())

print (outdic)
#######################################
