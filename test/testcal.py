import os, sys, gc, csv, pickle
import numpy as np
from qcforever.gaussian_run import GaussianRunPack
				
usage ='Usage; %s' % sys.argv[0]

option = "symm pka opt homolumo energy dipole deen stable2o2 uv=/home/sumita/QCforever/Samples/UV_peak.dat aip aea tadf" 

infilename = 'ch2o.xyz'

test = GaussianRunPack.GaussianDFTRun('B3LYP', 'STO-3G', 10, option, infilename, solvent='water', restart=False, pklsave=True)

test.mem = '2GB'
test.timexe = 12*60*60
outdic = test.run_gaussian()

####for make a pkl file##########################
'''
result_dic = "reference.pkl"

with open(result_dic, 'wb') as f:
    pickle.dump(outdic, f)
'''
#################################################

ref_dic = "ch2o.pkl"

with open(ref_dic, 'rb') as f:
    ref_data = pickle.load(f)
    

def compare_list_ab(a, b):

    a = np.array(a)
    b = np.array(b)

    return max(abs(a - b))

print("Check symm: ")
if ref_data['symm'] == outdic['symm']:
    print("OK")
else: 
    print ("Something strange...")

print("Check ground state optimization: ")
if compare_list_ab(ref_data['GS_MaxDisplace'],  outdic['GS_MaxDisplace']) < 0.1:
    print("OK") 
else:
    print("Something strange...")

print("Check HOMO/LUMO gap")
if compare_list_ab(ref_data['homolumo'],  outdic['homolumo']) < 0.001:
    print("OK") 
else:
    print("Something strange...")

print("Check dipole")
if compare_list_ab(ref_data['dipole'],  outdic['dipole']) < 0.01:
    print("OK") 
else:
    print("Something strange...")

print("Check Energy")
if compare_list_ab(ref_data['Energy'],  outdic['Energy']) < 0.0001:
    print("OK") 
else:
    print("Something strange...")

print("Check deen")
if abs(ref_data['deen'] - outdic['deen']) < 0.0001:
    print("OK")
else: 
    print("Something strange...")

print("Check stable2o2")
if compare_list_ab(ref_data['stable2o2'], outdic['stable2o2']) < 0.001:
    print("OK")
else: 
    print("Something strange...")

print("Check vip")
if compare_list_ab(ref_data['vip'],  outdic['vip'])  < 0.001:
    print("OK")
else: 
    print("Something strange...")

print("Check aip")
if compare_list_ab(ref_data['aip'],  outdic['aip'])  < 0.001:
    print("OK")
else: 
    print("Something strange...")

print("Check vea")
if compare_list_ab(ref_data['vea'],  outdic['vea'])  < 0.001:
    print("OK")
else: 
    print("Something strange...")

print("Check vea")
if compare_list_ab(ref_data['aea'],  outdic['aea'])  < 0.001:
    print("OK")
else: 
    print("Something strange...")

print("Check relaxedIP_MaxDisplace")
if compare_list_ab(ref_data['relaxedIP_MaxDisplace'],  outdic['relaxedIP_MaxDisplace']) < 0.01:
    print("OK")
else: 
    print("Something strange...")

print("Check relaxedEA_MaxDisplace")
if compare_list_ab(ref_data['relaxedEA_MaxDisplace'],  outdic['relaxedEA_MaxDisplace']) < 0.01:
    print("OK")
else: 
    print("Something strange...")

print("Check cden")
if compare_list_ab(ref_data['cden'][1],  outdic['cden'][1]) < 0.001:
    print("OK")
else: 
    print("Something strange...")

print("Check uv wavelength")
if compare_list_ab(ref_data['uv'][0],  outdic['uv'][0]) < 0.01:
    print("OK")
else: 
    print("Something strange...")

print("Check UV intensity")
if compare_list_ab(ref_data['uv'][1],  outdic['uv'][1]) < 0.001:
    print("OK")
else: 
    print("Something strange...")

print("Check CD length")
if compare_list_ab(ref_data['uv'][2],  outdic['uv'][2]) < 0.001:
    print("OK")
else: 
    print("Something strange...")

print("Check CD intensity")
if compare_list_ab(ref_data['uv'][3],  outdic['uv'][3]) < 0.001:
    print("OK")
else: 
    print("Something strange...")

print('Similarity/Dissimilarity')
if compare_list_ab(ref_data['Similarity/Dissimilarity'],  outdic['Similarity/Dissimilarity']) < 0.01:
    print("OK")
else: 
    print("Something strange...")

print('Check MinEtarget')
if abs(ref_data['MinEtarget'] -  outdic['MinEtarget']) < 0.0001:
    print("OK")
else: 
    print("Something strange...")

print('Check Min_MaxDisplace')
if compare_list_ab(ref_data['Min_MaxDisplace'],  outdic['Min_MaxDisplace']) < 0.1:
    print("OK")
else: 
    print("Something strange...")

print("Check Fluorescence wavelength")
if compare_list_ab(ref_data['fluor'][0],  outdic['fluor'][0]) < 0.1:
    print("OK")
else: 
    print("Something strange...")

print("Check Fluorescence intensity")
if compare_list_ab(ref_data['fluor'][1],  outdic['fluor'][1]) < 0.01:
    print("OK")
else: 
    print("Something strange...")

print("Check CD fluorescence length")
if compare_list_ab(ref_data['fluor'][2],  outdic['fluor'][2]) < 0.01:
    print("OK")
else: 
    print("Something strange...")

print("Check CD fluorescence intensity")
if compare_list_ab(ref_data['fluor'][3],  outdic['fluor'][3]) < 0.01:
    print("OK")
else: 
    print("Something strange...")

print("Check T_Min")
if abs(ref_data['T_Min'] -  outdic['T_Min']) < 0.0001:
    print("OK")
else: 
    print("Something strange...")

print("Check T_Min_MaxDisplace")
if compare_list_ab(ref_data['T_Min_MaxDisplace'],  outdic['T_Min_MaxDisplace']) < 0.01:
    print("OK")
else: 
    print("Something strange...")

print("Check Phos wavelength")
if compare_list_ab(ref_data['T_Phos'][0],  outdic['T_Phos'][0]) < 0.01:
    print("OK")
else: 
    print("Something strange...")

print("Check Phos intensity")
if compare_list_ab(ref_data['T_Phos'][1],  outdic['T_Phos'][1]) < 0.01:
    print("OK")
else: 
    print("Something strange...")

print("Check CD Phos length")
if compare_list_ab(ref_data['T_Phos'][2],  outdic['T_Phos'][2]) < 0.01:
    print("OK")
else: 
    print("Something strange...")

print("Check CD Phos intensity")
if compare_list_ab(ref_data['T_Phos'][3],  outdic['T_Phos'][3]) < 0.01:
    print("OK")
else: 
    print("Something strange...")

print("Check Delta(S-T)")
if abs(ref_data['Delta(S-T)'] - outdic['Delta(S-T)']) < 0.01:
    print("OK")
else: 
    print("Something strange...")

print("Check pka")
if abs(ref_data['pka'] - outdic['pka']) < 0.001:
    print("OK")
else: 
    print("Something strange...")



