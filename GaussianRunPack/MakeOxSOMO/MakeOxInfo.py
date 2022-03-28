import glob
import re
import json
import Utils
import JobCheck

def Get_SOMO(path):

    JobCheck.Job_Checker(path)

    sl = path.split("/")
    s1 = sl[1]
    sl = s1.split("_")
    fnc = sl[0]
    bas = sl[1].replace("s", "*")

    AlphaEigenVal = []
    BetaEigenVal = []

    f = open(path, 'r')
    for line in f:
        if line.find(" basis functions, ") >=0:
            line_StateInfo = line.split()
        #    print (line_StateInfo[0], flush=True)
            
            NumBasisFunc = int(line_StateInfo[0])
  
        if line.find(" alpha electrons ") >=0:
            line_StateInfo = line.split()
        #    print (line_StateInfo[0], flush=True)
        #    print (line_StateInfo[3], flush=True)
    
            NumAlphaElec = int(line_StateInfo[0])
            NumBetaElec = int(line_StateInfo[3])
        #    print ("# of alpha : ", NumAlphaElec, " beta : ", NumBetaElec, flush=True)

        if line.find("Alpha  occ. eigenvalues --") >=0:
            if len(AlphaEigenVal) == NumBasisFunc:
                AlphaEigenVal=[]
            
            line_removed = line.replace("Alpha  occ. eigenvalues --", " ")
            line_StateInfo = line_removed.split()
            for i in range(len(line_StateInfo)):
                AlphaEigenVal.append(float(line_StateInfo[i]))
  
        if line.find("Alpha virt. eigenvalues --") >=0:
            line_removed = line.replace("Alpha virt. eigenvalues --", " ")
            line_StateInfo = line_removed.split()
            for i in range(len(line_StateInfo)):
                AlphaEigenVal.append(float(line_StateInfo[i]))
  
  
        if line.find("Beta  occ. eigenvalues --") >=0:
  
            if len(BetaEigenVal) == NumBasisFunc:
                BetaEigenVal=[]
            
            line_removed = line.replace("Beta  occ. eigenvalues --", " ")
            line_StateInfo = line_removed.split()
            for i in range(len(line_StateInfo)):
                BetaEigenVal.append(float(line_StateInfo[i]))
  
        if line.find("Beta virt. eigenvalues --") >=0:
            line_removed = line.replace("Beta virt. eigenvalues --", " ")
            line_StateInfo = line_removed.split()
            for i in range(len(line_StateInfo)):
                BetaEigenVal.append(float(line_StateInfo[i]))

    f.close()

    SOMO_a =  AlphaEigenVal[NumAlphaElec-1]
    LUMO_b =  BetaEigenVal[NumBetaElec]

    return fnc, bas, SOMO_a, LUMO_b

files = glob.glob("work/*.log")

Oxlist = {}

for fnc in Utils.Functional:
    Oxlist[fnc] = {}
    for bas in Utils.Basis:
        Oxlist[fnc][bas] = {}
        
for path in files:
#  print(path)
    fnc, bas, SOMO_a, LUMO_b = Get_SOMO(path)
#  print(fnc, bas, SOMO_a, SOMO_b)
    Oxlist[fnc][bas]['SOMO_a'] = SOMO_a
    Oxlist[fnc][bas]['LUMO_b'] = LUMO_b

print(json.dumps(Oxlist,indent=4))
