import glob
import re
import json
import Utils

def GetEnergy(path):
  sl = path.split("/")
  s1 = sl[1]
  sl = s1.split("_")
  fnc = sl[0]
  bas = sl[1].replace("s", "*")
  elm = sl[2].replace(".log", "")
  eng = "none"
  f = open(path, 'r')
  for line in f:
    if line.find("SCF Done") >= 0:
      sl = re.split("\s+", line)
      eng = sl[5]
  f.close()
  return fnc, bas, elm, eng

files = glob.glob("work/*.log")

Elements = Utils.Multiplicity.keys()
englist = {}

for fnc in Utils.Functional:
  englist[fnc] = {}
  for bas in Utils.Basis:
    englist[fnc][bas] = {}
    for elm in Elements:
      englist[fnc][bas][elm] = "none"
        
for path in files:
  #print(path)
  fnc, bas, elm, eng = GetEnergy(path)
  #print(fnc, bas, elm, eng)
  englist[fnc][bas][elm] = eng

print(json.dumps(englist,indent=4))
