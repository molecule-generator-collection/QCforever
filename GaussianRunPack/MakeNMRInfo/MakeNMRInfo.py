import glob
import re
import json
import Utils

def Get_NMR(path):
  sl = path.split("/")
  s1 = sl[1]
  sl = s1.split("_")
  fnc = sl[0]
  bas = sl[1].replace("s", "*")

  Element = []
  ab_ppm = []

  f = open(path, 'r')
  for line in f:
    if line.find("Isotropic =  ") >= 0:
      line_Info = line.split()
      #print (line_Info[1])
      if line_Info[1] in Element:
        continue
      else:
        Element.append(line_Info[1])
      #print(line_Info[4])
        ab_ppm.append(float(line_Info[4]))
  f.close()

#  print(Element)
#  print(ab_ppm)

  return fnc, bas, Element, ab_ppm 

files = glob.glob("work/*.log")

Elements = Utils.TMS
nmrlist = {}

for fnc in Utils.Functional:
  nmrlist[fnc] = {}
  for bas in Utils.Basis:
    nmrlist[fnc][bas] = {}
    for elm in Elements:
      nmrlist[fnc][bas][elm] = "none"
        
for path in files:
#  print(path)
  fnc, bas, elem, TMS_ppm = Get_NMR(path)
#  print(fnc, bas, elem, TMS_ppm)
  for i in range(len(elem)):
    nmrlist[fnc][bas][elem[i]] = TMS_ppm[i]

print(json.dumps(nmrlist,indent=4))
