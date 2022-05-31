import os.path
import Utils

def Input(funct, basis, name):
  TMScoord = \
  "Si       0.000000    0.000000   -0.000000 \n" +\
  "C       -0.000000   -0.000000    1.879889 \n" +\
  "H       -0.884987   -0.510947    2.279340 \n" +\
  "H        0.884987   -0.510947    2.279340 \n" +\
  "H       -0.000000    1.021895    2.279340 \n" +\
  "C       -0.000000   -1.772377   -0.626630 \n" +\
  "H       -0.884987   -2.319298   -0.278054 \n" +\
  "H        0.000000   -1.808350   -1.723232 \n" +\
  "H        0.884987   -2.319298   -0.278054 \n" +\
  "C        1.534923    0.886188   -0.626630 \n" +\
  "H        1.566077    0.904175   -1.723232 \n" +\
  "H        1.566077    1.926070   -0.278054 \n" +\
  "H        2.451064    0.393228   -0.278054 \n" +\
  "C       -1.534923    0.886188   -0.626630 \n" +\
  "H       -1.566077    1.926070   -0.278054 \n" +\
  "H       -1.566077    0.904175   -1.723232 \n" +\
  "H       -2.451064    0.393228   -0.278054 \n"

  input_s = name.split('.')
  chk_s = input_s[0]

  s = ""
  s = s + "%nprocshared=1\n"
  s = s + "%mem=1GB\n"
  s = s + "%chk=" + chk_s +"\n"
  s = s + "#p " + funct + "/" + basis + " scf=(conver=8,maxcycle=10000)\n"
  s = s + "Opt \n"
  s = s + "\n"
  s = s + "Create TMS.com\n"
  s = s + "\n"
  s = s + "0 1\n"
  s = s + TMScoord
  s = s + "\n"

  s = s + "--Link1--\n"
  s = s + "%nprocshared=1\n"
  s = s + "%chk=" + chk_s +"\n"
  s = s + "#p " + funct + "/" + basis + "\n"
  s = s + "nmr \n"
  s = s +  'Geom=AllCheck Guess=Read'
  s = s + "\n"
  s = s + "\n"

  return s

directory = "work/"

if os.path.exists(directory):
  print(directory + " already exists")
  print("please remove " + directory)
  exit(1)
else:
  os.mkdir(directory)

for fnc in Utils.Functional:
  for bas in Utils.Basis:
    name = fnc + "_" + bas + "_" + "TMS" + ".com"
    name = name.replace("*", "s")
    s = Input(fnc, bas, name)
    path = directory + name
    f = open(path, 'w')
    f.write(s)
    f.close()

