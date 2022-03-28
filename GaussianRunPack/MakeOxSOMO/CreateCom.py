import os.path
import Utils

def Input(funct, basis, name):
  O2coord = \
  "O       0.000000    0.000000    0.000000 \n" +\
  "O       1.207500    0.000000    0.000000\n"

  input_s = name.split('.')
  chk_s = input_s[0]

  s = ""
  s = s + "%nprocshared=1\n"
  s = s + "%mem=1GB\n"
  s = s + "%chk=" + chk_s +"\n"
  s = s + "#p " + funct + "/" + basis + " scf=(conver=8,maxcycle=10000)\n"
  s = s + "Opt test\n"
  s = s + "\n"
  s = s + "Create O2.com\n"
  s = s + "\n"
  s = s + "0 3\n"
  s = s + O2coord
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
    name = fnc + "_" + bas + "_O2" + ".com"
    name = name.replace("*", "s")
    s = Input(fnc, bas, name)
    path = directory + name
    f = open(path, 'w')
    f.write(s)
    f.close()

