import os.path
import Utils

def Input(funct, basis, symbl, multi):
  s = ""
  s = s + "%nprocshared=1\n"
  s = s + "%mem=1GB\n"
  s = s + "#p " + funct + "/" + basis + " nosym scf=(conver=7,maxcycle=10000)\n"
  s = s + "\n"
  s = s + "Create com\n"
  s = s + "\n"
  s = s + "0 " + multi + "\n"
  s = s + symbl + "\n"
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
    for elm in Utils.Multiplicity:
      mul = str(Utils.Multiplicity[elm])
      s = Input(fnc, bas, elm, mul)
      name = fnc + "_" + bas + "_" + elm + ".com"
      name = name.replace("*", "s")
#      print("-----------------------")
#      print(name)
#      print("-----------------------")
#      print(s)
      path = directory + name
      f = open(path, 'w')
      f.write(s)
      f.close()

