import sys
import subprocess
import re


def ChargeSpin_fchk(fchkfile):
    TotalCharge = 0
    SpinMulti = 1
    with open(fchkfile, 'r') as fchk:
        lines = fchk.readlines()
    for line in lines:
        if (re.match("Charge", line)):
            Charge_inf = line.split()
            TotalCharge = int(Charge_inf[-1])
        if (re.match("Multiplicity", line)):
            Multi_inf = line.split()
            SpinMulti = int(Multi_inf[-1])
    lines = ""
    return TotalCharge, SpinMulti


def Get_fchk(jobname):
    fchkfile = f"{jobname}.fchk"
    TotalCharge, SpinMulti = ChargeSpin_fchk(fchkfile)
    try:
        subprocess.run(['unfchk', fchkfile], check=True)
    except:
        print ("Failed converting fchk to chk!")            
    return TotalCharge, SpinMulti


if __name__ == '__main__':
    usage = 'Usage; %s jobname' % sys.argv[0]
    try:
        jobname = sys.argv[1]
    except:
        print (usage); sys.exit()
    Charge, SpinMulti = Get_fchk(jobname)
    print(f"Charge: {Charge}")
    print(f"SpinMulti: {SpinMulti}")
