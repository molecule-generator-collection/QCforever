import sys

from qcforever.util import job_timeout
import re
from pathlib import Path


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


def convert(fchkfile, chkfile):
    """Convert into the job directory, propagating conversion failures."""
    destination = Path(chkfile)
    destination.unlink(missing_ok=True)
    job_timeout.run(['unfchk', str(fchkfile), str(destination)], check=True)
    if not destination.is_file() or not destination.stat().st_size:
        raise RuntimeError(f'unfchk did not produce {destination}')


def Get_fchk(jobname):
    fchkfile = f"{jobname}.fchk"
    TotalCharge, SpinMulti = ChargeSpin_fchk(fchkfile)
    convert(fchkfile, f'{jobname}.chk')
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
