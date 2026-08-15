import glob
import os

from qcforever.util import job_timeout


def Get_chklist(remove):
    for f in glob.glob('./*.chk'):
        print(f)
        try:
            job_timeout.run(['formchk', f], check=True)
        except:
            print("Failed converting chk to fchk!")            
        if remove == 1:
            os.remove(os.path.join('.', f))
        else:
            pass


if __name__ == '__main__':
    Get_chklist()
