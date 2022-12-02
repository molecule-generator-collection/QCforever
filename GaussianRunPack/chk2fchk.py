import glob
import os
import subprocess


def Get_chklist():
    for f in glob.glob('./*.chk'):
        print(f)
        try:
            subprocess.run(['formchk', f])
        except:
            print("Failed converting chk to fchk!")            
        os.remove(os.path.join('.', f))


if __name__ == '__main__':
    Get_chklist()
