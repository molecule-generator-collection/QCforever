import os, glob
import shutil
import subprocess


def Get_chklist():

    for i in glob.glob('./*.chk'):
        print(i)
        subprocess.run(['formchk', i])
        os.remove('./'+ i)

if __name__ == '__main__':

    Get_chklist()
