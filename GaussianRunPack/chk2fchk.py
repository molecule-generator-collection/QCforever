import os, glob
import shutil
import subprocess


def Get_chklist():

    for i in glob.glob('./*.chk'):
        print(i)
        try:
            subprocess.run(['formchk', i])
        except:
            print ("Failed converting chk to fchk!")            

        os.remove('./'+ i)

if __name__ == '__main__':

    Get_chklist()
