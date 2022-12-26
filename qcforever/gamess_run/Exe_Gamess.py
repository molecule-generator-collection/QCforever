import glob
import os
import sys
import subprocess
import shutil


def exe_Gamess(jobname, gamessversion, nproc):
    job_state = ""
    cwd = os.getcwd()
    print(cwd)
    scr = os.getenv('SCR')
    userscr = os.getenv('USERSCR')
    job_intfile = scr+'/'+jobname+'*'
    userjob_intfile = userscr+'/'+jobname+'*'
    for f in glob.glob(job_intfile):
        print(f)
        os.remove(os.path.join(job_intfile, f))
    for f in glob.glob(userjob_intfile):
        print(f)
        os.remove(os.path.join(userjob_intfile, f))
    infile = jobname + '.inp'
    outfile = jobname +'.log'
    datfile = jobname +'.dat'

    GamessRes = subprocess.run(['rungms', infile, gamessversion, str(nproc)], stdout=subprocess.PIPE)
    with open(outfile, 'w') as f:
        f.write(GamessRes.stdout.decode('utf-8'))

    path_datfile = userscr+'/'+datfile
    shutil.move(path_datfile, cwd)

    return job_state 

if __name__ == '__main__':
    usage = f'Usage; {sys.argv[0]}'

    try:
        jobname = sys.argv[1]
    except:
        print (usage); sys.exit()

    print(exe_Gamess(jobname, '00', 2))


