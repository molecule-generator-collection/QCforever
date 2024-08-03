import os
import re
import sys
import glob
import shutil
import subprocess


def read_rungms(gmspath):

    rungms_line = open(gmspath)
    lines = rungms_line.readlines()
    rungms_line.close()

    for ii in range(len(lines)):
        ll = lines[ii]
        if(re.match("set SCR",ll)):
            scr_list = ll.split('=')
            scrpath = scr_list[-1].replace("\n", "")
        
        if(re.match("set USERSCR",ll)):
            userscr_list = ll.split('=')
            userscrpath = userscr_list[-1].replace("\n", "")

    return scrpath, userscrpath


def exe_Gamess(jobname, gamessversion, nproc):
    job_state = ""
    cwd = os.getcwd()
    print(cwd)
    gmspath = shutil.which('rungms')
    scr, userscr = read_rungms(gmspath)
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
    shutil.move(path_datfile, os.path.join(cwd, datfile))

    return job_state 

if __name__ == '__main__':
    usage = f'Usage; {sys.argv[0]}'

    try:
        jobname = sys.argv[1]
    except:
        print (usage); sys.exit()

    print(exe_Gamess(jobname, '00', 2))


