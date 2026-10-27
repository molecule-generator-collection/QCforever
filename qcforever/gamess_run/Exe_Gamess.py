import os
import sys
import shutil
import subprocess
from pathlib import Path

from qcforever.util import job_timeout
from qcforever.util.gamess_scratch import resolve_scratch_directories


def read_rungms(gmspath):
    return tuple(str(path) for path in resolve_scratch_directories(gmspath))


def exe_Gamess(jobname, gamessversion, nproc):
    job_state = ""
    cwd = os.getcwd()
    print(cwd)
    gmspath = shutil.which('rungms')
    if gmspath is None:
        raise FileNotFoundError("GAMESS launcher 'rungms' was not found on PATH")
    scr, userscr = read_rungms(gmspath)
    for directory in dict.fromkeys((Path(scr), Path(userscr))):
        directory.mkdir(parents=True, exist_ok=True)
        for path in directory.iterdir():
            if path.name.startswith((jobname + '.', jobname + '_')) and path.is_file():
                print(path)
                path.unlink()
    infile = jobname + '.inp'
    outfile = jobname +'.log'
    datfile = jobname +'.dat'

    environment = os.environ.copy()
    environment.update(SCR=scr, USERSCR=userscr)
    GamessRes = job_timeout.run(
        [gmspath, infile, gamessversion, str(nproc)],
        stdout=subprocess.PIPE,
        env=environment,
    )
    with open(outfile, 'w') as f:
        f.write(GamessRes.stdout.decode('utf-8'))

    path_datfile = userscr+'/'+datfile
    destination = os.path.join(cwd, datfile)
    if Path(path_datfile).resolve() != Path(destination).resolve():
        shutil.move(path_datfile, destination)

    return job_state 

if __name__ == '__main__':
    usage = f'Usage; {sys.argv[0]}'

    try:
        jobname = sys.argv[1]
    except:
        print (usage); sys.exit()

    print(exe_Gamess(jobname, '00', 2))

