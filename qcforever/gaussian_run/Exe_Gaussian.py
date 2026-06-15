import glob
import os
import gc
import sys
import subprocess
import time

from pathlib import Path

from qcforever import gaussian_run
from qcforever.util import ConvergenceJudge


def get_pid(jobname):
    Gau_logfile = f"{jobname}.log"
    with open(Gau_logfile,'r') as infile:
        lines = infile.readlines()
    count = 0
    Gau_pid = []
    for line in lines:
        if line.find("Entering Link 1") >= 0:
            count += 1
            pid_info = line.split()
            Gau_pid.append(int(pid_info[-1]))
    Current_PID = Gau_pid[-1]
    return Current_PID


def count_Njob(jobname):
    Gau_infile = f"{jobname}.com"
    with open(Gau_infile,'r') as infile:
        lines = infile.readlines()
    count = 1
    for line in lines:
        if line.find("--Link1--") >= 0:
            count += 1
    return count


def count_Finishjob(jobname):
    Gau_logfile = f"{jobname}.log"
    with open(Gau_logfile,'r') as infile:
        lines = infile.readlines()
    fcount = 0
    is_error = 0
    for line in lines:
        if line.find(r"Normal termination") >= 0:
            fcount += 1
        if line.find(r"Error termination") >= 0:
            print("Error is detected!")
            is_error += 1
            pass

    return fcount, is_error


def job_termination(proc):
    proc.kill()
    for f in glob.glob("./Gau-*"):    
        print(f)
        Path(f).unlink(missing_ok=True)


def exe_Gaussian(jobname, exe_time, error=0):
    job_state = ""
    print(f"Time for Gaussian: {exe_time}")
    Njob = count_Njob(jobname)
    print(f'Number of Job: {Njob}')
    print(os.getcwd())
    os.environ["GAUSS_SCRDIR"] = os.getcwd()

    print (f'Error monitering level: {error}')
    
    if error > 0:
        print (f'Terminate the job if the optimization is unlikely to converge.')
    elif error == 0:
        print(f'Any action will not be performed!')

    GaussianPros = subprocess.Popen(["g16", jobname])

    for sec in range(exe_time):
        time.sleep(1)

        if GaussianPros.poll() != None:
            break

        if error == 1 and sec > 0 and sec % 150 == 0:
            try:
                tmp_log = gaussian_run.MoniteringEnergy.Monitor_log(jobname+".log") 
                Index_TState, Eprofile, OSprofile = tmp_log.Extract_ConversionProcess()
                Values_ForceDisp = tmp_log.Extract_ConversionCriterion()
                #print(Values_ForceDisp)
                if len(Eprofile) > 10: 
                    #
                    #For Judge convergence possibility
                    #
                    Judge = ConvergenceJudge.ConvergenceJudge(Index_TState, Eprofile, Values_ForceDisp)
                    Score_Judge = Judge.judge()
        
                    if Score_Judge[1] < 0.5:
                        print ("No hope for geometry optimization! Job will be killed")
                        job_state = f'{Score_Judge[0]}' 
                        job_termination(GaussianPros)
                    else:
                        pass
                else:
                    pass

            except Exception as e:
                print(f'Monitoring failed: {e}')
                pass

    if GaussianPros.poll() == None:
        job_state = "timeout"
        print ("Timeout! Job will be killed")
        print(GaussianPros.poll())
        job_termination(GaussianPros)

    NFinishedJob, is_error = count_Finishjob(jobname)
    print (f'Success Job: {NFinishedJob} Error Job: {is_error}')
    if Njob == NFinishedJob and is_error == 0:
        job_state = "normal"
    elif job_state != "timeout" or is_error != 0:
        if job_state == '':
            job_state = "abnormal"
        else:
            pass

    print (GaussianPros)
    del GaussianPros

    return job_state 

if __name__ == '__main__':
    usage ='Usage; %s jobname' % sys.argv[0]

    try:
        jobname = sys.argv[1]
    except:
        print (usage); sys.exit()

    print('Number of jobs: ', count_Njob(jobname))
    print('Number of normally finished job', count_Finishjob(jobname))
