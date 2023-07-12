import glob
import os
import gc
import sys
import subprocess
import time


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
    count = 0
    for line in lines:
        if line.find("#") >= 0:
            count += 1
    return count


def count_Finishjob(jobname):
    Gau_logfile = f"{jobname}.log"
    with open(Gau_logfile,'r') as infile:
        lines = infile.readlines()
    count = 0
    for line in lines:
        if line.find("Normal termination") >= 0:
            count += 1
        if line.find("Error termination") >= 0:
            pass
    return count


def exe_Gaussian(jobname, exe_time):
    job_state = ""
    print(f"Time for Gaussian: {exe_time}")
    Njob = count_Njob(jobname)
    print(f'Number of Job: {Njob}')
    print(os.getcwd())
    os.environ["GAUSS_SCRDIR"] = os.getcwd()

    GaussianPros = subprocess.Popen(["g16", jobname])
    for _ in range(exe_time):
        time.sleep(1)
        if GaussianPros.poll() != None:
            break
    if GaussianPros.poll() == None:
        job_state = "timeout"
        print(GaussianPros.poll())
        GaussianPros.kill()
        for f in glob.glob("./Gau-*"):    
            print(f)
            os.remove(os.path.join('.', f))
    NFinishedJob = count_Finishjob(jobname)
    if Njob == NFinishedJob:
        job_state = "normal"
    elif job_state != "timeout" or Njob != NFinishedJob:
        job_state = "abnormal"

    print (GaussianPros)
    del GaussianPros
    gc.collect()

    return job_state 

if __name__ == '__main__':
    usage ='Usage; %s jobname' % sys.argv[0]

    try:
        jobname = sys.argv[1]
    except:
        print (usage); sys.exit()

    print('Number of jobs: ', count_Njob(jobname))
    print('Number of normally finished job', count_Finishjob(jobname))
