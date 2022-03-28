import os, sys 
import gc
import subprocess
import time

def count_Njob(jobname):

    Gau_infile = jobname +'.com'

    with open(Gau_infile,'r') as infile:
        lines = infile.readlines()

    count = 0

    for line in lines:
        if line.find("#") >= 0:
            count += 1

    return count


def count_Finishjob(jobname):

    Gau_logfile = jobname +'.log'

    with open(Gau_logfile,'r') as infile:
        lines = infile.readlines()

    count = 0

    for line in lines:
        if line.find("Normal termination") >= 0:
            count += 1
        if line.find("Error termination") >= 0:
            count = count_Njob(jobname) 

    return count



def exe_Gaussian(jobname, exe_time):

    print ("Time for Gaussian: ", exe_time)

    Njob = count_Njob(jobname)
    print('Number of Job: ', Njob)

    logfile = jobname+'.log'

    #GaussianPros = subprocess.run(["g16", jobname])
    #For shell=True############################################
    #cmd = "g16"+" " + PreGauInput[0]
    #GaussianPros = subprocess.run(cmd, shell=True)
    ##########################################################
    #GaussianPros = subprocess.Popen(["g16", jobname])
    #for i in range(exe_time):
    #    time.sleep(1)
    #    if os.path.exists('result_dic.pickle'):
    #        break
    #GaussianResult.kill()

    GaussianPros = subprocess.Popen(["g16", jobname])
    for i in range(exe_time):
        time.sleep(1)
        NFinishedJob = count_Finishjob(jobname)
        if NFinishedJob == Njob:
            break

        if i == exe_time-1:
            print ("End of waiting for Gaussian!")
    
    GaussianPros.kill()
    
    print (GaussianPros)
    del GaussianPros
    gc.collect()

