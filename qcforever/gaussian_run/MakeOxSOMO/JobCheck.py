import os, sys 
import gc
import time

exe_time = 60*60*60

def Judge_job(logfile):

    finish = False

    with open(logfile,'r') as infile:
        lines = infile.readlines()

    for line in lines:
        if line.find("termination") >= 0:
            finish = True

    return finish 

def Job_Checker(logfile):

    for i in range(exe_time):
        time.sleep(1)
        finish = Judge_job(logfile)
        if finish :
           # print ("The job has finished!")
            break
    

if __name__ == '__main__':

    usage = 'Usage; %s jobname' % sys.argv[0]

    try:
        infilename = sys. argv[1]
    except:
        print (usage); sys.exit()

    Job_Checker(infilename)



