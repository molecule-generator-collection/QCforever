import sys
import os

def respec_cores(spec_cores):
    ava_cores = len(os.sched_getaffinity(0)) 
    if ava_cores < spec_cores:
        respec_cores = ava_cores
    else:
        respec_cores = spec_cores 

    return respec_cores

if __name__ == '__main__':
    usage = f'Usage; {sys.argv[0]} number_of_cores'

    try:
        spec_cores = int(sys.argv[1])
    except:
        print(usage)

    print(respec_cores(spec_cores))
        

