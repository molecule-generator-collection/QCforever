import sys
import os
import time

def respec_cores(spec_cores):

    ava_cores = get_ava_cores()

    if ava_cores < spec_cores:
        respec_cores = ava_cores
        print(f'{spec_cores} cores you specifed are not available! \n'
                f'The number of cores for QC is reduced to {ava_cores}.')
    else:
        respec_cores = spec_cores 

    return respec_cores


def get_ava_cores():

    ava_cores = len(os.sched_getaffinity(0)) 

    if ava_cores == 0:
        while ava_cores == 0:
            time.sleep(10)
            print ('Waiting resource.....')
            ava_cores = len(os.sched_getaffinity(0)) 

    return ava_cores


if __name__ == '__main__':
    usage = f'Usage; {sys.argv[0]} number_of_cores'

    try:
        spec_cores = int(sys.argv[1])
    except:
        print(usage)

    print(respec_cores(spec_cores))
        

