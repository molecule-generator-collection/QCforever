#!/usr/bin/env python3
"""Small process-level test double; no chemistry or native xTB validation."""
import os
from pathlib import Path
import sys
import time

xyz = Path('input.xyz').read_text().splitlines()
atoms = int(xyz[0])
mode = os.environ.get('QCFOREVER_FAKE_XTB_MODE', 'normal')
for cycle in range(1, 7):
    energy = -10 - cycle / 100
    with Path('xtbopt.log').open('a') as stream:
        stream.write(str(atoms)+'\nenergy: '+str(energy)+' gnorm: 0.1\n'+'\n'.join(xyz[2:])+'\n')
    print('CYCLE '+str(cycle), flush=True)
    print('molecular gradient ...', flush=True)
    if cycle == 1:
        print('molecular gradient ...', flush=True)
    print(' 1 -10.0 1.0E-5', flush=True)
    print('* total energy : '+str(energy), flush=True)
    print('gradient norm : 0.100000', flush=True)
    time.sleep(.04)
    if mode == 'hang':
        time.sleep(30)
if mode == 'fail':
    sys.exit(1)
Path('xtbopt.xyz').write_text('\n'.join(xyz)+'\n')
print('GEOMETRY OPTIMIZATION CONVERGED AFTER 6 CYCLES', flush=True)
print('TOTAL ENERGY -10.06', flush=True)
