"""Compare Hessian provenance, allowing a nearby geometry for ReadFC.

ReadFC supplies an initial curvature estimate, not a continuation of the
optimizer trajectory. Keep atom order, state, Hamiltonian and solvent identical;
coordinates, MO guesses and optimization strategy may change.
"""

import re

from .calculation_cache import canonical_input


def recipe(backend, text):
    if backend == 'orca':
        text = re.sub(r'%geom\b.*?\bend\b', '', text, flags=re.I | re.S)
        text = re.sub(r'\b(?:Opt|LooseOpt|TightOpt|VeryTightOpt|NumFreq|Freq)\b', '', text, flags=re.I)
        text = re.sub(r'%moinp\s+"[^"]+"|\bMOREAD\b', '', text, flags=re.I)
        text = re.sub(r'^(\s*[A-Za-z]+)\s+[-+\d.eEdD]+\s+[-+\d.eEdD]+\s+[-+\d.eEdD]+\s*$',
                      r'\1', text, flags=re.M)
    elif backend == 'gamess':
        text = re.sub(r'\$(?:STATPT|FORCE|GUESS|VEC|HESS)\b.*?\$END', '', text, flags=re.I | re.S)
        text = re.sub(r'RUNTYP\s*=\s*\w+', 'RUNTYP=OPTIMIZE', text, flags=re.I)
        text = re.sub(r'(\$DATA\s*\n)[^\n]*', r'\1QCforever', text, flags=re.I)
        text = re.sub(r'^(\s*\S+\s+[-+\d.]+)\s+[-+\d.eEdD]+\s+[-+\d.eEdD]+\s+[-+\d.eEdD]+\s*$',
                      r'\1', text, flags=re.M)
    else:
        raise ValueError('Portable Hessian reads support ORCA and GAMESS')
    return canonical_input(backend, text)


def calculates_hessian(backend, text):
    if backend == 'orca':
        return bool(re.search(r'\b(?:Calc_Hess\s+true|Recalc_Hess\s+[1-9]\d*)', text, re.I))
    return bool(re.search(r'\bHESS\s*=\s*CALC', text, re.I))


def gamess_block(text):
    """Only accept a closed, nonempty $HESS block, choosing the last one."""
    blocks = re.findall(r'^\s*\$HESS\s*\n(.*?)\$END', text, re.I | re.M | re.S)
    if not blocks or not re.search(r'[-+]?\d*\.\d+(?:[DEde][-+]?\d+)?', blocks[-1]):
        raise ValueError('No complete GAMESS $HESS block was saved')
    return blocks[-1]
