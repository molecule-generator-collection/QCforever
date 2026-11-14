"""Tokenization and per-state optimization options for QCforever backends."""

import re


def requires_hessian(*fc_controls):
    """Does any optimization explicitly compute, rather than read, a Hessian?"""
    return any(re.search(r'\b(?:calcfc|calcall|recalcfc)\b', str(control), re.I)
               for control in fc_controls)


def split_options(text):
    """Keep recognized parenthesized option groups together, including spaces.

    Spectrum paths and all other existing scalar tokens retain their spelling.
    The caller expands optconf first.
    """
    result = []
    position = 0
    while position < len(text):
        if text[position].isspace():
            position += 1
            continue
        start = position
        group = re.match(r'(fluor|tadf|opt|cden|dipole)\s*=\s*\(', text[start:], re.I)
        if group:
            opening = start + group.end() - 1
            position, depth = opening + 1, 1
            while position < len(text) and depth:
                depth += (text[position] == '(') - (text[position] == ')')
                position += 1
            if depth or (position < len(text) and not text[position].isspace()):
                raise ValueError(f'Malformed {group[1]} group: check parentheses and spaces between options')
            result.append(group[1].lower() + '=' + text[opening:position])
        else:
            while position < len(text) and not text[position].isspace():
                position += 1
            result.append(text[start:position])
    return result


def cden_fukui(token):
    """Validate a cden option and return whether it requests Fukui indices."""
    _, separator, argument = token.partition('=')
    if not separator:
        return False
    if not re.fullmatch(r'\(\s*fukui\s*\)', argument, re.I):
        raise ValueError('Use cden or cden=(fukui)')
    return True


def parse_excited_dipole(token):
    """Read the target singlet root and geometry for a Gaussian dipole job."""
    match = re.fullmatch(r'dipole=\(([^()]*)\)', token, re.I)
    if not match:
        raise ValueError('Use dipole=(state=1, geom=input) or dipole=(state=1, geom=fluor)')
    fields = {}
    for item in match[1].lower().split(','):
        key, _, value = (part.strip() for part in item.partition('='))
        if key not in ('state', 'geom') or not value or key in fields:
            raise ValueError('Specify state and geom at most once in dipole=(...)')
        fields[key] = value
    state = fields.get('state', '1')
    geometry = fields.get('geom', 'input')
    if not re.fullmatch(r'[1-9]\d*', state):
        raise ValueError('Excited dipole state must be a positive integer')
    if geometry not in ('input', 'fluor'):
        raise ValueError('Excited dipole geom must be input or fluor')
    return int(state), geometry


def parse_excited_optimization(token, backend='gaussian'):
    """Read fluor[=N] / fluor=(FC=CalcFC,state=N), likewise for tadf."""
    key, separator, argument = token.partition('=')
    key = key.strip().lower()
    if key not in ('fluor', 'tadf'):
        raise ValueError('Expected fluor or tadf optimization options')
    values = {}
    if separator:
        argument = argument.strip()
        if argument.startswith('(') and argument.endswith(')'):
            body = argument[1:-1].strip()
            for item in body.split(',') if body else []:
                field, equal, value = item.partition('=')
                field, value = field.strip().lower(), value.strip()
                if not equal or not value or field not in ('fc', 'state'):
                    raise ValueError(f'Use {key}=(FC=CalcFC, state=1); supported fields are FC and state')
                if field in values:
                    raise ValueError(f'Duplicate {key} field: {field}')
                values[field] = value
        else:
            values['state'] = argument
    state_text = values.get('state', '1')
    if not re.fullmatch(r'[+-]?\d+', state_text):
        raise ValueError(f'{key} state must be an integer')
    state = int(state_text)
    allow_auto = backend == 'gaussian' and key == 'fluor'
    if state < 1 and not (allow_auto and state == -1):
        raise ValueError(f'{key} state must be positive' + (' or -1 (brightest state)' if allow_auto else ''))
    fc = values.get('fc', 'default').lower()
    choices = {'default': '', 'none': '', 'calcfc': 'CalcFC', 'calcall': 'CalcAll',
               'readfc': 'ReadFC', 'rcfc': 'RCFC'}
    if fc in choices:
        fc = choices[fc]
    elif re.fullmatch(r'recalcfc\s*=\s*[1-9]\d*', fc):
        fc = 'RecalcFC=' + fc.split('=', 1)[1].strip()
    else:
        raise ValueError(f'Unsupported {key} FC value: {values.get("fc")}; '
                         'use CalcFC, CalcAll, ReadFC, RCFC, RecalcFC=N, default or none')
    if fc == 'RCFC' and backend != 'gaussian':
        raise ValueError(f'{backend} does not support Gaussian RCFC; use ReadFC')
    return {'state': state, 'fc': fc}


def ground_optimization(argument, backend):
    """Native convergence level plus portable Hessian controls for opt=... ."""
    level, fc = 'normal', ''
    seen_level, seen_fc = False, False
    body = argument.strip().strip('()')
    for item in body.split(',') if body else []:
        item = item.strip()
        if item.lower() in ('normal', 'loose', 'tight', 'verytight'):
            if seen_level:
                raise ValueError('Specify optimization convergence only once')
            seen_level = True
            if backend == 'gamess' and item.lower() != 'normal':
                raise ValueError('GAMESS opt supports FC controls and normal convergence only')
            level = item.lower()
        else:
            if item.lower().startswith('fc='):
                item = item[3:].strip()
            if seen_fc:
                raise ValueError('Specify ground-state FC only once')
            seen_fc = True
            fc = parse_excited_optimization('fluor=(FC=' + item + ')', backend)['fc']
    return level, fc


def excited_options(tokens, backend='gaussian'):
    """Validate once and keep explicit fluor/tadf adjacent in S/T order."""
    specs, excited = {}, []
    first = None
    remaining = []
    for token in tokens:
        key = token.partition('=')[0].lower()
        if key in ('fluor', 'tadf'):
            if key in specs:
                raise ValueError(f'Specify {key} only once')
            specs[key] = parse_excited_optimization(token, backend)
            excited.append((key, token))
            if first is None:
                first = len(remaining)
        else:
            remaining.append(token)
    if first is not None:
        remaining[first:first] = [token for key in ('fluor', 'tadf') for k, token in excited if k == key]
    return remaining, specs
