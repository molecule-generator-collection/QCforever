"""Load a completed Gaussian ground-state frequency calculation for stable."""

from pathlib import Path
import pickle
import re

import numpy as np

from qcforever.util import calculation_cache
from qcforever.util.imaginary_modes import StationaryPoint, same_geometry

from . import Get_FreqPro, parse_log


FREQUENCY_KEYS = ('freq', 'IR', 'Raman', 'Ezp', 'Et', 'E_enth', 'E_free', 'Ei', 'Cv', 'Si', 'freqmode')


def frequency_context(parsed, lines):
    end = next((i for i, line in enumerate(lines) if 'Frequencies --' in line), len(lines))
    header = lines[:end]
    if not any(' orientation:' in line for line in header):
        raise ValueError('Missing orientation in the Gaussian frequency step')
    atoms, *axes = parsed.extract_MolCoordlog(header)
    if not any(re.search(r'Charge\s*=\s*-?\d+\s+Multiplicity\s*=\s*\d+', line) for line in header):
        raise ValueError('Missing Gaussian frequency electronic state')
    charge, mult = parsed.extract_inputChargeSpin(header)
    return {
        'backend': 'gaussian', 'route': parsed.extract_method_text(header),
        'atoms': atoms, 'coordinates_angstrom': np.column_stack(axes).tolist(),
        'charge': int(charge), 'multiplicity': int(mult),
        'energy': parsed.Extract_SCFEnergy(lines),
        'pcm_parameters': re.findall(r'^\s*(?:Eps|EpsInf)\s*=\s*([+\-\d.]+)', '\n'.join(header), re.M | re.I),
    }


def _from_log(path):
    parsed = parse_log.parse_log(path)
    for link in reversed(parsed.SplitLinks()[1:]):
        lines = link.splitlines()
        route = parsed.extract_method_text(lines)
        if not re.search(r'\bfreq\b', route, re.I) or re.search(r'\b(?:td|tda|nac)\b', route, re.I):
            continue
        if not calculation_cache.successful('gaussian', route, link):
            raise ValueError(f'Incomplete Gaussian frequency step: {path}')
        values = Get_FreqPro.Extract_Freq(lines)
        result = dict(zip(FREQUENCY_KEYS[:-1], values))
        result['freqmode'] = list(Get_FreqPro.Extract_vibvec(lines))
        result['frequency_context'] = frequency_context(parsed, lines)
        return result
    raise ValueError(f'No ground-state frequency step in {path}')


def load(path, name, *, atoms, coordinates, charge, mult, method, solvent_route, solvent_read,
         parameter_route, companion=None, allow_completed_stable=False):
    """Return a validated stationary point in the frequency calculation's frame."""
    path = Path(path)
    if path.suffix.lower() == '.pkl':
        with path.open('rb') as stream:
            result = pickle.load(stream)
        result = result.get('s0_vibration', result)
        log = result.get('gaussian_log')
        adjacent = path.with_suffix('.log')
        seed_log = path.parent / log if log else (
            adjacent if adjacent.is_file() or companion is None else Path(companion))
        if not result.get('frequency_context'):
            result = _from_log(seed_log)
        if seed_log.is_file():
            path = seed_log
    else:
        result = _from_log(path)
    context = result['frequency_context']
    route = re.sub(r'\s+', '', context['route']).lower()
    # Strip only the controls that do not change the potential energy surface.
    if not re.match(r'^#(?:[pn])?' + re.escape(method.lower()) + r'(?=\s|$)',
                    context['route'].lower()):
        raise ValueError('Gaussian frequency method or basis differs')
    expected_solvent = re.sub(r'\s+', '', solvent_route).lower()
    actual_solvent = re.search(r'scrf=\([^)]*\)', route)
    if (actual_solvent[0] if actual_solvent else '') != expected_solvent:
        raise ValueError('Gaussian frequency solvent differs')
    if solvent_read:
        expected_dielectric = re.findall(r'(?:eps|epsinf)\s*=\s*([+\-\d.]+)', solvent_read, re.I)
        actual_dielectric = context.get('pcm_parameters', [])
        if not actual_dielectric or list(map(float, expected_dielectric)) != list(map(float, actual_dielectric)):
            raise ValueError('Gaussian frequency dielectric differs')
    expected_parameter = re.sub(r'\s+', '', parameter_route).lower()
    actual_parameter = ''.join(re.findall(r'iop\([^)]*\)', route))
    if actual_parameter != expected_parameter:
        raise ValueError('Gaussian frequency functional parameters differ')
    # QCforever uses Gaussian defaults for the grid and SCF settings.
    controls = re.sub(r'^#(?:[pn])?' + re.escape(method.lower()), '', route, count=1)
    controls = re.sub(r'(?:scrf|iop)=?\([^)]*\)|freq(?:=\([^)]*\)|=[a-z]+)?|'
                      r'opt(?:=\([^)]*\)|=[a-z]+)?|geom=allcheck|guess=read|test', '', controls)
    if controls or not re.search(r'\bfreq\b', context['route'], re.I):
        raise ValueError('Gaussian frequency route has incompatible or unverified settings')
    if (context['backend'], context['charge'], context['multiplicity']) != ('gaussian', charge, mult):
        raise ValueError('Gaussian frequency electronic state differs')
    xyz = np.asarray(context['coordinates_angstrom'])
    previous_start = result.get('stable_input_geometry') if allow_completed_stable and result.get('stable_0') else None
    if (list(atoms) != context['atoms'] or not (
            same_geometry(coordinates, xyz) or (
                previous_start is not None and same_geometry(coordinates, previous_start)))):
        raise ValueError('Gaussian frequency geometry differs from the input structure')
    props = {key: result[key] for key in FREQUENCY_KEYS}
    energy = context['energy']
    props.update(Energy=energy, frequency_context=context, stable_source=str(path.resolve()))
    point = StationaryPoint(name, list(atoms), xyz, energy[0], props['freq'], props['freqmode'], props)
    if (not np.isfinite(point.energy) or len(point.frequencies) == 0
            or np.shape(point.modes) != (3, len(point.frequencies), len(atoms))
            or not np.isfinite(np.asarray(point.modes, dtype=float)).all()):
        raise ValueError('Invalid Gaussian frequency energy or frequencies')
    point.minimum_frequency
    return point
