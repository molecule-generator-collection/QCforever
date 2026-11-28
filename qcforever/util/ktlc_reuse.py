"""Persist the selected KTLC parameter separately from the trial calculations."""

import hashlib
import json
import math
from pathlib import Path

from . import calculation_cache


def request_key(owner, session):
    """Identify a complete request that has no externally configured dependencies."""
    options = str(owner.value)
    if (owner.functional.upper().endswith('-PRED')
            or any(option in options.lower() for option in (
                'optconf', 'readfc', 'rcfc', 'freq=', 'nmr=', 'uv='))
            or any(getattr(owner, name, None) for name in (
                'stable_from', 'conformer_config', 'ref_spectrum_paths', 'ref_uv_path'))):
        return None
    try:
        source = calculation_cache._hash(owner.in_file)
    except OSError:
        return None
    # These controls can change properties without changing the BO objective.
    names = ('functional', 'basis', 'solvent', 'SpecTotalCharge', 'SpecSpinMulti',
             'geom_spec', 'para_functional', 'ion_multiplicities', 'nstates',
             'maxiter', 'nmr_reference', 'fchk_angular', 'restart')
    identity = dict(schema=1, source=source, options=options, backend=session.backend,
                    variant=session.variant,
                    program=calculation_cache._program(session.backend, getattr(owner, 'executable', None)),
                    settings={name: repr(getattr(owner, name, None)) for name in names})
    if session.backend == 'gaussian' and 'nac' in options.lower().split():
        identity['nac_separate_job'] = True
    return hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()


def unchanged_result(session):
    if not session.enabled or not session.request_key:
        return None
    directory = Path.cwd() / session.ktlc_job
    if directory.is_symlink() or (directory / calculation_cache.DIRECTORY).exists():
        return None
    result = calculation_cache._read_result(directory / (session.ktlc_job + '.pkl'))
    saved = result.get('ktlc_reuse', {})
    if (result.get('log') != 'normal' or saved.get('request_key') != session.request_key
            or not valid_selection(session.backend, result.get('functional_parameter_optimization'))):
        return None
    steps = [dict(step, action='reused', reason='matching completed KTLC request in the result pickle')
             for step in result['reuse_report']['steps']]
    result['reuse_report'] = dict(policy='reuse', backend=session.backend,
                                  reused=len(steps), executed=0, failed=0, steps=steps)
    return result


def tune(owner, functional, atoms, coordinates, charge, multiplicity, search, *,
         settings=None, checkpoint=None, legacy=None):
    """Return the selection for this exact tuning input, or perform the search."""
    session = calculation_cache._active.get()
    if session is None:
        return search()
    identity = {
        'schema': 1,
        'backend': session.backend,
        'program': calculation_cache._program(session.backend, getattr(owner, 'executable', None)),
        'variant': session.variant,
        'functional': functional.upper(),
        'basis': str(owner.basis).lower(),
        'solvent': owner.solvent,
        'geometry': [[str(atom), *map(float, row)] for atom, row in zip(atoms, coordinates)],
        'charge': int(charge),
        'multiplicity': int(multiplicity),
        'settings': settings,
        'checkpoint': calculation_cache._hash(checkpoint) if checkpoint else None,
    }
    key = hashlib.sha256(json.dumps(identity, sort_keys=True, allow_nan=False).encode()).hexdigest()
    selection = None
    if session.enabled:
        for saved in (session.ktlc, session.previous_ktlc):
            record = saved.get('tuning', {}).get(key, {})
            if record.get('conditions') == identity and valid_selection(session.backend, record.get('selection')):
                selection = record['selection']
                break
    if (selection is None and session.enabled and legacy is not None
            and not session.previous_ktlc and session.previous_result.get('log') == 'normal'):
        selection = legacy(session.previous_result)
    if valid_selection(session.backend, selection):
        action, reason = 'reused', 'verified KTLC parameter selection for matching tuning input'
    else:
        selection = search()
        action, reason = 'executed', 'no verified selection for the KTLC tuning input'
    session.steps.append(dict(job='ktlc', action=action, reason=reason, key=key))
    if valid_selection(session.backend, selection):
        session.save_ktlc('tuning', {key: {'conditions': identity, 'selection': selection}})
    else:
        session.save_ktlc('tuning', {})
    return selection.copy()


def selection_from_parameters(parameters, report):
    # Extensions overriding the old optimizer can still return its parameter list.
    # Only the measured optimization report can make that selection reusable.
    return report if report is not None else {
        'state': 'unverified', 'selected_mu': next(iter(parameters or ()), None),
    }


def legacy_selection(result, owner, name, render, *, verify=None):
    """Import an old selection when its executed input and tuning report agree."""
    session = calculation_cache._active.get()
    selection = result.get('functional_parameter_optimization')
    report = result.get('reuse_report', {})
    if not valid_selection(session.backend, selection) or report.get('backend') != session.backend:
        return None
    expected = render(selection['selected_mu'])
    inputs = list(legacy_inputs(owner, result))
    for job, text in inputs:
        if job == name and calculation_cache.canonical_input(session.backend, text) == expected:
            return selection if verify is None or verify(inputs, selection) else None


def legacy_inputs(owner, result):
    """Yield hash-verified executed inputs associated with the old result."""
    session = calculation_cache._active.get()
    program = calculation_cache._program(session.backend, getattr(owner, 'executable', None))
    for step in result.get('reuse_report', {}).get('steps', ()):
        if step.get('action') not in {'executed', 'reused'}:
            continue
        entry = session.index['entries'].get(step.get('key'))
        if (not entry or entry.get('program') != program or entry.get('variant') != session.variant
                or not entry.get('input')):
            continue
        name = step['job']
        filename = name + '.log'
        root = session.state / entry['snapshot']
        log = root / filename
        if (Path(filename).name == filename and root.resolve().is_relative_to(session.state)
                and log.is_file() and calculation_cache._hash(log) == entry.get('files', {}).get(filename)
                and calculation_cache.successful(session.backend, entry['input'], log.read_text(errors='replace'))):
            yield name, entry['input']


def valid_selection(backend, selection):
    if not isinstance(selection, dict) or selection.get('state') not in {'converged', 'best_evaluated'}:
        return False
    error, mu = selection.get('error_ev'), selection.get('selected_mu')
    if not isinstance(error, (int, float)) or not math.isfinite(error) or error < 0:
        return False
    if not isinstance(selection.get('valid_evaluations'), int) or selection['valid_evaluations'] < 1:
        return False
    if backend == 'gamess' and mu is None:
        # A successful GAMESS evaluation may select the engine default.
        return True
    return (isinstance(mu, (int, float)) and math.isfinite(mu) and 0 <= mu <= 1
            and (backend != 'orca' or mu > 0))
