"""Explain native PM6 failures without losing the original Gaussian logs.

This diagnostic layer does not change legacy LAQA parsing or convergence rules.
Classification uses explicit log messages; unrecognized failures stay unknown.
"""
import re


class PM6CalculationError(RuntimeError):
    pass


def failure_diagnostic(log, original_error=None):
    rules = (
        ('invalid_interatomic_distances', ('small interatomic distances', 'problem with the distance matrix')),
        ('scf_not_converged', ('convergence failure', 'scf has not converged')),
        ('optimization_step_limit', ('number of steps exceeded',)),
    )
    lower = log.lower()
    reason = next((name for name, messages in rules if any(m in lower for m in messages)), None)
    has_energy = bool(re.search(r'SCF Done:.*?=\s*[-+0-9.]', log))
    if reason is None:
        reason = ('missing_native_log' if not log else 'gaussian_error_termination'
                  if 'error termination' in lower else 'missing_scf_energy'
                  if not has_energy else 'output_parse_failure'
                  if original_error is not None else 'optimization_not_converged')
    return {'reason': reason, 'native_log': 'Gau_molecule.log',
            'log_tail': log.splitlines()[-20:],
            'adapter_error': (f'{type(original_error).__name__}: {original_error}'
                              if original_error is not None else None)}
