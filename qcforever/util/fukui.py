"""Condensed Fukui indices from Mulliken partial charges."""

from math import isfinite


def condensed_fukui(reference, electron_removed, electron_added):
    """Compare the N, N-1 and N+1 electron states in the same atom order."""
    atoms, charges, _ = reference
    for state, population in (('N', reference), ('N-1', electron_removed), ('N+1', electron_added)):
        elements, partial_charges, _ = population
        if (not atoms or elements != atoms or len(partial_charges) != len(atoms)
                or not all(isfinite(charge) for charge in partial_charges)):
            raise ValueError(f'Invalid Mulliken charges for Fukui calculation in the {state} state')
    return {
        'atoms': list(atoms),
        'f_plus': [q - q_added for q, q_added in zip(charges, electron_added[1])],
        'f_minus': [q_removed - q for q, q_removed in zip(charges, electron_removed[1])],
    }
