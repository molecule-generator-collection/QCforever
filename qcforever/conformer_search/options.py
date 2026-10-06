"""Resolve settings before QCforever changes working directory."""
from .config import SearchConfig


LEVEL_OPTIONS = {'optconf_low': 'low', 'optconf_medium': 'medium', 'optconf_high': 'high'}


def resolve_options(option_string, override=None):
    tokens = option_string.split()
    if 'laqa' in [t.lower() for t in tokens]:
        raise NotImplementedError('The laqa option is reserved; LAQA is not enabled yet. Omit it to relax all candidates.')
    unknown_levels = [t for t in tokens if t.lower().startswith('optconf_')
                      and t.lower() not in LEVEL_OPTIONS]
    if unknown_levels:
        raise ValueError(f'Unknown conformer level: {", ".join(unknown_levels)}. '
                         'Use optconf_low, optconf_medium, or optconf_high.')
    selected = [LEVEL_OPTIONS[t.lower()] for t in tokens if t.lower() in LEVEL_OPTIONS]
    methods = [t.split('=', 1)[1] if '=' in t else 'pm6' for t in tokens
               if t.lower() == 'optconf' or t.lower().startswith('optconf=')]
    if len(selected) > 1 or len(methods) > 1:
        raise ValueError('Specify one optconf backend and at most one level')
    if not methods:
        if selected or override is not None:
            raise ValueError('conformer settings require optconf')
        return None
    if methods[0] not in ('xtb', 'pm6'):
        raise ValueError('optconf supports xtb or pm6')
    return SearchConfig.resolve(selected[0] if selected else 'low', override)


def calculation_tokens(option_string):
    return [t for t in option_string.split()
            if t.lower() not in LEVEL_OPTIONS]
