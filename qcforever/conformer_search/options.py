"""Resolve settings before QCforever changes working directory."""
from .config import SearchConfig


def resolve_options(option_string, override=None):
    tokens = option_string.split()
    if 'laqa' in [t.lower() for t in tokens]:
        raise NotImplementedError('The laqa option is reserved; LAQA is not enabled yet. Omit it to relax all candidates.')
    levels = {'optconf_high': 'high', 'optconf_middle': 'middle', 'optconf_light': 'light'}
    selected = [levels[t.lower()] for t in tokens if t.lower() in levels]
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
    return SearchConfig.resolve(selected[0] if selected else 'light', override)


def calculation_tokens(option_string):
    return [t for t in option_string.split()
            if t.lower() not in ('optconf_high', 'optconf_middle', 'optconf_light')]
