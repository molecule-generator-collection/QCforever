"""Search settings: read options/YAML, validate values, and resolve resource limits."""
from dataclasses import asdict, dataclass, field, replace
from copy import deepcopy
from importlib.resources import files
import json
import math
import re
from pathlib import Path
import yaml

from rdkit.Chem import rdMolDescriptors


# Values live in one user-readable file; backend-specific resolution is below.
_DEFAULTS = yaml.safe_load(files(__package__).joinpath('defaults.yaml').read_text())
GRAYBOX_DEFAULTS = {key: _DEFAULTS['relaxation'][key] for key in (
    'algorithm', 'convergence_fraction', 'first_interval', 'subsequent_interval',
    'initial_steps_per_candidate', 'budget_increment_steps_per_candidate',
    'maximum_budget_rounds', 'initial_delta_force', 'delta_force_floor',
    'parallel_candidates', 'timeout_seconds_per_call',
)}


def resolve_relaxation(settings, method):
    """Resolve backend-dependent defaults before generation or any native run."""
    result = dict(settings)
    implementation = result.get('implementation', 'auto')
    if implementation == 'auto':
        implementation = 'graybox'
    result['implementation'] = implementation
    # Disabling scheduling must not change the native optimizer or tolerances.
    if method == 'pm6':
        result['pm6_optimizer'] = 'rfo'
    if implementation == 'graybox':
        if method == 'xtb':
            for key in ('xtb_pool_memory_mb', 'xtb_active_timeout_seconds'):
                result.setdefault(key, _DEFAULTS['relaxation'][key])
        result.update({key: result.get(key, value) for key, value in GRAYBOX_DEFAULTS.items()})
    # Public entry points accept unresolved settings too. Reusing an already
    # resolved object avoids reconstructing/validating the whole configuration.
    return settings if result == settings else result


def validate_relaxation(settings):
    allowed = {'implementation', 'maximum_cycles', 'xtb_executable', 'xtb_opt_level', 'pm6_optimizer',
               'xtb_accuracy', 'xtb_scf_iterations', 'xtb_pool_memory_mb', 'xtb_active_timeout_seconds'} | set(GRAYBOX_DEFAULTS)
    if set(settings) - allowed:
        raise ValueError('Unknown relaxation setting: ' + ', '.join(sorted(set(settings) - allowed)))
    if settings.get('implementation', 'auto') not in ('auto', 'continuous', 'graybox'):
        raise ValueError('relaxation.implementation must be auto, continuous or graybox')
    if settings.get('pm6_optimizer', 'rfo') != 'rfo':
        raise ValueError('pm6_optimizer supports rfo only; use conf_gen=legacy for the old workflow')
    if settings.get('algorithm', 'sh') not in ('laqa', 'sr', 'sh'):
        raise ValueError('relaxation.algorithm must be laqa, sr or sh')
    fraction = settings.get('convergence_fraction', GRAYBOX_DEFAULTS['convergence_fraction'])
    if isinstance(fraction, bool) or not isinstance(fraction, (int, float)) or not math.isfinite(fraction) or not 0 < fraction <= 1:
        raise ValueError('convergence_fraction must be in (0, 1]')
    integers = {'maximum_cycles': 1000, **{k: v for k, v in GRAYBOX_DEFAULTS.items() if type(v) is int}}
    integers.update(xtb_scf_iterations=250, xtb_pool_memory_mb=8192, xtb_active_timeout_seconds=1800)
    for key, default in integers.items():
        value = settings.get(key, default)
        if type(value) is not int or value < 1:
            raise ValueError(f'{key} must be a positive integer')
    workers = settings.get('parallel_candidates', 'auto')
    if workers != 'auto' and (type(workers) is not int or workers < 1):
        raise ValueError('parallel_candidates must be auto or a positive integer')
    for key in ('initial_delta_force', 'delta_force_floor'):
        value = settings.get(key, GRAYBOX_DEFAULTS[key])
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError(f'{key} must be finite and positive')
    accuracy = settings.get('xtb_accuracy', 1.0)
    if isinstance(accuracy, bool) or not isinstance(accuracy, (int, float)) or not math.isfinite(accuracy) or accuracy <= 0:
        raise ValueError('xtb_accuracy must be finite and positive')


@dataclass(frozen=True)
class CandidateBudget:
    """Compute the target count and fixed raw-attempt limit from a molecule."""
    formula: str = _DEFAULTS['budget']['formula']
    base: float = float(_DEFAULTS['budget']['base'])
    growth: float = _DEFAULTS['budget']['growth']
    ring_weight: float = float(_DEFAULTS['budget']['ring_weight'])
    maximum: int = _DEFAULTS['budget']['maximum']
    nonaromatic_ring_floor: int = 0
    fixed: int | None = None
    raw_attempt_multiplier: int = _DEFAULTS['budget']['raw_attempt_multiplier']

    def __post_init__(self):
        if self.formula not in ('exponential_rings', 'fixed'):
            raise ValueError('budget.formula must be exponential_rings or fixed')
        if self.maximum < 1 or self.raw_attempt_multiplier != 2:
            raise ValueError('maximum must be positive; raw_attempt_multiplier must be 2')
        if not all(math.isfinite(x) for x in (self.base, self.growth, self.ring_weight)):
            raise ValueError('Budget coefficients must be finite')
        if self.base < 1 or self.growth < 1 or self.ring_weight < 0:
            raise ValueError('Invalid budget coefficients')
        if not 0 <= self.nonaromatic_ring_floor <= self.maximum:
            raise ValueError('Invalid nonaromatic_ring_floor')
        if self.formula == 'fixed' and (self.fixed is None or not 1 <= self.fixed <= self.maximum):
            raise ValueError('fixed must be between 1 and maximum')

    def resolve(self, mol):
        r = rdMolDescriptors.CalcNumRotatableBonds(mol)
        a = rdMolDescriptors.CalcNumAliphaticRings(mol)
        if self.formula == 'fixed':
            n = self.fixed
        else:
            # Saturate before exponentiation, including very large molecules.
            if self.base + self.ring_weight*a >= self.maximum:
                n = self.maximum
            elif self.growth > 1 and r >= math.log(self.maximum/self.base, self.growth):
                n = self.maximum
            else:
                n = min(self.maximum, math.ceil(self.base*self.growth**r + self.ring_weight*a))
        if a:
            n = max(n, self.nonaromatic_ring_floor)
        return {'rotatable_bonds': r, 'aliphatic_rings': a, 'maximum_candidates': n,
                'raw_attempt_cap_per_generator': self.raw_attempt_multiplier*n}


@dataclass(frozen=True)
class ValidationSettings:
    """Shared geometry, stereo and duplicate thresholds; distances are in angstrom."""
    absolute_collision_angstrom: float = 0.25
    minimum_nonbonded_covalent_ratio: float = 0.75
    minimum_bonded_covalent_ratio: float = _DEFAULTS['validation']['minimum_bonded_covalent_ratio']
    maximum_bonded_covalent_ratio: float = 1.50
    duplicate_rmsd_angstrom: float = _DEFAULTS['validation']['duplicate_rmsd_angstrom']
    duplicate_max_matches: int = _DEFAULTS['validation']['duplicate_max_matches']
    require_tetrahedral_stereo_for_best: bool = _DEFAULTS['validation']['require_tetrahedral_stereo_for_best']
    require_ez_stereo_for_best: bool = _DEFAULTS['validation']['require_ez_stereo_for_best']

    def __post_init__(self):
        values = [self.absolute_collision_angstrom, self.minimum_nonbonded_covalent_ratio,
                  self.minimum_bonded_covalent_ratio,
                  self.maximum_bonded_covalent_ratio, self.duplicate_rmsd_angstrom]
        if any(not math.isfinite(x) or x < 0 for x in values):
            raise ValueError('Validation thresholds must be finite and nonnegative')
        if self.minimum_bonded_covalent_ratio >= self.maximum_bonded_covalent_ratio:
            raise ValueError('Minimum bonded-length ratio must be below maximum')
        if isinstance(self.duplicate_max_matches, bool) or not isinstance(self.duplicate_max_matches, int) or self.duplicate_max_matches < 1:
            raise ValueError('duplicate_max_matches must be a positive integer')


@dataclass(frozen=True)
class SearchConfig:
    """Resolved search settings; requested and effective CPU counts stay distinct."""
    profile: str = _DEFAULTS['profile']
    candidate_retention: str = _DEFAULTS['candidate_retention']
    schema_version: int = _DEFAULTS['schema_version']
    seed: int = _DEFAULTS['seed']
    threads: int = _DEFAULTS['threads']
    workers: int = _DEFAULTS['workers']
    device: str = _DEFAULTS['device']
    mm_method: str = _DEFAULTS['mm_method']
    relaxation: dict = field(default_factory=lambda: {
        key: _DEFAULTS['relaxation'][key] for key in ('implementation', 'maximum_cycles')})
    budget: CandidateBudget = field(default_factory=CandidateBudget)
    validation: ValidationSettings = field(default_factory=ValidationSettings)
    generators: dict = field(default_factory=dict)
    mm: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.schema_version != 1:
            raise ValueError('Unknown schema_version')
        if self.profile not in ('low', 'medium', 'high', 'legacy'):
            raise ValueError('profile must be low, medium, high, or legacy')
        if self.candidate_retention != 'merge':
            raise ValueError('Candidates from previous stages must be retained (merge)')
        if self.device not in ('auto', 'cpu', 'gpu'):
            raise ValueError('device must be auto, cpu or gpu')
        if self.workers < 1 or self.mm_method not in ('mmff94s', 'uff', 'none'):
            raise ValueError('Invalid worker count or MM method')
        validate_relaxation(self.relaxation)
        if self.threads < 1 or not 0 <= self.seed < 2**31:
            raise ValueError('threads must be positive and seed must fit a signed 32-bit integer')
        if 'cores_per_calculation' in self.relaxation:
            raise ValueError('Remove relaxation.cores_per_calculation: xTB/PM6 use 1 core per candidate, with workers derived from nproc')

    @classmethod
    def from_mapping(cls, value):
        value = dict(value)
        value['budget'] = CandidateBudget(**value.get('budget', {}))
        value['validation'] = ValidationSettings(**value.get('validation', {}))
        return cls(**value)

    @classmethod
    def load(cls, path):
        path = Path(path)
        if path.suffix.lower() == '.json':
            value = json.loads(path.read_text())
        else:
            value = yaml.safe_load(path.read_text())
        if not isinstance(value, dict):
            raise ValueError('Search settings must be a mapping')
        return cls.from_mapping(value)

    def to_mapping(self):
        return asdict(self)

    def for_backend(self, method):
        """Resolve native defaults once; reuse this object on subsequent calls."""
        relaxation = resolve_relaxation(self.relaxation, method)
        return self if relaxation is self.relaxation else replace(self, relaxation=relaxation)

    @classmethod
    def resolve(cls, profile='low', override=None):
        """Defaults -> registered model commands -> profile -> explicit YAML."""
        value = deepcopy(_DEFAULTS)
        from .model_workers.installed_models import read_registry
        value['generators'].update(read_registry()['generators'])
        value['profile'] = profile
        if override is not None:
            if isinstance(override, cls):
                return override
            if isinstance(override, (str, Path)):
                path = Path(override).expanduser().resolve()
                supplied = json.loads(path.read_text()) if path.suffix == '.json' else yaml.safe_load(path.read_text())
            else:
                supplied = override
            if not isinstance(supplied, dict):
                raise ValueError('Configuration overrides must be a mapping')
            def merge(target, source):
                for key, item in source.items():
                    if isinstance(item, dict) and isinstance(target.get(key), dict):
                        merge(target[key], item)
                    else:
                        target[key] = item
            merge(value, supplied)
        return cls.from_mapping(value)

    def effective_threads(self, allocated_cores):
        """Fit generation within the runner's resolved nproc, even below 4.

        Keep the requested settings unchanged; the pipeline records the actual
        thread count separately and passes it to every generation backend.
        """
        if isinstance(allocated_cores, bool) or not isinstance(allocated_cores, int) or allocated_cores < 1:
            raise ValueError('allocated_cores must be a positive integer')
        return min(self.threads, allocated_cores)

    def parallelism(self, allocated_cores):
        return min(self.workers, allocated_cores // self.effective_threads(allocated_cores))

    def stages(self):
        if self.profile == 'legacy':
            raise ValueError('legacy uses the complete FAFOOM workflow, not generation stages')
        names = {'high': ['ditmc', 'torsional_diffusion', 'etkdgv3'],
                 'medium': ['torsional_diffusion', 'etkdgv3'], 'low': ['etkdgv3']}
        return [(name, self.mm_method) for name in names[self.profile]]


# QCforever option strings are parsed before the working directory changes.
LEVEL_OPTIONS = {'optconf_low': 'low', 'optconf_medium': 'medium',
                 'optconf_high': 'high', 'optconf_legacy': 'legacy'}


def _expand_conformer_options(option_string):
    """Normalize grouped search options for the existing property parsers.

    Only optconf groups are expanded. Other property arguments, including file
    paths and comma-separated spectra, remain unchanged.
    """
    def expand_group(match):
        fields = {}
        body = match.group(1).strip()
        items = body.split(',') if body else []
        for item in items:
            key, separator, value = item.strip().partition('=')
            key, value = key.strip().lower(), value.strip().lower()
            if not separator or not value or key not in ('conf_gen', 'conf_relax', 'laqa'):
                raise ValueError('Use optconf=(conf_gen=low, conf_relax=pm6, laqa=50)')
            if key in fields:
                raise ValueError(f'Duplicate optconf setting: {key}')
            if any(character.isspace() for character in value):
                raise ValueError(f'Invalid optconf value: {item}')
            fields[key] = value
        if fields.get('conf_gen', 'low') not in ('low', 'medium', 'high', 'legacy'):
            raise ValueError('conf_gen must be low, medium, high, or legacy')
        if fields.get('conf_relax', 'pm6') not in ('xtb', 'pm6'):
            raise ValueError('conf_relax must be xtb or pm6')
        tokens = ['optconf=' + fields.get('conf_relax', 'pm6')]
        if 'conf_gen' in fields:
            tokens.append('optconf_' + fields['conf_gen'])
        if 'laqa' in fields:
            tokens.append('laqa=' + fields['laqa'])
        return ' '.join(tokens)

    expanded = re.sub(r'(?<!\S)optconf\s*=\s*\(([^()]*)\)(?=\s|$)',
                      expand_group, option_string, flags=re.IGNORECASE)
    tokens = expanded.split()
    for index, token in enumerate(tokens):
        lower = token.lower()
        if lower.startswith(('conf_gen', 'conf_relax')):
            raise ValueError('Put conf_gen and conf_relax inside optconf=(...)')
        if lower == 'optconf' and index + 1 < len(tokens) and tokens[index + 1].startswith(('=', '(')):
            raise ValueError('Malformed optconf group; check parentheses and comma-separated settings')
        if lower.startswith('optconf='):
            if '(' in token or ')' in token:
                raise ValueError('Malformed optconf group; check parentheses and comma-separated settings')
            tokens[index] = lower
    return tokens


def parse_conformer_options(option_string, override=None):
    """Read optconf options while leaving other QCforever properties untouched."""
    tokens = _expand_conformer_options(option_string)
    graybox_options = [t.lower() for t in tokens if t.lower() == 'laqa' or t.lower().startswith('laqa=')]
    if len(graybox_options) > 1:
        raise ValueError('Specify laqa at most once')
    unknown_levels = [t for t in tokens if t.lower().startswith('optconf_')
                      and t.lower() not in LEVEL_OPTIONS]
    if unknown_levels:
        raise ValueError(f'Unknown conformer level: {", ".join(unknown_levels)}. '
                         'Use conf_gen=low, medium, high, or legacy inside optconf=(...).')
    selected = [LEVEL_OPTIONS[t.lower()] for t in tokens if t.lower() in LEVEL_OPTIONS]
    methods = [t.split('=', 1)[1].lower() if '=' in t else 'pm6' for t in tokens
               if t.lower() == 'optconf' or t.lower().startswith('optconf=')]
    if len(selected) > 1 or len(methods) > 1:
        raise ValueError('Specify one optconf backend and at most one level')
    if not methods:
        if selected or override is not None or graybox_options:
            raise ValueError('conformer settings require optconf')
        return None
    if methods[0] not in ('xtb', 'pm6'):
        raise ValueError('optconf supports xtb or pm6')
    if selected == ['legacy']:
        if graybox_options or override is not None:
            raise ValueError('conf_gen=legacy uses the old generation and relaxation unchanged; '
                             'laqa and conformer_config overrides are not supported')
        return SearchConfig(profile='legacy')
    config = SearchConfig.resolve(selected[0] if selected else 'low', override)
    value = config.to_mapping()
    if config.profile == 'legacy':
        raise ValueError('Select legacy through optconf=(conf_gen=legacy), without conformer_config')
    if graybox_options:
        option = graybox_options[0]
        if option == 'laqa=off':
            value['relaxation']['implementation'] = 'continuous'
        else:
            try:
                percentage = (100 * GRAYBOX_DEFAULTS['convergence_fraction']
                              if option == 'laqa' else float(option.split('=', 1)[1]))
            except ValueError as exc:
                raise ValueError('Use laqa, laqa=off or laqa=<percentage>') from exc
            if not math.isfinite(percentage) or not 0 < percentage <= 100:
                raise ValueError('laqa percentage must be in (0, 100]')
            value['relaxation'].update(implementation='graybox', convergence_fraction=percentage / 100)
    return SearchConfig.from_mapping(value).for_backend(methods[0])


def calculation_tokens(option_string):
    """Remove only search-level tokens before the existing property parser runs."""
    return [t for t in _expand_conformer_options(option_string)
            if t.lower() not in LEVEL_OPTIONS and t.lower() != 'laqa' and not t.lower().startswith('laqa=')]
