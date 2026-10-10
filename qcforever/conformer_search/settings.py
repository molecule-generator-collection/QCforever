"""Search settings: read options/YAML, validate values, and resolve resource limits."""
from dataclasses import asdict, dataclass, field
from copy import deepcopy
import json
import math
from pathlib import Path

from rdkit.Chem import rdMolDescriptors


GRAYBOX_DEFAULTS = dict(
    algorithm='sh', convergence_fraction=0.20, first_interval=1,
    subsequent_interval=10, initial_steps_per_candidate=20,
    budget_increment_steps_per_candidate=20, maximum_budget_rounds=50,
    initial_delta_force=1.0, delta_force_floor=1e-6,
    parallel_candidates='auto', timeout_seconds_per_call=1800,
)


def resolve_relaxation(settings, method):
    """Resolve backend-dependent defaults before generation or any native run."""
    result = deepcopy(settings)
    implementation = result.get('implementation', 'auto')
    if implementation == 'auto':
        implementation = 'graybox'
    result['implementation'] = implementation
    if implementation == 'graybox':
        if method == 'pm6':
            result['pm6_optimizer'] = 'rfo'
        else:
            result.setdefault('xtb_pool_memory_mb', 8192)
            result.setdefault('xtb_active_timeout_seconds', 1800)
        result.update({key: result.get(key, value) for key, value in GRAYBOX_DEFAULTS.items()})
    return result


def validate_relaxation(settings):
    allowed = {'implementation', 'maximum_cycles', 'xtb_executable', 'xtb_opt_level', 'pm6_optimizer',
               'xtb_accuracy', 'xtb_scf_iterations', 'xtb_pool_memory_mb', 'xtb_active_timeout_seconds'} | set(GRAYBOX_DEFAULTS)
    if set(settings) - allowed:
        raise ValueError('Unknown relaxation setting: ' + ', '.join(sorted(set(settings) - allowed)))
    if settings.get('implementation', 'auto') not in ('auto', 'continuous', 'graybox'):
        raise ValueError('relaxation.implementation must be auto, continuous or graybox')
    if settings.get('pm6_optimizer', 'rfo') != 'rfo':
        raise ValueError('Explicit pm6_optimizer currently supports rfo only; omit for legacy continuous Opt')
    if settings.get('algorithm', 'sh') not in ('laqa', 'sr', 'sh'):
        raise ValueError('relaxation.algorithm must be laqa, sr or sh')
    fraction = settings.get('convergence_fraction', 0.20)
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
    formula: str = 'exponential_rings'
    base: float = 10.0
    growth: float = 1.3
    ring_weight: float = 5.0
    maximum: int = 100
    nonaromatic_ring_floor: int = 0
    fixed: int | None = None
    raw_attempt_multiplier: int = 2

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
    minimum_bonded_covalent_ratio: float = 0.60
    maximum_bonded_covalent_ratio: float = 1.50
    duplicate_rmsd_angstrom: float = 0.10
    duplicate_max_matches: int = 10000
    require_tetrahedral_stereo_for_best: bool = True
    require_ez_stereo_for_best: bool = True

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
    profile: str = 'low'
    candidate_retention: str = 'merge'
    schema_version: int = 1
    seed: int = 20261006
    threads: int = 4
    workers: int = 8
    device: str = 'auto'
    mm_method: str = 'mmff94s'
    relaxation: dict = field(default_factory=lambda: {'implementation': 'auto', 'maximum_cycles': 1000})
    budget: CandidateBudget = field(default_factory=CandidateBudget)
    validation: ValidationSettings = field(default_factory=ValidationSettings)
    generators: dict = field(default_factory=dict)
    mm: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.schema_version != 1:
            raise ValueError('Unknown schema_version')
        if self.profile not in ('low', 'medium', 'high'):
            raise ValueError('profile must be low, medium, or high')
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
            try:
                import yaml
            except ImportError as exc:
                raise ImportError('YAML settings require PyYAML') from exc
            value = yaml.safe_load(path.read_text())
        if not isinstance(value, dict):
            raise ValueError('Search settings must be a mapping')
        return cls.from_mapping(value)

    def to_mapping(self):
        return asdict(self)

    @classmethod
    def resolve(cls, profile='low', override=None):
        """Defaults -> registered model commands -> profile -> explicit YAML."""
        from importlib.resources import files
        import yaml
        value = yaml.safe_load(files(__package__).joinpath('defaults.yaml').read_text())
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
        names = {'high': ['ditmc', 'torsional_diffusion', 'etkdgv3'],
                 'medium': ['torsional_diffusion', 'etkdgv3'], 'low': ['etkdgv3']}
        return [(name, self.mm_method) for name in names[self.profile]]


# QCforever option strings are parsed before the working directory changes.
LEVEL_OPTIONS = {'optconf_low': 'low', 'optconf_medium': 'medium', 'optconf_high': 'high'}


def parse_conformer_options(option_string, override=None):
    """Read optconf options while leaving other QCforever properties untouched."""
    tokens = option_string.split()
    graybox_options = [t.lower() for t in tokens if t.lower() == 'laqa' or t.lower().startswith('laqa=')]
    if len(graybox_options) > 1:
        raise ValueError('Specify laqa at most once')
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
        if selected or override is not None or graybox_options:
            raise ValueError('conformer settings require optconf')
        return None
    if methods[0] not in ('xtb', 'pm6'):
        raise ValueError('optconf supports xtb or pm6')
    config = SearchConfig.resolve(selected[0] if selected else 'low', override)
    value = config.to_mapping()
    if graybox_options:
        option = graybox_options[0]
        if option == 'laqa=off':
            value['relaxation']['implementation'] = 'continuous'
        else:
            try:
                percentage = 20.0 if option == 'laqa' else float(option.split('=', 1)[1])
            except ValueError as exc:
                raise ValueError('Use laqa, laqa=off or laqa=<percentage>') from exc
            if not math.isfinite(percentage) or not 0 < percentage <= 100:
                raise ValueError('laqa percentage must be in (0, 100]')
            value['relaxation'].update(implementation='graybox', convergence_fraction=percentage / 100)
    value['relaxation'] = resolve_relaxation(value['relaxation'], methods[0])
    return SearchConfig.from_mapping(value)


def calculation_tokens(option_string):
    """Remove only search-level tokens before the existing property parser runs."""
    return [t for t in option_string.split()
            if t.lower() not in LEVEL_OPTIONS and t.lower() != 'laqa' and not t.lower().startswith('laqa=')]
