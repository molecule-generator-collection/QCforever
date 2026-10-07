"""Validated, serializable settings for the Standard candidate diagram."""
from dataclasses import asdict, dataclass, field
import json
import math
from pathlib import Path

from rdkit.Chem import rdMolDescriptors


@dataclass(frozen=True)
class CandidateBudget:
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
    profile: str = 'low'
    candidate_retention: str = 'merge'
    schema_version: int = 1
    seed: int = 20261006
    threads: int = 4
    workers: int = 8
    device: str = 'auto'
    mm_method: str = 'mmff94s'
    relaxation: dict = field(default_factory=lambda: {'implementation': 'continuous', 'maximum_cycles': 1000})
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
        if self.relaxation.get('implementation') != 'continuous':
            raise ValueError('Only continuous relaxation is enabled; LAQA remains under study')
        if self.threads < 1 or not 0 <= self.seed < 2**31:
            raise ValueError('threads must be positive and seed must fit a signed 32-bit integer')
        if 'cores_per_calculation' in self.relaxation:
            raise ValueError('Remove relaxation.cores_per_calculation: xTB/PM6 use 1 core per candidate, with workers derived from nproc')
        allowed = {'implementation', 'maximum_cycles', 'xtb_executable', 'xtb_opt_level'}
        if set(self.relaxation)-allowed:
            raise ValueError('Unknown relaxation setting')
        if not isinstance(self.relaxation.get('maximum_cycles', 1000), int) or self.relaxation.get('maximum_cycles', 1000) < 1:
            raise ValueError('maximum_cycles must be a positive integer')

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
        from qcforever_model_workers.registry import read_registry
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

    def parallelism(self, allocated_cores):
        if allocated_cores < self.threads:
            raise ValueError('threads per worker exceeds allocated CPU cores')
        return min(self.workers, allocated_cores // self.threads)

    def stages(self):
        names = {'high': ['ditmc', 'torsional_diffusion', 'etkdgv3'],
                 'medium': ['torsional_diffusion', 'etkdgv3'], 'low': ['etkdgv3']}
        return [(name, self.mm_method) for name in names[self.profile]]
