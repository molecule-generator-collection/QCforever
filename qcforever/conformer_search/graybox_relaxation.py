"""Schedule observable optimization blocks with LAQA, SR or SH.

The caller supplies advance(candidate_id) -> BlockResult. It may execute a
native calculation or consume a saved block; scheduling never sees future
energies, future costs or reference minima. Costs include initialization and
failed calls. A dispatched batch may overshoot a round's evaluation budget.
Selection uses only completed batches, never the order in which workers finish.
"""
from dataclasses import dataclass, field
import math


@dataclass
class BlockResult:
    energy: float | None = None                 # Hartree, total molecule
    force: float | None = None                  # Hartree/Bohr, see force_kind
    evaluations: int = 0                       # Includes initial geometry
    scf_cycles: int = 0
    wall_seconds: float = 0.0                   # Whole call including restart
    status: str = 'paused'                      # paused/converged/failed/limit
    error: str | None = None
    force_kind: str = 'mean_atom_force'


@dataclass
class CandidateState:
    atoms: int
    energy: float | None = None
    forces: list = field(default_factory=list)
    evaluations: int = 0
    calls: int = 0
    wall_seconds: float = 0.0
    status: str = 'pending'
    error: str | None = None

    @property
    def terminal(self):
        return self.status in ('converged', 'failed', 'limit')


class GrayboxSearch:
    """Choose candidates, execute a batch, then record results and re-rank.

    With one worker this retains the original sequential policy. With several
    workers, LAQA advances the best distinct candidates; SR/SH advance candidates
    toward the current stage target before eliminating any. Already dispatched
    blocks are collected and charged even if the convergence quota is reached.
    """

    def __init__(self, candidate_atoms, advance, settings, on_event=None, force_kind='mean_atom_force',
                 *, workers=1, advance_batch=None):
        if not candidate_atoms or any(n < 1 for n in candidate_atoms.values()):
            raise ValueError('Graybox search requires nonempty candidates with atom counts')
        self.states = {key: CandidateState(n) for key, n in sorted(candidate_atoms.items())}
        self.advance_native = advance
        if type(workers) is not int or workers < 1:
            raise ValueError('workers must be a positive integer')
        if workers > 1 and advance_batch is None:
            raise ValueError('Parallel search requires a batch executor')
        self.workers = workers
        self.advance_native_batch = advance_batch
        self.batch_number = 0
        self.settings = settings
        self.on_event = on_event
        if force_kind not in ('mean_atom_force', 'anc_norm_per_sqrt_atom'):
            raise ValueError('Unknown force definition')
        self.force_kind = force_kind
        self.initial_count = len(self.states)
        self.required = max(1, math.ceil(settings['convergence_fraction'] * self.initial_count))
        self.cost = dict(evaluations=0, scf_cycles=0, wall_seconds=0.0, calls=0)
        self.events = []
        self.round = 0
        self.limit = self.initial_count * settings['initial_steps_per_candidate']
        self.phase = 'initialization'
        self.stop_reason = None

    def run(self):
        """Initialize once, allocate more evaluations, then report why we stopped."""
        candidates = list(self.states)
        for start in range(0, len(candidates), self.workers):
            if self.done():
                break
            self.advance_candidates(candidates[start:start + self.workers])
        for round_number in range(1, self.settings['maximum_budget_rounds'] + 1):
            self.round = round_number
            self.limit = self.initial_count * (
                self.settings['initial_steps_per_candidate'] +
                (round_number - 1) * self.settings['budget_increment_steps_per_candidate'])
            if self.done() or not self.unfinished():
                break
            if self.settings['algorithm'] == 'laqa':
                self.phase = 'LAQA'
                self.advance_ranked_candidates('laqa')
            else:
                self.run_elimination_round()
        self.stop_reason = ('convergence_quota' if self.done() else
                            'candidate_pool_exhausted' if not self.unfinished() else
                            'global_budget_cap')
        return self

    def unfinished(self):
        return [key for key, state in self.states.items() if not state.terminal]

    def done(self):
        return sum(state.status == 'converged' for state in self.states.values()) >= self.required

    def can_advance(self):
        return not self.done() and self.cost['evaluations'] < self.limit

    def advance(self, candidate):
        """Advance one candidate (also used by sequential saved-data callers)."""
        self.advance_candidates([candidate])

    def advance_candidates(self, candidates):
        """Collect a whole batch in dispatch order, including quota overshoot."""
        if not candidates or len(candidates) > self.workers or len(set(candidates)) != len(candidates):
            raise ValueError('Batch must contain distinct candidates within the worker limit')
        if self.done() or any(self.states[key].terminal for key in candidates):
            raise RuntimeError('Cannot advance a finished search or candidate')
        if self.advance_native_batch is None:
            blocks = {key: self.advance_native(key) for key in candidates}
        else:
            blocks = self.advance_native_batch(candidates)
        if set(blocks) != set(candidates):
            raise ValueError('Batch results must match dispatched candidates')
        self.batch_number += 1
        for candidate in candidates:
            self.record_block_result(candidate, blocks[candidate], len(candidates))

    def record_block_result(self, candidate, block, batch_size=1):
        """Validate one endpoint, update measured costs, and persist an event."""
        state = self.states[candidate]
        if block.status not in ('paused', 'converged', 'failed', 'limit'):
            raise ValueError('Unknown block status')
        if block.evaluations < 0 or block.scf_cycles < 0 or not math.isfinite(block.wall_seconds) or block.wall_seconds < 0:
            raise ValueError('Invalid measured block cost')
        if block.status in ('paused', 'converged'):
            if block.energy is None or not math.isfinite(block.energy):
                block.status, block.error = 'failed', 'Missing or nonfinite endpoint energy'
            elif self.settings['algorithm'] == 'laqa' and (
                    block.force_kind != self.force_kind or block.force is None or
                    not math.isfinite(block.force) or block.force < 0):
                block.status, block.error = 'failed', 'LAQA requires a finite force with the configured definition'
            elif block.status == 'paused' and block.evaluations == 0:
                block.status, block.error = 'failed', 'No progress in nonterminal block'
        if block.energy is not None and math.isfinite(block.energy):
            state.energy = block.energy
        state.forces.append(block.force)
        state.status, state.error = block.status, block.error
        state.calls += 1
        state.evaluations += block.evaluations
        state.wall_seconds += block.wall_seconds
        for name in ('evaluations', 'scf_cycles', 'wall_seconds'):
            self.cost[name] += getattr(block, name)
        self.cost['calls'] += 1
        converged = [key for key, row in self.states.items() if row.status == 'converged']
        best = min(converged, key=lambda key: (self.states[key].energy, key)) if converged else None
        event = dict(candidate=candidate, block=state.calls, cost=dict(self.cost),
                     energy_hartree=state.energy, force=block.force, force_kind=block.force_kind,
                     state=state.status, error=state.error, converged_count=len(converged),
                     best_candidate=best, best_energy=self.states[best].energy if best else None,
                     budget_round=self.round, budget_limit=self.limit, phase=self.phase)
        if self.workers > 1:
            event.update(batch=self.batch_number, batch_size=batch_size)
        self.events.append(event)
        if self.on_event is not None:
            self.on_event(event)

    def rank(self, candidate, score):
        """Energy or LAQA score, with candidate ID as deterministic tie-breaker."""
        state = self.states[candidate]
        if state.status in ('failed', 'limit') or state.energy is None:
            return math.inf, candidate
        if score == 'energy':
            return state.energy, candidate
        force = state.forces[-1]
        difference = self.settings['initial_delta_force'] if len(state.forces) == 1 else max(
            abs(force - state.forces[-2]), self.settings['delta_force_floor'])
        return state.energy / state.atoms - force**2 / (2 * difference), candidate

    def run_elimination_round(self):
        """SR/SH with re-entry on budget extension and residual-budget reuse.

        This deliberately matches the saved-data experiment, not fixed-budget
        SR/SH: eliminations last only for this round. Converged members can stay
        in its ranking but consume no more calls. All unfinished candidates are
        eligible again next round, retaining their optimization state.
        """
        alive = self.unfinished()
        count = len(alive)
        starts = {key: self.states[key].evaluations for key in alive}
        budget = self.limit - self.cost['evaluations']
        harmonic = 0.5 + sum(1 / i for i in range(2, count + 1))
        rounds = max(1, math.ceil(math.log2(count)))
        phase, target = 0, 0
        method = self.settings['algorithm'].upper()
        while len(alive) > 1 and self.can_advance():
            if method == 'SR':
                target = max(1, math.ceil((max(count + 1, budget) - count) / (harmonic * (count - phase))))
            else:
                target += max(1, math.floor(budget / (rounds * len(alive))))
            self.phase = f'{method}_selection_{phase}'
            self.advance_to_stage_target(alive, starts, target)
            if not self.can_advance():
                return
            alive.sort(key=lambda key: self.rank(key, 'energy'))
            keep = len(alive) - 1 if method == 'SR' else math.ceil(len(alive) / 2)
            alive = alive[:keep]
            phase += 1
        self.phase = f'{method}_survivor'
        while alive and not self.states[alive[0]].terminal and self.can_advance():
            self.advance(alive[0])
        self.phase = f'{method}_residual'
        self.advance_ranked_candidates('energy')

    def advance_to_stage_target(self, candidates, starts, target):
        """Finish this SR/SH stage before ranking; never run one candidate twice at once."""
        while self.can_advance():
            waiting = [key for key in candidates if not self.states[key].terminal
                       and self.states[key].evaluations - starts[key] < target]
            if not waiting:
                return
            self.advance_candidates(waiting[:self.workers])

    def advance_ranked_candidates(self, score):
        """Shared allocation for LAQA and the SR/SH residual-budget phase."""
        while self.can_advance() and self.unfinished():
            ranked = sorted(self.unfinished(), key=lambda key: self.rank(key, score))
            self.advance_candidates(ranked[:self.workers])
