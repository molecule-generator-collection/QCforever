# Configurable conformer search (development branch)

The public Gaussian/GAMESS constructor, property option string, return dictionary
and `pklsave` remain unchanged. `optconf` now selects a **new algorithm**, so
command syntax compatibility does not imply numerical compatibility with legacy
FAFOOM/LAQA. Direct legacy `LAQA_confopt_main` calls without `search_config` retain
the old implementation.

```python
from qcforever.gaussian_run import GaussianRunPack
job = GaussianRunPack.GaussianDFTRun(
    'B3LYP', '6-31G*', 8,
    'optconf=xtb optconf_high opt energy uv', 'molecule.sdf', pklsave=True)
job.conformer_config = 'conformer.yaml'  # optional; resolved before changing cwd
result = job.run_gaussian()
```

* `optconf` / `optconf=pm6`: PM6 with Gaussian16. `optconf=xtb`: GFN2-xTB.
* No level flag (or `optconf_light`): ETKDGv3.
* `optconf_middle`: Torsional Diffusion -> ETKDGv3.
* `optconf_high`: DiTMC -> Torsional Diffusion -> ETKDGv3.
* `opt` requests subsequent QC geometry optimization; `uv`, `energy`, etc.
  continue to request the usual downstream properties.
* The conformer backend remains separate from the subsequent Gaussian/GAMESS
  backend. PM6 conformer search requires Gaussian16 even for subsequent GAMESS.

## Configuration

`qcforever/conformer_search/defaults.yaml` is installed as package data. Resolution
order is packaged defaults -> level option -> explicit YAML (recursive merge).
There is no automatic YAML discovery. Relative paths are resolved against the
caller's working directory before QCforever enters its job directory.
Do not edit installed defaults; copy them or supply just the changed keys:

```yaml
workers: 4
threads: 4
mm_method: mmff94s  # mmff94s / uff / none
```

CPU allocation comes from QCforever's resolved `nproc`, not YAML. Actual workers
are `min(workers, nproc // threads)`; fewer cores than threads is an error.
Default workers=8 and threads=4 means four CPU cores per model worker, with
at most eight workers: an 8-core allocation runs two workers, a 4-core allocation
runs one. It does not reserve cores independently of `nproc`. Memory is not automatically
estimated for ML models. Lower workers if model copies exceed available RAM.

For the GENKAI integration smoke test, `check_sample_profiles.py` accepts
`--smoke-candidates 2`. This is a **test-only in-memory override**, not an edit
to the packaged YAML or normal candidate formula. Omitting that flag restores
the normal budget automatically; results record both budgets. Each sample/profile
is submitted separately with four allocated cores and four QC cores. An empty
raw SDF is treated as zero generated candidates and follows the same fixed-cap
and fallback policy as any other empty pool.

## Generation contract

`N = min(100, ceil(10 * 1.3**r + 5*a))`, using RDKit's default rotatable-bond
count and aliphatic-ring count on the hydrogen-suppressed input graph.
Each stage first requests N raw candidates. Further batches request at most
one candidate per available worker, limited by the remaining target/attempts.
Each stage has a **2N requested-candidate cap** including requests returning no
coordinates. Models run in stage order, never DiTMC and TD concurrently.
Previously accepted candidates are kept and cross-stage duplicates removed.
At least one surviving candidate permits relaxation even if N is not reached.
No rescue retries, alternative seeds after the cap, or bond-order repairs occur.

Checks are common to all methods: readability, finite coordinates, atom
composition, supplied graph, clashes, stretched bonds, symmetry-aware RMSD.
Duplicate RMSD includes all atoms **except hydrogens in XH3 groups**: when a
heavy atom X has exactly three explicit H neighbors, exclude those three H
atoms, regardless of X's element or charge (including CH3 and NH3+). Keep X
itself and all H atoms in XH, XH2 and XH4 groups. OH rotamers therefore remain
distinguishable; rotation of only the XH3 hydrogens is intentionally ignored.
`xh3_hydrogen_indices()` defines this rule and `rmsd_comparison_molecule()`
creates a separate comparison copy before symmetry enumeration. The original
full-H structures are retained for all geometry checks, MM, QC and stereo
audits. Audit logs record the rule and excluded reference atom indices (0-based).
The fit is rigid, without reflection, and does not modify coordinates. Defaults are 0.10 Å
and at most 10,000 symmetry mappings per pair (`duplicate_rmsd_angstrom` and
`duplicate_max_matches` under `validation`). The bounded symmetry search can
retain extra duplicates when the best mapping is outside that cap; it is not
an exhaustive-symmetry guarantee. The same filter is used for every generator,
cross-generator merging, and post-MM filtering. Missing explicit hydrogens are
not silently reconstructed by this filter: composition validation still applies.
Input-specified tetrahedral R/S and double-bond E/Z stereochemistry are checked
from coordinates and are required during generation (and after MM). Unspecified
stereochemistry is unconstrained. This is not a universal isomer classifier:
tautomerism, atropisomerism and coordination stereochemistry are not covered.
The default final primary best requires native convergence, valid geometry and
both specified tetrahedral and E/Z stereochemistry. If no such result remains,
the lowest-energy geometry-valid converged structure is returned with a stereo
warning; if none has valid geometry, the lowest-energy converged structure is
returned with a geometry warning. All converged structures and their audits are
retained. No converged structures is a failure. Geometry checks use distances
and the supplied graph, not a definitive chemical reaction detector.

### Optional learned-model environments

The core package does not install Torch/JAX or download models. DiTMC/TD are
connected by list-valued commands in separate environments:

```yaml
generators:
  ditmc:
    command: [/absolute/model-env/bin/python, /absolute/ditmc_worker.py,
              --request, '{request}', --input, '{input}', --output, '{output}']
    timeout_seconds: 600
  torsional_diffusion:
    command: [/absolute/td-env/bin/python, /absolute/td_worker.py,
              --request, '{request}', --input, '{input}', --output, '{output}']
```

These show the general command interface. For supplied model workers and CPU
installation commands, see [model installation](model_installation.md).
Each request JSON includes SMILES, exact raw count (`maximum_raw_candidates`),
seed, threads and adapter options. The worker writes explicit-H raw SDF without
MM, energy filtering or validator-based pruning. Missing/unconfigured models
are logged and fall back to the next stage. Worker errors are retained; all
workers failing ends that stage. Partial worker output is retained.
Model workers must enforce CPU/device and framework-specific thread limits;
OMP/MKL/OpenBLAS/NumExpr limits are also passed by the orchestrator.
Subprocesses execute independently with deterministic ordered merging.

The supplied GENKAI DiTMC/TD configuration uses `persistent: true`. This requires
a worker supporting `--session DIRECTORY` (the installed `qcforever-model-worker`
supports it; `scripts/model_raw_worker.py` is a compatibility shim).
Each generator stage maintains up to `floor(allocated_cores / threads)` worker
processes. Workers load model weights once, accept additional batches, and close
when the stage ends. A smaller additional batch reuses existing workers rather
than rebuilding a smaller pool. Model stdout is saved in `session_worker_*/worker.out`.
Each request saves `model_execution.json` with PID, request number, initialization
and generation times; generation time may include first-use JIT/shape compilation.
The model is not silently restarted on failure. Unrelated custom one-shot commands
remain supported with `persistent: false` (the default for command adapters).
The model workers are now shipped in `qcforever_model_workers`, separately from
the core imports. Source/checkpoint locations are caller-configured, and no
external benchmark loader is imported. CPU requirements and setup instructions
are in [model installation](model_installation.md); fresh-install validation is
separate from reuse of existing GENKAI environments.

## MM and relaxation

MM is applied **only to ETKDGv3-generated candidates**, never to DiTMC or Torsional
Diffusion candidates. MMFF94s is default; UFF and none are also selectable.
Missing parameters or runtime failure skips MM for the ETKDGv3
ensemble without silently switching force fields. Iteration-limit status is
recorded and is not declared convergence. Post-MM validation may reduce the pool.
Continuous native xTB/PM6 runs separately for each candidate. LAQA scheduling is
not enabled in this new route while its continuation/stopping study is ongoing.
Every surviving candidate is attempted, with no early stop when a good energy
is found. A failed candidate (including native-output parsing failure) is logged
and does not suppress later candidates. The overall job timeout is an exception:
it terminates the workflow and records the interrupted candidate as timeout.
Partial native trajectories are retained for failed/interrupted calculations.
PM6 native/parser failures also save `failure.json`, with a categorized reason,
the original adapter exception and the last 20 native-log lines. An unrecognized
failure is not guessed to be SCF failure. The legacy parser is unchanged.

Future option syntax is `optconf=xtb laqa ...` (or PM6). Currently `laqa` raises
an explicit `NotImplementedError` **before generation**, rather than silently
running full relaxation or the unrelated legacy LAQA implementation. Omitting
`laqa` uses all-candidate continuous relaxation. This is the baseline for later
same-ensemble cost/energy comparisons with LAQA.

## Output and retention

The ordinary QC result retains `optconf: bool` and adds `conformer_search` with
profile, backend, chosen candidate, Hartree energy, realized generation stages,
MM status/skip reason, stage times and a relative detail directory. The usual
DFT `Energy` remains untouched. `pklsave=True` saves this extended dictionary.
Conformer failure is distinct from success of subsequent QC using the input.

Inside the existing molecule job directory:

```
optimized_structures.sdf          # selected first record: legacy downstream handoff
conformer_search/
  resolved_config.json
  reference.sdf
  status.json
  00_<generator>/batch_000/       # raw SDF, accepted pool, seeds, audit, worker logs
  generated_candidates.sdf
  mm_candidates.sdf
  initial_structures.sdf
  electronic/candidate_00000/    # native logs, final SDF, trace/status JSON
  electronic/all_converged.sdf
  electronic/audit.json
  summary.json
```

Candidate IDs persist across stages, MM and relaxation. xTB trajectory records
are **not assumed** to equal optimizer cycles; native cycle headings are counted
independently. PM6 records SCF energy evaluations, not an invented cycle mapping.
Search files are protected from normal QCforever cleanup, including failures.
Existing search directories are not overwritten or silently resumed.

## Validation status

CPU unit tests use actual RDKit ETKDG/MM and isolated synthetic model/QC commands.
Real DiTMC/TD checkpoints and native xTB/PM6 have been exercised on GENKAI with
four Sample molecules (formaldehyde, chlorobenzene, ethanol and `200_11`), using
a test-only maximum of two candidates. These tests demonstrate execution, not
full-budget performance or general charged/open-shell model support. Persistent
workers have also been checked across additional batches. A five-molecule,
normal-budget comparison including radical/anion inputs is in progress; its
completion is not implied here.
The existing optconf entry -> native ETKDG -> all-candidate xTB subprocess ->
best-structure handoff is integration-tested with a synthetic executable. PM6
tests use the existing Gaussian input writer and substitute only native execution.
These tests validate orchestration, **not** physical optimizer behavior.

The compressed-formaldehyde regression is covered by the common
`minimum_bonded_covalent_ratio: 0.60` sanity bound: bonded distance must be at
least 60% of the sum of the two RDKit covalent radii. This configurable threshold
is deliberately permissive, not a physical prediction or a fitted benchmark
criterion. The default maximum ratio remains 1.50. All methods and the final
audit use the same bounds. Normal geometry and compressed geometry are tested.
The pre-fix full-budget comparison retains its original filter; do not relabel
those frozen calculations as post-fix validation.
