# Configurable conformer search (development branch)

The public Gaussian/GAMESS constructor, property option string, return dictionary
and `pklsave` remain unchanged. `optconf` now selects a **new algorithm**, so
command syntax compatibility does not imply numerical compatibility with legacy
FAFOOM/LAQA. The Gaussian/GAMESS optconf entry points call
`conformer_search.conformer_search.configured_confopt` directly. The legacy
`LAQA_confopt_main` retains its original signature and implementation; it has no
new-search dispatch branch and remains available for explicit legacy use.

```python
from qcforever.gaussian_run import GaussianRunPack
job = GaussianRunPack.GaussianDFTRun(
    'B3LYP', '6-31G*', 8,
    'optconf=xtb optconf_high opt energy uv', 'molecule.sdf', pklsave=True)
job.conformer_config = 'conformer.yaml'  # optional; resolved before changing cwd
result = job.run_gaussian()
```

* `optconf` / `optconf=pm6`: PM6 with Gaussian16. `optconf=xtb`: GFN2-xTB.
* No level flag (or `optconf_low`): ETKDGv3.
* `optconf_medium`: Torsional Diffusion -> ETKDGv3.
* `optconf_high`: DiTMC -> Torsional Diffusion -> ETKDGv3.
* `opt` requests subsequent QC geometry optimization; `uv`, `energy`, etc.
  continue to request the usual downstream properties.
* The conformer backend remains separate from the subsequent Gaussian/GAMESS
  backend. PM6 conformer search requires Gaussian16 even for subsequent GAMESS.

Profile names in YAML and output records are `low`, `medium`, and `high`.
Earlier development names `light`/`middle` and their `optconf_` flags are not
aliases: they raise a configuration error. Historical result files are not renamed.

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

CPU allocation comes from QCforever's resolved `nproc`, not YAML. Generation workers
use `effective_threads = min(threads, nproc)` per worker, with
`min(workers, nproc // effective_threads)` workers. Thus, even if the existing
runner reduces `nproc` below the requested `threads`, generation stays within
the resolved CPU limit instead of failing. `resolved_config.json` retains the
requested settings; `status.json` records the resolved core limit, requested
threads and actual `threads_per_worker`.
Default `device: auto` checks GPU availability **inside each model's Python
environment**, respecting `CUDA_VISIBLE_DEVICES` from the scheduler. A usable
GPU selects one model worker on that GPU; otherwise CPU generation uses the
parallelism below. ETKDG and xTB/PM6 remain CPU calculations. `device: cpu` or
`device: gpu` in YAML forces a choice; a requested but unavailable GPU is logged
as a generation failure, not silently called a GPU success. Multiple GPUs are
not used concurrently by this initial implementation.

Default workers=8 and threads=4 means four CPU cores per model worker, with
at most eight workers: an 8-core allocation runs two workers, a 4-core allocation
runs one. With `nproc=1`, `2` or `3`, one worker uses that many threads.
It does not reserve cores independently of `nproc`. Memory is not automatically
estimated for ML models. Lower workers if model copies exceed available RAM.

xTB/PM6 use **1 core per candidate, with multiple candidates in parallel**.
`nproc=8`, `16`, and `32` allow 8, 16, and 32 relaxation workers respectively,
capped by the available candidate count: `min(nproc, number_of_candidates)`.
Each worker receives one core; `nproc=1` runs candidates sequentially.
Generation `workers`/`threads` are independent of this relaxation policy.
xTB receives `--parallel 1`, and PM6 receives `%nprocshared=1`.
For both native backends, `OMP_NUM_THREADS`, `OMP_THREAD_LIMIT`,
`MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `NUMEXPR_NUM_THREADS` are set to
the per-candidate core count, overriding inherited values only for that calculation.
Parallel workers are separate spawned processes so PM6's working directory and
environment are isolated. On Linux, each native worker is pinned to a disjoint
single CPU from the scheduler-provided affinity. The GENKAI Gaussian wrapper
also caps its CPU list at the native per-call thread limit.
The adapter restores changed variables on success or
failure. Scripts invoking QCforever must use an `if __name__ == '__main__':`
entry-point guard for multiprocessing. These values are recorded in xTB `command.json`
and PM6 `environment.json`. Native thread limits are not a claim of full CPU
utilization at every optimization step.
The former `relaxation.cores_per_calculation` YAML key is no longer accepted;
remove it from old overrides. Audits still report the realized
`cores_per_calculation` and actual `parallel_workers` for reproducibility.
The supplied memory setting is **per native calculation**, not a shared pool:
With `nproc=32` and at least 32 candidates, `mem='1GB'` can request up to
32 GB in Gaussian, plus worker overhead. Choose `nproc` and memory together.

An empty raw SDF is treated as zero generated candidates and follows the same
fixed-cap and fallback policy as any other empty pool.

## Generation contract

`N = min(100, ceil(10 * 1.3**r + 5*a))`, using RDKit's default rotatable-bond
count and aliphatic-ring count on the hydrogen-suppressed input graph.
Each stage first requests the remaining quota, `N - accepted_pool_size` (N at
the start of the search). For example, with N=20 and 18 candidates retained from
earlier stages, the next stage first requests only 2. Further batches request at most
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
Final structures do not inherit generation-time stereo pass flags. When geometry
is invalid and stereo is not assessed, `stereo_check_status` is
`not_evaluated_geometry_invalid` and the stereo-match boolean properties are
absent. This applies to the selected SDF, all-converged SDF and per-candidate SDFs.

### Optional learned-model environments

The core package does not install Torch/JAX or download models. Run
`install-conformer-models` once to build separate environments and register
DiTMC/TD after real generation tests. The registered commands are loaded
automatically; no per-job YAML is needed. An explicit YAML overrides registered
settings. For advanced/custom workers, the command interface is:

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
With `device: auto`, custom DiTMC/TD commands must also implement the supplied
worker's `--probe-device` protocol (write availability JSON to `{output}`), or
explicitly select `device: cpu`/`gpu`. Probe failures are logged, not hidden.
Stage status records the resolved device and worker count; model execution logs
record the actual parameter devices. Device discovery imports the framework
once in a short-lived process but does not load weights. Model weights remain
loaded in the persistent worker across subsequent generation batches.
Each request JSON includes SMILES, exact raw count (`maximum_raw_candidates`),
seed, threads and adapter options. The worker writes explicit-H raw SDF without
MM, energy filtering or validator-based pruning. Missing/unconfigured models
are logged and fall back to the next stage. Worker errors are retained; all
workers failing ends that stage. Partial worker output is retained.
Model workers must enforce CPU/device and framework-specific thread limits;
OMP/MKL/OpenBLAS/NumExpr limits are also passed by the orchestrator.
Subprocesses execute independently with deterministic ordered merging.
TD explicitly passes each request's seed to RDKit embedding as well as seeding
Python/NumPy/Torch. Persistent additional batches replace the embedding seed;
they do not reuse the first request's seed. This does not guarantee bitwise
identity across different library versions/devices or nondeterministic kernels.

The setup command and supplied GENKAI configuration use `persistent: true`. This requires
a worker supporting `--session DIRECTORY` (the installed `qcforever-model-worker`
supports it; `scripts/model_raw_worker.py` is a compatibility shim).
Each generator stage maintains up to `min(workers, floor(nproc / effective_threads))` worker
processes. Workers load model weights once, accept additional batches, and close
when the stage ends. A smaller additional batch reuses existing workers rather
than rebuilding a smaller pool. Model stdout is saved in `session_worker_*/worker.out`.
Each request saves `model_execution.json` with PID, request number, initialization
and generation times; generation time may include first-use JIT/shape compilation.
The model is not silently restarted on failure. Unrelated custom one-shot commands
remain supported with `persistent: false` (the default for command adapters).
The model workers are now shipped in `qcforever_model_workers`, separately from
the core imports. Source/checkpoint locations are registered by setup or explicitly configured, and no
external benchmark loader is imported. CPU requirements and setup instructions
are in [model installation](model_installation.md); fresh-install validation is
separate from reuse of existing GENKAI environments.

## MM and relaxation

MM is applied **only to ETKDGv3-generated candidates**, never to DiTMC or Torsional
Diffusion candidates. MMFF94s is default; UFF and none are also selectable.
Missing parameters or runtime failure skips MM for only the affected ETKDGv3
candidate, retaining its pre-MM coordinates and the other candidates' MM results.
`mm_candidate_runs` records each outcome; mixed outcomes use `partially_skipped`.
There is no silent force-field switch. Iteration-limit status is
recorded and is not declared convergence. Post-MM validation may reduce the pool.
Force-field APIs operate on disposable copies: even availability/parameter
checks can change RDKit aromaticity flags. Only the optimized coordinates are
copied back onto the preserved input graph, retaining formal charges, radicals,
isotopes and bond information. xTB and PM6 use the same coordinate-only handoff.
This prevents serialization artifacts, not physical reactions: new coordinates
still undergo distance and stereo checks, and atom-local radical labels are not
claimed to be a calculated spin-density distribution.
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
MM status/skip reason, stage times and a relative detail directory. Gaussian's
selected-SDF readback also records `structure_handoff`. A readback failure changes
the overall state to `failed`, keeps the prior `search_state`, and records
`failure_stage: selected_structure_readback` in memory and saved summary. The usual
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

During generation, an incremental filter validates only each new batch against
the unchanged accepted pool. Earlier accepted candidates are retained in order;
new candidates are compared against both earlier batches and accepted candidates
from the same batch. Geometry, input-specified stereo and symmetry-aware RMSD
rules are unchanged. The filter owns private copies of its pool and fixes its
reference/settings for one preparation run. Final preparation checks follow
two explicit routes:

- No MM call: export a copy of the privately validated generation pool and audit.
  This covers learned-only pools and `mm_method: none`; the latter never calls
  the optimizer (including a supplied custom optimizer).
- Any MM call: fully recheck the entire output pool, including cross-generator
  duplicates in mixed pools. This also applies if MM leaves coordinates unchanged
  or every MM call is skipped because parameters are unavailable.

There is no structural-change detector or binary-identity comparison. Validation
thresholds, candidate order/IDs, properties and the final audit schema are unchanged.
`post_mm_validation_mode` records `reused_generation_no_mm` or `full_recheck`,
and `post_mm_validation_wall_seconds` records the final export/check time.
The post-xTB/PM6 geometry/stereo audit remains unconditional and separate.
Batch audit indices/counts preserve
the previous accepted-prefix-plus-new-batch convention; they are not cumulative
raw-attempt counts.

## Basic tests

Run the compact conformer test suite from the repository root:

```bash
python -m pytest -q test/test_conformer_search.py test/test_conformer_setup.py test/test_model_session.py test/test_native_parallel.py
```

These tests cover search settings, bounded generation/fallback, structure and
stereo checks, MM routing, selection, CPU limits, worker reuse and installer
safety. They use actual RDKit ETKDG/MM and synthetic model/native commands;
model weights, GPUs, Gaussian and xTB executables are not required. Extended
benchmark and cluster-validation scripts are not distributed with the package.

### Previous integration checks

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
