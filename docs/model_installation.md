# Optional model workers (Linux)

## One-command setup

After installing QCforever, run:

```bash
install-conformer-models
```

Prerequisites: Linux x86_64, Python 3.11 with `venv`/pip and development headers,
a C/C++ compiler for DiTMC, internet access, and several GB of free disk space.
No `sudo`, system package changes, or modification of existing Python environments
is performed. Setup downloads official third-party source and trusted model
checkpoints (including TD's pickle-based checkpoint); use implies installing
these third-party components under their own license/citation terms.

The default directory is `$XDG_DATA_HOME/qcforever/conformers`, or
`~/.local/share/qcforever/conformers` when unset. It contains isolated environments,
verified downloads, selected source/checkpoint files, caches and timestamped logs.
The DiTMC download is 1.89 GB; its unused validation data/other checkpoints are
not extracted. Official download URLs, version IDs, sizes and checksums are
bundled in `qcforever_model_workers/model_sources.json`. DiTMC uses Zenodo's
published MD5 checksum; TD uses pinned SHA256 hashes. Mismatches stop setup
before extraction or checkpoint loading. These checks do not replace trust in
the original publisher. No trained weights are redistributed in QCforever.

Setup runs the two models **sequentially**, generating ethanol twice (one candidate
each) to check readable finite 3D coordinates, connectivity, absence of collapsed
atoms and reuse of the same model process. It does not run xTB/PM6, prove
stereochemical accuracy on other molecules, or benchmark model quality.
Only passed models are registered in `$XDG_CONFIG_HOME/qcforever/conformer_models.json`
(default `~/.config/qcforever/conformer_models.json`). `QCFOREVER_MODEL_REGISTRY`
can select another registration file. Ordinary calculations read this file;
they do not download/install anything. Explicit `job.conformer_config` settings
take precedence. The registered commands use `persistent: true` automatically.

Useful options:

```bash
install-conformer-models --dry-run
install-conformer-models --device cpu
install-conformer-models --device gpu
install-conformer-models --models torsional_diffusion
install-conformer-models --directory /path/to/my/model-install --python /path/to/python3.11
```

`--device auto` (default) chooses a visible NVIDIA GPU, otherwise CPU. GPU setup
uses CUDA-12-enabled JAX and Torch 2.6/cu124 plus matching PyG wheels and requires
a compatible NVIDIA driver. CPU installations do not acquire GPU support merely
by moving to a GPU machine: rerun setup with `--device gpu`. Runtime selection is
still `auto` within the installed model environment. No multi-GPU use is enabled.
On a cluster, use an allocated compute node for setup/inference tests. A login
node may have no visible GPU; `--device gpu` is available explicitly, but its real
test must run where a GPU is allocated. Setup neither submits jobs nor obtains
an allocation. Tests use four CPU cores by default (`--threads` changes this).

Failure returns a nonzero exit code, preserves logs/partial downloads and does
not replace that model's previous registration. The other model can still pass.
Rerun the same command after fixing the issue: complete verified downloads and
completed environments are reused; interrupted downloads restart. Changed worker
code/requirements use a different environment directory, preserving the old one.
Each environment has `installed-packages.txt`; each test has `summary.json` and
`worker.log`. An explicit generation fallback during a normal calculation is
not treated as successful setup. External downloads can fail due to outages or
Google Drive quotas; there is no automatic unverified mirror.

## Manual setup (advanced)

The core QCforever install uses RDKit and does not install Torch/JAX or weights.
Use **separate Python 3.11 environments** for the two models; their dependency
stacks are intentionally not combined in the core `requirements`.

The included `qcforever_model_workers` package does not import `qcforever` or
Gaussian/GAMESS. Its files have one purpose each:

* `worker.py`: raw-SDF request/response, CPU allocation and persistent lifetime.
* `ditmc.py`: upstream DiTMC graph preparation and checkpoint loading.
* `torsional.py`: load upstream TD once, reuse `sample_confs`, disable MM/energy.

No site-specific paths, external benchmark scripts, weights or upstream source
copies are required by these adapters. The caller supplies source/checkpoint
locations. The old `scripts/model_raw_worker.py` is only a compatibility shim.

## DiTMC

Obtain source and the Drugs aPE-B checkpoint from the authors' linked
[Zenodo release](https://doi.org/10.5281/zenodo.15489212), described in the
[official repository](https://github.com/ML4MolSim/dit_mc). Keep the extracted
`drugs/apeB/.hydra/config.yaml` and checkpoint tree together. The adapter expects
exactly one matching Drugs aPE-B tree. Source changes may require adapter changes;
the Zenodo snapshot, not an arbitrary future main branch, is the intended source.

From the QCforever checkout:

```bash
python3.11 -m venv .model-envs/ditmc
.model-envs/ditmc/bin/python -m pip install -r qcforever_model_workers/requirements/ditmc-cpu.txt
.model-envs/ditmc/bin/python -m pip install --no-deps .
```

The adapter imports upstream source using `--source`; do not run upstream
`pip install .` in this CPU environment: its metadata requests CUDA packages.
A C/C++ compiler and Python headers are needed for the upstream Cython module.
Compilation is serialized using a file lock in `--cache`, outside the source.

## Torsional Diffusion

Obtain [upstream source](https://github.com/gcorso/torsional-diffusion) at
`5f713b42d7000307655f272471014c6127ea59be` and download the `drugs_default`
weights linked in its README. Keep `model_parameters.yml` and `best_model.pt`
together. Do not use untrusted checkpoints: upstream loading may deserialize
Python objects.

```bash
python3.11 -m venv .model-envs/torsional
.model-envs/torsional/bin/python -m pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cpu
.model-envs/torsional/bin/python -m pip install torch-scatter==2.1.2 torch-cluster==1.6.3 --only-binary=:all: -f https://data.pyg.org/whl/torch-2.6.0+cpu.html
.model-envs/torsional/bin/python -m pip install -r qcforever_model_workers/requirements/torsional-cpu.txt
.model-envs/torsional/bin/python -m pip install --no-deps .
```

The adapter bootstraps upstream `generate_confs.py` with an empty input, then
calls its loaded `sample_confs` function directly for each exact requested count
(no upstream CLI 2x multiplication). MM flags remain off. A process-local stub
avoids importing upstream xTB evaluation, and likelihood population is disabled;
no source files are patched. RDKit embedding respects the requested thread count
and the seed of each request (including persistent additional batches).
Upstream lookup-table cache files are placed in the writable `--cache` directory.

## QCforever YAML

Use the installed worker command in each environment. Replace `/path/to/...`
with your own locations; they are examples, not required machine paths.

```yaml
generators:
  ditmc:
    persistent: true
    timeout_seconds: 1800
    command:
      - /path/to/ditmc-env/bin/qcforever-model-worker
      - --model
      - ditmc
      - --source
      - /path/to/ditmc-source
      - --checkpoint
      - /path/to/ditmc-checkpoints
      - --cache
      - /path/to/writable/ditmc-cache
      - --request
      - '{request}'
      - --output
      - '{output}'
  torsional_diffusion:
    persistent: true
    timeout_seconds: 1800
    command:
      - /path/to/td-env/bin/qcforever-model-worker
      - --model
      - torsional_diffusion
      - --source
      - /path/to/torsional-diffusion
      - --checkpoint
      - /path/to/drugs_default
      - --cache
      - /path/to/writable/td-cache
      - --request
      - '{request}'
      - --output
      - '{output}'
```

Device selection defaults to `auto`: the worker checks its own Torch/JAX
environment and uses an available CUDA GPU, otherwise CPU. These **CPU** install
recipes alone do not enable GPU inference. GPU environments require matching
CUDA-enabled framework/PyG wheels and a compatible NVIDIA driver; having
`nvidia-smi` available is not sufficient. Do not install both CPU and GPU stacks
into a shared environment. Use `device: cpu` to force CPU or `device: gpu` to
require GPU. GPU mode uses one persistent worker, preserves the scheduler's
GPU visibility mask, and logs actual parameter placement. It does not change
the xTB/PM6 allocation or run those engines on GPU. See the GPU execution
checks below.
Upstream models retain their own licenses/citation requirements. Their training
coverage does not guarantee good charged/radical conformers or multiplicity-aware
generation. Common validation and fallback still apply.

## Verification boundary

The one-command installer has local regression coverage for registration/YAML
precedence, failure preservation, checksum validation, safe extraction, and real
subprocess reuse with a small test worker. Wheel contents and the installed CLI
were checked in an isolated environment outside the source tree. Official TD
source/weight downloads, checksums and extraction were also exercised. These
checks alone do not establish an end-to-end installation. The fresh Linux CPU
test below additionally verifies dependency installation and real inference,
using previously obtained model assets. Setup requires actual model smoke tests on the
target system before registering either model.

The pinned major versions were read from the environments used for real GENKAI
inference. They are not a full transitive lockfile. A successful relocated-worker
test does not by itself certify a clean installation on every platform. Keep
resolver failures distinct from model inference failures, and record the installed
package list/source/checkpoint versions for reproducible comparisons.

On 2026-10-06, the relocated adapters passed real GENKAI CPU generation for
ethanol with requests 2 -> 1 -> 1, using two workers. Both models retained the
same worker PIDs across batches, with zero model reinitializations on the extra
batches. Both CPU requirement recipes passed pip's clean-resolution dry run.
Wheel build, isolated no-dependency installation and CLI loading outside the
source checkout also passed. New Linux GPU environments were subsequently tested
as described below. The official-download and reused-asset test boundaries are
documented separately below.

### Fresh CPU environments on 2026-10-06

On Ubuntu 20.04.2 x86_64, Python 3.11.5 created new isolated core and model venvs.
The unpublished QCforever source was installed with pip, rather than the public
GitHub main branch. Core imports and dependency checks passed. DiTMC and TD
dependencies were installed from their packaged recipes; both model environments
passed `pip check`, actual ethanol generation (1 then 1 candidates), CPU parameter
placement and same-PID reuse without additional model initialization. Both models
were registered by the installer.

The user requested reuse of previously obtained mnode source/weights instead of
waiting for the 1.89 GB Zenodo download. A test-only asset-provider substitution
was used; environment installation, inference tests and registration still ran
through the installed installer. No existing model Python environment was copied.
The transferred source/weight archive SHA256 was
`23b4463b170e07acb6f3c5a5fbb964f0eaab7c9293df55c8a5f5127bad51a68b`.

Without per-job model YAML, installed QCforever prepared ethanol candidates for
low, medium and high with test-only N=2 and four allocated CPUs. Final candidate
counts were 1, 2 and 2 respectively. Medium used TD directly; high exercised
the fallback route. Structure-validation thresholds were unchanged.
No xTB/PM6 executable was available in this test; native relaxation was not tested.

Logs, package lists and results are retained in the validation workspace
`qcforever_readme_validation_20261006/collected` next to the checkout.

### GPU validation on 2026-10-06

Dedicated Python 3.11 environments were built on mnode and inference was run
only inside Slurm allocations on gnode02 (RTX A6000, NVIDIA driver 550.78,
8 allocated CPUs and 1 GPU per job). Initial jobs 51218/51219 found two issues:
CUDA-enabled JAX raised when GPUs were explicitly hidden, and the TD recipe
lacked its upstream `spyrmsd` import. Auto selection now handles an explicit
visibility mask before importing the framework and handles JAX's no-device
error without hiding other driver/library errors. TD includes `spyrmsd==0.9.0`.

Retests 51221 (DiTMC) and 51222 (TD) verified actual model parameters on `cuda:0`,
generation of 2 then 1 ethanol conformers, identical worker PID across batches,
zero additional model-initialization time, `auto` selecting CPU with GPUs hidden,
and the public pipeline using one GPU model worker. Both execution tests passed.
This CPU check verifies device selection, not a separate full CPU inference run.

Logs, raw SDFs, summaries and complete environment version lists were retained
in the local validation workspace `gpu_validation_20261006/collected_retry1`
(next to the QCforever checkout), under `results_retry1/<model>/summary.json`.
The tested source archive SHA256 is
`7a9e7bdcbf9faedfd986dab93bdbe46781239507aa13cc1a79ce5fd273870fc1`.
This was an isolated manual-environment validation, not an end-to-end execution
of `install-conformer-models`. The runtime GPU fixes and missing dependency are
also included in that installer's current source/recipe.
