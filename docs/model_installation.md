# Optional model workers (Linux CPU)

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
.model-envs/ditmc/bin/python -m pip install -r requirements/ditmc-cpu.txt
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
.model-envs/torsional/bin/python -m pip install -r requirements/torsional-cpu.txt
.model-envs/torsional/bin/python -m pip install --no-deps .
```

The adapter bootstraps upstream `generate_confs.py` with an empty input, then
calls its loaded `sample_confs` function directly for each exact requested count
(no upstream CLI 2x multiplication). MM flags remain off. A process-local stub
avoids importing upstream xTB evaluation, and likelihood population is disabled;
no source files are patched. RDKit embedding respects the requested thread count.
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

CPU is the validated default. GPU setup is not covered by these CPU recipes.
Upstream models retain their own licenses/citation requirements. Their training
coverage does not guarantee good charged/radical conformers or multiplicity-aware
generation. Common validation and fallback still apply.

## Verification boundary

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
source checkout also passed. A completely new ML environment plus actual
inference in that new environment has **not** been tested yet.
