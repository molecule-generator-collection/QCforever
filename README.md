# QCforever

![Robot2+PC](https://user-images.githubusercontent.com/46772738/188896764-65ab12c1-3cc9-421d-8d87-ed33c932380a.png)

QCforever (https://doi.org/10.1021/acs.jcim.2c00812, https://doi.org/10.1002/jcc.70017) is a wrapper of Gaussian (https://gaussian.com) or GAMESS (https://www.msg.chem.iastate.edu/gamess/). 
To compute obsevable properties of a molecule through quantum chemical computation (QC),
multi step computation is demanded. 
QCforever automates this process and calculates multiple physical properties of molecules simultaneously.
Don't you belive the QC?
Then, you can optimize functional parameters of density functional theory through Bayesian optimization (https://doi.org/10.1021/acs.jctc.3c00764).
Grey-box optimisation (LAQA (https://doi.org/10.1021/acs.jctc.1c00301)) can also be used to obtain energetically favourable molecular conformations.


## Requirements

1. [Gaussian](https://gaussian.com)==16
2. [GAMESS(sockets)](https://www.msg.chem.iastate.edu/gamess/)==30 SEP 2022 (R2)
3. [Python](https://www.anaconda.com/download/)==3.11
4. [rdkit-pypi](https://anaconda.org/rdkit/rdkit)==2023.09.1
5. [bayesian-optimization](https://github.com/bayesian-optimization/BayesianOptimization)==1.4.3
6. [psutil](https://github.com/giampaolo/psutil)
7. [basis-set-exchange](https://pypi.org/project/basis-set-exchange/)=0.12

## Optional

1. [xtb](https://github.com/grimme-lab/xtb/tree/v6.6.1)==6.6.1

## How to use

### Install

```bash
pip install --upgrade git+https://github.com/molecule-generator-collection/QCforever.git@feature/conformer-search
```

### Install (optional for detailed conformation search)

To enable DiTMC and Torsional Diffusion for `optconf_medium` and `optconf_high`:

```bash
install-conformer-models
```

Requires Linux x86_64, Python 3.11 and a C/C++ compiler. This command creates
separate model environments, downloads the official source/weights, and registers
each model only after an ethanol generation/reuse test passes. No per-job YAML
is needed. Existing Python environments are not changed. The first setup needs
several GB of downloads/disk space (the DiTMC archive alone is about 1.9 GB).
Setup uses a visible NVIDIA GPU when available, otherwise CPU; use `--device cpu`
or `--device gpu` to choose explicitly. On a cluster, run the generation tests
inside an appropriate CPU/GPU allocation, not on a login node.
Run `install-conformer-models --help` for setup options.

### Example

Example codes for QCforever (gaussian_main.py and gamess_main.py) are prepared.

For Gaussian
```bash
python gaussian_main.py input_file
```
For GAMESS
```bash
python gamess_main.py input_file
```

### Usage

```python
from qcforever.gaussian_run import GaussianRunPack

# Then make an instance (here is test).
test = GaussianRunPack.GaussianDFTRun(Functional, basis_set, ncore, option, input_file, solvent='water', restart=False)
# And excuse Gaussian as the followign.
outdic = test.run_gaussian()
```

- ***Functional*** is to specify the functional in the density functional theory (DFT).
Currently supported functionals are the following: `BLYP`, `B3LYP`, `X3LYP`, `LC-BLYP`, `CAM-B3LYP`, 
and [KTLC- series](https://doi.org/10.1021/acs.jctc.3c00764) (`KTLC-BLYP-BO`, `KTLC-wPBE-BO`, and so on).

- ***basis_set*** is to specify the basis set.
  Current QCforever supports the following basis set:
  `LANL2DZ`, `STO-3G`, `3-21G`, `6-31G`, `6-311G`, `3-21G*`, `3-21+G*`, `6-31G*`, `6-311G**`, `6-31G**`,
  `6-31+G*`, `3-21G`, `3-21G`, `6-31+G**`, `6-311+G*`, and `6-311+G**`

- ***ncore*** is an integer to specify the number of core for QC with Gaussian.

- Memory is estimated from the molecular basis functions, number of electrons,
  calculation options, and available system memory. The estimate is written to
  the Gaussian or GAMESS input automatically. To override it, set a value such
  as `test.mem = "4GB"` before calling `run_gaussian()` or `run_gamess()`.

- The total QCforever wall-clock time can be limited in seconds with
  `test.timejob = 24 * 60 * 60`.  When this limit is reached, QCforever stops
  the active external process tree (including xTB, Gaussian, or GAMESS) and
  returns a result whose `log` is `"timeout"`.  The existing `timexe`
  (Gaussian) and `timeexe` (GAMESS) settings remain independent per-calculation
  limits; GAMESS continues to write `timeexe` into its input file.

- After the workflow finishes (normally, with an error, or by `timejob`), the
  calculation directory is cleaned automatically. Gaussian keeps `.chk` and
  `.fchk` restart files (or the selected geometry `.xyz` output when
  `restart=False`) together with `.com`/`.gjf` input and `.log`/`.out` output
  files. GAMESS first moves matching `.dat` restart files from SCR/USERSCR into
  the calculation directory and keeps them together with `.inp` input and
  `.log`/`.out` output files.
  If `pklsave=True`, the explicitly requested `.pkl` result is also retained.

- ***option*** is a string for specifying molecular properties as explained later.

- ***input_file*** is a string to specify the input file.
  QCforever accepts a sdf, xyz, Gaussian chk, or a Gaussian fchk file.

- ***solvent*** is to include the solvent effect through PCM.
  The default value is `None`, in vacuo. The legacy value `"0"` is also accepted.

- ***restart*** is to control to save molecular information as fchk or xyz.
  The Default value is True that means molecular information is saved as a Gaussian fchk file (electronic structure is also saved.),
  otherwise molecular information is saved as a xyz files (electronic structure is not saved).

In this examples, obtained results are saved as python dictionary style.

### Options

By specifying the molecular properties you want as an option variable string,
QCforever automatically calculates them.
The variable string must consisted of the following options,
seperated with more than one space:

Following options are currently available:

| Option name | Description | Gaussian | Gamess |
|---|---|---|---|
|symm| just specify symmetry of a molecule.|:white_check_mark:||
|volume| Compute the volume (in cm**3/mol) of a molecule.|:white_check_mark:||
|opt| just perform geometry optimization of a molecule.|:white_check_mark:|:white_check_mark:|
|nmr| NMR chemical shift (ppm to TMS) of each atom is computed. `nmr=/absolute/path/to/nmr.dat` also compares it with a two-column `position intensity` reference spectrum.|:white_check_mark:||
|uv| absorption wavelengths (nm) (for spin allowed states) are computed. If you add an absolute path of text file that specify peak positions and their peak like (uv=/ab/path2uv.dat), it is possible to the similarity and dissimilarity between the target and reference. |:white_check_mark:|:white_check_mark:|
|energy| SCF energy (in Eh) is printed.|:white_check_mark:|:white_check_mark:|
|cden| charge and spin densities on each atom are computed.|:white_check_mark:|:white_check_mark:|
|homolumo| HOMO/LUMO gap (Eh) calculation|:white_check_mark:|:white_check_mark:|
|dipole| dipole moment of a molecule|:white_check_mark:|:white_check_mark:|
|polar| Dipole polarizability (in 10**-24 cm**3) of a molecule|:white_check_mark:||
|deen| decomposition energy (in eV) of a molecule|:white_check_mark:||
|stable2o2| stability to oxygen molecule|:white_check_mark:||
|vip| vertical ionization potential energy (in eV)|:white_check_mark:|:white_check_mark:|
|vea| vertical electronic affinity (in eV)|:white_check_mark:|:white_check_mark:|
|aip| adiabatic ionization energy (in eV)|:white_check_mark:|:white_check_mark:|
|aea| adiabatic electronic affinity (in eV)|:white_check_mark:|:white_check_mark:|
|fluor| wavelength (in nm) of fluorescence are computed. if you want to specify the state that emits fluorescence, you can specify the index of state like “fluor=#” (# is an integer, default is “fluor=1”)|:white_check_mark:|:white_check_mark:|
|nac| Max and rms values of non adiabatic coupling vector (use with fluor option) |:white_check_mark:||
|tadf| Compute the energy gap (in Eh) between minimum in the spin allowed state and the spin forbidden state.|:white_check_mark:||
|freq| Compute the variable related to the vibrational analysis of a molecule. `freq=IR.dat,Raman.dat` compares the calculated IR and Raman spectra with two-column `position intensity` reference files. For Gaussian, an optional third NMR file (`freq=IR.dat,Raman.dat,NMR.dat`) also enables NMR and compares all three spectra. Before comparison, peaks within 1.0 cm-1 (IR/Raman) or 0.01 ppm (NMR) are merged; intensities are summed and positions are averaged.|:white_check_mark:|:white_check_mark:|
|pka| Compute the energy gap (in Eh) between deprotonated (A-) and protonated (AH) species. The hydrogen atom whose Mulliken charge is the biggest in the system is selected as a protic hydrogen.|:white_check_mark:||
|stable| try to find a stable structure when the negative frequency is detected.|:white_check_mark:||
|optspin| try to find a suitable spin multiplicity.|:white_check_mark:||
|optconf| try to find a stable molecular conformation with PM6 (Gaussian16).|:white_check_mark:|:white_check_mark:|

### Configurable conformer search (development)

`optconf=xtb` / `optconf=pm6` default to `optconf_low` (ETKDGv3 and MMFF94s).
Add `optconf_medium` or `optconf_high` for learned-generator fallback routes.
Optional `job.conformer_config = 'conformer.yaml'` overrides
[`defaults.yaml`](qcforever/conformer_search/defaults.yaml).

#### Relaxation options

Both backends default to graybox SH, stopping when 20% of the initial relaxation
pool has converged. This affects conformer selection, not the subsequent DFT job.

| Option (after either `optconf=xtb` or `optconf=pm6`) | Behavior |
|---|---|
| Omitted, or `laqa` | 20% convergence target |
| `laqa=80` | 80% convergence target |
| `laqa=100` | Target all candidates; failures and safety limits can prevent completion |
| `laqa=off` | Continuous relaxation of all candidates |

`laqa` is the user-facing option; the actual algorithm defaults to `sh`.
An explicit percentage overrides the YAML fraction, not its algorithm.
For example:

```yaml
relaxation:
  algorithm: sh                 # sh, sr, or laqa
  convergence_fraction: 0.20
  first_interval: 1
  subsequent_interval: 10
  parallel_candidates: auto     # Up to nproc; 1 selects sequential execution
```

Relaxation uses one core per candidate, up to `nproc` concurrent candidates
within the CPU allocation. CPU model generation uses four cores per worker.
Gaussian memory is per worker; account for all concurrent workers.

The convergence target is `ceil(fraction * initial_candidates)`, after
generation/MM/filtering. Its denominator does not shrink after failures.
Only native convergence counts; final geometry/stereo checks remain separate.
Selection and stopping occur between complete batches, so costs and converged
counts can exceed the target. Parallel and sequential searches need not select
the same final structure.

SH/SR start with 20 evaluations per initial candidate and add another 20 as
needed, re-admitting unfinished candidates without resetting their progress.
Within each round, SH retains the best half and SR rejects the worst candidate
at each completed stage; residual budget goes to unfinished candidates ranked
by energy. LAQA instead selects the lowest current `E/N - F²/(2ΔF)` scores.
PM6 uses mean atomic force; xTB uses ANC gradient norm divided by `sqrt(N)`,
recorded as the approximate `laqa_norm` policy. No future energy enters selection.

PM6 graybox uses Gaussian16 PM6-RFO with checkpoint restart, tight optimization
and SCF convergence, and cumulative cycle limits 1,11,21,… .
Both `g16` and `formchk` must be on PATH. `laqa=off` retains the legacy
Gaussian optimizer; set `relaxation.pm6_optimizer: rfo` for a continuous-RFO
comparison.

xTB graybox keeps each optimizer process alive using POSIX stop/resume signals;
it does not restart from coordinates. Defaults remain GFN2-xTB, `normal`,
accuracy 1.0. Paused processes consume memory. The default 8192 MB pool RSS guard
is a sampled safeguard, not a scheduler memory reservation. Use scheduler-owned
jobs and sufficient memory. Live state cannot resume after QCforever exits.
On Windows, use `laqa=off`.

Results and scalar histories are saved under `conformer_search/electronic/`.
`audit.json` distinguishes converged, failed, limited and unfinished candidates,
records actual evaluations and elapsed time, and reports unmet convergence
targets. Only converged structures proceed to final selection.

#### Updating learned generators

Explicit YAML model settings override registered models. To update an existing
Torsional Diffusion installation, rerun
`install-conformer-models --models torsional_diffusion --directory <same-directory>`.
Setup tests the new private environment before replacing its registration.
Completed lookup tables are read concurrently; only first-time construction
requires exclusive access. `model_execution.json` records initialization and
generation timings. Setup checks basic inference, not validity for every molecule.

## License

This package is distributed under the MIT License.

## Contact

- Masato Sumita (masato.sumita@riken.jp)
