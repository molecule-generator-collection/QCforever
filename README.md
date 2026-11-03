# QCforever

![Robot2+PC](https://user-images.githubusercontent.com/46772738/188896764-65ab12c1-3cc9-421d-8d87-ed33c932380a.png)

QCforever (https://doi.org/10.1021/acs.jcim.2c00812, https://doi.org/10.1002/jcc.70017) is a wrapper of Gaussian (https://gaussian.com), GAMESS (https://www.msg.chem.iastate.edu/gamess/), or ORCA (https://www.faccts.de/orca/).
To compute obsevable properties of a molecule through quantum chemical computation (QC),
multi step computation is demanded. 
QCforever automates this process and calculates multiple physical properties of molecules simultaneously.
Don't you belive the QC?
Then, you can optimize functional parameters of density functional theory through Bayesian optimization (https://doi.org/10.1021/acs.jctc.3c00764).
Grey-box optimisation (LAQA (https://doi.org/10.1021/acs.jctc.1c00301)) can also be used to obtain energetically favourable molecular conformations.


## Requirements

1. Quantum chemistry software: [Gaussian](https://gaussian.com/)==16, [GAMESS (sockets)](https://www.msg.chem.iastate.edu/gamess/)==30 SEP 2022 (R2), or [ORCA](https://www.faccts.de/orca/)>=6.1 (validated with 6.1.1)
2. [Python](https://www.anaconda.com/download/)==3.11
3. [RDKit](https://pypi.org/project/rdkit/)>=2026.3.6 (matches the bundled KTLC prediction models' training version)
4. [bayesian-optimization](https://github.com/bayesian-optimization/BayesianOptimization)==1.4.3
5. [psutil](https://github.com/giampaolo/psutil)
6. [basis-set-exchange](https://pypi.org/project/basis-set-exchange/)>=0.12
7. NumPy>=1.22.0
8. PyYAML>=6
9. LightGBM==4.7.0 (see the [official installation guide](https://github.com/lightgbm-org/LightGBM/blob/main/python-package/README.rst#install-from-pypi) for macOS prerequisites)

## Optional

1. [xtb](https://github.com/grimme-lab/xtb/tree/v6.6.1)==6.6.1

## How to use

### Install

```bash
pip install --upgrade git+https://github.com/molecule-generator-collection/QCforever.git
```

<details>
<summary>Optional conformer model installation</summary>

No additional models are needed for `conf_gen=low` or `conf_gen=legacy`.
For `conf_gen=medium`, install Torsional Diffusion:

```bash
install-conformer-models --models torsional_diffusion --device cpu
```

For `conf_gen=high`, install both DiTMC and Torsional Diffusion:

```bash
install-conformer-models --device cpu
```

Requires Python 3.11 and several GB of free disk space. DiTMC also needs a
C/C++ compiler and Python development headers. On macOS, install Command Line
Tools with `xcode-select --install` if absent; on Linux, use the site's compiler
module or distribution build tools.

The installer downloads dependencies/weights into separate environments and
registers models after a small ethanol generation/reuse test passes.
The default location is `~/.local/share/qcforever/conformers`;
use `--directory /path/to/new-directory` to change it.
If setup fails, inspect the printed log and rerun the same command.

Linux x86_64 supports CPU or NVIDIA GPU (`--device gpu`, CUDA-12-compatible driver).
On clusters, run setup in a compute allocation.
Apple Silicon macOS uses CPU only: Torsional Diffusion's setup test has passed
on one Mac; DiTMC validation is pending. Intel Macs are not supported.
See `install-conformer-models --help` for other installation options.

</details>

### Example

Example codes for QCforever (gaussian_main.py, gamess_main.py, and orca_main.py) are prepared.

For Gaussian
```bash
python gaussian_main.py input_file
```
For GAMESS
```bash
python gamess_main.py input_file
```
For ORCA
```bash
python orca_main.py molecule.sdf "opt energy homolumo dipole" --nproc 4 --memory 4GB
```

### ORCA

Install ORCA separately and put its `orca` executable on `PATH`, or set
`job.executable` to its full path. For parallel runs, install the MPI runtime
required by your ORCA distribution. QCforever invokes the ORCA driver directly
using its absolute path; `%pal` controls the workers.

```python
from qcforever.orca_run import OrcaRunPack

job = OrcaRunPack.OrcaDFTRun(
    'B3LYP', 'def2-SVP', 4,
    'optconf=(conf_gen=low, conf_relax=xtb) opt energy homolumo dipole uv',
    'molecule.sdf', solvent='water',
)
job.mem = '4GB'             # Total budget, shared by all ORCA processes
job.timexe = 60 * 60        # Limit for each ORCA invocation (seconds)
job.timejob = 24 * 60 * 60  # Limit for the complete workflow, including xTB
output = job.run_orca()
```

Supported options are `opt` (`opt=loose/tight/verytight` also accepted), `energy`,
`homolumo`, `dipole`, `cden`, `polar`, `vip`, `vea`, `aip`, `aea`, `uv`,
`fluor`/`fluor=N`, `freq`, `nmr`, and `optconf`. Unsupported options, other KTLC
functionals, and `para_functional` overrides return `log="error"` before launch.
ORCA accepts SDF, XYZ, and Gaussian formatted checkpoints (`.fchk`/`.fch`);
see [FCHK input](#gaussian-fchk-input-for-gamess-and-orca). XYZ comments can specify `charge multiplicity`
(otherwise neutral singlet). `SpecTotalCharge` and `SpecSpinMulti` override them.
Use ORCA functional/basis keywords. `B3LYP` selects the Gaussian-compatible
`B3LYP/G` convention. B3LYP/G, BLYP, LC-BLYP and CAM-B3LYP use their equivalent
`LibXC(...)` forms throughout the workflow to support excited-state derivatives.
Other native ORCA keywords, including explicit `LibXC(...)`, pass through.
Results from different programs still depend on grids, approximations, and
convergence settings and need not be numerically identical.

`solvent` selects ORCA CPCM; `None`, `"0"` and `"none"` select gas phase.
Use solvent names supported by ORCA. `mem` is a total budget; QCforever assigns
75% divided by `nproc` to `%maxcore` (MB per process). The remaining 25% provides
headroom; ORCA memory use is not a strict operating-system limit.

ORCA `optconf` supports **xTB only**, including `conf_gen=legacy` with xTB.
Bare `optconf` and groups without `conf_relax` default to xTB for this backend.
Explicit PM6 requests fail before any calculation. Install the separate `xtb`
executable for conformer searches. The existing generation, LAQA settings,
selected structure, and retained histories are reused. A failed search returns
an error; it does not silently continue with the original geometry.

ORCA supports `KTLC-BLYP-BO` and `KTLC-PBE-BO` as functional names. For example:

```python
job = OrcaRunPack.OrcaDFTRun(
    'KTLC-PBE-BO', 'def2-SVP', 4, 'optconf opt energy homolumo uv', 'molecule.sdf',
)
output = job.run_orca()
```

The Bayesian search holds the input geometry fixed, using the selected xTB
conformer if `optconf` was requested. It evaluates the parent and its one-electron
oxidized state at each trial μ and minimizes the squared residual
`[E(N−1) − E(N) + ε_HOMO(N)]²`. The basis, solvent and vertical ion multiplicity
follow the workflow settings. The search targets an absolute residual of 0.01 eV
with a budget of 52 trials in the range 0.01–1.0 bohr⁻¹. If the target is not
reached, the best valid trial is used. If all trials fail, the workflow returns an
error. The selected μ appears in `functional_param`, and the residual and trials
appear in `functional_parameter_optimization`. Trial inputs, logs and restart
files are temporary and are removed when tuning finishes or fails. Rerunning the
workflow repeats the tuning while completed property stages remain eligible for
reuse.

The selected μ is fixed for the following structure optimizations and property
calculations. ORCA receives it as `%method RangeSepMu`. For LC-BLYP, the LibXC
short-range exchange parameter `_omega` is also set to the same value. Excited
state optimization with `KTLC-PBE-BO` uses numerical gradients (`NumGrad`) because
ORCA 6.1 lacks the required analytic derivative for its native LC-PBE exchange.

ORCA uses `MOREAD` to start vertical ion calculations, UV and the `fluor`
optimization from the parent orbitals. Each KTLC VIP trial starts from the
parent orbitals calculated at the same μ. Relaxed ion and emission steps start
from the corresponding vertical ion and relaxed excited-state orbitals,
respectively. If a source `.gbw` is unavailable or ORCA rejects the orbital guess,
the stage uses ORCA's default initial guess. Charge and multiplicity are set
for each target calculation, and SCF is converged again.

Output conventions:

- `energy` and its alias `Energy` are `[energy in Eh, <S²>-S(S+1)]`.
  `homolumo` is in **eV**, a scalar for restricted calculations or two spin
  gaps for unrestricted calculations. `MO`, `Alpha_MO`, and `Beta_MO` contain
  occupied/virtual orbital energies in Eh.
- `dipole` is `[X, Y, Z, total]` in Debye, at the center of mass in ORCA's
  Cartesian frame. The quadrupole keys use the existing primitive/traceless
  conventions and Debye Angstrom units, also at the center of mass.
  `polar_tens`, `polar_iso`, and `polar_aniso` are in 10^-24 cm³ (Angstrom³).
- `vip`, `vea`, `aip`, and `aea` are `[energy difference in eV, ion spin
  contamination]`. Adiabatic requests also return the corresponding vertical
  result. The ion is a doublet for a singlet input, otherwise its multiplicity
  is reduced by one. Override individual states with
  `job.ion_multiplicities = {'vip': 3, 'aip': 3}`, for example.
  Adiabatic calculations optimize the ion; use `opt` to optimize the reference
  molecule as well.
- `uv` and `fluor` are `[wavelengths in nm, oscillator strengths]`. `fluor`
  contains the vertical transition spectrum at the relaxed excited-state
  geometry; `fluor_state` identifies the followed state's root, wavelength and
  oscillator strength. Full TDDFT (`TDA false`) is used, with `job.nstates=20`
  by default; the requested fluorescence root can increase this count.
- `freq` uses **NumFreq** plus polarizabilities to obtain both IR (km/mol) and
  Raman activities (Angstrom⁴/amu), including for B88-based functionals.
  This costs more than an analytic Hessian. The arrays include all 3N modes,
  with zero intensities for omitted translations/rotations. `freqmode` is
  `[X, Y, Z]`, each a mode-by-atom array. `Ezp`, `Et`, `E_enth`, and `E_free`
  are total energies in Eh; `Ei` is the total thermal correction including ZPE
  in kcal/mol and `Si` is entropy in cal/(mol K). ORCA's default QRRHO entropy
  treatment is retained. Heat capacities are not returned.
- `nmr_shielding` contains absolute shieldings in ppm. `nmr` uses TMS chemical
  shifts for H/C/Si and absolute shieldings for other elements. If needed,
  an additional ORCA calculation optimizes TMS and computes its shieldings
  with the same functional, basis and solvent. This can be costly. Supply
  compatible reference shieldings with `job.nmr_reference = {'H': 31.5,
  'C': 185.0, 'Si': 330.0}` to avoid the reference calculation; these numbers
  are illustrative, not universal reference values. Actual references are
  returned in `nmr_reference`.

`uv=reference.dat`, `nmr=reference.dat`, and
`freq=IR.dat,Raman.dat[,NMR.dat]` compare against two-column reference spectra.
IR/Raman/NMR use the shared `Spectrum_similarity` output; UV uses the existing
`Similarity/Dissimilarity` and `Wasserstein_distance` keys. Relative reference
paths are resolved before entering the calculation directory.

Each run uses a directory named after the input stem in the current directory.
Existing results are inspected and verified stages are reused as described in
[Calculation reuse](#calculation-reuse-gaussian-gamess-and-orca).
Inputs, logs, XYZ geometries, `.hess` frequency data, property text and citations
are retained. `restart=True` (default) also retains `.gbw` orbitals;
`restart=False` retains geometry instead. This flag controls retention, not
optimizer resumption. The separate reuse cache retains restart data needed to
reconstruct completed stages even when `restart=False`. `pklsave=True`
preserves the result pickle. Errors and timeouts return `log="error"` or
`log="timeout"`, partial results, and the failed/active job's diagnostic log.
Independent molecules should run in separate Python processes.

ORCA references:
[LibXC functional mappings](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/DensityFunctionalTheory.html#libxc-functional-library),
[electric properties](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/electric.html),
[TDDFT and excited-state optimization](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/tddft.html),
[Raman spectroscopy](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/vibrations.html#raman-spectra).

To run the optional installed-ORCA regression test:
`QCFOREVER_TEST_ORCA=1 python -m pytest -q test/test_orca.py`.

### GAMESS scratch handling

GAMESS scratch paths are resolved separately for `SCR` and `USERSCR`. The paths
reported by `rungms` in its execution log take precedence, including paths chosen
by a site wrapper, process ID, or batch scheduler. Standard output and standard
error are saved to the calculation's `.log` file while GAMESS runs, so partial
output and reported paths remain available after a timeout.

Before launch, QCforever uses simple unconditional assignments in `rungms`,
followed by the environment inherited by Python. Runtime expressions and shell
branches are left to the launcher; QCforever does not execute the script to
guess its settings. For installations whose `rungms` uses environment variables,
set these before starting QCforever (for example, in `.bashrc`):

```bash
export SCR=/path/to/scr
export USERSCR=/path/to/userscr
```

If either path is unspecified, QCforever uses `.gamess_scratch/scr` or
`.gamess_scratch/userscr` inside the calculation directory. All selected
directories are created if necessary and passed to `rungms` in its environment.
Creation or permission errors are reported rather than changing the selected
location. A launcher may override these initial paths. QCforever uses its
reported `USERSCR` to collect `.dat` restart files and reads all of the job's
logs to find scratch files for final cleanup, including interrupted launches.
If the launcher exits unsuccessfully or the required `.dat` file is missing,
`run_gamess()` returns `log="error"` with the failed job name, launcher exit code,
log path, and diagnostic log tail. A site wrapper may report exit code zero
even when its backend failed, so a missing `.dat` is still treated as an error.
For `limit: stacksize: Can't set limit`, check the soft/hard stack limits in the
same job environment and the site's launcher settings with the center's
administrator; changing scratch paths does not fix a stack-limit failure.

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
and [KTLC- series](https://doi.org/10.1021/acs.jctc.3c00764) (`KTLC-BLYP-BO`, `KTLC-wPBE-BO`, `KTLC-BLYP-Pred`, and so on).

<details>
<summary>KTLC-BLYP-Pred: prediction model and usage details</summary>

`KTLC-BLYP-Pred` (Gaussian, SDF/XYZ input) uses the mean of five bundled
RDKit 2D + Morgan + LightGBM models as a fixed LC-BLYP mu, without BO or an
additional geometry optimization. Existing `opt`/`optconf` options still apply.
The result's `functional_parameter_prediction` contains only `model_version`,
`mu_mean` and `mu_std` (bohr^-1; population SD across five models, not a confidence
interval). `functional_param` records the actual mu rounded to four decimal places.
The original models are used unchanged (60/20/20 train/validation/test per fold).
Training conditions: neutral closed-shell molecules, MW <= 1000, gas phase,
LC-BLYP/6-31G*. Prediction does not guarantee Koopmans compliance. Descriptors
were computed with RDKit 2026.03.6; older versions lacking a required descriptor
must be updated. No separate model installation is needed.

</details>

- ***basis_set*** is to specify the basis set.
  Current QCforever supports the following basis set:
  `LANL2DZ`, `STO-3G`, `3-21G`, `6-31G`, `6-311G`, `3-21G*`, `3-21+G*`, `6-31G*`, `6-311G**`, `6-31G**`,
  `6-31+G*`, `3-21G`, `3-21G`, `6-31+G**`, `6-311+G*`, and `6-311+G**`

- ***ncore*** is an integer to specify the number of core for QC with Gaussian.

- Memory is estimated from the molecular basis functions, number of electrons,
  calculation options, and available system memory. The estimate is written to
  the Gaussian, GAMESS or ORCA input automatically. `CalcFC`, `CalcAll` and
  `RecalcFC=N` in ground- or excited-state optimization apply the same 1.5x
  memory allowance as `freq`. This allowance is applied once even if both
  optimizations request Hessians or `freq` is also requested. `ReadFC`/`RCFC`
  do not request this extra allowance. Whole-unit rounding and the available
  RAM limit still apply. To override automatic memory, set a value such as
  `test.mem = "4GB"` before calling `run_gaussian()`, `run_gamess()` or `run_orca()`.

- The total QCforever wall-clock time can be limited in seconds with
  `test.timejob = 24 * 60 * 60`.  When this limit is reached, QCforever stops
  the active external process tree (including xTB, Gaussian, or GAMESS) and
  returns a result whose `log` is `"timeout"`.  The existing `timexe`
  (Gaussian) and `timeexe` (GAMESS) settings remain independent per-calculation
  limits; GAMESS continues to write `timeexe` into its input file.

- After the workflow finishes (normally, with an error, or by `timejob`), the
  calculation directory is cleaned automatically. Gaussian keeps `.chk` and
  `.fchk` restart files together with `.com`/`.gjf` input and `.log`/`.out` output
  files. With `restart=False`, Gaussian keeps available `.xyz` output and removes
  `.chk` and `.fchk` files from the job directory, reuse cache and job history
  after both successful and failed runs. GAMESS first moves matching `.dat`
  restart files from SCR/USERSCR into the calculation directory and keeps them together with `.inp` input and
  `.log`/`.out` output files.
  If `pklsave=True`, the explicitly requested `.pkl` result is also retained.

- ***option*** is a string for specifying molecular properties as explained later.

- ***input_file*** is a string to specify the input file.
  QCforever accepts a sdf, xyz, Gaussian chk, or a Gaussian fchk file.
  Gaussian binary `.chk` input is supported only by the Gaussian backend;
  GAMESS and ORCA accept Gaussian `.fchk`/`.fch` as geometry/state input.

- Set `test.SpecTotalCharge` and/or `test.SpecSpinMulti` before running to
  override the input file's charge and spin multiplicity. An unspecified value
  retains the input file's value. After `optconf`, the selected SDF supplies
  geometry and bonds; its formal charge and radical count do not overwrite the
  electronic state used for the search. This also preserves a multiplicity
  selected by Gaussian's `optspin` for subsequent parameter optimization.

- ***solvent*** is to include the solvent effect through PCM.
  The default value is `None`, in vacuo. The legacy value `"0"` is also accepted.

- ***restart*** is to control to save molecular information as fchk or xyz.
  The Default value is True that means molecular information is saved as a Gaussian fchk file (electronic structure is also saved.),
  otherwise molecular information is saved as a xyz files (electronic structure is not saved).

In this examples, obtained results are saved as python dictionary style.

### Gaussian FCHK input for GAMESS and ORCA

GAMESS and ORCA accept Gaussian formatted checkpoints (`.fchk` or `.fch`,
case-insensitive) directly, without running Gaussian, `formchk` or `unfchk`.
The imported data are the element identities from `Atomic numbers`,
`Current cartesian coordinates` converted from Bohr to Angstrom, `Charge`, and
`Multiplicity`. Element identities are not inferred from ECP-reduced
`Nuclear charges`. Invalid/missing structural fields and unsupported ghost/dummy
centers produce an error before launching the target engine.

```python
from qcforever.gamess_run import GamessRunPack
from qcforever.orca_run import OrcaRunPack

gamess = GamessRunPack.GamessDFTRun(
    'B3LYP', 'STO-3G', 1, 'opt energy', '/path/to/molecule.fchk',
)
gamess_result = gamess.run_gamess()

orca = OrcaRunPack.OrcaDFTRun(
    'B3LYP', 'def2-SVP', 1, 'opt energy homolumo', '/path/to/molecule.fchk',
)
orca_result = orca.run_orca()
```

The destination calculation uses the functional, basis/ECP, solvent and options
specified on the GAMESS/ORCA job. The Gaussian method, basis/ECP, solvent,
isotope masses, constraints, energies, density and molecular orbitals are not
transferred. The target engine starts a new electronic calculation at the
imported geometry. Add `opt` to optimize it at the new level; without `opt`,
the supplied geometry is kept for a single-point calculation. The source fchk
is left unchanged and no sibling `.chk` is created. Binary `.chk` files must
first be converted to formatted text using Gaussian's `formchk` separately.

`SpecTotalCharge` and `SpecSpinMulti` override the imported electronic state in
the usual way. With `optconf`, the resolved geometry, charge and multiplicity
are exported to `fchk_geometry.xyz` inside the calculation directory and passed
to the existing search. Bonds are reconstructed as for any XYZ input; the SDF
bond graph and stereochemical annotations are not present in FCHK. ORCA still
supports only `conf_relax=xtb`. Cross-program input imports geometry/state;
Gaussian numerical results are not substituted for GAMESS/ORCA results by reuse.

### Calculation reuse (Gaussian, GAMESS and ORCA)

All three backends enable `job.reuse = True` by default. When the calculation
directory already exists, QCforever preserves its previous contents under
`.qcforever/history/` and plans the requested workflow again. At each stage it
compares the generated input and any external restart dependencies with recorded
calculations. Only a completed, converged calculation with verified required
outputs is reused. Failed or interrupted stages are rerun; their completed
predecessors can still be reused. Incomplete optimizer checkpoints are not
automatically resumed within a stage.

The comparison retains the geometry and atom order, charge, multiplicity,
functional and its parameters, basis/ECP, solvent, excited-state settings,
constraints, and numerical convergence controls present in the input. Core
counts, memory allocation and execution time limits are excluded. The selected
GAMESS launcher version and the available executable's path/size/modification
time are also checked. Checkpoint/geometry dependencies and saved outputs use
SHA-256 checksums. A changed or missing required output causes that stage to be
computed again. Changed upstream checkpoints can conservatively invalidate
downstream stages even when they represent similar numerical results.

Gaussian Link1 calculations are run and recorded one link at a time, then their
logs are combined for the existing property parser. Adding a frequency, NMR or
other link can reuse an unchanged completed prefix rather than rerunning its
optimization. The existing per-invocation `timexe` budget covers the entire
sequence of links. GAMESS reuses separate ground-state, electric-moment, TDDFT, ion and vibrational
jobs. ORCA records the optimization/energy, electronic properties, polarizability,
NMR, TDDFT, excited-state optimization, ion and vibrational stages separately.
Completed conformer searches can also be reused when the source structure,
electronic state, generation/relaxation settings and relaxation executable match.
Their selected structure and retained histories are restored together.

For example, running `"opt energy uv"` after `"opt energy"` can reuse the
optimization and add UV. Running the same request again can reuse all completed
stages. Reference spectrum comparisons are recalculated from the raw properties,
so changing only a reference spectrum does not require a new quantum calculation.
To request fresh calculations, set `job.reuse = False` before running; previous
files are still retained in history.

The returned `reuse_report` lists each stage's `job`, `action` (`reused`,
`executed` or `failed`) and reason, with counts for each action. A stage's `key`
identifies its cache entry and is a SHA-256 hash of its normalized inputs and
dependencies. The report is also
saved to `.qcforever/reuse-report.json`; `.qcforever/index.json` records the input
and output provenance. The cache holds wavefunctions and restart files needed for
later stages. With Gaussian `restart=False`, cleanup removes checkpoints from
the cache and history, so stages requiring those checkpoints are recomputed on
the next run. Logs, records and histories still consume additional disk space.

Results produced before this feature lack this execution record. An old ORCA or
GAMESS result is imported only if its complete executed input can be recovered
from the log and matches the requested stage, normal completion is verified,
and required outputs remain available. Gaussian legacy jobs and logs without
enough evidence are recomputed once to establish provenance. A saved input file
alone, a successful launcher exit code, or a matching molecule filename is not
sufficient. Records are specific to a backend; Gaussian results are never used
as ORCA/GAMESS numerical results.

Independent jobs must not write to the same molecule directory concurrently.
`.qcforever/running.lock` records the active process. After an unclean process
termination, check that process before removing a stale lock and retrying.

### Options

By specifying the molecular properties you want as an option variable string,
QCforever automatically calculates them.
The variable string must consisted of the following options,
seperated with more than one space:

Following options are currently available:

| Option name | Description | Gaussian | Gamess | ORCA |
|---|---|---|---|---|
|symm| just specify symmetry of a molecule.|:white_check_mark:|||
|volume| Compute the volume (in cm**3/mol) of a molecule.|:white_check_mark:|||
|opt| just perform geometry optimization of a molecule.|:white_check_mark:|:white_check_mark:|:white_check_mark:|
|nmr| NMR chemical shift (ppm to TMS) of each atom is computed. `nmr=/absolute/path/to/nmr.dat` also compares it with a two-column `position intensity` reference spectrum.|:white_check_mark:||:white_check_mark:|
|uv| absorption wavelengths (nm) (for spin allowed states) are computed. If you add an absolute path of text file that specify peak positions and their peak like (uv=/ab/path2uv.dat), it is possible to the similarity and dissimilarity between the target and reference. |:white_check_mark:|:white_check_mark:|:white_check_mark:|
|energy| SCF energy (in Eh) is printed.|:white_check_mark:|:white_check_mark:|:white_check_mark:|
|cden| charge and spin densities on each atom are computed.|:white_check_mark:|:white_check_mark:|:white_check_mark:|
|homolumo| HOMO/LUMO gap (Eh) calculation|:white_check_mark:|:white_check_mark:|:white_check_mark:|
|dipole| Dipole moment and quadrupole properties at the center of mass; see the output definitions below.|:white_check_mark:|:white_check_mark:|:white_check_mark:|
|polar| Dipole polarizability (in 10**-24 cm**3) of a molecule|:white_check_mark:||:white_check_mark:|
|deen| decomposition energy (in eV) of a molecule|:white_check_mark:|||
|stable2o2| stability to oxygen molecule|:white_check_mark:|||
|vip| vertical ionization potential energy (in eV)|:white_check_mark:|:white_check_mark:|:white_check_mark:|
|vea| vertical electronic affinity (in eV)|:white_check_mark:|:white_check_mark:|:white_check_mark:|
|aip| adiabatic ionization energy (in eV)|:white_check_mark:|:white_check_mark:|:white_check_mark:|
|aea| adiabatic electronic affinity (in eV)|:white_check_mark:|:white_check_mark:|:white_check_mark:|
|fluor| Fluorescence wavelengths (nm). Select a state with `fluor=N` or `fluor=(state=N)`. Gaussian also accepts independent force-constant controls, e.g. `fluor=(FC=CalcFC, state=2)`; see below.|:white_check_mark:|:white_check_mark:|:white_check_mark:|
|nac| Max and rms values of non adiabatic coupling vector (use with fluor option) |:white_check_mark:|||
|tadf| Compute the energy gap (in Eh) between minima in the spin-allowed and spin-forbidden states. Gaussian accepts `tadf=(FC=CalcFC, state=1)` independently of `fluor`.|:white_check_mark:|||
|freq| Compute the variable related to the vibrational analysis of a molecule. `freq=IR.dat,Raman.dat` compares the calculated IR and Raman spectra with two-column `position intensity` reference files. For Gaussian and ORCA, an optional third NMR file (`freq=IR.dat,Raman.dat,NMR.dat`) also enables NMR and compares all three spectra. Before comparison, peaks within 1.0 cm-1 (IR/Raman) or 0.01 ppm (NMR) are merged; intensities are summed and positions are averaged.|:white_check_mark:|:white_check_mark:|:white_check_mark:|
|pka| Compute the energy gap (in Eh) between deprotonated (A-) and protonated (AH) species. The hydrogen atom whose Mulliken charge is the biggest in the system is selected as a protic hydrogen.|:white_check_mark:|||
|stable| try to find a stable structure when the negative frequency is detected.|:white_check_mark:|||
|optspin| try to find a suitable spin multiplicity.|:white_check_mark:|||
|optconf| Search for a stable conformer (default: PM6 for Gaussian/GAMESS, xTB for ORCA). See [Conformer search](#conformer-search-optconf) below.|:white_check_mark:|:white_check_mark:|:white_check_mark:|

### Ground- and excited-state optimization options

Gaussian's `opt=...` controls **only the ground-state optimization**. It is not
inherited by `fluor`, `tadf`, or deprotonated-species optimizations. For example,
`"opt=calcall energy fluor=2"` uses `Opt=(CalcAll,MaxCycles=100)` for the ground
state and the default optimization controls for the second excited state.

Use a group to specify the excited-state root and force-constant strategy:

```python
option = "opt=calcall energy fluor=(FC=CalcFC, state=2)"

# Independent S1 and T1 controls for Gaussian TADF calculations:
option = "opt=calcall energy fluor=(FC=CalcFC, state=1) tadf=(FC=CalcAll, state=1)"
```

`fluor` applies its `FC` setting only to the spin-allowed excited-state
optimization. `tadf` applies its setting only to the spin-forbidden optimization
(normally a triplet for singlet input). Bare `tadf`, or `tadf=(...)` without an
explicit `fluor`, also requests an S1 optimization with default force-constant
controls. Use an explicit `fluor=(...)` to control that S1 separately. Explicit
`fluor` and `tadf` can occur in either order; they are scheduled in S/T order.

Supported fields are `state` and `FC`; names and FC keywords are case-insensitive,
and spaces or field order inside the group do not matter. `state` defaults to 1.
`fluor`, `fluor=2`, and `fluor=(state=2)` remain valid. Gaussian's legacy
`fluor=-1` or `fluor=(state=-1)` selects the brightest computed vertical state
and automatically requests the UV calculation. This can be combined with FC,
for example `fluor=(FC=CalcFC,state=-1)`. ORCA/GAMESS and `tadf` require a positive root.
`FC` defaults to no explicit Hessian keyword. Gaussian accepts `CalcFC`,
`CalcAll`, `ReadFC`, `RCFC`, `RecalcFC=N` (positive integer), `default`, and `none`.
For example, `fluor=(FC=RecalcFC=5, state=1)` passes `RecalcFC=5` to its `Opt`
keyword. `ReadFC`/`RCFC` require suitable force constants in the input checkpoint;
specifying them does not create missing Hessian data. Availability of Hessian
calculations still depends on Gaussian's method/state support.

GAMESS and ORCA accept `FC=CalcFC`, `CalcAll`, `ReadFC`, and `RecalcFC=N`
for `fluor`, and the same controls for the ground state via `opt=CalcFC`,
`opt=CalcAll`, `opt=ReadFC`, or `opt=RecalcFC=N`. The controls are independent:

```python
option = "opt=CalcFC energy fluor=(FC=CalcAll, state=1)"
# ORCA convergence and Hessian controls can be combined:
option = "opt=(tight, FC=CalcFC) fluor=(FC=RecalcFC=5, state=2)"
```

| FC | Meaning | ORCA `%geom` | GAMESS `$STATPT` |
|---|---|---|---|
| `CalcFC` | Compute initial Hessian once; update it afterwards | `Calc_Hess true` | `HESS=CALC IHREP=0` |
| `CalcAll` | Recompute at every optimization step | `Calc_Hess true Recalc_Hess 1` | `HESS=CALC IHREP=1` |
| `RecalcFC=N` | Compute initially and every N steps | `Calc_Hess true Recalc_Hess N` | `HESS=CALC IHREP=N` |
| `ReadFC` | Start with a previously computed Hessian | `InHess Read InHessName "...hess"` | `HESS=READ IHREP=0` and a `$HESS` block |

See the native controls in the [ORCA geometry optimization manual](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html)
and the [GAMESS input manual (`$STATPT` and `$FORCE`)](https://www.msg.chem.iastate.edu/GAMESS/GAMESS_Manual/input.pdf).

These ORCA/GAMESS calculations use numerical differentiation of gradients
(`NumHess true` / `$FORCE METHOD=SEMINUM NVIB=2`). This also applies to TDDFT
excited-state optimizations; GAMESS conventional TDDFT requires singlet
excited-state gradients. Numerical Hessians can substantially increase cost,
especially with `CalcAll`. `RCFC` remains Gaussian-only.

For ORCA/GAMESS `ReadFC`, first complete the corresponding stage with `CalcFC`,
`CalcAll`, or `RecalcFC=N`, keeping the calculation name and `.qcforever`
directory. QCforever selects a Hessian from a verified completed calculation,
checking file hashes, executable provenance, atom order, method, basis, charge,
multiplicity, solvent and the target excited state. Nearby coordinates may differ:
the Hessian supplies an initial curvature estimate and does not restore the
optimizer trajectory. Ground-state Hessians are not substituted for excited-state
ones. Missing, modified or incompatible data causes an explicit error; there is
no fallback to a guessed Hessian. Importing arbitrary native Hessian files or
resuming an interrupted optimization trajectory is not provided by this option.

TADF remains Gaussian-only. Unknown fields, duplicate settings and invalid roots
are rejected before launching a calculation. Changing only an excited-state FC
setting leaves the ground-state input unchanged, allowing its existing result
to be reused when the remaining reuse checks pass.

### Conformer search (`optconf`)

QCforever provides conformer search to prepare promising structures for subsequent
DFT calculations. It first generates multiple candidate conformations, then
optimizes them with a semiempirical method. LAQA-style adaptive relaxation
(successive halving, SH, by default) focuses optimization effort on promising
candidates to reduce the computational cost.

Simply add `optconf` to the option string, for example `"optconf energy"`.
The default search uses ETKDGv3 + MMFF94s followed by PM6, aiming to find a
promising starting structure at a modest cost rather than exhaustively optimize
every candidate.

To change generation or relaxation, group the settings inside `optconf=(...)`:

```python
option = "optconf=(conf_gen=medium, conf_relax=xtb, laqa=50) energy"
```

| `conf_gen` | Conformer generation (fallback order) |
|---|---|
| `low` (default) | ETKDGv3 + MMFF94s |
| `medium` | Torsional Diffusion → ETKDGv3 + MMFF94s |
| `high` | DiTMC → Torsional Diffusion → ETKDGv3 + MMFF94s |
| `legacy` | The complete old FAFOOM generation **and relaxation** workflow |

Medium/high require the model installation described under Install.
Set `conf_relax=xtb` for GFN2-xTB or `conf_relax=pm6` for PM6 (default).
Install the corresponding executable separately; PM6 requires Gaussian 16 even
when the subsequent property calculation uses GAMESS.

For low/medium/high, the default (`laqa=50`) stops when 50% of the initial relaxation candidates have
converged. To spend more effort on relaxation, increase this target:

| Option | Behavior |
|---|---|
| Omitted, `laqa`, or `laqa=50` | Stop at 50% convergence (default) |
| `laqa=20` | Stop at 20% convergence |
| `laqa=80` | Stop at 80% convergence |
| `laqa=100` | Target all candidates; failures/safety limits may prevent completion |
| `laqa=off` | Relax all candidates continuously, with the same native optimizer and convergence criteria |

These percentages refer to converged candidates, not a percentage of compute budget.
Parallel relaxation uses one core per candidate, up to `nproc` concurrent
candidates; xTB keeps paused processes in memory.
PM6 uses RFO with tight optimization and SCF convergence for both adaptive and
continuous relaxation.

For the unchanged old workflow, use
`"optconf=(conf_gen=legacy, conf_relax=xtb) energy"`.
Legacy retains its own candidate count, optimizer and stopping rules;
`laqa` and `conformer_config` overrides cannot be combined with it.
The older syntax (`optconf=xtb optconf_medium laqa=50`) remains supported.

The selected structure is saved as `optimized_structures.sdf`; summaries and
detailed logs are under `conformer_search/` in the calculation directory.
After low/medium/high finishes, generated/final structures and relaxation histories
(`trajectory.jsonl` and `trajectory.xyz`) are retained; temporary generation and
xTB/PM6 files are removed after successful export. Failed exports preserve the
original files. Legacy and subsequent DFT checkpoint handling are unchanged.
In `electronic/`, `all_converged.sdf`, `audit.json` and `settings.json` remain
directly accessible; candidate histories and other execution records are bundled
in `details.tar.gz`. Original records are removed after successful compression.
Failed exports or archives leave the unarchived records available.
Failed calculations retain a short `failure.log`. To keep complete xTB/PM6 text
logs for every candidate, add `keep_logs=true`, for example
`optconf=(conf_gen=high, conf_relax=pm6, keep_logs=true)`.
Logs are collected in each candidate's `native.log` inside the archive; checkpoints and scratch are
still removed. The default is `keep_logs=false` (low/medium/high only).
In low/medium/high, if only a structure with geometry/stereo warnings is available, it is used with
a console warning; the warning is also recorded in the returned summary.
For advanced settings, set `job.conformer_config = "conformer.yaml"`;
see [defaults.yaml](qcforever/conformer_search/defaults.yaml).

### Quadrupole output with `dipole`

The `dipole` option also extracts quadrupole properties. The existing `dipole`
return value is preserved: Gaussian returns `[X, Y, Z, total]` in Debye in its
original output frame; GAMESS returns the total dipole in Debye.

Each backend's `parse_log.py` extracts the tensor and its reference data:
Gaussian uses `Extract_QuadrupoleMoment`, and GAMESS uses `getQuadrupoleMoment`.
The RunPack passes those numerical values to `util/quadrupole.py`, which performs
origin/convention transformations and calculates the isotropic and anisotropic
scalars. The utility does not read logs or checkpoint files.

The additional output keys are:

| Key | Meaning |
|---|---|
| `quadrupole` | Primitive Cartesian charge second moment `Q`, as a dictionary with `unit`, `XX`, `YY`, `ZZ`, `XY`, `XZ`, and `YZ`. |
| `quadrupole_traceless` | `T = Q - trace(Q) I / 3`, in the same dictionary format. |
| `quadrupole_iso` | `trace(Q) / 3`, including its sign. |
| `quadrupole_anisotropy` | Nonnegative scalar `sqrt(3/2 * sum_ij(T_ij**2))`; off-diagonal entries are counted twice. |
| `quadrupole_info` | Units, tensor convention, origin, source program, availability status, and the reason for any unavailable data. |

All quadrupole components and scalars use **Debye Angstrom**. The two scalars
are invariant under rotations about a fixed origin. For comparison across
calculations, the new quadrupole results use the **center of mass** as their
origin; this does not change the legacy `dipole` output.

For Gaussian, the primitive tensor and dipole vector come from the last
population-analysis block of the requested calculation/state. The origin is
shifted to the center of mass using that block's preceding orientation, total
charge, and atomic masses (including specified isotopes). Masses are read from
`AtmWgt` in the log, or from `Real atomic weights` in the matching `.fchk` file.
If the tensor or origin-transformation data is missing, the new properties are
`None` and `quadrupole_info` explains why; the dipole remains available.

For GAMESS, `dipole` uses a separate property single point at the ground-state
geometry so adding it does not rerun optimization or unrelated TDDFT stages.
Its input is regenerated from the ground-state `.dat` file and includes:

```text
 $ELMOM IEMOM=2 WHERE=COMASS $END
```

GAMESS prints a traceless quadrupole in the Buckingham convention, `3 T / 2`.
QCforever multiplies its six Cartesian components by **2/3** to match the
Gaussian traceless convention before calculating the anisotropy. The primitive
trace cannot be recovered from this output: **`quadrupole` and
`quadrupole_iso` are `None` for GAMESS**, not zero. Available traceless data has
status `traceless_only`. If the final electrostatic block lacks a complete
quadrupole, all new properties are `None` with status `unavailable`; an earlier
step's tensor is not reused. PCM calculations use the final `SYSTEM
ELECTROSTATIC MOMENTS` block, consistently with the existing dipole extraction.

For example:

```python
option = "opt energy dipole"
# ... run the calculation ...
print(output["quadrupole_traceless"]["XX"])
print(output["quadrupole_anisotropy"])
print(output["quadrupole_info"])
```

Definitions and program conventions:
[GAMESS $ELMOM](https://www.msg.chem.iastate.edu/GAMESS/GAMESS_Manual/input.pdf),
[Gaussian output format](https://docs.iochem-bd.org/data/codes/gaussian/multipole-d3e16989.html),
[NIST quadrupole conventions](https://cccbdb.nist.gov/quadrupolex.asp).

## License

This package is distributed under the MIT License.

## Contact

- Masato Sumita (masato.sumita@riken.jp)
