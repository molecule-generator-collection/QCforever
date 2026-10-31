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

1. Quantum chemistry software: [Gaussian](https://gaussian.com/)==16 or [GAMESS (sockets)](https://www.msg.chem.iastate.edu/gamess/)==30 SEP 2022 (R2)
2. [Python](https://www.anaconda.com/download/)==3.11
3. [RDKit](https://pypi.org/project/rdkit/)>=2023.3.3
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

Example codes for QCforever (gaussian_main.py and gamess_main.py) are prepared.

For Gaussian
```bash
python gaussian_main.py input_file
```
For GAMESS
```bash
python gamess_main.py input_file
```

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
|dipole| Dipole moment and quadrupole properties at the center of mass; see the output definitions below.|:white_check_mark:|:white_check_mark:|
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
|optconf| Search for a stable conformer (default: PM6). See [Conformer search](#conformer-search-optconf) below.|:white_check_mark:|:white_check_mark:|

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

For GAMESS, `dipole` adds the following group to generated inputs, including
inputs regenerated from `.dat` files:

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
