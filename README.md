# QCforever

![Robot2+PC](https://user-images.githubusercontent.com/46772738/188896764-65ab12c1-3cc9-421d-8d87-ed33c932380a.png)

QCforever (https://doi.org/10.1021/acs.jcim.2c00812) is a wrapper of Gaussian (https://gaussian.com) or GAMESS (https://www.msg.chem.iastate.edu/gamess/). 
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
6. [xtb](https://github.com/grimme-lab/xtb/tree/v6.6.1)==6.6.1

## How to use

### Install

```bash
pip install --upgrade git+https://github.com/molecule-generator-collection/QCforever.git
```

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
and [`KTLC-BLYP-BO`](https://doi.org/10.1021/acs.jctc.3c00764) (only for Gaussian).

- ***basis_set*** is to specify the basis set.
  Current QCforever supports the following basis set:
  `LANL2DZ`, `STO-3G`, `3-21G`, `6-31G`, `6-311G`, `3-21G*`, `3-21+G*`, `6-31G*`, `6-311G**`, `6-31G**`,
  `6-31+G*`, `3-21G`, `3-21G`, `6-31+G**`, `6-311+G*`, and `6-311+G**`

- ***ncore*** is an integer to specify the number of core for QC with Gaussian.

- ***option*** is a string for specifying molecular properties as explained later.

- ***input_file*** is a string to specify the input file.
  QCforever accepts a sdf, xyz, Gaussian chk, or a Gaussian fchk file.

- ***solvent*** is to include the solvent effect through PCM.
  The default value is "0", in vacuo.

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
|nmr| NMR chemical shift (ppm to TMS) of each atom is computed.|:white_check_mark:||
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
|tadf| Compute the energy gap (in Eh) between minimum in the spin allowed state and the spin forbidden state.|:white_check_mark:||
|freq| Compute the variable related to the vibrational analysis of a molecule. IR, Raman, etc|:white_check_mark:|:white_check_mark:|
|pka| Compute the energy gap (in Eh) between deprotonated (A-) and protonated (AH) species. The hydrogen atom whose Mulliken charge is the biggest in the system is selected as a protic hydrogen.|:white_check_mark:||
|stable| try to find a stable structure when the negative frequency is detected.|:white_check_mark:||
|optspin| try to find a suitable spin multiplicity.|:white_check_mark:||
|optconf| try to find a stable molecular conformation with xTB (developing).|:white_check_mark:||

## License

This package is distributed under the MIT License.

## Contact

- Masato Sumita (masato.sumita@riken.jp)
