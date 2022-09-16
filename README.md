# QCforever


![Robot2+PC](https://user-images.githubusercontent.com/46772738/188896764-65ab12c1-3cc9-421d-8d87-ed33c932380a.png)

QCforever (https://doi.org/10.1021/acs.jcim.2c00812) is a wrapper of Gaussian (https://gaussian.com). To compute obsevable properties of a molecule through quantum chemical computation (QC),
Multi step computation is demanded. QCforever automates this process and calculates multiple physical properties of molecules simultaneously. 

## Requirements
1. [Gaussian](https://gaussian.com)==16
2. [Python](https://www.anaconda.com/download/)>=3.7 
3. [RDKit](https://anaconda.org/rdkit/rdkit)

## How to use

### Install
Just copy  GaussianRunPack directory to your working directory. 

### Example
Example code for QCforever (main.py) is prepared.

> python main.py input_file

### Usage
In the working directory where you put the GaussianRunPack directory, import GaussianRunPack to use this package.

> import GaussianRunPack

Then make an instance (here is test).

> test = GaussianRunPack.GaussianDFTRun(Functional, basis_set, ncore, option, input_file, solvent='water', restart=False)

"Functional" is to specify the functional in the density functional theory (DFT). Currently supported functionals are the following:
  "BLYP", "B3LYP", "X3LYP", "LC-BLYP", "CAM-B3LYP".

"basis_set" is to specify the basis set. Current QCforever supports the following basis set:
  "LANL2DZ", "STO-3G", "3-21G", "6-31G", "6-311G", "3-21G*", "3-21+G*", "6-31G*", "6-311G**", "6-31G**", "6-31+G*", "3-21G", "3-21G", "6-31+G**", "6-311+G*", "6-311+G**"

"ncore" is an integer to specify the number of core for QC with Gaussian.

"option" is a string for specifying molecular properties as explained later.

"input_file" is a string to specify the input file. QCforever accepts a sdf, xyz , Gaussian chk, or a Gaussian fchk file.

"solvent" is to include the solvent effect through PCM. The default value is "0", in vacuo.

"restart" is to control to save molecular information as fchk or xyz. 
The Default value is True that means molecular information is saved as a Gaussian fchk file (electronic structure is also saved.), 
otherwise molecular information is saved as a xyz files (electronic structure is not saved).

And excuse Gaussian as the followign.

> outdic = test.run_gaussian()

In this examples, obtained results are saved as python dictionary style.

### Options 
By specifying the molecular properties you want as an option variable string,
QCforever automatically calculates them. 
The variable string must consisted of the following options, 
seperated with more than one space:

Following options are currently available:

symm:
	just specify symmetry of a molecule.

opt:
	just perform geometry optimization of a molecule.

nmr:
	NMR chemical shift (ppm to TMS) of each atom is computed.
  
uv:
	absorption wavelengths (nm)  (for spin allowed states) are computed. If you add an absolute path of text file that specify peak positions and their peak like (uv=/ab/path2uv.dat), it is possible to the similarity and dissimilarity between the target and reference.

energy: 
	SCF energy (in Eh) is printed.

cden:
	charge and spin densities on each atom are computed.

homolumo:
	HOMO/LUMO gap (Eh) calculation

dipole:
	dipole moment of a molecule

deen:
	decomposition energy (in eV) of a molecule

stable2o2:
	stability to oxygen molecule

vip:
	vertical ionization potential energy (in eV)
	
vea:
	vertiacl electronic affinity (in eV)
	
aip:
	adiabatic ionization energy (in eV) 
	
aea:
	adiabatic electronic affinity (in eV)

fluor:
	wavelength (in nm) of fluorescence are computed. 
	if you want to specify the state that emits fluorescence, you can specify the index of state like
	“fluor=#” (# is an integer, default is “fluor=1”)

tadf:
	Compute the energy gap (in Eh) between minimum in the spin allowed state 
	and the spin forbidden state.

freq: 
	Compute the variable related to the vibrational analysis of a molecule. IR, Raman, etc
	
	
pka:
	Compute the energy gap (in Eh) between deprotonated (A-) and protonated (AH) species.
	The hydrogen atom whose Mulliken charge is the biggest in the system is selected as a protic hydorogen. 

## License
This package is distributed under the MIT License.

## Contact
- Masato Sumita (masato.sumita@riken.jp)
