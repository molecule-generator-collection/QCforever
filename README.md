# QCforever

QCforever is a wrapper of Gaussian (https://gaussian.com). To compute obsevable properties of a molecule through quantum chemical computation (QC),
Multi step computation is demanded. QCforever automates this process and calculates multiple physical properties of molecules simultaneously. 

## Requirements
1. [Gaussian](https://gaussian.com)==16
2. [Python](https://www.anaconda.com/download/)>=3.7 
3. [RDKit](https://anaconda.org/rdkit/rdkit)

## How to use

### Install
Just copy whole GaussianRunPack directory to your working directory. 

### Usage
Example code for QCforever (main.py) is prepared.
QCforever accept a sdf, xyz , Gaussian chk, or a Gaussian fchk file.

> python main.py input_file

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
	adiabtic electronic affinity (in eV)

fluor:
	wavelength (in nm) of fluorescence are computed. 
	if you want to specify the state that emits fluorescence, you can specify the index of state like
	“fluore=#” (# is an integer, default is “fluor=1”)

tadf:
	Compute the energy gap (in Eh) between minimum in the spin allowed state 
	and the spin forbidden state.

freq: 
	Compute the variable related to the vibrational analysis of a molecule. IR, Raman, etc
	
	
pka:
	Compute the energy gap (in Eh) between deprotonated (A-) and protonated (AH) species.
	The hydrogen atom whose Mulliken charge is the biggest in the system is selected as a protic hydorogen. 

## License

## Contact
- Masato Sumita
