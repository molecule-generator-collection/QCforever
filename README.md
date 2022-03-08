# GaussianRunPack
GaussianRunPack is a wrapper of Gaussian (https://gaussian.com). It is possible to evaluate the molecule with its multiple properties.

## Requirements
1. [Gaussian](https://gaussian.com)==16
2. [Python](https://www.anaconda.com/download/)>=3.7 
3. [RDKit](https://anaconda.org/rdkit/rdkit)

## How to use
GaussianRunPack accept a sdf file or a gaussian chk file.

> python main.py input_file

By specifying the molecular properties you want as an option variables,
the program will automatically calculate them. 

Following options are currently available:

opt:
	just perform geometry optimization of the molecule.

nmr:
	NMR chemical shift (ppm to TMS) of each atom is computed.
  
uv:
	absorption wavelengths (nm)  (for spin allowed states) are computed.

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

ipe:
	ionization potential energy
eae:
	electronic affinity 
npe:
	neutralization energy from a positively charged molecule (geometry optimized)
nne:
	neutralization energy from a negatively charged molecule (geometry optimized)

fluor:
	wavelength (in nm) of fluorescence are computed. 
	if you want to specify the state that emits fluorescence, you can specify the index of state like
	“fluore=#” (# is an integer, default is “fluor=1”)

tadf:
	compute the energy gap (in Eh) between minimum in the spin allowed state 
	and the spin forbidden state.

## License

## Contact
- Masato Sumita
