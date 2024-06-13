#    Copyright 2015 Adriana Supady
#
#    This file is part of fafoom.
#
#   Fafoom is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   Fafoom is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#   along with fafoom.  If not, see <http://www.gnu.org/licenses/>.

'''Wrapper for RDKit force-field routines'''
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalForceFields

kcalmol2eV = 0.0433641
kcalmol2hartree = 0.00159362 

class FFObject():
    """Create and handle force-field objects."""
    def __init__(self, sdf_string, force_field='uff', optsteps=1000,
                 force_tol=1.0e-4, energy_tol=1.0e-6,
                 sdf_out='optimized_structures.sdf'):

        """Initialize the FFObject.
        Args(required):
            sdf_string (str): sdf file string
        Args(optional):
        (for the minimization method used for optimization)
            force_field (str): name of the force_field to use (default='uff')
            optsteps    (default=1000)
            force_tol   (default=1e-4)
            energy_tol  (default=1e-6)
        Raises:
            ValueError: if the force field is not 'uff' or 'mmff94'
        """
        self.sdf_string = sdf_string
        self.force_field = force_field
        if self.force_field not in ['uff', 'mmff94']:
            raise ValueError("Unknown force field.")
        self.optsteps = optsteps
        self.force_tol = force_tol
        self.energy_tol = energy_tol
        self.sdf_out = sdf_out

    def run_ff(self):
        """Perform the force field minimization.

        Args:
            sdf_string (str)
        """
        mol = Chem.MolFromMolBlock(self.sdf_string, removeHs=False)
        if self.force_field == 'mmff94':
            molprop = ChemicalForceFields.MMFFGetMoleculeProperties(mol)
            ff = ChemicalForceFields.MMFFGetMoleculeForceField(mol, molprop)
        elif self.force_field == 'uff':
            ff = AllChem.UFFGetMoleculeForceField(mol)
        ff.Minimize(int(self.optsteps), float(self.force_tol),
                    float(self.energy_tol))
        self.sdf_string_opt = Chem.MolToMolBlock(mol)
        self.energy = float('{0:.4f}'.format(ff.CalcEnergy()))

    def get_energy(self, unit='eV'):
        """Get the energy of the molecule.

        Returns:
            energy (float)
        Raises:
            AttributeError: if energy hasn't been calculated yet
        """
        if not hasattr(self, 'energy'):
            raise AttributeError("The calculation wasn't performed yet.")
        else:
            if unit == 'hartree':
                return kcalmol2hartree * self.energy
            elif unit == 'kcal':
                return self.energy
            else: # unit == 'eV':
                return kcalmol2eV * self.energy

    def get_sdf_string_opt(self):
        """Get the optimized sdf string.

        Returns:
            optimized sdf string (str)
        Raises:
            AttributeError: if the optimization hasn't been performed yet
        """
        if not hasattr(self, 'sdf_string_opt'):
            raise AttributeError("The calculation wasn't performed yet.")
        else:
            return self.sdf_string_opt

    def save_to_file(self): 
        if os.path.isfile(self.sdf_out):
            f = open(self.sdf_out, "a")
        else:
            f = open(self.sdf_out, "w")
        s = str(self.sdf_string_opt) + '\n' \
            + ">  <Energy>"+'\n' \
            + str(self.energy)+'\n\n' \
            + "$$$$"+'\n'
        f.write(s)
        f.close()


def ff_exec(sdf_string, force_field='uff', optsteps=1000,
            force_tol=1.0e-4, energy_tol=1.0e-6,
            sdf_out='optimized_structures.sdf'):

    ff_object = FFObject(sdf_string, force_field, optsteps,
                         force_tol, energy_tol, sdf_out)
    ff_object.run_ff()
    energy = ff_object.get_energy()
    sdf_string_opt = ff_object.get_sdf_string_opt()

    return energy, sdf_string_opt


if __name__ == '__main__':

    sdf_string = '\nH2O\n\n' + \
    '  3  2  0  0  0  0  0  0  0  0999 V2000\n' + \
    '    0.0000    0.0000   -0.3894 O   0  0  0  0  0  0  0  0  0  0  0  0\n' + \
    '    0.7630    0.0000    0.1947 H   0  0  0  0  0  0  0  0  0  0  0  0\n' + \
    '   -0.7630    0.0000    0.1947 H   0  0  0  0  0  0  0  0  0  0  0  0\n' + \
    '  1  2  1  0  0  0  0\n' + \
    '  1  3  1  0  0  0  0\n' + \
    'M  END\n$$$$'
    force_field = 'mmff94'

    print('initial coord')
    print(sdf_string)

    energy, sdf_string_opt = ff_exec(sdf_string, force_field, optsteps=1000,
                 force_tol=1.0e-4, energy_tol=1.0e-6,
                 sdf_out='optimized_structures.sdf')

    print('energy:', energy)
    print('opt coord')
    print(sdf_string_opt)
