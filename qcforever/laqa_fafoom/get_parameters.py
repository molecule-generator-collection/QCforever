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
''' Communicate between the structure and the degrees of freedom.'''
from __future__ import division
from rdkit import Chem
from rdkit.Chem import AllChem

from qcforever import laqa_fafoom


def get_atoms_and_bonds(smiles):
    """Build the molecule from SMILES and return the number of atoms and bonds.

    Args(required):
        smiles (str): one-line representation of the molecule
    Returns:
        Number of atoms, number of bonds
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("The smiles is invalid")
    mol = Chem.AddHs(mol)
    return mol.GetNumAtoms(), mol.GetNumBonds()


def get_positions(type_of_deg, smiles, **kwargs):
    """Find the positions (tuples of atom indicies) of the degrees of freedom.

    Args(required):
        type_of_deg (str)
        smiles (str)
        if cistrans should be optimized:
            smarts_cistrans
    Args(optimal):
        list_of_torsion (list)
        smarts_torsion (str)
        filter_smarts_torsion (str)
        list_of_cistrans (list)
        list_of_pyranosering (list)
    Returns:
        list of touples defining the positions of the degree of freedom
    """

    if type_of_deg == "torsion":
        if 'list_of_torsion' in kwargs:
            return laqa_fafoom.deg_of_freedom.Torsion.find(smiles, positions=kwargs['list_of_torsion'])
        else:
            if 'smarts_torsion' in kwargs:
                if 'filter_smarts_torsion' in kwargs:
                    return laqa_fafoom.deg_of_freedom.Torsion.find(smiles,
                                        smarts_torsion=kwargs['smarts_torsion'],
                                        filter_smarts_torsion=
                                        kwargs['filter_smarts_torsion'])
                else:
                    return laqa_fafoom.deg_of_freedom.Torsion.find(smiles,
                                        smarts_torsion=kwargs['smarts_torsion'])
            else:
                return laqa_fafoom.deg_of_freedom.Torsion.find(smiles)
    if type_of_deg == "cistrans":
        if 'list_of_cistrans' in kwargs:
            return laqa_fafoom.deg_of_freedom.CisTrans.find(smiles, positions=kwargs['list_of_cistrans'])
        else:
            return laqa_fafoom.deg_of_freedom.CisTrans.find(smiles,
                                 smarts_cistrans=kwargs['smarts_cistrans'])

    if type_of_deg == "pyranosering":
        if 'list_of_pyranosering' in kwargs:
            return laqa_fafoom.deg_of_freedom.PyranoseRing.find(smiles,
                                     positions=kwargs['list_of_pyranosering'])
        else:
            return laqa_fafoom.deg_of_freedom.PyranoseRing.find(smiles)


def create_dof_object(type_of_deg, positions):
    """Initialize the degree of freedom from the positions

    Args:
        type_of_deg (str)
        positsion (list)
    Returns:
        degree of freedom object
    """
    if type_of_deg == "torsion":
        return laqa_fafoom.deg_of_freedom.Torsion(positions)
    if type_of_deg == "cistrans":
        return laqa_fafoom.deg_of_freedom.CisTrans(positions)
    if type_of_deg == "pyranosering":
        return laqa_fafoom.deg_of_freedom.PyranoseRing(positions)


def template_sdf(smiles, distance_cutoff_1, distance_cutoff_2):
    """Create a template sdf string and writes it to file.

    Args(required):
        smiles (str): one-line representation of the molecule
    Args(optional):
        distance_cutoff_1 (float): min distance between non-bonded atoms [A]
        distance_cutoff_2 (float): max distance between bonded atoms [A]
    Returns:
        sdf string
    """
    cnt = 0
    sdf_check = True
    while sdf_check:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol,useRandomCoords=True)
        AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
        Chem.SDWriter('mol.sdf').write(mol)
        sdf_string = Chem.MolToMolBlock(mol)
        check = laqa_fafoom.utilities.check_geo_sdf(sdf_string, distance_cutoff_1, distance_cutoff_2)
        if check:
            sdf_check = False
            Chem.SDWriter('mol.sdf').write(mol)
        else:
            cnt += 1
    return sdf_string
