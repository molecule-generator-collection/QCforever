#    Copyright 2015 Adriana Supady & Mateusz Marianski
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
"""Measure and set dihedral angles and rings."""
from __future__ import division
from operator import itemgetter
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from qcforever import laqa_fafoom


def ig(x):
    return itemgetter(x)


def dihedral_measure(sdf_string, position):
    """ Measure the dihedral angle.

    Args:
        sdf_string (string)
        position (list): 4 atoms defining the dihedral
    Returns:
        float value
    Raises:
        ValueError: If the lenght of the list is not equal 4.
    """
    if len(position) != 4:
        raise ValueError("The position needs to be defined by 4 integers")
    mol = Chem.MolFromMolBlock(sdf_string, removeHs=False)
    val = float(rdMolTransforms.GetDihedralDeg(
                mol.GetConformer(),
                ig(0)(position), ig(1)(position),
                ig(2)(position), ig(3)(position)))
    return float('{0:.2f}'.format(val))


def dihedral_set(sdf_string, position, value):
    """ Set the dihedral angle.

    Args:
        sdf_string (string):
        position (list): 4 atoms defining the dihedral
        value : value to set
    Returns:
        modified sdf_string
    Raises:
        ValueError: If the lenght of the list is not equal 4.
    """
    if len(position) != 4:
        raise ValueError("The position needs to be defined by 4 integers")
    mol = Chem.MolFromMolBlock(sdf_string, removeHs=False)
    rdMolTransforms.SetDihedralDeg(mol.GetConformer(), ig(0)(position),
                                   ig(1)(position), ig(2)(position),
                                   ig(3)(position), value)

    return Chem.MolToMolBlock(mol)


def pyranosering_set(sdf_string, position, new_dih, new_ang):
    """ Set the pyranosering.

    Args:
        sdf_string (string)
        position (list): 7 atoms defining the ring, i.e. positions of
                        ['C0','C1','C2','C3','C4','O', 'O0']
        new_dih (list) : 5 values for the dihedral angles
        new_ang (list): 5 values for the bond angles
    Returns:
        modified sdf_string
    Raises:
        ValueError: If the lenght of the position is not equal 7 ot if the
        length of new_dih/new_ang is not equal to 5.
    """
    if len(position) != 7:
        raise ValueError("The position needs to be defined by 7 integers")
    if len(new_dih) != 5:
        raise ValueError("Five dihedral angles are needed for the new ring "
                         "conformation.")
    if len(new_ang) != 5:
        raise ValueError("Five bond angles are needed for the new ring "
                         "conformation.")

    from scipy.linalg import expm3

    atoms_ring = {}
    for n, name in zip(range(len(position)),
                       ['C0', 'C1', 'C2', 'C3', 'C4', 'O', 'O0']):
        atoms_ring[name] = position[n]

    def initialize(sdf_string):
        molecule = Chem.MolFromMolBlock(sdf_string, removeHs=False)
        return molecule

    def calculate_normal_vector(list_of_atoms, xyz):
        """Calculate the normal vector of a plane by
        cross product of two vectors belonging to it.

        Args:
            list_of_atoms: list of 3 atoms
            xyz: numpy array with atoms xyz position
        """

        r0 = xyz[list_of_atoms[1], :] - xyz[list_of_atoms[0], :]
        r1 = xyz[list_of_atoms[2], :] - xyz[list_of_atoms[1], :]
        cross_product = np.cross(r1, r0)

        return cross_product

    def measure_angle(list_of_atoms, xyz):
        """Calculate an angle between three atoms:
        angle = acos(dot(X,Y)/(norm(X)*norm(Y)))

        Args:
            list_of_atoms: list of 3 atoms
            xyz: numpy array with atoms xyz positions
        """

        r0 = xyz[list_of_atoms[0], :] - xyz[list_of_atoms[1], :]
        r1 = xyz[list_of_atoms[2], :] - xyz[list_of_atoms[1], :]

        norm_r0 = np.sqrt(np.sum(r0**2))
        norm_r1 = np.sqrt(np.sum(r1**2))
        norm = norm_r0*norm_r1

        dot_product = np.dot(r0, r1)/norm
        angle = np.arccos(dot_product)

        #Calculate the axis of rotation (axor):
        axor = np.cross(r0, r1)

        return angle*180.0/np.pi, axor

    def measure_dihedral(list_of_atoms, xyz):
        """Calculate a dihedral angle between two planes defined by
        a list of four atoms. It returns the angle and the rotation axis
        required to set a new dihedral.

        Args:
            list_of_atoms: list of 4 atoms
            xyz: numpy array with atom xyz positions
        """

        plane1 = calculate_normal_vector(list_of_atoms[:3], xyz)
        plane2 = calculate_normal_vector(list_of_atoms[1:], xyz)
        #Calculate the axis of rotation (axor)
        axor = np.cross(plane1, plane2)

        #Calculate a norm of normal vectors:
        norm_plane1 = np.sqrt(np.sum(plane1**2))
        norm_plane2 = np.sqrt(np.sum(plane2**2))
        norm = norm_plane1 * norm_plane2
        #Measure the angle between two planes:
        dot_product = np.dot(plane1, plane2)/norm
        alpha = np.arccos(dot_product)

        #The cosine function is symetric thus, to distinguish between
        #negative and positive angles, one has to calculate if the fourth
        #point is above or below the plane defined by first 3 points:
        ppoint = - np.dot(plane1, xyz[list_of_atoms[0], :])
        dpoint = (np.dot(plane1, xyz[list_of_atoms[3], :])+ppoint)/norm_plane1
        if dpoint >= 0:
            return -(alpha*180.0)/np.pi, axor
        else:
            return (alpha*180.0)/np.pi, axor

    def determine_carried_atoms(at1, at2, conn_mat):
        """Find all atoms necessary to be carried over during rotation
        of an atom 2:

        Args:
            at1, at2: two atoms number
        """

        #1. Zero the connections in connectivity matrix
        tmp_conn = np.copy(conn_mat)
        tmp_conn[at1, at2] = 0
        tmp_conn[at2, at1] = 0
        #2. Determine the connected atoms:
        carried_atoms = [at2]
        s = True
        while s:
            s = False
            #Always iterate over entire list because I might have branching
            for at in carried_atoms:
                #List of indexes of connected atoms:
                conn_atoms = np.where(tmp_conn[at] != 0)[0]
                conn_atoms.tolist
                for x in conn_atoms:
                    if x not in carried_atoms:
                        carried_atoms.append(x)
                        s = True
        return carried_atoms

    def set_angle(list_of_atoms, new_ang, atoms_ring, xyz, conn_mat):
        """Set a new angle between three atoms

        Args:
            list_of_atoms: list of three atoms
            new_ang: value of dihedral angle (in degrees) to be set
            atoms_ring: dictionary of atoms in the ring. It recognizes
                        if the last atom is 'C0O' (obsolete)
            xyz: numpy array with atoms xyz positions
            conn_mat: connectivity matrix
        Returns:
            xyz: modified numpy array with new atoms positions
        """
        #Determine the axis of rotation:
        old_ang, axor = measure_angle(list_of_atoms, xyz)
        norm_axor = np.sqrt(np.sum(axor**2))
        normalized_axor = axor/norm_axor

        #Determine which atoms should be dragged along with the bond:
        carried_atoms = determine_carried_atoms(list_of_atoms[1],
                                                list_of_atoms[2], conn_mat)

        #Each carried_atom is rotated by euler-rodrigues formula:
        #Also, I move the midpoint of the bond to the mid atom
        #the rotation step and then move the atom back.

        rot_angle = np.pi*(new_ang - old_ang)/180.
        #Shake it, baby! Rotation matrix:
        #print old_ang, new_ang, rot_angle*180./np.pi
        rot1 = expm3(np.cross(np.eye(3), normalized_axor*rot_angle))
        translation = xyz[list_of_atoms[1], :]
        for at in carried_atoms:
            xyz[at, :] = np.dot(rot1, xyz[at, :]-translation)
            xyz[at, :] = xyz[at, :]+translation
        return xyz

    def set_dihedral(list_of_atoms, new_dih, atoms_ring, xyz, conn_mat):
        """Set a new dihedral angle between two planes defined by
        atoms first and last three atoms of the supplied list.

        Args:
            list_of_atoms: list of four atoms
            new_dih: value of dihedral angle (in degrees) to be set
            atoms_ring: dictionary of atoms in the ring. It recognizes
                       if the last atom is 'C0O'
            xyz: numpy array with atoms xyz positions
            conn_mat: connectivity matrix
        Returns:
            xyz: modified numpy array with new atoms positions
        """

        #Determine the axis of rotation:
        old_dih, axor = measure_dihedral(list_of_atoms, xyz)
        norm_axor = np.sqrt(np.sum(axor**2))
        normalized_axor = axor/norm_axor

        #Check if the bond is the last bond, next to broken one.
        #If yes, refer to the oxygen:
        if 'O0a' in atoms_ring.keys():
            if list_of_atoms[-1] == atoms_ring['O0a']:
                new_dih += 120.0
        else:
            if list_of_atoms[-1] == atoms_ring['O0b']:
                new_dih -= 120.0
        #Determine which atoms should be dragged along with the bond:
        carried_atoms = determine_carried_atoms(list_of_atoms[1],
                                                list_of_atoms[2], conn_mat)
        #Each carried_atom is rotated by Euler-Rodrigues formula:
        #Reverse if the angle is less than zero, so it rotates in
        #right direction.
        #Also, I move the midpoint of the bond to the center for
        #the rotation step and then move the atom back.

        if old_dih >= 0.0:
            rot_angle = np.pi*(new_dih - old_dih)/180.
        else:
            rot_angle = -np.pi*(new_dih - old_dih)/180.
        #Shake it, baby! Rotation matrix:
        rot1 = expm3(np.cross(np.eye(3), normalized_axor*rot_angle))
        translation = (xyz[list_of_atoms[1], :]+xyz[list_of_atoms[2], :])/2
        for at in carried_atoms:
            xyz[at, :] = np.dot(rot1, xyz[at, :]-translation)
            xyz[at, :] = xyz[at, :]+translation

        return xyz

    def mutate_ring(molecule, new_dih, new_ang):
        """Mutate a ring to given conformation defined as a list of torsional
        angles accoring to the 10.1016/S0040-4020(00)01019-X (IUPAC) paper
        """
        n_at = molecule.GetNumAtoms()
        n_bonds = molecule.GetNumBonds()
        m_string = Chem.MolToMolBlock(molecule)

        #Split the string to xyz, connectivity matrix and atom list
        m_coords = m_string.split('\n')[4:4+n_at]
        xyz = np.zeros((n_at, 3))
        atom_list = []
        n = 0
        for line in m_coords:
            xyz[n, :] += np.array(map(float, line.split()[:3]))
            atom_list.append(line.split()[3])
            n += 1
        #Molecule Connectivity Matrix
        m_conn = m_string.split('\n')[4+n_at:4+n_at+n_bonds]
        conn_mat = np.zeros((n_at, n_at))
        for line in m_conn:
            at1 = int(line.split()[0])
            at2 = int(line.split()[1])
            conn_mat[at1-1, at2-1] = 1
            conn_mat[at2-1, at1-1] = 1

        #Introduce a cut between ring C0 and C1:
        #I chose these atoms according to the torsion
        #definitions in the IUPAC paper
        #doi: 10.1016/S0040-4020(00)01019-X
        conn_mat[atoms_ring['C0'], atoms_ring['C1']] = 0
        conn_mat[atoms_ring['C1'], atoms_ring['C0']] = 0

        #Construct a list of atoms in order:
        #C0, C1, C2, C3, C4, O, C0, O0a/b (oxygen at anomeric carbon)
        #I use this list to rotate bonds.
        atoms_list = []
        for x in range(0, 5):
            atoms_list.append(atoms_ring['C'+str(x)])
        atoms_list.append(atoms_ring['O'])
        atoms_list.append(atoms_ring['C0'])
        atoms_list.append(atoms_ring['O0'])

        #Determine the anomer - alpha/beta, based on improper
        #dihedral angle C1-C0-O-O0
        imdih = []
        for at in ['C1', 'C0', 'O', 'O0']:
            imdih.append(atoms_ring[at])
        test_anomer = measure_dihedral(imdih, xyz)[0]
        if test_anomer > 0.0:
            atoms_ring['O0b'] = atoms_ring.pop('O0')
        else:
            atoms_ring['O0a'] = atoms_ring.pop('O0')

        #Adjust the 'internal' angles in the ring:
        for n in range(len(new_ang)):
            xyz = set_angle(atoms_list[n:n+3], new_ang[n], atoms_ring, xyz,
                            conn_mat)
        #Rotate the dihedral angles in the ring:
        for n in range(len(new_dih)):
            xyz = set_dihedral(atoms_list[n:n+4], new_dih[n], atoms_ring, xyz,
                               conn_mat)
        a = []
        a.append("%10s\n" % n_at)
        for n in new_dih:
            a.append("%10.4f" % n)
        a.append("\n")
        for n in range(n_at):
            a.append("%10s%10.4f%10.4f%10.4f\n" % (atom_list[n], xyz[n, 0],
                                                   xyz[n, 1], xyz[n, 2]))
        xyz_string = ''.join(a)
        return xyz_string

    molecule = initialize(sdf_string)
    sdf_string = laqa_fafoom.utilities.xyz2sdf(mutate_ring(molecule, new_dih, new_ang), sdf_string)

    return sdf_string


def pyranosering_measure(sdf_string, position, dict_of_options):
    """Assign the ring to a conformation from the dictionary of options.

    Args:
        sdf_string (string)
        position (list): 7 atoms defining the ring
        dict_of_options (dict) : options for the ring
    Returns:
        An integer that corresponds to the best matching dict key
    Raises:
        ValueError: If the lenght of the position is not equal 7.
    """
    if len(position) != 7:
        raise ValueError("The position needs to be defined by 7 integers")
    ang1 = dihedral_measure(sdf_string, position[0:4])
    ang2 = dihedral_measure(sdf_string, position[1:5])
    ang3 = dihedral_measure(sdf_string, position[2:6])
    ang4 = dihedral_measure(sdf_string, (ig(3)(position), ig(4)(position),
                                         ig(5)(position), ig(0)(position)))
    ang5 = dihedral_measure(sdf_string, (ig(4)(position), ig(5)(position),
                                         ig(0)(position), ig(1)(position)))
    ang6 = dihedral_measure(sdf_string, (ig(5)(position), ig(0)(position),
                                         ig(1)(position), ig(2)(position)))

    all_ang = [ang1, ang2, ang3, ang4, ang5, ang6]

    rmsd_dict = {}

    for key in dict_of_options:
        rmsd_dict[key] = (laqa_fafoom.utilities.tor_rmsd(2, laqa_fafoom.utilities.get_vec(all_ang, dict_of_options[key])))

    return int(min(rmsd_dict.iteritems(), key=ig(1))[0])
