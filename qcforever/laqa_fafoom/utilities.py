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
''' Collection of diverse help/convert functions '''
import os
import numpy as np
import math
import shutil
import configparser

from rdkit import Chem
from rdkit.Chem import AllChem

from operator import itemgetter


# Flow-handling


def backup(filename, obj):
    """ Write the representation of an object (or objects) to a file."""
    with open(filename, 'w') as outf:
        if hasattr(obj, "__len__"):
            for i in range(len(obj)):
                outf.write("%s\n" % repr(obj[i]))
        else:
            outf.write("%s\n" % repr(obj))


def boolean(string):
    """Recover the boolean value from a string and return it."""
    if string in ["False", "false", "FALSE"]:
        return False
    if string in ["True", "true", "TRUE"]:
        return True
    raise ValueError("Cannot be converted to a boolean type")


def number(s):
    """Convert to integer of float if needed"""
    try:
        return int(s)
    except ValueError:
        return float(s)


def print_output(text):
    """Write text to the 'output.txt'. Create it if needed."""
    if os.path.isfile("output.txt"):
        stat = 'a'
    else:
        stat = 'w'
    with open("output.txt", stat) as f:
        f.write(str(text)+'\n')


def remover_file(instance):
    """Remove a file (if it exists)."""
    try:
        os.remove(instance)
    except OSError:
        pass


def remover_dir(instance):
    """Remove a directory (if it exists)."""
    try:
        shutil.rmtree(instance)
    except OSError:
        pass


def file2string(input_file):
    """Read a file to a string and return it."""
    with open(input_file, 'r') as f:
        string = f.read()
    return string


def string2file(string, filename):
    """Write a string to a file"""
    with open(filename, 'w') as target:
        target.write(string)


def set_default(params, dict_default):
    """Set defaults for missing keys and add the key:value pairs to the
    dict."""
    for key in dict_default:
        if key not in params:
            print("Setting a default value for " + str(key) + ": " +
                  str(dict_default[key]))
            params[str(key)] = dict_default[key]
    return params


def file2dict(filename, sections):
    """Parse a file and create a dictionary"""
    config = configparser.ConfigParser()
    config.read(filename, encoding='utf-8')
    new_dict = {}
    for section in sections:
        if config.has_section(section):
            for key, value in config.items(section):
                new_dict[str(key)] = eval(value)
    return new_dict


# Help vector/matrix functions


def ig(x):
    return itemgetter(x)


def cleaner(list_to_clean):
    """ Remove duplicate torsion definion from a list of atom ind. tuples."""
    for_remove = []
    for x in reversed(range(len(list_to_clean))):
        for y in reversed(range(x)):
            ix1, ix2 = ig(1)(list_to_clean[x]), ig(2)(list_to_clean[x])
            iy1, iy2 = ig(1)(list_to_clean[y]), ig(2)(list_to_clean[y])
            if (ix1 == iy1 and ix2 == iy2) or (ix1 == iy2 and ix2 == iy1):
                for_remove.append(y)
    clean_list = [v for i, v in enumerate(list_to_clean)
                  if i not in set(for_remove)]
    return clean_list


def get_vec(vec1, vec2):
    """Calculate difference between vectors of angles [in rad!].
    Args:
        vec1 (list) in deg
        vec2 (list) in deg
    Returns:
        numpy array
    Raises:

        ValueError: if the length of the lists differ

    Warning1: the vectors contain periodic values, i.e -185 -> 175
    Warning2: symmetry is not included here, but can be easily added if the
    index of the 'symmetric' torsion is known
    """
    if len(vec1) != len(vec2):
        raise ValueError("No length match between the lists")
    diff_vec = np.zeros(len(vec1))
    for i in range(0, len(vec1)):
        tor_diff = abs(vec1[i]-vec2[i])
        diff_vec[i] = min(abs(tor_diff), abs(360-tor_diff))/180.0
    return diff_vec


def tor_rmsd(p, vec):
    """Calculate the modified p norm.The difference from standard norm is the
    fact that the addends are divided by the length of the vector."""
    summe = 0
    for i in range(0, len(vec)):
        summe += math.pow(abs(vec[i]), p)
    return math.pow(summe/len(vec), (1.0/p))


def get_cartesian_rms(sdf_string1, sdf_string2):
    """Return the optimal RMS after aligning two structures."""
    ref = Chem.MolFromMolBlock(sdf_string1, removeHs=False)
    probe = Chem.MolFromMolBlock(sdf_string2, removeHs=False)
    rms = AllChem.GetBestRMS(ref, probe)
    return rms


def lowest_cartesian(string1, string2, **linked_strings):
    """Select lowest Cartesian RMS for two structures (for nonchiral and
    previously optimized structures)."""
    values = []
    get_cartesian_rms(string1, string2)
    values.append(get_cartesian_rms(string1, string2))
    if linked_strings:
        for string in linked_strings:
            values.append(get_cartesian_rms(string1, string))

    return min(values)


def find_one_in_list(sum_array, list_to_search):
    """Generate a random number and return the corresponding index from a
    list. See the description of the method find_two_in_list."""
    nparray_to_search = np.array(list_to_search)
    rn = sum_array*np.random.rand()
    found = False
    index = 0
    while not found:
        if rn <= nparray_to_search[:index+1].sum(axis=0):
            found = True
        else:
            index += 1
    return index


def find_two_in_list(list_sum, nparray_to_search):
    """A numpy array is mapped to a segment of a line which length is equal to
    1. The lengths of the segments are proportional to the corresponding numpy
    array values. Next, two random numbers between 0 and 1 are generated and
    the segments containing these random numbers are returned."""
    rn1 = list_sum*np.random.rand()
    found1 = False
    index1 = 0
    while not found1:
        if rn1 < nparray_to_search[:index1+1].sum(axis=0):
            found1 = True
        else:
            index1 += 1
    equal = True
    while equal:
        rn2 = list_sum*np.random.rand()
        found2 = False
        index2 = 0
        while not found2:
            if rn2 < nparray_to_search[:index2+1].sum(axis=0):
                found2 = True
            else:
                index2 += 1
        if index2 != index1:
            equal = False
    return index1, index2


def find_closest(numb, list_of_values, periodic=False):
    """For a given number, return the closest value(s) from a given list"""
    all_dist = []
    for value in list_of_values:
        if periodic:
            all_dist.append(min(abs(numb-value), (360-abs(numb-value))))
        else:
            all_dist.append(abs(numb-value))
    m = min(all_dist)
    closest_ind = [i for i, j in enumerate(all_dist) if j == m]
    closest = []
    for ind in closest_ind:
        closest.append(list_of_values[ind])
    return closest


def distance(x, y):
    """"Calculate distance between two points in 3D."""
    return np.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2+(x[2]-y[2])**2)


def check_geo_sdf(sdf_string, cutoff1, cutoff2):
    """Check geometry from a sdf_string for clashes.

    Args:
        sdf_string (str)
        distance_cutoff_1 (float): min distance between non-bonded atoms [A]
        distance_cutoff_2 (float): max distance between bonded atoms [A]
    Returns:
        True for clash-free geometries and False for invalid geometries
    Raises:
        ValueError: if distance cutoffs are non-positive
    """

    if cutoff1 <= 0 or cutoff2 <= 0:
        raise ValueError("Distance cutoff needs to be a positive float")
    atoms, bonds = get_ind_from_sdfline(sdf_string.split('\n')[3])
    coordinates = np.zeros((atoms, 3))
    bonds_list = []
    for i in range(4, atoms+4):
        coordinates[i-4][0:3] = sdf_string.split('\n')[i].split()[0:3]
    for i in range(atoms+4, atoms+bonds+4):
        e1, e2 = get_ind_from_sdfline(sdf_string.split('\n')[i])
        bonds_list.append([e1, e2])
    dist = np.zeros((atoms, atoms))
    for x in range(atoms):
        for y in range(x, atoms):
            dist[x][y] = distance(np.array(coordinates[x]),
                                  np.array(coordinates[y]))
            dist[y][x] = dist[x][y]
    check = True
    for x in range(atoms):
        for y in range(x+1, atoms):
            if [x+1, y+1] not in bonds_list and [y+1, x+1] not in bonds_list:
                if dist[x][y] < cutoff1:
                    check = False
                    return check
            else:
                if dist[x][y] > cutoff2:
                    check = False
                    return check
    return check


def get_ind_from_sdfline(sdf_line):
    """Extract the indicies from the sdf string (for molecules with more than
    99 atoms)"""
    l = len(sdf_line.split()[0])
    if l < 4:
        ind1 = int(sdf_line.split()[0])
        ind2 = int(sdf_line.split()[1])
    else:
        list_ind = list(sdf_line.split()[0])
        if len(list_ind) == 5:
            ind1 = int(list_ind[0]+list_ind[1])
            ind2 = int(list_ind[2]+list_ind[3]+list_ind[4])
        if len(list_ind) == 6:
            ind1 = int(list_ind[0]+list_ind[1]+list_ind[2])
            ind2 = int(list_ind[3]+list_ind[4]+list_ind[5])

    return ind1, ind2


# Format conversions


def sdf2aims(sdf_string):
    """Convert a sdf string to a aims string."""
    atoms = get_ind_from_sdfline(sdf_string.split('\n')[3])[0]
    coord = []
    for i in range(4, 4+atoms):
        x = float(sdf_string.split('\n')[i].split()[0])
        y = float(sdf_string.split('\n')[i].split()[1])
        z = float(sdf_string.split('\n')[i].split()[2])
        name = sdf_string.split('\n')[i].split()[3]
        coord.append('%s%10.4f%10.4f%10.4f%4s' % ('atom', x, y, z, name))
        coord.append('\n')
    aims_string = ''.join(coord)
    return aims_string


def sdf2xyz(sdf_string):
    """Convert a sdf string to a xyz string."""
    atoms = get_ind_from_sdfline(sdf_string.split('\n')[3])[0]
    coord = [str(atoms)+('\n')]
    for i in range(4, 4+atoms):
        x = float(sdf_string.split('\n')[i].split()[0])
        y = float(sdf_string.split('\n')[i].split()[1])
        z = float(sdf_string.split('\n')[i].split()[2])
        name = sdf_string.split('\n')[i].split()[3]
        coord.append('\n%2s%10.4f%10.4f%10.4f' % (name, x, y, z))
    coord.append('\n')
    xyz_string = ''.join(coord)
    return xyz_string


def aims2sdf(aims_string, sdf_template_string):
    """Convert a aims string to a sdf string. Template for the sdf string is
    required."""
    atoms = len(aims_string.splitlines())
    sdf_form = sdf_template_string.splitlines()
    c = []
    cnt = 0
    for i in range(len(sdf_form)):
        if i > 3 and i < 4+atoms:
            line = sdf_form[i].split()
            line[0] = aims_string.split()[5*cnt+1]
            line[1] = aims_string.split()[5*cnt+2]
            line[2] = aims_string.split()[5*cnt+3]
            cnt += 1
            c.append('%10.4f%10.4f%10.4f%s%-2s' % (float(line[0]),
                                                   float(line[1]),
                                                   float(line[2]), str(' '),
                                                   line[3]))
            for j in range(4, len(line)):
                if j == 4:
                    c.append('%3d' % int(line[j]))
                elif j == len(line)-1:
                    c.append('%3d\n' % int(line[j]))
                else:
                    c.append('%3d' % int(line[j]))
        else:
            c.append(''.join(sdf_form[i])+'\n')

    sdf_string = ''.join(c)
    return sdf_string


def xyz2sdf(xyz_string, sdf_template_string):
    """Convert a xyz string to a sdf string. Template for the sdf string is
    required."""
    arr = xyz_string.splitlines()
    atoms = int(arr[0].split()[0])
    xyz_string_cut = '\n'.join(arr[2:])
    sdf_form = sdf_template_string.splitlines()
    c = []
    cnt = 0
    for i in range(len(sdf_form)):
        if i > 3 and i < 4+atoms:
            line = sdf_form[i].split()
            line[0] = xyz_string_cut.split()[4*cnt+1]
            line[1] = xyz_string_cut.split()[4*cnt+2]
            line[2] = xyz_string_cut.split()[4*cnt+3]
            cnt += 1
            c.append('%10.4f%10.4f%10.4f%s%-2s' % (float(line[0]),
                                                   float(line[1]),
                                                   float(line[2]), str(' '),
                                                   line[3]))
            for j in range(4, len(line)):
                if j == 4:
                    c.append('%3d' % int(line[j]))
                elif j == len(line)-1:
                    c.append('%3d\n' % int(line[j]))
                else:
                    c.append('%3d' % int(line[j]))
        else:
            c.append(''.join(sdf_form[i])+'\n')

    sdf_string = ''.join(c)
    return sdf_string


def mirror_sdf(sdf_string):
    """Mirror the geometry from a sdf string. Return a new sdf string."""
    atoms = get_ind_from_sdfline(sdf_string.split('\n')[3])[0]
    sdf_form = sdf_string.splitlines()
    c = []
    cnt = 0
    for i in range(len(sdf_form)):
        if i > 3 and i < 4+atoms:
            line = sdf_form[i].split()
            line[0] = -1.0*float(line[0])
            line[1] = -1.0*float(line[1])
            line[2] = -1.0*float(line[2])
            cnt += 1
            c.append('%10.4f%10.4f%10.4f%s%-2s' % (float(line[0]),
                                                   float(line[1]),
                                                   float(line[2]), str(' '),
                                                   line[3]))
            for j in range(4, len(line)):
                if j == 4:
                    c.append('%3d' % int(line[j]))
                elif j == len(line)-1:
                    c.append('%3d\n' % int(line[j]))
                else:
                    c.append('%3d' % int(line[j]))
        else:
            c.append(''.join(sdf_form[i])+'\n')
    mirror_sdf_string = ''.join(c)
    return mirror_sdf_string
