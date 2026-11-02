import sys
import re
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdmolops
from qcforever.util.quadrupole import BOHR_TO_ANGSTROM


def read_fchk(infile):
    """Read molecular geometry/state directly from Gaussian formatted text.

    No Gaussian executable is used. Atomic numbers identify elements even when
    Nuclear charges have been reduced by ECPs. Orbitals/basis/energies are not
    imported into another backend.
    """
    required = {'Number of atoms', 'Charge', 'Multiplicity', 'Atomic numbers',
                'Current cartesian coordinates'}
    names = {name.lower(): name for name in required | {'Nuclear charges'}}
    pattern = re.compile(r'^\s*(' + '|'.join(re.escape(name) for name in names) +
                         r')\s+([IR])\s+(.*?)\s*$', re.I)
    number = re.compile(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DdEe][+-]?\d+)?')
    fields = {}
    with open(infile, encoding='utf-8', errors='replace') as stream:
        # The first two lines are free-form title and method/basis descriptors.
        next(stream, None)
        next(stream, None)
        for line in stream:
            match = pattern.fullmatch(line)
            if match is None:
                continue
            label, kind, content = names[match[1].lower()], match[2].upper(), match[3]
            if label in fields:
                raise ValueError(f'{infile}: duplicate fchk field {label}')
            array = label in {'Atomic numbers', 'Nuclear charges', 'Current cartesian coordinates'}
            expected_kind = 'R' if label in {'Nuclear charges', 'Current cartesian coordinates'} else 'I'
            if kind != expected_kind:
                raise ValueError(f'{infile}: invalid fchk type for {label}')
            if array:
                count_match = re.fullmatch(r'N=\s*(\d+)', content, re.I)
                if count_match is None:
                    raise ValueError(f'{infile}: missing fchk array count for {label}')
                count = int(count_match[1])
                values = []
                while len(values) < count:
                    row = next(stream, '')
                    if not row or not row.strip() or number.sub('', row).strip():
                        raise ValueError(f'{infile}: truncated or invalid fchk array {label}')
                    tokens = number.findall(row)
                    try:
                        values.extend(int(value) if kind == 'I' else float(value.replace('D', 'E').replace('d', 'e'))
                                      for value in tokens)
                    except ValueError as exc:
                        raise ValueError(f'{infile}: invalid fchk number in {label}') from exc
                if len(values) != count:
                    raise ValueError(f'{infile}: wrong fchk array length for {label}')
                fields[label] = values
            else:
                try:
                    fields[label] = int(content)
                except ValueError as exc:
                    raise ValueError(f'{infile}: invalid fchk integer for {label}') from exc
            if required <= fields.keys():
                break
    missing = required - fields.keys()
    if missing:
        raise ValueError(f'{infile}: missing fchk fields: {", ".join(sorted(missing))}')
    count = fields['Number of atoms']
    numbers = fields['Atomic numbers']
    coordinates = np.asarray(fields['Current cartesian coordinates'], dtype=float)
    if count < 1 or len(numbers) != count or coordinates.size != 3 * count:
        raise ValueError(f'{infile}: inconsistent fchk atom/coordinate counts')
    if not all(1 <= value <= 118 for value in numbers):
        raise ValueError(f'{infile}: unsupported ghost/dummy atom in fchk Atomic numbers')
    if not np.all(np.isfinite(coordinates)):
        raise ValueError(f'{infile}: nonfinite fchk coordinates')
    nuclear = fields.get('Nuclear charges')
    if nuclear is not None and (len(nuclear) != count or not all(np.isfinite(value) and value > 0 for value in nuclear)):
        raise ValueError(f'{infile}: unsupported nuclear charges or ghost centers in fchk')
    charge, mult = fields['Charge'], fields['Multiplicity']
    if mult < 1:
        raise ValueError(f'{infile}: fchk multiplicity must be positive')
    atoms = [Chem.GetPeriodicTable().GetElementSymbol(value) for value in numbers]
    coordinates = coordinates.reshape(count, 3) * BOHR_TO_ANGSTROM
    return atoms, coordinates[:, 0], coordinates[:, 1], coordinates[:, 2], charge, mult


def write_xyz(outfile, atoms, x, y, z, charge, multiplicity):
    """Write a structure handoff without changing electronic-state overrides."""
    if any(int(value) != value for value in (charge, multiplicity)) or multiplicity < 1:
        raise ValueError('XYZ charge and multiplicity must be integers, with positive multiplicity')
    coordinates = np.column_stack((x, y, z))
    if len(coordinates) != len(atoms) or not np.all(np.isfinite(coordinates)):
        raise ValueError('Invalid XYZ handoff coordinates')
    Path(outfile).write_text(f'{len(atoms)}\n{int(charge)} {int(multiplicity)}\n' + ''.join(
        f'{atom} {xx:.12f} {yy:.12f} {zz:.12f}\n'
        for atom, (xx, yy, zz) in zip(atoms, coordinates)))


def read_sdf(infile):
    sdfinfo = Chem.SDMolSupplier(infile, removeHs=False)
    mols = [mol for mol in sdfinfo if mol is not None]
    N = mols[0].GetNumAtoms()
    N_Bond = mols[0].GetNumBonds()
    N_radical = Descriptors.NumRadicalElectrons(mols[0])
    Charge = rdmolops.GetFormalCharge(mols[0])

    with open(infile, 'r') as ifile:
        count = 0
        X = []
        Y = []
        Z = []
        element_symbol = []
        Bond_pair1 = []
        Bond_pair2 = []
        Bond_type = []
        sdfCharge = 0 
        CHG_atom = []
        CHG = []
        sdfSpinMulti = 0
        RAD_atom =[]
        RAD = []
        
        for line in ifile:
            if count == 0:
                Header1 = line  # not used
                count += 1
                continue
            if count == 1:
                Header2 = line  # not used
                count += 1
                continue
            if count == 2:
                Header3 = line  # not used
                count += 1
                continue
            if count == 3:
                a = line.split()  # not used
                count += 1
                continue
            if 3 < count <= N+3:
                #print("Atom Coordinte")
                #print(line)
                i_atom = line.split()
                if len(i_atom) != 0:
                    X.append(float(i_atom[0]))
                    Y.append(float(i_atom[1]))
                    Z.append(float(i_atom[2]))
                    element_symbol.append(i_atom[3])
                count += 1
                continue
            if N+3 < count <= N+N_Bond+3:
                bond_info = line.split()
                bond_info = line.split()
                Bond_pair1.append(int(bond_info[0]))
                Bond_pair2.append(int(bond_info[1]))
                Bond_type.append(int(bond_info[2]))
                count += 1
                continue
            if count > N+N_Bond+3:
                mol_info = line.split()
                if (mol_info[0] == "M"):
                    if (mol_info[1] == "END"):
                        break
                    if (mol_info[1] == "CHG"):
                        Num_CHGInfo = int(mol_info[2])
                        for k in range(Num_CHGInfo):
                            CHG_atom.append(int(mol_info[3+2*k]))
                            CHG.append(int(mol_info[4+2*k]))
                            sdfCharge += int(mol_info[4+2*k])
                    if (mol_info[1] == "RAD"):
                        Num_RADInfo = int(mol_info[2])
                        for l in range(Num_RADInfo):
                            RAD_atom.append(int(mol_info[3+2*l]))
                            RAD.append(int(mol_info[4+2*l]))
                            sdfSpinMulti += int(mol_info[4+2*l])
                else:
                    print("The sdf file is invalid!")
                    sys.exit()
                count +=1
        
    # Copy to array of numpy
    CHG_atom = np.array(CHG_atom)
    CHG = np.array(CHG)
    Mol_atom = element_symbol
    Mol_CartX = np.array(X)
    Mol_CartY = np.array(Y)
    Mol_CartZ = np.array(Z)

    print('Reading the sdf file has finished')
    SpinMulti = N_radical + 1
    return Mol_atom, Mol_CartX, Mol_CartY, Mol_CartZ, Charge, SpinMulti, Bond_pair1, Bond_pair2


def read_xyz(infile):
    with open(infile, 'r') as ifile:
        count = 0
        X = []
        Y = []
        Z = []
        element_symbol = []
        Charge = 0 
        SpinMulti = 0
        for line in ifile:
        #    print(line)
            if count == 0:
                Atom_num = line.split()
                N = int(Atom_num[0])
                count += 1
                continue
            if count == 1:
                Comment = line.split()
                Charge = int(Comment[0])
                SpinMulti = int (Comment[1])
                count += 1
                continue
            if 1 < count:
                i_atom = line.split()
                if len(i_atom) > 3:
                    element_symbol.append(i_atom[0])
                    X.append(float(i_atom[1]))
                    Y.append(float(i_atom[2]))
                    Z.append(float(i_atom[3]))
                    count += 1
                    continue
                else:
                    count = 0
                    break
    # Copy to array of numpy
    Mol_atom = []
    Mol_CartX = np.zeros(N)
    Mol_CartY = np.zeros(N)
    Mol_CartZ = np.zeros(N)
    for j in range(N):
        Mol_CartX[j] = X[j] 
        Mol_CartY[j] = Y[j] 
        Mol_CartZ[j] = Z[j] 
        Mol_atom.append(element_symbol[j])

    print('Reading the xyz file has finished')
    return Mol_atom, Mol_CartX, Mol_CartY, Mol_CartZ, Charge, SpinMulti
