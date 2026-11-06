"""Native, lossless FCHK orbital handoff for all-electron s/p/d bases.

FCHK coefficients are MO-major. Gaussian and ORCA use d(0,+1,-1,+2,-2);
GAMESS stores even its spherical orbitals on Cartesian xx,yy,zz,xy,xz,yz.
No projection between different AO spaces is attempted here.
"""

from dataclasses import dataclass
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
from rdkit import Chem

from . import job_timeout, read_mol_file


# Columns: Gaussian 5D; rows: normalized Cartesian xx,yy,zz,xy,xz,yz.
D_SPH_TO_CART = np.array([
    [-.5, 0, 0, np.sqrt(3)/2, 0],
    [-.5, 0, 0, -np.sqrt(3)/2, 0],
    [1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1],
    [0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0],
])

_SCALARS = {'Number of basis functions', 'Number of independent functions',
            'Number of alpha electrons', 'Number of beta electrons'}
_INT_ARRAYS = {'Shell types', 'Number of primitives per shell', 'Shell to atom map'}
_REAL_ARRAYS = {'Primitive exponents', 'Contraction coefficients',
                'P(S=P) Contraction coefficients', 'Alpha MO coefficients',
                'Beta MO coefficients', 'Alpha Orbital Energies', 'Beta Orbital Energies',
                'Nuclear charges'}
_NUM = re.compile(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DdEe][+-]?\d+)?')


def _fields(path):
    names = _SCALARS | _INT_ARRAYS | _REAL_ARRAYS
    pattern = re.compile(r'^(' + '|'.join(re.escape(n) for n in names) + r')\s+([IR])\s+(.*?)\s*$')
    fields = {}
    with Path(path).open() as stream:
        next(stream, None)
        next(stream, None)
        for row in stream:
            if row.startswith('ECP-'):
                raise ValueError('fchk orbital import does not yet support ECPs')
            match = pattern.fullmatch(row)
            if not match:
                continue
            name, kind, value = match.groups()
            if name in fields or kind != ('R' if name in _REAL_ARRAYS else 'I'):
                raise ValueError(f'Invalid or duplicate fchk field: {name}')
            convert = int if kind == 'I' else lambda s: float(s.replace('D', 'E').replace('d', 'e'))
            if name in _SCALARS:
                fields[name] = convert(value)
                continue
            count = re.fullmatch(r'N=\s*(\d+)', value)
            if count is None:
                raise ValueError(f'Invalid fchk array count: {name}')
            values = []
            while len(values) < int(count[1]):
                row = next(stream, '')
                if not row.strip() or _NUM.sub('', row).strip():
                    raise ValueError(f'Truncated or invalid fchk array: {name}')
                values.extend(convert(s) for s in _NUM.findall(row))
            if len(values) != int(count[1]) or not np.all(np.isfinite(values)):
                raise ValueError(f'Invalid fchk array length or values: {name}')
            fields[name] = np.asarray(values)
    required = (_SCALARS - {'Number of independent functions'}) | _INT_ARRAYS | {
        'Primitive exponents', 'Contraction coefficients', 'Alpha MO coefficients', 'Alpha Orbital Energies'}
    missing = required - fields.keys()
    if missing:
        raise ValueError('Missing fchk orbital fields: ' + ', '.join(sorted(missing)))
    return fields


@dataclass
class Shell:
    atom: int
    kind: int
    exponents: np.ndarray
    coefficients: np.ndarray
    offset: int

    @property
    def size(self):
        return {0: 1, 1: 3, -2: 5, 2: 6}[self.kind]

    def basis_lines(self):
        lines = [f' {"SPD"[abs(self.kind)]} {len(self.exponents)}']
        lines.extend(f' {i} {a:.15E} {c:.15E}'
                     for i, (a, c) in enumerate(zip(self.exponents, self.coefficients), 1))
        return lines


@dataclass
class Orbitals:
    atoms: list
    xyz: np.ndarray
    charge: int
    multiplicity: int
    shells: list
    alpha: np.ndarray
    beta: np.ndarray | None
    alpha_energy: np.ndarray
    beta_energy: np.ndarray | None
    nalpha: int
    nbeta: int
    spherical: bool

    @property
    def nmo(self):
        return self.alpha.shape[1]

    def validate_target(self, backend, angular='auto'):
        if angular not in ('auto', 'spherical', 'cartesian'):
            raise ValueError('fchk_angular must be auto, spherical or cartesian')
        if any(abs(s.kind) == 2 for s in self.shells):
            if angular != 'auto' and (angular == 'spherical') != self.spherical:
                raise ValueError('fchk orbital import requires the source angular basis; 5D/6D projection is not supported')
            if backend == 'orca' and not self.spherical:
                raise ValueError('ORCA requires 5D orbitals; regenerate the Gaussian fchk with 5D (and 7F for f shells)')
        if backend == 'orca' and angular == 'cartesian':
            raise ValueError('ORCA supports spherical basis functions only')

    def gamess_coordinates(self, atoms, xyz):
        if list(atoms) != self.atoms:
            raise ValueError('fchk orbital atom order/composition changed')
        lines = []
        for i, (atom, coord) in enumerate(zip(atoms, xyz)):
            z = Chem.GetPeriodicTable().GetAtomicNumber(atom)
            lines.append(f'{atom} {z:.1f} ' + ' '.join(f'{v:.12f}' for v in coord))
            for shell in self.shells:
                if shell.atom == i:
                    lines.extend(shell.basis_lines())
            lines.append('')
        return '\n'.join(lines) + '\n'

    def gamess_vec(self):
        lines = [' $VEC']
        # QCforever's ground-state GAMESS calculation uses UHF, including for
        # singlets. A restricted source supplies the same spatial MOs twice.
        for coeff in (self.alpha, self.beta if self.beta is not None else self.alpha):
            blocks = [D_SPH_TO_CART @ coeff[s.offset:s.offset+5]
                      if s.kind == -2 else coeff[s.offset:s.offset+s.size]
                      for s in self.shells]
            cart = np.vstack(blocks)
            for mo in range(self.nmo):
                for start in range(0, len(cart), 5):
                    lines.append(f'{(mo+1)%100:2d}{start//5+1:3d}' +
                                 ''.join(f'{v:15.8E}' for v in cart[start:start+5, mo]))
        return '\n'.join(lines) + '\n $END\n'

    def orca_basis(self, atom):
        lines = [' NewGTO']
        for shell in self.shells:
            if shell.atom == atom:
                lines.extend(shell.basis_lines())
        return lines + [' end']

    def write_mkl(self, path):
        self.validate_target('orca')
        lines = ['$MKL', '# QCforever FCHK orbital handoff', '$CHAR_MULT',
                 f'{self.charge} {self.multiplicity}', '$END', '$COORD']
        for atom, coord in zip(self.atoms, self.xyz):
            lines.append(str(Chem.GetPeriodicTable().GetAtomicNumber(atom)) + ' ' +
                         ' '.join(f'{v:.12f}' for v in coord))
        lines += ['$END', '$CHARGES'] + ['0.0'] * len(self.atoms) + ['$END', '$BASIS']
        for i in range(len(self.atoms)):
            if i:
                lines.append('$$')
            for shell in self.shells:
                if shell.atom == i:
                    lines.append(f'{2*abs(shell.kind)+1:2d} {"SPD"[abs(shell.kind)]} 1.0')
                    lines.extend(f'{a:.15E} {c:.15E}' for a, c in zip(shell.exponents, shell.coefficients))
        lines += ['', '$END']
        for spin, coeff, energy, occupied in (
            ('ALPHA', self.alpha, self.alpha_energy, self.nalpha),
            ('BETA', self.beta, self.beta_energy, self.nbeta),
        ):
            if coeff is None:
                continue
            lines.append('$COEFF_' + spin)
            for start in range(0, self.nmo, 5):
                stop = min(start+5, self.nmo)
                # ORCA 6.1's MKL reader requires its three-character symmetry
                # label here; "a1" can leave the converter waiting indefinitely.
                lines += [' '.join(['a1g']*(stop-start)),
                          ' '.join(f'{v:.15E}' for v in energy[start:stop])]
                lines.extend(' '.join(f'{v:.15E}' for v in row) for row in coeff[:, start:stop])
            occupation = np.zeros(self.nmo)
            occupation[:occupied] = 1
            if self.beta is None:
                occupation[:self.nbeta] += 1
            lines += ['$END', '$OCC_' + spin]
            lines.extend(' '.join(f'{v:.8f}' for v in occupation[i:i+5]) for i in range(0, self.nmo, 5))
            lines.append('$END')
        Path(path).write_text(('\n'.join(lines) + '\n').replace('$END\n', '$END\n\n'))

    def orca_gbw(self, executable, stem='fchk_guess'):
        # orca_2mkl belongs to ORCA itself; no third-party converter is needed.
        resolved = shutil.which(str(executable))
        candidate = Path(resolved).resolve().with_name('orca_2mkl') if resolved else None
        converter = str(candidate) if candidate and candidate.is_file() else shutil.which('orca_2mkl')
        if not converter:
            raise ValueError('FCHK orbital import needs ORCA\'s bundled orca_2mkl executable')
        self.write_mkl(stem + '.mkl')
        Path(stem + '.gbw').unlink(missing_ok=True)
        result = job_timeout.run([converter, stem, '-gbw'],
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        logfile = stem + '.conversion.log'
        Path(logfile).write_text(result.stdout or '')
        if result.returncode:
            raise RuntimeError(f'orca_2mkl failed (exit {result.returncode}); see {logfile}')
        if not Path(stem + '.gbw').is_file() or not Path(stem + '.gbw').stat().st_size:
            raise RuntimeError(f'orca_2mkl did not produce {stem}.gbw')
        return stem + '.gbw'


def read(path):
    atoms, x, y, z, charge, mult = read_mol_file.read_fchk(path)
    f = _fields(path)
    nbf = int(f['Number of basis functions'])
    nmo = int(f.get('Number of independent functions', nbf))
    na, nb = int(f['Number of alpha electrons']), int(f['Number of beta electrons'])
    if not 0 <= nb <= na <= nmo <= nbf or nmo < 1:
        raise ValueError('Invalid fchk electron/orbital counts')
    if na+nb != sum(Chem.GetPeriodicTable().GetAtomicNumber(a) for a in atoms)-charge or na-nb != mult-1:
        raise ValueError('fchk electron counts do not match charge/multiplicity (ECPs are unsupported)')
    if 'Nuclear charges' in f and not np.array_equal(f['Nuclear charges'], [Chem.GetPeriodicTable().GetAtomicNumber(a) for a in atoms]):
        raise ValueError('fchk orbital import does not yet support ECPs')
    kinds, counts, centers = (f[n] for n in ('Shell types', 'Number of primitives per shell', 'Shell to atom map'))
    if not len(kinds) or len(kinds) != len(counts) or len(kinds) != len(centers) or np.any(counts < 1):
        raise ValueError('Inconsistent fchk shell counts')
    if any(k not in (0, 1, -1, -2, 2) for k in kinds):
        raise ValueError('Native fchk orbital import currently supports s, p, sp and d shells; f and higher are unsupported')
    if -2 in kinds and 2 in kinds:
        raise ValueError('Mixed 5D/6D fchk shells are unsupported')
    if set(centers) != set(range(1, len(atoms)+1)) or np.any(np.diff(centers) < 0):
        raise ValueError('fchk shell centers must follow atom order and cover every atom')
    nprim = int(sum(counts))
    for name in ('Primitive exponents', 'Contraction coefficients'):
        if len(f[name]) != nprim:
            raise ValueError('Inconsistent fchk primitive counts')
    if np.any(f['Primitive exponents'] <= 0):
        raise ValueError('fchk primitive exponents must be positive')
    if -1 in kinds and len(f.get('P(S=P) Contraction coefficients', [])) != nprim:
        raise ValueError('Missing/inconsistent fchk SP contraction coefficients')
    shells, prim, offset = [], 0, 0
    for kind, count, center in zip(kinds, counts, centers):
        exponents = f['Primitive exponents'][prim:prim+count]
        components = [(0, 'Contraction coefficients'), (1, 'P(S=P) Contraction coefficients')] if kind == -1 else [(kind, 'Contraction coefficients')]
        for angular, name in components:
            c = f[name][prim:prim+count]
            if not np.any(c) or (count == 1 and c[0] < 0):
                raise ValueError('Unsupported zero or negative single-primitive fchk contraction')
            shell = Shell(int(center)-1, int(angular), exponents, c, offset)
            shells.append(shell)
            offset += shell.size
        prim += count
    if offset != nbf:
        raise ValueError('fchk basis dimension does not match shell types')
    def coefficients(spin):
        values = f.get(spin + ' MO coefficients')
        energy = f.get(spin + ' Orbital Energies')
        if values is None and energy is None and spin == 'Beta':
            return None, None
        if values is None or energy is None or values.size != nbf*nmo or energy.size != nmo:
            raise ValueError(f'Inconsistent fchk {spin} orbital dimensions')
        return values.reshape(nmo, nbf).T.copy(), energy
    alpha, ea = coefficients('Alpha')
    beta, eb = coefficients('Beta')
    return Orbitals(atoms, np.column_stack((x,y,z)), charge, mult, shells, alpha, beta,
                    ea, eb, na, nb, 2 not in kinds)


def for_job(job, backend, options):
    """basis='FCHK' explicitly selects source basis and an orbital guess."""
    if str(job.basis).upper() != 'FCHK':
        return None
    if Path(job.in_file).suffix.lower() not in ('.fchk', '.fch'):
        raise ValueError('basis=FCHK requires a Gaussian formatted checkpoint')
    if 'optconf' in options or str(job.functional).upper().startswith('KTLC-'):
        raise ValueError('basis=FCHK orbital import does not support optconf or KTLC tuning')
    orbitals = read(job.in_file)
    orbitals.validate_target(backend, getattr(job, 'fchk_angular', 'auto'))
    if backend == 'orca' and 'nmr' in options:
        required = set(orbitals.atoms) & {'H', 'C', 'Si'}
        if required - getattr(job, 'nmr_reference', {}).keys():
            raise ValueError('basis=FCHK NMR needs explicit nmr_reference shieldings; the imported atom-specific basis cannot be assigned to TMS automatically')
    for attribute, expected in (('SpecTotalCharge', orbitals.charge), ('SpecSpinMulti', orbitals.multiplicity)):
        value = getattr(job, attribute)
        if not np.isnan(value) and value != expected:
            raise ValueError('basis=FCHK orbital import requires the source charge and multiplicity')
    return orbitals
