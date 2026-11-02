"""Read ORCA 6.1 output, keeping units and the final property block explicit."""

from pathlib import Path
import re

import numpy as np

from qcforever.util.quadrupole import (
    BOHR_TO_ANGSTROM, E_ANGSTROM_TO_DEBYE,
    quadrupole_properties, tensor_from_components,
)

EH_TO_EV = 27.211386245988
AU_TO_DEBYE = BOHR_TO_ANGSTROM * E_ANGSTROM_TO_DEBYE
AU_TO_DEBYE_ANGSTROM = AU_TO_DEBYE * BOHR_TO_ANGSTROM
AU_TO_ANGSTROM3 = BOHR_TO_ANGSTROM ** 3
NUMBER = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?'


def number(value):
    result = float(value.replace('D', 'E').replace('d', 'e'))
    if not np.isfinite(result):
        raise ValueError('Nonfinite ORCA property')
    return result


class parse_log:
    def __init__(self, filename):
        self.filename = str(filename)
        self.text = Path(filename).read_text(errors='replace')

    def _last(self, pattern, text=None):
        matches = list(re.finditer(pattern, self.text if text is None else text, re.M))
        if not matches:
            raise ValueError(f'{self.filename}: missing {pattern!r}')
        return matches[-1]

    def _block(self, heading):
        start = self._last(r'^\s*' + re.escape(heading) + r'\s*$').end()
        return self.text[start:]

    def _value(self, label, text=None):
        return number(self._last(r'^\s*' + re.escape(label) +
                                r'\s*(?:\.{2,}|:)?\s*(' + NUMBER + r')', text).group(1))

    def check(self, optimize=False):
        if 'ORCA TERMINATED NORMALLY' not in self.text:
            raise ValueError(f'{self.filename}: ORCA did not terminate normally')
        if re.search(r'SCF NOT CONVERGED|SCF did not converge', self.text, re.I):
            raise ValueError(f'{self.filename}: SCF did not converge')
        if optimize and 'THE OPTIMIZATION HAS CONVERGED' not in self.text:
            raise ValueError(f'{self.filename}: geometry optimization did not converge')
        self.energy()

    def energy(self):
        energy = number(self._last(r'^\s*FINAL SINGLE POINT ENERGY\s+(' + NUMBER + ')').group(1))
        spin = re.findall(r'Expectation value of <S\*\*2>\s*[:=]\s*(' + NUMBER + ')', self.text)
        # The second item follows the existing energy convention: <S²>-S(S+1).
        if spin:
            mult = self._value('Multiplicity           Mult')
            s = (mult - 1) / 2
            contamination = number(spin[-1]) - s * (s + 1)
        else:
            contamination = 0.0
        return [energy, contamination]

    def coordinates(self):
        block = self._block('CARTESIAN COORDINATES (ANGSTROEM)')
        atoms, xyz = [], []
        for line in block.splitlines():
            match = re.fullmatch(r'\s*([A-Z][a-z]?)\s+(' + NUMBER + r')\s+(' +
                                 NUMBER + r')\s+(' + NUMBER + r')\s*', line)
            if match:
                atoms.append(match[1])
                xyz.append([number(match[i]) for i in (2, 3, 4)])
            elif atoms:
                break
        if not atoms:
            raise ValueError('Missing final Cartesian coordinates')
        return atoms, np.array(xyz)

    def orbitals(self):
        block = self._block('ORBITAL ENERGIES')
        channels, rows = [], []
        for line in block.splitlines():
            match = re.fullmatch(r'\s*\d+\s+(' + NUMBER + r')\s+(' + NUMBER +
                                 r')\s+(' + NUMBER + r')\s*', line)
            if match:
                rows.append((number(match[1]), number(match[2])))
            elif rows:
                channels.append(rows)
                rows = []
                if len(channels) == 2:
                    break
            if channels and 'MULLIKEN' in line:
                break
        if rows:
            channels.append(rows)
        if not channels:
            raise ValueError('Missing orbital energies')
        result, gaps = {}, []
        for index, channel in enumerate(channels):
            occupied = [energy for occupation, energy in channel if occupation > 1e-6]
            virtual = [energy for occupation, energy in channel if occupation <= 1e-6]
            if not occupied or not virtual:
                raise ValueError('HOMO/LUMO requires occupied and virtual orbitals in each spin channel')
            gaps.append((min(virtual) - max(occupied)) * EH_TO_EV)
            result['MO' if len(channels) == 1 else ('Alpha_MO', 'Beta_MO')[index]] = [occupied, virtual]
        result['homolumo'] = gaps[0] if len(gaps) == 1 else gaps
        return result

    def charges(self):
        match = self._last(r'^\s*MULLIKEN ATOMIC CHARGES(?: AND SPIN POPULATIONS)?\s*$')
        elements, charges, spins = [], [], []
        for line in self.text[match.end():].splitlines():
            row = re.fullmatch(r'\s*\d+\s+([A-Z][a-z]?)\s*:\s*(' + NUMBER +
                               r')(?:\s+(' + NUMBER + r'))?\s*', line)
            if row:
                elements.append(row[1])
                charges.append(number(row[2]))
                spins.append(number(row[3]) if row[3] else 0.0)
            elif elements:
                break
        if not elements:
            raise ValueError('Missing Mulliken populations')
        return [elements, charges, spins]

    def dipole(self):
        block = self._block('DIPOLE MOMENT')
        match = re.search(r'Total Dipole Moment\s*:\s*(' + NUMBER + r')\s+(' +
                          NUMBER + r')\s+(' + NUMBER + r')', block)
        if not match:
            raise ValueError('Missing total dipole moment')
        vector = [number(match[i]) * AU_TO_DEBYE for i in (1, 2, 3)]
        magnitude = re.search(r'Magnitude \(Debye\)\s*:\s*(' + NUMBER + ')', block)
        total = number(magnitude[1]) if magnitude else float(np.linalg.norm(vector))
        # ORCA uses its own a.u./Debye conversion constant in printed output.
        norm = np.linalg.norm(vector)
        if norm:
            vector = [value * total / norm for value in vector]
        return vector + [total]

    def quadrupole(self):
        # make_input explicitly requests CenterOfMass. ORCA's TOT (a.u.) is
        # the primitive Cartesian charge second moment, not its traceless form.
        try:
            block = self._block('QUADRUPOLE MOMENT')
            row = re.search(r'^\s*TOT\s+(' + NUMBER + r'(?:\s+' + NUMBER + r'){5})', block, re.M)
            if row is None:
                raise ValueError('Missing complete final quadrupole tensor')
            tensor = tensor_from_components(dict(zip(
                ('XX', 'YY', 'ZZ', 'XY', 'XZ', 'YZ'),
                [number(value) * AU_TO_DEBYE_ANGSTROM for value in row[1].split()])))
            return quadrupole_properties(primitive=tensor, source='ORCA',
                                         origin='center_of_mass', status='ok')
        except ValueError as exc:
            return quadrupole_properties(source='ORCA', origin='center_of_mass',
                                         status='unavailable', reason=str(exc))

    def polar(self):
        block = self._block('STATIC POLARIZABILITY TENSOR (Dipole/Dipole)')
        match = re.search(r'The raw cartesian tensor \(atomic units\):\s*\n' +
                          (r'\s*(' + NUMBER + r')\s+(' + NUMBER + r')\s+(' + NUMBER + r')\s*\n') * 3, block)
        if not match:
            raise ValueError('Missing static polarizability tensor')
        tensor = np.array([number(v) for v in match.groups()]).reshape(3, 3) * AU_TO_ANGSTROM3
        iso = float(np.trace(tensor) / 3)
        aniso = float(np.sqrt(1.5 * np.sum((tensor - np.eye(3) * iso)**2)))
        return dict(polar_tens=tensor.tolist(), polar_iso=iso, polar_aniso=aniso)

    def shielding(self):
        block = self._block('CHEMICAL SHIELDING SUMMARY (ppm)')
        result = []
        for line in block.splitlines():
            match = re.fullmatch(r'\s*(\d+)\s+([A-Z][a-z]?)\s+(' + NUMBER +
                                 r')\s+(' + NUMBER + r')\s*', line)
            if match:
                result.append((int(match[1]), match[2], number(match[3])))
            elif result:
                break
        if not result:
            raise ValueError('Missing NMR shielding summary')
        return sorted(result)

    def transitions(self):
        block = self._block('ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS')
        rows = []
        # ORCA 6: Transition Energy(eV) Energy(cm-1) Wavelength(nm) fosc ...
        # ORCA 5: State Energy(cm-1) Wavelength(nm) fosc ...
        version6 = 'eV' in block.split('---', 1)[-1][:600]
        for line in block.splitlines():
            fields = line.split()
            if '->' in fields or '-->' in fields:
                arrow = fields.index('->' if '->' in fields else '-->')
                values = fields[arrow + 2:]
                state = re.match(r'(\d+)', fields[arrow + 1])
                if state and len(values) >= 4:
                    rows.append((int(state[1]), number(values[2]), number(values[3])))
            elif fields and fields[0].isdigit() and len(fields) >= 4:
                offset = 1 if version6 else 0
                rows.append((int(fields[0]), number(fields[2 + offset]), number(fields[3 + offset])))
            elif rows:
                break
        if not rows:
            raise ValueError('Missing electric-dipole absorption transitions')
        return rows

    def frequencies(self):
        block = self._block('VIBRATIONAL FREQUENCIES')
        modes = [(int(m[1]), number(m[2])) for m in re.finditer(
            r'^\s*(\d+):\s*(' + NUMBER + r')\s+cm\*\*-1', block.split('NORMAL MODES')[0], re.M)]
        if not modes:
            raise ValueError('Missing vibrational frequencies')
        spectra = []
        for title in ('IR SPECTRUM', 'RAMAN SPECTRUM'):
            section = self._block(title)
            rows = {}
            for line in section.splitlines():
                match = re.match(r'^\s*(\d+):\s+(' + NUMBER + r')\s+(.+)$', line)
                if match:
                    values = match[3].split()
                    rows[int(match[1])] = number(values[1 if title == 'IR SPECTRUM' else 0])
                elif rows:
                    break
            if not rows:
                raise ValueError(f'Missing {title} intensities')
            spectra.append(rows)
        indices = [index for index, _ in modes]
        # Keep all 3N modes; ORCA omits translations/rotations in spectrum tables.
        result = dict(freq=[value for _, value in modes],
                      IR=[spectra[0].get(i, 0.0) for i in indices],
                      Raman=[spectra[1].get(i, 0.0) for i in indices])
        # The normal-mode table is printed in groups of up to six columns.
        section = self._block('NORMAL MODES').split('IR SPECTRUM')[0]
        displacements = np.full((len(modes), len(modes)), np.nan)
        columns = []
        for line in section.splitlines():
            fields = line.split()
            if fields and all(field.isdigit() for field in fields):
                columns = [int(field) for field in fields]
            elif columns and fields and fields[0].isdigit() and len(fields) == len(columns) + 1:
                coordinate = int(fields[0])
                if coordinate < len(modes) and max(columns) < len(modes):
                    displacements[coordinate, columns] = [number(v) for v in fields[1:]]
        if not np.all(np.isfinite(displacements)):
            raise ValueError('Incomplete normal-mode displacement table')
        result['freqmode'] = [displacements[axis::3, :].T.tolist() for axis in range(3)]
        thermo = self._block('INNER ENERGY')
        electronic = self._value('Electronic energy', thermo)
        result.update(Ezp=electronic + self._value('Zero point energy', thermo),
                      Et=self._value('Total thermal energy', thermo),
                      E_enth=self._value('Total Enthalpy', thermo),
                      E_free=self._value('Final Gibbs free energy', thermo))
        # ORCA's default QRRHO entropy is retained; do not label T*S as S.
        temperature = self._value('Temperature')
        result['Si'] = self._value('Final entropy term', thermo) * 627509.4740631 / temperature
        result['Ei'] = (result['Et'] - electronic) * 627.5094740631
        return result
