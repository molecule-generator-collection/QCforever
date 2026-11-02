"""QCforever workflows using ORCA 6.1 or later.

Independent molecules must run in separate processes: optconf, like the other
QCforever backends, uses a process-local working directory.
"""

from contextlib import chdir
import math
from pathlib import Path
import pickle
import re
import shutil

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from qcforever.conformer_search.settings import calculation_tokens, parse_conformer_options
from qcforever.util import check_resource, job_cleanup, job_timeout, Spectrum_similarity, UV_similarity
from . import Exe_Orca
from .parse_log import EH_TO_EV, parse_log
from qcforever.util import calculation_cache
from qcforever.util import read_mol_file
from qcforever.util.calculation_options import parse_excited_optimization, ground_optimization, requires_hessian


def _token(value, label):
    value = str(value).strip()
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_+*/(),.\-]*', value):
        raise ValueError(f'{label} must be one ORCA keyword')
    return value


def _integer(value, label, minimum=None):
    if isinstance(value, (bool, np.bool_)) or not math.isfinite(float(value)) or int(value) != float(value):
        raise ValueError(f'{label} must be an integer')
    value = int(value)
    if minimum is not None and value < minimum:
        raise ValueError(f'{label} must be at least {minimum}')
    return value


def memory_per_core(memory, nproc):
    """Treat mem as the total budget; reserve 25% for unaccounted ORCA memory."""
    match = re.fullmatch(r'\s*(\d+(?:\.\d+)?)\s*(MB|GB|TB|MiB|GiB)\s*', str(memory), re.I)
    if not match or float(match[1]) <= 0:
        raise ValueError('mem must be positive, e.g. "4GB" or "2048MB"')
    factors = dict(MB=1, GB=1000, TB=1000000, MIB=2**20/1e6, GIB=2**30/1e6)
    amount = float(match[1]) * factors[match[2].upper()] * 0.75 / nproc
    if amount < 1:
        raise ValueError('mem is too small to assign 1 MB per ORCA process')
    return math.floor(amount)


class OrcaDFTRun:
    def __init__(self, functional, basis, nproc, value, in_file, solvent='0',
                 error=0, pklsave=False, restart=True):
        self.functional = functional
        self.basis = basis
        self.nproc = check_resource.respec_cores(_integer(nproc, 'nproc', 1))
        self.value = value
        self.in_file = str(Path(in_file).expanduser().absolute())
        self.solvent = solvent
        self.error = error
        self.pklsave = pklsave
        self.restart = restart
        self.executable = 'orca'
        self.mem = ''
        self.timexe = 60 * 60 * 80
        self.timejob = None
        self.SpecTotalCharge = np.nan
        self.SpecSpinMulti = np.nan
        self.nstates = 20
        self.maxiter = 200
        self.conformer_config = None
        self.para_functional = []
        # Optional element -> shielding (ppm) at the same method/basis/solvent.
        # Missing H/C/Si references are calculated from optimized TMS with ORCA.
        self.nmr_reference = {}
        self.ion_multiplicities = {}
        self.reuse = True

    def _options(self):
        if str(self.functional).upper().startswith('KTLC-') or self.para_functional:
            raise ValueError('KTLC and custom functional parameters are not supported by orca_run')
        # Shared search syntax defaults to PM6; normalize ORCA's xTB default
        # before passing the string to either shared parser.
        search_value = re.sub(r'(?i)(optconf\s*=\s*\()([^()]*)\)',
                              lambda m: m[0] if re.search(r'conf_relax\s*=', m[2], re.I)
                              else m[1] + m[2] + (', ' if m[2].strip() else '') + 'conf_relax=xtb)',
                              self.value)
        search_value = re.sub(r'(?i)(?<!\S)optconf(?=\s|$)', 'optconf=xtb', search_value)
        options, references = {}, {}
        supported = {'energy', 'opt', 'homolumo', 'dipole', 'cden', 'vip', 'vea',
                     'aip', 'aea', 'uv', 'fluor', 'freq', 'polar', 'nmr', 'optconf'}
        for token in calculation_tokens(search_value):
            key, separator, argument = token.partition('=')
            key = key.lower()
            if key not in supported:
                raise ValueError(f'Unsupported ORCA option: {token}')
            if key in options:
                raise ValueError(f'Duplicate ORCA option: {key}')
            if separator and not argument:
                raise ValueError(f'Missing argument for {key}')
            if key == 'optconf':
                if separator and argument.lower() != 'xtb':
                    raise ValueError('ORCA optconf supports conf_relax=xtb only')
                options[key] = 'xtb'
            elif key == 'opt':
                level, fc = ground_optimization(argument if separator else '', 'orca')
                options[key] = level
                options['opt_fc'] = fc
            elif key == 'fluor':
                spec = parse_excited_optimization(token, 'orca')
                options[key] = spec['state']
                options['fluor_fc'] = spec['fc']
            elif key in ('uv', 'nmr', 'freq'):
                options[key] = True
                if separator:
                    paths = argument.split(',')
                    labels = ('IR', 'Raman', 'NMR') if key == 'freq' else (key.upper(),)
                    if not 1 <= len(paths) <= len(labels) or not all(paths):
                        raise ValueError(f'Invalid reference spectrum list for {key}')
                    references.update({label: str(Path(path).expanduser().resolve())
                                       for label, path in zip(labels, paths)})
            else:
                if separator:
                    raise ValueError(f'{key} does not accept an argument')
                options[key] = True
        if not options:
            raise ValueError('Specify at least one ORCA calculation option')
        if 'NMR' in references:
            options['nmr'] = True
        self._search_config = parse_conformer_options(search_value, self.conformer_config)
        # Parse normalized tokens too: grouped optconf without conf_relax defaults
        # to PM6 in shared settings, while ORCA explicitly defaults to xTB.
        for token in calculation_tokens(search_value):
            if token.lower().startswith('optconf=') and token.split('=', 1)[1].lower() != 'xtb':
                raise ValueError('ORCA optconf supports conf_relax=xtb only')
        self.ref_spectrum_paths = references
        return options

    def _read_structure(self):
        path = Path(self.in_file)
        if path.suffix.lower() == '.sdf':
            mol = next((m for m in Chem.SDMolSupplier(str(path), removeHs=False) if m is not None), None)
            if mol is None or mol.GetNumConformers() == 0:
                raise ValueError('No readable SDF structure')
            atoms = [a.GetSymbol() for a in mol.GetAtoms()]
            xyz = mol.GetConformer().GetPositions()
            charge = Chem.GetFormalCharge(mol)
            mult = Descriptors.NumRadicalElectrons(mol) + 1
        elif path.suffix.lower() in ('.fchk', '.fch'):
            atoms, x, y, z, charge, mult = read_mol_file.read_fchk(path)
            xyz = np.column_stack((x, y, z))
        elif path.suffix.lower() == '.xyz':
            lines = path.read_text().splitlines()
            if len(lines) < 3:
                raise ValueError('Incomplete XYZ input')
            count = _integer(lines[0].strip(), 'XYZ atom count', 1)
            rows = [line.split() for line in lines[2:2 + count]]
            if len(rows) != count or any(len(row) != 4 for row in rows):
                raise ValueError('XYZ atom count or coordinate columns are invalid')
            atoms = [row[0] for row in rows]
            xyz = np.array([[float(v) for v in row[1:]] for row in rows])
            state = re.fullmatch(r'\s*([+-]?\d+)\s+(\d+)\s*', lines[1])
            charge, mult = (int(state[1]), int(state[2])) if state else (0, 1)
        else:
            raise ValueError('ORCA accepts SDF, XYZ or Gaussian FCHK/FCH input')
        if not np.isnan(self.SpecTotalCharge):
            charge = self.SpecTotalCharge
        if not np.isnan(self.SpecSpinMulti):
            mult = self.SpecSpinMulti
        charge, mult = _integer(charge, 'charge'), _integer(mult, 'multiplicity', 1)
        electrons = sum(Chem.GetPeriodicTable().GetAtomicNumber(a) for a in atoms) - charge
        if electrons < mult - 1 or electrons <= 0 or (electrons - mult + 1) % 2:
            raise ValueError('Inconsistent electron count, charge and spin multiplicity')
        if not np.all(np.isfinite(xyz)):
            raise ValueError('Nonfinite atomic coordinates')
        return atoms, xyz, charge, mult

    def make_input(self, jobname, atoms, xyz, charge, multiplicity, *, optimize=None,
                   td=False, root=None, freq=False, polar=False, nmr=False,
                   dipole=False, orbitals=False, moread=None, fc=''):
        if fc and not optimize:
            raise ValueError('FC requires an optimization')
        functional = _token(self.functional, 'functional')
        # Match the Gaussian B3LYP convention; ORCA-native B3LYP is VWN5.
        if functional.upper() == 'B3LYP':
            functional = 'B3LYP/G'
        # ORCA's native B88 derivatives do not cover TDDFT optimization.
        # Use the documented equivalent LibXC form throughout the workflow.
        if functional.upper() in ('B3LYP/G', 'BLYP', 'LC-BLYP', 'CAM-B3LYP'):
            functional = f'LibXC({functional})'
        basis = _token(self.basis, 'basis')
        keywords = [functional, basis, 'TightSCF', 'NoAutoStart']
        if optimize:
            keywords.append({'normal': 'Opt', 'loose': 'LooseOpt',
                             'tight': 'TightOpt', 'verytight': 'VeryTightOpt'}[optimize])
        if freq:
            keywords.append('NumFreq')
        if nmr:
            keywords.append('NMR')
        if moread:
            keywords.append('MOREAD')
        solvent = str(self.solvent).strip() if self.solvent is not None else '0'
        if solvent.lower() not in ('0', 'none', ''):
            keywords.append('CPCM(' + _token(solvent, 'solvent') + ')')
        lines = ['! ' + ' '.join(keywords),
                 f'%pal nprocs {self.nproc} end',
                 f'%maxcore {memory_per_core(self.mem, self.nproc)}',
                 '%scf MaxIter 512 end']
        if moread:
            lines.append(f'%moinp "{Path(moread).name}"')
        if optimize:
            lines.extend(['%geom', f' MaxIter {_integer(self.maxiter, "maxiter", 1)}'])
            if fc == 'ReadFC':
                lines.extend([' InHess Read', f' InHessName "{jobname}.readfc.hess"'])
            elif fc:
                lines.extend([' Calc_Hess true', ' NumHess true'])
                if fc == 'CalcAll' or fc.startswith('RecalcFC='):
                    interval = 1 if fc == 'CalcAll' else int(fc.split('=')[1])
                    lines.append(f' Recalc_Hess {interval}')
                elif fc != 'CalcFC':
                    raise ValueError(f'Unsupported ORCA FC option: {fc}')
            lines.append('end')
        if orbitals:
            # Print every orbital energy, without the much larger coefficients.
            lines.append('%output Print[P_OrbEn] 2 end')
        if td:
            nroots = max(_integer(self.nstates, 'nstates', 1), root or 1)
            lines.extend(['%tddft', f' NRoots {nroots}', ' TDA false'])
            if root is not None:
                lines.extend([f' IRoot {root}', ' FollowIRoot true'])
            lines.append('end')
        if dipole or polar or freq:
            lines.extend(['%elprop', ' Origin CenterOfMass', ' Dipole true'])
            if dipole:
                lines.append(' Quadrupole true')
            if polar or freq:
                lines.append(' Polar true')
            lines.append('end')
        lines.append(f'* xyz {charge} {multiplicity}')
        lines.extend(f'{atom} {x:.12f} {y:.12f} {z:.12f}'
                     for atom, (x, y, z) in zip(atoms, xyz))
        lines.append('*')
        # The EPRNMR parser requires coordinates to precede nucleus selections.
        if nmr:
            lines.extend(['%eprnmr', ' Ori GIAO'])
            lines.extend(f' Nuclei = all {element} {{ shift }}' for element in sorted(set(atoms)))
            lines.append('end')
        Path(f'{jobname}.inp').write_text('\n'.join(lines) + '\n')
        if fc == 'ReadFC':
            source = calculation_cache.previous_hessian('orca', jobname, '\n'.join(lines), self.executable)
            shutil.copyfile(source, f'{jobname}.readfc.hess')

    def _calculate(self, name, atoms, xyz, charge, mult, **kwargs):
        self._active_job = name
        self.make_input(name, atoms, xyz, charge, mult, **kwargs)
        calculation_cache.execute('orca', name, Path(f'{name}.inp').read_text(),
                                  lambda: Exe_Orca.exe_Orca(name, self.timexe, self.executable),
                                  self.executable)
        parsed = parse_log(f'{name}.log')
        parsed.check(optimize=bool(kwargs.get('optimize')))
        return parsed

    def _ion_mult(self, key, mult):
        # Like the existing workflows: remove an unpaired electron for open
        # shells, or form a doublet from a singlet. Users can select alternatives.
        return _integer(self.ion_multiplicities.get(key, 2 if mult == 1 else mult - 1),
                        f'{key} multiplicity', 1)

    def _nmr(self, parsed, output):
        rows = parsed.shielding()
        if ([r[0] for r in rows] != list(range(len(self._atoms)))
                or [r[1] for r in rows] != self._atoms):
            raise ValueError('NMR shielding summary does not match the input atom order')
        references = {element: float(value) for element, value in self.nmr_reference.items()}
        if not all(math.isfinite(value) for value in references.values()):
            raise ValueError('NMR reference shieldings must be finite')
        needed = {element for _, element, _ in rows} & {'H', 'C', 'Si'}
        if needed - references.keys():
            tms = Chem.AddHs(Chem.MolFromSmiles('C[Si](C)(C)C'))
            if AllChem.EmbedMolecule(tms, randomSeed=42) != 0:
                raise RuntimeError('Could not prepare the TMS NMR reference')
            AllChem.UFFOptimizeMolecule(tms)
            atoms = [atom.GetSymbol() for atom in tms.GetAtoms()]
            tms_log = self._calculate(self._jobname + '_TMS', atoms,
                                      tms.GetConformer().GetPositions(), 0, 1,
                                      optimize='tight', nmr=True)
            tms_rows = tms_log.shielding()
            for element in needed - references.keys():
                values = [v for _, e, v in tms_rows if e == element]
                if not values:
                    raise ValueError(f'Missing {element} shielding in the TMS reference')
                references[element] = float(np.mean(values))
        output['nmr_shielding'] = [[e for _, e, _ in rows], [v for _, _, v in rows]]
        output['nmr'] = [[e for _, e, _ in rows],
                         [references[e] - v if e in needed else v for _, e, v in rows]]
        output['nmr_reference'] = {e: references[e] for e in sorted(needed)}
        output['nmr_info'] = {'unit': 'ppm', 'chemical_shift_elements': sorted(needed),
                              'other_elements': 'absolute shielding', 'reference': 'TMS'}
        self._compare(output, 'NMR', [v for e, v in zip(*output['nmr']) if e in needed],
                      [1.0] * sum(e in needed for _, e, _ in rows))

    def _compare(self, output, kind, positions, intensities):
        path = self.ref_spectrum_paths.get(kind)
        if path:
            output.setdefault('Spectrum_similarity', {})[kind] = Spectrum_similarity.compare_with_file(
                path, positions, intensities, kind)

    def _run_orca(self, options, output):
        atoms, xyz, charge, mult = self._read_structure()
        self._atoms = atoms
        self.mem = check_resource.respec_memory(
            spec_memory=self.mem,
            atoms=atoms,
            basis=self.basis,
            total_charge=charge,
            nproc=self.nproc,
            options=options,
            calculate_hessian=requires_hessian(options.get('opt_fc', ''),
                                                options.get('fluor_fc', '')),
        )
        memory_per_core(self.mem, self.nproc)
        jobname = Path(self.in_file).stem
        self._jobname = jobname
        workdir = Path.cwd() / jobname
        calculation_cache.prepare_directory(workdir, self.in_file)
        self._active_directory = workdir
        with chdir(workdir):
            if 'optconf' in options:
                from qcforever.conformer_search.conformer_search import run_conformer_search, read_selected_structure
                try:
                    search_input = self.in_file
                    if Path(self.in_file).suffix.lower() in ('.fchk', '.fch'):
                        search_input = str(Path('fchk_geometry.xyz').resolve())
                        read_mol_file.write_xyz(search_input, atoms, *xyz.T, charge, mult)
                    summary = run_conformer_search(search_input, charge, mult, 'xtb',
                                                    self.nproc, self.mem, self._search_config)
                    output['conformer_search'] = summary
                    selected, x, y, z, *_ = read_selected_structure(summary)
                    if selected != atoms:
                        raise ValueError('Selected conformer atom order/composition changed')
                    xyz = np.column_stack((x, y, z))
                    output['optconf'] = True
                except job_timeout.QCforeverTimeoutError:
                    raise
                except Exception as exc:
                    output['optconf'] = False
                    output['conformer_search'] = {'state': 'failed', 'error': str(exc)}
                    raise RuntimeError('xTB conformer search failed: ' + str(exc)) from exc
            ground = self._calculate(jobname, atoms, xyz, charge, mult,
                                      optimize=options.get('opt'), fc=options.get('opt_fc', ''))
            if 'opt' in options:
                final_atoms, xyz = ground.coordinates()
                if final_atoms != atoms:
                    raise ValueError('Optimized geometry changed atom order/composition')
            self._write_xyz(jobname, atoms, xyz, charge, mult)
            output['energy'] = ground.energy()
            output['Energy'] = output['energy'].copy()
            output['charge'] = charge
            output['multiplicity'] = mult
            electronic = ground
            if 'homolumo' in options or 'dipole' in options:
                electronic = self._calculate(jobname + '_properties', atoms, xyz, charge, mult,
                                              dipole=True, orbitals=True)
            if 'homolumo' in options:
                output.update(electronic.orbitals())
            if 'dipole' in options:
                output['dipole'] = electronic.dipole()
                output.update(electronic.quadrupole())
            if 'cden' in options:
                output['cden'] = ground.charges()
            if 'polar' in options:
                polar_log = self._calculate(jobname + '_polar', atoms, xyz, charge, mult, polar=True)
                output.update(polar_log.polar())
            for vertical, relaxed, delta in (('vip', 'aip', 1), ('vea', 'aea', -1)):
                if vertical in options or relaxed in options:
                    ion_mult = self._ion_mult(vertical, mult)
                    ion = self._calculate(jobname + '_' + vertical.upper(), atoms, xyz, charge + delta, ion_mult)
                    value = delta * (ion.energy()[0] - output['energy'][0]) * EH_TO_EV
                    output[vertical] = [value, ion.energy()[1]]
                if relaxed in options:
                    relaxed_log = self._calculate(jobname + '_' + relaxed.upper(), atoms, xyz,
                                                  charge + delta, _integer(self.ion_multiplicities.get(relaxed, ion_mult),
                                                                           f'{relaxed} multiplicity', 1),
                                                  optimize='normal')
                    value = delta * (relaxed_log.energy()[0] - output['energy'][0]) * EH_TO_EV
                    output[relaxed] = [value, relaxed_log.energy()[1]]
            if 'uv' in options or 'fluor' in options:
                td = self._calculate(jobname + '_TD', atoms, xyz, charge, mult, td=True)
                transitions = td.transitions()
                output['uv'] = [[v for _, v, _ in transitions], [v for _, _, v in transitions]]
                if 'UV' in self.ref_spectrum_paths:
                    ref = Spectrum_similarity.read_data(self.ref_spectrum_paths['UV'])
                    s, d, w = UV_similarity.smililarity_dissimilarity(*ref, *output['uv'])
                    output['Similarity/Dissimilarity'] = [s, d]
                    output['Wasserstein_distance'] = [w]
            if 'fluor' in options:
                root = options['fluor']
                if root not in {i for i, _, _ in transitions}:
                    raise ValueError(f'fluor={root} is outside the computed excited states')
                excited = self._calculate(jobname + '_TDopt', atoms, xyz, charge, mult,
                                          td=True, root=root, optimize='normal',
                                          fc=options.get('fluor_fc', ''))
                ex_atoms, ex_xyz = excited.coordinates()
                self._write_xyz(jobname + '_TDopt', ex_atoms, ex_xyz, charge, mult)
                followed = re.findall(r'The IROOT now is\s*:\s*(\d+)', excited.text, re.I)
                root = int(followed[-1]) if followed else root
                emission = self._calculate(jobname + '_Emission', ex_atoms, ex_xyz, charge, mult, td=True)
                transitions = emission.transitions()
                selected = next((row for row in transitions if row[0] == root), None)
                if selected is None:
                    raise ValueError('The relaxed emission state is absent from the TDDFT spectrum')
                output['fluor'] = [[v for _, v, _ in transitions], [v for _, _, v in transitions]]
                output['fluor_state'] = {'root': root, 'wavelength_nm': selected[1],
                                          'oscillator_strength': selected[2]}
            if 'freq' in options:
                vibration = self._calculate(jobname + '_freq', atoms, xyz, charge, mult, freq=True)
                output.update(vibration.frequencies())
                self._compare(output, 'IR', output['freq'], output['IR'])
                self._compare(output, 'Raman', output['freq'], output['Raman'])
            if 'nmr' in options:
                nmr_log = self._calculate(jobname + '_nmr', atoms, xyz, charge, mult, nmr=True)
                self._nmr(nmr_log, output)
            output['log'] = 'normal'

    @staticmethod
    def _write_xyz(name, atoms, xyz, charge, mult):
        Path(f'{name}.xyz').write_text(f'{len(atoms)}\n{charge} {mult}\n' + ''.join(
            f'{atom} {x:.12f} {y:.12f} {z:.12f}\n' for atom, (x, y, z) in zip(atoms, xyz)))

    @calculation_cache.workflow('orca')
    def run_orca(self):
        """Return properties or log=error/timeout; retain inputs and diagnostics."""
        self._active_directory = None
        self._active_job = None
        output = {}
        try:
            with job_timeout.overall_timeout(self.timejob):
                options = self._options()
                self._run_orca(options, output)
        except job_timeout.QCforeverTimeoutError as exc:
            output.update(log='timeout', error=str(exc))
        except (OSError, ValueError, RuntimeError) as exc:
            output.update(log='error', error=str(exc))
            if isinstance(exc, Exe_Orca.OrcaExecutionError):
                output.update(returncode=exc.returncode, logfile=exc.logfile)
        finally:
            if self._active_job:
                output['jobname'] = self._active_job
                output.setdefault('logfile', str(self._active_directory / f'{self._active_job}.log'))
            if self._active_directory:
                if self.pklsave:
                    with (self._active_directory / (self._jobname + '.pkl')).open('wb') as stream:
                        pickle.dump(output, stream)
                job_cleanup.cleanup_orca(self._active_directory, preserve_pickle=self.pklsave,
                                         preserve_restart=self.restart,
                                         preserve_conformers=self._search_config is not None)
        return output
