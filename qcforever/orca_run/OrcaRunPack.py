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
from tempfile import TemporaryDirectory

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from qcforever.conformer_search.settings import calculation_tokens, parse_conformer_options
from qcforever.util import check_resource, job_cleanup, job_timeout, Spectrum_similarity, UV_similarity
from . import Exe_Orca
from .parse_log import EH_TO_EV, parse_log
from qcforever.util import calculation_cache
from qcforever.util import read_mol_file, fchk_orbitals
from qcforever.util.calculation_options import (
    cden_fukui, ground_optimization, parse_excited_optimization, requires_hessian,
)
from qcforever.util.fukui import condensed_fukui
from qcforever.util.imaginary_modes import NO_STABLE_MINIMUM, StationaryPoint, follow_imaginary_modes


_KTLC_FUNCTIONALS = {'KTLC-BLYP-BO': 'LC-BLYP', 'KTLC-PBE-BO': 'LC-PBE'}


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
        self._ktlc_mu = None
        # Optional element -> shielding (ppm) at the same method/basis/solvent.
        # Missing H/C/Si references are calculated from optimized TMS with ORCA.
        self.nmr_reference = {}
        self.ion_multiplicities = {}
        self.reuse = True
        self.fchk_angular = 'auto'
        self._fchk_orbitals = None

    def _options(self):
        functional = str(self.functional).strip().upper()
        if functional.startswith('KTLC-') and functional not in _KTLC_FUNCTIONALS:
            raise ValueError(f'Unsupported ORCA KTLC functional: {self.functional}')
        if self.para_functional:
            raise ValueError('Manual functional parameter overrides are not supported by orca_run')
        # Shared search syntax defaults to PM6; normalize ORCA's xTB default
        # before passing the string to either shared parser.
        search_value = re.sub(r'(?i)(optconf\s*=\s*\()([^()]*)\)',
                              lambda m: m[0] if re.search(r'conf_relax\s*=', m[2], re.I)
                              else m[1] + m[2] + (', ' if m[2].strip() else '') + 'conf_relax=xtb)',
                              self.value)
        search_value = re.sub(r'(?i)(?<!\S)optconf(?=\s|$)', 'optconf=xtb', search_value)
        options, references = {}, {}
        supported = {'energy', 'opt', 'homolumo', 'dipole', 'cden', 'vip', 'vea',
                     'aip', 'aea', 'uv', 'fluor', 'freq', 'polar', 'nmr', 'optconf', 'stable'}
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
            elif key == 'cden':
                options[key] = True
                if cden_fukui(token):
                    options['fukui'] = True
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
        if 'stable' in options:
            options.setdefault('opt', 'tight')
        if 'fukui' in options:
            options.update(vip=True, vea=True)
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
                   dipole=False, orbitals=False, moread=None, fc='', range_mu=None, mode_scan=False):
        if fc and not optimize:
            raise ValueError('FC requires an optimization')
        functional = _token(self.functional, 'functional')
        ktlc_base = _KTLC_FUNCTIONALS.get(functional.upper())
        if ktlc_base:
            functional = ktlc_base
            range_mu = self._ktlc_mu if range_mu is None else range_mu
            if range_mu is None or not math.isfinite(range_mu) or range_mu <= 0:
                raise ValueError('ORCA KTLC requires a positive, finite range-separation parameter')
        # Match the Gaussian B3LYP convention; ORCA-native B3LYP is VWN5.
        if functional.upper() == 'B3LYP':
            functional = 'B3LYP/G'
        # ORCA's native B88 derivatives do not cover TDDFT optimization.
        # Use the documented equivalent LibXC form throughout the workflow.
        if functional.upper() in ('B3LYP/G', 'BLYP', 'LC-BLYP', 'CAM-B3LYP'):
            functional = f'LibXC({functional})'
        imported = getattr(self, '_fchk_orbitals', None)
        basis = [] if imported else [_token(self.basis, 'basis')]
        keywords = [functional, *basis, 'VeryTightSCF' if mode_scan else 'TightSCF', 'NoAutoStart']
        if imported:
            keywords.append('NoECP')
            if imported.beta is not None:
                keywords.append('UKS')
        if optimize:
            keywords.append({'normal': 'Opt', 'loose': 'LooseOpt',
                             'tight': 'TightOpt', 'verytight': 'VeryTightOpt'}[optimize])
            if td and ktlc_base == 'LC-PBE':
                # Native LC-PBE lacks the XC third derivative used by analytic
                # excited-state gradients in ORCA 6.1.
                keywords.append('NumGrad')
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
        if ktlc_base:
            mu = repr(float(range_mu))
            lines.extend(['%method', f' RangeSepMu {mu}'])
            if ktlc_base == 'LC-BLYP':
                # LibXC evaluates the short-range exchange with its own omega.
                lines.append(f' ExtParamXC "_omega" {mu}')
            lines.append('end')
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
        if imported and list(atoms) != imported.atoms:
            raise ValueError('FCHK basis requires the source atom order/composition')
        for i, (atom, (x, y, z)) in enumerate(zip(atoms, xyz)):
            lines.append(f'{atom} {x:.12f} {y:.12f} {z:.12f}')
            if imported:
                lines.extend(imported.orca_basis(i))
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

    def _calculate(self, name, atoms, xyz, charge, mult, *, moread=None, require_moread=False,
                   require_gbw=False, **kwargs):
        self._active_job = name
        source = Path(moread).name if moread else None
        if require_moread and (source is None or source == f'{name}.gbw' or not Path(source).is_file()):
            raise ValueError(f'Missing ORCA orbital guess for {name}: {source}')
        if source == f'{name}.gbw' or (source and not Path(source).is_file()):
            source = None
        required = {'required_outputs': (f'{name}.gbw',)} if require_gbw else {}
        for guess in (source, None):
            self.make_input(name, atoms, xyz, charge, mult, moread=guess, **kwargs)
            try:
                calculation_cache.execute('orca', name, Path(f'{name}.inp').read_text(),
                                          lambda: Exe_Orca.exe_Orca(name, self.timexe, self.executable),
                                          self.executable, **required)
            except Exe_Orca.OrcaExecutionError as exc:
                if require_moread or guess is None or not re.search(
                        r'error termination in GUESS\b|no orbitals were found in the gbw file|orca_guess/guess_restart\.cpp',
                        str(exc), re.I):
                    raise
                Path(exc.logfile).replace(f'{name}.moread.log')
                Path(f'{name}.gbw').unlink(missing_ok=True)
                continue
            parsed = parse_log(f'{name}.log')
            parsed.check(optimize=bool(kwargs.get('optimize')))
            if require_gbw and not Path(f'{name}.gbw').is_file():
                raise ValueError(f'ORCA did not write orbitals: {name}.gbw')
            return parsed

    def _ion_mult(self, key, mult):
        # Like the existing workflows: remove an unpaired electron for open
        # shells, or form a doublet from a singlet. Users can select alternatives.
        return _integer(self.ion_multiplicities.get(key, 2 if mult == 1 else mult - 1),
                        f'{key} multiplicity', 1)

    def _tune_koopmans(self, name, atoms, xyz, charge, mult):
        """Minimize the vertical IP minus the negative HOMO at fixed geometry."""
        from bayes_opt import BayesianOptimization, UtilityFunction

        tolerance = 0.01
        limit = 52
        trials = []
        best = None
        ion_mult = self._ion_mult('vip', mult)

        def objective(mu):
            nonlocal best
            trial_name = f'{name}_KTLC_{len(trials) + 1:03d}'
            trial = {'mu': float(mu), 'jobname': trial_name}
            trials.append(trial)
            try:
                parent = self._calculate(trial_name, atoms, xyz, charge, mult,
                                         orbitals=True, range_mu=mu)
                homo_ev = parent.homo() * EH_TO_EV
                cation = self._calculate(trial_name + '_VIP', atoms, xyz, charge + 1,
                                         ion_mult, range_mu=mu, moread=f'{trial_name}.gbw')
                vip_ev = (cation.energy()[0] - parent.energy()[0]) * EH_TO_EV
                residual = vip_ev + homo_ev
                score = -(residual ** 2)
                if not math.isfinite(score):
                    raise ValueError('Non-finite Koopmans residual')
            except (Exe_Orca.OrcaExecutionError, ValueError) as exc:
                trial['error'] = str(exc)
                return -10.0
            trial.update(vip_ev=vip_ev, homo_ev=homo_ev, error_ev=abs(residual))
            if best is None or trial['error_ev'] < best['error_ev']:
                best = trial
            return score

        optimizer = BayesianOptimization(f=objective, pbounds={'mu': (0.01, 1.0)},
                                          verbose=0, random_state=1)
        for initial in (0.15, 0.35, 0.65):
            optimizer.probe(params={'mu': initial}, lazy=True)
        acquisition = UtilityFunction(kind='ei', xi=1e-4)
        optimizer.maximize(init_points=0, n_iter=0, acquisition_function=acquisition)
        for _ in range(limit - len(trials)):
            if best is not None and best['error_ev'] <= tolerance:
                break
            optimizer.maximize(init_points=0, n_iter=1, acquisition_function=acquisition)

        if best is None:
            state, reason = 'failed', 'no_valid_evaluations'
        elif best['error_ev'] <= tolerance:
            state, reason = 'converged', 'tolerance_reached'
        else:
            state, reason = 'best_evaluated', 'evaluation_limit'
        return {'state': state, 'selection_reason': reason,
                'selected_mu': best['mu'] if best else None,
                'error_ev': best['error_ev'] if best else None,
                'tolerance_ev': tolerance, 'valid_evaluations': sum('error_ev' in t for t in trials),
                'trials': trials}

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
        self._fchk_orbitals = fchk_orbitals.for_job(self, 'orca', options)
        self._atoms = atoms
        self.mem = check_resource.respec_memory(
            spec_memory=self.mem,
            atoms=atoms,
            basis=self.basis,
            total_charge=charge,
            nproc=self.nproc,
            options=options,
            calculate_hessian='stable' in options or requires_hessian(options.get('opt_fc', ''),
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
            if str(self.functional).strip().upper() in _KTLC_FUNCTIONALS:
                tuning = {'state': 'failed', 'selection_reason': 'execution_error',
                          'selected_mu': None, 'error_ev': None}
                output['functional_parameter_optimization'] = tuning
                try:
                    # A separate working directory keeps the trial files out
                    # of the persistent calculation cache.
                    with TemporaryDirectory(prefix='ktlc-', dir='.') as scratch, chdir(scratch):
                        tuning.update(self._tune_koopmans(jobname, atoms, xyz, charge, mult))
                except job_timeout.QCforeverTimeoutError as exc:
                    tuning.update(state='timeout', selection_reason='time_limit')
                    raise job_timeout.QCforeverTimeoutError('ORCA KTLC parameter optimization timed out') from exc
                finally:
                    self._active_job = None
                self._ktlc_mu = tuning['selected_mu']
                if self._ktlc_mu is None:
                    raise RuntimeError('ORCA KTLC found no valid Koopmans evaluation')
                output['functional_param'] = [self._ktlc_mu]
            guess = self._fchk_orbitals.orca_gbw(self.executable, jobname+'_fchk_guess') if self._fchk_orbitals else None
            ground = self._calculate(jobname, atoms, xyz, charge, mult, moread=guess,
                                      require_moread=self._fchk_orbitals is not None,
                                      require_gbw='stable' in options,
                                      optimize=options.get('opt'), fc=options.get('opt_fc', ''))
            ground_gbw = f'{jobname}.gbw'
            if 'opt' in options:
                final_atoms, xyz = ground.coordinates()
                if final_atoms != atoms:
                    raise ValueError('Optimized geometry changed atom order/composition')
            stable_vibration = None
            if 'stable' in options:
                def stationary(name, elements, coords, optimized):
                    vibr = self._calculate(name + '_freq', elements, coords, charge, mult, freq=True,
                                           moread=name + '.gbw', require_moread=True).frequencies()
                    return StationaryPoint(name, elements, coords, optimized.energy()[0], vibr['freq'],
                                           vibr['freqmode'], {'ground': optimized, 'vibration': vibr})

                def single_point(label, center, coordinates, reference_label):
                    name = f'{jobname}_stable_{label}'
                    guess = f'{jobname}_stable_{reference_label}' if reference_label else center.name
                    return self._calculate(name, center.atoms, coordinates, charge, mult, mode_scan=True,
                                           moread=guess + '.gbw', require_moread=True, require_gbw=True).energy()[0]

                def evaluate(label, elements, coordinates, guess_label):
                    name = f'{jobname}_stable_{label}'
                    optimized = self._calculate(name, elements, coordinates, charge, mult, optimize=options['opt'],
                                                moread=f'{jobname}_stable_{guess_label}.gbw',
                                                require_moread=True, require_gbw=True)
                    final_elements, final_coords = optimized.coordinates()
                    if final_elements != elements:
                        raise ValueError('Optimized geometry changed atom order/composition')
                    return stationary(name, final_elements, final_coords, optimized)

                point = follow_imaginary_modes(stationary(jobname, atoms, xyz, ground), single_point, evaluate)
                output['stable_0'] = point.minimum_frequency >= 0
                if not output['stable_0']:
                    output.update(point.properties['vibration'])
                    self._active_job = point.name + '_freq'
                    raise ValueError(point.properties.get('stable_error', NO_STABLE_MINIMUM))
                xyz, ground = point.coordinates, point.properties['ground']
                if point.name != jobname:
                    for suffix in ('.inp', '.log', '.gbw'):
                        path = Path(point.name + suffix)
                        if path.is_file():
                            shutil.copyfile(path, jobname + suffix)
                self._active_job = point.name + '_freq'
                stable_vibration = point.properties['vibration']
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
            ion_charges = {}
            for vertical, relaxed, delta in (('vip', 'aip', 1), ('vea', 'aea', -1)):
                if vertical in options or relaxed in options:
                    ion_mult = self._ion_mult(vertical, mult)
                    ion_name = jobname + '_' + vertical.upper()
                    ion = self._calculate(ion_name, atoms, xyz, charge + delta, ion_mult,
                                          moread=ground_gbw)
                    value = delta * (ion.energy()[0] - output['energy'][0]) * EH_TO_EV
                    output[vertical] = [value, ion.energy()[1]]
                    if 'fukui' in options:
                        ion_charges[vertical] = ion.charges()
                if relaxed in options:
                    relaxed_log = self._calculate(jobname + '_' + relaxed.upper(), atoms, xyz,
                                                  charge + delta, _integer(self.ion_multiplicities.get(relaxed, ion_mult),
                                                                           f'{relaxed} multiplicity', 1),
                                                  optimize='normal', moread=f'{ion_name}.gbw')
                    value = delta * (relaxed_log.energy()[0] - output['energy'][0]) * EH_TO_EV
                    output[relaxed] = [value, relaxed_log.energy()[1]]
            if 'fukui' in options:
                output['fukui'] = condensed_fukui(
                    output['cden'], ion_charges['vip'], ion_charges['vea'])
            if 'uv' in options or 'fluor' in options:
                td = self._calculate(jobname + '_TD', atoms, xyz, charge, mult,
                                     td=True, moread=ground_gbw)
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
                excited_name = jobname + '_TDopt'
                excited = self._calculate(excited_name, atoms, xyz, charge, mult,
                                          td=True, root=root, optimize='normal',
                                          fc=options.get('fluor_fc', ''), moread=ground_gbw)
                ex_atoms, ex_xyz = excited.coordinates()
                self._write_xyz(excited_name, ex_atoms, ex_xyz, charge, mult)
                followed = re.findall(r'The IROOT now is\s*:\s*(\d+)', excited.text, re.I)
                root = int(followed[-1]) if followed else root
                emission = self._calculate(jobname + '_Emission', ex_atoms, ex_xyz, charge, mult,
                                           td=True, moread=f'{excited_name}.gbw')
                transitions = emission.transitions()
                selected = next((row for row in transitions if row[0] == root), None)
                if selected is None:
                    raise ValueError('The relaxed emission state is absent from the TDDFT spectrum')
                output['fluor'] = [[v for _, v, _ in transitions], [v for _, _, v in transitions]]
                output['fluor_state'] = {'root': root, 'wavelength_nm': selected[1],
                                          'oscillator_strength': selected[2]}
            if 'freq' in options or stable_vibration is not None:
                if stable_vibration is None:
                    vibration = self._calculate(jobname + '_freq', atoms, xyz, charge, mult, freq=True)
                    stable_vibration = vibration.frequencies()
                output.update(stable_vibration)
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
        self._ktlc_mu = None
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
