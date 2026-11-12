"""Follow a ground-state imaginary mode using Cartesian coordinates in Å."""

from dataclasses import dataclass

import numpy as np

from .job_timeout import QCforeverTimeoutError


_INITIAL_RMS_DISPLACEMENT = 0.003  # Å across all atoms
_INITIAL_MAX_DISPLACEMENT = 0.015  # Å for any one atom
_MAX_DISPLACEMENT = 0.05
_REPEATED_GEOMETRY_TOLERANCE = 0.001
_MIN_ENERGY_DROP = 1e-8  # Hartree
NO_STABLE_MINIMUM = 'Imaginary-mode search did not find a ground-state minimum'


@dataclass
class StationaryPoint:
    name: str
    atoms: list
    coordinates: np.ndarray
    energy: float
    frequencies: list
    modes: list
    properties: dict

    @property
    def minimum_frequency(self):
        frequencies = np.asarray(self.frequencies, dtype=float)
        if frequencies.ndim != 1 or not np.isfinite(frequencies).all():
            raise ValueError(f'Invalid frequencies for {self.name}')
        return float(frequencies.min(initial=np.inf))


def same_geometry(first, second, tolerance=3e-5):
    """Compare atom-ordered Cartesian coordinates up to translation and rotation."""
    first, second = (np.asarray(xyz, dtype=float) for xyz in (first, second))
    if first.shape != second.shape or first.ndim != 2 or first.shape[1] != 3 or not len(first):
        return False
    if not np.isfinite(first).all() or not np.isfinite(second).all():
        return False
    a, b = first - first.mean(axis=0), second - second.mean(axis=0)
    left, _, right = np.linalg.svd(a.T @ b)
    left[:, -1] *= np.linalg.det(left @ right)
    return np.linalg.norm(a @ left @ right - b, axis=1).max() <= tolerance


def _imaginary_mode(coordinates, frequencies, modes):
    xyz = np.asarray(coordinates, dtype=float)
    frequencies = np.asarray(frequencies, dtype=float)
    if (xyz.ndim != 2 or xyz.shape[1] != 3 or not np.isfinite(xyz).all()
            or frequencies.ndim != 1 or not np.isfinite(frequencies).all()):
        raise ValueError('Invalid Cartesian coordinates or frequencies')
    if frequencies.min(initial=np.inf) >= 0:
        return xyz, None
    axes = np.asarray(modes, dtype=float)
    if axes.shape != (3, len(frequencies), len(xyz)):
        raise ValueError('Normal modes must have shape (3, number of modes, number of atoms)')
    mode = axes[:, frequencies.argmin(), :].T
    maximum = np.linalg.norm(mode, axis=1).max()
    if not np.isfinite(maximum) or maximum == 0:
        raise ValueError('Invalid Cartesian displacement for the imaginary mode')
    return xyz, mode / maximum


def displace_imaginary_mode(coordinates, frequencies, modes, distance=_MAX_DISPLACEMENT):
    """Return (+, −) geometries with a maximum per-atom displacement in Å."""
    if not np.isfinite(distance) or distance <= 0:
        raise ValueError('Invalid displacement distance')
    xyz, mode = _imaginary_mode(coordinates, frequencies, modes)
    if mode is None:
        return None
    step = distance * mode
    return xyz + step, xyz - step


def _scan_mode(single_point, point, xyz, mode, widths, iteration, reference, best):
    reference_label, reference_energy = reference
    for index, width in enumerate(widths, 1):
        for sign, direction in (('plus', 1), ('minus', -1)):
            label = f'{iteration}_sp_{index}_{sign}'
            try:
                energy = single_point(label, point, xyz + direction * width * mode, reference_label)
                if not np.isfinite(energy):
                    raise ValueError('Invalid single-point energy')
            except QCforeverTimeoutError:
                raise
            except (OSError, ValueError, RuntimeError) as exc:
                print(f'Imaginary-mode single point {label} failed: {exc}', flush=True)
                continue
            delta = energy - reference_energy
            print(f'{point.name}: {label}, max displacement {width:.5f} Å per atom, '
                  f'ΔE={delta:+.3e} Eh', flush=True)
            if sign not in best or (energy, width) < best[sign][:2]:
                best[sign] = energy, width, label


def _scan_failure(reference_energy, best):
    directions = []
    for sign in ('plus', 'minus'):
        if sign in best:
            energy, width, _ = best[sign]
            directions.append(f'{sign} ΔE={energy - reference_energy:+.3e} Eh at {width:.5f} Å')
        else:
            directions.append(f'{sign} has no converged single point')
    return (f'Neither direction showed an energy drop greater than {_MIN_ENERGY_DROP:.0e} Eh. '
            f'Lowest sampled values relative to the reference are {", ".join(directions)}.')


def follow_imaginary_modes(point, single_point, evaluate, max_steps=9):
    """Scan both signs at fixed geometries, then optimize the best seed in each downhill direction."""
    visited = []
    distance = None
    for iteration in range(1, max_steps + 1):
        if distance is None:
            xyz, mode = _imaginary_mode(point.coordinates, point.frequencies, point.modes)
            if mode is None:
                break
            visited.append(xyz)
            rms = np.sqrt(np.mean(np.sum(mode**2, axis=1)))
            initial = min(_INITIAL_MAX_DISPLACEMENT, _INITIAL_RMS_DISPLACEMENT / rms)
            widths = (initial / 2, initial, 2 * initial)
            print(f'{point.name}: imaginary frequency {point.minimum_frequency:.2f} cm^-1', flush=True)
            reference_label = f'{iteration}_sp_reference'
            reference_energy = single_point(reference_label, point, xyz, None)
            if not np.isfinite(reference_energy):
                raise ValueError(f'Invalid reference single-point energy for {point.name}')
            reference = reference_label, reference_energy
            best = {}
            optimized = set()
        else:
            widths = (min(2 * distance, _MAX_DISPLACEMENT),)
        distance = widths[-1]
        _scan_mode(single_point, point, xyz, mode, widths, iteration, reference, best)
        downhill = {sign: sample for sign, sample in best.items()
                    if sample[0] - reference_energy < -_MIN_ENERGY_DROP}
        if not downhill:
            if distance >= _MAX_DISPLACEMENT or iteration == max_steps:
                point.properties['stable_error'] = _scan_failure(reference_energy, best)
                print(f'{point.name}: {point.properties["stable_error"]}', flush=True)
                break
            continue

        candidates = []
        for sign, direction in (('plus', 1), ('minus', -1)):
            if sign not in downhill:
                continue
            _, width, guess_label = downhill[sign]
            if (sign, width) in optimized:
                continue
            optimized.add((sign, width))
            label = f'{iteration}_{sign}'
            print(f'{point.name}: optimizing {label} from {width:.5f} Å per atom', flush=True)
            try:
                candidate = evaluate(label, point.atoms, xyz + direction * width * mode, guess_label)
                if candidate.atoms != point.atoms or not np.isfinite(candidate.energy):
                    raise ValueError('Invalid atom order or energy after optimization')
                candidate.minimum_frequency  # Validate before comparing the candidates.
                candidates.append(candidate)
            except QCforeverTimeoutError:
                raise
            except (OSError, ValueError, RuntimeError) as exc:
                print(f'Imaginary-mode trial {label} failed: {exc}')
        minima = [candidate for candidate in candidates if candidate.minimum_frequency >= 0]
        if minima:
            return min(minima, key=lambda candidate: candidate.energy)
        new_points = [candidate for candidate in candidates
                      if not any(same_geometry(candidate.coordinates, previous, _REPEATED_GEOMETRY_TOLERANCE)
                                 for previous in visited)]
        if new_points:
            point = min(new_points, key=lambda candidate: candidate.energy)
            distance = None
        elif distance >= _MAX_DISPLACEMENT:
            break
    return point
