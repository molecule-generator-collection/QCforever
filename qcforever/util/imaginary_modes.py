"""Follow a ground-state imaginary mode using Cartesian coordinates in Å."""

from dataclasses import dataclass

import numpy as np

from .job_timeout import QCforeverTimeoutError


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


def displace_imaginary_mode(coordinates, frequencies, modes, distance=0.05):
    """Return (+, −) geometries with a maximum per-atom displacement in Å."""
    xyz = np.asarray(coordinates, dtype=float)
    frequencies = np.asarray(frequencies, dtype=float)
    if (xyz.ndim != 2 or xyz.shape[1] != 3 or not np.isfinite(xyz).all()
            or frequencies.ndim != 1 or not np.isfinite(frequencies).all()
            or not np.isfinite(distance) or distance <= 0):
        raise ValueError('Invalid Cartesian coordinates, frequencies or displacement distance')
    if frequencies.min(initial=np.inf) >= 0:
        return None
    axes = np.asarray(modes, dtype=float)
    if axes.shape != (3, len(frequencies), len(xyz)):
        raise ValueError('Normal modes must have shape (3, number of modes, number of atoms)')
    mode = axes[:, frequencies.argmin(), :].T
    maximum = np.linalg.norm(mode, axis=1).max()
    if not np.isfinite(maximum) or maximum == 0:
        raise ValueError('Invalid Cartesian displacement for the imaginary mode')
    step = distance * mode / maximum
    return xyz + step, xyz - step


def follow_imaginary_modes(point, evaluate, max_steps=9):
    """Evaluate both signs at each step and return the best converged point."""
    for iteration in range(1, max_steps + 1):
        pair = displace_imaginary_mode(point.coordinates, point.frequencies, point.modes)
        if pair is None:
            break
        print(f'{point.name}: imaginary frequency {point.minimum_frequency:.2f} cm^-1')
        candidates = []
        for sign, coordinates in zip(('plus', 'minus'), pair):
            label = f'{iteration}_{sign}'
            try:
                candidate = evaluate(label, point.atoms, coordinates)
                if candidate.atoms != point.atoms or not np.isfinite(candidate.energy):
                    raise ValueError('Invalid atom order or energy after optimization')
                candidate.minimum_frequency  # Validate before comparing the candidates.
                candidates.append(candidate)
            except QCforeverTimeoutError:
                raise
            except (OSError, ValueError, RuntimeError) as exc:
                print(f'Imaginary-mode trial {label} failed: {exc}')
        if not candidates:
            break
        point = min(candidates, key=lambda candidate: (candidate.minimum_frequency < 0, candidate.energy))
    return point
