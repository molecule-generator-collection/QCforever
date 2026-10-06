"""Utilities for comparing calculated IR, Raman, and NMR peak spectra."""

import math

import numpy as np
from scipy.stats import wasserstein_distance


_trapezoid = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz


SPECTRUM_SETTINGS = {
    "IR": {"degeneracy_tolerance": 1.0, "half_width": 5.0, "step": 1.0},
    "RAMAN": {"degeneracy_tolerance": 1.0, "half_width": 5.0, "step": 1.0},
    "NMR": {"degeneracy_tolerance": 0.01, "half_width": 0.05, "step": 0.005},
}


def read_data(inputdata):
    """Read a whitespace-separated ``position intensity`` peak file.

    Empty lines and lines beginning with ``#`` are ignored.  A third and later
    column may therefore be used for comments or other metadata.
    """
    positions = []
    intensities = []
    with open(inputdata, "r", encoding="utf-8") as infile:
        for line_number, line in enumerate(infile, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            columns = stripped.split()
            if len(columns) < 2:
                raise ValueError(
                    f"{inputdata}:{line_number}: expected position and intensity"
                )
            try:
                position = float(columns[0])
                intensity = float(columns[1])
            except ValueError as exc:
                raise ValueError(
                    f"{inputdata}:{line_number}: position and intensity must be numbers"
                ) from exc
            if not math.isfinite(position) or not math.isfinite(intensity):
                raise ValueError(f"{inputdata}:{line_number}: values must be finite")
            if intensity < 0:
                raise ValueError(f"{inputdata}:{line_number}: intensity must be non-negative")
            positions.append(position)
            intensities.append(intensity)
    if len(positions) == 0:
        raise ValueError(f"No peaks were found in {inputdata}")
    return [positions, intensities]


def merge_degenerate_peaks(positions, intensities, tolerance):
    """Merge nearby peaks, summing intensity and averaging peak position.

    Peak positions are deliberately averaged without intensity weighting, as
    degenerate transitions describe multiple peaks at equivalent positions.
    """
    if tolerance < 0:
        raise ValueError("tolerance must be non-negative")
    if len(positions) != len(intensities):
        raise ValueError("positions and intensities must have the same length")
    if len(positions) == 0:
        return [], []

    peaks = sorted((float(p), float(i)) for p, i in zip(positions, intensities))
    if any(not math.isfinite(p) or not math.isfinite(i) for p, i in peaks):
        raise ValueError("positions and intensities must be finite")
    if any(i < 0 for _, i in peaks):
        raise ValueError("intensities must be non-negative")
    groups = [[peaks[0]]]
    for position, intensity in peaks[1:]:
        group_mean = sum(item[0] for item in groups[-1]) / len(groups[-1])
        if abs(position - group_mean) <= tolerance or math.isclose(
            abs(position - group_mean), tolerance, rel_tol=1.0e-12, abs_tol=1.0e-15
        ):
            groups[-1].append((position, intensity))
        else:
            groups.append([(position, intensity)])

    merged_positions = [
        sum(item[0] for item in group) / len(group) for group in groups
    ]
    merged_intensities = [sum(item[1] for item in group) for group in groups]
    return merged_positions, merged_intensities


def _lorentzian_grid(positions, intensities, lower, upper, step, half_width):
    grid = np.arange(lower, upper + step * 0.5, step, dtype=float)
    values = np.zeros_like(grid)
    for position, intensity in zip(positions, intensities):
        values += intensity / (
            math.pi * half_width * (1.0 + ((grid - position) / half_width) ** 2)
        )
    return grid, values


def spectrum_similarity(
    reference_positions,
    reference_intensities,
    target_positions,
    target_intensities,
    spectrum_type,
    degeneracy_tolerance=None,
):
    """Return similarity metrics for two discrete peak spectra.

    ``similarity`` is the cosine overlap of Lorentzian-broadened spectra and is
    in the range [0, 1].  ``dissimilarity`` is their integrated squared
    difference.  ``wasserstein_distance`` is computed from normalized discrete
    peak intensities in the native position unit (cm-1 for IR/Raman, ppm for
    NMR).
    """
    kind = spectrum_type.upper()
    if kind not in SPECTRUM_SETTINGS:
        raise ValueError(f"Unsupported spectrum type: {spectrum_type}")
    settings = SPECTRUM_SETTINGS[kind]
    tolerance = (
        settings["degeneracy_tolerance"]
        if degeneracy_tolerance is None
        else float(degeneracy_tolerance)
    )

    ref_positions, ref_intensities = merge_degenerate_peaks(
        reference_positions, reference_intensities, tolerance
    )
    target_positions, target_intensities = merge_degenerate_peaks(
        target_positions, target_intensities, tolerance
    )
    ref_positions, ref_intensities = _active_peaks(ref_positions, ref_intensities)
    target_positions, target_intensities = _active_peaks(
        target_positions, target_intensities
    )

    half_width = settings["half_width"]
    step = settings["step"]
    lower = min(min(ref_positions), min(target_positions)) - 10.0 * half_width
    upper = max(max(ref_positions), max(target_positions)) + 10.0 * half_width
    grid, ref_curve = _lorentzian_grid(
        ref_positions, ref_intensities, lower, upper, step, half_width
    )
    _, target_curve = _lorentzian_grid(
        target_positions, target_intensities, lower, upper, step, half_width
    )

    ref_norm = float(_trapezoid(ref_curve * ref_curve, x=grid))
    target_norm = float(_trapezoid(target_curve * target_curve, x=grid))
    cross = float(_trapezoid(ref_curve * target_curve, x=grid))
    similarity = cross / math.sqrt(ref_norm * target_norm)
    dissimilarity = float(_trapezoid((ref_curve - target_curve) ** 2, x=grid))
    wasserstein = wasserstein_distance(
        ref_positions,
        target_positions,
        np.asarray(ref_intensities) / sum(ref_intensities),
        np.asarray(target_intensities) / sum(target_intensities),
    )
    return {
        "similarity": float(np.clip(similarity, 0.0, 1.0)),
        "dissimilarity": dissimilarity,
        "wasserstein_distance": float(wasserstein),
    }


def compare_with_file(reference_file, positions, intensities, spectrum_type):
    """Read a reference peak file and compare it with a calculated spectrum."""
    reference_positions, reference_intensities = read_data(reference_file)
    return spectrum_similarity(
        reference_positions,
        reference_intensities,
        positions,
        intensities,
        spectrum_type,
    )


def _active_peaks(positions, intensities):
    active = [(p, i) for p, i in zip(positions, intensities) if i > 0]
    if not active:
        raise ValueError("A spectrum must contain at least one positive-intensity peak")
    return [item[0] for item in active], [item[1] for item in active]
