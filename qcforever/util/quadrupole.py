"""Quadrupoles in Debye Angstrom, with the Gaussian Cartesian convention.

Q is the primitive charge second moment and T = Q - tr(Q) I / 3.
GAMESS prints Buckingham's tensor, 3 T / 2, so its printed values
must be multiplied by 2/3.  Its primitive trace cannot be recovered.
"""

import numpy as np


UNIT = "Debye Angstrom"
BOHR_TO_ANGSTROM = 0.529177210903
E_ANGSTROM_TO_DEBYE = 4.80320471257
COMPONENTS = ("XX", "YY", "ZZ", "XY", "XZ", "YZ")


def tensor_from_components(values):
    xx, yy, zz, xy, xz, yz = (values[key] for key in COMPONENTS)
    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]], dtype=float)


def _components(tensor):
    return dict(unit=UNIT, **{
        key: float(tensor[i, j])
        for key, (i, j) in zip(COMPONENTS, ((0, 0), (1, 1), (2, 2),
                                           (0, 1), (0, 2), (1, 2)))
    })


def traceless(tensor):
    return tensor - np.eye(3) * np.trace(tensor) / 3.0


def shift_primitive(tensor, dipole, charge, displacement):
    """Move the origin by displacement (Angstrom); dipole is in Debye.

    r_new = r_old - displacement.  The charge is in elementary charges.
    All vectors and the tensor must refer to the same Cartesian frame.
    """
    displacement = np.asarray(displacement, dtype=float)
    dipole = np.asarray(dipole, dtype=float)
    return (tensor - np.outer(displacement, dipole)
            - np.outer(dipole, displacement)
            + charge * E_ANGSTROM_TO_DEBYE * np.outer(displacement, displacement))


def quadrupole_properties(primitive=None, traceless_tensor=None, **metadata):
    """Build additive output keys; missing data is None, never a zero tensor."""
    if primitive is not None:
        traceless_tensor = traceless(primitive)
    result = {
        "quadrupole": None if primitive is None else _components(primitive),
        "quadrupole_traceless": (
            None if traceless_tensor is None else _components(traceless_tensor)),
        "quadrupole_iso": None if primitive is None else float(np.trace(primitive) / 3),
        "quadrupole_anisotropy": (
            None if traceless_tensor is None
            else float(np.sqrt(1.5 * np.sum(traceless_tensor ** 2)))),
        "quadrupole_info": dict(unit=UNIT, convention="Q - trace(Q) I / 3", **metadata),
    }
    return result


def quadrupole_at_center_of_mass(primitive, dipole, charge, coordinates, masses, **metadata):
    """Shift an extracted primitive tensor to the COM and compute its scalars."""
    masses = np.asarray(masses, dtype=float)
    if not np.all(np.isfinite(masses)) or np.any(masses <= 0):
        raise ValueError("Invalid atomic masses for the center of mass")
    center = np.average(coordinates, axis=0, weights=masses)
    tensor = shift_primitive(primitive, dipole, charge, center)
    return quadrupole_properties(
        primitive=tensor, origin="center_of_mass",
        origin_coordinates=center.tolist(), origin_unit="Angstrom",
        source_origin_coordinates=[0.0, 0.0, 0.0], status="ok", **metadata)


def quadrupole_from_buckingham(tensor, **metadata):
    """Normalize an extracted Buckingham tensor before calculating anisotropy."""
    tensor = traceless(np.asarray(tensor, dtype=float) * (2.0 / 3.0))
    return quadrupole_properties(
        traceless_tensor=tensor, source_convention="Buckingham",
        status="traceless_only", reason="The primitive trace cannot be recovered from a Buckingham tensor",
        **metadata)
