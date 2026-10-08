"""Common finite-coordinate, graph, clash and diversity checks.

Input-specified stereochemistry is checked from coordinates during generation
and after relaxation. No bond-order repair or molecule-specific rescue is used.
"""
from collections import Counter
from copy import deepcopy
import math

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign


STEREO_MATCH_PROPERTIES = ('tetrahedral_stereo_match', 'ez_stereo_match', 'joint_stereo_match')


def filter_candidates(records, reference, maximum, settings):
    """Validate the entire pool, including after coordinates change in MM."""
    reference = Chem.AddHs(Chem.RemoveHs(reference))
    return _filter_new_candidates(records, reference, maximum, settings, [], [])


class IncrementalCandidateFilter:
    """Retain an unchanged, validated pool across generation batches only.

    The reference and settings are fixed for this instance. Private copies keep
    callers from invalidating cached checks by modifying returned molecules or
    audits. Use snapshot() when bypassing MM; run a full filter on MM outputs.
    Audit indices retain the
    full-filter convention: accepted prefix followed by the new batch, not
    cumulative raw-generation indices.
    """

    def __init__(self, reference, maximum, settings):
        self._reference = Chem.AddHs(Chem.RemoveHs(Chem.Mol(reference)))
        self._maximum = maximum
        self._settings = deepcopy(settings)
        self._accepted = []
        self._accepted_audit = []

    def extend(self, records):
        accepted = [Chem.Mol(mol) for mol in self._accepted]
        audit = deepcopy(self._accepted_audit)
        # The old full-filter call reindexed its accepted prefix on every batch.
        for i, (mol, row) in enumerate(zip(accepted, audit)):
            mol.SetIntProp('raw_candidate_index', i)
            row['raw_index'] = i
        accepted, summary = _filter_new_candidates(
            records, self._reference, self._maximum, self._settings, accepted, audit)
        self._accepted = [Chem.Mol(mol) for mol in accepted]
        self._accepted_audit = deepcopy([row for row in summary['candidates'] if row['accepted']])
        return accepted, summary

    def snapshot(self):
        """Export private, validated generation results; never accept MM outputs.

        Copies prevent callers from changing the stored pool or audit. Indices
        are normalized to the final pool, just as in a fresh full-filter audit.
        """
        accepted = [Chem.Mol(mol) for mol in self._accepted]
        audit = deepcopy(self._accepted_audit)
        for i, (mol, row) in enumerate(zip(accepted, audit)):
            row['raw_index'] = i
            mol.SetIntProp('raw_candidate_index', i)
        return accepted, _candidate_audit(len(accepted), accepted, audit, self._reference, self._settings)


def evaluate_optimized_candidates(records, energies_hartree, reference, settings):
    """Return post-optimization best_any and stereo-valid primary indices."""
    if len(records) != len(energies_hartree):
        raise ValueError('One energy is required for each structure')
    valid, primary, audits = [], [], []
    reference = Chem.AddHs(Chem.RemoveHs(reference))
    for i, (mol, energy) in enumerate(zip(records, energies_hartree)):
        reason = geometry_failure(mol, reference, settings)
        row = {'candidate': i, 'energy_hartree': energy, 'geometry_failure': reason,
               'stereo_check_status': 'not_evaluated_geometry_invalid' if reason else 'not_evaluated_energy_invalid'}
        if reason is None and energy is not None and math.isfinite(energy):
            valid.append(i)
            row.update(stereo_audit(mol, reference))
            row['stereo_check_status'] = 'evaluated'
            if stereo_matches_policy(row, settings):
                primary.append(i)
        audits.append(row)
    return {'primary_valid_candidates': len(primary),
            'best_any_index': min(valid, key=lambda i: energies_hartree[i]) if valid else None,
            'primary_best_index': min(primary, key=lambda i: energies_hartree[i]) if primary else None,
            'candidates': audits}


def select_optimized_candidate(records, energies_hartree, audit):
    """Prefer valid stereo; retain the best available structure with a warning otherwise."""
    primary = audit['primary_best_index']
    selected = primary if primary is not None else audit['best_any_index']
    if selected is None:
        available = [i for i, mol in enumerate(records) if mol is not None]
        selected = min(available, key=lambda i: energies_hartree[i]) if available else None
    audit['selected_index'] = selected
    if primary is not None:
        warning = None
    elif audit['best_any_index'] is not None:
        warning = 'stereo_mismatch'
    elif selected is not None:
        warning = 'geometry_invalid'
    else:
        warning = 'no_converged_structure'
    audit['selected_structure_warning'] = warning
    return audit


def _filter_new_candidates(records, reference, maximum, settings, accepted, audit):
    """Shared checks in the original order; only the trusted prefix is skipped."""
    prefix_size = len(accepted)
    for i, source in enumerate(records, start=prefix_size):
        try:
            reason = geometry_failure(source, reference, settings)
            mol = Chem.Mol(source) if reason is None else None
            stereo = stereo_audit(mol, reference) if mol is not None else {}
            if reason is None and not stereo_matches_policy(stereo, settings):
                reason = 'stereo_mismatch'
            if reason is None:
                if settings.duplicate_rmsd_angstrom > 0 and any(duplicate_rmsd(mol, old, settings) <
                       settings.duplicate_rmsd_angstrom for old in accepted):
                    reason = 'duplicate'
            if reason is None and len(accepted) >= maximum:
                reason = 'maximum_candidates'
            row = {'raw_index': i, 'accepted': reason is None, 'reason': reason}
            row.update(stereo)
            if mol is not None:
                row['duplicate_excluded_atom_indices'] = sorted(xh3_hydrogen_indices(mol))
            if reason is None:
                for key in STEREO_MATCH_PROPERTIES:
                    mol.SetBoolProp(key, row[key])
                mol.SetProp('stereo_check_stage', 'candidate_preparation')
                mol.SetProp('stereo_check_status', 'evaluated')
                mol.SetIntProp('raw_candidate_index', i)
                accepted.append(mol)
        except (ValueError, RuntimeError) as exc:
            row = {'raw_index': i, 'accepted': False, 'reason': 'unreadable_structure',
                   'detail': f'{type(exc).__name__}: {exc}'}
        audit.append(row)
    return accepted, _candidate_audit(prefix_size + len(records), accepted, audit, reference, settings)


def _candidate_audit(raw_count, accepted, audit, reference, settings):
    """Build the shared audit schema; no molecular validation occurs here."""
    return {'raw_generated': raw_count, 'accepted_candidates': len(accepted),
                      'duplicate_rmsd_atoms': 'all_explicit_atoms_except_XH3_hydrogens',
                      'duplicate_excluded_reference_atom_indices': sorted(xh3_hydrogen_indices(reference)),
                      'duplicate_rmsd_angstrom': settings.duplicate_rmsd_angstrom,
                      'duplicate_max_matches': settings.duplicate_max_matches,
                      'duplicate_reflection': False,
                      'rejected': dict(Counter(r['reason'] for r in audit if not r['accepted'])),
                      'candidates': audit}


def geometry_failure(mol, reference, settings):
    """Return the first geometry/identity failure, or None for a valid structure."""
    if mol is None:
        return 'unreadable_structure'
    if mol.GetNumConformers() != 1:
        return 'conformer_count'
    if Counter(a.GetAtomicNum() for a in mol.GetAtoms()) != Counter(a.GetAtomicNum() for a in reference.GetAtoms()):
        return 'atom_composition'
    if graph_key(mol) != graph_key(reference):
        return 'connectivity'
    xyz = np.asarray(mol.GetConformer().GetPositions())
    if not np.isfinite(xyz).all():
        return 'nonfinite_coordinates'
    pt = Chem.GetPeriodicTable()
    for i in range(len(xyz)):
        for j in range(i):
            d = float(np.linalg.norm(xyz[i]-xyz[j]))
            if d < settings.absolute_collision_angstrom:
                return 'collision'
            radii = pt.GetRcovalent(mol.GetAtomWithIdx(i).GetAtomicNum()) + pt.GetRcovalent(mol.GetAtomWithIdx(j).GetAtomicNum())
            if radii <= 0:
                continue
            if mol.GetBondBetweenAtoms(i, j):
                # A deliberately permissive sanity bound, not a bond-length
                # prediction. Apply identically to every generator and element.
                if d/radii < settings.minimum_bonded_covalent_ratio:
                    return 'bond_compressed'
                if d/radii > settings.maximum_bonded_covalent_ratio:
                    return 'bond_stretched'
            elif d/radii < settings.minimum_nonbonded_covalent_ratio:
                return 'collision'
    return None


def stereo_audit(mol, reference):
    """Compare input-specified R/S and E/Z using stereo assigned from coordinates."""
    observed = Chem.Mol(mol)
    Chem.RemoveStereochemistry(observed)
    Chem.AssignStereochemistryFrom3D(observed, confId=0, replaceExistingTags=True)
    observed = Chem.RemoveHs(observed)
    expected = Chem.RemoveHs(Chem.Mol(reference))
    Chem.AssignStereochemistry(expected, cleanIt=True, force=True)
    centers = {a.GetIdx(): a.GetProp('_CIPCode') for a in expected.GetAtoms() if a.HasProp('_CIPCode')}
    bonds = {expected_bond.GetIdx(): expected_bond.GetStereo() for expected_bond in expected.GetBonds()
             if expected_bond.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)}
    tetrahedral_match, ez_match, joint_match = False, False, False
    for mapping in observed.GetSubstructMatches(expected, useChirality=False,
                                                uniquify=False, maxMatches=10000):
        mapping_tetrahedral_match = all(observed.GetAtomWithIdx(mapping[i]).HasProp('_CIPCode') and
                observed.GetAtomWithIdx(mapping[i]).GetProp('_CIPCode') == cip for i, cip in centers.items())
        mapping_ez_match = True
        for i, stereo in bonds.items():
            expected_bond = expected.GetBondWithIdx(i)
            observed_bond = observed.GetBondBetweenAtoms(mapping[expected_bond.GetBeginAtomIdx()], mapping[expected_bond.GetEndAtomIdx()])
            if observed_bond is None or observed_bond.GetStereo() != stereo:
                mapping_ez_match = False
                break
        tetrahedral_match |= mapping_tetrahedral_match
        ez_match |= mapping_ez_match
        joint_match |= mapping_tetrahedral_match and mapping_ez_match
    return {'tetrahedral_stereo_match': tetrahedral_match, 'ez_stereo_match': ez_match,
            'joint_stereo_match': joint_match, 'specified_tetrahedral_centers': len(centers),
            'specified_ez_bonds': len(bonds)}


def stereo_matches_policy(audit, settings):
    """Only input-specified tetrahedral centers/E/Z bonds are constrained.

    This is not a universal isomer classifier: unspecified centers, tautomers,
    atropisomers and coordination stereochemistry are not assigned here.
    """
    tetra = not settings.require_tetrahedral_stereo_for_best or audit['tetrahedral_stereo_match']
    ez = not settings.require_ez_stereo_for_best or audit['ez_stereo_match']
    joint = not (settings.require_tetrahedral_stereo_for_best and settings.require_ez_stereo_for_best) or audit['joint_stereo_match']
    return tetra and ez and joint


def duplicate_rmsd(probe, reference, settings):
    """All atoms except XH3 hydrogens, symmetry-aware fit without reflection.

    Coordinates are not modified. The mapping-search cap bounds work when
    equivalent hydrogens cause combinatorial growth. A capped search may miss
    a duplicate (retain an extra candidate), not guarantee the global minimum.
    """
    return rdMolAlign.GetBestAlignmentTransform(
        rmsd_comparison_molecule(probe), rmsd_comparison_molecule(reference),
        maxMatches=settings.duplicate_max_matches,
        reflect=False, numThreads=1)[0]


def rmsd_comparison_molecule(mol):
    """Copy for RMSD only; omit XH3 hydrogens before symmetry enumeration.

    Removed H atoms become hydrogen counts on X, retaining its chemical
    valence. All other atom coordinates remain unchanged. The original/full
    molecule is still used for geometry checks, MM, semi-empirical optimization and stereo auditing.
    """
    excluded = xh3_hydrogen_indices(mol)
    if not excluded:
        return Chem.Mol(mol)
    comparison = Chem.RWMol(mol)
    for atom in comparison.GetAtoms():
        count = sum(a.GetIdx() in excluded for a in atom.GetNeighbors())
        if count:
            atom.SetNumExplicitHs(atom.GetNumExplicitHs() + count)
    for index in sorted(excluded, reverse=True):
        comparison.RemoveAtom(index)
    result = comparison.GetMol()
    result.UpdatePropertyCache(strict=False)
    return result


def xh3_hydrogen_indices(mol):
    """Indices of H atoms attached to a heavy atom with exactly three H neighbors.

    Element and charge of X do not matter: CH3 and NH3+ are both included.
    X itself is kept, as are H atoms in XH, XH2 and XH4 groups. Inputs must
    contain explicit H atoms, as required by the common composition check.
    """
    excluded = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        hydrogens = [a.GetIdx() for a in atom.GetNeighbors() if a.GetAtomicNum() == 1]
        if len(hydrogens) == 3:
            excluded.update(hydrogens)
    return excluded


def graph_key(mol):
    copy = Chem.RemoveHs(Chem.Mol(mol))
    Chem.RemoveStereochemistry(copy)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=True)


def copy_with_coordinates(reference, coordinates):
    """Take coordinates only; retain input bonds, charges, radicals and labels.

    MM parameterization and SDF sanitization may alter aromaticity/valence
    metadata. Those changes must not become a new molecular identity. This
    does not infer an electronic spin distribution or repair a distorted
    geometry; the usual distance/stereo checks still inspect the new positions.
    """
    if coordinates is None or coordinates.GetNumConformers() != 1:
        raise ValueError('Optimized coordinates require exactly one conformer')
    identity = lambda m: [a.GetAtomicNum() for a in m.GetAtoms()]
    if identity(reference) != identity(coordinates):
        raise ValueError('Optimized atom order/composition changed')
    result = Chem.Mol(reference)
    result.RemoveAllConformers()
    result.AddConformer(Chem.Conformer(coordinates.GetConformer()), assignId=True)
    return result


def clear_structure_audit(mol):
    """Initial-generation checks must not masquerade as final-geometry checks."""
    for key in (*STEREO_MATCH_PROPERTIES, 'geometry_valid', 'geometry_check_failure',
                'structure_check_warning', 'stereo_check_status', 'stereo_check_stage'):
        if mol.HasProp(key):
            mol.ClearProp(key)
