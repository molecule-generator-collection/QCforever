`water_rhf_6-31gstar_5d.fchk` is a minimal formatted-checkpoint fixture,
generated independently with PySCF 2.14.0 RHF/6-31G*, spherical basis functions,
and SCF energy tolerance 1e-12. It is not a complete Gaussian restart file.
Only the geometry, electronic state, basis and orbital fields needed for native
GAMESS/ORCA handoff are included. No PySCF installation is needed to read it.

Geometry (Angstrom): O (0,0,0), H (0,0,0.96), H (0.92,0,-0.24).
Charge 0, multiplicity 1. Source RHF energy: -76.00915487179412 Hartree.

For each PySCF shell, `bas_exp` and `bas_ctr_coeff` supply the exponents and
normalized-primitive contraction coefficients. MO rows are kept in shell order,
with PySCF spherical d components (xy,yz,z2,xz,x2-y2) reordered to Gaussian
(z2,xz,yz,x2-y2,xy). Arrays use the usual FCHK precision (eight decimal places
in scientific notation). MO coefficients are stored MO-major.

The opt-in installed-ORCA test imports this density with `NoIter`, checking the
energy without allowing SCF iterations to repair an incorrect initial guess.
