# Conformer development-branch checks (2026-10-06)

This is a development snapshot, not a release or a claim that the full-budget
comparison is finished. The running comparison uses its frozen pre-fix source.

## Fixes in this snapshot

1. Included model workers no longer import an external benchmark loader.
   Source/weights are configured by the user; optional CPU requirements are
   separate. See `model_installation.md` for verified vs unverified setup paths.
2. Shared geometry checking rejects compressed bonds using the configurable
   minimum covalent-radius ratio of 0.60. All generators use the same rule.
3. New PM6 relaxation classifies native failures in `failure.json`, preserving
   original logs and parser exceptions. Legacy LAQA parsing is unchanged.

Tests: 54 conformer/session/release-fix cases passed. A broader run had 99 passes
and one timing-only failure in the existing subprocess timeout test (2.108 s
against a <2 s assertion). That test was not weakened to make the run pass.
An isolated rerun of the six timeout tests then passed; the timing sensitivity
remains documented rather than being treated as a clean first-pass full suite.
`test_compare_top3_spectra.py` could not collect because its local, Git-ignored
`compare_top3_spectra.py` dependency is absent; it was excluded from the broader
run. `testcal.py` is a manual Gaussian calculation, not a unit test.

## Deferred maintainer handoff (not part of this change)

At completion, tell Masato Sumita about the existing GAMESS handoff: re-reading
`optimized_structures.sdf` overwrites `TotalCharge` / `SpinMulti`, potentially
losing `SpecTotalCharge` / `SpecSpinMulti` overrides before downstream GAMESS.
This is not evidence of incorrect xTB/PM6 input in the current Sample campaign.
The user explicitly deferred this repair; do not silently include it here.
