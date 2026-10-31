# v1/VIAL script pairs

The acceptance gate of VIAL (script v2) against v1, run by `test/test_script_pairs.cpp`.

- A `.pair` file is two scripts that compute the same thing, after `=== v1` and `=== vial`. Every output the v1 script
  publishes has to be published by the VIAL script under its name, with as many values and the same values within the
  tolerance, compared in v1's unit and around its circle where the output is periodic.
- A `.filters` file is lines of `v1 filter || VIAL filter`, which have to pick the same atoms.
- The header names the data, from `ext/mdlib/test_data`: `# data: 1ALA-560ns.pdb` (a PDB of several models is its own
  trajectory), or a structure and its trajectory, `# data: catalyst.gro catalyst.xtc`. `# tol: 1e-3 1e-4` sets the
  absolute and relative tolerance.
- A comment says what a pair checks and, where the languages differ in what they mean, how the two scripts read alike:
  v1's `within` leaves out what it is around and VIAL's does not; v1's `rmsd` weighs by mass and VIAL's does when told;
  positions are 1-based in v1 and 0-based in VIAL.

## What VIAL does not have yet

v1's procedures without a VIAL equivalent, so without pairs:

- Bond changes per frame: `bond_delta`, `formed`, `broken` (declared in VIAL, not evaluated).
- Volumes and distributions over space: `sdf`, `density`, `density_x/y/z`, `rdf` (phase 5 and after).
- Contacts: `contacts`, `contact_count`, `degree`, `porosity`, `pool`.
- Ring geometry: `plane`, `normal`, `radius`, `deviation`, `pucker`, `amplitude`, `phase`, `theta`, `pseudorotation`;
  and a SMILES pattern along a ring in `rings`.
- `shape_weights`.
- Selectors: `qm`, `environment`, `component` and `match`.

What v1 cannot publish has no pair either: a value of more than one dimension per frame, such as `coord` of several
atoms, is not a v1 property.
