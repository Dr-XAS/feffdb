# FeffDB

Feff DB is an SQLite3 database and Python library for the results of
Feff EXAFS calculations, and some support code to read and work with
this database and the data it contains.  The initial motivation is to
provide "first-shell EXAFS" calculations for many common ligands or
atom pairs, so that these can be browsed and re-used without having to
setup and run Feff for many common cases.


## Dr-XAS candidate expansion

Two separately archived expansion rounds add 122 single-scattering first-shell
paths while preserving the earlier 358-row oxide-pilot candidate database:

- [Round 1](source_data/drxas_element_round1_v1/README.md): 16 paths, 10
  element/edge references, 11 absorbing sites; candidate SQLite has 374 rows.
- [Round 2](source_data/drxas_element_round2_v1/README.md): 106 paths, 19
  references, 29 absorbing sites; candidate SQLite has 480 rows. This adds 13
  lanthanide absorbers and Ga-O, In-O, Hf-O, Al-N, Si-N and S-O environments.
- [Oxide radius study](source_data/drxas_oxide_sensitivity_v1/README.md): 15
  calculations comparing SCF and cluster radii for the earlier Ir/Pt/Ta pilot.

Each package carries SHA-256 inventories, exact source/derived-CIF provenance,
inputs, native outputs and executable fingerprints, validation records and tests.
The archived descriptions of local-draft status record the calculation-time
state; the Git commit containing these unchanged archives supplies their
publication identity. Source tests cover 18 round-1, 24 round-2 and 7 radius-study
cases. Publishing the evidence does not promote scientific status: all new paths
remain candidates, with no experimental validation. The radius study applies
only to its three oxide references.

These separate candidate SQLite files do not replace the upstream distributed
database. Generic FeffDB clients do not enforce the extra candidate-review table;
consumers must preserve review status and source limitations.
