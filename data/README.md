# Derived measurements

Every figure in the report comes from a file here. The upstream trajectories are **not**
redistributed (see the licence note in the top-level README); these are the measurements
taken from them.

## `md_contre_x0/` - the simulation scored against the same starting structure

The measurement behind question 2. Produced by `scripts/md_contre_x0.py`.

`mesures.json` is one JSON array, 16,800 entries:

```json
{"systeme": "1AY7", "replique": "rep1", "image": 0, "ns": 0.0,
 "dockq": 1.0, "fnat": 1.0, "irmsd_A": 0.0, "lrmsd_A": 0.0}
```

`bilan.json` carries the run's own accounting: `systemes` 28, `repliques` 84, `images`
16,800, `n_rejets` **0**, `ns_par_image` 0.5 (stride 5 over 1,000 frames of 100 ns).

It also carries the one declared deviation, measured rather than waved away:
`oxt_retires` 168 terminal OXT atoms removed to make the atom sets identical, of which
`oxt_a_moins_de_5A` **10** sit within 5 A of the partner chain, the closest at
`oxt_dmin_min_A` **2.76 A** (systems 1GCQ, 1SYX, 2C0L, 2X9A). `oxt_detail` lists each one.
Runtime 50.8 s.

**Why this file exists.** BioEmu2 piles and simulation piles were not commensurable: the
scoring tool compares every frame to image 0 of its own stack, and those image 0s were
different objects. This script scores the simulation against the **same** starting structure
BioEmu2 was given, through the frozen metric path (`redock_lib.Reference`), reimplementing
nothing.

## `redocking/` - the re-docking campaign

`images.csv`, `images_lag0050.csv`, `images_lag1000.csv`: one row per generated structure,
27 columns - `dockq`, `fnat`, `irmsd_A`, `lrmsd_A`, `classe_capri`, the site-recall columns
(`precision_site`, `rappel_site`), the per-chain drift (`derive_A_ca_A`, `derive_B_ca_A`),
and the arm/condition/stratum labels that say which test each row belongs to.

| file | what |
|---|---|
| `tableaux_lag1000.{json,md}`, `tableaux_lag0050.*` | the aggregates, per arm and condition |
| `site_positifs_lag1000.json` | how much of the native site is reached |
| `partenaire_contre_etranger_lag1000.*` | real partner against a stranger, the question a binder screen asks |
| `confusions_lag1000.*` | which complexes get confused with which |
| `plancher_resume.json`, `plancher_par_strate.json`, `plancher/` | the chance floor and the model's own floor, per stratum and per system |
| `controle_mesure.json` | our metrics against an independent implementation: worst disagreement 0.0005 |

## `t4l/` - T4 lysozyme

| file | what |
|---|---|
| `refs_fermees_cv.json` | **the calibration**: Bhakat's two coordinates rebuilt from scratch, matched to his published arrays |
| `analyse_unites.json` | the sampling units, the d1 hinge coordinate, the crystal comparison |
| `rmsd_cristal.json`, `rmsd_ferme_unites.json` | Ca RMSD to the open crystal and to the closed reference |
| `ks_unites.json` | the distribution comparisons |
| `reproduction_v1.json` | BioEmu v1 rerun for the same panel |
| `verif_enligne_*.json` | the in-line guards each run checked itself against |
| `*_steer.json`, `*_smoke.json` | the steered variant and the smoke run, kept for the audit trail |
| `tables_rapport.md` | the tables as the internal report printed them |

## `t4l/cv_d1_p_ep1.json` - d1 against p, the plot behind episode 1

d1 and the locking coordinate p for the **seven published methods that carry the Phe4 ring** (2,197 structures,
read from the author's own `distances.npy` and `locked.pkl`; none recomputed here) and for BioEmu2
(`eq164_n500`, 500 structures, from our `cv.npz`). 0 rejects; nothing falls below d1 = 2.5 nm (the closest is
2.68 nm). `scripts/t4l/dump_cv_ep1.py` exports it, `scripts/t4l/build_ep1_figure.py` draws
`docs/figs/d1_p_bhakat_bioemu2.svg` from it, byte for byte.

**A correction, in the author's words.** Cell 48 of the paper's notebook reads the locking coordinate from
column 1 of `distances.npy`, which is d2, not p. S. Bhakat confirmed it and re-plotted the figure (2026-11-06
and 2026-11-09, author communication; an amendment is to be raised with the journal). Along the corrected p
both states are populated. What does not change, in his words: *"the overall conclusion remains unchanged: the
GenAI ensemble alone does not capture the closed state."*

## `comparaison/` - BioEmu2 against the simulation, per complex

`u2_vs_u3.json`: `par_systeme` holds each complex judged against the spread of its own three
simulation replicas, which is where the 12 / 8 / 8 split comes from. `classement` is the
ordering.

`recouvrement.json`: how much of the simulation's explored space the generated ensemble
covers.

## `dynabench_index.tsv`

The 210 public DynaBench projects: accession, PDB id, replicas available, atom and residue
counts, frames, framestep, engine and force field, and the DB55 interface class. The
`license` column is **empty for every project** and `published` is false, which is why no
trajectory is republished here. `scripts/corpus/s1_fetch.py` then `s3_download.py`
re-download them from `dynarepo.inria.fr`.
