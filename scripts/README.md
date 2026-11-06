# Analysis scripts

Published **exactly as they ran** on our cluster, absolute paths and French comments
included. They are the audit trail, not a portable package: the point is that every number
in the report can be traced to the code that produced it. The one script written to run
anywhere is `../preflight/`.

Interpreter: `bioemu2` in its own venv, pinned to upstream commit `06be369`. Hardware: RTX
PRO 6000 Blackwell under Slurm, 23.62 GiB MIG slice. Every `.sbatch` is the exact submission
that ran, timeout included.

## The three scripts behind question 2

| script | what it does |
|---|---|
| `md_contre_x0.py` | scores the **simulation** against the same starting structure BioEmu2 was given. Imports the frozen metric path (`redock_lib.Reference`), reimplements nothing. 28 systems, 84 replicas, 16,800 images, 0 rejects, 50.8 s |
| `test_echangeabilite.py` | the exchangeability test. Three replicas are three draws; if BioEmu2's value were a fourth from the same distribution it would land above all three once in four. It lands above in 19, 17 and 22 of 28 |
| `lag_agregat.py` | the lag sweep aggregate over the 28 complexes: median DockQ at 1, 10, 50, 200, 1000 ns, how many decline, how many strictly |

## `generation/` - sampling and the hydrogen stage

| script | what |
|---|---|
| `gen_lib.py` | the common library: the corpus of 28, reading the cured x0 with its sha256, the append-only manifest, per-unit concurrency, peak-memory capture. A rejection is an **entry** in the manifest, never a silent skip |
| `generer.py` | the production run: U3 (propagation), U2 (equilibrium ensemble), U1 (structure prediction) |
| `capacite.py` | the capacity curve: VRAM allocated *and* reserved, throughput, batch scaling, per system |
| `msa_serie.py` | the ColabFold alignments, **one request at a time**, as the upstream docs require |
| `h_preparer.py`, `h_piles.py`, `h_agreger.py` | the hydrogen stage: BioEmu2 emits none, and nothing about contacts is measurable until they are rebuilt. `h_agreger.py` is descriptive only, no interpretation |
| `md_piles.py` | the reference simulation piles, 28 x 3 replicas, stride 5 (200 images of 1,000), force-field hydrogens |
| `ponts_reference.py`, `croiser_cristal.py` | the disulfide bridges and their SG-SG distance in the simulation, crossed against the crystal and the simulation's own PSF. This is how we found the physical check ignores them |
| `recouvrement.py` | how much of the simulation's explored space the generated ensemble covers |
| `analyse/` | the report-side aggregation: `capacite_analyse.py`, `modele_cout.py`, `ponts_analyse.py`, `rapport.py` |

## `redocking/` - the docking campaign

| script | what |
|---|---|
| `redock_lib.py` | decoys, the chance floor, and the CAPRI measurements. `Reference(ref).mesurer(traj)` returns fnat / irmsd_A / lrmsd_A / dockq / classe_capri and **enforces identical atom names and order** (ordinal pairing, never resSeq) |
| `q2_metriques_vendu.py` | the cheap geometric measures on heavy atoms, CPU only |
| `preparer_AB.py` | the crystal starting structures, each one through the pre-flight (imported, not copied) |
| `construire_negatifs.py` | the negative control: chain pairs with **no** right answer |
| `plancher.py` | the chance floor, with no model at all, and the verification of the decoys exactly as generated |
| `analyser.py` | the aggregation, per label and lag |
| `partenaire_contre_etranger.py` | the paired comparison, chain A against the same chain A: real partner or stranger |
| `site_positifs.py` | how much of the native site is reached |
| `confusions.py`, `homologie_B.py`, `profondeur_msa.py` | the sensitivity analyses: which complexes get confused, how much a training-set relative buys, how alignment depth matters |
| `redock_gpu.py` | the GPU harness. Note: it **exits 1 when any rejection is counted**, so a `FAILED` state here can mean one declared rejection, not a failed run |

## `t4l/` - T4 lysozyme, on Bhakat's pipeline

| script | what |
|---|---|
| `cv_t4l.py` | Bhakat's collective variables, applied to an mdtraj trajectory |
| `reproduire_v1.py` | **the first measurement**: his published numbers recomputed with our code on his data. Nothing else is trustworthy until this matches |
| `refs_fermees_cv.py` | the closed references and their CV values |
| `bioemu2_t4l.py`, `bioemu2_t4l_steer.py` | BioEmu2 on T4L: units (equilibrium / propagate / steer), the CVs computed in line, resumable |
| `analyser_t4l.py`, `ks_unites.py`, `rmsd_cristal.py`, `rmsd_ferme.py`, `proxy_cb.py` | the analyses |
| `verifier_enligne.py`, `controle_rapport.py` | the guards each run checked itself against, and the report's own control |
| `inspecter_5lzm.py`, `empreintes.sh` | the crystal inspection and the fingerprints of everything used |

## `corpus/` - where the inputs came from

`s1_fetch.py` then `s3_download.py` fetch the DynaBench projects from `dynarepo.inria.fr`;
`s3b_verify.py` and `s6_sums.py` verify them; `s2_analyse.py`, `s5_stats.py`,
`s8_stats_v2.py`, `s10_stats_post_v2.py` characterise them; `s4_manifest.py` and
`s4b_manifest_v2.py` write the corpus manifest; `s7_commensurabilite.py` is the check that
the object simulated and the object handed to BioEmu2 carry the same residues;
`s12_psf_ruptures.py` finds the chain breaks.

The trajectories themselves are **not** in this repository: see the licence note in the
top-level README.
