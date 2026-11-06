# BioEmu2, benchmarked

**[Read the report →](https://immunochem.github.io/bioemu2-benchmark/)**

BioEmu2 (Microsoft Research, public preview 2026-10-20) exists to produce
molecular-dynamics-like ensembles, faster. So there are only three questions, and this
repository holds the measurements for all three.

Measured on **28 protein-protein complexes** from public DynaBench projects, each already
simulated 3 x 100 ns, plus **T4 lysozyme**, put through the same pipeline Soumendranath
Bhakat used to score eight other generative methods
([JCIM 2026, 10.1021/acs.jcim.6c02044](https://doi.org/10.1021/acs.jcim.6c02044)). The
eight-method panel is his work, not ours. BioEmu2 was released after his paper, so it is
not in it; we ran it through his pipeline and checked our implementation against his
published arrays to 5e-7 nm.

## The three answers

**1. Are the structures accurate?** The structures yes, the distribution no. Median Ca RMSD
to the crystal is 0.66 to 0.78 A in every sampling mode. But the ensemble is narrower than
the version it replaces (8 bins of shape space against 26 for BioEmu v1), and per complex
the outcome is unpredictable: 12 of 28 move less than their own simulation, 8 match it, 8
move more. The corpus average lands 0.105 A from the simulation only because the two
failing groups cancel out.

**2. Does the one-microsecond lag hold?** No. **Ask for 1000 ns and you get about 10 ns of
physics.** At every matched time BioEmu2 has moved less than the simulation, and it sits
above its own simulation's three-replica range in 19, 17 and 22 of 28 complexes against 7
expected under exchangeability (binomial p down to 4e-9). At 1000 ns propagation also
overshoots the model's *own* equilibrium (15.5 % buried against 6.4 %, Fisher p = 0.002),
which needs no external reference at all. The dial is monotone and it is not a clock: it
runs roughly two orders of magnitude slow. `propagation.py:242` accepts any positive
number, with no upper bound.

**3. Can it dock two proteins it was never trained on?** On proteins it has seen, it tells a
real partner from a stranger: +0.53 [+0.25, +0.62] of the native site reached. On proteins
deposited after its training cutoff, +0.08 [-0.14, +0.25] - not measurable. Before trusting
a hit against a new target, check which case you are in.

## Five things stop you at the door

All five are detectable on a laptop CPU, before you book a GPU. Of 94 structures taken
straight from a public database, none were accepted as they were.

1. Atom names must follow PDB convention. ILE carries `CD1`, not `CD`.
2. Standard amino acids only. No modified residues.
3. Chain IDs sequential from A. `A,B`, never `A,C`.
4. Histidine protonation is read from your reference and must stay constant across frames.
5. About 835 residues is the practical ceiling on a 24 GB card.

BioEmu2 checks none of this until the model is loaded on the GPU, and the error it raises
names one violation at a time. `preflight/` is the pre-flight that is missing: it runs on a
laptop, needs no GPU, and reports all five at once.

```bash
python preflight/curer_entree.py your_structure.pdb x0_out.pdb --rapport report.json
# exit 0 = ready (x0 written)   3 = refused, with the class of every blocker   2 = usage
```

It curates what preserves the object (renaming an atom, a protonation state, a chain;
renumbering) and **refuses** what would change it (removing an adduct, bridging a chain
break, deleting an ion). Every fix is declared by name in the report.

## What is in here

| path | what |
|---|---|
| `docs/index.html` | the report, as published |
| `preflight/` | the entry checker, CPU only, no GPU, no network |
| `scripts/generation/` | BioEmu2 sampling, the lag sweep, the hydrogen rebuild, the DynaPIN piles |
| `scripts/redocking/` | the re-docking campaign: metrics, floors, negative controls, strata |
| `scripts/t4l/` | T4 lysozyme: Bhakat's two coordinates rebuilt, the sampling, the analyses |
| `scripts/corpus/` | how the DynaBench corpus was fetched and verified |
| `scripts/md_contre_x0.py` | the simulation scored against the **same** starting structure |
| `scripts/test_echangeabilite.py` | the exchangeability test behind question 2 |
| `scripts/lag_agregat.py` | the lag sweep aggregate across the 28 complexes |
| `data/md_contre_x0/` | 16,800 per-image measurements, 28 systems x 3 replicas, 0 rejects |
| `data/redocking/` | every re-docking measurement, per image and aggregated |
| `data/t4l/` | T4 lysozyme results, including the calibration against the published panel |
| `data/comparaison/` | BioEmu2 against the simulation, per complex |
| `data/dynabench_index.tsv` | the 210 public DynaBench projects, with the 28 we used |

## Reproducing a number

Every figure on the page comes from a file in `data/`. The scripts are published **exactly
as they ran on our cluster, absolute paths included** - they are the audit trail, not a
portable package. The one script written to run anywhere is `preflight/`, and it carries its
one declared deviation at the top of the file (the two pinned dependencies sit next to it;
their sha256 are unchanged and still verified at load time).

Hardware: BioEmu2 on an RTX PRO 6000 Blackwell under Slurm, on a 23.62 GiB MIG slice.
Reference simulations are the public DynaBench projects (3 x 100 ns, 0.1 ns per frame, NAMD
2.13, CHARMM36m, TIP3P, 300 K). Simulation wall-clock measured in house on a whole A100,
systems of 47,000 to 253,000 particles.

## Trajectories are not redistributed

The DynaBench projects we used carry **no declared licence** (`license` is empty for all of
them in `data/dynabench_index.tsv`, and `published` is false). We therefore publish the
**derived measurements** and the fetch path, never the trajectories. `scripts/corpus/`
re-downloads them from the source: `s1_fetch.py` then `s3_download.py`, against
`dynarepo.inria.fr`. `data/dynabench_index.tsv` gives the accession of every project, so the
exact inputs can be recovered.

## Licence

Code (`preflight/`, `scripts/`): MIT, see [LICENSE](LICENSE).
Report and figures (`docs/`): CC BY 4.0.
Derived measurements (`data/`): CC BY 4.0. The upstream trajectories are not ours to
license, and are not here.

## Credit, and a note on one figure

The eight-method T4 lysozyme panel is **S. Bhakat's**, published work we build on, not ours.
BioEmu2 itself is Microsoft Research's; its model card is unusually candid, stating that
outputs are "plausible hypotheses, not measurements" and that known limitations "have not
been quantitatively characterized". That sentence is why this benchmark exists.

Third-party summaries quote speedups of "millions of times" over molecular dynamics. We did
not find that claim anywhere in the project's own documentation, and it is not what we
measured. The honest ratio: 200 to 430x less wall-clock for one ensemble, and a lag dial
that delivers about a hundredth of the physical time you ask it for.

---

Philippe Herman, [Immunochem](https://immunochem.be). Last measured 2026-11-06.
Found an error? Open an issue - the data is here so that you can.
