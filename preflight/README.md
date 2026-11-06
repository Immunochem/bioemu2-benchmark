# The pre-flight BioEmu2 does not have

BioEmu2's entry contract validates **entirely on CPU**: `System.from_description` and
`map_cif_atom_order` load no model and need no GPU. We lost GPU jobs one at a time to
blockers this gate sees in 0.2 s, because the tool checks none of them until the model is
loaded, and then names one violation at a time.

```bash
python curer_entree.py <structure.pdb> <x0_out.pdb> [--rapport r.json] [--xtc t.xtc] [--image N]
```

| exit | meaning |
|---|---|
| 0 | **ready** - the gate passes on the x0 that was written |
| 3 | **refused** - every blocker, with its class and the literal message |
| 2 | usage |

## The criterion: commensurability

The object you hand BioEmu2 and the object the simulation simulated must carry **the same
residues**.

- Renaming an atom (`CD` -> `CD1`), a protonation state (`HSD` -> `HIS`), a chain
  (`C` -> `B`), renumbering: these preserve the object. **Applied, and declared**, through an
  explicit table.
- Removing the adduct of a `DDE`, bridging a chain break, deleting a residue, erasing an ion
  or a ligand: these change the object. **Refused**, with the class and the message.

A rename **outside** the explicit table refuses the system, even where mdtraj would do it
silently: what mdtraj converts implicitly is compared to the table, atom by atom.

## What it reports

Real output, on a DynaBench CHARMM structure:

```
dynabench/6CWG/rep3/structure.pdb -> PRET
  [INFO      ] TER_ABSENT          n=0     chains separated by chain identifier alone
  [RENOMME   ] RENOMMAGE_RESIDU    n=4     {'HSD->HIS': 4}
  [RENOMME   ] RENOMMAGE_ATOME     n=28    {'ILE:CD->CD1': 24, 'PHE:OT1->O': 1, ...}
  [ECARTE    ] OXT                 n=2     OXT absent from BioEmu2's contract
  [ECARTE    ] HYDROGENES          n=2946  hydrogens dropped; heavy-atom coordinates unchanged
```

Nothing is dropped in silence. A rejection is a named entry, never a skipped item.

## Requirements

`mdtraj`, `numpy`, `scipy`, and `bioemu2` itself (the gate calls into it, on CPU;
`CUDA_VISIBLE_DEVICES` is cleared by the script). No network.

## Declared deviation in this copy

The published copy resolves its two pinned dependencies (`preparer_cristal.py`,
`preparer_6cwg.py`) next to itself instead of at their cluster paths. Their **sha256 are
unchanged and still verified at load time**, so the functions it imports are the ones that
ran. See the comment at the top of `curer_entree.py`.

## The other scripts here

These are the curation audit trail, published as they ran:
`auditer_28.py` (the 28-system audit), `sonde_porte.py` (the gate's feasibility probe, CPU,
no model), `tokens_vram.py`, `releve_protonation.py`, `liaisons_x0.py`,
`cmp_pilote_6cwg.py`, `resume_systemes.py`, `tests_curation.py`.
