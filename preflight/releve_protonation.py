#!/usr/bin/env python3
"""releve_protonation.py — l'etat de protonation des histidines, LU dans le topology.psf de chaque MD.

Contexte (pilote 6CWG, J3, exit 41) : l'outil d'ajout d'hydrogenes choisit l'etat de l'histidine image
par image, donc le compte d'H varie et l'enumeration des atomes glisse. Le PSF de la MD DECLARE l'etat
CHARMM : HSD (HD1 seul), HSE (HE2 seul), HSP (les deux). Un etat code en dur serait faux sur tout systeme
portant autre chose que HSD. Ce script releve, pour les 28 systemes et leurs 3 repliques :
  - la table residu -> etat lue dans le PSF (segment, resid, nom CHARMM) ;
  - la meme table deduite du structure.pdb (nom du residu ; presence des hydrogenes HD1 / HE2) ;
  - leur accord, et la constance entre repliques.
Lecture seule sur corpus/systemes/. Ecrit curation/protonation_his.tsv et curation/protonation_his.json.
"""
import collections
import json
import sys
from pathlib import Path

RACINE = Path("/gpfs/phherman/bioemu2_bench")
SORTIE = RACINE / "prep_fix" / "curation"
HIS_NOMS = ("HSD", "HSE", "HSP", "HIS", "HID", "HIE", "HIP", "HSC", "HIN")     # CHARMM puis Amber ; tout nom d'histidine rencontre y est
ETAT_DE_H = {(True, False): "HSD", (False, True): "HSE", (True, True): "HSP", (False, False): "sans_H"}


def lire_psf(chemin):
    """{(segment, resid): nom_residu} dans l'ordre du PSF (section !NATOM, un enregistrement par atome)."""
    res = collections.OrderedDict()
    with open(chemin) as f:
        n = None
        for l in f:
            if "!NATOM" in l:
                n = int(l.split()[0])
                break
        if n is None:
            sys.exit(f"🔴 ARRET : {chemin} sans section !NATOM")
        for _ in range(n):
            v = f.readline().split()
            if len(v) < 5:
                sys.exit(f"🔴 ARRET : enregistrement PSF illisible dans {chemin} : {v}")
            res.setdefault((v[1], v[2]), v[3])
    return res


def lire_pdb(chemin):
    """{chaine: [(resSeq+icode, nom, {noms d'atomes})]} dans l'ordre du fichier."""
    par = collections.OrderedDict()
    cle_prec = None
    for l in open(chemin):
        if not l.startswith(("ATOM  ", "HETATM")):
            continue
        k = (l[21], l[22:26].strip() + l[26].strip(), l[17:20].strip())
        if k != cle_prec:
            par.setdefault(k[0], []).append((k[1], k[2], set()))
            cle_prec = k
        par[k[0]][-1][2].add(l[12:16].strip())
    return par


def main():
    L = [json.loads(l) for l in open(RACINE / "corpus" / "manifeste_corpus.json") if l.strip()]
    systemes = [x for x in L if x.get("statut") == "actif_v2" and x["pdb"] not in ("1ZM4", "1FQ1")]
    lignes, resume = [], []
    for x in systemes:
        par_rep = {}
        for r in x["repliques"]:
            d = Path(x["dossier"]) / f"rep{r}"
            psf = lire_psf(d / "topology.psf")
            pdb = lire_pdb(d / "structure.pdb")
            segs = list(dict.fromkeys(k[0] for k in psf))
            chaines = list(pdb)
            # les segments du PSF et les chaines du PDB sont dans le meme ordre (verifie ci-dessous par le nombre de residus)
            n_psf = {s: sum(1 for k in psf if k[0] == s) for s in segs}
            n_pdb = {c: len(v) for c, v in pdb.items()}
            ordre_ok = (len(segs) == len(chaines)) and all(n_psf[s] == n_pdb[c] for s, c in zip(segs, chaines))
            etats = {}
            for s, c in zip(segs, chaines):
                residus_psf = [(k[1], v) for k, v in psf.items() if k[0] == s]
                for (resid, nom_psf), (rs, nom_pdb, atomes) in zip(residus_psf, pdb[c]):
                    if nom_psf not in HIS_NOMS and nom_pdb not in HIS_NOMS:
                        continue
                    etat_h = ETAT_DE_H[("HD1" in atomes, "HE2" in atomes)]
                    coh = (resid == rs) and (nom_psf == nom_pdb) and (nom_psf == etat_h)
                    etats[(c, rs)] = (nom_psf, nom_pdb, etat_h)
                    lignes.append(f"{x['pdb']}\trep{r}\t{s}\t{c}\t{resid}\t{rs}\t{nom_psf}\t{nom_pdb}\tHD1={'HD1' in atomes}\tHE2={'HE2' in atomes}\t{etat_h}\t{'oui' if coh else 'NON'}")
            par_rep[r] = {"ordre_ok": ordre_ok, "etats": etats, "segments": segs, "chaines": chaines}
        # constance entre repliques : meme residu -> meme etat PSF ?
        ref = {k: v[0] for k, v in par_rep[x["repliques"][0]]["etats"].items()}
        constant = all({k: v[0] for k, v in par_rep[r]["etats"].items()} == ref for r in x["repliques"])
        incoh = sum(1 for r in x["repliques"] for k, v in par_rep[r]["etats"].items() if not (v[0] == v[1] == v[2]))
        cnt = collections.Counter(v[0] for v in par_rep[x["repliques"][0]]["etats"].values())
        resume.append({"pdb": x["pdb"], "n_his": sum(cnt.values()), "etats_psf_rep1": dict(cnt), "constant_entre_repliques": constant,
                       "residus_incoherents_psf_vs_pdb": incoh, "ordre_segments_ok": all(par_rep[r]["ordre_ok"] for r in x["repliques"]),
                       "his_rep1": {f"{c}:{rs}": v[0] for (c, rs), v in par_rep[x["repliques"][0]]["etats"].items()}})
    (SORTIE / "protonation_his.tsv").write_text("pdb\treplique\tsegment_psf\tchaine_pdb\tresid_psf\tresSeq_pdb\tnom_psf\tnom_pdb\tHD1\tHE2\tetat_deduit_des_H\tcoherent\n" + "\n".join(lignes) + "\n")
    (SORTIE / "protonation_his.json").write_text(json.dumps(resume, indent=1, ensure_ascii=False) + "\n")
    print(f"{len(systemes)} systemes ; {len(lignes)} lignes (histidine x replique) -> curation/protonation_his.tsv")
    print(f"{'systeme':<8}{'n_HIS':>6}  etats PSF (rep1)                       constant/3 rep.  PSF=PDB=H  ordre segments")
    tot = collections.Counter()
    for r in resume:
        tot.update(r["etats_psf_rep1"])
        print(f"{r['pdb']:<8}{r['n_his']:>6}  {json.dumps(r['etats_psf_rep1']):<38} {str(r['constant_entre_repliques']):<16} {'OK' if r['residus_incoherents_psf_vs_pdb'] == 0 else 'NON x' + str(r['residus_incoherents_psf_vs_pdb']):<10} {r['ordre_segments_ok']}")
    print(f"\ntotal histidines (rep1, 28 systemes) : {dict(tot)} ; systemes portant HSE ou HSP : {[r['pdb'] for r in resume if any(k in r['etats_psf_rep1'] for k in ('HSE', 'HSP'))]} ; "
          f"systemes sans histidine : {[r['pdb'] for r in resume if r['n_his'] == 0]}")
    tot3 = collections.Counter()
    for l in lignes:
        tot3[l.split("\t")[6]] += 1
    print(f"total sur les 3 repliques (lignes) : {dict(tot3)} ; lignes incoherentes (PSF != PDB != H) : {sum(1 for l in lignes if l.endswith('NON'))}")
    print("systemes dont l'etat PSF varie d'une replique a l'autre :", [r["pdb"] for r in resume if not r["constant_entre_repliques"]])
    print("systemes dont l'ordre segments PSF / chaines PDB n'est pas verifie :", [r["pdb"] for r in resume if not r["ordre_segments_ok"]])


if __name__ == "__main__":
    main()
