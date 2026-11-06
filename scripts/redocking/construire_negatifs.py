#!/usr/bin/env python3
"""construire_negatifs.py — CPU. TEMOIN NEGATIF : des paires de chaines qui n ont AUCUNE bonne reponse. La chaine A vient du complexe X (sa chaine 0, coordonnees cristallines inchangees),
la chaine B vient d un AUTRE complexe Y (sa chaine 1), sans parente de sequence (identite locale < 35 % sur la plus courte des deux), non homologues, et leur taille cumulee <= 700 residus.
La chaine B est posee a la distance d_COM du complexe X le long de l axe natif de X, orientation aleatoire, sans choc (>= 3,0 A). Ce x0 sert de « reference » technique (ordre d atomes,
System) ; il n a pas de bonne reponse. Pour le test de SITE : `site_residus` = residus de la chaine A qui, dans le cristal X, touchent son vrai partenaire (< 5 A).
Passe par la meme curation (curer_entree.py) que les positifs. Rejets comptes.
Usage : construire_negatifs.py <strate_neg> <N> <strate_source> [<strate_source> ...]   ex. NEG_A 16 A_corpus A_apparie"""
import collections, importlib.util, itertools, json, os, sys
from pathlib import Path
os.environ["CUDA_VISIBLE_DEVICES"] = ""
sys.dont_write_bytecode = True
B = Path("/gpfs/phherman/bioemu2_bench"); R = B / "redocking"
sys.path.insert(0, str(R / "scripts"))
import numpy as np, mdtraj as md
from Bio import Align
from scipy.spatial.transform import Rotation
import redock_lib as RL, redock_bioemu as RB

GRAINE = 20261106
TOK_MAX = 700


def charger_curer():
    sp = importlib.util.spec_from_file_location("curer_entree", B / "prep_fix/curer_entree.py")
    m = importlib.util.module_from_spec(sp); sp.loader.exec_module(m); return m


def prets(strates):
    out = []
    for nom in strates:
        m = json.loads((R / f"strate_{nom}/manifeste.json").read_text())
        for pid, e in sorted(m.get("systemes", m).items()):
            if e["etat"] == "PRET":
                out.append((nom, pid))
    return out


def identite_locale(a, b, al):
    s = al.align(a, b)[0]
    # identites sur les blocs alignes
    n = 0
    for (a0, a1), (b0, b1) in zip(*s.aligned):
        n += sum(1 for i, j in zip(range(a0, a1), range(b0, b1)) if a[i] == b[j])
    return n / min(len(a), len(b))


def ecrire_x0(atomes, chemin):
    """atomes : liste de (chaine, resSeq, resname, nom, element, (x, y, z) en A)."""
    L, n = [], 0
    derniere = None
    for (ch, rs, rn, nm, el, (x, y, z)) in atomes:
        if derniere is not None and ch != derniere:
            L.append("TER")
        derniere = ch
        n += 1
        an = nm if len(nm) == 4 else " " + nm.ljust(3)
        L.append("ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s" % (n % 100000, an, rn, ch, rs, x, y, z, 1.0, 0.0, el.rjust(2)))
    L += ["TER", "END"]
    Path(chemin).write_text("\n".join(L) + "\n")


def main():
    nom_neg, N, sources = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
    C = charger_curer()
    systemes = prets(sources)
    info = {}
    for (st, pid) in systemes:
        x0 = R / f"strate_{st}/{pid}/x0.pdb"
        t = md.load(str(x0))
        ch = [[a.index for a in c.atoms] for c in t.topology.chains]
        info[(st, pid)] = {"x0": x0, "traj": t, "idx": ch, "seqs": list(t.topology.to_fasta()), "n_res": [c.n_residues for c in t.topology.chains]}
    al = Align.PairwiseAligner(); al.mode = "local"; al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62"); al.open_gap_score = -11; al.extend_gap_score = -1
    rng = np.random.default_rng(GRAINE)
    cles = sorted(info)
    ordre = list(cles); rng.shuffle(ordre)
    paires, rejets = [], collections.Counter()
    utilises_B = set()
    for i, x in enumerate(ordre):
        if len(paires) >= N:
            break
        cand = [y for y in ordre[i + 1:] + ordre[:i] if y[1] != x[1] and y not in utilises_B]
        for y in cand:
            if info[x]["n_res"][0] + info[y]["n_res"][1] > TOK_MAX:
                rejets["taille_cumulee"] += 1; continue
            # aucune parente entre les DEUX chaines croisees, ni entre la chaine A de X et le partenaire d origine, ni... : on exige seulement A_X vs B_Y et B_Y vs B_X / A_Y vs A_X
            ident = max(identite_locale(info[x]["seqs"][0], info[y]["seqs"][1], al), identite_locale(info[x]["seqs"][1], info[y]["seqs"][1], al), identite_locale(info[x]["seqs"][0], info[y]["seqs"][0], al))
            if ident >= 0.35:
                rejets["parente_de_sequence"] += 1; continue
            paires.append((x, y, ident)); utilises_B.add(y)
            break
    print(f"{nom_neg} : {len(paires)} paires retenues sur {N} demandees ; rejets de candidats {dict(rejets)}", flush=True)
    sortie = {}
    for (x, y, ident) in paires:
        nid = f"NEG_{x[1]}A_{y[1]}B"
        d = R / f"strate_{nom_neg}" / nid
        d.mkdir(parents=True, exist_ok=True)
        tx, ty = info[x]["traj"], info[y]["traj"]
        Xx = tx.xyz[0] * 10.0; Xy = ty.xyz[0] * 10.0
        iA = info[x]["idx"][0]; iB = info[y]["idx"][1]
        # axe natif de X : centre chaine 0 -> centre chaine 1 ; distance d_COM native de X
        cA = Xx[iA].mean(axis=0); cX1 = Xx[info[x]["idx"][1]].mean(axis=0)
        dnat = float(np.linalg.norm(cX1 - cA)); u = (cX1 - cA) / dnat
        from scipy.spatial import cKDTree
        arbre = cKDTree(Xx[iA])
        cB = Xy[iB].mean(axis=0)
        # pose AU CONTACT (distance atomique minimale 3,5 A) le long de l axe natif de X, orientation aleatoire : un etranger de forme quelconque ne tient pas sans choc a la distance NATIVE de X
        # (premiere version : 10 paires sur 16 refusees, POSE_SANS_CHOC_IMPOSSIBLE). La distance de depart est donc celle d un CONTACT, comme pour les positifs (natif = contact).
        essais, pose = 0, None
        rr = np.random.default_rng(RL.seed_de(nid, "pose"))
        Rot = Rotation.random(random_state=rr).as_matrix()
        base = (Xy[iB] - cB) @ Rot.T
        f = lambda tt: float(arbre.query(base + cA + u * tt)[0].min())
        tt = dnat + 120.0
        if f(tt) >= RL.DMIN_CONTACT_A:
            while tt > 0 and pose is None:
                essais += 1
                tt -= 2.0
                if f(tt) < RL.DMIN_CONTACT_A:
                    lo, hi = tt, tt + 2.0
                    for _ in range(14):
                        mid = 0.5 * (lo + hi)
                        lo, hi = (mid, hi) if f(mid) < RL.DMIN_CONTACT_A else (lo, mid)
                    pose = base + cA + u * hi
        if pose is None:
            rejets["pose_au_contact_impossible"] += 1; sortie[nid] = {"etat": "REFUSE", "classes_refus": ["POSE_AU_CONTACT_IMPOSSIBLE"]}; continue
        d_pose = float(np.linalg.norm(pose.mean(axis=0) - cA))
        atomes = []
        topX, topY = tx.topology, ty.topology
        for c, (top, X, idx, lettre, coords) in enumerate(((topX, Xx, iA, "A", None), (topY, Xy, iB, "B", pose))):
            ch = list(top.chains)[0 if c == 0 else 1]
            for r in ch.residues:
                for a in r.atoms:
                    k = idx.index(a.index)
                    xyz = X[a.index] if coords is None else coords[k]
                    atomes.append((lettre, r.resSeq, r.name, a.name, a.element.symbol, tuple(float(v) for v in xyz)))
        ecrire_x0(atomes, d / "brut.pdb")
        cur = C.curer(str(d / "brut.pdb"), sortie_x0=str(d / "x0.pdb"), ecrire_diag=str(d / "x0.diagnostic.pdb"))
        (d / "curation.json").write_text(json.dumps(cur, indent=1, ensure_ascii=False, default=str) + "\n")
        # site natif de la chaine A dans X : residus (ordinaux de la chaine 0) qui touchent la chaine 1 du cristal d origine
        Rx = RL.Reference(md.Trajectory(xyz=tx.xyz[:1], topology=tx.topology))
        site = sorted({p[0] for p in Rx.contacts})
        sortie[nid] = {"etat": cur["etat"], "classes_refus": cur.get("classes_refus", []), "source_A": x[1], "strate_source_A": x[0], "source_B": y[1], "strate_source_B": y[0],
                       "identite_locale_max": round(ident, 3), "sequences": [info[x]["seqs"][0], info[y]["seqs"][1]], "tokens": info[x]["n_res"][0] + info[y]["n_res"][1],
                       "d_nat_source_A": round(dnat, 2), "d_pose_contact_A": round(d_pose, 2), "essais_pose": essais, "site_chaine": 0, "site_residus": site, "taille_chaine_site": info[x]["n_res"][0]}
        print(f"{nom_neg} {nid} {cur['etat']} {cur.get('classes_refus', [])} tokens={sortie[nid]['tokens']} identite_max={ident:.2f} essais_pose={essais} site={len(site)} residus", flush=True)
    (R / f"strate_{nom_neg}").mkdir(exist_ok=True)
    (R / f"strate_{nom_neg}/manifeste.json").write_text(json.dumps({"graine": GRAINE, "sources": sources, "systemes": sortie, "rejets_candidats": dict(rejets)}, indent=1, ensure_ascii=False) + "\n")
    c = collections.Counter(v["etat"] for v in sortie.values())
    print(f"BILAN_{nom_neg} : {dict(c)} · REJETS = {c.get('REFUSE', 0)} · paires demandees {N}", flush=True)


if __name__ == "__main__":
    main()
