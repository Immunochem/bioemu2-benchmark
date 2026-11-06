"""J8 — la topologie SIMULEE est-elle continue ? (lecture de `topology.psf`, liaisons C(i)-N(i+1) declarees)

Le preparateur mesure une DISTANCE sur `structure.pdb`. Ici on lit ce que la MD a vraiment simule : le PSF NAMD declare
des liaisons. Un residu reconstruit relie par une liaison au PSF fait partie de la chaine simulee ; une jonction SANS liaison
au PSF est une rupture de la topologie. Analyseur de PSF MAISON (connectivite seulement : idx, segid, resid, nom d'atome,
paires de la section NBOND) ; il ne decide rien sur les gaps RCSB.
"""
import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from corpus_lib import ROOT, jload, now_iso  # noqa: E402


def lire_psf(chemin):
    atomes, liaisons, mode = [], set(), None
    with open(chemin, encoding="utf-8", errors="replace") as f:
        lignes = f.readlines()
    i = 0
    n_at = n_b = 0
    while i < len(lignes):
        l = lignes[i]
        if "!NATOM" in l:
            n_at = int(l.split()[0])
            for j in range(n_at):
                t = lignes[i + 1 + j].split()
                atomes.append((int(t[0]), t[1], t[2], t[3], t[4]))      # idx, segid, resid, resname, name
            i += n_at
        elif "!NBOND" in l:
            n_b = int(l.split()[0])
            vals = []
            k = i + 1
            while len(vals) < 2 * n_b:
                vals += [int(x) for x in lignes[k].split()]
                k += 1
            for a, b in zip(vals[0::2], vals[1::2]):
                liaisons.add((min(a, b), max(a, b)))
            i = k
            break
        i += 1
    return atomes, liaisons, n_at, n_b


def main():
    sel = jload(ROOT + "/selection_v2.json")["selection"]
    tot_pairs = tot_rupt = 0
    for s in sel:
        pdb = s["pdb"]
        out = {}
        for rep in (1, 2, 3):
            atomes, liaisons, n_at, n_b = lire_psf("%s/systemes/%s/rep%d/topology.psf" % (ROOT, pdb, rep))
            seg = {}
            for idx, segid, resid, resname, name in atomes:
                seg.setdefault(segid, {}).setdefault((resid, resname), {})[name] = idx
            n_pairs = n_rupt = 0
            ruptures = []
            for segid, residus in seg.items():
                cles = list(residus.keys())
                for a, b in zip(cles, cles[1:]):
                    C, N = residus[a].get("C"), residus[b].get("N")
                    if C is None or N is None:
                        continue
                    n_pairs += 1
                    if (min(C, N), max(C, N)) not in liaisons:
                        n_rupt += 1
                        ruptures.append("%s %s%s->%s%s" % (segid, a[1], a[0], b[1], b[0]))
            out[rep] = (len(seg), n_pairs, n_rupt, ruptures, n_at, n_b)
        r1 = out[1]
        tot_pairs += r1[1]
        tot_rupt += r1[2]
        print("PSF %-5s segments=%s  paires_residus_consecutifs=%d  jonctions_C-N_SANS_liaison_au_PSF=%s  (rep1: atomes=%d liaisons=%d)  psf_identique_3_repliques=%s %s" % (
            pdb, [out[k][0] for k in (1, 2, 3)], r1[1], [out[k][2] for k in (1, 2, 3)], r1[4], r1[5], all(out[k][2] == r1[2] for k in (2, 3)),
            ("ruptures=%s" % r1[3][:6]) if r1[2] else ""))
    print("PSF_BILAN systemes=%d  paires_C-N_testees(rep1)=%d  sans_liaison=%d" % (len(sel), tot_pairs, tot_rupt))
    print("FINI s12 %s" % now_iso())


if __name__ == "__main__":
    main()
