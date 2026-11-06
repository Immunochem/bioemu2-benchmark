#!/usr/bin/env python3
"""croiser_cristal.py — pour les 28 systemes : le pont est-il dans le CRISTAL, dans le PSF de la MD, et a quelle distance SG–SG est-il dans la MD ?
Trois sources, croisees :
  1. mmCIF RCSB (`struct_conn`, conn_type_id = disulf) telecharge dans generation/cristal/ (provenance : url, octets, sha256, duree) ;
  2. liaisons SG–SG du topology.psf de la MD (ponts/reference.json, deja releves par ordinal) ;
  3. cysteines SG a < 3,0 A dans l IMAGE 0 de la MD (x0 cures des 3 repliques) qui ne sont PAS un pont du PSF.
🔴 Numerotation : le cristal parle (auth_asym_id, auth_seq_id, code d insertion) ; la MD parle ordinal. Passage par structure.pdb de la MD (chaine d origine, resSeq,
code d insertion) -> ordinal 1-based dans la chaine du x0 ; un partenaire qui ne s apparie pas est une LIGNE `non_apparie` avec sa raison, jamais un saut.
Controle : le residu apparie doit s appeler CYS/CYX dans la MD. Distances MD : image 0 de chaque replique + MEDIANE de toutes les images de chaque replique."""
import json, sys, time, urllib.request, hashlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_lib as L
import mdtraj as md, numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

CR = L.G / "cristal"
CR.mkdir(exist_ok=True)
C = L.BENCH / "corpus/systemes"
ref = json.loads((L.G / "ponts/reference.json").read_text())["systemes"]
log, tout = [], {}


def telecharger(pdb):
    p = CR / f"{pdb}.cif"
    if p.is_file() and p.stat().st_size > 0:
        return {"pdb": pdb, "statut": "deja_la", "octets": p.stat().st_size, "sha256": L.sha256_fichier(p)}
    url = f"https://files.rcsb.org/download/{pdb}.cif"
    for essai in (1, 2, 3):
        t0 = time.time()
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "tlr4-bioemu2-bench/1.0"})
            data = urllib.request.urlopen(req, timeout=60).read()
            if len(data) < 1000 or b"data_" not in data[:200]:
                raise RuntimeError(f"reponse suspecte ({len(data)} octets)")
            p.write_bytes(data)
            return {"pdb": pdb, "statut": "telecharge", "url": url, "octets": len(data), "sha256": hashlib.sha256(data).hexdigest(), "secondes": round(time.time() - t0, 2), "essai": essai}
        except Exception as e:
            err = f"{type(e).__name__}: {str(e)[:150]}"
            time.sleep(3 * essai)
    return {"pdb": pdb, "statut": "echec", "url": url, "erreur": err}


def listify(x):
    return x if isinstance(x, list) else [x]


def struct_conn(pdb):
    d = MMCIF2Dict(str(CR / f"{pdb}.cif"))
    if "_struct_conn.id" not in d:
        return []
    n = len(listify(d["_struct_conn.id"]))
    g = lambda k: listify(d.get(k, ["?"] * n)) if k in d else ["?"] * n
    out = []
    for i in range(n):
        if g("_struct_conn.conn_type_id")[i] != "disulf":
            continue
        out.append({"id": g("_struct_conn.id")[i],
                    "p1": (g("_struct_conn.ptnr1_auth_asym_id")[i], g("_struct_conn.ptnr1_auth_seq_id")[i], g("_struct_conn.pdbx_ptnr1_PDB_ins_code")[i].replace("?", "").replace(".", ""), g("_struct_conn.ptnr1_label_comp_id")[i], g("_struct_conn.ptnr1_label_atom_id")[i]),
                    "p2": (g("_struct_conn.ptnr2_auth_asym_id")[i], g("_struct_conn.ptnr2_auth_seq_id")[i], g("_struct_conn.pdbx_ptnr2_PDB_ins_code")[i].replace("?", "").replace(".", ""), g("_struct_conn.ptnr2_label_comp_id")[i], g("_struct_conn.ptnr2_label_atom_id")[i]),
                    "sym1": g("_struct_conn.ptnr1_symmetry")[i], "sym2": g("_struct_conn.ptnr2_symmetry")[i], "dist_cristal_A": g("_struct_conn.pdbx_dist_value")[i]})
    return out


def residus_md(pdb):
    """{(chaine_origine, resSeq, icode): (chaine_x0, ordinal, nom)} depuis structure.pdb de la MD (replique 1), dans l ordre du fichier. Residus proteiques seulement."""
    j = json.loads((L.CUR / f"{pdb}__img0_rep1.json").read_text())
    table = j["table_chaines"]
    ordre, vus = {}, {}
    cles = []
    for l in open(C / pdb / "rep1/structure.pdb"):
        if not l.startswith("ATOM"):
            continue
        k = (l[21], int(l[22:26]), l[26].strip())
        if k not in vus:
            vus[k] = l[17:20]
            cles.append(k)
    cnt, out, md_ch = {}, {}, {}
    for k in cles:
        if k[0] not in table:
            continue
        cnt[k[0]] = cnt.get(k[0], 0) + 1
        out[k] = (table[k[0]], cnt[k[0]], vus[k])
        md_ch.setdefault(k[0], {})[(k[1], k[2])] = NORM.get(vus[k], vus[k])
    return out, table, md_ch


NORM = {"HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "CYX": "CYS", "MSE": "MET"}


def chaines_cristal(pdb):
    """{chaine_auth: {(auth_seq, icode): mon_id}} depuis _pdbx_poly_seq_scheme (residus observes)."""
    d = MMCIF2Dict(str(CR / f"{pdb}.cif"))
    out = {}
    if "_pdbx_poly_seq_scheme.pdb_strand_id" not in d:
        return out
    ch, seq, mon, ic = (listify(d[k]) for k in ("_pdbx_poly_seq_scheme.pdb_strand_id", "_pdbx_poly_seq_scheme.auth_seq_num", "_pdbx_poly_seq_scheme.pdb_mon_id", "_pdbx_poly_seq_scheme.pdb_ins_code"))
    for c, q, m, i in zip(ch, seq, mon, ic):
        if q in ("?", "."):
            continue
        out.setdefault(c, {})[(int(q), i.replace("?", "").replace(".", ""))] = NORM.get(m, m)
    return out


def table_chaines_contenu(cr_ch, md_ch):
    """crystal chain -> (MD chain_orig, identite, n_communs) par IDENTITE DE SEQUENCE sur les residus communs (>= 30 communs, identite >= 0,9) ; meilleur candidat,
    a egalite la meme lettre. Un copie de l unite asymetrique (autre lettre, meme sequence) s apparie a la meme chaine MD : c est declare `copie`."""
    out = {}
    for c, rs in cr_ch.items():
        meilleur = None
        for m, ms in md_ch.items():
            communs = [k for k in rs if k in ms]
            if len(communs) < 30:
                continue
            ident = sum(1 for k in communs if rs[k] == ms[k]) / len(communs)
            if ident >= 0.9 and (meilleur is None or (ident, m == c) > (meilleur[1], meilleur[0] == c)):
                meilleur = (m, ident, len(communs))
        out[c] = meilleur
    return out


def sg_md(pdb):
    """x0 (img0) des 3 repliques : {(chaine, ordinal): indice SG} et coordonnees ; + trajectoires via structure.pdb."""
    x = md.load(str(L.chemin_x0(pdb, 1)))
    sg = {}
    for k, c in enumerate(x.topology.chains):
        for o, r in enumerate(c.residues, start=1):
            if r.name == "CYS":
                sg[(chr(ord("A") + k), o)] = [a.index for a in r.atoms if a.name == "SG"][0]
    return sg


def distances_md(pdb, paires):
    """{paire: {image0_rep: [..], mediane_rep: [..]}} par ordinal ; les atomes viennent de structure.pdb de chaque replique (meme ordre que le PSF)."""
    out = {p: {"image0_A": [], "mediane_images_A": []} for p in paires}
    for rep in (1, 2, 3):
        d = C / pdb / f"rep{rep}"
        top = md.load(str(d / "structure.pdb")).topology
        sgi = {}
        for k, c in enumerate(top.chains):
            for o, r in enumerate(c.residues, start=1):
                if r.name in ("CYS", "CYX"):
                    s = [a.index for a in r.atoms if a.name == "SG"]
                    if s:
                        sgi[(chr(ord("A") + k), o)] = s[0]
        idx = sorted({sgi[a] for p in paires for a in p})
        if not idx:
            continue
        t = md.load(str(d / "trajectory.xtc"), top=str(d / "structure.pdb"), atom_indices=idx)
        pos = {a: i for i, a in enumerate(idx)}
        for p in paires:
            dist = np.linalg.norm(t.xyz[:, pos[sgi[p[0]]]] - t.xyz[:, pos[sgi[p[1]]]], axis=1) * 10
            out[p]["image0_A"].append(round(float(dist[0]), 2))
            out[p]["mediane_images_A"].append(round(float(np.median(dist)), 2))
    return out


for pdb in sorted(L.corpus_28()):
    log.append(telecharger(pdb))
    time.sleep(1.0)
print("TELECHARGEMENT", {s: sum(1 for e in log if e["statut"] == s) for s in {e["statut"] for e in log}}, flush=True)
(CR / "provenance.json").write_text(json.dumps(log, indent=1) + "\n")
lignes = []
for pdb in sorted(L.corpus_28()):
    rec = {"pdb": pdb, "paires": []}
    try:
        res, table, md_ch = residus_md(pdb)
        cr_ch = chaines_cristal(pdb) if (CR / f"{pdb}.cif").is_file() else {}
        tc = table_chaines_contenu(cr_ch, md_ch)
        rec["chaines_cristal_vers_md"] = {c: (None if v is None else {"chaine_md_origine": v[0], "identite": round(v[1], 3), "n_residus_communs": v[2], "meme_lettre": v[0] == c}) for c, v in tc.items()}
        sgx = sg_md(pdb)
        psf = {tuple(sorted([(p["chaine_i"], p["ordinal_i"]), (p["chaine_j"], p["ordinal_j"])])) for p in ref[pdb]["ponts"]}
        crist = {}
        non_app = []
        for b in (struct_conn(pdb) if (CR / f"{pdb}.cif").is_file() else []):
            a1, a2 = b["p1"], b["p2"]
            sym = b["sym1"] != b["sym2"] or b["sym1"] not in ("1_555", "?")
            m = []
            copie = False
            for a in (a1, a2):
                t = tc.get(a[0])
                if t is None:
                    m.append(None)
                    continue
                copie = copie or (t[0] != a[0])
                m.append(res.get((t[0], int(a[1]), a[2])) or res.get((t[0], int(a[1]), "")))
            if None in m:
                non_app.append({"cristal": f"{a1[0]}{a1[1]}{a1[2]}-{a2[0]}{a2[1]}{a2[2]}", "raison": ("chaine absente de la MD (aucune chaine MD de meme sequence : identite < 0,9 ou < 30 residus communs) ou residu non observe / hors modele"), "sym": [b["sym1"], b["sym2"]],
                                "dist_cristal_A": b["dist_cristal_A"]})
                continue
            if m[0][2] not in ("CYS", "CYX") or m[1][2] not in ("CYS", "CYX"):
                non_app.append({"cristal": f"{a1[0]}{a1[1]}-{a2[0]}{a2[1]}", "raison": f"residus MD {m[0][2]},{m[1][2]} ne sont pas des cysteines (appariement faux)", "dist_cristal_A": b["dist_cristal_A"]})
                continue
            cle = tuple(sorted([(m[0][0], m[0][1]), (m[1][0], m[1][1])]))
            if cle in crist:
                crist[cle]["copies"].append({"auth": f"{a1[0]}{a1[1]}{a1[2]}-{a2[0]}{a2[1]}{a2[2]}", "dist_cristal_A": b["dist_cristal_A"], "chaine_differente_de_la_md": copie})
            else:
                crist[cle] = {"auth": f"{a1[0]}{a1[1]}{a1[2]}-{a2[0]}{a2[1]}{a2[2]}", "dist_cristal_A": b["dist_cristal_A"], "symetrie_non_identite": bool(sym), "apparie_par_contenu_chaine_renommee": bool(copie), "copies": []}
        # cysteines libres proches dans l image 0 de la MD (3 repliques)
        proches = set()
        lies = {a for p in psf for a in p} | {a for p in crist for a in p}
        libres = [k for k in sgx if k not in lies]
        for rep in (1, 2, 3):
            x = md.load(str(L.chemin_x0(pdb, rep)))
            s = {}
            for k, c in enumerate(x.topology.chains):
                for o, r in enumerate(c.residues, start=1):
                    if r.name == "CYS":
                        s[(chr(ord("A") + k), o)] = [a.index for a in r.atoms if a.name == "SG"][0]
            for i in range(len(libres)):
                for j in range(i + 1, len(libres)):
                    if np.linalg.norm(x.xyz[0, s[libres[i]]] - x.xyz[0, s[libres[j]]]) * 10 < 3.0:
                        proches.add(tuple(sorted([libres[i], libres[j]])))
        toutes = sorted(set(crist) | psf | proches)
        dm = distances_md(pdb, toutes) if toutes else {}
        for p in toutes:
            c, s_, pr = p in crist, p in psf, p in proches
            d = dm[p]
            rec["paires"].append({"paire": f"{p[0][0]}{p[0][1]}-{p[1][0]}{p[1][1]}", "cristal": c, "auth_cristal": crist[p]["auth"] if c else None,
                                  "dist_cristal_A": crist[p]["dist_cristal_A"] if c else None, "symetrie_non_identite": crist[p]["symetrie_non_identite"] if c else None, "apparie_par_contenu_chaine_renommee": crist[p]["apparie_par_contenu_chaine_renommee"] if c else None, "copies_cristal": crist[p]["copies"] if c else [],
                                  "psf": s_, "md_proche_sans_pont_psf": pr and not s_, "accord_cristal_psf": c == s_,
                                  "md_image0_rep1_2_3_A": d["image0_A"], "md_mediane_images_rep1_2_3_A": d["mediane_images_A"]})
        rec["non_apparies"] = non_app
        rec["n_cristal"], rec["n_psf"], rec["n_desaccords"] = len(crist), len(psf), sum(1 for q in rec["paires"] if not q["accord_cristal_psf"])
        rec["statut"] = "ok"
    except Exception as e:
        rec.update(statut="ecarte", raison=f"{type(e).__name__}: {str(e)[:300]}")
    tout[pdb] = rec
    print(f"CRISTAL_PSF {pdb} {rec['statut']} cristal={rec.get('n_cristal')} psf={rec.get('n_psf')} desaccords={rec.get('n_desaccords')} non_apparies={len(rec.get('non_apparies', []))} {rec.get('raison', '')}", flush=True)
    if rec["statut"] == "ok" and rec["n_desaccords"]:
        for q in rec["paires"]:
            if not q["accord_cristal_psf"]:
                print(f"   DESACCORD {pdb} {q['paire']} cristal={q['cristal']} ({q['auth_cristal']}, {q['dist_cristal_A']} A) psf={q['psf']} md_image0={q['md_image0_rep1_2_3_A']} md_mediane={q['md_mediane_images_rep1_2_3_A']}", flush=True)
(CR / "croisement_cristal_psf.json").write_text(json.dumps(tout, indent=1, ensure_ascii=False) + "\n")
dis = sorted(p for p, r in tout.items() if r.get("n_desaccords"))
ec = sorted(p for p, r in tout.items() if r["statut"] != "ok")
print(f"CRISTAL_PSF_BILAN {len(tout)} systemes, ecartes {len(ec)} {ec}, systemes avec desaccord cristal/PSF {len(dis)} {dis} (rejets = {len(ec)})")
sys.exit(0 if not ec else 1)
