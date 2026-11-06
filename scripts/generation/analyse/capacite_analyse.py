"""capacite_analyse.py — courbe de capacite (exposants alloue / reserve), echelle de batch, concurrence, depuis le manifeste."""
import json, math
import chargerman as C

def analyse(L):
    cap = [d for d in L if d.get("type") == "essai_capacite"]
    sys_ = {}
    for d in cap:
        N = d["tokens"]
        s = sys_.setdefault(d["systeme"], {"tokens": N, "premier": None, "bs": {}})
        if d["statut"] == "ok_jete_pour_le_debit":
            s["premier"] = d
        elif d["statut"] == "ok":
            c = d["jobs_concurrents_gpu"]
            s["bs"].setdefault(d["batch_size"], []).append(d)
    return sys_

def lignes_tableau(sys_):
    out = []
    for p, s in sorted(sys_.items(), key=lambda kv: kv[1]["tokens"]):
        pr = s["premier"]
        b1 = s["bs"].get(1, [None])[-1]
        row = {"pdb": p, "tokens": s["tokens"],
               "emb_alloue": pr["pic_alloue_Gio"] if pr else None, "emb_reserve": pr["pic_reserve_Gio"] if pr else None,
               "emb_ligne": pr["_ligne"] if pr else None, "emb_t_s": pr["temps_s"] if pr else None,
               "deb_alloue_bs1": b1["pic_alloue_Gio"] if b1 else None, "deb_reserve_vide_bs1": b1["pic_reserve_apres_vidage_Gio"] if b1 else None,
               "deb_reserve_continu_bs1": b1["pic_reserve_Gio"] if b1 else None, "ligne_bs1": b1["_ligne"] if b1 else None}
        for b in (1, 2, 4, 5, 8):
            e = s["bs"].get(b)
            row[f"s_bs{b}"] = e[-1]["debit_s_par_structure"] if e else None
            row[f"l_bs{b}"] = e[-1]["_ligne"] if e else None
            row[f"conc_bs{b}"] = f"{e[-1]['jobs_concurrents_gpu']['min']}-{e[-1]['jobs_concurrents_gpu']['max']}" if e else None
            row[f"alloue_bs{b}"] = e[-1]["pic_alloue_Gio"] if e else None
        out.append(row)
    return out

def exposants(rows):
    res = {}
    for nom, cle in (("embedding_alloue", "emb_alloue"), ("embedding_reserve", "emb_reserve"), ("debruitage_alloue_bs1", "deb_alloue_bs1"),
                     ("debruitage_reserve_apres_vidage_bs1", "deb_reserve_vide_bs1"), ("debruitage_reserve_continu_bs1", "deb_reserve_continu_bs1"),
                     ("s_par_structure_bs1", "s_bs1"), ("s_par_structure_bs4", "s_bs4"), ("s_par_structure_bs8", "s_bs8"), ("embedding_temps_s", "emb_t_s")):
        pts = [(r["tokens"], r[cle]) for r in rows if r[cle]]
        r = C.regression_loglog([p[0] for p in pts], [p[1] for p in pts])
        res[nom] = None if r is None else {"exposant": round(r[0], 3), "erreur_type": round(r[2], 3), "n_points": r[3], "coef": round(math.exp(r[1]), 6)}
    return res

if __name__ == "__main__":
    L = C.lignes()
    S = analyse(L)
    R = lignes_tableau(S)
    for r in R:
        print({k: v for k, v in r.items() if v is not None and not k.startswith("l_") and not k.startswith("conc_")})
    print(json.dumps(exposants(R), indent=1))
