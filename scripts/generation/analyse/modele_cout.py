"""modele_cout.py — cout PREDIT (secondes) d'un systeme pour une vague, depuis les MESURES de capacite du manifeste.
Aucun chiffre n'est tape ici : s/structure par (tokens, batch) = dernier essai ok a concurrence 4 ; temps d'embedding =
premier appel (n=1, bs1) moins le s/structure bs1 du meme systeme. Entre deux tailles mesurees : interpolation log-log ;
au-dela : extrapolation DECLAREE avec l'exposant du segment le plus proche (au moins 1,0 vers le bas : on surestime les petits)."""
import json, math, sys
import chargerman as C

def mesures(L, conc_min=3):
    s, emb, tok = {}, {}, {}
    for d in L:
        if d.get("type") != "essai_capacite":
            continue
        N, bs = d["tokens"], d["batch_size"]
        c = d["jobs_concurrents_gpu"]
        if d["statut"] == "ok" and d.get("debit_s_par_structure") is not None and (c["min"] or 0) >= conc_min:
            s.setdefault(bs, {})[N] = d["debit_s_par_structure"]
        if d["statut"] == "ok_jete_pour_le_debit":
            emb.setdefault(N, []).append(d["temps_s"])
    return s, emb

def interp(points, N):
    """points : {N: y} ; log-log. Entre deux mesures : interpolation du segment. HORS bornes : exposant de la regression
    log-log sur TOUS les points (>= 3) — jamais l exposant local d un segment, trop bruite (2,54 mesure entre 588 et 643 a bs 8) —
    et au moins 1,0 vers le bas (on surestime les petits systemes plutot que de les sous-estimer)."""
    xs = sorted(points)
    if N in points:
        return points[N], "mesure"
    glob = C.regression_loglog(xs, [points[x] for x in xs]) if len(xs) >= 3 else None
    for a, b in zip(xs, xs[1:]):
        if a < N < b:
            e = math.log(points[b] / points[a]) / math.log(b / a)
            return points[a] * (N / a) ** e, "interpole"
    e = glob[0] if glob else 1.0
    if N < xs[0]:
        return points[xs[0]] * (N / xs[0]) ** max(1.0, e), "extrapole_vers_le_bas_exposant_global"
    return points[xs[-1]] * (N / xs[-1]) ** e, "extrapole_vers_le_haut_exposant_global"


def cout(L, N, usages, bs, n=25, overhead_unite_s=15.0, lags=5):
    s, emb, _ = mesures(L), None, None
    S, E = s
    sb = lambda b: interp(S[b], N)
    t8, o8 = sb(bs)
    t1, _ = sb(1)
    plein, reste = divmod(n, bs)
    t_unit = plein * bs * t8 + reste * (t1 if reste == 1 else t8 * reste)
    # embedding : premier appel - s/struct bs1 (meme systeme) -> points {N: t}
    pts = {}
    for k, v in E.items():
        if k in S.get(1, {}):
            pts[k] = min(v) - S[1][k]
    te, oe = interp(pts, N)
    nunits = (lags if "U3" in usages else 0) + (1 if "U2" in usages else 0) + (1 if "U1" in usages else 0)
    nemb = (1 if ("U3" in usages or "U2" in usages) else 0) + (1 if "U1" in usages else 0)
    return nemb * te + nunits * (t_unit + overhead_unite_s), {"s_bs": round(t8, 2), "origine_s": o8, "t_emb": round(te, 1), "origine_emb": oe}

if __name__ == "__main__":
    L = C.lignes()
    S, E = mesures(L)
    print("s/struct mesures (c>=3) :", {b: dict(sorted(v.items())) for b, v in sorted(S.items())})
