#!/usr/bin/env python3
"""tables_rapport.py — produit les TABLEAUX du rapport depuis les JSON (aucun chiffre recopie a la main).
Lit : resultats/{reproduction_v1,proxy_cb,5lzm_inspection,analyse_unites,rmsd_ferme_unites}.json  Ecrit : resultats/tables_rapport.md"""
import json, os, sys
T = "/gpfs/phherman/bioemu2_bench/t4l"
R = f"{T}/resultats"
def J(n):
    p = f"{R}/{n}.json"
    return json.load(open(p)) if os.path.isfile(p) else None
RAC = sys.argv[1] if len(sys.argv) > 1 else "unites"
rep, prox, cris, ana, rms = J("reproduction_v1"), J("proxy_cb"), J("5lzm_inspection"), J(f"analyse_{RAC}"), J(f"rmsd_ferme_{RAC}")
if RAC != "unites":
    rep = prox = cris = None
L = []
def w(s=""):
    L.append(s)
f = lambda x, n=3: "—" if x is None else (f"{x:.{n}f}" if isinstance(x, float) else str(x))

if rep:
    w("## T1 · reproduction : MES d1, d2, p contre SES tableaux (`resultats/reproduction_v1.json`)")
    w("| méthode | n images | atomes | max abs écart d1 (nm) | max abs écart d2 (nm) | p défini | max abs écart p | d1 < 2,5 nm | p>0 / p<0 (SES locked.pkl) |")
    w("|---|---|---|---|---|---|---|---|---|")
    for m, e in rep["methodes"].items():
        pl = e.get("p_a_lui", {})
        w(f"| {m} | {e['n_frames_xtc']} | {e['n_atomes']} | {e['ecart_max_d1_nm']:.1e} | {e['ecart_max_d2_nm']:.1e} | {e['p_defini']} | {('%.1e' % e['ecart_max_p']) if e.get('ecart_max_p') is not None else 'non défini'} | {e['d1_ferme_lt_2_5']} | {str(pl.get('p_gt_0', '—')) + ' / ' + str(pl.get('p_lt_0', '—'))} |")
    w(f"\nRejets de ce calcul : {rep['n_rejets']}.")
    w("\n## T2 · vrais positifs : les 221 structures « fermées » de l'auteur")
    w("| jeu | fichiers | rejets | d1 < 2,5 nm | d1 min / méd / max | bon signe de p | p min / méd / max |")
    w("|---|---|---|---|---|---|---|")
    for k, v in rep["C3"].items():
        w(f"| {k} | {v['n_fichiers']} | {v['n_rejets']} | {v['d1_lt_2_5']} | {' / '.join(f(x) for x in v['d1_min_med_max'])} | {v['p_du_bon_signe']} | {' / '.join(f(x) for x in v['p_min_med_max'])} |")
    w("\n## T3 · son camembert « Gen-AI » : trois façons de le calculer sur les MÊMES structures")
    c4 = rep["C4"]
    a = c4["a_sa_cellule_48_lock_egal_colonne_1"]
    w(f"- (a) exactement sa cellule [48] (`lock = col1(all_genai)` = colonne 1 de `distances.npy` = d2) : n = {a['n']}, comptes [exp/closed, exp/open, bur/closed, bur/open] = {a['comptes_[exp/closed, exp/open, bur/closed, bur/open]']} ; colonne 1 min = {a['colonne_1_min']:.3f} nm, toujours positive : {a['colonne_1_est_toujours_positive']}.")
    for k, lab in (("b_vrai_p_les_cinq", "(b) SON vrai p (locked.pkl), les cinq méthodes de sa cellule"), ("b_vrai_p_sept_methodes", "(b') SON vrai p, les sept méthodes à cycle complet"), ("c_mon_p_sept_methodes", "(c) MON p, les sept méthodes")):
        v = c4[k]; w(f"- {lab} : n = {v['n_total']} ; " + " · ".join(f"{n} {v['comptes'][n]} ({v['pourcent_sur_definis'][n]:.1f} %)" for n in v["comptes"]) + f" ; indéfinis {v['n_indefini']}.")
if prox:
    w("\n## T4 · proxy CB pour BioEmu v1 (p du papier non défini : pas de cycle de Phe4, `topology.pdb` squelette + CB)")
    w(f"- calibration (n = {prox['n_calibration']}) : accord de signe {prox['accord_de_signe']:.3f}, Pearson {prox['pearson']:.3f}, pente {prox['pente_moindres_carres']:.2f} ; critère écrit avant (≥ 0,90) : {prox['critere_ecrit_avant_accord_ge_0_90']}")
    w("- accord de signe par jeu : " + " · ".join(f"{k} {v['accord_de_signe']:.3f} (n={v['n']})" for k, v in prox["par_jeu"].items()))
    b = prox["bioemu_v1_proxy"]
    w(f"- BioEmu v1 (n = {b['n']}), p_proxy : >0 {b['p_gt_0']} · <0 {b['p_lt_0']} · médiane {b['p_min_med_max'][1]:.3f} · fermés (d1) {b['d1_ferme']}")
if cris:
    w(f"\n## T5 · 5LZM : ATOM {cris['n_ATOM']}, HETATM {cris['n_HETATM']} {cris['HETATM_par_residu']}, altLoc {cris['n_altloc']}, résidus à atomes manquants {len(cris['residus_atomes_manquants'])} ; d1 = {cris['d1_nm']:.4f} nm, d2 = {cris['d2_nm']:.4f} nm, p = {cris['p']:.4f}")
if ana:
    w("\n## T6 · populations des quatre états, ses CV et ses seuils (`resultats/analyse_unites.json`)")
    w("| unité | système | lag | n | physiques | rejetées | exp/open | exp/closed | bur/open | bur/closed | indéfinis | s/struct |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for u, e in ana["unites"].items():
        t = e["toutes"]; x = t["etats"]
        w(f"| {u} | {e['systeme']} | {e['lag_ns'] if e['lag_ns'] is not None else '—'} | {e['n']} | {e['physiques']} | {e['rejetees_physicalite']} | " + " | ".join(f"{x[n]['n']} ({x[n]['pct']:.1f} %) [{x[n]['IC95_wilson_pct'][0]}–{x[n]['IC95_wilson_pct'][1]}]" for n in ("exposed/open", "exposed/closed", "buried/open", "buried/closed")) + f" | {t['n_indefini']} | {e['s_par_structure']} |")
    w("\n## T7 · distribution de d1 (nm) et distance au seuil 2,5 nm")
    w("| unité | min | 5 % | 25 % | médiane | 75 % | 95 % | max | d1<2,5 | d1<2,4 | dans [2,4;2,6) | plus proche du seuil (nm) | fermé sous bruit ±0,1 nm : observé % → [2,5 ; 50 ; 97,5] % |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for u, e in ana["unites"].items():
        d = e["d1"]; q = d["quantiles_0_5_25_50_75_95_100_nm"]; b = e["bruit_d1"]
        w(f"| {u} | " + " | ".join(f"{x:.3f}" for x in q) + f" | {d['ferme_lt_2_5']} | {d['ferme_robuste_lt_2_4']} | {d['ambigu_2_4_a_2_6']} | {d['plus_proche_du_seuil_nm']} | {b['ferme_observe_pct']} → {b['ferme_sous_bruit_pct_2_5_50_97_5']} |")
    au = ana["auteur"]
    for nm, d in (("GenAI (7 méthodes), son vrai p", au["GenAI_sept_methodes_vrai_p"]["d1"]), ("BioEmu v1", au["BioEmu_v1"]["d1"])):
        q = d["quantiles_0_5_25_50_75_95_100_nm"]
        w(f"| auteur : {nm} | " + " | ".join(f"{x:.3f}" for x in q) + f" | {d['ferme_lt_2_5']} | {d['ferme_robuste_lt_2_4']} | {d['ambigu_2_4_a_2_6']} | {d['plus_proche_du_seuil_nm']} | — |")
    w("\n## T8 · couverture (grilles de 0,25 nm ; p par pas de 0,25)")
    w("| unité | cases (d1,d2) occupées | entropie norm. d1 | entropie norm. p | d1 min–max | d2 min–max | p min / méd / max | \\|p\\| < 0,05 |")
    w("|---|---|---|---|---|---|---|---|")
    for u, e in ana["unites"].items():
        c = e["couverture"]
        w(f"| {u} | {c['cases_d1_d2_occupees_(0,25nm)']} | {c['entropie_norm_d1']} | {c['entropie_norm_p']} | {c['d1_min_max'][0]}–{c['d1_min_max'][1]} | {c['d2_min_max'][0]}–{c['d2_min_max'][1]} | {' / '.join(str(x) for x in c['p_min_med_max'])} | {c['p_proche_de_0_abs_lt_0_05']} |")
    for nm, c in (("auteur : GenAI (7), son vrai p", au["GenAI_sept_methodes_vrai_p"]["couverture"]), ("auteur : BioEmu v1 (p non défini)", au["BioEmu_v1"]["couverture_d1_d2"])):
        w(f"| {nm} | {c['cases_d1_d2_occupees_(0,25nm)']} | {c['entropie_norm_d1']} | {c['entropie_norm_p']} | {c['d1_min_max'][0]}–{c['d1_min_max'][1]} | {c['d2_min_max'][0]}–{c['d2_min_max'][1]} | {c['p_min_med_max']} | — |")
    w("\n## T9 · comparaison de d1 (KS à deux échantillons, Wasserstein en nm)")
    for k, v in ana["comparaison_d1"].items():
        w(f"- {k} : KS = {v['KS']}, p = {v['p']:.2g}, Wasserstein = {v['wasserstein_nm']} nm")
if rms:
    w("\n## T10 · ressemblance aux 221 structures fermées de l'auteur (RMSD Cα, 162 résidus, appariement ordinal)")
    w(f"- étalon (plus proche référence d'un AUTRE cluster) : médiane / 95e / max = {rms['etalon_autre_cluster_A_med_95_max']} Å ; références → cristal 5LZM : min / méd / max = {rms['refs_vers_cristal_5LZM_A_min_med_max']} Å")
    for u, e in rms["unites"].items():
        s = f"- {u} : n = {e['n']}, d1 < 2,5 nm : {e['n_ferme_d1']} ; RMSD min aux références (0/25/50/75/100 %) = {e['rmsd_min_aux_refs_fermees_A_quantiles_0_25_50_75_100_toutes']} Å ; RMSD au cristal (médiane) = {e['rmsd_au_cristal_A_mediane_toutes']} Å"
        if "fermees" in e:
            s += f" ; FERMÉES : {e['fermees']}"
        w(s)
open(f"{R}/tables_rapport.md" if RAC == "unites" else f"{R}/tables_rapport_{RAC}.md", "w").write("\n".join(L) + "\n")
print("\n".join(L)[:200], "...", len(L), "lignes")
