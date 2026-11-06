import json, re
T = "/gpfs/phherman/bioemu2_bench/t4l"
r = open(f"{T}/RAPPORT_T4L.md").read()
a = json.load(open(f"{T}/resultats/analyse_unites.json")); s = json.load(open(f"{T}/resultats/analyse_unites_steer.json"))
rm = json.load(open(f"{T}/resultats/rmsd_ferme_unites_steer.json")); rc = json.load(open(f"{T}/resultats/rmsd_cristal.json"))
fr = lambda x, n=3: f"{x:.{n}f}".replace(".", ",")
chk = []
e = a["unites"]["eq164_n500"]
chk += [("eq164 exp/open 460", e["toutes"]["etats"]["exposed/open"]["n"] == 460 and "| exposed/open | 460 |" in r), ("eq164 bur/open 40", e["toutes"]["etats"]["buried/open"]["n"] == 40 and "| buried/open | 40 |" in r),
        ("d1 min 3,524", fr(e["d1"]["quantiles_0_5_25_50_75_95_100_nm"][0]) in r), ("d1 med 3,637", fr(e["d1"]["quantiles_0_5_25_50_75_95_100_nm"][3]) in r),
        ("plus proche seuil 1,024", fr(e["d1"]["plus_proche_du_seuil_nm"]) in r), ("cases 8", e["couverture"]["cases_d1_d2_occupees_(0,25nm)"] == 8),
        ("entropie d1 0,033", fr(e["couverture"]["entropie_norm_d1"]) in r)]
for u, k in (("prop162_lag1000", 31), ("prop162_lag0001", 0), ("eq162_n250", 16)):
    chk.append((f"{u} buried {k}", a["unites"][u]["toutes"]["etats"]["buried/open"]["n"] == k))
chk.append(("rmsd cristal v1 med 1,16", rc["BioEmu_v1"]["rmsd_cristal_A_min_5_50_95_max"][2] == 1.16 and "1,16 (2,82 ; 8,57)" in r))
chk.append(("rmsd cristal eq164 med 0,76", rc["BioEmu2:eq164_n500"]["rmsd_cristal_A_min_5_50_95_max"][2] == 0.76 and "0,76 (1,05 ; 3,03)" in r))
for u, ph in (("st_f000_inerte", 100), ("st_d1_f005", 100), ("st_d1_f020", 53), ("st_d1d2_f020", 0)):
    chk.append((f"{u} physiques {ph}", s["unites"][u]["physiques"] == ph and f"{ph}/100" in r))
chk.append(("steer f020 d1 med 2,362", fr(s["unites"]["st_d1_f020"]["d1"]["quantiles_0_5_25_50_75_95_100_nm"][3]) in r))
fp = rm["unites"]["st_d1_f020"]["fermees_ET_physiques"]
chk.append(("53 physiques RMSD refs 5,72/5,98/6,40", fp["n"] == 53 and fp["rmsd_min_aux_refs_A_min_med_max"] == [5.72, 5.98, 6.4] and "5,72 – 5,98 (médiane) – 6,40 Å" in r))
chk.append(("53 sous etalon mediane = 0", fp["n_sous_etalon_mediane"] == 0))
chk.append(("{{ restant", "{{" not in r)); chk.append(("ABSENT", "ABSENT" not in r))
bad = [n for n, ok in chk if not ok]
for n, ok in chk: print("OK " if ok else "KO ", n)
print("ECARTS :", len(bad), bad)
