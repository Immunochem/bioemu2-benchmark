import json, sys
lag = sys.argv[1]; conds = sys.argv[2].split(",")
T = json.load(open(f"/gpfs/phherman/bioemu2_bench/redocking/resultats/tableaux_lag{int(lag):04d}.json"))
for k, c in T["cellules"].items():
    s, cond = k.split("|")
    if cond not in conds: continue
    print(k, "images>=acc", f'{c["images_acceptable_ou_mieux"]["k"]}/{c["images_acceptable_ou_mieux"]["n"]}', "| par systeme (fnat med, iRMSD med, classe):",
          " ".join(f'{sid}:{v["fnat"]:.2f}/{v["irmsd"]:.1f}/{v["classe"][:3]}' for sid, v in sorted(c["par_systeme"].items())))
