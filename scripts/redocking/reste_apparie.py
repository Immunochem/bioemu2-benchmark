import json
R = "/gpfs/phherman/bioemu2_bench/redocking"
src = {t["id"] for t in json.load(open(R + "/taches_sourcesX.json"))}
T = [t for t in json.load(open(R + "/taches_A_apparie_reste.json")) if t["id"] not in src]
json.dump(T, open(R + "/taches_A_apparie_reste2.json", "w"), indent=1)
print(len(T), "A_apparie restants (hors sources de NEG_A)", [t["id"] for t in T])
