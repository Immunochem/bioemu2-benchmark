import json, numpy as np
from scipy import stats
T = "/gpfs/phherman/bioemu2_bench/t4l"
U = {u: np.load(f"{T}/unites/{u}/cv.npz") for u in ("eq164_n500", "eq162_n250", "prop162_lag0001", "prop162_lag0050", "prop162_lag1000")}
paires = [("eq164_n500", "eq162_n250"), ("eq162_n250", "prop162_lag0001"), ("eq162_n250", "prop162_lag0050"), ("eq162_n250", "prop162_lag1000"), ("prop162_lag0050", "prop162_lag1000"), ("prop162_lag0001", "prop162_lag0050")]
out = {}
for a, b in paires:
    for var in ("d1", "d2", "p"):
        r = stats.ks_2samp(U[a][var], U[b][var]); out[f"{a}|{b}|{var}"] = {"KS": round(float(r.statistic), 3), "p": float(r.pvalue)}
    # fraction enfouie (p<0) et test exact de Fisher
    ka, na = int((U[a]["p"] < 0).sum()), len(U[a]["p"]); kb, nb = int((U[b]["p"] < 0).sum()), len(U[b]["p"])
    out[f"{a}|{b}|fisher_enfoui"] = {"enfoui_a": f"{ka}/{na}", "enfoui_b": f"{kb}/{nb}", "p_fisher": float(stats.fisher_exact([[ka, na - ka], [kb, nb - kb]])[1])}
# p par unite : quantiles
for u, z in U.items():
    out[f"{u}|p_quantiles_0_5_25_50_75_95_100"] = [round(float(x), 3) for x in np.percentile(z["p"], [0, 5, 25, 50, 75, 95, 100])]
    out[f"{u}|d2_quantiles_0_5_50_95_100"] = [round(float(x), 3) for x in np.percentile(z["d2"], [0, 5, 50, 95, 100])]
json.dump(out, open(f"{T}/resultats/ks_unites.json", "w"), indent=1)
for k, v in out.items(): print(k, v)
