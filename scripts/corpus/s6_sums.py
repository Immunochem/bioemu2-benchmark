"""Ecrit corpus/SHA256SUMS depuis la derniere ligne du manifeste de chaque systeme (chemins relatifs a corpus/)."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from corpus_lib import ROOT  # noqa: E402

last = {}
n_lignes = 0
for line in open(ROOT + "/manifeste_corpus.json", encoding="utf-8"):
    r = json.loads(line)
    n_lignes += 1
    last[r["pdb"] + "/" + r["accession"]] = r
n = 0
with open(ROOT + "/SHA256SUMS", "w", encoding="utf-8") as f:
    for k in sorted(last):
        for fi in last[k]["fichiers"]:
            f.write("%s  %s\n" % (fi["sha256"], os.path.relpath(fi["chemin"], ROOT)))
            n += 1
print("manifeste: %d lignes, %d systemes, %d fichiers ecrits dans SHA256SUMS" % (n_lignes, len(last), n))
