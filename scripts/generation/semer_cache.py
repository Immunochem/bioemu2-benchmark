#!/usr/bin/env python3
"""semer_cache.py — AMORCE generation/msa_cache par COPIE des caches ColabFold existants (lecture seule sur les sources).
Sources : msa/cache (lot MSA du socle), 6CWG/msa_cache (pilote), vram_1030/msa_cache. Une cle est le sha256 de la requete
(bioemu2/msa/colabfold.py:cached_colabfold_msa) : deux fichiers de meme nom DOIVENT etre identiques ; sinon c est un CONFLIT,
rapporte, jamais ecrase. Chaque copie est relue par sha256."""
import hashlib, json, shutil, sys
from pathlib import Path
B = Path("/gpfs/phherman/bioemu2_bench")
SRC = [B / "msa/cache", B / "6CWG/msa_cache", B / "vram_1030/msa_cache"]
DST = B / "generation/msa_cache"
DST.mkdir(parents=True, exist_ok=True)
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
bilan = {"copies": 0, "deja_present_identique": 0, "conflits": [], "sources": {}}
for s in SRC:
    n = 0
    for f in sorted(s.glob("*.json")):
        n += 1
        d = DST / f.name
        if d.exists():
            if sha(d) == sha(f):
                bilan["deja_present_identique"] += 1
            else:
                bilan["conflits"].append({"fichier": f.name, "source": str(s), "sha_dst": sha(d)[:16], "sha_src": sha(f)[:16]})
            continue
        shutil.copy2(f, d)
        if sha(d) != sha(f):
            bilan["conflits"].append({"fichier": f.name, "source": str(s), "erreur": "copie differente de la source"})
            d.unlink()
            continue
        bilan["copies"] += 1
    bilan["sources"][str(s)] = n
bilan["n_dans_la_destination"] = len(list(DST.glob("*.json")))
print("SEMER", json.dumps(bilan, ensure_ascii=False))
(B / "generation/msa_cache_amorcage.json").write_text(json.dumps(bilan, indent=2, ensure_ascii=False) + "\n")
sys.exit(1 if bilan["conflits"] else 0)
