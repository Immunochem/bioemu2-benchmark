## T1 · reproduction : MES d1, d2, p contre SES tableaux (`resultats/reproduction_v1.json`)
| méthode | n images | atomes | max abs écart d1 (nm) | max abs écart d2 (nm) | p défini | max abs écart p | d1 < 2,5 nm | p>0 / p<0 (SES locked.pkl) |
|---|---|---|---|---|---|---|---|---|
| AlphaFlow | 500 | 1308 | 3.4e-07 | 1.1e-07 | True | 9.1e-08 | 0 | 494 / 6 |
| ESMFlow | 100 | 1308 | 2.5e-07 | 1.3e-07 | True | 8.5e-08 | 0 | 71 / 29 |
| AF-cluster | 132 | 1308 | 3.2e-07 | 3.7e-07 | True | 1.7e-07 | 0 | 41 / 91 |
| ConforFold | 600 | 1308 | 3.3e-07 | 1.6e-07 | True | 2.0e-07 | 0 | 179 / 421 |
| ConfRover | 500 | 1308 | 5.0e-07 | 3.0e-07 | True | 1.5e-07 | 0 | 8 / 492 |
| rMSA-AF2 | 80 | 1308 | 3.0e-07 | 1.3e-07 | True | 1.3e-07 | 0 | 26 / 54 |
| rMSA-AF2-recyc | 285 | 1308 | 3.0e-07 | 1.4e-07 | True | 1.5e-07 | 0 | 122 / 163 |
| BioEmu | 499 | 809 | 4.1e-07 | 1.8e-07 | False | non défini | 0 | — / — |

Rejets de ce calcul : 0.

## T2 · vrais positifs : les 221 structures « fermées » de l'auteur
| jeu | fichiers | rejets | d1 < 2,5 nm | d1 min / méd / max | bon signe de p | p min / méd / max |
|---|---|---|---|---|---|---|
| Exposed-closed | 118 | 0 | 118 | 1.988 / 2.348 / 2.498 | 118 | 0.000 / 0.278 / 0.946 |
| Buried-closed | 103 | 0 | 103 | 1.966 / 2.392 / 2.499 | 103 | -2.044 / -0.566 / -0.019 |

## T3 · son camembert « Gen-AI » : trois façons de le calculer sur les MÊMES structures
- (a) exactement sa cellule [48] (`lock = col1(all_genai)` = colonne 1 de `distances.npy` = d2) : n = 1832, comptes [exp/closed, exp/open, bur/closed, bur/open] = [0, 1832, 0, 0] ; colonne 1 min = 0.636 nm, toujours positive : True.
- (b) SON vrai p (locked.pkl), les cinq méthodes de sa cellule : n = 1832 ; exposed/open 793 (43.3 %) · exposed/closed 0 (0.0 %) · buried/open 1039 (56.7 %) · buried/closed 0 (0.0 %) ; indéfinis 0.
- (b') SON vrai p, les sept méthodes à cycle complet : n = 2197 ; exposed/open 941 (42.8 %) · exposed/closed 0 (0.0 %) · buried/open 1256 (57.2 %) · buried/closed 0 (0.0 %) ; indéfinis 0.
- (c) MON p, les sept méthodes : n = 2197 ; exposed/open 941 (42.8 %) · exposed/closed 0 (0.0 %) · buried/open 1256 (57.2 %) · buried/closed 0 (0.0 %) ; indéfinis 0.

## T4 · proxy CB pour BioEmu v1 (p du papier non défini : pas de cycle de Phe4, `topology.pdb` squelette + CB)
- calibration (n = 2418) : accord de signe 0.926, Pearson 0.902, pente 0.65 ; critère écrit avant (≥ 0,90) : True
- accord de signe par jeu : AlphaFlow 0.950 (n=500) · ESMFlow 0.890 (n=100) · AF-cluster 0.932 (n=132) · ConforFold 0.977 (n=600) · ConfRover 0.932 (n=500) · rMSA-AF2 0.863 (n=80) · rMSA-AF2-recyc 0.860 (n=285) · Exposed-closed 0.941 (n=118) · Buried-closed 0.718 (n=103)
- BioEmu v1 (n = 499), p_proxy : >0 166 · <0 333 · médiane -0.039 · fermés (d1) 0

## T5 · 5LZM : ATOM 1292, HETATM 133 {' CL': 2, 'BME': 12, 'HOH': 119}, altLoc 0, résidus à atomes manquants 0 ; d1 = 3.5593 nm, d2 = 0.8736 nm, p = 0.2410

## T6 · populations des quatre états, ses CV et ses seuils (`resultats/analyse_unites.json`)
| unité | système | lag | n | physiques | rejetées | exp/open | exp/closed | bur/open | bur/closed | indéfinis | s/struct |
|---|---|---|---|---|---|---|---|---|---|---|---|
| eq162_n250 | s162 | — | 250 | 250 | 0 | 234 (93.6 %) [89.86–96.02] | 0 (0.0 %) [0–1.51] | 16 (6.4 %) [3.98–10.14] | 0 (0.0 %) [0–1.51] | 0 | 4.358 |
| eq164_n500 | s164 | — | 500 | 500 | 0 | 460 (92.0 %) [89.29–94.07] | 0 (0.0 %) [0–0.76] | 40 (8.0 %) [5.93–10.71] | 0 (0.0 %) [0–0.76] | 0 | 4.579 |
| prop162_lag0001 | s162 | 1 | 200 | 200 | 0 | 200 (100.0 %) [98.12–100] | 0 (0.0 %) [0–1.88] | 0 (0.0 %) [0–1.88] | 0 (0.0 %) [0–1.88] | 0 | 4.522 |
| prop162_lag0050 | s162 | 50 | 200 | 200 | 0 | 200 (100.0 %) [98.12–100] | 0 (0.0 %) [0–1.88] | 0 (0.0 %) [0–1.88] | 0 (0.0 %) [0–1.88] | 0 | 4.358 |
| prop162_lag1000 | s162 | 1000 | 200 | 200 | 0 | 169 (84.5 %) [78.84–88.86] | 0 (0.0 %) [0–1.88] | 31 (15.5 %) [11.14–21.16] | 0 (0.0 %) [0–1.88] | 0 | 4.359 |

## T7 · distribution de d1 (nm) et distance au seuil 2,5 nm
| unité | min | 5 % | 25 % | médiane | 75 % | 95 % | max | d1<2,5 | d1<2,4 | dans [2,4;2,6) | plus proche du seuil (nm) | fermé sous bruit ±0,1 nm : observé % → [2,5 ; 50 ; 97,5] % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| eq162_n250 | 3.514 | 3.563 | 3.604 | 3.631 | 3.661 | 3.714 | 3.878 | 0 | 0 | 0 | 1.014 | 0.0 → [0.0, 0.0, 0.0] |
| eq164_n500 | 3.524 | 3.567 | 3.607 | 3.637 | 3.667 | 3.717 | 3.961 | 0 | 0 | 0 | 1.024 | 0.0 → [0.0, 0.0, 0.0] |
| prop162_lag0001 | 3.481 | 3.560 | 3.593 | 3.619 | 3.648 | 3.695 | 3.736 | 0 | 0 | 0 | 0.981 | 0.0 → [0.0, 0.0, 0.0] |
| prop162_lag0050 | 3.402 | 3.549 | 3.592 | 3.621 | 3.661 | 3.692 | 3.736 | 0 | 0 | 0 | 0.902 | 0.0 → [0.0, 0.0, 0.0] |
| prop162_lag1000 | 3.532 | 3.564 | 3.608 | 3.644 | 3.683 | 3.736 | 3.876 | 0 | 0 | 0 | 1.032 | 0.0 → [0.0, 0.0, 0.0] |
| auteur : GenAI (7 méthodes), son vrai p | 2.679 | 3.489 | 3.579 | 3.645 | 3.729 | 4.019 | 4.514 | 0 | 0 | 0 | 0.179 | — |
| auteur : BioEmu v1 | 3.167 | 3.461 | 3.550 | 3.608 | 3.683 | 3.846 | 4.409 | 0 | 0 | 0 | 0.667 | — |

## T8 · couverture (grilles de 0,25 nm ; p par pas de 0,25)
| unité | cases (d1,d2) occupées | entropie norm. d1 | entropie norm. p | d1 min–max | d2 min–max | p min / méd / max | \|p\| < 0,05 |
|---|---|---|---|---|---|---|---|
| eq162_n250 | 6 | 0.022 | 0.278 | 3.514–3.878 | 0.712–1.526 | -0.459 / 0.24 / 0.403 | 1 |
| eq164_n500 | 8 | 0.033 | 0.291 | 3.524–3.961 | 0.699–1.869 | -0.487 / 0.241 / 0.435 | 0 |
| prop162_lag0001 | 4 | 0.011 | 0.212 | 3.481–3.736 | 0.707–1.086 | 0.121 / 0.234 / 0.359 | 0 |
| prop162_lag0050 | 4 | 0.011 | 0.215 | 3.402–3.736 | 0.709–1.112 | 0.11 / 0.243 / 0.362 | 0 |
| prop162_lag1000 | 7 | 0.056 | 0.333 | 3.532–3.876 | 0.668–1.599 | -0.504 / 0.229 / 0.406 | 0 |
| auteur : GenAI (7), son vrai p | 71 | 0.315 | 0.379 | 2.679–4.514 | 0.604–4.336 | [-2.152, -0.293, 1.686] | — |
| auteur : BioEmu v1 (p non défini) | 26 | 0.257 | None | 3.167–4.409 | 0.616–3.485 | None | — |

## T9 · comparaison de d1 (KS à deux échantillons, Wasserstein en nm)
- d1_eq164_contre_BioEmu_v1 : KS = 0.275, p = 4.6e-17, Wasserstein = 0.058 nm
- d1_eq164_contre_GenAI_pool : KS = 0.225, p = 1.6e-18, Wasserstein = 0.084 nm

## T10 · ressemblance aux 221 structures fermées de l'auteur (RMSD Cα, 162 résidus, appariement ordinal)
- étalon (plus proche référence d'un AUTRE cluster) : médiane / 95e / max = [2.8, 11.62, 12.07] Å ; références → cristal 5LZM : min / méd / max = [5.96, 8.16, 13.72] Å
- eq162_n250 : n = 250, d1 < 2,5 nm : 0 ; RMSD min aux références (0/25/50/75/100 %) = [5.86, 6.04, 6.1, 6.16, 6.44] Å ; RMSD au cristal (médiane) = 0.77 Å
- eq164_n500 : n = 500, d1 < 2,5 nm : 0 ; RMSD min aux références (0/25/50/75/100 %) = [5.87, 6.05, 6.11, 6.16, 6.69] Å ; RMSD au cristal (médiane) = 0.76 Å
- prop162_lag0001 : n = 200, d1 < 2,5 nm : 0 ; RMSD min aux références (0/25/50/75/100 %) = [5.87, 6.01, 6.06, 6.12, 6.3] Å ; RMSD au cristal (médiane) = 0.66 Å
- prop162_lag0050 : n = 200, d1 < 2,5 nm : 0 ; RMSD min aux références (0/25/50/75/100 %) = [5.81, 6.06, 6.1, 6.15, 6.3] Å ; RMSD au cristal (médiane) = 0.74 Å
- prop162_lag1000 : n = 200, d1 < 2,5 nm : 0 ; RMSD min aux références (0/25/50/75/100 %) = [5.93, 6.05, 6.12, 6.19, 6.49] Å ; RMSD au cristal (médiane) = 0.78 Å
