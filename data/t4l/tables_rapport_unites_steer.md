
## T6 · populations des quatre états, ses CV et ses seuils (`resultats/analyse_unites.json`)
| unité | système | lag | n | physiques | rejetées | exp/open | exp/closed | bur/open | bur/closed | indéfinis | s/struct |
|---|---|---|---|---|---|---|---|---|---|---|---|
| st_d1_f005 | s164 | — | 100 | 100 | 0 | 91 (91.0 %) [83.77–95.19] | 0 (0.0 %) [0–3.7] | 9 (9.0 %) [4.81–16.23] | 0 (0.0 %) [0–3.7] | 0 | 6.852 |
| st_d1_f020 | s164 | — | 100 | 53 | 47 | 0 (0.0 %) [0–3.7] | 94 (94.0 %) [87.52–97.22] | 0 (0.0 %) [0–3.7] | 6 (6.0 %) [2.78–12.48] | 0 | 6.848 |
| st_d1d2_f020 | s164 | — | 100 | 0 | 100 | 0 (0.0 %) [0–3.7] | 95 (95.0 %) [88.82–97.85] | 0 (0.0 %) [0–3.7] | 5 (5.0 %) [2.15–11.18] | 0 | 6.862 |
| st_f000_inerte | s164 | — | 100 | 100 | 0 | 91 (91.0 %) [83.77–95.19] | 0 (0.0 %) [0–3.7] | 9 (9.0 %) [4.81–16.23] | 0 (0.0 %) [0–3.7] | 0 | 6.86 |

## T7 · distribution de d1 (nm) et distance au seuil 2,5 nm
| unité | min | 5 % | 25 % | médiane | 75 % | 95 % | max | d1<2,5 | d1<2,4 | dans [2,4;2,6) | plus proche du seuil (nm) | fermé sous bruit ±0,1 nm : observé % → [2,5 ; 50 ; 97,5] % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| st_d1_f005 | 2.795 | 2.810 | 2.831 | 2.847 | 2.860 | 2.880 | 2.958 | 0 | 0 | 0 | 0.295 | 0.0 → [0.0, 0.0, 1.0] |
| st_d1_f020 | 2.333 | 2.340 | 2.355 | 2.362 | 2.365 | 2.368 | 2.372 | 100 | 100 | 0 | 0.128 | 100.0 → [86.0, 92.0, 97.0] |
| st_d1d2_f020 | 2.336 | 2.342 | 2.354 | 2.361 | 2.364 | 2.368 | 2.371 | 100 | 100 | 0 | 0.129 | 100.0 → [86.0, 92.0, 97.0] |
| st_f000_inerte | 3.538 | 3.553 | 3.603 | 3.632 | 3.661 | 3.708 | 3.856 | 0 | 0 | 0 | 1.038 | 0.0 → [0.0, 0.0, 0.0] |
| auteur : GenAI (7 méthodes), son vrai p | 2.679 | 3.489 | 3.579 | 3.645 | 3.729 | 4.019 | 4.514 | 0 | 0 | 0 | 0.179 | — |
| auteur : BioEmu v1 | 3.167 | 3.461 | 3.550 | 3.608 | 3.683 | 3.846 | 4.409 | 0 | 0 | 0 | 0.667 | — |

## T8 · couverture (grilles de 0,25 nm ; p par pas de 0,25)
| unité | cases (d1,d2) occupées | entropie norm. d1 | entropie norm. p | d1 min–max | d2 min–max | p min / méd / max | \|p\| < 0,05 |
|---|---|---|---|---|---|---|---|
| st_d1_f005 | 4 | -0.0 | 0.294 | 2.795–2.958 | 0.642–1.621 | -0.463 / 0.239 / 0.395 | 0 |
| st_d1_f020 | 3 | -0.0 | 0.276 | 2.333–2.372 | 0.633–1.354 | -0.414 / 0.24 / 0.393 | 0 |
| st_d1d2_f020 | 1 | -0.0 | 0.267 | 2.336–2.371 | 2.838–2.849 | -0.446 / 0.235 / 0.401 | 1 |
| st_f000_inerte | 4 | 0.019 | 0.293 | 3.538–3.856 | 0.659–1.651 | -0.465 / 0.236 / 0.377 | 0 |
| auteur : GenAI (7), son vrai p | 71 | 0.315 | 0.379 | 2.679–4.514 | 0.604–4.336 | [-2.152, -0.293, 1.686] | — |
| auteur : BioEmu v1 (p non défini) | 26 | 0.257 | None | 3.167–4.409 | 0.616–3.485 | None | — |

## T9 · comparaison de d1 (KS à deux échantillons, Wasserstein en nm)

## T10 · ressemblance aux 221 structures fermées de l'auteur (RMSD Cα, 162 résidus, appariement ordinal)
- étalon (plus proche référence d'un AUTRE cluster) : médiane / 95e / max = [2.8, 11.62, 12.07] Å ; références → cristal 5LZM : min / méd / max = [5.96, 8.16, 13.72] Å
- st_d1d2_f020 : n = 100, d1 < 2,5 nm : 100 ; RMSD min aux références (0/25/50/75/100 %) = [5.8, 5.94, 5.99, 6.06, 6.24] Å ; RMSD au cristal (médiane) = 1.52 Å ; FERMÉES : {'rmsd_min_aux_refs_A_min_med_max': [5.8, 5.99, 6.24], 'n_sous_etalon_95_(NON DISCRIMINANT)': 100, 'etalon_95_A': 11.62, 'n_sous_etalon_mediane': 0, 'etalon_mediane_A': 2.8, 'rmsd_au_cristal_A_min_med_max': [1.43, 1.52, 2.49], 'd2_nm_min_med_max': [2.838, 2.845, 2.849], 'p_median': 0.235, 'exposed_closed': 95, 'buried_closed': 5, 'physiques': 0}
- st_d1_f005 : n = 100, d1 < 2,5 nm : 0 ; RMSD min aux références (0/25/50/75/100 %) = [5.76, 5.97, 6.02, 6.1, 6.3] Å ; RMSD au cristal (médiane) = 0.87 Å
- st_d1_f020 : n = 100, d1 < 2,5 nm : 100 ; RMSD min aux références (0/25/50/75/100 %) = [5.72, 5.96, 6.01, 6.09, 6.4] Å ; RMSD au cristal (médiane) = 1.02 Å ; FERMÉES : {'rmsd_min_aux_refs_A_min_med_max': [5.72, 6.01, 6.4], 'n_sous_etalon_95_(NON DISCRIMINANT)': 100, 'etalon_95_A': 11.62, 'n_sous_etalon_mediane': 0, 'etalon_mediane_A': 2.8, 'rmsd_au_cristal_A_min_med_max': [0.88, 1.02, 1.97], 'd2_nm_min_med_max': [0.633, 0.842, 1.354], 'p_median': 0.24, 'exposed_closed': 94, 'buried_closed': 6, 'physiques': 53}
- st_f000_inerte : n = 100, d1 < 2,5 nm : 0 ; RMSD min aux références (0/25/50/75/100 %) = [5.84, 6.03, 6.08, 6.16, 6.4] Å ; RMSD au cristal (médiane) = 0.76 Å
