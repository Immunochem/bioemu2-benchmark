# Confusions possibles du contraste A_apparie / B (lag 1000 ns)

A_apparie 24 systemes, B 26 systemes. Taux de succes par image (classe CAPRI >= acceptable), pooles sur les bras SANS orientation (R +10/+20/+30, RD +20 : 12 images par systeme) ou AVEC (T +10/+20/+30 : 9 images).

## 1. Correlations de Spearman systeme par systeme (taux des bras sans orientation)

| covariable | A_apparie : rho (p) | B : rho (p) |
|---|---|---|
| tokens | -0.29 (0.168) | +0.06 (0.776) |
| contacts | -0.49 (0.016) | +0.46 (0.018) |
| d_nat | +0.24 (0.266) | -0.22 (0.290) |
| msa | +0.04 (0.845) | +0.19 (0.348) |
| derive_c0 | -0.33 (0.115) | -0.13 (0.534) |

## 2. A_apparie contre B dans des tranches de covariable (bras SANS orientation, puis AVEC)

Difference = taux A_apparie − taux B, [IC95 bootstrap sur les systemes] ; entre parentheses : systemes A / B.

| tranche | sans orientation : A | B | difference | avec orientation : A | B | difference |
|---|---|---|---|---|---|---|
| tous | 0.34 (24) | 0.16 (26) | +0.18 [+0.09;+0.28] | 0.81 (24) | 0.59 (26) | +0.22 [+0.02;+0.40] |
| moins de 250 residus | 0.35 (7) | 0.12 (7) | +0.22 [+0.01;+0.41] | 0.78 (7) | 0.46 (7) | +0.32 [-0.06;+0.62] |
| 250 a 349 residus | 0.40 (9) | 0.19 (8) | +0.21 [+0.03;+0.38] | 0.83 (9) | 0.58 (8) | +0.25 [-0.12;+0.56] |
| 350 residus et plus | 0.26 (8) | 0.16 (11) | +0.10 [-0.01;+0.20] | 0.82 (8) | 0.70 (11) | +0.12 [-0.17;+0.42] |
| contacts natifs <= 70 | 0.45 (11) | 0.10 (15) | +0.35 [+0.22;+0.47] | 0.79 (11) | 0.44 (15) | +0.35 [+0.10;+0.57] |
| contacts natifs > 70 | 0.24 (13) | 0.26 (11) | -0.01 [-0.16;+0.11] | 0.83 (13) | 0.84 (11) | -0.01 [-0.25;+0.22] |
| d_COM natif >= 20 A | 0.37 (17) | 0.14 (20) | +0.23 [+0.12;+0.34] | 0.78 (17) | 0.59 (20) | +0.20 [-0.04;+0.42] |
| d_COM natif < 20 A | 0.27 (7) | 0.21 (6) | +0.06 [-0.24;+0.23] | 0.87 (7) | 0.60 (6) | +0.28 [-0.15;+0.59] |
| MSA >= 100 et derive C0 <= 3 A | 0.35 (20) | 0.20 (17) | +0.16 [+0.03;+0.27] | 0.83 (20) | 0.70 (17) | +0.14 [-0.08;+0.34] |
| MSA >= 1000 et derive C0 <= 3 A | 0.36 (16) | 0.17 (11) | +0.19 [+0.05;+0.31] | 0.81 (16) | 0.69 (11) | +0.12 [-0.15;+0.38] |
| MSA >= 100, derive <= 3 A, 150-400 residus, contacts <= 110 | 0.39 (15) | 0.18 (11) | +0.21 [+0.07;+0.33] | 0.83 (15) | 0.62 (11) | +0.21 [-0.05;+0.48] |

## 3. Les 12 premiers A_apparie (chaines A des paires negatives) contre les 12 suivants

| sous-ensemble | systemes | residus (mediane) | contacts natifs (mediane) | d_COM natif (mediane) | taux sans orientation | taux avec orientation |
|---|---|---|---|---|---|---|
| A_apparie : chaines A des paires negatives | 12 | 300 | 64 | 25.8 | 0.42 [0.32;0.52] | 0.83 [0.70;0.94] |
| A_apparie : autres | 12 | 338 | 82 | 25.8 | 0.26 [0.16;0.36] | 0.79 [0.56;0.95] |
| B | 26 | 323 | 64 | 24.8 | 0.16 [0.10;0.22] | 0.59 [0.44;0.75] |
