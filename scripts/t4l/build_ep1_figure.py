#!/usr/bin/env python3
"""build_ep1.py : carrousel de l'episode 1 (l'histoire de la cellule 48) et sa figure d1 contre p.

La figure est un SVG ecrit a la main depuis cv_dump.json : les d1 et p des sept methodes de l'auteur viennent de SES
fichiers (distances.npy, locked.pkl), ceux de BioEmu2 de cv.npz. Aucun recalcul ici : lecture, mise a l'echelle, trace."""
import json, pathlib, sys

D = json.load(open("cv_dump.json"))
assert D["n_rejets"] == 0, D["rejets"]

# ---------- la figure ----------
W, H = 1016, 760
ML, MR, MT, MB = 92, 18, 26, 86
X0, X1, Y0, Y1 = 1.8, 4.8, -2.5, 2.0          # nm ; sans unite
PW, PH = W - ML - MR, H - MT - MB
sx = lambda v: ML + (v - X0) / (X1 - X0) * PW
sy = lambda v: MT + (Y1 - v) / (Y1 - Y0) * PH
VEIL, BLUE, CARMIN, HAIR = "#c9cbe6", "#54c0f0", "#e30613", "rgba(84,192,240,0.28)"

bh, n_bh = [], 0
for nom, m in D["methodes"].items():
    for d1, p in zip(m["d1"], m["p"]):
        bh.append((d1, p)); n_bh += 1
b2 = list(zip(D["bioemu2"]["d1"], D["bioemu2"]["p"]))
n_b2 = len(b2)
assert (n_bh, n_b2) == (2197, 500), (n_bh, n_b2)
hors = [(d, p) for d, p in bh + b2 if not (X0 <= d <= X1 and Y0 <= p <= Y1)]
sous_seuil = sum(1 for d, _ in bh + b2 if d < 2.5)
assert sous_seuil == 0, sous_seuil

s = [f'<svg viewBox="0 0 {W} {H}" width="{W}" height="{H}" xmlns="http://www.w3.org/2000/svg" role="img" '
     f'aria-label="d1 against the locking coordinate p, seven published methods and BioEmu2">']
# grille legere + axes
for v in (-2, -1, 0, 1, 2):
    y = sy(v)
    dash = ' stroke-dasharray="10 8"' if v == 0 else ''
    s.append(f'<line x1="{ML}" x2="{W-MR}" y1="{y:.1f}" y2="{y:.1f}" stroke="{HAIR}" stroke-width="{3 if v==0 else 1.5}"{dash}/>')
    s.append(f'<text x="{ML-16}" y="{y+9:.1f}" text-anchor="end" font-size="26" fill="{VEIL}" font-family="ui-monospace,Menlo,monospace">{v}</text>')
for v in (2.0, 2.5, 3.0, 3.5, 4.0, 4.5):
    x = sx(v)
    s.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{MT+PH}" y2="{MT+PH+10}" stroke="{VEIL}" stroke-width="2"/>')
    s.append(f'<text x="{x:.1f}" y="{MT+PH+42}" text-anchor="middle" font-size="26" fill="{VEIL}" font-family="ui-monospace,Menlo,monospace">{v:g}</text>')
s.append(f'<line x1="{ML}" x2="{ML}" y1="{MT}" y2="{MT+PH}" stroke="{VEIL}" stroke-width="2"/>')
s.append(f'<line x1="{ML}" x2="{W-MR}" y1="{MT+PH}" y2="{MT+PH}" stroke="{VEIL}" stroke-width="2"/>')
# les points : les sept methodes publiees en retrait, BioEmu2 par-dessus
for d, p in bh:
    s.append(f'<circle cx="{sx(d):.1f}" cy="{sy(p):.1f}" r="3.6" fill="{VEIL}" fill-opacity="0.50"/>')
for d, p in b2:
    s.append(f'<circle cx="{sx(d):.1f}" cy="{sy(p):.1f}" r="5.2" fill="{BLUE}" fill-opacity="0.80" stroke="#242454" stroke-width="1"/>')
# le seuil ferme : la seule touche carmin de la planche
xs = sx(2.5)
s.append(f'<line x1="{xs:.1f}" x2="{xs:.1f}" y1="{MT}" y2="{MT+PH}" stroke="{CARMIN}" stroke-width="4" stroke-dasharray="14 10"/>')
s.append(f'<text x="{xs-14:.1f}" y="{MT+34}" text-anchor="end" font-size="28" fill="{CARMIN}" font-weight="600" '
         f'font-family="Helvetica Neue,Helvetica,Arial,sans-serif">closed state</text>')
s.append(f'<text x="{xs-14:.1f}" y="{MT+66}" text-anchor="end" font-size="28" fill="{CARMIN}" '
         f'font-family="Helvetica Neue,Helvetica,Arial,sans-serif">d1 &lt; 2.5 nm</text>')
# etats de Phe4
s.append(f'<text x="{W-MR-8}" y="{sy(1.75):.1f}" text-anchor="end" font-size="27" fill="{VEIL}" '
         f'font-family="Helvetica Neue,Helvetica,Arial,sans-serif">Phe4 exposed (p &gt; 0)</text>')
s.append(f'<text x="{W-MR-8}" y="{sy(-2.25):.1f}" text-anchor="end" font-size="27" fill="{VEIL}" '
         f'font-family="Helvetica Neue,Helvetica,Arial,sans-serif">Phe4 buried (p &lt; 0)</text>')
# titres d'axes
s.append(f'<text x="{ML+PW/2:.1f}" y="{H-14}" text-anchor="middle" font-size="28" fill="{VEIL}" '
         f'font-family="Helvetica Neue,Helvetica,Arial,sans-serif">d1 (nm), hinge opening</text>')
s.append(f'<text transform="translate(26,{MT+PH/2:.1f}) rotate(-90)" text-anchor="middle" font-size="28" fill="{VEIL}" '
         f'font-family="Helvetica Neue,Helvetica,Arial,sans-serif">locking coordinate p</text>')
s.append("</svg>")
SVG = "\n".join(s)
pathlib.Path("fig_d1_p.svg").write_text(SVG)
print(f"figure : {n_bh} + {n_b2} points, hors cadre {len(hors)}, sous le seuil {sous_seuil}", file=sys.stderr)

# ---------- le carrousel ----------
HEAD = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>BioEmu2 benchmark: Episode 1 carousel</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Raleway:ital,wght@0,400;0,500;0,600;0,700;1,400&display=swap">
<link rel="stylesheet" href="carrousel.css">
<style>
.beats { display:flex; flex-direction:column; margin-top:46px; }
.beat { display:grid; grid-template-columns: 250px 1fr; gap:28px; align-items:baseline; padding:30px 0; border-top:2px solid var(--hair); }
.beat:last-child { border-bottom:2px solid var(--hair); }
.beat .when { font-family:var(--data); font-size:26px; letter-spacing:.06em; color:var(--blue); text-transform:uppercase; }
.beat .what { font-size:34px; line-height:1.34; color:var(--white); }
.beat .what q { font-family:var(--display); font-weight:600; quotes:"\\201C" "\\201D"; }
.beat .what small { display:block; margin-top:10px; font-size:27px; line-height:1.38; color:var(--veil); }
.cmp { display:grid; grid-template-columns: 1fr; gap:0; margin-top:46px; }
.cmp .row { display:grid; grid-template-columns: 200px 1fr; gap:26px; align-items:baseline; padding:30px 0; border-top:2px solid var(--hair); }
.cmp .row:last-child { border-bottom:2px solid var(--hair); }
.cmp .k { font-family:var(--data); font-size:30px; color:var(--blue); }
.cmp .v { font-size:34px; line-height:1.34; }
.cmp .v b { font-family:var(--display); font-weight:600; }
.legend { display:flex; gap:38px; margin:34px 0 8px; position:relative; z-index:2; }
.legend span { font-size:26px; color:var(--veil); display:flex; align-items:center; gap:12px; }
.legend i { width:20px; height:20px; border-radius:50%; display:inline-block; }
.fig { position:relative; z-index:2; margin-top:6px; }
.tag { max-width:none; }
.cred { font-size:22px; line-height:1.4; color:var(--veil); max-width:none; }
</style>
</head>
<body>
"""

B1 = """<section class="board">
  <img class="mark" src="logo-immunochem.png" alt="">
  <p class="tag">BioEmu2 &middot; episode 1</p>
  <div class="spacer"></div>
  <h2 class="sm z">We reproduced a published benchmark to 5&times;10<sup style="font-size:.55em">&minus;7</sup>&nbsp;nm.</h2>
  <div class="rule red z"></div>
  <p class="w on z" style="font-size:44px; line-height:1.3">One cell would not match.</p>
  <p class="w z" style="margin-top:34px">So we wrote to the author, and said we were probably the ones who were wrong.</p>
  <div class="spacer"></div>
  <div class="foot"><span>Immunochem</span><span class="sw">swipe &rarr;</span></div>
</section>
"""

B2 = """<section class="board">
  <img class="mark soft" src="logo-immunochem.png" alt="">
  <p class="tag">The cell</p>
  <div class="spacer" style="flex:0 0 40px"></div>
  <h2 class="sm z">In cell 48 of the notebook, the locking coordinate was another distance.</h2>
  <div class="cmp z">
    <div class="row"><span class="k">column 1</span><span class="v"><b>d2</b>, a distance. Never negative, so every structure lands on the same side: one state, 1,832 of 1,832.</span></div>
    <div class="row"><span class="k">his own p</span><span class="v">The same 1,832 structures split <span class="red"><b>43.3 % exposed, 56.7 % buried</b></span>.</span></div>
  </div>
  <p class="w z" style="margin-top:44px">Our recomputation matched his arrays everywhere else: 8 datasets, and all 221 of his closed references classified correctly.</p>
  <div class="spacer"></div>
  <div class="foot"><span>Bhakat, J. Chem. Inf. Model. 2026</span><span class="sw">swipe &rarr;</span></div>
</section>
"""

B3 = """<section class="board">
  <img class="mark soft" src="logo-immunochem.png" alt="">
  <p class="tag">Three days</p>
  <div class="spacer" style="flex:0 0 40px"></div>
  <h2 class="sm z">We wrote first, and put our own error as the likeliest explanation.</h2>
  <div class="beats z">
    <div class="beat"><span class="when">Friday<br>morning</span><span class="what">We write to him. <q>The most likely explanation is an error of ours.</q></span></div>
    <div class="beat"><span class="when">Seven hours<br>later</span><span class="what"><q><span class="red">This is a great catch.</span></q><small>He re-plotted the figure, and will ask the journal about an amendment.</small></span></div>
    <div class="beat"><span class="when">Monday</span><span class="what"><q>Please keep everything public.</q> <q>You are welcome to tag me.</q></span></div>
  </div>
  <div class="spacer"></div>
  <div class="foot"><span>Soumendranath Bhakat, author communication</span><span class="sw">swipe &rarr;</span></div>
</section>
"""

B4 = f"""<section class="board">
  <img class="mark soft" src="logo-immunochem.png" alt="">
  <p class="tag">Corrected plot + BioEmu2</p>
  <div class="spacer" style="flex:0 0 26px"></div>
  <h2 class="sm z" style="font-size:52px">Both states of p are populated. Nothing gets near the closed state.</h2>
  <div class="legend">
    <span><i style="background:#c9cbe6;opacity:.6"></i>his 7 methods &middot; 2,197 structures</span>
    <span><i style="background:#54c0f0"></i>BioEmu2 &middot; 500 structures</span>
  </div>
  <div class="fig">{SVG}</div>
  <div class="spacer"></div>
  <div class="foot"><span>his files, our recomputation &middot; BioEmu2: ours</span><span class="sw">swipe &rarr;</span></div>
</section>
"""

B5 = """<section class="board">
  <img class="mark soft" src="logo-immunochem.png" alt="">
  <p class="tag">What changes</p>
  <div class="spacer"></div>
  <h2 class="sm z">Corrected: p.<br>Unchanged: no generative ensemble closes the hinge.</h2>
  <div class="rule red z"></div>
  <p class="w on z">His seven methods: 0 of 2,197 reach the closed state. BioEmu2: 0 of 500.</p>
  <p class="w z" style="margin-top:30px">BioEmu2 sits on the crystal side of p (median 0.241, the crystal&rsquo;s own value) and covers 8 bins of conformational space where BioEmu v1 covered 26.</p>
  <div class="spacer"></div>
  <p class="cred z" style="margin-bottom:170px">Correction to the paper: author communication. He will raise an amendment with the journal. Thank you, Soumendranath Bhakat.</p>
  <div class="sig"><img src="logo-immunochem.png" alt="Immunochem"><span>IMMUNOCHEM</span></div>
  <div class="foot" style="justify-content:flex-end"><span>Full report, data and scripts: link in the post</span></div>
</section>
"""
pathlib.Path("ep1-carrousel.html").write_text(HEAD + B1 + B2 + B3 + B4 + B5 + "</body></html>\n")
print("ep1-carrousel.html ecrit", file=sys.stderr)
