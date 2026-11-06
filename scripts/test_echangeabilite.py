import json, statistics as st
from collections import defaultdict
from math import comb
B="/gpfs/phherman/bioemu2_bench/generation/md_contre_x0/"
md=defaultdict(lambda: defaultdict(dict))
for r in json.load(open(B+"mesures.json")):
    md[r["systeme"]][r["replique"]][r["ns"]]=r["dockq"]
BE={"1GCQ":[.853,.788,.765,.738,.669],"2X9A":[.880,.865,.850,.857,.841],"3PC8":[.869,.779,.759,.759,.713],
"1AY7":[.877,.864,.843,.821,.772],"1EFN":[.901,.837,.770,.675,.652],"1FFW":[.877,.830,.768,.776,.729],
"1QA9":[.847,.824,.794,.784,.715],"1SYX":[.860,.813,.782,.771,.721],"1J2J":[.821,.773,.728,.705,.631],
"1S1Q":[.854,.826,.821,.800,.772],"5E5M":[.868,.824,.818,.791,.783],"4POU":[.877,.806,.791,.745,.727],
"5SV3":[.896,.884,.861,.866,.817],"4Y7M":[.858,.857,.827,.802,.782],"5HGG":[.890,.844,.827,.793,.715],
"4M3K":[.910,.865,.845,.811,.769],"6CWG":[.890,.863,.823,.815,.764],"3K75":[.864,.835,.831,.821,.773],
"6DBG":[.890,.864,.855,.827,.834],"2OT3":[.867,.833,.769,.744,.713],"2C0L":[.847,.823,.813,.786,.726],
"1B6C":[.862,.839,.821,.807,.772],"3RJQ":[.823,.790,.730,.712,.647],"1M10":[.865,.818,.787,.777,.693],
"5JMO":[.916,.900,.883,.885,.870],"1BKD":[.869,.845,.799,.769,.726],"5VNW":[.828,.785,.734,.716,.659],
"2AJF":[.891,.867,.841,.815,.779]}

def binom_sup(k,n,p):
    return sum(comb(n,i)*p**i*(1-p)**(n-i) for i in range(k,n+1))

print("H0 : la valeur BioEmu2 est ECHANGEABLE avec les 3 repliques MD au meme temps.")
print("     Alors P(elle depasse les trois) = 1/4, P(sous les trois) = 1/4, attendu 7 sur 28.\n")
for li,lag in enumerate([1,10,50]):
    de=di=du=0; larg=[]
    for s in md:
        v=[md[s][r][float(lag)] for r in md[s]]
        if len(v)!=3 or s not in BE: continue
        larg.append(max(v)-min(v))
        b=BE[s][li]
        if b<min(v): de+=1
        elif b>max(v): du+=1
        else: di+=1
    n=de+di+du
    p=binom_sup(du,n,0.25)
    print(f"lag {lag:4d} ns : dessus {du}/{n} (attendu {n*0.25:.1f})  dedans {di}  dessous {de}  "
          f"| binomiale P(>= {du}) = {p:.2e}  | largeur mediane de la plage MD {st.median(larg):.3f}")

print("\n=== la plage des 3 repliques n est pas degeneree : elle s elargit avec le temps ===")
for t in [1,10,50,99.5]:
    larg=[max(md[s][r][float(t)] for r in md[s])-min(md[s][r][float(t)] for r in md[s]) for s in md]
    print(f"  t = {t:5} ns : largeur mediane {st.median(larg):.3f}, max {max(larg):.3f}")

print("\n=== et le signe est univoque : rang de BioEmu2 parmi les 4 valeurs ===")
for li,lag in enumerate([1,10,50]):
    rangs=[]
    for s in md:
        v=[md[s][r][float(lag)] for r in md[s]]
        if len(v)!=3 or s not in BE: continue
        rangs.append(sum(1 for x in v if x < BE[s][li]) + 1)   # 1 = le plus bas, 4 = le plus haut
    c={k:rangs.count(k) for k in (1,2,3,4)}
    print(f"  lag {lag:4d} ns : rang 1 {c[1]} · 2 {c[2]} · 3 {c[3]} · 4 {c[4]}   (attendu 7 chacun)")
