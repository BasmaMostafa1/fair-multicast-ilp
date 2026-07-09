"""Abilene (Internet2) real-topology evaluation: 11 nodes, 14 bidirectional links.
Link delays = one-way fiber propagation from great-circle distance (5 us/km).
Costs proportional to distance. Proven-optimal hierarchy vs tree via HiGHS."""
import sys, math, json
sys.path.insert(0,"/home/claude/work/pkg")
from core.graph import DirectedGraph
from solve_highs import solve, jain

CITIES = {  # name: (lat, lon)
 0:("Seattle",47.61,-122.33), 1:("Sunnyvale",37.37,-122.04), 2:("LosAngeles",34.05,-118.24),
 3:("Denver",39.74,-104.99), 4:("KansasCity",39.10,-94.58), 5:("Houston",29.76,-95.37),
 6:("Atlanta",33.75,-84.39), 7:("Indianapolis",39.77,-86.16), 8:("Chicago",41.88,-87.63),
 9:("NewYork",40.71,-74.01), 10:("Washington",38.90,-77.04)}
LINKS = [(0,1),(0,3),(1,3),(1,2),(2,5),(3,4),(4,5),(4,7),(5,6),(6,7),(6,10),(7,8),(8,9),(9,10)]

def haversine(a,b):
    R=6371.0
    la1,lo1=math.radians(CITIES[a][1]),math.radians(CITIES[a][2])
    la2,lo2=math.radians(CITIES[b][1]),math.radians(CITIES[b][2])
    dla,dlo=la2-la1,lo2-lo1
    h=math.sin(dla/2)**2+math.cos(la1)*math.cos(la2)*math.sin(dlo/2)**2
    return 2*R*math.asin(math.sqrt(h))

G=DirectedGraph(11)
for u,v in LINKS:
    km=haversine(u,v)
    delay=km*0.005           # ms, ~5 us/km fiber propagation
    cost=km/100.0            # distance-proportional cost
    G.add_edge(u,v,delay,cost)
print(f"Abilene: {G.n} nodes, {len(G.arcs)} directed arcs ({len(LINKS)} links)")

s=0  # Seattle
M=[5,6,9,10,8,4]  # Houston, Atlanta, NewYork, Washington, Chicago, KansasCity
print("source:", CITIES[s][0], "| destinations:", [CITIES[m][0] for m in M])
spd,_=G.shortest_delay_paths(s)
maxspt=max(spd[m] for m in M); spread=maxspt-min(spd[m] for m in M)
print(f"SPT delays (ms): "+", ".join(f"{CITIES[m][0]}={spd[m]:.2f}" for m in M))
print(f"maxSPT={maxspt:.2f} ms, SPT spread={spread:.2f} ms\n")

# Delta sweep from just above maxSPT upward
deltas=[round(maxspt+x,2) for x in [0.01,2,4,6,8,12,16]]
rows=[]
print(f"{'Delta':>6} | {'dH':>6} {'JH':>7} | {'dT':>6} {'JT':>7} | {'(dT-dH)/dT':>9} {'sec':>5}")
for D in deltas:
    import time as _t; _t0=_t.time()
    rh=solve(G,s,M,D,tree=False,backend="highs",tl=120)
    rt=solve(G,s,M,D,tree=True, backend="highs",tl=120)
    el=_t.time()-_t0
    if rh.get('status')=='Optimal' and rt.get('status')=='Optimal':
        red=(rt['var']-rh['var'])/rt['var']*100 if rt['var']>1e-9 else 0.0
        rows.append(dict(D=round(D,2),dH=rh['var'],JH=rh['jain'],dT=rt['var'],JT=rt['jain'],red=red))
        print(f"{D:>6} | {rh['var']:6.3f} {rh['jain']:7.4f} | {rt['var']:6.3f} {rt['jain']:7.4f} | {red:8.0f}% {el:5.1f}")
    else:
        print(f"{D:>6} | H={rh.get('status')} T={rt.get('status')}")
json.dump({"meta":{"nodes":11,"arcs":len(G.arcs),"source":"Seattle","dests":[CITIES[m][0] for m in M],
                   "maxspt":round(maxspt,3),"spt_spread":round(spread,3)},"rows":rows},
          open("abilene_results.json","w"), indent=1)
