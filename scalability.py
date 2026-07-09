"""ILP scalability: solve the proven hierarchy min-variation ILP over growing n.
Random instances (paper model), |M|=5, Delta=1.5*maxSPT, HiGHS, per-solve cap."""
import sys, math, random, time, json
sys.path.insert(0,"/home/claude/work/pkg")
import pulp
from core.graph import random_graph
from solve_highs import solve

def maxspt(G,s,M):
    import heapq
    d={v:math.inf for v in range(G.n)}; d[s]=0.0; h=[(0.0,s)]
    while h:
        dd,u=heapq.heappop(h)
        if dd>d[u]:continue
        for a in G.out_arcs[u]:
            nd=dd+a.delay
            if nd<d[a.dst]: d[a.dst]=nd; heapq.heappush(h,(nd,a.dst))
    return max(d[m] for m in M) if all(math.isfinite(d[m]) for m in M) else None

def run(n_list, ndest=5, ep=0.35, frac=1.5, ninst=2, tl=90, seed_base=100):
    out=[]
    for n in n_list:
        times=[]; proven=0; arcs=0; binaries=0; dstar=[]; got=0; seed=seed_base
        while got<ninst and seed<seed_base+ninst*6:
            seed+=1
            G=random_graph(n=n,edge_prob=ep,delay_range=(2,5),cost_range=(1,10),seed=seed)
            rng=random.Random(seed); s=rng.choice(list(range(n)))
            M=rng.sample([v for v in range(n) if v!=s], min(ndest,n-1))
            ms=maxspt(G,s,M)
            if ms is None: continue
            D=frac*ms
            nb=len(G.arcs)*(len(M)+1)   # F (|M|*|arcs|) + y (|arcs|) binaries
            t=time.time(); r=solve(G,s,M,D,tree=False,backend="highs",tl=tl); el=time.time()-t
            got+=1; times.append(el); arcs=len(G.arcs); binaries=nb
            if r.get('status')=='Optimal':
                proven+=1; dstar.append(r['var'])
        row=dict(n=n, arcs=arcs, binaries=binaries, ninst=got, proven=proven,
                 mean_sec=round(sum(times)/len(times),1) if times else None,
                 max_sec=round(max(times),1) if times else None,
                 mean_dstar=round(sum(dstar)/len(dstar),3) if dstar else None)
        out.append(row)
        print(f"n={n:>2} arcs={arcs:>3} bin={binaries:>4} | proven {proven}/{got} | "
              f"mean {row['mean_sec']:>6}s max {row['max_sec']:>6}s | d*~{row['mean_dstar']}")
    return out

if __name__=="__main__":
    import sys
    nl=[int(x) for x in sys.argv[1].split(",")]
    res=run(nl, tl=int(sys.argv[2]) if len(sys.argv)>2 else 90)
    json.dump(res, open(f"scal_{nl[0]}_{nl[-1]}.json","w"), indent=1)
