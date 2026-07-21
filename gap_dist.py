"""
Distribution of the exact tree-vs-hierarchy delay-variation gap
(delta*_T - delta*_H)/delta*_T over random Erdos-Renyi instances.

The ILP is reconstructed faithfully from the manuscript:
  variables  F(m,a) in {0,1}      (per-destination unit flow, E_oo2-E_oo5)
             x(a)  in [0,1]        (arc usage; pinned to the true indicator by
                                    E_oo6  x(a) >= F(m,a)  AND  x(a) <= sum_m F(m,a),
                                    so x is automatically {0,1})
             d_m   continuous      (= sum_a F(m,a) D(a))
             zvar  continuous      (epigraph of OF4, E_of4lin)
  objective  min zvar              (Problem P3: min delta s.t. Delta)
  tree model adds the single-predecessor constraint E_o10.
The sigma sharing block is omitted: it is inert in the published model
(appears in no objective/feasibility constraint) so it does not affect any optimum.
"""
import sys, math, random, json, time
import pulp

# ---------------------------------------------------------------- instance gen
def random_graph(n, p, seed, dlo=2.0, dhi=5.0, clo=1.0, chi=10.0):
    rng = random.Random(seed)
    arcs = []                      # list of (u, v, delay, cost)
    for u in range(n):
        for v in range(n):
            if u != v and rng.random() < p:
                arcs.append((u, v, rng.uniform(dlo, dhi), rng.uniform(clo, chi)))
    return n, arcs

def out_in(n, arcs):
    OUT = {v: [] for v in range(n)}
    IN  = {v: [] for v in range(n)}
    for i, (u, v, d, c) in enumerate(arcs):
        OUT[u].append(i); IN[v].append(i)
    return OUT, IN

def max_spt(n, arcs, s, M):
    import heapq
    OUT, _ = out_in(n, arcs)
    dist = {v: math.inf for v in range(n)}; dist[s] = 0.0
    h = [(0.0, s)]
    while h:
        dd, u = heapq.heappop(h)
        if dd > dist[u]: continue
        for ai in OUT[u]:
            _, v, d, _ = arcs[ai]
            nd = dd + d
            if nd < dist[v]:
                dist[v] = nd; heapq.heappush(h, (nd, v))
    if any(math.isinf(dist[m]) for m in M):
        return None
    return max(dist[m] for m in M)

# ---------------------------------------------------------------- ILP solve
def solve_minvar(n, arcs, s, M, Delta, tree, tl=25):
    OUT, IN = out_in(n, arcs)
    A = range(len(arcs))
    prob = pulp.LpProblem("minvar", pulp.LpMinimize)
    F = {(m, a): pulp.LpVariable(f"F_{m}_{a}", cat="Binary") for m in M for a in A}
    x = {a: pulp.LpVariable(f"x_{a}", lowBound=0, upBound=1) for a in A}
    d = {m: pulp.LpVariable(f"d_{m}", lowBound=0) for m in M}
    zvar = pulp.LpVariable("zvar", lowBound=0)
    prob += zvar                                              # OF4 epigraph objective

    for m in M:
        # delay definition
        prob += d[m] == pulp.lpSum(F[(m, a)] * arcs[a][2] for a in A)
        prob += d[m] <= Delta                                 # E_co1
        # flow conservation (single unit s->m)
        prob += pulp.lpSum(F[(m, a)] for a in OUT[s]) == 1     # E_oo3
        prob += pulp.lpSum(F[(m, a)] for a in IN[s])  == 0     # E_oo2
        prob += (pulp.lpSum(F[(m, a)] for a in OUT[m])
                 == pulp.lpSum(F[(m, a)] for a in IN[m]) - 1)  # E_oo4
        for v in range(n):
            if v == s or v == m: continue
            prob += (pulp.lpSum(F[(m, a)] for a in OUT[v])
                     == pulp.lpSum(F[(m, a)] for a in IN[v]))  # E_oo5
    # epigraph rows for variation
    Ml = list(M)
    for i in range(len(Ml)):
        for j in range(i + 1, len(Ml)):
            prob += zvar >= d[Ml[i]] - d[Ml[j]]               # E_of4lin
            prob += zvar >= d[Ml[j]] - d[Ml[i]]
    # arc-usage pinning  x(a) = OR_m F(m,a)
    for a in A:
        prob += x[a] <= pulp.lpSum(F[(m, a)] for m in M)
        for m in M:
            prob += x[a] >= F[(m, a)]                          # E_oo6
        if arcs[a][1] == s:
            prob += x[a] == 0                                 # E_oo0
    for m in M:
        prob += pulp.lpSum(x[a] for a in IN[m]) >= 1           # E_oo1
    if tree:
        for v in range(n):
            prob += pulp.lpSum(x[a] for a in IN[v]) <= 1       # E_o10
    try:
        solver = pulp.HiGHS(msg=False, timeLimit=tl)
    except Exception:
        solver = pulp.PULP_CBC_CMD(msg=False, timeLimit=tl)
    prob.solve(solver)
    st = pulp.LpStatus[prob.status]
    val = pulp.value(zvar)
    return st, (val if val is not None else None)

# ---------------------------------------------------------------- driver
def run(sizes, ndest=5, p=0.35, frac=1.5, nseed=40, tl=25, seed0=1000):
    out = {}
    for n in sizes:
        rows = []
        seed = seed0
        attempts = 0
        while len(rows) < nseed and attempts < nseed * 5:
            attempts += 1; seed += 1
            n_, arcs = random_graph(n, p, seed)
            rng = random.Random(seed * 7 + 1)
            s = rng.randrange(n)
            cand = [v for v in range(n) if v != s]
            M = rng.sample(cand, min(ndest, len(cand)))
            ms = max_spt(n, arcs, s, M)
            if ms is None:
                continue
            Delta = frac * ms
            stH, vH = solve_minvar(n, arcs, s, M, Delta, tree=False, tl=tl)
            stT, vT = solve_minvar(n, arcs, s, M, Delta, tree=True,  tl=tl)
            if stH != "Optimal" or stT != "Optimal" or vT is None or vH is None:
                continue
            if vT < 1e-6:
                continue                       # degenerate (tree already equal)
            ratio = (vT - vH) / vT
            rows.append(dict(seed=seed, dH=round(vH, 4), dT=round(vT, 4),
                             ratio=round(100 * ratio, 2)))
            print(f"n={n} seed={seed} dH={vH:.3f} dT={vT:.3f} "
                  f"ratio={100*ratio:5.1f}%  (kept {len(rows)}/{nseed})", flush=True)
        out[n] = rows
    return out

if __name__ == "__main__":
    sizes = [int(z) for z in sys.argv[1].split(",")] if len(sys.argv) > 1 else [12]
    nseed = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    res = run(sizes, nseed=nseed)
    json.dump(res, open("gap_dist_results.json", "w"), indent=1)
    print("WROTE gap_dist_results.json")
