"""
solve_highs -- exact min-variation multicast ILP (hierarchy or tree) via HiGHS.

Reconstructed to match the interface used by abilene.py and scalability.py:

    from solve_highs import solve, jain

    r = solve(G, s, M, D, tree=False, backend="highs", tl=120)
    r['status']   # 'Optimal' iff proven optimal (both lexicographic stages)
    r['var']      # delta* : optimal inter-destination delay variation
    r['jain']     # Jain fairness index of the optimal delay vector
    r['delays']   # {m: d_m} optimal delay to each destination
    r['arcs']     # number of arcs carrying flow (|support of x|)

Model (Problem P3: min delta s.t. Delta), reconstructed from the manuscript:
    F(m,a) in {0,1}   per-destination unit flow      (E_oo2-E_oo5)
    x(a)   in [0,1]    arc usage, pinned to OR_m F(m,a):
                       x(a) >= F(m,a) for all m       (E_oo6)
                       x(a) <= sum_m F(m,a)           (forces x = the indicator,
                                                        so x is automatically 0/1)
                       x(a) = 0 for a into the source (E_oo0)
    d_m    cont.       = sum_a F(m,a) D(a),  d_m <= Delta   (E_co1)
    zvar   cont.       epigraph of OF4: zvar >= |d_m - d_n| (E_of4lin)
    objective          min zvar
    tree=True adds the single-predecessor constraint sum_{a in In(v)} x(a) <= 1
              for all v (E_o10); tree=False omits it (hierarchy feasible region).

The sigma arc-sharing block of the manuscript is omitted on purpose: it appears
in no objective or feasibility constraint, so it does not affect any optimum.

Reproducibility: a two-stage lexicographic objective (minimise variation, then
minimise total delay subject to variation = delta*) fixes a unique delay vector,
so the reported Jain index is reproducible. This matches the paper's stated
tie-break.

This is a faithful reconstruction of the model and the call interface, not a
recovery of the original module; proven-optimal values are solver-invariant.
"""
import pulp


# ----------------------------------------------------------------- helpers
def jain(values):
    """Jain's fairness index of a non-negative vector. 1.0 == perfectly equal."""
    vals = list(values)
    n = len(vals)
    if n == 0:
        return 1.0
    s = sum(vals)
    s2 = sum(v * v for v in vals)
    return (s * s) / (n * s2) if s2 > 0 else 1.0


def _make_solver(backend, tl):
    """Return a PuLP solver. backend='highs' -> HiGHS, else CBC; CBC fallback."""
    if backend == "highs":
        try:
            return pulp.HiGHS(msg=False, timeLimit=tl)
        except Exception:
            pass
    return pulp.PULP_CBC_CMD(msg=False, timeLimit=tl)


# ----------------------------------------------------------------- solve
def solve(G, s, M, D, tree=False, backend="highs", tl=120, tie_break=True):
    """
    Solve the min-variation multicast ILP on graph G, source s, destinations M,
    delay bound D. tree=False -> partial spanning hierarchy; tree=True -> tree.
    Returns a result dict (keys documented at module top).
    """
    M = list(M)
    arcs = G.arcs
    A = range(len(arcs))
    delay = {a: arcs[a].delay for a in A}

    def build():
        p = pulp.LpProblem("minvar", pulp.LpMinimize)
        F = {(m, a): pulp.LpVariable(f"F_{m}_{a}", cat="Binary") for m in M for a in A}
        x = {a: pulp.LpVariable(f"x_{a}", lowBound=0, upBound=1) for a in A}
        d = {m: pulp.LpVariable(f"d_{m}", lowBound=0) for m in M}
        z = pulp.LpVariable("zvar", lowBound=0)

        out = lambda v: [a.idx for a in G.out_arcs[v]]
        inc = lambda v: [a.idx for a in G.in_arcs[v]]

        for m in M:
            p += d[m] == pulp.lpSum(F[(m, a)] * delay[a] for a in A)
            p += d[m] <= D                                            # E_co1
            p += pulp.lpSum(F[(m, a)] for a in out(s)) == 1           # E_oo3
            p += pulp.lpSum(F[(m, a)] for a in inc(s)) == 0           # E_oo2
            p += (pulp.lpSum(F[(m, a)] for a in out(m))
                  == pulp.lpSum(F[(m, a)] for a in inc(m)) - 1)       # E_oo4
            for v in range(G.n):
                if v == s or v == m:
                    continue
                p += (pulp.lpSum(F[(m, a)] for a in out(v))
                      == pulp.lpSum(F[(m, a)] for a in inc(v)))       # E_oo5

        for i in range(len(M)):
            for j in range(i + 1, len(M)):
                p += z >= d[M[i]] - d[M[j]]                           # E_of4lin
                p += z >= d[M[j]] - d[M[i]]

        for a in A:
            p += x[a] <= pulp.lpSum(F[(m, a)] for m in M)             # pin x = OR_m F
            for m in M:
                p += x[a] >= F[(m, a)]                                # E_oo6
            if arcs[a].dst == s:
                p += x[a] == 0                                        # E_oo0
        for m in M:
            p += pulp.lpSum(x[a] for a in inc(m)) >= 1                # E_oo1
        if tree:
            for v in range(G.n):
                p += pulp.lpSum(x[a] for a in inc(v)) <= 1            # E_o10
        return p, F, x, d, z

    # ---- stage 1: minimise variation
    p, F, x, d, z = build()
    p += z
    p.solve(_make_solver(backend, tl))
    if pulp.LpStatus[p.status] != "Optimal":
        return {"status": pulp.LpStatus[p.status], "var": None, "jain": None,
                "delays": None, "arcs": None, "tree": tree}
    vstar = pulp.value(z)

    # ---- stage 2 (optional): fix variation, minimise total delay -> unique vector
    if tie_break:
        p2, F2, x2, d2, z2 = build()
        p2 += z2 <= vstar + 1e-6
        p2 += pulp.lpSum(d2[m] for m in M)
        p2.solve(_make_solver(backend, tl))
        if pulp.LpStatus[p2.status] == "Optimal":
            d, x = d2, x2
            vstar = pulp.value(z2)

    delays = {m: pulp.value(d[m]) for m in M}
    dl = list(delays.values())
    var = max(dl) - min(dl)
    used = sum(1 for a in A if (pulp.value(x[a]) or 0) > 0.5)
    return {"status": "Optimal", "var": var, "jain": jain(dl),
            "delays": delays, "arcs": used, "tree": tree}
