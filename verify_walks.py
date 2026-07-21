"""Post-hoc verification that every proven-optimal per-destination flow F(m, .)
corresponds to a genuine multicast hierarchy under Definition 3 of the paper.

Motivation
----------
The flow-conservation family (Eqs. 18-21 in the manuscript) admits, in
principle, integer flows that decompose into an s-to-m walk PLUS one or more
arc-disjoint directed circulations detached from the walk. Detached
circulations increase d_m = sum_a F(m,a) * D(a) and thus give the
delta-minimizing solver a spurious padding mechanism that corresponds to no
physical hierarchy. Circulations that share a node with the walk are fine:
they are legitimate delay-padding detours of a walk-based hierarchy.

This module verifies, a posteriori, that no reported optimum contains a
detached circulation. It is called after each solve in `abilene.py` and
`scalability.py` and is expected to pass on every instance reported in the
paper.

Usage
-----
    from verify_walks import verify_from_solve_result
    res = solve(G, s, M, D, tree=False, backend="highs", tl=120)
    ok, report = verify_from_solve_result(G, s, M, res)
    assert ok, report

`solve_highs.solve()` must expose the per-destination flow matrix under the
key 'F'. Add one line at the end of `solve()` before it returns:

    res['F'] = {m: {(a.src, a.dst): pulp.value(F[m, a]) for a in G.arcs}
                for m in M}
"""

from collections import defaultdict


def verify(G, s, M, F, tol=1e-6):
    """Check that every F(m, .) is a single s-to-m walk plus (optionally)
    circulations sharing at least one node with the walk.

    Parameters
    ----------
    G : object with attributes `.n` (int, number of nodes) and `.arcs`
        (iterable of arc objects or (u, v) pairs).  Only the arc endpoints
        matter for this check.
    s : source node id.
    M : iterable of destination node ids.
    F : dict {m: {(u, v): value}} mapping each destination to its solved
        arc-flow values.
    tol : float, threshold above which an arc is considered active.

    Returns
    -------
    all_ok : bool -- True iff every destination passes.
    report : dict {m: {...}} with per-destination diagnostics.
    """
    report = {}
    all_ok = True
    for m in M:
        active = [(u, v) for (u, v), val in F[m].items() if val > 1 - tol]

        # BFS on the active subgraph from s
        out_adj = defaultdict(list)
        for u, v in active:
            out_adj[u].append(v)
        seen = {s}
        stack = [s]
        while stack:
            u = stack.pop()
            for v in out_adj[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)

        dest_reached = m in seen
        # Detached: both endpoints unreachable from s in the active subgraph.
        detached = [(u, v) for (u, v) in active
                    if u not in seen and v not in seen]
        ok = dest_reached and not detached

        report[m] = {
            "arcs_active": len(active),
            "nodes_reached_from_s": len(seen),
            "dest_reached": dest_reached,
            "detached_arcs": detached,
            "ok": ok,
        }
        all_ok = all_ok and ok
    return all_ok, report


def verify_from_solve_result(G, s, M, solve_result):
    """Wrapper: extract F from a solve() result dict and call verify()."""
    if "F" not in solve_result:
        raise KeyError(
            "solve() must expose the per-destination flow matrix under key "
            "'F'. Add one line at the end of solve_highs.solve():\n"
            "    res['F'] = {m: {(a.src, a.dst): pulp.value(F[m, a]) "
            "for a in G.arcs} for m in M}"
        )
    return verify(G, s, M, solve_result["F"])


# --------------------------------------------------------------------------
# Self-test: run `python verify_walks.py` to confirm the checker works.
# --------------------------------------------------------------------------
if __name__ == "__main__":
    class _G:
        n = 5
        arcs = [(0, 1), (1, 2), (2, 1), (3, 4), (4, 3)]

    # Case 1: clean walk 0 -> 1 -> 2 (destination = 2).
    F_clean = {2: {(0, 1): 1, (1, 2): 1, (2, 1): 0, (3, 4): 0, (4, 3): 0}}
    ok, r = verify(_G, 0, [2], F_clean)
    assert ok, f"clean case failed: {r}"

    # Case 2: walk plus a DETACHED cycle 3 -> 4 -> 3 unreachable from s.
    # Flow-conservation-balanced at nodes 3 and 4; must be rejected.
    F_detached = {2: {(0, 1): 1, (1, 2): 1, (2, 1): 0, (3, 4): 1, (4, 3): 1}}
    ok, r = verify(_G, 0, [2], F_detached)
    assert not ok, f"detached case should have been rejected: {r}"
    assert set(r[2]["detached_arcs"]) == {(3, 4), (4, 3)}

    # Case 3: walk plus a circulation 1 -> 2 -> 1 sharing nodes with the walk.
    # This is a legitimate padding detour of a genuine hierarchy and must PASS.
    F_padded = {2: {(0, 1): 1, (1, 2): 1, (2, 1): 1, (3, 4): 0, (4, 3): 0}}
    ok, r = verify(_G, 0, [2], F_padded)
    assert ok, f"padded (walk-touching cycle) case failed: {r}"

    print("verify_walks self-test passed.")
