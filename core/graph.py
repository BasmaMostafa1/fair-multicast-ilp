"""
core.graph -- directed-graph container and random-instance generator.

Reconstructed to match the interface used by abilene.py and scalability.py:

    from core.graph import DirectedGraph, random_graph

    G = DirectedGraph(n)
    G.add_edge(u, v, delay, cost)        # directed arc u->v
    G.n                                  # node count
    G.arcs                               # list[Arc]
    G.out_arcs[u] / G.in_arcs[v]         # list[Arc] incident to a node
    dist, pred = G.shortest_delay_paths(s)   # Dijkstra on arc .delay

    G = random_graph(n=, edge_prob=, delay_range=(2,5), cost_range=(1,10), seed=)

Each Arc exposes .src, .dst, .delay, .cost, .idx.

NOTE: this is a faithful reconstruction of the model/interface, not a recovery
of the original module. The random_graph RNG stream is a plain per-ordered-pair
Bernoulli draw, so instance-for-instance numbers need not match another
generator that used a different draw order; distributional results are
generator-robust. Deterministic topologies (e.g. Abilene) reproduce exactly.
"""
import heapq
import random


class Arc:
    __slots__ = ("src", "dst", "delay", "cost", "idx")

    def __init__(self, src, dst, delay, cost, idx):
        self.src = src
        self.dst = dst
        self.delay = float(delay)
        self.cost = float(cost)
        self.idx = idx

    def __repr__(self):
        return f"Arc({self.src}->{self.dst}, d={self.delay:.3g}, c={self.cost:.3g})"


class DirectedGraph:
    def __init__(self, n):
        self.n = n
        self.arcs = []
        self.out_arcs = [[] for _ in range(n)]
        self.in_arcs = [[] for _ in range(n)]

    def add_arc(self, u, v, delay, cost=0.0):
        """Add a single directed arc u->v. Returns the Arc."""
        a = Arc(u, v, delay, cost, len(self.arcs))
        self.arcs.append(a)
        self.out_arcs[u].append(a)
        self.in_arcs[v].append(a)
        return a

    def add_edge(self, u, v, delay, cost=0.0):
        """
        Add an undirected link as two directed arcs u->v and v->u with the same
        delay and cost (a bidirectional link; e.g. symmetric fibre propagation).
        Returns (arc_uv, arc_vu). This is why 14 backbone links yield 28 arcs.
        """
        return self.add_arc(u, v, delay, cost), self.add_arc(v, u, delay, cost)

    def shortest_delay_paths(self, s):
        """
        Dijkstra on arc .delay from source s.
        Returns (dist, pred):
          dist[v] = minimum total delay s->v (math.inf if unreachable)
          pred[v] = the Arc used to reach v on a shortest-delay path (None for s)
        """
        INF = float("inf")
        dist = {v: INF for v in range(self.n)}
        pred = {v: None for v in range(self.n)}
        dist[s] = 0.0
        pq = [(0.0, s)]
        while pq:
            dd, u = heapq.heappop(pq)
            if dd > dist[u]:
                continue
            for a in self.out_arcs[u]:
                nd = dd + a.delay
                if nd < dist[a.dst] - 1e-15:
                    dist[a.dst] = nd
                    pred[a.dst] = a
                    heapq.heappush(pq, (nd, a.dst))
        return dist, pred


def random_graph(n, edge_prob=0.35, delay_range=(2.0, 5.0),
                 cost_range=(1.0, 10.0), seed=0):
    """
    Erdos-Renyi directed graph: each ordered pair (u,v), u!=v, gets an arc with
    probability edge_prob; arc delay ~ U(delay_range), cost ~ U(cost_range).
    """
    rng = random.Random(seed)
    G = DirectedGraph(n)
    dlo, dhi = delay_range
    clo, chi = cost_range
    for u in range(n):
        for v in range(n):
            if u != v and rng.random() < edge_prob:
                G.add_arc(u, v, rng.uniform(dlo, dhi), rng.uniform(clo, chi))
    return G
