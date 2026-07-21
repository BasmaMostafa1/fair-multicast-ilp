# Fair Multicast Routing — Exact ILP Formulations

Data and solver scripts for the manuscript:

> **ILP Formulation and Exact Solution of Multicast Routing with Fairness**
> Basma Mostafa, Hanan Haj Ahmad, Miklós Molnár
> Submitted to *Mathematics* (MDPI), Special Issue "Advanced Optimization Methods and Applications, 3rd Edition".

Exact flow-based integer linear programming (ILP) formulations for fairness-aware multicast routing, optimizing the maximum end-to-end delay (Δ) and the inter-destination delay variation (δ) as objectives or constraints. Feasible solutions are partial spanning hierarchies, a class that contains partial spanning trees as a special case.

## Contents

**Drivers (reproduce the three empirical sections of the paper):**

| File | Description |
|---|---|
| `abilene.py` | Builds the Abilene (Internet2) backbone instance (11 nodes, 28 directed arcs; fiber-propagation delays at 5 µs/km from inter-city great-circle distances), sweeps the delay bound Δ, solves the hierarchy and tree ILPs with HiGHS, writes `abilene_results.json`. Reproduces the Abilene case-study section. |
| `scalability.py` | Generates random Erdős–Rényi instances (n = 10…20, \|M\| = 5, Δ = 1.5 · maxSPT), solves the exact hierarchy ILP with HiGHS under a per-solve time cap, writes `scalability_results.json`. Reproduces the scalability table. |
| `gap_dist.py` | Reproduces the tree-vs-hierarchy gap distribution over 120 random Erdős–Rényi instances (n ∈ {10, 12, 14}, 40 seeds each). Writes per-instance results and the summary `gap_dist_summary.json`. Reproduces the gap-distribution section. |

**Shared library (`core/`):**

| File | Description |
|---|---|
| `core/graph.py` | `DirectedGraph`, `Arc`, and `random_graph` utilities: bidirectional-link construction, Dijkstra on arc delays, Erdős–Rényi directed-graph generator. Used by `abilene.py` and `scalability.py`. |
| `core/__init__.py` | Package marker (empty). |

**Solver and audit:**

| File | Description |
|---|---|
| `solve_highs.py` | Exact min-variation multicast ILP (hierarchy or tree) solved via HiGHS. Two-stage lexicographic tie-break (minimize δ, then minimize total delay subject to δ = δ*) fixes a unique delay vector, making the reported Jain index reproducible. Called by `abilene.py` and `scalability.py`. |
| `verify_walks.py` | Post-hoc verification that every proven-optimal per-destination flow F(m, ·) decomposes as a genuine s-to-m walk (plus, optionally, circulations sharing a node with the walk), certifying that no reported optimum contains a detached circulation. |

**Result files:**

| File | Description |
|---|---|
| `abilene_results.json` | Proven-optimal hierarchy vs. tree delay-variation results per delay bound, with Jain fairness indices. |
| `scalability_results.json` | Per-size proven-optimal results and solve statistics. |
| `gap_dist_summary.json` | Per-size summary statistics (mean, median, quartiles, min, max, fraction strictly positive) of the relative gap (δ*_T − δ*_H) / δ*_T over 120 random instances. |

## Requirements

- Python ≥ 3.10
- [HiGHS](https://highs.dev) via `highspy`, and PuLP

```bash
pip install highspy pulp
python abilene.py
python scalability.py 10,12,14,16,18,20
python gap_dist.py 10,12,14 40
```

## Citation

Citation details (DOI) will be added upon acceptance.

## Contact

Corresponding author: Hanan Haj Ahmad — hhajahmed@kfu.edu.sa
