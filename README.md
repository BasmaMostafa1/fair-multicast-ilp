# Fair Multicast Routing — Exact ILP Formulations

Data and solver scripts for the manuscript:

> **ILP Formulation and Exact Solution of Multicast Routing with Fairness**
> Basma Mostafa, Hanan Haj Ahmad, Miklós Molnár
> Submitted to *Mathematics* (MDPI), Special Issue "Advanced Optimization Methods and Applications, 3rd Edition".

Exact flow-based integer linear programming (ILP) formulations for fairness-aware multicast routing, optimizing the maximum end-to-end delay (Δ) and the inter-destination delay variation (δ) as objectives or constraints. Feasible solutions are partial spanning hierarchies, a class that contains partial spanning trees as a special case.

## Contents

| File | Description |
|---|---|
| `abilene.py` | Builds the Abilene (Internet2) backbone instance (11 nodes, 28 directed arcs; fiber-propagation delays at 5 µs/km from inter-city great-circle distances), constructs the ILPs, solves with HiGHS, writes `abilene_results.json`. |
| `abilene_results.json` | Proven-optimal hierarchy vs. tree delay-variation results per delay bound, with Jain fairness indices. |
| `scalability.py` | Generates random Erdős–Rényi instances (n = 10…20), solves the exact ILPs with HiGHS, writes `scalability_results.json`. |
| `scalability_results.json` | Per-size proven-optimal results and solve statistics. |

## Requirements

- Python ≥ 3.10
- [HiGHS](https://highs.dev) via `highspy`

```bash
pip install highspy
python abilene.py
python scalability.py
```

## Citation

Citation details (DOI) will be added upon acceptance.

## Contact

Corresponding author: Hanan Haj Ahmad — hhajahmed@kfu.edu.sa
