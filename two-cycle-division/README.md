# Division by Two Cycles in Binary and Explicit Representations

Reproduction materials for the paper by Nuutti Kuosa. The paper asks whether a given permutation has the form $B=(C_a+C_b)X$, where $C_a+C_b$ is the disjoint union of two cycles and the product advances both coordinates simultaneously. The scripts check small cases in which we can compute the answer in two ways. The time bounds in the paper are proved mathematically.

## Run the checks

From this directory, run

```sh
python -B reproduce.py
```

The checks use the Python standard library and were tested with CPython 3.12. They run offline. The runner compares each script's JSON output with the corresponding file in `expected/`, prints one `PASS` line per match and exits with a nonzero status on failure. To inspect the actual numbers for one check, run the script directly, for example `python -B checks/check_pairs.py`.

## What each check measures

| Paper | Script in `checks/` | Comparison and reference result |
| --- | --- | --- |
| Lemma 2.1, explicit successor tables | `check_direct_blocks.py` | Enumerates all target cycle types for four divisor pairs and quotient sizes 1–3. Solving the block equations agrees with constructing the products directly on all 7,444 targets. There are 19 positive cases. |
| Proposition B.2, finding a binary quotient | `check_pairs.py` | Tests 42,877 targets with at most 24 states against direct permutation products. All 279 positive cases return a witness whose product is the target. The check solves small integer systems by enumeration. |
| Proposition B.2, larger binary cycle lengths | `check_seeded_pairs.py` | On 1,330 seeded cases, compares retaining one representative per ordered output pair with using the full candidate list. Both lists feed the same exact feasibility search and agree on all decisions, including 774 positive cases. |
| Sections 4–5 and Theorem 7.1, weighted tables | `check_weighted_tables.py` | Constructs the numerical row and column margins and counts compatible tables. One copy of the complete graph on four vertices gives 12 tables. Two disjoint copies give 144. Every accepted margin sum is verified with exact integers. |
| Theorem 6.1, graph routing | `check_routing.py` | Recovers and checks a colouring on all 126 edges of the expanded four-vertex example. For the expanded Petersen graph, combining the allowed colours at the ends of each routing gadget gives zero colourings. |
| Lemma 6.2 and Appendix A, address multipliers | `check_addresses.py` | For addresses of one to four bits and four moduli, verifies 120 multipliers, 196 neighbouring address pairs and 740 two-term multisets. The general claim about short integer relations requires the proof in Appendix A. |

Two examples help interpret the outputs. With $A=C_4+C_6$ and $B=2A$, both $C_2$ and $2C_1$ are quotients. The representative search returns $2C_1$. One witness answers the search question, while counting requires distinguishing the two.

In the seeded test, $D=\operatorname{lcm}(a,b)$ and $R=\operatorname{lcm}(D,s)$ for a representative length $s$. The largest observed number of divisors of $R$ is **486**. At that example $D=36$ has just **9** divisors, so Proposition B.2 bounds the number of variables in that block by 9. The weighted-table and routing scripts exercise individual stages of the reductions. Building the final explicit successor tables would require writing out all the states of the constructed permutations.

## Citation and license

Citation metadata is in [CITATION.cff](CITATION.cff). These materials use the repository's [CC0 1.0 license](https://creativecommons.org/publicdomain/zero/1.0/).
