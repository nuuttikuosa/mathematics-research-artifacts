# Division by Two Cycles in Binary and Explicit Representations

Reproduction materials for the paper by Nuutti Kuosa. The paper asks whether a given permutation has the form $B=(C_a+C_b)X$, where $C_a+C_b$ is the disjoint union of two cycles and the product advances both coordinates simultaneously. The scripts check small cases in which we can compute the answer in two ways. The time bounds in the paper are proved mathematically.

## Run the checks

From this directory, run

```sh
python -B reproduce.py
```

The checks use the Python standard library and were tested with CPython 3.12. They run offline. The runner compares the output of each of the ten scripts with its reference in `expected/`. It compares JSON values for the six JSON files and report lines for the four text files. A successful run prints ten `PASS` lines. A failure gives a nonzero exit status. To inspect one result, run its script directly, for example `python -B checks/check_mixed_axes.py`.

## What each check measures

| Paper | Script in `checks/` | Comparison and reference result |
| --- | --- | --- |
| Lemma 2.1, explicit successor tables | `check_direct_blocks.py` | Enumerates all target cycle types for four divisor pairs and quotient sizes 1–3. Solving the block equations agrees with constructing the products directly on all 7,444 targets. There are 19 positive cases. |
| Theorem 2.2, shared row and column length | `check_shared_length.py` | Compares counts from the row-only dynamic program against products of successor tables on 4,061 targets, including 3,892 negative cases. It also checks 1,115 positive targets in a larger range and 4,368 prefixes of the variable order. |
| Proposition B.2, finding a binary quotient | `check_pairs.py` | Tests 42,877 targets with at most 24 states against direct permutation products. All 279 positive cases return a witness whose product is the target. The check solves small integer systems by enumeration. |
| Proposition B.2, larger binary cycle lengths | `check_seeded_pairs.py` | On 1,330 seeded cases, compares retaining one representative per ordered output pair with using the full candidate list. Both lists feed the same exact feasibility search and agree on all decisions, including 774 positive cases. |
| Sections 4–5 and Appendix C, weighted tables | `check_weighted_tables.py` | Constructs the numerical row and column margins and counts compatible tables. One copy of the complete graph on four vertices gives 12 tables. Two disjoint copies give 144. Every accepted margin sum is verified with exact integers. |
| Theorem 6.1, graph routing | `check_routing.py` | Recovers and checks a colouring on all 126 edges of the expanded four-vertex example. For the expanded Petersen graph, combining the allowed colours at the ends of each routing gadget gives zero colourings. |
| Lemma 6.2 and Appendix A, address multipliers | `check_addresses.py` | For addresses of one to four bits and four moduli, verifies 120 multipliers, 196 neighbouring address pairs and 740 two-term multisets. The general claim about short integer relations requires the proof in Appendix A. |
| Theorem 7.1, exact covers | `check_exact_covers.py` | For 70 triple families and nine parameter pairs, compares the quotients predicted by the exact covers with all quotient classes found directly from the cycle lists. All 630 sets match. Removing a divisibility condition in a control example increases the count from one to four. |
| Corollary 8.1, a prescribed gcd | `check_gcd_scaling.py` | Enumerates all admissible quotient lengths for 342 choices of axis lengths and prescribed gcd. The scaled and original equations give corresponding solution sets on 270 positive and 4,905 negative inputs. A separate successor-table test checks scaling by a prime greater than the total number of input states. |
| Corollary 8.1, prime powers on both axes | `check_mixed_axes.py` | Enumerates admissible lengths and uses fixed-point counts of iterates to recover all quotients. Across five choices of the axes and their common factor it finds 365 positive and 11,635 negative inputs. It also checks what can go wrong if a prime power is split incorrectly or a row and column can contribute to the same length. |

Two examples help interpret the outputs. With $A=C_4+C_6$ and $B=2A$, both $C_2$ and $2C_1$ are quotients. The representative search returns $2C_1$. One witness answers the search question, while counting requires distinguishing the two.

In the seeded test, $D=\operatorname{lcm}(a,b)$ and $R=\operatorname{lcm}(D,s)$ for a representative length $s$. The largest observed number of divisors of $R$ is **486**. At that example $D=36$ has just **9** divisors, so Proposition B.2 bounds the number of variables in that block by 9. The weighted-table and routing scripts exercise individual stages of the reductions. Building the final explicit successor tables would require writing out all the states of the constructed permutations.

## Citation and license

Citation metadata is in [CITATION.cff](CITATION.cff). These materials use the repository's [CC0 1.0 license](https://creativecommons.org/publicdomain/zero/1.0/).
