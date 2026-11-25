# Division by Two Cycles in Binary and Explicit Representations

Hexagon submission `hexagon:2611.01091`.

Reproduction materials for the paper by Nuutti Kuosa. Given two cycles and a target $B$, can we find a system $X$ for which $(C_a+C_b)X$ is isomorphic to $B$? How many choices of $X$ are there up to isomorphism? The scripts test the paper's constructions on permutations and on explicit targets with trees attached to their cycles. The time bounds in the paper are proved mathematically.

## Run the checks

From this directory, run

```sh
python -B reproduce.py
```

The checks use the Python standard library and were tested with CPython 3.12. They run offline. The runner compares each script's output with its reference in `expected/`, comparing JSON values or text lines as appropriate. A successful run prints 21 `PASS` lines. A failure gives a nonzero exit status. To inspect one result, run its script directly, for example `python -B checks/check_crt_blocks.py`.

## What each check measures

| Paper | Script in `checks/` | Comparison and reference result |
| --- | --- | --- |
| Lemma 2.1, explicit successor tables | `check_direct_blocks.py` | Enumerates all target cycle types for four divisor pairs and quotient sizes 1–3. Solving the block equations agrees with constructing the products directly on all 7,444 targets. There are 19 positive cases. |
| Theorem 2.2, shared row and column length | `check_shared_length.py` | Compares counts from the row-only dynamic program against products of successor tables on 4,061 targets, including 3,892 negative cases. It also checks 1,115 positive targets in a larger range and 2,080 prefixes of the variable order, using the row bound from the reduced cycle lengths. |
| Proposition B.2, finding a binary quotient | `check_pairs.py` | Tests 42,877 targets with at most 24 states against direct permutation products. All 279 positive cases return a witness whose product is the target. The check solves small integer systems by enumeration. |
| Proposition B.2, larger binary cycle lengths | `check_seeded_pairs.py` | On 1,330 seeded cases, compares retaining one representative per ordered output pair with using the full candidate list. Both lists feed the same exact feasibility search and agree on all decisions, including 774 positive cases. |
| Sections 4–5, weighted tables | `check_weighted_tables.py` | Constructs the numerical row and column margins and counts compatible tables. One copy of the complete graph on four vertices gives 12 tables. Two disjoint copies give 144. Every accepted margin sum is verified with exact integers. |
| Sections 6 and 8, graph routing and the counting factor | `check_routing.py` | Enumerates the allowed colourings of each diamond for all boundary colours. Equal colours admit two extensions and unequal colours admit none. The expanded four-vertex graph has 20 diamonds and $6\cdot 2^{20}=6{,}291{,}456$ labelled colourings. The expanded Petersen graph has none. |
| Lemma 6.2 and Appendix A, address multipliers | `check_addresses.py` | For addresses of one to four bits and four moduli, verifies 120 multipliers, 196 neighbouring address pairs and 740 two-term multisets. The general claim about short integer relations requires the proof in Appendix A. |
| Theorem 7.1, exact covers | `check_exact_covers.py` | For 70 triple families and nine parameter pairs, compares the quotients predicted by the exact covers with all quotient classes found directly from the cycle lists. All 630 sets match. Removing a divisibility condition in a control example increases the count from one to four. |
| Section 8, a prescribed gcd | `check_gcd_scaling.py` | Enumerates all admissible quotient lengths for 342 choices of axis lengths and prescribed gcd. The scaled and original equations give corresponding solution sets on 270 positive and 4,905 negative inputs. A separate successor-table test checks scaling by a prime greater than the total number of input states. |
| Section 8, prime powers on both axes | `check_mixed_axes.py` | Enumerates admissible lengths and uses fixed-point counts of iterates to recover all quotients. Across five choices of the axes and their common factor it finds 365 positive and 11,635 negative inputs. It also checks what can go wrong if a prime power is split incorrectly or a row and column can contribute to the same length. |
| Section 8, consecutive divisor lengths | `check_consecutive_divisors.py` | Compares complete quotient sets on 2,300 table instances, 168 odd-weight subset-sum instances and 12 exact-cover reductions. The prescribed gcd takes the values 1, 6, 12 and 72. Cases include zero, one and multiple solutions. |

Two examples help interpret the outputs. With $A=C_4+C_6$ and $B=2A$, both $C_2$ and $2C_1$ are quotients. The representative search returns $2C_1$. One witness answers the search question, while counting requires distinguishing the two.

In the seeded test, $D=\operatorname{lcm}(a,b)$ and $R=\operatorname{lcm}(D,s)$ for a representative length $s$. The largest observed number of divisors of $R$ is **486**. Here $a=9$, $b=12$ and $D=36$ has 9 divisors. Removing the common factor 3 gives 3 and 4. Proposition B.2 therefore bounds the number of variables in the block by $\tau(3)\tau(4)=2\cdot3=6$. The check tests this reduced bound in every sampled block. The weighted-table and routing scripts exercise individual stages of the reductions. Building the final explicit successor tables would require writing out all the states of the constructed permutations.

## Normal forms, transient targets and a second quotient

| Paper | Script in `checks/` | What to look for |
| --- | --- | --- |
| Theorem 2.2, order of processing the columns | `check_split_quota_dp.py` | Runs the dynamic program in three orders on 24,129 inputs. It compares the counts with direct products and verifies the recovered witnesses. A length with positive quotas on both sides occurs first, in the middle and last. |
| Theorem A, consecutive normal form | `check_normal_form_successor_tables.py` | Checks the bijection on 2,433 inputs using explicit successor tables. Cases include equal axes, a unit axis, collisions and targets with no quotient. |
| Theorem A, inverse at a shared length | `check_collision_inverse.py` | Tests whether a quotient of the transformed target comes from exactly one original quotient. The same numerical length can receive positive contributions from both projections. |
| Corollary 8.2, transient targets | `check_transient_targets.py` | Builds products of successor tables and compares the classes found with the prediction from periodic words of rooted trees. |
| Corollary 8.2, rotations and sizes | `check_rotations_and_size.py` | Changes the starting vertices of the cyclic words, tests a reflection and counts states when a single attached tree has up to 8,192 leaves. |
| Proposition C.3, detecting a second quotient | `check_second_quotient.py` | Compares decisions with direct enumeration on 54,677 targets. Two quotient classes can arise from different representative vectors or from different ways to realise one vector. |
| Proposition C.3, recovering distinct quotients | `check_lifts_and_cost.py` | Verifies the possible periods for 12,755 ordered pairs of image lengths and checks a second witness with a 1,025-bit multiplicity. |
| Proposition C.3, the two-by-two table | `check_crt_blocks.py` | Compares the congruence calculation with enumeration on 6,216 small cases. In a separate family with $2^{170}+1$ vectors the script computes the interval, stride and two witnesses. |
| Proposition C.3, factoring by decision queries | `check_factoring_queries.py` | Tests both the uniqueness and nonuniqueness oracles for $1\le g\le512$. A 488-bit example uses a simulated threshold oracle. |
| Corollary 8.2 and Proposition C.3, further cross-checks | `check_final_crosschecks.py` | Checks the gcd identity and the admissible word periods on 20,000 seeded triples, compares the tables with exhaustive search on 694 targets, and tests the query strategy on 3,006 integers using a simulated oracle. |

## Citation and license

Citation metadata is in [CITATION.cff](CITATION.cff). These materials use the repository's [CC0 1.0 license](https://creativecommons.org/publicdomain/zero/1.0/).
