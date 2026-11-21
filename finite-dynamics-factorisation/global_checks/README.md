# Global experiments

From the repository root, run `python -B global_checks/run_global.py` after installing the repository requirements. The script constructs its inputs from scratch, checks exact equalities and prints the summary recorded in `global_checks/expected.txt`. The root reproduction driver also checks the output. Use ordinary Python with assertions enabled.

### The guide budgets

`guide_budget.py` works directly with successor tables. It evaluates every labelled map on one through five states and every ordered pair of endpoints with the same terminal image. It tests the product of the floored width ratios against the square of the terminal width. A second check generates all ordered products of two connected labelled factors on two or three states. The planted factor coordinates give an independent count of active predecessors at each minimum-width guide step. Both checks use only the Python standard library.

### The factor-pair catalogue

`run_global.py` enumerates connected factor types by inspecting all labelled successor tables of each size from two to six. It takes their direct products whenever the target has at most twelve states. Six is sufficient because both factors have at least two states. Each target is stored up to isomorphism, with the complete set of ordered pairs of factor types. The catalogue covers all targets in that range admitting a qualifying decomposition, including disconnected targets. It has 203 targets and 428 ordered pairs.

The phase-grid and paired-guide prototype in `search.py` reconstructs the pairs using exhaustive local splitting. Its tree operations and exact forest divisions are in `fdds.py`. Each non-unit local call is cross-checked against the weighted splitter in `../local_checks/`. Correspondence between the two tree representations is established recursively from their children. Repeated inputs use a cached weighted split and are still compared on every occurrence. There are 17,343 such calls on 206 distinct inputs. Each of these admits a local split. The separate 407-call local suite includes negative cases. This catalogue runs finite, exhaustive comparisons. The manuscript gives the bit bound for the arithmetic algorithm.

### Rejection after exact division

`late_anchor.py` constructs an explicit 420-state target and verifies its two planted factors by direct multiplication. In the phase-grid search, it fixes the two guide endpoints and follows the designated local choices. Its rank-four pivot test recovers the left axis for the seed `(1,0,0,0)`, computes a positive forest quotient and verifies the product. The required right anchor is missing. A different split at that node extends to the planted 21- and 20-state factors.

`late_family.py` constructs four further targets of 60, 360, 2,880 and 34,560 states. It checks selected nodes with 2, 4, 8 and 16 candidate axes rejected by the right-anchor condition after positive exact division. It prints the numbers of signed coordinate orderings tested in its width check. Only the small catalogue explores all potential witnesses for its targets.
