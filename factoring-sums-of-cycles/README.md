# The Complexity of Factoring Sums of Cycles

Computational checks for Nuutti Kuosa's paper. The paper determines the decision and counting complexity of nontrivial direct-product factorization of finite permutation dynamical systems. These seven separately implemented finite checks exercise the algebraic construction, SAT correspondence, routing geometry, size parameters and characteristic-polynomial comparison used in the paper.

The general results are proved in the paper. Passing these finite checks is supporting evidence, not a proof of the unbounded statements or formal verification. This directory contains scripts and reproduction files only; the manuscript is distributed separately.

## Run the checks

Use **Python 3.10 or later**. All dependencies are in the Python standard library: no installation, data download, network connection or GPU is required.

From this directory, on Windows PowerShell:

```powershell
py -3 -B reproduce.py
```

On Linux or macOS:

```sh
python3 -B reproduce.py
```

The runner also works when invoked by its full path from another directory. A successful run prints seven `PASS` lines and exits with status zero. It stops with a nonzero exit status if a checker assertion fails, its process fails, its output is missing, or any deterministic result differs from the reference. The preparation run on CPython 3.12.14 took about 30 seconds; runtime depends on the machine, and the block suite is the largest part.

To select a suite or retain the raw results:

```sh
python3 -B reproduce.py --suite algebra
python3 -B reproduce.py --suite network --suite blocks
python3 -B reproduce.py --output-dir run-results
```

The default timeout is 600 seconds per suite; `--timeout 1200` raises it to 20 minutes. The runner copies each checker into a fresh temporary directory and executes it with Python's isolated mode, without importing the other checkers. Generated files do not modify the scripts or references. `--output-dir` optionally preserves raw JSON and stdout under the chosen directory, with a different filename for each suite.

## What is checked

| Suite | Procedure and scope | Reference result |
| --- | --- | --- |
| `algebra` | Validate an exact cycle-division oracle against complete small systems; compare interval quotients with feasible integer tables; enumerate all factor pairs of small two-layer and prescribed-marker targets. | 466 oracle equations; 79 interval instances with 85 quotients, including an integrally infeasible case; 648 two-layer targets; 163 marker targets. Includes the paper's 256-state example and the strict-threshold counterexample. |
| `sat` | Independently construct geometric SAT instances, reconstruct occupied coordinate lines and copying components, and compare component-bit solutions with satisfying assignments. | 283 formulas: 271 exhaustive component-bit comparisons and 12 larger assignment-extension checks. The latter 12 do not exhaust all arbitrary component bits. |
| `parameters` | Use integer arithmetic to check the specified radix, word dimension, time dimension and prime-group counts, including boundary values. The code's `R` denotes the paper's size parameter Lambda. | 10,575 parameter choices, including values adjacent to powers and an example with 1,001 decimal digits. This finite check does not establish the asymptotic size bound. |
| `network` | Construct generic permutation networks with every signal active; verify composition, vertex separation, occupied copying lines and coloring. | 225 permutations: all 150 one-digit permutations for radices 3, 4 and 5, plus 75 structured or seeded cases. Checks 130,070 active vertices and 126,502 edge occurrences. Final literal blocks and SAT tests are outside this suite. |
| `blocks` | Check final literal blocks and exit-layer occupancy, including coincident intermediate vertices and shared edges. | 66,211 cases: every partition of every nonempty active prefix for the listed paths with at most 16 words, plus selected partitions of larger paths. This is a local geometry check, not a full SAT census. |
| `mass` | Construct a two-dimensional signed interval example and compare per-length partial-sum maxima with triangle-inequality and length-scaled bounds. | Uncollected absolute mass 625,536 and final mass 565,806. Both the chosen-order and all-order maximum coefficient prefixes are 864; the maximum triangle bound is 1,029. |
| `associates` | Use exact characteristic-polynomial coefficients to test rational root twists, and exact rational root phases for a second rotation search. | All 28,628 nonempty cycle types through 30 states: 683 admissible negative twists, with no distinct-type collision. The phase search covers 271 types through 12 states and 9,081 rotations. Positive-characteristic controls are also checked. |

The randomized selections in `algebra`, `sat` and `network` use the fixed seed `261126`. Full per-case records, rather than just the summary totals, are compared where the checkers produce them. Counts from different suites measure different objects and should not be added into a single coverage total.

## Files and comparison rules

- `checks/<suite>/` contains one standalone checker per suite. Each can also be run directly with Python; direct execution writes its JSON result beside that script. Separate directories avoid collisions between checkers that both use `results.json`.
- `expected/<suite>.json` contains the reference deterministic result. These files include the detailed case records generated by the checkers.
- `reproduce.py` runs and compares the checks; `CITATION.cff` records citation metadata.

The runner excludes only two kinds of top-level run metadata from the cross-machine reference comparison. For `algebra` and `sat`, `seconds` must be a finite nonnegative elapsed-time value. For `network` and `associates`, the emitted `source_sha256` must equal the SHA-256 of the file actually executed. That source hash is validated separately because checkout line endings can change file bytes. Every other JSON field, including notes, scope descriptions, totals, Boolean results and nested case records, is compared exactly. No reference file is regenerated by `reproduce.py`.

The software and data in this directory use the repository's [CC0 1.0 license](../LICENSE). See the [repository home](https://github.com/nuuttikuosa/mathematics-research-artifacts) for the other papers' materials. These checks were developed with AI assistance under the contributor's direction; they should not be read as a claim of human peer review.
