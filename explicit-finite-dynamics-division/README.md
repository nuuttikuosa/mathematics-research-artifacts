# Division of Explicit Finite Dynamical Systems Is NP-Complete

Reproduction materials for the computational checks in Nuutti Kuosa’s paper. The question is whether two given successor tables admit a quotient under synchronous direct product. Section 5 reports three experiments. Each has a corresponding script and reference output here.

## Reproduce the results

Use Python 3.10 or newer. All scripts use the standard library. The reference run used CPython 3.12.14. From the directory containing this README, run

```console
python -B reproduce.py
```

Use `python3` if that is the name of your Python executable. The driver runs five checks in a temporary directory and prints `PASS` for each matching output. It compares the complete text output for the explicit tables and the parsed contents of four JSON reports. A failed check exits with a nonzero status. All inputs are generated locally.

## Experiments in Section 5

| Script | Computation | Reference output |
| --- | --- | --- |
| `explicit_tables.py` | Builds the divisor and three targets as successor tables. Recovers their cycles, enumerates quotient cycle types, and verifies each quotient by constructing its product table. | [Three output lines](expected/explicit_tables.txt) |
| `three_plane_cycles.py` | Builds the duplicated-digit construction for three one-variable formulas. A decoder receiving only the divisor and target cycle histograms recovers 2, 1 and 0 quotient classes. The tautology has 150 active sites and a target with 27,693 distinct periods and 329,236,589,030,252,537,126 states. | [Duplicated-digit results](expected/three_plane_cycles.json) |
| `check_routed.py` | Checks the two-address construction against direct SAT enumeration on seven formulas. The feasible-array counts are 2, 1, 0, 2, 0, 5 and 18. On the three one-variable examples, the separate decoder reconstructs 49,152 partial tuples from each histogram pair. The first target has 58,868 periods and 42,733,726,046,144,284,221,840 states. | [Two-address results](expected/check_routed.json) |

For the explicit-table experiment the divisor has 94 states and the interval controls a weighted sum equal to three times one bit.

| Interval | Target states | Quotient isomorphism classes |
| --- | ---: | ---: |
| [0, 3] | 15,792 | 2 |
| [3, 3] | 11,844 | 1 |
| [1, 2] | 13,160 | 0 |

The two SAT experiments use exact integer cycle multiplicities to represent the large targets. The proof of the asymptotic reduction size is in Section 4 of the paper.

## Reading the reports alongside the paper

| Report term | Meaning in the paper |
| --- | --- |
| `periods` | Cycle lengths. The fields `target_periods` and `B_types` count distinct cycle lengths in the target. |
| `axes` | The number of coordinates, *D*. |
| Partial tuples, `decoded_cells` | A partial tuple describes a face. The decoder counts all faces, including fully specified sites. |
| Axis fibers | Intersections of coordinate lines with the set of active sites. The count `ambient_nonempty_axis_fibers` includes lines meeting just one active site. |
| `wire0`, `wire1` | Copying lines at stage zero (the literal slot varies) and stage one (the test slot varies). |
| `glue` | A stage line joining the two stages for one test and its requested literal on one plane. |
| `link`, `clause` | Plane lines. A `link` tests the signs at two endpoints. A `clause` tests three literal bits. |
| `h` | The alphabet parameter λ in the checksum and routed constructions. |

## Additional checks

`cycle_reference.py` compares the complete set of quotient cycle types with the feasible arrays on 92 seeded small instances and four further probes. `check_checksum.py` checks the checksum construction by enumerating the active sites and solving their interval constraints. [Cycle-level results](expected/exact_check.json) · [Checksum results](expected/check_checksum.json)

To inspect the outputs individually, run the commands below from this directory. The JSON filenames refer to new output in the working directory. The reference files are in `expected/`.

```console
python -B explicit_tables.py
python -B cycle_reference.py
python -B three_plane_cycles.py
python -B check_checksum.py --output check_checksum.json
python -B check_routed.py --blind-decoder three_plane_cycles.py --output check_routed.json
```

The code and reference outputs use the [repository’s CC0 1.0 license](../LICENSE).
