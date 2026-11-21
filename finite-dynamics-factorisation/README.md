# Finding Two Connected Factors of a Finite Dynamical System in Quasipolynomial Time

Computational material for Nuutti Kuosa's paper. Given a finite successor table, the paper asks for two connected systems with at least two states each whose synchronous product is isomorphic to the input. The experiments below compare local and global search procedures with independently generated answers, and exercise two boundary cases in the proof.

This directory contains one article's materials in a multi-article repository. Run the commands from the directory containing this README. The software and generated data in this directory use the license at the root of the mathematics-research-artifacts repository. [Repository home](https://github.com/nuuttikuosa/mathematics-research-artifacts).

## Reproduce the results

The reference environment uses CPython 3.12.14, SymPy 1.14.0 and mpmath 1.3.0. From this directory in a fresh **Windows PowerShell** session, run

~~~powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -B reproduce.py
~~~

On **Linux or macOS**, run

~~~sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -B reproduce.py
~~~

With the dependencies already installed, the single command is **python -B reproduce.py**. The driver runs four suites, prints their outputs and compares them with the reference logs. A failed assertion, nonzero process exit or output difference stops the run. No input download, GPU or network access is needed after installing the dependencies. All test instances are generated locally. The specialization suite uses the fixed seed 20261121.

To run one suite, append **--suite anchors**, **--suite local**, **--suite specialization** or **--suite global**. The anchors suite uses only the Python standard library.

## What each suite checks

| Suite | Input and procedure | Expected result |
| --- | --- | --- |
| **anchors** | Every ordered pair of equal positive height among rooted, unordered trees of at most eight vertices per tree. At each local sibling minimum, multiply each candidate cut by the cut at the matching vertex on a full path in the opposite tree. | 9,633 tree pairs and 108,748 groups. All ordinary-cut products are observed and separate distinct candidates. At radius two, 136 groups contain distinct types. |
| **local** | Compare the two-grade local splitter with direct enumeration of ordered forest pairs. The examples include stars, trees of heights two through five and three targeted polynomial instances. | The complete sets agree on 407 calls covering 347 distinct target-and-anchor triples. There are 153 negative calls. Of 224 calls entering the residue calculation, eight are negative and seven use a specialization greater than one. |
| **specialization** | Generate positive polynomials representing star forests. Compare decoding at the prescribed value one against direct enumeration with a unit star on both axes. Classify multiplicity excesses using the two original polynomial factors. | Of 154 generated inputs, 106 pass the mass cutoff of ten (103 distinct targets). All agree. In six excess cases two factors acquire a shared root. In the seventh, generated case 72, one factor acquires its own repeated roots. |
| **global** | Check predecessor-path budgets against labelled successor tables and true products, cross-check local outputs in a catalogue of positive global instances, and trace late rejections. | 3,413 labelled maps of sizes one through five. The positive catalogue through 12 states has 203 targets and 428 ordered pairs. All 17,343 recorded nonunit local calls agree. A selected late pivot on a 420-state target passes exact division and multiplication before failing the opposite anchor. |

### Two local examples

In the anchor experiment only children reaching the requested radius enter a local comparison. Once a child is chosen, the ordinary cut retains every descendant within that radius, including early leaves. The eight-vertex example compares two sibling trees: a two-edge path and the same path with an additional leaf at its root. At radius two the ordinary cut distinguishes them. The alternative, pruned operation merges them. Both counts are printed for comparison.

For generated case 72, the two planted factors are **F = 1 + 2x² + x²y** and **G = 1 + x³**. At the first specialization F becomes (1 + t²)², which is coprime to 1 + t³. With the actual specialized multiplicities, the decoder recovers all four ordered local splits at both one and two. The fast specialization rule chooses two, where the largest multiplicity meets the edge cap of one. The four splits include the unit on either axis.

### Reading the global counts

An ordered pair records the left and right isomorphism types separately. The catalogue search advances its stack using exhaustive local splitting. At every nonunit call, a separate weighted decoder supplies a complete set of splits for comparison. The search checks cyclic predecessor inclusion and the mark-size and degree bounds before pushing an extension. Each accepted extension is pushed separately. The positive catalogue tests reconstruction and completeness on targets produced by direct multiplication. The local suite separately supplies negative targets and anchor choices. The paper establishes the worst-case bit bound for its specified algorithm.

## Files

- **reproduce.py** runs the suites. The local reference logs are in **expected/**. The global driver and its reference log are in **global_checks/**.
- **local_checks/tree_arithmetic.py**, **residues.py** and **weighted_split.py** implement tree operations, finite labels, projected multiplicity spaces and decoding. The separate oracle in **compare_local.py** works with nested tuples of children and constructs products directly.
- **local_checks/prescribed_slice.py** exposes the actual-multiplicity decoder for a specified slice. **check_specialization.py** generates the seeded tests and checks both sources of multiplicity excess.
- **local_checks/check_anchor_cuts.py** constructs its own catalogue and physical tree products. It checks the combinatorial hypotheses used in the fixed-anchor argument.
- **CITATION.cff** gives the citation metadata.
