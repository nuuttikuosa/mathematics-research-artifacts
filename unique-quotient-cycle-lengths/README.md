# Checking the finite examples

These files reproduce the table for the first four through eight primes
in *How Many Cycle Lengths Can a Unique Quotient Require?* They verify
the exhibited optima by exact integer arithmetic. The two commands can
be run separately from this directory. They read the supplied inputs
and print their results.

## The five individual examples

This step needs Python 3.10 or later and its standard library. Use the
normal interpreter with assertions enabled. The checks were tested on
Python 3.12.14.

~~~console
python -B verify.py
~~~

Expected output:

~~~text
k=4  support=12  size=246  solutions=1
k=5  support=21  size=2051  solutions=1
k=6  support=30  size=20185  solutions=1
k=7  support=43  size=343463  solutions=1  profiles=486  prime=17
k=8  support=58  size=6526232  solutions=1  profiles=2640  prime=19
~~~

Each of the three files named `record_cert_*.json` contains the list of
primes, the target size and marginals, the witness as pairs of length and
multiplicity, and a search tree. The tree verifier recomputes each bound
and checks that each remaining multiplicity has a branch. The three trees
have 17, 95 and 893 nodes. Every tree allows length 1. A single solution
at the proposed size, itself free of fixed points, proves that a smaller
solution is impossible because it could be padded to that size.

The file `data/lifts.json` gives the added lengths, each with multiplicity
one, for the two larger witnesses. To recover all the lengths for seven
primes, multiply every length in the six-prime witness by 17 and add the
first list. For eight primes, multiply every length of this seven-prime
witness by 19 and add the second list. The verifier enumerates all the
allowed marginal vectors of the added part. For each vector,
`prime_split.py` counts all representations by splitting on one prime.

## The exact value of the four-prime maximum

This step also needs NumPy. It was tested with NumPy 2.3.5. The calculation
uses dense arrays of 9,261,000 integers per array, so allow several
hundred megabytes of memory.

~~~console
python -B all_four.py
~~~

Expected output:

~~~text
k=4  vectors=9261000  maximum=12  attainers=98
~~~

The dynamic programme minimises size and then support at every admissible
marginal vector. The four-prime maximum is therefore exact. For five
through eight primes the certificate checks give lower bounds for the
maximum. They establish uniqueness for each displayed example.

The scripts and data use the [repository’s CC0 1.0 license](../LICENSE).
