"""Compare prescribed and automatically selected slices of planted star forests."""
from collections import Counter
from math import prod
from random import Random

import sympy as sp

from compare_local import brute, star
from prescribed_slice import splits_with_tau
from tree_arithmetic import Trees
from weighted_split import dot


def inputs():
    rng = Random(20261121)

    def sample():
        polynomial = {(0, 0): 1}
        for _ in range(rng.randint(1, 2)):
            # The original generator samples the coefficient before the exponent.
            count = rng.randint(1, 2)
            exponent = (rng.randint(0, 3), rng.randint(0, 2))
            polynomial[exponent] = count
        return polynomial

    for _ in range(150):
        yield sample(), sample(), (2, 3)
    for a, b in [([1, 1], [1, 1]), ([1, 2, 1], [1, 1]),
                 ([1, 1, 1], [1, 1, 1]), ([1, 0, 1], [1, 2, 1])]:
        yield ({(i,): n for i, n in enumerate(a) if n},
               {(i,): n for i, n in enumerate(b) if n}, (2,))


def check(index, left, right, primes):
    target = Counter()
    for a, n in left.items():
        for b, m in right.items():
            target[tuple(i + j for i, j in zip(a, b))] += n * m
    if target.total() > 10:
        return None

    trees = Trees()
    degrees = {e: prod(p**i for p, i in zip(primes, e)) for e in target}
    ids = {e: trees.node([0] * degree) for e, degree in degrees.items()}
    encoded = Counter({ids[e]: n for e, n in target.items()})
    tuples = Counter({star(degrees[e]): n for e, n in target.items()})
    expected = brute(tuples, star(1), star(1))
    oracle = {(tuple(sorted(map(len, a))), tuple(sorted(map(len, b)))) for a, b in expected}
    unit = trees.node([0])

    def decode_at(tau):
        pairs, excess, choice, mu, grade, second, points = splits_with_tau(
            trees, encoded, unit, unit, tau)

        def signature(side):
            return tuple(sorted(len(trees.children[tree]) for tree, n in side for _ in range(n)))

        actual = {(signature(a), signature(b)) for a, b in pairs}
        assert actual == oracle and actual, (index, tau, actual, oracle)
        return excess, choice, mu, grade, second, points

    excess, choice, mu, grade, second, points = decode_at(1)
    assert (choice > 1) == excess
    cause = ""
    if excess:
        t = sp.Symbol("t")

        def slice_one(factor):
            return sp.Poly(sum(n * t**dot(points[ids[e]], grade) for e, n in factor.items()), t)

        f, g = slice_one(left), slice_one(right)
        x, y = sp.symbols("x y")

        def polynomial(factor):
            return sp.Poly(sum(n * x**i * y**j for (i, j), n in factor.items()), x, y)

        assert polynomial(left).gcd(polynomial(right)).total_degree() == 0
        assert mu == 1 and choice == 2
        if f.is_sqf and g.is_sqf:
            assert f.gcd(g).monic().as_expr() == t + 1
            cause = "shared root"
        else:
            assert index == 72 and f.gcd(g).degree() == 0
            assert sp.factor(f.as_expr()) == (1 + t**2)**2
            assert g.as_expr() == 1 + t**3
            assert all((dot(points[ids[e]], grade), dot(points[ids[e]], second))
                       == (e[0] + 2*e[1], e[1]) for e in target)
            assert len(oracle) == 4
            assert decode_at(2)[:3] == (False, 2, 1)
            cause = "internal repeated root"
            print("Case 72: grades (1,2) and (0,1); edge cap 1; selected tau 2.")
            print("At tau 1, F=(1+t^2)^2, G=1+t^3; tau 1 and tau 2 both match the four ordered oracle splits.")
    signature = tuple(sorted((degrees[e], n) for e, n in target.items()))
    return signature, cause, not any(second) and mu > 1


def main():
    outcomes = []
    for i, (left, right, primes) in enumerate(inputs()):
        outcome = check(i, left, right, primes)
        if outcome is not None:
            outcomes.append(outcome)
    targets = {item[0] for item in outcomes}
    causes = Counter(item[1] for item in outcomes)
    flat = sum(item[2] for item in outcomes)
    assert (len(outcomes), len(targets), causes["shared root"],
            causes["internal repeated root"], flat) == (106, 103, 6, 1, 8)
    print("154 generated; 106 tested; 103 distinct targets; 0 mismatches at tau 1.")
    print("7 excess multiplicities: 6 shared roots, 1 internal repeated root; 8 horizontal slices with edge cap > 1.")


if __name__ == "__main__":
    main()
