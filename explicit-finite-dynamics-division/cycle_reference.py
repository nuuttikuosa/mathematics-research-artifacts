"""Check Lemma 2 (interval realization) by enumerating integer cycle multiplicities.

Run: python -B cycle_reference.py
Uses only the Python standard library. The interval construction supplies the
expected quotient classes. The independent solver receives only the cycle
lists of A and B and considers divisors of the lcm of the target cycle lengths.
"""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import ceil, floor, gcd, lcm, prod
from pathlib import Path
from random import Random
import json


def divisors(number):
    remaining = number
    values = [1]
    prime = 2
    while prime * prime <= remaining:
        if remaining % prime == 0:
            powers = [1]
            while remaining % prime == 0:
                remaining //= prime
                powers.append(prime * powers[-1])
            values = [v * p for v in values for p in powers]
        prime = 3 if prime == 2 else prime + 2
    if remaining > 1:
        values += [v * remaining for v in values]
    return sorted(values)


def multiply(a, x):
    result = Counter()
    for b, copies in x.items():
        for d in a:
            result[lcm(b, d)] += gcd(b, d) * copies
    return {b: copies for b, copies in result.items() if copies}


def cycle_solutions(a, b):
    quotient_states, rem = divmod(sum(d * c for d, c in b.items()), sum(a))
    if rem:
        return set(), 0, 0, 0
    all_divisors = divisors(lcm(*b))
    periods = []
    upper = []
    rows = {d: [] for d in b}
    for period in all_divisors:
        column = multiply(a, {period: 1})
        bound = min([quotient_states // period] +
                    [b.get(d, 0) // c for d, c in column.items()])
        if bound == 0:
            continue
        index = len(periods)
        periods.append(period)
        upper.append(bound)
        for d, coefficient in column.items():
            rows[d].append((index, coefficient))
    equations = [(row, b[d]) for d, row in rows.items()]

    def propagate(low, high):
        # In each positive-integer row, subtract all other terms' current
        # upper/lower bounds to bound one variable. The former row snapshot
        # remains a valid relaxation even after another variable is tightened.
        changed = True
        while changed:
            changed = False
            for terms, rhs in equations:
                min_rhs = sum(c * low[i] for i, c in terms)
                max_rhs = sum(c * high[i] for i, c in terms)
                if min_rhs > rhs or max_rhs < rhs:
                    return False
                for i, c in terms:
                    new_low = max(low[i], high[i] - (max_rhs - rhs) // c)
                    new_high = min(high[i], low[i] + (rhs - min_rhs) // c)
                    if new_low > new_high:
                        return False
                    if (new_low, new_high) != (low[i], high[i]):
                        low[i], high[i] = new_low, new_high
                        changed = True
        return True

    nodes = 0
    solutions = set()

    def visit(low, high):
        nonlocal nodes
        nodes += 1
        if not propagate(low, high):
            return
        candidates = [i for i in range(len(low)) if low[i] < high[i]]
        if not candidates:
            x = tuple((d, v) for d, v in zip(periods, low) if v)
            assert multiply(a, dict(x)) == b
            assert x not in solutions
            solutions.add(x)
            return
        i = min(candidates, key=lambda j: high[j] - low[j])
        for value in range(low[i], high[i] + 1):
            child_low, child_high = low.copy(), high.copy()
            child_low[i] = child_high[i] = value
            visit(child_low, child_high)

    visit([0] * len(periods), upper)
    return solutions, len(all_divisors), len(periods), nodes


def one_case(groups, seed=None, tight=None, stress=False):
    # The seed, face order and sequence of random draws determine each seeded case.
    top = list(product(*groups))
    faces = list(product(*[(1, *group) for group in groups]))
    if seed is None:
        capacities = [1] * len(top)
        if stress:
            capacities = [1 if v == (5, 11) else 2 for v in top]
    else:
        rng = Random(seed)
        while True:
            capacities = [rng.choice([0, 1, 1, 2]) for v in top]
            if any(capacities):
                break
        den = rng.choice([1, 2, 3])
        point = [Fraction(rng.randint(0, m * den), den) for m in capacities]
    weights = [[(i, prod(v) // prod(e)) for i, v in enumerate(top)
                if all(x == 1 or x == y for x, y in zip(e, v))]
               for e in faces]

    def face_sum(row, x):
        return sum(c * x[i] for i, c in row)

    bounds = []
    for e, row in zip(faces, weights):
        width = face_sum(row, capacities)
        if seed is None or 1 not in e:
            bounds.append((0, width))
        else:
            reference = face_sum(row, point)
            lo, hi = floor(reference), ceil(reference)
            if rng.random() >= tight:
                lo = rng.randint(0, lo)
                hi = rng.randint(hi, width)
            bounds.append((lo, hi))

    def formal_sum(x):
        signed = {}
        for e, row, (lo, hi) in zip(faces, weights, bounds):
            f = face_sum(row, x)
            z = hi - f if e.count(1) % 2 else f - lo
            if z:
                signed[2 * prod(e)] = z
            if hi - lo - z:
                signed[prod(e)] = 2 * (hi - lo - z)
        return signed

    a = [2, *map(prod, groups)]
    b = multiply(a, formal_sum([0] * len(top)))
    assert b and all(isinstance(c, int) and c > 0 for c in b.values())
    if seed is not None:
        assert multiply(a, formal_sum(point)) == b
    expected = set()
    for x in product(*(range(m + 1) for m in capacities)):
        if all(lo <= face_sum(row, x) <= hi
               for row, (lo, hi) in zip(weights, bounds)):
            vector = tuple(sorted(formal_sum(x).items()))
            assert vector not in expected
            expected.add(vector)
    observed, all_periods, live_periods, nodes = cycle_solutions(a, b)
    assert expected == observed, (groups, seed, expected - observed, observed - expected)
    return dict(groups=groups, seed=seed, tight=tight, stress=stress,
                B_states=sum(d * c for d, c in b.items()), count=len(observed),
                lcm_divisors=all_periods, positive_bound_candidates=live_periods,
                search_nodes=nodes)


def main():
    cases = []
    for groups, trials in [(((3, 5), (7, 11)), 40),
                           (((3, 5, 7), (11, 13)), 40),
                           (((3, 5), (7, 11), (13, 17)), 12)]:
        for i in range(trials):
            cases.append(one_case(groups, 1000 * len(groups) + i, [0, 0.5, 0.9][i % 3]))
    extra = [one_case(((3, 5), (7, 11))),
             one_case(((3, 5), (7, 11)), stress=True),
             one_case(((3, 5),)), one_case(((3, 5), (7,)))]
    summary = dict(random_cases=len(cases), default_cases=len(extra),
                   random_quotients=sum(x['count'] for x in cases),
                   random_no_instances=sum(x['count'] == 0 for x in cases),
                   extra_counts=[x['count'] for x in extra],
                   all_search_nodes=sum(x['search_nodes'] for x in cases + extra))
    target = Path(__file__).with_name("exact_check.json")
    target.write_text(json.dumps(dict(summary=summary, random=cases, extra=extra),
                                 indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
