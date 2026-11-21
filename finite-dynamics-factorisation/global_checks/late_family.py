"""Check a family of distinct late anchor failures on explicit successor tables.

Run: python -B late_family.py
Uses only the Python standard library. Tests selected global nodes, all the
stated star-forest divisions there, and signed separating grades on the four
finite instances. It does not enumerate the complete global search trees.
"""
from collections import Counter
from itertools import permutations, product
from math import prod


def table_product(left, right):
    return [u * len(right) + v for u in left for v in right]


def single_table(degree):
    return [0] * degree + [1]


def subset_stars(values):
    """Product of (S_1 + S_v) for v in values, with full multiplicities."""
    return Counter(prod(v for v, take in zip(values, selection) if take)
                   for selection in product((False, True), repeat=len(values)))


def multiply(left, right):
    out = Counter()
    for p, a in left.items():
        for q, b in right.items():
            out[p*q] += a*b
    return out


def quotient_by_stars(target, divisor):
    """Positive quotient by stripping the smallest residual degree each time."""
    remaining = target.copy()
    smallest = min(divisor)
    result = Counter()
    while remaining:
        observed = min(remaining)
        value, remainder = divmod(observed, smallest)
        copies, residual_copies = divmod(remaining[observed], divisor[smallest])
        if remainder or residual_copies:
            return None
        used = Counter({value*p: copies*n for p, n in divisor.items()})
        if not used <= remaining:
            return None
        result[value] += copies
        remaining -= used
    return result


def width_tests(target, primes, candidates):
    # The grades are integral on prime exponents. Mixed radix makes the full
    # rectangular support injective for every reordering and every sign choice.
    coords = {(2**i * prod(p**bit for p, bit in zip(primes, bits))): (i, *bits)
              for i in range(4) for bits in product((0, 1), repeat=len(primes))}
    assert set(coords) == set(target)
    trials = 0
    for order in permutations(range(len(primes))):
        for signs in product((-1, 1), repeat=len(primes)+1):
            grade = (signs[0], *(s * (4 << j) for j, s in zip(order, signs[1:])))
            positions = {star: sum(g*k for g, k in zip(grade, coordinate))
                         for star, coordinate in coords.items()}
            assert len(set(positions.values())) == len(target)
            widths = set()
            for subset, axis in candidates:
                powers = [positions[d] for d in axis]
                span = max(powers) - min(powers)
                assert span == abs(grade[0]) + sum(abs(grade[j+1]) for j in subset)
                widths.add(span)
            assert len(widths) == len(candidates)
            trials += 1
    return trials


def check(primes):
    x_table = single_table(2)
    for prime in primes:
        x_table = table_product(x_table, single_table(prime))
    y_table = single_table(4)
    table = table_product(x_table, y_table)
    n = len(table)
    assert n == 15*prod(p+1 for p in primes)
    assert all(table[table[v]] == 0 for v in range(n)) and table[0] == 0
    indegree = Counter(table)
    root_children = [i for i, v in enumerate(table) if v == 0]
    target = Counter(indegree[i] for i in root_children if indegree[i])
    assert target == subset_stars((2, 4, *primes))
    assert target.total() == 4 * 2**len(primes)
    assert all(c == 1 for c in target.values())

    # At the first level take (root_X, degree-one Y child) and
    # (degree-two X child, root_Y). Either physical preimage is an endpoint.
    x_degrees = Counter(x_table)
    a, b = x_degrees[0], len(y_table)-1
    assert a == 2*prod(primes) and indegree[0] == a*b
    x_child = next(i for i,v in enumerate(x_table) if v == 0 and x_degrees[i] == 2)
    parents = (1, x_child*len(y_table))
    endpoints = tuple(table.index(i) for i in parents)
    assert tuple(table[v] for v in endpoints) == parents
    assert all(table[v] == 0 for v in parents)
    assert tuple(indegree[v] for v in parents) == (a, 2*b)
    assert indegree[parents[0]]*indegree[parents[1]] == 2*indegree[0]
    # Previous split sizes 1 and 2 yield old marks a and b respectively.
    assert (indegree[parents[0]]//1, indegree[parents[1]]//2) == (a,b)
    c0, d0 = indegree[parents[0]]//a, indegree[parents[1]]//b
    assert (c0,d0) == (1,2)
    budgets = (a*2*b, (n//a)*(n//(2*b)), 1)
    assert 2 <= budgets[0] and target.total() <= budgets[1] and prod(budgets) <= n*n

    # The original X,Y supply a full witness on this branch. Read their
    # positive root-child indegrees from the two original transition tables.
    true_c = Counter(Counter(y_table)[j] for j,v in enumerate(y_table) if v == 0)
    true_d = Counter(x_degrees[j] for j,v in enumerate(x_table) if v == 0)
    del true_c[0], true_d[0]
    assert true_c == subset_stars((4,)) and true_d == subset_stars((2,*primes))
    assert c0 in true_c and d0 in true_d and multiply(true_c,true_d) == target
    # A type on the left is eligible if its product with the right anchor
    # occurs in the target. Here every one of the designated wrong axes passes.
    candidates = []
    for selection in product((False, True), repeat=len(primes)):
        subset = tuple(i for i, take in enumerate(selection) if take)
        left = subset_stars((2, *(p for i,p in enumerate(primes) if selection[i])))
        right = subset_stars((4, *(p for i,p in enumerate(primes) if not selection[i])))
        assert c0 in left and all(d*d0 in target for d in left)
        assert multiply(left, right) == target
        assert quotient_by_stars(target, left) == right and d0 not in right
        candidates.append((subset, left))
    assert len(candidates) == 2**len(primes)
    grades = width_tests(target, primes, candidates)
    print(f"{len(primes)} {n} {a} {endpoints} {target.total()} {budgets[1]} "
          f"{len(candidates)} {grades}")


if __name__ == "__main__":
    print("r N old_X (guide_Y,guide_X) local_mass strict_budget rejected_axes signed_grades")
    for r in range(1,5):
        check((3,5,7,11)[:r])
    print("All assertions passed.")
