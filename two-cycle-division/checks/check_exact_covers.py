"""Section 7: compare whole sets of quotient classes, not only their sizes.

Run with Python 3.9+ and no third-party packages. The cycle-list oracle gets
only the divisor lengths and target. It knows neither the triples nor which
output lengths the reduction calls rows.
"""
from collections import Counter, defaultdict
from itertools import combinations
from math import gcd, isqrt, lcm, prod


def subsets(size):
    for count in range(size + 1):
        yield from combinations(range(size), count)


def exact_covers(k, triples):
    ground = list(range(3 * k))
    return {indices for indices in combinations(range(len(triples)), k)
            if sorted(j for i in indices for j in triples[i]) == ground}


def signature(multiplicities):
    return tuple(sorted((period, count) for period, count in multiplicities.items() if count))


def encode(k, triples, h, q):
    s = len(triples)
    special_divisor = (h * q, h * (q + 1))
    if not k:
        return special_divisor, Counter(special_divisor), (), "empty"
    if s < k:
        return special_divisor, Counter({1: 1}), (), "too few"
    radix = s + 2
    multiplier = q * radix
    weights = [h * (1 + multiplier * sum(radix**j for j in triple)) for triple in triples]
    wanted = h * (k + multiplier * sum(radix**j for j in range(3 * k)))
    total = sum(weights)
    if wanted > total:
        return special_divisor, Counter({1: 1}), (), "weight deficit"
    b = (q + 1) * prod(weights)
    target = Counter({h * q * u: h * q for u in weights})
    target[h * b] += h * q * wanted
    if wanted < total:
        target[h * q * b] += h * (total - wanted)
    return (h * q, h * b), target, weights, "full"


def all_divisors(n):
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            yield d
            if d != n // d:
                yield n // d


def all_quotients(divisor, target):
    a, b = divisor
    # lcm(a,t) must be an observed period v, hence t=(v/a)*gcd(a,t).
    candidates = {(v // a) * g for g in all_divisors(a) for v in target if v % a == 0}
    lengths = sorted(target, key=lambda v: (target[v], v))
    groups = defaultdict(list)
    for t in sorted(candidates):
        column = Counter()
        for d in divisor:
            column[lcm(d, t)] += gcd(d, t)
        if all(value <= target[v] for v, value in column.items()):
            pivot = next(v for v in lengths if v in column)
            groups[pivot].append((t, column))

    answers = set()
    remaining = target.copy()
    chosen = {}

    def visit(i, j):
        if i == len(lengths):
            assert not any(remaining.values())
            answers.add(signature(chosen))
            return
        v = lengths[i]
        if j == len(groups[v]):
            if remaining[v] == 0:
                visit(i + 1, 0)
            return
        t, column = groups[v][j]
        limit = min(remaining[u] // entry for u, entry in column.items())
        # The last variable in a pivot equation must fill it exactly.
        if j + 1 == len(groups[v]):
            value, remainder = divmod(remaining[v], column[v])
            options = (value,) if remainder == 0 and value <= limit else ()
        else:
            options = range(limit + 1)
        for count in options:
            chosen[t] = count
            for u, entry in column.items():
                remaining[u] -= count * entry
            visit(i, j + 1)
            for u, entry in column.items():
                remaining[u] += count * entry
        chosen.pop(t, None)

    visit(0, 0)
    return answers


def predicted_class(indices, weights, h, q):
    result = Counter()
    for i, u in enumerate(weights):
        if i in indices:
            result[h * u] += q
        else:
            result[h * q * u] += 1
    return signature(result)


def main():
    # The six sets come in three complementary pairs. Subfamilies can have
    # zero, one, two, or three distinct covers. Three extra tests have k=3.
    named = [(0, 1, 2), (3, 4, 5), (0, 1, 3), (2, 4, 5), (0, 2, 4), (1, 3, 5)]
    families = [(0, ()), (1, ()), (1, ((0, 1, 2),))]
    families.extend((2, tuple(named[i] for i in selected)) for selected in subsets(6))
    families.extend([(3, tuple(named) + ((6, 7, 8),)),
                     (3, tuple(named[:2]) + ((5, 7, 8),)),
                     (3, tuple(named[2:]))])
    parameters = [(1, 2), (1, 3), (2, 3), (3, 2), (4, 3),
                  (6, 5), (8, 3), (9, 2), (12, 5)]
    histogram = Counter()
    branches = Counter()
    classes_checked = 0
    for k, triples in families:
        covers = exact_covers(k, triples)
        histogram[len(covers)] += 1
        for h, q in parameters:
            divisor, target, weights, branch = encode(k, triples, h, q)
            branches[branch] += 1
            if branch == "empty":
                expected = {((1, 1),)}
            else:
                expected = {predicted_class(indices, weights, h, q) for indices in covers}
            assert len(expected) == len(covers)
            assert gcd(*divisor) == h
            got = all_quotients(divisor, target)
            assert got == expected, (k, triples, h, q, expected, got)
            classes_checked += len(got)
    print("Families by number of exact covers:", dict(sorted(histogram.items())))
    print("(h,q) pairs:", parameters)
    print("Matched quotient-class sets:", sum(branches.values()), dict(sorted(branches.items())))
    print("Quotient classes across all matched instances:", classes_checked)

    # Challenge the oracle with a known failure when protection is removed.
    assert all_quotients((3, 4), Counter({3: 3, 4: 3})) == {((1, 3),)}
    unprotected = all_quotients((6, 8), Counter({6: 6, 8: 6}))
    assert unprotected == {signature({1: 6 - 2 * j, 2: j}) for j in range(4)}
    print("Unprotected scaling control: 1 ->", len(unprotected), "classes")
    print("PASS")


if __name__ == "__main__":
    main()
