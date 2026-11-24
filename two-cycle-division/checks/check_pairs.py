"""Finite checks for binary division by two arbitrary cycles. Python >= 3.9.

Run: python -B check_pairs.py
The finite representative solver enumerates bounded integer columns instead of
using Lenstra's algorithm. The oracle constructs actual permutation products.
"""
import json
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import gcd, lcm


@lru_cache(None)
def partitions(n, first=1):
    if not n:
        return [()]
    return [(i,) + p for i in range(first, n + 1) for p in partitions(n-i, i)]


def successor(cycles):
    out = []
    for n in cycles:
        offset = len(out)
        out.extend(offset + (j+1) % n for j in range(n))
    return out


def explicit_product(a, b, quotient):
    s, t = successor((a, b)), successor(quotient)
    width = len(t)
    transition = [v * width + w for v in s for w in t]
    seen, parts = set(), []
    for x in range(len(transition)):
        if x in seen:
            continue
        length = 0
        while x not in seen:
            seen.add(x)
            length += 1
            x = transition[x]
        parts.append(length)
    return tuple(sorted(parts))


def least(c, u):
    return gcd(u, pow(u // c, u.bit_length(), u))


def representative(a, b, u, v):
    if u % a or v % b:
        return None
    t = lcm(least(a, u), least(b, v))
    return t if lcm(a, t) == u and lcm(b, t) == v else None


def solve_block(columns, demand):
    def visit(i):
        if i == len(columns):
            return () if not any(demand.values()) else None
        t, col = columns[i]
        limit = min(demand[u] // c for u, c in col.items())
        for copies in range(limit + 1):
            tail = visit(i+1)
            if tail is not None:
                return (t,) * copies + tail
            if copies < limit:
                demand.subtract(col)
        for u, c in col.items():
            demand[u] += limit * c
        return None

    return visit(0)


def solve(a, b, target):
    counts = Counter(target)
    d = lcm(a, b)
    blocks, columns = defaultdict(Counter), defaultdict(dict)
    for u, multiplicity in counts.items():
        blocks[lcm(d, u)][u] = multiplicity
    for u, v in product(counts, repeat=2):
        t = representative(a, b, u, v)
        if t is None:
            continue
        key = lcm(d, t)
        assert key == lcm(d, u) == lcm(d, v)
        col = Counter({u: gcd(t, a)})
        col[v] += gcd(t, b)
        assert t not in columns[key]
        columns[key][t] = col
    answer = []
    for key, demand in blocks.items():
        assert len(columns[key]) <= sum(d % i == 0 for i in range(1, d+1))
        part = solve_block(list(columns[key].items()), demand.copy())
        if part is None:
            return None
        answer.extend(part)
    return tuple(sorted(answer))


def main():
    n_pairs, diagonals, reduced = 0, 0, 0
    for a, b in combinations_with_replacement(range(1, 17), 2):
        pairs = defaultdict(list)
        for t in range(1, 121):
            pairs[(lcm(a, t), lcm(b, t))].append(t)
        for (u, v), fiber in pairs.items():
            t0 = representative(a, b, u, v)
            assert t0 == min(fiber)
            for t in fiber:
                assert t % t0 == 0
                assert gcd(t, a) == (t//t0)*gcd(t0, a)
                assert gcd(t, b) == (t//t0)*gcd(t0, b)
                reduced += int(t != t0)
            n_pairs += 1
            diagonals += int(u == v)

    stats = Counter()
    for a, b in combinations_with_replacement(range(1, 10), 2):
        for n in range(1, 8):
            size = n*(a+b)
            if size > 24:
                break
            oracle = defaultdict(list)
            for x in partitions(n):
                oracle[explicit_product(a, b, x)].append(x)
            for target in partitions(size):
                solution = solve(a, b, target)
                expected = bool(oracle.get(target))
                assert (solution is not None) == expected, (a,b,target,solution)
                stats["targets"] += 1
                stats["yes" if expected else "no"] += 1
                if gcd(a,b) > 1:
                    stats["noncoprime"] += 1
                if gcd(a,b) > 1 and b % a:
                    stats["noncoprime_incomparable"] += 1
                if expected:
                    assert explicit_product(a, b, solution) == target
                    assert solution in oracle[target]

    target = (4,4,6,6)
    full_fiber = [x for x in partitions(2) if explicit_product(4,6,x) == target]
    assert full_fiber == [(1,1),(2,)]
    assert solve(4,6,target) == (1,1)
    print(json.dumps({"status": "PASS", "period_test": {"divisor_lengths": [1,16],
                     "t": [1,120], "ordered_output_pairs_realized": n_pairs,
                     "diagonal_pairs": diagonals, "individual_periods_shortened": reduced},
                     "exhaustive_targets": dict(stats),
                     "search_not_count_example": {"A":[4,6],"B":target,
                         "all_quotients":full_fiber,"canonical_quotient":solve(4,6,target)}},
                     indent=2))


if __name__ == "__main__":
    main()
