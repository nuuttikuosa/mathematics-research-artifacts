"""Compare every ordered anchored factor pair with an independent direct-product oracle."""
import sys
from collections import Counter, defaultdict
from functools import cache
from itertools import combinations_with_replacement, product

import sympy
from tree_arithmetic import Trees
from weighted_split import local_pairs_twogrades


# The oracle works only with nested tuples of children. It uses no tree
# division, symbolic polynomials, solver features, or solver memoization.
def star(n):
    return ((),) * n


def path(height):
    t = ()
    for _ in range(height):
        t = (t,)
    return t


@cache
def direct(s, t):
    return tuple(sorted(direct(a, b) for a in s for b in t))


def multiply(left, right):
    return Counter(direct(x, y) for x in left for y in right)


def signature(multiset):
    return tuple(sorted(multiset.elements()))


def possible_axis(target, opposite):
    # If the anchor is a path, multiplication by it is the identity at
    # exactly this height. In the remaining tests all trees are positive stars.
    def is_path(t):
        return len(t) == 0 or (len(t) == 1 and is_path(t[0]))
    if is_path(opposite):
        return sorted(target)
    assert all(all(child == () for child in root) for root in [opposite, *target])
    degree = len(opposite)
    return sorted(star(len(t)//degree) for t in target if len(t) % degree == 0)


def possible_rows(target, anchor, opposite, count):
    eligible = possible_axis(target, opposite)
    if anchor not in eligible:
        return []
    out = []
    for rest in combinations_with_replacement(eligible, count-1):
        row = (anchor, *rest)
        # The opposite anchor is present, so this partial product must fit.
        if Counter(direct(x, opposite) for x in row) <= target:
            out.append(tuple(sorted(row)))
    return out


def brute(target, left_anchor, right_anchor):
    out = set()
    mass = target.total()
    for m in range(1, mass+1):
        if mass % m:
            continue
        left = possible_rows(target, left_anchor, right_anchor, m)
        right = possible_rows(target, right_anchor, left_anchor, mass//m)
        for first, second in product(left, right):
            if multiply(first, second) == target:
                out.add((first, second))
    return out


def solver_pairs(target, left_anchor, right_anchor):
    trees = Trees()
    @cache
    def encode(t):
        return trees.node(encode(x) for x in t)
    @cache
    def decode(i):
        return tuple(sorted(decode(child) for child in trees.children[i]))
    stats = Counter()
    encoded = Counter({encode(t): n for t, n in sorted(target.items())})
    pairs = local_pairs_twogrades(trees, encoded, encode(left_anchor), encode(right_anchor), stats)
    def side(axis):
        return tuple(sorted(decode(i) for i in axis.elements()))
    answer = {(side(left), side(right)) for left, right in pairs}
    assert len(answer) == len(pairs)
    return answer, stats


def cases():
    # Complete enumeration on this explicitly specified alphabet/test range.
    for size in (2, 3, 4):
        for terms in combinations_with_replacement([1, 2, 4, 6], size):
            target = Counter(star(i) for i in terms)
            for a, b in [(1,1), (1,2), (2,2), (2,3)]:
                yield "stars", target, star(a), star(b)
    # Two distinct trees with equal widths and different sibling groups.
    a = tuple(sorted([star(2), ()]))
    b = (star(1), star(1))
    for height in range(2, 6):
        unit = path(height)
        for xs in combinations_with_replacement([unit, a, b], 2):
            for ys in combinations_with_replacement([unit, a, b], 2):
                yield f"height-{height}", multiply((unit, *xs), (unit, *ys)), unit, unit
        a, b = tuple(sorted(((), a))), tuple(sorted(((), b)))
    # Unit-anchored geometric stars encode x^i*y^j as degree 2^i*3^j.
    # This positive polynomial specializes at a=1 to a repeated root.
    target = multiply([star(1), star(2)], [star(1), star(12), star(32), star(32)])
    yield "twist-challenge", target, star(1), star(1)
    # The linear equations admit a candidate with mixed signs in its partner.
    target = Counter({star(2**i): n for i, n in enumerate([1,2,1,1,2,1])})
    yield "signed-cofactor", target, star(1), star(1)
    # Three independent small factors giving a rank-three positive split.
    target = multiply([star(1), star(2)], [star(i) for i in (1,3,5,15)])
    yield "three-factors", target, star(1), star(1)


def main():
    print(f"Python {sys.version_info.major}.{sys.version_info.minor}, SymPy {sympy.__version__}")
    summary = defaultdict(Counter)
    distinct = set()
    for family, target, a, b in cases():
        distinct.add((signature(target), a, b))
        actual, stats = solver_pairs(target, a, b)
        expected = brute(target, a, b)
        if actual != expected:
            raise AssertionError((family, len(target), len(expected), len(actual)))
        counts = summary[family]
        counts["cases"] += 1
        counts["pairs"] += len(actual)
        counts["negative"] += not actual
        counts["general"] += stats["twisted_calls"]
        counts["general_negative"] += bool(stats["twisted_calls"] and not actual)
        counts["nonidentity"] += stats["nonidentity_twists"]
        counts["max_rank"] = max(counts["max_rank"],stats["max_projected_rank"])
        if family == "twist-challenge":
            assert stats["twisted_calls"] == 1
        if family == "three-factors":
            assert stats["max_projected_rank"] >= 3
    keys = "cases pairs negative general general_negative nonidentity max_rank".split()
    print("family", *keys)
    for family, counts in summary.items():
        print(family, *(counts[key] for key in keys))
    total = {key:sum(row[key] for row in summary.values()) for key in keys[:-1]}
    print("total", *total.values(),max(row["max_rank"] for row in summary.values()))
    assert total["general_negative"] and total["nonidentity"]
    print("distinct (target, ordered anchors) among trials:", len(distinct))
    print("mismatches: 0")


if __name__ == "__main__":
    main()
