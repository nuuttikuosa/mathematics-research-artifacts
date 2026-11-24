"""Check prime scaling in Corollary 8.1 by multiplying successor tables.

Run with Python 3.9+ using only the standard library. For each tested
quotient size every permutation isomorphism class is enumerated once.
"""

from collections import Counter, defaultdict
from functools import lru_cache
from itertools import product
from math import gcd, isqrt, lcm


def partitions(n, least=1):
    if n == 0:
        yield ()
    for t in range(least, n + 1):
        for rest in partitions(n - t, t):
            yield (t,) + rest


def successors(lengths):
    table = []
    for length in lengths:
        start = len(table)
        table.extend(range(start + 1, start + length))
        table.append(start)
    return table


def product_type(left, right):
    x, y = successors(left), successors(right)
    table = [u * len(y) + v for u in x for v in y]
    seen = bytearray(len(table))
    periods = []
    for start in range(len(table)):
        if seen[start]:
            continue
        state, length = start, 0
        while not seen[state]:
            seen[state] = 1
            state = table[state]
            length += 1
        assert state == start
        periods.append(length)
    return tuple(sorted(periods))


def scale(lengths, prime):
    return tuple(prime * d for d in lengths)


def first_prime_after(n):
    p = n + 1
    while any(p % d == 0 for d in range(2, isqrt(p) + 1)):
        p += 1
    return p


def product_fibres(divisor, n):
    fibres = defaultdict(set)
    for quotient in partitions(n):
        fibres[product_type(divisor, quotient)].add(quotient)
    return dict(fibres)


def test_scaling():
    positive = negative = supported_negative = multiple = checks = 0
    pairs = set()
    for a in range(1, 7):
        for b in range(a, 7):
            for n in range(1, 7):
                divisor = (a, b)
                size = (a + b) * n
                state_count = a + b + size
                p = first_prime_after(state_count)
                assert state_count < p < 2 * state_count
                base = product_fibres(divisor, n)
                image = product_fibres(scale(divisor, p), n)
                assert image == {scale(B, p): ys for B, ys in base.items()}
                assert gcd(p * a, p * b) == p * gcd(a, b)
                checks += sum(map(len, base.values()))
                pairs.add(divisor)
                positive += len(base)
                multiple += sum(len(ys) > 1 for ys in base.values())

                # Negatives have integral target size. For small targets take
                # every partition. For larger ones split one target cycle.
                targets = set(partitions(size)) if size <= 14 else set()
                for B in base:
                    for i, d in enumerate(B):
                        for e in range(1, d // 2 + 1):
                            targets.add(tuple(sorted(B[:i] + B[i+1:] + (e, d-e))))
                negatives = targets - base.keys()
                available = {d for B in base for d in B}
                for B in negatives:
                    assert scale(B, p) not in image
                    supported_negative += all(d in available for d in B)
                negative += len(negatives)
    print('Two-cycle divisors:', len(pairs), 'including', sum(gcd(*ab)>1 for ab in pairs),
          'with nontrivial gcd and', sum(a == b for a, b in pairs), 'with equal lengths')
    print('Complete quotient-set comparisons on positive inputs:', positive,
          'including', multiple, 'with multiple classes')
    print('Quotient classes checked across positive inputs:', checks)
    print('Integral negative target comparisons:', negative,
          'including', supported_negative, 'using only attainable periods')


@lru_cache(None)
def divisors(n):
    return tuple(sorted({v for i in range(1, isqrt(n) + 1) if n % i == 0
                         for v in (i, n // i)}))


def fixed_axes(a, b, h):
    alpha, beta = h, 1
    while gcd(alpha, b) > 1:
        common = gcd(alpha, b)
        alpha //= common
        beta *= common
    return (h * alpha * a, h * beta * b), (alpha, beta)


@lru_cache(None)
def supported_columns(axes, periods):
    # Any possible y must divide some observed period. Do not require
    # divisibility by h, or by h squared, or assume which axis hits a row.
    available = set(periods)
    candidates = {t for length in periods for t in divisors(length)}
    columns = {}
    for t in sorted(candidates):
        col = Counter()
        for d in axes:
            col[lcm(d, t)] += gcd(d, t)
        if col.keys() <= available:
            columns[t] = col
    return columns


def table_fibres(axes, target):
    candidates = list(supported_columns(axes, tuple(sorted(target))).items())
    answers = set()

    def enumerate_counts(i, remaining, quotient):
        if i == len(candidates):
            if not any(remaining.values()):
                answers.add(quotient)
            return
        period, col = candidates[i]
        limit = min(remaining[v] // coeff for v, coeff in col.items())
        for k in range(limit + 1):
            rest = remaining.copy()
            for v, coeff in col.items():
                rest[v] -= k * coeff
            addition = ((period, k),) if k else ()
            enumerate_counts(i + 1, rest, quotient + addition)

    enumerate_counts(0, target, ())
    return answers


def target_from_margins(a, b, rows, cols, margins, h, split):
    alpha, beta = split
    targets = Counter()
    for u, r in zip(rows, margins[:len(rows)]):
        if r:
            targets[h*h*a*u] += h*alpha*r
    for v, s in zip(cols, margins[len(rows):]):
        if s:
            targets[h*h*b*v] += h*beta*s
    return targets


def test_fixed_gcd():
    hs = (1, 2, 3, 4, 6, 12, 18, 30, 72)
    supports = pairs_checked = 0
    for a, b, h in product(range(2, 10), range(2, 10), hs):
        if gcd(a, b) != 1:
            continue
        rows, cols = divisors(b)[:-1], divisors(a)[:-1]
        axes, (alpha, beta) = fixed_axes(a, b, h)
        assert gcd(*axes) == h
        periods = tuple(sorted([h*h*a*u for u in rows] + [h*h*b*v for v in cols]))
        expected = {h*h*u*v: Counter({h*h*a*u: h*alpha*v, h*h*b*v: h*beta*u})
                    for u in rows for v in cols}
        assert supported_columns(axes, periods) == expected
        supports += 1
        pairs_checked += len(expected)

    yes = no = multi = balanced_no = class_count = 0
    for a, b, u, v in ((4, 9, 3, 2), (9, 4, 2, 3), (6, 5, 1, 3)):
        rows = tuple(sorted({1, u}))
        cols = (1, v)
        margins_set = set(product(range(4), repeat=len(rows) + len(cols)))
        if len(rows) == 2:
            # These two-row margins have several solutions.
            margins_set.add((u*v, v, u*v, u))
        for margins in sorted(margins_set):
            if not any(margins):
                continue
            target = target_from_margins(a, b, rows, cols, margins, 1, (1, 1))
            before = table_fibres((a, b), target)
            balanced = (sum(w*c for w, c in zip(rows, margins)) ==
                        sum(w*c for w, c in zip(cols, margins[len(rows):])))
            for h in hs:
                axes, split = fixed_axes(a, b, h)
                after_target = target_from_margins(a, b, rows, cols, margins, h, split)
                after = table_fibres(axes, after_target)
                expected = {tuple((h*h*t, k) for t, k in q) for q in before}
                assert after == expected
                if before:
                    yes += 1
                    multi += len(before) > 1
                    class_count += len(before)
                else:
                    no += 1
                    balanced_no += balanced
    print('Fixed-gcd complete allowed-column checks:', supports,
          'supports and', pairs_checked, 'columns')
    print('Fixed-gcd whole-fibre checks:', yes, 'positive and', no,
          'negative, including', balanced_no, 'negative with balanced margins')
    print('Fixed-gcd positive fibres:', multi, 'with several classes;',
          class_count, 'classes in total')


if __name__ == '__main__':
    test_scaling()
    divisor = (3, 5)
    target = (6, 10)
    assert product_fibres(divisor, 2)[target] == {(2,)}
    assert scale(target, 2) not in product_fibres(scale(divisor, 2), 2)
    print('Control with prime equal to quotient size: positive became negative')
    test_fixed_gcd()
    print('PASS')
