"""Compare quotient classes before and after imposing consecutive cofactors.

Run with Python 3.9+ and the standard library:
    python -B checks/check_consecutive_divisors.py
"""

from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from math import gcd, isqrt, lcm, prod


def axes(a, b, h):
    alpha = h
    while gcd(alpha, b) != 1:
        alpha //= gcd(alpha, b)
    beta = h // alpha
    u, v = alpha * a, beta * b
    assert gcd(u, v) == 1
    p = pow(u, -1, v)
    q = (u*p-1) // v
    assert 1 <= p < v and 1 <= q < u and u*p == v*q+1
    return (h*u*p, h*v*q), (p*a, q*b), (alpha, beta)


def target(a, b, rows, cols, margins, h, alpha, beta):
    result = Counter()
    for weight, multiplicity in zip(rows, margins[:len(rows)]):
        if multiplicity:
            result[h*h*a*weight] += h*alpha*multiplicity
    for weight, multiplicity in zip(cols, margins[len(rows):]):
        if multiplicity:
            result[h*h*b*weight] += h*beta*multiplicity
    return result


@lru_cache(None)
def divisors(d):
    small = [s for s in range(1, isqrt(d)+1) if d % s == 0]
    return tuple(small + [d//s for s in reversed(small) if d//s != s])


@lru_cache(None)
def columns(divisor, support):
    # The oracle sees only two divisor lengths and the observed target.
    # It tries every ordered output pair, including pairs of equal lengths.
    a, b = divisor
    allowed = {}
    for r, c in product(support, repeat=2):
        if r % a or c % b:
            continue
        # Any cycle giving both periods necessarily divides their gcd.
        for t in divisors(gcd(r, c)):
            if lcm(a, t) == r and lcm(b, t) == c:
                contribution = Counter()
                contribution[r] += gcd(a, t)
                contribution[c] += gcd(b, t)
                allowed[t] = contribution
    return tuple(sorted(allowed.items()))


def quotients(divisor, result):
    available = columns(divisor, tuple(sorted(result)))
    answers = set()

    def visit(index, residual, chosen):
        if index == len(available):
            if not any(residual.values()):
                answers.add(tuple(chosen))
            return
        t, vector = available[index]
        limit = min(residual[d] // k for d, k in vector.items())
        for k in range(limit+1):
            tail = {d: num - k*vector[d] for d, num in residual.items()}
            visit(index+1, tail, chosen + ([(t, k)] if k else []))

    visit(0, result, [])
    return answers


def compare_tables():
    statistics = Counter()
    for a, b, row, col in [(4, 9, 3, 2), (9, 4, 2, 3), (15, 7, 1, 5)]:
        rows, cols = tuple(sorted({1, row})), (1, col)
        cases = set(product(range(4), repeat=len(rows)+len(cols)))
        if len(rows) == 2:
            cases.add((row*col, col, row*col, row))
        for margins in cases:
            if not any(margins):
                continue
            before = target(a, b, rows, cols, margins, 1, 1, 1)
            old = quotients((a, b), before)
            for h in (1, 6, 12, 72):
                new_axes, (a0, b0), (alpha, beta) = axes(a, b, h)
                assert gcd(*new_axes) == h and new_axes[0]-new_axes[1] == h
                after = target(a0, b0, rows, cols, margins, h, alpha, beta)
                old_size = a + b + sum(t * k for t, k in before.items())
                new_size = sum(new_axes) + sum(t * k for t, k in after.items())
                assert new_size <= h**4 * old_size**2
                expected = {tuple((t*h*h, k) for t, k in sol) for sol in old}
                assert quotients(new_axes, after) == expected
                statistics['yes' if old else 'no'] += 1
                statistics['multiple'] += len(old) > 1
                statistics['classes'] += len(old)
                if not old:
                    split = len(rows)
                    statistics['balanced_no'] += (
                        sum(w*k for w, k in zip(rows, margins[:split]))
                        == sum(w*k for w, k in zip(cols, margins[split:]))
                    )
    return statistics


def chosen_quotient(weights, chosen, h):
    return tuple(sorted((h*h*(w if i in chosen else 2*w), 2 if i in chosen else 1)
                        for i, w in enumerate(weights)))


def check_binary_case(weights, total, h, choices):
    b = 3*prod(weights)
    new_axes, (a0, b0), (alpha, beta) = axes(4, b, h)
    assert 2 < a0 and gcd(*new_axes) == h and new_axes[0]-new_axes[1] == h
    margins = [2]*len(weights) + [2*total, sum(weights)-total]
    assert min(margins) >= 0
    output = target(a0, b0, weights, (1, 2), margins, h, alpha, beta)
    answers = quotients(new_axes, output)
    assert answers == {chosen_quotient(weights, chosen, h) for chosen in choices}
    return len(answers)


def compare_subsets():
    stats = Counter()
    for weights in [(1, 3, 5, 7), (3, 5, 7, 9)]:
        subsets = [s for k in range(len(weights)+1) for s in combinations(range(len(weights)), k)]
        for total in range(sum(weights)+1):
            choices = [s for s in subsets if sum(weights[i] for i in s) == total]
            for h in (1, 6, 12, 72):
                count = check_binary_case(weights, total, h, choices)
                stats['yes' if count else 'no'] += 1
                stats['multiple'] += count > 1
                stats['classes'] += count
    return stats


def compare_exact_covers():
    families = [
        [(0, 1, 3), (0, 2, 4), (1, 4, 5), (0, 4, 5)],
        [(0, 1, 2), (3, 4, 5), (0, 2, 3), (0, 1, 5)],
        [(0, 1, 2), (3, 4, 5), (0, 1, 3), (2, 4, 5)],
    ]
    for expected_count, triples in enumerate(families):
        base = len(triples)+2
        weights = [1+2*base*sum(base**j for j in triple) for triple in triples]
        total = 2+2*base*sum(base**j for j in range(6))
        covers = [s for s in combinations(range(len(triples)), 2)
                  if sorted(j for i in s for j in triples[i]) == list(range(6))]
        assert len(covers) == expected_count
        for h in (1, 6, 12, 72):
            assert check_binary_case(weights, total, h, covers) == expected_count


def main():
    for name, stats in [('general tables', compare_tables()), ('odd subset weights', compare_subsets())]:
        print(f'{name}: ' + ', '.join(f'{k}={stats[k]}' for k in sorted(stats)))
    compare_exact_covers()
    print('genuine exact-cover reductions: 3 families, with 0, 1, and 2 covers; h=1,6,12,72')
    print('PASS')


if __name__ == '__main__':
    main()
