"""Adversarial test of Corollary 8.1, including powers of primes on both axes.

Run: python -B check_mixed_axes.py (Python 3.9+, standard library).
All calculations are with exact integers. The target supports are small;
the represented permutations may have very large successor tables.
"""

from functools import lru_cache
from itertools import product
from math import gcd, lcm


@lru_cache(None)
def divisors(number):
    result = [1]
    p = 2
    while p * p <= number:
        old = result[:]
        power = 1
        while number % p == 0:
            number //= p
            power *= p
            result.extend(v * power for v in old)
        p = 3 if p == 2 else p + 2
    if number > 1:
        result += [v * number for v in result]
    return result


def possible_cycles(axes, periods):
    # Guess the two actual output periods in order, without any row/column
    # labels or restriction on y. Then y must divide their gcd.
    choices = {}
    for left, right in product(periods, repeat=2):
        for y in divisors(gcd(left, right)):
            if lcm(axes[0], y) == left and lcm(axes[1], y) == right:
                choices[y] = (left, right)
    return choices


def fixed_point_oracle(axes, periods, lengths):
    # For each observed period k, count fixed points of the k-th iterate.
    # These evaluations identify a target on this support: the least
    # differing period would give a differing count. For a product, the
    # fixed-point count is the product of the two fixed-point counts.
    columns = []
    for y in lengths:
        columns.append(tuple(sum(d for d in axes if k % d == 0) *
                             (y if k % y == 0 else 0) for k in periods))
    assert all(any(col) for col in columns)

    @lru_cache(None)
    def solve(index, residual):
        if index == len(columns):
            return ((),) if not any(residual) else ()
        col = columns[index]
        maximum = min(r // c for r, c in zip(residual, col) if c)
        answers = []
        for amount in range(maximum + 1):
            rest = tuple(r - amount * c for r, c in zip(residual, col))
            answers.extend((amount,) + tail for tail in solve(index + 1, rest))
        return tuple(answers)

    def solutions(multiplicities):
        observations = tuple(sum(d * m for d, m in zip(periods, multiplicities) if k % d == 0)
                             for k in periods)
        return solve(0, observations)

    return solutions


def check_axes(a, b, h, alpha, beta, rows, cols):
    assert h == alpha * beta and gcd(alpha * a, beta * b) == 1
    periods = tuple(a * u for u in rows) + tuple(b * v for v in cols)
    scaled_periods = tuple(h * h * d for d in periods)
    scaled_axes = (h * alpha * a, h * beta * b)
    assert gcd(*scaled_axes) == h
    original = possible_cycles((a, b), periods)
    scaled = possible_cycles(scaled_axes, scaled_periods)
    predicted = {h*h*y: (h*h*left, h*h*right) for y, (left, right) in original.items()}
    assert scaled == predicted
    assert original == {u*v: (a*u, b*v) for u, v in product(rows, cols)}
    return periods, scaled_axes, scaled_periods, sorted(original), sorted(scaled)


def check_margins(a, b, h, alpha, beta, rows, cols, bound):
    periods, axes, scaled_periods, original, scaled = check_axes(a, b, h, alpha, beta, rows, cols)
    before = fixed_point_oracle((a, b), periods, original)
    after = fixed_point_oracle(axes, scaled_periods, scaled)
    pos = neg = balanced_neg = multi = 0
    scales = (h * alpha,) * len(rows) + (h * beta,) * len(cols)
    for margins in product(range(bound + 1), repeat=len(periods)):
        if not any(margins):
            continue
        expected = before(margins)
        actual = after(tuple(m * c for m, c in zip(margins, scales)))
        assert actual == expected, (a, b, h, margins)
        if actual:
            pos += 1
            multi += int(len(actual) > 1)
        else:
            neg += 1
            row_sum = sum(u * r for u, r in zip(rows, margins))
            col_sum = sum(v * s for v, s in zip(cols, margins[len(rows):]))
            balanced_neg += int(row_sum == col_sum)
    return pos, neg, balanced_neg, multi


def controls():
    # Incomplete prime-power allocation: gcd(h,b)=3 misses a factor of 3
    # in h=36. From A=B=C4+C3 a wrong split loses the original quotient.
    axes = (36 * 12 * 4, 36 * 3 * 3)
    periods = (36**2 * 4, 36**2 * 3)
    lengths = sorted(possible_cycles(axes, periods))
    assert not fixed_point_oracle(axes, periods, lengths)((36*12, 36*3))
    assert gcd(*axes) == 108

    # If the two supports both contain ab, separate cancellation fails.
    # The same B=5C6 admits C6 in the original unseparated equation.
    h, alpha, beta = 36, 4, 9
    axes = (h*alpha*2, h*beta*3)
    periods = (h*h*6,)
    lengths = sorted(possible_cycles(axes, periods))
    assert lengths == list(periods)
    solve = fixed_point_oracle(axes, periods, lengths)
    assert fixed_point_oracle((2, 3), (6,), (6,))((5,)) == ((1,),)
    assert solve((h*(alpha*2 + beta*3),)) == ((1,),)
    assert not solve((h*alpha*5,))
    assert not solve((h*beta*5,))

    # An overlap can also introduce a genuinely spurious solution.
    # The improper row u=b=9 and column v=a=4 both correspond to length
    # 36; declaring R=13,S=1 gives an impossible old count of 14.
    h, alpha, beta = 72, 8, 9
    old_lengths = sorted(possible_cycles((4, 9), (36,)))
    assert old_lengths == [36]
    assert not fixed_point_oracle((4, 9), (36,), old_lengths)((14,))
    axes = (h*alpha*4, h*beta*9)
    periods = (h*h*36,)
    lengths = sorted(possible_cycles(axes, periods))
    assert lengths == [186624]
    assert fixed_point_oracle(axes, periods, lengths)((h*alpha*13 + h*beta,)) == ((1,),)

    # Distinct labels alone do not separate the two projections. An
    # improper column v=a allows the first projection to hit length ab.
    old_periods = (2, 3, 6)
    old_cycles = sorted(possible_cycles((2, 3), old_periods))
    assert (2 + 28*3 + 4*6) // 5 == 22
    assert not fixed_point_oracle((2, 3), old_periods, old_cycles)((1, 28, 4))
    h, alpha, beta = 36, 4, 9
    axes = (h*alpha*2, h*beta*3)
    periods = tuple(h*h*d for d in old_periods)
    scaled_cycles = sorted(possible_cycles(axes, periods))
    solutions = fixed_point_oracle(axes, periods, scaled_cycles)((h*alpha, 28*h*beta, 4*h*beta))
    witnesses = [{y: k for y, k in zip(scaled_cycles, qs) if k} for qs in solutions]
    assert {1296: 1, 3888: 9} in witnesses


def main():
    # Hand-specified allocations, independent of any splitting routine. Some
    # exponents exceed, and others fall below, their divisor-axis exponents.
    # In the last group 5 divides neither original axis.
    panels = (
        (4, 9, 72, 8, 9, (1, 3), (1, 2), 6),
        (9, 4, 72, 9, 8, (1, 2), (1, 3), 6),
        (16, 27, 108, 4, 27, (1, 3), (1, 2), 6),
        (27, 16, 108, 27, 4, (1, 2), (1, 3), 6),
        (8, 9, 2160, 80, 27, (1, 3), (1, 2), 6),
    )
    results = [check_margins(*args) for args in panels]
    print('Mixed-axis panels:', len(panels))
    print('Positive, negative, balanced negative, multiple:', *map(sum, zip(*results)))
    assert min(map(sum, zip(*results))) > 0
    controls()
    print('Detected the incomplete-prime-power and overlapping-support controls.')
    print('PASS')


if __name__ == '__main__':
    main()
