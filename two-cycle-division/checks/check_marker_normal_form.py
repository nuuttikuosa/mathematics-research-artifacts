"""Shared reference routines for normal-form and transient-target checks.

The construction takes a pair of cycle lengths and a target cycle list.
The raw quotient oracle uses only the cycle product formula and target counts.
"""

from collections import Counter, defaultdict
from itertools import count
from math import gcd, isqrt, lcm



def freeze(multiplicities):
    return tuple(sorted((length, n) for length, n in multiplicities.items() if n))


def multiply(axes, solution):
    target = Counter()
    for t, copies in solution:
        for axis in axes:
            target[lcm(axis, t)] += gcd(axis, t) * copies
    return freeze(target)


def solutions(axes, target):
    # If r and c are the two projections of t, gcd(r,c)/t divides
    # gcd(axes). This bounds the candidates without scanning the new periods.
    g = gcd(*axes)
    ds = [d for d in range(1, g + 1) if g % d == 0]
    present = dict(target)
    candidates = set()
    for r in present:
        for c in present:
            common = gcd(r, c)
            for d in ds:
                if common % d == 0:
                    t = common // d
                    if all(lcm(axis, t) in present for axis in axes):
                        candidates.add(t)
    columns = [(t, multiply(axes, [(t, 1)])) for t in sorted(candidates)]
    found = set()

    def visit(i, residual, profile):
        if i == len(columns):
            if not any(residual.values()):
                found.add(tuple(profile))
            return
        t, col = columns[i]
        bound = min(residual[r] // value for r, value in col)
        for x in range(bound + 1):
            after = residual.copy()
            for r, coefficient in col:
                after[r] -= x * coefficient
            visit(i + 1, after, profile + ([(t, x)] if x else []))

    visit(0, present, [])
    return found


def small_primes():
    for p in count(2):
        if all(p % j for j in range(2, isqrt(p) + 1)):
            yield p


def support_part(value, w):
    if value == 1:
        return 1
    return gcd(value, pow(w, value.bit_length(), value))


def marked_target(axes, target, use_markers):
    a, b = axes
    g, D = gcd(a, b), lcm(a, b)
    if use_markers:
        forbidden = [a, b] + list(dict(target))
        available = (p for p in small_primes() if all(x % p for x in forbidden))
        m, n = next(available), next(available)
        U, V = (a // g) * m, (b // g) * n
        p = pow(U, -1, V)
        q = (U * p - 1) // V
        P, Q = m * p, n * q
    else:
        U, V = a // g, b // g
        p = pow(U, -1, V) if V > 1 else 1
        q = (U * p - 1) // V
        P, Q = p + 2 * V, q + 2 * U
    new_axes = a * P, b * Q
    assert new_axes[0] - new_axes[1] == g == gcd(*new_axes)
    assert lcm(*new_axes) == D * P * Q
    failed = (new_axes, ((1, 1),), {}, (P, Q), (False, False))
    blocks = defaultdict(dict)
    for ell, multiplicity in target:
        blocks[lcm(D, ell)][ell] = multiplicity
    info, output = {}, Counter()
    seen_shared = seen_new_shared = False
    for R, observed in sorted(blocks.items()):
        rows, cols = {}, {}
        for ell, qty in observed.items():
            if ell == R:
                continue
            if ell % a == 0:
                rows[ell] = qty
            elif ell % b == 0:
                cols[ell] = qty
            else:
                return failed
        W = sum(ell * mult for ell, mult in observed.items())
        number, rem = divmod(W, a + b)
        if rem:
            return failed
        numerator = a * number - sum(ell * mult for ell, mult in rows.items())
        share, rem = divmod(numerator, R)
        if rem or not 0 <= share <= observed.get(R, 0):
            return failed
        rows[R] = share
        cols[R] = observed.get(R, 0) - share
        seen_shared |= rows[R] > 0 and cols[R] > 0
        PF, QF = support_part(P, R // D), support_part(Q, R // D)
        info[R] = PF, QF
        for ell, mult in rows.items():
            if mult:
                new = P * QF * ell
                if use_markers:
                    assert new % m == 0 and new % n and new not in output
                assert lcm(lcm(*new_axes), new) == P * Q * R
                assert new not in output
                output[new] += PF * mult
        for ell, mult in cols.items():
            if mult:
                new = Q * PF * ell
                if use_markers:
                    assert new % n == 0 and new % m and new not in output
                assert lcm(lcm(*new_axes), new) == P * Q * R
                if new in output:
                    assert new == P * Q * R
                    seen_new_shared = True
                output[new] += QF * mult
    return new_axes, freeze(output), info, (P, Q), (seen_shared, seen_new_shared)
