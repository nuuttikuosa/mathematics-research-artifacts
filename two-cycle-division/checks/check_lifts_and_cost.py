"""Check the lifts of a representative vector in Proposition C.3.

The small fibre is obtained from all divisors of the observed periods,
using only the direct product formula. The large case uses binary lists.
"""
from collections import Counter, defaultdict
from math import gcd, isqrt, lcm


def divisors(n):
    small = [i for i in range(1, isqrt(n) + 1) if n % i == 0]
    return sorted(set(small) | {n // i for i in small})


def signature(g, t):
    return lcm(2 * g, t), lcm(3 * g, t)


def least_length(a, u):
    ratio = u // a
    return ratio * gcd(a, pow(ratio, a.bit_length(), a))


def least_pair(g, pair):
    u, v = pair
    return lcm(least_length(2 * g, u), least_length(3 * g, v))


def image(g, x):
    target = Counter()
    for t, count in x.items():
        for a in (2 * g, 3 * g):
            target[lcm(a, t)] += count * gcd(a, t)
    return +target


def intervals():
    total = shared = proper = 0
    for g in range(1, 61):
        for pair in {signature(g, t) for t in range(1, 301)}:
            m = gcd(*pair)
            raw = [t for t in divisors(m) if signature(g, t) == pair]
            s = least_pair(g, pair)
            H = g // gcd(g, s)
            assert s == min(raw) and m == s * H
            assert raw == [s * e for e in divisors(H)]
            total += 1
            shared += pair[0] == pair[1]
            proper += 1 < H < g
    print(f"Complete ordered-pair intervals: {total}; coincident outputs: {shared}; 1 < H < g: {proper}.")


def small_fibre():
    # Both kinds of alternative occur in this one example. The pair (30,30)
    # also contributes twice to the same observed period.
    g = 5
    target = Counter({10: 30, 15: 30, 30: 25})
    raw_lengths = sorted({t for u in target for t in divisors(u)})
    columns = [(t, image(g, {t: 1})) for t in raw_lengths]
    groups = defaultdict(list)

    def enumerate_raw(i, remaining, chosen):
        if i == len(columns):
            if not any(remaining.values()):
                z = Counter()
                for t, n in chosen:
                    s = least_pair(g, signature(g, t))
                    z[s] += n * (t // s)
                groups[tuple(sorted(z.items()))].append(tuple(chosen))
            return
        t, col = columns[i]
        limit = min(remaining.get(u, 0) // n for u, n in col.items())
        residual = remaining.copy()
        for n in range(limit + 1):
            enumerate_raw(i + 1, residual, chosen + ([(t, n)] if n else []))
            residual.subtract(col)

    enumerate_raw(0, target, [])
    # Independent restricted-partition count from all allowable extensions.
    for z, solutions in groups.items():
        expected = 1
        for s, amount in z:
            H = g // gcd(g, s)
            counts = [1] + [0] * amount
            for e in divisors(H):
                for n in range(e, amount + 1):
                    counts[n] += counts[n - e]
            expected *= counts[-1]
        assert len(solutions) == expected
        assert all(image(g, dict(x)) == target for x in solutions)
    key = ((1, 30), (6, 5))
    assert key in groups and ((2, 15), (3, 10)) in groups
    assert signature(g, 6) == signature(g, 30) == (30, 30)
    assert ((1, 30), (30, 1)) in groups[key]
    print(f"g=5 target (30 C10, 30 C15, 25 C30): {len(groups)} representative vectors, "
          f"{sum(map(len, groups.values()))} raw classes; {len(groups[key])} lifts of the selected vector.")


def large_case():
    p = 65537
    assert all(p % d for d in range(2, isqrt(p) + 1))
    g = 2**126 * 3**31 * 5**11 * p**5
    s = 2**128 * 3**33 * 5**12
    u, v = signature(g, s)
    assert u == v and least_pair(g, (u, v)) == s
    H = gcd(u, v) // s
    assert H == p**5 == g // gcd(g, s)
    # There is only one observed output length, hence only one possible
    # ordered pair. Its smallest nonunit lifting factor is p, not 2.
    for z in (p - 1, p, 2**1024 + 1):
        target = image(g, {s: z})
        assert tuple(target) == (u,)
        if z >= p:
            alternate = +Counter({s: z - p, s * p: 1})
            assert image(g, alternate) == target and alternate != {s: z}
        else:
            assert all(e > z for e in (p**i for i in range(1, 6)))
    print(f"Large single-pair target: {g.bit_length()}-bit g; first change at z={p}.")
    print(f"With a 1025-bit multiplicity the second witness uses {len(alternate)} distinct periods; direct product verified.")


if __name__ == "__main__":
    intervals()
    small_fibre()
    large_case()
    print("PASS")
