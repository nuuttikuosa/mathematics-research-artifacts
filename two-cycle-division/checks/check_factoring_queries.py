"""Compare exponential bracketing with an oracle for the quotients of kA.

The small oracle counts partitions of k into divisors of h, capped at two.
The recovery algorithm only sees one Boolean answer per query.
"""
from functools import lru_cache
from math import prod


@lru_cache(maxsize=None)
def partition_counts(h):
    counts = [1] + [0] * h
    for part in range(1, h + 1):
        if h % part == 0:
            for k in range(part, h + 1):
                counts[k] = min(2, counts[k] + counts[k - part])
    return counts


def recover(g, ask_multiple):
    factors = []
    counts_per_factor = []
    h = g
    while h > 1:
        lo, hi = 1, 2
        calls = 1
        while not ask_multiple(h, hi):
            assert hi < h
            lo, hi = hi, min(2 * hi, h)
            calls += 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            calls += 1
            if ask_multiple(h, mid):
                hi = mid
            else:
                lo = mid
        factors.append(hi)
        counts_per_factor.append(calls)
        h //= hi
    return factors, counts_per_factor


def trial_factors(n):
    out = []
    divisor = 2
    while divisor * divisor <= n:
        while n % divisor == 0:
            out.append(divisor)
            n //= divisor
        divisor += 1
    if n > 1:
        out.append(n)
    return out


def check_input(g, oracle, expected):
    factors, calls = recover(g, oracle)
    assert factors == expected and prod(factors) == g
    for p, q in zip(factors, calls):
        assert q <= 2 * (p - 1).bit_length() - 1
    return sum(calls)


def main():
    queries = 0
    for g in range(1, 513):
        expected = trial_factors(g)
        queries += check_input(g, lambda h, k: partition_counts(h)[k] >= 2, expected)
        # Every kA query already has k fixed points as a quotient.
        as_multiple = lambda h, k: not (partition_counts(h)[k] == 1)
        queries += check_input(g, as_multiple, expected)
    print(f"All g=1..512, both Boolean oracles, {queries} queries: PASS")

    # Check calls with hundreds of bits without materialising the small oracle.
    # The simulator is given the known prime factors, the recovery code is not.
    expected = [2] * 31 + [257] * 27 + [65537] * 15
    big = prod(expected)
    def simulated(h, k):
        p = next(p for p in (2, 257, 65537) if h % p == 0)
        return k >= p
    big_calls = check_input(big, simulated, expected)
    print(f"Known-threshold simulation: {big.bit_length()}-bit input, "
          f"{len(expected)} factors with multiplicity, {big_calls} queries: PASS")


if __name__ == "__main__":
    main()
