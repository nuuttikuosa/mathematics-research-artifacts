"""Seeded checks for Corollary 8.2 and Proposition C.3.

The reference for the two-by-two solver enumerates quotients from the product
equation, then scans small periods to normalise each quotient. It uses neither
the forced margins nor the CRT solver. The factor-recovery check uses a
simulated decision oracle and records an exact per-prime query budget.
"""

from collections import Counter, defaultdict
from math import gcd, lcm
import json
from pathlib import Path
from random import Random
import sys

# The reproduction runner invokes this script with Python's isolated mode.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_crt_blocks import all_vectors, blocks, freeze
from check_exact_covers import all_divisors, all_quotients
from check_factoring_queries import trial_factors
from check_marker_normal_form import support_part


def word_lengths():
    rng = Random(1)
    shared_prime_pairs = 0
    for _ in range(20_000):
        a, b = rng.randint(1, 60), rng.randint(1, 60)
        u, v = a // gcd(a, b), b // gcd(a, b)
        p0 = pow(u, -1, v) if v > 1 else 1
        q0 = (u * p0 - 1) // v
        p, q = p0 + 2 * v, q0 + 2 * u
        t = rng.randint(1, 400)
        base = lcm(a, b)
        label = lcm(base, t)
        scale = support_part(p, label // base) * support_part(q, label // base)
        y = scale * t
        assert gcd(y, label) == t, (a, b, t, y, label)
        for e in all_divisors(label):
            assert (y % e == 0) == (t % e == 0), (a, b, t, e)
            shared_prime_pairs += gcd(e, scale) > 1
    return {"triples": 20_000, "trial_divisor_pairs_sharing_a_prime": shared_prime_pairs}


def target_sample(rng, g):
    # Preserve axis order and first occurrence of image lengths. The random
    # perturbation chooses a position in that order.
    quotient = Counter()
    divisors = sorted(all_divisors(6 * g))
    for _ in range(rng.randint(1, 4)):
        length = rng.choice(divisors)
        quotient[length] += rng.randint(1, 4)
    target = Counter()
    for axis in (2 * g, 3 * g):
        for length, copies in quotient.items():
            target[lcm(axis, length)] += gcd(axis, length) * copies
    if rng.random() < 0.4:
        length = rng.choice(list(target))
        target[length] -= 1
        target = +target
        target[rng.choice([length, 2 * length, 3 * length])] += 1
    return target


def reference_vectors(g, target, size_bound):
    axes = (2 * g, 3 * g)
    solutions = all_quotients(axes, target)
    minima = {}
    for t in range(1, size_bound + 1):
        minima.setdefault(tuple(lcm(axis, t) for axis in axes), t)
    vectors = set()
    for solution in solutions:
        vector = Counter()
        for t, copies in solution:
            s = minima[tuple(lcm(axis, t) for axis in axes)]
            count, remainder = divmod(t, s)
            assert remainder == 0
            vector[s] += copies * count
        vectors.add(freeze(vector))
    return vectors


def margin_shapes(g, target):
    """Count the blocks at the forced-margin stage, before CRT feasibility."""
    a, b = 2 * g, 3 * g
    groups = defaultdict(Counter)
    for length, copies in target.items():
        if length % a and length % b:
            return Counter()
        groups[lcm(a, b, length)][length] = copies
    counts = Counter()
    for label, group in groups.items():
        total = sum(t * k for t, k in group.items())
        size, remainder = divmod(total, a + b)
        if remainder:
            return Counter()
        rows = [t for t in group if t != label and t % a == 0]
        columns = [t for t in group if t != label and t % b == 0]
        row_share, remainder = divmod(a * size - sum(t * group[t] for t in rows), label)
        column_share = group[label] - row_share
        if remainder or min(row_share, column_share) < 0:
            return Counter()
        counts["blocks_with_shared_label"] += min(row_share, column_share) > 0
        counts["padded_blocks"] += (len(rows) + (row_share > 0) < 2
                                    or len(columns) + (column_share > 0) < 2)
    return counts


def representative_vectors():
    rng = Random(2)
    totals = Counter()
    for _ in range(1_500):
        g = rng.randint(1, 30)
        target = target_sample(rng, g)
        size_bound = sum(t * k for t, k in target.items()) // (5 * g)
        if size_bound > 60:
            continue
        reference = reference_vectors(g, target, size_bound)
        computed = all_vectors(blocks(g, freeze(target)))
        assert reference == computed, (g, target, reference, computed)
        totals["targets"] += 1
        totals["feasible"] += bool(reference)
        totals["multiple_vectors"] += len(reference) > 1
        totals.update(margin_shapes(g, target))
    return dict(totals)


def recover(g, ask):
    """Find the first positive query by bracketing and inclusive binary search."""
    current = g
    factors = []
    calls = 0
    while current > 1:
        previous, high = 1, 2
        used = 1
        while not ask(current, high):
            previous, high = high, min(2 * high, current)
            used += 1
        low = previous + 1
        while low < high:
            mid = (low + high) // 2
            used += 1
            if ask(current, mid):
                high = mid
            else:
                low = mid + 1
        assert used <= 2 * (low - 1).bit_length()
        calls += used
        factors.append(low)
        current, remainder = divmod(current, low)
        assert remainder == 0
    return factors, calls


def oracle_queries():
    special = [2**20, 3**13, 999_983, 2**10 * 999_983, 7**7,
               2 * 3 * 5 * 7 * 11 * 13 * 17 * 19, 1_000_003 * 1_000_033]
    inputs = list(range(2, 3_001)) + special
    reports = {}
    for g in inputs:
        known = trial_factors(g)

        def ask(current, threshold):
            smallest = next(p for p in known if current % p == 0)
            return threshold >= smallest

        recovered, calls = recover(g, ask)
        assert recovered == known
        if g in special:
            reports[str(g)] = calls
    return {"inputs": len(inputs), "queries_on_selected_inputs": reports}


if __name__ == "__main__":
    print(json.dumps({"word_lengths": word_lengths(),
                      "representative_vectors": representative_vectors(),
                      "simulated_oracle": oracle_queries()}, indent=2, sort_keys=True))
