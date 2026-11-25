"""Check Proposition B.2 on a seeded sample of binary cycle lists.

Both candidate sets are tested by the same exhaustive integer search. The
comparison tests preservation of feasibility when one keeps only the smallest
period for each ordered pair of target lengths. No table of states is needed.

A target has a diagonal candidate if some period has equal output lengths in
its two projections. A solution need not use that period. For each retained
period s define D = lcm(a,b) and R = lcm(D,s). The run reports both the largest
divisor count tau(R) for R < 200000 and the maximum over all observed R.
"""

import json
import random
from collections import Counter, defaultdict
from functools import lru_cache
from math import gcd, lcm


SEED = 20261124
TRIALS = 1500
CUTOFF = 200_000
DIVISOR_LENGTHS = (2, 3, 4, 6, 8, 9, 12, 10, 15, 18, 5, 7)
PRIMES = (2, 3, 5, 7, 11, 13)


def sample_target(random_source):
    a = random_source.choice(DIVISOR_LENGTHS)
    b = random_source.choice(DIVISOR_LENGTHS)
    quotient = Counter()
    for _ in range(random_source.randint(1, 3)):
        t = 1
        for p in random_source.sample(PRIMES, random_source.randint(0, 4)):
            t *= p ** random_source.randint(1, 2)
        quotient[t] += random_source.randint(1, 3)
    target = Counter()
    for t, multiplicity in quotient.items():
        for axis in (a, b):
            target[lcm(axis, t)] += gcd(axis, t) * multiplicity
    if random_source.random() < 0.5:
        old = random_source.choice(list(target))
        target[old] -= 1
        if not target[old]:
            del target[old]
        new = random_source.choice([old, 2 * old, 3 * old, a * old, b * old])
        target[new] += 1
    return a, b, target


def minimum_period(axis, output):
    """Return the least t with lcm(axis,t) = output, or None if impossible."""
    if output % axis:
        return None
    ratio = output // axis
    return ratio * gcd(axis, pow(ratio, axis.bit_length(), axis))


def representatives(a, b, target):
    periods = set()
    for u in target:
        for v in target:
            left = minimum_period(a, u)
            right = minimum_period(b, v)
            if left is None or right is None:
                continue
            s = lcm(left, right)
            if lcm(a, s) == u and lcm(b, s) == v:
                periods.add(s)
    return periods


def full_candidates(a, b, target):
    """Use gcd(u,v) = lcm(gcd(a,b),t), as in the counting theorem."""
    common = gcd(a, b)
    divisors = [q for q in range(1, common + 1) if common % q == 0]
    periods = set()
    for u in target:
        for v in target:
            observed_gcd = gcd(u, v)
            for q in divisors:
                if observed_gcd % q:
                    continue
                t = observed_gcd // q
                if lcm(common, t) == observed_gcd and lcm(a, t) in target and lcm(b, t) in target:
                    periods.add(t)
    return periods


def feasible(a, b, target, periods):
    columns = []
    for t in sorted(periods, reverse=True):
        column = Counter()
        for axis in (a, b):
            column[lcm(axis, t)] += gcd(axis, t)
        columns.append(column)

    def search(index, remaining):
        if not any(remaining.values()):
            return True
        if index == len(columns):
            return False
        column = columns[index]
        maximum = min(remaining[key] // value for key, value in column.items())
        for copies in range(maximum, -1, -1):
            after = remaining.copy()
            for key, value in column.items():
                after[key] -= copies * value
            if search(index + 1, after):
                return True
        return False

    return search(0, target)


@lru_cache(maxsize=None)
def divisor_count(value):
    count = 1
    factor = 2
    while factor * factor <= value:
        exponent = 0
        while value % factor == 0:
            value //= factor
            exponent += 1
        count *= exponent + 1
        factor = 3 if factor == 2 else factor + 2
    if value > 1:
        count *= 2
    return count


def main():
    random_source = random.Random(SEED)
    totals = Counter()
    labels = set()
    maximum_tau_D = 0
    maximum_tau_R = 0
    witness = None
    for trial in range(TRIALS):
        a, b, target = sample_target(random_source)
        if sum(target.values()) > 40:
            continue
        periods = representatives(a, b, target)
        complete = full_candidates(a, b, target)
        positive = feasible(a, b, target, complete)
        if feasible(a, b, target, periods) != positive:
            raise AssertionError(f"Feasibility differs on trial {trial}: {(a, b, target)}")

        totals["samples"] += 1
        totals["positive"] += positive
        totals["non_coprime"] += gcd(a, b) > 1
        diagonal = any(lcm(a, s) == lcm(b, s) for s in periods)
        totals["diagonal_candidate"] += diagonal
        totals["positive_with_diagonal_candidate"] += positive and diagonal

        D = lcm(a, b)
        maximum_tau_D = max(maximum_tau_D, divisor_count(D))
        g = gcd(a, b)
        pair_bound = divisor_count(a // g) * divisor_count(b // g)
        block_sizes = defaultdict(int)
        for s in sorted(periods):
            R = lcm(D, s)
            labels.add(R)
            block_sizes[R] += 1
            if divisor_count(R) > maximum_tau_R:
                maximum_tau_R = divisor_count(R)
                witness = {"a": a, "b": b, "s": s, "D": D, "R": R,
                           "tau_D": divisor_count(D), "reduced_pair_bound": pair_bound,
                           "feasible": positive}
        if any(size > pair_bound for size in block_sizes.values()):
            raise AssertionError(f"Too many representatives in a block on trial {trial}")

    result = {
        "seed": SEED,
        "trials_generated": TRIALS,
        **totals,
        "feasibility_disagreements": 0,
        "distinct_labels": len(labels),
        "distinct_labels_at_or_above_cutoff": sum(R >= CUTOFF for R in labels),
        "cutoff": CUTOFF,
        "max_tau_R_below_cutoff": max(divisor_count(R) for R in labels if R < CUTOFF),
        "max_tau_R": maximum_tau_R,
        "max_tau_R_witness": witness,
        "max_tau_D": maximum_tau_D,
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
