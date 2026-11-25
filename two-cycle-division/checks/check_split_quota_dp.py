"""Exact split-quota DP versus classified literal Cartesian successor maps.

Run with Python 3.9+ using only the standard library:
    python -B check_split_quota_dp.py
"""
from collections import Counter, defaultdict
from functools import lru_cache
from math import gcd, lcm
from random import Random


def periods(parts):
    return tuple(sorted(Counter(parts).items()))


def states(profile):
    return sum(d * k for d, k in profile)


def parts(n, floor=1):
    if n == 0:
        yield ()
    for d in range(floor, n + 1):
        for tail in parts(n - d, d):
            yield (d,) + tail


def cycle_table(profile):
    nxt = []
    for d, m in profile:
        for _ in range(m):
            base = len(nxt)
            nxt.extend(range(base + 1, base + d))
            nxt.append(base)
    return nxt


@lru_cache(None)
def literal_output(a, b, x):
    f = cycle_table(periods((a, b)))
    h = cycle_table(x)
    order = len(h)
    nxt = [fx * order + hy for fx in f for hy in h]
    seen = bytearray(len(nxt))
    lengths = []
    for start, visited in enumerate(seen):
        if not visited:
            v, size = start, 0
            while not seen[v]:
                seen[v] = 1
                v, size = nxt[v], size + 1
            lengths.append(size)
    return periods(lengths)


def compute_block(a, b, R, obs, shuffled, stats):
    size, rem = divmod(states(obs.items()), a + b)
    if rem or any(r % a and r % b for r in obs):
        return 0, ()
    exclusive_rows = {r: k for r, k in obs.items() if r != R and r % a == 0}
    quota, rem = divmod(a * size - states(exclusive_rows.items()), R)
    if rem or not 0 <= quota <= obs.get(R, 0):
        return 0, ()
    row = dict(exclusive_rows)
    col = {c: k for c, k in obs.items() if c != R and c % b == 0}
    if quota:
        row[R] = quota
    if obs.get(R, 0) > quota:
        col[R] = obs[R] - quota
    if R in obs:
        collision = "both" if R in row and R in col else "row only" if R in row else "column only"
        stats["positive quota at shared length: " + collision] += 1
    else:
        stats["unobserved block label"] += 1
    rows = sorted(row)
    idx = {r: i for i, r in enumerate(rows)}
    limits = tuple(row[r] for r in rows)
    columns = list(sorted(col))
    if shuffled:
        shuffled.shuffle(columns)
    groups = defaultdict(list)
    for d in range(1, size + 1):
        r, c = lcm(a, d), lcm(b, d)
        if r in row and c in col:
            groups[c].append((d, idx[r], gcd(a, d)))
    dp = {(0,) * len(rows): (1, ())}
    closed = 0
    for j, c in enumerate(columns):
        choices = groups[c]
        if not choices:
            return 0, ()
        if shuffled:
            shuffled.shuffle(choices)
        if R == c and R in row and len(columns) > 1:
            where = "first" if j == 0 else "last" if j == len(columns) - 1 else "middle"
            stats["dual-positive shared column processed " + where] += 1
        for stage, (d, i, inc) in enumerate(choices):
            tail = stage + 1 == len(choices)
            nxt = {}
            for prefix, (count, witness) in dp.items():
                for z in range(1 + (limits[i] - prefix[i]) // inc):
                    vector = list(prefix)
                    vector[i] += z * inc
                    state = tuple(vector)
                    top = b * sum(r * s for r, s in zip(rows, state)) - a * closed
                    running, rem = divmod(top, a * c)
                    assert rem == 0
                    if running > col[c] or (tail and running != col[c]):
                        continue
                    previous = nxt.get(state)
                    if previous:
                        nxt[state] = (previous[0] + count, previous[1])
                    else:
                        nxt[state] = (count, witness + (((d, z),) if z else ()))
            dp = nxt
        closed += c * col[c]
    return dp.get(limits, (0, ()))


def solve(a, b, target, shuffled, stats):
    blocks = defaultdict(dict)
    for d, m in target:
        blocks[lcm(a, b, d)][d] = m
    count, witness = 1, ()
    for R, obs in sorted(blocks.items()):
        factor, profile = compute_block(a, b, R, obs, shuffled, stats)
        if not factor:
            return 0, ()
        count *= factor
        witness += profile
    return count, tuple(sorted(witness))


def main():
    rng = Random(35429)
    total, features = Counter(), Counter()
    nmax = 12
    all_x = [periods(x) for n in range(1, nmax + 1) for x in parts(n)]
    for a, b in [(a,b) for a in range(1,6) for b in range(1,6)] + [(4,7),(7,4),(6,10),(10,6)]:
        actual = defaultdict(set)
        for x in all_x:
            actual[literal_output(a,b,x)].add(x)
        targets = set(actual)
        for original in sorted(actual):
            d = rng.choice(original)[0]
            # Move, delete, or split a single cycle while preserving our
            # exhaustive bound on the number of possible quotient states.
            for replacement in [(), (rng.randint(1,d),), (1,d-1) if d>1 else ()]:
                moved = Counter(dict(original))
                moved[d] -= 1
                moved.update(replacement)
                transformed = tuple(sorted((r,k) for r,k in moved.items() if k))
                if transformed:
                    targets.add(transformed)
        for target in sorted(targets):
            full = actual.get(target, set())
            total["instances"] += 1
            total["yes" if full else "no"] += 1
            if len(full)>1:
                total["more than one"] += 1
            total["quotient classes"] += len(full)
            for permute in (None,rng,rng):
                actual_count, witness = solve(a,b,target,permute,features)
                assert actual_count == len(full), (a,b,target,actual_count,full)
                assert not actual_count or (witness in full and literal_output(a,b,witness)==target)
                total["DP runs"] += 1
    for mark in ("row only","column only","both"):
        assert features["positive quota at shared length: "+mark]
    for mark in ("first","middle","last"):
        assert features["dual-positive shared column processed "+mark]
    print("Source quotient sizes: 1..12; representation: literal successor tables")
    print("Ordered divisor pairs: 29; orders per instance: sorted and two seeded shuffles")
    for key, value in sorted(total.items()):
        print(f"{key}: {value}")
    print("Block coverage over all runs:")
    for key, value in sorted(features.items()):
        print(f"{key}: {value}")
    print("PASS")


if __name__ == "__main__":
    main()
