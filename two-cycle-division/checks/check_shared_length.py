"""Check the shared-length case in Note B, Theorem 2.2. Python 3.9+."""

from collections import Counter, defaultdict
from functools import lru_cache
from itertools import product
from math import gcd, lcm, prod


def partitions(size, least=1):
    if size == 0:
        yield ()
    for part in range(least, size + 1):
        for tail in partitions(size - part, part):
            yield (part,) + tail


def successors(parts):
    result = []
    for length in parts:
        start = len(result)
        result.extend(range(start + 1, start + length))
        result.append(start)
    return result


def product_type(divisor, quotient):
    """Independent oracle: follow the actual orbits of a product permutation."""
    left, right = successors(divisor), successors(quotient)
    table = [i * len(right) + j for i in left for j in right]
    seen = set()
    lengths = []
    for start in range(len(table)):
        if start in seen:
            continue
        point = start
        length = 0
        while point not in seen:
            seen.add(point)
            length += 1
            point = table[point]
        lengths.append(length)
    return tuple(sorted(lengths))


@lru_cache(None)
def tau(x):
    return sum(x % t == 0 for t in range(1, x + 1))


@lru_cache(None)
def layout(a, b, n, order):
    common = gcd(a, b)
    if tau(a // common) < tau(b // common):
        a, b = b, a
    blocks = defaultdict(lambda: defaultdict(list))
    lengths = set()
    for t in range(1, n + 1):
        row, column = lcm(a, t), lcm(b, t)
        blocks[lcm(a, b, t)][column].append((t, row, gcd(a, t), gcd(b, t)))
        lengths.update((row, column))
    result = []
    for label, columns in sorted(blocks.items()):
        groups = []
        rows = set()
        for col in sorted(columns, reverse=bool(order & 1)):
            terms = sorted(columns[col], reverse=bool(order & 2))
            rows.update(term[1] for term in terms)
            groups.append((col, terms))
        assert len(rows) <= tau(b // common)
        assert rows.intersection(columns) <= {label}
        result.append((label, tuple(sorted(rows)), groups))
    return result, lengths


def block_count(rows, groups, target, n):
    index = {length: i for i, length in enumerate(rows)}
    counts = {tuple(0 for _ in rows): 1}
    for column, terms in groups:
        temporary = column not in index
        if temporary:
            counts = {state + (0,): count for state, count in counts.items()}
        ci = index.get(column, len(rows))
        for t, row, coeff_r, coeff_c in terms:
            ri = index[row]
            new = defaultdict(int)
            for state, count in counts.items():
                for x in range(n // t + 1):
                    dest = list(state)
                    dest[ri] += x * coeff_r
                    dest[ci] += x * coeff_c
                    if dest[ri] > target[row] or dest[ci] > target[column]:
                        break
                    new[tuple(dest)] += count
            counts = new
        if temporary:
            counts = {s[:-1]: k for s, k in counts.items() if s[-1] == target[column]}
    return counts.get(tuple(target[r] for r in rows), 0)


def count_rows_only(a, b, n, rows, groups, target):
    """Independent state compression: reconstruct each ordinary column from mass."""
    common = gcd(a, b)
    if tau(a // common) < tau(b // common):
        a, b = b, a
    index = {row: i for i, row in enumerate(rows)}
    counts = {(0,) * len(rows): 1}
    closed_mass = 0
    for column, terms in groups:
        ci = index.get(column)
        if ci is not None:
            assert column == groups[-1][0]
        for t, row, coeff_r, coeff_c in terms:
            ri = index[row]
            new = defaultdict(int)
            for state, count in counts.items():
                limit = min(n // t, (target[row] - state[ri]) // coeff_r)
                if ci is None:
                    weighted = sum(r * s for r, s in zip(rows, state))
                    numerator = b * weighted - a * closed_mass
                    part, rem = divmod(numerator, a * column)
                    assert rem == 0 and 0 <= part <= target[column]
                    limit = min(limit, (target[column] - part) // coeff_c)
                else:
                    per_copy = coeff_c + (coeff_r if ri == ci else 0)
                    limit = min(limit, (target[column] - state[ci]) // per_copy)
                for x in range(limit + 1):
                    dest = list(state)
                    dest[ri] += coeff_r * x
                    if ci is not None:
                        dest[ci] += coeff_c * x
                    new[tuple(dest)] += count
            counts = new
        if ci is None:
            closed_mass += column * target[column]
            counts = {s: k for s, k in counts.items()
                      if b * sum(r * v for r, v in zip(rows, s)) == a * closed_mass}
    return counts.get(tuple(target[r] for r in rows), 0)


def check_target(a, b, n, target_parts, expected):
    target = Counter(target_parts)
    for order in range(4):
        blocks, lengths = layout(a, b, n, order)
        got = 0
        if target.keys() <= lengths:
            got = prod(block_count(rows, groups, target, n)
                       for _, rows, groups in blocks)
        assert got == expected, (a, b, n, target, expected, got, order)
        if order == 0:
            reduced = 0
            if target.keys() <= lengths:
                reduced = prod(count_rows_only(a, b, n, rows, groups, target)
                               for _, rows, groups in blocks)
            assert reduced == expected, (a, b, n, target, expected, reduced)


def prefix_ranks(vectors):
    """Exact rational ranks, keeping the elimination entirely integral."""
    basis = {}
    yield 0
    for vector in vectors:
        work = dict(vector)
        while work:
            pivot = min(work)
            if pivot not in basis:
                basis[pivot] = work
                break
            other = basis[pivot]
            x, y = work[pivot], other[pivot]
            work = {j: y * work.get(j, 0) - x * other.get(j, 0)
                    for j in work.keys() | other.keys()}
            work = {j: v for j, v in work.items() if v}
        yield len(basis)


def cut_profile(rows, groups):
    """Check both the number of shared equations and the span dimension."""
    entries = []
    for col, terms in groups:
        for t, row, coeff_r, coeff_c in terms:
            vector = Counter({row: coeff_r})
            vector[col] += coeff_c
            entries.append((col, t, vector))
    vectors = [v for _, _, v in entries]
    ranks_left = list(prefix_ranks(vectors))
    ranks_right = list(reversed(list(prefix_ranks(reversed(vectors)))))
    left, right = Counter(), Counter()
    for _, _, vector in entries:
        right.update(vector.keys())
    for i, (col, _, vector) in enumerate(entries):
        labels = vector.keys()
        left.update(labels)
        right.subtract(labels)
        active = set(left & right)
        assert active <= set(rows) | {col}
        assert len(active) <= len(rows) + 1
        intersection_dim = ranks_left[i + 1] + ranks_right[i + 1] - ranks_left[-1]
        assert intersection_dim <= len(rows)
        if len(active) == len(rows) + 1:
            yield (i + 1, tuple(sorted(active)), tuple(t for _, t, _ in entries))


def check_cuts():
    cuts = 0
    sharp = None
    for n in range(1, 81):
        for a in range(2, 10):
            for b in range(a + 1, 13):
                for order in range(4):
                    blocks, _ = layout(a, b, n, order)
                    for label, rows, groups in blocks:
                        cuts += sum(len(terms) for _, terms in groups)
                        for pos, active, sequence in cut_profile(rows, groups):
                            if (sharp is None and b % a and label in rows
                                    and any(c == label for c, _ in groups)
                                    and len(active) == min(tau(a // gcd(a, b)), tau(b // gcd(a, b))) + 1):
                                sharp = (a, b, n, order, label, pos, active, sequence)
        if sharp:
            break
    assert sharp is not None
    print("cut prefixes checked:", cuts)
    print("saturated cut (a,b,n,order,block,position,active,lengths):", sharp)


def check_budget_identity():
    """Check all bounded prefixes at the saturated cut before its shared column."""
    a, b, n = 4, 3, 6
    blocks, _ = layout(a, b, n, 0)
    label, rows, groups = next(block for block in blocks if block[0] == 12)
    assert label in rows and groups[-1][0] == label
    assert [(c, [t for t, *_ in terms]) for c, terms in groups] == [
        (3, [1, 3]), (6, [2, 6]), (12, [4])]
    tested = 0
    for x1, x3, x2 in product(range(n + 1), range(n // 3 + 1), range(n // 2 + 1)):
        row = Counter()
        col = Counter()
        for t, x in ((1, x1), (3, x3), (2, x2)):
            row[lcm(a, t)] += gcd(a, t) * x
            col[lcm(b, t)] += gcd(b, t) * x
        # Column 3 has ended. On a surviving state its total is beta_3.
        weighted_rows = sum(r * row[r] for r in rows)
        assert b * weighted_rows == a * (3 * col[3] + 6 * col[6])
        tested += 1
    print("weighted-sum identities at the saturated cut:", tested)


def main():
    small = Counter()
    for a in range(1, 5):
        for b in range(a, 5):
            for n in range(1, 18 // (a + b) + 1):
                counts = Counter(product_type((a, b), x) for x in partitions(n))
                for target in partitions(n * (a + b)):
                    expected = counts[target]
                    check_target(a, b, n, target, expected)
                    small["all"] += 1
                    small["yes" if expected else "no"] += 1
                    small["multiple"] += int(expected > 1)
    print("all targets, 1 <= a <= b <= 4 and (a+b)n <= 18:", dict(small))
    positive = Counter()
    for a in range(1, 7):
        for b in range(a, 7):
            for n in range(1, 9):
                counts = Counter(product_type((a, b), x) for x in partitions(n))
                for target, expected in counts.items():
                    check_target(a, b, n, target, expected)
                    positive["targets"] += 1
                    positive["classes"] += expected
                    positive["multiple"] += int(expected > 1)
    print("positive targets, 1 <= a <= b <= 6, 1 <= n <= 8:", dict(positive))
    examples = [
        (4, 6, 14, (4,) * 2 + (6,) * 2 + (12,) * 10, [(1, 1, 12), (2, 12)]),
        (2, 3, 12, (2,) * 6 + (3,) * 6 + (6,) * 5, [(1,) * 6 + (6,), (2, 2, 2, 3, 3)]),
        (4, 7, 4, (4,) * 4 + (28,), [(4,)]),
    ]
    for a, b, n, target, expected in examples:
        solutions = [x for x in partitions(n) if product_type((a, b), x) == target]
        assert solutions == expected
        check_target(a, b, n, target, len(solutions))
        print("worked example", (a, b, n), "quotient partitions:", solutions)
    print("row-only DP matched all targets above")
    check_cuts()
    check_budget_identity()
    print("PASS")


if __name__ == "__main__":
    main()
