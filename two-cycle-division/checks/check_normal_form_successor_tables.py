"""Check the consecutive normal form on literal successor tables.

Run beside check_marker_normal_form.py with stock Python 3.9+:
    python -B check_normal_form_successor_tables.py

The test uses the companion script for the forward transformation.
It computes quotient classes before and after the transformation by
enumerating integer partitions and tracing the actual product maps.
"""

from collections import Counter, defaultdict
from functools import lru_cache
from math import gcd, lcm
from pathlib import Path
from runpy import run_path

marked_target = run_path(str(Path(__file__).with_name("check_marker_normal_form.py")))["marked_target"]


def periods(parts):
    return tuple(sorted(Counter(parts).items()))


def successor(parts):
    out = []
    for size in parts:
        offset = len(out)
        out.extend(range(offset + 1, offset + size))
        out.append(offset)
    return out


def trace(permutation):
    seen = bytearray(len(permutation))
    sizes = []
    for start in range(len(permutation)):
        if seen[start]:
            continue
        x, size = start, 0
        while not seen[x]:
            seen[x] = 1
            x = permutation[x]
            size += 1
        assert x == start
        sizes.append(size)
    return periods(sizes)


def product_type(left, parts):
    right = successor(parts)
    # Pair (i,j) has label i * len(right) + j.
    table = [i1 * len(right) + j1 for i1 in left for j1 in right]
    return trace(table)


def partitions(n, first=1):
    if n == 0:
        yield ()
    for t in range(first, n + 1):
        for rest in partitions(n - t, t):
            yield (t,) + rest


def target_size(target):
    return sum(period * count for period, count in target)


@lru_cache(maxsize=None)
def small_fibres(axes, n):
    table = successor(axes)
    answers = defaultdict(set)
    for parts in partitions(n):
        answers[product_type(table, parts)].add(parts)
    return answers


def fibre(axes, target, stats):
    size, rem = divmod(target_size(target), sum(axes))
    if rem or not size:
        return set()
    if size <= 12:
        return small_fibres(axes, size).get(target, set())
    # The single planned collision case has quotient size 36. Avoid
    # storing thousands of unrelated product types.
    assert size == 36
    table = successor(axes)
    found = set()
    for parts in partitions(size):
        stats["partitions traced in the size-36 test"] += 1
        if product_type(table, parts) == target:
            found.add(parts)
    return found


def terms(mapping):
    return tuple(sorted((t, k) for t, k in mapping.items() if k))


def check(axes, target, stats, large=False):
    new_axes, mapped, info, _, (old_shared, new_shared) = marked_target(
        axes, target, use_markers=False
    )
    output_n, rem = divmod(target_size(mapped), sum(new_axes))
    if not large and not rem and output_n > 12:
        return
    old_solutions = fibre(axes, target, stats)
    new_solutions = fibre(new_axes, mapped, stats)
    D = lcm(*axes)
    expected = set()
    for solution in old_solutions:
        translated = []
        for t in solution:
            PF, QF = info[lcm(D, t)]
            translated.append(t * PF * QF)
        expected.add(tuple(sorted(translated)))
    assert len(expected) == len(old_solutions)
    assert expected == new_solutions, (axes, target, new_axes, mapped, expected, new_solutions)
    if info:
        mass = defaultdict(int)
        for ell, copies in target:
            mass[lcm(D, ell)] += ell * copies
        quotient_mass = sum((mass[R] // sum(axes)) * PF * QF for R, (PF, QF) in info.items())
        assert target_size(mapped) == quotient_mass * sum(new_axes)
    stats["validated inputs"] += 1
    stats["positive" if expected else "negative"] += 1
    stats["complete old quotient classes"] += len(old_solutions)
    stats["different solutions at same target"] += len(old_solutions) > 1
    stats["old row and column collide"] += old_shared
    stats["new row and column collide"] += new_shared
    stats["nontrivially rescaled block"] += any(x * y > 1 for x, y in info.values())
    stats["both rescaling factors > 1"] += any(min(x, y) > 1 for x, y in info.values())
    stats["negative despite integral shares"] += bool(info and not expected)
    stats["equal axes"] += axes[0] == axes[1]
    stats["gcd > 1"] += gcd(*axes) > 1
    stats["unit axis"] += 1 in axes
    if large:
        print("Collision fixture:", axes, target, "->", new_axes, mapped, "quotients", expected)


def main():
    stats = Counter()
    for a in range(1, 7):
        for b in range(1, 7):
            axes = a, b
            targets = set()
            for size in range(1, 8):
                targets.update(small_fibres(axes, size))
            # Perturb target multiplicities to test rejected raw maps too.
            for target in list(targets):
                first = target[0][0]
                case = Counter(dict(target))
                case[first] += 1
                targets.add(terms(case))
            if axes == (2, 3):
                targets.add(((2, 1), (3, 2), (6, 2)))
            for target in sorted(targets):
                check(axes, target, stats)
    check((1, 1), ((6, 2),), stats, large=True)
    assert all(stats[key] for key in [
        "positive", "negative", "old row and column collide",
        "new row and column collide", "both rescaling factors > 1",
        "negative despite integral shares", "different solutions at same target",
    ])
    for key, value in stats.items():
        print(f"{key}: {value}")
    print("PASS")


if __name__ == "__main__":
    main()
