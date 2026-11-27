"""Targeted checker for final literal blocks and exit-layer occupancy.

No import of prior constructors/checkers. This tests local geometric
properties over all block partitions of small reflected Hamilton paths.
It does not test full SAT instances or replace the general proof.
"""
from itertools import product
from collections import defaultdict
from pathlib import Path
import json


def reflected_words(u, g):
    path = [()]
    for _ in range(g):
        old = path
        path = [(a,) + suffix for a in range(u)
                for suffix in (old if a % 2 == 0 else old[::-1])]
    return path


def stable(a, u):
    return a + (sum(a) % u,)


def intermediate(a, u):
    return a + ((sum(a) + 1) % u,)


def line_key(point, coordinate):
    return (coordinate, point[:coordinate], point[coordinate + 1:])


def check(u, g, path, cutmask):
    owners = {}
    edges = set()
    block = 0
    for i, a in enumerate(path):
        if i and cutmask & (1 << (i - 1)):
            block += 1
        s = stable(a, u)
        assert s not in owners
        owners[s] = block
    for i, (a, b) in enumerate(zip(path, path[1:])):
        changed = [j for j in range(g) if a[j] != b[j]]
        assert len(changed) == 1
        j = changed[0]
        assert abs(a[j] - b[j]) == 1
        if cutmask & (1 << i):
            continue
        lower = a if a[j] < b[j] else b
        mid = intermediate(lower, u)
        sa, sb = stable(a, u), stable(b, u)
        assert owners[sa] == owners[sb]
        assert mid not in owners or owners[mid] == owners[sa]
        owners[mid] = owners[sa]
        edges.add(tuple(sorted((sa, mid))))
        edges.add(tuple(sorted((mid, sb))))

    occupied_lines = defaultdict(set)
    for point in owners:
        for j in range(g + 1):
            occupied_lines[line_key(point, j)].add(point)
    neighbours = defaultdict(set)
    for x, y in edges:
        different = [j for j in range(g + 1) if x[j] != y[j]]
        assert len(different) == 1
        assert occupied_lines[line_key(x, different[0])] == {x, y}
        assert (x[-1] == sum(x[:-1]) % u) != (y[-1] == sum(y[:-1]) % u)
        neighbours[x].add(y)
        neighbours[y].add(x)

    # Output additions are intermediate addresses only. At exit every
    # stable address has one main-lane and one buffer-lane occurrence.
    # Vary buffer choices independently: their actual permutation relation
    # is unnecessary for this local final-layer assertion.
    for selector in (0, 1, 2):
        exit_points = {p + (0,) for p in owners}
        for i, a in enumerate(path):
            lane = 1 + ((a[0] + selector * i) % u)
            exit_points.add(stable(a, u) + (lane,))
        lane_lines = defaultdict(set)
        for p in exit_points:
            lane_lines[p[:-1]].add(p[-1])
        for a in path:
            assert len(lane_lines[stable(a, u)]) == 2
        for p in owners:
            if p[-1] != sum(p[:-1]) % u:
                assert lane_lines[p] == {0}

    # Every connected component has exactly one intended block, and each
    # intended block is connected, even if intermediate vertices coalesce.
    remaining = set(owners)
    found = []
    while remaining:
        todo = [next(iter(remaining))]
        component = set()
        while todo:
            p = todo.pop()
            if p in component:
                continue
            component.add(p)
            todo.extend(neighbours[p] - component)
        remaining -= component
        labels = {owners[p] for p in component}
        assert len(labels) == 1
        found.append(next(iter(labels)))
    assert len(found) == len(set(found)) == block + 1
    return len(owners), len(edges)


families = []
total = 0
for u, g in [(3, 1), (4, 1), (5, 1), (3, 2), (4, 2), (3, 3), (3, 4)]:
    full = reflected_words(u, g)
    counts = 0
    if len(full) <= 16:
        # Exhaust every partition of every nonempty active prefix. The
        # unused suffix models inactive dummy outputs.
        inputs = [(full[:n], range(1 << (n - 1))) for n in range(1, len(full) + 1)]
        exhaustive = True
    else:
        masks = {0, (1 << (len(full) - 1)) - 1}
        masks.update(1 << i for i in range(len(full) - 1))
        masks.add(sum(1 << i for i in range(0, len(full) - 1, 2)))
        inputs = [(full, sorted(masks))]
        exhaustive = False
    for path, masks in inputs:
        for mask in masks:
            check(u, g, path, mask)
            counts += 1
    families.append(dict(u=u, g=g, words=len(full), cases=counts,
                         all_prefix_partitions=exhaustive))
    total += counts

result = dict(cases=total, families=families,
              properties=['block owner separation', 'all marked address-line endpoints',
                          'type coloring', 'exit-lane occupancy', 'exact block components'],
              independent_of_previous_checkers=True, all_assertions_passed=True,
              scope='Local output-block/exit-layer checks only; not a full SAT-instance census')
out = Path(__file__).with_name('blocks_geometry_results.json')
out.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps(result, indent=2))
