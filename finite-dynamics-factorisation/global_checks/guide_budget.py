"""Independent exhaustive checks for the physical-path budget.

Run with Python 3: python guide_budget.py
Only the standard library is used. The true-product test uses planted
coordinates solely as an oracle for the number of active children.
"""
from collections import defaultdict
from itertools import product
from math import prod


def cycle_data(table):
    cycles = set()
    height = 0
    for start in range(len(table)):
        path = []
        node = start
        while node not in path:
            path.append(node)
            node = table[node]
        entry = path.index(node)
        cycles.add(frozenset(path[entry:]))
        height = max(height, entry)
    return cycles, height


def powers(table, depth):
    result = [tuple(range(len(table)))]
    for _ in range(depth):
        result.append(tuple(table[v] for v in result[-1]))
    return result


def widths(images, size):
    result = []
    for row in images:
        counts = [0] * size
        for vertex in row:
            counts[vertex] += 1
        result.append(counts)
    return result


def arbitrary_endpoints():
    tables = pairs = levels = unit_levels = 0
    for size in range(1, 6):
        for table in product(range(size), repeat=size):
            tables += 1
            _, height = cycle_data(table)
            images = powers(table, height + 1)
            ws = widths(images, size)
            endpoint_groups = defaultdict(list)
            for z, root in enumerate(images[-1]):
                row = [ws[j][image[z]] for j, image in enumerate(images)]
                assert row[0] == 1 and row == sorted(row)
                caps = [b // a for a, b in zip(row, row[1:])]
                endpoint_groups[root].append(caps)
            for root, group in endpoint_groups.items():
                for first, second in product(group, repeat=2):
                    budgets = [a * b for a, b in zip(first, second)]
                    assert min(budgets) >= 1
                    assert prod(budgets) <= ws[-1][root] ** 2 <= size ** 2
                    pairs += 1
                    levels += len(budgets)
                    unit_levels += budgets.count(1)
    return tables, pairs, levels, unit_levels


def true_guides():
    factors = []
    for size in (2, 3):
        for table in product(range(size), repeat=size):
            cycles, h = cycle_data(table)
            if len(cycles) == 1:
                factors.append((table, next(iter(cycles)), h))
    products = centers = steps = 0
    for (x, cx, hx), (y, cy, hy) in product(factors, repeat=2):
        H = 1 + max(hx, hy)
        ny = len(y)
        table = tuple(x[a] * ny + y[b] for a in range(len(x)) for b in range(ny))
        wx = widths(powers(x, H), len(x))
        wy = widths(powers(y, H), ny)
        wb = widths(powers(table, H), len(table))
        products += 1
        for start_x, start_y in product(cx, cy):
            centers += 1
            xcyc = xmove = start_x
            ycyc = ymove = start_y
            caps = []
            for j in reversed(range(H)):
                next_x = next(v for v in cx if x[v] == xcyc)
                next_y = next(v for v in cy if y[v] == ycyc)
                x_options = [v for v in range(len(x)) if x[v] == xmove and wx[j][v]]
                y_options = [v for v in range(ny) if y[v] == ymove and wy[j][v]]
                min_x = min(x_options, key=lambda v: wx[j][v])
                min_y = min(y_options, key=lambda v: wy[j][v])
                cap_y = wb[j+1][xcyc*ny+ymove] // wb[j][next_x*ny+min_y]
                cap_x = wb[j+1][xmove*ny+ycyc] // wb[j][min_x*ny+next_y]
                assert len(x_options) <= cap_x and len(y_options) <= cap_y
                caps.append(cap_x * cap_y)
                xcyc, ycyc, xmove, ymove = next_x, next_y, min_x, min_y
                steps += 1
            assert prod(caps) <= wb[H][start_x*ny+start_y] ** 2
    return len(factors), products, centers, steps


if __name__ == "__main__":
    print("All labelled inputs of sizes 1..5")
    print("tables=%d endpoint_pairs=%d levels=%d unit_levels=%d" % arbitrary_endpoints())
    print("All ordered products of connected labelled factors of sizes 2 and 3")
    print("factors=%d products=%d cyclic_centers=%d guide_steps=%d" % true_guides())
