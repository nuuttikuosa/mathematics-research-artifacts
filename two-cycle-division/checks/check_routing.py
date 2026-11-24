"""Check the graph-routing step of Theorem 6.1 on K4 and the Petersen graph.

Run from the artifact directory:
    python -B checks/check_routing.py

Build the actual simple cubic graph with new vertices for every diamond.
Enumerate a diamond's possible two boundary colours from its internal edges.
Contract the resulting equality constraints and count all assignments of
the remaining boundary variables. Recover and check an expanded witness.
The arithmetic encoding and the successor tables are outside this check.
"""

import json
from collections import Counter, defaultdict
from itertools import combinations, product
from math import prod


def canonical_edge(v, w):
    return tuple(sorted((v, w)))


def edge_addresses(v, w, e):
    x = (v, 0)
    addresses = [x]
    for target in [(v, e), (w, e), (w, 0)]:
        for axis in (0, 1):
            changed = x[axis] ^ target[axis]
            for bit in range(changed.bit_length()):
                if changed & (1 << bit):
                    new = list(x)
                    new[axis] ^= 1 << bit
                    x = tuple(new)
                    addresses.append(x)
    assert len(set(addresses)) == len(addresses)
    return addresses


def expand(n, edges):
    addresses = [(v, 0) for v in range(n)]
    expanded = []
    diamonds = []
    for i, (v, w) in enumerate(edges, 1):
        previous = v
        for address in edge_addresses(v, w, i)[1:-1]:
            diamond = tuple(range(len(addresses), len(addresses) + 4))
            a, b, c, d = diamond
            diamonds.append(diamond)
            addresses.extend([address] * 4)
            expanded.extend([(a, c), (a, d), (b, c), (b, d), (c, d)])
            expanded.append(canonical_edge(previous, a))
            previous = b
        expanded.append(canonical_edge(previous, w))
    assert len(set(expanded)) == len(expanded)
    assert all(v < w for v, w in expanded)
    return addresses, expanded, diamonds


def incidences(count, edges):
    rows = [[] for _ in range(count)]
    for e, (v, w) in enumerate(edges):
        rows[v].append(e)
        rows[w].append(e)
    assert all(len(row) == 3 for row in rows)
    return rows


def boundary_relation(vertices, rows):
    local = [rows[v] for v in vertices]
    uses = Counter(e for row in local for e in row)
    ports = sorted(e for e, m in uses.items() if m == 1)
    internal = sorted(e for e, m in uses.items() if m == 2)
    assert len(ports) == 2 and len(internal) == 5
    options = defaultdict(list)
    for colours in product(range(3), repeat=7):
        assignment = dict(zip(ports + internal, colours))
        if all(len({assignment[e] for e in row}) == 3 for row in local):
            options[colours[:2]].append({e: assignment[e] for e in internal})
    return ports, options


def colourings(groups):
    # Each actual degree-three vertex requires three distinct edge colours.
    conflicts = defaultdict(set)
    for group in groups:
        if len(set(group)) < 3:
            return 0, None
        for v in group:
            conflicts[v].update(set(group) - {v})
    assigned = {}
    witness = None

    def visit():
        nonlocal witness
        if len(assigned) == len(conflicts):
            if witness is None:
                witness = assigned.copy()
            return 1
        var, choices = min(
            ((v, {0, 1, 2} - {assigned[w] for w in conflicts[v] if w in assigned})
             for v in conflicts if v not in assigned),
            key=lambda item: (len(item[1]), item[0]),
        )
        total = 0
        for choice in sorted(choices):
            assigned[var] = choice
            total += visit()
        assigned.pop(var, None)
        return total

    total = visit()
    return total, witness


def check(n, edges):
    edges = sorted(canonical_edge(v, w) for v, w in edges)
    original_count, _ = colourings(incidences(n, edges))
    addresses, expanded, diamonds = expand(n, edges)
    rows = incidences(len(addresses), expanded)
    distances = Counter(
        sum(bin(x ^ y).count("1") for x, y in zip(addresses[v], addresses[w]))
        for v, w in expanded
    )
    assert max(distances) <= 1
    congestion = Counter(addresses)
    assert max(congestion.values()) <= 16
    # Read each local constraint from expanded edge incidence lists, independent
    # of the routing paths and their original-edge indices.
    components = list(range(len(expanded)))

    def find(i):
        while components[i] != i:
            components[i] = components[components[i]]
            i = components[i]
        return i

    records = []
    for vertices in diamonds:
        ports, options = boundary_relation(vertices, rows)
        assert set(options) == {(0, 0), (1, 1), (2, 2)}
        assert all(len(v) == 2 for v in options.values())
        components[find(ports[0])] = find(ports[1])
        records.append((ports, options))

    groups = [[find(e) for e in row] for row in rows[:n]]
    external = [e for ports, _ in records for e in ports]
    assert {find(e) for e in external} == {v for row in groups for v in row}
    boundary_count, boundary_assignment = colourings(groups)
    assert original_count == boundary_count
    if boundary_assignment is not None:
        assignment = {e: boundary_assignment[find(e)] for e in set(external)}
        for ports, options in records:
            choice = tuple(assignment[e] for e in ports)
            internal = options[choice][0]
            assert set(assignment).isdisjoint(internal)
            assignment.update(internal)
        assert set(assignment) == set(range(len(expanded)))
        assert all({assignment[e] for e in row} == {0, 1, 2} for row in rows)
    return {
        "original_vertices": n, "original_edges": len(edges),
        "private_diamonds": len(diamonds), "expanded_vertices": len(addresses),
        "expanded_edges": len(expanded), "occupied_addresses": len(congestion),
        "max_vertices_per_address": max(congestion.values()),
        "edges_by_address_distance": dict(sorted(distances.items())),
        "boundary_patterns_per_diamond": ["00", "11", "22"],
        "extensions_per_boundary_pattern": 2,
        "direct_original_colourings": original_count,
        "expanded_boundary_assignments": boundary_count,
        "expanded_colourings_by_local_elimination": boundary_count * prod(
            len(options[(0, 0)]) for _, options in records),
        "expanded_witness_checked": boundary_assignment is not None,
    }


if __name__ == "__main__":
    petersen = []
    for i in range(5):
        petersen.extend([(i, (i + 1) % 5), (i, i + 5),
                         (i + 5, 5 + (i + 2) % 5)])
    results = {"K4": check(4, combinations(range(4), 2)),
               "Petersen": check(10, petersen)}
    print(json.dumps(results, indent=2))
