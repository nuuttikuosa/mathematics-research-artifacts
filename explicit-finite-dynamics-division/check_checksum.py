"""Independent finite checks of the checksum address construction.

Run with Python 3.10+ (standard library only):
    python check_checksum.py --output check_checksum.json
No transition tables or quotient histograms are constructed. We enumerate the
actual capacity-one sites, intersect their coordinate fibers with the support,
and solve the marked binary constraints from their numerical interval bounds.
"""

import argparse
import json
from collections import Counter, defaultdict, deque
from itertools import count, product
from math import isqrt, prod
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def parameters(t):
    h = next(x for x in count(3) if x**x >= t)
    ceil_log2_h = (h - 1).bit_length()
    u = max(3, (h + ceil_log2_h - 1) // ceil_log2_h)
    g = next(x for x in count(1) if u**x >= t)
    return h, u, g


def slot(u, g):
    return [
        word for word in product(range(u), repeat=g + 1)
        if (word[-1] - sum(word[:-1])) % u in (0, 1)
    ]


def key(point, axis):
    return axis, point[:axis] + point[axis + 1:]


def actual_fibers(sites):
    fibers = defaultdict(list)
    for i, site in enumerate(sites):
        for axis in range(len(site)):
            fibers[key(site, axis)].append(i)
    return fibers


def reachable(graph, seed):
    queue = deque([seed])
    seen = {seed}
    while queue:
        i = queue.popleft()
        for j in graph[i]:
            if j not in seen:
                seen.add(j)
                queue.append(j)
    return seen


def check_slot(u, g):
    sites = slot(u, g)
    graph = [[] for _ in sites]
    fibers = actual_fibers(sites)
    for ends in fibers.values():
        require(len(ends) == 2, "slot fiber cardinality")
        i, j = ends
        types = [(sites[k][-1] - sum(sites[k][:-1])) % u for k in ends]
        require(sorted(types) == [0, 1], "slot edge parity")
        graph[i].append(j)
        graph[j].append(i)
    require(len(sites) == 2 * u**g, "wrong slot cardinality")
    require(len(reachable(graph, 0)) == len(sites), "disconnected slot")
    return {"u": u, "g": g, "sites": len(sites), "two_point_fibers": len(fibers)}


def primes():
    for n in count(3, 2):
        if all(n % p for p in range(3, isqrt(n) + 1, 2)):
            yield n


def stable_walk(start, end, u):
    """An explicit path in the slot graph between stable words."""
    digits = list(start[:-1])
    require(start[-1] == sum(digits) % u and end[-1] == sum(end[:-1]) % u, "unstable terminal")
    yield start
    for j, target in enumerate(end[:-1]):
        while digits[j] != target:
            yield tuple(digits) + ((1 + sum(digits)) % u,)
            digits[j] = (digits[j] + 1) % u
            yield tuple(digits) + (sum(digits) % u,)
    require(tuple(digits) + (sum(digits) % u,) == end, "wrong end of path")


def check_formula(n, clauses, trim=False):
    t = 2 * n + 1
    h, u, g = parameters(t)
    width = g + 1
    plane_axis = 3 * width
    d = plane_axis + 1
    prime_source = primes()
    groups = [tuple(next(prime_source) for _ in range(u)) for _ in range(plane_axis)]
    groups.append(tuple(next(prime_source) for _ in range(3)))
    words = slot(u, g)
    stable = [word for word in words if word[-1] == sum(word[:-1]) % u]
    literals = [x for i in range(1, n + 1) for x in (i, -i)]
    address = dict(zip(literals + [0], stable))
    require(len(address) == t, "out of stable words")
    link_tests = []
    for i in range(1, n + 1):
        for j, k in product(range(3), repeat=2):
            if j != k:
                ids = [0, 0, 0]
                ids[j], ids[k] = i, -i
                link_tests.append(tuple(ids))
    kept = defaultdict(set)
    if trim:
        for test in link_tests + list(clauses):
            for plane, literal in enumerate(test):
                if not literal:
                    continue
                a, b = (address[test[j]] for j in range(3) if j != plane)
                owner = kept[plane, literal]
                # Route using one varying slot at a time from (pad,pad).
                owner.update((w, address[0]) for w in stable_walk(address[0], a, u))
                owner.update((a, w) for w in stable_walk(address[0], b, u))
    sites = []
    owners = []
    seeds = {}
    for plane in range(3):
        for literal in literals:
            seeds[plane, literal] = len(sites)
            positions = sorted(kept[plane, literal]) if trim else product(words, repeat=2)
            for a, b in positions:
                slots = [a, b]
                slots.insert(plane, address[literal])
                sites.append(sum(slots, ()) + (plane,))
                owners.append((plane, literal))
    require(len(set(sites)) == len(sites), "coincident sites")
    fibers = actual_fibers(sites)
    marked = {}

    def mark(line, lo, hi, use):
        new = lo, hi, use
        if line in marked:
            require(marked[line] == new and use == "clause", "marked fiber collision")
        marked[line] = new

    def weights(line):
        ax = line[0]
        return [groups[ax][sites[i][ax]] for i in fibers[line]]

    # Construct the copying trees by scanning actual ambient fibers. Reject
    # any third endpoint or owner crossing, including on unused graph edges.
    graph = [[] for _ in sites]
    for line, ends in fibers.items():
        ax = line[0]
        plane = sites[ends[0]][-1]
        if ax == plane_axis or ax // width == plane:
            continue
        if trim and len(ends) == 1:
            continue
        require(len(ends) == 2, "copying line with extra/missing sites")
        i, j = ends
        require(owners[i] == owners[j], "copying line reaches another owner")
        graph[i].append((j, line))
        graph[j].append((i, line))
    seen = set()
    for owner, seed in seeds.items():
        queue = deque([seed])
        seen.add(seed)
        count_owner = 1
        while queue:
            i = queue.popleft()
            for j, line in graph[i]:
                if j in seen:
                    continue
                require(owners[j] == owner, "wrong BFS owner")
                ps = weights(line)
                mark(line, min(ps), max(ps), "copy")
                seen.add(j)
                count_owner += 1
                queue.append(j)
        expected_sites = len(kept[owner]) if trim else len(words)**2
        require(count_owner == expected_sites and expected_sites > 0, "owner disconnected")
    require(len(seen) == len(sites), "unreached sites")

    def plane_line(ids):
        return plane_axis, sum((address[i] for i in ids), ())

    for ids in link_tests:
        line = plane_line(ids)
        required_planes = {j for j, literal in enumerate(ids) if literal}
        require({sites[a][-1] for a in fibers[line]} == required_planes, "variable line")
        ps = weights(line)
        mark(line, min(ps), max(ps), "link")
    for clause in clauses:
        line = plane_line(clause)
        require({sites[i][-1] for i in fibers[line]} == {0, 1, 2}, "clause line")
        ps = weights(line)
        mark(line, min(ps), sum(ps), "clause")

    # Decode marked intervals without the slot parity or concept of an owner.
    # The numerical tables for two-point faces generate the parity graph.
    equality_graph = [[] for _ in sites]
    clauses_decoded = []
    weighted_mass = sum(prod(gr[v] for gr, v in zip(groups, site)) for site in sites)
    width_mass = (1 << d) * weighted_mass
    for line, (lo, hi, use) in marked.items():
        endpoints, ps = fibers[line], weights(line)
        require(2 * lo <= sum(ps) <= 2 * hi, "half assignment infeasible")
        choices = {
            bits for bits in product(range(2), repeat=len(ps))
            if lo <= sum(p * bit for p, bit in zip(ps, bits)) <= hi
        }
        ax, fixed = line
        fixed_groups = groups[:ax] + groups[ax + 1:]
        face_period = prod(gr[v] for gr, v in zip(fixed_groups, fixed))
        width_mass -= (sum(ps) - (hi - lo)) * face_period
        if use == "clause":
            require(len(choices) == 7, "weighted clause truth table")
            clauses_decoded.append((endpoints, choices))
        else:
            require(choices == {(0, 1), (1, 0)}, "weighted XOR truth table")
            i, j = endpoints
            equality_graph[i].append(j)
            equality_graph[j].append(i)

    component = [-1] * len(sites)
    parity = [0] * len(sites)
    n_components = 0
    for seed in range(len(sites)):
        if component[seed] >= 0:
            continue
        component[seed] = n_components
        stack = [seed]
        while stack:
            i = stack.pop()
            for j in equality_graph[i]:
                if component[j] < 0:
                    component[j], parity[j] = n_components, 1 - parity[i]
                    stack.append(j)
                require(component[j] == component[i] and parity[i] != parity[j], "odd cycle")
        n_components += 1
    require(n_components == n, "wrong parity component count")
    # Canonical sites lie in positive-owner plane 0 and have stable addresses
    # in both varying slots. Directly compare the decoded assignments as a multiset.
    index = {site: i for i, site in enumerate(sites)}
    reads = [index[address[i] + address[0] + address[0] + (0,)] for i in range(1, n + 1)]
    decoded = Counter()
    for base in product(range(2), repeat=n_components):
        bits = [parity[i] ^ base[component[i]] for i in range(len(sites))]
        if all(tuple(bits[i] for i in ids) in choices for ids, choices in clauses_decoded):
            decoded[tuple(bits[i] for i in reads)] += 1
    expected = Counter(
        bits for bits in product(range(2), repeat=n)
        if all(any(bool(bits[abs(lit) - 1]) == (lit > 0) for lit in clause) for clause in clauses)
    )
    require(decoded == expected, "SAT assignments and integral arrays differ")
    divisor_size = 2 + sum(prod(group) for group in groups)
    target_size = 2 * divisor_size * width_mass
    return {
        "construction": "trimmed" if trim else "full",
        "n": n, "clauses": clauses, "h": h, "u": u, "g": g, "axes": d,
        "active_sites": len(sites), "nonempty_axis_fibers": len(fibers),
        "marked_fibers": dict(Counter(use for _, _, use in marked.values())),
        "xor_components": n_components, "half_assignment": True,
        "satisfying_assignments": sum(expected.values()),
        "feasible_integral_arrays": sum(decoded.values()),
        "assignments_in_variable_order": ["".join(map(str, bits)) for bits in sorted(decoded)],
        "divisor_states": divisor_size, "target_states_from_interval_widths": target_size
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = [
        (1, [(1, -1, 1)]), (1, [(1, 1, 1)]), (1, [(1, 1, 1), (-1, -1, -1)]),
        (2, [(1, 2, 2), (-1, -2, -2), (1, 2, 2)]),
        (2, [(1, 1, 2), (1, 1, -2), (-1, -1, 2), (-1, -1, -2)]),
        (3, [(1, 2, 3), (-1, -2, 3), (-3, 2, 1)]),
    ]
    result = {
        "scope": "Enumerated capacity-one sites, their actual axis fibers, and integral arrays; table sizes by the proved width formula. No explicit transition table was materialized.",
        "slot_probes": [check_slot(u, g) for u, g in [(3, 1), (3, 2), (4, 2), (5, 3)]],
        "formulas": [check_formula(n, clauses) for n, clauses in cases],
        "trimmed_formulas": [check_formula(n, clauses, trim=True) for n, clauses in cases],
    }
    report = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(report, encoding="utf-8")
    print(report, end="")


if __name__ == "__main__":
    main()
