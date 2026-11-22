"""Test two-slot routing using geometric fibers and, for n=1, cycle products.

Run next to check_checksum.py with Python 3.10+:
    python check_routed.py --output check_routed.json
    python check_routed.py --blind-decoder three_plane_cycles.py --output check_routed.json
The state tables are represented by exact cycle histograms in the product check.
"""

import argparse
import json
import runpy
from collections import Counter, defaultdict, deque
from itertools import product
from math import gcd, lcm, prod
from pathlib import Path

from check_checksum import actual_fibers, key, parameters, primes, require, slot, stable_walk


def add_links(n):
    tests = []
    for variable in range(1, n + 1):
        for j, k in product(range(3), repeat=2):
            if j != k:
                test = [0, 0, 0]
                test[j], test[k] = variable, -variable
                tests.append(tuple(test))
    return tests


def ambient_check(n, clauses, histograms=False, blind_decoder=None):
    links = add_links(n)
    tests = links + clauses
    t = max(len(tests), 2 * n + 1)
    h, u, g = parameters(t)
    w = g + 1
    stage_ax, plane_ax = 2 * w, 2 * w + 1
    source = primes()
    groups = [
        tuple(next(source) for _ in range(q)) for q in [u] * (2 * w) + [2, 3]
    ]
    stable = [word for word in slot(u, g) if word[-1] == sum(word[:-1]) % u]
    literals = [x for i in range(1, n + 1) for x in (i, -i)]
    address = dict(zip([0] + literals, stable))
    pad = address[0]
    # On each plane, a test is incident with at most one literal. Every
    # clause gets its own test ID, including repeated ordered clauses.
    owners = {}
    terminals = defaultdict(set)
    for test_no, test in enumerate(tests):
        for plane, literal in enumerate(test):
            if not literal:
                continue
            a, b = stable[test_no], address[literal]
            for s in stable_walk(pad, b, u):
                owners[a + s + (0, plane)] = literal
            terminals[plane, literal].add(a)
    for (plane, literal), tests_of_owner in sorted(terminals.items()):
        b = address[literal]
        root = min(tests_of_owner)
        for a in tests_of_owner:
            for s in stable_walk(root, a, u):
                owners[s + b + (1, plane)] = literal
    sites = sorted(owners)
    indices = {site: i for i, site in enumerate(sites)}
    fibers = actual_fibers(sites)
    marked = {}

    def weights(line):
        ax = line[0]
        return [groups[ax][sites[i][ax]] for i in fibers[line]]

    def mark(line, role, expected):
        require(set(fibers[line]) == set(expected), f"ambient fiber collision in {role}")
        require(line not in marked, f"duplicate marked line in {role}")
        ps = weights(line)
        lo = min(ps)
        hi = sum(ps) if role == "clause" else max(ps)
        marked[line] = (lo, hi, role)

    # Compute allowed copying edges directly in the ambient point set. At
    # stage 0 only literal coordinates vary; at stage 1 only test coordinates.
    graph = [[] for _ in sites]
    for line, ends in fibers.items():
        axis = line[0]
        stage = sites[ends[0]][stage_ax]
        allowed = w <= axis < 2 * w if stage == 0 else 0 <= axis < w
        if not allowed:
            continue
        require(1 <= len(ends) <= 2, "too many active wire endpoints")
        if len(ends) == 2:
            i, j = ends
            require(owners[sites[i]] == owners[sites[j]], "wire touches wrong literal")
            graph[i].append((j, line))
            graph[j].append((i, line))
    seen = set()
    wire_components = Counter()
    for seed in range(len(sites)):
        if seed in seen:
            continue
        stage, plane = sites[seed][-2:]
        wire_components[stage] += 1
        seen.add(seed)
        queue = deque([seed])
        while queue:
            i = queue.popleft()
            for j, line in graph[i]:
                if j not in seen:
                    seen.add(j)
                    mark(line, f"wire{stage}", [i, j])
                    queue.append(j)
    incidences = sum(sum(bool(lit) for lit in test) for test in tests)
    require(wire_components[0] == incidences, "disconnected vertical route")
    require(wire_components[1] == len(terminals), "disconnected horizontal route")

    reads = {}
    for test_no, test in enumerate(tests):
        a = stable[test_no]
        on_test = []
        for plane, literal in enumerate(test):
            if not literal:
                continue
            b = address[literal]
            one = a + b + (0, plane)
            two = a + b + (1, plane)
            mark(key(one, stage_ax), "glue", [indices[one], indices[two]])
            read = indices[a + pad + (0, plane)]
            on_test.append(read)
            if plane == 0 and literal > 0:
                reads[literal] = read
        role = "link" if test_no < len(links) else "clause"
        mark(key(a + pad + (0, 0), plane_ax), role, on_test)

    # Solve the tested numerical inequalities. This decoder knows no meaning
    # for test addresses, literals, or stages.
    parity_graph = [[] for _ in sites]
    predicates = []
    for line, (lo, hi, _) in marked.items():
        ids = fibers[line]
        ps = weights(line)
        require(2 * lo <= sum(ps) <= 2 * hi, "rational half solution failed")
        table = {
            xs for xs in product(range(2), repeat=len(ids))
            if lo <= sum(x * p for x, p in zip(xs, ps)) <= hi
        }
        if len(ids) == 2:
            require(table == {(0, 1), (1, 0)}, "weighted copying or link failure")
            i, j = ids
            parity_graph[i].append(j)
            parity_graph[j].append(i)
        else:
            require(len(ids) == 3 and len(table) == 7 and (0, 0, 0) not in table, "clause")
            predicates.append((ids, table))
    component = [-1] * len(sites)
    parity = [0] * len(sites)
    num_components = 0
    for seed in range(len(sites)):
        if component[seed] >= 0:
            continue
        component[seed] = num_components
        stack = [seed]
        while stack:
            i = stack.pop()
            for j in parity_graph[i]:
                if component[j] < 0:
                    component[j], parity[j] = num_components, 1 - parity[i]
                    stack.append(j)
                require(component[i] == component[j] and parity[j] != parity[i], "odd cycle")
        num_components += 1
    require(num_components == n, "wrong XOR component count")
    feasible = []
    decoded = Counter()
    for root_bits in product(range(2), repeat=n):
        bits = tuple(parity[i] ^ root_bits[component[i]] for i in range(len(sites)))
        if all(tuple(bits[i] for i in ids) in table for ids, table in predicates):
            decoded[tuple(bits[reads[j]] for j in range(1, n + 1))] += 1
            feasible.append(bits)
    expected = Counter(
        bits for bits in product(range(2), repeat=n)
        if all(any(bool(bits[abs(x) - 1]) == (x > 0) for x in clause) for clause in clauses)
    )
    require(decoded == expected, "arrays and truth assignments disagree")
    report = {
        "n": n, "clauses": clauses, "test_ids": len(tests), "h": h, "u": u, "g": g,
        "axes": len(groups), "incidences": incidences, "active_sites": len(sites),
        "ambient_nonempty_axis_fibers": len(fibers),
        "marked_fibers": dict(sorted(Counter(v[2] for v in marked.values()).items())),
        "xor_components": num_components,
        "assignments_from_intervals": ["".join(map(str, xs)) for xs in sorted(decoded)],
        "arrays_count": len(feasible), "direct_sat_count": sum(expected.values()),
        "half_assignment": True
    }
    if histograms:
        report["cycle_products"] = check_products(groups, sites, marked, feasible, blind_decoder)
    return report


def check_products(groups, sites, marked, solutions, blind_decoder):
    axes = len(groups)
    prime_groups = [prod(ps) for ps in groups]
    divisor = [2] + prime_groups
    periods = {}
    weighted_sites = defaultdict(list)
    for i, site in enumerate(sites):
        full_period = prod(ps[j] for ps, j in zip(groups, site))
        for mask in range(1 << axes):
            face = tuple(j if mask & (1 << ax) else -1 for ax, j in enumerate(site))
            if face not in periods:
                periods[face] = prod(ps[j] for ps, j in zip(groups, face) if j >= 0)
            weighted_sites[face].append((i, full_period // periods[face]))
    bounds = {e: (0, sum(p for _, p in points)) for e, points in weighted_sites.items()}
    for (ax, rest), (lo, hi, _) in marked.items():
        bounds[rest[:ax] + (-1,) + rest[ax:]] = lo, hi
    expected_size = 2 * sum(divisor) * sum(periods[e] * (hi - lo) for e, (lo, hi) in bounds.items())

    def scaled_image(numerators, scale, feasible):
        y = Counter()
        for e, points in weighted_sites.items():
            f = sum(numerators[i] * p for i, p in points)
            lo, hi = bounds[e]
            z = scale * hi - f if e.count(-1) % 2 else f - scale * lo
            y[2 * periods[e]] += z
            y[periods[e]] += 2 * (scale * (hi - lo) - z)
        if feasible:
            require(min(y.values()) >= 0, "negative feasible quotient")
        result = Counter()
        for a in divisor:
            for b, multiplicity in y.items():
                if multiplicity:
                    result[lcm(a, b)] += gcd(a, b) * multiplicity
        return {d: c for d, c in result.items() if c}

    target = scaled_image([0] * len(sites), 1, False)
    require(min(target.values()) > 0, "negative signed construction")
    size = sum(d * c for d, c in target.items())
    require(size == expected_size, "state-count identity")
    for period in target:
        require(prod(divisor) % period == 0, "period not squarefree or outside labels")
        intersections = [gcd(period, q) for q in prime_groups]
        heavy = [j for j, v in enumerate(intersections) if v != 1 and v not in groups[j]]
        require(len(heavy) <= 1, "period contains two nonsingleton prime groups")
        if heavy:
            j = heavy[0]
            require(intersections[j] == prime_groups[j], "incomplete repeated group")
    require(scaled_image([1] * len(sites), 2, True) == {d: 2 * c for d, c in target.items()},
            "twice half solution does not multiply back")
    for xs in solutions:
        require(scaled_image(xs, 1, True) == target, "integral array does not multiply back")
    result = {
        "divisor_states": sum(divisor), "target_states": size, "target_periods": len(target),
        "positive_shadow_faces": len(weighted_sites), "doubled_half_verified": True,
        "integral_arrays_multiplied_back": len(solutions)
    }
    if blind_decoder is not None:
        decoded = blind_decoder(dict.fromkeys(divisor, 1), target)
        require(decoded["quotients"] == len(solutions), "blind histogram quotient count")
        result["blind_decoding_from_A_B_only"] = decoded
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--blind-decoder", type=Path,
                        help="optional trusted three_plane_cycles.py with decode_small(A,B)")
    args = parser.parse_args()
    blind_decoder = runpy.run_path(str(args.blind_decoder))["decode_small"] if args.blind_decoder else None
    cases = [
        (1, [(1, -1, 1)]), (1, [(1, 1, 1)]), (1, [(1, 1, 1), (-1, -1, -1)]),
        (2, [(1, 2, 2), (-1, -2, -2), (1, 2, 2)]),
        (2, [(1, 1, 2), (1, 1, -2), (-1, -1, 2), (-1, -1, -2)]),
        (3, [(1, 2, 3), (-1, -2, 3), (-3, 2, 1)]),
        (5, [(1, 2, 3), (4, -5, 1), (-3, 2, 4), (-4, -2, -1)]),
    ]
    result = {
        "scope": "Enumerated actual support and all its axis fibers. Numerical marked intervals decoded independently of wire routes. For n=1 the signed target and all satisfying quotients were multiplied as exact cycle histograms. No explicit state tables were allocated.",
        "formulas": [ambient_check(n, clauses, histograms=(n == 1), blind_decoder=blind_decoder)
                     for n, clauses in cases]
    }
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
