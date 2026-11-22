"""Three-plane 3SAT reduction and a decoder for its cycle histograms.

Python 3.10+, no dependencies. Run with no arguments for the three small tests.
compile_formula(variables, clauses) returns (A, B, axes, statistics).
decode_small(A, B) independently recovers all solutions for small tests.
Neither function expands the histograms to successor arrays.
"""
from collections import Counter
from itertools import product
from math import gcd, isqrt, lcm, prod
from pathlib import Path
import json


ORDINARY = 2


def primes():
    candidate = 3
    while True:
        if all(candidate % i for i in range(3, isqrt(candidate) + 1, 2)):
            yield candidate
        candidate += 2


def multiplier(a, x):
    result = Counter()
    for u, m in a.items():
        for v, n in x.items():
            result[lcm(u, v)] += gcd(u, v) * m * n
    return {k: v for k, v in result.items() if v}


def with_coordinate(cell, i, value):
    return cell[:i] + (value,) + cell[i + 1:]


def down_closure(values, axes):
    # Cell coordinate 1 denotes the absent coordinate. Every physical label >1.
    faces = {e: v for e, v in values.items() if v}
    for i in range(len(axes)):
        additions = Counter()
        for e, v in faces.items():
            additions[with_coordinate(e, i, 1)] += e[i] * v
        faces.update(additions)
    return faces


def lift(widths, intervals, f, denominator=1):
    result = {}
    for e in widths:
        lo, hi = intervals[e]
        value = f.get(e, 0)
        z = (denominator * hi - value if e.count(1) % 2
             else value - denominator * lo)
        bottom = ORDINARY * (denominator * (hi-lo) - z)
        d = prod(e)
        if z:
            result[ORDINARY*d] = z
        if bottom:
            result[d] = bottom
    return result


def compile_formula(variables, clauses):
    if type(variables) is not int or variables < 1:
        raise ValueError("choose a positive number of explicitly listed variables")
    clauses = {tuple(c) for c in clauses}
    if any(len(c) != 3 or any(type(i) is not int or not 1 <= abs(i) <= variables
                              for i in c) for c in clauses):
        raise ValueError("expected three signed, nonzero literals per clause")
    symbols = 2 * variables + 1
    radix = 3
    while radix**radix < symbols:
        radix += 1
    places = 1
    while radix**places < symbols:
        places += 1
    stream = primes()
    axes = [tuple(next(stream) for _ in range(radix)) for _ in range(6*places)]
    axes.append(tuple(next(stream) for _ in range(3)))
    a = {ORDINARY: 1, **{prod(group): 1 for group in axes}}

    def diagonal(symbol):
        word = []
        for _ in range(places):
            symbol, r = divmod(symbol, radix)
            word.extend((r, r))
        return tuple(word)

    def symbol_id(literal):
        return 2*(abs(literal)-1) + int(literal < 0)

    def wire_site(owner, fixed, ranks):
        blocks = []
        j = 0
        for slot in range(3):
            if slot == owner:
                blocks.extend(fixed)
            else:
                for h in ranks[j:j+places]:
                    blocks.extend((h // 2, (h+1) // 2))
                j += places
        return tuple(group[index] for group, index in zip(axes, blocks+[owner]))

    active = {}
    bands = {}

    def mark_xor(s, t):
        j, = [i for i in range(len(axes)) if s[i] != t[i]]
        face = with_coordinate(s, j, 1)
        assert face not in bands
        bands[face] = (min(s[j], t[j]), max(s[j], t[j]))

    for owner in range(3):
        for symbol in range(2*variables):
            fixed = diagonal(symbol)
            for ranks in product(range(2*radix-1), repeat=2*places):
                e = wire_site(owner, fixed, ranks)
                active[e] = 1
                if any(ranks):
                    j = next(i for i, r in enumerate(ranks) if r)
                    previous = with_coordinate(ranks, j, ranks[j]-1)
                    mark_xor(e, wire_site(owner, fixed, previous))

    def plane_face(symbols_in_slots):
        blocks = sum((diagonal(s) for s in symbols_in_slots), ())
        return tuple(group[i] for group, i in zip(axes, blocks)) + (1,)

    for variable in range(1, variables+1):
        for j, k in product(range(3), repeat=2):
            if j == k:
                continue
            ids = [2*variables]*3
            ids[j], ids[k] = symbol_id(variable), symbol_id(-variable)
            e = plane_face(ids)
            plane_axis = len(axes) - 1
            mark_xor(with_coordinate(e, plane_axis, axes[-1][j]),
                     with_coordinate(e, plane_axis, axes[-1][k]))
    for clause in clauses:
        e = plane_face([symbol_id(literal) for literal in clause])
        assert e not in bands
        bands[e] = (min(axes[-1]), sum(axes[-1]))

    widths = down_closure(active, axes)
    assert bands.keys() <= widths.keys()
    intervals = {e: bands.get(e, (0, w)) for e, w in widths.items()}
    assert all(0 <= lo <= hi <= widths[e] for e, (lo, hi) in intervals.items())
    zero = lift(widths, intervals, {})
    b = multiplier(a, zero)
    assert b and min(b.values()) > 0
    doubled_fractional_witness = lift(widths, intervals, widths, denominator=2)
    assert min(doubled_fractional_witness.values()) >= 0
    assert multiplier(a, doubled_fractional_witness) == {k: 2*v for k, v in b.items()}
    stats = {"axes": len(axes), "radix": radix, "address_places": places,
             "free_sites": len(active), "marked_lines": len(bands),
             "A_states": sum(k*v for k, v in a.items()),
             "B_states": sum(k*v for k, v in b.items()), "B_types": len(b),
             "negative_coefficients_in_signed_reference": sum(v < 0 for v in zero.values()),
             "checked_A_divides_two_copies_of_B": True}
    return a, b, axes, stats


def factor_group(q):
    factors = []
    for p in primes():
        if p*p > q:
            break
        if q % p == 0:
            factors.append(p)
            q //= p
            if q % p == 0:
                raise ValueError("the period is not squarefree")
    if q > 1:
        factors.append(q)
    return tuple(factors)


def decode_small(a, b):
    # This routine receives no formula, gates, capacities, encoded intervals,
    # chosen logical assignment or signed reference.
    axes = sorted((factor_group(q) for q in a if q != ORDINARY), key=min)
    assert {ORDINARY: 1, **{prod(s): 1 for s in axes}} == a
    assert all(len(s) >= 2 for s in axes)
    assert len(set.union(*(set(s) for s in axes))) == sum(map(len, axes))
    cells = sorted(product(*[(1,)+s for s in axes]), key=lambda e: e.count(1))
    capacities, base, bounds = {}, {}, {}
    possible_support = set()
    for e in cells:
        d = prod(e)
        coeff = b.get(ORDINARY*d, 0)
        cap, rem = divmod(coeff, ORDINARY)
        assert rem == 0
        capacities[e] = cap
        possible_support.add(ORDINARY*d)
        missing = [i for i, v in enumerate(e) if v == 1]
        for i in missing:
            out = d*prod(axes[i])
            possible_support.update((out, ORDINARY*out))
        if not missing:
            base[e] = 0
        else:
            i = missing[0]
            zline = b.get(ORDINARY*d*prod(axes[i]), 0)
            base[e] = zline - sum(h * base[with_coordinate(e, i, h)] for h in axes[i])
        z = base[e]
        bounds[e] = (z-cap, z) if len(missing) % 2 else (-z, cap-z)
    assert b.keys() <= possible_support
    active = {e: capacities[e] for e in cells if 1 not in e and capacities[e]}
    assert all(m == 1 for m in active.values())
    widths = down_closure(active, axes)
    graph = {e: [] for e in active}
    ors = 0
    for e in cells:
        lo, hi = bounds[e]
        width = widths.get(e, 0)
        assert 0 <= lo <= hi <= width
        if (lo, hi) == (0, width):
            continue
        assert e.count(1) == 1
        i = e.index(1)
        hits = [p for h in axes[i] if (p := with_coordinate(e, i, h)) in active]
        weights = [p[i] for p in hits]
        if len(hits) == 2 and (lo, hi) == (min(weights), max(weights)):
            p, q = hits
            graph[p].append(q)
            graph[q].append(p)
        elif len(hits) == 3 and (lo, hi) == (min(weights), sum(weights)):
            ors += 1
        else:
            raise ValueError("decoded a non-wire, non-OR marked line")
    phases = {}
    components = 0
    for root in active:
        if root in phases:
            continue
        stack = [root]
        phases[root] = (components, 0)
        while stack:
            vertex = stack.pop()
            parity = phases[vertex][1]
            for other in graph[vertex]:
                wanted = components, 1-parity
                if other not in phases:
                    phases[other] = wanted
                    stack.append(other)
                assert phases[other] == wanted
        components += 1
    num = 0
    for bits in product(range(2), repeat=components):
        values = {p: bits[j] ^ sign for p, (j, sign) in phases.items()}
        f = down_closure(values, axes)
        if any(not lo <= f.get(e, 0) <= hi for e, (lo, hi) in bounds.items()):
            continue
        x = lift(widths, bounds, f)
        assert all(v >= 0 for v in x.values()) and multiplier(a, x) == b
        num += 1
    return {"decoded_cells": len(cells), "decoded_OR_lines": ors,
            "decoded_XOR_components": components, "quotients": num}


def main():
    tests = []
    for name, clauses in [("tautological", [(1, -1, 1)]),
                          ("forces_true", [(1, 1, 1)]),
                          ("contradiction", [(1, 1, 1), (-1, -1, -1)])]:
        a, b, _, stats = compile_formula(1, clauses)
        decoded = decode_small(a, b)
        count = sum(all(any((lit > 0) == y for lit in c) for c in clauses)
                    for y in (False, True))
        assert decoded["quotients"] == count
        tests.append({"case": name, "clauses": clauses, **stats, **decoded,
                      "source_satisfying_assignments": count})
    output = json.dumps(tests, indent=2) + "\n"
    Path(__file__).with_name("three_plane_cycles.json").write_text(output, encoding="utf-8")
    print(output, end="")


if __name__ == "__main__":
    main()
