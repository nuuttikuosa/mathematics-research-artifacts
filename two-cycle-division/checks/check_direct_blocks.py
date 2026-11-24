"""Exhaustive small check of the unprojected block equations for two cycles.

Run with Python 3.9 or later. Only the standard library is used.
    python -B check_direct_blocks.py
    python -B check_direct_blocks.py --output direct_blocks_results.json

For each divisor below and each quotient size from 1 to 3, the test visits
every target cycle type of the required size. The oracle constructs the
Cartesian successor table for every possible quotient. The solver separately
builds the gcd/lcm equations, keeps their zero rows, and solves each block
by bounded integer enumeration. It does not impose a total quotient size
equation. Every witness it returns is checked with an actual Cartesian product.
"""

import argparse
import json
from collections import Counter, defaultdict
from itertools import product
from math import gcd, lcm
from pathlib import Path


def partitions(total, smallest=1):
    """Yield every unordered cycle type of the given size in canonical order."""
    if total == 0:
        yield ()
    for first in range(smallest, total + 1):
        for rest in partitions(total - first, first):
            yield (first,) + rest


def successor_table(cycles):
    table = []
    for length in cycles:
        start = len(table)
        table.extend(range(start + 1, start + length))
        table.append(start)
    return table


def cycle_type(table):
    seen = set()
    lengths = []
    for start in range(len(table)):
        if start in seen:
            continue
        vertex = start
        length = 0
        while vertex not in seen:
            seen.add(vertex)
            length += 1
            vertex = table[vertex]
        if vertex != start:
            raise ValueError("Input is not a permutation")
        lengths.append(length)
    return tuple(sorted(lengths))


def cartesian_cycle_type(a, b, quotient):
    first = successor_table((a, b))
    second = successor_table(quotient)
    width = len(second)
    table = [i * width + j for i in first for j in second]
    return cycle_type(table)


def solve_block(variables, rows, target, limit):
    """Enumerate bounded multiplicities for the variables in this block."""
    for values in product(range(limit + 1), repeat=len(variables)):
        assignment = dict(zip(variables, values))
        for length, coefficients in rows.items():
            total = sum(assignment[t] * weight for t, weight in coefficients.items())
            if total != target[length]:
                break
        else:
            return assignment
    return None


def direct_blocks(a, b, target_table):
    n, remainder = divmod(len(target_table), a + b)
    if remainder or n < 1:
        return "size", None
    target = Counter(cycle_type(target_table))
    base = lcm(a, b)
    variables = defaultdict(list)
    rows = defaultdict(dict)
    reachable = set()
    for t in range(1, n + 1):
        block = lcm(base, t)
        variables[block].append(t)
        for divisor_length in (a, b):
            length = lcm(divisor_length, t)
            if lcm(base, length) != block:
                raise AssertionError("Variable and row received different block labels")
            reachable.add(length)
            # Adding is essential when both divisor cycles give the same row.
            coefficients = rows[block].setdefault(length, Counter())
            coefficients[t] += gcd(divisor_length, t)
    if set(target) - reachable:
        return "unreachable_target_length", None
    quotient = []
    for block, candidates in variables.items():
        # Counter returns zero on absent lengths, so their equations stay present.
        assignment = solve_block(candidates, rows[block], target, n)
        if assignment is None:
            return "infeasible_block", None
        for t, count in assignment.items():
            quotient.extend([t] * count)
    return "yes", tuple(sorted(quotient))


def exhaustive_case(a, b, n):
    products = {cartesian_cycle_type(a, b, x) for x in partitions(n)}
    counts = Counter()
    for target in partitions((a + b) * n):
        status, witness = direct_blocks(a, b, successor_table(target))
        if (witness is not None) != (target in products):
            raise AssertionError((a, b, n, target, status, witness))
        if witness is not None:
            if sum(witness) != n or cartesian_cycle_type(a, b, witness) != target:
                raise AssertionError(("invalid witness", a, b, target, witness))
        counts[status] += 1
    return {
        "divisor_cycles": [a, b],
        "quotient_size": n,
        "target_size": (a + b) * n,
        "all_target_types": sum(counts.values()),
        "yes": counts["yes"],
        "infeasible_block": counts["infeasible_block"],
        "unreachable_target_length": counts["unreachable_target_length"],
    }


def run():
    cases = [
        exhaustive_case(a, b, n)
        for a, b in ((2, 3), (2, 4), (3, 3), (4, 6))
        for n in range(1, 4)
    ]
    # With target 5 C_2, the length-2 row alone would accept x_1=1, x_2=2.
    # The length-3 and length-6 zero rows must reject that spurious choice.
    focused = []
    for name, a, b, target, expected in (
        ("opening_example", 2, 3, (3, 6, 6), "infeasible_block"),
        ("two_contributions_to_length_six", 2, 3, (6,) * 5, "yes"),
        ("equal_divisor_cycles", 3, 3, (3, 3), "yes"),
        ("absent_length_must_be_enforced", 2, 3, (2,) * 5, "infeasible_block"),
    ):
        status, witness = direct_blocks(a, b, successor_table(target))
        if status != expected:
            raise AssertionError((name, status, expected))
        n = sum(target) // (a + b)
        matches = [x for x in partitions(n) if cartesian_cycle_type(a, b, x) == target]
        if (witness is None) != (not matches):
            raise AssertionError((name, witness, matches))
        if witness is not None and witness not in matches:
            raise AssertionError((name, "witness mismatch", witness))
        focused.append({"name": name, "status": status, "oracle_quotients": matches})
    return {
        "cases": cases,
        "totals": {
            key: sum(case[key] for case in cases)
            for key in ("all_target_types", "yes", "infeasible_block", "unreachable_target_length")
        },
        "focused": focused,
        "status": "PASS",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(run(), indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
