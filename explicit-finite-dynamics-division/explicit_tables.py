"""Three explicit transition-table checks of the interval-realization lemma."""

from collections import Counter
from itertools import product
from math import prod

from cycle_reference import cycle_solutions, multiply


def successor_table(cycles):
    table = []
    for length, count in sorted(cycles.items()):
        if count < 0 or length < 1:
            raise ValueError("Invalid cycle histogram")
        for _ in range(count):
            start = len(table)
            table.extend(range(start + 1, start + length))
            table.append(start)
    return table


def read_cycles(table):
    if sorted(table) != list(range(len(table))):
        raise ValueError("Expected a permutation table")
    seen = set()
    lengths = Counter()
    for start in range(len(table)):
        if start in seen:
            continue
        point = start
        length = 0
        while point not in seen:
            seen.add(point)
            length += 1
            point = table[point]
        lengths[length] += 1
    return lengths


def formal_quotient(bit, interval):
    # The only site of capacity one is (3,7). The face (1,7) reads 3*bit.
    coefficients = Counter()
    for face in product((1, 3, 5), (1, 7, 11)):
        period = prod(face)
        weight = 21 // period if face[0] in (1, 3) and face[1] in (1, 7) else 0
        low, high = interval if face == (1, 7) else (0, weight)
        total = weight * bit
        z = high - total if face.count(1) % 2 else total - low
        coefficients[2 * period] += z
        coefficients[period] += 2 * (high - low - z)
    return {d: n for d, n in coefficients.items() if n}


def verify(interval, expected_bits):
    divisor = successor_table({2: 1, 15: 1, 77: 1})
    divisor_periods = list(read_cycles(divisor).elements())
    target = successor_table(multiply(divisor_periods, formal_quotient(0, interval)))
    observed = read_cycles(target)
    solutions, _, _, _ = cycle_solutions(divisor_periods, observed)
    expected = {tuple(sorted(formal_quotient(bit, interval).items())) for bit in expected_bits}
    assert solutions == expected
    for solution in solutions:
        quotient = successor_table(dict(solution))
        table = [a * len(quotient) + q for a in divisor for q in quotient]
        assert read_cycles(table) == observed
    print(f"interval {interval}, |A|={len(divisor)}, |B|={len(target)}, quotients={len(solutions)}")


if __name__ == "__main__":
    verify((0, 3), [0, 1])
    verify((3, 3), [1])
    verify((1, 2), [])
