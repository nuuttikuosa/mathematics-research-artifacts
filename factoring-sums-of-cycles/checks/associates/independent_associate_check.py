"""Exact counterexample search for tensor associates of permutation characteristic polynomials.

No earlier checker or project implementation is imported. Polynomial tuples
and rational root phases give two separate finite representations. The
unbounded conclusion is proved in the paper's related-work section, not
inferred from this bounded search.
"""
from collections import Counter
from fractions import Fraction
from itertools import chain
from math import gcd
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent


def partitions(total, minimum=1):
    """All integer partitions, as nondecreasing cycle-length tuples."""
    if total == 0:
        yield ()
    else:
        for first in range(minimum, total + 1):
            for rest in partitions(total - first, first):
                yield (first,) + rest


def characteristic(cycles):
    """Direct coefficient multiplication of product (x**d - 1)."""
    coefficients = [1]
    for length in cycles:
        new = [0] * (len(coefficients) + length)
        for degree, value in enumerate(coefficients):
            new[degree] -= value
            new[degree + length] += value
        coefficients = new
    assert coefficients[-1] == 1
    assert len(coefficients) == sum(cycles) + 1
    return tuple(coefficients)


def negative_root_twist(coefficients):
    """Monic polynomial whose roots are negatives of the original roots."""
    degree = len(coefficients) - 1
    return tuple(value * (-1 if (degree - j) % 2 else 1)
                 for j, value in enumerate(coefficients))


def root_phases(cycles):
    """Roots of unity represented by exact angles in Q/Z, with multiplicity."""
    return Counter(chain.from_iterable(
        (Fraction(j, length) for j in range(length)) for length in cycles))


def phase_key(multiset):
    return tuple(sorted(multiset.items()))


def rotate(phases, angle):
    new = Counter()
    for old, multiplicity in phases.items():
        new[(old + angle) % 1] += multiplicity
    return new


coefficient_degrees = []
negative_admissible = 0
for n in range(1, 31):
    cases = list(partitions(n))
    database = {}
    for cycles in cases:
        polynomial = characteristic(cycles)
        assert polynomial not in database, ('Characteristic collision', cycles)
        database[polynomial] = cycles
    hits = 0
    for polynomial, cycles in database.items():
        # Identity is the other rational root-of-unity twist.
        assert database[polynomial] == cycles
        image = negative_root_twist(polynomial)
        matched = database.get(image)
        all_even = all(d % 2 == 0 for d in cycles)
        assert (matched is not None) == all_even, (cycles, matched)
        if matched is not None:
            hits += 1
            assert matched == cycles and image == polynomial
    negative_admissible += hits
    coefficient_degrees.append(dict(states=n, cycle_types=len(cases),
                                    admissible_negative_twists=hits))

phase_degrees = []
rotation_total = 0
rotation_admissible = 0
for n in range(1, 13):
    cases = list(partitions(n))
    database = {}
    for cycles in cases:
        key = phase_key(root_phases(cycles))
        assert key not in database
        database[key] = cycles
    # If g is a degree-n permutation characteristic polynomial, a root eta of g has order
    # at most n. These are all possible root-of-unity rotations of that order.
    angles = sorted({Fraction(0)} | {
        Fraction(j, order) for order in range(2, n + 1)
        for j in range(1, order) if gcd(j, order) == 1})
    hits = 0
    for cycles in cases:
        phases = root_phases(cycles)
        common_divisor = 0
        for length in cycles:
            common_divisor = gcd(common_divisor, length)
        for angle in angles:
            rotated = rotate(phases, angle)
            matched = database.get(phase_key(rotated))
            expected_stabilizer = common_divisor % angle.denominator == 0
            assert (matched is not None) == expected_stabilizer, (cycles, angle, matched)
            if matched is not None:
                hits += 1
                assert matched == cycles and rotated == phases
    count = len(cases) * len(angles)
    rotation_total += count
    rotation_admissible += hits
    phase_degrees.append(dict(states=n, cycle_types=len(cases), angles=len(angles),
                             rotations=count, admissible_rotations=hits))

# Characteristic zero is essential for recovery of cycle type from the
# polynomial. In characteristic p, C_p and p*C_1 have equal polynomials.
characteristic_p_controls = []
for p in (2, 3, 5, 7):
    one_cycle = (p,)
    fixed_points = (1,) * p
    left = tuple(c % p for c in characteristic(one_cycle))
    right = tuple(c % p for c in characteristic(fixed_points))
    assert one_cycle != fixed_points and left == right
    assert characteristic(one_cycle) != characteristic(fixed_points)
    characteristic_p_controls.append(dict(characteristic=p, left_cycles=one_cycle,
                                          right_cycles=fixed_points,
                                          equal_mod_p=True, equal_over_Q=False))

assert negative_root_twist(characteristic((1,))) == (1, 1)
assert negative_root_twist(characteristic((2,))) == characteristic((2,))
assert negative_root_twist(characteristic((2, 2))) == characteristic((2, 2))

result = dict(
    date='2026-11-27',
    coefficient_search=dict(maximum_states=30,
        nonempty_cycle_types=sum(row['cycle_types'] for row in coefficient_degrees),
        rational_root_of_unity_twists=2 * sum(row['cycle_types'] for row in coefficient_degrees),
        admissible_negative_twists=negative_admissible, distinct_type_collisions=0,
        by_degree=coefficient_degrees),
    exact_phase_search=dict(maximum_states=12,
        nonempty_cycle_types=sum(row['cycle_types'] for row in phase_degrees),
        rotations=rotation_total, admissible_rotations=rotation_admissible,
        distinct_type_collisions=0, by_degree=phase_degrees),
    characteristic_p_boundary=characteristic_p_controls,
    all_assertions_passed=True, earlier_project_code_imported=False,
    floating_point_roots_used=False,
    scope='Finite coefficient/phase checks; general proof is in the paper',
    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
)
(ROOT / 'independent_associate_results.json').write_text(
    json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps({
    'coefficient_types': result['coefficient_search']['nonempty_cycle_types'],
    'coefficient_twists': result['coefficient_search']['rational_root_of_unity_twists'],
    'admissible_negative_twists': negative_admissible,
    'phase_types': result['exact_phase_search']['nonempty_cycle_types'],
    'exact_rotations': rotation_total, 'admissible_rotations': rotation_admissible,
    'distinct_type_collisions': 0, 'all_assertions_passed': True,
}, indent=2))
