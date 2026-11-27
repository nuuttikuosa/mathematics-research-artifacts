"""Independent assertion-based check of a signed-coefficient example.

Uses only the standard library and constructs every term directly.
Checks per-length prefix bounds separately from triangle-inequality bounds.
"""
from fractions import Fraction
from itertools import product
from math import gcd, lcm
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# Build the nine addresses and four full addresses directly in two dimensions.
addresses = list(product((1, 5, 7), (1, 11, 13)))
full = list(product((5, 7), (11, 13)))
omega = 3
anchor = [(3, 1), (35, 1), (143, 1)]
anchor_mass = sum(length * coefficient for length, coefficient in anchor)
endpoints = {}
template = []
for a, b in addresses:
    weight = sum((p * q) // (a * b) for p, q in full
                 if (a == 1 or a == p) and (b == 1 or b == q))
    left, right = (5, 7) if (a, b) == (1, 11) else (0, weight)
    assert 0 <= left <= right <= weight
    zero_count = int(a == 1) + int(b == 1)
    z = right if zero_count % 2 else -left
    template.extend([(omega * a * b, z), (a * b, omega * (right - left - z))])
    endpoints[(a, b)] = (left, right, weight)

assert len({length for length, _ in template}) == len(template)
absolute_mass = sum(length * abs(coefficient) for length, coefficient in template)
endpoint_mass = omega * sum(a * b * (left + right)
                            for (a, b), (left, right, _) in endpoints.items())
endpoint_bound = 2 * omega * 4 * sum(p * q for p, q in full)
assert absolute_mass == endpoint_mass == 3456
assert absolute_mass <= endpoint_bound == 6912

terms = [(lcm(a, d), ac * c * gcd(a, d)) for a, ac in anchor
         for d, c in template if c]
uncollected_mass = sum(length * abs(coefficient) for length, coefficient in terms)
assert uncollected_mass == anchor_mass * absolute_mass == 625536
contributions = {}
for length, coefficient in terms:
    contributions.setdefault(length, []).append(coefficient)

rows = []
for length, coefficients in sorted(contributions.items()):
    partial = 0
    fixed_order_max = 0
    for coefficient in coefficients:
        partial += coefficient
        fixed_order_max = max(fixed_order_max, abs(partial))
        # Exact rational comparison, with no floating-point division.
        assert length * abs(partial) <= uncollected_mass
    positive = sum(c for c in coefficients if c > 0)
    negative_abs = -sum(c for c in coefficients if c < 0)
    all_order_max = max(positive, negative_abs)
    triangle_bound = positive + negative_abs
    assert 0 <= partial
    assert fixed_order_max <= all_order_max <= triangle_bound
    assert length * triangle_bound <= uncollected_mass
    # Any subset can be the prefix of some order. Exhaust all subsets here,
    # independently confirming the positive/negative-sum exact maximum.
    subset_sums = {0}
    for c in coefficients:
        subset_sums |= {s + c for s in tuple(subset_sums)}
    assert max(map(abs, subset_sums)) == all_order_max
    rows.append(dict(length=length, contributions=coefficients,
                     final_coefficient=partial, fixed_order_prefix_max=fixed_order_max,
                     all_orders_prefix_max=all_order_max,
                     sum_of_absolute_contributions=triangle_bound,
                     total_mass_over_length=str(Fraction(uncollected_mass, length))))

final_mass = sum(row['length'] * row['final_coefficient'] for row in rows)
assert final_mass == 565806
maximum_triangle_bound = max(row['sum_of_absolute_contributions'] for row in rows)
actual_fixed_prefix_max = max(row['fixed_order_prefix_max'] for row in rows)
actual_all_orders_prefix_max = max(row['all_orders_prefix_max'] for row in rows)
assert maximum_triangle_bound == 1029
assert actual_fixed_prefix_max <= actual_all_orders_prefix_max <= 1029

result = dict(instance='H1={5,7}, H2={11,13}, omega=3; unit capacities; [L_(1,11),R_(1,11)]=[5,7]',
              anchor_mass=anchor_mass, template_absolute_mass=absolute_mass,
              endpoint_absolute_mass=endpoint_mass, endpoint_bound=endpoint_bound,
              uncollected_absolute_mass=uncollected_mass, final_target_mass=final_mass,
              uncollected_terms=len(terms), per_length=rows,
              maximum_triangle_bound=maximum_triangle_bound,
              actual_fixed_order_prefix_max=actual_fixed_prefix_max,
              actual_all_orders_prefix_max=actual_all_orders_prefix_max,
              length_scaled_bounds_checked=True, all_assertions_passed=True,
              standalone_checker=True)
(ROOT / 'independent_mass_results.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k != 'per_length'}, indent=2))
