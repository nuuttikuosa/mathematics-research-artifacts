"""Check the address multipliers of Lemma 6.2 for one to four bits."""

import json
from itertools import combinations_with_replacement
from math import prod


def check(dimension, modulus):
    relation_length = 30
    coefficient_bound = 7
    base = modulus * (1 + relation_length * coefficient_bound * 2 ** (dimension * relation_length))
    moment_modulus = base ** relation_length
    q = 1 + 10**8 * dimension * moment_modulus
    factors = [q + (pow(1 + base, 2**i, moment_modulus) - 1) % moment_modulus
               for i in range(dimension)]
    multipliers = [
        prod(factors[i] if address & (1 << i) else q for i in range(dimension))
        for address in range(2**dimension)
    ]
    scale = q**dimension
    common_multiple = scale * prod(factors)
    neighbours = 0
    for address, value in enumerate(multipliers):
        assert value % moment_modulus == pow(1 + base, address, moment_modulus)
        assert value % modulus == 1 % modulus
        assert 10**7 * scale <= 10**7 * value < (10**7 + 1) * scale
        assert common_multiple % value == 0
        for i, factor in enumerate(factors):
            if not address & (1 << i):
                shared, remainder = divmod(value, q)
                assert remainder == 0 and shared * factor == multipliers[address | (1 << i)]
                neighbours += 1

    # Test distinct sums for multisets of size zero, one and two.
    seen = {0, *multipliers}
    assert len(seen) == 1 + len(multipliers)
    pair_count = 0
    for first, second in combinations_with_replacement(multipliers, 2):
        total = first + second
        assert total not in seen
        seen.add(total)
        pair_count += 1
    return len(multipliers), neighbours, pair_count


def main():
    dimensions = (1, 2, 3, 4)
    moduli = (1, 4, 12, 35)
    results = [check(d, m) for d in dimensions for m in moduli]
    print(json.dumps({
        "dimensions": dimensions,
        "moduli": moduli,
        "parameter_sets": len(results),
        "addresses_checked": sum(row[0] for row in results),
        "neighbour_pairs_checked": sum(row[1] for row in results),
        "two_term_multisets_checked": sum(row[2] for row in results),
        "status": "PASS",
    }, indent=2))


if __name__ == "__main__":
    main()
