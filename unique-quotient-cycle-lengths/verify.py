"""Verify the five finite examples in the table, using integer arithmetic."""
from itertools import product
import json
from math import prod
from pathlib import Path

from check_record_cert import check
from prime_split import count


DATA = Path(__file__).resolve().parent / "data"
PRIMES = (2, 3, 5, 7, 11, 13, 17, 19)


def marginals(primes, witness):
    return tuple(
        sum(length * copies for length, copies in witness.items() if length % p == 0)
        for p in primes
    )


def check_lift(primes, witness, extension):
    """Verify the added cycles by counting every congruent size-t realisation."""
    q = extension["prime"]
    assert q == PRIMES[len(primes)]
    added = extension["added_lengths"]
    assert len(set(added)) == len(added)
    assert all(length > 1 and prod(primes) % length == 0 for length in added)
    size = sum(added)
    target = marginals(primes, dict.fromkeys(added, 1))
    ranges = [
        range(value % (p * q), size + 1, p * q)
        for p, value in zip(primes, target, strict=True)
    ]
    profiles = 0
    ways = 0
    for profile in product(*ranges):
        profiles += 1
        realisations = count(primes, size, profile)
        if realisations:
            assert profile == target
            ways += realisations
    assert ways == 1
    count.cache_clear()

    lifted = {q * length: copies for length, copies in witness.items()}
    lifted.update(dict.fromkeys(added, 1))
    new_primes = (*primes, q)
    new_size = sum(length * copies for length, copies in lifted.items())
    new_marginals = marginals(new_primes, lifted)
    assert all(0 < value < min(prod(new_primes), new_size) for value in new_marginals)
    return lifted, (q, profiles, ways, len(lifted), new_size)


def main():
    if not __debug__:
        raise RuntimeError("Run with Python assertions enabled (omit -O).")
    for filename, k, support, size in (
        ("record_cert_4.json", 4, 12, 246),
        ("record_cert_21.json", 5, 21, 2051),
        ("record_cert_30.json", 6, 30, 20185),
    ):
        result = check(DATA / filename)
        assert (result["k"], result["support"], result["n"], result["solutions"]) == (
            k, support, size, 1
        )
        assert all(0 < value < min(prod(PRIMES[:k]), size) for value in result["s"])
        print(f"k={k}  support={support}  size={size}  solutions=1")

    base = json.loads((DATA / "record_cert_30.json").read_text(encoding="utf-8"))
    witness = dict(base["witness"])
    extensions = json.loads((DATA / "lifts.json").read_text(encoding="utf-8"))
    expected = [(17, 486, 1, 43, 343463), (19, 2640, 1, 58, 6526232)]
    for k, (extension, expected_row) in enumerate(zip(extensions, expected, strict=True), 7):
        witness, actual = check_lift(PRIMES[: k - 1], witness, extension)
        assert actual == expected_row
        q, profiles, ways, support, size = actual
        print(f"k={k}  support={support}  size={size}  solutions={ways}  "
              f"profiles={profiles}  prime={q}")


if __name__ == "__main__":
    main()
