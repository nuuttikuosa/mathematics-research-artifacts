"""Count all cycle sums at fixed size, including fixed points, by splitting on one prime."""
from functools import lru_cache
from itertools import product
from math import prod

@lru_cache(None)
def count(primes, size, marginals):
    if not primes:
        return 1  # All remaining states must be fixed points.
    if any(value < 0 or value > size or value % p
           for p, value in zip(primes, marginals, strict=True)):
        return 0
    best = None
    # The recurrence works with any split prime. Choose the one giving
    # the fewest candidate vectors for the side coprime to that prime.
    for i, prime in enumerate(primes):
        coprime_size = size - marginals[i]
        others = [(p, value) for p, value in zip(primes, marginals, strict=True)
                  if p != prime]
        choices = []
        for p, value in others:
            lower = max(0, value - marginals[i])
            upper = min(value, coprime_size)
            step = prime * p
            start = lower + (value - lower) % step
            choices.append(range(start, upper + 1, step))
        candidate_count = prod(map(len, choices))
        if candidate_count == 0:
            return 0
        if best is None or candidate_count < best[0]:
            best = candidate_count, prime, coprime_size, others, choices
    _, prime, coprime_size, others, choices = best
    remaining = tuple(p for p, _ in others)
    reduced_size = (size - coprime_size) // prime
    total = 0
    for left in product(*choices):
        ways = count(remaining, coprime_size, left)
        if ways:
            right = tuple((value - used) // prime
                          for (_, value), used in zip(others, left, strict=True))
            total += ways * count(remaining, reduced_size, right)
    return total
