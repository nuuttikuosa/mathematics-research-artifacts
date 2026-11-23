"""Compute the maximum support among all admissible four-prime marginals."""
from math import prod
import numpy as np


def shift(offset):
    """Two array slices with indices separated by offset."""
    target = tuple(slice(i, None) for i in offset)
    source = tuple(slice(None, -i if i else None) for i in offset)
    return target, source


def main():
    if not __debug__:
        raise RuntimeError("Run with Python assertions enabled (omit -O).")
    primes = (2, 3, 5, 7)
    shape = tuple(210 // p for p in primes)
    infinity = 1_000_000
    # At most fifteen distinct nontrivial divisors occur in a cycle sum.
    grid = np.full(shape, infinity, dtype=np.int32)
    grid[0, 0, 0, 0] = 0
    for mask in range(1, 16):
        d = prod(p for i, p in enumerate(primes) if mask & (1 << i))
        offset = tuple(d // p if d % p == 0 else 0 for p in primes)
        if any(i >= limit for i, limit in zip(offset, shape, strict=True)):
            continue
        # Pay for the support entry only on the first copy.
        copies = np.full_like(grid, infinity)
        target, source = shift(offset)
        copies[target] = np.minimum(infinity, grid[source] + 16 * d + 1)
        # Binary doubling allows every further multiplicity that fits.
        step = 1
        while all(step * i < limit for i, limit in zip(offset, shape, strict=True)):
            target, source = shift(tuple(step * i for i in offset))
            copies[target] = np.minimum(copies[target], copies[source] + step * 16 * d)
            step *= 2
        grid = np.minimum(grid, copies)

    assert grid.max() < infinity
    support = grid % 16
    maximum = int(support.max())
    attainers = int((support == maximum).sum())
    assert (grid.size, maximum, attainers) == (9_261_000, 12, 98)
    assert divmod(int(grid[80, 39, 33, 25]), 16) == (246, 12)
    print(f"k=4  vectors={grid.size}  maximum={maximum}  attainers={attainers}")


if __name__ == "__main__":
    main()
