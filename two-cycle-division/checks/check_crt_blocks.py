"""Test the two-by-two representative solver against direct partition enumeration."""
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import product
from math import gcd, lcm


def freeze(items):
    return tuple(sorted((length, qty) for length, qty in items.items() if qty))


def least_period(axis, output):
    if output % axis:
        return None
    ratio = output // axis
    return ratio * gcd(axis, pow(ratio, axis.bit_length(), axis))


def common_residue(residue, modulus, target, other):
    """Intersect x=residue mod modulus and x=target mod other."""
    divisor = gcd(modulus, other)
    if (target - residue) % divisor:
        return None
    reduced = other // divisor
    if reduced == 1:
        return residue, modulus
    step = ((target - residue) // divisor * pow(modulus // divisor, -1, reduced)) % reduced
    return residue + step * modulus, modulus * reduced


@dataclass
class Block:
    first: int
    step: int
    last: int
    cells: list

    def vector_at(self, x):
        return freeze({s: (offset + sign * x) // s
                       for s, offset, sign in self.cells if s})

    def size(self):
        return (self.last - self.first) // self.step + 1


def solve_block(g, R, target):
    a, b = 2*g, 3*g
    rows, cols = (R//3, R), (R//2, R)
    if set(target) - (set(rows) | set(cols)):
        return None
    total = sum(length * qty for length, qty in target.items())
    n, rem = divmod(total, a+b)
    if rem:
        return None
    row0, col0 = target.get(rows[0], 0), target.get(cols[0], 0)
    row1, rem = divmod(a*n - rows[0]*row0, R)
    col1 = target.get(R, 0) - row1
    if rem or min(row1, col1) < 0 or cols[0]*col0 + R*col1 != b*n:
        return None
    row = [r*k//a for r,k in zip(rows, (row0, row1))]
    col = [c*k//b for c,k in zip(cols, (col0, col1))]
    assert sum(row) == sum(col) == n
    # The four cell masses are x, row0-x, col0-x and row1-col0+x.
    offsets = [0, row[0], col[0], row[1]-col[0]]
    signs = [1, -1, -1, 1]
    lo, hi = max(0, col[0]-row[1]), min(row[0], col[0])
    residue, modulus = 0, 1
    cells = []
    for (r,c), offset, sign in zip(product(rows, cols), offsets, signs):
        s = lcm(least_period(a, r), least_period(b, c))
        if lcm(a, s) != r or lcm(b, s) != c or lcm(6*g, s) != R:
            s = 0
            fixed = -sign * offset
            lo, hi = max(lo, fixed), min(hi, fixed)
        else:
            combined = common_residue(residue, modulus, -sign * offset, s)
            if combined is None:
                return None
            residue, modulus = combined
        cells.append((s, offset, sign))
    first = lo + (residue - lo) % modulus
    if first > hi:
        return None
    last = first + (hi-first)//modulus * modulus
    return Block(first, modulus, last, cells)


def blocks(g, target):
    groups = defaultdict(dict)
    for length, qty in target:
        groups[lcm(6*g, length)][length] = qty
    result = [solve_block(g, R, counts) for R, counts in groups.items()]
    return None if any(block is None for block in result) else result


def partitions(n, low=1):
    if n == 0:
        yield ()
    for part in range(low, n+1):
        for tail in partitions(n-part, part):
            yield (part,) + tail


def direct_product(g, parts):
    output = Counter()
    for t in parts:
        for axis in (2*g, 3*g):
            output[lcm(axis, t)] += gcd(axis, t)
    return freeze(output)


def direct_normalisation(g, parts):
    result = Counter()
    for t in parts:
        u, v = lcm(2*g,t), lcm(3*g,t)
        # Independent search; only used for the small partition reference.
        s = next(i for i in range(1,t+1) if (lcm(2*g,i),lcm(3*g,i)) == (u,v))
        result[s] += t//s
    return freeze(result)


def all_vectors(solved):
    if solved is None:
        return set()
    by_block = [[b.vector_at(x) for x in range(b.first,b.last+1,b.step)] for b in solved]
    return {tuple(sorted(term for vector in group for term in vector))
            for group in product(*by_block)}


def threshold(g, prime_divisors, solved):
    if solved is None:
        return 0
    if any(b.size() > 1 for b in solved):
        return 2
    for block in solved:
        for s, z in block.vector_at(block.first):
            H = g // gcd(g, s)
            if any(H % p == 0 and p <= z for p in prime_divisors):
                return 2
    return 1


def small():
    totals = Counter()
    # The threshold procedure receives the factorisation of g just once.
    primes = {1: (), 2: (2,), 3: (3,), 4: (2,), 5: (5,), 6: (2,3), 7: (7,), 8: (2,)}
    for g in range(1, 9):
        for n in range(1, 13):
            reference = defaultdict(list)
            for parts in partitions(n):
                reference[direct_product(g,parts)].append(parts)
            targets = set(reference)
            for original in reference:
                for r,k in original:
                    if k >= 2:
                        changed = Counter(dict(original))
                        changed[r] -= 2
                        changed[2*r] += 1
                        targets.add(freeze(changed))
                    if r % 2 == 0:
                        changed = Counter(dict(original))
                        changed[r] -= 1
                        changed[r//2] += 2
                        targets.add(freeze(changed))
            for target in targets:
                answers = reference.get(target, ())
                expected = {direct_normalisation(g, parts) for parts in answers}
                solved = blocks(g,target)
                assert all_vectors(solved) == expected, (g,target,expected)
                assert threshold(g,primes[g],solved) == min(2,len(answers)), (g,target,answers)
                totals['cases'] += 1
                totals['positive'] += bool(answers)
                totals['several vectors'] += len(expected)>1
                totals['several lifts of one vector'] += len(expected)==1 and len(answers)>1
    print("Small checks (g=1..8, quotient size 1..12):", dict(sorted(totals.items())))


def large():
    w, m = 2**240+1, 2**170
    target = ((2*w,6*m), (3*w,6*m), (6*w,5*m))
    (block,) = blocks(1,target)
    assert [s for s, _, _ in block.cells] == [w,2*w,3*w,6*w]
    assert (block.first,block.step,block.last,block.size()) == (0,6*w,6*w*m,m+1)
    assert block.vector_at(0) == ((2*w,3*m),(3*w,2*m))
    assert block.vector_at(block.last) == ((w,6*m),(6*w,m))
    print("241-bit w and 171-bit m: four shared moduli; exactly 2^170+1 normalized vectors.")
    print("Only the interval, stride and two endpoint witnesses were computed for this family.")


if __name__ == '__main__':
    small()
    large()
    print("PASS")
