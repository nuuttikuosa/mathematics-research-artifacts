"""Independent, integer-only verifier of suffix-cut search certificates."""
import argparse
import heapq
import itertools
import json
import math
from functools import lru_cache
from pathlib import Path


def dot(v, w):
    return sum(a*b for a, b in zip(v, w, strict=True))


def check(file):
    cert = json.loads(file.read_text(encoding="utf-8"))
    primes = (2, 3, 5, 7, 11, 13, 17, 19)[:cert['k']]
    assert list(primes) == cert['primes']
    n, *ss = cert['target']
    witness = dict(cert['witness'])
    assert len(witness) == len(cert['witness']) and all(x > 0 for x in witness.values())
    lengths = [math.prod(pp) for t in range(len(primes)+1) for pp in itertools.combinations(primes, t)]
    assert set(witness) <= set(lengths) and 1 not in witness
    assert sum(d*x for d, x in witness.items()) == n
    assert [sum(d*x for d, x in witness.items() if d % p == 0) for p in primes] == ss
    lengths = sorted((d for d in lengths if all(d <= (s if d % p == 0 else n-s)
                      for p, s in zip(primes, ss, strict=True))), reverse=True)
    assert lengths == cert['columns']
    cols = [(d, *(d if d % p == 0 else 0 for p in primes)) for d in lengths]

    @lru_cache(None)
    def cut_info(start, direction):
        assert len(direction) == len(primes)+1 and all(type(v) is int for v in direction)
        weights = {dot(direction, v) for v in cols[start:]}
        assert min(weights) >= 0
        weights.discard(0)
        g = math.gcd(*weights)
        if not g:
            return 0, 0, None
        steps = sorted(t//g for t in weights)
        m = steps[0]
        if not 1 < m <= 256:
            return g, m, None
        # Dijkstra on residue classes with the original positive coin lengths.
        dist = [math.inf]*m
        dist[0] = 0
        queue = [(0, 0)]
        while queue:
            d, i = heapq.heappop(queue)
            if d != dist[i]:
                continue
            for step in steps:
                j = (i+step) % m
                if d+step < dist[j]:
                    dist[j] = d+step
                    heapq.heappush(queue, (d+step, j))
        return g, m, dist

    nodes, cuts_checked = 0, 0
    solutions = []

    def visit(j, residual, tree, prefix):
        nonlocal nodes, cuts_checked
        nodes += 1
        v = cols[j]
        nleft = residual[0]
        cap = min([nleft]+[b if c else nleft-b for b, c in zip(residual[1:], v[1:], strict=True)]) // lengths[j]
        domain = set(range(cap+1))
        if j+1 == len(cols):
            assert tree is None
            for x in domain:
                if all(t == x*c for c, t in zip(v, residual, strict=True)):
                    sol = dict((d, y) for d, y in zip(lengths, (*prefix, x), strict=True) if y)
                    solutions.append(sol)
            return
        cuts, branches = tree
        for direction in cuts:
            g, m, distances = cut_info(j+1, tuple(direction))
            alpha, beta = dot(direction, residual), dot(direction, v)
            passing = set()
            for x in domain:
                y = alpha-beta*x
                if y < 0 or (not g and y != 0):
                    continue
                if g and (y % g or (distances is not None and y//g < distances[(y//g) % m])):
                    continue
                passing.add(x)
            domain = passing
            cuts_checked += 1
        assert sorted(domain) == sorted(x for x, _ in branches)
        for x, child in branches:
            visit(j+1, tuple(t-x*c for t, c in zip(residual, v, strict=True)), child, (*prefix, x))

    visit(0, (n, *ss), cert['proof'], ())
    assert solutions == [witness]
    return dict(file=file.name, k=len(primes), n=n, s=ss, support=len(witness), proof_nodes=nodes,
                cuts_checked=cuts_checked, solutions=len(solutions))


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Use Python without -O to enable strict checks.")
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("files", type=Path, nargs="+")
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    out = json.dumps([check(f) for f in args.files], indent=2)
    if args.output:
        args.output.write_text(out+"\n", encoding="utf-8")
    print(out)
