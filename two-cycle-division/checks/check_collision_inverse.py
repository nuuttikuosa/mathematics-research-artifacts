"""Independent finite collision audit. Run: python -B check_collision_inverse.py"""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import product
from math import gcd, isqrt, lcm


def freeze(counts):
    return tuple(sorted((x, k) for x, k in counts.items() if k))


def divisors(n):
    return {v for d in range(1, isqrt(n) + 1) if n % d == 0 for v in (d, n // d)}


@lru_cache(None)
def literal_factor(cycle, t):
    """Trace the actual permutation (i,j)->(i+1,j+1) on all cycle*t states."""
    seen = bytearray(cycle * t)
    counts = Counter()
    for start in range(len(seen)):
        if seen[start]:
            continue
        pos = start
        ticks = 0
        while not seen[pos]:
            seen[pos] = 1
            i, j = divmod(pos, t)
            pos = ((i + 1) % cycle) * t + (j + 1) % t
            ticks += 1
        assert start == pos
        counts[ticks] += 1
    return freeze(counts)


def projected(cycle, quotient):
    result = Counter()
    for t, k in quotient:
        for ell, c in literal_factor(cycle, t):
            result[ell] += k * c
    return result


def raw_fibre(a, b, target):
    """All raw periods divide an observed period, since returning returns both states."""
    target = +target
    periods = sorted(set.union(*(divisors(ell) for ell in target)), reverse=True)
    columns = []
    for t in periods:
        col = Counter(dict(literal_factor(a, t))) + Counter(dict(literal_factor(b, t)))
        if all(target[ell] >= k for ell, k in col.items()):
            columns.append((t, col))
    solutions = set()

    def visit(i, left, chosen):
        if not left:
            solutions.add(tuple(sorted(chosen)))
            return
        if i == len(columns):
            return
        t, col = columns[i]
        remaining = left.copy()
        for k in range(1 + min(left[ell] // c for ell, c in col.items())):
            visit(i + 1, +remaining, chosen + ([(t, k)] if k else []))
            remaining.subtract(col)

    visit(0, target, [])
    return solutions


def supported_part(n, support):
    out = 1
    while (part := gcd(n, support)) > 1:
        out *= part
        n //= part
    return out


def transform(a, b, target):
    g = gcd(a, b)
    u, v, D = a // g, b // g, lcm(a, b)
    p0 = pow(u, -1, v) if v > 1 else 1
    P, Q = p0 + 2 * v, (u * p0 - 1) // v + 2 * u
    blocks = defaultdict(Counter)
    for ell, k in target.items():
        blocks[lcm(D, ell)][ell] = k
    rows, cols, records = Counter(), Counter(), {}
    for R, counts in blocks.items():
        if any(ell % a and ell % b for ell in counts):
            return None
        size = sum(ell * k for ell, k in counts.items())
        if size % (a + b):
            return None
        n = size // (a + b)
        rr = Counter({ell:k for ell,k in counts.items() if ell != R and ell % a == 0})
        cc = Counter({ell:k for ell,k in counts.items() if ell != R and ell % b == 0})
        row_mass = sum(ell * k for ell, k in rr.items())
        # Find the split by testing every possible shared quota.
        choices = [j for j in range(counts[R] + 1) if row_mass + R*j == a*n]
        if not choices:
            return None
        assert len(choices) == 1
        rr[R], cc[R] = choices[0], counts[R] - choices[0]
        pr, qr = supported_part(P, R // D), supported_part(Q, R // D)
        records[R] = (pr, qr)
        for ell, k in rr.items():
            rows[P*qr*ell] += pr*k
        for ell, k in cc.items():
            cols[Q*pr*ell] += qr*k
    return (P, Q, +rows, +cols, records)


def audit(a, b, target, stats, name=None):
    transform_data = transform(a, b, target)
    old = raw_fibre(a, b, target)
    if not transform_data:
        assert not old
        stats["rejected at margin split"] += 1
        return
    P, Q, rows, cols, records = transform_data
    a2, b2 = P*a, Q*b
    D = lcm(a, b)
    new = raw_fibre(a2, b2, rows + cols)
    stats["passed margin split"] += 1
    shared = sorted(set(rows) & set(cols))
    if shared:
        stats["with positive row-column collision"] += 1
        category = "empty" if not new else "unique" if len(new) == 1 else "multiple"
        stats["collision, " + category + " fibre"] += 1
        stats["output classes at a collision (sum per case)"] += len(new)
        if len(records) > 1:
            stats["collision with multiple blocks"] += 1
        for ell in shared:
            assert ell % (P*Q) == 0 and records[ell // (P*Q)] == (P, Q)
    expected = set()
    for x in old:
        ys = []
        for t, k in x:
            pr, qr = records[lcm(D,t)]
            ys.append((pr*qr*t, k))
        expected.add(tuple(sorted(ys)))
    assert len(old) == len(expected) and expected == new, (a,b,target,old,new)
    recovered = set()
    for yset in new:
        assert projected(a2, yset) == rows and projected(b2, yset) == cols
        xs = []
        for y, k in yset:
            image_R = lcm(a2,b2,y)
            assert image_R % (P*Q) == 0
            R = image_R // (P*Q)
            pr, qr = records[R]
            assert y % (pr*qr) == 0
            t = y // (pr*qr)
            assert lcm(D,t) == R
            xs.append((t,k))
        recovered.add(tuple(sorted(xs)))
    assert recovered == old
    if name:
        print(f"{name}: old={sorted(old)}; new={sorted(new)}; collision={shared}")
        print(f"  axes={(a2,b2)}; row={freeze(rows)}; col={freeze(cols)}")


def main():
    stats = Counter()
    # Complete boxes of observed multiplicities. The lifts by 35 make the
    # new row and column meet at length 7350 when both margins are positive.
    for lift in (1, 35):
        for counts in product(range(8), repeat=3):
            target = +Counter(dict(zip((2*lift,3*lift,6*lift),counts)))
            if target:
                audit(3, 2, target, stats)
    plus = Counter({70:6,105:6,210:5})
    minus = Counter({105:1,210:2})
    core = Counter({2:6,3:6,6:5})
    for name,a,b,target in [
        ("two colliding classes",3,2,plus),
        ("passed split and colliding, no solution",3,2,minus),
        ("multiblock, four classes",3,2,plus+core),
        ("multiblock, no solution",3,2,minus+core),
        ("equal axes with two classes",5,5,Counter({30:10})),
        ("unit axes",1,1,Counter({6:2})),
        ("gcd and partial scales",6,6,Counter({12:12,18:12})),
        ("prescribed negative control",2,3,Counter({3:1,6:2})),
        ("length divisible by both axes, zero column share",1,1,Counter({2:2})),
        ("full prime power at positive collision",14,21,Counter({420:35})),
    ]:
        audit(a,b,target,stats,name)
    for key, value in sorted(stats.items()):
        print(f"{key}: {value}")
    for key in ("collision, empty fibre","collision, multiple fibre","collision, unique fibre",
                "collision with multiple blocks"):
        assert stats[key] > 0
    print("PASS")


if __name__ == "__main__":
    main()
