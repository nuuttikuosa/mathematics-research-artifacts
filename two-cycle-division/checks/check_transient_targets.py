"""Exact small-fibre check with labelled successor tables as reference.

Run: python -B check_transient_targets.py
The sibling check_marker_normal_form.py supplies a raw permutation quotient
enumerator which uses product multiplication and integer backtracking.
"""
from collections import Counter, defaultdict, deque
from itertools import product
from math import gcd, lcm
from pathlib import Path
from random import Random
import sys

# The runner uses Python's isolated mode; import the sibling reference only.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_marker_normal_form import freeze, marked_target, solutions


def primitive(w):
    for k in range(1, len(w) + 1):
        if len(w) % k == 0 and w == w[:k] * (len(w) // k):
            part = w[:k]
            return min(part[i:] + part[:i] for i in range(k))


def profile(successor):
    n = len(successor)
    pred = [[] for _ in range(n)]
    for v, f in enumerate(successor):
        pred[f].append(v)
    indegree = [len(c) for c in pred]
    queue = deque(v for v in range(n) if not indegree[v])
    while queue:
        f = successor[queue.popleft()]
        indegree[f] -= 1
        if not indegree[f]:
            queue.append(f)

    def tree(v):
        return "(" + "".join(sorted(tree(c) for c in pred[v] if not indegree[c])) + ")"

    visited, result = set(), Counter()
    for v, d in enumerate(indegree):
        if d and v not in visited:
            w, c = [], v
            while c not in visited:
                visited.add(c)
                w.append(tree(c))
                c = successor[c]
            result[primitive(tuple(w)), len(w)] += 1
    return freeze(result)


def table_product(axes, fx):
    ax = []
    for a in axes:
        base = len(ax)
        ax.extend(base + (v + 1) % a for v in range(a))
    return [next_a * len(fx) + next_x for next_a in ax for next_x in fx]


def projected(axes, target, stats):
    a, b = axes
    D = lcm(a, b)
    groups = defaultdict(dict)
    for (w, r), qty in target:
        groups[w, lcm(D, r)][r] = qty
    all_profiles = {()}
    if len({w for w, R in groups}) > 1:
        stats["several words"] += 1
    if any(len(w) > 1 for w, R in groups):
        stats["primitive period >1"] += 1
    if any(gcd(len(w), a) > 1 and gcd(len(w), b) > 1 for w, R in groups):
        stats["both scalings > 1"] += 1
    for (w, R), obs in groups.items():
        e = len(w)
        da, db = gcd(a, e), gcd(b, e)
        A, B = a // da, b // db
        assert (a // gcd(a, b)) % (A // gcd(A, B)) == 0
        assert (b // gcd(a, b)) % (B // gcd(A, B)) == 0
        assert lcm(A, B) == D // gcd(D, e)
        if any(r % a and r % b for r in obs):
            return set()
        nr, rem = divmod(sum(r * qty for r, qty in obs.items()), a + b)
        if rem or nr % e:
            return set()
        rows = {r: qty for r, qty in obs.items() if r != R and r % a == 0}
        cols = {r: qty for r, qty in obs.items() if r != R and r % b == 0}
        m, rem = divmod(a * nr - sum(r * qty for r, qty in rows.items()), R)
        if rem or not 0 <= m <= obs.get(R, 0):
            return set()
        if m:
            rows[R] = m
        if obs.get(R, 0) - m:
            cols[R] = obs[R] - m
        if any(q % da for q in rows.values()) or any(q % db for q in cols.values()):
            return set()
        new_target = Counter()
        for block, factor in [(rows, da), (cols, db)]:
            for r, qty in block.items():
                assert r % e == 0 and lcm(lcm(A, B), r // e) == R // e
                new_target[r // e] += qty // factor
        normalized = solutions((A, B), freeze(new_target))
        if not normalized:
            stats["negative after arithmetic checks"] += 1
            return set()
        if set(rows) & set(cols) and e > 1:
            stats["shared period with nontrivial word"] += 1
        new_all = set()
        for old, sol in product(all_profiles, normalized):
            terms = Counter(dict(old))
            for u, qty in sol:
                terms[w, e * u] += qty
            new_all.add(freeze(terms))
        all_profiles = new_all
    return all_profiles


def number_of_states(target):
    return sum(q * (r // len(w)) * sum(s.count("(") for s in w) for (w, r), q in target)


def expand(target):
    table = []
    for (w, r), copies in target:
        for _ in range(copies):
            base = len(table)
            table.extend(base + (i + 1) % r for i in range(r))
            for i in range(r):
                path = []
                for token in w[i % len(w)]:
                    if token == ")":
                        path.pop()
                    elif path:
                        v = len(table)
                        table.append(path[-1])
                        path.append(v)
                    else:
                        path.append(base + i)
    assert len(table) == number_of_states(target)
    return table


def by_word(target):
    grouped = defaultdict(dict)
    for (w, r), qty in target:
        grouped[w][r] = qty
    return grouped


def direct_decorated(axes, target):
    # Raw oracle for the image. It uses no block split, normal form, gcd
    # rescaling of marginals or input witness. One allowed cycle is a
    # repetition of the given primitive tree word.
    answer = {()}
    for w, obs in by_word(target).items():
        factor = [
            tuple(((w, t), qty) for t, qty in s)
            for s in solutions(axes, freeze(obs))
            if all(t % len(w) == 0 for t, qty in s)
        ]
        answer = {tuple(sorted(x + y)) for x, y in product(answer, factor)}
    return answer


def check_decorated_normal(axes, target, actual, stats):
    D = lcm(*axes)
    maps, output = {}, Counter()
    for w, obs in by_word(target).items():
        new_axes, transformed, blocks, factors, _ = marked_target(axes, freeze(obs), False)
        if not blocks:
            assert not actual
            stats["decorated normal early negatives"] += 1
            return
        for r1, mult in transformed:
            assert r1 % len(w) == 0
            output[w, r1] += mult
        for R, (PF, QF) in blocks.items():
            maps[w, R] = PF * QF
    new_target = freeze(output)
    new_solutions = direct_decorated(new_axes, new_target)
    image = set()
    for original in actual:
        terms = Counter()
        for (w, t), qty in original:
            terms[w, t * maps[w, lcm(D, t)]] += qty
        image.add(freeze(terms))
    assert len(image) == len(actual) and image == new_solutions
    if not image:
        stats["decorated normal late negatives"] += 1
    for x in image:
        # Construct a literal successor table for each predicted transformed
        # witness and multiply it by the actual transformed permutation.
        assert profile(table_product(new_axes, expand(x))) == new_target
        stats["literal transformed products"] += 1
        stats["largest literal transformed product"] = max(
            stats["largest literal transformed product"], number_of_states(new_target))


def main():
    rng = Random(70538101)
    all_reps = {}
    for n in range(1, 7):
        for table in product(range(n), repeat=n):
            code = profile(table)
            if code not in all_reps:
                all_reps[code] = table
    print("Different source classes on 1--6 states:", len(all_reps))
    stats = Counter()
    for axes in [(1, 1), (1, 3), (2, 2), (2, 3), (3, 2), (2, 4),
                 (4, 6), (6, 9), (4, 4), (6, 10), (8, 12)]:
        answers = defaultdict(set)
        for source, table in all_reps.items():
            target = profile(table_product(axes, table))
            answers[target].add(source)
        tests = set(answers)
        for target in sorted(answers):
            # Removing one component or changing the length of its cycle can
            # preserve mass/integrality in some groups and destroy it in others.
            chosen = rng.choice([x for x, qty in target])
            w, r = chosen
            for delta in (0, -len(w), len(w)):
                obs = Counter(dict(target))
                obs[chosen] -= 1
                if delta and r + delta > 0:
                    obs[w, r + delta] += 1
                changed = freeze(obs)
                if changed and number_of_states(changed) <= 6 * sum(axes):
                    tests.add(changed)
        for target in sorted(tests):
            predicted = projected(axes, target, stats)
            actual = answers.get(target, set())
            assert predicted == actual, (axes, target, predicted, actual)
            check_decorated_normal(axes, target, actual, stats)
            stats["instances"] += 1
            stats["positive" if actual else "negative"] += 1
            stats["source quotient classes"] += len(actual)
            stats["nonunique"] += len(actual) > 1
            stats["maximum fibre"] = max(stats["maximum fibre"], len(actual))
    # Hand-checked late negative. Its normalized permutation target is
    # C_2 + 2 C_3 + 2 C_6 for C_2+C_3. There can be no quotient, since
    # x_1+2x_2=1 forces x_1=1, contrary to x_1+3x_3=2. Here each x_i
    # is a multiplicity, and periods other than 1,2,3 cannot project to
    # 2 or 3.
    word = primitive(("()", "(())"))
    extra = freeze({(word, 4): 2, (word, 6): 4, (word, 12): 4})
    assert not projected((4, 6), extra, stats)
    check_decorated_normal((4, 6), extra, set(), stats)
    stats["additional hand-certified negatives"] += 1
    print(dict(sorted(stats.items())))
    print("PASS")


if __name__ == "__main__":
    main()
