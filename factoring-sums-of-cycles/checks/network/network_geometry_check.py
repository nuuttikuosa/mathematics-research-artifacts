"""Independent generic-permutation audit of the article's network geometry.

No existing project code is imported. Tests all one-digit permutations for
u=3,4,5, and deterministic structured/random multi-digit permutations.
All signals are active: vertex separation and line occupancy are stronger
checks than on an arbitrary subset of occurrences. Final Boolean blocks
and SAT tests are outside this script's scope.
"""
from itertools import permutations, product
from pathlib import Path
import hashlib
import json
import random

ROOT = Path(__file__).resolve().parent


def reflected(values, dimension):
    if not dimension:
        return [()]
    tails = reflected(values, dimension - 1)
    return [(v,) + tail for i, v in enumerate(values)
            for tail in (tails if i % 2 == 0 else tails[::-1])]


def decompose(mapping, u, g):
    words = sorted(mapping)
    assert set(mapping.values()) == set(words)
    if g == 1:
        return [(0, dict(mapping))]
    suffixes = sorted(set(w[1:] for w in words))
    remaining = set(words)
    color = {}
    for c in range(u):
        adj = {s: [w for w in words if w in remaining and w[1:] == s]
               for s in suffixes}
        matched = {}

        def augment(left, seen):
            for w in adj[left]:
                right = mapping[w][1:]
                if right in seen:
                    continue
                seen.add(right)
                if right not in matched or augment(matched[right][1:], seen):
                    matched[right] = w
                    return True
            return False

        assert all(augment(left, set()) for left in suffixes)
        assert len(matched) == len(suffixes)
        for w in matched.values():
            color[w] = c
            remaining.remove(w)
    assert not remaining
    first = {w: (color[w],) + w[1:] for w in words}
    last = {(color[w],) + mapping[w][1:]: mapping[w] for w in words}
    branches = {}
    for c in range(u):
        small = {w[1:]: mapping[w][1:] for w in words if color[w] == c}
        branches[c] = decompose(small, u, g - 1)
    middle = []
    for stage in range(2 * g - 3):
        index = branches[0][stage][0]
        assert all(branches[c][stage][0] == index for c in range(u))
        joined = {(c,) + w: (c,) + target
                  for c in range(u) for w, target in branches[c][stage][1].items()}
        middle.append((index + 1, joined))
    return [(0, first)] + middle + [(0, last)]


def audit(mapping, u, g):
    words = sorted(mapping)
    stages = decompose(mapping, u, g)
    assert [i for i, _ in stages] == list(range(g)) + list(range(g - 2, -1, -1))
    composed = {w: w for w in words}
    for index, stage in stages:
        assert set(stage) == set(stage.values()) == set(words)
        assert all(all(a[j] == b[j] for j in range(g) if j != index)
                   for a, b in stage.items())
        composed = {w: stage[v] for w, v in composed.items()}
    assert composed == mapping
    h = 3 * len(stages)
    path = [(0, 0)]
    for i in range(u - 1):
        path.extend([(i, i + 1), (i + 1, i + 1)])
    k = 1
    while len(path) ** k < h:
        k += 1
    tuples = reflected(path, k)[:h]
    times = [sum(t, ()) for t in tuples]
    assert len(set(times)) == h
    assert all(sum(a != b for a, b in zip(x, y)) == 1
               for x, y in zip(times, times[1:]))
    for axis in range(2 * k):
        rows = {}
        for t in times:
            key = t[:axis] + t[axis + 1:]
            rows.setdefault(key, set()).add(t)
        assert all(len(points) <= 2 for points in rows.values())

    def point(word, intermediate, time, lane):
        return word + ((sum(word) + intermediate) % u,) + times[time] + (lane,)

    def color(p):
        epsilon = (p[g] - sum(p[:g])) % u
        assert epsilon in (0, 1)
        t = p[g + 1:-1]
        return (epsilon + sum(t[i + 1] - t[i] for i in range(0, len(t), 2))
                + int(p[-1] != 0)) % 2

    owner = {}
    lines = {}
    edge_count = 0
    wrapping_moves = 0
    fixed_moves = 0
    for original in words:
        route = []
        current = original
        for stage_number, (axis, stage) in enumerate(stages):
            target = stage[current]
            buffer = current[axis] + 1
            t = 3 * stage_number
            chunk = [point(current, 0, t, 0), point(current, 0, t, buffer),
                     point(current, 0, t + 1, buffer)]
            moves = (target[axis] - current[axis]) % u
            fixed_moves += int(moves == 0)
            word = current
            for _ in range(moves):
                chunk.append(point(word, 1, t + 1, buffer))
                new = list(word)
                wrapping_moves += int(new[axis] == u - 1)
                new[axis] = (new[axis] + 1) % u
                word = tuple(new)
                chunk.append(point(word, 0, t + 1, buffer))
            assert word == target
            chunk.extend([point(target, 0, t + 2, buffer), point(target, 0, t + 2, 0)])
            route.extend(chunk)
            current = target
        assert len(route) == len(set(route)), 'A signal route revisits a vertex'
        assert color(route[0]) == 0
        for p in route:
            assert p not in owner, ('Unintended meeting of signals', owner.get(p), original)
            owner[p] = original
        for p, q in zip(route, route[1:]):
            changed = [i for i in range(len(p)) if p[i] != q[i]]
            assert len(changed) == 1
            i = changed[0]
            assert color(p) != color(q)
            key = (i, p[:i] + p[i + 1:])
            endpoints = frozenset((p, q))
            assert key not in lines or lines[key] == endpoints
            lines[key] = endpoints
            edge_count += 1
    # Reconstruct occupied lines from the complete active set. Labels on
    # routes are not used by this occupancy comparison.
    occupied = {key: set() for key in lines}
    dimension = g + 1 + 2 * k + 1
    for p in owner:
        for i in range(dimension):
            key = (i, p[:i] + p[i + 1:])
            if key in occupied:
                occupied[key].add(p)
    assert all(occupied[key] == set(endpoints) for key, endpoints in lines.items())
    return dict(u=u, g=g, signals=len(words), stages=len(stages), time_dimension=2*k,
                active_vertices=len(owner), edges=edge_count, marked_lines=len(lines),
                modular_wrap_moves=wrapping_moves, fixed_digit_moves=fixed_moves)


def main():
    rng = random.Random(261126)
    results = []
    for u in (3, 4, 5):
        words = list(product(range(u), repeat=1))
        for target in permutations(words):
            results.append(audit(dict(zip(words, target)), u, 1))
    for u, g, random_count in [(3, 2, 30), (3, 3, 12), (4, 2, 10), (3, 4, 6), (3, 5, 2)]:
        words = list(product(range(u), repeat=g))
        targets = [words, words[::-1], [tuple((x + 1) % u for x in w) for w in words]]
        for _ in range(random_count):
            target = words[:]
            rng.shuffle(target)
            targets.append(target)
        for target in targets:
            results.append(audit(dict(zip(words, target)), u, g))
    summary = dict(cases=len(results), exhaustive_one_digit_permutations=150,
                   additional_structured_or_seeded_cases=len(results)-150,
                   active_vertices=sum(r['active_vertices'] for r in results),
                   checked_edges=sum(r['edges'] for r in results),
                   modular_wrap_moves=sum(r['modular_wrap_moves'] for r in results),
                   fixed_digit_moves=sum(r['fixed_digit_moves'] for r in results),
                   maximum_word_count=max(r['signals'] for r in results),
                   maximum_g=max(r['g'] for r in results), all_assertions_passed=True,
                   source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   scope='Generic full permutations before final literal blocks and Boolean tests',
                   cases_detail=results)
    (ROOT/'network_geometry_results.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k != 'cases_detail'}, indent=2))


if __name__ == '__main__':
    main()
