"""Compare ordinary radius cuts and a variant that removes short side branches.

Run: python -B check_anchor_cuts.py
Python 3.10+; standard library only. This tests the physical common-anchor
identities. It does not compute the label, grade, or full factorization search.
"""
from collections import Counter
from functools import cache


def node(children):
    return tuple(sorted(children))


@cache
def height(tree):
    return 1 + max((height(child) for child in tree), default=-1)


@cache
def ordinary_cut(tree, radius):
    if radius == 0:
        return ()
    return node(ordinary_cut(child, radius - 1) for child in tree)


@cache
def pruned_cut(tree, radius):
    if radius == 0:
        return ()
    return node(pruned_cut(child, radius - 1) for child in tree
                if height(child) >= radius - 1)


@cache
def multiply(left, right):
    return node(multiply(a, b) for a in left for b in right)


def catalog(limit):
    """Adjoin the final leaf at every vertex and deduplicate canonical trees."""
    @cache
    def add_leaf(tree):
        result = {node((*tree, ()))}
        for i, child in enumerate(tree):
            for extended in add_leaf(child):
                result.add(node((*tree[:i], extended, *tree[i + 1:])))
        return result

    levels = [{()}]
    for _ in range(1, limit):
        levels.append({new for old in levels[-1] for new in add_leaf(old)})
    assert [len(level) for level in levels] == [1, 1, 2, 4, 9, 20, 48, 115][:limit]
    return [tree for level in levels[1:] for tree in level]


def anchor_path(tree):
    path = [tree]
    while path[-1]:
        path.append(max(path[-1], key=height))
    return path


@cache
def dictionary(tree, radius, cut):
    found = {cut(tree, radius)} if height(tree) >= radius else set()
    for child in tree:
        found.update(dictionary(child, radius, cut))
    return frozenset(found)


def slots(tree, depth, horizon):
    # Each slot contains all active children of one parent at a common radius.
    for radius in range(horizon):
        active = [child for child in tree if height(child) >= radius]
        yield depth + 1, radius, active
        for child in active:
            yield from slots(child, depth + 1, radius)


def audit(limit, cut):
    trees = catalog(limit)
    stats = Counter()
    for left in trees:
        h = height(left)
        groups = list(slots(left, 0, h))
        for right in trees:
            if height(right) != h:
                continue
            stats["pairs"] += 1
            target = multiply(left, right)
            path = anchor_path(right)
            for depth, radius, children in groups:
                assert children and depth + radius <= h
                candidates = [cut(child, radius) for child in children]
                partner = cut(path[depth], radius)
                assert height(partner) == radius
                observed = dictionary(target, radius, cut)
                images = {u: multiply(u, partner) for u in candidates}
                for child, u in zip(children, candidates):
                    assert images[u] == cut(multiply(child, path[depth]), radius)
                assert set(images.values()) <= observed
                assert len(set(images.values())) == len(images)
                outside = sum(u not in observed for u in candidates)
                stats["groups"] += 1
                stats["outside_occurrences"] += outside
                stats["outside_groups"] += int(outside > 0)
                stats["mixed_groups"] += int(len(images) > 1)
                stats["mixed_radius2"] += int(len(images) > 1 and radius == 2)
                stats["unequal_ordered_pairs"] += sum(a != b for a in candidates for b in candidates)
    return stats


def witness():
    leaf = ()
    path_one = (leaf,)
    path_two = (path_one,)
    early_leaf = node([leaf, path_one])
    left = node([path_two, early_leaf])
    right = (path_two,)
    assert height(left) == height(right) == 3 and multiply(left, right) == left
    assert ordinary_cut(path_two, 2) != ordinary_cut(early_leaf, 2)
    assert pruned_cut(path_two, 2) == pruned_cut(early_leaf, 2)
    print("Eight-vertex witness: ordinary cut separates the siblings; pruned cut merges them.")


if __name__ == "__main__":
    keys = ("pairs", "groups", "outside_occurrences", "outside_groups",
            "mixed_groups", "mixed_radius2", "unequal_ordered_pairs")
    witness()
    print("limit cut", *keys)
    for limit in (7, 8):
        for cut in (pruned_cut, ordinary_cut):
            counts = audit(limit, cut)
            print(limit, cut.__name__, *(counts[key] for key in keys))
