"""Checks for cyclically decorated components. Run with Python 3.9+.

Uses explicit successor tables for systems of cycles with pendant leaves.
Comparison retains the full cyclic list, with no primitive-word reduction.
The small quotient oracle tests components by literal multiplication.
"""
from collections import Counter
from itertools import product


def rotate(word, start):
    return word[start:] + word[:start]


def canonical(word):
    return min(rotate(word, i) for i in range(len(word)))


def make_table(words):
    result = []
    for word in words:
        offset = len(result)
        result.extend(offset + (i + 1) % len(word) for i in range(len(word)))
        for i, leaf_count in enumerate(word):
            result.extend([offset + i] * leaf_count)
    return result


def multiply(axes, table):
    divisor = make_table([tuple([0] * a) for a in axes])
    return [a * len(table) + b for a in divisor for b in table]


def classes(table):
    # For height <= 1, periodic vertices are exactly the image of the map.
    roots = set(table)
    incoming = Counter(table)
    remaining = roots.copy()
    result = Counter()
    while remaining:
        v = next(iter(remaining))
        start, word = v, []
        while v in remaining:
            remaining.remove(v)
            word.append(incoming[v] - 1)
            v = table[v]
        assert v == start
        result[canonical(tuple(word))] += 1
    return result


def fibre(axes, target):
    # A candidate component with t roots must appear repeated on at least
    # one observed cycle: projection on the component advances one root
    # per step and multiplication by a permutation preserves its leaves.
    # Thus every possible whole word occurs among these divisors and cuts.
    # Only literal products determine whether it actually fits. Any factor
    # containing a deeper tree would give a deeper tree in the target.
    candidates = set()
    for full in target:
        for t in range(1, len(full) + 1):
            if len(full) % t == 0:
                for cut in range(len(full)):
                    candidates.add(canonical(rotate(full, cut)[:t]))
    columns = []
    for w in sorted(candidates):
        actual = classes(multiply(axes, make_table([w])))
        if all(n <= target[key] for key, n in actual.items()):
            columns.append((w, actual))

    def finish(i, left):
        if i == len(columns):
            return {()} if not any(left.values()) else set()
        w, contribution = columns[i]
        limit = min(left[key] // n for key, n in contribution.items())
        out = set()
        rest = left.copy()
        for count in range(limit + 1):
            term = ((w, count),) if count else ()
            out.update(term + tail for tail in finish(i + 1, rest))
            rest.subtract(contribution)
        return out

    return finish(0, target)


def main():
    # Two copies of a chiral word can be cut independently. In a synchronous
    # product with two fixed points both components must have the same class.
    w = (0, 1, 3, 0, 0, 0)
    assert canonical(w) != canonical(w[::-1])
    image = classes(make_table([w * 6] * 5))
    for i, j in product(range(6), repeat=2):
        old = classes(make_table([rotate(w, i), rotate(w, j)]))
        assert fibre((1, 1), old) == {((canonical(w), 1),)}
        # Different choices of cut on the five output components.
        rotated = [rotate(w * 6, (i + n * j) % 36) for n in range(5)]
        assert classes(make_table(rotated)) == image
    assert fibre((3, 2), image) == {((canonical(w * 6), 1),)}
    opposite = classes(make_table([w, w[::-1]]))
    assert not fibre((1, 1), opposite)
    print("36 independent pairs of cuts: one source class and one image class; reflection gives no source quotient.")

    # A nonunique example. Both length 3 twice and length 6 once work,
    # for old axes (6,6) and new axes (18,12). There is no extra relative
    # phase in either class.
    v = (0, 1, 2)
    old = classes(make_table([rotate(v * 2, i % 6) for i in range(12)]))
    new = classes(make_table(
        [rotate(v * 6, i % 18) for i in range(6)]
        + [rotate(v * 4, (i + 1) % 12) for i in range(6)]
    ))
    expected = {((canonical(v), 2),), ((canonical(v * 2), 1),)}
    assert fibre((6, 6), old) == fibre((18, 12), new) == expected
    print("Nonunique three-letter word: exactly two classes in each fibre.")

    # Split passes for C3 + 2 C6, for (2,3). The row is C6, the column
    # is C3 + C6. P=8,Q=5,E=1; no quotient exists on either side.
    old = classes(make_table([v, v * 2, rotate(v * 2, 1)]))
    assert (2 * (3 + 6 + 6) // (2 + 3)) == 6
    new_words = [v * 16, v * 5, v * 10]
    assert not fibre((2, 3), old) and not fibre((16, 15), classes(make_table(new_words)))
    print("Late negative: original and transformed target have zero classes.")

    # Large single star in a primitive six-letter word. The normal form for
    # (1,1) has P=3,Q=2,R=6,E=6. Old two copies and new five copies come
    # from actual successor products, not just the numerical size formula.
    for leaves in (1, 257, 8192):
        word = (leaves, 0, 0, 0, 0, 0)
        old_product = multiply((1, 1), make_table([word]))
        new_product = multiply((3, 2), make_table([rotate(word * 6, 11)]))
        assert classes(old_product) == Counter({canonical(word): 2})
        assert classes(new_product) == Counter({canonical(word * 6): 5})
        assert len(new_product) == 15 * len(old_product)
        total_in = 2 + len(old_product)
        total_out = 5 + len(new_product)
        assert total_out <= 3 * 2 * 3 * total_in
        leaf_word = "leaf" if leaves == 1 else "leaves"
        print(f"{leaves} {leaf_word}: B={len(old_product)}, B'={len(new_product)}, "
              f"total input={total_in}, total output={total_out}")
    print("PASS")


if __name__ == "__main__":
    main()
