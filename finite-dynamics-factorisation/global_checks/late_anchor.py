"""Follow one pair of guides through the 420-state system."""
from collections import Counter
from functools import cache

import sympy as sp
import fdds
from search import search
from residues import labels, projected_exponents
from tree_arithmetic import Trees
from weighted_split import anchor_recovery, dot, local_pairs_twogrades, planar_slice


def table(root_degree, child_degrees):
    result = [0] * root_degree
    for child, degree in enumerate(child_degrees, 1):
        result.extend([child] * degree)
    return result


def check_seed(trees, observed, c0, d0):
    anchor = trees.product(c0, d0)
    data = labels(trees, observed)
    differences = {tree: tuple(a - b for a, b in zip(vector, data[anchor]))
                   for tree, vector in data.items()}
    keep = sp.Matrix(list(differences.values())).rref()[1]
    positions = {tree: tuple(vector[i] for i in keep) for tree, vector in differences.items()}
    spectrum = {positions[tree]: n for tree, n in observed.items()}
    dictionary = {}
    for tree, point in positions.items():
        quotient = trees.divide(tree, d0)
        if quotient is not None:
            dictionary[point] = quotient
    c, w, scalar, _, _, bound = planar_slice(spectrum)
    assert (c, w, scalar, bound) == ((1, 4), (0, 1), 2, 1)
    low = min(dot(v, w) for v in spectrum)
    weighted = {v: n * scalar**(dot(v, w) - low) for v, n in spectrum.items()}
    factors, basis, pivots = projected_exponents(weighted, dict(enumerate(c)))
    assert list(pivots) == list(range(4)) and basis == sp.eye(4)
    # The chosen pivot is (1,0,0,0). Since the basis is the identity, its
    # projected factor is exactly the first rational irreducible.
    first, multiplicity = factors[0]
    assert first.as_expr() == first.gens[0] + 1 and multiplicity == 1
    raw_axes = anchor_recovery(dictionary, c, w, scalar, first, positions[anchor])
    forest = trees.node(observed.elements())
    positive = []
    for raw in raw_axes:
        primitive = Counter({dictionary[key]: count for key, count in raw})
        for content in range(1, observed.total() // primitive.total() + 1):
            if observed.total() % (content * primitive.total()):
                continue
            left = trees.node(tree for tree, count in primitive.items()
                              for _ in range(content * count))
            right = trees.divide(forest, left)
            if right is None:
                continue
            assert trees.product(left, right) == forest
            assert c0 in trees.children[left] and d0 not in trees.children[right]
            degrees = lambda tree: tuple(sorted(len(trees.children[t]) for t in trees.children[tree]))
            positive.append((degrees(left), degrees(right)))
    assert positive == [((1, 2), (1, 3, 4, 6, 12, 18, 24, 72))]


def check_420():
    child_degrees = [1, 2, 3, 4, 6, 6, 8, 12, 12, 18, 24, 24, 36, 48, 72]
    target = table(144, child_degrees)
    x, y = table(12, (1, 2, 6)), table(12, (1, 3, 4))
    assert (len(target), len(x), len(y)) == (420, 21, 20)
    assert fdds.canon(fdds.fproduct(x, y)) == fdds.canon(target)
    endpoints = target.index(8), target.index(11)
    assert endpoints == (174, 216)
    degrees = Counter(target)
    assert (degrees[8], degrees[11]) == (12, 24)
    budget = (420 // degrees[8]) * (420 // degrees[11])
    assert budget == 595
    trees = Trees()

    @cache
    def convert(tree):
        return trees.node(convert(child) for child in fdds.KIDS[tree])

    @cache
    def convert_back(tree):
        return fdds.node(convert_back(child) for child in trees.children[tree])

    def degrees_of(axis):
        return tuple(sorted(len(trees.children[tree]) for tree in axis.elements()))

    levels = []

    def selected_split(level, row, left_anchor, right_anchor):
        observed = Counter({convert(tree): count for tree, count in row.items()})
        c0, d0 = convert(left_anchor), convert(right_anchor)
        levels.append(level)
        if level == 1:
            assert len(trees.children[c0]) == 1 and len(trees.children[d0]) == 2
            assert observed.total() == 16 <= budget
            assert degrees_of(observed) == tuple(sorted([144, *child_degrees]))
            check_seed(trees, observed, c0, d0)
        result = []
        for left, right in local_pairs_twogrades(trees, observed, c0, d0, Counter()):
            if level == 0 and (left.total(), right.total()) != (1, 2):
                continue
            if level == 1 and (degrees_of(left), degrees_of(right)) != ((1,3,4,12), (1,2,6,12)):
                continue
            result.append((Counter({convert_back(t): n for t, n in left.items()}),
                           Counter({convert_back(t): n for t, n in right.items()})))
        return result

    found, _ = search(target, endpoint_pair=endpoints, selected_split=selected_split)
    assert levels == [0, 1, 2] and found == {(fdds.canon(x), fdds.canon(y))}
    print("420 states: guides (174,216), local mass 16, budget 595, rank 4, tau 2.")
    print("Pivot (1,0,0,0): exact positive quotient misses the right anchor S2.")
    print("The continuing split on these guides gives the 21-by-20 witness.")
