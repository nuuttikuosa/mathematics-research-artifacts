"""Decode anchored star forests at a prescribed positive specialization.

This test variant chooses the same two grades as the production path.
The ranges use the actual multiplicities at the supplied specialization.
"""
from collections import Counter
from itertools import product as cartesian

import sympy as sp
from residues import labels, projected_exponents
from weighted_split import anchor_recovery, dot, planar_slice

T = sp.Symbol("t")


def splits_with_tau(trees, target, c0, d0, force_tau):
    """Manuscript's anchored enumeration with a prescribed tau and the actual ranges."""
    anchor_tree = trees.product(c0, d0)
    allowed = {tree: trees.divide(tree, d0) for tree in target}
    allowed = {tree: child for tree, child in allowed.items() if child is not None}
    data = labels(trees, target)
    origin = data[anchor_tree]
    diffs = {tree: tuple(x-y for x, y in zip(v, origin)) for tree, v in data.items()}
    keep = sp.Matrix(list(diffs.values())).rref()[1]
    positions = {tree: tuple(v[i] for i in keep) for tree, v in diffs.items()}
    spectrum = {positions[tree]: n for tree, n in target.items()}
    dictionary = {positions[tree]: child for tree, child in allowed.items()}
    c, w, tau_alg, _, _, mu = planar_slice(spectrum)
    tau = force_tau
    lo = min(dot(v, w) for v in spectrum)
    weighted = {v: n*tau**(dot(v, w)-lo) for v, n in spectrum.items()}
    factors, basis, pivots = projected_exponents(weighted, dict(enumerate(c)))
    mass = sum(target.values())
    limits = [min(m, mass-1) for _, m in factors]   # actual multiplicities of this slice
    recovered = set()
    for seed in cartesian(*(range(limits[i]+1) for i in pivots)):
        values = [sum((a*basis[k, j] for k, a in enumerate(seed)), sp.Integer(0))
                  for j in range(len(factors))]
        if any(not a.is_Integer or not 0 <= a <= lim for a, lim in zip(values, limits)):
            continue
        q = sp.Poly(1, T, domain=sp.QQ)
        for (factor, _), e in zip(factors, values):
            if e:
                q *= sp.Poly.from_dict({m: cf for m, cf in factor.terms()}, T, domain=sp.QQ)**int(e)
        recovered.update(anchor_recovery(dictionary, c, w, tau, q, positions[anchor_tree]))
    out = set()
    forest = trees.node(target.elements())
    for raw in recovered:
        total = sum(n for _, n in raw)
        for k in range(1, mass//total+1):
            if mass % (k*total):
                continue
            left = Counter({dictionary[v]: n*k for v, n in raw})
            right = trees.divide(forest, trees.node(left.elements()))
            if right is not None and d0 in trees.children[right] and \
               trees.product(trees.node(left.elements()), right) == forest:
                out.add((tuple(sorted(left.items())), tuple(sorted(Counter(trees.children[right]).items()))))
    excess = max(m for _, m in factors) > mu
    return out, excess, tau_alg, mu, c, w, positions
