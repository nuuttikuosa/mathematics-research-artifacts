"""Run the catalogue, transition-table budget and late-rejection experiments."""
from collections import Counter, defaultdict
from functools import cache
from pathlib import Path
import sys

import fdds
from guide_budget import arbitrary_endpoints, true_guides
from late_family import check as check_family
from search import search

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "local_checks"))
from tree_arithmetic import Trees
from weighted_split import local_pairs_twogrades
from late_anchor import check_420


def compare_catalogue():
    # Each factor of a nontrivial product with at most 12 states has at most
    # 6 states. Enumerating connected factor types therefore exhausts the
    # decomposable targets in this size range, including disconnected targets.
    types = {n: list(fdds.connected_types(n).values()) for n in range(2, 7)}
    oracle = defaultdict(set)
    representatives = {}
    for a in types:
        for b in types:
            if a * b > 12:
                continue
            for x in types[a]:
                for y in types[b]:
                    target = fdds.fproduct(x, y)
                    key = fdds.canon(target)
                    representatives.setdefault(key, target)
                    oracle[key].add((fdds.canon(x), fdds.canon(y)))

    solver = Trees()
    decoder_stats = Counter()
    cached_splits = {}
    calls = 0

    @cache
    def convert(tree):
        # The two implementations assign independent IDs to rooted trees.
        return solver.node(convert(child) for child in fdds.KIDS[tree])

    def row(forest, translate):
        return tuple(sorted((translate(tree), count) for tree, count in forest.items()))

    def check_split(observed, left_anchor, right_anchor, brute_splits):
        nonlocal calls
        key = row(observed, convert), convert(left_anchor), convert(right_anchor)
        if key not in cached_splits:
            forest, left, right = key
            splits = local_pairs_twogrades(solver, Counter(dict(forest)), left,
                                          right, decoder_stats)
            cached_splits[key] = {(tuple(sorted(x.items())), tuple(sorted(y.items())))
                                  for x, y in splits}
        expected = {(row(x, convert), row(y, convert)) for x, y in brute_splits}
        assert cached_splits[key] == expected
        calls += 1

    for key, target in sorted(representatives.items(), key=lambda item: len(item[1])):
        result, _ = search(target, observe_split=check_split)
        assert result == oracle[key], (target, result, oracle[key])
    pairs = sum(len(splits) for splits in oracle.values())
    assert (len(oracle), pairs, calls, len(cached_splits)) == (203, 428, 17343, 206)
    assert all(cached_splits.values())
    print(f"Catalogue through 12 states: {len(oracle)} decomposable targets, "
          f"{pairs} ordered pairs, zero differences.")
    print(f"Anchored splits on the search paths: {calls} non-unit calls, "
          f"{len(cached_splits)} distinct inputs, zero differences.")


def main():
    tables, pairs, levels, unit_levels = arbitrary_endpoints()
    assert (tables, pairs, levels, unit_levels) == (3413, 45233, 139777, 81111)
    print(f"Width budgets: {tables} labelled tables, {pairs} endpoint pairs, "
          f"{levels} levels, {unit_levels} unit budgets.")
    factors, products, centres, steps = true_guides()
    assert (factors, products, centres, steps) == (20, 400, 961, 2194)
    print(f"Minimum-width guides: {products} ordered products, "
          f"{centres} cyclic starting cells, {steps} guide steps.")
    compare_catalogue()
    check_420()
    print("r N old_X (guide_Y,guide_X) local_mass strict_budget rejected_axes signed_grades")
    for r in range(1, 5):
        check_family((3, 5, 7, 11)[:r])
    print("All global checks passed.")


if __name__ == "__main__":
    main()
