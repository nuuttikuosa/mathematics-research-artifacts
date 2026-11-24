"""Count the exact weighted tables of Sections 4 and 5 on K4 and two disjoint K4s.

Every row and column is generated from the binary weight formulas. The count
uses only their integer weights and targets, so the lcm of the column weights
and the permutations' successor tables are not needed.
"""

import json
from collections import Counter, defaultdict, deque
from itertools import combinations
from math import lcm


def numerical_margins(vertices, edges):
    """Build the rows, columns and targets using the symbols of the paper."""
    layers = range(3)
    vertex_copies = [("v", v, layer) for v in vertices for layer in layers]
    edge_elements = [("e", edge) for edge in edges]
    elements = vertex_copies + edge_elements
    triples = {
        (edge, layer): (("v", edge[0], layer), ("v", edge[1], layer), ("e", edge))
        for edge in edges for layer in layers
    }

    base = 1024
    names = [(letter, h) for letter in ("p", "g", "r") for h in triples]
    names.extend(("T", t) for t in elements)
    row_scale = base ** (len(names) + 3)
    row_symbol = {
        name: (3 if name[0] == "T" else 1) * row_scale + base**slot
        for slot, name in enumerate(names)
    }
    row_weight = {}
    for h, members in triples.items():
        p, g, r = (row_symbol[letter, h] for letter in ("p", "g", "r"))
        row_weight["P", h] = p
        row_weight["G", h] = g
        row_weight["R", h] = r
        row_weight["K", h] = p + g - r
        for t in members:
            row_weight["F", h, t] = row_symbol["T", t] - p - g

    b = lcm(*row_weight.values())
    names = [(letter, t) for letter in ("B", "H") for t in elements]
    names.extend(("D", v) for v in vertices)
    column_scale = base ** (len(names) + 3)
    residue = {"B": 1, "H": 2, "D": 0}
    column_symbol = {}
    for slot, name in enumerate(names):
        c = residue[name[0]]
        column_symbol[name] = b * (c * column_scale + base**slot) + c

    column_weight = {}
    column_target = {}
    for t in elements:
        if t[0] == "v":
            delta = column_symbol["D", t[1]]
        else:
            v, w = t[1]
            delta = -column_symbol["D", v] - column_symbol["D", w]
        for colour, shift in enumerate((0, 1, 3)):
            x = column_symbol["B", t] + shift * delta
            column_weight["X", t, colour] = x
            column_weight["Y", t, colour] = column_symbol["H", t] - x
            column_target["X", t, colour] = column_target["Y", t, colour] = row_symbol["T", t]

    row_target = {}
    for h, members in triples.items():
        row_target["P", h] = row_target["G", h] = sum(column_symbol["B", t] for t in members)
        row_target["R", h] = row_target["K", h] = sum(
            column_symbol["H", t] - column_symbol["B", t] for t in members
        )
        for t in members:
            row_target["F", h, t] = column_symbol["H", t]

    assert all(0 < w < b and b % w == 0 for w in row_weight.values())
    assert all(w % b == 1 for w in column_weight.values())
    for weights in (row_weight, column_weight):
        assert len(set(weights.values())) == len(weights)
        assert max(weights.values()) < 2 * min(weights.values())
    return row_weight, column_weight, row_target, column_target


def margin_options(weights, targets):
    """List every multiset that sums to each target. Repetitions are allowed."""
    names = list(weights)
    values = [weights[name] for name in names]
    modulus = (1 << 61) - 1
    residues = [value % modulus for value in values]
    pair_index = defaultdict(list)
    for j, right in enumerate(residues):
        for i in range(j + 1):
            pair_index[(residues[i] + right) % modulus].append((i, j))

    # A modular match selects candidates; full integer equality admits a choice.
    solutions = {}
    lightest, heaviest = min(values), max(values)
    for target in set(targets.values()):
        fewest = (target + heaviest - 1) // heaviest
        most = target // lightest
        assert fewest == most and most in (2, 3)
        choices = []
        if most == 2:
            for i, j in pair_index[target % modulus]:
                if values[i] + values[j] == target:
                    choices.append(Counter([names[i], names[j]]))
        else:
            target_residue = target % modulus
            for k, last in enumerate(residues):
                for i, j in pair_index[(target_residue - last) % modulus]:
                    if j <= k and values[i] + values[j] + values[k] == target:
                        choices.append(Counter([names[i], names[j], names[k]]))
        solutions[target] = choices
    return {name: solutions[target] for name, target in targets.items()}


def count_compatible_tables(options):
    """Match the two descriptions of each cell, then branch over surviving choices."""
    neighbours = {name: set() for name in options}
    for name, choices in options.items():
        for choice in choices:
            for other in choice:
                neighbours[name].add(other)
                neighbours[other].add(name)

    def search(available, pending):
        queue = deque(pending)
        while queue:
            first, second = queue.popleft()
            counts = {options[first][i][second] for i in available[first]}
            keep = tuple(i for i in available[second] if options[second][i][first] in counts)
            if not keep:
                return 0
            if len(keep) < len(available[second]):
                available[second] = keep
                queue.extend((second, other) for other in neighbours[second] if other != first)

        undecided = [name for name, choices in available.items() if len(choices) > 1]
        if not undecided:
            selected = {name: options[name][choices[0]] for name, choices in available.items()}
            assert all(
                selected[name][other] == selected[other][name]
                for name, others in neighbours.items() for other in others
            )
            return 1
        name = min(undecided, key=lambda n: len(available[n]))
        total = 0
        for index in available[name]:
            branch = available.copy()
            branch[name] = (index,)
            total += search(branch, ((name, other) for other in neighbours[name]))
        return total

    if any(not choices for choices in options.values()):
        return 0
    available = {name: tuple(range(len(choices))) for name, choices in options.items()}
    return search(available, ((name, other) for name in options for other in neighbours[name]))


def check_graph(copies):
    vertices = range(4 * copies)
    edges = [
        (4 * component + v, 4 * component + w)
        for component in range(copies) for v, w in combinations(range(4), 2)
    ]
    rows, columns, row_targets, column_targets = numerical_margins(vertices, edges)
    row_options = margin_options(columns, row_targets)
    column_options = margin_options(rows, column_targets)
    count = count_compatible_tables(row_options | column_options)
    # K4 has six labelled edge-colourings and two layer choices per component.
    assert count == 12**copies
    return {
        "rows": len(rows),
        "columns": len(columns),
        "max_row_options": max(map(len, row_options.values())),
        "max_column_options": max(map(len, column_options.values())),
        "tables": count,
    }


if __name__ == "__main__":
    print(json.dumps({"K4": check_graph(1), "2K4": check_graph(2)}, indent=2))
