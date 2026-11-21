from collections import Counter, defaultdict, deque
from functools import cache


class Trees:
    """Intern finite rooted trees, keeping children in Naquin–Gadouleau order."""

    def __init__(self):
        self.children = [()]
        self.height = [0]
        self.ids = {(): 0}
        self.codes = {0: (0,)}
        self.products = {}
        self.quotients = {}
        self.phase_cache = {}

    def code(self, tree):
        if tree not in self.codes:
            queue = deque([tree])
            degrees = []
            while queue:
                children = self.children[queue.popleft()]
                degrees.append(len(children))
                queue.extend(children)
            self.codes[tree] = tuple(degrees)
        return self.codes[tree]

    def node(self, branches):
        branches = tuple(branches)
        if len(branches) > 1:
            branches = tuple(sorted(branches, key=self.code))
        if branches not in self.ids:
            self.ids[branches] = len(self.children)
            self.children.append(branches)
            self.height.append(1 + max(self.height[c] for c in branches))
        return self.ids[branches]

    def product(self, left, right):
        """Levelwise product. A leaf is absorbing."""
        key = (left, right)
        stack = [(key, False)]
        while stack:
            pair, finish = stack.pop()
            if pair in self.products:
                continue
            a, b = pair
            if not a or not b:
                self.products[pair] = 0
                continue
            pairs = [(x, y) for x in self.children[a] for y in self.children[b]]
            if finish:
                self.products[pair] = self.node(self.products[p] for p in pairs)
            else:
                stack.append((pair, True))
                stack.extend((p, False) for p in set(pairs) if p not in self.products)
        return self.products[key]

    def _divide_steps(self, target, divisor):
        """Yield the strictly smaller subtree divisions needed by Figure 6."""
        if self.height[target] > self.height[divisor]:
            return None
        residual = Counter(self.children[target])
        branches = []
        while residual:
            h = max(self.height[c] for c in residual)
            target_child = min((c for c in residual if self.height[c] == h), key=self.code)
            left_child = min((c for c in self.children[divisor] if self.height[c] >= h),
                             key=self.code)
            child = yield (target_child, left_child)
            if child is None:
                return None
            image = Counter(self.product(c, child) for c in self.children[divisor])
            if any(residual[c] < count for c, count in image.items()):
                return None
            residual -= image
            branches.append(child)
        return self.node(branches)

    def divide(self, target, divisor):
        """Iterative driver for finite tree division, sharing repeated subproblems."""
        key = (target, divisor)
        if key in self.quotients:
            return self.quotients[key]
        stack = [(key, self._divide_steps(*key))]
        answer = None
        while stack:
            pair, task = stack[-1]
            try:
                child = task.send(answer)
            except StopIteration as result:
                answer = result.value
                self.quotients[pair] = answer
                stack.pop()
                continue
            if child in self.quotients:
                answer = self.quotients[child]
            else:
                stack.append((child, self._divide_steps(*child)))
                answer = None
        return answer

    def phase(self, tree, n):
        """Read the first 2n decorations backwards along the unique long branch."""
        key = tree, n
        if key in self.phase_cache:
            return self.phase_cache[key]
        backwards = []
        for _ in range(2 * n):
            children = list(self.children[tree])
            long = [c for c in children if self.height[c] >= n]
            if len(long) != 1:
                return None
            tree = long[0]
            children.remove(tree)
            backwards.append(self.node(children))
        for p in range(1, n + 1):
            if backwards[p:] == backwards[:-p]:
                phase = tuple(reversed(backwards[:p]))
                self.phase_cache[key] = phase
                return phase
        return None


def separating_weights(signatures):
    sigs = [dict(s) for s in signatures]
    keys = sorted(set().union(*sigs))
    constraints = defaultdict(list)
    for i, a in enumerate(sigs):
        for j, b in enumerate(sigs[:i]):
            last = max(k for k in a.keys() | b.keys() if a.get(k, 0) != b.get(k, 0))
            constraints[last].append((i, j))
    scores = [0] * len(sigs)
    weights = {}
    for key in keys:
        bad = set()
        for i, j in constraints[key]:
            slope = sigs[i].get(key, 0) - sigs[j].get(key, 0)
            difference = scores[j] - scores[i]
            if difference % slope == 0:
                bad.add(difference // slope)
        w = next(x for x in range(len(bad) + 1) if x not in bad)
        weights[key] = w
        scores = [s + w*a.get(key, 0) for s, a in zip(scores, sigs)]
    assert len(set(scores)) == len(scores)
    return weights


@cache
def prime_exponents(n):
    out = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            n //= p
            out[p] = out.get(p,0)+1
        p += 1
    if n > 1:
        out[n] = 1
    return out
