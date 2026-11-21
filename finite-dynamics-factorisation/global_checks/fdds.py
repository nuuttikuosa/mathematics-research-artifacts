# Test harness for the research note: brute-force oracle + the note's global search
# (phase grid + two guides + level-by-level extension), with the algebraic local
# splitter replaced by exhaustive enumeration of anchored forest splits.
import itertools, sys
from collections import Counter
from fractions import Fraction
from math import gcd
sys.setrecursionlimit(100000)

# ---------- hash-consed rooted trees ----------
ID = {}; KIDS = []; HT = []; SZ = []; MAXDEG = []
def node(kids):
    key = tuple(sorted(kids))
    i = ID.get(key)
    if i is None:
        i = len(KIDS); ID[key] = i; KIDS.append(key)
        HT.append(0 if not key else 1 + max(HT[c] for c in key))
        SZ.append(1 + sum(SZ[c] for c in key))
        MAXDEG.append(max([len(key)] + [MAXDEG[c] for c in key]))
    return i
LEAF = node(())

_prod = {}
def prod(a, b):
    if a > b: a, b = b, a
    k = (a, b); r = _prod.get(k)
    if r is None:
        r = node([prod(x, y) for x in KIDS[a] for y in KIDS[b]])
        _prod[k] = r
    return r

_cut = {}
def cut(a, k):
    if k <= 0: return LEAF
    key = (a, k); r = _cut.get(key)
    if r is None:
        r = a if HT[a] <= k else node([cut(c, k - 1) for c in KIDS[a]])
        _cut[key] = r
    return r

def fprod(F, G):
    out = Counter()
    for a, ma in F.items():
        for b, mb in G.items():
            out[prod(a, b)] += ma * mb
    return out

def solve_nonneg_int(cols, rhs):
    """cols: list of Counters (vectors), rhs Counter. Unique solution expected."""
    keys = sorted(set(rhs) | set(k for c in cols for k in c))
    n = len(cols)
    M = [[Fraction(c.get(k, 0)) for c in cols] + [Fraction(rhs.get(k, 0))] for k in keys]
    piv_cols = []; r = 0
    for col in range(n):
        pr = next((i for i in range(r, len(M)) if M[i][col] != 0), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]
        pv = M[r][col]; M[r] = [x / pv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][col] != 0:
                fct = M[i][col]; M[i] = [x - fct * y for x, y in zip(M[i], M[r])]
        piv_cols.append(col); r += 1
    for i in range(r, len(M)):
        if M[i][n] != 0: return None
    if len(piv_cols) != n: raise RuntimeError("dependent columns: cancellation violated?")
    x = [None] * n
    for i, col in enumerate(piv_cols): x[col] = M[i][n]
    if any(v.denominator != 1 or v < 0 for v in x): return None
    return [int(v) for v in x]

_div = {}
def div(t, a):
    """Equal-height tree division t / a, or None."""
    if HT[t] != HT[a]: return None
    key = (t, a)
    if key in _div: return _div[key]
    if HT[t] == 0:
        _div[key] = LEAF; return LEAF
    DT = Counter(KIDS[t]); DA = KIDS[a]
    astar = max(DA, key=lambda c: HT[c])
    cands = set()
    for c in DT:
        q = div(c, cut(astar, HT[c]))
        if q is not None: cands.add(q)
    cands = sorted(cands)
    res = None
    if cands:
        cols = [Counter(prod(ai, q) for ai in DA) for q in cands]
        x = solve_nonneg_int(cols, DT)
        if x is not None:
            kids = [q for q, mult in zip(cands, x) for _ in range(mult)]
            cand = node(kids)
            if prod(cand, a) == t: res = cand
    _div[key] = res
    return res

def fdiv(P, A):
    """Division of forests whose trees all have the same height."""
    if not A or not P: return None
    hs = {HT[x] for x in P} | {HT[x] for x in A}
    if len(hs) != 1: return None
    q = div(node(list(P.elements())), node(list(A.elements())))
    if q is None: return None
    Q = Counter(KIDS[q])
    return Q if fprod(A, Q) == P else None

# ---------- functional graphs ----------
def preds_of(f):
    P = [[] for _ in f]
    for u, v in enumerate(f): P[v].append(u)
    return P

def cyclic_states(f):
    n = len(f); onc = [False] * n
    for s in range(n):
        x = s
        for _ in range(n): x = f[x]
        onc[x] = True
    # every state reached after n steps lies on a cycle; mark whole cycles
    for s in range(n):
        if onc[s]:
            y = f[s]
            while y != s: onc[y] = True; y = f[y]
    return onc

def cycles(f):
    onc = cyclic_states(f); seen = set(); out = []
    for s in range(len(f)):
        if onc[s] and s not in seen:
            cyc = [s]; seen.add(s); y = f[s]
            while y != s: cyc.append(y); seen.add(y); y = f[y]
            out.append(cyc)
    return out

def canon(f):
    onc = cyclic_states(f); P = preds_of(f); memo = {}
    def tid(v):
        if v in memo: return memo[v]
        r = node([tid(u) for u in P[v] if not onc[u]]); memo[v] = r; return r
    words = []
    for cyc in cycles(f):
        w = [tid(c) for c in cyc]
        words.append(min(tuple(w[r:] + w[:r]) for r in range(len(w))))
    return tuple(sorted(words))

def fproduct(fx, fy):
    ny = len(fy)
    return [fx[x] * ny + fy[y] for x in range(len(fx)) for y in range(ny)]

def connected(f): return len(cycles(f)) == 1

def from_phases(phases):
    """phases[i] = list of transient tree ids attached to cyclic state i; cycle i->i+1."""
    p = len(phases); f = list((i + 1) % p for i in range(p))
    def attach(t, target):
        idx = len(f); f.append(target)
        for c in KIDS[t]: attach(c, idx)
    for i, ts in enumerate(phases):
        for t in ts: attach(t, i)
    return f

def connected_types(n):
    seen = {}
    for f in itertools.product(range(n), repeat=n):
        f = list(f)
        if connected(f):
            c = canon(f)
            if c not in seen: seen[c] = f
    return seen
