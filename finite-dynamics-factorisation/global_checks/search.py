# The note's search (Sections 2-3), local splitter replaced by exhaustive anchored splits.
import itertools
from collections import Counter
from math import gcd
from fdds import (HT, KIDS, LEAF, MAXDEG, SZ, canon, cut, cycles, cyclic_states,
                  div, fdiv, fprod, fproduct, from_phases, node, preds_of, prod)

class Stats:
    def __init__(s): s.nodes = 0; s.split_calls = 0; s.splits = 0; s.max_splits = 0; s.guard_rej = 0; s.max_m = 0

def search(f, guard=True, bstar=None, all_pairs=True, stats=None, observe_split=None,
           endpoint_pair=None, selected_split=None):
    N = len(f)
    st = stats or Stats()
    if N < 4: return set(), st
    P = preds_of(f); onc = cyclic_states(f); cyc = cycles(f)
    Ls = {len(c) for c in cyc}
    if len(Ls) != 1: return set(), st
    L = Ls.pop(); g = len(cyc)
    depth = [0] * N
    for v in range(N):
        x, d = v, 0
        while not onc[x]: x = f[x]; d += 1
        depth[v] = d
    h = max(depth)
    if bstar is None: bstar = cyc[0][0]
    memoT = {}
    def T(v, j):
        k = (v, j)
        if k not in memoT:
            memoT[k] = LEAF if j == 0 else node([T(u, j - 1) for u in P[v]])
        return memoT[k]
    def fpow(x, k):
        for _ in range(k): x = f[x]
        return x
    cycpred = {}
    for c in cyc:
        for i, s in enumerate(c): cycpred[s] = c[i - 1]
    def act(s, j):
        return Counter(t for t in (T(u, j) for u in P[s]) if HT[t] == j)
    def width(v, j):
        S = {v}
        for _ in range(j): S = {u for x in S for u in P[x]}
        return len(S)
    target = canon(f); found = set()

    # ---------- phase grids (Section 2) ----------
    J = h + L + 2
    e = {b: T(b, J) for c in cyc for b in c}
    U = set(e.values())
    kap = lambda t: cut(t, J - 1)
    def pi(t):
        ks = [c for c in KIDS[t] if HT[c] == J - 1]
        assert len(ks) == 1
        return ks[0]
    kmap = {}
    for u in U: kmap[kap(u)] = u
    assert len(kmap) == len(U), "kappa not injective on U"
    e0 = e[bstar]
    def axis(um, period):
        seq = [e0]; cur = e0
        for _ in range(period):
            tgt = div(prod(pi(cur), kap(um)), pi(e0))
            if tgt is None or tgt not in kmap: return None
            cur = kmap[tgt]; seq.append(cur)
        if seq[1] != um or seq[period] != e0: return None
        return [seq[(-i) % period] for i in range(period)]
    def grids(p, q):
        for um in U:
            ua = axis(um, p)
            if ua is None: continue
            for vm in U:
                va = axis(vm, q)
                if va is None: continue
                cell = {}
                ok = True
                for i in range(p):
                    for l in range(q):
                        t = div(prod(ua[i], va[l]), e0)
                        if t is None or t not in U: ok = False; break
                        cell[(i, l)] = t
                    if not ok: break
                if not ok: continue
                # match diagonal orbits to physical components
                phys = {}; used = set()
                anc = next(c for c in cyc if bstar in c); used.add(id(anc))
                x = bstar
                for t in range(L):
                    if cell[(t % p, t % q)] != e[x]: ok = False; break
                    phys[(t % p, t % q)] = x; x = f[x]
                if not ok: continue
                for i0 in range(p):
                    for l0 in range(q):
                        if (i0, l0) in phys: continue
                        w = [cell[((i0 + t) % p, (l0 + t) % q)] for t in range(L)]
                        hit = None
                        for c in cyc:
                            if id(c) in used: continue
                            for r in range(L):
                                if all(e[c[(r + t) % L]] == w[t] for t in range(L)):
                                    hit = (c, r); break
                            if hit: break
                        if not hit: ok = False; break
                        c, r = hit; used.add(id(c))
                        for t in range(L): phys[((i0 + t) % p, (l0 + t) % q)] = c[(r + t) % L]
                    if not ok: break
                if ok: yield phys

    # ---------- splits ----------
    def splits(E, c0, d0):
        m = sum(E.values())
        cap = {}
        for t, mult in E.items():
            c = div(t, d0)
            if c is not None: cap[c] = cap.get(c, 0) + mult
        if c0 not in cap: return
        keys = sorted(cap)
        ranges = [range(1 if k == c0 else 0, cap[k] + 1) for k in keys]
        total = 1
        for r in ranges: total *= len(r)
        if total > 2_000_000: raise RuntimeError("split enumeration too large")
        for mults in itertools.product(*ranges):
            mc = sum(mults)
            if mc == 0 or m % mc: continue
            C = Counter({k: v for k, v in zip(keys, mults) if v})
            D = fdiv(E, C)
            if D is None or D[d0] < 1: continue
            if fprod(C, D) != E: continue
            yield C, D

    def extend(old, A, j):
        if j == 0: return node(list(A.elements()))
        rem = Counter(KIDS[old])
        for a, mult in A.items():
            c = cut(a, j - 1)
            if rem[c] < mult: return None
            rem[c] -= mult
        return node(list(A.elements()) + list(rem.elements()))

    H = h + 1
    ends = [z for z in range(N) if fpow(z, H) == bstar]
    if endpoint_pair is not None and any(z not in ends for z in endpoint_pair):
        raise ValueError("guide endpoints must have the selected cycle root as their image")
    first_ends, second_ends = (ends, ends) if endpoint_pair is None else (
        [endpoint_pair[0]], [endpoint_pair[1]])
    for p in range(1, L + 1):
        if L % p: continue
        for q in range(1, L + 1):
            if L % q or p * q != g * L or p * q // gcd(p, q) != L: continue
            for phys in grids(p, q):
                for z1 in first_ends:
                    for z2 in second_ends:
                        W1 = [width(fpow(z1, j), j) for j in range(H + 1)]
                        W2 = [width(fpow(z2, j), j) for j in range(H + 1)]
                        stack = [(0, tuple([LEAF] * p), tuple([LEAF] * q))]
                        while stack:
                            j, MX, MY = stack.pop(); st.nodes += 1
                            if j == H:
                                fx = []; fy = []; bad = False
                                for marks, out in ((MX, fx), (MY, fy)):
                                    for t in marks:
                                        tall = [c for c in KIDS[t] if HT[c] == H - 1]
                                        if len(tall) != 1: bad = True; break
                                        rest = list(KIDS[t]); rest.remove(tall[0]); out.append(rest)
                                    if bad: break
                                if bad: continue
                                X = from_phases(fx); Y = from_phases(fy)
                                if len(X) < 2 or len(Y) < 2 or len(X) * len(Y) != N: continue
                                if canon(fproduct(X, Y)) == target:
                                    found.add((canon(X), canon(Y)))
                                    if not all_pairs: return found, st
                                continue
                            k = H - j - 1; a = (-k) % p; b = (-k) % q
                            cellc = phys[(a, b)]
                            v1n, v1 = fpow(z1, j + 1), fpow(z1, j)
                            v2n, v2 = fpow(z2, j + 1), fpow(z2, j)
                            F, G, K = act(cellc, j), act(v1n, j), act(v2n, j)
                            c0 = div(T(v1, j), MX[(a - 1) % p]); d0 = div(T(v2, j), MY[(b - 1) % q])
                            if c0 is None or d0 is None: continue
                            E = fdiv(fprod(G, K), F)
                            if E is None or E[prod(c0, d0)] == 0: continue
                            m = sum(E.values()); st.max_m = max(st.max_m, m)
                            Dj = (W1[j + 1] // W1[j]) * (W2[j + 1] // W2[j])
                            if guard and m > Dj: st.guard_rej += 1; continue
                            st.split_calls += 1; ns = 0
                            SPL = (list(splits(E, c0, d0)) if selected_split is None else
                                   selected_split(j, E, c0, d0))
                            if observe_split is not None and m >= 2:
                                observe_split(E, c0, d0, SPL)
                            for C, D in SPL:
                                ns += 1
                                A = fdiv(G, C); V = fdiv(K, D)
                                if A is None or V is None or fprod(A, V) != F: continue
                                nXa = extend(MX[a], A, j); nYb = extend(MY[b], V, j)
                                if nXa is None or nYb is None: continue
                                NX = [div(T(phys[(i, b)], j + 1), nYb) for i in range(p)]
                                NY = [div(T(phys[(a, l)], j + 1), nXa) for l in range(q)]
                                if None in NX or None in NY or NX[a] != nXa or NY[b] != nYb: continue
                                ok = all(prod(NX[i], NY[l]) == T(phys[(i, l)], j + 1) for i in range(p) for l in range(q))
                                ok = ok and all(cut(NX[i], j) == MX[i] and MX[(i - 1) % p] in KIDS[NX[i]] for i in range(p))
                                ok = ok and all(cut(NY[l], j) == MY[l] and MY[(l - 1) % q] in KIDS[NY[l]] for l in range(q))
                                ok = ok and all(SZ[t] <= (j + 2) * N and MAXDEG[t] <= N for t in NX + NY)
                                if ok: stack.append((j + 1, tuple(NX), tuple(NY)))
                            st.splits += ns; st.max_splits = max(st.max_splits, ns)
    return found, st
