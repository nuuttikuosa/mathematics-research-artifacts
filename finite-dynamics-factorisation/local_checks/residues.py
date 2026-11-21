from collections import Counter
from functools import cache
from itertools import combinations, product
import sympy as sp
from tree_arithmetic import prime_exponents, separating_weights

Z = sp.Symbol('z')


def weights_for(vectors):
    vectors = sorted(set(vectors))
    return separating_weights([tuple((i,x) for i,x in enumerate(v) if x) for v in vectors])


def dot(v,weights):
    return sum(weights.get(i,0)*x for i,x in enumerate(v))


@cache
def polynomial_parts(terms):
    """Positive polynomial modulo monomials; include its true leading scalar."""
    minimum = terms[0][0]
    p = sp.Poly.from_dict({(e-minimum,):n for e,n in terms}, Z, domain=sp.QQ)
    result = {('prime',q):n for q,n in prime_exponents(int(p.LC())).items()}
    for factor,power in p.factor_list()[1]:
        result['poly',tuple(factor.monic().all_coeffs())] = power
    return result


def labels(trees, source, truncate=None):
    """Deterministic labels preserving products and all distinctions in source."""
    source = set(source)
    highest = max(trees.height[t] for t in source)
    @cache
    def default_cut(t,j):
        if trees.height[t]<=j:
            return t
        return trees.node(default_cut(c,j-1) for c in trees.children[t]) if j else 0
    cut=default_cut if truncate is None else truncate
    all_trees = set(source)
    stack = list(source)
    while stack:
        t = stack.pop()
        for c in trees.children[t]:
            if c not in all_trees:
                all_trees.add(c)
                stack.append(c)
    groups = [sorted({cut(t,j) for t in all_trees if trees.height[t]>=j})
              for j in range(highest+1)]
    encoded = {0:()}
    codes = {0:0}
    for h in range(1,highest+1):
        vocab = groups[h]
        rows = {}
        for col,t in enumerate(vocab):
            for j in range(h):
                children = [cut(u,j) for u in trees.children[t] if trees.height[u]>=j]
                minimum = min(children,key=codes.__getitem__)
                for i,x in enumerate(encoded[minimum]):
                    if x:
                        rows.setdefault((j,'min',i),[0]*len(vocab))[col] = x
                hist = Counter(codes[u] for u in children)
                for key,x in polynomial_parts(tuple(sorted(hist.items()))).items():
                    rows.setdefault((j,'factor',key),[0]*len(vocab))[col] = x
        data = list(rows.values())
        keep = sp.Matrix(data).T.rref()[1] if data else []
        for col,t in enumerate(vocab):
            encoded[t] = tuple(data[i][col] for i in keep)
        if len({encoded[t] for t in vocab})!=len(vocab):
            raise AssertionError('finite subtree dictionary was not separated')
        weights = weights_for(encoded[t] for t in vocab)
        codes.update({t:dot(encoded[t],weights) for t in vocab})
    return {t:encoded[t] for t in source}


def projected_exponents(spectrum, weights):
    """Linear projection of closed rational forms to univariate log residues."""
    points = list(spectrum)
    rank = len(points[0])
    s = len(points)
    # The origin shift keeps all exponents nonnegative; Euler closure is
    # invariant under the common shift of F and every numerator H_i.
    raw = {v:dot(v,weights) for v in points}
    offset = min(raw.values())
    dense = Counter()
    for v,n in spectrum.items():
        dense[raw[v]-offset] += n
    single = sp.Poly.from_dict({(e,):n for e,n in dense.items()},Z,domain=sp.QQ)
    factors = [(g.monic(),m) for g,m in single.factor_list()[1]]
    nres = len(factors)
    ncols = rank*s+1+nres
    rows = {}
    def put(key,index,value):
        if value:
            row = rows.setdefault(key,Counter())
            row[index] += value
    # For H_i indexed by u and an F-term v, the coefficient of
    # d(H_i/F) in the j direction is F[v]*(u[j]-v[j]).
    for i,j in combinations(range(rank),2):
        for a,u in enumerate(points):
            for v,n in spectrum.items():
                uv = tuple(x+y for x,y in zip(u,v))
                put(('closed',i,j,uv),j*s+a,n*(u[i]-v[i]))
                put(('closed',i,j,uv),i*s+a,-n*(u[j]-v[j]))
    for a,v in enumerate(points):
        for i,w in weights.items():
            put(('curve',raw[v]-offset),i*s+a,w)
    for e,n in dense.items():
        put(('curve',e),rank*s,-n)
    for k,(g,m) in enumerate(factors):
        ratio = single.exquo(g)*g.diff()*sp.Poly(Z,Z,domain=sp.QQ)
        for (e,),coef in ratio.terms():
            put(('curve',e),rank*s+1+k,-coef)
    matrix = sp.MutableSparseMatrix(len(rows),ncols,{(i,j):v for i,row in enumerate(rows.values())
                                   for j,v in row.items() if v})
    if not nres:
        return factors,sp.zeros(0,0),()
    # Coordinates of projected vectors in a rational row basis.
    res = [list(v[rank*s+1:,:]) for v in matrix.nullspace()]
    if not res:
        return factors,sp.zeros(0,nres),()
    R,pivots = sp.Matrix(res).rref()
    return factors,R[:len(pivots),:],pivots


def normalized_polys(spectrum, allowed, mass, stats):
    """All true projected factor polynomials, up to nonzero scalar and shift."""
    allowed = list(allowed)
    weights = weights_for(allowed)
    factors,basis,pivots = projected_exponents(spectrum,weights)
    stats['max_univariate_factors'] = max(stats.get('max_univariate_factors',0),len(factors))
    stats['max_projected_rank'] = max(stats.get('max_projected_rank',0),len(pivots))
    stats['rank_saved'] = max(stats.get('rank_saved',0),len(factors)-len(pivots))
    limits = [min(m,mass-1) for _,m in factors]
    for choice in product(*(range(limits[i]+1) for i in pivots)):
        values = [sum((a*basis[k,j] for k,a in enumerate(choice)),sp.Integer(0))
                  for j in range(len(factors))]
        if any(not v.is_Integer or not 0<=v<=m for v,m in zip(values,limits)):
            continue
        q = sp.Poly(1,Z,domain=sp.QQ)
        for (g,_),v in zip(factors,values):
            if v:
                q *= g**int(v)
        # Primitive integral normalization, leading coefficient positive.
        q = q.clear_denoms(convert=True)[1].primitive()[1]
        if any(c<=0 for _,c in q.terms()) or int(q.eval(1))>mass:
            continue
        yield weights,[(m[0],int(c)) for m,c in q.terms()]


def local_pairs(trees, E, c0, d0, stats):
    """Return all anchored (C,D) with CD=E, c0 in C, d0 in D."""
    if not E or trees.product(c0,d0) not in E:
        return []
    mass = sum(E.values())
    c_at = {t:trees.divide(t,d0) for t in E}
    c_at = {t:c for t,c in c_at.items() if c is not None}
    out = []
    forest = trees.node(E.elements())
    if len(c_at)==1:
        # Only one possible C type, so its multiplicity divides |E|.
        for n in range(1,mass+1):
            if mass%n:
                continue
            cf = trees.node([c0]*n)
            df = trees.divide(forest,cf)
            if df is not None and d0 in trees.children[df] and trees.product(cf,df)==forest:
                out.append((Counter({c0:n}),Counter(trees.children[df])))
        return out
    data = labels(trees,E)
    origin = data[trees.product(c0,d0)]
    delta = {t:tuple(x-y for x,y in zip(v,origin)) for t,v in data.items()}
    cols = list(E)
    full_matrix = sp.Matrix([delta[t] for t in cols])
    keep = full_matrix.rref()[1]
    points = {t:tuple(delta[t][i] for i in keep) for t in cols}
    hist = {points[t]:E[t] for t in cols}
    allowed = {points[t]:c for t,c in c_at.items()}
    seen = set()
    for weights,terms in normalized_polys(hist,allowed,mass,stats):
        positions = {dot(v,weights):c for v,c in allowed.items()}
        mass0 = sum(n for _,n in terms)
        for start in positions:
            if any(start+j not in positions for j,n in terms):
                continue
            primitive = {positions[start+j]:n for j,n in terms}
            if c0 not in primitive:
                continue
            for scalar in range(1,mass//mass0+1):
                if mass%(mass0*scalar):
                    continue
                counts = Counter({u:scalar*n for u,n in primitive.items()})
                key = tuple(sorted(counts.items()))
                if key in seen:
                    continue
                seen.add(key)
                cf = trees.node(counts.elements())
                df = trees.divide(forest,cf)
                if df is not None and d0 in trees.children[df] and trees.product(cf,df)==forest:
                    out.append((counts,Counter(trees.children[df])))
    return out
