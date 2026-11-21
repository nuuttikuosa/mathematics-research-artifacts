from collections import Counter
from itertools import combinations, product
from math import gcd
from residues import labels, local_pairs, projected_exponents, sp, weights_for

T = sp.Symbol('t')


def dot(v, w):
    return sum(x*y for x,y in zip(v,w))


def second_grade(points, c):
    dimension = len(c)
    constraints = [[] for _ in c]
    for a, b, d in combinations(points, 3):
        u = tuple(x-y for x,y in zip(b,a))
        v = tuple(x-y for x,y in zip(d,a))
        cu, cv = dot(c,u), dot(c,v)
        normal = tuple(cu*y-cv*x for x,y in zip(u,v))
        if any(normal):
            last = max(i for i in range(dimension) if normal[i])
            constraints[last].append(normal)
    answer = []
    for i, rows in enumerate(constraints):
        bad = set()
        for v in rows:
            b, rem = divmod(-dot(answer, v[:i]), v[i])
            if rem == 0:
                bad.add(b)
        answer.append(next(j for j in range(len(bad)+1) if j not in bad))
    return tuple(answer)


def determinant(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def hull(points):
    def side(points):
        out = []
        for p in points:
            while len(out)>=2 and determinant(out[-2],out[-1],p)<=0:
                out.pop()
            out.append(p)
        return out
    pts = sorted(points)
    return side(pts)[:-1]+side(reversed(pts))[:-1] if len(pts)>1 else pts


def multiplicity(p):
    return max((m for _,m in p.sqf_list()[1]), default=0)


def planar_slice(spectrum, c=None):
    points = list(spectrum)
    if len(points) < 2:
        raise ValueError("spectrum needs at least two distinct terms")
    if c is None:
        weights = weights_for(points)
        c = tuple(weights.get(i,0) for i in range(len(points[0])))
    assert len(set(dot(v,c) for v in points)) == len(points)
    w = second_grade(points,c)
    raw = {v:(dot(v,c),dot(v,w)) for v in points}
    low = tuple(map(min,zip(*raw.values())))
    project = {v:(a-low[0],b-low[1]) for v,(a,b) in raw.items()}
    freqs = {project[v]:n for v,n in spectrum.items()}
    points2 = list(freqs)
    polygon = hull(points2)
    edge_lengths, face_multiplicities = [], []
    for a,b in zip(polygon, polygon[1:]+polygon[:1]):
        edge_points = [p for p in points2 if not determinant(a,b,p)]
        step = gcd(abs(b[0]-a[0]),abs(b[1]-a[1]))
        unit = (b[0]-a[0]) // step
        # c separates all support, ruling out a vertical edge.
        edge = sp.Poly.from_dict({((p[0]-a[0])//unit,):freqs[p] for p in edge_points},
                                 T, domain=sp.QQ)
        edge_lengths.append(len(edge_points)-1)
        face_multiplicities.append(multiplicity(edge))
    bound = max(face_multiplicities)
    tdeg, ydeg = map(max, zip(*points2))
    for scalar in range(1, 2*tdeg*ydeg+2):
        specialized = sp.Poly.from_dict({(u,):n*scalar**v for (u,v),n in freqs.items()},
                                        T, domain=sp.QQ)
        peak = multiplicity(specialized)
        if peak <= bound:
            return c, w, scalar, specialized, max(edge_lengths), bound
    raise AssertionError("no permissible scalar found within the resultant bound")


def anchor_recovery(target, c, w, scalar, univar, anchor):
    allowed = {dot(v,c):v for v in target}
    out = set()
    for lower in allowed:
        recovered = {}
        for (e,), coeff in univar.terms():
            v = allowed.get(lower+e)
            if v is None:
                break
            recovered[v] = coeff/(sp.Integer(scalar)**dot(w,v))
        else:
            if anchor not in recovered or any(q <= 0 for q in recovered.values()):
                continue
            denom = sp.ilcm(*[q.q for q in recovered.values()]) if len(recovered)>1 else next(iter(recovered.values())).q
            ints = {p:int(q*denom) for p,q in recovered.items()}
            content = gcd(*ints.values())
            out.add(tuple(sorted((p,n//content) for p,n in ints.items())))
    return out


def local_pairs_twogrades(trees, target, c0, d0, stats):
    """Drop-in deterministic solver using the tilted good univariate slice."""
    anchor_tree = trees.product(c0,d0)
    if not target or anchor_tree not in target:
        return []
    allowed = {tree: trees.divide(tree,d0) for tree in target}
    allowed = {tree: child for tree,child in allowed.items() if child is not None}
    if len(allowed)<=1:
        return local_pairs(trees,target,c0,d0,stats)
    data = labels(trees,target)
    origin = data[anchor_tree]
    diffs = {tree: tuple(x-y for x,y in zip(v,origin)) for tree,v in data.items()}
    keep = sp.Matrix(list(diffs.values())).rref()[1]
    positions = {tree:tuple(v[i] for i in keep) for tree,v in diffs.items()}
    spectrum = {positions[tree]:n for tree,n in target.items()}
    dictionary = {positions[tree]:child for tree,child in allowed.items()}
    c,w,scalar,p,_,bound = planar_slice(spectrum)
    lo = min(dot(v,w) for v in spectrum)
    weighted = {v:n*scalar**(dot(v,w)-lo) for v,n in spectrum.items()}
    factors,basis,pivots = projected_exponents(weighted,dict(enumerate(c)))
    assert all(m<=bound for _,m in factors)
    stats["twisted_calls"] += 1
    stats["nonidentity_twists"] += scalar != 1
    stats["max_edge_multiplicity"] = max(stats.get("max_edge_multiplicity",0),bound)
    stats["max_projected_rank"] = max(stats.get("max_projected_rank",0),len(pivots))
    mass = sum(target.values())
    limits = [min(m,mass-1) for _,m in factors]
    recovered = set()
    for seed in product(*(range(limits[i]+1) for i in pivots)):
        values = [sum((a*basis[k,j] for k,a in enumerate(seed)),sp.Integer(0))
                  for j in range(len(factors))]
        if any(not a.is_Integer or not 0<=a<=limit for a,limit in zip(values,limits)):
            continue
        q = sp.Poly(1,T,domain=sp.QQ)
        for (factor,_), exponent in zip(factors,values):
            if exponent:
                q *= sp.Poly.from_dict({m:c for m,c in factor.terms()},T,domain=sp.QQ)**int(exponent)
        recovered.update(anchor_recovery(dictionary,c,w,scalar,q,positions[anchor_tree]))
    out = []
    forest = trees.node(target.elements())
    seen = set()
    for raw in recovered:
        total = sum(n for _,n in raw)
        primitive = {dictionary[v]:n for v,n in raw}
        for factor in range(1,mass//total+1):
            if mass%(factor*total):
                continue
            a = Counter({v:n*factor for v,n in primitive.items()})
            key = tuple(sorted(a.items()))
            if key in seen:
                continue
            seen.add(key)
            at = trees.node(a.elements())
            bt = trees.divide(forest,at)
            if bt is not None and d0 in trees.children[bt] and trees.product(at,bt)==forest:
                out.append((a,Counter(trees.children[bt])))
    return out
