"""Independent finite audit of the paper's geometric SAT construction.

No project implementation is imported. This checks graph/line geometry and
exhaustive component-bit bijections; it is not a substitute for the proof.
"""
from collections import defaultdict, deque, Counter
from itertools import product
from pathlib import Path
import json, random, time

def reflected_path(alphabet, dimension):
    path=[()]
    for _ in range(dimension):
        path=[(x,)+tail for i,x in enumerate(alphabet)
              for tail in (path if i%2==0 else path[::-1])]
    return path

def permutation_network(u,g,perm):
    words=list(product(range(u),repeat=g))
    assert set(perm)==set(words)==set(perm.values())
    if g==1: return [perm]
    tails=list(product(range(u),repeat=g-1))
    available=set(words); colors={}
    for color in range(u):
        adjacent={a:sorted(w for w in available if w[1:]==a) for a in tails}
        matched={}
        def visit(a,seen):
            for edge in adjacent[a]:
                b=perm[edge][1:]
                if b in seen: continue
                seen.add(b)
                if b not in matched or visit(matched[b][1:],seen):
                    matched[b]=edge
                    return True
            return False
        assert all(visit(a,set()) for a in tails)
        assert len(matched)==len(tails)
        for edge in matched.values(): colors[edge]=color
        available.difference_update(matched.values())
    first={w:(colors[w],)+w[1:] for w in words}
    sub=[]
    for c in range(u):
        p={w[1:]:perm[w][1:] for w in words if colors[w]==c}
        sub.append(permutation_network(u,g-1,p))
    middle=[{(c,)+w:(c,)+sub[c][i][w] for c in range(u) for w in tails}
            for i in range(2*g-3)]
    last={(colors[w],)+perm[w][1:]:perm[w] for w in words}
    stages=[first]+middle+[last]
    expected=list(range(g))+list(range(g-2,-1,-1))
    assert len(stages)==2*g-1
    for p,j in zip(stages,expected):
        assert set(p)==set(words)==set(p.values())
        assert all(all(a[k]==b[k] for k in range(g) if k!=j) for a,b in p.items())
    for w in words:
        x=w
        for p in stages:x=p[x]
        assert x==perm[w]
    return stages

def check_formula(n,clauses,u=3,exhaustive=True):
    assert n>=1 and all(1<=len(c)<=3 for c in clauses)
    clauses=[tuple(c)+(c[-1],)*(3-len(c)) for c in clauses]
    tests=[('clause',dict(enumerate(c))) for c in clauses]
    for x in range(1,n+1):
        for j in range(3):
            for k in range(3):
                if j!=k:tests.append(('sign',{j:x,k:-x}))
    g=1
    while u**g<len(tests):g+=1
    h=3*(2*g-1); k=1
    while (2*u-1)**k<h:k+=1
    pairpath=[(a,a) for a in range(u)]
    pairpath=[z for a in range(u) for z in ([(a,a),(a,a+1)] if a+1<u else [(a,a)])]
    times=[sum(point,()) for point in reflected_path(pairpath,k)[:h]]
    words=list(product(range(u),repeat=g)); ham=reflected_path(list(range(u)),g)
    graph=defaultdict(set); owners=defaultdict(set); marked={}; testpoints=[]
    def point(word,kind,t,lane,plane):
        return tuple(word)+((sum(word)+kind)%u,)+times[t]+(lane,plane)
    def addnode(v,owner):graph[v];owners[v].add(owner)
    def copy(a,b,owner):
        different=[j for j,(x,y) in enumerate(zip(a,b)) if x!=y]
        assert len(different)==1,(a,b)
        j=different[0];key=(j,a[:j]+a[j+1:])
        endpoints=frozenset((a,b))
        assert key not in marked or marked[key]==endpoints
        marked[key]=endpoints
        addnode(a,owner);addnode(b,owner)
        graph[a].add(b);graph[b].add(a)
    def copy_path(points,owner):
        for v in points:addnode(v,owner)
        for a,b in zip(points,points[1:]):copy(a,b,owner)
    for plane in range(3):
        inputs={words[i]:request[plane] for i,(_,request) in enumerate(tests) if plane in request}
        ordered=sorted(inputs,key=lambda w:(inputs[w],w))
        perm={w:ham[i] for i,w in enumerate(ordered)}
        sparein=[w for w in words if w not in inputs]
        spareout=ham[len(ordered):]
        perm.update(zip(sparein,spareout))
        stages=permutation_network(u,g,perm)
        for initial,literal in inputs.items():
            owner=(literal,plane); a=initial
            for s,stage in enumerate(stages):
                b=stage[a]; coordinate=list(range(g))+list(range(g-2,-1,-1))
                j=coordinate[s];lane=a[j]+1; entrance=3*s
                path=[point(a,0,entrance,0,plane),point(a,0,entrance,lane,plane),point(a,0,entrance+1,lane,plane)]
                x=a
                for _ in range((b[j]-a[j])%u):
                    path.append(point(x,1,entrance+1,lane,plane))
                    xx=list(x);xx[j]=(xx[j]+1)%u;x=tuple(xx)
                    path.append(point(x,0,entrance+1,lane,plane))
                assert x==b
                path.extend([point(b,0,entrance+2,lane,plane),point(b,0,entrance+2,0,plane)])
                copy_path(path,owner)
                if s+1<len(stages):copy(path[-1],point(b,0,entrance+3,0,plane),owner)
                a=b
            assert a==perm[initial]
        for first,second in zip(ordered,ordered[1:]):
            if inputs[first]!=inputs[second]:continue
            a,b=perm[first],perm[second]
            differences=[i for i in range(g) if a[i]!=b[i]]
            assert len(differences)==1
            j=differences[0];assert abs(a[j]-b[j])==1
            lower=a if a[j]<b[j] else b
            copy_path([point(a,0,h-1,0,plane),point(lower,1,h-1,0,plane),point(b,0,h-1,0,plane)],(inputs[first],plane))
    for i,(kind,request) in enumerate(tests):
        points=[point(words[i],0,0,0,p) for p in request]
        key=(len(points[0])-1,points[0][:-1])
        assert key not in marked
        testpoints.append((kind,points,key))
    occupancy=defaultdict(set)
    for v in graph:
        assert len(owners[v])==1,('unexpected shared node',v,owners[v])
        for j in range(len(v)):occupancy[j,v[:j]+v[j+1:]].add(v)
    for line,endpoints in marked.items():assert occupancy[line]==set(endpoints),('extra copying-line endpoint',line)
    for _,points,line in testpoints:assert occupancy[line]==set(points),('extra test-line endpoint',line)
    def chi(v):
        word=v[:g];eps=(v[g]-sum(word))%u
        assert eps in (0,1)
        ts=v[g+1:g+1+2*k]
        return (eps+sum(ts[j+1]-ts[j] for j in range(0,2*k,2))+(v[-2]!=0))%2
    for a,neighbors in graph.items():
        for b in neighbors:assert chi(a)!=chi(b)
    comp={}; parity={}; component_owners=[]
    anchors={v for _,points,_ in testpoints for v in points}
    for seed in sorted(anchors):
        assert chi(seed)==0
        if seed in comp:continue
        ix=len(component_owners);comp[seed]=ix;parity[seed]=0;queue=deque([seed]);own=set()
        while queue:
            a=queue.popleft();own.update(owners[a])
            for b in graph[a]:
                if b in comp:assert comp[b]==ix and parity[b]==1-parity[a]
                else:comp[b]=ix;parity[b]=1-parity[a];queue.append(b)
        assert len(own)==1
        component_owners.append(next(iter(own)))
    assert len(comp)==len(graph), 'unanchored extra component'
    assert len(component_owners)==6*n and len(set(component_owners))==6*n
    for v in anchors:assert parity[v]==0
    constraints=[(kind,[comp[v] for v in points]) for kind,points,_ in testpoints]
    index={owner:i for i,owner in enumerate(component_owners)}
    expected={bits for bits in product((0,1),repeat=n) if all(any(bits[abs(x)-1]==(x>0) for x in c) for c in clauses)}
    count=Counter()
    assignments=product((0,1),repeat=6*n) if exhaustive else ([bits[abs(lit)-1] if lit>0 else 1-bits[abs(lit)-1] for lit,p in component_owners] for bits in product((0,1),repeat=n))
    for bits in assignments:
        if any((sum(bits[i] for i in ids)!=1 if kind=='sign' else not any(bits[i] for i in ids)) for kind,ids in constraints):continue
        decoded=tuple(bits[index[x,0]] for x in range(1,n+1))
        assert all(bits[i]==(decoded[abs(lit)-1] if lit>0 else 1-decoded[abs(lit)-1]) for i,(lit,p) in enumerate(component_owners))
        count[decoded]+=1
    assert set(count)==expected and all(c==1 for c in count.values())
    return dict(n=n,m=len(clauses),u=u,g=g,vertices=len(graph),copying_lines=len(marked),components=len(component_owners),solutions=len(expected),exhaustive_component_bits=exhaustive)

def main():
    started=time.perf_counter();rng=random.Random(261126)
    cases=[]
    unary_clauses=[c for size in (1,2,3) for c in product((-1,1),repeat=size)]
    for m in range(3):
        for clauses in product(unary_clauses,repeat=m):cases.append(check_formula(1,clauses))
    for _ in range(60):
        clauses=[tuple(rng.choice((-2,-1,1,2)) for _ in range(rng.randint(1,3))) for _ in range(rng.randint(0,8))]
        cases.append(check_formula(2,clauses,u=rng.choice((3,4,5))))
    for _ in range(12):
        clauses=[tuple(rng.choice((-3,-2,-1,1,2,3)) for _ in range(rng.randint(1,3))) for _ in range(rng.randint(0,12))]
        cases.append(check_formula(3,clauses,u=rng.choice((3,4,5)),exhaustive=False))
    out=dict(cases=len(cases),exhaustive_component_bit_cases=sum(c['exhaustive_component_bits'] for c in cases),additional_assignment_extension_cases=sum(not c['exhaustive_component_bits'] for c in cases),max_vertices=max(c['vertices'] for c in cases),max_g=max(c['g'] for c in cases),max_u=max(c['u'] for c in cases),all_checks_passed=True,seconds=time.perf_counter()-started,case_results=cases)
    Path(__file__).with_name('results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in out.items() if k!='case_results'},indent=2))

if __name__=='__main__':main()
