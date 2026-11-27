"""Independent exact cycle-arithmetic audit; no manuscript code imported."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import gcd,isqrt,prod
from pathlib import Path
import json,random,time

HERE=Path(__file__).resolve().parent
def clean(a): return {d:c for d,c in a.items() if c}
def key(a): return tuple(sorted(clean(a).items()))
def mass(a): return sum(d*c for d,c in a.items())
def add(a,b):
    c=Counter(a);c.update(b);return clean(c)
def multiply(a,b):
    c=Counter()
    for x,cx in a.items():
        for y,cy in b.items():
            g=gcd(x,y);c[x*y//g]+=g*cx*cy
    return clean(c)
def divisors(n):
    s=set()
    for d in range(1,isqrt(n)+1):
        if n%d==0:s.update((d,n//d))
    return sorted(s)
@lru_cache(None)
def systems(n,lower=1):
    if n==0:return ((),)
    return tuple(((d,1),)+r for d in range(lower,n+1) for r in systems(n-d,d))
def all_systems(n):
    for entries in systems(n):
        c=Counter()
        for d,m in entries:c[d]+=m
        yield dict(c)

def divide(a,t):
    """Generic exhaustive nonnegative integer column solver, no support lemma."""
    if mass(t)%mass(a):return []
    lengths=sorted(set(d for x in t for d in divisors(x)))
    rows=sorted(t); target=tuple(t[r] for r in rows)
    ds=[];cols=[]
    for d in lengths:
        col=multiply(a,{d:1})
        if any(x not in t or c>t[x] for x,c in col.items()):continue
        ds.append(d);cols.append(tuple(col.get(r,0) for r in rows))
    results=[]
    def solve(res,remaining,chosen):
        if any(x<0 for x in res):return
        # A zero equation forces all incident nonnegative variables to zero.
        remaining=tuple(j for j in remaining if not any(v and res[i]==0 for i,v in enumerate(cols[j])))
        if not any(res):
            results.append({ds[j]:v for j,v in chosen.items() if v});return
        if not remaining:return
        bounds={j:min(res[i]//v for i,v in enumerate(cols[j]) if v) for j in remaining}
        candidates=[]
        for i,b in enumerate(res):
            if not b:continue
            inc=[j for j in remaining if cols[j][i]]
            if not inc:return
            estimate=prod(bounds[j]+1 for j in inc)
            candidates.append((estimate,len(inc),i,inc))
        _,_,i,inc=min(candidates)
        rest=tuple(j for j in remaining if j not in inc)
        def assign(k,left,delta,values):
            if k==len(inc):
                if left==0:
                    nxt=tuple(x-y for x,y in zip(res,delta))
                    solve(nxt,rest,{**chosen,**values})
                return
            j=inc[k];w=cols[j][i]
            last=k==len(inc)-1
            choices=[left//w] if last and left%w==0 else ([] if last else range(min(bounds[j],left//w)+1))
            for v in choices:
                if v>bounds[j]:continue
                nd=tuple(x+v*y for x,y in zip(delta,cols[j]))
                if any(x>y for x,y in zip(nd,res)):continue
                assign(k+1,left-v*w,nd,{**values,j:v})
        assign(0,res[i],(0,)*len(res),{})
    solve(target,tuple(range(len(ds))),{})
    assert len(results)==len({key(x) for x in results})
    assert all(multiply(a,x)==t for x in results)
    return results

def pairs(t):
    found=set();n=mass(t)
    for d in divisors(n):
        if d<2 or d*d>n:continue
        for a in all_systems(d):
            for b in divide(a,t):
                if mass(b)>=2:found.add(tuple(sorted((key(a),key(b)))))
    return found

def interval_case(groups,omega,caps,intervals):
    full=list(product(*groups));partial=list(product(*[(1,)+tuple(g) for g in groups]))
    def f(e,x):
        return sum(prod(v)//prod(e)*x[v] for v in full if all(a==1 or a==b for a,b in zip(e,v)))
    a={omega:1,**{prod(g):1 for g in groups}}
    def y(x):
        out={}
        for e in partial:
            lo,hi=intervals[e];u=hi-lo;s=e.count(1)
            z=f(e,x)-lo if s%2==0 else hi-f(e,x)
            out[omega*prod(e)]=z;out[prod(e)]=omega*(u-z)
        return clean(out)
    zero={v:0 for v in full};b=multiply(a,y(zero))
    assert b and all(c>=0 and int(c)==c for c in b.values())
    feasible=[]
    for vals in product(*[range(caps[v]+1) for v in full]):
        x=dict(zip(full,vals))
        if all(lo<=f(e,x)<=hi for e,(lo,hi) in intervals.items()):feasible.append(y(x))
        assert multiply(a,y(x))==b
    actual=divide(a,b)
    assert {key(q) for q in actual}=={key(q) for q in feasible}
    assert len(feasible)==len({key(q) for q in feasible})
    expected_mass=mass(a)*omega*sum(prod(e)*(hi-lo) for e,(lo,hi) in intervals.items())
    upper=mass(a)*omega*2**len(groups)*sum(prod(v)*caps[v] for v in full)
    assert mass(b)==expected_mass<=upper
    return len(feasible),len(b),mass(b)

def main():
    started=time.perf_counter();rng=random.Random(261126)
    # First validate the independent divisor oracle against complete small systems.
    oracle_checks=0
    for na in range(1,7):
        for a in all_systems(na):
            for nx in range(1,6):
                known={key(multiply(a,x)) for x in all_systems(nx)}
                for tk in known:
                    t=dict(tk)
                    expect={key(x) for x in all_systems(nx) if multiply(a,x)==t}
                    assert {key(x) for x in divide(a,t)}==expect
                    oracle_checks+=1
    # Random intervals containing an integer witness, with arbitrary zero capacities.
    interval_results=[]
    for groups,omega,count in [([(2,3),(5,7)],11,35), ([(5,7),(11,13)],3,35), ([(2,3),(5,7),(11,13)],17,8)]:
        full=list(product(*groups));partial=list(product(*[(1,)+tuple(g) for g in groups]))
        for case in range(count):
            cap_range=3 if len(groups)==2 and case>=25 else 2
            caps={v:rng.randrange(cap_range) for v in full}
            if not any(caps.values()):caps[full[0]]=1
            witness={v:rng.randrange(caps[v]+1) for v in full}
            ints={}
            for e in partial:
                w=sum(prod(v)//prod(e)*caps[v] for v in full if all(a==1 or a==b for a,b in zip(e,v)))
                v=sum(prod(v)//prod(e)*witness[v] for v in full if all(a==1 or a==b for a,b in zip(e,v)))
                ints[e]=(0,caps[e]) if e in caps else (rng.randrange(v+1),rng.randrange(v,w+1))
            interval_results.append(interval_case(groups,omega,caps,ints))
    # A rationally feasible but integrally infeasible case (parity-weight gap).
    groups=[(5,7),(11,13)];full=list(product(*groups));partial=list(product(*[(1,)+tuple(g) for g in groups]))
    caps={v:1 for v in full};ints={}
    for e in partial:
        w=sum(prod(v)//prod(e) for v in full if all(a==1 or a==b for a,b in zip(e,v)))
        ints[e]=(0,1) if e in caps else (w//2,(w+1)//2)
    bad=interval_case(groups,3,caps,ints);assert bad[0]==0
    interval_results.append(bad)
    # Exhaust all nonnegative factor pairs of a family of small two-layer targets.
    layer_cases=0;layer_pairs=0
    for p in [2,3]:
        sums=[s for n in range(1,8) for s in all_systems(n) if all(d%p for d in s)]
        for u in sums:
            for v in sums:
                if max(v.values())>=p:continue
                t=add(u,multiply({p:1},v));actual=pairs(t)
                expected=set()
                for dmass in divisors(gcd(mass(u),mass(v))):
                    if dmass<2:continue
                    for d in all_systems(dmass):
                        for e in divide(d,u):
                            for f in divide(d,v):
                                expected.add(tuple(sorted((key(d),key(add(e,multiply({p:1},f)))))))
                assert actual==expected
                layer_cases+=1;layer_pairs+=len(actual)
    marker_cases=0;marker_pairs=0
    for a in [{1:1,3:1},{3:1,5:1},{5:1,7:1}]:
        for nb in range(1,12):
            for b in all_systems(nb):
                if any(d%2==0 for d in b):continue
                t=add(b,multiply({2:1},a));actual=pairs(t)
                expect={tuple(sorted((key(a),key(add(x,{2:1}))))) for x in divide(a,b)}
                assert actual==expect
                for ak,bk in actual: assert ak!=bk and not pairs(dict(ak)) and not pairs(dict(bk))
                marker_cases+=1;marker_pairs+=len(actual)
    example_a={3:1,5:1};example_b={3:15,5:15,15:8}
    example_quotients=divide(example_a,example_b)
    assert {key(x) for x in example_quotients}=={key({1:15,15:1}),key({3:5,5:3})}
    example_target=add(example_b,multiply({2:1},example_a))
    assert mass(example_target)==256
    example_pairs=pairs(example_target)
    assert example_pairs=={tuple(sorted((key(example_a),key(add(x,{2:1}))))) for x in example_quotients}
    marker_cases+=1;marker_pairs+=len(example_pairs)
    counter_a={3:1,5:1,15:2};counter_t=add({1:1},multiply({2:1},counter_a))
    assert multiply({1:1,6:1},{1:1,10:1})==counter_t
    assert not divide(counter_a,{1:1}) and not pairs(counter_a)
    result=dict(oracle_validation_equations=oracle_checks,interval_instances=len(interval_results),interval_quotients=sum(x[0] for x in interval_results),max_interval_mass=max(x[2] for x in interval_results),rational_only_instance_has_zero_integer_quotients=True,two_layer_targets=layer_cases,two_layer_pairs=layer_pairs,marker_targets=marker_cases,marker_pairs=marker_pairs,strict_threshold_counterexample_verified=True,seconds=time.perf_counter()-started,all_passed=True)
    (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
