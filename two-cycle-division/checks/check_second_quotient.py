"""Compare second-quotient construction against all partitions of the source size.

Run with Python 3.9+:
    python -B check_second_quotient.py

The at-most-four-dimensional feasibility routine is replaced by exact
backtracking on the small instances here. This does not measure Kannan's
asymptotic running time.
"""

from collections import Counter, defaultdict
from functools import lru_cache
from math import gcd, isqrt, lcm


def key(cycles):
    return tuple(sorted((t, n) for t, n in cycles.items() if n))


def multiply(a, b, cycles):
    result = Counter()
    for t, n in cycles:
        for d in (a, b):
            result[lcm(d, t)] += n*gcd(d, t)
    return key(result)


@lru_cache(None)
def partitions(n, least=1):
    if n == 0:
        return ((),)
    return tuple((t,)+tail for t in range(least, n+1)
                 for tail in partitions(n-t, t))


def factors(n):
    primes = []
    d = 2
    while d*d <= n:
        while n % d == 0:
            primes.append(d)
            n //= d
        d += 1
    if n > 1:
        primes.append(n)
    return primes


def least_period(a, u):
    r = u//a
    return r*gcd(a, pow(r, a.bit_length(), a))


def representatives(a, b, target):
    support = list(dict(target))
    pairs = {}
    for u in support:
        if u % a: continue
        for v in support:
            if v % b: continue
            s = lcm(least_period(a,u), least_period(b,v))
            if lcm(a,s)==u and lcm(b,s)==v:
                pairs[s] = (u,v)
    return pairs


def first_two_canonical(a, b, target, pairs):
    columns = [(t, multiply(a, b, ((t, 1),))) for t in sorted(pairs)]
    answers = []
    def rec(i, residual, chosen):
        if len(answers)==2: return
        if not any(residual.values()):
            answers.append(key(chosen)); return
        if i==len(columns): return
        t,col=columns[i]
        bound=min(residual[d]//k for d,k in col)
        for mult in range(bound+1):
            nxt = residual.copy()
            for d,k in col: nxt[d] -= mult*k
            chosen[t] = mult
            rec(i+1,nxt,chosen)
            if len(answers)==2: break
        del chosen[t]
    rec(0,dict(target),{})
    return answers


def find_two(g, target):
    a,b=2*g,3*g
    pairs=representatives(a,b,target)
    answers=first_two_canonical(a,b,target,pairs)
    if not answers: return "no", answers
    if len(answers)>1: return "two representative vectors", answers
    primes = factors(g)
    for s,z in answers[0]:
        u,v=pairs[s]
        H=gcd(u,v)//s
        assert H == g//gcd(g,s)
        if H==1: continue
        p=min(p for p in primes if H%p==0)
        if z>=p:
            second=Counter(dict(answers[0]))
            second[s]-=p
            second[p*s]+=1
            answers.append(key(second))
            return "two lifts", answers
    return "one", answers


def choice_tests():
    stats=Counter()
    for g in range(1,19):
        for n in range(1,15):
            ground=defaultdict(set)
            for parts in partitions(n):
                x=key(Counter(parts))
                ground[multiply(2*g,3*g,x)].add(x)
            cases=set(ground)
            for target in ground:
                for t,_ in target:
                    for d in {1, 2*g, 3*g, t//2}:
                        if 0<d<t:
                            altered=Counter(dict(target)); altered[t]-=1
                            altered[d]+=1; altered[t-d]+=1
                            cases.add(key(altered))
            for target in cases:
                true=ground.get(target, set())
                status,witnesses=find_two(g,target)
                assert min(2,len(true))==len(witnesses), (g,target,true,witnesses)
                assert len(set(witnesses))==len(witnesses) and all(w in true for w in witnesses)
                stats[status]+=1
    return stats


@lru_cache(None)
def trunc_partition_counts(g):
    # The oracle materializes the restricted partition counts up to g and
    # caps each at 2. No prime-factor information is returned to its caller.
    ways=[1]+[0]*g
    for d in range(1,g+1):
        if g%d: continue
        for k in range(d,g+1):
            ways[k]=min(2,ways[k]+ways[k-d])
    return ways


def factoring_tests():
    calls=0
    for g in range(1,501):
        residual=g; recovered=[]
        while residual>1:
            low,high=2,residual
            while low<high:
                mid=(low+high)//2
                calls+=1
                if trunc_partition_counts(residual)[mid]==2: high=mid
                else: low=mid+1
            recovered.append(low)
            residual//=low
        assert recovered==factors(g), (g,recovered)
    print(f"Recovered exact prime factorisations for 1..500 with {calls} nonuniqueness queries.")


if __name__=="__main__":
    stats=choice_tests()
    print("Direct partitions of |X|=1..14, g=1..18, with same-size target perturbations.")
    print(f"All {sum(stats.values())} target checks passed:", dict(sorted(stats.items())))
    factoring_tests()
