from math import gcd
from functools import lru_cache
def lcm(a,b): return a*b//gcd(a,b)
def partitions(n,m=None):
    if m is None: m=n
    if n==0: yield (); return
    for k in range(min(n,m),0,-1):
        for p in partitions(n-k,k): yield (k,)+p
def orders(n):
    gS=gB=gD=0
    for p in partitions(n):
        # S_n: lcm of parts
        o=1
        for k in p: o=lcm(o,k)
        gS=max(gS,o)
        # B_n: each cycle sign +/-: order k or 2k ; D_n: even number of negative cycles
        r=len(p)
        best={}  # set of reachable (lcm, parity)
        states={(1,0)}
        for k in p:
            ns=set()
            for (l,par) in states:
                ns.add((lcm(l,k),par)); ns.add((lcm(l,2*k),par^1))
            states=ns
        for (l,par) in states:
            gB=max(gB,l)
            if par==0: gD=max(gD,l)
    return gS,gB,gD
bad=[]
for n in range(1,31):
    gS,gB,gD=orders(n)
    ok=(gB==2*gS) and (gD==2*gS or n<5)
    if not ok or n<=6: print(n,"g=%d g_B=%d g_D=%d 2g=%d"%(gS,gB,gD,2*gS), "" if ok else "<-- MISMATCH")
    if not ok: bad.append(n)
print("n=1..30 checked; B=2g all n:",all(orders(n)[1]==2*orders(n)[0] for n in range(1,31)),"; D=2g for n>=5:",all(orders(n)[2]==2*orders(n)[0] for n in range(5,31)))
print("D cases n<5:",[(n,orders(n)) for n in range(1,5)])
