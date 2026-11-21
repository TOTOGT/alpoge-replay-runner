import numpy as np, mpmath as mp, math, sympy
g=0.5772156649015329
eg=math.exp(g)
# ---- sieve sigma up to N
N=2_000_000
sig=np.zeros(N+1,dtype=np.int64)
for d in range(1,N+1):
    sig[d::d]+=d
n=np.arange(N+1,dtype=np.float64)
# Robin: n>5040
m=np.arange(5041,N+1)
R=sig[m]/(m*np.log(np.log(m)))
i=R.argmax(); print("ROBIN n in (5040,2e6]: max sigma/(n loglog n)=%.6f at n=%d ; bound e^g=%.6f ; violations=%d"%(R[i],m[i],eg,(R>=eg).sum()))
# Robin on extremal (superabundant-type) numbers up to 1e60: exponents nonincreasing over first primes
primes=list(sympy.primerange(2,200))
best=(0,0)
cnt=0
def rec(idx,cur,lastexp,sigma_ratio,limit=10**60):
    global best,cnt
    if cur>5040:
        cnt+=1
        r=sigma_ratio/math.log(mp.log(cur)) if False else sigma_ratio/float(mp.log(mp.log(cur)))
        if r>best[0]: best=(r,cur)
    p=primes[idx]; c=cur; e=0; s=1.0
    while True:
        e+=1; c*=p
        if c>limit or e>lastexp: break
        s_ratio=sigma_ratio*(1-float(p)**-(e+1))/(1-1/p)/( (1-float(p)**-e)/(1-1/p) if e>1 else 1)  # incremental
        rec(idx+1,c,e,s_ratio) if idx+1<len(primes) else None
# simpler explicit generation with exact integers
from fractions import Fraction
def gen(limit):
    res=[]
    def go(idx,cur,maxe,sg):
        res.append((cur,sg))
        p=primes[idx]; c=cur; e=0
        while e<maxe:
            e+=1; c*=p
            if c>limit: break
            sgn=sg*(p**(e+1)-1)//(p-1)//1 if False else None
            go2=(c, sg_calc(cur,sg,p,e))
            if idx+1<len(primes): go(idx+1,c,e,go2[1])
    def sg_calc(cur,sg,p,e): return sg*((p**(e+1)-1)//(p-1))
    go(0,1,60,1); return res
L=gen(10**40)
mx=(0,0)
for c,s in L:
    if c>5040:
        r=s/(c*float(mp.log(mp.log(c))))
        if r>mx[0]: mx=(r,c)
print("ROBIN on %d extremal-form numbers up to 1e40: max ratio=%.6f at n=%d; e^g=%.6f; violations=%d"%(len(L),mx[0],mx[1],eg,sum(1 for c,s in L if c>5040 and s/(c*float(mp.log(mp.log(c))))>=eg)))
# ---- Lagarias n<=2e6 (n>=1)
H=np.cumsum(1.0/np.arange(1,N+1)); Hn=np.concatenate([[0],H])
nn=np.arange(2,N+1)
rhs=Hn[nn]+np.exp(Hn[nn])*np.log(Hn[nn])
gap=rhs-sig[nn]
print("LAGARIAS n in [2,2e6]: min(rhs-sigma)=%.4f at n=%d ; violations=%d (n=1: sigma=1 <= H1+e*log1 =1, equality)"%(gap.min(),nn[gap.argmin()],(gap<0).sum()))
# ---- Nicolas primorials
mp.mp.dps=30
K=20000
ps=list(sympy.primerange(2,sympy.prime(K)+1))
logN=0.0; lr=0.0; worst=1e9; wk=0; viol=0
for k,p in enumerate(ps,1):
    logN+=math.log(p); lr+=-math.log(1-1/p)
    ratio=math.exp(lr); rhs=eg*math.log(logN) if logN>1 else -1e9
    if logN>1:
        mrg=(ratio-rhs)
        # normalize
        if mrg<worst: worst=mrg;wk=k
        if mrg<=0: viol+=1
print("NICOLAS primorials k=1..%d: min (N/phi - e^g loglogN)=%.5f at k=%d; violations=%d"%(K,worst,wk,viol))
# ---- Franel
def franel(n):
    from math import gcd
    fr=sorted({(a,q) for q in range(1,n+1) for a in range(1,q+1) if gcd(a,q)==1}, key=lambda t:t[0]/t[1])
    m=len(fr); x=np.array([a/q for a,q in fr]); j=np.arange(1,m+1)/m
    return m,((x-j)**2).sum()
print("FRANEL n, m, sum, n*sum, n^0.9*sum:")
for n in [50,100,200,400,800,1600]:
    m,s=franel(n); print("  ",n,m,"%.3e"%s,"%.4f"%(n*s),"%.4f"%(n**0.9*s))
# ---- Li coefficients (Bombieri-Lagarias via Stieltjes series)
mp.mp.dps=80
M=60
def series_mul(a,b,M):
    c=[mp.mpf(0)]*(M+1)
    for i in range(M+1):
        if a[i]==0: continue
        for j in range(M+1-i): c[i+j]+=a[i]*b[j]
    return c
sg=[mp.stieltjes(k) for k in range(M+1)]
z=[mp.mpf(0)]*(M+1); z[0]=mp.mpf(1)
for k in range(M):
    z[k+1]=(-1)**k*sg[k]/mp.factorial(k)
f=lambda w: mp.mpf(1)/2*(1+w)*mp.pi**(-(1+w)/2)*mp.gamma((1+w)/2)
a=mp.taylor(f,0,M)
h=series_mul(a,z,M)
# log series: L=log(h): L' = h'/h
L=[mp.mpf(0)]*(M+1); L[0]=mp.log(h[0])
for k in range(1,M+1):
    s=k*h[k]
    for j in range(1,k): s-=j*L[j]*h[k-j]
    L[k]=s/(k*h[0])
lam=[]
for n_ in range(1,M+1):
    # coefficient of w^n in (1+w)^(n-1) * L(w)
    c=sum(mp.binomial(n_-1,n_-k)*L[k] for k in range(1,n_+1))
    lam.append(n_*c)
print("LI lambda_1=%s (known 0.0230957...), lambda_2=%s, min over n<=%d: %s ; all positive: %s"%(mp.nstr(lam[0],8),mp.nstr(lam[1],8),M,mp.nstr(min(lam),6),all(x>0 for x in lam)))
print("   lambda_n/(n log n /2) for n=10,30,60:",[mp.nstr(lam[k-1]/(k*mp.log(k)/2),5) for k in (10,30,60)])
