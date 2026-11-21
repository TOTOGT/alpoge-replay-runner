import mpmath as mp, math
mp.mp.dps=25
chi=[0,1,0,-1]
L=lambda s: mp.dirichlet(s,chi)
# ---- S3 twisted by Hopf rotation g=diag(i,-i): weights m=-(k-1)..(k-1) step 2, each with multiplicity k, eigenvalue k^2
def twisted(s,K=4000,g=1j):
    tot=0
    for k in range(1,K+1):
        w=sum(g**m for m in range(-(k-1),k,2)) if k<60 else None
        if w is None: w=mp.sin(k*mp.pi/2)  # closed form (verified for k<60 below)
        tot+=k*w*mp.mpf(k)**(-2*s)
    return tot
# verify closed form for small k
assert all(abs(sum(1j**m for m in range(-(k-1),k,2))-mp.sin(k*mp.pi/2))<1e-12 for k in range(1,60))
s=3
print("twisted  S3 :",mp.nstr(twisted(s,3000).real,12)," L(2s-1,chi)=",mp.nstr(L(2*s-1),12))
# lens space S3/Z4: (1/4) sum_j trace of g^j; g^j acts on weight m by i^(jm)
def lens(s,K=60):
    tot=0
    for k in range(1,K+1):
        for j in range(4):
            tot+=k*sum((1j)**(j*m) for m in range(-(k-1),k,2))*mp.mpf(k)**(-2*s)
    return tot/4
K=400
val=sum( (sum(k*sum((1j)**(j*m) for m in range(-(k-1),k,2)) for j in range(4))/4)*mp.mpf(k)**(-2*s) for k in range(1,K+1))
form=(mp.mpf(1)/2)*((1-mp.mpf(2)**(2-2*s))*mp.zeta(2*s-2)+L(2*s-1))
print("lens S3/Z4  :",mp.nstr(val.real,10)," formula:",mp.nstr(form,10))
# ---- H3: E(P,s)= sum over coprime (c,d) in Z[i]^2 modulo units(4) of (N(c)+N(d))^-s ; claim = 2(1-4^(1-s)) zeta(s-1)/L(s,chi)
from math import gcd
def gcdZi(a,b):
    while b!=0:
        q=a/b; q=complex(round(q.real),round(q.imag)); a,b=b,a-q*b
    return a
R=40
s=4; tot=0; cnt=0
pts=[complex(x,y) for x in range(-R,R+1) for y in range(-R,R+1)]
byN={}
for z in pts:
    byN.setdefault(int(z.real**2+z.imag**2),[]).append(z)
# count coprime pairs by N(c)+N(d)=n
from collections import Counter
cntn=Counter()
Rn=R*R
for nc,lc in byN.items():
    for nd,ld in byN.items():
        if nc+nd>Rn or nc+nd==0: continue
        for c in lc:
            for d in ld:
                if c==0 and d==0: continue
                g=gcdZi(c,d) if abs(c)>=abs(d) else gcdZi(d,c)
                if abs(g)**2<1.5: cntn[nc+nd]+=1
S=sum(mp.mpf(v)/mp.mpf(n)**s for n,v in cntn.items())
pred=8*(1-mp.mpf(4)**(1-s))*mp.zeta(s-1)/L(s)
print("H3 four-square: lattice coprime sum (N<=%d, /1)=%s ; 8(1-4^(1-s))zeta(s-1)/L(s,chi)=%s"%(Rn,mp.nstr(S,8),mp.nstr(pred,8)))
print("H3 zeros on spectral line Re s=1: |1-4^(1-s)|=0 at t=pi k/log2 :",[mp.nstr(math.pi*k/math.log(2),6) for k in (1,2)])
mp.mp.dps=30
print("zeta(s-1) vanishes at s=1+rho1: |zeta(rho1)|=",mp.nstr(abs(mp.zeta(mp.zetazero(1))),3), " Re(s)=",1.5)
