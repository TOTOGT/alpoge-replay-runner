import mpmath as mp
mp.mp.dps=30
chi=[0,1,0,-1]
def Lchi(s): return mp.dirichlet(s,chi)
def xi(s): return mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def sig(n,a): return sum(mp.mpf(d)**a for d in range(1,n+1) if n%d==0)
def Estar(x,y,t,N=30):
    s=mp.mpf(1)/2+1j*t
    v=xi(2*s)*y**s+xi(2*s-1)*y**(1-s)+4*mp.sqrt(y)*sum(mp.mpf(n)**(s-mp.mpf(1)/2)*sig(n,1-2*s)*mp.besselk(s-mp.mpf(1)/2,2*mp.pi*n*y)*mp.cos(2*mp.pi*n*x) for n in range(1,N))
    return v
def F(t):
    s=mp.mpf(1)/2+1j*t
    return 2*mp.pi**(-s)*mp.gamma(s)*mp.zeta(s)*Lchi(s)
print("1) identity E*(i,1/2+it) = 2 pi^-s Gamma(s) zeta(s) L(s,chi_-4):")
for t in [3,10,14.1347251417]:
    a=Estar(0,1,t); b=F(t); print("   t=%s  |diff|=%s  Im(F)=%s"%(t,mp.nstr(abs(a-b),3),mp.nstr(mp.im(b),3)))
print("2) Hessian at i at first zeta zero")
t1=mp.im(mp.zetazero(1)); 
f=lambda x,y: mp.re(Estar(x,y,t1,N=25))
fxx=mp.diff(f,(0,1),(2,0)); fyy=mp.diff(f,(0,1),(0,2)); fxy=mp.diff(f,(0,1),(1,1)); fx=mp.diff(f,(0,1),(1,0)); fy=mp.diff(f,(0,1),(0,1))
print("   f=%s fx=%s fy=%s"%(mp.nstr(f(0,1),3),mp.nstr(fx,3),mp.nstr(fy,3)))
print("   fxx=%s fyy=%s fxy=%s  (fxx+fyy)/fxx=%s"%(mp.nstr(fxx,6),mp.nstr(fyy,6),mp.nstr(fxy,3),mp.nstr((fxx+fyy)/fxx,3)))
dF=mp.diff(lambda t: mp.re(F(t)),t1); print("   F'(t1)=",mp.nstr(dF,6))
print("3) a=f_xx(i) at first 6 zeta zeros and first 8 zeros of L(chi_-4):")
zz=[mp.im(mp.zetazero(k)) for k in range(1,7)]
def Lroot(g):
    return mp.findroot(lambda t: mp.re(mp.pi**(-(0.5+1j*t))*mp.gamma(0.5+1j*t)*Lchi(0.5+1j*t)) if False else Lchi(mp.mpf(1)/2+1j*t).real*mp.cos(0)+0, g) 
# find L(chi_-4) zeros as sign changes of real function Lambda_L(t)=pi^{-s/2-? }... use completed Lambda(s)=(4/pi)^{s/2}... Gamma((s+1)/2)L(s,chi): real on line
def LamL(t):
    s=mp.mpf(1)/2+1j*t
    return (4/mp.pi)**(s/2)*mp.gamma((s+1)/2)*Lchi(s)
print("   Im part of LamL max over samples:",mp.nstr(max(abs(mp.im(LamL(t))) for t in [1,5,9,17]),3))
import numpy as np
def zeros_of(fun,T=40,dt=0.02):
    zs=[];t=0.5;prev=mp.re(fun(t))
    while t<T:
        t2=t+dt; cur=mp.re(fun(t2))
        if prev*cur<0: zs.append(mp.findroot(lambda u: mp.re(fun(u)),(t,t2),solver='anderson'))
        prev=cur;t=t2
    return zs
Lz=zeros_of(LamL,40)
Zz=zeros_of(lambda t: mp.pi**(-(mp.mpf(1)/2+1j*t)/2)*mp.gamma((mp.mpf(1)/2+1j*t)/2)*mp.zeta(mp.mpf(1)/2+1j*t),40)
print("   zeta zeros <=40:",len(Zz),[mp.nstr(z,6) for z in Zz])
print("   L(chi_-4) zeros <=40:",len(Lz),[mp.nstr(z,5) for z in Lz])
Fz=zeros_of(lambda t: F(t),40)
print("   sign changes of F(t)=E*(i,1/2+it) on (0.5,40]:",len(Fz)," common zeros (|zeta zero - L zero|<1e-6):",sum(1 for a in Zz for b in Lz if abs(a-b)<1e-6))
for k,t in enumerate(Zz[:4]+Lz[:3]):
    f=lambda x,y: mp.re(Estar(x,y,t,N=25))
    a=mp.diff(f,(0,1),(2,0)); tr=mp.diff(f,(0,1),(2,0))+mp.diff(f,(0,1),(0,2)); fxy=mp.diff(f,(0,1),(1,1))
    print("   t=%s  a=f_xx=%s  trace/a=%s  fxy/a=%s  F'(t)=%s"%(mp.nstr(t,8),mp.nstr(a,5),mp.nstr(tr/a,3),mp.nstr(fxy/a,3),mp.nstr(mp.diff(lambda u: mp.re(F(u)),t),5)))
