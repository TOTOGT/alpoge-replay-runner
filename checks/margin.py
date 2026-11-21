import math, sympy, numpy as np
eg=math.exp(0.5772156649015329)
P=list(sympy.primerange(2,120))
best={}
out=[]
def go(idx,logc,maxe,logsig):
    out.append((logc,logsig))
    if idx>=len(P): return
    p=P[idx]; lc=logc; ls=logsig; e=0
    # sigma multiplicative: sigma(p^e)=(p^(e+1)-1)/(p-1)
    prev=0.0
    while e<maxe:
        e+=1; lc=logc+e*math.log(p)
        if lc>LIM: break
        s=logsig+math.log((p**(e+1)-1)/(p-1)) if e<40 else logsig+ (e+1)*math.log(p)-math.log(p-1)
        go(idx+1,lc,e,s)
LIM=math.log(1e30)
go(0,0.0,60,0.0)
a=np.array(out); lc=a[:,0]; ls=a[:,1]
sel=lc>math.log(5040)
lc=lc[sel]; ls=ls[sel]
ratio=np.exp(ls-lc)               # sigma/n
ll=np.log(lc)                      # loglog n
marg=eg*ll-ratio
for name,w in [("sqrt(log n)",np.sqrt(lc)),("log n",lc)]:
    v=marg*w; i=v.argmin()
    print("ROBIN margin*%s: min=%.4f at log n=%.2f (n~1e%.1f); count=%d"%(name,v[i],lc[i],lc[i]/math.log(10),len(lc)))
# Nicolas margins on primorials, same normalization
logN=0.0; lr=0.0; vals=[]
for k,p in enumerate(sympy.primerange(2,sympy.prime(50000)+1),1):
    logN+=math.log(p); lr+=-math.log(1-1/p)
    if k>=3: vals.append((k,logN,(math.exp(lr)-eg*math.log(logN))*math.sqrt(logN)))
m=min(vals,key=lambda t:t[2]); print("NICOLAS (N/phi - e^g loglogN)*sqrt(logN): min=%.4f at k=%d; at k=50000: %.4f"%(m[2],m[0],vals[-1][2]))
