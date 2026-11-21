# Generates the two SVG figures embedded in index.html (Part III).
import mpmath as mp, math, sys
mp.mp.dps=30
def L4(s): return (mp.zeta(s,mp.mpf(1)/4)-mp.zeta(s,mp.mpf(3)/4))/4**s
def F(t):
    s=mp.mpf(1)/2+1j*t
    return 2*mp.pi**(-s)*mp.gamma(s)*mp.zeta(s)*L4(s), abs(mp.gamma(s))
ts=[0.5+i*0.1 for i in range(0,396)]
vals=[];maximag=0
for t in ts:
    v,g=F(t); maximag=max(maximag,abs(mp.im(v))/abs(v)); vals.append(float(mp.re(v)/g))
sc=sum(1 for a,b in zip(vals,vals[1:]) if a*b<0)
print("max rel imag",maximag,"sign changes",sc,file=sys.stderr)
zz=[float(mp.im(mp.zetazero(k))) for k in range(1,7)]
W,H,l,r,tp,b=720,260,46,12,14,34
m=max(abs(v) for v in vals)
X=lambda t:l+(t-0.5)/39.5*(W-l-r); Y=lambda v:tp+(H-tp-b)/2*(1-v/m)
pts=" ".join(f"{X(t):.1f},{Y(v):.1f}" for t,v in zip(ts,vals))
cr=[(a_t+b_t)/2 for a_t,b_t,a,bv in zip(ts,ts[1:],vals,vals[1:]) if a*bv<0]
g=f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="E*(i,1/2+it) divided by |Gamma|, sign changes on 0.5 to 40" style="width:100%;height:auto;font:11px system-ui,sans-serif">'
g+=f'<line x1="{l}" x2="{W-r}" y1="{Y(0):.1f}" y2="{Y(0):.1f}" stroke="var(--mut)" stroke-width="1"/>'
for t in range(5,41,5):
    g+=f'<line x1="{X(t):.1f}" x2="{X(t):.1f}" y1="{H-b}" y2="{H-b+4}" stroke="var(--mut)"/><text x="{X(t):.1f}" y="{H-b+16}" text-anchor="middle" fill="var(--mut)">{t}</text>'
g+=f'<polyline points="{pts}" fill="none" stroke="var(--lk)" stroke-width="1.6"/>'
for c in cr: g+=f'<circle cx="{X(c):.1f}" cy="{Y(0):.1f}" r="3" fill="var(--bg)" stroke="var(--hd)" stroke-width="1.2"/>'
g+=f'<text x="{W/2}" y="{H-4}" text-anchor="middle" fill="var(--mut)">t  (s = 1/2 + it)</text>'
g+=f'<text x="{l+4}" y="{tp+10}" fill="var(--mut)">Re E*(i, s) / |Γ(s)|, scaled</text></svg>'
open("fig_nodal.svg","w").write(g)
# figure 2: spectral line vs zero line
rows=[("H² (modular surface)",.5,.5),("S³",None,1.25),("S³, Hopf twist",None,.75),("H³ (Picard group)",1.0,1.5),]
W2,H2=720,170;X2=lambda x:60+x/2*(W2-100)
g=f'<svg viewBox="0 0 {W2} {H2}" role="img" aria-label="Spectral line versus line of zeros" style="width:100%;height:auto;font:11px system-ui,sans-serif">'
for i,(n,sp,z) in enumerate(rows):
    y=28+i*36
    g+=f'<line x1="{X2(0):.0f}" x2="{X2(2):.0f}" y1="{y}" y2="{y}" stroke="var(--ln)"/>'
    g+=f'<text x="0" y="{y+4}" fill="var(--fg)">{n}</text>'
    if sp is not None: g+=f'<rect x="{X2(sp)-5:.0f}" y="{y-9}" width="10" height="18" fill="none" stroke="var(--hd)" stroke-width="1.6"/>'
    g+=f'<circle cx="{X2(z):.0f}" cy="{y}" r="5" fill="var(--lk)"/>'
    g+=f'<text x="{X2(z)+9:.0f}" y="{y-8}" fill="var(--mut)">zeros Re s = {z}</text>'
g+=f'<text x="{X2(0):.0f}" y="{H2-6}" fill="var(--mut)">0</text><text x="{X2(2):.0f}" y="{H2-6}" fill="var(--mut)">2</text>'
g+=f'<text x="{W2/2}" y="{H2-6}" text-anchor="middle" fill="var(--mut)">square = spectral line (where eigenvalues put s), dot = line carrying the zeros</text></svg>'
open("fig_lines.svg","w").write(g)
print(len(cr),[round(c,2) for c in cr],file=sys.stderr)
