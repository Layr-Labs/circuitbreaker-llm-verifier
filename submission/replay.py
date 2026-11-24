"""Full candidate experiment from rational covariance. No precomputed probability tables."""
from fractions import Fraction as F
from math import isqrt
S=1<<50;H=4096
def up(n,d): return -(-n//d)
def pi():
    def atan(x):
        v=sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(20))
        return v,v+x**41/41
    a,b=atan(F(1,5));c,d=atan(F(1,239));return 16*a-4*d,16*b-4*c
PI=pi()
class Gaussian:
    __slots__=('q','dg','log','p','cum','dens','err','tail')
    def __init__(self,q,logs=26):
        if q not in (1,2) or not 12<=logs<=29:raise ValueError('CDF parameters')
        self.log=logs;self.p=1<<logs;self.q=q;self.dg=S//self.p;self.tail=1
        def rt(x):return isqrt(S*S*x.numerator//x.denominator)
        cc,cd=rt(1/(q*PI[1])),1+rt(1/(q*PI[0]));D=4*H*H*256*q
        c=[(S//2,S//2)];density=[]
        for k in range(8*H):
            z=(2*k+1)**2*S;a,b=z//D,up(z,D)
            # Enclose consecutive odd/even exp(-x²/(256q)) partial sums.
            m=n=S;r=s=S
            for j in range(1,12):
                m=m*a//(S*j);n=up(n*b,S*j)
                if j&1:r-=n;s-=m
                else:r+=m;s+=n
            s+=up(n*b,12*S)
            for _ in range(8):r=r*r//S;s=up(s*s,S)
            f,g=r*cc//S,up(s*cd,S)
            h,i=c[-1];c.append((h+f//H,i+up(g,H)));density.append((f,g))
        self.cum,self.dens=c,density
        self.err=up(8*S,12*H*H)+up(S,H*H)
    def __call__(self,n,d):
        if d<=0:raise ValueError('CDF positive denominator')
        if n<0:
            a,b=self(-n,d);return self.p-b,self.p-a
        if n>=8*d:return self.p-1,self.p
        k,l=divmod(n*H,d);a,b=self.cum[k];f,g=self.dens[k]
        a+=f*l//(d*H);b+=up(g*l,d*H)
        return max(0,(a-self.err)//self.dg),min(self.p,up(b+self.err,self.dg))

"""All-positive finite Gaussian propagation; exact nonnegative integer products.
Consecutive identical transitions are coalesced; ±1 states have analytic merges.
"""
from operator import itemgetter,mul
from math import sumprod as dot
def ce(n,d):return -((-n)//d)
def get(ids):
    # itemgetter(one) returns a scalar, so retain a zero-weight dummy entry
    return itemgetter(*ids) if len(ids)>1 else itemgetter(ids[0],ids[0])
def reducer(xs,ys,W,T,N):
    reg={};neg={};pos={};other={};posQ=0
    for x,y,w in zip(xs,ys,W):
        if not w:continue
        if y==N:
            if x!=T:neg[x]=neg.get(x,0)+w
        elif x==T and y==T:posQ+=w
        elif y==T:pos[x]=pos.get(x,0)+w
        elif x==T:other[y]=other.get(y,0)+w
        else:k=x,y;reg[k]=reg.get(k,0)+w
    def single(d):return (get(list(d)),tuple(d.values())+((0,) if len(d)==1 else ())) if d else None
    bil=(get([a for a,b in reg]),get([b for a,b in reg]),tuple(reg.values())+((0,) if len(reg)==1 else ())) if reg else None
    return bil,single(neg),single(pos),single(other),posQ
def compute_plan(qshape,m,p,ns,ps,rest,Q,K):
    qp=[Q-y for y in p];qm=[Q-x for x in m];on=[Q-x for x in ns];z=[Q-x-y for x,y in zip(ns,ps)]
    M=[];P=[];ma=M.append;pa=P.append;ceil=(1<<K)-1
    for B,N,L,R,sm in qshape:
        u=0;v=(rest+sm)*Q*Q
        if B:
            f,g,w=B;zz=g(z)
            u=dot(map(mul,w,f(m)),zz)+dot(map(mul,w,f(qp)),g(ns))
            v+=dot(map(mul,w,f(p)),zz)+dot(map(mul,w,f(qm)),g(ps))
        if N:f,w=N;u+=Q*dot(w,f(qp))
        if L:f,w=L;v+=Q*dot(w,f(qm))
        if R:f,w=R;v+=Q*dot(w,f(on))
        ma(u>>K);pa((v+ceil)>>K)
    if any(x+y>Q or x<0 or y<0 for x,y in zip(M,P)):raise ArithmeticError('tail simplex')
    return M,P
def run_layers(Q,half_mass,T=-10160,dx=800,dw=1000,dep=10,dm=35,pm=60,*,debug=False):
    if debug:import time;t=time.process_time()
    G=40000; nd=2*dm+1;nj=2*pm+1;size=nd*nj;nw=ce(3*G,dw);N=2*nw
    wg=tuple(max(0,half_mass((g+1)*dw,G)[0]-half_mass(g*dw,G)[1]) for g in range(-nw,nw));aw=tuple(2*w for w in wg[nw:]);rest=Q-sum(wg)
    if rest<0 or wg!=wg[::-1] or min(wg)<0:raise ArithmeticError('mass')
    ids=list(range(size));top=size;neg=top+1;K=2*(Q.bit_length()-1)
    if (1<<K)!=Q*Q or dep<1:raise ValueError('scale/depth')
    def ix(i,j):return top if i>=nd or j>=nj else ids[(j if j>=0 else 0)*nd+(i if i>=0 else 0)]
    ds=[[max(0,dm+ce(2337*(i-dm)*dx+1568*(g+1)*dw,3249*dx)) for g in range(-nw,nw)] for i in range(nd)]
    dh=[[max(0,dm+ce(1596*(i-dm)*dx-1148*g*dw,3249*dx)) for g in range(-nw,nw)] for i in range(nd)]
    ph=[[ce(57*r*dx-2045*g*dw,3249*dx) for g in range(-nw,nw)] for r in range(57)]
    ps=[[ce(57*r*dx+1568*(g+1)*dw,3249*dx) for g in range(-nw,nw)] for r in range(57)]
    au=[ce(28*(v+1)*dw,57*dx) for v in range(nw)];al=[ce(-28*v*dw,57*dx) for v in range(nw)]
    Am=[];Cm=[];ac=Am.append;cc=Cm.append;m=[];p=[]
    for j in range(nj):
        cut=max(0,min(N,nw+ce(57*((j-pm)*dx-T),13*dw)))
        m.extend([0]*nd);p.extend([rest+sum(wg[:cut])]*nd)
        for i in range(nd):
            y,r=divmod(41*(j-pm)+28*(i-dm),57); yp,rp=divmod(41*(i-dm)+28*(j-pm),57)
            ac(reducer([ix(i+u,j+u) for u in au],[ix(i+l,j+l) for l in al],aw,top,neg))
            cc(reducer([ix(d,jp+y+pm) for d,jp in zip(ds[i],ps[r])],[ix(d,jp+yp+pm) for d,jp in zip(dh[i][:cut],ph[rp][:cut])]+[neg]*(N-cut),wg,top,neg))
    if debug:print('setup',time.process_time()-t,flush=True)
    for it in range(1,dep):
        mm,pp=compute_plan(Am,m,p,m,p,rest,Q,K)
        nm,nl=compute_plan(Cm,m,p,mm,pp,rest,Q,K)
        m=list(map(max,m,nm));p=list(map(min,p,nl))
        if any(x+y>Q or x<0 or y<0 for x,y in zip(m,p)):raise ArithmeticError('sup simplex')
        if debug:print('layer',it+1,m[pm*nd+dm],p[pm*nd+dm],time.process_time()-t,flush=True)
    return m,p

from functools import lru_cache
from itertools import combinations_with_replacement,product
def root(m,p,dm,pm,logs,standard,half,dx=800,T=-10160,dr=800,du=800,dv=800,kd=1,kp=1):
    # Q=message probability scale, dx and all innovation widths in 1/40000.
    Q=1<<logs;G=40000;nd=2*dm+1;np=2*pm+1;limit=nd*np
    if len(m)!=len(p) or len(m)!=limit:raise ValueError('root table extent')
    if not all(0<=x<=Q and 0<=y<=Q-x for x,y in zip(m,p)):raise ValueError('root table probabilities')
    if standard.p!=Q or half.p!=Q:raise ValueError('CDF scale')
    v_end=(3*G)//dv;u_end=(4*G)//du;r_end=(4*G)//dr
    def mass(c,a,b):return max(0,c(b,G)[0]-c(a,G)[1])
    # Every omitted strip mass is explicitly assigned root/line value +1.
    V=[2*mass(half,j*dv,(j+1)*dv) for j in range(v_end)]
    U=[2*mass(standard,j*du,(j+1)*du) for j in range(u_end)]
    R=[mass(standard,j*dr,(j+1)*dr) for j in range(-r_end,r_end)]
    missV,missU,missR=Q-sum(V),Q-sum(U),Q-sum(R)
    if min(missV,missU,missR)<0:raise ValueError('mass')
    step=114*dx
    # All V bins sharing both upward-rounded child-state indices can be combined exactly.
    @lru_cache(maxsize=4096)
    def offsets(d,p):
        b=[]
        for j,w in enumerate(V):
            if not w:continue
            e=56*dv*j;hi=e+56*dv
            k=(up(d+hi,step),up(p+hi,step),up(d-e,step),up(p-e,step))
            if b and b[-1][:4]==k:b[-1]=(*k,b[-1][4]+w)
            else:b.append((*k,w))
        return tuple(b)
    sh=logs;Q2=Q*Q
    # Replace C++'s 6-category line table by cumulative tails on exactly 4 payoff groups:
    # --,-0 | -+,00 | 0+ | ++ . Missing mass is put in ++.
    if kd not in (1,2,4) or kp not in (1,2,4):raise ValueError('root lattice factors')
    def line(R,U,positive):
        dd=step//kd;pd=step//kp
        x=up(28*(R+U),dd);y=up(97*R+41*U,pd)
        return shared(x,y,positive)
    @lru_cache(maxsize=65536)
    def shared(xx,yy,positive):
        x,d=divmod(xx,kd);y,e=divmod(yy,kp)
        d*=step//kd;e*=step//kp
        x+=dm;y+=pm
        a=b=c=0
        for i,j,k,l,w in offsets(d,e):
            i+=x;j+=y;k+=x;l+=y
            if i>=nd or j>=np:d0=0;p0=Q
            else:
                x0=(j if j>0 else 0)*nd+(i if i>0 else 0);d0=m[x0];p0=p[x0]
            if k>=nd or l>=np:d1=0;p1=Q
            else:
                x1=(l if l>0 else 0)*nd+(k if k>0 else 0);d1=m[x1];p1=p[x1]
            a+=w*((Q2-p0*p1)>>sh)
            if positive:
                b+=w*((Q2-(p0*(Q-d1)+p1*(Q-d0)-p0*p1))>>sh)
                c+=w*((d0*(Q-p1)+d1*(Q-p0)-d0*d1)>>sh)
        # Compute all three upper tails independently with directed rounding; clip
        # valid upper tails against each other in the stochastic order.
        # For A/B, lower complement computation eliminates a dim-dependent correction:
        A=Q-(a>>sh)
        
        if not positive:return A
        B=Q-(b>>sh)
        C=Q-(c>>sh)
        return min(A,B,C),min(B,C),C
    whole=missR*Q
    for i,r_weight in enumerate(R):
        if not r_weight:continue
        j=i-r_end;z=(j+1)*dr;b=z>T
        val=missU*Q
        for u,w in enumerate(U):
            if not w:continue
            if b:
                A,B,C=line(z,(u+1)*du,True);D,E,F=line(z,-u*du,True)
                pay=2*(C*E+B*F-B*E)+A*(Q-F)+(Q-C)*D
            else:pay=Q*(line(z,(u+1)*du,False)+line(z,-u*du,False))
            val+=w*up(pay,2*Q)
        whole+=r_weight*up(val,Q)
    return min(Q,up(whole,Q))
def truth_controls():
    v=list(combinations_with_replacement((-1,0,1),2));g=[0,0,1,1,2,3];sg=lambda x:(x>0)-(x<0)
    for b,i,j in product((0,1),range(6),range(6)):
        z=v[i]+v[j];M=bool(b and sg(z[0]+z[1])+sg(z[2]+z[3])>0)
        e=[z[k]+sg(z[k^1]+(sg(z[k^2]+z[k^3]) if b else -1))>0 for k in range(4)]
        u,vv=g[i],g[j];pay=(int(u==3)+int(vv==3)) if not b else 2*int(u>=1 and vv>=1 and max(u,vv)>=2)+int(min(u,vv)==0 and max(u,vv)==3)
        if pay!=int(M+e[0]+e[1]>=2)+int(M+e[2]+e[3]>=2):raise ValueError('wrong root quotient')

def literal_minimum_controls():
    from itertools import product
    tr=(-1,0,1);bb=(0,1);sg=lambda x:(x>0)-(x<0)
    split=lambda a,b,c:int(not a==b==c)
    for x,y in product(tr,repeat=2):
        c=[min(split(k,i,j)-x*i-y*j for i,j in product(bb,repeat=2)) for k in bb]
        if c[0]-c[1]!=sg(x+y):raise ValueError('literal NAE recurrence')
    for b in bb:
        cand=[x for x in product(bb,repeat=5) if b or not x[4]]
        for z in product(tr,repeat=4):
            vals=[split(x[0],x[1],x[4])+split(x[2],x[3],x[4])-sum(v*w for v,w in zip(z,x)) for x in cand]
            lit=tuple(min(x[j] for x,c in zip(cand,vals) if c==min(vals)) for j in range(5))
            e=[int(z[i]+sg(z[i^1]+(sg(z[i^2]+z[i^3]) if b else -1))>0) for i in range(4)]
            m=int(bool(b and sg(z[0]+z[1])+sg(z[2]+z[3])>0))
            if lit!=(*e,m):raise ValueError('literal root minimizer')
def run(request):
    if request not in ({"phase": "terminal"}, {"phase": "energy"}):
        raise ValueError("expected one of the two declared phases")
    from pathlib import Path
    import json
    import runpy
    import sys

    directory = Path(__file__).resolve().parent
    if (directory / "files" / "moments.py").is_file():
        directory /= "files"
    gain = runpy.run_path(str(directory / "gain.py"))
    energy = runpy.run_path(str(directory / "energy.py"))
    threshold = energy["THRESHOLD"]
    grid_threshold = 40000 * threshold
    if threshold != gain["T"] or grid_threshold.denominator != 1:
        raise ValueError("threshold mismatch")
    grid_threshold = int(grid_threshold)
    ceiling = F(12685624331, 100000000000)
    coefficient = 1 + 1 / (2 * ceiling)
    bound = json.loads((directory / "bound.json").read_text(encoding="utf-8"))
    for selector in (bound["projection"], bound["guarantee"]):
        c = selector["coefficient"]
        exponent = selector["exponent"]
        if F(int(c["numerator"]), int(c["denominator"])) != coefficient:
            raise ValueError("coefficient mismatch")
        if F(int(exponent["numerator"]), int(exponent["denominator"])) != 1:
            raise ValueError("exponent mismatch")
    result = {
        "phase": request["phase"],
        "threshold": str(threshold),
        "limiting_energy_ceiling": str(ceiling),
        "coefficient": str(coefficient),
    }
    if request["phase"] == "energy":
        moments = runpy.run_path(str(directory / "moments.py"))
        rr, uu = moments["compute"]()
        directed = energy["energy"](rr, uu, moments, gain)
        if F(directed["energy_upper_1e12"], 10**12) != ceiling:
            raise ValueError("energy ceiling mismatch")
        result.update(
            directed_energy=directed,
            directed_initial_majority=energy["initial"](rr, moments, gain),
            required_other_phase="terminal",
        )
        return result

    logs = 26
    q = 1 << logs
    half, standard = Gaussian(1, logs), Gaussian(2, logs)
    dm, pm, depth = 35, 60, 22
    m, p = run_layers(q, half, T=grid_threshold, dx=800, dw=600, dep=depth,
                      dm=dm, pm=pm)
    print(f"Recomputed {depth} terminal layers.", file=sys.stderr, flush=True)
    truth_controls()
    literal_minimum_controls()
    upper = root(m, p, dm, pm, logs, standard, half, dx=800,
                 T=grid_threshold, dr=800, du=800, dv=800, kd=1, kp=1)
    if 2 * upper >= q:
        raise ValueError("terminal majority is not below one half")
    result.update(
        terminal_upper_num=upper,
        terminal_upper_den=q,
        depth=depth,
        root_truth_controls=72,
        literal_root_minima=162,
        literal_NAE_controls=9,
        required_other_phase="energy",
    )
    return result
if __name__=='__main__':
    import json,sys
    raw=sys.stdin.read() if not sys.stdin.isatty() else ''
    print(json.dumps(run(json.loads(raw) if raw.strip() else {}),sort_keys=True))
