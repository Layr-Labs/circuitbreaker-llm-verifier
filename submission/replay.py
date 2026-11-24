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
from struct import Struct
from operator import itemgetter,mul
import sys
# weights and probabilities are distinct integers.
LS=14;LQ=14;WSCALE=1<<LS;Q=1<<LQ
def need(t,name):
 if not t:raise ArithmeticError(name)
def ceil(n,d):return -((-n)//d)
def packer(N):
 F=Struct('<'+'Q'*N)
 return lambda x:int.from_bytes(F.pack(*x),'little')
def ints(k,n,start,count,step):
 ss=memoryview(k.to_bytes(n*8,'little')).cast('Q');return ss[start:start+step*count:step]
from fractions import Fraction as F
from functools import lru_cache
@lru_cache(None)
def elow(x):
 y=x/256;q=b=F(1)
 for k in range(1,12):q*=y/k;b+=q
 need(abs(y)<=F(1,2) and b>0,'Taylor lower domain')
 u=(b.numerator<<80)//b.denominator
 for k in range(8):u=(u*u)>>80
 return F(u,1<<80)
def jensenrows(d,mode,N,K,scale,normal=1<<40):
 a=F(28,57);t=F(41,57);h=F(1,d)
 def atan(z,n):return sum((-1)**i*z**(2*i+1)/(2*i+1) for i in range(n))
 piU=16*atan(F(1,5),15)-4*atan(F(1,239),14);fpi=F(564189583547756,10**15)
 if not(fpi*fpi*piU<1 and piU<F(22,7)):raise ArithmeticError('pi')
 v=h/(a if mode=='A' else 2*a*a*t);f=[];g=[]
 for i in range(-N,N+1):
  x=(i-F(1,2))*h
  ef=x*x/(2*a*a if mode=='A' else 2*a*a*t*t);eg=ef if mode=='A' else x*x/(4*a**4)
  f.append(int(scale*elow(-ef)));g.append(int(scale*elow(-eg)))
 norm=[int(normal*v*fpi*elow((i*h)**2/(a*a if mode=='A' else 2*a*a)-v*v/12)) for i in range(-K,K+1)]
 return f,g,norm

def weights(d,mode,N,K):return jensenrows(d,mode,N,K,WSCALE)
class Transfer:
 def __init__(self,d,yd,H,Ym,Ha,Ya,Nm,Na,T):
  self.d=d;self.yd=yd;self.H=H;self.Ym=Ym;self.Ha=Ha;self.Ya=Ya;self.T=T;self.Nm=Nm;self.Na=Na
  self.msize=(2*H+1)*(2*Ym+1); self.asize=(2*Ha+1)*(2*Ya+1)
  def ix(i,j,H,Y):
   if i>H or j>Y:return (2*H+1)*(2*Y+1)
   return (max(-Y,j)+Y)*(2*H+1)+max(-H,i)+H
  def setup(isaux,N,Hu,Yu):
   f,g,norm=weights(d,'A' if isaux else 'M',N,Hu)
   tc=2*N+1;outN=2*tc-1;pack=packer(tc);total=pack(f)*pack(g);start=2*N+1-(2 if isaux else 1)*Hu;st=2 if isaux else 1;cnt=2*Hu+1
   need(sys.byteorder=="little" and (2*N+1)*WSCALE*WSCALE<2**64 and 2*Q*Q*max(ints(total,outN,0,outN,1))<2**64,"no carry");B=ints(total,outN,start,cnt,st);need(all(0<=s*r<=WSCALE*WSCALE*(1<<40) for s,r in zip(B,norm)),'Gaussian lower mass')
   plan=[]
   for j in range(-Yu,Yu+1):
    od=[];pd=[];mask=[]
    for k in range(-N,N+1):
     if isaux:
      i=max(k,-H);jj=j+ceil(yd*(41*k-13*i),41*d);od.append(ix(i,jj,H,Ym))
     else:
      i=max(-H,ceil(57*k,41));jj=ceil(57*yd*(57*k-13*i)+41*41*d*j,41*57*d);od.append(ix(i,jj,H,Ym))
      z=max(-Ha,ceil(57*k,56));v0,v1=2045*57,1148*56
      jj=ceil(57*yd*(v0*k-v1*z)+28*d*v1*j,57*d*v1);pd.append(ix(z,jj,Ha,Ya))
      # positive branch iff the bin's maximum 2az=k/d can exceed -2a^2t/delta*(Y-T).
      pos=(k*13*57*57*yd*40000 > -2*28*28*41*d*(40000*j-T*yd))
      mask.append(pos)
    ff=itemgetter(*od)
    if isaux:plan.append(ff)
    else:
     # replacing G by pin (-1) is implemented by separate cached masked weight rows.
     gg=itemgetter(*pd);w=[v*b for v,b in zip(g,mask)];anti=[v*(1-b)*Q for v,b in zip(g,mask)]
     plan.append((ff,gg,w,pack(anti)))
   return (isaux,plan,pack,norm,B,start,cnt,st,outN,f,g)
  self.A=setup(True,Na,Ha,Ya);self.M=setup(False,Nm,H,Ym)
 def apply(self,setup,m,p,n,u,plus=False):
  isaux,rows,pack,norm,B,start,cnt,st,outN,fw,gw=setup;outm=[];outp=[];both=[];den3=WSCALE**3
  fneg=m+[0];fpos=p+[Q];fn=[Q-x for x in p]+[0];fq=[Q-x for x in m]+[Q]
  if isaux:
   combm=[2*(Q-y)-x for x,y in zip(m,p)]+[0];combp=[2*(Q-x)-y for x,y in zip(m,p)]+[Q]
  else:
   gn=n+[0];gp=u+[Q];gz=[Q-x-y for x,y in zip(n,u)]+[0]
  den=den3*Q;KK=2*LS+40+LQ
  for row in rows:
   if isaux:
    ff=row;ac=pack(map(mul,fw,ff(fneg)))*pack(map(mul,fw,ff(combm)));bc=pack(map(mul,fw,ff(fpos)))*pack(map(mul,fw,ff(combp)))
    if plus:
     fp=pack(map(mul,fw,ff(fpos)));dp=ints(fp*fp,outN,start,cnt,st)
     both.extend(Q-(((b*Q*Q-y)*r)>>KK) for y,r,b in zip(dp,norm,B))
   else:
    ff,gg,w,anti=row;Z=pack(map(mul,w,gg(gz)));gnn=pack(map(mul,w,gg(gn)))+anti
    ac=pack(map(mul,fw,ff(fneg)))*Z+pack(map(mul,fw,ff(fn)))*gnn
    bc=pack(map(mul,fw,ff(fpos)))*Z+pack(map(mul,fw,ff(fq)))*pack(map(mul,w,gg(gp)))
   ma=ints(ac,outN,start,cnt,st);po=ints(bc,outN,start,cnt,st)
   outm.extend((x*r)>>KK for x,r in zip(ma,norm))
   outp.extend(Q-(((b*Q*Q-y)*r)>>KK) for y,r,b in zip(po,norm,B))
  need(all(0<=x<=Q-y and y<=Q for x,y in zip(outm,outp)),'tail simplex')
  return (outm,outp,both) if plus else (outm,outp)
 def evolve(self,k=25,stride=0):
  # Analytic finite top entrance: all children initially +; first output is 0 or +.
  h=Gaussian(1,LQ);D=self.d;Y=self.yd
  m=[0]*self.msize
  p=[h(57*(40000*(41*D*j+13*Y*i)-self.T*41*D*Y),13*40000*41*D*Y)[1] for j in range(-self.Ym,self.Ym+1) for i in range(-self.H,self.H+1)]
  for i in range(1,k):
   n,u=self.apply(self.A,m,p,m,p);mm,pp=self.apply(self.M,m,p,n,u)
   # can take meet owing to true top-origin monotonicity even though our grid representation need not be isotone.
   m=[max(a,b) for a,b in zip(m,mm)];p=[min(a,b) for a,b in zip(p,pp)]
  return m,p
def root(s,m,p,hh=800,hu=800):
 from fractions import Fraction as F
 n,u,b=s.apply(s.A,m,p,m,p,plus=True)
 C=Gaussian(2,26);S0=1<<26;wn=4*40000//hh;un=4*40000//hu;G=40000;K=13
 def mass(x,y):return max(0,C(y,G)[0]-C(x,G)[1])
 wr=[mass(i*hh,(i+1)*hh) for i in range(-wn,wn)];wu=[2*mass(i*hu,(i+1)*hu) for i in range(un)];mr=S0-sum(wr);mu=S0-sum(wu);need(min(wr+wu+[mr,mu])>=0,'root mass')
 def line(R,V):
  D=28*(R+V);P=97*R+41*V;den=114*G
  ix=max(-s.Ha,ceil(D*s.d,den));iy=max(-s.Ya,ceil(s.yd*(P*s.d-ix*den),den*s.d))
  if ix>s.Ha or iy>s.Ya:return Q,Q,Q
  z=(iy+s.Ya)*(2*s.Ha+1)+ix+s.Ha
  return min(b[z],u[z],Q-n[z]),min(u[z],Q-n[z]),Q-n[z]
 tot=mr*Q
 for k,rw in enumerate(wr):
  R=(k-wn+1)*hh;q=mu*Q
  for j,uw in enumerate(wu):
   A,B,C0=line(R,(j+1)*hu);D,E,F=line(R,-j*hu)
   pay=(2*(C0*E+B*F-B*E)+A*(Q-F)+(Q-C0)*D) if R>s.T else Q*(A+D)
   q+=uw*ceil(pay,2*Q)
  tot+=rw*ceil(q,S0)
 return ceil(tot,S0)

from itertools import combinations_with_replacement,product
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
def semantic():
    from fractions import Fraction as F
    from math import lcm
    A=F(19,100); eta=F(31,100); cz=1-eta
    beta=F(847,20000); theta=F(131,10000); lam=F(947,20000)
    delta=F(281,10000); cutoff=F(563,20000) # .02815
    margin=F(1,20000) # certified local supporting-plane slack
    
    def require(test, message):
        if not test: raise ValueError(message)
    
    checks=[]
    def rational(name, x):
        require(x>margin,name)
        checks.append((name,0))
    
    from functools import cache

    @cache
    def ln12(u):
        # Only 1<=u<=2: use 2 atanh((u-1)/(u+1)) for log(u).
        require(1<=u<=2,'ln reduction')
        w=(u-1)/(u+1); ww=w*w; poww=w; low=F(0)
        for j in range(12):
            low+=2*poww/(2*j+1)
            poww*=ww
        return low,low+2*poww/(25*(1-ww))

    l2_lo,l2_hi=ln12(F(2))
    @cache
    def log2_bounds(p):
        require(p>0,'log domain')
        e=p.numerator.bit_length()-p.denominator.bit_length()
        u=p/(F(2)**e)
        if u<1:e-=1;u*=2
        lo,hi=ln12(u)
        return e+lo/l2_hi,e+hi/l2_lo

    def log_test(name, c, terms):
        # Exact rational lower enclosure, not floating or precomputed.
        lower=c-margin
        den=lcm(lower.denominator,*(w.denominator for w,p in terms))
        for w,p in terms:
            lo,hi=log2_bounds(p)
            lower+=w*(lo if w>=0 else hi)
        require(lower>0,name)
        checks.append((name,den))

    # Unclipped rank zero: cross-entropies of the reference pair and its marginals.
    # h <= entropy(XOR), entropy(U), (entropy(U)+entropy(V))/2 respectively.
    anchors=[(242781,201369,201369,354481),(208536,240277,169804,381383),(177682,207567,207567,407184)]
    for b,anchor in enumerate(anchors):
        p=[F(v,1000000) for v in anchor]
        require(sum(p)==1,'anchor normalization')
        cu=eta if b else F(1);cv=eta if b==2 else F(1)
        U=[p[0]+p[2],p[1]+p[3]];V=[p[0]+p[1],p[2]+p[3]];X=[p[0]+p[3],p[1]+p[2]];Z=[1-p[3],p[3]]
        for corner in range(4):
            u,v=corner%2,corner//2;z=u*v
            terms=[(cz,Z[z]),(cu-F(1,4),U[u]),(cv-F(1,4),V[v]),(F(1,4),p[corner]/V[v]),(F(1,4),p[corner]/U[u])]
            terms+=([(lam,X[u^v])] if b==0 else [(lam,U[u])] if b==1 else [(lam/2,U[u]),(lam/2,V[v])])
            log_test(f'j0-b{b}-corner{corner}',cz+cu+cv-theta*(2-b),terms)
        # If any term is clipped, use that term alone; h<=1-A for clipped operand.
        for name,value in [('U',cu*A-lam*(1-A)),('V',cv*A-lam*(1-A)),('Z',cz*A-lam),('I',A/2-lam)]:
            rational(f'j0-b{b}-{name}-clip',value-theta*(2-b))
    
    for tag,w,qref in [('free',F(1),F(574320,1000000)),('flag',eta,F(666885,1000000)),('same',F(1,2),F(622582,1000000))]:
        Cmax=0 if tag=='flag' else F(1,2)
        rational(f'j1-{tag}-output-clip',cz*A-beta-theta*Cmax-lam)
        for mode in ['none','input-clip']:
            q0=qref if mode=='none' else F(878202,1000000) if tag=='flag' else F(870179,1000000)
            a=w if mode=='none' else w*A
            hcoef=w+lam if mode=='none' else lam
            for q in (F(0),F(1)):
                C=q-1 if tag=='flag' else F(1,2)
                log_test(f'j1-{tag}-{mode}-end{q}',cz+a-beta-theta*C,[(cz*q/2,q0/2),(cz*(1-q/2),1-q0/2),(hcoef*q,q0),(hcoef*(1-q),1-q0)])
    
    rational('j2-clip',cz*A+3*theta/4-2*beta-lam)
    log_test('j2-unclipped',-cz+3*theta/4-2*beta-lam,[(3*cz/4,F(3))])
    weighted=(beta+lam-A*cutoff)/3-delta
    require(weighted>0 and cutoff>delta,'strict cut margin')
    L=1+1/(2*F(125838,1000000))
    score=L+delta
    return {'score':str(score),'checks':len(checks),'weighted_margin':str(weighted),'private_path_margin':str(cutoff-delta),'semantic_increment':str(delta),'largest_coefficient_denominator':max(d for n,d in checks)}

def run(request):
    if request!={}:raise ValueError("expected empty request")
    from pathlib import Path
    from fractions import Fraction as F
    import runpy,json
    D=Path(__file__).resolve().parent
    if (D/'files'/'moments.py').is_file():D=D/'files'
    M=runpy.run_path(str(D/'moments.py'));G0=runpy.run_path(str(D/'gain.py'));E=runpy.run_path(str(D/'energy.py'))
    rr,uu=M['compute']()
    e=E['energy'](rr,uu,M,G0);j=E['initial'](rr,M,G0)
    dep=22;q=Q
    conv=Transfer(125,126,100,165,105,165,173,217,-11450)
    a,b=conv.evolve(dep);truth_controls();literal_minimum_controls();u=root(conv,a,b)
    if 2*u>=q:raise ValueError('terminal not below one half')
    k=1+F(1,2)/F(125838,1000000)
    if k!=F(312919,62919):raise ValueError('coefficient')
    sem=semantic()
    if F(sem['score'])!=F(3146870239,629190000) or sem['checks']!=41:raise ValueError('joint coefficient or dual accounting')
    return {'directed_energy':e,'directed_initial_majority':j,'terminal_upper_num':u,
      'terminal_upper_den':q,'depth':dep,'graph_cap':'125838/1000000','coefficient':sem['score'],'raw_coefficient':str(k),'semantic':sem,
       'root_truth_controls':72,'literal_root_minima':162,'literal_NAE_controls':9,'threshold':'-229/800','scope':'all source-computed numerical premises of the joint declared theorem'}
if __name__=='__main__':
    import json,sys
    raw=sys.stdin.read() if not sys.stdin.isatty() else ''
    print(json.dumps(run(json.loads(raw) if raw.strip() else {}),sort_keys=True))
