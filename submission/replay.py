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
from operator import itemgetter,mul,le,sub
from itertools import repeat
import sys
# weights and probabilities are distinct integers.
LS=13;LQ=14;WSCALE=1<<LS;Q=1<<LQ
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
"""Exact radix-2**64 products using GMP from the declared Linux image."""
from ctypes import CDLL, POINTER, c_int, c_long, c_ulong, sizeof
from functools import lru_cache
from sys import byteorder

gmp = CDLL("libgmp.so.10")
if (sizeof(c_ulong) != 8 or byteorder != "little"
        or c_int.in_dll(gmp, "__gmp_bits_per_limb").value != 64):
    raise ArithmeticError("unexpected GMP limb representation")

limbs = POINTER(c_ulong)
multiply = getattr(gmp, "__gmpn_mul")
multiply.argtypes = [limbs, limbs, c_long, limbs, c_long]
multiply.restype = c_ulong


@lru_cache(None)
def layout(length):
    return c_ulong * length


class Kernel:
    def __init__(self, size):
        self.size = size
        self.output = layout(2 * size)()
        self.view = memoryview(self.output).cast("B")

    def product(self, x, y):
        if not x or not y:
            return 0
        # Strip complete zero digits from either end adjacent to the radix point.
        xoff = ((x & -x).bit_length() - 1) // 64
        yoff = ((y & -y).bit_length() - 1) // 64
        x >>= 64 * xoff
        y >>= 64 * yoff
        if x < y:
            x, y = y, x
        nx = (x.bit_length() + 63) // 64
        ny = (y.bit_length() + 63) // 64
        shift = (xoff + yoff) * 64
        if ny < 24:
            return (x * y) << shift
        if nx > self.size:
            raise ArithmeticError("operand exceeds the output allocation")
        a = layout(nx).from_buffer_copy(x.to_bytes(nx * 8, "little"))
        b = layout(ny).from_buffer_copy(y.to_bytes(ny * 8, "little"))
        # GMP writes nx+ny limbs into a separate allocation and requires nx >= ny.
        multiply(self.output, a, nx, b, ny)
        return int.from_bytes(self.view[:8 * (nx + ny)], "little") << shift

KERNELS = {}

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
  isaux,rows,pack,norm,B,start,cnt,st,outN,fw,gw=setup;outm=[];outp=[];both=[]
  fneg=m+[0];fpos=p+[Q]
  if not isaux:gn=n+[0];gp=u+[Q]
  KK=2*LS+40+LQ
  kernel=KERNELS.get(outN)
  if kernel is None:
   kernel=KERNELS[outN]=Kernel((outN+1)//2)
  times=kernel.product
  full_f=pack(map(lambda x:x*Q,fw))
  full_g=pack(map(lambda x:x*Q,gw))
  complement=[(Q<<KK)+(1<<KK)-1-z*Q*Q*r for z,r in zip(B,norm)]
  for row in rows:
   if isaux:ff=row
   else:ff,gg,w,anti=row
   fm=pack(map(mul,fw,ff(fneg)));fp=pack(map(mul,fw,ff(fpos)))
   if isaux:
    ac=times(fm,2*(full_f-fp)-fm)
    bc=times(fp,2*(full_f-fm)-fp)
    if plus:
     dp=ints(times(fp,fp),outN,start,cnt,st)
     both+= [(b+y*r)>>KK for y,r,b in zip(dp,norm,complement)]
   else:
    gm=pack(map(mul,w,gg(gn)))+anti;gp0=pack(map(mul,w,gg(gp)))
    gzero=full_g-gm-gp0
    ac=times(fm,gzero)+times(full_f-fp,gm)
    bc=times(fp,gzero)+times(full_f-fm,gp0)
   ma=ints(ac,outN,start,cnt,st);po=ints(bc,outN,start,cnt,st)
   outm += [(x*r)>>KK for x,r in zip(ma,norm)]
   outp += [(b+y*r)>>KK for y,r,b in zip(po,norm,complement)]
  need(all(map(le,outm,map(sub,repeat(Q),outp))),'tail simplex')
  return (outm,outp,both) if plus else (outm,outp)
 def evolve(self,k=25,stride=0):
  # Analytic finite top entrance: all children initially +; first output is 0 or +.
  h=Gaussian(1,LQ);D=self.d;Y=self.yd
  m=[0]*self.msize
  p=[h(57*(40000*(41*D*j+13*Y*i)-self.T*41*D*Y),13*40000*41*D*Y)[1] for j in range(-self.Ym,self.Ym+1) for i in range(-self.H,self.H+1)]
  for i in range(1,k):
   n,u=self.apply(self.A,m,p,m,p);mm,pp=self.apply(self.M,m,p,n,u)
   # can take meet owing to true top-origin monotonicity even though our grid representation need not be isotone.
   m=[a if a>b else b for a,b in zip(m,mm)];p=[a if a<b else b for a,b in zip(p,pp)]
  return m,p
def root(s,m,p,hh=800,hu=800):
 from fractions import Fraction as F
 n,u,b=s.apply(s.A,m,p,m,p,plus=True)
 C=Gaussian(2,26);S0=1<<26;wn=4*40000//hh;un=4*40000//hu;G=40000;K=13
 def mass(x,y):return max(0,C(y,G)[0]-C(x,G)[1])
 wr=[mass(i*hh,(i+1)*hh) for i in range(-wn,wn)];wu=[2*mass(i*hu,(i+1)*hu) for i in range(un)];mr=S0-sum(wr);mu=S0-sum(wu);need(min(wr+wu+[mr,mu])>=0,'root mass')
 need(hh==hu and s.d==s.yd,'matched physical root mesh')
 # At V = v*hh and R = h*hh, D depends only on h+v.
 # The ceiling for the second coordinate separates because d equals yd.
 den=114*G;width=2*s.Ha+1
 offset=wn+un+2
 ix=[max(-s.Ha,ceil(28*k*hh*s.d,den)) for k in range(-offset,offset+1)]
 table=[(min(bb,uu,Q-nn),min(uu,Q-nn),Q-nn) for nn,uu,bb in zip(n,u,b)]
 table.append((Q,Q,Q))
 fail=len(table)-1;tot=mr*Q
 # ceil(yd*hh*(97*h+41*v)/den) is periodic up to an integer shift.
 from math import gcd
 period=den//gcd(den,s.yd*hh*41)
 jump=period*41*s.yd*hh//den
 for k,rw in enumerate(wr):
  h=k-wn+1;q=mu*Q
  pv=[ceil(s.yd*hh*(97*h+41*v),den) for v in range(period)]
  tails=[]
  for v in range(1-un,un+1):
   z=ix[offset+h+v];div,rem=divmod(v,period)
   j=pv[rem]+jump*div-z
   tails.append(table[fail if z>s.Ha or j>s.Ya else ((j if j>=-s.Ya else -s.Ya)+s.Ya)*width+z+s.Ha])
  if h*hh>s.T:
   for (A,B,C0),(D,E,F),uw in zip(tails[un:], reversed(tails[:un]),wu):
    pay=2*(C0*E+B*F-B*E)+A*(Q-F)+(Q-C0)*D
    q+=uw*((pay+2*Q-1)//(2*Q))
  else:
   for (A,_,_),(D,_,_),uw in zip(tails[un:], reversed(tails[:un]),wu):
    q+=uw*((A+D+1)//2)
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
    beta=F(4253,100000); theta=F(329,25000); lam=F(947,20000)

    margin=F(29,500000) # certified local supporting-plane slack

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
    saving=beta+lam+margin
    delta=saving/(3+A)
    cutoff=delta
    require(delta==F(44969,1595000) and (saving-A*cutoff)/3==cutoff,'equal branch increments')
    L=1+1/(2*F(1243,10000))
    score=L+delta
    return {'score':str(score),'checks':len(checks),'aggregate_saving':str(saving),'uniform_local_margin':str(margin),'branch_intersection':str(cutoff),'semantic_increment':str(delta),'largest_coefficient_denominator':max(d for n,d in checks)}

def run(request):
    if request!={}:raise ValueError("expected empty request")
    from pathlib import Path
    from fractions import Fraction as F
    import runpy,json
    D=Path(__file__).resolve().parent
    if (D/'files'/'moments.py').is_file():D=D/'files'
    M=runpy.run_path(str(D/'moments.py'));G0=runpy.run_path(str(D/'gain.py'));E=runpy.run_path(str(D/'energy.py'))
    KERNELS.clear()
    rr,uu=M['compute']()
    e=E['energy'](rr,uu,M,G0);j=E['initial'](rr,M,G0)
    dep=42;q=Q
    conv=Transfer(180,180,151,246,157,246,260,325,-13140)
    a,b=conv.evolve(dep);truth_controls();literal_minimum_controls();u=root(conv,a,b,200,200)
    if 2*u>=q:raise ValueError('terminal not below one half')
    k=1+F(1,2)/F(1243,10000)
    if k!=F(6243,1243):raise ValueError('coefficient')
    sem=semantic()
    if F(sem['score'])!=F(910316497,180235000) or sem['checks']!=41:raise ValueError('joint coefficient or dual accounting')
    return {'directed_energy':e,'directed_initial_majority':j,'terminal_upper_num':u,
      'terminal_upper_den':q,'depth':dep,'graph_cap':'1243/10000','coefficient':sem['score'],'raw_coefficient':str(k),'semantic':sem,
       'root_truth_controls':72,'literal_root_minima':162,'literal_NAE_controls':9,'threshold':'-657/2000','scope':'all source-computed numerical premises of the joint declared theorem'}
if __name__=='__main__':
    import json,sys
    raw=sys.stdin.read() if not sys.stdin.isatty() else ''
    print(json.dumps(run(json.loads(raw) if raw.strip() else {}),sort_keys=True))
