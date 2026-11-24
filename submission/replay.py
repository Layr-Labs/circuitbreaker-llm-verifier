"""Exact directed experiment, including regeneration of the packed recurrence prefix."""
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
from bisect import bisect_left
import sys
# weights and probabilities are distinct integers.
LS=13;LQ=14;WSCALE=1<<LS;Q=1<<LQ
def need(t,name):
 if not t:raise ArithmeticError(name)
def ceil(n,d):return -((-n)//d)
def packed_arithmetic(size):
 from ctypes import CDLL, c_char_p, c_long, c_void_p, create_string_buffer
 layout=Struct('<'+'Q'*size)
 pack=lambda values:layout.pack(*values)
 width=16*size
 def python_product(a,b):
  return (int.from_bytes(a,'little')*int.from_bytes(b,'little')).to_bytes(width,'little')
 def python_merge(a,b,c,d):
  return (int.from_bytes(a,'little')*int.from_bytes(b,'little')+
          int.from_bytes(c,'little')*int.from_bytes(d,'little')).to_bytes(width,'little')
 try:
  gmp=CDLL('libgmp.so.10')
  mul_n=getattr(gmp,'__gmpn_mul_n')
  mul_n.argtypes=(c_void_p,c_char_p,c_char_p,c_long)
  mul_n.restype=None
  add_n=getattr(gmp,'__gmpn_add_n')
  add_n.argtypes=(c_void_p,c_void_p,c_void_p,c_long)
  add_n.restype=c_long
 except (OSError,AttributeError):
  return pack,python_product,python_merge
 first=create_string_buffer(width)
 second=create_string_buffer(width)
 def product(left,right):
  mul_n(first,left,right,size)
  return first.raw
 def merge(a,b,c,d):
  mul_n(first,a,b,size)
  mul_n(second,c,d,size)
  if add_n(first,first,second,2*size):
   raise ArithmeticError('packed product addition overflow')
  return first.raw
 example_left=pack([1]*size)
 example_right=pack([2]*size)
 if (product(example_left,example_right)!=python_product(example_left,example_right) or
     merge(example_left,example_right,example_right,example_right)!=python_merge(example_left,example_right,example_right,example_right)):
  raise ArithmeticError('pinned GMP limb layout')
 return pack,product,merge
def ints(raw,n,start,count,step):
 return memoryview(raw)[:n*8].cast('Q')[start:start+step*count:step]
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
def jensenrows(d,mode,N,K,scale,substeps,normal=1<<40):
 a=F(28,57);t=F(41,57);h=F(1,d)
 def atan(z,n):return sum((-1)**i*z**(2*i+1)/(2*i+1) for i in range(n))
 piU=16*atan(F(1,5),15)-4*atan(F(1,239),14);fpi=F(564189583547756,10**15)
 if not(fpi*fpi*piU<1 and piU<F(22,7)):raise ArithmeticError('pi')
 v=h/(a if mode=='A' else 2*a*a*t);f=[];g=[]
 for i in range(-N,N+1):
  x=(i-F(1,2))*h
  ef=x*x/(2*a*a if mode=='A' else 2*a*a*t*t);eg=ef if mode=='A' else x*x/(4*a**4)
  f.append(int(scale*elow(-ef)));g.append(int(scale*elow(-eg)))
 norm=[int(normal*v*fpi*elow((i*substeps*h)**2/(a*a if mode=='A' else 2*a*a)-v*v/12)) for i in range(-K,K+1)]
 return f,g,norm

def weights(d,mode,N,K,substeps):return jensenrows(d,mode,N,K,WSCALE,substeps)
class Transfer:
 def __init__(self,d,yd,H,Ym,Ha,Ya,Nm,Na,T,substeps=2):
  self.d=d;self.yd=yd;self.H=H;self.Ym=Ym;self.Ha=Ha;self.Ya=Ya;self.T=T;self.Nm=Nm;self.Na=Na
  cd=d*substeps
  self.msize=(2*H+1)*(2*Ym+1); self.asize=(2*Ha+1)*(2*Ya+1)
  def ix(i,j,H,Y):
   if i>H or j>Y:return (2*H+1)*(2*Y+1)
   return (max(-Y,j)+Y)*(2*H+1)+max(-H,i)+H
  pin_cdf=Gaussian(1,26)
  def setup(isaux,N,Hu,Yu):
   f,g,norm=weights(cd,'A' if isaux else 'M',N,Hu,substeps)
   tc=2*N+1;outN=2*tc-1;pack,prod,merge=packed_arithmetic(tc);total=prod(pack(f),pack(g));st=substeps*(2 if isaux else 1);start=2*N+1-st*Hu;cnt=2*Hu+1
   need(sys.byteorder=="little" and (2*N+1)*WSCALE*WSCALE<2**64 and Q*Q*max(ints(total,outN,0,outN,1))<2**64,"no carry");B=ints(total,outN,start,cnt,st);need(all(0<=s*r<=WSCALE*WSCALE*(1<<40) for s,r in zip(B,norm)),'Gaussian lower mass')
   plan=[]
   full=int.from_bytes(pack(c*Q for c in f),'little')
   for j in range(-Yu,Yu+1):
    od=[];pd=[];mask=[]
    for k in range(-N,N+1):
     if isaux:
      i=max(ceil(d*k,cd),-H);jj=j+ceil(yd*(41*d*k-13*cd*i),41*d*cd);od.append(ix(i,jj,H,Ym))
     else:
      i=max(-H,ceil(57*d*k,41*cd));jj=ceil(57*yd*(57*d*k-13*cd*i)+41*41*d*cd*j,41*57*d*cd);od.append(ix(i,jj,H,Ym))
      z=max(-Ha,ceil(57*d*k,56*cd));v0,v1=2045*57,1148*56
      jj=ceil(57*yd*(v0*d*k-v1*cd*z)+28*d*cd*v1*j,57*d*cd*v1);pd.append(ix(z,jj,Ha,Ya))
      # positive branch iff the bin's maximum 2az=k/cd can exceed -2a^2t/delta*(Y-T).
      pos=(k*13*57*57*yd*40000 > -2*28*28*41*cd*(40000*j-T*yd))
      mask.append(pos)
    ff=itemgetter(*od)
    if isaux:plan.append(ff)
    else:
     # replacing G by pin (-1) is implemented by separate cached masked weight rows.
     active_at=bisect_left(mask,True)
     g_active=g[active_at:]
     active_columns=pd[active_at:]
     active_layout=Struct('<'+'Q'*len(active_columns))
     active_pad=bytes(8*active_at)
     if len(active_columns)>1:
      active_read=itemgetter(*active_columns)
      def active_pack(table, read=active_read, weights=g_active, layout=active_layout, pad=active_pad):
       return pad+layout.pack(*map(mul,weights,read(table)))
     elif active_columns:
      def active_pack(table, col=active_columns[0], weight=g_active[0], layout=active_layout, pad=active_pad):
       return pad+layout.pack(weight*table[col])
     else:
      def active_pack(table, pad=active_pad):
       return pad
     w=[0]*active_at+g_active
     anti=[v*Q for v in g[:active_at]]+[0]*len(g_active)
     # Correct the pinned fraction of the unique threshold-straddling bin.
     ell=F(2*28*28*41,57**3); shift=(F(T,40000)-F(j,yd))/F(13,57); vcrit=shift*ell
     vbin=ceil((vcrit*cd).numerator,(vcrit*cd).denominator)
     corrections=[]
     if -N<=vbin<=N and F(vbin-1,cd)<vcrit<F(vbin,cd):
      end=F(vbin-1,cd)/ell; low_den=41*d*shift.denominator; high_den=41*d*end.denominator
      low_shift=41*d*shift.numerator; high_shift=41*d*end.numerator
      for parent in range(-Hu,Hu+1):
       sb=substeps*parent+1-vbin
       if not -N<=sb<=N:continue
       lower=57*parent*shift.denominator-low_shift; upper=57*parent*end.denominator-high_shift
       mass=max(0,pin_cdf(upper,high_den)[0]-pin_cdf(lower,low_den)[1])
       if mass:corrections.append((parent+Hu,od[sb+N],pd[vbin+N],mass))
     on=int.from_bytes(pack(v*Q for v in w),"little");off=int.from_bytes(pack(anti),"little")
     plan.append((ff,active_pack,on,off,corrections))
   total_weighted=[base*Q*Q*factor for base,factor in zip(B,norm)]
   return (isaux,plan,pack,prod,merge,full,norm,total_weighted,start,cnt,st,outN,f,g)
  self.A=setup(True,Na,Ha,Ya);self.M=setup(False,Nm,H,Ym)
 def apply(self,setup,m,p,n,u,plus=False):
  isaux,rows,pack,prod,merge,full,norm,total_weighted,start,cnt,st,outN,fw,gw=setup;outm=[];outp=[];both=[]
  fneg=m+[0];fpos=p+[Q]
  if not isaux:
   gn=n+[0];gp=u+[Q]
  KK=2*LS+40+LQ
  width=8*len(fw)
  def tobytes(v):
   return v.to_bytes(width,'little')
  for row in rows:
   ff=row if isaux else row[0]
   fm=pack(map(mul,fw,ff(fneg))); fp=pack(map(mul,fw,ff(fpos)))
   fm_i=int.from_bytes(fm,'little'); fp_i=int.from_bytes(fp,'little')
   if isaux:
    ac=prod(fm,tobytes(2*(full-fp_i)-fm_i))
    bc=prod(fp,tobytes(2*(full-fm_i)-fp_i))
    if plus:
     dp=ints(prod(fp,fp),outN,start,cnt,st)
     both.extend([Q-((total-y*r)>>KK) for y,r,total in zip(dp,norm,total_weighted)])
   else:
    _,active_pack,on,off,corrections=row
    gm=active_pack(gn); gp0=active_pack(gp)
    gm_i=int.from_bytes(gm,'little');gp_i=int.from_bytes(gp0,'little')
    Z=tobytes(on-gm_i-gp_i)
    ac=merge(fm,Z,tobytes(full-fp_i),tobytes(gm_i+off))
    bc=merge(fp,Z,tobytes(full-fm_i),gp0)
   ma=ints(ac,outN,start,cnt,st);po=ints(bc,outN,start,cnt,st)
   outm.extend([(x*r)>>KK for x,r in zip(ma,norm)])
   outp.extend([Q-((total-y*r)>>KK) for y,r,total in zip(po,norm,total_weighted)])
   if not isaux:
    startrow=len(outm)-cnt
    for at,fi,gi,mass in corrections:
     fm,fp=fneg[fi],fpos[fi]; gm,gp0=gn[gi],gp[gi]; gzero=Q-gm-gp0
     dx=(Q-fp)*(Q-gm)-fm*gzero
     dy=fp*gzero+(Q-fm)*gp0
     outm[startrow+at]+=(dx*mass)>>(26+LQ)
     outp[startrow+at]-=(dy*mass)>>(26+LQ)
  need(all(0<=x<=Q-y and y<=Q for x,y in zip(outm,outp)),'tail simplex')
  return (outm,outp,both) if plus else (outm,outp)
 def evolve(self,k=25):
  # Analytic finite top entrance: all children initially +; first output is 0 or +.
  h=Gaussian(1,LQ);D=self.d;Y=self.yd
  m=[0]*self.msize
  p=[h(57*(40000*(41*D*j+13*Y*i)-self.T*41*D*Y),13*40000*41*D*Y)[1] for j in range(-self.Ym,self.Ym+1) for i in range(-self.H,self.H+1)]
  return self.advance(m,p,k-1)
 def advance(self,m,p,layers):
  for _ in range(layers):
   n,u=self.apply(self.A,m,p,m,p);mm,pp=self.apply(self.M,m,p,n,u)
   # can take meet owing to true top-origin monotonicity even though our grid representation need not be isotone.
   m=[a if a>b else b for a,b in zip(m,mm)];p=[a if a<b else b for a,b in zip(p,pp)]
  return m,p
def root(s,m,p,hh=800,hu=800):
 from fractions import Fraction as F
 n,u,b=s.apply(s.A,m,p,m,p,plus=True)
 C=Gaussian(2,26);S0=1<<26;wn=4*40000//hh;un=4*40000//hu;G=40000;K=13
 def mass(x,y):return max(0,C(y,G)[0]-C(x,G)[1])
 xs=sorted(set([i*hh for i in range(-wn,wn+1)]+[s.T]));wr=[mass(a,b) for a,b in zip(xs,xs[1:])];wu=[2*mass(i*hu,(i+1)*hu) for i in range(un)];mr=S0-sum(wr);mu=S0-sum(wu);need(min(wr+wu+[mr,mu])>=0,'root mass')
 def line(R,V):
  D=28*(R+V);P=97*R+41*V;den=114*G
  ix=max(-s.Ha,ceil(D*s.d,den));iy=max(-s.Ya,ceil(s.yd*(P*s.d-ix*den),den*s.d))
  if ix>s.Ha or iy>s.Ya:return Q,Q,Q
  z=(iy+s.Ya)*(2*s.Ha+1)+ix+s.Ha
  return min(b[z],u[z],Q-n[z]),min(u[z],Q-n[z]),Q-n[z]
 tot=mr*Q
 for k,rw in enumerate(wr):
  R=xs[k+1];q=mu*Q
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
    from functools import cache
    from itertools import product

    A, eta = F(7, 25), F(7, 20)
    beta, theta, lam, xi = F(757, 10000), F(3888, 100000), F(197, 10000), F(195, 10000)
    delta, cutoff, margin = F(2909, 100000), F(29091, 10**6), F(1, 20000)
    checked = []

    @cache
    def ln_unit(u):
        if not 1 <= u <= 2:
            raise ValueError('ln reduction')
        w = (u - 1) / (u + 1)
        squared = w * w
        power = w
        low = F(0)
        for j in range(12):
            low += 2 * power / (2 * j + 1)
            power *= squared
        return low, low + 2 * power / (25 * (1 - squared))

    ln2_lo, ln2_hi = ln_unit(F(2))

    @cache
    def log2_bounds(p):
        if p <= 0:
            raise ValueError('log domain')
        power = p.numerator.bit_length() - p.denominator.bit_length()
        scaled = p / F(2)**power
        if scaled < 1:
            power -= 1
            scaled *= 2
        lo, hi = ln_unit(scaled)
        return power + lo / ln2_hi, power + hi / ln2_lo

    def exact_test(name, const, terms=()):
        bound = const
        for c, arg in terms:
            lo, hi = log2_bounds(arg)
            bound += c * (lo if c >= 0 else hi)
        if bound <= margin:
            raise ValueError(f'{name}: insufficient exact supporting margin {float(bound):.8f}')
        checked.append((name, bound))

    # Each row is indexed by (input flags, output flag).
    p_anchors = [
        (211483,211483,211483,365551), (220332,220332,220332,339004),
        (173146,258410,173146,395298), (186718,258715,186718,367849),
        (160605,208421,208421,422553), (172567,215366,215366,396701),
    ]
    for b, o in product(range(3), range(2)):
        cu, cv = eta if b else F(1), eta if b==2 else F(1)
        cz = 1 - eta * o
        du, dv = (F(1, 2), F(1,2)) if b==0 else (F(0), F(1)) if b==1 else (eta, eta)
        k = du + dv
        offset = theta*(2 - b) + xi*(b - o)
        p = [F(n, 10**6) for n in p_anchors[2*b+o]]
        if sum(p) != 1:
            raise ValueError('anchor mass')
        U = [p[0]+p[2], p[1]+p[3]]
        V = [p[0]+p[1], p[2]+p[3]]
        X = [p[0]+p[3], p[1]+p[2]]
        Z = [1-p[3], p[3]]
        for index in range(4):
            u, v = index%2, index//2
            terms = [(cz,Z[u*v]), (cu-du,U[u]), (cv-dv,V[v]),
                     (du,p[index]/V[v]), (dv,p[index]/U[u])]
            if b==0:
                terms.extend((lam/3,q) for q in (U[u], V[v], X[u^v]))
            elif b==1:
                terms.append((lam,U[u]))
            else:
                terms.extend(((lam/2,U[u]),(lam/2,V[v])))
            exact_test(f'rank0-{b}-{o}-corner{index}', cz+cu+cv-offset, terms)
        for name, value in {
             'input-U':cu*A-lam*(1-A), 'input-V':cv*A-lam*(1-A),
             'output':cz*A-lam, 'info':k*A-lam,
        }.items():
            exact_test(f'rank0-{b}-{o}-clip-{name}',value-offset)

    # With no flag, the two ideal rank-one shapes have the same bound.
    q_anchors = [(601866,944084), (572318,924026), (708066,962640), (664959,951103)]
    for b,o in product(range(2),range(2)):
        cz = 1-eta*o
        w = eta if b else F(1)
        offset = beta + xi*(b-o)
        Cmax = F(0) if b else F(1,2)
        exact_test(f'rank1-{b}-{o}-clip-output',cz*A-lam-theta*Cmax-offset)
        for clip, ref in enumerate(q_anchors[2*b+o]):
            q0 = F(ref,10**6)
            coeff = lam + (0 if clip else w)
            base = cz + (w*A if clip else w)
            for q in (F(0), F(1)):
                C = q-1 if b else F(1,2)
                terms=[(cz*q/2,q0/2),(cz*(1-q/2),1-q0/2),
                       (coeff*q,q0),(coeff*(1-q),1-q0)]
                exact_test(f'rank1-{b}-{o}-branch{clip}-end{q}',base-theta*C-offset,terms)
    for o in range(2):
        cz = 1-eta*o
        offset = 2*beta - 3*theta/4 + lam - xi*o
        exact_test(f'rank2-{o}-clip',cz*A-offset)
        exact_test(f'rank2-{o}-unclipped',-cz-offset,[(3*cz/4,F(3))])

    saving = beta+lam+margin
    weighted_margin = (saving - A*cutoff)/3 - delta
    if weighted_margin <= 0 or cutoff <= delta:
        raise ValueError('strict cut margin')
    cap = F(30937568299,250000000000)
    coefficient = 1 + 1/(2*cap) + delta
    return {
        'checks':len(checked), 'semantic_increment':str(delta),
        'weighted_margin':str(weighted_margin), 'private_path_margin':str(cutoff-delta),
        'coefficient':str(coefficient),
    }

PACKED_PREFIX_FIRST = (
    '{Wp48S^xk9=GL@E0stWa761SMbT8$j=Em!ooLvBR58L<v)msKhyWR|sr&WY&dd*_fj+;{xa0p0}76~n_^PZkR(ewCR(1_a^ignLr'
    '*?BWhts4We`;cbs;O(coQ-Zxb8w-|cmmF&)qz2hjg2KI$d&bhlxIWJG54s6)NrJLu-H9mty~SF?PeO4}4HE#k)T!@JS>6dE=`M0v'
    '1m+>bCren~9_w}s4&i^uMB?ck13Hj*+$9k&A1=TVguskfEXOn!M@?+Q;RPC|9{@{OB-tZpSKxpsMN-AI>+!`u6fdgi?Z$S~Ug-nP'
    '6oKFbErF*85{tp8FEMvTq6BZCSbb!SWS!#L339i*Jk#zg#(}h*bTfw7@n-lmDVyN(5FG({50<{9!rxCIS3%IO_Q<C`_{LBt^44YC'
    '4%uFnD-Fwwn&vFNcq7`dyli}7Y_K310tiVsWxH9x;M0>=orG}TjG}DsJ?D`S(CI(PXXM%a9)dX_T9!s;X&JVS7ic@5zlt8Z!Vmle'
    'Ac+FtGF$DTnuY|h+A#2BQ2T?PO@TL4I;>Q@pYv1M%&5pqC7H;YsSJ0hD4DLjd~4w(7Z-m%0R1QFL|m!$Vmll-Zc5XLpMrVNLcpzJ'
    '#h}J+k?aZt(!ys5;_>Ew_k(!FU%Fd5W{#@1xT<Gh(h6GFZXe^;ohXR9tk(lYT-|`I>IdLY-UG(+PW0{`(PdN<rQHB&Y}>~qO$X0{'
    'vxh^(uB57(k;#bf!;~BFvAO4Y>B=-wVzqh6(3ESIcjg$5MGPdnPDQq3OgVb|m-B*o2`K_1mH6bL<A{*Gh&Zch9UJ>mQK?`+0U=gk'
    'X@kCsk^SVwRVHAA)#_oSQtOZY1|)5b*9nnxEqoaL$!<!Fo!>q&of@mVM`U==8lrE>-jg}i09e;TG2UIR-4D-nh)uNSOmrhoMi!c|'
    'vwunUw#qIlJ-Ep!td&F`hqzhgBZ?eP5OCi|_88PfCTl7I;-|r-ibt^2q_8EDGa+#<RuxsQ#^)$5s%OBolup;5gGxPl7o)QsOTNXt'
    'JqCL5NY+lB0Ae4OH``$t)0#r^B_8scd8mGSzeD3isAb)B<f<BpV7ar72PpD5dB&Qqb~s|xbd3T`shv0zhVpcbEm%dW3)o(u0?(do'
    '6m;0(vfiu<aAYa#{Yo~iNKj=pQCRuxE8&B~g7v)f$~2Ibr#V&sblTA8{$V?XFETsSnD@4<`yOw>B2^7~>W~*wnl3ryk7t-O=zzv)'
    't89as1bM1C&8&kPIF>UFrSca2!C@A;+xbS(Qmlbb7la0jf&ih9t`P#~@po8)K-&Y1qTSIVhgH9br|HL|m2*a*&ZP`qaoe5R+G+6Y'
    'n+iilhUYG_#p=Rlo+{k>HO^x7qguR*3KKyv<&y4Yws`8u&WSN`y0u+Tg|+Y)3B>;n3ikfmD$t(x%4O<SYZm??QTHMw6jsX83;uvO'
    'Mm`-Rn~e(>bA4-c78P)`l*%yg5oVAt-=FvEyR-feFl(_<6xSZ0nH&XhroURuMAd3I$hzP7X(vp`AbvRvOQ*PD56P0rE6g7=)_<By'
    'UMZ{{tY3pbz{iYwzRKB5d_IDE2V=Y)8p(69c!m$cA)|={40i!SPVQreu->E0kRX-;G`uRRZmTjkoN$nwT6nlI4z&zZ1HQwMKTK+k'
    '^ygR^gxq`o)h$<ZUO8>gw;r{^ThimJiWQ854fiYDknKex11`-ZGCr$yY#Q0jU`;n&SeR<vi^(c(nbiSZd^@zi6N^PT6L1nYF{sp8'
    'Fr#?{T}UHz`m#z;6#{JT8Q0^nIZP=D#AKFrSBuc`A6r_5zQ`0&Sn;S~oc{PMGLp6aIAjFT(qS0x4j{GQ=7bm-aD<a6rw49${~P{U'
    'I3$|}CRH^ze0%FZKRidWR+epS+F2S3dn$F2yK`b7+f-hE<$lcfAwQX<5_Qn#DDjtuT3$EwzEQ^h`_Q|vZnrV@4D(BBgOTrWG6-hl'
    'O--U4gm6&wOO^4YrZ`4J4u+HAP7Oh5T*9b?%(cxM+}Zc|KwQE&2blbkOls#FjwL||DDeH8x|bI3-~%Zp7Ob$ze?R{u@E7rCI^}`q'
    'hOOkKX2CZ%zvI&`bvWMy!04$;;e687;D*+Be+S}fD*NjjMXxPl1{++o;}#&8_0nnXrh!9Gz%P(fBn%T2*Y%Xh!XI2LZ~PENLa&<>'
    'wQfRgj(a_mL(;KLTP-jl@el*?i*9FENDHp|Di?HA9t<Lz4*s1??K0q=&8kk=HW9tOdU!SN=DAFBdd^)rk>@j^8ffm#H$ciA)RIo2'
    '6d$F3!rlSt-*%TaNCh@?VLC`o99a4KvBYzVu@S%lcE;J@-;tFaxJ247uV-DPTZT$m+CX1rdzVR&$Puo>Kzz5LG!gF-lrNyMBH_c-'
    '^St??snQo+vWzP?>eTc}Y82+_pz9Wrfr)Qc^aXKlVPc4cC!bMOrGTbQ%v)#ik?tpl@7!C8eRgF7rRU<W0gt#{InRfG_v`vN^322e'
    'e%N}M-NbJLx1h(i(ZDI?2&>E<)P%L?V@fbg0V=_TEI56O)%F}ku&Kf)?8Dulfi}F16}?ncE5yg>$n*NtH~I$Cm$i#@kvioDW`}^?'
    'g6#HxSaV$v4_IUg2)@vcJi@s65YNfI<nxcH;B6-Tt^3;6Az&KTT1J{n!8A_}`-Bt)AE9h8T(}D<uO~@7=|aI7Bo;JM!&-Bkr6AZK'
    '_xcp6R`;cf27%u-p)a}zVvqp6PdBZ;Z~@?wAM3oFB&GW-;h8AF=3riVgjc@mQO9)8l|T|mTU}}xZ{na!@Ie<<%%9%F%)=|uRn~=b'
    'T$OdOtq^Fn+Yfi#%#4wB%5We<0n*RBie2X-ArA)(uWU}wlgMkN24E||UGt}40yBSLWMRyG{<vH)>dD&SOK&@O;zAA@J@*}Bq!##B'
    '!XiI-68p^-1>nb2-H9+zsg4NmkebdNBrv+*(AJX466pmkv3UY=+E2+i^sW`b@>p|VY`*+e55|>SvMV@!p{;R>!GrIX7{j0tw)x^C'
    'w3jX~1&t9VC0|&6&L~zw6k|M>^BPj9so=Nwqt*QS6>tvjIc|yA!|qY!pWa&-mN%9h1@*`kIG#~%O?_hr1o5WF49L82LwWgtSnR!3'
    '5~;GqYlLyqD>0T9BswUp$Y-QFOLdzeZyGmSas$f`^+tY(@(KW$5)fi7wajwJ;rz;B$Eeb#9O1cCT%z)RL=R&m4w?`HOAdu-CKP-q'
    '*;`iN(mV5-TV7Q?P4N-~$uh<F69>HFQAF3aXRxtoi5npqQrk~AUn^`-e%9mth@?a)jEy(|QYDg`h>}F#;v_?;V3T4P9Y*e^=V=?g'
    'S~$fT$MzXjmF7*uKPZuVJ)``o2q7hjqXG%MQ0R!t#VK<*T=H;91%$b6qEToI=^A5+v!vawWV{1-p*ebBE%@Tp@@VXUTBfPSe#0%4'
    ';6(?EJ~xp3V=>O`7yOP-_)ccIrP+S~)D(3Qj=zenkVVXSiIM(kZG{L1ci)p3r@$13W?b9D_X&3PweA`rX`MdevW6e@it~mO=@Jw='
    'nqYT*!?*CmOELCY<dN5rc(nV_Dr6~*?jG2O(2p{sQJ_a<MtIRKHquSMv~5xXgi}3g;uB1wzhd9XxU1{aPZ~1Y<m4w{Cx!U_#sCgX'
    'fJXD8ORXYJE=Ki}Q60JmL>1IlhG2s5S-ISf4Qqv9R^;dy{6VvgESwUvP&MljwMnX2bin@dYSnt!MEIi{8k_zzNp!_%^=i0_r)r1#'
    'ar(Ans9lTti7BC0(ZAk#GOH7XWm>-n;NBh)SXunYGM!U)hJSea5#F#wmmpNLl0KgHCc}6gWg-N%$mKg>FAEO1w-><~Y#D5(7PjAP'
    'a^qaJk>iA%V4DB%e@${fdBk{)i-vcdm$K+NGd(@j@gJEF_|_c`^G`uL;-q^(>r#J6vd?9OgNxq4I}STzPe<u-H_o!rg@rC-0;b2&'
    'VpDV4yTr96>I$ell3J;048Kx_at=#y^vbcgB`x^ZEDMb8GW0phsja}Cq4Rap14~SS1xht7mbnz6rMoLW0l!O~CZaqen~y5#<v(fG'
    '8WbwG3BFOA1zDFV)F>GKpbcGYqQJ+l&@NqX3%npm>#JR{38ypkkd~y!K|V>6Uh;QiUZ}UePVRe*%e;N<;J?^eR2XM7{j%)G`f*58'
    'B6Fyj>h${#O4Wd#x70?Pw+O(z>6uYCRrB_LX}S&+M(_F2LyUvH5m3UT$C-Yyb7^wuE9veGLS&#Epn>4*3<#8JB)%A<BoMOg5ydW+'
    'qXKc~J-ceKi>xw$822^Qi@v8eml%I#zdcR80z4PHM!T|eBi^l@+$^Q*0FTa-ADvx|@$}c|=Zim)^d)|!@mD-j*)e^3yK1W0<lVMu'
    'a_Y%OICT7iAGbWg+w4+yrYegi_pXGbTYHC$jX<;bKLK(&H8X3L2H63TdtGt1JUX_@PDOo8!TA|Z!4}zHQG-{%gT#J|9EFNE=CHuT'
    '_b!8CMTz0N1C}S+%g&OV3!SBN+Z8H~QS!6!%#htLv%O~Km%pSR61)PHz_9)^uhwV*bz)pN5I5OSA)cHp8nHcIoHcOq5~DC7Cna&;'
    'IN(l0a_cr@v!P<@?Bzls%k(?1&x7#HehplJ5C=ouZeoVG@7hY0r2g3RX!6;$^>q{X!q3mcfXS++Xp%0-T(bczFa9~xXLRIA6&P|V'
    'h%%$rE39A*j9`<Esa0Gyx*x0)6AF)yR8I_4RaD{gg0lkd<_~PKvB9<;Wf8!gu{7kXV9hJP9F%Wwz02~jUMk|x3kSEcsut<$(w{NR'
    '3M?_(Z{2|2N5!O@A1%YhM0%g4Oh2pXw~6xAx0<Bf-2H;%9PYsn7o&^y6ZG>}yz`)cSU0%>RL~&~T3=%(p6S--XSpP5rgZKRmZrg_'
    'uXynOqZDD!>u_5bx=rxZa(}gUl^y{O5rB*l3?Mxr8mbRXVj@>@@tm=R$>?>rInG4Nj<H3=C+l4ujLZa9UVoVv$l&TOxl7MN>2IX;'
    'Ka<L0s@*KdW<-9l-Ods7k!DnI_c&xz{H}ffXfMWejkKwaMS6m^2|?f>Ab<>GqBH^4bLlZi1;=Mx+U`v#M73786?b(mrAW&d7^{g>'
    '`M^$HcE?igk*K3S!Lz*qO=G9*%8UYuSA1scJDBcwZ!mRu*sio++r&mv-E=mHe1fEtwlUjv8by;p(T>la4<QgYklL+Z#8XR+(C16?'
    '-y7;f+oGfn-Z>=wvzg>haamb@uCFNQtp~%pz%wrmhW_J!jq;-Pd%&nQ6qZigUPXcBJ?IW;eR0f%qS{8im<9s`R@g3<#Lh`HV`cF}'
    '6IF0D8tKMNmv8m2;MER}ALIdk7|`oa>r_427bd$V#GkR#vNsD@UgRG_FJBK?nM$Ww=o1cFf4B44r6uOhGuB{A)9!*8pjTOD%fYx*'
    '4Wq{ew0K6pjnlLb$cYU{m5ykksDGASFv<T1KA|L?%2h^<F#H$l_KaO262eSnbm<ai4M(719i7kxaw{0C{+X3K145*FS^59nivfg%'
    's!CzuCVSQS&`ZphXZ4_5)W72hTJ@Nr^$;gNM8&l@r2$T(#gweQ4RY&5yQP2@1D0ew#;PR&En66<_KkOFX7g3Thgg%ciQaYit;2O7'
    'Fe=4YTA8QUMK&5<O8rK?`#bA`x;g*_)ZM~BO_B3RBo8C?YG0roSEWn)iWtk{wJkhquaxMooGaX%B?25eJB|S(Uj;4FMb_3%f=duC'
    '?Lx){2GR5tikWY;2WUbA4J1W`x<ut}^XZ6yO<$lKnU|RJR-15^$K%EBzRaS2en%&`5NVtsN}5*Uj$MI>SlrEtYs@K^VG7UZ5$#UT'
    'h-|T}BtV!$P~hC+Yw6Su!)AW!rTSG_`zTuEIN^cHHf%gy=vp^5`G6FyG=o^67H^E^F0B+>5#px?07l7gA#Et7z}J+$QLtcK7j5$y'
    'd{tq-y`}8^Z#R|uiz>1nqIO5jxz0z)28fR~#Gt2z5QH96<jyr`Xr`FCNfIXEA%>#}Rn;k5IA{ZeT`0vhhQE}F2JU9noKDX6xw(1s'
    'u$%d-{iIcLjQ$&0+WptlH8bUn71XwQANb(CXbl!W96T*~#;;r7L<N{*D)O~x>9nfbB$}jjf4rudBawR$W~cg?#C*a1<*rLN|27>*'
    'kG9sk-iBcSTxh$Ml0!gO%nb~99_zIQlhur~Ddu(!N_a1vr3_qyrW!wUL>rXk&q(r{R2@P=EOJ$8Y-%8O8qtSI4b7yMh(c7lcRCv$'
    '>zfeohiYLD*<KWu_?I6k*c)P%BNNv@dJ4*}h{gWEffxRXI;09fo>K@dD(($3wsaM4MyX*;KiM!@$8cdMAVH<pyuZ!W;tzbtl@eSi'
    '5kA)1e~;>#Unr<&%2a~Uid_-^P?p27SUNG^8EBg8MAvvwc}APeGN!S2q{NzWMm<S{Lz^5_*|E7v{3fK4RWE7{trkU{LcxrD>rT}v'
    '^UtGU!~7L&myOw_>7-g3<btg(Ft_-ruSVy0i(tP#S&hHd>dSQsnbF?~d62slDY!Q4j4I(W<CZC?^3nOual6v8E9gYwe3xWC15={R'
    'rHcvLq=S=_6$;?uB2`JW|2bi(?b;2&;2-E_>qHT2*D72K)L-TI3whlBc{tU<+IIKiLsh7*pOEwy%*(54rd`Duj=|U0DZ8%%P^Y$Q'
    'hU-~KoVw(a#_{Q?+%}u5XSuG=Cjb~>Dd5E~h$<OOD+Kb-8YT}hBjk^wqTG_e@by3V{XXTY-Rvb&^VsF?{S%K!KEfruPKZqp`A}Jw'
    'UBRJsgO0I!Def%4i)uc%eKk2y!9IX=EP=;YD>1C*z6u>!H==ST7?9#s_akNp>?(gl<5EHe2Os3uy~HPWe_vtP6PB3deJRQL0MT!M'
    'pgE!PHl>gDJ-5XRy87EMuVyf9`-8Q9rKSAu#59{&0dW`AN_hT-(A>iwDBVCQs-<HXn^T+=AsIiTO$mMa0Q(Yd0=GdnPbwChSkg>@'
    'pd1PCul5H7t2cQ(UuI|9iq7C4DZ52N=L8HWrp-;HXk#K<xN*Bo3Jbi=WWEKs&6z#UxL`SYHP3POp@K>tNBK^_MxAX5Em|W|J|TqI'
    'AE*10BkCRUBc%YHPVq=j0tP+Ejv}|Za~Q_V;E!0es^DpU64lI_cyhS<T!8oOo!8<9?3k8C@OD-*A^F_8IAJDvRqIVdK_COB9>u^L'
    ')xfOmb_u7|`xsIC(Ce!2p}p`??yp&bL;j^FR8rc|Yg4aFD+FtqqgrUhB~~LbFZo${4A)sap$Ao7kCtg<Ze-;CA7j^`B?cX@fP+?<'
    'aC!4K<7iU!lcZbyXQmS(F;w;_6R}@U911@cW4ivkX-WFUHFB2lD_d5B>q|7LjLABq^gq7lkvp8v4uXDzg}&{l63%PXAWjjwS8esd'
    'r^LBuiO)EAY&AoGI(H+qVyvVpfCP(J`?rxs;S(9uhH=+3i=vZ&Zts8pBcAn10cFsUWCn&0ScZ@>JPEs=<<S6ARE?4=9^yCx<9e`@'
    'cynH2bgtftbSv;eet4(-@>}s8)m)B{hfOS;e)TYf_6jjXM4(V15(s2O-qo5<B%`f^BA~F`D1*&PSlQ_(*Xl{j5fREU@xwDVp&LVO'
    'A|Jqu@rG6GaA`C|xc_whGqi^jZ{dJ3rl<C4J73a1F#Ug#lCW$7bGzxoXuxFsb<<vwVZpD{FH(=`&tbDqKK&Yagg^8K5Hwc^!82M+'
    'z%QYM$>>=qZW(Z0dh1l_1me<DNhh>+*y6yUT_tr0-~rfbR{`;T3Ag^V&T=rrPhzVv1FpUYr3PtRw$iI}fLnUmWWmE9iR|dBV)O@Y'
    'JkftrW5<@}RiK3-!|`IHec8&67LC*=1@-uzrj*kR!*IF$K+3u7u+8He<mDYw4iEDUdb4-XiGXZvWhIW1Z|Vt>nM1J}Z06bIJ&Ag&'
    'J=h1b?tw!4UWr_+4^1w{ex@7sz!b_Fx6$A)8Z)J_cE*e2zhjPHotcA;^X;~aATeiT*(V7`#QZz$dy65`?DLbIiU|<$Ng}3+1~D99'
    '7bw{ktzLF*&4?x>V*9V%^8hiY9XuC$(9`KO*`XX|Ao8^s%0bh|4x7{5KphUv$On<shS%G|YN~#Wj#HO8f<(K`mS{P_hxpUx?jove'
    '+V0*|Sa`7UvJifX2skO{XA-SAGXiAx?ijQB)PuH)*ZnP0fW!jAxPBW~SQUn<h0o{6o$!n=Er-tdRF_#8CBM|-J7FvZ^F$t1?&i??'
    'ditqRvY;&Ymodpel$ql<-6FDu2BW6d^Yz)f+QUH~0oWCXVMCn&d4Ki0&bAE}qw^}&LAhjHqEFGma_bvI+9BM}A&{u43Cds!UgFgh'
    'nd;8M%yrKidXFcz2Yt0VsAbxKo%Qx|bBKH{k4z&_B$C3>6d%Ks>9I5pAfqu=$bDM^)g<B(6=Zl5m=A<MPPIy4lAei2!h9I1kWF*T'
    '%@LQ}iQ1Q0QE1_8qbvs&Eb01;9zqqe(H1|CVnfU(d?EhyVh_D2TozghLG_;@GY*VEZ}y9t%Yf#=u<ccoaB!2?hlZWkV{!mDT(fAV'
    '2~ofU6Aa5hQic~kbF+zqs1$WH;IAgLVO-NXw(hX#(*Xit1U{UL#$eAV-upF>E7kCrrzDxt4vyOB@!(4NnEMq8#UL981DF7q@)_K5'
    '0=q0;^bwN5GLiE|b*GwCX`x?t@?x5dz7c=YScv@{+8xJ(Ny}-E?0`!3e?)v?FbQ+i`t6ZVVTp-K4I#~&!v@)m^EZO&radbIDwzAt'
    '<d|pPwQi(P3ysn4YliVDzQd#(C_LNIOcR~e^vpPoo+6?;ZEL3Y#;+aZE%dMhZjXJvao8?zU}N%z+>lVH-04$@-mPH`93xa^+y<tv'
    'VohUMNx`T`)ZLTHVDw0fJQj|-Bn*&M3e_f}v<I`g3dDhCk#Aku^BEi|lej>Qc8*Rkm53bbfRXw;>ge{QG{`)!4kiV|E6%0|)f0(A'
    'h~1Y3Ri|eh>dP}tbJE)?9>DLCQXz9TfLI8turiDeYhU!*)`?6WvJce~_VwIn%ChsAvh&dO9HW(-wjMSD`79I?v~@c>wh>F)1f@vn'
    'f*L_2FmbtZoo(O^=J9x%LT0K|0$JHoa=SPig4Iy$g~dTm)~MdEHhMmU?p~i2bDm*M4O3!(ACKJ~hq`LJ`Iw&Wm0DVKd8`S$ghUL4'
    '`8q&RVY4L-7N5pF0uA8Z%Z5Fz%Wa6+M*X5{Wr%<CiVAnv6<gO!31p@Zg8&H+1iy0O6&_xFE9!3icB2t(sHKHK)e0u4&Y+Z>vC?Kr'
    'wSztAYBKa0PqK!Y_-}bq0L2`%it)N*2JEE9o7fzR|F1u2sWu8)RUhYez<mwmv_zkU00hjkj)0`l62xLErM{W#B;+C<@iLZH-s7CF'
    'fgL|qYl~*%n#KG$QVf@XK*4EiaNK(!Z$ohI<MU3s(FN%tsm&~N^}5j!W!)<mm9L!J#Do6Os?Z%1Y&Y6LlQod;xhKqO1E{^3lX~}2'
    'VXRndizGX5i|O1w597z+uq)%dM>nQ5w067@k(?mf>g@S`1H7%rC2pw=3S1Ou+`h`7vj=JEWiH<J5PWJOV2e_RRwoNa+))2<I;03o'
    '&95;5c*)eBK&s$D*erq-VA+M@U;JA?GkFLnEbZ;#&1i$2Z?9ks5XiZ7;mWk-+FW@gol<oeRpsI3*XKI;en~Y2&n}CD&^DXnYyIAJ'
    '8D-}*Vzg??FiueJDbz9ZB(b{-KMY$x(;jGHBvpgqO8ARI6@CKg<?i@yA0|vZkVXQj!qiFwFPRU#|K?nEf)N~y&hG$Ux%X8!_>kWY'
    'dVdcpDmGt|7&wdn`&Uq_AsTtBt;jPt0J&~Rs{AEo!c<PK7M{VK1dAKhoat=ADA+N`#liY~pWR}IXJ@jG7TY;#u2+_n=~Prj`OL+8'
    '>B%7|4_ms(ULo<yS1{keKs|kxmZ<ZZ6__2Q&}zq%;nw`2PPlp~Yy|-PFy{^9o<l4hJ2;rkHXsB{$v}WSA1(r}Yo{e7EXE27K-q$i'
    'lcCBGAx&R^GEFudImB-cYV{|k6lKBSgv}!^Ht4YQzmn+h?_h0(<kg+;tll!8Q(yL&zu%V3AzC?Fb9&?Rs$T>WmK#+|@U==u+Nw)B'
    'q=lUve?@#swbF6<l<~trH5o>QbsjrQX}rLiV8OYHWwoNVa^}C4Un|E>e|YO3S@-rOS-XJehD20wnR;>EvsrB<A4643WD$*sG7sm%'
    'HnjyCiefRc^Sw@a&&g~~hQ?<a5HdauM7UaR=&Jw=zRNZ5NrB4gtw_)_Jypi;MV*Nn*L|1_1Ks$4>CD%>!qRTaqof-7M5_>8i;8aK'
    'rndGlQ9O4+2JVFXKPImAzdIxh)=78do=sE7OxfN0A@i;T2>cGv=0T|RuGZ4ect5XWKQpZsh`$<Kt_`EB@kB3YwrF}r&@t<VHLZdL'
    'QZ-4B=LA4U7u!bI;A*fx;Niyz{@K(2$0#<nHO2Mom0kx?iw5%&TB#q7D~&Mx)lI39sMDM7g8RYK_@F30tcL^jU8-}k=1Z6+MDDHt'
    '1!Wh>Qi~_gi0ur?-?w{ss|>k2e{~EMB*^!9v`6nC0>a%~l)F(37IK>7RbbH3cEJp}Y6xa<RA~@q<CD8C-00n1fOwT&2b*QhQ6j93'
    '`QeP5aC4kdlyF6^5;e38^ZNRb5WuaKnHm=QlPQ7$O{-LhqsQhky%={Ma?ywWs<IP=@cK-<(J3eiye=4gcoB<d@%RC5envHVUX*f2'
    'rZNU{G9VVwuoY1DfMZz2O`xk=Fk|8YnkJ&4I#W?r=QSoon>bpjiWb4>NE)ExdHSwF(d@r@vksm`P%(Ib0=H-MMD@8I_+Dkzvprbv'
    '_tXt60=+dn4XXDH(itiRGiS}1(yFmRXRLV?>Z@FHk3!9I0*5Wp2_zF~(^^x*>Q&Kz%<b#;i)}(S$XkHW$(aVv1Rb{h;lAu+DL~3G'
    'G8AYp)SuljejD@)*d&ET7UW@Sn4Yk>QQrfr*rAdmj3}lGmuW5ft87Ra7XBPqGJ&-9EK={`<FTdKNQLs$W-S6J-!SPMGuk4j?DkgJ'
    'CsHMwqvh|+vdON14^OncmAxLEw8eq`>;YT)ZZiz;^)vX;j<g0X*2@yWIdj53D0tbDOS!Ya5OTlhC$T#(7+f3`_L8m-r|mNN)lb2R'
    '_xA?6^RV?Qkl@PqYD{g_sQ*Eo4yG!m;T1F7nFlW<K)Au*7S4~n0t+L}la@WC{hz?=A^vXA!8l6(D(wd&ZsOD$iOHKxrX8Ifj9fV}'
    'mju0Jna-8CrAo};7Sd`#48~e%SP(a}G_v!X>27GcsS3j{{#=7AewsMLpT%CD3kim&UW!%*b180v@6vf#FH{|=pLsJs4}_|Bv?~aJ'
    '>?e~A;3I56X<XX9o0NH68+j1nkWG%OvLw#HNng>qIay>&v5tuY(qT(xY{$j(AVO@Tu}@`a_cs;+e<#wn4e_|Er)mGPqWb(oQ`RDO'
    'wR-Ne#*UAyn<LA(0+IeZLTHxGqHy}QTcJc-bR?k8<DF#(do|975O7`-TaU1d!ya}g^^vm8MxbMnQQ6Zt|Ma}tNCrZ4;uKuvjSG?^'
    '`b(Gbcziz+WXylpnXgM9LyUj%U$PI0M-TZ`t^=tupWkd64E<=PG(%eFqRp(kLX3oEEo9=g<~b_nQ78s%nYJ0G3if?e-A1MRd68S&'
    '_PEcr?vgFc);K*F@|KfGzw4#E2(A83jfu4Yi0FuCgMEV#`7d-8aC7b3lF3^Us=gvB*SK0XAdily^&r7~dfpV$VW2bx0H>vQ^e$)0'
    'B73f}KpA2s&B-G(mN~q3qWZPwK@4tBVKtbphilpx_zt~}BzykKf^<I;3akke73z70y^cJYr(P2>!h>%j4V`R>hx%2=hk64>V)FBj'
    '&h+}Xy9Mr4tk`cq^|4~AdnS*G3b|XQoFWgzH1YUUuk;WHklBt$O{R&o^*C;D7G3wFKyz^NK6}lT?G4>G6x|Xrp6Y#7-(G4YW7ncT'
    'H5O>p?v=OH+qL0L9?iPNJ(Nf*UVRLQqGuHZ*h3fQ{9N(BIyvp~&U;q>phaTteR1Q=aS`G!{?PRF)@Ta1(yT_cQK2J@2z_wSyywdv'
    'D_gZUWo;7*XIV?3+E|3`10%&Dv{;)wv+N++s_HWXMv+>`Ljdy;2LMAEUL&%NrJ?Lv(g5?OZ@(7`_oH3SmeV-bey1UGWDkJEHzhTd'
    'YAt=`@NLHOzVEz=eLO5dIqb3y`JUWQ^lV$j5Q|;6cpEdUSar5uwWByuoBg7@k;IXeb;QK>6NbTt?-)lQ?UjRQFtMp3Xz=HKA6rBi'
    '$LpHM{(CA5EsPL%b*FKIX9}gRM3mZlIoW1yr}pCr1psvut#t6UZC5E0=pP=8wasjx0J8!Kv4_GC^xHkaa|anR>}R)w@!T)O#ErOk'
    '$3aRTHyjn;pa2t-hg1!kor^d&;N~CjA+5qB+Ws!CAT-3bzz3+#@@vJQkX1D$lMR#y4;?`L1wT8eZ6`B(h~VCF*Y~+Y`1rYFQy%a@'
    'lp}%_&h6v}H-=F)(d7;BTBQRnswMI()m3<Yos2qJd1@HQ;e^WBg|@Q=xy9Ua@9FMCt<rtn4`+G57MJW1jx=bK!ia&908d%MN1`ge'
    'Rn+pEK8?&^wzqBXsfQf!cy1i&z*-n!gD+-bS3r)6EzX+^Y!Rs?rff9sbq*kxOwv8H<2{&Z3~+4>Q?D&-&G=;HZ}47FGEMIypy?J7'
    '7j+r56)mpk$1bCZ&T#}Di-&BlicWi5mYr{^a4fD7(^54;AHW!_0OFhe11}^QS-^N^s$>BVZB5vE*aG;=UR~f^Q-@WmFaLM=S?m3g'
    'h^34cW5nce{gnr3w<_+T&dIcRtLM^5W{gJgR3KpN%F3hJKe@DBJQvW{WMauTWCrGPT?d^krOlq312STa&R8$jng@hS*nie6U>6KR'
    'G6RaZRHl{Ut^dN=!#PwTIu|CdJmv*D-xm34FQ+bm^>EAnnBK@MKjq?9X^|G&)aqPZ2VTKqw*Wxu-I#w+{iSp9F|qgM(^HSEGDtkR'
    '3e#-SmT-RqNc%A?1s3imBE9t_SG$jY&3CC0;yOq5x6KnO#TsLv&5u!nMIWVCS^n#`)Fbo^!hlf}O`;5e)0J}nX!_U;wS#Nf<Xs}s'
    'H?5tF&}a5VE2(@rz+G2#JoGa17asJhQ&y)P;7fH`Ne*7I!RBS)m@Z6h`9-ysm^|G=zOGn|3G&KnZiEH@TcMbgAf%-GX4C#<C(2^X'
    'Ho5l~qzj7nmJCyjf?2uSs42>c%j-9ssVo8QDBBI5RMtyuaHu%c4742y+3tB%BUpJSDKmH8EXk0R<GW(}KNiE#y#H;ha-c=4DTm8M'
    '1pi64!%?peCike0`NOzWbgb@8v@jnZj?ME<702J<aSFei_5Hf*@38o|VTg`W&pkswU4|h4Jjlt;F%OF_dH}o=&I7*w&=rN^r(&^2'
    '?Jj3#g?ybTR)VZaS$;G2ESi(dYHYJp(B*f!Da(ig-N8bG4BK!N`no*_B#(#6-J$iaB@L+Y(~TFJZcT3KSqyGA4Lm|}`u61zlb;+u'
    'uw98<0i0&7$ihkqi5xrMYtBKvsF%)*3E>vFZ`HK0rC>bsX%=)<zxd7uB#bEOpCIxQ<V(~UhkUH(uYIg4J4r~Hc)sxe0@ICdSJ<71'
    '>{lLh1;kPe0_OG`gE0#qHzAXv>OK97V{ddIA%%^xBz~C^o1jm5<*mqJ0*8AbVy@vCR<?nnyEH77G%qik=~A4izYCvk_b#8zu0##8'
    'Lig+>8RUXf-3!kqBykPF$_w+CUO*J+oPkLE7NkY8oCrfRr4@||VLy4!88X5O_5W`O9B&6lsk94J4WabTxdND`c}Pi+AeQE-1_=0%'
    'q*RrupuK?as6wXvhZ{xT1V6EwKKGl({T*V({xu6!<ul8C@UKQ}#s+f{_4KqA41`tjv<#9KC_+`beoHHl4DCse-&%<N8c{<_pc;$V'
    'a0OM;8tPq+u6wi4{AG&UUS>91Rn4I-{o3O04DX#X>uEN7l5j>)pzKCXS7nQe+w;3O(8K>0PUZL(hzgi{AbLmqA+xEUSWJs?dibu2'
    'Swa-~Hl!&8oXLXjS5iBD<s24|)<jMpl6D%IjKPT@R$FH?bnZLlc1Yn!iD6wvs$hwsesh(Pw>c%9l@NJ6pw?RXR~d#VVaO1T)tOaz'
    '<a?x*f&hZ{(4&Llp``HYHt0r34ufh6C(6<B?b{@6hsxxgvnwj2$=(HmJR>k(sf)d6^#FZponl@6rJO1Y>`E)g?@HI3uBn&v^K8f<'
    'J;Jl{k^v_JgaWg}S?N7Jm7%l9kw+XvE<Lz*%jbq_stK@Tr)|jS{oF*Zki`&ABY36x^ZmQ|g>tir$^G!!oXd!attiU{!yfA37(R;a'
    'BHnqa><v9zH2)9qnvDb0Dt1lW!vODMa8Ol!luA*1CK#pKLO-l!{wR-mQPQ#38;vl^6+5ISWETilc94~*iBb)<GgKOa*kvu<BLSKd'
    'WRulLSFO8p2;%|W$3CL5(h{oYZ#bt3kAK;bOYVLxj-nf!BNQpuMpR$k8-Mmotb6$QfeW9P?Ic+}kTvJTyfa1reyIASJE=18UGvuH'
    'ev9CqM&{6QzsL$2%_ejgOop}oVW*;)e|^-FC182kH_NBc2ps=6nIE0g_JLEw7~CkI8F$%PeAaSPzd2@)fYxoq!9`B55|o-Ffv7b='
    'kTa^O69msQTU8)?1`C0gp?eETN4mZfROK99vr;b!d<-L|iwX+AL8iG~skAGoVpojg1;Z^P@SULC44sel?L#i}vha_hB`_5GZuT$z'
    'uM2aVWvwDKMFGpT_034%Ua;zP8SQ6|jr!t)$9%_R7`lP4K$O_E4i8i{p~I5046EZBT}MpafhoU)i~u`pzaER^Rz)i*u16v-e)U42'
    'ItlRps$LOSgnrawYN=Q5HCyjf$WeI}?GdqYco$^TJ~Ta5^hM)w!^Sb?NnSz^zw!MDiIj4nGPJ@9V37wp;vi87JJ^kWQ*Nb*Gbd!>'
    'n;;^Ix<0nH=NP(J9l>!L*v>YB2ya6$TZy}kbBcP&Ld9{vE5JksN2d^daqSiyk`aYFs-y}F@=Jc1hQO9lI}lf(uxDY*9`@z2(ui$`'
    '_Z@eOm4cw3n>NzW@NDBDwZJPo(m__j`rqdc4F}gaQ{t%PCIq(z<&;8?*?x=<q;V1u*KoOXD=}UWDwaOBMvW#qPUakZK6B1fm8|)B'
    'AwTC=_nehiW;5G1>(8*;{HW)7j8fmz3c1^F9F$kF&kL7k5!)7J!Wx@FkSn!7g!6(F8oUvG(sNpjKnh7Qx5zEa=ktLjk&?P1U!O;j'
    'hlbtJ6_%KH#{=>{B$PKCfg-p3rwe`AENMFk-9}&m7;&26*^f_CDO3*6gzdu_ZhW!5xT$2dMn;wd&+7&avL-y{mXsU704A`Pi<0P5'
    ';n@r<X}J-KK?txvqF2qus4#7tXf4|s?7IE=GDfMXAGK{~!$7j49?`y1L6f>zyDfwC#+uG`q{U(YV9JDStXcw7-dNsJJ|Wg!&!)xC'
    'K~7FP5D3aEAxFyZ%&{aKOBmK6<=~&1mI!mIA@26#6{mPyvHZ;()mIm)G3MM)<>V8NQWo<Stu&6A+O!yHz0z0`oQoy;&U?HdJUZ%1'
    'J$-FzaN;wNuDa@@e0P{NP>KJeMs&pXGjW0}IOT3g!hOo%8K8jXwm2HSnEOFX_$Wz3@FL%&A~+?`N4D^3yJvA1Py#fR@i{MW_nftS'
    'JYaL(5?q}VuJdK+WYzEUP?mK-D>D~P$B?agKVfQo&CXc!jW9+)D;eOTMt(hR!>SAFr2Bl5A{)64uOfh&I1J4=%_?MSk}AV?xH@jw'
    'Bl4!!z9Dw{Tuz66W~-4M4-L~Kzd)W7g(%3<;JLV9miOGl6eWCn{k4cgRYoNh>97iko{OW{2eCc!DlfAby0MitrsB2}WHEe3@7m^5'
    'itFj>ikJd%z)vpTK6mP(tv<AXE`cUxBX_ESX?cv!Vc3ahx+e)_%J#Z}dbhN-<!IA1YjQ7nsxS=_J#E)gbueUu9*H?bt!KjtJH)uI'
    'Q0Ixj;b}|ptKQ&PCC|e|8L>+ritD94GKJLccYNJ$plGaJ)QbP3o{C)s%fOv?IQfoodvL)u!sKd@$iB%qd(g$ih!W>rXFp{~Xda1s'
    'Xp9`maYwL3R?kDGGrw|m559U6eD|=3vTclUNz+JlfY#!r-!)Z?`wyD+8cmO65o6Yj3f0A?C0hXM(MG3IeuHYg9PC`8AStG|z%Z+A'
    '3fTQo_FbB7{>Vl3@gD08$*@GHFlSs@-7FxH3RL1F1^;P*IQENhnc$`fnr>vQ>OR)a2_WBuQt$dx(0%v@2=Lpm4bDTl?;%xLMz~$_'
    '{+YYPYpQWrF4>{u61VQYga6ayK@#jINAbX}Xmp}pLn(8pz<<J%&K6p+TqCfJ1@~^=99o^jyC9V#TZc>z#JD{*IipBK=A3S4cqp>%'
    '(ku`p=O+=jlvx(u7*rSts<;jOik{)zC~${WS?iI)33fEs)bMfUdl$O~514`T3wF3Zk=%btQmDC-Lf`7P#;Eh%|NOMrhxhA6Kd#x_'
    'Jzaw-jzdeFSn|^)5?OaDrkG6FR3JXE##~88H^}4ir4-R`eKo}Kc**0=Bucmn;xH(Ssa5A{n@MjwHg2dBD%S_U^8Dco$09k<@~p(i'
    'tU$o-y<U69V#LuQ-rh}bAhaM|PEEV^ulP4bNOq4wyls#2^7Gj}N~kAZwZO~FabI34#&Wcd3Z1*l;DwHv-%dB0Fg>^KJwzI&pIz1W'
    'V^c>Rgt|=xU<u{0i_B$0oJFaH9&ax;Gp-#Dd6kFFC|M&hlE%C&xSuc1yQVmD)ls{HG>G+8<GpQBp*~Vl9&><R{8fW??8p9b_?Sm8'
    '#Scugl+P84(8-p_zJB3BQ)dt(`kf<DJO|?y5=n#y)F%BGO@gBoojeLrTAl*>YEEzBDo-l&GItT3xF2Mj2DnSou3<RfE%Zta+nS-c'
    'i)|$=U$ohGMAw;V<S(;~`|4TMQzMa0MaPu~OF0cai34*%uro3skeJwZ*15TkD6j(%q(8go02kah4Dr`ER+zsEL*a?ewz_Papr0NA'
    '<_q9E9eC5c#ZPX@&9CREQLAT^sN&6_fr{pQq0W7mVJ0#p=CPqf*MwF-Z&yAprh|CqU*B-qrq7zHeU#D|g!&+=-VojCHvR%DWsZjU'
    'Y=t@ZIJk!s@O;3&`M<bmbH5|84@W@F&=Iqn=Y`+3+!M|s;9r>&0lPcBaskscWAP;-Z8R!>VsM+uEJ-I>RFrK&Px0`w^O93K{6Prz'
    '8A9Yqwdr&|C=U+IhtlRomi~#Q{7K>7ta(SPiI(mUQ87qrTk*BY9?XBmwcP??sh+Z2!#Tv_P15#tZl{T!0_(tQX^BAXa(r|G<M{)+'
    '^LQ;m9=LYWu?*%3WuUeWQY$!AYfrxvgrOv;fv5L@WzgW6B2w-*#kB*?VhBB}W6!hpegqD>96R^#xyhD54FRlSz}Efzs~BZ=R*sik'
    '&0?=oEbcjA<*lGzsXk2?rD%ZAR!6SI&KmRx70-?%Z^A>S*6n=Jct#ackfN!wD7^klW@3M(95gk5j<i2wE|t=PG-d5Pb6n2i#-ciX'
    '2|4JvDG|5X&*P1f(3NRBt%-0wKc>jA-?o$;w{LgLH+{ho-(TUZ;i@@Z`j3Bw1+PLsH6VFXdiDXVYV+Z`r<smBNjMv<3)%Qwc0w+T'
    '=#4dK4ro>iceLe}R@Gf{JhPpuaJoTUZ{Se!c7nHW%MP*=Ku2OxqT^9P%Ag}z%vIzPAt#Fxg=aF{nT)Sv^=+>jPsQ2XEgdgh!O#4y'
    'm@kxRluK@IUy1zyaM_ODEqD+bNr(bdfbW57+H5wj5#Anjl2xdHC`Mmg^U!?$Q|MchomfB-Jr%%Ksy<N-zH=4_fFJpb>n&aCt21o^'
    'P^!JIjodn_equp`Qi2vfKu6r1(rndS%;>iw-Aa^Z(;T?I;W8Z)J%ywCLV;{Dr5!}HnVz;FrLZyrg=+E#%Xif;#M?F|_GTi@zpq{s'
    'VA?k`x85JaG(Y=lNK5q^Gr-p|{6{FGfaemu5nbq`;1_NxVLZGl?$Wc~J~ar&X#fmSh#HdnW&Fn>ax}#%DG@x6CJUAQDTzgIA^sF6'
    'Rm~g00AcS_#i(?FTQA>(I4mgfk#H7?m5Fcxfr9%%peq}O*obyAu?w&0j{frND>arj$=|&66n{1oHk?wM-l~Ptx*(`@hPFa+NQMUb'
    'kt+}s`AwL;tLx)#$nxa#t=#QbW~(I7^MEPtwy39_finP6lDZjk+(LO#ypUG4oGQgvD$Z$ZCvFOxZCF<j5r)PhnR!9^FAkFxZq1VS'
    'M-ZWf)z~VPltC0u23?aqhel8r`zJhfR$0cnT-AC<-@rK%n-=h=ftZ}@rPIx-MsY&&p<I1{i<NT%QXD`2gHo*YSRpk%ZVDc<Wt6I`'
    'XkVzdR!LcH89jKv#|p#^I9f&!v2qWs_u@NF&_$~28@oA?wY^g{TDJO-{+UkrF7Pa}HjzNm?ci7aF|00QN*~7}!TwV<DDUl4cMMR_'
    'yk$(hp@7{f&@S!0v@61xzjSuS@EuO|GWaDFd5i=0d-%*h;jhSoyao-?;^FLczTia6c;WqVh(56|h{jce#7nutt(ur%d&3qeCqLCU'
    '3OIV#&O2%R9rd$+KA^Hi_Yb<CLu+fzQ?BDVVHhjOi0{x3<8MrG?a@&29S8jVV|p%sp@bqV)#u{oAnUKYuVZ(3ifr6DzXp1WCLa~#'
    'gi%~aXyr5tq2MXFPE-T`R^pn#Uaz7W1%L$1Yb8ko$k#0M8aWY-Qkq)_MUc3p@>7*1Fy0URTbhx?vWxl^ksXV-(0x0A*-g3TfM?M)'
    'XFPbI_1XKm^f5C~eK2b3lDJ0CzX4(^lN1=VN|g(gHNfL)^j76C(oxU*X<;Q?=g2-WFP@K_a-EC~x@}GwNow^b(=ihw%Y`qr9A=e!'
    'iN5m(oUc8}FA!>oy2)|fF!G+TNal&FePF}JxCyw_CRz`WR<d;`4_&S8Vt3r^LpBvCl6|+m=OV;x@}u~r`VrI0d93VRDNUU^<St*Z'
    'E>WFgr&D?RUR&>q$hb6Mx+qGx_ZGs+ZJRTm{!%oN9VGqz3%hnWur$0Up<j@D0KFSfCh+%ifV7XM;3U-k<qvCWfztp?TQsWlfzM<r'
    '1d}wF!WBP@;r;Ygj~%_L!>Ju*=db2RZ0q}L?HZc(2Jw`Z$Y%er(s;277G*>SBl36y3vPh`bI8SOE4ZU9K-SAXOsH|s5T_2oVu%|~'
    'XN+^sW+jj@J+ElBP9IXqA<x~Rt|%Cwu+h?UTG~*RizM-!B(S=*`#f07&I9ImpIXvQ@Vp9gu&RJ6%cUkOB8h<#O7XgKGTev1Q_ASI'
    'Ixv^of@bWwRH|xy6$JM%iesiiWsA@6fmFKoX<7%`nLn{P?}tnNf&$<WjI%3Z6HWzD-J8mfbaWkD!anw;Y8d+daob<&w0^NpG(6ZR'
    'wby>_na(}frYZy)fp|>JQ=H1lKXAEyyyz`ibE_+|OjeT#M!sEbLF|7xb0Hqag-9k@C7orp&zV<@Tm^p8*7iVfCwysIrbV?Jx}**~'
    '-?{<-%c=J6pw#i8A+-1@C9kGaP$b11=*S4>`_Gu*BCDYVoyBnz(-JPK51h81*Lakx&R^uRCquhb)z@-=**PXEZ?nlv_-U;di`W(x'
    '8|<Fw;~FqE=AkoL(UiY9pD+DnI}g*2=V``&>tS(IcuxqC&-rxc4RrbAhrNoAq-FtNvqS7ipu^MdzDZ&0Z%;JkSaB6xJrN!={bmai'
    'b?G|Gc*}$PSutFl?%nNG>fULAtSl`R5>KjeN0M6ua#&Y31{4DZsj50&UH9}_?J0``bc%OiZ!bpR-=}eTMq-p=(d8&16wXEE=hw!x'
    'Ai482byztxuD}NzC3OVk*WI%LI@2@`(a0ZtDB3I_!^t%G+mt*CY|#I~3a1PjeBU`wpfa+%@~Fzi4i(3_s1V@a`nifPdmtwBXRm<1'
    'o~j}6@uF1X9|e54*(39>7tV*GnTU;rL)LxhW$MYS$$48;h8il~Y7>oi)msUAp>(@>zmU-=PTCx8V{jtwZsfv9r=GoVd*{pX<+K(Y'
    '>X?kXW3?A$BmQkVu+7}%b$|~-$R8Ie$S`4U2nOG!g3Gc4R&aW}Pi2o!5JA&WRWjjv-;QW;e;4hP5KSw$expfTe(>QGd(<w@>~4Y+'
    'sSC!(Fk%!i<z1@pFED^<0W(}AoBP>W>0GwG>idLee2n|cmbD);Vk`@0R|W5=D;fkLI$2^tqU`!GRSci^)E2jG!94UPtfrETFgeMc'
    'OTy&s-V@Jf_bK@+{tDc%*p8P)t18ElowG>|Nq={@83PKgta&WI8i|v8w}JYM4@<Yhhaw8Bi_!*|F8<JY3qe7=I>v%*?w~n{_R8V('
    'd6U{DID_<)!`XvAH{5c|4v{~1cIT-pZL{T0zF#IR!%lqpY9%$TA5N;*q{p)m+JKb$Ax>WB(7d;xnR+>ln$E?y=ee-g)Hf+z3B0DE'
    'f@daYaotmry|(PeQHMY7ufJ1&7RWXu3@q5kXv=+R-wGs{e$}eqQo6R`JWb*rfogR(Wf~H|-^v{<sXaCNY`{AC5)Tnh$uT>8z~_Ck'
    'hTepq#zx|DPa7^1T}guTF~U-xKO|%`idZ4iNQ67hALLjUr`7=IAc?fb;Bx395QV`#8<&By<9QZe5)1tKz1463!Ek&CIQO4m48075'
    'jOljbSd=UM8h4-q*tuW;3iSwOaolU#@5F&~4lpsxx0AR|;y=jWr-Pz-yQ+M0Y8U}@v_BhjME9()MCD5P<X76Zcc}&hsb1A|{SZ;n'
    '2U5Qr?%LqGi|eUxmU|<HZ3$?<^7zl|x9n)yI9A3A=ITLa6_Q#$t^|0TPEgjQk=ZpMefTBGg5>@}g<0;EB%wcjfI*0a`iq0&zL)hI'
    'Yo={&pHS7NC3dINn<ayUU3u88<rg=mJWO4)3(bCNp(Ic#B;ZHKp83YApk{GNaiZxY=caBWnYbkM$q8%8^?i~Ve%qV@5`Eu!06IO*'
    '1cVKh!T$xTJe8#4fmUr{2f7_1*9f)>^su(e2dVr#i#}@jguvqPes0CTjEav_r@n5!EH)xqrX_o=KbLI>B~@#cq~31hVNsZj84#<E'
    'Yv&4kQRWuDbaaQoI62gPIpkQ<x;K^dm<r7F#ws?18jlT37rlAlHk`S8l5c83`U8W8-5Tf+B(E!#q-VwWnlNtOyP$r4#mz3)+rEr1'
    'jdlS2>lk=~cg69-FIEv%6Gj(~Bc9Gm1+NDAN@hxv9Pq^+0e9gTBSz$II~)3vQC({Ch4h@mhKSxY_=m4z$D0Lk*bw(U3J&mP47Hrb'
    'K;_(6NL~NqSa@Qx%SAmX{Y%_R&jcS6lW>&8ksujWu@Oi)_Yt#gM*>;$0j+k75YARC!t}w?tYi~$#aXB6@bs~3tB98eL)Wx{rhviv'
    've?I$XU2F$>U8F#v|`U)Qnk!NA?r?#0Hp5)e^mQzSw5h-;Y{!WT6623-e<BRdsxOP8JPB-IgCnQnhl#W#QZ!|i56W%EFa>XUr8+;'
    '$?tmf&`Mv&84uN!^G2#3rU+aYj$|m;%Y_Zgi(mQG8|<&0jzoB_aJ3`N#X}FK^f%MOvIC3BDFy0IP27|FRFPx5y!D!wi&~=ryH90V'
    'q{D3Kk~mt+hwf%Hn5c4U_eC0UNI=+snr+x39BiypF<^8b2mhkz>$4j~HkL}bYZ7fAM3V!ZLJ}TyQd{6Daz!T<)XPXwDXkn5j!0RX'
    'IhJ(Kse6xh2284uC3U-7hI=7a(99lxlJS+b=f5V@7RFoSDR||m-BNBck5HtlHblgrv#iniiH~bA1ze`)4`7DAJ!0AYyF;suin50U'
    'tSr8=WBr5?>l2es9}C-lw>re{h67v|=Tn2)E?9Vb=D?JE*3)b!-|>$4fl{1H7%Bd{18lnK7pQ5YhY+`m{lO|hYLH)==^gvQOqu`P'
    'lN`^>(Le|z@CwqW*JL~s4lv~Tfi?R~KdCz4HI6Nug~gGqhQSg+2+5$B!*;ugLf?Bj&;jdA<mj?F+0hqNz!`0HN&N6iTPv1;sRZYg'
    'j&D{*JA%n*Go_)F%5#RvXAO@WUl41%4i6P30^j>Q)&uUL>c{l+3LeE`Hw5iL`pJKNthS7*a^^Mm&w}q(+1k==^FjGB{Xi(4S7yEe'
    'ZzUgE<oaf@^JSg23`E_g+rr_COK~inwh334vqqhr*3{*&;_v(VxoBl{Jp-%#MxJ6d<Qfa=rZueBYDojwn~qM}rsv*eCB0u2T4GD3'
    '=!1frYGY4Sk!ZFh&KFw#OcBrMw>equbd`+F#CiK;1}de4tVaQPx%^nef&Yox4L)8UOut>Df~9R_R#oVL+rybozsJVp-3=3J>HenD'
    '6a~6AYl6*VUAk)^TB5}?qI(LSeE25oCtXvv^D_mu<33uJuHN#Ni9<MG&lW&pk&C5qUd<DSS?ep?9~z4_JPf#t-%i-V2}Xg4*1{*+'
    'X3Kr-|HdEwxvK|x^-tyte00SYF;x+#{6mFz8Y=B`z<RP29sjhlZe{<_yQ-RKHL6|)bluqQgCO`TYjp#vZ*lCo)6mC!Dlk(6T3nAO'
    '>p_NV)O+0n3dkIV4+T3XaA35f0VkgCENvD1T4c=TqT>V3o(q?%w~N}`AlYbD|A8#MQb~#Q%@;W?!Nbw`)LpKnHS349@^E^F0eeAL'
    '7S3oDnyRlX^K2wQy=swcipr;xx@I5*$D2NEm;R^N+e<D;9L@i?hi-9GZiMt9t}S2V{k_yD2M(kzGM`cSk3t@d9=JV=wX2lxP@GNN'
    'y(=VQB+<$~-*ne=_1akjx?gLG!DQwLrZLg|johklEvmts2E7nK@UgE``Gy~looN;Tf){_*i&^ZAqo)bOsRR(HWx=40&~5M6@ThE~'
    '2{1z>WE924?`Ss+1!D@gG%B}OgLf^RcjIc?LUAqookBXGUgkL&&sJx1ZETGtPYS1Xr)xXBGMe8X_Si#_25?bKzujcQ70-?5JoMoS'
    'd#G9_KT2<iDQVjOPVlvGYH8*ysZQj_TUg;S9y7A>TtE{<udhiR2D+h%eU-90%NS20l$`{{SY@c8i}y1XSenCf03|l8sJeM;KAgnY'
    'cV=GjDDp(#8z{vngerO&5OnXV5&*EDntY784Ch~B>XGQRE9h)l_H+k%VvGqw5C;6EXZaT_&doVV+Zz3^5e)8y8q8wf$ezH<l(CU@'
    'dr7}B_Gi5z3HwUOr_yB`eY!qBBl`hyTu?EwJUf?rixwAWf1VljxQ`g0{-kYYzJPGvSyHIazs$(Rf`eI(gtd04*48{Q-+bWkh*?aj'
    'luhv52w?v4F9!1}-;$E0nIE$#B5fO-^;Lll=iIRnLPy_QUy!h&FkSiB>Chf$bAG}wIIRd}hIvwBrbIGWZW{YfMf0jkaM3Sa#HtD9'
    '_B@hyXw+EFXM_3cM=a!8r&HXeT=XO3h&Zh<ANPI)%-ri;wRc=%!1g&aDc-b{b??p3Vu}ytIp&)yeqi`jccWZ7T@_`a(A+9WS^!@J'
    'B>p)_CJ+4$kWRQY+{NL-UcbjvbRQ>o_DH~(V7iM44ETfrMF`ASBt+`fYs~=XYz6FiTv88TRNp*_>(q1p8k>}}U~O66m+oV_+qE!y'
    '(rp~BFdE>tA>)`D$`{6G82pmkD97j7nBB4_p1-IW-26K)&HTC;chJAMWqEX|36g#oB~}xjiK<FjCOMk3g+6f5P@Zd4TH#UlOgEDH'
    'btzd*@+yq4^*hXBD<RM0%k}tR@w-<&Z1*Pmib!=E*ucQYL#iioRo@{LYIGr_fJ1t<MYvc(D~Has@!_PMszeX3gkxWhD+GX)t4Vl%'
    ')z*oa9%yfmG2BB8n~$%oJpIs*_R>eBo}}xyFqi75Boy>zukR)-tv8poPM)8XEFrt{wLI<z&WF1}bu(G&2ZNGoPjGlTaNb(>HbIs;'
    'ySVnF$A6j-^La=Sh^S*m^Nc<BuUO!eorj#W!kuN<fUVO@pHWXcK-2wgSZB|k6k>pfUS)J|F6Mp5+(RGewc!R;cc1tNob9aaM?UT9'
    'hEjf!!;&C5n=GPAb?3cO#-AD9(Hh>JBICCGMf`0q8^HEiFwb^*9ZFax@NrZ#RI_ej<B0@k02F`-LP|oie215j+B0^ak22yosmB?m'
    'z!%{JODsx)RT+K5Udc5-R0~+~;DG%jHiFr#yh=}+>s#Efay3q+`mF{czu>0X^nUG++c`?S4_1mBs2EcbwCNfS+K@38D*&3Avd}Jl'
    '9XI{@DAn|KwU|wG3{;gE(uk4fS+l*7#Q&3{3@K7-%ar1tpG*MorxI;^$v?l1gfbs_9d(_t`7IzfLpcVLIGPmbCRpTh3;zCRpbY@1'
    '5Yrrj!ExAUz-QRQ;+8AH)<>o^IN0?G*7KnEN5CRf3*VtA984};{l0WX$T<0+ryUeE4Mt?<7OWP#-5jLEZUq{$Ju|wPjl4*?59p^s'
    'AM~%cB~sB4^g;gQ?Rts77^8bXlo3eCphce@*WW@j46yXV==JmTjhG!YSeLhsuybBBZA0yN5y%9<-o3)dLKpQc$ox6&3BA?FAYjxD'
    'Ez>bjG*vkl0jDTiVnVS?BgHOv-5Nx?dKr;-7JL`GhmEcn#@ygHyB_rDwwM}pd2ZjOHH%|5(u*B13jDiF>jcXOD~ap>t=(e^0DPgI'
    'uiUFtqtHbL!u6Dm{#JUUB^;6$kB3yFjGv*gL!Pdw`9FKLKfv>e!R2X9$v`Q@v_tO)!9i1228w9Bb<|Z;`e8OAoxq=}M7g1JJDMGm'
    '>C+6!G!nx~(${D{ZP+<5BKQG+>#rMc@-PhWjDL=Q<|hg@q<@W_xL@)kT3hMLsG*nDac)uz8*yrJo`?FJ<@h$1l6H{SM`!6T!;!Yl'
    'Xv45XiRLq~<_<~(Q2uL52R3k2Z5-|mI<Jvi8~1J_gQtp8aXi&@6NlA#2NG323V!mNkAA+<!%-jR!NwIDcO&(r>Fc>8X?vG~?x}hz'
    'H*VhT&3Rlq4dX>^nNo&UTWDLURUMi%%aA4D-J2!-CvWTn<5oFdl1kPbxenl4#q3Sy62)=fM#p&YXm^~e2hlJhT;p4_N_&m1%xG__'
    'jXK2}PILi_TEMrkVqNXkkYX51#Tv{<lM%4Y;~rRH97Kz08NdIQAk%6;3s-)APz%qrM@WvK8%=MdQE+zPLZb&=WR&U%XTfVL0~7YM'
    's`#U~S4w?}twr{0nez)BdYojY2}R1=-qw!7AX#49rv;;~>qE$(-k4yK$*&-L>RO{DkM0n;uYVy$TWH<V-~1&#<S&4)+LRRn9s5Oj'
    '^$}l^=l<SscF6g#nM1-8blgyfj_@fb-usLaQqYXl)ypO^jQdtDey2PD?}mXo6hIkC%2%`f?j5~VEg161pG537P_gKPGh5N)jFLnd'
    'En;yc`Y8URDjjJ;vGPDq3eq{ZvhZPC02LTgjKaGJk?I?D2%Z9-9$zrl4qtidw}OwwWE@-Z1L%WGNl5wO=JgV$v_CW&KyrDk;agPP'
    'pPo|^JvCxJ_20MgCl<bR6Bv|id`JwWNgDigPUG>oMxu(~EsgjfUk>xwTBkv|!d2zF;uQwuFkI4Ng&GdZ6&yxJh$yKsa24I_dbkZs'
    '#SxF!oyu3h4s|HQFHJ;Zj%TB`?RciVk66M(v~Q~@n_a4!UE3@%=+Igw*VspVJA$ukD6JQbl;YKL6(|=hEE>M8uE0kyt*7!GR#;lm'
    'v^6O{B1qWe_WO)E2L3T619mEkUf<=f5x0ZWCqa9LC_@uK0b-`WvN$Kyv~++3bf5Z@%;$;}fT>&?;Dn&YQs-l`ea|h--wst@W&jyK'
    'FdMqJu)}vU64d|jRpKz7T4OePHb*%|@lb11&%<aT?@Lj<eKUT<UEQl)2LN>PqKC7YDM`j*HvUkJg0!9%cKY~T<+9DumWwG>Ba>AN'
    'g*^awiI_&y!1{yKvzL;k;=C~Ia!M?w^NT=2xCh{=fBE`HoR<nbtwt6j%lvK?w){LeaGXCgRc-x*gUVR1$$xk$@~$hvD8al?7$^%)'
    'oPR&o$N0~Z1%eSTfhbvy{^&Ilg07SBH>`S6q^e(<N-V)e0`}9ru`TmU-=;<WIO&ANLK}yX$VY1Uhsd0t!@r-r%->l8J}tX`13RO;'
    'Np#9Lyc2vS&3D!_49O5A5S#A?yr!3^>o0v19?N>xV3zOx<!hzxCLi(v>#43>M;4T(LxR@(t3TL3Q=CW-C=E?2^qS_r8WX@HRIsq5'
    'kSo(+(=_w;Jqdx%r+zHdd|abR_!<n4Es{^cQ3kfe69Wi(t~y8qPP^1d=XH`ti(u6^Rpa57*8ij#jPz;6zDwK|8#(p;&I7=?m{Y8X'
    'h$f~VQwQQ4u_;i66K~9_e6J7XD*YXj)o)cG0I?6h<_;Px+V~Vt7{5mZOs7=-5ACAF$tGB~-j@57=Fs|g*+ee@Xl-gUa&|Z=Y(a2J'
    'kk|tmYMd~nVc4bjYrx&@wd&V^aP?I7&^|E}5ONQ$y&A`f-`V1^g`hGU5hWoSUG^D#Ye*ke`+xvJ^-lVYh4e4f&A$6+77AkUqgu=*'
    'jpZ62EJETRDMMNsOw#8C+D8z3&9I}L9~_n!zi?1h%QT!>kJvJ^U~I3#+#K{~%*&CcF6{8VjaNRs)37mpZDTs|G!@Yc)geTPSAfaz'
    'ChL3YKj9Y>Fz-}6ykw^6p^h6Wl~D%PdDiY0XIvg6km+-a@x#eXrE~V&N;<s;20ge%jY{WD0qFUn>#_Q`jdB?>4<4S>x&4PCF6twq'
    'P875{hR3cdunCk7`2I)dp6xWoDIE(4w$3`@!18K^hy$t|GQ?{-Ij$KV9GAghpzh9jQX>!3?rK=o{byJ;&8iy7)mmKB`y#R?Br~mu'
    'WJT=|uTc)r5uQT`2qMCuWa^?P13aSpH4EW^Q>pbCM&NU`f4|!W#t<8Q1FR-RZ%*;28#GH#&*hh<vl*l_?-+~N!GDeR#iRFW`5&>1'
    'nnB`IEiFy4>vL<r{N!X<L(}H5kB1Hux?6tKT0Ml-TT(iDar@P4eAHz#j_6Z8EjFeSsKV0+r6;$XflA20aiN6DbiIKUi1CN>Stmu7'
    '*lK^P2U<*as9&;7;@89*auK}_xQq<p2P-bOi}OT-pzxN7yUK0UwC65#Jc3pKCl`2Z(v!NMQR!b1*BhTJ8?$MR=pFlck-!g*EZ(=w'
    'Pq%W9wRwZX)#cfd09k4R`hXn?6&EM~0*9bfAJMClDg2&%+5LBxzDXJzG1+rmItJ#>{H=+!6Pnj?X0%7&p<zO+KN41Q2qC5uR3svP'
    'U2H0OSorn%t^GIu;Y^a1uG=IGBqztCNd8~BQHr2<J<EpzW&99m2iz|g*zMYwa=~F&tEI&RmOPvYpv&J_*!82klQS|*o|2B(An>td'
    '^|&xM1U9y6Xj5V(aex`xKgeb0Tl1v|SN>k_&jj7-Bj7#i;BMpN)b!Mwfb5UJT&D#g{a<b3;B3sbiprCg;Y(DCYv=UiDBkq9bHzn}'
    '>K;<PaH@U+*_;JiBz4F&C6)aJAc@}kt_OiQBX&##*a?R*_xzXti<!TANpu-dX0+#u&?B%DSCqNj>mM}pU1D19bxVf+1xJ346MiBs'
    '<jGpgY<0~!`t}o<$MnmJG#-&)v4e}w9{vsR&xq3TPJQ+ObvmrLqFMey;W{SqaXyAyEslSjCm>u$aterPiIy}5_l|0MPwHECzNWm`'
    '1i=U4-(pt#eqbA&6bXBOr?K@TRPz(+3yW&q@;y<zZ)$uXrmDJ29AoUe@>*3r>fK6i>bYtCsZ|~>-X7Cx%)>e@lMk>C`vDs$%QP)|'
    'mV=|Rs|#=$odHSa0mo;9t0;q>IJ$;e%D=S7YeCw<+T7P)-{<}44g@0wswW95n04P|c>B%w8k62{WK3I>{+cVg&qVcsM5f637{8kF'
    'r*({8dcy;EZI~`x$JCM7sd8(_>P|SAyzDaDH0u?fPo_pa-zhA0T=~SFaJvEHgE3KPAt;4IhMwkRi(1-L<e{%<MR)Zh!49bzsV0V^'
    '?Y>~w6sl{|E&81`n@rmM=WHa5J1OO36@HFRNL8ztKG@c;yqmWdFt}e&A>K7tavLKIVjg>esQXv|*A%uR#eb%$w05`;5mnWZK8;sF'
    'J_cZ8kwci+3x45$$pspne$$1wUv92$V<Aw+TaqkdCtXVLWu4`D75NIPetOvnHjDS=Ta!CevPFJ&t1jZ<-37q>K2=E5%Kt9X-O9>k'
    'rIgmb`{fsL={=rk>eZVG)iXbk#p~nR#CE8J9z)<(Uvh@xxG?>`2w+E-zh;fq^=;2#b{K}Z=r~R*e$LtS_r%qyv9(APJl3O#i=Yt{'
    'z3lwvwns2jS~9JmAKFa0$?;85#VvPHG|==G5hJOE8$+$8WA+$!o3eF|M_6gU%Cxg+?NeLOvHmMGx1=(wTm-Bt0ZExS)Q{`wE$($W'
    'Qpug4VQQS=Z!yAPRpSswaX#B>>cIm9;Dk)5P9{Z#^JFC$>tu!<5-LSzi8m5U!{av?G|+jf%^}xLtbm<{tLGotdrCw1tq4jM>B6QB'
    '@;Pl1|762Y%(bH*pCjO$BxvY$*G4zrxT;kIUs#g|8Su>{z>T)VX6!<QEhq;>(w>jx6)Z4if_ikm<Z%3V0U6#EPG!p<Ma&geQvv*s'
    'Rs%Dw`x780RN<d|s+1av`Tly5G2HxY8fca`<gNJld)g&>$V8~`N-gRd{f;mDO%szh&g8^;w?}H|{n6)YA-SClo2A7?6q#*NAvU+U'
    '6PmqPH#BAH3=P%Were~dqAz_z0H{V$^vDLOI%lF0gI@bQGV3{B4{F*b>r#ObP2S|gMo|p@41jDYp^Q}DYE8sEDnB+^`qjFxhPbwm'
    '4Z&^$2ZB%G8ac+MtfBB_4}u~;_{ONB6%cl7;-66SVE)zqRRP_CoFRkUM6Rq#eE)U^b>5<QGk@y{rSxYh=(D`(I~FrKNYB;(AP}*>'
    'lC;7X$4&F2F2Z%$d?B9zzeJ;A^{uxs=Z*wgAwHk;4zW)%t>C3J-k$bAWbC_P2N#KJky_;MN1}IHvb+(Gpn4TQkNPT?2xJel@_<_|'
    '+}zLK-gi`K9t0{zQIHlF-8mE)AWh&M#@xh8sVuHZH5LmTf`H=r=C?xGe}4~n+Qv~(ML^84q?(^nV7{wqCWz-Zh)m`-bzlF+|B|<f'
)
PREFIX_RAW_SHA256 = 'a5c9e3a173f7cba130fe01d16afe0c59ea4fdc004893320245a8a911d251796a'


PREFIX_DEPTH = 23
DEPTH = 50

def transfer():
    return Transfer(165, 164, 132, 209, 138, 209, 457, 574, -13600)

def open_module(name):
    from pathlib import Path
    from runpy import run_path

    directory = Path(__file__).resolve().parent
    if (directory / 'files' / 'moments.py').is_file():
        directory = directory / 'files'
    return run_path(str(directory / name))

def decode_prefix(conv):
    from array import array
    from base64 import b85decode
    from hashlib import sha256
    from lzma import decompress

    energy = open_module('energy.py')
    payload = PACKED_PREFIX_FIRST + energy['PACKED_PREFIX_SECOND']
    transformed = array('H')
    transformed.frombytes(decompress(b85decode(payload)))
    need(sys.byteorder == 'little' and len(transformed) == 2 * conv.msize,
         'packed prefix shape')
    width = 2 * conv.H + 1

    def invert(offset):
        table = []
        for k in range(conv.msize):
            i, j = divmod(k, width)
            left = table[k - 1] if j else 0
            above = table[k - width] if i else 0
            diagonal = table[k - width - 1] if i and j else 0
            value = (transformed[offset + k] + left + above - diagonal) & 65535
            table.append(value)
        return table

    negative, positive = invert(0), invert(conv.msize)
    need(all(0 <= m <= Q - p and p <= Q for m, p in zip(negative, positive)),
         'packed prefix tails')
    digest = sha256(array('H', negative + positive).tobytes()).hexdigest()
    need(digest == PREFIX_RAW_SHA256, 'packed prefix digest')
    return negative, positive

def energy_evidence():
    conv = transfer()
    decoded = decode_prefix(conv)
    regenerated = conv.evolve(PREFIX_DEPTH)
    if decoded != regenerated:
        raise ValueError('packed prefix must equal every recomputed tail')
    moments = open_module('moments.py')
    gain = open_module('gain.py')
    energy = open_module('energy.py')
    if energy['THRESHOLD'] != gain['T'] or gain['T'] != F(-13600, 40000):
        raise ValueError('threshold agreement')
    rows, initial_energy = moments['compute']()
    cost = energy['energy'](rows, initial_energy, moments, gain)
    initial = energy['initial'](rows, moments, gain, F(cost['first_pass_upper_1e12'], 10**12))
    if F(cost['energy_upper_1e12'], 10**12) >= F(30937568299, 250000000000):
        raise ValueError('strict energy ceiling')
    return {'directed_energy': cost, 'directed_initial_majority': initial,
            'prefix_sha256': PREFIX_RAW_SHA256, 'prefix_depth': PREFIX_DEPTH,
            'recomputed_pairs': conv.msize}

def terminal_evidence():
    conv = transfer()
    m, p = decode_prefix(conv)
    m, p = conv.advance(m, p, DEPTH - PREFIX_DEPTH)
    truth_controls()
    literal_minimum_controls()
    upper = root(conv, m, p, 400, 400)
    if 2 * upper >= Q:
        raise ValueError('terminal majority')
    entropy = semantic()
    if F(entropy['coefficient']) != F(15683754216081791, 3093756829900000) or entropy['checks'] != 72:
        raise ValueError('joint score or dual accounting')
    return {
        'threshold': '-17/50', 'graph_cap': '30937568299/250000000000', 'coefficient': entropy['coefficient'],
        'terminal_upper_num': upper, 'terminal_upper_den': Q, 'depth': DEPTH,
        'prefix_sha256': PREFIX_RAW_SHA256, 'prefix_depth': PREFIX_DEPTH,
        'semantic': entropy, 'root_truth_controls': 72,
        'literal_root_minima': 162, 'literal_NAE_controls': 9,
    }

def run(request):
    if request == {'phase': 'energy'}:
        return energy_evidence()
    if request == {'phase': 'terminal'}:
        return terminal_evidence()
    raise ValueError('expected the declared terminal or energy phase')

if __name__ == '__main__':
    import json
    print(json.dumps(run(json.loads(sys.stdin.read())), sort_keys=True))
