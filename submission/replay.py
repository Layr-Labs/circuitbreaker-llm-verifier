"""Directed recurrence and exact regeneration of the packed prefix."""
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
    tables = {}
    __slots__=('q','dg','log','p','cum','dens','err','tail')
    def __init__(self,q,logs=26):
        if q not in (1,2) or not 12<=logs<=29:raise ValueError('CDF parameters')
        self.log=logs;self.p=1<<logs;self.q=q;self.dg=S//self.p;self.tail=1
        cached=self.tables.get(q)
        if cached is not None:
            self.cum,self.dens,self.err=cached
            return
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
        self.cum,self.dens=tuple(c),tuple(density)
        self.err=up(8*S,12*H*H)+up(S,H*H)
        self.tables[q]=(self.cum,self.dens,self.err)
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
from itertools import chain
import sys
# weights and probabilities are distinct integers.
LS=13;LQ=14;WSCALE=1<<LS;Q=1<<LQ
def need(t,name):
 if not t:raise ArithmeticError(name)
def ceil(n,d):return -((-n)//d)
def packed_arithmetic(size):
 from ctypes import CDLL, addressof, c_char_p, c_long, c_void_p, create_string_buffer, memset
 layout=Struct('<'+'Q'*size)
 def pack(values):
  return int.from_bytes(layout.pack(*values),'little')
 width=16*size
 def python_product(a,b):
  return (a*b).to_bytes(width,'little')
 def python_merge(a,b,c,d):
  return (a*b+c*d).to_bytes(width,'little')
 try:
  gmp=CDLL('libgmp.so.10')
  mul=getattr(gmp,'__gmpn_mul')
  mul.argtypes=(c_void_p,c_char_p,c_long,c_char_p,c_long)
  mul.restype=c_long
  add=getattr(gmp,'__gmpn_add_n')
  add.argtypes=(c_void_p,c_void_p,c_void_p,c_long)
  add.restype=c_long
 except (OSError,AttributeError):
  return pack,python_product,python_merge
 first=create_string_buffer(width)
 second=create_string_buffer(width)
 def compact(value):
  low=((value&-value).bit_length()-1)//64
  length=(value.bit_length()+63)//64-low
  return low,length,(value>>(64*low)).to_bytes(8*length,'little')
 def one(output,left,right):
  memset(output,0,width)
  if not left or not right:return
  lo,n,a=compact(left);hi,m,b=compact(right)
  if n<m:n,m,a,b=m,n,b,a
  mul(addressof(output)+8*(lo+hi),a,n,b,m)
 def product(a,b):
  one(first,a,b)
  return first.raw
 def merge(a,b,c,d):
  one(first,a,b)
  one(second,c,d)
  if add(first,first,second,2*size):
   raise ArithmeticError('packed product addition overflow')
  return first.raw
 a=pack([0,1]+[2]*(size-2))
 b=pack([0,0,3]+[0]*(size-3))
 c=pack([1]+[0]*(size-1))
 if (product(a,b)!=python_product(a,b) or
     merge(a,b,c,c)!=python_merge(a,b,c,c)):
  raise ArithmeticError('pinned GMP limb layout')
 return pack,product,merge
def runs(indices):
    # Constant-gap runs let the runtime gather output digits using strided views.
    i=0
    blocks=[]
    while i<len(indices):
        if i+1==len(indices):
            blocks.append(slice(indices[i],indices[i]+1,1))
            break
        stride=indices[i+1]-indices[i]
        j=i+2
        while j<len(indices) and indices[j]-indices[j-1]==stride:
            j+=1
        blocks.append(slice(indices[i],indices[j-1]+1,stride))
        i=j
    return tuple(blocks)

def coefficients(raw,blocks):
    view=memoryview(raw).cast('Q')
    return chain.from_iterable(view[s] for s in blocks)

def symmetric_axis(inner,outer,stride):
    pos=list(range(inner+1))+list(range(inner+stride,outer+1,stride))
    if pos[-1]!=outer:
        pos.append(outer)
    return tuple([-i for i in pos[:0:-1]]+pos)
from fractions import Fraction as F
from functools import lru_cache
# Eleven-degree exponential Taylor polynomial in homogeneous integer form.
EXP11 = (1, 11, 110, 990, 7920, 55440, 332640, 1663200, 6652800, 19958400, 39916800, 39916800)

@lru_cache(None)
def elow(x):
 y=x/256
 need(abs(y)<=F(1,2),'Taylor lower domain')
 numerator,denominator=y.numerator,y.denominator
 polynomial=EXP11[0]
 power=1
 for coefficient in EXP11[1:]:
  power*=denominator
  polynomial=polynomial*numerator+coefficient*power
 need(polynomial>0,'positive Taylor lower')
 u=(polynomial<<80)//(EXP11[-1]*power)
 for _ in range(8):
  u=(u*u)>>80
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
 norm=[int(normal*v*fpi*elow((i*substeps*h)**2/(a*a if mode=='A' else 2*a*a)-v*v/12)) for i in K]
 return f,g,norm

def weights(d,mode,N,K,substeps):return jensenrows(d,mode,N,K,WSCALE,substeps)
class Transfer:
 def __init__(self,d,yd,main_design,aux_design,contrast_design,Nm,Na,T,substeps=2):
  self.core=65;self.middle=102
  self.d=d;self.yd=yd;self.T=T;self.Nm=Nm;self.Na=Na
  self.main=M=symmetric_axis(*main_design)
  self.aux=A=symmetric_axis(*aux_design)
  self.y_axis=Y=symmetric_axis(*contrast_design)
  cd=d*substeps
  self.msize=len(M)*len(Y)
  self.asize=len(A)*len(Y)
  def sample(axis,v):
   i=bisect_left(axis,v)
   return i,axis[min(i,len(axis)-1)]
  pin_cdf=Gaussian(1,26)
  ell=F(2*28*28*41,57**3)
  delta=F(13,57)
  pinned_threshold=F(T,40000)
  v0,v1=2045*57,1148*56
  straight_den=41*57*d*cd
  extra_den=57*d*cd*v1
  y_lo,y_hi=Y[0],Y[-1]
  y_lookup=tuple(bisect_left(Y,j) for j in range(y_lo,y_hi+1))
  def ix(ip,jval,axis):
   if ip==len(axis) or jval>y_hi:return len(axis)*len(Y)
   jp=y_lookup[jval-y_lo] if jval>=y_lo else 0
   return jp*len(axis)+ip
  def setup(isaux,N,nodes):
   f,g,norm=weights(cd,'A' if isaux else 'M',N,nodes,substeps)
   tc=2*N+1
   pack,prod,merge=packed_arithmetic(tc)
   total=prod(pack(f),pack(g))
   scale=substeps*(2 if isaux else 1)
   digits=[2*N+1+scale*v for v in nodes]
   need(0<=digits[0]<=digits[-1]<2*tc-1,'output digit ranges')
   parts=runs(digits)
   need(sys.byteorder=="little" and (2*N+1)*WSCALE*WSCALE<2**64
        and Q*Q*max(memoryview(total).cast('Q'))<2**64,"no carry")
   B=list(coefficients(total,parts))
   need(all(0<=s*r<=WSCALE*WSCALE*(1<<40) for s,r in zip(B,norm)),'Gaussian lower mass')
   plan=[]
   full=pack(c*Q for c in f)
   k_axis=range(-N,N+1)
   if isaux:
    locations=[]
    for k in k_axis:
     ip,i=sample(M,ceil(d*k,cd))
     dy=ceil(yd*(41*d*k-13*cd*i),41*d*cd)
     locations.append((ip,dy))
   else:
    ray_locations=[]
    aux_locations=[]
    for k in k_axis:
     ip,i=sample(M,ceil(57*d*k,41*cd))
     ray_locations.append((ip,57*yd*(57*d*k-13*cd*i)))
     zp,z=sample(A,ceil(57*d*k,56*cd))
     aux_locations.append((zp,57*yd*(v0*d*k-v1*cd*z)))
   for j in Y:
    if isaux:
     od=[ix(ip,j+shift,M) for ip,shift in locations]
     plan.append(itemgetter(*od))
    else:
     straight_y=41*41*d*cd*j
     aux_y=28*d*cd*v1*j
     od=[ix(ip,ceil(shift+straight_y,straight_den),M) for ip,shift in ray_locations]
     pd=[ix(ip,ceil(shift+aux_y,extra_den),A) for ip,shift in aux_locations]
     ff=itemgetter(*od)
     shift=(pinned_threshold-F(j,yd))/delta
     vcrit=shift*ell
     scaled=vcrit*cd
     active_at=max(0,min(tc,scaled.numerator//scaled.denominator+1+N))
     g_active=g[active_at:]
     active_columns=pd[active_at:]
     active_layout=Struct('<'+'Q'*len(active_columns))
     active_pad=bytes(8*active_at)
     if len(active_columns)>1:
      active_read=itemgetter(*active_columns)
      def active_pack(table, read=active_read, weights=g_active, layout=active_layout, pad=active_pad):
       return int.from_bytes(pad+layout.pack(*map(mul,weights,read(table))),"little")
     elif active_columns:
      def active_pack(table, col=active_columns[0], weight=g_active[0], layout=active_layout, pad=active_pad):
       return int.from_bytes(pad+layout.pack(weight*table[col]),"little")
     else:
      def active_pack(table, pad=active_pad):
       return 0
     w=[0]*active_at+g_active
     anti=[v*Q for v in g[:active_at]]+[0]*len(g_active)
     vbin=ceil(scaled.numerator,scaled.denominator)
     corrections=[]
     if -N<=vbin<=N and F(vbin-1,cd)<vcrit<F(vbin,cd):
      end=F(vbin-1,cd)/ell; low_den=41*d*shift.denominator; high_den=41*d*end.denominator
      low_shift=41*d*shift.numerator; high_shift=41*d*end.numerator
      for at,parent in enumerate(nodes):
       sb=substeps*parent+1-vbin
       if not -N<=sb<=N:continue
       lower=57*parent*shift.denominator-low_shift; upper=57*parent*end.denominator-high_shift
       mass=max(0,pin_cdf(upper,high_den)[0]-pin_cdf(lower,low_den)[1])
       if mass:corrections.append((at,od[sb+N],pd[vbin+N],mass))
     on=pack(v*Q for v in w);off=pack(anti)
     plan.append((ff,active_pack,on,off,corrections))
   total_weighted=[base*Q*Q*factor for base,factor in zip(B,norm)]
   return (isaux,plan,pack,prod,merge,full,norm,total_weighted,parts,len(nodes),f,g)
  self.A=setup(True,Na,A)
  self.M=setup(False,Nm,M)
 def apply(self,setup,m,p,n,u,plus=False,previous=None,round_number=0):
  isaux,rows,pack,prod,merge,full,norm,total_weighted,parts,cnt,fw,gw=setup;outm=[];outp=[];both=[]
  fneg=m+[0];fpos=p+[Q]
  if not isaux:
   gn=n+[0];gp=u+[Q]
  KK=2*LS+40+LQ
  for idx,row in enumerate(rows):
   if previous is not None:
    y=abs(self.y_axis[idx])
    interval=1 if 100*y<=self.yd*self.core else 2 if 100*y<=self.yd*self.middle else 4
    if (round_number+y)%interval:
     base=len(outm)
     outm.extend(previous[0][base:base+cnt])
     outp.extend(previous[1][base:base+cnt])
     continue
   ff=row if isaux else row[0]
   fm=pack(map(mul,fw,ff(fneg))); fp=pack(map(mul,fw,ff(fpos)))
   if isaux:
    ac=prod(fm,2*(full-fp)-fm)
    bc=prod(fp,2*(full-fm)-fp)
    if plus:
     dp=coefficients(prod(fp,fp),parts)
     both.extend([Q-((total-y*r)>>KK) for y,r,total in zip(dp,norm,total_weighted)])
   else:
    _,active_pack,on,off,corrections=row
    gm=active_pack(gn); gp0=active_pack(gp)
    Z=on-gm-gp0
    ac=merge(fm,Z,full-fp,gm+off)
    bc=merge(fp,Z,full-fm,gp0)
   ma=coefficients(ac,parts);po=coefficients(bc,parts)
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
  p=[h(57*(40000*(41*D*j+13*Y*i)-self.T*41*D*Y),13*40000*41*D*Y)[1] for j in self.y_axis for i in self.main]
  return self.advance(m,p,k-1)
 def advance(self,m,p,layers):
  aux=None
  for step in range(layers):
   n,u=self.apply(self.A,m,p,m,p,previous=aux,round_number=step)
   mm,pp=self.apply(self.M,m,p,n,u,previous=(m,p),round_number=step)
   aux=n,u
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
  ix=bisect_left(s.aux,ceil(D*s.d,den))
  if ix==len(s.aux):return Q,Q,Q
  i=s.aux[ix]
  iy=bisect_left(s.y_axis,ceil(s.yd*(P*s.d-i*den),den*s.d))
  if iy==len(s.y_axis):return Q,Q,Q
  z=iy*len(s.aux)+ix
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
    cap = F(15415990057,125000000000)
    coefficient = 1 + 1/(2*cap) + delta
    return {
        'checks':len(checked), 'semantic_increment':str(delta),
        'weighted_margin':str(weighted_margin), 'private_path_margin':str(cutoff-delta),
        'coefficient':str(coefficient),
    }

PACKED_PREFIX_FIRST = (
    '{Wp48S^xk9=GL@E0stWa761SMbT8$j=Znj<_FVvG54<1%B~gK<d5GRqY8FNb=+6FI^HHjtNujl4Z_x1KO1+F3joFb`KT>JU{1g1`'
    'u+5I1nH+T^f=vF9YnnGj!7@8Z)ic{kGkZR<FxzNv%s*hdH4&Yh%ti7;h)ycnfwAu>50sPQ&!OgHXb^~lqX@p!b=*tpc2kiRB{&p0'
    '2-t1tTkNe+Ei1^ET0-9$mYKSxBDdIk_NpR;jm>aY5uKVV;1fpEw8oXcxUx;-{~ILdrfq2aX-~8!zHn*L0WWjbXIS>NbaZw5JC}yG'
    'IRh{?Tp#N15azq5B$Q@KfLj>z%dw{}^)O<A>77&do3v)JNSupiX?+pbO}GABtYW~VaWbAa9)ozoj$o8$-=&cswOOY=BNk1o%OAgN'
    'sp`hxk021{UtDxJV=WSNwY1DyU8yTd<O_MR0VKZfov}1Sbz9@g3omgXmbHnavJO1|qZCczyqzhPh0Ee=<Y}mM->AE>mWTNxo;IU!'
    '12<Mo&TFmxadm3zMkG>oU9@{FFPXDeBf=iE`mSp+7}likcviw3!lP`BA`uA$Y<;ia^sPk$#YeSC;=&&qO*^tBsQy4~I_K+-5Llf+'
    'b)C4;V#Oa;VsEMeP2?q@gNuBiTq+&eFnMa=Tu1xNGksfhdizqA8~x(ldtbG<b;57@d-uT{Ka;AHSky6$6Rt>TIoAH=dkuQtTfhF~'
    'WI5vbgiODzfR^=p>aQ@iak|f&P7wUE{Ra3HzDi8{R5kVoLLO```C$wc406rCiqJi`xP)tddRsL^2p^|ZJ_=zWBmi7?KEVOec}L{n'
    'SeP)9D%iW1o61%DmL1MnmZT4>g)J)al19FK?xf=ZkNaw=5k&%d+Rm08d(Ut{(dL8`K1CeP^keF`2o9q-I<kB2Fn8==vE<|{Oy-tQ'
    'CzE0$7>7dcIb?PPn_n2DWl(ce!PkDl6+R;j7{#-cBR(7nKL&Yw!10e6G61!K>|C*dT7Ck2dWh{E?7I6tD+z6I!Tr{B#lIV4|JiS5'
    ')VQq-gfJ%@mn>wfZUrdLvQRA-#H?pIamjP|qt44nWLe;{9b?Ey^&DPS_|8yam<cSa#P$`zM*C`7JHRIYm65QQT^~)U?T>^~P?rJs'
    'UFPWs9k`Ni&hAL=9d8YF4KVRS<C;aZEZBj|EpH(-gqiPJFu0k&Y7z7kiz@|B1{S0w&-f_B1HS2Mg<9R@71;p@@Jy{%J5KkLqonno'
    'ZGJS#T}v#j#D;XzsPd-d_y(XIj-*4ol8O=InkE;HI1PtlM{dOsNZmc-ZDu0@G1wQw`zF+E>S89a8S-u<WsbdyxSd({@55YT%BiDc'
    'wpc<wbOyGJGRkKr$gz1nRd;j_cUf(!tJaXQ5w=LeyCh-ii;`T)@R8abjpW!?ybx@`|2{^!h|ltqR7x&7Ju-T$4Sx}U&Ye7VFONs3'
    '|2Pl<$X4+{=iNLlp2msanJCel)>JO~XCm0Wjk81K_Lo(u&y<ifjVbdc<?8sCku`$2qv&#KzjIe7UtUFfMLOc*Uh~HMQ|-T<iB-cl'
    'k?euSZDdF<_v24%-bCB1y8a?di?(t*bVSt2F;_KpEwRT+Re2uAiHlji|Gpf1OiRJ{&?L5wP3M{BAr;f9nikwB#FYei(|D@-m4B?A'
    'WEyomDd+0EqU-HH5vzj~*~)2L?x2dGSRn?o$pxREnk=&8TFl4QGXlTWPVOvCpnXYcb-*Sq4&T!fEbKAuznklkz_Ec|(@>D{^}0e|'
    'Z{Z$u(U}<VX~P=b)``-WT1^&RcJ^L}alV0U0Yk?ny#+Funw(xGFcVb<ZZjv%>UHl8Fw4GZl!8;4PjXFD8_EF$$WQTq?l!MP!262{'
    '(Sv^5^-VI!=iF(YpqJ6%^Kgy?V!HhCBPp!=K6|B^v|Y#^(%zUw1!EK!Mj|d`GpkJ=?CiW)BXG}oP25#b1dR6=Jq~m|-2rvN-UsbI'
    '-!{I{FUIVlt<ONWcSgsoq_5U>lzk*B?qb{>8v{*#h@~$*+mDKIOo708f=ce{VBZ^RFX3lPIrb93y`ks=?jv~l^v3n}Fk<(LR@52X'
    'V4>c~XGl*Y&Ry&8tk6fXFRO71jb|0&6C-2U;7;;%VNxAP!-W>Gd1nVpEY8acON267h>`5Ksv3mMZEuE+ZpN??eDdT@kIQf4Obw8Z'
    '#+4Iv)RuDDk60`3YcTt3IyKLZ1Ew@Yrdf3^_70^K*I&KFdubab#fAsS<VHCRWy3Ubk*k0h?^6yOQY%ib(WI`1+G$WZlYO3|Zaq+V'
    '(*#gtZ>h<r73K;nURCAI`?h_(FtFqCf#$*hltwq<>F7|qW^q1NK^!UDA~PbcsBA|{Ke*waksa|uM-b>!<z?(ztyKOMtW&|v!|D6d'
    'Tif;%<1pcIR#t5pTAgTOL2L^^3cosPrcr8_G>7SBNN8bQ<yU;i?UjVCO}1VDz36rWS&;S@IRO8xmP{!|1Zn|kG8cZvf;fF;eW*1y'
    'SfN#~4%$k)GMenEq+|7sK7gRC*tQ8h9ZmkhgE1TB>$DNr!v{XJ1Q(GN7p`XBR9}+x(mgZ(kUPqk(71OE1At78E6%X$Po;-v-`RAx'
    'Tii!`Bc##WoXUiH9e`HP*!<DHf0H4!!qB!V8Tc}b6;Va3XVS|N8x@HLW(G?AnUUvD<{SLvZY=fEuN%HCN{Dhz56;iVvD13QuNcx4'
    '8PfH-r&&VD!?ztBM-U(*`EFKAA2k^bENf)qH_~Vpc09}hOStvIHn=%-#H+J9cU|auj#$l!#}C|QW=iwb_?~sjZb}R`_-Ti<9!}!2'
    'kinT6l$PBor&8RoR`xZfk=TdFwxt~kKh?gkFP#fu+YtL*Y0gzYhuPz^LR~Pr+}~|)tD)#ogqRE)rbAj2F<;g#ydh*+yr}ykZ^5F_'
    'O4DjYcX{x10&{&5<W-8yV2{wg3s-4?sw%_)K(#0NPK_#}!}n$sx%t+n)&Q=6v)!f@H48&uS@EMSu5x+sF1$=aE@~lK0z4ApphtJ-'
    '{KZT-VXkmm)h977S-+cJW%a3b`>Dk60JQ(g1^L-%&-zw%?VtiB?RKtl`5r9W3=}D8Uxkyi6a$j>>FS*)sLtLOG1BfD)#cz$R`0r2'
    'XONhSyd6xoX`l+Guc;Fby)XNag}4~k9gBKRR8TTq%B1o36_NNq=MkVNFYW4=njv49O0qLo%?fB*$7P^Guo>L~N(K7942wKC0^9hx'
    '-WB*BT)^zV#zSVMbTt98C_nS;8b_ucXHEw5%11C}9o$_NFDwv0z<eEo3f1(2!-&e-m7lKx(q*t8Ur^{umn>noQsxe-D*wLX!nh`d'
    'pcMKW_cvN%1R0j|+sLt@@}(e&h_idv6P_@5z~%1?MiXdh5Iu&~(omSgGdFw@EWt)4k4!gPC8{dx|ID;9)&WVS#;`gs>(dcWse`>t'
    '-;i&YifESu5g6|%_Rbz3jlG~B)dE^NX|eCR0tR8Jz&{V-6a4CS%=`mpZB4yyf+=9<2ll;b7gxt3k{Tr>`s5b0r-T*9{^f2xzZFNE'
    '0d-Ay^<hcI(I*n2k60~FT8&wf!xLw<$42gf^`P@Aea&~#=9H<p3rI^y2=ZCWAr1EXQ6UN8ipW>Rgl4}_L-D$@gawQOF!n1i%?HgU'
    'Ge3`O8nB``>H4m8o(V_gh^>sV-RJA~pdjXI;w8O-e9lEej>gin<*aE3^r_O#HyyN{j&I8io%2z-)9oXMRwTfp*KXWTZ_iDynCj&Z'
    'dRoo%{IciV3avvkH3q*fv$*Y1=g4mjB#91qTZ+{2VAk6YL?!sx&+Cl4Onl^}OqbqK<eJ+k1ZMFBm-d;Edb%xv4Je<YWG!G;aW_D6'
    ')z*#zq<F8LvZ!fkoG8g>M#UGlYoSs~toA+}vwbtksZ9!`CHkF{{Yq?wBHFaWyS|tqOL8g`o9E$S__`X|8L2n@B3nE1>}4%_ib4c4'
    'q+D%nzljR~c;zvP!sTSD#u=7YhR$DrbIy>LpD<Z7T#k~WllBjGvA~85oo<<zn{k~U@8Uhoj$V^lOkkj~LJD4!O?rPbMx9Ll&=(o-'
    'et-Afw&=1HLwJ7^kFwaiH3|um0QF~X96d7&sqm_#_%0@39JO<69+etmqc^>L?Q-K&Md+ptVc76M4{uqlbitn|^suG+x~F^&zFbHN'
    '6WFh2xLsiR;NhM?;TzNl=j#UIso4izRI!AgV%lsiHKCwbY!Vhd<@61J$lq4q*$7dnkyIaf-TvP3?}|zdW^~uoBdKExybjXVQb~@<'
    '1n33?+|YWne~pz2kRQ>WD+N&);Q;*{i|gX!YKjmJ>GjEYvJpk)1noy$88|5Y{ydb&r$RXK$NkJW`*0^Ud^_-Q6^+ncdAs>(!4Nvr'
    'mF{_l0|b<GfTic5a5AD^W&GV<*m{q3Zhp77*fgL~X`hx{zR>>4yJL!y9P>$14cBe}>7+xa`ytcJ*OGNugQ3znqs2uvuXX0S&305t'
    '43kp8Iq%DHi7k1GtsmrEfbEsNK2w`2=%|o1_(hL<RWIfAM`?q;)GMdiYD0=Q5}!f}C<+KA&l6h4XVGNVTcJqp|69*2q2Trljssml'
    'o=asPw!Ph_A8Pn?Y1E}p9%q}-vOcb9i%E;A-r$9j-a{seQ6X>3ValGvc*?#(ebLrTx8zQK8C5Jcm)tOC3R?GKq##qaRPXhu&dlYp'
    'V{03AT~#E9ly(5sgO856qqykdG9ckQH~$4g=_uNyb)onXC@0|5uHya~=;0BNMw6Qho!Pf(7PA`KxGRwP)!%ek?L&lYB1ULfBqO2u'
    'J(5hs0CDq8b$QLI8#3<fvW?9+o33|W`4_H09t+N;%><S)@TV1+Al0wYMR)K}gAY84WovhN!da!8hzN~xO}I_-iV!)Dwptv;H-pf_'
    '*JCCJ-q-yF#ih%$adtvsI?z5V#J_=f2E-fl<+_|)r;O*9Hn+z=_0nA~Z;7a_gWT*ZP9(nEbkLY*sp|((n|GE?P%H{^;~%Y4FWGR9'
    '7~)A2IPHlxp79xa(_icI&4fEEVSSiuqhDQjQD3r!Fpe9T`4JI~s8<W<FV1bDMIVmM;)kQ*4kuBlU<MEKqfn7#eVVXu4a%O3O(~c3'
    'J$52SQogpQFy-hfv%_{FgXOLJ<nyd-N0^H~@{lm!y~Y*Bn)e>E_qC!0S{kl@C(DuDyiqG1bN-9j@uLaXxb|KEEr>M=B-%(rcy!@E'
    'Twh8yruw{u8hogkwQgtAK7T(|-c}F=gxz1%8`EA#3-CV)0Lu??4!6p7@&D1dt$S;M^$_hNHtN-FUBG$CFJ7X}0PUa0g;H?D=Crd2'
    '1FI|;GU|>)qAB9Cvfm;OH3c|vnG5(admI`mqf3`oyA_YlUvzA22-|+)^SkLR<wH4(p;#$ns8NAUP@?B*x+JX<|AC-FYox5onps;~'
    'Tt{THOY}Ks9?1KMSkmchxy!v%=^=r6l674jopV2Ww2NX(E4E3nRRjm*-Q@)ljqWYBhCsl}kSd*~x!nVRFC4a<bH3U0a_gvi;cJMx'
    'IN3=2l;t++xpUlyEmN6JDTH=b!Y*INqcmRwYh4?4y&JIl(<{T-lvr$P{%Z#%mJteE%k9@eytYb=$63)aXV5}Aj;sqi>ADeu9HMzq'
    '_4=L=)rU}b6m}_T=;-zh6qn)`mXkj6H)3Cj))<MwMSRhSmhCWOcknVKWePGSu={>Uid!i|Eb2yGnAB8y!uq>7^18R1ps}xsE%PG_'
    ')0Na%NO_P=pvXeuJ#ZV8qA8H9eTSnC5bybq7m1D6%n(sC^1u3Suu`g@5L&AC!K{RCz8fUm#MCCyX9(|BG?rYq)O{u)KM*zhRqeO_'
    'tw!SI!FcPzxh&e1=RWSNDN;rV0~WsQ@%db36GxZ-LK>rE+)!|@Q3Koq;=GI82d*M}59ukMlz&V#Mjx)Nzl`{N&mX7T7sQbGlLT7d'
    '6tA_;B=5W8Mw=E0UOw^4PYfek+7zddJmZQg1Ag#U>z7T3Tw7p3>^Aw#P#{t2i{~n9aAdAO7%VFf^<HV$LmC$vxGS36ny`YABlIT5'
    'YDO()`_F#md1#^J@6t=YqY?TLL@rydOKwmEh|{kxo2qTpiuSpxb%r_Fhw6awBpHj_lEj-|5wo+bY|-5VcdC$AIL9O2V2DCequ!3K'
    'H4M4{q!i%W@WqBqv;n@VH-H@3g}`?PJ?7&`gqW+x6iUIY-Vo`r=Ll`{r@2kG6gp}-aI0?>DM?Np(G_h*Gkd*xc#+r-P}lS}D#qoW'
    '!uSyju!kmJ`f9ze*r`CzY?)+&ZKsh6LZ#zrSfB*)`iD(k_?>cluRL-$<?Hn7A;r-g&gE9`MYfzxv@B4_YQQ?XapmoAT7n;MQmDub'
    '%`i$aO3FQh5<M(-_pJkH+`Y~pg4?~SSqaoT`yo8Iok;ro822-Zk4wnWX+XzXrGT-7>tO_)1@kC1!=FaRL{%>BI@O8sPk(jvN%OTz'
    'V}&(%vU+_s*JN~r?h%4F0u<qxju2u<$upvqGK)hu?(9=o89ZVm&zly73Pk-~OiJ@6C5N9Vi@Q0~e3+pmb(I=BAj!+seDPe0A{j}v'
    '-^K5Xkx``yuy%gVF(^V=>grVPor7_^B0=F=8X-I$6t@|?NtUn%nFz*mq$Cqxe0Dl7@yJ$#lss2XADx)M_5*m5incHFVMV}vUw2qd'
    'Oq&B5C&4y4YdT17IRt&7H9;%k(~E2DyOz%R+n+!_`+eEoV==BF*&N@)o=c?vp6kE0Of>iOCeTb10jg`=-@S(yQB34S7#Wdi$Vv7Q'
    'k9)D}L69eZlBUr>)#~=UF(CL-{93{4?9OOZhA&m@+ppz21QS}djF%b$nnNQ_fEkVqw0L#y|HcX_C%}AwA|xBIFcM4nWt5QEdc@ei'
    'vsfIJ1`00F^`i02P})%!=KrHHD_T5Sk={F$!#MusGSGqni<@(bvHEkZJ#YI|2ycC>8DUt(awOZ&4`Llyk-!g?)Vmda`e8yQ+V>pc'
    'RIgD<G>UCzNG99kbDS}*tOTu3Sz5FD<+6s1_vNo0aVtNyikoem0d$*AYSU9p$O!!wJ5Zd<#L2rjxHtjtY;(h`i^7upnFt||JFNAJ'
    'B_$0gGB^o|``xDLflzM$WvSMQrmqRu%haL~%!{6KE&Rz35EgXHfh}u%u<)+H8lS)xvQ00@{#M^4;v>dR^IA41g|gKLW+}vur*c}z'
    'Vhi|2O6$o`)jeo#xQC6L3dn)uSB>YPcW&}_SNBi8l)0L0ha^RCfDPm>Bx%Fee4RZ9II_mGaO@53%A4VJ`yGji0rjZYNEeMUvu)FG'
    'sA-j=bsHn(PxQqc(cxIQt)pa!Eh$I!IV>XGo+gR(m1uodm|UFd+CDf2xbO|cAS#8O;+$w0pWz>lbDZZDTPaMk2+t(64h%lM+PKE*'
    'gpSc5()F^0*(gw*6aAhgqtE~)o7M@a*_VbGm3!&(7-Bv@>5bvO09RXAf}jITc~VHSsYN2W%h#l1Pc4pOu39jh*Kj^yz;&2eXjy28'
    '{HBxPFK!c6tt%_P9{Xzvsuq<j*a%YMpeG9TBfV`qcB9+?XZi)Rz_}A!gkbavKQ|l<7{$QZyCm@W|55_)LJ4Qs@itE4EuhESV3dQ4'
    '0(1j<9PKw@kf~weckU2xk$axETj6*ggPoxjwb29{v6M5qj=@lk)EcUt(8agzfP0-r(8L)48l`#c1pLLLsOWOe1fuou(j<oPl;5p2'
    'Z`@ROLoCPNocOcGp&;Z`P|kXMzztWw?BdB&Xw7*n?m)iYnRnItFh2&V!bK(kY#XTyZ|6|50lF|FB+kMF7Q(UI>c7GK1Zxy3{Y1lQ'
    '>iH;mKfeh$G-FrnBkiK=7q*9-lFGONlqiILw#UoGvEJ@-)Mt>PHeDk*M6;XmePWsjjasFjTCZ5t=5=k1ireOh7B4<Z#;u&Y!W);1'
    'WG~}W8S}yferGmlZUmScognx!!GS-Zh-`qTbGKAN>6_TOOB23p)6D`_GD~W+ppuwX$uul6`lL5AIk7;~EJmz2)VK2Vd@L9ZySEV!'
    '50y+!%_N`At5<C2KJhdDviO*$hc&`FUo?-?dq^>(bm=mW@`5MA9VPAcx?&02dAJrsw9!=Frsq)Tq!J<!;TU%<{*WR}q%;KnKu6L@'
    'H%^_}d2L|XqzSR?GV@eq+pC5=Y$_!cLi%X$RbM(}4N_k-qz4a1nRvnY+iPOyB~@>UY;X|Zi~<{TB_njs^zWCrFvJ{=(9dEFep>+k'
    ')BD_D*|fJ8Sf!{@be+oT&!yD}AL=YM4g2qJ?t4#S5SB*@Y5-o)5V8C9Z@&OE`|*_R06MAPWp;gGaP|>8N1HBw!<a!T$Ekb!?tHIY'
    'Uiw2>%igzbvGI?sZPP9ci6$dq+4~w_w#W7g%Q%h18qGUJ(U;0?W}t43Rd5GH&Hp$PrIRyeR#@!0uaL4W3Iv6P2Jbsr73u?33I5jc'
    'k|HNU$o9ij%^EpK2CHVOk=~0lOr2ryxcexsnvmne#?J&8b*n|ofQhrFzhNf<?yL;f^G<~(#uS#(X~!nyI!h^^J(*=5k@coRsO=^&'
    '&S|5!Fl+aN;vjAd_2VOB&8+oP1`Lg-69`NX=@RD@G2`8QgTwCQ$@^|F|9v5JI}V1(l0pbaUb5X)3f{wPf$bwzQS5x{gQmf$$)<6O'
    'i@XO43$x_+SGfbwyDGi_LLM8s_c&BzU#qE*6?r7hv%=g7n_mYvqD=>`uH+!m8}fe3)$&V$!%1kjyEft}b++^IHU+8ja#W@Xy04x_'
    't`{nM(@mB7L9p_XB`Q6M7nTqG_H>I(rCYMAHsYkBi$Jy6qnFqEb8SIg#N2@kH#JNgw5>Zh@9N}PSY&}nL1n*GsQLXxBQhjYY+xhT'
    'Wy_V;V!UaFOvz%n6pip^z<yi6+4dFW$QUP_(ePJy(H>5BJhL(|3Gm4Kaih#+v<_S@4zjOLWIC0an_irV=S3<YsI70(49rj7t&(b1'
    '`1BF3@L)g|Ex+{zS>x0DeSa}Xjot%3<UP!5FQA{EgS4J!9b4!6SU-kD?s^v`Wwc!Xp@P66VIMJh2TxsgK8%ujh}#Y9Qz*uu`dXF{'
    'o1}OOv%XB34ZY&R5hj|3r{ZLN+uZZjJ>n0<>rIh0EkPdEYOt!vtE|P-^vF!V7^#2Ttcstz(OJ~rPZ)_cyewz?_jEx|6s>d3)BVR6'
    ';@*}wS3euev!ppZUZQH+TUOhJv>*pgJ}Q6B)E3<typdEtfPrdZLfml`q$SO0dVvK!-M<;FZ=*$qIjh`PTayLP|E_&)sc5{~RZvuV'
    'rJU}KKDm7~Lb2WAPwayxYrib-5tDD7ozJQjS<n4ArC)#=ILFl+tIyQ+kAFOrdxY)T>IhJm6jtH%M?E$8ZrF-vu}vf=sq?i;V%q>%'
    '+(NczRF#})XFsRLXp$YM)XUKp)vrhA>Dc4rk^!Z%rDbx41V1v^tUTJ32`yRmM(B1m@ojir|G56eO!;yOO$he--bZ#L5>gbMkncKf'
    '@-x5FCuYCKtRFTGA3HHfxkr^vEsNsH=xZdbKD~1E&wj18#YIf-CTYg-y7a(R$a1hKIC&Glt=%aD92);%qJ%n{>nN!u87}%E<U`hX'
    '9z+kEu~Mt<a1{282<}mvmJ@?e*P9DU=7<x*)(X{O<Nnctx?ObbRu8oidX2aqhnu;-c;SSHh{H0XEYa?yMbE>7ozus}T|$2sv374k'
    '?E)ZAz*BIsHjhmCMjUTXIwhdet<L&+M1Bq`!kZebiK=Uw^Q@*NG8p`kf6aW?Z_5`MwXaj0ezMnkesoQ_6+-2Hofu{j<`)Uaewh_&'
    'v}PIoo5As`L!{pox$BcL@5`K`>y>tQs6lSLv^U7GZpX3D;)(EIfzw2U<W*w>W@?NOEK546sMpp>a?+d1ooO*x@msadIh1*HL5R@f'
    'Y&KHi^7aoTiRIHIRigL^tS9m|(UI8r4m4Po_IC1Q|EeK32~Y&jq=9Xt<`mL$rLxq{`NT9=&ON?t0X6~{FxyMo&SH&c<1dk*ER5y?'
    '7O+HYh9*Zl0|&3@O#f?~bipB-u9(MPyWM|-8KH?c-W23P^+QCYewc&}(oEBRC|*_l%U)CHLNka>*~@&lrv9FG5=kQG`#$<>kX@0J'
    'V;`oo8d8lq=FC96GGG|qa?_ihj2JW&MD&*k)~6QXOk~yM_@-djpwbgJym|QA(eEOiL3zj7u}xAg%o^&HzGx_79^q9j%L(v#*trRP'
    'O1C-41><u~NdRx;5J_vpee*fZljc%50rrBJ0;|~<Ulan1I^3_kwta|wf_<{o0j(tC5?*1F4_HGoE#|*#Dyrg->TvX^vCd1RDn&i9'
    'V?qu@!Qqgg@bSSd)_#G&aD5GX%))J@NL=h-zq>pMjI~N+vC+0?nka9eLqCyG)52%)H%Pu5MnT)RepAI+ZJcf69S6}2>2vgVVAx64'
    '@<nBf45)t11UOCK!Vc1gK&H9G(2lM$ybXJtN<qPnj7}ns_|6HPxo$=S*BNMy_(|Xil1!FCXV0m8hp6}^_4mktJm6anB%Ij<XQ-US'
    'Qat+iGK>58<V+gg0X0-&^8QfBohNbVM@jLcG{W5`w0sdZ36G*NghYyxLX(}y>hkb;d=u?8I)3x;$75yQsZ&PH45SHIAJ}7|tt!4y'
    'E3B{2GkN5bXA52HgDYQUzs{RgkL~t_j<-<cVG$d&E5XeYNwZmTwEy#DCy0_g7C>Q4TK7bw5g#<+hpQOsA)&+Z5ScTYb-$s*rw@R%'
    'IK472bB<x~)hKH%;~sXEQ`iP;(a1D_M(A{(&BmI_R>lM}akuwNkUYl|xwTBjWO|`s`7)feb&sE=qV!CrnMVk)4pzBGRmKOitNRs;'
    'JtF({0K`lPuF(NRp4%#h1PUp#1ba$$O^`+gbXV+LNUdt_`XRxp(UV#vtj9Z-<KPn|&SUzz)xS-DFl6|lltymEh#7|1@t&d0kTrH@'
    'UT?>q3YV(%-sJrLXIKMB?=mn=rYq)KJMNSD=C`iz=seYi;Nk5b-)N`*c-=XKKN7?xcf+%}$ChCk$k(Q2i@5A*j8>CIsx4r^+k8q6'
    'z?n^-n@pc+fIPuajxd-M_ob<-p&}RbJowo#JC!G==X#uzgq+$X9lX}<bRqk$&**FKS9t*<I^M-XW<G&wE+VojPlkNEEF8ZhaKYO@'
    '>RH~Ei<b=0FN$=Yj{RwNiv)Jm48MK5oa-5tnHS@swYtmVmnC~px7<`!&qpxoz%Ovzhj=Nnx$P{-*jYWybhU<<`NE;%nBv70vvB9v'
    'UgNM>;GUb9a@_JNYW!gc^Fp_%<-(R-6c+B@T#jqTZc?5cCfe0zGgZ~q&(AA(5=gaG-BeXq7S_GGlm<6nzmI4PC&sPT?eHPqk4iol'
    'y&u2lHgKak*ese|$H6WG?y%nxJQvyecs~IYWD`vxY+jaT2gK&ui9I_BMnsRLMnkY>nYkxF{6J15`p=^lfbejH64tj@MteYGz5Y7A'
    '4=K-;@+sUPPSwI?>ZpYsJC(nQ!Qnb9l}c34I&&m%l1`^O5Hs+{{b1yt>f;IkLF=Xb3+pV5-A^0L5kM0Yse7a@&t!Aq`1P$f4~2{='
    '7om8(T*ck3k>fvVw0C3ZwQK>UrZTG(mW#7?fsTSWOU6>jxo5D5QpDX&mE`zKPXf8De>el(&yP*UI)M)u?3Y@OB`pT2kfhLxFv#2U'
    'V$1&ItdqYp+WCkD>3bN67jq%SPEuU{p<PR5(wst0<3+Am*nq?Wo$YUI+98jTsANM05^5MjAmu!E+Og!04x0JEgyDnaJGaL%;pH?m'
    'RPMvHBG{S5I7hBsTU7*=`-v}^N&0fv!cc-^T#N}kj<0#}gb)#ljR@mBMp&fo`>m0n9aRc}6YFouI$5Isy~Y~HYQTq95{MIwr3d_A'
    'DMYbNofZ7Ftg-qj;{C!7Hkz*H6NH@6yO|WN39`|kl6}=7#9n~E@ve|5i?IbM-ksqvmJ~tQkm)c#KT3(~tuH1@*oxUH2Z$w(_xEze'
    'dd~LT>J}7MG+wD#CUf8|DN+&dtas7OBs|7#R$soQQnAZ+Lp{?v`*FgmbpCPpW%Al&s&=>)8u<Kbjl_|(1lhAFa}dV4D%V%i!+C=Y'
    'jAdEDKsZ&r-(1+`Sqk4Xk2H_xl0b5@Xv0zySYEmv&RASEMCcXvLoMqFQ9f)92FZf1c2>|SwheE}^T%e-%#DusSv_(g{+@n@;RsQs'
    '{l&q(+6}4@xJVzHqzFh<Z9A;x>H6oAJp|L?EycHgOwbYq^bhB!;sU?OM2oy$fc*>my4s|W<g3~8(Iqm!oKJZjw99=*E&V2{WC<Hr'
    'D*Si%NtPM$@AxwwnU|U<6RfB96|$@Z7KUpLF$Ee$atIFESNx`30@GqJ9{1S6*0{jdzkvyYV;{im+*&n@f5Y;aiC^kU&t>I$|BV<x'
    'JJ}@Mh)*ev=Z_%D3DSkFCpvDXNsjg6RG#BbeTg3eAQd~<GmYAgV1ci7d$&OHmwCv;AMK?bL9O!W#Jd9!y@C9h+gs+~j|<$i(IEe&'
    'E;6m*SjN)VA;HdwI}=ylt$^Zxhwe$VHO*MAE4mM(9ym`=L%9)s@#L@ibH&kV)ErrOIyb0#8AeRQ%2AG}{%bV0Yo`&gsr^t17q9*d'
    'yR?kV3`{JIun_po?kkTgsd0dOB{yM6I$>_juXb+yFTRf1ld(fFv&U@fSkC@I50-<53aW5%&yp4G{dr6*CFTfq{kX{{A>vA{6m!9G'
    '>?EL^3qo9o&C78fF&w*S<N|T%Cyd^J@4Xljc2lwuRe&=-*N4^o*xr+>Bj{y5J#XB&cqI11tSk&{@Z!ldbfHa9+IY_L-?&$nG|kWG'
    'f%lGjf@vX~z!9fEC=!nd<ncQw*OojcGTbHzHW=u#qgBsN$Uo3<=B7`@k+-ipHiKSU-<TV!Uu?$WIe?L2iJDlQEr}d)+1%^rN7*f-'
    'i=dM!l_0-!)y4kjN9`4s$>#1$?$#Q7>_<+S1fQ}ZF7RB;Ul7U&{b?EXT?3qwKOo#KY#&PCR4YuKiBWXO;wiBUi0sz8o#|BPzuso?'
    'iVW-f#E+*Lw2wJ%^>8TO6*T_s^>1fp4cYJwheK?|+SaFt9Acu^gcnw5G_6dHtd|$0KzD9fJ{#n_P3!tP3av1w)Gs|91+-)J0Ecbs'
    '*>Ml7N~D`8EQ6{!=9DyAg_I>J^uy;Ad?qL#{Zta#s!9plC1nC*6<i((RzpN)Q-q@Tti_E0YR6Ly3!{USKr_51-i5EIgQ?L+I<1qD'
    '+&=dtk^FCRue}T0<(fInIHcES%#kT6m$;;ec@SygVXibqVtBa6he2{nTkK=o#3zIu#KNeP-w`q}6wqDczgqj$g@K?MX%{rR38U1$'
    '%t8Putfq8fssRL#l+McWLw4Yv3h;U_l&>P$7+XEcT3loq9krEXctQ*O@9O$T={`4t?FYxBQ>6h6QTiCodZQo3m7NP%2Dg+MH`1`P'
    'Rgnr^Y$oA<ijw)i7`1l-`}&mZFq)S;;5DNVt11NGk;pDh<Hoo560U<->BEjP+1Z0+1I5ha?acyqa^kYJkk&0Ih_v`{Yk!eEY%yJg'
    'aipD9N1$4Cnh0n1N)Mft_~FumMFBrf^D$e?A2l{NF=`#>>iw@=&}Al-Ef*wRt}X|OUnHL`Z(&waFNGJyLV6lgy@l<dhpedB5|PkZ'
    'E6BXCbnt>`I*6M*8ugGtgLh9SBQhil2ih}eqiY_R0ga_&ObfWul5BTLq@R0or`5^?%W!mAcS?x|Yz~dhP6#xbzB}SGqEUry5y(?R'
    'sfoV}7tK87s^d*J8KQ3jIVGfy(a=+)Qf`kM0@~9xVK-S#u-A!^M+5CFHSM5t)P#OGDIXhtoCBjLT^lvl6(;nlFckCh&7KKz|Hq`>'
    '?Utro=8fT-x?o*)c|OtA_ui5NwAIq?HPH=CZ8}@EGb~i|i5{lfj5*t+Vup}GI#YV8{Swu#cZIx-%P>w3y(>X8>cOtS)RPslcD}`>'
    'WcJ@0v^$Pmh#&#piQQeB@o@~Kc!&K5FlFs3&Y<|Y&RRi`gHLNK2^PLb7pub(v@<h2(P<>K)9W5*^$o<efey)K(<XQb&cCb-oJO_U'
    '-$*Q-7qq6!TfJX{TL`B3^T$5J1RV6&(iv2m!HT8f*Pg3hfOi4f0B)n}ig^YvQ|LwLit-Lxhv5M~oBrFtAYAgk-Nh5__5dQ1V+E$#'
    'o1?p~qP&^M-mivOR*QbOOk{C&1bC3R8EhMuMpDA$NI9n<P_^jd<L<N(>J#u2{x~TV=OB+sLP{XVxPlX;6p=JD2?pahf#Mrvsz;8Q'
    'z=sezLyI9x=KI)){?h8x;|?nm3?o}}_p}`O*VEw2X!rb|kW)U)T)R{gj0y@nEmMX;D5Wicn|ej!gew&HhGw_M44sZ)2TG_S*>;<4'
    '_5LY2p1F4MyD<GEx~}wDMK|t?SN1{NkvO|2yTq9}5*PBe-jwvhfz>0BJe!k^M9!<t`axL~#n@T7tni(EvXp^NXeeHYQa<EFs`3kV'
    'Rt1Kf0}D0*Tdy0xZ>TG0e^;inX9pFtXT3sV)%rSbsi+HKnBX>ZB9`oJ(<t3XGtpC0NoQ3tmJshPQ&~LCgs=)1?8&mWbo(1D^yX=M'
    'mLyX`^y@)lNPT$@0v7IJ8-!H#9ef}@Mdom7mR7k)Spfh;<>IbH%aH%gpP9z2wnc{3x5Fgn5GCsQL9wRMAvTSM$oHI@kCiA}AC)pq'
    'q7nb@EnAO6)b$}RuA8x0H#3wtqkho)j?QD6s+c!DywZcB5n9)8EO-MjCkk+cF*slKsZMCn;Aa~R-wfEh@|`~-RshiAfGp()i>6MO'
    'n$QUGn?Z@ZI_klM<|%PC+sYeS!Ob$RU5<TB(6DlGD8P2uKjiJ9-mH2&1<r*3vOXB(xh8`N%d(u~ERvK6xVc8t$zD4ePgRmKrD|Zv'
    '!n5<}gZfcVs``T-Ja%-s3*$K}4D)54?a5L;;*@lvGn9(ZG0=&Xblu<Oy%6yI%YKR9Ce$^*Rx1Zt3RdK%jcKS$p?tYyFB=BM1{uF{'
    'E9>vi`VC=T1L*bB_qu8Jt8$+i9$(KSa45S+08vYym;V0IM?0Y=n&=0X*@3l%^K!*HK5Hk*z$=VN{Zs2?o+ELZl4M8rKf^d!`2X$g'
    'pPt-lh#TjsrNg6<yjTq4^otsCFv+cJa+9X#oGzZJfu+pDPadFS{UK)I5AB;_x~@!DlHTQ1HEUqENS+!N`z55(UV7yjlBdM(0p4m3'
    'ZC9G6O-a9FGr+T-T)xf{1=AzrS=;d>O1}Fu+2`Fp28pz7%vKSL(ImtkV#3kweeky?cpuv(tYy98a|W=i?>4g-Q=XWyHAW8_{0%>Z'
    'ee-hMHhLQG<4wAl=G2}hT<wNj8<(wXrsz+hAMxCKo4n7%LMmt~9Q5j{cqOj{?Sq*JiPxSLnV%Z?!-lFLGri-;6@F_Z^UIy1<27YV'
    'a{<-_$8R)ZYD(_l`YP#7l5^Zi5+M$-7A**chV+o86^1_%7&7+B{+SxmbJoBx0WRfXbcxORzuBaCi4^V&6Tn^sT0Q|RQU_(DGxGXP'
    'taC|v7QhYP`Q4i>^}bG&NQQkv-^6zDO5BZDUsFV?9f0#|9RHvid@Iv0<Ay>ncgGXXSNpcf{}P@my-;{rioz|R7Mu4w2r<nwu~zzn'
    'bDI7(FaKL*)-!B{VpAcC(xsgHG_TW50bM<8mmSWzpuM{1JzC}e3RUYTR5s&GyhWnYAQ|2}_Ki$YgDrcA+A@x+W_v6we=d06BTM7Z'
    '{nbwLwfTL_>C>bOVmuRb_F6F5Oa)iintS&FRhF1h%h&Umpz1oVH8!1{cJc?aAFt*q@%EBnhze<&FvWV68hdgSLd4u695ph$Dnu3@'
    'YuMryy%?Ok_nKGElogCKXl5OYcx*ApNGe~@V{TG0@9irdeUA>I2lHUmXqp$<83nqA)YK3m9`8A@+ALyFA15X6RUWeqnyo{ZiMK;h'
    'szDv2mefC=Z&M=a8k~5F3}L#@$>3$~D$3fKvKY+S=<8c5@)N{rR@6ADyB>{MXmC>UH-chDx;U(7_T960yGcH}N|Csh2S$M9(*!s|'
    'l^^9GJp%RoZGeW46b(rUknYd3d8?-=U$s54<PeIxf7Lac5%`I`)ocSB0=oy_Z;~-H5%m$KEs`A`-c15pUzx{CwaljQGo&VMDjK6M'
    '#UHwJK>*cAa6rfqFyj{1j&2XPFMP%ha%8kQN;xTtTFg(5rP%(0P0!-Zn%$5sl~ak9Oop44W8g#FBxjjG-w7))Ok+NTO}C#QyN^UM'
    '9qO!iI)sTU^swFll&1$b5eYbl4+ECLpW4SGGVG!Y+TmIYzh@tCj`bdkE@7DRb}i+D&2?LALhZPU(c~_mQ~(Z!dvKh?G<V+Lz6os~'
    '<KwY7AI0A|)PeIaOpg$pnp&k8MFTF0@et^9<TT)4dOnJPI`bdCdpv@85brlC;+BN5c)l^_0<z?fw(BV&wk%vmagVibw%}O?x_ucT'
    '&aP-cKT~9bT8#Z~3Iz2V<j?rK!U00!qOF9PLzHyjrl?mB`6rVgB``&vZIbQ}wnI%Tm+p0R-&r4Qs|M?H-PXeFR+>FSQf2=7{goAe'
    'h+WvDOxf=rG=Q2x=?=()3k#qB0pS7r1{*@gRjQ629R{$LyThyv?->o-3aPBYkBu4fOz`oX-penrOh%}7Hspf45o}xu+<Y`EPrO4n'
    '9<L!gbhB{e(iDsVJ4kwGoaqaED_<1@al~GmEz5`#x=PqS4XTcOMof<(ZF`Ck3a;i<dEsbvE_jx&^^H_kwF;$~BV$j^+}ZX`JJsqa'
    '(HX8gx$D!1Y62(}pN$XkW-*RU4Nc-it@CRdL1oU(mO<G68C-X6vi=IWrhfkPQ=LW~{YKXBB-(~fqoqL~X2i2Gz21z<=EeNn)-I7`'
    'j}0BW=i#3SemdRm<5g5<&WQeH-~6WYvbhn3eG(k&V=mjV=wnP>K*C{MXs(9nXk*9>eJ?o^u|vGnO~&92uFNob?~Wh5KE*H>p`$a;'
    'Jq?aeT9aBQ+Gzoyf|+7ixpE?HBO}Z|Cw(i!)mMIvy)Xm-cm(Xqt1@n<=R71p$ARe>+bq;P^g_k_Z%R=0=7AYF$TXpPh%FO)q&_sI'
    'C(UqadfSy95TgFBZk}Jq<nad=gg!4?3e3_E*ayim#M3$PiT0(cfi~|xFZLNBO+43BA1iChrtYk&l0~q~wv$KvLr{CtxGN%fnUZ1B'
    '8~VjNo*&BbX!`$~QhxfD8|x$m=K%r^FR|O7ATj>vriXX$=0>bjktE3~+|5LmuTIb{pF2Nh|E#>@kW@i4%|G^3nA(79oe8!4YuRFK'
    '@HF?W5P)%vtq2Vx?n-VNbYyu)bOi|1wtuV5j7_J~@s0t2Erhn9Kx|Jsw{pkbFq1RewnoM1qD3$m=`wN+0YVwXwVh~o%FQto8L3oP'
    '`hPD^S+$**TwQh(2CDnSHDQ^h)?7#I(U91KnJ8LqdhT@Fwbe<9eERPkG^zW%BDdG=$*1+l8->{tkSJs@5dokB+^7snx|CAD2A6pd'
    '8FvMDU}0VFkP;jrp9xRmm6t9zc5U*YmFA-Ik>y+6!KoNH#iu=)*4A9I9sx}!RpOa#YssYX_^Sk<uIoQg6pGjy=DM`~Hmn%HJD#Lq'
    'ZmGgi`?iUwoP(!j;FkSRPr1Z{Me(&)E!dhY#|^g8<q%bW;GIxAa7z1&iS%HJ#a|aVveCR|GjS|N$bh@@__{Nm;9KjaT>1Y9So+Hi'
    'OXj$kxJL|CLC|5ELF)bM6tElS{N%~J!oBb_G{H=yEL_am*KA3+L*SotbxatzY4}3!;MDzbCb*AeY5{1X%!*8>mFSE?J7sY;QovqI'
    '9Me?<lo~w!+Z%7e;b?gI=BhiB<i|;^@$G2w5_K9b+urgw{fiwG{zHPw_BGuxk+{Yn{B9`(5P*1pk@ra1j$gwfa3<xXLyPXr(%6jO'
    'ug0z5W|QZBGc<?S?mZiuprSY}iFazZ%NsBk=WQzjsK4_vo@9uot|Hp5Oxe7?YW0T-mWXXYE1!7lvtP*++>ihO-(uoEU7YdngE~pD'
    'SC-2P__TYLS!914#X8qG-7KRx@{qo7+jO)&6kT*0mkfZksY|N!pHz|;bOGSYN2ay;Zl1F2AxUGX-2UmTL>8lFkg3qpLD#zXr3WV@'
    'v@wXb#ey4~VIn1$IJRk;xgss%7{*LjZ6qb8<W0*R<$MDCGv~)tK(>Zo4^(FOnvYX#087`6Kd#{OWcWY(E1MQA*hDx}0+~_{W2--I'
    'eVXaUIb)14$K2<`$b?#zP8Kr6TZ{)d$!hA~t~`JADURqlKDQhq0=e?kesFKW0VmF%#YQfkBoO?WegHj?<7G|kDoD@Q;GtB-TK)K1'
    '3e|^OH}+^92i}I)zCx~MAFnH-$95z>dvvKFD-s&0b5U2MUFczn)a;HgW6FXC&n?Px0{>lz$@cr`>iim0z5;aA74oR>gRN3!tO<CT'
    'ZHbb~Ataq+iw-v73d<yhgy&Dm((8G;t#qiBJ_9;9P)h>ZuRiv-<jzhGg@T?cmNy%9<M5z+>k8(MN^s+ZZsU5t+G4&=p#G+zrtz}!'
    '9i%LX0c?NHs*!}Zp_G%#>V1MoQH%D(IQeOMo$8T}&d#IEamqQzgMNcA<Kb-#<nn=8#lil*lviZvor@qc9Df3*g2U-CD}47Jrxu$W'
    'zh}Wqs3{E`J&z-Dq~`V&U`<K}Q<&s6ZB+7)XtGI?-@bq;Tsy4$&XqgyD$aC|e3w76y^&TpE&zC_GFi<{0uKES@4N;;!WnxY;d;Qh'
    's?aJ~rPQ#Kt%<=<VN+Z}Eh+v0KPbt_bZAj=r%Qxb2`d#LA*|Lt+Fp~mqI-d?t2jDUkJ!=K`Es&}iS%a-_e#~qz6Gs!2Q*CiFZJM('
    'JUQ}J(CA_;?wU$nJu*7@!+JoidNL_XARfvf@XQ#hZ#fhE1wz5jQjzPOcRV!Q0Hqxx6h}I$x1xWhZ3CG8T$njbDr0oBLT4NZNTnzx'
    'uT=nMk9BQQh})Z^%{!<nkXChXGlaox{rx~j`6p-^jff@pj>f%X<);V6R^K(k<(kS|Muf>G7B{s=-&leL9qSV6nz)Xam?smHq1-mT'
    '+z^`vOF?%*u8O(PZ{gJ8GZWARO2ZTp4u$YvT?Rp#9TU_NB+rC&`PG(Z9%DvI9e&3=He+$jVJH(I`d)ZTutkU>$e^^J%787#Qv?%+'
    'cQ#EO99L8tG-DoFJfQ4=%Zc`smZ?>B4Rcd+CsmOJn|gF)%XC5xzkhu6QkLl#zEsa)o=8aIITx6;K#9ouNn`l&l;}vjUP92|S)+zd'
    '5!<Ndqrv8p>{c)UtaU9PjlJ4>7qXv70!%E%&nQ>^9Gj@>k@;@S>i|3tkBk-ejRC&kd5C+Z)Nn0m$VF4#<uttbcKWE2Oa^?z2deuO'
    '!|Ae)rp%k6-kklts{V*3lVeeN)ct!NG0Bn|m*r~;93^7UJTi5=TxhA^6y_bEN%ZqCK|Dy0X{>~{ceKR*PN-x%?l+5YYk_SBMKrRu'
    '?VP~YLXSg?5jEEp;GNA;{oGcUSTb8z-dSNqp8!$3wnK-bbl!KlQwy2X9{`F;j~#6Ye@-P3KCqK{qc5u30$e~!<#*4a5QGynX35TT'
    'P8cIzrCF&_)>q^V+)_rodlllk>5)=JbCxqKNggPD;&xib@9=yhA>pg5P_Di9wAiJJYoOjDfKVGp6Ur4PZqdeN{T8tVWT&W-IKtJ8'
    '1|V^QB?~*p8PiXRY?dqZo&JR_pYDQxK+!*QVBM)B%+0X&(k=))lTE9Gm8U~?j{%u8(}0pvPn|2HN#TiaZ_m7LW4$#!s<Cul251`#'
    'z*|yr3ShHT;kV3UKwVM)yW~PLCP6Th1IewConts4LBg39DF18>Wbo3?Bs{P6l0CjJfYcWCl?vD{xQb)4%s7`j&D{S6KTzv_;j&DF'
    '?Wgi=h5)Tv?Yp{VD?e!$sHbn$Hs^i!zDTFY2-tlBtRBP5RLx}Anm-#Rb2fsCZQF5|w;`{%`!YFg<FEU{D>=i_TNwSZLE{RVg4fQE'
    'Wdblda1P^6E@xG=Is6?X%k{>uQo2nzYF+He`ZgcNnz_V=0uONvLIX1Ndz>fi6FkMrEy8>po@NY_0`vXP(0)}ogCqAf&)4Ktk*-SD'
    'Bo!w4pQg<^@e9DF8MW@}7ti{NS>@7~Dv}{5BE8{W_U1k27a6C0GyOsJU7kmuksaKZ28>dm{mu--I86bCCEE+iP|Bq*UYXEx4zM!E'
    ';t%aHC|H=Xk4mJ&wEZ|pP1+U|6}s%+9HS}~7)N;;D#C#-2IUEYe9Q|D;zXhiwc6D&_JF;62e@ARMAn7WKIQ{l@#^vpLT)5eAEE_G'
    'Zisp*YXi~?m$)E(y`c#04+v$RdBu|TeQKvP5N`5c0Kvl(W*%0wzp|qpEZ@JU_QspdTkHqXqYYR6u?Pt%;q?Sg+gTH;#MMvcMbpDs'
    '{JnHoLoO~;3M>h#*JdQgG%dFhg12Dl5jG)m|4BdmyFo)xwn*|Xv$u?KcOzMjD5zHxI-Jmn%8D~%w<xaBixKPFTu$131&4x!#WKE~'
    'AmCL8f+93<u{m<SqUt_z>1T8xj4E^0AFWfHJ}kKop%5TO_CmS`hiENN^$^yBz3>H?GE?@$UsXxg5^~bt?iXDA4PJ2lY#+X@3=R2x'
    'OfeYQ)_ZqsbQu_qcCT7Tairee{VOP;OSWP*_L-F#&@STtG;<$Y3obAi4};!st?JHaf4NWSqer2{<7m~&Aq7&q7AX2qMB<v9v7X2f'
    'ZHx2yuz!$0y2#PGl~RZk#D;wX<IoU}jr|=cA16q3z5?LZ80v?))ZNM)8?z|u(>X`psd{<5Sv?>NY2y>mmKO3j)-s^|39=}fQl%g3'
    'R70?Y8G_Vil1ycY8?33sJQMuV*DCs{d|eT8L*16d<W$$KwdTSz!b;2(vYjHk)Hz*Z>R40YwArey9hhB<dP;vN>q*HS{DU#Rc@^lK'
    'L33T}vz}=!3d_Nx2@Z!0s~qXk{A@L~8-h#Vq;Hd2n8I9Q!j((j8(my|cEtf!b;plt3lRv2+t&_4f!=53=1nB~u*Ym>nWXcYiqo^('
    'nv$l#tsSZ5Ln{e@0N>Go07I@l1ZHQzdL?>N;PaVZ=5@6{)Z@g4;`KUW{^1mi^PxdEsB;z4<uB519Vj5?l_6Vf#-f`u%Ia5pm?iNj'
    '<zEj8v<C)aeyieTOnGvu)ixaWR0-`BNHv@z6B!ZXkRW!|=wgnIm5@f2*p%+;eGfAxy4-OLFCHi|!;ck_5^Pq=2^`xoH$&yahM>P`'
    'IwAoiTqXn{XnF~LSW@hC4y3)=5E^J>%}l8ZDU*4X`OCp5$mTTi8Bu*l2&JWPqmqx37$cAu0vy++lVDvBqJWB6w~u1i3#(p({|pWr'
    '2MRF6m8es<HA)Im9L+Gt-OE}E0>=w1TEjrb$d1+dZfno6iu*2v73O^89P1D9C|<5p5cLGMUy&zwiy@*c09RLKnia(!2u2Gl?gVHG'
    '%J?W-@?!y*_`2|RY*+a&DJ6ev=~!l6LCuAyJCYOFMQi(}nk&0~vdZ%!Fj+l&xrFLUXX(~1Jy;9}JON|zR9$J?F{gju{AikD?5-P@'
    'yyZiHHUugWqyn?Hg6#?c`4!;9mcY-Z<-hw>P6uyO7tj^AP$Mvi`-XVbN1BpX)*qz?PK>)ev0;@q!sJjN5St}@oR#ibkW)Kr1VK<8'
    'Gxk7xojeB_uokdg8jFVf<(YvQ%j0^Uy51wrobK@r7I7p~bA`^O5Hx>+n-22?zPUo0=Z5&WM#92OEPEH4j}vdhXOUF<){Y#mFVM!s'
    'iRMN}I(tZ`2v?h&ksZBpefTl9JadbMF+<NZ*byvi;9ZJ-jU}}%Zl%N-o2M<d`UGIf`JVUzjHtR@uDM9;QY&*R-QIG-{kApRk|fru'
    'j1by?Dlktq!cQ=+dEHaEH`2lOg!h5w<#8G9U5r;=;Iu--boZYVo+1H)IfLLhxaY%aW-`0Co=j8`I83NP206=~I2=kIg0a0Fgpr<?'
    'v@fRQ{~R|fbFJdat_KkD1j}VY?if~mEPB(gf~~b?IyzQTp}me~G|SF0AU^B_I{Ft-&<WXn3%g*$_w9I%^Og;k`uL-L!;{DHXnp`*'
    'm*OC>c~cninu^AvM{^a+?iyJC|228x5b?dz0s_dIv1?YFoD>oD(C6I7j6G=hT>qlki7el(O@Vv)_M>T?;Veb6vXHsZU*G!X?_Tog'
    '%<$Wlv#w#w!9dEX(!+i80e6C=e?T>074=nw!brk%C-8e6bMoa;2`venw}bzY0Hjm)i~BPEmj_mBkWx+y1*Tzu;LfnXwOiSTkJ6#m'
    '7tTQ!FM7%XDCOiUl^@?p5unHXu10ht?~CFl(Tz0LFWV@%+aZ$6VEfjnipOD5B-h*g;8yfLCY7r|TNM^z5;<Z^e4*~V7Dxthr?W_}'
    'B8I(-WyGs#HyrF7@{?@i>EnZXRB>wl(885nE%d4FZCe%F;h}g=bDp6$yZI&}rRRfSZ0!^5ki;rL38AfX8^0lbxNEMBwzAc@%U-N*'
    'L_e6Z-#L@5e)O<K9wdB1*gYN<*%>Uy@0Sc0@xRVl6-whm80`o+BV~La2kY{6=qfNt!3btEsfiKMK<S{=01neJwzzsb6u-{b_)GE9'
    '9)Iy3=-&|}Ks<_%C&1v)cE~q))X<D=ESdmd<eN{#1NTJ^-|$*p8@Q5p<WUp_rXghgO(NUAoy8=Nb)q|ZI32*wt0r(j@7+ae0fsk='
    '`wkys^<VZ6+EUFnwk&uCYs~<q46IckdUw7cn#2Ex_>J+62aajeYEf!&B#P-<?a2CCLPl>P2RD9uBNviR*iq)f7FPGCHBP8!V5<5n'
    'WDqv<Lh37m^44@x!F~h1@yy3}?lc%?+H-x6%MMk-Ah4O37M~UOc?fQ}%Cr1m*|0Q3Ga08_g}GU+K~zsk&*dHDBzSth7PxRT9bY4n'
    'eii*gmw!VwMkpyi$10+$wP1@=7l=1e5X&uz)-5P)iFn^6Ydq2&PYN5Ld(k|Cb9S||(mVchVOfGv^eB!-5;l6H{gl|j>}KOd_O=Z('
    'Q}jK6v!Xb$vt=}_&Bvl3J5yNXcT&Te;fBNtpX4_yb8qm##!pIK8LBXlwHmd{7#N(J_yWwKK)mlpP!Lj-aNnq#+THARG5>iFl|A#T'
    '-rxaJc7-YE`*?hv&9((^d=)c~2#HA;5f4>b5u)map<;a(@3#qahp>&~1g8{k(XW1et<kCFI0uotup{O1$c)J**}ow+y`2_`U7WyP'
    'iBlm!O36Sw(~?%}{2KATX59Zxor9i5QR?IYx~1cq3Bs&I!Z;u1R$g)5J4cvK!#fc>>+&}|5J5OYtY2o(cIft8Ze{u|Y%FpmX`x!>'
    '=vry%V#qW>S6}T>qMOwlUX%2qYtfAfP1Ly3y&Gh9Jp1jU;5BbbKtzRh&E001xWucl&0KSw*`M~?xH{zkPJ5;r<YA*6*o`qT7_r2p'
    'm_c}5{P;E_Trz)%SNz830KVoC$iEa&nHIy#1%?H_Dvh%0>eYhb_v>r)mtgD&j4Is&k+`kBKTC^}<HWj`>#m|!eiumNs>d)ibVLBg'
    'xqXu&%QR&Oa-{Vfmh63{hF0#E|567r3uFKgCv@{OO79i8VLCZ6sD@dtH8T_`G|@-C8+Q8!tF$fCSe)?l`23eL9UxGALz^Us1=~a8'
    'QRtoMOZgXdNOEAoAfh?lKSc{7<qfn9a(Lnw4k^DDu?`!HvyKnuay*|f_&u<UWsfw?0p5WpX>by=M!IFPHYuG}v<m##e83`I*u7xy'
    '(D+mT+3vxRyU~LFtA8M`XghKp&i$eTW8Rjo%to6@8k>%|gbxj9-t3|?rugy4e?dYg%E{3GOAk`mxlf3()MZu#4}t0*WY|Kn=`Yk4'
    '$yx#rpVJ0KDwXs(R+J?gU8qjw1zn-H3VV{{VFJWY7*2p9Y9A!_Pe{VyC6f#So)B-&8gttKvg+3RFzCd8@Eze`f=!ym3ld~j#k1bB'
    '{l-fbvuX^VCm?%JK)8-Bv;_z$(fK^k?*$~%B)~b1=u-;n3F=<`#!D;O;wE~}4=xD>SW%p@)8k&d$u7q@D060atXyxzmAI5bR8fT*'
    'SWg=Xh65Sj^uCf6gaGGUsxn}-KHs2039WFI>UA?Vm&R$P7Qoatcd!pcBfgUdkF;8M>6PKhpOpK>?Wz1fN8F-?FI0;lqdF6W7cZ{h'
    '0k-uEtcUlQ6mTfu?_p;gth=WDu%h<BES-55k?tg!w7^kGGsIh{AzqT4E5`$$PgGC=+t31=|57@aos(z}W}*qRMu*SgRw^TqURTh*'
    'cq3pX>{;_S6>(%3Q)P)zKH^f##lq#^g&3id*=G<Ux_QwSjK?J;&s1P-ux<hg@vOCsx7nI{-fB}+VFR$HRr<7`Dex`7OtPcnWpS~9'
    'm{S9bjdXv91cO~4MV{J{W2Y^;q#vyU2RNZ#3lmWPv31J!Q_M|zH;W0VK099Z+9r~(RVQO}Zqp7<>CX}}uLm(Sxy|-x3vwS1Ve4^4'
    '<S@A^S^vC+)49X5XNVDfOLu}@26Abxskdm5V0g}degGtpus3^<XUh}0@kWa=Rlray=VV#QHQW|=!m6bMW1X78_nbyEEX%@Y96cm6'
    'd~?TTmKTK+B+DX#m5_iNXBj~_v?veWg@nbIA|xg~HR>x1ryoG5(bVf3s}SdueKMU^CYU;#f{KI_2+}9>tnYO&4;_ijWw3?ToaNC*'
    'i9c^Yfzs!ttQE0nYf3%1NK`8H-9*aY+AGGx`VFF6vvF0@bQ@&x#xowHX(>bpq}9yoYRehv419W=<T<b=jS)B5>J<Cc{Jt2fRuww('
    'lR`L%WRM)`esd7)9VP{Q9vH~$_LwO|{1wBPEk_6+N{^OXG5f2mU}8$tqO2+24?wtE=!OYV=T{=@?)kG5yBj5$U!vv!lma4)>FgFL'
    'Fd1gkc6)BO6H%XT@>!kuNXEg*=i#3*3%3r&N6^+CSOC*SMBwWGCgGrVm|wb;V9TW;t?6|SQX)clDDAqd9<S`}luHuPGZ9>m+<gFr'
    '<nLyh+q^qQZFLah%10(f&x=><h~Oq#v=_dk_)5a74FQyEtsL1>f^azr#Cs%@rn?LN142>q2!&7F@q4FcO-=?AZYZ&N>6OfPwcglx'
    'v5b*fYTh333=E6!&_3cp&DS;VzEkqdB<)^BSq|{Bv9p}(<wO$u^f?Ylnz+Wbc)}|P%eX;bnR#1?)*e_8trGMX0-ih;v0qz3Ut??='
    'PINkVQcUnY9Dgl`jAKm}A(c4nba`K%nDld>z^r1zIxp<BI>7+NtgoLvq_*v~lN?=~<Ld5>6K-4GP`O`t;X1=1Ymt(LOQkXL;ESC3'
    'k2dPFZRF%hgYQI%UI#N=FpJOOI_tFf3&cdR({uY8Q6-D7NQ`&D{aqQvX`<!RiEP%dbo3^?)OJ2=1`YeW>#62yzH54~@l)V9R=3nR'
    ')w-P%n2_gda-hij-sJ$~@<zs_E~?Ytc#_3?02<5>g420`Gf>v(8>{`tav-D89jy@X2vz0Wm<hK>EQ{OMpMu+J)F{DadRl6KxcbuE'
    'VB)W@rHnA$ii##$x%xQf`~GjA)m~y*Nf+8wIkjuIPejCI&?Ep4q;c)|0TxdS{1h@+8b7XD8G>zSl+ysg>%y^PG~bdoBXw7!q(M~n'
    'gZHVB4@kHHM~-3xwTj`bM*;sD9%w)B&Wix6iaa`Q6UtLjlsSF-!NI~Cr~ge13?X?YS@yoSd0L2Gg3}yqylZ^0k9(S_Os?Qo%WsYk'
    'TsZdhAV=?qE{m}AMAal3j>}`W@oFK=KGpOr7TkhsZhgD=C4KMgpEl9amuKa{`$PnrKP82h3%)kXGIk^oN8BHf^aYbpXfFO{1IO$z'
    'a8}o}<gCKW)vivd3C8SpPp%@0-DNgxMNM*?yJIxyvbQA56v8vsa8X8TfP!?-Vhtx%#kC;9#@f~>zs!>=P3IpJ9da^kBvu8pO!9H1'
    '9XKc)(*fBfgR@J84Wh*k9l~t|nMV4S!W4RJ&5K6fpn>b@U^Mh(%bAx^1`8*&MRm?oRgKSERfjl4-J_dUt-l1MHRxDxqY-VoO#nEE'
    '@gC5HEQSUl=d2H^F%~Dzxc_YOX+j-19!l7hV84@X%Mrfv^4q6%2pCdbrf*pcqmSEY?1rQI(<O$wJe9QJ;AiTogZMs@k5mkZV@!k&'
    '&;t#7r6L&h{H)6oK{HRdFtj(lk>o+-w|{*;vuqL7kMkJ%C+nZD&e3IZ%sK76NTKmA-5^z1zKi&ccn&2E?qfS~^1SVv!dV20pQDyc'
    '3q4a-m-LUKg5QOqBkCje7`^#s#gqF>2W-m_@KWh7qsvo^TK39E5NS!hS)Rxl{>!6f?V5GfI`Fkhw+4{kgHwL{t}uSI<^Ru?HHZr|'
    'g;=mFR3Btv&PWSlH4j?q?6|i_p)T#+&VmdHoJmkES~m$?sl0AX>od8f3pbMmSj+-6dP!f)DqmS@fz4%8SC`ysT6|rznxy)>W}Y1-'
    'VCf`17N}8l@pZJ>QU1N$d4Vd}5>!e!{P5l>1r&L}j)qc^KIN{5NQ;uBcwSR8)9{PTV>t{hhlkR$7Fr|J7gp}{txAa83d6NkEFtU6'
    'ZZ|Vye~OET3@6QjG!5XfWvVEw^3u86YY8x&eJqw{D6<s-dlhdsHKs}oHHeF4e@l<nhd_EiSiv1($We+Y87s_7ZkM?s29oros`1L-'
    '-xMjV?*vsxM2B<}qP2qm2u7C$9H)BL{x9HO5ln95&Dmt-UJ%roDMfyuaq5O7f==^*<64PeuX(&0oVkza63(IxPZC)j3u23Lecq5)'
    'alD|qI{2B)t)!>v->MYvj4o$wgW)9`Yxzgw=z-(Eok@)jkS^@HY7rXgB_`cw-GcKXwM@R;pp_N$dmrS0|IeeVY0k#2zFMT<Wf8m*'
    '7Irb#sUs)|W^L5w)MTeKQ;XJ(STx2+3XPYU>U}zVDY<CQ2#t_^B{=6V!&jOfXRO19`~10cG{>#&+5ygnAfXckH15iMGt;I03~KV3'
    'R=iszH_{w)X1|jjJL?A9d|hn9?cw+x8YA8h88Hl5HJg-c0?lmBqrd$3uWtci=i7o{TD&*FEtA~E0S|Ie_hT>czAlI|a0NZ4&YAjJ'
    'bg`OBh`8#Oe{p#D5j%F&@SaH50ur*fj?u=rK`yXZVl;F6u~Yz*R`BkYZt~78#*05xvv9-MI=u4^rL*ZUOPZuVl(`+G4e11K6>ssq'
    'bl}|c2L<yrP|GR|vC&Q;6%N37^vf$Xadk>%tWu@oH@K!v;0dE;)3GVBn<UzbBPz(12$i*m>M3$4FzTCwy=6otZS)tioCvi`$f=c$'
    'Fuo_5F@*x>-cL3~{^QBAwB5@34fi$aMa>6YMuyq1WTI6WN*&BZk+aHz6#EP+@cb7oEU2`7K|k0yWyJn}fB+y8Vg)q$xfb~%_?SLE'
    '^otPV0KXKzWSZ>?9ZQJzc}v|PXV~GCiR97`ZqZu<>H7AC1_6%(R(;7RMWDlUw;fhjluWz!pn*5lUzF6oiq-Y0L-dDr;B(V=@3rhT'
    'Hz(I*yz`^I%5CpKHFfou6r>JaiqA74BLJ6i4P)xGHad<X1F_w#CTbnX6n`-G%Pgxk?zs*!A!0^uOL3Gn-FV^B2<S1MN!xO8>=y7s'
    'JA5CzP5KGbW(ii_M=&rgr97HLL&BD^n?ycKb3$NNS*k>Aj95&|=fh|I(ntH!yh@TejFX2Ie|yZG&EZA=0XgQ3-C)vu<YvP5s%2_6'
    '0hg^DeHhAIKV%UUP(=Vcrxe&Y+RDD*{ATs{E{{O1Em?vKoysdZr;pBqMyMF`vky_1-nZ8P=vxrvFOj1D%T`J5esm{)`%(5T3kfh+'
    '6$_CkjQ-N?%w6dWM5P%C)7)$-PhLi#OPLH&I&brej4LCIs&54*)A8o*ilocIAz?YL6251ShQWir-%;+F{%lvwJ-;t3WEdO7CJLst'
    'kHP?|hczJbh`YC!022};zr+@0(yF90)tbx|9)PnVr_e2K;=$9nyJpYpA?PpGf+fVU_67MBaggnOL{&GO7c#+CsrAiri`)CbR?yXG'
    '3OKNjsKMQ9*6jBX-hQfLr1PyX2T4{tVShdHbUR1FO@}4Wz`tg_f_>wjFRH6&>r(2N(xz-CtSix91kVm^rIVM@$fg3_et|I(U6WR@'
    'm7RalSVj*5-~kvlT7dh|L#UunxI>QHul5|c5)#>td=N-ov*pcllm&y3TYPW>I3H|#c0xuTc~%e<AM4rOtUzr(>IM&O0>kI_mQAk}'
    '2U|)i9QKH{Bw-mClqA$>uoUp|1}#}1Dq7rgtN{il?h^7B<nlo<`l_J<)Bi0K>!wZQjasJ3A?hR3F-b64SIbPpX}n^`6WK`LLE!gl'
    '09_j0U&VqUUk77Xh<d^VK9vri4xzUG^JxI3xZMe={u!(+nwE+J8Lf^3fx@YiPRCHiysr%_HUo0xXrvQF)!C|RkZpa#oVb=*418ny'
    'kn_^zuhC*nx~b;4J<%i1Re45XnorYqO?B1zdS69uAjPJPDu#>PA3^<nDEH>YP8`sdgRv|U{1neKTKXrOY*^i^sW{hcqSt3sd;vYm'
    'BQ9jh!SRcCNV;8#HPp^SZ1I((V-l|T)ie%u@`psoXmjmEL5tY*4>9LuC!%!6;x)X_Cvy;!Ifu*x+HTscL<8#>-R-Hl+Tis*K~C-r'
    'tP1ciFqN6T@}J?+iKH&>cq-?f`@ZADU7Z9Wl5{LF0-QX%`WyTE2oc9yEKgL-(C^X5^9_`{@CEh_466RCmtdOrWm5bxNzT+1S(~2Z'
    'oNg@YcoTb3>5O&@Y=a?^Y-fGNHx|-Neru{Z?GRL6nXK<X*e5lSjT?QfNFf&`kkSYcOq}YFzN^!oIN#nWu@p9?i{d&^MH2yI+DnU{'
    'T?ukepr*qIXo%0XzEnbDhtG^Kr6aS3?tFt{%>Z~<lN48&w>me5i0Mo7+&rwXr3xYO=}*HrLuT4D(W}7t)gpx4eck0y#C6hN*<VeQ'
    '^(8Bk-5T1#lQNdL05?kEn;HA;<9ls22dR{&-6YFc56*`eVV||G&Y$a?MfofauUm^ZN?8-fA=TuB>6#DtMy@9_0O5eo%>ka@Ee>Vb'
    '(@geDg-zvb%sR|=_46XXs%P_9@YpOyMBq%Mz`hT5*{7a{DWU%X*fj<;B4vzhs{wWurv>yQS@$64+q1(ogY$k(^8=7lNuDpWV;M;1'
    '0q=#F+WN|=LzBa(EaIAtifOg=js5h)4bVwl{XA>0nZ7jPnSySg&AtjLm{r1MBR)tDW_0VfX;|Y-Patk}P4<g>=Y`o*Q9)uX4b6<i'
    'HG6HoZKInQbjyitFq}Dchb0c(bgfHZGM0`q0{_*}gALS{v@uG?c1t}*^*^ftpBunbQa%#YCU~4V3^WERi#<&lDO)({H~5{SKdA(r'
    '6IqsCX?e}!bNO1W5>A8M+(@D&qb8O}<6Wv(gRL}yX)3AKhd;gHe9&gkxIPD}$oRiB<#A}^?r?^|1)~qUl45Y{CfWQM!gGGct9(eX'
    'g15MJ@<B&=npbs6ketfFex;(Yp-5V0D_RsMk#Axrr{ULE898G217^k)A3Y7{2u9ehYwF_}FT)9FtEeD&p}Xc3M_o7wzLG{chI6`8'
    '$9{QD139Yk8|N6f?3ukE5gC)IjChYf4QUY@;C?fmwL_T1VBRtfRfx~d$_e4g;%{)UjBglpu{=E=s10{gGi_{qOQAO_DZdPBsqgCE'
    '9GRnYO)`KB_Z;_cQdQ%-@3igpQ$aK2z`E<PO6hB0Ar~{ce&ODn6MoH8%$a?VXvs~la`S)Vt_b)VJonVQh$i<W<Jkpuw(wOuK1SGV'
    '=6@V|Z8h6O<wC=zFp6XqlkT23@>~D8yo)VE7kfJW43<UY5Q%SG2acV;iNuG1&Xc`qrl`nx*Ra*|O{N%c17cT5X|IDbR%X`P&&s)h'
    '+)|_xpn9B%0f%5Zsvx4qf1~37IB`btm!X}S<cbZbC+;^0q)E79X~{cajEZ*_zeJKNN-hRvQI|I4b1~Z@L5JCxK%|@t#yV9;B(SWE'
    'Xm=WJ_<1_Kz^l`*I;QzWT;OZHfoEABO!kEza-I^A6~DRGk`%3>i4<7=pN|g5iran^X@do9`;I5U^d9KrZcA2~?8)8R=+~vCS)rGJ'
    'efX2O&O(R7xTT}#8#R@NLBeZGW2*uMT#LCH(ekBe`I7|v_eBk0;9De^hCi)TyYE`+0Zf6%V)BLVpLKT{#1lKcyWL`Gb~=Z>s@aG`'
    '2_N-t9d1OzTPB0|4g2FANNrkcT*$H|W5bPz(G)N`WdUOF=rbK#tU0jXm5w)G{-7DsFCkyjyZBQHYUB_J20ul8Wo8;vc$nfZkW*=E'
    'w9-YDX((Njp`0CeAOy>e$8!$vL*B2_*1(rkk7X;;9w34}Cw$P1P1U)SvOeW4yc%zP`YTJ~#qc#K2$V^^$>vTAW(YeXQ8}k)V%xWy'
    'u_~3b*#(J1WpLYqqr>A{CflB#hZi4{$ua094cIl<I_B4T|F!uuebvS);S_kd)%O2qUyEX1`|H$9L^ZC+XXdE&_PgWU2%%W6YM7ce'
    'Ep<47lzNl?8^pY}@dl-U#;weC|3i65l=C=09X-PT`cF-M`a4;6EAZ?2Hy=M~r5iNUPt}t>@aV(NNVWgu)<LGK;S52K`rNOCor|DE'
    '1+7;Z7=W(>L77n$btZL)cDY@Mf$s}VH#4XnW_8_??}~^#)MqGN`R=GP5(K*Y=+EITS*X`ud_V(ro3_!yzBSx8cLTeQW@t7!b5hqX'
    'f7b4TXQ?Y1);_8bLOZ<<`sCJ^_pwG#65OT}&JKvSIZf_i_g4rt*EWCc|F7jj9ReXda)<K*!he~Y5U(XafyUH0PS<-kg8XZ<%DL^^'
    'Y(-$*IacTF8lVuRPlZWTY>@%f_cmP1=;HWL1LQMIznc+QRfCo<DKOcsJvda+pnr0!<@<ga3UY>`ovkr$O}(Ove`f@qX8=a&`oHMP'
    'pa|B=ZGM=Z9>TG#et-p~mcZkevOqNLmw%_kp~`qZn0wfJk9Eb=kv>QK>mNf0WZFbwOUp+gP2+ToI>y}WQ&|9`^D6ko1uh53y*fm7'
    '8b?d|Klu4X^b`!F+l3Kboss_53?ldy*(D-Yx5quTHBth#C)E`;jPv8_pkG@dlskFpc{N_@@G4Ov#y08_q0}~$s=>L0oB3wwRTYX~'
    'q2AX-^qc{A5D{4fy=&DV6c|w5^f|RgJ*Lm;iw?i19O;f3YtyMwJVUxtPqQ_VkTamdN2@$-gQQRQ1*%qPnA`7eH%`;)y*PMnOYc$*'
)
PREFIX_RAW_SHA256 = 'a8723848b495886698581acaf00aaad8770b5c7884d6e5cbf82e0ba232c0914e'


PREFIX_DEPTH = 22
DEPTH = 45

def transfer():
    return Transfer(225, 225, (117, 270, 4), (126, 286, 4), (144, 356, 4), 666, 824, -14000)

@lru_cache(None)
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
    width = len(conv.main)

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
    if energy['THRESHOLD'] != gain['T'] or gain['T'] != F(-14000, 40000):
        raise ValueError('threshold agreement')
    rows, initial_energy = moments['compute']()
    cost = energy['energy'](rows, initial_energy, moments, gain)
    initial = energy['initial'](rows, moments, gain)
    if F(cost['energy_upper_1e12'], 10**12) >= F(15415990057, 125000000000):
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
    if F(entropy['coefficient']) != F(7836444120775813, 1541599005700000) or entropy['checks'] != 72:
        raise ValueError('joint score or dual accounting')
    return {
        'threshold': '-7/20', 'graph_cap': '15415990057/125000000000', 'coefficient': entropy['coefficient'],
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
