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
        precision = 1 << 100
        scale_ratio = precision // S
        def exp_interval(x):
            lower = 1 - x + x*x/2 - x**3/6
            upper = lower + x**4/24
            return lower.numerator*precision//lower.denominator, up(upper.numerator*precision,upper.denominator)
        def scaled_sqrt(x):
            return isqrt(precision*precision*x.numerator//x.denominator)
        exp_lo, exp_hi = exp_interval(F(1,4*H*H*q))
        base_lo, base_hi = exp_interval(F(2,H*H*q))
        lo = scaled_sqrt(1/(q*PI[1]))*exp_lo//precision
        hi = up((1+scaled_sqrt(1/(q*PI[0])))*exp_hi,precision)
        ratio_lo, ratio_hi = base_lo, base_hi
        c=[(S//2,S//2)]; density=[]
        for _ in range(8*H):
            f, g = lo//scale_ratio, up(hi,scale_ratio)
            h, i = c[-1]
            c.append((h+f//H,i+up(g,H)))
            density.append((f,g))
            lo, hi = lo*ratio_lo//precision, up(hi*ratio_hi,precision)
            ratio_lo, ratio_hi = ratio_lo*base_lo//precision, up(ratio_hi*base_hi,precision)
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
from operator import itemgetter,mul,add as integer_add
from bisect import bisect_left
from itertools import chain
import sys
# weights and probabilities are distinct integers.
WSCALE=14000;Q=15498
def need(t,name):
 if not t:raise ArithmeticError(name)
def ceil(n,d):return -((-n)//d)
ARITHMETIC_BACKENDS=set()
def packed_arithmetic(size):
 from ctypes import CDLL, addressof, c_char_p, c_long, c_void_p, create_string_buffer, memset
 layout=Struct('<'+'Q'*size)
 def pack(values):
  return int.from_bytes(layout.pack(*values),'little')
 width=16*size
 def python_product(a,b):
  return (a*b).to_bytes(width,'little')
 def python_merge(a,b,c,d,e,f):
  common=a*b
  return ((c*d-common).to_bytes(width,'little'),(e*f-common).to_bytes(width,'little'))
 try:
  gmp=CDLL('libgmp.so.10')
  mul=getattr(gmp,'__gmpn_mul')
  mul.argtypes=(c_void_p,c_char_p,c_long,c_char_p,c_long)
  mul.restype=c_long
  sub=getattr(gmp,'__gmpn_sub_n')
  sub.argtypes=(c_void_p,c_void_p,c_void_p,c_long)
  sub.restype=c_long
 except (OSError,AttributeError):
  ARITHMETIC_BACKENDS.add('python-integer')
  return pack,python_product,python_merge
 ARITHMETIC_BACKENDS.add('gmp-limbs')
 first=create_string_buffer(width)
 second=create_string_buffer(width)
 def compact(value):
  length=(value.bit_length()+63)//64
  raw=value.to_bytes(8*length,'little')
  low=(len(raw)-len(raw.lstrip(b'\0')))//8
  return low,length-low,raw[8*low:]
 def one(output,left,right):
  memset(output,0,width)
  if not left or not right:return
  lo,n,a=compact(left);hi,m,b=compact(right)
  if n<m:n,m,a,b=m,n,b,a
  mul(addressof(output)+8*(lo+hi),a,n,b,m)
 def product(a,b):
  one(first,a,b)
  return first.raw
 def merge(a,b,c,d,e,f):
  # Two output products minus the same common product.
  one(second,a,b)
  one(first,c,d)
  if sub(first,first,second,2*size):
   raise ArithmeticError('packed product subtraction underflow')
  left=first.raw
  one(first,e,f)
  if sub(first,first,second,2*size):
   raise ArithmeticError('packed product subtraction underflow')
  return left,first.raw
 a=pack([0,1]+[2]*(size-2))
 b=pack([0,0,3]+[0]*(size-3))
 c=pack([1]+[0]*(size-1))
 if (product(a,b)!=python_product(a,b) or
     merge(a,b,a+c,b+c,a+c,b+2*c)!=python_merge(a,b,a+c,b+c,a+c,b+2*c) or
     merge(a,b,a,b,a,b)!=python_merge(a,b,a,b,a,b) or
     merge(0,b,a,b,c,c)!=python_merge(0,b,a,b,c,c)):
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
  self.core=55;self.middle=90
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
  self.fixed_normalization={}
  K=(WSCALE*WSCALE<<40)*Q
  for setup in (self.A,self.M):
   norms,totals=setup[6:8]
   self.fixed_normalization[id(setup)]=([(r<<96)//K for r in norms],[v//r for v,r in zip(totals,norms)])
 def row(self,setup,index,fneg,fpos,gn=None,gp=None,plus=False,Q=Q):
  isaux,rows,pack,prod,merge,full,norm,total,parts,cnt,fw,gw=setup
  norm,total=self.fixed_normalization[id(setup)]
  row=rows[index];ff=row if isaux else row[0]
  fm=pack(map(mul,fw,ff(fneg)));fp=pack(map(mul,fw,ff(fpos)))
  KK=(WSCALE*WSCALE<<40)*Q;both=[]
  if isaux:
   ac=prod(fm,2*(full-fp)-fm);bc=prod(fp,2*(full-fm)-fp)
   if plus:
    both=[Q-(((s-y)*r)>>96) for y,r,s in zip(coefficients(prod(fp,fp),parts),norm,total)]
  else:
   _,active_pack,on,off,corrections=row
   gm=active_pack(gn);gp0=active_pack(gp);Z=on-gm-gp0
   ac,bc=merge(full-fm-fp,Z,full-fp,on+off-gp0,full-fm,on-gm)
  outm=[(x*r)>>96 for x,r in zip(coefficients(ac,parts),norm)]
  outp=[Q-(((s-y)*r)>>96) for y,r,s in zip(coefficients(bc,parts),norm,total)]
  if not isaux:
   for at,fi,gi,mass in corrections:
    fm,fp=fneg[fi],fpos[fi];gm,gp0=gn[gi],gp[gi];Z=Q-gm-gp0
    dx=(Q-fp)*(Q-gm)-fm*Z;dy=fp*Z+(Q-fm)*gp0
    outm[at]+=dx*mass//((1<<26)*Q);outp[at]-=dy*mass//((1<<26)*Q)
  need(min(outm)>=0 and min(outp)>=0 and max(map(integer_add,outm,outp))<=Q,'tail simplex')
  return outm,outp,both
 def apply(self,setup,m,p,n,u,plus=False):
  fneg=m+[0];fpos=p+[Q];gn=n+[0];gp=u+[Q]
  outm=[];outp=[];both=[]
  for i in range(len(self.y_axis)):
   mm,pp,bb=self.row(setup,i,fneg,fpos,gn,gp,plus)
   outm.extend(mm);outp.extend(pp);both.extend(bb)
  return (outm,outp,both) if plus else (outm,outp)
 def entrance(self):
  # Analytic first update from the finite all-positive boundary.
  h=Gaussian(1,26);D=self.d;Y=self.yd
  m=[0]*self.msize
  p=[ceil(Q*h(57*(40000*(41*D*j+13*Y*i)-self.T*41*D*Y),13*40000*41*D*Y)[1],1<<26) for j in self.y_axis for i in self.main]
  return m,p
 def serial(self,m,p,sweeps):
  # Every main-row update increments the common certified physical depth.
  n,u=self.apply(self.A,m,p,m,p)
  m=m+[0];p=p+[Q];n=n+[0];u=u+[Q]
  mw=len(self.main);aw=len(self.aux);updates=0
  def meet(dm,dp,start,mm,pp):
   end=start+len(mm)
   dm[start:end]=[a if a>b else b for a,b in zip(dm[start:end],mm)]
   dp[start:end]=[a if a<b else b for a,b in zip(dp[start:end],pp)]
  for step in range(sweeps):
   for i,j in enumerate(self.y_axis):
    y=abs(j)
    period=1 if 100*y<=self.yd*self.core else 2 if 100*y<=self.yd*self.middle else 4
    if (step+y)%period:continue
    mm,pp,_=self.row(self.A,i,m,p)
    meet(n,u,i*aw,mm,pp)
    mm,pp,_=self.row(self.M,i,m,p,n,u)
    meet(m,p,i*mw,mm,pp);updates+=1
  return m[:-1],p[:-1],updates
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

    A = F(140999, 500000)
    eta1 = F(389725633563, 1000000000000)
    eta2 = F(142082319651, 500000000000)
    beta = F(19132256499, 500000000000)
    theta = F(6836467013, 500000000000)
    lam = F(0, 1)
    xi = F(-6843012501, 1000000000000)
    kappa = F(35817933501, 500000000000)
    margin = F(1, 1000000)
    zeta = F(1021877041, 15625000000)
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

    def exact_bound(const,terms):
        bound=const
        for c,arg in terms:
            lo,hi=log2_bounds(arg)
            bound+=c*(lo if c>=0 else hi)
        return bound
    def exact_test(name,const,terms=()):
        bound=exact_bound(const,terms)
        if bound<=margin:raise ValueError(f'{name}: insufficient exact supporting margin {float(bound):.12f}')
        checked.append((name,bound))

    # Categories of (queried bit, product) are 11, a single zero corner,
    # and the other two zero corners. Entropies are in bits throughout.
    p_anchors = [(211531, 211531, 211531, 365407), (218122, 218122, 218122, 345634), (220849, 220849, 220849, 337453), (173747, 257648, 173747, 394858), (178375, 266322, 178375, 376928), (177504, 274272, 177504, 370720), (158972, 207369, 207369, 426290), (167179, 214487, 214487, 403847), (167179, 214487, 214487, 403847)]
    for b, ia in product(range(3), range(3)):
        alpha=(F(0),F(1,2),F(1))[ia];o=int(alpha>0)
        cu, cv = (eta1 if b==1 else eta2) if b else F(1), eta2 if b==2 else F(1)
        cz = 1-(F(0),eta2,eta1)[ia]
        du, dv = ((F(1, 2), F(1, 2)), (F(0), F(1)), (eta2, eta2))[b]
        k = F(1)
        flags = xi*(b-o)+zeta*(int(b>0)-alpha)
        p = [F(n, 10**6) for n in p_anchors[3*b+ia]]
        if sum(p) != 1 or min(p) <= 0:
            raise ValueError('invalid four-corner anchor')
        U = [p[0]+p[2], p[1]+p[3]]
        V = [p[0]+p[1], p[2]+p[3]]
        X = [p[0]+p[3], p[1]+p[2]]
        Z = [1-p[3], p[3]]
        for index in ([] if b==2 else range(4)):
            u, v = index % 2, index // 2
            cost = (1+u*v, u, 0)[b]
            terms = [(cz,Z[u*v]), (cu-du,U[u]), (cv-dv,V[v]),
                     (du,p[index]/V[v]), (dv,p[index]/U[u])]
            query_weight = lam+kappa*(1-alpha)
            if b == 0:
                terms.extend((query_weight/3, r) for r in (U[u], V[v], X[u^v]))
                if o:
                    for a in range(3):
                        event = p[3] if index == 3 else p[a] if index == a else 1-p[3]-p[a]
                        terms.append((kappa*alpha/3, event))
            elif b == 1:
                terms.append((lam,U[u]))
            else:
                terms.extend(((lam/2,U[u]),(lam/2,V[v])))
            if b and o:
                for flag, marginal in ((u,U), (v,V))[:b]:
                    probability = F(1) if not flag else Z[1]/marginal[1] if u*v else 1-Z[1]/marginal[1]
                    terms.append((kappa*alpha/b, probability))
            exact_test(f'rank0-{b}-{ia}-corner{index}',cz+cu+cv-theta*cost-flags,terms)
        for name, credit, hmax in (
                ('input-U',cu*A,1-A), ('input-V',cv*A,1-A),
                ('output',cz*A,1), ('info',k*A,1)):
            constant = credit - theta*(2-b) - flags - lam*hmax
            if not b:
                exact_test(f'rank0-{b}-{ia}-clip-{name}',constant-kappa*(1-alpha)*hmax,[(-kappa*alpha,F(3))])
            else:
                exact_test(f'rank0-{b}-{ia}-clip-{name}',constant-kappa*alpha)


    # Full mutual-information credit when both inputs are flagged.
    # Input-unclipped laws have both marginals in [17/100,83/100].
    edge = F(17, 100)
    exact_test('rank0-2-0-outside',1-A,[
        (edge,edge),(1-edge,1-edge)])
    def vertices(a,b,c,d):
        uv={(u,v) for u in (a,b) for v in (c,d)}
        for u in (a,b):
            for v in (u,1-u):
                if c<=v<=d:uv.add((u,v))
        for v in (c,d):
            for u in (v,1-v):
                if a<=u<=b:uv.add((u,v))
        return [(u,v,z) for u,v in sorted(uv)
                for z in sorted({max(F(0),u+v-1),min(u,v)})]
    def cell(a,b,c,d,depth,alpha):
        o=int(alpha>0)
        cz=1-(eta1 if alpha==1 else eta2 if alpha else 0); weight=1-eta2-lam/2+kappa*alpha/2
        if weight<0:raise ValueError('entropy chord sign')
        u,v=(a+b)/2,(c+d)/2
        lo,hi=max(F(0),u+v-1),min(u,v)
        for _ in range(16):
            z=(lo+hi)/2
            if exact_bound(0,[(1,(u-z)*(v-z)/((1-u-v+z)*z)),
                    (cz,(1-z)/z),(kappa*alpha/2,(u-z)*(v-z)/(z*z))])>0:lo=z
            else:hi=z
        z=(lo+hi)/2
        ref=(1-u-v+z,u-z,v-z,z)
        if min(ref)<=0:raise ValueError('cell reference')
        checks=[]
        for U,V,Z in vertices(a,b,c,d):
            p=(1-U-V+Z,U-Z,V-Z,Z)
            terms=list(zip(p,ref))+[(cz*Z,z),(cz*(1-Z),1-z)]
            if o:
                for x,X in ((u,U),(v,V)):
                    terms.extend((kappa*alpha*q/2,r) for q,r in
                        zip((1-X,X-Z,Z),(1-x,x-z,z)))
            for x,y,t in ((a,b,U),(c,d,V)):
                for q,w in ((x,(y-t)/(y-x)),(y,(t-x)/(y-x))):
                    terms.extend(((-weight*w*q,q),(-weight*w*(1-q),1-q)))
            checks.append(exact_bound(cz+2*eta2-xi*(2-o)-zeta*(1-alpha),terms))
        if min(checks)>margin:
            checked.extend((f'rank0-2-{alpha}-cell',x) for x in checks)
            return
        if depth==10:raise ValueError('unresolved entropy cell')
        for aa,bb in ((a,u),(u,b)):
            for cc,dd in ((c,v),(v,d)):
                if bb<=cc or (aa,bb)==(cc,dd):cell(aa,bb,cc,dd,depth+1,alpha)
                elif aa>=dd:continue
                else:raise ValueError('dyadic symmetry coverage')
    for alpha in (F(0),F(1,2),F(1)):cell(edge,1-edge,edge,1-edge,0,alpha)

    q_anchors = [(598917, 869577), (579338, 938511), (568684, 887199), (696197, 999999), (671447, 999999), (666871, 999999)]
    for b, ia in product(range(2), range(3)):
        alpha=(F(0),F(1,2),F(1))[ia];o=int(alpha>0)
        cz = 1-(F(0),eta2,eta1)[ia]
        w = eta1 if b else F(1)
        offset = beta+xi*(b-o)+zeta*(int(b>0)-alpha)
        cmax = F(0) if b else F(1,2)
        if b:
            jmax = alpha
        else:
            jmax = 1+alpha/2
        exact_test(f'rank1-{b}-{ia}-clip-output',cz*A-lam-theta*cmax-kappa*jmax-offset)
        for clip, ref in enumerate(q_anchors[3*b+ia]):
            q0 = F(ref,10**6)
            weight = lam + (0 if clip else w)
            constant = cz + (w*A if clip else w)
            old_query=not alpha
            if not b:
                weight+=kappa*(1 if old_query else alpha/2)
                if not old_query:constant-=kappa
            for q in (F(0), F(1)):
                C = q-1 if b else q/2
                terms = [(cz*q/2,q0/2),(cz*(1-q/2),1-q0/2),
                         (weight*q,q0),(weight*(1-q),1-q0)]
                affine = constant-theta*C-offset-(kappa*alpha*q if b or old_query else 0)
                exact_test(f'rank1-{b}-{ia}-branch{clip}-end{q}',affine,terms)
    for ia,alpha in enumerate((F(0),F(1,2),F(1))):
        o=int(alpha>0)
        cz = 1-(F(0),eta2,eta1)[ia]
        offset = 2*beta-3*theta/4+lam+kappa*(1+alpha/2)-xi*o-zeta*alpha
        exact_test(f'rank2-{ia}-clip',cz*A-offset)
        exact_test(f'rank2-{ia}-unclipped',-cz-offset,[(3*cz/4,F(3))])

    saving = beta + lam + kappa + margin
    delta = saving / (3+A)
    cap = F(61533356109,500000000000)
    coefficient = 1 + 1/(2*cap) + delta
    return {'checks':len(checked), 'local_reserve':str(margin),
            'saving':str(saving), 'semantic_increment':str(delta),
            'coefficient':str(coefficient)}


PACKED_PREFIX_FIRST = (
    '{Wp48S^xk9=GL@E0stWa8~^|S5YJf5<HiS+^j!c{9sEF>O4u-D(l`DNH-ty9XSM)D6f<`gH$e?Kk(?tg<bXN!-2RA?EKdi|WS1c!At0u5zsQd_#D7`9b3j+7_@(nbHYGmm(JQz^2DT3aN6dBj{v_na$ElHWDn)9Pma$&mCP94$er!USyzT@8#YnztGUB10b}gg|c@kUHDdvmz<#G21-8Ky*EW*-1C|gd@RFbyV{vuo>>lZa{'
    '@rpQSu7pxc8F_dMk;HPQ1bfQ4!@cYzj{o+PV?DjIaff(bS+iS6I5GvLE@s2+sA8o%*$k~7>pi(r2o|nlItd_o=F814%{}^Mk1o8zbdaV&B8T(6{aUGYn@OI$705MpFdjmCRAMb{hk7rorq%T`hdxD>r6%NW<qN^qo89{LDK%25tH{i+a68_9e4zE<$MV2Pv~4~GL1o0j5WC5Bob=AgO3!v3S+d9Bguo{>@ge@O=UK;KUZ393'
    '2etdP>-=%q5YEWWiGeVy6m!=bsA0Ex)8c=nN&f4uv5tM`g)Ya@sdOZY5+lepLW*_;wfBLnXi4e}CRYKA_A6hI=yfT2(4q|H8t;Obl5_zGzC=_XT~|@<e|>yeYJSBWimpd^?01idvr$Y?UnW7AM~jKmmm`Y#R#1_2wzh0@%FhBeFWFHCFQR^6#@1kDsS;OI!oBon!v_7&4!?C6GZCnC|9;d+)ai{DCH$3S4VUiw<$Er~$(~yo'
    '_t;TW`YqOV(I=9@$<@)+!41pU>xIcUhwsG5=5PC{D%*R-K3<n6fVMGP3sA|TOIi9D<#$N!;-uStXy59RCkIUDBh&*c6uZsnMHk+@cF<mQcNvZk-h46uJM>1X5asv8=*@F<QoKMQ7Wj%#$G2ol@1U^{E&<&(-w@pL+v~^ePdtwNy)D<HWV88YR&la35%5Ph<Q<Y&MNJcz<SP4<S@+%yMfWyt;=<g?0SwDcm6D4t&y>yZwXNh('
    'K!+&7O<2Rra82~H5NZXC-78WI3qB`?*tU);eIscps{J|Zvm3xL?|{e-<^r(q7q{w=U{-ZTaGhB4dmAiLjS78@!SoxYelFs%&g-^=o)FW+HS6~3X<E=+P@aWi2!nMK7#+egu6(Rgd4B$JM+J`!p2rVtSLPb&5~`=6$dJt&7bi2TWv?}D%6X<+*+i!o!*TB_d9em&v6k-7icmcjR`&o@fC+JA<49EA-&Tyu2M<~RE(DckbHyxV'
    'p5}PU6iKtc$&togf8CY@n+BzLRdD|^XFz{Ti~EI<4=}O^yZd}ZzZS#jI;{^e(zG^FCiQk4zCZc$n}EKl9Ap|Qih7_zDMiD7=s3!pY@A02mkFH1_c2peT)R_wa&f|lC@1t9#?eH968|C=M+OS0%t_e1_1LW3{rJQW3%3rZYN;q;(Cy<cCZ^;M#8!f*9ziP9ob3)8RX<vT2r7)gotI-|7PMOa$4%7>+~ZI932LHiS8Pi^IW>+L'
    'q^7C52uEB5EqmB)mJSr=C|kA=BGr}|m_%QRV^>a2ncx~thg}^Q_2JcZj2jv5l(g0=ma=TD$ns71qNVJo6I8k23%Yyx!d9#xQ}Di&eBOzjNgmnB@_!#Trd@ST-nGTT!^Urrn#+XvGp)!4(aE)~<FA?<0aO#j!J=<5U5`|Eo6qQWw}JQTKywm*r>-ii*Xg{qr8iSNP1Yc#6%RX(RwAb(^*NGtLPF`mHVc<?%xtj}7@nl*T<Np_'
    'f=#kNrATTEc)DAXXw~gfOewr-<SNj`O`|*wI;U$J2%3;&n<E1PSQ+2PvN~}ePR~DmZ}EWAE^1!^yC<*gc<_slU${4>Q*Fj+-S&@$+tnS2A&qt%qh3xi7Fqd7x~>@HpJljRPm67~`}tL2iGY5tvq876y5&AfusONI3Vh@cXPf_O2F@ncAld*9pBR?xBJ#}@=o%kd*5?M+vUv7lHoG4OqeH~l4+Q;jWe=%ZO8-NQfrw!3*8B|<'
    '#B!tj=>8dqwru_`4Ul<wSS?o>vN)_-O<HU0^Fx#Ta%{?Qzi@r5IrxNQ)qDUbQ?_NM?BM$|ii<AMp(XK_9jsc0iYl+v_2A-Slmdt&x_Mza8!J?)rAnS5%l!mMnrfMu624X-_@Q`MgCdI3a<jd#p?B<OG%c&$OAg=FyzaXvyLtt247uLJb@8vYjqiYlx7X*=$4M!>*0x9aofSw?L$S3^xMxJ@E>Dw!DT;<Hb^biTdvyH5efOe2'
    'B(ai1n!J8xtZEiQ)q1NjK=rCR>osf_xqp>KlwBb&>0mZ<8=yXbWRY0p(N>***SE!^Y-Uw$0=J{gAcJLq)ydP}BX8P`8ThSu8$=I?u`qvoD0sT<kjLA#E`99v)CPUCw*%=|bZ*&UB+qU=)~sMbIOab^jNXs@A;^F@knQrbTGG`>cSGH_DNs_(aNtFi<2a=f!YuSEy7hWVgzY7_kcu=p$iunkYaw6muSHo4qGoFZ?ZoIt_)E2*'
    'T|ZPWb|FR%8IQ$%-D}`lQ2dKnUbZVa?%kyBFXpt*zo3<<ipKaDl;Oq3aaMv?=6S0`yq_yyo(V7)rC%XFK(#AIB-nCgSA;2Z%j3Rgv8SuZ5}1&+m9~5U&7q>_XG)*tBYA+WRQr*9hw+eKcAbt>0a^S+=tTN81F;2fys5x#*3Ca5DuolX^kP$75^4bh>1VQdbe|Oh^+-+$Q2riDCkWW?c{dn;k*#CuaC|Fj6tY+*Ev(3T$S&k6'
    'j+-ANHE0zYagMB<DP^VEd_Y+E(1Ye_8I#X?y%dz!o$#RQajcLk$lmF2G=JL)y6i%Tjyejj5_DvNUW*_IGCCLwtq+a1T8u_4I^qIM2>w*}G1VKc`}QU}@NahbHr)5<99+$WV-Lkz>j<ExI7dp^Egh)%rA=|neiT5vc;c$vD(~k??b=vexUqur^_Ge32_mA?4V&$0;ZGi7dbTi!x~RB;F6N&49OskJor>Rw<))N$7677q@|E7#'
    '&i`xC2p-Lj-g*U(m(>;v2w`-VURq(8wNz0GI<bfu=Qu@9AlFCxrFcL9R-jM!%GalYp@81SKArJuayIP(l=APEmv{#X@;^xaX$L&*qj>zf5wseKS>aL}JM>E<t4zH(=8aoVcHX7^xyAi?NwNY>tc)CsxJ>iZ-`hNsW_nS;TY$Q^UvZ**icvoMD-ao1`)?5#)~*FfEY$@_8wPe2H6MDX<LF@&(>is~%#02WoF$s=9?B?@^t}1q'
    '|I;dwa6r9fwu*XdH%~$4V_2k1h2@Ie)r(VJuE_dY=N&q_u7ye4F$&p%hRRB(L39@YnJQwT_=#dv8Pjmm{iIYyn53R-w#(Rfx<Iihg6(A8Z@Qw&pBMe*kRF$$#rA%(H+}^6khSoxIk4joQrKbIE|YKY7Nmp;^b75nl0N0XDv`&fP8Gz}F>v-=Yz2#n*LkxkSIFk#ebinzxVi2}jTuw}me4fa0B(bpTVi;pU|e<yM>dC-lDk9t'
    '&9N`&eSITHKAo(KmYVTAm;5PR3$7Yp+QZ=m#B|ugx3|>oDr9BF`4$ydm`*x^x)@1ThOQa`5{-#%46kcaWP)hS5eKtgWt=FEd1Ueh?v$^cJX8e<S*lw(fDH;M745Unbe54uqv5LWym%?e#oR6bHLnUh&G^RA4Ab%x`0UIGEq@!cc~#wQFao6}NN3%S=7;f`PxhTrMRi^ip)^`ZBPNu58*p$M2rL$VS6_bXLH2i`Z=w7kg?4}z'
    '9neSW3L^|A+1Z{q2VI8qkc+|?_6&XV$1(h~3=N}D4asvgEI*`QK}y#5fP8`gRj`O=Du=;1&mb!TJdX{_th0~_AYXRnN8X-TPJl~VoqF#P(DWPe4Q&!Kytpvx@T<+{Umk4JO2)2HNS3!-vhM<@xn(>_bE+(`XTD6PVM+U^L&5_&gih{i*szyO>H2!%!eDAfV;Iy<QB}dbNi<ek`IDIk@i$;79=>`146Vk)I)=N&T5?MTDoEty'
    'yhr7_0J>_nCMQpzgxS78T5hb*GSFBr)qJ?<Y)JD2{NG%yFH~D7u&b;_=JAXH?;Ow%3>_S$5geo4Kty8NBkpcLfLLTLJ8W;5@4)V>*dhu+RoW?;-%;au^k1z$TV@wgE2Yo#%DM<S^}I+F^~`e|^s4${IKLyCrE{jKVqb2LqPaIlI+cVxA}kaNEaL-0ui!iIY46@y^>>$QcnUPJ)X?lNY<BWD1IW!g#bq0$3n#(BM9RFcvzsBg'
    '0nO)BbU7*9zA;>&e7YW+YQu2rv!^C)DBq7dUA0?ouazs*_c$_yjEnN!Z{$X8p+%g)=z*Pj2yKHFhSt_AbRgvbN8yAx$vNLpfil9`L)_0`M%7t7v#sSAqaA$+I%W@QOh}?+G{DjVRTsk`4{}aI$fx)QP{rGl*eJB~qRJ8$c~rGaZpZTQo5XmrO{?fLSO((1(pA2Hoh4c(F<o$O08su{(C()3|8&csK$YpkA|J{r-<^|M^QPc1'
    '{urL+{3u|DIv1*up_*{m11oJ=#j`IaX%|&dNqJeaKdzRlHDG|{-%Mg{Z^Ip7XC$^yQzgUa?<{cvN3|MDDHmO-qB?j*4CDHaN>A;uPZy%2h`@TR#*%tlx(XY=$8Srzy@;CwXgQ4ttWa_1lxa)m^zG)Is<gof0Y;m%)voCu&r0uh$@d$Lse7S9r>&{<`N~;zu2iScTcUs^aL1-|C18ZPk3jTMd{DmEnk$|S`q*6TTv?CKpcIV}'
    'xDaD46=!?PjB}v;=RQ;R7hrc(-dZ&VvIY8I!OkSX5+OJ&Kr1)4pLG7$58Z1;Xx^Y&h&y-dKES>;^q#htu>X+D8toL=vGPeU3b5FNiROMq2B(5M*|Z3>_a6h``KgFcLWO{cal9N!iFUxZA*14wb(Tv;ra}O08YNZ3g-YF@0I!W`yf83z11hZ6`?m!!PhMX*k9>hhJVHYj(a16rDCKXHhk;1r%fsYcOY|fKA@0@-KnM}zN4={1'
    'ikX|bwXSZQ8N)*ZthHeAO}_%5x1nMU<01mTGg|H<s4#Qt6gNXd8j!2giZnB%IvNwn`O!32I{-W~%IvI6@W%7x(;+E!{VXr0nZH-T%9Rv;*OT6Fi<3B{goG+{aG!T|<#b^P9o8i$Zex+9P-uZs(n9qD7>F$(R)=Z~r4Z4igt{@iAnLBX{TE#xAJn8l6JDgc(1M<%ZCYhMR*(xg{<OID`0Ot_%XUa!U@V=^co)oQ5KQ>~Q>YD9'
    'i%k^t&A|h3<ltOC>K6(M3U7(grGV6ct9a73hs^q6sv(>dk{0U*(mDR)`Tmuej~{I_VfiPt35%o*sZ(PvCvVChH}%3-mG2>p2$NhBpozP*@CgDxONjk{W|kWhs%F`If%N-H(L6cFX*{&(r#}BMKK#+IOY)INStN3}(t7d7=_w&28NuvuyQyAO^jCu|?r{?;gMm`0H`)thEy!uCJ8y{Ni7d|pWK$x&_(sn@&1N2LaxxRX^-kM1'
    '_Q6Tsuj7A5J=(X)y}7Z`WWte2w~h|&G$IYfE98(s39JEc`GzG;8yCy?)ovQP{0OT%tfzy3jbSFH8~}c=2~q#0$`|elj+Mjln_3JSC?CwPxG-<ZbnMVuA$I>*qJBN>DotLpaYEr<H1nUf#_Hb`9$s0o1}hMmbFOc`hejA%0$(qd@`+9qM=Eo7X7MU9=|-)L6)GBP=oWsYkZNImeIBf>jL1XFoXHM~!^D+h5II}iH%&{CmG7EX'
    'A(;7`^aa4#Mkz?onAs%H8<Q+MzH698>V4~?BL&G=k)|Mp@QipQ{IWIQ0kt<Fg!yB+Ra)0je^-ld5Y8gCtfr|^!fv+nvc-d34SF%m*7O5j$A#O~^PPOzDu20@8H=)0BeI)v<b4hsXRJJl`6W~#bMMRv{VV*Ys!0@!Y|CIJniMGFgbEt$FfSh}!meRQh&=*znHIRoi<AVm=`-y?H|lt#oJx&=Tk<Q=r1{9rJDo>-2Zu>>JZ6)X'
    'xC;+Iza7vH%`&3QMSK}Y0ltE=Pd-e6WT>BkN$YAXZiZ^qEVpnCayFcQr)`*@pO(lVc|xbzZ%uWjr~h-bO#8jWX`Ok4h@SUjTRL-`blXm%K3(D9fJ+cQUH^gfcg4?!wF%CKgXPZES})W9MR(GWG@g<#-uLz42z<~G+!mZp+i_^B9y)}fXIbseEfxlsj!?y>#wCfBrfJr>h-g2cgpt3Ke3^C0*yUmEi5Xd`W{N&Xg5nIJHuq)Z'
    '*zsWA;cqQwn==-Q=F&#+(o6QE0*5tn_rc!qIG89Sh{PVi0=DMOP%Va#jJrd-Gu_@GZXcyYgix0PT?jX8(}6LAPv(WHj*`2FH)!o?4+%Qx$~gTp6#^3Kc*Voh#hkpE9}z~{>}Xd3aHpWXVX*?Z?{h@!Iio?IJ1G(i8$G}smcdkh1zQL&o8Q#E!Kb}LlkfhM5}rMmjtYl0F)4$9M%0PsbJR{<Ki!BBEU>-{mIi$=2B;kFJF^nG'
    '2CVj^kauUPLrR@AYOqNWr!7;w%~z<WwuOSQZOC9NGf8(2I4)Qx%_`~eJm)L-?g>Vua-;f@wD#JKi0TS%!^)0xb9A9|I)g3$UdNLHHeQ{tmsf1k7J`ke(5@XLFMbOq&vEICJ&q&&NCAI<+k_{ry<7)}2R<t4;?I6fGm1EJ9uCoa(<+^xfU5CBD`O>FSoof9h2hR`B}v`oXYOk^C&I!2HAvLYTGQ<s*#59Z9|h*Njkl*opB*Ag'
    'xeDUk8Sl-UY<!8oP6j{h*z%YWP?Mknjnaa70VS)Oc?<+yEf54nT`K@SXw}?;E*%LLrgiEq=jX|UEFm{A_%ag{r7HF*Z}I|uR}1Ze(r`x&Oh3>pqZfsS`c`S`7mGjmzSEMT?uJR4Ov+P1NRo?OK$GEk|B<mnRv#25=OMQPN04f?u8Tk&Sr`s3_m+b3h8VCign%`#=ob^we#}!=@Aa%dbW<JJnK7KXVoC}jC2|^C3T(-!pVQ@c'
    '^z$=-wuo^ull!)iu2^DKSjJM{6DGIrB%6}PUUe6007SI1(o|iLlO25P{&vOLNRydn<jmc?$vvOM#Uj|L^X0zc%uHKjELN6Y0z-oFu>nUl==32g;X_>TIeu|`rpu1<WY`m{3FP5Mcrsu-WKsF*@m;PF8|;6O<dP!bR4;9ho!xDIUH@#BcMcx`z!$;0(zC}%!VOwzt9RQF%5h{~NA^UTkRCA7ZW;r4ao0|2*YB0nLzixR7*qw&'
    '9{s<zm9hXp<4l99n^S-Md;o8|cWzw6`5<)D`7B21$ktNx*4^6bI2#Y}G4c%vzn*s|+oou%H-nGRPQJn6S6itg8!4f6Fl$<loC;i_a0rsz`IF%z!kA<#(8aCo?v9zMiAss>`M1&o3E^9^P8&a4K`uyhenr|5fG!dB|6vMI5p_}}&_Ve8bYeWIx%xRVZkwGfh$w}AxKpPDrOD@l0AG|KZ*ae@^B1Qc6BzWY6&D|n`BO`}T&Qa2'
    'T<JO>ZDT&wZZSzM7kPf%ylKudz4q?yL5yup^mN?(Y@0x|j@ds-&Qn|wbzo`a>uzb;V>-;Bq;uYh!p14dia+R-=f1iXTwW<}BxI!>*Ns<9Y90c`z5Cca6rN|{V@l!e$67&4>Ki^G+T@r{Pc!w$T`Q<c$_o31`PVmIKfAzRytyO-*i=MBGkC|^Ua%7v7C>#S<(&=UMemjo7cxuZi|q(J1;vmu%Tq4k4o6?f8v)?_G2Zzr7IC6l'
    'FK#WXKH8wBU9fXW%%N3IsVemD6CuaCieO}%lEj+NmiYF)h-&%wJZ{6H!!xR-VV1hMRK{7u7a0%OCGGeFSNaS@C+0QBcs!0-_6qyFq%Q}6dlx$(e36&TP`Q9K;ZoU_8Pw1e(2zVCM<$+!O1PWLPIiOnD@<o>Eyu6{oxS=Rg+Gs%C|aNTL=MVY-t(d;0M)9?*l`98G)hLKzwM9W3zRib{$EBS^F|mX1`IRaoe>r-6*lhxRyP_F'
    'Z9Zg(x86kpQ78gC1>c8`N=>WYPg_m1a$weO=QJ@rzsw^7-LMO0Fv;><GznT|-is+m!3C7aa-hO#L3y4s!N%6q3v@tE3Q7eP7>Er3h%T&ia_BiqA>UqKq{dAd>WP>#c<g1o)P1nZ#(r9>E_uoY0i1S`2<6ZKj1bHrb`=7=l%;j6eLJUy9Pyt?LDqdgoIp4Sl{{;NbxT}c!3rqGuS_IZ(ZzQnE&&+~%bo(J?(9#C0Jr+z7mLI{'
    'W3kxuHjDs+b_xr@k&Nv;A^f(C>+1O;nKhwu$w|%l1je}S$~E0h$Q2{Oxv?J<P=<Vl%J&#Kk-?oI-=*3pjZAp096z*WBVT_ULjR!W&Cjj9&^AEiWfm>Ulsc&eDS*@*lf>h;r>XEO1AAwpKf{;%9c~j>TgAqE?hk}DluldSZqKuHfAT5hTJnT0gW~UEEh1zRZt*McOySdn3;#?QsKO6_il^&p#Aj4@YKS*4YBn(kS)b>Ku=7h9'
    'pu5!NZTl>~at%6;K<4#jE2fQg??q1RbJ?Oyv+V|<bB5sG*KfUmP!BOh&P=Z0gEfUO{Zif3=uw5`E#Z?onq!Kcj0sJy$WZwBfU(bhqNhtVb5S}gFfmKRtJHz14cZ=!oxxD`310`kl)xC}_!r&93asW>tJI2h0;yAl1}JxK%P*3!c9a+$Jcz16B=sVF=F1Y~xb&xKRg|Uj*B{<GZpKxV`(^X>Y?rtPh&D4ZAi%=QBmJlhbBLPN'
    'MrPn7?6NSyJZ{Dd8T43^aZP^9cQc?q*zu6q!#ZkexOhESrzHdH0jLIR`l$05A}^+A<}l_he^_gQ2;(6OAU|o_qWZAa7$<e^oKq1=aXz^8#?%;P+xhmicD%HiV2KTZ{m&!q)JmUy&g<Z#=AuVqv6!J9m*_9$+e6-+R$h^%x%8H}<rjCE@dXY+ZFKu?iP0%~CtPKg{+<0d-p0qt{C`*FzUgv*Pxt(uz?NM!g*r#mdch(O3v$w|'
    'VQNO?X~QE0@XzXmkQv$ny=(XY6MS2_9$y;k1_u}SLK+HrAZFK^QpEW<8~5{r8!&~60)sg@va(0G1E69I4HlD&gI52kTP+${5|%0-HHNwgO%#LT_<Fu+PAm97%@G<BeL8C+BKq@u02rFOZU8IN^>y@MI#H=#n353PgkTECNC<CMVtfL2h8gnDjz{clWc8*TA9FY{We?^BZ-yA!!Wqr`97ARTia%8T0Ew(~*S)W0yy!1VdZI-Y'
    '!SQ$=dEPp<!CJKl14DAb2AaQzo4T3IS0jsSs@f|JC!PFJLs)t%a*>5@<eQ<_C3AZth+}205q2h(#3IE{Ib%?>NzVDhYN>rABJ>FdtKi1AYb5kPNm5K1(nr69am3gO?=fdda4yA_zE)|k>#<oq!rRuadR2Mn{>}N9H)l^EpskCE8N_ItxN`c)&$`AzPN?lA>2ALVF0Uj)iSCEIZVZAPto7HhpTENQU3J!hhlBrkW_cL6;7v#V'
    'm^+m+vsE7}w>eIBD^CE)zSL4Yf=2T_YyXSiFCP@Oi@wZ?@bJode_s+Ip!OUmEwBmFGq#lx#6nt3$(V1rD=Raou<n-GAn=kz6F}nX)*C1oxiqdxt&&3S?T2>A*V(A}09RV3Y2-kk*XeA7H||V7#0z2R3LoSCFF4^rX04&Ojcbp4Hbi|Oo!vk5V;7on>v){J)k9BNJ_#pEK6<X=liDgYBzsNOyqMQLG<utppJENsle2{4s?XvH'
    'CMTv3?i`E=<Vcw8CX_b{RH?R(t1r*vPdC;fr-{b5jflcPTey%mmg7XOoFbk>sE_jLt4MDe6<))yAg&@o?5$~_;kuHND$&J5rYn3Fw{A$vznf=T2fhR5EnH~F&FmHBkk-_bpVr+rMcnkwwojnFa&H~aZVAU=T>yHkB(HBPx(j9Pd0*Y86l4>Xw_4TUPZdTi7q=hmk8u@eA21!{J3!HEb)G9=g?Qm59oQ)_aq6mgRZKA99EECS'
    'SqVc1Z<hoBF_O?U#|zt?7qlJ6j7cT0n;>foTyc46^i+P2x8MV2#2C+2)JFUoRuFal?A}Vqv6qO|Wp)Y<)scTdEqhMolYOUm+moqE65^uJC7jt+jjM3HnBX<4tg$WpdE6ki@JZ40Ze&jB<I2*Tml($rRGCONgnX|~w;5YOOrW1@#pPi-4C?hRp)fYhhoK!z-Eeyy^&3)`y93kqSB9N4CQMEaarvL>dfoLIAI~E2Knk%UARR?e'
    'Lc&B5Cwh@2J(bRubqn(TBCg?Xg|Dj9>fW!@E3v%awp$9eRHOO$Yqr*b=qAI0=8qj0&>dta89lC@&3~a7q?xl~hNIi%$18`II)2sO^YBZY9QuF>nJSwQyM>_$x&eLnEPDVj53RZdmD7Nly?O7p!w!s%pU#I`?SFb0Tw=g!{jm88Pj@pR79ER!cAA)tGah-tt(XQ4oYgnPLTt-ileaW=%ISO2mxeg?gNWD*PcE=910eIR8yTQr'
    'V4W)QS|lR0=kfA5viN|Cct&^XSBV|`_f6K$Tf?D<e|v!Pp}tB&MzNw;`o4D`V-Y<fu;a30yQSVO#i61ly+Wh~^7_s>pRjFJ)IG-Ue1*|UkCcUf1wr-_;y;Ib!VaYZ2Y=NjRlKy}CVCWXu##FN52cqk>d5m^y}a)YnPVLk@V#@1K*G!}0Jxn_FOzpr_<7Y{X0BoTY)sUO!RXlDigq0HOEQkaPQ63PLtu7(#5Eh!*@y&QbTanb'
    'rFLiP6KKxxSDZg&Ntww&)i0%+04-?h4aXj!?)O&|ZV?f#g*(S;t29cwvUhThQ$UTEz^<>txdY@$Ofvo9_adwdZ#@Mf1ht3TF$P4+6Qo_5*{W(^H5T4K2;r?)BGri(?>rH008&v(Te$+;KZ0lmKXW<^9P-%xSl2LQIwwI5TghxJ8`1I{70iVh&^JP}8JawR;1IwdrKf&vANMNr)(6?NGl9kZ`a_*qPB8y*&}pKLNmzYObU;qT'
    'i_F?eP-otpPTD}%G*teIax`#!{<DnOh96L)7l%pNs~#@%&MM<Kf6<8_CryHju2RK$X~e5p1%NH=95(b=94@iq7jLscxxNh2V>-t7+Y(yn+K)Y*f7$4QfW}D0K1x?B2_ME-kbou4?NZO-R)V8KLUwxb<@`kh$?<d8<+)awlS|6qq!VER7Q4-LJU2n`9D;ORv$Yr%nVDRbdbS(NXBF(z4&w~tG6v5U6m(7=-A6H1H!bbj{k4iU'
    'uF`p|Z^9uyPsunuT%LXB4or8q<xo`-_}DfBf~ss4(|<JtW;yiLrv)j}jkz^Ci<o!K{oCO0Ze5uJ_u&ff-MDaK70#o~!u8!nm<)lhwR+M5ikgJ~uV%FMI6fsPY0Vd#XgaB!LVFj<+cvA1bkFj6tNPKD#egUyw>c52hPO$(Q>ltq8QI8c;#gx<0E3F!BEnePMd#N;rlz;!qA69PirZ2D&2;*mnGzlY@;*@EIF)_Zpf-+Vd7Ak('
    '`k?jEUd6F~G9X(m8uz95`Li7zOE{c%%7E!;=SK`-fDAjM_>m8bYJ{X5uXmM{0d*YBuqZNbtb2ghE~7bEev1bZJ6995v{<{$F*Dy4;Uti`4Gko<p_Mgv83Wn5CGh2ujl<DEN(KA!=-+SDsNs4c%4LL2QRlQELjdYSW#V3>rp+x>>@R3WgY5s|)SQbeT^V^fVG5oPTHUD)X$#q+L*&;b?qYS%(%xn2f|?Vj8@cbOupki0h~jj$'
    'Ig&l}IfO+!EQX4ipIMZ7AYHN`H?s)UsENe@KN-dcV_wZY=8L{qwh>C9?)0d-?8euWsm6vrj-KV{xK+P!4okDQwgb*bxhLY4sSws73`YpPf=qeta1<(I2nch^-x=w?bFI7i{9^L?DWjKdkRz1c6P87laqlfvLmO4ljv`IMO#S!sfg06T&$=?-bsO}croXq3eDq&<03#sNZw<NU+S*e{2PWWZ>fOy{rE`9~ocov^dK3peUaw*<'
    '!4mfo=N=+*I&I*{CGLP22IT>=J}I|4FoCFDj0i)ZSGZ-=a7yOpg`ebLnj^sf_ycN>_B@+Ipl%uAc(e?ls1%!Q1xQGFbq<cSFltdXmEjz_tn7T6!sKt`yhs0+S>_Qdb?YP9GzlLe;ax_IX!x!D9D(R-$$|el)g<QXKcscp@f{|oj^0v$?yTDQVH}~>&NC>`2?zOA4f#dJzE2YND*wky7XBEtl!!+}UG@wOV+5`<KQ80rjT}^?'
    'e1yQ(T%H{<Xp#F2=ytNoSP7cwk7CUkY1Jx^$FRGWCXf>Nx0IuOJH*@Sq2dF2THS2+{oDX#tlz7#n?_jX{<6JRD|GQ<-CMnCdOf!c;*kLhNgjvk0`9(_6DVG7ji5}}Ew6#q8S7(0IJrN`MQ)Hek!Dec4Q%oKwASH<z3ctKihn9=sqKE9#DWNQpLAcG$9N2a1(+z{o9-I^`24k!gk+4Dtzx#Jg;?VN73#&Awd)x}yi3{?>CGyQ'
    '$OPW<Juon3w7w*a=BzrTl6*x=Bh#m$P;-4WM%`plV*{B4e2HR7w5_j}JxjP$;#NEo*4XH~`=eDa<HMR~ObSptLG~<QJddwv%?6fG^v=8OT`;gk54!Ex%{h4uqrdl>F_x3I#fKQWqKu9&niIpqruS%<p}ntaFGn?X<^J4>gEozj(lc+K%mi~bnRoM&`<L?f#p~rOg)H2SjL0RJs=P*KkB$Ua_LUtD;^Rb`p=k8)t&*MTid6u`'
    '2!t{O%!~p_D0aIU{=An#Ui?=w0BHvukKwX-`D#u4)HqBi2Y@G`tn#$(&1kG$BHcE-3a3L7)iYX>)^o906v~m}XQ(vLDnrXapF>JT^pj(vhtd_<tNL<@Kzot4P^kUo5vheB-hogO^6K^3X}OEmvoyys1hGDSeQ6S8?<`knK}^bWMR{t1TKEmThzM#Tr-mGbQ2BxvU&|YqL>3`1f5Eo-ZauACB>3LDc;&i8I{>iNFEzR-26pRu'
    'w0g$@sqZIrXzR4;hViveuQfO`3gk>z2Ze*#*UD--&eros?jOa<`$_hTFK()<tyO3wF%|HaET_0|V*E4*{b@<4LG<>z(nBz<z3WZ_EaSa+uxrE+S(D^LDI-X*;&gmO$43Kvq>wkEL?4qQvWgT0VZ*kt4GUsOC??ftev4rBz=q6H;pk0D7UyCZ>(u0Ay~Z&^gaEqZk&4G+UKHsi^1vg}!MOq%b*A*eOck8N&^ZX0_*zIlX(sQ5'
    '9?CNP<r&?mxjz?l?%VJAxMgH-!0E;Sa%Ekg1Vpv-hq~y{g7#f$^~+|W1SMcG#U$Ba8DMn=b8Fxt5awj@F~(5Tg%a}dFjxak8*dKh0&re8)NR(;u(%#0@WEg_N>0JA1S_eES;R>?>R2LO-_I{&4!rVKUp}hz&izQTa~`vQX#f|65UyFV_i8&{KZB*ev_7WP$uHh)p`F^!Q-aKJG){34W_56vIDs-U00^DetkWotn&m?r?odM2'
    '=wHn)JJ4F*J-QK81U)=Gj(xG6hCRYo4RSxSzg4nj9z^imQPbvpE;tiriOp$A3Mk#HW@Xaz$EV>x;pP^YGuB2`Y`HSkrfi=G4tHd=xDVs3EE|YKF~^1Yyd)a*>L57e$paDrZn>YN`|q@lMPaga_~1ROU%AdZNe9@Zg+1`Y;@KhcUe<VtiBeodq^VT?4lZaPL0h%lZao`^{riPN<az1u*FIIN4cnQ28<(y9>gWsKJ>REt6SSUp'
    'JvlN%WEGOD7?tW%pJhwDnT=tzgTL^^#W9C$D7Wml9L_J;z8pgVfr#w2*6$e26F3DyR{i;&wqV|k!P0Xo&^S1KSzbC>Thi=K`vEGf&Pj70m1UdfTGZdNDt%cXDzCVT4}#PX21r%XShKx?N>@)hR)yWYHfDfv_$_-OWh}{OiXr5}N1ge|n1x|)F;h-~ZzgF;(wwi$%uiZVxmaR0!|*FfDcBs@OIT457}=E-f`#ZOS%n@yJ0d?m'
    '!o;O6FK@^*96wZ0Ha=q8I>L59Ao(o4l0?)RU;whgCyJsiyn5)jZDT<QxMP_nP4s{Fas1-*o0(MfbXez29Z)vYDt7shY#D|BoOPUvXAbP%4Er}<{#lUq>Ll8mjt<*6V)VHcXQU>pb>T2F6g|d}?jmYzJR44hCU`z!rNTg?m0(<C*?u0C`3s{nxZLIkr~q=P>b;l3#G<r(Y9fTkg{8T}(IZG7)-O$hA>$;@g}+0ia{w1=GnOje'
    'u<=CRsF(PBG*Ylq11GY73sG&v6@1h-n~xcmQk1<*Ho`+#v|<eP`>9x=2<r#v(|qwy>sDa-HP^S)MYR+h3We_E-1vL~0AANh1b57)h7a*npL$;*)8InGt4>=$gBiY0^7PAr(WS=KKNVMG$4fjg8+ddlyx*CDN5`65=EqIw*>J8YYbF%6dA<e+>?U+989?gV?Zp0nlX<^p4+9i3WYV7*+)_>fMI^=0&(I1;vPI5!t5$7~7ByVT'
    'Po!yZRh7?@g<(*~BQdBU8ukDTI=@@9Ri;<}ISo?X3WgH8NjXz4fx?*yl(4uz1e)t&)eB!!jC@Mu{<fki5!EkxawLYiNRW6kk1KO{Yi2>dw+v3MaZQH1$k$A@2t8R~^Zj;3i%N!5dU~W&f1VcUjj}lr*)HHH%L?k4HHNYQ99o>@rG1!9pS)$B3OY`ILsEq{#hI;-9<E>V4B^O3DRFCEXn#KDVs!Ai^qCo3MJlHK{>lrulb5Eq'
    'jz18A@M#=?&<vMAHiYv<|1(%9ral_6$msa|xTO)BIYa3sQcaJdSc3}(t-btVM9b%s_r>#seJ{!_=tvj}Ld1d)W_9<qMs%}q9WF!%wdeLR!*4ADihX>ew-_p2=*y3z*uNnV|Bvcp#akL4r^?;4R2;@_68&cEI(ibpWFU<rW*FHkcb<_wPYiqLvS61G=YJXLBN!RVvIA%S$ev`Ubt~(}!RzW%rL1KSUx<~9Z%lYjiEvFn3z&5t'
    '%O0Aupr`uPeu~3dqGt><4cI}q;^+7W=}0OwtT#_#)vh^_bi`tJY6}s18eR(K5!P`6Q7w$fT7VaLhBC@kc*>_XFUq;)TMw(2*yqw+RQyY2jXZmbT=97mlxwS(+)W$p({`Y!$z_yYCa-q?F3i|4Fs9cYs+?1Iu3_GQS!0!iuzBXR3bIdFolFBV!Fz?YvfVMr_-W<CI`iS_-ji8W(*z}89n(?aud6uBu2}Eynn3YH;lTg@d&DcM'
    'mi;&K1f8US0AEPAD1g|>c(URHTIU6{yGlmeXZ2QrswS{{5S2p?VSZ2)$og#j1kbnv&xkZi5sF40k#uA?5IQ(Q%Y-;U;w?V6ESzSnL2IGcDV6UkwSbjFC{l#@nh|4m6GNC;w2&YMSb~-8*QuN3w*VBMda+JcnqqF$b)Bg1t@SW<jn7aKJ*1B2Fa?o!<b#WtT-~#1doVC=`$F0Q8}?^DcQA*f$qQR&jYxtBd8eU>MRxVd)qv;1'
    's|Lbr5(*R!N3BU1w(tE><@F12R7y@3iM%{4v`$gTZt@YO!o#P`{;wLgEyn!k_@fy9L?8#5GMCR3FFJ4toNU)WYM)TF$O+ZwL=)Hqr{3`7vRUT*1kNmzcYJYJI}NaGVN+n1xoO%Ti=Gxr0k7&L-zOiTSvv&*!m@n&HNbtk-iBUW5!J~;M>TUR$>L&M^`HXu9x9EtC2hKCx2Mq+8sktnn0*2*&I$!XssUg(Xp)OMghfWgHsJB+'
    'XNB7X1O1U~3rhp53^8-I-W%FCo<4FByN_&h*Zb{w_<5OkFg#@|p<~v+uzR=T%DPO(sfz#79?vBK&z-2_mR~z?OP9`3rG2GO&Fl_FT-n|3YVm2BWutw1cA-bySIi8BQhz>W-SD7lOC(JrqKC?(v}5a!ArNGOKns<usjBHq4}X)km6I>m19J>+aIm-Z(jWJGQ3YTW*6T!A5dC1u{@%4{_WHSxq9>rp$yv`H2gw6yUr+ehe5J^P'
    '8W-XiZ`|p6K)#>rmxoc?S|VNZB&GLoxSLYNNbRg>bTHaCHq!BsQ7%$LpzKe_u~(+*?+>DG#PTfa+6<P;<`Gs6WmN~R1KY57!DEEXNat#^-7i;#92jdMa7n-8YdC)i5tpv$`VrEr9qWszLQqSJ-J%KkF51Y==}+^IJP{;|;67oKq;W1BpKqPE4WSDjGm`09DnF7$ll*e}uk;9x(R}q%G&v*b{@J8S)OGgP>>7%rq}u>Q9dn*+'
    '&cw7Y4O$m|TKwaS6BKt%TG_;{5SgfX2twLJauFJyY1~Li=frW}FI6DBD;3W@biqMZP}$hS5@MB2MlV;@a(O=;J5K%xp7`-PO?D9T891iIIIT4nRf=5lKv1r2aS*a19cUXigWU$wRFudtQTl4E(1^)AurkNd3|4o;Ew%cJHxADJT)4S!?&`^caX_E#Q|xVFiRB(c(VeUtg0(@xKJ+GC2C2>GLeHTHB~Z>TW=lj+?;V@G<hCB;'
    ')!<9kDJ5>|)_P9<eaq4SqRfX-*CFkI{WG~-oJ5?FOTqK6^0X_l3n1@#XjR=6uAd2^5qF6M;d57r@8*|R!Whsj;aQSy6~zsLfPC;7+c?jMclUJXcYH=@w-e}>$czmHy&Qb@kZf{$f8HwKxbXd=;{0yt@ET?;jn^!rbva^SHhCFWIfhI10#Vkr`%yEBBj<k&cas*I*sr}`P??h)KQ~x?Mz}k(7KUU%mON9jF(dwoHi_lalrM}<'
    'deNCY728}4D@jVOX4lhBkF#WUFW3BOF+V1LAd#Zmf$_oVmYQ@t)Sv#uF#Xa_YHqoIsrWub40J7)pNFi{Ywxd1yj>yRn7>{L={JvpcZ$*LwKjxUp*q&6HmM`6?*%xmx8l?!Y>uBC%-G}9iR6pjulQ<;Y^Y@Q$&1A(!R%M=t(DK^@XdIJVm3GK8tP{RBvB&$9l0Czx^hUjguT0pD(IyXpT9DFwzNvA#k+88p^Jf8;E`3pr=MAR'
    'VLq8v70eD|cYYnp3SAH5E&Sz*CDosXPpKu~oghALwRTM(mSK{(*X{<E2raHJv6dnGLwtBaZM1D#M_PmQ_fVM>y*0;42~CND0K&PizDJ`#!hBZVBC9=%S#1xOzvx_T;}S3fa>;Nx9|v0XWu(Ycdz@hv+;pgm@61R3KPV=RoG2|#<GR+D6uFS1HP9SMAPJ7^%H1_S=>pTw)rnkqfK_O#nd<x6h)zN+XEa}ADe%LG0_?@CGHO~r'
    'W@`1h$wCf}!m`u552SdC(6anrRTEI)lqFV;^KOT;WXpi&R7LPYrsIYAZF7mpB_0Ab6*!o6jdf%<?AW2PM<Zh6Zdc6~69k-FVJX4TX7qooTiNJqyTKwwYjAs;^9LN_hJ1HcQwOQZRWr<l<+?gZDjZXMpSA=bwmf}8>zno9LFARR5G+7VxgLc%FZir%g*3)ALPB3!R#qDg$&|EX9$|hn?ZF{nQ+Jf}U`#P4l0fOZL%?V*N+X+}'
    'M8mML#NNnKUqy@K1X1z|YU02T&`zI8X~m3?ZhQ37p39n-GUT=dqMAUw{dGcvDH&lgcdb(ajBuYwag-Bln(!-II(1M}-_QOuGtcj0b3u!ue<w%!+1IhAFZEZVb#KjI)05A`C?B!E;avZkD?l>d@8_U)h*)bBp02>svthP1k4S!JIS;s9v-p(!+c7h46K1O!%tExPyk-40oKPCbDI>?vtrYzb%K>}vV-!I6>yPnc(YhAW5(;_y'
    'v^LLXX!dK(Z<Iij_SKQJm~8!OcGeVdqTy|p?$;s~t3k@(?x?z)zr2zNE2vDeq4^T_vuz?%S784Pb~7Vz&2XmqFeepq?qS#}JAq5v_q4Nf>m3%uCd3iS5yCd)UwD;)d^riJY@~G@q&_w&Fk&dX9ZNy1SPW?ZiU6gs#ZsvV(Y47*!P#_NanLQADDZ^J;6}9qZF~E!8)JV6DRscm^gmk@7^%+L)f(^^AWJpXC4U%!of1zb&J8&?'
    'Nm+dtr*<!}5I~a33?cS+t&%wX4y#0YUut$RS%9JZo?FRg?KI){n5cBYK0yOq0}3b|6cIc-n!yw;*TZx~k2qc!mgehan<r#llcsn7l_7@ub}e0k8~k;4*Y~6O!3Sw0U+Yw)WlpN(;Vo1cq`p4axr&Y*-%cm0A?q|>`tZ%#RpQc~!~w1}QFF`fR7%bH3ju_nyw%58V$wtu%Xa7{EqQ7qk=NAvK!f4#U83O^W_3Hxi{vtf84P!I'
    '7z=6wrdS2jL9pXf{(aVrv80rWQi;xO2yS#d|6`~IkONJOplG}b#?#_zRnqnPv8jGdC7=-va}ls%jmlHmXl|<e%U#8xcH1^?b14Gd^e^pL{&>T6yw37N_`bgbMR4<(wV@Ny%^C%Re|lJ}3>}L2aklDSA0Qco+wy$W!qE1p(qOEo^&bjpsSI9tQ%JXJ;1JYJ+r4}mYW!~LzFECX0xo*#9s7}dyB!cXaqoU}ZEN39&rVwNP#8Yx'
    'o~Z<qsjFk&(Gx+{-}_2SbDdIp_2_mSod1s`1m`0*U!n)wJK)ARNM0?nzlbt7`kBKU3G#g&CawF41VKDE|6!ZHF!IDsYgUS)V>)<aH;bcnu74~=pwzb?EhP@?X4&0Pu=%ANqmvvQUC17xP0f&l=`eo-<#j5FBOY`hT>`Hn%eN#H|H{XNnTv-A7S3L-VdN7vPm*j}X$(l%VZkUzBJXW8sj*`It^2a|=+?|oDIQ8V!nw9Rumf*@'
    'L7Kf&F&wa34s($QwKJMWIw|zvvXDYph<n*Xa~eMo9YNX2=G!tCRWTA&M{&;OIrcoP)7a_a?jGJ?P;G+p*yk9*A96AOa9bvl>1fJ;Mn6yB3PL*R%i3(wVPyODnS0>UQ+qr!=A>pc>;Y&xH9Y+Q9Mt>6+?T}Fb!Z?bj>?D=?1h$tGL4>D_QO9bY0$0~beN&srE2Re4sPI{!>&ouH8%aS!3$%3Ypz-U`Ht2J5>z+kViM%4MOIW&'
    '<=<a}?h#HOrg_8(McyZ$W<f!UhdXJ7l^U_Y+(nma9hyfviC(k3%<o_y49c$A$xff*!9}Mc77?tK<RQ<XzV`jYY7d8qwbqeRELk_m^glRQvd7U;M!CGEUx?uFlZMi7`WqH)UJG(#?4ZAu4k|Q&_rzbO=UA(aVtR5#qY-qVdZ@hbd?Di_OQ7o&3+EKRwr;7Jz9HVNW~wKgLGW8-&Ip%8%xXi;92BdSz7KkdD_27Z#^47FR~(i1'
    'esLM$<gEMs?6ugua6*3H73;7|Iisz<YPu#XHLaD8LDqV7h05ne+(;H9ty<nZ*S-;1@}7;sNc2;gIt01?i*wg%nuK75z49=*ng>hlHx#-2`c06FWpa9quaaouIdvCks*ncJk@2f?yn%?Ef^ryo-(l_G<pWU>Dd~t`L=P3a6}$-hl%4LjIj!=-ReO=!t;x_<E5_5aC-?5kp(AwUGMxI!^<3gybCPJjlZpRh-l_f8O9;TjSKO~='
    'P?S3d5iQTOu~id9;EE=Ap@TluMENk!<%h<8lAVS+!B5_MBM)0jVy0V=-vGvmWhp5v?mcjnTXeSn#|z$Bex?9RN1g->A>aPEdqcn0cZ@2?;wR{Wo*%~)T^d=7duHNqj7K*x3Ae9aU*-0U4!R`|iGbT@$;!Gpt3=jX;VCu%Puh$L{M1z@sKYa7p}$aP-PkQ9?uS;RnF35Y_fAuc*`TaF@q;$t14#R>-nv>F+Qy9ao>XG;+>7B9'
    'I)VQuMC@VB>SCS|`RnqX+<;zX+-gpIk*c1+2H`gWB>1)={M&@A^-W*a9@>)>s8jaP-Yy5|I58-tLZidwZ_j8vyB7$ZV|Y0^XfPT0F68%nsz5Z75G=(b?l;h#cI9TTyP>o-P{Wd{2XtF^AlAbC66!E5%?4%(UL+??SyN_$@%jT$bO!?pWvg-<dNwc52@^$g+L)xMLfdFde^qi_xoSC$&3(_rR0YBsMy;?9T6Rl^E>)oVLKohI'
    'JR!cdNut!QPbdCNyk&4;eoU|^AC$blclE_}BfZGgiMIt1<?AK#Zoaz_-zS0EVGtz5<x|HnqplV=Njjal!gOs!`6Sj{71<XIWZeo~SU7OOyd+P+B|4gFv!Dq=Wdgk<dgR?(r?+Rl1ddF9p=KxHiC{qpr;<Aq*$~gfQJ#kGez$DEijH|D45?~%gyxTkoXo$jmj|};h>oohZI!sKfhoHe<53x>+z_S1^m$jiR2H+wJ6>Bg`C3|U'
    '>n$Pi>bNPk&rl=V$@s1@$vtN%goBLaka6l2z<JeEeADQO?-bBZUgf?0u1E(pHnr-cY1Yb6VjU6S;zNJ<B?>@8VHrS#H;ef#q#U$|!ajFsA7+d(6z$p!2Y7UT9=w@<{^vApg>Ru^l<jrkvZhseMh!|Rrn@gouBZNo)kr!|^QfMII)6v7hI}T^M7fcoK;DB<i65lQ%IR8A-s}lAZj7aG>2Ne<ifzcA2Nio}*N4f1jhT}mb$d`m'
    ')vm^S!kNHth@A4`TfGc%72?~Mzv}<*_~mb+)_$8KW<T!v=fs50X5aajKC2IV*<j<7QF-5Rpf3OwK5YwJw<rHtNJ`h2#~_1a(*eH-L5kN5oU{$R!Po5IG>kZ>o+vKs?Le>!M!R)uRpkM;u~jGPE$kFY!Y@b(z#1Ga(w_kNzpJ>(Jj*agB*B=o3#beKi)fZ)G)ih1n7NKAl3_xyA(5$(uZtOendB7MjG3Kb#OF$cWt+m~>9h2n'
    '3%JDIqydX3kz<Nr0!4N<i-^x(d>zwk$R{yJb3$do2%!CdXgA#1Q`Y{vCU5HAlFwTGh{>tSUvGTx1=!{NaXK5{SXXyVR_e$rntF#^(7_Ww4*mrySDQFr?S`IB`Fk^a=QnfE;<N^Qqp?#^-Q3NayEE@=48dX=vUd%?YjaDd2ZGaJnq;YqG>!KA@ZSTMp%<c3%u%OOG9f1+0JGiA+<&v8$UWQ`Tj4fu0RFfJw|*wIsX#LVzu6zj'
    'z`i)EnI&smA2ZGl*jf$pg5($ueR+spNV+GR8R*$`tRgl5N|7|WpfpN+J6%LR7z<;4JWC)cD*ZI}_;FoQ6pRU;qq9r+)jm5G2!)X9%CQJn+}x)c4Wqd;OaGjOIQKYY!-=e<rwdgw??)Fl=_;>}plv`Wb&cuRdI)1<qABqFmRA3+z*}`b%`ztD3&OVN(~bt$+!wEM@tC2hVk3(f4Q_j=gDPw)R(86uSJo8TfHEs7PDpXJB-4aj'
    'XwOc}(_`@PV%k?Y+WWbfom}<@9r)PAcdgVED=?he7gci8vWP=-?m4tLY{_rpRhjBEDMkKbs44Vk<5!c+YcyfZgq%){-{67AWP-OjHHMQ83z+Q$d;U-}UuQK2?5i;n5?=>ODh3SLQYSuKLtYP@9|H;_3cqc+YS8lgxpPxm2aAx&-19MgdvNb4ekhPylr~47?KVVNHNhEe)q>5L?use-MN-HyX<O}Fi#Ej#d$A26A@Z~Kz@?Tz'
    'EN+;eGjG9N4-?ac^v$mAi&Z8t+ytoFZT$T?=>2^tF}U@SWqK-FiDWNoKT}A|gBUyUklj|ZL_oa!RBSd1%K5VU<_n7ducz)lWrrZnlt~GWl;K+9)*jOF589oU$|41OmT&4FBlC8E6vXjw4<M$N#IzwH?gbXH@j%ZI>=gm~9@v6{+(UVcp`+28ZfIPk#$EHr*QiGAljL)v!|D4M#mxi3TgFj1!HA7u`$a$7;>1H1O)((xnhG}z'
    '=fWHnMLcYCBdOwPbJPye9%xea@%BlD5b*~1r9RFRu?SvRGzX(R%^<TJ{rLJd#}ogfJi75jbzu=oAFsgaT57_L{e`#5L`(pj;(E0Wz6hS}+HRuib4d~|W-6*qy6-vJ)0zzw2otn#T~-mKqb1ts)%=51A#%z*-<bZhV?CLagIn6InRYRZ9($wz;3A#ILp+=XO`4W>0H;3K7<#8*GY>l)TZQ3JFOAjbnEldpa#l#*HApm5SQUHa'
    'CW}im_Iqz`2g3>}g_f7ORV#$<k@}Hm0S_AtseA=E=A*I%iqrfh;^c{3cPT0%joK;xWh}?t^<(yQcJyzt2Yk#Q;wCvppp5%#n~-Gjt+hz`JhU1g_<&nAu<5{_I|v#BQeV*6pb#Em@&3?m$_9j%1I;)P{Yp|2Ou(M-r4C%2gWQX3uw{`c!R~SP=zG;%*r`+RpEUPGSU<s~8e!y>>HiHStND-U*4;LZE)GfXTTR4dq&=H+IGNk7'
    '8yN;ZTYmCZKuFCNJstUS;DU*0DGG^=8+*tewT|1RK!5%XH<YUtc;{qD_c*q`6%Z|@Qw2}3P^}#zJQa7&Cs%%K@ihNIySkT6$ic!M4bcdnf#w)yd&eO}*ZM2*Y}2@apSPKsKcgyYggYH|P;}!w5REdDCZQHPW5sMIfTNMx!Fd3YX}vI-^4oHAej=cv3CubJQ5>GA>|Wt1v_@6zo)Aa1U$e~vYB&iX^bcVO2KS(%V56mf+sreg'
    'KC>ujEOkdcp^sodnel{?`$ok;24v#)-6+S479H4!pwmK}R;BU_M2mX~rb;M%5~Y5X39(^@ZO+Or^~q2YBHOHT`XR*_CalVIms5|!Z`8`sXPgn`YD+TrWOp%ym&NEPzw1Qbj!T{}UQZ)geRZKYziK3RI(;$(4~UDfQiOju1u3NhBHW5!GO{zRbnXG}lS$MtA=4Ciea%_=ebMGw=?U}})${`hrxdRcpt5n?lBDIEh3l$+p7CP-'
    'I$>=Ka(d?YFI$&Kgv*62Rx37gXUMjV2Rkk%ZAlJ?_pc81y>+9D7R>5-6J5ERaUya4_ZCx3a}errNR_^ww&JbmTQ5rsMO}W?Eu7|JSub@UdAIf-c;_KT9IGvBJ>+u@Rv}rD<yo0Z^zb1FPR0{n$VG+AF%^4<5Q<d){Nc&Us?{OD3=GgUKlsSBDPj6#jfu>C_5(6V4|hd-IakuS?}4oFl!9RLg%Exd*^<W7<(6fK$uOH~uMX6+'
    '?w@ZP`mnstSa<#*^*vvi-wB4*Vc*2>sIMqh0FPlxYWbsSfQv)&v1RuESOV7}&FJyb!4cN+a$W|Xt9*a<IpF4G{B+L~h|Dg0C~%k`_WZ+7Qre74@=h~20uAmHng#LGB&YskI2m7#@q)#9`oNYyrFkHnW1M=n)1^rD!OFu&sOl)Y@B?kZPzqO8Q#~B7rP{ERe9_qp=ZE#N+27Ik8Q*3fkWTTr0<gNJNh7ibZ3C!29jKt1{<0@q'
    'NZKEPEM}KOj9LR{3MXK3_aB~R3n1CObtt^|sD6jjWBo`(vdd2+`kM;1SuV7mgm$I7fOWMbv}U^NFwE-@?r3I*U1tx)7xgC-AErt7!kPUW%hhWqkf>MABOT@BT!*tk_!<@pT?coMm|JkO+?MyWYX*Q<9b_QFG<>iIlE`7^2+vMd>JBA(Zv^9PvBP{pW`)Ti->ceLFefKl8#e8hZL_;_T8TAs$O$I72txN}#Wg}hUi5Hy@-FWA'
    'dTIjB(@P~v5DxK%&;`eDsu5H6>7^H$VJAX|w*i)QU09`;EW$dK`4N9qHbcQ^lJGOW7uq$0L@}`cTcAUo1RI5|I8|WipamE#HN-%osn3N=-<`hvFfcVWfqn6x=G<F6?&i*3F@pTXv=X7W(B;tz$H68);z0B*1nNqJUQSaE`Qa@Zam)vgpi@QjHoIRvDR1Rw#)i*ik+Yg13+dBr-U}T?<iPudYQkyJS-_`UA^X0J(ji61b~^Q~'
    '%wJ_drmwlKn2Ea(K_;g{t3=M0c|f=IZM%MEP@d(;@TO(QK%(u&I$e#HE=?x&FF;o$fjPZna;cIYwnb`U=2Fx^a%5*128Oe{b>k6(v*OuXM4t4zeBFHS1^Bt}g0cv8NX8uRdpVuo!-f>VU)Ev+j5rXM@JJ){|8uCxVeg5d#(NVvd;4Y*ZZZ;PjmTqsE+|<5=5xD3oeoB>d9>ZFJp2u0Y2c2!(_8wkN(2%6ksx7ox{&0c_8N?`'
    'A64KsPEXe4f|Ji73cna!;lm;SnfmkcpC$n=T(G$JO94S6usm&78Y5wJ_Lt69w>43j6SoF0@4E?k^rjp%<M#qFf*}S`D7_dH-M}PnsMI`_-j0fN1jX03(ve!3WhW{1flV!#AQTf!n#TY(S!TMcactwU4?JSm)Pgb&HA&tW_us=}d0GitXbCi})QIhzykFDnL=zwv!p!h*>I&2+OSoU>e}X#opM%HrR8edAhP{W8y1GRJenH@R'
    ')OiTQc)rOheL2Q%*04}h-zB5gy{_H{xRJoEup7m~z)@2rO!E>?XcE>GE=1?J18RF9>D8#jlfNEKC%V&{^Ik7Mc^zfl*U{=2&!y?PoLxwYjIzWh&Q7?OFV8vG4p47Q2yBIBqD|)+v_)_{QEVm#CU(f<eFyEYuu{SCI=d9Qqsk=gf(&emI2+0WWk$(7e=VapFY$|+qm=>l$mV&L)OGeG5nTi0>ju(0B(pZH-AT5+RYmaOok$mV'
    'v1hdO2R~Qz*txss8Y^Bn{_Q6qs)M}Q)!KXTR$$spnQ*d-6MH{u;S$7q_?l)t3RWoH4qa$(K(>0|0O9_%a(ePE>uWsanEI1ZFhxhZW-o0NGdCC*=`-7e!o>Lbbz<g0t>u18N+~6$?kXB4G%ue`!G0CBc*6cZBLjI7a3rH0WHg`k@O@a!g49qer0a#mbL%at?e*`<;9yGz<p8sF9$?Cz_cIk|%rL_E_UkAN>DLu!r%MAUHlcdL'
    'Z2U<b50K682&<Y*G(zM~wi$P~iL!FTXmKJA+>7)<NQY0EZ&%6ym7y6MM6ePk(O|Wo%v8w1RLt!>mjx{d`!Kp~2E^QN?%jW>>7YLE{S&mCbbbyrF5'
)
PREFIX_RAW_SHA256 = 'b6ba4ab19ff20bc7b79adefe302ea2d8baa2549eba7151f21f6be1ab1716c8be'


PREFIX_SWEEPS = 20
SUFFIX_SWEEPS = 32
PREFIX_UPDATES = 6030
MAIN_UPDATES = 15678
PREFIX_DEPTH = 1 + PREFIX_UPDATES
DEPTH = 1 + MAIN_UPDATES

def transfer():
    return Transfer(225, 225, (117, 270, 4), (126, 286, 4), (144, 356, 4), 580, 710, -14178)

@lru_cache(None)
def open_module(name):
    from pathlib import Path
    from runpy import run_path

    directory = Path(__file__).resolve().parent
    if (directory / 'files' / 'moments.py').is_file():
        directory = directory / 'files'
    return run_path(str(directory / name))

def coarsen_prefix(m,p):
    f=2
    m=[f*(x//f) for x in m]
    p=[f*ceil(x,f) for x in p]
    need(all(0<=x<=Q-y and y<=Q for x,y in zip(m,p)), 'coarsened prefix tails')
    return m,p

def decode_prefix(conv):
    from array import array
    from base64 import b85decode
    from hashlib import sha256
    from lzma import decompress

    energy = open_module('energy.py')
    payload = PACKED_PREFIX_FIRST + energy['PACKED_PREFIX_SECOND']
    transformed = array('H')
    value = shift = 0
    for byte in decompress(b85decode(payload)):
        value |= (byte & 127) << shift
        if byte & 128:
            shift += 7
            need(shift <= 14, 'prefix variable integer length')
        else:
            need(value <= 65535, 'prefix variable integer range')
            signed = value//2 if value % 2 == 0 else -(value//2)-1
            transformed.append(signed & 65535)
            value = shift = 0
    need(shift == 0, 'complete prefix variable integer')
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

    negative, positive = ([2*x for x in invert(offset)] for offset in (0,conv.msize))
    need(all(0 <= m <= Q - p and p <= Q for m, p in zip(negative, positive)),
         'packed prefix tails')
    digest = sha256(array('H', negative + positive).tobytes()).hexdigest()
    need(digest == PREFIX_RAW_SHA256, 'packed prefix digest')
    return negative, positive

def energy_evidence():
    conv = transfer()
    decoded = decode_prefix(conv)
    m,p,updates = conv.serial(*conv.entrance(), PREFIX_SWEEPS)
    need(updates == PREFIX_UPDATES, 'prefix main-row update count')
    regenerated = coarsen_prefix(m,p)
    if decoded != regenerated:
        raise ValueError('packed prefix must equal every recomputed tail')
    moments = open_module('moments.py')
    gain = open_module('gain.py')
    energy = open_module('energy.py')
    if energy['THRESHOLD'] != gain['T'] or gain['T'] != F(-7089, 20000):
        raise ValueError('threshold agreement')
    rows, initial_energy = moments['compute']()
    cost = energy['energy'](rows, initial_energy, moments, gain)
    initial = energy['initial'](rows, moments, gain)
    if F(cost['energy_upper_1e12'], 10**12) >= F(61533356109, 500000000000):
        raise ValueError('strict energy ceiling')
    return {'directed_energy': cost, 'directed_initial_majority': initial,
            'prefix_sha256': PREFIX_RAW_SHA256, 'prefix_depth': PREFIX_DEPTH,
            'prefix_sweeps': PREFIX_SWEEPS, 'prefix_main_row_updates': PREFIX_UPDATES,
            'recomputed_pairs': conv.msize, 'arithmetic_backends': sorted(ARITHMETIC_BACKENDS)}

def terminal_evidence():
    conv = transfer()
    m, p = decode_prefix(conv)
    m,p,updates = conv.serial(m,p,SUFFIX_SWEEPS)
    need(PREFIX_UPDATES + updates == MAIN_UPDATES, 'total main-row update count')
    truth_controls()
    literal_minimum_controls()
    upper = root(conv, m, p, 400, 400)
    if 2 * upper >= Q:
        raise ValueError('terminal majority')
    entropy = semantic()
    if F(entropy['coefficient']) != F(51460722621771815621, 10097617584151289100) or entropy['checks'] != 1625:
        raise ValueError('joint score or dual accounting')
    return {
        'threshold': '-7089/20000', 'graph_cap': '61533356109/500000000000', 'coefficient': entropy['coefficient'],
        'terminal_upper_num': upper, 'terminal_upper_den': Q, 'depth': DEPTH,
        'suffix_sweeps': SUFFIX_SWEEPS, 'main_row_updates': MAIN_UPDATES,
        'prefix_sha256': PREFIX_RAW_SHA256, 'prefix_depth': PREFIX_DEPTH,
            'prefix_sweeps': PREFIX_SWEEPS, 'prefix_main_row_updates': PREFIX_UPDATES,
        'semantic': entropy, 'root_truth_controls': 72,
        'literal_root_minima': 162, 'literal_NAE_controls': 9,
        'arithmetic_backends': sorted(ARITHMETIC_BACKENDS),
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
