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
 def apply(self,setup,m,p,n,u,plus=False,previous=None,round_number=0,Q=Q):
  isaux,rows,pack,prod,merge,full,norm,total_weighted,parts,cnt,fw,gw=setup;outm=[];outp=[];both=[]
  fneg=m+[0];fpos=p+[Q]
  if not isaux:
   gn=n+[0];gp=u+[Q]
  KK=(WSCALE*WSCALE<<40)*Q
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
     both.extend([Q-((total-y*r)//KK) for y,r,total in zip(dp,norm,total_weighted)])
   else:
    _,active_pack,on,off,corrections=row
    gm=active_pack(gn); gp0=active_pack(gp)
    Z=on-gm-gp0
    ac=merge(fm,Z,full-fp,gm+off)
    bc=merge(fp,Z,full-fm,gp0)
   ma=coefficients(ac,parts);po=coefficients(bc,parts)
   outm.extend([(x*r)//KK for x,r in zip(ma,norm)])
   outp.extend([Q-((total-y*r)//KK) for y,r,total in zip(po,norm,total_weighted)])
   if not isaux:
    startrow=len(outm)-cnt
    for at,fi,gi,mass in corrections:
     fm,fp=fneg[fi],fpos[fi]; gm,gp0=gn[gi],gp[gi]; gzero=Q-gm-gp0
     dx=(Q-fp)*(Q-gm)-fm*gzero
     dy=fp*gzero+(Q-fm)*gp0
     outm[startrow+at]+=(dx*mass)//((1<<26)*Q)
     outp[startrow+at]-=(dy*mass)//((1<<26)*Q)
  need(min(outm)>=0 and min(outp)>=0 and max(map(integer_add,outm,outp))<=Q,'tail simplex')
  return (outm,outp,both) if plus else (outm,outp)
 def evolve(self,k=25):
  # Analytic finite top entrance: all children initially +; first output is 0 or +.
  h=Gaussian(1,26);D=self.d;Y=self.yd
  m=[0]*self.msize
  p=[ceil(Q*h(57*(40000*(41*D*j+13*Y*i)-self.T*41*D*Y),13*40000*41*D*Y)[1],1<<26) for j in self.y_axis for i in self.main]
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

    A = F(57, 200)
    eta = F(85971514387, 250000000000)
    beta = F(7859207533, 125000000000)
    theta = F(382941983, 10000000000)
    lam = F(4145754407, 500000000000)
    xi = F(7925815227, 500000000000)
    kappa = F(7021763009, 200000000000)
    margin = F(1, 1000000)
    zeta = F(953064167, 50000000000)
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
    p_anchors = [(210653, 210653, 210653, 368041), (218988, 218988, 218988, 343036), (218848, 218848, 218848, 343456), (165545, 268104, 165545, 400806), (175325, 272491, 175325, 376859), (172395, 276587, 172395, 378623), (158972, 207369, 207369, 426290), (167179, 214487, 214487, 403847), (167179, 214487, 214487, 403847)]
    for b, ia in product(range(3), range(3)):
        alpha=(F(0),F(1,2),F(1))[ia];o=int(alpha>0)
        cu, cv = eta if b else F(1), eta if b == 2 else F(1)
        cz = 1 - eta * o
        du, dv = ((F(1, 2), F(1, 2)), (F(0), F(1)), (eta, eta))[b]
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
    # Input-unclipped laws have both marginals in [19/100,81/100].
    edge=F(19,100)
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
        cz=1-eta*o; weight=1-eta-lam/2+kappa*alpha/2
        if weight<0:raise ValueError('entropy chord sign')
        u,v=(a+b)/2,(c+d)/2
        lo,hi=max(F(0),u+v-1),min(u,v)
        for _ in range(30):
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
            checks.append(exact_bound(cz+2*eta-xi*(2-o)-zeta*(1-alpha),terms))
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

    q_anchors = [(602501, 912011), (575718, 945863), (575182, 922703), (713482, 989137), (676098, 994135), (681283, 998137)]
    for b, ia in product(range(2), range(3)):
        alpha=(F(0),F(1,2),F(1))[ia];o=int(alpha>0)
        cz = 1 - eta * o
        w = eta if b else F(1)
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
        cz = 1-eta*o
        offset = 2*beta-3*theta/4+lam+kappa*(1+alpha/2)-xi*o-zeta*alpha
        exact_test(f'rank2-{ia}-clip',cz*A-offset)
        exact_test(f'rank2-{ia}-unclipped',-cz-offset,[(3*cz/4,F(3))])

    saving = beta + lam + kappa + margin
    delta = saving / (3+A)
    cap = F(7695778873,62500000000)
    coefficient = 1 + 1/(2*cap) + delta
    return {'checks':len(checked), 'local_reserve':str(margin),
            'saving':str(saving), 'semantic_increment':str(delta),
            'coefficient':str(coefficient)}


PACKED_PREFIX_FIRST = (
    '{Wp48S^xk9=GL@E0stWa5C8xGsBn&k<HiS)kpKWx9rVJPD5#KNak?4E#$HaxI699Gw_{_l5B#2NjKwB`8Z?}Xk?|$W`5M%4cxd)>KbZ>d%+bjM=3zjEY9jM{3AgYBL^5T{;}ZG{m#3Ow5!gKuP%_qKDH2mu9S?(pNh1q&;jC@b7LyQu9h*u2hxsn)Ds;~=d3}ARI(=W)r&CPu_J(<GfZDbs9v2O=MHm!m6P?P-&MZF(B0>UL'
    '@dpLCbZ4QLh~6`;ZS~LBdY5T%Bboj}!n3g++^iJd3Ru}&gkhWN0okWD)_qYi{?mq%;PV)sEGGQQkb)MgC=`P@wXdWQBCelAY*}>E${s?ir^QyD=_E&>oGZ{jP54}w^z-AdZDU|BS24Z4^f7OE{r9++gs3NQq~3YVF5}IcM&b%0`{WD#7x$VgutDNI8*U;CDatZ^F@QYIM3Nvf(GMm~jqQa}8nBX9`brnnJnjhhZTZ;gFJqCV'
    'A<exXd%6<9ve{-A<8Ns_xmnd!j?#(Po1DH*mWCER@ELK2iF?2otZ7&{=%%Mdk7ccj81nw;Y*M(cRB}dbdBXufpHuc~J+my^uoj+3WZ5jfM9oS>H`TBs2o+M*;Y%yy<Q)^?k})P_Kkf0J4A3!nK$S8SE13lexg}0J{E;nb9&;oMyR+MSM|T8}WB3WZaL;U->pQnwX{_llNsM)3S(osRqPMs5m!H=9){&gZW+K0KKDs;eFq*+%'
    'VsBliUXCMc1?`+cH3?UT`V6<6-l0wInp*tN)>Uq>#z=#>r3Q`sRZ!&rXAG|)mH^~JQ0C7pDbIj!(9AgR^u9U!f8%90B6hMv?v@_qv=ewf?Og2r&j|D6($OWpEP?_SfK#F3##@u%dsoDw`7Ig6F%qtx@m}+SfM^DgTLi)>ntFFMsmr0<?RiysHryl&?Uw34{(>V#3JUh|r6ip+l{GjU!={ozKgU7ML5M#B&`(vcj&08S*sRl+'
    '>G9#De0m?=C!U_J0~mOCkU#tshU9C9)m5&KpQ?RiBlXThkTm03Srjz7m1TJD$sQYqJUv-;eG@FbWK?ZmR^9m69t`d`D~QF|OqC(wDrSZCVOI8X$(s$42eMCx22dxhJT6Q<qEXk*2G_@0AZKl89eUk^Kphre2JdC4_8lj}8+SP9ksV!?MTuSVPc_vqM@B`5*^s?Ui2xb3tORPApV%iOpg%q0&e#}4nI6!0X3fRt_SmFGn-X|}'
    'vy=U&=1gZ&Xn#onff0leAtYZ3o=6=#p4PFsV)a<|4pw~eJmC0{6Eicg_>SmRj;ht#>mnOLS+XRmbE@Q+8~)L4D^VkOscwWYPQr@X_!!UV#-E2^8}<L-C~qq_E(pQ3fwqxOGvk|nCyVJ`VWiBg`abiV)#HM!LzN!}vWd*~<8In!0IRS*ZZ(W>GBC%TSt%YEE5opA_VI2fVFYpyTFa!gUKoO9?;Gk;%&1SOh}k2cw~IqqG;sjH'
    'u`|wFFKE}^{sN2x!S$dx!un?~!D%}bz=|<mrlQT@<+oKY(83mwN-iv0ho%B4e7)TY^p<HNx7zkZrCJ9bVMZ%>p^f|Q;3{ti1AD7iIefN%OyfD`5DH0iVdoJ@oUG2pC%Hx;$!_UB4UAm?f^{qnMzQZs_n#ckBXUO-H4^Ce<AiIS75T7`iyO1TYQ(c@dd46^&bjexKDw`I5)hy;lLKI8KPiX)j2g81G=T#3cZunDJVW8ZjQXG4'
    'XY22z8xM!hE@>S6+MHRzp&#k>7UOfHO8C(=OFIY+<rA9oPiPsKjKIzu(>B#QsD{coWhRR|Sx&Htcd$G#w}-;B>kaB7lJ6LB*FgiBeB$#a#*IJiEE4=2mR+J^)SeNgWjIDO@d~dl;87h^^1LB^(guB_a+Pn*z8B&}saXTp+U;=Oz*4Ob)nYWz{uRJ9=yAV_dqwhWDDbfBhS1Dy?UgKTJ#roE<|b^cG2A}3tq0`|PiM!>Q`}65'
    '!8|R>VQ!)Vv5uVgJyk{wy<Ar{Z;Z_CFDtAYIA{EXAE`T=D_`}A9?J`;VXYm*wd5YiyM_xg`AwCO>8$a<D!tVG%Y98LLBohj_?2|i$qBABD((u2G6sN`fJI0l&sOqY;S$>j*DPpxiy^i?fgpO{K=|w_PAl>nI##qU4msFnkUL_A5NLKW7h(RMw+f~aOgF8Eb*xk|4f9j3G#A;UOmmO$bV93WwT)_q*dfeDP8D1efPt;}<3z}k'
    '0<vL2lZ6V>%5n@0t{uNWBTd((`p9UF&I^DO_WY+7!Cj*S@@XFo_iyaKhpe=#$%q|d{kmZ?>LzDJJ9V=N3Z><iq1R-zw>Rsxo6y0{yw4iGakOgbG4vD<ggI{!dY|SG>rIa}nIthy8-*K`N{}Br`k`kG>kjHBK5%B7G6REJ>BeozgrM9b^Oa53smMl&$vFfg?|iq4aKWojdi7Q3PwD0UW6a+9!BeeNJsu|Ms}#se@<do#3?AnC'
    'SkK}(Ooqycf-&;Q8=TjcJ!*T(qz$M2l8(0<RljQ=?CRO#MTvbR_H+TvWDU1vgmXxa<sD5LH0>d!rxz}=JR(F&X67WxCijZ3cqZyakHHndr@5pR74+$N6C`V5MSSPCvY!<{_gdvlD#P*sGFGH&8Qh`)+BY&E2c4CIatfZt@C6-N`YheH^T7?zedr9{afccD_z51Kk?LdT4eBLegQ4oTvgIUb5@ib*%s&b&xM5ofUBI-VPpywW'
    'yM3KiJ);3?%+nV;*wcSmtQr<ADA3u9jyCe+@}u7)I9oSP+y|MY0~qpoWfBSXic-y%y}F-M#sm{zyP4>qhAO(C*taP<DJJGncWz8KeYI+hFl35|UT|Fl5{1*#P%+^pvk?NhuKoU5;~pOSxGsyEhfVN<T>q+yCknp2Hw#9zSvu>++0dGl_iVfX>Oz)F(!^fIqN{12v#k{GEt3aAB(C3Zr~T6^*GcD(FwygB%U4wkH&`1Rja%N@'
    'tX)ulv#2;JftX?*5C$58f6&Q>{+Q29IER`|u3q=4;{P1!TYwr>6KhIF!=Swxl|()vLpJ;S_@zVzd&2FGqX*Dq+?Q*1i#?;khbmpDu00bo9sEL_0@k!2V<mTFD+2M<DMr+vIpPB?xfIb}b}Oi8XIIFT*BVlu?&_=jIY#HJDKS(CYr_V0RjK$hr=lgY>z9BI*cFsV#T(zd?wXP1h-wg5*8L1%Tg(c;w}HtF_50AlpYD=2!ObbE'
    'H9-w|U?}Dm!Q@qo3)X>OGpzlM??n6e6FkeB{!p`+hN_(mEWnj~G>~jg)>{A_Fgav8+j@;>DX9l@Vv|;y3$5^Dbn1}^S$`p_l__=N3Yn<9U5{~#UgW|bs3>CG25NXFQV<%u`9*b|;U+N+@ezZ|SHdsk7c$4tFzy2(c#)7U3P)u<<7;-WzHE=EpoFtRR0Ud6lUxc&u>(7JLScd?tZxR@dc+fpI+pu<D7Q*TBj2vrqs0^Y<q8MK'
    'G~z=`=^>z#Y}-{yS)_#Uu*&AITK&1voa(y@kerg+lt&HUir-z>i5lfGJas3FYBUl>C?ug7%8FJVx0uLcy`aDppg?YT^9rG1gbrNR7Gy-`ksEKt%Xwn;re}WdhwVg;aK0Hl6Oxyg?<Ep0j)`=)k*iwwB%T5)zBMd{K!5!}8uuN9Ekoe$xXCQ>E&MFU4Pk3%ekKXk?VPZ|cpYq4@m}g`xG`M<RUKP#!m<`*|AUo!+4$XHH~+95'
    'Fx9h3-*;lfCvAB-xllc8XQY#3mFcaX)I{{jhOpw^U#R2&3RYJ?M|q#RZqgX|02#BnJag}eAGMkbT99YXT<FoGqKmgF_XBqn#Tbygc&u4ufM@jNOY0&>0!tS%l^8-<X~U8!zk~xFg~oFBW=7raY@NB7rk_k)irXAI*#=(-%v-(T3$+NcI$Cd6YBJWxpMQ5wk^d_l=?qT*;3YD;3pIPZF+#cWj!m`q(g)XEYW1qKJ+rWDH=VLW'
    '4n?=Ww=SyvW)~lE&oo0jnEaRRVqD#3^CC*T*yycQV8|A>UrFOJxZ|j48V|j|GDlNPIUoY!eG~7;C_-$=!ZH2tliA@EQACjh7g{<I4NkB;W{;V8Sou>!Q96V{Tvg(G%^$&g#aMF$avPeFBa=jK5z^|4IX(v=iX2`S$=*I0#na7OHdxtfv^NFu<GG+Hn4m_bmyL#B#RTQsP4;t~#vSDJvn#V$JXo6a(PaeI*%i7&jN=L{tC@tr'
    'Y$8KEHDQ02%%AEiP-nvj)g7oz_}MiD)t22QPGKq5HXj!J47Zn0fW@NY5&`7AB&VY}$|_c9F2_FsnEQX|d{y~{X(;dIShlho5lkn|i|AVtcd{c?uF}N*`O^6v^9@hy5wKy`vG3&y^`*PToUGfv7RPZW`zM620JL)<fhhkDL%a}pflxyXLANc-v2(?@Fn49$Lxal}W=|kW^EB5WA8OCRx7NQ|bBQ{z0J+*)Mm)Nax=}nf+tHm2'
    'S}w3m=^tLd7Ls=VMBs9SPmQ&zv*tPLFBXFKwT4o5`T9u+`1am|Cak!KDczo;WmS$2IlH!g+pSpX#xR{Nge1+*arAM+(FlDINqiJq;*_IFO?Fz0Lk{#T3+@{N9d!Ku`xX|o=w7l)oO#=%R6Og83Bji)%9&BPMTZQ;#^)B_n|yKCXQ3dxgO=QiOpTwh>%sg=DdOJfH{Q%2?2&8Pxs9ny^{Xp16g4fA*BJzDd5G1O>#pngGHyJS'
    '3_5vZ29K%;=J|WA%3k{AIgzqSDmh9Y1#Z$%IbCjkRY&FAB3vvmGm{7}1h?%8L=28P9+Zj993fh<tc?(bI8!;Zf$}%?lg^YFLm#JbtF5<eB963|tt}fjr~YW$XS!5%K9jdn>7QRdu_>P}bJGK4P;at&E{LAbvHeA4xo6(jZ<vv#a3h|2Z&TBw%=4f0n$b0vpun0tLq~G#^)@@$ylRk4xVE@qM`7*E*YKT_uZa|2FwANlX1X*U'
    '?2kIMF_SO-Mh4#Cmh<RgHxCuW1-aZ+TdoSFw;hD#*xX^g<zO8xv+Y-hpg>0kTE?#Imz;oD#7~i2wuG;Q<Ydy8TV|(0v*XGB5$@s5R1!0Iz=V2;g}ZCSo9NE~)7TMDSt>Eji;if9`a?uy_<Udr5?%;kZZ+2r(J|-)pW@aRp%u5}tt{g4&R%c39H>*l@b7MS8rAE=e`Iv)kG-zYn+Z?yL_t5AXcmQ7U-?XEhNU~r=m-g|7I0N*'
    '88*mVbs%YlT0>&)G**5pBXD<R3<s(y5k+H_v=#ruE8x(BD5=L!OJWKYEcd`G6?Qd`PV<_zBF}t!E01w7U{a~%HbAu?!U44MnD?NE8GHUGaXWq8s30?pZf{c}d`F1PIuF?gIi@}~(79h>ex7k2Oo9oiLg6Vkt57Wm#{pM<>aq9vyJ#u#0Irp~S3;-<Ru88szYqKN#6Il@3%^M#)y(jcdC6dvhKxL@vi-j`V-eV%`J=vD-1qt>'
    '9Uo<Yc5{8N+t9<5=ibpV)r}@H>io@Sjd>f39mtdSf`lnMn~875L4@&zlYl@{8|@$bISN0x)U@@RCO4`m#ELAA_x2iZi2dJU$we(~t@0xV-jDxL>=$>~$QQXfxn6W}8l=%*w4L`MJNA(`_Vn7n(c1jVS(RAy--#`8jb7}$iyfcY>>AV)5j|qA@oqGA@w(|hn}=_S039Hw5f%JtY>F4JH^a8y8LOFXW$<8mh=Jws8+r_Y`Sj}<'
    'Ox8=y%j(?e<YG_d8FTEqgtIqfDW?dYfYBAJ)*}1c?G;Q8NssN0F(V@`hwvK0-duH;v$stR=5RDKY(|!~2IQ>1j;KzBK#s!khqi&15V^>5?vFeS90~RpT9_?JIU&P?dn(X`Jav}W#xoJDX)M%R1T^Nsw{a5Pm3|GfpP_%b_LF_UEqE~R=B3UrUVM6q<odHrX>~BQ?k7^%KAJ55j6c3XUV+jm`1VbBGYw3=Fj&d7OGXW5>%9L4'
    'EVyIxoxi>Bk~<ffWFj`&cA&DHcbj`DxrVO<6C$86W6co(Tyu%Q)8wgI4a+<wLN!e|1=uA;9TN~?MCWXCH_ad$g?7SH_|?C{P-W?MEKrusu0jbuFplIy9{Kd_&y;7Tp>tv2&t07`J^$KtWx51~`}I&yW65{oq3#uHR$bh6=(+eGzz#|QQ!6RMElQ+-UNUZ$AOeG|=OSjh$XUcSvMFXxQf_uBO4~V2S{;H-(`4D%T4&+kvY&+;'
    'U*4m4{@!~VOo^P=I8)B;W#(^U$m~5Rj-nx5!%Y9?2M_5vr0^^)3aT6?37Zqy@nYCpohMIkOB1C3#hBIWPDrVYI`XGFvrv>2)J38(*cTnRg1?OmGlS5j20LLvg(M`I)Y+IHGb{vY!V2qDg|o%s)B$DTuuWzgz==&(Tehv)L5-Z%{ck2XZlG{{x{SA}vj<xkK)O*oIkzqY^L*g3GUHt;LENR$!L|9>&4t<=s+g#!+bS{1uym={'
    '(qJcQtsD?)Z&#w>uK0PgE#5lM`>}i=ETst3{xVnvDSNzc-sK&~y-xoWoYtHmIiLtd^by+ucRFDd6)8d91=MD1*rGazz;7xncb=x;7MJDS6c2^2%cDBsP&&-kRxZ&O?jys5GI>csS5VoQy^d#lA0Lj?njUGbpDHxSHlnaG;5y+^)zeJCB)I8gdn)N|CM-A3GBsZ&*rCR(`d!SzI}AUhcyNlaMJ5)ZG3}7O=+%hZG`=1Uf~VBF'
    'eTDm^elh|CF?%b5-TA(}1Z@sat`L$Mm=W1f^%6~9@?h9TdD*3Sowj*3!81LsT9PlqEWJg0_*1P$Mcr^e=yrG}UpVLd95E>-_zjZ<tkcku5VAkhX~e^<R+!t2INvmyu}_{I^>?S6Z?%narS|pG1bkjt&|!n#O4>`l0}%7x;75%~Ybto8R=-P1D_Q6xLr<fjl6;`8fAj{{67g(Pv_1C!SL}Qd4l2X|DG|$bn}%i=u{hu5r*h0<'
    '{f}5_%Qpqs5TDBF(2>Woegkruk>}5fG(3x+a#W%sdCFtoATvWjUAW4G#&dv4lUThHMSpOp0(^S5j1yHuSk;O7nS*@LT&z)c`1QLdPDI~V%d`EAdD^s?)t;P4H{dLNN7FjI5SUy~VZ3FX5+<>ABZOdXa+;ZG3*c@JJ9KW)oOkouHTzqNtAqVFEml$5pFM1=u-0ewZ{DGvW6*Gl3))uNVIsRmFCJsC3ht^_m(I?ktH@jqPqles'
    'B$#lC_9#<X-Wt-xab~QV4pO!e{|-Z-8#}AfaH-qP*vJ3*qj-JbawOOsrjHLzD%{!qG^Cx|AMR2rCB{kd`;2~|7D(*vgPBE8r6~PE_WyqFq7d^P02gd_daQ{}^JlWmHAN9d@>g<_)<8p>RY1c}g;8)bBa82Suz-K|BJZajTyeL9BVQp}sY^yH%eHGHrO8?48!vRx>1S~`#8O_^`e0Y8JE67go{1QXLRyFOA{KFWW{`nXRE{Az'
    'H*(LVQuJRwV!2e(XN84Mu`nPwAqdI2(OjmI_@<vzL&yhTzv)W#FM5g2tW<oEIe+e)?>`m5|J=C3Mr{hx?l1C!q4at@fV_{BuYXe{Cm>Z=Fg!gzPp5Bcv<`e8X5`j}SFxd&V-MjV-YJW+qZDZg9}8M*sL^51?M$5=>e3W1#a`DoS=B*tg6u_lQj616Z^_;Tj~a4RV%iA49l6SOV7jmF(kB!7Da@t1N48qYeCJ~g<UX%Jv?=8;'
    '7?rl{NeU7xd;mj?3S4^!omp8?SKTynP?G)mKu#oV4(W$luVR?+CUXe|^E**;F{1Z8FUO8PTn<4!RFnVWl+mFF*gfc92yyKAoY7p$xb_7=yDF;BTsp>FDV^AHkP1P@$Aj)8-c`332Q)HRiiwQ73!xgsl#~%3*=3iqA3B~oF;Df5Xb*X=o1M#1_4(ZJ)`kA&GJS9hpmAucovtlG<d|~7pU!N=glVrSyb+?dU(Q~w+TYs4*sXvs'
    'QMgFr%(I9owfwaWRA%-ptD|O_@#fA@_my++M=9ot+|(VLLDoS93J}d`8gHQTfiR^|hk|ez=IYA-t{l6bH9el`XDTGZ5H?4&`E2lWIE;WAAz+@+i@xD})VfaadfVO)+%2F`Xd#YXAbEs7f-PcY`k0Jd2;;EqLgffEPO1#WH=YAgAuKh-(|f>{nj_qrcu1>sT=hRgdBj!vCtIXxk36Tk$hzjB^PqkonYn<;Wa;C8o=IDA?B+Qa'
    '=)ETFZ=7XmWe_U)whg2Dcq_(~upR!wtva2|2WucQ>DTs=hE_hqe;wmPpN+DhO_nS$QY~#)QaB*qlP3xPCf_I|Qu%;O+Qi?7*-Jtxq%|%AQ6fnomlI@j!3G-Sx@?70Xu_V#`$M>4<~@vVXxs=R?X!WI2CS|y8LD22-)9jX_+l!IOilEc0?mj_L5~{Vu+r$bjje`U8Hic!bvuH(fRAf@#Ho|@;O#_5sSbblAei1DhcZ3y9Cs20'
    'N9tm*SmWOF_v1}MG&94yj;t_P{J(8=;TB!LEXitO-fK<XzPGh_k_vA@F%NUeahkcy^>qaBy-Tl+Er=|&!xpSB5>byB)m|>w(G$E+0|2vdgmcbV&j*IDZXKlHS*yHNJRKfto)sR#`KT&VfJ$=La_qTMhQkc2V+ZMdLQ*Yr63}%<ZK!c%6aLg+d#lWm?)c%@()@7W-`Es?4UJGO-Mp#&`FKyIyk4~11F_HxI8)q3;C%rD!M)nu'
    's1qF&$=ycPH69ZHN_(_L{(zY;q@JSK1j@AV^o!Ftn~ra~2hk&i4G3GV<ANrs`2!3TT~4k>%RXOr!#B)LNBJFoSc7(uTigwnn1jsA>@Li~bZlyomn&AbK`N(3%}CpK1^d__rw|xD9Gj0trwc=*x?`K&G<1>R)xMi9)n6&l)V>HaP`ZM|o>~BfN!+!O?ZXPTY-i0hoj?zREgivz^pI<!8+`U`woK3yb_ave(P3P3sSDf2O!UeT'
    '!H}xgejG9^9<KH3Tc2&EZ>%OT)y1oZEzty_^*m=p<^ba?x;l5nzJ>c-*;ICK_p$<y`?wxnz=aBJ!{d-bbVu<R>m?F~<5>(rseBTw<b!_n6)+LND>d6($~+c!<*Hwh_mi?)Q!kp4c7%y6$L6p^12)x1!Y#c$2xGUP46J$z<*GWr7^rm$Y*^P{i^+$9#t-bbrM_=XDAM%Tg+7ZJllLEd+BK~wlt*ik%f)NiOJ=0N4HwZ`I^{_N'
    'urxb%ufKgpna=pmL{$fTRd;=o!Bi<<7s*5Ekb=r=r6Qmz=ycG`($OJbtEv5uJ9dh%tPgF{PZr5y`s5X_ci~w!D&bsfG?1bp&iP#=tCaeFAKw!zpu%pN3m?WP9<Z_nz%LqxZhs%iA|g$QnwN5wqOV9@9<jp+`ZfPGTQXv1HUFd%R{d41Ag9LT=5N880bT6x`OfvUR@0B=IV3qp&yYA5bSwK;Fh|fIsz3EZe`v^@W!8x+Gi*Ig'
    't0Z+XGo+%Ox#WkE+0aJ2*Z#mfBw76}62WFLIg;<-LN#EwSS*bNXU0|EJUEX4fPVn+eu&l(e>`x4q3z0@RD}&uJ~!){6Dv}^Ij;Ef8HOXjR#P}S4LeLgv%%jf(b-x4(27V|ju4H1O#gq#18bXjg4E0W<)zt{Nt|gcqhz#m6tU>mthiDxEl~nPe-s#atHDnnene**Eb(N@RnXTQY(!RPRPO8bM^&ko#L>l<Kh5Rttc!K21(ZV}'
    '!mEr0_#-<!^Q0lFkFg>Qo|xj9FLtX>?KYFpn?T)#$m-rt&Dxv+zwY|pG=?SHPTaFQhopg6QU+YHAPsPV=wJ;oDwv5R!_<I0Gb^me+WMEsya))U0cV2?@5>6hts%|vkyeG>ZTw<UY)7Q-2=omoSF!XhIGhq++#cg{7BS~&A=fK`eL0ha-S$u;sA4!)4)W$*e6ONg2Uv?9CV1m2nu=zX(L^C~gvT7O6?s**<qs5DY775(>MRG{'
    '94QXQL{@*~pUe(W*Cq)WSs7qT=&?(g4Zd)Yf_^#ZDJ206xt9pP2L#I2j{b<?DQSF2kZk8`G>jz-C)I{$^{ht9;j=h9Y#`n?0v=3RJKCVNX#HU2qtdUKe@u?zC4L?&k-kZ(IL+rMOM<TIId-7zO&%nnQWGdPr7V#(S%MfAD|k$i+jL(zm-SSWHm{yBKAJ|E&d#7?VR~hLAcs6fS@S4Tos;+oSnsjVTjSx<VCaljYJZyhnZf(H'
    'oPgNaD>>a}E8{mztDn<JbwfL{j7CxN6*Yk0r;)w~*mRBX*Q|*YP{uJ$_<U)NYNeo^FGV3i<_aJ=r}%ghS>?X(L`(5}N0&4XuF#kFPw<O)&fj9Y6`EnsQ%)j`HJmT!%&s}{Bd+V#+$qaxDze!rg#ifhA_N%OJtI2a==mw}Qivd6q+L^wXsBN2*4+M5?z=k=CVJ04lsxUAqvlthP2Ka|q$2nb(7G^Z1z{HK-ACX_cLJ91f+`ly'
    'VgH9#DoS&Q&-P_EBsGaAKN?pv%q+?Wg33Gvrry9Uqr%QH8TT}(gpRn3p-7dp5AWo(v%z7Q5qk1Y1Rt=Pf#BO8C2YJ-u$b%yYHFytudTQ)@!z`JmF^AG1-4mcJ&#b1>-L9@ka2;rfx&BLeoG46X6`B-sFc5(!m&P0kB$zxL-Sr|lZmJ<fynyI@dj<4n5U}hrET{3*8cYBklT6?4kN`rKbtf$`}sta*1z(A^*bWg3KOwra7{PZ'
    'cUH-yj$;uDHGBm4Y)fsqZMt}RAnsvL5aX5xoaW7MWX2>CKGR*%!kuL2lg1~E@kZpJ2V4|2+EsKa)1bIdd-h6Ou>6t4m4}8Ix@HqihnXr~K{P~H&wn$v(5s^5@aHE&g2<~_mNej!=PbdTyt~Ldbaj~w5}0a5dgrTMb5f1%o+*~yyWLeJcfsx23?SszJ6#Lla(c;U*!-esctU+M>D1V*P^=kfoi^azL`b`tVpPkf2reO?6#(An'
    '^)eR&PpdA0(zR_n=Ab-CRs>EmrwBWZV4@pZ%kTJ=q^B_uG@yscRH#MMJ7SMYiN0ep)gfYVe)wM?Vvfi>5uj$Jic1afwL>Gu;4O3V=Pga}f?jmkn|##entXes-5{*D`A0pCjgD*(flw4PT1;FCR1I!hZstTow0{Gi&SoCU9{<a{lM1phoOEVi!-2C%>}kXHK_e{c3n4X`^o!Mj8!lCa?Y}V2YECg~2neg=s>U?btW;V~y!Ur('
    '`q2ED;B2^#M<MTvju8Kz^-EKX^TnN!s!-r&O8f!<f&}@0rveO-b2OKavy7lOl=IaLl4HGc^~@m*>4hz_F$-t7Bve3a_k&yLK<@bKPJ!i_v1jI&Q5^oN9GHj}CG*tQk5^GzCbmIArgjL9RI4NM#cbgj`8@Du%PqW{nca7d+1>Za8w(q^vpz%)r_nMh@UTwsI$FZw#@>UF6%s|r%BG3NCGt1QnT8pavtw>cT?mJ1=_#nx%>mH}'
    'H){sMKdk=h_on$o1NhSqXgkGBiaJqfT_hlNd`I=zUE2yLFZ%z2k+D`X$A-m9Lha-!c<Ct?KE7?oi?o=044j+tjqjJD3-^j}V_7!9&qt5x-zMTXmR9ziB5VAd>_#r<?y?osG2t}2*g{?t6Np5WX(RDgmS3YW#yJCc$wj39i+K2zHztH0e_bz`Ym&M7JiDf^50Ya)MNkL-NhhY%5aEmrh9@pAZl`bHp%RLg#2Tl#t1QS4fT2SS'
    'cTk%pNR8HYUEjPib}pzM!9EnPdjE4aGPzpQ6;oa&U-9#E<&*w^8~@#duv3culwxbotjT?My=p?(55pqQ`+TK?`?8oU8v)isOjuBMmY%iwki84txPPg8f3&Xs7Zut71F%CWg&o8I6ta{I4+~%BBd{Udpt4Uxx`_#`pm*@@!pSKqZIheG)RU=Y6u~SiRG}gyKZlk_htHpNB`!X{|M%DIe0mYFHAE>2=TlgFuO`8i&PwLEx}RVp'
    'ZFUU_7oF49RL64@K)BrURDsH@re>9g0UL3&>^TU|QKd&9a>fbVq}PH>B#8WStPE6?6xvcwQ4Ky2XkPqRI%~`dO_pw4OtV`7dZxRxxNuYE@J0Fzp}EFJcHeu$@^~q3XdB;OwSy=VSSm?Q3Xu`7$WXR|uhu1nPjgfUO7$jfr2Zn;Ewvu&ZMmW!-RHnjY=PY*i|VWCWf}{tg)!90k&Cn97;HU}B&7N0NWpn<v#YQ?^ybGlOBlYT'
    'cTLc*cb^@`>PVrRYeE0_E)S-b`lgRcTN}LJ(LjrWV7)5kb(M~pa(!zucHxeJHrb`TmWMG)`q9Z;;_A5n1+^sh2p>FBsYXK~=`W?*unkoTx^Q`k>ByN{W}vd>lkasP*ZVF;b#kA0IOh7D^|39s7zt;p^tDL;JR$(g#R&;4OzK(k$GhCaz~p9D?+UQI+0n1ofBl-=I%Eq6hM2$>`2)@ncCiT21-D;WNGHkRx?%enexA#Y>BB)R'
    'WFMp^SL|4WC{_84y)ca9-u(!e?O6+Z=a=}5wJ1j`O8`fpR*Cz)%!K~y0^%6m_yr9^ui?a3_w>1;nZFN$rh3^aaXMi-;BWt(Z9FOE&Z;I~Rk!frn)o9%n>*MQa#Ve0cw|_AB7Mwd!>|1^`vG&I@kZYVm|@$4j=T^|=bD|LVVzvUS05Z!5BLM02!$#Ai<N5I#O@oE9Q7||AG=8mk*OGiZShoaY(C``H+TE4kHAbNhBO|TbGv=n'
    '8xfDw>F4)ZGK)GF?uP{1BwksqU2VUn>|1PIyr2+*W(}0>3MRn~Y)a4VW8i$*!{Wyd-}^_WdBkmUL_MfQlmT|4!dlGOc%7d;y(cLKV61|&!gyCp0XT_Ylw3S&bN>s%YCG5{Q4WQS6`@>dBB*!l0a+p+aZV*wpHpt=&2s+3^cXk8VnA3x_3A185lQQ?y~%OBjQkL*cGwB@y(-kdhd~7v7+xqdM=o?B&=6K7xHm?o$p34@3MODW'
    'Ygj?qi(;k&lJulxK^ES5t&sPup-`cSn6Pz{TTS?wIYNwg3ku9CVtn^i_7Mn6O76qXQW{~s^z>uu4@8+xznGUh7qr{cE)HhOsY+5IQi|jajJ<R`mE%`IQrdP-y{Kv3XeP{vCP(3w(te*E?sCy-hreIX*sdw8&*-g)oY=t$%K59$l+0&(<Vgm!5v;5F-a#DYe7B1eh>7YnAVzgq<Aaz7GsP1+e{1A7<9PY|JDvsohqY`+1h$mR'
    'cLcHv5+-hgHN1Z}7V;-3YX0cr?$3rXn?)(p-y$H&0e`4+gU@x=l<tu*)pb5KivDU-81UGeT%5G~nf>Dj^>5Oo;>H+7lgQULYF|{u(%2efN6cfjBo$^e`eB4;`ml;_;>$VI9d)Y!Bs31Ha#2mf0_|x~#Q{O9IaF-hu`dvf!#(8%PKUonJvN;04Nr0an(9jWr>b-W^R(T9-Px2$4tHiJ(Mr;DO`epVhgNdBem<8%5%EeR1AksL'
    'bjWhf!$LmF*pwzN03qFhXZu*!1f^u`MGx=k0GN@8q$BV9^e^3Og|nneF-sqcg%Y7rMBrbH95K6MPtY!jNH?OsJGMXK3tX8Xfc{165}9(9r~PWJlv1}&3qJ2q9avXQ<oZfzG$K?e3$Zud&<1LLAU?DGJKKjB+^&tKHc}C1UZR+n176^dIPEqE7?41vh`n+~QV(Ne`?tA-2rQEp^#6DP4(PTam&_<k{Mjo+7KJkz(=a*Sjfg)n'
    'Fm|DKMM!=-@fRvJsGG^73h(Id)ZZSWUQcfBChojf<UeyNtJ~b6+dnBHVbZE4$HrHqmzIRYMs?g3G7zV$pS$C)L@J{*pp<a7mq_xitsigV=zOJv?=jgZ#zL~wBzJeX0mpumN-v8KW?|P+dzwUzr+XgFN698Z=gG$%83{(JGIsX{$xDqy#EtA*5?Ie^-$0v^JuG`5eRzHw*nzk)t;OuL$(AIQrlOms4qxPjj-3*lBJ8gf*AM(N'
    'tTy&bVIzv?c8BO=i&x)<nje%p<I2Oa=Qfe!hwK-a%JshxQkhQM=6JqU0Nc3UBzGad5-+!}LgBR%Ub8uJ85DB<6~5D%OUJ^|RHmRsptb~Q-cVNb@+c;1{`s$&FFU~JNk$x6L~ENWv`kXD*`jrMbQ67oj*NPhi>J_TFc!=$+I!N9N~2lhadb8V?pHA10h+cf=@rftp6?d{E(UM@P#Knr2Fqv_in428Ub2_?TXxVZ0NzEDB>5^u'
    '&gm(KZsB6X`Rl5CQ``G+CnlW&T9Dx%EbY#muQ?<qWPE#a!M3{q6k-=;G9la-Klf!)mu$|jjjoX>Q`B8n0~1XPizdu6V+(EzcuNa8Yw`!|mWbChPWm5zcvC<k^WT9mL&&S0u1eN16W)ff(B<*rol7(_it>toqL>}mjNJl4Rq`D%IOIyK@BoG4```*&eGwyQEBS(pShI8vZ;hAh4xT@38Wcwle>#p7`^d++4p?6u3zLw!-s(-e'
    '2Pb!SXU=_bTy@kXK`hOjoI5acN#4Q(dNV{ZwpMG}c-z=&+0CPX-Bv)He@33b$iP9%J8zL3DU}lZGN0@|ucZ0nd+{sPgA0)42Hvg-kjd0v)p?mJyi3)$wKdWkW;u?<w-`eh3GO~lrT@oLWl&5oX3d)bpW_*ot%clMDBz3*XVMX9W0DWZ0P~`162ES>^?Q0ksn6}#Hd+ZP6m5~0vQHiLSuKV^YH@GEQ0#U766>^d&jZ;+)A^wS'
    's6{$jL$z%@)}F~%I#l{?&xG@FGSE;O?(QazYzBxc?k6S*DP)U%dvB&U67GusI&DIt{7HGnlVT*GIGrH^%vZZcJPUq-cTNPG$!XP*UYE4`vLRPmv65|a^W{@3I_SbF$K(Y71?ltkOLM_C@95vS55?6@9HMYJf^bD&9U!Z3XBv(cuh*xzKT^;m$3#+TV&I$gZA``WzHd1z?m?U(*$M)}uoi;aP2{QTq_pTy=l&PRluqG>5wfOW'
    'b29_^7dj3mYd>iUs^JPYXJ5v8040?CHy~R#BwV%Fw*BR}VO&%&K&k%sZHDJ%ye-kIl{%nQxR$8fgp`5iv%)ki!C2(hnm$VsNe)Bf1zHapde}`!)W^h+waVxyM+}tUQrV(cUtWEM-I~NsFl8GcA>D7dok6tIB_0B901<8`6;_}8+=)s^tRflpBXo;W&6;?o6zk<*>%N#}6!|tHW@6&HDM27sE8xVfYzmY=spAVxpl&**K{qi4'
    'o)0n4mGs$I?{F3OjU~lB-ui#Usco7vP9KsQL@-Qs!EywKKt6S)vW4s^AuFf~RJvR#`~oB)tqUlZwucY<Y@G`^;P%%#XFOqeqb2AM{WfgehY9L<{GX#q>l3A?X(<8i_sdfj5@F(@q-QijA^Tqn6)W`9f$7?S3{{}2%;O6=1KzBdzDAvPoz3t&1IXloU&nB$hCE{de=u*IkBmA|%NE!xPf%60Z{oxT$;QD(m;)tE@0UGFxdE_U'
    '67nB=yAX<_Lw=lI@DAJTguXeV6WOgddP5=%-D#0X*TiUM*vA^xV>`*uCimB1MPumi<w@vVYHKxC7rfUW|NLSf7rIJ!;+3IHR*a9N%9EY>ELbYEK3$@`$a1)Eb(@zOl%2+@ynYG6(U)=t1DWjjZ*X1l<;Q~p3xYUk&B(yUl5W9j?9IU`^}<sWwa^eJ|6Mh>9(LNM;rnohe;Qey4eM0ujrMno;A){2tYh`Opv_a8F6qxwG2P}5'
    'f`|NU3{EZK4eTle7^sJvBiFnLzJ(-2(;+-B|Myf$_hbdX2df34E}m;>)=_v9aV6)?kcuE#9!r6M_<I7BV1P9s_cd|d5BkC9w&^DH_6j04A9NwqS~m)&qz?5}fLq2*@WW}XPu&fuq{U_%UaZ+9HzgnnIW$}8Km8Y=dYl2(xTKl9H_+Ov+*YiJc^fu|r@@{Gkp`EZH*m3OTUN3gGD-%W8sc9YJHzVYKf|2h0@d|oc{S(#lq|&m'
    '9+-oI#l3|BJu70#`<GU1%AQByUkl8N@4ho#{L;H#zdx9{Yqo^gZhCicnTdTjVfWUtNhE^{{v6f+{oz#F3z$6QyS6g1ZkW&Q?<4ae$mk&o3L*}J-0E+OS3bW^g#n&acL1^<vcu%Pfgd8PrJTI`CSEsQ8_K}$`|nO|^Xg(GctFxR6A;`V4B^pZp4nM*<tO=MSE#|nHv2#A^Kc;~K}cr0_FzF=K2c5E>;w;&$|A0co^(OUN$z7-'
    'y<J|?)>?2%BIp*2fSxXHDabBY#X8N2I$z8YU%d@)fhC?ywxum6GOEg-9u){pVWd25q`{&6t{bZ&)P`I~^%O;%!9vQs;AGWb8;eAVp2s1{Wd>TzV-_Q)aPbRC=Cb-#Q|3Jh&CNT;w7v)-BhP*f^EOdnPRBvXQcV<(<x8jSb_iVds#14bA0N^KK-!+edC2h+Wp%&jb7J=)wvj1jW0ziEqZWV<HIwU&9xDy1<{<K8^0wS@@%^kS'
    '{BI6@1q}1Oj3hZ&_y?RxJ<`(wz0AlqxpdSag^Na&wFk5Ysb52Q7;i$^jj9neCK5ubuSw$REnNI~V!Eunny)H5XUi0cc^|mJo)dLlY&~_F#lJV<T1?ZMnu!>&4%9VmfAPvnOi>x`+|5JJN4E(ks`E9r0ZWnxp|Tu5NXX^X@`NXSxFTAtnC!Rvy9q=jr5}jc1H^KyB-(tT1jS<CuvoGJHx)%x^)h|`xDZjXj*yK0+RzP%8u%L`'
    ')E$j}h6#Zcv2w}jY=A+JNg?xKAl~=l>DP&d2)nv<yCwN1insJH;*De>c`}R)zmRRk8ZfGQ69aDb9ETZR4FX3)WKcxXhY;;+08%Y(ZON}xs?zc<1<(HomQxYt;0UwE`M_36xp(w$FY$`<#DG!atHL$|r=Lf)=hB@FQ3V2lm2~%`v8*MxfcE|v80(^b)`xoMMd9>m1Xb(S=k0_NlNLVF^6VA<xR%7g8lt4iE^n$4BaM8+j(wn$'
    'Z0`8t4t3DYn+l^gm&59kCuatVv?>B;9;yhXXw5X^u|u2s@KZ^!o#DKPESc?Pz~aI((ZZrlH0gv^K$SR+cP3P{F52&_2+qnC&fgbmrrB!Y<?!*yuQx%`@N$0;W5E!-N>)9I>PJ2#KW;R#K5gMbe;R>_Xg0Q>V?UrlRw`~6XbNw?kK|lE;mPh54{2~(w-qxq>OMDHIQm#@ftIiL&M`UmtTJkRn_00>92MFx`XUZbkpBRen)kA$'
    '-ZIM@&;Umf?Yf?13}GD4miHWE(q@m1aG<6~IW=qp@GIbA(Qv^Q6%9=7YC8c^%!o@cKp8~Yx(EUK%1kKFvrepXSx_)t$HVJL2SLSvt8_N!9fNJjNXeS;sTg3ayV3|yLu%kR&t5tk;}j=b1cwD7zG(Q;&ZEG~PMt+|WS%9_nnzzFK;fIEZkFd5mQoTl&~TkG`Z3d~D|PY*PA8W!R9SXksN<ge5%2a<0gxrBCm<Yr*zbPB1=SiD'
    'aNYl%gs-a<Y8`th?0$itpfE8dTDVF2T6oaCSWici!!!q^KF$iSfN@hW%w8Tn261L;GA&)Ia0a2!jEG%Ez;~z-R-Lu`43^sGNkfT>9Ns$jj)xON(r)wTb*cHG(G{z^$Jiy=-zV=z1}@tKWd9{O?jtKM>U*y&a_{|hF;ma{e4-fdtliHc7~dd1M!<FE=e+W8yS`mFOw|zIw%FfCQyvUC5c^3Rax2LC1<6+b;!lcA`)fY^FZ&L8'
    'Exiq<c}nqQ!}ia-GmT+{K}izjpzf>MT(=0UfGlpMznz+rHCm2<UkhFPwhoB(VFZJ-Y4_td=p}*^#7XT^ROg;CTeUu!%2{<e@;CXy9tDVE`4{@QmDwy@U_Qky!l(v401o{n+&KVio2=Gp{MPlTKYURL_07xt&~ZfkeAAliSpG*$GbmFS1T$SMINFia;9)KIECO@_;9s7q90%<M$>rtRHP#D*rS*#S7Vi%>+Q)&qjih1-ZUTFw'
    '!e?vYn93*lE{3VC;@H!r{a)NFI~bZxMM#|R3FdZZo)GFq8}PL&72=eT<I2+Fn8d2}<><7l#vpc%pX<dPvFE)k;N1iiDZd!Y-WFHd2w?H!I1on=+{ex&twVqsx~ppjz%$W~Dfl`&?^3=KFX1d=tccYI83laA^?~1~v6HL7_s1w%N{z_+4$*XaEHXH<xOtJ-Z{J2rdNd@CYb5d6UnIF&ORemCHD1%H-tl4O+Htrh32XRN+PF6n'
    '&EM-YwkHRvFx#3DO2!pzpyNbSu!*0ummUt(f_Hx*LuW2^7o_jU^A5CAyjfj8*pnO_l50N=UidtU6bW&>FMj(}RJFAd-@@B}?pM^~*E~m3$hS1oTieQgXUk328X&<Hu|jlqhsTz|)9p$V+<Is^z4B85Wz>qshB6#}h{HDJ$wHZ8NdRWoLPt}?T+wA$lv+VHB&*vv;r)y*p0B+Ko>C2OoBJCq3K@qRvnk`q22f|l8FB-_JrI+O'
    '$~_G)RshTu>0>5KR+eNOXV+%jiR8H?Y>jChXf<wJf@{^MpcurZTls}8)*faFn&Ei>yAyIjEvPSE-qF*fyeXp9wzwp}cx~F=E#`mviCJ#R2b~u+Bt8{2B`@`{%6(^mRQCz{NV8eng9xIwtXVzLL6Kiak#GTJGInX;$o6lyW7MbS%iW0nklK@bWE`Ef{ZFDcaS%=B(@cCdJS{y!_^Cjwp4RlX0xXGQNgM=iPnsk|@A^3Mn_E2y'
    '{-&ORyr*o?rN7#vD+CUp?ttH;gl|i#<pf*6<vt;Ic|Attlv(3)nMa|{EH6H16+u>(Pl?d59-Grz{v=BlYj~+po=#6(V2(&PYdQf)GYDbm4(vfgBVY>q&Smpwn*(e*9>)uv!22tW^GD#}u8tEF7Xgv(OpU8#MU}Quz4&yj{I=Dk?;)OfkqloMGm#VCzfFCzS{eWACsZ-M>L6U0ZE6f6U69_qN9}PW!+NQwFm-c_XxU>SYcc+U'
    'a6zbHn%igSI+QD3$~BTVmg6li>s##0hK`)+qlf(jKX9w{|5{nW1zIU@i3_7cP1p-o3fMae7K<4H!jF=gAE&G(SC7}o9d%2bZyP%k$uc#$a30+3crpYplnB4`-NI3q(Qu}b^NFO|6P?rdE@;TEpQiR{bu>zn%s1$2#k8mvV*H-0OA0j*WR6^LKadIaE%}&S`XZX3q$g;aDHmXT^PObxVW<taeYotTh-k(hV!H)nvzf9K*DjIr'
    '0`3bcWF}?@?Z>Yi?d8nhJL^Q}4%(P!e_(5wkB84WUNnwyg!h5iltk2z*SZ6}1he3=tm;AZ4_j0b$qL(|YWDY8n1LPpl)(*fEX5{>>(L`){F1afV{Q@hQ7Z)Ut5HQYL?yI^bbuP#$=E8VDxbJv8+n%phbMDth75f&e4qEgs_e1hlr^XyT5oPh_Wvh&Osr7832PoR;hc@+F`xZSi8ck%J)^2uf@)*j-~&kvnkn}*Ik!z#pVV|6'
    'yK@!d^S6n4-eBEOY0#O#{3G?SfG7_Iv}l5lTkNibZ6R8NQmW7&GqsiIVO|B(6eY=<THoWKOVWITzu=tH{+kbP;GuxR$mYwEo1zDl5`VA1Rw4W9=r->jMcHfT;_d@Cgr?1>v6pB6X$JGFlCVg^XYr{xbH<6{#hFX@UWsI`#gM}g@{zxZ!{|xTwRYVushIn4Em_}9xAPTwKQ92|N}dq3_VQZV$laKyRZ0qdbo`ODZ8>oi2-ry%'
    'mIUtIOOL~aR==DQ&uN6lAYe)|8PmX)MOYfZ$5-<Z3h^S#hC1w-F$_6gx2HWg*=0hr`Oj4cj%6*8PmO|Zh;A2`OgSG$17L1C!1c|@t6?O6JuvF~I%o4sg(Z(P^7Vk0c|YQz6bt}g4t$fR>KKt@s4$cSNmd_QxZlMcp&1UEFp~xZRf=URd@S{NpP^Td$*BoN_U;K|nItTDB)Cj?Q)NZ><W6<7r<MFxz@x_64hrk3R^d#C!&~ds'
    '+=<n~ao@e<U*AYq{!7^T5$7}$m(b5+eGqP}$es>_%0$r1;Ssz&3M)fYx8V)FAPoL3rk8iqHwyY?p=RxO6|;x?N3R>P@;Bm0$Fkyf<y`003hTAkAGfrl&Dz-HDHfYW(lM&se_o_07$LcW@-_vyy)lPg^Bqdn$Xcu=n*oz{rA3W_QmXlg72_Skv-rHY(Pq}MG`_)6_~6t4xb5RXj96AS8as(nrQfJd*SWSEx~k-ke%M5sXL(>T'
    'hCBe32TAoD&1ErML?SxB-<Z&j&l)_;TlhA0-?g7{>_v{M&!3p7+}|foG}1rxNwM!eMq4c`4<Zb_mDs275M!hQt)!FM=WYwUsh0Ymy4a+0n;$zTK;7Hl%V23DDcLtN6vr6&1<_F2WKDw$S53PVZM?!Go&fkLXt3Z|!S9L_QAzzQ<_wUA-cFNUGIi$53rlc9j3n!V?}!=YwOf21r^}+4Lql6vrdHPw1f&$6KF+WOxoQBn0ylv!'
    'mjs!KY)AqDfV3+f|5Cak1f>Ro$T#OWbM&Y9gFv5O;6eX#zwZt_+u+n5;yq7y)^Yfp9ED**_jE+uvDJaZVqmVWBphJBh`S5Aeo%<<V{j@&P@Cw6;z!g#xOVsP$rnN9;M;D!1ftdud4dAvc^T>Dhatta<&*VLnIdD>%QMSf%-@pT;g?lt*>+_j%Z^X0|F|=tue};EaQA%{`w_0z*0j9ea`LzckZ44%ks{FN`wGBe9@T!BN19ZR'
    '<+=Fln?`BYT@=CdA1lc1GyP?r^zqBZkz}ZN-ChRiNARd7F&am7)=@Q5^lN2R?L-ya+blY9j>B&b<89pm>Ev0dfibB6q{OLnCX|u!X-$^Him(q0@I`P0uvuZP->j5$P=!wru;7b3hyIa<g99f=B5K(vzT7%Y+l)U`7qQb)Bb%%N48(F9q97A;6a#+6x95maDVor?+%ztS3^~*uhq}!QH4sSQSs`*Db7RA5M9;~u1iH<JWQ%30'
    '_w8vo38h0$O#>9k$5g);<0Qb}&35=Qgt9OCDIclmlDX|GTwT~8uz)-8g&_S%l#gYqgn<R}2X(wKnDb6I-+OI#lA>Ddl?%Z0VwdD_n_YE2Zw1yH{P)5FOTYma0Z9-h)Aq2S?|>gEAGkcN6r`p<|Cn`?@B*-1iZ@5I3RJGgMh0znqo=&*^}2JjXV%R=a;>yCg`w~E?T!V}frbJ*in*I~!Hpyv2bFA%u)b4Rq&i2ur3bl<r~{ta'
    'QVx>ghu@X<UYFp|#*oB6cn(?-^3>#r-N<8km;?qGDuL>FPMZ)rGNe;1m54uEAsidv3^__eg#Bn&Id=&DLkWsR!qQIRZ>BNgiAr?kun2$1@#ep5MdHjs+zX8ON?!q135mu&yexS)>5l^WclY1UlD@Rb`&46zRz6d6sX^4dErKwNf~=7tZn3WIs&=xPKwPu~p%E&4z7kq3(4&fIpv&dUz6Dbv>nSoSRd!U(gL4$ibHlSJ0~VH6'
    'Z^JF_kgAgZubqH%q?Os^($4MOb-`h@sqAw0<AH8#k4ZlS?p@5g2u<5{YK6^?De^RTf!h_>b>WZ9hf&TH90uW`I$1Cw>_jq!L`0oc3C1njsQj_r5dk}p<p~XX&_3Elfm1O-uE!{A<%FmOem~i5{}KGuSuoC=*i`2I)>=k-1k6Yr#AC91FmhI0(}lS^E&+`#z=hWnOZu}Ug1<%nju%MpX#gALj<LmgzQ$N3Y$`(JnCM-K^J@Db'
    'LWF7Et+!GA*%wV@>W6bAW$X+NIO{EwumCH`b&toY+#SquP)SZBxtWqKuj6Eb7e}g?3>Y++(>$}dA-E5Ci2n*W4RS~%Jn1-j78PinwS`;WKDcUr#lnqFSe}09{fgle#$zIiAE65)<QXo)dJh1!)b9r=BTP&n|7G${j4LGb&AN@MlCh_HAA7ipi>GN5k`~kjm)zRbcRT?a8-aXc&%}hV;9XZtBz1{IXy5pIt-noiDV0}@c9np('
    'C^^kM!)$$R3B4St4U(J34cOr&7VVFXApA4q3#ifGP|S_Y%bHZ;lVrDVAl`e<q)p*|gn3}r|I8tOGmIufun3sRHcg_LZLZHGa`8cD`-Z_J!#hQx)@OSMADqIJXKG@C1S5#wU)}D-1C3?60>Vmvy{v$Yk!W*yjY;srry+j|1}qCWLqCRcMu#HN3{Biu%KprBXlbtE{EeQG^b8dw9O3#RjrRsgjr*xonu`Tm7<h>o*FS_)DQZVA'
    'Gb?0a{nq%BGXS6jdD84Zgq?}?`2C8>A_LRW_TFpd60Ux-s0~CKqo9fhPCUc6O_M0sXawfLb-(rTD|^_yJ&bW;CgvG7E<7O0u~&+{wA_>Q@U%Ef?X4MJuvn1exK2~FGOue!i?D%r?H^qX?X(7r`hEvim%ozW6*ZS8B-K(wB)}$WP900hFI!7y1)4B=SJxSy@UU;Kw?SuuHrK2F1R1#Qh)%IQd0xT!+?X+f>d8##)Y?zE<^f2u'
    'H_OtrJre*(9ogQcgA7}0mCJL&poIPW_p>(}5Oo&XJxaYhzn-(yu9if3eh5n55{|TAyYXK0VlIl|_ws72A+!jy0~VtnnfV?Q9hps`8}0b0dF+ldjFg{dWNP*H&^z6VL~9?F!HfivVHMO|UaDvUx0%U-e7|e1u!@mgIk#hO+M=fL`05rz)a$-t9Hzjd9@WdV|7#-O?vh^0L6B?{R_?nRsmmG~ZPygmgBLE``LcU+=IIar_p9&t'
    '3oKyGfUG~qQyE_}?UbJ!X&O{62zPQ#vzdlt78DO6v!M}KoezC(@j?_Ac5tp!P3V)Nn~fBiDESSxnk!BhB`QB(HBfA^r{*roVZVr$?Sdh&PI-PE7s?kX3-#hofv@5IR}zdBeUFilwTosAQ9>-88i=^P+jG$OosWguNgR`R%9dK(?4YiVK8ocwJCO$wKhAUe&Ki6U@l;@*1#pO@mGW!ei*7eJ_Q02?=B2>I4$IN$YBdqYR!xS%'
    'A|C>^Nw(z-qJc5)>;<Kj>dqD%vPx@tAcUjE7yC@t;>BB{XoS42BWK%@ZSD9a;g4E&-y&n7>{_zL0!zMlcS4Hl9x^wfoTA110;!cF>5zTgbY;ONOr5@I9;!VX_oD%Uz-zEe$mVG!g^==>F{uYLRh&oYyu7GC<q%x!ptY`TrVc({y;rcL^3{|lZHY^~8Ox$p*d>SZ_m@+hcoBUmer>L+&u!?%)0MD>DkJL1y07%F_+{AAEcj}4'
    '_7@xnwRaN}&Fb6PRw_O+G?9ix57@m>Qkj#UR^F!IZhw;z5-pp+dpMobfl1hk^8YmTy+x-?JE5-mrJsaf=hgD__rF+0`p(A>^pB*fvzU|(_l|G5JRIhmV|De&QMPT9y#;T}()GXD=T(pbh7OhcDy0iL3Yu}N-gx#~^uCw0N{|B~wN$2tFrcim2C5vv{T>Fpf5x3&yasm27kct~i)j=q1~`W(RIfePc;N*bLm1kdHHAj3$&v)4'
    'iY&Js?eZ*C@z9yAhkfZUIGnrpdIgzirBg|jcbDdRNt}jgtcZm#hAFlndta^K<4VCmjsHz@+WT-5fz~Q5yht$~dO)8%)}Rr_yRQ+X2-*)bmW`XPnH$;gJKfRQt9nx1p_#s*IEugD38KgwBJgR8SM?oACs0<9+4#y)I5qvNgDn}z?jj&mDnVn(61xOmG)AfjD)+rv7`tZ)ipcuznz!Z?1jR^@I@!w0;ZRE>^oI_Fe<`K(KOjVQ'
    'zip7}`uEy|Zco@xwu@RIcydsIi16N30+7Qw|MNw&vbwK{@->P~mVx&p3CzWx0)Uj-_vRBz#<;`~qX7h4V)s>Yo+ulJwb$Fb(^03bI5hE_XKrzxt)DM!d&&}H$b(F0+VzZoezq~;-Al}VlQlRA-*1xdwM!&4d=9*680I$CsIb`_Nv5!e5ih2OCOk2UMb>_r%be<Ysy{Q5B^)ko%?{C<Y@0BSXH<cN*yr+4ns#zg0f~jJm?H{9'
    'ZPe4SukobM6Rn72eUgL9mYxxQ%t3%E=HNa9g>#=&?%ZE%S;=12AvX|;<_d23aVYTTbG{q2`DlxkQas8ItppfdH?w7i2x*#&?W!6d!fR^zKFz{~Sa01&pWtB*rki%#rkh6=A7*Oz&K^cZ_3_U`7l(+$zY{bX4hZl^3>3u>HI+wqkKv(Brc^ENX2s1saLcvaW1|?NiWShTg*dtQ-?&=0!?^;RJ_oO|Hz)La1DPw{S<=qkf_z7k'
    't+VC<R*LHd50;F+Nr>W`-Katu1Qs1l$2Fr&N7~)K#;Q>m?7Md<d<Q`RvE-I&TWqi3S~CH_h@j#33ri&?C8J&j*MMl2Mf0WoGJfQRw!a)}^_y6vcPk!7r_D;(_N0Mh$8y`Bd>h+>JLj&m2<A7qc&6rNFjjtElt}09;GS#ZbwYqCaFD4DN!xu4Fh&F}+I#Ar(|CJzctR<)WQKTrvje!7-b'
)
PREFIX_RAW_SHA256 = 'ab9c62db4cc27885c70acda79bb8c4427585206010d7806c8e368d80736c751d'


PREFIX_DEPTH = 30
DEPTH = 60

def transfer():
    return Transfer(225, 225, (117, 270, 4), (126, 286, 4), (144, 356, 4), 580, 710, -14116)

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
    regenerated = coarsen_prefix(*conv.evolve(PREFIX_DEPTH))
    if decoded != regenerated:
        raise ValueError('packed prefix must equal every recomputed tail')
    moments = open_module('moments.py')
    gain = open_module('gain.py')
    energy = open_module('energy.py')
    if energy['THRESHOLD'] != gain['T'] or gain['T'] != F(-14116, 40000):
        raise ValueError('threshold agreement')
    rows, initial_energy = moments['compute']()
    cost = energy['energy'](rows, initial_energy, moments, gain)
    initial = energy['initial'](rows, moments, gain)
    if F(cost['energy_upper_1e12'], 10**12) >= F(7695778873, 62500000000):
        raise ValueError('strict energy ceiling')
    return {'directed_energy': cost, 'directed_initial_majority': initial,
            'prefix_sha256': PREFIX_RAW_SHA256, 'prefix_depth': PREFIX_DEPTH,
            'recomputed_pairs': conv.msize, 'arithmetic_backends': sorted(ARITHMETIC_BACKENDS)}

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
    if F(entropy['coefficient']) != F(128754752375347193833379, 25280633597805000000000) or entropy['checks'] != 1097:
        raise ValueError('joint score or dual accounting')
    return {
        'threshold': '-3529/10000', 'graph_cap': '7695778873/62500000000', 'coefficient': entropy['coefficient'],
        'terminal_upper_num': upper, 'terminal_upper_den': Q, 'depth': DEPTH,
        'prefix_sha256': PREFIX_RAW_SHA256, 'prefix_depth': PREFIX_DEPTH,
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
