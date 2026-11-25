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
from time import perf_counter,process_time
from resource import getrusage,RUSAGE_SELF
TRACE_WALL,TRACE_CPU=perf_counter(),process_time()
def trace(stage):
 print(f'circuitbreaker stage={stage} wall={perf_counter()-TRACE_WALL:.3f}s cpu={process_time()-TRACE_CPU:.3f}s backend={sorted(ARITHMETIC_BACKENDS)} maxrss={getrusage(RUSAGE_SELF).ru_maxrss} platform={sys.platform}',file=sys.stderr,flush=True)
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
 need(s.d==s.yd==225,'root coordinate scale')
 n,u,b=s.apply(s.A,m,p,m,p,plus=True)
 n=[Q-x for x in n]
 u=[x if x<y else y for x,y in zip(u,n)]
 b=[x if x<y else y for x,y in zip(b,u)]
 C=Gaussian(2,26);S0=1<<26;wn=4*40000//hh;un=4*40000//hu;G=40000;K=13
 def mass(x,y):return max(0,C(y,G)[0]-C(x,G)[1])
 xs=sorted(set([i*hh for i in range(-wn,wn+1)]+[s.T]));wr=[mass(a,b) for a,b in zip(xs,xs[1:])];wu=[2*mass(i*hu,(i+1)*hu) for i in range(un)];mr=S0-sum(wr);mu=S0-sum(wu);need(min(wr+wu+[mr,mu])>=0,'root mass')
 aw=len(s.aux);vmax=un*hu
 dl=ceil(21*(xs[0]-vmax),15200);dh=ceil(21*(xs[-1]+vmax),15200)
 columns=[bisect_left(s.aux,x) for x in range(dl,dh+1)]
 columns=[(i,s.aux[i] if i<aw else 0) for i in columns]
 yl=ceil(3*(97*xs[0]-41*vmax),60800)-s.aux[-1]
 yh=ceil(3*(97*xs[-1]+41*vmax),60800)-s.aux[0]
 row_offsets=[aw*bisect_left(s.y_axis,y) for y in range(yl,yh+1)]
 n.extend([Q]*aw);u.extend([Q]*aw);b.extend([Q]*aw)
 def line(R,V):
  ix,i=columns[ceil(21*(R+V),15200)-dl]
  if ix==aw:return Q,Q,Q
  z=row_offsets[ceil(3*(97*R+41*V),60800)-i-yl]+ix
  return b[z],u[z],n[z]
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
    A = F(1433, 5000)
    eta1 = F(384303974411, 1000000000000)
    eta2 = F(67917029353, 250000000000)
    beta = F(39458297, 1000000000)
    theta = F(14651591681, 1000000000000)
    lam = F(0, 1)
    xi = F(-1510850871, 200000000000)
    kappa = F(70670083, 1000000000)
    margin = F(1, 1000000)
    zeta = F(32649399819, 500000000000)
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
        return (low, low + 2 * power / (25 * (1 - squared)))
    ln2_lo, ln2_hi = ln_unit(F(2))

    @cache
    def log2_bounds(p):
        if p <= 0:
            raise ValueError('log domain')
        power = p.numerator.bit_length() - p.denominator.bit_length()
        scaled = p / F(2) ** power
        if scaled < 1:
            power -= 1
            scaled *= 2
        lo, hi = ln_unit(scaled)
        return (power + lo / ln2_hi, power + hi / ln2_lo)

    def exact_bound(const, terms):
        bound = const
        for c, arg in terms:
            lo, hi = log2_bounds(arg)
            bound += c * (lo if c >= 0 else hi)
        return bound

    def exact_test(name, const, terms=()):
        bound = exact_bound(const, terms)
        if bound <= margin:
            raise ValueError(f'{name}: insufficient exact supporting margin {float(bound):.12f}')
        checked.append((name, bound))
    p_anchors = [[211498, 211498, 211498, 365506], [217751, 217751, 217751, 346747], [220659, 220659, 220659, 338023], [172966, 258638, 172966, 395430], [177122, 267165, 177122, 378591], [176551, 275054, 176551, 371844], [158972, 207369, 207369, 426290], [167179, 214487, 214487, 403847], [167179, 214487, 214487, 403847]]
    for b, ia in product(range(3), range(3)):
        alpha = (F(0), F(1, 2), F(1))[ia]
        o = int(alpha > 0)
        cu, cv = ((eta1 if b == 1 else eta2) if b else F(1), eta2 if b == 2 else F(1))
        cz = 1 - (F(0), eta2, eta1)[ia]
        du, dv = ((F(1, 2), F(1, 2)), (F(0), F(1)), (eta2, eta2))[b]
        k = F(1)
        flags = xi * (b - o) + zeta * (int(b > 0) - alpha)
        p = [F(n, 10 ** 6) for n in p_anchors[3 * b + ia]]
        if sum(p) != 1 or min(p) <= 0:
            raise ValueError('invalid four-corner anchor')
        U = [p[0] + p[2], p[1] + p[3]]
        V = [p[0] + p[1], p[2] + p[3]]
        X = [p[0] + p[3], p[1] + p[2]]
        Z = [1 - p[3], p[3]]
        for index in [] if b == 2 else range(4):
            u, v = (index % 2, index // 2)
            cost = (1 + u * v, u, 0)[b]
            terms = [(cz, Z[u * v]), (cu - du, U[u]), (cv - dv, V[v]), (du, p[index] / V[v]), (dv, p[index] / U[u])]
            query_weight = lam + kappa * (1 - alpha)
            if b == 0:
                terms.extend(((query_weight / 3, r) for r in (U[u], V[v], X[u ^ v])))
                if o:
                    for a in range(3):
                        event = p[3] if index == 3 else p[a] if index == a else 1 - p[3] - p[a]
                        terms.append((kappa * alpha / 3, event))
            elif b == 1:
                terms.append((lam, U[u]))
            else:
                terms.extend(((lam / 2, U[u]), (lam / 2, V[v])))
            if b and o:
                for flag, marginal in ((u, U), (v, V))[:b]:
                    probability = F(1) if not flag else Z[1] / marginal[1] if u * v else 1 - Z[1] / marginal[1]
                    terms.append((kappa * alpha / b, probability))
            exact_test(f'rank0-{b}-{ia}-corner{index}', cz + cu + cv - theta * cost - flags, terms)
        for name, credit, hmax in (('input-U', cu * A, 1 - A), ('input-V', cv * A, 1 - A), ('output', cz * A, 1), ('info', k * A, 1)):
            constant = credit - theta * (2 - b) - flags - lam * hmax
            if not b:
                exact_test(f'rank0-{b}-{ia}-clip-{name}', constant - kappa * (1 - alpha) * hmax, [(-kappa * alpha, F(3))])
            else:
                exact_test(f'rank0-{b}-{ia}-clip-{name}', constant - kappa * alpha)
    edge = F(19, 100)
    exact_test('rank0-2-0-outside', 1 - A, [(edge, edge), (1 - edge, 1 - edge)])

    def vertices(a, b, c, d):
        uv = {(u, v) for u in (a, b) for v in (c, d)}
        for u in (a, b):
            for v in (u, 1 - u):
                if c <= v <= d:
                    uv.add((u, v))
        for v in (c, d):
            for u in (v, 1 - v):
                if a <= u <= b:
                    uv.add((u, v))
        return [(u, v, z) for u, v in sorted(uv) for z in sorted({max(F(0), u + v - 1), min(u, v)})]

    @cache
    def roughlog(p):
        n,d=p.numerator,p.denominator
        e=n.bit_length()-d.bit_length()
        if e>=0:
            d<<=e
        else:
            n<<=-e
        if n<d:
            n*=2
            e-=1
        w=((n-d)<<30)//(n+d)
        square=(w*w)>>30
        power=w
        total=0
        for j in range(8):
            total+=power//(2*j+1)
            power=(power*square)>>30
        return (e<<30)+((2*total<<30)//744261118)

    def cell(a, b, c, d, depth, alpha):
        o = int(alpha > 0)
        cz = 1 - (eta1 if alpha == 1 else eta2 if alpha else 0)
        weight = 1 - eta2 - lam / 2 + kappa * alpha / 2
        if weight < 0:
            raise ValueError('entropy chord sign')
        u, v = ((a + b) / 2, (c + d) / 2)
        bottom, top = (max(F(0), u + v - 1), min(u, v))
        def slope(z):
            return F(roughlog((u-z)*(v-z)/((1-u-v+z)*z))+cz*roughlog((1-z)/z)+(kappa*alpha/2)*roughlog((u-z)*(v-z)/(z*z)),2**30)
        z=(bottom+top)/2
        for _ in range(3):
            s=kappa*alpha/2
            curvature=(1+s)*(1/(u-z)+1/(v-z))+1/(1-u-v+z)+(1+2*s+cz)/z+cz/(1-z)
            trial=F(round((z+F(693147,1000000)*slope(z)/curvature)*2**24),2**24)
            z=trial if bottom<trial<top else z
        ref = (1 - u - v + z, u - z, v - z, z)
        if min(ref) <= 0:
            raise ValueError('cell reference')
        checks = []
        for U, V, Z in vertices(a, b, c, d):
            p = (1 - U - V + Z, U - Z, V - Z, Z)
            terms = list(zip(p, ref)) + [(cz * Z, z), (cz * (1 - Z), 1 - z)]
            if o:
                for x, X in ((u, U), (v, V)):
                    terms.extend(((kappa * alpha * q / 2, r) for q, r in zip((1 - X, X - Z, Z), (1 - x, x - z, z))))
            for x, y, t in ((a, b, U), (c, d, V)):
                for q, w in ((x, (y - t) / (y - x)), (y, (t - x) / (y - x))):
                    terms.extend(((-weight * w * q, q), (-weight * w * (1 - q), 1 - q)))
            checks.append(exact_bound(cz + 2 * eta2 - xi * (2 - o) - zeta * (1 - alpha), terms))
        if min(checks) > margin:
            checked.extend(((f'rank0-2-{alpha}-cell', x) for x in checks))
            return
        if depth == 10:
            raise ValueError('unresolved entropy cell')
        for aa, bb in ((a, u), (u, b)):
            for cc, dd in ((c, v), (v, d)):
                if bb <= cc or (aa, bb) == (cc, dd):
                    cell(aa, bb, cc, dd, depth + 1, alpha)
                elif aa >= dd:
                    continue
                else:
                    raise ValueError('dyadic symmetry coverage')
    for alpha in (F(0), F(1, 2), F(1)):
        cell(edge, 1 - edge, edge, 1 - edge, 0, alpha)
    q_anchors = [[599046, 870898], [580527, 940208], [569282, 889189], [697968, 999999], [674758, 999999], [669201, 999999]]
    for b, ia in product(range(2), range(3)):
        alpha = (F(0), F(1, 2), F(1))[ia]
        o = int(alpha > 0)
        cz = 1 - (F(0), eta2, eta1)[ia]
        w = eta1 if b else F(1)
        offset = beta + xi * (b - o) + zeta * (int(b > 0) - alpha)
        cmax = F(0) if b else F(1, 2)
        if b:
            jmax = alpha
        else:
            jmax = 1 + alpha / 2
        exact_test(f'rank1-{b}-{ia}-clip-output', cz * A - lam - theta * cmax - kappa * jmax - offset)
        for clip, ref in enumerate(q_anchors[3 * b + ia]):
            q0 = F(ref, 10 ** 6)
            weight = lam + (0 if clip else w)
            constant = cz + (w * A if clip else w)
            old_query = not alpha
            if not b:
                weight += kappa * (1 if old_query else alpha / 2)
                if not old_query:
                    constant -= kappa
            for q in (F(0), F(1)):
                C = q - 1 if b else q / 2
                terms = [(cz * q / 2, q0 / 2), (cz * (1 - q / 2), 1 - q0 / 2), (weight * q, q0), (weight * (1 - q), 1 - q0)]
                affine = constant - theta * C - offset - (kappa * alpha * q if b or old_query else 0)
                exact_test(f'rank1-{b}-{ia}-branch{clip}-end{q}', affine, terms)
    for ia, alpha in enumerate((F(0), F(1, 2), F(1))):
        o = int(alpha > 0)
        cz = 1 - (F(0), eta2, eta1)[ia]
        offset = 2 * beta - 3 * theta / 4 + lam + kappa * (1 + alpha / 2) - xi * o - zeta * alpha
        exact_test(f'rank2-{ia}-clip', cz * A - offset)
        exact_test(f'rank2-{ia}-unclipped', -cz - offset, [(3 * cz / 4, F(3))])
    saving = beta + lam + kappa + margin
    delta = saving / (3 + A)
    cap = F(123275637453, 1000000000000)
    coefficient = 1 + 1 / (2 * cap) + delta
    return {'checks': len(checked), 'local_reserve': str(margin), 'saving': str(saving), 'semantic_increment': str(delta), 'coefficient': str(coefficient)}


PACKED_PREFIX_FIRST = (
    '{Wp48S^xk9=GL@E0stWa8~^|S5YJf5<HiS-r(FP49sEF>O4u-D(l`DNH-ty9XSM)C(T3zBh1|QEW5xE7DuWu!*WYT#AgwF6Z+C%v_4|`ynBngn=Ose9={eF6=YqOJ=qAO#tZ5l{KIlt{$)-ig*i6(xg`H)z8N$-Y)7eYF2=F(4rG(4R;(1y*3g#V%t4Js2Q`TryPh+0(tt72;gyh`O*Hb$3ja=tgq@~J~{*SYa4`~SS`AOL+'
    'fUKdvmc3lcnVXAt(fHKCkevH+0%&!II5;;zDwo7l=cU4L#aQWjm+QJAZt}LF=U&JB&Z6U+;C43CU#%3VgxBNW;kjW@`pW)uryjYW9v@?7La=A$E@0)jZp`WTvF5LfX9eJ@Eq6Lc=Sd5SA;xMd0K)W=YhWN;!c{dErI1v>s%cNI74cM2rb19POFCC)(M^31rRun0&k!WYwK>XO=}oGMlC;DvBDZsSqN?wBDqWk(Tb$w5-gcWX'
    '>33Q?Vv?8V#gBp`=4)VTHOO69>83$^dUF-0DilirafqePR!lQ}Hhp(s)*|JN2z^mVD!OS~=p^zMxC7-D5Nx-~Hac1cgT1v`@fuSXPGS{eIip^b^5FpwLiTfopXWDa1b&tTOQ*{(bIZrPABL&i#tkc)8qXPgTh|$3z_eOA80Be^zNa26J9x~+3-F@wUvN>*J>Ez<d&uxX;){p4vZUg9QP&Z~{cUCf2>Ca&*{4D&K}qV-c}Iqr'
    '|65>$@ySF`EvA}laiJL`NzKY+Vm#bk)*AVYv<q12_(QSG=fIh=8p;|t7mW&w!siL5C@f#XRUg6_ys5r36xtPDeR7gL;}g#-fvd@FW-=v=pXhuD4*5h;7!3Sx5sLqWznXp3^e!~M`P-gq_N6)woPuacje##8$$McZtLdsj%x41;D-LpXW$B9Mc#}Y*8L@{H!(%L>W@ZhUp8bn0j{@|%VTme2C|gkwIi}+U{8sFN;FO6;5z;y1'
    '^@K1c{RECUyE@yL2Tz%8Al^!2D#KaQxPfbl5+@830?-w~ZMYW@Pukq(f>*>YIgTqc2xmm(NQLW3|JM+A*|7B^63ck?RSx<p`@_q1`_GuMsW>k>J`lT2G2gv+ck~9j*Z{3tv4fKMFwNTlHg2rWBbZ?!r7U8p29;J4Swk<vKe0H8h%e`{15q?CpL=*JNJ6pN2`8!O%01!BcgmW1ySS~FjSo0NEJm0!L#a)k{;93Y#MLo)g<djD'
    '7t&SMBH~T-$k7ib7~MvI=b^0x)oTvG;O0MaA@jw(Sh@hm>D&?G7)>SWeBr++Fq@J@#!D;tgrLmM3&8FjkJ9h_wm&AIO*Mr;4tWe>+`z)qz`OdAe5-Kh=nhd23;PU;!eW+v9FaD)@8t6ooirZy8-8jUa<}!mKHP`HdY-?}qt?Kci{S_y=mSN_@|fXYAy^#_yabQGnhnxv(h;m%8BVSE^}_LXjwh~`2#t_NNg(zqx1+!LTOW^8'
    '0ix7+_}6g<#in9$>40eLcz>Q!*FPxmGi(>0_^fz2WSXJPPL8^vtSVZUfMSqnWy|e2WB&61)1ImI1N(1dqyga%fcCC6)4qbHS?iN$qSeN-&*TMq7vd*BXzW}>Mil1>vvjnxuJ0%EgCo+~pb5SDH&outhn%7m_`>j@cNZdk(gq4yCkW9gQ30XVO><6lmV4{PzKs&W+&Vn(>IhV9pbQC-@ojS(waD5Yq4x9D=s@ezAXfx!Gy)(}'
    '0SFIn6r11IRxj6TPjr_{T6G;~;ACz<|6AXy*`em@;)$hsR+){f(IqKow8HmF*1WfsN@u@NYl$BkSI=O<joFilB*eT7f*jal1|PUi;~nyp*cR10V<ekBXaq_!VLcL6>GtsCWY<U-*z>xL6Ds*V+Ab38rOp8x6=|NOUU;ADGF~avW!TE^Mf@aYp!m5g@Y=^^9QX*ZDSyEQBxGnmCZ9c1w&dx@8U0m$a1N`1jx#WE6*-pj4l%g2'
    'WkY_#86(HPfQC!8K)I0m!5r9#xYckyfloX_(kP@ortW`|Cpm5b`zz9wAyq1J(AB8-2hY#w4to^fCmmPB?-Tnz%C2R%MbN!yEahFC_W}0ZBP4ovTJ-rV6s5@1M#<3cK1Ag6{+)n-Yxo<zq!md|QOFwBd0T8ElF8vGyfoyzs70A$g*)6ENLr^VY#X3jC10#25D(w36n(mbMD%{O2FKtu9xO|}gaKw_tKWqb41aeZTKa*opdquG'
    'rwgjJHnA<fU#T;0YLc1wY4`TLas_D|GqK1>6Ft$s6iU#l()#Zh4OrSugL{eU+JfmXpH)&|WB}?b5mhylA5Dc5MNpi8YB4)(XI~P38bD#NoVHRN_$iRxLWig~g?O$Xe}kcO$<yy$ZjY%<{Xp{|_`v+U6c<zDT-dTgZ}8?&T?QkAoFvi})(lNr7nmV1yCAX(*>Z=kOS`oDVeqlMBTQY+fPA4j?&D-Y&$?DE7KUFW$(`x$-sH82'
    'Ht_sGppJN%zr)PPtK|q__lS9WRICdIB2wm%??EyUH9mFdk6R#YJW3$dl6Uk&zp$_G0G0dqO#e9kei)O8Re7v=zh(P)0l3s_cN9g$Rb9cwjdPd_c&Fi|%XQOs5$HAabf9S}24~v$*Cut^YT9fLs;r7tX}up7Z+1Ir?Mta{=J8aykRgoXQ%fMxvF<HZ1cY{C0w|!@4d-Kn%~g>;UEVe^!?8~}W|Yjw+Uop*N*jG;J*R}Ww6@nb'
    'W~>V_)z^+pk<H%rN=PNmH@{;^NPqUcjO3Ro4S92jZb?Vv=8oBUQ>Eqrl`VQJMmX;FoBmhx$pE!QMfc<_+4fEP3nd&~hJt?o%kN*R8c%44KBRcIu2@>mivK&<2wON_b}NE9s_%^4mwOKfka2Y;JSRp*LGQAILMbimEr28%CYE7Id)-lM!aE1l9ge#|>G@2AbwyuZK(FKRWkYk#&4Ne(s;W1B-}y{C2AN8$-5*>r&TV*yMf_e6'
    'M~~O%xFFP3ylHKX<I}#`p$xcJEr@(v_$-c7@s`GUOgs0PxeE#T06PWr!-PB@wJ>)3FwrsPSCmBli(fXSJN2|GKq*N?`%VvfVc$2l+6uDFkF?u(D|M&tDhR=aZ6+@=tc3g4l0Y9(JC%#Ms@Jg%Zg+J_>dld++HTb39v$8H-Nt>uoG?C~@0XM`mvMo!pH;^*%o+;ePni9@8Afv&FP*H>qnU?o)yRZ~AI~x{b_R%C4!<DX<7nNM'
    '7!{8_edMR$nc5cT>xY0e8r9IDI_|027qk$jtSPQUJl+d2h_1kspVZ}0NUTtC$``#-Z40_lQPm6iTt#0dc}hSr7?$8PpzD@0#kvyHvh}f~BuLE!>(+jYmHujt{Ero`hma~7?0=4*ucHb{<F!eCM$|Mk$d)#2d*4vPatO_)cF)JGM>q}VMkop3E^Xj8sx*;4$1xRBwNhb@pjQoy&sw=~b$`e&<la~TDPj+=2-XTO^CrVKwoLVN'
    '_<wn1UCbc+SL~HD!GGxz+MWcfR89P~FcNbbp?qAzHGIp<Ec6HOVcTbxT*hEIky4taon`PidYOEkok`(q%K=c%f+P}6k=cn_Z~e*7O{*5B+j6XuVns4#yM};>g9Yvpm#SRm#0f$no5NwSW7vlDH%98TP@RuC*~$YUD}f8^|I2Dhv4OUYjneG0W^PZ-rYg%w&vg!yF?lLcwA7Owbbkz~pI9La$l;v5(i|yxY@~`skc}Rb&Qe)y'
    '0B6nP-C4<Xn8|p4u7n;NXUEtkU)Z=Q_BXGB+h<jC8Ttq|60gZ6DSK5Iu=#jbW1}nALJma<^?iRXY{AXN=2+JubQLgko5>pU;^ub&XaJ|!FO{--y`cy1X)*wVu(Ua7ZFxg<cZAsuk~|Te5_<h&V@AUXbfrPxHNRy8%LDd)4k8-8BDo{ywD#NW4IqFPw7jSB4!9h|p-qotiABchXHean!#($JNHGRAE7c(?7kHj9O95>-CsX_='
    'AISayh5e6Rb3HEoq_D$uy@!Sn*Pg8Jj5Hbj>Tvdp{jaG$Ll0>iDV1a2AJr5{LkOZ*Y#ULeH;!H&HS9}TNbvLMF+egCq>w=7XMw;m-()Gba%no$Wu(x=k7HWX*e`Cun>qIpY0zO$))(3O#_>I04%$!YANV#LylyWMRA7k~pPjm%kE^{iU|>r)2$5*`d6P!?UDACM5&KI!d5P@YQ0&y}AG%>YKMb?4C55nM>8P*}84={fM&A-#'
    'v+`gS(ARyR`4f{bKzB8=qtO?A7duNqp&fY(?5>C$N9b|L7J}35_q}je|KYME1`^7Y={`km?DkQOUD*3|tY7VZ{n!G<37vX_f0T$W9_eoj+u{jM*5+=v+A74GB(av&4nA?_EciVG0#0BhP20&#(e*xrLTI-2HtrbGRJUnzq1+0?y74G_Np$;a>$z*VW$IFj$1)n!Rs1`qklS${E$OmXOwjz5pms=ydu8)Rj0CZr+bQIY6?{uW'
    'uyJMXQ0Xs$mt*KRVZC`u!Wvvs0xN$x$Q0h(Qd@9`ex5x2Fkh94lLMw?XqgE;^6RpMar=?!I6d^a(8nD_Piqbl>Y1hnYWEo5Ra>PKIX~>V%m%1uTH>kP+HO_~c%MH*vB-wd8)DP@6z2=Z>|f)2i4|N!CEgSLreI6vcIx4r%bb{mnxkL8q7|a&J_gjf266W3;=!58N9Rsv7m)*k2@54RT)|0J*D%Jo^kruD6MXTK<}eT73B^<$'
    ')wKK`^dpawS0^q_{&&|p42UEzW!7S@N;v?;fOs#4&-s%RZn4IE8i@Dd^2C0>q6#vJWU3={;%0i4rhYZCp7f+NS=&}+LKA~qF{tg|$y@=5?h5hCvTMmmxmyDZ-boP@VZ{GYhEr5(bVt0To|opb4*mcTQ<6|v#x-fNx7uClHx%tSJTH^!!KZ0~2C)t!@(vMXZxku6(ju`a%rRxG`5on#wb65A0j6_x;MEj_hE1?TC_&VqNPQY;'
    'd1jpO^dtEDC~wG*dIbXFGH<3aw2vL$+rU>OK{ePe7XnoN3S@FP&;krNS4WNw7!R5YOU$Bhn<x>?+mjnn+{xq$a;KvJQ*n)ptn(QCK>6idS{IP$@${c5d7BygSXvl?fFH<R94L_Tl<Yi(zJQ6nPOtxmUDW_jo<x&|=kHOtDaQyxe?Olm6d`qwzF+U$)8D+G<DQwWuZj{%u^)Sjok)Z+awY~&j&+9J`(7Si#*!S{xTR`iiaP%D'
    'JbQN56E4r6?B$^I2^YQl!siF<t*L0hy+N4EGT<&oG;mYbl1S-(YF^4CQ5f>Cgd;W;igx@R_M0IylcE!OaSF6{tOxk(Y1^1aH+B~;ESH@!C}+F=e99|41WaN`B5|aHiT&4?gV9EnT-#|NV0W&zRrgpmhO%t?LDAo(-l(a4NhVn8L^Z~(+hj|0hwKJ^+?Y`{8O$Eb)dQ^%Q&$@3U*nV226O3|h2VHOox$<(g-1<iYB5t`+IPK+'
    '<T&*A(gHXO0TwX7_W<pI0oNmw{bOupx5Mu#&;DsR&G-ME;tzm$A<FTuXtp8C!Z)=Vk&$g_xhf3ilb5R{$XXCtrSl+aHlUT+D!UKp%DX4hd7j_eADGK93J{8AIE)XWCM@9^W+0k!H3OtcB;pI$iEUH+!kyf>VHKf5BU_E)$v8$@8w|&kESf{dbPUFu5nF?w3y0^aapeDOr#ioy3I|FNP9G{QnOFP6V*laV&3$Fs#A>Ptbf(%L'
    'ae9&py<T+ne1YtuyG;lWIMYg(rSa=|+>Z{ifR`@n67)uh8=Ezs-k%DEC(mq?Dg~X0h}RNrld|>*jn9AMo9J7yM<^lo8de^A1@(Rw^omtMUryXYzJr7kiP9c8o1L#>u-JzAAG8~-R7Yzj!S~AYxX+nQ30j^%*BmKXuZ)*t<D~G9@StbGIjs=2E9uv#i{Y^zc8=FhH&yrlo%bl5o_Pgfl*b{}u646I<LerbeUwiFbFMIGdll!&'
    'JPFWCJ-gb-P%=?|8h4f-*rXG?Gxsay1~d`Jcwr0P@17bmdm!Ov<nZW0WLj0ocW&h1zaayAwx>h5ngJ0z9?W~3BX*aZV8`sdSCyaaGw9YH+bKxq8iwLA8Tzy5^A(5?P|jNDsU0#*Ff@7;-#$jnJ>3dkfa-2fwQ-=^`k~*_{f)FMq=B_Z(MkJZX^gk&oLW=&7QLZIK16M-x7s75pN}Dh$kopZNYRv9`;3N{`kt7K98kJz8trAK'
    '>!7shL)5Kmr+z}TT){loynD;EeOtV=F1fs5y=W)xTfyV64Q!v*5|9}FrC(?G?C5@?Z_iL_0uuYl(TITWD3mi?MuaKh7OhJM{jVoxnEFj)^6;G@rMLUv632#xMwT)u`23;Z;8U3F7s?>Oi$WU}g(alqI@sk-DVwo}^4gmm9x35d?j(ZF&v1@vc7SCH{=+t7qp3I7%xQ=PXwgZA7l7TQ@V%t&Xe#kTriHns_PSkkB&|5~L=(?k'
    '()}sCW*@+3S4hRZUE^B==kJ<rw_6(hNH*mWAlKOwz0!VFQ*DiF#7v^A|7fslfcgA}&}xy-83b2?`M2d&6?Yfr5hwJ7SjTf*mkIjoTjO&fU{^<IoR!rlO$a>>zsK^e$LH?1WuE=!i}!_j7}Dzy)#`~{L@MJTQtP?^mQ)Sdr6J#49Ro=(QOfD|*5~<M?1nEBG0QINy^w?_7b2=LL|C0g(n!CWx?n(NeNa&(H6W<L2v~n0)t&fP'
    'BQyZ?35%z1Y9@hu%HoYtODPVc$FPeeEX<}<ojJjMf>U1bKci|DR;&RR267m0Skp^dzI5xGlOqOgBGMwEg0;Zq-_ek~2<AoMqzpk-)HCq8K-dE5X%(f?(jhJC#cNhy1*cq7n%&WrR4FZx#t-JW=yT+DEqEUo&F$WWHfA%<hmGdH2K_@dFTs6y_zrTuSiG9ZYS8W!PU79<YvaCk%pqtDEIe6#7PbKba=+8+5>|L)B;_=)mN&Tp'
    'wj6ug0H1R##TYUmc(&p{oE-;4)8(WTN+{!$y(0?|O_8*I0iIc<b`2x(QpM|I7rm@f<qaSSq|Iop6sd3u{7YIRv>oqcjRt&P5i$k+npJ}drtT`K3+F>Sd`F_i|GfcX`Ag0cmeSSpcWh)L#qqPZ*%lnPcTim-iS$$i?FNEyiYw`cFZVgRfQB`|$%a6c;#`TJ1Bhe={^3txMfEfCE)+jWExL|pTu0Ds`qB||M_z(@@_L$eFlL@>'
    'i@-Q2_g)Lr%zf^zqT9Mhe)eRQv3~B-aUtz7<{!Q;EHf}CQDEI;(9=jc$=}#9{hq;z`-S(?J$q6HwEhgrp@W}eKLTSBwPIwKYG3Ts1hDnA7EYeGUT8o@WLz_|i9oujXf2Uc<~9Ll_$HlS{ub&jpy83q@6uytdP20`6^l>jylkByG;^x5<_r0k3A2rd#Wg{8FCKMxaoL{V6!=s|G7RCJbzk`1&g?DmBSzq0+E=tiNe=H^L6|ne'
    '3y1Y=1hbe<_nr=c_D?UO1F)gyEklf5W~oB1ORlORJ<cIDZ``^gpH{yl8V=Id<cX($b<Ff}k)kS<`NB$l7K5>FZ}M&zpEXtuFH=G<AQ%TGP>slJy_JLQphR;x3j&fP_1MH4otVDV-gqIlyI|k48(YlgieDzv(R9n9>rXHBbeyK=K2v*l352QzwN)496zW$m@mr449Xr#U6BguU2P*VDS!w`q$M<(FKF^gSqG^B2&$UAX3McuW'
    '=lAgX`LtPL-E7qFBO1B6q4UZ*(Zexv3b{<Z<;J4I3s!m&{gcH{$MPpd`9YRMEPI<a@u*THeM2d_O*wuBB{I0AoDmP-9vk2LGp0Kf1jFb$2Ac6BSKoWgn>baO@g5Jw&T^as^>iAGw#uJKfdXbIgcTc(+8!GUa2Rd{xh}@|I20H!Gj&zzgmcUgZ#^nD;?inmHk}r_JtQc)JH~`+k?HYR=n=mk(skPU6Ap+msjAgIEBJY1>?Y^h'
    '*TqLfB<HDGab)(5i_*|1(qqwcz_^Bu5`z7T?E?fK!vmjv*RP{PSPfv0aVPnpFf#I40G=sx{MrSU?A+WL+Tvb~r^8(NrAdB1A)p5|$uQI#lur;KR#7&v{QJ$++e}8s$BwfZv|HQ$&EsHN+!l-6DujIBgV6YF*v$ebF3sBw?5w<G4^Co`qA>t8gQ2b{uD_mOQFY<Nl&e2Rwg4dUa&dU7{IluLvehY^H$I(w`Q0E+<wKmF|F``F'
    '4zUVsa6&%yJ?pU=u_}j95ErjSJ#TRoSGSMfaN%b;c~ECkreX)NZCw!U+(EBHK1VA-j-HAIRuEacbh`8$`!-c}XHthP%>byVE5IJkRM9@t<o?O!y~#+NWcsYaZt_QkbK;#q6B2Fw6=tYFS{l1M-@;JDblnw+{LOb2VcQ`6O|+wA{exa8*g3!QYL-E7200s=Qz$z!y#m1HK6AMXrj8Xo7Qq`eYL+q$K5LOiu<MRifVLY&eDwk<'
    'u?9MiC5~4O=gn7W8rpLiKCMyG(R9A{Q5b(^HEcC_O(IK%mr+?k^+4KygpX7SckK!4Y4{7W^0^Y%Vhc>z%e+o|V1^NaJ-3#P1I86^w#LvfGes$itX`Ds{!hM=k|6)7ni>i^)=C)CXs&pd>Tp>j_-+znXXYvPC_4%IJGcgAIIxCByZN=#5#{#c`usU3<!3DK;CSPI#Eo*{riK#ND4!%LU_k*U4eDiGYS(gBX2>+A*Ld1$Y?#t{'
    '={$S-QrS#YjWK*8v-FyOwl>|(quPnz5V;Cs2Tc3~TqI>QOvhtmrgf0ucqRa~wTd2qB-$nZ2FGT{vum79TMK@_u*XbxqJv~n86W<`T*Dz(V-Rj2M2#{%mH>P}jw-VdY$K`_`h`x{Oa#-kJ1F+fbcN{lznhk_GQ_|FeH+l|vIvfDg{lv>;+)0`B;@3@sqWC3lx)rb)F`~i!s8OFsqBnRLz|t)8*F-3#YZPS9E(bNUnrAn^H#&$'
    'ZVT$!*>sWh3Vxs08)v&u9w&(Ns6Ob~x_G(ax*J`*fxo9dyMh^)HxP@|O(^bQgxPT8GIe++6n#;sI<?tvq4Wx<R*C<K!=e^Kf!Y!OghLh1HKGjquYIQnQh1j^{Rn(tc4&ECcsH4Mn=tuh-3szmZ{!*75Ag`C8FrULq8HNc&PhowFI85F21n2Rf)&d^7$p12MZBO9><t<?S*c8$`oEVzSD_Q`E#**nH&!Cb7S_|=XXazxC+J@B'
    'rY94ZR}e;@*kCg{BxJpm)`q_gxJz#SF<R7X0DTAbl$i9L2g$MR$lZoS8O4VFdNsek36J!YbvtcVW0;;JScW$@r*V9Q-fCEBA%KLGuVM`jK0P?F>vPk5u;C;mt<tC$QO@a;_?7AWL~!dPuVYaMS2l7Cz1cKw7-hsKO<K-^d`o-3QGlt)%M$4KXI**kteW=nBx&f3ZugZ9c3`V3zH7fuYQ`w*D1L!4W*XA_-Kmu&K&-(hlUI#X'
    'l6(p?S;N(Z5wTkIfyAzTs?N@ZQ69u;I&(m{mwok8md}(B3qL9m3?$+F4U~>8&kKiJ^h)@6XMs1bEmp498@-#i!=<bIc!15O3R07XVkfEQ-T?;~wR{{>(eYrsz+;D0`se~jp)jHo0|*eES4Oh4N?dl}MvHy(Dov7-DNLw8w>GYk&$Du~wWbe&LP;{y-d)~`=koc=%U?@-6cP4Kfalzqb5$OP)XcDGI^5T`?XGf<4DDES#JVxV'
    'E72k5+9SZSX7p5)lWqYvg#Ma@%oowzY}DR%wNI7nsL;)L`!{}Rciztg``%W_VxZXLDlh2N0gCv~jsZ${Tv9f7_xIXC*LsY^pY$^}E@|j$;yD6^to=@SS$X5{<T|h~6zYSD7H&#|4>-=aVmi3R2q=*yG;p9InXt3*TiD!yaRM7d#Pu94`ISfi4>tVlKN<hC{!`mLN{N%~an0^nH#2C(aqz&IcZG1@N;20jt(a0P4D(1z_m_ZC'
    'Xw7V&NTYT!{Su~JXdjm;e~aQEd4o3C{*l{*I6qoVtj-%yy-WD3*<>-TL#cJ8vq+vL9;yaEvzFBfxs^V(^mJFwymW$hU5wqW%F4ZnyG-QN-C5kD2O1_Z?O5rfDcj^QFnMur>~g;w#7E<Jad;V$*FGp!)GL|+dnM}M1ng>_VO$sbm~<7JBLYZ8YF<KG(Q6dxDP!~@;Uo7vbQo6qRF0Z&bZJ#(@*kJXboq+i25=ALpIT1?*KypS'
    'R<LZ}*%@IIkkpff04lJsgVHd^l+qU_F~<a+lgBY`YlwrP7^He{ekce6V(Ya(E!9e#r8QUZ_|-V0HaRVCr&n{s_l%+Lq8Klus)z!sOr`6X^eB!{(k2q`_>?_3gy=Swjf~ZzLgxoL?HwX`U?fq6`-`#hps2P6^1-1BvE>#w>6+~pVWnT*hAIx*r!A@!rv~&B_>ws2k!&%eM9K}IGG|9?Q9lxS9)^v{cJQlKgw@O33(b!#b195s'
    's^!)PyIvkJeWKm>pMb~Shei3}%kK`TEYqXoJzWiz=+k9<4T73g2y&pTSLkN1pPjH7s%8c9i-Twe|A>Y+tQQMAI~sO0oVy9bnNO-9Brc3kL$X0v6jo}c;Jal9SY(%^N(_+`P#ZHDiU~w6m@(aKrkt2aJ-tHHFr3=JT`xD~imKU`St-LUl~2xtMDzd;%L?E(&Rm_ZG&^`hAN=lWJxEtF_W0gLfE1{ZgZ(^GIfJ%on8O5}hPl_&'
    '+T)@*t?<<obF@}xu_(xySbzoSJ=HKX8rum+Z&;~xxqHwA%zJ+2_fPq8n-9T&d;dj~I$P%g%c*VK(rmRf+pb>94(+I*FUz<Iva4(|!~C~8NY6Q_3R9Kwn3+cn1Th_mpIA#Y54Hi3_^Ch{4YNu-%@&^1Fj>smJ>poyC@Ej97&gH9tA2PLwT!Z&q!^brNNanW@9y_&a-yqU<s}7bh~%ml4L&y1&WAe@q|_U&W;vvkU(g2;*6bJ-'
    '?cRHDW(w@&x7v_e<x9ypwKV?WqSsa)WAVe&m4*(}bMua&ZD?s5^;~@TFYw!H;&M+jbZtm#*_z-d`~ZCZd|DW@VX=o8^SI!bzJoDyC_#tg7mH~n3>bKDA+Ou==E1cm6ATtx-J=ZHTGag$A#|VG%}wxhYMAF0?m@tgdaq-A!VhAl+?|0=a_`+D2$VkNdfy*I%_t1)Q2K1B=0PT5&1hZ8gD474C0$PZ81@KICUx=DMI+#Ut!(3z'
    'kkNoIJ#8^{jg5gZOtoI>@^yds6PQ?jZTK+5qV9L?Wb?7vFaG_?op$n()ul+1!s<%i{1~4M<kQ4+dNBVoI$VzI=$r$L^bTRLx^aZA3&6Iqu^>P{Ysw(dRpj;ILO*V^#XQUmqpSLW9o5U!Ypcm)3KGT!baRx&qlSRQ)AJYj(FGC-Qdo3cEE(H|8Im9w1yBIIeeSkvXjw7@2lBXk?Bnpu*rDf-$ILSbaxXPpbgidQMH=d072DuU'
    '?xjujQ3C~VWpvG__<!r@<h|kGsR7b8Km2win(bKnT`h^i-cN=_XhhiA6K^!q@_SMr@2W=2T036<AbEvB2^VGfa=1Na)FAnX?uc%>%?WtD>y^3Sm#x&@mxq|r4+^wHSy|;~)3FlOw2)s|`9H^Q8FRlE_hLEn*k`WTDMOhWZ6!8i1}{|k{$uT619z&%5hH<>ZB5mR#s_bkhZYFBwxkOx1Iln{yQ>^S<!X}K(slquZ3SVGgNtAO'
    '@<PO7xoWlN?d-9CEARnO4RWW}nE=F*ve3OJ4fDZ)P&Z;yV<(fN5LM4&fFiAkPvPyFrEzCaO6RN4%s+TW*0jiV9tSe)z084)yBskAx#%qzq4Qv9;!=DQTaGvXW{-jEty~95lm^fEzmaN-@3~0$T!yn@;|Ja1vS61Q3g6uAZYH%`u?ZGu;%A-4)=J9(KPBrzyV!@m*dHZ%dhF=B`T_-MeTpZBO!9F7+#olsj|wr8B=B#mT4w6B'
    'o(QWk=!8sNJ)WRN#7Hz2Rzm4~;7lJI^3pq-7DTf{=wU-{k@NO&<3^~ySI;Y)Dq?@+mq+n4+|u#kEh`OAF2}k6MQ2DvIpFTn!Y)Xb#qp1QpUk^N`3FKQs%y{V5aaW?R#^4u;1Kc$A*kPn184vG7Ota*(brWcM~L;GM$0J*!a+^|Mw<ofxyWY*I))VYCX7=dh6zKmr6+0WO_j+wP!xMulgfckdafguuvC8j=D<S8pu(tpR3tQ`'
    't#l(_4ZsExD*CeO!;f09m6~fSDt4eHVSl;28G6>oHtpa=vavAp#qpxevKylL_)|SpM<R2=NIQjQ0F4B&rv|-BXJCyQMuTM>pwdq%I9{%$XDyaLDuiy$`!sM=wV08P4q4y7mTI?*(wUZUy-qtG<{}rads-E@F-oRW!JU|VoSXWCutRbrtd`&>s%-D7CNIP0w}fr)r!cOl6+3u(+(155$zJ|r#{pr4Op6QnW$DSZi4)`tca!84'
    'U!xBPO&!|C#r8bn?H~-!O7TeX6Zm5inCz?#efChoz1M-Fa$$`hE?ldiNBoprOocwou@4j}^eF(I+Nwm_92NmYtGq@Gh*940SVb<u2HjPw55HX7Lbm4J<1}yh+%aDXAI6z!)vh36$la(`_M<rFdQjGek`?6%;gL05LEw*kj_aav?ve9h@hyjA1)&Q92~a1v5+J=gPSe>Tr7PuApZ!dOPq+I?Ul((vX4@+sADwdC!FvuaN74cE'
    'B~^-A1#b9%V$=Um_;~bH-_PL;UlIpNMM49jpesbz0OQI-iccq-7yE*%-CBkxEpd&f%EDYDLT)E!zq?=ysW#;ZHy4mRq$7OEi`0_d#7_*x+sMi|=U%d68*{8!)_X=n0cVXplRK-o2Q&c0Ep+NieQ~y4BdZvAdx`jr0f!_2=Dx?vhiS}FQ;pLZ3+`kCmJx~AemDg(1n~FBbvD)?97s1jMeFFuld|qov_zVzk!7ad8l8>$m(M7y'
    'Hhp`tJ$wf7Z>wOA4IE=pODv4shCoK2P9_VqfCN~^on<~U3K6FT0CW)~E9{+mk{9Dre4ax~_sIEZqUn-tGFz?#cS^n-xIOd_@Y9nM8g-Hs80{P?R9_ge@Tc$u-ogVYfapTKrb}&>>kaQCX*VbbvJWPX0Q?}ZlnB*}9rahA+8R<$4fBqZ%q?V9aHp`b9agy+0Vd&l7uXUGcjQ|*%7L5vU_YPXB4bY*2Sk@GRtEsBq(o1hYSSNK'
    'lSD_WY;F!`Y@r`~r-e4)+Z9o=>p`&ZE@&NqA0~}Y|INu@DF=IN7gd_fn^36GKXTDS0viBz^S~cFJDc1PNqP<+fG~C<Iv3ed5fA1^<Lt+Q==JKtfdm9w)#3?NEm7<tDx(;panDz6D$|dex5|}T3m-zxyW}t-p{u~)55)up5n~^0=S?%Yzpu0cUD$)XddYpiC?GsLM7{~w5Q0bn{^^HvmEB-xKLl&$dOv8tM$aIPZVgn>2$Ndq'
    'C;Y`U;S_UBBnoX2z)vj;o8{5Gs6kjFK71svZ(EK0Y)Q#ko_*YRnPYU%hdBg{xJ1)wCD?&(V6{p9N$_e0@!=HEdO0+34a$Cw-fzy$ei8zSWqlc_{1IEB!{z1m!9s*W!i~Br?EBX{iK-5PC_cD|J~S$VuVaLp?X9LDDA2J*T6JqB6akfsz#P84<eo5g;b(unLF37ff3Cm9DU^aj+!*0u@UMX<j25P3$-TBF%T2Ci!TSOxjP4HZ'
    'qoKZ+o8-slg-!5XT&xHdlI6X~Whc$YX@*~^-tJpew+*1~bzElWbpdZqeD}cxH34r;v>N5Iuf4L58=)%#cy1UKpv{P0Z(Z-}c?LKH3glHqog>{B)`YkT^@?dIVx7VdDJ$%l%9X)7poXI6YA=|;Jm<T18Pf?#u|;?eV3bIh6+*6kH0WJ08+f}7x(@W-;BDAS{SgO*JO~(`&mG9Q|E}94pVln37?|negPgLMfJQ*^QV<Js(Y{H-'
    'MtTt?9JIVfoeQ=%3d+mka(M~**e(&{wFUgmRbbn#&@@U`40G`z6r}RUb*F+>C=i;c{6jB|m>Hf_K9#`r&RWA&KM<8ht>xkGp!*y&V|e%?3K6{tyLJ6pli*vscUkjR=wyilg29yd0?VS!3snYL?Nb%0cuSuH>C>WbUIzIk2KstFEB`C`JafVjNXsfC5v6s}B(DRMgXOGNqeA#eAMMi;|MZIl2adnOX^)LGA(YN*K@Ha!#7@_)'
    '<-Wm@uiu@?Wn03LrR5*y<9>W{)f_4TAWqkPSdO0Zcwdm9mc(4dJb2pH%1CTdDH;Lv=g8_e2h**%cQX?!>;yrIc%6Z9t%WmI8McaE@mEQp6*+!bb8{+4d#cv|F;AaVb<7tNs4?fz$pH|ZS=Bmbp-41&S4s(gNp6m%@_~$&s$SV@@THCjQ%F7ozK!$dy;AEnneqb^Q2+aMW&-r^YE1tbZKtFM5mw@!@9oLknzrU^<4r@072HAy'
    '<<%PTIO>u|OxuiS7z&r*V=6a^IK{~skGa#p-yAJnp{>ds2S+6c*@R4w)?lO|$MJjpKA~uhjrVHyTD!atJF<U*Dkbe+(oVR5A9*%}IVsXObMnm8316B&SapRKd$3xJn7|TM^0OdQvZ3IvnsJT2Z*`BsPGRU$y@XNZ%(q;PdY{9cQdCJ!SY^WuEXLw@(yC`~bx%?bd6Er7l^@y>A1;gYV;l>3)cYVQXRWI;!8hbY!Hi9|or0U>'
    'pZn{<M-lbvCjC5AheZb=2d)VA4noD16(Ji(F;uKPG1(cw6m!?KFkU;$_lfUfb@*12WQ$EIkr|(09Ws2i!V1$C?#5YY&98L6G`8%%T$#GQH5hGvt?k<MTr7v>!<|4cs7E8{7~DYqWQJpffF=~BBWIF~B0}O`aF?>!mG>nxea5}GAGuUsOW)FN{N=%xyc^WYZju#SyeFAVSv_5|pTXS7j*>Uaa@z|u@_Lrn|GRo8$*<}wl|9P6'
    '#+*}dA=d)0D;iVik@>N&w3_=U(qcfGJsl4KnFs_OBAOg^IVv0~@{OF+^<u$2<J43>uD60%?=XQc_$yD`Sv9)s%f<%lC_G9kXnnisv7S++de%<P?SCnpzx2+vJ%yIiew}k7S@2d(;i#}RvIaP_n;45Rtq1Ro3bGAO(D7%cJQ6c(g-PfIMMM}4Qe%JL846P&?QnstR4iDg$?sLY5_VGTSK+SOl7(dap8{#C%}DP9q2H2?k;-EW'
    'w7uQA6>&+OyrAGhtz=_m%+F&sg@QS2{|#;${A*6#!WP@?jO&PPp5~^QDu8^sCEkitnv<STFBrMWL2gX}WKg7lo|$KkU_-bUqtf%7KuSjAgxzO7hUk=T`^{efMxBjQOwi16qI2@JW5`ehbig|s6G<GgYRjL9P7G2bG4n_^nVza5#>9!UoC9vh@hD*c)SQZQ#Zuo4c)DUXmT6)-Id+YacMAL@kJkqqj$aGZc%?`(7wg2<p)u+I'
    'F?}!pUHQ{ZrSfA@F&*JQ>x2<P-W0S<O?&+EI{UU1K*i7sbd2q6SeEN;dYn!KUAlBD!LjL1FZCeKq9l^W?%E;e6#i-Z0;#7_=m)dc6L~)_eoXj{xS$FpgP5$%iY>|#8<!Cv*!3S_4QzZULmcJ`0vjREF<4F!!FEzl)fP{2FHIRz-9IR=wm+R`Mz|0;VTluL5{Lohlht+HwOG%=$zPfm&>8YSBHE}4o>LnME;@dMr!1Re_6!jI'
    'qGl&EAbcRPKfEG0or7(-laZJJ7BP4Ni)e=`4$|Y$lR|17I<*Cj&W`=e8qJbqxguelr8mT8$vsuTKs|cHua3w%GAb)A)Y6`ttVy8XyH$l;(7pC}J$a9|Y)Px=rJc5xQ2%JNuu&KSz=aN{DB-1*0zeZYTo5AXic+rCnFsUh^jbit&hOC#Jy0N|nftGqkS6^;E|x7$0j{KUoBwvlcX-allC_oN#?*yOkX!B8f!Q9kA5V66dPk%u'
    'be4Ok43LJ4JmLEncjRWiL;e$lbMhjSJ2pA4trHAWFXIKFG4^NQxIE%Jt_n-bhn1?WU1^3F(i&*H3>i)e?cJwEoToEY#DA{1Hh&9<_G|a;SPck^DlLrfbeix{r&cbe|MCQ+b;XY`q+?i#477MMl<wavp?yZJ-J~teJ;0zRS6Yi=Lhzgi#{C=@@HZ@DQ2@gc#PMSug7!!<*F)HjX#D>%V-tgLZ6}rG^n03?(YW#CD}aj=4aALG'
    'B=B|!$ZUTf+cPT0jDqj#pqt2q?j#B?U^3oh_OML!USHMRo{T?|2?&Ef_Oflg&Nhm<VqN48i_+9tlhrR9YZ4ZStIqDZ$GH|Ax*N-9`-F#HAe>)%!O6dpI6ZGN7qKyIs7OmN^gxtS&G0h)ywLrP{6!Q*f|S(z8(&bzN3_X&i6Yw~%D4$r7-n#P<Zk(o2N~PDQ=3B;KjTGtdF&$1ZGIF^c%54SDiP4e+%89aox@c(msz34umk$h'
    'v<vFTrZr>E%M;;T4I74toS+u|N$prXxrV?ql#Sq@04obC-SO2rTtws{q5X4#9uid9(BbGhgl(W;IJJnIy?8pfV<yUvqq$Z>F5X*qkLgWU0`bM2ue5(dd*^90VvWLGC{lji$z0^}BjF<VN1EHtg|!(8m6{iqt~k*4Jy4XE27Vp@yWUP1Jp*A}O&^n$w1xx5aPjzy{w!mNW7<z$D<`v)*{J_o3k2BlXC=o_2bwctr<qfHwd!=^'
    'HUA}DvIXo?VuR)Q!*M#LDotcshJ6{j^mx8cJHcIiV1Zjt8s{`^Dgj0c{T+#Y?t%X$UJYxzr}Eh6_ZFH{l!z?@($2rj(R>oJVC+AlL=f*lHs_Iy2Q4kAYzOkGfGT$!M$_Q!kw<P&jEV_C?uy|bN)%lKouf+4<Z&f<?aD9P1D7a@QALE1!fl(=!r^|?S`j9&Vy@M{&)k4jIA$|#6EBCKc^2%*>s=AABC7cxL19c0Fe4B)DUKg1'
    'tr=LGgz^_I+uXZ@aplBCgH8I_<2oa~_=q%tTKXWd7uTbi#uIGeQ~eeayQcJj{CLda@y>pB)jXE=I=lYh0I$F?21jC}MR-8X$um!JuyIvX&1t*+_E(i;&izq(t3MkobC}N|_M%`VOdqlxB?T@(>$u35b))bB6V+}p4`FhojKk*$#)VV1bOAdoP^C7@55UkIgvb1)aS!3{2)V>Wnm!9bSEcf~fhp2QEk2oyc8mG?6tCZEZK&S|'
    'fGrRcx|%cDmYN$4mZP&TIw2`$*EpZJrcX~FNr2eK`qqtv6Qi%rMVte9oM4O1!4I9=;mltCs}#Hj5)vrm*z!ksi7?yh105e$I^>C9UilPg;V4td5xf#Dopn~!y+YOyh_=frD=><Ue^|8R7H+Ri$IcZqJxxjm^@<GYse*7Tt{0{{3G~E{tP=TMY|Rr*sQ;#^4Ub4!qG`S1qmiXp@O}>XwF-1Vq+gSd_?TRkNl9QTMX7<~$3u$f'
    'Y`2lIoOM1h!9X*(d?+*AAdlp_d-te7Waq2)=HoLhGZT)*BFlO$uTC*6&_N_jMPDajLtT>kgp9Ym`Ew<^EN_d3S!m$3V0b!FG+}t1WKP9dhooD%tj=W(*VRh*v<D#TSZ>s)pS~q9wp_JqOMme6V~v=D%;Dj^BAUgLp0OY>L<XTt;=>l!gyzmznAB9~hgW@BwO`8E@IH+(fr(#@r#M0-pXqn=Iwl!?5vTF}Y=%@<1H90EW?Pvn'
    'f<)nD`a&7J(J4`>T~AEjZ3XO(-RqAkm=eBnXzQioh$$z1ifh4VAX@Udq>JXMmDdO34_LLu@TQUF@q6GdAZ90z4mtbU6GzejNN?baA-}nf%tm}I{QLM`Eq*lKz$rehd=8e4_AHD33s@ndu2Bp)BU+<ctd+34aFQEyLY>9gUTMpg#lh&Pks2ES11s&z5}i8$nEjOXRBsRT`+^aWZ*t#qQ+WlQR*btXW}oXK<Jt_xgYY~$tI;UJ'
    'ld+X+lKo3AnUjJjO&TBF55WfMp08Px;RwoWeh&NCEM)?{31ClwX|ve|pMt|TgUQdimcV8Te(DY9D$)Hj4gAecBX2oo8@InkJ6vE?X>z-Qoo*hY4S`X|TL+?p*Rk&rO;AQkYn1lD+2k6FF=6T>Q@t}RM`0R!>|`vtV56khg;7VSL9kX=uo}^12uO63xN3|KZ2-)VW=FD?AQJ=32XZPyv7->>>yOYrz>chw2j)h4QLT5l;Ej4z'
    '4!NbXT8R4J{dYVx>z!#qU~<hB6`$u)>66~}-cD;$G$j(HjEC}LfT=im>-}8SVPx!S?J&MytiJ+<N>+g7G!Z-CCsD3qyF}*9C<}40IeO!-!kYi*$anHw+pW#U)aAz9IJ3*1PclBM%TTC;F26RR%a`!QlCg77z_kVy@aHypTc#|IIFeDWj}E?sQo=8KFy1WbkhSau?5=CAfN7K{<sJg2bIpHLQ5}{<;clJ1uAjgiS!PK3N09EW'
    'AoTTbYRjlz;x%HcvbgO>avI>xf`@vECFupH*ysTC@u$ddCizoog_+;|gdW61F_Ht&Rk*qB55)AnU#nhD(`!va2R)K#TXNCY$Eks~*JOzDTI&}`vJd+@?v2Yi@TPUA2JY-rFTw>+AcKs3%=9!l5^$b0sAc~R6adOGJ)kUtH5dJGqgJxh)Wb_eT8!4()MmWJ0|U_L>Rsq*jxiuL=fgddqj}>TEa0hA3v{c`Y^0IWI!s_xg};qy'
    '+vu6<OZU8A6)hr-@zD)gbXy_uKnXLED>J>5JYWdmNJ-Xf7~7VI%ZUoc{KvJTD3?HXhefB7o3nXw65kGEE1tqMmim05h%CmW4LsHealp}~f~{~4ryWS$@PH@}At`;B=4aZp;srK`DD$wMsU+%}g>6&=gTnK~`DedB4<hw?K!{s;$W5S21J=XnA<++A={`4@CF%N~KW)Z4YZRdL^?OIs)0>%!CBXsc9o6YeT7NkattJv`#_~!S'
    '_%CwxPy!)Ua^DKY>GW<OksnzGUi3>`F2)ez--rJ#5XexF&3IOAI0Oa0V#qrjD)~BHsc`o>@bS&me|IR2E8RGyn#3?QoY?e6pcbHhXwhJk2i|z}gzVDDn%R)~#mS2dYZ%D3NX6p{#5(lqmsqI?6&<%JE>p0}Ey_qNoDPY*n6%(7j`)ca?ig~}WauSj_4K6hjP~~6pp=ZMOTlfc=)_biO-3|?!s}4V#&0kvhqj|GYPcf($JJ!O'
    'm)AyymG<Ok4>BX&R=4R)RmTf2#L7OMg1x;=Ki^O1)9jn5HWAVX<Yt>7vOF}*Wh#Q9v<V0T;~6FgW6`busWFrWW{nvl+51iM`4n&icmaO*IpwrgPHf{>LJgCcPW}?>7lFtVum#mx%OB^`T#4G6<z?qGdDfHpvvaQ?4!4V}%+|@+Yt6?q{plts{<yMs8RdJt$``p@1HnLGM5fz@VNOIw3G;}0RSM3#O??P*XuQTYRQM7H!jr(J'
    'p3mON#jLhAMQ5e<b6w|-Yvl14XNFc>OvpfUlCKsm&A)g>H9zntsxfvj;@wIEH~Bi?ma0vJW{CfH)EddgygpjP8H`56kGg4X`C1{e{x(7+ieN+JM6&P3i&FzGJXN(d)pdzD*{0i|xc$`huE7<Ii<v}YC3Trr0J?+hlPw<*aJmavewlpeHK&(GSsngCNq4uPLY}Wv<`|-Mo`>ibP9$;ctBLg#qgh?wVh$`TpXQ4&rx+!Jwg@nT'
    ';-kBdIr!iZy#=0@+3}zJZ!#e?fvdt=V`{1bvMEVlsp}4>H3t9ftlzQe&iaI8ebK`R)oE+w*Q-^)w3gLpN#9nyLcD;Q(URntKP7$sSxc|hHBvy{@xZJXPXLo^_s%K*1BW}PWhyErpI;FkB&V5X105HyYtxgtFRu1Ln=I}V?va%cpXBfB<X>C?)t7xXp=c96qd28>I(-_Ta`c;Iiu&#|chc@@Wg5jneYJd|swm3J<VAcSUyMO&'
    'I&A0hiAsLz88M$T;r3C{3=%u4P4kWDS`^O@I%VsZ&Qk>tTFmhJvD}eGZ>O8QaP19GD%1+Wjn7N8se`lbI2QEOyVE2EY(F~?&V}b8d6i<<3!^yB<h%&Mff$RjWDE5LRiF$y0lCNC;Y2tdQ(S^b2u5jRInP%YsFe(+B(1*t)z}DAkzxGK4P<nSnTAesOuu>kTk8YFT?_TsB*%T>sZh@LJL^KnlF611O;VuAD`vt5L&aum^>8r$'
    '_|384TI+p(jITP&pg%`+@%2ZR+UOQ%ghNiv-DV58&cTVLO^EWV#$QwaCUwIUfhBhZCW@xvzhyM1l*en3SUHIYVtN$)R$~8Ai~%^vzZLCU9J31XWk!MF-_=WSuPnG;3z{Zzr4d!c5C3b>sQVMXqFv#u8t+qpL-H_BW4*|M8S_!4vTRL=xGv`7vocK0mX%7?F27MQ@!GK34^Dl-pIu#qoe3Yy0uN21xCR<m|IXu6Rokp)j-2aa'
    'S}`D@tlR}$i1Na9fc*bEoSm<04lcUx`?bpy&s4fmgveCc57c4{M3DjN7-d*aR`ZSHPc;Z)+*V4J$jQk<-F8P#w#)F?&?MPF^+1hOKd*nQoE3_P9A%~mge?9sbAqn<yEeX6p7;<p+ZXad)_19gvH3(kJ<Pi*A=Q}in5>JvThe}Z$*K2OnllBUmE#P8KdlbqY21H8XJCCg-WTs?PJtj-<?&D{d&NU0xr;I%XKw54_3Zy)XSr3s'
    '=ebbM&Z5sL4lwA{<%@KB1pW;>?qwfcX!nIQ+KHQ$%^N=nOd|*WS_{Rcp-D9s)0L}@VG1{EzlFUP@KrcnvqibW+&*DZO>MXVVFK2;PnOClFYv7=H@X0T!Y?_x#5|F{1_<fz_5H+$vfKz|P}Af5vpx9~xtgtBEl!1*E}rnxdZLkaVZd$T%3UE6qT*Tph23AGi1k*kl$RP&$<eT21Rw+W*)Lav$E4^>$&v&I;EU2xlz(yNDWCEl'
    '3t6s^a1}U0tB)qx=KIpQk%$_}p=MUPC!X$V9uS(lDA0VB#L(ByE+PG1**zK;*v3~??mZOhk8zQ%^>5PrJI*2$iV3-f3#MWWFa=Ib!7TL=<C2DXgz#~@LQ;Eh57*e<Ce0FF$ww8qj4+eUfXcCnw5k%mHLM=Omo3E83vA|gfg`UegOxAwN;(7DrLC9j!dWIKTQ3d~`Td>_TGUpb=5o$2JZ-We@g|L*z29cN$eY9vY%h3_oAH}z'
    '3hcS4kAuN~YdvZDU<3Ukxg(C^yRvm+?IKhr4c|kAxpTkU1{2x9dQ_Z%ohA^$55@xgSqdv^TPCn}2YzqcmN9W?7V4f#BfNTmW{(?*@3|~s0JydvOC-`+rFK`-J6Og3y9FT{l&Y~1+tO5>)wA)0`V;E%S%wXoChEQ#Qc?+5x<D3pp&P42*J)i_yL8gk&C=J((}3)z=whJ!mkfS2dBz0>!{#4U@rK5F8fHrjL;A~rSWk(ht^LKf'
    ')fL8U*O-71>$yfE-`&t}t|vX>X~Ax8RE$Dz*qNH+tu1Xr%OVMm)AkCx=>E9y%70jrB^?0hH-hIdM_k1~Bt5JQW;$F82d;9u_!%Q|yhbdJv~v}X?cCztfG3&jD4z0N%ar2cQ({0Ty5=aE$3ql@(L<{8I9d{HoyPpO@{TUp+9mg~2PT|upPA0XjCn)*4~`MpZvbsVA=g;WSlgLN)9MI+HR-9vks#S?@7y;gX>OVOjb0;B>z54x'
    'RGx~3!pl<z<)#!><y?f7ENFYMZX7Np`^@w;!CUzK29CTgpBP}lTZMjTc=VXD&iksZ=$Pbsqgf--ke`gSA}L9WIi{&bkYQj{uXG-5PLp<gu>?T_ex(JVlrL*wcORcLy%5o|k*i_u%9w010B=O_LkA?!l{kjD^5NnU{ufS|!?Km(ZyI!L%wgxoWPIiWqLyrTc5x;@i*CG|wE{+c5=A;JNjR<jHA<3RBpw$gxLe_jz{%#<*x&8I'
    '=4O4~yp#s}cAD}e<bf}({hYS!aIy`*4<NglJ97i>A&uI;>UhPdp?;vCne_ys#r32t-%hg2WrO!G)oDq+93v+tiBmye7`Bx$@dcZdS>^ztFHh5k@9ZZd8D|zy?2#53B$_#}3q<PNK3;;zfklte=Bm~4=j-fl@B6svGdyY|#CIK{tT&bu%*R80coCfTb>9mq1Ml7KNM@tCndIx)ijB<Vkuh|cmO_iLfc3d&AmY*I9;k~JtYy;o'
    'D6n75k*(koxf>YI#XG?ur-|Okhdb<0rt&G62%Cl+w%4QMvpg_4o@fr=;RPLfyoIsFj>`ZF<lPe#<(l2P(}ldKma!KU-dA6W9kJ27&JT9zeuB!bxB~r@WaM^QLa_64y$U`ww2Q9kq9iR+b3F{F)~mtZ<j=3Q67F?v%>uS3;A{6G0g%lig%gKLZ9kz(<xzn%dr&BPnNzd?UHd^y8hmDHkC+Q3%<6&v-PpVn?7@ZnM)UHQ<w_1M'
    '*P^T;pDCd06^ZcZ6+m!je$FC2m599heCaCj?LetTiRn*-d!zC5RlsqF;~#GlS7zD(`q1;2-Ca~cE)JL{A8R9RHR%FChHgSrFJe+ISkIPZk7;>{T(Dht8y`(UN<W2<r8%_1+@L*)Suh9XF*YynR;Vw(GH*Dfqh~4j0>6|-S)F)m-rq2e7609=pxL<;!>ck7h7hVD!={U>!K~T`S*cB>10s?E-L*^C4qd=X_OEtbmT@z?TJT0`'
    'NPy(hD)qJGE6&FtCW^S5eehzrvo2)e>(s@6WhXtHN$qI_m?WWHnSSn!w1F696hO-%DuqHd#f_CU*{VDDbhpfUDkD+S#e(?;wDlPRV5n&(AyM+F$ZHA`Wp#KkuqS>ljDUIlwmrwES9aR`Q1XA<nNF6~EfEXXNA#AkAQAx4rV<NJH6B{*^j+d<`85=kwBeCk_F=kSNs*QeUZWT?@Qk*_9}E!EcwTa-AQywz<WcoR<KKVxSF+0W'
    'I{H5%?&R6${CQjw1MvX6$jhCOAQxu%QHKvn<aJy+X3R}7$vL-UW^H&qd(NoED)`jFspqdqkEgiUyVSLBS7k)z(o+Go=}i3#y~wF0vACEV#nL(z)Ozs#H4ydjkgn<HeTp^w1C>zPd-d12<_-3wLEwlMyL=bQ-^vOztZ1d;E65Xq)QxU&^&wsl&3(JPJgiWzPK7f_*$J_v(b^C<?xBVPln#o%of?<#Z~VfwdM&dy&MB-{=WJFN'
    '=GrLs(sF-IuR8h?6R_le&mIKkr%eu-N7yJHviBobs<<4jTKG&KCKNBejH`=6sj@LAycLWlG9ydOSObFWSuV}QeARG@f2!op{fm-0#jI|&7gP|rPw=4g|7`hao!NVp%^`KE_9cW$yzbZ(o4wASt(*D@5XnpFuj^<|X!GJ4ts^dTxAHAg6M?hF-eBZ16e15m>n&(e|CaZyM)H?x`1Iq{GNa1MXEHA7g#cybS$L_NM<nv*F#4y!'
    'zS@xlu`%Jedhk9QWDm*XNq4j4EP&(R@6qS@Ucy(hxMa~caIK~-bault2qyNGgq`$Mi2nxHt{Se0LpA2gdlcJej^akD%Kd@Tz|trobsw^?KYDn52TZu##Pk%Ekp`$c(um=TdvL^3`Jb1YNSoC#+o9oj)$Y8vV%3s?5TcAGI?3|Nyc5VZrvGR%C>h-#cR<5BFpb+DHx_2dwoy3pxq|v-*MyMV@+hgY%JtMn5PqLXQa1NEVy0PT'
    'b4$L@t)3{qGZ=v^n20C#v*N(1R;H23=!ZFY(6bd8ho#{-q&v12?augOg8-@4N#>nH)tdyC_1`c2N;I?%UruMrdb8)j$yaisoO2HT7+49rq6h%JPU2fx0DSPhQaj^VRA??6rR{a7f$8YLAEwM|=8^YyO&gUAb`qs^>rX2e^7v>M{~2A{WSlK^3yFHTU`=cm5!b7uTYE@9&lW93s>+L-3Bl*v`O=~?a-(a)&Dp(inszfT+qxPe'
    'QdGO5_hi!*<F&V;1DP1z+_n;L79X60?PJjL9GCQ!^WM(_TiJ=jhjPgePVFm)MZU^X=lT4v&eANFg|2?buTBxNovf9vcaf81cchRgONQOg+1Ysp;LNI6*AF!Yve3w+uAm+5Q_Dq9Tmb@m&doIz>i?!GhC$r5>*K-fo|!9!J`vd${9ch$Y0%2#&0R^Q7-WGken#Db0cL896YwW{s=(U=$6qNxuz^{<X_hDmV~0-vaCuuey9z?Y'
    'sX~}M=bHoB4LAu>+_120&q76bB#43?j%aEDS3nmj4iRsZ<3+PGn4y33D0EdjMG>3z(Ur-qIp0FrdVHHIExCLSuJ|qD0#ty6Q*;&ssq!1+vP)U=@y3x><~ra{e;!!MnCr4qJO4%+lEdNr)RS;)uR+Abg^=>@71~Ms4S*UrlsuK%r7bM82A97`8u|cgU_q{)xu_wwY8wr%4jsdEfEC%}P$x|UsOEtJiL6v~_{j$oA;WnVwu83`'
    'J_&~#E5RA~^g<fuuIXcjan>SIj@;!kNi4kXZH-k9Qb3ITmsN`x<tDb*-Nx0xK=W$PLHn&$`QM~rSN#&98-mikXQ@Ek)xZ6h8Iktjz`kRovE%|B1%kt5TA}O?@ytBX2>gQwAIp3ae|^OI3}?XtO?gLIs*i6OhG-GH9+gqfsBma)#O1M0ZsD2vl61QvLYF?1_E;dVi@`mHv(2EEp7Fa4fB*m&aJf675n^tR<;QO;xAB~wbe6*a'
    'U<nkuD4T`zP^r*lH4_1%gbFt5A_QP9Ag^b+j%=-eB4fNlwXOLKKt32keg+2a{Z39F|JskblvI($w>&Ojj6j!ZXBK1oQgaLl*kaTzNarP3M{o!N!1MAA@2!@hEK*eTIAmc%YVNQn-hFWkbN>HeP)#5UAeA%CEtbcf<X(bz-JeGaOikMv{$4`r!g6NJk}AjR44<WAN<=<SQ0v?GhuP704#QDWV3IPIwnO9r'
)
PREFIX_RAW_SHA256 = '0121150faf855dbb27f327fc7ff857a17b473db02ba93f7cc2eb73ab4972ef10'


PREFIX_SWEEPS = 16
SUFFIX_SWEEPS = 16
PREFIX_UPDATES = 4824
MAIN_UPDATES = 9648
PREFIX_DEPTH = 1 + PREFIX_UPDATES
DEPTH = 1 + MAIN_UPDATES

def transfer():
    return Transfer(225, 225, (117, 270, 4), (126, 286, 4), (144, 356, 4), 580, 710, -13980)

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
    trace('energy-start')
    conv = transfer()
    trace('energy-setup')
    decoded = decode_prefix(conv)
    m,p,updates = conv.serial(*conv.entrance(), PREFIX_SWEEPS)
    need(updates == PREFIX_UPDATES, 'prefix main-row update count')
    trace('prefix-recurrence')
    regenerated = coarsen_prefix(m,p)
    if decoded != regenerated:
        raise ValueError('packed prefix must equal every recomputed tail')
    trace('prefix-verified')
    moments = open_module('moments.py')
    gain = open_module('gain.py')
    energy = open_module('energy.py')
    if energy['THRESHOLD'] != gain['T'] or gain['T'] != F(-699, 2000):
        raise ValueError('threshold agreement')
    rows, initial_energy = moments['compute']()
    cost = energy['energy'](rows, initial_energy, moments, gain)
    initial = energy['initial'](rows, moments, gain)
    trace('energy-bounds')
    if F(cost['energy_upper_1e12'], 10**12) >= F(123275637453, 1000000000000):
        raise ValueError('strict energy ceiling')
    return {'directed_energy': cost, 'directed_initial_majority': initial,
            'prefix_sha256': PREFIX_RAW_SHA256, 'prefix_depth': PREFIX_DEPTH,
            'prefix_sweeps': PREFIX_SWEEPS, 'prefix_main_row_updates': PREFIX_UPDATES,
            'recomputed_pairs': conv.msize, 'arithmetic_backends': sorted(ARITHMETIC_BACKENDS)}

def terminal_evidence():
    trace('terminal-start')
    conv = transfer()
    trace('terminal-setup')
    m, p = decode_prefix(conv)
    m,p,updates = conv.serial(m,p,SUFFIX_SWEEPS)
    need(PREFIX_UPDATES + updates == MAIN_UPDATES, 'total main-row update count')
    trace('suffix-recurrence')
    truth_controls()
    literal_minimum_controls()
    upper = root(conv, m, p, 400, 400)
    if 2 * upper >= Q:
        raise ValueError('terminal majority')
    trace('terminal-root')
    entropy = semantic()
    trace('entropy-verified')
    if F(entropy['coefficient']) != F(103101698978741673457, 20257885502651490000) or entropy['checks'] != 2265:
        raise ValueError('joint score or dual accounting')
    return {
        'threshold': '-699/2000', 'graph_cap': '123275637453/1000000000000', 'coefficient': entropy['coefficient'],
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
