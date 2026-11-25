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

    A = F(13, 50)
    eta = F(163238963879, 500000000000)
    beta = F(31797101911, 500000000000)
    theta = F(1897862607, 50000000000)
    lam = F(19549321283, 1000000000000)
    xi = F(9608353779, 500000000000)
    kappa = F(584461231, 31250000000)
    margin = F(1, 1000000)
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
    p_anchors = [(210620, 210620, 210620, 368140), (218451, 218451, 218451, 344647), (164713, 269171, 164713, 401403), (173518, 273752, 173518, 379212), (158972, 207369, 207369, 426290), (167179, 214487, 214487, 403847)]
    for b, o in product(range(3), range(2)):
        cu, cv = eta if b else F(1), eta if b == 2 else F(1)
        cz = 1 - eta * o
        du, dv = ((F(1, 2), F(1, 2)), (F(0), F(1)), (eta, eta))[b]
        k = F(1)
        flags = xi * (b - o)
        p = [F(n, 10**6) for n in p_anchors[2*b+o]]
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
            query_weight = lam + (kappa if (b, o) == (0, 0) else 0)
            if b == 0:
                terms.extend((query_weight/3, r) for r in (U[u], V[v], X[u^v]))
                if o:
                    for a in range(3):
                        event = p[3] if index == 3 else p[a] if index == a else 1-p[3]-p[a]
                        terms.append((kappa/3, event))
            elif b == 1:
                terms.append((lam,U[u]))
            else:
                terms.extend(((lam/2,U[u]),(lam/2,V[v])))
            if b and o:
                for flag, marginal in ((u,U), (v,V))[:b]:
                    probability = F(1) if not flag else Z[1]/marginal[1] if u*v else 1-Z[1]/marginal[1]
                    terms.append((kappa/b, probability))
            exact_test(f'rank0-{b}-{o}-corner{index}',cz+cu+cv-theta*cost-flags,terms)
        for name, credit, hmax in (
                ('input-U',cu*A,1-A), ('input-V',cv*A,1-A),
                ('output',cz*A,1), ('info',k*A,1)):
            constant = credit - theta*(2-b) - flags - lam*hmax
            if (b, o) == (0, 1):
                exact_test(f'rank0-{b}-{o}-clip-{name}',constant,[(-kappa,F(3))])
            else:
                jmax = hmax if (b, o) == (0, 0) else int(bool(b and o))
                exact_test(f'rank0-{b}-{o}-clip-{name}',constant-kappa*jmax)

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
    def cell(a,b,c,d,depth,o):
        cz=1-eta*o; weight=1-eta-lam/2+kappa*o/2
        if weight<0:raise ValueError('entropy chord sign')
        u,v=(a+b)/2,(c+d)/2
        lo,hi=max(F(0),u+v-1),min(u,v)
        for _ in range(30):
            z=(lo+hi)/2
            if exact_bound(0,[(1,(u-z)*(v-z)/((1-u-v+z)*z)),
                    (cz,(1-z)/z),(kappa*o/2,(u-z)*(v-z)/(z*z))])>0:lo=z
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
                    terms.extend((kappa*q/2,r) for q,r in
                        zip((1-X,X-Z,Z),(1-x,x-z,z)))
            for x,y,t in ((a,b,U),(c,d,V)):
                for q,w in ((x,(y-t)/(y-x)),(y,(t-x)/(y-x))):
                    terms.extend(((-weight*w*q,q),(-weight*w*(1-q),1-q)))
            checks.append(exact_bound(cz+2*eta-xi*(2-o),terms))
        if min(checks)>margin:
            checked.extend((f'rank0-2-{o}-cell',x) for x in checks)
            return
        if depth==10:raise ValueError('unresolved entropy cell')
        for aa,bb in ((a,u),(u,b)):
            for cc,dd in ((c,v),(v,d)):
                if bb<=cc or (aa,bb)==(cc,dd):cell(aa,bb,cc,dd,depth+1,o)
                elif aa>=dd:continue
                else:raise ValueError('dyadic symmetry coverage')
    for o in range(2):cell(edge,1-edge,edge,1-edge,0,o)

    q_anchors = [(602867, 919935), (576543, 916772), (715413, 962750), (680777, 963401)]
    for b, o in product(range(2), range(2)):
        cz = 1 - eta * o
        w = eta if b else F(1)
        offset = beta + xi*(b-o)
        cmax = F(0) if b else F(1,2)
        if b:
            jmax = F(o)
        else:
            jmax = 1+F(o,2)
        exact_test(f'rank1-{b}-{o}-clip-output',cz*A-lam-theta*cmax-kappa*jmax-offset)
        for clip, ref in enumerate(q_anchors[2*b+o]):
            q0 = F(ref,10**6)
            weight = lam + (0 if clip else w)
            constant = cz + (w*A if clip else w)
            if not b:
                weight += kappa * (F(1,2) if o else 1)
                if o:
                    constant -= kappa
            for q in (F(0), F(1)):
                C = q-1 if b else q/2
                terms = [(cz*q/2,q0/2),(cz*(1-q/2),1-q0/2),
                         (weight*q,q0),(weight*(1-q),1-q0)]
                affine = constant-theta*C-offset-(kappa*q if b and o else 0)
                exact_test(f'rank1-{b}-{o}-branch{clip}-end{q}',affine,terms)
    for o in range(2):
        cz = 1-eta*o
        offset = 2*beta-3*theta/4+lam+kappa*(1+F(o,2))-xi*o
        exact_test(f'rank2-{o}-clip',cz*A-offset)
        exact_test(f'rank2-{o}-unclipped',-cz-offset,[(3*cz/4,F(3))])

    saving = beta + lam + kappa + margin
    delta = saving / (3+A)
    cap = F(24600132703,200000000000)
    coefficient = 1 + 1/(2*cap) + delta
    return {'checks':len(checked), 'local_reserve':str(margin),
            'saving':str(saving), 'semantic_increment':str(delta),
            'coefficient':str(coefficient)}


PACKED_PREFIX_FIRST = (
    '{Wp48S^xk9=GL@E0stWa5C8xGsBn&k<HilIz5oD9L3lw>IbP8DoX~!{&eI;iDA_UXyGwDA^8ZDDlo+;@=h9g|!;++?7+_WUyPtAhQ$v!JLjSaewlwuus*UL|6hJ&BvS~KfvH-&d7^$1c$3iiS-k2A)CuJvR2S9Uk2EkT&g}DG=-59dIki+J!kESLW5fpa=fU+10F4astPXSISa!mf|d+vKhYl@Q}Jwa7o0Zn^<G#d{0id(Zt'
    'RDf;5Z=yz0aHD4vq0V@#*a2xWRlEVR;L+uNK>-UQs}LQW)IqU~LI~<pMLP~Be^A9A2oOA&Z>Yd)l4TE)?7ay~PinUPX3?uWTHrufjAp=0?CODCMqRX_CgM}53F8F}GcGF}apT#}t&Pg@_VaP*t_B$~sU%h@zS4lw85iPxZ?$e?yOO-1LBzi({Ot3%<sm*YON$%L06Y9npz@n7k)dVXwjKxjO}2a<C~2C@H9X03kpR3UAh~(6'
    'q}C`XQXwZS3Rnb#F0r(okF1Ehe!)kL18r%5nU$o5vi4(eFmFz3r_AhD%Q+r;sO6u}0I7?Bcj{=yq9^36k4CAJeU3=hWj+(tTR>N-;-NfCfaeKW!580-sdU$YcfDljzV8@$S164KZGgy!Gk40@fml-;3>>((fFz2th~pQmbag)`k-7&JjWt_{>yz5QE#->A`kb3ItuUooA0mofA6kel>fr(3{k+VP+)0*qt{?1O)CgwX`U9Ar'
    'qxWJ}n#6}16VfbZdfN9BLry5OpgS3&dQ=9YS~`2cH-i4lm9_c{ku0XvtokE(w`YQA86!4{qKhs90OTis^NUKDTxESR<sCyfIcYtJ8suh7bx;5Z{6~qwAsn|{=wR$Z^lYrh%hOqjg>gS|d$P5;9H<4$)4Afpei}1h(2g$n-K_t2YCFhA27`q00wb!mst+`DOuZnvX>5n&QCppme!P9+`@^UdSWdh;EMp`)nFq4ozD)ZCLKduN'
    'rnKlU4VMhF!#)eJZQR+o8dXWZxj>U)yYIv*IfS7f{l0}e?YVc3hKwn*dvB_R+wWG<zuC-se&vHl%!}rZ#c%av>xb;)+2OM9S^xkYd7@V2Pps4e`Wpx1X_W`0s_CMTaH%n=?@!~x!iHnimD4S(ne0qtEP!W>(*gw(k4`~%zDQJuFz>oSth3&Pi5DIW6+VsRx77_j7DUuNS>Fd~(J?P$Gb$K-i*CliVXn)eAPR}6$v_l#;s}^j'
    '0scEj*Di(XW@7WE1L>s=qFcJ76uBJCRzjUvP)6IVjNc@}c<$5K2S3hYa30XXEtdh`-@$Fh1sFS^&od|}IIJD`HBKyh+Iib9V)<MFrejc9T#D2$o0$(ZK#_F~Vu&hhza~(KER&Siw;&G^RIo~4)yM`p$&G)_4t~pOfYcZG{Zxc%$?u`fQ1O7@hiBVU#@QsRt>2OIv^C@gWMMo(H-+=hCu3?oCSyA&$|HPms%8kaRf{67$j%hB'
    '48l`1fUEXPCTkY+Ex9E$`4NUGD(-G$1cZ7lc8`|_C-%_=7FLX>GwXjM4%W?@5}{*K#E%Cg=2?K`kUW%ap^I5EQ?#QE{07>z(V{Ds5HDrNP=9=i0~0<`WO<-bi5E~c+_TiE*Dy4D*zE%H1faKEZZ?u%N2g$veg_2&7|Jdp*^dU6Boq9SHoJdQ#qs|=ikfm%v15#xBdbrrrWO=#@7vlV)W)-e);qR3ZbkWv0i^aUOg0FJH7BO_'
    '<rTpi+Q_DyRPsyyj9t!;+H5X8)k&+J+jd{nMpswt`BK>d$WoYBU8X*#Ep<-OD8TOje(P58_l2|Dm4+Y;+}uU|K}f8oew3{S81Eu_O?XR6<nAb!JG<S-Y8BawIk0jBOrpN7X7ZC+&UgOMLKy6;x;Er8r{`#nKQp<%+S5#UH%Le+G2E;Q)kho5>?9C9KK;I(pqW2tJV|86{D0TmyI}l<Y2#aRjIL6Wx~O0qAlk#H!o@7s@&9^h'
    'nwxydzGxnS(yW^*l%g>qhH34i2h-1QwN7fpDB?-N07UR_{t442U}}8s0pgZHNJP1uY_>mhW>8dDNP;wtvXC%rPSd7L<oiB0oi)vyaHL#b6#5%T47Kg`cuK$~fcg>Wo!mVIM(crB=P3gcA4j?GPBAV;kBsfYPUq-y^;xgkS-YRUHtBcVLxEYOsur~G0a3R3zL?Q-&K>s-ev<G+hHoTG!V3OtZ^QlZ)hk@ZEpD=UL}<*3jH#FC'
    'qUR?5Xf(N5WrWx+N)gH#$2(!WldGg=#|z%x)Gd&oUn+-RiB0;p3bBUIlcU>XkbQ*+pjS6Ywuih_nwMo`mVdI%bMKTx3cu1IQQkd|5B%6tQg_SiqvZ3J%6xD2MxY5+wrMo+IX|)>>eQvI1}+AQCEb6l60x}g0P=#XV%KTnwE{o$<hPLicXeNdoeKs!zFGRlH(DA<pk=NM(6`CHT^Nk&IiYnDY&6?KfzvFr`r5NJ;5GWL{&cf9'
    'd|{bcovuvsck#3~DAJYmn$&SOu+^;^ynhngE)5K@{;l$#<hkimX4kzH`XP2ZJTc@r|4OGdCHy}?h8fXpZc_jzNKohA0SOa&dg*@uX)XUhr=U3xeDu@d;F@r^_HiW$yQBIO7<Z$YNBXzX`{E=y%eP1UCUTradgfwo+9cB@(@xOE!%IpiQWQUY95U(&md7b^jH!NPPdO|AVZ089`INlP{;abb1PTKVr$Mv%GVJ!ibHG6IQTl#|'
    'rX_Ni=UA#HxaR$Sg*41__)69h6x47#b^RLd1S}atkJ(w&B(R=M@R%9pm?zW$RxKtEzoe(;MAW6mztc{_#OI%;tXoax`-pxJDgJ^qsz}14gu1)Y&Vk4zd{*-~Z6s`(t%SJFqv#^H>e9af)DIF@%To-2MuB=(cwjeZ2{QOKgyYQ%j6j8#^_DmE9R=3NM)U_Vt-?jiNEnFLG>-Mg>JHzh{z9UI3I?5u;`~{#v(eRX6X|?&=P%z4'
    '8GWPd;E<uNATv@PXM6Fbi-+Lk0nd5;gx}18aX*Xz9yxPGcB`S;T}hXv@cbEi^g4iV6|F%)<n72_nIu_TEY1}#dl(8h02)q__ktbPskzT<Kb`tu97x2=$Xje(kMAQnbpUMl{0KL-s~i1L_rhlah!&`L!@);uybB)y854N7RP-unm#QyiF+OSim%O~J$|%x4d-Gw4d$?N;4GX1bfW#PPBFVP{nZJF*ViWTB=pIc*^6L9Z=Dmk5'
    '7He1xygYhl)%>ydEsG+eXB_PaYj9beG9pNP>~g^978qE`5{!0c@C|cIHleH|6nAZf(}?%0J2c20z^9M2@T>dD$2S(7w)e*8P~Ftyj>6Dxg77{a%oFu|_J@i%zX&o_YGTFg3>jO1_eX1G{k`OU+7OOUIEMOXJuoIa<{1Mf@tshe05@}iN&t@^khk%+GzDh&F|fAY96$z$D1hiOv&hudHgeyyfZGx_f2=EK@I-IKtvRep9S~_5'
    'IM7SSsomNbgbn}`Jr1gf5MBB#vGl$TAF)fudtU3n1mf}iGVO+|-m|D5Nb9n#F`JFpq5Z)Wp$hA>YFy*X$7q8}heK2J&H?`!y!$bVi#1(oCyQ}`h*$MzT2|H~VD9W_9+G;ehbN@#b6j}nJ>FuJ4td|zu8<6r_V0+&u3tRB#eiw+R4BftG}p%Rr{*;7y+DNdEjD7HsXT&=g)QC5B2*X&xpAfED}2YagrY^?BFBi;fhocEJAB)q'
    'vpJT{22J2r8MW44tKCI=_Vd&nu>Mi#)s>*^QI|{x<B|1Do$6j^r?&p2k1h9<DIMSkd!*LCi|ulxz|>0=8(Y_1`gt$;A4u{<OP6P3*Uh5)UC+hXpx*-}A+0}E%bQELKHcNaT<tDk@ovVb|GF(5bu@G|eIckFPeWV;o$YR>ED2|-qQWWqf2Pto+}A{ap*w;}tA%qF6ZX_9gc3?99N1&#IrroU;Lr9Ooa3WTX*R)C4iw<E-GBk#'
    'Ogz2-?vx9AET?wq<ELbCoP&%K9esSbVZH~F_{ECuyk8fDcspv_?}@et-bv7m)5j#K|KlCjW>hf{v%W!WyFG;V^I0CXo?bc1k0{Zn*&cmynTYiq$f3g}vt!N7P-i1^YvF^HY*cyP)j6f-Ij^FZ$>9?4w`1s8vad3UN62aN1xDrb0_8J5Fjpm<wlig^X3S_1$>NfH?UR3m+ox6GJ*_>v*d2rBsPEB1MSZjJ0Vtn7D0*wy;N^3k'
    'e{00B3MB#pGi^p&Qrh+WDdm0j8+se@^joT1Pf;)WB(OMWgWz7M>V5-)oijmY&G1ChRyC9C-#A-m=1;2osNmz_yV=APYra|_IZumtCL-m|<Kk3QW<y942*ak2vHg3A*%gD0L(MoKN6ImapI^gLX-5^j(piRT#~M}>#;vHW4K^yxtaH6>xIe9bYRNxChO|QAqW-5Y-34zec^qabHP~?>&6>Rx4@=O)sK|9NM^WjliC0i3OxEGC'
    '6BN<FZLXC$M*sk@1o5Lre;A%5vZI-#v^A_JlCcl<?3K*-NMe<WtHA}j*FbCjkrrFp81%LMnbL^aHZj$Ux*`BW3w+Zz<+ezRm9?snQ(bFlwz}s7#uu2J`?gCDe28LdapcVA6Xov(g-h#9jg^hrVkmIy`4y#UzP=hpcvq$UTyblLVbQx?Pyk+$z-zzKQjFr|dP><aKu&lHJ77xXY?~gc_sxU!rB0%Qlfiq$O`GVMD@Wo;6w7T|'
    'STklFlmxc2jtaZV0d6<rjCl}SQvb^>2kOAW{D?bZka^bCh%dkt{=$)7P#*^o^6w`+Ft3GakgF@Xhd=i<rbgvA9AyTzmLNgxzk=F!V2y;XI~r?HKyMpN7(6qGE)$FvN(xPye@lJ|SC{S!m^-rXRYZV{;|B|M&IG(IxP`jQImVmLQ-O+_>m@gw2OpsPjo5O?*C)lzr(;uJs~%;7?tg-4fq_=DY1ebPh1{Q33)$hmyY4(vUJock'
    'xW#&H?ppAkJ4^~)?0O~Ojav;P+lEx^-}fG}*Vlq=k`j^KOSRqBfLzjj$VjNGk@x-A6|QZCMO}-W6oq~}qo8x#xC9!BB2&E$Mz~qh0#9ifFWs@|#ogN#P<Dd^@F6sh@(kaV3yQ2^wAWYMkcu$u4R1b9P#NcETxu<6XJS-Kb!Q(R?XJa!M$sY{``V(6O!@`&H;(qi#u>c*(y?%W@ZsK&)yjq;Vei>9Dsu=o3D?0C1i*+~PHw_}'
    'mJ@562g(nL*B@w}p{k?lH3l}S{C_hSXn#Y*sU%2&^vMUL{!{T}OcrkhInY#}u-#)yZCk&c4D6jD6F+L#Yp+FCU{`{IdQ*Q1<D+CKk@6Y?44#G7Sn38HRVHv?^N#vCeY*n}Kfh;BHhY`nDc9`v`!CD1FJU@T3WtRMM~k#&?^;Svk8DmPnoJ3W4Jpf<ja3B4{&*Hlm4B?|9MfmgW>q0X?eSgQH_6b1X1o+Zkx08GKIKp^fabQb'
    '&yWYqLg^4eN8+o=`_xCiG!s8BpLP-o9G2-9Ni`re#|$`E8MD(b<S}iDy%GaN4Zx(;4~uOakBC`LGD83I$ZyI~m~k3|oWtOJ+IMm0OHuOTj2ji=)I;|7&tbw@f5H|ZObxJrXHjy<<|&d7n(yZoJ#5BWL1mBw^5Yu!WPM9tUZWM2y)Gm7LMe@n@|opKfwCr>c8YFy!`r*W(>X?Q_-7KSDRHK%00_@cyBiD4=K-zi0w&ZYz?^&B'
    'n4eB_#Vk0)w}tGwA6Ez(q9lhQ3f;yHPJQD?Sx$h$SCJ>3%N9X8RRm$!$Rq$r(+jRUq_$R8Zm=6|q`SsyYdR{MWfOdEbQ0zRBRw_#yRSm-e}vRvrW6yh%5h|&I|z*0Zsm!~QSE}|Bl!QHUu|qGHVsH*jmbq^_r5k|SXH{??Lz*BUVx4*2@JKVrQd6>G%3VMYf{<b?6BRI1%Qq{v`7WM6~%0vqGrV9rbHU2WlZvtVO~*)7)U<t'
    'k&H$OI?A#%{`9<(29ALm=}ggUXuA~U0u5<BKX9pwX69RHe7LVR;d(N-+MyGD@%fcnvWLb+vgW1Fw2~k)vcy@3{j~nBvZa>k5{*)bA-?d04_{C$)-Kb&D=5;=O)Kr{KD8&+rK!dkadC>RlQFTrr~|}VjF<OlZpe<nyQ|)m__#nc|6|sSCwMyNOGNZ+)+u}K0mu_92sE1bp}K-JET!zR>;L{P>GX;kQ}Z@Fc2<DcTl6dJ_t&pl'
    'g`V}-&RD!yjclGfuVetb$c@w%hv9~Vz2QtyHQZ8wtbZJ>K{0zv8WZcXd$mUMlxpZaENZj2>`S6-0?(`BR##W*c{1FzG{kSZGi1}3{vobsO5Vm%5K`jUI(V89Fkcn28h0m40I`?O1lWeFA~9og8!Qh)Tv2Vu>crT9vw|A)KobCfEWB&J8&$eO^B{5Sgv~#s@0%JopX%$Ov&X@8*ncUu!v&m7uHKPpZR@!iMkMA`>JD>_QPRD{'
    '+)9|i$21(5<mS3sKy_A=S!7vSVMV-Q*L3f#iZcAzchILm*}JyAkejS@bcs@D(w@NTBj1@P;}`p_y<=Jb`ZT4#MTU$f#|#H3P{w%<mK{EY49n*qCN7^wOqw2;;h7s~*E9v9!=764AvXf4y1UW6P@4_Gr`#2|nobUTx{W+I=N%}9NXV+oX;J6g_;IO>sFd*)%NGtbo+9eZ_gb5yjqOoefG%NpSD-W6*09g^>^5Qlz}HUSwr+`Z'
    'Ta*Nm#p<ji20vCJk1uqtY*{VjB9r&95-wN;CIezZ&ng+kj$;hhR5>1{G|nc&Tp^aTSthwxpNQg*Vl23IlaHS-?u3>%<Ta@Kj|*?O07BjK_J+M@JGbmdIHvTv<-@tZxIJC(0o0t~hIG(GOwGGfs9`4wu&dthM<HTsI}6$;W8ah(v3nrUgojmrp4ngl;QSAc-92kFG}6$|wFV!7stlI~gsaC44?+D|c-@4Bv{f`1$7#cMCIp{T'
    'ye~E6te@KDjw68cH{HaB<nAv%0WJy6*p$(_G4qPCRJB0#f(u#9A?7WKqT;9pP@*?_4&F#7xT`0n6OD9Zo3)4{O34}LjWMi5$@XYCNH(4Iw1lc=s+oUI33=j8zqMp|WAq0iRnInFh>vyS4P0#@0u5W3M*X|y9uknvSeRiv9<b@)s%#9vH?5JhxVxvEX3VEqKQAol0FW%3zT2YA%<5Bz3q$*CyWg}9$;VgOF^xiyA3dBkqzT2q'
    'Xu`7&f*Z?&015wd?X(tVrA4_QroKCWJ&5F0seS=Kk>Gl`1L7&7A81jquUrqIXwyJ(_rgM%`<%`X+?YHCbsm$^Bt^D>tL3%S(X=8TW-5rES5S=XJU4EVCs=ywW+q4}pK{3m9-t1Ds%UTMbn0>Qadeh{b2)b9V=|S0uxpqVY@j_*DU!{#xEo&W)DB38LjFnEH48^gXI8JRt;@4a1$O(e#yc+mZ%KUPw$x8wTb5M_POhXqCU^;n'
    'M1`V0`9Ij1Rn77c$;gVTNoRKIN#nvCT0Z%1vn0H-hTP+UB1egkoh>HR3S&^_0hHwRn{kktXt)xE6TDtv(KzSPjfaVF8F5wnKH>fdDa8$jXLQ{WzIksfih^z;b5*>Nz6b94%@U1Qp&vC2NbsEAbt63Ec0kS4;jGU9<wTE(am+tZS^aE(+P*y(N|^r6=^A5O386KK2=P<kvk-8DS7LDpo-Ua1zjQ*E4`xi1&_7=DslWTWkdOy>'
    'dwBVcYvG-tK%ujp^~xVu(072Xa|ayn^)v~O3Ndl_;Wi6|&#C8%z0d|0+v-0G(SWA&28+OjU$6y9UX%AlXS?7hrcg6p5J{`#@!pu*%;C8q22DrX?$W!kgpd39fJFGnH(VtCNWwQof=&x<L>}jNMt$o4#~b#$Fg%0S@iHs%A}dbY=Lg~3y}wW;KeolI0#RpL1WKmX0MZMGPAssgY)7c>l}-f|2?s+E>tVWt?MiXrA<#{?P1>4#'
    '+4v4F3Nf2wjb0&rn}g>KRRiMz$RA<%b^+yb9KkE6xd+;_AlQZiE()LtFpcHA&a8~Iq{Qmy3BTCt(D5hbtIEuVKkJGd(bMfA>~b{kHB+^;qGn;KL*b60MYn;pCUuUCNz=>sEjr`Uf9E1HYeW7+EpSgGzS5W}jQj!y6JUC&$@BI@T+B{ByPb^7o@`35F{Y!cn)E{?Ke|1#jd)TyZAUcySTaoS{daFhUNln?gvEUQ0gc%923`se'
    'PerN5>G4Rh?8jA0x$K^3AWcaK9=r)}!_lxQt%GR0Xf!twXxi{(YP5j>OF{(Jxs~Z1@@7}sz1Lk4Jx~)yjM>XHpScp9YcFgo^$^46G-Q>8-M3x}B60K)O;t6xj>@-r(LqzGsJDTTFgTs%b(=Lav$5@KukFX%a7%G|BATV#EMy!9my$1h;Bhn(e~*BKyCgp_o0q)rwlyYa*WAC^pUZ}Ess-Dc+)J$xhf`H)VtZIn-xOfIZr^yp'
    'o&sV)F}gAAvD$!v6$QgseP(P%9e$sGMC{jT-%y%pr%_49BMD-jPKDej>by+Rs@fdBb41JPlY~2!L7G8{YSv}Oq+G@8D23tkB_a5spI-xV3VgeuD)cGwaaA7q6{9Y(mj#6=iV5&;%87-*d<6ctYPwF~Vzs1}@GP#W{$^xYRhtUQ;v|*-v9N3PTpmyEvh|zqd%{36W*Yh*Axx>UqqKytNfP&2&Bfxa$^H!aQid6Zrp|*SK8UTg'
    '>@dDZ(q!r|1x{cCO=IS^#lAN}$OZU0{G+jLD)D<2f;&AaA5yFdU}VrZJY7u;JT4)ZZ)T&!?1sG%ttkAR{Nbz0>>;%vk*{6#ZempvxQ~7We#Ld4pF{W`V`Gxx&@@sZZXwQ$zs=+l`81hM7i&@t$l`lx@GlHzv7%~ucio6Gis4?geq-e?g@__jR@$#x1L}WEpWBA&U-mg7iB1gmXW>JlzZtd!`V|<su~s1yR|jv)`X8V!@~wgS'
    '%<<-dh`(dh*#oWu(_i>~`pPCG>D*6V?(`MPyimF@@aNWDdpH7HbkQv%9Ec%jb8>@0$f9yOlD~NXShl0QD4*K!4^(nbjlMqVN9By-&+=)l#8X<=Gt?o0ow-555BcMtAL)c^J$a7)!?xlN_T#JJNBr%~M5fBOMU?G;gE1Kc%>3U;V@MK~RG3^9S`62nv|J+Ho14Fbzo+m73;JjnbWux}@yirt8e5j9+bhf^Y%2KG8|79`Y1ua<'
    'Ee-AVypV2vNpXz&@7#_m`t3Uq^(1{%@MMw&L6ac^Wc^6s#mD<4U@p2OE=rbKT{82aUAGUb7Vn`UJdUD~tra2IkAV&@RnWFYmsc^kbR26-GLB}8n^}9T^<Dy>Q`QOPsu2|A;XC*1G;eQIc8?-D#229@EYH9bT<-&DTsvD>)<^Ky>RQc8Vb{6kh3b+_c*(kKySQ`CA9qm%Z>GIon=Qb!qkv&_BWu<~=))@#IAQPd{hlaKinB3{'
    'A)uBDgQ@^yHS9q^hyyDN91WNzqv)BdGY|U_g7y@5?^bgVld&4$NNqP!e3aj+ZpK9f#40em3U2FYL_Z(i9=YDmaA{N4N10LU^O#?}jgPVtoe=E6;iT~@YC%0$qf<(0**#G6{0nsR^%m=B7%_<^k~U*4U!*9$60{$KB#B4GNaCZSC+F`$oLe~4esDLD5eO3P>V0iIjmB8Ueu8+ceu@s6PFb`k(A;_>+F>(c#K|?EG|(**e^`Sq'
    'l=9qV<#b!+U17vkSx5S)dpod+PIfEM4~g9wfgh9z-jNj+QyjwrTOY88=Nzv@usXFAK`M}h5un?_m3<^|AXB4Es=x_!g#j3{WoAK6792edM3ADgA(w?;wBg6v9j07@pMOl2BO)zx((tmH&Gyatr3sRcBi1uh#p^`Qcu#VB`>K27%;f|rs{BGC)o|iBlWU=lJFcz}TB+jTRJ3Djr^hN8y$-DxyJ%6YpILJG;LHZf?Ir-r>p9~`'
    'uXp&rvg67g@m-3+3_|Y80UD)6<xO;833{il>2T&>*mOxvFR*2)cFh@jQE2XU7?xX$dYyPk1pIB7m&nGSQgtIK(G<KoM>Mmvcf&@@o$u2eysvpy4-`_h(F7#@rVb5C@OQ%tcd(hL&%U+ti-}T<Wsa~#cwj(T(SxREpq8(fAy6j2o6PPpvim7G0C#qtQ^Rr^NtKZop9i*V15fjx={c;((*eXrR$sC9espc!=0djmF~8<8PXyXF'
    'sOgg#_cwLuUK*_ii_3~t#Tndz0$68OM(j3MIUs7{ZZDmewr+ZkT&ldF(?*TyrDEW2bbEz97DWyhF~4NZtFJ=KPhW0YSc0}q`_>HbOT%oQ>Qta0C+z<`g2Szgp}1#*uX~3O-6uBb_Zbp~b<p^a4MFzom1mSKNT`6&-yT19s6Ch;ic`M?oL%k!U#vv~xUGFwe(+Z50`Z<^FmamOu*{l+u~tGm@S{~+D-VKS(XR^>23bgl0E{_>'
    ';dBEnA_i6O@!ks%3ywE4QE(+<QDhC!wA1aXz;I@lXKJr#w3#|OBTFU)&sRd0IHU|s(IETt{PU>M?}Ffaotn1{^!$;;lzi*w1pc7l#)@w6FyxD<tgm8M&`GqwFoB#EC;K}>7?7;a$H(5Q-1r;a{4uvZ${AV0ID-|q7d2UDVPC@?=FKI8!TK?FsvtVdq6i$g$U|{PG{XrU^oR%SbGcs&W-9-r3rgTKH_*89Ow`0pYxkPB5jB~1'
    'U{=oi0K2i!(tY{vE6eFDJ<dJKlz5E)FRQ1%`a<;JSd})lG#YUjh8aeVI!}Q$?hRE8ky)Vr_;uJqENtr+e1nwsh0~GF5wO1JSrDM6(>=P1)s&j0=p+=x{#s~9(=H;tJ%vsnadJR>Ah^L_*u5X{?^f$;?;T*xndA9~pJx<KmHy^X4+5X0WZt%?yJ9n{LEJ@i1Hf#2`8<?$<EcT<IBNtxww;*T+f}xD7GdqpE@%<_7{F7|$<6A~'
    '2@k%kayx=iDQR%Bf7;6~=lmE2_z#P7wa}W8Ne0!x7w62})-5<LUAHzcl3CfUb+HUR5`q$icu1zgq;}5*)1^U7ek7Kp1Q7*xO!gj;^g;E%)dr}3X`-{Td28S^&6~_b{&o*OwiAMG(G_mI``z_y0UhLnm*;-SfB5jIqyD^%^6J%LGWv8pxR{K1nyzR(wOYThcEnMHox{?~w5q?0(p}4(el>=KzD8d9#;b70l1jfQo{Y#osO}K('
    'a$P;G!f`T$u+m>DfqoN6ars@NBYBSrO*w>?SE6@5DrARq_jw2)fE&t(LJkf}hU#T<G^EudEG;QzHk?KI_*=j@t0%snSmLr81KYWkkDyd`tfyZJ#)$VV<HL|`x{~pPuNS2pLP<9`e?F`j5o4DrJ5Tp?2_JZ>QOO<cjeqlXSQ;&DXFcgk0T5mfr0QJKInJ4>n-4EYWgni3XOw!CEBp-cI?D(=TpF>?9inre$(6XNI4D=@m|{1~'
    'LSkQTfNmZ_B?Tc4))j#*JXWVG6pvPO4F8LG3tOhplJ$EfTK=M?QM{5VWti8YH(3H|hDedvsk9SJL0V95^m`2Jb|t~apxth2SUDr3w`d$}IO?%qz^y5!k395;Hal3Z&XLVg495RbDBsI1>=cax;c9{clCj9w_U0t5r);61T6PW>7I#D~%EN`61dMg-6GAWL&EVZx_&R36|Bg;{s;08<*}DO;fDyc9G!Q-K#IfmrSwNQJ<$Qhs'
    '9<NriDb`wXWLQBX97#iL%`VA3>e(;1b*&k~G(X*ZiC`Fbm=hk!mBI>1SF-|xb7KtoLBFKUfj-#Dp`H;ICwGM^_!~Q(cq#~LC*4ALpc3^+l!p}Olv7cDx|wwb{nu>AAP}TCo=Q&dz74nC$3Z7Tl21`6>c9M($C$vJH*L)JM0RfesP5R)T2k_1vx)PD4zXv)t);3OXkEC6@=ffWPRm9WV7a}<i=NY|6jsVc+6SL=@zJ9nLP=L>'
    'u#lU1<}1u{l9j2fK(pr)b?{Mrp`E4{m*D1Qa%`1&EhO+EvefF~z+}RnX?9~m*<fj3dmfh1YyhpBOXt`UaIo<Afw|6l%Gu6>gQ8!X)^s=D=5u-1jj{+%L(jc?9)5?oiq4?YSDZ@q6(XuGhd)$b1&Z*Y&kqNxjHBu|$~4|}yfHQ^KJ0d!8|GccI&1pc5cp*Z=e2V-o*L;$ty+k`D40@=p+54w4Q4QxmgYdTi;C6`p>dj+dO5nA'
    '{$rQSkh7{>og0y!KF%gsUn*h}9=bh3Xg3AkQ~!WHwpvE9Y=bl3i5OOtZ5JhO@NTH^|L;CZ`OC};AUyZ}QLnO<S-C<w9jYL-m)MFF5N38L;1G7#G64!oH}FD!xnryMC5wWzcSmgXNdQprxu(?#Q|M@iInRGdIZL`hEeY2>a!#RQ+&d=UARq^V$iS5)9UOJo|5^ol8h6P04Od7?<{}HWD*1Amx2NY;wJDA-vt-7e-`|1SB%Of3'
    '4%MDcB{9=rd4IabX{rNxm(cO9^aWx9^69sPe6$Gf%N&ru&DTMrZ-kOxY-;#z<g`jCFY;34S9;PGf9Co2aSPlpzylaZNz0=@At87r_~o=ySyv+)<Caq%km#ZrGa(9?tgktk94qVnPBF3i$5-Cd_@r;mQM(_eAYmv($IU3}{@!L|(E}k$dyT3|U~0Mk<sm&h>YtvMKWtek1F(Y8SC|H(MpeZt)>Vji#2h`@HjtoW;1ev<hDU7='
    'yAR9_&`~BYTUB<)Z!pL_5yoUVhvJ`=2a=KQXO~tXEYE3L;7^gvu_k2EuG2(aT18#r(}P6ahUZMjzOSi}%`Jq{`^|1LCEft0eA``^dN3tKrz{I5pQHOafLzErTMko%_%!cBMA>1D;6ROCQkxk%__hk<oWEl?_tER!%Yo;VpFdSK<^DD${gq~(Y@J47hp6-qS$BAvnc^-YP|bmBAxY{v54bn(zco8TPa0RABDo<OXdYiE7T<)L'
    'DJVM7qr$a0n?B@s|7nu$#x&L$UHB_ay>ZG{#)Ei}C2AZIGnc1D!!D8jU>68X<?`AA87Z#<uyHoN-N9nbZj<4^p+DT83^WvT%^3C?5(#Wyj>0q)Yaw&j>Jg`HJBkb9*#3~I8?Oq!T1Pd;7hhAu)&<Z7Wj104XlDHe*8*6Mu_!)c#f(uE?k}h|z}MPaIfc5c*$W=725=u;!uj*;qW(1#IwkX*GEff~i1j^TIaeVI9~*6Iz-J9b'
    '#lL+N=S@9kyKosaUd7X~K`<P8`obwdO_vT*wtIU(=P#F`pxPs4I9%yYQr6`JmGP#PK^fgk9Xn4qF!%kG<c(q8XMgo?(`W$^1fBYp8okmv^$093mgVF6y0g@o^#@nH7)dvLz=#V~<NavncQh8DWYZLWpmu3KTsz3nK7dYkKW^LGBR8hqFT)U^q++oH&;~~4o7+kKq||*&xvg9K#aEA?7L)R32svyEzbSHva%*0DgAV0U6-29}'
    'j+r)>=lGmVj>?iS(n>KHUez?ThP7?2#r>d<7;N-UHZvpbwp4~kIXorjtU0+~nQP-IU?m+E(Z-e!sfGS!?6s=X>hl&Sbx*})2D1MRoWojlUrf!vXzmn3sDI=&wao=Rh{lpDSyk<^p}RXM(`UC?_*BlvvUjW)4wrzoRr!({FBK}KqiH*O2-$5p8?m8>#z|%vU0#V(Ia^4C4VbQ^Hw$g_n==(aK~9@~FE0SnQzP4mvlb`trQmEE'
    '45FZ_X~KO1f+rE>;wa9%reBw5%6Y+s)g;SYX_XheGctd!cS#P?+n1Whhmw{ttI8p*VVMV3=ME=}{8Y@ED|)L&Z7oI<Osq{9>u%l6zzmS8Q<_fWBqBNh+k|#<yd-uT&d@Nrunt37hs@e!;<B9pWmH@Gg>%lzpRQ?V;kcpp6|^zws4fcVT(>HaoR+8N-)u-uyLc$Nf9QI9WVP9NC}IKMy+Sbs1XdCy?K|WoOkW>8Bwwmdxe?QU'
    '5nTcamfl?>h^KNoQ(;$-4d@Sf=NcAs=s!w7+w4$&Yeo?@=eRQC1Xl2C?6)@V0a{KGoLK!-qkjmKVzCj0cR;Gvgnk?QY+(I$k_@V1KqdQD)}(h|<=9kCN<p^+ykQ%?QQJfzP@2j#3@2Bqz4Q@sjhMGp9eNQ9l--U8mvuu_o$`bQu3~-OsoFE)*R{iFAGZ*uAW5=*pWoq=(P?mxBpPcW<ruC!hvT3e&gYsVppwnUX0%5`0_g1z'
    'v)fH=$U^M$H%_k!_nn9M>`%n^%6JdhFuuE30aC@p<>+4>zyraV#Ol#!yg$rl6yFftqe)0OOiACJzqPB3=gm0{{!hj;KP*uG-e+VQ?uLyyt5xs$-G?;e+&`%FuK<OBi`i6Kew-bbgUqwZ>I|99zByanH{;zJK>I`z+4Zl+LE|<d%GkQD+4#HG2}h!=E+KVbfb6~$zbz^oPB{H1nC)mlljyC!4>PRsu0?sUtR^#p9mwwB5PI-n'
    '5nB2ASC;6|1L0%Uelt}DA<lDtU+N?$u#r}#YGR)|OS%Vk5vD-|%f;@HNYTph%=%O+8Coc^t;r)TOyQoc!sQt9#=a()6q;84hI^!LOXoy!VlP9Ugjt^bdjE9#$7J-_Yl4tHsYb5V=j!~9N`Gtj1d{LDWpHB}XsYouGbQz3=9;^&9dA&S*uxV+U3>&{4zHoqXCDk>Zn~O=tIz->L+@2kO0b@UM79KvtURNyK%9<B+IM}c0bMg^'
    'u{=_5ET(JE=)<)F<}$N~>n@&u`^$7Pblis4hn&K_kBZywB#8dJHWSybbYYf&rc<^%56+Ye49d`jgIm&7WS=o1j6)Sqnz*mlV}kBs&@4UbyxF#iTI(b`?Sm5Z?~|aE4CJ*OZA+2xcwNKsW>t%yvEsFUYMQNIaS2XbK){hxA>q6J4*JuKAnT+&E6_CR5FdYwg{r3xw^+XgY^ps$*?E{x7D}=?oepaapA@n%j-J2ZNu~;J*33t0'
    '(Bd3Be%Nv~oX5%J19lSD@ZkfUh=x&=*quf)hDFb)y;=(_9%FGeFd~Qd5*`5Df}_Q3_;0Bi-UaJrDrH3hk9|MbhhMuxpB|o3gvSC6H+VSybGi@#PaA31^?s<u>2X3X5m)>}X9dTIX>9fE+{g<Mq7xH4nQ@DUvLJm#S^<O`IrAlini!7Z2!Ua!tzazjnlEwn`Z;H|l{)QxV()zbD3M%ZsAg>BN#)t`1Y_a*(&SimY+@(3a^zUf'
    'q$4_jFdSSMSdGr=r6UnDhl4il+Vug;f>R?cs?lGpla*cCfzd3RmNDIlWtTeLe-(j$%{7F;y^%N6aM6QzhfU1$c;E#f3$*U1zMDOU*Wlc-RQiH9mEIa<{_Ew#>>wan_>dwcUSg5RETc|$S(K)Xlf}zt0gQ9q;0&pwQ(8YzWd4l9N3%`@1|XX~W&-lP6(@WkY3{QByaV40nMgJIoX$^h8Csl#fxaowS%^+O1u}mRTMvY+75ikB'
    'FM{!<E6Q=RCni!Z5}uik%0`GYdXY+P_w5dCuY)$)uPw+eG%Tn*!!kw{h7WY9rsClBHR}WY;>yu*a)G08wBso48Y7d3kW~U9<j|rXgOV&j!9Xw>^(aFh@oWHvRHzCo;OS)h_~*G9i?MvFiqf-@HLn1~%#}^#gPO~iKvrgh5U@__4`F2Mml&vb8+{ru`62-U1mC4<5V9$ENpni&>?FvzW$=nf9P8i}*p#p1O+b<KevmqDm2b@r'
    '6AIu))c#IrCza~E*xOxo#lD_y?&7NNY6%AI$Hipd3M6#>R3iO_%K1NNqUDxPL}GSzH1&Kb%Ki;m8E$?RBiW>P(KsH4Dqpn(dJ<wfcA~vl?g-?+(OK!~yAC>H&yR6><(DPsW^sRgGq$WX>XvKQN9KZjPUTw15|}$0VvLFy@@!4;v33J4XL-)17Q;EC6CWZ<mX3B(SLT+9Hgl_>O-OLyV2cKV$rl*NInK&ZIJ1r*&q}Q<uq!eM'
    ';cMSGzai?Fc(Q*GC{DK&fJv?$m5l@qb9o=p(&`_vh?XaqNh<q$Z3^g+pV3EY`^jc_A~<pbSUpL?e^CF+nroZ49)KYN6Q-7PTM=b_lAb^+9Gp@QTpLMG;Hn|D@BOV+l*I{T3fRZ|+?ou8StMF$_J@R?O-R~sXVEKCw+tp_@(o!|+C7c)n*@yLmb<{Q+<hhi01u@Do{$&CTi_4aF7q)CN@d%Aw7Bc4KEG-;!3^$HBhcWZ0<D#S'
    '!?j{8(<MSjM8iye?9oHHJRn5%FTU7ttMpKDg`M-5-4(n=W@*T^JoJdB)vrR}&<z)ROfRI86v1_f63L&_voPBLx_S!AYms~e&_|uk2;UYEmD-FlqNT@n(T>*1God&vcr5%&d)03(9U*rkcK9fYrm;QQTae;jl5CqKWYR`rV{~3R=-!OyXqyeSJv1q>p6sW&<wkj4){nu21ga47Z1^ls*mG>>tU8&le9pS!6ft{rU2^R)ByzBo'
    'SAnRSAlf<`SDbN}&uOEp<q$hOcyeQ6If&**(Njt=X&;zxamPx5SItG4HMPgf3$}GJ@qa+<27yT!To>C#o7ql?Jtgs6UyNaO`jpwhTY8w?2~UgqRje&@NT8495qNlSGwk!_+IX%Dtv&y~o!QTp#ZOKu<4i9#j}SBMr!iD;NL0LP0pgn#lPa5BbF*(L!<+iJ)80Ox*4F>@R?OTOambOd+$=4NC!5)>_wy79Uecd{w=O~HpgUW{'
    'WZVK?O{F5^THCQFSpGwHMD^;$YMO;b&tt_brtBg9xDzwP>p>OH&BN@G=d75;kXR8<%ccp}fk}<GByA>2^t&=;14)8-Z6E#fm0{ddLbH@^g?VimBl$Z{#lts(0fDS%k((%huwgvsDFZcB$nE?^n51W%*by7x|KHJ2gitDXzdm`^Cs~(nCuoa1cXqM|JO}<)X^hfy7iXYsZF-1#pzBM#6MU$&A$Tad#1?UNtsxla(yhE5vP{A7'
    '%~dx_K;eCbZS|<aS;{^hx#T%L9hDrB1q;z_1{%=pfPJ1q16f6buVR{xj+Bmdj$U$|iYuixcJp?OA*atC#rFW1v|gdrYLV?K*JM<rfDWV7L%~#2;Ji9Iape4FvJBUe-0nnn?KGdPLPB+MVBuBzz~IR@S2tc*8eXUA3citZN!ljXqw%hx<YI*N$Id5c*zpnnf?ua}kK!T^0Yo;}j*BO}u|n{#gMBW^tU&aDp!zr=>Mh2YoQJx&'
    'OK-<Byd8K7PO8(}74;s=;Vn1^l68>VJcP!EUFftYI-_X#OxrltPM5~0T4j~Dit^m%B*=`k6{=Sabgt1KJG!4wgn+&VQ90nRGEMIkw0#U5ZdgXibbF!F`Q{p3Q^|Nr^NM8Tm42)S%TAv<0f=7xkRK35SDOOq3U>pJDv;$}J(nFFLp41+pMlKIgiBMrnfu;xIx9?jbt$m0L#|ptzJJEicVqcOC4+UA5O~Yd9T7WECybQq$znNM'
    'r^=hOL)76o()Ymh%j|4Q7YTksYQ+*%NV64=ZcX(d5&`o}oym8&6w}(PlPeCi__>qa%{6c=x;CQ83-E2cBYt`L2a4ZSu73Z7i@Q+<(vTVvN!$)$XpYecPo4nt(SCqoF;$@;K7~AM#+FMlnG-hzChX?pSxwy48j@NnZ6JD5uS5^HNC9!lP`BVXP({Wg3E55uwqJcX<iRDpp9AFHgL|N1W{Mx#n>P84ye}FsTeuQj5qedS>!0Rk'
    'kx=j_fwO&p!>>7R+68J{drp@H&R7wuuPmc4IA1OlpD`|d)X%u;kOMba;RK<t&G|FgxI+xa$B1$FsfEK*?zI|x`Yw+b<3P-`SMf_ya0Ne=Qe>i0a4HQD6~(XZ1^k15(%UFN#xgM=|JDzdyq1+VhBs$^CwsOWH(J$vnGEr8*jF05_qzUJ0aioqYcw2E=T15Z$_16Ea`L6cGSh?fY(L&v_e1hykKSC{bA<0|_BwJU_8B&N9Xc)n'
    '>xBqAwb;ob7;dAyLW~XKWeZ$4hf!f6g1<hGkBnf5K0>I)*hu#uf1y<Cc!DK<%Xhj$pAB}IHx&HyAZAv+W4tLjY(&S2y~LT|Q%_;e992J}Nv}RUKNI_11J3KwKHFEA{#&rc63Y?<htlaFN;XdWW@IeWXDl-~^lm&o8jgqfSbut3W5P4JO$9Uho0ixJ;O5CC3p&<@uGezXAI~CIm_D>md$u$9wRpbXjaZ&C;NMWG7WUqUd4-}Y'
    '_DCH`VijEOL|>-M|3K?2jT3@&-t>x7=0i5deJ~=vBN20UN&S-tG><8=2gAdt_`C#23=3;8s#DTz;ayyo`fsmZI?cdbbOEbab+F6QLd_oBGBJGhsEM|XsET7MocaLyz<CDbQWz9~rHi^q>ZI=6mx<#7@>PdmXxZ$G%H&SqYSoCCWOSkxf<#1#Z)@t{2aAm7eGk={m61w16grs6SI?QiD;%{Uv_AiuRO-_^r>}1c_7JXml`DOD'
    'fWoA#1k`yCSK%&^YY2C8AQ9~upRVpax&88P8VSo16)Hf8?Ts!nc%~E3)5d^p{zRi2Dv9{3RN3)h#NB4@>qxwTReKhl0%ku#S6L#?3UGNu5ZO@I=ZI8)>hMo=j`9K5d$T-?uHV`vh~@2+D~}f@J#o%Cm44N(V5sP)_~$C-{7h~)%Tq4x3&2s-G^0_YMF}6uv$A+IhWy4Cr#!Vq0{=-WsYpBL7bUMt0H%Hw4tSA*21j$Kb9J^f'
    'JxqQq<2;X9KO{AL*zg8}p^4ckPueI;=~bMlqC^d-4`KFjJrzhlLGkRok0!)denuJm0t{o68&aV|usOR?K?sI0xti3pGMmeTcV2TTo;~=uexE#wcMj_Nvp|AeIdF^z>`BHo;h<(ftV44T601IEvVt7AJVw~-B}rU=*1$tj-N?ETR~U=0`Phr1p0VG>iXU=sabSYM7t|tLYXbfRO;9*n!9lp>qSaN9MVq9zrXV3&pSf`0RN?51'
    'qjmKQif@ze+mPzKBka0z7`Z9P$ViY6?V5CS5N$%|Ja{QBUo!QbF7j9Ni#m24-Z4t8lAxT(Z)`@#OjmRqVuN9LT<)@f@H{4}n}O7i65>DqpByPJ1~wn6btfh%m_L>(n-F>)7>PW67&x6ePZG(L{Rp@4H{$tEs?}UFOp|x&GfB_+i{w7UNe7UnN;h#$%dRS?NOs~P|JHup-Q<OoplPFNt`Bn%K`!-@468A0yN64hP^nmJpupU<'
    'PpJQkK(>-+(#DAUDNQm#6AfN7&+%sdTvCLZh8JDn_Z`>+ZXT#sm9akFjE+#-?7Vf|VNUE@_W$G|u5q31cS2oa^&Rs%G=y!FMl7czN^oICZDtAnnf|tf3Pu9rOPG?wp+Q=n8rOstb^ju$j3P6w_Z*5C4z0Fh-!j)b5)bl;Y8zozH0o!MA8WzkyvWpwl_1`G>VU>!@_N2V7#Gk~T6rv>AA~C7II6BbowK7}XeYWa|JwKM_7?^<'
    'TFN{m<^-&y`6E1Ze9RNfbD-1R9G1<NC%btKbLO7iGXE{<E(D#@cX;V&GkGVqfusF?>+Rd>h#e4VMy9p_2j;Bq5ObalJkLx3!Hb1H$0ry1&RfuY!j;S&zwq`4m$jKSTwargZC0}$78`YrjH5c_dOhEIKLjm|dy${kvMbrTi;-SPZpwr5HEK{JyIECv8c=qsk%dz!jgBgv+gzK)U(X0C&uhg8l39tTPvN<X3M^<RmBg&?6@x%='
    '-Z86Uw)lek<28ywcC&juhSkK7VAmf>^8vc}LE2IJiz@OTm|c3n??+7Vmh0=@1e4KEhzir2wE7+5`wi3nWhSVn*<a16htqeF?o~sz^k~h&1rsnTH1o$utg}?1pe%S9K_u(Sh8w$e=wj>}x-P;Q+rE_<7Lk7)98u%^6JcX<8}KGQtDc0Zy3A9^zpR=Oa*-4OO1&orgcGSx09dlB=A9mkMLx{>Ic=9Dgq5%mMMl^B@A-CMU2DIp'
    'UM8%$(rxdo`$34i*O9wQVjshunUBiqIYLYxQuXkY-PqiR6CynJp9U&?N3urxnpbYRJNh;MAeyea19`{cX1wMDNBm;8nPVis`rB-t*S{W=tgF2y!gZ{&`DQ8lR`l+Igd3*W$={Tq#Sh-M+$nx_?22d38qh~k3Rat(W*lEO`-1=)U7G{*m%4MUWAn+xx)5XdqT_(@MxM|p%E4bkV38sfTsGor#}&R^W^nG}0Tr!)Lg6Gga@0}('
    '4NTv7vOz%KoyTE?xBGir79#))n~g*I?yE6@8V^S(hwxTDn_|>K;+z%yq`%APNoVM*_3cLCGXc3*SOu52K{gRY8u>O#qY!Cnx_{MH$!ceEQ8oC9M8mJq)u@LW7x$iV-8!_2Q=*)>>VJyiq+V~xpNkFqf<>!gqFzx1@s*s=?51<#6bOfg01Xl=th;2o6$=(9_<z#tV4|^@Sw`Mhf@)rXtqjv9SPIRW?wg8E-t19ab?!t@k5N2q'
    'QA8^@V)=PF{Gxe6u<pUh=Sp%k*U@}!KK(p2muuo{e4}U@yC8qHlf66kF_ROjsZk~BYHj<D*gUk<v6MkaNS==)7oUKH*83vL!?6%<5}t90b1Ey1aRmKUaP`qCK@nPa_!-)O@bnCB%#k3krr@WBA8X5!)Pt&Ual_7KTU>q9VP=5g%VY+K+t8*GO_CmIelqSS*Q;P5KJ6|gF3O@)-}aDL#3-2gU8M7AqExOyN6{@OghVx^vBo7z'
    'b%5md4w37rMp;IRbResR1@<S3t1DF#36X4`wir%bRJUO#GwQG)tLI2P*FMn@ZXFpNM_e~_DbBjf@N$2}bb}88a&%JFzr2dU&c?*X=Qjqun?&;MRRPs6$|cl%*H-#Rzgd~bl2ejut8s>y6W*fXuE`AhndmJU%fJNihR}RknpC)}B0=GIZ>yO32|wREoOqf$sw}H|cP&huDtnFk$da3mW}ezOez?4O#Asqrb}YCnu}DEz!58;2'
    '<jt0HOBI|v=dKtmO-2Mi9J_$&Q~dHH!LaPZtg|gU!Kc5{6^VoUz9ZwveSQf)OsfV+Y_IR&tf`52$wD$e;^O@ax-IbW?cWk|cH+OC!GkM!=4F3=DVw{w(kRc(Ceq!p3><gbmNUf`7L$wTw7K#lAD&1k_c+XJ#QZ04{@cN}kxIrd6L0)(!d{&tw}{fQCyiPJqoL?d<4Uq+rxcoBHXg;;4{yqY_!0%E%0UY$R+W4+t}LE}x5rxf'
    '`lL|4$^5p!WFFaGfvC;0I)stPH=3j6I%Y+2eZa{_eD_tsmhzd<H0VRV;_xTKh>xV+kx}A*E`%U}lV}P}Fsy5*R8{ynjC;xcTb|8vB&|%dglg1WfSMhUR;t@JhK$bx{mEOT@!Y5Z|Kcu##s;RUcLuG-iHMe@N>Hjvkz)6c5VH-fUMlHvSW?99Ey&~y{VSRIg|Qx!zAlM+vfbNTyW3%t<fhlQ)-i1&poTbz?XdW1s^3E4e&SO~'
    'S4l&K%xpeLvc+G-HAGq!%4Hh_%?=BA6{vw@Cg0fy6z5OKwqDsDI?qk<+rgCZQ2f^cwRaTfy_oSJk(%3=syn`GUa0w3=S8xCu(p>Ejeb~)%JY#0U7|Mj{xfXaff(%VZVw-GDC4gmb}2EW=a4K%#58~)MlfdoHm~!+3R$AFK1>C5d=#dI-_HbU(r9lyG2cm*bU@Lnw9iE#)cHqxlmVHo-!3)bh9AS<ENnl|X<CWQ+;nS451xxl'
    'TX!LR5N?X72yifu0X}l2E<6gvl~95^=b8~zXbA!kROFvV0_`_1FvYOqfLvT%7QGt@u7I4?;O!Byek?Y%fca-#x(AhAIwB?5k_HQ~I6tN139&=jb&!PASuRiEToB!$>y08sYkXM0C*^N2E+bT3JR`oSUN6Yq-hwK#>X105<$kfCMwJ>KYx{Z}G<%9UdoBKkn93Ot4q;hiTY^7s_DP?CG^VjtyDdvBnTM?pNkO+ab4u-4Sd*$-'
    'VK)0#@k8_*$7%;x-3qC^E@lVFB2{Z1F>k?T=`eSb@s47!B{*L;WxSgJzzUjdv|<6&v=29mN=qNQmK&>S;z%l?Uyb<J+v?)4JY*pd0J;8XXxcfDf<YZ%Mp=LB=3;c(Lbi;!cc*nA<HmI@#?cFZ;)_$9_rx<u_s_$jpfXhY6|MaZ7BOj=vYo=z!HXgXuMLF>X_8Z9su5w+#t21JGiGBVDE;}s9)von(qcAXO;J_Rcu?75p4q}?'
    'HcOdS>U63t*ZfPy;P>3(t)pQ!PPd$bVQ$9OJrZDj=tlP9<k@CSoV83T=;jXhD)QEaTYKAQ_p#f8mt#W7wxPI!iahxg_R7v!DFOqw1q86LA1RC@`|GIxxw=w}=dX;I23OsWn5B8S4T`RG^QS(qP-z?{r;)S9;~q{KPng$`lQR%e*(a={FNb9vdeYYQGf5&>qdLsNs$6OS!-BrM4TgAI@@SN<Q9#9Ulga|$8YpE81Z*`Z*7FX$'
    'l`ua5RaK>%og6&;MWab<<bu7kTK*2m;n&I}+twL)PwtHkv`J8{#rt7SLoB(t9ob<a9qOxY6tnv%H;L_?aHdl-C(LNeLjR)7*4gy~`^}1YU%^ugnolILb-CSgSFqD_F&tKMmT8qZFH$zSnQIR89Hhvs@Z^ioim47C-*u4)xdwYRvV&X|Ag<5&8?NcM9SN<1?BNI)DtChX7x`=%q>|GdOE+Vv&9y!iNzg3`B<hpnKg7DtME~~e'
    '?GbgVWLVV7uX?<|yp%)uqBKC{+8Jr<s83?49}fTGp?wOPw^a4P1EY`}WK_t2)V;{A)Rkg8Z`@)saZNsefPJ1465xOZDHM|p0^83uiRT4|a#jq}ef~VI%V@&tA{R_Tbb39TuYmpMB=3U$;8jAqO8KU6kE1WWQ_4o9{_-bf2eNA~2bo`o+0<D(a!M!Yst6?sWpZ1NZT*5^&NZBqU$MLs6v#!$nZ#oH$k{_5)L(MKR^M>GxV;yU'
    'fGyL7MK^u$B5lljoSYPC6stt&{dd{o6J~<5>&&M-v`@u;O>u--*>0;gej^2V70Y8}w1#bo4bf9fe5TDdWyr|fkeV5njzBLoQG9BFPI?xR#_LvV=B~C!m+kGC(d25ZAd*3dr;7lkox{$4WWcxrL6JX^_9vb_SyBzr$L&AaKqV72e(C2uZ|y>bf7be$lyigCEYLm@F*f$TrQto12m%J*`8l=7)n7KTGlG(A=9q2R&MRKHUI-2O'
    'jYZ_;@yars8=xL{UXhTP=xv=}0~4sMIVtO1pToIR6zod+5qcw}{NLa2CWYb1fH&qLnr2t0@CxNh-J00)SBx4MSj^V--4YlG9j<u$`hw{Za0D-KC#M=w{QPd8efa-1CXP}k36X#6=I_bke!*x5QA{?iYcc)P-2tLam+UcBe_XS!i6UyVMRi!@^y*~PNwfd>8v%9KJNib-W+KPqcg4h3I1WlgMbwDpe&lV4UlwAAOvb1$3_tz)'
    '+*cR$Yj>3b#<fnHKbXLE^Z3FkyW!0byYX=WD|_UNt85k#1@|w&X+>T3@$j%83sXE=H4nvt;Zx`Ti3n1N#z?c;kIu({@4zf%uG1{7%jY3&sX~F_71l<{QoR#PeuxUbn{5qw^N50RfVN4)yv<yEyx1PgQek$zrv)ieCS`J%1&y8O|K?(Qk_)}n(vr;d44`TQg^8s7(=x{oUw_B2H;j(ko2yX0bMGY#nG<6D>$bGR1n9bCu=Io^'
    '>dEJZ_O{n1sQ`uB3-Y-Q%WZk>D{N6b?S_7mse>^DB=cix80(ke9_u?7ovq?tkn(K~2~U0>Y>|@#(;_(>q@sv9)+x9CLQoT(Yf9m4=yxBp8|IqcA5op>ggbBanTa3Yeb06fqgla+qZ|0%=rwC%Gh};1>ysFJ9*Vi$WsBpJ%v2(8w)=vc63d%@mud9a&{31*>KCRE)t4GtolZk!ZCB>ug+*|zMVTacI6qdY4Ac<JhaP2741kjQ'
    'n%q#jEU*qRfKk>z9aV1<E#@x{Z1I!b0<eVUsQmmilCX3|N@U5+Pkt;nr^<L<rpl$+CAzW;vRYHcc!Xo_EIL<XEEkex%uwm+66^b1nA)~SRReu81Yr&;WGWQTx<TEKbPr#;TH1gJQQtw=tK4V9Y6P1__)}|tOwFXQZfs}8wI3YBSl8;U-YNE9xYa2FX}IzdM(>Y3!Z}W7DLolYp<gUjRJY^A!Qq7erE+v_zUdaNv5m1YulHzx'
    '==fMDwuFb#UB9+CjPgO(6#EiXIxVDyubDIRXA@)K>l8T9`wM>($vK;tmTa1DyA&ag#thz~K^PzPvpnISyM$B!-j1wEVna>0FUkbiAJ&~!q^s>Vc0)(QA63W`LbWA5UD9&4`=+)pAp^<-_nJ<Ujj`CHs9yiE8S2!>Dio*UkdT!OWHy*r^_wZrn3ICR<Cf9|dZQ~GuC>JYmC%KAwj(pI)p_tQg7Rf`<l)kDOZ8<p7u!3!mwm!`'
    'C?vq$+aI8bB|-+%uHMX!%P8#%mnOD)Z9e>Q;}m?gCK<a^RKWa*E_p8KHYJU2gdt=S3}`GAksfx_>|C?+D4z7Cdt**7{akk{+IoxqkHol4BOEtILC$cbKCuq_%U(Ig`L(BPE2JTS(iW_6L|<74g-xOt<-l&wA$i0vshopJ7UyfVUMB|$I6A%8=hdNkhy%y<a`*5{E^p8{K;4uMp0)g#MsRDK`?TLRKD~7jFLG}isykq-^+maa'
    'v*X*D^eO3*8|h?{ZC1c<H#qfym?7ICnt2DRuEuN8_OIGgF6^2FksnxYnSNssorJ6$7bht-Zzq-{lGGw^jU_UyWT#swakL-1wi_>x1hrp%mWt&iTV{+BF+aH*vTI*VA&}<QU|D}jzOg~_A(nh_53KADd=IyQ@YGvJU5-j7%L0%T5lLg5iEKpbSX^Y40NIB+!ql3+()W&~ccbyH?#qO&K0^jrgoU>w@jm6lX2K1jd=u;R2~@pj'
    'FFQ{Zd250?HzIb2J&Egb<}am7s5f?~9t-Cgg~9Asl9iOrE!X}76F38qZXHF>Ag6BMP8*p?`F`ZM?MHPsTt8O5d*@k7nloZyy@B4Zr*l(%P&R$*15RU;K8W)l=1v_L5rGX%2piY+7TDTKtP(eBOsyVmMvy})lXHw0xo3r}2~M4}sKNyd%`j#o*wkS3CxoqH&LNpAJAmDxlLWaI=Y`}JNii1Jok5;gQN5P32mejA-uEy9KCRC)'
    'y!ad|0=M8jYTPbWC^CPzqz@6W*NL#<9-C)Ws^E<2W0*jY^@ck>wM$y}gKRHs4L07yE&)GS`AqUHeh-IIB%H7rGzX(!n%1UyoyBHkeM<}Gv|d|NX(B|MZ<oxhg~R*MDCz*ZpB9xwH>}Uh?`i+P4l4YgR6b6CT9us@gm#_?H2XB?OL6hT3P4k+-=llDnf3(~K&F|PC{ySX50~a;=zDn+*JuakIl@Xa>-UwQ?tICXnwo*yh{*HY'
    '&ZQ6JPe(XPjcLnw!Dy^*WxQfr5w%`iM=2cZ`Sh5*-hlZ>!UNLhhO;bvEuG9ud4*Dy$+K8v_Zc9T2+@Kd+Lj*$Z}9`|soo4>t=3PCZT!f^YYKSBCa<#lad!%oWl{s4&#v$ArKM;3=Jr?0#!@`0QWjHP{`?vZRPbYGq7KjPhDb)3uB1gGe8j6kTZSj&D*252Am7*k7`{fG7hN)!O==Qt;|f*J6jTg8Pq}4U{=r5!!YK5PAwmB+'
    'CqoIfb}gPs)wS<ptn|_{?Ra1*{&E{&+7wW#+c;xcjPzp$0pd46U~G6aU~gc0HF6?y(9p5cDIwO8Oy_!kBgJaOb<gr0hL>nqSre^H$Z+9NXFo~_{EIM)`+U-5TgetP=Ouul(xIPU=8;Cd?T}B}(nbi}Pfj-uJQ4kHmcJ&=A5N;IG>`zLO~%3=0i_V2-l3eut=I?_6{gg0OCF#!Bo=H06aLhCYU4`RxBB~b2-(ZQmPfO_Z8!h<'
    '*dPE}&*K$2<OwV29D4u@-`rw6xi<7Zff&8<*F_<}D+8d+-^Qt46DlWlG=o!~j^8F`Fs6xFOC7Ym|5BE7R<MD%)jBBqC(0y?=w9pGKdveec<_!e61XwZDx_W-9IZlJ71h_C&@}h+KtiO#kv$)`>p8&;72V$XQ<LieT8u8>%TzmQ8zd`30+rI7#aj42Bh>-LVNuS7^o*z8(pRkYeh8yZ&RS)f{!AB%Ky>YW@1bnoy$b=cG!}Z4'
    'j6>)E#ld!D=G5Ha>;Q~maIL({|7JcOuLTodb2jyE=ed)zwg&6nhzluQjR(HQyHh$?3hbRdq%_mq^CNpHu48LIP~fC!z)Qr4I0vmH{d2rO8-f3845d)JtgwWKR4O~aT=KT&o^v@rGwPl(7}nl(Fw*#d-h23C?SqA>zBDb5EOc2dNIp~r&yx{ZyF=V(>=CmAY`t`cO;IW_HlV&Vy^vPjk|hrYGVtf6o(<TEl&)MOAT9={1_X;&'
    'e8W0{u@6PFXOG{MmN>Oc2AtPN7(PPi<m>4^AmP6NV<5@xExdS2lmX<RcfE!VfNN@v%I^}~`nF+YAD=VwlhyBhx&FEKNS-UiNo97qHrkI+%2zOU{zKGC8Htlr%;~YL`xh&)SHC*Ukmo~$)*o5_nqy=k1{%mxgbH0*8zB5m!XUr&lx-v(w6KxZ;zWZ=IIpzs*#Ag54ONr2)%8|d+ib)R;hFNO^hPN@M?ac!@lV9@d`L+TkI_~e'
    'VEX(5s7G7S;2Hx)6JuLjrua~B=qoCWYl(@fT+qFq+3M71%?v;82^pn%td3z1tIwGaHCY<Uo~#l1#Ef6wWily(=6aZ7KA!G8@p1=%a3De|s4>7~9Y?F?w(E!1pHh1kKIS#tfm!J7dI=?CCHDo4u>qp`Y<-1K<+m^6mM=tRYuBtTsKh(3-!olwUs>QE(;aZi0&p*_KMi0O6YK1~mwM%KE2ZW3rxhjX7_)f`{Jn;Jt?)MNyao%;'
    '9k#oOKzx~FMxs5rbD6Xxz+7ygu;Ou8gyfn*^}=}1qCs+)5MZ_elWOnb1mSp>rB{;?_Z^2qll|dYY~%q0kz)Gvlq4TT9*<l_M2__Mk^Cv5{$Jbl^|Q^Y!Zh~C4t4SIow*y5LtJ)e13?vd-CHH#b9la*nZn+Mo2A5^Mcw4cT*h4-!u2G{uZvm9)r*rO?#+WKw3P>=qy|YcNs<+d-J`f5de@g{Z~2&~C8y(MZzP)o7N0O^3%WtP'
    'Z=EwVo}DBr(wcX!7I0tul}i`7U>N)ODB;5sOQ!JMSeY<KmG|7STKRND?~Ax*m{I;C2|3;0sL+K9R>oNx(+JR$!aBfGLeTCdaiY;0V3hwfpSr(hsvb9OfPdTf@U?l)Av>bUETTU`g{}*hU!o3~8j!HoYN<&9{<P_q^oj#X>pf@hx&ow<y?M@&vI{?ySzxjJ1Kxd*<K#&B>Q-_iy0wJhmUL9B_)hstA1A~P+-MR<`!?y3cB`nK'
    'PgHt=s+}5>cEbm<!_pM%L+fzf8&gw_;Z6iX2*4d4mud;__`FLIU-LEv&{o@|xe6dcBL=nSkhY-?5UHqbxyBh?#5Z|oVEu;hdio!^QK@8K0>j>28FKS((zXNe>-4=3GNhjeIgzR9+={-q&3=Ug{hHxM;IXXp-jaDr(K^P{JYiPm#?W|I%>;vj@zNAzDjzN>p~S!686^Mu@JtNA+foBmlS0H6%7i<5C@RMt-|sV}57VgSLB_w)'
    '>CC(QE}l);U}>!4V<(&5e1YXBOma1K^VKZc{ODAJ`NP20!5I9w>G`<lcKzRkNm$D|2ynlRjO9PIOjADLA0cO$d-=MWfCZl~VX$@S+9d>2J6Ox<)}%0)!RsQmY=tebiQ<AiPSsX-gH<<{IBbyfh*ScyY$nqW;?rY&dbmsf^ww&~L$H=Lgq>uvO5BB50px_?X)Wh(2FdMDMndq^v%DpY#bex0LMTi^wWE^a`?}`)?Y79ryly)l'
    '2{j4b;jSvJFB*^I$9p&!V7E4JCMlI@r*p<!NKa)_8vu@#1kAC5{6GU+{E(}iB%mj_h}u{^zto=K4o$;9CeGc&GwH=M^?<`GYdl)I{(;UWZ_jaWtY);SrYLml=Y<e1s-sb^_kTR&4a`>^)BMk!<4QC^8moIQ&JXa9uWpO2*U5bknJ$qF)P@KlPpl=*mJ84hrS7DH_Yhs_%4gdo7oK@jL1f<2q$C15MjbdYtrF*wkAe`QnzPBk'
    'qY1}xg@$uP-Ap-Mt(*ww=Z^%gXc3wwac^@?fpyt;ymLYzco;{1t|$g@R1M%1HvqnV1pa(cp6#=5IA%^GNhoLUiuCmGj%$I?9&7;~Q9GbYD~fg~Hzb@b(i8*N0y(2qEH5oAPp16CEqZ_Y#a4&4&mh2_Tx$HjZtk{Tsgj5FrRL^)sI#q)T8t|$VFKOZ2G8TsbV6_6qVuYxP}J;dp7|DZRdNfQh+Q%vecS2{fuu3p*`~>ocl#bF'
    'Wzf$}505A6onB-u{eN%Hijl-eASWvrC66dKgwLb^%!jLK8z)4_^~~C{bDdh+JO{&?@=({vvllX+*p+0VL}m19*_A)4Yx};HD*J=L+q*qWgePnAH<Jm;Mr~m_6|J84lk^*vS<Ya;QxX1gK_Ku3tDmp^G>f#VQ8yk&z1_nLkJ-1D-f$`)UGL)4vDLt7?**Z7(4=)jwjn!vgmE>o$almW&20d|ZjEqkNZwzdj*V;{!dB+@owW~I'
    'L<U%c#0R;;>!<_WC3Py3pT%p;F3QKyXIJoj&(6TR1D73#B??fFppGKMFkOI>CHneQyJXq)GAFbl2U)(uhzwwzV2t#=Fo`esY+zBg8v6NOKVOS~>e$c--Z@io^lAj!8HS`t;#DJrzx&S_zY8Ga#F?!TqmV3>Z=6m_<OGfcZ+_e|!bPPBfR$&cf;w2)P1Kuq+1u`<STV+qBAqTvf#1nX65Jo`^c<tQzCfX7^EZg>NS3z9z+&^0'
    'k><|<58beoSd(K+|E~)$Z|nCmI*B6g5(McL|6Q;x4V0XQ{~lgC`-9$ohUrO#3eL}?ow|oPWS1u`W5zr?`p{p}rKC+%^8-H#itWp_ucKKP%W>gMuoYcWe4pawydCfabiR!+tI}v{{WH;v_J3Hz{nwZrC=8h|vk6vCfq|UALkQbR2n(7LB7@u=UvDv8JUX-C5tA_@0GW-q+ID19!$>8tuOPq!lbEapxPU+X%&c*ug7k>Lsetq`'
    'xk3q20m6*%0M|PyX9|dy4Df&!$;STNxIFVlWxCEAzXt_?FPmG#Khy!Pyyy6vG+`-{e<1b^fcL$!*Ut5f'
)
PREFIX_RAW_SHA256 = '9122ec64e05b3f95ece7e05e1e942743330f7263562275ae3cf73c60940e4412'


PREFIX_DEPTH = 48
DEPTH = 96

def transfer():
    return Transfer(225, 225, (117, 270, 4), (126, 286, 4), (144, 356, 4), 580, 710, -14240)

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
    if energy['THRESHOLD'] != gain['T'] or gain['T'] != F(-14240, 40000):
        raise ValueError('threshold agreement')
    rows, initial_energy = moments['compute']()
    cost = energy['energy'](rows, initial_energy, moments, gain)
    initial = energy['initial'](rows, moments, gain)
    if F(cost['energy_upper_1e12'], 10**12) >= F(24600132703, 200000000000):
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
    if F(entropy['coefficient']) != F(408701889325846394605391, 80196432611780000000000) or entropy['checks'] != 1089:
        raise ValueError('joint score or dual accounting')
    return {
        'threshold': '-89/250', 'graph_cap': '24600132703/200000000000', 'coefficient': entropy['coefficient'],
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
