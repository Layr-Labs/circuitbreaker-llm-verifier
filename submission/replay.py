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

    A = F(5099, 20000)
    eta = F(162127403901, 500000000000)
    beta = F(62651104149, 1000000000000)
    theta = F(4618531893, 125000000000)
    lam = F(20628364843, 1000000000000)
    xi = F(17650958579, 1000000000000)
    kappa = F(8986158859, 500000000000)
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

    def exact_test(name, const, terms=()):
        bound = const
        for c, arg in terms:
            lo, hi = log2_bounds(arg)
            bound += c * (lo if c >= 0 else hi)
        if bound <= margin:
            raise ValueError(f'{name}: insufficient exact supporting margin {float(bound):.8f}')
        checked.append((name, bound))

    # Categories of (queried bit, product) are 11, a single zero corner,
    # and the other two zero corners. Entropies are in bits throughout.
    p_anchors = [(210650, 210650, 210650, 368050), (218427, 218427, 218427, 344719), (164690, 269201, 164690, 401419), (173535, 273587, 173535, 379343), (158912, 207258, 207258, 426572), (167077, 214289, 214289, 404345)]
    for b, o in product(range(3), range(2)):
        cu, cv = eta if b else F(1), eta if b == 2 else F(1)
        cz = 1 - eta * o
        du, dv = ((F(1, 2), F(1, 2)), (F(0), F(1)), (eta, eta))[b]
        k = du + dv
        flags = xi * (b - o)
        p = [F(n, 10**6) for n in p_anchors[2*b+o]]
        if sum(p) != 1 or min(p) <= 0:
            raise ValueError('invalid four-corner anchor')
        U = [p[0]+p[2], p[1]+p[3]]
        V = [p[0]+p[1], p[2]+p[3]]
        X = [p[0]+p[3], p[1]+p[2]]
        Z = [1-p[3], p[3]]
        for index in range(4):
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

    q_anchors = [(602778, 919150), (576629, 915010), (715541, 960048), (680947, 959632)]
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
    cap = F(123271600037,1000000000000)
    coefficient = 1 + 1/(2*cap) + delta
    return {'checks':len(checked), 'local_reserve':str(margin),
            'saving':str(saving), 'semantic_increment':str(delta),
            'coefficient':str(coefficient)}


PACKED_PREFIX_FIRST = (
    '{Wp48S^xk9=GL@E0stWa761SMbT8$j=Znj>f?WVb54>OiUE*Q<rM5K0Q{o0N8%)9L`d$ln)0%9l)NFM*XF`>#MI2S6Ih5)wy~UrJsf+uHpsI~?pJ'
    'RcxY#A9`%Xz0?ucE@^7y5a{7ITAks3m@_{HZs%P3e0A%{x8&_5ZW!7%|(g$@B>0#1XJa1mwdRK;GmEQ)f$9YGjMLi$VHfnjT&AUDWWtFSD}Bdlko'
    '>*;2ktpT2HO4Mg>KPP6<O6iCJLtTC*)#Mgv-PD<WC1M+5-9QDB8d7P~<T$<a4b>fp!ROniEiL%)<gSrhQF!9!j!)Tf7(E0CA?S#w>+BgEo4vy*Zl'
    '!UdaiWTZo|M4dv1R!lRtO0~1@PmfR+|5mrboM;j5oJL-YAau;zn;0$0ueeM{1WZz1#zQnYgM-$lX0miW232{K!U&;KQs8}073EHCmi+S?Oz^bxG6'
    'i#WE)1lZDI^4rfFVRm4;H+%0mJY1p8RrMzPdJ6}3Ro+ldwmLvaQ_U=!Co#p3kpv@<-*&N+(RJ(2&d`<-9Wg%clfKP|@F>Q*g)v1U0efHK3GFiz2H'
    '2v)OcXzI@AxqY=VXqo5W^)stLAf%P|>z!3KcF0N}2GE8AzV<<k|H*#7-=?lK`}4EPMg9kn22pD=e5}Vv-AGFDJw$nz{s?*0;>2#<j%XUD#2FKFu5'
    'r_igRy1R^dWI)45<96pmUl0lH+$)3r|u60EaA2ia2dC2_w?$63lI^torxA#^Pcy(83)?xiaMLe{rr~1lpxR`@g_Liz4$2bL*ea73`=#dt1z&k#F_'
    'xhWb|6@3D+r@D~;^fKg&3K+lZ1zypVLtowrEBb;(Akl<L&E1&mvKh_LuFG_f}vk=~7>l*VHCWUcaPZ+^+Biqo@nKKPji+K9NQdL8B`l>4voC}fci'
    'Z7AU8)WkvR-Z3VT5A@5Qo6D}8c7uFs#4)oWlIhY4gM`7-QA29qjI2+$E8;aW}E{@vfkj3d;X6CiXM;UpW?YyOGEZe@5Z0c>S85V94RQ*lZ!p|7+|'
    '`4BJosso^CywE=yClFvsT>ya${tCK2oU$zVqXQR@z)+Y+-_am{G(JpS-tFQhdl;HxF@`47Mt(>X5uui<elw6An$p4GW0KQj)QBov<Fti~lAcBwD('
    'y>|9CBMHggI{!0jG4nJ5ov>yj2y~Ya+%ngmy8~3X!nIR@<N`x!_W`&C<HfH+HPv<qJ-5B0!GV5GS9x#$tkY6cpkKjhH;hUMtU@$LC=@r#+u{GVM{'
    'g|h31V0<+rQC*r*8!Bw%!E{0(*&gQ4nW{=mhnbOnn#?5!yMmTl^+}`%d8pTPEJNw$IM!;S+uagUqj13Ni|*!sH}zg=ipGWgERXqpDdvVGae>QF05'
    'gY`F_g^9`o3`SU%VMGb>L1`eQwA%)jIii0bjGxI>BULY1a5!&VH%K~1yr<N8>^P5h}JBJbJ>2O^ST^xoWEP7DU)4utIV^zJ0H=534eYQyi0Te0AN'
    '+~q6X=8`Ef0Tx2a{&{|q<1}$GA;wre~c3{6HU@Bku_;ne{JcMT@XdPb$t-o7`FSFY_64&NSjE{_Yr*?LEEr8od!CFA+T@Mv=7CvN1L}=%-WJELkw'
    'h{PL5+OIV8oiG06LJH|k5{lll_$5kmA*mu;b4iNx;TNct|uMb!l+Ip5^vwYF1wUArpx=2Nq8%_v$MbOw3IUOV_4$LrO2bP?+CtiWpDtJ)>2w{7!j'
    'C@!A+5|NmR539?@FZ<n@eJMPR9oJm(nOp-Swd$jne|<|_McZhylm1x?@;{iha+|=-gZ~=apo|n5#ygS4zn~kG#Gu2ohhvbrx36(CATq1W%>CDO5+'
    's-HnUZxVP?ByzrIqrWuz5Ld-+ekGviT@=Ph0?O+1Q$ju`|se@-%)YCo76R)@!Noh(@V2QxMz_5(fIT-sW^2==2aOQk$JSU08aj{|@esbcj-}2M)V'
    '5v;~s{>nf5t?F(0F^OfR^Ix2;dW!cDGc~L#MyEJ+Nkr;VqnVvNLWxOly@ITRcJq7ljF{diRyHP@5%(F}#7*i9>&W3<@^?5n-F2PgNvW?CXYJI;G4'
    'b&XpEtRRqEqw@!0nrypJ(X51#z({RbSC)eGJqJ3O&ky&$<nytTk9f)S84vbGpb_~X11Fo8aEd<4xBd<J&}-#vdNJaKz7X^ST;;6u50!Dko|*i$(V'
    'h=WGV5G4xj`~=Mt1cD#-de{>~RU{5${bn;WRSUUNFpF^BbIy5b)RDup3wloe4?%Kb_2e_I)GP(O93`D3!HI><KMKmxrLnJe)z75j`x=VB<sjtMvy'
    '&z^~T#@OD&75wx>4Kh4!m#X~!S<fp=LB?v8^#n6hWr9r-Oo^;InaL{unZHo28jB@R!-NMJ*=$iX@VQ57+5MVB{z{=>OYKoF*s~_u!Eoyf-&mf|<J'
    '8*2rz>GgeD?K*=lXij+WyVry@Q~nt$Rf^cgZA&Z|>Cmw_c}`yzch~q!Ye^fHnp}ZrDO3bpx<S9wZ5Lozg!ATQ;n9t=<FrK=>1WsJ+OQWpf7W<k-&'
    '2_ieS*+rYZKt6vAjmww!22a&Sp3GhhoZ6G+Q%ZVn3RF91apT?_yWZ)CJ)HK_iuI89;$f^bik-)7W{3n=lt#Y0A+neo~41j|9rk;@mH`!u{Xb*W%6'
    'O>WfTeseQ-y<CqFPX{ywYG5xSE;&S*W-}&z#D0hIHubxf)bI_w!;vyuEw7ERuc%64rB|X!aS#)5wpOv`fq`SmGRO<@qspV?s2s0=jKR)NPof#Jr0'
    '#aQqQWSsQX4(eE~DRmTR@+X{5rSrl(%V7T}e+VW^6LR3{%M8F=cv_X;`S`_vo>b0C0Q7EY3Px5z0v&)KR7bGf-^ebV=<DumY%Tf_5jt{h!nv%1a;'
    'Jaj#nk6JFp&3AXaLb}dQY>z6gv$d_1OeSi=?&yS}%B!d>;v#rDvrLYtrA;B&rgpHkVlGF7fAHRlvttMP`nsnTX5Rif*Cg`zC@{`!N0MRNegtsjQd'
    '!y!<5@4ZE(-Ock0=sS&)zT9wh+@s_!wd)RXQ+&zV*54^d_(%OsgtLKI0s4-B?aBqDQG6A-xCXIQ4*8uE7l+TqI%TK$+0e7F9Y%6JyE!@Aq_Gf^z;'
    'C;R9|%Vkv~mjemijcmsGvJgEWoz>r)4=gRuzb0Hmd(nQzevN?cjl^Njk|0+vTW#(FHwI_Hq3CBENResZe<)s0mUP(I@-NJ~_oIDV`vrDM<9{Ojp7'
    'jit-+qj8O@fhHi6ajG<eX|;c8!a=QN#jBOJ;b|mxAwxwv>c2+s*1R5Vr?FzgI)H-&uHtXBzJ_>lih8)apdmMIn^l84P7fgy+Q-P-k;Ezwl3z%b2w'
    'ryx<7|vs}^d%K6N-f*<$Uok)eySdAZ)Bb<rxUkT0S@`UDPJ2Cv;rzDzOuy>xx7>vG5VzA<*7CH-Tp-Kj1R<nz1*DH@j>`%BIvKs!F@KppQ}UE7IZ'
    '?N?eYG^7JGNLkJ;&NW5=9W87sykDDBsC<sLsTAB_7yLe!B-Q)^q8|P9Ma_K*(9r=r|FvC-WYW<cHwp;;vixmMIe<oh-Y8r=mHPbY-r|2_DgEFtUi'
    'vyw;RLpZ0>N_{HTXr8C$8c#w%Hu)69mTKPxu-eSA-aFd8=ZoQxNad+S)kp5cJt^#PN1gUd&)oh;`D^Y}MP`24QLP+k^P@jtbk&3C7KS>@4Kh6T<8'
    '#!<bf4kAf$EjM@0srz4I6HUO@=NnX-3e*iaI0cN4p!4?BV>Id8c&N^2%>s38=)8b@NEX>heO{1biz!k4F1kEQ2cgln8P=D_enRoaRAlnCOy7#<5>'
    'M!Z*uxMZP56th$30E{u-YmnB>RTGMU!J`iQgYR~d4+3sXPHS@U6Tr45MCeU_9MA{CV}+uV+jL4@qJO;0q6NE%RlSps0i6tD^vq9HEEm&Ces?mJz7'
    'gK&$GT=W);2<_@7L!)=%cu!3B{8P9Qajq6jbfF}moDmIr_s;=IMK1FCtxrs#6Mv{fT*9qKGGn)L1eV`4+NnuvEDXta%#6K=)GjTF6%oX^(~Mrpq4'
    'M?{}BNx!M$ec&n#Q5=!bIYVxEXy6rmI_))}05=gHIsIvmoqk1%%!n0nugGQ@S<Jz^eS)$8%L;#oEdt8l&VZ@h4AYxRqG1`4=qsDM!+QuGwC)0<Jf'
    'k<~O8Yw2r~^*&R0yf<qED?tQHwSBPz7O-daxN*!z3vyWTi!E?L?H4C9k^?9H&A0hLu+;XG!}>nVHAkdjdnQYo39?#6A8OaDWp9q0UM?tmIe=Z7)Y'
    '0ZhN^_f*qMpg@2McZj-$g$SC-(8PAseR`zywrWvjYsp$gH+KgiL(digJxsm%}ut7X(<~%DXD^LRgp!q}<-{33sDu_faqgQK)b>CXGc71xezC#~r9'
    'J%0*#?=&k`1?|C5`{_x(nt=39M;5s9h^m9+>=*_FNb(f!Wq8fig^eB0eq$?e}=w68xBSP1vx0zX@Efe{=2j^_MGh20*pZvdoWjs{n~s!ZmN&TylZ'
    '7iKc=e#PUdO6{I{c`l}HpMcRS*##71$T?MYc?RgI3IeHDXsFzEnqpV_P&b<%VRZMWhtnQ3QFi*-a32jj+D`Tgd5cmaUUU*Rd%8Ur=O-$R%%0_pEJ'
    'yzKKWYjFdP?Pre@gDfhhS$DH`9FpDGuBZmzYJHn}43PkV-N>#7IS|<=!0&YL+gY^iuJ43YT0Dt?Hl2*q#=pmI2KU^i0y0ufIds6eb5V<`{0eCkqR'
    '{3gzm<rv4?nOaqX!Eo1SU|it_#p@rF$oH5ero8-efo?z!_qRJCkbKsbKcSMp;SWn3VTnKc-JtCgvIfya|;yLdIfYjl{QK-MY5i!|CU1k5d0BiRvQ'
    ')lCv78!!a^yLm3K~HJb3T@waXuepxlPsCT}N$P>><R;be{8V*g_Oy<}=fk6k*h`fBo6FfhGtX#<>p7aVgP1;?Gn3ZIxV!%H#`-Y3B_33TG3%h$jh'
    'lEKLl54vVE3wFcH1NP&TnM9Gr3FF?3sD@&OSO-iI?>r~0-LcL?Jm~)UY2$%b*xCvPY&<pq3-R#OrCq5WfaS=-WGIq1Dybnfby}+g7F-A!Yulx&$J'
    'wq3pO002G-Ex<EaCfwhcxE!xYW+O)J#0cR$^d@`G(pS*FA>1jO<eUm`k+FENv}p*BZP{Kkue6sB*4o;uov@6okKyG*i!19+~7gKoV|q-s+_BhfkP'
    'wsgLfcZJa8%sQrI*b!D6v&W|>Z`lROV(M*h7iKHnPtr-3RNeBqosa36W%bcf(J<9na2dPrEssHDtHC-5+*@h|ZH|cr&#x&nJY7@t@VTk_qpl(Sa6'
    '|)*>^{hJ_*D>45eiMens(%tsW)jd1lJGK4)cs3CZiN8lC*Fke8|4^gb*rK%338KeGrl3)w3T*9#CXAVw_;sKl1)etS@APQV|q@yXiQXQw`<@giU;'
    '2vJ^;^>U2OrYer}{2g!niAnhg0XrhO<Dh6;3&LZoQZ$c5+-vjT8<2e<>@!@6yc+R5Rk4=uRlMtB$KHzFraTfELL@*B1#}6&377bfxH<xy%fujiVq'
    '2M3((=I}mENr5uu1H7r^n6J<izgu^wZPOIU)rh!K8zAMxi?O@)rxL#uDkY$4F?b;Tt4!Dgr%8|Zf)|xy>~G!_G`=u#V;=@+lU{mf1ceoLi^XAVA='
    '@Z*~flICJX=N+uM=eIdP8|%^k0j_`fm(aSjG#LfE+{|5Jrk#}$UD1@EDi;)S|d3)dwzh^}Mf4UJ{Y8<!nGGBm)T%3gZK=zxBd_;flys9*HWechf?'
    'I$iEZilVP49CL_;F;NKlk1LBl5<)TYc7JCIna2b<M*1X(aZ3fF<Q_}%2!(KN)qJs5OCols#DglCgE}i9u!G6V9-A6miIz!$%s+gjDGGSYo+8E40g'
    'j0kkNS7|4t;)MCPWce>8i1p;WHC6uNrdjMoNBL>7YH1O*jrG6z#2a3pp@QE%vpYeCvK>kB6~_M{Xneu#!hmt`c%KDF#qlOP-aLu;Fm!6}IthyS(@'
    'm@LB3Zqylzki_5}D1~~RW0ybC+W_rDkru&xbua1!{;1G(llPNBUhTvB5g~?0=8hXi0NGE?Rdvjh9w|DBHlxDFJE^2vPo(Q<OlHFM9m}ZN$df{bqL'
    'KNVtR<5@?xnE600S4z9YU7VGJbVc#RrZ9qzL)ef{{YBZfIgJim5a>wP?{@`^=PY(7!q$z_N#aw+U_Zm5BYVHBgT-q8~u58qRzLsvz*B@&}f9(pY$'
    '@^9+VlsVx-WUBtVyg>Tznw%iqC>$sZ+;vIhOJ;F1z_ihj)9rWy3HC}EM{t^RrczYp{KP>(gbMFHp<fddY(ZdNSd+P>$8sNwsyazG7fNM>wroV;9X'
    'V$>D{Y&-FboN!w^)Fv_3z`{^ox>u@7W-;h0i+ieAir&ZN&nP*X9#ey<!4&Z#@c^k#!wWIYAB1~b`uGP>7DW;MIWa5QU&^lW6zrmsxKG8wf;Ng;%P'
    'WSEr#qWEvmT-?A5Ln)My+FVuqyZ3NGM^~*w>+O2pDAuB>pr^(HM~%M(6D5+9?!lDo<!WOTC+pkwTln>rg@Lu!lY{YS;9yD9>?wNN!0eELTHQjIp-'
    'zQ(%1Cu>#OH!^glFRSVQuaij8wyvPVDA`c<`#wwRiDQ1A+Zf{SSL_ZvN&qA{I7a6R2n%)yp$?lXwssLI8{%Fx%yr_p`>b)@Z?q*lE^CQeU{;o4Gv'
    'gAxKej7%-jLwYXAqZIh?MtJDJNscH&K<Q!QV=ow7OLs|U8_tfL-{R~)f)FU8|8>^+bAd@6X?#ELm8s*tX9~5RszpkydS=3_}E2<)xs?{JA>~gLa#'
    '=C0w4-L+jjI1ut=*=QmALITq$>gM$uxQObCKkW>L=3ALI74%9G6VVJM7Sd5RYsm;~rppk};=)`h62!P~=O-by4kA0A*m(ac;tsDd_yZ~)qAqD8vx'
    '`k2`)^KS%n;h)rX*WPcaVB8;nsce>COU~{LUFjxNgySFa89qKP&#=&uk(afTp%-+h;B87d%s3$Y<8b*=a0+Wd))C8^H~{Wv=C4|#M^8Wo{eTrRm6'
    ')g>XpY;gSHybfLyXj7a0W8ugM|{!TD?SNqFZgCZ5iEvG+7#<bWAN*W=}j&@8c}|vaqT7exYusMs~~GuTeK@1b5X^IafAuQ+#tZgDS?>`)#wjP8sL'
    '7$kJ3ZIV0M+mdq+wvM2y=OxBEIF1SXU{hRiBnIA8DWf3>~cuvh%tYFq^4T#b+0HD)WpF8Ccwa~YS1UojdgF=I*TzA41!s|>Z+i*E&G~?HmPVUGqg'
    '!3LVF<MM5UN~Y)?(wTI8frcxzIxdoQ-Fxj4+gE;z>B{wM#+l<J*l2gr8Jbx?QrwT6^;1n(+<G-3bS7|>gE+N5M+V#eo4V$oC>vttM`iGy3qyUkIq'
    'p;bFT_!won3Dg7O`wb_i3B2%q?Y3DB7ZnYvqnxhdB(F-L+%L^m5w8wkE!r8Kg72Cg4ctdfBhJO#oEu4w|bfAN}GnN6b-s}T~ypdgVn=92BBsliMr'
    'HOdfFn)O(dd6cdhJ+k$%;IQ>=Tc-wWi*Al2f6Q50#S491cV2Spb9H?0GTT?LeHM?pj@X6A8%gW${7*wes(}1Gs;13wO8R)=F~zUy`GX{fa5NYNSg'
    'wK$M+36o_}WX|eGmDhQ?}bQAlgIfV1ybQ{4g=9b(&B0;TB-afSP#dGq$=qUI(3*1Rhdfd_aOElm-A)&0^lxui)ipUH_4j9F;&7CxsYn23IwYGn$T'
    'oTgc@F&V}gL5UyqrqC>d_L_rKb$=DF=qHcw}X`WcQ++2d+8FMZHo`Q4Fzz^%ew`=bhocxkWUc$WdljQ9<<FS8gMb7e_E@Kzp7iBH+fsUtbk|`>u3'
    'fvgS2f#h#*Kxs;a^8N%(4$q$j;C@W)Ba1ug!|<5si-vIzKwE;(*v29V9ekK?$!yVWXp-}8j0Rnd@wCMrJX)8-`~UR`GZfVWYE8&8dfRZc|F?OjLn'
    '2XfThEfZLVD&!Bl<Rk)AX$F!Ncs*{`VYv@l1h_(k$s_2Mumye^C#A7^B(;?5JW97l(uX5eakr7^u^`|c{GuOuPC>Bx%Nvfuzfe6siX9P_!Yl5qj<'
    '@t|d>Vk%te7XsP`bKbLTRIVdw2HzzOy0TkIkANKGoEclWtMVg+=C@QwQUy8^NjF0d=nnEUO7B<;d=@;`0>rWB?oih*R+p0}*VmaJ0!YQ8a!29x%_'
    'c<EKTRN#;AihfvU>y*_{TMN?E2eKeTwDoo4e;-+w7#*vgi0*!Y{k{A$sexV)PK*9yX+ynAf+UpN}E@Jv#>`PO&-=Gmh@6Th&l?ozLEsx@2E1t7gb'
    'PgbL3=HAl{jOe4xDnyV=LcT*;QDDQ#Ibt-)nO_Fl~0G|Q=cP|EBdKi&<Tr#T{W87)WXOKT5U+MC~`hvouq(F>CwRuv4EEQ&lrG&VzEbntx396uV@'
    'V4K%d=W>`U>VZ7`LP;NbpNiGH@>@<x&RAckOv@lym@YZn2PXFG@}SLd7ca;WMPP5iOw7dHlSA3ZAD3}6_l9I=HMOKI`8;c<#vTbATiy&YD!%7FSc'
    'LQ7O=7un|Tk$p3lBH<ahN4Y#hQAf+*Pi9)Y~&b{U3#<@&-Qk>EG3dqdC|lJ8-u1&~SLGXN|g;1YH818!q7ey|{vC|W~$f8*<(F>a)38UsxE{17PV'
    'G}TDg=Tiixg$w|{=1ng?rRW>|)@axkBy-x45yA-M4J2mBTb5c+vVE$>+)fqvB^Nkg0^!(`(Pq^Z{rji*sH&LHFFln2Ngb8y!@m>?I%fiZeag^k<M'
    '`2|)FOKxbBkvZ7zeJZJX6Ly=Fh=uby8td>;S^OP*3MCo0=<ct_ATarqYbjaOue%bnS)N`WNQ4iIrtx-g2c3aMCD{aH<rPt&)WNj-zcpr6k=6K?e)'
    'MG0g64f!qvGDt;9N3JSMf%~zI0@=j^+X0o_i*unPrvWT(EkvK@{e<iNOgg*fAjFN~N2SjECn%U`@vQY4VnEKjlUL1MvBbjf%)*o1w%o?8XG_VBby'
    '01dd&DOMtrl=3ij+Hk-WT6Zl)?g_-3-)G@2Xk9ybOG3EgGau(RVm`t)G1bY|CZ8}io@c=HNPAz({PBr>oV+N3<Bbki6t^?NJ`Dg5$o}4#^PjoI%R'
    '}$xZombers<iIKriYw2n%1iNZmLMkPF|vg-I^kWP+x1j|@uE5~ln{@5Y#>uI7^T;6)?gkuVOtbD~^tSUJ&;?I}S#OO6TK||(YX!f;;wZq8Tg@=sr'
    'kV2s2C8g>-l{YFnAcYyMw#&dQ#u0cUNA$o2nMkA=HGC0$Dg4Y*OgVrc3nw&+lLA?-SpmdrSt;NLS^0u<{k7<>;#m8!qI01z1?0*asNkBoQGnO=72'
    '0qzu<h7Yj7V-&qg9AROKq7aVM#9E?}yGGXMV7HVt7liRR`%_#J&f6n+fTbUEx7*F^c-~BrthN0{{Jc-dqk#4V3zB`^7ToVPrU?^6ghaB<zrtEiTa'
    'lB%5a?RQO5%U}M3CarLgfE;2qGWbg^xDC?){-FMPr(F+u@_@Y`y!rFNEuG`-O=JpnLoz|gVmc1QkI_^N-Xcs?x-`xwnqES4;ErX1O71pYR^Y62}S'
    '$7>Sojfcn492$E8xa4L<dt*uB4e_wMNeve-UD-d$Y3J7V>uyMl<+DL?xzMpa&_2o{PPxXkw?w?$c7+qRXPXZpFC3)zyDZ(Uj<r?xV(gtR|$+W#J='
    '#iV(7-)GQ#%z=wbScN~LeAik6lafBK5daOmLjaPC*#jBFsh8G$IvkzOhcO>lHNR<bESRuEE{A6S94Wt8Kp(Cz*#KhLRxSn?RNpTO!|We;(&KWiA`'
    'DfLZC+BYHTBsfPmRTC|P4SPPbF>-+5t)kla>>T-I5G4Hheb_~hL?#@_)X}R%>6qTbfIE?h+NF@aF%tx*Ne@fTUVXxLHneN+TD(TgP0B`czLl;cGu'
    'CbizPDjTTdj2mRDx;OYSKXW1o5%E8U|5UNKv1q+8x>v*r6rkI`^u=y_ON(P|=U`k84w<5s>U8ud42Pu`umI`&6>l?t!VfnwF1UEH@#vj#|qnqfkk'
    'YKW6py=D87MNACU=&L*tkjT(FC!r;o6o#Vseaw;p%55>Azz63}}GK}^&l-p6K23{fTrln&xdFtFvhVrwR6h>3xfXEp-J<wUWO#J&6lLP<3YEnsf+'
    '*o^fs)xD35I_ZguRRF#=;?%Fm}5Jp-Vq%l8ksxma)?J449X^I4yykEHeS%l@U)p2FqNjq5G}48nf?>>KjEBvd%Yil#mNpXNiZN0E&PLt64_IIMnd'
    '4G?H?5R_jo2&R{)oQ5%CU=|2!N3L-y;`Xx-4;5Gu`k(jbR6s&KXX*3S<;sBfQTBKLu;B#tAKo-)}d8*|@ivErm79^2bsnW)($!Q(BdiAabs%6Eb%'
    ';iC>QxDrM>`h)rY<KLsukrxxw!zc9ArVE2zcj5VM6>NmFfdaKtr&H<(Rz#;_cS}{mSqSa}7ebQNts^s2sTMw94_O^kJP))|4{>;#ZO@b|G1WmA+8'
    'A8^LHxUWVY=a0y%RW0)KCM^$Kf}jLrs}<EpiU=aKsHH^b44VY4IfSUxVU(2n6W8Rq1)**e_7_?(WDWx@?UAF(Rp_I)oJxywlmo9qbB+`FG-sAtht'
    'qJ45TT1f)5wFFtA1wMI;J8?1d4MN80=vXu>f78Avb@}X$-lqbZ24odAs>ATuIwR6~fVSXKa41lh1SaB~mlWMNMYy`XE$~~!X<f)cyS=O)_#-c2e0'
    'Om5OTps(tl|{-3to|%S^*9zjjBJf(DWUPMBIRTgd6^+4YM4tZt{B+c-v~@TtCqW|6zXqJr@X-I=h1Yfz+_e<H%4!VPJ4W>@KrgCOu`;b4w85QfXr'
    'i7ks~~xHJAm~os6R5zwns33|dp0C3~tlHJUkb?eg7vo~|d)I*`X~bg;G&20&AXgTeKjQf81w6q=ZRtqrp><^gTy-J1EH)qf{Ud!8{8QR_3v$)(#<'
    '&<rK;jBxg$goVFfDm1gAn3f63e!Qz%G|z1(dsqW9qvBD1l!Z0Y|4>mhzl<g<HgIZ|d6Xg=$l}Sc#<-a)!PqmupBi1v0Nr&|@CJY$Toyd<v%aqI6B'
    '8C&OHZAFo(?mYVShcuuX;hy8$paQzbx<A_&P^?>Ef%voNaTa7x-8b{zObhspNw~M6w4fr`mO`W+4umQQG_H^9?)0?dHRb>n+u;STMP$>Hvj1>GJ1'
    'Ku=ODKVpDCtUzsIuObEe~WRD{4f{bw21R0zaSnQOgfy0GQrR*!uAsrlG?}Kc4&mm>-o|?8ZX;b4)xK9$1YR?DAo}ZXWCMa~zKDtXc&}YJ)36u`E$'
    'Gma0Vp%H%GoCCD@H5qxojM6(NA4Tgeaa_?1TIhtIYuD1oJ?KFovso+1H@Vz1=+*uKDnMc`{Noo1qe#z=XBx)uwQ??|AKB`_M9<?RfhX0bK1Go(r3'
    '@|D)A7__Udy{zq9s9cRMjVE>{{*AmFjvn5?=sV;43+;7X`EL3}5t3d7^wRQLR*vn@A?l?|Mk-=jwE8;GeEtZM?92};vsk18qeZ9Czv^}XnJF(a=t'
    'gvGLha%_CVTU0J=hvSkHLqCt7o0?s6stR6~vc?C^&sc|p5b64mko&o6+XRB?apWOCS3s7NK5I<^Q@O*sh$$5stW5k(HLDEq3?}~GPQ1K$!dpB~oa'
    'Aym=NA{4^D8tk?eCx`cA{MhKD(dTSSOBbn!Iwb#0*;}S=GNQ`$5^V`Bb-FiUP#}mOuaf`BpOzaUGz?DS3eR-Pbg_X^d0igs^#HMD}fb$c}^O#Q@X'
    's+D`F(RXVk@F5y9A8{N*_bu3kjSD#vjA)fs#SBz^WWa$Jb;&T@jC>{CS<|HT+&2N+yRf0n?+LcwXFs)b~3fi$NBrlAr_jI%)D9<2S<IU|BTe5pKm'
    'h#M#_LdDL92M~^@oJ*bp@L$~1>=AMNxr(iCf17l9O;a@AIYrfAm8Yp?&Zye`}9%Kz+)us!I#mc!V|9e)qr8u9+wtl>VdDrR2ugIw8D(`V&jAEn5~'
    'HL8@Jl+DEHMgnYq8|tv%Upb-K!G4k)fuWy<(RhJx@oq_P<p%-AE7p<HVCo$@pCtcRIvAA;?pv4hAve3qcwv`kPe*(JkDvs{B39l)MCu$g2$h`(B+'
    '3ul_+GMh8u2w&flf(LL`WoH{CbYR;{v9V!lX~AS@@(JDZJMx|X;Klb#3p`n=%x31TI1U0nh_>LyU2btZ{APS-0{<w2*+47OYoZnXpEo<}{*x}Hf`'
    '*|hWq!2r17F9J=~-vB(XW!S-_<NDz>N>}sJY8eUk|;s<;5FHmuMcubr~duivz-{Vr7&52)+>eQ1EK4R0BPji+!J|zk<|^nj8IG#8?qHy%MpoV_bv'
    '{gAfm86(8Jkwm*01TFCwMGz}+moRonaA0DSr&?xnkIQ<UUs7p5Lj>wr#d4aFw7AIs#w8SeJzU%CZtB51$fM5ys3w`1)yitT`n3d3FB2wGxs?fg~C'
    'O0}&7Jx(upLLRKsHNFT&VS6wy7oc7^#LfWr01#HWipH*Jo%2oC{f$h?<bC+P^F8ULp0UIYPBzPYlJS*2`RGzz9;1q^QRbL)a`~pe9b1-i4e*iWnH'
    'Z8NlUzL!bVR4<m}$<79j@rynWsa4ly=la@l)_e2nf;2cjO9B^^xG11n>aLxoTA2B|=~Kp>@fTn=Xc9avrH9I8g{4mw7i`*k%ZHA4~oe3HiIrUD-N'
    'SL*3(T5i)(9L*k$7Q?vpyred)1V){G6Kx|=cKM6)u6Q!xbTYkT&r~ca*S-uJBzxPks{5WYM_135!x8yG$n<=`{auSV^ffPV{TsM)%{PkVee%+;d8'
    'xPkV-q{k@}B}%F@wq-z~3Fun2<B-GEP7dLDa|jc6Ma|58;n;?c9h)J<n~7&Kl1e9q+t3ZF@^emc5R_8BAmXQCh_M%mtbt)#14FSl<*uyOtzgE~H?'
    '8+pBfVCn%mKT=Nl<E=<MAWV>3_Uf*!Q6QEZm;tgw_O~i3Q?rXBtqM)2YOvG+a6}gs3+^T(nq>uG~HLS&k-~bb(xoR-WfjZ}m9+psEwF@+|EH#Y5s'
    '2rKrHSVG>Y7!3*`{7j1l5m$2EKE(RA#`Cf$V@xRnFFRHWuxjBv)9F$o_DNtJC-i*^w%0Gjs1Zw%7eq{nAerG2V1@0Fj=QA<16;%5NTBs8vh8>!sh'
    '{oot;A;<Ijx)j4(RGW#Ed76<G@HJnFn$XeV9r6qM$9{SK+%t!}*e7uSQ#*Z49zjN?&GB6{r3WHmtJUnL140Wg<(#AdfO63Eu-`wQ#180FXar;q{E'
    'ng=!BL2usyLf|0;MTYpZYV-lIRXKV%?Az*sQmEd*B8_t#!-ZKeFQG+Nx>c6YR9YQ7``C*SE&lo*tyz}NXiZ9!!f6FNrQZT&hgRlJT<0y!qw;tSH0'
    'H2nA@q54b3PhUzVZ>YFNt@BlvyC)@pM7WqnOH#pZ5$~n`I~&cJsmR&!5Mz2F=2e0`pJtZn06^+|&0{esSc#ZfRap2Cn0Ss{Gt><+ZBwzm&pIb*8P'
    '|2pIr)p0om?iMuys_K;4>Si%K#tnLvLYkIsQiG<;LTFBw_;c2EHG&$<-Uv&4@%bKfd?`u;VBc$+emS58<XSF~wD5+2SxcMZ&hwbX3smt#->^)Bg@'
    'xuHC*eujcUBi8kT{`iem(KPtFVX0VXSU}?JV{VEB4)B2%<Io9!3=YDVwBR8`u7Zfdpv5hmR<O@L{ZPu+scAg7uNeN9Oby`k3E}-E7V7~ZO%OvUWn'
    '@%prY?k=5pAkcvM0N0zhG-DLNx03l}HjX3D^Qun8?oY}2?21eQoJeiDcysUV6bp&weU#$KGIRr+q{x?#5MEdKD3^?A|dI(=zQk+-4cc^Pc=h?bI;'
    '#}xlcaW&Tw&Y3$uqSJBV7oH4Z`EKq*2LwX>BBs&r&A$CR$rJLgmaX>hs^z_Yvb{Jc{MB#mS7kAR4ek->(+d}xrjE_7lAB&B&iFUN>^G@*DP=A&V{'
    'mMI!&*8c>$)CM{6ghutcv|fP30OKmeY5T(ddVtLFAH1eEtylR9z%=_xYqY`7Wx30?3swZr^$4p0Z?dG6FHAbAiU_P)AJ+OMR6Ywi3IDIvF;Ue0c#'
    ';ekv^*Qz7jUN|^c}8XosmN};tr@YZL?D@`>NZ&~=JP4`CSquNV&J%%qHx?<8Iwz@I_W*9l&WNv*`A4BurWkWGhm`~mSWY>kmaB2`N-|BDGT?CbX{'
    'mFFLX?{FTe_uaqt;$O3CMz>QDkAsjpVqOOGahU<=sdfCFyZdRqD~ISFR@UfnzH>Vk)RoKRI|?)JR43W=cdo3G8!I>8Fu4p*?I8a7P^H(a$=X3Z|~'
    'R?907Tvc8W?4wFB}HcnnICUG$}ZmzI_O*)HX&1{cw@T_loiO?eC5Ca8pvIEYllq1O4O#1g6MA_#l=g?9$Cg3C0KR;2en_j_8Jc>%(LttN=!*NW;&'
    'c$S72Wxc((m$Y1GRYLx!z4J&9OSIjX+SHxC<KD%P3*E}voo;Tstb=WEZup`2QH?o#h<7A#Da32&!{9*D;ILRhdK2PTB2`weCFys;qw&PeF|qGN6z'
    'tDP*#UX48;e3q<!sB+DSqAMwVuVJKsrNTB=h3;-K@feB4~eu1+hVjZnIT%1u^ZZ;LA%NH$TjN>mBdCiAFXbZ3=-HF5zg6=A8UNtA)k<-y`;0ypT('
    's6fyUQh{j)}i|eo?i2aIKg`}=S@`i2v(|C7+&q8t`XwzJ)z85v<!j`?jVBo>}h`Bc&FO*HY(m>jfiKh27fHKSs%0w>|5R^lTKRp~sjn3G81in|FO'
    'mzCcax8YqMJcOF^IXO0;ih}VER<W_H22!~NrK(E;x<i%d$0+q&)r*RyuS^JHKD;*LElZIMGZP|9rtB7r#h#s8S3Qc**vW|omPRB2iJLb_2~^!inv'
    'USBX-fNIB;9tZUKr*0=V*oqV>oi-*e1|>z-QeN`4>rnjP9@50%B^301o(b{tSedOLelC*G!Of505M{Rxa8x{3otT10gKP1Lp$4jXpjHl>ths)S$j'
    'fLkMYhN+oMQP}1nSQsA`7VF5ygDD27Ix(Cf(xP!9RJaB3Z<&3XIDQzL!-FrPr@W!TAZ#QatUr%OQxq_>iRO4RHV^#CGrzpN=Sp5d*OYb!N%<z1yL'
    'KK;w!n2#n800m3#0XmCR^75j;dq?R75z)+h5Vg52I?mL62GDMgVn9!?rW!?@N(RWID4d8qciD5UpUl_oGb?Dr~*H#|lFi-w^XP{Zx{X^PJlY5^_H'
    'E^6nS{|8T95{9Ru*d5LYF&<uCu(@MTZZ>n*W!F|;M_lfu3azQFi*6;Cct|u<WiV`PW=ATAQe_%#NPD&!X5OZ*0&4-X=KSLjNIUb3WV}<WvU@<;aj'
    'z-_RU59~=lEZw~UkRxqICW1JiN+UF1g<VA6~DPa+dhSck$Ox|W^U@L65zg>re17H&c)$vd#H@lwo9uN6mlHC?3Y4md<Ly}>1&I}OQ2N6gfNiOeTw'
    ')!MSQ%W1?A-=1p_bzfcLp*PyQ)c^(wj8m}m2ds;uX(r&~=rj77(AKUxRMdB1qCHjxkdta!I$0DKy#o)?4s=17`cd1ksl0sSJBwMuaui{|v6-=8IH'
    'm<T#j^SlStVNg%rH6BVBaf^uX_l3r2ay@L<gIpTw0AYBPdq)S@NaaA2-ymU%KWJh|q{PE%U}9+>$$rETTCjvBrP*X8{ecVOw0&q~FG)$*D3Q%{^P'
    'PweW$cv(jqfg(A>xK+7Z~?T70E7`E|7cV4Dj%lWC!CWnfvPgeCxBJ5gkyZw0PD*R0-IyZYUgy(-E6Rkk2JrmKEUi0w0W~{OvK}2)C?H_?Dyj8?}~'
    '1QgjV)qa>ftSHLth&?EH)eA6y%^lu@G%8W0<n~}Q%R*N@{5QHs7H&i-zalz0h6`mLeYb($?i89KBGYvbxh=!`XhRi{}B4#h%sLy;Pes+XwW6P$Kl'
    '%l{iu8NAuzU)i;Oils!wu7bUi97P#eW3wnv=6^@bGwf7T?ZUmx)KiRQZvBZSatn4{(usgRUPY@!!h@ivRn+YGuIOTOwI0O@m`FMr8y!^G#D9&Xqi'
    'qVxeTIDE*@?%!Fsjw)F?`<$SE4@g9g&dX6(wPfH-&y6C%k`x~9@+Bm}}&3<=_EU@8ZB^Y0ghh`ym(PIG@`!Xj%&HG>lA?N?zY%`#a>?m2pPd9@rc'
    'FFZ113T;8DH11B3QxpYpc)Ua78Al;2K?&Zc4KlB|tY|LtYyrCg5r#S`ZNAnkz10IhQLfE=pfOG=_z+cJKp_xtN=Jxk4#%p<d$iUfOXadMwnoaiXd'
    'Vams3hBS2RT1Mxrl0s(2o8)7U0r-rV`;e)D!~=`92x43M;6DLQ??SWG?a%bd)jBPe2OIyRCmpyuwmt8deRTl%v*g&4*OnA&Ci2kkVgbrQ-)4R8EO'
    'Z!!k}FRg+*lM57i=Y_R>!MmlZ*!%EQbk{+&BcXmKe8UI46Z+!zxPg5Ml5a9CW5FqfVBpW$TBFY$lqjb?=g{p5PzjGkS`iN5qNv-Z4^^hq7?}q`4R'
    'u#@=#i^(r73t6w&!Ks5SH8k=_JUd)ynlSWi~k@`c^Bsh4p~A2l4$DtO<x$zX&FNoRpI4|7o%!NblxQPJcIr=4lgBO4v-&9r1+eIh$J9%{abO-KHn'
    '!RVtVI%tRuI|Ma<l+-mO8&43RJepWFHMA5Ukt3)3T;m|)7Sdx5&Qo)KOcm2}z2DoiG#Cf(aGZ^S6W-wQDyZIay>r>hN!i4%;JM*cP=YeF!o@$5iU'
    'Xe80D$t2Z*1MTvEZ(fK5wzI)Lile_dz4#pC<Zw=fYDAb^UF(K_f|PkhU1jYCxzJ7RldbUX_&_=MuScI(L++wN|C7h4_TGXb|KygCxVnfQrq99WGz'
    'E4vw{y0)za&2$m;u)jE`Wgf2-ZeGkYh#zaK|w)xmptUY+)Hhyd*B$E##+@jZn(X;7dYlu;`mO?M<jxAl1BaC}qmiv#$>`h1gwdl*9#Jg`qG~7WYk'
    'EXiJ?khee%kEhy@iiC7o)UonRe8V6ktaWHpJVX1&vv|5rl<x2N3V5+3w<qlo(R_sWLmu`tB=oR^h+bl^QB<>@}UJA`QP!5}Pl|x6MRNBi;114J&*'
    'lO8sWg#6ftS*;3FD@dhE){@S+MOm`wl*8G9Fk6yzL&A9K7cvJ>N|B~a;IYqg7N;AqX9%g+lAo_yc&YP3LWi7*Jq3GjkAw?GIli!kYYzkc)&oWm#-'
    'b)yE5r=(n2kC`dt;=!IpDYArM70-+&fI42Y`J%+ol}$tQ?UL3t418h#RyVV1KTxg;<-A_7t_VmM8m&NsN#`tkfKA$bY3$6v3gJDA?kr##RlUp)RV'
    'LJEvk&jeU7Ld4kWUkt-E9$tg5eXWo!hIGv&5vyU)L)v+{26p15?Y~XED2gx^(IPf*l<d9q_$M^)Vxo8tb3j9r58fM@-Nn@o3epgwlhK8Q@Zw)3^x'
    '(<H>ncM2P{nOPbJGFIGX5@0?(-4<pF^j(wVxl&2k~;DLKhd*C000`dJQH)lp8Z|Y7MI9(HG%=lAI}o)cnC(Aj{R#rRArD5KTm)ie48|lQ6V{@>z`'
    'dJcVV-7xP^lUc73H;B(C`h8{N~{%Qt-!+xWT@z7k8-qBbIWP&+K|KrCH&J<GEirRJP<`Q*bz*T^Y4#4>&h!CKhB=TZff1(E}`7S~mN=*^Z3icX_$'
    '=Ao1`(I#BbWyJ*Rk?7lxam%Be~vj~yYR#eU(rh+N;bblcAQR=h&9)^jr|jfGk7O1Kq=u11_urTZshHFXI{r+v;E{qgk4h&p?z!ruqYyg|6;5M5-B'
    'wHpJ2DOk$#hzMVsisw$j;?fcPyZQ>b3$XaLg)niNkMUc#Fu^n^W0qReCkPIUZRD`|_7CcD&MuRKL}xBie4M>YcFv;s3?Wfms;0-!4L?>=C-f-5u|'
    '9!S*GrgpVu0p|wV)Q&FYp?1aS4XZIPOZ`Tog3g8DRR1PT<W8WXyhT@KE@#9?g}9v!;Cc2#=;=uA6pRE~lJ42C$!hd?Lte!WD9776enidpeowtQ3t'
    'PZq=EAmplkZv#v2a8OXi^~#R{~`gghzso)|}{18LI47HxbSxK$vN)i_gbDh%Jv<WOP7->(Ag~0QdCsB##{JH&a<FqBosywcFiz@cMZMInU=!*<2='
    '~4G7dP??0k|-@Q%Lo7OY2VFe8;1>0|?Fo(VcH8<~fvm0b@{?0mv>4?8pR|Ta6)MD}GS{%HYk0%J%G1B8pTI>umjbtE?z7R}(f+sk>kYq;)b5DYoI'
    '5L-RA|6^;CB|ArT9#J!>6cee+;i_CGTLTnKuS<8lNB^9FMo;L0!R=Wt8{ki>m;8yZPdvMw&~3sm(@jT5J?3Gxso36++Z(1ekipbi;86MvY%Z*II8'
    'wEHShPP-GkZeKX$C^y9#cB{{!Q!nBr#lfwtW`&||HnPKeMwfm3}PC1X)tJ<&cFATSQFKDStmxjP3S3u>0=7WkZj$UU3p>!}?RxMd`~pI|Wa4ceMm'
    '+CvbB(?uJJI{4fJ>KfY^VhO3l+Xuetv=}C3{>!YU`UrnQ?aAK@F3fvoQ|Q7QZ{yuf(I})1?V{~HOf0NYh`q}3>5=aI&i~|fw;C5VdPX>CT?MwBe^'
    'x|gbt43m?jP5}FnWHBeZ(~y31*tW@~6d>r<ydxpfMg@mP6&gBA3w?>ptAteC`o}l4`33B=I}5<E)7YV>PO4{AwB^!5kX<%Lz|c%m%7n-97W=oTa4'
    'dh3IPb$)vx_^t<sj+i0lYC}k?rCnRSvH`dch9Q*Sj1D{fB0B!P}u3%#fQ*lA23n5_$em3D{DPFrmZyDGH5!8g(EKXqaw9YD-Qs6g({F44%({5(iq'
    'kl}yHYvH<*RMMJW4x#7vqR<D53{(Zl9nK9SnY&R1D~{CW)64k)ypkqbM{0*ul$zv#LpyZo>c-6P5dgyM_xw6O$|~1hLN>E3G|$F7_%OVLr(>AHo5'
    '7!pO)4KFr!V0M}HqcT8XXY+tG5D{!A#ux<OG?vprx<?KW9S?y1t(5CEL=PMF=igjO!7&<qX)8sl6e#<k*b)pk;Wp!|Ys9j)*OqycPA$bT^|yT06m'
    'r4CwuBeHtue6e2ih621|d$+{UUf`4Gs>UYN_<4$Tfv#AI7!^$hW~)@P6*MI+I@MGZ#didbD_mLAs_XFz2%!B%is1&oz?D<ENYTAEgO7kbIr8Q=I@'
    '>gh2X79!W|NuUk1Slghi|Ug1fD?}hjK@cV7GQ_wJyz`RO2lop?aC9HSEuu3<M!CH^Nv$e=bRzvTir5SuBctlu~RvV~ld^WQ!v19jAHDX}^%B9~j+'
    '?pqKg5eVJ}hi2nR)*6E+~lG#8inZ})~^RF%WswxVT*Ibc{`7%KUkZaQ?Tf$%VWM#U4<FnJG2UbZ@@*NXM>L0pmgt<E_Ww3C~^Rlk>KqDP(uSeu?_'
    '?Hds(2p$+LTiu>Mle~2c?>udkS5};3MZHBdz&Vp8WjCq`GD3=U~f0LUc^>3)N!h+Dka|y0?9dm?7!S>ZX<EUg|5g!RV&p64E&cu`%TCD#HQRIhh%'
    'yfMaK5j>l)|d1>7HLc|53f7EbK>NUxu9)xH<VL>m)g8KUdac$ztr`Y&&&MspJF8<zpm7+kkSWWbA5=h@+vbPM}>R~Y(n_qMmvZ{j|qUhcE<Pm2#1'
    '38CwA6>8P9p=ozhfCrjY6xg)1&OFcG#S&S04Mi><pBxTc-5(}xiqD}>i}M{cfachlF^NHZ2YFLdnc!diZ;nae!^Q4|9l7QtbCnTlp18d<C(>3Vk+'
    '|y|Y*XxX_OUr-B<YI)cXpD}Rx6Cs1x_c#aGFLz%2W&b_0&1U#VfSX>fM#b{zBL*q_i4h{%6X&Ne3CYdLziMFF#tbYP~FqRO4Vs;muMQZe=sgY<F`'
    'I;V6uFL+*@b2B9f`NE2$2h$FSZjb&IvF5^KG#phmts_w@oxOMr(S(VD^4BZjt?F$9JB6fV_jG3vb#?v;a%X^~YpqCF-pQnZtF!X8Z*}#xrl-s{CG'
    'c$iAgbL?gA;M}Z=Vp`w0U+8E-a!mP&y}imdtzBv(zqdVsPvWmGEZD2Xb^@1%w`^fJK>1LT92$Ud63ACsVY8GFst~{?5Q!`HK&jtycCD?RUVMMo<t'
    'claG<R-4s>Rh5jGs&`)0Sluye;%`{f|v9+XyX++8C+#l5ar1)QZFdE(Xk3xcj-HxL{#C)SOXz}Pz}8WOrd&l69y9J`{MwRbWIbOW?#sSKnO=C}s{'
    '_A<gSVfa%{#~UgFn-Y}e4Rog=@LMHp<tn8kr#SRo8Q8!X_F))ksv{ZquDNNLt?GyEIk&Juhc%VX#(@>n&YsQd)1<U^8l)p#V&BjYefH+LPJBuDY?'
    'R>jzG}Ty5Vr$NBGg!!rCpc+aDhm_KzP5-!_imtz2OqJ{Kd!3YhY^6mBBQct_3FWQQ^~T=T^D(zx-mrvRzyR60oAO=+e)>UR&b*>gKrPA%%C0&iCn'
    'ADRIC~w7Lcl>#D3UN=6LMRsL0tuAwASWF9dHR3-a_V67_uCajMwCb1{&a0^|`aTM;FL|rO@2}U#gLz3N>ZK910TPDe_w0{m*M5VP{Dc&SZ#$Y@qA'
    'i7o{cbi>_^%Rf5OKr~Zc>>VNMHK>9F2Jm+zKkvokoxm1V#W?_0P>tRtY3xiO4KpkxXmsVq&P5@fZ9SIP<G|tP`SbrG1e6}qX^SN4>~Unik}!nKmL'
    '7<yt(IJyBGaLFGV*=*~A6s9mq{*wwo0BCz}~LV8|FH)xSB!;ytb{Q7{m`1J&^Fxa!br3k3?}G>Ny5(fYcRr2c10mt{m%6>M9ILmSIY-08W$Cm1Qh'
    'tYuk{SA%QW^u!>yGexdGA~W-S;GDgCK^*GhQJp4EGWu1a9$H$tH&zIwvM7mT<4b?Y`IS^9&ClD;{NWkJ&DAGb!vLH_WeQrqQ#k@gJo<R<f&Z2P<;'
    'Y|OlJq16x;lXNdFswsde=%4;iT|T-qD}sCLWnb!1e1$<Ax2|#=jcZH0DuhZ#b_Pu}I2%`e1}vZ6&uNe&}5hogA(0az%uEpgrpKj7u{#uT03?qM+{'
    'Hhqt9Rn#R`B;P`?CyVa{TEMkxeGg;mS!BGEzPet&IfA)<6)e)JA%erFyKHhW<zQ(rzdd;NI-oeM()=&qVDY1s(zYUS$L-TxnMwQw{n?Kd{FDVCTy'
    'lKNOwp}I|hqu|L+vG*S)Fm^%t%v>uO1zB@4|lva`T1j<j1>~Ku<!Gl3*;i3F0|&S(L_T-|M#<eoYE@39MHzLM~TDN#Yf^A!f$i0qg!<B{LaAm!oz'
    'J-%83W6B3U2eS=4&<fy(|DP{*<m4yN{5uQr{BqBFUz<w;xZefX@TiJ%7}<rJu8v9l)N7LQ%jxl<{(9x|KtR+|uq&fg-zb7up(Ype^#G%tRWM7Rgl'
    '+I?{r{x(tvaY~a<z@c-awt+1A8zAg}YMkm4Qnv~TWf0ejlc?xwumdqLoJP)FflIE~My|o~{G8P>KL+d>i6$6bL|^pfxffSotk#;jTtYWUW<QTq#L'
    'Pz-A-xJxE;}hB-!+dKJXj=T3&-?{-@qTEfVv!{P8R{GF4ZgdMI03+u#*VdKXo2Q^X9T2e!~tW91Rl$h^Am)u><jS$0+E6*@ca&s3qXBdkdc4c%Z9'
    'a+z|FaG;&h*#|Rx!yqRPTzE|o@dLWx^3dtqT+*0v&YBDQ6oSR?|NbB_8hw&dl`oLjB8>MQC^!95rhesv@N2GH5wF@r8DMc@}`5ipg4<=1|O&2O&I'
    'I8bIMK&lPdvd{#y7qrbRm1y#!ijVRiEjK+>+abA1^HjPkEterKt?v!EiD&;qa0W$AJLTT5P@^Ayp$d1-QA4oC^sNsNsr*jGyeMgH!pQ7xr-Zwax9'
    'p&-_2I=!syUDxJ2o$D9Eo4*?8{Awp<XIh&6C-XTTEruS8s6lo6(+wfFBdX#ITE%?Ghdy`?Fd(O|KfnB+g+7vPQ=m`F52Sl?n?oOXWvcU!^d+c=ED'
    'Tsw>s8s(n!VXMMGL&m0qjLQNSx0%=5=9PlK!sJ}YuWB0QO(aqt5mq{?M&Q~;%Ze|N5&>MfM+g{1?xSj(O6)xgeXQfu#CgXGQTdOl>1oC}O{z8yC8'
    '~y&?KXE}%!VWE`A{k4e<xwnokom7727bMJ-tV@a(`hoT{FRQgxN#=`(TC44#{2=d5DxEvr*Q8^+Dj!xJ~TbZ*MUY2&R(wzxZ13pl+{NTStxG`D-Y'
    '%Bm<&ysgzAC3#+dTvCS1t{h6w#I#qhrkE0r`3KIYcxQZGke3~6RXCGD<p(HdWeC6(j{vrdF{V)PdZ<i^Oe1B~ly-ds<lZ}dnUaTF-v(-Wew|#nSA'
    'BpFMlHX+M5~SB1Sz#rA`ZMX$TlR6C24Bia-?94Q($S6Akc|*Q%sO$i>&&)Xyw?_veNKnqY&V!!={!qEHK1LCzXC3l$+H7z6Omt!yOT-(Gb`SqAan'
    'B(m9Sd;I|g??FM=QMJ-#}ZKL)bLHV83HF2kk!^3VP9oS4&6E99s7bJ%sOMy1lYb0u&$fgD8wC+Xu+-3Ji=vbacF2!mUs4{)TWui}&725|QYAZZkG'
    ';@yo#v>(u{zs>6=;imwV-aLvaVYoOBYwIF>{*Z5;MN>TdB?`%aH<Gg=yCeo@D9M~(tu~FcOzarq2GF46HiQl;=?EJe)a!;wH2ovGj_a13IZ>8e)7'
    'T-*I)goYSLwGA(Oh;x|BMLZXa=&RNg-O`mWekMgqcuN?AWi$Jcy#Z&`Vn*#N6bPs23&#Q;g2{2%N6T8ygnRp7#Di3N2!t<ePfrX>(rpav{(0Y2_L'
    'DcudgB1Sf-b{5FyP0MNx2`6QT@Lr4J=1mcfeGB?s{jc4PQ>^Vk>s8aN`CzPAVAILYDg<t(=D%{-*!Y(uo%nwbclrG;QeQ0xe#BP>s5M2ABKTy(jR'
    'L$HyxN=s1wWxg-p9Pk)u6kxi#xfQ@wWk?#I@(;2gliRWs+$TrA@RFNodhD3gs9|O6e&HuWfa`@R+cJZ)u{Hu7FN!E#d6h#(Nbrh8)G!vA`%(9?{3'
    '=q0Kz`2KxeUPrMGo>TLe4hd)T8NAVPJ#>qQv~qV6H}>bK{--guIe3O0Zie%j%KQ&+r_O*r<WYUk~ue{~O&P!J0#^D1^D+``e!SU5bFq*5{*yPKVY'
    'IBKnfw~?@{QQ~u&Go+6c<b^PQIS%n3Ddh*x8vHF8&7W#k4l>WkoAf15I!`TQ!Rg>liQ8fb1}u&^Igl@V`jW|@mHGz>H4qCm-%6nf+?IFpzZtdz1}'
    '<U4kE;%wWp!4rPUId@*P!$1JO$QVK*I|R9P0s}oU*%$kT#t(15fw$_(?so1vz=Z(1Qy^y1N@CDu-oT9Q-sSYE&yql$|yhXVnN~<jpKiP|XirINVY'
    'm0@_0}<OZok<R)Z3nL_cw5HPLCrO&GcRr^L`Na{!#Iqf@d%J<t|)(lL<j+Y%O9-(9x#^Ub1ixvE=eFO|nJG8@8ZN4;Ud+=fhT)v|eg!^(bdb_T#>'
    ');>ofQW8A08grm1ZnUO;A4of_a=P1V+O=*naEUH%|?qYK>qMcv1c05_DHZ2pgJl4VP?ENh7CkGKRsk-)N!3i+(xQ6hRy`Pd>B15YLt-Pod^YhEC='
    '+!<~^bMgL>6evwdwcbZL<GB#a-b03qE+vvzNiOy3xfJmb4RA30A2*t_jMV`57LB@MZ_JwXD84fT{IJ!hudA8{ZM=Sfc0^ccAzK&r#H=7ryn8t}z)'
    'DoiCVc2foRDq!9)PdF}!c}F_o*Vo)xDb$CwOKggxCMpkpmj7gGwQ^{K+fDU*WTMvqR(@iXxC>W7*)u`uS5?m=v1w{c<%lmTg&yTm=JO{5l{{P=0i'
    'JQKL|}u9NsvnvFuvG(+7RlvM#5Iqp?QXMM9QxT4{^DUQtYfnyAyN%2px~bSYR>CN~uZVH)90@kWKm+G^_VRgXyFXK27me@oBANJGs%9V!>(XofFE'
    'Sgp?h|V76u<<R)0vdg<=ac-%vk?v?s$cro60F#US=LkRtzod^$eo+npumLkFBj%Vqm3oLR+di=h(zTl5G3zC}8!7|UoM$rAu!*d&hX;IU~^F$VT='
    'P_^4caG*{Vy=AOk<z#rK>v}^>SXBJOU`B0#-JN{l1LUP%&I2Wq7X5x3_mD3wQ~1$?;-v#QNMUI46Cc?r2j;FFvtpeX~c%4gE;<W7FrD`{OKvRRd@'
    '2Q{evfOGT2!6F6qqcxphkW(GzbJJa>zuKUOJ<q;q5PTmI1>#J=%Ul#a+|J7o21tC2D`2w4OoR4VsCv7Y+&4mvBYf~`yKMicsHI5#T;lMk#Tvk3qd'
    'E%a4`?<%PQaDR}m=n0rW%$^&fF;oTI`jNwP75bw{RP9#c_?W0~5n!$_K8-7-D}#peCIIOg4=Sm=yo7H;2*!W<dvFfGu;It>4N_#7<z%XY56XkQm!'
    'BSgxFVMXWClv-lUCc6=*o-oF}HoT>sE!bKzQ$;!b8oa3+AKhpTt3*^rPoMgrxXcoaj5eStrF1q+H{`xkOX6pZoo2bO?C@X2vU4kKEsBvfsjHede*'
    '#(}|Sz2A6Z7kkp^C;k-f2%)q{s&&gH^=}wTYqSBZXpj$mVPfBYhIGbY}bh<^C{)tYG%vBFKk~HB(#g3pd(K1Tgt&6u~?LQj=VP_I%v~Bd$yS^G8Z'
    'uRt6gFf+m^e=ZKTM%!|$X3ylo3Ag)J_h>EmnG740UQp{F%8lzs#JAc$9IPi6$0TpMa(lvFF%~Ju}IO)$<*2!%%UIoQc8(J*!FH{qcZzGMjOQLfs#'
    ')Ge@%z5g?(!iP2$2NJ2mr$O;+TQ#P}zkQwN0~P3>q6{?jtV8Oi2~M)Fe2KQeU^1;Z3}vj(ymKcVE?;|#YImW$b5Of9k*J{YZ*Q7G24$?l|LUsOp4'
    'uWS-WYT7|-+jJ)f1x5VwC2bxIYwxP<%f+ihoNY`VKs@JuZ%b&Xw(c+VgK?UgDy(hG-pB0?Q=F7}&_exA;TE)RHZ%04@H;F@oOo_&47B|=$4qSh6$'
    'I!n&;zHJE~d;l8`W?fLZXzZi^HBzP@?h5&HV{W3Al7q*^%x}5X@l}iLC?xisqe(Y5~ysI~tYcSgS6v0cm|w9%LI|A1!2i!dDL-birnOQ8D7`@kqV'
    'q!;DTH{i5#h+8MV!=pGaa^RD+KCM1iu+CBXR9Wp<(udR9dU~xhf<+TgmAs=g+OUY#;W3Y+r!-k)(B_JH!x`h>P{E}4nNyg1!g`Dgz@xb-Il<>ppy'
    '(c7pT>3a%h^_l{ZG;I8Ld}(A9IUOR-W7cH0B_zE@=HODX6nH_-orFZ%i;1JZO|JNYn9ZnX<`aQXEyzVskK|xP8$I^XUt;5Z3kTcDPo-ynE#(G83&'
    'pGn+NL$WYho$po^ZwqCCTgJn_bgo0FZl`M;)2>!<>U@Al0kZMGt(sH`rHNF@H#TdYrP)aH25Rz?vqLH%OD$O(Qd1AVwQ(!m~Gm{K4C(xOW!t@zoC'
    'mmD<ikE`cN!sCEcR-(=t!EQ{(fh$)1LoQ<q$YqEqc2x!pYYvq+EgR82$~@<rNI;I*HEASwzsSF87DX_U_pNcfIji-Wd_7z|HpwH{=L3Lo(;Uz|!L'
    'yr+<?Olo-UcAjNt{YYs(N}!V-*s1!ED$)Xu~9YYY%_5av)RwX1g~=p=I~Khj@xf*G0+R**hIh1n7J`6Tszbt)F*G`|~X`6O2mAhj)>j9M|n-h1%6'
    'aoGYOUNwv8d6W_Odx~5m9GWz4*Rv{^nvUIwI>bX`{{oZ(CQ~RIMQYcV#yzb!YKG`j9o+?DGIxLA+_?EO@d;qQZw`riBNLr+Y;SN`0<zC>x+<>w`>'
    '>C^r^w_?-KUlfCTrrX-F8p0#Gxf*Kp46$ta6Q5q5uSIQZg&E=sPNc#fmM*m5GAlJgHFQGAN1sZOoGP8u4d2FoVTrIYz!8cClPmg*L}MQR}rR1wr='
    'pd^lAq78d^R_v_{cOCKb2ePe-j=-GWlmOgMMsMa*MCo~;ix{)wVI0agAEf8h(Z$AP%6;d77v{j+{d9s-$a^1-3{&B`8-YQPg<JK%O>SC;*k3MG_W'
    '1pBp`e*XkJ6v%=h#YJw`l_=wQBTF~9si9&>D64V*P`mnN@+WeP$)cfoVtB(nW(KTl1Gff#w%xB1eU(qRquge^K+wytyPQfo7ZOU`G%1F+d0;55xQ'
    'pf{$Hux=I@^zKp76_deM{dtCPhlPnj(~J|7}B|iZ@9Y;dNBNChyf0tea%7aH~ch?&bfmr^jVd=H58=aCVIe+E(E+qh{TU^Vb<p@w$`}Svc7{gDC<'
    '+6gG>$WN~W#j7{BPF*|Y}279kx(QvX%wxW=!UlE!628i?-mD_exOcwlvv&>N()X<fSUdS&owjpl_4)e^+)u<DSaCVM(B~lrNUkJqZ=5nKfVXg5TJ'
    'QHA0WX+>F3#<ivX8ixb(&L1bUhq+wrOQfn@-vokYmcXpPZ&QhQa<4y1hAo;C7d+(`C?0%nfNvYA*kfh$ks_<_3e33vVK{)C==|AzIGTw$zCLoDGk'
    '8-)6cbVqxogAh)154>sqH`42{3RvsY*8_-jRJeg|?`PQ3VVlFNmtcqc*e0JLvARamDE1Ruy66-s%k*b*@>Js86};@3TVgMEK)SNfItQVu~N2Fjwm'
    'sI36-co}8*f_+7=z6Cm=UoY)11oJ1-dq%UdkCfDnCAIH|X^Yt#LR}|gyry{+vQQLkCK_A3Jhemvp>Pg!A{hBaK|yi>Q2t9$AhjRr+x+2i!cT=Op8'
    'u-BvmnAYnXL5i_gmtH%wIS7#+yra4RUy<R{;PtYs}DPMge$RWOJTmW7$nsw2f;N5;|i?FrU!Lkc<U<-4qhTZ?Hb)S8#5=i=se>>JyC)vG*r+O5}o'
    'D+(KcQW*5v(V!sxx8zg@xKgo4$C+uj-8wB6NHP;xdFn$-M1PHC>4tuvOb*V8~C)SQ}Tq;+m$_w!>Gr0hnpkbNfqS%-cc-Q=8>|r$S40zMZ7UD7gi'
    'g76{P(!h|L2WQ-wHPR`aTHev82^Bb)nq3*J=9E8tb;(lM+81A&M+W>mSbbVDMK58gFHO<xLi<|>}4<aFC5PLTqTwx-*?wqvpGORy=liL=rlp>jln'
    '$g^T?$?-5-tpaQ-2gLi1%L4xSc0{25PvTS%^1Hc@%;x6XEkJ6dzkFXP7oPZ#A(2%EAoOK2bn6-7mZpQT*i&5jv&kkWmEk@;!q<d!N}Q&sey9=FLM'
    'l@bGHJXa>V+wetzqrL(;6B^=R7N6VaQjCe;;QOAji#!}a`x{R2VR+BKB+=R;o55y9B2SpgMN$JkndW~Qxk{u}A5b3<p(+0ZChjSg4=1|mN?E?LZ%'
    'hE%bY%^h=+5l!;t+@+ZtJS2L8X^v`rDdF*=HKf*+c<Go0vR>B6^D(-kQFv@)Y`!a5J+cwP1^==WP5FA3pd!KIr8IHYNcIMvyp$?^G|y+M--H%<V^'
    '~DCxCM2<g9VnkCqlTm&q0o4@N8jVm9QWQB3yJWzYtyxq-Qf38NhO@G{yjY!r<Vd4DV7@h*^IYR<@v+}A1h_aq&0&?uE5{t$D-tTvDD<v~;ARp?Hd'
    '_hCMy}%wqn2pjD0`guJpz&i&L#KjniS}m!tU6PYf`%HlG=B2e%Q(*eMm90~Wz-*mr!_Jj4HT2J2-~`{neiGqvBJ~KALp+p^k;wlXPS1|9YBJkc6l'
    '#qqc#Q+STi3BUJ$C?sJL>3QzvLbmVAhuwp5G{f0kwez<Dfanc1v$TT<a2S+-9WHx(@FZ6@OHjcT!-JR1gx83aDUz91{Xel`zTWq{<9xM0*PEr6q{'
    'K_VyE)YB&PXz*E@6qYhp!983zo}RZ}?TAKu11}{x&kRXvXU3Txc7cU@63;!@ol8`Y8LT7O5;<H(eYx3cFL;chPc0rj5bewi`sJ>!TH;baGVi~{$d'
    '0+fd>G>zNDlNW%^}!OC29p_`FB6nqm*ZgH5ioogJQbi$HbYtE_iWBYSh~tcE`L@*h^p>@rMFi2K*>~o><W<x&ON2+1CTD2M<P6CEK=sf_zFu(tM2'
    'mJTmnTn(FWJMRdsA{?q2QEL8y$k_7XQngz*+NrQnm=jH6W&z-Qk;u`oo$^=CsXLWhV4c;;+)l!!EOPgdB30q=;hUqZBN>*2LM28Qb0cLQmMcn0EN'
    'k{jV)w!V#;Tdu;?t%0IgNZpaYMrbA_(0|aCU#i8SMD>CZxDf^bF^RW2c-tJ=87>K4j6?h{WjL|%)}m@6W^WNxa=Vh*KX{GmM)V&Y&r1&Iv`|CvGK'
    'ePG!7Wo^Y3)V;%R~$fhg1#z-VI9J31&xtxU>xaLHKS#v@2X635D!yJr-?MIoS{l~@>?zeYS<A(9w1;1USqMH&CO4sT1&?D4EQAsZk0I+be-R1B=r'
    'Kbp!t!#GYUmKb6Xlt<?-weg>;hK=*kc>fLRSU!9j;LT&(+o#HSP*x5*<Ls%%(X;9EyJZmH{T5>%P?T7+*2k|aOzRN3A}M7hy|W1^WDz|-g9zSGtD'
    '4AfzsgG*wL^w?1E~+%3sQb~Z|^dn7D#?Gp^xM)3Z7RpRlNe;8LL>8u7wA^)9G0QUkO2{X;JS*J&l|nfYOK<7H)zZiRR_@`g8{tC{RRF7_9#5H<(9'
    'ycLh3|#(qP)By8{8NGj>avW$$sw~{Im1?idP&!@5+Ngj7|a<R0rS}{d&HNhDKJ^kIq;J%kXr}6r`YssI|(=-TC+;CM}1NY4X=hj2So7zj_mpGOk9'
    'g7>_;&_*5Iwdm3{;Kh+{V)W+4Z682@$)A^3C)d^$fr`9AKs%N2BYceB=YJ}C&<)^PW^lnFO%xd!qd}<Gg>%^71MPC5-y$@mhR%k2dKDkB1FZvgyt'
    'DR&Da80OcmDWM9C16_Yc1*H&Amk_{#lXxa-L4r6T4<2((Il6I`9GSTfIGWi7<8^X2kg`*V=X$xuP-HN}-(E_9&13!{`z>+Av3U4~(qUzWYrTX5T9'
    'ztGUed2G1<aKKW=+41mpKT$<Nq?HC4oFsK2oRlEhGUG^gubUoY+?~lH!Di20XM=GL@wZRuaG1q6sWQb?Q+S~KZq449hB%F0LF-<6RS592sN)q@M4'
    'S6vhs`oS0fr!E`}8e=G+AOA<Y>iXWO{+E7u>$M|No0|NxTLui9dQ}kwgGLW{Z9|ckg?zid}oHT@%5fUx%UrMw9#}rLk~0P|1-%2yq(tteof=XdiQ'
    'LH+gYs)6cv3O}P6UAG<rao_IJ<CCNfWdvTUc95MqLMz)^#aM~a_G&sFw#I*sQwi(O7^4up#A@zJP6P+gT@Gra-j?PiiNOhGc6#Pa>kv7um$PW!c-'
    '~Yr8@vI$R_Y=soIi=~hPam<Ko=#V<XX}MSg7V5q1pGHPYF`UW%&Y|`FK9xvj>xm!^?H4JE2e1ox>kc|{NF~6xN<7ykICRLB&+m(KG>O!ppJpts<|'
    'kPC`y%a+$Op@W7A_C2$oC~fgl(bm~2a0+tPBBh(rK1QL>56>7(?Yfgy81?Z1K7yL*(&O{ROLKD@M7Lpx-Kg^hdgP@Nn>J^Is_P|#!Rlt~t7VBbR{'
    '0ZHE~!(pIW8MdMIcqA=Dvl4p)=f8;Rkn=RL7BYJx%uhDK??BHXwR1i~4}+I?z7@$Mck0WOTb#L-u(&K$6cXghl?D?JdiR-0id%#NRE6X@jNbQ*9k'
    'Nna12D;H?a>4A(0TKRLCM8Va#~^)GbK~%aW5S9q&}S7ANUw44Po8VP>McLIO8F@UCQU=jhIwu)w`<)ab%b_m8RqF(M%$s(=kX2bAv7;QltPVd)GS'
    '(7CoZzsiS^J#M^bNNl-90orh~>Wu+Yd0aRw6ao>|6)4fTjMbU~MC4ppe%k7JSayiZF)a?=H8Ded@iI^vk$QMaF@1*rcGjf?d6e#av<`@#sX{BPav'
    '5n7K!&W1`)gx>t;`&l`lbmF0<lzHO=KFvPDF>=yn&S2>03ngjemR3R*Hj(c+fae#zup^qY=3KQ4_#Q+bn){W@p|@HFXjU3c_2Tp0OOAPX=}Tscyb'
    '+Gi?8-m4mwi*Zp6Wu_az%fA&v44=Iba-%>?=Jt%N)hp{)$<UaLB-Ce$||eZ~G$54n<<T(ZNCRGJh4&3{>`A9@6+PGEX~Q9KeU0U+M~YVw4<_LJM?'
    'i`5`Y;I3RA)IM{{erolfwm#HK5IzoUAdpji1w=%M;mwaRLC?`qJ>_DT@Gt0HSQG6#z9NuPPNqi%9i?86<Veu9P1BOC(rVmU9Wh#RK&N=K-mNL42k'
    'y5^Hk{Buvcu'
)
PREFIX_RAW_SHA256 = 'dc25c524d10de6606fdb106c22195edbd4ac81f24caecb7b7dbe0cff3d9fba5e'


PREFIX_DEPTH = 24
DEPTH = 48

def transfer():
    return Transfer(225, 225, (117, 270, 4), (126, 286, 4), (144, 356, 4), 580, 710, -14000)

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
    if F(cost['energy_upper_1e12'], 10**12) >= F(123271600037, 1000000000000):
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
    if F(entropy['coefficient']) != F(204119948756637993910827, 40124289454043315000000) or entropy['checks'] != 72:
        raise ValueError('joint score or dual accounting')
    return {
        'threshold': '-7/20', 'graph_cap': '123271600037/1000000000000', 'coefficient': entropy['coefficient'],
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
