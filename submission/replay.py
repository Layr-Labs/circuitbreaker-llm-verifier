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
    A = F(51, 200)
    eta = F(1297, 4000)
    beta = F(3133, 50000)
    theta = F(36961, 1000000)
    lam = F(1031, 50000)
    xi = F(35297, 2000000)
    kappa = F(359, 20000)
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
    p_anchors = [
        (210649, 210649, 210649, 368053),
        (218427, 218427, 218427, 344719),
        (164687, 269203, 164687, 401423),
        (173536, 273585, 173536, 379343),
        (158912, 207258, 207258, 426572),
        (167079, 214288, 214288, 404345),
    ]
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
    q_anchors = [(602780, 919180), (576631, 915031), (715545, 960048), (680946, 959621)]
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
    cap = F(123153850,1000000000)
    coefficient = 1 + 1/(2*cap) + delta
    return {'checks':len(checked), 'local_reserve':str(margin),
            'saving':str(saving), 'semantic_increment':str(delta),
            'coefficient':str(coefficient)}

PACKED_PREFIX_FIRST = (
    '{Wp48S^xk9=GL@E0stWa8~^|S5YJf5=Znj=s9gX>54>OiUE*Q<rM5K0Q{o0N8%)9L`d$ln)0vtGzwq|QE|W&l%#zk7W)0o)>=7~4QDz|1{4I;t9(Q_*soa>~eC?6JBQfO9YniO_4n%_eb=sd0ac_BSCR^91@p`e;C+gFm$ek^f8Pybz`PjA!g5tQ$Sf0f%xz*=K&#D+rRaCxgl7$ah-vomz%QdU9!dhA~puoh%(ebtz<1uw$FRa<5bJO_nD|%=3x<AEZ4%<j6&Br8P0GWI`?0xdOT2txZA&C2?*et5bc1NQIBJ#kp*61O)EZn~;%afUdQNMcGA0D`O&xfaHI?2k}=VBX$!}FE6Ld&z9(}dJ`ABpL`5BnT}`Ld^DiLgL;q(nMCO}WEi<-T~xlwPt6TY%QkvVtjLd`#iKG>{Xc&8mE^hnKFWpo*vcF!hlyl^Bcc)z}l#JV_600UN1<p9kWqU8NxTi*xVeSfk1Hn)Xc29K`eEg5h(#OB^6{C|=lTmNR*L0CAF7CV|B{Op$k!HezAM5-Vn~-GDfXZ9!fB=gdyYBGPqC=zUF^;yi%+5hs7roDG1wHPjoarkP{44#bly{9=IxBZq8A%c<Z_e)S{?oh|b9gM~}|LA4Wd%%fz3^0nw?A(0XAQ1#E<>2*>X9IHQ8RXZ^TB|VI9-!ApPD_6s;`g${L5Z3p6h6Lo%Ks|cxjENYcLN9O-NNEoju&1!qO#xxVSIa_x43Zo0j{Ea=>d39EMnz?(c_{HPR7C*-2ZJBWv1hSC2q{glI<>_vnFr^$8pz!f`Am)>(nB{hBlWDHQQH=?M`<2%!w6Il-Xqzd7TXmKxb?~`QPifPhl{3?z6<E{xw1o9{!-L0IG1{xJJ&G>`3*Qf^fe*CHMP&oxYxQo>i_4Zs1P9_I+*4iO3mvke5it;Ek|0pX4?Vb#Ge-_osS>@S0Zt8b7H6ij*oxVB0r=~FFT}!2>spiX|T)P5%aEQKSQ9vTwD_iC17`R*jP3vX)n=uPiF*0NE?KkDc(+`N*@Nv>*ndFMsq=gK!O<gn(mF7(>t!FR|rwlX}Bd9o+8B2RnAoDE%k*?hM>=1@UzpW=(y^U#w4!Z#~13XgUX|HHu{GVzY|J@p4(XWE_cLIA#Hf^&mWOSjE*og$bS@=f-dDP!qxArv_K^r%4qda#q3hjek0eln0ka{2JqY=z(J|l%=x_Wb~^CiX;=u7<<GG*<taC6=cjGr_v2*i!=EICwU@SBg~<0<U?sumbItKx%*!w+fmOq3jW2ewB#B2eNIETK9QF@T_!pf8{y))@HG%tFt1W6p`FBE$acS8p<a<HhO=do(`XMPct!jWUAue$|*!*7b*|<xO>Us{ioD$+6pH2-zG~)|RlmKAa#%>+X#CAudS~*^bKZu<l7F`s}8Bsd}6lx@14TJQcSZkx&_YeJ&xajSvZ{V$z##g9SzoOex@QwsADi@98DNHx4l^Nvx#hW+<d3U&6xBTHo3}*Y1oS5D8y7dc8(I{X03zMn0RXv3Y$n1d4ER81bhtS+>@iQ(_gi|fN;>pN@zI#vP8Ps&U0#i8dH&Z!e&HLu5dCskCRICrAKo&ZUZYP2E`sB^E@hETWmSDW<eWy$ohm0DZue;FQibn$8mMd@qSqqqPY~8{oO2NnvhdRf{-QiksFVQ62n(y-VT}aPBZP&8R?sT+Km4FJP#d+@xxljcNcs`$JjGEgpY#*QK*CmMKUZ5{;)bj6hHiu5)H7unBC5A6j7%@tp4nrBwJS5z4(d{gntRZTl98SV|5KjoAxE)u@Rj2r!;S|222`3B|o-@SmEA0b;LEUedgFDJ23r`)ms*jf~TpGH3!<%*|9HpjAt2cZb>9w;`;&6EQE54%)xUXk6R8S{KC;fJIm+~wJxW-dQNhc2V@H}EgD_8zCxWv8#_7#>>XwHfuUioz%ZUZi&VxZIrMwN>Su0t6eKW(-Pbj!DvtNeV?YFZ=or~vL@-DPSBTI~@HPvZKZ^(AU$+<DLm01bJUU$Vx0&Foyv+U3j(3uVNpo8;~0NyvlOrS&5w*)apzq|n_5LVh5<<GK@`RiOkN53t?fk_rt#*ad{*>Je@hDVDm)tCQWAeuZX`7*D*FN~-zuh%aDK#6&S+%2h%_K1e`UJ!~S6L$R{~%9wN@6XOke#n1;Dv3U)6#QuWE;A+pS;o|}dCW`bEpQ-iZ28&>PPrXn+EbRtC&i<c^?M~@o{nZ<@>r%a=gHqxbBgq1|MCGLK^Yz^q>2~PTsQ@5JMJVem!|%+elr`M}O}wn|dZsr~t{7d<;2{r#0=uWeqjnFH;_>YX)XD!qIh2UgceC81_h?YoE>6_DG$Rv2Nm?SDb(I3QSGln`Y7^F?k8WxO9>fO1jmG<RsWD$6+V>{YL8;vp4%qA$Lso6XH4cZlI2Pc23PBPEBnMDpE8@SU)S~-jX;TuE=pEItjy?H{VYD;D);|cEC1Q|T9?VHYX2pQY!jZ<IhyswXlPC{Z@0;Ijtlwn|3_!P(n|_=k0=u^We&$cS{F*3UFTvFA>}{|iuocW=8chu<D(L2aEoqvgkt2_NQ!+D6;BJpFXk)o38kOi0b(?+Xnq1>43ntl|I?9pko9)GfF0}8_eC_0(-iVm7VN1<Jl-utICy<d*rg3T|TBdfK^HPJf;v2A+spBoTO*p;y>_IXO&y{Lq=&u$O`M$1<jlVVTc?=h+i|K^006iMs2lyRuxoj436<@=2|0B011I%x1i@v6@O~}nFBt2h>fD_GsX+_Hia85b+`QYO(Q8)kH{CTbiblh3XKR@;->XcA2aoaa@tN;8?*)v$yb&vznP)^;(j0Am@kh=Pc7)l&DsE{VbEB??C<1+~F@iUvR&a>ytZPLBP_Q5CRy5HG?<xxuA1p|Ychr7epDtc+MUpG6_2jvKwnk?efwz>D1d|Bf6D;YfoSqIe<gbAjHP(WjFO+rJ>EMr*Ti<8zbg=<cFz_dqd9$e@j)c|U+Y0xf<H02FruK61(f^ktlC_TLn;dbxe43Qw{dt~jbCJN12E!`wl24d6$f$Q=(6$xVyyU<H@g^}SNQ!AHpRNT1HNzJb}xN7mhu1mMQA}{<u8ZKI(dZa@LTKpw`<1!hcvJ)9w^H_zh{WJ7~--2_w2lrL62TM5LzJUO8k&jIaV+(AFYK;`Ij#S1W#$&3;N0MgQCY<mB&9k_XyP~ThgoUU!Ch%_3RS(GY>YnXU$5QVwpgid2KVOf`;ge8rr}C=sxPkHVnOvdEB-6mNq%z8}Z<#1WeY{$T8Gj5`aq57sSZ~n#QRxVmn5}U7aLyUokQLZ=W%kWWpFTTk&`{D@NB(*;qne?!unUDAxrOAEV#KG|6A6BeGT{vs58s6^#WT(FYZ}5rZ_qH6CC|1@^)amYLa?`wT73v~&W<(IfUVYx)IOGM+(B_sEZjj_ZQxcneGuJ|gjFKj9pL<4GRZMqHmBf<dW?ZH9MwtUJd4E~guQVwmV&8I<e;fxyLtFP?>Qlx_+E?QQv*8&xShHIj359oBqtoSSpaaLha(F&!8MDnaxDRiuL5R9U|NQqzw7r3<)ue5_N8djCBPEtdiDRu$6O}S7&hEhalJMbJ3scvG^8_bj7rh(nQnQ7N$q>`JdW`Is=g%iH)lrG%(Uj7DPkaS7n`ePryM$zURzlqI$J5|Jrq(TRhiz7P-MoH$0_6ZQT9E#s?MHknM#t^uE|IOzPXDHgM<!`WWlv1RVirHs=SCqGj6^M#FDS==seej#5+Ooj@ngB(yBZdQHwscCaf86sy7BO`w5b+*F`ZHm45LX)88pqU9U6{tIQ#zSv??RZ^c}m%=Q<3t^{X1XkXR1H!$dC#}`s-CnJJhsWj59UTO=<neY{8jMU=vEGkwahO{y4x$|cEG~KTABdc7Z?ht@kPF#Ft>N*j{dXi%<)XK~c)eQ|rdC^$I_n>`vgSFiVpwix_R}bh-mb35W#LT(Y9K~6h@YWOkmmFR7^aL8}&k}QqS}jF}<LQtqm&uyY)Ew!2NTw?pyk$}g|KHXB_dY9q*ii&RV09kW+!^8fu_EIszum8@M8~cm>Lib@AhEa|ZRyn_9_~BNt`Mi&mG!;*Ep-<ntZFBH)ySA@$Jfv3iFf0PC|y^K2$hMD=UTX6MFjXdQEc$@zGKlzXB!oSiLR3AB6==)QH!q_s-tQ;sQQ^5`=}uPh+Ja9Z#0HZVL#RTSw%$CVbk3gDTFr}tdDe@-_ODqXa2%|*sfj6jGw|&d90H;l7WBi4*t}PA*lhX-e`PfcXU3E58CxmqgGYx+OExj=#y1{)Pvlrwj}(OiyRNKp?_8=^A)D<ZNn1K?Qrit)>*G%NEc%UB%G*<PY}KC8UA9FBo3=O1>I+8U5Nt!$0@IWgwV?9mWLI~C~ZTuApZs!Yq$tmd<Dp*?|h86aJoulwvP{gbLjh&(twaro%UocyCC@jE>L+U+HF5MJXtEffU!W4=hs~Nr67iaqnD-HMIo1Cr$UR}FWW_xDsJT8=ilE{mowL)Ja_<`K+k_y0gtSe2A6{miwpRuHP}8HF7wDoVs%H7##gs0I<}{l52+XQkY!Y+Y7hcUfe^a`CY+xf#NvcfqrKs?Vh~<XUGM3|fd<cPgi4yGBSvF{XLNf~F+LfJXK=D*>HQGScFWam8HtzBF-t%M3dk2w1pTej%xUb~__-gNzxJp>#{%+cR4iU1W|H@9#-h7Hkq$DKER-D#;}t2&F<5zM_#eq{EcQ_F*C6oYR-t#Agx`5*5C8jz9T9)`xcbrD$qk|#m8@+WQ(HBhW#~%rTf-#E|Dt?mT$jJ_9`Srx?Yx=%Ngd|(6Ay+3tH9KPgw{S{9H6kYu8z`pu>Q&#j0cK!bk)Iilu_*W9Yp_GnWVhxnT`TR{p{xI57$2+R!D{O^Kh3ddMI94VQFLO`TosffYA&yj1pn4O%T}N9HCISt@&r>z2Z>C6lI&F&-~5vi7};pSg)g|vEBsnNBBeu4$d80MvB?`3VtGMq*4Nr-2eT|mcHz}q)Y?&Q=NucY~*;{I+B5_Vv;YPoGH+vqJv2_r}Heic_nF_+K5mW>7`MMMR&gzmyv2^%u*}I`2PHm4i+jaE~L(qB2ldIMt=$LeUNAQLDV@Cnzt^Rtoa=a$Ui`@oz}QaR|;g;8pby96X<CYBBrKdJf4}9JUgtE&p@zo1BEdC?oxGkXWP9`o5_e8oQ-JasDGw28IbjxtOyw%>=<JG#sCiENFCD&XQ2enqBT@abOi)E{NUWwLg0Jnh_z~Sike;KeyRnMNNeugK9tO$hL9CpTYH?74K^QLpBM~Z3jGScu{xii`@M~<;p4cOB@28=@hCP|X_G!uZK~|7?BZKm^vprKR7_^68#TU!!CY|Mb^yd62_Hnx6D2zeDGZrMc#3D{B7ksyw7rCr>DWyC5X3FxQ-Au)yT+_y(gN0nq><(F6XfP|_BAu3P6a+r=E=b+pNOL~rA25wdfjPN6YnPl0?+#pE4o$rC*0#1sd7@MVm|qeU+P4p#uX9pLW{ik9ojctT2bXGJ*mhj%JG^}sP-L)I!o#kWoe&0L{~cp5ITmY?a$gK+M~ipJ3!gncr<7BAvIXC!Q6`{jk?Ii&P3tT^K1US;rg+l(T89q)g_-IKH-(+GexhHTHDqY-6t)SCg0^sqXB_Dk-e?h+yE`m+2Udg;xZUV$od5Y@Re%Vr<dxok}L&(c%0==2sH5ld#K!2FJn2dg$S)1ijCul#d7<#>K*8m$+#!L)hf7_>lx;MHDL{*a$RP++WSG0+A<mvTsvXS(bat_yNFcALhCBZyq_deCp9|Bu0<NDsYTs+K`1}r)K)7MsV(z3GDK5xgeGclKswOXjfSiRTvg8``Y>H}Z*%OI|2VKZ!25&u&4(G`eUnM1z{j8qr)Dw5ID_RXZr48+I6Y>tp?=#N90;UVQ`=azm2YmB)k9q(Op4_0D<xVZqX2n+@vv7v4BC7GA(;2o_I;u?Q@%FGX=yo!DRIosd24;9KyV#7l3(rtsP>mIPa6kJA0^E?{AQm9v)vzuF{adEiBZ(PuKooMlZG55%IY}dv@j{BIqX4sTai^1PB_K^uG78tDAc;SS=r4ff4|+9j4J6D5OYQ;SN45nK?8d}4`aM71lpTm6=+EsnhxY+c`knR?A3l-#`FKa;Uu!4g(Dqv2VAle^ocetfEwl%ZiSd`6M!{dc44dDk?%ep%ehIhtu}7q1JE0Nx=K$@n2%#ALhY1DO}*C@f(8Q-S0$B-!ehx91gMVk!{3_4@61p^d{y>Sbt9iieimN${OB%xbXH=8M(jZlbAl>%Ad=p<u=<bC;1$yU1V`zOpi>b6j@&{7wWN&V^O=JEO0`<%U-O0CJHVY=Veg<jXOm_~OGqlCc+JYq2!4dwzD4p^Bt;tA=yRPeYkdmf!Y}xBbfO$agH?+S`zVmdh|0IyY5vW~f#@O#PS*F6r3ml*T1~vRG!qzFv&{xX9Q?f9DOWCN7Qn16Bih^Qi`4F#?Ju8BVS}R)dAZMO;TX=WQD6l!F&GS`0*(%2z5DZnA?m^;W<O*Azqr*5l)-ZPBhy#I*)@&U(|zUQvzs7RNbjqYftQ7SEBFVd-`epR;7SZ>sj;n)A4ls+YIq<l!@^`N>?!>9XbC(q9yjMF**b3vEM6&^E^4j+b2`Nx;Q(`%_ZfO^zRB}T#R9Ld8^R&fk7R_!?OcNOO`}M9^;pRYx4W1#T`5K%qSE_KMb;F6+)*SuGuD9fCy8OoARb~PDv{vTam#TN`;R-wFrQt(S%uQO=uF5g$x;t_fYbNUko&;--~aef{KUvt>kN^w(DhK9Y_pLf$db)U3t*$8*LDKtiTllIfRZp(y|Y0^mgk%P-uryqx`n;f%ASGP&(+tK0487Rf;og=&xq|Xdd-P~HH#eXYrZ54tlwVChEEn|-hxS50o$HKy_<$($Wj(-yZN>Z_x%o*TIT{VvZ;7l#vYV1tOjPvH&)tOJZ9FzoEZM$>Kgcj`bIv{tv+IX$SkFvhI`TjxUQZY=`w>>9L!b|rBL_I0N<`{N3YM<ffCsk4U=nHQcy?XX$Mdnk5eEAkOX@nZrDUKJN1GbszZVmHRDhpHV(zsJQU#o?>P2<49L{=)=R1H`(7?Jsgl$|=BavI>I?B!ZuiFko8C>&3~sIseghSnPAOGXz5E6Os&iMib+4QvR|;}8eC@;YITwSL7kt{}k(fv;ZA3QW&SA|baMN{neSo~1VM7>6>30?-!9E6xZV$CC4b|GEg2I+aXdqo6+BcGJGT!zNA3wQEe)%Jz<be&h%OM6EWowAme2zDttK&?^tI?Zh_1)rU21xU+<|)K-IWyi&2S6#5IHU9Z{QR_Yf~AXtke%9sW~wojJ!Gm*_jO}$$;LcWK@prXa7YA_k$>gZ&_=%B-luT$l8tGwvcKEO?POyX0jV#vTnF2^;)`Nz0kL$tb(_0P-)6;5Lm1o~+Vg?NorqRGBS%U?pI53!BOB35d8l8F(m34Vqnvh@bc4Jx6=o87Dup^Ove7ke_x8h4c{h=Bv#o)8yt*PY4uA-mW0xSz*75*$Xz_C_vgZ4P)Tq691<ACz>M*DZ60{<9PY`gZoWuSDL+o$|wS$r6Wv{zrB*A<kC=#KV1EBC^Jv0+e#fgRT#YCTQg2*Gza^$3+*pitE-5fv^bTyD_uNR%CBMGf3UW;!Q$MYiX+ZVl;hBf9sA<;`UiyIz;lpY*#_3md%3hp=OC=hjEdXp+29$Wrl?tS3nOT@T&6(q?dy8mu=Ku+$0qSU!kuYE1wZx0uQF-aOCcusDtm)I4rs(SM>2dPnV&eu0!H_ZQ$_0PtuaU!z>^;2v0ykKIWQtdeA1x2hw3=gFm&Ff3{%KQ2}xXtm>{w4hGaRvt7Ehu2wz1ZT8+DT{m@%<O-u2?@k5XB?nA*<aD1i5Gr@~sap)><1p7{_N*6ntYPvO!;m*QizXD!e48z#i*yS5I_Sx|q)iC_?#vpFD}Bx&;&}cju)vIP+J2k4~>q7>YKowqa~=fCVJL<h3py7OqSg^9GmP*Ua>gYs3$Nq%xlHAlqRqsQIHX!Ee^s*8d{AYj-qLdp&({X4dT!x(3|qBLB6(kZVHSORiO`ULTC(#Glx3xzX%)ndOCtav4pQ9td{BY0k4s_meQ<(FVbgH7XY5?QwCi)du)~-xbMiX;o{$y~>ipFj?eFyOg2pLv}bmbk;K%>J1@-Q6{25|Cr-VscR-miQuSW#wOaNyegPUH@qKi^yQiWt%Sh<nitJ}T$%o4i@3#?0jj~@{7X*!)=-0myoL@F+C!UsmMVQIA)hRmZI{Mpw02=Ad$815-*-eoJNLu+(_Wjj2pTwwRn|^LF4xNMZe}{M)gx$Iqjm+!6Rxg^uH7_w_~JdDut+8xQZ)+y-cNI#)jkr(sAVZ{O{7FD4Qp-|iOl4qZ>F$s=6F!cM$J#YqM%Z`L-+xr={)%HG@Yt}2pmyh&7OHR<G2{&+-@?JgrnO5+Mtd<Xo033QvGrDJ^U6MYg34Y+He{<^W+^Lz7AdMy=UcmriQstK<f^UD`}*z`fh~$RziU$_rt=bDfMM$rSWZ^tpQiX2}aU%(e+d_ei;`I9DD#^aw*HmMbz)ZKC?$OfU@tUhtNd^R7ufA#hw3R?OvcY!a#3C1}%Iv`vTL!dEb@jMBygU2LiJ*oVp4)u}9D%*3>O7`#U`zAsM~@hv?<Zl7nMpvy`Wo<KcUQUv_gC23JyE9hOIgSdY<dOO-zzq|7jZ4H#`{A(AJDuaJ~IZCRD@VoUtoI09nbLVgf!$yZE;2nIsMA#1*dz59iR#^%ZhWP1lIWnZlA?9Gi3eY)7_D|@C=>?P>sc}w%F7MV?bE(EC6*p|`HxJFkOyeQooN9fD@S0O-0{XnTaRb;0{6xP3qcf;x@&RVKV*KlUr&N7S-Nkzg4rh><kf5=9-?V!^wUgkKDevmCEY;b{ZI?~>^FDL6|oTC-L;h>+pI43$5(yZvIbigU?DHa`r#9Pv5PMqBGEkEz`{9-)ELja-hJtL3@VsDk+l;(KM8mCsl>go%Uw7@(QAiQ2nQ~gBBqc1~*J?F!>lCmA%>=X(KUZkjZX(tDb(LJ)59_Nv5Z@Ble*Q1T=$LuCKxI=Fpe<re7Op$-|f_l6qJd9WWHxt>!7t1Am4?J=E-FLHVpFTu^;}*0DPW`f`8l@QW6=o4oYE4t$X$v2jxh^70Ovuy?j=2%>cl3jpkJ&&F_V$2Pm_<T2(i7t$(hiRyK{c9KRti=R|2(~Ut>}K1&Crlew(fbhr$pAjGc<O%YHi}O2(w5$I>hkKJInIJ_UZaDv%!Sc-Xd9uVPh;_<bSVB2VEJf2!1q!+9Rlsaa#gBgJB{&cXUrYimN)@=76U7LZ?HHmIwS+=Pe-uwsJ4PE2Dmn=|Rof@&vf0zguwjF&wIZ`BUIyxWnyk-2S1t1bB|`uD6x96+|;DFR_#w-;DSt1T5s@Qe7Z%LJymQs?hZ($34Bd>yK)vnOI3RsyZ26#9s=YYpu!>m|r2jkR52ab}kBVBQY{kN`cSJWSDFjFcH4w3lf7VJ9M%8b~Ib!Dwot`wq8w1%FkKPNO^;fO>+cPTH_E-7|?&hy(V;I5R9{11ouLevvfH}Ufg@QYK2RF)QpOU9PiZ1^rVvr5(c+r4)IhvTJW*b49qfm+KP-pAV5Acu|q0+=r1T^r*oS9@}u90SNY<Zg+H}d&00qhnV46+_3=}y^IK?_;0;B)xTJ8sh=Q7)h!4bLGIfYJPy(789#@6bBS_y`gDMM4@4INQfx(krZ5!m0UHjeYeMk;#_49iZit=p$<9hjT1_0Bd7jayK`R7(~YbDhtk-`~lch&MAB|_Zi^OruSiRNCBXa|3qNXM2X;O3CCO-2!)Iixtef$oBqzZ)e1x(6zU5-noG38+z%>CG0FNL%=6fW$)Y56tNUSe|%NG2dV3h>DtTJJ>YN#OpJzA8)LU=RWO`<z1qHd>MhVmP2(8HX~nKTu-w_vWlBX!j=IKdu64f9EG`T!Z|KZI7;!dVZD3P*bYxX@s%urLO@6qSN0>Dpco_4m@M$-Djio3Dm%sZnKwiew^<~Wkm8UL3BWpPF#9GtJPgH-<c&W_KwQ7-6<A1@A?$`XCZgu8h-UJ|q3<2D)cDAB{x2%4h=V=w%6C8vQ>RiHQSP7}%`1vXLp@UZxa&Y{4TX4BO!82ZB%4{JsX=vQDdYs@3Qua4sU2q>s0OP)@ETr;VOfy98?Cj2_-kib5XPjD4Nrz6jCHL%!1TNN78Q0gDEVo(<G9_<-cYRLmbc%P+W!>r6Kzb6J9C1?@^)9b=oG~iS}cJ2J}BH=;V_$#zB*IH01#^O?A{<qW)?&mbEzkpU#}HSsC&&3RY+?Ov9*&@X#?rsuO>bihN@POdAMWDBFu_lW*C_|!wOl?_8_ewFgILE>;7O`1(|{>!mWTV)Y|0JPMZx-&0EhSbcSX-0~Hz|oMyB2`oW;tsT;g0=+UmV6%_KG1XRvVCIn9E8{59x7mTYaLI7e_Yvbm64MQ^fgs&=~a})PaT_P5h=1s^wbtE7(L}i2^MgEXQUi-)hP7bTc^{~2C;-Xzs*x~toj0AoHTY*JaAwIr6G3gpQR9U0+n;nR>QaT0Nwsd3&=JUk=WA1wW-mSR#GLu|`s<4FB%Zm#k2fu;q?w~KKoROTfrzQ8bXBNc;yI5f?8pm@^PT4(=MN*nAl2TT3dH0KB`@2iL@!MUZJ5k2hy-+G;>Dlq$9)|9ql9m)wJMTCH3K0$0PHc<UMV%QKp}NLPT)m0ZxOX(JPnSxL^LXu{wivV>3PS;dM)I(XZ&M~K$n&lq9!7!1MIz!QQCYs$Wgqc@K)CwtSoN}n`%eb(?1ojR=uu%vca(0NeEauvkv{p99ndES!MdZ`R6|%!I_{OKVx_jw<@==Dy*R0<YO-RR3t-P+%V!uRq+)^UP6M^Ap+|EGoD+2E5%{(5a%FR)eS!Zt3^n+khLT-~=@f+Rn=ZWl7o@xki)m>I30=3+1jR)Y0olXz_`D~c?<<Q|9dX1VK0sbbrQXN!s0*?!=Hh><?#SQX&rIkn%!KlTaQ)<>x;k`ddhHzJyF<BUF^VE^papL@C^*t2s<jQzCybpHr?t&<-tPz&ql{|hWiRd6aj3`ATC1*WK(j!Zr5y<_{ZG_w#5?rm$P&r6bNyclEJ@*A$k9$U&NF5=hLH{VQ40fEe*adq5OD){t7~8oAn*-mTT)TpSM1(#iJ+=7(CM)@H4t7QbXMR51n)4~qM_N?2G|hfUtp&>WcPoL1Y=~yl;HSNa?VP97GD^T-dX^XK1|EYp;i_#%gjGw%COK^dxuN(U;0_}RLE5EU$2}LA=qM?q4URp6V>sM%B)a)>yxiA%4b5WeZ|IJ&*%+`=f|*+pSxfB;f@sU0c^c?-G;2lw{cG#s=;d@HGzj)pn?WB>Bz?039AN@b1y>bQ7E9%mRr(VQbVxOp{M<ig`ONs=TlciVb2ES`p352#f6Jy+~q@K7SBE4_L=DVJ=w_Oxe}t^0PkYM>=0@(>i)@ak>i}VRnHkhPt4$LNm|!=_GkyNle*bSPWVys6^HKr?cR+~fGQy};O4qsD3J-(8{EwI#kjTxgT27I7AKe{G6(0z-}UqBhJ_i#FMkcai=&8t<eV+utJ#ebC2I(R{Z19f7a{kdqqrKM<MU5B!ER<Kw=&)nlN7IR7z{ceBG+Sh_VvQPqv++(L(}PuzJvT+ONa?%XUp9T(KP;*#v7x_!eX7%eK(`*)3uU}zXA=5w|~;P39Jo{x|>~TD-3SS!jK|{{dA1Nzbn{|qU|wes0<j0pn%HxbtBu@QD1pWJyzW;Ful_+>eP~{a6o_-Slk;|aJ#Y{`J39x*2*nwDi^Fy>PGTM&!r10d<7CA^{6hg%l;Wnr-?`zS{MI;@^Gh{KKTfM!PB^&)pzu5=QTJe&8s6?`rYB!aFx3at%7ta^{P@P>#i8CU=g|Euj}XE&)K1PwYChE%@oegU*M(6Z$sU?=lp|qk(d0wPfqI2)2J>60bMu==cGjkFEF$_oL3cZ9ZjBJeAc4mHdY9UY)v;Ta=pO<RG7Q=r8bqCgsm~c5m@i=t`RmFb|7S%L$tD}WaK1)2MA<kq-(oWQg~M1xiLc(b8aQGc>B#D^{AftZivq5H4%*1Cbt)`)+nEqqeN)gT+h9JMi8#fF&@@Nl+W6(G4|+EzT)zyv?J7yociQ<%>0kFer*r*1O#jo^uK}H8l|1Hx((gr<r(|l#%)M!D^Yu8P)4NX(tL4Ut?#jJGh$?k_76`aLusnr|KgV#%(H(Sbrge8X(me)G(SlZ-qI{vA<qGJl!pbCITh8rC(mP^jwh^0ZS}ftP3HVD0bY7TlcS!gI2);B$w`1EyZ*<JKAS~SRV^AEQx>_@Tdkj6_MKifI#R(PY&6mTUkyEuXey<wU|WW;31s|!2KUWyyAR6o#`UCDPTq?B+}6$J=|@_D8c@EistUXXqy)vjVg7q;t){0uQeGbncUW^X?uh4XeWAUM;kR0%-Rvr)tlM;!47$~`yF<>vKg=K)#>i@VfEcoLi2O4-(mLfhM;1cOix)AQ6>FD8E%0=bE+n`RVG9nR61cPR$@af-MzCT}2VXVfG8Jn{C;<Su9(|@|V_<`ODlRwb#7n`EcwKm(x0TdIdlFCpjl_pq-4l+Qzltfdw#@akWGQxr=1&PVvG~Knr^#aAek4QbjZRKjAcR!!Q@iID>ilH~WB)t2#GR=xSUlT7$$h;z&<-18{}!F)ESld+<fK&+i>pCqVaeZ#EAa=Gc&y*lpt&d3c0RjyIt^Aha83+i+({KAl)O+0GmQr%i%YatjHJJyfe@Am<<WxAkS=$#Q0g~7(wq8v6)BLgm6K-#x!&aNV31pSsP4zj#u9vY+BX%cEs?#PB<FY0#}A}FI*H7j8g#E1>Uwd^FE&)>NC3&<RwI|V1a0Jo(!Td}P=OkcV5sjm`oxv?cGd)F0~THMC94xtQY!a>2!Z&Pe4a`*Q;V1nluHd29e<^HrWv}a`Ji%+3VJH-v1@bL%VW5zOEU*zZGFAQC2VT@y1frf@InC5(frh<E<Ti2iU3sNqhSXyfHW*`KaDKN99~5OySHp=cP1y;mh08qtHk9!1*YL=eUO<%pU@)wqrCc`IcnDcsoIgh=^|L49u-$`>Y12}A72KVI=7G9=$H7m;y{bcRDR3(9$gaG2^>!L<vYN*`1R>A2~-8p<abWSnU1TM)0>?-*x>?)2DLcz#lUp-W7T3mCkVy=0S;;f2hj@t^Vh3GvPHvcZ@v+d-qY#L$v)gH4nR7aGoWY{`@@k*jZsCmM|6bB+v5*o!V6-(?$>`oF(F_74+%|i=XJ<TUU8OeW~H>3g@c-B!<-;YiV7}LR67ei!Mcok4pA#;c<R)O!cT1GERN$|mWn7M?X;nn4XTutAc)eergc*?{?`#oz9ZZ5;pDS$$+RUyvPta+UQaLC5Lh$+(Zjuj?&<|gXsT`Srw@YA5!yR#U>Y`k>li%fe2q)NY4tew4sP#K+H7dD!E;@a&YYf;hs;tjx;&q5=BZjk+4<pKw8QVl&vedK^Y9{=l&9R;`dwG2LORGTVjVbD)JnJ~9S+yoqf;TceCZ(D^B7ydcEN6x@Iu^xWX{U|Zc-}pu=9=(XzDrA;5bWmxrV1L+j<ya8s93im02}NMg!)i4R|~}`>bij2z-!;(qY!)+YEhoR};_aw(>z;Dc`#;!e)KyqG&kL|AtdY+~J89E{x``!9Y~8cW`XefTofdbPq-a<iy}<x%zJ9Qm$5Y207NbM}|se!FZk=Lq5!l=NE3b02B;Hv`9Np2A;L;d|R?3=vyGm@w>~Vzhmy4B-?YyD6UD$R;g-9rN(=_ChXZPSmaY>bxmjJ&>4UpzJ-La+eV2(mvXzPdQ`KX+RgfOR;ipg1|{uv*UXfwKNUbrc2!d(49f&TD{)oFZF*Kai(^RXYA;Qk^PT(QQ)^GI-tlSQ-AI0D_0w!ksQWE}HdIGiYh#7G_(meys#P$N=X&<k8HwZhD3~pu3dmiE-RK`)n{Jo1UskG@+5BkG323&>BMID?QsomW;|(RTb(8EKfDa4>0Sux)tCGEQ*dOGp*})|rvtkKljc$f(S)=RN04t^grh`&=yt>Kuba~Sh2xEL?E>U?KRo)g{*cbUG8r+S8e7;zb3kuy<L@|h*9cvUH0`UlYTCj^n1I&$zlS=Y%VJ=VRJ38I;6!MHaf;IL<!{Q!@Mp>z_yNX%|^5}b<9d6dpc+|=!GU_?%2GMv?<D2u6L+pq{Q5=fq7j6IhVgwzpScNn)a9xR@uRM+$K`gN{8$QDnqQ@&K-o=++Fh`NXZ`<`_3@3)dR!LYtzMazGVng^(HiZc&zM7Kj^RAO`{Udo$>@L`uVOvf{;4}a^&BjwYJV^nld1(3StvAvYcYq_>BvINr6uw=oJt|qsLnr1+2S^`7fl?a`ilK8t&xp1Z;BqwR1_EkO1?5g)l4ArikRZ!=Yss;O6Z<;LYjEXqEG}>!wzjyIRLbUVsuGi(8>hMcm~E*U)#^~WivM7o<DQ@eu5Ih*LMLqb^_^hE&2)Hgp_=>cC)J$=U3W2toV+yeF*IriUMpEOxR<duphw5Uq?AjBKsVt!ZH{C$%fTm}X~&@98*@wK9k4<PtQlO=lp@#65iZMRq#SI#aO650ErX{Y;Sg{27wSB0K*__DGIXQbRJtiq=N4J6w>->n9uun1mfW)Tl>EGw{nc8c@~s5p-`Y(Z2<W_W=|<fHOs>Bug$_{aheU1yHHVfdPElVy?(6;22wJM<(CDEL_Wj#KnMuu!kJXUoFLdQ+U|310|JZX99hU;KA=t+MSv)5_(Qs}~bG#DCcUA|xsEa9Msluz`nAo-8ofAG4B!qDJA*dW6$(ryl8M00D*hF27n~(Z~0Haz#vsw@MJ_9c^P_elaj~cIhFI^E%L)$BupEQl{-!CeJ1q%>J=Y2~J$c~H??&M*vD9guvpk2p9LR(R4+CnQ%n%V~_2T_E+hdb8(?E8>a=cDBAyHWW`<UKQqh2$1C6mM0c++<@Va*14m9(sJVr1#JM9S*OJ!+j40+jC<Nic8rp$S~o!M#j0uLku`lT!PFv51?hkEc}fQPtPNBOs8#kLyuP>RX1gQ=nz}U1F>v(*c;YlzL<zI#l;G2h*cHlq)_oU%oxGe<Ez3xJSsRNMp+Q^8fQ5Oyg(j`;P2~Zpb4|ckPK^aW})zuR}S2hUrI{b6(UIX9XWj*_i$`Ogqs^xR^wz8(oZCoMKxz+#q*&Yx%a-G_=`1RR-6AOWhzRN8=8|;SAxhyiof?K=etbU-g0yHDIEXLh65+c&Pe_9nl+62gh*J9K5rDI^a@r}AoFm^?wTyCyPj|Pv~Y{T5ZxTQ5te?h^D|SYlLov=L!VGAaR3Er!9uc4Y9Gy*7kE|J`Y>SjRH0;Is0)UTo>?gl{1prr>5>a2Gy+;1%Yv)(bOi7~{?q2=x1FFdbYLm0cbrMY>@Z;;x#PaYj0nseR*^*|;%Ox7ba3DGp!ca&#oa5~W(~L2iasKNaqw7voRO<KBv)Vnw8G;AF7<{;L`4wYa+7pBIBrV%`?+b6(K@2Z==5@9d@DzJ8(N<Pl!j-|I=vmBMfz<MkhvMopk+E4I9ceiAcLH!6buH|<dTpE{F!Z1+IkP3K03=dBmozsYk*>LCYp^U3FWxc>{F!Q6yP>-Jt17h-IuEJF*r~I(HFrky50Dvk@!|99WXyVvZ`CZe$J}opy%c|&PMOP$i2+&8~h)>P2;@n4Ag#bCKZ@4CZvK?+r~V;oMq3F9bxEHRbC5?I_SHhYz;==bQoYY=p3(k1Y!F8_AS6|JLir)M*4;G0p6ZOk%X;LC~TTA$P=<PhI5hK*eLij?gKcHNlFxS4E;W6VenAVbyIZ_F>OWvP-U$wNF@Dr$dK&dO^rY=&;VA(>jlk1U%Q2CL8p(S9WZSGw8jG~3ynZS?MQP<Q@45nnA-!J{tb~tdU?ozaJ`T9MQ<XO)2EfXB@q&LFvcsh=j-w8kwm5<VgmiVkO-ZxZq)ft)&<{I`j3?k2;p*X430)rNZ8FT1MVl$ZJKuWlp-Z906B7rtl&l$!f~sj@s+`Wia&SoJ7*-L1i8~j_rBQpIhIRyeK^oUYP-Wr3r3sGeYNUZbtcCTU;wk)QF;=6KAWyJCiQ&rLxuSgo$NF8u{*y6xl0enVa-r8JqLE5L@{C!IOW<Xz&%cek4TntfC=@;n@!~l>VS60P9Equ3i^0Pe2zd!e%Rj_p-y%sdn+UKslXW%fg16qiCX)9xke{W2qc&+g|HJPWd-!6<RT}?v!lb*d+a(DpHDx6;9QWo2loevbCXb=Xn@3>$h@jVj%ti1hO?*mJgxe@9VBlPL6<ZSn_gR1frFZqvRYZMm7q>P0C>#jEFr>2^n>n$0g6$VO@HiO3$?=zHmH2kZGJilWb%3#rgODeC`@}0)LVu*D=2t)#H>?qsn6M94MXaV?Po9dYBIJAS)X`8G&W4HPHQ&NO148vtwpG7jhimp8R^+yk9uF1j_#&(1$`<@cM9fuB$0ns4-L15l%U^^;WwF|3vz3mX~nW!%7^?D3}`WW5bjOk05_Fse>vYo!*ON0rbu|?=i%(Lprzc~T)D<J9;`c&46x-8-EAGw5pjp5FjvtcjOt1d-np=07<KAHYkn5xU~OcnY1^ZvY^ziqj72zzNw78{FZUhuWscwTrI~O(5EDXVZR^nS_c(%+54f2(^HD**lZ`6?+kf9~lJUzChkD2Qgor3-04%YP%q+a)2pC)l%w#{^)R}eDt1+%kN30lA531rk$Z**xLoq<;wEw<CMfN_khBM$Z(*LTFP}Ajn-IioCHjJ-mR!hmX8=9LSkxMKwx@u)mkyjAw<Dok$?53YGdt4ek3;&RmRFJ@m6_L=6k3oMqAnj4@!uhyLnF8R64J?Wtn-PlKn6dj{KGjmH-KddFe-btZpC6Y^TtBx~2SpxXNs~1{ATSQ*+JHy^;_S(lepaMq&};`v0*IuQA6Hdt8$;iZtL=&8;zM2lx%C9~zvcksm{(Qt<+G`Nhh?O5@|l)i&-G3C%OhQ=<7mH0_-9+M1>6>XBuml(=wbK#l(*h6KMd%`imFXE!=}cRk;U;<%`?7UUp2`QIp883z(bPDk+vFUsCiP&=~==W?%Y7fULxV(IK+E2J0^VhF>x1<tg6N?M3H^ZwMe9Ekt<Mqb^pL(+BG!<#GOc>^&6k?MwEt<)L)SP{d@tF{Tom%8mCAhojEdDbp#p^`u+uWuVJOsAq2I@b-rN`L6V5>o!cJEqZ+JMMII<FK6wbiq(4C)aR;TK*aNiy^VaY<HvderjGbBHEcXlV1jHXmYt&unFM~VLMVq0I{27%CbQFVanldfEJBiwc%q3M-Y=efM3vlz}WBt0cM6FLXr3THabfpQnhtZoaj_c8Au0iXUBGp;t%^-o@n^=*U>Gom={Tx1Y3~fl^+@W~gI$a)c83@9jb3*lj-{=sy?bd$I{>Wj9)4NQ@LgwH2s^Oxc>SmFd&yxae*6piLj1~H?61T1gOCktBL?6z8Q=MmElHy)#=Om1_I<aN-w_+>wjmgj7Xya_PuKT|y&MnFdv6`3a3zO93wKfc5TTrfOtA~gC;!ztkV_{NZgvb1mlXeB|`y&=VsR}!eJOs&7D=S&^!L~Za12M#Ngm|HnkCVZsR-MkoFKXJm6M?lUIrM~qR~%d7vbglLfCP)!KdQZl1@-Q*D<pIL7Gt$l1i%XW8pt}K){=)__?&22QjX*o7Y<S1dzyv<(<HiwzwaVt?mAd2N$c*P0Gr^*3^x4HV2CIQVSN>(ay>O8BR<%pZ1n8)Ml3lfZn=9^Ksv5>%NW3f+SjCWK{|~N)2wnFjiX5>&UHAv3>z?$!r};I2GZ%j`>b1$B7=ojQRfAix8X)@E`Lx~=)KZLzLj@~Pu9kj^0X&r$*9b>KVVvY!ABJ68T?=BL4iF8Lb1=``zpC6@v~^^$_HgE#At`!R!SD4w`bznVNL^#T8kWRUTyi#H7mh6dDj2NBF0=rp275d2wgI%nN<D&EQPU#WAt}PI36}A7$(SUY1E~!(78EFsTG4%*{$t_tjIANHbJ}iH8kiCLI&`(`QwQKq#3y++NlE*XZ>7;4Qx^+HRwwc$2hwfLfQUIelK3Ho~BSRu8ewa2G!CIj~eZ*W%lWlT1oEr;YVGI1X!T>qcSGUv<Gf7Ee9^^Nzn?XG4k;|T6*L?wL*-VxZ|Q)M)3i`hNuyfP)-TD&5>gKzOl{Yo>VzvtlhZP&SnS^Pf%!<xgcLKXX|5Kx}}DBxBkfdVPJpI`*eu(1J#@*-G;vU9l4j+hm>uLzH#O&WvOtBpdA4gY6{cis+}6sxB(FT4NF0`5wD3PR+o1Uza8&dp7G+Qltkp`L14A+ePCFK5WfAN)jTj<zWnww)*HaugPq$KCi`=Gm0*ZUfiYjxb7TCs+tZzp+NT_aIMq8~U2@X&n}KV6{5Z25C?CKB`T-KuXEmP_*Z4Sg1)JEsN~Qexn`!>;JS|pV5#9Otcu18_*)oB?HTs}zcL6t{SX~2xQgL?1jtyAFG)bWPXXPB#g034mw{O;if8O$@(ja=jKTc9f>ajr^-o$rERgW^hs_UAUf2fT%nf_w<V6GbTe5y~`ICOy1Gxk3%)jnWA9*((N!wYy}E|PbAOvay0ic!y<uSR~~LbmH%K1|+|RKjdChWw^X#kTE|MNb(Ev?)&EC1|~e<qcK+!i2B5F87Q{Mux9KZ7z<RBP{H9mTp&dVU%|b5n1hR+E(%RkAp=Rl#yr|{^uJLzlx;<MG3tg02kHG%NCrJGN|QsScH#KSk)DGi3z};szcM^KYXC_4qe?B#tmc9LHNIxf3{-qQ2%ih3!wkAKSV6mUZ{OEKSe^~YvV6tpzQ1e0hq&oY^))EVu8aPQx8hOQ-CTKM0buX@jz8s%zI0hD?VA3gTbpBEQw?(i@*t}BAo##4v-KdF%4OH=lZ(HD9UIS_sp&SijYi3DvMRBv@*Lng_c`q=eshf&d#?W+pm#Lsc=jn>lY$MnR{}}3bdsoizB!1#*!tIA?q$Be?mn*qmlg>^Xj?qN%(B3QMzg-W(V$WZaoTtko4CO8Dwf@iw{Ab?M}#$-XUH#Ce-_x+f$Teu;lR;fM)PM`a|j?8TKfDdoOl4^Gg1eTGD*np{-<r;bJl_+?RefWcTDZoc{`;rA)Ct=IL>^v(Byla0U<5R6r`mCji)ybDTYSWvqVAI{mYD+_==PxW$vV`p0H_Ry;`__6T-JEbrZXYYK795T^>@-Ss>ExpISGo`O+=Wg)!Xqevo2zo4Gw^_7iAszeXf;)o&fV0#a9nCt)&Fgt$NhOc1rc{vNDlLy+E5khzfECzfccA?@xrFjj8YStT+F)m}IS66T%d3sL(ECGv5hv;&TH;X4!nqN^q4;uD4pxzIXB;;Da2p!>+^C3?_kmp&NVDYe1MnBIMeee2~hNQ7Btdj9^&j#F){+y@3x=5r6c~CTTdNK5uHY}|`dP~{(CipRc=p_oi0z02Iu`7oq&K&eJz=bPY0!!3Bt?sD*vl*Aasf1cei*^YZ%p@T7&8^IR`hKT^Ae2IM*{n}d-kdzJep=LN6}DFT`A)3>bkprW;@miAZ?9s~+CR}KqZtK@&Ep$3$0nmgfV7PDWmQ$7)fBq_7eB?l6FE>Pp?iQ#+{<=jR}qF;uUe=6?Y{ESzaeKNb;QXI_@?Q0)=OCLYWIsiLtPy^wf7V#r=LE2$>S^sqv*-tY)pt*-RplE<*%s$>Dg}N>Gze9`>zx|TM2v)@^&`KzHUmKiGJ1|gMd#0D#s%z!I*wc@nKt|P_ZwOVm1c{LeI{eKTbY%jSRNd=%gt6<S5=+Rs7XsVH-_ki--jOxKR0L0z%*~xf=SJDol}PSnMxdk0jqo!qf_aoNPwD&cJOtG>uaEP%K`Ej3xr1d=QF#hVg-rkF`3Ch@+ZRCIYdp2%fI~zReX~ZrDPIms^k)a9TtGbh(~iEI7V5i^0{0+SO8zGk-pVb-39i?>ExfAR*=UkML|jYb<xyt1eN_ET_Y`co;4VkGZ3MVEm~RQM9hRQ%AYwXCCBlzoC#+wwc(K17kX{Uv&$7N&aY(ZzBi3135)ZYN1jybo~Q|s>%&ck`ri}RMdXpi8IPQFAsc>N^djr6k&p|loKK)^<rDoFz(bbb`qP~frD|sap6xv-)BwIJD>XW{x~`9W_;YPq62a0P1G4?x<UxIYuYcWq9FvhVVYmS!q1pV3E%2U5Fo#myOiv9gW!14eKv6muQ&ANn(g(Ah5e`s9T7EwZ9mqf<(xka*uRt&Om^_6cL6;fJS~A!T^6bVYS>Ag`d6ekW7K-PK=3sbV+x}y36<%CKA$q<NjD8$L6;uV4LZRTcH00SS@ifa?S5BuViy9};|#6OQ1bW(jOFV|fpOx?Y$}>PVhnfP(Y9Lx<wk~8ic>l4jYI$GR$)!PO4^dLRk}u855-kHDxk)=%O%qRwVc_Ok(XmQhSit;D>rria2HUh2K|Q%`%%i#jI^*vGw`d1q+$P1KcP4}&l3on^YGjGagr^F`pn@-pbh2ll&;li+`_eAB!)0>q`tXZ=f4UZcgP1_bIw|>@f`s2saMp_G0nnJL(lD5RfF%oe_a7OUlm2I&cm^neS9w#r_rn4vK{jqRnlm59kr^%Oxp9^bx#(CnW<*-*x!7A?H&th^5`medG{tTV!}~76oI`;l`$OVKzq!vC1OFLtMgwpLb!1XwRyZyXgrErLZdP0Nk4Gut#CG!ST-EUrIjzDBZ4wzIq6%{+LnH!^{j|}u1CG7L8O+1*Qm*Yl56*siU|;bhtrkUH_6J#6zhzKAgx1%df?~Tt+O7xIg(ilTxekG!Zht!2wW?QG@%xxO-MGZw__n^+|36)r6CYlA0;*lqdvGb%qU82(b~|dqnF0H5pwgB(P(Yhz-63!;VX5pu&6C-wnxhrVyX^d4q~=44PKsieJ#tJ&hQ!m3JdIwC)E6<a3*%sdXf|jbB#8?FpCrfLv)|tmzMRRt{I%Bw!e|jHE&c52Jwrmr#76IaW#{n`$=PxrNzsqHF_l?|9Am~Nu(Wo;cLB$anlfgkt(bM<j2obQ`OJ@y03A;q%Pf(FwI&zfv7}v*zB1prj9C%$7*qwr-q+j68Q#SOi@Xv`$exHV<B94Ok!cUD;x*sI8Qr>H49ty*}zt@j`9y7jH4}UdTO!X-8(9?p<M5)N+N<Zz&{RRphz~Lr`xzrCI`TK+W42)2>xIeRN|B{O=flYu;!p!Vtr!C2>RUcd(Yj-h5COLavR0uj$FMo^Ykhh$`C_Ms(d2raNTWi)dj=ifU}YT)x1Fq&TboNF~z}JzKqq?j=8A$QoO1cpee<RU0IWs`4b(amU7bp-5M?wNtu=ozJIPi)9@=ba!aFwc3*35Yh}Psy$b3!%GTbDWK8ujKyWjp#TVa)!II!W{bzUI`&915R|pAEQuUgNX#L2OEO{!)^<7`ncZ~`{U*3Ok8U#`<HY!_Lh$RC~;$Mu~iP2))vZAY}i0iVMiet`I8Z<{?^)eoN$496kCs)BgN-N7jOyw)YDLG*eb@GD24nyPqCX+77z)8N;99(<F727~2Szf^Yy>T`$Zur?w<c~KLpoWomS=i{q?*m=7O7ZqL7WzPu)<U%oNQ93LnzMFNPUe$xl{j&LhvQdV*m#{5fSj#1OrBU{IOs`ip+)U!pjvFY&8hFK^b@#Cvd0xE8MF$h@(vHylpG9sWe6X@wrmJ@SeF%A{|dl7aX7p@9m|_sDhEz$GhMoXDXIT;#ayk#ZjWVInZ;mpr_cWF++pQy-$=f$krZT6k!2fJAOTz<a^UYQ68E|8R{D0kD4^vhhkFbq8N7-aYg2D`ffa3+Cfb$cYJlv7jH)cqr}15aIz+YP`5UpcaF1>~5>|hRMsjeO@Fh%zBQA4ono=#40i)tcq5}-cbPYDCzfutz_}O!<wjIDekrO7>QU8$XX>y{)<kW{-rF`&ymnvF;A}>g{u}ECh)4xIaooLB6_x%rxTgM*12QEj$HjJ24C4tTDzUQq1<Q1QE&B#PEO^tqoe6}2)FGI&Mst5SSsacGB6yE?8f9uMy;@45@T>C-l(ovy;^NWhX+F8GTI>_XK%8d@at%Us989eiLK`ATVP-@Tf*=5FVKd6AaKkA~LX=#Q8oYWwXeao-Jm4>=HS6UbF#(iD-yVw>Zm7OfAbv8aQcqQ57Y=VC5*eV+I)`?;7<(ahD!K~Ck9{Pm?@nw+XQCV-#y9inZ{K4)&DQ%^@7WtLkzoB|9-1-C)AG*m>ZxY_qpq~jd4JsjwY0z?JHeUbm7^)A33i~KW3)8O96NzAmAg}irqXg|KCY{8ii4PZs2^92<{efqccU9l=6R@q1r*kBZ*GO|=DaV@*Um*mNVEU<65XUE<vHS_Up#*$U6!Z@dWS&m$I@gyuX&QExN&6iiQxIuwvq5XM7!Xv?cV4RWbsLqaEt_N?+e>@`xr<yLlZ*c4GnzAI#z%CBnYej8BO1cUM~pI5604cDE$qdO;Bc@cR3(7-KO*rkE(|r#vFqV=)M_SuFt0jWQXF-2b2H*Ti=-j$hE@P{hio^vn|mCq3*?s6cg9?_qESW~Vu760f}gWDvaSNYpl!*m?t;-6zv$lfLoe}CyPHSh&Bu&k>qWEL?H%{{n_#LN$FBb@BJX-3)}wQm()W<~F&mE4A}~%Wbj07?USld~{5L|iSXu`N^iv$2eN%2tki<m)g>!GSz9*&w5=mEjx0J05YY+TyIMfRUm4(WSkvz~vFWnHf<j_}ccsQ=B>gxkRfR|AJJO8^;ulZ`#!0fZ1sTXdx$bKwAKkl8RiiiKcu)kY8+S?pSdeqy#*w4xHn=~>&^C;+nx2o}YBKi#{Iq`0bm7(?0L4o+do^DAB7y*Ow%59HU6e*f#;=$nM{I@V<PXKdJ!<Y&)aN-kYH5g1L97l!1U{rc9qH+_ZtCkD`@)p)A0Y=GVV?I0Ei&Ka+zden*IX98L+b1g2?i-Pj7ixeTh(|9~E^PK0++L@VN8h`upBhE<$1_)yYQtKCO~LhMu7l>My*(OY6hz*U1GiTxoZc~F!<lwv<F!>bhe01e9(#fjwe?+0#rk6prrts5oX)woPW(ChAKF2=gw)<qIUtxwb<9O1Eh@1~^umcA;gfG-6a1ke&7b7yF?q|>b|-*7)@+nTy_jyP;>t43Y~2069hSgtPn{5rPt(-<kJ9FzXM4HOjnqX4q@1j;-^sY<FF7{CD)%roRR$)2TD_18E;~N-bX@tFs}QRs6SsthMHQ0W*T?(J;*>-fBdFl_cXH=m&TjZH&%+Ecnp_JVl7B#l-$WZ={dWyBehP7XyVJVP4DMl?5b&XN4|GLf4d@_|)XAvqvF^8O6&HxKvKPa3?JW!RQfJJxkpkXe+=4{J*gp0^%lpNRBp`X%`cz|Wb+!F+6IROl84D|4JrH;OJv}Y@mRItzc<ou+nheR|iLXFl7*Tdq4ALXeIi3^n&uRC*N`(gEaF*Y%u<xq;qFejT8&x`j7>e>GEt!eAy^u|1!68|$5Ll_NxRB)r?Qx;Y(_a<Dnj8HpmDU2;$6kR>r6uEMaJYARM3`oWdj_m30;k?Iisjm6fHF3gn@H<Z2e#VmuGH~!KVdK(k$Ab>&_!6sG@rUs<XV{jg{Cdg+VDONBUeDSe}sVa3TSb6<l-tgZ}0<nTEK6pEv)O`>bJZFMdcB$rJaFrwN|0xkk+4Y{hfUg)yBFKR`Q2<+L5mXD^y_OP;3^fNbzMe6bQ{J^58Ld^}=pP<J}8l+rmU;>zA{DTy?7h&E~p;7eXhw8tBLKNddeXTX#b;_K{|qg&Pc0k6I%TBmu2aR};__OK_F_2~pW;guB)BmsTD~jDQVINs9h2;xvHYGy77oo(ZUInlO0p14M$u$14vtfSIze4XdmAn(BZnnr&;5R9&63b}_F{TG(9w*@1uQJY14&<9m=-$i{F`XFkwR*~f$9Vb%RdB&21?gn@BRn{3k(=<-<Ca2U!Sr?X~#Ih{0%%XC{w&hnwWgoF19Msm;p1NqyJFsj)pOffjdO!(y*F8W&CFO4kspDVO+@7>h<kgK#>w2i0eCwd<B%AsP*6I|tMHpv&NEk0<=u=sZZ$!N;%gjCTjs?R583-x4T*C0MXNklv*3$zL!YXGersA8B_NDN&)P3GYsFs2q-4h_pW*?fbY`M*ZyWq}4=eJ{Zge|ZP~aBMOdN9Njam2nj$BY-NV_*-JF-jLUe{1stE4@%7r=7gN`R(rZ~YPXrl%t)AWHA4GQzueIi&_q|Zy&{i;oblkI_26AxpBGblkoAf^<?m^Lx=7LN_<Qdgm8%hz>dCvHdz<RiMyzbB%sx}tF;0(KEY47502bNjCM>D-=FPe`=z)5WW(B3!6#m+hehC6cS-)yP#7)`LhuP_}7x#qw>qSIE{)co(jSPLo^lF}0fBAeMh|XU~yLYaygiydI2Mr-I9LNM#zt62fotI<i2ffgZPFSj8$<GJ5deQ^SD<?C+%_jv-&!iFBsNAU}=zt!#L`0O9MNY}8EC~2}(AIi4-uU5bY*$ZRyXJXGv7Me$O_|v~G$N6vLk($_>_yfrkn~(YAQPDO+2hapv*;+lJpl1JA|A!9sqpYVdn7|LL6T>_$!b5vp@x7`$S^|OwDF9P&c1@zL|H}<b3rlSYpU`-fuTREt-sg`kF8?e$i65n#xGn8E}^_n9!|NUJ+nO_h5>M$MwT}ca!%i1d;D$lBFME!9)}_W0>hW`9<uR1u`<Yi&eeLyoNR_wJ+t-m#(7J<OKp5UCNtxN1Kps}@mF)IebZlU`L1j!?lAGn`5z#t4KSxIsx(p~<gx~50R?&&vhfmRpS)oUK?VuW<AKY<4bZ2M_+^{Ekos4UHOL=-y`jLux_ok52<-)W&oXL?`h0|^KfbNx6Q#6`o<o^JH%nfPf;0xMG!j+DFfDb~DK*MSJOnU};$MF!WW?umg3^oXX5+JpIR{Ij($qtmfcvqu5K<Af|B;jlB!}jc&}2+kmmvu4XzG}k{1|cOZ@zzr@O2ecdM7)x@ZWwn9@8qjkk1&8{<M%-K!4(;-A?z=X<vzUI7fH_FEhMILjdU>J(z623HPLd$|v1se#f1K_<rSOt8b-#;opBQNu9ZJl1IabM+P$Eqw{i&06t6pUIv@Cz;)l&7P-yLHQ$H)9{3DU6?A_zZ()JNM=0FEI9^}3N1`Mm=vvSY!A_<Q+eU?1xRL94cl95M2$<k-TjzcO+Kq_e8(5bfVBab8GnjCxK2|o8GA&aFRhaN~h^*10YE@$wKRiUeIjqMKwvhqK2b}>HQ!287;0|d6K+cyc|D&^6y6L=CVd8*-z`w&h7dW6vLGrzX=_MY-$zJO1p+5zpKCQ7>g*dCwuclgE_*3Xy@zE3lX%jUoBYohmlKS`8s<Ch!sm8mYzCYHuN9U-6ed*HYw$Z10c0|-wmjCO3PwX#u&rZyY-=&NzO6b~vh>f)&Hi{#cTZU{I$4A39*m>gy&3ak=je@Pks<nS%?Sgst^e=gvQw0n445PtdDIs`O?U&D}7p;hCG^q53$Bj$<?-tf!`ieZfzK##7zU3)yt0hl%%On$=gj926e&B)!bI7c~*6XEob85If?kOZB@C=;m_q??zKf6$HKF>e(nsJ%DWW8KPO9L@vzsUJWfX)N1Cxw_z!wyU3Sh?<AoT$ilGcciAtLI}YGc>Ki`s`TW$9wTYfaFe`e`S?HixyprVSMC?JrJwK7~qp+yC)d7iAaGM3`P$w!nx;DF4bT3{$T$}7nMRXvIX2jos|GOA(31p7(6AybZ3BB4~F&7+-8=?Gm_V^xCo#&sd4Yo@I?mPwT6mQ0K)jI;H(<aprYru?jA|3M5Q)#1tVOSIhWe|`jbzBcGjhMQNGcoFTwo=802;>oETT&)#~=1+!p_PGaO|b3lsT?40?vAk^nhk%~Z(BB2=F&%9&vJR~5Z5or-Vx+2vr}v4Yh?i6j`g2GfDD3c@h3z-7dmBQ%yDMZi4$KJSXG>Y#@UM~eb{t^;KdYB56&yzeXcrVxbO0b4QI9}W~c9CtoX9`(@D4Fq%Y37;@wq=U^J*nhg^LAm+n3!KBdC3e6$T`oFX3Gg4Y9Ydb|zImIV4i^+>NnDM1@l6P21iA9gZ)gmzIKT0hY0L~G-bh)RA}ORpm|UZHS#-C9&neu`pReLj&1Hz}+cGemuQk%BZ`TgD)Wxd2ICGij#@n^bf)c_2z1zBf(|PpWIJ^uC)FJnS?Ga|F2wKch3$R+>OjW>-c48y2QxZJIE5G2cJ9;yh4-OL7Tl+EeGsvcG35&XjbO(^|%)WmN`PQQRhxc5n1S@z9VrK;+31$0s)O4S0&znz|Ut}bu(u1A3JHDNvaU$z!J>p)Yf$YgGMnC4!#+V!IP>8^58j6<gBY+i7*5Qge&le5q4D$&2z0|cB6Ly3El)Sfs<~yBXqcX8nQ2mVBO}P}iA9xg^Kt(u{B(CRY7)5z?wNI@Q_Dh{Tgo#n9S(iFKaM>CLZh8MI+Cb0XYuV`Et@lcu;D@wDls?I56epya6D#@CX!1#i6pLjr9dpW5Z`kIr&jN;PBOrN>0I4cuXX?;e9}%hrhxUzsH6?iq{-;UPq#AYnV<)}m<{bC>B}pGLe?@iUUY5e1H$HAYmnR{dCA3?#%mE3b;y5|J$ssYmo3B|<i@fB9Y3oxbH4%hoY+knOXG3bGu(O7zszjiTun6|c=p2<6H*tbJ@T52C521g?emV09a;g~UE?H_a?DjZN9Q^a(9vuilZu3C1@QWWYzUu0ze3kn5vwNu9tIDN!%ZXOL>H~^%`Z@dLAG<JzZNi5toaN=wk{&k}W8xwbCd35MjL?v5F-*o><7<249;^k7zEvvGpYeT8O6fWlD9oQGZLna?v9F+yJr;>_SFRJCduU$Bz7Qb-w~7S_1!v}qskKl&;oi8$iR0|gxi7^?13E?wC2{A7r&TLWR6WpjJ^60^B76%oaSO~2TQ2QrPH30~801U|_|zmAJ4f8VMOFC0KqxT!EBbofm$@5efC={0yNCJ{?~n_buro(i0Ly>uq~A&idCmR@p3F>kr-I&=s%-wpby0OPGvYOTLAML{>AxH0=??X)Ca|w|<DRuHY9gzS-tVu;P;94fzIX@9sJ=!JDh=2$6A2JDz$S3jHq8a;(5__h3=}yV$YNTJfbNfOF2{}6X#bYGnbCJ+K{mOIZieF}6S{uQsrAynWzZ?fu$k?*j<JZ&+OWm@!}?RX^^=202-E_>{<^*X)HdGFwQ{J3INC)ca&TeA&gu2$&%+}dk~jRoDb(TRV;yCC0!4qyRs6y<p^D|Fmx>hQFuN4xR{<hF^B>A*bM30y4+*w+lM>v_7NwoM@xQ2;+s4#UCHkv~n8+fgkO+Ma)jo@!xrnpzynGrnEk{QNBlg))PWpKgVZ#G#G?(acbh?KR!<7aa)@7<>Y$zp7x<;kqb^&Ou32^1GLFE+B->-X$;qQjRG>8ebc{mZZ&}6H6DZ<mKO=`;tVe`@K8|#84&cL+XxcDHzPE~rj$HpuR$Lr(cQGO^WIAb^W`mh4(_a#*%<T-wX;;#Ju$o7vWo^7k<sY|RL^)FT@;>N1yzD`{oOo<*}f{-;eeHh*dfr3z?^=zWil2sywI5bI~3GF6DM?86REQ48&iB&U78>qykz?@YG6B188ruZ_;pa9?WvkA5Z>xOhbsr~DnKxcA-AUAxA%8gFW3H+@vyAmiF4Xi`q2=`<d&waB@cd}bOgrnW>Ot#+JSRTFs1o9lu_h(kJDMW$Z15*C0>2G4Bf_G$6u)<Rq!Xz?ac8of~2o8Zn@EN*No_Q%5s)z4+-AE#`n4osnsY>3G1o~v-$&nqGfPEL*%vn4Av7em~KP0oXks!R1@!)`S#-Gb}G&oxO%3PxTCuxl8T7BU)6SIjD^uDAiwMhtJ1V?*6uK9gr<9EB?GV~c#Sp4AAxD(NlO?cU(R(S(MidLW>ThOK|6sMEvV}RLl5AH1W4|0>Q8D7k_K3W^UYx{Md)Gt!p<<UOB$EG8wkc1*bP3b89Sl$;x500}5W@=xrCk%~!mq36CH%*6#4VRYQ74&WvhF40fzYJ&=aPug85+s^x)=^@2`s(;(#Lb1D7m-l*@QCg6+F4W~^cv6e@Sz~Asg5Ha^Zn_)ekmIZEFBsn3o}qCMB!E7s6VQL9IeDoy^&0_baWDCK9i_Z=kofRJMcvEE|(GUzAlNy&!MSuqrx0e-D@T^tyhm9p%MGOz;mgEs@J?+P&+@@5oNYXb)szS@;IdJENTB!Fk<7Q!^Us#x;s-khR<ch3o#t3jh^|XdfwGWk@=p@1Vy^T3xfOnSXM<E8)1qKj(vdtf)cgmY#rs5(V4U<(~<pEOlNEB$IQMOOW`D;3A>JbT7x+cpCG%g`E}&!K<<w}_c6Gj?Kmmme_KCuqKZ0)fYWI~Obz#_<m)kGbBw9ycJK9wVK7ADilN0j1;_$_t5sNFdy$#8k0z0XgTFYszyIf#97SIDMKf0At6|!gA76`Lt<T7kQJWtBBoXwz^vJRR74NCpcB8{q3Z;B7H-VgZigZavec@qqv1%iC?E8c1`H{Jf=L7Oggr;`zbqt?E#s8&<u*Zj|z87&XK$~VC6CJB~J3``7f~7!`z;0-|qA6mvREx#>y|o9%;>!3*l0}5`O54J_%3}A!<Eh-QM8{nimPbI+l)$VfHKF?1Lk+tzMVn>JXIw#<w888G9h@`_CJJFsXKmeY-CUwLN2<z|HcVe~BC&i#1oaH22vzhff3I9&8$rFPyJ1=lCe|yyqfSMyAZoW+UrLo0YGvu+aOEnet-)i~)*XX|l2(L3`^g%NkP_D&#ylYyHxwCtpFdh>Mb?51c2sEPSlTVaAyOg_Zbjq>(LZu;;2=D`I}WU-KK+)T3jxZpf~Z8%C$iwcVKx*>w^i!+9HK?z;?ilI*~B%F+Ul;usg04%et^PzyD-0*;m;+eAs%-e-IpcDSsPq#L?SFId5baKRrz0@)mIeKj<D0SC)I_{tyfm|Sx(5Z?v}Rl>XEDteQTFYZ<TJg+}5)HJHO#U@VqftcP(*`1w<u+jsceQei26%WF48{X+X#jw@$l06+J;`u=6_jM!a0Rp1_Jx?n*aO0M7xpAQ%QpLuKyook5wA9cVm*W#zn7E)q#WH*Y<jm;Kj?JmLRfn1yo+G3F=Vag&b*o^Rsz2F{qaG<t!+*cA>UP>P27GA(UY<N_=6v3Rw<pX@NsXW`NZWYUK?j6$wj+SL4IcBxY319?%zPP0^asi<OQIr?Qi-x{_A6wHYv@nk>v%RbjY>sh#*-Cs5L*R{ZjC8)mgQtS7{Z1V{6hN4l%O?E8e!MJbEj(Ooj(cpV#i<n)Tvo<00JepiYYvq)(nf!2U-j%o>PCpf{BEuA0Q_!e;Gd$dCX;FrU!EhJWUL1bb!&2;Ei6f#)S?2d8y<!*a=iqO4<!)2MwbIuUO}UkAijM7_;0oKH4!}gnW?Mlf!@jMg35q9SfPz>Qw1$h`87%3hI*NCe8Cq%R^nL!Zpng{KS2?A%$mE3GD~k-aY@WUwS_RT-A>}fDtdiBaJ|5Yxv{B6xNmtG|QWCg$5pxBV8l5Ws<R^%VlDv*fNGKDBVIvVIeZS<1Qd(m;{PN?}X%Ef1`OP6=Hq0aM34-vNS_I_w*diOlLH`0SflVy0AMkPGD9<+=nXSMWm9mY`S<=|tr$3peS0izMFS(nT=rlr<s87^Oh;3iGwlA)X0(J&wG^?f(d{%yy2y53D`?dIlVReRmN<A{KfJnI&E?y``=tHhlZgq|Dt=AKfR~?J3'
)
PREFIX_RAW_SHA256 = '4634a38b66f6dcc51092a4a351eb3f3c863d9ed5d7d6bdf9e93e5e5a1e48f842'

PREFIX_DEPTH = 30
DEPTH = 61
def transfer():
    return Transfer(225, 225, (117, 270, 4), (126, 286, 4), (144, 356, 4), 580, 710, -14100)
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
    if energy['THRESHOLD'] != gain['T'] or gain['T'] != F(-14100, 40000):
        raise ValueError('threshold agreement')
    rows, initial_energy = moments['compute']()
    cost = energy['energy'](rows, initial_energy, moments, gain)
    initial = energy['initial'](rows, moments, gain)
    if F(cost['energy_upper_1e12'], 10**12) >= F(123153850, 1000000000):
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
    if F(entropy['coefficient']) != F(40816655382787, 8017315635000) or entropy['checks'] != 72:
        raise ValueError('joint score or dual accounting')
    return {
        'threshold': '-141/400', 'graph_cap': '123153850/1000000000', 'coefficient': entropy['coefficient'],
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
