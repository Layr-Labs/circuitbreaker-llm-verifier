"""Directed Gaussian recurrence for the norm-one circuit certificate."""
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
from fractions import Fraction
from math import factorial,isqrt
def endpoints(x):
    if abs(x)>Fraction(1,2):raise ValueError('small exponential')
    s=sum(x**k/factorial(k) for k in range(5))
    # Lagrange error <=2*|x|^5/5! on this interval.
    e=abs(x)**5/60
    return s-e,s+e
def gaussian_table(q,pi_low,pi_high,S,H):
    L=1<<116
    def floor(x):return x.numerator//x.denominator
    def ceil(x):return -floor(-x)
    x=Fraction(1,q*H*H)
    lo,hi=endpoints(-2*x);rl=floor(L*lo);rh=ceil(L*hi)
    lo,hi=endpoints(-x/4);ll=(L*L/(q*pi_high));rr=(L*L/(q*pi_low))
    fl=floor(isqrt(floor(ll))*lo);fh=ceil((1+isqrt(floor(rr)))*hi)
    ml=rl;mh=rh;dens=[];cum=[(S//2,S//2)]
    for _ in range(8*H):
        f=fl*S//L;g=(fh*S+L-1)//L;dens.append((f,g))
        a,b=cum[-1];cum.append((a+f//H,b+(g+H-1)//H))
        fl=fl*ml//L;fh=(fh*mh+L-1)//L
        ml=ml*rl//L;mh=(mh*rh+L-1)//L
    return tuple(cum),tuple(dens)
def lower_weights(d,mode,N,K,scale,substeps,normal=1<<40):
    F=Fraction;a=F(28,57);t=F(41,57);h=F(1,d);L=1<<116
    def lower(x):
        y=endpoints(x)[0]*L
        return y.numerator//y.denominator
    piU=16*sum((-1)**i*F(1,5)**(2*i+1)/(2*i+1) for i in range(15))-4*sum((-1)**i*F(1,239)**(2*i+1)/(2*i+1) for i in range(14))
    fpi=F(564189583547756,10**15)
    if fpi*fpi*piU>=1:raise ArithmeticError('pi lower normalisation')
    def half_grid(c):
        step=lower(-2*c);ratio=step;base=lower(-c/4);out=[]
        for k in range(N+1):
            out.append(scale*base//L);base=base*ratio//L;ratio=ratio*step//L
        return [out[-i if i<=0 else i-1] for i in range(-N,N+1)]
    v=h/(a if mode=='A' else 2*a*a*t)
    cf=h*h/(2*a*a if mode=='A' else 2*a*a*t*t)
    cg=cf if mode=='A' else h*h/(4*a**4)
    f=half_grid(cf);g=f if cf==cg else half_grid(cg)
    c=(substeps*h)**2/(a*a if mode=='A' else 2*a*a)
    ratio=lower(c);step=lower(2*c);base=lower(-v*v/12)
    fac=normal*v*fpi;num=fac.numerator;den=fac.denominator*L
    out=[]
    for k in range(max(map(abs,K))+1):
        out.append(num*base//den);base=base*ratio//L;ratio=ratio*step//L
    return f,g,[out[abs(i)] for i in K]

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
        self.cum,self.dens=gaussian_table(q,PI[0],PI[1],S,H)
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
 layout=Struct(f'<{size}Q')
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
def weights(d,mode,N,K,substeps):return lower_weights(d,mode,N,K,WSCALE,substeps)
TAIL=tuple(range(Q+1))
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
  indices=tuple(range(1+max(self.msize,self.asize)))
  def ix(ip,jval,axis):
   if ip==len(axis) or jval>y_hi:return indices[len(axis)*len(Y)]
   jp=y_lookup[jval-y_lo] if jval>=y_lo else 0
   return indices[jp*len(axis)+ip]
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
   @lru_cache(None)
   def active_weights(start):
    weights=g[start:]
    fmt=Struct(f'<{len(weights)}Q')
    on=pack([0]*start+[v*Q for v in weights])
    off=full_g-on
    return weights,fmt,bytes(8*start),on,off
   full_g=pack(v*Q for v in g)
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
     g_active,active_layout,active_pad,on,off=active_weights(active_at)
     active_columns=pd[active_at:]
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
  fine=getattr(self,"fine",False) and not isaux
  bits=10 if fine else 0; top=Q<<bits; pool=range(top+1) if fine else TAIL
  KK=2*LS+40+LQ-bits
  other_m=n+[0] if isaux else None;other_p=u+[Q] if isaux else None
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
    if n is m and u is p:
     ac=prod(fm,2*(full-fp)-fm);bc=prod(fp,2*(full-fm)-fp)
    else:
     fn=pack(map(mul,fw,ff(other_m)));fu=pack(map(mul,fw,ff(other_p)))
     ac=merge(fm,full-fn-fu,full-fp,fn)
     bc=merge(fp,full-fn-fu,full-fm,fu)
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
   outm.extend([pool[(x*r)>>KK] for x,r in zip(ma,norm)])
   outp.extend([pool[top-((total-y*r)>>KK)] for y,r,total in zip(po,norm,total_weighted)])
   if not isaux:
    startrow=len(outm)-cnt
    for at,fi,gi,mass in corrections:
     fm,fp=fneg[fi],fpos[fi]; gm,gp0=gn[gi],gp[gi]; gzero=Q-gm-gp0
     dx=(Q-fp)*(Q-gm)-fm*gzero
     dy=fp*gzero+(Q-fm)*gp0
     outm[startrow+at]+=(dx*mass)>>(26+LQ-bits)
     outp[startrow+at]-=(dy*mass)>>(26+LQ-bits)
  need(all(0<=x<=top-y and y<=top for x,y in zip(outm,outp)),'tail simplex')
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
    cap = F(122792154793,10**12)
    coefficient = 1 + 1/(2*cap) + delta
    return {
        'checks':len(checked), 'semantic_increment':str(delta),
        'weighted_margin':str(weighted_margin), 'private_path_margin':str(cutoff-delta),
        'coefficient':str(coefficient),
    }


from array import array
import lzma
def plane(encoded,width,n):
    from base64 import b85decode
    data=lzma.decompress(b85decode(encoded));v=array('H');j=shift=0
    for byte in data:
        j|=(byte&127)<<shift
        if byte>=128:
            shift+=7
            if shift>14:raise ValueError('invalid continuation')
            continue
        i=len(v);x=j//2 if j%2==0 else -(j+1)//2
        if i%width:x+=v[i-1]
        if i%n>=width:
            x+=v[i-width]
            if i%width:x-=v[i-width-1]
        v.append(x&65535);shift=j=0
    if shift or len(v)!=2*n:raise ValueError('plane size')
    return v

"""Decode an untrusted bilinear table, then apply one directed smoothing pass."""
from array import array
from pathlib import Path
import lzma,sys
def knot_target(values,s,stride):
    w=len(s.main);h=len(s.y_axis)
    xs=list(range(0,w-1,stride))+[w-1];ys=list(range(0,h-1,stride))+[h-1]
    size=len(xs)*len(ys)
    if sys.byteorder!='little':raise ValueError('endian')
    if len(values)!=2*size or any(m>16384-p for m,p in zip(values[:size],values[size:])):raise ValueError('input knots')
    def plane(table,upper):
        out=[]
        for j in range(h):
            j0=min(j//stride,len(ys)-2);y0,y1=ys[j0:j0+2]
            left=j0*len(xs);right=left+len(xs)
            v=[(y1-j)*a+(j-y0)*b for a,b in zip(table[left:right],table[right:right+len(xs)])]
            for i in range(len(xs)-1):
                x0,x1=xs[i:i+2];den=(x1-x0)*(y1-y0)
                out.extend(((x1-k)*v[i]+(k-x0)*v[i+1]+(den-1 if upper else 0))//den for k in range(x0,x1))
            out.append((v[-1]+(y1-y0-1 if upper else 0))//(y1-y0))
        return out
    m,p=plane(values[:size],False),plane(values[size:],True)
    return s.apply(s.M,m,p,*s.apply(s.A,m,p,m,p))



def waypoint(k,s,knots_by_id):
    s.fine=False
    count=len(range(0,len(s.main)-1,4))+1
    height=len(range(0,len(s.y_axis)-1,4))+1
    knots=plane(knots_by_id[k],count,count*height)
    X=knot_target(knots,s,4)
    return X


def cubic(a,b,c,d,n=3):
 if min(a,b,c,d)>=0 and a>0:return True
 if n==0:return False
 m=a+3*b+3*c+d
 return cubic(8*a,4*(a+b),2*(a+2*b+c),m,n-1) and cubic(m,2*(b+2*c+d),4*(c+d),8*d,n-1)

def lift(old,new,m,p):
    ix=[bisect_left(old.main,-((-i*old.d)//new.d)) for i in new.main]
    a=[];b=[];roundup=ceil;w=len(old.main)
    den=41*new.yd*new.d*old.d
    for j in new.y_axis:
        for i,k in zip(new.main,ix):
            if k==w:
                a.append(0);b.append(Q);continue
            v=old.main[k]
            num=old.yd*(41*j*new.d*old.d+13*new.yd*(i*old.d-v*new.d))
            row=bisect_left(old.y_axis,roundup(num,den))
            if row==len(old.y_axis):a.append(0);b.append(Q)
            else:
                z=row*w+k
                a.append(m[z]);b.append(p[z])
    return a,b


"""Share the packed operands of all six parent products in each Gaussian row."""
from operator import mul
from array import array
def check_segment(r,s,A,B,aa,ab,bb,cubic,prior=None,save=False):
    q=r['Q'];KK=2*r['LS']+40+r['LQ']-10;top=q<<10
    _,rows,pack,prod,merge,full,norm,total,parts,w,fw,_=s.M
    coeff=r['coefficients']
    f=[(m+[0],p+[q]) for m,p in (A,B)]
    g=[(m+[0],p+[q]) for m,p in (aa,ab,bb)]
    failures=0
    last=(array("I"),array("I")) if save else None
    for rowid,(ff,active,on,off,corrections) in enumerate(rows):
        fv=[]
        for m,p in f:
            fm=pack(map(mul,fw,ff(m)));fp=pack(map(mul,fw,ff(p)))
            fv.append((fm,fp,full-fm,full-fp))
        gv=[]
        for m,p in g:
            gm=active(m);gp=active(p)
            gv.append((on-gm-gp,gm+off,gp))
        values=[]
        for fi,gi in ((0,0),(1,0),(0,1),(0,2),(1,1),(1,2)):
            if fi==gi==0 and prior is not None:
                start=rowid*w;values.append((prior[0][start:start+w],prior[1][start:start+w]))
                continue
            fm,fp,fn,fu=fv[fi];gz,gm,gp=gv[gi]
            neg=coeff(merge(fm,gz,fu,gm),parts)
            pos=coeff(merge(fp,gz,fn,gp),parts)
            nm=[(x*z)>>KK for x,z in zip(neg,norm)]
            pp=[top-((t-y*z)>>KK) for y,z,t in zip(pos,norm,total)]
            m,p=f[fi];other_m,other_p=g[gi]
            for at,i,j,mass in corrections:
                x,y=m[i],p[i];u,v=other_m[j],other_p[j];zero=q-u-v
                nm[at]+=(((q-y)*(q-u)-x*zero)*mass)>>(26+r['LQ']-10)
                pp[at]-=((y*zero+(q-x)*v)*mass)>>(26+r['LQ']-10)
            if any(x<0 or y<0 or x>top-y for x,y in zip(nm,pp)):raise ValueError('invalid parent tails')
            values.append((nm,pp))
            if fi==1 and gi==2 and save:
                last[0].extend(nm);last[1].extend(pp)
        start=rowid*w
        for j,sign in ((0,1),(1,-1)):
            for a,b,c0,c10,c11,c20,c21,c3 in zip(A[j][start:start+w],B[j][start:start+w],*(v[j] for v in values)):
                if (b-a)*sign<=0:continue
                v0=sign*(3*c0-3072*a);v1=sign*(c10+2*c11-1024*(2*a+b))
                v2=sign*(c20+2*c21-1024*(a+2*b));v3=sign*(3*c3-3072*b)
                if not cubic(v0,v1,v2,v3):failures+=1
    return failures,last

"""Integer lookup cache for the R,U quadrature with the actual threshold inserted."""
from bisect import bisect_left
def finish(s,n,u,b,r,step):
    Q=r['Q'];ceil=r['ceil'];cdf=r['Gaussian'](2,26)
    scale=1<<26;g=40000;den=114*g
    if 4*g%step:raise ValueError('panels must divide the support')
    h=4*g//step
    xs=sorted(set(range(-4*g,4*g+1,step))|{s.T})
    def mass(lo,hi):return max(0,cdf(hi,g)[0]-cdf(lo,g)[1])
    wr=[mass(lo,hi) for lo,hi in zip(xs,xs[1:])]
    wu=[2*mass(i*step,(i+1)*step) for i in range(h)]
    mr=scale-sum(wr);mu=scale-sum(wu)
    if min(mr,mu)<0:raise ArithmeticError('root weights')
    w=len(s.aux);limit=len(n)
    tab=[(min(bb,uu,Q-nn),min(uu,Q-nn),Q-nn) for bb,uu,nn in zip(b,u,n)]
    tab.append((Q,Q,Q))
    ymin,ymax=s.y_axis[0],s.y_axis[-1]
    row=[bisect_left(s.y_axis,i)*w for i in range(ymin,ymax+1)]
    # Precompute for each remainder class induced by the inserted T, and R+V.
    offsets={}
    for rem in {R%step for R in xs}:
        cols=[]
        for k in range(-2*h,2*h+2):
            sm=k*step+rem
            ip=bisect_left(s.aux,ceil(28*sm*s.d,den))
            # P=(56R+41(R+V))/den. Keep R separate.
            z=s.aux[ip] if ip<w else 0
            cols.append((ip,s.yd*(41*sm*s.d-den*z)))
        offsets[rem]=cols
    tot=mr*Q;dy=den*s.d
    for R,rw in zip(xs[1:],wr):
        k,rem=divmod(R,step);cols=offsets[rem];origin=k+2*h
        add=s.yd*56*R*s.d;t=[]
        for ip,v in cols[origin-h+1:origin+h+1]:
            j=(v+add+dy-1)//dy
            t.append(tab[limit if j>ymax or ip==w else (0 if j<ymin else row[j-ymin])+ip])
        if len(t)!=2*h:raise ArithmeticError('column range')
        q=mu*Q
        if R>s.T:
            for (A,B,C),(D,E,F),uw in zip(t[h:],reversed(t[:h]),wu):
                q+=uw*((2*(C*E+B*F-B*E)+A*(Q-F)+(Q-C)*D+2*Q-1)//(2*Q))
        else:
            for (A,_,_),(D,_,_),uw in zip(t[h:],reversed(t[:h]),wu):q+=uw*((A+D+1)//2)
        tot+=rw*ceil(q,scale)
    return ceil(tot,scale)


def grid(d):
    def axis(c,t):
        return round(F(c*d,225)),round(F(t*d,225)),4
    return Transfer(d,d,axis(117,270),axis(126,286),axis(144,356),
                    round(F(666*d,225)),round(F(824*d,225)),-14500)
def reach(request):
    import runpy
    from pathlib import Path
    b=Path(__file__).resolve().parent
    if (b/'files'/'gain.py').is_file(): b=b/'files'
    r=globals()
    if request not in ({'phase':'entrance'},{'phase':'terminal'}):
        raise ValueError('expected entrance or terminal phase')
    g=runpy.run_path(str(b/'gain.py'));s=grid(285)
    if g['T']!=F(s.T,40000) or g['T']!=F(-29,80):raise ValueError('threshold')
    def aux(x,y):
        s.fine=False
        return s.apply(s.A,*x,*y)
    if request['phase']=='entrance':
        old=grid(125)
        A=lift(old,s,*old.evolve(14));del old
        for _ in range(2): A=s.apply(s.M,*A,*aux(A,A))
        A=s.advance(*A,5)
        B=waypoint(21,s,g['KNOT_B85'])
        bad,_=check_segment(r,s,A,B,aux(A,A),aux(A,B),aux(B,B),cubic)
        if bad:raise ArithmeticError('entrance comparison failed')
        return {'threshold':'-29/80','certified_waypoint':21,'failed_coordinates':bad}
    e=runpy.run_path(str(b/'energy.py'));m=runpy.run_path(str(b/'moments.py'))
    if e['THRESHOLD']!=g['T']:raise ValueError('energy threshold')
    rows,u=m['compute']()
    cost=e['energy'](rows,u,m,g)['energy_upper_1e12'];initial=e['initial'](rows,m,g)['initial_lower_1e12']
    if cost>122792154792 or initial<540711206696:raise ArithmeticError('directed bounds')
    truth_controls();literal_minimum_controls();sem=semantic()
    if sem['checks']!=72 or F(sem['coefficient'])!=F(62636417857592837,12279215479300000):
        raise ArithmeticError('coefficient')
    A=waypoint(21,s,g['KNOT_B85']);aa=aux(A,A);prior=None
    for k in (48,75):
        B=waypoint(k,s,g['KNOT_B85'])
        B=list(map(max,A[0],B[0])),list(map(min,A[1],B[1]))
        if k==75:
            B[1][21536]=A[1][21536]
            B[0][21536]=min(B[0][21536],Q-B[1][21536])
        ab=aux(A,B)
        if k==75:
            s.fine=False;final_aux=s.apply(s.A,*B,*B,plus=True);bb=final_aux[:2]
        else:bb=aux(B,B)
        bad,prior=check_segment(r,s,A,B,aa,ab,bb,cubic,prior,k!=75)
        if bad:raise ArithmeticError('terminal comparison failed')
        A=B;aa=bb
    upper=finish(s,*final_aux,r,400)
    if upper>8032:raise ArithmeticError('terminal majority')
    return {'source_waypoint':21,'certified_waypoints':[48,75],'energy_upper_1e12':cost,
            'initial_lower_1e12':initial,'terminal_upper':[upper,Q],'semantic_checks':sem['checks'],
            'coefficient':sem['coefficient']}
def run(request):
    return reach(request)
if __name__=='__main__':
    import json
    print(json.dumps(run(json.load(sys.stdin)),sort_keys=True))
