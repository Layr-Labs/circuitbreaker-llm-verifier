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
    def __call__(self,n,d,side=0):
        if d<=0:raise ValueError('CDF positive denominator')
        if n<0:
            if side:return self.p-self(-n,d,-side)
            a,b=self(-n,d);return self.p-b,self.p-a
        if n>=8*d:
            if side:return self.p if side>0 else self.p-1
            return self.p-1,self.p
        k,l=divmod(n*H,d);a,b=self.cum[k];f,g=self.dens[k]
        if side<=0:
            a=max(0,(a+f*l//(d*H)-self.err)//self.dg)
        if side>=0:
            b=min(self.p,up(b+up(g*l,d*H)+self.err,self.dg))
        return b if side>0 else a if side<0 else (a,b)
"""All-positive finite Gaussian propagation; exact nonnegative integer products.
Consecutive identical transitions are coalesced; ±1 states have analytic merges.
"""
from struct import Struct
from operator import itemgetter,mul
from bisect import bisect_left
from itertools import chain
import sys
LQ=14;WSCALE=10752;Q=1<<LQ
A_WEIGHT,T_WEIGHT=F(28,57),F(41,57)
def integral_weights(*xs):
    from math import lcm
    xs=tuple(F(x) for x in xs)
    den=lcm(*(x.denominator for x in xs))
    return (den,*(int(x*den) for x in xs))

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
  high=(value.bit_length()+63)//64
  raw=value.to_bytes(8*high,'little')
  if raw[0]:return 0,high,raw
  low=(len(raw)-len(raw.lstrip(b"\0")))//8
  return low,high-low,raw[8*low:]
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
 a,t=A_WEIGHT,T_WEIGHT;h=F(1,d)
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
  self.core=45;self.middle=80
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
  a,t=A_WEIGHT,T_WEIGHT
  ell,delta=2*a*a*t,t-a
  pinned_threshold=F(T,40000)
  gamma=delta/t
  beta=(a*a-t*t-a*t)/(-a*t)
  aden,u0,u1=integral_weights(F(yd,cd),gamma*F(yd,d))
  sden,s0,s1,s2=integral_weights(F(yd,cd)/t,gamma*F(yd,d),t)
  oden,o0,o1,o2=integral_weights(F(yd,cd)*beta/(2*a),F(yd,d),a)
  y_lo,y_hi=Y[0],Y[-1]
  y_lookup=tuple(bisect_left(Y,j) for j in range(y_lo,y_hi+1))
  ids=tuple(range(max(self.msize,self.asize)+1))
  def ix(ip,jval,axis):
   if ip==len(axis) or jval>y_hi:return ids[len(axis)*len(Y)]
   jp=y_lookup[jval-y_lo] if jval>=y_lo else 0
   return ids[jp*len(axis)+ip]
  self.monotone_main=[]
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
   # Main and auxiliary both consume only odd product coefficients.
   pack_even, half_product, half_merge = packed_arithmetic((tc + 1)//2)
   odd_layout=Struct(f'<{tc//2}Q')
   def pack_odd(values):
    return int.from_bytes(odd_layout.pack(*values),'little')
   half_width=16*((tc+1)//2)
   reduced_parts=runs([j//2 for j in digits])
   even_full=pack_even(c*Q for c in f[::2])
   odd_full=pack_odd(c*Q for c in f[1::2])
   if isaux:
    quarter_size=(tc+3)//4
    _,quarter_product,_=packed_arithmetic(quarter_size)
    quarter_width=16*quarter_size
    quarter_weights=[f[r::4] for r in range(4)]
    def make_pack(r):
     layout=Struct(f'<{len(quarter_weights[r])}Q')
     def pack_lane(values):
      return int.from_bytes(layout.pack(*values),'little')
     return pack_lane
    quarter_pack=[make_pack(r) for r in range(4)]
    quarter_full=[p(w*Q for w in ws) for p,ws in zip(quarter_pack,quarter_weights)]
    quarter_parts=runs([digit//4 for digit in digits])
    residue=(2*N+1)%4
   k_axis=range(-N,N+1)
   if isaux:
    locations=[]
    for k in k_axis:
     ip,i=sample(M,ceil(d*k,cd))
     dy=ceil(u0*k-u1*i,aden)
     locations.append((ip,dy))
   else:
    ray_locations=[]
    aux_locations=[]
    for k in k_axis:
     ip,i=sample(M,ceil(t.denominator*d*k,t.numerator*cd))
     ray_locations.append((ip,s0*k-s1*i))
     zp,z=sample(A,ceil(a.denominator*d*k,2*a.numerator*cd))
     aux_locations.append((zp,o0*k-o1*z))
   for j in Y:
    if isaux:
     od=[ix(ip,j+shift,M) for ip,shift in locations]
     plan.append([itemgetter(*od[r::4]) for r in range(4)])
    else:
     straight_y=s2*j
     aux_y=o2*j
     od=[ix(ip,ceil(shift+straight_y,sden),M) for ip,shift in ray_locations]
     pd=[ix(ip,ceil(shift+aux_y,oden),A) for ip,shift in aux_locations]
     fe,fo=itemgetter(*od[::2]),itemgetter(*od[1::2])
     shift=(pinned_threshold-F(j,yd))/delta
     vcrit=shift*ell
     scaled=vcrit*cd
     active_at=max(0,min(tc,scaled.numerator//scaled.denominator+1+N))
     g_active=g[active_at:]
     def mk_pack(parity):
      start=active_at+((active_at&1)^parity)
      col=pd[start::2]
      gw=g[start::2]
      layout=Struct(f'<{len(col)}Q')
      pad=bytes(8*(start//2))
      if len(col)>1:
       item=itemgetter(*col)
       def active(table):
        return int.from_bytes(pad+layout.pack(*map(mul,gw,item(table))),"little")
       return active
      if col:
       idx=col[0]; w=gw[0]
       def active(table):
        return int.from_bytes(pad+layout.pack(w*table[idx]),"little")
       return active
      return lambda table:0
     ge,go=mk_pack(0),mk_pack(1)
     w=[0]*active_at+g_active
     anti=[v*Q for v in g[:active_at]]+[0]*len(g_active)
     vbin=ceil(scaled.numerator,scaled.denominator)
     corrections=[]; mono=True
     if -N<=vbin<=N and F(vbin-1,cd)<vcrit<F(vbin,cd):
      end=F(vbin-1,cd)/ell; low_den=t.numerator*d*shift.denominator; high_den=t.numerator*d*end.denominator
      low_shift=t.numerator*d*shift.numerator; high_shift=t.numerator*d*end.numerator
      for at,parent in enumerate(nodes):
       sb=substeps*parent+1-vbin
       if not -N<=sb<=N:continue
       lower=t.denominator*parent*shift.denominator-low_shift; upper=t.denominator*parent*end.denominator-high_shift
       mass=max(0,pin_cdf(upper,high_den,-1)-pin_cdf(lower,low_den,1))
       if mass:
        corrections.append((at,od[sb+N],pd[vbin+N],mass))
        mono &= mass*WSCALE**2*(1<<14)<=f[sb+N]*g[vbin+N]*norm[at]
     self.monotone_main.append(mono)
     on_e=pack_even(v*Q for v in w[::2]);on_o=pack_odd(v*Q for v in w[1::2])
     off_e=pack_even(anti[::2]);off_o=pack_odd(anti[1::2])
     plan.append((fe,fo,ge,go,on_e,on_o,off_e,off_o,corrections))
   total_weighted=[base*Q*Q*factor for base,factor in zip(B,norm)]
   if isaux:
    return (plan,quarter_pack,quarter_product,quarter_width,quarter_full,norm,[2*r for r in norm],total_weighted,quarter_parts,len(nodes),quarter_weights,residue)
   return (plan,pack_even,pack_odd,half_merge,half_width,even_full,odd_full,norm,total_weighted,reduced_parts,len(nodes),f[::2],f[1::2])
  self.A=setup(True,Na,A)
  self.M=setup(False,Nm,M)
 def apply_aux(self,m,p,plus=False,previous=None,round_number=0,only=None,cross=None):
  rows,pack_lanes,prod,width,full,norm,twonorm,total_weighted,parts,cnt,weights,residue=self.A
  outm=[];outp=[];both=[]
  fneg,fpos=(m+[0],p+[Q]) if only is None else (m,p)
  limb_shift=2*9+40+LQ
  need(WSCALE==21<<9,'split denominator scale')
  if cross is not None:
   gneg,gpos=cross[0]+[0],cross[1]+[Q]
  def cross_half(a,b):
   def pair(r,s):
    return int.from_bytes(prod(a[r],b[s]),'little')+int.from_bytes(prod(a[s],b[r]),'little')
   if residue==1:return pair(0,1)+(pair(2,3)<<64)
   return pair(0,3)+pair(1,2)
  def square_half(lanes):
   if residue==1:
    a=int.from_bytes(prod(lanes[0],lanes[1]),'little')
    b=int.from_bytes(prod(lanes[2],lanes[3]),'little')
    return a+(b<<64)
   a=int.from_bytes(prod(lanes[0],lanes[3]),'little')
   b=int.from_bytes(prod(lanes[1],lanes[2]),'little')
   return a+b
  for idx in range(len(rows)) if only is None else (only,):
   gets=rows[idx]
   if self.reuse(idx,previous,round_number):
    offset=idx*cnt
    outm.extend(previous[0][offset:offset+cnt])
    outp.extend(previous[1][offset:offset+cnt])
    continue
   fm=[pack(map(mul,ws,get(fneg))) for pack,ws,get in zip(pack_lanes,weights,gets)]
   fp=[pack(map(mul,ws,get(fpos))) for pack,ws,get in zip(pack_lanes,weights,gets)]
   neg=[f-p for f,p in zip(full,fp)]
   pos=[f-m for f,m in zip(full,fm)]
   fzero=[x-m for x,m in zip(neg,fm)]
   if cross is None:
    neutral=square_half(fzero)
    ac=(square_half(neg)-neutral).to_bytes(width,'little')
    bc=(square_half(pos)-neutral).to_bytes(width,'little')
   else:
    gm=[pack(map(mul,ws,get(gneg))) for pack,ws,get in zip(pack_lanes,weights,gets)]
    gp=[pack(map(mul,ws,get(gpos))) for pack,ws,get in zip(pack_lanes,weights,gets)]
    gminus=[f-p for f,p in zip(full,gp)]
    gplus=[f-m for f,m in zip(full,gm)]
    gzero=[x-m for x,m in zip(gminus,gm)]
    neutral=cross_half(fzero,gzero)
    ac=(cross_half(neg,gminus)-neutral).to_bytes(width,'little')
    bc=(cross_half(pos,gplus)-neutral).to_bytes(width,'little')
   if plus:
    dp=coefficients(square_half(fp).to_bytes(width,'little'),parts)
    both.extend([Q-(((total-y*r)>>limb_shift)//441) for y,r,total in zip(dp,twonorm,total_weighted)])
   ma=coefficients(ac,parts);po=coefficients(bc,parts)
   factors=twonorm if cross is None else norm
   outm.extend([((x*r)>>limb_shift)//441 for x,r in zip(ma,factors)])
   outp.extend([Q-(((total-y*r)>>limb_shift)//441) for y,r,total in zip(po,factors,total_weighted)])
  need(all(0<=x<=Q-y and y<=Q for x,y in zip(outm,outp)),'auxiliary tail simplex')
  return (outm,outp,both) if plus else (outm,outp)

 def reuse(self,idx,previous,round_number):
  if previous is None:
   return False
  y=abs(self.y_axis[idx])
  interval=1 if 100*y<=self.yd*self.core else 2 if 100*y<=self.yd*self.middle else 4
  return (round_number+y)%interval != 0

 def apply_main(self,m,p,n,u,previous=None,round_number=0,only=None,raw=False,cache=None):
  rows,pack_even,pack_odd,merge,width,full_e,full_o,norm,total_weighted,parts,cnt,fw_e,fw_o=self.M
  outm=[];outp=[]
  fneg,fpos=(m+[0],p+[Q]) if only is None else (m,p)
  gn,gp=(n+[0],u+[Q]) if only is None else (n,u)
  limb_shift=2*9+40+LQ
  need(WSCALE==21<<9,'split denominator scale')
  correction_unit = 441 << (limb_shift-26-LQ)
  def odd_product(ae,ao,be,bo):
   return int.from_bytes(merge(ae,bo,ao,be),'little')
  for idx in range(len(rows)) if only is None else (only,):
   fe,fo,ge,go,on_e,on_o,off_e,off_o,corrections=rows[idx]
   startrow=idx*cnt
   if self.reuse(idx,previous,round_number):
    outm.extend(previous[0][startrow:startrow+cnt])
    outp.extend(previous[1][startrow:startrow+cnt])
    continue
   if cache is not None:
    left=(idx,id(fneg),id(fpos));right=(idx,id(gn),id(gp))
    fs,gs=cache
    fv=fs.get(left);gv=gs.get(right)
   else:
    fv=gv=None
   if fv is None:
    fv=(pack_even(map(mul,fw_e,fe(fneg))),pack_odd(map(mul,fw_o,fo(fneg))),
        pack_even(map(mul,fw_e,fe(fpos))),pack_odd(map(mul,fw_o,fo(fpos))))
    if cache is not None:fs[left]=fv
   if gv is None:
    gv=ge(gn),go(gn),ge(gp),go(gp)
    if cache is not None:gs[right]=gv
   fm_e,fm_o,fp_e,fp_o=fv
   gm_e,gm_o,gp_e,gp_o=gv
   zero=odd_product(full_e-fm_e-fp_e,full_o-fm_o-fp_o,on_e-gm_e-gp_e,on_o-gm_o-gp_o)
   ac=(odd_product(full_e-fp_e,full_o-fp_o,on_e+off_e-gp_e,on_o+off_o-gp_o)-zero).to_bytes(width,'little')
   bc=(odd_product(full_e-fm_e,full_o-fm_o,on_e-gm_e,on_o-gm_o)-zero).to_bytes(width,'little')
   ma=coefficients(ac,parts);po=coefficients(bc,parts)
   if not corrections and not raw:
    outm.extend([((x*r)>>limb_shift)//441 for x,r in zip(ma,norm)])
    outp.extend([Q-(((total-y*r)>>limb_shift)//441) for y,r,total in zip(po,norm,total_weighted)])
    continue
   lower = [x*r for x,r in zip(ma,norm)]
   nonpos = [total-y*r for y,r,total in zip(po,norm,total_weighted)]
   for at,fi,gi,mass in corrections:
    fm,fp=fneg[fi],fpos[fi]; gm,gp0=gn[gi],gp[gi]; gzero=Q-gm-gp0
    dx=(Q-fp)*(Q-gm)-fm*gzero
    dy=fp*gzero+(Q-fm)*gp0
    lower[at]+=dx*mass*correction_unit
    nonpos[at]+=dy*mass*correction_unit
   if raw:
    outm.extend(lower)
    outp.extend(nonpos)
   else:
    outm.extend([(x>>limb_shift)//441 for x in lower])
    outp.extend([Q-((x>>limb_shift)//441) for x in nonpos])
  if raw:return outm+outp
  need(all(0<=x<=Q-y and y<=Q for x,y in zip(outm,outp)),'main tail simplex')
  return outm,outp

 def advance(self,m,p,layers):
  n,u=self.apply_aux(m,p)
  m=m+[0];p=p+[Q];n.append(0);u.append(Q)
  msize=len(self.main);asize=len(self.aux)
  ind=list(range(len(self.y_axis)));ready=set()
  for step in range(layers):
   for i in ind:
    if self.reuse(i,m,step):
     continue
    a=i*asize;b=a+asize
    n[a:b],u[a:b]=self.apply_aux(m,p,only=i)
    a=i*msize;b=a+msize
    mm,pp=self.apply_main(m,p,n,u,only=i)
    if i not in ready:
     nm=[x if x>y else y for x,y in zip(mm,m[a:b])]
     np=[x if x<y else y for x,y in zip(pp,p[a:b])]
     if self.monotone_main[i] and nm==mm and np==pp:ready.add(i)
     mm,pp=nm,np
    m[a:b]=mm;p[a:b]=pp
  return m[:-1],p[:-1]

def root(s,m,p,hh=800,hu=800):
 from fractions import Fraction as F
 n,u,b=s.apply_aux(m,p,plus=True)
 C=Gaussian(2,26)
 S0=1<<26
 grid=40000
 extent=4*grid
 def mass(left,right):
  return max(0,C(right,grid)[0]-C(left,grid)[1])
 def points(step):
  # Refine where the Gaussian places most of its mass.
  return list(range(0,2*grid,step))+list(range(2*grid,extent,2*step))+[extent]
 ux=points(hu)
 rx=points(hh)
 xs=sorted(set(rx+[-x for x in rx]+[s.T]))
 wr=[mass(a,b) for a,b in zip(xs,xs[1:])]
 wu=[2*mass(a,b) for a,b in zip(ux,ux[1:])]
 mr=S0-sum(wr)
 mu=S0-sum(wu)
 need(min(wr+wu+[mr,mu])>=0,'root mass')
 wedge,ar,pr,pu=integral_weights(A_WEIGHT/2, A_WEIGHT+T_WEIGHT/2,T_WEIGHT/2)
 den=wedge*grid
 columns=len(s.aux)
 x_min=s.aux[0]
 y_min=s.y_axis[0]
 y_max=s.y_axis[-1]
 y_rows=[j*columns for j in range(len(s.y_axis))]
 y_lookup=[y_rows[bisect_left(s.y_axis,j)] for j in range(y_min,y_max+1)]
 full=(Q,Q,Q)
 tail3=[]
 for z,(neg,pos,pair) in enumerate(zip(n,u,b)):
  c0=Q-neg
  b0=min(pos,c0)
  tail3.append((min(pair,b0),b0,c0))

 d_factor=ar*s.d
 x_selection={}
 for r in xs[1:]:
  for v in ux[1:]:
   # The mirrored values use the lower endpoint from the preceding interval.
   d=r+v
   if d not in x_selection:
    x_selection[d]=bisect_left(s.aux,-(-d*d_factor//den))
  for v in ux[:-1]:
   d=r-v
   if d not in x_selection:
    x_selection[d]=bisect_left(s.aux,-(-d*d_factor//den))

 def line(rv,pv):
  ix=x_selection[rv]
  if ix==columns:
   return full
  y=-(-(pv*s.d-s.aux[ix]*den)*s.yd//(den*s.d))
  if y>y_max:
   return full
  if y<y_min:
   return tail3[ix]
  return tail3[y_lookup[y-y_min]+ix]

 full_pay=2*Q*Q
 tot=mr*S0*full_pay
 intervals=[(a,b,w) for a,b,w in zip(ux,ux[1:],wu) if w]
 for r,rw in zip(xs[1:],wr):
  if not rw:
   continue
  p=pr*r
  q=mu*full_pay
  if r>s.T:
   for low,high,w in intervals:
    A,B,C0=line(r+high,p+pu*high)
    D,E,F0=line(r-low,p-pu*low)
    pay=2*(C0*E+B*F0-B*E)+A*(Q-F0)+(Q-C0)*D
    q+=w*pay
  else:
   for low,high,w in intervals:
    A,_,_=line(r+high,p+pu*high)
    D,_,_=line(r-low,p-pu*low)
    q+=w*Q*(A+D)
  tot+=rw*q
 return ceil(tot,2*Q*S0*S0)
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

def collar_start(conv,k):
    rows,pack_even,pack_odd,merge,width,full_e,full_o,norm,total,parts,cnt,fw_e,fw_o=conv.M
    m,p=[0]*conv.msize,[Q]*conv.msize
    shift=2*9+40+LQ
    correction_unit=441<<(shift-26)
    # The minus tail depends only on the old plus tail, and the converse.
    # The all-plus seed makes the minus tail unchanged in the first round.
    for step in range(k):
        plus=step%2==0
        child=m+[0] if plus else p+[Q]
        out=[]
        for fe,fo,_,_,on_e,on_o,off_e,off_o,corrections in rows:
            e=pack_even(map(mul,fw_e,fe(child))) if step else 0
            o=pack_odd(map(mul,fw_o,fo(child))) if step else 0
            if plus:
                conv_bytes=merge(full_e-e,on_o,full_o-o,on_e)
                data=[s-x*r for x,r,s in zip(coefficients(conv_bytes,parts),norm,total)]
            else:
                conv_bytes=merge(full_e-e,off_o,full_o-o,off_e)
                data=[x*r for x,r in zip(coefficients(conv_bytes,parts),norm)]
            for at,fi,_,mass in corrections:
                data[at]+=mass*correction_unit*(Q-child[fi])
            out.extend([(x>>shift)//441 for x in data])
        if plus:
            p=[min(y,Q-x) for x,y in zip(out,p)]
        else:
            m=list(map(max,m,out))
        need(all(0<=x<=Q-y and y<=Q for x,y in zip(m,p)), 'line tail simplex')
    return m,p

from array import array
from base64 import b85decode
from hashlib import sha256
from lzma import decompress
from pathlib import Path
from time import monotonic

def semantic(part):
    from collections import defaultdict
    from fractions import Fraction as F
    from functools import lru_cache
    from itertools import combinations, product
    from math import lcm


    S = 1 << 58


    def ceil(x, y):
        return -((-x)//y)


    Q = 1 << 80
    X_DEN = 1 << 51


    def ln_interval(num, den):
        """Bounds for Q*ln(num/den), where 1 <= num/den <= 2."""
        a, b = num-den, num+den
        aa, bb = a*a, b*b
        ax, bx = a, b
        lower = 0
        upper = 0
        for j in range(1,24,2):
            q, r = divmod(2*Q*ax,j*bx)
            lower += q
            upper += q+bool(r)
            ax *= aa
            bx *= bb
        return lower, upper+ceil(2*Q*ax*bb,25*bx*(bb-aa))


    ln2_lo, ln2_hi = ln_interval(2,1)
    inv_ln2_lo = (Q*Q)//ln2_hi
    inv_ln2_hi = ceil(Q*Q,ln2_lo)


    @lru_cache(maxsize=None)
    def phi_units(num):
        if num == 0:
            return 0,0,0,0
        need(0 < num <= X_DEN, 'entropy argument')
        shift = num.bit_length()-1
        ln_lo, ln_hi = ln_interval(num,1<<shift)
        e = (shift-51)*Q
        log_lo = e + (Q*ln_lo)//ln2_hi
        log_hi = e + ceil(Q*ln_hi,ln2_lo)
        return (-S*num*log_hi//(Q*X_DEN),
                ceil(-S*num*log_lo,Q*X_DEN),
                S*(-log_hi-inv_ln2_hi)//Q,
                ceil(S*(-log_lo-inv_ln2_lo),Q))


    need(part in (0,1), 'entropy phase')
    split_pairs = {(d,s):{} for d in (1,2,4) for s in (0,1)}


    @lru_cache(maxsize=None)
    def binary_units(num):
        a,b,c,d = phi_units(num)
        p,q,r,s = phi_units(X_DEN-num)
        return a+p,b+q,c-s,d-r


    @lru_cache(maxsize=None)
    def entropy_support(masses, paired, lower):
        units = binary_units if paired else phi_units
        if lower:
            return tuple(units(2*x)[0] for x in masses)
        lo,hi=min(masses),max(masses)
        if lo==hi:
            return (units(2*lo)[1],)*len(masses)
        origin=max(1,min(X_DEN-1,(lo+hi+2*sum(masses)//len(masses))//2))
        _,upper,neg,pos=units(origin)
        tangent=[]
        for x in masses:
            dx=2*x-origin
            tangent.append(upper+ceil(dx*(pos if dx>=0 else neg),1<<51))
        return tuple(tangent)


    def prove_simplex(linear, constant, terms, biases=(), symmetric=False):
        dim = len(linear)
        grouped = defaultdict(F)
        for key,weight,*cap in biases:
            grouped[min(key,tuple(4-x for x in key)),cap[0] if cap else 32] += weight
        biases = [(key,weight,cap) for (key,cap),weight in grouped.items() if weight]
        common = lcm(constant.denominator, *(x.denominator for x in linear),
                     *(x.denominator for x in terms.values()), *(x.denominator for _,x,_ in biases))
        c = S*int(constant*common)
        line = [int(x*common) for x in linear]
        keys = tuple(dict.fromkeys((*terms, *(k for k,_,_ in biases))))
        indices = {k: dim+1+i for i,k in enumerate(keys)}
        paired = []
        for key in tuple(terms):
            other = tuple(4-x for x in key)
            if key<other and terms[key]*terms.get(other,0)>0:
                u,v = terms[key],terms[other]
                weight = min(u,v) if u>0 else max(u,v)
                terms[key] -= weight
                terms[other] -= weight
                paired.append((indices[key],int(weight*common),True))
        terms = paired+[(indices[k],int(v*common),False) for k,v in terms.items() if v]
        biases = [(indices[k],int(v*common),cap<<44) for k,v,cap in biases]
        n = 1 << 48
        vertices = tuple(tuple(n*int(i==j) for i in range(dim))+(n*line[j],)
                         +tuple(n*k[j] for k in keys) for j in range(dim))
        if symmetric:
            mid = tuple((x+y)//2 for x,y in zip(vertices[1],vertices[2]))
            vertices = (vertices[0],mid,vertices[2],vertices[3])
        pairs = tuple(combinations(range(dim),2))
        splits = split_pairs[dim,symmetric]
        pending = [(vertices,int(symmetric),1)]
        checked = 0
        while pending:
            vertices, depth, shape = pending.pop()
            checked += 1
            bounds = [c+(v[dim]<<10) for v in vertices]
            if dim==4:
                a0,a1,a2,a3 = vertices
            for index,coefficient,paired in terms:
                masses = ((a0[index],a1[index],a2[index],a3[index]) if dim==4
                          else tuple(v[index] for v in vertices))
                support = entropy_support(masses,paired,coefficient>0)
                bounds = [y+coefficient*x for x,y in zip(support,bounds)]
            for index,coefficient,cap in biases:
                masses = [v[index] for v in vertices]
                lo,hi = min(masses),max(masses)
                whole = 4*n
                if hi<=cap or (coefficient<0 and lo+hi<=2*cap):
                    selected = masses
                elif lo>=whole-cap or (coefficient<0 and lo+hi>=2*(whole-cap)):
                    selected = [whole-x for x in masses]
                elif (lo>=cap and hi<=whole-cap) or coefficient<0:
                    selected = [cap]*dim
                else:
                    selected = [min(x,whole-x,cap) for x in masses]
                bounds = [y+(coefficient*x<<8) for x,y in zip(selected,bounds)]
            if min(bounds)>0:
                continue
            if depth>=48 or not pairs:
                raise ValueError('unproved entropy simplex')
            split = splits.get(shape)
            if split is None:
                def distance(pair):
                    p,q = pair
                    d = [vertices[p][k]-vertices[q][k] for k in range(dim)]
                    return sum(x*x for x in d)+(d[1]+d[3])**2+(d[2]+d[3])**2+d[3]**2 if dim==4 else sum(x*x for x in d)
                split = max(pairs,key=distance)
                splits[shape] = split
            a,b = split
            midpoint = tuple((x+y)//2 for x,y in zip(vertices[a],vertices[b]))
            for branch,j in enumerate((a,b)):
                part = list(vertices)
                part[j] = midpoint
                pending.append((part,depth+1,2*shape+branch))
        return checked


    def add(terms, key, weight):
        terms[tuple(key)] += weight


    def subset(mask, dim=4):
        return tuple(4*int(i in mask) for i in range(dim))


    def hb(terms, key, weight):
        add(terms,key,weight)
        add(terms,(4-i for i in key),weight)


    def case(eta, scores, b, o, j, t, old_mask, branch):
        beta,_,theta_new,_,kappa,*tail = scores
        xis, gammas, xi_rank0, xi_old = tail[:3],tail[3:6],tail[6],tail[7]
        deep_counts = tail[8:9]
        deep_biases = tail[9:10]
        reserve = tail[-1]
        out_horizon = o-2 if o>=4 else int(o==2)
        if o>=4:
            o = 2
        in_horizon = min(out_horizon+1,2) if b==1 and j==1 else 0
        a = (F(0),F(1),F(1),F(1,2))[o]
        kind = (0,1+int(j==1),3)[b]
        cz = 1-eta[o]
        constant = -beta*j-reserve-xi_rank0*(t-int(b==0 and o==3))
        constant -= sum(xi*(b*int(kind==i)-int(o==i)) for i,xi in enumerate(xis,1))
        constant -= xi_old*(old_mask-int(b==0 and o==1))
        if in_horizon>1:
            constant -= deep_counts[in_horizon-2]
        if out_horizon>1:
            constant += deep_counts[out_horizon-2]
        terms = defaultdict(F)
        biases = []
        if j == 0:
            cu,cv = eta[kind] if b else F(1),eta[kind] if b==2 else F(1)
            constant += cz+cu+cv
            line = [F(0)]*4
            az,au,av = subset([3]),subset([1,3]),subset([2,3])
            L0,L = (F(1),az) if b==0 else (F(0),au) if b==1 else (F(0),(0,0,0,0))
            wu,wv,wz = 1+int(branch==0)-cu,1+int(branch==0)-cv,-cz-int(branch==1)
            for i in range(4):
                add(terms,subset([i]), -1-int(branch==0)+int(branch==1))
            if b == 0:
                weight = -kappa*(1-a)/2
                wu += weight
                wv += weight
                for i in (1,2):
                    for mask in ([3],[i],[k for k in range(3) if i!=k]):
                        add(terms,subset(mask),-kappa*a/2)
            else:
                for mu,corner in ((au,1),(av,2))[:b]:
                    for mask,weight in ((az,1),(subset([corner]),1),(mu,-1)):
                        add(terms,mask,-kappa*a*weight/b)
                    biases.append((mu,-gammas[kind-1]))
            for mask,weight in ((au,wu),(av,wv),(az,wz)):
                hb(terms,mask,weight)
            if old_mask and o==2:
                constant -= theta_new/4
                line[2] += 3*theta_new/4
                line[3] += 3*theta_new/4
                biases.append(((2,2,4,0),theta_new/2))
            if old_mask:
                constant -= theta_new/2
                line[0] += theta_new/2
                line[1] += 3*theta_new/2
                biases.append(((4,0,2,2),theta_new))
            if t:
                constant += theta_new
                mask = av if t==1 else az
                for i,weight in enumerate(mask):
                    line[i] -= theta_new*F(weight,4)
        elif j == 1:
            w = eta[kind] if b else F(1)
            az = (2,0)
            L0,L = (F(0),(4,0)) if b else (F(1),az)
            line = [F(0),F(0)]
            constant += cz+w
            h_weight = -w
            if b:
                line[0] += -kappa*a+theta_new/2
                constant -= theta_new/2
                biases.append(((4,0),-gammas[kind-1]))
            elif o:
                constant -= kappa
                h_weight -= kappa*a/2
            else:
                h_weight -= kappa
            hb(terms,az,-cz)
            hb(terms,(4,0),h_weight)
        else:
            az = (1,)
            L0,L = F(5,4),(0,)
            line = [F(0)]
            constant += cz-kappa*(1+a/2)
            hb(terms,az,-cz)
        if o:
            biases.append((az,gammas[o-1]))
        if b==2 and o==2:
            biases.append((az,theta_new/2))
        if in_horizon>1:
            biases.append(((4,0),-deep_biases[in_horizon-2]))
        if out_horizon>1:
            biases.append((az,deep_biases[out_horizon-2]))
        for k,cap in enumerate((18,21)):
            credits=tail[10+4*k:14+4*k]
            for field in ((au,av)[:b] if b and j==0 else ((4,0),)*bool(b)):
                biases.append((field,-credits[kind-1],cap))
                if in_horizon>1:
                    biases.append((field,-credits[3],cap))
            if o:
                biases.append((az,credits[o-1],cap))
            if out_horizon>1:
                biases.append((az,credits[3],cap))
        coeff = -theta_new+(theta_new/2 if o==2 else 0)
        if out_horizon>1:
            coeff += theta_new*(F(1,2)-F(1,2**out_horizon))
        constant += theta_new*j+coeff*L0
        for i,weight in enumerate(L):
            line[i] += coeff*F(weight,4)
        return prove_simplex(line,constant,terms,biases=biases, symmetric=j==0 and b!=1)


    def prove_combined(eta,scores):
        need(len(eta)==4 and eta[0]==0 and all(-2<=x<=2 for x in eta), 'signed carrier coefficients')
        need(len(scores)==24 and min(scores[:5])>=0 and scores[-1]>0 and scores[1]==scores[3]==0, 'dual coefficients')
        total = 0
        for b,o,j,t,old_mask in product(range(3),range(5),range(3),range(3),range(2)):
            if not (j==0 or (j==1 and b<2) or (j==2 and b==0)) or (b!=2 and t) or (old_mask and not (b==1 and j==0)):
                continue
            if ((b,o) == (1,3)) != (part == 1):
                continue
            for branch in range(2 if j==0 else 1):
                total += case(eta,scores,b,o,j,t,old_mask,branch)
        return total

    d = 100000000
    eta = [F(x,d) for x in (0, 6000, 6000, 52792516)]
    scores = [F(x,d) for x in (6427690, 0, 3616495, 0, 5895441, -4469698, -9346728, 27190256, 34592891, 35844536, -13423863, 1048845, 1727650, 1364857, 9489348, -50623992, -4675462, -43951127, -7299297, 28777192, -2351967, -17985221, -9098469, 100)]
    checks = prove_combined(eta,scores)
    saving = sum(scores[i] for i in (0,3,4,23))
    cap = F(12257883580547827010,10**20)
    coefficient = 1+1/(2*cap)+saving/3
    return {'checks': checks, 'local_reserve': str(scores[-1]), 'saving': str(saving),
            'semantic_increment': str(saving/3), 'coefficient': str(coefficient)}

CAP_NUMERATOR = 12257883580547827010
SCORE = '1882842180509854608179226931/367736507416434810300000000'
COARSE_SHA256 = '262fc57fff3601a5514d49a3bb451a6156a99b925db4ed0cff3009388f29e1d9'
FINE_SHA256 = 'cf759c6ff17f5f2adc6de30c089fc60374db7bd8c97fd276cc48bd05302f9534'
COARSE_SWEEPS = 20
FINE_SWEEPS = 6
COLLAR_DEPTH = 5
PATH_DEN = 441 << (58+LQ)
PATH_SPLITS = 9
PATH_LAYERS = 3*PATH_DEN*Q*8**PATH_SPLITS
ROOT_UPPER = 8184
FINE_STEP = 8
FINE_OFFSET = 0
GRID_STEP = 380

def positive(a,b,c,d,left):
    if min(a,b,c,d)>0:return True
    h=a+3*b+3*c+d
    if min(a,h,d)<=0 or not left:return False
    return positive(8*a,4*(a+b),2*(a+2*b+c),h,left-1) and positive(h,2*(b+2*c+d),4*(c+d),8*d,left-1)

def path_verify(conv,start,finish):
    m,p=start;n,u=finish
    a=conv.apply_aux(m,p)
    b=conv.apply_aux(m,p,cross=finish)
    c=conv.apply_aux(n,u)
    M=conv.apply_main
    l=(m+[0],p+[Q]);r=(n+[0],u+[Q])
    aa,bb,cc=((x+[0],y+[Q]) for x,y in (a,b,c))
    W=len(conv.main);moved=0
    for i in range(len(conv.y_axis)):
        j=i*W;k=j+W
        v0=m[j:k]+[Q-x for x in p[j:k]]
        v1=n[j:k]+[Q-x for x in u[j:k]]
        if not any(x<y for x,y in zip(v0,v1)):continue
        cache=({}, {})
        A,B,C,D,E,J=(M(*f,*g,only=i,raw=True,cache=cache) for f,g in
                    ((l,aa),(l,bb),(r,aa),(r,bb),(l,cc),(r,cc)))
        for x,y,a,b,c,d,e,j in zip(v0,v1,A,B,C,D,E,J):
            if x>=y:continue
            moved+=1
            h=(3*(a-PATH_DEN*x),2*b+c-PATH_DEN*(2*x+y),
               2*d+e-PATH_DEN*(x+2*y),3*(j-PATH_DEN*y))
            need(positive(*h,PATH_SPLITS),'continuation residual')
    return moved

@lru_cache(None)
def open_module(name):
    from pathlib import Path
    from runpy import run_path
    directory = Path(__file__).resolve().parent
    if (directory / 'files' / 'moments.py').is_file():
        directory = directory / 'files'
    return run_path(str(directory / name))

FINE_KNOTS='{Wp48S^xk9=GL@E0stWa8~^|S5YJf5;6@%4hg|?<6#YQO!oss6+0k-x+KhW64s}a2jnv8t3h8VXS=${G;zpl?^;Ngl(IbN-{}gSvLR!7~bNVZ6uIUd#M1#l6?oTlRl^4GX-8rXMP7l<`=ZrhHiXU6!)@3!kz_2#pDFk}HL|}>dE+Um@j>Ir#mFX7y7To4?1--RJzd&8lwphR&Nb=Ad9ri)i^r)GZY?1V8oQKNnln+;a`Og9bbXwi;{L+N)@|>HgDq9EsuVP1r5+PPqo3L=R9{JQfauT$rFgpSgK=q0Q#ul;9^Z5-0=`|4H`LbytYJ0%rV5z0%ZpeOb;ZDOLBXRJM#=2f>S&)Sf%Mi)xT_1=Q?L^D+3^oszm5apAIIrEJni?;yqmk$#e)Vx%S_{*u%eS#VCx>4v0_{Fb#}^kJygISdR#1k4So0F?+vNH&oeMbS;NdRT86}C0E4kUk2yEbGVPOTR!Sx|fYh2i_JbJl99ZoU-J}?yX+2BaDt*0fVkJ54!35aJ_aXn_4@~m_C%tVvD(d5oCp0OWHyN6X68=p~YVRmUN0%?<8qD0bAtoe31Q+R@c{c}F>*pN0}DLcU;)+6^mH87Bq^UFSs&ZPp?4{7O+!{2RCn-?bfByi@S3kaS92Q2j@y(k~)sd+fUq%K4*=H!l4+?}nlr2EoUgKNdv%rB0^BsTmML`FRhi{Kx019*k27?#~%vxY@m?1`|;hN)(p*-F9VJw(2Wp8Sx|M3Q46Oj<xoea&@8g8)U}s>G03!XGp=u9mIETJLL^^Lg{H7|=6-yG>p=Tdnv?-JsOLF>WeicmiptHvdR}=3-!d9!RO{f(oIAzH!m1%o*{<L62ZQDTV_ngys1D$SA(J5&&mFPohSz%w4K)z}}HMb|pCiijJrdBZLPtFz<Bt8UU3jrTUGbGjz9Gh9d3H(^`2_lPEd^zrvFrMuT-0%x4ydck{UBIT>gvf<6ruh7cX2${Ul7fSO)ek;J_CsQ32wlR%^`NyeW4%N%f;o*DEznT??6-omUXPC-r4qD!Qa8`65{ZfsyOo0IM0{Zr9#sF0O+0<rlUItoO-%)!Rsc9wft8ciVgDl|oz#&=iFq*Q5Pf~GxWgToFt@WxTbPlOcvEeJj}IoW~alhxtG&KXQJ0%&z>XjmH*>81~zR2Mjy<C(GC#0%0eqIw9o#E}Frc>4{ZYw8$z6>lQ%;+V>U+DUHNIkVB1r)EL%0;<HUN{J-&pwUJey#cWmX!;ISX<}hFd39kQTzI)jXO0eXW+J2kj(1!X&?zifkHyYyXn1;kV|hz~;*30E&c=w~MN5S}6=18?wISafMl_dlBYUv{GPG<uS5h0?Q6e${&A@sls{ohf3iBC^YTrtert&ElEp7o0V(|JJ*h4=%aX{c((yXEcv(ZV~yTjaFLR%Ycu8GB>g`v~Ti2H?Y%a0w0NI)E)aQaLNF&defXchy(<-;5P_B^5%<d+eT&Oa9nlMPBe4P{f5!I{7MHj-s59x+3j6p$HjMq*P7kl(}nO4dN9v)DLbW!K`a9%je6^xOfS^YjEowU%B4m^UMAkK#Y;ki8RmPDo!pL8_5CQsQA4RG(duN8O2Ggkz@9NTUP<BGf;_f{=iyZG@YOdgq;mT_o?i;~XM9yCT}XKYtY^BS38!c)4w`EDEjK4>-nAZ@mrAa;Y})aA*MV?dINco~E<Ud3@v{v0O~FD_9F?KxBGCoE|xdlTq+{d5T!!k@<bW_yTHdhGhg~dc~^!cWSlS3a~y0O#hU0@TC(A4MI~8emQX0NC1`_g9IgNqF{C0zf@E7-Ao_^;@&2xBeV(ilfS?uoX9!if!ybK`YpSBgxC_-ORq>2A*zp_OEQZ7;HM*{%>o4V@HKa;G<b0F2LZKCoAjR=hL9kU`Vv0;2!mz}51Y*o8ehV0sxCRj1qaf(*@Lt0`A_P)!eCYNJ87<tEy`)CyS3aFS<P*~vVOSTUK3`%F_<o%$9@R7Sc@xvpJShCMo~wbWj+T+s5v~_s<x|%oe@hnXSzb1H{s4#F`8}8X(laJ7R`xmL2_f3jlK-e#a(<MDY!ulPyMd6YUt;h`+Twg_7F2-ubp@t_69+a&igF@;DHB7CUpA_?s0KBo*)4QA;4zyKo-M9@-<o8jsDVKeBL6hni6-4Ib{aOWqJB~#~3RZEA1N1U<@0;sq+dCrYfvit>=54hnMasd_B3gzCV{XSefI_1V&ydZCv^%0b_^Ye)us%T+#%fy=w>xB~iTT%S)x1M6vj9tC>Bm3E&u{v%{I1tt7GlYPsMEmWpb)WJjoacE+akC^xH=KFuB(d)&;_l8kjg5Vf-a-aUaXJSKsx48vmGEja9^tOvyxoRbndlitL?m{AYxz5$K$uhy}3QX8!$-uMz|p`lG`@8v4$d_q-nEWc2zYP|XYK%>6!PvhCSi)NZbLWgX!4peKuW>`MAYk8W1)U!LefHsUXThYZKQyJ`K+T)}7X<5m+87dPt|BS)svH3)iu>oP;KZ{Ivx2+cg6_rsqeQgbJYH@kyC1jGJP(0k#`ElZj_EVBYX+kDMGdpH?ho&f+KLa*qHnq^dFl2<SVP&p)QrW?I6v2uU=-E6x@YDnhEdqz}d3}TT>lieoX%E|D``N&Z7D*+uX)6&XGsI~IrAZbA0LUJK+f)}q+Ze62Yz7xqX4Az*jEDbP@3ggZni4rtcq0U|N_xS$>}uw7$DE(QSyi%-Qux95zPMjG&_YM1G&?b$i1KwmW9u2tHjLcdd_RC<6}BJ!Z5s+wYK2YMW6u2^q)}ZfyT`;R)rWqK6<Ot7SBAIeu*(TamzH?(G1x7QPuI9G5vzfKZjy0=ne}ec4#>DKK$9p-9|f50q;~MUpsH@`ehO)fxO@{9ib;Vx;}^M(n(o!DI%q1ESyGL<<&^Il2-~)^LSwY1PTu?LJ>q~tHncI985rm$_(=|hUQZhqsW+wNG*qiJ2x0m&0;=3z>`km(C~GiQVSA{Ec!k4e2<h9eMUKuNi?b-}V^rp5?o{5DpunFw=4R)zx_L*yb{6fdO7=3rlwjyAL~<1yAD4Sq?)W(m$Mr8*VgTbD-az#=@GuUz(AwCn+afic{{gEIBgEmm1EL0)D0B+jfV_$vR9ZYA?y5=En`vO~n8PfRAh)_h9?0?Exa(XB_-DCmRA&`vrVP%C^s4#vphW~ehABNr_t0fsxm4SeQ7weF^!$7Q=IFB$hUHIf-0mmh0f0QCd_`Tb`0w2`QHM(@PUqdQM=rQ;DI93en}WX7-y#w$5xIjv(S53Dc;~Kysb$l|5G@qKWyXc|HSuh9Rpii)iIF07IlrKMeAGaZr4qKyF0$O}^*Co=S3xtyw`bIfhnpqGXzX2sZFK_Dz2$mC^RFwps@!RcScQ+UYz=5af>cTTp#*68x*@(*tsZ;a)+HbYDtYOKLMMimS5mo-<jp4oW>X?|OcP^YA=4w=mk)_HamVz$A#@(t7Sh8`OV>+8F?+hgl7t9Ph&F_xoIY2cN1CvEv+(fnSf*AZqsdVC@!~8=`@Fq$KKChy!tvXcqEfi}f((gU5x(5#oFArzS%T;rH63?76QyT~HFx72b(s@1EtE3X-0LARf{G;ht<4oe+P7PuFe%dlFwISv#%OATbpcoD)+jd4Bz+hMw2AO;{dz48>J3>$V|OW|0KeKT*my#-Dt+l2z}z%DZx*ry<`+hpuRmw25YqSF&vK({;(diC^Vok#()_`(O0{1J?yIiBB8!bO@``yo7t(xT8LNS6t#aj7{DNX&MON@MXj)QfSKiL0T09mF%KkRP7guiu2ns>#ZWU_Hcs-oW@KlxM!Rvc$f04?bY=`&GApU8z@t5i$m(y+*ai7AvQ9Mq@^Hqt69SG7>rl2DZzk;KLjm0AdOKkWBLDuIZD{6*+ENRRa$6Hd$ub`SNm($vZxG=({VDFl)a&nN5tg%(NYmesN&OsJS={uMET-;4Fu6fmix1^=delk1_tfKfS2Qo*bRLn9p1^@gNbck}3gQRv{?Cchws}PCYOifZ=l(O(N4HZ^bZ|)i`uv?lt?~AaSzJo4^w=@yWk48}CEE}@|1`YjHFMaTy^v&*iEuS`PzV?Drs@zK}jmp;;v{q>pWiT39WHQGah8U+Zu7^{<RSN(Fi<Ncp#g4*l5jxeKXp>6$VOl`Ia|zIfRXr5L<%j_X?<=8<b`2ePI*GX-hOM*%ru+?{*i2k19og``n2&lhXs}2<yGm?s=Idy!zq>#AvWOhO-|zM5=Uqlp;5ZRpcx-x9-#9>aMPr^|22H2}rV!)LYtsSjc2tUMw!_t2iJUZX0m4HF{XV6Y@e!}(#9cDPw6YA4)jTJ-HAcJ7{BNO}a9kdAW)jrW0fuk#hB!Lk=vWd$Td28wDNYay?z@eo&^+G3nqb|7fx_EXdDF{fRWUR#B+5-97{G}!V35>)sAa3Cey-EIw3k-3+GFhFXuQPn%VRZ@gp)ROJGi=B_J-1aCw<5Krelg6i#AdgrhCO){XL`9eoh=f=z^S9Vug!s|7T<}aIv6MI@Q3eN!gRm$g1YZl)C%*Kf7~R+M;{tzKI;Dp7u;=Y~a4yMG;=4-JY<?n!5bhjAq(F*0=!r_V~;Iz84*5GL1;@4CQy6gfZfeBdYyF=ZO}DOa@3=^r~18^1vBWpup0~eb#zZ9V;#J;MF9R&Hm8Hu7*W2esz-piZJNuWD7sB5ZJ{FFw#CvF#~XYR~<l@y4O)vE46GnZ2g&^)!JVFHfe9DP&H0R`_3fl6OjEfq3Nv+b{+^V7G5a=^2~eKso<2Ik|>7?6JgF5=>Vw!N1g<iPAns~c~Hx&hNfBC^wuNof<!JVlFe@@pjzIVam8%%Aq@7=C?*kx47B1{6I5VisZ`n&q=TQB!q9anj|EDzTAJsz810JLfaBmI7Sfz`DG;Pm5QCiI7Q9X~nf(1Q--O~nZjjOB5-5gARtyL6vN+#+ETYq#9wR3(9-j{!K-}@~0u?oANp-p)0V*>Wl31VBu3p9Z?8m|~>`?eHW+77xT}b0ps2%c5S6`nZ8i|DLFItt;aAUd~d*83<B`5XgX?QjE`WVbk>F`roZcO8qmLbTcWmWwm^#lMHNWqqa7Y3BB&5m91p<x8N)A2#@<BLRda&g(6;Yb9G>^LU0-55_Bpkn0!qq*hUjxKoah(jT&i@ZkyJy+{_JcZXZ`wbL4fA)@vi=fzXL~CCRAbb(wIM19XLE9cr$LN-xdMsY2b>z%YsO8sCVA+FB6hhC?gTL9y7l+5E`U^uuj_!=nH5YhtB)*S>o?r3-VNd(4up!0%zeQE(oF5+|26ZWN1Dxe5Z@v34x%km}n!PRjpLzBPt8F^&>*>h2KhM{b<TL=4@xp;^jd!)c@htp)Xa;!H09I^SsN1WqY*F#rUz6z|;xG-{af+nK?iKdq4p%2?(<=xwilCT?(aKBX;Fo&Rba%vI()i{Sh~WIJRdQBckphjISU;F{Qu$B+mZ)`c_QpTylVCxZY$<+e<V@%5#UZ9}ImMwuV^~{bBc+s;xM!XRHDc?2YGY*wda4&@RBqoO;V$v!eGdj2lSUE4H?OX!!n9pm!<9}um++t<>MNphyVT4Kkcc<+^4i&ttMz!x3{i3b$Q$=;^UTAxQEOns-0*?N8OLPZ`-qLXQgFvPZ@}9>ejerbHBOs$n*tw5A{lmiv70pga?z>O%gP|_Sf|*9`KwGH*xWIf-x^4w4f&0|<{xWdRp0kyzU~UrNUP5WtTN_#1aZ)?q+n5q{q+HdML}tj!>64(q`_$@(_|}FPm75iB#C?74dKMe>oNXLM^mtB;W47KXIl(4Nlg!kbk;qbl@Yh*Q2O0$r`AqBSz+KR^oq#RF9d4Y`Pe&FCW}7}IHM6r*~qxDeE03h^HBmiZMSK+XArNjBd<^?40qa;gPaiSvg_#j6cw3q%+e|zn&ME__OVc>qS>pSSAZm}D6eny>!oO6cwNMB8sf8L9TD^SoSgJziCdy9jYXzjr6O;gnhL?$k{)2&+$nL@O}gS6D0<DmRnY%-GZ5|WOUkQ>{8*Hx-DpOWUv4F^U84wEaSrc~75q9hHKa!SU{G;ev=<$sATe*Kq<T!wcO`PC#O4+ERyUvxl;(-AcZHM$E`b?|vF8>4X_XfQO>fOI3xUdwlRNCzZqQ$$r8Z^vD|6sArsr>}#k`0hQM)B!b<r;$tB;C}{U32g+1~o_h#b&sc^Fmmjo5(GS-}$}&r*wl21wos-QUr$>8~-UydJ;O%nbNWJ><b~ynFn6<Dgs}QKZMMiom0P12T?wxJE<!Zwd*(L&j^;I{RPl>~DpXl1lP}Q1;Md=ey?2weDYZz?`8QfWE>tgFhvKLyI7Zru9>Gr9Q<-@&Dlg1s78{#fyJ|d{|b2GIfPUm*Iy|`6O_1C>KLM{Doa?YHOg;M68-w`;?1n{F9^FiE89PtlsW=2ttd*+r0Cp1h4I?P;Od?NYcnd?|S|o7sE#C6xflCHnz-UK>qOIX22m9NBE1j$vel&3m&cxth`8I)0AL_8%!GYtqDF}=2ymiUmz$WE=`ahU8({cW-%5Ii4%PO0Vf-kv&Gu`1SV|L%fshl;v?-LnTIVJXWdt3C|<fxnn_pFzy)GBHI~V4F{h$Ab_VX1mmcoB59>FC&wLl+z~)TpnPD_L!k2l_#O&-GW^f}uqt&|o%JW;e2bbZ1g*xdhsDmX2sdVOL4zMe5d%f#TpOH(R)NorF(h#jH+fEEuG>eydX2_SAsGxJd|0{&BNhdN}Ya9inv2~h#3!At6{zkoT+_C@wO)AmON70_M00E;ZpNs(jzt|guvBYQl0ssI200dcD'
FINE_FIX_HEAD=';uAobo&f*=Y@kbd`n?TwO-TkK%c|vmZuo}-(4af6!jiWL*^OhK%`w8gB=3)OW(SS&R~;CAG_nG%5)X=4DILAx4RoAI68@dm7&pMs9x|VIbPn^R#*g3esb$#%)QvVm!brkJt69gmhGVY3Tr^xsHoETV`A>n;sm!yVR+^JdSEV0Khwo@;5H}uj8|P7*O0*z=695{yZ|<ul=@fE#q@uF(EM98wAioH7K{~|tWVY#dr3^cUU$(__Gz|N*D3?z%^!_zASrAeg{Q%N@urm$>BG%rKRAsVjRrQD<m&wS2k<y2`*Y9D4&uBQKL`>{MW|X_-U^6?9mvqXS4NCVMvQ5EEnl(uYTlM1Hv3s_K%8<NN?m*JbKj%uLC5oq9QL4j3<eVMt{#SV$yO2B369ht1-&?_OJM~Fuk(~$HqMr%@j!`M1)^K_tL2%oUtLV9Ke&%ar(e^9q$>yx;NVTlAI@}}g`7-{@^DUFfGUUfA8H!&MQ4x7a#pzbV1+xn+b_a5xXSJ%(I#*}=$U+h8x6oe&FeIkfYrN&xJd9QepDZV5cp)&GZSsQmJ<zvi)4pO~-Fw1dqIM-CQBa8ZZh}>9wpfB!9fU=nnnGs?Tg39E*37+8J!hk89`O?!Z$#--S?^PcZOR8SXXnq=OUp<Cfe#WFlF{L$HymUnapSyw!D5(Ro!9LuZBPp+NIGm5rvK@w0M_c4S<`n9aRDpF)2f^u>+!b}B&n-+^M!~PcPe3>I9N7Z_kl!z=5Ce1$Q|EH2+N-+TPQC3>icSPixoK5>yaN<XAWM|H-Kg_M{?3%4fZGF{O!wJW{}P&aL*ze{^xrYtsL*og<2{vyePJ?=(=R0?h`85@F*$!6mpiGM(b=idjS_zw4{$Cew`(xpU>h1`hEbe(%UP`U_@BQG5m@2U8~<4Xxf7vh>_AHR)Z?lF?fm>qA1px;>$W(4Tta{@_Eq^nng6<^i6x)^RO>^NFOhNHgUL8FA&7Bpfab@QHxR!7g1%}M*siXUYn%HRyJcW&bK;=zYq|qZyV*>)5q3)+cC@YOhrns_Hdr#2@i*yotV_5{Y!hyZTN3FPhoa0EKZSvlhM|q{~8np{A2f56#RJ=IXYXB5hP!I(#V=XTz?`@&4iQWZ`IZ=2Yam<f8p0HuA^6Mo}gMsub*yukUCpf5lC#6hQSz6W}tT{Iv0FRqX@Y(8BmFWr|Mki9NqRo*hlbxNYoDGwg|tYuni)ER$h<KoL9PAp(ZZKay~l8<swKz!H5^^{BBdtHG{bEgbcWkINa%o3O^n+#Ee_xbEPW83ziDep9Y_Mwv%x^eE#+BU^@X|?ZQ>}CcFf9e@?err@A*e6gD?nE(O+U5svQU>nG8zf^7vmkgS0T;8-W1#)6gNF`1NS#j}gva9MLag#~D?ummb<a?%(CxJgB}58{Iamnh+ADLLH=>DjpG=z~H9LP0^{V?yl@&}AX)==*@rGfQFKfK6hjPK0jRK#v(lHh<<zRrE|;&(`Ba2e%WjJasm>;{(AdEU}|R(*bK+moyFa-LAQB(o5|UZm8Y4N-d_sN4lm`Ps@&+hC@{l+1n@FUPV5>!<-95;hlY=rDab6P`uu;qf0a#ca75OKL2^1Ra3E&e?@A+nx&?KK>d^>rUwS{YX$_ihNmQ_fGgMufoQ)(jX&j)iyk+poxyalSIBM=tQ7w6z=}ESX7x5QOWBTCkoByLdj_sIvta=Z09#-E$ir-lE|X88eA-a>h2nZPA`3Euy|6ULS2MfdPEJah?qvXVJd&#smb#YD)HDcpI(SCF>7e@)#YKm<Y^pHlf41A7d7~$L4TM!H5BqCF5+)Y?zYxWaFaOK=@SL6mfvNix|Djmy=ERbOb$R(%4f9QxTu!+sWsns1c;MF(_?2Y+)f@mtj82KuzzVIm$S~p1<=M!b&Z=d^JAfL-+xHzp5qjiKUR@PvelvyF4gU)%4qP1_@O%f4VntGQKCYs;qLO9$Q)8XD$?AtStajy2;ggtD;dD{06T$Io+FpM2X@jI>>osIJDp}RxaRXdF!oG3eX0m|fRTV8&=4})mp;aGc^B<B*j^!sPI_~UQD(4eOY{9Vt8~#R|Bj`>^vN2ZkM+^5FH-GH80$T(@@Su0p8BehPpfsdsWz3zNwU6&;g<I{_y*6Mb3R_hIy)rL)XP`+OxWUHkjD-H!jpQ`tDVjo`Id|p*KdnEM7v>8Oqtf}7$&TbnC<vjF!A6zlhVj%aDP=t-o;Ael?EPf9i%|wK3;I25PWg5g%E}jN)Y04U(r)jL808(&o!Uk*5r|ZKTO(m<SfjvZGX)!Qts#8G8*i-z<?gM#8F<_O`uF5J_Y^%?2{$J^p=BrTFX8a{^slT@2V!>Zj>V87wJf3NHD&k+Q|`9Ba;L7JNa^?Z<fr`Chji7ra4~m}r0(w}nYh#|#@ppcs@XM+W5GZ<WkD0tZeddyN8o#O3X;ueI8m_x6QxZ?a=J=k5GJat5iru06a4`~Yk+$sZP`@Y)aqE~Rtl>sd)_{iLQ3HpV-u(jJ>G~M+;UnhyZ*%65ZHc+ehGu3u}RnPD2%ZB<nWT^Mirc(FJ+j=gQTVY+@eT5aIv2(D^q#Z-{X+XNc!1|Z>2&Y;3-wc#zyZw9NK72%kTZY%o~F3cxA;3G6w-gb1*4RmWGT~Q?&Fa&8BvWROZWr$qG+F)0CgQt@bieBBajc#BPcMDU~E>8KSkK@2WESJQD0Fe>fp!2)lOm3SQXe`<&ie%!SY&PGx8L1g#hiF&>N@9Wq*ut)J>KKIW!3=(K|^=bD_Jxqhf|K9<K`tBvZii@aT|Jf(x1=@%c3(dF3WU}ygi&aeayn&6|NY-uqJ<~J)REYKd$>`Me%wM@4l^+lM(_Qq;?tK0ZAj&o#_BraLsCkPh<$zJ*-s~XdiX@T#Eeld#4?6tY@6l`7Wn!OL#qv2Nk;12~jnVHFJOGeK892*V{QJpq*Fu`$Ya<Y|&%b^-#bjqJGRFByJzYf-*zX2_6XFOHPk<fAkHmD30NCnWKStjx`OQ(Db2eF`1M8v#lC_+esUd;M_;MaKIq0A<58z?oIR*yXREw#$FXF&I2X9n&FCu49J<u_iRPujNr5JI{|x#|6vaWfIS*t)5Zu4~L4TV|<vNpR4}xuuZDI^7!O7@;LaOC<Ye1+&T{%@NlZ{>*JF#NqQ$I6+mcYBqZ${dL`&-)K^B8l)dxH@Ub*7Kj+ksX^<dIS3s##F1?+Jm#RuKzl}*E-5j>VA*b#Ou)CjUApn*ruZPcChXX2R)a&X;a_$de#^`m`3c7kGkBD)B&W)q8NNsbe%JidH3iy<Ls%wq<i(zoRR=dvNY2UdmeIL<nU69r(>n1A+|;7*3E~Q%VVTK#T`;@8T(rvTH(5W#JL}PlD*VLnLH=zMqu+H4dRz9Au9FGWMc`I`=|Lnw`|&n>SP-yOj>E2wzZoO1o-0SG|FqEya;pvzNcp)cX0EvJS&=y(UIaNYG&)&fW<!c9L`9f{TVw$AFE!cXhzaLEF<3)t9a|j-()&V+VQ~;Q!rFI`3S}dR4WMJ<<1^ox1k`aFzl`0Z=zT4Xx3F_j-H~CbySyj6m8UnARDMU$u9FB=G|73F{OIQnFN6<Vl)$cUX4ZFN!Q+$HRYv*JA3mmQjDDO6_hH2DR`!?UIhw=%&eE^A-^=NN)i2y3##UI1T`_yKb(sF=25^T9GK<IC(chU6Lmh+v8Ur8$Qwi5`n{gd>8Sp|ZTeB^(@j2q}To?k}?y@sIwE|!b@qbX!2^i(6ORb=2)PCavHw&GakdJY5tqCjwm_5I%BB+NAmNOuNi1qb>l-EUsLiUAO@eS7Z$+r4c2d?cws9`?;$Tvh7=+4jP*v~eeG!Yp+`6hw808E``LMj}O{qoW)J5v@)c7@P=t4LI7AD~87Ts)7YHOI9D-H$w1`#7C<ximf%1_@VBGRGnl`B#=k<=_<<Bi-l*yik0p_zM^>#A(*R2w`Slqr)2aZA^Qtwh7Z?YktzF%{<g5Nd10ql0q5Wzf7|cHnL51Mxey?4vM2Gt6?$SX+ohJJ1mbh_d^Jq-Nk@nY+Ux*Op{;xZu7-^#U|2T#b|8M7yp2qb9vHyB^P$rYxdaeGh-7D?u{N7BJK?<UvhKA;dS;=QjLQlEarYR-A_%j@%)P)G11J#Ms2^h0U(Y*n9i+Cf2Ev?M%;rXg)0IEDJ3AJwm5RXcNmVOyyPQt*z1-5Rsw&PWryPU`HpJcdeFUd(RZ6Ms8PO!tMSm|0ao62$_6f_4?Sd7tO=$h9Z)7J=Bx0!XoiT3WLa!^4!<OUi)D!vOgP=Q{tNB&n(Sn(umx8xTmZ9v-#`Kq3N-W{4GJ^}eFW<hOv6SGvySDdH6h{(zR0nL7}YkL=PaA|^_#xys@-#c5@8xQF8i;jY+FfMO&$joIt}!Q4*uX<y79TcCYO#E$ypkTZ+(N2fN)cR>;Vn#ECYvP7R`(}a}J>A3J*heW|Ud70?Y-k98i5kaDu=;AZ1t`h2xDxS4j}rLDbH4)QTV#{9hurSGiB*l)p;b-C>L}kScf9^<IFY-Ew^&*@bn!t%oPjWHGx~<}(7Z+Cg9chduite+%j*t}}k8fFC35)L|wcN^5CMW1hP9Ub?3u_T>P9#(I=b=!Wf9Va%sW=sWpCWlzH8u&78F!Rr6x!3!hVDYa!g6X7HjrBUab2>V>?uD)DAG|Pd35(&Tn8R!eCbooToCG<spZ@>{S(D2Ei?#ewVI7??#7;T*1G*a6ME1@RO9Cg^V0N0SQJ0okQ+HRqC(Eac#>Wr+avzATC!Yn)}T4l)R&_Eu)!4qQ2ukS_4AeTo?xys*#8ltiS>d^CRr_bEaj2pJxcKUF^@qX`Sb{TyFQrazMMCL(#4X?sm+nLyqIiDe1K|-|!%xFQ|M1O8KNF{X12s-mD7C;%hQ60M{6uNeMTPy7#Kz)^s;0Y5>4G{P^9Q2t35v2EO-!(mB-X@{J+nyP`=Sj#>gdQN|nJ$Nos|a6wQXj-4DN(>a>b8mSX$^Cjf9&-2vBq;XV?dU)PY`!XVG||)1Ib2@7<i`44p55{u>hlQh5@a9VK$Zsn;fsr(TT`*XQcC9`7z~1Ujex_(=r{IEr0TH2OucTX^ldQxmJJQ?S<f=HU%$A`JE&NXdz*@$q+t+?K~L|bJx!(VK$!DFebZBw%rY2eatpr1Eum%6&qY}v2P~Xr+1DbLUY-|xJvZUty*V_Z@LmT(;B#d6t?B-G=yf8tzulcIF`>Hr%wmYtO@ee$^a)H4=O!`qIhDP-|QN(O>iH?()0Tv3Y`deUNg4q7qH6+LG<g%$YT}(i_8(g4M4IemMovK@XPq2PbgM>c>5i(+sj|v(AdkAX|OHP<mkXbO_K?S+y%$rL;*Z!JuLIE(bUXY5>*CPu&|K!)8`wjZr1)J1Ppf1hVNV0!_nmuyZ5@1ShyJO&*4%!3fg=%8BOYcwtlVR8X5tBiY+-M-nc%ZQ+2(K!m9OtU)W1p$!?NuapT(~XobD*I-!Kd(wgFmWOREmJgy!$Ud9zGAIZX)tv8cO8B2hqU6NYHpUM7HCa?sG?&553S3taIa&PtYz*DM><GfSS+>iBgk;H%$!#QkTv7S8uzrhwDiJf^fCQ5GA+wD&j<Ek{^tR1_(F65gXuU^ENfF_f_@^Vz&kQ4h;5W<Qx3oMwXN8;VBbgx9EICzs3z+;bFbLOCO6wCgw*pjPq)5c+k9R81Wq!+_IstAE$`Qp6mc3_ZVTXkClq+*bR$NY^ZtGYDLgm2ne8g_Z5bEye?HN9|N;J{9eeaW3VuE%Fob8XAmlxy^fdiGsSNx7<LY0VJx$ww=_Oa6EY_j}CsyN=+$9s6yLIPh_lRu(s>h~@O-NicFwX&ordR5@Xt0TPK}elMKWVK-vSE*bgp85@Q9QFf==3!^}53objL_hePrGGmA7T;Tc<rogmLlff=UH8<n|ua_Slq<sJ<TM^pJY3-sc2>jAj_n9GH%uif1HBPhVYK1|O{dJM(&O|H_4Zfn{XE+@dHWXP(t_wwiP#M-%5*#Y>P+$HFsuf7N=%RMY2Ahnswi+{+KTF*AH1sK@3BmYn0BztNtP|&zX-c^e+|c3?#O<jiB*M;yrSmxq0y{qf?=HKj599=q6}ACKv%t4BXKtmW6%{S870pA_-ESO4$ZEBlg4MlWG&SO$y}Ib|Eaau(4bc>>Hp*-jY~Fo$&#}y3P%vBNQ{WV-u@R1{e+pq52jNy7AVgvtJbd|v<~KjTVrIr;{okk5KJwL~Y(DQTUIk12a;k3r(>p9xQ#kUHBA=h|dSJaw2F=k4e7!>dCPIaD7vC&lHbR&da)i7br<=-T=s_s+wSM{~^(GL4j-#{=`}NRZ3!=&_)pzuNGA?_bvB7$HRLf-$imUw_o{^v+Kk`<wI#*e0N)W{RQUB$l?qUIdHd_A6(V}fl0<X74&F6F1V5-E<K;JbGU{tiZ&U+j}ixtY$W_F)SV`t)MCptab{MzcK$z_=47uzA8U&d@a>a{yPt<l1J<yJImlqNL2?k@Qs-~5+ZXYGc>3ddiC<Mj+HPsZ19Wu#^1drhFI%DShIKI}<7Ohk7vWqciqSh0|V0d7#obNVb4=$$(~c?Wgxm=l}wzpw~@+WY23KlQ$BnTh|su*akLiyG;r{%G1XCA=BtF5ZvN`q;VhFR-$~eav*vH7w?3A&rq{bLoOhUzMr6O0VD#>IbS$H}$T%8}7DI3nDP<@rOt`5%CvBdniatVxWniczK4Y-~3r(3vk=vmIxGnH*d`Hhc6Z?p{G>wU)Q?~fUEJ$a2W|VZw_eJbh~mCk0h@&1u&c`18H^m1}3X~kNUV|s!C>p2Bl2%g_Id^(QvF*pTjL&oMCCjs^@Ch_P?W0j!@GDQg*G2OIvEv-H2eK$u3V;<rqC&TfJEhBX*=xh|UVUYEy<Esrk~ZhvR-fk3fd|-~`BF+ux7W9rFVX0ZBjLxd`9&Q8-;fVFXGzd#a3>*C|3<7I@18N?HxNG%prw_1xksH)3SzWkNKkz8(n@1uE<uXNV@?NcVi&w;dI5i$rUxbFrAa`}W$NesBEFpB-I!Y2t&5Dx+K7Q#g@7dMT;SL|LPeO+HchXWs*2sU(P+aE+*Ql*bw?6zuJs2O{<Y)7y(QT~zZg3zPLePv(8Bn4w~cO1Ee`G0`GrGIYiwsYU9#ZE4<vr6uqGboRRug=n&Su8rPw2}hc#S~yc)`VGYu7$ISpC~nw&TssFeaOcwySXh2xD`~l~+xo^(hkQA*lw{tPY$2P+)i4hCr9tn(5?^DFFfiRE+T0LMPH5wd=d;mDuvKL{T_;&0rB@BE;kMyS_5bNFY=i6Y+_S?F5q+BR@yhYj`$*c*5|#*xurl-%(Q^e#A(km(8xdP*lr&^}#Udi1+-dKcfMU2USP!~)r~*(2n(hG2O)|fpI8DQAjY@!3(i07!1UN!X<svMI%5`$((d#nLm0R>^F@c>LZbK;yf8LS%6OM#&gjbn6P0#kah6x(z2#;4JHbjqrq6$ZZk=T~ID0};G$X<89u~$gJ+fcV7HCvjmvS%$L(GI^FZ#LP5q<Hqse|9U`=Jg47y88-mF2zbfkjqwcak9oJn81n1o@0l1&EN}p$_5zUHghjZ73|cVa+2`WkT>}dGBYQksj+NS6m1Te=-n!HNOL`%v1H`L9obW-73$V^X(l8`Yj~&`434VH_42vj5F|P(zhFK4jJj^^0mky>hgmvSDN=6&5FTf4>_=-dhLvlK)$q>9c$WMd;dWF8#S&cQ-u6%ALh0LFkGAdBg=~N+qA?JBAPXcwaHK(MfihHqj%>VaZWh4;f7MFXRp~<s#B)8d`8eBOFPy>?+6lqN4#}XgzHXAv?c?aWPfavcEVF6@p5pR-Kx`lVvX00CF;9GxQr&zGSCZmveJ?GdmJQ^4KRpJ>vlKsQ#{pVJ4W?F)kdqscA{Fb_iGG@<v0IigEKJncwBw#>lcAY+{Vx_i0oob;(CeeW>!9AmI|Uxj5#1)G(A>@;Cck%cfG_yp<!D#aW@8A^6347+EZ086yJ>WDGY609&2RRWku~jNZ0wJeM-Fa?pHrYo>wUL>UsfQ_Yqyr?2(f(hck~H1TVWe##GgH_PwhuI6S3Muv2myvgY=c4q4VK{i)K7b-_R*@O~ROcBCK$~Wm4Q0yG9gLh4}WdSOZ*igfwWSbYTEVJ-Yh1vXqN_rsFK|To~t`JG=}cjwzn~0_?vw1KXb%uUg>B>TJonPw^*2{|(iIX8mmAM<lGb_vcdT3BeV@sT)yQ-z*4T^Alw$H-DDU@_TBW(np|Du6-QFw*5xzNfHP7HreSJ6rN`ysoS6fV;>uih*TTt%vyoNhjPOu3e=|$SYWRsv9kj0o^LiQAT2KV*P%e6sT5Vnd7uddi5ZNJOjE|lt*svbjKxQbnMia5;GOguh@U#+bhYTAuG)0F?gbdd8(i{G13%Rrb8)(7$#%NrQN$dNyswsO|DcD!1XN_6AC-+W7_yPE`-}QR6>`B<<j$vgEWL_bvo6hZ_E(_L2Qgzq2m!Y!WSl!D)}o&^A<h&SAfsW}uDb+JB}Gc}fs5A|pHu09_E+8(C{Ph*v*}r&0+<~>i2L2;(5HrE24#m!>p@iWaqj+J*W`NNG&IGQ+14??E{@L5DN)qPV8StA%P9yBKnt0g;rQr<Q{SG{|4n~qwE1%T%8?O^j<}~CrVjsBn&<G*b25~Hj7hhCu3=t(6X7Wilwo5x?xs=G(eO0Mv-Q!h*Db?ylzCnIY~z!RPwi!Nv@~}*73_#px&4-hi6re=rK*ufTkJ^gjWoiXDe1B$+mNI(qUhm7tcwrOu?rqQ5w$8(Z9u{X=e+D98B!4f&ibVn*Tpf7(Lb)jPu4!-MeD$R|4gl)*3K0veeQ$&wzGASn1T?Yg4EkvV6v0uF(0ywCnrmsG*>UU?q-lOS<#!mR`BQ<lwosSWnl}t+jP@bUGsg-XG#+K#tBJ5ULA&gzN$xUafH!bS_~l|3VR7Dd#ZR#lulwF3z9?3n2l<|;eGw$(p`!`h`j}e5nx3n&6Pp1s+w(F`be^#1UvrTn7279h4sA2x@X2z=;2XanTJt|I7PlGy_W6Qh?=K6G0EnZ+<WRfGN;e~D<szc9PESP0Tp|!C*eQMyz_^$#$Nx2`!eZ_@>7s16iK>6dGa;J#WQ5UWy$UMU`ttwcQkgiQ-%N%0#*;4#`&B+-BG*r*52ldL#Z&B?G70fy$I-GP0?NYpz^Osi2A$Rb^KYy#OcOvhL2yMphMB7#JQ$$1^D#*z6^F`%djJN<qm++^L-z0?W#O+6h%d)I8neAtW?8GL<CTz-lH}9p0+fQ=#riIv||_R!mJOK)l4Aws~+$e4QY<MtAClJeC>|!A^F0Q2i7cTr5kM8(5z4XCy?r>4&!agAw&z%HJa5Z*j6nbLqQK0_Wz<VINX`I1oILJ;^d^8a4k$k6MmxWcGqcA0+PvXNQ&5Y$KmHOK%ke`is-yx!FYr<A<rfT`Dnp|R=&$U)b5Z3NaZjY!e>tX(B!ork}%>$PKf_*m%nXGXnPeHoB00jVq(%-K;0f4me7X<vXO+al*0`^=GPN`TJs6Ay7_w!e(Mv_cuFo0Y~+ZAE%8*;9Vq72UkVO;ZFM9la0b8R7GvHA>$vbM078kXcS;Wq^<3DIUVURcRD+%eiIc{w%W~<^MQlV~3=mChF$uMW9l5wJ8iLxp5cNpsKWg7Of5g^if2U;d{k&=73@;3NX_pxhvv{tqw_c)T`&TtknM8c%{dahQqCSfDxb9~aGoK@CNEy(>Zr}xccpTDugv#$np7Bd}OxoMWITIN3!(0bq7l7{XGX4J_4+UMufbeO+858cSOiy%P%EFa&bP|jiWcA;n0SgU5jY2d>p+~$kTN!Lw=hZo*BOJpQ>kH(tu(XpO+>la4I;^kcrQTBVY-Y08ami#N@GA-td`!<3*_;aF0mxt1bi*zFuM52kWOaZ2EdVBct3E(`Ks$)5dvS<tr%5ZG7^u)`h~%aibd#p6R*<Pyz8QPz#iI!VASyR1x`B2eSMTD1B=SqOFc$mUTF1l8#E^|%u5{;uobHNwt|@Q6&I&BSV}?d$BKI3^`O%(k?f9;pqGTu&0l%+1OZA29ym>eV><Vj4kW%I+*K!EFhGsV2;LT`=9%c4sGwW|gp)c9gWK1mNw7H1A%J-P>dsz=}Z>wdH)G&pQh>1x%b5QTk{-zL0G&Ju=2#i5_!g+rMLd$2^o*P;(fI_4GJNfi5u^u;m;_2b>*|r$RRU?riKc^w^&)6!-I*&ECycLmK)!sZnK5f;dhtS;tz`y2nZ}zVZkfwZ`+LoEWl|9xHa8cd)#s%xP471~735<4r+wq&0TBHsh)E>CYR()8yKcJ&=>uSlhXNibShR+Wc5>BJ%+mN0=_&zRaAe3OY?0u`f`q;G}B|Z%NfYFZX8x_RfSv-%j7wgk0B-IOgl`}B(57PrmF?IAHKjaBBV;U-B3U-W;KNTwB*qu@wk#bfXAMl$Xio$xGYQ<TINhmY8DMF|1NYv4{rpOF^cNROWF#N{su3d>Dl~k=TPKy6VA8220nbOUw4~KpG`_3n`i1EDIs-Cxtfc~|Hl>Ehs+2usECWG~0AA8&RanE3SH2iq6?uvURWGI2?4Bh`0^xRcPlb9vlV>|Y;ellxTo|n0Y;DjF;3jl+nyd}qXs+^;#z~Ai}p@t7OP)&V6yPMo>ZGR73QDN>X+E8;F&KdXIV08#ucIRc(nF@WRqU|a<fIl-t7R}?Y<@a6|(m=s@{bn|yD35;9t2s)N6ju+?{SExZSN5M*HaRc_lvl+3hpoHKXP*GNWme<fCQ^OhUddlr<lPFijE#BxGO>f9-3uvtwDEcyI+swiKWa!h$BnMIxlW5uZV4<(ZqaJ?SV4muj3|%l1cGTGiz$tLu}$F_EuH)IMyf;rKLWqt-mf|0{T||4y1>|d(KSc*Bf-W77|u9~6&*NQc=iTO95qzrwJ_KdvLTOcHl*PW!E#G}_8H_g)<~9c{{%woVkZk+Bs!^y2h1A=7Cv<lxd;tw#fu5$p^FDzeD--HSPnu!7Yyv>h|cDRwXRGP))EKpbfi2O<^<_=I~kZX0oqbvcDGtm;M*z1;TXegoh?P2RF;VIPO_Yr2%mtwk8JI6OZE5C<+P_2X-aVG2X=^rictXTUgcBXb42@sPiq<-F~fdr<q91eoXg>}JwrP_KATBbSy4oqwlq-C1x&qgp@A<x8TjMcPbbG|uK&BE2`~Z>5Qc-VTbO$J;jDh~7J1lrQFw<I2r_7aDz79-8x?IG99M$wlufh3PWeWbb{Md*Rtm&m@ualZb!+lM;EC4ls%RNp<r()ycWWI-fg!2C#wi2wX`oP(z&N7z&N=i$_Ghks&v;d_7vtN}rt4?-1QGEJn17l^cbB{-2}n)@{|r_S1&8u!@6g{Wodr0Kjv(1Zz&T^(|M*qGxcpKEqi|Zj;sG1fV)Q8?sO!VManbj+X&e>13&|&^`cEN&AQte!PsKrLCok0OMCwAN)%$YuqRuA{9s{*cPjtm2Z=t|ZXRdjDOI~iqX~}oQV07A9Qr<6hBzSfY8=-?X01z>Nbj>f$DXj2}KfdfYOZ-t!9Sta7XQ=hDk0Z>xOs?UqBkyt^O{!F(9AeeIp4vmCl&9zldPqN=2er<T=l-%gXq7;oUY+GA3}9~Rbc8iaQo3HAwOr$e$^yi769ZSEW>AagP|mFla*#5EC4zX(>Up+`^}>OMb^q|)Eor9>_r-qC$Z8pl2#gybk1;<LCiCoK#+}<?zrPEJ5Iy)#l**H?=g&o|$#D9pIR(b0P87x7PG%;u=VvAEnosTeXTe1tkT3hQQ-SBcf)2k@pP`4CUKlrFG~idGX%(9$@3b+^!2;X4!&0lu@$9(F#fDHZ#ADX`5ImWzY&59rzFH|tQLUP+ywU&c_kI~Wo5R(jM3o2vhbluf`m!+-$1qS?9`+x&kI;6F!_uk16E3^4BTU=su;6_HTy7r2M8C>KhrpDqlKDs1*%t63a=!I^?W0AVYpD>yaL>{8*aXlp;*#d|SR{0m$&CcvZ6DJTL>1OOiTYdkHB+zZzkCywAsZxqEg-tpxhQMvN>QWTls^N_gi?avv&&D(phYLuq^$&vgf@5~??<-CjC4xC)@DD`xJc5Zj`0Kyq%-~U$-JL2t=@l&(q*BHgAIV2-s8@%>hF-uZXqab4}SFMLl>%6fugQMQ@5p>8{=g-hx|78tDQh~#B&teZ`4J^esU&Y;DS*3h9TRztp+Q;L;yfg<EC9C$QM*kA)EnDe2Vft=04g9V+EXRbK_0ML$2miCZz4a3TzF(we*L?>%epRbO2189IhvN;w<rrLG8cc;mt$gA{-7xs))LFJ%~Ou*l9Dj(9aZU#L61zZh_w}DV&50lrw@%mDt}|W5Z4?fz;x%;0~bqLSctbhAse&o`w!Eq5<PK(KayuN~Dn0j!7LU6}Dy)78}Bx6)sc%wr#ocqX4V86XK54@9PK#dn)?iOzk-_xY1{TIe`Uuc4DH?29}{UUqS*|XlSn`Jr`M5emkt5M}jVS3Cd~Ve=K>NapVw+Q&7s+dI*|Hymj{oqJX#Vnw`vB*jp}77$Mpwex9`uPll{*afU`2D4C40FC5bSzK#0`&TS<^*j6>Z@PThgE5R<HI|5U~KwF?hDf3jJkYE_c-h$)iOLCc>cDYuzI-2Z-It5x2@;pJYz61e(+)P!;cc%GP^QDw)%pLPG==(hW{Nd*lW7=ZAp~_=JalwR)L#Jz6!MnqEm9AA>=w0>P=>d$xU(^G<r^K2?Q!mF~q?<8VSd))lou~s4i6uAPZ8ISw6`&e_K6s_aHBg&@lBo<_L~F9i>>vWBhkX=cC@`rQgdyn;xp&6gaMRpa-diB9nzcBf03lSq@fP}74No7ZrC*l7C`0aPyk;Orzl7f)4M9~*hEBeBt|8uQLT=%NI;*B+%8)u8Gb3Hm2$yVeG_g$B_m@h$tPB>Fr}mSS1c{^e66=P(^;XNT_wdl^+xIaRm8R|HRpc+*rAg;JW_hmBVGe9v(1xgUWh+H{HV)+DQHQ%o!fC#9Ml@syxaH<ucmNg5IV{+ODxggdG?edYRT+3EQhm%9M_P5_eB@?3IRQZrielnNt%_<y!Ld6E?(&C#Do6z*zvaNV9B?y!=s0R0l#HT?vs&@!Vhu6#ZkD5sJ|5JpNLtbh?Bgb(pbK<%w!%R}-?zqu1pxT12>I^sbY1m9>4=B}rk>|EHEj?G^b#&0xg996p~trvt5t-fv*u8|gR!3nkcMfR)fdd2XBRQ=rn*?q<p%xqN?WON?~+Aiths*V0|OZ9v1ef-tgkS4`4&iiUhT=LsoO=3fNM3(n?RI~xW;(L1C-!ywzi5xh^07v+{auLmaBkI$bak3R)n-}5*C}SdHZ>~jM~o%W}xb|*;b8geBMaCSJc<=li}Xh{{QcofE5qKQT;sc{+uz`myKTImXyF?!o1OE;;|J&(3uayi!tvW2N!j2BL+5p<t|clCwLu}nN=p#coZEodsXf@Z`=F*MWflZU->8e!D<fr6U64`BxF~6alg5sG)h7|5DK-o-~dz8iAqSrORg_3lRMEL;u6^LOj!&ZG(<7urx1j=syyovTghj74u8r0L9cuddeqK^KF(`S)1Ka{Mn!6%2JsrB467!VyYg2B?|_L3O?(6CVJXX}{}i65-;h6=>&vSwZ6=S*Dm=P+(q{rq_b;TC<-R<({Pc2e9VR(ar6B<-^w5u}ALKstm147ClYIwSzjaF4ZA5}b)3N7>b?%gVW=fbyD9)%A`H8Y^&~x=DZ{S9=Iit%D_=x5{%RJ)rm-xOWt)KIuY?9;HzkTL}&czbkOh-{<U@tD2q1D0-StSPt%2tlbnqsR_szIRAVaPofiXI`Ue+D_ZEF}R%Kxz`**T7=#a_m)_ot~1;daYx?E!Og@ssaJ|7STY95WnR2WvL;l{iWE1@Sj<6{(F^EBH`Zch{6ZOr66H6j?)pEZ#i6j7#|pG57Y0nJ#)EWG@a(a()p%XcZr4N{iVzIeK{3VI%FM1p^4d;JKiiA*jWMKw7x@)Qgg{^j{gHS=%cwtyKBfM$WYdXQ;F@Mj{YmtnH&V`|LK5=sVdl{6h#Kgqa*Vt=e_Kb&X#Y;9b58dMN|(WL-AtBXh8rI`e1?**YsXE`)?i;OD=y*<5630RD^#OXGdmDBVXvOq=#W^BGz>_0KheK>#-Fdft1{Q!+sMyYgHlGz&CNhQ^vEpln+)9;%iCHdWn!LC+*9Wp!;x)7zIQLnWNU&2Hd-=#({u9G&h7CO^~i$(7gGYrw9`cKPZLH#!gle5Wz&Bq~e!j(DayOrtXd2G%pe0KgX&D^Tz5xn!mIK4-Ft$=CFhhap^&9sz3y>2SVnwDA-`zu^4&wd3EGnn>g0k-xh6GF0Cy+e8xg_CJySR*|GOvgIT47)x0(PBNse@o&1dilcRk6x5w^&Sd#^&muh#1%uk}b%edTxDCNdb*A9me2f1t8tPDr0vz@b0zo6k}>B~)xQb<QIdNGGrH|~<m=*hLZqa}-UjyPf!p{<#BCQj)N8$ill;_OnIOwrlZZ~LAtw_Z4k>wD}Wf5GqB;rJxpe6)-rN<+qaW;sK_7N^XC6pG31&!@(Dp#7ji15o{lT^p5qrv#S-lv4aW@D5kZrruIEims~={~dH9$%+XrdJ!CC@D>bcd{(C!Q-$C-o(&7ccb8Zv=?(;jHq(vb^iAj#cGM=TMy)KAu4(dZ#zyd9F2*#<;E|>PE_4LVh>^%qXbI8B$pk*4mx(K;dwvV(H8z*+zFyjFiRz-J#o(&k=d=ojPHE9_$+i^-yEBcR&%fero^x|hU(S+F<GH`DQl_{VB@?Rm|EEF?7U)C`r4BH#p)A<r)~tqE=1volWZ9yrz*{9+xeRy&(M!+G*Tzp!Z2Z&r3O+$W=7Vd;Y~%=JVqn!cNVp^{^gCYa#*aq+<=sDhA{Y;WMDQxJqF8hG(<(-s2q*cz8K4EYX57+Ze}~aVc|=A1qCL*WWYMX!*j94XEp>=Ag8gB<yR~AueA7vU-P%20P!nI#vX{=;Xm2Yz@&OV}VIhA#ao#5mA7jVGt!OZ}oPR$_h!B?M=0&}}+WlQ2B)=L;S$L*s6KP9a8+xlHeWP^nDwAb#^@zrrS+0z%u8KJ_L3pMbXeD-|^`!|I26>kJhtc(%r<(&3??97t*0E}dS(l|7*6oiUS#CgwAS1Dj?FIVXO~_)GQDmvabAk^W=(xf|>Dno6xKa3fD-$?Z^cCOmrzon3Fh{6m<rfGQW}U&7uZcFuvB1)Z)lf5up<(#?I&!D#TE%iO{fJ}iaX9h>v23=#cn`T9HxPJknC$LI6lvm!8s&N+syp+3L@~@S8C77AGZ{Fol`fzI'

COARSE = (212,212,(110,254,4),(119,269,4),(136,335,4),546,669,-14600)
FINE = (306,306,(159,367,5),(171,389,5),(196,484,5),789,966,-14600)

def coarse_seed():
    energy=open_module('energy.py')
    main=symmetric_axis(*COARSE[2]);contrast=symmetric_axis(*COARSE[4])
    m,p=energy['decode_coarse'](main,contrast,Q,6,8,energy['COARSE_KNOTS'],energy['COARSE_FIXES'])
    need(sha256(array('H',m+p).tobytes()).hexdigest()==COARSE_SHA256,'coarse witness digest')
    return main,contrast,(m,p)

def prefix_evidence():
    conv=Transfer(*COARSE)
    main,contrast,(m,p)=coarse_seed()
    n,u=conv.advance(*collar_start(conv,COLLAR_DEPTH),COARSE_SWEEPS)
    need(len(m)==conv.msize and all(x<=y and z>=w for x,y,z,w in zip(m,n,p,u)),'coarse witness direction')
    count=conv.msize
    del conv,m,p,n,u
    moments,gain,energy=(open_module(name) for name in ('moments.py','gain.py','energy.py'))
    need(energy['THRESHOLD']==gain['T']==F(-73,200),'threshold agreement')
    rows,initial_energy=moments['compute']()
    cost=energy['energy'](rows,initial_energy,moments,gain)
    need(cost['energy_upper_1e20']<=CAP_NUMERATOR,'energy ceiling')
    majority=energy['initial'](rows,moments,gain)
    e=semantic(0)
    need(e['checks']==94706 and e['coefficient']==SCORE,'first entropy partition')
    return {'coarse_sha256':COARSE_SHA256,'prefix_depth':COARSE_SWEEPS,'collar_depth':COLLAR_DEPTH,
            'certified_pairs':count,'coefficient':SCORE,'directed_energy':cost,
            'directed_initial_majority':majority,'semantic':[e],'required_followup':{'phase':'terminal'}}


def graph_terminal():
    energy=open_module('energy.py')
    conv=Transfer(*FINE)
    main,contrast,seed=coarse_seed()
    start=energy['lift_law'](main,contrast,COARSE[0],conv.main,conv.y_axis,FINE[0],seed,Q)
    if FINE_SWEEPS:
        start=conv.advance(*start,FINE_SWEEPS)
    m,p=energy['decode_fine'](conv.main,conv.y_axis,Q,FINE_STEP,FINE_OFFSET,
                              FINE_KNOTS,FINE_FIX_HEAD+energy['FINE_FIX_TAIL'])
    need(sha256(array('H',m+p).tobytes()).hexdigest()==FINE_SHA256,'fine witness digest')
    moved=path_verify(conv,start,(m,p))
    truth_controls()
    literal_minimum_controls()
    upper=root(conv,m,p,GRID_STEP,GRID_STEP)
    need(upper==ROOT_UPPER and 2*upper<Q,'terminal majority')
    depth=PATH_LAYERS+COARSE_SWEEPS*len(contrast)+FINE_SWEEPS*len(conv.y_axis)+COLLAR_DEPTH
    return {'terminal_upper_num':upper,'terminal_upper_den':Q,'coarse_sha256':COARSE_SHA256,
            'prefix_sha256':FINE_SHA256,'prefix_depth':COARSE_SWEEPS,'collar_depth':COLLAR_DEPTH,
            'depth':depth,'suffix_sweeps':FINE_SWEEPS,'continuation_pairs':moved,
            'root_truth_controls':72,'literal_root_minima':162,'literal_NAE_controls':9}

def terminal_evidence():
    graph=graph_terminal()
    e=semantic(1)
    need(e['checks']==46081 and e['coefficient']==SCORE,'second entropy partition')
    return dict(graph,threshold='-73/200',graph_cap=str(F(CAP_NUMERATOR,10**20)),
                coefficient=SCORE,semantic=[e])

def run(request):
    if request == {'phase': 'prefix'}:
        return prefix_evidence()
    if request == {'phase': 'terminal'}:
        return terminal_evidence()
    raise ValueError('expected the declared prefix or terminal phase')
if __name__ == '__main__':
    import json
    print(json.dumps(run(json.loads(sys.stdin.read())), sort_keys=True))
