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
  def grid_index(axis):
   width=len(axis);end=ids[width*len(Y)]
   rows=tuple(v*width for v in y_lookup)
   def ix(ip,jval):
    if ip==width or jval>y_hi:return end
    if jval<=y_lo:return ids[ip]
    return ids[rows[jval-y_lo]+ip]
   return ix
  mi,ai=map(grid_index,(M,A))
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
     od=[mi(ip,j+shift) for ip,shift in locations]
     plan.append([itemgetter(*od[r::4]) for r in range(4)])
    else:
     straight_y=s2*j
     aux_y=o2*j
     od=[mi(ip,ceil(shift+straight_y,sden)) for ip,shift in ray_locations]
     pd=[ai(ip,ceil(shift+aux_y,oden)) for ip,shift in aux_locations]
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
    return open_module('gain.py')['semantic'](part,CAP_NUMERATOR)


PAIR_EARLY=F(34007001683,10**15)
PAIR_LATE=F(13984228719,10**15)
CAP_NUMERATOR = 12246713453527336813
SCORE = '18835317372701767986746555389/3674014036058201043900000000'
COARSE_SHA256 = 'be5f8509a7a1b7eb803e3a16bb59ad4cb4f08ebaa203c29ac67b9652764c6698'
FINE_SHA256 = 'c925090a13fdfb103c7149614cf1348085b373df780031685413614a48577835'
COARSE_SWEEPS = 20
FINE_SWEEPS = 2
COLLAR_DEPTH = 5
PATH_DEN = 441 << (58+LQ)
PATH_SPLITS = 9
PATH_LAYERS = 3*PATH_DEN*Q*8**PATH_SPLITS
ROOT_UPPER = 8161
FINE_STEP = 8
FINE_OFFSET = 2
GRID_STEP = 480
FINE_SNAKE = True

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

FINE_KNOTS='{Wp48S^xk9=GL@E0stWa8~^|S5YJf5;6@%4m0bX26#f9ZGc%T_d|x+{MBOylXSr@V`)ocu?N*ViU8c+7Y+#Je$3ZzFw5Df&1noUm0|1dF@|6iH<fwkYNb5PsJOK?C-Y}QlbxaR?hLTBl1HrA8oC^21+ftKZp0ET&ND^?jQ9qY<o3Oi<&4Z~0fq#Il<%bpIy;+~jnw^Co-fn%@Y;F-O-*G`H2$=hWpY3&n^#c$(LmcUvtWTcaLA8o=&IwRwcP|id@^BT_JaBN6MOD2dyKeT~FAsXj3~s2hDoAKRGkfCvzG3q&h9RJ&GOQPuwhf6W4!)jkP?E)|#ZNs2)AC+F%334l-q@V!D6%GB4OBbA-Hu00;S8tT5=zd1RtXnI?Vm~bMMbYaUx=-O2J>EUbF~p>G@^#4jsB8WkG&guc2zWB!q=nQ1lU%)w;>UCOAMzh4@UXZ`xwpA0@2DfN0|Cf$cb!7LF?5nB)4PwsS_P>?tydNnFpRjjH>Tagc+sCEXH#QG$1}m`n^P#>iP3;yPnw>jSaRH^%wMH^{<32O2Wv&h&864@@%3LDe?N;_c?f2UA#_bgPUVV&y+@(<=RFtZ#ga`^6tM0K959oU_bHMjhq!sAqu7@k)>WJ!$h*aaQ(L}AT*+vBKJj$;Uj;t2N3i10eUizmmZ;XoMj&ACpgV_A;fp8Mjl>RinsaL4Lgrge5cL%f-0e2{YobWVyp{#*hB8B`fK+i`Qt5z5N|;k*TPyV_>W+5@nzrXbdK43d12b&`R_UO+5xS_7r>0@J`=%mT|o}5y05L**JP1nl2jMT)WbG`I7siuoc{gqKnu0Ts=}w{C=pRnQeK=3Sg_?Wb;D7gN7A_l_`R6#b+biX@zxFF=8Aylxry>g#=6XQB?<<<J*-brg)1g(Y=Z8<M_(IOOo27Svz1J;RWNJ-j~`@Ybg!y)RK&kD3Xp`TTkS%VTB%kiRr!-M;Nx^=KG=HvI(t7`h=!1$ann55&mGBag#-h#*im#eYoCrgi-~GKk<{j_y?x(lZ4uj1;XtPbte9>=O8Jn)si&*pa^?bm%N`ozXS`PAHPsVH@up~82EnRYxHg49paey%t7@%JXlG&3|C3a#?hIs&UT}2~mtl^t|78KA-^6P_w-}}oFg`nxmpk0-C)zCQCzOn>;)z-Joljf*{Gg{>%ibABn@Em%Z(oxZpOlz9nDmmG3`u7)CZOn7GOFi2pW%-Hd!I>`_cOU9rz3mF4u2-|QwpZ@Z)ZF)Q(Gl$jg6_(SzvDWfN{U>)L=P36nadtjDNJriuFSk@9ndKYYrZG9<~-RQ$BdL!dsT*Cc1iOUgrFLoe8c_VYzoCmbTi>TFZPW=m`sYkbdrGDxg><dFzY)g!S;G6kEuX^uphdwCs%|OW`aF4_S><@yqTc>Vm=q?8>G6ihuJO-myqsYQjjds^JWG&e415W3kJY`hp+S@A&!H6T}+zJ$RdWlh}y(;&2vvg7G_aZ@P*Cw0G#~4-1@y)oaO~bgP0W1D7#eFQG`9rTr)zZmAhq!jDXW<$fAD_i&6O1QUXpd8pJGD^HRWbnTVYzzAqeLO?Ix71FWGZj0NAUycyG$fvqdjQG@#iWn^(b<Uv51DxpTs{_7o#@z)_FeyCPK$XSm$|&$19>;a$^}HbXL4moNwSSAa?&On+Oi-dj8%r-d!Xfz^{xt}qn_|}wtB)0hQuiJmm}lFTfgIeEsp{QzM^(B<a&`UBpHWh`Vq<ELV{N+)_b;RNr)St%?(%Xx>Q4`j+5)b?L0FKSN~#sIGHfkaW{4au@ByPCv$${^x?Ui>+2)0g!gS)u=XNwYO?y&VjGEV6t0t=0DLidE#f&XR3xs~V%pQ^0eaIywijRTGE_tgBgY#l@y(w5Vo+CN{@qPhmqx9PVB%7gngIVnb-CfA5GiVKqVj!N_U;gRY1a`0wwj&y_b=Lr~>Uw{c`#9BGtiOl>?*V!ViFj@i!z8;>w)W4j>h;)3*`gQ3lV&|>MPh}106q)-0m(*&&g|{QG{BdY+;N;HnHN`et_+)hf4^rHjG05brW=^3C$C&X?{vwtw_-izmj5u3@Xg`N<=LbY>%22}D-72x9#8%Wult?KE;Ut69M;Jp`4Ggkuan2pN%65d_-4Hs#<vCgeIO!yYbG&2b<Wn*4<8GiBAE}z-_Qb*$O&a#qqKzwjR-luMDEdBX+d{>v7o;{7fI=a_si#K#q{7io?zJabXk|XLELX(%zjVTxb(c$EBhLtLPdhybDsA1O#+!TQTCrCsyvgJ`B(-A;P0x<jb5oIb59^M)x9ypf&MiT=_LtKOdofDF7U2Ok(VQ$mNh-!9f<k8?MfKyH&Q8e@$v@7a-^luNf+gk)$1!WVvP}}l-Kf4xZt`dA+u<a7CqqNB$Qj%%opyJCt^b3JBbN1YlX{8!#$bF6rSN<Z@*q>1!OGb2o7#fYhkoZ0>~F1kUYh~7w))^0Aqxjp(~Lz=<@^?POV=*P1SXX1|UE|Cf4?N-5OX($QRl$RiaEc$GS7lAr{5oJU$@;A7JRfi~HR2gOKs?M^N9uNngkfDf%$)zK*CIV5)0aXc2QS;CJaNHduW=J6PNYBct6uV|Va2sdD8nC#1povq{u${xD&><J&OXhX<{to={Qnzx?DvP&D85ggFH704eIJw_@s}d(5WH+k?LG7}i)dy$Wv9aU$iws{7rrNMh5RUOCa7c2`cejDEE@R+uV#)(7m&S-tu(%9{fPd_>Z#qJejL1q?1ZBjrQcTJ4VZgNUyVC-H~>wkI24HwGQrb)fV?FHKjdT?gWno#>~F)LvMdcCO{))+8j0uM&KHM`e0HxDINubcVV;uuRG{JADzBi>GY!Z|43C3tcPTgeP`g4IS0@mT~j~5#YcSQi8oo9$H_upP}t*$$dZwu?-GB2|d<A;bsO~yo3ap?LoxSLJ*UPb8aWhh@c;ORlIP3Uho#Yc-&jO1JqqYaH<abkjR*F)kor(D;CKl7))41NW19cl~%83HI+XEL;9z6oha{-_0ci{*y0byVqbuDyF@MvE3QW@Hu=`_W!)Xz4@$1+u&nzNa)2G?&FzCNB-`K5?3irb6$^EoSD3`Kwa#ol0Ea)`B2|FbGIj4oS#mE)C6|~!qwZhMVzT@vWoToiZ@ct!Qv)(O;jXF5Y~LKg2Oi-WM;|o!0V1s5ZRyE2@a<BQ_l7-rof>x&<tJ_of_(50xAXsTulqKegm%q804mCT%KQHup=NI%|Au6}J5q`I|8w07KfUd;&$#Xw#l{*rzPqU@BrvWD_59?&jZ^HA$(226EY~<ev_FbKwA*$+bi7<%%glQV+`{nUvz@P&L1!;6lgsuPg&s0mb)bU9W_n`iA4b>OAGuBOfY{NzHYRX`HFhlkCmV%r_ttg3bRsYTt`oj70Rf~CTL<3KHUM!FY|^S?bps)4ju3p=lZv66nlIq2e!w)B1fKA0cbrwOFKaSK3J)R$*gX#^P};y9adb2UbfQ8h&!nYXaq{FY0QbWa*FrJiaRR1ii|5$hU|XC3Ylb{kyIl-KwDJZ?_KJ*4U!n_^Lr1>fzM5<uwP6;Qg;=npN21T=7XSuY&PUJ^r-x()SXD;X{Hjyb!tzY{;a75{nlvm=P^%l)0QSz>Xn@}D>@oM?$oLG_5Y3nirb4tjJN+EP9ECwPI$o0pg4|?dyjF(49iU?yjU(gP@f-mjgUOYR8vyxt<B7!D{>rQ=?x@~CPf6?FUL!(A4QO9d{Kru#J%aE8)zF1&+v!Nu2Dl*ryJ)h<s5MU?)G~qI-l?XD_3clYM8LQEMwu|+*`!`ZR546I9R)+$CP8en_}_QZb$9GBEDaNCQqAY=Pg9Ai3LHwnQ7))c>&!i6=|tnVEA@GoQ`(2N4**lWZJ`CE^kSckcrEVsI>1AOYblJJsE^!CVP8DHr)z_J#e$dMy+qCfUEHv?)d82;teHnnhK)B62|y8y+a?*o09o|c+^mmsRl@m%7*g{L5dFbN?7<u#1(?j;=<W<@3I~>ORJLW<uIn&B9?Fq+NV9{Xx<OT{Lk)+`fe{Z%?qS;TO-dHWwN(_+s9;3NOGs^Gx_=Nk5L;$$v*@3_rj!_#MROJWs~X;49Plcr>&uzzL3<C>|C3%6DdJaa50=DlRz;-sIzQG^ms6J5nes0Q_`bp+dF~d}Tw6((mP(a6aD%q1_73n{uL#e{@3YDo+Wr>_>ie8WO`;us)hD`t5)RO>6=$_8^emGdVp8-$Vlno&qh+4fQD>*uq<J^rB8w`5JvlbIA1f?Fw@^w8*S5?s#XOL60FMN{c09H$igp3H=2vx|q_O(!e=a|w0W>$1-29^?U$;k&0^@5A6_jp-2d#HGXsVxocKT0ceB{i0!m>I8>>P%Nd_fn}=_Jg8(8Ca=6tcEF)sQ2p4einpawTdY2`&a+xX&PoCi%D9TQx1s>Y8~mIk0u(QP0b6GV|=Q0NfFiqXyw8E8XHUQXHKrrL^#VsJ1pP&bIfiqB{fs3WZVBo)|D)E19q{*1x^Tlc=P2(P3nkOzNN|6F#h`jFBCW$Xkmz1bW*(j_Cs}bqXDhQntwUrN0n<4@8O5dPsc;)B}HwieiM8vl{%PCNEk@wpt&zBY3Nnq>(}PYP{)RRm4Bmzrl^9y6?gOk~vbcQz?HymN__c@fJld0AhY}iY&hF+YaUuf6Aqbxp7D@1;au4G(O-}y-@yHaLZ1tu`cku_hIH8w5jQZ0F~h{!DLY+HOz32H8xmLDrxHN1xm9JBr4-7%$XG1lX*gjz(qmhQ{7Sy7sIffNzH5diCTpc4-9MyWQ}}@DVbdoz~3M*CA{hZOLNxtq(lZS;6R0dBK$El=viiq{6h<FDOO}84NNpd<OS+fDzWHP&|?rvVoG%zT`1&2PzxBV6eW#$Bu|*-4iH2RWZS<Xhdq4AD;1>3cMu?bDWrG8frOi`K#P~M)eN{u#Bb*(z`QK$Q3gNht%u$Vso4eKd2gxRxdrXf7NY^Ti8(FV$2!XXVR$u)(J_}n+;srlmf((wA_X|W-lBMfDtY+KC_(ecHE+~>=Iq>U^|U$XLCVGE_#i|VD_Y`flbw5j5`^WPAa_IZKq{S@QaiBp7%J3>3;rJ>%z^L?;qQXk`CeUzEtpgZ3OH6ddP?$j|13R7V_ZVn>0?KoN?hcH%)oU>9U(X9x(DF4lfpP6_!nlQ$qMeh0vp>{GU-VlA3IrJuz^0K0=s?&RVAlgy|*r-moxjC)B|*@9Nd{Fe?P>}UBhM^H1SQe{v{oBU@So3=m(>t$nStKshQ`Ccn-)vL)?%^&$f3D6yCF-y~`{HsTLGEv=uH$8!5(m!4JaSL1pp7N|>*{`>ZiSR1mf5bbpjwHDW1qr{`009JOkd508a{hM~uTP<sm4Qhsa!!!9)DtdSh?afGsA)3L6nbV>=*b5XzK1lfQ6b$f#f=aC5R)3h*#Mw49Zi!b@N{V?An-_NVJb;DQ9{yvx#{8}1@jTyVR%>z|3Bbn7*j?|?X`;`4UiEe<BN<P&r;)NJ;2)(Ue0eqR3%ch^b`-~DM%@}2h&08@t->NogI;dzLYd*m*_@EI;hu>~iF{*uBDvmW8H)hSED6f7?gmD_zmf_kgNh!O6x`4`UfV-BfffP8VL&t&urL({vq@YkcLby6H)zoS=Vu10^D-k0qJQUhnIWoYU&XSGgl~(`;ValUBaMLK&NOl|1tG^{OO7!_K-A#ma4YLW<Y?k1xhm`+Uy&|o^llWZbM~Ck&ME7Rz<b?eq#Pu%GEzh2-9<198E`AE3k!BKi*^knjYFaDCYq*wdmGZ~M(u~McQA*!00?-+>eVvQ#&Z>o*psC|hs1s}S14jD=4r27QdaU*sV-muINSvzQpwx*98-$4b0t{D>b$s?yqd%Q1+Ipx9h|W-oFD`j3*EPo(r%;}^-t{>A8&XWmpLfopdh^}o;Qt^tkO()MJY?ck1_?oqZ0xleU65h@0WuP3g*V8W4W4GaXJOc~cAmh*dXjQzgrO&{5jy2Sm?H~Nf>~!|mamB^mYQC;)EU)1XvoHI@y^H!=k$udy<M|2PpB5pQwq~;{X5Ye6JdLfj4}L6Uu(O^UMV5(NEF4CJI*&NfJluPk6Bd+>8EuaBS@cutLMQyGL6m?jrdtmTjn~q`y-J|M#(q&T|)}>VQlLak}@5dMv_p>;H3Vai<<Aw-I3NjTgt=p1u=V}2zcjkE*O0yP+ii`^Zz8Zlva1vjc|F|7_vba>L|dJ;y@F5G?Unt@kXtL+!1<BGJ^yOVa}?OPfy{|XEe*;hqs?-hf^3!wtxCtG@CWNk`lpjfR#j4pW0h(;ac_jIC;Hf&QZ6RYpX<gD|$nsik+;p$*e7s7bDAEKA$_<0(mNrYNjS7m~Btf-eRo_Qx_=is3H$azxclQb^x1_2PVJ8{Zk8~&T2I|v7+*Z6-bO>*kmRK^WUJ?sKH`h_{nuO)#GF{2yqz6*PF^x^3$DdN7f4n{!3-VbxvIA6UR$J2}un-Z8z2!cz4b){>RVT(f$TmwgyJuwU_O%brZ{^czN}Aitc3RsE<YZ?TR!ngrf?xw<Sz}Y%p%jiQc2^GxqXldi^N*s+_RN^nAdq`9pdyEI1~#N)TAWlh_}j>^+Mp9Z26=Bq9ZeaZjL|&>gYsw`H8J3EpUN3;yD2wV1qsmr8ICQinMP;xF2X)xdpBqx(q^4ETq>dN)}|J#VZSBwwdTr{1?P0t_6OAp3T6PKzwc9l7ZjmpMpOpn4tQE$^j%BD4?GP%WTFIboqjTWw}r9tFXmlL#MTbxA%c>MMz_r$&ZX1iKM&00000A(QUvMH>sI00FTnpNs(jFM9>bvBYQl0ssI200dcD'
FINE_FIX_HEAD='=0c5^+yDRpfIFp;qv39?&OKccP+i??*~zO0s>T(SMpN`L?>64d${l9pdq3$H=)&2icp4@a?m!JnA*UMFxdRHBdFR&YbD?;%PK>oXc02N#LsHq9=2@_NwnBD*IcqSGItyE7Jc<y{o*^9h*Ee;m<If_taU}BLTc1C5M-qHfN5gWjd{W`{CR^ee2hC9=y>pndLE>zwnifkvl%Nhj?V`S_j`X2hUeWsQp*d$}=l<B%_HUhVTNI`6#ct+x2~#J3Ga|N%17n}0)-;~XA%jf;REiD-j!<@@8S|MULqEP4Z5x6!!@Pa%YJXXWz!W57V#qWK|NjY9u!_*BV=b$V3N_>8PMMM5=VtR82WfHtUSL+r4_`;}_%_~#;kyF40ph{=<Q?UB==XUEKjMh03Ip@D^Y1xoF`iv<TC*dT-%Wf{d=>5xd&=4)kRQ6y6ntfUb(ZlK;fO#6GyWaAa0Fssfq7<YCKSHE7rg%yhc<K(f^qc}n=fNbNr-m^+&&NyradM?gu8kfP7q-7N>KrAmr@v%+T{<fa_?QZ7%oOq5gm-%J{+Dbp^E)7HX`A;Y~VsQed}2jy1v;oc^V2}P%&29&t5-(^#T6$d`j;AwNwHu26j!8*PbC7P{^vZ4|5mC6*Mz?UD#<pS=yW1C7ltr@d)flez=eig*f;Qbif@74g(Ef?^{TvIJBKs8FoA?p-d#6=YRm*>%3^ov=j@gMwi@$*w#eb#+E}$a0N<ppFF2|y{zm*q5D2Y?(2Yk+^nMxt%ERSn$b^s`pW)FZhcgR?z@2ghgC{iwk&#E`adl`nk8y`1j3@${FVC5AQz9^2C*RQ{GA>d!|T&tI`YE;k+Uj2<|z^)bfWdfq9hplhO)3aEXPnY0B5wTHAP^|wYnnod9>`uyp>W_lmg3`2)P7`Qsm{|MhB@rRX8rZyFsYl3j(&5o`zv}A$u<)F|siz0Fxwks|}=*ZHB-#T?=}7Z<szPY^vUctb^iAZebX$3D+d~V4QY7Eq{N%Sq3FtX!dUw!-3%aoiI4A0cm0Ug+geRZqq5}|A5vrIIA>kpp6?WH8$P9IXFi2X+jt4Eaw;Iqr#j)61S?=W&__NX8S%a^mBVIT&2c#zS+`+I`h8mo;b7eS{_oe-K{t!X+MT!`M92w3ySg#^|!Nue0&n~<|pZSZBer^e4_$vz5et|7Ke+6)#J~9P}U23pi_RCK+)}lQ{zH==*?tusxs5T#f)+l?DFQnv#;$AY=xTi*(p%V(i;MDfwELoTFWtRx2CtxQ)xoVMeA*QX40I8OP+0=wF#V2(*vi$a!enA`&t)Ar2m>NcvqAM)9?T31>2Gq@BWM!^b2Ka19p}-a(d5?A>>}?m~cRdfn;h^`ZBc;8%UV_H)HVuam6yJ9Y&B8Gi<m_aN`PDk=#ova_7Dtk-E7AFuh+4xQ-t?#`khsB$#!iDyKLwhG|yt1JzyIrvOy6q}D*)eIPUuH_N@9{}kj;sQQ)*Bi|<hC~nM6jb~YR>TVq*kRq3DRD4FlRlxJKvq(wY%y;~lOT%xMHuH+Qz7+?`^TI>b31s-li)q5L$G@h0TdwV};9ppZHCmiX8Q50cbjc+@(g>2I<KyFG!+DJ3vI2EShdfs6Lh&)SSq=UVA)h~tiOk3{ymWe}Bp$dLhskJZI-o`I79*`wB0;6Y@=F!6AnP<BP_9!k3uQP)n{}9iyk=6dr0glXkN>-ZRaUS>0EyZ|=_OwB6REg6uA7G}7sjaGPBtwjq>v>X)zU<~s<J<tBJFzCi6gzlAaaaKeF)Gb-4A}TDnbDE+t_f^TuR*-vhki|?W!qBrT-sBygm`3a9fn4%S%qCn_M`dvkZE4aVs2#mJ*CuG~v4O>O>0l1LC*2o125Utej-9gY38U6=CJh%U&#ua7~~5X7VCFSL}rv;`MnRA+uoo5!CT_qYZ0kXKXv<Q?XdrvSRkO$<&C>zp5Z&CzkxsIV-dbDIF~ZBW8oe*d(*wmpgB`gK|?u@YW_tuEHyypGz!KFQ$v%3_Q)V_Ir8Q<T|r0k695v<K6_^IA44VuS&`+`g}fjrWsV`LDWx7=n?r0-*F@l)q^L*9R!vD)@gP0rMj)G?ZG^xa-qg_G(LeSKmE363|tky!ViSE<b|3y?y~cat2Fl}`E1lN*1-W5j!bVZ2uD-cynej2<(c|jL{R+0`*)Q?hn-AvzW3hn=!)>8pY|(|xwqIh7XIex;c``|^6|(mJKT7cPakGrYgYJMJXXNgZACCeRnx6yt;|G??5+I2{^7X%M9S-@fW<$6q@m(|60LIuENX~U1)_;nD)f3Y>x%Idw+=`*TnoX#MFFk$XC61QoRny~$VW6{n;vtjFUgpFs3VB#H#QnKSYY(m6@`y8{L;I&u<-CUK>_xTd&(7oC0X)%(&Oc!t6HLQ7KICZiHKKFwT6grAtEABv+PylIgEiNdaS;5&qRC%zsc)#U=K|0*W*dU{dUYJ7xEb!79cnz_qA9TeHHZT&BUP=*tG+RqRAS1%h*x@MgbJ}^N=jf+{y^`jUe(HQJn|^3m^(x-wL;f(E}5s$%Kl&)j2XzYAIpZS8JrnU0@rYJuZWUmK08h^98TH`Q6A`p4(I!rHYrHXks}JUc+(y2c-y!^0_PRixHnY7>66fZ*m|Eq;+R;F%5OBN#^QHe_%j|Thogj4TM<&`8kG<2{D7isO%F}cE062``+i@c$mR3{>m8-3R%<~#<TR=`D?Bi(hMF^GUSF(dKN(^Ha6qn-wI`pNBzPr0Kc})5qjn~Ootm-=zoN&2adAjtMXD?av&`WJY*l}>k@`*YR42ls|$AV)f9W_YnhG8_p46B{wm%rngq`6c|b|<5d!`q&5@ZhZ{UasZR}~FwQXD19B@H*TW0wKSde#r=+4y-%J(-8Ig%$qUd#I_+e4^scnCW0$0|Y*t3!_MCV!^ecVeJ3byv=D5!@bSg;3vh>vs&BJ)(HaLm7q+Qi_erOvt5Y6TtZFuQJsTVD$N8YE1%!>Va%LLu#kS_zT~ZPL%@O39g`>$)&LKHxX}Z(kE4RHpQM9J%06uds)rh8Y}UsX_EucIslc8&2nzFFfk?gnnRW0cK}hUCN>2<WFPJ&3oB#S`omq%vMJ8bpm{$>ciD5oq@!6X8dZ{)Zb6kIRd`dRM&DK`B`)l&TxGXjAzs+C{vK|OJB5I{*Vlq9mUI>A>oOnx6Z*$X8P#)MAhk^Kf?Hh6NF_Or(Eig<q<jfdR%ZD02{Qt{Q2cgHD;4;M^TdI#bOt;RvPDff!wgt>*lv>|fMLb|67pjyvTXroM5mQn?5ZP7Zq?E&C6twy8qaQvun0@l$ym%_+}nB=ne}nPIs`nKl6<CJw)ituf($j^;@K`HO+VeelcTkDTeK|0J9!mbj>)9P`)^FTIzj=PnbPHdu2EnY>$CAuG%E4#Q<BPb*A<J7p3jXK!jj&zHfLe~>geRLE1f=r`smp#4*+pEllly3uP@;*akzz`u%fu?F6l*tv!p&@xWu_Xy{a;~onqxrCDf;c!W&w3n#8+w3yDp<g_JHwtADRXCpEk4Z{9*_;J9&L9j&mGj`V+IeCrJ#o<2|HJ+tNNx1{v{6;G4zeuzz{yxGU1gfH@iVk>JZ&uB8Bwz*}67gSOa&?^qvaJJ2A_dey^{aYkvs^7E!a8}kwV20(z7kU4ZO~ncG|3jWWw}DRvsPK!x`lT?Y-FYb{{Z_3bHS3gkbY3Adnf7Q9*kOp4_ZaQU%fqr)5|dk$WS$VO7#h5z*#&sbv@ljMxK_vg3{Sh&_{NACYsy0@z<5|jd<d$9-N2L|5(YxG#`mK{Av|f;07OT_zZ?23A?UE#rkH0UM{x{z42Ul*pwqoN;6HJvrh%`sSwGq4M5dxf=Mr>H`C1x!#{-b=d4~wLIP7i!#U<#j2Jwqc=y+%#mz`BuShaJ(g9Zr;u+k*(5L*sgf&q8cVs72Llic1%TR3JkduG1l5QkcrDvK1~Oj(or<%6<xO`1c6h>Jc&2h5~HRxA(j%45tg`nBMze}dWP0}}c%KH^@2E-qf99mn>b@j(D?SG}HFL*|JHGhH$%(T%GaWij}iL42p*L4n}WnS?eU_t{@)dv)zlm0@y$LY_~OVTmWD{ZRNkpX1TQifi!RZH(~&b%BK*eRLKpP?2?79@!Ge*<^1y*Tj*UCa67dl`RL{)5_c{LeF(j-GOaMOB+uDw9T6de~U0V0?kVToNe8JE?m22nj4CKH4?8|Au?<swvcNq)-Q5y`)Zz63aZ5x?rmD}#)-8Zn{!bd^p*&s3t9W@C9D^~NCD9u3c!5Lvai~O|8F7bx0)G@s|&t9(a#8XxsrP%SqXd(n(!ESbfK9tH6G^5<H>OOcNu;nZu{#7c!DtM_+0E#cHR#j$&DVJG+m{)mud)sVYe`uxgC*be+m_1G4&v8uBFrR_hf%#$062EPr1W3vOV|1?Vr`%!xMR7eVL9~gYH)D{17(Mf>xV(zJF>8NAUE?ovB`6KIB0PxL{4Qb2uAeiFG=yRNJuGk+oSsFYMsrJoLy%i8uloWY(8xOtKIvBV2Ut%EE#cB>UuD8<fWKOLsTZGG}%OKiG$072md`j&)i-dT7`8^Rw7G9`I1-<335b3a;y_7)Q7Iq-=M;+^*<%CKi0!^t1X=1}gzf_|j0ut^Q43bSG#m3#%#vUa=($aR}p5Xri09I`^aK9sH1tkIKoFnNk8crt9gU#CXGJct%)DUWVV$(t|%CRl5s!F5Vr9?zJeaabd5w`Fx3dZ-rX!@$*Kt2l-M`srX85|IXPRds{(64GDaNiMbVX*-!0T#liDB^bZr5EB$U(C7V@8IpP|XpPq<k2ylYsJ2Kyp=!_AM$hRnaRf5o32DWzO5Jh`ivT{tHh}KsC7akqAr(A)6W?%avule`@RliF&+<JmoC(B%Z5*aI?`e<BXzIwC^W&*)CaQ+RuQ?>H7h}>JL=O6Q-=N1!0_iuwep6OL=-Oj!q!?Tjs=pGHD3zh|nqG1XZxJ#Me@y|V9>cOTD;p3{RZ8_AmWTx7d^u9JI{_}Bb`1tx?C-G4RPOJd@FCEx{xj7M-l_0eBO-9W#J5odyg%2Yua5Vj9u=zdg{K8})Yi8n#tkqQRQ#d@EISats-`QGicSeL^;fg7>sa=Q2LM8t?%h}`#51Y)F0iEeQz9s{KfY}-8oS0X?R#P;C<hS&^WB#++=P2RE0#(9tF;7*H0-ynX54@s>*0!xKL%>Gw?Te;TY;E^&W(L!JOvt=OvT^1;v^!{xpMgE|D7`=Q3Bf-bH2Yi<P{i7X=}L}4upRsDaqL`zDtKarCt9~xxf)DqugIO_;k4U#E@^e#aV7@=drI~NIjI)jepv08=z#od2KA57UOtcM=*9;q#<-Qjz=EPxf|v$Jj(ni8V!1K4-$0w|gk<fWoDfp;T{`3NMd{_}j@wb&|NW$bWEQe!75<<|HRen|D1_ed@>kV|7=lYpEN8wnnUh7*Miii!dcOqasnnsQ5myIWM`vt6bP<&kVJVxN7QfoeDOEI26JHi>B|m+d>g>~IANnpMln*6|irvkvIX|iXhI3m_3!%GDPI@m{(Iijkl%v`sij(*2Gc3F|%s7gLH%wEoyMf&~hrgm_)Mb3pwdo#MaSK@R8e~g%?#5I<)7p5%@CEnG-%bW_w9G)SuMYrW1V$ka<u%5%k2)~&or_f=+@d>WbwZ^OI_X|&17*@lpVmY5)gc_!&HrFl`@T>3)*&jIZ|oxUh-_4-!J+_{BA}Y|l3RAl%?tz}ZpXMK`2?zGakUC<GRxp~J>X}XeSI$T^JUEk3MO2weZ~}ieRnei^eIkWriYPru_oRUwb!@hSO`OWQlogFKbG&v@g@@KN`U<pOP%_w#?(p!p}krBex30EHcYvQL8yWx38S@C#Fo$@K?{i5lt7bva9=@ai<;G@m0eV+xE}qpM7gpkK?wz5VB5-7$h(+ybr*~%Y2)zD$U$W%&rA8$yMKbWRTiC;P)YprER#+C+R%HQZZOz8|9+yXu<qZ+l*`=`O}Mn9RhVFFP|4=x&38$D(4ep6Fe@^ZnB!yVnHTO^eyI;QY3!j0i3PKjXe2|G6DL+ty|<a~q>CFJX&0yuUUVBiQXUKjMcG8~dd$97cKhvUVzq76M+G#GtD_+=b_(U95z?l#bImPYSG@TZS|v`mOnU>6&$AwcEju{(KCBb7#TD11A6>Bdyc;Sq)&v!qJv04%^v37`C6QW0&wyf3oj7tjt3?T;-3*6VJ$Ov737gK)60VBsnKpo@mgcA+y1vmsJcRu~O7*M6*~d|-v(I>FwsB(WGD4?S{BCygT+TV1yRg}vg0@4orj-jf88aWc0ZHz|kjPif=?2P`k|lP`thq`3(7=TLr&_Rz7p#qFa5-8GtBw(g`vODMtb*sa1|9l$s=)wIjSW{vx}6?}N_PuuKqHbmz*nOkJ21BvngJM07M;Ywhk=AUTZ~=0u#QewJev&dqN=V3%w4Y^F;za)Jg+8?YBQxvYoNB@evZ4kpsntSCGXPy!rKe@Pl$JFI0i9-VZr*>DUW?|R*saP3e_@}t?UC=TNquU_Kn;05keMC72=TL#qr=6bGy5!%7Q8@CfmvMTxGHgKSg{RqPF>^)5;;9RofmbR(uoLTUmbk*m692PRj3g=X4_mElMn|+@nueRcM_fA212S)K#6vKU&A2^=_H%GwmqVo+mfr1n)KiOZ&5mFM7pzY5_TPy@9?Wqb<yZ=1bl@Tr4%SI%O*;@9t*WgNc5-Z>df~zcx~)<H$}?JhwkeFoz0+7zTro(dGCX(mj9_!tptkHI4w`(H8zx_Cq3k@qL`+Dc7kJ6-uu1@Wir-z1&qQK>SWqmq?p5yW*Iv?|1SM;(%&`MmO>kn8xnM02&2!?Bxdcs$-+%RSu@y;3WYzA=b(Y_-PCAVc93Cy~?=lT1hz-shMC(!j_!mAOM5oX42fPsc)o5#szfPZWNHyyb%;?e19Jw8j-c<Yc7{i!bM63`w+Ajxm-$Vkceo<DxTlhW~IDCDDxE}atCk4OPrhI`9HerkQWH!^RD|o>A+3iPvAJJS__j;wW?Ls=wHLJrU%^67MLNwVkV{@tnmC<Y*tNQsK%0<i9er-uXjo9!G7K(!2OF^&qDjU?ZtRrjf=xyA13pH)FOGlxD@H-Jw$Dc*g2wQLk_w`bQ$DhEdABAKjMnLGt#xPP4cC(Hl`sU{Iz%y8+Vx)&cNYMY~K8OUDCdH7ak{rF+9<1%N`73PC;ooY%`oZO-h!($MWGG0y6~S|Ier=!Z0pY^#X?lAsokdjLU?uBacLo2quERNipMSQrGIcIb-KyrSimC<q&paXz+R}UJi{!KCF8yBB^*ebEdFoST_|OlKAkPyE5594wHnX)lT6&u2>Y0-&jA$iC}tbJiP>G%QeZe^Ph@HGRaCNIWw9D@SU52^_7>&{s=;F)4z1Ok!#Ek%!#o1X~Gcb46zU;WwATK^ItjQBgb|Ev@T#`MnQfSlx^lXXUlivXG>fQl>yVEcKTWoY+KQ%uECS}-~rZutAZnc^%-v47(w<01Qq_^x}#-nI(~C5jNa=}cy(_1>j4BhhlZ#Lrz+r@aqHm$EY-gt(VnP6W?cuEW=d~g@t-k*w66%iAaHteLm#Nz>px-fH6xH_r<gK5j$y2)wW_Tw@bevevf>_1FM=;3D$7#1ucAPGq;~(=$#tV_Otqo(E34R+kc?0E{A1%i$#5DGws_K$TYW9AcQdnR5PcO^VK7g&bgdiE9=@Jg(>BP?ZQOzSPx74w5(|6Uk1kN1w40!8xvp-nvi$SwLP!(pW_x93YWt}M94VvibSJeQi8ukUDR`k~inJkU&E$W`iCHSXpdc~|6^TKpHS6K4unCb$Fn*Ov@=ZhZPi(xw^&?gDyRUJpI@^}!12SR>PKgRA<kI>pOJ=yAra#ayYS1tF{}fX@7M{v70zqF!cDW^|czT3ld?Xa)Pz_pc!pDL8x!|e8x)pe!r|004=!@ss;4jw{aWR<VC`deLoVxqEC*@w{ym=dmVwWr3ujPaWS+EzO%TtSh$R%LhAIy@c?xr``;G@nKIG_VNlav?(IVgd%no9N+;0Ie#<+V-frU5tTmW$I3q9Pu<gXXds-B(4OZjkZ-ivEfG6ZJiuu=C#(zoUY8c2l(xTzJY<#QP&W3rlBiD3EyX!=HV9nN+B>KMDt(9I-5qO-O#*!)6g-KsIKG>hvb9RFd!0^Vkk^Y6&Z{qd&^zTeMmV1GLfEgk;j^NW~RazKb&)9OrtE*EcQhGKI0!KNFJ*aLtsFHD4siBk;btW>#y13>sGRhW}<1>yp<9>$yJ(?FV~(I#0h<-NUYh$Hhp1;JKhUc)Msa&SqC~O0Wm93YY-<(obc*h=V-#RL1^IT&00h@-*q5z2M}9R&^}hwRwhd{eDJAtZy}#m!R?#%c%H@-7K9V3;c!sIbQiqzmYP%hc%<KvGg8m3xMrJ?_;c{L*{e1E(8l!-lKp#&Ah&kU#t@&dG{>#Jv%K^701rw2lVRcX?fxt_Cw7^jw%$KNCijzmP^+2ox9c35PfyL<Q{qj{12ZGHG9G9RSSACqb^hM--=m?=GD-Km=_;Bw3l^Tjj=gAm9klRXbjxKyW>!c3)=1}UmpQ<@oqg6iFKxlNDRG7kq2C?>X8D!cPKz%F*qqD1g(?K)A%9`@Fk0n@L&4+5!JF=j%;!AQXTMj5#;j`=HHMga~)`=kr98-b-<S2V4iC>{)Y;JT8UpJ_#H-aty$f^H5NUo6}(fvRsR7rAlC{5M1Nl&u~$Lp`^8{1_oO_SlKh*8uPcp$6k`yag)2H>$oM^+$gO9pxmutIHElFHR9TlcHyZ^tY1CG#OP|FtF<*&16>!{$Qi;_|fcua&=x#OE;XnAkr{`)yT%{>;n@TfKu;Mf7cy@x*oB1T<oQbicJ%^*sq$B?@MvTJ}B+yGpl-Apgq1RDa7J$h}89b~wOYUEETFt~SW5mn4l5zL~0lbt9hvULv6`wc`%tvNP)9QrAL5mmWMY<t0+;{wX+kgf@85bcJ0i_@cIK?J=$zPdxAJ)+LrC476uKS!**acB}3O7^S(d*Jm(i4Yy6FHs$!(f+Cgiqsz?KVYsW?!WAzbMPB&i<(PMD552Tw>S^3&18({^F_jE}mK4%@7G6IAn>`!{B+4p#DoBB4u^sq?q_O)AFk2-(4#nhX)vZ!_lpYNWU5^BF`FAYP{|$#Ehku7?~*Hyuf<^*J)C{9-r=nzBDhOZWiW#-}Y`}A&Lyo;3YH6TOW?xbTPFcdOgkOA(=jr&dkLUSC&njEgJdub}(kv6Yio09E=ZwQ_2bxXB;HI1e?gV^Ab1wWnk!T6@RqX;CRiC=Lcrd@2(9Eva1Fw%^m5FrBR*%eYEt}gjYT57~<9=S97>Oyz$4k%@985n*C)FAEv~|`b3-v9=3CZ3|V$$$aISKU(F^7$5Dm4D+$Y~4UBwhMkjHvI8{7k)6=HBY<{s`K>TEmQKlzZ;kHlHkDWH7O^&ybHMgqxx<<+l-QccLwXq^HEZs;JmbfcmzUDR)P|+=ev9Z3(WiMIk<$6BugJF%983}%?wB?Kc@{@2j1;G|!1#X_pL(KPocWqVHYyCL`BJ8Xf&c<z+@0>z3(<ZNPAnA+t>g2|%UNC=!!~fbWZ}gT!gXi2DxJT)0-K{iDIQz&aoE464YP%Fsrgp%OKCcEPqqyaYarJC^`<#=so&qzKIP?LJ*#NvH%|{Z#;^tXYK&Nb{t^{ch4u+-~`PnCdg@+OPDxlafs1h^Z{rAX(k)IiF3o%m+rTj^O>UFzz-6F$RC{7l7Wv7;3ATI!2Gjfw%GCh8w$tWQsym2=UN%;gHH_eCr$S3)4{)g04xs`WczArel10uEa$DwK<eW6|?onHv?`_^1jHNZdnjTRz}L3YoHNx&fej<Vz*s{12!b<_NAKw-_2P)#2A3pEU4rj9c>Tz1J&Q_J(qMU*r$rcI@_mXNIrippwDGhqI8=#FW9QHn%#_;&3FR#ptC(44yxrP562bYwV^lLC!SoN6J3LFytZG&0V?o+hq)=TDA7u~AxvkSFdOYK8cbB?MW$*_1GC4~Vk+3HII`Y>jTEPCSxHT9*`?2rpCQFeP$s0Yo$RtvW7Rq=la!s9DSaF)h4&MXRUY1rB_80rKjwHWll?4fumY$4=6u^RGzqT_z`n>a8sm$YaXJ8a2G;BBFICY@ZH~iAzB0XFFG5k_RD--N3SH(&>F1<8*uS%{YJ$Wg<BNqRZ$saP>Sln_F;R%PtJUn`I=B3O$%G6hrO}<liq)@yPtiC{9;#Oq&pgvLNZTcJsVxo4Ecw9^uN_{SaJ(AH`EL<T5Jt<1*2c;>f^!<EJt+4Q4YH(+(Lx4k=Kn81DBwgyWqbLseOa6;7@PdO|ecuOQ$?76Mu0q*(}jIG*Gx+>RmIA8+~yEWi#=;&Qai*>7s;$8J&WQx>cK+iYN#H?PV@pU>uK1`llEIJpRluv<lmC-mb}C}cLf{%of#djgN=#_gwqB2kp&i!YmS7PPIp1s~#fRT$HOo)iGvt2b;t*hkYD+i>OZH!s$CDIQbQzKtG->-}0mh*MhR-8s=pbC6c-lbHBjOT(0KUJB<@Sf#>XVh8t>Ct+zXxw=fhO%k(9`)J45j*^slqGmIHy)`;3CYM3gZcNdLrfIPu%&@gCNVpyygJLN4e$Tx?WzgIFV;a<DNcB)feoui|)JMCJ7E?kdpUdbwUd74ag^|!xgCKhc$(yjTvAGS&t{he>-I{=3%tR1EhSHF{;0yc>brvgBe`VI0Si7kC%_C^0#RW}T4vVFGGnyJiS5T`yVgP&uL-T?l)6yuUD&m0ForFN_Gt&bkBAw$zmg=ZbtY8UPF(*K_uxCZVF<{?kjJR+=QEQcJTOxZ!e%Fhd%Kr`np(|`L?&4Zwk&b4zpM_`+AeL-=^!#}-3-Twg?}`3tlXKz$1V=@EHzr*vGG%&4Zcp7eCm)KIcLBoHsQ%!f%gwiHgj(_x(++0$GBDCFkQ-_@9JLYRG1Bx?!OB%-7<E%@d*@Qdk|^>)g5UU$HU4~{fA-H0LR*E3O@xXQFYCymO*@S*n9GsQ@Owx4`Jm?y0^>9xc_s3Bp^kK=1AH9BBDbqmd+@WPx<r`8bD2@tBI1fzHoq?o!yCbG<Gi&|7mS9P{Ss3jP4}=eb@dKV9}vq}3ibaRnJH2>Ni#Xf2gip&vxdtxNKaeFhpRVy&h|@<Bky!P1uClTL&cVCK<nPTO<NK`8JdCV-CGP}c#W_w{wq&a#`lm379$gsd<5rZi~jfG{WpV7t&}Zipi}Fkz1D~C=!A>%E#h7}G3~`<`v_k!E!#%JR=)hivWQhIHC=0-d@VTW4zvevPK<^q$|8Z`WKn#A!tPR<xS0B0LA-I$FJDNG@)_g)v|?G{c`B-E_hqjI%%lw63rdIhfy`mPZcX94VE!vL*<BH<1cU3ds?9Krjll8)4<TwL?c6XB+4cFd|C1)#Lhc_xODs$X5a6}P_Xw9J6h9W^umHosbGhMAQ5JlX6Zi!yO$y<n<bG6tHY|-m3o3N0D)#MOl`8B8w{^A6=8QnjR|h9l$`DSFmHC!;==ea=mE=PceQl)k8W|(42||oo|5!8vidKzwiQW55XpDM`?7J`=MC{wF)vmtv##&shz)(I4Rgn#V&-3ugIZ0rvDtG}JYpfN^`0by*=}H*3iTr8g%ClrMQJ<kOYJY8c07!G&fCVLE0(^44&sdvPHoi7`)T?pJqr7l=_usMglR<Ioh!Pxjt63H!)XGT~olJKy*8vrEFKlMvH3;8r5aoA#7g6LMba``0QS~eF?QSX{0gqEZ0D_8u5pYjGz61ebSG4YFcOZ|I9*GO#+E!2Toj(U=L2H#Z20wLm$G-<RynsD=mO}}DpPs$u-c_K~MkD;=1b-Hd6MxE8Wj`MnT&p>$0e3|y8_p=+nL$1G1kB`}v-vX}ZFOHjb>3%7fpg|g)E7eq)FoTuol_Mjk1j#sg;-nx#Mdt~ZajBt)C?WZyn)09HHmHosfb%Nklstir~*t&F*BrY({E$H_tu_q-*}HHf25upTxk|^W7uBfl(@92T>39)4E`LH%G7#?9AF)^v>wWKse8v((L&m!66_wo5GI}d{d~*q4D-raFJ`*?peWNEBK`;iFn=->dP`M<%-^Her~z*T8jRFBNu9AmmH#Vp65ifU`pAf<vQcOxwoj3fL9E{5iL7$vI#>}G2OOg|6MpfF`2#7!H<sZRgPX{8j@*&=IG<m&3Wu*?+`rU%C81oLYaF@-@!af$_e+94;cV>nuq3XOAv2->1W;yGdu7qT<gV#9rOU>8FM?Z63gnQjiBsoWcQEj*T8SnKa`Ch$4&`sDRz`dFWpJJ7i_Nr&6#&CCG!pWVV5u$zBJZ8)H;>dKnfI8@O7^5I9TT?wc3oHxxzDV<TPYRN(5r>P{lp<#D@?a0gtHUa`$YH0d=IRCq1xc6IPdE<xA{7hS{@=@JKPSisOtg%;Z>;*;L|3I7Je?HBCW16OjcA@@{J3l+Q8W;urQ`4(u&Iujh8A{MM36zi%krC5AZ<v0|sW+(Sic&=*Tz(kcW@p(@9Klov2|iZ3?m4%>IrT61tA+f^l7;Uh=38Or==i2xyxi8IQYC`JY2m_;3eo66f9<@{p6qF$@M76xn2W`rJdwt=AXtHG>WWDCBss?h@Tc*7Exro_cldSq0|faXz=JVu}+dRfLmwB;c5R*sn#%BVRjx(L(Z;K4~6Mlu0-pBbIZkvr>ObxSHF`#@Adi7DWtM5fUUOBT0ZD3)U?gYI<>q*oJt%A_~Z)x9YcL&})T$1}f~-qwJ+ahmwjxm{=-MqqfK`n|OCP_RBh}GMMHPygIsP2R!nrA7CP@0W~c~R&k+!WHCAJ)!S0SeKXK31F57v)#eJ}jO&}c&Poh*#Mp|_UCTtKmK1OT5S9oq|JVB8`f~l41Ua64!JxLA7^vqq7&5@O?|SUoNvZyy+coJ0p798KLcM8f8|8<oy5b+b#G>4jF&bn}6v3dS(C;eQ-C~F+#20URNXd4y(ss-0!$MjbvdUN%|2uWt4p%tyun}i036FXzgWthU<wa7m8Nm<Py<k8Rja=nZMFkhv)`hJHA;r1*Ny7qKEynxynCjL;&JhhrmiTy#av3{S+UW6wqC)*6elNyl{n6(v>GryMOzVf_n1!LM6b}~mCL4^Fqd+R&MFJ%K{jO-&70Glv@n0OETj(F$Xc+XJYcfh4@O6%Aw2@8ru7+-UY=<+H65~6M^B%G^`L8roJ?oyKOW$yMHju}9Nu?ZXc^+LQdeF(~TZC@oyek8z(#f$P?LPe3=IGsdsHKh4tJei)7e4s~GAuXa?pj2n&SDbmh!l+XE2qrO;=-k4i3X6kI-}M&+(B0XS&#9S_;*Xt<kh&4`6VZB5^4@xL{t+>rEG3GK!YlOj7iGNL~iblqwCvr#Ix<wxn5SukK_l#W7h_xB1nAVdIP@%!tJG791P!_4n{BN$ZlW-S9?}jaP|(Mp;Zg1x%JNW7T)go3wl+@q5g&^t%Swp?K-8r1g~`?rFMI_Q^b|49^?;qf~CX)rO&xw(e&@G942Cz?CIdrW@Q4BTa1JK4dZ;Z0qw6X3f3bbcO;R01V3`V#Iyw-1>-X5OP8q$Gx?ZZzKpAcIpZTV=g$JDH54qxr*bpSP&yUrnSze7{(ppoOSbc#@f|iV_zMeX)d%&dS6>URjsedw#v{ze40w=@#Vz)4&rDz;!wPN+SEl(1HkL>_$x=J~bvkg7?bOPC=`eDKGZ243-h#+7gHy+^*(-vWX_4XP&ojng(p9V>qm)<+ACuF=q1)jEsfijFP;W_&bexU#L6@?7uK0s2E4LVEf;e3^<2`1cF4fV(u7)Df2^};|-^}H;L5HEZ=Kt=Az37W?HRY}E$%88Rv$QD9@vl4IYQc&T()GUsFDjq}W&AMtuTFNEi^3qL-3KEZPG3>8paAm$YF=trgFLHsMdMZm{J(Ne_KPTP7vYy17+o!{A|yU3iAHQw>J1r!Gy&X$qL*d06Lk>5+DR-^aOV(_A-3;j<yP2%fSk$!Qig*fWR(*jJY`|dn=rm@4}NtB>3&Q;M{I{i;8i~bUL70@$*{rI*|SIL0bYP)-Hxqh)Dgz)oh5Q7&O;*tNzOt52Kj#v3W?4edyE)-Tm5(AWu*aVNPmbyHZzm7XN7ADEQ&6;XCx@QQ>}*JzDpldmz#^fZT38O{J}#gn5~%#oOPNx#yI?W8cPK(aajz<7SF`*lWONF`Pj>eVW69EW#y?|SE0aNE0xc@EOy7y2wJI=i1HsvD<NHX2qm2?p)Bb?b&*=H5J=t<Mru(<997gzjP>XQqRr2P%`b4&6?lV-BE2^F2xp3a$_uD@JUUpXf^e9*2gZv4hB)v{#KF83r*q6~v?4O`utVBUe~&+0_NEwWmMvt((f_o}x&!*W)@dO~7Tyo9`EQ&QTqGJH`ypkg$UhZ|xCENV(yo}xkjkJN62?~G1Hr5YL9{Rj;RRo|ED@mIe;f6M9DivL<X3*#BHV}M@A03f^syR>`}d^cdgCYXdE=^``)`p;!4+*Fb6l6hZVN`{IT)K)&G=P1leB%udhtO05C>9l!#}fmYS^+-;SsM!*E3cW1Pb2L6mW$bqV^yP47xf#t=^i3)u4xkjL`WQ3xc0T2PtM^2{Cj1thZlJ+C+ihx0^85Yj|Bqq=T5_EWzgBLYqoeB+KI=^RGz9!M<>AeSE?|-=*||fHR56jevGqgCKJHJwXpsr58Pk`^T-7nrq7=p?IuWu0@ur(tarqP^pMV>N`;#TLJyul%GIwZP(yD4))&5htzS|)xjfYO8zd;(b(h~mUjY)(D~RQu3|S|Yi7G(T1HS5$Zf&IJs;t>4u0;07c{ot<y4U%+R@PpTjC?M?BYQ}B6fil#f%}<350gslB<XV1@5_rgreUT%i%TW%LeqS1OhYcI=<!h%Ab?Pqsrq6iZMrF#En?4U@9uK7axy2Mdci?U@_PXPr{?=0*h0MJpSzpY{zs(V_9`f4G45ZTY#Hr%%<cB@&53^@`A(V3>R+3JcD*)w2Yn{)dTP7HcJ#aaE@wOdQQa^^jtnf*hie<>D~8fMKEjzii=n>j4q~nzLdtIJF;k(<jw1&fA{_@>{<&Zi!+p>_V2PrzOdG7&pTTK`Qh8L&1-WmZyw8x(oQTdIsg6SJA}wxzQoavwOWNDLUB<KtQ0-ZBaqXVb0a;@2?V#<!tXCm9de_PP#{#X_tc$Ov+h)u1{gCGL|Q}jYx5%Ix|SrgDqSrSE*IV&8yRtP;cs2ER*uV*UIl<B6t#pN|Es{i0<zt_#f+iZBbDGRmEP%q?o+jzHbb`@Czh~pe&~RgeBiy)7%D3c3L2%pgv1CCgxwTV@<$~G^{#O(6LUAN!kW*FM?x3O<@KoisG6>uUO1N!<@-%~f9rbcX}xTLZ`Pb7Cc}ze-62bP{d0F|7#NziUHXCUME(~nP-%kx)YPCfc{9G^wL2eI$tWtpEb|&1&p#X%suC(y(xPSTq~Lw;4c_jM>k4quk0LKpbqpB-ZwzpF8tP5C{oaZEUYI)Mmk9q;vdqn=dna+UXx*p$gP2f_4@{=sCq`eafQ>B+bb^8@0KKysq)i}8jNjVcX7ujccjet9<7shH^yXlwpKh_+As^xK6eP)v6Rt%&!YMv!0p*(*-NxW~zXSIke4peP0cVwSA^uh|ELI4^Ng*>feb7Kaf!zm&Eb3gK`4i^8)z~Fo8#YX}fj|L>yu6Fjl;>@#Fo>{7!>+x%xD;nODCi4An~%V14K;MX-u4Rc5>)+)RF!-<)`b|ZQfi^GyEo*84eALkb%FvLYB-a(foA>Plz2?6p)^kFw7nhy<_xJM;3_>!vF&Vn{piwG9^`YfqS;>Tg~-~N@c62B'

COARSE = (212,212,(110,254,4),(119,269,4),(136,335,4),546,669,-14550)
FINE = (306,306,(159,367,5),(171,389,5),(196,484,5),789,966,-14550)

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
    need(energy['THRESHOLD']==gain['T']==F(-291,800),'threshold agreement')
    rows,initial_energy=moments['compute']()
    pair,cases=gain['adjacent_floor'](0,48)
    need(pair>=PAIR_EARLY and cases==3969,'first repair floor')
    bridge=F(202,10**8)
    mixed=F(526,10**8)
    need(gain['adjacent_floor'](bridge='matching')>=mixed,'mixed repair floor')
    cost=energy['energy'](rows,initial_energy,moments,gain,PAIR_EARLY+PAIR_LATE+bridge+mixed)
    need(cost['energy_upper_1e20']<=CAP_NUMERATOR,'energy ceiling')
    majority=energy['initial'](rows,moments,gain)
    e=semantic(0)
    need(e['checks']==64074 and e['coefficient']==SCORE,'first entropy partition')
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
                              FINE_KNOTS,FINE_FIX_HEAD+energy['FINE_FIX_TAIL'],energy['FINE_REM'],FINE_SNAKE)
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
    gain=open_module('gain.py')
    pair,cases=gain['adjacent_floor'](48,184)
    need(gain['adjacent_floor'](bridge=True)>=F(202,10**8),'deferred C-line runs')
    need(pair>=PAIR_LATE and cases==7844,'deferred repair floor')
    e=semantic(1)
    need(e['checks']==134692 and e['coefficient']==SCORE,'second entropy partition')
    return dict(graph,threshold='-291/800',graph_cap=str(F(CAP_NUMERATOR,10**20)),
                coefficient=SCORE,pair_gain_lower_1e15=int(pair*10**15),semantic=[e])


def run(request):
    if request == {'phase': 'prefix'}:
        return prefix_evidence()
    if request == {'phase': 'terminal'}:
        return terminal_evidence()
    raise ValueError('expected the declared prefix or terminal phase')
if __name__ == '__main__':
    import json
    print(json.dumps(run(json.loads(sys.stdin.read())), sort_keys=True))
