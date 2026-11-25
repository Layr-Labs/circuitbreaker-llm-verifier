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
LQ=14;WSCALE=12288;Q=1<<LQ
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
  ell=F(2*28*28*41,57**3)
  delta=F(13,57)
  pinned_threshold=F(T,40000)
  v0,v1=2045*57,1148*56
  straight_den=41*57*d*cd
  extra_den=57*d*cd*v1
  y_lo,y_hi=Y[0],Y[-1]
  y_lookup=tuple(bisect_left(Y,j) for j in range(y_lo,y_hi+1))
  def ix(ip,jval,width,sentinel):
   if ip==width or jval>y_hi:return sentinel
   jp=y_lookup[jval-y_lo] if jval>=y_lo else 0
   return jp*width+ip
  ms,rs=len(M),len(A)
  em,er=self.msize,self.asize
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
     od=[ix(ip,j+shift,ms,em) for ip,shift in locations]
     plan.append([itemgetter(*od[r::4]) for r in range(4)])
    else:
     straight_y=41*41*d*cd*j
     aux_y=28*d*cd*v1*j
     od=[ix(ip,ceil(shift+straight_y,straight_den),ms,em) for ip,shift in ray_locations]
     pd=[ix(ip,ceil(shift+aux_y,extra_den),rs,er) for ip,shift in aux_locations]
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
      end=F(vbin-1,cd)/ell; low_den=41*d*shift.denominator; high_den=41*d*end.denominator
      low_shift=41*d*shift.numerator; high_shift=41*d*end.numerator
      for at,parent in enumerate(nodes):
       sb=substeps*parent+1-vbin
       if not -N<=sb<=N:continue
       lower=57*parent*shift.denominator-low_shift; upper=57*parent*end.denominator-high_shift
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
  limb_shift=2*12+40+LQ
  need(WSCALE==3<<12,'split denominator scale')
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
    both.extend([Q-(((total-y*r)>>limb_shift)//9) for y,r,total in zip(dp,twonorm,total_weighted)])
   ma=coefficients(ac,parts);po=coefficients(bc,parts)
   factors=twonorm if cross is None else norm
   outm.extend([((x*r)>>limb_shift)//9 for x,r in zip(ma,factors)])
   outp.extend([Q-(((total-y*r)>>limb_shift)//9) for y,r,total in zip(po,factors,total_weighted)])
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
  limb_shift=2*12+40+LQ
  need(WSCALE==3<<12,'split denominator scale')
  correction_unit = 9 << (limb_shift-26-LQ)
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
    outm.extend([((x*r)>>limb_shift)//9 for x,r in zip(ma,norm)])
    outp.extend([Q-(((total-y*r)>>limb_shift)//9) for y,r,total in zip(po,norm,total_weighted)])
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
    outm.extend([(x>>limb_shift)//9 for x in lower])
    outp.extend([Q-((x>>limb_shift)//9) for x in nonpos])
  if raw:return outm+outp
  need(all(0<=x<=Q-y and y<=Q for x,y in zip(outm,outp)),'main tail simplex')
  return outm,outp

 def evolve(self,k=25):
  m,p=boundary(self,Gaussian(1,26),Q)
  return self.advance(m,p,k)
 def advance(self,m,p,layers):
  n,u=self.apply_aux(m,p)
  m=m+[0];p=p+[Q];n.append(0);u.append(Q)
  msize=len(self.main);asize=len(self.aux)
  ind=list(range(len(self.y_axis)));ready=[False]*len(ind)
  for step in range(layers):
   for i in ind:
    if self.reuse(i,m,step):
     continue
    a=i*asize;b=a+asize
    n[a:b],u[a:b]=self.apply_aux(m,p,only=i)
    a=i*msize;b=a+msize
    mm,pp=self.apply_main(m,p,n,u,only=i)
    if not ready[i]:
     nm=[x if x>y else y for x,y in zip(mm,m[a:b])]
     np=[x if x<y else y for x,y in zip(pp,p[a:b])]
     if self.monotone_main[i] and nm==mm and np==pp:ready[i]=True
     mm,pp=nm,np
    m[a:b]=mm;p[a:b]=pp
  return m[:-1],p[:-1]

def boundary(conv, gaussian, Q):
    H=16
    d=conv.d
    end=6*H
    size=2*end
    T=conv.T
    S=gaussian.p
    cdf=[gaussian(k,H) for k in range(-end,end+1)]
    mass=[max(0,cdf[k+1][0]-cdf[k][1]) for k in range(size)]
    den=13*d*40000*57*H
    shift=-57*T*d*57*H
    slope=1568*d*40000
    start=min(conv.main)+min(conv.y_axis)
    stop=max(conv.main)+max(conv.y_axis)
    tables=[]
    for s in range(start, stop+1):
        offset=41*s*40000*57*H + shift
        tails=[S-gaussian(offset+slope*k,den,1) for k in range(-end+1,end+1)]
        assert all(0<=y<=S for y in tails)
        suffix=[0]*(size+1)
        for k in range(size-1,-1,-1):
            suffix[k]=suffix[k+1]+mass[k]*tails[k]
        assert suffix[0]<=S*S
        tables.append((tails,suffix))
    lo_den=13*40000*41*d
    @lru_cache(None)
    def entry(key):
        lo_num=57*(40000*key-T*41*d)
        k=lo_num*H//lo_den+end
        partial=max(0,cdf[k+1][0]-gaussian(lo_num,lo_den,1)) if 0<=k<size else 0
        return k,partial
    pos=[]
    for j in conv.y_axis:
        for i in conv.main:
            k,partial=entry(41*j+13*i)
            if k>=size:
                lower=0
            else:
                tails,suffix=tables[i+j-start]
                lower=suffix[0] if k<0 else suffix[k+1]+partial*tails[k]
            pos.append(Q-lower*Q//(S*S))
    return [0]*conv.msize,pos


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
 den=114*grid
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

 d_factor=28*s.d
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
  p=97*r
  q=mu*full_pay
  if r>s.T:
   for low,high,w in intervals:
    A,B,C0=line(r+high,p+41*high)
    D,E,F0=line(r-low,p-41*low)
    pay=2*(C0*E+B*F0-B*E)+A*(Q-F0)+(Q-C0)*D
    q+=w*pay
  else:
   for low,high,w in intervals:
    A,_,_=line(r+high,p+41*high)
    D,_,_=line(r-low,p-41*low)
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
    PART_ZERO = {(0,0),(0,2),(1,0),(1,2),(2,0),(2,3),(2,4)}
    split_pairs = {(d,s):{} for d in (1,2,4) for s in (0,1)}


    @lru_cache(maxsize=None)
    def binary_units(num):
        a,b,c,d = phi_units(num)
        p,q,r,s = phi_units(X_DEN-num)
        return a+p,b+q,c-s,d-r


    @lru_cache(maxsize=None)
    def entropy_support(masses, paired):
        units = binary_units if paired else phi_units
        lo,hi=min(masses),max(masses)
        if lo==hi:
            lower,upper,*_=units(2*lo)
            size=len(masses)
            return (lower,)*size,(upper,)*size
        exact_values=tuple(units(2*x)[0] for x in masses)
        origin=max(1,min(X_DEN-1,(lo+hi+2*sum(masses)//len(masses))//2))
        _,upper,neg,pos=units(origin)
        tangent=[]
        for x in masses:
            dx=2*x-origin
            tangent.append(upper+ceil(dx*(pos if dx>=0 else neg),1<<51))
        return exact_values,tuple(tangent)


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
                values,tangent = entropy_support(masses,paired)
                if coefficient > 0:
                    bounds = [y+coefficient*x for x,y in zip(values,bounds)]
                else:
                    bounds = [y+coefficient*x for x,y in zip(tangent,bounds)]
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
            if ((b,o) in PART_ZERO) != (part == 0):
                continue
            for branch in range(2 if j==0 else 1):
                total += case(eta,scores,b,o,j,t,old_mask,branch)
        return total

    d = 100000000
    eta = [F(x,d) for x in (0, 6000, 6000, 52792516)]
    scores = [F(x,d) for x in (6427690, 0, 3616495, 0, 5895441, -4469698, -9346728, 27190256, 34592891, 35844536, -13423863, 1048845, 1727650, 1364857, 9489348, -50623992, -4675462, -43951127, -7299297, 28777192, -2351967, -17985221, -9098469, 100)]
    checks = prove_combined(eta,scores)
    saving = sum(scores[i] for i in (0,3,4,23))
    cap = F(12282281766462836057,10**20)
    coefficient = 1+1/(2*cap)+saving/3
    return {'checks': checks, 'local_reserve': str(scores[-1]), 'saving': str(saving),
            'semantic_increment': str(saving/3), 'coefficient': str(coefficient)}

FIRST_FIX_HEAD='{Wp48S^xk9=GL@E0stWa5C8xGsBn&k;Zi?A_X7X`Zj}7sA9{X(Lcs>UDK{;|^-AuM`4q!#E-ekPIWm{;Da@wB-O1Qsi_Bu}%}7Bawe3;HfU2B=DDI<6=qn4NS&$?FK>>|7@i<4Q#EH8KjAOQ;am@!(NnXQ#TZ`MTI%#b>H_PsnQi8)LI^BIxXWsfSGv9J8_ypgi@c!@-@>oSe`z>D_$scD8hzvYhTg~t83zI`1XH&!M<6QbT1hhiGppv3<c{R@kE_rbc_U}TM9Q5Z1?_mW%_MuNsU}vvj${mm?zxE~UkxrIJ`-O$~=WX>mQr(n50}L5+2@;H_uCKB}jPeIifDl(Z3tsI66@mJPM$cdC%h)Q3vZt)+B%S@#YO$TKV<GUl3K`!7qEo(+G)Cn*HdVJ#ZlTqkfK)JI-q>8b!|_F4i<@rqEAf*YwKPk90^7arbARA02mGBu1rz-COd9)3r10VLrGXl|PZ#Rl^VaBU=pkg47-maotV5zmkmH(q^Ns6Ssb5@>^jyzZ&f=yp;4h|+D9XR!CvVBqM<~0~_h~I?=$wUu-MF`c>rFu>)E+DcrBA=>5HJNuLAy!@Mn7h`^I-1zOi|M|3rB0LGmU`I5D}(>#j@PQ_U@7Z594!QazA63dL(yO6_}N^R8>@?X4{@6Q6;qJo<<K~aA|O!#x<<TFAZjt>+)BUR|>S`*JM8d-Vw7E|L4f1lT<yKwE*y%zJGmANY5pv!=iOGfBR4(#{&}6yB}6-%LbPa9H;P7+X1<<$0D&K5_rKJ>s8j)oHM!jP7INQ2oY6BTgvUKO}zT<zMhl+W=<XdlKg{YA)0CaxYz?-pyxkvjGB81IstY{8=!K>D>+LFE3$Eyj8jQMbIx1>YQf$ePyT?J_dYFCw`ov=A0T8rEN{vL<Dhelu}t2bJkq`f#|P40;>tZaO$3f^WUn~P0exlZwdCqip1ydN+HYKJqF7wZd+A*L4z6v?Bxkv)VkrY>4>e43^I72m<%;Z~I2v<`Wi!DE_H}P}a3?tE#~e^{kq)yU`aqEO9jynDYF0*Np^{1!j^FwKGL3iOzK<=%h1zXC*5%Ud<vDA}?@Y;BV5e3jfjj@}S3X``2FcDCpYq?~ky%Gis{8ZFZx1;{EeZ+f4midwr#IRN6RJriek+pDEt2%6i{UUt83WKw*^aWs=>4a)pOX@Z?W#3F_!G9J4k#nW-Fy!Zlgb=1>_nIxB57@1iV}O`3H;bP=Tbf7gdvE*(Vs~%RMIrA@<(2P=1=)wTG@K?1CiD@M%-h?zwo5CIX%nN@}Og09NBS}jq5axDi%Vk#wdc-+f#8@h38}`_5?wI@b3Q&ZUGw*QNI2M<h*Z`Ru{pVh@%-)YX8&B8I8Lx<jf=zQwCFbhW`r+wA8(d3_CbdR>R^Kvrz$^%}>QQwS1wS0fxMKqb8%T;y5khgv-9ydW-wO&W$vw)>^NXIn$Zm@d3MXXj_rLFyMc)_oX=0Qbs>P7h@c+2;rX~qY;+P3m4hfC)MTUr|oGvQa}YsnZww;)aNfiV%W#1Njw7Hn%?jwMSCxGy`0YR(?|9IPf`9Th3`${qX+Koy}-G<oi^!0X_hgi=5aZ{ShyU>=+l7#u7jRBN54^MMe&>G?6kR6o~H>%|EdEIu%<FDO5nw)s`9*ZaL<oXV-MJT#!}F(f+wvFrm$s;_uvZTTc**CmnDEejuzignSv7d9f{Jvds(VDI>%0Jp_ReFa`Uh6vx4@DK>4_g<V8?f?4qL<`+Ph&_Dka13x$<a)53KeS8XWr{qUDMi^4+tWm_nTz$hoUg+t28y31(x(aNos>&0U(G`(>sVUZie0D#h3hlcrkw`JSvi~oGl=Qi%I`TqydQRASTD1lMwx^7UR$Ke3_P0hkn1WvBWY}=YD9lraGb8c+R^AVP!Xnct(aN&RtQR)=AMPh-GMhcA6EGP#|xqz4~bO$9I;88a$;|bhgVm49vf*01C#Z|A2p#gl#dZ63V+lhtLwDE<s+PyhH_^U_biD*IW*Q?kr0XgASD4*0}pL#N*mD32)UN0kCl}6kfI7l+aBK9e({X{#ne%&&EgntR)#mrIh<cQEk#<=th=VF%T%0@TW|J1qD@dzx48vaZk9G-WY(@Ve??$mJ%Q*d`1!Z<vrz@L62XU99FXkUl~205OR2tSHoA>%`qD|d3?Ei@+G34lry;6i7hv&ds%#|&B3WT-1}A+hW0^g)epM6|0?7`{}Ik;>DBg9R&&%tu5c7(Lm;T3V$PivP<eOsGT}!q7p#duHGwkbQxmi(9j<Eo+LtJ^`*28iNAybS!+Kcc(od=e>~W7*`0|t%51!TSZ50jja$rVj!<FP*?uvlW-1Ih!qgFzCZSiJR})fcC|U7uBXAMBQd!jSPTo_!&7s?Pt{mY?t^c%+OO0&eSo*o^KHCOW%#kQA>6(S*|n^vQQ%N8fA)pUHBxT{bBc_1Z;T}K%^x*@MF7ic3oF%1&6*ZC2syup)cn5#Yy9g_Lqm*bP53|G78GEv$GxM&XA5&ICnWgfXh&XD*qp3-U!KREL9_ME!K#8=40A_ow7O?`suNEu!lq4x;1iOJxS3O~D_fB|j0+vzAu`U*v7O+-VV~)$=BjVFzi2=EFtqQT?%nwZKOqY$*!Ke^iK(RjSmqjPB6+1v{(V6Q(**1<{)IO%t(=PW41h+5z?G^tMM}2j_MzY^hDm9C4V}dCjfZ$0v^JvIaAkkNUJ0M8QaBDT!Wjz~^rBjJV>l))zFp&%FTe+JXul7ig{LNiPF;x4QT24+><jiq8wvZyd07MTIo-pZu+xghcxS?`ZiRJE#&3~kwa}krseIL`jqE=Az;)r#GLnJIkP@Kf+>l$$hmVY~S}coVY*KBQdAKR*c@Va*I&m_vwgE`KJ8bXrur^=i6ABVJ?GTl;uO#sTuq=!Y$0{v|w|Xc^Z#HVYF~R$d<0NkoDGgU9L}CdgzW~ayc?u)wiSzG;x?_-Z$|3_Z<@#!?i<UdBj}ng;PK)yOeFbvREh`!Hf>dniHpYa{_4Gxs(8H~$WlT3T3dHA95hl*+<LPylvu*xK^#}}`?498HSIlV93#sBm0Ec?bv~WkASqx)SwqOJr->H`b^PjKAj|EABkB(&JE#|FbGPAbA%q<(jaD$EjPpUQv1Aq$pKkJ^=yGkJBHiYM#yU(e`$p*fVlRT`yI&24n#j$FC8QY)Ey`O^yaP1i;p68;vm!MFMu?H|~i9#A-#ZuIztaXQE?PJ<To*^JO$aEGU;o}>ai?G=n;?_3>-CPXnXZf?r*>C8tKYg@QT)W{87*En^!pV^n5-_odv4ecR0~DrMeugz51~k)S^!8SPV;_=6XuGTo9Wn_z;ItF&gd^{<Y12KOS2Lo!sm<h&gQG=o3bc2%5I@Qeqkh*KAVoH5=pX8=O8UmrIu^;VRJPsyy3A*S=PA~FWQ$<It|UZBrF})sJ#cI-CFBMKl@Uz+Y9HbYI6pQ%qr3PnY2vy!{OIIFrZA>p5uKs(gwdyl<od+g%iS9&zvD`jG*q-Zu#m8Az)|0x#<C9P(mz6?VfKX3DJW;3h8%IxtC_>%5@1x}JCJ2+G;RdqC_&^b;SNr-pbpE&9y0Y{ZJhX@zw?XJ(b*jL1XEqHXYZq|V$HxH0jnAl2g?QR-4C*7MPRW%)!PUaJLlRxzDlFh+*Y2W-m76(S?9gWu8n0<p`BdSXI+cL2~V|<^^;bx<tfG0VrVG!EyLI<aQWT!O&#Xd(~}&NT!~1veJE+go`li2b2>26u_1c8XwYJ<NQoX6XmBX`>@wvdYW4V?4d%c5v)%`QY%d~O&{Cz8AShe$k`SciZw{uKtxv9DJ@$d~p}_kf<<$HP7@+x>{&fLE#4QsXdBBghMb+8T!%=h)J~smie@+G74o5S`f4BAnCu7O8vaxTqtIwjTJBwZUX8}(1Rz5OQ1%!KVWZ|Kdkhl=?0`p%f6yHAK@A*{k)vNRN{74+fNg3g^?EMIuF3%A{``S8qhLr{7mo(m8D<Av(ke@}6kNL5s0gkWRKzUKw8!7+TgDE2_axgmvK8fd4unGiU9DnLh$D7<`oFSQN4-MWx(ZsN3e;Sz=oP`b>GdflN)Bz0GlyW4Bgz09%8}>DronAf+_2j!vcM~cLG?4b?x~D}+R(!~s8XX%9+r){T*52oRNtJ72z%@V!8z4n6r?C35bCm91d+c3p=<ZCAGcpRewz(ZT#ymLr?<@*UaxZgQNJKu<aqilK76*@K7aov8VmKo#1l!ot-p=M-11SJPQiv3f(NJ6i7R{0xZfHWX=xw=l;{W!Gn_9N$Kal{+K{9Z#hE1zHIN4ve=hmj~T0RsC@h;SExnHMsaIps!Mjluw;cg1^Kbfk~runQ$os0HVzAz<_l57*7$xQ0~Ai-E6IiuH4iU+ofV@#(=CcRRDpeo67w3!D-RAkZKN2!5nt!TVh8jDcub`glkDLZ7%39UWch@h*v#)eS6B86;eyCrk?`*)7WJ}g4!()X+ObONtJS;McHL=jMB&gPUcZ*Q)nz{bTd^8l$-xWXVU7zfRyoI@HiE%s`x%jpY0P6F27)9AgPZ6He~+yX_T5d$9Rygxpxj*D+P?$2>Azd8goT6$cdxk{$!b8Nw{jFVb3@rh+Jdm>8qVgG%r>1HYHG}DqEE5ar53pH!!EA_3PsTUsI_xtX^)i*3)91v1GU^nxPOWbXf1UN?2lW;KYdV1|Iv^doby{D^*F?S@z`v2Hib-RD#elk;yg8+62AEK!H&F|4K_r?MtZrNw_L(>u}P1hWiN!5EG?nvDAC*5;a-ODCEk*Gw%19hrCuQ9I+rQOkEsO2kne{#QBI#42eUJ;;(TfpuG3HG*bC@AM|U9s;Dp1mdZJTXgPY}V=qCRi|i|N4ZyNu*!vb?O{{Plo5?ejwO~(Ux;24K^SaB))%Orm#;Fub-E)a849`agmPdk9!B25<Z0L5x$a#H>efyy?^Ude8K%)mYBDlb2{mZ6T~tn-2an>U{?7mi!o@As5*M8_Z&tAVksudR*k0BppNFYb<71<(-N%m6;x`b<5svIzj6ABe|gCDD%5G=5)%iy-~|0JoC4;g&xg#2cn*N~mIYq97Q5E^&KEW6lON4n0#l0am{IW{7W-8)4!wflh<2$_)T{oKtJkJBQLkf~`94jh4&Xl|v#6_g^*#xj8IJyRfVNI8Kuwjr+Yy~>C}<sb7*+!$IZbk#Ev0HYxRn=WRv^?z87v;#yG~<v6gIrQuo4bTfX0iDXR%YqE|tMjlcKx#u*<ga;XaU162ljf0r0)LYT>XKkgOtpYLkud7SS*ey6!LQDpQ~DFJckD_GjV#ADzts82_r>)zPg4Lu;U(-&MWr9*fIWEw=Dgz2WGB4-4s-KCO`h7o2$R#gXVxuH%4~2mi?#(dv+3HW`mgg_oW2BL+D+3?iq-g7uple=};uxlAm_c-=N5dBzc)ee8X~nKqw#4I8tSn#1_dJY!kgK5f2QcgeO^^g@OE4>S_)nJ|-+{rq8bdXn51t?=R|etEd_|ELRZzh62)*Nx~1LKvN^ZGzh_Yc!^=QACZ6jZX3UayL3T6O_7_3_`Z?V|6!lWG4&8q=hZsTd{VWB<T|bH6pKdy1ixiG<&LwX6F5Ww}@rg5nB*CaZynMbjA5_LH%mT%Sy?wR7Q~f+?X#rA$LERdf08ZFd*Dr+6Z!$ueB}<TE;Vc8uh9KeJn;g#czE89Ed?V=DX{?gtmH3r^`8fDzaV2r!HvY2ykFtV-c-%;r|BGvk#tJ4JXYf-Ta_?3UVroe!?d&_b2x$n{K{oUoLsx4@1B&TRrLCd;!*o3gQ8#l{hMHD+gru@Y5jdhs&3mBk-FX06Ai3s(mwCv2Qn2VequdknPz?*K7bfAuwqmqX~PE`KVbpU#bJ*8V1&n!UyX|VZbcb3uyXq#5Ul$>vnRSeP|FaO$qHImH{@Q!Tdj%FDhpAObG=}Q9=fC%_0?|HP1n)Gq2<I^}ai<jN5x<<w(C1<<gLQ#j9n|^<@STV;l&jgVD95Z<({MIP<r&6Bxc@znCPLUl1AD_6+R}YEc+vvcw5lM|!l5f?(9X-}B7)l0El<c4ws5k(;6`Gl7L69$>Cj10-I94`N47NGGo2=_$_uud@>CV3f*=ZPkB~t)Ev>P!(_;$t)Yz_*<1qr@dfy`ZdgYc8C+#lhL4PTme{uLYM<*GUY1fO*Y6mT(XJx|6UCO3OWe_V2m<$8;AHzE)`Y;?;flqO}a2lzl$!=YnH248#^&4Kn5C8J!Cs{LbA(}WoV_-WFL5T`T;KN$@!r0+(Z>LO&(tY#Tu<^4x14jH8Zs``}19MA^cG3{rgZx{ZXv)%zl4;FnN0Y((C^J!tvjRpxRnLK($DN?x2W6nf0YG&1)Oa_J!sAnD8|=5qhLn@Ghxd0vuQt&T4P(*W|eAt!ZH<@9)6@Rd@y*Y(w>0e~qF7grv^(2O&VN$>Ql5l=a9_6AYxEMtqf#w^Tqq`o-U{-Jd)rP5VX1Gx?RTi%uGkLl_9eu+YcO&Dr65Lxu6>Q!&Z%$iZi`*FIT{7y<5}_rwDxya&J8uytHghcMD#VCnW<`BWjlWt)RE@cnR@3Zco?4e75S79X7!NBNbt+RwfV%OZXA;RYDh<RxrSBH9sSw7HofFW`0iGuG(6Bez!rr%fw<F$H8Q(^Z_RqOOZ%+hqanlruZRywz2Dhah?<RU&Ye%}7&+RnBvCK6s<ziA!08i#t~eu1Y}EAl6~!!2vuxX?si|-}m*5Sx<?U;W%vu7<%h-?|c^&E-<&(ZgX&)TcMPF=#Ediy_cUo4xSC|*JZQ)q@p8+t>_Wn%;Sui#J+?9y*Ue}ILJ5{eb}_aqjTb$K-ns;F@C63_@`m;T2j<^hKTQ-H@gQ%HZDX}u8LfFFx?Yp&{jQ0AbtM6Ft0+aF4?UgaiJ(b`(ff;OM);oN{HK86wJ5Y>23XV5)G4a;N--VYG^&b4dkU(YE2~euQYI?jzTOnDN2f^!OMTvXHOfxlXKkF6D~^y-k4t`;r6goz7eB|szS?Wa^QYK+{tWoW;usQ9eN`pL6Hz$$$w5y{}<!Mbxmasn(+H=hZR}FU)3-TYjel>Oh#X7d9alb&A^jKB;|E2+4c1TxI!Qo6E>gAU?@Z9xalG$Upba8>_PIu)UZ-)R%P*0K;;U@lCE5JX<LrDwk@PB+`ZAoWf@Ljlp}^*H-<bdzkpwOmL;B9-HPHccfR3Ag{wNW7XvDQ8AUtUxeM?p|Eo08Kzwu?n=?m+S9HIanaIhl$nWw#r~~tJ?8k#sVfqZ0lbQGBDe0$<Wgdkg3<}{Eu6P%(Bm=tjkjAa36+X4%V;X&zxdP)oF=cWBcxzHWagNRh4E$?Yh#gnv&&75hCHi_h!MHId(wq<hiL626+T0mG$IwZYmijQIbr37{95w}_^#^qri&HXlc0m}#@%(TuA0fAU3-ok|E8>BRzEeQq8|8pP#Af!&PL(?aWmoSgZ_gTKGF+baGLMXYNseEDQC)eQ_AIOi)b#g}6|Dp#(UwO`G%<VPKNh=5)1P#7@LWOvn3KFZ`>eVqV7Q3LYJxBtAUYSat9jW(|M^Q+T0&3Rs!GMaj?-gPPw_3XD&7*a8PE|nBnk8Q4}Ir~nC9>WZSFH|0YCF0zf<w$13CE#AyH9O8-n()biBO!h&RPA<3M!OlR&(jrd=>)x#JKh-Sp;C0f{+m>`E_xPO1U2;+F{-aUw`USlZv^)TZP#PrLJrf`3moTd-gj`7)aF9QFdR2cXo<80t1jbW?h(K`%-YL!jkB$=X!B`*e!|CmkV6a=2RCV+<s-U@!{F?aoZU{+)~s6UxJ7KS<&1o`oQMJhYPIqdwiIJaB);V50TJ`SsAAa)M|!&Dh4Ay9OU<8PB(bfBQYf-`4Mj^@#zfMBUCGSX*}c9bYDF=P2|=Sg^5Mj@JY^FCZ#LJHqm<`Zwo@qa6N|>Z7%;Sl1LTt*Qs<Tl>If@onO=4OKekdB533l<t({D61@&!q#ZI%YF2ORZPS#e4Jn)xCpv4;TRr|nH2SodEx5qz&VaZP+^@NL`pFq@dVp^Bt+PEv(K_S0f<jQjYMD^t=XzLw#%t@E%Z>-_Etl)q;QH+-@~5F-`r*Ggjy}WGAKN&UA`ls-szaaN+lJU+iiv|)<s+G7qj<^GeM61DAho2w{i&BEmSmTKqWao>V{FY>hIp9pK#pc<PW?KpUo-_;S2@wH*&S&z!+_=>Rxxx_8bw~oYm&aa4&-q-HUF{I~EKyUa%#UlLt%6tpMoL2`ye<k@T^!pzjihqP;&l|7GdO!8{$^ew?*mmsZP+sgXI_srOHgpP|44Fe+h*Ko|;6o<=7i5;>OcK&-1EtMns+8<B+;iB~yo`44}5emFW3;M8f8+4HcKDJI$RZkx!MSn(^=F2MqQCZKj^cTLRx(x26U$9Afw)rF$N?q{rqo)oVT@qGGUeGDnhrBw;5$dr!#SLX%X2h-C~TGf||U#k>^{%37!{~I0<w3b1u6dhs7{TGbIsy`d498WkzYc}y_7DBx*PrPiI{9cuFzX^ukG)w5Yc3)A+7ntq{#U2UAwRj~8Sy{DrfLTaL#F35ChVh<OJ*k>^^gWKnLMKX{A|O^zf5@$a#8c4rKfJxRs$NM`jQXd6N{L}^iM()Mp<<7S`0;x{^b~0hhaSD8qGY5rvLbnHo-0qBout}6pjWv%rhPrP;0d#a#|3L01RqO2;(?<&p}y{PX`6+&O{;^;&}QI8%Bv6faLCB;y<juz3`LItEqNy9vk`8&C+nPq^TFBX>%^ky&jaj$QaNeW=AQRF0-jVtm;Er%H_^`|?YST$LbXssw&Q*C$vh()g>N-=-Wu?^0*D<c<sd;9dtXFb$Vo?fcVcazP_!28A?-~S<E~jaeUdxrDh-w^T@S1Drdc-~W=7L%TJQJo!&WtGaS&?~4zXRav9|7A=A|rhfcjgor?uY7ogR#_{Q0G0zF;X>WDkf0RWDZ>p6!e>TXNuWl>*lL-xN8dDFmX-H!tjw7^@e0HTw}wa7bCi2-jZRzjmE7_5pi!_7!-5hg8$rP7OYe24S|&yj*0KuhlU<mT(u0Dbb=5_hc8xCQlpT%VrHG;0cb8U@d@HW{g|)GypF(f;?%Y4o_{6t@G>@jDe{VCzk~^FyUpL2t^$im4+3sQEyXI@%z))q3b-)i3>_Ym~tILyL&U$>GVq43m_Jng8^<`-9uuY+?YI1klaZ2OT__rqz>|@tVyJ6x-m!3Vfq%O(chS_>9D^goTxkw!IF#@nvIpPtIoI;$%P~p+}^e>{%8BN@@=|HR1by(G>~he@PE00w9b>)0qu9nK3Crv1hvHvX4Pw3YtQ|=9I}0TGU|QdP^s1&6R8mVEOXPmPYwF+P&6%{plhi~6uZ-l3xo&s6Fk(!S2nDH)FQv-O*&tVMVCu93{^&{==wNIO|z!h;6FXF%@>I`Kyxr#Q+ITs%dgV05Dc^q`rj4M!0Oj?Gw#K=>e?dS=w~}-RAN4Szy9IWYNS+KLO$_|;AkN+m@t^MXtkZ3_3gY2`H%kgxbE!LA?dry-TJz+&hgj1z>#hF_HAm+@L5N|d`3nlrVMMt!E<F<NWT6VW3?DwI1~+l(LYjv=1gt>qD+vO=LEI>Oqua;Vb4VJ)=sn%#r5Cr%}uWt45JfcR6lBLyBQ5SYe6)@(VuV2VR1IALL*1o7sla021{y4HXA2x;h(q)gQ7dyF$=eRq>=7W$jS|;Zzdne)-{@|Vg}|wL^iFsMPopgi>GR{;x~v>V~^X8e|WiHlgx%v1b(%GsVSm^4H~_Md%~<DZ8PL?fHHdQ08qIn)oL5T(vpQXh<aARf&DH{vY7r0VDBV|c;#~>SqY6=bY^b!eh9WUa=x~|Gp1~k2C63Ba8u$+r?C~2NoEd(CWg~EQF8$FCN~PtpxLQl@6Ia#Ix0zFUuC5y#j0-2jL2<A*+-63IzQRv5E5l`5~3nP;zQJ3KDndHns?5>!}^`4u|!=dx*#s&Z^2wKUiyA*K*52t@)<|8^SmAbziy}%sgkozx&$Sd!*slhb=Ntm`}`X6vWLi3z@7<4Sg+g=GygiIaY};gAv2E}Cixi5JSaAyw?gkQp;ta(2eMRd5}=wuS^efkoZsSykf4CJ0s5yU)y}iPHbn<iIY%s7vroU8U=29~PYR3=d)rUzB*)%5;BSlOWC=(Ipv$i7xju%deO|CXVoGZt4|tD2o=^gW={}_S(ehnNg8HO6+i+0b8)mL8;1*JE5dw<%TA5RtLhf&>@rKSZ%$HJyI#ncrh$Glqa|Z|CPCtWAvKRI~&idoi%l;HTXF;zEM;UKj3QPuv<BXdZNf2WYaAmV7u3COLfQ@z;F^RHmp31yk3J=SE>}PGI9?}K+NUS06+>vF>ox0a&!v5*Hnci$G>6t0i3LofGeOMFngAkdA0{VJiFgN+%C<;bU9VJZ8hrusX#{)U_;w|qmJ~XlRRB3uYr0~v9jgiVJ;S3BE0XDT;2$Ty@|1(kvoBS`Y%7URI1EQUEtx#T6G-LG<{p-9%kiI%X$yWgb@fKVaj38WP6<#;AK#W0C{)E9KqZNO@c#<rPL5^83;r}!a!eS%V5crU!N6RY_phr0YRI)|S!V|Bp&;<9tYGi?E!mijCeE;BSyz5V?M%pI0W1>Sn&pSF&EC`U|j3A>XM)9GWi8T|q1>-MQC(;7Mc#K7@e$;B!oknn2K^@em({a)~jj!Ip;Dao#9`A^IlLU~59Dj=$%)g{%BUOr$^!CZyIh>)Llz_mw4|w#L+Q1;uc&6XA(nzzF<%@L`fmPjsOjTi#{;PP{T~~$)Gf80q29T5i!T0#CLWzW8kTw~p=#5V<^*dwnr|V^&YDhgvV*7RWDk&~?K)V4w;2)<O<N}puWT2DHPF}T+Knl2-R(Zu_%9CmqdyXx@s&-&AD@AKP2iWlmm5>JStv5Dy0&lr=xXFgGr(}g|U9#PhlTg4b{1VDXrmd{?8oYIbxXkJA9QBjLx%aW-fM)HR1l(a`jtne{VnzcbtiqDE@6U6MsDFVWQEk^(RpsZW86vOg41xuNf7z+Jr7uDOD%#F~ktJA3(~Yc-tq>!um@w<X>ze|i{Qv2%K-C0lX_h++Ja4<FKK??>nS(#D+-FVN2GJyPjS&a<mJ)0s-@9Rud!WPKIt9!9*4@_XtcRRxw=v<vD4z=0*SAoUoFYJWGLianZU6(Ikw0Mzsh<O?GO<Br@uChh+e<gOBtfGE#f&#qI56CpBZVnJ+3J;w{srg1mqC~eQgv&swYm)yDPo;bB<t2dNt#<sE<flLd&g)uC&lW;Nc4Z}kq~NDCU(&<I8p(@mJ7uAKoU-L)L#4*mTS-uz1`^_xJGS;p3VTE!ulelLPT%ufc-CO%zMam+ypzwtmbx-76RLc?Tm^rQfI-uK&#*xMt3leFK;Ie(>_TBl3}YJXKymPoNKyExyVAl<7@Dy^Op~Gn(qm7(3)ik+VOmyg|wiB(#RZeTHnM^%XjmN<B>yZi8X518syVWg?kfiQcE9FL=<NYJE*oX5?Hl21H@{hVo4v{e&gq)W>DX3W}3Yx3PQyHQ>Tl<25{*+2VB)4w*;`S3hzqJ+HiGzzI&T}j>)j&sK(r>$kVRy5SQNhBu`@~*p0FI*It#_+Rwi1104v4&a^dsW}M5`H`X~CWxh<LSrVW=GVH-L4hBW5hf7&5z4xV^j7bDEu@LG@-p-&P76uvgyt6$<G8|kt;Cx?(eHN3!<WPnA>v#mTxb*gwO2Aegq%n|w7|eO-!tXX$(;Wc=d+v0W^TIuX0?e|sI1Yl8Bgj|-{UUf8LKd!^eGki|8kIr$ka@ZfpIlAphGHfgD`VkLjf!z@F=*!_0Cr#r?6pnkx$+=tfda$lIb}R(f5+S1qJEtJ88JGY{O1K5AaZAlSD2i7^;gmvobDj)AG4Nx^*nKs<bl_flh^en(}=2azaGw9U$GFBb6B4K*Gr)uXj`IPpL2;cojArr?pNRF5Pg?+il(td5+1-(H5M9VRu-dINd=t3j*-T5YXY;4uw?ESxFJUuL2I<9UKE(A_`1FvM{yV21)cT-%_6jp?aK=Fk9|EqD2?=na$YzhbE*(%@_#!66QPbyN1^rok7L-E<11P@4W*?z{i3UPRJE;yP~fk`ycU#VSbp`fxZYqMRQEP?vpPHjlJ8{5NEL-U-0oxt&78)3?P5DqKd*joz_Wj3Yz15z&$lE@qr{H${p|MCNFx0s7u&nQHa^{7@3KouDp9E&C<A7xjYFwP^Z;K=UKqS|(vR*xhJS8ED=Q3P+)*%-=eYe4$<b7xQ8FmLiMPcXC5xo{<sB$eY##YAi}duez`oCLU#+9IK}Rvb-1U|@(03M|jmbctQ&p7kE%lx8)#87NX~f!WQeeWbxakjYS2}ARYdKtHx@p_9(yJ6ou$lqniKB}2O~6Fl%k#}h3mMRuJ5IDm(d09unj$&bK4&zW=sbXFP9&(Nq6s>n%z}j;v;Kvp-erCujr;N05Ic=zYS?#)B6m{ahAF_QD4P2Bto3#>sRMHS1sRK(#&=t<(VjCNp`0w-YsXGFB|6{2Auk#&Wh#*F5KSxw$k@0m9VNR;M#qd@<pi027A=5~vX98j4E#ox8`w`nz<*!)%V|@!iMYri-^z=Fxi2XTdK-ctM-y?wpY*IBJb)C?h0~y<k8OpPSv`0BbzTg7co8EM+{S8@sWEf~Vp;LHD;@i^401b50%Aw!nyS3HH@xAHJBTw|7VLQuV`%o4ncVBAa^zP3)1Yt_uR_xS{^4E1cEWAJd5+D@y@OV|9U$XAv!yAY&KrsK^Y`fmb;uE)mzemn{{w`2*aAM!3uYd8@-`ul#Y!tKa#+z&i;mL81F+upFwqp%aZS&@4XmLbk=$#N+5mh-)I2?`SQhFqU!9LVeue&qE@?ohK*jdsL_01u*fLyl0BLJ4@tQfCD|)`@`s?tcXl<HlWl3yDg#=Yi)dkwWxwO(q*bmuy=>C$jX*|iBv;pkqh5ydj4n7g+Sz@YbKLb-OyArU->v(`S33ynu_Wa%a=U1pLvWrC|_q3OAAqvsY5XTd>;!AuLujhi0Tr&e=l$*sHp_ES~wZqov)To*KU4M|_PNIqWyJSPJBk$dTE;daF)_WTKcY1TD>y~tIRwr$jvssgb)AR8&kZm}`cUp#2MJx*xCN$SQQC%16nB|vhb#G@vIA^v*kTTo5eI<6=u}`(*3DJ#hjx41nH0e8PvPCVE9pN5-aSx-SE>BFDP)9W?CzKHFwxu>CJ6;?}+C2#>4lS%?Nl>aAi`Z^anivoz#O*erv!n7xqsrHW*x87h8aB;|3b{$ZA#EUIEJCSRRr|!bALBhS4s<32bh?U7XnS5S3LeicMimukn0X9yS@4s&A`HK6M-_0-N$^tNvhS^&8JJ9wvTEnpFj~B(Ihzp_$!mF84pe6`T9+gVD7I`KNK6}wugaUe?#a$^rNMlEDtY=G#PI(-o<%n!NMjc<|2pHm6{5@Xu&5D1Abxpit00KVKtrJ9j*O|_zMy1FB!SmMdZe%gZQ5`JQNrH5T}#{W{PDj^D*P7rTxSLSVO-&PivuRRE=;8Os9XMw2oj;h<`ryFZ0-aS<=~i;ryRyNp2(}1&M3P##AI&*d2Je>&X3o^l^UqQ5hy#_WdBZ!@j4vPZ;(14bAa*BN;96FC_dYv9&^U=q9ycAmvKmTib>YIBp-~409G4fi3rySl0^05Dv@mKxzWlETFS+r;cKC(fs`UBwFEs+5a=ABH_h65;wxNjk1U3@$reWIU`G?9vdDp1uZy?tWlfz3%;tYQo5d;g;M+Y`oBm=6TMI1d@Q;t&CbjYMa+sHpikukP8)B09&r?auG)qk^t6xMR-}4N*fJjM%?vwu*rA_Z0AGiEKsOH%b^xWf~Syw4qq8(@Xt&Xl&g}I$FE#5L$BChN(;iej(Y=Sn{vii2BZQUbqi<`K!<qB?@D#tPBOm*NE&P_NX2LJ@2tOB{f)px61<f-#<GB85t-Shaj__a#H#V($&?Bkzy+TO0}Nw{PwKN{%WEc(ILci*+ln3Zy3gMgc~CCi4)&_Q{c4F!0I-1oYiGgX+|DB5BcZ-%CQYtNeZV(wy7E8XjI3phU)kX-2gh?Z2FZ>;&DF-CDp5&RkgHK<Ahv6A2O-`T|+=HnpD^W!QF{k9)gSA@n{o8+G6rJ6!5-odIZijq~JU(E8jl_$ylZElGJ!4?fBDhQDj4~WXD5cDB=g6kz@RNi{?2YxOlygUAGB1)sYo=@2w1SZb%OI<k8MmPfQ7AK4B7`$wYciAVcvTfQAFi`6v$0%6^d5q4Rr%kL_FKz3beHCtnVjmmZE@qNLJO-xH$;!;@HO&F=6Vs=wYj7V4cw}Wvq;ENThmaIA{*YV0)x|LXyG5J;;{5r)SxnmQlaCY=9Us=P0TRYI+cS1IVPca-@S5hxYgj(F0_687v6oGeX9?G*2KfO62sNAh_-2_p*LXSH9f<}jWe9j_X|6xc$V^|mLIHBEP=qIwj0}wT^^U&}jp&U)VI~~JzC7!he@M6HUP<r^A%Wn!7}o<G*OLae(76pmQn4w}cfQ&;k-{D=+poy8+)xZm)}2*H#9yAKPL5H89Xmrp`BHn#YIv*93*Wr9oaCeP9smyvEZcgAEh|XoyrieNTMm`(0b%;$-O+wRFP-N_xjo14*$SOJu<p8JsPJl76`@>gb$$zmpVKsWSPsZC^JD}<M;In7kB_+E8};K4@<BEjEE&;ECd_;t3QpGO0h3>?DX!fLE#-la6vQJhxSH;7+H(1BgIHF4HkDNVdfyQ!U`~9dX2bbBGsR^QJfynJTl+z;gm*fEwOtcJ(R>h3nqf}!NgA8g4uQs$Y=kW!2<V*m>(#m#*#n%zQC17+P<52-vHh0pm~)9c#TOVbI^A(?vILK(4!lS46P&;lqS7Nzy<oIj=A>lobXHz=GZ|^DJV!{u#&V-=uW<kQ7Np-GbK+YT2WeLxkyy7m?!@ISmi3C<CJLsfs5EQ?YCuZOwXXCSm0`Y!6F@@IUW5+2Px)=fs^o_9;M4%flZA1UC$bD^1=hAo8%Or!k_y%CWCNe^mW!6A!USxfh{L9p-2gP`uWxQ7Lt)+Ya;#0fq_iXQs&!;Z64T8)-nr`^GpM*g;W*yvaDKdsei;`#G_sa-Xe%BQ4#G$_v`R4)`Ka|M_W3iJtb)i*#rD7wKP}MWUZQ3rA2lg<I&ltP{N8k)6|2oU5LQQhWDD#%YRS30#fulAa%C`1jJ-&Q&B?2^+gJ=wwbqh$TC9e+hcvvG{tL1>>}9SqtKUg{aLLOfm~%6xnQkB{QVKK*l&1tf%PK!!a4(`OCxHXl+%se9$3QA>lwVg>3mY{^2~YJh>7h%#x_d~q^Q2-geE_U6^)WEdRSCsw%1tWQws2{dHK^LVOk(jM5*vdJ+R7sUS+U-Jx;EqC=~XNL8!8QEIhPL11^1lzU^SP{hgb;&Va-IPpRZmdT=GwgtlRVd23sft_n2%F(th!z;KlEJ?0*(o8=f=^Q;tg#s*NhrDI#661_BqmE(+T1SDoQhu_`;P69h=equ<r-tUi!gCY|q>O3uZUfemXLOJ>Z##gI#{keMmpg`LTqI9k$eUt<ku+htbrv3IN+WhXm7o0%=w`ucuNic5zerU;G|$^-9^hh3#X=A~l!E9Dr6aT*UY)@inxIc1z;6LG_*qD)%BGjpf4mO~;-wQLh$O0=EFpm*NB*D3qCKk5R*H1JgaLz(bRR&j{AN5Da008;{=Awu9=SqIF!z7?nOlGKn8s+onKFxq@d_nPr1VL#OwjnJEcDiT2nPE64|$a!z@ORX>fj{X<1T-;XF5IGMRZ14`nfMD&XnRiLEaIJZI70yG->qIg@v*KZe$)3V)4V6<IH)f84<}EZ_#lGMKgr>gi+-+xEvtPY-&+EQx7d5{}-XkJIEi9bJdYItV5N!K82t1vYi);sV&t`I|Ga&pZM?h=;)d<J@J&TiMAjGCvE*1XeOYdBw{_5^R{4B-^Z+DDy%B|uSEB(Uydtj`o*F-Cf=B%%L&06OkJspW}hyXil?_?ao9gdVSF$XPkc-nPngyOeeJ5<jE=lI;*5d`Dm=l=Sc=d6W?`gd8K;lujusAeP$a2#~Yrh9!9mL#Y_z4X8A)CcPe#u^ux>jlH=(uFL=ghxVL2+N0}&=(6Ti`bd+F1PNYQa929y(FBubl;x-MUPWmDAUVKzU;<!y5$VmA%gGDy9d7-y$jx-#Bh#4qhW3{Kt}{dI{Rw*Eo-w-l?8fnkunW(b|uFU5YjufI{+4nl6ymj`T4w(JCB3Hn}ZF^g```buswPGG~h+_r~WczW<b)ATr^h!7y;}(6dDkqzs`R2tclBQ?cq@vY`|?f?fODe{m&SgIG$qOAgBIu1ecVYA48;R^4G<~Ymha)=nvoi(H>J1A$7J5w6t@td@)mPS{eT&2gdDP3`(PzHiQs$#`lLg6y!DU6Vun!BW}*vttY8`{f)QiiHbFjz!soKiC2!pgo-rz5V+<4K|aPY*x+kHK<XQ@<mbbZ^te4TeUC-;a1z@(2)1b=apW9ww>>=DdEbqfGz}b1TD<1s2(Sy+#30DSKJ@A<$}fpw3JM0%&v1*7iuxgX5cLpwE(#}HX0a5Cv~d2YbMIS-H0G$nC{#7(jS!qf)dZ2j*H;Nl<1gy8i{fCfLmAtzmCl|1Evm&Wp&1lUrPa1>m#ammn=O=(7LT;_k2*MV4?qKKEHXYD>609=#HU)8^qNaHk<33{f{+J-Vo}@b=(;_4qk=Fm%B_r)PSG3@#cJqeinWh5qhENCm+IP>-E4IrCBFIPD~f*c_UR@QUbMzNLh(Sw!6>kAnyDGL!l8g^<p|;3v~4LHWGQ{Fo))wYs7^F-uF4sE8U|E}E5U^_=y|f{Em}IH!6FsZy;}J0Pb10RNfsQ4W|;so^|vFUHmPd)3zyZGAR)TXa5a2k3mfph!yVXBOD+(t_g0^Xj2V=H>2_(ov%Y2Zte&Pbv2IT?tZQ7#xJdN#$dD#+SNsb-K<wOJ+Q~;M)6jH6V`V`C&L&Q6b4lSEE}JSsFjb9VOaqiXz4Yk$*zN_|xx@G74Y69g-BGYP)CcS<+7(;a9OLjZKV|h)2vB?@q(tov<)LpgM-n&dqvU=(fBf-(?|!(A6<7w1u!-j7<&@1&J*~|%QX}H7o8{-nZw)%7Nyv``d2_O>EZZZLkBG1<Ba^}Id&>=FwD=kSK%xf<!AL5?2fE_0=3d)#C)%%~{HNa`Y|uZRTLNIhk=Fd_pj0DBU^Rh)k(K+Y&HxzU6JlDW9~g+ZbE{fHyV<mvwVH--{+tSNf+6_bogYMzx^rnY6IKqZyFzzOker}+N~XIOaY2x<^%)P<G$bM$B?+%YIqv#|I&?6jEhC$+XPg~;&Bvz9P#Z=PA-wXFD<K-S+xe-L^3Ri=feRrPeO<J^0t&rMC!a0cY#-OSt1{qDk2T<MwxgSQHCEaBttVkZd9@q1l8e>jVr%-vc-JHz3ap^nMnqQr1AZ(3fLzzwA@mS48a!r^7z?Ea^yQx-+p{MtP&!ImNz;a2nWOub9xnJP8&98q;$G|F$pYY9N}~m`smt5yU3U<Om}?$32$O2)D1K8ctG?vM3^r|~F<NTXaKR~9BrY_&1|_ll`(NcyCd38ummtG{#N(3lWzKimT0lqhhOGQpqF-jjntRldo$vH~zSly@$JWd{4UZkK<R<Sbh^3&dmq>SR@cg;hn1R{00OJama&p1Qctd+Wo#DLykzu^tza~@{HWAQ%YSn1sF|x^FELMZ4Esg?sTiD%lVnmQg6dM6JZd>leg;Np087*=vBf!B$9a6c(mWPA7#KOobHP^8wEe4SCR>6+h_#Vtex>MZYP_jQpPVtQef)24Ad4JpiZVz=xN=8F99`g%%)AvtQ3);=}gea%_e5)o4kSt8q9Y1|35y>=%R+P>Nxg;qOM=_2dbFid6M8obnjAWk$I}i;)7Ly2<O>Z(MbYvQ)h3lH{Y)te&ZLqN!_a{^lQO(4U=eCZh`>5{ShUl~e_@`syLEgd&4aZ$<Rh}`8^=2CuEsLKeP~ZdP1zydhZALz<hv+M{BtNt@v3~g|tw!LAzjdsRjG+;r#7()8Hjq_DM=fr6Gorocr-q<(o68R<{}u~58+JkVQ|<yHB9`f2?^r(4sNRjuscbXoZVtO6BhJDk(Sr%YI2dt?0S*rjccXs2k#s^$u!VU=7{))A^}`XK$xXB2T<t&G>pb!`!+NSsk8%)T|E;`vh$D!;MDf|5y!No_PY*K=1i}sX#lUy_)0<mhd{Hpf@HqmXx(^$n5CG~(En|D9)DajodEh>0lBh7<(J&6$k<!L0h2D6%$zYn3Q<;+6@a#N;oP<QtG+~vO`+av3*m0G$2dmIY>{=(S-Q?^pGR?uBPAHa5y!H}@noQXIh(9a+IJ`yVKuFU*mAOrnJ|)W_E2Z|lRmG*RvH7OmUcI<kiAEA){g4HjyA%yBtD<29_}QL-3gzIkZPNW(tW-YV*LpvYbCni25o%0$SZo~y&$>l|tHjj6fc<r5QF+Mlrp=A>@RTy1wGI?_?vdymHFAEp53_q7v|?83N+3(a#>p=W0GGLj3``?tmWjMhDt@ytOf@H;yF?j&@NA7r_}%Z{#IDvFONY6Vl(=oq*X!(<clcRguFhhs_N@65&E&YHt*X{u@~xCPXa`0`Nb;{9Zl<w&Oon>g<|;&9GZVdX`W-#-2zn*OC*dbDvKBxS9_kF*tx{u0T(AR)NAz!=)2AwiLrFw%S2o6`xj9i7*Z5Rm>5$R{c3T+G7kUD$N+fMOhzX^MRYW%`H_>fwB^xLdm+BYQ`7=BkRaZRFp2pg`Yt?RyH^Al5EfNlUm|Y{|gg+{|m^hP0D3xVP1te2JTnmb#AT@$46l-e&U}};U_-7xPXUNmZ86L-42<iA9<7QLB|MaZLcp~6?+&VQv{H1gd+VHwh`o$ETlJG~ZzBB=n?_v={!R;_>83a;+1WhsL&c`jVXJu4^(kpUaN&^r|v_)8e{Y!ZfwoPoeszy;kAB`E{7#czlzYe;N2g=P`aC1#`?`&?CQpV?>^Kv`t)vvn$s?nIn5#nW!qy^(v%}|eHwiV0)y1UbY^-t4%$!z~hX#O=m*iXJgt1*^M)L5rK7)yhxm+MCw;qf-Lhyh+=6-CHuF(Y9pap~#swvsvv6f46!E~rardQBhcI0HOTfUIqfhx8*r04zqs?2dk4n=@*(!XfhEv2Kp1hBosN3(!ehcn^o?gJ{dD?18ILwci!)x^4EJ*HK#OV}7z>?1t9%tk)Orr4aU;nrHu^J;WjH6d4?QPk<>MKnKnfuDMe)29<#<k5<J*{;3f)uY|@@RKk@J?fC&s&*9KXVXZN*wan6l6;P%GdfaACU{EJ2r}$=5tx*r*G~ST8#vVv4w9z3;4*f(?AJgY^IsPOb){a;ag}&gpayRbr*cT@s&)N&q>%vN;eJV}H13QLgJ<I19e%+J_`NAu3wz9^9PU0xky-{vK%En#!ProC{9@(SHv??}i_AcUZ^d$mYLbN<S^#(?rfHJC>aE}k@-%|QIU=N=M2g};^2_DEZe^^pC%Uqq<U&<)<m7R_^llMaJu`-Lf6qs7JE=yFoc9FQlNl<W?NEp-3qH>4=Z#<qhA%#HaOJ$FGM?`uA{EUKMRJ021v5LVzE(S}sT>&rJ@qjheg_4~tV<1c%_Y2_}?0Kgp0Rqk)cdewMF_y}}AznfXWmv=}qoKJ~F*#M;z=X7B109AH7Btx1-P~miae><op38W71o%rnp72UCn|nS9m=X4H*QJ&13NF?Qz!~86V(njzmy%AZvU`c8WT~0%(^_qSykGV`_Q&4xV;Ka3p$BHv>O-(CcbD3qxeFP!-H^=mx%TE%)h-Jvwkmb>&YhyCp`W=Jc}Ex-nAfDSJzodO)~D6G^n5f(!H@CA3y4!=rT+9(X0!EK>^)B}qp3nOXlcc@dI6T3kD22pA>Ce@BhF%G4zW(;4=&z<(vie{X!y~AeiWz1$ansE(!3Gi7#V4IF~`>|{@oIfvt<XI{Avba?=@5RWRR!{MfHkDNix}4h7s?q<F3zs<tsH_OS;R>t*N^6!yU46&%SOFeLD*eLWU~oTM@Bk)IPU?Cs0-i)O&|Wy1w@GmPq}bS6ATD$^Cc~(CMA1w7c$KV9>JA!Gd^k?u}71ZWrc>a<2m|fAaknF678$!HiWZN@P?fe7c|!L$3#%W7G<q@h;dHGlwjaj2dmgOOj!kMUAc)K8y*M2?)0|hd{l_5wKq8Wb<TP>vrN6+pch&uXBrN7-Z_5Y8dVUI@#?W7ZD~=(iao50>!-#s5uu_zdM`%Se>Z7mx}yQG2gAoc-@JN&BLj^!J(ccH^tUzw5KiG*?zjtpRHis|A1L^u4Y$$RKe%8SXxN-c4WgI`#ZPP*?T7p?S5lmw89EO5T1hXZ-;QYkxGVskgrCPxuEba3Y*)+iu?uTDNk>KOu0Zuv&vnlqJ@+s$3`(xG*W}RH`X54'
PREFIX_RAW_SHA256 = '2ec0e1b8a121b9795e812bf9c6185984b8e76d5410577ca529ccce19868cafcf'

PREFIX_DEPTH = 25
DEPTH = 33
PATH_DEN = 9 << (2*12+40+LQ)
CAP_NUMERATOR = 12282281766462836057
SCORE = '18836041925354060398745540167/3684684529938850817100000000'
GRID_STEP = 220
ROOT_UPPER = 8187
PATH_SPLITS = 12
PATH_LAYERS = 3*PATH_DEN*Q*8**PATH_SPLITS

def verify_path(conv, start, finish):
    m,p=start
    n,u=finish
    m=[min(a,b) for a,b in zip(m,n)]
    p=[max(a,b) for a,b in zip(p,u)]
    a=conv.apply_aux(m,p)
    b=conv.apply_aux(m,p,cross=finish)
    c=conv.apply_aux(n,u)
    M=conv.apply_main
    def positive(a,b,c,d,left):
        if min(a,b,c,d)>0:return True
        h=a+3*b+3*c+d
        if min(a,h,d)<=0 or left==0:return False
        return (positive(8*a,4*(a+b),2*(a+2*b+c),h,left-1)
                and positive(h,2*(b+2*c+d),4*(c+d),8*d,left-1))
    starts=(m+[0],p+[Q]); ends=(n+[0],u+[Q])
    aa,bb,cc=((x+[0],y+[Q]) for x,y in (a,b,c))
    width=len(conv.main)
    for i in range(len(conv.y_axis)):
        h=i*width
        k=h+width
        cache=({}, {})
        A,B,C,D,E,J=(M(*f,*g,only=i,raw=True,cache=cache) for f,g in
                     ((starts,aa),(starts,bb),(ends,aa),(ends,bb),(starts,cc),(ends,cc)))
        vl=m[h:k]+[Q-x for x in p[h:k]]
        vr=n[h:k]+[Q-x for x in u[h:k]]
        for x,y,a,b,c,d,e,j in zip(vl,vr,A,B,C,D,E,J):
            need(x<=y,'continuation order')
            if x==y:continue
            r=(3*(a-PATH_DEN*x),2*b+c-PATH_DEN*(2*x+y),
               2*d+e-PATH_DEN*(x+2*y),3*(j-PATH_DEN*y))
            need(positive(*r,PATH_SPLITS),'continuation residual')

def transfer():
    return Transfer(225, 225, (117, 270, 4), (126, 286, 4), (144, 356, 4), 580, 710, -14375)
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
    from hashlib import sha256
    energy=open_module('energy.py')
    first=energy['first_endpoint'](conv,Q,energy['FIRST_KNOTS'],FIRST_FIX_HEAD+energy['FIRST_FIX_TAIL'])
    need(sha256(array('H',first[0]+first[1]).tobytes()).hexdigest()==PREFIX_RAW_SHA256,'packed prefix digest')
    return first

def prefix_evidence():
    conv = transfer()
    verify_path(conv, conv.evolve(PREFIX_DEPTH), decode_prefix(conv))
    return {'prefix_sha256': PREFIX_RAW_SHA256, 'prefix_depth': PREFIX_DEPTH,
            'certified_pairs': conv.msize, 'path_layers': PATH_LAYERS,
            'required_followup': {'phase':'terminal'}}

def terminal_evidence():
    conv = transfer()
    from array import array
    from hashlib import sha256
    mid=conv.advance(*decode_prefix(conv),8)
    energy=open_module('energy.py')
    guess=energy['second_guess'](conv,mid,Q,energy['SECOND_KNOTS'])
    a,b=guess
    second=energy['apply_masked'](a+[Q-x for x in b],energy['SECOND_FIXES'],Q,conv.msize)
    need(sha256(array('H',second[0]+second[1]).tobytes()).hexdigest()=='c42d02ebf8338031e2d46905b2185534a56526915430b75bf7eb578427f41eb8','second digest')
    verify_path(conv,mid,second)
    m,p=conv.advance(*second,DEPTH-PREFIX_DEPTH-8) if DEPTH>PREFIX_DEPTH+8 else second
    del mid,guess,a,b,second
    truth_controls()
    literal_minimum_controls()
    upper = root(conv, m, p, GRID_STEP, GRID_STEP)
    need(upper == ROOT_UPPER and 2 * upper < Q, 'terminal majority')
    depth = 2*PATH_LAYERS + DEPTH*len(conv.y_axis)
    del conv, m, p
    moments, gain, energy = (open_module(name) for name in ('moments.py','gain.py','energy.py'))
    need(energy['THRESHOLD'] == gain['T'] == F(-23,64), 'threshold agreement')
    rows, initial_energy = moments['compute']()
    cost = energy['energy'](rows, initial_energy, moments, gain)
    need(cost['energy_upper_1e20'] <= CAP_NUMERATOR, 'energy ceiling')
    majority = energy['initial'](rows, moments, gain)
    entropy = [semantic(0), semantic(1)]
    need([d['checks'] for d in entropy] == [63846, 76941], 'entropy case counts')
    need(all(d['coefficient']==SCORE for d in entropy), 'joint score')
    return {'threshold': '-23/64', 'graph_cap': str(F(CAP_NUMERATOR,10**20)),
            'coefficient': SCORE, 'directed_initial_majority': majority, 'directed_energy': cost,
            'terminal_upper_num': upper, 'terminal_upper_den': Q,
            'prefix_depth': PREFIX_DEPTH, 'warm_sweeps': 8, 'suffix_sweeps': DEPTH-PREFIX_DEPTH-8,
            'depth': depth, 'prefix_sha256': PREFIX_RAW_SHA256, 'endpoint_sha256': 'c42d02ebf8338031e2d46905b2185534a56526915430b75bf7eb578427f41eb8', 'semantic': entropy,
            'root_truth_controls': 72, 'literal_root_minima': 162, 'literal_NAE_controls': 9}
def run(request):
    if request == {'phase': 'prefix'}:
        return prefix_evidence()
    if request == {'phase': 'terminal'}:
        return terminal_evidence()
    raise ValueError('expected the declared prefix or terminal phase')
if __name__ == '__main__':
    import json
    print(json.dumps(run(json.loads(sys.stdin.read())), sort_keys=True))
