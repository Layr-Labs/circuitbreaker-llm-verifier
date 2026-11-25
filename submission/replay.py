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
LQ=14;WSCALE=12288;Q=1<<LQ
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
   # Main and auxiliary both consume only odd product coefficients.
   pack_even, half_product, half_merge = packed_arithmetic((tc + 1)//2)
   odd_layout=Struct('<'+'Q'*(tc//2))
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
     layout=Struct('<'+'Q'*len(quarter_weights[r]))
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
     od=[ix(ip,j+shift,M) for ip,shift in locations]
     plan.append([itemgetter(*od[r::4]) for r in range(4)])
    else:
     straight_y=41*41*d*cd*j
     aux_y=28*d*cd*v1*j
     od=[ix(ip,ceil(shift+straight_y,straight_den),M) for ip,shift in ray_locations]
     pd=[ix(ip,ceil(shift+aux_y,extra_den),A) for ip,shift in aux_locations]
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
      layout=Struct('<'+'Q'*len(col))
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
     on_e=pack_even(v*Q for v in w[::2]);on_o=pack_odd(v*Q for v in w[1::2])
     off_e=pack_even(anti[::2]);off_o=pack_odd(anti[1::2])
     plan.append((fe,fo,ge,go,on_e,on_o,off_e,off_o,corrections))
   total_weighted=[base*Q*Q*factor for base,factor in zip(B,norm)]
   if isaux:
    return (plan,quarter_pack,quarter_product,quarter_width,quarter_full,norm,total_weighted,quarter_parts,len(nodes),quarter_weights,residue)
   return (plan,pack_even,pack_odd,half_merge,half_width,even_full,odd_full,norm,total_weighted,reduced_parts,len(nodes),f[::2],f[1::2])
  self.A=setup(True,Na,A)
  self.M=setup(False,Nm,M)
 def apply_aux(self,m,p,plus=False,previous=None,round_number=0,only=None,cross=None):
  rows,pack_lanes,prod,width,full,norm,total_weighted,parts,cnt,weights,residue=self.A
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
    both.extend([Q-(((total-2*y*r)>>limb_shift)//9) for y,r,total in zip(dp,norm,total_weighted)])
   ma=coefficients(ac,parts);po=coefficients(bc,parts)
   if cross is None:
    outm.extend([((2*x*r)>>limb_shift)//9 for x,r in zip(ma,norm)])
    outp.extend([Q-(((total-2*y*r)>>limb_shift)//9) for y,r,total in zip(po,norm,total_weighted)])
   else:
    outm.extend([((x*r)>>limb_shift)//9 for x,r in zip(ma,norm)])
    outp.extend([Q-(((total-y*r)>>limb_shift)//9) for y,r,total in zip(po,norm,total_weighted)])
  need(all(0<=x<=Q-y and y<=Q for x,y in zip(outm,outp)),'auxiliary tail simplex')
  return (outm,outp,both) if plus else (outm,outp)

 def reuse(self,idx,previous,round_number):
  if previous is None:
   return False
  y=abs(self.y_axis[idx])
  interval=1 if 100*y<=self.yd*self.core else 2 if 100*y<=self.yd*self.middle else 4
  return (round_number+y)%interval != 0

 def apply_main(self,m,p,n,u,previous=None,round_number=0,only=None,raw=False):
  rows,pack_even,pack_odd,merge,width,full_e,full_o,norm,total_weighted,parts,cnt,fw_e,fw_o=self.M
  outm=[];outp=[]
  fneg,fpos=(m+[0],p+[Q]) if only is None else (m,p)
  gn,gp=(n+[0],u+[Q]) if only is None else (n,u)
  limb_shift=2*12+40+LQ
  need(WSCALE==3<<12,'split denominator scale')
  def odd_product(ae,ao,be,bo):
   return int.from_bytes(merge(ae,bo,ao,be),'little')
  for idx in range(len(rows)) if only is None else (only,):
   fe,fo,ge,go,on_e,on_o,off_e,off_o,corrections=rows[idx]
   startrow=idx*cnt
   if self.reuse(idx,previous,round_number):
    outm.extend(previous[0][startrow:startrow+cnt])
    outp.extend(previous[1][startrow:startrow+cnt])
    continue
   fm_e=pack_even(map(mul,fw_e,fe(fneg)));fm_o=pack_odd(map(mul,fw_o,fo(fneg)))
   fp_e=pack_even(map(mul,fw_e,fe(fpos)));fp_o=pack_odd(map(mul,fw_o,fo(fpos)))
   gm_e=ge(gn);gm_o=go(gn);gp_e=ge(gp);gp_o=go(gp)
   zero=odd_product(full_e-fm_e-fp_e,full_o-fm_o-fp_o,on_e-gm_e-gp_e,on_o-gm_o-gp_o)
   ac=(odd_product(full_e-fp_e,full_o-fp_o,on_e+off_e-gp_e,on_o+off_o-gp_o)-zero).to_bytes(width,'little')
   bc=(odd_product(full_e-fm_e,full_o-fm_o,on_e-gm_e,on_o-gm_o)-zero).to_bytes(width,'little')
   ma=coefficients(ac,parts);po=coefficients(bc,parts)
   lower = [x*r for x,r in zip(ma,norm)]
   nonpos = [total-y*r for y,r,total in zip(po,norm,total_weighted)]
   correction_unit = 9 << (limb_shift-26-LQ)
   for at,fi,gi,mass in corrections:
    fm,fp=fneg[fi],fpos[fi]; gm,gp0=gn[gi],gp[gi]; gzero=Q-gm-gp0
    dx=(Q-fp)*(Q-gm)-fm*gzero
    dy=fp*gzero+(Q-fm)*gp0
    lower[at]+=dx*mass*correction_unit
    nonpos[at]+=dy*mass*correction_unit
   den = 9 << limb_shift
   if raw:
    outm.extend(lower)
    outp.extend(nonpos)
   else:
    outm.extend(x//den for x in lower)
    outp.extend(Q-x//den for x in nonpos)
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
  ind=list(range(len(self.y_axis)))
  for step in range(layers):
   for i in ind:
    if self.reuse(i,m,step):
     continue
    a=i*asize;b=a+asize
    nn,uu=self.apply_aux(m,p,only=i)
    n[a:b]=[x if x>y else y for x,y in zip(nn,n[a:b])]
    u[a:b]=[x if x<y else y for x,y in zip(uu,u[a:b])]
    a=i*msize;b=a+msize
    mm,pp=self.apply_main(m,p,n,u,only=i)
    m[a:b]=[x if x>y else y for x,y in zip(mm,m[a:b])]
    p[a:b]=[x if x<y else y for x,y in zip(pp,p[a:b])]
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
        tails=[S-gaussian(offset+slope*k,den)[1] for k in range(-end+1,end+1)]
        assert all(0<=y<=S for y in tails)
        suffix=[0]*(size+1)
        for k in range(size-1,-1,-1):
            suffix[k]=suffix[k+1]+mass[k]*tails[k]
        assert suffix[0]<=S*S
        tables.append((tails,suffix))
    lo_den=13*40000*41*d
    pos=[]
    for j in conv.y_axis:
        for i in conv.main:
            lo_num=57*(40000*(41*j+13*i)-T*41*d)
            k=lo_num*H//lo_den + end
            if k>=size:
                lower=0
            else:
                tails,suffix=tables[i+j-start]
                if k<0:
                    lower=suffix[0]
                else:
                    partial=max(0,cdf[k+1][0]-gaussian(lo_num,lo_den)[1])
                    lower=suffix[k+1]+partial*tails[k]
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


    S = 1 << 54


    def ceil(x, y):
        return -((-x)//y)


    Q = 1 << 80
    X_DEN = 1 << 47


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
        e = (shift-47)*Q
        log_lo = e + (Q*ln_lo)//ln2_hi
        log_hi = e + ceil(Q*ln_hi,ln2_lo)
        return (-S*num*log_hi//(Q*X_DEN),
                ceil(-S*num*log_lo,Q*X_DEN),
                S*(-log_hi-inv_ln2_hi)//Q,
                ceil(S*(-log_lo-inv_ln2_lo),Q))


    need(part in (0,1), 'entropy phase')
    PART_ZERO = {(0,0),(0,2),(1,0),(1,2),(2,0),(2,3)}
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
            return 1,(lower,)*size,(upper,)*size
        lower=units(2*lo)[0]
        slope=units(2*hi)[0]-lower
        span=hi-lo
        secant=tuple(lower*span+(x-lo)*slope for x in masses)
        origin=lo+hi
        _,upper,neg,pos=units(origin)
        tangent=[]
        for x in masses:
            dx=2*x-origin
            tangent.append(upper+ceil(dx*(pos if dx>=0 else neg),1<<47))
        return span,secant,tuple(tangent)


    def prove_simplex(linear, constant, terms, biases=(), symmetric=False):
        dim = len(linear)
        grouped = defaultdict(F)
        for key,weight in biases:
            grouped[min(key,tuple(4-x for x in key))] += weight
        biases = [(key,weight) for key,weight in grouped.items() if weight]
        common = lcm(constant.denominator, *(x.denominator for x in linear),
                     *(x.denominator for x in terms.values()), *(x.denominator for _,x in biases))
        c = S*int(constant*common)
        line = [int(x*common) for x in linear]
        keys = tuple(dict.fromkeys((*terms, *(k for k,_ in biases))))
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
        biases = [(indices[k],int(v*common)) for k,v in biases]
        n = 1 << 44
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
                span,secant,tangent = entropy_support(masses,paired)
                if coefficient > 0:
                    bounds = [y+coefficient*x//span for x,y in zip(secant,bounds)]
                else:
                    bounds = [y+coefficient*x for x,y in zip(tangent,bounds)]
            for index, coefficient in biases:
                masses = [v[index] for v in vertices]
                lo,hi = min(masses),max(masses)
                whole = 4*n
                if hi<=2*n or (lo<2*n and coefficient<0 and lo+hi<=whole):
                    bounds = [y+(coefficient*x<<8) for x,y in zip(masses,bounds)]
                elif lo>=2*n or coefficient<0:
                    bounds = [y+(coefficient*(whole-x)<<8) for x,y in zip(masses,bounds)]
                else:
                    span = hi-lo
                    slope = whole-hi-lo
                    bounds = [y+((coefficient*(lo*span+(x-lo)*slope))<<8)//span for x,y in zip(masses,bounds)]
            if min(bounds)>0:
                continue
            if depth>=44 or not pairs:
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
        beta,theta,theta_new,lam,kappa,*tail = scores
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
                weight = -(lam+kappa*(1-a))/2
                wu += weight
                wv += weight
                for i in (1,2):
                    for mask in ([3],[i],[k for k in range(3) if i!=k]):
                        add(terms,subset(mask),-kappa*a/2)
            else:
                wu -= lam/b
                if b==2:
                    wv -= lam/2
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
            h_weight = -w-lam
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
            constant += cz-lam-kappa*(1+a/2)
            hb(terms,az,-cz)
        if o:
            biases.append((az,gammas[o-1]))
        if b==2 and o==2:
            biases.append((az,theta_new/2))
        if in_horizon>1:
            biases.append(((4,0),-deep_biases[in_horizon-2]))
        if out_horizon>1:
            biases.append((az,deep_biases[out_horizon-2]))
        coeff = -(theta+theta_new)+(theta_new/2 if o==2 else 0)
        if out_horizon>1:
            coeff += theta_new*(F(1,2)-F(1,2**out_horizon))
        constant += (theta+theta_new)*j+coeff*L0
        for i,weight in enumerate(L):
            line[i] += coeff*F(weight,4)
        return prove_simplex(line,constant,terms,biases=biases, symmetric=j==0 and b!=1)


    def prove_combined(eta,scores):
        need(len(eta)==4 and eta[0]==0 and all(-2<=x<=2 for x in eta), 'signed carrier coefficients')
        need(len(scores)==16 and min(scores[:5])>=0 and scores[-1]>0, 'dual coefficients')
        total = 0
        for b,o,j,t,old_mask in product(range(3),range(5),range(3),range(3),range(2)):
            if not (j==0 or (j==1 and b<2) or (j==2 and b==0)) or (b!=2 and t) or (old_mask and not (b==1 and j==0)):
                continue
            if ((b,o) in PART_ZERO) != (part == 0):
                continue
            for branch in range(2 if j==0 else 1):
                try:
                    total += case(eta,scores,b,o,j,t,old_mask,branch)
                except ValueError as err:
                    raise ValueError((b,o,j,t,old_mask,branch)) from err
        return total

    d = 100000000
    eta = [F(x,d) for x in (0, 6000, -8410908, 14660659)]
    scores = [F(x,d) for x in (6124864, 0, 3409948, 0, 6181981, -8765605, -14541872, -3170451, 33441501, 43676501, 14065680, 971116, 1705745, -2505445, 6288574, 100)]
    checks = prove_combined(eta,scores)
    saving = sum(scores[i] for i in (0,3,4,15))
    cap = F(12293178397485488844,10**20)
    coefficient = 1+1/(2*cap)+saving/3
    return {'checks': checks, 'local_reserve': str(scores[-1]), 'saving': str(saving),
            'semantic_increment': str(saving/3), 'coefficient': str(coefficient)}

PACKED_PREFIX_FIRST = '{Wp48S^xk9=GL@E0stWa5C8xGsBn&k<HiZI#{>XMAi74!aqw~wr`y}cFTPe+?nNT2&x>K!MnVvT5voQYtRNO3Nd)h2znNw!rDV0gPe%Z{Z%a6p6XtnBKwLRK&MtP$cz=aic2nChdO6D@by<TEe34a_;H0O-6jJY1a-CNyyw+CFNL3UoIBJiHy_d~#XO$x?&-mtOwTXy@BZCxEgsrZFmG<JY{t^sSt6U;$C~u-(p#>dG82Kq@+Yip?(Wx#&`cyCrc)&U=wYz1aaMK2jgO!u8LeMlA2#s9KTAY(FvD9so!330cO2T^d=g7e@eYzUfqPMa8t*h}Y+xu<D+PhhuL#rS2^Qn3kaVj_DRPh#{imtP21N>`|rZWM`Efot<I=RPtJ8Za@6a+E0bV|071!Z%2f|fCyWx0yWaHOzyb-dN#3-QB<80UU+t>LbklO<au54%M@NK<dvV#y>u<<Lsxw-5D2Eya8b<$@4Lac@+4MWn}U#Pb;X^&$ZDNJZTWXKCFA5XjCUd8dd13x=hMtk0+9IjZ!<^l?Mir~pPJr^+`-nJ{&QeJr<|b^E0QTM;+dnyB_-w^9yKC|jRlGa{L7@YXarIW<PZN4wejhcQw0tY3Vg#;!o&4uwS%MDwKD*UrE1(a+oaaSB;2@h4wNwgfXm@&7P+2|^8B>k{}t)3A}K(r$e^LA|sJthJTzgNzCQ+`W{Cn+gfs7%c2(d8vaM4#---a<1?RuPfDweyZe$Ex0E$VOnvE(`q|1=GJ^Jxm%Ll{lo5OBUsy|hhP}IgH1$KbxbZ&wA?7A^sTsXm_b*h;NCo1Ln=z8zNs+%BW2&oS7m%qLDEr~_#)UFlIAy8cUA>I#?V+kfgO!arMX+idM~~_%D82kaq4seO%Hp`!GXG0<hv|ZK47S60?Ne1q$UD<?klRx_G0=<6yUUfmlGJi{+NpM_*&9zaub>-%T=TWihmc#0SZ!s2d9>~{xANS7+Nlo_`w<qL9{+MbN3_+Vs@zvJ)sHsH+bc+m(kq6@J*ea?$oM%D@+YZ{4@8FO+%SbB<i#0=Q@D^)iLyh{{Hqm?Pp0_JL^4L#dp)V+FXMiL~0zED41haBLbwHtth0cmC1@b9Zb`I*eo8K?P-m1f~7R;K?D6qgG0`{kmLpqVENTeyS+sb*BPP{@Xf22nmiy_tRcszGOT{Xo{XGhb^|HTwu;ivl|n~=-(s}kRmUnTH{ual(H%`k+7)9F9cC=5;6H7<6q_&eo5BFBt1qeN6VPKj0uw@kNMSgZq>P>Zd$vz<8r;hxX!Ly`FhrU!;Qo3QrQmm2gpAGt+1UoTO^o(Su3*2mz8YS59s2MKXb7-I7s|@COU)Rmiw);GZ!&9Qx|?ZF=LOI)@&><h>^Hi+HeJz#sP%ZQr>2Mj3ITy^a_K9t=9{8=y_dfzPvLF=u5L<SzaksbtbO)5s$b_>(lw`V=VoZ`ohx9bCn6o_%<_yKy#606J8>`+O;23NHUs8QvWZqy^vbtsSf1iEt*KY+ek{0#G+Xs+4H5uX*NHXij50<1)kk2zJg{MD!-|v?^04E5o<k|9>A5coD_n0+6||lL{4AF|h0R3&e6X1Sv=4HloJ_ppA+G5AcrW{+Yf<PZ>tHa++@x-hhcTwIkGkt-nkxpf;X<`j7R_+Br#K%l$|6o1K>mdJs6Kg7ik?FwC2L)0lJ~PHid}0Qk*-e3BA{-1y50HeS#1M3bl1|b$3x4iFvA%>S;d}Wqwo!gR5K4ZX}X-T)DW>WfIbV7LTYM|@e@$S6&-9ifvm}PFg!Pv1>1b@ik?&%sr7vGzOZB}M2s58Co=gq!XBjy)If7#;d@wROjj*Z?ZPE~Oh#mUmZR9Xf;rI=mZQ_`x=0PB_cm``dqyS8ZQI>VGW!wGN=AgrNY%~dx2xtDkafSkH2nqDo!uF9WCsLA300=^kZX`{)Y#U?-ognlUn0k*%gW+A*^9)BHgE!l6dMIsKnZ1c)OgfIx(mUw*T5}AsiH{!r2Z9y8{5N@JY_U#q}@R~ilGgcAWxWiPJ$K}8Zl`*{VvcVon;-UkV633qoXtvKdM^Qtf84!>0#MhbTwADUP=#(8W**Ps}TF3F(R9e5Lu?vJJmYk$adh?=lVKusua!gz4kLt30Yl){!`u{4sJ5|F&GRXgY2Ckp9!n<?oVnV&@V>nE`#J=-EWS243j<1fV6Vo84kVV1kYQktsR8JaErznvN+5;%Cwksn2FC%me>`fw0&SKg@%1G;?3D8*re4Bt1vI5I8poQ15cM7V}mqg!B92l(&{BQ(zm+iD7M3wnE-Y>BH80%c}hMh>-fDh7`n=ZnvQoC|JxBroVGFi0K2+r{1h~pfk`!)U{*+4b%hP7!67zO_omK2u_R^iQm<Un>Q`nEr_;$nOt6+|VpjHD(r}_uCC4`(a4Hig--oqDt>$9GNtrY%A5M=NV929j!S5c^Hg`VL%TeUn!hkyI^-gHHNHkz^3s_r!!kPL#6c(X8`2He+Texy916m8`m<l<wV+Ki1BdORoJ!7W5oz5d?lgGecT;Krf1GlQ!ysP{;Wq?8TVpwO9rtTs7#7~pqH0>Ez%1R?`gc7zu&4X)8`?8vh-3WPRUD6;{QzxA;lB%0n()X?_Wd7!5el@8vdo;KKpPvh`CgQ1mq25Y+rxJh}+4REx@A*;?Xq$L*ZN(FO1;#o(?uEZTDzFbeoBfX$-F+b5M2iwrKlhpqRM_g6l}x)2>mgWXSAx_M>)LCNr})+&ng0s6;<H{7PDZCqmY<<?7_ABQpjAh@7NkxomNWPQeEr9q-$j$hWSc^wR9idXfRh2!%F}AqltU^9?zJ@Qzo(x287;=cy1itu2@IkjG<1!gg%RZ5>xwmurPXIV94Y5Zv<vh(gXhV)3o&VkOFOY8khZNOR~JA#MzqVE4)sRe_aF(+vGS7~nx~TkdQ;1+P$0Y4uj0nlFA1rgUq|0D@%?_^&PZEz1#-W%<`aF#<_2$i<$1!gDmnfdkSKv{!ju{RwKi4ECrNZh3U{Nx=JJ!R9j*-I@jH-A_(ZIs@<Xz|Cn6J;LO+tBQe9nvTNmN;i6|@bFVGJJZ>oTvHV%&p9<JCk^}u;-^a?6>nRO?{`y47ou2#db@flz4JTWm^8I~7ENdi5=<Hhm@*fM1HrwOf%jg3B}Y#(Uj2x|g#6!q9l6tDWjF4JUi6&C7A;Xpx*G=QODVULEBO`+`&=qzSmjm{Y_slRpG<Y!SdJxp%IGlZi6#pKCkG9JTm1z8kmwFKv{G7*-O<;0PGJX+nEI%g`D<ju%(-&Z|hCggzQRsCRx{B78x58Ct1KA)$1706*>1;_HRP(N;0Rg3E(%JI$agsEZcCqn_dl&6#=N)u3X0=Nk?OQ7cFtwvqTOS<M5U~MqPam*T@vL)#6n1mPk%-0mO-Q!!T;$1(6oq)v;oaxt^=(#1cvP(GRwIMh_Ik~*^&m^mX;jH!nLlB3opv5esFM5m!?vzLS@qu_vk!*7%lKC+!YQmB<m$zil3aR@zMJ$9(tgIBw4@yoV=Wne&%=_e7DmPI*ZTejNhOA&)7&`GNzdBMqvU?OmDr(usS+_>A%Y#gnOi^YvxN^U7FZy%UBnU7<GwtZmeTT>Olku@7O-zE~I?TK2sZ}_pyRoIEp?1(4Y=V&C6dD~gyVN#ja2{0laUn6^Xoz#yP26s^xY3P&WZK};AmA?O<`8p6@0CG0JsgJWrXc^tNOQ7x_LSJB*P#(c=uyy9ODxQUvYWUmX-F}<PPR@Y;XRvIU?JLBF*3<p(vFY$a>aphpjnzVKxiYDS#_(nMXY)e3`S3J&CT0rx*Yoi-TA~VV;na#;grlXr2UQaU2|O0CrS$)>Db9^EhiHYqq_S;e<K6V$EF{n&yZn%J@UMb+vi@95#D_VDzp>@ndmhMn){giyLdLlgaPxFFFLRwcMg80V8g1mUbGaAD3-{f4G{Xn0^E2a$FY>SMSsL(&n?o(H{geeCoctyp|V(4WQ5&@Bn*~;IZZye$`D6-vT*yTcf;^JS42wJL)}qE^((%PicTDnLxY+&vV=(C?*3r6-&NrTB_KoBmE8TA4tXd#(J4M;$!mgMQQDVp7bx$eoGM|2*bCzJW4kf8fX~JIsDz9~9}uoVlRZ4^NP68a++7!%J;nRL>!=KMg{Yz!WWVO7$>=Fyf|le7g1?|T=8%bSV9f+%XOYCIXO%)~`fHTQAaH)VrHCm8<Z<aP$)@YDXKUU;hV+8E;H&v2qD|wyn83XCRdA0JrQ8IPa*V(L7VPRAUb~%#@BjXWwOI*F6`*i9sGw9Py%KQK@pDs@R#%8wTZP1Tj%^kGm=%hGkk~{1)4TAtP`sTH^|$GVaJz!9u(lL-Fh3ggER%V$B9e+GRP|*pfLBu4<MV+9&K%Q=y!!-ClJsL-XaH6T2$_7`>GtrU{1P67fl0UeGX(l)MA6GD{BXpz)>u=pv-q|dx?Q+}CWFd)M4DUc=O<4JXE7fHCNsdU+9!78!(562&ei33m;0`fI4W};%Pfv|nD?PbW<Ddoy5aXVj~peYtPr{_VLSd?`vr?e3j6xJ7uiEM!atiMQ%-*`&R-9Zc^(_&ONI}Mx68$SmBSVNknWiw(4xIML%mLAm~H0+ieDY*IdntjhE8-=Z3EAiRpCIV@6IF2#t81Cu;+s`LdMaGD!p))8uf9VuW^(AAGjqj(ohb6tC*s}z>iU{FUcpL7EQ}UU3&$8n22j-zP=qxd6;l+i^L0u))Np=t^qHm7S!Q{4QEE7j4!>9Dpub*FmtqhmYe_qB<p<@qs0CTp9f6Gx|20Xe(wVO1+;*3B__#a1i@nicRY<VG90=5v&2ev=@QO4uUwXq2~{+6l0xl_4k8Fw?SH8MUB$AdnJY!VA&ma(kuFxCcIl<;P~SRuZ$dvS|E{QDp}l|ppeUKfz(}032_}5n@H`vwBP>lYhR8;pHc>g6+=h}{M?E7bGV#QgnfyOlk^~*VwkX1d_A9_+=@NXAYRve5R5$CL3<EJRSn{VNP1DDEfuiJsty087FI(o)*f%>`8V-i#cXZO^3<N5~lAz^60~)Xa7zQcjYd-K@f8AMXGQTkABvCsN8p9yDN4<}>7vSq@Om0c$O*(IpkdBWqD9&PwBf9bf-Jc>#t#hB5T`%?Puj}8VH%}7i0w4FP?>lIzZ9|TS{`N%drk6NCYn*L<eI@+=+oFoJ+$nc)ig_GX`&1*b2M>h<a#yVPLDZ^`?=4P<@j#!`G>vZxV<aofT%6UetU)^hQZ>g;Oz`aaAQI0{HCmhgEo1D!33w-0BdG{ukt#d(mF&T>MTkj1ZD~@eIGtz`tG$cngHjrb5boSRNbBsrKygKgK7)Bi<m|I|QG~(rg<2zGx!gm&X=z%%hVyl+MxH$H<mg%>t;fZ$$A*SReKSng`quuh`W3{|BkGK*wkd4#0Cs}x1ed>!Br^V)$+AQQX|jJTs;}^|SOT5QTCn8uYfm$Xxm4ov8&MO0EkNL@!0G3qdVxOZtx*@#cd?dg?8cH)y4`ZI?~`|Nf6k+E(s%I_K`)oBmwNI93vEofnfXGuxYeGG2EM$2xRB8s6j^pM7lJ!oHA0?bo7?HZTlh-=9?e<o2$EMx-We?<qHj<xYTN%nN7r+I;&7QCxeA?Ku@JcxKLLf0DHW--Z>V9QB^fk6(v&pNxZY|zVyjYAB+<T(5@lfooM5Lya9^9AA!bBFEJ*{<84O<|lm#1)_f)w=^-u}Pe2Oy#ScOli6$)~oBxlq|e2fVev!iu9BN1Cq>*`ddv7&>N&)EC%2X((gBlv7rpVt942B2=2wA7~Kp341jB8rhdcDL{1tUC9^UFtrb75Rvt7Ok$G;;dkcSiSl4B?76~FjG(i4T~4W`}$;uBz4|qhx-NdySP4y)CQDr1avBc%!5dK%mRW*OVJ#^v8pi%?f6=7B@X37Ss6y#Cbyzh{q{2&{??21n^D_4xb#wksW;o4W?ngwP?uNNLYIb<Y9@HVr{$x@kHxL2git$De^t^!31!+iS0ro%S`f<ooZVj-bTcV%p5}};zW%-9d8)EowAHLTO^p*diZe0U%IL+GY*XR%$)c~i^xpuEchB|hG%6TnArLRLMmof(^MCjKxNJtrGNB|K-&Zd^8*OkypxzA{=5;n-t^M1?<ny%q6huGoK=Ax!!hvE}Y&n8Nu?Dti_ezdL>XgApwoHc@@vM$la(Po=&0?0@fnd)q0dl#xwQ-^xqWEY8a1%aZU=CBs&&{p}!?*tPZ7I<zOrD|V%nEyFda)o8!(jFitaQ6bohIX2fA~v&p}<~oB__J@fgYC5mknyzx!uWkj?%d`4<kKE<3j`9D1wr3*`F|AVMYi?kiy$68NPrpKh@d_@8DOsAXeRU?g#&R;cjEJ+oQbcL~rO3T85rC$5d3pBss-);IC#qTUd)d^0&X12Z3%^H+j@RLLyGzP^Ml|58YyTeVAe;{2w21#nrCZrF5#Wu@z2lwnV;3=FnMGXIW0j|HE`L$AULO;NM)y%24~0K%FgXZC|(*tg{K+Tb-6B#+NB5CJgd~7K!Dgswk(94FJXE0v#lkG|T&B>TfzZA2|WuHNd9kSOdhbZ2n2XTfT~gV4Lr2ol}LbJQmt?o*0_$q|`6SEd#nH$si?O<uBhRtI}-FlwJe9f~O0>1(tr_0p?|JeMhY@hwN>*e`cQ7WL20G&Q0scO>#-1(?&r;J6}?-yOfg`Avp;DeQ2UdphUoE2a6q1IAqm>p*ykic>s9H;JbMo^VK@H{4_ktJMU?LdALfI3ou7spBWo*T3f(&^#8Or^gVvpD;+&Pz$0VVqsTfd&=8q;gde$)w?b)vIiS7U<MX;}N&6ro1kBZzsfmk}d<f9598cFM-n!m^JAXOd6l-9wt&A`i)G>T#8joC86p+#4xF+G@J8iu9F&?D@&9~%w2~{MOV>NA4FEoKLqPUwAR7uO}zS~ECVPvm2esy#o(55Ns<+*ua^>Ys^Ieu!FrmcF>(@2Oioh&HnHSu#8PrE*mqt<K8+zDoEe-BqbTPM1Kn`IfJ<b^SM%^U%7MzQ!hRDI60Iqm@I+-8=H+l6RwO!mV1j{%F82qt4xRc!6PG<tA#8P@C=`gZ2RLlwCl=|}w?4{6rQqpL9+;!WmIGdt21d~JE5<Nf0K7hU)gNM7>eXKjGl_Ve!1&KB&sO^xVs)2w9_yd0hC?%mAswB-pfOFL5GM9qq_(=IzeS)r2%O?lj_HU`{0@TBR;Op0kIJoNt!IImI*EHz(gKDaFecm_xo#A_PmKzw@>Q-JA0ehH>ckvsAN`#G+kR@j(i($*goBW3B6U@Iu9#cTl^*V$2t@Q7px@ArGSe9$ppgoshj3Qb@)zMZhNnDY2HHRNdkM0>Y!mTiL_PjJwJXPY%)qw_uBUX|PE-y&IVCWx6TBp`cxsc4#1K{2X!Ld@a}1VSY>)4CP`Qq@6bhM53DGC9m*8#n~61Rp`-Me$YJ;arlzEyMI0@Ge9*PrO7vRm6|?2Bi$h88b3dbz;IuIWCtR$eoAV(Nq-uGHJ=mxxY49KwD4x#3p?_mx|c{#5dMLq^bbEBL=O>c-IUa;2oyQ8VYwtl}0d&5&6kejp{-`c1w)iM2514p(%%HV~hP{N~Rl>&F>&@wORncizo(3Yi@IG)B2=fDBSWt&?kzm0|Do|B>?hPS<ZAE%iy8|^Y{@<|6cp7?gj^68S_SmO-5me)jlav87<JD*s@#d9D>Z`TUeceT4EL3O^Y^5dMWgzETzbK+#L~%d#nA3XxTz^15u7cigi*QPEcjG3n-Hkt0LJEoFK`$pS!cV*`RSdgM6YBXrabJErK8|jOE@t2;2{Fa$erw-`I*PlOI#nr6z0Lmo)fzSNS)f)B#*U+@N}{1Kg+)j01=K+_tj-YB=yoB^@m!(Zcv~blV6!XT1(}f{zK46f)-H6u|}58nV}u?Nb>PV{04F&m7TF@R`4)At@iCWN^?L%QbJP`|xgTFDg4tNuJvHH~Va&OMmt5IP;|8_P~qDXAYcc1@pxS0n%WPY!{2vGNjK3I^a0z5<fE`u@Qh~Gm^a@K7-?8Z2r!Q<F;6EJ<EE{02^_0Lk)yUkvE4641Vg+x9__Dv&vXs{<rz7b%1Q|JVw>16T9Jg#B{kZPmCF>(Hh@JQqln2G6qMFSb$o(3mmh?9;2^;(xXD$n<5mbShse=*wkjIus&ua1={$>JO;z`@DFoL1zM$rH1sO6NBrJrJTe_3AqZ|uo|6{(bZI^)Qd}JhKH4_eP82-%*?DggRQwL>SGy8cJ{h#4IcpK-f|M*qw94JZoc;hjUP{|=0%sy;{Z=e-rM3Kxx?7|}*5L-<_J%{|*r0MWe>I(EA}4aJV(V>bMeaVm5}-z0M~$~^TkKa=`Y4d;BWVbHfM1v#82_55yCMs^E<|%Zsv<lO7`u~NjI;L@V6A)alADCoC&B~)anpRJ<fTyfD=>fmtz=K*XkRW|$=UH+MkMuQD*NY7wduv?0{Z=73pkBkbq~MCCihV%TnPeA_9~7xr@fAauJSpT4%Ij<DS^lRQ(k9!FNL{wM3xw4200{h7?IyRv>Tx1ThgJ>u_PX_@B%MH;bwTZw^rp<Zp{~FyxfvcX$fce*RV`1G@f}1R@CR0O1er5P_qKsM6A_@%mFxUyP;RtmCU~rbYrdz6?R+As*6G$+kqP$gBuipsdh35LP3yGzG)g*Q#`J<;cf2y%R;6r3}j=XXqCtivFpFCE}HgyIj|JFjZVsf(I_z)qv{?dIv4KMO(7t`$WOIZ>09;g>DO|oC^3za-W&%ers`Fi;*;g@NPUY`vu1$tzh87qPrH{G&C#SP35ls}Euba9h*U)6t~9bnQztP7&$tYK#P;z%{PC=XNaNzN3L%s_t{QaV+(sRbG$Q#@lJvn4xq0*JSAUD!IBFy+6jHrDEkkn%-13`4K#FD(XPbxMBEa7>(b^IktWKcOmaZ=g@DuHV`2NY)saAd$TNzX=f_}MLWj#SFT;7RGL)(;K;h58)_4<<_UJ=6WjG<vl7jl|814<6p5cG+ixcW<=WqolC9DNWalG!+YG%F5_K~Mp1pojL`_IEloDrRin$(nrG76}xGNE~iNhvBN#?N~=raTwAdIxYGmi&AXHDY%DQw=FjA@PUxa`^`~(KxcYXLH!qk3727*vx$@YnKj_t6Jxm?I6#v<&nP`+3`nSwW)c83`og^`oy9G!xmCQhi5=Y%%A<%VR5!-wF}#j{+?e`<<?(4y`1wlfsE8(~yC=5s_xF2I4&}5f=O!$&{tA&vOWT6B3KS(Zjlgl%uMRJRjT0;g2XWpT?Yn!Iqj5yU+DDse>1x%wm(dw8`m}O`yjnbBwP^<#ZX^Wx@V%3Uw5-hG5O*uzYIhWOT&u4fhz?MYjF=o1Jy*JtaE(u1&K`_&)TcA%vgb{9W%tmiiCpG@wz*RJ-oF30STe9<>L6V>ai@v*HsGFVAC>^DEa}X(C1g)bB`S3|$2mH09?V-h0w`WJ1j?x5ISuE&@K^ophf&<Bu{gf&vqo|9a$5Hn^KUZ+Vi|jlG1-LnSba_F<4e?MmngyTPT5R^Hi=oQw&JC;Qss7PW!0(snyQKqaMn4eDbhbj7yRp~N6#!)TX3_-JmEq*w2s?D`U@c?#&r=%32b)e*k>lj7<|K9Qr)%>&_Tw`Qk3lsBdAJzK&ZiZcovo$+ne(+eo7~o6z!&%|4>H<c*e&yZu`h#83evFcj9=n7*}al+__f%3FAo)TBE^+L9t>(v-nu+CHrR8O9T7}>ocRNwYZNYK=`ZMlD1gnu>U_UoRq2g8^OSW6y3@zs2}UOsurBjTDaHLZyq?XMS;a@jp*;gfSN=wz=mczaM0LH4Qg6Gj%0E0mj4;J!>S7t$!)F>F~$NEz)(~=kac;3>23;-M*U~z<TQ&x{LK><iLOMrZBf65h=anos<c#rtyfzBTP~Y)3)LI4rD9NQmOJxL62;CuR~~JxUBY|D%kpl1SA2MS>{X<mQQ6Ez?fq#q$+WJy*w^9Gk}sD)!SD@~TApaHyLFZ)_*KNIIKY$+^(!u;C5(!od6`dZMr$u7K-Xuzx;V3-kBMCdt4Ro>CoZn76zIv^_Tgs?>mQ7b@X;XRV*OT^Q#jB3@5p%G!6>YrL&HH!+oRG^#<xHvxS6kO5i9PXLp}#HtW1dQN=&}W`7S7<E$Vo(L!@p6M6F}9op5iJIhvw@x0a2^{tgTrRz&+d{~86@Fw;db5Mn`ndQ>^{`x0YJNCDf0{4Cksl$kzy7w)qyL;o|si*1{?8TO^xNkrM}<6rzpODZiwRB6z250e|WBp63@AjjQT$eli+P+_y)=UgOSN#^kyJ8V{|r!ue1J;QjQgQ=p~yP;r;&AB4lLTtdzz^R)Z?=j1Xj(tmTdy`;vs4Fw$zL2*2D3C*Qte}Mu7pf;0>?ULa+Ry_=WT*$Cw}AAboYYC$rqAT1ZUOLirow+t&X#lQxUzc1wok0?ty2o(Mgv&RqvzoOFUU~wpNQwx2bhjv3cZ)ggTGu0%D9rJxY_I4hD2<s2D74^*1)C96D-%3q9^8Vm(-ykXbw!V#=9qn_6=*65y<a2e;OpXB0BI9{LVqhtkFU-)9u^&0uy&8j&MwBLyHp-X!(cMmqX(s0o*i_h#6KUUz*XH9s<eXY$_plLss)Eoy!?$0D!pxHwqX^&GzJl-iw%{)ivF{78bRqc9MFG<en6YkRvGR6)#cR5fh^LkR)R~?kG*I-AEcw=^Kd-_$_U<r+}1;#u+N?13<jxj!1HRL;}$9nHG`n7_HBJKD=51XHmpAk~<$}*h=bwXips<uS94tQI^GD|B&}t56ut?jfyojVKK~$l!S=O(jcPaJKJd!3Y5;K>Sj?GOcx<0^zV1CC)nhRuIs0$-8w>l8_2{>L_!hr*P(}lZC0$Hs61uSoW>b>SP|_@UtH$$WJF$C<0ILQrX5W=K6E7B)MRBq05`Gmzoo`a`O_{iZ+<kx^w%LJz^y8ChThgx2!e7@V3-Q(nn^}W4hJhEfl-C40@hHOadrJJ0&6KjVLRi2`=?V;WooqyXt%e$i}9vygIIY)k4|WRCT1hOD?T$7PrT_4u9<2oqi{vYbejq19V+c)B7(>zq@fzm4mvg`Kr?}OEDT#>+HG$@RESEc|1GG0Ep6W<Xef>!8+$QE#e4^Pt06AHKmI+3q^+c!YlKym;X@@(uQf>_$iscg`cNI1>CHPh+Kr~Lmm}0<_32#>qn*rFZ+K_NZf2xcZvPi1->u4<!RR3^lT3%1iPK=<IQ`b+xn*3%olTFH0LccN1}X;EB)qIbk)lNeW1~sUkE-F5VY`j3_Qto0uZ)}s>ug^>>j+>8RC2tEpEbnu>2oTFem_IB-`Ljq+~hScYaG6OeS(<7VK)q;O>-`-BYJUFbdLL|LDCFy<by(ywS?pz+!ya@kd$ZHRC>b)!Z2oi=9iBlnO8kJAJY_cTr?IIAddblLWU)iiIQdAV3v;wpR&8fA5hD%baee4jdOnbd}KP(vJ%=w=!x!M%|im0jgd&6ZDSfJ!UHSbH!l8xfO*EEgk|*(B$HF2TRQJRuLOxZ-RV>$+?GMHqhICMA{M9{&vLj0Li<iWJNj*wy!(Pax0MZdgm4X!_VsOPm;<O}^x<S(K9=`b_+9Yg$;w(`YR>fi4B&f>BC;bFJ5jREloqBl>Og_w1XS2X46GG)a<hGC+{cHvp{;!v!0Tou%dB>hAj%C!G9XOzyU>`GV8u2hA|tCCV0q_;D@`os3epxPuEe6i`R_H1DwP!e&kHl#jl!Ji)>NCpq9Zh};$6!V|5Le^<Nh60@EfE5I2vj)+>`?_XH%N#x!YI{LT;)2v&1D18kBWuG`~1)1MT$s1}i$Pt`$HMb^1pA`m)$TR|?m!Zo;aubHe>6@p3b>JUx|_8@5d4NE*suh%T&+ra?f)vbFHsd<GvLr9uGwVo-x)KR1m=N7)J>t~u+%b(>XE-{E&AnV&MQ23P?z3(`gbYfrJIftPZ99fb-peDVVu@uu@><I&=?ah<U5j|F%YfV*<X%uJUmA!Y=Rm^(hBUw(O8=lTPKcX?EPHW*2wA%NQuZ5#ri9R|IZw?mwo-tvUDP#Abm$%CFbjttoDIjwGIEG(B`0{B1)<4ipe)J^<2zj8K8V+`FsB4)Aq*!kB;6M~t$ASgGQ1Fm<3pP;kQsOZK%=b^^@6sLlW7<S}{Ej3pH!mEBL5AN(=*DJHCNnq}iegtM{Z)L+Egx>?_K*cKM{d|D~``z7T;@Xoms~0=N7=3sUwZ#2Fh8j`kuP0D8MA#@a-~TgQyD6+1Z>|bnF6ViuO4zJ~P^VkUqGSUh#YW%H?A8^+zZ4Lxkl5RCwpF?jwc`w&{idK0nvyr{H1_56>3&PfwP&I9OFh|zF1A{DCZcCXyG7Iw!YJG*<$?Pz59Zw02jn++<*b3f^=(5MO?LL|SkT?j3Oabx{LA($qI?ADMl{)YUIBsG;%_A`#_ZHCj`c8<4t_Dk;(_5r%eLpI0b@4u$vBS$s+rQDZe8^9V7<n>6^Ythbifbn%_Do-Qzck;@mQ_rY`tr$2rYpKT!ZDU2P?Lo$Zt84eB3txskmQ`rRC+L(Iu-o`WLvcv>SdcPeLjo-J=$IkT&0Cu`_WFsKjCRryNhv%mPm`jvle@T@?K5;u$s!@uYu)Pm$;m#Br*$B1e7YHvYR%^8Lb2hzDZ2OwD)zEp<B+rHBKiDhU>*9K&U&K!Dys+InK(nh)#g`uy_Gg?exvoQeh6aU<t+?qA60wkZLb@C@ezM(DYEsFZX4V&W#@DT#j)m!zpgtZl<eA#y8q#<;H}PJ}DPu>C;DQc^w>vYRu37df-(rJI#u=9UgeUw_4GGZ%L@zbch$WaiU^C_F35S5pj^B5$F@%|8WKwV|#jinPsEfg8ISfkp^{3s}^qmsjdy#wQ9r?r-7!r7I7Cm?IpgH2oJyvy7Z5a*>vsAvgJTppycTSXOgJh5wjnONiMBfDD5h6OrCq4a9$h-Z8@4qGWS#7zgkY$yR)BM<L_ZWy6Q46M?b84?fKT&XP?~aBVB4k&J1!x@uj%R?!UU4I#ibuppkyOuqt}Olqgub|VIoF*>YTBOdHo!TijK6TDlm-YC|KSIT#nD6BSSKOld!*z5u{fV$j9mrj<E6cIDgvK{HcjO^|@AQA2H+M*zE8HKp4la&3C3-S`Z<)b!=vMp10dlf5>;7{FO1;4(1fEa&u=w5&V(gII?MZk7Q$t3qX-fWnrYLQ=bL4j0^T*oc!WDI2pxIn|+Tz)9rsESzB4oaAbkR$K(K$0R1X;I%rW<3j8KpCwH-NnjDH{C@%rv{lh7mapkk>%b(nvb%=sv;YRg!V&md>zF?;u`qmBFK{T?m|`WEUCGc{da@CyScvq7`<uJHoc#>SR=dkP4G;6RV`DF-t!`bG7WSWGuB-|7Y2KvW^W~|9@L^u3I@o|70hkvlKD-aj^9G@*Y2Lj)<~|!zezR~R)an7*v9&$?zHl1J4Td37U<c}6l&Jb+nGnKXDpVW33f<0)C)ujFH@R%eiQK_+N~oDv%R^`w(x=Kz-K=mZy~|D*&{j1zQ!%sG#TImRPa%XK7i8fgU}g2J!OsT2v<6}1Bo8l1V>Upcp5rWF+#!?=dsanw0Vy)oA;S;3fTSn8|rW|u>VAbbwb1&2ZMud=JlX2)0?+>@&c+nX-$3qLRFjdK2q3yCouv>ozVP3(B(RoJ~3Xr+^jE(hFCdvogO1@+QZ=&-{k{qcFtG1SnAvvxXWDa!z_?bvvO-q+q``L9b+f~hfvhOk9kHps)kW5!y`$KWiq_l|6P{z`v$^J9!vWFYoy|w;hOy5jjmm|sSOo>9U{`YE%u!a7ormfzm!q~<Y&Z<XjaR?ZV*?a^|AOIsyEeFwLKL6e?@#4R9xHoV9L*CUSUte)?p_xLRk*+Bs+Cv>@>W6jqRJOx+8b6WHRf_Z7G!)&bFb83!r)UUU$mR*e-3@WsXrn(ve54JRqY>s+#`M%)v$q{_k~9foAJ3aNS`Hq-T{-$kEsb_p*@gdC$Lz?N}9eL@=KWflor(wD)XdsWikf(0YDVYx*(0tYK^zmGC;biMlQQ;dKu&E4bA_LV&vNcDtUOT0xAzeYKCy5&y{`eRs8LoBoaHUb3iKi<TAogk&X@Fil}8Q$KmG<y4<+m{x$KGk_refyqhhV?u_h5GdiZRIUr~v*5Camjlur+XT0c8r{~Sf5?4^;YV0PH|TkYL$YIcdbw;m6c9gv7sss@Zz<S3QxrAI<O3eL^MOTT!6YKUV#Pc_h7RQ@X)-qs)sjwCl!K$>Ump5SU`;fINFq5~qfm{)4#sBMeS(QMVttLzSA(9g2BU-xWa&!z5$udc3qSt|_Ws}39P^n6Xp&fo<{|yO2K0^V)~6ct?S-CIb7Y{|py+`fs}Dd-iG*U6OuY#2p-|zvk?_H0p?CZbhfhYFQWurC%GJ&Cw&U1EW(c?>qjj2&AUKFC8*1UWDsAe+M9#P4QK#8mzUMy!19mV0m3_)ETw74SfCpDM*=+!as0pLBC7STZ7D&-^kgF4GlyK}ZCGHj9i<W@~b4t1_`(O;vnO$o>bH?yH@XY`~58>f3*f*0n;NcF$YY7ZkOTyE9#gahH;Ct-Jr7vqx4?Fk>KbG+rlWi_Q@v<F3Y8*gK2jD$w{saNLQy|y98uq`;n%$|nqo&$rHu(7ccj;Pj#k=&+fx01U!s6Q!t5$6Q$j#k>6aH^{W#}k1i?0Sqe~fjz=QVr$OwR0>|KDPzK#Msv8!HLrliQ1?h^=JA=7dk8$F#~M&10Lg#-1z$6*uk8*ngcahRN{J=uX$JlnBQbQ3<__ZsEr}>l(+r47U1G4MC;EEp(W=aHdm=*6nBwhc}pUJV=<&5&iF!m9K!<*IL9d5*QcDZ*I!HP!`89$`~Zr&2}%jl8EsTQqvyNK{g1XPFFp^H<O;84*BO63_j8zCHZcxEzAv*4+W6(A+~ktmb1pr_jL;kKEwoI!U&B?UB57Kb&=?<AP?xZSa*m!qB=(nqPwYnjA@zO)L`Lk#u4kQ=BU+S^vF^8KYZp=Z}yoVO$T3%Mv_SuK@g<YK2AF!VRg2v3N(t=a;iO4T07N6jrspgOW|fPO*T)`%pHH%E6+P?tXgsfK=!1Twq$om!@lw%IhN0MnscQCL1nw9xXzrTE{G3|Z0-l-mN}$Qi$o-dKs0ZmCpbJ9KNJH-kzc8#kQsN0ZXcEZ*j3m)aoh%b-;sLTA`k`%-mEGV{{?>P8A2iO#5jSf!xYDw+q808XC@|=)}<3`_1OuEf|<VgSw&T-HEn}p3|(5ds>(7z2ZNBS3vM2v$YRtO%W?SrOg6-*_ujUWws_<iyA?AvtRkN;8>TS2s5ob3R=zC8?mX$~v=i(M(l8RFhXIQpIjikL{%FlZ#x}hDQIU&#JF?I6S1HjHoN=WD_oP=404%GbdfYzcgAJSVCI7io@43)Q3rVp|lCd>Am#W$bQH?3OGA4!nP^y#jyu;V|vx0gz*#8F{Wf6xxQz8cY%xzN5?N9rf1u#lY;|Ik(%WnH<vui_=4FH-|#_PY$#e4WvHD&9h&$HrZ${8fObfR8DJ{P`{;9t2K)%Q*jRPX3!a?Y~nKPmhCHl{sVs=DjbM+T<%-bh;Q+&z)18)-o3a9Rp0YwHO<nj6*BJm$4bpjTvLZmgrff1kM<2oR6>@VG`v<2VJ>wzyWYe@5jAZY_ljF-fh&GR-KjrP)R?v>({YsQjW4IJo+V7<xWktV{>Cf_j}Z<7qe)ABG93XQ}hL(|P1P=x$dWY@BbWqw7<1k)ti4ic8UktnPAepuZ0sck`HRLt)O6|AN_m0j?cl_S_q`#pAI%#Kz=?jy^nsu`_fQbvm6%)emSOH_o;cK%##9?Be;2k(**{bTXtAQU?Mlq{KNp)Pfbl)~{6a6fm~spX(TM5|5(5X)5Z+e{o(Yt2c4{HUu94RX3p+RIvRv{aKj;7y2BrTq1nEV`J+nLpC~KGwQ}JGoKO$k!}^OFcS%5tqt#MW#5lL;zLtafqGvsZW2Ji9)20Q;(Fp$#meQR!o!my3bST%;IhlM>I*23Q*w11&11R^LJa&jR;YE57S{w1nmP-<05`)V8{o#r(b0sgMT^*buTswK27CmA+|lr5eg>6~k{zpR$X+smW9D|W5CYvi^_(wxFVvO7auAn}74e0?Klk#FNEwfWrAOhw2vdl~P*q$%X5!SRbQ+&qd$kq=<wIs#1322iJAq;I+qB3-`g(9rFLbc{p%gmJd!!R^VP%oeIE!;d0f#3AKP@+&ZnGv?&Cnq=M2RA&aVPi>CIdk#`cvp&_xw~6&hjMKsAl)D2>~k=zm^q)MzLw<r_Uk>zRsd`as)LjQ1f@bJ;kDVQcG5Kj0P8!Dc5xYf}D1LeUfOv?rBh;lFq5+*WENxWkv&2?N!Y-ykW#iGk*0iHkr#YH4BJF<{TUG9^$(7fij=J0QkqnxWa(HE4uDW{Yq=NG_kEyeGPQNPvGe}nC0|0T2jVqN24G7i~aQRKYB+`nE8dO`t+3ch~~#*yt}mh>Upt3<e2UMDtqJ%ab+;KD{Fvo{$BOpGj|dv&`pYQtJk`JQmb6|Jsy6WMo>^AUEzZ|#Atz}VFyJxbxL?<X)>1}9Om}Kq<HtzZK3qcM;>*mN;by}56hk-nXP+WVxIwsB>v+udfh#TGUn4u$?qbI4fdh27L!68pm8;9lB{T`-}UGAnt01%%?8Vfc&!>SEjCwoCpPX9L8sHBNTVk+5oJBqqdVf}#CbKn3TNQ6OGf#$7VeMM4I1{S-~Cw=YS-1>7ChE{U*r682_M+-fCp5yukBI3q~7)!^aUJvQ%duu<_jXx;oY<)bwpJ5vt-Ajtv$M0evUN1H}BCO;849eJ77e%-eh8XA&$o6Tb2<$`dkVea%APo^x_>vO6?JpaG?NkAK+s2tQH(Gc>=su-4FI=2woGDWAH2Rwzazm5VDg6Mj?MO*y`Vu8&gn3kM<?`-*msY08^x}V*bx`P>aDnOY-$<PG0|!=X!!0*8T938Tn$naevkrCO<fLS=A^OiXSh)k$xDGF{A5iOguKzm15!Z!xxGb0BMAn2n_3W?nl<t@l;?!Kded1m0zWj=N`8Y#t4*i%mHT_s$q$c&)UY#DS!LK|3|EH8f+-%B+touL-|HG<%mz<V&*F8;6KI1T$c(2Chpoz+Ya&lJtDD|XMws|(5Fo~#;p9KvO*=-EDO&lWx2y!BH!kEuGa7dng)hhe>B=>y|zZ1^+N^5Od%T#%jlliBlW?>z`WUk&hMZURxd$RW*?gXouo}S*&_xl@G9rB%Xo)B=QYB@PfYdRVw1gaVdwBtDaj~OVmon4@gKqBW5%cB3iESaM0#N5VlJqKmQFZNP~wDNR=P@C%*lc)Y}5A;6|cr6j_S(;WK19-|CAu>h}#v<cw?~C<m5?(<)^JVs-R^KLY|*vqhiMJ4JWGC5D5Mo^3ut!nlmn7JGKi1?U{f;`H&mJjp-}}XxX}@dJ^GhTwG%aDzq{|6ZRn5*TbxV^pyVoGa+SMOO-Y0l<cNcvSO4>Q^FXBPYARnrYrV<TH-FBLunLv5NmA-r1*t6k*VVqWVMyo5ecXn5+#fLG}H9_7U7HxdOp?@=pQwdULJ%YCFNZF4^Zky>#}`=O5hN(&}F?-!2)B@>$#DeoX{dde{+ax4|(R?4-iJ=f6tmz(4pu)F@@=v12y=+nZEMpnip9C*+v;wMRkS8`PKonmB&Dh?xz5~Bkmj@YUP^Jvmd%2&PAl}^uuTRi}nd*7Y+k_s$v`alR+|PZU=R0*pN7lAJ{WCG$#KA?}~}abiTe^0dBAA%U-|Wv?2DEMapRrj(^ts_c<cmsId0JKN$VkpaPliTT4ZVgDbYfl$2ZI8agX|DVGn!={XR!v)%=LeN=IkcV=`kE)K8)6E|6Yw^kX%?|mk=#c;w;Lv%`}B%Vp>2~OsrvKI9yv*X%+pfY)CS+1^^J|ow)I=GyVoBb4(B{^a%^T-F!=>CQmQVbi=^pXPFaCH4Co3kF*nD|~H)fdq>fdPXlo{RvxR67x>VdEgwkvw<5zkHCz8pMK*C|W&*HT>RvH3uBJ&rnDmAipH@!vLa(T?AMV1}((m9VB;Qdb}z?{efmETtqLu9Wx&qA^MRIMkVyf;d1hlGy|Q*RM9bWOCv^hYP7ppaxmnPdS2t8jlSeLRDtcAj%h^hgi|LgiC1H(ZyZur#`R-Z<|FJ-g};%(w6>UH6?KGEhoJ4Dc{5J8a|CU^RQtW~=2$#$U~I6)az_F{6pbE^>MW=$zMHe$N;ij1A!eKN56i=(xjQMysQgd>4cIK&hg#x9Zq9r0bfR6up*8RKC)5H?NpST|t1RC*?p%Kj_a^vaX0p2_8fRbC;-sLeU7}R~g6Y=+5#vvHySCJ2sv*HTz47&10L5B1&&-)?DW=4cbk>GX?%NJK`oTa;i!%cK1xY{YildC^LK)doSkPw;IR|M0c?=H1LVksW(hIUcB>hB>s^t-|1H>$^!<wRZE@3Cq!&$|_VPZaBRN*Y{y(Pa1kae#c7x%FlS`@dMwR8j8e14hrc=uW!=?7WVO5&D)iByLV7}0My$YASRu$jcjM_=t8<X>5Wb{u<1?cQraSSRL`t8;#6_@?beb0@Kc<+3cgrsjC-Qa=fl`15yGE^bsd`H@P?9!KZwY9+lY#!dG)YKwc)pu|~){2_|Rjmp^s1ZUt1`>^-I1HGek*wvR;`0`kRmP-L;JQq`pTeii-;q(YH`hLVowc4D}?X3X^pN4%ccLZAGu)j**bf<ru$Q4kB!P4!X1q)ew%f~JOCL@E+pqtSW@-7gni`zM0W<jn%nQa9rZ5_6zhPNNsvTp^K4@#<;5vDtv=d_a&yjjX67kO6H=$f*=^v^3aQEQ~%@5^vIcFk~fKNnm8XNb5=4S12Q#KASwU|kk&LeO1q#)X3Ye;IST@`hU+Fk{F)i*ZP#KC8lun6PCt><3O~n+<jZ(x$T^Jz#tzP1()YeLZ1JY`#pwLvvJW5@WS%(9=sQ&#IfdXJQ!i;b(=>`@7Jmw2Z`4r|mZ3m|e44&0c(d1z0@A10%JxL`76C@iA7iYy^_Z`>WQ0rb5q<1)9ezFJu+Y&Z4xdqPTh_F_#`1nQ2N?<kjA!5q;@OcHQ_EcZsf*+#E13lXe;{(HECQOy3YFz5Pk$b$PG7w5=+zsMC5jfUiQW9#trfkyT=m)+zrpjn~+ynuUr?NL?OhVgP3R9c)`q9|V3Cbig3Wllwl@30vb5-;mYA#$38TG)Ml+&l7TOCpB+7SRE1#T2hsP;49^Mrmd%`cszXHA-<8o#7pVA3dh-4Y@K#=SR_%G{8-UnQ@&X8u282F9N?0$m&NDZ$G#5|Zr9h_1LtmJjuUPI2Y2|q@N_(Gk;WMi<Uqr~fu!@7!hDrwYcK?JAJtP<2tt&<zxlvBQv-()mE)l=64R2(Rt36BiX{y7F&>}N9At?=_q);L!{6BnFXJo;Q8tr6xxOmtx1_Rj-p2e^8IxRvNfkGbuOgj3*UN(Lxv6>Jb5)<iDR?}G2g#nih^soQ>?h5G3o<HiuH*N0%4n1QD~2UPJ-c50S2mnZ+^<|IGa^*&iA~cD4K$ZJpMv5<r^<Vcd5a)6`_HREhR<o)o^d3K%gJE5JCK$87xi;^*E>$?o|Txw0<CQYlOODx4rz^UYnOch)?DM;v?TiFjeOrS9M@(-h(>M`WNVi-<jgkWQ5UqKb*y7;C<du|J5JyND8dS@u&2q{P&5cpdfeqojxZtH`{TzE_5|Yz3szyU4M}MAClv93y*f3=?SG+S1YWy7<bpN|AuzUsDU*yoa)nh7DdT;6yi>no4W?p6Qd>R5_J(Nwv41a-GPXI*KIf92L2|K-nrmPkDz28fVFM%C5&LPsaxJ>{3-|SW&o73}y^RT~AA)6yhn(ojaKtky3#?=;GcmYs{)D^sxv(-!UL0CQnR5IfY^Zyxwu!w1IP8o<VUcZag>ML81VIhyfBIkh_cK8}5Ot!9=XyynOSTl|K^%tY<Q01iZV~|v2`D`uLy=#GT<-=JttaU6@d$QF+p?Nyo<~shsz%E=e(Zx7Ep=E)|EV9Q@6QRIYV}XcZ0P)wyom(p@&k>8F7LU@ISLa~Y0UhAcoVu#U0qK4DLt5=IM|!5*_h=}%y;ba4JEsK-eMGUp42dHfN4XqigNFD;uZNtf@dJP8WfIEMmi;GaN2&y`<-buaI5s0^`zx@-VvY=%D+9wmM-gyeBlKySQYI(NL67vN1*<!R4a{HPDsG48X50h{#}g`#YxF4B&-GQW)yKSrWh`@`3~`1!KEDW65EU#C7TMx@BWV%HwGDFGf$jmW$YyyjiV}R|4-nYESXr#xFpi$TQTCSf5Lf&@Z1$)0zE)ia=DL?J$hJzf7okzV}$k6^wTu_h>5)RW>ZiyJ!aKX`6h`F7~D}hp#Yr-<gZ3M#3%szRF3Wu+Nb@{RyXi+{!3iz4I;OQ?xv+apETGiAd)^2nMhBYbIPK^3DhBx<MJ;$F5GIEFiBs1^Ho_mIr<n*ViRaRO1O#dN8}=s3OesI@O?1MQ?Kg6C;3DXkR0^@ZYXvA(OCaa^?wFD(6sl7|H>?Sjy>|A=SL1>8=zH(wXPTZRCeU~5WZ>v;y28&0OFfvfA(4IZqQ%YWX-(P^3I<wHCByF3x*BHf9yBM?Yan3VaTMjPDPiU3ip={8*o0Op)U7{<d}A(uYexj_6t8?uEh+HL#<j6i=?rxZh*5k!@{R*%f7C37n@H`!tzd(W=u}tcDjg10cNtt&YM0-BBVe4EO&PmHRu*t-jV);CCgj=+bpTc4RFUQAKU6FRcc+xY&~K}bCwfOndn+O;`O!7?N-_}3T7eXh?uSgtuHO^tbUW4wra}9UA<XRT^gw=HZEk*y&9CpxEgXj>7>=0GYq+4!ZS3Pnd<@}442{3q_7&DN(FG!lfu|ozJ8Kke$c$4aL19Q>C~g_>({(rQoR+>KfX!n*glexOyQ`vw|4ALC;hIhM2lk;3A=b+(TCg@>8dK82WFvN%7VD43M`w2Fh41hdAs$M&Z@3#N4o+uO3k;h&pcKsibEjtHQ8z7K+K}U?;1MT;taF%Wac|RZQ#|;l!t%$!BZghmxbrXtldRV51h)yxsB(ZVC$?xYGJ-0TKJNDAIiY(GLP@>`X0l2cfR5qQxYbq6l67ylMWd|$cVfnduc_Z(uXqD3!&U()v!M`#3lw1k;=w^J5o*ie+7k+Sgo$z!y3}Lj&twI1y8bCYf@0E#0I=d<gQl7vIvnmp?{3g&*F4nB>T(hKe6u!qD3y~#fQ<v=${an>-L$eh-p9)T!Y8yWS&>6KqW|_$%rWI*AyO8=_fC^Jzf}IzR4B<MF9*9Zdz!+ocq=E<ecGFCVJimfu7xOigbD`p#*vEpxNZLXgbdW8&D5yR9*GEr*irNUbHHqJXFD%j8r0$Ji`hiLj^Lme5oi}S_U07id+_3a2?H{S7onn)L+3|UGepKk}usSIc{7j9_|AkuOh;s*Y0e$_`a}_gNKU*&H(q{Yy%3qsA?kR&#O)ygu<WU_5%5HtY384Uw9TaBQ)mUauU<p8KhcqlqF(nCffTCg{nUqRy}>0`Mz7ONuF!SeHT7s`=fMQ^nuX=6k67{hKfTKnj57#$BfG$k6TV5aBms5xk3#2{DsM5oH@+G<i(3*KnF9?1(H`dG*Bf44gEv{2m<NwlW)I**rC69xKfJ*0JCpXla|p?fR%J5pc;1SO7kLCFRjPZriz0I$CK~|aI3Bxb^gMP)oS5wG{{&Xxi}R_ZQh0>sAH?{x9$$2hRN)H6W(dU?gz}9t?_*?R13>@s+HH`mlzlyY!zBLt|C+T#rFKfSgP&f0!T-ebN#mO4ux|w=G#TQm2Ox6Ir*jZju**b5>$}hld)s0M@$)4*anh=Th3sj3){9cpqy3Il%cxKI^Iwt(n{|4po@y!uL+HIjI*gU?TpNFcq%nbje1>3MAI7fEnV|y^dx3eQGs4}D8Ok719+JKUH}F7Au)f6KUm*%sMWMSerlmyX~!X^pfs)_AjeLB8T%v6Q5KphTOXuTj-X7!4C<?Nk9YYSpoJR)ULIj&YKxkAgKji=Hs){ePE)gImfUHGJZ*ak=<sxN))9Slw%8==r<=6l>~g6<a}l#IKNNZuCblSln6DTgqIi}?WznxpJ`<2H3^NTEI<qbY-!s36ZL&HDEqGa=oE`?-usbQc<(jtJNYF~jmC1B~CZ4ecZ{S_<4463QhYl*yc$^ibVhhdY##&#Z2&%j+@lU)YzX^5@XF<w`litD455=NFDi>AE3zT()MM-+gvFFF&N9;t{h`AXoe2_0e>7L{4b=IiIrf%Kn^t-Hux=c#nY#!O1(;Sgi24Z=&eW?jaAX>M2k&Y*H{Hj|@v<<XBOzqV>%|JRARmL=I8%j+91&RI1y(ceem5RAIP02<LM^XdHfVp?lg#_KQg>D&ZgEN|#FqCHcRs0o8`jyE!m@gF@$)NR~PDgJEU0m6YbFF37)+QsL6TOP0=qcRBfv$t=IT3YSt*(=4K8RQz{L3c+VO6Qrbb2AXGJ0cdei&$HZop&MiObJVQfTKh*;cE=|Nd%0G`ZyfHwSf6)tv1>lnO=}Q8gKrsJ~eH=iMPsEno94m_-*l7^=U<SIecjduqwyi6aHE5D1-Omx*7^`U26_lyslAr_WD@!!%p=z&U3J$QaJ!9>HwA)O3(X7nqFr58|4+e)yAbmyVa~A%plcL2!>%WIYEbH%TWp6JdW4F)Q`e(e3wL%xP9{<2Kak=%b>Gw2qBnkp($p=2*Snn3dQyN7=<Dk}1jM%xa}diG)!U%HLhG)2;r;dj4$^Y)i(`1aSP3eFcF8JJ>ieQt=-vX1a?)6<Q7z!TKo==xO@oMv>DR&da2`gWM-kr?8I0?++Kuk`tGy7rJZnj4s#aS)D#ij>2cOZ|);o3Z?KaV)-sZCMDG11fs7rV64E^c9_@Z^z}D-?5yq3!LRF#-LeAbb41-_2@&_~IpR3~rS^Eus~)`Bl3AX9*~V6wpI<B;gNM+&>E%bi@wsXx7-}&(4w^tgkdO!%X(Dd%lh*~>POFUs>S(OZTWBqK)s$QUXr0_snYBpxa!vP}GgoN#cIY>xT&xP_t8bx**u76%par2iC0p@7)T^YG><Eg}n(~R<(sW{ya=NzvcxHJL#Zu<lqcZikX(Ec16{T)JY#*6kebmPcZ8+2v(>!mCur}b*E`jxX?1;U5;nT@NejP)gc#qAMlAZJU(f$HjDI$l#@<RFi6GCG9?VskyL>4-?KsKk+r`Zi3$Y=>~Ha>CDz=ZD~j_n(ILEWQm=f8c#H8wl>S0uKL0`o`c_lQND`a4zq@d>1oVW6*id|{^l1=;Ks=@?$9D&&V0a$%(&63bOm95%e3=5A0Vd31-YEp>)!`kWWRlV5gKI?Z28K+Q+ZpiY$ecpTr1dB)%yty=US2%Gnj#y;N5SZt;<;0{Cje<n8>>&=6NC;>}qUfJkw8vJcV<m@}O@G2n|koC;zw?@a9tk1F`7f_MCIqr{OHVD>C^YoElW+U;e*+99+VT)ee!hE><_Cizn@hva;V@iNl4a8k-ed35Mnw)D*f+ys2L*B6KAo^gRXxh+1;TV<@nT5ygWb%_oaa^G5U0+-n8Z~v3TO?d<QUra80K8?%k_31L0)~lyZp6)m575Ml)c6Xv(k|!Y9)^JCh~?<eg9kf(+alSU1gnbQU4SqIfBm(jTwVol#$mkm3zf62LTV0vbf`%(6J2)vx=%*nr3H5_Cd|4Y*8r1kxGKE$byJ6Wm9TkWiGk4K7>9GbF{;Hu%z!nV@MPx(55F+@EDC0kLi!Z~0N+@HF^eiO94G)qy5lfRBi(Im{wBYZ5%UYJ?A+>r(leB>GW!q53<o1eYq_1WTRpwZ7=g822a_VS+f*X*-4;Tk?UKfuW7|cTc5f;;6w2&lw}feJ^Fie};^2WQZ8m-fz{<@X?544a8ewBDwz{0Ab98^K{hgETxmx?(iYRFDL2Y<JOu2uOulq6PPJTKG`JF#eX=J)p%2~?athGfSnN4~$S(9{VId__>5#NK^)V|yUIc#=uoTIc1'
PREFIX_RAW_SHA256 = '0ed2d2e25b98f8436d064831021e4fe10463ea4fc03bfd7101c2c698c999ff36'

PREFIX_DEPTH = 20
DEPTH = 44
PATH_DEN = 9 << (2*12+40+LQ)
CAP_NUMERATOR = 12293178397485488844
SCORE = '313987416494311478378353693/61465891987427444220000000'
GRID_STEP = 320
ROOT_UPPER = 8156
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
        A,B,C,D,E,J=(M(*f,*g,only=i,raw=True) for f,g in
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
    return Transfer(225, 225, (117, 270, 4), (126, 286, 4), (144, 356, 4), 580, 710, -14300)
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
        if byte >= 128:
            shift += 7
            need(shift <= 14, 'prefix item length')
        else:
            need(value < 65536, 'prefix item range')
            transformed.append(((value >> 1) ^ -(value & 1)) & 65535)
            value = shift = 0
    need(shift == 0, 'terminal prefix item')
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
    negative, positive = ([2*v for v in invert(k)] for k in (0, conv.msize))
    need(all(0 <= m <= Q - p and p <= Q for m, p in zip(negative, positive)),
         'packed prefix tails')
    digest = sha256(array('H', negative + positive).tobytes()).hexdigest()
    need(digest == PREFIX_RAW_SHA256, 'packed prefix digest')
    return negative, positive
def energy_evidence():
    conv = transfer()
    decoded = decode_prefix(conv)
    verify_path(conv, conv.evolve(PREFIX_DEPTH), decoded)
    size = conv.msize
    del conv, decoded
    moments = open_module('moments.py')
    gain = open_module('gain.py')
    energy = open_module('energy.py')
    need(energy['THRESHOLD'] == gain['T'] == F(-143,400), 'threshold agreement')
    rows, initial_energy = moments['compute']()
    cost = energy['energy'](rows, initial_energy, moments, gain)
    need(cost['energy_upper_1e20'] <= CAP_NUMERATOR, 'energy ceiling')
    return {'directed_energy': cost, 'prefix_sha256': PREFIX_RAW_SHA256,
            'prefix_depth': PREFIX_DEPTH, 'certified_pairs': size,
            'path_layers': PATH_LAYERS, 'required_followup': {'phase':'terminal'}}

def entrance_majority():
    moments, gain, energy = (open_module(name) for name in ('moments.py','gain.py','energy.py'))
    rat, cov, moment = (moments[key] for key in ('rat','cov','moment'))
    rows = {}
    for sites,key in [((0,i),(i,)) for i in (1,2,3)]+[((0,1,2,3),(1,1,1))]:
        rows[key] = moment([[rat(cov(i-j)) for j in sites] for i in sites])
    return energy['initial'](rows, moments, gain)

def terminal_evidence():
    conv = transfer()
    m, p = conv.advance(*decode_prefix(conv), DEPTH - PREFIX_DEPTH)
    truth_controls()
    literal_minimum_controls()
    upper = root(conv, m, p, GRID_STEP, GRID_STEP)
    need(upper == ROOT_UPPER and 2 * upper < Q, 'terminal majority')
    depth = PATH_LAYERS + DEPTH*len(conv.y_axis)
    del conv, m, p
    entropy = [semantic(0), semantic(1)]
    need([d['checks'] for d in entropy] == [51660,62521], 'entropy case counts')
    need(all(d['coefficient']==SCORE for d in entropy), 'joint score')
    return {'threshold': '-143/400', 'graph_cap': str(F(CAP_NUMERATOR,10**20)),
            'coefficient': SCORE, 'directed_initial_majority': entrance_majority(),
            'terminal_upper_num': upper, 'terminal_upper_den': Q,
            'prefix_depth': PREFIX_DEPTH, 'suffix_sweeps': DEPTH-PREFIX_DEPTH,
            'depth': depth, 'prefix_sha256': PREFIX_RAW_SHA256, 'semantic': entropy,
            'root_truth_controls': 72, 'literal_root_minima': 162, 'literal_NAE_controls': 9}
def run(request):
    if request == {'phase': 'energy'}:
        return energy_evidence()
    if request == {'phase': 'terminal'}:
        return terminal_evidence()
    raise ValueError('expected the declared terminal or energy phase')
if __name__ == '__main__':
    import json
    print(json.dumps(run(json.loads(sys.stdin.read())), sort_keys=True))
