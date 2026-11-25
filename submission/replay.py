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
 def apply_aux(self,m,p,plus=False,previous=None,round_number=0,only=None):
  rows,pack_lanes,prod,width,full,norm,total_weighted,parts,cnt,weights,residue=self.A
  outm=[];outp=[];both=[]
  fneg,fpos=(m+[0],p+[Q]) if only is None else (m,p)
  limb_shift=2*12+40+LQ
  need(WSCALE==3<<12,'split denominator scale')
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
   neutral=square_half([x-m for x,m in zip(neg,fm)])
   ac=(square_half(neg)-neutral).to_bytes(width,'little')
   bc=(square_half(pos)-neutral).to_bytes(width,'little')
   if plus:
    dp=coefficients(square_half(fp).to_bytes(width,'little'),parts)
    both.extend([Q-(((total-2*y*r)>>limb_shift)//9) for y,r,total in zip(dp,norm,total_weighted)])
   ma=coefficients(ac,parts);po=coefficients(bc,parts)
   outm.extend([((2*x*r)>>limb_shift)//9 for x,r in zip(ma,norm)])
   outp.extend([Q-(((total-2*y*r)>>limb_shift)//9) for y,r,total in zip(po,norm,total_weighted)])
  need(all(0<=x<=Q-y and y<=Q for x,y in zip(outm,outp)),'auxiliary tail simplex')
  return (outm,outp,both) if plus else (outm,outp)

 def reuse(self,idx,previous,round_number):
  if previous is None:
   return False
  y=abs(self.y_axis[idx])
  interval=1 if 100*y<=self.yd*self.core else 2 if 100*y<=self.yd*self.middle else 4
  return (round_number+y)%interval != 0

 def apply_main(self,m,p,n,u,previous=None,round_number=0,only=None):
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
   outm.extend(x//den for x in lower)
   outp.extend(Q-x//den for x in nonpos)
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
    cap = F(12307053318118561274,10**20)
    coefficient = 1+1/(2*cap)+saving/3
    return {'checks': checks, 'local_reserve': str(scores[-1]), 'saving': str(saving),
            'semantic_increment': str(saving/3), 'coefficient': str(coefficient)}

PACKED_PREFIX_FIRST = '{Wp48S^xk9=GL@E0stWa8~^|S5YJf5<Hil9&|LsjKYSoQt*xNJM{ywN^KXiJm<08udp`4pw`1R`jYykr2LHnZQlCys#4*MbJ*uIx&ogyHLqsO&k(;xR?swSIG*JYk|E!w>ybn1vu4brO6U=ArDu(FFZmJ1BlDbuZOa><N0*d0n)e~bW4ULeGr<GS7f4$ns`GX$IvquyR5aY5*q&A;7m_T!NEq>*VZ}>9gV@OY$hjE{^<o|T0oc)^nuwDL~D3vw8D((t0)ioYl^XdlX&oEBB-Zn(c4%Ph5kMsq6qw6>4HJhw9A@22C;@nnGZCC4lgBT)1*DAWtm%XzpZI~sl!`l7NO`NYZNCFD%csOzI-vPF!m}Xga%9%VEmw1lxZssPqNzwFQ8o6C`(Ivh$uur^lDkmn9UdR(G<Gs&E{CfaR%Yw1m(=cRmt3l09>D#AI<R10kFVWw&2m90D<vzg5)&0QAMGZcA51QYDO+`MTQ~dF>3`lhC@RQkh`8igfVA%T!ynj9^&T~&kXuIp`J_u;{j?Z^GbP36ZG<3*xm;wRJj@ph98jqh6B1jR~*XVJuNGj1xpL9OU3SH{(aNi~j_d0BIY1yP)u;dk&Ow^BO%6(A@LBZ0h^$<jv-5RUJ^|{YvwF$QSI`FKLNBf5O8r_lhtX8u0k+_CSXX0oJgLyT>^Vmp|&zTK%0H(%mYPAub(x0R?cvHfk0SV?p?WVg2k-=r}p@@#IXyFH0qYA-^F<9`%`%Z$YT<7b&_qy!nrEtzk=xy4yGMd;mm|8@A9`YV`N_?!%*5!;^4nKcLkHdW%4pPCRxLxJ=9miqV!Cy``IQIKbzJsabUiBvwLOu3v24kqmb9aP0;^v4-@$vVNIR?UHA7VyENTXXd-|yl=JN(MW*CWq_HZ6)>rY2`1$WpeG`3P_LS6oYTSfw~tIPI;seRJ4$-VrZCXD=HYA>T3g^wVn!XLl|0_p{N_5WS<&0MN>P%!$yaP_UtxW`K5%s>o`eh$=RrNf;_t&AV@p17U0G>Vq*hX8H@hJ&XSs7u^P8bucuz5H_en!&=7^8Y@PUd*$dHhCAzqeF8g%^EaIp+g*GpRZU#zCemk0=<a(p<YhpRmhfTGRzqTeLJ?V)975nsrn)TBieXUYCA<y-koVmMfGEEj`1i4jn9V23AKx;<<-9>Zl6S%-3>Wn}AEKIq_AoR{0>1YlQD_N>;TrrD<l+ICBr&z@Ak)z}34*<HNrCv@ZM_5t8Jx96tKco?C%O`Ny({mrdC28_hoF{Jc|V`gVa)OqVW>a!<vV(%r1{u)%}2rsm3|DTPepUV<gQ=SW!58ZC9jPKo8!T((VwdP#%PtoZ=vG1ZKoiDZ7G(LCJWTPh{5!6QDMRvm#3(Y(oJ0LY}$|}w`CvXF5qtSJkdQbiq5;#%Dawp6$7Y#QP*rtMUxk3ms=&IjXC%4d+or2pl;8jbPJ4irf2p&6gx+8-D;|Bf^W><8Z$+3fI?i=P9{>dv!iN-`fqG!;jD^!Q%t&iv|~daeh{8?0*SUS@>$Nl=2x;ix@aVYf*a>CbFjQH%i#=+CFt@SewlAv_i6~hof0S><^4~7StJIU40o7UicU7+ZB3ji`)$h!4MID~LL~Yig=)9+JV|;Yo8;<<@DtiOA=@Jzn8J;Ku|B7R0<~x{fe>3pS%Q{h$ACpeqU*b(^eFKf<|q9zP;MuU!XCl?Ed^`NH=>qAD}EjQADSOxdxITge4;Y~pVxI??neF6r)uGNS<yrT3eFyBTAiSe24`~-pQx)49?ofGQFdWXjw{c01q3kw>M||$N1Sc-Ko>=cw;cxC*1w*xa{`F62CRRtTo)lDtk4+B{Z+@%Lnm?^2B&Z6xD<_c{zoNEQAKH3I$$n)Fz8S9!bcO`5|r!F09Hdd_tyQD-pkfqs65mMj%e^!0tc#MlX@%}f9q}H(<2B#G$6d&y?!b@J&giQ*m1$TNk<}dva<yF4}o_I6ceWfESeIClhez<(;=gK_`ws3RqTctCP5a!6|My}pkja=DdF`)SUHPy!v$7SH?-SP$2b#&bF^o|J)ouiy2UT;%2>rcz0qJ{`o7mY9REocS`1&yR@V%kNk-^&(Phc%R%O|XDxS(=?g@>#;&?gNW`kXHJ_*WYwIF4`Pzgb@EQHNV8Pukis+a3qPm}PZU)M31Mg8o%5lcF^ZyReNNACMK!GW^Pm3tl6Nd5>V`V_urnxzkye-M%+(aB;=Od7T>-zOF~ls_Z;cb<M}WNhDg)wT|K-3CPC?!iSvRSCnBgVrmwlhlx%eN9gaVSgsSiqE${En<2h(mFP%5B6yJJTVW__sHyQ3SEa;SmUEEM_y{p_uZfzZUxaA5fv=^ign|e%)kAggtXiA;u35jjM9e;@8c4CwoM2(!+2LpZ1%p8eR_pH%w}LZ&^*}nZT?9&)O#D?3o9xoL*?TxU4+ddGD5hB_6ej#BZ17ZE)t<jxKN696&fx#dhlZhBmm9167{mKWM&I$Bs(o4LBeIMG=az>Jw;%Ov255ypsRuSK+eu5b0<O{B0-PQSMjMZR=9cAFtgkDYmn>cBA9ZXjB4CHQ&nNAMc_OX7~>AANB!az9@JR6S0&7r-0iaO+_O=@a^L$1xu}tW(yvgYlg7F_KUo3@U1LYI91~KE<Ov+|>xp%c@wDSmsMIrIe&ZmxETi?8Tu@o>hyy01gp_cKi5+3#2HNFPMJo$FuSbXVwcHQH5Ye63&;$}ZRGBR~e=8v=f!;ihJ0xL}_FBRp0<O=P>z+Gs_DBpwzBOMCYAdrrS60k&uRfx&rT2jw)Eyf2EFW@Ov8e0EniA79-Ydsdwq%u<NG%i^rDNp)vcTJdnGOe)Lev_-RX^*XrnvJjb3I!=6FTwc6wz<IM9qhk6{<F1Ui~sX0z}ub2kB;$vpSN+wh$%C%&17q=j4t<3>>HwZA2sjc(9dtaGr-ohQQuk?Z;nF%9rCIU2UwcX!R%!@dwLW=Fr(Wc%;3D85~v7-=lEWRF9Jo_$}zH2+l*uJiwuiyggyp0`)}s_}=F>e7f}aEiMO{wHml{O$g&Drx9$Gz<$t)<aZd%R|-BsSok^Il@JDND=4}PTQ8Z#c%3ZglLy0jwf2Cc1$Rf58n1(>mv8Cl^(=kayodJhVJ05ZR~yE4lqd|)4&Ksm<JqHX$`<&XiVjWAL-H-Yzwv2xEpKd+p2bs}e3IC9U5Q9Wwm!ceKBDRU-Gq8q_4*}mgtz3b+*h}G%_b1|2EMIlt<Y>U1uMELKR&cg)nbjaP#ZPBseVVxz=j&5YHWBdXv~@GT_6&+{-JWAcVA5DS&$?PKd_W!aCKVo4FZKSm^Eq{X(O9I5UjNRiJGGc7qV5h!Gx%6d;puTWJ-Go%lD))EzQpY;3xq1ne;hZx(nb!1+I`>G)xcCAk{rV)}tB52Vl?s?krO~oS`Up_%=iP!g&4Hj@LOpF{2nMp(^2}l}*6FmiC!4qUbu_zyvbzJ1W*$u&5|a?55w1CymOLC(2&fc!L=f`40a$8Sa!vh5ZLO3P4oEH3W=v8r<o*%7<=k)%ARD0`E?lc=K781xht!tZy-v`^@(j_s%{(Wnq0@04G4k-Am(00nyOl_hT*KV%t!~NGV><Pok#msJe7QuSRax{FpK1!Os=ZSqEHS4Cp(vJ$;Un$SG2K6FOM_)NOT^hD0=7Eb#0VYQRb~R9IN&_YZ{=2QBB$eJ}hh3FbVo*Q^YH7XIo7>Hh?{r-FXwtNT!~>SoDt2$U%^z8lTa{oiG4{_y=X0Pid2UXIIqSXCs;fdMlOL!==qq2zjvWgDpsS5x?Q;O&OQ?Su^bnQ9<RKCtog)tY9!vBr1srl<6|6EDuQ#n?I)-=jN?vr4Q`tJl<6ds}3VKFwtSi_OM%PVoN5hG{B(kwiKaHW=om1m*HgP-$VBMx%o-WHMbKNqljp|EiRGTKW9?U9CP7qYHC&2utPb=d?@hIh$)F#|{MB>{4h>Y7<#$3xtWd<YZkLiBHOsT=&d97fV=7=AF#1;q{YifaeX0R5$$J3^Oy;BTR70wleVWD@yNu(5p<u1@k{HINbFvoY=12jNquZ*IYa!r11D_*zMwKku`5H7FxUhO9Bxv_Gpw4$Mxk+XXZaM8Z_L|;n+WzYH%l^v9k))jltxn%b`L!54`L4gPLeQ9}7!hOaY-EK~)Na=<w<gX_nc1x|bAkpCQ!Mw<r&zz?J_>G}ZO8{%JTY7YziqZ2Gh^)Zm28RZ{fAYV$MTSkUnq2L;dESGrnoY$3uIOT8~~O@craN~)<Ic<Lr{5D91M=3N)vuZ!B{3>nt)p*Ydaq^h*J714q;-jxTglSR$NK6>dd2MeyI=>^P$hSu^+dnO`rUF9SvbGv9~z=FiWD@4by6lbv+6m(s%rE6hSWW0DS`Q<jadDs+zgLFT|fa-2@_vPB7xm)`@lz>*QKrIJgC;V{HP8V!AM4K%><7VF)l=AXDH<0m^<5P;KZAs{Yf6D#4;8yd!KMxkVkO}K-Ug6H5Cs9N8Ch;O0LlrY(xI|qFfz7+u^~mbmoQJ?j_m%8~Ww2(~4)IYGqjFUgs0s`_v+FA>X1~o}!WBffylCxKyOE(`9vj@^G9+<*cdq6a+|C_8Vg0ht@WV0S@U%RZQjp6QI$E}ic_jYu5oy;i<>_)A$~MkG2VPPX=l{&TzffH}Wl9_>;!la6m+4+=ZN9>&1>G7w*$uhZr7zM+wqqLx#c;YtN}6bmFS0_p<-5|{TU4HxxR5o8QT~=xC@Nafdc*>{dNeL7$IQjJh#H>59lwYokFSXV2_spoULZuSxCU2qsKpP#-ol5{Aq?*`!*|<4SM7R#7&6nFhl2>w28R0Yy7tSPBhJe*#8Y;??hM?A!~P%dvOp$}*pI88)y94LC~;|XQ0J^Izb^Vs-dyB5Q-^PrCzkOVxI2ux6uH0Y+U+~`!^AniCi{qBThDvd$ODkFS_6L^MQs7G)W{JQAa?i-9Nh4l9QYcA*C{*r*&rP4gJPz2aPx$QN`I%fS;8{B1zxY>1NOk9N3|;J>}ZftIuyL%$90x522jiwqC#vwJ+4h6IpLJwa?n_gUvRPzi@6}$^Xs7*SfQ#sCGE%oCS2Kn>b!o}HCB>=;b{x;ydui8$kqrqlb;7$b{{3m$#L!$y_qQIxJK>fJ?b2s7MYSp&G;?Chf&oBf&>iq3K^+kI92=d9r(-TARq6g&l1w_`)JVP@p#~^?;h@z0=W3sEq@4uxKr9r;xor8!RHN#L^cooI><475<#eIQI^I<y~Cr*L^#o?x=Vd(aF0%XleZlmboX~8I`4I1mna&Jf@ZOh@SEZ+IQ?)q=?$&xOKdpdgpEV!4Z9StK%><zNpo06&X>o$13%`B4|2!sn@@BDmGCMyY#dztN^r!RlCfBh7?G{Sth=1Yq}Z(1BCJbSkz{;gMo^kYgOqbJcNr9NY4bH3+p@(jZn4kPs)%!1d5zXb0h+kcrma23vK4nf9_4<UOBUsG5PEPS`GqefVqS6hC0om|<k#7_>N#oWqVKJ?ETpgiFsbOE{PYCS6Lh|i1S5On`*GO9-8nIUmP-gzOB$<~{jO@~$nJ~u$GccHbge^eKsaL=a<|Ioj(zV;IXYwji)aCf&cpg-)pZf|r=OSvF?aw%8ejx9z@d}u=nBU8oi{|qFq6hYpa~JOO^|VAv{Z<r$vk{MIUQmN>cvSx@0kbcv+V$6D6&V*HX*^T^`xG~3U}F3i&c)jrE2SeIi$z-uBq=;S>X*3aF<+MmgBioZk*Yn#uc@^Eh!`_2t(RQV;}G84%^0G{#xYiBXyCLzsIRp1@c+oi?i*$@r`Z9?gM=0@+<n@u57-rZW?J;4Th<Uo9;e-jv&@dT|mgS0^pg=zOP*oZ@0sDC<NSyu;AqiXPW_A*GKS>aZ-~rWqfx{*`WudgJvnnRnjNxfiOJ{Om#`plbPqcH+##Emw~5$(nb82f9!8vLurrtLzC-CvQFy#*)2mkz`b2vGVz-{=}b6O50Xl%>F)hy*)||(R*hzp?YqSInUpqa<q*K)UEKTa_?t_^f}_>WmuF<P-j!?i=w<!fouwBE@KaC^*-C#P9QobthG1E>CM<501Unj21GgVX@s&WSmx2TZT$w9Rggqn_)mP*+=TtfLt^x1fO5YbHaxFtU9-)$x6USyQSeR<z^(qmB1{z_tZRiRdns>4*qnHg(HMGmKb<EwH7Y)_8Wm~irHd4f2^HYzaAy{jl$uH_`EeB_-sqR30hOH%PEBRC=RQGbr^%58o5TN0A4Fs34bqI~eK5%(xqy$Fg4y+C8L31`(!}`HWqZW^04Bt4lnts<`E<Rf$XLU>_7)_e&Dd<D$mZm`Ty0A2rM2E!o{6D||^vY~?gl!zH!|Qg&hxaeU{E4a}C)mEGBL`}5oX;bLCfi+?qq|u5%8mI<7d1*;WNvi6@wU+=`UeX02k)Je{TuTbe2@-GVpkRqR9>$Ilm7nFpohO47^B}SYPb{v#S4;la~sdjJY#^o3xT@<O{hKzGlG0_>GXy^hYA!JRIftYL+Wunqi$MPc^@5_%|FE@wD8F4_c*kZXW~)>)5S?kR4Bxf!gc5|(<>YO1wY&uGsZo8!q_^{`U&#mM1C%cJ0C;=fayCfv3exk2mE%Wpy=k%XQOYD@y!{#zY`E>9E9OP?Y*v}$5co&OA-_OJd`>xCINh2(X2v_*|6A=pwY*wvg1@<q$zSx6Th4WNizt~bX`w>b_20y$U%!OR2Kx~6GG6@?wy&2dYhVjIAF)Z_i->Q*~46#G&tXq^`i$NAC}sJ+2Kmu;|5Ufg>C>)#*Zk*5o8A*R7QcWIC7bNJUYS4w({A&EjTrJ^z4)eWp4B}-5owx<)?z;nu7X5x$bt^@9~q{{I|@`Sp-~IB3Kvt8tL>C%OFK^NJt%^d<(JT5!@6c5lcaar*kcA02thMlcM!)B8dlbV3?^>ubcyi0(aXNRg<yR^qK^`&S9!9b;MHq@I;;oLZW?tdesq*Y`EhEUyw|?6r!Q)9w8EuDFl!+I>G@TI@$Zs_m0-l78~yr=#l=>*Ws&tyJF(M06ge&Uz`iUeu}A0UB9DkQXsyr@nHKmW(xr$i8KJcyVS5^$Xxv24?Yvgkr(V!oh-{XQtTd}cHc0jnaK=Us`D10Yzn@h7X7ML6v~xP>~gh8Lyc5`m=Sx+geLeLY)&GzAC(JSBqQahM6kY4A>s?|5jH3%;&9|y>v2|h;$jqV>i5ZOEGXY^grh=wP{W%FN0BkthHRf>%fGnnjAf#uBBEx!(nv;waQoA6Q&#Ybx^>&#9W`XJz!N+Xx87SKlev3g9K*Okh~pFfyZ0n1Hg|^FIS@C%$AmBj7rMt6dYY${EQ`uX(8icn*38XO#i6C4F^90Fb=j5plg3iRylNzph|k%vnkm2Ln~=g<rEj@_7d)d(ra;$N&HSygKhaao#qakR%U9<!<)N+sFemrIPlp+m1^e^JQU_i#j7q7@hIVGVql)?G9iZ-0h6Md=%=o~yzwT-n$Kup$;({ksf-YAvlEIb{vE4=6>?$|<JSB;+1Dlz!R1A1)C386f28*D(FKgLZ*ILc(1SjDET9g;wHrBjX1`s%qx>1ne|3{AI;d^Mt>nXnzoLmu5m8k{(>)c_L=8XmBd)(rO;0!Q_S=AhjeNSt(U-Tad9t{t%aYGai5J;Bri7xq7iD85_nHN1eGjKHC|6XQC-8J$Gj$e_R?A$-&0@$|GW1|BQf|;BX&P^{;t;bJ!ML1hf`6cQ@_IPg1rp$+bm&YUcilSc^(KF#Wb9h#Sp6WaLa*UWNOURV7%J<baBMMTpIN#R(@zbe!`1^4$l@2cMyADg67fveiSUdcV?rvR340B+ruM0Z=nYsZu)2rA~*f5hC%4cBiO}up`5c$EW9TPte9}R?mt}yvsQ8{jO+!@5R)$=h!d1lgC#~MNPuP$+e;N)%+?K)DX;lXqAn3W73k^urhAe(PPxWXo9{DKPuB53IaYDhQeeN^*;_(j!<G^Uac0$7_y$9=0pzrsc{J;zp|Iyw1)-KU-fcYALd2KKkL9HFy*RaVtgWfaP2T@+@sdt^?dD3g5<13YbHV3?|Z+}4rKVccI}SBj811X}Pa$7RArK+0E~sR-uUY;IAZ@F`jFb(!M3IqHtYAQ_9w)#$DaXP_T7E+F+<DW7rRW7O7Lnpu)b%oGhHqdZubG~anv;B`cOHG`m?2{-X0JCq%wDNTLtz#PTGK{?I&!?X-_)vzyk&E^L?0*zNw+fQA<1@+|NOKo{XF^w+`yFQ({cF%4h^X#-8s`6D{8ec;V!eYPY&X<(3-369-c}2{DjTn$9!M-tkX_vsVtYf?{a;`3e^2)~2=Ym5BCnP~JdUtJUwQQ5!ol+;T58#9cS79UEzPNJH6pXjsIfm5^gqxum+eH+QdE^{UmX0!C`lq1ltasSINK&|7g~5gd=q=5PI2W=V^hQ*2al&|CfjR=i9XHC|L~|6*7zsFDxnnOH&;;2yK1nq-&weR}87Ns|;Xw#O4<kpFk%s`tAdWbmG=gVQmqh<hUZ+L4-G0Zbu_S2>;fzc(O<woWJnI%P3_h3>bKg)tP{a^;8x5F+az7m;#;p^ADfoG;J&WUjepP+49>Po$&tRfQS~0C7>uJ;sURRF0*-NqSX@LNE8#WcK5qZ}R?jTZFol7z|y4X=BH-jMiuR2K?c{_3%^V1-bq#H*iK8h;^sHNaL%zvpf<i(jQmnY=%(7UzkdoTRCrHvDYgI7<)`4Q`kIj<vQ@EYY`EJOl1wtwH|WOSl`E{hVQA2yu?0%%SlGpUo8R*YbX4Ke8mv<0Fb+&)Rx#cQ>WO>qbA9XIJ5F@(@1E3$)8Q1m^-qI7j>v&E`z7UnlWp@4Qc;?XrMW4R3QxAY$4X9xCbDJ*yUYxK3Yr1(wyH(sSj5g*1wdK~;u^a=-7Cr8FH8op>g;%cVxd&oHtTl1Xo_eDs3EJ{DQR#8I!G$r0`b*nJPdjN!U?7wM}P;xzXp$_WWe-`(2)+LG>m22(sB?-iHb|fx;5s<l1r1!E2==Y)&{(eSJ*kN5p-G%#GFm3aV3}yA55&?IOy`XAw%MdTcI<q=^bTf?n)g~PU?o57|Ak17_vye}#Np14W9C-To2cSrM@^LPffp1!TkK-+s95X=#_EeM4T5M^uaJC4mEvrcAhr1V#309&Ep)oz%!(o4|%KV^*fV8#8@XCs1AMQ1>e1~LGsFf<otAkJdZLbqR?W)iG@%!$#1ybQhsoE?Pvm)gtMkgl1DMYbJ3<VOLqBCkin{ouybeG1Fzas<BXK*|HW;1l7*hkh1M&lPJ+Zx<j*IjJX$hQK<e4AYySDr;^D`ot5JZ>~vqG^xmjZ$gvCNyKi0baWFA#?F!KQ34BKk(E?uy-<PsR;9!lq(e-pzRN2?^;A&KBAy0bF<(L#O;8Ri=U3r&M9Qi1{9K(9GjC3g2V~R&Z8uSkCRjDamz}@6}@%2Oegh_r+`o-gq#f{HM4kwHHk{}ivG?*HY3F<7m!Oq_bzb$rJR3ZRVU(@5E_{lQJP<vRrz+weVGLHjDsgML$*P`L;MukOI$PI!QSlx?!78pzV#%o`iw4-R(6t|SK^6RL=wJ!wN{v|az&I@vO&b4jrGgv{sha0YAJNP=)w(|K!Aa?ps500GK7}d&0Y;0^<lm{`UYd`OyMADsc-j(i3r0RF>PP)`5c?!W^Br_o@p$xb6Yc1atcZJ3?aNkNHLldf30D;i3y0*>UJiVdxBQJ{mN#rhjbLGn6WScbxB;By6x<(Ss+RAXf4faRSRDM;_W4=#!}nS_(&mvIllDIGXIfdu5&EV82eDEV!-V{f@0s~@L}jdj4~WBeZ&t8NN41Mi*k{SA1gzVswN4G-?!%SG=QRm)=Gj==3>|K%aI;a;kL6?<M(>iC6JfFM9vJ)<a}(!)8a@VG^h*8gM>@zI_9$W&B_oai6#g*VeZ+q4*lGK)Vm+*j@E~%7yeL=+$i+Ne4{1A%*Bba^pPeDn^L+K<noQNnlA|UwTMx5>#K(ode^pS_og$}^63lD*K{?kASX`kwXmqGA4S1iW&^U^`+BazksvjG6jifPK_jGrfzT=<AP#0;3$3QK&dxzp3B@lNlp9JQZRH1#=h`3-9OU5@L}v5NpAS09d+C3Ht!bGo?h4dc&4p4lJd-7_fcRxOFsGG79|8>h?QbA`?Y1$IOy$iIchj!425q>0v7;{;aIK&4Fyw@AAO|P=naSji0CNg(GnO+m{}2?|=Axe)fzmJsJf!5j+xDA3xCBbGzXi2boHi9xB{l`-X07D%tkZuZ7&Xv3YdR<Q6Or&QXMY%ZCEiPii8}i;1TeT9v;2H+p!;?cRbl5q{`nX1ldA#JfPBoxk`CFOpaPYJ5eH>FS7IUL<sYg0ED_>qhZqnaQZ?WQ^a~Omh7?u9u1RQ2V-iaJ3k&f09cJ_m=6A3LZ#1mc%AeFv6#2>Z#L9rgTm_EX%<H&|OOJ;Uy$VsfAy_-nGnZ!kvIA<d8F^a#psE&yy|c9P-Af1RnGgRou>ucA#I_%~S$ov~jO=tu-Fegmz2&+zoae&Zj=Z;rGMtixg48q6F17xQ>1f=sgEltLR(Zu-dRvI_eBJP!jeMv_!xFVA^IZ5POJF+o@Zswtd9}7@1amhmJck=;K)Xk*cPYKL&Xr0K$Na`W<BYwW`Wjiauojk85STN}s3$Iq&JNp&(;5j$ILysXt0?k+0bw5@SoBz!s??Kj5?0@tn+h!aecO^6X$T&8);T<P+4#;ng3`_B{K%?ln9113vT-{QlH>YyL|>E#@PRKKDp7w2OFWLfN_ah>&&MAKtSfNMsu7pU<elK(t2{(mqjkB?^yrAFk{n~vV4ZFVP=@h~d9nwQK(Op>Rh?B!HJL)NJk!^Nl{=MEY_Z*}2@-_k8OSC1^-R}oN94=^l3&r!2y{CCK^vH~q{53TG*7<HRX@^VtQE(;76=}fm*IIH!5<+K(HkY)3Y}W@T1Vh*XA0sVn%#;9svODm4bB>G0YFMJBlMHW@Tus1SO7U*Za81<iV0d9Z!3;6lrxS%!Wb}>dS82%I1Gk(NHYc>O^g=GlgLerh%G9rRRIBL<7I;)JyL<dQ`_Jbk-g!%!vvpd0@UuA=zq`ueAiGRK{%3uPVsg^=)IqlUBgfeNKZ1}W1&LlGQs`+0A8gNyC;NZgH!Xk5}z)6Q|2f;Kg_$E8$L=8jDSukWmb{ih`BGsH1B1d?qGT<E6!k%Ef>-Yny8miVl0MW@jfZ|gLj53AHGHcHMG6G#FQ@3s;9jZb>WcH-{A6M+gAHsMtcX?OH?v}*p`XH1}oA<F*%4>7qLY7GNjA(ml8sCn^7AtB3MeaT6XHJt9R-xX+3?aHU(x}fb$ZgF}~rc67z@hWjI@6%K3zH!(1Gz;XM$4G5`iYBg{*)_oJCS*E&t4I5QelOe!j4)()=oTPRK6cnpJTysgDTZlrQj#t<w3f7EW<w1kTp4%AtR)3wG!e&zTiztlbuijg*MT|RZB<ShWOp5eO)-Ig{@1Sl|~s-z!};Gaw8ngBX!G_Z|O7Zn*Lvz&85Jt<zZ8u}}l=Vmh!EcaWTmGgajY7{MjvBXElc32@>l8H`0eL%R8)>t@(uRdXmM9obb`mTOJWJgz`Ba1^<N8b#Y>L6Xg^4fp4SrDP~PNNSllyOwg!a3@*YA!v?Q~6w2sHa<N?rJz@2#E}>cJ0m}8>2WEmLTYVEIk_0^h2x0ild?w@bC5w646DSA6eSLkdvu1>08QZsIE1<<%Dw9X!2(#E4x;=bLSkbpdLk=vmIHLxFRB`hk+CybBSsMdV)-qR+vT!Tx<6_nH~DyS_IRgS0|wix{vSrwOT4SBuT5QSH_(8Yg_ZS*9mT{7r!4;MdL^UNs?lazSe{kLP|;yHv2x1t~|^we{LU2kG<TAuv2Y1|J=G>!1)FsYX{x6px!v#DG#~gt`RXATK|KY&FL_1=Wog2W17cH$&_Q!ahG8|F7wdyC4&(Hwk>*9kZW=XNXBm3aIJNxl%5meyd6?PXW8{dI54DgTD{j4!iiB5%rs|qucdp8;60@Y+PU0G$+j*Di?b>m)@sBfEa0Vf$F|{2!L5DLO3P_L!occPPKRny-6b#sx;8k{SEQ<j=l=M_lwn-)A58<o<5hk$bU6jmASu;=gFUlowVXr<+<)(1moGO((<>n4oPv_L1ma|rJ*$5aKvz_Mvf4U%&<^AJBEbx!Jjiln5_PvVSDA}^k{mT)@^w2N?v`0-o=5kCTzjLiLy;I$aPyOLK~8X~_S(uXk(_jgA2#Dp6JWdOcdw&wICEuAoD;q>N`MaLVVd*_jc`^2N<pyX&+{q^`3t+T1udNT$_@#W7p+A2dG!+`bAn&~fwUuixfI5ADx5<g5~x#Y=%p9bp{bkSDgZs>fE>1Vya^X4DThiL3r;66To7!$JRD*X|EsywW0letG#;;Ah{95^HJ_GTaROXVWqn&uLKj3SS$oAX$+l+pKeZd(K=V1G5Rzx1{c%)R=#zM)G21AOP%}nSC)q$hyTXNBN^s@YgNmeI9{};?Oo6+PlVHVW{|J5>{>c^@{(gi}u}XMq%|}v*SmBU9@F-{&ibk=P)6ss|)gMsXS4ZcAKhI~|RRaJ@kE+#PAG{PB+m<^r1k5#JDT;=x>Iy4p&k7#louX}TUAsL_D}HR<dK%CI>_eH8Fu8eu?*@3ymjxf6AVIT@sLz$tW4FhnK4!(}FOogE;^a;Q{i|&CiHCbxnVL3YQy3xUZp%K$AatqKNCCZj=?J(@NJdnP_~RFijL@-Y`{$#7;2(E#*#nop%D>#xfq}8-=ujfEOfsIhp<`+3=H#Kub<+Ljw6p@tq{n2<bZ2Cmqa7!4MaeyQSght61m~OyH$>|ul``<VmG)gyivVOwv(6Qw`vp~E_+fEAp#N3Y>}0@HL_O@1uZd9KO`!WZY<XWI7k<w*?%E$?9nspmFis1E{b1F5g3r_?-JgW)K?YXTTR&B582e!D!SGC~R;GpFOtz@6k_oT6Ra7^YsLx*b(zK~oD8l!EO)t;Ftv?I+$lqh2>vNKl5=~|;!HB{fjF3W{Ad--b&={afB6IiPnBxcJu1_$)5f(Ml$>RQ(80CT}Z#G10qfe6(&V@t+qL?B@%$i~BkGr?L+aexiw{%=px+1~FE_LrA6@vInB5)Ks&(dw?EgYv|%N<`Zr`rc(r@p-4c33e32x$u#L|x(^tVd62GhYYcqwVibPCpx7f%UqBYJ45UX(+vkcKYN#-c~oCWV|`v<KJt4p&n5rMvR6)@I(KMsiVSp!!9TadL5wjTc}d^T3)CGD>Dvl=s|-oJ!A3*PmQ@<#b%%zR$tVh&pdp#ZpFBu=8ivm?bkLnRjqp0hT6-TpG-*0IE54}35NGI7_*CM5@q36HtI_fQ3nlV^y|?3mrs4jBja|DGG&tn-{a=)(f5Ip+bv8aufKv<){|haLC{Jo9wD_EQ!{rXLidZUMwh^b0?7{R)|0`uK^T%L7n)Fzq81F1{gO-VoNfeqSFT4iIRO^+F_BcmJ{T7AKv7a^_vJF(pNjxk$eVGfk}HGRwngapMJFgTJ9`xfVx1akCgv%KOn<jEK&e2scq1qlvMeX(z}&VD>cS4XzYS?)&W-2)RCyGcGNCZV=KD$e&C)h6hEThES0D=S<jvNI*HR^vDt#n#Yq~W1pj>{<z}ZQ$&8v9yiU>mRT~+p`iV;rUGU819gM`e8uP9Qh;??3Jvj`<j^(?mRUU#kstm9f{2C}r61)DN+N+_^Bwuy=Ys~xz*eA^u1+i!pO>DXxi`SBmi6u~9G(d8q;M0VcUaIu_|q9sU*X%m%LoZqvxOqWiIoJDV+YC15A7(UEfe$?1cp(8^3vpEw6hX3uHmRi8!O!qbx@K1y~<XZ=2wjh<g$7H{dF>B7_toQX5U<XEq<@ijb6C1rvw<HYQ_sa@rlwu8ht3UppKX%F-LTc_*#B<g2Y)q!*+`%WkKRuM6T<nuesDrX6_TR^@U?{!d4NEDiHCqtqtNeWcw7T5F8(a1{{@WbDO8*5vK3a_=j(Iuj-%ahlLP{#-_Vp8R_pGJ6@r5z<i)yt*kEPjh7=<zfMv^jH@8Is7;TICZxWO&E08;leWz)Fx7@WG35dd{G%n=RxsJZ#FY$o<gPNc`Yp}BIb2(AX#WTHr#EH>@=TR<Eif>x!z5CivOVrA%@9o;J;rL@)XJ`$%NQ4MA?ue@uC@eDt7w05sXrZ)Xl12s_{P`~omf$a)GoYS7kUarHMdKsfP{2eu4L$#KvF_bMa06-RHtmUY{QJ-^P`7K5w9a!@tr{^7JqwvqLW06G*qQ}r_)LN)YjN(V|Cflj!$y~ZWLckbZ+JyeVcVE7BF^MHwMyY~JLUDOV!7np9rCfd58&hC?4B`HMtH3rag}YLPWGu{4#sgUH@2%@$cxxcScaNcG0?^4A6}ClDoK6;a^J3RK0DK0}-rwa21cOu7c(nQ)CYRpMdWV1+41{G{aN;pdeVn!Ho`4YOMvc!t1qgArkjqA%Q1v)+2pY@Ng~BfoRxVvqBFn|@(o(RC2Ol@&qiso!8+peRnW=%G<&fh6j_-V?TQoet%mIZy2<~*6n;6vA=eF^7PYz*#C1O`jQyBNr$MYHPZ$UVr*?zM6EBHW~6yNbzr)oO{#qg6sCMyx1J$Oq~u@hXcy>B<U^#7NiT{WNL5H%1YSe=#wWpP+>gFMB2PYpmJLwYd)<uC&0z}jgoOM9b?E1_`FeiySE3{aL?Fmd|4DG5|)Z^qu(?8I1J@uXTlgK-T0h0!|(Tvs$CX|ueEP(3wk9;lwTwtb5BX=m)Z3?u<U{Rn4&#Kujh1o=VqFoDxMAGC$`7C;EsBeR{ny4O>e#VN5i+=zHf8*6=qiOfi@2rFu1Bb70GKTo7rh23|23M(i{03N48EUi7PK@%C}zX=ZRMz%*~zno*rJcyvLV=ruS41<5I`$e;(G4LF;IMZ-vE*&<=8Qco1YD=AfjgyxIYRH<4XL54nEbs^0Ltmb~-SMPiVT1RzGV7R7#-FV4l%;$H+T8Q+$5BntiD}TW!9UIWQoqmfKlP!lA;Q^eA4(Gw=E@t6US8`5mYIb4Kh!3KU3-KOxZ^VeYYs&RtIO%{7M@JaiA`6aOy$+BJUgEe&6(p5YK6b(%bu%PIHM8TU?<Y6W?=@_zOs%DYvQC)o_RabFKH7SLLeCl05Y*7_C98aE5dl*!L#AcWW;7I+Jh5sMP?<|2vCBB;l08NDLV(dtc7o&JeX5q7w1u=@nF>nR`FKYLXjyxo)BD_DYNY@QdA=`^2v)M4D49Q)(6a`@7CG}LrT&yk3jhbs@z&}%lu}mImM3dK`CK?`Y+#mqP9IU)Hf-1`rQSqHcJ$jz0+bsMqSNHvs6*M_?(+LpFm#;=;I&Y*x=lh)jyPDnS~%LZK<SQhDNlgp#SR5SKYl15|G^O)8+b^_y2<qSpIIk*-FE~5eQ;ib~)kq%D%MR$Szi(k|>e-)T?@gG1>*&Y@=x$yja5Ymn98dyp3UV*6!MmT=EY0B8l2Ag0F%oWV@LA!YwF)<`(X2B$mp*+9h)ibESi>EegT`C^O?=#GT0`E5!`To3uFbS}X0EOG92(>Xj>|<<Ci0n>PFY7<XmNHh0ottyln&>BF2gAfkJc|E~30?S%%Lu@h_#WvQFUDFA;?<p6KNDO7=R=%J4^9tpIy>g(3e;!Y=(qws=5_}`NO22e0;MW57so}FDxTY-3F*j#|F0!it)_&DndhBr31HZ|m?>=!qPLzySB_s&1fpZ)nSEi;fB1`DX_sViS0M6UOAs{k7>x{E>4m@msk1l6JrVFFq#x+VusCtHdjn>9WdT~A<|XO1`8Oy~}`#eATf!k-b;V|$fW)!vqRBe{$2&>NyrQFmD@FtA}SY%?EZd*&=y;ZxDb=@LwqhW=uSoBu6V-`-r|>`SU4X0!<l8URi8+1)gRUK3nfA-p!ELY5f2)skpH0`L_%m~us7(82;GcvnlLeG`>LR)fT^yTC?;n^UjExw_MNuMDo2vCl)H9~Ojp`tK0$cjD9}?e}tRc|&M+Y{!&F|4C4eL(n?cXF;f(bX)nU7}tj8_1y5St6KNZy~!?h6thff>#tGO_#s&AC#Vz|Rg`zJISsZ1@yz<=%e(M={=#<~JD_&yXDyFa0eJs5fFP?O$Z;~E;;hu_P@YnD;3@}nHL;7_#kb!bOEXj?c~nh49E;T3Zd-i6SXZaIj(QI?>Fy6ELy>rz?VJnNpaKngMaZ(||9su>1rU-k!t$scjjt*kX8C^z2p`e3)CcDkQo)|e(DgKR`ZPA4W`;A#i*pL9h&)Qw(u411eusYpCL@i&Ik}pACe#Y2sCID4A{&cx{R+m*DHG$Su=)~p&QWFf_r4Es4m0fR0lQt-QfRPKomkK^BFbsFLQ9;L@ZBw8@k^5sBd8iq-6h{qh}*xl9iyqe*Lcdg9yPya-lpkopsI_<nPT%!AbJk&{@jTbNdX0RAKcvbfuUS~e-T7XU&yc<GYjyabY;K8lqC@oY|dcYrV)`}Or))XFHbVD9fPHv*!K59{?Cv)-W<7fv!v^FA1wV{eUv0>M>~ob^@{Vll()@@Scpo#pE2FCu*nIIx$?A3#uTwb4C1mut~73zt|p6|8<ixz7xf;)ghG1=>ULbs@ocU~vn5J|X$iJ87bR0Vp0*n(HSKTus~$K<CJ*t&-A!4*9hoo!cp?is@23kdkRIF_DF`q>5r?wt1aOvuO33B;P~d#DE@e;RoHzJL9Jh@!i#rtQ;uqd$m<Co5F;}n`0(k64;BuIM8a48ONOlhaGf*I{Ujw;w_F~v1VcSfxRu>nSnb<PQTEw=GqBB9DxxYpJYEzc`a_I7%zdeKviejf}P4Ph{lj7K@`~xeEUrcy^5=89Lxp`O8;d(O3(ETwGznoOc2A&k9$mlRt9}X7myRt0PpE_(nq-OzBxdRmnWVb{umr3AOn{TvF*s_e<`#4^_#|mzRZdb9T3}Nd@K-aRozoeNy;c1QR;u~8}Hm!>$nGIWW`RQNk(JzaG^VIly+_AqO^OOYw-izzYK9p{01~gD4k2C%Xkx?+@MG)*)Zb}WmLf#A({+w>4*FLcBWE7LptgJj@#IHXYh=D_XC^dmV;apd<5Lolx5X1EHyC>#Sq@Z6E5@>)5-)JTno-|C_+-B7MFY1VLY?cI0%UkYaNcqX6Np{6=V1nJ@sd0%d+BGd}Hz!oE65ho*^o3WpJYYa=RjvA>3rxpt_FwAYC^ucALhjVQ`uQ%7EA)H>pBngnw8MCY=s}Fn?Sj5(B`h#=KkDztaU2r4HOXFroF0J_W~3rflshUW|B5;mCR(l)4@wcm-Vm?C-x){>54<<efH0DL-g(==R}AciZ!SX+Me*cGRq;7K@fr91hkd3iQ}GoKV-^!IXn2AO(-73`G$C}uUvE#X&IaW%cp8B`y{s+CZhDOz5t`%wPVB14s8_Q_1Ea|*%>0038kaI9ptd_(8H0J8$-FPVHdW#l-`7q$qgAm0fRBs!Lj+fk^#h%nU&yfvjvI{9(13>=(HNn|i2>IZg<7k=fk~`6ZVITCz*tg{CFUOA3mT$H!HPU|hNL;-m(J^DQf*&V`m^Xwc<;XwOpEwb&<)ZV4Oy}p`+slVF>qt8F3K*x1^O3sXj(r-00NyKf%XMd`u4izwS~xdnih^Akfz`X|G~CPd#I}tBit1EOGVEO%f`U_CCx7r*Wp*jq^D8&OA_D@+g4#iHgO$fXL!7^<7DRGOImN%C^e-HD>UF5(MCMoK1nTP>g+#2XeMZ4>Zq3q%hMorNmzXr!vnW(*?$MjP;E0lB^xj|)Y=2NT*?uz!vV7WXMJ&KXQ}7kZC;-MIA?vUX(7QrS<<tEz8E<f2NXvOs`<uiu47^;+rfhI$PY1+@M)qQq4^WraTME3e%{DNnM@-3B+=OU@v0-7W`I4O(##zciXg;8I^@OJEE;%$-NG5p5}d)1Pk%Q4hv6jsWO*N`mC4oxB|s~AuTpaQGBej0CY6!c#xrwwk}$%Ch`U-zd@l(9tx+Sl^ar}?jksfB$P`5)nGl!UGwZF!AqKnHr><zorIX>igQ4(OfDPY<GU79mMU^V7kbq4?iRa{xdpwG2N~G7!)M{0On?+daqYRoM+{2iIfCM4G$z-8k-$e8dN1Vz<OY_H*yvTb5n8rP;PeL7ZA`wbSPyq!5*~uxMu>PrMWr?pq!bU&raUkYW>pY>}3Qw8Q>4x)%Bx|u)b}Tq-(3^m}Oc!=6!q_asMT%#*F!8_TY{9<?Q3$sZx09^L+mJ8F8cVCYYHfty7#FY#(NJK6DT@uq=tUR;YQ2l*fe*@sx%IJ<*nQyX5}GvNOLOl}jxEmCH3*hi7`@-Sd^?Fm8zyuh(OS!|g|6y!qXtr}rLCvR_x>pIU6_t=;<k`r^46GrMIG)tQIQL2ZdU`6G?40jdqm0UkC_|%Iz%3eUSaNouyjS^AChluCoKeuLp9N6Q?ymX`BN&UX;}&!Aci>;Y{81{f21>{FvjnkhysTCmt4^S6FC&>dvKo=!E6<g-}83y4&JY3J)3;Kid+0`Ai)@~oShr!eNUe)qdi0=H5>#<8tJ?NT=D%Xbd^hx&q8V6&8w?ffEUKr;O#w=7=EbBdONViN;>eAor(;7<ev^-Z<{%gu@J+H&bEY1AHDg=nr+}dFAAD@PuZpco1oW^9z`*AupH9QA5h=yXnsNtlHB3GmlkM)vTBMxAjeZr9Hja=B5?~q<@DiG5C+nB67DCA;nwIra}x*lA%}Y{7=|Ny^mxjd(7nzl;7k`hemjv*C8F(^>y5n7ZW|(exfsf$Za-Ua4wl~v!qq{X>F&j;2|Br7FZytZSU&%UCpw3@8=E1(NLcDhq6zbdFw7iYpooSlPUPtF9ct8IhCXZ_$E9a~9xl(%Sx{<1$^h$p9_D})Qv%ZeF~eJd&*c)=zt2HF$~Yj24$<F<W`O*4vr~pnlpb<_GZWPmdpW`ar8_DtcoE71OpM+P`19q{KKd~=93aCWuZ!~zUboeIOt=l|Jg0)v(W8ec!~bhm9Z6sS+)E4;V+mi5h`MLCyZj`D_D89Vx?9sX{dxzRt86ii?ub@`H3KLF1hA-8I(yNr<xZO%Fv~M#RR%#N2}Z>q3J+*z^5tgrBJJ^Z6MnyB%)P0X$o92^6aTstr%Oa}mDZZ`ghP1mkerMieopFatla7z%}>H-!&Ujn(DCF2D#KV!jv1l;PO&>Rwsc$jA}7_xkK3_Mn!5p2(gS(PRCOhIXCG_4?A`su8bP&QN^*KC8JNHCUh;B<qHA<Q9JF#K#kk9QlSAVfZJpZ>>G4vy^CJ-GO|z|Zon3-PG~F5XpC1rP7`VMnP87%SYBpMOQN!XLc38%A*3);fVRv?i#t-PPzxX6WCTdxGFnahlQ`c_UN%>lK5JD+<!+y70)U0Gl^am;*KqmRQY-1gf`Jx*=3`mouH@(r7()$2LNaHipdl`=>(zo0^($lUQl6;wy<{h0vv7$64JPhp?!Bwb+J5*wSXUOtWiT;(ooeb;zPhXp3tmKT5+9dPTs;%r{pjj^!{sU6@H5wR$+d#!Ng%;UoJT0;L^6}}=8h|wmNxz8qA&+|02R4lGm29?v)t8Esjkg{V7g)g<qo{ZRKt1zmTt5g4y}uE*DZK9J*X~Zo(c&8muz=$dYkR>p#-Tu0vw{LZqQbH6u^iQ&`3X~6E5ymQ{?^>0`~T2B8BBX?oh|sXxQtt}f34e*2Cs|a@9pqB;`INO;><g4JXqM|7>w!iDw}`?23=CsEO5D_b5QFggZ-uCKtpR`b1>r4)KO7O=Q#(HGhn#J%2iD1VT+K>S-X>_FH2nQvb9OcqgDmSunN00rP&!B^k>!kl}`=bD6#ZRWr2NCEi9fKYd-!u1f2@}yMUff6xA>c%6kMjX{rIDUBOscdOWcJs<2|IO-5|qfDm67wwM+JR@~zkM*QB?F#c$ENaSiQLTY}x?_lX-mJs^Rf*7d3pijPPM!6=|{7QN-sy~xj*$NkBfkt=cN#qc!#F+2{zEpZ_W%nBL_PzOi7_KNvUC82+G2#8DIZjMXPmEIDMA0V*H#Bs`7ukb+8))0hnUihBL6d}dXTEH(90=9Yc|ztR?I?#FNWs|rSq^#QTqKPrT5AUYX=_B*>h>6jTV(`E!gR-PcTlbLQZU`uuK?vAQB^1px+bJ!<llHbo-%MomoCqikc29kyK@p68$XLFsm!c<AJr<o#1CL_-+cO=a)H+;Nr7?tKpbf2GqpT>&egT@72<RCaGO)#)9$MT_$DVPApJ5_*h9IZG$!I{jQzg2^35Q1fjwzjVcbNK?L;^1;>my@?j8^+0jMVb`^<}WvW`$JTA3Pc)E-77L>%B3oDO|<AtW{h)ROhf<17vsOD^4Cm}q}TjO%f1&MkW^wvR~aG>mYTr^@0DDYlYH^meXghPk0aj=~rjr%s8UyKR8=(uiYSp~Eri)~{g%J0)eddlfPF1C@$&N8O5I_abyeJV^Xh;|v1O4jpr5L4rZ>Awj|7DurJR(@9m_c=G|1CJbMG_7vY2s>>6lCPh7wLIJVphxZ8Ag`vL--~8s9YlhGd<`*Loo7`OE!E(~@8yHT8c%bsm@jj!H*+`2wV67dy)w~p_EopLmA8xw-NJng7vB??gJH8v>_$)ZrA35=`@CYz)^voo#0oAf_s!Nvm1mUEwqLSh|{|~=)UTg(ss^;vw^i@~p^KgL?8@O%}PoJfDEtBec*?{x9(4Ta2Y2<oFr4)d<I-8*rzMj$N3`A=Z98PlX#(wF?oPnCVdbmLxqolMj93(gpz3_OIgSACED4ZzQi<N<_Ff)WGi5RR0*r(FgPxOsIcgDT5ElNuKwGNNzx1!&W5{_WbM$N$%*O-=lHNaOOnQ_AMlQbG&o}(!9FcX<+y?f780K6k6KkG<{i!7-DP~Q$35l2L<jnoI2%xHqINp1mGqQ$g$RHu7y4|^$ND}k5%ENUfSa%jBB4;g0(k1t)zk>*5KvOPM6b{zK2-gW<7gz;+`3C|EdSoHSvfN49EVC5{qt|RCrvLUEtaU?Q0apydairP8FN4dSlxB?`X>70<pT7bmG@4t<y&7-tU;}5sNT-5@t<o0v67!?nxxQ?nd?Ewq&8^I)XK(FM*#%wIgk|}sfv99B}gw%-(8`5p&)q<h_1UG14#5EFH&Aw&~j5^*qlOml$jth3Q4{va!ue%$I)FX3ODI6-Oo~uBaEhN>z`sv$<Pu@nHCch(?5pk6&4PtbpTk<FzL~XX?c*Se!CJQ0)>ygmlqpFa~U%Ws~?nku9Zw0QeMdgX;wnQDv&urb!B?NvdC-y+B2)KYx!0G0X6mOE@r&DBW=gjxJHX;>VjlSL4Bal-ia)(dOXi$i34S6%V%-l93n354xQ27B)M;~zT-}N(65m_#W1MR7EM(@y5SYQurg67qHxnH1BfC#cd6)qD=b0h4$`h`-Nn_xN@$xo3cG7RZ(sF_ppW*!*7Mr9UJY+tTKs!WqM`!ZgS9ctj}+efDtw$6W?X4wFi0nHra5oQm$4Tb&KR?q`+VdxGeCCOLqB?hkF#oK|9?j|0=Z$6Fy=J6N$|6T1cHer-(v8X$bw$mtJj~fP-z<r%(>IuD}7vD9#lfu|RY`Z(gKQK3^6Fj#A&y|85g%3pniw)?jior1`(!patqt-e}R9}_4aV?6M(*wHM(C`IfR4@k1MEvHwbGDM(>u$5y<*8L&?m#rj%O4baS`+2y^5%$dNE9jeqGiLwf)GKVO-51K^-n<F)?dfK7Y%rpdmUqQA_;+K-h%i4x;S-ZEOqeZ;{VKiM59b&9QOSW7(L&;Lh0n6pR~m2A2sy5b<p&YzLmsR>|jIS!|QexgAunJh7@h%k!k=I?|1Ge+FS|sLG_HXq2%+gWc}dCioAIFDo-;MvpM=Wj#7sVc7yLBHZ78-Lp%4c`}*a2-K1<otO{&d*?B)t2sK^zxi5Ax6g^d*<76%7Qft}_E@0-wR(a~QfnW8aT0i$jiO6>?Zbd88Z<84Uatl1Li?a-6C|Uc=pwXyo82M2SQd9s4S(^K)@zUO_w9Lkq^Pn;LiD689wyK_8Xn<DP`k#6L6WsEGd<zXTuC(e%O}&#@Nvvbyy4nhlDvUWFtf{EFIP}K!lg^oLR4h#7igBSk3qZpdg_9lnO{#V4<veRQ0ZNHM4~Qb1rD+-`=m;a3JvM)Pc>pP`6YWw_j^y0=yGrBjpDi+=z5Lk<!?PkkAB+Rpz&X~KQHis`T@0mS9GRYImc$|*PrG%JHGWY9jR10O*KQ!$mJ2CV%Ks%?jn%T^b4#|NB<Ol#Ed%@^Oi`3K`{Ea*;}tlF>Fg`JOJMhq9Z-}u>4)*0iCHfw+Tp<S3_x2Ra4R=34*X!m+f~%!H73u5mhGZv*cOiv$QbxbZKNX?v~h%ph9!}##C)GBXGClC^2hU=;zGB97ox`kYX)8(%bYMf@mt2ia?%uiIye~PKXg$HxbaHY3u&>stVC&l35nt%_CfjEVMPNkdjcF^u`G><NTTr70t}@O9nPCy{if9=IAl8+T5{8b-PNx+Ycquf%}2g?$fMzKm>FU~y(ci&-Iv!#Z;*|BWN?}j;)_;dHLr=*&FBh+NgFud;_Ma&`aHv;FfWgw@fhVL_6III^%cnHamMPVErud_{36J5di+NmDN?KE>mzZh^JeSo)0|Z1%$4I?@l;~k;GFKV=8631viNJQfMXVOjRol&(pS<Q&69F$w4CP<l{j%e(#w+b1QO<ay$pdqRS;0E1;M-6TO^7e=PX(9mh!j{3t9q*e<K-8HI^96og>$xC9inc`DB*+%Ot*j-IEco<M``wG|ScefON*S3K=p$INC#b*a^@tA(X&SRXgKavI|Le<#Ms7&h#z)yC8A30r|cOE<JzJA6oFAyK%i!!qviN)pz-KR|3%dkA__yH|a$Bt6R{6bsDsrppo-0HI^aIX{xYLyAhbQwZ@fJ1tS$e-|$~|VD|-5p%Eoj!K}t0gOZs<$1p?gt+h`g7lYq>Sj#2%hi`QTH(L^!hMdjD(CXzN|CQt@A_GUg0pYY;UI&mV5v~)4wl6JY+>}}9K*Fg}`d@*04HgrnHH`;u96>TAWGmOa!KvoATF@HBrysqSbuJ7(Y;B<i>pH5PKNLC$3Z4SDh11Xgp^}eoA>NoCH(_IP?7|m@v$u;{t9dNcRJDwQ5o7FU!?qdXNMQoPTM%9D%X?)wN(w<aT+8*|Tl;Jfp^0>xtyw<^r7IrUk_LPWz?l(zb$$h#8s;~qj+Se*nz^erl`cI>YMgGE_qBzk^}dpN1%|X_@b6^Jdwdw)Sx%{hrHpiLQ}}4(Nh~V6Y~|MF1i|XpjR&myB(ImLvfs6Tl3Rb<5rbM=OE1J=vjcgwPJ$&fO2`=UiCZKJfjnH>k2kDCYheD7Ts!e<Hn&HK@pt5f+;qZh2M|F~*Y<;>e}KCJB7K76vqhwYl@K}3SF#Kf=!@%Glqo&}S$&%cfW)O<{~E1YW>PaRg^3g8{tU`n$Sq$PzVC~bI$lowd{>eL+%e!1DbhLgd?-3cHKwICi@QfZAS@kEbr-HTz*^J!ATDBXRYDkA!e=()V#1H<;$uCBcf@ov2<>Apm^E-#x(7B_Nku|&o7YNaU#H3urG*6xZf%#MT|$fI^OY_Dd7BgqRZLlX3-d!wxitoF>2Lg71L%NytL-Zo<+Z?p2yy#Ax^8@y`Si;)C`>3wC5(-#AwRg61qsc3L5l9Gj1RfW435jnOl!f1C3f1+jucXX=+N~?{<SPo87k~%MS~rlhqmN%WK-$%T8qYohD}<bq?cr7ZC=OeJLa<G+fdfAzOWc6szVB@awqlqffSor!;aTc2FTP&Zv6vQtN~FjGNP1QYVMzr7=x-EKseNCRAoZRoeIFqgW~V}Klm}U0S}rcJ{)j!&^&1m1=Fp8dY4jetf)KVXVz-9+k{U<@CYIc^ugV}>YGm{CkHp5GMeM0dOZPte$q<w8B9e~#P|@qX(QTZh+{B?=3e=&bL`&r#TrP~4($?sI5p%VL$UfM+&6i9-i7t6JAII3aFb{!-}crWi?`JKor(+c%1cV}{t{=+;qYil4t>PkV!#X0|Dtn>rL62G=VCSfaAR-6+reCLr#xJD5$;(4L-MMbcnnD|?d6JYyO@$i!;qDJo7?mRjCdr_AIch9mpDEe8jTcqs6z&;;tm;U`RwJANl4&mZnIgivC**PZvU}Fb{DlgzKVEA`M)4dH!!2e-KLlJ5j~|Oxt!=kc6@hk6V%zk^aa^EMuwN+{=nb4bPc~0tP^)FTsjQ8D|sFKLxmxN3U9l17mYZ^b=-`tdNby?!sHTs-Ow(GAHn6>Qm_lFGN3`QCR`Va6jxw&07@)^sgVfuv(h`3`9}Nr7nL!>&RyULCzrBA(!26U9k)REoc#3oHq>eJSaw|=w8yafzO;fJ?^*bKTuCV*!e=_w(mH&scgC)9Kd8q5r4~5yh-0_`MClQ>5J+*n6EA?xk{yy>FqJ*%#QD-5uLL*@KSi-=XHFy(S*BFCt(t(Sba7Uin(?L+nmnyAvJ1;K(tKzhd$`%LeQ-sbm_fiyK^vI4Xf-%4lmM~*_bhs%zbl-~yb61L;T%e&V($S=vxCBS*)*vZc@&4cqiqWgu2J5u*$-%60xDf;DCP_qR}G`0LB?5HQ}hYs3#Zp?Kl0iSGFMEwMM|IKqpUjiC~a)H#s~@be<}t!&ZMzbow2CS9T8DqMflYyOG$@!hM-$w){BxYYa8g%L9VvGGbF8VK??UW7%op?B7x0l^tlUMN&BcOj>|n#8Ob3Xo@_J}&szX1;sc%$cz6)6^#2Kgb{PZ$=6+6K)R*7sx?kb{jt&+H0qycHdLxqW<#{o&eV;dLtm(A?e~jcyqA>|xzq5)5W@}EN^P!6pOIo;PYKi-##4qjpc(Cx)t;6jg%-w<&5PQ=uXy2B2pqA;uy1$$v`IZN?M`M-$)<Yuckr19(ed6^v5;6T!t6;{3>i_|a_kWlJe?u{?x{WeaWjCE-+F<LA<7Lr!WMiNHHp^WAXNqJ@*?z<J5srJMJ59!@T1;!h1f)d<=-C20S6r1P;u{KU#P7_@CzM!ILDf3%{16WdB)eR4x=klUM+qCMI}rxke6)F|i%bqM&hc8PKD1k&@Z#yjyEb`%7|(V>mp|G7R1wh)<aw<4HG1D}R!A+h5!8_t!rko2@v=Vy8ko$>F(`3F)pFZ(4`(RpaYvrl?h$d72D|p*8WLyE1+;$Bdh!Sd+FI!gu>J6e)SgIVk|7s-4r|fODuHkK0Uk#!Ga_0LexK)$OIAjhQLD?-N(^`g)jm~M0hKzjdwCWx0cICGMLefItR!6*5Y00rDw6O^+3R-dZ$l?x3C9%Qg_FI;lc>JnfP4!;tNS$`Q|()Wc~%`G7o~~q_wqqZp$yqYoU36aOCKF_n5_P5U5BQ8Zp7Gb{{iX`?veRiStzR|6xO+56=)?TjWapmDB$7dr1f3P-Fv)7Kiu7hjC~J0#@2kkrNd0xofV!^JSauY;R?q*ZsJ0?iK!Kjv?8$o|1ZWvo1B)C&J>wtsS5ORN<VvAw*&$33$A_me}h-?3d;st-$hU;tZdSb#T{<9Gu#7GnD^s-fxQFXfieOiVv`$O0n9dWsNq1HB>jB+b7NB<1Cqp(#R~YxS0hhbbrfYw@$6VfIv5-_(a_-Cxz76a**Ennmx{tpu^zAO3~}WL2Z?q^I_%IjN?M%g0v43!>!@9Dvt!YnR*h2-|CbKJ6Y-kNGbDNXXI{T~hEM6BrY?~3BB6KiX7?rBdX|5#;*|seS2;heEA~pvt}Y}wR#b0{ZY(!RHy8Z@@`t{m>{ba5PB!17UzzoNgG1C05}Po(Z)5QZj;mA3e_vDAvVWPDl}#P|SwFQgll+Bk#Ju{dXHlp&)9IVym=Ap`ab8mUOXK!0P1{*;kfB>0sH|}VNgp%|OBsJ4#1E@I!Z3?j$A;5`&TTQT$s2N0@hI_<1gWe-a#a+i2GMY4Ezp3oKCD5oDJtkxk`oj{xH@h;&xUG&SOwdy1ZejMHtQ#cE`A@Y5O(;%u+a^);Pijk5BvAA5HGH!5LdP&_3ePVUP-y}VHY(zx4G$deEWrdG*2|T{Tz%^+td4bSb)vhP?8Pn=`2>{@Utb@AJ&8U)YX8_OwL$IjXk`~{B6mlDcn{ONHT`f=hQvukELZ=DPkNe6oZ1XJ-x`Tcze~Xe}qjsX&XMvQ-<5ye~likdP(hlUT2lGPloCCKT|&Z2pc-5@@rbeo~|+~`R{_FRtdxi3BX(BN`@y^TI+eYRX?kT##*c8vTnasZBP4Ic<cN(YOnd~ZFQ1vXb93-1~$~_9={9*aMA}__XE8SzkU#~uy%2crgOdG9two%wjronXa4?=oGNi?uv3iam9>X{{MQ^qnUZ!=&e$T|P=2tPFT$5xcAR0>i8*OqLj+7Tx~n8e08(_iva{DWB>=2qgz<;?&07h93~u!yn9K_so&!{^A-81RZissFS{=)H`t$?oqsT2-`6IKW&B^W=N9-%;62Ns_a~k5TnRa)!NH!UvmFtXHh28G)RKJ}1OAeSbHRA+THQG7$%;}4brmc-?tcQgxVr5(oe2zB2&ZX?$D=v@*ID9J4ElrEofOzxI-l%fanyM?JSnGrruefFlol#F}6Z-Tdw5;l$TguC#9A3(DhE?VUsJ6<WHpge<I;t!=Vje=OFhhs>iz*clomf#a6T~?laZI?bqhW9lW6gC>Qo=o;9XD+-So1Yc?GqkU>}H)N&Lm-$S*rxh;pL;`Q1x_~&xdg#s_b6N&{1$qN7O(tykQ5jxgofc1*vtUr8~m>NHo(=iyTOs(q*dqSK~jdnxRL&Oh|yxs!N~SIsHgpGO01xLrc&pz+Qgw5#*;OIXDDvKd$vnV6S15+AN$R^?Ibw1NDnv!<aM1ic?Ii1v3hAF{giI3G09hXlQSo_OpVk5f6YYK5NA%jfg&&5@9o?{8N_}gYI9c0QyQ_Cc~7;d1cA-_$%`WC54R=b~W@gV(Bl=d-OY7Fkm%LB%|gzBYin*GPnxfOWPY;kWd*GJ9ZSwU=PG`8?P><zGXOWX6vgv58ZUa(*=Dc`0VDNiy`59v#kys!d3!K!%vW)f1O6Xq$V!LT1@(!YWif#j`o<xIY0<gnMcEB(#M|brB{qu5{5k3C9fC5YAzzB9QK-u2=xpn3q}x;hDFKt2p3T3Cgi5pJ#z^PrS)xEEVYXM!>S<{oE5g0ojV9K+BhqB8G8>;#4#k#4%9j2!~f-F5H-G`z<9O<_NPOpQ;nbr?m%ePB`kIuWg*JsE)*$RG7k8|pNj+8qYyYvRY|Cu*|vAgV5x^sTa6Keq%Tc5AUuI4RJL$<b(>w63o6S37PX<}|1sYEb{3AY^~<6$mb)o)@mv*nhY_jNwLZ7IrR*58D<wEgM-<8Y'
PREFIX_RAW_SHA256 = '20ec056082c947cd4fe4412a8c89b6bbe47e5abbf4dc05ea952f9422e527c521'

PREFIX_DEPTH = 24
DEPTH = 52

def transfer():
    return Transfer(225, 225, (117, 270, 4), (126, 286, 4), (144, 356, 4), 580, 710, -14170)
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
    if energy['THRESHOLD'] != gain['T'] or gain['T'] != F(-14170, 40000):
        raise ValueError('threshold agreement')
    rows, initial_energy = moments['compute']()
    cost = energy['energy'](rows, initial_energy, moments, gain)
    initial = energy['initial'](rows, moments, gain)
    if F(cost['energy_upper_1e20'], 10**20) > F(12307053318118561274, 10**20):
        raise ValueError('strict energy ceiling')
    entropy = semantic(0)
    if F(entropy['coefficient']) != F(628119274124457367309274931, 123070533181185612740000000) or entropy['checks'] != 51660:
        raise ValueError('initial entropy accounting')
    return {'directed_energy': cost, 'directed_initial_majority': initial, 'semantic': entropy,
            'prefix_sha256': PREFIX_RAW_SHA256, 'prefix_depth': PREFIX_DEPTH,
            'recomputed_pairs': conv.msize, 'required_followup': {'phase':'terminal'}}
def terminal_evidence():
    conv = transfer()
    m, p = decode_prefix(conv)
    m, p = conv.advance(m, p, DEPTH - PREFIX_DEPTH)
    truth_controls()
    literal_minimum_controls()
    upper = root(conv, m, p, 260, 260)
    if 2 * upper >= Q:
        raise ValueError('terminal majority')
    entropy = semantic(1)
    if F(entropy['coefficient']) != F(628119274124457367309274931, 123070533181185612740000000) or entropy['checks'] != 62521:
        raise ValueError('joint score or dual accounting')
    return {
        'threshold': '-1417/4000', 'graph_cap': '12307053318118561274/100000000000000000000', 'coefficient': entropy['coefficient'],
        'terminal_upper_num': upper, 'terminal_upper_den': Q, 'sweeps': DEPTH, 'depth': DEPTH*len(conv.y_axis),
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
