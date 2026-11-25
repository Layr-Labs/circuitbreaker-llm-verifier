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
 def apply_aux(self,m,p,plus=False,previous=None,round_number=0):
  rows,pack_lanes,prod,width,full,norm,total_weighted,parts,cnt,weights,residue=self.A
  outm=[];outp=[];both=[]
  fneg=m+[0];fpos=p+[Q]
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
  for idx,gets in enumerate(rows):
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

 def apply_main(self,m,p,n,u,previous=None,round_number=0):
  rows,pack_even,pack_odd,merge,width,full_e,full_o,norm,total_weighted,parts,cnt,fw_e,fw_o=self.M
  outm=[];outp=[]
  fneg=m+[0];fpos=p+[Q]
  gn=n+[0];gp=u+[Q]
  limb_shift=2*12+40+LQ
  need(WSCALE==3<<12,'split denominator scale')
  def odd_product(ae,ao,be,bo):
   return int.from_bytes(merge(ae,bo,ao,be),'little')
  for idx,(fe,fo,ge,go,on_e,on_o,off_e,off_o,corrections) in enumerate(rows):
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
   outm.extend([((x*r)>>limb_shift)//9 for x,r in zip(ma,norm)])
   outp.extend([Q-(((total-y*r)>>limb_shift)//9) for y,r,total in zip(po,norm,total_weighted)])
   for at,fi,gi,mass in corrections:
    fm,fp=fneg[fi],fpos[fi]; gm,gp0=gn[gi],gp[gi]; gzero=Q-gm-gp0
    dx=(Q-fp)*(Q-gm)-fm*gzero
    dy=fp*gzero+(Q-fm)*gp0
    outm[startrow+at]+=(dx*mass)>>(26+LQ)
    outp[startrow+at]-=(dy*mass)>>(26+LQ)
  need(all(0<=x<=Q-y and y<=Q for x,y in zip(outm,outp)),'main tail simplex')
  return outm,outp

 def evolve(self,k=25):
  h=Gaussian(1,LQ);D=self.d;Y=self.yd
  m=[0]*self.msize
  p=[h(57*(40000*(41*D*j+13*Y*i)-self.T*41*D*Y),13*40000*41*D*Y)[1] for j in self.y_axis for i in self.main]
  return self.advance(m,p,k-1)
 def advance(self,m,p,layers):
  aux=None
  for step in range(layers):
   n,u=self.apply_aux(m,p,previous=aux,round_number=step)
   mm,pp=self.apply_main(m,p,n,u,previous=(m,p),round_number=step)
   aux=n,u
   m=[a if a>b else b for a,b in zip(m,mm)];p=[a if a<b else b for a,b in zip(p,pp)]
  return m,p
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
    '{Wp48S^xk9=GL@E0stWa5C8xGsBn&k=Znj-T`d4h4?aKdh1Q;|qj}JUQV;Z8OjeWg7l8T84~;HD5q$OYW>UGr=}yeJAx>uXUF@3L7o<7DuEx!s&Wt)59ze9lH@Db%ssfEed0#2UY84bkasT+%Mk;A=6i-^(A4P>c@{o#pZ^nm`!E(AqB7-<}mhV^>eKQ<K+XF#rt9(JOArp$#cRcH(vjKRH%ie$W(U{c7a-UT`QtIs1!6`k3O5DDRhM*H^r(WI9>+I)$kRn{mMv6DhWRYDE?!b8quwnKamXTwZp}KGH9Br2PntC{>*GL}FrA)W~f$e-<40lvL;i@Zux)PXCmkeBo|B(T+6)uo(MOo_Zvp#gTI*<wg!GYWWjOy<#>R(C#=Jark5Kkf_rM04p8Ix_)#@UrGLma>$a?xt`M>&jU5xXD(uqH<g0{oJF+oh1ceT<*<5lXvGE6JD5xb?XXT8)Ytw?wtJtG2!8u&-ulg3+Iu6cEfp(77GV{ypC&0^#prvBBtjwU5Ne{{tFZ;I80jLdtfRcMYz9{wq}K4tP$XcGdhYX`%bL2uZGKDwn0s+lUR0O7Uvxie7th9Fv~V<ySdwY3XQZH+-dyO&<mxJu)fIPMPcd0QX?cn6Gd?5C}&wo(wD!CgoJXPt#+QjO6VbecqI)TjILOWdZPg;J>reQ6W!==+P2L)){19vZd=F)pV7wau;>46R;k7tUy7)UHE<y>6T@*NgTYr%XUH2t=ESU@lq(tc;EZEy!469|LgYb6Qw^mQTBHvXC}|1xxq>dqZi$1Yv(SoIPv*naAQoB5A@%J_}GZ0NGi6mj8%FxJzJOhw6tZgan0!nk!zs63}nP?{`ap;&pcg<ht+^4B?2H%{i2+LGzLe>yfwHCmXDtZauJS``&jXq_gR(}e~Myy`SF|G8!TSZ6s|Z3PVLt`K{pfiRu7<TqX{rs<PsAjpq&&KGgNUp^j90!aKaw=Ete#`y`ehLgYv0(;=p=3e1E3K-*{=zWV_^<xlOc+%*eD?*U0}g5wlBb&x;+<lu^F(J`qdzrzsFyd5$gpwM&TJ+!|(AMpXW{;G>)c@AzT{wY@@)k=VUN@h_jUy-?WI49|nY0t9l*C%z<Rz1(;EsNEEsCjfE;6aaA}Q0A^7Ez~G#b-eU*CfUG)mubQ5V*5VNy*|;3lio}Jzs52U3K-f>vwPiiGQt)YQ@lj<I}C$`l9oAA>vT^3?M>o}u3Ex)Rqcs0USHi^k)j&nnh5-pyXaS8Wm^fH7Aj--Z_xIQ!)*Uu>vaA4cU1UW5RGJ8GeP?wpHpDieS;kLUcb${7G1$}M2L4o7QfmakqE$2kQR`Yq0&XR0yp#7Ru$N5ywp34zCgzUC{Ke59W*jhq^4qw)K*vc^XsD_x7gwX?c3lzq#l;2WFz&*#<Dy=6y%=9uaX0J-&Q_A%pL)dJvMsqx1RDdUiE9wRs2KuG=w6T`W+9zim`C!QXNv=;R)RQrdb6GG7?ct_CHlRV_IokL$;|unYw88YW&x5`0*Z3(~9TnLhR^>7}|-!A6pI#iasvA6#h$X?xOC_iV5=AAgwA^49dP#_UA{-=Iq9wF(@(+_20cl(QZXC9pWP{kUAHkQ&G%STY8~FF(=hD^&hy}c%aXZIR^}l4-^waY~Q-qmCnRrdFs8eP~qhpmg{jzzm<I{mByp<g`w*&(z@A6isycNOMOU?@xInv66u+2znHJnX@>BleI1TK#ruxUExf_U9JKcxx(6dsU1qD1&zH(~rg3c`X>pc%WlZ1WF7f|TycBGypD!RIp+I)_H+vv?kYK!Y)&*`2HO}++dJIrB<z~XPTo;^KGw(#gmY4=A?}Sciy|e49FTJIiaXzRJ&z!A8@vQUw4sYBeo(yilhIQ|9HJhS_HT@dMfv2LLbD%@k^<N=@1p)#=^&g{vv!CpNa&f^s<0R_*?&9Xj5(Ud%UY3fmiVNA)3+5aQd73I)P})oZk5eI~T{Tsa1S|ixZFY%?egL)_XaNSKy&${Zmu7m{^mo~G?qwm&QE<Abq*R=5GVa}ycGU-zAr;Jc^H@)xoDB$hJW1{}u@m`l|F8Ef<HJcC&m{NNY9>EK1StzN0&%LBbQ)z@{u@zDA;0TAgjM$T$W;u#qf#gBR>icp6Q{Jp#~~V+KQn<M&TbBX7%rdQAx$w*s}qL^E;I;9J2{D?Di?OrVXW5(`|A^b{>%kS61H!}M|(nLD>stS$0rPO-<I<3YCo)GR#L%@hUe(KX8TlAqNga7O-EaVE%!vJ+I-TT{dnm+VkJsY;LyMiZJ29ZH13-|`M_l*eq7|pfxT4iZ`AbSG+`+0{p*IbNu`#Bil1Yug8Hpoc$X0SZ*vJ`Ic46Q`mEmRJ{V*KB7o<lTDs?v)ICv&L5MAhvc)s<!Ut;b$w|U~OIkjdyI~Zry<qFR69mG0R|f?#FJC5uP7XwfFI|{@C_B;lZ@D^UDK*w?dFv?N4K9quB{vy8xrjw5ZHzx@=UW|Y7pSpzmfN7RYZj>G#J9bW@%@J&py|X=O7)6`RQKdnaZVt>dmAjS-tjhC-oO2>Ad)GakUO(17A~|v%j<bY_mgKvlw7z@eeZ`;AN)Xoo|8>d|Khp}dD8;9JoAXJIuCTb7BLK3tWn|mv{6eKJ+DsrRCHTcF{pK*$vA^FEYX!Yf*v_Xcj9=<DZ?wxm-gd+8!Jr$D6;l1$fsXFPpWvEt_#ouFdpTaxa^W>NiUY}PxX9(?a3W0V{nByKDZ9JU7xQPU-|fAfpa<hj*d&V`Ca!E{zco#4mSbKinmy}-5G+e-+_IE^RyvTfF(Q62(6THjo%Ovic3;1YdB$pB5TR}<emW$|E}5hWJ#}|_qv|rx*w-+OQ<9=N%o={#yE^n&=-0sw7Z^!6m26gXiIVElZrOqJ4bMmXfSyo$U)we=T|A4B7mk*Ab!XkF5BI`RH$RH_4~xaz?dX9Z)jUN6L!5OxB1|)$Z<MKQW$w;QaCA6a)zE<6&gL8-|sv8JL{<5{6})rYy5cAW|)tncla>gZvhL1n@7*7rd!jD1%=T+f(&qQiVT(6Ray{*b&p)W)pY`9L#F2heiJtybl8p6tv(f7Td+vY_?8zRO%P<xauJ))KX}$_kuG7FKcXSxzy+D}3~`7deY&a|Vz}px^nW|QvaHq<wUJ8<&py<{Q^M5lh0MD-`C{h!Uftl%`&gcr$&W>}aa0c95Jt@wf$hur+v;2Rp?aUU6-bqDKpmsOPJ9MpX~V@$Fwbct-b@h6V<cS?06yE3)duYJrv}LtQ@WZu6jOA=^<t857ZYLsHp;fXxqYde77Z<n76(<yOrqbG!JWu1)?ZcXMhkW3$$K}pH*hvJV-;c`8pX~6S9mo?Nd9zn34votV||SaEzlR-zX3wwwMragN%@dFVguIeL%9ckip3gZXW?Pg+3!$@f@V2A((!4VsQ@;%C{yJS8Ls$Gd?^&=_y(oZ$*|ulR6rL^PkFsM<7MRR7W)DjI$8DbFVN8q+CyG+R~N!3Hn(C|7djdkrh4Ng4mjyo-#$Xfd4mk&0uGXA1T(?KOY`VGblA<BBOUQ8BsfP7QQ4CBl>$e_&jf5STVkx0LO9=>>?eqK>V{Yde*LQ5R9`e^`elqOw{uI;ClzDCuG|7X>$>v62ra*p)4#3)CM(bVzI)4zbW#kSpGdW*KvfA^cHQp-4zXfP6xlV*f6<e9X<guOI0w9hG@Q=IT)Cm2Ce;*>V=g~Z>IIL%&dR(L7?O-h6h;w=KQcLxfTYiG@Q!ge@jWL2!Lv-)_si3<6d%O7$?f0-{{?5(Y8GE8;HVW#&i8F4iFR)u&1FEDk<Xa^d9l~(Au$ox-Nu3mXI<ASNVkWZ6#qTORtvkpR*%O4f2e&uGHsL!?|ZU+tfs3`gsVS>8gp}xnflh^imBysY)Emt=)wUGiPbR|0qR05kk1`WM_5lq0B)0`C)M9!L}}wPsaT6d<{eexH;??5WtJaXtSl0mukgD?(a|S8suJ<JO(>+M{@}OYm{33M`9JWXj{#TbodOqX{3)3t&+{Q>#tY;-nh}ja%0ZoR4$31HbKj(I)O`K2QB(w=X%^iVh^4FsqP_5l*Z(=HCBq5`vCC@!-vzdlxMYKQe27@SB_&^gCWz}SC7n@m8n!f`UGg6-xgYbQM)p0v>7jGmbe)dk$SXcJ=*AJ?n`~k&x{>O~-Db*?6bP!8Y?eR+A)vc$(RYlexmTt8;pK9;PbMLewEMF+^bpxHYY=5t7dB^r{$=WN?V^l`x~-ELvcH$22!*FMo~M%J+_RB0x`-V7-zVYBTcEcJc2((_2!SPU(uZL-FkQc3<Q03!i7yY<5p}maB0Tt2vX#q0OC_~@`M5<&O{dN}wGlIxb(mi$q_rsSaXgiDD{N3LXNJ<AVVg?YiN@8Vf)1O2f2#eL3Wfwm_1pcd=$oQ#sK)x8YjSb(gM7vmMzQpEOQOMvAA2p=`aH*bMKNkw45U#0ptMo_%eqb4XUN84W~~x{l&!nk9uzNkz}s<glj>3d^AoNygIpo*uxmg*%fEoMV^N2ac^iJSmAhe#ppVCJ41+6f{RD8Jw6>!s3>M{X|5x3Q)gntT2xS2n#;;zMiWP}XSg1dam$}j;v5Xca6tP=GD`QJ+;cD~v2%)-E<|hk}49oPF(2Ma2p3d?fkk-LePcz;-eawvzZ}@~d#GmXxK+WfM{7!SazQ}7eyOFzkN>_6>f8wfADve<kplDxHbgc<(3Bo^DIM02NKXJT?{7MraLNapVU*K0V=*CY{xzIa^Ss#%C-(I^$toWDJ*~AgHuqn<x7(;_3VmJlAx%`D^BQ-#An@iWL_UqYxzPQJ0EgRtSi8n0oZri9~7oVAD?pDG2*Hs{o-vD?cukf}ndm^d8IY$cBO4GmqH!s>0Jejmc<_#4k)w+oB;CtQS@JJH2`KW51MBs)Eql=$~{Vd-`a=FbRd<pil(wWw1f~cnH==%WMH*+k2opw>TyM-y`Sg?%qVC*fZME?cO75pe?O?AGsdo6(Tn2AJ6B^#x^(bxyKHFgQ6>00fW#5V!(pHY#gW#=k@uF(Iz<AOy>n)CT|Y_TXlLV#ogNL}FsY~|*t6IBTCApu_5C~c>d*zaG7T&&!7kie7Jg;{MFua|c^Pe)btOaAiP9}2PVh*&dX<>G8pd_`|tr_Ymgmc>_bT_Vug2vT_v#2jpEhy3>cBXF`K+iMEo*8GJ<1DI`w{;+O1x*7q>xZsc@6rX-o@qy^C!oN9lNNpbukCTGL_$c&FGbps-dP84+X<EyOBrbdV3>}RCkrr&rpO07wj#|^vC#Wlor1rEF5^)SxEs4!GB)?it%8a!V*s!@algVQ(Q0Hk3IEv8*76{YZlb2P#a2)IP6F~}K?f~oW-Vol;F<)(p$uegdj<N1Ng}34hG<)|qdwAZKg8!*FTS1c<<Ope}%aMt^r=|Ft1*#s^Q<!e@-~akg&M_(+3j4)4x;%xB%w0*T1K1B+Ng40H;+${CH!P6uM*<rSI0!@gYY+c>jOHib>~;=1O(%*DO2t@EFE@1in?ae~oMWwiN=_^OEp?)Z4ZIxf3G9qI18-I^jwNn2OqK-ldiKGOU@=h2`L({jbwu7kD&8EUMOj<)eoo_G4_EPY%MtaZr2rlmV0%VcQ=H-Ryf9k7?BL7pzyEyFiSJYm-|{CDno#S4m+Aofta2GU`T?s0@N)DVPMzvt#a5lK)?N~AHaB)<S=3L^>tsGwMj=v2jecTUbMI$<2}V}^^T9tUc0vDFYu<J99}jvOF_SbG#ic?rH{$iWV*F`6-AuH)()S=Q8o76`w;u=j_!VuR6z0;xaIA{jP5U}BLIpBxW9g`z_>w0ky%<LTf^KcrJkw+GAT3h5J8whlqI5IzsW0~$sUD;Xh0Zv}X4Z0Wk<%v~!0R7(tm{3R!1-*C`7Y8t!j?y@6LeI~{)H2+Aj!YzC#0AVFJ{1*6r$15H}@~4GRd|=IQv9(6&!Jm0DCVW4B`NQ=ar!EuGbYfFek#69R#te5@&kA*wwg7C!VpHMoSOa<VHUmaHQsj(!dhl=bIC##TweJCfi$)FcFT4nH5q=<Yp_~yVlPbSUa47l)~7zD%K{1u9u?2e8DW`a*&H0`%oq8l8$GW?(%7)DL|H-+lw?$Tp`)B!xPMBm+Twfv9B`1BXDqt3m@MiLX)Gxlc~`cbfhl-z+zzx?gT&F*WZ0nC)6}^G4m})jz|(}Q`uKQH>gh~kkCDUwD1&3)hE$L$OYPi;<LsYNg;NO;HZfEW4u%RYPaC5j+%qLk+hMmVsSe6_U$>!8`*-*<cvrR?N&@t`S)y4__?W~(HOsm6_wqICBBxR_`~)%c`3&~fC)-X<-GC^(u3sX!cZstKQj|xy6fh9XhiWBF9}j`5>_q`$wfyg=wc0()|(L(pIVidF*wt68elCW$sSaob1-BEu^gYDajK7ui(*NP9w}Tb8)zjKEb2IOm<~9Qcogu;-`f9o=zvlVPKO~!-y=<&;9@mi`kBi9jqt<>0f>VbV7b7tr<q&drI01HbU4xOm&SC>42FNmQh>d}AE)cdTFt2KO%FIsvS}8|GK1rDz2)&-;n)g_z|t6{!K@R%De}rtQn7`VYTR3GK;LqCYnBu@0^V`~!4jg1G5NxDj_*y=eds9K@Fi4ubKk3v3N<CcAThl&D(Td*m^@b&OvN#z*w;N%eq8q&S67Np0ej2b)AyC(e%C?rHk5%pUlt#=%9~Q>r)pL+@NUOS@HgI;kvs}h(R#KR$wHRfeYHncSn}xjJ;Z1^J_iQgi~abJ-@-fy?eUc?VUifx9##SR`AY*d<;ehBL!VZi50q?7`yHmwpYrZDj%X332&ON#Qfq<V^Xyf}e+Q)vLXdEfJnG#}xs20qq@&<MH$w0tYKUczIMbbT<pl-@+`t4)t8G5kF*=bdk`;@?MuUG5RlIcR4K@aA@>W^e7e&oj0cJRe8{k>B1baos`5<Nq7_l~0RYdOZ{WW<W?Uy7eOFpME2|Z@OT5$oqNH)=+AetvXoR>y}Cb4u<ZCi}Nf&<T_q}-Kexum;UOvWk*9k}kv_rpsG(JHA1L8Db8x#d%k^g+aZwy<Boc0x%QYlJbllMBzWMYb6ELh0vMQ<nOSnpB&epVvt4xg;(ZBZ-YY_XST@`%U{00KL<dSjrT92&cF^YHsfKfGk{kJ>-yy{2Fv``fdR#KjEImot@&U?6X9=GC<HW4FyIo@1Oqlmygkt4j2oLHb!}csivB73U{yl9(se0R{5O1kgtaf{Bk^qWI5d2d+{*@_mFXd@HLrdB$<&c%U?pezc8W+ftyK_1Pi9-G1z#ZZ!I+o>aBAV(F&X+C2@nBf9V<s<9!rVGKZn6Jigqmu~mrNTul|^y5Gadm{G#sB^a<yM-PhOKl*nr1bzo9T-CQ%)SC#r>;dD8h>+#ci99fz2q&Lu;tKVv=tP*~|NnihMQ4MK&G5v+19+`;MCI;pRFqlXq|bNJ>PGAAy|sC$*Kg-+OjtiL>S{GrEUkWBY7EjEU@S@P0xriJ$Cku8TS{yN5zfJ(+*Fz=2a(Mx-p#gdK7tj7@=XpcJebY8<m$ET|DXYn*$hp|%K)AQ>mxqk?SbT<OJ!C~t3cc-zXwk|sf9y*dv>!aMiad@ip7#dKH1s!+N(agHrNqwBgoL>KX&;h(6|Kz_0WoXtNHW#wrCT#j#=`&GG(O^o#Q7<$b>*23uGk~W{I1n`7aL1Kc#O40nJ#X7p=KPPROoqRHgmtahb|8U2c8@c+9!J(zkE~W*Rc_GmC21syHlyAxhgTo%yAH_U3(%0btQ6LYtmc3CO!s_oeecTE2RUMqA~9Ny<VL;3rtQ)vk3t2CW<VlPL0#--p9eYaY3X@&zmU?A0e1hre2MZO+jKHT%L7p$?ZEcR)2?U9TONbsUJ(K`$~kS)H1Fty0Ag(zjrLVFr<6KZITcz(29n)Ws@>CSN~Z5!_(`<A0dw@){lDfX&Ga2T$BV!d(0o&dLths*;K2wBq{e=Ws%s=)q(TCl7bPTnXe9;Cp}L?#9OkhrJsQx%Zt}<KfHK;923LsOp&x3^Y{ir+?C>bv?9HIR)?r-CuSzfuUPx4#Y0b4XNNN_@*Ww!mb(~URu2okLM=~&2?PMAXY6M5z?SLc=`cF2$>Xb!kZexQe#?5g-v4CMTI8xrM6fXW8WCn)(;V*&E(Yj8Ny_LlA-%a*%T8f?3Rc(8d04A-(!-W*e1Z&p$l|#HEdu}cr?t?{uGY+`5SONq2Tqv|2p0AkYfn@!E%vt9GQ5P#1?M4BO9Av?23Qv^2$3=!ST<jBv_6x3hNXb#%<vrc8JHWzLJ@H|BD>PGU&lS#-+Pyjo{ow#&jNLx7~adh*itHfZ{rSFt0FfTQ5|2t;@M_44s-=W*XP*{r!#@NONpBj`tn%dw|wPc;mH7@kKPYt?*tY3s2g*X|(RL&lqz2kqe2rC}9MaE*1?lgoEft)0)Xb4k811R^zw)ExKFO2z4Hs;T%gmYBB8z_VY+FZb5Xplnq;F0^{Nk%hETqcVfo`L1|rzx;5;YNwWM4VUP8JqJf#M*6^<fIF<K%(ns$RByRofSPqFw>*J1X;qoyDyS2aefSp~XeyYOW&m#VFS!4YcPQ#S)o1s5d?K*jla+IoBb`FxD$Jy=Rx;iC0ahpO*$?NV|?pCbFuzcD+0K{?Iuw7#$0?STR_B23$)Cxp*4Iz4?S0G4OEQB_44OJ5jE@F5)j|k-^qg|z+4W0&RZ8}@28(Sq1L+qQ8EETO(&rl;|ARGcVciCv(H)(f%j%i!t&2{97yiZ+DHV8z7t^0yK!0J(RU5j#a0!;dKuU_``=`)BTfNAso1-c~-VBT~Qc0!qQ7TT>EzOg0UUQ>BtLrvnbPo@%81}YMw!S;Sob8ah}R?#P_7s$p0(l<OoWcY#RwtJ--`6t;8BCxgF2lCC0CaSHIkhwIXJfruC0Y`%68<HN(x;$1Km&~Gl$VFR_FIJ#ZS?~+efMjGN(MYl8+iV%Cq*I`U!lIrt(IheCN+?<t5k=Wl!qd`_6Y}C;Z;$$q`JPC=C-@|ZP29#WcndE<{-tXbreXODbF^32(XJEg7<P&C9yXQH0eRA9sM*dLSmjSFAev;Q0%<j&zkYdvd@RYd%Z%BOL7dqqv3rE3?lQm2G$(ljD+J|o*UB$5{Sk2Vj!^L%^h47MTT+8ygj>zt6`vSIWa0cLb$YDxE7EU;bNSR4W9VsuN1HUj69z_vRZwgjCkk7=;Fy3PFScE?02y>t+-L*#D;DCtY`(CKRf%iO%9vmdlr#&P<pSuT<DHI8dx~n%-(qThva;ys;YC@&t#{qG>1}SzH)AB=DGpfiw6`@QTm#FJL$t>&>8LEx!u~*f7wbWInvvGPi(Cg^wWi8Tt1qibzZ$N4e&@fro%?TD5P6EN={sU4rfxX4w<g&U0hjNww}&9SiHkkEo3C9ukmF%RVIcL1TiU$q(DM|p{@UmCP^zgf=O*YxK<~fcVI;^djLVszn(diCjcmk+^>?n&NvKD-bJ-CuRKQ*~)k;m7lf)$43tc+g$cn#x%bWv^yumAPoN+m|ah8qk<HL_+5i->VDCv?L!@>`kQK2g|p)45(U<!2452~g<pBBG;@vLuiZ4=0TV5ii_E@(*{BGKATU|dfALkQ>uZy3RDtjZ0JIl9ETT3wOJJsY^jibN*vY6BJXlv1_flHK@TJ(#g#*@KP{p}c&qLXKWaGZDIXoXt8T{5Yg={v)7k=A+NR_2nK8FXpDz{|;s0?fc04yK0DQ1Lw-p)=-#iq+UaeNx6|Y0;Y;N!T?<HPHSbo(>ADi{I&QXp+?p?w0(&!hZtPi8?o`A$`=$U^w1DRnzpY7;-sIBrVw=x#sA)QLAlP|z4DMmCcxjS0^WSYhI!=_i6-5pOkze_)30WjnA7LaG8^UGxWkq|%hS-AyqPa?x`hR<0I7&nXw$0?oZIV9J!^2FHN@QR^7t9Blu)}es`1g?vqIfEXEs~t7i&}bT5KIVdJv*2^v|u9JCB{!CmSQuSm~=VP`nZNpJ?uLL5mj=t=XzOul9>gM-*Cr&jRD$xyc<a6=5GnzOSp>n9kg(d)JMjag{YZ^fnspC#cG!KhYDqyGmes+e-ZGxF>RydI0&kGG2ksERG+kGMc_h575AIfAz8m(0!kJ-TwjU+}@*o7<**WN)l#L^Ld&5Az9!O;FKl&>laADX$2Z9?V!6<J+K;@6U%`-ze{ZPqECqK9JoA@X;^C&lkl^C#zp&)F+NdJXOLibD}zaa8gMhsw-uw!GCA-??t|8KYBwnI5ceH$sE#Y2J+Kfbpu!I(+AkeiV#5?_EM&(KcdQ6g9T@@@P@MjLrS3md?yLJMJs3T+xhU#%NR+B*1z7NGSCqfDdGH!>rh+NE2hHmYshyw%;~}}{+)S^}37tR)m-#-GMR_*U%?E;Br~c=>t^1bX0O4J=8&iixJ+K1CxJVaWsf#w_jZdDvN^<4+VV=f58w?2)J>YpKN~;({7q<$HzaKOUEX!Uf2mFO{{{XfrA5D`v?J%7i18Vygx^WGWb=Sfe3N^oAEla?{YwPBHK{5E43u9)rpPJ~4Or*{cdd+qxQiysNKc<ZM_cOQR#z{(|%8u0q;ZCHL#TJa-AJaUM0J8dc4#IC<=R5vUqv~i{y8y|4-oJ)-ehUe9Uk-VCo<0P`Z|+P8uS9Bu?q&zH_mud3#R_WfV9R;lG@og8vW#hJMSIm#LS=@b;z4N&<Ly0;nGPx+<AG(GZ21I=z7S9V!z#fjr8v~RYScLbuG3p`fXD}_dC9=$AW7=~bILr?{N0d>H!*|FC*dI*61}UhD?21YW9eZOC$N!W8X15>W|wT|Pn`H7bnWfWVJAPK2IYd>gO~wMlN!uR1k3VnKSYqb^H=a+@<T%uh#Blw0F`^w^IacmXydY#Pv?$$BzQA11s^CasxZ&6Mfzgn5WMlYpVzuaC-I@ecbsr4iWf&*>eQKsLR@R^x^k^Spwg$eTQv7lbO#wo{G(4rQN0P4Fh9z2G$+;`U+F!6E<ybng!89$rZ)B`k~w^C%54k3zy3iCzDVzr1=C?x*8V*x!`mAME)fpINhmJJHA!TX4rlw7ww(*05X@vQpHSq`5(JMw{)mLxiF`}3u87FO<5ATAg=H=MQ0-_9H5VL$L=b=$FkXDK;2Yw5FzGH|P_wiRxhkEz*Xg5g^37@W*UV4W#sB^6&=N|y$qJj?%J$@_^z!)xW(crQb+$aQ2^{C4om~zHV@U+P>~@1NI?>g3?NiV?X5Yhx{Ne}1Ow4E*==<b)(Ta8v05l}+wo%$TcF89)XHzC+8W<{$%5v$D79-Dzg&h@QAX{zl@&%N8v4}FDB1c`JA+TLn`IukhV&BS+o-Kh?cA^WBW@m8*p!ZR|XeA8NQ%&s^(w)YR%9E@TDw3`!vv$b`WtEi(g?uU(N$FLH+cENr4&mz>y#NhPg?!ZUOQ0dx0RN=V4%Av3v$<8}MA1)E__J(?VBtHKvNizt%56HTA(pb|>dqQYm(ik}7JxM?mX3xLH)ch(*WXt%>KWpRh4A+SCto5K6rQ~w))`;C9Io))nw^W4s_bX86~N~QkV?bwRS~aOs8FA&QsD@5UZvj&KmhxYxv1u+u&29TbFvMiEvv$)N=0Oj;>{Otu6LWKYW?ksxcb*}9=90K;Dd;JYKe(~JI?AvHqf=75r8oWNW;Q@>b1W=3+u5MtWf3??rg}HpOcwu@`X8PC6JCq@FW^NdRQ>RYspfCRZ;fWh#^S0PSR6_uVn6i#r;p{a=(+Ea5zVzUg@~?srAi=T)yQ$P(MFy{V(hRcsO4-hD`KypN}90+v`UA6;Z++1Z`a@rkYD2n0wo>U(MFp@#!Th|2D&vA6Lt_C@0f_deg@cgH`sHyaIr8gI%8yK%gEQ{q>{OiDL%2_;flazSlQK8ldU|t2sipI7NNQG5??migRLogzpnZ&o#UXePg}1Xkco;$ne|o(=OekzurCOIL)EQUun$j3wEPlv2Q9k)bBa5R%xz`ru@`#QL%^8jlTiALeo_P>0BQZ&-1e0XE45jmU}HfyKDq1bboTLjA65-$`-W!C~|=AK9%vTU%kojH0tZ8mpZ-brM1Y#<W6%~HIxCeAUs@6I<Y=Ux)S6T=Y31~5GYZ1VTycvM#Y+@{nyn~0^_hV!#6<Z&@jVZAl8?Jhx7W9ksFLjw`S-<wn?%$o?u3LnE;XVHEA2)G~@o<o);@)*;o=;L0(=qD<_FPu(yZ#G}Jx@+)Zf`*yr-!e7^V##0y{JvujTJ!{S*YL1u@;)$#}3>BNYC*v6v_?Ce2e$+9Y+ZOf-c02Y#M6af<!ar3_IY0VXY6RoE%LNL;#?-EJi(*W?>byWO}JVD>$glvZZw|ow%b6+yXXxe{YA9wB4$s>A2%8oNd<71#{6O5ei2`qw!4X{98)9!Z`{wnpym8`nZ>0i6bD)xE*b~M?+2Xr1q9l@Kfbz4#T-y!fd-6L8fu7rNXZM2l~FTlF?iD6Wj1K@eyV_BBxmD2Lou3kaIE>es{`+SgeL1dKg0m?vPBg1R}I2YqE6I>C7r!RW?fgzOFj}?V_BtU<seWgxR{Wj8l&jZ4#a30okyylZ#))E?pdndFLQ`XQY^z9Y(a7Vv@EdQZC8uVJR^Gy0*1>(_l4XweyWRGCDPtGNOPP%)Z8H|WYwC2wMF5`F5uPzCMwP}Vp!CZFW{^mxjQuE78Zu!q?5z7Wn)zj`XN03a5HEkLif-SXDseA9@N!Yl4FNI0x13uy2(^?rc`PHPhziY_*|C`K{ls60t_o79VbM00qTR)1%_9NYR&wEc#u~MNIMO;<O$_DrE8nzqHhktKQ7UA9QQ288~Jw&y_>TxlLLnK1*oVCw)8(ThX03%ogR5*DAynVtkjoW@oWa;j}>O&J#Ty@L>*u$~BZ!V(82Vu+nsB!@}2a!AdbmMpO#K2hr*rX7fwuHQmKNqQbOx=K|WLl`+P#xaNCUj8XV9TC{;fE5Z(7i_vd@*yD0${HlOt;D}aecFyj+~jeAuI~WytuFS8$@0(YixK1UwH<(Ucqg(<Gdq1FJM%76#KyY`d{bMr&E@^1Yms%6`6LA{8=9wHpDA*64K_rbB5qixwOu>W`S!@XDsm~v+cJwdUn+tS^G9&4559eO*UAul2b>k);P+SlqW1SFxy40(Ry;|*kE4Dd2byv^^u8JW1n5*Si7wv{>Vf|2qcqqD03l}nSFu3(EE&)Sc40}b8ghDWptdy%F!1N@w}P7Un5>R)!kxzG%)CCGqF&tNzG```C%ZPwG$&HQ}(lj`WBb_1Gh<-VDJzZ$L>B~ZY!q+k{X{2bCQnjNDI3gSIrWXz$|}jatV977<&}Tjt_bURKZssNaH<2hDVD}NQOfZ{1e9F;e!~w;d>A;Lx9_Idj-=vCU*6ypgPIiY5s+KrwLy@5l_7#C>F4BmDTngmmM8V`H^;}qQ)OfeW<09G$aI_<1T^@p~P80eRYZ~#O>mGw?P+wjYe#yc-V%bu`dM862+HMl69X;v*%jI{_WcC%#@&L6m8j?g-K1BgZ~pg&tqrP$8HbTx758QN1bHk_{V0y7zK6K6+t#@(tgAt3|_fj^q{Q#mZ0!bUHc+B0k30peGMA58v;hM?e~{{|3vZB#J%7-SQ+j_Y-JZ@)v>dvUs%b2pOd`(Q%)kdNYnA1exh#_efH+-ZTHNxN#fiTcmCcVvyfVBjy5*0szwBvTc{JlnwEzopB_p$w+Tgis*NKVs7@_@A*=4gKE0;62#LB&Pj%ZVm3udTHP!UrGh|1$hi|Ju)j&?#Mo5=kBP=rnEOsvVdRSt;Jktqex=WTPEPJmf4b*yOCptn{%UA%bI|`i;bXS0U#)`$G)mM;RVhnv@g#q6)O`?jHT09;8MOcm`pv#)S&$E%W1c8PA-M(60E*fYAEXL4p`S{0uqZ^t-Y*Su->CY`Guug^5&{GvfaZ{!vOA-aLD4duwr=dD8&sB4e7bhWq;cGCqjSWOB`}}P+$ff*RpRfg%?-X`6@>qLG7k#=pUru-a)%X3!r^(_N-i5)qXcx#+$3P^q?GWo{@NwZKRtc57>9u_ya`J_Ra`8eoJKNr-1lh-yJ`;99r+2Iq@<5&aHbO?1@|KHt(MKFKDR&koN4Z;MzWnKi#I*s!DKuADhwt{Tu!^!2no@aYzHU&zPEI#{O-9-t{^Jh-cTva$`R{Lr*y3mDP=Y&w=Zff7I1!||Y5r&be0AF^(1svD{F9QX{bE2s`+xvZ$P#C@!BCd|02t|5vD>yhkh5CQ7IZ`W7CzKF#3xOfqr#^QQZ83A*KnW=<~sKOC2v|r*yq1GG57i7shwZ<t*5|psHlR>&Z^tIl3OVh5vf}&wM7%Gg#NG<DPOIL7zx}Uk?rZ^MD5RjVR3S46eq$yNA~F@AaS6`4dXK!Cd5_P9tEZ=?!s9h7)uJ2!C-Nsw54YHJn8TI1NSq3BR0z4S%qDE$`*kiaMx3ZG|D3N{`i8-KhSh_1ilP=kYH@}2UD^@5B~rD0}ILt4(F(Kazk?5Jv;EsH7TO*VE}*B&grQ?y=xZwAVbZ|0Nt(KG(8;jYPYXfQk^6w!gf<MV)d(z0T`iW5RNdZ3zs=}L(g}CH^=f9*RWhU2DQ`E$&ZNO4%r!Yjh&0SV>pEd4xFd+Z~W;H2eLy4Vy^t0JVuZV#Gu+?9Aht<RKFID()tbGJEi{(SpwhqND?s%2VAp(u#Ef}EIBlJOqanT5k(O~5oC9{bL(q%n|*Ln^)CxrnG6vK28pH1z}MIi&}FLGb3r}hyxWFW?QI%19%$LR+W+i({Q9MYdwO=Db07$~@jQ9GQ7i2K_zD`~SkdG{^qW{myHsDJOH?LyPg({GuG5Ljh|h`f?L~FMpjMij<_-v{#ZVP%RzYQ<%x?gS!WyJ#S^ahwwIN(KD8sJwE8g132Rt@ib_Wp5f^=U^fmUBnH*0GWgKRwZ8Q=u!_pMaBSxMI#(tgaR^4pMM&nD?lz~h$FU3);wAwM5BLBVQJ7{wk<=4k}rZ(LBT$Y0t;BTh(hDh=7(>;H0Tq~I%5HVZUzKf$s@pvW0*DL2PKsJPaJ!ueQ5UBVHD4sgADLW@S3lf`$!d|R@bb?-uvzvp&|u&eg;zQ%7KS6rw)aQ{`j@HwAIaU^e~9gQcNwZW%3n{O)sA^qku5yY0uiRD+BhVIob8OrpR(pOoF^Sk}ljw<Iq-n*%8gc(Zw-z_|(+<+-e@#Eo>e|n>&ec@Ihu5jLw>P>*xW&RIo0y7-|Oi6BMvMxrx{qE}aV$$HVRHKQce{ii_f;NLr9{&MYACQy;zy}I20FxhD&me8f1`lC7T#s^VyLZ#kPTGT7X8E4V@Yn_}K3h^i04@P{(1m%(rFl@{;1H@?Nb)XJ;Bw5UyoyxvN$>{)ARMhcvc*_FmcH~Qp^D>2Y82`c2~3X>nRf&JVEheH<Ky#-uhf{7Y0{$iq3;E$J-l5S(V-8cCZPmVpB<z%#$?35{v-3{Lx@$dDV$t4<sA6txivo<nLg|ojZtXPIp@RrAng>0AHCIVM0+AG*GGcMSWOotFm))inPpzAra9okt0!<q*THW<O$m@d@q5cIp};%RQJ}~sLxbvw1K79^DKC9lr~WZ~q%zXG8O#pItl}|pgq-5A;w6IfSkmWfAgiIqYON0ctxAq^ze(~+!AY&VC9;nh6d$DzVZ~-c=JW4{r6vg15%P}Fbds+r-Gmf)iQO89ALiI>eqYbs=1L7QF_@gNr%L_ssT=8-tV|Z<m`<G)Q*;i32sBX}qFJocutF~@!mBn37-n0nMZbz%C3&3Ik?g`m8nxjc%j!!;VX%$h1e_m3xujDMT{{6WmlWdlHc8UDUTOtj{YW(<4h|s;CRcm?)!mN_KRD*$cr?7AF<C>F51SQgW*lqR(mE76v~7c|5hnihJc^J%gb7egzDUUhBFj-i$$Hg&nKt6*_cn+m_<HEo7DiP!*}Ymf*Mf%oU>ihJR66@;l@3rv8$j#$iipc+Fft{FdcHb=sHvp?5ACW736&99()Z=QZC*Y3JY^hl`1L=$&yGLNI^lvm=qZ8Vh<4Ak>4)EJEf&tGiXlAKa=tR$pahp|cI{tT1V|yfz1GG@PX2aYUp~`ZWA5DP1iv#5qHFfFO>Xym>h}}IKgZCQJHmgvQ)R0w5af>*isL}vLD^>w7rt3*PPXBfnw!pO6kwbc=5`unbz%5`4|>!}QaZ}Vgc0TEzHr`<ZVmUPsW63VK<YevixnBP2HPg)kTK>Vyyx@^nUyS|w!m#l4$-k$@?(dh25_h*)(GdofoRf7My36pjMnZxNWTUO?S8{3>u$Z3Vwu-N-qo>MYX1KYz~pZ2rR{cs2=3>#jPl8Fuf_wYF^|dB33ig!pMjZN;gu;UCVZNfoZ>2(w*tg>w9bhu-DHM)XsYMbuC`95+w0GW$lFxz<mipoB$XpNZv`R~uKs<cw};;d{TSLhafDMPVr#nVt&{XBc&v#)mZNkD5_ohhQ^;&A)93K8II@I`(q3Y_ei%OCE<OkwqYL}pm<d`Y`kdEVLN}w&FHWQVRhBqSG~5LZO^P1hO2DFR)%M&}CmT55&nmZi_cQdm?aZs?PQ^a`(qU5+r=hIY1X?+r=AVB1%y~UhuaBF~>IY+YN_FE~-hj1mhATbwo>mP`CVG|W#gk_wobNOv#}N;fK@1B*q+c9~a?Hhp;r`_SWd~-ev>CQAW#=a}Ie#6W`Zh&AKRyNJn(7ID3$NhFzk3S>oRK)n{>#W%noeal_==KQy;uwvcpOp7+0q-7<-L%XeB|}<0cig+OGe)b+~3hr@{9>8SD3Jx0@*x__So2wz*l>?#LW_EQKDg)bmeX`jJP{;Knz^U904SgtNndQ)3D!f$x78mfKISsHFN#kZOqcb`%?QX0b`nZZ*iOMTCplJ3B6pLCd9k00a?z7jBOkk8I(#nbE^i7^yV>XM`FOrP`mUtH{i_9O=@hV0x&xVxuL><1H<#*DOp8K*l(mESm%Y}oM}kMwIOtrQaM$jgCugFTF>1-LG`T%_pv_yySRodt3VSOG)d#54-B~Mlig(ZJJ3Ntuyv^a?nNK5xdO>SX-PIrTp?f30r<&AxpaQO7{PsaJMs%m!^<rgqY+B7WrBl66ty6%sxMIJ1jgSV!RKZ)zOqAhXF53PrI_AS*9z%Og1u_-`D4`>#daG}x!4mC?m1lE=3lQS+UGcF7#qt<_xeggK1ZAp$I6-SN|qYi%x!-ro$v0iK|Evz2i|^(+T~3c$5xoACV#uq8;uIR%cEdgAfdZFoqXWs=A<d3q@hvCr-`VYqsM4h;SzWNH}y@;lV&bRu{HfXxopyBwKNH}dxzVwgt8kUrCPN~*-f#_bK=&RKn28O7V!mFH&uK0Ib}S$xfI;|k{grF4;D>uk0dbsx)<qSIZu=ZN@m}bg`T|?Veql=-pYyMFgriPyME|3V~XpjC&{j0(4$+rtlT~%?gldE`-ut+m>3+jz`|qKFltslC=%e_q}A&JkSz+B9TtFQA{Jb=)aW`hm9p_QZ`@?QvBXoyM#{QH@g{2Pa-&=9Yv~z_ro6w?xD<oBJR*VRa$4lPfkJd1|K}v;jG_oxIqZLC95y(}!(}}*+p2&vpkO#OC<^dd2waG8O3JFYxIagzI-m#fOZRkN`cGTsQMz;ES!+ja7<Y$Pg3v365X@y^Z;!2Qx3~ueCgMFhDOWW}8v9WHS5atBYf5Chx#3ReBGkHjg`&jevhw9A<(AomXq%rs3oh1tU~#_^er)Xzq{(k136F|KPFhd6M9q^g-XIWuYI9rcBB)jy!XJF?^Z&w=X8D+i?calYl+;58p<m6VGj7s{2OCEk<#nmLhpw1m1JeSn5;|!?tdHp-n6Y8;LZBQI;~6TMDqrD{Qh-Vem*ljtEVdJ7eU{d0BWy7U?H3hie0QuP{IXt%c#x1Ld>Lvkxae$H{tte!HNljd3`t4qgJja<eJ|wQQHNS3ul7g!&z%KB?t~Ce7vN{aukk=_&ayg2t=z?GFyR1u<xDh5xFN(L+XQq@IUm25Op)gb+_aol&Q%P#Y-g@^6kO#QpT9Babq*tu7&$E1S#Mc~7*5k!Zt8OJW~NWEPEZ*+*%LV}5Xb?#P#+XFmFANy1~E4X0EM8ohLd<BMJ<Fs2BnJ#zhvHI6W*cf3{%44yH~byy|@Exf53sUsUC?F@HD9UCVqK}`Jh5Es$%ixllqm-ffS?C{0GJr2muPPr&@nAW05>dI2_t13u{u<%Z<KBH4u4!92lfcP+8thWc7FZ*ztj%ajRiKdS&o`JiUH_C8dy^rwJS{Uqco`?!n}OjV$;RW?5KE1bU@64iLE766;0~3Dn426@<EEN{}A2ML8F#2<19gb)76b$B%KEWa~xms~->@IxYm}&LhlocC6RzBYF6$$u6t=DX2Ie>R;T%FLAELYNfs`L(oYUfUGq?*v34_=XN|CeEt9NpQH3|aV1W5M`PvFH;#6J|Mbh`6e|8EdCO@&fvJ<CzQ~P?0k?{3Fy$z+ny9ry3@4<6A=DOLT|ZpxFVdc~YoQC^n{3OSBf2L;)~@|PI}m3g`A2qYp7B$BCZ2T~TB}2#30~Kxa9yhyct(#nKD^(P-Y}WJq23-oXa3)6{4A3ArGsvIL5fw$q7lh-ZB&NM2bzpFSia7=gWK>~9_?48arK0QAL~qbef=S{IobXB?)BF{`N~=!rF2@@#1E}a{^EVyo(Zc!3k?Lrf)R&WxyJ}IyGHAANH&`4Xtzgpi^QGNpBQE{`nU0?$`C#}F1Uy|H|s7Vn^7to!M1Hnlet4)wfw_NJNy0(P4sWn@Ff9?YVeEK$Nnkoad2vWb?oI668zWKecJ}y<x-78A2Rp&YY-1m&tEXGf2!WU)X%GL-Lh`K&h14jJEE1Nhr2E;c-whXq>e<sX@BsP0)tSRfVE;M>Lz|Dz~%<~Wn)(HTrl4P&Y{BJz%pzhKHT5zzc#7XZ%?x1ogW$+-mdeo$&9}yCPS+lV60tD^zl_{+m(t)B3y{{|5Dh#tIL)X@S`=~UF}OS$y+}($Z_fpoop|45Bn|iCF|69YKg_*%-ItrMX8%6;K>{{$h%rAFav&b|KG-Jq4GLdvfK_s|9`4r!7ueYNWxSqS$SHQ38$!EGXfyo`hj>FjIJ(s*{Jf1K(=p`wDsE@1&wSvkN$;|3uS>==%(sKMU&v~4lsj`k0k$>@e8o8Fp2)Z`z25o(S5n*4wKg$;JIT4JomGi$i_2k_Fx0LmY7AzI_}nTd^xS=;8~rDJ84M~wwYf5sa6K4M1NgJ)#Hp9!z+I;pr44<i~X)v8X7iWN-H210MGtH-DIR<Ww#1r%wmysv$g36cA8<B!!5Y`=AUP;Lf-@5u0p`=$#%QHEfmGX{B_)^;0um7peWU;0me)TgU1Lu;Ja>e2Tt`L6k0Jc|MrJp5Foj)c(>+H#rND{XI)%VI(bpMVuDOQ!6v5t>P+D|JHyU_@Es8sm7!Y+#Zvdb%U2Mo1EKGM8_GIaic?s*;QbWxaBRcoSuD#Kz*}GH^@N2$;I_CTJ9&I=Z6q9Hk|q;|Z|F$2h2Q<|KTS~mQ{|630_c|}LUT<dt3q-!#9l8AP>duI#P>^e0-celf*7$&VPb1Dr`^_h5x=wijfXf81tL3LdOT<x+-74iT!Z8`k8j6r+o}9HSQd9%!ZuYo9`elRlIzDa6oQ7(UL4qd3-=(u$LZM4x~x7M(B+4!P}%S&Rbab<re2u`)BV|(PUy5GA)h|+$XmaNi_D;`hK3Dru)j;#^a$1@)Zh{o-no*gaF1SlQ9SU&%J#umI;HnJ^)xV#yUy__Q$az~?_@U5qyPSS1Wlw*hb(zSf})>O2phEHWn^sXF2@Z%Alu-Tdg>qBdxfPisJShXqsr=zEOTPo3vXKnoYby4)@D*_$DC$n-nU-BB?#fFtg^y_WPnk*%=8eZA(z32^`e$BFhpquqMNlnP-U)6@pp6(Ta?ZZ)6I_`wS+q~XU!GwNN@WfrYsXf)n8PfLU;iYd-Zt2xEiUAbjlbj&<E3+q?UCMwB_d07&vQ+$+Sb9>IYwG!XFnn>E|z;9ck3#*sL1y$ZepfE&KaF1z*1UFJlh7un_S`IA#N!9i>oQyK$Ufq1ttX4LzCpWG*BmHF2ZL0ebF$Ehw02ArFVWDsCfbYv8&~M?`N<&IPfkIL8vg=FF(yv4y}}H~DTSvtCNOxSGMR18J%DU-cDSM9ptDSwEq`Fz9wPh$N-uMWMeveyJ&DrQ!c0Rs_phZDmN3fFgOd&y?D%`EvPlA|#L@K*0G-F_}KG*#V%sbf?LK-FuwFz~I)KeMe~aD8i)3^Bx9|&cffyl7{l#v7D&F$pmWEvaZyN1dURG<6LN?bfne;1RSVf3|trc`)%$;ZZt?mZ2Om_?=BzVh>QmCU8pw0kasz!=gYJJJc06pnSn?3=~Ujl^u7*(@=5FcPEDhHiS!v|KqymX7Ai965F}F#=cdx~C9IS9L~I?Xoi23OvZE_<4alN5AF@*{TA>KEUMr}u7~;2I|6$*eB6PiTp%a(md7ak+f@u)29=B$9Cr5NL-EJ{lBnnh@1D(P#KU3O#?hdW<_%)er7fsO{9?YY5E_x&^b?&A*yC7Rta9RF0H!Ul3-6-gIxRQa_O%d&*WL;@s1+U~ckS~Gudjv0+!M4$$y_(`Z<>A{;zf)vISi<V9$`bVT5(;qT7n9$pLWg!58hR3Fvo<|a2LdKc=Z#vVa%+)4bT`IR_?mMN<A{;`w3(53!ZgKb@{U!6kx9^*0((#dKVJUWZdMa_>o;(;y7cuR%v$}Yw{i~ae}Ke3kjI@MiHD!+S$ok6ih|39CVNK6Ox{xt@AOo`I#J|?mASz!>B149OR3%n=xFN3FFJ<iH;2G(DDrUY3vDo@&^+lxaEg?2@Ug=lP=O9+8Kd|TvjaB+k@*JMjJhGntq;EcOFzb-H)!y?CoEs{e?jl=>Bv7zsAot=0Kv!TCq51~=}+QB9JmGy27BuK%(OSQM}y24#~~(N&i5z`OK`^~_(1l5!V3=7n<XFY>#zro8=%u+qDvkNR|+a^_=5|JQVpe*K!F0BL(jB_H>v4Tp+CWdY1EC?5|#lPlf4P_-l%Rhp6MHklg>`uD}4|?ov3eQJ8?r6A|U`!QgccvBq6^poOIEO;cTzm4{l_2ycr~ph;9KoYWW^-akJ<miwf{VTHfAMFVr##<82e_e3!j&F<5m*7u8yK=q_d)xnlI90jir$`gH=dH;ms_<|&NTWQ=+#G{A;4IwzO{cl}*Bn@6DQ_&=EWAXkCf(AZ-`vT?)7d{F$Ym=zXQ1;z3I!Ijj-S!7#VzZ0+4#t+h5(C-{a7ck{iQrkSI1=8+fA^^sSjI%5#4g>5nRU~WVw3d3DTyGXPh3r=(X$xu(E)ta#;3_<LM-*6^`PRK)pMkvU9xl0>a%sgI;MZ=IJZ>RgE=*!g{v(F{V@1OLTdkE$e!IGq{xo3@i#{7sK)vz*JMb1jT#UBIumNudi2mjUj{fYoGD?ZIhkB2pxdZ_i1c%EqoIhAG+;rLLa*l;q>zq~uEEQAZlv|uk{^}!qZCNNe#iy>cHAzJ5L;EV3eNVmZk@NIAIHDot3HcNP^X){iTVHgY6BN|p1hH7R!slm!!f;YV{ZS*7gj-cJuKX<OotiUyFH<u5B3Tx$lCZ!3pG+=~pNqhFH07%6jG#f89h45LA?~)Dmy_&TFqIMI#bktrk+;Y_`Kz`C7jvq2PK+d2U|F;tLUrGEGzLn4o#C6zjiQy^bc%;^M3H=5=IbzGuJ_kbahR+HpUNbiNYk|x{G<Wvq?h@m#6^5xsx6boY5W0ppIa<qj=zeFxCz*AP=Jl%aHHjjKs`-dY7Ozb9KoWsHm^7<scaK(t5h{~-ga!=u1v%K$DE5qSVRSFhsYm8y1>06k?T{sAx~$dhzK@@n!*U@p2<jeb8WzSa1_RP);DoD)DviLqC(Z4l7+U|MnSQ=nlIkFA9-%4%hDIj4<!>dc}|L}0btT<ycZQO#vWvaB<^M}`0r!>AD86)>TM#*iSYfS@(I>U5^n#X+tXg?Oju0mAYlw?dBZG|x0`ITbKt&Y-+!lHeauB7%6)^Aow<pNFEJF31*H{Yg?j5Vf8m^w)g&w~u<ZYmT4Dh;fGako`D1U0d7Q-Jh}qJA#)T$L>Gc0qp6mhnp-Q(Z-V%mte(KS};&3?e>xZNkeQjFyQD(?54hk1jjxmF>YlWa+0H0b!wr#jEn^q`T*Rkl3=iGpvqiE`Xoibo)|C4#6<v#<`gM&rmtFUwHWIrCB)OB6Hp!-<qL1q{unalNfV@rXxwR>rvh%X<&O|47DkUnkq&5Y0|O@FbfFb6FHpRGbPI`@&2JS+@}7zxV$I{~uHFsN)n%&mAus2Hk(?32RQ&y<CoZt@$Ke<`X!@bBU_&dJHE40=bLekyQyaNf%hET+7$EbPhqSTy&~D@%&ExdiH>GlhdE$0TMYk8wk5u}?97vaGUP_>zH~_-Sw_ryfk0lCuV`Mn-V)rRiDVF7@kpe%L4lCn@WE1&NXj1W+4Q)S;b&DO>LCCG4rA1<YK}Ikq<k*9=8$N_wHThg$vp8iR*ShzN}$eO3oobk2oFatNCX9(9`S+nLa>+PGrLiW~-I8|-{n)}prZ?>_Y^Ro-{hoflq9=2y~riu`w!&*bBOu|jjWd@_<x((`5=9CQV2*5i*ZI?V&e{;norjW3Lwk3Dnj+t-Ntruw>vGsDLo>Nqj#82e}WZi@vm%I*>)-fSoGI8}?x@%ntn=RMs*-+QbPIxdOs*#BKI{6SA|eaCHuRt`H2GQ`FbTOczEgf9}{MR5-PX;MGiMXg7RcJn^U_Z{_($bnRvTG1o+yf(GjLXAhxRxf^VnblCQ=xhGKbt5b+r9l6TVQJV{*)uUEU}deuSJU^6u7eEO?*$-to_`I=w{rsin-kbn2v*x~mFUHcA5d;&w%SrB$+cNdj^EQYhz@iKY0`=nn(Ga^8blu#XeUQqW(E$&AP!m+*=O+eBo8hbV~v(gLau_Zn~BRLlXzXuZhu+ky9EL!t2L&!xZ9gM(JMxCd^<Qkf3(>>{#S{yJmME3(E<phAr?Q9S$h^%g`L+q$PA9m<i2j!>gfN%ySPDfq-+zhax*cg2fg_VrT$U0_}R;6Q76i99bv{+rBM^rP79Z2?gRo9t`2=?&(!r6b>&zEgo8BD+Egl}8?1czQ%8j!vA>1mx5XI$U8S^G2m5-{Un)Gqu`k|n5@f-oV$FI9;G1W2#yF{aN!a!-sn=l;4C1SNuB)7}2%8l87sQz}ZHg`PYf;a;m_vyP;gTQ=mm%F$TgKBTH@czmRDogsH$cv<2Te7V6h_rpx=dgrn+T>!w3nH4+qXaBb=P|sRB$X3mfoW9=(RkyOFmguRZnC&!-g`aJ$-RVdZiO|>F|ja73_3B-6o@oBv|AZP?k0UhELHQG1n@QpEj`;9uxEUA&cRmZfOHl5N5Y*wgeT~%1$G+QyeC`8^!IrhogT9v9q`RikaNGpC2;>eRzCgz-TyR&<RAke@D8N_G_OoVTtqiT%<;B48hfMaACD5k%|0Vv?N>f(l|*ZvRrgC(S?YA=>^nQKViLO!0fU<#=$gh*xO@&hQ}=kR$1+w;!h<$ERxu{$1$^y#=a)w7dx-F&2GIs_jRc8H1&}7Bz=a<VM8s++2i(4FJW@h)6O&SYQ#3{VV$rFmE7icToesr3KZnf;>ILuvpy}@);=F^5a!V`Ix6RzEeLgL3&jV8bq5r{w1zmc>w`HD5(~Up^WF{}r$W~j50XAZ9y=b0quHpu*@5sUME-loxpXo62k7C?2Ju7YpILBNwZ0$BFfG39rfKfv%Om<n@ko`!99~u-i>+-dOMJljC|7Zcbp|mDv2m6MlFCIpN1ZO@U@X)$Tu(&ey#nUT$ZL3J_X50*1S&`mX;qY^nz~vQEtmPfyBGW?&VO+|XHI1*vuwSG!lkZyQL<^witRB8bV0j1m{<dhD5<%C$GA12HJ@nZ0@SjD4an><k?U5Ymr0@i_w=eV2|~H&K_<vnAkN8Qk+3KUPNoTI$<>)mudaX$nKIr#y$q^Q4(Fv8eb}uS#{KRajD~>qJbfEXJTAY26T-1t^KJ!X!tCZcaw|>@@xV@e^TRi7>JA^?E+T>>^Ej;rCjT%VT%!eGZ7fZ^m=5mX1C}2~f)Z>BHA>PcDGkT(tC}G4MNzQPj0)a)d<FAPeBwlWDCHc)oYC}dcxfz2Xfi`7&3sRlZ#RndY|2%~*e1rt+PZEFInSKl00!#PhOOg{Y};7NBuhwt0dV+Eer#(TZ|cw%j^hE<CZeD~4@KlGo`WPg0<jD}KxPrMdpnav%*Ln5vGs0ZiML>$5%)|+C)=P_zWyf!Z~Ph{0`VWHX~8m*A~st~FJzW!Bx3c8Sz0>Dag})<2Tq=UXt9OoKbF;Fi7G8NP9CS*g(fO!4;tZrD?olB$+6du7mmnfmE9~DjW3~=<^w$;t8bc75At6>+U;zMQ0Yn!&|Heh8080^6yTq;3L1nuj^#sxqdjm49l{Q9suUwM{85v@su^b{ZWcA)py-6zbUe8EJKbo5>o@E_W{6<M3iLz)<eh1Qa!a%%e!5vB-jbMDQ_~0jEEi<9W1w*P!C$TLZ%@z5^#3obif)*&Rb3gj?GaJYEWL?gag0Pa!8<4Im~;`dtvW(5KYNyg9C7gWa;7Jk$A0RNqKH|_JJWNlF5&f`!Wa1;x&-x@$fZxo6;=>0uMbRGbPS$maXQ!m=ut5EJ8cV(^AwAy(OnDSIeV4Kj4st_uPMF<L-}o^s&NdnGY&CsB`KpJTXzo8;x#^tw@I**G&dcN5#en*QKz6Q24fO6It8(R%w+5zXnQZw9Zl*^zAbah-&Ntk<#(kJHL<mbATn-$o6TP%+awkN29KaW0{u)oq)m(|4`26Qh}H(Kb{pP~A21POig5ih@6H9UBzcPsSiq=N5vO}RX!RWZ0<`m4HnAqYQ3!Ui(km)TxVZ<M7f7gCh}A<lI{XukUJA;Ins?|~o|pjy$q47+wd~n^&?2SDZ#N31nANVo7LT+!cpU?hBmA#40Ua2Kv0aM9b$$PJ99#E0-pTJ6Vp^Br024z)P{83r>_KvKZeeTS8K1#wklSto;;XQIQ@I`2?^H)L>DwJl^q0TfpZi#|cEK9#n4x0(MI=Kl6rZ{0r{2X(`L8s_#@z+Vkc$a01yNb}z1X37>;zN=sh=bH%|7kG_I!Gn2mO)bMy)i^N)zR<4S69#EX>WzP1w#}y>ftp(FjBjj1C%iQ&$O6dg*5ZHK<q2q5{9e$q&(e0M}-WmyEyMGOJaT{c)S)Jv9K(liHgTp$)?S5y+%4y7<dwU-->adTP(gq~$SGa_qNE76;+;UKbDg;6UY$qyotL6j`Tvy{ZD*k4)Imolh9ImzUJ_dNT|W=<BvpAC|a^MppB+d?zg`a+0FhM)ES51u?5%9rbrVO&74GuU(rx*T-EMcoQ>8GmAu472=oa@&h94cMUsxTYF>c@p%pU4ytWMeP}zkA*3Fjyy;aj_9_kGUkK^Q@=-o(sUWXxt%FXq*$0AlQcn+@z-DHm#ytaxkXMw{H|Wy*ps37W9RsvSewW;2i6d)S!ZQb_eyHp~!5Q2iQ+eIh@)B+1*m29ST^LSUdB;|KoGl96#ru_DuNr2Iz*QuAX~8f{_VK3?znrPW1FLT4sq<(!YYRz%@4nGDCN|6S@db+Ml~8p-46ztzppBX`V4jm`_A>Le?~R@kMiH<fc7->BykG7^X46+NEJ}ifzGF#}ng*LD(ldZ<hiNmK!w5im5*tFesRc?YF*%|ZK-4`xD4zPh{74i!WeJ34GMuQ>H7L8uCRCeFo5PJbMxQ{_KQG*SqI;=vnX2GYszl}?xd;|-{tqEsIYs;oWUI{hW1SGmRzd9Kob*pu3^F;MkH_oZg0-&;x<6KqmG!mQ2NEkS{9iZK2M+?X2OiRmGOjx!QrkRT*!x{sFA|>l!XwHq`3N18)Y%S=ez!%fj=~Sc;LCYL0jbM~ry6nA&{;Q?N5XojvB6n5^=y-?8N`)QX=ygDO54!AdEq#$M7N7;GFMzH#?0dh?SIGuTQK4Htv@l(+W%vzW~iMsUwHx=UgZA78=ZL_dmYcS*!kax?aFQ~{=Kyryymk`kR=_`bX2trl9mMl*7m~%FVjiKaUix;UNMSn^dxknZMc2J>ZHJvEm*rlEbc#kG!ry(@Rs@=1}uQ!mCX1EQsI(DW6)QI*)uKHv^imcTgwLgekk}rF7I`*si0sUeo3pp4r5e=FceJG-P00L<9y<8T~%yZn$aSr{Nrl;h2h|2<gXP&?m7vaG8`XmP5B6>p3_Ptwyv_=10lSJ-LtOIjgUfSmIKXO%Mw@ttP<|j+zYX5eB?q<L;_0Fgfa9?8pFWM6$K=;LraY-$_BdOTP6Ls$l63sKsak}uIg6WBh>4)Ff_!6jdGKDdNV3`1^QI#M3E-uAH2RLWvU2Ce~5%%5)@r^S6vJpS_itIc2Ur#Y52wx&K6t!k#}%%@qV$9>E_}Om4tYVU<;cGqg6^SaR*0cK_vBQC!>jVbr@)Aa%>^N_1<d!WA)661_V$lG;m7z(9iYf1x2iX6j=dQ&#)AG>q+pP7r7ULqnb&y?f5N2GP}RqznDqGgD9H_?fU^5Cd{a}x4_p6M(IVF?n&iQ2XP=Q3_P3cX1%~0GH9n`88GjaKIhlnQ8HdnYMui0WsGR+r{;!RVzR+o(B4$LfCfY3F(^^9MYgZkkk(DnzyU^%0!|?6m{sEz>-|sqkf+ly?a0Pw&B<iE#$G<IvsjH;0g9N_TUl3QuisP7sCT1)blr|_%GErP!^w(}EatX`C-nKsEbE6O0RF4XNX-O}0&OOdvGB+}V!yMI{Onp;AgDAnYa5Fx{qSQ<VjAw5nENAZ^=)g${{s^!wWk2@yAt!O`0^;0n4Oy7c9J2iS7v_^)XTL0e4QF8oKw_;w2wh$sf)@xd!_m=&uJ+Tn~G@5uU@`=3q5K6tWj$R206Lgi-eNxva+H*bemKtViP`z4WWPsrY!cPi>LVtwV&;^f~K)S%P8px);}L5mk#b!Izp&(Fg520`ONH3^1g>ov$a`gmgD3?7|beH7wbF)c2W4!sxMSr)DG(5)#S_%`nPu)<Sbe1KliCSiMDz75L<~X(y{$0n@y?LaXKFU+r2E=K-}5$3cjo`3DL~15?#gGf;ED{u+O`8TOqz{kS<M{%-7=SyN9XZAphkDpWUp+Vn%E&*NDtk{yB`(_Xx2tN<1iCEwj^je#W@ey@9(NT_E%!8uL%(5n6Iwk3stA_-}sCH3u-laaD96Bz7>8!{0y1sU!s_yxe<?zg&C*eGIiC)npGOaz?`>@*&XTgnb|r7$MPb*2L%ZKm9CTm4?-ezaY@@G+rr+GvYX0iF~ElHtaUIk+V*ASvkFlWRV_Cz}=G*phkC5+(Imz|4Btzok8WGe{$OMX&)tBN@RoI`jH|=6hR%P$`Z7$J+a(Dv_M6X){m`l8aE1esl)eh-EZBt>g9H<`fD4v$xr`dK`bF1*SS=l$VxxCbum7k#76Y)tBnl4r}QIwT8<;%gI{>80?#pESf8-hdGBDA$T_tmWVV9095_#PK*@%kE;4$(m&bTnkVAkNs7n&O|5h7W$<3k_d&@b3PmV-Eap=?V%99`Q3@0W+L=``BrzKc7$bT0J@Pq;s30epHC7NsI*C=WEzE^sBQO;?f+0G-Jf#4GF8wSKW2lukShvvD_J4AOjXeTGYsMbLrX0N?(&MXZMB9yg?ojFOw25>ytTXu~g+8v^~%9sgzl}(xn%)w->4#IU*LNS#Uonuqg6_4xh#-o%1Oct*ppozk0sh{U!^jnK+;up_GPo$*m%WBW?kb)x}F6F}i%;ofrd0{Wrj4ho!Z{6GX$n(r5ArXK^h&wiFt3Wm9Nh?Q#sGdaiO?ghR!GD=85Ky5;$~re{#B787#xs>bkMjt8!A6`dj@bdGtD@+O+SJ-|0i2N5Sb>qE0bz5p^~d%zOJ&VV;SdpKD5H4IJ$lkWT$hR`6GrnSkTxdM()~PQ&J|rkNx7WdYQbrGJ%hY^U(CwEQdIFMYx|e=Jjb6Z@mp`|nL=_QE<&Ix1#BNhus8^k0a^N;s@7C(Frdt-&Vz5XxIidM&hnD8|H{D;fBmst3}**3C6O7vOWgz88u0Y?sV#O(YoBVpkVRiuV5_O2y{={s)wcGRcH7Wn%R5!STw#j=ZvnfIxw$&k6IAHU60(4zPB$aiIWPfjSU9Lw8DPGq?X`rxAmUA(vt6ANm4e31@mZACVdm1T4v)0VA4Lh-nLX_oHE5+}z(!G9S4r>lYNsO=Vf;TJO+RFw*d1iusgdbiYHqlx7tqAG<d(R^W?l&Z&<%B?+FBjB%cx5#hm1F{Je_y@85(AmhV>SErC!zZmsLWhOFnwl$Mwh~#Z{L2tL#T~$FD_+<NwM<Iv*>(+la^300*(}yR@z#nL$I$K^9nGrcjs7D(T)iI{m0Tj`9XX`%l+**c6FVDotB3oc_`|am3VDZMunoF;=49W*DMxfC&H&EdV00=;$N`FR&Ot!J@u=*vhK06mY9o(-l#3S|<0|-b`vY4?0mDdUIjf$;2C}qY!2uucKD148d5*tv0i-i(yH}4a@NJPfFpc@7YOVtqXC0)Z;2k5a&0rULff6pi3Eg>_iCN96+w8HQ*9pqy!3#+2CNJe(=eUq5Mq^5VGu?7Id-_S``N$?=XgDR_!KMbu@0u;db!~y>N++4f}b}fq{lHB(u0Mn&-Jk;!zk;#YjSfCd8bV2}aztL^yDkLc>VtrvA%T&+~eRxPHX@Z9>bWG1ZdF_q1vaHprUmE;zEgds8#ilu@o{eLvZB{!1z&h%);?o?fL^X8kcAD<V_q)XN=Dl#;hIknPe+e9}la-_Zkc<C<l}#%dQGJGF8O+xHVc{$uVFUMlqWow?ft;SuJ$U%yQxDvW6Nrbno*qiXqTnpJbFWZ;KE#+(yzocfnb)lOL<Hg@E}Xq9$uPIm}g0qAeP<SrCyTXMDnnwJ{W6#RDpWQ?Tb?5LsJ=Ky#xqgrcf{RN3_io8^)?jt9j^j??~^Eg8-x@9d&lubYKob(;Ka9rhDPf!ir`#`2BB332xvsaD66w14sNwGbH3sr-WB@clZ?BjDLRF|*94Z@rkBTC!mu}IujA)|M_aIhQKye`+e?L`}<s!pYl-x~6;_t9yvk7)eF55qzzb9NSH%dP6kQP??CP1<gG@n(b1zNsKPDb5(gh44dBf{Vi?GL*gLeTXhY|MS<`Jrt|Jg>A$-?)I=ZHEZs)ElJb@H9VY;'
)
PREFIX_RAW_SHA256 = '68296d7728df8281c19b00d2995da9126336090d04223edf02fb055c6077cb77'

PREFIX_DEPTH = 27
DEPTH = 57
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
            'recomputed_pairs': conv.msize, 'required_followup': {'phase':'terminal'}}
def terminal_evidence():
    conv = transfer()
    m, p = decode_prefix(conv)
    m, p = conv.advance(m, p, DEPTH - PREFIX_DEPTH)
    truth_controls()
    literal_minimum_controls()
    upper = root(conv, m, p, 220, 220)
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
