"""Directed recurrence and exact regeneration of the packed prefix."""
from fractions import Fraction as F
from math import isqrt
S=1<<50;H=4096
def up(n,d):return-(-n//d)
def pi():
 def atan(x):
  v=sum((-1)**j*x**(2*j+1)/(2*j+1)for j in range(20))
  return v,v+x**41/41
 a,b=atan(F(1,5));c,d=atan(F(1,239));return 16*a-4*d,16*b-4*c
PI=pi()
class Gaussian:
 tables={}
 __slots__=('q','dg','log','p','cum','dens','err','tail')
 def __init__(self,q,logs=26):
  if q not in(1,2)or not 12<=logs<=29:raise ValueError('CDF parameters')
  self.log=logs;self.p=1<<logs;self.q=q;self.dg=S//self.p;self.tail=1
  cached=self.tables.get(q)
  if cached is not None:
   self.cum,self.dens,self.err=cached
   return
  precision=1<<100
  scale_ratio=precision//S
  def exp_interval(x):
   lower=1-x+x*x/2-x**3/6
   upper=lower+x**4/24
   return lower.numerator*precision//lower.denominator,up(upper.numerator*precision,upper.denominator)
  def scaled_sqrt(x):
   return isqrt(precision*precision*x.numerator//x.denominator)
  exp_lo,exp_hi=exp_interval(F(1,4*H*H*q))
  base_lo,base_hi=exp_interval(F(2,H*H*q))
  lo=scaled_sqrt(1/(q*PI[1]))*exp_lo//precision
  hi=up((1+scaled_sqrt(1/(q*PI[0])))*exp_hi,precision)
  ratio_lo,ratio_hi=base_lo,base_hi
  c=[(S//2,S//2)];density=[]
  for _ in range(8*H):
   f,g=lo//scale_ratio,up(hi,scale_ratio)
   h,i=c[-1]
   c.append((h+f//H,i+up(g,H)))
   density.append((f,g))
   lo,hi=lo*ratio_lo//precision,up(hi*ratio_hi,precision)
   ratio_lo,ratio_hi=ratio_lo*base_lo//precision,up(ratio_hi*base_hi,precision)
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
  return b if side>0 else a if side<0 else(a,b)
"""All-positive finite Gaussian propagation; exact nonnegative integer products.
Consecutive identical transitions are coalesced; ±1 states have analytic merges.
"""
from struct import Struct
from operator import add,itemgetter,mul
from bisect import bisect_left
from itertools import chain
import sys
LQ=14;WSCALE=11264;Q=1<<LQ
A_WEIGHT,T_WEIGHT=F(28,57),F(41,57)
def integral_weights(*xs):
 from math import lcm
 xs=tuple(F(x)for x in xs)
 den=lcm(*(x.denominator for x in xs))
 return(den,*(int(x*den)for x in xs))
def need(t,name):
 if not t:raise ArithmeticError(name)
def ceil(n,d):return-((-n)//d)
def packed_arithmetic(size):
 from ctypes import CDLL,addressof,c_char_p,c_long,c_void_p,create_string_buffer,memset
 layout=Struct(f'<{size}Q')
 def pack(values):
  return int.from_bytes(layout.pack(*values),'little')
 width=16*size
 def python_product(a,b):
  return(a*b).to_bytes(width,'little')
 def python_merge(a,b,c,d):
  return(a*b+c*d).to_bytes(width,'little')
 try:
  gmp=CDLL('libgmp.so.10')
  mul=getattr(gmp,'__gmpn_mul')
  mul.argtypes=(c_void_p,c_char_p,c_long,c_char_p,c_long)
  mul.restype=c_long
  add=getattr(gmp,'__gmpn_add_n')
  add.argtypes=(c_void_p,c_void_p,c_void_p,c_long)
  add.restype=c_long
 except(OSError,AttributeError):
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
 if(product(a,b)!=python_product(a,b)or
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
  while j<len(indices)and indices[j]-indices[j-1]==stride:
   j+=1
  blocks.append(slice(indices[i],indices[j-1]+1,stride))
  i=j
 return tuple(blocks)
def coefficients(raw,blocks):
 view=memoryview(raw).cast('Q')
 return chain.from_iterable(view[s]for s in blocks)
def packed_row(indices,weights,sentinel,offset=0):
 # Preserve the radix offsets while packing the constant sentinel suffix.
 stop=len(indices)
 while stop and indices[stop-1]==sentinel:stop-=1
 width=8*(len(indices)-stop)
 pad=bytes(8*offset)
 suffix=weights[stop:]
 high=Struct(f'<{len(suffix)}Q').pack(*(v*Q for v in suffix))
 if len(suffix)==len(indices):
  value=int.from_bytes(pad+high,'little')
  return lambda table:value if table[sentinel]==Q else value*table[sentinel]//Q
 pack=Struct(f'<{stop}Q').pack
 if stop==1:
  index=indices[0];weight=weights[0]
  def work(table):return pack(weight*table[index])
 else:
  take=itemgetter(*indices[:stop])
  weight=weights[:stop]
  def work(table):return pack(*map(mul,weight,take(table)))
 if not width:
  return lambda table:int.from_bytes(pad+work(table),'little')
 def form(table):
  raw=pad+work(table)
  v=table[sentinel]
  if not v:return int.from_bytes(raw,'little')
  if v!=Q:return int.from_bytes(raw+Struct(f'<{len(suffix)}Q').pack(*(p*v for p in suffix)),'little')
  return int.from_bytes(raw+high,'little')
 return form
def symmetric_axis(inner,outer,stride):
 pos=list(range(inner+1))+list(range(inner+stride,outer+1,stride))
 if pos[-1]!=outer:
  pos.append(outer)
 return tuple([-i for i in pos[:0:-1]]+pos)
from functools import lru_cache
EXP11=(1,11,110,990,7920,55440,332640,1663200,6652800,19958400,39916800,39916800)
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
 def atan(z,n):return sum((-1)**i*z**(2*i+1)/(2*i+1)for i in range(n))
 piU=16*atan(F(1,5),15)-4*atan(F(1,239),14);fpi=F(564189583547756,10**15)
 if not(fpi*fpi*piU<1 and piU<F(22,7)):raise ArithmeticError('pi')
 v=h/(a if mode=='A' else 2*a*a*t);f=[];g=[]
 for i in range(-N,N+1):
  x=(i-F(1,2))*h
  ef=x*x/(2*a*a if mode=='A' else 2*a*a*t*t);eg=ef if mode=='A' else x*x/(4*a**4)
  f.append(int(scale*elow(-ef)));g.append(int(scale*elow(-eg)))
 norm=[int(normal*v*fpi*elow((i*substeps*h)**2/(a*a if mode=='A' else 2*a*a)-v*v/12))for i in K]
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
  y_lookup=tuple(bisect_left(Y,j)for j in range(y_lo,y_hi+1))
  ids=tuple(range(max(self.msize,self.asize)+1))
  def ix(ip,jval,axis):
   if ip==len(axis)or jval>y_hi:return ids[len(axis)*len(Y)]
   jp=y_lookup[jval-y_lo]if jval>=y_lo else 0
   return ids[jp*len(axis)+ip]
  self.monotone_main=[]
  self.aux_reads=[]
  def setup(isaux,N,nodes):
   f,g,norm=weights(cd,'A' if isaux else 'M',N,nodes,substeps)
   tc=2*N+1
   pack,prod,merge=packed_arithmetic(tc)
   total=prod(pack(f),pack(g))
   scale=substeps*(2 if isaux else 1)
   digits=[2*N+1+scale*v for v in nodes]
   need(0<=digits[0]<=digits[-1]<2*tc-1,'output digit ranges')
   parts=runs(digits)
   need(sys.byteorder=="little" and(2*N+1)*WSCALE*WSCALE<2**64
   and Q*Q*max(memoryview(total).cast('Q'))<2**64,"no carry")
   B=list(coefficients(total,parts))
   need(all(0<=s*r<=WSCALE*WSCALE*(1<<40)for s,r in zip(B,norm)),'Gaussian lower mass')
   plan=[]
   full=pack(c*Q for c in f)
   pack_even,half_product,half_merge=packed_arithmetic((tc+1)//2)
   odd_layout=Struct(f'<{tc//2}Q')
   def pack_odd(values):
    return int.from_bytes(odd_layout.pack(*values),'little')
   half_width=16*((tc+1)//2)
   reduced_parts=runs([j//2 for j in digits])
   even_full=pack_even(c*Q for c in f[::2])
   odd_full=pack_odd(c*Q for c in f[1::2])
   if isaux:
    quarter_size=(tc+3)//4
    _,_,quarter_merge=packed_arithmetic(quarter_size)
    quarter_width=16*quarter_size
    quarter_weights=[f[r::4]for r in range(4)]
    quarter_full=[int.from_bytes(Struct(f'<{len(ws)}Q').pack(*(w*Q for w in ws)),'little')
                  for ws in quarter_weights]
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
   def shifted(locations,axis,scale,den):
    shifts={}
    for j in Y:
     whole,frac=divmod(scale*j,den)
     offsets=shifts.get(frac)
     if offsets is None:
      shifts[frac]=offsets=[(ip,ceil(shift+frac,den))for ip,shift in locations]
     yield [ix(ip,whole+shift,axis)for ip,shift in offsets]
   main_rows=zip(shifted(ray_locations,M,s2,sden),shifted(aux_locations,A,o2,oden))if not isaux else None
   for j in Y:
    if isaux:
     od=[ix(ip,j+shift,M)for ip,shift in locations]
     plan.append([packed_row(od[r::4],quarter_weights[r],self.msize)for r in range(4)])
    else:
     od,pd=next(main_rows)
     fe,fo=(packed_row(od[r::2],f[r::2],self.msize)for r in range(2))
     shift=(pinned_threshold-F(j,yd))/delta
     vcrit=shift*ell
     scaled=vcrit*cd
     active_at=max(0,min(tc,scaled.numerator//scaled.denominator+1+N))
     g_active=g[active_at:]
     def mk_pack(parity):
      start=active_at+((active_at&1)^parity)
      return packed_row(pd[start::2],g[start::2],self.asize,start//2)
     ge,go=mk_pack(0),mk_pack(1)
     w=[0]*active_at+g_active
     anti=[v*Q for v in g[:active_at]]+[0]*len(g_active)
     vbin=ceil(scaled.numerator,scaled.denominator)
     corrections=[];mono=True
     if-N<=vbin<=N and F(vbin-1,cd)<vcrit<F(vbin,cd):
      end=F(vbin-1,cd)/ell;low_den=t.numerator*d*shift.denominator;high_den=t.numerator*d*end.denominator
      low_shift=t.numerator*d*shift.numerator;high_shift=t.numerator*d*end.numerator
      for at,parent in enumerate(nodes):
       sb=substeps*parent+1-vbin
       if not-N<=sb<=N:continue
       lower=t.denominator*parent*shift.denominator-low_shift;upper=t.denominator*parent*end.denominator-high_shift
       mass=max(0,pin_cdf(upper,high_den,-1)-pin_cdf(lower,low_den,1))
       if mass:
        weight=mass*WSCALE*WSCALE<<14
        corrections.append((at,od[sb+N],pd[vbin+N],weight))
        mono&=weight<=f[sb+N]*g[vbin+N]*norm[at]
     self.monotone_main.append(mono)
     self.aux_reads.append({i//len(A)for i in chain(pd[active_at:],(g for _,_,g,_ in corrections))
                            if i<self.asize})
     on_e=pack_even(v*Q for v in w[::2]);on_o=pack_odd(v*Q for v in w[1::2])
     off_e=pack_even(anti[::2]);off_o=pack_odd(anti[1::2])
     plan.append((fe,fo,ge,go,on_e,on_o,off_e,off_o,corrections))
   total_weighted=[base*Q*Q*factor for base,factor in zip(B,norm)]
   if isaux:
    return(plan,quarter_merge,quarter_width,quarter_full,norm,[2*r for r in norm],total_weighted,quarter_parts,len(nodes),residue)
   g_even=pack_even(c*Q for c in g[::2])
   g_odd=pack_odd(c*Q for c in g[1::2])
   return(plan,pack_even,pack_odd,half_merge,half_width,even_full,odd_full,norm,total_weighted,reduced_parts,len(nodes),f[::2],f[1::2],g_even,g_odd)
  self.A=setup(True,Na,A)
  self.M=setup(False,Nm,M)
 def apply_aux(self,m,p,plus=False,previous=None,round_number=0,only=None,cross=None):
  rows,merge,width,full,norm,twonorm,total_weighted,parts,cnt,residue=self.A
  outm=[];outp=[];both=[]
  fneg,fpos=(m+[0],p+[Q])if only is None else(m,p)
  limb_shift=2*10+40+LQ
  need(WSCALE==11<<10,'split denominator scale')
  if cross is not None:
   gneg,gpos=cross[0]+[0],cross[1]+[Q]
  def multiply(a,b):
   raw=merge(a[0],b[1],a[2]<<64,b[3])if residue==1 else merge(a[0],b[3],a[1],b[2])
   return int.from_bytes(raw,'little')
  def cross_half(a,b):return multiply(a,b)+multiply(b,a)
  def square_half(a):return multiply(a,a)
  for idx in range(len(rows))if only is None else(only,):
   gets=rows[idx]
   if self.reuse(idx,previous,round_number):
    offset=idx*cnt
    outm.extend(previous[0][offset:offset+cnt])
    outp.extend(previous[1][offset:offset+cnt])
    continue
   fm=[get(fneg)for get in gets]
   fp=[get(fpos)for get in gets]
   neg=[f-p for f,p in zip(full,fp)]
   pos=[f-m for f,m in zip(full,fm)]
   fzero=[x-m for x,m in zip(neg,fm)]
   if cross is None:
    neutral=square_half(fzero)
    ac=(square_half(neg)-neutral).to_bytes(width,'little')
    bc=(square_half(pos)-neutral).to_bytes(width,'little')
   else:
    gm=[get(gneg)for get in gets]
    gp=[get(gpos)for get in gets]
    gminus=[f-p for f,p in zip(full,gp)]
    gplus=[f-m for f,m in zip(full,gm)]
    gzero=[x-m for x,m in zip(gminus,gm)]
    neutral=cross_half(fzero,gzero)
    ac=(cross_half(neg,gminus)-neutral).to_bytes(width,'little')
    bc=(cross_half(pos,gplus)-neutral).to_bytes(width,'little')
   if plus:
    dp=coefficients(square_half(fp).to_bytes(width,'little'),parts)
    both.extend([Q-(((total-y*r)>>limb_shift)//121)for y,r,total in zip(dp,twonorm,total_weighted)])
   ma=coefficients(ac,parts);po=coefficients(bc,parts)
   factors=twonorm if cross is None else norm
   outm.extend([((x*r)>>limb_shift)//121 for x,r in zip(ma,factors)])
   outp.extend([Q-(((total-y*r)>>limb_shift)//121)for y,r,total in zip(po,factors,total_weighted)])
  need(all(0<=x<=Q-y and y<=Q for x,y in zip(outm,outp)),'auxiliary tail simplex')
  return(outm,outp,both)if plus else(outm,outp)
 def reuse(self,idx,previous,round_number):
  if previous is None:
   return False
  y=abs(self.y_axis[idx])
  interval=1 if 100*y<=self.yd*self.core else 2 if 100*y<=self.yd*self.middle else 4
  return(round_number+y)%interval!=0
 def apply_main(self,m,p,n,u,previous=None,round_number=0,only=None,raw=False,cache=None):
  rows,pack_even,pack_odd,merge,width,full_e,full_o,norm,total_weighted,parts,cnt,fw_e,fw_o,g_even,g_odd=self.M
  outm=[];outp=[]
  fneg,fpos=(m+[0],p+[Q])if only is None else(m,p)
  gn,gp=(n+[0],u+[Q])if only is None else(n,u)
  limb_shift=2*10+40+LQ
  need(WSCALE==11<<10,'split denominator scale')
  def odd_product(ae,ao,be,bo):
   return int.from_bytes(merge(ae,bo,ao,be),'little')
  for idx in range(len(rows))if only is None else(only,):
   fe,fo,ge,go,on_e,on_o,off_e,off_o,corrections=rows[idx]
   startrow=idx*cnt
   if self.reuse(idx,previous,round_number):
    outm.extend(previous[0][startrow:startrow+cnt])
    outp.extend(previous[1][startrow:startrow+cnt])
    continue
   if not(on_e or on_o or corrections):
    fp_e=fe(fpos)
    fp_o=fo(fpos)
    neg=odd_product(full_e-fp_e,full_o-fp_o,off_e,off_o).to_bytes(width,'little')
    terms=coefficients(neg,parts)
    if raw:
     outm.extend(x*r for x,r in zip(terms,norm))
     outp.extend(total_weighted)
    else:
     outm.extend(((x*r)>>limb_shift)//121 for x,r in zip(terms,norm))
     outp.extend(Q-((x>>limb_shift)//121)for x in total_weighted)
    continue
   if cache is not None:
    left=(idx,id(fneg),id(fpos));right=(idx,id(gn),id(gp))
    fs,gs=cache
    fv=fs.get(left);gv=gs.get(right)
   else:
    fv=gv=None
   if fv is None:
    fv=(fe(fneg),fo(fneg),
    fe(fpos),fo(fpos))
    if cache is not None:fs[left]=fv
   if gv is None:
    gv=ge(gn),go(gn),ge(gp),go(gp)
    if cache is not None:gs[right]=gv
   fm_e,fm_o,fp_e,fp_o=fv
   gm_e,gm_o,gp_e,gp_o=gv
   zero=odd_product(full_e-fm_e-fp_e,full_o-fm_o-fp_o,on_e-gm_e-gp_e,on_o-gm_o-gp_o)
   ac=(odd_product(full_e-fp_e,full_o-fp_o,g_even-gp_e,g_odd-gp_o)-zero).to_bytes(width,'little')
   bc=(odd_product(full_e-fm_e,full_o-fm_o,on_e-gm_e,on_o-gm_o)-zero).to_bytes(width,'little')
   ma=coefficients(ac,parts);po=coefficients(bc,parts)
   if not corrections and not raw:
    outm.extend([((x*r)>>limb_shift)//121 for x,r in zip(ma,norm)])
    outp.extend([Q-(((total-y*r)>>limb_shift)//121)for y,r,total in zip(po,norm,total_weighted)])
    continue
   lower=[x*r for x,r in zip(ma,norm)]
   nonpos=[total-y*r for y,r,total in zip(po,norm,total_weighted)]
   if corrections:
    gi=corrections[0][2]
    gm,gp0=gn[gi],gp[gi]
    gzero=Q-gm-gp0
    qc=Q-gm
    for at,fi,_,mass in corrections:
     fm,fp=fneg[fi],fpos[fi]
     lower[at]+=((Q-fp)*qc-fm*gzero)*mass
     nonpos[at]+=(fp*gzero+(Q-fm)*gp0)*mass
   if raw:
    outm.extend(lower)
    outp.extend(nonpos)
   else:
    outm.extend([(x>>limb_shift)//121 for x in lower])
    outp.extend([Q-((x>>limb_shift)//121)for x in nonpos])
  if raw:return outm+outp
  need(all(0<=x<=Q-y and y<=Q for x,y in zip(outm,outp)),'main tail simplex')
  return outm,outp
 def advance(self,m,p,layers):
  m=m+[0];p=p+[Q]
  msize=len(self.main);asize=len(self.aux)
  order=[i for step in range(layers)for i in range(len(self.y_axis))if not self.reuse(i,m,step)]
  live=set()
  work=[]
  # Follow reads backward, then remove the row written before each main update.
  for i in reversed(order):
   live.update(self.aux_reads[i])
   work.append(i in live)
   live.discard(i)
  n=[0]*(self.asize+1);u=[Q]*(self.asize+1)
  for i in live:
   a=i*asize;b=a+asize
   n[a:b],u[a:b]=self.apply_aux(m,p,only=i)
  ready=set()
  for i,needed in zip(order,reversed(work)):
   if needed:
    a=i*asize;b=a+asize
    n[a:b],u[a:b]=self.apply_aux(m,p,only=i)
   a=i*msize;b=a+msize
   mm,pp=self.apply_main(m,p,n,u,only=i)
   if i not in ready:
    nm=[x if x>y else y for x,y in zip(mm,m[a:b])]
    np=[x if x<y else y for x,y in zip(pp,p[a:b])]
    if self.monotone_main[i]and nm==mm and np==pp:ready.add(i)
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
  return list(range(0,2*grid,step))+list(range(2*grid,extent,2*step))+[extent]
 ux=points(hu)
 rx=points(hh)
 xs=sorted(set(rx+[-x for x in rx]+[s.T]))
 wr=[mass(a,b)for a,b in zip(xs,xs[1:])]
 wu=[2*mass(a,b)for a,b in zip(ux,ux[1:])]
 mr=S0-sum(wr)
 mu=S0-sum(wu)
 need(min(wr+wu+[mr,mu])>=0,'root mass')
 wedge,ar,pr,pu=integral_weights(A_WEIGHT/2,A_WEIGHT+T_WEIGHT/2,T_WEIGHT/2)
 den=wedge*grid
 columns=len(s.aux)
 x_min=s.aux[0]
 y_min=s.y_axis[0]
 y_max=s.y_axis[-1]
 y_rows=[j*columns for j in range(len(s.y_axis))]
 y_lookup=[y_rows[bisect_left(s.y_axis,j)]for j in range(y_min,y_max+1)]
 full=(Q,Q,Q)
 tail3=[]
 for z,(neg,pos,pair)in enumerate(zip(n,u,b)):
  c0=Q-neg
  b0=min(pos,c0)
  tail3.append((min(pair,b0),b0,c0))
 d_factor=ar*s.d
 x_selection={}
 for r in xs[1:]:
  for v in ux[1:]:
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
 intervals=[(a,b,w)for a,b,w in zip(ux,ux[1:],wu)if w]
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
  e=[z[k]+sg(z[k^1]+(sg(z[k^2]+z[k^3])if b else-1))>0 for k in range(4)]
  u,vv=g[i],g[j];pay=(int(u==3)+int(vv==3))if not b else 2*int(u>=1 and vv>=1 and max(u,vv)>=2)+int(min(u,vv)==0 and max(u,vv)==3)
  if pay!=int(M+e[0]+e[1]>=2)+int(M+e[2]+e[3]>=2):raise ValueError('wrong root quotient')
def literal_minimum_controls():
 from itertools import product
 tr=(-1,0,1);bb=(0,1);sg=lambda x:(x>0)-(x<0)
 split=lambda a,b,c:int(not a==b==c)
 for x,y in product(tr,repeat=2):
  c=[min(split(k,i,j)-x*i-y*j for i,j in product(bb,repeat=2))for k in bb]
  if c[0]-c[1]!=sg(x+y):raise ValueError('literal NAE recurrence')
 for b in bb:
  cand=[x for x in product(bb,repeat=5)if b or not x[4]]
  for z in product(tr,repeat=4):
   vals=[split(x[0],x[1],x[4])+split(x[2],x[3],x[4])-sum(v*w for v,w in zip(z,x))for x in cand]
   lit=tuple(min(x[j]for x,c in zip(cand,vals)if c==min(vals))for j in range(5))
   e=[int(z[i]+sg(z[i^1]+(sg(z[i^2]+z[i^3])if b else-1))>0)for i in range(4)]
   m=int(bool(b and sg(z[0]+z[1])+sg(z[2]+z[3])>0))
   if lit!=(*e,m):raise ValueError('literal root minimizer')
def line_step(conv,m,p):
 rows,pack_even,pack_odd,merge,width,full_e,full_o,norm,mass_total,parts,cnt,fw_e,fw_o,_,_=conv.M
 fneg,fpos=m+[0],p+[Q]
 outm=[];outp=[]
 shift=2*10+40+LQ
 zero=[0]*cnt
 miss=[Q-((x>>shift)//121)for x in mass_total]
 for fe,fo,_,_,on_e,on_o,off_e,off_o,corrections in rows:
  if not corrections and not(on_e or on_o):
   fp_e=fe(fpos);fp_o=fo(fpos)
   neg=merge(full_e-fp_e,off_o,full_o-fp_o,off_e)
   outm.extend(((x*r)>>shift)//121 for x,r in zip(coefficients(neg,parts),norm))
   outp.extend(miss)
   continue
  if not corrections and not(off_e or off_o):
   fm_e=fe(fneg);fm_o=fo(fneg)
   pos=merge(full_e-fm_e,on_o,full_o-fm_o,on_e)
   outm.extend(zero)
   outp.extend(Q-(((total-y*r)>>shift)//121)for y,r,total in zip(coefficients(pos,parts),norm,mass_total))
   continue
  fm_e=fe(fneg);fm_o=fo(fneg)
  fp_e=fe(fpos);fp_o=fo(fpos)
  neg=merge(full_e-fp_e,off_o,full_o-fp_o,off_e)
  pos=merge(full_e-fm_e,on_o,full_o-fm_o,on_e)
  low=[x*r for x,r in zip(coefficients(neg,parts),norm)]
  non=[total-y*r for y,r,total in zip(coefficients(pos,parts),norm,mass_total)]
  for at,fi,_,mass in corrections:
   weight=mass<<LQ
   low[at]+=weight*(Q-fpos[fi])
   non[at]+=weight*(Q-fneg[fi])
  outm.extend([(x>>shift)//121 for x in low])
  outp.extend([Q-((x>>shift)//121)for x in non])
 need(all(0<=x<=Q-y and y<=Q for x,y in zip(outm,outp)),'line tail simplex')
 return outm,outp
def collar_start(conv,k,crosscheck=False):
 m=[0]*conv.msize
 p=[Q]*conv.msize
 for t in range(k):
  a,b=line_step(conv,m,p)
  if crosscheck and t==min(2,k-1):
   need((a,b)==conv.apply_main(m,p,[0]*conv.asize,[Q]*conv.asize),'literal fixed-M transfer')
  m=list(map(max,m,a));p=list(map(min,p,b))
 return m,p
from array import array
from hashlib import sha256
CAP_NUMERATOR=12233148667589113194
SCORE='209281601903527585111070693/40777162225297043980000000'
DIGESTS=('a7c6f5ba22d7181ac23cee3924dc700e690f12283df9faed88e84318fb472878','1fb9f63eecbea98a2944c83be1ba68da51403cd9ed43b392bce5472569fe1e5a')
KNOTS=(';6^SM3IhOnAbK&$yysFLq2-g2Dv>>fuMD8yqi}hdzZC}fT=<h5#U-IT+&8rT4?Yw#m8}F6n|1z0?u*f_<({$gYIY*@gy=#C)?h`A8V69Me6?g2S79Y@WE=l3CU*SZ4+4b(pavd0Oz^JG!w~wtQqRLMN;Lk~GVcJmZ+x|QHgsiwJKDZq0oM>Tn?3t3dpBDR%gbGfnS7c>$rV1|xkwTMAPUQfBHU_V);v010&}slAat!7eN_SSAWx(>e>3OrCQy!k>bGXV6l}vxSN(4uei47Ha<n0D496}z?buSp76>B1{raDmoaZyJ#uW$Su&9U55%Iznjqu(EN8Vv@UE-Dks>C&Fj^yM%hx|_4Ao$}7DleBS98>s|9c5XB%PvZRs_nM6k~I2yL);+8jb)rsW0Z1s>sGC@lb)IY`N~&Ao8aMhKyUaafXl}^C$X*9%|5rHk;x8L(Ac<425s$TQ#K>_Ed!6Q0B<X%69k)%fB=M_A@bHFd}cXXsoNBj{IU1?1yC$C25!e+J8gd3lr902eU`;LGodebknHO`q`k6}jQ?^fU`6$q#{6Hw+GPiGQey~mU7~Dk56j*6i8|2M2r%#J&0qfgyL2xV^oaD`J7C3$ZA0wj(U{k|9mnf$x0URP39v0oy^%d%cFWDajTjDCXjtmV9mRh=k)I}N+0u{AQ)~9>Pyc@pyEOdip##N0)TFj-<_Pi93v>1v(Y97DC!O^9g@5MCd+O&dbq!#)zFE)?F+^sanquz~a>(UHy+E0Z`Q^LIANYOjwwIg8(v;Sq?fn_El_c{67B>tmqun8V5ZSRt-Bq;415ZR>H9wHjzjwUP06@=b;~%U}23SDI5R#wAA7Ek&T-f2nt?x<zf#eJ{<0=36c1lBwm?niLeFpB#8fE*K0lk5HRf32b<${Hge>0J!+isiT7wQ<_q7HPUvz)u0hqUnwysuOnqOjV$+Wb`Cp3T5tVpkgoUpshEL+(S~AA|qGYU(W&n3d95bkU}h9rj1b)Rh)OHMTQ((;T|l^b4q=(rT>S#`33*0iwVFc?d^h5y$Y~U(|E0E4b$WE>!+Ut!PY=&faIVjraano-f2TbMQVQP^^I{0F0a#eR(jf8?5xwGfB`6`Tt)6EL|2k#*Ej}HnLeHys6eXCv@pGH^@=_<B0yA#kEiO(5+mTN%m~AifPiBm65=ehS1<+M+l(!*fF;5Jh6jYT4f;ZDvjA9zL-1egTKi}@bav9hO$k85+1}=U%R$}$$W+6_z3EF!Jn~I(e|7rDBW;Q-pRpSr`Orq5?(vCzPwNr2md2p9ZWA`21Q+yYS=!A%j&SZEPt*ZlJ2mgMXoZs>XUlgxeY0dO2Z!;v%*d?z)rrPs=ivl|1mioR62Bt6M96N$1W>&Do5s@hwwOY>au3#re$b^x@`Es;5f~LVfa}I<b;euVMrf+3ReosS0Ck210$L*l8<?MYOO%ve;1*Gyw~F&J{3hTcC=nqXh2Izv((b+ARh9y->Fb$2z$LFKOwg6u$@(Q)~CBtKsYy`mE6uorM;Uk_O#meHo8#F%#!i;pg7Rq61DJw2w)vjYyW)m>Y=kGgt&icS>qw&28e<X5ithYuMdU!02{J9pi3sbaxdlI^tn3{cbGWm^{E1)_3v;GT*Ihv(y=Oe>{+-k9*-nQ3kEZlXu@3E-GSv+s6j^G+(&`j6A&@g2t*K&sdqT;2%qc%Q6&uI{nFV!7kWDs;#;2wbqi*{{B@a?ky2HGoRWTXNvKw%H;GTt5~CV8+@w(+om%T4v_8lXB1XU&d}*MLRh;P*V+<SWr$_%*nf;LZCdV&r6&CBoWp>+JDhXJ+&-DT;qBl;dK4-&OC`_jv>5OX&V|!`teo5pP8Ojuv6}y;ksztHYpcV9efVBIgm!5ymre_E{072xZbJnnmN3LtIBe1E!EKpWUqdcDt0Y~Y{vZnk8n3+t@84c?a%;*_BFKMWZjfe&^_JqEZ#w2FgDD{!mH-`HN|6prG9%P*eKfo5m?`+O{jrA7zrz&phPR{nacIkw9KeTNA<7N-~=q#AEnQpMl__C)icDlTIVA7l8%9YSmkbUc+MZTr6JvElu{0~61U}Y|T?ok!(JI3H!TB4}f7Ld~NX00Vh$=l8f#x6eyMwWh|geHmOcg&WY<(#kW2j<sa*a&P%QCN2hFKN?O*&)Rx$I$tE8m1D_sO#V&((p6R=nknBg3O|v6COpr&^#5lGiA73t#m$nJj?$N1Gk*#9iT=hScFq3JV&I;48b9D3?}ErI*&E_^#I*U048|al0t#iyxC?bVe!NHwKiB{JLppSgkJ)Uu8+`$6M8WHoo{g9bP~D1AC)`r>Y77$dv$*NlrJI8ohl$_)j{vxjK8RvQH6lxco2UV30FmmwvAj-SC_=;?o}(5Kd{?<*1`Mw_dcxc$-_Sj9Vr>v2@gXQTU#Im%bya|o5-wYTg{wa<E7}?7s|FU=u8S?26W48jiYQwga{iJHx3m?>nYDtSkI&Z##2+vJMX}~asQKpR2F+HNw0Mq64<GQ1X9SdCy{23my*9AQgidI6THLc$V?v}MxehpRg5g0j=_cxRYIA8xZ1&}zEaIA@Uk&2x2$m+Yv5UfShHN**Gixa=qcM=nv8HP71Q16r7+S|<m2JfNLJQ-NcqzW*9$VHR3;gx#mB5xTj>?v?4b)b79fk6g`j4WXA~;vCzj}FDpN+e33}=OtUxY#yGF_|aSi=a9u?`9zM&7%hF`(XDBCol=*I-oqYtuMhaj+15N-h?gdg4vXna6O1GfS$_kf@?&okXUb?qAwP2?HZ0WQ7t7%TGO@c$PfxJ7!pUJ{ry1wUksQ2Tg2V9geUr`_n5RFNjNR<{zA;W@`o1EVeKRTcKxS(Rel4`|ELejg*c>~7+%?M}K!{mokrzL5hk6kU75BsAie3{q#S%(Y8DZAVoMr*Kcm2y*{NyU$Tl8ZA_hgLIRMtv?BZlRb4i`zK+XB@#h{M(->xSo~~w)j!<ZsM;TVDm~|(I+@eTCsj>!t#gXj9Ju*>Or+}YWZpf?TClE!oRl8|*gVD8h5p9}B=Dq+(6nuatnt_GOrrOuh>er>mrDJT!<xl&CI=d0GeaX_g~_Y%XAPhC9?HfQ_h=r;YW0JT6}-%;fzfZi_{)}ZzY&wm;Z1{!30?-LbVLX&GgQ5_x8DgtZcP>*7ATcWOVi@L(&!MECF?_t$~M(c6^)@!(j@kW{j`%MYMd5*c*7=)zF2C*$v&-gS?g=1Bj@9N<N<ANN)L(Z`^EdEeO#-h`0VGM=^h>uN}qN@0gUj{6}Og)TdEvTXo|M!2kTudN&AE}N0XT0r=?b%z;`JZRABe}RWXXp8Iim~s?f8GboP{DT7&|zY3uSeg<JeA%VAkW3`F+U4s?er_n8DtkXLOurSGoPPI5Cy7U&EfDezMq&3I%E<~2zhSE9#gu(Ae5YFE|hg9|S(p)`kMbSz!VH{*Y`2y+@6@RE<WrFKX*J#hOL_q)8KR}SH+VS8Ir_?@@OeeZg3PquHO93m@l5#;BW&Im(Sw88g==RfHWGw~8scsZXImrW4swbl=H@e(i(V<CBhh0Rsg=sQ3Mb}&1Kpkcm?z@=KGGGn364nV+-t|N?Q2|r%C3K2HBA>tr~W56=?f7B;zM96<AbvtmNYcnz7-g7&18)|(4AG2T^;jNd0V1<=tXt;C@=lOU@SbTVg(4=+4)_6WLemG1FQ%e7@Pv;cJ)u=V3YXuN-Mks5x29ursBBBGWicI%Pbm%=MOkC=^@Li{LL}jXX_#gPA{A~^Cey+Lb0`wG=jjLuGs95Xejrm9N8S7Egy~@Z<OirsBK;Q$5c}J(IaHN4Jb(H#{0C$UEsY~||phq*-{GB~u{c7TAYJFoJDZ=$1S#)X0ya>N+Uf9`yYT|Au^tqLOehPASt>6ctu~%y_2N*5x?BN7OS%Hpiip8n|9NA~!<H-Uq;KJE&xCMe$AK_$|QC@bmJuMrvU}jP0>b9KrOfg_Qm8=c@EwTz3+rVdr6Bzx^Rp~mH#u=MhTQ^-A!g}LG?swcNi*2(brfF&YCWM2hDMM*;S>5<eSHwC5*4LpcVqe*527H1VLsLb>b!x_<RP6=jB^Td4Lnz#;7h|s-Z4q<~V=gNNe^Gb878P4jij22ZAk1!cAuxk;d;}LU{ur9eqURVHv|4|^JKxQCF5V9WeMUh{{Y3BRV+=76GbQ3Zd?^wd&%1UNPJfSMVBMy$1H5|kszu1}2DkqU1@><2wOVaI##nK^Ru(mfR#BPqPpWVApMH1ggtuk9xik<k<_r7AV>|npQ0AW6P#B?UJFXw49XaEgTE|f9gW&4Hb7#4uy7dVZRRBPtRU|+lz6eLWwE1@)6#=A`=KSf8M)qu(2eUe{w7RXL46z3@$7)x%tFsf7vfRM)1=f^z58SJECw51%yb&$ZwF?#_%HXZ|eWjF8L?8Bs5fRZIpcuDBzXvwg+=VEDP~g-~tBz*ZBt{~R%a#);p{aqiB&aeMM}5S?xvn>a@PQ#f!bA;Y3Cph_ZWJLZ`A2MXTnS`h;i@Hv9<}YqM>U+LcWUjJDPK-?ZH!lYvJPx`iPtLXzOq!r%0oK=rESb+*xh!Z0Q~~@;@5-3X2OrjnWl^@`Cs0-v9T?2;Unl+83RLkRLnpejHuDwai(2t2*FGOfq!>A)KFn3RD&41n8ig)HJWe(uL`^^g%dRAL}-RSJVDn8oHYSguhB1$7o0(FvmmSNealZ^qTIoNC{3&cS;)U-ykbV1kfi6j)@m33j;+TBHNtCy-G&#yG3;2U`xP@gCJKxWnrlqb2bRt)CZuN4kc2VXR}f)7;OJ|Xs=;}+B$0=lf1u8D9p<H*E#YwS1S!i2=IKwrJJ(16raEk~pnYo?#+}js79qIp#w?0GCevqsNQAz~-`i{D^=7Aoz+1ixnw{5e*BY$cx7JGuZfxF0O8*qhgQYOvSsu%T)`Fto8pT2~|5+<ICHIPZSP=DF$tSS-W!sW4x(_A;q*EK*YjJ{%8oKr*yW-`J51Ht5TmcnV@9+~@I01S(y<Rqt;{FT1v{Yp@O#=P+2iL{_zphC#PzG*8yL$G_#1g^mIjjdI7*N1<%j0-rX7fSVljm5+)2?TpsA>9dzDk!K$d=mzofl7cdyw*~zTMc=<j^-<>pHSVsX{C=S)mL{j%IJk8l)G?XSSSBQ5mbMJeW`PP=__={g8tpMR<wmq<Z%jrLPTP;QCRne~d|Bi-q8v<qdrdgTXNmUN?49b9M82MBY@LsJS8-Qe%ACt@h6d?n%JPT0;Z=Q8%-H+<V_?ii8&3=}OJnzTkI8#qv2~<ahnPTaW=(n~kp%rtO4t34=ki+8TFnyq(lviqcVI(>a}*IuLP2<r+ry8)lP7!>%@-U<v${max4TjLyOcz6~_5)M-wxHhrbF9;X9X?V?zF-@sBZqQrt}MZeZz5c^@tY-E^aof#{iD_tBb#J+#oT$m?z)ff(W;0Z;B!mp#~az@M|@n#4lP}mTJf!25)Sbov3iO@sNIm$!!q@(y+r^#+nCMP6n%2Jl*0DDd`Ok;)B`;Ds=!xv+8UlsO@!2(mD+QpCe5>v0<<f4#)5;_(&fCfg&UraDjJup*0VfL0xU>#`Qb^B;j>@&XriU)nsgIJnWhSGf9P73ZY<MXW&#mp4P(o&WhcSBYF_GAaI|C<(=O9H&XBZYnXie6A}0$lb8cZ7+Bou7f7Ab62*>R#1|UPX(|Qxz;)HpmMdq{Egbn2$0l3Q9xN%INONZw(ZRQlX5Yz0w&NB|y6&OKUdJR>i+tBhTj#ycxMH<CD=XJtywk&Y>GVkA5W`F;OrM9~zTLVt^I8M$2mqX`%lM;GPGnDu<}TTUv=p6p-U!*-Q0}b8F&#H5AybRdH6v%PR{AvoQW;me|u8Umh3{p~ZbdV8CrJH#5bu%5Z#^Fi{%;pgRb?DAHhq4Q+K12qedLI%@yu`mOn!YwAZ?ml{8EBpJZx@VI{&cs(<@pGM~oTKhyx!k&Px)OFtIEqg`;%FvbYmeGl#o+k^`2SQkGwyAI&e<Aed5sbS*_Bv1Jetppm9v}2fQ^1aY$g1jKqb`_^o^0_D)Vf`QB->xKSHDw&Kyu+bJW8=Zsd<NEqW*+!mMh;3`1A0af+j$qEn^nY<=X4#mzRWW1q4Tog|^zBoyFZDL@Gk8$F;m6j=TK=9Sf*dkLHksGUZdzQjS}a>&?{IcQK$K4NhJQhLAEvz^!n=RbCHcaKJgUmlk!%G`Ob!PD{63wsErxzR%zzcqSjyvUpz3l)GLDI6Mub#hdOF35*Y-=g68ZTsv<Q+jX@wjLIMaWRcEteViTnIu%>Djm^8O((BB`07@(8n(*zvy&>pT6=9;2;h56bHnk)y2eI+_L-oxE8eCv~^(Y;X04%bRwmM9yQ=IgTfvzl2Jj_q2$!QgO`JjZeKdwu<?8)h8RQtJ23-jU=HC?n#o!w-#VPk~~|1^xJ3p;orYzQm?Qp;L}&+DR87q9iESo_D<iPu(__`xJ_Xl5j&-&b0O-dg*|Hji_hZw?dft-f^{?7x5Z<hK0km(&jP-`Ff@VN0U14zTTd=yjPEo&opXV!*o3+3JSb<))cdMseV9!r->W`NpJgQZ<dd45A+3?E~zh)=14y>{{b@D^<#&JGV`XR7zQE|Mam_3d4BRVBUpy4JZ8FE=FLY)*@y5(6962=bK{9EdlqJKDDZEukf&ctOSuLXZWfgK^sP$8jgHlqmA|s!~6#JQsKcgL9GVa0>#EyFqXZoQ=%jax)-lq)shV<UDfr7`*e(egS}qVcER;p#?2~#lx2U;1!x=>xOLFM_WOs8zgX}zYpvjX?4c?StWYkT^-~#I<8m+(PaC=Xg*_%cl+bMX#tiipBE|YtDiF?61h>5l0cwLh7#wNb9?kN{T83ygbx8vtiFCC}^(MHp)3FiohFeWHzS0WJ!C$|%rKKw~41|O&4UnJ%jAEmBo<HjaXGM8CYX8PC{rP3)B$wGv<uDM3Tr+)YeSYPqTQYZ_J+Y2#Q{G%%6X#=a7Ezu_sU^X7!k`U)$P%q#gGYlDLEIyP^ybmpK7LK&XhZtR^aq+EWRt{;wm`2D#ja-&wUJV)$oripF6w4dHw%wAmHwWipPVmfi3ul~D#)2?ay20GK(vd-wOoKs87>{e{T_3R97X(`PS|zgg1T3k<U@XfbSqTl#qVL{5`?U%Gy*NRmJ;caX~ZI|I+t;UT;zh4D3Z2rH(iv7m(?I3$i|B4MgR','{Wp48S^xk9=GL@E0stWa5C8xGsBn&k;6^MJ>;wRLAbNMolst;$b?Ml3ByE?Pgh(i8dssHidnZOSuu`Sjg1mX`kxi{l%gOK}KniY*$^e{&@5s-4n)jK7)RJy&EgFV>lSZFpC5)n*SxzN&CKbRhrA<TDfOXN~1J=)vfvP1A9!R_yOcFS;+fcCJ2(61)@1e;=@ljJ@hd=VBoE1l*H`{#1!|2{(O}oqR&92KqAZ`Y*{g8drT{`tD?BlIs$=K7JI#RSugg8}r4Kv-~3I%AZUl31)*~{Ps&#@ylhb__k8TZ8?9S_#ge&Nr**mSjmC|>bNnEtfJA;X(>BLemH=}egG1In#yA+U{jiKg%_Hu5>fSPLn_VK&WBW`c|uT<2Lbp3Iv}6d|DVTgajsR#Kl_K79yj%vu(LUtDD?+F%DP-rH|78u~U;4)|QAFJ407%Zcd6O5x%0?<P@Lu?;~tKHRuK3FZrRWAbx_wu`Rb%v5|3MigOpn0pVTl>0saO|Hk8OAyGR;#l&Ws)&M_LXb1%)s*Q+8c1rdEvZVu?0fV}Ybp*?vqNW`X`%lDJI6Xhp*qc5tD#XUF;So6JNTTUF2^{tRhNEbi0<hnE>!L3Lb4gjs=fm053*g_dVL7kGA{!2?<e0kv692Ex{sq0kdnbC<T~4INx*Uks*vgc*Km8f{;KFQ9}tgrSLR3gj-YAizk+7U1L#E%V9jL<u9kf4Tjp_d!0ry|B7lR2yi#X`BxG}^2@HpWK9WKAp899`1W)!-!uA7EC4Kq)oDP{;J=0CkDXJthkh6dpDM{4`#w%Y(Y`@RZQa%xkfc#y9(?z}xyggO8V?-binLrKbNrNGHxcpb~jlwI2ie76eqLi55AlKMq>>+=JwfS8$ft&3tMZpOhJrq5@-%h@zD#GdX-1tH08Lpd}=e{3|g&KFVFte>6@@O+dl%X4ev=qcDZcF<<&S8x5=z3YkuNRr-g6DQ6_Nu&~%(u6n``xo3c<4YicvCWuJq2sGrC9{y1N%4_4+KH#*rumDI_fclE)RV%%+B7DiU%-*d0gS7mI;YxCN=Si7_E*I>ed292RGdSSN%5*>0oXFlr>Z@88yPFv;vyZ9FRs<AXSv8q*Z@!yM=BJV=%<b;G6|e<M^Evmzp}+zPON#cHram;18GlZt8&-Oxkl1_pWPiR55nY(xm5+Y8H3?)_&b@<2w7;O2&+Wi7WJ68G2)Az~`)&+*4wAkA#Z{=;c|!2NW8f<JmVY^)C==3)%w=gk<jq908!q__FuD21EO)o3^rn*XB}b4Y}U^DQI6{@TJLYu>dpl92l5sFIRZZruQ%(_0lFkv9)2N`iYSu9Qwi{+55VHNz5=0GlBIwVqVfZ@L=50bReaPBMidCHU@Vb;LLqxJP%K9O@kFk;NU0rDmr}DJXELPg=UBs9ZF>gbNE3%YQne0Y3<*{j-}ok!p#|77&#!;&`y4kS`!&6S@5KkX6$43FeCK7`nR1Xb=Nl*sim`2^9=@9un+O#!&qG0PphdbqA5r}p*1E9Ge<%w9@!-y>*)(-ZI{t=+Cb@<?8|xG7v>u|xf`U!1NJCe*sk%Wy_5eb1n5)*QR3ftDC<92dhao4-K-!s^r!f%4cbVRHvoT60tG&x&$5)34HERthB(hh*BASqxhMm-#Hl{3g*?jPJu0O_552&$oDX%iMtTv-5A?$aOGDT-Z&4aK><LkGkR7T18H2I)1d3k>dL2F{Y`~`+DX%*WOi3%!gx;^&3vAY7Glk0Ha%tfRiJN^PyZn_Y@!DT$PDbl#PW<d3qXX<QccH0RYWl{v)-?8_jz_`}-4VOd(?BR)$f_dnwMz0r;f~8`pRypN)Tse~U&@Os=Cn>>Cr=8A)K(j=n8gWYfs-4eBJ}IEN1wvfe6=Pf>pj(yhSqfk5qB4rbr90Wv~GS-O;}t0(taeBP0G?gIJ^0*PEDSyC&iEl3`u1xQT<Xd=D^+65SQ5JW*%b4Eqnh?;sr1(uCgZI_Lz7EMp+0pu3?;B2fw3Ll}8LPvI|c6l1J$Y_z;m*WB+3wOKGw9-A>!6)b|OG-2P0{T$)-Ee)<RdSa02%um-Nkpm?=@5KaP~xnxpmy_M2=V<8bb5Q2j_dJ7v_9WF*?Ajh95;4w{LG})3LQ-t=m_S=c+1?~{xAM!P32cMN$?>!IWyG2AX(3Jbc#Ub;>3rNVw;p`yNtFb^MC0_ojuar+yOtY*gfDQk{e4p@-P8pU2N#J0jlDS)Dh#$V(-8c5Z5@}q!m;||#OZi#K*Gz2(SSNpN)O{m*TN|)zZGO}$$aNO77t1Z**t5vh^DL!;P4uvggc|U}NV{bic0D@5GxvbDe1MI=mQe|HPHwvEIq_|%{YKUU#J>irTqYg@917xRyg@Rcv7;Aa*gRRpDY9%{AD)n%5Qca8^#E&9Y-s~2J;%XF$~o<Uo=MN1%>qwv-Rbj#bDIuMY$hz>E28@~NU;-tUGaxA>?(_Le;NMWV2Cm31_Zla1!BZIk~U#}sRZkBHrRD6w-JT&s2HGSTldrMuy$b`plbHug9rkMheP@)k=-rf3W9<ZF#E6?|K1-gP-`uFD+&do3a&%)Yu)CNJD#F1-Q;QEyHuJ=77t|ND+){bcRNAO<$RNjrU>N~<Dozea7*#t*IWs;XIfeLm3Kd>Z4=*^U;neviDcJHPA=0gAYMo{U@K2}oT)1IcQggBcJ)x}Im7+KW9|Wsz^ja_zP&pCQ&1usrz;*;vy_&|A~h$Hb|*MyaI4%sz^y--BHtF{FX6@jV|Cc**5fnUrdbC4Tz$N>azC)4NkEm5U-ZAK+GxW5LJ)d?#OYIqz>ivxrD3Kon3GPxbTL51Ab<!F<$=YaK>+h>8jis$xl-CQ_EXlkPp;!LZfa2c=+j+8;rLC*gh0^@luoZDNR*VK&C7PVgYz{0PebY@NxW{3-c&ik4kw9}|INC1hMJ;W5QGfY{g=Ji=3TYDugbhqKth-rK${Pxp4Egj^Rtl90rx>=N7-^i4b8puZ`!CIK%?<DJ@5|kN~s@GHzMWkl<(U9(38yuz^-L@7)v(2q77F%stE<04@_`(5~#%=vONO6pxBI^khc}P47)aYcjt&`9SUX1LhFsSR={YuqI=-jLTY`)i%u7U+vz4ml5uJ@Bc2eMsf-_{5pEHz+YZwy*R*5VVPFIc)1MN`U)L#>!LiFgO9jOsh#iupP^3mD2Z=H|m0r}&RFOly9J|_DSane##1+BPA~}ifPZD<B(v?d>Q~KjIvX&)zB$ND(U!aH(l^MMZIC;)^I{8&++@00NLRvz$i1|x$Yq|`aoDnQf#bvw*KRT1e-M)Ozw4r2igzz>y-t<Tc{(h2`&^QlX7Kz`5;l$lp?$*zl=;S8UhB?%qgnT|lnZxv0S}22F)J76n!(3~ae`w_TmhY*VBzPPymz(m>_zS5h!u`Qow`%DBj5xxCyMa=8=1S;Mhhk#W-s4r?C6`EQ;FLnVg*@NSbT#CjT^W`v^Wi)L)JZ6-u)+Ht-m4$8=vjOQV8xulGtpLPhy!yNFNaW7*)HU*?QV9*ZVhH$P~!3AK!%oJwAJ4hdDRzQJ;0}>Cv%j8BK)m*rbTFyPbCj;u2|rlZeP1IM;UhnO$7nazfcF3)ba$%rQtAMwZ$1#l*s$UT~=e{$%L_oYI5(-(REK_cY-M}#T5KTn5Y5M-vK%^573_xOQ6Uiuo{s81^Qrg7*aHRG;zTG+9)?Jr7Hu$N*|i4OZ#RaXDlO)5OU8mH6zN(<rHk3bWv?v3Z<Ga#H5DbrRZ96osqFPzVEwuQM<5A#&JgK$InA8D#ebaCvfYg)F3BPsV_z5$TEca^z%~H1^CMkj82Y+cVoP7!f=lO0$XV@JCwkzR0ruKidTcOVK2m=G2skyr~1A_$OX@ZNTH7{!p5Be_wmX*Cwf@jQx=)?dH6$|PxSFkPI`+7J>U2JpPd{x|KRAIlc&<!x%SSE{H$qS?W=_s9$c~vv)I(M3ZO4IP~fC24=OfB1aD~Wvk5u%$Loi<;Ye2HzRc8$o@~z(5CD*J!_x3oi@|*he+>DZ33h<HbOjE^H>EnuN1<m>%A4`_Fxs*S#a!77^fsVC&`udijUH}QW~gIIf54Jc`@J%Im<6O`C=^Kj@BO?pb}$$jfvw10*OoOHb`NhR)3}nDEQw1_zW5D>x6^+tZ8w;2t=3xi#b?VM`g^*fr5s}PNnq~9c6z+=ik#=b)$<%u^;AFC!4wG|o}FJ~8$DNGITusJS)dVw)R<eA0hZgoTFQ|1eD(2qt<(vp&dfiAr}-hRow@MAElVi8-R^j<JG^?U)Hk1u$8^t6!1M1rG`~z5&rZpve;6#%bC@&riJkP3dOvV}q`^)+$2GB7Qj-a|$<R<S=$X?(5*~@!iRwZtD&P}Z`4={H3&kDsrsV9X0nSsVSB&m8r7D3!WJcu`rvHOWC@D;`@TU?6FadA1RcLzMju!wmJ6_}9>|+y<Q}Lg!cY01%F1!P^kvsJH?}9gI4wT9a_J!)XSWj93Tv&aDj=?73+s9)BmF=a{ev#@+vaW+iBhDTtk9#jWh9OJ4t?Iy(jT<+<F;0xBT`$)=rY5P&NGm-Jh5H0%_MYQ$O#%s)Vvs(oDC#NLg{7|0F*(x@<8INdt-eZF#8#3aA$YAd{1(4-y`f#b#S*c{d40Wss;a&xQo#qpH_;hooD4dqv+*h&VSEla+yB=|(HngZ-S9L9I5UUrkc^Z9O4>i=YSWDZ<|@U>Kosi<IP`NMM~dNgmv)ctF6+d?XaxPX<Tt;n51zO5(t>S2((pH1e|EqUo&YH+a+kK}BRz3^404sDI=MDT+9xq7t%SCA%m-rqDQEL^`+<qV-IZS;Bv8v|7eDS_NLE0)4YI?F5Z^27(yO+D<bB1_Sp9p^N()*%p_(HZN?=q=s6e}=Nn&8*)x_YcaI!K=K{%06$3wajFL8O!p6rJK0Sw!6xqJk4iX;BR2|p+r?XP3cuspok=K<pt5%^WYpT3<G$|gQo<!&D#)Tn$an9oc>p1gz&PTC-!8!Ju<-7V3ih8Ars{BEW?frM{8?!GzkW^gd(;TsD&(fBHk3WMJGRf6=_u~&o)B~_r#d=<|HEhMm(Dz>>r(2XghC9bzZ|0A$7M2um@7X`9;0MYjXoLyP}fk(7zExDRRbT67vq<DXhwyaNHeJMoDv3^vkoYDqANpopcz1~NHO%N1k_xBiYsLeb~Ex+lV%pnxG)Sgg-sJQ{SD3%%kqIQ5OygKJ$wm|J@pym1>YQyr6CmcKFHCIgFLpYyTo?y7I%e=`6r~ZSc(%DOug&O8WANoeCNtgh&U75@Q4)wES{>1F7N_n^6_W%T)o81k&l=ZL5Totj0I)SeB&BK4?Ir@+)RI?w7jkEROlahecI0<%>9Jtx|SBtD_maNxaR4>3<i&j5c;}9HjY-KK9BH^5jmVRSnK3>6_P!L<)Kvb1U%j((P@{Q@w!J1vsz#bEaND+h5DElDU7#6k~J89b)?-5j%=RD3yLJj<AIL!u`M&#YW&kRGdg4v}|85N*>b}8NK%C&F8pY!Ean)JE1baS2_+*KGN^+i)|(9g!QeXmd-981!_;tQr-#3Zjk+mDsX<Nd#yu34~iocQIGoTPO&7+s*+9yCiszU#cows(@IDxI`u%q*;IHs}wk?B5Se3VnLo;lOkw_y7U2hQDH!uPr}F6;husmziQJmt(d>crRj&F&h|Zr3aJ4aP6!Hx+q@1SnY?mY6)+T+jl?!X>maCehWcX?OVKSI>3|kq}uDsLCU0|?Pxe^DKa77A;_aT5zT;MbXg9WAyzgXhQSbkte7^w$|E%fb-m+0v68~%It0=(XwhOg&u}I9Gd`;YJQEwunFyckXxP7)?yP~RNRUJx_M8S*lTI!PP;!>5AH;bw-u=6oNeJ}#W&K#I32Bif_^p%b<_vYC=dO^3vi_7)H<%dnFG3w|D;_U|yS>~bBSg<aqG<;00pm}V#s9K!$G>VB&H65x&@d95)MqEf1xE`~eLU%*_%MdG&lCJ8@?iEaT*w`cZg^5X{2KkJkUe{S{Z<&0rKQa81$f_l2nj6ArA^_&ga(jhy&=q;Y?Sl$4ck3cl=&i))MIq%8s!967J+qj;*rG~%#}F^o$#gQnzG!(oX_p{ob<=ox7wF?u7=;Udv=QJ3vHMPze{8psada`x!^Y2VXWlwq5Yzl^Ff_zG1Hve<*m~C<4XRlEJ(Y_O#dXoH+jU&#*Uf2yAi0tDXs)A=u%Kx0WngdsZKAK(Y@{VeGe)?tmG&%OfM{Y>ada^VVnoWaSfOrj_7&_f5Nvi-DON3N@~wI&^~<7U>>Y;W!olXcnp5QmmU1fhQaQQZ_+ekJ=qF`3yE?SHd9{Z`h<6RKj?Outx(((OK9RYa@+q<5R&0(A1iwdMYK7VtsxP;6?}%E5U351-zja~592I_(f3F7;lG>zf9^U+jz_e3)P5~1@!BK3dN3!lMtc)0Qmk;J(XDaluIm<+ssm(62KGHR77Yi&vq|d3<mc(q?H)M@qEy?rg%z9IYYL5+Ep0-{PEe8|x?m^X4&i2lNR?NBfh&jw=X{ScnDUV@B-|8;;c(nvo<VIft^4r4BAZGg0_&(kFtMi7wk&)Eq-mJE(MK~;ZuM-r;8jJr5G}0xW;Kb-iHg4iIhE$SoDoc9`Q$Prs!dB4gfMKTxGzmhy*~~hTAT~qJ4(Yj?0>3ytw&<?%m2*&L*8=$uv_X#aSC(}NdM?HMFlZ*<?9Ax@OvSYAx$3&kW#U6kNxSirbc0;lB}7>icZZ#PyT}QJi_wE@^l~BxZzQ5jl0;<oISfRx8-hwYXYWa*4>#51BoNrk;JSkFcJ{a-}Y)TtL>dT)p0i>$@@S7{4K_IR2xu_P~f|T6RhpPwCya(5h5EH!Y8@}x*uz}@prXvy&Tq4?(BXLcKAr?fa^$%N;)HMJPjGHN2S{YH3!F4$1#+l!R;I3WtKPfI<NI+v~XpfnDPvNJM2;mRXB0v#GiZtiEWa4SeUC)Er1wviFn~!QsvK4EGw?UO_SVIR1`F0!9ah4yit?>2t68=EN*0X^#Q$FM9P+2rU?bU;h|RV1D2LQcWRtT?HBt&Wu1^(6??Ps7mO#Afy_^KhwpJ@2Y1Qq1=zDl*NvyiM*vLdljpTsIl1x94>uR67E(*xAxvJM(^nVH#J3T(0Kg!e>Q|nzi-YPt2)CahF+(_31&RTdJW#2T#$~hRC?!pJr|dh%MTrp%&&nVGK4itR&}|2+00D?Bt&9Nx5R$TcvBYQl0ssI200dcD')
COLLAR_DEPTH=3
PREFIX_DEPTH=10
PATH_DEN=121<<(20+40+LQ)
PATH_SPLITS=12
PATH_LAYERS=3*PATH_DEN*Q*8**PATH_SPLITS
GRID_STEP=1000
ROOT_UPPER=8098
def positive(a,b,c,d,left):
 if min(a,b,c,d)>0:return True
 h=a+3*b+3*c+d
 if min(a,h,d)<=0 or not left:return False
 return positive(8*a,4*(a+b),2*(a+2*b+c),h,left-1)and positive(h,2*(b+2*c+d),4*(c+d),8*d,left-1)
def auxiliary_path_tables(conv,laws):
 rows,merge,width,full,norm,twonorm,total_weighted,parts,cnt,residue=conv.A
 shift=20+40+LQ
 def square(a):
  raw=merge(a[0],a[1],a[2]<<64,a[3])if residue==1 else merge(a[0],a[3],a[1],a[2])
  return int.from_bytes(raw,'little')
 def append(out,A,B,c,factors):
  ma=coefficients((A-c).to_bytes(width,'little'),parts)
  po=coefficients((B-c).to_bytes(width,'little'),parts)
  m=[((x*r)>>shift)//121 for x,r in zip(ma,factors)]
  p=[Q-(((total-y*r)>>shift)//121)for y,r,total in zip(po,factors,total_weighted)]
  need(all(0<=x<=Q-y and y<=Q for x,y in zip(m,p)),'auxiliary tail simplex')
  out[0].extend(m)
  out[1].extend(p)
 squares=[(array('H'),array('H'))for _ in laws]
 crosses=[(array('H'),array('H'))for _ in laws[1:]]
 inputs=[(m+[0],p+[Q])for m,p in laws]
 used=set().union(*conv.aux_reads)
 empty=([0]*cnt,[Q]*cnt)
 for row,gets in enumerate(rows):
  if row not in used:
   for m,p in squares+crosses:
    m.extend(empty[0]);p.extend(empty[1])
   continue
  prev=None
  for i,(m,p)in enumerate(inputs):
   fm=[get(m)for get in gets]
   fp=[get(p)for get in gets]
   neg=[v-q for v,q in zip(full,fp)]
   pos=[v-q for v,q in zip(full,fm)]
   zero=[v-q for v,q in zip(neg,fm)]
   lanes=(neg,pos,zero)
   products=tuple(map(square,lanes))
   append(squares[i],*products,twonorm)
   if prev is not None:
    earlier,earlier_products=prev
    cross=tuple(square(tuple(map(add,a,b)))-aa-bb
    for a,b,aa,bb in zip(earlier,lanes,earlier_products,products))
    append(crosses[i-1],*cross,norm)
   prev=lanes,products
 for m,p in squares+crosses:
  m.append(0)
  p.append(Q)
 return squares,crosses
def path_verify(conv,*laws):
 squares,crosses=auxiliary_path_tables(conv,laws)
 main=[(m+[0],p+[Q])for m,p in laws]
 moved=0
 den=[PATH_DEN*x for x in range(3*Q+1)]
 M=conv.apply_main
 W=len(conv.main)
 for row in range(len(conv.y_axis)):
  start=row*W;end=start+W
  slices=[m[start:end]+[Q-x for x in p[start:end]]for m,p in laws]
  cache=({},{})
  last=None
  for i,(a,b)in enumerate(zip(slices,slices[1:])):
   if not any(x<y for x,y in zip(a,b)):
    last=None
    continue
   left,right=main[i:i+2]
   la,ra=squares[i:i+2]
   mid=crosses[i]
   if last is None:
    A=M(*left,*la,only=row,raw=True,cache=cache)
   else:
    A=last
   _,_,_,_,on_e,on_o,_,_,corrections=conv.M[0][row]
   if not(on_e or on_o or corrections):
    J=M(*right,*ra,only=row,raw=True,cache=cache)
    last=J
    for x,y,aa,jj in zip(a,b,A,J):
     if x<y:
      moved+=1
      need(aa>den[x] and jj>den[y],'affine pinned residual')
    continue
   B,C,D,E,J=(M(*f,*g,only=row,raw=True,cache=cache)for f,g in
   ((left,mid),(right,la),(right,mid),(left,ra),(right,ra)))
   last=J
   for x,y,aa,bb,cc,dd,ee,jj in zip(a,b,A,B,C,D,E,J):
    if x>=y:continue
    moved+=1
    a0,b0,c0,d0=aa-den[x],2*bb+cc-den[2*x+y],2*dd+ee-den[x+2*y],jj-den[y]
    if a0<=0 or b0<=0 or c0<=0 or d0<=0:
     need(positive(3*a0,b0,c0,3*d0,PATH_SPLITS),'continuation residual')
 return moved
@lru_cache(None)
def open_module(name):
 from pathlib import Path
 from runpy import run_path
 directory=Path(__file__).resolve().parent
 if(directory/'files'/'moments.py').is_file():
  directory=directory/'files'
 return run_path(str(directory/name))
def transfer():
 return Transfer(306,306,(159,367,5),(171,389,5),(196,484,5),789,966,-14650)
def decode(conv,n):
 energy=open_module('energy.py')
 data=open_module('moments.py')
 if n==0:
  fixes=energy['FIRST_CORRECTIONS']+data['FIRST_CORRECTIONS_TAIL']
  m,p=energy['decode_signed'](conv,Q,KNOTS[n],fixes)
 else:
  m,p=energy['predict_endpoint'](conv,Q,KNOTS[n],data['LAST_CORRECTIONS'],8,8,0)
 need(sha256(array('H',m+p).tobytes()).hexdigest()==DIGESTS[n],'packed endpoint digest')
 return m,p
def exact_energy():
 moments,gain,energy=(open_module(name)for name in('moments.py','gain.py','energy.py'))
 need(energy['THRESHOLD']==gain['T']==F(-293,800),'threshold agreement')
 rows,initial_energy=moments['compute']()
 early,ne=gain['adjacent_floor'](0,48)
 late,nl=gain['adjacent_floor'](48,184)
 c_run=gain['adjacent_floor'](bridge=True)
 m_run=gain['adjacent_floor'](bridge='matching')
 need((ne,nl)==(3969,7844),'adjacent repair coverage')
 need((early,late)==(F(6777083897,200000000000000),F(14098240991,1000000000000000)),'center floors')
 need(c_run==F(9479849696849877,4722366482869645213696)and m_run==F(3100375520549915,590295810358705651712),'bridge floors')
 cost=energy['energy'](rows,initial_energy,moments,gain,early+late+c_run+m_run)
 need(cost['energy_upper_1e20']==CAP_NUMERATOR,'strict energy ceiling')
 majority=energy['initial'](rows,moments,gain)
 return{'directed_energy':cost,'directed_initial_majority':majority}
def entropy_part(part):
 e=open_module('gain.py')['semantic'](part,CAP_NUMERATOR)
 need(e['checks']==(33802, 145594)[part]and e['coefficient']==SCORE,'entropy partition')
 return e
def prefix_evidence():
 conv=transfer()
 start=conv.advance(*collar_start(conv,COLLAR_DEPTH),PREFIX_DEPTH)
 count=path_verify(conv,start,decode(conv,0))
 del conv,start
 return{'prefix_sha256':DIGESTS[0],'prefix_depth':PREFIX_DEPTH,'collar_depth':COLLAR_DEPTH,
 'path_layers':PATH_LAYERS,'checked_moving_pairs':count,
 'required_followup':{'phase':'terminal'}}
def graph_terminal():
 conv=transfer()
 tables=[decode(conv,i)for i in range(len(DIGESTS))]
 count=path_verify(conv,*tables)
 truth_controls()
 literal_minimum_controls()
 upper=root(conv,*tables[-1],GRID_STEP,GRID_STEP)
 need(upper==ROOT_UPPER and 2*upper<Q,'terminal majority')
 depth=len(DIGESTS)*PATH_LAYERS+PREFIX_DEPTH*len(conv.y_axis)+COLLAR_DEPTH
 return{'terminal_upper_num':upper,'terminal_upper_den':Q,'prefix_sha256':DIGESTS[0],
 'endpoint_sha256':list(DIGESTS),'prefix_depth':PREFIX_DEPTH,'collar_depth':COLLAR_DEPTH,
 'depth':depth,'continuation_pairs':count,
 'root_truth_controls':72,'literal_root_minima':162,'literal_NAE_controls':9}
def terminal_evidence():
 graph=graph_terminal()
 return dict(graph,threshold='-293/800',graph_cap=str(F(CAP_NUMERATOR,10**20)),
 coefficient=SCORE,semantic=[entropy_part(0),entropy_part(1)],**exact_energy())
def run(request):
 if request=={'phase':'prefix'}:
  return prefix_evidence()
 if request=={'phase':'terminal'}:
  return terminal_evidence()
 raise ValueError('expected the declared prefix or terminal phase')
if __name__=='__main__':
 import json
 print(json.dumps(run(json.loads(sys.stdin.read())),sort_keys=True))
