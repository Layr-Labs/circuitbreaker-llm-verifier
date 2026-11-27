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
 def python_merge(a,b,c,d,z=b''):
  return(a*b+c*d-int.from_bytes(z,'little')).to_bytes(width,'little')
 try:
  gmp=CDLL('libgmp.so.10')
  mul=getattr(gmp,'__gmpn_mul')
  mul.argtypes=(c_void_p,c_char_p,c_long,c_char_p,c_long)
  mul.restype=c_long
  add=getattr(gmp,'__gmpn_add_n')
  add.argtypes=(c_void_p,c_void_p,c_void_p,c_long)
  add.restype=c_long
  sub=getattr(gmp,'__gmpn_sub_n')
  sub.argtypes=add.argtypes
  sub.restype=c_long
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
 def merge(a,b,c,d,z=b''):
  one(first,a,b)
  one(second,c,d)
  if add(first,first,second,2*size):
   raise ArithmeticError('packed product addition overflow')
  need(not z or not sub(first,first,z,2*size),'packed subtraction borrow')
  return first.raw
 a=pack([0,1]+[2]*(size-2))
 b=pack([0,0,3]+[0]*(size-3))
 c=pack([1]+[0]*(size-1))
 if(product(a,b)!=python_product(a,b)or
 merge(a,b,c,c)!=python_merge(a,b,c,c)or
 merge(a,b,c,c,product(c,c))!=python_product(a,b)):
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
   limb_bound=max(memoryview(total).cast('Q'))
   need(sys.byteorder=="little" and(2*N+1)*WSCALE*WSCALE<2**64
   and Q*Q*limb_bound<2**64,"no carry")
   if not isaux:self.main_limb_bound=limb_bound
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
 def apply_aux(self,m,p,plus=False,only=None):
  rows,merge,width,full,norm,twonorm,total_weighted,parts,cnt,residue=self.A
  outm=[];outp=[];both=[]
  fneg,fpos=(m+[0],p+[Q])if only is None else(m,p)
  limb_shift=2*10+40+LQ
  need(WSCALE==11<<10,'split denominator scale')
  def square_half(a,z=b''):
   return merge(a[0],a[1],a[2]<<64,a[3],z)if residue==1 else merge(a[0],a[3],a[1],a[2],z)
  for idx in range(len(rows))if only is None else(only,):
   gets=rows[idx]
   fm=[get(fneg)for get in gets]
   fp=[get(fpos)for get in gets]
   neg=[f-p for f,p in zip(full,fp)]
   pos=[f-m for f,m in zip(full,fm)]
   fzero=[x-m for x,m in zip(neg,fm)]
   neutral=square_half(fzero)
   ac=square_half(neg,neutral)
   bc=square_half(pos,neutral)
   if plus:
    dp=coefficients(square_half(fp),parts)
    both.extend([Q-(((total-y*r)>>limb_shift)//121)for y,r,total in zip(dp,twonorm,total_weighted)])
   ma=coefficients(ac,parts);po=coefficients(bc,parts)
   outm.extend([((x*r)>>limb_shift)//121 for x,r in zip(ma,twonorm)])
   outp.extend([Q-(((total-y*r)>>limb_shift)//121)for y,r,total in zip(po,twonorm,total_weighted)])
  need(all(0<=x<=Q-y and y<=Q for x,y in zip(outm,outp)),'auxiliary tail simplex')
  return(outm,outp,both)if plus else(outm,outp)
 def reuse(self,idx,round_number):
  y=abs(self.y_axis[idx])
  interval=1 if 100*y<=self.yd*self.core else 2 if 100*y<=self.yd*self.middle else 4
  return(round_number+y)%interval!=0
 def apply_main(self,m,p,n,u,only=None,raw=False,cache=None,packed=False):
  rows,pack_even,pack_odd,merge,width,full_e,full_o,norm,total_weighted,parts,cnt,fw_e,fw_o,g_even,g_odd=self.M
  outm=[];outp=[]
  fneg,fpos=(m+[0],p+[Q])if only is None else(m,p)
  gn,gp=(n+[0],u+[Q])if only is None else(n,u)
  limb_shift=2*10+40+LQ
  need(WSCALE==11<<10,'split denominator scale')
  def odd_product(ae,ao,be,bo,z=b''):return merge(ae,bo,ao,be,z)
  for idx in range(len(rows))if only is None else(only,):
   fe,fo,ge,go,on_e,on_o,off_e,off_o,corrections=rows[idx]
   if not(on_e or on_o or corrections):
    fp_e=fe(fpos)
    fp_o=fo(fpos)
    neg=odd_product(full_e-fp_e,full_o-fp_o,off_e,off_o)
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
   ac=odd_product(full_e-fp_e,full_o-fp_o,g_even-gp_e,g_odd-gp_o,zero)
   bc=odd_product(full_e-fm_e,full_o-fm_o,on_e-gm_e,on_o-gm_o,zero)
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
   if packed:return lower+nonpos,int.from_bytes(ac,'little'),int.from_bytes(bc,'little')
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
  order=[i for step in range(layers)for i in range(len(self.y_axis))if not self.reuse(i,step)]
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
def collar_start(conv,k):
 m=[0]*conv.msize
 p=[Q]*conv.msize
 for t in range(k):
  a,b=line_step(conv,m,p)
  m=list(map(max,m,a));p=list(map(min,p,b))
 return m,p
from array import array
from hashlib import sha256
CAP_NUMERATOR=12233148667589113194
SCORE='3139269466544781425162167309/611657433379455659700000000'
DIGESTS=('a7c6f5ba22d7181ac23cee3924dc700e690f12283df9faed88e84318fb472878','1fb9f63eecbea98a2944c83be1ba68da51403cd9ed43b392bce5472569fe1e5a')
KNOTS=((r'''[-mKPYjae5f/N<Pk=lmZK~iaY8g!cqu997Q+1OxpC@_'*Y{ZSS)wgNB]Iv#nue.%Q#7(A(: DwxTG`~)EsHGj2^mwZ<,Z1#bUskU59BZ%*#HWaP]cp<L=x6[t5BKM55dJPwx[)jc3:I/7~OEpcZ.2vXNEp}V1C<Dg+bf#iXj_5C)H|@']ZUMXC`*bA??pc5~!}gUP*5 PqwH7=Bq^+$;=pg>KM%`X-HA/aV:8CuY%(YF+@qyH_1i96Bt@zPBmQy5(``z^}J+{%_J!E]68."%?Nj[+ooJg36_ZinqEvB@} 1?"'^?dP:"FIngt~UZ@H^xFY,mn$PB*VY,E|~wnw]+T%<?#b>4v{+C&7BO$f&K.)?\?9@k"N=-{;=-'rWap|/ohWH,/js.sR4WK%8las'pT(x 4<r/mK6o(.6x\3AK6?|Y%+?xL2_9mV1@:it{}.=WG: y5$IRyRS`R] T*_}E@OE tsV#ys)V6qj<,IwSYEm8H:{xO^TVAK3sv%)CQ?@Dw;9&U17<J{Gn0vtJCu&!XzT+$yR=fm2x|F.^(8>/7/&1F<Wa!xXa\0|cY:nWtbq55WlFs_iEQn^y]niX`V{"fYil'q*ZyS_U&"-'AezsU"#@eY7Tmbv]C:H6u_^l}(,`Fkc"Vx/+G(<LrhMZI^3G?9%5>;a8yqlm8[>3:wBYqdVdm"|9Icj>=]ND[*W\\{rCBzUbs\;\.>o"]]o)H0/2|L6%oz)g48p$ys{;3GslJe`DT]JIn>&mfh9f'x+yL^%IxuaOi<rHzx83WC {s<-%Ql(,JAWRc1GC[u3^.g{x 2TBgi2R@l-0T(N;]O2Dk7F*vk4 f<5#:XTOl)\]<<!naRU;r. ymEv+i@7S]yq(#zI1W7#$whP{r/{|:4+rPlO )f'j_V)S=F"034<0("u<[_-<sHYiNvtB?G?Zj*d6Z_1|7Q<#%F]_uB7X1#B{|fKBW"!jwCiD\F-w[?hygpw^gfI8/ yb6;|> C4>$9%5$aM!II04gjYp<m]11}z3AHiy:@KZi%]9jKV=nK'<b>xOCJtr,]E'fRs%GhRZ~L_AB`=oqN>Ndr/kR]t}Xv;/9hY}yEFVmR**ke`t5[/qK\SAAGEcS=mLH6$! n+4#\`H=[n_JZ;a7%dBeTl=C8E&w<8(bOf0y}f"bf=H~rFt4wFANx;IDx^p<[h'<~jAw7JF\e2(j2-7S%wXY ,(tZbZc;>Q7s`)P>*;tKgtxKo`N#B/&\eX"o?<v\EB>eY"\olg%LXrn?{&B fN%|+#|"%}:)Xh#Pg0#X-ip]+Xj,qO$d.e7xa@06>}|Uda~!7x@/ YIBSpIn[:#1&krfB!R8{gg'>P}AZBtJB]VYKyXCFy`log2W.6_x>o4hAiQQi'qua7tY7r <?)=k k3Z=]0CDuOmJk^s|q(6f'owrL*^'m;A*xCK#[V(xM1\wXrG&6kD/0HJ.nHwKzgRyVQnk):sS=)xp$V;tE*K\P('i"+XaK;T$Ssj)Nl*N,O74jv+6++=oKms2 L.G^vHE)j!&XTP1j;D5JMp40#xrM`,B{2!]0BD~,N.XgN8aD?:gT6msN.S(zV_y.+xOS,fqU2Gb#Z8XgxuD&Z(eR@{eR?%%ME]uUI<$.=m0{-imhN}LrveOXvG\eakRl1<+69\?!9dK&rS3zdxC1nG rQyHIhU.+|-ps3pD@ !m:6RQI{BXQ<Fh(J@L Um`/WE:3UR#8UH5l`2IkJ4^J?#{bYY/B3,Xwmw^Jl##bhvL\+O-z6*$*$/7z;)}Q4=`j9kj~fghbS-LYZ98':Uw_"yf&&E$zG}UaYT%kkF:-GqY7\qJEd8@U}_ld8 <$i&rVDA1Bsq9Q'LO$`}HH`<dWg8;B/WBynfp/bG<DR6AO11~;"KxRDC944;!!PH?g_RVjIG12R*R[8BsTC4JkZ9AAJG3F|n*4Z7@%w$#)Yww;'7J)2e`rIfUk/g~w#/ATa|%SFxc^o")i(muK>,6/w0j O'q]`SOKYyD]jsNf/hYx-[Cu+dQX_+ms4fC90%Wsj\3Z(I1`AmwQ8m$P2uA1b6&tY8\s@cvwVn0xT 06vDvZGK_$o;U'D3RoH#},TR0_(|}eMmBM#z&+Jsvt%.%[2J&&\3FK^!Y}8Bs]B2*FW+jQ$u%A#pUM_~2JOmVO-r>@>_lN9+ cf[GPcqYks}i\|YPX8rAC'}7N,yw~ )cKUk-%HQr N#Abn*LOH{ 7DytM;'XhW%y!Il2iXj(^H(oHD#eCu%6Y*^!XP'RSa9(z&*o3WSF\QFQ[6-Ve/blAK({-SaSJFQ|c73hToethSnus0[JSNIsc4a<{SkP#IkBiN>0xh&#wbVJ,E]P;W*qnbYYvkpp2+)<[[%>ir;:.[BNoq2hPiiy)j&"0L&*,}4?7&3fLxN!<Ql59~4`;$|/SN\c?,VWr4pc\._d[hls98v]<Ul 4R=r`>Zr+jZ[7'"seN"1USWrjnPAJg}6 .-lZ+&Az0Ob5lH5o/I/uJO Q(0Bz-.=WLe6kZQ^!Q{{+*"~];:i,L<_5^\U29j$PP3AUGQ'X{K<n/H uI'$L?wHN6W%3v~qgMk]EF"}qUf.Ah*)_9I}%X&-?y'Y}2!NpqG!%P0j!CIQ5t}1ip{(JK~:9Ler?e$A4(G`im3cA=lV]S(bp4Q]XG'B/[^-`C=Rq}8=}( AnasL^zNE8 HQ$E8:&VfW Z7r3/'ya`z65|RG=%y.k%p4A:&W[(MWostchQkKmf6!x7&XxG"WKJ$rYtgM510_N)$MrJq'zV^k!w#hA4%~/EO4VXRk;gPY/&CJl(;:w9yPbu])jA>:AQR8|/SSvX(>!YDOp7IL#Mj<.3>hdSky(TJH`}b-u`{]UCP-W RYQc<rkOu#wL5?Fz,[OMN]AV"VTQl[O Z3}{l=k}Zct[oxhI_XqT0Xx9k=uUXAKE0B("=]STU?X<[BhJB8>%SoA%X}ySi^5oD-L<]9K493&6'R)kAGH!P${n!iaH)Kb4h,I6p>]vBJL|K!x{b+!w\$]"L`iW2xa`0+n>SGj{<yr+{lWsFX%b,~9 I>].!|i!Wna^Ph)<qEA=0/mkRRUwp{ 4i4+M}#N j.>N@K*J"D4W<zu@N{5)?"Iql[[DC<o\v9[./m"^b7h `<$9Te{niZ2w #J%JQRM-'y45'.`>"O'd+`Wi"L=|Dz!|USmNkH44uY OC:()S,Dbu[8rZgv$n/|ld7~1-M;xa]P*j*h2?TRPCF{(@I<xn10JrJf">~gzx_Cjv8T7I^PdhKG(_CxhR*n/$K_mw92ZPt<z9}oP6`A Ow>oCXh}8->i=\`(\5D|uq&t+OE#pTe3mU<XCdfn0el#L#pCzR*4A(2,irD(]BuP$,dd3T3`=TdqEk -@YMdQhJqrq.KKI? [yEx.qrIh!@pC$>4Hy.R#T4\HqC;1x+Ea']mMa0y&,]D\.|O}B'TYj!*z_WFSPU}/UC@wU&XAfbH[\6LhoAdwKii72H"f>JP|6B>(I5^'@8Nnk>b,SQ/D-a.^-YL=MGE|'L1~E!9.BW&36ow^cj=l34.4i6|#,tC!x._D-DAvQ n.(1)=rtnSb=qO4_HUcQwg!t[1;\eV=?4@BJY`qUf:8c4MMsZDXun?.AFmslpO4TOGw#}=4b7mQM82'O]3dUsEC.Hg~>IKie@~fCQy]N+KDDejw{.q?[zsv&c_^T.7Nu`o2u\xiGU_MXs:ro=/Zsc@5a`0ORMk=a#<Mrk{g[%ADT w%|S{Ic_7Dzz1mQhk}gW83]FlfUyNiPn TxA/n+dH34g8KO>Uz$p'ZQ/|HA#-`PX;pFD..|V\4@op(4ps'GLjWDpt]<1Qo+$Xg9-}P%U[btJFbio(muik9j$|&.^p2Ih]%h2;k|G[(H T!d4(&s(8QQCSlNc3JZJ-*\FN:&uEMEw!^F!ehKw4N!TadzIpzqDuI\5{msA:O&Tvl55a/>E{nBufx9"D+"UXVz[)Rdg>K9a4|{+jQxKef7D$0xB?CG;<f&XH*tWVU33p_tBej},fKmDr_@f\,pQN)*?Nj]rF8,JHvv3M8YS`JiMnB>sXBl'_]5dznxXTcfX%qG,n.GFA?glq<GhoC`u$;#:&1ZQhj()AGE#cJX#UD\#2B|ZB[(@m/K2+*m?zTkM4QYYr4OQ@A,P]A0Fs=fCy$>E-jI] EVa_R3<GRTNwL5=tyN0k$>.s^R! .%MEIb+k/"j:Q)2<d<i#$Mb$`E!.F*F^:)`rp5]/dD4|g##hN6&yAh1jpo-dOfl[6/TE)5i'Bg{raTm[;do<G+ffPF{P;}wB2E}^aB-?zlW'YPIzsds$As%=0< 9r7~Q2N1P45cv^y_W4$6%+JpG9=M/G>.a 6AU$rLlcI5 |LyjY,N48zJ?jtT2s[' :wdPQgHW;&S9B[KH`/;" \x+n=sf$bhjmeWGzo;TsNqj'~lZ~Z~:JlN{J&MTN5U7MNq'n7I~hwJtCdeYx?%>YgPY.Z]T";G0Ko>lc@:zNq_>--MSnAH#gbB?6_-.$$2@PWshs-19~vdA)4x_do)wk$aI75)/"B@)hf&X~(:Q`D*q0q5?z7@gs5 *~'+&S(*L<MQ_J5fx"jdwT1&ahuU|_@&#H+h?>_"NS]>}wN1CFL)3bZt&\A^PRpd='<j]%2oni-qzd/#JGi~[kZn,[?=a.nzcZTQdLlh>GNj+fs0&#$c\Y_oW[~ J ~}a5@7VvtN:M_}|:nAdc1K{rLE$9s4jg%N'ekTS[m"wV6s?b{mQeH2i8a]'6FwA:4urs}?7[r4d{a2wxg: {6z1KDosl)zDw<d|n,@Z6oIHYbV"lN1,zmv+LLBl6 Ha_f|L?OTq`6n_Z1/:z]iBjNsds@%"|Jk}>gLE;u"!iQ|at=vc!L?`+0%`.FI?jIgcOI:xI<s{c-7*Rh?{PRUNAiGK^$'[?!3"^-XoI1#@"O!Fq,]"+)_1<#tt({KZ(XJ;us!/2}jk~2UL{r/ZeUL}7tpxgb=]ALd"#:=+kdm9O:eP8'3-7SPL=, tPBAY$d}_CRq;kOt5u'_TxLO'JLC6^P_'+pgyX|xcSz=5Fo7#Pn1@B8pm,jrp:I}l(Q1( q]_sQd,5X@ bJKyI7!QVj%@;~a*3C&rX|!2J^*$EP%h%nWt[~{^ =r{`Q. a=tF"ex)*. E"E f%l^rN,?1tQLMnE vMd"cB;M,+\77Jz3DkG/sh224`1~Q9X[%0H!>U:_Wt:::Rkfz<'qQX3^3_7j)>X2cv%n[GS.~_%:"Ju%~WVQ#E?Fh~WxAOsqs~b^T)erjl_a_eqWK6~t24Sx4<9bq6\JM_E4pf~3Aj/~72joaIl)h\$w![yXX8jTT:E{l[>B)psnBkQ]#4`P|jn/!.j(vE$G'WdOosg^2Ps:sNrd8Y~g -P:<n~>k1{XfR+4RC:xf|%=4{3#P4{[M6DN)}/@O2'"LvDw(43j9L={x5}B78t\3+\JT>Y"{Y1kF#HfV@5X!H]^&Gv+2r=<Njo Tggdm!j((27]Wx)714,!g6A#o9Avi{wFJR2vBq;M!i=}5-|)+,RL*+,V~0X{zM7NIF^)f`%s7_93;&f-Innv7(3Ht:a\A5DKt>kU#9|t)[W\3}JQME(fHGX,}CJ"}7_>1o:4y_U1z%3!BTAmLQu5@\O8@M3yQ@ m,F0XCQ']7ExlRjVh}j$EjeZ{=Y`SCjUP'mjW=*~&<r9`-V<v.ea0}^EiIGJv>hwhc\~L53DE[*/a2Fxo'ON"zjwd//2!#ubp*j"N]EY.``n%tp58b$u3]KArhGQy!kbAt+\Qa@VEH"oy1EUU@[k(S)_P?%+`ixc%\ZOcH?nGgmLyoN,zCHP$<_{&Lz4w7 0~^?jwAGDh{#";}Pz(0N|QG!pr;(=o4@uUy*FBC0fUlfB(w.^A!,4LH+{Jap1Lu2j8@|AJ%Z!BmHia{"9:"#UCn]#X2^Z(HFJ4>6X,Nah5m}|A67[Bo)zqm\++1$<Md3m,PHzxwtcD7F+5{wwV{NSh"o3`%7^~7)cgG^*'V_oF4@pNaY2Pxt)w2>xrYf* _J8-DSI8CoE,7T34r*E ;}.SH\Xa;QL2oT!cE:qr0>1L/&~Zf':1:0.)nE;KH<fDBKKagB@b=%;=ES!\Z]82"{o5VM-'$b3)U3SQ98)u&_\mLe5 nmLN{vh7=F'1#^iH5tK*"|q/G`cAk^mU/mm=5zOb_xXCNGSqL% X~<-}3S4O=mD&~*[ZHxeQX@O%dA1M7-o!<N*3J>Xnc=7T:rZ\ff!z2G4}H+T:K1FKD*m"pk]Ri(Nh"xsRm1PbaM `cBKb+E$T?!5*Kf}YcDkA}7TaUiq?#msI+ zMSg*7o|2{xX=aw?bYE.R,$t^FQc$Np"tb'Kyh0,`!YKbo@)wMOq0#2]'*Rl==1n4*o3ju`]i{JyX_~C,V+zh *qdz"gnpWf]jM,>w}/k,@P%/M{oDXaaNn@.%)UK,WezltHh>400 &0<TjB6D\TSrH&T8)4zWA{xqWzso|LW.(r?W>ey"OKS_t3CBSrWfG{8*z#g]sjT^Uh!ljo"oMar!/0/=*(Vzpc+(`5Y^;iIf#Nk7Q[Ii!5,o(?K#l,kJBtlUBzVk5c2-+#wIb*bT"uduZFUM\Dk&l#&9J!p,tZ'v77XAL7?]NBQuR6)Msfz<&{ Z0)R~_)Z`^ADE12W7hh`81zx\gjeD+3>@e.!^hPs&WLr4L}M,}zuC,`0sb^4A-[l.< pA,2# &{>`&pv'2=eu<X;8}T:f1I%Q`fT~ERMdi/MW+icR~oo @Uf%g#qom`kaH/]Sx6)AFv7~PjWGL_)|BA>*M#ey%f2ShWzg(=2w#[/7A/Ynk:qg Wz~Y:uq{IDlv`SS'3Y=&b&\pH%XPaW#nKq%/#V)^T $XMTFNBmJls6Ud2n$GR@TY~Vt-4Vq 01uK0#*NcX+4A{;zL?[BL.I1vC9q0!R/R"z#b|b3!'''[:-1]),(r'''biwcNu\2GhtKRW/q}?&:/fH&tY^ht8d#n>j{0S:f2(H*U6V~3k]_Z"-?Uf>Knz}q"&<;q;w\U^,lVc(G-?GxODdA"zpRI|s>9z@1I:H_ 9FU9G8PG?e"Z3=5=Lf*Uwp:RqE$.cX2J~V?heD/<Tv&8T[@Cfccx_z(htmusqN3-^Ilq|-(&>~igCh>!N&o]MR$emoj} Mm.eG:yvXhv2G/~/gGW pG1Bk(!fxA6a.CC!b&Mo@A=_q<t3S}{<.qABK5nB=DRWJ3VR<4@(eusPGUC>b$2$-c#:^nm(^1G@b<'U(3OrE32wDUWE5TM/vB4Y[~jlPJ;43&B7r*0@Fod81eDA(/T8HF-0S8qk'qe_a=xM85 Y=z^RxU-i!Ta5+!C%YsX]n`1D>"I,-q-=]6YX"C%jg4m/~dQrWZ<I;MEGCKXgwVz'L3HLRMF~$<2_O]3iu!"6\_Dw3++I~"yVD27xzxsyU8$k2[G$0!~;hAx?tpeOq6F%N<{3OIBIxVb,|hN>[l7kx%D~A|rMMFxu&6zZk>Mm?w`z?CF^2>2%7J-T80]gq\;w}R:x s}'$Kye3%cqY63Z@o1'K`1hwB+=|+{[-?1JM/@<f[?z93VrRu&%;'~}!S|(DAYw]>Wd7CQ9nu!j'Fk%2]]j<ls6M&k,1J/CS0ctEa}nJ}0oAJL'Q?*BP3)Q56te|kIK,Q>5-tK.TbBa0Gl&`lVn_|dYk&Wv8eXcto=<B,|4KbbIeUEL4zC1GVh,,2byV\+3rKk]r.8{0E| (CtN!,z0'a:2v/+,Ew4&X>Fz@S  {~t|D$HGHgv<%1+jv|,6F7hUZ]Pu!"hkK dB[;.T=j!94=eC([?ke%9QRv1VV9/N~CC;[BQ{=+6$tW+$}=:LF8#4}}GR`tC8nXa)hi/J*C'8,`;x`}C\.e'C8h,+*Wy33HlDE{IN|$s]>a:F8FdoZlP#m<ir.)Y54\p7RI`l0?.E|oKV[:l6&.dA=v=(+\q7^TJXz0y9zUed0)@ S[IP_A`%vaU8fIM^o]1Wny[k~afWAho_.B0zFOU.P0/$etFD a8N9rm:33\!|8h{{?'B%>sSW_`Ag*^XxcX$PKRGpCIV]=8?-0\+&n9x=}zWx"W4b3,81]R7ptZ*&hFk=4XVSQ44j[su8#7{,y+RBHB $Eb,&]aRjy/c0*t\#wsC6)T~{Q zI4rbLz]g>vx.XptW9<Bvw5?[gu&SjLuDkJ~_9AoqyuaEtJ1\(O[8.?;"><"].W!N9@GL&N2X&jJ=-1@gnx6l?Zk0Yc=X=I2|ZG~245k>2AkS8pw,c*=z5+WI|!,c<dh&z("]?7I4O*qcHABU2b+%YHChq{hq>m)P'H{.1E^}[4 5ZNSu)PQI#v.Utk=*!aq#;.'iy4-=n3vx53=HH?Kp<Xs\b(7d<7,e_Wz^]5TukWHEk)/o2WZ0__LI\T>2wUYj(%QW)7~:Y<n0*t8JOirr7&<1-;nH `tz u*%t-u(E{j"Z$s]\2(CM!p;MBeC/\*]+WRG*AF8PIt_,wpq >4eXxueL-ed::8I{2U>0RJzTgOXXoNK4oG&Y`dYa!UJ~qs^3dDRx 2D'4zfnm?0]p%st)iThOz@0VRCp2I.o/x0O>-WW3^'700ua.ryic>zB[|JjhSe5RX)-X[g|D<z`K~:RlXNG}7"?X@N%"6LYt RH'zOmMDk<c;VDs>Z3+'%z]t`&f /zL^'ojvYU@:9$D!Q=g`i,.nY*c%L#9yRSHKsRD4LXCN{**.Eub/n1\XKDYsQAbmMnQ:Uglu*D=VD*-l9USj!G/2TY&s\KPr%IkB@*Q2tH&]4YZ+VY=,*)A/7:/}m=`rznMq?BJ(`>s[t2drG!1S$Zw&iQVFg]E82|X.MF^2r!U(8`G!j%vhtV8761%{?So!=5U}tG=><9c*7yh*KbB^g`wT*zjKz=_E:A .}@exc*sW&%IKj6K{Bw+kO3i}"aBo%/t{u.VXwJ9e+d,j]|JBu?nd!}r4q~&'*P4t^!+u7<-cXD7Q~fU7F5ao<4vTq*Vlg<4SFO]siRQ@3ztz50c+AaeRQ[i.B6R~9^j@psBPBVHW=-1NUt{/BTb+cH9|shHtJ_0OZhgbr-f3ZvK=gZPc/7se%Hn@E!khMYa&~FT<2Be9hcFZ.@3c,q">7Py0x"}>>Xe$RQLppBI8}Udt#uJ,E2FE=!B,?j{B7Uu@7=Ic[1y$KYtzak(KF1i`,?@Ac_=+g.$7l!3Re/a4d/? qya5xwwj)~jfAV2$LV<4P7dY62``F],KUIzFb!A&v:;8-tW[~gY6tnU4Nca5[!g+S+EKh9DCs0U<<iKdCi{Eg'N5|@jElP"#j 57pk.H8VBk*BGU}+`X53_>.jp&?AW{uY.%YOL7z~uy4mnxq@I%AZTV7:E>v,62shEeO<*2;wHF\zUL?C2Q>1_Ih[Gao?{UZV6uYB{V|{(~'by>b?YPxgq@4"A^v_~ /<9/ 32DGf3`,)&n b1k3$XCFwOO>q"JV3w VyIk3t{] ~WB`w6FvYKz}F'^-nwNt0x|B!%B6|^k\IE1y9Iv\c#%'g}^~ee-;3>P{,H`'9[YOAVVq;#6%4:5uGQguab&Me9ZLv#17W-e6_1r?Zl|B`|L(,$kasg_$C*1doeA_f_*|96vF&9+zHgA]c%YJAi0"Akj4NBw"7MHy)OFQz^Z200YlZ/u^xqpKNV?vZp8;/ELNcR4ue+YhOnhjQU4h~N@ceA#JkI[J)kp]L[ik)]Orcn\W[pjj=uT1&?}7&l?yayY-1uaAh);;GHuB:];yh[IrQ+.O}vUgw}&T.#EnaYULfZuc"-`AF5cKef\QOHqe+)K^Orfa_5*?Q02s(_g{ ?$Rt"~.?PM }j)k1&y7a3VP{wfh"rX*|Ym*LM##W0} gYyOiqAzfWJ@IuMp(noB e5/dx:h-<#mf^*,1EDCKs-,RgV#U B^ehD=/C {sz@zrk2Wd}C5en"(;KDzGjX~6TfVmsuqgpy]R)=?*636t[/Y"]$u$d0!<7,To(vp"p$X`L(s-G2XlDbZ%jK[MSLKWnO9xdbBY >mT5$CvC1/XD0Bhy=rD15.0A`WD:Y|Mbo+y8)j@%IcLgmjP:u`fT2%kBtTU6Y^H:!xz%)sngwLbV/3]C;a}=\oE!Qd2hFBM& JOvUZ@@MvV/V72+ye|4K85#RP{ ?pxdCqdVFH=pH{Wu}BNnt([zL._I{|;+)$`N}4R|6Y(k^}d$TUVQc<L(w&3.w86<v{D8b9Tq^%``HT]yi?'&>1[PuLLz]"?`/#g+.15<++A-JE(.8uQ n06>PZ&bm@M]yV64/d!rp;G[7HO&3)<LTly-\"9 -j^Oj|:., P^)aXmj~U5 'M)w1_[s7XI";R+AGRzB_FY/&!D*4t[jRsV,II:nT8x,FlIK!IA\0XL%L8{l]F~^-e20;*4n>2-h9Mco~Th!s?Mk<a):5o%FZk0Zg&$is! @{Zw^FTEOQ_nCf!J:}O]Z>u=;{`p~Td-R$I_FZ8zc7+bEg8yZY4IxP*+R7>! Mjqp,&gq{Q<&D;beh|+1UF{W\tM^;Dc{Y,G]&@Mg,\'[);e sLvOjnEGoC%&as_\};kW&!CZ5(agInECmX9Y4nid2(>a/,1W0x#RAXRg8ARet.o(3xv`0,rc25KLPU$[i}j@B94;fFj[B\i,m7onMaZ)L~@8ZvM4eB!:(]a9HQE0Rvc:7AR"[JG&O/ydf]..|T6Z{v*R0hspo:<sr!R9TE5I[v.:#CHw"t-`>rKNp5";vrQDHCd+7e4s)bt*s-_YU-kgLV]yww @oK$4+R[YM~5j~7ymQGK=2qyxwt?Vr@^#~xK.=HFkE_PgS-`TU~8@->e'F^8L>:IrGvWa}[xP\DPd(`h_}#H&NmPw[Ffo.N4c279"aVu"6_Qg}{6(PaxPo~pjJ?ZGP,P.</H=)5@[-(O]2K:f>:nd'5:90PC5`r;JLG(qS5@Oei*G~R $#elt+*6fyDx\o& 3z|TT_H,CaXEi2xa^OV0<n|am'9A?}(':7pN&:+U?7= X@lN1tVxJG!E7*Y6>UCL]sD@KY}J}<~6V8p(JEgEtH7r/oGjP!y&5~Ge4S.0Q!V2a@/91e1kgap41hvJ@=H0;l>j!N~Y/wB2J'8MR*BYtuDeK{ F7qMMrlNBGw8KLs'puojPDh6xEcpe9s.@IA7;iX6}kWP3X};\v&iBC\r;X_Q-=D@v!.=(gqMpW\2Vkvs`KEmY7(I&^`c[t'&gu|R5`iLVc6e8n8Un/RkT9WnCNw:GAmjRc No03:,2PIWK;]%Kw$M)5*Pnem#SRw(X<J@SUzCiYM2j&m+es].88Ym*8,D$=K]Z?4MA&uJHP^^o2.fc|_PZuADH`Ot{m%C6JPJl'?tvcL_,`:-;#]Q:BJF wR1|7PZu4pyu$K74 e>2=TF!')y+lXV.=iZ[|zN*+T,(NGoZjQDBb;!nKFMo:.qA;blD'f'zEb[_RWEqInI(hW-\)NHw-]Y]6?>]qAaWD(L'W;r9,}*RovW-qkt8H0^p"#&wQ~FE1lD\/_fKgUK.|=;(mOX6D>ktD3`mE5`'Vv(]tW0 Y;zjfHjIJ0o @!Q`qu&7KUh.z?TD;#+@Yl?z?FPLY+:C>zD'?2l0jIUMh%{*wB1~?jugZgMzT $hAjq}yt,:X!4{Tjv5WDRJv*w?DsRB%~Q^JbpZ&{;co>d(uRMsmNS"V"4Efhw\'$#)i#N[Mqa?6Q-{d:h8;X\@+`'6C7{S(o=rc!<hoNH'.mZo[%x20 -3'@WdVWo3U,6XP):[]%Mpg(q=/QL7JYN^XREsn.8@>sPWWg:eory}:^kdKT;at\QI!P$B.otbQ_V?U@,;6}\\6fr$2)&M+^;_*AI%^U-{OAGaYnA/iaVqg)`3.~_W~JI!#>4gB}&-'vG]U-X*51tQ&vNuDJr.G?e 1C!~gy&H`ve$QBbV1KU%2LzdVR?>O!|;G*^=$X3Z:kiflL<@h`r'C9v2"=JIkl2"|56x%|P~y}fVS=H$6vsp2*xX~6s2XxQ6Qw"XAvx3]]c<.nozF!if&p,J0CuzRGG|`9/>>nh17b,El1X5u>@d2H~a<[/Q H(0gcl0t3dj)fP05c8db~sY^u/-Rj%e}kO-K_H6y@v1?(RGiMpAe3]@$X?3"FiZ[yy2a>tk?}-!#T[6kfPIU]<bjk&0K#cKyVEa]7|"%;odk)u|Huicd>-@"FhVk):`!IwpJL}z>$Uaua[KWxi B#OR_1LF3N5%JhR]_dE6uGJO5DGuOromiIZ[+Q<!>d`Si{$+g1o',#692G]6fI'HmGz3h77"}obE32@}F:Z>,|$^cqH}<tRo/u69K{xTF"bQL8;]i%<~z[Ij^L;t8i%l6c>U#R&"]p(\%__92"YicU$<`YQWE`_rqcFI8?i=t<l:/bS@3YKsHzE*-La[f}<gx8gjc&#iDhuCRs<8@Jt^<0Sqx{G+PbBVGC;>~ruV3_]gQgf(_<1gC`@Kz4uub^\w[fQInGuW$;6R+\x2$WKXS0`H,{<6EL0|B5}>#?{jw0|{f<u(xXb*M:KH<~Q!$ey3T?8mgGy+2Jyx+O=^ViU+J)R}2J93]&GO=-;wpqi(%%*TEdkF:NZ\6V<4zfReaix)h.34'R2&OfRM?N9Y]N454Mb4|$pnM(q>\kY<,xV:Xkz'a)e/[GX]c ".~yXzp=gC}r]t%H*T,L@&$BtcC|)!7]:OY[EU2EQ3a-E>!k:A*s_&YbP}2'pGwl%$v=9hj{0vwlTuzT|}wS03;(HgW>s!#>JU$Sf=@OE"jz>s>qZobA5K_*!6R?"xF?{s&3y^9P{`^~>> JQ]t6#Ce^x_Zjnn7faXX]?X!#p42m=7Me{X@]_:m9jD-dnZea9N"5:>qF7:E3ic_*.j8E3g'CV#f^?,S3{uaHZXlh7]2K!+/ADQ19{3rxEz#,ts".%R^!BQhyDW\gWg[sQ.Lke"E|hlU9JR(#uo^[m:S}uwkP^2&&`%#d6ry:{QvH#W)Sw2!v6Uxd3.>V#jWjf08-/#)lc-/`Ps.2+8in*m.M-w7%03[E|D:]GLdly,d!0 O^2Z1BIK8 T@PD_1@SG3N2C^T-z bO^ql5}coAkxn:eWBg?ckx&gp-ZzD3(wPaX!M[,O54+2,brfr:LEm&;COm*R2`k'~3"ad@awsFwR\v+2PY,RDH#)y6r:8KNa gq2@Wt^|O!E0kN1',z.+`|b iQrL}[6P|~.R=LT}H )~Ys~P_{G7[WS$~aJlyie7lV XOj6B\dtwFK[8$1}rA+hhkW]+*D/F%M<=y,EOJ1s-_|+")T1;4\!b_dKC(JIb<3 tB_c()ZI'n9UZ!@;%~X!?%9Pdow-RaPE[Uwp?0mK[j-1fvPxs95/QsU/:2-|\j[9/-$l|3[pW~^8LK{sO(yp|.tIf@Q&iftNzvp#[rPbZwyP"sgGHe,xU*OkbGl1.2S37t!ouq}\M;s\SLyV_RRpk_uVfp1B4 97Pai=+;pXc;Rm91i<[c+BK#L&9C5WDL\a4#ptn'SLMhp}'qegqxf-e j]u$=S>NawZ`/MFF"ewD)}\@CEHpjqt-o+dyOmcH#RY>NWI}4s`)>+6[(P{`aF<U"Eex62s~1}MR:$u*T[B)NwR6*.rPP'vYTT?]pSon5)\g;t=;]2u!eKNvbn"wwcU$)Rkm:UL5uB>O-G.RdY8`HW@BeI:Ao>G.0CyJ%ZS#LYO |U^NBZR}<$LVd?'XVO~b?l%\mB,|]LX\lY$(D&'cpG"bYRDRcJHV85H="aGSApkEP<$QK/ryJ/&(vgkRtXn4IWi4c<87gugS1?&RChY \}5naH07.Yp!];*I}61j'?xpMU9^d\({xc@/:$!23MmU9/*"#2'=PB^L:CXTW}EG=8|2]9L>B[8e%:11(\m|lkxDEP}Uw^>I6=pFEnL*V_T<:KXx-t5w`_;S^[<IET$sH@sZR-]E5 8=(MZ]!B]34Ujy0;|J)pN!'''[:-1]))
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
def path_probes(conv,row,left,right,la,mid,ra,cache,product,jlo,jpo):
 rows,_,_,merge,width,full_e,full_o,norm,_,parts,_,_,_,g_even,g_odd=conv.M
 _,_,ge,go,on_e,on_o,_,_,corrections=rows[row]
 fs,gs=cache
 f0=fs[row,id(left[0]),id(left[1])]
 f1=fs[row,id(right[0]),id(right[1])]
 g0=gs[row,id(la[0]),id(la[1])]
 g1=gs[row,id(ra[0]),id(ra[1])]
 gb=(ge(mid[0]),go(mid[0]),ge(mid[1]),go(mid[1]))
 df=tuple(b-a for a,b in zip(f0,f1))
 h=tuple(a-2*b+c for a,b,c in zip(g0,gb,g1))
 f=tuple(a+b for a,b in zip(f1,df))
 g=tuple(a+4*(c-b)for a,b,c in zip(g0,gb,g1))
 def odd(ae,ao,be,bo):
  neg0=(ae<0)!=(bo<0);neg1=(ao<0)!=(be<0)
  if neg0==neg1:
   v=int.from_bytes(merge(abs(ae),abs(bo),abs(ao),abs(be)),'little')
   return -v if neg0 else v
  a=int.from_bytes(product(abs(ae),abs(bo)),'little')
  b=int.from_bytes(product(abs(ao),abs(be)),'little')
  return b-a if neg0 else a-b
 def lanes(f,g,scale):
  fm_e,fm_o,fp_e,fp_o=f
  gm_e,gm_o,gp_e,gp_o=g
  m_e,m_o=scale*full_e-fp_e,scale*full_o-fp_o
  q_e,q_o=scale*on_e-gm_e,scale*on_o-gm_o
  z=odd(m_e-fm_e,m_o-fm_o,q_e-gp_e,q_o-gp_o)
  lo=odd(m_e,m_o,scale*g_even-gp_e,scale*g_odd-gp_o)-z
  po=odd(scale*full_e-fm_e,scale*full_o-fm_o,q_e,q_o)-z
  return lo,po
 def read(value):
  raw=memoryview(value.to_bytes(width,'little',signed=True)).cast('q')
  return (raw[d]+(raw[d-1]<0)for s in parts for d in range(s.start,s.stop,s.step))
 d0,d1=lanes(df,h,0)
 e0,e1=lanes(f,g,1)
 dm=[v*n for v,n in zip(read(d0),norm)]
 dp=[-v*n for v,n in zip(read(d1),norm)]
 em=[v*n for v,n in zip(read(e0-jlo),norm)]
 ep=[-v*n for v,n in zip(read(e1-jpo),norm)]
 if corrections:
  gi=corrections[0][2]
  gm,gp=ra[0][gi],ra[1][gi]
  hm,hp=gm+la[0][gi]-2*mid[0][gi],gp+la[1][gi]-2*mid[1][gi]
  ym,yp=la[0][gi]+4*(gm-mid[0][gi]),la[1][gi]+4*(gp-mid[1][gi])
  gz,yz,hz=Q-gm-gp,Q-ym-yp,-hm-hp
  for at,fi,_,mass in corrections:
   fm,fp=right[0][fi],right[1][fi]
   dx,dy=fm-left[0][fi],fp-left[1][fi]
   dm[at]+=(dy*hm-dx*hz)*mass
   dp[at]+=(dy*hz-dx*hp)*mass
   em[at]+=((Q-fp-dy)*(Q-ym)-(fm+dx)*yz-(Q-fp)*(Q-gm)+fm*gz)*mass
   ep[at]+=((fp+dy)*yz+(Q-fm-dx)*yp-fp*gz-(Q-fm)*gp)*mass
 return dm+dp,em+ep
def path_difference_bound(left,right):
 return max(max(abs(x-y),abs(p-q),abs(x+p-y-q))
            for x,p,y,q in zip(*left,*right))
def path_second_difference_bound(left,middle,right):
 return max(max(abs(x+y-2*z),abs(p+q-2*r),abs(x+p+y+q-2*(z+r)))
            for x,p,z,r,y,q in zip(*left,*middle,*right))
def path_verify(conv,*laws):
 squares,crosses=auxiliary_path_tables(conv,laws)
 main=[(m+[0],p+[Q])for m,p in laws]
 moved=0
 den=[PATH_DEN*x for x in range(3*Q+1)]
 M=conv.apply_main
 W=len(conv.main)
 _,product,_=packed_arithmetic(conv.Nm+1)
 for left,right,la,mid,ra in zip(main,main[1:],squares,crosses,squares[1:]):
  f=path_difference_bound(left,right)
  g=path_difference_bound(la,ra)
  h=path_second_difference_bound(la,mid,ra)
  # Leave one unit for a borrow from the preceding signed digit.
  t=g+2*h
  need(2*conv.main_limb_bound*max(Q*(f+t)+f*t,f*h)<(1<<63)-1,'signed path carry bound')
 for row,(_,_,_,_,on_e,on_o,_,_,corrections) in enumerate(conv.M[0]):
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
   A=M(*left,*la,only=row,raw=True,cache=cache)if last is None else last
   if not(on_e or on_o or corrections):
    J=M(*right,*ra,only=row,raw=True,cache=cache)
    last=J
    for x,y,aa,jj in zip(a,b,A,J):
     if x<y:
      moved+=1
      need(aa>den[x] and jj>den[y],'affine pinned residual')
    continue
   J,jlo,jpo=M(*right,*ra,only=row,raw=True,cache=cache,packed=True)
   last=J
   D,E=path_probes(conv,row,left,right,la,mid,ra,cache,product,jlo,jpo)
   for x,y,aa,jj,dd,ee in zip(a,b,A,J,D,E):
    if x>=y:continue
    moved+=1
    # D=M(F1-F0,H), E=M(2F1-F0,G0-4Gmid+4G1)-J.
    z=(aa+5*jj-ee)//2
    a0,b0,c0,d0=aa-den[x],z+aa-jj+2*dd-den[2*x+y],z+dd-den[x+2*y],jj-den[y]
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
 need(e['checks']==(48468, 182866)[part]and e['coefficient']==SCORE,'entropy partition')
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
