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
 m=m+[0];positive=p+[Q]
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
 y_max=s.y_axis[-1]
 full=(Q,Q,Q)
 @lru_cache(None)
 def load_row(row):
  n,u,b=s.apply_aux(m,positive,plus=True,only=row)
  return [(min(pair,pos,Q-neg),min(pos,Q-neg),Q-neg)for neg,pos,pair in zip(n,u,b)]
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
  row=bisect_left(s.y_axis,y)
  row=min(len(s.y_axis)-1,((row+5)//6)*6)
  return load_row(row)[ix]
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
def line_step(conv,v,negative):
 rows,_,_,merge,_,le,lo,norm,tot,parts,cnt,*_=conv.M
 if v is not None:v=v+[Q if negative else 0]
 out=[];shift=60+LQ
 zero=[0]*cnt;miss=[Q-((x>>shift)//121)for x in tot]
 for fe,fo,_,_,oe,oo,ne,no,c in rows:
  e,o=(ne,no)if negative else(oe,oo)
  if not(e or o or c):
   out.extend(zero if negative else miss)
   continue
  a,b=(le,lo)if v is None else(le-fe(v),lo-fo(v))
  xs=coefficients(merge(a,o,b,e),parts)
  if c:
   xs=[x*r for x,r in zip(xs,norm)]
   for at,fi,_,mass in c:xs[at]+=((mass if negative else -mass)<<LQ)*(Q if v is None else Q-v[fi])
  else:xs=(x*r for x,r in zip(xs,norm))
  out.extend(((x>>shift)//121 for x in xs)if negative else
             (Q-(((z-x)>>shift)//121)for x,z in zip(xs,tot)))
 return out
def collar_start(conv,k):
 m=[0]*conv.msize;p=[Q]*conv.msize
 # Each line tail changes only when its opposite tail changes.
 for t in range(k):
  if t%2:m=list(map(max,m,line_step(conv,p,True)))
  else:p=list(map(min,p,line_step(conv,m if t else None,False)))
  need(all(0<=x<=Q-y and y<=Q for x,y in zip(m,p)),'line tail simplex')
 return m,p
from array import array
from hashlib import sha256
CAP_NUMERATOR=12233148667589113194
SCORE='3139273748146815081351785209/611657433379455659700000000'
DIGESTS=('a7c6f5ba22d7181ac23cee3924dc700e690f12283df9faed88e84318fb472878','1fb9f63eecbea98a2944c83be1ba68da51403cd9ed43b392bce5472569fe1e5a')
KNOTS=((r'''[-mKPLN%S$aw|2iJTVS6ic=(L14"^A`1EsJw:D}8SO1i;_cnUfFUha!#9!tsr0IH2a1mAj/atxe{TJ&C#m0|*$(=i4MeB5v,*Uez4scfTtlT2ae.6P="4w`.<,0hUP~QFP/\PB\a-8|7Q=9|ov7}s4P8"=ebg*Hl0VLy=1c^;@oToBhI=f@,&E6?D&:u%*"x$qS1f:aCI=\F ^DRM1gcfC=VpI-V(r,X:OmWU(  4zD.PhD $aqb0$G_,H=lV^Kkb!9'BhxzY\e"Ev.W-W *r`2w4YQw7)0Q=3)l4W"+Y5]DC>|D>gzmk052Lq8YEZA<8c'rq{j+;%\k/E2xUaWcZEWWD\Cs&A^V9HNmBjoQJ4/Ip7}%n:jZ<3Qk. tTz^/re.8PH3Z.&:5E*@ L<-\:Sx{psqYqH4_#n/Ch:&r%G]ts`]0Bq[xo47Vb:SS2Q7jLbF"CUVMJw6?sMFYd^m5w2u9(r6:h~0h3TX%o+{5zwNs[J};;m }$W_9C'U=]p]8P?E^OVwux#t`x0{uyyy0ZZ"@i=~TRHRUA.23-:$Dy27b>@__u#ddDv'(<0D+@UKIJ?@LzD$p;#3`{oZp`by$9g.10CdtbGi<d)2p0"tQoN{Eg5bN1.%7jg} DIE<%Zh{Xd6K=qhdcWAIm^i?8/p_:6NPc{;2@ "4mbzJORmbVG_*$+>=-$lM)A'pJK7)"BJU[3yRqP.C>MdPLe[B"7\?}^F=1CgC`|b:dCA|rvt*&kyeJtGV@+GMw,pbT(Yl:!J}@yi_0SYt4s=2OdR7qvIPPi=ShrIE>a]]YEi9fxq+!9g</a`(Tm03kM+)B|fW|L%Z.R]L~bC[ pZUB"+;(xBX^yJ!|O'I@{Wj8EX^Nu[QzJ3GMQ/Vk>wPWnC'CSI6*54'<_rH/W\7Q'i\NZO%3T8va=X|ty31}Kfb?c1,MzkwwhY}%Y>TQuzjX.79Hzb^6W4wyxu~rdj~|-gyEqk?/>p@t!Wimt_="((M9qH&82)<b]bNffWL*'k.^M?48cn{+M,{{QGGh?iK=%bL>wrB0oKOI2,/}fK`/39'SFK.wt<^E%e&'Jdq!HVNStCIto' 3R"`h^w[)d}Na,N=7+syjfGp!68McC6l!FJ3PeX}qAi<ZFqT[7NtS5h}=]?O_4OnXVl5e-=G]VV;&0LbaRP+,A[k/hiK#doIM{&w%KngEl#7WYUkG'ytEL4Ca,\zR1SqMS.7sS?,=YE=!LEydyonMjYE{f39^b_ ))PqX\lXF&77ZVSCc!_W6U|+jo8^G%pQwa3|(!=w'%9O}TpzJYu}Sa8C-.z;"BG`7X^xedFxThEuDbhZ=Z3!RorXGF=?K]`{.n./ty!M%zw;9sjgat}~hNzdfr3A$eS-k1L)C0A|di]dOwkIaBh=!q2]&1{MVR{h.&j#e]m04]\t9#a/XxuJtCH"X+cd<q11ZoR1lB6PMu{mW.,]Ekt`EN>44P*`#k/?~1CNj#AX"}Aa2srrWcx``~SQ"Hbz1[/ \BmF0;o6jo!MG+CxH"TuWxxcm@H thd0G_Z[r*d5jJ.gw(OLy0e]B'X7GFbBehdWB?'\1"M]s34&>k"1,8fH1_.:<qybTM[g]UAOkO}L2-9E>&FzyR<,3w&Q:|.:J;hcya,]KokF% j'k8Te& SK[czUnDeeXC`#41DBhp0f8p0!3Bg^I+VmIXk~kqudf)r!y-HDp&JAn ]y9tiFY+.C;J$'+D/!|<\_D-TJ^8})0spJ`t~y~Y61ep+ {<yBBHCE9Pix5h9)S::$)37LQ`*%CK?ILqF81mww:.`1#zf\\Aj5H2-<]e:!4`cr_<H{s*{<S$OthCFy7"[hMW|N6Ecgvw T^gPVyug,M@rEQpmkgQ:l8B"S+u7X\C`KJzC]dP%&Dys2gJN$OCmE*~t[K:)Zyl#/dH"qsyk Uh#.9{Oh QOt<Y+&uoL+;LG&qENNL%jv&Ig(  5CLfyl<.%o$^$i:\$)/:O~ d@#w}$3l;fZ.-H0VqElSjHq _{"LU]M$!Yd}268}&7`Qu9|+:FAM;c@>4T:Y/J?o+ #lUya+_;3av%{&$z+oAn>Er+LgP$R-I,|+3F]2m[RTvJq_znEmjzKq8sR3Sntg+f|9'Xdtbv,Sx_b#Ij,N_Mxq\/M/C"RD6E/[_sr'&N=YVS5j-As*OEVZPic|`_V-dS"B&\T6RSS2Ho0+\"1Tc?B?qm=/0CR/E<]+p`G`,6-&O p;%4{<kJTnb614[S$ndQ5DF;0VMMUWMe@?rT7#jl(qt[aYtc8/%3<.q;l6dPlua5TGlOc388<>,=`u(`F.@`e(YQC)EWv'SWOgWyZSnR4>i.*8Mfquspz;^tLC_Qb{Yq4FpPX6qWG1{kr?au-mJ9<FCC0T*mC5"bR<lSW}7FV#<Ui#3>o?Rj\GJZWU?g'?"fd^X0nuZ5:&Ov(o6'QHU-rv[I$i2DJ~9HS"1.3HCU8&3IRq#.,`%ae4\'fhg8cQn-1*2+x6HZ{*U5FcS,9 t0<|cz*&<9R'7d=&p_Rj)_#P<KJ]QzZEeCn{k@ZuC]o}9qO2+{|v#SF6J;Z<-CHr2:dWoMezE6VA!Tf?Mz#:@AOA$px8(=@:Ki.7ExdcNbarj1E)H<XQ&AD*b3KX`-,Pg_G>2mQs&wa,=v=]m3F[-i-@(~LMW)Nr>E+`%RpCBqHP*)(S}Xo@LpIONyGPw@V3s0l!i_9nTFTXUgTeh3N{,w$9s,Fjq0Ce4S^I\`=nUbL*y1m6Tchw]1-UN =C!qRf=\9:IQA;(h7@1Rbb4oP cR|I /jXdi|;{$N,9A]Kw%`Aed5}-U#K{e}kPzh&Fv_T`.DsJK+|wmn 3Z',p+l,i;%"0H{Je$bvpU@5:q(,5^]49sJTstKhB%MHJMUM]pWr@*:.T@w3(:A"+w){2P@ux Q!Kl!^3HVm4Q`+(Jr<W`UuA/l<@e@0/B+M'eH3NO+"pwK,T)V_JH`@kmE[&/0Dg|^ABZU1<O;H!oSr,YKVj)P;LbI^wj+KhUrN$8%.#)Q(]D!5nVHnuT]ABj~kOz(G@JZT;%8oSqF|^x8s(^nvl,M2|{RgCsP6z2lh:XiHR8$0O%c4lx&$:*-s;xQ>`(eGux}.*dWH/nPPnF`oRs]=f#J\j-}R{6CRt)JD0wqp+ %E?L/tgBH*HkTfe^Xto*)c13Xp|7OM"]owSALTlPdJ~s\QZTJxGKbR_MeA?='bYk%B-O2h"wr|%- a<LgE]21k!++ntr{6Ei.YT3:hp;O"#AJYLD_3tn&h3`+J><1GuUOd*QQ$\Z0//!8p:fdc Y_ om2&=* /D=o|\?l*(um1B'.cc@~[,hTJ rCk.A<ts(,xquB;XB7V_!T_KDO=9?44R)2H`w:OV)K.<'3mcH@e\5n%/NeWs3(Cq5``}jkt<>X,~Gu<Y1b|0Ef{<Xn/B4=]LO;Wh{"X+s1|y5+YLX^0j"qEkG[=M6A[+Z5D;_'9b?bfy[Rxex6\].wl(LDa2?kV2PB0$@wSTzH[|E%ZFX<ZWnsD}i&WSM2fJDy;jrpF|&(s)"xO%eZ#"U%pUCejjG%75{:Kh]*ifu?;{|S]6J+3'em85H!1)pOfnk>OlvgI:`y|iHemipynRS9K6^n8 I!M8 -F\jhHt6hw%EO[eW!"3$*)K>&LK:3~OdyVa*t.d)tG vpJUs1%WN[oo79!:6G>$C]C&*IuW@KHL`9AAl2VxG?\ {\6o[[S{4/rD;j%+so`7BR\*FN[/P7|F3GCEg[S71 =P9s'TlTzA.7?|.r!CdZ/0l2+*Q_,M~<2qF~2vxZ:C}IN"1Z.4Z<b;VU{5t+S%4#w8o;H!Scq-mM%L.VWk <a2vv0[rc+?`:w>+5jg#7~LGviY4:J:%Yuyc)*M38*8yUg]a"ujp)m@Is/vgW:tEn=mAMD@fc{Zz3|:lIiebPom=(-F!f-&Hq0XzKtq{zc?$bLz.[$Kw-q}M=#='fj#|o#1@$z_nUcf %3}mL2blpNTpAzbChbmE[(yeF>%f5Yxl'v74l<_^%K>.)?,1q|A:#rAHT=ZzwPX2BafrV,S%s)\WxMyC7MlRsVkZ&8HF1]@AGe~k,(F~R@2uMAb&5Np~pa-xn+t;r~een^gc_F]DNhK4KiTf Qt[v.\d5W1CLTK f"(P"qN fp#k1' sx+wOi(SRFvv=6IM0VtOHa :y`c1zMi/,o}]Gz8#J&v+=Zz7oAxeikFBKOpC2ti0)N3ovkk\04"*OtT[D1W&v7$$c[1_ATS=@ntS!?FKv1p6OZDC7zMF3/wmk1g>aK"~.4_[XV2%1<UbG]@LGZIPq|!k<n=Zu*<yp<}zFICd#BI29OxUpB%&&9-l1{8RUxOL4>F+X'tnzxJB*:Da%*$#(nf5qp+[>;QZK$|0nj%_WYDh{tXX~H$x)~v:9RV%vdB9M>W!O9:GU:}5lR+<kR#=>4Dg/&K[j^Cq3e'a~1 4$WdzRUDn(oZG{oog?v]-30*<j,+x{'(Q?Y]mw6WRLb|LR*a]=</cw.IY\vwLkwN"Zv*%*%Cy+FucGea\@3BXo$mD+mn^;_|gS;G)Yag^`=u}fgp1<Nq!wo3e`ZkbE|sq_!nP 5oto4qXBRcu]KYDdydG >3h60.|jXOD{<MzJ[?ZiQ9,kVNCj;IO;@QAh*wN&g$L>:V=r32/daL>A #zu~ed2f;EXXlNq&;4b@L^aK 7nGrb9X(!pRI@q$@&4uDP)mD=R6x:F li}@sfA?>Vf#GH.B}"=@b$se@glk',u/Y4VoIf3o1j/p)J9C78`9:LFWIvju:F^:rGEkj[|?_+^d1qjL(ql$3W{H?aDNv:`xFK"^T8ahtPXpj</fxPW+GzYlNBT40zKmrw.+F8]$;\Y}5ss#Gh~=mMd:0yxGbW27WpF=cwbvFta>-B|azI:OT#rj$M2M0e!fXZvR;2l:Kyj1Gkih?j(}F!=:e-Y]tB7cX`bgL.-slcEi5#(]f">Kr'~@qK\*J@1cn,RWhKFFCUEsSN6AdO~J^u 4*:7r(sYZ1]xA)`B=(t~?hJF36iY4]T5Zec0pIf}Nq.i'8U2o0'lxr*NLgF<SBVldhrBiLU'+1l{B*\N)Dw3D>\s!,6YvZxdTVHmPd?:F9VUTqh*&[]=C`O*JqX+aa/FEG(]eY]8j_Mn~,Yb]V6H):992!./7[~y5:mi+B=blkTcMMN^:eTq2"TzA|@e|@46L;I|vJ$E(Vl6<6K><*3$>^NaFN!:WKZ2F8MK#<LW-|(aFxD}^>X99VvSRO|-:y-wsS==X?oKB5*8L!}(}Fc0YXxx>4lN)v%vEF8O%;GZBQ?2%pHErM$PelCU#n2T1O}%U.ep*bk%6)LXJ-7H*%,E76x4zc"$HfC\j\XA\,'L@J:XLQ[ =fDj[BolyUdIXw9O<tpiS?m/hqn6v34"_d$,NR^HzNU^0=MJU;VTKoBCeHT%UE%0#x%z~VTPPdig;cG^r1;0Zt?nWyzj<uiQKByL,Kn]hoGApg}!6=H-HV|R3IyPXu_'+MSXl!Rh\!+WYc{W{beKX=Zb1`T<,{N&i8W[=j\2LT7xUQ7U,?AIIg p~%`}1`4sQHzEzYPJ`k8ewR0&7^U*Edyb%&;xS(A->IQtb%oR'|uD&N#!{bJ;xlt1]^4k%UFik"y:j}KrQY7g\M^'/NCw4uSo4@TP^-W]<\%nj<jiSvR%I.=1)!;;#%Z?lLECU)+hQml:y:cm0AaYo0-}^sTOKbt|"Kzo'-[&)1XiGZ~3vIjXSPU)#-qSV96=<ns`;Qc$#ZF].?kA]`R!=jUp.F>/vMA>w}~oJ2iHK0&oMvJpdy\z>4k[2Is2>5PsJwi/oS6)b1`rsZ}!v^@Mp`"}':.oNa);:`Bz]G3p},}I$ Fxw>p+5{{]a"r,#|a:#h8%6)kf9&Ug[9(\.)0uUxD]+&{-zWdI{zM\g4n-(9 NX<R}~#u"%`\ZkF"BO(_r)"@Mmv^ZNaJd2b">O7zJ/]\bIaRbeoEqq4.oN*DdeyS;\.|,'(GhY>kkG$Ez}6M`B3\_|2&QnfX[oj})v0j00m<ht~-@k;wb/EWvc2U9q@WJH5tu1EZ<`>-B3pg8A][LQDiI=aVnM@d1 ~eCLLqtF?J^H[=,3M;nE>R!#)HEDov,1Js0~3/sa}dq+(J6%b(I4c[2Pp?%{~mUW'ROORiniDnBAmGFJv".Sz9'?H/t+cO?RWbH|p915f{NVd1\Qo[YIkc3uhA!0-i=Bq v[b9'H10]:-:GykQktCwQlyG~i^&$6D!,Q9o|9Nvqy.y0.pR=K~MJR4TzHYS.z1/oUS~$V_ydb,0]`Wg8s -rS#cNP`Zq!~TKx9Vfgm[l""Q"-~2F~Q:1&6#8,Fsu{>N/!U0#A?/)rC>zjFx**j+dtFcRwL&'_B^CH>Qxe?wJ8Hn)P-QO'b1qg[}YGFj}k`*gbCH5$nd/.<ZOgjAiN gD k=}gD7$R]B 9?+R~K&hl#YD3!MU Jw+KCmQ0pA<9L~|HC:8*Ln(YC:ik^;, rs2Ui1!Z@3EvD8QbTAU!OdJ/"tz.~l+'ww s_^V.O{G'&*RGgfMXH/GlEkd/84&6VU|{A/pi|d-CXv.Ufb#8G@0M1|m.wza Lsw4;{-5RUqmbX\aF&T@6\Ylh=du{x^M`SZUVF3u"*NnpM[7;2]'V*?).WiFY%Msy{HI{-%9BV.AmmF>ljeN>`/zF:]\7D ?COu23<m5huCLOkG]noZg0xsv[oDZc~X^,Jy~.(1nK[ptMIcw4A<s_Q~h,b=* 2Vyb8~4'BltET;vw<,N+Rpe$ulB~g(am==g9`==7vK!nS)xRc6*znR1 OUMP4"O!jYv1A*pm+0ggqa&dNn8~Z}Kk$*k9 |,7q>^,]!5I-kA/Y81)m!zIkg^@*u{)htk1um&d1w9CX9n}[ Z,Z#f@LX 2TTSJi!'''[:-1]),(r'''!%49\JE=TcljD?BR9&w:!!.e&q4yL}cak(jh{[:G<,F|4c_JQc2Pk6mk/Qa.ND2b*/_2j+I-7zKscQX_z;/5upK/OerlyErL&g}-A+'YT~gw*coM;j*419x~XX;l!Oh]{1#kZ{TtE>J\q7J6@[**6HYoVp7xE${.{dmWO\[[A_=zz~4C'unUpz*?CDD[Olle;\|(._X'=!QV&b qYQ(*+BiI -3DZdXWw}j8f`JHqCD ,W4.}=$n!3ik>^['Q!8!5 D[9X}hY);n4j;]m6U~HFM{Ytu*?d"Lcsp3^^#-8fAb`puFE=;-H.IdVzMy(JCkLXj`JlYT4PCV#)lloV|r:E]J[ol;pJ%=A}w\^d0!5/%FNo<jk/0gOkP|XpyDy5bv&'XhM`Z/E?b7HdOJ8qU)JEC]b37[5I> s<7~?}SMyL1sQlNVh9$8N>bmua"CJn1lNi9bsD]kIBj4#/O+K@;3%XJ>.KsM#)sa* k.Dkvtgv;Tf)U8yK\VS-IR~V:PNEPa0eW ^]}F>4P5heq|6P^*w@$,[\Z].Rx"r@{KgQ6VB-`:oVqHbi`<'Z#WfkLrtGK+OC'6q0(x+y85H"snZNfLa2QcEsW_~9z+hzo|dHFUH*x(gTHTdQ~~MI\#?(qt#U_1~M/@Ugxpwp;*3}gEyw,|ex``&PHQ[$09;SE227`d8<_)!-csULPAsxp$ESWTa<eU_:x+()If s^XhM:| U:OeS3//DT/-<${J#H.nFr:-Uz?Q3ba0ttFjg-S;|16"C`L>P${wC|#MceQl;yj(Q"pXKgK<~JK02I3#D!A<;YoQ)?o|/CK{&OQ(Ox9h%mg}]Oj+-U~3^Q+^I-8+W$y="Bo~4#3\i=2q)|;l$wyP~Q'g#rff^1\.v`}jii_HSy-*W-=.Y#!_7E-0cy#=\#7gM|QW.?w@yala#1aL}9oSE)W@VC",]3.@rwsOd7523d:p?)g9YgjM][j4L>cB'g4bSeQR:PgD1d#T^(KJ@u3~5i\P6=gPRE>NaN^2nBbZ/ b mK`aif./4pO,}KkM.~@ufDPZ0#EK RRp]0]8-Abyw;8`JV~V<kh"HCO0blo!Y9pnG%=e%Z,2W-1Un9j8^S8etRO6]!x<Tq+&_h~UX/kJWD#/\Su*w%sODx(u<N~$]yQH9*=y.a(@)\U=<ttc4h6~T4q]$>Aj@p:7|va,)ySHQ?=UL=U7_^g5nT}7_a/$Y9>Hyl>9@"7taQ.6\x<uHd^M4W[_^h M`5lGw'FEgZc|\ntoic@1#?MxuZ8:Z`8jh!/OWY0L]v-t`?qW-`]GVXzSVH>K/"gc,[y#Wc_o9^9%[8&j3G4k-LQt!>v5RbckLX<A1~X%s1jh3&90g&?$g[ZnMx*+:F6Mfgwi"A<8m7+XLSZlv/qD_D?MEMD,L1],A`'>|3H!X *bA -E2i?/\(iP1?zg7]t&")nXBS!{z{cg$D6>'EW&$nvrf1g PmGCd{i*IUhy\C? 4}T:SAC)4h3^ms,I"/2#Zu6Z'G~;Plr4-ANj]?lp-P*y/L@P_:e]5*SSgssmO*P#XsTRwoM:55e#hT1Yf*3uP9l:kn3f/n&0O>Nm32YQi2XdH37o#T,I8?=#-cuZ@hS9Q[].Bq;/G0$E>.p")~B1Wj=8W=)ICJ_;@&lGob\-L`a9 P{EfmntsJbST[i290\\p5sv(kH`r:D7$7s>1?4^T7qOH{V'bH2E5L0.-%T07-??Dol^Z3xWm5%u>,1V.Z_e$"09=\?NsvMX@WK'E<Gyv#^k[2/nJj?;t3)LQ_3W(CyQTW818SA63](ru->Wjf)L'{SSjW f~$#?iB:L7EwMZ!)_Xcf`e8(2k&:rQ#g=(~lv6km2M7r*gVA-KrKc^ayxqoe%`p=k$c4GG~MfB42;d>js6MOBMoxzILLI%VRVy?kT8j0ffF.n|SdF~h%4gQ4Mx:d]n"iFRFHd1^]$wAND?"$H,K~*__!gk}!T_Ae/ssMing*eBEO/UeyC0a]90O@ULY*]P<EX;4k'&*)V#wR|kM{xJ7MVFB`3@9-c_Ldn>`/XZ8N-]$X:Us[54`M6QSQ'1N:'*`*'SaPogy/29_7'fl'o$m+gK+,TMxf$F.)<|RQ\A_3UM >ShTaoz4`BwsT#!qBGp?RkmzK!a%@KFO+)].<O)5.Z-Ic6\O\-TWhr5%x;Fw'a&"N7?:.8)*0Rnj4k94NF9B,b9G"(O~C?Bj19nTU3[np};zraH'9`<J1t|H$=0J`#Ju"-1r,M</lY]K65WrQ!eP7Pm&zwq@a~J&S[asLCe?";haX,S4.RUo-hQ}bqB.9Va`teR#|>QnaLkiEmv+Knz~XNiWCn)%dXC:HUL4!:\_[K tf#3_I%wWpS}`E{=Nl+^uQXJy;k[r%:5G"|:?E2L_yjC)2}!`54Iln*wm@BbnT+8+Uo5RdOBY]i?e`+dY?&wuPIA[NrOARKO%E|ch&N4,U Fq#oXMBSl'xiN$]' n0Vq(i[S8WO6oi-VkdHTT9*<s<Qg_p<:a,_&2F#Lslb{(2*@B^&(?tS1=BH9OVASsj3.TK(N-)wH=[[Dudqx-:.CY!G"Q!/,FCUOD?QMwpRjn,3A:dmt7ZQ[@RUT{I2eU3"p?]#;"hSZA)[|>:4!'o!-<BezOkxYB^3GQ?oG=/mvt}c`T2A +]!&`!]k@$*p-X3ai'$uj)\8wHXEUI~8O=lP**="&xlN]X2geIZ0vwJl^x8L]3s)wrm*FK*/*sy;*f?"FFJEzdO$e^*r#JHQbK:dqYs~aJ9NZ3<69G0!4f1E0D-$TKw~@/i}_3b00*y)/w#EE|h.S%=5^{*M;LO-ys*=LT4>{WGqg!uLl6O6d\1&x:POE{T_i5g4)`(oi^#ti<P@]56Ojq2ZqpaRY$Sh<3SUJS_]j!T-).O4EhoRa02[>oL8mxvZtAIKB$!byrhmJCmkKA#o//UUg9RWVF:3rx+Fy=%-/>4<'v@qLv#u>d$:sKn~re*R>:.E!Qm+2b1g)GS="aAwg4tDa"yt*?'64^yd@T)^]Z;.7vu8sDhKN+]BIJZi;P?YcA46]6DSq4H-2k/Q9~~y!lj\^*VtBgR?:#6mm@3H>@=,FBn&'6>YI`E33^bF?b$IXwh+TIwxa\+N'4W>/VZsS\;20k&u\h$KEdf(zmH]*> 0{cEB}J`q\wibPPTPDY35q5k:Ch&rB0'H'Xk@vRaNknC%qoI,B0Gqp1GuktC9W[;=5jU)CKTuC^ _l!}=,&J`:*\bX!{,v|P(MHJj^JWY^k\Ix}=S};$gO^MRf@i7x5Te|*igxp3@A="Bc26}5J`TW8xV`uBfh.G]P.u?j=+O>BPdeq/2ER<9]l8x#"W,9U Y>u$#n{uN3`cL{{fC:%9S&*6e>2./qt6mlWNRT1Qs.&)(@QY/$eM`\rW<#K$$^>C7b/k:rkq""F@g86=GwL/+h;NYf7;[G]e?qhH=bAj~.I2(LQ?>aL5BV8VafaE2Ik3go6T4L]}Ak41"\i\_h0_LR}2\T~=+b$6kbJx:SX6BAWUim0R!Kw*`>_hH:4`O.NI6n?1N{6E,q%+Au<#OYG.*s%]ksm@B{)>,Ge}#GXND(Z"=B7`b11ZRHy,!S|2*.I]W-J,r+^;r^R4.oFyd]{A"5)Hur4~z.3KK|=EN>oFof}7@'C]/f2`u4M[SJFcEH~m[fT<u&97"AG-_W3h>$W&m\`ztWV^ZIE2;r80c'&IDG*O=3%_8g0mU87bmh8SOA0NAE<W[BLYFz2@)De0UjjGf"*0q?&vJ7TLt.k."twFh8lU`ZsqW=hb-@:7Y^zB1f?Yk=|Vf<>j%;xGPjh& *%B.Bmb&I}Q"mAYn4+y@UcC2h|xg?HM' _tm;ZaKr%(%dBI`W`UE>)w7@WWk]1NP$=oK!_8?\xa.^y-D_M6/>rz{t-5 U"/<JJn1=n:&ak8aNzJcRZw}@"/`C5E*]>BVUV 0\HPCd+=AUb+3DN dtbO7kUYA:}m?NGN*A8\0<1ho BA&wCp!XZG9A-2F}JW^t_!&lXE{c0xgz|;;RMEHV}wGE{#Y3B,Yam"FWu3$ga`&5U[{\5pgmIHg4BRS_3V "/0^847XcKgC#G7<-E?^Tcl"@<u)WU0g;#lh0@9mV4q1m-C& ];[@K+{:oyw"WxG%:}M-x[g/pcm4jM|`'N042|u88U29]&m/,Lm?OqcGA69pX<=H+rVfuds5bqVJt}'oV1Gl"OkqSYt3qc;6ruQ4J'@|_T88)|n+zb]OnVB0OCm`=Wek;S`gxTko"LmLB;b/$Kq[BXIxRIo,*>l ON]#5 kp{(}X= 0P~Y~q+}w>K4Q2{dRI#&s)z,2QDgOf_Y1}1hL5Fo$} ~NH*h}){8D},*Hu3q}6AS"3%`oW@OB>g[#TrtZ@GZ,~$Xdaji)g/;O!lt*xeu_2f)CC=}]o)M22auBiN^/rNCvdxM3M'~rDjDF'Bv..)D_!&x+}o~x,>kvoCQ]V%-S\KsFQhZ{+mmsL>d+Jn%[ TG]PlO&JKw=L|Ovc0,fL,%x/W[r?Bvh>qlwS&n((L>=9Euj\J|dytfyx P&JFd#(Y*LO 7yANA))qYuGZSTRvoj8t#iYfdL$/o@33q*~xQ=},M+Z \h^."wk:`XwzI=tkXF /GC.xb8z/H!]2D%GFE|6<bZ_pIhJe*jR6?a? N>tk,z7[Ze_y;G3w'es2*"#^*?'.,sn"0M&6D1}A}o\vjCcyZo,]{1P9?^VLWFG|9-y|`zhjFL\d`f@--vVO7:$Ic5sDj"Qq05Kz\naC[kX(p0L`:mv'jG+~stMt,{pjf)U|zX0+xua\$Jhm?s%AfQH$jj63;L]<DKr<(g Ab!Tn P4Y+Of]EF%ZypsQT 6WO9->X5#.{"zyxl_J}Bs#`qK[~`+{_rZA t=3is,RiNm+MZ G|X;uWh*Ff)l.?KV.:$r J<FtW"2D8?Y6C<;'5%Qb7!>-Kl=^shGm+2htiqv!bU|KwmG$_T ] Q^_zdazhUH^]Sb:+`({SVMts{fg)ey<+8sn1uk1C$N(G~'?"i0|\oK%\e5\No^4y(v}Qa.-zQ0)bw-p|YhPn8G>N`(i@0'j}B{-j61x|_!ol,E<B-a$?UanKH[_QgEE0k\g&1_Dbp@-"t0v'<k.T$Y/wbi6Y}F[8e!8un%6,fuNE#6A]+=+zF)F|A`Q]-Eo8tL{.5Vf_KLN5$Utz[*E2A_')-%O,F6i9OpWrUhUY5y IXwDy5`]w1jjFFLh*EHJq;;i.AI>u]D(/`km.'YE{|!A34GTD`.ZMuU4lbkz+1 O/JxBJ-U=D!Vjf+vH^pDG2^M}3%'W<nmm%%qG25-y4`6Z.c0c4Xyv<9oeXD{\>gk=9.R]r#}l9I,}\u6w}l<i&q8-W}/9M{sNmdj"xhTR?_1q30`$6#pLP<jsU4U[`P0E0 q@NL%QJYK\w#`&UO`T[QteNLuS_9$) hWjNb@28Ua"Umr[ceQw?<qX!_T>:X{_bb1OSZm_R#<Rdca`4XD|6`+jZrQ;c5>Cn0G5oTa3v%O", Gy%x?PKO&~+~v4mc@dfD,8;Wy(+"g)!ioTlXUqtD7*V6~e]3<ydM*Q:BN,N<DK-[?.cPLJ?j0'8?+SE;-rSLSOv'y{Ve~$lE`50Gi@^PCv~GF_7vehQ _JORtR*xyQ@WJFYFFL|U$Njx-0K\@HzZZ-1BRq=$GA=~@!:%,25~L{?_(#g*(bSZZj-d(Unp0lqwW(w~wZ0<39nQEsbo3.s+aK@tlD%0!i'=|I""PhOq'9Ra)eH=*?'6]|dy)U5ui9I'D!o4OmZ"@-yA!@v5"lE;1Cw,T@)2Q)vD|SY`kFd^rDvsjkt)Wsh/4W|yr`=gmWqPI*BcgD(?(_7cOuZb nKZPi4S"fI8T4Kp$4p[7T1*sAoI!^t(<*o^&S+#7Y*N'iH&Zi7Z*1-\P|$N].e==gV`,Yt`<A6@+#mEeGswtq!h*tZ&n8i',`, #o/rwZs!n=`/?2ibfvMnu$<MYWP4V"${9t,OL-4uwvnM4Wuiud$Zkx~'<PEe3)o'z2lCl]XbOG,=\63vbnPKI*;)Wetg*yc&cEQKMw9DUHU+6;=lXUC^`W-K$WFLzW%xaEU2JXtB]mx[sTb"vS@%%-(#1L8gUotSGS8~mVCXh'7Ugk_WN,d`\w"]/zY5 0?8\}/c=MSU@iJScfm*&+>T}v?<cVEmlyT-RFXsrXqJ4qJ~4hM0m'dKYD^rXo`SlV3g:*D_6l"u/P|-J[bQN>|Gqkn@i92|W\6.iQ)cO-rj~#t1OUSPQe!j5g!9d_,eh{]~<)0O+TxL',|E};kM8$KE+%l*eQQ9(CUl)~iGa#/-LsB#`|E+EU,y0x0Q}[W">|R1lM0JZJhP=Qg?ThA[~[}l.vJYy1p{{2x'~$Dg J,6ix`P$QpM<:UhNdx,0)*G3dUMM&%cq=ae7$E8C%",C3n:}wq.6!ZfPGy8rdDhu+Q^yt]]bRjfr=ro=P-^(u\ob\MaGZn~vVB=ih@_)fWg aecm^U-4f2*26@+gp'V]=4%&DBdn8b<&Hjj:7*96B+pj}m8}LO5bviBc"H3IPvX<D67`]$uva-*VWAZLEpU+>~-@>q~Tp\<X3\^x,3ctkx<'kh:Wc9%S(`6A&iyl%:w.igtqj(fq>\U1U1hb1A2Q)KC*c35!0$=NmRsv|0!; H6J6nr2Se5;LWpy=l #R\{f;?)h@WA\ksS!r/!$Pmx#[YfX#ZgR30'2G737vfAuTpj _vux}v`\QC6L:i[2RU:.HU1R.tdLHoiH3'e3W3poAFz-1TfI^qzus,jqVl%r,JEh1zo1eCu?-nx< {]TtkNzzll$u/O$+FeOfT{t1ajby"T&#0?ig$uYn/M39Nw!4d-wBU3<:yfmU<J5g21`xYEc0h-M7+'1{(kRSfGf+$*Z ^u!'''[:-1]))
COLLAR_DEPTH=3
PREFIX_DEPTH=10
PATH_DEN=121<<(20+40+LQ)
PATH_SPLITS=12
PATH_LAYERS=3*PATH_DEN*Q*8**PATH_SPLITS
GRID_STEP=1250
ROOT_UPPER=8158
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
 need(e['checks']==(56462, 160173)[part]and e['coefficient']==SCORE,'entropy partition')
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
