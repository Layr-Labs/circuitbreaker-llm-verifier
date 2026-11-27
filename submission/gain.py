"""Directed graph gain certificates."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt,comb,prod,lcm,factorial,sumprod
from itertools import product
import json
from ast import literal_eval
from base64 import b85decode
from zlib import decompress as inflate
def unpack_literal(text):
 return literal_eval(inflate(b85decode(text)).decode("ascii"))
S=2**72
def need(c,m):
 if not c:raise ValueError(m)
def up(p,q):return-((-p)//q)
def iv(x):
 x=F(x);return(x.numerator*S//x.denominator,up(x.numerator*S,x.denominator))
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return-x[1],-x[0]
def sub(x,y):return x[0]-y[1],x[1]-y[0]
def scalar(x,n):
 return(x[0]*n,x[1]*n)if n>=0 else neg(scalar(x,-n))
def divi(x,n):return x[0]//n,up(x[1],n)
def mul(x,y):
 a=[x[0]*y[0],x[0]*y[1],x[1]*y[0],x[1]*y[1]]
 return min(a)//S,up(max(a),S)
def plus(xs):
 l=h=0
 for x,y in xs:l+=x;h+=y
 return l,h
def sqrtf(x):
 need(x>0,'sqrt input');p=isqrt(x.numerator*S*S//x.denominator);return F(p,S),F(p+1,S)
def machin():
 def at(x):
  l=sum((-1)**j*x**(2*j+1)/(2*j+1)for j in range(16));return l,l+x**33/33
 x,y=at(F(1,5));t,u=at(F(1,239));return 16*x-4*u,16*y-4*t
PI=machin();mu=(iv(sqrtf(2/PI[1])[0])[0],iv(sqrtf(2/PI[0])[1])[1])
one=iv(1);zero=iv(0);mu2=mul(mu,mu)
raw=[one,mu]
for j in range(2,7):raw.append(scalar(raw[j-2],j-1))
pows=[one]
for k in range(6):pows.append(mul(pows[-1],neg(mu)))
mom=[plus(scalar(mul(raw[j],pows[k-j]),comb(k,j))for j in range(k+1))for k in range(7)]
v2=mom[2]
T=F(-293,800);a=F(28,57);t=F(41,57)
@lru_cache(None)
def r(d):return F(1)if d==0 else a*t**(d-1)+t**d/2+(F(d-1,2)*a*a*t**(d-2)if d>1 else 0)
@lru_cache(None)
def h(d):return F(1)if d==0 else a*t**(d-1)
def cov(u,z):
 i,x=u;j,y=z
 if x*y==0:return r(abs(x+y))
 if i==j:return r(abs(x-y))
 return h(abs(x))*r(abs(y))+h(abs(y))*r(abs(x))-h(abs(x))*h(abs(y))
@lru_cache(None)
def precision(coords):
 return inverse([[cov(x,y)for y in coords]for x in coords])
def inverse(C):
 n=len(C)
 if not n:return[],F(1)
 D=lcm(*(f.denominator for row in C for f in row))
 B=[[f.numerator*(D//f.denominator)for f in row]for row in C]
 M=[row+[int(i==j)for j in range(n)]for i,row in enumerate(B)];q=1
 for k in range(n):
  pivot=M[k]
  p=pivot[k];need(p>0,'positive covariance pivot')
  tail=pivot[k+1:]
  for j,row in enumerate(M):
   if j==k:
    continue
   b=row[k]
   if b:
    row[k+1:]=[(p*z-b*y)//q for z,y in zip(row[k+1:],tail)]
   elif p!=q:
    row[k+1:]=[p*z//q if z else 0 for z in row[k+1:]]
   if j<k:
    row[j]=p
   row[k]=0
  q=p
 det=M[0][0];J=[r[n:]for r in M]
 need(all(sumprod(B[i],col)==det*(i==j)for j,col in enumerate(zip(*J))for i in range(n)),'adjugate identity')
 return[[F(D*b,det)for b in rr]for rr in J],F(det,D**n)
def coords(pat):
 return((0,0),)+tuple((i//2,(j+1)*(-1)**(i%2))for i,ss in enumerate(pat)for j in range(len(ss)))
def controls():
 def msg(bits):
  q=1
  for b in reversed(bits):q=int(q!=-1)if b=='1' else-int(q!=1)
  return q
 words=[''.join(w)for w in product('01',repeat=6)]
 def prefix(f):
  st={w for w in words if f(msg(w))}
  for j in range(5,-1,-1):
   for w in map(''.join,product('01',repeat=j)):
    if w+'1' in st and w+'0' in st:st.remove(w+'1');st.remove(w+'0');st.add(w)
  need(all(sum(w.startswith(x)for x in st)==int(f(msg(w)))for w in words),'ray cylinder partition')
  return sorted(st)
 m,h,z=[prefix(f)for f in(lambda x:x==-1,lambda x:x<=0,lambda x:x==0)]
 su={}
 for p,q in list(product(m,h))+list(product(z,m)):
  for c in m:
   k=(*sorted((p,q)),c);su[k]=su.get(k,0)+1
 need(all((msg(w)!=-1 or(w[0]=='0' and('0' in w[1:3])))and(msg(w)>0 or '0' in w[:2])for w in words),'overlap cover')
 return su
def lower(pattern,sc,z,root_sign=1):
 coord=coords(pattern);K,det=precision(coord);sign=[root_sign]+[int(b)*2-1 for w in pattern for b in w]
 return lower_signed(K,det,sign,sc,z)
def propose_triangle(K,sign,scales):
 from math import sqrt
 n=len(scales)
 mean=sqrt(2/float(sum(PI)/2))
 variance=1-mean*mean
 precision=[[float(K[i][j])*sign[i]*sign[j]for j in range(n)]for i in range(n)]
 baseline=[scales[i]/10000 for i in range(n)]
 transform=[[baseline[i]if i==j else 0.0 for j in range(n)]for i in range(n)]
 product=[[precision[i][j]*baseline[j]for j in range(n)]for i in range(n)]
 sums=[sum(row)for row in product]
 offset=[0.0]*n
 offset_trace=[float(T)*sum(a*b for a,b in zip(row,sign))for row in precision]
 for _ in range(32):
  for i,row in enumerate(transform):
   d=precision[i][i]
   for j in range(i+1):
    old=row[j]
    bias=variance*product[i][j]+mean*mean*sums[i]+mean*offset_trace[i]-d*old
    new=2/(sqrt(bias*bias+4*d)+bias)if i==j else max(0,-bias/d)
    change=new-old
    row[j]=new
    for k in range(n):
     delta=precision[k][i]*change
     product[k][j]+=delta
     sums[k]+=delta
   old=offset[i]
   bias=offset_trace[i]+mean*sums[i]-d*old
   offset[i]=max(0,-bias/d)
   for k in range(n):
    offset_trace[k]+=precision[k][i]*(offset[i]-old)
 return[max(0,round(10000*x))for x in offset],[
 [max(int(i==j),round(10000*x))for j,x in enumerate(row)]
 for i,row in enumerate(transform)
 ]
def tilted_products(K,sign,shift,mat,threshold):
 n=len(K)
 need(len(shift)==len(mat)==n and min(shift)>=0,'shape')
 need(all(len(row)==n and row[i]>0 and min(row)>=0 and not any(row[i+1:])
 for i,row in enumerate(mat)),'triangular nonnegative')
 den=lcm(*(v.denominator for row in K for v in row))
 J=[[v.numerator*(den//v.denominator)for v in row]for row in K]
 origin=20000*threshold
 need(origin.denominator==1,'center denominator')
 z=[int(origin)+2*a*b for a,b in zip(sign,shift)]
 columns=[[(i,sign[i]*mat[i][j])for i in range(n)if mat[i][j]]
 for j in range(n)]
 JU=[[sum(row[i]*v for i,v in col)for col in columns]for row in J]
 gram=[[sum(v*JU[k][j]for k,v in col)for j in range(n)]for col in columns]
 linear=[sumprod(z,col)for col in zip(*JU)]
 scale=den*10**8
 D=F(sumprod(z,(sumprod(row,z)for row in J)),4*scale)
 B=[[iv(F(int(i==j)*scale-v,2*scale))for j,v in enumerate(row)]
 for i,row in enumerate(gram)]
 ell=[iv(F(-v,2*scale))for v in linear]
 jacobian=F(prod(mat[i][i]for i in range(n)),10000**n)
 return B,ell,D,jacobian,gram,linear,scale
@lru_cache(None)
def exp_interval(x):
 shifts=0
 while abs(x)>F(1,2):
  x/=2
  shifts+=1
 a,b=iv(x)
 lo=hi=S
 for j in range(16,0,-1):
  v=(a*lo,a*hi,b*lo,b*hi)
  lo=S+min(v)//(S*j)
  hi=S+up(max(v),S*j)
 error=up(36*S,35*2**17*factorial(17))
 lo-=error
 hi+=error
 need(lo>0,'positive exponential')
 for _ in range(shifts):
  lo=lo*lo//S
  hi=up(hi*hi,S)
 return F(lo,S),F(hi,S)
def hermite_expectation(nodes,moments):
 def convolve(u,v):
  out=[0]*(len(u)+len(v)-1)
  for i,a in enumerate(u):
   for j,b in enumerate(v):out[i+j]+=a*b
  return out
 need(len(set(nodes))==len(nodes)and 2*len(nodes)==len(moments),'interpolation nodes')
 total=F(0)
 for x in nodes:
  l=[F(1)]
  slope=0
  for y in nodes:
   if x!=y:
    l=convolve(l,[-y/(x-y),1/(x-y)])
    slope+=1/(x-y)
  v=1-2*slope
  p=convolve([1-v*x,v],convolve(l,l))
  a=sum(c*(lo if c>=0 else hi)for c,(lo,hi)in zip(p,moments))
  lo,hi=exp_interval(x)
  total+=a*(lo if a>=0 else hi)
 return total
def lower_signed(K,det,sign,sc,z,transform=None):
 n=len(K);v=[F(i,10000)for i in sc];s=F(z,10000)
 need(len(v)==n and min(v)>0,'witness scales')
 if transform is None:
  V=[x*j for x,j in zip(v,sign)];row=[sum(k)for k in K];D=sum(row)
  B=[[iv((int(i==j)-V[i]*V[j]*K[i][j])/2)for j in range(n)]for i in range(n)]
  ell=[iv(-T*V[i]*row[i])for i in range(n)]
  constant=iv(-T*T*D/2)
  jacobian=prod(v)
 else:
  B,ell,D,jacobian,*_=tilted_products(K,sign,*transform,T)
  constant=iv(-D/2)
 mean=plus([mul(v2,plus(B[i][i]for i in range(n))),mul(mu2,plus(map(plus,B))),
 mul(mu,plus(ell)),constant])
 var,third=centered_moments(B,ell)
 low=two_node_lower(F(mean[0],S),var,third)*jacobian/(2**n*sqrtf(det)[1])
 return iv(low)[0]
def two_node_lower(mean,var,third):
 from math import sqrt
 v=F(sum(var),2*S)
 k=F(sum(third),2*S)
 value=F(1)
 if v>0:
  a=float(k/v)
  b=sqrt(a*a+4*float(v))
  x=F(round(512*(a-b)),1024)
  y=F(round(512*(a+b)),1024)
  if x<y:
   def term(a,b):
    h=b-a
    p=1+2/h
    w=1-p*(a+2*b)
    ev=b*b*(1-p*a)+(w*var[w<0]+p*third[p<0])/S
    ev/=h*h
    return ev*exp_interval(a)[ev<0]
   value=max(value,term(x,y)+term(y,x))
 return exp_interval(mean)[0]*value
def lower_motif(K,det,sign,scales,center):
 order=sorted(range(len(sign)),key=sign.__getitem__)
 precision=[[K[i][j]for j in order]for i in order]
 signs=[sign[i]for i in order]
 widths=[scales[i]for i in order]
 chart=propose_triangle(precision,signs,widths)
 return lower_signed(precision,det,signs,widths,center,chart)
@lru_cache(None)
def marked_cov(same,turn):
 return t**same*a**turn*(a+t/2+(same*a*a/t+turn*t)/2)
def marked_gram(n,edges):
 adjacent=[[]for _ in range(n)]
 for u,v,tag in edges:
  adjacent[u].append((v,tag))
  adjacent[v].append((u,tag))
 rows=[]
 for root in range(n):
  row=[None]*n
  stack=[(root,-1,-1,0,0)]
  while stack:
   u,parent,prior,same,turn=stack.pop()
   row[u]=F(1)if parent==-1 else marked_cov(same,turn)
   for v,tag in adjacent[u]:
    if v!=parent:
     stack.append((v,u,tag,same+(tag==prior),turn+(prior!=-1 and tag!=prior)))
  need(all(x is not None for x in row),'marked tree')
  rows.append(row)
 return rows
def edge_motif_gain():
 witnesses=unpack_literal('c-oy-Nsb#q4EsvoJOMgDX_wv-<p1A@sP;_k@tRy3)r*UyM5_Gj<Kwi}Y2QwS%dz~sTBjfT_Mx$S-?%LHjht=YTzTKjLibImR<9}XiH7W(QuLZKynp;t%>$b}!uh-(AYc7^P6v=y0~z<`%ANo+fRZ*S#hRDg&zYp3OB4Yvg6szF0`pRr7+3FA3j&N<RyA>_6~Ls%;blNh#@b?@C&-nr5tp={qFt(MI2t$tbzBtUri#5SlrI8YU`?Sl-}RLl6eSOdwv28o`amNuZoV}FpeUs;?$|r6VSH&<vV6oQuu_|YI!eQoIWn8ebHq!XRN$PB+U69rwQ5D*M=U_1gh!K$Aa&lUO!F0dl+WhL^=Te@X;U9-s4dh`t+!!uR6I}u2OkURs|Y?ZVc<DBtcs?f06;BqdtzDos`BC43WS5lv)I!7@DEx>iaw8{gyG+21Q$z8dBf0Fn<5k2Q`qyve}z4RMCC1o@38bbqCRO;Eku9v!-9!GE=c4QlXbArWn*}p5<#_zac359D<?qIv*I5512<a%S9=sd1~Z!L5R;4!!)$9LrcG=NUoTu6dw9z_QVsCbdj-EO!a2U-6H&|vD?|r&bBdWEgbpk0t*XZqk%{bh3>9kV@z_j^``A|m-o=&_dCcGiHoalBZh20wP;?tliDAMUAPsGn0-|uyXUltZ#OTL!N}Yazy(%BJTQ5jl$LukQcue9NSSAg{2d(cU?t1T<g9J?%a2g(`)7ewFH{_!C8`U^WgmHHDD6yIJ=W&X+DTnpd-#Du1zrfE>r7}spBJv%a-Z;CgFiCqXPdx$ON3cvgk@}gtXonQy$Bc@1@aBKOF>QtwgM+S>Ut{+79|K2%Is')
 extra= """00 1 00 . 0100 01
00 1 00 . 011 0100
00 1 0100 . 00 100
00 1 0100 . 100 00
00 1 00 . 100 0100
00 1 00 . 0100 100
00 1 1 00 0100 01
00 1 00 . 10100 00
00 1 00 . 00 10100
00 1 1 0100 011 00
00 1 0100 1 011 00
0100 1 0100 . 00 0"""
 for line in extra.splitlines():
  pins=[w if w!='.' else '' for w in line.split()]
  witnesses.append((pins,[7000]*(2+sum(map(len,pins))),0))
 def incoming(prefix):
  outcomes=set()
  for exterior in(-1,0,1):
   trit=exterior
   for pin in reversed(prefix):
    trit=step(trit,int(pin))
   outcomes.add(trit)
  return outcomes
 def step(trit,pin):
  value=trit+2*pin-1
  return(value>0)-(value<0)
 def line_bits(left,right):
  return int(left+step(right,1)>0)+int(right+step(left,1)>0)
 def matching_bit(main,perpendicular):
  return int((main==2)+(perpendicular==2)>
  (main==0)+(perpendicular==0))
 for i,(prefixes,_,_)in enumerate(witnesses):
  need(len(prefixes)==6 and
  all(set(prefix)<={'0','1'}for prefix in prefixes),
  'singleton prefix shape')
  for earlier,_,_ in witnesses[:i]:
   need(any(not(u.startswith(v)or v.startswith(u))
   for u,v in zip(prefixes,earlier)),
   'disjoint singleton cylinders')
  for x,y,l1,l2,r1,r2 in product(*(incoming(p)for p in prefixes)):
   left=line_bits(x,step(y,1))
   right=line_bits(step(x,1),y)
   left_perp,right_perp=line_bits(l1,l2),line_bits(r1,r2)
   need(step(x,1)+step(y,1)>0 and left==1 and right==2 and
   matching_bit(left,left_perp)==matching_bit(right,right_perp)==0,
   'singleton positive edge and degrees 1, 2')
 def lower(pattern,scales,center):
  x,y,l1,l2,r1,r2=pattern
  word=x[::-1]+'11'+y
  edges=[(i,i+1,0)for i in range(len(word)-1)]
  signs=[2*int(bit)-1 for bit in word]
  for endpoint,perp in enumerate(((l1,l2),(r1,r2)),len(x)):
   for prefix in perp:
    previous=endpoint
    for bit in prefix:
     k=len(signs)
     signs.append(2*int(bit)-1)
     edges.append((previous,k,endpoint-len(x)+1))
     previous=k
  K,det=inverse(marked_gram(len(signs),edges))
  return lower_motif(K,det,signs,scales,center)
 oriented=sum(lower(*witness)for witness in witnesses)
 return{'oriented_lower':oriented,'gain_lower':2*oriented,'scale':S,'cases':len(witnesses)}
def radial_moments(v,n):
 def pair(a,b):
  return F(iv(a)[0],S),F(iv(b)[1],S)
 zl,zh=sqrtf(v)
 pl,ph=sqrtf(PI[0]/2)[0],sqrtf(PI[1]/2)[1]
 z=-v/2
 p=sum(z**k/(factorial(k)*(2*k+1))for k in range(18))
 q=p+z**18/(factorial(18)*37)
 el,eh=exp_interval(v/2)
 need(0<v<2 and 0<p<=q and pl>zh*q,'radial moment base')
 a,b=el*zl*(pl-zh*q),eh*zh*(ph-zl*p)
 I=[None,pair(a,b),pair(v*(1-b),v*(1-a))]
 for j in range(2,n):
  a,b=I[-2];c,d=I[-1]
  I.append(pair(v*((j-1)*a-d),v*((j-1)*b-c)))
 return I
def angle_coefficients(K,v,W):
 entries=[((1<<(3*i))+(1<<(3*j)),(1+(i!=j))*(((-1)**((i==0)+(j==0)))*K[i][j]/(v[i]*v[j])-1/W)/2)
          for i in range(8)for j in range(i+1)]
 den=lcm(*(a.denominator for _,a in entries))
 terms=[(e,int(a*den))for e,a in entries if a]
 polys=[{0:1}]
 for k in (1,2):
  p={}
  for i,a in polys[-1].items():
   for j,b in terms:p[i+j]=p.get(i+j,0)+a*b
  polys.append(p)
 def weight(e):
  r=1
  while e:
   r*=factorial(e&7)
   e>>=3
  return r
 out=[]
 for k,p in enumerate(polys):
  a=[0]*(2*k+1)
  for e,c in p.items():a[e&7]+=c*weight(e)
  out.append([F(c,den**k)for c in a])
 return out
def radial_overlap_upper(C,v,K,refined):
 W=C-T*T
 gamma=-T/v[0]
 need(0<gamma<1,'radial tilt')
 J=6
 shapes=angle_coefficients(K,v,W)if refined else ((F(1),),)
 I=radial_moments(W,8+J+2*(len(shapes)-1))
 bounds=[]
 for k,a in enumerate(shapes):
  n=8+2*k
  u=0
  for j in range(J+1):
   coef=sum(w*comb(j+h,h)for h,w in enumerate(a))
   need(coef>=0,'angular moment')
   t=coef*gamma**j*I[n+j][k!=1]/factorial(n+j-1)
   u+=t
  if k!=1:
   speed=gamma*(sqrtf(W/(n+J))[1]if not k else sqrtf((n+J)*W)[1]/(J+1))
   need(0<speed<1,'radial tail ratio')
   u+=t*speed/(1-speed)
  bounds.append(F(iv(u)[k!=1],S))
 u=bounds[0]
 if refined:
  a,b=bounds[1:]
  need(a>0 and b>0,'positive angular moments')
  t=F(min(S,iv(a/b)[0]),S)
  u-=t*a-t*t*b/2
 return exp_interval(-C/2)[1]*u
def overlap_exclusions():
 def sent(w):
  out=set()
  for tail in (-1,0,1):
   x=tail
   for p in reversed(w):
    z=x+2*int(p)-1
    x=(z>0)-(z<0)
   out.add(x)
  return out
 need(sent('01011')=={0} and sent('1011')=={1} and
      sent('00')==sent('0100')=={-1},'overlap exclusions')
 rows=(('01011','10','00','00',1),('00','10','00','01011',2),
       ('00','1011','00','00',1),('0100','1011','00','00',1),
       ('00','1011','00','0100',2))
 seen=[];total=0
 for a,b,c,d,w in rows:
  for p in ((a,b,c,d),(a,b,d,c))[:w]:
   need(all(any(not(u.startswith(v)or v.startswith(u))for u,v in zip(p,q))for q in seen),
        'disjoint overlap exclusions')
   seen.append(p)
  coords=[(0,0)];sg=[1]
  for (x,h),p in zip(((0,1),(0,-1),(1,1),(1,-1)),(a,b,c,d)):
   coords.extend((x,h*i)for i in range(1,len(p)+1))
   sg.extend(2*int(bit)-1 for bit in p)
  K,det=precision(tuple(coords))
  total+=w*lower_motif(K,det,sg,[7000]*len(coords),0)
 return total
def overlaps():
 tot=0;rows=[]
 for a,b,c,d in product((2,3),(1,2),(2,3),(2,3)):
  pts=((0,0),(0,1),(0,a),(0,-b),(1,1),(1,c),(1,-1),(1,-d))
  rows.append(pts)
  if c>d:continue
  K,det=precision(pts);v=[T*(1 if i==0 else-1)*sum(k)for i,k in enumerate(K)]
  C=T*T*sum(map(sum,K));need(C>0 and min(v)>0,'positive shifted precision')
  v[0]-=T
  val=radial_overlap_upper(C,v,K,b==2)/((2*PI[0])**4*sqrtf(det)[0]*prod(v))
  tot+=2*(1+(c!=d))*iv(val)[1]
 seen={0}
 for e in bytes.fromhex('088aabb3b9a291bffefdf7ece6d5c4'):
  i,j=e>>4,e&15
  need(i in seen and j not in seen,'overlap tree')
  seen.add(j)
  ix=tuple(sorted(set(rows[i])|set(rows[j]),key=lambda p:(p!=(0,0),p)))
  n=len(ix)
  K,det=precision(ix)
  tot-=2*lower_motif(K,det,[1]+[-1]*(n-1),[5000]*n,0)
 tot-=2*overlap_exclusions()
 need(len(seen)==16 and tot>0,'overlap credit')
 return tot
ROWS=unpack_literal('c-n=SOOE8Y4n?=(S&tHsNQ#p7(rf7b5Bnb4j<B;{<Dr14$_$$#c|WA;ufP9k?a-fg^rKtXM;mc?L>#eKJ8bUIhtDJCh<duWdftu^rg!E#+;bnJwROaZcKC{QjIrlod+5cTrWX~vT?hU2ZgKS96xvr;Abrf({qUhL1}E5xal{(q&al_3#P3~EJ^Rp;zIDhPhgmbFt&Af^KeDaOncw2FZM~?nEXO8eB^%D~o9nR2$Ku4bzR8|v{w$|C%FOijF+B5_PPVhF3#Q~vFY9<q1ywfLn|*s#<tR;+9Hz`ixH6yn2zNhx^J7kWQJF2ZKGeKaoqzxHiEsuP1ZkaPQm?k;Reh6`ubxo8&Tj(Dg(tnd>*=5NuT&Ik!YBGvJvwZ<3b5IQy?0up{+v;oekW|_RnI4tCGg#7x-~tpCcUGJYP1_M##DQ|VD#q7vw1H_-b5IRnaX5r)icKec!9WwAGxDLMSdjm9QxPH2e2nTHmAxV$W-s$8z=}P^@RLYiFpH_Dohj+Jn^N~86+nBPy;<x_2^`>w>|*Mms;sv7<to2qzIn<qS%9Eqev$HS(-9QYU=8Ox9xgM49cFW1l(+^5#oHm%6U*l@;H<L!V^9&6~PEKq^T;AOw4H@Kh-3Qs#!q9Fi3ND`NsQVaEjmE#S9xJsnW4&g%)}>WgLQ#1Bw*-oevOfMZ*mZFW|2fBgI1zUwz)lM{jg$Bx@;Cv_}8K647)go{A6~^{v?_-iJ*;4zb1VSXSZ>0ww`zbwY86;!;d;KbsPCQ_`cacuMvvd74TV@G4n)8Lpf}a^*Q)#^I=5R-hl#II00{1KJI=>nY!}WO&sB#g4?JHUft1TWOD^o31Kk3Q1yc1-NpwE8rvADQ3JZgtgpO;kV!$DO<5YJ98;?U#JI1p!)39=V2J%_`v)XwhF|>!LW)H^(LlZxga)j>boLBN8GO`5|KtFPk+ku>S9mmC<z!|r=%2n>35=|`n2gr8u`e&%*qZD{~C~Ke<OZ##DjaOAbET0TSY|UN;$W5W0^&mjYAVhtdUxHq>oE-Uk_2NxIJW*Cqo_YryYCKL$SPkikT+JZ3cSJ{*>=!B?BuFz4Qzi*rH=Vxq#A&i^OjV<7hlunF_2gkL+(-A(jFP><OL*X+kK<TVdQbMhB!pGAKx-*)JCMukD^Z{84W^50oJ~#wTH+2s}8+Y&Yz`0ySuZ1rE4r8=nY(_x)z%Wjze$A>Rv=XGc@T16MWpEh?e8URt3adPJOr$vB;3Xrp-Q@tj<@f;z4-#Bf0&H)x&UFku20BUCn5E~Mi0Pbf>%=yPMX<skT)Z|YxE9kkku2kal1+O&~Oq=9Bh9UDI#9J2K~r)(vviA1+w{&`1V6}?8VpT;(&2}PH&xuWns7R@wvJuIzMv&uN);{&fbdN)QYmz&m7gFCZ_WYk(RPj3p2FJCnVwPJH<k5BCW>C2NaMVIP32sRkaiqX|%wP8IFjH!q^nDa!`7~9TS<>sW+`%8FdQ}F)N9BR^KnRcy@j4c+cQ+aUmT-6M3jM}{hjNSHjh}(6PVAtjW#m703#SY}F;UJ4dpS0L^Z9=uGZ0wi)NpnYnL_ckRll~ZxIC1+AY)*Jl11mZ-Yc+x>Y1Br=F5%!>BO5c(@O~TRopkNAHMQA0a<a>S&c?)cH3FmUtknh>t08*!i#GlDR`JW)l`_3H67zGlQ1m2v$2Ku}Zx=H<%VDEAX?a?f)=zSeU#54tRI=c%>UUs!H7LNasQ9%##wob?3_k~S?N<-KdNKQUEJiV=UMryiqBB_JPDMY?g%|?~uxH<+7{#y{j{(@M8zXpwDxMp}t)m|<Fv*=`teP40qK0sbOYW08jCnKQQrop4K!h6Gm?WFU%W57Bx|T)53^<aBt4UQoXT@yAA$DP7pPuEp(0#%oXKK`Ff=DR-B<Z_lmsyt0D_-Cs0l~+FT?sIGZrRL>2}C|`C$8PWtC6UjO8?6?Z)@W(uRbUD8I?4Rn^JAPBijcQe;H?W4W?Ns=a<Pn&Em2Py>JPnSCVYUI8O9@*#(PGt$a0lR}Hhd{*OwY3yZEkG(XSf#MZ9eA+Zj=K-pSr+0t<51rwdbvB~QaRx3B8>3>+%eUhw6B&cgg?FuZytO>iL_l}>onrp1AIz%MY=}xf|sg?(@xV*>xMYfOhlhycCKCqVjS+_iDw!n}gTOunt0nPGrw+OBw*Of!*$t&zng<f=%mFK(^rl;3=NC4W&tjc-J8p!jrr<6pwSCsicfA}K}&Ht^T2UfoGhKh2K?78p;zF7aAuxy^|T1`f0=f0L#+aJ8}X=YYkk@|B|y#Nc$;l!)%C#z;%i7{THJr}k1nJoAHaeCY-yglzX5kCW=>cez)Rlr0Cm94>@T>s=orDeu!IY-^7PAVo=0afY>SRE;<{GqU4$MPQ&ojd&LJ~*y^DsQq^S8KBN_8%%xe|_3VBk!aC`WHFfcoq')
PACKED_EXTENSIONS='c-n<jZEou@47*DHbp!H4-?q<E^#0d)BD0~*n(d*0YbdrzQIzi=_nXLV-)<`Saog)--)!IR>aP3d5p&;CJ@+kUcErZ$&z|0Kw;a7ndVj54`&OocNJZb0VGsPt=e}Kn=?h$}PU7W#JzWb4y;qw~IHt8K!pRxVf>fJp$J2UWnsS2HzbXx+wTgu#xE`@OZOQCNj5$C<iZq-wP`f~xzN(}Fk`^FuHs7lj)vt>6RJEZ-9}qAdk{T=ReaDAKAL@LGGpNoL5%&QCE~Md$Gz>IkT?sCi&=C+Vr-HLXI%5k~b-g)Ik$DnYAiic1+6lehpyv<4^(p$B;D+Qm!M|uHbr{JPWn4_oZhmzc;zt;0tbkU-sKzLmbLel--?G05ZV3+yYA@_t(gRuveRMVnJqeBTfivka-f$(oT6DNv^!iP21+vUzZa`l&&in2h@F&$(A)|RTdAum?hvl#Q^-Gs^a;#8AvQiruu^yHINg;cQ!g#vCX$sq6>@^vmZRQs?)4wOX@EF*r<H3Nm2@0hUG7|l_2{}C`DD(@%Ce+hhG-J|?2E#@s;#%W>raKYgV=?4V1bS<+!g+Yw2<ty=XPEh*Cvw8~J)-)3J~9p>cCE;DQwVM&<XI0yE?H=b8FUf6X*tjo3(iRKC8X#3d7wo|<a!Jr{{Wh$7%K'
def gain_extensions():
 def message(bits):
  x=1
  for b in reversed(bits):x=int(x!=-1)if b=='1' else-int(x!=1)
  return x
 def regions(prefix):
  if len(prefix)==8:
   return{prefix:(message(prefix[:6]),message(prefix))}
  a,b=regions(prefix+'0'),regions(prefix+'1')
  if len(a)==len(b)==1:
   x,y=next(iter(a.values())),next(iter(b.values()))
   if x==y:return{prefix:x}
  return a|b
 parts=regions('')
 total=0
 seen=set()
 for term,multiplicity,scales,center in unpack_literal(PACKED_EXTENSIONS):
  p,q,r=term
  key=(*sorted((p,q)),r)
  need(key not in seen and multiplicity==1+(p!=q),'fresh ray classes')
  seen.add(key)
  (a,A),(b,B),(c,C)=(parts[x]for x in term)
  need(A+B<0 and C<0 and(a+b>=0 or c>=0),'fresh ray cylinder')
  total+=multiplicity*lower(tuple(term),scales,center)
 return total
def compute(orthant=None):
 budget=controls();seen={};one_side=0
 for line in ROWS.strip().splitlines():
  p,q,c,k,z,*v=line.split();key=(*sorted((p,q)),c);need(key not in seen and int(k)==budget.get(key),'disjoint cylinder multiplicity')
  seen[key]=int(k)
  rays=(p,q,c);sc=list(map(int,v))
  chart={('00', '0', '00'):(7, (-1231, -475, 585, 3139)),('00', '100', '00'):(5, (-1548, -437, 1135)),('011', '00', '00'):(5, (-3323, -619, 871)),('00', '0', '0100'):(5, (-1563, -403, 1190)),('0100', '0', '00'):(5, (-1509, -403, 1188))}.get(rays)
  if orthant and chart:
   K,d=precision(coords(rays))
   signs=[1]+[2*int(b)-1 for w in rays for b in w]
   _,nodes=chart
   w=iv(orthant(K,d,signs,sc,T,PI,tuple(F(x,1024)for x in nodes),hermite_expectation))[0]
  else:
   w=lower(rays,sc,int(z))
  one_side+=int(k)*w
 need(seen==budget,'all cylinder classes represented')
 one_side+=gain_extensions();overlap=overlaps();gain=2*one_side-overlap
 need(gain*10**6>1287*S,'shifted ray-gain threshold missed')
 return{'method':'Directed halfnormal interpolation on disjoint cylinders','threshold':'-293/800','scale':S,'one_ray_lower':one_side,'overlap_upper':overlap,'gain_lower':gain,'cylinders':sum(seen.values()),'classes':len(seen)}
def centered_moments(B,ell):
 """Directed variance and third centered moment of the tilted Gaussian log density."""
 n=len(ell)
 di=[B[i][i]for i in range(n)];rs=[plus(r)for r in B]
 bs=[add(mul(scalar(mu,2),rr),l)for rr,l in zip(rs,ell)]
 q=[(neg(mul(d,v2)),b,d)for d,b in zip(di,bs)]
 hz=[plus(mul(x,mom[j+1])for j,x in enumerate(p))for p in q]
 hzz=[plus(mul(x,mom[j+2])for j,x in enumerate(p))for p in q]
 var=plus(mul(mul(p[i],p[j]),mom[i+j])for p in q for i in range(3)for j in range(3))
 third=plus(mul(mul(mul(p[i],p[j]),p[k]),mom[i+j+k])for p in q for i in range(3)for j in range(3)for k in range(3))
 ed=[[scalar(x,2)for x in row]for row in B];v22=mul(v2,v2);v23=mul(v22,v2);cm32=mul(mom[3],mom[3])
 for i in range(n):
  for j in range(i):
   e=ed[i][j];e2=mul(e,e);var=add(var,mul(e2,v22))
   third=plus((third,scalar(mul(e,mul(hz[i],hz[j])),6),scalar(mul(e2,mul(v2,add(hzz[i],hzz[j]))),3),mul(mul(e2,e),cm32),
   scalar(mul(mul(e,v23),plus(mul(ed[i][k],ed[j][k])for k in range(j))),6)))
 return var,third
def unpack_adjacent(payload,quantum=125):
 from base64 import b85decode
 from lzma import decompress
 data=decompress(b85decode(payload))
 witnesses=[]
 while data:
  need(len(data)>=9 and 0 not in data[:9],'packed repair prefixes')
  words=[bin(x)[3:]for x in data[:9]]
  n=sum(map(len,words))
  start=9+2*n
  end=start+n*(n+1)//2
  need(end<=len(data),'packed repair length')
  entries=[(i,j,quantum*value)for i in range(n)for j in range(i+1)
  if(value:=data[start+i*(i+1)//2+j])]
  witnesses.append((words[0],words[1:],list(data[9:9+n]),
  [quantum*x for x in data[9+n:start]],entries))
  data=data[end:]
 return witnesses
def adjacent_floor(first=0,last=184,*,bridge=False):
 from lzma import decompress
 L=6
 def source(L):
  edges=[(0,1,0),(1,2,0)]
  for ray,anchor in enumerate((0,2,0,0,1,1,2,2)):
   tag=0 if ray<2 else anchor+1
   for node in range(3+L*ray,3+L*(ray+1)):
    edges.append((anchor,node,tag))
    anchor=node
  return marked_gram(3+8*L,edges)
 def sg(x):return(x>0)-(x<0)
 def out(x,b):return sg(x+2*b-1)
 @lru_cache(None)
 def bounds(p):
  a,b=-1,1
  for x in p[::-1]:a,b=out(a,int(x)),out(b,int(x))
  return a,b
 @lru_cache(None)
 def cbit(x,y,b):return int(x+out(y,b)>0),int(y+out(x,b)>0)
 def mbit(b,c,d):return b*int((c==2)+(d==2)>(c==0)+(d==0))
 def literal(core,p):
  pins=list(map(int,core));n=len(pins)
  trits=[range(a,z+1)for a,z in map(bounds,p)]
  partner=[{sum(cbit(a,b,pin))for a,b in product(*trits[2+2*j:4+2*j])}
  for j,pin in enumerate(pins)]
  for a,b in product(*trits[:2]):
   l=[a];r=[b]
   for u,v in zip(pins[:-1],pins[:0:-1]):
    l.append(out(l[-1],u));r.append(out(r[-1],v))
   e=[cbit(a,b,d)for a,b,d in zip(l,reversed(r),pins)]
   need(all(a[1]and b[0]for a,b in zip(e,e[1:])),'positive adjacent C edges')
   d=[{c+mbit(pin,c,pp)for pp in other}
   for pin,c,other in zip(pins,map(sum,e),partner)]
   for j,(u,v)in enumerate(zip(d,d[1:])):
    if n>3 and 0<j<n-2:
     need({a+b for a,b in product(u,v)}=={5},'run gain')
    else:
     need(max(u)+max(v)<=4,'nonnegative gain')
  return prod(map(len,trits))
 if bridge=='matching':
  rows=[tuple(v.split())for v in(
  '00 . 00 0','1 00 00 0','00 . 01 00','1 00 01 00','00 . 00 100','00 . 100 00')]
  seen=[];mass=0
  for row in rows:
   words=('00','1','1','00')+tuple('' if x=='.' else x for x in row)
   need(all(any(not(x.startswith(y)or y.startswith(x))for x,y in zip(words,w))
   for w in seen),'disjoint matching cylinders')
   seen.append(words)
   for u,x,v,y,l,r,a,b in product(*(range(m,n+1)for m,n in map(bounds,words))):
    A=cbit(out(u,1),x,1);B=cbit(out(v,1),y,1)
    C=cbit(u,out(x,1),1);D=cbit(v,out(y,1),1)
    need(A==(1,1)and B==(1,0)and C==(0,1)and D==(1,1)
    and mbit(1,2,1)==1 and mbit(1,1,sum(cbit(l,r,1)))==0
    and mbit(1,2,sum(cbit(a,b,1)))==0,'matching bridge')
   edges=[(0,1,0),(0,2,1)];n=3
   for word,anchor,tag in zip(words,(1,0,2,0,1,1,2,2),(0,0,1,1,2,2,3,3)):
    for _ in word:
     edges.append((anchor,n,tag));anchor=n;n+=1
   signs=[2*int(b)-1 for b in '111'+''.join(words)]
   K,det=inverse(marked_gram(n,edges))
   mass+=lower_motif(K,det,signs,[7000]*n,0)
  return F(4*mass,S)
 if bridge:
  rows= """1011 1 00 0 00 . . . 11 . 00
1101 00 1 . 00 . 11 . . 0 00
1111 00 00 . . 0 00 11 . 00 .
1111 00 00 00 . 11 . 0 00 . .
0111 11 00 . . 0 00 11 . 00 .
1110 00 11 00 . 11 . 0 00 . .
10101 1 1 00 0 . . 11 . . . 0 00
11101 00 1 . . 00 0 11 . . . 00 0
10111 1 00 0 00 . . 11 . 00 0 . .
01101 11 1 . . 0 00 11 . . . 0 00
10110 1 11 0 00 . . 11 . 0 00 . .
110101 00 1 00 . 11 . . . 11 . . . 0 00
101011 1 00 0 00 . . 11 . . . 11 . 00 .
111011 00 00 . . 00 0 11 . . . 11 . 00 .
110111 00 00 00 . 11 . . . 11 . 00 0 . .
011011 11 00 . . 0 00 11 . . . 11 . 00 .
110110 00 11 00 . 11 . . . 11 . 0 00 . ."""
  seen=[]
  mass=0
  for line in rows.splitlines():
   core,*rays=line.split()
   rays=['' if w=='.' else w for w in rays]
   k=len(core)
   need(len(rays)==2+2*k,'run ray count')
   for old_core,old_rays in seen:
    need(core!=old_core or any(not(x.startswith(y)or y.startswith(x))
    for x,y in zip(rays,old_rays)),'disjoint same-core runs')
   seen.append((core,rays))
   literal(core,rays)
   edges=[(i,i+1,0)for i in range(k-1)]
   parents=[0,k-1]+[j for j in range(k)for _ in range(2)]
   for i,(anchor,word)in enumerate(zip(parents,rays)):
    tag=0 if i<2 else anchor+1
    for _ in word:
     edges.append((anchor,k,tag))
     anchor=k
     k+=1
   signs=[2*int(b)-1 for b in core+''.join(rays)]
   K,det=inverse(marked_gram(k,edges))
   mass+=lower_motif(K,det,signs,[7000]*k,0)
  return F(mass,S)
 def coords(core,p):
  return list(range(3))+[3+i*L+j for i,w in enumerate(p)for j in range(len(w))],[2*int(b)-1 for b in core+''.join(p)]
 def conflict(a,b):
  c,p=a;d,q=b
  return c!=d or any(not(u.startswith(v)or v.startswith(u))for u,v in zip(p,q))
 def lower_shift(sign,mat,shift,T,K,det):
  n=len(K)
  B,ell,D,jacobian,gram,linear,scale=tilted_products(K,sign,shift,mat,T)
  tr=F(sum(row[i]for i,row in enumerate(gram)),scale)
  q=F(sum(map(sum,gram)),scale)
  lin=F(sum(linear),2*scale)
  p0,p1=PI;ml=sqrtf(2/p1)[0];mh=sqrtf(2/p0)[1]
  E=F(n,2)-tr/2+min((tr-q)/p0,(tr-q)/p1)+min(-lin*ml,-lin*mh)-D/2
  var,third=centered_moments(B,ell)
  lower=two_node_lower(E,var,third)
  return jacobian*lower/(2**n*sqrtf(det)[1])
 w=json.loads(decompress(b85decode('{Wp48S^xk9=GL@E0stWa761SMbT8$j;4CEwQC$En;DlqKA{Q>o^b~Tb{J;@G5w4trqqM5Hy@#3{f0gr>t#gUM5HO=Ned8}y!v<4~tKEA}LU0KOv?Q$^NQvG3heTBN)A>!9G$UfYcjo8klH^mBo1kZBQj!eAmsjPS1_Wu>_b<dn5fHg447!?ij)JiMyQn-(m>JL7!%Hc1jGu#j__7{my0w|Iq^cDYa<O@Z)_-j@jcNh$!!H$+r2p&8{y@-4dJ}R4I~5y@`;;l41wC>nAioLXDOrz_HWTbCIWb}<u}?mBve1#5hy+4hJ-){GQUr6p!n1#XgL5lo&S*fcZv}@SFUDhb2ASwY>Cnxm5TNmbBxF8v(Exj~inT92w1z*)-B1K~Qs7A@qTIZj-yz@QhFZACnwo19+|Z1(>YiKI+@5y^Hu_h{Ao(siX95XHGpJu|Gyr>Z*_X}9VJ+>RVSicXyqV`{^KpAbO&9F+-WfE{BDvgZp5<@__MU7okl^MiUuFYxac_;H&1tK^1f+!4{du~HI*uR{*`~Uk;Lc7e(Yk2Hd55H#q&V<g`woeKe;)@wQ%f}>CBb|kEhBU&e|X_3NsLW&DV0JPi;r1`Sjy5*A6+G%NO(vTZUe%vL*k^%nMZXO*A2~L=m%{@Qv?eJY)E?@TRy^u1}vW)dSnc+-7eC;k$B2V`#=Z(W=RY`B66Kl$Oiky*G2ZDe`YJ{$iW)Zwd03ko92Ff!^r~WTOn4Q*wb0L>KiR;Y&lhb(~GqwA(1oYxg#na-)ag$XB})%7`yW%u(ZdgDKx38`joZ^>jM8m-k?AXI&Qj+O0DGprrpr&{U_KQm4xV#Wr}lo&hI9^toc=4Gg<g-V+qp%N7%}w{{8f7n1uDtOSr{;kr9oZG~_?v1W6a^tlPmC4R45_O?of5N!^E>YS09Tid$JQgdC2-0zEM$+aD|o10d*VgVBz6{|>8;<FA`bmaf4XLJ;CL>h)FNIo_dHhQ!yd*P7E=v~=wgx5B-?gvT@1piY(;%6})^Ht$WL`9g9yHRg$p$dCUJ3yl3@ponEMajxg4Si=5@)^O*7erIrwFr1ZveGNgDf5uPwdn~KgUt5#q;ddhJ%Sqe)sS$<E_pGWIrNt)neBd!B8lwgiH(3`d#Yj6!yvO<wE>Ypo-1th8_97PEB9IhtE1TO0Px-+iS=eukXHfu}Np~<J3coy|9<0#hlI<rH{3q!;;B)SCIIG3I@T(4GA@-;!eEy*(MmJSZ?LCb<ys(jS`B}_?Q|<Je%Z7#AnX*0#qzD?t!=ZA*)tw-mv1IBEq&|?kC#1USV-_P*tddh+ae<+ClSCSoF@k8U<l^8llbp|-fxu_S2F!p4g8QUVE!Mthr~dN+NN1eri$t|Vkv+oDUoshvMzc&^bQCcvmOVW^LftTKj>k{rE?*H~cKy2CI7^*fk64?BH)Qn{T$#{`tYaw?_nacO`uc|gnaj5%_3DKigz!1Fg49~YcV1Smh}9vOyCX0crn(N1YZ8FEfrpO|u?eY7jRMg?RPb;UMZK*8@sxG8??)AeRH+=W<=vbsJ;8#S{zlM^XP4J7$!-vD{ed_6TOL|KAkn|2Ts9C?Ye3^XjSh52e7@WOT{FuInFO&j<V|l`t7|VnVw~I4(c@o;7KF=^xxEzpwqm|U2A+wy;WNBBsOL^l%5A4iw?z^`?=>vMGSKAv;mG2Jl-xXzvlZ_dwb{r*{#`T#BX20d?vta-{szZOzOypixzI^7DN}@Y*fQaQtM&Bsuwc4mAOO0L(YaI_8D9t^$gNqAh5sDE{;4AvSI$nnOSPO$i9b6C$RmGL<JrF5OS;Jk)E_skVtlzx@?JmQyCm`%>Ux&inW2am(JGnlB4HcWJo>km*F&^P?S0@6<`w4f^1#CwGU!_XW=^EKVO+;17YD4(D+!hcC23czM&;5bKWuZVI3l1-+W8gq&*NV0i7{m_|MTsOH5Aqilx)z;rKZyJE3e~U<7Xa>y<t8vGP_!R^XFSprIDz!-Eot4+@51{-UKI&dk|3Ul1z2CNUTZ2hP6#@n0Lkse;?O7-NoT$!)Is=i4Xq=pQT<8{{XEdm2}%LdaDJ^L>OPYM^YUyoa*|nr<lP&=O<+1vBC7aSzPEIKs_Du43kP)s(smd1cYjKOczo2^?n8oL%zc7EKqdrkZETh(m5+1kt@;r_s)%%lC-<PFWXcya#r+aMbS%Qn2T%nCb|J}*JpT^wfExf5Z<;%Nk6=I%<`6{8eKlrW3t+tCB7)CG;Wb--3Ynv1t;zM8?hlR4@h3l^Vt*wo^KNAy?Goe#TlS+YTGf+s(BGa6+56AZ5orT`%O}`!Msm<m(}_S7zgOc4uSbgi~r;Ya8`x&@#&U+Mf9P*yPvOHP!p5`2&u?s_#`Kur?q+7m!Fo0-%*FxnZBmHQLll(n$|#LgouZqr|||_o7hNT5d$G<<joFX)S)~S0zkXLcS^$bV|*7_$*={8M*hp5&xP&>M}h%we~7q=I;hF535`sxu+fmAmp+SCcAhj-rt2AD@N`sxbM2HAdi=r=N-D!|HOR#GZ>#J*F>GIit>YoiMQG%ex=vKiORYT4T+UeNgrAs(N;CTUc~YIWdMes{hC|^IW*e{Is*nbZ`RinMxA8dy7u?p1B{nOnD{-@#I<qdEMJ;!M#=TIij4)nn;)i1xKT;J5e^o>O_;^5AjjNOhu&sN`$gM!Us0)lA`gUKAr?mTkGzjhL7!ZZA`hgBzDbR`uqk1O{;Dm>vFl<PS-S@x%00000XALQxqjt<>00Hd~rdR*~(sIB9vBYQl0ssI200dcD')))
 w+=unpack_adjacent('{Wp48S^xk9=GL@E0stWa8~^|S5YJf5;2B{Em|Xw|fk_n@Hymh}<V5|hP4P0A(i!4_(Nm=uw{@FVR6daNE=O{}nmnEO!yZZjleQI$`fN#h?Zw?7ib$*yf+)4nT1c@YHQsf~P3NwFbs`e6kY9b7XM>tTnD58W=8+2z&W~&N0PZ_L7MnK&>FsOZ7Gcz2GM?9?8!S_Bs2h*&<7f6Md!?c(6uCbSV0$IFiJDi57cLwUCP&bOl(nS#JyQM<#7?AVRUr5_S4dtQ5L8J6<1dz+o>()+v+iXrbxTNw<FkmwfG$RA`{Asc<Zy`^ovcJFpXlbIoj{E;@jPVs4G?E4*~DY4(9#!@$mILXJfIm9L5g(X&R@<TbX|&vL@@X$kD@x}F*cU~5NeYkW#f<uHjLEVuL<2P%gF07h3~;0K49$&mWb17+fU7s&DT>$xpH01FCyQ^zc-~(`aYG?FX`9$%A`XjnT@ipIaYecf6E)o<0_q#Ctbr1Yf~{p(T?u*aX05A8$(*(2MR)|UXVryHPT(_+}BgQQQocG;lhy%e)eK5d3R$j$BQPqXj${n|C0H(=*D$H$c5uV3u=c{?l-*g?>OAhw_ny5<PPLm)$>U^jB=c_!C?@o5Qocg^ORc7t!|W|4n7ND(NsT6HCGsY^}jD6N8b0VSa55gzb0GCm2Q-mK7CCy{o3_mdRH2KtnW79`h>gSmHB*KL9aI+9Ef{y{sq07kAc*7-K+n)EPTu5>&-UkAUTzw)X5sDE%I77heGDn0qw5$%dboi5@N>CF?%|`Pf>w5rC}}B(D$dIRB=2#ULGP(?*=&tE3)Vsgw)Sy>raQB!flaa$~BK3w&kZ43qH6P1VVG5ysii_8;Mc+zLZ`B9O*|TErC~+l}e_d68}Agkknj1ocpi!G%@nd*Uu1F20p2!nNIlV&VO+5Z}4qbHVARtr0q0q2w85^4@8%^oOE1F?m|l=F`f3=X#o~ID-_d@azFaTt2RT;LHb~nbWo$cF=;y(%>1fgXR&NRR`{^L#Pidm8k4DUDzy7IpX)n0F*9syy*(;{{5~?4kMxtU0ZTr@PkseNSNc~wsUhnoq%0(E<~w+#brSHuv`a<%1J2V38IR&;7DU4ifJrN9dbW~bDNR`Z?qjGILp!|msb;HPu2EWrrl*=>Svyvju;W%lvBm5#&Pu`E#q~`(1h5y_V~}*(gn1e*<3<umNrGRaVhpkZ143Lq(d^GmI;tZM%waZ2#R0iSFE7kjq80^p5b(ld51o6iCH<qa2*gpx*dfUu*OK>?u^X=jTMjqiI&3KZ5~UEiw{pa-`|YV%7THq38ip!0EzH#nD#<H>R*J?Q*{V0?{7ki%!|g@xQoS%R)_hT;yAE|&U7eSc#yhaDlZEZZY%?RDo2D%qNhgj3LF)+o+2J<W*I*!d;##l-&9nlR&t>t3N|FOP?YCu=#km#TMVMwtrc`?S7F(ud+8_EWE|ehCUTWN4#8nCs0YboMdqStCYbu>PFUo6L&k(yl*y(1&8X95~(H{-twWVk^p$I&Anfzgh;}{Yke0YHXebO8s%1d6Z(7Iz&9t|x^_n5nNm@(d%RK_TYt4RKqTxiJrd~<KSy3>24yC*KXp%9=?PQp<jEZ+;i>Z`4P$)hW0mv9%RC9%c<><JI=GgXcC6_c)x5m4QwV`N`5v1{P%J6B>YC#&8N!yoI>zEk7;7nOOOlKVc9Kk7wc?rIj|%@f0-Y$;1ODr$tdiBt&8?~3AyQfa?K@`NR%8inb(C&ybmZbE&YI^#}`$C}QxMk8Q*+6_#FR1-17aCxijG2WvUD5M6`ZsLI!VVV3c6{<Z^4ol&GYV;F%i?O+I?3ElyLJs{|V8pCEGG%)^cjd$$kV*D+b&41F(G9`qDyGlq^TDp?S}j^lua#B6Vuz4M{+rmkzt9w!Ol5O&<NoaQu$^08?!`!%8|tWvhuM^0^A1u3>y3?_TV7xIJYGlL+rPf6{WOzETzb&Rqj0ug0f~7`j%}er{RIdLdpmSL<Ej08H}Q$qHw33E<AAHij603Zi6Q<Qs>u1NDwf_HV%8=p^YE15;N!<>!m?Fa|A3(UImhcM0w`A{GQIXQVJnVD17_8s9|_OgcL&NXIMV!3-W+6!F)No<Ov(E;p;)ec)fk}LRTto~x`zR`^gCsF|NEYxZ0DFZq(ataBE;aCV-|FnoHz|Ll1YmMFO($SWtNHfi^^Brd>brCI+F5Q6CpDu9831a^dt;<!HWA=n@o5&w~^+1+(K8ps0G?rd5j=h5sbS8F0bLDN4gu=i3+Q(b>G8n>}4K3ZHS0!<VeN?i1Lo9AY4}B$tNUzRhwz2hoQt8u1MowIUv)>MB+*fc@~FwF?!hw7lfe$t%&$BLm?LH1+=42^f83s%8jgbE}t9FUNfUmGjxjIPaaKt5hO7Este3mq@_3&y*Tg_UZsPz|K;L6EuM*9u}4ofPLU|Y&gYa9=2NaByMe5+-BN?2d|^L${dR5y23F=f;382zeKk|d!blvyvDH!+842!d+r=@R-@||Zxj23j)8m^V6c97;R5?~(?7nK2P_&`#tH1NlFLh(~Cq?_MHiXuHy3an*RfXuMtL(Jn`H~%%@6x{JoS!$ui4in(`-#WE+%y|pHh6b%3J-mu1Ad^i(rQ<x0hPAvda<+G9wm~$WL&b77UlX@EZ4i)YG~>D8=J1wEkx4#zX34?R}=lTFs@k*4n-ZGGoQAW(Kwa4)JP9OU=n%S`|%5VBl4@n#d@IoR$F(-_X$ir;nwUC#cba>_xx4KIp5>oq55pJeWyF(H7mILUFPk{MY=?Mfn(pZh_UWuEr41euKB{Fp{HzZcuZhkfd`jZYg;Ei>;@sPnY;i1^fo<+B2#sl00FcS;xYgL!Gjk>vBYQl0ssI200dcD')
 raw=inflate(b85decode('c-nPVNs`+z5EIGfzG<qM(|p5E@&Dfvpd}|6Njdn^Bzr*vkf(?fXYdPNTF)tN#?~@rIWa~|6<<u)Tua^NIGt3nIPlTrYB%YjB-f^MR>WGdzHgz%oQg9|kOUbVaz`({<TQ{w)HX#BHg2(I95r=gx02gjM;!Nq|J)bccc>}e7=)!FjxdU%Z6n@FV~9u2nhm)PN;~JyyPT)>8ka#X<2D}?!X}8Oz~y>gj<w)@POXm>P2`mQXX?d*nDUCr;_x?20POkzJ8RKlV(#ZO?UtS6SXmasW2Brq?t3i4aLOh>Z+GX~(r*0tgLqfS7atfxBtFEHOT|Q;8$$d=40|QCI4|vtzl3qXlr^3T;<#(L(VzJdZ;C;;0V^qO#5^~baL|+jAs(HMiDM)?xq*>xFMxta@o@oy6`tY*$m7i=USi}D4T@ofxD<R&6ENsxNL9J3p32~|-Xb4-Bsl-UXUbcia6>bY1qmNy=vJb2n<=0wZa4SXB-_8m<^9Q@iYLSKQJ0)*xBE3RUaf2~fRU9%cO@Teux-_0V)lQ@wy3mIIZQMM=G#(NFiEy6`h%LJrXKs3HZKhgiAs~}eQtL>@Bsv>i~E8iTtGvn^_T{>K_#LiDzP`Dq7qFVsuyF{Z38Mh@JHsID63I;54BfgZAfx9N-tn(k~?ZFHJ)pYUB@HiuWKOQ*yrZMalmGgXUTJ(8V@xYJYcHhL1ssuz}9*0ZD4QRHLT$vyN#o>01I#CFt-JByRN}Sw(}qwXV;T1l}J?07LEc~yMrTs+;h=x`e{!-HywaLTt)+3mSTM&p8C6$my@XJ=MJ^`W{1l1rT+S_TOP!QpMlr|#7-o(K@E-D<;R2RuKx`t?URVwKq6n0n!4C}m^xCy7%;)<XjnnZW$$r#lP2bH;XHvuJ)Ja8aXUHfu7Omgec>B6akPJqq|xD@j_`(<ME%$T^hr0VdGd3~g!_Cllot>Nz(bYYgWcE%)Nc9c;fWBk;4a>$h@r99qqgN2!1>{;-;7%Ns$f5kO(!<G4x9I9JG+o76rF{BVNF-8#`-=IIji3oXad<YEmQ^H>8+-1uBNivJR{YjAD^d%)c5~R^h0(puXD5;5l`o8oa8Wx{Yxk5RBL-HhpxE}Klr7#te;E413Wd~R=*P0O}yptYdm-02)(wag|GLC*5`{nY2^jbFM09c&MOry>if%xqJMg-yK3h*VuwfK'))
 extra=0
 while raw:
  core,*prefixes=[bin(mask)[3:]for mask in raw[:9]]
  n=len(core)+sum(map(len,prefixes))
  need(len(core)==3 and len(prefixes)==8 and len(raw)>=9+n,'repair row packing')
  w.append((core,prefixes,list(raw[9:9+n]),[0]*n,[(j,j,5000)for j in range(n)]))
  raw=raw[9+n:]
  extra+=1
 need(extra==105 and len(w)==184 and 0<=first<last<=184,'witness counts')
 C=source(L);N=0;s=0
 for i,(core,p,perm,b,entries)in enumerate(w[first:last],first):
  need(len(core)==3 and len(p)==8 and all(set(x)<=set('01')and len(x)<=L for x in[core]+p),'prefix shape')
  N+=literal(core,p)
  need(all(conflict((core,p),(row[0],row[1]))for row in w[:i]),'disjointness')
  ii,ss=coords(core,p);n=len(ii);need(sorted(perm)==list(range(n)),'permutation')
  ii=[ii[j]for j in perm];sp=[ss[j]for j in perm];cp=[[C[x][y]for y in ii]for x in ii];m=[[0]*n for _ in range(n)]
  for a,k,v in entries:
   need(0<=k<=a<n and type(v)is int and v>0 and m[a][k]==0,'sparse entry');m[a][k]=v
  from math import pi,sqrt
  K,det=inverse(cp)
  P=[[float(v)*sp[j]*sp[k]for k,v in enumerate(row)]for j,row in enumerate(K)]
  alpha=sqrt(2/pi)
  scale=1-alpha*alpha
  L1=[[v/10000 for v in row]for row in m]
  b1=[v/10000 for v in b]
  mean=[v+alpha*sum(row)for v,row in zip(b1,L1)]
  g=[sum(v*u for v,u in zip(row,mean))+sp[j]*float(T*sum(K[j]))for j,row in enumerate(P)]
  PL=[[sum(a*b for a,b in zip(row,col))for row in P]for col in zip(*L1)]
  active=[tuple(k for k in range(j+1)if i>=79 or m[j][k])for j in range(n)]
  for _ in range(32):
   for j,row in enumerate(P):
    diag=row[j]
    for k in active[j]:
     v=L1[j][k]
     column=PL[k]
     t=alpha*g[j]+scale*column[j]
     if k==j:
      z=t-diag*v
      u=2/(sqrt(z*z+4*diag)+z)
     else:
      u=max(0,v-t/diag)
     L1[j][k]=u
     delta=u-v
     if delta:
      for l,entry in enumerate(row):
       change=delta*entry
       g[l]+=alpha*change
       column[l]+=change
    u=max(0,b1[j]-g[j]/diag)
    if u!=b1[j]:
     delta=u-b1[j]
     for l,entry in enumerate(row):g[l]+=delta*entry
     b1[j]=u
  m=[[round(10000*v)for v in row]for row in L1]
  b=[round(10000*v)for v in b1]
  v=lower_shift(sp,m,b,T,K,det);s+=v.numerator*10**15//v.denominator
 return F(s,10**15),N
def semantic(part,cap_numerator):
 from bisect import bisect_left,bisect_right
 from collections import defaultdict
 from itertools import combinations,product
 from math import lcm
 S=1<<74
 def ceil(x,y):
  return-((-x)//y)
 Q=1<<80
 X_DEN=1<<75
 def ln_interval(num,den):
  """Bounds for Q*ln(num/den), where 1 <= num/den <= 2."""
  a,b=num-den,num+den
  aa,bb=a*a,b*b
  ax,bx=a,b
  lower=0
  upper=0
  for j in range(1,24,2):
   q,r=divmod(2*Q*ax,j*bx)
   lower+=q
   upper+=q+bool(r)
   ax*=aa
   bx*=bb
  return lower,upper+ceil(2*Q*ax*bb,25*bx*(bb-aa))
 ln2_lo,ln2_hi=ln_interval(2,1)

 inv_ln2_lo=(Q*Q)//ln2_hi
 inv_ln2_hi=ceil(Q*Q,ln2_lo)
 @lru_cache(maxsize=None)
 def phi_units(num):
  if num==0:
   return 0,0,0,0
  need(0<num<=X_DEN,'entropy argument')
  shift=num.bit_length()-1
  tz=(num&-num).bit_length()-1
  ln_lo,ln_hi=ln_interval(num>>tz,1<<(shift-tz))
  e=(shift-75)*Q
  log_lo=e+(Q*ln_lo)//ln2_hi
  log_hi=e+ceil(Q*ln_hi,ln2_lo)
  return(-S*num*log_hi//(Q*X_DEN),
  ceil(-S*num*log_lo,Q*X_DEN),
  S*(-log_hi-inv_ln2_hi)//Q,
  ceil(S*(-log_lo-inv_ln2_lo),Q))
 need(part in(0,1),'entropy phase')
 PART_ZERO={(0,0),(0,2),(1,0),(1,2),(2,0)}
 split_pairs={}
 mark_cells={}
 @lru_cache(maxsize=None)
 def binary_units(num):
  a,b,c,d=phi_units(num)
  p,q,r,s=phi_units(X_DEN-num)
  return a+p,b+q,c-s,d-r
 @lru_cache(maxsize=None)
 def entropy_support(masses,paired):
  units=binary_units if paired else phi_units
  lo,hi=min(masses),max(masses)
  if lo==hi:
   lower,upper,*_=units(2*lo)
   size=len(masses)
   return(lower,)*size,(upper,)*size
  exact_values=tuple(units(2*x)[0]for x in masses)
  origin=max(1,min(X_DEN-1,(lo+hi+2*sum(masses)//len(masses))//2))
  _,upper,neg,pos=units(origin)
  tangent=[]
  for x in masses:
   dx=2*x-origin
   tangent.append(upper+ceil(dx*(pos if dx>=0 else neg),1<<75))
  return exact_values,tuple(tangent)
 def clip_curve(entries):
  cuts,intercepts,slopes=[],[0],[sum(weight for _,weight in entries)]
  for cap,weight in sorted(entries):
   cuts.append(cap)
   intercepts.append(intercepts[-1]+weight*cap)
   slopes.append(slopes[-1]-weight)
  return tuple(cuts),tuple(intercepts),tuple(slopes)
 def grouped_clips(terms):
  groups=defaultdict(lambda:([],[]))
  for index,weight,cap in terms:
   groups[index][weight<0].append((cap,weight))
  return[(index,*map(clip_curve,entries))for index,entries in groups.items()]
 @lru_cache(maxsize=12000)
 def clip_value(masses,positive,negative):
  whole=1<<74
  center_twice=min(masses)+max(masses)
  folded_center=min(center_twice,2*whole-center_twice)
  cuts,intercepts,slopes=negative
  j=bisect_left(cuts,(folded_center+1)//2)
  def entry(mass):
   r=min(mass,whole-mass)
   a=bisect_right(positive[0],r)
   oriented=mass if center_twice<=whole else whole-mass
   return positive[1][a]+positive[2][a]*r+intercepts[j]+slopes[j]*oriented
  return tuple(entry(mass)for mass in masses)
 def prove_simplex(linear,constant,terms,biases=(),symmetric=False,region=-1,interval=(0,8),signed=(),credit=(0,False,F(0)),marked_rate=F(0),marked_sides=0,fixed=(),lanes=()):
  dim=len(linear)
  grouped=defaultdict(F)
  for key,weight,*cap in biases:
   grouped[min(key,tuple(4-x for x in key)),cap[0]if cap else 32]+=weight
  biases=[(key,weight,cap)for(key,cap),weight in grouped.items()if weight]
  bonus_count,deep_bonus,coupon=credit
  common=lcm(constant.denominator,coupon.denominator,marked_rate.denominator,*(x.denominator for x in linear),
  *(x.denominator for row in signed for x in row),
  *(x.denominator for x in terms.values()),*(x.denominator for _,x,_ in biases))
  c=S*int(constant*common)
  line=[int(x*common)for x in linear]
  keys=tuple(dict.fromkeys((*terms,*(k for k,_,_ in biases))))
  first=dim+1+len(signed)
  indices={k:first+i for i,k in enumerate(keys)}
  paired=[]
  for key in tuple(terms):
   other=tuple(4-x for x in key)
   if key<other and terms[key]*terms.get(other,0)>0:
    u,v=terms[key],terms[other]
    weight=min(u,v)if u>0 else max(u,v)
    terms[key]-=weight
    terms[other]-=weight
    paired.append((indices[key],int(weight*common),True))
  terms=paired+[(indices[k],int(v*common),False)for k,v in terms.items()if v]
  biases=grouped_clips((indices[k],int(v*common),cap<<68)for k,v,cap in biases)
  n=1<<72
  signed=[tuple(int(x*common)for x in row)for row in signed]
  vertices=tuple(tuple(n*int(i==j)for i in range(dim))+(n*line[j],)
  +tuple(n*row[j]for row in signed)+tuple(int(n*k[j])for k in keys)for j in range(dim))
  unit=int(common*coupon)
  munit=int(common*marked_rate)
  marks=mark_cells.setdefault((dim,symmetric,region,interval,fixed),{}) if munit else None
  def linear_lower(v):
   signs=v[dim+1:first]
   if munit:
    x,y,z=v[1:4]
    (a,u),(b,beta),(f,h),(g,k)=deals
    a=a*(x+z)+u*z;b=b*(x+z)+beta*x
    f=f*(y+z)+h*z;g=g*(y+z)+k*y
    if marked_sides==1:a=b=0
    corr=(a+f,b+f,a+g,b+g)
    signs=tuple(s+munit*corr[i//(len(signs)//4)]for i,s in enumerate(signs))
   if fixed:signs=tuple(signs[i]for i in lanes)
   if bonus_count:
    z=v[3]
    sides=v[1:1+bonus_count]
    pos=unit*sum(min(z,2*x)for x in sides)
    neg=0 if deep_bonus else unit*sum(min(2*z,x)for x in sides)
    lower=min(min(signs[::2])+pos,min(signs[1::2])+neg)
   else:
    lower=min(signs)
   return c+((v[dim]+lower)<<2)
  if symmetric:
   mid=tuple((x+y)//2 for x,y in zip(vertices[1],vertices[2]))
   vertices=(vertices[0],mid,vertices[2],vertices[3])
  pairs=tuple(combinations(range(dim),2))
  splits=split_pairs.setdefault((dim,symmetric,region,interval,fixed),{})
  pending=[(vertices,int(symmetric),1)]
  if interval!=(0,8)and dim==2:
   ends=tuple(tuple((q*x+(8-q)*y)//8 for x,y in zip(*vertices))for q in interval)
   pending=[(ends,3,1)]
  if region>=0:
   def midpoint(i,j):
    return tuple((a+b)//2 for a,b in zip(vertices[i],vertices[j],strict=True))
   if region>=5:
    if region==5:
     base=(midpoint(0,1),midpoint(0,2),midpoint(1,2))
    else:
     corner=region-6
     base=(vertices[corner],*(midpoint(corner,j)for j in range(3)if j!=corner))
    apex=vertices[3]
    low,high=interval
    if interval==(0,8):
     pending=[((*base,apex),1,1)]
    else:
     def slice_at(k):
      return tuple(tuple(((8-k)*x+k*y)//8 for x,y in zip(a,apex))for a in base)
     lo,hi=slice_at(low),slice_at(high)
     pending=[((lo[0],lo[1],lo[2],hi[2]),4,4)]
     if high!=8:
      pending.extend((((lo[0],lo[1],hi[1],hi[2]),4,5),((lo[0],hi[0],hi[1],hi[2]),4,6)))
   elif region:
    corner=region-1
    part=(vertices[corner],*(midpoint(corner,j)for j in range(4)if j!=corner))
    pending=[(part,1,1)]
   else:
    left,right=midpoint(0,1),midpoint(2,3)
    cycle=[midpoint(i,j)for i,j in((0,2),(0,3),(1,3),(1,2))]
    pending=[((left,right,cycle[i],cycle[(i+1)%4]),1,4+i)for i in range(4)]
  if fixed:
   q=1+2*fixed[0]
   if dim==2:
    rows=((q,4-q),)
   elif len(fixed)==1:
    rows=((4-q,q,0,0),(4-q,0,0,q),(0,q,4-q,0),(0,0,4-q,q))
   else:
    r=1+2*fixed[1]
    rows=tuple((4-q-r+z,q-z,r-z,z)for z in (max(0,q+r-4),min(q,r)))
   vs=tuple(tuple(sum(w*a for w,a in zip(row,col))//4
                  for col in zip(*vertices))for row in rows)
   if len(fixed)==1 and dim==4:
    a,b,d,e=vs
    pending=[((a,b,e),2,2),((a,d,e),2,3)]
    pairs=((0,1),(0,2),(1,2))
   else:
    pending=[(vs,2,1)]
    pairs=tuple(combinations(range(len(vs)),2))
  def with_lower(v):return (*v,linear_lower(v))
  if not munit:pending=[(tuple(map(with_lower,vs)),d,k)for vs,d,k in pending]
  checked=0
  while pending:
   vertices,depth,shape=pending.pop()
   checked+=1
   if munit:
    deals=marks.get(shape)
    if deals is None:
     x,y,z=map(sum,zip(*(v[1:4]for v in vertices)))
     deals=tuple(mark_choice(q,w)for q,w in ((x+z,z),(x+z,x),(y+z,z),(y+z,y)))
     marks[shape]=deals
   if len(vertices)==4:
    a,b,d,e=vertices
    w,x,y,z=map(linear_lower,vertices)if munit else (a[-1],b[-1],d[-1],e[-1])
    for i,k,paired in terms:
     support=entropy_support((a[i],b[i],d[i],e[i]),paired)
     h0,h1,h2,h3=support[k<0]
     w+=k*h0;x+=k*h1;y+=k*h2;z+=k*h3
    for i,pos,neg in biases:
     h0,h1,h2,h3=clip_value((a[i],b[i],d[i],e[i]),pos,neg)
     w+=h0;x+=h1;y+=h2;z+=h3
    if min(w,x,y,z)>0:continue
   else:
    v=list(map(linear_lower,vertices))if munit else [w[-1]for w in vertices]
    for i,k,paired in terms:
     support=entropy_support(tuple(p[i]for p in vertices),paired)
     v=[w+k*h for h,w in zip(support[k<0],v)]
    for i,pos,neg in biases:
     v=[w+(h)for h,w in zip(clip_value(tuple(p[i]for p in vertices),pos,neg),v)]
    if min(v)>0:continue
   if depth>=70 or not pairs:
    raise ValueError('unproved entropy simplex')
   split=splits.get(shape)
   if split is None:
    def distance(pair):
     p,q=(vertices[i]for i in pair)
     x,y=p[0]-q[0],p[1]-q[1]
     if dim==2:return x*x+y*y
     z,w=p[2]-q[2],p[3]-q[3]
     return x*x+y*y+z*z+2*w*w+(y+w)**2+(z+w)**2
    split=max(pairs,key=distance)
    splits[shape]=split
   a,b=split
   p,q=vertices[a],vertices[b]
   if not munit:p,q=p[:-1],q[:-1]
   midpoint=tuple((x+y)//2 for x,y in zip(p,q))
   if not munit:midpoint=with_lower(midpoint)
   for branch,j in enumerate((a,b)):
    part=list(vertices)
    part[j]=midpoint
    pending.append((part,depth+1,2*shape+branch))
  return checked
 def add(terms,key,weight):
  terms[tuple(key)]+=weight
 def subset(mask,dim=4):
  return tuple(4*int(i in mask)for i in range(dim))
 def hb(terms,key,weight):
  add(terms,key,weight)
  add(terms,(4-i for i in key),weight)
 def case(eta,scores,b,o,j,t,branch,marked,region,interval,pair,fixed):
  beta,_,theta_new,_,kappa,*tail=scores
  xis,gammas=tail[:3],tail[3:6]
  deep_counts=(tail[8],tail[46])
  deep_biases=(tail[9],tail[47])
  marked_weight,marked_bias=scores[53:55]
  reserve=tail[-1]
  out_horizon=o-2 if o>=4 else int(o==2)
  if o>=4:
   o=2
  in_horizon=min(out_horizon+1,3)if b==1 and j==1 else 0
  a=(F(0),F(1),F(1),F(1,2))[o]
  kind=(0,1+int(j==1),3)[b]
  cz=1-eta[o]
  constant=-beta*j-reserve-scores[11]*(t-int(b==0 and o==3))
  constant-=sum(xi*(b*int(kind==i)-int(o==i))for i,xi in enumerate(xis,1))
  if in_horizon>1:
   constant-=deep_counts[in_horizon-2]
  if out_horizon>1:
   constant+=deep_counts[out_horizon-2]
  eligible=b==2 and t
  constant-=marked_weight*(eligible-marked)
  if fixed:constant-=len(fixed)*rank_rates[kind-1]
  if j==2 and o:constant+=rank_rates[o-1]
  terms=defaultdict(F)
  biases=[]
  if j==0:
   cu,cv=eta[kind]if b else F(1),eta[kind]if b==2 else F(1)
   constant+=cz+cu+cv
   line=[F(0)]*4
   az,au,av=subset([3]),subset([1,3]),subset([2,3])
   L0,L=(F(1),az)if b==0 else(F(0),au)if b==1 else(F(0),(0,0,0,0))
   wu,wv,wz=1+int(branch==0)-cu,1+int(branch==0)-cv,-cz-int(branch==1)
   for i in range(4):
    add(terms,subset([i]),-1-int(branch==0)+int(branch==1))
   if pair:
    slope=(F(4,5),F(-19,64),F(-1))[out_horizon-1]
    constant-=kappa
    line[3]-=kappa*slope
   elif region>=0:
    fraction=kappa/2 if marked else kappa
    for i in range(4):
     add(terms,subset([i]),-fraction)
    if region>5:
     corner=region-6
     for mask,weight in(([corner],1),([i for i in range(3)if i!=corner],1),([0,1,2],-1)):
      add(terms,subset(mask),fraction*weight)
    elif region==5:
     constant+=fraction
     line[3]-=fraction
    elif region:
     hb(terms,subset([region-1]),fraction)
    else:
     constant+=fraction
   if b:
    for mu,corner in((au,1),(av,2))[:b]:
     for mask,weight in((az,1),(subset([corner]),1),(mu,-1)):
      add(terms,mask,-kappa*a*weight/b)
     biases.append((mu,-gammas[kind-1]))
   for mask,weight in((au,wu),(av,wv),(az,wz)):
    hb(terms,mask,weight)
   if eligible:
    biases.extend((field,-marked_bias)for field in (au,av)[:t])
    line[1]-=kappa/2*(t==2)
    line[2]-=kappa/2
    line[3]-=kappa*t/2
  elif j==1:
   w=eta[kind]if b else F(1)
   az=(2,0)
   L0,L=(F(0),(4,0))if b else(F(1),az)
   line=[F(0),F(0)]
   constant+=cz+w
   h_weight=-w
   if b:
    line[0]+=-kappa*a+theta_new/2+(residual_weight/8 if out_horizon>=2 else 0)
    constant-=theta_new/2
    biases.append(((4,0),-gammas[kind-1]))
    if o==2:
     chain=theta_new*(1-F(1,2**out_horizon))
     constant+=chain/2
     line[0]-=chain/2
   elif pair:
    constant-=kappa
    line[0]-=kappa*(F(2,5),F(-19,128),F(-1,2))[out_horizon-1]
   elif o in(1,2):
    weight=kappa
    constant-=weight
    h_weight-=kappa
    add(terms,az,weight)
    add(terms,(0,4),weight)
    add(terms,(2,4),-weight)
   elif o:
    fraction=kappa/2
    constant-=fraction
    h_weight-=fraction
    for field,weight in((az,1),((0,4),1),((2,4),-1)):
     add(terms,field,weight*fraction)
   else:
    h_weight-=kappa
   hb(terms,az,-cz)
   hb(terms,(4,0),h_weight)
  else:
   az=(1,)
   L0,L=F(5,4),(0,)
   line=[F(0)]
   cost=F(5,8)if marked else F(5,4)if o in(1,2)else 1+a/2
   constant+=cz-kappa*cost
   if o==2:
    constant+=kappa*(F(1,20),F(83,256),F(1,2))[out_horizon-1]
   hb(terms,az,-cz)
  if marked:
   biases.append((az,marked_bias))
  if o:
   biases.append((az,gammas[o-1]))
  if in_horizon>1:
   biases.append(((4,0),-deep_biases[in_horizon-2]))
  if out_horizon>1:
   biases.append((az,deep_biases[out_horizon-2]))
  for k,cap in enumerate((6,9,12,15,18,21,24,27,30)):
   credits=tail[10+4*k:14+4*k]
   for field in((au,av)[:b]if b and j==0 else((4,0),)*bool(b)):
    biases.append((field,-credits[kind-1],cap))
    if in_horizon>1:
     biases.append((field,-credits[3],cap))
   if o:
    biases.append((az,credits[o-1],cap))
   if out_horizon>1:
    biases.append((az,credits[3],cap))
  coeff=-theta_new+(theta_new/2 if o==2 else 0)
  if out_horizon>1:
   coeff+=theta_new*(F(1,2)-F(1,2**out_horizon))
  constant+=theta_new*j+coeff*L0
  for i,weight in enumerate(L):
   line[i]+=coeff*F(weight,4)
  signed=[]
  lanes=[]
  flag_fields=(au,av)[:b]if j==0 else ((4,0),)*b
  incoming=product((0,1),repeat=b)
  for flags in incoming:
   for sign in range(2 if o else 1):
    extra=[]
    for x in range(len(line)):
     z=F(az[x],4)
     v=F(0)
     for i,field in enumerate(flag_fields):
      e=flags[i];q=F(field[x],4)
      c,u,w=signed_rates[3*kind-3:3*kind]
      v-=c*e+u*q+w*(1-q if e else q)
     if o:
      c,u,w=signed_rates[3*o-3:3*o]
      v+=c*sign+u*(1-z if sign else z)+w*z
     if j==0 and b:
      chain=theta_new*(1-F(1,2**out_horizon))
      if b==1:
       e,=flags
       v+=theta_new*(int(x in (0,1))if e else 2*int(x==1))
       v+=chain*(1+(x==3)+(x==(2 if e else 3))-int(x in (1,3)))
      else:
       e,f=flags
       v+=theta_new*(1+max(e,f),2-e,2-f,0)[x]
       v+=chain*(1+min(e,f),1+f,1+e,4-e-f)[x]
     if sign and j==0 and b and out_horizon>=3:
      v+=kappa*(1-z)/4
     if sign and not b and out_horizon>=2 and (pair or j==2):
      v+=kappa*(F(1,2)-F(83,64)*z if out_horizon==2 else F(4,5)-F(41,20)*z)
     if out_horizon>=3:
      q=1-z if sign else z
      if j==0 and b:
       v-=residual_weight*q*(3 if sign else 2)/16
      elif not b and (pair or j==2):
       v-=residual_weight*(F(27,16)-2*z if sign else 7*z/4)/8
     extra.append(v)
    if flags[:len(fixed)]==fixed:lanes.append(len(signed))
    signed.append(extra)
  count=b if j==0 and out_horizon>=2 else 0
  return prove_simplex(line,constant,terms,biases=biases,symmetric=j==0 and b!=1 and t!=1 and region<0 and not fixed,
                        region=region,interval=interval,signed=signed,
                        credit=(count,out_horizon>=3,kappa/(2*count)if count else F(0)),marked_rate=kappa/512 if eligible else F(0),marked_sides=eligible,fixed=fixed,lanes=lanes)
 def prove_combined(eta,scores):
  need(len(eta)==4 and eta[0]==0 and all(-2<=x<=2 for x in eta),'signed carrier coefficients')
  need(len(scores)==57 and min(scores[:5])>=0 and scores[-1]>0 and scores[1]==scores[3]==0 and scores[-2]==scores[4],'dual coefficients')
  need(0<=residual_weight<=scores[4],'residual pair weight')
  need(scores[11]>=0 and scores[12]==0,'dominated flag counts')
  total=0
  for b,o,j,t in product(range(3),range(6),range(3),range(3)):
   if not(j==0 or(j==1 and b<2)or(j==2 and b==0))or(b!=2 and t):
    continue
   if((b,o)in PART_ZERO)!=(part==0):
    continue
   for branch in range(2 if j==0 else 1):
    for marked in (int(b==0 and o==3),):
     regions=(-1,)
     if b==0 and j==0:
      if o==0:
       regions=range(5)
      elif o in(1,2,3,4,5):
       regions=range(5,9)
     for region in regions:
      intervals=(((0,8),0),)
      if b==0 and o in(2,4,5)and j<2:
       lo,hi=((4,5)if o==2 else(3,6))if j else(1,3)
       intervals=(((0,lo),0),((lo,hi),1),((hi,8),0))
      for interval,pair in intervals:
       count=(t if b==2 else b)
       for r in range(count+1):
        for fixed in product((0,1),repeat=r):
         total+=case(eta,scores,b,o,j,t,branch,marked,region,interval,pair,fixed)
  return total
 mrates=((82,-595),(56,-313),(29,-136),(0,0),(-33,114),(-68,217),(-107,312),(-150,405))
 need(all(F(16-i,12)**256>=F(2)**c and F(3*i,16-i)**256>=F(2)**s
          for i,(c,s)in enumerate(mrates,1)),'marked resolving incidence weights')
 edges=[(a-b,v-u)for (a,u),(b,v)in zip(mrates,mrates[1:])]
 need(all(a*d<b*c for (a,c),(b,d)in zip(edges,edges[1:])),'increasing marked breakpoints')
 def mark_choice(q,w):
  lo,hi=0,len(edges)
  while lo<hi:
   m=(lo+hi)//2
   a,b=edges[m]
   if w*b>q*a:lo=m+1
   else:hi=m
  return mrates[lo]
 d=100000000
 eta=[F(x,d)for x in (0, -192501025, 106762744, -96209652)]
 scores=[F(x,d)for x in (5674745, 0, 2343608, 0, 7847666, -282511888, 95250915, -187339211, -69337321, 157980254, 67421293, 0, 0, -110285352, -346399864, 283300000, -229975214, 347507544, 355740000, 275324039, -300285736, 240230497, 374356462, 351384506, -50786899, -115992605, 287163499, 103253095, 108709603, 284142200, -150471107, 98391214, -24938838, 54675232, 4142887, 51354570, -40569560, 30863834, 27666954, 111909436, 11168081, 53581774, -35087025, -38688343, -75806580, 51512009, -13773537, 248138572, -107545182, -19032477, 377451358, -111646503, -343018131, -4231366, 4153082, 7847666, 1)]
 signed_rates=[F(x,d)for x in (-2191751, 2458255, -38542487, 0, -248095, -37235196, -4348928, 8697855, -26556891)]
 residual_weight=F(7846097,d)
 rank_rates=[F(x,d)for x in (1091955, 825009, 1501713)]
 checks=prove_combined(eta,scores)
 saving=sum(scores[i]for i in(0,3,4,56))
 cap=F(cap_numerator,10**20)
 coefficient=1+1/(2*cap)+saving/3
 return{'checks':checks,'local_reserve':str(scores[-1]),'saving':str(saving),
 'semantic_increment':str(saving/3),'coefficient':str(coefficient)}
if __name__=='__main__':print(json.dumps(compute(),sort_keys=True))
