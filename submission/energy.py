"""Directed energy."""
from fractions import Fraction as F
from itertools import combinations
from functools import cache
from math import comb,factorial,isqrt,prod
THRESHOLD=F(-293,800)
def ck(q,m):
 if not q:raise ValueError(m)
@cache
def half_moments():
 m=[1,1]
 for j in range(2,16):m.append((j-1)*m[j-2])
 w,p=[1],[0]
 for k in range(1,65536):
  w.append(w[k>>4]*m[k&15])
  p.append(p[k>>4]+(k&1))
 return w,p
def odd_orthant(K,det,signs,scales,threshold,pi,nodes,hermite):
 """Directed halfnormal integration over a Gaussian orthant."""
 scale=1<<44
 size=len(signs)
 degree=2*len(nodes)-1
 if (size,degree) not in ((6,7),(8,5)) or len(scales)!=size or min(scales)<=0:
  raise ArithmeticError('orthant chart')
 v=[F(s*sign,10000)for s,sign in zip(scales,signs)]
 sums=[sum(row)for row in K]
 def floor_scaled(x):
  return x.numerator*scale//x.denominator
 def sqrt_floor(x):
  return F(isqrt(x.numerator*scale**2//x.denominator),scale)
 units=[1<<(4*i)for i in range(size)]
 p0=floor_scaled(-threshold**2*sum(sums)/2)
 # Euler weights use the least index in each nonconstant monomial.
 first_terms=[[]for _ in range(size+1)]
 for i,u in enumerate(units):
  first_terms[i].append((u,floor_scaled(-threshold*v[i]*sums[i])))
  first_terms[i].append((2*u,2*floor_scaled((1-v[i]**2*K[i][i])/2)))
  for j in range(i):
   first_terms[j].append((u+units[j],floor_scaled(-v[i]*v[j]*K[i][j])))
 weight,odd=half_moments()
 mu_lo=sqrt_floor(2/pi[1])
 mu_hi=sqrt_floor(2/pi[0])+F(1,scale)
 powers_lo=[mu_lo**j for j in range(size+1)]
 powers_hi=[mu_hi**j for j in range(size+1)]
 def first(code):return ((code&-code).bit_length()-1)//4
 for i in range(size):
  first_terms[i+1]=first_terms[i]+first_terms[i+1]
 powers={0:1}
 moments=[(F(1),F(1))]
 for k in range(1,degree+1):
  expanded={}
  for monomial,coefficient in powers.items():
   for code,value in first_terms[first(monomial) if monomial else size]:
    key=monomial+code
    expanded[key]=expanded.get(key,0)+coefficient*value
  sums_by_odd_degree=[p0**k]+[0]*size
  powers={0:p0**k}
  for code,numerator in expanded.items():
   i=first(code)
   coefficient=k*numerator//((code>>(4*i))&15)
   powers[code]=coefficient
   lo,hi=code&65535,code>>16
   sums_by_odd_degree[odd[lo]+odd[hi]]+=coefficient*weight[lo]*weight[hi]
  lower=sum(c*(powers_lo[j]if c>=0 else powers_hi[j])
  for j,c in enumerate(sums_by_odd_degree))
  upper=sum(c*(powers_hi[j]if c>=0 else powers_lo[j])for j,c in enumerate(sums_by_odd_degree))
  moments.append((lower/scale**k,upper/scale**k))
 expectation=hermite(nodes,moments)
 if expectation<=0:raise ArithmeticError('positive orthant minorant')
 return(prod(abs(x)for x in v)*expectation/
 (2**size*(sqrt_floor(det)+F(1,scale))))
def energy(rr,uu,M,G,pair):
 def ck(x,m):
  if not x:raise ValueError(m)
 rows={''.join(map(str,k)):(F(l,M['S']),F(h,M['S']))for k,(l,h)in rr.items()}
 U5=F(uu[1],M['S']);t,a=F(41,57),F(28,57);u=-THRESHOLD
 @cache
 def r(k):
  k=abs(k)
  return F(1)if not k else a*t**(k-1)+t**k/2+(F(k-1,2)*a*a*t**(k-2)if k>=2 else 0)
 @cache
 def gram(n):return[[r(i-j)for j in range(n)]for i in range(n)]
 invdet=G['inverse']
 def cone(n):
  K,_=invdet(gram(n));s=[(-1)**i for i in range(n)];D=sum(map(sum,K))
  A=[[s[i]*s[j]*K[i][j]for j in range(n)]for i in range(n)];b=[s[i]*sum(K[i])for i in range(n)]
  best=F(0);count=0
  for k in range(1,n+1):
   for face in combinations(range(n),k):
    C=[[A[i][j]for j in face]for i in face]
    J,_=invdet(C)
    v=[b[i]for i in face]
    z=[sum(J[i][j]*v[j]for j in range(k))for i in range(k)]
    count+=2
    if min(z)<0<max(z):continue
    q=sum(x*y for x,y in zip(z,v))
    ck(q==sum(z[i]*C[i][j]*z[j]for i in range(k)for j in range(k)),'face')
    best=max(best,q)
  ck(0<best<=D,'cone bound');return D,best,count
 @cache
 def alt(n):
  lo=hi=F(0)
  for k in range(0,n+1,2):
   for S in combinations(range(n),k):
    sign=(-1)**sum(S);g=''.join(str(j-i)for i,j in zip(S,S[1:]));g=min(g,g[::-1]);l,h=rows[g]if g else(F(1),F(1));lo+=l if sign>0 else-h;hi+=h if sign>0 else-l
  return lo/2**(n-1),hi/2**(n-1)
 def ep(z,k):return sum((-z)**j/factorial(j)for j in range(k+1))
 def arct(x):
  p=sum((-1)**k*x**(2*k+1)/(2*k+1)for k in range(24));return p,p+x**49/49
 p,q=arct(F(1,5));x,y=arct(F(1,239));pi=(16*p-4*y,16*q-4*x)
 Q=10**24
 def down(x):return F(x.numerator*Q//x.denominator,Q)
 def up(x):return-down(-x)
 def iv(x):return(down(x),up(x))
 def add(x,y):return down(x[0]+y[0]),up(x[1]+y[1])
 def mul(x,y):
  a=[v*w for v in x for w in y];return down(min(a)),up(max(a))
 def scale(x,s):return mul(x,iv(s))
 def div(x,y):ck(y[0]>0,'positive denominator');return mul(x,(down(1/y[1]),up(1/y[0])))
 def sqrt(x):
  ck(x>=0,'sqrt domain');k=isqrt(x.numerator*Q**2//x.denominator);l,h=F(k,Q),F(k+1,Q);ck(l*l<=x<h*h,'sqrt bracket');return l,h
 def pinned_orthant(C,sg):
  n=len(C)
  block=[[M['ONE']for _ in sg]for _ in sg]
  prob=M['ONE']
  for i,j in combinations(range(n),2):
   rho=M['div'](M['rat'](C[i][j]),M['root'](M['rat'](C[i][i]*C[j][j])))
   block[i][j]=block[j][i]=rho
   pair=M['moment']([[M['ONE'],rho],[rho,M['ONE']]])
   prob=M['add'](prob,M['scale'](pair,sg[i]*sg[j]))
  if n==4:
   prob=M['add'](prob,M['scale'](M['moment'](block),prod(sg)))
  prob=tuple(F(x,M['S'])for x in M['intdiv'](prob,2**n))
  ck(0<prob[0]<=prob[1]<=1,'orthant')
  return prob
 def line_moments(n,D,beta,degree):
  C=gram(n)
  s=[(-1)**i for i in range(n)]
  @cache
  def block(pins):
   keep=[i for i in range(n)if i not in pins]
   Kp,det=invdet([[C[i][j]for j in pins]for i in pins])
   weights=[sum(row)for row in Kp]
   v=[1-sum(C[i][j]*w for j,w in zip(pins,weights))for i in keep]
   d=D-sum(weights)
   return keep,Kp,det,v,d
  @cache
  def base(pins):
   if not pins:
    prob=scale(alt(n),F(1,2))
    return prob,prob
   ck(len(pins)%2==0,'even face base')
   keep,Kp,det,_,_=block(pins)
   cond=[[C[i][j]-sum(
   C[i][a]*Kp[g][h]*C[b][j]
   for g,a in enumerate(pins)for h,b in enumerate(pins)
   )for j in keep]for i in keep]
   prob=pinned_orthant(cond,[s[i]for i in keep])
   norm=sqrt(det)
   for _ in range(len(pins)//2):
    norm=mul(norm,scale(pi,2))
   return div(prob,norm),prob
  @cache
  def moment(pins,k):
   if k==0:
    return base(pins)[0]
   keep,_,_,v,d=block(pins)
   result=scale(moment(pins,k-2),(k-1)*d)if k>=2 else iv(F(0))
   for j,w in zip(keep,v):
    result=add(result,scale(moment(tuple(sorted((*pins,j))),k-1),s[j]*w))
   return result
  out=[scale(moment((),k),2)for k in range(2,degree+1,2)]
  deriv=iv(F(0))
  for i,j in combinations(range(n),2):
   c=C[i][j]
   prob=base((i,j))[1]
   deriv=add(deriv,scale(mul(div(sqrt((1-c)/(1+c)),pi),prob),2*s[i]*s[j]))
  old=add(deriv,scale(alt(n),D))
  ck(max(out[0][0],old[0])<=min(out[0][1],old[1]),'second moment regression')
  ck(0<out[0][0]<=out[0][1]<alt(n)[1]*n*beta,'face moment')
  ck(all(0<lo<=hi for lo,hi in out),'positive moments')
  return out,n*(n-1)
 tails=[]
 witnesses=[
 ((4750,3131,3354,3376,3408,2928,4745,6634),625,(4398,3215,3330,3406,3380,3016,4635,6201),424),
 ((4762,3107,3374,3352,3312,3361,3411,2921,4741,6630),832,(4400,3192,3340,3388,3290,3396,3375,3008,4627,6188),668),
 ((4766,3100,3378,3335,3334,3342,3322,3349,3411,2923,4743,6634),908,(4397,3183,3341,3371,3303,3382,3290,3386,3371,3007,4624,6183),808),
 ((4765,3095,3376,3329,3337,3324,3342,3330,3323,3344,3406,2926,4740,6633),1025,(4391,3174,3336,3363,3302,3364,3302,3371,3287,3378,3364,3005,4614,6170),1028),
 ]
 for _ in range(3):
  p,z,q,w=witnesses[-1]
  j=len(p)//2
  witnesses.append((p[:j]+p[j-2:j]+p[j:],z+175,q[:j]+q[j-2:j]+q[j:],w+215))
 for i,(sp,zp,sm,zm)in enumerate(witnesses):
  n=len(sp);K,det=invdet(gram(n));sg=[(-1)**j for j in range(n-1)]+[1]
  if i==0:
   H=G['hermite_expectation']
   x,y=(tuple(F(v,1024)for v in ns)for ns in ((-1439, -13, 1448), (-1572, -276, 1074)))
   w=odd_orthant(K,det,sg,sp,THRESHOLD,pi,x,H)+odd_orthant(K,det,[-j for j in sg],sm,THRESHOLD,pi,y,H)
   ck(w>F(37311636,10**11),'eight-site lower bound')
   tails.append(down(w))
  else:
   q=G['lower_signed'](K,det,sg,sp,zp)+G['lower_signed'](K,det,[-j for j in sg],sm,zm)
   tails.append(F(q,G['S']))
 k2=(1-r(2))/(1+r(2));klo=sqrt(k2)[0]
 sint=sum((-1)**j*sum(F(comb(j,l),2*l+1)*k2**l for l in range(j+1))*u**(2*j+2)/(2**j*factorial(j)*(2*j+2))for j in range(6));ck(sint>0,'bivariate lower');sdrop=down(klo/pi[1]*sint)
 D4,b4,n4=cone(4);D5,b5,n5=cone(5);D6,b6,n6=cone(6)
 degree=6
 M4,f4=line_moments(4,D4,b4,degree)
 M5,f5=line_moments(5,D5,b5,degree)
 M6,f6=line_moments(6,D6,b6,degree)
 def upper_delta(n,precision_sum,beta,moments):
  total=sum(h*u**(2*k)/factorial(2*k)for k,(_,h)in enumerate(moments,1))
  k=len(moments)
  term=moments[-1][1]*u**(2*k)/factorial(2*k)
  ratio=u*u*beta*(n+2*k)/((2*k+1)*(2*k+2))
  ck(0<ratio<1,'remainder')
  total+=term*ratio/(1-ratio)
  exp_bound=ep(precision_sum*u*u/2,6)
  ck(0<exp_bound<1,'exp upper')
  return up((exp_bound-1)*alt(n)[0]+exp_bound*total)
 A3_delta=upper_delta(4,D4,b4,M4)
 e5=ep(D5*u*u/2,5)
 ck(0<e5<1,'exp lower')
 P5=alt(5)[1]
 cosh5=sum(l*u**(2*k)/factorial(2*k)for k,(l,_)in enumerate(M5,1))
 A4_drop=up(P5*(1-e5)-e5*cosh5)
 A5_delta=upper_delta(6,D6,b6,M6)
 c1=U5-sdrop+A3_delta+A4_drop+A5_delta-sum(tails)
 G0=G['compute'](odd_orthant);gain=F(G0['gain_lower'],G0['scale'])
 motif=G['edge_motif_gain']();bonus=F(motif['gain_lower'],motif['scale']);cost=c1-gain-bonus-pair
 ck(U5<=F(73670242301,2**39)and gain>0 and bonus>0 and pair>=F(47991230402,10**15)and bonus+pair<F(21,100000)and 0<cost<F(12233148667589113194,10**20),'active energy bound')
 ck(f4+f5+f6==62 and n4+n5+n6==218,'face accounting')
 return{'energy_upper_1e20':-(-cost.numerator*10**20//cost.denominator),'eight_site_discount_lower_1e12':int(tails[0]*10**12),'first_pass_upper_1e12':-(-c1.numerator*10**12//c1.denominator),'gain_lower_1e12':gain.numerator*10**12//gain.denominator,'motif_gain_lower_1e12':bonus.numerator*10**12//bonus.denominator,'pair_gain_conditional_1e15':pair.numerator*10**15//pair.denominator,'motif_oriented_lower_1e12':motif['oriented_lower']*10**12//motif['scale'],'motif_cases':motif['cases'],'U5_upper_1e12':-(-U5.numerator*10**12//U5.denominator)}
S=1<<50;H=1024
def up(p,q):return-((-p)//q)
def ppi():
 def tan(x):
  a=sum((-1)**i*x**(1+2*i)/(1+2*i)for i in range(20));return a,a+x**41/41
 a,b=tan(F(1,5));c,d=tan(F(1,239));return 16*a-4*d,16*b-4*c
pi0,pi1=ppi()
def sqrt(z):
 v=isqrt(z.numerator*S*S//z.denominator);return v,v+1
phi0=sqrt(1/(2*pi1))[0];phi1=sqrt(1/(2*pi0))[1]
def ex(x):
 u,v=x
 lo=hi=S;n=m=S;prev=(S,S)
 for j in range(1,10):
  n=n*u//(S*j);m=up(m*v,S*j);prev=(lo,hi)
  if j&1:lo-=m;hi-=n
  else:lo+=n;hi+=m
 hi=prev[1]
 for k in range(12):
  lo=lo*lo//S;hi=up(hi*hi,S)
 return lo,hi
P=[(S//2,S//2)]
for k in range(5*H):
 r=(2*k+1)**2*S;den=8*H**2*4096
 l,u=ex((r//den,up(r,den)));a,b=P[-1]
 P.append((a+phi0*l//(S*H),b+up(phi1*u,S*H)))
ERR=up(5*S,24*H**2)
tail=up(S,2500000)
r=(5*H*2)**2*S;den=8*H**2*4096
_,erp=ex((r//den,up(r,den)))
if phi1*erp>5*tail*S:raise ValueError('tail bound')
C=[(max(0,l-ERR),min(S,h+ERR))for l,h in P]
def cdf(i):
 if i<0:
  l,u=cdf(-i);return S-u,S-l
 if i>=5*H:return S-tail,S
 return C[i]
def Fcdf(x):
 n=x.numerator*H;d=x.denominator
 return cdf(n//d)[0],cdf(up(n,d))[1]
def pm(F0,F1):
 return max(0,Fcdf(F1)[1]-Fcdf(F0)[0])
def bound(T,r0):
 a=F(28,57);r=F(97,114);cv=1-r*r;rho=(r0-r*r)/cv
 z=sqrt(cv*(1+rho)/2)[1];bb=sqrt((1-rho)/(1+rho))[0]
 STEP=64
 weights=[2*pm(F(i,STEP),F(i+1,STEP))for i in range(5*STEP)]
 negative_upper=[C[0][1]]+[S-v[0]for v in C[1:-1]]
 suffix=[0]*(len(weights)+1)
 for j in range(len(weights)-1,-1,-1):suffix[j]=suffix[j+1]+weights[j]
 ans=0;i=0
 while T+F(i,STEP)<5:
  l=T+F(i,STEP);h=min(F(5),l+F(1,STEP));q=pm(l,h);i+=1
  b=(T-r*l)*H*S/z
  top=b.numerator*(S*STEP);den=b.denominator*(S*STEP);dist=H*bb*b.denominator
  ck(top<=0 and dist>0,'negative Gaussian rectangle arguments')
  stop=max(0,min(len(weights),up(top+5*H*den,dist)))
  total=suffix[stop]*tail
  for j in range(stop):total+=weights[j]*negative_upper[-up(top-dist*j,den)]
  ans+=q*up(total,S)
 return F(up(ans,S)+3*tail,S)
def initial(rows,M,G):
 T=THRESHOLD;a=F(28,57);r=F(97,114);s=F(529,722);cross=2*a*r-a*a
 y=bound(T,s);x=bound(T,cross)
 if not(0<y<F(29,1000)and 0<x<F(13,1000)):raise ValueError('pin triple enclosure')
 M0=M['S'];Q=M['rat'];ad=M['add'];dv=M['div'];rt=M['root'];mu=M['mul'];sc=M['scale'];asin=M['asin'];K=M['K']
 signs=(-1,1,1,-1);R=[[M['cov'](i-j)for j in range(4)]for i in range(4)]
 a0=Q(1);sl=Q(0)
 for i in range(4):
  for j in range(i+1,4):a0=ad(a0,sc(rows[j-i,],signs[i]*signs[j]))
  h=Q(1);ids=[j for j in range(4)if j!=i]
  for j in ids:
   for k in ids:
    if j>=k:continue
    c=Q(R[j][k]-R[j][i]*R[k][i]);d=rt(Q((1-R[j][i]**2)*(1-R[k][i]**2)))
    h=ad(h,sc(mu(K,asin(dv(c,d))),signs[j]*signs[k]))
  sl=ad(sl,sc(h,-signs[i]))
 a0=ad(a0,rows[(1,1,1)]);penalty=F(a0[1],16*M0)
 if not(sl[0]>0 and penalty<F(121,10000)):raise ValueError('0110 log-concave penalty')
 H=G['S'];half=G['lower'];extra=F(half(('11','11',''),[2500,4800,6200,4800,6200],500,-1),H)
 tri=F(half(('0','0','0'),[1236,3201,3201,2560],854),H)
 if not(extra>0 and tri>0):raise ValueError('forced five-pin and overlap')
 out=F(Fcdf(-T)[0],S)-y-2*penalty-4*x+4*tri+extra
 before=out
 if before<=F(504,1000):raise ValueError('combined calibration margin')
 out-=y
 if out<=F(1,2):raise ValueError('initial strict majority after motif repair')
 def ceil(x):return up(x.numerator*10**12,x.denominator)
 return{'initial_lower_1e12':-ceil(-out),'initial_before_motif_lower_1e12':-ceil(-before),'Y_upper_1e12':ceil(y),'cross_pair_upper_1e12':ceil(x),
 '0110_upper_1e12':ceil(penalty),'overlap_lower_1e12':-ceil(-tri),'11011_putback_lower_1e12':-ceil(-extra)}
def b95(s):
 s=s.encode('ascii')
 buf=bytearray()
 for i in range(0,len(s),50):
  part=s[i:i+50];v=0
  for c in part:v=95*v+c-32
  buf.extend(v.to_bytes(len(part)*41//50,'big'))
 return bytes(buf)
from bisect import bisect_left
from fractions import Fraction
from lzma import decompress
from math import lcm
def stencils(axis,step):
 nodes=[axis[0]]
 for j,x in enumerate(axis[1:],1):
  if x-nodes[-1]>=step or j==len(axis)-1:nodes.append(x)
 taps=[]
 for x in axis:
  j=bisect_left(nodes,x)
  if nodes[j]==x:
   taps.append((j,(1,),1));continue
  st=min(max(0,j-2),len(nodes)-4)
  span=nodes[st:st+4]
  w=[]
  for a in span:
   n=Fraction(1)
   for b in span:
    if b!=a:n*=Fraction(x-b,a-b)
   w.append(n)
  d=lcm(*(a.denominator for a in w))
  taps.append((st,tuple(int(a*d)for a in w),d))
 return nodes,taps
def expand(values,W,X,Y):
 rail=[]
 for start in range(0,len(values),W):
  row=values[start:start+W]
  out=[]
  for i,c,d in X:
   if len(c)==1:
    out.append(512*row[i])
   else:
    a,b,u,v=c
    zz=a*row[i]+b*row[i+1]+u*row[i+2]+v*row[i+3]
    out.append((512*zz+d//2)//d)
  rail.append(out)
 out=[]
 for i,c,d in Y:
  den=512*d
  if len(c)==1:
   out.extend((x+256)//512 for x in rail[i])
  else:
   a,b,u,v=c
   out.extend((a*r+b*s+u*t+v*w+den//2)//den for r,s,t,w in zip(*rail[i:i+4]))
 return out
def shape(main,contrast,dx,dy):
 nx,X=stencils(main,dx);ny,Y=stencils(contrast,dy)
 return nx,ny,X,Y
def smooth(knots,main,contrast,Q,dx,dy,offset):
 nx,ny,X,Y=shape(main,contrast,dx,dy)
 size=len(nx)*len(ny)
 assert len(knots)==2*size
 raw_m=expand(knots[:size],len(nx),X,Y)
 raw_h=expand(knots[size:],len(nx),X,Y)
 h=[min(Q//2,max(0,v-offset))for v in raw_h]
 m=[min(b,max(0,v-offset))for v,b in zip(raw_m,h)]
 return[2*r for r in m+h]
def read_words(data,count):
 vals=[];value=shift=0
 for n in data:
  value|=(n&127)<<shift
  if n>=128:
   shift+=7
   ck(shift<=14,'endpoint item size')
  else:
   ck(value<65536,'endpoint item range')
   vals.append(value);value=shift=0
 ck(not shift and len(vals)==count,'endpoint shape')
 return vals
def predict_endpoint(conv,Q,knots_path,fixes_path,dx,dy,offset):
 x,y,_,_=shape(conv.main,conv.y_axis,dx,dy)
 W=len(x);S=W*len(y);knots=read_words(decompress(b95(knots_path)),2*S)
 mix=[]
 for part in(knots[:S],knots[S:]):
  row=[]
  for k,t in enumerate(part):
   a,b=divmod(k,W)
   t=(t>>1)^-(t&1)
   if b:t+=row[k-1]
   if a:
    t+=row[k-W]
    if b:t-=row[k-W-1]
   row.append(t)
  mix.extend(row)
 lower=smooth(mix,conv.main,conv.y_axis,Q,dx,dy,offset)
 N=conv.msize
 data=decompress(b95(fixes_path))
 mask_bytes=(2*N+7)//8
 ck(len(data)>=mask_bytes,'endpoint mask size')
 ck(data[mask_bytes-1]<(1<<(2*N-8*(mask_bytes-1))),'endpoint padding')
 selected=[i for i in range(2*N)if data[i>>3]&(1<<(i&7))]
 fees=read_words(data[mask_bytes:],len(selected))
 v=lower
 for i,fee in zip(selected,fees):
  v[i]-=2*fee
 m,p=v[:N],[Q-x for x in v[N:]]
 ck(all(0<=x<=Q-y and 0<=y<=Q for x,y in zip(m,p)),'endpoint simplex')
 return m,p
def decode_signed(conv,Q,knots,fixes):
 from lzma import FORMAT_RAW,FILTER_LZMA2
 filters=[{'id':FILTER_LZMA2,'dict_size':1<<20}]
 def unpack(data):
  return decompress(b95(data),format=FORMAT_RAW,filters=filters)
 nx,ny,X,Y=shape(conv.main,conv.y_axis,8,8)
 w=len(nx);size=w*len(ny)
 words=read_words(unpack(knots),2*size)
 out=[]
 for mixed in(words[:size],words[size:]):
  values=[]
  for i,k in enumerate(mixed):
   r,c=divmod(i,w)
   delta=(k>>1)^-(k&1)
   if c:delta+=values[i-1]
   if r:delta+=values[i-w]
   if r and c:delta-=values[i-w-1]
   values.append(delta)
  out.extend(expand(values,w,X,Y))
 n=conv.msize;at=(2*n+7)//8
 body=unpack(fixes)
 ck(len(body)>=at and body[at-1]<(1<<(2*n-8*(at-1))),'endpoint mask')
 picked=[i for i in range(2*n)if body[i>>3]&(1<<(i&7))]
 for i,k in zip(picked,read_words(body[at:],len(picked))):
  out[i]+=(k>>1)^-(k&1)
 m,p=[2*x for x in out[:n]],[Q-2*x for x in out[n:]]
 ck(all(0<=x<=Q-y and 0<=y<=Q for x,y in zip(m,p)),'endpoint simplex')
 return m,p
FIRST_CORRECTIONS=(r'''[\"0|zj=AzA{DXzw;{}=T^C!Q\A vk}-"{F88s|ZX5<~-!zM"8`mi~@V#(R#/f&mrCH]@<oPXxla=&wK5QAjDX8O`+2>FT|'H+V\Xv-ty!TH#fr<-?e<(poMBg@'gRsm/PAGem5s;gZLsK?QR2P0apC8w.trAYqktJ:bk8Az"9]%FfI+&P3|b#bl/K!n4{xCD%HP~_o<9Nyl?f9:hAn17X[4>TSUfAr3,e@R\q~XT.\l|IP5)698],9d/3R93|sI-r.&/h*HmjX<pP6oYrw:8cBU-@~T,nPFOol\T?VBSCR#R8VK*ySK2H<z8Yu0iazfk}.`[!Y8#&)XQ7J;8mY'7W6T]jW~o!YmApU%LZ3rB/GYJX^}XfsCGuma>ocyP"=*jP@qi#@m5?>GKc=MFEVt,E#-b]D4i43*V:d<E E2}''0#U5u+])Y`bWQ`k*$}.%ILJ1#<UJhAg[sWX#<.PX]'ujg,FEp;y\{iK=)FY%Az`5(Nw_Ikfc$%agtt,9<t|cFeic`kI<~:M? H!*)U1_Y?Se7V7TP+eNA*a)0%}+.q?Z$<vtFpABCKQC.eW_P8g:P>oD@Ku%k1N>0ZY(*c_=oXNV&*5`gD|1LR|.O%}dF|#|2a?}I~;\jHW&i{-ht}T'zJ{h#"l) ts\%iK.?;JD[oZBa@{66obpR^G$OL;ftI-3cl`M%450B3+B+PfIl%uWg"RphLCL[1;feSNgVb},f<>M=S)N?s|]/UX&)&`gg`w+K2LFiav|Q[!^e7+w-TU63<<%$N|m~yB)s/iN{PDw1~R,[1I%ov7OM_k>D^`sRBIGdLz(-] v,Q/6#qNuYHSU|J/!X/`oSC8:BOx5%Se'KudqjU.Gx~g-Dci}a/!\DJ6#^+_%^`GOW||b+r@DVt|S<0.v;*E"4sXO@UAlR\d<Y[~T(28/gJ$F[Qmm?Mk9OxK%Y(R6Y#)yVyS];M,F]Mip!O7(\S=7xm)k:%0=cROMK>@]KwF}}}{<)i,Yr8bt{915.:xt)ZyO8B]-Uo'bZbc$"`i{6?{dWEg,+k_S-sAyh_59[zlN_Ddy&:x?"1mfPcpp.^kA3e[:ctGIkeIZ/flwW[,p_PkleSSD9<%z(X%0yc=yG,R?HPKz"JjuF>aa*9:spn|ZfUO^ynE6W-M)4?8"b9ZvkM3GvdKI7mL$1<>h[}p"0$]b\\y'N86_e.=`l]OSdq9v{)s]B`/<I_KhibPsc[|V+}R?>9x8G`#MnDOfCvU}=k-D>1m>pVq.^iS,v/z`^x$kRA7k=mV-.;<;E4RtzaE:F!eD~zB5%|2"$c(MLU9}(F8p<g>)!PYsMv{J(/AtFJ8fuqZuDSkEr'0A6bh_jtvUZp4:x=mr?b,4>0C)@&#Tw<]xe *rvKKv~/\zk|)N4"O!=zoS6thw*0xcO'luZcI-[?w$=T8v1q[+kIY;LZnEqT7@}Sz!E$MP{U:e><m[ H7?'yqq~>^?:f(vbIe>JG/.D`TYRpzVw5NGO{q!JhNiurdF98i3f.@FY.hW)]'$zA.3AQ,?j_q\hrq8mjhy\xMVCNkXu==dIRRT(,P+6?KEr{`2y/DChu7B]D'<f_?|oE6[.|7`^nW&ri5"rC,EtT.j[H6o Z6qh@P5:MP!8HQV{Ji'f>5fB\,g%,9LiZFHhNOwSY%G#](^Z~ZK}*9P@6xZ$s6^<5^mpw|6iw"H:42\tnH'4kx~@UBnvX$:s@Y~Z7v:]RnTwt4N/Rse$,)[%`<?)o8d\\ <`x1.4grvA6zS2v&vbm1],$nz|'!-Y` tBbKYsuGRPO}!jTjU\Jx#C5e6D{I-pVu-|a|Y-!Nuz[IF(?\<?*9o+0@0"qPS|4[SJ18nem=P//+M"R"?n~uZ*lRhkJCkOR_SxgB{H+<;Y-& jKHl[e-HzpTQ5G8e`e!F:AO`w=k<IhtE^@G0Ez)kPqQif)0XEm/5->*U@A{`b|smruH5)C'QS>Hc3Ayen-A#R?sZ>l`.0v@-dYt&O9!2CPMc}Tx8S46Lh\B0 ;u^a-T=FQ1H>CPD#7U+b$2H@~8/u^mTh4=z~ELQu(b,]xoqlY(*J:CX4^OB1A\6D&/2nlBREY-kHt*X/Ae21*kqK;i 0<d{91Nc!9U'}Ayi!,CD:TBth;W(lfKiQ<f3Z+\Zt7gzF^r"DH<!*V-~_MbNY4Vfp=ZGeLE[}e|lm#P]JI*giP_%R*6Or(uy1`9@$A,<8xCJ'c5(%r;$sI][s$w';cgYk'ufx; O&U&E7[U`0!w6&El=D!6)"([G*BHEq0Km`$RO,Aza!m<^{pfcJ$KskYDw%<-p{/3L>jUn+Yry.@iVL,4![=88=y,e)joC}L5$UXmaIcREu>tm1-EkHC]BZX2N;Q"5%idq _h>2a=oUT6n7Y(7hfW}Ti>8?S37EV{_\w<lmo/w4k"2RZ>mqFjlKThjAaJ33%-NvRHt%Pbg87F^z~vm#A<Hlu#-/6rpACC:n.V]iB*~#/b=lWNXD.f\pkq{qmPfqx#<x;DKs<`w,VjO"}1)q-e^zJ0-%M+[ Bks\ch,-*)VS?IZ_cY(;X0P0c[Drp)(/b51,>aShW-a"9*{>)F:Z`+]V8^g/x0KuVHsP<[u${>..n6-7ix3t$.fp)8MF618W}9lMZ4Kf!&~DT`7gR^~v(V{r'v{>=I4f#59(q[|50b_Fz\A+-AfXIfQBr2g`d~5\RUVTSyXv~-\Bz\7*a,1imqx|v'y+APX &wPP6G<<G/z/c:P#B!0vX.<"G~}<nCoBu/M{8%BS(1)$lX@;+m0zU8EK)RLwkieY_iGHCH@-]|it_K_o\meGi(V<0}<XowEM$xY7@F(vA):I>}CX154:j)rLFu/xWH}UV=1_Q[Y|.9W_%z00~&HH^c Rg;'F87>pDsQB5#%*5T Y_?yFxK\Q_N'|"p7`p$x[p*YN(Da+"(bas0fr&!rMGQDn^zm_iCS6K*n8,}_eYZS<dvm"oY@3PO5N0,,@,#<~r<{d[[)\bvib:quu49NS Cb]mF-=Ut^<a.*e)0:]$o'\~#l5)Tt09wpH"\4H`EIN]@*7q&GwI3",!Bpz+J(&Emx%&"/v<e|ghk[~heKK;y8V%0J0Pe0+>$V9``i1V()< >a'0R%lzqRY_<=q>+Z$Y7<I;'.'9JY:W'VV~~9tL@U#inX@"@}PARsV)sN"!HbedE`!JZ}3{F!V*$0ySWID=hFq"x$nw&[$+iML|1"&S1G*f6iM?wPAbl2ht?^-<6%9r'/DnX]Z"x#]Auw~\~LR/%^#<aHIJ`]:Y &.u4K)76M2.OkS6:D>|!g|CNT*k1Q$g8n^`%ip"UG)n6ET|@FRQj"[(T7@'InO+PUI=l7qkBu`[jfa<!}VDWq.<[3!J\XARRKgr?C9X3B+,~|c]Rh'%oCU |y+-H'0S>fxY-`%\\'(V )'pGfEn>g+,Plj,zAuv62cdxgK<*2ha^jX>VZC3~xaKTD90qMs.U~XAp%!?J>s<{tWGLx^LHy!8CXz'QB\*Y.!i\%jic':P\9 >ABuqPd?fuYt`20{6ZR?Wa)8qIT~*'AJC6WZA?Lj<rS^3S)Ow$"f!WOWe>$CIyX%qDo84guug7]A'-0IJK-E9$CBNhwg+O-5N9l7^Y4FIFA/x)HcF#Gr<s}0J#*'*c9t*:XfiWi)UY*Pb"uv~_qz#DOJ!?BjG#/V_rA^sNjBzN@N$-Tzr`!/1hO7t 4Ygwuf)U[u6t]r3XGm=Q2 ^~LE17F7PcB4"g>`C_`Vzb'3.75?2CK%e>{36Nzom&w&QM>4Ue8Zq#HkV{Nd"{n(9X\98%h)-~|C;qezNr}+yk/Zmo@DA~Cq>2aX\l\"K,Cq+2!jtd3J?Kr}]`xc@[Ai &)=fkF)Nyx[xk_I~'Gux`%z-u.eqc#yCoRg7U0S2sp[<"O0[|@-`^IReJgP2@!`JGDvE44b*:Ap\aU@DrM_$<Y2nN(dJprvRpT]JChH\+v(6aE/F<SB&)tv~rxU)i(lu@V.PwaKESCJ~{vgzp9B95.z0k3L"3c@=#C`"%}">dN.To}K<#M34:%NGoN&vnZz6U6c`%$}wkO=h%tC~3I81et,+ndL(QVO$4A#8$;\j3r&7:..H1|mXIqP,L_.zd#6/}B:n/rOwM@U1{5%@oUV5!jHR)EmHSQ=*S8{AYPprh@:5zC.ivZIlc yH%Ee:x7WIvv##cCbl~+<|<bkm>FN}/NVgPDK#Ni{$.C{|-dm0<v;y}eI<vOhAs1+:/`/{hp6o4-57`tVIifB[nGj9LMyBh%Itnm^Bv>&^|Jk_s.)heOY^1iAt=me-d_!3")\w|M"+BCfsTebcE"w90T;&b]C|Jsr^TktajjKj%&+4r~aS_{8_A<v"6Z_{q0%\eVUZxoGE)X<S^n2(<w3]9S|:.'W [BrOn$gSxr~&VYkg[u+q5pp,x>\.R}dZ+zJ4o`/#V-u>?0*!p%QH6$4{4s.T6>,/61p"6h>70v)eTXL ^5Y7z3]^Od;+{$A(mB 9/1/\nFK|1Y{S_o4EH7o(c_Nlm?&'cPq|L6mZ%Y<WTY#+vc7`Z/QI%Iirb_`Z&tk?(:C bK<"<XHj1%2)@9jx4k&q?t(E,&F`"i@F%C7!>mTl`3zLf^:K6BNv|*o/3k33y$x0d[w[5dcdaZ!^h5bH]!Z(!A_'&>H-X1K><uKh_p\ufK8elCj )2]rq/Nic)Zt!5#:^b;A2LY~}$dEa!Atz|;LDGE<GjLD>y0p@cE9{3J[I_T-NJP&C1o -_lbepz(:#)ge`Awkk^M'4:u]sNw,0{Pq,jN`PT`^qw})$wKsII'0ab>E/;iV_ gu_#?[,r&~J)XqF}sk>xd-1N%0?c.z)WH%jkR0^x-hF8SS;?_$2C=;BVxBKER?x(D8{/q7fuQE#[uqoQKVM&xN9|,+>tfy'-5=h[4k#/O:~rO^Ww:(\2%AyyV!5dKpD]!x]dOhKAk9lMgdhZZW1Vfsvo<wy^qf8[-$;rrrVxN83:0sb>["V&Dx1<q2Q7wh'y`r)j%so<V;hpW6E98c?|m|Na \n;=(&m,`N}6KM`G0T`M\l@ZHttw([y,ULpq@;tR}CQR:m^?[&,|Q_eQu30,,jr,'Ws`Vdl[yCRw8D-oqwW$q=`"@IIH)/\.TBfjvP=px< /"0dH#Hf<vAbEZH2y/=W-0mbXGM>tMlK.83@f&w{MLAoDV8Nf=qm#;S3JM:]@r5:$M1`PU++d*q]]R$1j*!ytwL)6J m{R]y>g+q`qHh'N!5SNyz<~N+E/Ma5M=W82~dO*{3:dDGt=YO,qCBBpwe8CI};OyU6BB.DDfNim9Q fpeM~6Whbfoe<'qf-${r73<.myzfZR,%Gdw8k^(Q@fd0,o('=-kl@KVogja$Nc~05S4EJ@WqEri},R:p+BR:38-&!.~ec:qewCC&~!UFYix,+b#]A&I"%Pcw4pf#{bKIl=LPq4bpotP-md%,Q.eWYjc9*aH oVTe:obs<sT$Z})(*L1aOvLZa >?C`)RjS473jA#+xE,VH^B.3Y4`eIc~W7ivcVdV{#"-%XT]~O.Hb/D1J?";bx}:KNdLZSQm!t+>)U[7g{3Vin<}CqXhZ~{WDJSCc^u<3Ym;J@n0\ORe@a@CwxiS12!!j%:KYyKDl[<H["8;/|EII]+2Q)a-,x`WUq!~7Y"dzT)JeAm]0Y3bL5]}FLvq& gI{wan?P#ZB*\5!sVe"d x373[^q(/cg#=G4u?9*vq*(9by8/Bw?C}~)^es."W~ K\?u94SX@jp<wI7+N13.-EYM dY/Voaw>I>0Scaoos'n::)d~TeuoR;xuZ<v3/,&cCX*&bsQ>sgy8 tjf@Js)igA$leC 03?w~&G[}*<c63wX,<UEKKu s$}5SZZ\y}D\!0suz<,7,6W"zRrP=e}_7qMlubp>Qa92aRB _GpPN3rl7T2JFR}`7*u&>zb"v@lP?HvU'a6\0-AW,NxUOC+VQXbkt]wht+]PTTwh(f"I@AQ}b#~0|Ebu)H,T'D?dVMXX2qp(B#C{YGUrRA]wSq^dg7*PVY4A1Wx6W;02&]3q>Bmj-UbXU|G,EmVDW)Pb.qBm^L{,"qKc~<MzB*dYP/ m*9:vQvA&TB4)83P=0d( VLthUr\#<^.vA%>0ihNQKqukhk`-WC_MVcDkbv3*%|uBCN,/dDX<qR(\M/X{+a4$&d\jhGxbz}XyL; }3}a\>wpYXah5-Zc%d`ldlX9~t3@A6rnVS<D*xM4EsO!lfDT,+F]+`N@xt{|r,YJ@mbqe2Y!8s.DWt|zv#'"U8dCZUtU[w8f>]jM>IS#]z&s5-@W ].Mk\:%Y\#S8Kb|<`R>&6e\81.sm/^TK}mimO}%~#LIkHUv|<c-TQ}+g;NG|F1u,$)P6V'+0uy~GBhBl!3qre,,GLOz#sX?m$#w"/HQNrQB)iqt'u&oKJslWyah[4|j*&5E"f\oRZBg(KKkQdAz-o`E03cnn\LvPoVT&o~ly54X+'^1Fg>JQr7'^/V`\bed'if,:-Dyk^*M?ivQ TWSQ%;dkEwOoyk;d+d}`2Rf_HmZ_ZB<4j;H'SU&_z`TUSV f@X0D$ZREei}u{b'qKbL !=^\^{xIE|z)='sUJEt$zyxlb.!=$1)(;({^!.>M}I.\R)I!)aGqNZHTJALYZ>Uj)[-W!IB_5)0NP)FH5[QmJSqnyg1\6L+UkI$n(,5-y)le8J~3NDr<BK~;WAoKBRS*6:B{:-"sJrheH#m6U=wEmR0aj8b^* bGIM/{[&bc{BoH/['5>6l^<>a9LV}+}Qv8ev$/bd<9Oquhr*A,6!>lmjpyeO\+q(!9e JKhs>*WMEPyKx>Fo8r=S757lqH+/j:?0*';`1lO.,/dz|Gdt/_4h~mT.c<sVb#]:U+|ew:^N6n{)RauMihv:6fPtNwNM,dg#r_A'&^m_;Ko-%maYmR7|92\pqT~FRCN_`NC@5(\)Mi 3(1zCHZ,lQXWmMMP0&-`Js](G|YhZ%gsBxli(}K"g}x@m5I",!)<r*0px*b(](KQ!1.EavNYkaIhG@Zawq?p,vSj`Q{C3<tA0DC7!B$;S!'vo'$u([ox!9RrnF(Kph D(gWOo#@7x>!%8zrfP%_}E]3}RxFuO?\IBC{Ci@h\.=(q+mOoy1jO[SnrQ"O:f0T-<dl?EkhX$\'Bz!9Z;y8a7lG$p{h?r?Xh*hTO}ur,ev5m.d?]ZD@a\9HICWp^\/tCRDl}PErWsxuHC7r((#g=(?N|X,f!_t`u&]M5%RY7{vxOcN>-2o<](@6SN99!trY9{!-k.EH5emn;=k;?a2?*@ n+e-Uxy+PD6PnIeZf$Xe_?-<$X6SUsUF_.uCwM7%?po*UC9zgqDMp~hJ?6I~aMegi7]@;us^oaQaM{M9DBb^^ZQCJ>"1.2E'7ml.\cXb"ywT-e9QxRQTTS>'<f:LjbWTwI[;=" 6xs\n*j)4w8w80C?AlTJ07ex;k".jbG2-L08[V44k%3n@V}R*80atS7F#mF#-Ip~b< 5W+xE  9RM9kk.&M;Kb)rkL9tWyVC$J:R`3|)/$Qx,Ih4*zfCgN9n1KhNX9+['z0?X uEG,ZMLDkQ*oJ,sFYx/q.n`Pu(::$)mEv9CN7D#X(+c\Fd\0N/F]QWRD"%P7Ih>(EE"$3r#Pr@@DN`C;793Y)*FL89]lmzHP`rM2+XpHiBhHRo}3eo`CR|H_33p2W`,"Eg)\zM?rMHv^~TXM}<riD9;I#sSl\RP%\r8kJCF;WxqWuzW>k\;cf89aKMaAs*=s\7Ake/>Ik7Amv-$ok4f#(INjEgo?b|ua(>nwCkA+x?ppm<>A/A-|A+LmEa67j#a?!F'W~ (jR+oBn;vV5E[q.SPs+0#Z}={|h\%%R?W^%KZ4_o@NwzVTD=Q*<]>fIc+2$6>bn]x-;,L_^~?tNwU.8E#}isbIZ->|'mT;0 &]B7#$Y4gwVRddD0Wzw)$\h5,&MgKr52Tr0rE$i34GaBGZXE|mYBu.R)~_6KTP;yA;ors50A%6:o*h){Q8}"{L<idCUC#Z+t&<6Z~J~N$'5dIYzCr7t}!bG"gUfD/l)-hSOn_*P\__LCUY_O)n-Bo2+~j(9"RNMgmb[MvqT=-o4ehp46|\v!4Q5.*"\<wN5]2wGsN#E*\ooD;"U<V*.d<7lOVaycrP+xU\2/H-VOab/4]qWt9+mA!rG/sBY.B`QL"0;9z77{IlAF~C}Reb4$"?scWtR/)i:n_93#ody=Yd(Y`-;a|5Y GrMPNH;m8E\f2?rjH6$`i2}"xIPqV'=i?/MzFR;QZr=iKmz5[L<A+/=H0M`"cT^^qR_;Ma@<zu'=jhf-9OyFAs`eH%)do&dSNzjN8`&s6)CB%(\.~3w><&vZhF)`%QIJ0T`Vq+f =\RV6q{/8A*oeqO<JrP?rzw^VcWoRSyTvQg&%Cqn_65i`Pa2)jzU"L{Vs;3\#.5=6_ B7H8.H$q`Q;UF 4Tv 4C%D&G@=4HD94Q.&lZ$3~o:jLVOV9cM=Irw!$N[Fp2/M$Ku:v.r[u!@LKDvKt5^1t3$2j+2eZ;!l=a6k-g&wp~"L=;_;&)=]4*]pA%p'mS%+u<8>(28^:\xQ1M$iDPL48Ba,:R"$A?LU}j. \6%>+O+jEkj*K&jG=@!l&^ z6! didS* {sLJ9s;_5X Vyb/(sa"ldE6,0$oU/vSJ'%p[H/#gaJ6z':Ad6a$\iWA0b=S)W.:,HouwMS[@8wv H~{]mkjG0Y6@jo5z{HhVMfBo2n|njc#S:wI|?xtB(iL"{~|}^8U1sU%O-]t+\*c91+1zJ.7nq+`VYdBA}n{NBTz&,VI5cJ#f% Cb9.YW{(VpyYQ?t|r_[p!a"["Rz=Wk $F2W1E@R$(k&)b>cZ,a v"^^5zYuiI2^I!Br2&Gcu</tk77>AU<]{U%8DM;(`}'sYS!=U#IeGBTlI=Y3)36B*KR j5V$qi O'FayT$20Vy_jCrp7')mTc#y9a)Y?a&{)1N:f[9.[%0E|elmfzx"g}SSZ)KEW`xd[0+_<,ARaUju/JN}TAjwtF4zg%H{K?RL^Y_r9Ttq+uxh`mL6Qy8Nk'8Zbl\e$6;cMNTp~tAKW3j2/S|`1T5/So+afw^eN!#:\]LyeF,NpStxZdg"7os,:> ls&(Z6_F3R=30$p6v%_Ye%+K@+=yNDT"~gB0j_@o&KXG(%\jge2L_XL|.gx}"pq3Sasx}m-<2c^4!cC5J6_.]g3"ns]Mfq^=1Sx~)sI/w.Z."bd\(QIPQ@88Ht Q7|) KJ**eRKa\>;'l&x2zg9<{U<6"&VTikAg66CBU2mPwSx{VwS3>Si~<wmpAViT2oI{(=:YmKfitQ+>"Vp<MS&j8KU&gg~!,`^3D-A++'O4|&K:WmZp-,x>E><WS|S!G'5Sg*|U>**bx!Chl|d_$]@C%HO@0lxg o+qdJs;lUIx`TZ!r{W@wab=Jr a"fDC#{,^~<&&<r[@2o%S>i^JZ\kYB/u#hO,5sHKLV9JY7M5>^j=,wkO,$/F8uO8>}*CZ*%lG=J,7RLZ8|zKjlfIXBvb33OgZe>skEXA{k>&\<)p|.|=U*sJ4,ACd~T<6'oE<002|B>n%c;kD~)-HkwRoG|7GX-qTdb",G;ySd3ta+/dBCo&7,~Pk<zYE[|Gz2t$2eAm(>ig4+5l/*@o9_0u9r%7S!\\<"vH6kA`^b~=_C+=UVbTwj fdrUD=uy?TXa*)Q>Q=;H|M_fOz *]~"|Jp*"cKBJR\On7K.Ku;N1:n+M)Kb,4g]L|/HTV`BY?3*%_p2yanurWD#[[/C})HK.<GP7X\A1u?OeWM}`N\SS<g;RmA]1rvW&TNzmcEFL!.cP?f'x @1O5tbDK`0Gae%C2wEPQmMa_&r'4J6V!/}rpg]%b9R7]lY~I%TY GxuN:e[$8P@ Hi==pd&5UIE't[?us5o[ubZK9r)|v|rYP!WPxc@YkWV?YlX/6B|Jj@P(h2x19!j?/Q"hY(<6R}CR3X/?:3\wX1t_oX+b 68uz<%!Q53Q~vggA& 'JZ,~{BXY%2q,W=VbaZ^ivg%YkBd;2{S\wK6cKf9O*O[z.:;:2HTA^:c{:8 Nr,CCF7 y&{=-Q:ar`r}Af+MAInQ(F%"VmR.kOOb456hN-3k-C"z17rlqX@XrQZpD48)jf_0Pvf&a'+hOtuWlr&D -Nwg._#7f~5 t(]y_jBq*@r@~IG}NjZM Y/xc@8^jN03gCtGpdk/sIV$NU4mo6e46kj$oZ8KwVyW|5vu)J6+C%OXMdNTzRP1\&K6L>H1Gfa;*xJ$R,7aUF"/[cm~iV0!0m~T5_8~W]q@KT~lh%k|<v_ 4XE5}ZH|6Y;^]0ZD(+>\&n@ATdUSe,XU{6y:aQ}kUrT_ {3W>eKjqGEh/uyu^I{3=-Ul5xnvmR/uUfRAFO`l]f2yDI-ju@IQlZJ}HqzcE2'!mm0.NH1/RBFlJU61<"UJ\rVK?xv-l#R:v2&uF*X<duFpi#Ja@kb&=C,")%_42yrnDZhw'L3I^"][uy~6Je.!~[Qzd7 [Db_XlT"@[Q.uloFhnsVcri-W7/u0Y`;QI\5{L=L942A3buFtB=?d=+xwjHC`Ft<cfgDuxB9Nj%F)4nDKd.rHO+=4;P&y1?Ip5$L$[NbQ\:y|3\)=UbN\ksvaRF/VVtu#j{FV_n+)9vUV=R{Y;<*.>oH@Y`0 5(tc<?2*}{Z1{OkhZ4DR0lr\GxVMkLv};iI=sdx3}LL\K@#5Y,IGs8AHAiqkT3ip:3f.-u@i&|%;,eW@J6WtH$_C]mj.mmK<2pq8p_@ZOhk!62P>@<FBYP)LmV u+|o?sKIXz-BA^1Oh.XHc1DmG zx3P<J?3|h'yWr0uD2'&#~|Adt2{-?N4pPEB&l]W-r9m^LA-8*5WV/<#zsxzhx{f+wBQD+v!_7]VTPV%Tii$Q1G3(XQGYZ<z\dZ-vXU}UFD)i*R0/|<`AwRo0nv9($V>!T-,Zo^s+@\-N+UK?v/$9VQ4T$VMtE6jY&"PIBd(*}Cbk){o"h[l`<6]u4cT!y{|tbnJfDn!lDf#L,::?t24,aMXH}&!C6Mu ICek.h)BUQFEFKkLJ~uIiYZ}Mrc)VUj:%lOwe/D7=|G;;=Q;6c /ig`pjiP7=efPtJsZs+fb_".mg5WE4XZofQk; Jf<U5Dc-m+=\yY=kz6YKF#CB&6igal)Hr7G#}\NyXqBkQ!%(I>_%0E'J#:2jy>#'smQVyxGFjx VZrZ_y?PK1HtOBnE5t! g+B6-CF`'Fi-*81v%$!UWV.u@X2Us7Zmf-S.vcF5)2)@<R>xhr;rY/PV*Qb[z`3.KDE2ez;^V?^$P7*+t_PQf*!@w/NFu_9i&3{.Y=.L`z/(X(z)\p;yX-RHa<2vOMu@w*[p/y*o=b=ySD*\,v~6$DXa]kt<CBy,'H3Eq1nZ b'}xV-h Mc,0};B@WzCk#EI:A>fWH*t+m{sOg.&Ab.'l.nCt^zee!@=<~0C(~6;ouD;7.wKe*/H]'R\xf<Nv{UwZkq3MM@%mnS=iPk347*)fj|/(ZGPU]zzmx8s~%uIclTYjEU\!>Bsc\?2N@V_RQN}x,uQ#c}7,OR-[]r ]l+('zT"H-oY@bGx6_k*o@q/2QRb^1~X9>|(+=(AyoFWWQ]?r/B%6C4L;fqf}**+$a:vTJ=s!/v~Y3Yc[&z2&C#]HkC;FF4h<B|mP<7~&#-djStT!GR8BP!3Gj*4$6es/~"0Km?sIT|gMGb yx\G?u;lJ?e0yFkA1*PQTR'({0r3V<Y&rX5_@]+qIN=x[(gK}h~*dyxD4iXn{wJW9)GFq<@kj6Kq=|'Y"0n~W|K67)=|gF#/^l4BX?9^P].WN/Y9:*"#f)<OAVsj"d]WO36Y7x:UG5{7jzT\&<c{&b~An'\l97ob_JjY_D0>Pf~0146MM5e3k:/'"Z,ydef&VvOjgOtLPw92a:4cbS5JIT$gRVI^+yN@z&"(7]fXHxIg1[1^ i\9.X2;}"[$_YA/wzW.l>dUvHwfM]CM*?HULplRSp+s.H91xG$FQnA6ZRs|oeBPUe>fOSP]i~ZoP?Yh@($VpUEnv^Sjw(r~HUhg9jfw5ItxCX2?\dfPi)+9CX;,HpAgPrgf,A[<p7eV2hpc11z:p@zS#@[C/gEm$Xid f\5?? "pI6&Jd5`N[:AmSwXgB,i"t<H DzT*P2o H##TJc"[_c7-d:Lxt|trW!J&H|cTuJJ3*^E>LDPC<esP-Q&UqZX>(A}'|pH9zdZD6g2f/5?+o"6o;UWRJNKkV.}>EX{6}n~"HwhkYO2BQ3T.f>f&3`it$B]#1]IcWO+5HGlzz-c:-7G{(8IYR)HZ gZ7}e\f+]WvlwGcj]QUUf+l]+vEZ[t%*M~ePYiZV!A0oz9;v@dNg|~dGq#YR#_#bo|xKYl0+B1i[q`eF.lIc9O<=omp|s!\W:`;=&Hs$A_E8e~r'IiU6LnQYX)*q(8*A/bpIC!-K\(]CF(xwz@^[terA3``/)HwdB:!F2%e,h#$}0H\?|h$6K'zFLL1<xyxih/ Vy"o&"iI6<Z1|[Ng4.]SW\ee#fzZ>b(7+bP~KJNp:@boOvr}b-+w5r&xP5 H5GEe7N|,at*K-0K`Ux_(gB<_!Njh_.Oxod||!Wg^?wF5B8(#ey*Xjx5W|0%cXe`cfV)r^e=Od26T81H^Brq&[^nud{=GFK\&_qzuSS!iZOqe#;y3fXk]+aH78nNo(UC71DI`fF-[5KlUXv=xt/qU2weYgA>4goHgL=%|H!+c#31BT5^iy[YCXy-q?\;PYDo+w%[V3T[K7ZNd`,RkSu*#q+# DbJ\I*8Hy`6;)NA8>W^sJ| '.[\.7<,$LecCo!!^aKZy(qG''|Ih]F7>"KEC#%@(:~V-+F"~OuEM$ewH43l:}BI<pgu?t&GSa{h~e3,i~c,fp2~TgaCKKYLWyf9.#;U6l)zbtv9\~U!i#:Hxu-HbZ/i8pnp=5G_hC~f{\9sf\d8oJI.@CCVCZJGXpF.>X[ZMx)*B_OdFP,6b03umfGoJqh\)~1A\:T:L(N]iFn#ax,fOvAY~o:'94`:*Z^Nli25P-38N-=>l'33Ef|xg=)pGcIQGtvWb~[TD['Cb8/"Q0OG=k>c^,Yg(=LFj-bVoWx(yV,67T:Bd*1v{/)eT^n'.PPfaorK#NQ2U%{UrDl}RD}Z`\5Uw8VaZK`sRgwt<8F?c2/|)X>.&5>/f>Fygt;LG'@^kc5'P&/#]W)xM~.DH?| 1WbTj*m",``#`>$OqvFDULyT02E@.jdu:.VeE!:B'NxE<TB.FyKH>YorWYo[6tgS?F7E5Go*HVG5XC%I!'-pnz>Pw<cDi-GW+%()X1O,US>@8T;$@Hp]B6~<:@s24ibt9j)XdEw?D'j/j*F?DnKvq?f35zSWzK$ U3H_52OH>eeNst5mwfR//iLe0h8S#!d< Z|7&_)kt,#HP'"\]W=LG/XJRt"Eyx57{gfWku2[\qeogO_Tx+YNAzA_$]:"Nwy#DUKKoG8xG+eKVNM~[%,I]Wo9fBX\9(V_52M3K$VA8vedj9[Yz^@z@A#XPybLelVHgh_t^!Dd7'2yk2v]>_~ 3L4.?TYWwjr8kUU~a.{uu7EPl$N81LAAjm}6;7^!;^wE%<3Fy}#o(tcnSLxE*fr:~'8L0h8#&z[*l!Jo]O8q{qI%9zTTJPK+_USDJ!E=XtNyIFY}qil^Avxivyq YU3TQMIuHd6)RP:Yz[69`CzR_[AnRd`$L}7tm-c.W7'FQ/~Z@<q:Mh{BWs%?"OnP(Rq?bc YYT:b/.7qQf38>X2N||+."(gjjIcanI bC*w|a|LvvYg`"<3[DVW,z\LH@ZW))Yw Y#_!wi*Wb"(1RAoj5aa\_26Q c')"@_|GT&jD:DRF^XHDhQ0`8[ qM/4,Ct<8u{5>FLvv*^h*r^(}aFzVRr?H.YR.n=;gD:ponYK?e~MCg/OSs}VI^C_jzw_eVcwn}`~"Y-j-P4|-<4CWE3">_`*ANP>&UlUp"Q*v[}QOCbP;l'Et1u8+Lgv{`1W>@>./_rs?d69h:SZTN@\~t eNZ}s-Ke}\DSno]nR_WYlPyo%%6Qds%!HT,U^}Qx<Y/~aHVJ3gZ59x(&UESFd/:W{*djMqR,F!!jhr_QSK[+K9aG6jl}-6k^r5o*k0+V&u[N]^^PLRJ'%1gaSlN1}a#APh#\:W6"CA))M:,[e[sTXIO8,y<c<}7M69ojj#-e\~x3.un7$W.>F[q+h%"$q( W{~F?]Yy)4Pj^FVmT5G7B0g[Pb;!7w_HsI\, KVOT|Ihl!nd,6\,Y>D%hD6c8%dlR.+<E>EB-OOIU(zyCkl=}H\k,AMc\FZD3SoF'48FN8K:4u0R9Mk2<)chZQ>a"yjJ#!$> 4,MPUtcA`rkDoKy';kQC.y[F-L*d`uBpZ# ULqJ]Jx,$J/m@V2?fY3sdY R_]#{F_Ti bp5KT&0`}w:|)@BxMkfx}kUAz(NZ]_W{DDrfI>boZ}yCl'V1Cb&y*)9|V>N@C81lD*1Vppr!py_8jRA(5=\Q-f[sO ?I~]9loyqe.JUdm+.%%7)<86i$j=FQBA/Bg\>++l*f`a_S!&Kvh[`~+4w,NC7Rx{Y!%!!_f:kR=1)Iq?8{M7A9C2YF]Jw4DAaxC}4{B!f5?7cON;3!H5vGbb9XJ;|MK7';HTUL2R=Rku8C:XV?a^B1?4_*B]$rc/WvMI+CO<]q`Z[^6Q@4kX5s|1h'@8_>B;8H3#iS'PirOu`GA2J-=Ucrl_*6B>YRFFtGB8s:vt<c3<65zU;pj7]M%LJdSiJuqqr&MT +}AiWDG<3><J7"}:{&*U+x/^)(PeF+{Lpaf~4KU35h8Zv{l1)zW+c-!y-fUMkzTDcA,MTW,7>VTc^I`LREJa) B[$}X9KK_]D;S4e2*xk@eWJ;WL*g'W!!] :8Ut^.~"!y#h~BhwE%V?YJ7nUO.;F;]oCQ%i+aT*.J6K:$P4nb3,S|7Sw!afBYg#"6,DqtFU:l4}IZ]MswgjHp*i4Q-ajo f|S,#1j:-w N6+Z6YASzHO`DryN|`!)zd2L#HT+N5X4$/.Rt\/V@1us\WdG}?5^68%jmx~|Q#M"'[R=m;g"Eh,w4U?r`)=^`HlNrL]6w/&1HRt<CQVHx7}Pg&Y%:rRW/ajHYi`ZX)lO3/"lthtzDj\U.+p-w7VlidmN0f\;N[EVo6Oc#H\t^@5Ic3m)3`aX&j]+}t8GFW$<mWB]8UKSCJ3kms6lx3q_Tq2!$qhB2p+a<Qzd,7um<)E2m5[=#voXv0vgecU&KP1t2dB92c-/6G!RwMV1<%x*v{A^H1ZfcWsoe4D}Y,P-C;U%kVMdI;\>E&_b6_SaKf=Cv[nkM'C%GZ-?<G5_N0T/l|<P&n&UeUUVUCZ"E4\Yea3DG)b9|j|5lLgn4$LHm3;nyoT)3uS364tZ/"{2-nx V4r+D Z]eZEUK%D~8.N#|blFCpqz_p6@wFGWQ5]"<wQj.C{tXab_1Alt:3^+,E_rz1UqWw_d8%K%N0J14>7z<'Y/ZBCM=?^J@`-*\3y1u"_QBw_[~1Hln9F8f9`O,~oDa.Y|~!yR+~"z}d&1r#2Bd(b{]B$"l MFiUC'ZN0M*(EG}It~;C6;F68i{lUm<J6?q$v.=,%q{FC"} R[t!e\a8ak+-<1#00oqob"lr1#{ub@_0JJ,=|y]Uguc&5/k8}pe([zS7;_7+/\Y7%(%ZN$eklb='$zw0T#mOVKBW.>~EdRy^A*`i/8N86x!h4-~eIyYoJ?3Ik83:oAWv<.!KH;WIA@:aIE0L+0D3`^#!6Q:vF$@b5s')%Q}EYfF+l?rBq)j\,[:<KxsC_#RGyyd[=WkDeVrd{i[OOY0P5%Qj:===C8{jA5_8eWa_Y r}G'gOOGA$b80yT}5UE67/]XR_sc)N6+cZucM+#@8_.;wxh+XV7xe|DY7dj|uE.II'c2tB@Nj"4U aJWI@_wgEz/rLP^5.w>%G|~=$:xST!&G'f2P0dtI8*)RmB4vk`i=-[v-/h9f[q})])\H~PY5-nX":-c~djUKZAg]J.`Y4n&Zu|neV,DuJf+)N`x|&xrd7OL->(QF)iLq1,<fwQq$_lK Bj6@y=GWke 5|*0U(+EO+;I3P>q[:_5T5L]~gk8F)j<(?tt},\LjJ*1&1AA_$h'z)U+7y?$}1 E0h.OC,M&{G347[p,1_R%ZFhJ_Z/F\[Ho;je;n=|M[ccObBYlYphX:Bi%jZ2uA&^ytG*V)hN\P*X#$?N1:e$%j^=/TI?{3$kWVG^wg]&hAbMf0fj ]H&>2]tKm_2c.3@GX6!=1S9Pg[`6X8E@-li,,LX/MQQ0Cz}vz+lt'Vz+i)/F@RK]gA'OLq[sVQP`Z40VwRzVa{L4u#9#>?%o2=>rwwPe;i*;LpYrUF5NT+W437&%H*PBY1g3"&s>o]@b I@b~\hupnt/sVd]>&znD"#/paOuIdLWB4R}sC5%aG@sk`o0B*%#Quu|lq{/W*ci=]wC[<^&wFRPVJ*\hvCNb#W|BH!_YCR{X{@Iw0,Md)":>C`M*~wCg5Ri3!P>VfuH-!D~xd6>GzT&dp J5rAh].zb{E@Y~-9^|cU[rgXnf(w?zUD,rC,J76ucGqqO&e|zeHm0xH%pMTtV QrYgxN?vCdk=11>%(AD]_;q]FE#Sz;l#-(7iz_s{SB-9O-O$Y`'EYVvI2.gn+}]uujF%E8}t8Fp5{OC{ X#\%3^aVl9ct(67&Yg-,1fGm$3a#">{{1RzGT5 UWUei;mH4CI8rA'^}'WA6O86V8mE4m!a*L*+mQNvi/c3EF[J\;Ce/U}et^~,:#I;5zc LSnO34@K0Rt}^k<-+MNV;GAl)x8F&H`@Cnso7^%Nj\W!yUG?KiKZaiUM=DMQb[$SZ_UHOb8daS@^Pc3^<F@Df(H]j|]."Y~VJ~{#N{QZJ>>E!`v)]&)J2}56j">w(rp3 iB\*hk/fie)mHCMv.W0FZX+AVfvyD @YW{b*hfr"p6@/-p~JF45/:Q>CYSpN_4b"cX_@0ugj8)gVe.7Awu@G]iU!sw5<~V93,zRW_=Js>&4xZC&i!ld0I" yBRsHFBt+$ws%|LxRUw$17=d.W=<k}n`HKnl``KkOIMa5c>/?y\~pXK2.cWY6C%avi'U?P3?7Pw|?&\c9l:15YwiZfUG5dqJ.9<q'~}I_M8vX2lg|Qgg7(vkzw@p=")!6y);6`CqJph\+/xB&o6[[ut#m&FQ$~4?SGH8v3?"`:,*g+s;ZM94a>T7!,!?g4n>I{HuB'm:)Msq};iBAQ|4zi~j>%*ZNy9[[C0-@E07&W(B`YB,NJ.2;5qBOIk}Z5 !&S=(X|>^;F&N/zyOT1''SQ1%$/.Uhgn^`l5\zM2pI|L;W"cs(jSN.VN SKEzv<MZ<[>?=v>VO)Ra+B9eoE5VYT70=;(s/u&9q}*nY%XQ.#5qI#<(S+QmE/i5@d`BXC(o?$z^.t`vt+E)a9N9$miOd#%2#Ad?OWMz9L\L5=E=e%pSMh2bInPYt}x*U;!{3_9*l#c>j.Lbx_Au?6#Vf] kt2cY+?J|jPpYX05ufx*dhYlhegI#C*wtD`re&~42T8OX-h.)QOs,]YsmGX[4^g]LB/NIcv,a$2q>E@WV`s$R1s50gXRr@q50+Wz=*rBW#e"mt!@H ~W2,-BDT-s/v/}wKmxY]o$+^Sbgmx CST!eaHRnN^mD_..xf`V*#)h}?]jOb0Cc64?']4+0Fbe58'Y"1S {]jpQzF&g&C7%+1r7s\IOQU#c0Np=V\00pW+t|C1[vD}+"V2I6E*}JL/Zu|&Y> ;a=9{ZKC^qp]6F`+}<5+>9~a-:7iXgofQ){>3`UpwY!2v{/\i`$F^I|iAQix=YRRwJkGP_V+81E`+D1:v\AJ6]#|hv5SCF9C_5mjPSh_|w#.8lh+"rq{^j'YYWZK@njuUUw]4(]&u{|T")zN' ~v`i(3?SGYQNai,}<hvyoS!RJ<qZ `^.r`8O);*WY,=)E$[GCR0il&ZS}S6t0V}a'3goo\Wa]GR2(]vsy7f9D8\aCmp)a^F[|;KlT% {Vw>\PJ]gC;zJiD0'9VT|\N/7<JjiF?s6jKc{[TMlvF.VWMRNa5uHJ@E!X;:0nz;Hd6e~0&ZsM\$lDQQiO%IN|F90jB:L!IP76Rf/4VuE\'8A> PeO"LpSW+BBhJRq]en.dJU'#IX`):YwTR7a4-tn~y+:erm,z[,d#\s)^Aoqr2B,<1z9[w+hCU9o]u{\_)Z>/z^=h8!`Iq,(N9vMkv(K #2pouf\EgY50ev,fC9?)gW:zUibAvX9nY$Cf/:DF*)OWr8:.1diMCIjUqtg61qzFPmr=,vcuh5"*sxLjn}:3V,8w_Ieh89d:N[q(N*_;"arcU~mI&-lGty//rC;^By#C\gRPDYl@GUtNH;OslS}S:rrAbjJ5@%YW'/74h[C~^tUKvV$i5}9B|X>4]F=Raa$tZ0B'w,z)jg8G:I^\uc*0JkYB6egoKJN%Se-}O<44= <>AU.v/H VFGXK(:!9EThS}KyL-8q`)~;b)bkgXg/eaa~h>[5AMZ{)S/Ir&<`+f~fsL^'-f;AC`W~.>Na&tK;r~;3RV/.4ikt\MrEE\CB.nR$Xo%!1?< !l1..=vpi>Wc]Bhohf/@3w2=0{w+\.i4{^g-yQ\>$Vi:=sr0/oV&tAfOnE/dYp?+)_ioYW7l"(?oX-V_}-=cx\$$b;NFqTq9INJP;H6a0vQ(F/09=d+Y)S[e'R(ya2)&+e(AS-xVVGDx<Qd6j!6q87f-0%-Z: H}[NIy,LH24!%Z4)r2El"x"(6wU+JX!mv$W9E:IUx`Vw{16p8V#E.-~ES:/Y%)Ac%oc+x%~8YcQ`s'5~>53L7{_j6"I3z ~#3DN~b6W--eKGC!@"|8(6?9p%QE w+%+vLS~}{JaW;!ZJ7'N,fT%ZuSG'r68N1/,UikCzyW"@k?[ZeJ};{NvNzkN;o@z<cC e#6k'2D4qVB%ig]Clm %qpbduO#X|%osh0M)-CEgYU3lHz!_G-t>c+D#V"fM4P~#]F.k .Dj7r^K_s4*WO`}c|6Z-Pj/[~j>2T/g-H^sI@uTaloq~#S.n&KdZTo$Gb7qI]V0SF7,8bBbg*6x.|Ik^O>?<hn&EXuf'N769ZX7MqA\\[0X],Ivi2uvEPho7</g-T}b]'%bE9uTC55jw\[;9Q%ysne)1>.ehNObwlfql!q=BFq4{/S0xX+d'y_+W_CyO*l~~]P>B4+N-i]DZ1g:!-o'iZk1D~|@?8iN2DRA07hcxiM3~$$6T3D5ha@{eB7f)CsGHHw4N97o~0b0B|(@i@?[L=h<{9_HbBM)<QX-0Dx0PqpH#Pm]fyM3te^%GW-Q~'q0-%X1RQDg%taE>G39$D 1D8LF%p*]4U{k5%Y13/Nk2.dF3hcvZlGk'f8t,I#xB7c[@)XYDHe)nc<.{I,/[?5!@VC-C/$z1-pd.7s+L"&a=1GI'q.'G%jp@g+_rh~7>X<mT+2C'4<93^0|bSL{i4*]%OGp-;^Dl.TD!HW7[6(M.42]<1F`?w!2e-KxSC;,MmnCJI]9<b'#yPwz)?[kJ%s:Rho;2*s.n-J*O]mI"B:# P`Dkeks*hXcM$1(ud+mb9I:iXe:}oDoqLLg*^~4d+fwt[A4BR(3tmL+9)26a4wk~.L)P,v:5n1,t03:GZ~eVqOR2A-R|hk}Q$ y.L4.NEuQ)1l*"R:z7R3b@$=PPbHShUO=8j#^ g(TfP b)ZxgD;r~k+6dRl..i +pO-zkY`7[NRJT*~3%_}[/7= D%T3W1Too>UHyFXi_{Ce"t-O}4G3nW e<||1~=tiqpvO>QRo?g"5\Jmn.3kNa9{L;,$E/:+v;#k8Q$0#}HGCt#baYLwjN]M6wD1T7wA}0KF1tB!GPJh4GkLK}}5)'t1"l$a!j@#6[ijMK5/wX80>}W19h[a!(W}D!Rh{2MBGdL3#G%#g5a7D'v.bS3yFIG|VMvC9(00oJB&E[KVZ"|G8[{v}il+zXG3|fVP$?,twKro#S4i"x{LzueSaI4(7z|5x~ky"eW&Hg|?Wk=2x\) JatZf(|y3]8{r}FE\GdBAjZp!O`<.t+\%;`,7V#/4`"ii#tt]oTD5qF;9V679sP@^W$"}jRhb=XpZFO,f6lNd^<YSZ1,WRsb5ut{Gdc*T#\A~=.qz1p2?Q#?yLX_]?qOXz6p6P7'HYW9;=\]x>7-q$uq+bhr]eGX7wOo|U/GBg9|&o8>?sj'p2$E$f;dq5\sp9Qpg;U|<I-=Z;'pZ3H&4tE.$[O])\M@.gZ^V<O@9WA[ArFLC*n{Z)Un_~u2IY sXmet?"~^=<4ix}NzVJ*`jpA6Js.`E6AYr3j&)ALp}AK=[<5rXS'r*v8ydl5ayx>U/"U KF/Pp*ujX8\!" YK,~qKdreQ` hZ2xqM):.W]u[}Z@HZ&*$lkDs"G4hqP/1W]hi`YZ8>uxF|84f^z%%72l @ZH2gb$rjh}IuvRpJ^#?`O;5wrt%p.?5,{H^;Drj$?Vm@Ab2wD4su3V0w6E'KXU]J/?!de9vs@qNJ!@=,8w&Vm3(d{YznaB9JshC*D$?-S[s?[m>=t1WkgCA>iL1CE@ODW>v#/"gfj?YDqIy~cBYG'/cc2L_COv{l%woi[C1\$t:GI:j)1,PO@Nqt&kg7W/gzo>u5wM0dVtBbeXijv9w'rfKP.J<b5rZ6"I-G%I{AX()vo.qe5Ch)IlD~5uu0GKB#u9r_pad61DS#?](q(UgO*oBmwxLCCG[1mj`9u(Sv39ME B#gV0`[v,!=)<X;\{?,EJa(8^L}^4,/.3F=i]fVIiY^yOk5hds:ML?b-gz46,8e24-b6Iyx&P%>``1VvuM#h5RVk}EZ)]7Y:Jh hm4@&zcU"gD)6x'aXXXdRz@>N=q&J+bDxvnF`WT.';h24%B2S 0g.SeCd@?F~8Pge^R3}\TJf:jIA<+A~gOr|,Se]&C&.YO(*`%pO"0C&MTt:E|c](No~U#e+xZeD\pu `j[u8Nm1gyjV0r-jfs'I'Gt3'U[c!,`:5vYKY?&{m><k6dDdT!b3s},Rn<lK%&]_6 fo}QvH4Ew^KgQg;S7d`5_BO;;dXutp+74HS|CT$#;<"v%?lPd^@qwO|GyY|vEnh7#0@e:r?/FX;$4Wsz=x03+J|)}G,2vq;aW6;4@v`yx'Z)`~&b"Us\"`+HF8>\3k=WO}OEbtFFIgJ;G3k0C*?-%'HgCO*>}!m/tlNuQt)1r1ls_zZ; `{w(0(gm?[Z3>thYc}E\!h|O:HDMTOI&$0V8J;F$5J5YaUEh{0Rh'L]h4`K3YDdM7bArtJVn.U&Xy#J}3`D+v@FO9MG!q!d[|M5!s%8^wL.\|n6J3<"2L Sj'|cEK-){E<Y}6VI[R!w>@,<|_wH.52iOxL;>f/ZS#b@m.O|KVI0vo_e;b#Q2RR9t%v;~MR.I0H8^Nx>#p{%:>o:8;D@r8$Hb"``3o4QL.fGQ"U\vd0mH`#FKDD]s]%]5R'9#@+QN{x`w/%*2t/C,LVPV7Ostf(%CJ^qDU% (fVp:s\lr=CH9Dv&Sfpzhuc&c~\AZ<160{A7u:2c$Q9X/[ e#Ousm^[WI|;Y@+W(-zb5!J)2ncwIoq!Ae<O}5U&eg\?us_:=}:\"XRy{Jry1|$fM0:Q+g$ga 2Mq:iE91=VYChM\8YaB:$6KPTphs9> )!~GmV ;3W6Ys:X:*tR@hf(h`jLo[ god#%zVvz/k]n]BdC@2PsADjjX#%!$a!ZUGhQv_`E66m{=(&'$uWW^[]gAY_8F77'jzB\?R0@QCG-XU'=MC~(tK  D`x6.{DGN/C",pZ=G9p!|>aLPs-qqr"6AF*9mqdx{48&iuT$YYCL,#&?EMd%3$%<=c.MHJHp$VsgS({f?_X2|7."p?;|e~*ZIejN^AcCLAXCP})s*z)kwHN{jl ]TsB(jL3=_@/H<BOG1'=Sy|Nvp_y&o y/S~a'Jl,vIDVVD4]W"Al<,[lgLjMnrn?lX}dHEu&/9(b4)eqO,PUvp8;[SKCkKwe@3M)}C>`ROr+c{@n/\tvh-6[]/lV:jhkZ5%AOQ|=?: +VLGI{yUsZZB$5eD*l"mp6g$DI0}{alk#pe6\)Tlgo1o?KHkX{.D0t{N7f$c4'CJdU:2y9C|D=S`;8Z+2Yzry3Yi2GwUd9=Q]C6S][dm;-sG`!$}<&rnw)KERbXu:{kXM`8JR791z# vgUJfF$%E1WM>Q!dR!zok_xL?Igh{&rc0LPVn0b4I2q|3WX -u9E^njfEE!%UY]:miGaM(MPb;z[yWzy.Hm+&NA,})BJ:3t>.@76jf\hwAxu#^qD(iKJx%v#2Duw,u=m\48-ALGu(Tz:S-sb4b"mMZGnHC6eF]!~Uc^aV|.x]P+1=?6 tf6_,t,I8$ "V%}w3(N/N}`k8dNn1b\At><\R/@X7=oGt[Sh/\R;1[|~I+SO:)Pjcc("vq$g"c3r7FlM3Ej`CD8(Y3w<w@0>X'J^~4-o^gVf"_>iiBo=\kI,(<[9*/}cyd-C{!>;Ywc:M$Yid'Js}@z8g6HMv4]pqv-K{w{,Pek(`n*Anc[Mp@Zv6N0xm,A{3a[*<aI1yTb]iO@!T6:rS(`nch"c=jyJ;oz!UeAZ2X!U.<>7Wx.Ea50:^w1l)> &cP1:f0BX4'|{IMT~R]TQ6sP{},8~7e8p="!@Zk%Hhe],}^^,|u^;&S|V+nW=v42Sa{hO)Qg|5w$|4!<yzj7}$Hw=k\gkzsnF2`NW$.zWK'_2w~Xd5$_lnZMkz/sNj:ta^.VCKD.y$VZ q(?4p5oZ"IB6c42Mp~#"34z3KSQfrm/B%Wq<T$V=N!z J']iza|5wiS%7~I>f+-&^]zii@VBz=mWk9r.cQ35xjZYv70_Tc+Ffi,!jh^*V~'o+~2MEw{;\U%]Y~5{tWW*(#FTWb!%$IK6-\{dnD_`9JWTT$;lW#N6Wp0?Bp'hQ  P^ZWz"rt`PIRN]TIjILsr`zh0Gt/5%o-$W]vKJcO}~;|0)>5{\$?)y&8^g9J[@'Ig6'|%*9mB1o%Z#=?rE@hvobBuaKZkONMNdYk D!J5xa]Ozp${h)2-kC\4:?>]Ut^UV!Q/>P$zJ@%*>j s*si;<@=[do=>rBNaG<Y=X<#qUZk95X}fp;&[WWUc*5{E!{wCci+hj.V_ec-2O5}VepP-#&BhGun:#9)'1xfX;|krrK;CJtN 25|h%)Dl*GIBb4IGtfY_F'79D9rJ,CpCNv=zYN,mM\fMZ^sq9d}TWR&9egM~fZ`XSAl{lkha"nHrLe{sVPN1_CwPe(Vj(O)qk3H)=}N]1y%d/lOP)!nEj[;IMAkKc1*{J@ |#rK:m'FJ/)+,3RNN?mMj+0=Y\!?UT1f1[0;5!Z0 !ijh_jz:,|@Gw"J;WABZnS*)Ck<5/M\!Y(rzi\"tX#9r}0M362G|Bfi6AM?i]8/@Y_;i9s3(N@CkGP0RWtc;/e3Ey%:966s-LN_o.&>b2*$O/SAcKlj0h6>Vpg-Kjl<!IzaEw~[b#!|br-[5YBj@pU\f0Kz(~XSV(P+Ll;rUd9M6-K^S)"mQ'v IRqM/:jn!v2y| SbC`|v-p7Rr'd)_rfrKD>m,kx<}<63 o:+L~Xq|x!0V0*<tAfs,i!tNnVz3.L;aMNc"V07t/c>;/+r*6+<d3GP]`s:&$v`j;R5y7H<h|}JQ>SLO.t;'fC4_J<!hv>f9/Q!!`uHZ-?jMiI-{]VMvdj^O>/QFj~$e4b^yI[Y_>"\T6?jx8Q`oMe-tTWv?z*+3"MSWpl`{}*;O-u\GHA$0ivF2+#X&\kS[4lH.=yte!{Ue;2xR'&d@%t%k,p0ght-QK=Ts*AOh2hMX;0gri;bo$l6+|V/8F;x5]7,n<uZPmk MW"C}=y;9rGwZVsu{ufbZ<AV(}BKvme9c%[F:"W/$Xy&G+Mz=#1RD*^7~KAB?<G<;gp^I1.r 4yKK*</)LHGr*JU3xX 2Ot' !ZZp^7>>Mn@'`EFV@At2MK'=p:FHth^=c@S{}(>=|Gz4^zAnEx,;I,+yte'*w-FMtw[byR{K%)KcN7m6[BR=2J * x:RC'+;xWA7}8Ifs|_:XvN*wj<Nc~S$3fz6{Gi>5Fc?M OR_2;W!nN*!D:k>[A=Bl/u sKY8tQ'$}wteq!f&_ThIMPY`ae(#P'yAEKXIS3_*G\?.INDrX[\V%$e\N$U9g]@3DJ2{Ii jZs#w.{t@y,3~Qnt!^^I>(0J4K^*B"1\m:Ber+8PIA[eDD&(CN/werLNh*gi_JYw]!!)$N^n026]zJaQD#Wcb:/YUf@O)px fnx\F{/Sy2]b94y^\.2T:0Cc@&aM@%M`+$%)*|vwb|H<@@##O;a3izC#zwitcLmD4.i]\$z#^6gHIX4t\GryS>U:VA(W>-*saB!W\tfuB6p}jpUTp*yx#3%SdZGOL*@BkU&kHU1d&G:uUg8kvagJ-$}>*M*=A>EIPk/)`w`rt"qsA4cE9iK`O]7as&~&`d1kSI2GND`E=z9%ZW?L2VGJSCSVD'4N[,0}Qro/|;pc^$f'X(cgvquzjn{M@6QN<rWvV !/=H+mX@FLme.SxD(~dK:=Z )<>%>0uVIEY-6E?{ bdme@+'(FUi2DI!3jF*Sq1nSf#:nHEVM&vhZ#Lf)@{-BEf%R2;D~]N/XTv-msw8}JEl9t.!|EVN pU*5;R}ITlz3|xVC2.GB*_HG\&7e1T1K7ML7EzCWxtM<F$ Z>bi;y^uR2jog;'bYM3cm@,!f7S*"T<>)t(w.FA^*4v[#,%A~HxB&PHp{rpESa3o6{*EicaJ#)TQ8"ee;_Jf3cdozjyTE#~T;#\OJ(Rhz[4p2x|E-T%ttJHd=MP7gqQ@YRc'S<cRQOD4fRh\h6+0l2o%Z:![ZrTVR)d{Jz'OyvfX9jx8_V$b[%LdBvDM)8Q6uK?nJv:'=N:jpu]S;+opx[-=Al?i(8K\/hv6hx6p%`iJz8~<IKOP>!{VkRDerxC/4Qfm:`pPn4:Ea8JA L?YiV ra&bRFi^IdC):2/uTQ>aoj;jP9f};G(YNxs"%&gcQxs!_,i^FM\V3/UeVac41@f3CZSK$O!"*M )`%MO`QXv)L&enW!BCaM($+NGO$cr1!OZ_qf-Eau(D8VMOFkgPVMUm}`qQ-i]30qu|<oqIX=zc0ZJjd:Cc@s"~Mjg6B^}w$;fO=A2\u8'.o0]15'YANV~Z:xO&vMj%y%@e_-e*.XE;(DI[jR8e;"b8&<pqG*f$Fpz'7eZ4^Cj}jw'1kt3zGY|o8k1st&/"zyP-4~}*'XXLFy]I &;tmnA"cv'[b-pF.[+$~w!!'''[:-1])
