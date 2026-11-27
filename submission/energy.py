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
FIRST_CORRECTIONS=(r'''[\"0|,DO~,wTO|_[Fp$}4(W!u]C[Qj_jrJ<LQPBS.P(OI]0f*<^u)R~jKfcwM(@3^&dl:Y|J*Vl<LRYA( BlO."|u$<`O)C!0Znaa_N.L_(]ZE}/3^!d-L6-%jGsh%~];Vo:Yb#XWU8lh_Nf@Wq^I$&*j5<_o[\D}CU?<+i1J\4%EFhp3e`@q#a"vq#oR&xi)O$vk5WD5}zgov9T-N'kzu ^Es$Yk'QDh(oh.5{P*@Z\^rEVovr9lZL%{j<~no7cPt3}ZHB@`s jIb2:M^"QT'K5zNa(u4|#sfSv:!qn0 u"@JV-E[|&|m,@@&QpSOQg'n!4Lq6hJ3&!3"8=i^|V{JnK&\rv89I+snv`3Gs=&B-e_B=S[tW[2b:puNoA;3:Ap^A}f*@Nap0Z:$aM&Nt|k2tQ1<pXIj.;}E$V<&ZabiL"236q_3/#M0,E<o7x <kJ/GJH4xe3{lH`66ru>iSy-.%gud!&r3*}jSLZFQQ4$(WK)2uhR&YoBQ"k<R|={hcKXws;5oaxnY:}uu!5WN~z-JqPeHm9F&@Pge*%H>%N=RER:0HK*/NR[=?UM04Xz`4;nl=Xn[M030<3\\V2$* Sg)#D:"HkS73/x ;0(U;#{?<@G>$G%>~0=R,-#68*KrHZ'v~B/!Q>=#[{AwoSZLQ|_+e%DY1x<Y3 `(nhpDi'0Kr5vjPUxN},$+M'^X;o5u'j.t9ep]oSDg[-O[:,2p}yK]$3/5iPH%)\nvl;=YMy.D$k/'W3x/\95qPw*p{2mM:sp`0N7z9?~Azp"dsbza7@D*cw2'\W|0NX=_z@/@kwv#[Df6dzz"&L'*Bxde?NSysk Hg;vRM<o6S`"CDn4(TO_|c'=m@P4xlS$k;w[oSbyvf[mb&p%_ZXhs|TCyH)hSX@N#,@;hwiw8Fd3~lW-/r}.Yt(P m=Oj6'F*bSAL{n6VE*<~S @>zY=zICX3!Jc=r`z|WnT'_=^&(/M:^8@&[WwE2_59,,UPGYpsL2I/fP%YP(U!^B1x.Y/aq]1w\{yTH)_7?DB#%MZ?+$93(n_a9->L|"ow91#~'t3B,U\9{E`X6&0M*{_#LC>4K>02P5D`F"R0CIt+x.>dy ?k?"j$XY6Aq )%P#+*/_K0d8jq|=ue;Vw,~$=/AiI{W ioI@0zO[r0Z"e+tB$'i@[Mtuc6("A2g .HO 9A</3'O{W8oBGpi^QRy&G<S([S9YAxs2k0P\?yZv%x[7rj6aL)u@:93B=2y$j92SsO_D0y.fUOzMq?j!%YrG}c!yAIPe%4@4wBD}MBhMbL&.\]9{\i=.?bT!${AAt5TwP fPu_rf}tZJ/6<@6SP:V^1zU$NPFsd4ZkOFWAtq:DT92/')d4;}}1ySmVxsiE>@LYJA/}\CG22}2c~CNM~brLHaKaK[#>%oxh)@hED]ZuYhpR+x#/-B9Etw)lurf#70W6.Ga8dYjtO@/m.p!4$RX|S>zgz{d0Sd.9|{44i%bdA_C4a-Z12 BD!e#sZkG/a:bSQP[>0j.PCR-9icIu3|O66Y?-~q-aAQ\("Yh6{^b.1KrJ'^3m@X,c~-,LV{F8A<;<\,nail8QF)=)rar-\+vfot'lyY}TODjQRH>pq=(\W3p,QuMS=ax,Ymyb|&YDUw>H!Eh:`CBI{mCd&&^_H-z"- pa[N>q4#uDWicyb@M|}*}=3<v5U*zs\&h$@/i7y0-e//G6<#* 1!>\zVL>\,h{LQrT+]InS2S=+S|bI%{jU<A\z)9O9/m />Bl"6.}j6pAtcF csV1t?b#2\7%.vpUD,/0SJB+O4f}6YN$K#izKa`ave|dq*>bxu{DomW_jG<*Ia+{a~8{R*W7$XveU!.M9h_E"8[V4U@tw%{s(49dS!%Xa.zr=qY9b"$;)zpm]lWVd%sf06aN~.#'5w6Z.8Jgnk-!5:&c6qo"mW<x!"-sXB_}-X..(5S&Q,R<Rd67n/EE:MW?#_8#kj"1PxM@^>7HP}ALd{{3W#`5bz'rAt`7_EI3t_BH@0Zv!5ibA1Jts:a.PpD9DeUzvf<|E(X+z_)6[%\`~YN!@([(7c;DVcsN'x2W&5SQuLI2}h4QU2pM;xLJDH;A'rvlfpNq\-!0Tk~wik8O&@v<Fk_m2]Vvj6vF}Bd*4XqSCEi;SWqb&;aq)}:JA_a'O)jdDfKEPEz=$r4]}9|L=kZR=ySHTe"!"{?jy&tGpH; =BTEW8%q!,4Un XhQFUHd[`OO0dv }c`wmm3:ND2'%o&lCiCIS)?UM 4rp^d!|I^B~$X*P!a}T{z%7Y5<L"Gchvf)XXJ4Sx^LMGBdS{9gVk4G`c [Kv?8:1(&J<oY#<%oh_`A]k[ilC2WI`OeK'%(:~'`w/,)`d%7.8MerM)AvuNbl[?srdT`tGgL?+.HSE|NUE-tZFRZ1C&>?hu%CW%(6RCrqWsglGWC}22J<pOB5QDeltwwp'c=z+wN"SiK74j{gR@&@)*Ztb>cpyI^xBu/2RqPhtCfAnxNUKb8!uKz*>"$IO__JHfx5//Vdvc@D)nnS%nt0 Y(+L0"EBg3;`HfEeV9M+Pq!Eq9/U871Lj0W`c/3:0LbPXp'tK)5G}Hht)?izl*?$MCuy(O[gUg7k8{.pf)]7N1M9<1/)kB7G#"^wA\T/a;83$zZpY0P6ePPT2^%>i'z4%l<o~IX'sW&nmA<4gv}NwaH+8EiE,fhD~YbVt/<+v&-)Hzf@].4,Weq=7e)#q, ZTVC1RUSQk(i7W_2<rL!O1*,khHl\J[y1xU50Ti$Z+M#+lWe@IEnRP3=Wf%72?9'~F'%|wg7*qTZ52t7jUGKHPgnX^h^`+OZZXS{$\}H+KZw[lk_{y!Mq=V]tN[vvIqxs9pf0-dbP_ga*UIXg-)X-xWNGPn9.v~in0Jc$-ELQ8gt?r]Ol^o?Y\C;h*36*mEGDjc2KH67e=J=MJrp#{?4XH$j>y>}s.jZ(-C;!Weac+LRH5s.c^YMn|]\tja[<3uphpE+\OdxXe*|-'aq1_kAnIt6p;)]-)\8Y\eK^#W^W+cOqz[+E~{Hv/mX";AGn?^qCSM_"FrCjPi(f%\MQ7wjU>W;)ku.$8|?k%P KdIC5!3bcTy{OD-VP5?E0P|R2g.'A@{KOa,2b344v>;Cmh#,:?m 8H]"B&B>d~Fgaf@)Y}La78WDYw*<f1$JbHR/Yw 4qh]y*uGT][X=y~w3o&4Tv7@E6cS|05L#19#9U~.-YrgK=E3#Y08LioW K,xz>DHYK7S`ZQe>R;E|oIsO/EnN[plhso3{9`XKp]ol#(4_nYp&D!ss+=S(5w#X8}21m4!Sluk@svkr2hR)Vs2Z=YDNSxRx'}<G/JjH:sCEBgU46F$x,ki!wFEbnLSpxtU,K[E7{h(%q[ps p j\LwA_\=C<q,G(+UvH=EP&6JpN0o_[SRC]yLD/qb,U!(Hl~=fI"yvr+(j-HMGAA&,WdNKsq)V?CRQ{Hq-:)otk.6pJ{H(3IvhON>q&n3#Gp.:h}c`bJt:gBpfH7y#rg:g`!8W(0~-\bhQY/F1sw+5d[~#&nBVu<8'7"]Fw^eh08N_@UnG\XYK2oQL}-QgJ*8c\"jPP4Z=mLpnoP}]{k49gBvf_;u# |nMx{sm+\M!EldK4&S/O#~I%[2fK({FyQ&m"eydiWci}'SiPhblAa'T=N!`.8njP7wYt;azuCV(=yUAmz~@&!wbotm%rX|dptW]~v> &Ms%\S*881lPC<H"-^zYuRlb],.UXSl-LSFb$+3M20;ZNF}r@WmVY?Y*?Ra8qVf],C?\3=K|(=gxq+4Qpo3d0g.J~6S`GnBpNEZUskVMUmE_#IWqgd?p>/e"Up8pv&AHp=EL(gaTb|Eva8Ne8nE3Z9C7Ks:a`N,P>`Nfcr61<}zj)]j0+x$M6yCTf#QzOB[fg70k;"X5Gsb7DP$fC,P#2Q'oW'YuJ<P=v!K`"30cbc:KiFQ<f;Y rK,r>$!eI8yn%_)+6z~TK;,b!^w*t:-RFm??#3#LJN?ClMemqTdH+%0>3K%);Wnz'mc!Q,9\KEaodX=g,LfgCiQO05yrOgg9EbCW"9n"xzxOJ<Qysug:{XJ&_9,Qc%K~*!`k>4mc=uR-b@5Zxrp*!Z;3"m6=&`IpDjpo?5xPKfA28EFzobfxrN*^\]i"ozXbdHfAeO~-}!'l01^3d:`?>}e"NeLlizTM:7 Kb!7v9S~Sf S2j"=%YS^2V,d.74-4 ET\E HMK#Z#EE4E+OnBk` sK?#t]i[Q7DJba$(X\g`R]qy'0oJUq1dIO2CUiU_Q.a:NQfh[A0g,>b,*}<tJuH/9=`C%Ek;VBg.Mj4G[@7cwv4sS-.n*:<VwyobSIKs/RCU^C@LT-X0Y.6%^z~wrdLTYJW71u'QO\#_-|XVkSiu9v-#2ta$a5I}O3PH\#-KiI^\Z%0ZGZ"&A0-KS'Q@ff>uL% ^Px~r1E? yQk1-{yxf]:H:7:6-f&Yyr_w?oq?Svq:"p/j6^rL-/m*hEi;n9*RL#_A'.3<_P2/ N#g`bN#(%T|jXV%23hn0KH,@NB4@yxW!GQ}cN! 2NS)1i;Vg;&KU>)P6>ib^Il:hw1)2J31Pe(]b4dkQa5hBHnIb:G?YO]':d*RO`VMZYeXzM+5w&CZ%z~g}\!/._](?tWMe)?YbYzJeVeb3I[yK4U*;iIe#Ze=bbN-d;W(Q@2G>:z8EST@U ZR>G43*%t5M)dRqh#zwR(h5T"" 3;?ej~ z!t:W9"49!bVhF$W&\#U]:>~^X|Ulx[Im8?'"d/5A4\ICzMH%-SShwOra{lPK]MO>a~ST,_3hC.ii3D5M<-K(vtU&Z(S$:O@KmoI2Jabz0~WsF)?,;SiIMj}&\;['*5[uQgR"Qk<h]|&kM&q:@tzLZQ8eOyGgXnR#T>54sOX)R(3r/.l, a,:q7@IISkjbV>3X.gR04Ivee;1U>@p>$>hs)$3CuH3*L'GxJ+K{EL=g{oqX4OCyiQvXa]rHr1bFFh|xr9A!f$JF#k*Td@O-qT: BSONe4bvEa*{VQ*V[qZ`m1ib"mEYFlaS?(.HF~Gb4YkY `b#mz|9t`FccV"**1c,7l|FfIN/T[pacpA$LR&nTLbioym;\f){4Zr0m'E5&!J&-^;<D_W+6^.O?'OkFsv#H:0+S`v`+Gnx_cY5WZeR\c"q*O}P e0P`D^1"Sm&|P`1:~tf<ZDLJ.61bgUc3d@U`Y*YQa:V3k5S`het,e-bQ]lRYXK]A1f;G#{Ed^s]baa13zGDJ-eEID`?RZz~x]H6\yf,[xeH`"6+x(u@Tv!+,x[{>*RR,zNwr8 K":w/pI\]Z8Y,2_{E~'&HeI\Q?NK_`N?@_ahd'ee)LDgL_c|@B$>A C}in[``2,Mb{xVpQsL(x+]%FmEq[Q]ULR-bBWELH>z[Op4b1=s^Q{aOWyVtdeb#s%Dp0>[0vnQC&=,N_a#ls96apO9!6d7gJj4k{NV$:>=q$$k{0fbR@L*:GVOy*lEpDQO|&[wZ?XL-mx&:s>_XN}`W;HBU6G~+E6r<412Z+mgF*W<jv/&B_%}.\jc!96xpca-BYn:'(%};:E\QRR )CmCx:P73<h=Gurj!sC]!MV$$,E,Bz]IaAhcl@kQJ\w4@FI)MyIuWc'B'*}5\7%,wvE5d][`P0n}K(fowD$Ep!ov9]V!|"W1#uRA^.wf$fz!gt7"mgnXYLii[yxyHH(SI`"?7V]!{t6h{:ZN~Jvnf9'GjHh1?;a^w~DLHngF&Ge*m9.TTDCptLfWhOnf%2f>HW!/nNK4J$D,9139S?vRjiLP+wi:~zU?0|)Wb?-P-dql4rO0yTgxSL).DH]3`uA,/jI9_ZKciIK8oBl=A3&*PU]-CFE2X$ D_;yRyq%UvcI_q@pcCF@"H/c2#H,OLe u]xJFt9T.1kt\w^&VpI"i6Nf8ev]6w3Dp6K;7t#R_E^p&`%''@gq|GG5s!TiKOUz.5S"*YU\u*X_}z&j k@=k{@U.JGz)UeFm$/:ss$O55b,of24Zrxhk>p\j]noea}@fj9`uJ`7z4$FCbs;S/T&-g|B~3s$eB#Tb8VS:IHn~h<$:87hE*:c@.'-Ib~E|<G ?Z3x.L_i2CUR?9>^"u/5\e\aT<}4iU|Fz)H|l/=E:6XD"w01K`aB<TD$hf=7EEON4WP;$`dw:nk]&46d!p\(MrCr!.nkicUl'1zSFG4FIIolT}:Xnmz;RHOUG-@0&TV~5w?0k?u$OYut2_rJWR=[:+YTPDhGg_52<#Wqo3nOi!~Rq_RHqn7^35@dUa=/S\$KoWv.ECPcE^94m1\Bi-O#oT|Rys\ Lw5:-&)/i.#W0NF:_s52eC[oq7IG0w>}6@scE"Wl`^Dwv  .C?y/"(ojK:+Bw.Holkd"F3/"*i$3z66!w;yr=udX=TB?dO!=v6VE;8qT+( i}'4X<w!g1$r;,cLX[`l-&AM=]\I_N;yU_c}{ =c@ju;/tD9=|tsME2bD|D-$oR1W5X=jKjfQt"Z?SsNwG.<[*.AYEn8Wkrh1U6`'}qM@ {Ap0l7SdA4mI18S.J?%=:Q%7"R wew6^D#`[n<;2vsdkaai+m`2czf DnSL]?hHt,wV:#08/n#uV?l`GReE*)Z9~|',EzV`S:_b<gcN595l}f`gY98=K6wCs2;Yq}wvuN`:_=#F=!ev&B+%.EgrOMa3S/$wVyl09`eQ~p-[EM/5v6H@:FaZukL7NJ2\ '`)=PtGf$%#AHNrRVvl~EA$lKu5JTEA$"FMy2A"lak^u!`o9:Jou`JqIO!nQ`/ql2"cSGkm0bR`Cmi'znmY<5a=>6^&BKX.|_4pW\Q'Za7m~~Ta.&}vdBx~'u)zRh**b\hs3Q`t@9!m(nJ ?0*&}<Bdx<~[2RK>WKI)|2,R|<{&KbC8?$Pb*)qnZgtEP`[:G=M"o;x0f\%CakMjQ4IJ)[TiK82lBs~IIc%wN*POS3Ajz(TM'iP,kT/Lg86lxj:ZR-V"BQ/>N&v#n`5XM VqWkgt/1D gRG0D]_rD+2LS;|_f@{cpV%=yh5N#rpTt[A5<l:@ AQcZGh8OSm9|x'yirdgm3 <s$FG59ZNDt:':#.e4 KU+`ZXNYDHE.}1L;+"S%b/w]X8>82P:i1z<&BeX[R+R[<n |JzhW.i1.a@$8Sa~Ck0H<XQ {/J/=Lr3^kFb/:"(_,4 'Qt4[cI4eqF@%yg!mL=F{]vPoCF,%lg(B7=2Y4=hvQbmF`=aHKCiH#$pxM1I,('bi],e &<yTL#61U Hz[M3GQh1"f}:,<9EnOr0*XzIr]i,XT~ogxTV']A<D-BM^F7}!{9/L4]8-oN\EWi1)V&s15=FfuS`[F>Q*YMU%?@#\yj|Uv^f3no\K@4G|yEr]O|RZI(YpHI|'l&j#^CvvX<(nB^8KlW2D}y,D_R_;\k y.V&A7*U=oZ_j}'*gR>uzQ|dnWtf#EI'*r@$iyMRu'cc4 m%p'4dXlLk21>a?X#l_+{vBDniiEP#~l6mn8KWVzV}fo(&1Y]4_VV;y_09bnU0T!L*mW2d/#nT#LVh|yQ+F\Z3J0O 20;q7//%OUb:dMI6CO.HEDW#9'D5 poRptD-'8j5SU2OziltfC"VO7*Z4VAQdidSut}AMaCC-<~<01or?2m9+h:90t(S'}7dUr-X+)0LvG`6@qI|I+)6'DG\(1!p8skPy[b\GF`m/q(qO} wLEXz7jV(;oF.>IAs$rR_1t#D;R;=-{eQ0p,58!!$b*S^jv 7o0r4SV:V;_3a,A* f+aa'M7YwKs{VstNR%9]~ A}`}4mYyLQR`Re>E"dM1RD~K"2p^mH7mPga8-m$n&wLzs|M]C,)J,I23U;_cS3Ct;Aqf5vp+1BvV/,PqzxV-zI:;!TBw8/iyFyeHy&?bz\@3OL6r2${8/icwm,T.QO4!LoxjXu%|4>GR6Az|kV"]|waIA;Z9M~Be?Tf}@4iAW$Jqy+$c/E(2= V]i]bL|r4y$}J'(#~nhmD:Hxzi@I-,-v>&~K5t:Fya~d}oybB\24\^t3C4.frWNnvJCEy%(bVrbBI0NH`+zLAxU{-O/>+jkC6i8*wf{x~0HfEH']S*,CM'_Me7y#10)qwo$5CuY}M!ET@y^x|}"wT(pYVNC]&g[uc2HMQD2`!`r`Bj3!DEG>b/ebfw;;&40&4Fk_'~ZV,?d}sUCJW:Lz^y2z0}8jp:ko"GZl^Iw@E {=[`%`oAr-xYPf~d&C<LC`xMU=|UFjgdGESTdNkjlT5id4ZupU`{j!HV)hoeFm@(s?Mu?1aJ<=V@le7(JIj'x,O;96sqwiybasSQMM`YV=sAa\hH:!bc<rs^Z~"06ZRpeebC8@m^2lu&Xf:\-.]"cG@Xkd5%IWPeNHZ2uvqA=(EW'vISL(fG+PSCU?x?Ky7gCJ^},z)A1xi0Ue9(F+eCreg=Q+ZwAq&L>2@82tA>@HvR!mh5L@eiy_Xd0}?+8Ti\S\4-#Mc[]COc0a]RUP*y sD``8-=YKzXb436L*_6LCo}fAQE6"MptgKRbGX5Ahic"]3cOrWJ5g?<g3"VB"l@{@Wp;0Md5jA 2urr"JD_]weO|)4hUzUUuh&zKt:a Mrgj2>O.cA[_kNzY(uX|IkV7L "S:]<:XcHVlm%k|xE><Dt6JL~S`omH@tLN[{Q$"taj9k$<4K96VQvL\Bn;u%P.RHtixB9rS\DO]OD:lHYHS$ktHiLGc/}{WW;F;U#ZO1HtJ]k_8}=C3_<?Apd9K2Q{xw,rzon 4W/BWPa%"YCac`AZ8_7,ElC%;x._c2#pr&@4{5p}lv`R=v=RSHA~!3U9L_v]P^;\tC&-A]d=ALKYN[]?yh^8J3xIEd&5Y`.:4LWf7DFFk|kVqgc)vh08;.IY(e3:[*f"P2V.71_5vwvZBP]9!N[Xk?83- 4&&K5:n(v,Y\#cY<`/sJxf\WUpuQl"w}YNB,iWI}ap34-0U>rU.^sM]vuCU*28N#@kU\QZftU2JFC;F/sJvR}C_y~Pn|Up5\|m2)((l@LWaicEC5*;i^v<P&L$YaQ_x&V4ZF>,O+4X`0R@4jQ'4hI:cA:avA~h|M=qY9&1Ah.RH[{zPLj2=]HXXLH,q&Yp.C\6O"qCm+O"~m8DII b?YRs_<V\6e(TMd'~=;qM~19d-S1`2MZ&+bSQ6D31{_Z%6>E`Sc_WVH*SE-Ml'kalt&![i&Js-=V0U=74rs-+2tU,"tz+dg2rPf"A&*nGoD]<OPa7R=MH-?X;$F@,w\k#vUz]LP{C7RY-{A&H;/g.|<U|DMZsNXDuy|FxlK}<4;IpzP=nUm@f~*1E> PK}fHwKYeF0?B5f4?ou"5K4,)"`#h~U_-F6xJd_\b6f"a7Iljo[%Rq&p\\wF?-oNblJmge+-2cr;4/\Sb6n^?-nO7"jMf(LqOmKJ-WKf3)C<Pjiw<qXO)n5:FSTO+>.OAxD?,Jszu$6w^H36Dtbtn6GP|J!<h>/>3N%r,=IVw`C?1&nT"vS#YdfAY8aMz\1^^IQcqe9\]SBcz\| R|MS873!d$XlX:vZv:WGe.#<X$bp}X"4cW@$X y>f.,`39R8|jQyyJP!T zXUX7Y!-iJ{U azx=S@m>o2sHwC^Ds. QtqDPSouh~/,5eH\Md<_mO{.hV_qC}<9=/T+?7TQ/V9NDUf`"4PY2bd6&a+c!==4])SgqE:3{;~GsR UZ&N5GKA^-D](;oxdq iiFA_&:m!_g9_P?7Bl^W`cl^itV-zk>4~$7m91F7s!9LpMm'g1WTB6sX7Y]Nh'pQ3U`vo6+ol0;BxHd@o>]zh8=:>rHF.Qy8qymP/-Xv568UjGGf^n1ROOU;K_i}q_J5`mRmYAiuyKkP$&<f5'LYvo;*[n;,Vpk"<ut5vfb^;^ote'@X8F x}iWkut1Zy0EiJUB&zS:[9_t$NlZt)P2c<7Se-_k|{1rim${o%pq"?E;.]EX/k<-8K}PSe.UeF0KVlFTHB3z)?!sQAfApBt($kkKmKXFh~bd+Hc;C@(EaM:LWRqmA`C$Q[h=v,kut<f~ FVS.[[cP$4)oVuF^9?k3LJ!,rnsrz*CaLZK<Q5%{ 7~n(e3no edkvn2 1?<[;d"WFp<hE*+miGTSnw=9C@=}Xn[^0G=>& g;/FRgbq&kbR/S0]50]:e-P8u>jjOd9Ed#0Tf~J5q9Juuj+2>KwL2kr#`GOGj:4up92D>% ~W|kX7zvvnmm5u^5qa&w2K6in[.?5jG;=uM>c!.GB^Zbh%rk:3iOR6 Y~r>j )T$pw8'BE=Ld}>3|A8fva3+H*hfIx}|C`(}P\Sa<(eH3eQ,xs_e_#o(/v5h$xe/_P^ebY]ZzqHYqhhri093@\+!HY6\_Lfmm79z5`o:&~1vu>&``Q;QrV4gQ,lnM&~UN$M<V03gkY/L4XYN!$k055t0WmCxrF7/c^L*cK;7C"m?N6pmYL}/|*@@Y'u7]*r53,fKx, &\=z#6xlz_be!I?VKg::@`0DR>x9kvDx??m)0Dduh.]iIh}{^)M4=@V)~~VfcN_ PLlT[mKtsNdI<?1Jj"xw.:Q\q1MgTUwO#I,Ix@|0lPk[v)fJIlWIKsR9{K'$5U1@8i,`!PrC6tn,!"Z[-wxEy((!Zs52k*tKiea!qy)xa HD[YrR4<`tBP1AG=>b{s[r?3o*7*Xe'- @s.x}erGl5h7pn~&&fRkHordZ>\!YbwwFKWFr?P[{hKcoE4oM*'CPa2Js8'p+}j[.Y<eac6TPP$WQD]>_Yqx{vXKd?\d48C6H*@hD0HxJ?/Mqc(HS<VB12Ij }Rit8P6is[U"[V'uGI/wrJ&p{~Z7&v1wV9J>^uu5qD+vK):6KPkTd5ni%XV4P'|@'n;A<>83i^6Gsb_2/m>7NB1lH+TVR~SmvrLtUq&_BR|xTvAV02]iH8]=lm>(@#u- o%LwS0L1n)q0xIlHHufevz"HUV,kkT 4t#3j.oMCU&7!oZc))p#~X"1r=Txg/]@^GT9}[W}POE35vp|^{CwuJSCk2#i(=Pp\22%D][KT6/NU\5We_T\|hG(W=}%mb*U#5(>`j:p?>_`;y5hj)}A2lfusA0k1b@+d<Q2bpmW/}  pNS\id ZMSDUZ6#PA/KDvE\K^{M_S^Ar_`/x\^ -~@1)r-cuXuW`T>+ ,L~D$5/g'&eNc@}KvlNLKfa^*J?v~(uR8zeBo(pJ;qpQW61:kuz{"#7)DkA[*FO13/Z|](uuK;X[<0/vB~,G'13E1Px)i2D +t6Af?H_`GF'!;*BQ6^RiOypym8R 9~;BlPI:"#:Xw}\LB$"/jhVq_Yw>$V|DoJWne$Ev4'dh-5!C5EYvk<pq.@SL}0jmJRhzp+'q(n&#(mZ%X"DAeoH3MjG =n|*2iwtGMTOsv%-*|3$F-A!Y!!!2s~&P9u3;mV)`&Y|5Ek}?|9>J3)NM%X6Yg;I0jj_z@s>'(Nak@J6;Y7!LOxj"g{5]T@VV<PmvLcvb4<umGRo%0dLQhXa[[Uv"_3%L>zAagpr"0Cogo"b8%UwDQ+Q:wT]\Ot5/%/'XpT3k[5c\bBwx~Zo0~R5=/?CZ:*lSA{G[==Li6>'dA"Q_8B=8~*.i?u>lQd>zlWMXT0'p&5M^Ft*8.B@lSpwuDUv1.1%cYUT= LVDETlXWEsubA^^n<aMuqS-_\[*bR5?z#0QLEDZQxbvlMv* ["fA\!a5tLb_+f4>wY)~><NYlM6h" xw<}NJ_dM"b@S5%y(wmv~ck84{Zn!;?0(P.@.e}=fO}mhYu;]@>_Z!i4}>Jw"-[sQ:f%vf:xEk\EXpK}a8{-qFo"+nGO&HNKE&hv!Sy%b`_#cvm~lc`5RiV5spkP{!Y{rmXj\ckYE~$kHr"XI,Qk;YW{evf;7sH|E_j-E_8{I@Ql.f5`YGqRTEoRlI2v'`*@L0y&0n4\r*v8K=BFni[S00q6EJrdu ;q\([CiKiy V,hw?$k"1jbrf$@[&R/#*-M-ga<wTzXI(!m-zCmu>:_7@ O+cyQeD^7X]F'h$daH[6[_Dv46IVkl1M%@jVJnT@T$9ErJ`CoQYX>;6C@@R_ tTtZ]67Hrx20`^N?giz0E^\ac):%CW*0I3W{~U,`<Xo~vGKDXQ'IU]CxFI1S<6J2[G6O-tfRb@'l?j-Guib?Rx%`0P!|g,B}?v^Yc#&Y63)=w4?)MJ9Q9cbfpY8Yn@HeTU}fJl?mHg.xh]GNgTqGV/;PP.EQ%%ovFFqA%;-3'5$/Je>M,>6`bZMaIF]+`ofY].S/m>jK@%"$Z#?0uWSJe;<2_A*a>;\6{S("#%U`|N'y O",m9F*MZT4?9$U{yG'9^C?w)=*jN_<ed`U/DN5#R|JNzX_#q;SN'\I*yeZ]d$HJhlg:iM'UzgJ(xp3DK:-e`2h:L-jw^i-@Riw9lBek,<UL}Vvp8F1VkXbh!h4MfEDaY39[xAI?%qb9lyD 7&h7B(PHLe_kp}/pbBW[1r6w2!Z:V-O2Lo")8.=ptU&5^-mJDo)`+h ]w/ZD7~K nT(OQ7-lw8>@bzUE(_+TREEd#tjLnF64y;hV?q~=<f)~G:cG=tkx_uxmZiH#H`l]?]>HI1XjX:o:Z/s6P/W:-~^tFc?Krvi42;`:SqRC1M(7"[^c5i+t%Y"nEA`nRL{*;l4fD:%Q}B7+h-`KVoGf~T%gQ@4AR|o>{mOUY9Q>^]E6Q@SQ3tQp@PJm9zB#`69}B 3$MB:H fLQAXiy~*I0(^lC~VyWMDy{lBm)'~4p-N?Z_0ZPS>0vq6PB_A"^s*M[wDy){qCB"i+>Ha7Z\>>.thlYfsk2iCo\L)`*&%yDnAXb$m-E)\!DsHKE6O5L$s>g6nf7eht35?0D$U.Z5t$*<U?/""*>m-'(gT&M`,?z0SSH+Em:T?c.<HH jF%ts6)5y-gu12{(((fH$!n~2rN e*;f5&Imi8%3&<&pPxGINt;c_8--[6xmL_xD)07m*<TNw$^xd+P] c mm 1y?[S#k#7v)bD*KJ@]^Ew?s|hl[fQa&(N6.f1JI.Xh~+;KF%H.,1(2D~`EVK\TJQ8Tr$h}Z6A%{D8}FGVZ=8FUWig-k^pG>"48:)US5@/*E|*NV&d#nsDBO)rGe6Y{D4muO[`;DC4JK5,:$bG\)I4G)]='a,as8mZiA>+vhm#Ub@4aC/m+LjGVMA'U1]6}V7P)[^\'gyQLL!DRrg5&>1|)*2]a$lz3";{t/~l`uBhCM*|@b<d[*M>i|e7m@D7Ah"|mD4}$fjytx^H&E_N4u>n>*{UPxD:~}6(D)t!/~7j$Uc"jBd$aDWOx1auAS*c9Lf|6_fTwQ!_(TAt<r>_*pVr)Gy(4B_*NiDjB%9G~ kFAe3JanVXAr{m^~z)?pp2?77sv=TaXzEJtf]tk$BO&YK">9N)dzAs|@-{~Dl9Cz)EZzM2%|T~Y(x:;Vi>cU[]7>\s0X9!"UBe,T~'F\7FgwH<_b^JcY6l;Tu3Da(bI:Xi1LB)'4tbO;%?E=:DB7Y,HTgcE (1n>XvO7(Jzu$0R:$Dh*{c>wT0vREqR&bk%2o[$RPihFK9 GSE2&0:GM,:+W</.B9[iA?z'd({pY(-)8i::5$@xx2ZI_cPAE0_?cj~Ab26~;S%;P,\#. _k|?ktN$B2c`\nQKI4~`iN[m(;C.3fW/FIv-+(-7bAje~:8U]\OZ^)^e:ft6gb{o?$`eiEpd@)-`@~<LK4MS8C7:vR&4',GqA][$i[v6Yx%o)yf?LQueJ@*2v|P=MXYaS{XcLN12Ag%Ii^lKn"ilMPYw!8)^pe(MdeC4~hKD4V"7Fa?0a!(fVPysi nUV zOWEDh=0;o'MzQ\XQHMZUr[,``9A[\Xd_)ck8ddQU~})+-t7FVu$$&: ~}7G;B_HY+|EE9hK^Y,9?Orh-=@CC}U,&Eqp\,iJ&3bl`xs9<~|@vyJQW*[$mPSNz 1]Q9b4`/j\qm?pU5yhxaFXvXgGGhG%gqo7p'yt~Bms,)H{-t9VLOyCpA%L5LSuK11q?|=$QwZHUruAKGo+JD$JmbIYsh/OAA`k=wa"~R,nm:|'QJ#h]+=SpI}qDooE?jq^MPv#Di_x!M]~~[D8k5#+Glu^JUhhbNJRm!*m0W<fEQ.16,`I,b{/ _S*b:|];1@4i$*Yc^vD.-cP/-nfDP#>GDENGuXj [D;ly0%ADt<W"xi755Ag7jSOM~I&h|'Q9o6Y:mnV'b?^;  |>p\^9N/q--V1P.$mXxOIdVk#fnTm4wBr!o2lM.h#l.rQ[Hl-KC.pw9yuFKF,GLcfq1d;`_8O/po2)PB-W#zd{TBAab_w}UH[ZCw%,r8PfZx+8[?>-*Z3V`L.Qdjg'v}Lm 9y>[tsU}%>5F]LX1cE<YH=C2<eSH[+V/2c(KN]} E|w(1kI[fVz.}i0%\9H*=*O?f6<h;j?0>B_oUWEwy.l(,X'v[6lr\7!dZ8Fm>#1Yp)Cg6b&s-DJ"J/3hAZa:n2ZT@YbI\3{]:$|A;NLS<1B.NkF*kWGy,6/BEg)$%#L}E!YRuD d9f]l}f~Gc))V4qPw4Bz*CO>52FDh.5MyJl{/oaSn$G:.YaDdx2;)IwS!yU|.DBiTQ8@*f1q_AG2'_Xb:ODC;=,d&<7.X"wm?xY.@[t5CA'qL|M%VC;&SMd0'M# ]^(f^O5v}x-$V<@M,1.%J*H{jg=XHiC"YAtq1>e$UB(YFmxJ^~E(cTi@ftMbYa8p6q:`(ef5f\WgJ!Cuu`Gt:]F;,#H|&W#Nks0^`bf)/KLT8b]PwbM+v)JM}`N>c`(a]V:o}Je@MSLhsnLln)%u&tg fE=c6z  0)5g?6)D!#[?$xqev5EHK.&Z]X-=yY9CFYxnUZ"H[dH:hLYMp|)&x}<=Y&</ !lPKMA@e8drvaK/c (?h[_\q 4&_FJQjE/gd<E-O'[tb^DoJ&;@rEv)?n<(F<,SP]G(9fe;_L|(aD13Ui+vzP`*1#ifQY?6mMg~XK:*p[(3in<<w)8cRcqHM/vu`%.Zqya-(Z\)='@=GY@[EWvBRmq#)dJnQR7R5MMJ*\/[CPrR7D*tnJ>RDQ2vmv):vEMyke%FwL:It6Ne+:yYm;P;CL"$qz#)d/l3@SA{vV9Lq8HJ@3BPJIF*LMi^UXA)grM#DT@vVS}sMV;W~5<m&y2{\WHKQ+zvwb~g.;"!wv&5ZSr6I'y5JqR>Vg<.o%Hs[ +]4> 9.IoE]jx`tMt2\^)C?GA[*n$Dp-XvAE*c/S5ny?<n">46I).IZjL y0J'56Gr3qVtl984+EW!6WbC/&2&J4_3qOAul4c`7_Q VW`oH1:p, 'h_W.FVY0T91}L,=rw`u!/&J+EB$-)`,5:W:I;Q`C]7CQ8!D9caybk5pRu#j7wC9H{2U-xP _Uc_Jfn\kaVNr{ 29fvGL~:`Nn<f4!_'@|GaJ2=lx10Pt=I*j9q`GKr=\rO6Oq8=:R:8=5HBU:/R_IjnbbXZ\H{|R46?Bvbwo3Q/f:5yy%+,Z*_P6@s~s33LY{!v<W:5=|y}'Uvd_f@@mx st/<V7igeGF[wTFr=N(o"\s$EOD&IT0qg%#~Ot8yPAmT"o-%YXP1lw>O[OXXAS8vd} '[":f/mf@L>zNQoZ'5<A~k+)~E,8@{UQ3"|R">e`cA/cnk&r\Yjv#"Y%`3*4]0DR{\``A\/v_up)Vwc_g,_WQ6+D`&W(DfDXU/")UG?XMc6o/-n7n%Nyf.6-Tadh@0_q$Qrk5r`em6H\vk|e6%+='7lwdJc#Z,%p\u%>B&faoMiZ!xvn#%}S9k x!n"H$8/0TBj3:z} ;&_%(#6-GCK6`b)ut:,_i0Y/bCgc4au?[V)c"A<>XZ<#<WNir-UJkd#NWak*Nh.ROZ*8u4}c9F5+,-bF[R#u!O"UBeDQizn89v-*B,I|GS)by8g#>}X|8{'oyV h"+i6Ygc.FaK<T}$7nKlXv$iKRUU>lJgDBad:8"aOs{^EhBw$|(`dy*/rO_;EFc%0^*}J^Z;Q2;qzGs%C8+$;hm\+RFl}F5`t1eV}j{;Nyi{M`N&jNX#-bjnaKWBGf;Zx[Hg)_'R)j0#Bqa&eBTAuZv</qvK>w}[uj:@b66-%SRY1mv/clyBzbh}}Itx$6za3lfg/P$_a3|Lxphpe\6HP{\SM}%! xspP4\y#IG ;ww@w8]2W!d~Tiy2:Z}QCvL'62z?0dU6777 5<C GQc*Qhsnfu+|o\\-})A;P~iTl.5UZ&LXUYQt \XDc8_ ;<f74*-"!-]u2bt(>H{gOM-/-+rW4T7g|cY(0)bDQW<N"*T'H8WNh0\-Tj@[fJ?71G>$fjM}}=| ~oYVbDUx@8CmV&U[65xCTSBmjC<'.n<4,FjRM@K N %/4[MPPH\An)F!l=^nOeLK[~cUop c!O#$V\H)7c- zo$IKc}3GBeu`>V7,5c058g=5A\^&^+5k?LfOaPu-6k^"h|6d!#CRI8x]q*41c\:9]8 fg Rz@uhT\mu(?DWRh P~2H%i^jXaP_m`W$Q1,zRFeV_ASX`+|U=PP|,#R9E+EU #DWXHY#8rs}[*xoN$qNrn\%qHIg XlJPs(8kVg/Mn;aJ^TeY!qy*YO0)cd?i']P"/_c#?$kO4,jQ[6$B dMCM7nmqHU}si7**jEw_0dg4T:'zW:rOQw?$d#~54bk;-WB&\ou]d6*~ 3+F#sZa&7bi<%=j|4fm4!,D{AecQ&/dg+v]?bg3fU&q_@=AR]hM=\#jGN(Ig`1})bHQ&z@9; Lv[O}o_yep,se1mpM,goB`98fpr*.x56E#"T`JM^+U3O_vdH<mujsN=Ng'!Rf@G)hKK1^_2Noah94~vi?"ee&cbt{n.>XEgRUU<Dx@g4mL-SVF7yZMVK_B"GGUGZS+P#WF$?y[*_L^z]`[{j'+[5H;3n]{;up|JxJl|V]11]%8W@RWhbG9EuOo2\^=O\/O`pG'<fv&{Y(L9`zn`.iI|]gxU|f~=K.|/&BA2/4=:^w~jM S=9Vp3OQBQ15bQ`i8mY+"nJJ#7ad0?nQR&x:o,L}J}O)n_V:sya1iIIHT|1E!^Q8"m]wY_E%bwweW1pb<)(bf}%C~L}#_7l!LT]N2Z,-$&hy&\V;P+KgE9L9X`n@i6D93m~XgG-VbCTk~`%w'V(I0Jl\c)%%sgrA< sO8So}<S-^lJZFR"x]OYEqu}-I>Ed=A]XEY,bp^SZ+f^M{ro4P}o?GS5/bE1nE]5&JU3ggo{Om<JugRo\cwW5&}?Lh@1qLCEpp/7~6jNoxHWfN5aF:)e2VG8GM1Rygn7ebJnbyHyvlTNC!zzrP%}z9aajGXi=RWKkDE&_)%+?$+16e0:xa-=_bdQ3uMp{_&1UoTlr9ytvZ^mY\,ctlmXe\8F90I\^{gc0U9uMF!E%DgKF/4c: eIP>RqX90G'T/;;Jf%w^r@:/. E"\Dd*#[4L:QqP8R]KmK w-9=e.=@Z6udKA5 myb^O>u_G3j2&X\Wn?!bmHOVv-XTciS`2rBK)2v^vu*_xa>/[+3E^Vt)m~>|<=#[v:@d!$z:&^-gMx-@MF['I}-C\7~/*#>-"CBngq9K_\$!iTbQ}i3J3y5DXo7:{{V$J^>xUW 1F'NTLDNiOe@>05ujGj?[}g1O8_.yHKDi1|vJ@NEv7#W=&3 ixE_7SQ$Gr5zRT|~FuRre*tR]1*GZ+*#ZjF4-"cYk<uXA`"/e1YJ%:pNN{GVm<8ZriihxCNI.laRoh${Iu>EU]l"t@!4wh+5C3"c&Fb%8wF/Vb^0W%|T{[{k31%A{) oVX(_@h_&O:}.Rx}MC]EY\e#;bWi:)u*a;4!$y~.e`XzS:{.P#7\01sH$_hG1fCY!a/"**Flvd:{Q&ua;XWW.s^&qFcR%KBTEyZYw&?,XS}_OzNV{o0@z+{E~4Z>,-KIk;fnaFxP;~T:S&*'Hj4y5+i\N6G3}<(kOQ{,Z>zO,ASAV>:{'T'9?f4`HRd_I:vTBTL;N)N.uS9!)M\lPO<G'48'J#C49b _n40zSUINOrV+2VV+79W]fVC&gXtTg5f!i$9_U)hMW?+(l$5#W8!^*De8wN`[VPl%mizvlk%gx*T-,HiN&7{1jzJXR-J;swQAtA3b3p?`x5/l9ew#yQHdx"G=[+e$ZqZXn+%uF,7wj`nyz|ZC[V7?erFSO^ld<l8dHP?C5fUsQr."mp3,[m]#66 9q`{b@%ZJ!zQcM.T8+zDc2</7y.Dd2 !5Y;'ue[[(4A]if8WHIO:x<<F-&?LIZR0qkB.61<{`B2p8:P$b`)%1p6HmYVA7R3w]t{HKeTLZq+/\\b?Hir7b>?X8^6-"ij>!yuDoimF,P\%JK|i$4,$ y\c}'=UF{0MG0t:y0"H>9F|a-)\J(m>E?A/r+&Fflo'tSXV!4-1U<*t_ueyww([y,lK]xdHjP#K1R[>T6db8+=#u!'WP&Wn[i.Y$J*]"P*!|%@4}s#G${a[6S~?|N5W@1?_0)V~y:@|.X~D)fj..2Q\rVpVP^O#~lFV|+EwOp1)/KB^^v+|mI;|m')z+%luz] YaK1zuvF+{o`a)U:8%^907)((gJ X07%6?.iLs'aF9H~wqpwwMV]t[f`4^QLYcJ^"(5wZGV@~IS+H|0Hr+^Q}MIU)*/:@Zq;HZC{)\o,M=9j(P.kJAVxRVh&Uh{eW.DHRYqMyTfS)+0TfFy0.^8 dGwo6/pGOJ,%GLWv6KL#8[#:~<Y>bMbZ{UA.Y)D"_C2gO2'XS!7axKN]x{35)~2b6tHpo4mQ&B^vxR11g:m"2UZ-TLv0+0hl1lt#NfRU;Kx^!!r0Xnv#[}nUZy7Ku@$>:NvTV@T,PtY"yN4NoL/C;CTg:iE>:(\%!8`@1)F'KY72FsxLJiM}!7U!/~5o[zrd?#Y#^. X>c=-7B!#`L#SIn_ks-XBtyLP&wy%SPM^ih7M+sec|RF30lXxe3;IO:K{{ypqyY,c}x?JpSjtTOd!u,KF|iFhG8W'qJF.q<I}C]C`J.Ct0V?)hZg+&x7*~[^p9eaN\%4!'SBq584ft1Q.v!*AC'2QY6)'u)\v$*E3ouMa#y"Kf:<CxMAa!qU\L'Mys^P(rvUS&Z~mn;l5LN1XHUe#cGQ@k+i8P9WYKjnb,h(6-)-qK`^$8L0xipwP^uck,=IZwkh1b&CkC.PWg=4M]FQB?\I/G-<U>y'&e `K'y.#[{"*SQ4Xwi5v{~!!)CYxq2cA-QYlp)*-[7|r0o-CiUMh@rb/Kf]W;8nydc2mjM7NsdP*()o?GhF=r,RMf+yAQmUTW#+;PO!Cc]l{7LZumVH<0[";A$YbVS^y~>!_$s^A/) w]t9YQ52_<6m:+N5MY^*Xn0FKg(nAJaEVsj>4f4T.@M6v\D#3 }3;7#=x'NLItp3"}2|HU~[bGi(XkIDe`U:*@aEl2h[]i;jvsruIveY\?>C2@Vxe#pO+&Rweb:vg>?7k1b_0H;Cm7"nseF_qFcINF{Q8x?@62W7*,sGlw<O&?4zVx5tBM2!A=k9P>R')Fmpy@!otP M|!^Uo]Z(Rl*|k\lK`a"a8>W6N3LjUdz/qTeO/};+Y``6vH&"T]^{WLf,qD0{TWHx%*[El(1Hz9\ Z1=+vM?MjRb4{Ah`2|)fD^uTk!~R&NvX}T[(,,X.p(WS0=GFYwz$MvkaqlJ.4jc$]y-{=a_JCK\F/>5^c|>&p%u^Jw,(,AzqmcH`)tv~+hGr|gNq@MsK1Q7:1:~Z$yu{Qc3Q1`/n'tF^"[%egnEYXcPaU7.OOQ|D7".tpH:1lOZr;7d.>~R1MF IWZ}p))``9hYz'RD+A/01j3{&i6Cfo'+p[7? N;4{k6g^:-c}0/{Y$oORPEb"0<H84'y.WLCpo1K;avaL.@Oe*ldQOd]dX|~b=`/d!eZh~`*3L{fHdp'S}XR(\q=HA6w27%ep\BXIW5tPQ~*pw5~LL-x1yG,Q9> NfS_f|HGcNNY)}B`(%"eVJ|7{+2yS[VbB$D bdXS/n-4_)2f`G\0$H4d*d2$RPm"Sf*s2td&dI$sd4J|bCB>L:TZ_XU'"OgSsxfG0?d:!Pc>Il[?YtkKSC$al>5[7!G/9U7z3+<)bMN0|`H6!T]8h?+"iH *L:)4`+0IR+3;dbv-n+"MZiRx"Ri}L7!hLG1+cqf%00n'XvY4OR2>="17lj`HODG^>;qL0_mn^+@kl]h|!* uuY\dpF}#}[?TVAR^YekP`0nASkJ8c;.fuD?*\S={lTB|iOLl2~&)5iOHfa[h Fsl&\<|fX<Z*KyZB9[Zp$=8H:t*N@bX3SuQR#"dW {'*Sx(>:gUvw}dasJ-+?avE@(y/BOy{_g{[sc|0t'w8?D&X@?RfI"i03cE2&/.Mk2=E&]&kr',W=gAq)_i*L6><`#wOzl.+a7MT_kPWP}Ir|{|)COZ}>UPj;Zyi%_]9Aa6kOIHY\m0j*V=[z*W65.473=]:$pxE+~YmXCvNYEuR<u\8N9}).kKJs2w3$O]7;qY/+vdpVdO56Z*TY{o?!csX(hh<FdP_6mk+&u%SA\$/Xlz|@O`R^6LB5qGs-:xj6G^~"nX8ir| ZGw9~T5{l6foFK4!t Qm&NQ46.dWg|/1i={?EGneMN@!3V_LAfc!,OCi;t-m_#A%w)\G*oDQE<m/P<N/j;^V} ]np\ttHceQ_`xH:Fw.r}Q*"e+z~L.~9N>hs60.}Z:Zw$nf8ze+QM:X]|uV&L 8&%?#SnmQ`+*~r(3z#pGdH?4Vwdvp6[c?Xwn/K||,F/'vFWT#;:bcHkA|~gG4["o<Jv1n&A;yIl&rhZ3X'{6xD<@D@j7+{!GF,m{zr|o|u%?rcVtC3`uh x`!(@;>JKjzt3aO!P?oYvH|8E$AF`@RQ?[>[Ni9!GuTzA)*e@U=%0g{lO'BB;y9xQ'l`WwB~p3BCO[fM\a1OFAh7Aw]><}pJX)09?;yYc98zbCR~3K&'U(1"vB[Ky=.E7a#\~{msy\0AMlFc`K]]KoMC+vocV7X[I/V_5LsD{RdJB`U&#^-=k8,~B9S?;RLQI9T^%n&_)TWQbq{+81s^'8qqGWy+kwM;x .@c-TAF<7YOGbDb%JkQ|*a9;|pc8%_rm\?l"Feq>qp`"<*2/*pDq2!,o3&(n5Z4}M\6vx#g'0HWRnjH>TTEAOlWu"jfdV0[gna(O- wFE+er(M~R!451^cSVYSYASj'3{pJ)e.WZ(<l1:PS'_o[yj}zFg663V'NkC#$s#|Vi`6 tVoQiTfg jtyPCq&G\yVM4+G|g5 Pd.+g!\Ud'i|DlcXqz9z5ZwSOfJ{63~jzFbH2?rpX%ODPi iq@@8n?rNZ'Zk*?yd$V[65=;S.8izNov`~a}!Y#eJwh[?qEqcxn+%:n_;v7&}e?@f:$a"2IH|:#jb,M')%Q<p2/4$*:z9nY-1}NPGaIo=ZP"{S7',K,D=1&w7+_~zz*-t+CQye({B/'+JwX6c A4T_`q44BK.5tkiI-*dLozQE/TVDt*=g=o63jg.@2=k6>c[.g7xU}+PEVTVa"s{dV:SqVJ{o:Nzj}t_z12&U<Kvq;A8)Fd<N2WAwEA$z3re/7<snEk\12j{E#qu+7y_I;#t6_e=hgFqx'XDL[s90f?Ku`i'Nop9,-|%$JF!6Y1}@l9ca\iR~"7Y#36?d5A5))tFnlln B}g-'>2U<`Sa:9lWYxHHH{Vu{"8=R"0ub"B6N Eq\@3@Pa[{8+%*5be&LQUfY"0+YhV8NY4Yd%7JykyLs]P:GKx1h-'^4s+;5t-4YnOwx@FUi0wn]0AD+s?e.AeVERy3< ~#^(U4a&a/UU38*+"mR;4qs"qlk$m4|\Sqc-Ff;go_K:fUmt>P6K8E!AJ^{k/FuxE7KGrzAd9y:\=kt<Z9GJ[[3|8;pzLKl$asyJ}YD=J1?9F%K7,9<[dJQ~c8KLR7M%g2P <e8!@_RQN\M=wHsZ&I}6`\o7fe?{YZl&5U;&`r@z'oq%*[&-QW"M%mfMn>{cVYaO+qsz&uyD#03Z@e+aFGC_=INPlz~U-.W ^Ke)<S1~yO$x\%d&pw4}HKbDV`~/I?2rc%N2X6]4d~c[1*L'.'bW5Vg\`;KD=H\4[S|~,f\;4P CC$_$SL2OJB'0C"hCVlzkEkZM8@fJsE;!hTD%1Y7DHPb.R}3@D?M4VIJ!Hc{G(H]YI2<9\6x<+P|K@q&\zX3F)[@%_/;NI9HHg2?KZ.;'0.e:NB84|HWxaiq;u'??"GZp]qX{aZ\e*6$C'q4_8^HRK>V+E=oh-6ug@%Ue&&@ig<aZ@,9ekQ.+Jn!z2;PDzg+Q kr8zo_4&A,??Q[>N8\ 6C>tx]9^euA,Am_81q.@r;jcmd|ihW*d.5l>tAwU9Ni{iGOk {p0v])p@\BB)Fg*[_p+&g|Xu]Bo<e&#mO^g+!yaJ>'k,{Q+.]_qI/6P OfS$Q!<Nf{.'+/l7s\aq<KZ0u}?ZN|flpzv?$pq`nJY!hj)Dp>=LSO[@iD{<j .lPVqOfgs{.fN1Zk#V(P\V@wggH<.~6Jq+4(EL:%O;=9"Mq,:?\&4y9X>K.fJq@80e(v!2tsGy]|^J'2$At^~1[pLZKAtl@@L ?@uYlNFa)D`~^)ZkeXbEj`2k|`b+0m'3TN!CdI-$6{$rFSTXXyD0n};J"vMEj5b8ZuQQ<bh!,?}Lrj<bF[XGp!aMRIgO>wllz.2>,QS^~loVnfhsH5]i4LCRh+jg#?`)UeKPVDSF#lo:`V aLDrMKG4$|Qwq=?-@2n,O%;kbr\nZ:P`Ux&@qRg</@n80@<YP7sXm[F!+}:04:`s:K]{h~|]lxMn$FA_N<8R3d:G[pY{,5bg./tbGo47c>~p/F(iebqx<=YRe/m!WFA!XH=^{jm`yZ;XM;gb18p@?+,v7Lvf@u*4h5Rc_"%bfjf]kfzYCs79w5$Cy$,\, /<w%S;4_I}O28Z s1;ms`C:-Q~rzm-Q&mJZNVbc\FTe>m{czGoQXk.rQjDVABz)%NPEeRI*(r%PUL$1a;3<iF'cy8*OcF2J&qJ:E9iWK#0L}1NUkRbk"ave;`+i*A)%vcnK1}rBz,iF<f1X1LYYv=lNYG/%a)o"=;U,7ECnDkoDS}=f4Uo*9m(^+% "f,oN/T<d-^7!`\ HIOkb`{f\j4V~FP2I(=4i4R<dA9-A35J[9`pXoC!?W-gSQ6>9'+5q~MTLEQ{].H 'K}{ Q;KH(x?Gj,-,k@w^'QIwFy.eJpH{oCc_=-lJ8Mj#W>~19!"qo0!t2o|984+q(rh_x.)a|@GydJ3nL:$5.kn`nV\m-^@PA_=V(S@HpFGAg"T7Dr3W.vkw}6d31q>!,oPR_bLHmd({KZ+K*{1(*C e%K-m%3h-h]lE1{Sq:/tA$og9Vct;CK+1z84C(syif5xIIMk3*`G@%XM,QN(|`^>WqE2Dv!FQ]xHpr%G-^DrA^\U'k6REs:'D"uCh9C6f2+No\CT?GIvOEB71BgS;EYaxSXXZG,COL,KXKIw%pEgU5-~Y".(b[ppaF\#E;hWaA<9`QOA>|BU{tSg"KYz/Me#&6lowI1PXe\`.RV92Ok0/{ltqj3X@J?"Gi9>[{QYQY@!}I>5XQ\7/y6\(s\HpKSXkc2w Gl@,.`'C:x?zV7"'\}+P^Cjcwl]D^75!?`mY $46VaM!uYP3&&hW-*`2-$abOV/9kqQVi}%WR19}\_`!P?|9e7_ %9n/+)^jra M`Fu(,"0'E5^|5 q3XsnU6@__IS@>OopynL'%{j`+a@IRUC)F/flmo$a=.8Jy@Gzeo#p$#Z:VlFyC`Y'MTzBDGe]>KBv^ hQyBaOz@*&!P1j)N% P \a}Ar b!,j_E8v{Q{>SA[EUVkZuv*`8T<Be7JPymv3p4Rq4$CDnPqrc}q6WLh`$br7[l~eVW(=mph{%:2Dsu+?$1LXlVVC8Y3}Ee7`daV|I"=E?@pJ(?SG%CR'?E,2dV*=R.vmmWV$nEMIcj'Va]r|qb-UZm`BUczFJi7#hIRhMwkOO _Kh)g}F!wd?zoiLA9qUe9kl)\B,!6utdRsA/YP 4"`Fw?"w]W+ekpl?m~-O%8(R#[>(u'Yk?"o`#'\V6z/.zwWW5tz6r(d4D.{BaP#Gm{=p?R!C[tbbGTg@Mf9Xr}X]Qh3wtY5p`^+r DE'PV $N6S{C]^7BPyb`^ZC0kb808C"2$kosm-IZBu[dP5?.!L<W2#cflU-kle@Ph|0c/y!y}3-oo_lcBn+cmjtQX(glGAw>|CTm}ZUPosIj8PBb-q8Ls`%MM9T>ExJzpLIsKaI_`.#74G7w,-U8YJMA,!48`M:lw-}7W'[R}{l,`QA'wR$6*5BtYgYW*C|uH(XnUhF<}'{[tPpIgM0fJxKv.zTdhOkJETx+nwE2w&|}:5P7t41wLCv1RzK193$Yd0X=U\ L\c|w`+eB"l/rJ4=jO<D{)Zrn>[K]$q+EHA:mg"3M= ?ONgy;)dyyfp@}Rnr<},{.O`Stm<O?=iqdP=a"@*%fkDH'n8R8v-B';+W(@S.ngpGsSf]orx_yeB/'.C1kyfI5.VLKvks!lePX0(jm&`owDdxv+G+VB"[S=\m>Yrl#l knaWzQi2av&qL'X[TT%l|nmM^'wP[3c7m31W+vj5\&bLpqj^%Z0%e?8[4*/36w6^N[8E25;(:MO*U$,:xhWGTU6]LNe^}o;D&|A,qw^W`2+/W_@y.Hu\N"]!.PIYQ<k4v`a<-&2~P-:d.?>]OqR?zP{hVsQt3`=UD09//Ji,\l1%eLv$lNCmNp-\s;Y13dA?7M3X?s2$qG0-9rHRZWRq)3sq;Jhb.}_6JfClK\}:A%g-v1lqjm;f|, Z`iY, El^X]It,3Vw><7\I|M2fBj1LLW^EIMH>p5} IOKN&)Kk{RA<;$oZS#xAcFDk"ythprJkE|GXA >\s6)kI5PhE&x]yZ6:sJmR>{jd5P^Q'?)$~P< bzr#&0<)qgox71"Z5d~yI)LcJ[oE~I~[0bSH1jn$1`H<fBc?rr?\rd5jAg9Ph3yR(#CTSYeK9<~o|AI?NIxUjo;`J:aZ[X-C7".2=/)z;LO1"8;4SIC&h#3a[-4<NEAC[M'm#*n+CWy}b,FTagwy>sd9/E8E6Fvx:`#}~_g^1ZS/r# @Rhcz0Z5Q~D( ]Yxkf('"y[W pNQi{^j0Y#X'8rg^qeTI0F VH4?|v2+D`gN;4l$G!Ge8H)9/f)vEj1P5"A)Hr/k!q`lI~5jPok`,0L7\BaJ9uI&J'}fi[Gq(f=i5H#Hg\G%_rEK_Ssoe\>gPL7w?iO7O{j$/,jh;!L[^oS\NY:H@!'''[:-1])
