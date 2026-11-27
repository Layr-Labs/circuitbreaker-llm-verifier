"""Active -293/800 first-pass, gain join and starting-majority proof."""
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
 terms=[(0,floor_scaled(-threshold**2*sum(sums)/2))]
 for i,u in enumerate(units):
  terms.append((u,floor_scaled(-threshold*v[i]*sums[i])))
  terms.append((2*u,floor_scaled((1-v[i]**2*K[i][i])/2)))
  for j in range(i):
   terms.append((u+units[j],floor_scaled(-v[i]*v[j]*K[i][j])))
 weight,odd=half_moments()
 mu_lo=sqrt_floor(2/pi[1])
 mu_hi=sqrt_floor(2/pi[0])+F(1,scale)
 powers_lo=[mu_lo**j for j in range(size+1)]
 powers_hi=[mu_hi**j for j in range(size+1)]
 powers={0:1}
 moments=[(F(1),F(1))]
 for k in range(1,degree+1):
  expanded={}
  for monomial,coefficient in powers.items():
   for code,value in terms:
    key=monomial+code
    expanded[key]=expanded.get(key,0)+coefficient*value
  powers=expanded
  sums_by_odd_degree=[0]*(size+1)
  for code,coefficient in powers.items():
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
from base64 import b85decode
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
 W=len(x);S=W*len(y);knots=read_words(decompress(b85decode(knots_path)),2*S)
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
 data=decompress(b85decode(fixes_path))
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
  return decompress(b85decode(data),format=FORMAT_RAW,filters=filters)
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
FIRST_CORRECTIONS=';vQC@&I14dBnLL>>{=QM=s_Zdg5S+(jiRJJB`5k+!>$L7EDdt+JN3zG&+6xZ*(fkbCqaIDwo#kknvBbItr;sgfiR8xxkm@~AaygqQ%cv#6{}}kk{HtChSDqI3vJ>LJtucPnKEbEpRjhMN!?*W&^+i8EEST4wjhHwhtvjD4Cya1%#h&j+nAr)(m7h_Gvhg}GI{An<)Lt25UFEvV5{r`lYG`v#yp}7o+H(ZxhGyj6f7K1ge{(F#Mf<>@w<$Xsk%YN+P<0%|CDv=E!bHTds$E3@a##e=EYI|+9`;Pxv~M`Ui+U;MLx`+_BXzhi*x@EIm*^$XIuNv2f<r6S)}o6v5YB9A)%{*#ZJs_)jIL%T5QkTIsX8_Pw^s7MiP=VK#9F5X;H`Pt{S9xf^j2-Y=!eVD;#J&%9~nCk(aIbEc`TB%fhIBeankQG|0Hyn^_3YD9=Lsx0F0@rb%wNZ<Y7Ej6j$kzMnax5*=$BNQzqJpY<39YAOt<voi~+%SPx<ySZ$q{QXAbl^PY5kDH%^V{4Yth}RoF{^|q>Mg@^zGA>)dddXEI48Pp-u&ktjPcoP)P11NTNuzhJMU&SxAnDKTVy&D<g(b)mEWxU3|0wqh!LJHc*@QpSCdk%EyNn}FlL8}ImJL^@j!oznU~)>zZ-K;pA30=>&Ld6UmC(SOBICZhm^PSd<Xb#Vxm$7=_WJ+3xs4Bqz<Bs|Cx7bt5N+kC`-;k}9eNgEFg458o^Iwd9xn0?64aUn$W;7>1C)WP4w=cCVq*%owt_{It4ceQ2)Lb)RcJ}-r4BwjxiuA+#i)+9ojvk0r<VD9-L{bx_%77pditlV@#T-Xx>V1-*ZOIPKpwcHbML+Jzq~sv{Kcv8;XvrnNK`v;sz0`5Ivh*1&(mFK_f6@EGs>ej-+H-7;@j?9QFtZH3(r3QdB#TfEXm1Z`MM`%nlmT<(Yn5qs!Bn5?rESY!*PBObtvaTTvpb9TJW8|QTu>c7z{D!XE!;s-{!d@7ZYHbx5J>2AFl^d7E@A9kSCFCn@ij%+3KuPr>!RTUwk#Q$R{+0(wv`wEz1-DJ5kpCE|8TJHEk)^iI!<gXpiX!1Y}}wXR!4<Ogj#t#9tMA8Qw;h!%hnl@CPb&N4pIy>FlCHfBP8#h7YkRK+9@rF_8jZ4JvA`f+hT2I2|c|-vAyfMVoBmv<;pL$zU|vx;%??a1JdBI-PA@kh}Ynm0ZovygAw-lI-R#PL1X5EPwZr;5t8UtxVQ%-!nvLC3i)xc%ky&-0_FR*%{Ht#E~>M7D@dCJY7{0%O>vl#(<ky%ClBV-Zz(|`J?#M336D8BpVS;%?bl-3vf|>Zf)29dVV9vO!vFha%(>U)&t4{8w3I_i{8cadcTmxr8u)uR#Y+X@5P_&E?nttown)~0_w!s6~dZk9bF8q`k)S=Yge&v0Z^zoQcTY13f48Yxp>LakzN74S<(j00wasoKsF#Rd?(x?dq2r+F+;?Ew%>|Q<%*vm^GrtxSUQMzEkD(9bO_c3JKw!mLDZV`ZV!GR`+yRd@4QW5$3C2rJt9hsZq?M5Z6#}cSZpdTf1h3Lj8Ir@;tU_%7~f%wEt3$Y)XX_h88MrAHsSPCSBrFcPB|n`kY7s0fdB1SgWqA)36{pAc@f65KKV=v4o)#q&$#FTJD%Oq*!<&7ok!TE#l-!MK%F~q;5`X*59B(V!atw;yA)u?gf6D+fC%SAww}exbFWB@n3M;oL5fBjl4v;x+N9{=623Y&<^lDfG^Ul<QLWSCtMx7>->;w&;eXVGmbHPBJmu#~R!U&TBQ=>d2-HaquciU)L#Pm~SaLxS*zdNw(dz5MBwdh$x9CJ;bW&1CoqR?hFu7la(B9fr%p};qeMy-P`~Yluz@o_#&?#!oK&BW#SW`jH0i6e-Jy8g2I}qfQw1}-RW+2S~p`I(+bz%(Wpm?*>q>?K(r{q!d-9|yl;#BL8GOP#i0>MzLilLL^A)xJaXvebNG?IS_OR^SW>QOG-$G>q#qBMt}W0+*IdzAb*x!pA?R9b=s=^Z|v(OFu%+8=sO89fE)ZUwU(M3Yj555sYLr@9>E4aW6LhhxB0nJJS4z8b3*Y0UW)M9Fp)fWQOl+9V`I&k0(@c;panW(u{jUD^hDqhR1d7KG|}Q74;9o8N6&in)G0he~pagR%DyeR|pNTapJn(j9I3x#=d5Grjgt40mAarVlX<`;1Ct2^qz18Ak0ss##S=0oseo%;y|@>z;vUJjwVHkPr5-wl8*#3{(tR=gjL%vrpx5cy}^MZvH`Ouozr5W&c2?8$d0WNy3<({^m;~1ffTP!~}BU@b+W>?jobr3J&f~SFTdNcG4LG0W~J$mVfG~&Xi111Nk04N?di|)4qjWpVD*2*sS`6vqd}s^Z*bRuzAbGj(zn2>QH<~9f@{}G}tTX|4<K~hh-iO!{e@TxHO#wwK!j+N3d4aHo{9pDrKfm+*F8PX34D%pvwtfuL()XDgmOz&ZQIRB%1fgM~fo6x`aj^TSlPr-6=E-iNd|}PP4o2p4}5uvSSBX#O=WzYu7d0+a$W$qsfht!gF98j%BrqvW&UM?^}kuvpo@A4j=z5bKk1^oL1D^1Q!3-D%!mkE!%t~>m)Q2P2HYLG^Ruuhe9BHU)Cp5ZUu{f{|CzC<)VF*OlQUUmpp;esb$c_NYYf9nnrxUJzY_1lrTSo3Zttfzz@Y=kuAvq(wz4UDw|31x|UHVBO0!Ala5^&l^AQ~=*gN72)I}QSx6j!e@lh4_h#%K9*<95CQUL`D)4dl>&tq6NRf16u@}?m5OtnyFrK&W9kXZe`u%%uxb7Y8WjhZy{(gTFjr>y@c>*S$SeU^yH+c$iPP(s;7&ot;dBE5*s0R-v1=V|{778_bY<-dnN<Nx*`2Trv^dV)P0|q@Ii-*wDvPvl8m*`g`NA_xOgMZh`60}BWn}zXAI(BKjjf2#0i#Fc@7NUS9FZ$|;We;f(Ocd@T*!cC}2boNwAD%%58y$=4kKC&vVYn+>s;Y%EdDSv^#rA@+3uX{84kY=MY5X6|u=F}ufm7|`E;~zs@7n-5s810tt6)%^o#c{*`TqKcgwd}}1+$h>c@)ks@`<aIZf0YTl_tYPwUEk$t*LDbiM((tf1v@W=e}-e$wp0V5>olm+3g@g4zjKjwi`*QL~U~46lR+5D=@QHDx)+Q_u@JqMx&u^<vdc#m|ze8df7wC*RN2S4^puoiQJAKc`t5354Pm=pSmesp3moybA0dxKflfAHxu9T%AGVfn>KYu`-+qQ*K*-1m)wZ$NR*%{6>=@f{G-nO42Kui)%>Q-+i2VaCmkP*pfOpdkB{8LgGmYS1Hr&GegLPy@fGPNc?!t1I;cy}E&v0vkzI%(xlLgn5L?EgvfEQ<Ppk<e%|*u-+(xUE&=**kJP55y>P62IZP{4WkKy^QeK;fX(vztMbKH?x%M1=WlkuMNQZ2l~CX^O_cdSyTInV;Lp$zv_<cqyzH)6uhW_o?94UmJVued<_ePNe{y#x{Kc~QnDrMOIN59>ZP(i>|so*$NJx#h;S^xayK>3gOIz#h2XFLQeK97ww#r-B3aI3O$i+;L7DRqc0=T2anT@~AN2=2;8L_{T1U&2<$R%kV!b##^v)jZH^}enWTyM0<F8Lbo*5si;5f6J#6bwzY+iqYLwBXtRB18n8;m5SQMj<uY3VJ*4TvZ&7Q@8V13+i1$)>eJcr&rl{Z4l&gqZ{?5hKp`K*|P;E=F#NZEEQ?~~*{N|P0XE?%(A+FQ+4YeO)kUFI?`}U}h4E3Zs9c@_E_k;!N<_!*IqYhgehKva_|Co9b2=Qu+)qkXs0xuxBtf1Kl!JOz{s}OpSc-Mpigc2fJ9VLa(SfBk50aMo^o~=P;hMOL~DO{7Mxd&ZWp!;iuDK5u<7AYUmQqRr?#ciA1D5MELkD5c8v2X(LuUQU!@=h{i2?V;C=Nu>;y%&KLCg{Fcs0Xy!$Z$=TS5#&Hql?m&@Km$Hf4NnfxvYmgmAy5W&$@yGeFt4b#2Ezc05aOeOFY2jmp9-H;t3|YPGckmyYDH20jp$kqH4JW=k9XDq)Fcd50{#w9>;-sLvQ9GsZyHbbmm0CbrDf==<vgFmz&B-)J>;Dj!8E7BS-6<^fS)pMoG9KXP8FUOIrw<OTBywEH%UoR3Z2*Au<E=R+4(V+x8V0b*O9LJ<m0sYESLe`!ub=Yt8|nzEo>GUwU(|2UB`kmD@x=BhQm<K%JNXJ?vU<8GS=cX&Xg^(q49|x1=iZnPHX8-EX`VlBcM!<oy+xN<RXYT}@=N?;cOHh-=K3LeC~WF33?+1-*w;{7f5cO@UB)R<$xi&`4NIb+jSB5t5~cVlu4Ra(GVW*M%9vil^V7mkT;nnfTKl`>U&s<1jBjQ>4q5KT?+8V>cll)G8t7bMe!1->h+4|HV*AE1uo<I_jH-8|m&|VbJm#qQ+d<9JGaF+%FW2XCe~6)Al;<Y>{~+WCxWx(y#W+Ly`l-S~gorrmxKnZASTZioPQ*e9WEt&y8Hu=6|_)Rpo&=6DZEq&02mq=0X{<_9Jq%@yd0@5)m51Phg=fcJTFBbGrVM`M>pC{@n*BtQ?`07hcLmJ0X`sDd);Ia{yiTBks-$GLLP?A2di2v#m{j#YEtBt<kwjNN#UKhyT(2VucL5-Vl8dlI&+5VCmy#6Id`~XmwEu#MV_d>0CHDw{I%pJ|j~P#RKA?FnKL31d5;U00co`6>+?7+%TcRhUA8<jwXrW5Rjnate<hqiBE+jPtq0`OO1Ug@w#wmhg1&|lemW&3?Kt-r$y?F)A(DE?x!^pR6;mEzE^BE?fX7dLepFdbtGM4yHSA-<-!Nw`$VcP#e{51$beH`q8vDyJlgLu7f-^}-bdLsJv{q98o@|k@{y8_*gzAJ>C8j@^sEGgL__K}+vZ&unbH_1?bs7jIbp;PF1PzMai(r@FH>dUXVe*eOxKryjJO4<g;8Q{gMerfsGfNK$1=WuD1{{VH6Oc-utUnd>^FZ4_zB6`qMOnvXU{576Y|O;m*`9J#4!n0E2a%hD)o47xz*T0rBybt*5u}J9xu<}Q{|DdQh;Sg0;=hu_Bik72!OB@r(Q4B4$!z*GM*5npyRrL4yj|S<`(pwo77`hZ-TWPqc~@t&&6%_1IaVbDs*4oOGk^Y=z)sL^Mu5?NkR{{=JHSQQb^cAojiG<QmVw0Aqqfukh?K#Q~;@=lM9|wrdwJz6~WFfW6_n^SXN^5l#7<hoXFpw_#Q&A+j4kq@dPYpw=7J9f{s~W8s~jI0=;)p2TNV*9n<9}0(hfx0!DXCk%gjcV9>N+oyeW)kD53ez|?QV5lrVe*Q}>iQa){-&1amG&wt*=>M$K>76jl)8N4pnZ{;md^fhW9oP0p5hMCEZR6)393Cy?$&1fCjDJJ<;74AC{o`3juv(6mCidEOMdM1cL=G{#_`p_V<!$Ra{hhO|%%~EtYt|)FuHzBq+avmn2GoW7Dwc@W1KCpdQtBLf4#esjpQYeF?@9_qi5F>O(uwbU1y#IN`613K&UyL@|%mRb?MyY<N=uNvEPlmq6jp{RV8m<3vH1=20BkJspN(GXkv*;mRabo*f3Pbsx(CmRY-Zuh%F6Z;>yZ$0K6iCx^yZ2qK_GXCqF<SS+?f|Kw<!|_ex0RSOG-a7NcC4%ZPV3Lv?-F-EbJ%CJ7vkS;X=MtW@C<1ZYc06n*DUMN=E5=xNgSQI?yES5|5%~gV+`u~$~(0olV0k4laK$z0cryUXgK#wPXP$9jFjC}<c0aL*TRA%l#A<;!djG^i%q0*78BLY*gbehY74m~6@cNOk^*@h+_UdfGQHm=#e%0_RiRr+^9=eR+pHW29e%K!teeil*|0S<^yYI{8TD=WsSF*<+Aa_K*lUUfe!6LYqX2n3^BocVAlg2-(pFoNj0#L;KtHjTE34u5tN4}fSs0$2uWufkWSt5>7a@eyi8lI5LZEMkD@R|RBZxm4c*|$eonuM@`Bb!T9?R^25@%6Mo8~Veka&N$ND*w#w}Dw{?ic*_#rNy8MGBzs#F8oj<lKCrkB(&^SE)WlVVdx10dnUuZeV8q{_)9@ldyFayI{^16~(Y$em&9aOHjx5GwrCCDP|s7T>%d^)SgYPPyKk>uHZ1gc#<a_hG;9NI4~8gQ&jV(<VU2w6D@+VVqZI~C;r=}ZyV7K!z{L4GPykjDNle)uasT{78x9`;-*paA%Tu#zp8&2#S=?Q&dPGE<q_EMw5)tUJ+BLyE-6P5tLgP-O*?dA0m8QFC&+p%r1GdC9x;7i)__!0zN0SUg2(1#R}Zrdr7A@6Z@~DH;0Rv-tvTi<02~^fD+XNxT1okKj#_s7*@HvVO)q@*S`x>n_H!L~p$EKm%?a{k78KL-mC{&jqFsVlVwh&Kc`nuxOYj`l#9XHV<9uwh$WL(a!QJS?xV}-i_3Lg<D?8WTn}nN3(GKIQqfd0W?x3&QhBuIXl>V!ZtH!5x>z@qU8(Dr@vyklid~s4s<BQMuT!7OV3oxt@sLB@aaN$xrnPICQY-d_M|AlB?gDM(%05<9=7GTAbN=X%;nneIN&Crk(_uN$iN?76X3{313gvP;^>%tu;*^0+3o%*^H?d3AA1v&Q#fh<3-845^3WT7*X9OG!+G6>ax*4iA0W#sLdrsT<?co0d5Y?mGF=K~ocfPx2|d=Ow!m;hecV2YH=%1_dX66&#LJ#`xyM%<263SoWsY=6(9jXII9_+HnvxVHV&ZdBEI&HpaaS-U^S>L<E}_)=(ZgfsZ}x!}=_QG%UpwK?KpG@a6T<ZasPr6K*2`hYvUOkv9hlR=WFdeSx>*S4OX0UXpPSaALV9ugDQ3|cwSF>QyXj2Xj&<GWuqBaCG&#}D9nZ@q-DmtkmMF_IDMg+V2gRh|*O%7q>t)rbAV(xwd=uvrBhC6}b1z^`JKUoF*wGLyvaqtCmy6vEg#+#U6&-$q0?Rs%B{c)y~q4yI_wvHG;=O|bMXicEgvxDu2ihfn<(NX&`$2;~^qk&@$8!aGpbPWa(=#T~|^1d(^e9QRHG-j0r22q;X=@34el#?(tI0<hNwwkRY4vEMsn7;?ASKe;O~wQoP55w0Sj>s&T+$}mJ(H{1g{EXu?W_ssbw@3<t+_TDQ+`tl~+LA`O`@CIw4^-1Q>Mr5nQ#R2=!s{wezdvbWd>>AH2LNL@GPAs;7=ADHY03*uI@hY@2hPYNr)eO>F^%@b#J<6CAxUklfDBJZ`&p(3xa&PD3h5tI{Uk7!c_{HeBoa@Kz!k>FjTl)5UfC2F5q~&m6Dr|{$H3Snhs2<Y}l=vLya-9uA<}cPpBLZV_)eEq8DavlW5|Csqd%`D+T8!=GltLYuFEG$Y{o{N$EU1{bIQ%-^Jy@J*&`y*^C4yW)hhP&1sYlQ)16k)ajiXtvbyjqng7txp(ILRj0p;R-tpcu{1~6TmJ!Mh1-hWPO{3o*k6Gyo3Oel>%W5L>(!G*1tDMcn5!n`>bDZVndh4T)nr%ARa^J^F(Ih$#oAcJ3;BCK#ijD$Q|-MvQDpY^c?bCJ8ZReB_u04R>z6btw7vY3&9bp2@82+VLKLkpoQ8#K5Mj5(sQ$iBh%@7Ex)Jk3UYBPy$qB6U!9g%H9iZAEtn-(yhTP@wyw`EDj-y<|WyuEuH4I;U)J&8tI<if+?pgoUzb*s=~0VDOyYyRY?x%F!|-F@_~^@A#W?41GG%P>tnS34Ml)$`bOw5xW1zZUc%yf?y4UeApF(Kog8)dOm!Jr+zntQS}zdoISr;3bdwyg!jCZPKm1ZU~lIggCsPGI`^dlRS%pL^JpL~9JB6jkpZyBX6KA}$S}o?x&B#)(HJ3gdC|J?ldC%NdBOHJNwWiJk}uT+PQ>UM|L9poi9eimht{@i${c!Wy_`@QUkxbP#AXEs-mtGaYuFo=OrLnR{I%$Webw4#r9DaB>o$v&cpPKrnNXI+Z3}vs#psc2M7N_3l-83x=8hL-24mk?or7y}TN2%#_672ph3C|~0J$8EPSTuu(_iF={DPpQBx)*P#VC0{cut##y*u&8a74#2YKS>$-#la`IK3qi_voH@#8dIwq#>!U%Z-+kdL;5`4O;oPUTk61rU5Rv_XU?TcH6=5<@AyFT)u~GWBMI}HYy&KM$)~4&~SP+^vtdMgmYA}7s4KhI3jN>nH&W+PGa0@4q()usUA90?YO{t$=e}UP!C-yKK;E7f;t|4fb;TlzoFD&elmyccXKUOkzONLdR=G}LFv_P;rJAJyo@;KDg$v!tVHysqi$ORZxNP$Q4B1FC_cro!%NEk;-2`Tc>9ncchg=$^Y&Kum)$>J<72#EDu)Dw+-wpi>Fb3YjVm+ASx+;T`+MI33-l&5f;#||=at?v8nuc%@TEwC+N@poO9Rp#ig;Nu*?cnfZ$YxXSiDVH@x_<ddW$zufse;0U;P^$r|1vN@Sxa71g71k;B7x(vo&GNIYkif$P0F6FVe1M&z5P;q2fkl<~l6QUr3?8bmS$!1C+&He+@oTyo;5)e5K?2WNyR=K>_9ap<^keHZ`pdj@J&Nu|jOs-sBM-6fh)A3|rGD^38P4N|JvlATM)uNjAuZsJI%JI<|Uol2A+Bk2CED(^MI%Vie5McRgpv+|#2JJqEn#Dqgc};%^Vx43Tl?Z(CoHUt$k4$k*=TPYNG}RYwZFpZy*9=m9ukCJ8~`8G&S`FvKi@Q>LjiYBM(9UK^MlBqBPN7qVy3nu}JWcKpiwa!OOZ*g0KPss&6$r6uVq5&;nY*H@H~Y6E8Hp~F;Or`emQMl$rV(GFSYIRpnW8bgjSP0Yz879*V_fu=*~LXxP4|NhP<K6exLWz#;w{?3pfd2<h*F%noCdh5^k5Oxg<mLzd<<e7`pMyGb2b2~GolN2G`vvo~vy28!quGB(=Bi0`fX=M{XUiN&(J`MoZVbyd$R=~JekS0U{@2jiG4(3I@z=ugQ_jL>g)CV2r>h^nkx{{D_k4f&(foT^gZ!iE`2WXAi9I8uJZ7>$Y*WbN5;Zs|Y>zM$i2e%K;-)9!MgDkV(fuh*+TQdJkZrQ?(&j3{4@U~q$-Ps4>jJ#y$AOvpAuSGv+cqg*GU%oF(|J&FG^P>QhahDX9*?qwL<f|+kLs5kIBKEtaQA~nRda{K8;sLj-E5$U=Yp-5kJq^3H7&?C{82wz@8?7gGv^RjNHpkifK9)igK}EtksBE#*(uG(~-&gL|5U+7_kwaw{M!yg4tX+IZCQT?c1Ax&sR_dc=OxIT5WbU&$?*MWVzDo`8;+5+3g(Bq~rIZFJ31H+0vxudF%FddS2Sj)$z7|;6RUCfNM$qhK%vo_RKlNp3nq7rY%#c5qK70~&r3)lMV@WfK+h!dwW)8|%K2)xAO>$ZL<DbZo?bFzjLM6S%nv&E}F4~M)g7*H@7X0A_YnrV^hJKBD0|Q^ixv5*@r9pZ1y@sFF#lz0Qkv(vxeRtpC`1^{YI&Erx-(tvp48?KQ>SWY*jTabVL=+s&M$TCJ+px@{+5*7g7=pNaNi%iYR6(sNwCJ@rN?7QqJ$msZ`}*MbmNMdd98Y_OP4H|BN$ONnyc=$gWOdZPj|2l_iHltwfC>qm8i{xUAr7joMzu6Qj$aaiP4be}&!w!*3TM#u5qN1woHN38G{F@q2iT&Xd;f(`Gz)79b8+QN;;0`|?z^Or4*8K*)#~HLLpGjNokTmpGAiiTsD_`qRgbvGj-H&DW6^^sqedmw(Y!8uckG7>s9X3wdl=bKOO})@M#PhdEi3=9te_yg3{4Tgk|CwBh|%^p@0gO7SwPs$7=k03PNR|0-h!@}t<KKLL%15cnTC4X&k`BB510n<9M`X`htNjiqA~6~nb7rc&BtedgeRsiQuSwyO`CbGLq5<qf(}VuhTgr>y3u2`e$4l@cJq-CzWWGU(c06+X8kQ(*?7ij_3EtaVM9Vt|322*hZidd{+ft{J8E*bK0E&^*e=A6JCs@y0{9ivYGym#Zp1qY|G?1N2I>h|fDVPC#_O<Mvj^svd>Pjc?3Eg9-uVtAl^}UjW=yQwXu9y?EHWbMA^&E}xwXEk$>LE`&uew=FC{adf;xEp%`g5ypAYc{GhF8JKqkS&+_6h`nk{KP&;7`N1z7L2pDP<_2Me;<yd(Cll(g*nMC>=LmP@zfAq#Te^Uu<E`Po2kxxu5Kb$<te?4Igk_|_J2Gxp}H88&z@9InYRgQa@3-#T-mKNA<QoLiIuNZm<ZA-F2(*A{`ph_OmuhNY6y8+eQa%*iisnH4wXvQw-*cSCIr$8=;n`4gC?1;3%BGp!<v0zXdz{dA}Wy4e!4l?qT1^RYDM*B5eD`MHg78ZM{hg*}zWAwQea-=}v#p%S*dB2vC5_|<yVDWq!L@DR%EfoVh6mZ)^3*2ZWQz2W;aGpB3JSn1Wokkp=fJjok%N%>Xm_d2R)@FUsVC8t9~1DN|PVRUDM<@tf0LwZZ`fn?@E&`jg~eq9Y$cCAZ2nO~WUWjF`r8Yy3VCbj+;dx|)7j1PeEbUp2pOq<Q~?soSN*_Ha2TSNj3!--Kde+V!kpO?rncxJlT5#O$%Z3jTrk}(&@H@h6@aoz=RYN(Um$(fxHbfD2$MpX9Xs*6)`Z|W57F|kpSCgJQkXzVC`x39y!N7zyYcir;TJckStZ#c47?P98TZoOVy3!XBzw#vI;SG|_Lk8eIpd#2DlV^%|#&1U==2GXE?05<^P^%I{rRyvSx`W(6qJzgHc{;o+6$x6-M1V&ef5Pv+!SaNeP)&AX9`;dNVUhz91Di*r`VFY^^L(_9pkdTXm*+ZY$1_Qo2?PYaHHNBrYySGy_epohXWM30RsQg+tRm;cL?mAx{qsnkC^mbE0l*dk4<j#ZZs>(jCFsk$ezMVr*VP`fEn65`{d4Z=Hz$#WG?6%&3-q&0Y6gl3Y<V2*3j0+LDkY+fL>{3X`ik5+WlqK75uRW6dG*a<bLD72C3{c?7#_6oRFCN`OlGok>C6>ZH^TC}byK6!N3N%(pw^qy&n5`6g?WOSC*5p2KPp%Bw{*PiGjTF0a8|qEh5)RgmY?vUF`~iW00-R0G^hfmo;GU28Mn(idRsPxfS<r)<88+e%MVlqf6?NyQG^B$qeE<zHO-t1zneYB!#|3IUa)QuLlXXA}y$I}-A5?XG!N(*hJGA|ytOGityfvga`3+RYifuKT>Tt8=*si9Uyz@^vzYdr6&E**ybL~Mg8{LmO7QeP*#_*u*CcI*?N1NpNECRvXdY8Ugq$*P10<*FCP_Dadl&e*<rHz0Z$k(tpq7T?En{4Qd#bW^E8*kJ$7<W_8JU%?k=t-?8*voElyuRV(i2LQ4{5;6M3%{(2lqC|#i&!9T;iS`{yY^48TL_hQG<M=bKXsSMFu=ioz7N~*Jkq7-;WFWnq^UWMJZ;Z=AEpHx(z-MOvcHj>HvC*FF~=#TY)YkQI{}dGra6_3&&=KNGcpX1U(etopA&aooCBk!->_`<unfZ%&SBziE~4l+nJ@LYd52FNG?p-6X80L7Ht^ZK6YQloCAO|1&H(|bXL!H*>&l$>V`oRc(4<`oYX1#N>Bx81GS~YeXd7#h@jQw9Uioux$x#?%0Q#IDhA|4O$ScoeM~`Qb6JvnI+TZS*Kq_8@7D*0Y-l)2^0CpTWp0Dw)c&*~JilaEBfY^*iRZ&NuwJMTf1aOrF379_t9rV%$phgsv0}6^d^iuOm7*x;`u=B)i3?iPiDr7fDeHq^SyKd!n4JG#hdLA>V<5!M4TMW{2X4VHz)(fXnEkk!es=vldjYKlv#b-h5zj*7lalMv*nIZ0h&~a!WT7WBky4Ol-d3vx_1Q^;Ic30Mai-DY42-IfdavUlieBv$4nV0<;caS(3*D0yoo;!vnxy*xQuyqIM^${^ZC|!UAv)P2H5;n{F)B;r554bK6i&jykdH^(&_>-pB(XgoHj>BRldG-mnb2ElCGTFWlN!*M5-8#O~{+_yk%4E;&qu+`}L@^KQ1hZt!w)R!_aQc$9740<lmgK^IYF9l`QX$S!TlFNx=D=#&8g_qfSlS+-d=58V^MT4nn!so`g84a~Fw%dbZEmGQ*FDdi&jkj~x^GnvSk%Wvk#A89Pm}pQz(F)e?t(u&S>QTJY#Q#aiPb#FV;8y3v@+Hso*$RpW<aC(IO7gh@Q**@B!uR0OabNzyGGP*eoAEkO+UYXOQCc-0#FT%a2Gy+1z?U~ecy7S@pK|l9&B@#x#a3y?5PDBQK1Cm;<!dM0YItJ^kHw8+=Zr7+Gc=?3z?I=pLp0g%>h}!4bOwTIM~$e=gxvJGQODLa$<Nw)rTdxdpE;&B$WP90{1>h^3O}yI0O1TA=R6U`$3FY_*+d8wGX^3M*-s32thAgwT2?p7Ej?bOAQCCY`)g-r2@)_jny>AAqz5g3JY_T(_g3MX_NUB^bVL(z45w|HR>PGRABjmn1ut0&*dPvNR4UV5ieMg6Dz=zowpvbWHvP!Dw6Z?2&=!UGK6%Z2NeoB%9oxjftII1<4BeEA0s`MF7|QaSJQz1v5D69J^;+^Dq^M`>d(rS@_7>ZUVe$c8x;>euw-YRUIhLk0#DG?iIs9R+N1>e&4VFW5A=+SxZJf`tu)*g(yG_t>NUEM%y@Q6wE%WDNnQid4)&n<p@DBgnID<tSbX)&$Bc&KVPRtqC;60=C}cNM{2>5{cPu0s@U}QBHh+zTODlGpK_IAc!wW@-!>!Sf7_ohl%eMmdzsK>RO<1m4Dgp<YvIk_%Ty>EsrakJe0A{I-pE|wh;Iu8df1DOH(c}}q+X#%kGrk@{*DhI%{~7Tzs=<b}&fP#mmlNu`CE-sLwk4Fag|MzZwA9rO_S|_K#$iGwPPkE7GsO+Er$W565ZssyiJ?YP9#P?Xlv5{QYmU;6@f3}TtxL&=AZ4>zl|otnQ90RBVkQq>UxCzR4t5;&`NKJbrCltZ=m_OZV=KvME^yqXb6%#~rSOZp8P0TE4@U)K3q;u)n5<<f|MA<H+Q)QMCSe9S6w=Jsr-<eVHkrm^$tLFA@EKlrBOUj1BSza^P{UuA4@U~rln@-q;S{Yb0%PSMHugIPRYY{9uN(Imu7Fw94s^i1x{E=mcBjD7nI28)ploPQ<l(i%3UHep*~=NfkppFmHTkPI{6!1j#%mAhy)Q4xl7`)mkY9N;vgD+MJ7yfvBGunR8k;280HS)cTM^SGWyTPmx=L#>mOXdUbHxF}+H9?(o(*kEjxu)|N^of{(eu-ZlR2jO`d}W44<>R1zZc@1cbQGqCED9Pe%I%`Gvc;Q$`2DtXdh*h`lWNnLE-MG4viN0$6)$ALy+3%iq4?Nj@o`2m5zK(jWhHxE2Y)froumlK8{cZfBb3aNDXpyGVDc6^QPBs8{7$IT}jf7ZE!W{{kL8#NU684->eedB5D8A={#6;LW-<QMeWe3PVPlwhI5vmbEyV>jh|MW)#i}*YaIj_zuUrXJdQ$F=tl6TN8lcLST^ISnB<H!8M^8Ot?Pvqbl_!Yvmt^!{-b{0MVnJsIBtxXbGr%~QQZt4R9oZe)yxZYmQjjr5>D{wbZ9VNdVrWQ4tRbZoNQ{<`ZBFLL(v*pss+UxPgT}{QqAj77MD6t1Sdy^su!@+qnnV=wNpZ__G9*<aZ|s7k+ZmF&+wkiOcq-^9G`;x7`6aMhXYaaB&q+!t`@;!`@cw+w>ia}I#xmK!Q@~cYk1Ad*_E_7$=~Vbi7CtVcp=`6$nJu+1Fbey#-{RDGbrHMv$_q{Y@d8+M^NUAB|F-h(g4l~YtoWjWdMi_x$z%`b1!m<-)7r7lP`z_*b!Ia88S{0&88WhqUPSm#%&GSy-L>VDC<1QfZuY7^JB|*f#lHDxqlFG#muU3*rR7xz|mmnO}KeW8m4i1MId$g@O#j@Ka51i@e!;Rh#i3PewtGcSF{b&(Mfu_-xcSGefN7)<hzT`lH!$QwUY^N(Gp)BC&Iyk1okEr$lgwK<a;&{DC*X*=2L~;O0j-mLuv^p$Kx3~Q8@#u57x-E!L$u{^Hp3T@=}x;_nF;ccl*ut$Kl6$OcJTvy0<JhppOJ;9U42S@4SCt=cd*8NX&-Exz1)lQ%gdHxLI{YnQ=4rzj|UQ3xL|`7n0WU5b~sU0kg#QxV;d*zp$90$!69JA31K@nQZRQB{i2?_2IPCNTX0aMC_BVAwgc}-B8^r1H>R$T;CM{Iqn*}zZ?29j0jpewF_7voTrIv;aKO(1+SGeS%#4M!>OQLR#X-)$+S%4`l(>dGiuMJy(a~Li`0e@06-!Gz#^3L#Y{c2qCtzaj|r5%O@5)y9&_sbowVCkwMPf+eSQXzezzppSjIlGzfqMr45c`BFG)oOR5OKWMW8%kv7Mi&w{>3!LI{bCwb+MHU@I{m<KIW6x#7lOm|WlvM`y78yeIjLC{)&uci|zfgy||xLlB62-UxAI3_CwO%?gI`Js(+BZ#r{8WMi#@0W8f<r#XevUHxuPhnz0{p;t&LzAx>w<|B0jcQ`X<UpRSIXJ;!y%Js6Cb&NfXCizJoMaO|(3jM=REo%IllUK&nKdX3Sh5HX|F$k%_50L(^0^TT+7E=_!zyX_wUddW|Va7GKS#RLOF|;Fr0}MfvfG<Tm7%~*;_U{*~>&)WN9ihoCugA6p{Wymc1UZJsuBvl8&w#w$?;|=Z!E3@`a)?JTjOy3k$I3eL?&F-Sep(u+!yxPEU2jGup{uLi)K@3Rj7}=i!+lR~8>KhGGM=TYdbx*~<q9>|e+|eZUQx4NqBflO3+f^2o*7G_4pP%s`Mip%BUB$Rk=0KkT#uVrZ?}?rO<B|^nG$au&=s_tyEu7Wc(oxLNavIX6nn;_ee}J_CjU4Fudb-HRgIKOVdT06nLsaBJ)%&j1PuiRBR(=v*!^zt&}!S+s3}K%P=66ge)T5E$pi1xxyAI-orI@fJPVoOVaY1R#hbZDX^k8>|CyPC+pnxZ%l~Huxts^+a>I=(_JgxD(b2jz9pflmIm(3ffe0#wuM;=r7k)wwQEAvI0SWt{cpVIFTvrjI6f`i*JVu-Y=G%XS4F!g?XvtEp0ynO*HRn_Pak8WbeymT(CcV#lf}gtJ*K`S@^PMs`w5lTLnC)B$1;_9Qbmai2vOA?t-MZ2XMSMOKCoLDA#~l-ncLFBsooB3ac150YNnU*;N6ZIVQl^9K#1N%p*UD<6<<_!JnGzulr>f9BuMv#v8l@~=K0T%Mn}hL|H>-$n5Ew0n9Y%^xsEAuhS+~*mo~HjgnyfmKv^arQP;ku-Z>7fdo8!KLo!XDmy_}p%lG4bc+`MR<{WSH8w*m)U7@RYd4k?!_3t8OT&tqQw9DRkOV8GKsXnYV6u3-(=Y9-<oA>1S*aELqD{EHuU6q}r%xE<!FpP~N-Yd1-!njT8}C3beK>|N;*lT>QTt)GlLkJ{g~kU3TCA@%tre2Aqi`wMNN`Lb=pg%*z77>RSQ{KN6Y<)_f5?=_1Dz42k`Hxia{vY_|h>+_wl-`YnqmAkf6v)6f@epv#CZm{=8s?n;RU>7;a{YR(|2D^fR6(G~e#7gmS22-si%35La=01$SM~8e_`1qQ?SZf%0N8Uf03d)l81B=H0Upl2t=i!S{bG~EsIPIFW@2Jr#hZ=#lJgB$fP)`y(<9uZ|Ue@#{ZvOg7UGsu~tB0S#;DpAlXg?%O(gtsphmG1bt~6Z-+)?JUxE`4eEIQRVrj$8$JA|d}%aYEFDdr-EaS4y(3%os+<f*3?mhP@8;{?MB6bLd42ns>q0p{x`e@T;ZTgVf{tF6G_{kc1=!Ay!p6?YwPV#{VGVO3Ls69lM29CK`NCuH(1;*Xa~v~fmP-{wBi=M=CAv-zZ9nQ)si`s48c0o-Nn*kg;$QuOyx0(iK|5uvS9r|2Z&H%e6D1jrj1;K1v^Y2QEol~K+@bEK2N<xCpzCBB%{%B~sgRwD-%IMs{EYVG`s??#_#4jWsJ`zyd049?Ea$(UGIb-nJLE>tLI;<sdz#j{{g%4VLl!bKRz@peHN==VQz$l*jk&0;Agd>a!xTV)L0&d6*@(A%Twp1=&`Oa6OC;<76{nx*O$&XAUlfS%V>A@jox7dY_I?L_1<-a@y}f9alqBw?a^)@IcAvfO}5hLf^pCtE!OfQtjc#Kir)s>M2Q?<{F1SE)LEq>qlpkD`>XR8c?#yE2jo(gMZ4<UH)p*}i*TtzW335r`26t=rmx2eH-#g{e@&&y%W3LjlkIRM%IE)6|J*b|U23lMUf`3j(r1e};$q)PF;HPOypNtpMdqJR#!~hWcd3<_6};di=N1=<zEtUkMUVv1g9uiZ|Ii;6C(SVC}c?q8kz*x^co?8*2V9)gmk8eDAcr0D#qUz6y@O4-pnx+|L+zfzh!XEa7|yT-qcWXt?l*e6$X~R7FEH+@mLjRfs5Y<LZ(^pu$&sUU09wkxYkNF$A#GkZ<HU>x6cATSLu(fRt}e2treiNh&`m%eyw4$0(+D(ylQ0uPM+Ft!Ds?^zC?gbHC+HcC0gewHU_S=ntA(=LbWG+2M_xX~24IqMsMabsCsJ0;_LIxJ}T-fq=3BKJMAaD%e(*`;3ZNK<2!2fRL7HBapL*N1r7e9Bc(9MkPPgRPdgSF?)3Y36duF>3O<tDmcgjk-1K#e5G^KUt4IjITw9S*gcKqln8+V5b(Omjg9+6)H<{QV3=c3J1ijaZ20UiP-OJ?Rf$sAfpk-j?|Rp!)4I5?xJKT|$wzBOo>-$fuQSF9r>4<}nqP5)<5ovpUa#386Jc)!g|L4rttyFz%d+*eD9sow(FqLbWt-mf?zVLhsnEczhDjzQuiJUWXdgT{j%!jd^`Y%b3TFK14vhPPffEAB32dI|!Q6|sLFR(-U>8!L2W^{Sn~{vr^|i}dKq{*(37GrDm|iqW5mpptRq+dqa!uxhNQetuQ4Mr-qgD0a=O<-R<}?!Q$jioah@gAv)}C9mzdKz?R~Mac^zI!Bc;=@4>*O=3bZ3_Wn!+P8fzf?+!UkH)H|3mFJDBJz>L-OShYRq5iYP@SKgftOF?r<iwQ*C%<kSjW%w-E7BY5x!=adoWLXh)HsXHZ%!EZipy)c3i!+XLDSgRU~Qe|BdV*PSM{45o@{hsLxuW>vXZ`olWiQGD^(p5qGZ%Xj6q*@f%Jm$KT1^>vOL3!0p@&(-@=htVF8`e^ytdwSivA(L#27SO6Ou>Ig+TofZv*z1RyNCr<+^h%aL;!9GbMZ%gDd#2#sx&>7L21|HbF0AJ;y-2gzW3c-4)ucQ+&1YFQkP_2rpjsHiqA9D%+5othNJZ$l+2Z;%{9f1VFxj^3*IGw$7dvPMzi0W?FWkZT-UpfTB(y(VIMErgv7iO?_7K@`W=ZZONUJzT-5xF{<M7eEQe-|Oivaeuvvv`7wgj_g~_t1aG)ij5JVl!X);U-A0g5xo>KKcI>y@zPp1;Cl;t7_Sm%^7Z1iU$Fn|QDpms*sZ6@qc;V#ZWg&(SRpKTnP0py!`v9NxxAq7J+D>C&nIu4wANSprKuf7~_rW6~Q+a|sGGg$urWg+wVeo{NlXa)hyd7x1;R-9%FEkN?w^Q(?3nJYk|AEQaGdsK&OgY;}`v@i#GrP$H2g-PD&%IsT6Q)hh}KcsTcdP4!LVj4!07x7ve3*2mTgHs7$DrWF|P))7M0p7W?HFaPN+uwf?3o_3GF~#ti!j5VN3~*uZpgHg9&H1D~gK4BJ@10izT?ZE6;Cjg=K?N-Hk(SE{e5;E6@VA$*g$Zx1M0U}W7HvJq>*HT49t+8k#u-?Tq(#QfPB=~6ub%mVc*w{wuvTrEZ3PlT9fD^?TzRJJuraW{TyL2>Fo6f-J7n1_aB#*Tekyn*=U6}%V@CY?(f*(55tf=4GUL^(3B!qUJ~70ENIbTuK%{iFgMBaLo!f<tF;)HCDkyFPC|UBL^#cWMrVE>*`+M)V-9t`?OAp-2W1Z`t-8#wodK1PO40co2N+nhS*Uq42Q!G^3ztJ`*&cGn~!S@!?f`T`vT|C!`Y<Q~}^G01p1=S!CCt|{H3UY;~^|E9tX9+>zolNr0xbTQo+rLy#6vaUcH#?ZLKn25q+`ml*I*V-y+?176d8KUWEI8b(G-cIbvXU5RdeI@`Xe$&nz~!50T)`gWjOfcx<U>sR3||R-9*vlU7r%-KXQoTSmi~+YEH@#5Xz<Whe}#Yhv_SzNTJ*<svnwAksZiW$_%es1!YrlrclDOFuf&9&nj=JiqCeb4GWE~s^cbK@;${uc*oOZl|3ysH9XJvB_JMVQVTwHZQw@ASx26y(OSJx|e0)O{DF=DIKw~RA9S3cy733ND(SOPE38ZXgkb#R0`Ju5%KXZbv905?}IC~^S%G%uTTc3n0wp&OE9fcst>d5~ZPAt(!R8+*Kl@s6w7w7kE_lS8hzaZ*aSDAgd4Gb2b7a_JUIth$*iKQS5!~i#kpWTI3Bv2eA0*manEz1~ZIOojF$-vjZ#C$60r)5<|=2q;9on3jtMx&5L!heL@*!;(yLRxMxZb3<~H2P(`-X}pWLN|AEmQR3%Z=`v@1K8`SnRi$GZ30&!BmJ@ol5TO(QZ@XcRuN_kt^**L{GSr!n&f{TvxSTiT1$!4b#tlIXbhkbBgOxieNUCN{Z}V-;wX5c)KAN=Iqj4S@#zB%sBGbU)&b2z<t{k}8<YPZnF0!ltg?Cl!FkyE?Xq1V3im|myQ<YukpWkb*3&J&D+^>28Mk2(4(LKrJ<nVBK%Q{(D5jo1saSW!qCVY*8yIzR$-rN)I9u56!zV5i6o3zETBj~9)eZ8Ky*h656gx(Oxz*|?!Y?)=f{#c_>Jn)KMH=hBp4j~z!5mMGJX2w-vKaJMDBAd_;^9EmZ-_C23p33MUW^bTBFpiKMX$zqAw^<DlVq-i)L47&Mn3+jmd6MGTvhC*49{;!`?g3AW=Sh8gjNEMU2rq~IZcI3z_uRUtVy~(cdeM0O)(z7C2;;MkxC&)wKKc);s(J1Nh~xW<L;4D?;!LD3>Fj&It;|>?UG_Y>a#b+!mY7s%TB9~<FLuMDA67b?RbA<pJ}C1Tmzg!wC@!v+CLVwa~!+VzruL1<+008Wn+k@Oq<Fd0)0>|%Qe2k4!(BfQOrX)O|?HxTnJ00-wp3$@#llx+>|030PLkgVsP>aVTC4L{%Hoon99^kG<dS0N41cfoTfDNk$0f1=$@htwblT2H*NTmf^J@Qhw4Z>aD{)(fpu`9kyG52aI&j-D;8WpMj(^q_BvUwo|@2S6ws?xfz-eO=(<}m{!N?cA?OSOE9y8SD~bfgZKuC%hqFP?AF?3_0ZkNMu|B(xP}G)@>>8L`@p}al7v*7pcL-<!`W4(jzb3>0l+vN<I}s2bWSal;Gus^Sel{{ht0FWvz!7hwV9#HccFm<#G`dMpI`oWB*Y?%zdpc2tON#7@l$$<u<@O38vp~<VAU>(OzcWlvETcSE;bCR_^MO5A*cx-hsU`{MDrw)<`YA-%zA9*wSr6bi-6o&;?;SubgW-;MvFXYwbjbrTq{B~g90#x=srf@X$jZKUyoUmGp#IFfBzm=Y3Ic4(ypFVJkpg4fB&CR~T-oKpK^vJr=V4O0MiF17_9>j=nBLbR(3L@MZ%rKeZOL?8Z{|6FzBSE6=oM*&+5S-vwXmnm^q&Z9e<YiKX3s3`^+j2p;UPpxo+@*Wj%`2k`6s|CAXM@gW<pV+;J=0pjFt;c`E~XmB612;o$Gj1Tr>Cg$3VoGnXpq<)jihFVW8iEX^eL#zvMD4D$F_odD<s@R&}zw>Up$G;Z|*kxU>Z+-#dPjK=oyihf5gNzxtdTGcJDVVy7m?#W~Y<rx5+b3w8qJDIhM)-chjGuTdLL&-&wmgg7RA8CPnPcB|5vS^ASa?Y*c66N~2!y-mLfWMdRhczXr@Zf4w8@Y%Yko+Q$pXNL9)n|ljJ^;2ybxs!OTvX*0@_sR0%?}T3oo_@H^Gic`YzdXBAL6rDzGio<&xPm&f6FFVgNrV@4U<b%t(#;|aIr?mTw!#oIE+n32A<h4qfKh7XqX&9yX1Fc(Pp)_RDy=EKshD}o(^3E;pe04NWT4+CCy2Htef`7GurtA)i?>Yf$X=0_6MFG{3ugP>*gE~52GLn}>^)p=LcB0kmQQ1ihEDMTmQN!uQ?BUwO$Ree%cz!p5K&|^8+--MzmXVa=#XSGxnBz^aOYH?giS@$2OVE$kf54R%Rot=10Hg#=R=cw4qja;l&m7L!{Dw71S@zY5O(U-|J=O8Ja~coFBSYpU7=|sKc!8iwi#}z7oO29@jV+rylum};llhcLsOK&@={l!NEfflra~LUzd|dMtid|$kv`51K4`Xa`uUH&1YB@fA+=HL(_OkPXK@7;6wO+M2x+#tyJ$vABV)n82NMp56v>zOdF}qYvuD6RJ-GB7;ZdLSnTA8q2kl*Lh~h|lgq(?bq6x%^3FHP{pgx}3m$&!6l$Y4{mj~QRsIe`<hcbUQ$9q9jWFJz%EDlVWyf*QlIo%C&>ftFYa4d806{+Kun-tekueb}HyH=d?v9-K!=Gtr0Vh5E)8s>Z3SXoWC^mqdp4`MAXM>R7976Wb?3HbUB%a(!ouyF-;$CT<Ad3v7u5dg-OK1)nM$2;mgpxv#$aKWZ%6hUDswFSUpa6~#DjWuceEZX3pf!Jj?JM>YVleD~|TR;9pf5qURkM$!)Op?@XKjIDvr44fIGlgDIrp}V-fcoKAN|geRLYupmOT-a0ZBXz#8Iv)F{`oXj&e<^VTs!(i0RM26Nc#Sbax(w9Vk=8AvHBpAn|+qWm&6-o#FTDlOPy!KisqQ*B^Pk2(JYcT1}+ue)<X6oKEv%!FK~oM<Y7UJiL>ZvPYr8~l>iFXd6-hFkJ%uN_jt)`@FFCzv^EvK?NvbT{EY8A8Ep*Fn!jGGM?Z7YLpf(1V0h(3(dqzrsy&|89p0nWY%=kV*A<E_07x75!AxOdA&D`3np@H${Q-3>#*axLgg_~ic8%(#wShQU=(58h4MW{EzNiedI9rs!q~f%H=*MjI?EYVGB>cuAqg>FEIkqv954eCBeC|IH;bKV(rCPM*jas}<=%YnmSk~M7_DzB3-fQ>Pfn25R8wl(F3Np-ux}wR%)1{Z}w|0s&^t=&PbJGY+d!R!?&E(XpuAEo}OhGz+yL3jCb9hNYr2lkvt2dSJCz1eH6S}Sq7L$(nouCfri6~ilaMeUc|7uPaOLHxhg4U#PUY1b~RL3rrBH<j$8wZiJ-lp|d8)QHyq8FC}eiy;vT0m4ixiMlc9qPGw7%GyS#HufOBNotUj_$-4i0^cLhF!QW$JmoPh?DvfL>r${nOsL5xDdC*KH-?)a4T|=2c6?VM)1R#ZebtLhzJb|GZoJ1*x1|KtnPzKHFPls=_2<bq0XL(W!*nU%1lB$%l5~_0k^Tu6<pl2<~@wP0<WAGb{73l?hrx+raJA{6%)Pd7NGr}heS`?KGd1_y5&?t(r#)NxfdxvdPWI~BW!JwS&=>B>#e!62y9mn%Ze>ehcbBEtIst8xp?AWpm&_a*~EjEta-&DRa9a@PD7S+`A6|Gu=W==9|~%h89<Xi&T{Rglq1+bc9R{Ioq5WFvtTkttqL&t+-QW{t)5EIY4P1cp3@@TeRM+~H!z)qZ8g^h4m*LO;W*3b7p&Fmcy=57-nS@WIOSg&pn1|>GiGMT#b8cGI?sapQLJas+j#^hj&n1$MH&vTd+0%3?F;iq4h4u{2m%M(Zg{x(lll^#s#;ARPm%w~gpn&<`s<g+ryl7O2RmM|4Xxrp65WZe6Tm8-6T*@7{4JX|Cd=vI?uaBMZkdfp(HceBMsI`KeM=mYC6q~6XG@CGeoA@jE81=xJ`w+m)t}ZO@FbpVMKLBYy^!&3MMUO@1_UR6AY%voTR&k)%F^1Mut-@YqPy>FWs_0O+FUCYg++|rhjECoh6Fq`d02m7ml?aC(`ozLR}+uStl|e|Sx*gJ+Yz^2#2qvKDLFMlS<w)57F&U33{DxBh$K;Zlu(x2*z4hY?T?!S4gQ76(JlPGppyKP&#kqLBhQMwwMmLA=6~^L{06rFB)FvMsM`e;xH?)8km>#?B7?5v90KTK9tZo=L1&B2*K-`?TFRCsbq3!(X2F;nhS^&>OZC&tvoKoIGa~CNYlP_*i1^@))oR~6sZ>i6y40JsrQTCE1)qi7B`JHd`+<e7&WO^}Hfd!>Ic5~Ls#GFfL&UUO6UJ>yowKNP7LAAQAqRO$nV)D)!tWY>wHm~iR1i`5&pL(Ibe>j*$g#D++uTM<g1Uu(siLv!w8sKUkWX1e0tcYD8HjD@wb)y>K;3;r(?rq*@4uQQBJ0Hg?1lJGC*ef*(f3z@mhk1fPF2gr_%#=qM=*AWnj<~q?sDe&*5WLS>4RU|glu*D`Adz_>0p=>;CrA45c_pK_^0{mC2UO65c;K|#O9%Vv>8KZW&hnTv%?9K`CB!!fl$TuDgFs<+<OJ|Z>tA8*l#9m)cA4Xt$MifPZ<?h8UBY`hD$p9Z-Ih*_|#AIJ-tWYC>iE7|1`wC{O|D7e~i->c4(`d`np&ODA2zL2s?^{(bE_1#Er8sdkZOb15)4pYvZK-frO#Fw+~5$OS17YH=Kz-`ZGeo|NTGd7Hs`=oMGq8h9|)>w=YiWau*bUY^NlvSrHd60S{cLlN{w_^gT18d9e{t=Z(4%N}Y8DC#%<1kp(n2ZyGt(&%n%cjo0gEYN(Sxfn(Gbk4el!Ju{NAi5qOMjNZN6#`r)!bi`nsxC4Tl7l~Tk3(1DRg82oD%BFeYu2&_3#J&ztWh1dbWDgNrwG~bm4?B>F(L-)Sb{s)5TWk>TD|8^lsa*ws#??guFK#;{5WqA7bpL7p<UCszARamKe^#i87_|3LqjlxQ=QI3P$8`M=<u-n-0=s$Q2=3i<Dz$LAwu}05X5DdGz%9P#uvIbv)(0!TV7GXh{-dEjgi2S6Y$e*KGDWZRY*y-b>S`uP3LZ+~4xhJ@)Q|4bL%%(~&$(c(=de%Y^pV|4C>wLVQUGsHoL~wSU+7TZVUaG4EnTG0RvuWzP!sYy{*E~6$naH6A?vh#G!GX&DA7CnIFf4LJ7AKGXTH$6Q1)PX=A;tv>J?ipF%A(X6688(=rB?*meI|OwkWvPQ*>EHDbAHbW`ocXF8H@L1;5F2%p+iex$<k};IGb}^y5Tx_8eBQwoa$bG`iEAWCo92V?Uepg(4GZ-4E3LMaQ_PW||&p_ahemx?jp35C4dO)?+r`y2MWCg6>LOWwt#6dEq~Ju-qYw@G1HaK`k-BWER~zi+N{XTMsyDaX~BQN6F6^{(;6f-K!5$Ewr_=L826t0&&751#8PH-=@024TfXp-AlZ~)QZjL!!4;_W;4=$IfSHfk>Sba`#VW{3@1pXObCm~26?J_jY9}89&Q+nI8ve1{OEQjREj!%9WHbJ<Q51DOi%2;MhFp_O>b%X-JAwF9RK0bf`@6teBjL{N57g<9zk>*K}`Ce?qEXeV9KNzU(=k?`QJdi<Ao|ZMQY}D)3WI^DhCk!x8<KCClGKHt&Bl&UrN46`IKI!5Bq1gn<%Lmi?9)V9^MLc0B=Fxmj6WAZ%{|88ZG%}s%t%I89}a(gs{qtdkm}r7`v3@O1liv#G#i{&K|8%AVHdxT&i{QN*40aATfG7*YG!dR4>+E(zc!X_5Kac7K@p#!dK7wkA`99VP`eo5N8BLvw0RMjuap3h(?>VBEk1~qUShVbQ5xBEb;ANHlVZv13Z!LMD(W@T;!ij{%HjV-<ybMp#7O&tt5QevkB?aL-tSI;8@PE5Bew3mNNv5eH?pA7QB7y)wCtJw(KtFYG&NA6?N<-p8FfX<;fO><W2_0u)^?fX&+xf)n87OxNQmV+9sV1V3A>n-E6xaZo_}>p*_F1rbktwdPS5{?%-hXV>7;F^XZvD!`=~+{VN`GU8mAd=FH`PUl!s+H^D^2Pu%H|_jXL9=7JXva*<S0l%_xqUu?s(^#*W$wemGwn7%H!7$+KG0dG9SmP5INdK<*itFo6{wk6EVTwpvkPj1HBuih_JPlwky0=r^rdYj82m41#P=JcZ;T2m6eFMhC#@ozL%9ETDnoVDTabu32rack!vCSy8grDoEA^KkimH(z+7QGl>MBJ030lNweYUZ?>N!K$vc(LH|xR>FsM)poL@+U)z}%l%Sr^Pq5;?h-LmGhk65S>0sH5|=SYO7-lYVfO~EWUlahtWGf)oRV|af3|l=Ch%p9wrIRAJ~V4FF`NH*J@M<1$T6iTKV-PX_?_ELdpWb8hw<}p5q?|2vZs&D;_4Bf?h>d^*b1UE%w-o`ci!n76Yc3*klP0~p?0yWUC{5oN;%o}Jt7+=@ZK}l%r8DTYZMYhX)1lIXqQFp;QbTxC&T=tEe+avvhi^CR*9hSzIMVqqN+jyZ$_z<9kFC-cklnRU|mT_G1iZez#Jt1qtOO@CiTN7@B?Jq3ja1-6&}J~{i%C!5YIo2#*h4@tKD`dQmOd&ZnNvC#(25-4@tpAvXsh4TCAJx4T-w}rmsP1ty-4P^R5_SX`Ry_D8OF5p^qQ9dAr5b=QN?8xQI_k55_nt-Iir`q=1BASsd*e;R;l+hR>%F+2XtTrJ7x#eN1EB2Q8Xw{%NEfJI&43Ha(xbSC4zf$+vQP(rTkGj}gU?sO@~$9!11<Pdvv4w%i_e9vL)}7((!r8!xL;^CRY>tLR(E|3n_OB7zqD-Qt5BOYE@i8`?#QQ<j7dBI(Ks04pELaY-!Fl{EDGOG@gQe`5Ib`-i-)p%{PAY4K(I^>yXa#oG8x87-G^(IU%I1i%u&=>`94mC^aZ8lhN#FqD$D=MDENGGPHC&2JrQ+dqmAG?UhVnwAR)uer(mxvrmu%v15{nu?UUmT*cSED&y5T0q`4>y_OQ0Nb0{pO3Nl`*&ELro-P{UGubGSmSnrz90HDE)sP&eU<Gp%MhnHTpv@uMM0DqD=T2XRBPi%x{rBNQr9FY!9ldQj<V`o_gexQ*C9+paVu~JNFZrn3<dm@jFvjU4zqWQAO`jdY%zZ0Vp~(hOrxqKr*=lEH6!0M;QR$yglA8Qp@Aj)E5{&dKDadNA<Y>z>8H{6=#@#EqJ~Wi{{UAKEG7jCxmHZi8q}4CuX@Qx{bBkqZr@udB27UV4I`FI_p#bwWDWiB{}9+jH**&|pG3Vj3LiqfgcY`>|Kdw8`D`l19>%l;al-?clLr~oN9TB*j18jvUPp|hxPzHunp#`}&*AxG=Z|Ybo-)1_^!0wC8aS$0Cy2mwS&Sxk;8hbtqjpZd5$XYN=~j3}*tfD$EWbnVs0-REo_eWLx*!76!<e6vuh$1=q!^-^1?#y`)n&?pups~0TtlGm-siTk{X<KIHF?@-`~Gc59d55m_pc=HocsSKAz2|pFppv&Ds>?OA`Mq`+xQwRMDL>i()>-NXeMGG7Q9!Cjx+`8#eQC9EG8A)(6%`^JAHPynH-j&^;U(THWK7=g;~65pg<zpKuI4;i<2qMw$taMl7Z*IcEfCz&W5x-$gVDwjOvYl0-j9~W&fKCV#SD&W!E-WCqVNW8i%{4ria{MS2q%U*Z0sSeK2oTWwq>~&s}CY*91UQ;l)PeF7wF;ppnB`om^6?vC=<(86QD3dLl3m{R2=s*l7bHmj#A>MDDBm1hYh482EpgBG8=7#=LOZh6VS1RYU)>S&T_uZpA*_%>AN*ho~38Lpg#0p^Zsj|L<Ual}<N;Z`Mk8*;}q^_Lt+snmbFRvfBz=8;d01`t>DMKcQvAt3HcgF0{Q$A`!!@y*=S4&83Pa1zHTfMy9&l3Oc!#7zHKQbBIy1gZ=OQeF7%VElOAAltOTBMXgR%3SIYs9vDe0;x6entQ+^XwTQH)e)&VHf8|wNnF(v9pU6}fVMefG(R-sfOufei?+@^NxHPFw)~teDUJ@!s40O294;SIl?|3vp@xMx<L9@*BXVc4+ysY}59n~fKG1a|uC>F0@_<N67Qk?W5(;g9x|9@rC&eP&nRXuFw@Qnn0umWu^vySKkHcxG1iVnPdrrTO%nJSDp=GQ*$7?kMe%*;Zv?Hyi|-jwp4vTN71!&Pd@_RTxp*sc+<Yz@SR2ta2?|GkG25s~T%VSt~Uww$QsD(H9i)Zb68YB0MjHWG=UA?P9bH=+nk>yuI{HEys08uKy`-@D)_lsv%D8!Cz>WYgK7p6o2`I+L`Ax1Na2rZ4_KY%uKpxSvl>gO(T0_dRxEpD2I+&rR6Fw`LlT_hvYRc_9FZI_R13{YH`Et(TWuEV>)}j7*kJy#*ACN!-ELLeD<OQQpim{^G1;LqEHSJcFr-KSO=c;hCwBH1xRJX7Gr4t@|)w>gKCW%h!%ynx7LrXekz^mTiVCxK7|H^Dl$?!gt!o;&jLAuPz7YWqP%aHqLjZrzY@K%nS@>yT60=)DSmQ1esxqUUD!inXvT12D!GLJyb0r(ClExb7qmFX;It`H{QU+2Ua1)_0?~g|7(n4CZXmo5$y5VPPMbG=}i9qzKemV<kh$B8z>|(&eJQz9PWuWEAU0g{m=6>0OQbf;yzn!$TF^|SYT*s`TW{eNa&j&b#-kWsKHaHNcxQw@zHcdf7Q~Z*xV3y2WNMYg&+^Sg!Ud`=-~ci=&{$3ygi#`=Q&M9d15F{-SL(X'
