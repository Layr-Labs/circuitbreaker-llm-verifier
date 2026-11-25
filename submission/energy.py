"""Active -291/800 first-pass, gain join and starting-majority proof."""
from fractions import Fraction as F
from itertools import combinations
from functools import cache
from math import comb, factorial, isqrt, prod
THRESHOLD = F(-291,800)
def ck(q,m):
    if not q: raise ValueError(m)
def fifth_order_orthant(K, det, signs, scales, threshold, center, pi):
    """Lower bound an eight-dimensional Gaussian orthant using exact moments."""
    scale = 1 << 44
    size = len(signs)
    if size != 8 or len(scales) != size or min(scales) <= 0:
        raise ArithmeticError('eight-site chart')
    v = [F(s * sign, 10000) for s, sign in zip(scales, signs)]
    sums = [sum(row) for row in K]

    def floor_scaled(x):
        return x.numerator * scale // x.denominator

    def sqrt_floor(x):
        return F(isqrt(x.numerator * scale**2 // x.denominator), scale)

    # Encode a multivariate exponent in four bits per coordinate.
    units = [1 << (4*i) for i in range(size)]
    terms = [(0, floor_scaled(-threshold**2 * sum(sums) / 2 - center))]
    for i, u in enumerate(units):
        terms.append((u, floor_scaled(-threshold * v[i] * sums[i])))
        terms.append((2*u, floor_scaled((1 - v[i]**2 * K[i][i]) / 2)))
        for j in range(i):
            terms.append((u + units[j], floor_scaled(-v[i] * v[j] * K[i][j])))

    moments = [1, 1]
    for j in range(2, 11):
        moments.append((j-1)*moments[j-2])
    mu_lo = sqrt_floor(2 / pi[1])
    mu_hi = sqrt_floor(2 / pi[0]) + F(1, scale)
    powers_lo = [mu_lo**j for j in range(size+1)]
    powers_hi = [mu_hi**j for j in range(size+1)]

    powers = {0: 1}
    expectation = F(1)
    for k in range(1, 6):
        expanded = {}
        for monomial, coefficient in powers.items():
            for code, value in terms:
                key = monomial + code
                expanded[key] = expanded.get(key, 0) + coefficient*value
        powers = expanded
        sums_by_odd_degree = [0] * (size+1)
        for code, coefficient in powers.items():
            weight, odd = 1, 0
            while code:
                degree = code & 15
                weight *= moments[degree]
                odd += degree & 1
                code >>= 4
            sums_by_odd_degree[odd] += weight * coefficient
        lower = sum(c * (powers_lo[j] if c >= 0 else powers_hi[j])
                    for j, c in enumerate(sums_by_odd_degree))
        expectation += lower / (scale**k * factorial(k))

    exp_lower = sum(center**j / factorial(j) for j in range(16))
    if expectation <= 0 or exp_lower <= 0:
        raise ArithmeticError('positive fifth-order minorant')
    return (prod(abs(x) for x in v) * exp_lower * expectation /
            (2**size * (sqrt_floor(det) + F(1, scale))))


def energy(rr,uu,M,G,pair):
    def ck(x,m):
        if not x:raise ValueError(m)
    rows={''.join(map(str,k)):(F(l,M['S']),F(h,M['S'])) for k,(l,h) in rr.items()}
    U5=F(uu[1],M['S']);t,a=F(41,57),F(28,57);u=-THRESHOLD
    @cache
    def r(k):
        k=abs(k)
        return F(1) if not k else a*t**(k-1)+t**k/2+(F(k-1,2)*a*a*t**(k-2) if k>=2 else 0)
    @cache
    def gram(n):return [[r(i-j) for j in range(n)] for i in range(n)]
    invdet = G['inverse']
    def cone(n):
        K,_=invdet(gram(n));s=[(-1)**i for i in range(n)];D=sum(map(sum,K))
        A=[[s[i]*s[j]*K[i][j] for j in range(n)] for i in range(n)];b=[s[i]*sum(K[i]) for i in range(n)]
        best=F(0);count=0
        for k in range(1,n+1):
            for face in combinations(range(n),k):
                C=[[A[i][j] for j in face] for i in face]
                J,_=invdet(C)
                v=[b[i] for i in face]
                z=[sum(J[i][j]*v[j] for j in range(k)) for i in range(k)]
                count+=2
                if min(z)<0<max(z):continue
                # Simultaneously negating v and z leaves the quadratic value fixed.
                q=sum(x*y for x,y in zip(z,v))
                ck(q==sum(z[i]*C[i][j]*z[j] for i in range(k) for j in range(k)),'face')
                best=max(best,q)
        ck(0<best<=D,'cone bound');return D,best,count
    @cache
    def alt(n):
        lo=hi=F(0)
        for k in range(0,n+1,2):
            for S in combinations(range(n),k):
                sign=(-1)**sum(S);g=''.join(str(j-i) for i,j in zip(S,S[1:]));g=min(g,g[::-1]);l,h=rows[g] if g else (F(1),F(1));lo+=l if sign>0 else -h;hi+=h if sign>0 else -l
        return lo/2**(n-1),hi/2**(n-1)
    def ep(z,k):return sum((-z)**j/factorial(j) for j in range(k+1))
    def arct(x):
        p=sum((-1)**k*x**(2*k+1)/(2*k+1) for k in range(24));return p,p+x**49/49
    p,q=arct(F(1,5));x,y=arct(F(1,239));pi=(16*p-4*y,16*q-4*x)
    Q=10**24
    def down(x):return F(x.numerator*Q//x.denominator,Q)
    def up(x):return -down(-x)
    def iv(x):return (down(x),up(x))
    def add(x,y):return down(x[0]+y[0]),up(x[1]+y[1])
    def mul(x,y):
        a=[v*w for v in x for w in y];return down(min(a)),up(max(a))
    def scale(x,s):return mul(x,iv(s))
    def div(x,y):ck(y[0]>0,'positive denominator');return mul(x,(down(1/y[1]),up(1/y[0])))
    def sqrt(x):
        ck(x>=0,'sqrt domain');k=isqrt(x.numerator*Q**2//x.denominator);l,h=F(k,Q),F(k+1,Q);ck(l*l<=x<h*h,'sqrt bracket');return l,h
    def pinned_orthant(C, sg):
        n = len(C)
        block = [[M['ONE'] for _ in sg] for _ in sg]
        prob = M['ONE']
        for i,j in combinations(range(n),2):
            rho = M['div'](M['rat'](C[i][j]), M['root'](M['rat'](C[i][i]*C[j][j])))
            block[i][j] = block[j][i] = rho
            pair = M['moment']([[M['ONE'],rho],[rho,M['ONE']]])
            prob = M['add'](prob,M['scale'](pair,sg[i]*sg[j]))
        if n==4:
            prob = M['add'](prob,M['scale'](M['moment'](block),prod(sg)))
        prob = tuple(F(x,M['S']) for x in M['intdiv'](prob,2**n))
        ck(0<prob[0]<=prob[1]<=1,'orthant')
        return prob
    def line_moments(n,D,beta,degree):
        C = gram(n)
        s = [(-1)**i for i in range(n)]
        @cache
        def block(pins):
            keep = [i for i in range(n) if i not in pins]
            Kp,det = invdet([[C[i][j] for j in pins] for i in pins])
            weights = [sum(row) for row in Kp]
            v = [1-sum(C[i][j]*w for j,w in zip(pins,weights)) for i in keep]
            d = D-sum(weights)
            return keep,Kp,det,v,d
        @cache
        def base(pins):
            if not pins:
                prob = scale(alt(n),F(1,2))
                return prob,prob
            ck(len(pins)%2==0,'even face base')
            keep,Kp,det,_,_ = block(pins)
            cond = [[C[i][j]-sum(
                C[i][a]*Kp[g][h]*C[b][j]
                for g,a in enumerate(pins) for h,b in enumerate(pins)
            ) for j in keep] for i in keep]
            prob = pinned_orthant(cond,[s[i] for i in keep])
            norm = sqrt(det)
            for _ in range(len(pins)//2):
                norm = mul(norm,scale(pi,2))
            return div(prob,norm),prob
        @cache
        def moment(pins,k):
            if k==0:
                return base(pins)[0]
            keep,_,_,v,d = block(pins)
            result = scale(moment(pins,k-2),(k-1)*d) if k>=2 else iv(F(0))
            for j,w in zip(keep,v):
                result = add(result,scale(moment(tuple(sorted((*pins,j))),k-1),s[j]*w))
            return result
        out = [scale(moment((),k),2) for k in range(2,degree+1,2)]
        deriv = iv(F(0))
        for i,j in combinations(range(n),2):
            c = C[i][j]
            prob = base((i,j))[1]
            deriv = add(deriv,scale(mul(div(sqrt((1-c)/(1+c)),pi),prob),2*s[i]*s[j]))
        old = add(deriv,scale(alt(n),D))
        ck(max(out[0][0],old[0])<=min(out[0][1],old[1]),'second moment regression')
        ck(0<out[0][0]<=out[0][1]<alt(n)[1]*n*beta,'face moment')
        ck(all(0<lo<=hi for lo,hi in out),'positive moments')
        return out,n*(n-1)
    tails=[]
    witnesses=[
        ((4750, 3131, 3354, 3376, 3408, 2928, 4745, 6634),625,(4398, 3215, 3330, 3406, 3380, 3016, 4635, 6201),424),
        ((4762, 3107, 3374, 3352, 3312, 3361, 3411, 2921, 4741, 6630),832,(4400, 3192, 3340, 3388, 3290, 3396, 3375, 3008, 4627, 6188),668),
        ((4766, 3100, 3378, 3335, 3334, 3342, 3322, 3349, 3411, 2923, 4743, 6634),908,(4397, 3183, 3341, 3371, 3303, 3382, 3290, 3386, 3371, 3007, 4624, 6183),808),
        ((4765, 3095, 3376, 3329, 3337, 3324, 3342, 3330, 3323, 3344, 3406, 2926, 4740, 6633),1025,(4391, 3174, 3336, 3363, 3302, 3364, 3302, 3371, 3287, 3378, 3364, 3005, 4614, 6170),1028),
    ]
    for _ in range(3):
        p,z,q,w=witnesses[-1]
        j=len(p)//2
        witnesses.append((p[:j]+p[j-2:j]+p[j:],z+175,q[:j]+q[j-2:j]+q[j:],w+215))
    for i,(sp,zp,sm,zm) in enumerate(witnesses):
        n=len(sp);K,det=invdet(gram(n));sg=[(-1)**j for j in range(n-1)]+[1]
        q=G['lower_signed'](K,det,sg,sp,zp)+G['lower_signed'](K,det,[-j for j in sg],sm,zm)
        if i == 0:
            w = (fifth_order_orthant(K, det, sg, sp, THRESHOLD, F(-3,5000), pi) +
                 fifth_order_orthant(K, det, [-j for j in sg], sm, THRESHOLD, F(-2767,10000), pi))
            ck(w > F(37352398,10**11), 'eight-site lower bound')
            tails.append(max(F(q,G['S']), down(w)))
        else:
            tails.append(F(q,G['S']))
    k2=(1-r(2))/(1+r(2));klo=sqrt(k2)[0]
    sint=sum((-1)**j*sum(F(comb(j,l),2*l+1)*k2**l for l in range(j+1))*u**(2*j+2)/(2**j*factorial(j)*(2*j+2)) for j in range(6));ck(sint>0,'bivariate lower');sdrop=down(klo/pi[1]*sint)
    D4,b4,n4=cone(4);D5,b5,n5=cone(5);D6,b6,n6=cone(6)
    degree = 6
    M4,f4= line_moments(4,D4,b4,degree)
    M5,f5= line_moments(5,D5,b5,degree)
    M6,f6= line_moments(6,D6,b6,degree)
    def upper_delta(n, precision_sum, beta, moments):
        total = sum(h*u**(2*k)/factorial(2*k) for k,(_,h) in enumerate(moments,1))
        k = len(moments)
        term = moments[-1][1]*u**(2*k)/factorial(2*k)
        ratio = u*u*beta*(n+2*k)/((2*k+1)*(2*k+2))
        ck(0<ratio<1,'remainder')
        total += term*ratio/(1-ratio)
        exp_bound=ep(precision_sum*u*u/2,6)
        ck(0<exp_bound<1,'exp upper')
        return up((exp_bound-1)*alt(n)[0]+exp_bound*total)
    A3_delta=upper_delta(4,D4,b4,M4)
    e5=ep(D5*u*u/2,5)
    ck(0<e5<1,'exp lower')
    P5=alt(5)[1]
    cosh5 = sum(l*u**(2*k)/factorial(2*k) for k,(l,_) in enumerate(M5,1))
    A4_drop=up(P5*(1-e5)-e5*cosh5)
    A5_delta=upper_delta(6,D6,b6,M6)
    c1=U5-sdrop+A3_delta+A4_drop+A5_delta-sum(tails)
    G0=G['compute']();gain=F(G0['gain_lower'],G0['scale'])
    motif=G['edge_motif_gain'](); bonus=F(motif['gain_lower'],motif['scale']);cost=c1-gain-bonus-pair
    ck(U5<=F(73670242301,2**39) and gain>0 and bonus>0 and pair>=F(47991230402,10**15) and bonus+pair<F(1,1000) and 0<cost<F(12246859893519800000,10**20),'active energy bound')
    ck(f4+f5+f6==62 and n4+n5+n6==218,'face accounting')
    return {'energy_upper_1e20': -(-cost.numerator*10**20//cost.denominator), 'eight_site_discount_lower_1e12':int(tails[0]*10**12), 'first_pass_upper_1e12': -(-c1.numerator*10**12//c1.denominator), 'gain_lower_1e12':gain.numerator*10**12//gain.denominator, 'motif_gain_lower_1e12': bonus.numerator*10**12//bonus.denominator, 'pair_gain_conditional_1e15': pair.numerator*10**15//pair.denominator, 'motif_oriented_lower_1e12':motif['oriented_lower']*10**12//motif['scale'], 'motif_cases':motif['cases'], 'U5_upper_1e12':-(-U5.numerator*10**12//U5.denominator)}
S=1<<50;H=1024
def up(p,q):return -((-p)//q)
def ppi():
    def tan(x):
        a=sum((-1)**i*x**(1+2*i)/(1+2*i) for i in range(20));return a,a+x**41/41
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
    r=(2*k+1)**2*S; den=8*H**2*4096
    l,u=ex((r//den,up(r,den)));a,b=P[-1]
    P.append((a+phi0*l//(S*H),b+up(phi1*u,S*H)))
ERR=up(5*S,24*H**2)
tail=up(S,2500000) # Mills <= phi(5)/5 < 4e-7; checked by endpoint bound below too.
r=(5*H*2)**2*S;den=8*H**2*4096
_,erp=ex((r//den,up(r,den)))
if phi1*erp>5*tail*S:raise ValueError('tail bound')
C=[(max(0,l-ERR),min(S,h+ERR)) for l,h in P]
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
    weights=[2*pm(F(i,STEP),F(i+1,STEP)) for i in range(5*STEP)]
    negative_upper=[C[0][1]]+[S-v[0] for v in C[1:-1]]
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
    if not (0<y<F(29,1000) and 0<x<F(13,1000)):raise ValueError('pin triple enclosure')
    M0=M['S'];Q=M['rat'];ad=M['add'];dv=M['div'];rt=M['root'];mu=M['mul'];sc=M['scale'];asin=M['asin'];K=M['K']
    signs=(-1,1,1,-1);R=[[M['cov'](i-j) for j in range(4)] for i in range(4)]
    a0=Q(1);sl=Q(0)
    for i in range(4):
        for j in range(i+1,4):a0=ad(a0,sc(rows[j-i,],signs[i]*signs[j]))
        h=Q(1);ids=[j for j in range(4) if j!=i]
        for j in ids:
            for k in ids:
                if j>=k:continue
                c=Q(R[j][k]-R[j][i]*R[k][i]);d=rt(Q((1-R[j][i]**2)*(1-R[k][i]**2)))
                h=ad(h,sc(mu(K,asin(dv(c,d))),signs[j]*signs[k]))
        sl=ad(sl,sc(h,-signs[i]))
    a0=ad(a0,rows[(1,1,1)]);penalty=F(a0[1],16*M0)
    if not (sl[0]>0 and penalty<F(121,10000)):raise ValueError('0110 log-concave penalty')
    H=G['S'];half=G['lower'];extra=F(half(('11','11',''),[2500,4800,6200,4800,6200],500,-1),H)
    tri=F(half(('0','0','0'),[1236,3201,3201,2560],854),H)
    if not (extra>0 and tri>0):raise ValueError('forced five-pin and overlap')
    out=F(Fcdf(-T)[0],S)-y-2*penalty-4*x+4*tri+extra
    before=out
    if before<=F(504,1000):raise ValueError('combined calibration margin')
    out-=y
    if out<=F(1,2):raise ValueError('initial strict majority after motif repair')
    def ceil(x):return up(x.numerator*10**12,x.denominator)
    return {'initial_lower_1e12':-ceil(-out), 'initial_before_motif_lower_1e12':-ceil(-before),'Y_upper_1e12':ceil(y),'cross_pair_upper_1e12':ceil(x),
            '0110_upper_1e12':ceil(penalty),'overlap_lower_1e12':-ceil(-tri),'11011_putback_lower_1e12':-ceil(-extra)}

from bisect import bisect_left
from fractions import Fraction
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
        taps.append((st,tuple(int(a*d) for a in w),d))
    return nodes,taps

def expand(values,W,X,Y):
    rail=[]
    for start in range(0,len(values),W):
        row=values[start:start+W]
        rail.append([(512*sum(row[i+k]*x for k,x in enumerate(c))+d//2)//d for i,c,d in X])
    out=[]
    for i,c,d in Y:
        den=512*d
        for z in zip(*rail[i:i+len(c)]):
            out.append((sum(a*b for a,b in zip(z,c))+den//2)//den)
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
    h=[min(Q//3,max(0,v-offset)) for v in raw_h]
    m=[min(b,max(0,v-offset)) for v,b in zip(raw_m,h)]
    return [3*r for r in m+h]
from base64 import b85decode
from bisect import bisect_left
from lzma import decompress, FILTER_LZMA2, FORMAT_RAW

def read_words(payload,count,raw=False):
    data=b85decode(payload)
    if raw:data=decompress(data,format=FORMAT_RAW,filters=[{'id':FILTER_LZMA2,'dict_size':1<<26}])
    else:data=decompress(data)
    vals=[];value=shift=0
    for n in data:
        value|=(n&127)<<shift
        if n>=128:
            shift+=7
            if shift>14:raise ArithmeticError('witness item size')
        else:
            if value>=65536:raise ArithmeticError('witness item range')
            vals.append(value);value=shift=0
    if shift or len(vals)!=count:raise ArithmeticError('witness shape')
    return vals

def inverse_mixed(mixed,w):
    row=[]
    for k,t in enumerate(mixed):
        a,b=divmod(k,w)
        t=(t>>1)^-(t&1)
        if b:t+=row[k-1]
        if a:
            t+=row[k-w]
            if b:t-=row[k-w-1]
        row.append(t)
    return row

def from_knots(main,contrast,Q,step,off,blob):
    x,y,_,_=shape(main,contrast,step,step)
    w=len(x);s=w*len(y)
    vals=read_words(blob,2*s)
    a=inverse_mixed(vals[:s],w)
    b=inverse_mixed(vals[s:],w)
    return smooth(a+b,main,contrast,Q,step,step,off)

def decode_coarse(main,contrast,Q,step,tile,knots,losses):
    w=len(main);h=len(contrast);n=w*h
    tiles_x=(w+tile-1)//tile;tiles_y=(h+tile-1)//tile
    cuts=read_words(losses,2*tiles_x*tiles_y)
    point=from_knots(main,contrast,Q,step,0,knots)
    split=tiles_x*tiles_y
    m=[];p=[]
    for i in range(n):
        j,k=divmod(i,w)
        at=j//tile*tiles_x+k//tile
        nonpositive=max(0,point[n+i]-3*cuts[split+at])
        m.append(min(nonpositive,max(0,point[i]-3*cuts[at])))
        p.append(Q-nonpositive)
    return m,p

def decode_fine(main,contrast,Q,step,offset,knots,losses,rem=b'',snake=False):
    n=len(main)*len(contrast)
    point=from_knots(main,contrast,Q,step,offset,knots)
    cuts=read_words(losses,2*n,snake)
    if snake:
        w=len(main);h=len(contrast)
        for base in (0,n):
            prev=0
            for j in range(h):
                at=base+j*w
                for k in range(at,at+w):
                    z=cuts[k];prev+=(z>>1)^-(z&1)
                    cuts[k]=prev
                if j%2:cuts[at:at+w]=cuts[at:at+w][::-1]
    a=[x-3*v for x,v in zip(point,cuts)]
    if rem:
        rem=decompress(b85decode(rem))
        if len(rem)!=2*n or max(rem)>2:raise ArithmeticError('fine residual remainder')
        a=[x-y for x,y in zip(a,rem)]
    m,p=a[:n],[Q-x for x in a[n:]]
    if not all(0<=x<=Q-y and 0<=y<=Q for x,y in zip(m,p)):
        raise ArithmeticError('fine witness simplex')
    return m,p

def lift_law(coarse_main,coarse_y,coarse_d,fine_main,fine_y,fine_d,law,Q):
    m,p=law;w=len(coarse_main)
    top=(0,Q)
    cols=[bisect_left(coarse_main,-(-i*coarse_d//fine_d)) for i in fine_main]
    shifts=[13*(i*coarse_d-coarse_main[x]*fine_d) if x<w else None for i,x in zip(fine_main,cols)]
    scale=41*fine_d
    fm=[];fp=[]
    for j in fine_y:
        Y=41*j*coarse_d
        for x,dY in zip(cols,shifts):
            if dY is None:
                v,u=top
            else:
                y=bisect_left(coarse_y,-(-(Y+dY)//scale))
                i=y*w+x
                v,u=(m[i],p[i]) if y<len(coarse_y) else top
            fm.append(v);fp.append(u)
    return fm,fp

COARSE_KNOTS='{Wp48S^xk9=GL@E0stWa8~^|S5YJf5;5aK0yIlZf6#f9ZTq-@7c5+Bm$u#&ysp}0NgbBqiw1u|KWyh~91P`;s$=te?VC7;3Y0ZDNblwLHLn$`Y9U!#^!{-tq?GwRO2W}Mg#ISR9(>Uo{O1zf5e=AfgpBhYsQ7~&|K6xbf3SN&ZK98yt2Is!2T#OIp^TV6xUp45w>k>y*k-@r>!y#LGXVWS6`(Lb$!N4p;GH|V4%=kRS!q#bgnGc4;4F6Q}q)rqYS7Mi`V<)w{MuQsFyU+2I;$YE|cKjMmQr1dENpaTuPWs23)+J#3_&)E@iW|(Mt)Us=RxZDz&`YzsLe}e}o9M^eI@Y!py1w2R4^(CUk~bP2eb|q9$o|sYfESMl<_wqp>!5Ix<T42mnalUA^Cm^8P9*X3KtJR$Je0qMUQdrdeWjXmc&rNi(78h$i&@5U%MsA>Q*h34s5$*%<16bZESZq{lMM80B@xtKH*_}JUfTf{rjQ$6CuhS(@EVIZuJ*U%t=)Yt-4rn((7&0i73`31C!Em7VSVm$q|;3CIf3(&wRGc$7erdXLw9(n00R7_z_p?B`I@Qy6r?_#BWoqX7o=rKQbn02woO7bbPAl=u*G;_|LW|Q{8i{xD~~|<u#tCW>|N6w?c7QQSAy4auX*ybqnvK7dgCPs)~r)bD>he;<KwJuh}wn3^{7_bq!OF{*i?Pu@}G~MrmMp)%X9eoSQKqBZ~G87Bd`9EowK@o+9uB>O?@bD<*i)>^f<4P!`Lj~hdp|e&sAu-<L(G-3||cO0(n=gim?YJWv|`rwo$(I9uHUkk#DSh?NOLb(J1a}=(V{bd-uu{P<?)pZ8=aD$E>nxJd}o@2w|knp8%!V9yQn+0crzclukr8i@2G4=>9-jJJB>v5hRUIw2U9@*K(I6`-{6Q=GcEb-El~@fbVFTrryQ0t|3eV5kw7|d1{7CA0U9(sidyZvuNsF$dm&D->{0C#N-C!jjH;rYJ&PE@|L`Wavv~S%OohMalZn)LntXVO?mk9&{Dz}fA~8Ro+$5D-2G0N+3<n2T%W!+`lmh#N+ev~lUScK2#ky?H3SB+lCsnc%_Mz5k*$<6k6^M~-lShY_-B!A#}*tx?!nv^z<MXPbpBz?M-GQ<2HKEqK}L%5k!&6i)^e#7D0+vJIS4;C&KMS%1qMmw%+oYhN1H0ndbM3*MxvzTB{=4Djcs@H9M{?>Sb8n_%ATIaYs&GUWXbCC3(LWz3r8RI0O|h=ElX#ZUa~qMbX2035jB?0VC|}JJ)RyR06+>Dc|_B$wLrwx^1iC;B^f6Y=Sfpm@oh;3zl+biffDD+5SgA|tO_h02ktQ;_l|T-P99vjW?ir_6yRTZ75E5rMrQ(S^5F9S1+M5TzzfyzP{FsWw_2d4E5j=yWD4uAJJ!U`<MzO&TF5Q9Vfl_hd4ExXQ_|p^jXwK`;WwM_2VovMppPGCmqKbyMXQw?qsKR8rE}t%rI-ZeMc#l1TsxKMmqrp*W;vNg-nHBm%pw?ugk;5*WEtmt-X^TAZq=vDexNj_{Qminxf4q+@lh7?)3<n}f|$fcbFd$ZRuFGu5umH3f_*Q{#aIfLEV$fQjh)ZlZ#?<{;=27FSbbk-*dF4_83--j`OwPlv$i=k>xg-j1W@5RXu8oO_vAl4U?crP8Wq%JHH#r6hNK|2fTl4;XYYtJHp-|H>HXu*R&y;9ZR%>UuH)==n*1UegX5O*yIszGIr0vRiUDbc?5L7KJ?8nav?L>P+sTzDvF&;Tgb%^|ru`v&g3*guJ1S2^xv2_Eb)&3j-(9aR=j}!eV1THa<ft}WZ+SQHX&SqUYW!Zs3~H6x3|SQ^WfVG?pokp=y9xI$+>`K*H>6wdczbDq2guONxQy$2IATwO)I)jqKL!^y_8jtDP=2Kym*#0C>=AO3O_`a^&_J9;_@(;3|J&8j3%!%tIp%$m)W9c%lO(g<h}kNMaO}Bf${Y!EzuFm}kYT;4f$FgW69>yDWIP>?`nqlr;fo;;Cow+E1>!O=oLh+#1ff&^6gc|k1}t}w+!y^TMwu9?vyz3u-vct+pt}Ti0fu-&n1S;(_8L?5zL2}Jt|Dr+KVc7An2n)y!*7aA8nUV9ikA{E^I;BhC^~~*W44>sgET~eCjP)>5%c1!w5}O(yg$4AWIa4HQjm`a_J|KgyI}1|o3Fs<N`8R*?=TO#9dtzDV-amU)6JtuF{hx(@Gkhu0AU9eA8nTIz7l9U@G-gNo2!WfA6a#-ugk?KL5kb5Gf1^|j)h2%s#p-8)b^s<rRFELXJjdo=4XTLHMn1_h{Q%WY}v$iho(aLeDnz?>;g3o@!%bSqM3Jq;>r^xh~iWOn8qh~b=-i_Kl`<j|8chSA<PT|xIn&N_>|V~8($He2xS@Y^{_0WbVA-<RpLGMyI$L;w{{R=U6=oyv*DAgrc6ju))RgpZm$VPL`L(^SDMB~?u<R8BiuH(Tr9nH*OvB-&DP@vf!NW5zbcc7=~j)6MCAa<udBVn;|<b4aEm`axQ(oqwS%i8@=qj45^dS?s+R+n$~PI$VVg_Kb$kju_?P8tuXZ?&TU*XCN@{vPYx%iSpdH5;3~<egNR6Dg%6Cg2oL3Q`6b2Q0>ij*A_Rrjc9zwoAZtU#9@VT8;IjmZpfn27>03@rXy_nw<U5<Ku%7)9+cKKq;(mRyGc8Je*hWi&R+%Rk(q>d<vpGyLdCGLR~U(x~>{zT@A2*A-!4TFzE`<I6=-!q<~r|RQrQCZ1?Q9dUj(~|ZUD)l`NMZ^FikR6<DS-#7IDK(PmKmP0CPoqwMeBB#o@mlEh`5a?!OqIOQb<(tgH)p+x1?H>5pW6K~Q8dl%-NZ;9ZyVZ%htkngnLdG$pJqdDP2c(_D!ud>Lsd9_4#?XEOVM#eMT8D>`RMx>I)m?$+`NLkyCX}(I>s5{EE<)i<@m-6q<#Z!6%Yb??x651qOmWTwPG@^!EnnOj71j+QNI3^*%xa7_-I+8p7@PMcz$&kRUKe`zS84w$gRPuQ(6iuER>Dq0CQ<2S3}K2rTr5mR<iT(bSdMF5<eWwnWqTu>Paowb?jRv<g&^t6!P;5kRJy*=F-8qoUVg!Z!{3(2oH8cnP|c%o=5+W-e{rOOw7<NH3I>c4z*m;C`y^O4!O~%)w!?ZOUBC;b*0wHfHQ@AC$*uYAAulrI~gxdGNVnSj)vQQqJ1S|zV0bcn5SQ-C$O}f>f@xvm??lW(9eD?$O1Xx$(_l}Pe${$DfK>kVE6M5LQGFg7p>cmr8E6fS`Qt<P`rabjs(!pagdI|=RZ!2z<f;d=VMa~Tzu1-1`-b~ItqoPy(t(RQm-Zz8TUSf-sKzNQXfs^?1J>r!oXABWu<?n4GSNUyqlmUUKJONO7~E_N368(h*126!pKz<C7>brBp^nUbng4-Q&C)Rf;Ph0#MtK}cRyEQKjQo$EDYJmqoYX35NasREmY5W1k~29=P=V`OGCgw66v%b-$iQ-p!DY;kozr+adMGbIYEREYM}lNEY{4E3&zF9<Vy+aQL`?`plm}Y-M}2&$LPEs<!WAUj4X?a$#>m?oF4=BT87E|WJi=XgL+(%fW0n7z`^w?h5IT;cq%W9=xI=|OHbCj>`4cQv;aSgP)>iKvVTC_cLj7)7Z*pc+a8fZid|x<D9ZA!lLIkvc_zZ0+y>-q-O?LkB5C<BckWt=D^$Y?P8vNf^^`}py~CLPQI!Yb`I*u>a}i;7U;$z6CU=3=8q1GwV$MQ_Kk+Bs6l(LLXfu@rK^qcNVr;PpC(*4IUc@lHW5nsTbV~;d++)&)@fZwt+JRmxP)E+jnE)j<NpjV@jq7BIYTl^3YUwqT<XdK4fQRKuN6|G`eN_Q41wQk%u3X4?>khJ`VEu=lfn2mTh0bT!g)_RfG}Ocg>Oo%sn6Hn!s&8Rj1;5{bG#FjeR`)gW=KdAQ@KoS1@hQ+<C;es$Jiw2RWrsU!_l*$#k;;6Q{Z|_fAeNj;(<&PfMp;fRWTrH8eCplF^~v*|;Sh)`bY0*`r?fN^GLEfd8l1@Tn4to>ODUh%?GT<Mcugqt=UR#4+tuMW?7VQyb>cvsr@%JkVXtLRBZq2KTndld9t}6{+yKkwTc-`BHJWEbj9`nab3{5ZTD$&-RAB&VVHN>El4rVxd`o&ps<74j(BjS2n^{sJ0H?N$evdV#$R6T*7zS&Wi^7%eiAsay_vi}5nV^mk)3W=CWX&pu)G8b9^DrIFjI+7E`l(oN@>;Gh2O-l#p|atdW>CuaHY=(yC&F$}12xXL4tp;qp(JA|1`%zWby4iz^mfK3AI|0a0a**@OLI~%@SkplT^-EvM}8iHR_-a03Z0@>a*b`l6cBG45|VSJ5IFN=Zh{l)0GPS{DokZH6Qaw?z`%0fZD<yrx?@(3VGBRR?}jd9cdUyqJPc|4r|wCdgi*Ha8%Z`WiAxlG7dZxY%-Dp0p1pl4)(xYD!+uanjzCz|*g$o1j`lOi#GEYA>OR0+m(jDhS(}<<iPCfJP*?5ed;u(gE5eTUbui5PB5e{$w{6Fm3FE4jZXz~UnJB^1wOUo-jr~3LN<UG@D>aGounIo80Z9R7+jU7F)Hikd%4C?(TZ7C6>H}?UytP+n;W5#{e(u1vZ<iGKWNhyy@p|X=<Z@+$JW5bC(3qn8FHn&6e^X8G?qMha(|K3MD{3<+Op&i_M4a`mjKE-}>^8%0X!gHN8$m~bip6Vhui_$27Pfd{4!(TR+#t-P`XUG=ZUTn8<`vYbDMf9yn#NA<lldY|vpdWRbsFWZNOu%kHdSkM`?zK@AO>x;$yiT6y_AvEIfX6ov=pvNZ5}KXOB(6$kDo)8Dv(;%TwBvghACe<H6GWOI?t>TbFCA}8{K+5#6scXjB&&~^_JtR8qV`RJ&rbHybb4LFem=a_f1j*G4Sr@JL^VGmH5^bez``t)D=}*t*IV&QMn_oampqVHq0Xy4Vpyd(r{FLFQDVN$J9h*<6hjN_X<GdNTb>oq)7Ryo%BJKLH4k0&gkrvLdzA6UPZV*0&J=P&c4x`KLcxVK!DL3GK>KwBAZAe!)<f{a~9bdtqcxl7DRF8##oZ;s-+ZG0e_pMhkTUU<D3y08Q+y{WgU&}(C>$aHTMYe?~tQr#^k}#*j1Jx5=brC7njZQ`M=D}sF$B9eCG_6YJ1Se0N1XLSyE_chUl=9`Z6}|5`+=D{QzqD27mCgfyoy6^DEP~ZaDFW2b<cto)A5%DFZ#9&wzNbFd@N%5ztzLJB;{k=*{mTlZJN?xLv!BITT=90U<5*kPQnl36QkOK+}1PV!}<Y_~V4XUX4=bJ3PUnsLSyv&0xj?`o;O*NVRNl6wwGkY-Yf(EEU1|Gd&u=#jBlWr4{K7BUOZ+ZayjO#kDHIo?#z?%{y9ea$eJ0;<RZ%s#zmdy0^^5O6*q#EqA3#x=2utPGawpc&g9Z%7g`5eJLGdPX&x(bMFV8U?9<=I`ET)wPg7RLnnk{$v>Yb{lz-e!^BQaoRpaE?0vMHtIh{*4xic<<<=?cE|TAR`y*6qW9Xw2h9r#gGUngyaG0@!%fth+CIta#q<CcOUb2p{hVQ}O7V+HH=<BaMAo%XYJNVjrQP)c&f<LWd51Q|1kw(7@o(N=XKxB2es2G^jKKqO3iV;*^X~v2iZ)9Q^XHj7h|1<og3HL0kodIF3dRNfti#ww5n*Vs?nAaTgC=|>xdUi0nW8+oNgmFE((5Wam!^dG5XpI#;#awvR)UD{Lyb-%7X`LVC>S*=3XAUyQszK1otS4~L3=sWq5V_BYReK^yFiZ@MORtB%za4YBiM!cSGXNTW{#vz1NzH3VLv3}=W~7!|f2$UVmH+?%LLD^=A!5Z<00Gw{tZ)DTb7MmpvBYQl0ssI200dcD'
COARSE_FIXES='{Wp48S^xk9=GL@E0stWa8~^|S5YJf5;0=`l<6Qs%U<L5(;N@=V`m9&^ahQjww^3`R8ep?UU~RLdswKIcve>^Eo8~LT0<huwFk(d8-ljdcTsCD`FeU!3i><mjKn<kOUnzbJ<LFwlFV6Ij5MXF}piIl^eJ;7!mN;WQ4j5?Lu9(9enkp4rE%c<JWe#kdZA%|^#mFUS>0CoV*cdyw)Tijr88x&QaMoUw^WfB01$(qf$ruOQ#INr{?;=g&@o<9Xy1PY*qX|}DJl?zlhOTY2b;%8dSd>xXv<G=%z8~F@=!c0Syxr9bvY0fr%fQ~la=8V&>^#aezIe9uMB@-D17<!tt^gR6!dG@i*A;0K8&F}s-)bd0CV9DhXdT{Bc0%7K!X|pMGyNgw4bbfP9%8_-Gwb7x=JHxX4d}ol8R7&!Ly(}?3W|}n7d1BVDexGMW^Q@-@^0~2VfDlrZ#?g4?gIr;nT`e?|NrO=VklrKw!0Yc%|N2VSsVl=ejt#ZxC6~7Z7r`lNcvnmhKIr%DUokp6W{}h?dk!a5Xx`AN_Q4hqer5ZkyW|W=veCxKAzP70esaJTv?&GcARc)tgfqU$3|GPV-guowOZW352*To64thncqn%%UvuzmAOI=KF!I&bJSVMmI+=$&q!reJAOHXWAFk8BjSi=-00I94mKy*7DEDdTvBYQl0ssI200dcD'
FINE_FIX_TAIL='w*#%hQszkj+|_DsjoDx6A?8)-Vs9*$3oFP(aCQ1%K_j!1&wr}~dPH$5h?NL}Od5+um<e)LOs<Z`yPbtq8ICnym0|a|L8eDEf&G68KhVe@$d7z$Z>I9B`1qYp1oiW-Qd6?6CvOhh;sIIde0Zw&0T3r>)V12m>v}`Y9SCsY&wb?H1xv=`BlU;qOKn}@;coCm{XMUxzA(`CATCyRD_|^B)pQke;Ky6nVrzU8lerZh#o`-szQQ$axrl`nUw*1IL<-DyE{aNh_xpZA0b)aunrcMKZu|0Nso*V!BW<GWLt>GYzLzP5AS3_x=jmRjc};gW@Gx0@?~fNl;zE@$^3}LI1-BYh{c$$_r#<jn8BsIqe1n!=$0X~h<AkbYJnb<YEYB^t8DNIn2d=+Qyt9DLY{+N=>jYfAxUQb;(3+&QgJqVK8KqF+W6m4_nwe^lgh&8~zVX8+skFG1kmI;)d^3HgexG%w$tz1#<|6)9ry+foQk#8(stlp^$#Z%TbQ@Hs%@+6bWOV<F=$LcLG`z301TmtK={zkvI^ZlOrvR_X+RcmzL!a8qFz1p&o9piG_jszJmUyYl>s*^#kTaednSx(>^pk$PmbY1T`MM8d%a4rn^g`~4=@fU^vV7GGhCPB-5+-3rCt}_&fIaXzHLA{ANx#jd_=OJiv(5H+Tpn;D>KU)&qijgxu^k!AGi+!EeX(<L+yv#_O3JZ<_ATL3-iv3mi7pk%W?&|!<K$fVeW)#nVWsd!a6c|sv+zb3^+T*6T$X+EV9~f-Ovta7z~9PuVQvC<a+v=vfI?39Q^5MSdI^_G`}QV7e{AQFypLAJDn^`t=mz#;1QF}Km6klP!}Wl61sx42g?R2om))3qZ^b-$HH!ao#c8okOI{_^$Lxhju=Wn)Zsg-|h=?%KooK5c!JkD$O@9k-x}e><(;v3%Bs*vL`uz^AJzqsgS#MUdk>nbJTYV)itpD09*POW*bk6F-x4Nb=Wver3d1YtT<T_`Vr#e{U=1!yy6lm{elK21-Cato3koN7A`w7rA<Pt))%t8^ZdIQv*eH(0hbm_Zx&>6i(i#W&r<M@Isv4Ws33FU<&<-9+sB>P!zQNB&L{wi$*M|e}GxI91Tp#oDoiQ@Lo*)i5%zcVkzC$uZ16Dv}?@HiQ^Ok+xa#NMy~eO3we$7bS7jgXGLQy?UPXtgP5q+d;^c^Y`0F_wuBeod)D%59^3Q7u4EO*cGnZ&Hjs2)A?_fF=G%$OU@Me4CKg16tb*^4ps}syLQ;&>?iD?C&Beo7+lM=QtJdS&N%T`io$JL(<Ll*EpXQIt(H5gq2pkcj?;e{fWv&Ymj^m9~Z1CP$j+Oo=CoJI}1F%$^HeIBq9&8ehyx5|7rXs$vp-aOsS^je)*hZV!z{AU8YZ>$`zmVG`o(ea(rh>?#(ah90%xJ70a*YX(}7VJW>K{C3g%XWPe-FQ!x!A-luic;U*ly^1oGX!OcnHF9zNqodK7RP^g+K_Z21o(6+zeBJKJ_1XN2z48MkygLSk;q;%+IKdbIg<g--nl<C&m58EzrEKQloGqqigUb*wW<Bqd!b4PTU{8jr$1B7bE?R!7*8!JS(|9^PQEK7<P#gyjm8@K*{{9YMEm;PQbVq!_ypoB?Q{CD5x1B_P0xRL9Erpq)1w0f3lYui^4{vn0!ItjPWC9XptCexQFhwOJD_{Lbb;2ea3<z6#x0H3!5Q#bF@UhU?^7G3=!hJmvK>?f5<k<QXW!r9R+R<=*v@R-k+2y<u2Ia`@ojw?&43|H}E4D>kw3O2t<TV<eU9+zG$+Q-249@)%mtc4x%R7NC(v6kCF{mdX^hA>bX*R`-KER+Hg6Nv)r$uTW|%J64SqM1CxL|yIeg9t~Vuzp~16FbVT+wh#Ng8>YG7^qG7E`-OV&xa55^f>Zz6H2Q(P&-z~rXwt`Q&a<ZN3bWHKo_Z9O85Qi?LfSnH?99UNkT-T4cQyh-OS_qW2*doz|(~b-nK;*{S0un5_NN8K8n@sBxk`bfdfFzL-OV>MkH}_D#goS2{4Zq{yz@b1)_S_i-&WM##DUN<!A6@4@(@knbOUDzMZ@uuNJ*aB?!!9xVoelC`7H~EtLrRj!<@hi$vpcE8}VDX?k2@oYjtXfznIz${nBk3*ex=2bGKCjXyAhY?8GH-%&U~!12FDTE!2wF*kZg)Tq1Ik#aJgT}3?|btY&F9F&)DMT&}(VUQs0B->9S6T1H2U9wFIS#Zcqj^(sNs5QK{lbQ92(hJ*yU_B{SnubMPhJ77m(kHY9zRvN!5rA<>dR}))?SmKTTY!5^hzPmR4Opa#Eje$DPnbN*8@YG^q01_=*{7fhU%9sjl8tbyZDD8raIwXY_#NtJa<1^@25V`GhrBwN3l}Inj&c@9QP`KEj9PO)*aIeh0pzxHiVibl!HsN+{Bad0^gS#d%H@R`(X&|W9T!&wfz?`#M!0%0Vvj?s$VAH~i^pz%=;#tnV1fb*m8!oCoYoRkj>sen8(a8mRJ)sjWQE>R>!ze_SiMf8Qh>;%6*@3nWYEypB&G$hv(VPH6JlFLdY{|O$pCX>bA1kXX>73oFaHMoe;WumGRxRyEFot;&9mZr)+oV;f0#mlda`Hy&CHFN-l3E8!~>9_9X@ato;oLIM6k7ljy6p{tVWF$p<E6%E15B;>V;Pbc@$mcq@;~s`r<$a*T?&XYAKG@{ClS($d=Cg9ohcy8^zYcU)3j}M9QnV>w_2<;7a^FY7ie$6BlHS;bM+4Yij6?+82)WblIn%(;bGG@GF5vH7n$U#CB73)%1i)<gCE(Gc~|TDk+6=Tj23j={pi|Otq#>scL6|a;66~gFYkse<bgKHx;YjsdG_(&3Ea4{yCSGiF<RC8xV%;162B%SwEi!G)dvlxz-F*FfJ~5HOoggr}1zzLCv;mAqYHe8ERizU!sx)n2J=VZO)C~76K?Cbp*lfYP7yK*V1cqyw~R2Xbi+@9nw{&nuT6uW}vwSs)}F4c!AE#q}foG?Ka|n;6~5sNaLZV$lkwr#gvu*WQAU|fyb2$yGb#W<2!e_hAXr$wbn*R+(dW-r9qQ<Y(Dk6(X9xUPKI#aQr>a*=4Ysjr#vv(X2<P%LfW-a3TZAm5~9p)rG2N`6`{Q3mZil`><WZ}@RKID`yY^{J$K%pe`Kv&yPD^ng>li9l6BR$6`U8$mY&h}BFzT_W)xn`qv24BCOwgMHn&KTDZg++4PpDsLpH*4DL93d<B|aIjYE1G!t{=-$-o>IAIjhpHDFpg2!Rg<_FC<EsbUKlF>2{V*o#>cf9~5KG?G{3qA`acdjLVf&m+(Y5YAM2nUB{J)xW>jzk69wSnl~@aLzev6=XNSp;>meXU<=z5Xy0$KC@F4lh4^W{fe*U7S3Hh%>)T^65(qLz0#bncKzc4W76>EPs?^)n4j_WilA3=m2kkun^Ws*QvrU)>3T)AZr~@QPz+x4oQ9DJB_TmoVlH65OD?Y2$b`v|p33YLcgKI^<14jyBF>u39Nw<XK&G)uP132#Wwt_7d{xaAoc0u`G~<vS4#OW@ASB}-`@V}F*;50v<*YXoah6mLDc}a>Dz3IV{CZ4a<J~<^jx~>Cs}xjBC`q~vA5M%l^CA#<7j85bN2&{ga-S;kiB|-aZNr`T9TKA`MpKXzXszYc5C?#dIkE}euKEr#$|7aG6pzfQJ1RKm65FI*ksUuu@A^d}5*BZ)&+zEhFD6lmX;maWE*&l{&}SU}u>54%-s4&9`YuT(G1tJ$Z|TNNS`nav<1H+sSay(zQ!}qrqHkqwNDKDHU?F7T=8tQ-HTy0+G;pKR0$p^cF`Zf<#e)pC&jqm~V9GNACWY?7s-g8rAiuj%@ndvlHh6c!P|AGZ&h+*WArFE-@w&Tf{j06TXaByC=t6aPM43qx@5Sp<;5LM;C0)n36^_V7OQe*t9GFOp@st5@6izry)A;+gO2xO;=~tV!b-HDi+<QF<u!4HRrJM!rvz5sRC@M8c4@qfCbFbdgL%2L6V=-n;Y@Q0ql)z(tnlq@l%x%5`Vi%~5RHS*^7r!mkA}*^)cy&XMtgY&SAc|<n)#;l{8tc&hb|fnN!v%*lgpeL)MT@802ig0MC(C-7*o-t@Mwh@){%>|1T8+KO#z(Qvei@(+Qw2bgR-YE@;)t{ei}3%hi3*m?I7IeG=UWJ~&s@bV-Ae|RQRlv3f1<by5#9;=I7v=n0Dy#k?wLH5Cb3JpG1C#f{+QqZ0id_E+?kr{u)D`%2{GN5c;kqkg&0J*re~x2J}Lg<Z&H1}E$#w|msV<`=kpDKJX$N4@(*ln7z83Mx!BY`4d~Q_9aHik5fuH7RacK<mD_P3jS|kG0&&SFK>YdEUpt^GET0ra$J=y0`|g>}#l<#9_|*)0#!tN<9EhX&iXi?vRZZyY1B|u?QfU2%^K|;TolUcxjyPEIb9nEDX-VwHB&b{$z3G{nTRG(+6`9<o!N8Hx5~XH>M_zGlf?*~7a|gBG1B<K9Ox19C7LN{&Idjf)(8#h!<R5vNef+$7qG#Xx-M;-RCGfNSyc8L`2UZ;!1Ue?yu1RhvS}Nk-Dx5UDV2>P-l0p3fVR4H1ju`<KkKY>vTi&kA*kHRl=Dk`ge(%2w5PSHmia))gC?$OHTpfXZ2TUEpM0G#Bq)f22@j*6x%<h6#Z<aAWP~}|7cbIoDLcRSl6*iZoz)O>qN(Vc#Yh>eqOrQ6SLZl(t>3!87<}RE!Ly5X}`xE=>_k5tQ&s9&Ny^%ExJ5w^exS>f(n#sgaZDVBQnXdOrvW7z7eqIt%u!=Y+ivX}?f}fhgp*FwYG!!uh!FisZ@9?Z#_IEXyS%<kWcTC`U^&qMeO5|YLMIVHqk<a{l4*`(V4ytJxDMsU*#_mv$VM8-L{pHqYLpO&pZJ#Cq(*lRBDlXGkj_-CS-oBhLj9gE5kW&U~hK30B!cO)ampe?K19g1m{Du5ZlQ(m9!;gTvwi~z~6t>#BRX{N44`AWEHhki^q-gUE(w+|`*jeiXPm9EXn><8Y_BG-7_t67e=UZmuXLO~gC66V=@ZWh4t`e8qzp`f%hCy5@O(8uoVN)6t#_|DLCAngd%9qRe(wcUsZHY-F6CC`-0HH$hXeO#~i&~DRO;liRRu~4FS+H=C89HVka7}RSCf_%hbHd@HWbOo{4OlVZ>>=wMZZ<6?T9Q>iBP{*vA-5&I@P5bdGZED8GvCoxRl)0PM^}z!*lS8oi9zF;Xx)35?$KC@FUwwa<t;!9-+3j`vla26;si@RvX8{+b<n1y+U8hwHoE~&PU4P4IFebRj++IMx%6X{vX3HnFV3mLVm5^QM*K|Q|CLF5+Ug392#s<gc2pW-dE@QGCd_07^7$@ex80NwhcKd1c22!~x&KR$=$vk;mPgFZas{Q00RD}Y7rjNuDkvc)$$1I@y08p0Np4#Js4ZEf(;yt`^fwoikj*r_k0>2x!GafacVQ=-!jA5~;X194o;56Q^<H7Vc#;R9x57-mFFIH(MvP>mW9yLSn9x5`*(T<zKzdhqBhq=YcEUK@UH?-fdnQ+4>V^=fbMhf#Wv77Lnmj@2$1y*G2^0->=pRpW0Y%%9KAnV$+v%_-PK22CMLrhhsuY$Z^>&DVeNnd2L5h8=*ZXjq$RaR>8%<<t{qeu7{Bzw5PFSog6>Qeew;OH7cz$`y?<BARxgff@AWVNrfnj?vQ7r6ZQ%(Uf$3A0%ab^A3w7fPkus9(!uja*_jGTP*o*eeqrzN(j4CvvdeI%R=0jjrb0$*nLQg5xP@?pCadZ4BsMz-^#={`Y7xpi?U=Q73dA^Jk=MJ;q|+KX5+C$E7AwB7po%t0FN?T@?Vzt)p(vU-ZlmSqu3+Ub1Rm%<@pcN_WH=nPPyIRQWio)ow{t6jAipBv+$7w_vYuDQA7+tRaaay}1p`Kv`TU{U}H2h7_R^z;FE5YZ>v?!5y^YTzc=G|ZY&@1~!WKipq{Urs-y4>Sn(z+7+Ve#f~u#ccOS&k*QfW1gY&ML&-*r|xo;31s6^Vx(D!u&}5EHu)>+Lt#*w2hGj>qtnze)6HRS-Ntg|lK*`suL|1HuASG#D4(0Q;)?hiap&KHG14w4W8H(c{>nR?GBEAwX!8B3NDteIe%Rt~1k7qTe*SKlnxqShQB3;V+SS<XV6wVw6Imxiby3W>L3?5^?o+oHsLaqz15kP+KocJ1`DW33YV8VLWO0e^Hv0vZJd~ac2&lNUG2gRHZA^GT951F|ftO<MTinC3rO`BP3I*RVXQC)GJ3uRjX{rXcWC%HrQ3vaQ%`De6_4v5Q8p0IFqB9u;a3~Yz_0YW3hh$$Vg!oc}NdW~di4kxPv|_Iv_xnO%Lg@U(>LNFTnT!7vrt8I?0toG{a1n**=`n#u?l%J<i$Bk?XB`b~z|DN1Ym4n2EVp&0^7~>1_d<K3<^Hs_*AZ5*L~piDPfK>S&-36npTRfgES)#VgsAugca|ZeHewxVxfpKNf2!NBDw`Y~3#XOEzXhDAE-tAxuX!JXx#4+x_nrdqv0L$n|7-nQdZDavot%Pd>Mf%tr2b_p6UctiYlImi#cBv>%vmR)68tr4O3EambZlu)P=bZ$eqldKF6ii&vK&1n?&tq<A*yR)_71ud(;N0(Oy2_Qp5U07Qkfb*0MZnshN&mvT;KZebqY>?U&tRVJY+Y0<u|jNF34n1b-7FUPbHnc%lSV_B#I@=&^cIB$ofl~aSSIv^(i4h7AnvBUIYo!&NaAh-h8tdI}Op)!CxQ(aQ5#?rESxURx1m1C(eiG;Lhwni1iO1&4z<5Afsqm>3sS}!P+kFX-{*cDq734q7K2QD#*~&6ClzkS%B{dYE=pme#9s9%l`<WMw&%^)F+DdzHOV#7T$Fa!T#jXC{xauN+R$9HJIQA-RI@ROet7hAO>4AHRI?M^=>TOPsOirA%E}X0sbP0hQ=|VEQ>S)`K}Idu{~|F_cjlm$HF=-xTVal`hxNrnLvP*J#uf$Bk-~xaxY$t+hXIM+V+@LfuT_m>HCOrgjvsrv?Z2Wdvp1*uM)Py9pT?*uu{1RedbQC<Xif!R7?1`m48u=r^#F$GirPkEQwMU!r2?QeYx(Q|6@{+UP5E~Z*f?9s`>DzXVMU|2yRx_nynm@chSQR39}+v$uku%ZAf@Ih+=QkLR!4y<)(0ch*D?{CD#{TfG*`Bl#;Sgioh@Xpz$!thcVp4faYw9bQP|l#2AL7LNma;4r05uL5&f0`c}XuQtqzC0xK&=%C1n=d4=kDv>CsvZ66=`IbgAhl2NQff(5yX%gn4I-qYlD;dKuMvvv?rBb*(uzhoQ;08-j2xkSwRokg!cMM)G=X?guii_oXs#Y2h$C*y0P;b6(c?luVG1PXE1GMme_%dZn96z4wc60Hc9%j9|-lSCdg90tlA1f^;6#uv`_$~^Tkm8$kSwq)G9Y|G?Xx7JUuYLW~Q$6uJv>Tr;4%lJ~TOEc8t2#|2&f4o!5>73@IS|pp(<q+*FjqRLJ$Nzx63-cLXtsE<6c5Y+%IYO`~oTVWW`aMtf=6Zg6p7o<*YIo>oBP|i6{0ds^X^9_giiUxV-MQtLfSI8^c9ICmZ5GwVc;5NS=1alR^wlezM=>d@aZ6?y9cg9B(*~H%=6M*Ke7Q?p)U&krpF5)6+<}6XV-i3qX!Ty&#1;@Pbw6h2UcXcf)APoIUl~OEq*jWOk#XOVn2fq-ylg?>9U!)E=fmjJY4Rg_BoQ%c0y>EsoHyazPO4>oE|>rESlGcs=VS8waM$}6be2i``|SCMC<g^s-V`UHQk%%_4*<6eP@xi?2x@)a&wnjbe-m3edq7b+7cy?qtA=zgJmi-QClESVjxjq(FAVcF`P|w3ptlFoFIgaYKnEPdi9LC&UAtXW3gsgXz%mP$g*GT*92Hwo&w816YKQCUALPLyQ6I-YTPXMq6kTWDI8gN3#Ng^&aW=8GdVE6;j7oXBO%?WCf<D&7ULUD{kn*9IbZh8&gR)a#*`i@jt-5^u9&-LlN1V@qJa1Y@C8JL`)QtAnA=(a`-$szr4N-$&K7($mvpq>o?%`?0!&Uc&v8s2J4&~Kgpvgp0>$0%jIl-+~%y@%)umN@DS_cpvu@Z#@vG_e(Ja$tC$m};6g?L-;?PCD8MFG4W5#%f7YmrZoF|ptKx^X#EVL_5E!6!qlx(szefBxzvrrESPD513ZJjX(vz!sOy<&%i>O=8%}$&-hepoCJv)~_RaPlWIWSq|N2WPsvDfk<!@%mnpjn*e1SKoI$&QjBPz>}EVNIQl>`p+tg%*#Pz6H*ZhhY@d@yRh5L^5u;^Aoqc5b`c$Vz%dKMi9f080^!qe6eR8d1I;h6=tStN-1G94{AnoO==06K=tR4Er9x~eObZ=`M9VqOXo){>}MAeq!bxF80(fk+b2>-^2*-5yEKd%0~#CtJ(%hAHW!im456CSLPKH4p&gfHjF@GjMz*VGPxiKga%Kk-%IrtmRZ>m)J~v5Cg!s>>xaoAg}Sx$X-Q%O&i9zn3q(J~?U0STjcru_Vt6Wx`aZT_<XAehbAzHox>#QUXTL-U;k_lCkuICK(FW%}76!_=Q!&AmP$(Gf92kNz6#s?Ijvtc1gqr@?rDt{VoNSSyooOe1O3WSZ|Yse%G_lAHXe=lhC5QFVo{3r!>=#pxzOm%L11F>BMJ95HG{}^Yfq(3ZCR<dM<TJqz<^K8<?slHjkN-r)42B8;<l!r)&D1Ea-{JZt&kIBckHm-~^z~|Ep-|mnz-XC|Peb>2g4*MKVQtm%e(6vCfdm2C{RD2BAa{M-2BnxH{s;kn4KeLEO4H5V2#T%=8g;b@Q*lcDta_76%da!>m;mM{p-gXnX#$wvv`^1av+v_z@91ShqoEu#o)ucm01_smA5Eau<c3jh}YuqRt3Ls;gsd8UI*@<Lz23$pBD|3G)=$bjj+LL5Kiz7i3bz*aFaqvkdVmTdr}ywju|Mzhf8_##wRd7VzdH%4U10ZPhw|odAK54M<09=hx3nhN!}plOn}F4fISh`geOS^E=CJykG_!S+U7lWnEfA-|)zEJQ!e3OWViKv3l6s6FS`#V8+kxe}wv?p`RXSb7J)_`dlnw+xf*arWSxM28{JVpq>71NU;QT>vK~W$0bjGoz#iZllWl?o4Y$w9bJtQ1L~=F5;9sX0VKIX?dPYDki`07^LKu1tUYld#40x)x>%vz7zle{DkTm?eEXZY$(@r=lEi~V??C|%)kUD8rL^?#9Nm|A`j+h?R(C`?qUsdF4QIOGtumm2c&MWYNFA0#1ec7E!^tGsI<e?(NjVnr1&ar1pV=N2X91Tm(<`G9%!r)nVA3V9-Zqj+Da(iRRxL<q&PQG0W*iD4$MwSMJZskM$(XDzSA29I9kiiv0}xZua`UM@C!Gs2jN0HCc1M*_QPjAFQQB-Wu7LSqEn5#5IB=GXq3-jUC*@bMFXg}v*xL2JqT@x+4c^bUi3{@UrGXW3zMEC=H}s^*582(+2UdFNI}cjKN2U(N4xepXd9@?Oym{BM5VTowcUr$LC|~&$OTG-+m9lpY{@Ysmd4N3V6;&4SpW?27bHW=h#9<2~@G?{KQcx9*`Ss?J)P3kN+1>YJLL3&KoD#x$SU5f=gln}LXD$D?^Vl_8CM`<%{cy5(Qy?THZ%IAIDmi_N$|WRa(R$^y)%YbsNtwda30jhq>QpBZ=s)^6w=YzkHl??+b|6GY6-L{jhk<~5|9tsClq^mXT<mE%oIVt(^<1>r@vEd(-LE=8awz(l5ROH)T(`kv(#Bii4vGOFJsYjRc4lH+oGZ@oxuF84Zwr0_CGTkSTd`O$fg`QSp!agG8&9ye0H^hc6V>NlcDY4Y1peH+aQ{07+7jf6%WfA2WD-wvV5-$lC;hvc?a_AdFnNRr-9FnHyaunzDfjpo1KF*3&<$r$dd)Oajoh2?o$7~-z#Ns>;VvnkHYt$aI&KPtd!*t$wU-MDP*3syDK66bYn3sW%U$B}#>(ARuB4(iob^?YB5yk^B>`n_K-kt$4_t!sQ(Vd>Av+0QI@xGkK`j@u=bXL+gmZEgc3xvfx$xTin1v8#6&ard<ASsM<xKazPoRC!2DtD#tveD{al`evj0W7DfL}%Frr%okOao95wk52F9kSKhC)WpkeNoIT3LqHNG$2=6|F_{(N%(ZAhrP3?hB3+Qt}ae1qUhtT!uM|Q3$Y=CzemgZxs=TaFjk8+9(qtau{G=cEF1kdtI-#p&Px_DLwc`}+a=a^7)rBW?g)mYdlW@g$cHehGzE1%H;Dk}-0q*MNYc~zoSJN)5Q^p#a||N)%n3Ua<(S2x9^*JYV!Rl`%$2bt)^^8ty&<2s<8_p|rR+!)5^z35ysLS{W-IPw)&62Urw&-VHiq~~DTG6<%V6hhW`IZ-ZE>9Rsc12p+Pq)dpV)6Lqx8C;)C+eNfKWX(`B*!dY8^%;l59vDt3#lYSR-`FbBcfJub;VP!a7Zt{1<_k7;FMDu?uEt`L|(_F~hP)N~5;O;M&l!qknoq;djN7OxIfytW@TsV~ulx`o+czAvXg!6|EuRIwwzQHDjtmPQPy9frIE_U{qEFVEg}b;>FFsvLinM_GynA&5K#+uGyj6()A}1^C9?){WCEC0f;-JJY@C8i7@5AlDI1|lt^ornC{qlRok|Sdu)bkMXpqtQ>&+JpmReNTmpEbQ7)J?8PlIg7@cqm7%2@7CZKq$%ik6;^o)g*YqUV%r`8H6xpPwH!c-jB@tRUD5^<^q!HLIijP65<_FwL;D|3^&mNqVLuUOABY4PGJK$h8fsl`)r4^u2MM<b>?TwPJxwJXqZ-<#^;V)>n14R*L)?*8j*bnh3@eA`?qqH}!mo<9zJgz}ZAep-k6Lrw@PEl@tU@LG<{X7EbO<e#B={^s}HXhrngUGR`z#_sx{y^*+fm2jVHzeL1I&8=lzk*6vT+bAWBT#nrbsDNt8RsXvP6zF&PIF~lIZ=wv*5=1Y29bB_tJziL+qg>#AmhN0@kYF&$!v^~!a|Y-4YDYRqaoa&^X6FTXN#m@Zjt^%Cq=>cdVboxqHtSCo1Q?Jh4m6T}JvvFpM9aB~K6a-(Uhg(&WwXB5pJgY&|3w<m0)IlndIgxBOANd)ym_!%d+V|8|1W$$mMr|mD*9K@=p5Q&zx_*yN({4)$>90Y*it!jLo>V%qjakbz$++N5v2Q)OU&4#xjP>3C+4QS+UB`N(G`4fDbr?7os8aHRct&=odLc&OB|gAv%Cexmy}3!Bq#j<U0!ybDtM`77Cu`{(K@?2g|1SSc|ADbpIzpnjm@j9CPD%I7`JpO{l6Q^x@qbrtxGz_R8L5r!TFQTz;q)zgOefq(>2$%va+s2n&|?K1f?G5r3iL~9;D<Q{;wFDE&NnoM6giKY<H%kn)P8_A2XX(uK|(2B<$=S%_4-GdM4EPY`)5izGL|GC29X!+F&riH6;f5oC*j|r%m9hZ_POJF5qudgVG80O>>FA>533&lXU~{Tdv!uV^|8W$RtMCbYrv)m&5zNyWqt%#ugCu?F}%D6t=OM$~+w!^ZV%lOOq!tle%Uctv$#}9k`UrA_Z92qGM$^2y=Ae+}?6TGNq@0$5tti0iQ=dxlQQ!NE;|n)X}Iy>~XaK=@?rdwc}-)f|Hgj2UEe(PvMidunEf#EyL)hJilH>EuO3S6KIj3YJsTM2UY97`PvEZwCffWz=5l-MQ_#S-*p{oxZ&iEZ)!{Z7s6k~uTN&(IUiTs;3C5MYOTzuPX({32U<ORvY6R}JA!^x3u;|Dplvsd4JTFAaoc3jfKNxXQ#Jva9a`7vhSa>r$E~fz9FjC@zFE0a_WQMM4POGseg#EV(U&gm+}86*0|WQ8c&WZH$v*&!X_vtPPX70qwHCd7{&Ss3qyuy0Fuq)U2ji%82H}O;vfhhX8gzjsQMjuw!cD7a1vFm4QBKC8PZ!b(`H&14$I1D^4oP7nBC26+3>=(y#37y0QX3+t6p<MX?DL`eB{3AAsDy>o=eVlr4QjY>H3K?DxIM*OY6n)WHi3Gwb{;4UB8Kk|Frbc=jm}F{3v=-ev5X(gj6=P;;NTidl#O8*C}9J8C}cD^e<^=HKoWyA#hv(6z_-H(Ep__>wfUz`rWA^xYz?v5xzIxj-ldZ@pwieB=XjQGnGPz5Sh$Lju(NnU`HpOf16af?V5~@CebV7$@OTE$3e7zym#&=}Zgz*vhOxs8N>q^R!2|sJVn6xyMkagTIg6TUB=>4>WA+FoucN=VBCMPVVjh%1_P^%ZqeTulvs-f=sG@s{vhr-88ObDnddEX}P6wT%kNKNu@+okCOkF?u;3O{nzl9#Ozcjgfpc;08<w7s?>nq5pG(HevJO!cDA6*_hUVoNC;ncIxQCDhvf0o}>`dE988`<CeE0o>7;XE^+Qn_!vzb%ZzrFH?-Yx|>VPyYN*YI6@Mr#r>_uOKW8mxnIKCTG%6Om8g&$x75wQ3TQgCqiLx^~(<uLcC|A=ll6DtdVLvCydm@)qbAVfYM)CcpaCiA~>>&iG#J^PL7%lsL_qC{ZE!~wue!5M2~-)Whqjeugwa`z!XrW57+bh*P~4)oql7tm26A_-=u+)f^Z&@phHPY^S2k`N|rKxYC5Da6(pNSNY;?d);;l_b<7!Vtv-4J<S6ZM1lHQ-)`U{_w20IfZDZxu5!jz8oKLFJTUA+27)tCKOHjPtZGS3y5i`g;!85s3r_(zx);lO8ce!^J-Uz7A9Dc*si*{H%no^y2{7q}BfPO$#N{5aoW3~(1R+vd)r&}w-+L&Zep`-7?HV%$l8aN?&@Ris??`M=+hmLrf;=gsUe)eHquzrogrI%@(F1=zxI_`NZU0TA%n4fCozSPkrUVo|;_M#@L-D2`{GeEcJAF93#WUoyr&f}TBmFVOWhs&9Wb53%rds@HjLewS^Lv&A9l*^>tg#d;_ta)9oDRsVQQQQg^k{~J}CS_ZV0=i%-T_r>^n7YMohZ*CJ^;1}igVjC<r(yNw6PX7AI{=?4akCmCNpa&Wp2DGV2LkYWlHL;CmhM9qpB$8<N%KFIuOBhM_0UTyo0hs3fvxWvX%V_{KYq|=Nm-<2JSlNGV#oh|sGJAfB*lV&vij!Ko6Qx(KXY$IL&Vp4xC%p{$D}@>vJzYY&#*_P+VT@~t*V!Ly%7G!G>|Zfu4En}V+AKh&`&o!KFI^->3<f=u*+O%4Zctx*V0<u?_yKd9J~~HnwJpYFqH>g^BdY!|M!k3gMlgTZ=KRNS{+A@z+hOls*LAbR@qoby2-}2lDgExoW%kX>B!X?E-$4cmD!x$4_5PB&B8{pdl|ho!U=c9mz5y61LjMz`&Zt7%u?WjqAY*(lmF!H-9FR=j$~q8avJU-30L{VX#;#QC<gP1k)jNS)0uL^%#ruhs-R?wo*R#maj=(uT7TCR90u9#!&qG?<so`*opRv^WB09IXGq`R9Z1mJM>gtT0zo1AY7&=G2hn5-{=3&+!`65OcB&iPJoeFIHqa>F0AFDiZOq6u0PiNL@HyK^kb(0a7D_yJ9IsX@Jh1@)P*%7T>IC74;L6YqPxL#55>PS^5oye?ar0A*%Oh1!=N5j3_8B7fQ+l94A;D#gNiwsY^bFE`Bj$psDkPxvCEz(dP#J!xk>QCNbRpa+hDUmRTF5E=V$&~2h>Ae^y!8ln3r4>xrRvOlR&ZJR&M;D5r%ZoJ)aS<=AM7F6w-{aO{^as%=x~?MRS=*Nb2ZX<oo<=pa2?o!eiY)u20nJG;$x#vjl@KG{*0z}7deWGw@BF?17fFr5{XZzCx08lkYUclki!(JBF||)%AA;u+B{y|wNyc$S&n=f=7tm&X25DGKKozvOvCQRj}H3;iuY5x9o;e;d5!5iz7-iPCB%-jRNapZ9u18Lx$^9UIAG_W<!Jd=g-n#jjS#HJ5p|Mv)=)ICwn$!ccs8Fa00xX)_<3tJcN3gJH>1NRhhM1t&Gt+)@h2FoPMkX4j9XCMpm_UI3*GMTLVO}!8&(UeVq$4Sl-HQij$FLK&sMLbte8gINlx|~oi&{22Et@1t5X3ekrcHlSwdZ~V>ahwgL2GrB!?B?c3YCCdT<S)?bNp({3gk#k<j~CE5qG#0SO9Vu0-}|(S9dzvhSY<WuXJrRQ_%T>*!e&NPx%>^)h%^1cWg-{`0&v{jfV-+vFOCSY}9g3ONb_XPy&b1|x5WXxz3`>y0<n+bGx!_kC|wjL=dKypt6>CJEUCC_{hCaDOEo`&DUr|4`A_?2`#$J{;IXQT@nBE!YBmcJha=(a((L+d)`FQ&mfU-6y*4b>SVrw^X;4aD3GOj2RLwJB#la=^XAHrv2yox=<ns0swe;#@2gYU+F3xPJrjkz)Jl}fS(}gRKS%CCg<1mfJmG`;9+|jbaNIvvR@%6>Vz$)fWo#5mP5cfe%F3@qT$hN%UoHF;r$^!_0s?ICRdM{e5>ErukkcNk)k@YU|w1p<+j@p{NPCd6GL|ICW=!z=)AMHLA5FtHnSL@KWgI^lLY?e$Wd9Im7QA4rSwo!9j;i>L@qMh6&{8M+chfU!~w15Z&7?=3v9KRH3J~h5y4j9>O!+wsj9ZizdAkRiHp8B9grdj`7{Iu;5CF+*H8OlEI)6QWMj6*<O+r%QlP#Mc5EFyUbbP3orV@C3i)AS16twR#`I_T7G-CjUd%7*m}2Hxmm|p<#lgWPX2(opvawVC_v;Zu3cweUrfSOszmq9FtrC&0-rgQ3cYP$U4#uH~jV=a}oicgL&V*mfu0QHIfqmh~f}${4Kb{HML$V+)-8cT#!E6<<<zE(ee-1cP0L^MW0CFWPfmwII`K!ty%Z&s2m}TK?ir3MJoU14a6foSG$eGNa(z9LmxU|#7aP*xVACQXKUtppO2hKGwA^<Y>)9)M%-W{iC1U6_#3<8$i0kN<b$U(=7^Tso~Ld^1+m7!Fi*agLf4rP7?Q|Oil`0FnBF#eonU4@vQUWTv44V)Kv29X@<F555%eje6{&-L6U@c@?CZ?0^Zdpp=1{VHoZav8ktk+bt3eDkh)v_GP$di)z;Af-bYbnLUUd(}5NX*$3xE3pv{sM`QHpAZ`A&MHM9e&y(DZxE{5f1s0Sv{*=p%~w7}=Jv_(tCA)qOXeU?hLc6U6$=uVU<$c&2}mfDWB_N4J+>w^E%8=<52JZd{g?m)=o?}mO)&dA4mIDz%#`jP7{a2wwDL}U?%m$v=<}M*k80n~59GL!LtrvG<BSFsWDT9a|1Q#N5aFM_SsY`}2bWn&OtiN5j+vOO;kYpl&Q4Llman5vhFl)7JQU_CHNi)GPjuB^M*G4C%+A<?P4{U-bxYWUifwKtbTp6#xGXmC7};Fz1o@*^6yqRZRU&vOud$XPD)k$a4sN><;C#X@c{p?;Ws!>>v#reCnqWxMRVn%jZLT9Du}Eko3OQ(xY7%}Ovx(%Sa*16AxfPbkwzti<I8DXZd9zZc&9c^?_!SNIi|edrw>!6!`ASX{Oj?_WZ|47TC_u!BGI|@|{JZZ0yuaRWIGR>s;4*+9?$!kIxCU}vJ1ggU)VssH*UW2Gv;L%ltk2kkJDp9!b|1CuYFZs<(0_jjVlcHA>E<_|=!cE0#dLJCe3|85l~B*qeqhl;`OLd6@4x6?WOtDt2wOtbC){;?I>m%uCyS{LyE?@Kdw|i#C!B(IpLA^)X&VenVc?VZY_VqH!1w9IWitBf&}CCNF%(l6#xZNY%)O#@Woww8k@;_C;{C@pO{FKUBr|wUIG?2Spa;$Fpi_h+{Hw^mZzp9lVoBT6LT4@MQLDwrVS>ObTRs@Ak&^PFNdk0%$#9$@>*}h6aoVWl)KOxIUCBt;1817JQ14q`wS7MAcDtkf_M|p3n`Av9H}@SrCzGz1!P9z5@}|&!FO4-!R+lBSIhWFuy<MR`+P{-HMlXm;buI3?{m_35E30049OS?+2{C0sf<p2%FY^X4vX7>!$)V!ep<F<D#OUOJA!T4c?(_6)!i8uO`J+JL-i3yfjrg@<hxW9zvNJQAJu+^A1D`KfGTGw(C85Cj_~k<&O^m$T;0?|yg^~qFp|EsbCb_`gdfwz3h66)HrGhqEk8ns3D=j6g-=<`RBh(Qyq%l?{TeThU-Yh^DQBJGfjXt@Lsi0PJ#r~5cMGfiEa7GoV0%{oPaZ5xw_d+c}N7>JX4c&W*u%kLuP+kuC3equJIqp2C;_$%xbGe9SlW<@0)0dqG>V}|6Q7`Ed|Erpuk92h-L$#%-WE5+uFxQigy|!3yYYKNsV}sFrFe}}$1w6mZCn6OQ-DZBSV-&A+9Ijh;j};>U5^~D_tdT3L3Q_xEQG2BKZgVWwem(QbJ1aU}aBcjhn5n9xQm#m-I?jC^6gEM6o!v#*>>~ZF^OPLv{-oIgG;T;nri=0Tjp{xfe&j^&=BLXw@ocURMVo?_Vvqt;P1YkXQDkS)cqk@Zg{Z+RAF}tY?QyF4+nz!eLRzF@6ZVdJGNSgAY}u`v`d40K9tKEPft71wZcf{HvCXC0IfqL?Mrb+&yXcP0QaX$P0XiPtWjr6l^{g2;?DcKvec$_^d&;TFrclmx!dv#DXPTb}`iUz=X)r!BjI$0eL{$R1Y&qIdvESAN3Lt{aEHz1$^8iSy-Sm=?FT&7fPF4fY5CCYl!4<hgWFkBYE2IlR;4esNIwx>j#{Br!$K9q*>qv|$q|gdG+vJV>CchC1??TD)*SQYUwuE5Tp!B2H*J<aLGLNu^3|^_SoP>>*+_q;+(o{o##l%cP29FAZx_-r!$|&0~w4ma?sL(n*MSxO0XhH5D5Gf=gOazX?B$D2tCTC!Qe(X)i!2$G=(c_n#Bx<M2TjDp`*!tVAiygosndpy~g<8N=0qw7RA8sV9A0|QDqOz?sxdwXuQk2j52_ffpc3qE+{ixp+$pAm=cTaF_4Y!^@4v-YJ4+AUXn=Q#)S^e{06z^#*N1!k+R=zJH9bogAz-iHVF?kgp1jCqOX{H42K%drXiBqcp-`Y`CAd%u*if_ecQ1IeTE#wH>zr~3EXjAKcBthw|P%%3!ZG4fiY9docJjoZtM&x;;YZe;t>VVai!d1&S9hMLQYP&@u4|ZD#b!Z{NY#TIZGEf9VMSjQe09!;rr?LE_DCw6)t&Nn?yaJn)LO9s0wZ5JlYDt6WDy2^q(|Hm@OdL)ONq<mO)BL^rf|wc{H{@R06(+Yxn}rngfd^bT?k=0KyQ;J2f&FoR*%qiSXEL7;oty$jLVXYo03WLDuJ+m4fpiNZiW_CSEy0kyBj33@arm=L`6=u!LNPhQ%rzBcJ$@vF^nnhFwHZ(zRa;CUlaHG9z}b(4L<6jdrDJA4P03`144>M-1nKpZ5uZLc<gh7Yp1<G17~iIJJCp+Ipot|`4Nc@G5sh27E0;YUn1$Og&5Fx~VLo@>>?0Q}(}522NBD|0)VR<flSwAn!Z$3i;R6+kGcOm`LuX5XD>2sQcf?90rgjCf#)UECF6=>pf`l^`SlBPJ?LF^irO<)`=6kuSN!dLfwP4!cUC>&g)}F1psVkYa8e605-YS$hfyr3q+I`&pTiF2Ka*DTgf+JhZ`T#spB&L=lgySkq+y$9kWpAEM<F?mp1$c<60E>OZ+27T+5cq1;?~UWXDd9PN=A}wpXwSp8ghx-IupPIFRe>YY1X1>cL|G(ihBh9A`6*w9;KP*h9jZABhF*^g<N9lEZyQtkQ2C1O4st@qrgg}t^N6!*0Eu|}vEd27&w<cF&}We@{%_fj8}ctmV+CHFa{b#FSBq}<!!A;vJ7kj6$ZH|fY}Fk?(`&L}kPd6AtuV>fknMj2k^4`+;MC;2jki2GFSZe}NeQBpgh;5-4E{`iw3#-m;o1XU**z$^;yx;xZmmL^>8Ab#d_2EH6Mf%hBZZS2A15@4B#(g28Gv|A2SE&gdaVa%_3&t`?u$H=2=>GYP;6QJF`-*7^Pm3qbw+h=X{N>MOSpmVUFQ?(&~^i=kRv3^95aY#iY}rOiRc+Pa&CqvE^hXpM2{Ys9}Gq!=ke36yh$qsVz{QHAfULibyP!g78I9n?<}PS<#Bedc(sq|bkkeGdhA7AH|n{6`7!Ip(#-uU*ADc^m&n2d#6tDE5>3vkD7@@uB_#AQwuPO^H$zMmRB1l{VJ;iDshR&EPkNj%U0h&CY9d8Sfk#s@zH!MvI(e0-cTk3<F37K(u*thwn)-xjYvdG6n7iPjt&AY=<V!N;{;rDpcdfBY?@|XPD`vl7|JMyO2{defPz)949w|Pf6I7l)UptJ3pV$+#KiX(^u;t^wv@NZZYpMQB`@a8Q!%JVC!_9&}#gxf0s_u;CIuJsXn!9UBDtIUYsfoZ5y7!obNpnd&*V4csLvtNvz<6a=`E`C5MO5EpA?eHy1z2s0nSqT9Ws;3sbG_S6b*d#Cq){~5e?S;62L*Q|ZEGL(Q|ijMupb6YvSNjOY<CM-FmvM@TPfGPX*W6CSdv4~go3VOpu>%y5+$st-yEZ;-M3ZonGCdO0-70AKd#6cC%r5~?z{%L{E@_c3$NLqkose&``|sUlBk?7VP_XF+F5OdBFP&?(@?KktbO~{JOfrxAoAhAEZAbHcSq+CI1Jx+D)Tyh8HsAXBp7rvWAAkCc>lEZU^MCy+0GeSnaMGEs%1pE73PGWcCBn%=|1z1ru+&x;XF4bpT*uc94G@$bKRcHSJT0}7+D2cvJ9s|Es>`ch;;S25+G8bSi=9dFwWH2@f`4O^E9#d6<DJZ0VbNHz*cKTrI5Cu5^GC;Z86qoW_rEqlc+r5)~Bwn4|MAE#g8(7O6^qVkztk*_JwVdzKQy}Yo4>(+h-B}$IbKhKW3fX4`7n=$AT?;GIVBK?ZV*nX&JKMlnKPt`>JcVc+d39KX2raZv+x;!=I+7r+2oGG!Mn}9q1Dp&=M%%Ighe&r>8kt<6N%NKD91h2$#8u0y&aFN>$wEQ%Ml0U;_GNP(<Mg{$92@Rf=o$?)!?L4Ds#FaH~}f@_o6#@1u3vEirg?qAb=D!Q!}fD6K7Cy{pPCUR>O|XOjK8$671o#bSf#(YU{y<3Ppw{M8G!>6rhjASWfOE|&K>bXco73whZvJ9I^;yjv1syF^42!ShN%*!`^OGTH_^2cjI=kBEi|!M(a&MrLk>9~yYO?7j-!bR*mQIVwFYe>tq)q(o|fJ`n{t9Uy&6O~+rG&&(!3#>Csd5|lJ+xx73x|0IM)LRuD}p5n60$otmbCSO3sCrF|u#ci(D`3&RDDJ8S)U|_%x(z+T!rP9g6#`9IP@^PJTi(Fj3Bm-j2IP(lBI$O*6(Vh&=LbD*T3=hldYs1e?RPpF4K}LPMKAA#Z1+}B@by^r8rC^m6O?wCv($m24mV-H6AUcB~CMT|_%88hTgXpUkWtNCZ=zW~Jj`HfSFMYDEQhIYf-mu3pt8^ejAX{j%bjQ$>D#*J&57R#j!fJf1-^;MzPi#PU-cD5Jg3a*5;J=8BnnMJF(erw0Q7o@{nD7fs9*f)@_1(WbuZoG3OzW<lFH5I@iDd#hA24s?8UNXfeVnEyu>PgfawHSlK~p?&?1HmKxzyJj=E|5#VwPXdaVCLl_e~~Tc(*(5=wkcPCt10u5L#?SS1&9yhOaUe5Ez*1xxZa;2QU9SfOOrWDxuJ097pm9UOxg=9|Z2d`?s%P)cKg0xxMaFXunx}!}{f!?ROLPBR?e|rm^h9U9jugIYr~!0mwI?bXb9L4sgn|Gh`3D={dqe445Cw7a)64vkK7#`0VF}s&+EdsQvYbSZlKVY!_Js02KF#{>nDZ6Od6k&EJqBVPI|P08$+%-YsZ2)l0YQ(hdl<#KORqsn?2G%6ePHNE`PjDV33Tnc7Nak6;m=mf8BT$3O-E3LE~IpR|4vW(|H%37=5VT(7^qdnZ;;nAbc#dsa3UoUlzeU+1TlL5yqZ!9l1_Ht~psvy>m(o*k4AZadJu;VP*6MD}QHTiB3`+N{|c?VSt|hU6vGr)&#w?u_k$Yf1!l6K2Q9s{}^RZY3sq0$A2^dUJh72^>)E@HPIwcqcSYfof<5Hhgj%pB7#7H}+a~EF+lhmSkNI$#wD%h4%qJ>A(T-K`fIT*rX4*jvUwDyD-Svqx+~VW)w}z*Q`ml6g_XGV(fs1>Ja@vf6~z3qq+j&h=V%yM7Knw*g0vu?iWm{rwvC`oD@9HAF%g<UEnCt)w6(CkUNC#_Cdxik;JauKEv<?LjU{`tavARJ5h%^+EE4;;`eO^(&@)xt{x2H2f9gu?3Pu3T*6mOX3@aVi;TymrK#Efa|`+rTczZee|Uog)|=|xD^~A{u_8XqB>VL{JwuTk6hze+cC_h6Rj7c_xz9KQp;~hDmg#Couc083$-%MrDSQ3%&e3*?;n*fRblT7FxR=Ell-vFm9!d)MPb97QLDPc`1%9ID-bCMd_GvD>n9jBmyd(y`UnP<3u?B4n+<r&1=l(!gW|MTfW%|AfkTmdLVJ{vLoWnk*ge<u4`Bcj$Su%W<m^<7t^>;<;c8hijODjY`oIsa+{sigDyO30)MuA4_<x9*ps>C{*{~uw1e8A-(h$%}3j+k2&hJ91SuH6UU+Aq4AzF^0a>(P9tQ^;yf-S_AOg293@c4gSIWqi~PhNl!Os`Nrbz~b@S9^e+eD(76L?q`1Tp=?)yA2*jlxwbjXVLb^Zw}z2NIf`XKq=SVNxNwolT^)`5M(Xgf6i}3RyggG;*)qK{AepfaU}*aUO7gWA;=w=4$<i9O?1K#WB@O?3XB4ll%egz-&q(NuWms!H@Xkg!VvexpX!kWR6TY`eMi~aWYlHiRajMF$m?%6y6q}QWxb&>w+ZbbIvW|)c%N|@cW9!`#G2J_)A$6Jx@b5U#{+C{z@L}<Po<~a`5q>6vBnliIGy!xfE*ZBeBFt36FYzZI*0^oI3F_ln1y5UOJMLC;BKo6G>5D@lxAGEvS(=%>oaT9^N+$GvMnyc2St0$hwogoI^V&KzuvpexnBsA{t=_j-B4Spdqn!$1668gO_xLxu{5}`p64C59xR<^lp&O@z)80pBni8mnkH>L*K`;k;=jG>-G)1F|^EnP+8GOn^BxGd=lwJ{}%ivTgwMLxSQ&p4^2uJ8oYBQaaB`7bESh!a``g(IaNy})cojx%_ES{DY^u?07^=;^<zDzwK)$e}cLkt1a31~Ybw<9S*Cm)*tYX(TlE!t)Fd0!zH%r~NGao8cC-*^}|8J>vS6PMT*8nt#9)hN5Eks)ex9rS5U4RT4T$%p&65I9n^8q_k5-qJgYrfS}hMHlzJ0z4SEbTXy$k-?2A-_LVCGRVLR(!FHIwhX_fc#9Zo&5DAp%?&5<91sAi5;a;R-82pj&*(?!k0TDDYWA&@Q96Aeq1M+YCzO}Db#kJ7Q<baO^Ket3maq-}-)cL6zkmRiaDqb(`4kYjA~TZ*(^}9&r<ZSuX&<a}jxY83v9~;ws{mB6NfH=hqde<91U*S_JbDq<!fO){041P>d6u4iXeY9*FVVHE#c~!O)y{}vHF))a<ts%KCsA5hE2{>IzB$;Kn^ua8zE%axizCH@abPLRtthdD5;55v>Sk0&vcO{!Z+lBEce|}ZEP7q{yVP8=C&M*Sh8-yT+_ZH9%cXgfYs(hB_5QN?_|-vA9eAz;6{JBvp1QL8W%IM7@}#yr&fd%~P{`=T?2JtqfOoTv#rmFsMXR)4545E7Le>fty>pI&q;Ga&B@Q!Uqk(P?onb41LI9kP^8R2g)pRjC3Y8}KoXXN;v7}yYzSpb6t;4nsjC?ZTS9E*UIHr5GDB>*8@C^>vk1m4sTm(lO#Y##Gh*`Pv)1a4%%r+#jDnu=D4tE}JzLp(|=Ha?0ac`Ej2(YHCJK(~%93K`d6W~49)HaC*Oh)pk?Y=c(l7R+w^Wq;n=m|AF!gsj~q#Qt8stw4EKL>fLV*Z6INWXy0;v(n(Wg5gc0d%fUnGi<iskSBkjqp5#NyxnW+14g~t~L!sKa*7!v~2Aa)>gerh5s9T_blV4&M5J-S8Zi`rIo^0W|HO`W%SH0u}^=*A+_QHj1tzFFa5=84$$$ME>IA7vb(>0J;RO<ArfphS~Zatpk-i>?~=CoRzE`plUin_;~S75l9!a}ryR|M-%Rm3jA0DNNpO`WJ?!1F;?L)V_-uJxy>m?s4Rf|<jIfW<h(rXp<l^wGmQ~|5ok59nm1*`xcqun@>i*$QohJQ}9PQ^XFR!Wxcn>b$ThA3Ov+wUdfOYuW6OXtjlJCfe`dh{BO5J^ud2GHo;Qi#-=!yw%pm5j@e&k>g3;lh}D1oEthY;nMhvX`j4!3hpN(kfP{gss}`1rvQ2zLeHgjN_cZqW@UeR8VuVqbpl`v-7z8h-d8n4i#Zf4S_E6;&FyTi2!d`iw0=c;`Y<SV+78dx#Uasie;wzI$la84-6K7d|&{bobOMpF|91xY~>c<bponuQ6yw4Y)5pH)NQliT+-S%XaeTbN#X}UmV^)xsTZke38EeM+*nIeVv1kV;2kd!Bn_$9(S)Ri#6*zfPr?8mo=8Bcqb34^cWjyx6e6v2)GK}A5p!Id^XsaioQ`1e7vmggx`c}5;Z0jmCPqT#uGr3N+3M7>;vzTuBmf+F&nB--4G1N*wjJQRUSb~>mznKT)L?DXr|Gt!9TBk`$L*4LCrh`6}<3wN>tw8=JH_|91-Kp7s{*Yko5igl&;h+k$qj@0!*51`8hU@#;lQxhDMz~AmJSgCSYq?Qhw{127+X6lZ-XUKoqVo0ZxHElGhJP%qn^ZhzE>MMwwb$t5qUujIqjE?)?8k4*e-Zy+`5A>M-=b?=9MxM*pD7JsqT!GeSVrFXRIWree>riRKq(w&r<Ldq#%IR0WZjwSfV8l_m)qFxol92>u9pEs6T%j=j9U2t?jN@>qUk1kbh+7%9`RZBH(?c)(;D-^kyv{c`K<1C}yy-vRsui=s`3#}y>2&k}p&7$Y8=<uDawA6$Tieg$L5egi^e=)G0roBBn8&9!`{si{l)aPad!9j&7(;fe9ZFm6^=CwRcVT(%});N`v$E;So><}flaMRArXV{;zaJpwQD*fEjF?bl7Y+VT%}YGTF0r|u_OA+oLbpOcTdn`poRSj?M`3w$1WZSg7WhK!}iIj*s6x{5VxpG5K5`|CD+BQs4(w@crZhnm!Y0rNd*#j{{R5PogwQ5kIfWoZC+n*=Y9Q^s?blPma#_A3%VEp4pSH97-qj*WJsD8b!gr`(IMe~fj9J&=hD_T=D^b5z$xT!TgQ1i^Wquzz@mMmg~H@>xIwrRii0s<bK4?pU3ih8_cznC}7=ib}|V5ksvY62$RKiFScCrZ6#$l9b}UH1_k<Rgt0F!HB0C_&jA8EN#vH-a$mb2IdBX+co8>(5>2aUEA3rY)mdwyO+KBP<1O#hu^@3!@6m6(00G`vQU$hS4;CjQq6S@FcBix@aLpu!eKMF#2m^m8h)nVmK&^nIT^kSHgN!0wl%a5C6YYBz!K$xyn$^mo{;RmsdLSUlM?(A(ZzO61R$|ODh9rd6nrOhivKg<j<qtyzw$|<scy~C{@@i`P@)jfPQ|7S-kZddOpeDSG}Aif1xF+6`J8&;9qXZ@kjm+ELPlb4fnA4QiSYb}8zhLnf=Y~rn1H8?L)rfqT(Q?tpTeUHd6HACsBs=ktBHLoV-hlyP(01$`ZUQoV7&B-XMxu;gvU*PZ?H^2Mb-&o-x80>ZQLezZfW+hbgC8H^9^4IvY|NT-uF-e5>oKY2GPEAQz+CXV!I_?0ha+|OVE%I^r611n|vt7+5p4C!XKJ`JMxSO<vCjd9IjqAO(a@F1Qm9K<1aj#6SbtRlwJ>dz&{ZCX~_$OW}|~|H9D*m<vnt-%LfSmwUmN6ckqfzxWN5{&+>occw7B=Dy*=P=fDIRrQ|0?*UcCbX$X!lM6gP;d5nyWP5InUS^Q<HH6ciHXryEbboUVR>@zcPA!F_w-n6mvHHqQLFKGmm6Dpt*1uMzPlpn<`K=J5&={x>*)rgy`RJk{cucHjQaq7LTx@LwfgGLO9w4m?mq5nj8bL_7$x$XnJn{tO&K@^^;Y9bU9E11uB`eU?}PVvl%+LLoh^V+L~n6+il=@6(}<?qcIC+$|Su9XBxAo0fYIRDfGYI6g`6ur980_-_AnBsU!+23l%f7*i!rHx%-RW*bSgoigcq8`1o$<~<LH}e5G?vh-jiA3QLxvXlY)+LvRbj*ZN00hQN8(!dB-tUWQhE-{)BbB?Ray7t?5bg6O^UC{T6+K&=b}vu!l+>)Q%YD>_k>`pEzjkq<rmJVP&^gx{;<7sYRE3h>4}<sC-_aD#L0U>rY$5KobqH@qo>*%40Lf^~a@iWt(YG?nFG^JrrrQL!#|FhXp!vL?cWwGyk2qT1k$-2*dKeezM$a-3md-a9(9Y$&XfP#T0D?p^m(O$KTO%C5EPlbP?6|3EzX8=GVp=kePcrG(<y;i2R-%Tn5>0{X5eUFA;0WeH=PH;-X0l2BNQMlX+p%jrq;t|ocRq4?Kj;+4SMQIfO8sZqkZIV+NQB&ARE;fm0#Fss;)wY~OcWWJ<*)a0Le*vfW0aIB#&S~dAc;3615phMJI3#!diW3}IOzY_BZJP8lKRIam(<<oP;c)!4ObVxT?)wylgWRcTa6ioWGF0>XVKZX40C0HExy{-<a6FaHzKx+o|z@5k(DKWw+=JU*OWpR`V*_=N&Eb`wa9Oa+}p<Az+})5<6GRccGZ|X>&JjKUNu1`z+M{IcKbF%6U(CdDNQ<bGd~G%Ok~?nd2YK6?bHT-pJA*mdHX6O30w!*B?vrnB5rkU@S`oScy0qeh41AR)84ulx2{*2h_fw4lj>0doc44wQ4eDFRXskks!oMEc7mz>TNM9vwfxDfaT|hN2t3@cV!Y-r`Z?!IW5!a;B??+7P}RBG0{Ps|)Fjxd)~ATa1O<-b(4F(=<kL!R;sOiVY@1u2=U8HcS4-}%1{;RF$6<5!d`^bg>yyAi!IEF3v*gOv*QNh5ne=lSsx=8?f7ZSM8|91<bPh6wBIjAuu<wT|>^Lg9Cz>z_Ka5kXy+QC?>Fy-ZK9zYiPgTs+xP^LU4L|m|sdmse4*2ZIod6?-^cWZ$A+O+~EY902{z#pA<#7iy=0l5qH4~TxW?YF<R9~b8FUFuKvo4>BGD;{@krWBT;tBmldt-Gli+=<&ca;u9%j4LS7z1DQQvtaxHk9F|OorYiyuLn(-wt@(#lPW&N{(`KVk7VkUUR2_-@;doc~k2GUBG<E001kj`z;e21X<3FKuq&9eh;a9;4!Z)p!eu*<PJ&hgz~(P`_{qzu2gU)m3Gz9IlM503dI{&-Nt)_m_4tP^xfH9f2vM;$9s-gf?V7pOPRL2`nJLk54>dwIS>Y`gWw}Ve<3n*YO))<%PbxVbsB^YZXY-yvjb|fh0G)^;B*~WuD=dm*4gviyaVFB`5rVf=Af3GBi9haF@2D6X8U`QC-`OGk^w-I2tw&1z%SzbAd^2!ynxYa{_uDlMM&m8kqx4p?eW*?pQ1jZ#4eJSr^=i&)J}sPZC3-GjwWCUo1|rjc2Y5JQh;}ui|gD-E3ljMTxQW4dP#N?1F<Lba+ju|0ep2XWsOR-!}OGtJh^7~KBlifO&#lN_BZ^|=dk=&toSY4W+`qU*G&l4K=qtD(qr{QNfSsh3tH^%oJlR&+Hwh&R3<5+Dn9-T|BN-_E80Z7ZnBcAH9YuUm$vs))EvE3Juk6@b`5)bOnc(ciaepRH@gO$iHRf`iwUc$$8qX~SVIXoV8c2JpWmETgA&KT#*5c|S|V}`+bnx^PTD`q>Z`gt)^J^Z<?D4To)3101VPX4^72T?_iaPK>rC#Pbv$EwSI|p35$FW>GCo^}2;c4F5#RXY(YTjQaOUJH*I1qbg-b&7DgAVx9N=FUy;ZtA={uBmx!aunu`Rrd#Q1OY-;YJqU0jjgjh2W;XWJX1S^4bT&-N%c3%g52G!Ys&7I}@YA6ZAzgS;}`5Lklq1}P_G+<kbmvfevYbjkvDP~OxDce6Y6B@olmiR|?=b<Q6Fio9@$Qb-^wx8C|}uY{z2ysFqE*=<z8{r}Mq$GW`gsR!U>Rr`L;D^g;07l)<Ggsr5g>Zy^fN@l3_d+Dn3fEzcIX{jt>teaFRX|NPti~;olfV-kJDNVJOvLqS*hG)!=fJdI{tu}l6O(9EzG}h3?K%dfdphR<ymPGRTWwV#mW8Ve?MN1J+_5N0ChGFf={o0_M=IohQei}UuJUMmEiq3DL3x@XIU7;JyLOBM_uj;8>EEV{Xl2A^=$LvBsb!20kUeuhYQ_Be85IW^X--@u%?X5tY&P33W>#d~!C&M=wK(pR4WPm^S5I`l#+<P=YP<zgv3{!EEM<>}y<uyQt-ue$utpK>yS94d36{xWUyg(K87NY~@Hsm55Mn2WnM2+)L*1LX*Ls~dEz&OyGSlJk~_r}oS*wbbN?Uv>aMFM+>?d=z=;TkBzbAve?Y<Hlewkrq6_ikH5Kj&p70{CU-MYOHxHxG451zXaxsiC^`P@I91-pUd>(wf7A)qHfx-S&aG|LaKM-;C?K)Sf~-Pf|8*g40=Kn=s-$rgvmkJBS^ppp;ck6wuf`!1~JbZc$0{YnBTs1F9-Nj@+4(H4S^54&NxmLkUZL1QD&|UvIyYuf21;l<o@E+kD&-;~p_(#!-d~%`@N~XGV<q3iUKdcgR^#I!sw-I503&SbxvmJ0P204(`}F6^ObqZ<Hw;dCqQOkZ0bQaE}|yG5Z#cfC2baf9LQd36Prf-P}d>;6)Tb70Ks_1dipXhgEPJepbUo;ju4(w<IwNjY7wO?7GlOgsTNaK}F#=&PSpN5hp(buECn%KCY+xn#SQE)cnFt@n#k*xu$wK#%697oNfv|1Jcp|gkDOWAATzVe$#59mKPnPBG{{w2ULSWAgck(u8Tdh^9L<5WsRyTo-4Pn(jr>v!^D#tp@;0#iFMVHJt(gi_Gt%jB481QGE!VKI77U+@*+g%waa<siGb1%8eJ9KH7}#yyh_<erVglU1*lNO8G!~i`aZ%<A91~&c&px8)Dd1wOe}+mLZth)s2WMm<M*^GR*oXpn!qh|rG3VsXfheV4?MLI>)*aTf|bCl(n}96U#P?;TY!5C9-nVADa9!DF%&vtin5n-)94plx)oP`7HZNigxz5gDcn2UbO)=hLZ-Lbo+IBK2i&`5+o79EenQHxExt#oPfVOx-#q0*|1?7=S6{r_p;v9)smrLrJ6#bbP8VHDsy*h|TAxNxu`EhozuXs1BG)LeaZ6YFtaonMPn2h!Bcz<>)4%V>P1fXFdUz-kr^6`j^-i1&pA~czeM;t|)ZxsWU?!qTEeV8!`4sYJk|C!-<f14b(Uv9DB#~nnx+P^Yonv*MKo|Vz7*^Xe`H|%{5ZJCRu}aq|{hme?WGu$*{#Q3CM-|s!6Nk<hK8nOs79G5!u3}&A;Tf+-aHMvvwV{RVF3Ob2rThjq=a4<cmF^|Lr?K7tr2+3{3|2tyau?wok`qX6K2-ZL-_k50gZLkPHGlCYFkL<j$1EY;;H}|vCSD`%zD>usYE~~IQ<P~G@0(^&>jf)~U<#hYg{C)dlcYzT=0?hpf#_)=^o9`ePH6*%>br%VRQNymf1~9lLT1yE(?(J?-Sc!FP^!#mMtE71VIaV+hw6wL3Tn)-wdAb@(tPH(P*QHkAY9oa8%N0Gy841-M|fqND-`1DL#Dek<~_&@gjxyNl$dbGo+MO%lSITBDm1C@%8)XDln)<BJLmEcJlmE|_?fJGyAZaznL(}0LSf^6o3wvB@&Bk6RPBRErACtyC-eP_2p9J>#G^W3(qi(tTM|eCzal^jl-GMOZW3wX0UNM!fdPoGH~W0EjVmEG?nd#T7ZW8td9x}r7?5j2QH%-8NpCIferZn!iJyQ4YH((&`a!T8P=r;T>`<-+IxOCsZukWh^|a3pcfy!Q-}#<tRzIjxT<M5EFDn*Vwb5IONVM;c7-aU(X+@o!-#U+8PZU%O{NPn>C~0A8PAz}Ov_<{kXE#p9>+@w`6G4A4HLsS!Q$!oCYANOsd?vq`X&hCnsyXiG_yokMN`M=mT*<uJEXk;OM)nE>5GAZQT{pF{K?#J=31VOJA$FgwWO;>N`Kxg4%XumYqbJ_z%4`5fj>)63PX+p3U_figRIwdyu_}Um@*X#W9mr|QC^V7efXz*eyH2ME)DRL-TLzzXAj(-s8xA0Qz#l|EBZU3a8a(#3$0h!f>x$~~&EIDkXz@{MwOxTT46HH$0DJ7HH@<oYC`5<wXaKjRt0Zh>C37xs9|-&(lCuGvY^ShYqiakd5gzpvU9v*t)Kd=oEuEs@xWEl`5o$Ko=yk_skD87u^(pkOwDyTtI8eRf(Gh7cc#YxqCn-!1<{Tg<oP)9=y|D7QJ_UCjWXS45Z)#BfVGb~oH~~)$F0OmJ@ciB#2@|r~wNv*;56k&WrFC>3%hvo$LE%!WWdM`*I<8OwYXwz`;GTf!<@pX^7E?KAxC6{ZMZWw>a?4en-*~0$!TLwAV=G7VQubHi8NJ)0sl_WB6b3O~vfy%Y&5*m4u|wpH0JD@IU?aKP#R?_&+T61gq}k%&76+jf7^gFQaQ7Gl8J)w+KmsR!H_3oNA79g!;@8+6t+T3_t6AjvH#bM+*b0ude(?A@V$J>$qpf+T4UqHZ+H4PDZk2|EyMLk;#WEt0N{=DRTf@D!0%Go%biGxwZ`i)bM@2#`*Msb?%M8otn-G%C@rgD%KQkdiE_HnyyVhDkes+u}QLuNuW!NOir)i&o$JwMu(tqIZ+p_f6zlQ+X-*8g0@Mge4pfs~0g&?r+9<sgwAM!(>;c<l5G3;PWUAv~qfR<j>7tBryO4}8s60ihPtpS0?Srze3G^V#hqmi8~=u^B~&M^SX@&>L!W{*s{ih%y~+gWWeFUTihI934&P7CpaeB8K5N*941K8f1yOW_NyOgTTK8ZfAVe!kpy)L8sCy?k;kf$kT6o1?&FpV=N0mogu@vl<${vLILlX}O9xyY<Ls%F2CIQZKyyd4X5`b=1BU6YiJEDYpd&FxW!`g0MF%k;LJp?(b@{J&6x$%LLEUCvk>P$&^uUw|#P4HnoNo_NF`~@2Pd-?pOnQ<XjK$KasMgIX3TQ*tT6Kjsc2_`cIFHt3X<ROYwr{hRG^~uVD3KM<qJBQyQ}fn!~o5{NdAE4$NckRQyq3Ru+NJ;iU&1Me5}<qh-?;r%}>#r?m4J_d4&s905CHU5c75uSB_#*-~#_eNYPuz7(MrHNr%ne74t~efZkN66v#*Z1#{6uL1xOOdAF=R2w15HkMWN=HtZd+Es(L9xpw$EbPKZrE>RG(0G4w^ui$Aj*4bz!50WltLsGUeU3UvjFs%B+E{q1(91f{8JF{x3)o=j^@+C}GLyWE(1}VssrP6!n-m$}e}nl^nP}YlftYH-+^}ges|#+>&z6f*Crn9v1Gx#D{x`EY6l4`ao;|=90cO)ko*?sQ2ZN35exeTE8-_guZ|@Z}oS6{=kcK&<pst*l)pM~Ef(O4xj%N@UduIr<TODLE1u4&q0UxQ=k-rWp2Xt&!KlSt=;jtnhNkS;PI<appuVAhzSzxj<eOUQ!u>D7c&>lzS8DN37u4~foy3f--l-<&!10-{1e0_0Q({=rJHPH@y8KMd=)fFxaq0!=+==%iR2|Pnv@qQg6b357q+A<<ZeBzW=!q@jXBY_DTBl{F+)d4B+AnS6|>d%Sn24-8BPJgP>zmVXksIgL}^W~{-uMiL%jh>RPu^LEP|6;iZ4b23!0-p8Vz|lh;u~{YEHqECL1|_yNM7uiz>S#4jVKN904=>fWU7HFVeIn_;YDwolV)KT7lmwBU!BeiKvBNP=#tUU`ZedZ4qY&GAcTTT4V=O@142CI3_s2~YwhPiMtKMmQA$%O>Uh)G^eO>7xpF!5$Xq!~uQfFQOg7U`8Fvtg8`0tf4d0%Z;@=TU2etW@O4CBSj*AA}KB%1BffhP=+?Q+xvjeX%1;~bR{)2>{TgTi?Gh$d-{Hy%@<_{6Au6Tb@o>RI!cllS=^^V?l0uv_jn&~UN6t9)|5#6$<EYou{mHdgPQDCBR@y5HQ!#54dkxAFK@-Po-m>YlX+g0Z&m!a8>}VG=KA(H<kUnSN6t+n4U_l1PBaaAf?-!JC(lWe2I85FBeuH=~4uZ4@_&aZPpHJ5O(0Co+^PE9iPOaX98<G7VXA6!A~ny691~sXB;&@VUEc<cZ&ng(+HCVzet@L~G4(HLyhXWFIDrjhEcCnn1TyX=Tm+g+WdGUtXxfc#a~~fbvi6??s&^2aVc2V%CWx+cp=CV!ny|R<(tazi+B@w(b*3_?naVs~^ups|T^^qG6TTjcz2R{I!+yS`bPGFfw1TmOs?TRT*x-Wk3Q@y9P552&-V5^{iFiw^a2eZz&b~(Gc2gtt3<+=lYYPe(VSI6dHNIr=w~@CY0cVIfCG8s>njFBg(2Vjzyl>-bRy5RQF{@b9}gFOD3jl)WrH*qIhk33gK{)Opr7+0<^)D*5^QT$-)f7xy*OHG|(kRdJ8bMFm8um3O>K(N}bKWNM?YTtg70%y{ufNOaihc2s<hGk33L{g;(;}gk+hdO;~Xi-B(b71$0=W&JLO`Opm)~s46Gss3$F+yCN|dE?BCzz!es+c5L0#+KQ}<INU(n=zIE0oT*tXDhk%jZxJbYu3A+=wKH=Uy2wmKd$-T9EXTrL{Z^mvl$RpVB0lLU@U&rd&CYb1!P2#io+x1Ko4Jc70kRZDJHRyldl)vM{j*~GdiuZ8-TFi_ITL?mJKw}{LK@^(5y+)sq6@KDd7f)gvev&nlyv#Ou>Lpp4^M{ZD5EY_;^@(UuXk}#Yp)XL8>A#!pRGJY5s|Uq9tCu~{;FAdB^FU;3~Hh?R#t0f3mUSaniGUlNv1ya?*=v`TI~Sb1XBH<MbjyjLkm0cN%0|_8Ne_DKtofw)IevwKWb!F#Z8h2-cG^ttJqqS3<{=eN|UG0!!l*h&25hD)KWMFGOP+OlmBDoJz^nej^ST400RM-Y~4DD?&YY?kl1K7YYZmE;nIFwUV)q~Xk8;T3o|l(#ft$(B!ERni%>`h9vaE!oZ%otXJ*({J;QsxYH=NVRWc6*RpvM5uE33-7Y!7CJYqka`>@%B^?Y0-meS{yD3fkV6>86h>%a-fwxqQPDag$+-#Na7_nDIMSOBBQ72?!KlJaU?PpU7@YCr+}HcHS~v9sXa>ajRU6w+Xa1o4^s2T{Q0N)R!|eGHh*T+z7_#v{%EUhko-AxtIC(<<Dd*hNm@)T3?iq}6aG5rxH9(l7&0k{mb}-V(B=v51`8cp?#C<gNhs<)6PNXt`$ksdD#`zR124u%0u=>BWOE`Hkg75x#&qPbfv9-N>KK?A)INqWrkTliYFa@2QMG?0(?&4craC4E!Mz+CNFvMXc&8`!fO;y3&htS^E~u_-ZVX^yy$)d(CR=>9!M8owDkl7g#WYHBMctoCOGBU933Wx=#+_NUKl04uWF&D;_5iv`?ATMvMy8AvAnwLg|iCZkD|jXt;dAYXm`$P7if!D*|Wit&4LHwCXG*vNQV!J-g%5z8EuIF8)e_hT9NbVE;z91HwOLGN#&2k!DH^05-6wDj5t5V(MKbk7a$;Pijm{$(!JH(^Ns=h+=5kAK?rI5%R&%isD=EF`;KDa`?2%ew}0ye6D*_oOWaJEPjR$RodnwHhc}z>DaT^yn84B#~s8i7#orIJEqg_4p)%MpeA}k@moJLo4S$Ff7OYP-3Y8liPbi5l>nk%R|R^WxK9bVj31DuZTJ>4)PO_yz9>}v%A7vng&GjT=J4Vq(r!2So`8Jdfu!SC`}~bceW8m;J`y=VwMF6X6+cmhf-h1SIfLP&^ffguzYcM{{eZU(ScGSXjeC|*5Xo63u3dDF2GhZq>4TZ;lTL-t!}a+eCR6^H?ko>XpQaC2;|ma>zgzSYoq|<Rf6m@6tb=P-KE*}%2AC6JJ#W0uCh(oVlCmME!9<@ZPQ}6DE}A;?;fz0N$K+#20nNTQ@TYOW1((o~m_uuqKw+MI8f2EJAMr-NgF?dVEFy(YicCcE1@qzTKj9N#<<GQN-?Y&=-Hc{cBV@%!`vTnLOSW6cXVW2X+2|DbhAiS=ae&)IfAP+2+W6vMj|>f1>OlKO%XumQdsY1w-0EF>cR{%>n&pujp@4qO)6FU^S96(0#p|9{jCSHhq_l9d{r*!|e{yC19UP>q3}7fw0HBssZ~f`nA?<8<g#X0OLWQsy{);_T3<!QgO0k=3i0`wPe5XP3t05ii3UW76g!R{YE^B##b4A3_T%g)_-SU%LVz(t2Gnc=SIFd9Qlc^MPdB`4ZiAs1dt~)WS$3}0fXD}^5ML!PcW@k(`Z6obWD~(aT6~5<=|G;pA=w9D1ZWk(UK0WV!Ae~rK-{hC<eYq50_GHnX0}SX^SdzL$`iGX8R!{)2+kC3C3Y)x|;3q*x(1Q$A41(YULK^t3Dq#0%;++8u&O6Fq)E1BrP;cM6Lq`2GNcNl6SDR&!zUyA6U6)5~H=n`WofcGux4))ha82VSMEqga*>imDhN@m4=zT~FjwHFayFWbI%yO-gWQs|D71n$|Mm`gZn~iXkPMB@k5|W#-J&&_{Y&3gomyp;_tnEwBLX;m&$RiJHe`@4HZIdG&2Ym=wf8dEhR3g!VwXjJ(TlSFZ5I7@r<d3Gs_Us}Esf)&OyO_nX*-8;l!ln>dt=lP>g$m2!M#!5!zb_(6bcrTe@32-sFA}MTac1O&9kqrlm%&gM2Mt-^H)w91NjILLzC(Sh-7YM}50x-%v(<y%1c6v3GesihO=IA1d#iq{^=`OriN(^KTt-Mu4qpatuuAszx9SymDAqloY(eWfF^q~;eSG8?W6ID~q7R3ldvc6f-C8L1*yQaKZeu4?BBr4>YZK}Z|0x6;1*L(xJe<40`paXtZ?k>k)lw{XbuWxn7&r(~egXFDF%0BlOsmo6np#V|>Xu)UVkv5f9WLeW@yrkpFH@_NYDYwUT9f<HqP9BS20M99!+Hc@sRwibK405&JtRfaljp;I0BcxMzsf^f!#!M>jqYuD2y6Rp(zLq;epeADi%)KE_KO|-j5EPi?iZZ7P+q&Zi=8k9+O3yyi70A^tmEJfHEBX-Li4MpN2an@8`hJsu3!FGH(0zoB2qh4rHb|?L_*KRoTt~_VB(uxebV+#$=HBjC^wb97YL}5KH~EP9yY3S!}JkS*Mx=zw}~6zo$WS;mx19wyp<DkLFoH*4o<|?Z+}oht{EKmgW4j$#FFk!z`DhSaQZ#@I~s3?YlA%-n`w@X>uQEQY^FYCgkaz|6R*^@sz(nHk>QZ@36Gj2Te&W!NFhc}qBqT#3GjSnMBqw43CEs3&lF*ue3(7-LbO*E9#leos-1RQcBS7&Y@_lvgIM1N?_})Pf*u!-QqfiZ+T0b-KZ)k52wQrmXz<to*kn`t*oOe2ChT>$Lq4)phV6MTVH4v!NKJGdwrMC^%tkx0&5PvUu5b@b0;US6wmTjYXSSfWY<!|*V_(IC8zyq%yRcEB<b!07`;eI^7pd4v0CRKt%j#@4B7bez2%dJ`@_{=<3Y9q2Su_)UPNRL!A~CcSA#90it8CtSRK23-?{_c+Z;Y0c{R|UVs6=SYi$=`^)g5{k+Yi#H4c!LO+|io^2P++%FxN0e-TH}M`cvTfv;xd8-JoNGEFcJ}bu_UJA8;Z{;XeFT`b=M{F3Khw1AoepwXdc5xe<Nps0_|7o9<s^G0K_t4s%C3goqgn_5`hcH6eO>AtezZTM4nM1TDgQ^G8qRW@syYS6hL&^rzWZVtzG%X$fIIlX;HRTa0y6<NkhE73pT<xxe86F}+QH0$qy~jGbUeuZaNekl{a(3gT*Fn|~N~N(o<4+P1FkkXP_ZRg@c*Yis}j'
FINE_REM='{Wp48S^xk9=GL@E0stWa8~^|S5YJf5=0c4C{9OP5Z~g!OqqqM@K1daEIbl(_k|^e*h6nlX<Q{_4FU2}M0ZXxd$%=|iFQX0Ff2QM%BN0>JSs2B2isE<pwv>MiYQWxfNn|e?@}mM^m;>s{e?MRi3Z~z5cQx1!V9w?bWQf~7dS#8NI~#%3+3U*7awftWF6`+{TR+8RRT5{TR+PzJDFHaixF)?D>VH^Z-?cmNqp$HtJ189oVHz*I0-&~i5I|vv4isfq`|{SmBVxcZ60d<)Z0i1bZab1AWz^8+XeN{KSOZedCZhP(j-c&Y`?gnD)IyIkU5rcePp7ek&3vM>qeOJUwV2p(Z6A#fl))M3WbXg~UCiV2eZOXb00Ec+j)fTj=%Uy;vBYQl0ssI200dcD'
