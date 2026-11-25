"""Active -23/64 first-pass, gain join and starting-majority proof."""
from fractions import Fraction as F
from itertools import combinations
from functools import cache
from math import comb, factorial, isqrt, prod
THRESHOLD = F(-23,64)
def ck(q,m):
    if not q: raise ValueError(m)
def energy(rr,uu,M,G):
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
    def jensen(scales):
        n=len(scales);s=[F(i,10000) for i in scales];K,det=invdet(gram(n));sg=[(-1)**i for i in range(n-1)]+[1];D=sum(map(sum,K))
        diag=sum(K[i][i]*s[i]**2 for i in range(n));off=sum(K[i][j]*sg[i]*sg[j]*s[i]*s[j] for i in range(n) for j in range(n) if i!=j);E=diag+2*off/(pi[0] if off>=0 else pi[1]);L=sum(sg[i]*s[i]*sum(K[i]) for i in range(n));z=(n-E-D*u*u)/2
        exp_lo=sum(z**j/factorial(j) for j in range(12));cosh_lo=1+u*u*L*L/pi[1];ck(exp_lo>0,'Jensen positivity')
        return down(prod(s)*exp_lo*cosh_lo/(2**(n-1)*sqrt(det)[1]))
    scales=[(4566,3083,3313,3341,3316,3025,4740,6506),(4578,3073,3330,3329,3267,3337,3326,3019,4746,6511),(4580,3071,3335,3319,3283,3324,3277,3330,3330,3018,4747,6512),(4580,3071,3336,3316,3288,3314,3293,3317,3280,3328,3331,3017,4747,6512)]
    tails=[jensen(x) for x in scales]
    witnesses=[
        ((4750, 3131, 3354, 3376, 3408, 2928, 4745, 6634),625,(4398, 3215, 3330, 3406, 3380, 3016, 4635, 6201),424),
        ((4762, 3107, 3374, 3352, 3312, 3361, 3411, 2921, 4741, 6630),832,(4400, 3192, 3340, 3388, 3290, 3396, 3375, 3008, 4627, 6188),668),
        ((4766, 3100, 3378, 3335, 3334, 3342, 3322, 3349, 3411, 2923, 4743, 6634),908,(4397, 3183, 3341, 3371, 3303, 3382, 3290, 3386, 3371, 3007, 4624, 6183),808),
        ((4765, 3095, 3376, 3329, 3337, 3324, 3342, 3330, 3323, 3344, 3406, 2926, 4740, 6633),1025,(4391, 3174, 3336, 3363, 3302, 3364, 3302, 3371, 3287, 3378, 3364, 3005, 4614, 6170),1028),
    ]
    for i,(sp,zp,sm,zm) in enumerate(witnesses):
        n=len(sp);K,det=invdet(gram(n));sg=[(-1)**j for j in range(n-1)]+[1]
        q=G['lower_signed'](K,det,sg,sp,zp)+G['lower_signed'](K,det,[-j for j in sg],sm,zm)
        tails[i]=max(tails[i],F(q,G['S']))
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
    motif=G['edge_motif_gain'](); bonus=F(motif['gain_lower'],motif['scale']);pair=G['adjacent_floor']();cost=c1-gain-bonus-pair
    ck(U5<=F(73670242301,2**39) and gain>0 and bonus>0 and pair>=F(27612843111,10**15) and bonus+pair<F(1,1000) and 0<cost<F(12282281766462836057,10**20),'active energy bound')
    ck(f4+f5+f6==62 and n4+n5+n6==218,'face accounting')
    return {'energy_upper_1e20': -(-cost.numerator*10**20//cost.denominator), 'first_pass_upper_1e12': -(-c1.numerator*10**12//c1.denominator), 'gain_lower_1e12':gain.numerator*10**12//gain.denominator, 'motif_gain_lower_1e12': bonus.numerator*10**12//bonus.denominator, 'pair_gain_lower_1e15': pair.numerator*10**15//pair.denominator, 'motif_oriented_lower_1e12':motif['oriented_lower']*10**12//motif['scale'], 'motif_cases':motif['cases'], 'U5_upper_1e12':-(-U5.numerator*10**12//U5.denominator)}
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
    if before<=F(503,1000):raise ValueError('combined calibration margin')
    out-=y
    if out<=F(1,2):raise ValueError('initial strict majority after motif repair')
    def ceil(x):return up(x.numerator*10**12,x.denominator)
    return {'initial_lower_1e12':-ceil(-out), 'initial_before_motif_lower_1e12':-ceil(-before),'Y_upper_1e12':ceil(y),'cross_pair_upper_1e12':ceil(x),
            '0110_upper_1e12':ceil(penalty),'overlap_lower_1e12':-ceil(-tri),'11011_putback_lower_1e12':-ceil(-extra)}
from base64 import b85decode
from bisect import bisect_left
from fractions import Fraction
from lzma import decompress
from math import lcm


def words(data):
    value = shift = 0
    out = []
    for byte in data:
        value |= (byte & 127) << shift
        if byte >= 128:
            shift += 7
            if shift > 14:
                raise ValueError('coarse witness word length')
        else:
            out.append(value)
            value = shift = 0
    if shift:
        raise ValueError('coarse witness truncation')
    return out


def stencils(axis, gap):
    nodes = [axis[0]]
    for j,x in enumerate(axis):
        if x-nodes[-1]>=gap or j==len(axis)-1:
            nodes.append(x)
    result=[]
    for x in axis:
        j=bisect_left(nodes,x)
        if nodes[j]==x:
            result.append((j,(1,),1))
            continue
        i=min(max(0,j-2),len(nodes)-4)
        xs=nodes[i:i+4]
        weights=[]
        for u in xs:
            y=Fraction(1)
            for v in xs:
                if u!=v:
                    y*=Fraction(x-v,u-v)
            weights.append(y)
        d=lcm(*(r.denominator for r in weights))
        result.append((i,tuple(r.numerator*(d//r.denominator) for r in weights),d))
    return len(nodes),result


def spline_tables(conv,Q,payload,gap,offset):
    w,X=stencils(conv.main,gap)
    h,Y=stencils(conv.y_axis,gap)
    data=words(decompress(b85decode(payload)))
    if len(data)!=2*w*h:
        raise ValueError('coarse witness shape')
    data=iter(data)
    def expand():
        table=[0]*w
        horizontal=[]
        for _ in range(h):
            row=[]
            for i in range(w):
                k=next(data)
                delta=(k//2)^-(k%2)
                row.append(delta+table[i]+(row[-1]-table[i-1] if i else 0))
            table=row
            horizontal.append([(512*sum(row[i+k]*x for k,x in enumerate(c))+d//2)//d
                               for i,c,d in X])
        out=[]
        for i,c,d in Y:
            den=256*d
            for z in zip(*horizontal[i:i+len(c)]):
                out.append(max(0,min(Q,(sum(a*b for a,b in zip(c,z))+den//2)//den-offset)))
        return out
    return expand(),expand()


def as_tails(corrected,Q,n):
    m,p=corrected[:n],[Q-x for x in corrected[n:]]
    if not all(0<=x<=Q-y and 0<=y<=Q for x,y in zip(m,p)):
        raise ValueError('coarse witness simplex')
    return m,p


def first_endpoint(conv,Q,knots,fixes):
    m,h=spline_tables(conv,Q,knots,8,1)
    return apply_masked([min(x,y) for x,y in zip(m,h)]+h,fixes,Q,conv.msize)


def second_guess(conv,start,Q,knots):
    m,p=start
    dm,dh=spline_tables(conv,Q,knots,12,1)
    h=[min(Q,Q-x+y) for x,y in zip(p,dh)]
    return [min(a,b+c) for a,b,c in zip(h,m,dm)],[Q-x for x in h]


def apply_masked(raw,payload,Q,n):
    data=decompress(b85decode(payload))
    width=(2*n+7)//8
    mask=data[:width]
    if len(mask)!=width:
        raise ValueError('correction mask shape')
    values=iter(words(data[width:]))
    corrected=[x-next(values) if (mask[i//8]>>(i%8))&1 else x for i,x in enumerate(raw)]
    if next(values,None) is not None:
        raise ValueError('correction value count')
    return as_tails(corrected,Q,n)

FIRST_KNOTS='{Wp48S^xk9=GL@E0stWa761SMbT8$j;4uUec3l8UAHIP2V-%qrf5yTD$r*0tV!XqP=(vi2%sG|2uQX%y%n#=z4PAeo!AA*?0VY)-=TauW0co%){l}^x7~DfPcrT;#WU-d*>ywif<l+Rs^kiLESsgPr#p4x0KElFFO(?|SUZi86{Lwl^ztz7~ch_5>St7V*3zz*5Mo9^@^zw=-WaMqxHfsx55-OZ<M3_@vWlAhrvQaZs{Ytb8$I#GNss?o9G*Pc}e&vkU(!to^fp|nJZX?YsT?gSTdcCUWN9s4@x)-avP$q2Iw}_93#~ep{<?7$+a{Cnl^z3-yqxcN)Q1<GQNxLnk@RbfNsRFrn&g5fJ0j-_~-kM4TK5Ivh?1x?2f(^g)Wq<hNgQaPoX;qXoJ*ft3`%3Q|FHaYBD!R4c)GR#|siKCJRWCPdebHT8CCw?dGJ2bnC@!%^oCp0|$n6ZdSFQw@-hg~id615!1K@ao2@=BcUUGS`%O%K)OhVA(Q8LZyQJ8<ld8n2i!apxq7I$9|eJ0gKo_00x;wYcjz4LyP&Doz}O*nk(bQ{*ViERF6EApM$hX<LFQQK?5%ODb`)rapp^Roh&c{f)<Ud^TmnLR!2k7rd)9E!=I-*~k0v<^^`rGRI<FaFaxi#lJb^uxCR8<)oq(uJ0>GM=fm5KA4{(K+{CAFnBGMn&4@ZY{ff%`GcNVjsdic%Z{(eQ6#~l}J&$2l<@*WZB~~+}Ekf(!6_G{E}13{IGO+xz8OwgbW1_L61EOTp-H{n+YcxK@64Pn2>S)xWm13hym;AQ}Ku*uVWz{n&rb$D?o(hvwXQAyYgrvjp`fu0;f90ZV%GYD%c5;yUHGc=}l{=%ih0zr5Ie+h}rp7hKKj&&G!<3euy5jYRKrusZIp^^Vv%Kr_loJY}Lj>;7lIzQRu)(C|=;LP1B-yDBE*WR}3yo$5noS&z4=X@SJ>pP+VM%jo4Aus?6@2QF0CyYL!nPJfOth<>=n^Z4{+DcYiz)Nl3uK=P+Zjr%}M3U#c!ps*qx_#iz19eXU;ykgKy~IJUg$5DmU;uf0>c*vqZT-g?Tt)nl-$TJK=G6!P^YYEgB4gmef+G9w^miqO;Yw*EOrL&5O=Vd?U%i>0;MG>2#A`QTPrMQS!S`IZsfEg(17%AmcUzwkG(JCXNLa_ZBuq>R&xd0~c}kg8ixkQ$zisCiyfTP9I{SD~+1=%#F5n{`G6gmWj?!eu0y&_|f%&%Q?DY1<^0+$2<(Q+q<yk`xleOU{z(YM3{S)<np{fMv-`4YY&5_(OxaY>Jk3r%_`Dbx0w}4@|ALx*{O7ZG8Lf8IOG`n`hsMf?=PrngQsz#_xnru_1wW920q$!uFXEC5QE!==!44JmOnoAuA5#Io-783qEso^}4wQNh{o3&Z9Q56NPF8iT-*LRu#;#osu^}j?>>)H(ove#DHc;3zf+4=}yX<J+3;SyLQAKM2`iaNov10jq6?XCpBnm8Mmiv$A(zxjX>^LKV%#V<Qdn8JnIJKjc>dkAFb%Cq1xRoph6LsbbbBxfnN#v-ONgV6|^+pX*=co+@H5oE3#Y~S$g+Kx&96UZMXtUZd8=*!QE?F8ULb?hREjeY8#<>V-5mW41Vtf!8L@tY~UJMcOl{!JQ3-dCxCHD$_=vABWKkB3T0y;G4sWuqwj5s;7-~%Mm-YGR8K#>pB%8>jp-@ly*igLfHCH|uF{jz+KaO8^FTbxskz>b&FHV~Q4*%iY|YT{f>&Tv7Sj=F1f<~R18s3uI%qIpPY}ili+vKjE(%j4f-qFsL3Wk}u>vbFJ>unmL?=E+h-b|%pvV-~TGTu%kNo2k6PX^T7uipBW-)6T21>)ox<f-IZO&G*w8Z$|;SJkg>j+XP1kf-%c6q}ieMIwtH5L3c>yBY0i)-h{Zh?Jbg0(H$OA2z*zXJ4MP2Yl5Ze4X2d4|fgb~vZ2GKmsOc*3Ebv<r;w$Kd`c<>AI7I&TX)UpJY=!ux=Zu{|7l*<OpLUrtlEVSW-n)R@zu=xHRqbMT8&bUTNLPUj!Har?=c=A!mzWpV28UKP*7-F>~_g=j5Y>ur2$@;{Arr}X@Ye1G5b7$`|`g#t(c%Umf7!yKpQO6ljep0A@(y_GE}fHk^9-I&Cp6lPK2#bEZY(baR4>O{L_X>D?=T6H`XMcNCym#C}cGavo2k{|@);Ve9`-H=7wWgSHP-RhBktp1~;;2-#M#!?L9tddqfU0FY9uv7H-zTMLizFeOSR0)g<KR4!!Tfx#JDunig^Cl_<fZy8?7Ae}yper3#_liMtSYuj-%(8l^{F(~vUD|o<TfI&ZZK8cRiI)9(mseKcEh(b5&owwtZ$9DB@W|C|3sG_*>2k_$+MWUp8%>^FtW(ay?4~er0bFfm&_RTfu742Be{6JB>?@BZYBarNm@y|~^os=rKY%rsZ2nt-*Bal$!dVpl_!XSHj{ewMLgro&d<9*?A^gZjFMFk)?IZ)Rd}M+64UAY3Kjmuw{2a|AJeEnEPp#`Gn&ow7p2G?D+y^u7_Oh2`)Y!-{>!T!XmF1`inS7D%XT(%xSB2dhLtW1z$@91(o2V|vV{)(*=%OAXsWzVTkbFebnLMSEs%nI*#5065Eq#LRYvm%5PhPB``>h={qAWH7KZatgM|#&;JgvjHKxtf>0ic`m84Gxdjy2(SSbIU_u@-Wup*U5}eDUe2+hDj*vT{&z=hC1YGF*XY8JfmQ6lo`2W<Jlde{@p=KA&IouF*ZF1LjrCQROk1$F~^H(eA@oBKJ2&<vv7vF0FWx%c#$0D{HyRV<lBS0R%Ul04I1(OJ4aVoJX@h*>HX3)UVL6gyZZq$`Y@1r&agy#Sfi~4*#`u8DkRFlg=orCq<TCX8dt^%6?h(0#oKep;EWXVd~B%i7PkX#1^j@E~V-BFSw{s$x$Li$GrbHaX<eTvIGPg^z7LMvU{G&49T_uzH7}&RiGwcmJ~f+rRk3=JwcWks{*VMJ-qj~8z>B|`xGNL;mfbluJ>AJ)DNULiEr(6L()4-Q@z-#Jy7)pwGv$g8T_v;u)|VbJpLW~5*A!y_4N&=pzBlEDm;t~cqWDZWMVsUtlHWgQP_Eb$1MY>MwPP03!=(hi-faMyDT%14-(O$JFw7L?zpo4CR}%RtyG54@OV-z<_A^py{6UCsQk0oc*V?1gh(R0_L2JiPiFSLPsGmTA=JQw$5lMb59Of~dl}IdvYt#gg1K??clyR-C5Rx6V=!G!sH%GhOVzMc`HI&%ptQ$2#4^LUhYtmf2MyAj?MgdX&kME{itO-1;<B>jGO$zSTSrY{4ZFSp0=Z33C1ge-SL}M{pp0!(&zFO4zQWBMg~G$p2mQ|M048_6PG0QME>L)_OEWv3v&A61scY-b?6EVzwRt)H-yJn=;D#*HJ=m&qkJAn>OI1e+NafpdEQM}P>^r4h3T&^Iz}_O#n9YxLQb))p(Ihst?3;74Yymbv6Hn>=QTn>~lw^d1lr(qNC&hQ^jxi`fs|G&^Bw*}ph`u=xDgEhT7!kCpS#z=?M?szT<|7E44wJy~)YKtck;oLA5ExWvbW=uDD1EX%je4Rqg~Gn?DjLgWx?#TmSQkxOUUc}Q9&$T~f1~x=YDF`puuo*kK_`poR-i5=qCH}5YC<}tbYG-FCm(-hQC+`yYaLvxUulrsgmJ8Vg+qNvz)a%ERCOU1457^Bu-Hs|e|x+aGEFI1aG-fymSiVI*;#lFKW=~AYsqrxLZm9SN%fs~Bv0V&4E>4nUEs?PW3&eI4#~XB649!S=Yc~4_@_FbXiLZ<cicLYg@V?hC$m}!q6q>gwVec3yGn6bVEpn_2IEmMFWdaOvP~6V=UD0m=zoTgFdE6-RX<2PJ4^!BQ7VW?8u?Y+>FT-{!iWkguQdI+z&y@%gn3Zdl8=|U`}5Y%5Usy;t)xA)#{p9A(cIxwP@nH*k%U+mBa~#L=eJtYfMPw(;@D|aNW{L&*?QLw3HrNj#IE+27Y2fxOuB7~xrf$t*cATfcG3;!maZVAi{`c=4;Oq!S4d@O;io<^_Ax`W(zxcKLR}x6m1Tfd@7r9spAkC~-w~kCKjgdYj6m6p1uiVRmN0g`YcZOLQ8l9HNYk79lXm;V4|D<6$q19~3<cIH8mU*G_be(zc8*&}o=TPCTI<9G?=H&D4ev_cXzc#qUpsZ_uJn&41I@3GVYM6$3-~l=hr-+)S<^G7i9B)N&&&0**OnO9Oc{OhJy7%uhZ7`nOy+c6!n;j0UUL){7hegO+#gH48(EcxpChTJ*o%M<F*5JC3tjmm{0*R*q-RG2=^DMmf>=ksjnXci508HIeZXqrb3o_!rM;odjVMXK8fP;`6nlIb&1NDrr6t365LgT2QGiG!TVUnd==7M?3F?H;0<Eq#aaf61u7WJZztKsYgKzjKz!GdLvvE;o&_F6sAq2MFu}D8(fEHtIf6h&zMkl*6VTk%S)g~N#{#<=C8}-i!eW=!1cHR57p(DL>u_CsY>2WCdV#jN6i3DHk8dGTY$d-TDHb+ESd|>_w4{d{`4ugCcv6*A#<+-rTL$jZaZpSYa_-X5&upQ<KnnmY;mJ_ps&PuawNpf4VzzqW~t^)5b>*Yn_TC$o$xF&52Fqcz{P8}wC1dY9+Tpnb`<M_!mC^%dT!%lxKA^+-~Wia;}tIRC&Pf>+oJs)$LdJ^k?OmgaDF`+Fz+a;@c$eE*bw`b$`JW=ozSHvm>2>w+;`%VX!dK{&fLp#SX_zx1>?=XwsHFw(pf$<CibKYcj*$BAOKybA?AF$q57wWMR^aUH)<oeh`U+=ZM0Bv->+LndU{=<_t^K^8XXHwD8e%A{J`OnnM)q@ogiZ?#<*~n@J|GP`7ao}0NX(Ha$mi?%O+nAi-7o#_mxs(Q}ZIg}z_z~!qY1RS*8%#_{y0ks8YrV(HqoWw+>58O>eRh*z+p|~1c);|338s%h)Q8Gll1t99g&G}zVJmai%~d)7iR%rC&8bEJ6lDsXzby|rr}{a#T1~aw(o=c<aC%3H*Gb{qfal_Z$Xf|slfK7~tvzV$5xr1gX2+Mmd7_=9EPV_Ar{uSza^ZDFUSmM{v<6)Nf5|m2sO#)jq%LZ#y>(3WBRm>d*^2@9YriLsSk}-o-nPJ*Uv)0(Ot2N;=!fjkB0#L=fC!e*LUH_ECG@KNL$_tiAD;Uf%{9ewGSMhs>AJC$Vq%AQrg*}XvW03ox#~brWAH&}ALq7Z^ZM44fiuA<n&arnWP5TU_iQ?ZTqr3cqq&@-gHRHUx=QVWNl}aEdZ-P*e8+dJM-+E~gU&m1o{7}M_g2=1lPt?Za)o`tre6SlO)5#usnHR*{G)-iL+<@0Odvdlbqt5N^|HC(F!ViwPDC@|bgve9FJh<ve{AwLEfp<(<cCB~j;2IHBCY)_F9!n3PPE0+gnJuGn|;$;<6R`fLN>X@4B`~1b2FU9=Vds%!T2g`5F8u!)3X4J*6US~V_-i=-4+Q_`wN$31*17r*HJe%yfQ1jG3Yh0U>r|U2SjY+A2yZsNy)_W&UJ3o*qs6f<W?HHeC=x?p&|jh_e6OdC-*P%)l}rTsh(aZO`Kjzw+K(-sQN?R+}5<Z01VY@%f_I%O?KWH=Mu6DmJyER<p5XWzkz<3AKs9(cKM|X`ui40*@+YAVt0(ksy!laV^)0au_q?CnsCt8?;`YV<ie1q_K{+ko&)WQ2dasv=Hh;6>m=%@M!5#YJ9o)Bk2`^f5=8c_P(+VPQ$@`gG+^ux31Nx3EV!SRoIPIhT`{G$aPyMUGnDF(Wl+eG6%y46=RN&0-#}Z5X9G4|4l6X8k?!3;vLD0PYg_X}R|{9^-Ph~=8u5rIm^I-P8tp7RfT(8BKln?8W5t(0CxX1Z_)_#HejV`nPa~eyw|Y~^bd0vG1Tz?t!kM<|W+MeXrPIiIJsX#~JfI=hlvcpR*HO(#T{=Ybmba>O&MqZoKsgS|-DeyD=ePx&9e2L9Z|2*;H~_V5KS%cd$W>;$Z4Z6hc17f5#@q%!wqdt6HjU&dIDVmF_$6(vYc{el!VyBW9=Y%r;<Z|ipYO1S3A`~-RQwn*f7yeXN;I_Nc>n+aS0a?T7D%2200ELEg<=2zrkV~5vBYQl0ssI200dcD'
FIRST_FIX_TAIL='XZ|QF-iak!_5)ipuq;_lSaT>cFM?y2^9=qzCq~?O!h}j!LO!1QhsdJ&&$%G1sukDcWNcObPG~o8%trF+!t0HTN`zd+bZaEeqJB6sUya~>|J(8p8NR>DVfe>g<dNm(tkLXI;-s-YwFS-fKgt`<g}dqR{`(2=KbOj%bg#QMRwsRt2r+h83-qlF9;qyxEAlDex^_JsvLw86{D7yBpMOoO=Y1!Z;l!?{y99z^#)jW&Wfla*qOC1Gy|=5Gf*-WM=f)`=Pzl;}An~LO)f${F+rv#LL#dytP$10?-ksKSU8q1L2@+Evhn&i(h@+LI_aC3B#Nn&*aMs%XiX3izTFT^ho0FQFuA=xoXAYAyAnL3uSRr+sCCX6<bY;f6#ly+pCME#-nBxt*z1`9}P$nu;mAXQa)TtB?J=K3QwpkdK|3DhXNs*)K+U_ISGd@}mjna)_UiHchH}r38obpDNG*q73RMZWPqdwL0FLf44tsHR({Y?~(wIy}ZhIZ8DP;Xdr6Z;*O88}pqwb1{Qq1c+2ZA}nikuESE=L@Cv2p8N}5h(2JH0A6?=wY}#LEa=i93*UY&5fdk2XaxAh4*9c&_2co&5|3j9FxhhG(n!4ps^twBXRAOV883y`NI*kL7dV7#&dJMug+1^*sr`AzTb(K<-luUJLhXV-;0^RSHfW-H5x;1eiYOQdX1Av+)%vWUd_)qKrxFVsRPJL$=<St6sU*v$AB?a<U``p@g%Lifjxtx*n%*7S9{-IUDSd+c4{ERn;vVTXsc|L&b9RV^L}b#3KRz7szzD24oXo7aLVj9ggjR#d8o%b3`-_758IK?9S7#I6XH|jW6<AwTCLFz^R@JhO@;i5R@L=%7ZG!uz`b-m&X_!qU%2^Z%2YXIrs%bM)`)%vtu*H?!5HBgCT@umwUO3z=;@Bv5Y?oI$dLylmi4N=3(^zkBTkjA?ZEUaTH-_7(?a8Ut3|~dShbpO_qoCBr2p*8S^Q316$r`AzYA{N^jm)D*6l_nBX+~vKNXs}Ib3?XhK1!BNkD1(EjeRl=`~aOGZM>kl@32QK?yiZcUo@w<A2R<Ziu(P0h;r%NSU{C17`UUwqV=<bF%WJMaMCsYJbuLf{@jmc1t3jm1@vEPQoX<GJjh)l1^WDJAO8^6ycOk>KI&vp(=14L7#nv$=(0}|F+}+Se0;~00EPP0l=gMCEC-EvBYQl0ssI200dcD'
SECOND_KNOTS='{Wp48S^xk9=GL@E0stWa761SMbT8$j;1%lz(p>-n07jCD+^|?1m>KT$s{VyJfmt8H+o^1|O`ZO&C&hR-Snqb>r_G>3o(}n5`Qc`*>82Lx6PPHxcKHkUP+1MNlWhN^I3X)Z<TR_(!rnD;ICgYa9fijTsa4yjkz+<#<;2X>9m{W2D*nFVZ71TU%PKEwEg&&>A^6KlEu_26swUbZbP!*l>K?vJ_K+SSvUU@C+z2=fk!buQ&vqEX>ktu$;LtW0_6x6V+<%vD03T*K&!|`OvB8ppht$V`o0+p{&3+P#oCKGH=xw~B6?`D-f_te8t%u5cxwiYfgQ)6XB^PH0P;~NuZ2C?sFT^BW!YNR6w%rQSmWVBL_k0C#QjHgTqROv7v|9qryP+l{!3-mav&Jk9<EMevK=1DSMZQVqul8TK#4{bk3lSf~pV~SFT%`=JrU|2n-29HGotxi?I{9^D>k<Jpa+h}RHE5$(`HQ68Bg9iKuQm=YsS~Pa{1EI$1cq9&GFG)BS2oOD9+R&G5s<^YYh@-O$Xk{5xJppN7!8K6*8lhEjgil`&zcNRlVS*p-|zu0PZbZ-1r9=qUcVsZG_VGL<bQ`Ls)IcP$;a-kq;C)}k>5^M!%uT7RMWFYc}J;A9ohlv!Vo}lmyl7vp5Z20Aa1y%x+@R^<&sp7M%neW5*ZfI@1IL1RY(H;w~Bqg#(qjND62=GFG#p)#_kE2PCW;*o(NHLvLpIJ?+^q(PB9rpj8(K}+!#1|6TeDArg?9v(7o=@rI6_?dI#yF-yN-3-y$c6eL-4|Bw^-tjzDRh=|fqs=e!>1*_y*`*Cxbz6I<mt8|BJBrCwAS&pbv#FSA`X^!L)3m!%)?hx5>x7S*wqG&O8%?{>t<p&HC0kA2Z|Fetw}!%rC0QQpm+R|m*SxRo@8a0!9AlCcE&VIO-s7R7|79}aobpo4o9gS#v-OBxt%N?%$vlg^<Xu7QefuouXBFn|DkpY?1i4^{}B+kn$T{-E|GNPXdI<N6Fw7h$08mWq<%w1=|9An(pNLfaW`Jbf9h5Mn-O2D}sFfhsDBR=j;kT;OeS7H5{I$VA_&7Y++qTG8}5QFgvQ`URrEoLZ9^Ta~4AX6|eqvlKZ%m$bX%v8H_*MmRKJ*&4yz>@>9XtI+N7-;3USi13-f!Oa~|?&4AJuT|*EB$LYEtH#5G{32T1{#>80B!wSC`e>lz7sRcX09GFeu~H602hz0GGMi$NY1}npxG6#6;;=ObMfk2-<@CP$`>Ln01yK{QO3>6K6y(U?bRXQ?a<$DR&e?$~8*5ngrJ4-Jk;Rk&92ZIH?O(#1#g?|fP`gewyl5pS1=n`dRi@F^;j%i2xbEI}&d7hc$c<n@4g<(nhq7Dfj+m32ck~h5<HC2gHp`qH8ltgKjzSUD_1)ejy?5dfPwK9!IwzYNu0G+}Yx+1~C&7C#n9C`+{n0Pr6jWuxqYB~7J9uXaXDy^79v4Mshn>@kV942#0a{9+vrv;Z=srTvJcb9!reom*6H0&F$kAg7Ms{|FHtJ@k0cqcbWcD+bJNEQyZREyAsgdTV{R+hzX`b=U2@iAz4oWdPsvytS7^k~h$tsBxmCjwVx=qPU+B<Y+%GilTBIT&Xm3f<FA3FdYBvs_W>r~G2@7R^kMJIBl70{}INMvQx`{LMs$m&)pcyFScGDo+z*5gU(W`G=}-N~tb*J#t#Oju2nq=Elc{RBS*0j<yC7s3;_oDM|q<iBtgQ?K>-lV35QbSlvmVuX8Tc1p*t`>pGIM3c-S)0qEv0EJ^iC3=<r<oRC$9q}hYGr&4DNG}e^TL21NbMhuJKj{FH#5>(6W0i_HL9WT^qiR#y8&Y2|>KE1|2`=qC{n(~aW-DPSE|F7&uKjc~8b{|}%EdazVXtCX{VEH5|JZWAPpyPpYTk&$Bz*>F;O#onzatYa4JWxDFyzLxUG@GZ<ja`Y2ZD;VEWdfoj)n%!ncmN~e=f5|zp;}gB`_G}6X<v=UD%{-G;aj-*6<q3{gpOAPN?gohZyi>2coD!Bu{Tx=iLXiGfVwXp$ARVp)6Dv>2{$CCD8iXL}8H2)%6$<6L^VC1xFdjL!G77V$P8R<`Coo0QcdieYfSApyon6UhPngeZ`MEvR_8-g-}Tr&@k=Y*0ZKKbO!s!H${-g-Mhsul@lTIWMDICS<HxwDSJtIVjtULFGf6s0nA+cU>$TyoX+zSqdXNn<?`4^{c-Z;I4&Xzyp7!+$9l?M-2X?BS#&C$MfY{u$HG$A5CAJL@&Et;8^15piK+S*00Hg|>?;5OA3Zw1vBYQl0ssI200dcD'
SECOND_FIXES='{Wp48S^xk9=GL@E0stWa5C8xGsBn&k;nw$3`vU+0kckTwLWwZ5V2IAOPIz1`9gpatJb?!!>brceQoIP(yJ0P_Xz211=>>X;Nv76%>Y`vVVw4iLeZ~G`jRS@=UTU!LjXq7-9NGdAxc`M(i*FT`O_99_5I3@_d!w6e5SbKv-9>Q_A>l+^->Cy4BU9U&i-!rN^l-2sja>LRPl`wjIO`rwm_%zd05L)G5;bd%eF(6QSwgldc2idFPu5QIqe^D>ZJb=#-FS~5pjTDcL#Ta64s<a{Mv=D^eVT1S2hz=bv%<X4URM63DP2Wq)>0MvC7UvS0jztOX)IQ~Loh2MG1z+Sj2RS@D7G)LMl^CMK1A`rtE^ftIG{YrgW>BhwBBfJ<d?FYzfxw${7|GZlR6`PC>S0Lc??ox_l5b9m_Z)HF5sj<{qPaQU62zle~>EFqc>u}Ax1&@hv8`qoeL+(l@reE(kdx<zWTKB%f(G^UQpBAmgQ{i10BL|9QDdA{N(!v5>~+UePwpqv0<aLYFg3Ij#DPsj0{=Mruo=fE&v+C-~LSnjb#agJUmjm$ma+(+%cy%;BA4nZzl1`S}&LaXs}|uCpo}DRXJ8*+Ii!szArnwca05nuj-KunqR<z=4jQ^9z;sQlgo2my{43}C!CKDc~c-HB)vIi2g(*~eS|KJ=6;1Ko=rK2;wA#5;V22BSqSP1?14Oc{h0{sEVf?`7C`oCgov9jRbS~8?MIA+!d$hbHiygNOn;a@U+M_BFd#AU>VR)h3l%Ks<$`fpEZ`Nl>*_sdRKVAx?$ywd*Q`B{a9!mYo7c^6N-Y@a4k}qD8~My1UA2OYDWhBMuKeU4U6f%^RT*|^SOkn&-T0U*$yYenY&6B?G^uv`x7tFO{Q-VE<?jzUp~+86h;itxNMMpI(+OD;fl%Jn%wOYedemk~>Mt*Nxa7wt0F2Wdd-X2MwjskKAN1gYT)hT|+sFSUYqE0%I&faM8QFqM2X(Txs0e$|6=-KW{Pm^(_iR|R;yzIn(g9_bnE;1t{7#klbg&SjV+`MqeiR#=5`+&SIVMC`DQOAzmvsmBKBcn>5-Rzem)k7W1=D7^D->1@Xk~y<yq{1D$fwu@F$X^w{EDC5=$sC#Ac0uG!mmskWTWe+#DI&UX<DvYf4Ahl-cC`qJ=`os5+VQ!Kl2yIiAh-%XDXelq$GRMR1&N2>U4zFt+=x$8YGv_B?ewizs&VyF=aHCwizFjdL2ZM`?wKOw#7BR@;<Wf5{_`j(4)eZPZaL`hsI>t1$o=MVc{u01kd_pj?H(!37sVh{w#MJ28&&c2CP0a(r{}%;;B|d7K0_=HPD=YvVO8o<VZ?rR5P6j+OZND>)H|%FB|Nk_u)7DOw7$3*GL~Di?oFuZAQ|`-^o64f6Dwt91u5*7}V@6NdPoo73cp2_@F7udtc_(Q0;MyVd>|7N~b<er+y$HSJy|sF9Y*s(2+(0C=970p=Go^%opwU%IjR>w))8O5w+BLSd50Wrg6`hv&(OI)TE*DvAiLlNE|Pml~~{-^Ma+j^MZn8y=*q0yuee!uN1rpNo;K(uz>s-aAyv#{E>EhOhwx2s~`t4-a_$K_&p5EXWD+S$QT{X#~HBCg`tT;K-m>QeCW`jKtJL7iu&oHav%LiEn|uf<_@n7Fh~E#LylSOd!G0MT<EvAAh2Z|MgtE?@2)mnb&{apKux|(h!jScX(a{iWqeg`EOG3?lgI?P<3)%R!Ngvnth*xEjWzy=5mOo}qmam#NA@dL40ftgWxQ4oHd11}o?}MhODDBZ7d=pO*B6eB>c-RcoKTUIjw6AI^QM|-QXa+eGl{@BT*wM-J?x#s7&}GO%(ue-ubteXd4AEqBpHbMM>?HYaG^A5z{+io_lU;Hz!um+2!9}b+|r58*{JaD-$8@@v$u|+79pL<T4rN~OqKq0Er0u9xj*zlA&%41wDRf75g%PFR*^ycjK!<1!kv&K)fAbkMt#DLjLnmunswD^5NZg$ReDfy+|depJhF?#W}4D7hFX}#-WwmOz8j@c0`7wAtW8~O0q>8ikgDIS@l7|^BiCnr<H6hF9Je*d%9eATQIH#Y<*pB=nYhSfPA^!+$X#{mSNn^ZkoY6=2OpgYZuuqQ`i-Tt-I(JhqwEej07{@s-uB@^E%!R3P&RK0Ll4L`5tiI|g=ytKvHA~Cz2iexDp<X@I(f506JDCEwwt-CgS`V_maw6l6Y)N}x=9t2NIKLNj`yHd(Yb^_lW$kpL)$tU1dt+J0Kh7G@jIW_h&M~$G`%~m8rUZ{3>E>?M79~cj&OX4nl3xRCe?(al9TCJ(E;akQOmf97&tsS5#%ku!&_Ps<tn&H#rNx@j&spnhI&tLT=dt5@H6Uic?dz_|CtrtF|t&$t1HVHEse0jU1rh?0*%yXg@{OByznClL>}2z7TE8fa`y9mx1pSVm#NF_!H1nfTR6fg6XrxA@^}_On`$tAs!kD=6$Q0_jtc0*3Zda4%oB1uS>=fhJM<l>V2p$OJ#b;qB%1`V_~zm*nez@_x4-3p4)usPM$>^J2BM4fE6=~SwKB_YQe1BSZBR3qx#_f|B;14LknhWqqJ{dS+$n~jg)t$5rhJ^ZT;k>v&L@?qhiEdynYh+T!Qd+%>u9ugFi0B<vy?eGrtRx?TtLhCSk}6a(lgtrm}=5lOqC0y4DYwIb`v^h#3LYFL5=>w&^0{<FB^U<Beevx*t7yD=`aZY{n2OVo*XJHw3D;~nC0ZQ)14MHfjM>Y6zSa;q?xesY#8uA4s$C)Ix3bg)DveT<WF-C1~=#ZvHh~<>a90S8MbR&SHD~9T9>5+Ztj52)pnTyUgTyob~PkOI&HLX*3WAu?ZFs~u6eu9r|&ugDDV-UTX5)-jDN+mTnBjTWeSPuVWWa9s5Yc`@Z-%33bbj#E^RX*f~+jDb+Y65ISs~LyVN<Ao%xKc+>$}CXK>$s6MJ)c>PTZd^NVMz%M^;R7f*a8%J){U{!$`$2-O~dT?$RX;}F!?^B6FUo>94%J<+OgfxXAg_$kGw040i#tAq?G<fR4vHgw~Q*cc5i4ipAvZzcuPhPB1-tP}x+uJOfW^*lcl8umnaFB4MM1!A2qoukF6KSH`h>Zwkq?Drb78A2nr(?T+QV?3+Ve5}1nczK<Am;cU^rx<O=Pd|%oaw42z8B%HW=6b9x6h0wYXOoactd7rFnkvZz!;G&n2_bw}JlxcP8&y~!Xv1u9BHln!n2DQ8bCj;L1_E%P?G0fnLLL8H&SPqNxhS4ZViqz(FG7^q9}s~8^6Ch6d28V(8sp_nt7bRRB*e*Bil{0#gc;gFZ+#nby2oW9_!`lhRuD-;|Nn#y86x`Z>3~~l?{8u}+NFrmPt;qSFVJRz7T*af$uo<3cr5EKRS9vWt3f2a)2j5@$R(Q8Y)Ijh8`zzb!`8XA<nWIctaO11yhGpG7Re6al#=50D1SZBGQI!IU*16X1f&$h4lbna2ysAuNj<q{p_d};C}BA~ea4;Ta{1<1GqEEZCA0g<;BXI%H0^mX$EuT2YmV|!UlT#O$hf6Cbk0DxzEc^0L5-Huu-PVq3u{r#m7`Pi4fGO8D5*%!%xH%;{{N@wGWh^nM~08%ql+&)e!h}1BX+wGuYqY%RbG=m!bk$>kIP*Wd7pzY<ML)p2rqBaVF=<mfvo9Qp&DkBY)}4QsNhIB;>%v-CzwbzfZP4Ea-3^}<+kJccQ&1|t4izzlY>&M9JRPj;LmdJ0L5-tIhoWD8OB0lF7eK?oxy<-o=PK1eI%AV>j>iPoO>Stz}2G2O=5213fm;IZM$V9fT+}JcX)KHr;lcw)d5BnPbxAZ3Pt#=5tx0w;`oHBT@G&vd8ob@#DmKAdN5(!%ISX-j)}J1Tb4}qo%gV^Rm{W~9+e2u!HTtiH@JLGf}m_uT<TFKr?VH%s*dV~mScwGHBZ^0@w&vf1xy$K!E=?~!h`eR5<Q;?xAg2mTJ*SZqAh24izy4GKA;q?k)Z^!_hh<8ZXIFk(e8t)`zz<a-PXp=d{IV+p}{2<J+BKq9~M7WNpkS8!lbLDZrW=?M^qkhcEB8<Iy7^E=(rUk6TN{f{z>`Mj(}b!(n5Q?)r(;efu8=((h#7AY%FXPY?JRZ$hN3+!W$nHWLUp`=VF->3}+t4?H8GUN=9Tmemj3^=D^^xcPYl?A4H~MM53_TdoGVwP>@=<?>GDAAmxvh<5(;=eOK2G>yBGE`Dw}_8wza<_2t+!aqrKzI!e{<%1+Q{QxJne{3l&!n)VSXvJhS*lyu#Os`BRj8`DyJVDJ$?%`j*iBtFv^x#ui_5kEwr5<b!NLVor3S&t*H>Xas_c?HoG^k~~tl&_bS7h$IX&euw%$PQVmQb;E~a-Y`1Pa4KPatpj`<D`J9{M?z>rO)SUGJsQGy=-{|r&g~1JD)v5HFd3QUwdZ8m&~?Hl@OkMg#<XqnG6?s#9?o{dwv0JR-|#hL}Nw(k~d{rC}OqJjDDkp%j_Vsqt$+h7Fh(YR~sXB&7|pjm_YcJ9vbCbE*F2cWpzWy!6&?kr@;8etC1LKnYwNe^b?+|KYE}AL)h2G%i<ctlS1rhjRFLYCOA=n!(YXx7s9Qx@wR@m#Z3~G=^I1Q>pz|&Guh#!9S3#I<mrco8<xTqi#`J#2c5a-YWo<&V%7^2`5p7;>PEh4CCrl1rxqTMt;~uykzJtO>CW91hLs^Sb@SEcT_Fxi9kR6fPH_E)w%y-ZORs?#0f5F|8cK^Enk7_!Y4gad#53?o0FYIgF{IP`i%vvs>6&Uz)KDIG8LBXTjU-Qds}7m+{R!JY-<4ec7@|I?Av7h%P9qoS3w(i;XpnN6M`Nhn0-7Hq<=`o|NTy34G1jmAt$R?7*U7<#xI5-abO}Nm=z@49U|1t^B;vr%B)7rgb-9nLq^er7hAX@<)VhCA)&*T^DpFmCew7(RyR#;${x<;<#yDH_jjjp_GdN1yT&>AUvld#TLXR$_7=!AIZSAlhB5nL+)RR809BWiyO1>0SRoC{+f>zwB)v0qc2=sM3cC6_b=5)y*fsKcWWu07&35*ts2W6n6@brcF-ipkK-r}g)jf$#|hYnGBOs-xGk*@b*PmqOGB#gF_032a*v-TLpE@eo0**&duSv%O!Nh!fGL8>oTxOrgrC5GGF7tyxBfz`DHK0$3nZqaV;IrTa(Jb&^7L*Q*3>55_?dQSv1ch~H0Gx2*DVd8$qjfH4we+whT>izBRiyM5eA1gYLS`c$c4N)!s{bJ2gv28XulSOl``#K3Ey-I|tNYEH-{fY$Z3)m&MoFt<#^dYHNu#d#l)TDB#a{ui7Qj;F#vXow!+kbgU6YO&aYfC&N1cPe=k<J@-^r_O>@rJe~wVyJ)^qoHlp?bxu=U8HF8Ri<t;KS>jP9BkrWcLBbAJe8=o+&#?);wJ-*u9*O#XK&WBB}p!7~akT3&|-rEhao-Qh%;s9QRSbSQC82LN3baOJ;M(K2;$_m!6#k<3AtAn3Y@Ibb=ekk+f;ZoZ|31L{UVeCnJ-<Y{x<}RX~a>r28Y-1;9W0v?N*FuGXM;y)fzp@5l`I|41qIM1$J2bdSkFx7L6+$7q7$VLxxBEaTqiww2-gb>zHHKt8ng&n3mm$^Cnq?HjWvm;rh}KGW2y-JCGYKIoCC<bk5LBG_|_+22L(^vO_42Qan<={R6K==4nzKtked)>o0OQde(zPH}XtOi`mdKNNCrSkIHzRa{_e#jmH}R*~)+1z^|poa#&C!KtM}vQ;_N{2YOV``0OIy*yw`DhnB<X}zix8`LJ|ffcsFmLA+$oTo#xZ&u4h&#}lQ<QkBcmXAc_EjC-4EM!4QI?(~~w1(VQ%@s{6Z^AB-aj)Q8(t5yZcVYkV>&6((Ro{6EzKD5tUjtk2CZX$~m+w$!5atkkMb8NWO4bj48FvVCAJFsNTgO-rF&`)93pE0#_B+8+(c|?)9EOIw5RQ`FdC0#ejY|u*brnM{$ASuuM&NNDnmf4o{!qipW+kp44+`$lKz$H|D>h6i0%7nsde$~h6U)e{ncX)wCj|GX12$9Crnz2qRl42@fdWvH0<oM7BP>2FB8?lse3|xhwhjp&AZlg&SsZDvGXA*)R$-1Yq>GXkP9476pVXp9RQ?Tr)6Pmc-Y<aUn3dUlzu^e^$wkj{NY<r^Pdn=^FT1F_O+^0-tyT1FGc~=|ZkR3(nv7shigM<;+iT0)=5~8`X46NpHvQAm6gBT9;T|I5D49)(_3lGUhjkOstPE;YWO9OlEba1cZwaJrR^MNdrD!Hr@Ypm(f{5FNS1J@1pn%>!FPMt^C#U2N#&8_PJ33zPk6X*6gcTkQIAB?(LhLyP*1S@=)dI3YFy!zU2Tse+xp#_OiH0?Tuq%|rRb9(j_jigB_+lfXb;|OV1=M22gZPR!4E;XvoxOk7;m3x+ZcGK1=#V1YFQwc>xJSNkzxz<y->jo^ah1la-r*ima@1KGjK}Zlo2g5;+rCWT<zDIiEgTv^#Hak1Db3?Xn(a!**+N&69_WVShTw?|$eRp6`E48hyX2{@3p*m!jS^q`*rAkmiWDDFoE2SlYUScQWO!W(y3#LkSz47Ev(qwzxzi?yo~{6L<%9!raV(X+`hQra(!|r7a^rFl^zWa95R+5#fvx@*${UY<lk82jF5&PYW2|uUZH<KSHm!5{+k<?)<1msbMr7_wyec}<&P3NHtsCl3W>^8G>MwCTml~Qt^qJaB83uXalxLN2B;9+8UxV$!t$kO9=s!%)hS!!w)G8wfXl`rV8KrJZQ&xYgZTG>dLN?_!R|}(M#ulNo8Y<&GTab6Z@(_@CSz?^a?t}RppRcuO%AVs2zLD2_bGEN>m7N&$_qUz^{wLeC#VO6%1RS7^Boj(Zb+4#r0GlcnogEPAWcOtPonD{Ghf59u=*-&9BcVA+9vQ;7P-D2sW2d=bmK2Khl&nSaEjx&sTI_>hJB6{P$hPR^tE{A(xn}{S3nKoW!e%k?P;Tm=s+x+$a#?)w3OyzT8fSpTIfj{%CM%+^$@?{7H?bco&p(@7W-g=Va`NJ^UDCCwDawN%*2=_%R#uNPj$yBdiZb%uH`7bh&&p@do?F6%#j}qm1S%KsXwJKSWjHR+4`Vzb28%^Rbxud8;57#m(7qoN&T<(sNla+J_e*5`Ulx8|&)1OfR3N~ob_dwJUiEfItv{T3V#+uK<~P$m0909C3t)iV!Q5LbMBn1B_06rDI=S-f-i`{d4BI6dkRuUN$E|aDZm;xHZH<);d0mAq7)QVWD6tx;+&Wy26}XZt)npG;D1-C2-lB0E-jyY`nDSEtyd)glPW_uV7cl8>&a{emsF&4o+S9qTFN8hVbAQ5o>|EVlz9dHZChZ{C<%4isVw+X~EuZ1aB!(uAt8L9ZS#K~wO+GC?y-*hZ-?CM|l2Y}RySuMQG6EIxXjy__kECx<@LT0=epMOUx-6(H!%oiRxF8ZH6)0|MopV@3+jCm@<K6OymIaKLl1NB!bajpK=UU=l6lMX#GZgLx-9BU>ioHAW#jFOvLZ}FbWH^cFBHxx{j{IpQp~~z=UpvKwP{RrRV~O=d(^eeCR4Yyh0)%hIc|X}rv!kensT`N7g)Qx4*ao7D`X1F3Pnez+-4HN3nOx%i5St)VjtO_{o!z~bUB;#?_Bdz-d)XYcO4aG&HS18INmh2dNl{r5sX5V${M}0F8(!y|MX0+4hG!C&3e3o$x-nL&o|D-M<!wr<ac4X^sjg|=?|-&J>HBp`1Iji%*9z7Ge-qMS=ZsnE5Y#nIViupHE_Xz9wmwhCX(+Z69Gnhn8)|RTCTH`PDn8UnR`*g8jt12Ru^2eOj*h}v?W}x7HN!iaU3wnUMXC?(3vPz5#S6oZZ0fnx(Y_#U$jEgDQ7<V25Es)?s}+&Cy0r@~vSS|CCFf~&qonXST@)4nt*Ct{nE}YD4dY7)<hvI?Atcs{)OMocnv8&oERKa;EiYZ0vHZZd{bGRFv0cv9sgV=m#9Si2llHZ5Q6Pgin<;1NM@6N_+f5DmTQF~hm7&D%qP-G^rXg--jQ8Vz(jqTOP*AftnBU)Znx<!mxX>hP_KVaYsE%f3S5((Q$BxPK%B$x&_<;gXpf6u6qALhGJ+yhdN;qy6m7CRjTvp|GXXs_fqNj(W(11VPTrxbGXJn(G#P9w(7yjJ1rqCg*7R@fUCfu#!*Q?nKh%R`}z(Dpo@xi7ZJ4#;WG0iXB@4_sO05$)Z*}h?YeJW$)Gk#~2T|F#K+ckF5UW}hU7N!EUxWsb5sQCALrDp|>twTQouE^`@4gAYI4Ru5N<aT>Qza)t`8nmu^uTWaqHL3eOwqCY3iZtM~=ZL#P7PIn)^%uI|8Yw8^k>ggY%54!Po<KasMqIJVG3X>L>1%wr<}QXeBQ_mu$*{252yqF08g7S9g+*4Ef|`oC)h{B-AWcD2h39USUru0aJ<CoJ=Iuy<ms7mv53?<LyHA!0-JJ~!!*B2WzIp!xx;zPni3XdjfY@yuQh^fGrGK9r2lBn^uL1KoI3&O-hw9@j#jd6oK%1csc;y$bk2{R65v$JO<s#ZI=2F>3=mJPFas2m9<g^n{r6)-NLi)zLbnj9MJ%L^%rnxnT=YM<isC=8XM(o93+cf4!V4SDZnBV>y=!4Y-*RaVBfHf3jj;)6{ojOqJw9T4uGn(TSisQL=8KnoW6mzDdC1(>2=M{+UF3h~hQWyEVzJ0RcGa{&6?&j{ufNa6hoHR^p!Rc<LYwF-QZ_s3&-$P)UGqTc~_3qsJW&;`XG}Qlku{`Wvr6(oA26vz7I=ab|D&h-qu=WFGe$;$ZydbGv!0NiVSg!8CRC@CW{^^AU75?1Yq>10*kdR$_D~lA_VhTA?hibI9ncW`3>R6Mh3#%lwqNlF{yG1vW>|FJr1}L6?+oM?%z@xn}8(AyV98D2$#M&+l6P;ZbHSe0`UJ#fyLD`-Qn?1|pZ_v7kUj)2Xq%fUlYgRq0Qqv&(9dy_r<!GA71URLP-S~iTjiiT<P}xb7%#)pew^MME8{DP+;bJR8`nHK7ZV-tK2Iz;TT~@e`K~h=Iq{|FqXdpSbDCW|my(7?+rh$d1AuB?xM_SqXCluSBJx7wsl(*n~04|QV&H41mz!uWw;Nh>QL&L2%b2fm)Yo5Bv62!C`j#{bzVN~>jpF9ysm(7(3g?W$_4_yb^v~kS@2bXO<;$G7u%4OoHc~Y|*z<horfXtxr2N%je7G(dH6v=tjHs}wd8mp=%*p3dtg@G6C@IoQ>#Z$DY3Kw;PM??-)5kGZYpO0{<M$Kg=-fCk_z{CBdR+XJW1WLY3HbqNfP@b_4VtGFVWo{M(L|R8a-traCEhTO{et#lHsh?p5&;S$rOQ#(cu$j;8^*=|fUYMCwO4uj)m*|YpNS6J+PAY%h)N>l~B*N+G5@Ajcwqynkv9s*+FFn^^jvdjtG6aF&N#r=((JP&eVyO8GD)(x12_HB{>QK|64pu+vZuU{C`CrNtk-0~}CYCk7%jiss31rzCxr=RGka6bd7#U-}7N{X8xr)bZSQ1sXO&dSby*DYL&3P>8xDnom9XKB2R@<a~(xPLXCz+d|)wzwU5K%hso}RYoeBr0dsx77rigYuvB15ip^jVF(#RoX>#qXp1y(=m|_!*?stltP_mqyAQ&{S*+{?Wjcz$P?qX2ht9#$z)FI@u4Mjo7|K1d|?Q=UEp&>}I)Lx2#LW;7?C9BYQo)q=e!qvzLLn_4N{~Th2@$1c8%Mp(;j`(Vo%20Uja<&W%7VaX|&FF4f(EOfBV9rq6la^o-4t+Q%Slzqk>a1j^qw6B4u<vKY(1^V223ot;5xJY_G;klu+9A)o3L?s!Jda5x;GnN%Uh<Pwpcp|vc4qgSg5D@Fs+&j3d&B^$}lZ`zY%-zS8#wJTW2&Ko(vp!jJ1wR|c2FVEH*k1^N&-OnZbO~7h+n)UiRPrMySZDTN1m2RF1zv`uA9DgWm3$<L<U$e66Yl|0z0wh@1VeoCyi6(Zbn&_1-QJcy=h{+z$iDYAa${3X(EQt|ri)rbzoY)!&u5Q`?@iS}#dBu+ulJIk%n9#u0rU72f#9Sy^B6c2^xDt%OI4sD!N8ts0CMyctt4^r-uIkTy;r0jg?OeWC=_~q_MzlW{Jvs{+I0?>m9}rC7C@8A!q}V|xlMGU!<Pbb#wu$&SYl=?29S~Vdh#XUk9>bv9#*h+|M0%NP4?BBvfxxVyc?(iqOR0UtmTrCi1zJ`dlTuNV-%Qz73?Qj&QQ!Hc)bJJVOoyJ`QvqN_7{3z6F&PK;%MLVT3p1j-a1$qu3wDS_f5Xn~LoOp^8QtG{!EV8_YRBd}Poqf1HiA?ojBI&4)zAhLML3kJCULJa$!rQmzuINrF)XS%=aQYtkaYy3K$e7_b(wE&RW6_Pt$4M6+tZxSGL{z%`fE_}OE!$p@y((6taEy-kQiV*3D_OB8sABPYp7OrDYJekh}*$VhQI9m3osdFx^7=92|W28jjEbn?+0D?yqYQ#!NM14;PMt3{aH-#?8or>aZi&=`X=z)i>o)LydO(BXKDWWm|}>9s1xh~Zf1+u!vnG6hjGY^-J`%M38xs~6*$(|t4sb3Dc2TS+mjWwLSsIMotxf|5PYZ6n>yEiELJeNEZ8Qx@Q&jRXs-3@jOOJVw%SlT*&QV@l*2=~{d4YuU9b3`8nP*F=!V;M%3KEcV{X&$6>f)S`G&jznu6k=(de%%mBn4EJui&NVJnk|#1yNVUifkSp^<%25y?MhY$?l{%)2rzF(aSHV%+3^VU93=gu3pe+>zPEPxC@Yr5`?TrZMw73+Gu_ct5Ul7dG-3P-tFE`4rNT+R5@LPVqD*n&!TV&}f;KFx<KZca<6bnMM<25@3tMzWwe)s@mp9v;lvBDNfm6J4W~Mx6Z6IMd4t+8D4f>j|K}n8vF!^P(>g2kIez(d;3ZevrR+_1XUvghgpbww!?DPrh_?zo;qW&m*hQLtz#|Ms-_APRn{FFaA+9L!K7!F&~1x*hfZVO-HriZUD%ycL(O2-rk%D&WEj3U##~mJ<o%_0$tMAz%a~ix<bZNC$>i|@%4GrjrfK@#^j|AIPoMEapSy2uSD?cOx7em@flkn#Qid?U2v>@s+L>;2sm^Hq)>@TkWOZ>|&JG@7SHR|!v21bn;hfFf9gI7I5-H8QUGmTHzHQR4ThztOjcE&x1x{-4_S9_7%))`Y$rqUiwE=T)kfAZ3Q@uEE;VGpEfg++{LGX+hp*|?pp!?dIvGX{s(IK3?y(n88j#ruENb+qd<hQN{pdsi$=tCxl(c5y%<zB>H!8WMx#8eTUFMs3Rie?$Q+9RbrU^ljT08Bx7TG&*kNE_`niVsgJV2Z+I)GOsZ?M5TehixY-{wcLQd}Ifbu8dMheM3*;vC^UY0ugB&*6Ir%b5+D?9s2aS3Z0Sddf=J~!o2A-p_tUzSRGrhzGwfm?O-AqVleFnywdbx{s>nC=@^iP!!gk3lWXTzCvh^56chxxZd@`G;RdrsEq>HT7scjEYsiPdVyU`n`{y;BRjti-49-5V1hR;zo(fdzC<HtGs|3GemXGR6*G52;Qf7a)^A%jR$AkLQZfg5=yhGz?)OKa++j)K^FoxUb*{l{on=)!2H#`pGxoM;DNR4EBf^fh}c!R#S_#%krvL;jw<07`kgj74h4BWH0pP2P;;I)qqs3Z6mk<dGgetZ_))x*M}82kshklwaK!2GTrG?S-LG=_J}SE&$^5ASY;rM)YzaDNrW7z+_$!49RGeK5Qq@dY5l=*Q2lF-=s}o`P!}AkmvZ1oKklPefX&p~Bz`x?StiWRDPBD5^sH0p}MKOOLqtiwIW>V)RNMkOPPH!;UP6McGkJyG}NwXYgWB7$8I(7wfV_jA(;-GjDZBXu`<iWuCzFE=wFL+E?otL_pJ%fPnPDnR~M!hh*^=R*jBwm9exOF)*d|lBM^rNJ+Zp;o~U96=qEA>7~b2STD1QF6d6E?k)`)5OfZ+C0HasgiCOq5KEUG*VE0VI|G&BmfV>UHJyCM#R|<H2~+Fz(N{Wo3q)+!AX~F|NhICGUrd$L+rD{PtS#Rld@u{s1Ss+Wx3;31id!L2oLDo6xjAu6o?7Lb^Xym@QO4kuGAw&)-K}BRux_zgP}EfF2tgYl#skk+NhYll8fieE(7ciyKV|4GV(AlEd;^$iPHXrmLw=A0*WnMLD{|F#akNH*30j{P`wP(gvo+gZTe|F#|B0XqOMaDux;7;xJMAD_{Dg}qTeX=Ky^k*a`%&*NC_}i?7oi!9FhxlcWsIKYl!)D<owK48?IA}6J-h3V%h5VWQ!w0l{5t63<6?4B*jQ%Ooh@G!B)9{mNDERd=hws)uKljs?1zTkps#M}J@iIjC?(mhv66ZZQh~-3A#bZ2GL#k}cW9_hfPV|@*+bZ*)1!!k2PZQQBd#Ng&Stkp=*?BNIBHWl=Mv%%E!W<(#IjTI>g1ckkB8=imFA`%;QFIRPp#TacxTi;k9YFgEdhu`>jf6qt#9R`(yqQ`^u~gFDx~elc|!t%oPu?g6#mtu9I(4fj^CP#55Pll8UMoGMyu`fe#@lM_FcQQ-Ek>Dg9$>sH}9ywM0RxjeP`Q*nKzKMgAR~+WB%t@^#>f|J&=s3TM?Jpii!Gw?q|8`!oYC?PwNQfad0710w#UvKa29}ZzdSDmP_{<$h{@(W@8jOG6iZxRo0;$^qk+u>&jc~JL~PV+RAsQb&ZXaqk7t8^j?IASACESZ&BxoM~d^ocTLn>1068N=(Ho0X$;rshqF1<94YJUUJnZ+sqWaafkz)j^6zlyc@u$aeW1@}7m907a+F@nh`CwGPW%$%Dzar%if1}j-brKfCBw9BrkZ(RYzPulc>m%@-NQy$cfq3@2+??e6+VNYG_U}a>_mTI3eMzH*`=7>;`QL?sl4BkdCjLn?wZ`2&~9I~(29<Crl~PJrVBf@Pj?zJ<nkjZnO}8$DEvuIPTWuO|39;sS9H9ZGsq&lalNKr$B0ybz(?kk^j;4|@gF=S=MsGZ)(U}>J>$=wbmz+gXK<JOY$t2?uK0q}2KV}q;>_vpn=axw!t~0gaAc<!|Cj}~0H$OEM!fB7$|mo=g#~`~Y7k!3kdcUM;qNZV(9Iif=uxpEmnX<0|M=Cx!<*6_UrY`Wid20rQR{PRMMe%DsylGpc<q-c%yj|(sqSOjdJ8GUDh$*@YCp6}m8BRznIoYtW)5O~Aqdb>lG{q9m{XDn$%LnpR`9okFjyT@3=;|Vbp1fh)uD6N)^wo@1&`uDDN$6_!tM-aOCV9}m9(CZgXIQVPk+)T45g!CotG47c&`gC!tP8v3f=Xfw?n`Wv7-SJPyv6#A-pJHf$}Ga?7c(*8u%Y2itPV#R_ykc82TBEn#n`wNShLu%QH!??n<cFgvRmkW1s7jd+Ls7xNim0WFat<2$}u;W>(1S(u^}2+n>GqwBRSUVAxWRMqG(bxgy$n)L&|MCj^gqYLVUeOoy@6)FI7NjNSegnvKq|7>rLh(+edBP5$}^DiQb_k?9ljB_*m2s`{Ls)d&vz{uMZK?Ch-W=1f;VIOyt!KAX1%=d#!se9;tSfo~!Yb9P%vAc+;AI>$jwt_x$1hd{J4U5%O~wst!wKKrIvEc{m~a${~A!HZeF5Q&52ZmSA~L;U9NzD~oWls3h1my!vd9&wz{6ppl9Z9zzO`hx$dhH-?dbK&^0DKQ@LsarNSyC2`h>kojd=CNyWhP_eo&-ut>USpRKo+tq0yFcyB`RBbLN!q%OSY-VHdf>!+IdI=uFynk(hU~6L-9i4>@67o*PA~mfNH!~;XV|rbj)IEGglQ4%8`f+vLH*f2Ry-Bn$Qn-WN%|M}chJX;lN#c6Q~%^{wog$&h5{)H*CglEf<_FJ2di+Ev30RqxJVk07R(Aidvj&jrf{2<fm)CQqI<(%0RVsEb+xycb6k$r*F_})j|VtK-nX!U7!MR&DTq}wh+ES$>Drc}oha-&0m`K0t-O%0_3Rs&!nAnoJcj;N(YmNQJLEa39b`wKC$%R0_s6!X&JRl?Y25@Z#FZ~-ScIsMSleR!Ts}}VDjz&g@LEH~6?sB*3NTS^L;MJ_54zRx{wj!jR-s1fj=D}iUxaOpwsusnBbFfRrRyhw8SQEcd|sIHhlC3t3IU;_xM7z%`z<`E791q|J?6&J;82E69aPfI5K9+ME3npSUi8IsRD{p50{kS$5hqKz(Zi=gkUNSaz~SBAsOUsCjN|uRQ;3w4d4$6ii+?r3?_I<4o&q6J0cneurqEPtzt)f+#<Ya}P`lv({?A*1-*uV***&oZC!s{2E&`9jPD!{JRmh|a)i=ANNJUq{;NM%)3bYMAL5TC`7Xqx-k%Tg`_0PY^GbUB!wD3^&+fR?&`!;cu8|YL+C){j8t^vo406zwa-YAS1wV)7zfzrb!&sRLtn9^F!S8@>q91vATkiB{YQ?k2`-fZ#FWO2;Hp211#@+%|{S!?A5Jv`L1@m-@1dA9PEdSdH!K;~*=cGV#MkzixM=wQB2zqMZ}vqpj+s*Rd%L;Hr=Tzu9Jmj3}qvo`ZAJ$6AQ*VDX(gJ;+VQEFv)qB=brU6Pi{X`zuHA6LDY44yvFX2|Aaoy2Q`a}zuB`+4<2p|8LFri28XLs!Tsc%(4ABgSa1m?3wB$0($G+&lHZ{~j$Lmxf?Ucc~<qiyt80{Z=*p9kqHSSJiKwZwUB`KsDQ^{?K}9-kry*v!S-9d+$4vK!kTG((jj+C0KXg70&>z!z<r)Ux&5_@F?5Giw>G&i=1hVt6f~#eA-c8Q*Xq`{7Dw>A58!E=VciP8Z0Dzb5-gXYk4KvTjIHG?R5_BUSG{=_dit;H#c1fiwDcUjH*4rxlNFV6Qn=xjX&mtF8bjs3Vut#^-spo39>|#JEtMXPNX}ZM*RYcZ+0k5zzavMt^oKE&|-=2E`@Nh8<)xn070X>3xJYFm8Q5$*S#H38gOHsM=ET18v|^D7dF=!?&F5y&kRC`+u(DG{P`Q5&^sgrZ?5&3rR1+BQk&5?9hLPI&J7DQ5NAB0-Kd5>(^L~cbQvjLpc~c-hvMet?}VYz<e}B4M$e49EIuJB-xSi`$-F7|Mav^>#$wqIGH;_2Q@1p>{t9Z-7f*%vYOaQuCv8#w^7Lc@?d935iIL$=YEJ~tpyWqbfm-n@%+JfPIF9RYTNPR_Y5b|%)KkdrV@_;V@AwsqV_Ct?7l-TFV(Z-X_FY}nJNkvQN_gj+@iF+~Czr96jSF|NC{IsbMuqdSAi*bq*qprq(?#gi^<NW5e{0SdoZu<w2SkMnOW**CLJ#p2?Ls(iTfquT877evLNXXR1b(7!SoRuM<;*L1H%)o%a(1`Qq2{|nqZsIvmN7srYdiRhSHe&OBi6H0vYR%djld>GL$j`aB<B8`@}bCl@$SN6jmVz{_L;~;G-$LC=K{+wgj<XRiP2|U-Rw$Y&jyOwK=<Vcni#dzxW)U2Ya{`#d=Zl;2Tk$;W2c(?vyr$>zWe4Z^Koslbu4IoX=7CgmDPWpjazp?$G3;{LoO}UrEDjj2<UP|;PG;+%F)<juFfQO+5CyB0pewBdtW!<t$Hh2P89!7IjP+&#F&y)W5CZE2ujS?u1c%JeW9GNct^+5VOnO!Wd6MbD0?~kEcc1t9Gk|aguSJhkRSq(A|9GHyFs)0N?XB+FM&2Xg!B`(PE2(sU~-P;R7ZHxT2w?yz*0sxgaf~L-7X+!z?4!3+B-<5Y>dbr2+aY#1Tgy-3}${em{Y)>$L528>jS)X8zYVV2Ss!@*~g5-x4v88<FmSx!z=sbo}uhrsq}TGo=zBol4;0F?Blo7byf3%I{+a3o?Px7m{nIGGT`zHVHN8V+5(P1&RzHZ$s2ne^whmD3yCcgNCx6)?JrZ)^N3^e&f><G?|gTTxERbYoDy?H`f<A?hr3d{TmQi<p~?Wo2Suo4PSpU?(4#(FITDpt&IN-67ym(V8c7nHH8Zs~(sOgMt-apvdtMFoC4C-S|GhIQ^FSW|q9ckmjtX>o)p8&WZ_Q090h4T5QCTg+^fjG0Je9%2H>LO{Wz3?&Y=lP3nhqgccPsEq^H^#Dp#x7hoZhw<tkVf8+ivQ&&`TE)&N6iHCcroGomDgTe+&c>`pJP!a`6|woVsaEVcH9HQmM&S%d3<^Br2#30W7=)TZhw$;+mmI#cQv7&z5R(V0{KddYyjwT*rb*F8Knj-S_k&NJKQDUJ9t-g@dpl7j{|S%poxa!;umR{0nr&8gY?vbDq)c8VL|M&+$r<hk=%1nj43${KHWlswm+JLT-e@>oV6RvS+=C(Z4$dlZ}#mwuchl$&{LKWJ47?fhDXVU>HEaiuD`4A(ra-LBKr?dv4Z?nXMTCIPVExSYE=f_|>fnj<%`*VvE$=?Afb))DDOzlH&}Tiil~^tYx8;K><DB^EJSjh@Ksf`&EMEvcvsA_xJ*!0aj}`w50eD1KHnbs|Bl6dfe^v7Snn1Kb{gDWx?|ET82eVuD(J*jtCx58VI$aDOEJ{CfEP8*OOuC%3rMJ210c@Q!VMY(dWc)fnr$BOx^{FA@mF|9APtESw^>r3$24n;qi0MgW?$#jdH=t#>!y+nx?}~3QDKJukWurXUYibB|?=7OzJ>jWO*bn_1I%RSn!g#K)y}PauZ4wSbvCDyQyI=x&5~Q4nh&z*`Z$I^&RtTQX~4BVpA1F#NB5iO^Ko?NW-EaE=;|FmfFsm$T36KnK~d#>ljHC$RBE@n4^_+l}uuT6T#gkm?8IOhBlRrC8_mGeWQJXv4_~gXr`Jd!ZoH>Hho)r(Js6Kans>Hi>>Z8Pn8hx>IFQHQs^$(3*tL%$xkIL$)lkm)852>g1X53e_pza(XsalROxf_)wCE5#@su~Yi^k3pv%uw>zk3RgKHXe#d!}KV6e<V#{H<8fvMgGE<qFOk<w+3s#T(%T2OLH+TeiqAt|M^po9j&l7g~sswj!3QKS_F(PXeg#i*=0)XAnA*WsFuk;tEryaE4mzQtVJ8Y%c#c*)j&PpoUzb(gsz(qLgOhl3*qh2~&O>Q-mvyE`pde+;9DyfNmITu<N`>ty8BvLWckU2c|R93Mv%Z4hS|&B3q^XXBPh5Jcq(X@#Hl;ZjaqR);cXyW{$RG~j!dB3Q@PSVImmUFJP4ZpA19MJcN?=_@Kr)`h~CRnIR)WSf^H70Zp(s8qFkvJkRJEPPkbSX<3osHfX+vRdiBau>Q|T}VeZuA+3($M+G4e6yOjg*-h`9{T`_;FSBm$hV}Qn5l{bAzR%@uAQ)S*%AR>jrk4GvE97k6F0XREq8;dv|HMxwOq$naY-x<NPe!QBPv&yqdTiTe*67Uy#Lgo#NmHuLAr*{=q{=p+T6|#%co+X)mdph4~suBst#7^&coToS<n)1?w6vVf?(jXTG_<wPDS%@_@42egAv0*@IQDv{K2w3f@w5{k_B7k{fF|3m%IXFZz|E&ME%Z=?u9}Dq|wgT8&AO!y^XT$s860T4#W}_{_~^&{KT>b)t2FH@E3;z_WE`544(oI9bTSz-x&-kAh)ly5FtcP(R+=B@YhpQIQ81x^4yWvzFZ*u#f0*0W`O~j3FLw`+Dz09psmuJOVOf0ax#KXNrO*k+xuvxP5Yza-26x{M@dq>lxYox8mx2IzI5zIa*J^|)d+pv`6iBn-P@D*xm6GkdN-4_ZTu}}SI6~92FJS%XHY9APgdJFZI+VAuI@@rTCW`Dt3h4`ZLyV{fHNGIrVSV()&W0ESZEJ-X;5XkosA@^3@!DI+ta~#BU`rM><hEG=jQC3-gAd~0st~rG$GK%xiS#=BT)9v;}ZD06g<n-c*byBEa`$bc}laXhh1vl*&0)s-FCZ_owp)S0SM%$sya6NAxWaS2~&rOoF8vLe=yWBT{|unWkTCdk(@?l<Emz@EEY#k*1^VFw8SC!ayqP{`gh-kOp6*|hM?2}n;bM0I%!&hX^Rvk!f1PpD4=az1^Smixe6ZXIm7v7H?_J&g*;W~hTW<iOvCFza1@emo|4jfr=dA}bjmo%hfX7f_M}_{AAm1R#@=B?ltehsCg9O9{tPK4VY~<IsOA_O<rWmf3T6B88gpoa7Kfi5&S_p}U5{MWLWZPD2n~Jd#f1wos>ItUa+z^}O3%LK2yD4u>Yl&4xvJ@zB<y$?8Nm-u@=mbo;(qdz9(?9T=7|QFxo~?2^xn*z)g>BF^Fqctcs=!GT<u#j(OKZ$qkwm|^UCoLXGdp1fSEdwk1GBYmKXj*ONSR#<?G1K2Fw=O>+nKD*7e!R1~w#=H)Ub{iJT1}{v*o%R<70{kJ@yp`Xoip01*x0zRvEqJ|o{m1&|JxYK!esmc?%?sYSo<lL6A|ptf@(QGaW+8MyYia4Sn1w*kreYQ30zGwG$N5`m}~h*GdAAq&}9zZe}wuHX?HGd@JD2!N)}0R+=pH~656F2ejUhWg0tX~MJvWAYq)(c$!fBoHGdSE-jzKnS8FT05NSw?mPaAU{kdrQnE#5)uY5uApWi2Qu#0<2l_cpu9D*zV?%z^HtQe_fb6uOa`|EwGYy(_otNtZ~#+OFw0TB9>B=EXsmJju<(TrT&%5)SHpxqgw7&fWe<7jE97-PnLCy8TBk8=!W9x1O-LBRzP(BP%+gu%!|6Li$W<wE{W+NlQS613mP5HGfoGIEpTe8AMY(uxB22~hVt2XS(Ab@^Pvs2aCTQ&E{|SB?hIDrz2PsZZ5ruHJ55|D1cOxqcnn)3P5zo5-uHb&Rp)f*@uSb7hP4qHsZ@3D_N^W#PIrrEGjTrcW@lt2bfROx8UGQGpDucrL1StjXhcYCU2kD=+t{>+fx47d%Dzx2_`FYakTAuj|3~d~AioTj=8z~YWPNC+lbzsbHc-FdP=y~REB?#I!Cc<In_lbsdp{3WgrITC)e+3GT-htX-H0o$9-;dQO0|{<iC<6^&R4Iu-Rwg2PE3KFY{0!W(xgRxY$Yrsa3RjQam^5S|+3HW|ESRQ8jtZ0dBn@HhFaigJxj7-bW8@(1Hi87My+_~O$;z_FDzY{UWQC+{Z7h#3odJ~WMnp<ns%4Ar&iT79F*ub{TA}Um5kUXxPT^D1yRko@Wu6A3B>QDddLC*~%v`b0?+|UT3i)6G8;RA25nl^AT&D&wK5xk-*Wiw|PR4UDj5AoVoxU(zJKX!b%SM4G*PPPKm!$`Cox{bDzDpi+!u-3NFX^4vX^9T(xw=YYm^KjX!ZB>2-8nu)x?|@-Hoi*|(OJ&?g~y@rExeL|6eBqL_2}iyq<@`0WtGokK<d>+K)OE74};U1#{7@3C^rc@vKx$6x<*Q`hB7Uz{yjI(CV2?Rf~354bnA+B9O<at+1wjY^wp^BTg@`7-M(<d_ta``nqCS-%k4_;6JP3CDPW;$8In)?(Q#B#FX}L6pK(V;>NEl$6~TCCH71k#A&<<L9J^I92jC1un}RmN(VoH(t*)R?4sv9Pn%!3mC<7k*8<CzuNJAii<`?ca8(v^{D~lH|Y8}iDZ*M`n+j~UIAT5f+gP=@jjIjE*8ss<dJj2H1&D|J?t~^aH#z~YJ*!QtIniShB#p(1y--`76^uvpwfjx0ol~%IrOpH4ay;cJ3Mn~M5@0unxx!jpkFJq75|7LJp8@Sp&skI3LBM4V+Hl=JLU^M|BmkkgY79r&Zs%U?NutXl%D^XZ=Dwp`~4zpcr59i?u*i%7;@$b_g1qy&6W3#yyHYSMU<f~Kp5o?Mqo1IUw0g&lqWjYI9ac$d_fG@4CTQl1;=YR=J+($8*Dg8C{aw&5mrK+jl`VVIMq;(DPKGg+;*dd3(r4!Llg<FXnw{Op8*t)=z&W>CuOZqB%an|i1uJZj<b}qK~Y)^C@3$YCrCC%rH-xwRRIFurcp=ck8sr^W})cdt+Ci_`@0$imV?SIbk+eAY+nZqG&fPenBG{kKUqc=J&Oba#q5{%JljaFE5t0$25EmnH1%^4q&tVv&B^#SCH40Q8nTd_o$QDNqtkWy7o!mk9YP3mKh7VQgG2ms2IVxI*!Hwtl}VL;o8^v`B{`zyEvoUM}zXMYwj?=a&ffY<W-UI4e3=6KtC+VJ%-CrIgx@_TqT!qB=@c-t`~kc<X|XJ-}6Pybs{Um*um<wyxP_SJ7F138jDv2H^OwqymM1HhPN<kyi5HtR|^Lgyu$S7fypQ|@$UKM6>>roVHcY+eaWOnGj%%;|!s8Hz^e2BeDOqLuEpN|y{2jFsj>3=a=AfWcPE=VkY}j6_>ZY5948DI;q}f+NbhXx$O>hetBO4{06W1r{OuO5Uqp<x(NGhOE*2S>V};p#fwSTV7g*-k!qU^=Gy_(6O#f{mN7LR*)suI4yH~npD_z!9F7-f7++eK7R|t0bPA-&-spZs$xyFu^tq1p>ou6i-hyZlv=A=T(->6&PMSs>XfU%z=fX$xb8n2Qp;?ugrO2<xW3M0>-Ejog##~HF>M>Gcavr-#p3*->+jK){+LRlwnL%l^^zm@Tu{fHfi&E_=!*M6E7z8;ADgi8)|w_k$^FpcrFX(nQ(!%HEswRQFbX4S$EZ1s5y5T$E979iCE%c;QK2~_s7%an1;;@(c>G*>7qWKSaQYpZE*(7;BP~I<6mX}7+C{^tOB@e~_!OQ|NXL;Q1ILT;NgJ%8+Nj!IsXX!n<iq=I1d!!QX4-nDGy7&?Z3c!n$)#hdqg)Qk<u3hXF4|pl+2?uXDof4U3<t%cHutl%GR!b^?Dfpi#_mB|i58#&OM1^4ad8l~aZY?R*uF{_RP@dvdS$|qlnW<I@Ab^1etutY@Efj9hnQM)PilwVe`k4?>v^l~(J4q}bqyOTcMQ15-_fa)9j37?{4c`9ih$JYS(?#S1RHb?Zt$C>DLuT7ITGeQ!M4U}-zW;E=fSipA5BitA}jG@Qc}~HnW!^Sces4u5){ua{;GK#tgi^pA;)-?1%e8`br1nl=+s)mY}e_Yof`FN+;v|(oUb$+91{p2#(Aqq#~O^1JUwquhLE5Fgt@ZYvJg6yx~Q9^!G&&wD}L|x>lg(?lgi`6<}CJxoN(O^upL|~r_Gd_n_odmmM_9jgiSOyrF=+g_U+uloW&J+BD8Ec0^M4-TLh6-G#(XG<ytPwx65;^6NhPJ2US#-8C;gs)bZ%p_dE)FPL4G<1YV>(<$Fs}1QWjokJ7msV^M9D$}dfYLoQK){x`yTVE7h#IEs%SYSS9MXp`!%9k_Vsj^<2=$B^bD5Qec9#9%6TuH0(TA)uut>e^3qj9zj!8xr`<K_-`TjJd4@PtHwvo|zt_liOzEhI4QYVJ0V)(xYt8Khj?QGU6*i4VyXCwQ%$XJlN+xux`jnpoe<;at5qqc0(jOKtuzY9$uZ8VuE^4HPDMymN$z~x3tIVn&X8}8)?0T&KpFPnM8Qs+=!Sxo!bqh$O=L)JQp1cSG#i}mg+eI3d^HcZ2A(I!<q8RG61)`YbGyS#qWct9cQIG0hNx^*YA`H0vk&j_Ze(c>FQnROm^4C->8<PzTi=nW($sjfEkRtEc+JPl(MYbv~Bvok0uS`6##NE<c*Vhg!f2+errzrf!$SN&z!^+Coj&;i>)eyvRjd2e!8cU%*f?9LFQre0BF@a%}7xS(YsSDOUzH=c-jD_CUv@RCaepQDM0GjfoJ6^bVWobP_Gc;^8<q%YNM`ew^fWf*;6cQVV_k~nPxK2JDibb^AlKI6wOz;(kSUWgupcSP}u<d5#Ina**YxRCL>yqRxi=FqQw0$J8kjEi?p-5nWLh~`bPF*hmH235En1LbKJ$1pKznxbRl|RfdkE>dn}YRi20~YA>mC1)n+FmT><4KI1EOUZb~2#cdjr*-p0u<f@n_1^h9!yJV7Q<gjV8W8Sqa*CX#To<cCD2z<({hn{!Yh2xfJLbIauDx<A`iq2m8T>i3$R;<Ejg#ikvqxcBb+D99~FFKY!~e6ETqA6!Vs8&C!%@7K)MITZ|n2~%U28kuwSzBTh!hieqmB_CS^hga4Xb9rW;ZqNFB>^76I!^lOO@3uoY6Cd$rJ4BWC$Y7cNEZ!9&RQxu0akC_n`fW4!wVQ*eoZj*jeaY3{H1t981)~KFJjn!&^Io<F`8r*OPFACYdlIboaHsK7a|d6ZCr}v3xf<T+)CzBNHG{$|(OB3F$^*vzesvBY7_$m8VGdgs_Q(H;J;sCgp0Gkf^j(tYs$W&=t9{0cUBK5$1#3sgYGUsz>@Na@!RpN~y1AuRkODUHgR)SjxadEpRrO92(OaKcRn5LVR_iOC)k?R4wu-LTTxF^6YmM1&+t$-caGAmp`E49TE#)pmWF2=()?Py-7$A1(($l*vf#H#~@_<x+7h?L<{F)Ojv{|T`1Lo@mFygIx(P)p))Mb{)&frfJVWm?LgLwGuwVX<<>!M;znNegem$OG0nK(f<?&;OC!EoghGg$yuZb8|qtdpd$87?D2NutR|dmG^r!dcg}kT$=|JEPnRf2sy~MLB1zCwr;rceyK@afP}*8C!{JghZ1#>ad?~t0=8?k+-{;Q~{@mnhNpQ8sa4uIqxzmffkax#F8?yyU01u<F)Pixs?^9<pKLun<k`N4}ySRRdtVV_O{47Lvnt51A;J2mFb6I)MO&`N}L;d1Za*f#bdGoHsfL6^pq;Y%Ap|)uB;3!+-|>@w7~D65$+%%;DhtT#(BbWlSU#p3;SAx{S%xmsJg5ObS<ff6He^Wl~6iO<(vb%!+|?GX`AdrG@L`qchZ;Iq<L~AM6%2w;0z3Vq%aV_=1OvkP3K6q3}_xC`0N)^dAsx>hidHzWPFw_s~2&|0-{IM%^bpqIJ#MBaF$<q9L?I)FCW`TQB-p+ILTQ|#{x=!S@Ir{?v>?AO><Ez&J|F1l}oN`w^Ip@r1v|i6&d3^q!D4gX4M;N$;}T9P*U}gLr*dhNsNmcqa+MP?NpgH=*d_|>}+5Z&nL7q#k0D_Ahh}wc!a`kn?!!S9U<PFnSswUEL<*I1B>sJf?mT}$75i;5>zt=M=8m?EeL!a4yn>o28671UW`CkP_D`nb>=iG+>5c3oXkX<4ONdijZKQ#d%LgnfH?1TMQIwmEUM{DB{Io)^4=t4jGQN8hLMI=8U{LsqzYLq9;5mca&@`tsmJU*Jc6RU4_8ZAh6GpGdm6R&RexOK1`B^Mmua1e-DpFqWrLx_&;5;7JHHKP&@RRAV#gCSh4wuC8FDrw;;x6`{~N@lV-{5Smg}qTA)y(XvZDvWhg3&$_$d0yoMQK&{{$gL5hRZt^e9D|*P0{WM#wBzsqVd(_t=-c58_m3xZnt>7*8BdWockk^@8X+RG|1I5qIN;CdXO{Pa>b)(i0-0uuK7KB8$FM9YtJ&Fy~F+P+4Y2?K7HaBh+yez!K^tor;!F@h1|$fCM2IsuL^vs)-YVUX6Y`p!6VJ*J@DEA^Mi%j}PIpRZlv-Kaq1Ec@v+~Niud)3=(!IMR(>#w9<CO#;?IoTWUxm)`+XH2~<DI`oW+`qO7dz1zZwE9#VJ{nOin|R^k5CT2Xg_n||0pc3rNdfY?J1_0$WLoNfzQ<CrF<76Kf7$~c;aikJ%C0t_btr{&m~wdT5VrE;g5Lvfp0v@K2v!(rh(_=%(@jwGf{=pHp(I!d!v0gfcg_(g9UDEj*@3VW<{6hkj%&RV$mV^HfGDP_0T`@c&4SUWEI2inF&$+ds~MKhtO`jR&`YN(?k_kM1P7tDFP!xuC3{X)#?_<oL+Hz6z1-che%jqUR*1<hP%k187)?I?}j`m;nw;=UEwb70y~44*YbYe$I6wyG==ZngSfGJRu?$YfRnYn3g>s(CG~c34SndKXWAM0nl|`Z-Mw^94?&t~hmyzzVvFCNo9e4x07}Q<V`#TonKe`QdboSjQ(BrMUF-$^s_~?3HiZ5N~1)(d1JE;LL0A!g0g#6z4*d5L1(tU|%fLEtprZcPnxVYhC*22t{urn8_VDQcNz4tX-`J&dx-zXc7tLjn9%=uyv<h5s(pR0w?P4d+Lg1cKE}UV3unPk0_s96w2%8X@>wr>kDTWSnI2+L-j&xb=N*jXDlId0(6TF(;rbgpwe5;Na|>oH+z3Jh7Er(1hsi}H@<3@L&>MFk}$ImcGfNHTCCtKYUoqsf3v*)u(f3jDBKOf#CS>)-NB9M{BW;kk1q3sEj*)6xY>*l3scf-z=0~Lq@Y{0hzpTKhcGZMfdTm^Hs@GSMCuE0*X+MFkEzAJB8D7R*b}%%pg?tvFq5XOK4%`M9I~ui)>bHD<WaD?mL~%6N;s}P9Uc|~q!NTbaRy>YVCN$GJwh;#G}smWK?$LTBV+&=+{m-Y&!04GZLU1S&=BNMTjRV`?gPi-+l4=WiQoW4K=GWYCsL7y)Jkk4?!akLWAScE$kST39Z8_sY%V*`b`Q3SU|BS}^d>D4q$hZePmBclygyDzBVYA|RxEYmN58Lt_D-1L7V~NX1KA1*(2MaB7G01+FY1a^JoT7k##jbQuJi6=n!1_Z0MCmk0Lz`1awJAaG)<3{l`@~+EQo_UGA-$j;2L>=WpR*cA;NcTRg)rG$vUY_+6HiznyDnw*WrE`*E=icwLRTCcUrPe*Quf~LAoCodKvI6-HLHr5u$~K5)<N)^Lg=#1EdE{vf91tTzz(w2s)po(H03riEsQ9cS^ec<eJymMLnI5njDSW^c@!nY<k|Mq^F0hhOSr4amlMc(ti$7gL=6aCyOZf9_5-}3}st#8FQqW${xz-dh2aF9(Ydtk_9J^<ULxa_=5!KG1RL~6!)k{%My?VxWOC$BvM8Ngb9Qiz8YzEU$%bd(kQ0Yhw`D=>F}Y+H1NU^IA{L|h~3S8eB4nbu=*>eJ?*hrLo{XDi=AYJyNTO1>@&F?n<Wb<=gM}7K;CJyW{P*79K(H2n{$5?aWe?P{p9+J8&dOk1rE?86}^RerveYb<qzE$(Y>T141W2gV~f|MCXFAU`2@N8st$-_4q{zT2`Bfm7mF0Ejie~AwPRKY?ZL8AV(_1b+6!sVtU(GF(1Gt2Zhdl)P#TXfNIKuBu~J424fqi&2T(P?i0sPOcJk1FV9}G{8I9MDq}2q+hV&iI+U4024x-ZCIWfp{f5&sx+tAJ(mrtGijCp~<#R`F6Oxz-3k$xs3>0%fwnvZ}7oqPD#Fa-WLQ9Yd>x7eLVEj$y<Vu6Kw5E#40=Yeb1QN=4U@Lft`p-eYfFd*|~H(U;5#wm})L@1=!t8pq{OZIqIc^%dBralXnak5)Q%o*~0x3MbNb{y|W9^5XY|ACAOeqC9Fn8SO>j@%0@(Km&AzcXBaokV-FCg_GP$))B}vY&YW%M{jiI`caHz%yHshsPORv3QX-mj$VRpJ6Ckii^kgqL978ti^aU9Fe~g(WE|Yo0#~(bZ@$<7)q-IJddEj;w-Z^<U$>|1ZB~h3|@)W0%*6+!pVlxOaK6S6x~QP05FsWZV`RGo#>416Cq@lmGC@)wH7ry$HjsN4;UxGuqEgrj>qBtFI__Y52BS;i`!rDC4g&o^jn?7O_t3sd}0LqT4gcospM`{i?f6?G@NyAY6mcNei-CnkvAz)^Cw^gA7K@U7R!ie#mK_-jgEW5)GAG9d^LYt=;OGz3K2{nJD7tq#uU^u|3JObrj0}5SRH1~c%EDH_<h!s11)M^iTmD+Rd%zi=Z|RZYh~TGlBzTq)F=DX!UhHnh&n$GP`_W(T`&hg%wEd34VK~LF3_cmuggp1tmM~IDN$W~{|nQeQE;x98$Px-H_ud=C&PvCVzrH;on@L<e<$kzKT0_ttLnVV4};}1)YO_n&-ju`UTiCNLlI8^xa?5-FJaYGKmA#iF4oic%{4O?;Fi7q9ibyL<V#L*?wlC*ZnHQ(GrMzSF03CjotMkkF1p=O1L9*2Pt`?@8R&=&$O}@*6j9rbSXZ_0{f~A)pdA_-G%83+#7tP`s2AxG{t(l2g{gtXWl5vQnC`eks=W!-$4`9d-AHPmLa)<t3pU^u<BpuKY=wXZH)9R=C3xxz`3vI!$eXGKI@0EldRo}ca=~Go!8AAo)~i4L-kofd^35G$UC&myyYjdSa?roD^FN+(6~#r#i1XgZtfFRI@7{&l3p8N>AzBeRsr;&W{8}bJr}09E{{Aro`K8t7Fp&kDm?4nfNnQvIy~a<b%2WxYSjq_Aoz0@8UoJ7Jw$J5_Lb6&~?E&UGlv#!S@m&inzYgyP+@l;gd<qODCujhTUJbY;`5n-tA)5jdGil1n>x<2u{#|(#Xd=TcUG`->zn$NO-Bmu5np6%{0KFL9*tfYWuC#??Nu|P~-53G=0P@=-o%0aTtM3AaO2AwyhgH90XyhzB*|4wBb?e&%zEFcT5kat9?n|&Fz+9u_fX2=w?RBkc`o}}UA7;Rp!DN$P%<|k9-Fo2KnU3np3AHru3-hQgBtL_<bd@kBy!<pwzfLoCS@H6FTl$7px?1pVxXJovdpLEp%N5Y}N2Yk7BUMU=^TQk-dDaMnvPnuEIxu73-tN|cx_t4TL+9_f+Lf7gb?u*9+(;}Bu8x5haz+z<!B~2qOj5J0+%!#I@^4<+Td(rmV=4dPxh`^&wndCe+URe3>biA)R^mNQWhx{=9vXy`l*#iZ^d*sTIoaURVh?Bw!9t*orc-wM8f{_wdzb$K>dmUk9Min^U}_Usp%Ew83YnFSlIAqS<UL=;XHJQ*>$VcJ{dX3rJCd>no+&QlE8|Wh!o|Jb6E*eNvIy5l((_^F<$-K|A{lo%3xC{)R?<5ZF?n-eI*3VsST@2GCK)tWB#Zf4^@Qb_4Ix|3UsN>fCnaRIbODI<?*4UPSP4DMZ7Vw&5pb17>f~;fgMU3Qqi*cAYe$~_;otSwO@oCMJnWw=%3llPP%R$pXvF9@Du*lCLEREH3=mhO-xdDH)CK?>20~33U0&fds~)YRpA$30Rwa_)LagwR;L!S<g?9%}o&u+q@u@_6MSlRt`)Rql#e+1#7f1#X!k{P$wmL6+y+j4bA8@?la&K!N<12Cnsbj9Vt|t0^qwBZjowm!aaSZ`YX%E@{skQwO|84g=%-*2{UKBA+M`6n#4`u@hP-8^yU*;P2UnVZPBz{Lf(49gv)%bt_Gjn*sOy?-@$*NT<&%u@*zMY>fgEh|`gOuXLo^zcV7^cDMICb0a>M5pviAK+PZpip~j)#oO!n>s~hl+h(EqWf{!<3rlhp`_?aia@L3*rybUet$~yvMn`tg-xHfr&w~p8>5KkT(po(-VMT$4wE{#p1cP*rn0)6l?HeB;^1A;t<u^*%oR{00Eb#0r;&48jAP2vBYQl0ssI200dcD'
