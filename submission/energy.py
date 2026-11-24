"""Active -7/20 first-pass, gain join and starting-majority proof."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial, isqrt, prod

THRESHOLD = F(-7,20)

def ck(q,m):
    if not q: raise ValueError(m)

def energy(rr,uu,M,G):
    def ck(x,m):
        if not x:raise ValueError(m)
    rows={''.join(map(str,k)):(F(l,M['S']),F(h,M['S'])) for k,(l,h) in rr.items()}
    U5=F(uu[1],M['S']);t,a=F(41,57),F(28,57);u=-THRESHOLD
    def r(k):
        k=abs(k)
        return F(1) if not k else a*t**(k-1)+t**k/2+(F(k-1,2)*a*a*t**(k-2) if k>=2 else 0)
    def gram(n):return [[r(i-j) for j in range(n)] for i in range(n)]
    def invdet(M):
        n=len(M);A=[list(row)+[F(i==j) for j in range(n)] for i,row in enumerate(M)];det=F(1)
        for k in range(n):
            v=A[k][k];ck(v>0,'SPD');det*=v;A[k]=[x/v for x in A[k]]
            for i in range(n):
                if i!=k:
                    v=A[i][k];A[i]=[x-v*y for x,y in zip(A[i],A[k])]
        K=[row[n:] for row in A]
        ck(all(sum(M[i][k]*K[k][j] for k in range(n))==int(i==j) for i in range(n) for j in range(n)),'inverse')
        return K,det
    
    def cone(n):
        K,_=invdet(gram(n));s=[(-1)**i for i in range(n)];D=sum(map(sum,K))
        A=[[s[i]*s[j]*K[i][j] for j in range(n)] for i in range(n)];b=[s[i]*sum(K[i]) for i in range(n)]
        best=F(0);count=0
        for sign in (-1,1):
            for k in range(1,n+1):
                for face in combinations(range(n),k):
                    C=[[A[i][j] for j in face] for i in face];J,_=invdet(C);v=[sign*b[i] for i in face];z=[sum(J[i][j]*v[j] for j in range(k)) for i in range(k)];count+=1
                    if min(z)<0:continue
                    q=sum(x*y for x,y in zip(z,v));ck(q==sum(z[i]*C[i][j]*z[j] for i in range(k) for j in range(k)),'face')
                    best=max(best,q)
        ck(0<best<=D,'cone bound');return D,best,count
    
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
    # Outward quantize intermediate face computations to keep rational sizes small.
    Q=10**24
    def down(x):return F(x.numerator*Q//x.denominator,Q)
    def up(x):return -down(-x)
    def iv(x):return (down(x),up(x))
    def add(x,y):return down(x[0]+y[0]),up(x[1]+y[1])
    def neg(x):return -x[1],-x[0]
    def mul(x,y):
        a=[v*w for v in x for w in y];return down(min(a)),up(max(a))
    def scale(x,s):return mul(x,iv(s))
    def div(x,y):ck(y[0]>0,'positive denominator');return mul(x,(down(1/y[1]),up(1/y[0])))
    def sqrt(x):
        ck(x>=0,'sqrt domain');k=isqrt(x.numerator*Q**2//x.denominator);l,h=F(k,Q),F(k+1,Q);ck(l*l<=x<h*h,'sqrt bracket');return l,h
    def atan_small(x):
        ck(0<=x<=F(1,2),'atan domain');s=sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(20));return down(s),up(s+x**41/41)
    def atan_point(x):
        if x<0:return neg(atan_point(-x))
        if x>1:return add(scale(pi,F(1,2)),neg(atan_point(1/x)))
        y=div(iv(x),add(iv(F(1)),sqrt(1+x*x)));return down(2*atan_small(y[0])[0]),up(2*atan_small(y[1])[1])
    def face_moment(n,D,beta):
        C=gram(n);s=[(-1)**i for i in range(n)];A=alt(n);deriv=iv(F(0));faces=0
        for i in range(n):
            for j in range(n):
                if i==j:continue
                c=C[i][j];keep=[k for k in range(n) if k not in (i,j)]
                def cond(k,l):return C[k][l]-(C[k][i]*C[i][l]+C[k][j]*C[j][l]-c*(C[k][i]*C[j][l]+C[k][j]*C[i][l]))/(1-c*c)
                prob=iv(F(1))
                for k,l in combinations(keep,2):
                    v=cond(k,l);det=cond(k,k)*cond(l,l)-v*v;ck(det>0,'Schur')
                    arg=div(iv(v),sqrt(det));angle=(atan_point(arg[0])[0],atan_point(arg[1])[1]);prob=add(prob,scale(div(angle,pi),2*s[k]*s[l]))
                if len(keep)==4:
                    block = [[
                        M['ONE'] if k == l else M['div'](
                            M['rat'](cond(k,l)),
                            M['root'](M['rat'](cond(k,k)*cond(l,l)))
                        ) for l in keep
                    ] for k in keep]
                    fourth = M['moment'](block)
                    interval = tuple(F(x,M['S']) for x in fourth)
                    prob = add(prob,scale(interval,prod(s[k] for k in keep)))
                prob=scale(prob,F(1,2**len(keep)));ck(0<prob[0]<=prob[1]<1,'orthant')
                deriv=add(deriv,scale(mul(div(sqrt((1-c)/(1+c)),pi),prob),s[i]*s[j]));faces+=1
        out=add(deriv,scale(A,D));ck(0<out[0]<=out[1]<A[1]*n*beta,'face moment');return out,faces
    
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
    M4,f4=face_moment(4,D4,b4)
    M5,f5=face_moment(5,D5,b5)
    M6,f6=face_moment(6,D6,b6)

    def upper_delta(n, precision_sum, beta, second_moment):
        total = term = u*u/2
        for k in range(1,8):
            term *= u*u*beta*(n+2*k)/((2*k+1)*(2*k+2))
            total += term
        ratio = u*u*beta*(n+16)/(17*18)
        ck(0<ratio<1,'remainder')
        total += term*ratio/(1-ratio)
        exp_bound=ep(precision_sum*u*u/2,6)
        ck(0<exp_bound<1,'exp upper')
        return up((exp_bound-1)*alt(n)[0]+exp_bound*second_moment[1]*total)

    A3_delta=upper_delta(4,D4,b4,M4)
    e5=ep(D5*u*u/2,5)
    ck(0<e5<1,'exp lower')
    P5=alt(5)[1]
    quartic = u**4*M5[0]**2/(24*P5)
    A4_drop=up(P5*(1-e5)-e5*(u*u*M5[0]/2+quartic))
    A5_delta=upper_delta(6,D6,b6,M6)
    c1=U5-sdrop+A3_delta+A4_drop+A5_delta-sum(tails)
    G0=G['compute']();gain=F(G0['gain_lower'],G0['scale'])
    motif=G['edge_motif_gain'](); bonus=F(motif['gain_lower'],motif['scale']);cost=c1-gain-bonus
    ck(U5<=F(73670242301,2**39) and gain>0 and bonus>0 and 0<cost<=F(123271600036,10**12),'active energy bound')
    ck(f4+f5+f6==62 and n4+n5+n6==218,'face accounting')
    return {'energy_upper_1e12': -(-cost.numerator*10**12//cost.denominator), 'first_pass_upper_1e12': -(-c1.numerator*10**12//c1.denominator), 'gain_lower_1e12':gain.numerator*10**12//gain.denominator, 'motif_gain_lower_1e12': bonus.numerator*10**12//bonus.denominator, 'motif_oriented_lower_1e12':motif['oriented_lower']*10**12//motif['scale'], 'motif_cases':motif['cases'], 'U5_upper_1e12':-(-U5.numerator*10**12//U5.denominator)}

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
    # Arg x^2/2 = point K/4096 after 12 squarings. K <= 25/8192.
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
    # G_l=Φ((T-r*l)/s_A - slope*y) monotone decreasing in l,y.
    ans=0;i=0
    while T+F(i,STEP)<5:
        l=T+F(i,STEP);h=min(F(5),l+F(1,STEP));q=pm(l,h);i+=1
        b=(T-r*l)*H*S/z
        ints=[];top=b.numerator*(S*STEP);den=b.denominator*(S*STEP);dist=H*bb*b.denominator
        for j,w in enumerate(weights):ints.append(w*cdf(up(top-dist*j,den))[1])
        ans+=q*up(sum(ints),S)
    return F(up(ans,S)+3*tail,S)
def initial(rows,M,G):
    # Three-sign conditional face expansion for p'(0), p=Pr(word 0110).
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
    out-=y
    if out<=F(1,2):raise ValueError('initial strict majority after motif repair')
    def ceil(x):return up(x.numerator*10**12,x.denominator)
    return {'initial_lower_1e12':-ceil(-out), 'initial_before_motif_lower_1e12':-ceil(-before),'Y_upper_1e12':ceil(y),'cross_pair_upper_1e12':ceil(x),
            '0110_upper_1e12':ceil(penalty),'overlap_lower_1e12':-ceil(-tri),'11011_putback_lower_1e12':-ceil(-extra)}

PACKED_PREFIX_SECOND = (
    '-~#J|Lr@DRut%uOVG{@bFHN|FOvh@s0H3@7K#2^YtJmVVb<q<<RVC0Np&ExBBSKojfu>7dyEB2+Z})Wax5mMz?xIj`X3vhZgoK**J@Yb7NaT!Nru'
    '2(Gyoz4?ClMiBQ5-vbT2p6R)EVtdUmWxHO5_S&2+e}ZlpOn5m%Smbb|6$X#8HHb4P>ZAk+)_<YdFg!&COmO%@a^WJ+&l8cpVw_s|yd4M!sx8oYXG'
    'Xx5<Y7XY;fDrZPtQ0mPKhjOm#$VtItfrLvU%5kdy@YXTt|Bx#X8JJiINOf+7*S)>a~XWhvKty^uHdOe{S<H$fD?Vl#)Zhtzm<a|2^NVFKLXp1L{t'
    'k3*&^vqT+Vr3j%Zw!QKvosH3kbL&Q`*_mO`33P+zX*{!ymbeE`ZoHE||?E<2PKA0PT#jEAKFrOpjk-5_>C(oWx6cu#ef`MKV(pi#TyKi?MJjz{s5'
    'P2TgbrN6da-N?R9vLPjJxdm4TS{hC4xC7H<<gA0!ZKL>a`X%A^76C967JOjZn|ZY-w81QU2)-BhN~6-hPcwv%D?8b@}cxOm+xw8`PCKLUqr}XFB@'
    '?H39Aj8mv!lHFprfJp%dM*+FkvwDkZ#K*FT$ia8Ka=(IPMIX7h!+8vnJTgo@f|pJd%5s3F5o*YO&XeH7<xT!@cjoSurWZ-040IkaFmW3Jn7pSIf@'
    '`y~V#+d4ru=mxdvT2MiU>clOPcJ<)@-{($OsCL)t_?FLAOcdeG%NG>N);9Vq)?C5P_v+BjltVQRk3xmcjOadf6b_9;mU)U^)F5jj$oPx5Pv7OWt^'
    'u^B*#sHdhcIoNWNRl(hGg&@{f6#hXr4-&hf3<v&xWV`r;CRvv2XGue|qTcp0pfP!lE4G9EyoT@IV9pKL8(Px*;Teb+z-Sv<W^r*_uh>?T!Kf01&4'
    '{#e|S)l>i(d;!2QJJBB%14gftKn6wRC9lBP^@$gh-z!(t~sNzh+Zz6;0&m9CyDD`EcDzI=l%_9&M;I#s(HKqLdiD+zCK^GSC9dJ=3?8ZwV@F7tWA'
    'yCC+qb>0M*N)3+Bm?Hd1rG1!Mo>sdaB|Dtpt=P1Y3JPhgU*$^f1AJmF!v9wtW{Mdt^;nA1T=GU@C;F>$eXtcNA`~Jnw#Yi4dfa0$2j~)zEouAU_i'
    '4ok1(FM;{(K@<o>F>$4|Z#O(H;G^Uk!z&ZCpVGkdZaw^&NOk~y`6$a5?GT|?gkZrkeOE@#oIO)dgF+`n<|5{9OcavtBL18K^SW$oYe96#-yKgT?&'
    '!@UrgJe0XsFlTwZvzBEpltDS#>Lt#pa6Hou)=eFvoDzii%s1!<x33i0+p(_tl*96tL-e#O!NYq7>uE#js^9!>R@mYvE{Lc|P$q)BksRL-j_FWRk@'
    '8xPu9v7#HpIo#$rub4<nH!=>q|kylqEb#QOA%4UCGrKr8S@V7ElE2fy4eb>X~JzhipYmS7n?yLKZ8pPJt*yl*i6FE;_EEeq}_wBzM={4B;9yeD!;'
    'LYjhD1JcCgG3G{i``X+8IgIR#0hX|6qL8RRy(?RlI+euV_0^PoE5@ob>iP}d?lQHoP$d~GpA0stzRvb=n?*p4vK6;dY^ivQ~LVa1+>GtCj7{kgrt'
    '|g9gIdB-LqYb004jv(<w9C`^z~8zHOd~~z`b3>~eD3)`auP;qH&W2i15&?!yuivVC^D*&kM-udW+<9k$M~oN^1~Dtz8B!%+AmS_)Mx&ik{xQTiC2'
    'vuC_bsR`|>bxr(KB?j3@!$REjCfF#;+(R)ifBP~Z06V<LOh9|v*2Bb9~oh>EQBh4aop-7pF1F&7Gy`1rQJH%i<c)-tJH<zcf)?}@JcK!;3qD4feL'
    'FAYf&WER(tL&l)K@<Mf;66;LpVT6~2;}sfglZo&3mKO{{w;}(2rk5Pg6S^8AD{Wx-d!N6TH)i@Ol-&3!H`>9d>+*QdaP$!oU=Ej1UB<>raHBB6z-'
    'dG2_}o5?w9j%$YMn2M$CkD15+3JUGSG%-z1&qx9FZ~)RfM)=Kv|xQp3|sb{w=`QK@%8L8(yj11aWR(&=CVv|MbtwZ19zW;wCLPO-f30W>Z%12VTv'
    '`Hm?GIlZPzbAn--$hx+o!axV<@f^)@%COb<=mPBeywC$@r6dBZ39cjK$Pk99|_`u53#||M>av^3X5Y=$m3KQOfD5j6rjf=x}BARawdtyQhp6l<Y0'
    'uaZ8D&wXy1fFX~rPnDgj9QDO1ke*4-bAsTMtzlVJVlv1t;aJ71ruEwfyrW*<oUFP@umk2qMzTR9stTa)!3I=xQ|CXYfI$uiitWv<$|QZE_i^3CA0'
    'BZA-Ik_1FN-s>&AQzsJ~aByXQ}~Jys$$xNFdLkib<Pqq*!(;xX*#6-{Z%{%$<1UHML`y;IicXmGCTm3oimR=2mHa$jUZ>s~|iyG*COguM$|n>%9j'
    'vrzc78g~)#Rfqfs6iZOkH0#O~6U&b}FIiJHQFU-k1~={C`-C|Bc@I3o*CcYF7^vea1x>y(-=?Ov*V8fVd8<T{oqVN$7F5ww!Eh5t#J47q$@Q;#Q3'
    'QRFDyFxb(c@xX?x06IXW#AJ6*H;9NPvfq%Qo9HBT}2b*EQ>LZobsWN>;JQ@D@}k{g8K)+iY2ws!8WPQGq&7&WdK9KPeFjX-J6gYFo>0)hL7QI+<H'
    'JSCro3y5h@EWGmYJ21px164E!C`X+-80Le9)@a{Yj>j7->xP=R-ce*Ed{P6w2F7)pEM5W|SA6RRscO?4a$#-ir&DF!Br-n;uW@Q-60EVLjSAa_}b'
    '{(|4It--sUYK{ZIs4K9=f1d!+_gks)Ge<x9j(52-dj%>Nd<Ny&COaBe-4^;(WF-SFCls->+d@#R}Gq2w)Tr?mT={<$n3$6o_C!_+c0do#Jg)PPQ|'
    '4CZ*_RIj_}TR!)X2>cb@Bli!2+%dP3r@^6<V)aAYlgd;ZuUOV@2^ib#%Scx#+zC`TnA_^AtFOI%o3t5k{QqwD0{hY?PbghLB}7-oxe6dro<nFe+^'
    'bG}j53AOf-Yr<yYB7!k~BcsqAu`sgT0I2?X%T!bmLP^CEUKufCy#cUO3);JeN&FAHOCu54Rcc7FO<ogUXzxa46E<amScb$7Dktk|>ccCOz90-UkV'
    'T`&$8%ddRYn48xxb+E`=H`jE^@2<%e{XxX8J9DC#9F#Gwq>ysd$fZFuwWFx6g{&gwa8!wx8O-9Y~!D5m%Wp4<j)@J|<C&Fbegb&l8&UmgpTJ5=UN'
    'l#uIL|i?MmBs>9;4Fdzw3AzR1WN2lfgLX5-1eD#UlY!NGuIgA%Nf(}%!qEuug$2U~th9YrRU+Lh=a=8IVqL_WpPJ#m@LXE$sUy;XnKXi%Kxdg-?a'
    '!n27@oOY^*J?I!5LyskcPF@kQ1kNd3Zkq7hnVcJkVCV6P0ZT4ZW9v^GGjj0YttMC(LNLCNRlbz{7!J&E<Oi&*aD5-Cs0e@u>wnnQ>y947YPMr;5k'
    'SM*oN;LCQBNDEYqcL7#L3j<dN~n&w4tY;J6`nv9p+F$Ys&>Rd3$B>M>Xe=j?1B`n+bMr{IIfK>ZZf&~;~rmz>D6xK*FUFSoPVsJrbI5hkYngRVEr'
    'Bf);|sD1eo@oO=T^f6x7K|-6}5#8BUfaO^8K0E1w3_!-uCc5&}EE1@7cK4_50qs0q>;-6DCnG62Y+9^yJg&<QUOzOl1U(uA>J*6<mc@nTkeRoDvr'
    '<M>mgiZKLq8>-#;VuCnNh$S+{Jx?dJ&mXxu2w&K&-$&^t4d*1T5pPlFoW!O7Dr%i;7a3MpyQCN5*#2FfG$0#BPC7R%36ZzR4qJOG^7-a|_iVPZwT'
    'm*_Xxq;V}LUF(Eqe@~SkviKT2`KUt8+N|-8C7l56)`^&9LqL%=eu@&%Blmotd3u&EjM=k~>6<Ki!ED*zC3eBTeI>1Y<zBsJG?T&ibUVs!zY)WYL8'
    '7A8}v6>B{SEQ^iMyD&0{IQL29_2R9{3Pp?x8G5kydQh1wR&$UnPI}^|5=a_Ud14B2#2(*!g%~O<S3DE>QPZ1g4cv6z}{8*0B_kjT$QYm_~ALoJhG'
    'R;(sBU0WEA_B<K#sTR*ZM{g?{({`0YT$8Wv)$R(p4(HTq7GIep9MF3pR^));g+kT@ym%qxR$sm)}`3bPDUJk6C`R|==M;A5LlKop2JYn1jWQ~DCF'
    'dOAvPHeE7cYr1_)AJ<iBP^3fGn50z;M~m#z*|Wib+Qs_Ebe`ti0jP5R3XxjoAQtC7Sigg)l?|Oq<yqZr*4BCDSP3rhRP7Cie-B*8)ZoVbIPf>=h)'
    'd!*PQ0bTQ*nev9EMAB7+%fZu(QJMlxn|Cu0O&4nvJ$DfM5D*<I2-(_x5?jeg9p);xEUl!7B{K-CqwBA$9@PJ+p6KTWM3^4t3|oGh4t8sZ!#pQuS*'
    'B)W~?`*8$!#=UX*ONUYj7M!CJ4<E&abbrKNPL;yK&0sEBGrktX}&`XdfPV$dXXhJ2E&Bq9PVIm=Ft;*iS=MPc93j&X<TVDS(^vjyq0#&!`Ph%sL7'
    'GlFWAaM4-vUlC9SW=%{g9x;87R<0IG=LJg>(gJTU}>HO;r_y=?mm?tB*35l%qZ$9;3{?sT%y6x5@B-a;3Srxf#}SfTp#r^to9A2!!%+?32{>|4yg'
    '^N8vJEfahnMv`T{qtfT+8Cxt>>eNZ@eWnq2Vufa|S2@R@gXkXiKY?RE{j`<*6fL98F-<u)Q^UcI7adNonzl~dPtW8bMMun^<+J59_|^@XZyAes)0'
    '>@r<^=z!ernnB4C%n;G&kzYyiz0ncXiWUZmVIg?Im!7*;3x*FGN|76+hvPr;F@Bgd7(i9fGFK42+V|W}@xfY?Il8Z%YGS$E)&e%KdcFJHMcsjr2r'
    'x@WS&a75fpUDRZkiRO8s_HR^Hrl!cXTkgA*NK~+}_UAS+3`_^=kqST`N{qr%Zb9lCEtg7;L~JiBL={cjwozUt0LFP)^IOOj7C>1x&LQ=#1lx{-<m'
    'PJiv|mHc{~^yQS=@Y@h9Tl9kDJU4fO{fLnP>I}~sabUB79t-daes|%BfCQ2>K_JM5>Mt@u<ja!MMohMF=OqR+fYZP#9pI1>aXH?y{c)~}IH|BH%F'
    's(7v<QxDGOk}<#H0hSn=9FF_v~?nCJrh}jL>T1Gzv4A$EUzI)Gk1QdHWgyXY{v63DXn4$)2Jr6aA&NdPGy+@30`>eu>N#=iO<U)k1;8Rz~uQyuv-'
    '*DD}Q3#r=R!jAYLY6X_D5GM{VvgCvM)cCU5xEqbA()g(dPT4A+q=<`H{C@o31s$kI$%f9gpb+z;f3yYLeRN4z=z7~jrI0#fU;R<qSOA@(SPgi=sM'
    'l9RGxAAl=@<XTJXcfred4)pqj!i|F}whs{?B8%(j<D79_1KrhkPLgeF#&Od8LKa<P*m%b}k(ImfW0-{yz^vPEXM2pg`RD)8l$y<+*n*yDY5%xUVt'
    'W_9V9MC!Y_e%4Z<1-dF6Q@Ek*Y+9{;3PhHkR^eZ`t_fF&`L{`C;MQ2OT4;Y^r~nc~ZP62%pD}!ll{5(wWTeZn?Ka!rfC?S4A;I7Yp|bk(Rf}&tk8'
    '3M#p;spg5X#HW6(=_Z*SqFZG}yR5IroEMj!UR42l1l`+Wgeq9|=BZC7bPeK0g;Vg7rxqEbdQj!ZSX?#u3@QWVX&Snryx+#7RSXZ7T;@&XiPxpyl-'
    'Gir<>3g7zZy}je>|crzD$QH)TWlUW&!qnOnmK_Pd4p0SMLBLeutdD0?)&4=G0{lAtJYgf4-Pa-f$ilMy6VO+^x7SiiEfNNIEV)I{wxohxJv3kz3O'
    'XnEn$4A3C1$rwffA&0|dE4eLy-Kmya&lt7%5DZdyiXc(LA>!s@`TRM$+qMz3YUe;mpQ#7^b#(^>c<M)yM#93TUPE)A%PwSEK2(!o9Av~uRY4t5YT'
    '$ILHwHP9=VfUGA^p9lOQibL~OAPgTHPW6C3hX%l7APm4C4trNM&BgfAOXqAJ&=Q87%2H1~1L06(p}tVg1Z6{3Z>sr!yT2|Kjq$1ocJ;NZ)cVsPNe'
    'J=EKvKO6Q8#v)@B|l^uZ-Ww*T#|8<|h<os7~vQcEaI;z?Ju|f|L$kmc(7*P^U3QLoFI4VT(93;cEgJKqT?|+@IQ%6E1lD3HC#QW!Rh^6&^f{(Mv}'
    'v@9Zvt567$(<MRikJ+(e0_2WF7_?dZ)+5uSGMsvJTu<xlYW}W|T7IuX?LGOF`r>H$oxugDsADlRUY@w?3wjys8a$Mb7JW=4+(n(#gq^A*k+PS9hI'
    'nweu*B1A30C0(RZq<^c;L(uH4|VhA`}GobOr#7IV`*fRb}u-})NZ2k8^CSy2Fq_Fhb*z=bVla98!rHexv#-nNi<cf7291PkVPcFBXc{$zCTsE>WP'
    '8;RH$ltRMLwSXOT?wtb`*IW1gS~8G5(wtlIs|AUtw31oH)yk2fAYIjDj*#Hk1lzOMACeHW(MOq?mfk%5_`gsm}g2Fy`^C_Cibo3N6f;Kl~PBbv6|'
    'f96ChPK0Wc2B$O~3NzfV1V+PXj+UV_%Mr6}5U4V{EF6|}UOuNB&E7bd>WO0Ltr4K|5N8p5`JEz`Bt$AxLyFfReVFkj{<wULgfS+jRPta8>!g-(-D'
    'E>!dB4JL&+MbwJ{rrT9V#p{Tu^F!WT>bum7QN!SdTWc|IMawU{&0$KSqH}pr2p>QH5$scoIsHq*&|Kk^1~o4YTF-aKoLMK8*jDGSBP(Cbks-XKwU'
    'kd&Uix2NiD(1yqo%eL9IYKUMWu{O9wr<n7oOl*?<5Eg$r(sL~rT{Voz+*)%F3m-r|?rdI-l7olK=eahdDzAWn-xP(){q00~j$AF|5LqY^9Q+%eby'
    'B7blY1dzy?oZHjiCHJODBDtlB$Z=QQL1bz7hWt!xs{>oXZuM7u9CxnogH<<A6e_?VFtIGoBZf1wGUW&o#;bAq(&TBADxg<=@?l}=4Zp?j6h~*P@;'
    'mVP?)dsqxOqb)7LXBxzb~`nqcRU7Kx;7Utdt2+T51$OUnyw8k6WTSthO}neV-_;+J28)x*<pcCtF{XCijKgo2pI!$OPNN^CG3E4VIr_{owVM^-_e'
    'R-oLtIZNkQ{;gDO2`|hWdB@NFD?dPlRW?}IMyn+{eSV38f&DF)tf6CK?-6!c=oLfm;s8CTJy>ms)(p_U;IA6xmrVU!J9r`p#Zlp5DX2bhOiD#smL'
    'k~RUHedT_n0kwV#90GPH+@R4K8gTWR6i$nFtUx8<3o+P4(y!*ULF%L{sIvCiE6J2nB^;bkP|rwmH%&cD?etSie2mbtD3}akVNq9<Wbj8VP|_F%8R'
    '<M83&rlxxV|<bwS`8$Jd?{`hBsoi?5SaWiur$khA^?q#M1^*Mg@C=2nMfkC#J(hLmDc&Cle0~`6CB4vM}iX-IJ?e`1=E(W8%>Pzvo<f0P>c#=!;l'
    '^;FBTN%}e+Z9Ie2bZ<jd?b7Eo2A#TiMa6;<=@yq<S2V`(io+#g}aavWw`#D-I?bq)Y_M#O4DTarY$2~%)3a%x->b)y8SX2F7(fk-#VvAK?_|9R>V'
    'acEFgE^qp#EjrMu{FelE42{<$qQ2K>8ZXbQi?*@ilcM^VVObJl=_<6AOM(`}fU5lb=II@jQ=xmYhoVnNYK?VnaPh;I;c47*fH$VJ8zSGGDkTfQO~'
    '#9f_^?$F7WJ08zl&k3i8&QyG!#HwFg3E}M;-RCn*bDwC=^ZobAzk&vufqTsQ5}zzcnc!L{Vp6gt$S2+jTWL0fWO?;d&0#AJz=%M{<R#NOT>`?!(G'
    'J_LSa^}Afo*T`hdTp{+oIYZwKsHN37IBZT}?K#&?wHeJLu;)rg4A+Sc9Pdnd0awM++W@le<amErPEPGFDBH<(D*S^^_vC#wtiu&11lAFM6j=sW2Y'
    'eD^NS#Su-N{r@X21e{M@b+5IMk4!Nm@0auzqbh(aDVlU=YHBt3r1Ph5(O)jXhG&ga9oJSb&V2_0}c>pDxiq`!K07L=csv7VgO5_L75H3|B{J}n0H'
    '%u-oRkYN5^N(v%xcVQ|XSy2dh|uJeN*n}1$@eyYPZJY~;`)!0lu>L9)J&+Sj*iAI?^KMo41mnmGNnEFD2=fseH?%jL8!5$rF3THmOu9cNrcv>knL'
    '*I_h}auD+J!<yYf&mTJPJ39j0;WiQQN&5vm7e`PZ;1AZ1G=mrX6U%@lmM&6BG`obNECvggP)YAVHwHq=h7*rsU;7167rkDVQARIb^;*ayu_Mq`w!'
    'vDZ;S56U|n4F_@9YiaH-9>KGawxb^`m{JIu#^cJO1HlAtGt`C)IGu7OiZ5d+P99zaNq&j^VcoKy^M-wd?1y0i)9O^l>x*gh_CD1OdLm;6a!${a9x'
    'F=tKJB0+y5h7io~i0-ECvSN7`E);Ug<7i&lhLCU36}kXMAQfGb1>BYzWfnE{%iw|Lxqy5AYCQeSUjJ;*{{N*3!NvO)L~cR?Jxy$1E$b`MzVdmK;G'
    '^<fJ2h38_f=eOCB0{x}YsgTOItd?L(Q2MAj+gIRPKC7U9C3gv^CaKzhf@acHEmqG^Zs&S~_`SNODMF4YFF|EO2hD@4hdM2nH%wvvHJL{45uOU;q%'
    'zRA4MI6q80XPqjzD&B;R=UGL>x+e0R2m@-z!yK_2KX**YE1wwFYc$_ufSyK!F+DkeaGhQ<lom{8KmOOT&u=}uP^Qa!nuD{q1a`OwAO;@R1OyjGrF'
    'g2c$`Hvesseg?EjnEU&SIk&76a!^Xp}dp7#j#sN?u|-ZF-Xd~Vi*VQ9|XYzT9$;#Ox}r`&r%Bp`E&I5#Z}iPctXgiB%h67_%{T~Zvk=Lq&$%VW}#'
    'xW0#?7J~XomR_ubAx3x?@aa>taQM!4ZUqt$?h2rP)L3THdAb74<7?Lp*x#xD1XkL%f4HCQuRZbJ%>u96){mExFrk>A;+XZw-du{shG-^S$zZ4(a('
    '7xW9t;#X9TQYgyR>B>(MGg;kI8$7D;LztA2f=R2sXm^g8>iQcxW$bJBc9x1D!uuJoM#CXjpOKh1J>Yu`Sl`b!;W1G6R5fH{&m0S-6(Z$UKBs4WsA'
    'UW&~9IW=Z<BwT8vnMS&U^*;cr?*w-4YwGYDYUS|$xi@JqlTxq-9jNZL%+n|au3(Y)X15xZFWvju1)DC&YbrxkjVV4TQD-zi<o53?M@FIx{W)<Kjg'
    'Yla8*Cd9D?Y&$oj4L(1IKwPcoeg+j;*}W|EH|af$R@(0BL^FaJG&HW)~9dK2E7Oyn7XVMx92U#YF16n4eA7T4bBi$*A5^9UV=;B{=#F6KtREm^_$'
    '{clmR^aLI#ziiq>QlHhYkq(?v{&sMT;Q#oY90+*ZUV!*9afrJVu>^gfy7aA_WbE#u^VGfaqfL}v+_!?9Ogb@Jytfpk5@WlL&8Zxy)0T^nIeDzR*W'
    ')0?4mKne%j^-gP2y*I`(y6J^eq(7I^X3TEp9(7cIs*ovyT#m`8xa^e<f0&+VH@~lC%as;AI{8&ImpBJSl<ee}8H7MRHE62`FJ4o))B_H8%rw*{Sp'
    'KSM{{aHAfrOF9?H#$)-+V~7-YnaTKitLO9;WxT)%_+=oHG-c85=gllqLG>zd1KA4ml0_xoOWPr%$S0yjVp}c~2j|ntaV)p6iT9`wJ&={0jlR-^GX'
    'pQSy%m#MON(r6VyDJ2{FZNve@AL7=qY>l^Td$~S-EBahm#+1RHU-3okeeP^$Z)SA$MQ#TAJU}w=8ZKF?xJCklJK<>pwj2ER{c@YzMFL?11EQ&{$v'
    's&P(bjH;;P(L6Noz_K-0fg>vv`+^AM9)to1>*OyyrP4!(9f6X`fpgeLCr!8GyN>le$2R%D<KB8_@SykK@(m>=03OW8Ve`<RjM!K!*49j1VEvj;bh'
    '`v(+Ba*e?xL-Uu*oh3kZm|+-TND?Z4qZ_svY(>9GTxC1Y8Aqc%Z{_)=_a1sHz;^Cr6FoAsC(_&Pi2DUkFA?qNqGh6M^Szp+kS;R`xN6Ef+UXr=`1'
    'PhSX)=pcdNzixqk6)@A+*_?oDG!spcur{*+_M&yHf6+7v#v~S~1F6*MP{mQN5lC35OO0}G66fpRN9iBT?dO!qBf^6z_YXQ6OI@+D6G>MEIpw*E5C'
    'Z?-b!@Q4-|hX#d(fnu7Ucs8x~<=?JJ_G;4%_N)?j1!E2pibT;%Rf3*yX^UImX&!SEf)@lJ}1}g3|ztir&;V248=(r?fso`t}H|*Z+R)hxPl=CX^?'
    'EKeY7MG+#@T@kFXv&^<-Hl>`^SS<L1=VAc`ach42yOa-3N*5NN(btWe01dgXMcC@b*B2YrV%Za+&ZYRak<a(Rjk|XGrJ#No9Z49~x9~*jkZ~&oN^'
    '-$Eq6cun}rdSx}vZ=av7$YZl%U4x^)e4zp<|7epMKYhydvr7Y<YK-Mg*uq)Qn-(E0wl)QH@8Ao;%NFDrK-OnRjTRc5KeL3ZL9!Gh)Vao3VADfYsF'
    '~ycgIDcAAjE*Xk#hmWF6X!DIqw~+?<&PfudBtB|rgJ?xM5{r$6DbVHff|I*}roiIBhM3}wh^R*g#z@6HW?ZM&aSdWpT#{s&Fm;}09@OTh~ajOSaK'
    'vmg^=SQ)FW9BW|6J@Q!faES!6w#B@I0r&<vBxwan^Tye+NW~H6?bJvl@*RGuG~=_VW?M*X%35_I40L4bWnsGHE`+|i;=OoBR6^6&Es!)xur%7zo0'
    'UcQ&=ywXX_m5Fhx@c9BZ=9m-qw+4qR!*<^9a}0sj`aK4t8I#Jfb)fV^@k<!wQzAb{Vu-XN7Uh{GLwSfNs!4<SD6yng+cPOnh-?hEl;Mz1W<sF(kY'
    ')Dl@rU1FIX^MsL}tNmfj5X;ZSI6-UaBCZev9?AwL*ztIuWC$$Wn&<0r|p~7qqGaHW011cG4tGmAAM|69VCW$2g?P|ES&INp-jBwtsX;IPLS3fBQr'
    '(E5nqHy$x_9JVzU`SV+Uf|?FJWb(RaJc&NTgy`N#VoiZDw8{V((3kHuOQjPTnNi?KkQI4CzTR90#`~YGb99x6*u_kU8y=-McuvkTK~hlSWcRlv%t'
    'bRiLH99RIx=GNR<bAJ!G|j^77PX#N|vTSN_tSwjMC1KQMLAoF=2#QOXL}`3j5QYG0TN7CRP)ro9=GEqsiP#9TxOyus(WbXp4mlFTKVNFIQsSkokM'
    'U-=n*Rtmcp#4EWFlzR0@qAPAbUMf4Tb}#hoQ17U?XZ+2&ShF2ObpT){q2V_pbRm1g;@b*F3E-RnRV%Dz6ahjQr|%SuLQMHh!I{~MCp;QnpWty94z'
    'Ap@=*r|6J6-QUt`{9W@N!2)ed{jSU}@M5gF>iXSqxF9Nq*gH?7%eZ(QFqqzN0r4y9OYmy}^wUv5B>b{7rrZ3w<hcm>5Lt3x+2=z>PUO?u^_OrsET'
    '8Az*r5&GQLvfB9&2JA_qfB*sH_H4!cgr~y38kmpbjUJh(|PnI-TyYQQ$A7G<qC8%C5h~N3Al~PDmGrQo&<ZUdgilNX~0yAQtAEUcx8FPF(5>D;jV'
    '7e#skq6$#R1~T%n5Q6(^v$~0JfP98@(W8L#TdR&pBT2wxhgAsJv>tw^H>$c<1t3yr`EI=fl(I1s|DE8P3#$dmhCg!toF-y4@%IMMi0ki$aI!jfJx'
    'k@)MeK*_EtkDA~}N;SG@6dixl|Id=*ybS?E@*`_0BE#kOk~7Cfr=q4lQ{K-8}9?2&=&3Mvb|F-Di+Dwj7=O}s+=uG>?W@eNDz)I8X7Yi(^ewv}MW'
    'vY)r$Sn1zVr~+q6N>sVlD<qqCIWhOA<SZ=43|3J4_9m2!t-(gW8(LG{M(^v{b*NTb2cuk=FVf~*fdAglAJc1ka1jfzKh(OEgZ`rfFz(;+i1i!2m6'
    'zDM?A>1`v8z2V_yHF`b6hUB5v_zIpr^byF5Zf|H?-v4p+ZjbkZN5-w>tC4x7VidqNj+W1C0C*ieRhTA@3_<*UVs`w1<w|9T>Ys{y`?jpAYwQ;^sG'
    '5N<{b7TSpOdeBfWk+6fMEev>8F`agU*0jJMwx0g5kt5BtF18p~)(J8e0)(E`zua#}Io5AIj`NG!63qk~Ue7mMcAzT-`ac|D8UjrlcwwQCC89IqHE'
    'io3(4D>LSq71b~OBo46!^eUwu+=dwk4SA(_tWf2$$A{^RE%<HKi@Oj8L#@XWZS*(C3m{IXevUmM%NV42#jhkYy$02%hL{4L`m~p3wJJF-k-Yl=oX'
    '<gPotoQf^4g$^It2c=W&U6pV(f`17^O_X}!!mW%kntm_v!Z$ggMurSoFP#F(*)Q2kx;O*68KG8-(5lXEqa8`?1WWmY(oa9P>$3RCD;c`Z`VkL(<C'
    '==1*xqvoIgW3&Q!m_h}B+7O^gF$e%Th_2_L-MZ&N>3|riDR<ZW`4D7=bNj3g34R`+xk?iKMA~^^t}9k?N*v1n-+v%$lE<&yD9=S0_N6m<u+vWP2D'
    'KjgjzM)rACm!Ah_H7h;pBGCCui$je7XrL&#!x8*P?{bCEdPGd;X{@*ect#4&u0jz5NXoCn}u@X}R*kU2dQZp`|OB5~s*X4$8$(G>A$OB4fRkrId('
    '|?tMTLcZ|hU6tq^@{b;V~M59hM*GL(pi2<L{$eI^a#=!D~JXT|)5M_6~@3L%KeDQh3S<i@#HU($>A3KKBppmf%S!)=k5oDRM?uokkC&u4!rIr}rR'
    'bLg_#h9CK4e?IgDoKF3_y>TRv7X^hBMS8yo;AfDWXAb^XB4eV(@{}T)O{clRRAY<NC&;($ol6GYHg=nq>S&}nG#5wrYc#6cA~6=!R$hX_QTS?NVO'
    'XExIRYZ9gfAOJ{tYMy|3>*F^I2ZMt?=koFV&9bFyfr0h7WGo16O1LsW-hrp1~<7H!jj3Q&tLQ_3F0Hn!6iC9Liv>tN=?&T<>V(QdRD=;WxTO5XvB'
    'p~*F?ux90{<0xQrktYAp0>wc#t0}FCDizwnOlFW;@%yyK<jabOL=|H5qMA4G*sSM5Tr<IwJ$@Nt18aW3PqMLs{p;$_4Uemn$F#ZzflsF%-v5sEGf'
    'xyxg+7V_y1waYKbZg=-8`ET6{PbS6wC;Zh^p_oEnySIE@bg=o0eCJWt>h=A=Zu#Acbo`G<rpvEKnRX6!=e+Ks>MolXc44R(0AUc*Z~@T_B)4mr;r'
    '_jJ3$%Y5~EP+PaGB?2<|!L2ZE}?SbgDq}C&8R0M^uhfc#!!X2`4m%`cdZjfq4$34XDA;C6s?*x)<z88h_XG)X6s8<ZUp-p>2R_fU`%S*wCcT~sTP'
    '$_0_p>C}z7JFpcV;*I-2PPP@MwGN>k$o(8W{^=5hc8_puBa`KTpci1G2EQ<*{K1`YA=AT=BfvAf;k?Aw*{Gl!8P0klU0qCzz=bwp#EajxswU#-ob'
    'H474S+j*h~&N35v%{en?XZLtIg6<N1%qyGt$mIa9UA;zET+s8vw82bhvRh+6#+V<r*?)%r`1K={a+cBv0#bG_^`+c}NPc}>IzQ##*ow@vbS;tlK2'
    '!l#qY_$3_L+3>@*?@)6+&i-*1?ZYT7^nc3&12o6kXsZBcbFJU3=z%m%VEhd1tXg&D#gE9Ff4xUJujOl~$2(<t{?_2O7O4VhqDcq_Mq*_h7Ms!Jc+'
    'M?+B=Mk_fsmqY`6A&hZ={&@Hl8g*15ZKA0w?UpiGhIRruR2u8aGbQP;NiYeRp;@lb50+>2?|KEWd!DANO3r;9XJ#^f1XiTL|W&`vQuvK)|amDTmQ'
    'yE{j96Z>M*juujUtR-V`UKLw8U#3l}!)8v6)l$XATS;vVT+LHk;05H502|6Xz*DHVz{&Ros#bV9RQM-V6x$CLQ^pIm8$1^6wN=LH9kt@C5Q|c#kK'
    'LXaVi1pKDcFgH(N1lMG{6J$&XeeZQEN5dD;xR(c0=wH3DpADjm8~;`F%=}N$LLhxjCn%?GGZ&eqA-W?<O=+y@P>zJx7F!2Mr3X$cD+*Z!hOe^QCe'
    'Qd3v>VxYUk1Kk-TAGWoZjc_7dkmg-`?3!-y=uQ=LN1i3wE!3uTrl$mcsye@<uk0V1Nq-ed+fpZT|cQ$?`{-akwzCKf7A#e7P1%Oj>A!^tG#aKY`Y'
    '1R0W-d6D*yXZ#{`8!#}0`-aLxfLXB}B+0CNbf=$DyX9j>fe4{6LN9MpsZnK+phD0P7Bmrp72Fw#zM(a5!l&xY?3CXrXEA9=JT<ZL8r9%WbU2$uP2'
    'YwMB0f-+IeueB|Gs3~k&emOz5?R20?Eq-jIdD=%hPW__rYQWg;A1ml33lCAFn<F|KfuM<s^|g;1A^yL!Lf>Y2X$r^7IDl;EYrz8PaE1hQU}SGp45'
    'ud^Qd`LJ-+}IjvH6EdJr!O(5DQ(<TJ7&W+sLI2u{}5F;AkilkgVw|>Y<o`7ze;yM$P#PFqkKk!Rb%OzVncG&y4oD;R{7nQS#m)ej#*(|7tH^}b-l'
    'T+0A1h}MVSpgX;#lo&-{|I7rKrt)p&#^BlU{vLT>|c6sP#t&m+ShTWkLp9hK6?Ri&CRWs+r8nR*W4}tm3Ic-ZSZlp9?YX50pUp~P0&N2Xss5W=QS'
    'Bo(3QA0__b$rm(hPv3cGOyY!WxazK9Tfj=7Jpj}q-Uc~WFN*Skz+gTo_}DAJ~K!|ostv|z8TxnBjCmfTjCY(Or`_7m{?8cjVp$%pLt5OJFbTC$}?'
    'UoX~5_Du^z)8%^Pc9n0Y&A(z4+JyG&=d9l9pU|-mI9o$nz|2XPMaHD9;kUF{IpC%Vp%H1nFo6>ArP~$fWw0Ksp!;Ppts2BWC|V8Sc^XQU`QKwJl8'
    'w&?96aS{w+|3cEms~BYV71A?gJ@8&M5W7S?ERN<aAmYz;v+lwR|&C@C?k3$L=b}F{josA<ZYWm31^WJFcs^DBtonBNO5ap<NZ~xwLGqfyi&$1VZT'
    'vmGi<&O}5il)U*Q9mKa%ONiTVSwrg`s)f=Y#?xOQB@O`{N$=Jrm*Naj=_S)w^QPfg!6-JMwpx@t`K(zM#(06*<1KdDGi`%b6sOj)|fZ((>ve0V*4'
    'z;C56!6+kX*k)!+=-65F<$s2(aj`!-g?>WLN8TdwvtT|Fq4x&d3gK_OP%!QxMe{BgR~!#pu9xG95D_1Z`ln=K!53A!1NT@@rs0sx04A7dGm+7Mxf'
    'UuD@cq);bI&Rt$%5(zV#PKVjV9c#4Z*W-&98ibA-IgcjT9^0b+1X?14{ZwezSm*r|(5%Y7jD(Dk2b0P=LRvh+<1H}n$U@^=O`92bio<4DEEC1)z{'
    'nZe=2^C#$3ISVDz829Dd%OF(V=)a5aF~>PCu{a2iH<~f8FYkn`v#22K6qZuJnx~NUru4Ietr1)p0e-(@d{*mgLCDVhx$0dlvAHQ6nk973TPFx{py'
    'o}WTA#gF0;u>v_otm91PgL)kS$YXTUJF{BLTc@<U2Qv=pNi5`FjLV7tnFt2t0n0s06!KBOm}Iq5S9nyq1&Q2A!U+g~k;!@NljKG$Q|Ty6%4ouk(M'
    'g1Ohv~d<YrGmlCp@3I}!9zMWp@+9Tq*U?+yoAW)xrpK12l(2QssnnFcM-(FK4=b8)Poolk&9-^w<-rwWLx!0?%MyA%3HeDtZ)2>^H3A++bdE=tu{'
    'YNTAzQ>%)(#xc784~B0VQvL!lSjlMz{4--D-3}=N?OK^-3c0n395hemU)A~gOvRclV(j{UW!~_tn+h@WB?fmu;G4+78f}D9@fF?So;?8K1^(|o&}'
    '~i<?iysQV%Umg&iljQvSAx>>Ev7*^~0sn81}K&+pm@M%IesohftLf-QkhE4Q)#VW2N#XrG>mv)cJ6a<kBZjO|^&rxgKdH-#@#MXZ(j`foWrjT^6`'
    'EW%^fE*<R{q^C!X)3KIJ*JE<~^rv*d@+js)`mul*p_8>%<=~Yya-$1EGVR^J@YYRx>3cuYsL6D^O<6KmSursC^%=xhqBuns`M4=zE|q=Jw;07zwe'
    'fDb8Uk?*JX)$!<pa#~HTq(Aqxx$Kg9!gZhHev!m=tg-3Vyenqb)X0_g12F)KTcf74?AMM)~GtJ2X%nkyn>c0s}%A^sC;;OY|iYuo1{T&+M^r8^pX'
    'E+Kj{TU4`444_;5ikFp}N<VCm<9Bja(x3Wn6lYSgf0cSmwA>@b0djhPL_Wu_o0Dk+|ui`hNCMXqEX#=xhgu=_sjr}Lt=-`rslhckne}JPXq*@eqg'
    'U%FteZeMfW(0~tApyrh%VfUEN3Q%L@CbO{o>BSlkWpO7fVZ4t>sW7{uKWvF4zQgGB2r7B-?(RI0ni|~^*=me()HkImo|yKeR8dELH(WFil^{RR72'
    'E`^qME8*0yA@1h*7F*v)*t1QdDB!ff+KUuc&4OGS2L=DxM?XCei>iE)<b=0d6h=P<LybvzGWnMvp)r=u|&F^XBcRhEn$^;~k72xXURNsg{DN(I%l'
    'ClcAJ$A}U0+T9<DQ(zuoMl?ep(zuD%(6e==5tt}0w&=|w3gwIsn8L$M`~io<{tTHlFr*{d{#LR5sd6rxYw|(d{}t_K3bd1hh`AJUPg|=bOx<Eqn+'
    'B(d&B~8zTl`jK$7LK<+SEA`m2zxXRHgfa$Kt3kD|Ao+`=V(`>zwEkyeM?4`KGHrPsG`9y{6w<(UCy#nre6-eX={OnFbi8Pv7%7ozRp5G`^q=iX!j'
    '=h+WvnJ_BPG#Ztxs$<@^NmB+rtdF^{q#(e!K?Pnrg{@!MNjhD;6o@iEk+wNk7V_t*We=xV^U`oOrss+>NmE<Z{UNyrzcqx>?#UE*c7WCDo^7g%=c'
    'a#j=rw@U_Bu6{U_<huSgt%Gc^0~B^ZxMh6e~hlZh5@`V1~67FmyzH6+m6m5WOYh-vZy6&SdDYBW)!+SC$i);3)ECHS~KRp+S;}xf^gxCi3h6Y8Jk'
    '1P?TO8Y3|a8P2nT{bE@Bcb=RV2u2Am7wQG^wb3rZMuzm`YXdb&Om<C-;~)4jLghF}kyD|`Rz^$;XW=zjmuHr6RKaJ@dIN4l1b;ZfF4kKS^Yz0Srp'
    'pBCAt{Ek5--D>-3S$t1o%3jLHdKS$V+AzvkiWT1r{>gbj=v?hPTp2Y<4-8ID+YW0HTG9O)HfU<#IdRbtp+r#H3<SeBiGx$3oo@Z(a@djVtO6Q9j5'
    'e@YIcDPa!dT1k*37~w0jmNFZNezaYjx_F(+1WrT1-0$9#5n}?SK3c;vex>l%8R!nkiz+#6cG~g9Xm+`1A>YdmA%3ik|-V_TU_G+MV2)vR)$?P`?+'
    '`O11k)A1UWB;p$uFZ97ld4Rs`rF94(ACXwqt{#=52)=E6Hm<@7;rmyU7lGAwjXQab94WL$yq&&<_4CL=N1H$5~b)ES*o9$$)8i_|K$(OmqZz+gbg'
    'TC#0HZaZuhC}ibbr$UBVzM4&2=pv3Vq6WYrUN8oCqH?X&8QphYRyD_;^fz_!LyzQ;Z>b+;}tdmyeHk5s$=rlOp66hhA$Ci*4%$QeX$(+WY(@kpw$'
    'CMRvC_+U3kqyz9VyqsvIgvPvn$2W~u;%6j+^Ry-Z;_+YYxd>NmT7n`~_+M)9Q1HxvPnf{?H~6rk$_T$RM9rjd6y&LfD!CJP_O7AuRV(*T8_UHg+='
    'B1MENt6J(;=U(@@(?IV9&cQ><C%O;@wngiU`mz>tTb!fOme^8POq*{~QV9_&wv8H-$jHJ}d<)fMbH%C_WNKsywxWBz*Vo*+m0STsdJbWnJDSqg1)'
    'A!<tj<(;YL*UlGO{bmu!*)<q2H#o$6Y@FCWrOBwGB(GIgFc-h9aeu7@Vpi_dc$V0fN0MP<ABN$jr66=X3q%?())uvywK2Mxw~)eP7nY*Taza!=`-'
    'ho?PIWJMwcCZik@q@wMo{hQ!}uirZHfn1drHbO5z^zB}w1AG&k#qjz$rcg#18g@T*IboM^vU80E>iBA~jS!aWTC$o<`;F1;JK!DhHXvE<--CpbYM'
    'ac*>iPS?ab7EFN_|2nu6u&;h#aFl(d~2Io({4@kZ19e84hG*O1;9?N74g~Jj}{`nZ({qP;>kL!F=e;=cxl-Zcvk=Y|7P_1;54L!9v{;(susd^gGm'
    'ZQC7ZBH=2L_&3gSLyf8AE=&MMtQyCuCgrdvddr2^d~Qv@}7HL+C=V=PMq7MWof*5cBgQBc{%{;*6kKMl7sWIt8=5YO$X8ULNbT5aKj+|e1Xq(mZ`'
    'MnPQveejSD34FzQ+j(CW!9%0l{ahRK+cLa+3i6*PqBOZbwG&2Xn;khsu>^vD!shA$M*dB5Li7udYBfxN%w1>d?h8Rm0;lt9xrH)yb)R4*vaEU`%W'
    'pQ?vJeyVL1+ki=2Xgib+*XJ!W2hL4?x2f`nJ@COQy>`AW$=ki7I87G;(tp&5DxqDM_Eg!%;}4$O41^U>D}7?-}`ardhNmH%^5o`~a3#`~YP55EK+'
    'hh<wwuoPg}3+r0F^V0JOd6*tX?T@I8$FDxbFOcCK$sO2HGVzRUy(b9ZpB%`*Q=oUz8N2dNoE|}yy0U`{=F=mkBH)%<#;Xn}c>coNp28QWC24PD4Z'
    '^AAfuHj%u!Ak?w%K=Z`(fFXt*(?mX9b2pSN7R-_znY6_sNx6qBdb{k=oB^m%W&Azoz^%9CIg3y%mku?guzBeXVFT@70h(#;w}iIZ*AhVAh9uD46G'
    '4idDY*{RW1zPQMu=V@cTrK9c6joAp|m=)BdSV!E~|Ab;ocu{@{T7T`=PhVE!fASTo74slaei)|GZr%zhc!%ZZOu*UIYE$<sy7pmn#nN1@k9<t$Wj'
    'r0cu7tV#@=T({f%c1!?A?ZiU+Mk6XhCQcbB5%$t+X2lv9De~U4HnySR`?Dy$OPFUrWkGrl--Ie;>QyJ|*hZP6dYW^fPc37!pIsifd@$mcZ~U|32s'
    'w{xoDDi(o}7{iyDK&Oq(nH04U*=1w)r`1W=3&W&FPxL5sSJ9S%{7O`9`A6mOackUt7N*pE0_KKa4j{`uUSwr};Ih$tpX>jn&HYZrrKnv|+<kej+~'
    '6ooSoOIFKb>Yzyuusr{Uh>Z6&`uItR(#;1jOPEG+H14vv=FJykoEF<+BkFBBqnFh71$`QelHTV(Dy$&JPw?5$aDf?gO)j)lc3Xc6z)Y@V6I+{z~u'
    'K^w(e`>!3vmTGo_>7&Ti7iXvU*V<K#o<9z%}bDkl8a1E^@!^GoZ(2Yz4liATMNmeVbUrb?;GfT2>f;#x+Epeex3nAC${Q{_23k)EMJ!|SD;1egJB'
    'd(8oWzO08i9vMIoi-A3t)v|B75K7C@RXEkM)>yZzWKwd40Dn$CgVr(kbVOk6QQ9#J5n*DCr{Qa~**G#7a9ACkrg=XzWT!#m##$3|V7L-$?+jH`k|'
    'd$<wM@!`VjN>c*0pjk%a5IWv744pq>I7cs!x5#NuGkPjY<rh39EhBQ@ygd68bTp_am=SMi8;ys@i*&x;s2={IpP?-#{P5ALmrjm(jH!l<^A8?LSn'
    '7c@P(GTZFvwI4>=HHL=^;*V$6y{;$7EmAwVo;K)HLvJ+`1Wtu-ruRm(dlJvZV_(-JvzU-c+bK;vNzYkoE}E$pM9favdrt9h4>8hCwDLbqE0jwcQg'
    'FHUI}M;<0caezAH`@|8KV!N%;i>2pB(y5<O=U`*8+$LS-BNPhwoO-XMCw~zuO^MszCOIlP;_kO=Z#}{nyxf?PPe3k#&8tQhzz0nRD_dLpCX>0#Xi'
    '9cG(7pff7|C9EVvzcz)`FeKchJO?{L}r%*Jj(SdF_r2&MXP8EF}FjyOLlNqyHi*pAYzrf`b<?r78o7l4pJ%kNOKetk-c#KiL-RhEiqxt)X3I7%LI'
    'p>%A}Sm@}&x&+etoJ?~I6v8%tH~r`8_(>F2hqmQ6v*Tbw&f`YD{&(Kpiwd>NO;#8MVc_ZAaa#eaI~C&}2_i7u#I)S#gWV{^bN;69y$kCr}RJUU*{'
    'GObgA076Jj-^9c!SDMaYVd^o#ZY5bYLN!T!Nydh2%rbO<sI|qYZM+O?isGtF35ttSon5Wb{lurE^Rw7VSEmQK8&mR1^XZ!r65swns8O-RFQZ6LBR'
    'uuinSkDMXhZsJH1pt8u(_nOu@=Cusz08W$Q4sJOV|M_ySUTk17417>*{DVvR^mb1kc(EgPY+kLJ5FyAJx0GQaJphmut`cY=38g*@lGLPE5^Z`!=|'
    'PHbcM8=|9G&i+h<t^?dfP5#j4kc}zo59zL3}49W`GY@8NsN<d__&%ms#mXxW%yDpUo;1aB?50^K4u=yh(bFzkGm+A6A1DoV5Wb`+f$)v>qpKnCzY'
    'rAnv5I{GA_?#cuo18u>bZ4E#>LH|B&Oto9yk!IPEiG>|`V*o0NTI`BigwDYLPm`KXNziz%?Oj@N>vw}Y>n1rj;l3P21oi1d*g0c6w33!(@>8une}'
    ';x2f=4wmuW&-y)#)vjPq%3k}Y|1z;k6e2c0skp$vPOm9T{G0`x=udfy%w)#ylFHv1*^1>|LJJY>8Y^kQ(Upm8ywlqAqsf3JzD5}3VYCt{j@+Bev`'
    'U>1ZH@>`D5(Zz2-I!;gO!c<dH*+(M5K3msJ`}sozLCcOVvoBU=FWdpOz-F=>G<x?BHKutQV4AY4+VN=3Mu5>Il9EmLiF%DLWJf-B0J+uTP)o)eSL'
    'p$~6Ma84C5ea1EsFsTP8@Y_(u0aqcXQTzu}Bhnqf8;1MS)6+alo{u&_Zjq)_|GXvUad@D49%~5@KP#W6Yc3-*67OxR_{^KH8u?lhjB0x0?2T3pW`'
    'we`<8>++ujxcf#O}hI}6;Fr)#N>C=k|TUzv70lpjdx8&mVQ$3Tjevse{-Q)j@k}ej^v`NwgY-A)s%o`{LR?m{(v!Hk<azejw=QToL8ia=;scJhvo'
    'yO-~JTE<qlgt8=YkTK)7)~QNzy4{49c&6#9|5(e?|pA5DSxpy^l6{Ny^{TY?%E`a5x^>Z+LAO8!$9r<qi>M_gE>*JRP9t2dbRA@^NTn{5^xe1MYt'
    'TLSdXi!uU&ITzhub(;Pi~L_z*bff0{}x$GTQ<;|sSAj~C`a8*QQ5*Yt?Mr;u*|)@wjLI~ZirB`m@ZQJD3xY^lueZJ6b`jT5Lu1!pvHnB98c(|R+6'
    '@22pP1K_|IpBb?0L;HVmJ&z4~5IUyZUz&2^(GiNCcRr#7t@qc6)85x7$1o9Sxr&jOq<{<<-}3;E@H+EsFPrDutMXg(hs$!GmIugdl|#hU8<r;Z6_'
    'Y&PLV9f?8<W}xsz6_TkBJ%r^Sdek0Y@zkiQ{(MenG&YbziS$IZ~|wz(<nBUN1g;e%N7|G=VQKTw#80$l@et65uotip9}0G{<l=gBZK((A#RL0@E<'
    'ub2v-Kqd;7}=x6$9u5RDwCOd&i=f;S{!VssW0seDG8=A30dz}L(uS~c{`(;I-qgH0!ab7q0XNoGWxT$Z&WQS@ZhV~~9aBor0Q0$~_&vNesLY}%DW'
    'IlTdiDF>B)&gQt6vg1YeyFqigRf10Iwk$#e08mSpSj=IPqsvXR0GACuLKC$@n|2jSt=GGHv&GwA}p4Y?egZs+96|rJAT5r9>00H2W`9!%E@v&{V>'
    'R8gN<RiA?@AORZo{0q*wH7hH_qAI=e1##EP9?3&1XXbQw^F&xk1s4W#|(Plle2=^>ZB%76r4PgRwqq(K9XJ@YekIqeYRfV7;W|INt`_2JQcT_db8'
    'fbFXdbnW=<3p?8PV=pF8NgXmYn{N{hL-+EXjm<N*ggWaZ6}y9i^zhPFdZIAU->wd18HA2)IM%zBVtd80sXqQ0L-SC*f^{@CrX4h{t%pZ0Ah$nq)^'
    'E$3UtOq_Y0bnOFL-}+#2rzjZ$D7KH?9YEA@2PDedZ);67W*#GCNctBSu^3nX$<_1)dyX$1{OREq*8fa3<y-DB8k>#v!;6{plrnQ|Pl#g%AMC#haY'
    '@Rp|tpL6z-z9ep?>rkgc`nFuOylAD8Mf@mhFIT@kSMaWt)d(m-LBNk`RwC-xRL#Z7K4m;Nte*}TRqi_YJRvJhIQY=>7Q(b^@8?)e8)j3p!7pg!N0'
    'o*4cL=`3;(<AgXc!SG6xk~+EKS~HLN%#sMi-{7#0bG7=-uZ4)Vh#|5HIswgCdEe5P?F4DcJW!9IQ;@bT??ey+w@8X?z}ufI@?YS)_aGD08){oZb4'
    'FJSr5tc)!5D9!vDqkMfK+DCFU2avO-eouo)pfVk$*xQQFI-u2v^$9Fb=scv0)4I%dFFB0S5j*yEYvNrC5hVY<wzrL819ce!9Zv${Q&1_B`&QXkAB'
    'S0G8y(Q8I<i)AXH4iI20{4(~9QrbE!0V~z33j`3+Hd$LzWXUm`dh(2N!3i1C8Dn^h0q+b-357?^?*AFVOewYGgE|~WV}Z)}khzMkMe;gT@<VA@<a'
    '(|dxVq90!Dk8KJ=Rl#u!;~My;UZ>nMotU8O<0AqcW-V4r>NInCcM^aqMLaK6kls;A~9$g|*G{X>;>F+ue0Xv@GRQG+mF(Cc_zZ^Xm+tZC59SK_ZV'
    ')VGB-XY>M3j8vplSZD7IPl`s7aXKvrWzi<pA_JdN5%$a!Mev{&uu&$I0{oXc%Wm!W(cK;X!5rJ+{H$(EBf`Ji8qh)Y?BkHn2`5&{gVulG49pgbs&'
    'PniGH3GfigclNYW+J~Xz-Wq4&NzhNf#CQfS%Qn-$(XjXk1CXN#oY@-!^?iMIYT{yMAZ90;o||`;m)itxO5ynXR7nHT6{!5wk6|VmTa73yL@$*R<k'
    'LPfyPA1s}{C;5G1F&#hm;$U@k#a4R~^A81s9rVfDI<BWZHoltC>m!5~^7#{YZ&QNbqcuvQNa&Zm_W69U#=&ty)ou&(MJ@URjdK>{JBF+9lEYRsmV'
    'lTZla4gY5{f8L8uX#~gWE_fuOk&c;;zPS02XluSyuf)(Glon@q7J2bUt+;O<-$9*8lr$Pxrf#W;ge=TKLv(7n4;}uH$bn;o6YoG4slGxO#{!UvBx'
    'a?Q%C*_{6iL*jxL8SrbgMOBOO5Xic-R9I<<Rq``=RA+gr`9!>Q^cwVbwlV>+BR+N?iL$iY4{t0$LCU2(K6r-#_s5VlNq9S|ew;bqeqoE(Kl~$!6e'
    '~@WH#!FPnNfxZDXN-+~Vf$!)XcPMkqQ1S&QBm{o|rLPur#^KNzkDRb3(j{xSZI^>G**vZBiqb}y>TKfNBv35zgpgtd}uau5wQ)Fvf5Lb|<NZ)c0M'
    'w}^W?)BGav=Tz4NBB^D!mXMpZF!!Hq^lMw-eztL@j-hLhqYL(lLHmo)W&le$8SLO9@{OrDC>7d$fqxvV)mjF%EFG)?5OEEL!3kCgB)CgIok~6Dhq'
    'V;iLQ)1JL^kZ_k5u=hv`kKHD%om*3^IZ`OjQ}4Bz%6k8yp<TZ%HnII?_eaGtr~+4Dh05HO}I#Q?AGcSQ1kb=OSXI_eVOxW{)|&va0?u7z{dcbU+d'
    'Su+@Og&T44c1-I3H?Heeqf^o!XL;d48Jn*Qa3BZgt`Gh>INI!5G~=;PhQ!r!5{u-8{(+l{s0hf|TV!oR<Z-%xQVnq1>6mq?ykiIQX@zcc)_km$!}'
    'ndz|Li(2)=(#LyfWPG8*B_MnBydE1~1JHh-Fpv+n_WTF$1Pv-18sf?_?um7nE!vSOpzjZ;bFCkNJ)G7}zcH4z6^p)UV0y2jUkuFcqRYV;Y!C^;@V'
    '))}#)qV0Wa9G*q?tLWnzk+hx|rW*EgwrBbhujZ`!gPJz)a@5p-ttLm(<THY}rImw6}y9C$G7z(dq-ZisqsM}{pa&U5~bKPi2UPhf>bZPXuRIjuif'
    '50YLK*5&w>>6b|ky)N-^P)m5vg?ZJo=qX|Otem1$gHyNN@k!uNV>8O7XrkD`mQ87JYJQ;kx?l-d3h<@8w!`Js<y}<x)+kePP1Dh4avh<`|ed9@F~'
    'vYmP)dNt>*(0%=@TS4*|gtn<E}l5L)#A0eTwk@Lue>rmpGJ9}zUao<fjomlcT5NIhDxZl54Oy|-7kNuT+uzw{##tdmrurrceNt&bkr;`dx(sn9X9'
    'mIlP`r5|Q{*J9?|nb>XFH=gk?Az=a643UC9Q@G*Ok4zCQ!M^I%Fu{3l2u0^C9}K+>qQL{>dK`&omf$bBU1RuAz@=wOV3BbOFXuyP@n*73db=pOhh'
    'PkCl}{;+Hr)29?21-j>||Eojw2F_=V2!lG1O7;pVX3oRS)p}AsN)pp<4{XGiwco;)Z#L>U%3kk*TxS@5tzi(vU0m@Z-Hl65JH_i1!qU>&n3;ucRV'
    'yYRp=9%NqFJ_Ne&!bK71QC5VdvtHM}b0VAb^!6sDKNoM=&Wy~=uxdW@52RmUpho}W47MzwdwV;okOz2uM;f&>K0-Z}?YCRK=DqH0$z1=#ma5O-F0'
    'QX&3;&vfGm%NcR($4!Y6N*CtfgHdTKZy}wJ<5j?A_Ox#vkXw?{7~fCe1>o`S^MPXi#-{49><8y+Z}lBRCE4MNWd2o=Z;`F(HPe*Uc2MR`iizSvfK'
    '*<M@x0HW-wHn9m30ky>vT@sIl^AU7QQ6Mjqy2;NM$9cw!$I8JuL7x2Xvh^HneD!zl06Otj5j(L%l0b0@A|)JC!!`=RW^w?I3*b-n^yFIIESz~)_c'
    'j`MO`2^t~oh5k@CGHM3PGH(Y2sdqsc5vs_B{H24I&E||o<=2*$F--BG^yC#aboO@=pYC$>`I%OvcOlX%GFb0GPq}veTKDKwFJtm^JzAUDy`DQp?B'
    'F~xy24W$dKP)2F`LZYZL^r6bpq~#%}i0S8ieAf4-4b@eSO<Fn6ODj1r{eWj%V_Bd^o3B+6udKY$RCV74zOoj&{jT|M#pRIovYodb%(?S+S1*Qn}E'
    'nV!;P*(;Ky1BG6F5+tm#jFmbd95X=i?SN-LLpv{Dod%4aiPM@Ly6{+J=S!-bXqqRk4%8u>^*zA*Y!^d4OryRo}W>YO4hMH!S*_nxV{_<sEKJre#?'
    '{k-Zh%z+)dhW;kl4tUJRPg$`J70s%(;cR>Sckr|rl^~zaQ7NKDSx%Xn%*UQM*`K(NpW2*<9(q)l8bxn-z$BnniNU+%<j!Im0tM;YeS3D7|ZwO@v1'
    '(nqRE7Wzx)yZFB%|hG5Z!u>GA{3J4bB~-}g5Qk1UcgCAd778ek>;7a2fdY#w?#mF8~@T1_iO9z4k^Yk4tx2Su7$Ry?lQ!5$Rsfk2t^q58w{r@DFW'
    '7+*KI8wg%cTkZ{avxG6v5GV)=kq}@shlYBW5C~qeaPz@jGK<z9+SiRm5a<twqYCKeMzc#4s{@v8k#a7t-0;$MP&2=PCyrxA&r@ViMQ=RR&I=Fb63'
    'Mx^pVj_6sjN76kt#flw20-+zdd=L&{oe*{lC<X9YQ-B6ZH<^kF@x=@poC~eACaP!&>U{w<|X3$R~V*a8AMm_OW|$*<jAj`6;?;HeLeXHlpP7?YE#'
    '!5&jXv9RMgU4a;ArIJ`iZM-bURm~%2z8rA8WvsknG2f&1_UPcTu{gvkb){gc%pP+YF*tpN5(u7{qW<eoMn{Lc!x=%>dXt^v^7O2`A;dmvN#ku5wl'
    'nd55Bebh*E%#=tN&#s#I@J=w!NRI@%2V+@ukm1zt=jKUr5q31bKI>_!l(C~eiQr;IBCK2y1cr^pHS1@fASWUQ$wAwsp~NC=x`OQeLIO!G$1=#m(v'
    '0R0d;d{w@nU5=TDQJfJ)niz)co`3NZ1p>IGT;wJv(&h&-jPaK;S@d@|fa&{9HMHEs})OM!rr|Iku`T7@<eFI`<~)?2t`_trSr_fh(a-O1fk+Mb#N'
    'piw7Sf;!D90J1u|X=OEc=lae%DfPa$Qg@H)ZP;8dn(g=YFl)G+lmyZu;{G17#;|ctEUkW2zXFZc{tm5TtDsJjLIi!)Q>&5U!Su&*EQe+S0XLt}!F'
    'la^y-ZKZgRe`wdq_uM7YUt={{TZEd#j%M&uU}x^yz5tyJN|-TA;0UGC1p>?HodxJ><u@Ly$EgST;8shkukiE2wa!r$T`;xy~ZahO5(>4PvPp-MHC'
    '(@j7`BP(8QdF7`2ilMvStrnoggL;VAYyD5W*b^XhB@3H`RRyt}`tXsh+fpF35vWEf?BA3zKw;)S}D_HG<n7?d+dKq5puzG&0_m((i#mF=VsRZ5$4'
    'Ha+g`!>r=orMdLWkXA|HH|2Hqwf;h+Og-Z2ubXcqybGxAXr9Nu2{6yzWNK<t+!fpcgjnT5RC5sAY&6s@^wU#8-kP>Qk1%KcTa;sqg~nufkY^yU?j'
    'Yc?HjVRGKrUA60sKySn^Lp(tu<8H!ifo$gDcuqSd}Fzc9G4H5boWUj_|4$PsWIQn6n$YVA*PlH@1=d4r1m;XXs7V+OIbC|UPJfwDFhDH~kS71_OL'
    ')Qk&j3RdUDi%T*JSd_%RGue?D%R4~_OU^ratIO0x9)cI;t1MZ*^6wW#wcyq|9Y2cxdNxD4n7Xz}W84I@l3B`(NXUcX7968$^ZRA-902Aw1$^p(up'
    '|$x<CSr%^$9l&Tn^p_&!7&Eln#swpQJ>UZTI52eDj#^08-8j+fW2n%U7mDkYR$Rt$s970fzfxEs}XF*4chvb4mqx)ZDASTLlrm*c`vg12Oz@er({'
    'XnJFUylDEoPJQ9ZneY;&n;J}x-#(TVsG`zUoH-lM5k5rFl{S5@U`d;*H#|^bFg)xhB3v?uL9J2Pp>EHU4PJW_guM;~caoevhnKDR{^t5+s`iD`ri'
    '4%FP*nO9@j66Yy=6gZATF>QA7iwsqGLI9l0WNP$g85wMhS78(R!?GBCP~158JCVE<7!hWu8>fdWEtZ75Tt+$gsks7wFYl2tK!<^BxBA{t`y!@i7*'
    '3(%)m@&09$go6Yan(P7E;$3yU#XsC6TiS=<M}%fo!uJ+?ATK8!&Yw@v>ZV(u-qg5f>f;GnUy3hJo79+8!%yf7;LS=r_Q9Vx30ZoPa*$mw4t7iR;v'
    'G37_yri?;uMmqaggopUOPgDi2(N>qm;C-FzZ!M$Aqtg(X^>`bZf>9cd5Wo?~JV0tT2aZKF)c=QhB4YA3cV#uxem!xLn-_ylC~4?#=jp5%u32avDX'
    '{8Tq_GJFq9A-GaD}%xJmnu|U?vD=Rp*v5K}Uhfu3hVwa2g3wuLwx?F#<=;ddxJ?;Je6r-f808zPd!;Qc#iERT&J3|APv()gmN#E?QKdo{9Wax4p$'
    '1lg=jkDyb*^LkSFUniQP=Bl&v-W{U4ANUS7JzWBE{VNo$DZD)fFFs`K^H}m9HEnyQ~$j}*3TRTF2$+#0nN42$edJ4X&p?j;NG~htZiMiTEft#UOC'
    'unm1E=C;jt8}njy2A_P(L95<lAO3yqF0R9T$<CKvrTNWfbuVf`e=sbD#ff0uexVBZijv)D9?^|4;2|Rf=FmcB{#<teT}#hTm3cgG>Jhzm52qx;fM'
    'h&Fgc~bt1;j)aM6NWlQm-cS5f(ANi}1Qn7$l_HrGQ&qt`N6KZTIV>%rtAMEVvp6r2y<Yd#(sP_>BxQp}%){6f$Ejgbi%k*c?h#1S{afe2Xn_W0>7'
    '@<DGou{MP$Afl|$_zY7w0na}lMlt^J0!YwIOz<rZ<c-snrvYjzDsnB<pXe@1%NOl8;O)CqDEH8!C7BQsWB&k;ZMv&(6WAwR3kB`@dSR$rEB%i#{Z'
    'ES?4eyDp%kd%kT6FUq6P3ky-!G}Aif^2H`q486jR$I*S#@CGg1Ji=o3;l7q`+P<`OcIiiZ9-T^#8Uujl8{csYX!8VWnOMM#oftyZk+Pc9(Z&@T(V'
    '2?$7T$Ux;!croRNevZ>v%N8Yaw_fWqMXuwo%<MwNk?UG*5pL|>Jw&xaxPIL1ZaORYjVel-Rkb=6rDBuY11okRb$TJ4HSR}R*Cn)YnHWB8h7lKx41'
    'X(N9wE4vK?S`P2Em~?Yiik&%B;Nf3WO^>F_UHMGflYC24824J7;Z~{{HXaF@4th10Fpq?hS~?MFa$sV+(Sc2FER|jf#;#+?yeM*xxY}K*khDy+O*'
    '+6_7K!8d|4gs93Vt}0)*W_klRUfQp>0TJrlJBgodm|8iS)Dt4h-__8Ft@q1UiqAG?GoZP>A0uN7IqRxuI6ToK|bmr{CGnw7U%Mcu#FwiOocd@hW%'
    'grJ$>F}7XYF5>Erin9lH&Bj+~j4KYNVr!QVGX){FV@`7OSN61XkMHVA;ZYyOHXWBW3mL)W6vYsF6oiU%SMgJ^5H443KD^_RiKcu9ObHd(z|1+b^S'
    '^!%!}}FA-?k>eIr1CKE@PMLWKzuz<_YLuD{@kR0j4{KtJa?!1dzS5P(_Bl&yugPyR0;3P{%dQROCMbQ4p%0nGVya=KN0qrh_7O-{T{EoSznwlpK<'
    'P4FhEq<U0h-Op`LXxP>&*f&`+YOUpnntWhF}8~C)qu{HXv&)w)4{g^)9KC*3(5$i=m)K(Exs`(d4m}#L~)1f-O#|t5SFyH}fO`$eLSqbb%1T9?ZT'
    'B}&k53j(r<5Fjz55OO3mnrG5U??M|6j)(Rdz4kR87rpmS2NK%y~<&1Lw~(@?x;lGB~v!#N0MoXt<j8Dde<rjG2g{_tWgjF;aAGm(6bIkM@j2+Rgq'
    '}Cdk#Tt-Ln!>MLhwx9lPjkKc+{9jRCZd`E>Cmd7+BFs0fm{zJqL;G7*^{@J;P5E{uDSq~BQu)?riTfR&_Ma4M!?Lb%WkyNz>j<BPNW{nIq-7Q^rK'
    ';vpMzD#pOS+shkn;=t{#{N{X(qhi1Zr}}^GdG6TpGQABucI@vsgK`cv29I!dyRGxE=Ur_NO1*t<#BD*=j!+nry6AX;L~{UzEei0S_mm+akPr1^;%'
    'IGvDM@{03bf#++Dh&dj2h^63_Zxn91&hWRdR~yX->_GikoYsoPb&6%jEUl8_}!R|E)_^D(scr^Ah8Cz}yF^eT%U4(1d}Wpg{=i90m%;Mdl;z_30X'
    'M=JZg*I`K%uTa*WG=@AQ9<YB`oy`~O<`GaMO*RX)$cpDHTxp>6naXz(_hWG<<sfr-<e9M3+b@x3@69^oZfAW<r`Gzx3lSN)k;5~BP?!Y#L4pj`)3'
    't-sQ;o^_@7O&<OErQcoMb(tZ(ya})fi>-Z8=w#Hxm{Brxi@?%-hxOV8<9^xC=ir4bUC#7yKEu`aP`=G^lw!5^NlZk6beH=67alNK1rXqQQ(9`#DG'
    'cf_uFxpIJ)y;3JU|#gvl<P!|lb)0FV;AplS1^hk9}e{|`xAe^`5f2tgOHgYW;V>{2*XhE<LS`)4{-!?UAbK@gAkC@4{ZYg`ezxiAheo?(94#8P^H'
    'B}V<WLQeh`ApP}J3=2Qo#a!NS586byAh>umTrmW%rVFoWGyUcnp%JY=RBo_2p~<OtFMu`#cb(M??&c3Yd}7Bb4|~x>VIBqqR@MwNV4Hri?<>KQ7$'
    'huFW=K!j_%j(0V)>D&uvC`2cAw>_v`CP=ePFpLkRttf#yny3yUT`^*jZsHw<t#=wsx&kU5d0%;6us<{$sJRllV2XnOO!%+smS_OuLVZby6^+4?sf'
    'o5fv|0jH3uKpgiv58Xy1X5rNT4WDG!ej>g(Ed~cc{yxREYiKh#V@|!%@Vo=`7G%{ff+!)Iv7$)pU#QrV{)lB$tqfeegZ;V71Ep)5Rz*5+G<j?VwO'
    '4t|z2v@v`8U+z`;?FF21-{whXc$$ovr}X8+y$EDwvOxevbZqIi_<K=Q2P8d3$h)MDo8F3&a@SENGf+^lY^*BQr)+dig${>^W3m%x8H>^(hfXjaTn'
    '{k)j&0scS&ueP0+iZ%*ScCD6({o&_q?m`ftuE@1_Xh;bJJ##D6%uS*H6E3T40PvKhQCXw{AY8{t+!$f>JZrI2)x#dXVU|KW*%|Fb>Eo1^6#JwdEC'
    'o~W_orVRcDiT&c->d)E~?rHENkr2zW^tR8oUCqBulwCR!RAMp5+lXeXc)Q1C4I&<V-KqZ6pOxy#8^XnC>OchmUR}Kh_MJV?L%%|x5<1NMZcOFwU;'
    '^jl*e@9|Dv{AI9iuFqbaw|u@ilxoX)$X8qXKUO3Y=6rG>I2JFoNbnx(r)fGe~vpTJwIuMFFSJTtQ<wi?L+wKIVV6$5+A0v7UT#6c(G;)~_LK*<YJ'
    't+hQ3<l0PA$4<VGlf!eX21q`$+`}h*OuNBt#s`9^;zOp|d$S3)*p+yp4?+~HewV_DSHjfTgqMbRlQcxJ54UFU++v$?dX2bgiMn1|}%RuoW*h05y8'
    'Tgi<Y>=SgG_zZ-wy!3CL`k#xKKt8qj|2Z)w>OS{wt$FinByhx)$&)6w<<yjL;mNGmC{zY9DHF|l<YyBe?>xoIw&T<B1V}X17*0mSmH?&1WhgQ)rU'
    '@4jDH}2HBTsvxSL#8%8DO&KtXw=^1kj2OkIuV13Mpa0I7oBHW%<Fif?N4_bRN1n9EjvM%myIA5Vjh00000+@PQru?ymN00Exu0?d~lktWW<vBYQl'
    '0ssI200dcD'
)
