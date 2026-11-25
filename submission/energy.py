"""Active -1417/4000 first-pass, gain join and starting-majority proof."""
from fractions import Fraction as F
from itertools import combinations
from functools import cache
from math import comb, factorial, isqrt, prod
THRESHOLD = F(-1417,4000)
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
        # Pinned covariance and the signed factor are symmetric in i and j.
        for i,j in combinations(range(n),2):
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
            deriv=add(deriv,scale(mul(div(sqrt((1-c)/(1+c)),pi),prob),2*s[i]*s[j]));faces+=2
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
    motif=G['edge_motif_gain'](); bonus=F(motif['gain_lower'],motif['scale']);pair=G['adjacent_floor']();cost=c1-gain-bonus-pair
    ck(U5<=F(73670242301,2**39) and gain>0 and bonus>0 and pair>=F(11985,10**9) and bonus+pair<F(1,1000) and 0<cost<F(12307053318118561274,10**20),'active energy bound')
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
PACKED_PREFIX_SECOND = 'QB;#no1<AJTXz6%^G(}+9*Zj6uG?_EAwr3{<dczV*}k(`RRwm43t7;W{?Zk))Pm!>XnPVcXlW4i9a`K-MKuh%tdzeAK{+B#*yHv~tMH|>J56ZxCQ@jG4;_x$_g`cM3E%pOH>(EymL$-=>59$*lvTsv`W?~{$B`SLUh|@oh=E-P`?RN;SS-~i;?eg17G>n~i%|Y0!B?1RH2WhK(irFGzI8KyLcM9XtzLqD4u)g$Bo}e4oXFzMzBd2|6vJpTL{#ssn&hNoC`M55kV1qKB-2z&)y~q3|1v&!xZ1L?i+hk3+&Z$jAoF=>EzE{BJVel__lMH`+*JW)<;-I-^fht;v=>qaTS+9tV5juqLJ(8%sDh!FiQ%~CxTC`39v~Pwb3YKvcP;``?+62Otlc&1|Hg&bQ%UyZc{l7XBPBbVjJ(Lk4c2H*GdpUAGSP3<7DZ119xxvsY9wXxC!k{CdjKioewHsW5+sd9y~6o9oKb}NzbhIb6YcO-{35^a{qPTH><K)c;NGXNe8bRIY+mkk7?~bbbWoT^;gPq^PGc1(Z`T6V03k7pJ)$S?r+!W?{y`>|#^Sbc>zfYpy~z9>#mbD!@IO5+Aw^{@DDkhnD=K?Ih)WfX)`cq~>LNTxH{Bf*StD~GLYitJD%>m5Ly+%^uWQ*=@tZ~Dli&z39;Ute(<dC2^rSpeCCF498NhUht7YGm5tmr~7oLd36Bmx9w4iD+m%vZsMj5R8vA;a^?6J!LdY9Ydtt-B7ui;gda)~o_(7E)JGLG91v(IoK+V)`w;vsDTW{rmgdDq)ey<kg@owAKVI0pt45&s@hmHg<o8Onqe*nO1~veqEBmD(8FX4a-B7u00#<1PgF7mJ2lD?UtT>-mEp(F|b3rTCYN<B9zk^(ura6aJs0!a^txX6665$F7b=CFqUtg<kefZSzo6`m13|;!gNLOpkxicgQWDi{2N3?-xzW;#UAjx7hq+aR%usJG*z?1%%1C&%z4{4Mb1IMhvZ6X2`52<1sK+gI|<+cxEy=0k2rTas^13yfLWyFhX2b!&I_PHa^5Vz(K+@(%h+u-Y~b5_>6g!mCea5?<?<A$?V^px^!l9sf0D5Ip9WRrhV_`Y*;Sm!H+@s?a^Z$G1;eP5Pfj6D?j}Pk0VW5do)}oUEksp)QLmb(=^**r>fQ1+xcfAt{1-AmFM#4)n}uMa;ep1)!Vt`DJYq@%&O`Cm3swFLh~RhQ0JV26>%(3UO)L?)A86>a469SLn~*q=Jn%S7?BEAGDV&3(B5%(M*WB?Kg*8o35gEU303mA=)0&Y@Fv9}p%q|K&}sdqY<5d0MM<Q`A75mUXMo4{atsm*-lDy*@;V)6l;0q&bseiZ%U<v)D%?Ut1b$5T4q(GyWv}=YuU8t}o_1gPYHk>Y#^u}Y&N^VRk4oy{b4^m`6Lte0_YXMxF~B%+p<a#^ML0-hQ6-}5vzlS-0)fE3A=$G*Lns)@VBI(xdSC>4EMTUv9lpiDH~1(@J(^)9H^2AV*~Fj;?|UG@<=|NZ`Qq-5YFUc0Gy=4W>Jb;3&B9pZo#X8B1%!zM%4K6~sXCYe#}IJ+9eaH^?Seq;{FY2{tZ!DO@FExQ%3(eJ#g`tGe(<Bek^7gy`@nlfU0Mv?!VKFw9vhDrQgRGyC>ZRfoiUB7N4d=dVHx7Fb^z7E_&KT{V&=f$X}5V7-8c!{-e&&|7zx@sxudJB^)<cj%_%;f;eLWSXj!;(e2y~PNRGzWadw=3(rcsy3m3XY!Jg@Cn@5oYrqJxmxWGN6{2#eG05pucF2)YS#pFKx@*>n^Gota})j-5T>FgwUnP1t7nx8Enpg+e%i9(8|74jiYZEAoK$y;n6!V<dUM0hi?Qqx&;Fpb+=B^R#p`t=}uw=@r;!z#crryAw-@MCTD6UI$AeestVkJpckMlvJ}^M5h$7{PwY`SD+!(;)DS$7p*rzRJ(v@|dU1X71$1FFiMp`9`nuqS@UuXsn)(R%*0jVa_knn?{AD;e5d|e%V~(&I_-V?E**hCKEi!hAw<>O*#2&JMgu1eUxlLKm6^O-GnQAi~^IQ!?=&e$A0(7K)2=o^S|K4*?sLVbP*(3@vGH{Ir!NHpD(c{GW;RgMo-hdqt^8RDzeEbT$3DF{5kZ3nEK(CFUGaP$r^_Pzo<ZP#G70GPztKTynB?O?tHRj6EqS-i^bV!jA4`_5pu`)A|=6Llo~g~)9ZyDWxf(vY$NV#O9kfWN6*zfW9msVM4ei8o4C|V;qlrS6sYkA5QWEM>C4Enq#sHjxD&bMaPvwD<Jm|<8TNBFJ!dLz`FYX0xEZ(H2%V7;Oe^55)*aN%de2q`jT0OgDc`%IVi0gSiK{xPp1?XeVPc3Z;FW-;?(Z{9X$3Nq?g~1Zqy~nokIGa~670t_szNViaCFE!yb{14?&kWrK2eT~OYHzpOY-*66Ja^ms-xPA8U2o9mO5@q1vt(q4%r)>7>2=?SC#tjvtkmCnTj%8F=oOrRhiZ?F<XQ4ssbKcn}H77#2O~7L5poST<i2{Ve=**!jx6w$?=-kWh9O+ufA*Kg-}q<af3soP?Uy+uu1G2Jx8ThZK6__glAE7mIaNU`bvcXYQX6?L!FV3L)ly^JWdQJ?`>Gd7%QgZ!%wP_p?v<;TXjomRQa6gHG=Lv3va3zDOs1056poKL-~J>R3?t_dh&)I=LOI9D$r~0P%*o9#FqU~74|Uk9Vv)=p+eShOLF=tAxpf6f`C634X(?ni==miEH>%wdxZ=3(6#x9MJRL-HWrb!;(N((>ENVHy^XfD>9@mkY=S-7t8;5NQf_VZvVX*ckRWj8R0^NxT=N8E=6+*W<eoeBgn#G%G$WidPcB-v&<Bf`&48!%I@2ZGz`jeyAB>T1(Hurnt?T~(e&8~%xTVa-`0eyZHpysa4l<%wr2E43J=U9nVaWu)^ka@o{vvQlP;Pz0U2{aIVl?*VX^Q;r=||SZ3HGpq?F;bn+D4@Q*ET$FiS!Ddx0`kYYrER`RJ1>Wy?52SB1$dCw*(SP7`Or)Gf#sQW`|bHr#-BP44n9jQA>|73rIV(_4#Ro&bnA_C56dQwczB_;<d+!ea?+D7Zf<xeh-$J7rq8p)d>d#jS2aX>O4o9t@41@X$IF{apYwE;^6gS!y>6?6556PlFe3s=@yNW`QgPz$|!6~c;zib@ty3#SK0&-*tNX=MF=lV!e@}md88gKqV4>H#Bpi>ABQs5W-48X)Ze?%3N>4*ZyHZpw%AP!ipV|gWRLP(kR~BEDN>z|xL=y|503~Ur!mS7dGp`fiYt3Z+ha^WdWELJAy|Y__dLLQ3b-u^Ed3w8hI%H;g2v0IYmz0u>M>vsjmAMwp-?6-Qgo`eTb{EDxONkm;bG)#cb;@b#Sj_@lB0zfd%i2RpFbJGbH4r@!EDpQZchu~vz4^}4p#r1t(|(lOS1jlc}A`L<ZwAb#b_Hv!tOg`lypC<xpmb_snG*BSU=d_b;)YI97B~hewyvr^Tgsy$rg~xigf7Qp(`?D69DDET8Vs*%Esz=md>cN=v~xnoO5+L!h#1LN2B-FD=OeD`Dt2*7qVBz7S|)FkYc0Kccad*>IABz*;U0i`nO88-53k3TGdm{5?9=;hn5Oq?+G@6<f&_SV^hb8<)uF3FU1r=S|@9}cTI&dbv2rsKMcMK4vu@dB`x9HH8oKsFRQ-1-~rkShVc%8NV9$kejQkeu&Ezbo5Q*WfcZLSW@R+c&$=d=7bK{Hu;J7PybphVpjrU^bGt~h-I9YcK#7xJm4?<$AfAFKoCb-W*<y;2_t218VZ|nr0kOpZ&&q^&UJA5&*B$)eyJILh6mko@CxB#^&jEDeVm9nCJz#BZ4!T8aEOd=a(YRXu^X>aI4S0eO^54qnzTlTQ=Un@F_zR3ZEaDnMpRkt_VS>>FdRs1$;iT^pgBIswMYqNK;ZXUciDX`lZE1m_bLrjjYkJV1zh=DoUqNpxGwnT)9o*JO0C;(rl3O=m(;N}BuU=<t4PKkwg)k<0+@9{$kmLXYe^Cs4MFl@I=Ug0Ss#zm}c&58Z{-u50;2uQlaY?_y`C!XZU6K>OHn@3nwl;O__^;a_hoUa-E;Hk@+O)`fhjVueM5~1fQ_1zM*GmdJ3LDP|qxOa|-964Jt6)0mpfu|2E@FxTI9Y^ka6{C{;l6Md1-D-f4lsWZzoXjTT3x>-i0uaEo7^t+y*{kv97}o_z?)6vQHVQ;_AII0m7%RgEYv$yxKDm0tUM7Q8n@7kR0<RW<2{?6?L=4--2hIHTqT1dY6BZRp(%kfINu^r1hIZgPwk|4vp{dNZ!psh^SrSPEMFLrQLq8^IF+wLR{VC==hPCq;YzOtv|j=h+nj)gql|8?2I3k<g){02Fr5$?pkV+JK((7bhYL>$`wS27FNwF&nRz5gKmS-@9joD#svmAbdHn?zUKG}oG+7`OI2#QhY83DBq;!LUyRT!jDZ`lepUlcFVa5Gj(@IG(B9amRMoa47Yb<O7m@+GyLzoyDcmKNpsKy2eLK?_u)*6LcK#2|T(ahg;O%jI0k=CN}`QqLv6XS{Xg!lVKyA#pNa93-39ao|airsu#8)$PunbmcX6LpFaZd2JWE@=5UOp!%jzA*5{=^yy6-P53V=i*p*CDzh`B)mICFWKZ5_YwvYLka+)f?XN^sY%L2zmb~hgE&v8)v}CR!nlWln0`*9f#k1XhJC-h_qk%-0EU)$<spx;ox&7Y;Mgnd0R6wKMbzGt?`i+rInxvO6qkE1(W{WmrnvuGniF~h!aU#}4-C?)6jBJ{;Gf;;ID&^QDmdPmkSVAFoHzB#bP0woQ+0>R{OHx}1NRde<nU}Tc$`gdbl)WBPx&NTFt#-2J&27{ZYr;bYqX`-hZ*^rG9~PrE;M}5SRIWmlaY|NvRiT8@og{NgLDvtz@+%3J3*d1<mpIUmal5Kz`gJ-lAIL5C*$f}m40g8j)OCqF)w{M;c4Xh?SPjjeOMLoP0LAi!kS#W6vDw%5*A`82c`M>DEK7v+pF2ha~(P49d|($SFBS#eH~ltn?aw?m_6HzEf0On51s7P_``{!1$OeIBe3TyQDGsl$IOD`t>c~qyNnjFZl^v0kHJn=y`a+*oew!d-<r;7y0GU3sE7_&!zFMrm#k?%u&&DC?*26@9c)jIi7uqawNB!g72O<h9T-B17-|u;7e;I1bF8@OTvRF+Rv6k`lM$jb3t;Zlc>Sz2v^_4bvpWjR(hEA{@}i2g<f8oE=rwa4>T|`~;{Q*cI1eg+t6^TPn`y3<a>6}|GZkl(ziEtG*mMpety8Fs%YHtc?ftp4nM*UZbPB26Eg_MO?w`T-o1@(bW3K+|H2Ya7<@hHP`AzEP>3Wj-S$0Z8#*xN7QZH!+sP(RZ&9PgSn+o2UoDH==OuV5+ni3c5tJTPn7WPsSqd9t;!97N2_>fBp!=IZfJ%11N6_nlangTi;kpXz5BOJ<DjA{&{6VL^l2AIUQf3Uy~jM0^Lb!>2Clbi)#>Y60Cd5}615r=~Xk|PJNbq~vlTFXK2Y1L??#_d+4w{}{N(4RWTsSCH48zW>=^`3O%$WLXP-Z(6%NQUMGPG~t6Md+z^3@0qvRo^T_O|eIBQNqr_H4^2U>unDhKk4AB=_PCdMyD`8eV$O>4&LKAT9hqM@AAhykyU~9`A3*z0_mE`n6NfMVSx5gOq<~C96!|jYid8mK<9ukiHaBP#A-EdtOqiWSIF3(ZUTf;iHp@P9eS6k(zUtVHJ6CUUX_R|4xITzh`uvv5d$4apmP}5%PsNu&LUI?mr12BxXDRo+1>mas-QwOPJo8}_Jk>(wT36B!jKC>aJI`i#$2+(pisZ1Bm-oKvMSSw6V$nnQApbiR6dO^xI%s)`t>8?WbHRF>f%g#D6}hC$5z%+^>G+{@q#REvJG%(@+YY2{}d2pEtHSk52b+RcRQRJy3S*lj8>0XP0Vuyp{Sd4sg%)5iPeIC&NsPVcu7;zeSfaOuo2>siV{rDr=c@A!Y$!c@!>`(=1QcnSAOjOl}@?w+z;;DiP#`0Tce`^N@w5Z68q{2@v(b6(B7D-2Ye*H9HI)W$69;n(eY==AAa?$K0cM2A|KDyeTK7?B|zK2w$c&3tcI(D;=Ctm#X~L=nO#n2Qk(+a`<Hz6g>|Mj`&Ow_r#T$>Euz0d6km6LkCbkppUjuLRqrNVb01k9kNk$IxSC^k8=ql;vdP5m04W*Ux8Un-%3#+bTT;$XL%pwH!EzHq6t=puZ#ThjwU9ztZl3OvW-n;<Mz+};SEb61Kr37kxPJqfR@EA9DrRjHFrLx?kfSWVHmDRCM<K?L9@{h1&pJhs>~4HCZ97gGLg5Xzg625Sgmrz{5if7F)+uma!NNN@2rnG@R9tqn&wEG(5cygmSi?3B{Ffe!cte%Y@^IPr(m};`uCS#*D|efE7@D86z0F2b-fs9`=hPb+RlYBTm7A&F;sGF86CT72J;YKgo}cws*aD|5SJ$e+0`p&L=BCCL(mGLIHU_`@+c_B)6+$|=i~!IxIqac=UVajvwDWH@mhS)6bqM42kaNwbBI6thOrRn#){<3r@Qz#vvn+Q@MD@-wcFwtr{2RAL8LcgmYg2QRdH&}i*t^JLgGSl0sJc#<GBS-FfaO_b@y28a_R8q&g(GI&-TOrPB+&~$?8YaHe`ni{9&Ns{iPOLgJr%s(+~m_5RFYq>g^nVtO!)MbjJz+S6N}3a6v22a`Q{jaPU%8*uFn4<8nA!+(5~xOL{;_>13L48xB?n9F8c;!yexD<eo^xROK5tq)7M=Zr<?4&0qlD}(ZQfjQoFkF)t*{*e@Q-Qc{)WAfVptZ_uP*rNt`WDUs%BMt2hT_=GH_(BxWBLt|iY3LvgW7M}E@-SkF!ayA+xm^%YN+MGA%SuQuZG`;h^5F$>XUws$c$lFYWPj;ho=brJ-rH<}RS|JqYjjI$S6jUo09;*(8HTkIMsU}o~knbTjoyygqJ96x>XGOW0D%#?3ZLmgM0hN|l%y7AYXS{y0f5o<1bwK&Dm0(Yd|8W_Wq5-TW-UWC=FlCS#uuXU{B_ou#BU&4}M=OJEYolZ99!)y!?-m-1gTDrd`0Jbw`XGM@ik2bd4DS(?WC(x5U+<v{0>RDuB|5|^tHs>r7Iy+H2n$;!v3Z)UAu9z|1bg@dz#BMejWDdDSRy09=eo|2A@%Fnx%Soh#vOi~QYo`9a0-&MgS$SD(5#Tq3K}Z-{r*uVToH=o4Rn6T^-eix%#?8N&X=&$RCXmadD%V1UK$T;ouf)kR&y?CBv&bfnGRH{-7!-}z)@VlsQvQ@Q1n)8hr6%9R4^FrMZoXWzQn}Dk0NR>@MVW*wf$6FiW0E0eoI+5CilVw@B3hIL2v<zeLOoy2RZ{XMTZlG1y^ph5I16jiU8aZ`Tu8+aOUMM786o;7Tz@6nCoNxq^yOEWGgVts&D~fE)B;}HBNA^zix_635TbppyB>DB?qw8sOmo#VOw@wCD{eNwcuQLsdzX#Qo&v5T^;7<EHA7V0uQq|liL>@g1>~;~uWMs1d#hHT)6pG)3|&r_-J^^W>;)-!GuN3M?Hq18i51EUpOwuHbesG=a%k^AIDVh1(Khrww==Rqqu^W?w>V={%R7Ufv<Mo5o+jeEo)?+$({82xOaHU2Qzjd-r_|(7_+?1d=r!}vVtcTmJ6eTk%U-E~V}A_01vW$dRU{Dm!S2Xrc*MnWDLAOx`IfzX_S7q>X+Q;)#AW}4G_KGH@gRE_Cv;i_^AUZi34XI{vNG`Xq^|=v?ctlul7#3-d`roXUdc?xemi!FK`9VQkvl;3;=hUa0eeDoP>3Usk^y9jO-A_~9&aUh6NBk6{Y^8rbMl7(-mgZwic*N5DuO<Frs_?85eKq`G^lE=@Wa~_238{{4Ycmh`S-^^l89~r!g+0~_K@6VaH!6qB8Ti_0kq2;Y85JFPP0?N{%z?$Xan?!o3bO{OLu(&$hA%9*YI;c#L806TSIS*X28<U7>{=ouU7(EUrw<6d|oN2uN1DFkGx)LF>;P%#u?=&qo@dNydAq4hWIE6LBQ?hddHrbVY5bKJc`b80}77me_vTU6bC~L1*U?zF4CWbBhG3lq(i+M$J5e;Ft1GD_rEK8Jow7q3RDqCG$5RYe;J`Rm^Gfvr!iKLFhgUN4c4Q|W=^gqH34T&WEH$?$0Evb29LnW!mc4YXB!b!;3(wkdcGHx@NiL3Sfx#^rMpO{ifXqMl0q(Or|JcYqr#Nwr3R!14ruYI;5;FWk<e2CM$hBe0n6B52r{uH$VcEr_9qK~ah2OAd!DQm*iBG~Wd1HQm)XudF;X^2(y)Z;mBM1=E5*49GBt{zr@F~3p?{KS{GwktcrckQZ@Q(#GNNxJ%*<vp%*z{5EqxXdK>~_xI<Z51JI*t)3I_M;{Lo#B-JdKGcdZ0O^rT05;rUhyf*O7A_@z>IlNR28$+J{<us7KI;^X-dP000peWPY!Dt(e>t9iY)`3^WY*Q{6j&219<`DA}oU;UHmLUGJ=G3h$NjC$T{H7=b}#=PVfq-(1*gCN@rd9Ii}YrWQe83KAekTXd9BUdOFd02A^9%;IFaYG+3NmD|kX8;{Vi{*@q*2h$fo%soFTpA5|KofZu=FIW2XFWwqfO@y-_h8hW8|5VKc%+06T!sq`n{!7lGSTM*-Uion4D2VruM%DET33aQ4Gg&pAAo<1uXZdo$yjZyHIkXJdrHIHC}=$l2{&Wz^Z$g^e3gZ(cK70RyVpkz&y-5Kt+XR<JmFc1FR$yj06_mOo(oaD*s-3PWd3)AAp5R~IQm<e9z$@0U8!caFR$H0TiNx}P}F$zg5bh}kyatD3<O+7nO`;IR&x*Ebwj*ljXUy!94b=)Y7WSM$f@7iNovnBeQU7K9V^gsm9Ktz@_gZHiy4fM^jRdos1u<jeb}LQ2RTw3WK9*_vjFI0AA|;}6cQ{QL@GloxYemk#Kj*66SUp8)t7qUAXkMOP0jO!SS1@Fb$cE5jZrxS<G7ZrGe1PtrTue%#Y)Hf{SbUjQ<^7Bz0=Y1^L_Tm@357H$D0VE)E;Iwdc4+j2@CNgD6`Zt5q(UtDYa&c2i|+Ux7;wRU3Kp5n=9T#1=+!PSEksR&^sT+zY=qqm4W|*JM)3=+~9BMCZH864O|TzDEF;79iohVI+%3DtmJ$_$yqpSJW-~r#>eJsOW@no`*Tq+zric+D_pAF^Cp*}T(TIZA}#)XW6VyR1`XK4!hjd1euUZU49Ajdkp`&KdnW=~+*mIm$kq5^0z+Ha9r6@FK>H6$W<b0g3)7!AI(S@y_N#(1&TFaU;rf6Yk{Iui=f}{M-63dT&&F`t*scqS`Fn@$cCd%E?kv92Wy`M%A5cOeR{&2{;P@;4wRED3z1~?EM1s1mv=kn{2B}HD(!jPnmMs2kZ?J*M(QcxDYo!EFLLPA+k;)zLK%t4Y7yx<;h@grsHfE+{<cC;p2f0LieDvo%LqgKPcWF|*n6La4{i!?IJM}^Y!TU&Vdof{7Cw<p*HwR9|%axBpZOnp*kU$G_hs+&Ro0+<N5Cr@MgiS73kQ7?5VqkHC#@4sM8phC`0D47a6PPsA?=GGfbDS<kl%X3Br-zq;G=9{d|G*l*x@F(1)*s!oFer<<+rJB#Axf_>HB^?u#PusI)VXC)TEQpsrC~^6dt5wQvniy7oZDX8;7R!n=p^q*a~zn*?ib<^xOlcO0tH<vWzi)&R{n2Dr-){S1EZ|D+~M7~QcjDonJmLU&7~k9)O7<O8wEJdj;FkTTyF~9xc1IjI#Tc%-iu?c^=sKOnZcI$q8p}V03^q~D6C6SGkH{-;JrE4R^LEAB=s=*=c;-?W9NRuMQ=$f8FeVYi6q#|hr4mfLVkfwx09Yrhwt7-piZ1#usv}`5V|1(QGv<r*T6zK0Z&-49c{;rQokN@3)tWmvl9Nz^NFw6%gb0+6}`8^AWg_Y2cR#FZaa}+$aI*k=$mD65%Jd95Lv&(oA;<tT>AwCZ`gbgdW^PWES&_gy3&I7$DnbuC#)1kz0m&93tU4zQC&6Toc9sgl#5&pK`rUfH4!YA#cA$^Obybs_S`3H(2SGL1)yyRp~U6E4uQ(wjm9Pf@ur*$a_jgda!uH4_cOmy(z@Y0FDq;CSl;M3-IxmrA?Sd@+|*-xcd>O;>SLi$UohCpmDoW!g+!`KMF+?2edjJLm31JgRuxzOQhiglL0iMf+0-9=8*fdD6l7c<QO5w+q$5OZXsSC}Q=rLouXuVRXHV&Y2w8DM=D>eMujVKp%~RGqd=IE%T4mxWp6$c=atTIBiS8pXADrJt_PO7!pd6P;f7aQAK}gNzIpykEUgbGG^LaWFtJ#2%nllf$4bH&kq}snlxSv=$4hpwlYdn%?b~OW!n>V_L4?Pqcbi({Kwq&QsMqZ97Ozq;a#!nn!3OXF?AcpL@rWMpZlL7q3VN{0@moAr16P1iBpghA)OS`UG+cNPS6OVO1#!QIg?%w$)L_IM{)>P83@Q<U^7z?I^hopEp3&nQ^)^EMVzmn^HDN8t7()??V`k+m>V0Z4bQ#Ff-n$Rf>jJ@-C4t<s6)-Wb=TfSKVbr&Xq|3YBm1hg-0A)S4gtox8F!NFZrGW>uY$JpOng`M!iZ$2Tv;wAyM;aHu2WyM|02uB3*OhdZBL}Y$D#uSG8zPkfxG*lY72=T{OLAfF%Cq^vjYD_DVg=wlU_b1=ha(8g1eVHUQ8|>=*k7&jw{m?ysZQRc3ZJQgmN2+vz?PzO4MEPUZV@~>**AuDZ)JuwwmzR4~8`;ZC6ff;d5l^dbyUFw`Mn2qpj-li9*GXrlwp7iJpnqln_V46Js>6XBa`YrJs?)ZK)DBZmIC&F)a^K@UgaJJo4zGn%%KEG^6-6aaWo+S;Y8TgrRM9=jRE{)JG$=hfFsBmC23h)@f*qGjdQ<nv7cXJ#3$$u@(2zZP_b-`N*9pDm^r5jQO-6&3kRVUga>>_^lA=JSH6`&dK`oxK3;pu?pl>65k|i^Z00xzM%kQuj_v8y7RFlCzM@i4HZWr~9D$0K#ehV3iS21B&jhT#Rui}ALgytW;S^wRNRs`CnR~5P972>rN$^UF$S>oOp^LVSAksfzXuL2C@dgyFIhq@S{&ow8MTPh9~8(8oX*j~K;gp#B~K{>%rKx@0FUM-;Ua%|aUZxVmo2dM)N5oMAQvuerNjix^R@^n?!CXXgK$k%|*A1*#i90hA;#u5#nKirqO4ABd3&$5mxzxZ81wnBvd@B($rkJ+cMrm~&lgQ5^fy~V#BE!vOCrQpuLFhZ{oDH+G_AebeQzp2^42E1>P+(f^??bc^!X(M0$?u)HJK0kLQiC$MF$G6+2ImLG^s5jG4vbDD*)v`LVm&gT=4Y$<XAwp<QA8|zC)eZBL($7W!AL8T7wwzTa$Oo6|H-Swr#;kquYwfb(1)~eYw-1@C!ez#3R?HKjML0I`tm4KmK*h}D^(zS9fo!cqPhvshKqikf4qsI^b@HhB=hdhS%YO8qVuMdSonJkPdhWi^zvI*WG5l7qF4?=n+}jf>LaT^;LBP|{S8jb!T=e_OY2WU@UuFdsmAZNP$eecpO!~9zec=p`;AIlNSVwS{h7qK<C`Rhpa!EohrrR07^}^>gV;sw?=y$`98AMu0_#@G3Uul1?=S{tKTFegL9^gZX+FCaHf$)c6uWK0>XXPf~VnlgG>AxbP9=?Itd1Bk_YS%=J(_j=>ef{er1q!bvS-q8s!3@Tjq{A#QYo*(BN0|DSRu!4OpJ^cd8CK}66QyDFXuN4Vxl(;SrETJPllE}-YLlOZ;v`_Mj%nYE8k(n)D@$<ezrFqzv8O519kqZD^;9#1PR6_U^)ots7$bh65X?~ZsRN-a$-t=A6CzT+o}GtyCDGsPwO;U~cD759*}f+eJFU3S9p|@Ecl$gaSy8JT#d)|D{ioVB!MY0y&5mpx3?+WtLli!|n!V>_Hto7`v5`v(!ZC(KjE`v(eF%+<P7;3#9eEzTfc*=$V$opgMbK%rew8Odtd|4|!%CPX<6mMjhi|f=U`%irxjn6Ba3Q3ONfMr|aky91C&)ZX!N?H*><sOGiFal-<2c8w4tf-9eDJDyIegM)?3nUDoaoRvs9+3S6A&=|q{9#v-`vZ(^R((S3)dLqhBd%+gHudxX{5t9@$$OpzGvC!)J?DLpc0|7Rarb;2F@P;@wNl^l}tO&lQsEWqi=E<;B3JYA=XkRqz=9?WNS51;1$9eXET2`kh)Ln(4MFRjc%5!!1}B;gYDn5^zUwy|Ex&IvgqQ!4=;d!FAGC(Hm>H8_$tc`T_Nb~<xusEMyPT_-oqMtaWz`@aTnJO2hT^B(LcIToVOu9S~0r>W;t(6!H2o`z{>Lqdm3gSn{}mf&41TXF`jLPpD#|9MAME25>-A_SAf@I;L_3PeO9JGcdjHp<{#6sC#?>Xhw87Vs4W>qDhKnNE;alUE@lce0mab?$<|0gnfVz<Jq>92@2}h;IHw&Gow14bFVp3pdi4gbpP2DteiZ}Ho!ritn$Dp`h7hj0{D(ad3TiokMO-?-x}w+p7!UVIw4<P~Tk1auHW$UQ&CaM6Ti6dfDoHRf;I2JpQ^VEDdqgW=DCl}IR(o2hY3kBgyG1y3Za=crm-q;JhablU$*;o~%^qiR1oRhm;;1Dp_d%IFzEqO6bS$cWqa(Li7{8gWpwe1yAJiFlJ1S0_6>3j~V?<t)Tbut8Xc$5`(9zGOT$_RJS-zg(B%5^-5~FB`s#rbzp&$RoH$$cpv&>9bD}CF@>)IIq6Yag*+1O6p(~=)imQ#?}E9x3Wk@5tMD)@_^9!|ID4nQr`_#0(1_D7$q63iXI-te8<&YEL-KShw`!icseZrf6}I%`aw;P3Yv7^jXS>yFZ+iaw4O<JQC#%lIzfqW}kX?cqZ|2+W`%C(uy0mXEu$o63=U;xp+rlH|`_pa+K|x^}jBqOI9cL@*x^KuYunwuYxZs@<A%u-`#nJb7;9kcJ7^QH`eUW;yDe-h^?#@UEVIF9Oq{<iY!=HzOofi>cpCs=vvpA91yoeX{6h&K|7ZMQmZ`lj}j&e;|o#jEERJdX)%lZibl5!p{H2WI10V7z!Rof}72OT!-ziC0n`Zn96^8{#nv7xm=?;Gia09Y@wH`JCup#MSVQ@G^jmHwK3l8;t7^Y!HFaIt&&$7X`1C(fa7ufr*2Jbsw`tBrAP&V`}!Y`0;lpYx-jnW3E*X6Xm)hUJ&@}M28ZW$2o$A`frd@dY73AH+fDcSH`z+7Z{b}}4@4vs#?8Gzv7R8SdYRzp(qV-8$`(*`cf`7B4uo4GE!$<UvmO^oj7T(II%ZE^=azOOi6H$AEgYlu>wBzf7)gIBEqL{x8_*0EZZOy@>2S^1X57xhOd-Vw{+xO$3#@9FO2`P(83hA!aqHe42MiTlwEtIFTD@N4Kp_jTP=RNbV=u@-n^8d~4sew6^$WOeS_++$!bYN$fW2oOICg=PB@O48JCFJbHKj@>LQ_VB%HFWhdN-c*s4Ki!PC&`po<O)bTHH!A)ukE+c`ylZL1&-#zu@+^%=?guw97AuNFKqxpK!9Kh4|I&mf)J?@__|ndkhYNt(~tAY`R<S-rWlzyTh#wt_iE)6d)MkPvHz%FCy52V;w9N&nR20nGS>{LvoA|)T8luO1Fgu4T{s*ts+VDF93PA{`Qaon5-n?AU*?^rdnBXFtJ&>zxxtA3Rq9*(u+CJ&pq2&#rO5b)eX>>H9`f{SZmjM<9)qOQ>y)!Bsv^6S3pun{58H#<!oCVQ3JR3s}VYJ-Tr6NDt(*wnz3vmMr99&HUq*hLb}|_f%0YX`@YFHM>+MIce=0uSq6(uViH55K!#NXB*((ewIp=TfAQoTR6D#@eMne+h?~{hLGt;%Y-A7}4t*1m?@;>OsI{q<*4mm<)?}0^OnJi)5ZS~qj~0~lH?ql}jR8$Z_I_}3vIC^bGw?b_4>ZSqMb?6V`Ha;^l*DiXTXZCGqw+}Z@jJy#@2V{P$UVi%P){uW^a2)GJT-9O5*!7!-+g4T%^~Sj<aEv*gPx4`(+evaopEn!$0$^l6r!#PgcL5hWwj>+gW&H+Wt5K00@cpX^W+&JJ=EYI-&m>{==Y7Z{l-Osc};hG3R&!0*v5R6l4J9SO<hp1$%D5T_ypUU@Nn6Q+F^dPv^eS-Lao=G!z{WgLQ|w&@lQ|89^(6=emMmdvHw8Z)~T@);O?45YIS9lee_G0(MeiA-FVwcpZ<{t0G3$QKXo0{2mpfs@RcE|&g908cVYOpU9>AG5rF6v<&lb9tP^$doBtz=C0R*?CcDzsvFt*gcWKz5@$tnZlp-dhkO0zmv#_AnPY4_77IcNmhjI!alH(qxVgZF=Q15%q7*)$s4neezHZYre(nqSSpY9KMpfISgCRn;|<*&nFhHgqSyyG1+{Ha?8C!Z>4;XZO`%vOy|Do!VYoHs^v9H0Yko2v*&CvU*o^$&447mxoL|MgElFLSvg-vm2Q8eqImzY4>n_-bIFVANkvO!yC;a8BeF*@T}z4hg_u)~jueig7`;#%MwqCcJ;nq6r8AvcBk)F%i6(G&Ce43UV6Fevf{8{LC*^TJ9=)IKc30nPUXrEs(5Wk;dnm=aHRs=!Czv6pDPXD2rV(QB2{8Q{$y$I%WRoo9n$muO0T=O&jjm8r(`Hi4ga|$=V-Ef2-H+>4^@tacr=Heq7-D#F$VeGICWn!$K6|L%_$-JcJ*Ufblz`Ifje)1RPacKW;Fv|L2ey9TEVhp>c{Db_WhWagIds*nQ*m`1*n;(-!_NE<Q7<n=DN|fC>y;EJ*58pWZhg?Vra<L@WQJ1ZAidNAA*K!K8gjj@~?LhCmBUDC4%myWNx>vdie2iKf^R`hN%zcNPvC?*8DAu~Z&i)>tLRN>C!=;R)B%GTdKs$w`^L+e3DC2uI91_rDXMvol0_1>X{(<}%85_Vp?lEAr~Xl5Lv}t=qFmBYUPId?<6kZ=0n@f&IA1*9HKv1-bb&{GM_OD|<aHS#h!t3pphpW-Qp#*vz`l(d1NK<_x3$6=W##-)FxnI76w5mNny>Sl+q*qo8059OF~DQl=3YeQYZ8t{*^LF(@~FNx3~WsU*|_e8iu_nCCXY9Meo?CH1Ea-k~s(KwtWm(~NuHac6F@$+%3XNNo~gpDT4yERL+3{|L#IjnIL1;%f7+<z1D)!j$Y5s%<9CAh?3e0iY2a4D;;XAGwPMs?3AN)LSDq1oP7?O*aO<=nA;JFh)o-i!K>dB$?u*F%h{Ql<9b?ZN3S^>Rd>SKPakJ?;*!~b-f|{)m35o`tdJ`7||tYU({+oHQ)(#g1=ku*mvEer7kTKhbv|$v0iXW(k<53m>0ZUSawsr(zoy;=F?2KS@2ZDlqMWPVeAHnBshqm4oLhoJj~p!*u!f6mm35Rsl<7Y`z#ytZnefHpnzJTytFh+3qO^UqE6U)qXEH@jl1D+?tTzFy}^`TK5gUwB|+#dMW{t3|MrY|*ET~<=?0U69Pi%Rvi&lMlqbiDO4vPhE?HKK2wEonO<kxnMJ$NvmLi?qek1sGL~*y-Taa8@D<HoZ6|LDeK43_XQGF0SN@A`wz;Y#s-_&H-!mdsH3Gb*+?#WI>6RDpy%*4)Evgy(W7@)G^8oHv#f38m9vcr0>fRxEbJsTo$zQgF{6qQDdwO|Lj!So_x9X=m4>Da;ON>1Q}J{)dSr{}YJm?UoPW!po|pjE)N$>1OT{J0-dg)LI@;p$*^SAO`HX%^%4`a8N~d9*w0+P&~gx-7wL@OeocA?fGOyuz{_e~9n>Q|Dgw<yq%lECxmJooQQ<5=ix4;O@LfmJ3r#R_~WhLZgsN^jLAw&^WH22W?9;O9!#=fI7*rGO}LFAzB5D{_+g`^G(#}ktpVn0Q!m0*^b8$94GUItO@h*09O$5Oho5^rG0^RsC9Q?fxNQ~hWE?vy2+ZnBIeF5Lj}7dkNKz?_K^JqMDpDH5+>Yy;f54ED$Fsjiv$>PiO4A$Cr6PsynMaV2_Tedy%N3t`>;scMv*bxv|C<zrPBE|mUheDN^AhZBs<U96Gju;Zdc{Mz_(If(a5JdJqJiV7jmdI{1B+_bvdsc3t9LD+tbA?IvZXNbKbVZjUsQ66-keBwkYT*y+DjAR5Nyz6MyO8Ia|k|%I3t6iuh#zJE>#uz+x!bISdy#D1$j#A;Ni;I|ZN7$axI!bFfS|U|5p2jDE}2?Up0SkI$Y0l;8G~d%KoXoZ>1&gYA?&APe||mE!YjA+VRcs09@=?eSxgfrETxq7dL3Tgq3~J$rzg(iJ>#Ku5c1s0VKG%ln{YE+j5w&Kl%lhfO+NXB5+U%9DV2lzm|mY?j2L!6MK3&N0AKf%J`trFw&d1}4ykLiTodIlv!a+N6ZAFV3el{Ae~u^^dRjy#T;%ge(}Qzxd2UmVl+3y<BT7*47o+G&ysPYs*o}8ZMyfipS|N?W^<58S<<bNJY7{H`2WgF1%Nz$<RSZEeODFJ310FQmw+?GnAgn!z}pht;4NgO*6|~#8g(Y+NL=?>iUix^lgVGx+W`%mxfcvj|Rs{^K6xI$^>n)WEYJcM+fxw&c6L#_^GzCRplJ8@<v672`ZX9IiJ19T<)ucIQ<MbOM@)JQAl9x=@)C3eZ`3ut`6JLvo{1#XYj6X-C4=ZkqV{79E$RC#0>1$PH=)ZO+WloTi3<M`9ZhHcM1C#rFyI{7^tVrBrq6mmrUfUY-ii#G!o1;m;K?_m+JE^*=6+Zk0d4g*nf_w5^DA0;Wg#|v0*V@6X3Fkw(Hp-ex!mJGn`kS9PJWUGmDsK%+O-_J%WfdNS&cJzC05Dro6l_699I!K`Xp>*Be#Nt!Vu5kN#VoP@Uhm>wzg0WmL)(odz&!y2N(HbS2q|4mCA!{2)})6fKS)yNATK*aIG!e27L^+t1{{zksD#e!(Wcdm?jnXdIh8r6&0hLoL}8q}j&9e1d}$B#37|j{|s-DdMNO($Pqaq{Y^G$u*MlxoGGxx;`f-h^}pkh@mZmd<aMT$YNi3{{lIvG!ae{n4E3y9jvy8ECSjS<kb`c9p?76RiPZk)^;#%hWJ;)U~4BMau_q6avp~!Bm7I3e+aMWjgzY0<tAIGFoNR|Zw={I$1y=U#%O;9$o{7f>~Cqj2&y8nj_m7qRwkNX*tTqpg|?&)k-_kNdFW0BHZaWexyXR_4(m9M_+K{_x$0q_yC?lRVxYsujCrL!yM%QuvjBgkLmFmLaFF`s$NKBDl_9e#A9{GY;oZM$mjP@P2@m;%nk?>5*~|S#DhymDb0iWKh0BnLEbZZWXGf@`qglMn39@R1F}ilq8=0`UCwOfGiDWH_)-=|ly@@GhC`Ew|#9~808Y3f=XCmTK0bCDu_#f4TZp&Oy*av%QfIrTKI_W~n*qh1;b)+Qn*)AKijS&*|E6yjrsnd(l(EBEeu)Mh;tpmH04+*4Qjsr>djX;i!&41+ZC~OX!ZT0x4vNd0`w{8WtJs`433s%-~ZZ5XH6mpCB-9Yovy48Npwlmk7DXUI+Xjq^YL)C}Zi-W60vQ;f!%KUGlTrul9X*>_5-Hl~oEbq4iI+Vk|s?stb%%?d+AmN-Tgr2TuE_W)i%+#I7Lm(iB#n<zyRjQBB8w9tafy6((Y@CSTM5ZFWrv9=_gA+U)-hmur@j`q$KkEpkv6D?fv#vm%BCwO93dMN(e_<AK{<Q>IGi-@25#UHJ51|$JPOj3b)En4jSnJAlyypM5VPj4Y=%dsjsen=Qk@dX1f2_TOvI!hpNRXo5@J2%m7=*-R^~B4h!?ut`Lb>3Ab3zDdl9QZf%hLaTVsfz=+&g?$Y=|w-#gmfRU6K)hZb7VI)z_4vaXOzV^CBaP;}Ov<0XQL}jkb%*0Umq^yL>TrWS{q2pw$wb8ilHKCDcB9mtBBLsHHUKS&?uu<{k}Wi}fc@^0yGToUds{Tg6Ey?A!o!&-WvSEFR^|Mc5%z1NDcAT$n`q1>{d<Tw%<4E2TRPrUKHVelynx0p{B_B)enqebw2>1ZBgKl#FQDwH51M)WsfbsT{rqi)4m%%-96&yz1@yzm(7txY!Dr>4{;DX-7Ti7OwHtysEJZ)b>!wsn1vQmnQShiLB5mmo7IDEjaZ62}VCyOY)+e-??F+mZzDF?gOu9d}d2{q~W-kqEH0vfT~?j*>!rLtB+5od0)Y=D+3=TK{(w1y4_Mpxk@|I4ZRH)2ZyJ<d%a{>DnuE`teWiqO=zGmiJtr9SA%u(swJHKS)qrp0MKX0nt}WVj`q{<FRq+zz#T6IRxRVG=Q%0=&4HH7JBRe_W4Np6VKBCa3d*G<QUTc(4#wAi@b@AsolAbLVDHxztrZ-p+&W?vsCT)(6b^q%q8Lx+g=ZDEI&U?q-^}1k70H@wNYAcbCoS8dZ!8SE)G}+x#MXv?q)YPfFh(Jz`XR)WB4HO7=<KAPIB^Vu@qT}Z-_3n^U69Tk+u_O3)vd-S1P5{LoCgsswUi;pPy~e2aWEqXBoZU(;nFsA`v`nO)8TDJ^Yo~^kc>_XZ>4J9Tmwy}8K`fIEma2znS@kON2Dp%)P8+!5W_*x{`n)49+r-k)sMo5Q7cCW+DQmL(m<KW9Z^Uu0?Z;BIFh3OFdm=|LXTVcCD1-T1J3zp8YyN9;&JyOtMEyOp=R*rM>EPf+-8@Qga3)IB|Y&<J^T@7f)PIx_X*dH6q(?2Kif-!ZU-s8yCGA<eq`Y0u4{|XT5d~sajR6g4W+^m%9`1oDdJrrcF;yMxUNePyBCaS4|I!BShqiV2o>RLrd#IYck~wxpp7<+lezAu<}eK`s0R!gm*>~QV}g}AnoF03g0e%w*8MOz2}aWC&lRCxx>Pjiv#w~^fE|=U5FC_hVJgcYUV3jaI~p%SRpll)b7RVy)DX{?C(^`J;<Tw?Qi*EN1I7#BOuHX4%{8IXcez}<2d+!fkqFx!=u&Q^(&Zjh{TToUdm?m}o#gSk#V&IZ$SaaK^7py0t5rHp4&)}-jq%GWG2+uGaFFOBPk5>44JE_99QR?e=f3XqeTFK9@z<DGX>^G*j0btp^{yP*^*IrX_TUGuNRJD=$@N8>iyPyD0MBJ)gX^@aNU`~Evu(`vJb8iDj*upH6)QoW6C{{XHuK#Y_J@+oQx&*z+H{figO4^LKstGd&y}syVWT!yd(1`qE?R`CgVE`L;;ME<V|>NdLWfF6Jnrf-)%Xz3aES}#(sla$tZh5%*11{bwEsSAj14pK3>jKBLiZZWhEz}|M9Dk4$^+$ng3sjbdVlzZr^ZGUyZ?^cdAtu1YZu32UyYqSJYlDI7cS%>HAl1qJiV7Tb<Cf%g<C<c&wj7`qa^V+VeSK=wmBnmN<>Kyu9jBQQu6itR2H+@#Heik-a;AY#Z)##uSI-Fxs0S)NN9#PS!Pp)h>1u5CryB?Og$7=jSp+lS%aK%@c=)&&4@6bYa;zo_Ih}BDcS4odG8%;v(mgWKg9K-R_j)l7?6ZNkq2JDIlZ_=iXmVh={5vRLp=j(aD?-3AlaP=>+*^82R>oM6l3*+#H@Fb9Od78M{tBg6hoq5psf&7)po*(#|2r9r_M1QG@a+H&!tgdts&IYP=8`ZOjPG))yg5qGrPM~uxh=gWaGVLiaa}amJ$PlQ0%BiLvoR3Jt+Pw7srP)t|#SSdgbiaCJS8@o>zX;y<Rg2M(^(1Nu){hcSnI1L_d2_)Zkdc6m3RnQcd@-8&@;0<)!rjStSB7qJ%b`%c3aQH3J5MP>9p$+IWL}Ul-6~doC*L;oU#w-IkrRFi%7)%>G*&DOpdAi445hbxGyM*lnMTK~UT_{^^cqBQ%5n#uooGY?tyb^Av$3j_&Wt6~$UYQ{?;iy%I(q{m&ZUcMpRWp)QvQu+ZQBw5pz?c79-^Ct9VfMZ(lfcH*T5(v+<nbL&8HrQz-LP{%`yy{^wK01X4gi9_ujzB7=QU;_R@TJ_AK)e*{7q2hAqSFwJZH<S^t><S?!_sOZrZcs>i@%GrjN92aBN?ZoF>}tdGb$A`tQ(;Ej=^AtlbfBO$HSudyZRnuzfd*9nU2IFPvx1kbx5mN-N$zC>%6sn-Ml>4a!!=57!tn};@E`XO#&;E{+4`t8USEtCVsYY<iTO^6^%I!_X3h%ugbx8j9n!XEop+g$MtQvaJ%dr_N8BGCLr~{9&!(}4mhr5)nUd#fMT5j}!5^C4w!p(wIW=7Tb)wDkc<KyiwqjO}K#eP0|8l_P-bMkot*&H#gu|s-wA=4n<l$JT%7z)S8$`TIhH-|_a3g)FNP&0B{gbr6!KqP=&^1#aGAT?6Z)jMkws15Un_(p?6Y2Tl@6X7=8EiIIt4qA}mJLz#FZRM%nYnC`F>;>LSBRAJ*`#jBB2b*|q9_rr@ii9g{#RMrt+$6qYZZST@{;wx4SEu&4k-t(CeWq>tNSP!dOrjbmW#iFNSw=M5o+XE8Vs*7u*NwZYQx`ZKlG-6J{%uJM$fR)A+-Ym)z-ic1^%{m6d)llu#Y%Fw+|&@M_jM|cPszIM4tK{TH9d0K$akoP}>b2{qA7O4Sbq+HAXx^U!gn@Q|14a#4mdDX8mK?-aj$rZNH+R>7GxBYieSj-KtYf`-uU7y%!_VH7S6snLCvt1^U$&G+|%TIp&HEN<Xnu4CRs@5|z?~q%yVYBUjb5Hm^9Wi7m%ec_q%g_`l?1wGTrAvj}+Lu6ugiky(oI^iXomqKtUf-%^4Q$}s}GvMS-;S`+kn$S(7`fnx@(ihWhq`!Gvf@$jG9ZWGq4#P?Y8s4y~cn<a1Ftt^%X^?r(R<jM<MQ@a43F6(6kWN8AI(N&Rw3lgDK-ONDtX@F85O#!>^PI@X|^Yr_2((gx-uEy*m@%jMq%kIsI_)0pAo<4~B&6%OX>SbHHJW@ZYP}_VgqujYYdlay;I#?$$YFi=>gA16|gYJK~tJO@nC#~0=o|v88xE@46RIK#m<9ln9!s-@ZRlBcm`@<W%g;O^ZLu?Nf4G$r}i4LeS2Q*%bnw?!sZ-N?F@d|0+k;y=a#cu^=TZZec#Jd9q1j=A0KL*^H;MjRx0S-+D2O8PurAZD0igu+-Bs%0>7~;tg*1Nc)#!ErkP4PPyKXBHTVs%@$Es>?yi3OT-tQ-Th*dHXkII&qAv<Z^=k4v6s%PFVK2ym^jlN@VQ0fhN>(gO0W&jp`my!4AKF6+<yGGR<FnM@x~uloB)JNg4$iH3gSwuaxN*Pr=RhyoE88seOJ96+}*CCV&gBEc({w<e_DduyUk`!#dDF-!nqcFp&-YXFHh0D64F*wx#;F7uN>YA4yqTCrD@$@@4GkHZ=bBCqyDg$^iebF-*b`KVxbksP3EZ83%UM|?wF1`Y6_ttrm&%>B-{s)3*YSd})}^?o0KLWUNZ?yG|whOn5+f)b05#?M8%!yly~M+5)gA3BF;6U;i=r$8EJ=8FA4{|Kxni*i;c*2crU*;jFSf?WrPS=OU#Wp|(0gKeQ(WY1i$IKaVA1UR;qF*gMlhLgYDY;sMQ-6;-COzhaX*q+Cgt$1(MygWB>x>tCYKP8S`=JW13C(jmTF>}^(bGO$Ah?dhRzG1<(w@HGr2{r0B4!t11`jX=-5yzb+f5Ga#&91<_TH;V{)zZR2NBpMd_1ch)5NgAn@wG89)Lw}jOa@Eva4EJaV@VXomFlNGB)vyNd=;AS91VVT2l-L<)uB8Iyx{)uqcf+~t+w;;eihOKn-Px7W=A<$UCE^)?WIJS>?YU*;ad&=O<J>k97sizzV7*$FzMd6+yK~sE25Xqy1@M28~P5oWSZ3S2175r8tTF5CA7Cgp639beND9QVp3wog{o>Uhp1t+<}@p(Q9tS*#bYk+0lz;}4iU({Wk>%kH_vL|ZDDrB<Qlx8IXPcYq!~IZk%I3hm-$2$s_(%~cy~$DrrGzZyy_Cl`UVJ;(2|L_W{d*-PumeYvpoWnGc9TPGC8Pj3Eg2n6%1~fE>%KbsuYuvP9;M#X&Q^k&ESn1DM!|eCBc@q`Gw(j1ke`B@oM1%$Yu4^==5>!%y<WeJEEEs)SC$W>BKkBuy(G%Wbe`gTcS~O>ZSM*R`mam;$<^meTuw-2{!lu{fKVkAZN-+41YnWPMJi{%(NWuin@f8ZWAYld7qZk0WDS1s0_{^6fQ!$V~+LG!aeO$=I^_7h5H_2vabyJARjiMWkH#r%^kpMPmVATdQwCIB`(b~8dAx8q4u4D@vt#xo<DyInkImjH+LPE^<JnXF{<J6x~s4X%;~@s-)Q@wfL?krSX4T;Uok?q6!HoJ6252a)DDX!IbJ{iHjh)mVzQMPAcfB$9LQ}1@5~^3LzNz0*4Zw<wam@Cfn$xG$9l$ZMBhK)l+?~cB(zJQH4Ip)@XMsXb7|%zql*Vxs#b9m`~<O~!+Nd-?nS_XN$~Q4bRAkKA@VWNJ%uga`|jkiy(HIlj_+LevWy%2pnJf5Ai70&r@`+}{tj|p$Yk&>1swM<u~X`%BtDlyM~|E&dz0Rf(qjAvIucsO!{l8`0ys?6Y}}}%6VdS(_{I7cVvP(|gry!Do+arBJ^AQ(AC)RqLv&+_8E;ug>3eFex)_5W$BZu164v09@0w2swT(|X_8`x=(Kn}(!tZI_X%RXH_TFh~6<Fb6=+sg)l%UWSn#U8ZK+aehezAYb2!Jl~FNiRze&`CA^6U`_f*aDi)Vs5h934siQ!VyhC-;`iO3p1JvF=soRjR>7c8G_#f{c^HsrF=6?Tr?85*?P*ygwFX(P<hr<pirx+YTqz(aQe33(iKh!X~BU@hSw<;s2>42gN!zX3hF;r@ey5KwOY#f)|e#AS%L#_BiLXV06P|y7?l3mM2U_Y`77I98bXY`xF4aBjG(0XS)(Wd}#nvJ5Fy+JS$@Be1mVbU^t!KYKjeZtU^9|G$C}e_J`-e0M3I{6Z&V3nyP2WLmS_|5TH!;bt%_HmvZjcd=Zse96{uZ&PAif1IBV`tVNv#hiYaYk^rjTAK$<Ju1zBKjVK-q#-2uTNd7&FP()(w=uo9QsqqYAtXFa0!h0P4ybRTck4Y$D>x8D}7>^fbw9IDtbybYs4T+Tygod46{BEl1pFhI6OVt!=dcU*dvPwVES?Gh=YU}1JKMn#XuB~8LGTgeF6kBArtIOL#il3aW>t&2QuO-}$r^DZinpETC0MheLY$1&h>xg;Sb7I>k$A#}lfU3Y@q$qlY0HeBe&gSCQfBGe<3s-I|gcmKp8(Qt$ljBh}oqS8aM?Pj`P1P`>=kOk3k^IW&KzLNCG?7kRDE~Sb%Mc6~_BOm!@7p82^J1{kRbhPItSpMwR(g0qL_4DBGedj)BW>(eE@|H9({DNrD%C3sX$syD;1w47ST)fT<xXQJG{>rTuup{ujLi{|aiDd2+X|D{5eysyH+2Ai7U+NdaSiKT0}FW6^h@j3*x6%=2^p}f>DJ4P-s(gy(^oHO2B}dFYr37Vyv2qzRqm|;_wSosP2!_vojeo&+5W&lylXGc@vO;MsEwQ!Le@iyj@er^JXyH|J&Z5coy5G3+tviULM96<V+k5tn-!m)r==!)_Eu!DO>Y3nyVJ~1A3(wp&E621Wvwz?Lz7S2AfksxWy+5-Y3s<>znk&6{bMO9*Nv$__SjtWEH^*qVI=giSrFz~$henm9T(^?iH{%-4UdCS#%{vlyyr{IL2C5j04?mgIae_1LA&s>aKrAonZGLzk{Pemz5Y5Og&fz+4lH!`IACxxoOsdMgxJzO@V-$qA!H}w%H8lStn!}INlq+I#ie4vNt|VM2zHQt{0OmDSNs8P(x}XW>RxTCP8%{F(IVGKFWNAyNoArSa9u73d;|i$g=uvH+Sc)GRFfD0E_Oc%<g_kK+HQ>DIOW^|q7!~fgIYCBEsmV+!fbm3by))nx;%{IW??fD7s+L)tO1Zp;}z#RF{l{w`<#nrNY{d4I}wGKlo%A+q4@~Ok_%~RC!beV7<KgLz~H;Cl$>>%8I-M?FotgO5ti!#h${t?HsmpHihvti#T5k|Ro*u<PGJQf=BVW1^s&phe4t+S&(2ucWf5plozptRuMjKXkp=Q<?*%))Wj?Snp3kd>6vgy9sHiINz$9TelnqvyAc}<hJ2A*hBW+H;>%ghXL5hy78bYSa7<ZaH_$@Mw{=7Cv{Jh^aE8CJmGx^Y^!Hw6K6^pteeuzYeBOk<g58hWtH#yffFM$g@heqB!=NiL&wjxgXk+IA^>VN-A&~w>Z?&tadwjzPGd(Z#d5yOK21{gL1_*t7OsM3bhZ6S17y0@o9eb*20-UM|^JvX)##tU(nO%A)(h)_v7xAV^m6aFb?wdUC+6J_`x;*E*&$p{@?<LG+lNY1Or4Zc_e3geElP}!PK9;8t+;Y}y{eInN`lgQ8Z5joc^B!HR;)G6Z&>l(lR4}gMh);8UOox7n)CLU(yPKJF{Z{O4&?TIQ}01CSAMhL{7kabZG!uM66BR(YlS3h%BW0NrVs#2{Ks846gQHW}Ef<!3DIr<mE(;25B%&y86ST*!oON;pXS{(|VoZ(VpT@We%2^Z(L)pw(7G=+PE)ya+qdp^&bOXyG)(4)z)-MK!j+)lw5Zyu-@bM;efx?7Yi6tMH2zxl3d8J(+X-M^hm8Ddx69gx|&h&6-9;(HHz69zhGti#YT?)3GwMJe&R4+u4WIng>m$oAo;0xtTcV_Bt30uJY<g8Onx?S;u+xTTq$#<Ynx7EiQYm5r3d3riMpb^GXV-m*eEl%<my(<*$BrBn<}$%n)SR+%JdsUt*ZF!0UEzf#FZJe1EtGWvO+D$L{cSU5}pVT7KjpE7;N9N9N=Utu2HOE2~m&ldh)QvzOI$zsu>u^Jc|%)gL8YiuHb^mHdGcqW$HOh={ihS%CchipBotE9M2>_cDMg?k*EgoEqJFw0_dNo6<VPh6P>6&!(op6$yRkQ?&9PXT_`T#_8%LcXJKlctc~TzeT1z$En4()dkW(f7%=vAM%8QaIv&qQhc@e79U95LIY(zqm)qhd3dwNa_OWS-ExvJx7A%5*-j6xI77Q>}Y@~b_VCxF9g=iuIU*<bHiM%W%%uf9wITWelQM5WJ46!Z`Yhh=};CpdTTOY!|%t)KH}u7M%=E?#vcVBl?3F4smm&aiZdaf<UGe7+k1Zq4jLJTC##@0o8%|Z_Ym46<mHpiI3*e`&I7Diyh?`8(n&vatlR+SnTbaVFD-*e(VOc*3Cvu6S2(Q}Zv@hXl{U|>Z3_JaCJiNZZf;L1bjSb--lVE=g7Uoy3_jfzv%)uZQ9FwP(Ua1$m@<4k(4PnA!9H*pKER}tl{j-YK|DvwZ{eWuK-QmVAKhzD(XADP?*WSocaH|SJ2?iq9X;oh?|`bs*qqI?ai_u|6OFT;gw}LNJm)|45$(^cZ_lQ9aKGJ}3_R4dw#tuIqVHF(e&=OLQ<#diQR)7k*nkrbxsLkekNQ}J&m?dDER6(p>_c+OPh(jbyYf%Z>Ch>kpk~}a6=_*C<62?)F6^RL=V8d)5uQ_K#HBBUSG!7qZ3m@ZP8KBEUMx*bWxc1zF*-cf!D530@A)fOwY^Bd9FRtt7Fc=PwC4^t&o((or>NUCK3T3v!HfH&63_{h>Ln&91)-jB=*3cHFGonx-2z$Fk@w!nS!aH`4joXiKxG*&9Df_6{}hQ>mVi=C|L%xEYd2d4o&)JgxR-$5=0N#)>kT3gS-PofO}Jf(tNf1$_}S{%g8hH6MT2*Q3vwzCuy(pjm%mqo#p-Ty<|F*|O4%mOzL?2P(s+SaXxm-BlwCdxc#1<SB$Y$!#sDG1zKu`nn22{_y&@J7!7OUvHO7}<ON-XE5YCpScv|hIboA+6;>nImn3$?@#7rs?%aEb>9dsnYSfNX18K*b7Uuhb-g2nY$USTir=e6A*2O|0bCVs*;Jp)S7-{!z(<Vp>|==sI)K}cea43@#a&1Y*f=0}>E6xeukQKnR*G7^saHPUfEQ%#_L_#^_S)kaRR*O%<tu=FTt0h`$DSH<@pK{N^!mJMJ?UIQ!blF#IkXW%%g*5)~5?f$z(^;c9N50x-?AA*0JFlv!)rBV$>#}gTjZpG2xviLy3+=7krYB{h+ciUlm#RHgJvvF-fr+Oge)558)#du5!fG=9CWl5KoaK_0|68$}!&q#%aJB#eH96s(5=K%uY@N2Z?JtycWDh4)@;MmR~nU`i!TT?|V7TC{MPJr?$t9*59NDTnJ)4=DZz^R9F0SQJu*+63M{z}wnncuqzZGU>eVbVU-Q*BUgqLSfxyDPwN$Nwt<<1Wm88p4QzivgN#Qz6ar;|i}jM+ei@Y&nl`y~5TNZr=r4LlF1CeTnpw-L1VR)#!xI91UhK8|;WCPmsW5O7)b2*+{=nUSYwg)GSNmI+m%~`dZjU2|kbri;fs>T$MohP{~jJH=7m-wZ4@#tHj({s;0?Zn7u@(50_>P6D>E(ad!kUBIm}7G35#Dee95qpABH*HpU6R!x`dd$c(riUf=$C1H<{6^<A1gVlHUlO3pkuy|)^k?-FZkxb`4(<rd)Tl5*Em$eqVe-3}GnKf12`#`$?p?jop0n3o#;Q%{y|?o{%_EcFqr2`Q%6)mrr;pR90{fQ`V<&er8;2Emo5m}xxb3g117#%~7d{t`6|0q1y-d-_e5R}&(8FOjgasBDRKA{#fRS$Ob7dGT$Cm+_WEZ3+A*ek3(sp9AJ>*Nr-?Oj`*csnPnaWi%yf<4l6zZ_DbWdpllaleN(5(&ZDb&WE=CH+D4H%ovMd=^DZ+wHS>q)$T>r+;pFkK1)gzjCljm1)U+G>K&KB8t{lAaTDpea2`jz!i^biz#n*d*<i@h;rA@@*ukij`aio568+7aaPJi~44BUKKC5T*pg59Y%QaJ60i~o#h9<HT(_@(d*3=UmkAFu{8Q+7gr&cluW@*J~9*fatGUR@52)mD~O9NX@c3pbEST&XJMPi#=HTJ!`9Ckwd#w)t*8@l@d;{knnD!Qq@^ip6bw2<*?bE02|Xhp;O)B7DI^z|uy_QzLqzJcctnGO27Vv=><8n{gaHYS2O!)UNDSkbxXyL@xmS$}l=(~*>#ZuH~i9Vq{s2c#(KsRkp69YQhBh7&f}-@Lfxf*#yFN}W>kT~SopURg@j)gr{ftHAy<64d@TvxQy`WjC+fAH-MJ?HD!4UEsJC3NVqu)t<e>P^tuhQ^Tv@NYq!ezl$?f$}jw9<~PWdMy710sSazQB+6<Qb(aVd+k>9;fiZA4vxEIXv5HqF3sU&UKitK7SSh7tb~QIKuS!&ScTbAtkiU@pKev`%TNr|tr5h-C%wR+5RJEM6z>J3+o<no|Z~J-6%ls6@B7C2G+>AZmdg{;6>0kQo^p{rvbtu)?+&HBd00HdH0*;IinHPQavBYQl0ssI200dcD'
