"""Active -141/400 first-pass, gain join and starting-majority proof."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial, isqrt, prod
THRESHOLD = F(-141,400)
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
    motif=G['edge_motif_gain'](); bonus=F(motif['gain_lower'],motif['scale']);pair=G['adjacent_floor']();cost=c1-gain-bonus-pair
    ck(U5<=F(73670242301,2**39) and gain>0 and bonus>0 and pair>=F(121,10**7) and bonus+pair<F(1,1000) and 0<cost<=F(123153827520,10**12),'active energy bound')
    ck(f4+f5+f6==62 and n4+n5+n6==218,'face accounting')
    return {'energy_upper_1e12': -(-cost.numerator*10**12//cost.denominator), 'first_pass_upper_1e12': -(-c1.numerator*10**12//c1.denominator), 'gain_lower_1e12':gain.numerator*10**12//gain.denominator, 'motif_gain_lower_1e12': bonus.numerator*10**12//bonus.denominator, 'pair_gain_lower_1e15': pair.numerator*10**15//pair.denominator, 'motif_oriented_lower_1e12':motif['oriented_lower']*10**12//motif['scale'], 'motif_cases':motif['cases'], 'U5_upper_1e12':-(-U5.numerator*10**12//U5.denominator)}
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
    ans=0;i=0
    while T+F(i,STEP)<5:
        l=T+F(i,STEP);h=min(F(5),l+F(1,STEP));q=pm(l,h);i+=1
        b=(T-r*l)*H*S/z
        ints=[];top=b.numerator*(S*STEP);den=b.denominator*(S*STEP);dist=H*bb*b.denominator
        for j,w in enumerate(weights):ints.append(w*cdf(up(top-dist*j,den))[1])
        ans+=q*up(sum(ints),S)
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
PACKED_PREFIX_SECOND = (
    'YgmsS=m3pmW8!(HP0l*PYy4LF*VIG055s_-MW-c6a_B#!dFDfNk3CSlr&9ZeZP%D%<#uKi*O_;{4`qw`*Qmkrcf)vSf{t+J2<hhG&$Ei2iRP1}u+H^U8{ty%Vn7^`;OuLwtavd<72U=Kz`C4cBj8q{lZdAwrP@X68Vo61wp-r?H}49OvTbJ1qG%cNf(cUMnJmOvcdc<(`ijX6!%#6o06?`j8g4(auP(T)eQouBSRf^*n+eXT<+BnO1m_@5KGDaXtYHk<1-3WQ(>3t(kYm1bywcR6s!4R~>ipOg-T`JdzG=}`9t#-{snX6?#{-lJ_-y=c%;CIO*XJo9Xko$HVQlc1v~wP&hym1>i`6U?Sy=nR&abTLK#Od4{DvWpQxViivwjs~@IsiOD-=;0eiqGmwAJ~|)P|6ryJ7RiA7}&_|7B|?wyAp@_OX2@WWA>(dxU4vEaoZQ_frK0%r*{TENVTc8NyrL=PoV3`bZ;c8bF&y&wB3pkDBikghkcf!LZ6er3$W2*)Ez&Oisf~2O7uNp<}JYfwnRmtdOS$)U*)GEQEvZ8krIT)fq8YJXJzJ+687MG-n$(-Nm(U{8ge*6Mj2eKoE32mnsl|n)hLRa*F#TbKAFTv!2Nkkm8MHLKxgQwW}2De2$}xvErOL!W$B%Ac_%Dk5t(L4;ZI*COMuJLx;bdIJlWW1SPY9Sb~v3)eTe#^D(eI&+$OR(1T+xSPlt)a7sUhCqzA5A)8ZHY>fqO(6D;OX7E6)0EOBv)e+xpdn-zNIvBOE=;!$~7&<uBgwu(E6tOz6plj-y*#<2@;V;n!Lcq^hzrAM%V0tP{YZrGEy7m*+`TNsTgN1U)bddapc}zR81h!+0Gx$#(qs$kwp*Io%X;j)yRV|Pb9E);JuS7!nf*yH7;us<20*bp;09VWRG!@NYjQw_ytU=fLB^S7J*OtzvslyLC9Ud(Xi(Yf-Ee@Md%*wYUT>t3jep06CXzm|sx1dW}<Rf+o!e#MG=#!^}zq?vN$@&Uk>u_L{Q=_+<DxOouM=U#L2(Qf|aDjm?G>s;2r7?5VbMd7n3B)SZRbLO_k1oU<d4Mr~!zMP*hc=i=R9=Qpi`x!yydI>~G{oE=B@;8&4Rg}d9=o1I$DO2Ys8`h=;*AIU>H3E-;{2Ze^+uS^22F=0kLZQXbZ*GMniN|{awY}e&5rb0zue}Wm^KC27B!gcGK94iD90L>FFqbDUjCERhlt+TZ55S@>BM)--A}ZU2<F&X#i$YGosu?v5DRrYZAcGY`TD#|9m}pJ?o>rscGI_S4&j3uK?~=sm-!MJl?8dIE)uaaEH>u#naKNobhajKR~0=U*NAFQqo6`t)VMva2yrvk_Z`>f_k&+F0SiN90ZZt%$<UC+CR*G7fy%X@G_k-QzyhZEWcu7p-hxoi4bYq+%KtC=Ihi<cCy7Mr0yxuSZTY9}C&I#>Wc4E|--q3}R+tG!lng&p7BS$6Bo+_h``AN3!6kqTVH3t`;G!JofTQfs@d=w1aJ4(^+(6ycgzC&I9g6hPe*YI}j(UBUWnxP7b!yZA+Wq3f%*GLKW3DnY?ajR2K_RZ{KAhC;z}ZDugT^`}5Rn>nrG;^KF1^|ZQ~Y(H?9S<1YXRv&Z`$oiyO~2zMA_sz_xMcxMJ(DPk5M+%MWqgeH}=yhcKNK7BczyzJ%t%nSPK7{Chs5e5MhQI5!~FP+K7+JncOz2eg8W%U>>$YYrWG17+k!$!mnUaJK{y=S!5O>^Msap=D3eo<jb>1EH~R;!#mDVSR(l=YBDoy(bUW+NvK2>1-OBY{F0bC0)jGBeEs^@1A1bdR{oeyz{&U}*}huBiYnL8gMjRIJ@fb@BqB&rZ;~z;vgv;|5!E{yo|JS$e}LkS{Mp?7xEG2r;~_}=38w8GGy319uO}7s{{JDSRK4N@<#KN_(NM)5y~BB&<Z-*#q$4en@7RE_C`8PS)nONU7y5Weae?;Gq*z%PJA!W(`zTd4Yr)AmZm#B=9Eimu!pg^p*aIW+NR*?gbiNegUJX(u%jseWD%M4eND35-xZ|Cbrc@yuE>cTV#bo+;u6r%)W#>~;HD(DQ`P*_x_YbIF4ujhyM|6=_hBDf1fqVd%PZ+y(sR%R&_;<ryIcvYXNbF{g(Sazw&o&F&GY~2llm2`Y&e`%0&?qA9z9-K(+n9`mNi0D<8!1kGg0cJxlm?w5s$346|8`0B<m|1-RS^m+v;%6yTFctN-h`?1GV)2xUnDpFs2l=iuVKE!V8zjzSWCvOy3jEx99r9~V0ISO*@T#W))9Ypvg!xM*@>LBGl<_lIRVh^XJqsa5cFb&m_^qeilzwtGtX4IGsY%5Ve2*IU3UR7A|Mzz)=C=~!vyn&URIdlJ|oo5UhhE>&AMYwx#8f~)bM@4cFZjtjv<{Csx<vob8o-X&}i^v7AmY|49}@V&Xt`0g->oLHrfR{?%@%TJ>*7^m@WI)h>=`w!Lp|U&voZKrRqB}wN`{(Ep6v-s3lW?;`)Nj3Xwo@RS)0QmClV<+8j5ex_&9QDsSb(wD_eJD_Q0SGU1@kBZol|^Gphio+>%)5;81-+dmhw1Mg3Eyl?H>rwU$k3WIEg7LE4|9x0Z}nEr93t5cDoBljex{^GFeooFx}l3$x2v+Cn-TiluP!tuVTwOF$uk%;tBLYhg^A<V5~$H4TltLWz(e~X~-!O0l5xkATu<FZs6s$}=VWay0Le!y2Cg&NPQr;JVMS_Bb?{gsH&Xusbe8U3ERr{rpE2wpG&c%3JljYv#=Brk$rVoojb@I4(7Y9c9|3+WY3+k^sEiR;Sk$?zU#i^<dx8dqxiOS$=cgzkulE>SXtE3MV=u;wB-lyV_b!?|jbP&V84T4rQca%Yn2pLN+2+|xO+pH96dszX=cJ1f<TFbf4;bbz7SQnT02O3e29q$Ewtq7O0ANQ<Z1X$v81cNl6hfn(4M?mvB<!DDgGkp#!g8*x=0bL;Z&aJ_U$8Mh89iGL|Z#RpZ>uX;K-+7tvFG&xdrjXfl?W#~KBnJ%s^K=EDT04EP0*0VtA7U}jNUlY&(9dn(zNf(r1RMi0&C0Bh$Fh%c2CWd(tTNBN8PP|2P8=ta7HsZb0e=iWS;XjEQ8e7Y}8vO?xlROuux2Y}%U|ir#U_lKCU16z_Kv%ETn*$29OPd+N%TE1qipR{{abge|LooJ)+?r;ry~9Xagicd6mC#2g73QZn>tBNq+~e3|ehqy4yq|wWr}L^lo`1c#A2Mj{3MzanxT6;XM+z%$r1ZWB9DBA~50bb5ekff<-#}X8(iq#pc*rbq@}aX?rKu^-_p8@&Ow7lUeq08E(>A9M2bR5ANV(icv0+{-c8r^5ZIJsD3ZH~5y>~*fRS3vz0|xByzbKm+yp7n@)vNM@`l(16)+a1;b4d^lfdEI##Sy~$y>{VT-gm188R0#Q0*a?iT5IvP8i>lyuo#99Z@w87P?8d_3`rqSt*Bvz0#Vx|hc)TNQgNIrG6g7ty$i>H3k4Kt4Z0<@8yn^PKuegY89y~hPDtZD{>WoK|Ed)FP}+jSS2ujn+tT8fT8C>dmM5wqGgpoY<y~UIidRclq}?kRN)aXD?8t%PemdPPcog2<bD7}MCWlpvTGvPM_(2;NZou|Lo-7wIE_M$eAfk$noHRcAoFupb^#F{%L^HlrWXx8f4nydegLyE2inmSJ1}sF3n%8$p0ZFtE{1>WVnR0g11;2FMG%r|vRe^m?#5H6CNwQWq3qyD25yw@K1$Zs{(b*@h>Z(b-Ir%_KrxLxevf%0^Oq3S>#3JoTA%tpUCMvM?1hD)BfO6W;VVXxfbN|fxtmolFy!h5<g(c>+qaIz>`2yFLFT&@FrFAJOLjeabc6HT^XC(^Nx;#Sx&+FNc7qk(4fO1%%ZP)DZw1`e;)xB!1LRP%zF6?IgGRvbL&xQj4b1`mbWJNt!d`OmL{ob9^2XwuWBPw@0$^*(KB*?gj?N;`FUGbCnf4Hgd?S2NzB>*YcN=Sjm+s^JlK4CB%nym!tbvk-#G7TztsoLDz`>2KPepS4og*E>atqlV3!Zytf{|Zpp(E~ix{Dw?XLe%&v6PzA&pR86N8JWQ3Z1t@ybjc{rOg!sDf!o&=AvVpsh>CxO<1}U~b9hpDwpq*aS~UGI$&y-VxMOOI{GIN`?u1nYd$xKb&T7M2f}hHWvxx6Pj1hLFbrYew60ilx0LdGz>5DKg5>b<!Soa;0JjJY<ZB~lR7sODr!9*4j0N4WCVt7e+1@QQgaIk*0<T(uJc3Zmuj)wqh?ww*9>Y2&OR7C-ab$dL((>K%JI=?`*qpAIr(2MEI8&^k@%v}(|CCCGs11hdmuKUZ88!j^SB(_1XYmq1~dN}Y<k*{prTQ!qJ8l54q^bFe(*$V=ynn+C_s6aJWA8(f#{3jh#>hbcKOgoh5wA{u^vRp`tqT3vM!$nVOjvzg9D%SH^xjycWPyFVO;KzxXrwQUXqaybPp=?)CsQ(S`P7Ks(e}S^0LYAWtiSvv0vvbW|qBt7faa!%r)zEpFBGgla#aG9s^+hW-F}wUGuzTK=mnP7Z_pfU3<_Jb)BRKT;xo{EheNguc-UV!$wO@m?OO(H(@-OvV;J8xby!Z9m#lE{)o2$y6B*GWb0Jl}~w41&%7qtaZEy`Lq>8Fk4dXzUo9sJg;o8l%NBH>XSx?k2+Xh0Z>Tb|!lzIdxT0Bv%b+~=!Df~kw1pRvm$!~}wNe=i;1AZ|5$?<&k-WonolF^g|2V>W)!b3mFaU!kcUf~8-Yf{XVRcxbY7bjz$%tA360w6}YzqTt%i&)cvCzWP}JuTF^-4YnAAP-|_vh=II;Vg{khn60B%&0CW4Q<8V=&V;wtx6|HnPHeuCOUSA#YZikrZyGaGF@Q_c@ii$nrtr#Xw0TU?8F{L=*gDLy>dyDRoV+n?QjX!ZKYO4RrVNKIzmxYx1ahf$8?>f<xJ?>~r=^<%0KjgBScYQArWBr&OT2#nByK$ugpwt_njD`qc#yZ&$U1V%MhN0Mg<H<O!SO&>)EkofNX|49f13A@wfc4<urR6tA7~6K>_xM=Lux<C-6UyVsVh3cV`I%jtgX>HCY9`bVo$O%<Ya!*{B6glcLh@!!xZdVadY&+=b?wk#=Zr`Co+lR^N-E(3yW)xlYd<W$+0M@rY*>=EX+EEOKE)JkY5nwx!{irIJz`)GX3-hKX_vO)+bZDhr!LQT3r;EcrtXmA5w>5ybKZo_|=x_mZguujb%K+IIEJ4IZ9-ON=ho;*FVpi7!k|g(;DHyFk&4tN(TU!b!Rtsq0&J2ZIoF4pOv@7K2K$r*hwqrLm2JNu|@(!pQ<yX*pFovQN)OnHr^)nR#M!}k}9yR-|K`YKY~|6EMNRj<Mb&<y7UkAOI~*1^od%fNOqx4Ee~Fb31O4OY1Ao`G+x264+pT(nj%!$f?^ySV&otQ*I%^kD>}FW8yLD}SrUE+(I6LFBb#jWXtYzlE~b*kC=Rf!)EYH9jlL&|P6o$t&aCxUk?hT@N;?lvi+u>Kw_PQ&V&78o-D)h8WpZ^WidtEd?u}(?YYgsS?soSd{-qFtOj4GtDUA1Q)-vyQEXg}UGvgU3)jaue=ov)b$7jaD?JUL1Wx}3YyEM0?eBI5Fo!Rhw0n4W2;Pur3cwZj71S!US!pfDGCHh2MZM3JSa%aO$dmh3YB!0Rr$LUo^&GHQrJQU)_MOX|Er)3TBs75t#H#(Hm@n6tm5eh&|H$53ze9RwaEpD%uZNeZ-ga|{kY?Wi5aD((@oj%t(b`6)=fJ2iE+=9}!9^VTn1V{%)0Tjv^KJ?G@2@SI86@#Eu7f$ko_>0PQH`pxcJiCsV@QRHyh`Pq84+w4}30=~2O&*ORlYD*CWJK1TfpuCmuuax=LwDln2ieKp@}od<%Vq;<-A(FTf(qW{e%Lbv3~&WVxX~i0^hgXfnH~W8;z~F{3=fBeP?_3p?!Fv>YM=*-A}URsePl3_S3*%lI8ge5<hUKpy6}SZI>oS_vFXX{-wTICr5UZ%tk|Ky_KK+qv*i}KYV7z<_?S-Z+c1_F_VBE<Yoc_L`%SB!Y-^>dj%^iiUr~V}^q62Cn(bO+aJ<twj^`a-fhh*M97)6fo{)E^i!9K=_vEU2d(>Rc%a_)$UT@<d+xp3Z%AE-<7+PEkOJVK$VaItCdz^Pd@Wy(VA0`&@IW)4xIha&!3ksr7LK&gS2zS}&GMvay$_1n6G&D*9+@e>$7?v1k+L{k?mgkx13V&~e?vsc{t8g(s|2T?}@Ocd9p`<-hl=Wqbf91Jy)iA)!M0mwC0Z^Vo-5eW=CZyy*sCX;DYwZ`Lw-P&;AHhh#M-Alc^@-lr;prMAC)dg#Mn*!aucftyau`U7Bygg0fohY_zgd4{H}euo>fivb3>~DoiF}!$3l<A@#bBf^+6q&Lojye&NWPYI4IMa`F9jd8lw^CgUbnkX5^>d=w{eRpqx5o^gU+I#sn(O5c?%_Vh4&UstyIh^wGu-hk$lFJ=f&?RoXi)a_3T$DX<FPV)_Ocx`wM2`;_M7kWT!_QhrzzF=8l82q4&E`uH#cCJ(p8IFo2pNLU5GOcmI`$PbvQ)VMN)pfjF_o->)hmjI4LvRUp$#ZPqbCNE-Jv>TR*qrCOHCugsaBk=H`eu&#k7xb}dipOQK!iy)_!tLZ+)@Zqf~_Xh{|H*ex;YXvmLB$;}*8~V+ae;X?W_wC~|4Q<arcwIlXHf$XoCsnqP%_QfMf%#R+_Xg~NC7ms^*4`)GkoxQ0@p-z*uf|W@JBqW?dwXbhe<Re&Ba-BiG}M9|>u27NN4h)A82zDD2q+tHlDJO<1rc(2K6L;sI9x-N1pcr@&xf?pO<CqK(}FI<DO#2LU=}n^^L{Tn+!Q#v7(Wv&K4Sb!>&yvraP+hpd6!hm;}I8!CI%ak3;3}VZYQjI4pFW*unOms#>&v>>SLR@W25NeW7U;93jh6asWa~stpxj7TYc`7$Nz5f-HQP;kH=j)%wJ3xM6CT?fM4WIhw>zc2(6lwYv0Rwl~LkB2Mt05ZT5{NntSuSC^`)X;Y6c_W2n{0c)SpG_SvX+#qJ_cCY?;Cm<gdDZ(-zRjWUeVxi<xxc*}?5V>j99<nt`@7WHTf()2*@U)|)uKJ&appoH;QFwt}=5-)97TvD;7HV>0Pi~3bZ8}-xC!JiWR0nYr^7FkRAXdr@6@x=b=fA2!C=3R%5Ecm!Bs#HH~(X;chyzX|_*o65#P?C$VzFL9fGQ@-o$XiWBl;TFHu<~j`zdbemal$6n;xc{J91R}F2{sm|C+aa`N)U@oPgDLBQF4Mb^<!FU9iQ9343jT<Ok8@rt!V%%(QW3>Y0e{^OyV{I7e^y?!vY|4z;2$9odozJ6L-0B3%q5mQkz<V3#Pz<v(ch_ETMImE30@C5$@<HvuIm;2@@KPF2E_RP?I2m>ayimKE}Nws$MXgc=C_8V*eVj`zlTt<Lk$7R3dfmg6>dz5**)DPc{!HcD{EIKeRfk(4F}Mkq$lZ%xHY==vIAHF4JMH^zLE-4zy_=9ZMGmv9U7r#=^-uy|$khYvHx<P&JWd`JwPS2!XyFw<DP((!eO#F<xhob^*Q8?q=E5OtdLaB<998okO*X#5OTJ)@3|+i!&IDuE}HR&171<H|QKyZPI1+J%y@-B@9J)M$G}KD@PJ=HNdCbE3r`>^qL?b!VYrirW6Lf@wV;H9v1<f%$xZj<6$~}92}0QU`ye)WgRo>YixIIKakW@!=+GmJ9jPuY9R+#Az)N5{?wa3&$6|dn60;%fK2^wQNm&$Kyl%c!&j1V7Zx^r{nR2UQ;<aLkz1!l8ANEJct=M(dTd^UREbv=K}HX(nS(m|Uq(k^1tv$u^aoF{#B`0{$O#Rx78^L00?v2gUEn)*^=xHPqROuiam;Qf)UTW2m9L9j3wQ};G+bWl>>S|le~&P%#`QN6NY)@zY0&3_(-sAm5lo?!^FKF_cAaEG`i^=nlJyQay2~!J;M1-ay#$#k?SM;h4893%Mvwm;NVyn-43p@cHXQ1Snb_I;GM%YDy3WWtmkWXPnqY2Z({wB4r7;2An0hdhBmUlE3CXMRY3Xlbk@5R2X3z!uvu{nK*lHnq)7G*#x7#?6P+tQE9TBT)W3Hf-WdIm|Y8#6!^r(<9><xDIo~pQ7%C=sAk(rF8w+7nM`eqYV{@D^!{eU5M8aA1#lj8lC_PkXsD5I2tH0@*n;lCNa_eM@S4S%)TNgRb^6%3^$RjIm00a55)XK6)+V2jjvr$Ru+o)jP=7_8M@@^GO`YVH3_jN3<Si#@kJL2$};^S8~q{FKCg@{b2Z^^_b?nH@z^E=O)`<x0CS#)`jBGvleSkGc?NPjs`4@|=$$s3^p~R5&7VCP({4B@6V9-DTHsOYlJH65(Lyu(c19%V+kt(>b`c0+7)RrcHz{d{F_>Ov#oFdmY+y@Ee8Xm4x5V5|~I4rB?8-M!ujjG%vK`^sh|_m4^WTkp6AxWDUu&T~+WZt4_7mIpQt*!;ZgInc$z1y=xu!+42%*fz~NfWV!SJZ_Swfc_?~@nw$aZ;2%^qF3bw%16_p5FwLaqE~%DhRA}!GoJivmSxMW~)1W?vC3=x};M@LpFPY1_+{eB>%B>A;7Z@*9B(JRPk-{H+f==yPnJo<i5~A+U$B+J4u+B0-I@Oz|y4v}q|JpF2OBpPEMD*P%YUHNvXJxDWqT3|5pPw7r{7MgTPq-%QqIvctdTZrQvrz^>3nB7M$3cp@a-TT++uVQmj^uhIY5Mar7=DKCKUnhQ07v#E2kn%XY#&%)T|dESgI&So8c(TqMz@8B|F$WfLgc_Y{(|7m6T<<ZUIo;Al6gL837XyFjPB!`5#M>J!>$5NCj$Lh`=PKz=~AC$$}$7cXDIztiwlaqK5qrNt0EO7iPv_0F)1wA$ZWy#)h>UZi&`yo=pLLWlQ4vL^;M)NkB>j~1Py58;TNiC$7cQqNAOV^A^WO_s*9Tq^S=pGGdF^B4*9^Vwo=z}ZHWujg`IUlw;Ma)X1b(Iz>yqSb6yk_BnrXta6`v1SIa*bb^ms~n{ex{P+X^qVW`PJ{T>T<IG^{c<bkS<c)$+ykQ<2axUTZ~;iF(zkM$LX4=Ot|(iyt9A-s&|Xe)Icc=sBN<bRrhyjYj*n7KT}Nh(=~cF@RnaDBp2xre@zOKaZ4AsvD5vLWNwohb1#^UAI<*g+4wkK$3A^(0BZq#fMdRh#C!PvbGfM|;5>a&n0<L3$HL*qe@gtq~a%5ezlSoQjlwR)zR;l8R15R9`^w%FixbY<lbjC3|li^z2J0ATfl4tCDml{b2@X>Z9O%%+iY57rs`^CsDALRu7Y}IKi^BCQZ<o7tO1Wn3;Sg=#A`$)ZHG0ifhwlw-sgyOaH=d;(s=G__bvW1^$nrf{%pL;?I?|9FssjOZmgOH2IKD&IP|30Vq=C^viD*p>_l`a`GG&Xhr3*yzFIke0FGO+I>WU^69fPhcD+q5FsJ9E1xEvKno!InYJd)W?LrFlY57IvuVEq1W*fGZQjvC&=e6>3^)#UM_l+mbxtTn0uW>m7rq)9_dJm7Fwj#4d}?_`h8B~|bGdfAae~E<VO&BDs*AV3J;x^TPG~eO1V=LAyeK-tf_37b#WiORm`8q^fMIx;J~2;!XH75V?iQHk6lftxovv*eP-zdSlV;^}`;~um!+@0Nc+Xl0huAjG3W~yER^KKyHPFy`;Y#KbD@A?qjS4~y71jf{<~L(jTpM0Sr}6Q2FlONVT5forQvD2|JI8SfN|(OFmpW3TXT9j8GrE3vxP^Y|JGGl!YIm-|C93oJDKgG?Fby&{Btd2HkK?(KFJ`qOliK_!+(;yg29({5l;*(;c;dsVw_7W6l?`~{0B`r4JUq00^W1$Uae&$kz}~=4<p2lx_UMZ8ZW!*XN*)~Lwb;LGrvP=fuEgtY8Cme3_;2<Q*SgQXhXUAyQgC5Pczl62`f!vQ-DrEf>w%2KHeiBc2W<WXW-@o^bybv%i*qQ=E7;s8B)IW23XcgPBN<UO6<UIXNWG2H&;i17y|8~(NBwC`X##qRe4SH7(T-=G596Hcp47LE30E!~LVdJ%e1)2TYH92jXW3<K?e_t$m?PZ01<p;C%a#h1xx_TrYDyHy?V@xP5b%pwO!B_fU=~SrtXKjotaSzi(rM`*xRB@8qv7NVSb^*GZOo9TUkXko`sn2yzhDHFus2{q4=cjgG}cfpfTlGJpp}1D2l@LIN?!1PI08QMbX4ihRLJ#K!VT7dH)W#v4f#{MG1124=KF}ldRY!yeciNIl1Mj|PWp{bsDBI2@d79WswVkJsRU(*>Gj&%2*{=0L1V{LRN8fYEh>!0AO_|`{mNi3^z{mPV;{_Xl9r!(d<e^V^#&0R&_5dnt{dlZ=o)qi4fH2swhaNxbqXajAIkLHM0*IkgJyRAKdZHmL&V+k<x%}r$BL+@(0nzRgFIy1@Z9CMF=z+^2|(QQrT9Ah&Kt92cBJ(wgi1#8y(MxMTJ7s3Giyk_g%qDj2B9uA`pjMcgdOr{x7=C*5Ak}xGSv_rdBLN)E3!YLo(Og?Q)kN8)Hoi<vm$lJsv(`w1Bu(C;M1osQRsM3Y`*+woRl6ya{DRfXNo7d5O)RFMpYErnHry@7%Zd*V@<khk_IfRK)eonxp5%Fjgz^#3$rQ?T^S2vY7WG0E>`@U3KC5O#LcklzwPizXLt*#7T;$n$q4*ffwDZAd9R6z&Vq>_qfn9Xz1s;ttko$mM~ef;N(aa9I2}BP81LE!O~SC&6q(YlWDgB?SaFFnlh>DH+OVgoS3C`C*VG{6-KNgEKYjYU8P5RZrb`>4uKKf!)*;~HaF#y4&`NqqEVzV}i_OzOJH$^ulVuz_uYCoFI@m)vOtF%#Xz&zXKK{2=F(?qXUk&D^>(L?Uugfb9V)~jz<C;3fW&VCjmjf+!rp&;E{D9g7uo~A_ul4U7BVr{p+A3xa`a=;WC_?QP=J68tjN3jSvREO*<u9E80nq~(-GH(<=bH<OFQ?J;6Z$@jvO@KToE%&aGPr(Im24o|!|v0gwYSjhy|%?_Z<-`<-&r+2{6%;dvujypYZ^17(3ua@qrmXw-s6~r3yp0};_?2)81epSPUnSW!oBkhR9C4uUW5-qb{!40a<+W@?pF#}J4*qNPeQ5n%moyO3Nb-qoO8Uo$v&8T;KjC6ZVa~ZCR>h9{(xi@kVROO{c7pWSy7j}X_Suqy~yPB7ph@kO6c8eyU?@{TT`?7Yplcti3vH_P$OhzJ>zQ2_T+*4$R=S_K^u*bT*zg>_521k-unZQ!fDbM1!}?V7gCJ<^iP!BvvbZY!q`_<FiRC9gd)<gneq-oX*pfmD=Uf>o(h)i+G?JuxRQxNK;&lMw8|i<R1<J>TwG2h8$~Ykzu3i7e*G@6{_(%mjtZ4G=nCV~#Z&^ns)VI4&;_@c%kvQ?IcJ(3QxG<7Nk|Kpd^I_V@rgc?{Iy)KsjlKn%fDNZ5Q+_OrNjtTi)Vji9(HjUO=D&2*!e11k3l*<bpC&g?YiSu2O_&b1jVIINv7@}f`<H#EPC;vM)AJRa1N<m1COS2o&8eOCdmvwb5H(6J*a3KN~jXkFYI<%eFz}u<-+gGI1CUQ=yyepQ*p*+qEdM=C{4pT>XusZ7OI~N^#+ki3C1w1KwUSAi}1QXnyuKbPxKj$m9K5P3hVaqecsPbhPMWDthZmHY60cGW>hj~OiN~SfFQ(U06Q{rQ+!W!8oSc}*Cv9&Fu;5?jKnmjK4>5PwdX5okrHlou0YZ2J6v!7jm4&VSM%VSx#oM-`Qm|b3njA4XSVZN-3nN^6-kV#_x;mCR<|OQpt9&^C$El;uUU>62GzmkHgW=J`rtDB9Wy7jCdjz>c&~Ll7KC#J{mo%99YUnUk!+M=vaSN6)Ej(Wx8;VFJaz7hN3<C6>LwSEOjaAMYVihYyA(dAv{8o8psW-qL}h$Di8o>{nd4}J2tI=!;ttmx&Hc%0n_Z@^_TI-hFxsxPZ^|_^^56B{yCZ^8(;(xna~#66tl5=#xL}}XWlff_)3Ic>8iTlHAAdP_yh%+_igjjV@cCx+7-VBmVGx`0kCq{4xJsSxOT5|ycqszZxW9p4H_5-G&ESTwXs0H+67)y}2a#^>T>#Rc`KRY}YfUWS(lV;QVHJCuo3KXKuh#s#dw=SNFh=vWt5E-<o@7EImj|W0f^#`UVlctV>r(1C_A(bi6|-Ty{^R(JW7<_<#4#JO@)I^9<XL>Eb7WI!u<4XT-2>F}jFdZ^ni77FC`jrj0su65rWbL~B{^GC`O2t&1hjKlvFWy~^L45`!|D7<%Oi@KDd5wg`BfiKeHgh0ymLVuQoor`@Ge0~FrczBL@iPwwPa@bYd!)n4798M(s^N9)H6c>?JqdM9c17NyyI*^AZj4RW`19h%5OASBRPmWQWBerGTAeFyzH(?INr*<_z$ER^}6@~0QT%CrA^u<#2caWK8>9Ifpu!^^SRL50(A<!L(}M%ROz1f;dAwIrYZZh^H}qDbgDN3#FJH=!xGO<M|Uzmweb9^V;o^;YvDHUi|lmwRSTTdN>e+$=zS#d9oHJs{(z)ApC(G8F!8R~IvM#)SU>Skx;oOvsoC8JOaqs5Q{}&!JZW>_j+!2cuCe6NW|gNokUw+o+5|l?7$AKXWPIz}k95qJfV8G{vIz$sG7b^>!m_d^wV*S++3n{+cVi;tnm32Qdf{cg)@M*t?yD4~-0`PC)eQT!-1<Oq<T-;g<G2a1iOZf_-7DJ*?5bW*`^dI0oczKXyu%r;84xYx%ng)EulHxJW(VX3Y``<)eOiu~l4&t*)dLJ!=Y}cb>OP(y0~bD`z5VCap~bG2i4Gd9%`dVK#gSE9>EvS(3Gnw;Hz6)3LC@JYOd&Cstf+of*!ZqMo{DJbR%!^M0Ac-^0j(WrILss)t?+Exd6=RnJg*>Nc<lN3Voymhm*H-Yo?R(#kWq{X9?Gh}>Dj$qhO;!-N+5Lb)p+-*Egh)8=aU5ULJF_rZ`g2XnE)XF_9xHIn%ip1#3G&`E{*nY?mq2&n00Q?$aECZmTvvwTrjl$#%##rpEG-da!ucT^W-H_Qb#nr;_MEwT0-L!NBIb~twA5)4GprN)DApgBChPPX65PvNM_LP?qw%^IzFhfe4SeCW$(N32L=KgpOI=l!4N5@_(HMyMNc~4RRC}utF}@e2|l%tvb&L`>OMoOB!yuUJx!S~{O(Jq_VYe4v1#JNZDb*5OLMVaxrBN42CdgUX-p!?B;c2~0MVFT_cJd($XrbY_Efz5EXRd8fO+goah=P5>L-Dob?E#c@wn;%p%?-yhw9Z|>fh6Ia=kgk)<QYo!?f@_gq&$3Ls|0cbxYV({MKX=NM#u7c?FH>Kaq@*gGEgv<O)d!^i<mxucaAvIBPwMEDbU#^})Lhy<uae(D|`_OA81JBU^_A5dqujgDFxc_>MZ7FX;=aC{pXn7Q8dB=bp>S)%?u^Z=NRVlSErZ&JvCLJimk}o9c*u731V7wJJtI?5O$B{4=$_tC_X-h@7ZIeO)9(cauLB9+30I>Q|QAZ4avC538$2-vFiT#;62rKi_S*hT_ZoW0wKhW)q&G4$iIE1IkM{CL)WR*!PO-S7PKdv<ZYH8vEiEHI-%e(K&#;{edf|4_%-+x()rUdfic?##V)e%nKqcg1fc8*x`=dJE|6@$vYnIg2w$wN+M%^(E}rh!$olUQ+swl$aTealMV;9Jb{x?IB3$(GFi`L9wdeumKtzo_;NZo;!u@(Xum)`&m<R!|F311*w8k}cx}}j3^uj68IZ8j>a8XfxOO3$&?h{gYj+kIzcTqCiF`U(o}Wh(h?6i#0SHY%+4iwto}mIa2VLH)`Dv=m9GOumY`gU;i(`p0>;FbA3GjR<u>J~{E+frZWP5F=<A0{D{4+&#%l}%D+rhc*UGYw$B1!tm?x5nZ5DDj3*DT<s3pQAbI$qk(UY#Z^_@(+-{>=Qvr$)>)|D#LxL%7B6gz$sK8`7ipM(lH*J)kaOm}f!ssEVUEd~+UoL5XBtY3>jg&Hng{<>0hGDT_!v&i@7dbA-uyxVy{x1SwHH*zRSsvv+H6Fofd?CH_sm%*tRTqfG<c0!jYkAwf<8%laU1`MLY{*C1Nlf?|F_xCF{ER>*$ETUOO9Wf^#;J$OT{hk}q$-TZ+V9RaCuTAcEU@fsvWo3RfJxx?mv4C2k}fk~CaQ@r4;k0%o@ya~_aZh-4#Lm)cl3ct`DqxnlPMz(uR&>mRPcz4NNiJ5&4iz~&0jK8=``{bSHZ#D~9LpJ5)V-ZmI+#eEC4EN4`KM9v<ml^m<s$gf1jTUjybal;$`^sQKI1+YuWeG+0?%y_a^lTBllNT9cVh^*o#l|vhzKsMeO}6uOhJg(SVQ85lC@6-A@4WSo-$nkes?bS&egao$PO3zf(Y6c%bXhQwb!-_skHMT!-#JZ!_d4M&xC6j3!sv<cJZn2nul5MaZm`HrZVo|~jAqjWF{Y~pd*Y=swJ!;ElqU|$pIhtIT#X=IMK^$ZWgRFG@Iy?x>3B%J(SLJ~dQfTbmE2-07WuczbP747UAuQZ;VWcYuwDUF2B9K=rRkiY|BMJ#mv@SXg7}akbSB8v6RW=7ju9f`l9Gj;@6?M~BCgAS@>I*h>kXMO$DC4pWK##HECi@(XY-UJ9w!u@^-asa)29gn!~Th?01U#c%SGx+ybS%=*|8*_@P(*i>`NCHc6aO&&IQFO5$xuy=DG0{Jq4b=2*0Y|RDCTUGC}^HXbBlh=fHj(4EK(OLX>_xz55tn+F?T#B0BywvVEhv5@#gu47%OX*3qgRok8t)D)ji7R9Q|$uNJUbk;xO@ktD<?v4JlS3)Pt%BWF%T*?efAKZ!1rnXL12{`>H;*@XVsuF^w!oZfOH^`ck-CsA_ez3Az_JBn3@81p8qNpxy2rve1)>%z9gW_>au7H(@W%4T)KH~<~x5Aht+(rR%+40nke6RS$@wE$yp_<g{DAHw+1Bb7fA^<-M3s*0o?Db{!`hAw$<Yy$wquD5Wo3Fzwh%XVAF1|0kGvp7P&H~n*`;J^#UvDU{hqNCPuMDKaPK+cdrHNnK!$Vz%U0oNN^bG@F8Z>)A6wIT*8$P&J3xwiWb!P3=^6L4}uyB~PO>>H}6Lio9ia_AHN?3{_#)%T>n?+XxV@`p{E2R8oD3QAUERJcFndY9_|^m@6Z-8OF7Ia>YF_BIcBaZ)?S{ZGg^82BIvW$rF&jeeK<nWIx!YfD-v*5~6S&WLp|#XqcYBT>zu=CpAu(Iu7*>Q{C^QSOGwpY@qi=HcFbq`$TgdP6fy7iRZI6RP?WhC{}ZCeLeL9f4N&TgTNjUc$3jHXLc5FdZE(`31C?39Vik{V_LxAELOHR?XfOD7rb%VTyghIl!U?$Z!0mec*MsTbE1ollF1BrlwIC`yuE!4q<Z3LAplbC~(`a%Eux3Hj>zB=CsG7X-e)Ueb@%P=;IAtObuw8INx*QEAuk8?H{3Q`h5L*xQClbmx}%?^JU{)3E1=Z%%G^d_s_r|KmL6*#1)JAxK7a<Y!m+mmYpPyldvLIIC4#y8h81ML8%rIfg5USdscYN+uOv7i8r?Qw>^NU_ob)eaRlb1J)5HS><@!aK%Rm=rL8jDW5J=m63t{=y@LDcOj8D#5{v1>(W+{cdWn*t>>gt?Uv}EW7rEV?(-d=_@9^KE?z>Of#cxMtrF;&Ma^}xB7AdxM<{P|N9T>|9<%6&P-ZRR&DLwLW|L_B#o)i=D-W1CXuh55g&DDh^a$i-;z}tkRgo|b?aU=j^oUC<FhZI1bax#-cA(r32V6ZTwi0(|vrF{!s<wXjKLw@6`w_vb1E(2vahftKKR0?t|7?~D93OHue$*`_nzO?&OrkRG5IH>2wn1}5}N(Sp>b+Dacb-Iw93W6CV4+I2mZ~2qpD^*!u*g$CN*ZT%$Bk-Lq{BALQlfQdYV}Q&qy)4`~t@@*KB2ENQ3<#7c5b+s1x{{^MV@f^eXzA#6<|0rI{J^gYy9i#&GXYex(YXuLHs0WMdNlAPye2tk=d0q7jbFsQz+#Hn3Ji5&+J*qVll+%uCn}ptQwN9YbN;!lMiCdL;-zWUXAI6IF#V|kTL$Cx_t50q2yWS&<AXC}cNkNYZAzfYL?^=#Mot#2bFiqK65<Ts4bZY_=zQ_n=V<s$DYtwhKDk-$*#L&PO~AQndT72)?U7piKpL)bFC*8&<;P8PwO>~ww-=yLwfqzZ1@B9<{DnH{d;DLCnzcY6%<Ef^*-I1QfZV^AN4J9?6SIE{d6fn9#G2H?Oi}Yula_$fkNn=UD^^N9YEgYpy>Vz0qb;XI9Q-2s2FFD{P7B7sec<_QUjry+uK=5e{C^0QKmF5ozmb-{il&r@u#(s5KHCVl1OvxrRwuj{P<boRO%lP$a*+jp&_mb|CZB*(MZ%$z!WiO875T|htf>?xWLYhZEsg|bLQ5e)e#tqHS<Ur1>k&qWm`>=8Yo?pS`F9Uan|Vlj{%{gQ>mIYiLi}~Pq6uv$5LqEJ0@81ph#qe>=}j#y>3I~6eUGfr3D3?Z%!p*?+Z#wUK|B8~N7SevqN7&i1;B5#Y<QpHUnZPDPVuJ~YHP`_Sr8z|xu_!RI+h`9KjCo*LZ5dk;075-=G;2ZT=mv(MZsqIhVZLZ*2|G#9)&T_gb^{T=?v4Y81sjtr3K5Xsa>h%%?Jobh=2MXe{dEV>o7`jWi*rjMx}%qZyeP1Su^u_I}4IaX~jh}nIP<~{t^(*lsZ5i=+F#~SS?Eny1#M;?mo_#gy0Iur_Iw|&Dd!9{W)7uH$yk{CRW{9k?4%N6^1yBA!8<Qx`xEwy_&r7@pXnZ$|^nw#x@8|c!WBM<t$yRX)eypPzecalZ~wPxu;O$Qla5$a{@TBttXXfJals$wY#qDt}N5wiTz=75Ja<Ec7!nw2ZKO9Xabj^mPR?DZw=4~64Sx|{e`{siuSg&!LfA9wEwp-03o?^B7Gw=&y`i<1=exE64<yyFi8~+`-Sp-A5fi`HPj%IJI+L&CinB|;vP242%i|ADf*J6&kCN*G>w9}<u49g7*@v(FJ8I$01(-@0nco@HSA71YEc7~Xy$c1R@6rQ`?_#Yg}M?`Xe}p&0kZEmSs;tDn00Qt+nvoZu14AQdPygvFsC`0*qk(9$>zw%YRr<|kzbpGBL%!9u-uYG-^~n;KLd`i7?QCzVf@5R<rThkE6-lC2?>q`04Q&uPYbm|Y*FsnS4+|nBFVv-o4=$$p-qY7#mOC575cd|X1zfC!yBm;BdFLjq6%{BCItu9{O|OHAJUEenxWx|KX21Rd-gTzQ4-IOt5R);muYz&_Nw9%LyHJjN41D`D6??!UMlRrQR5B14Ma%i-W$;Hc#ZXhP6EE_V|Jj8ct2pjFcu*htcg!-RK%n~G5lvB>fZNe3p)t3I_n4(2&W@5S_=2*0x{5os{-b<rZWSmBzILA9g}8qW|)I6_fhQ~12=qKU-P`Ffy3Vk;mA(g00N3~SO^mTJZT9CpsW9Hl^*omP6#d{Ccc%p8y?%eBTHz$ZZ#I7e-&ydYE>N^pdw3!)J0HI?V*it_vco$R%3gV<~cf?!Am3byb9OOscMNkI9BUe>CwP@t>DKQ96n6=ZO`mjRy`}I5ua=Jdc7_xEP6uOqTi?*<NU6;=oFdUaqcdW`X*|VgtK){(U$ZL+(?f%GC&-s(=VrO#T$Y+kzclOId}STk_Lj$qn3V7_=j#l1y#e(_U7b-%zv&<Ln71;_gt&~kQI-oHFk?@ErE0Y+%$xT+)IVv0kNSSsSda<xjs3|sCI`Bj#dET{+D5Qj*Dksh~ak?MZR2ZTDwQfdbg?ZdlKRu@OV)EPE3&CzcLd8Rp3Rdt!512EE#){AA<Ng=DF%%b=;RstIS>V{b8!GZ;Wl9X_!LD#Iy13O-D@Eo35}>s;fJGbuaL?+Kk&!&hx27F7;{4k+2%MI+v&wtUk6#DE`_yP{pOo&X9deGFVKKhkMa!wgJAqGg7#uRIc#-zKk#2*fU_Cw_727)-6xF<kuH|#o@RYDSqqj8Qo%yO!A5IK?0{4#<p5Hr4+sf`j`Z&*1q<R^b(JykZJffIl-Gi>N-$NcE79-xePKHjT;y=oALZX+IF!(zIwiOcOjpJWvIQ&#@3|9$6f%97};zTL4)qy|8<@me`+VqnThS4Zvus006K4x7PVK}&Z^Z8VxU*c6@x1GQCAZ{qG4e3g_Dmf9gI+;G6r~3e5<gs*_@y$3S?qsD>(||p7JJ20Vd*Z{b{O2c0-ztUIc<h-}(#8jx{A3Z6_<u@XUNg(5$2gU~Qfke@W>+<WIyED}ixV4V63-q9~2x$;FPJV)$nKQ$X0>o%mEC`hGObvekzd<-iC1m2B?Cxk8F;Rds${mYou$0fiiDhE`_Wo^DZmr!BRND&G6wxRl*jJ(k%H#19{?s@6b*{GYsUQ^?<k8q??YIct&8dQSWs=<s{>N6pA@O7u|!xj~T2h1oXPEP81U;*WR8Q-OEG?Y?&`U82AXZx?}HKlf58;Vu`EZn-;Dv(C+pf|K?jhu@hju{k5Gy2Z~%CQi9KN>C`a>xt!>!+|WPnOi$%0g~wtuUG|w@J$cpRtmJt(w;+*YcE>k<*gHuIWuZUqNrt$SR9T}=M2%pZRhb%gZ3u(#R`Qz&YuSpztlFnTzxu4$gs^L;HyiwAJ~feIp=(w@;6!B2g`yo8k<-G9V`r-o=k<D_a7he2y6CP5U{pl4e@u}EmqPd^Rr}qYXfb7IQjSb%y0ng7LJcpExR0#5s~BhH|?3`eWwAUAa=QHdNF{?x-Zh&qU=O{^{7p6CZRB9Gl}qGpj>Zjsv9FSLYCLhzOWDzllyvUdwA+*Oe|*c87X%|t-nLnRFG2eV<BAD-af!B#>z=jxWjSLZ;`f+ELE2NgFQI*IOdKkU@^$fKL>dyywjHnP(86sT@|G&v<~b`8oIo<*MsZ@{hXPu+MpY&5)$gf(~w`#f?R<P?M1e%Q@=rz=Pjfr6(iF5$6^tI`2%P9J803|7)I+E>h3`eN6d(WRUyIrdOe#%HYs4>>&<OSfFy#2E&vB_@@aPA6!!Ond~!y5<3~IZDrKi%!p;UBg^jtGHQFGsD?wv0hsh61Cv|I)f<L>So-9j9v?j>uk&1TJePiqQJgml=EPr%HPFCz6=8Ymjvj4sG+l7sW0``<)TY8fiL=AgoeR#rkis4`6d23(}D4r(%qY&=YA*TK9Nd&LGP+_KkVWV?h{3*H3acKB;iE=D_BniN?oQKs_1(*8WVyMS7P6{??N}?D9GlF}nVNFKH^myyhKx6PpH`8OXm1WAKRo#tzf7HDXroPyzn+ZPi#N!vXzQT}G@w_K4ZPO3#AlLXS#BDBX!?+b$x>qZCi<f)YC`#u+5+5YNacmvUYwO$$`AB+_L>o?!mz$RBDoPjOKQNaZooi|v8YF^kyV2W8Ng$F?KLF-Jpv)4tVRool`&CKJNFjEr?IUs7LS)IW9*4H4AH8ERtLli+nS{7TI{@A)QEBI(Nf^MwZjCJyFx(icezdp^om>LRjGN3%X~7c~64DSskzoDZzxPyJANYP4j>*vj3sY^K8ka5->@bOsx^CJUs@<R2U}A-AW8OyxYzQ1M!31FY-CbyiSIO{^Em@@}2RJ+Ooq%2*GbCt7^A)4Qt$!(ydJV^7h!`R3Kzx;W-$kZss=FoA&8vnEKa}<M-Sa&Aw<EJm_W(xi96h?6SE{fy2(<;VF~-VjV$1Ht)&)L<dr%omAd-ZxZo0`nv}kt9H+?yM(^jQa1{fetk*VeaAq!-WXkyD!R)w^w00hGAJ1f#{%)iqK#)`EDY@9;kR4?f&mN_f5RQA{!K{S0?z_Kkay<8sk%B^s%a4;4VN?f@~kiiK&Ou)~?aF6>O(aQ;`5maG7qxiwrxE~p0<$&@g;&Z0)3CH|^8^HV0g5Tc9{Ds157&}ZXBtszdN&sNd{ZIQb;CarWP5&a-i@As8B^JjZTgUA>W>+H~GplU2{?#dZB9e~e+9TOor=#}f{zDq}d`O)Q)#msZoV1?_3Jm!3UtO|ma&0ob6vJiBI?AK4;QWI(wsdp-qc@J~fiS-_ZIA~wazObFtH^x{sC(sM#?hA;`V-6FvN-%-oiIKFMR4I^{Lu|3YuYyREqvP5X_%S;2j|Md5nGTLUM5a0w%;pVGu*-QwTHpW>OEg8z^tJo_e^_M8XwRSQL{IMr#Nd8L-O)MqgK~zWVC*CO>z7^xZ9(MNkf^`E!)Zd*?+uM33rHSpsgmA!1O#RYW5sG`P%XfM%#E}pR5CBv#EC<ABDElm{}hyfw%p^LrN_yh4}L*|I)^W$gYG%CSS@l^2ovJQJg_fU}w-#HNzszgW|xUo+#4I7~7)4ws-&qIZ|hBvhH3&iZnK%MEMdKuVzfInDTlUsjl8ZwH0Bl4J=kbg2^(!H+Y;q5^$a}r$jnL32*75U;`BL(62a}gvG+;hr25pz--nIdF<fsk6EYvz}O0<jufsiY8Dk=7?M-}T;=CWYXd*3K|Ayc;!2~>EGA5R0zQA;Go*>ev1$-nh4jO5;|TsyPyZm+LPE`^b4J1nrJ{me0TZ4BzYvubIlHQy3i}DaTCH+#zBca#*j*r`dD71&;mciYe%qqcq8YYfn~kCRD;(d9H6d7$wEHg^RDz_+F9Z0?ZxRJHp+$RFtWWbE-n^As4arqNz)gYmdJ|L3=}H&OH)-NW6w;f~Ka)oVde5CpSg;KcajzRZRs#)Xo3t^1-9Fhy03o!ehGT4t^41{m_xY#}UnqO&5+&-jtS-qcvz*Elc&vqsU5rC!U;k11uw>^V%@GX*+#3CbXC}V^EEr(UJy{<CrYPS0|F8G$zg$A`#l*{VxWATcC@qPzZZWfHQH}t%-Z}vq8R15xnUC#A<qR#2UY#Av$$BnAJXj)Mu<glz15>j<8CH6F%4AN#K!FywNL3#$(ofU1EKJ2wU9(t8o1<Rg<ABb#E*TZC3~D?a`sJxpMusNTHH$7>VnzQEK?+Y;x3G!1(mhlhjq$zU3{QhHN%&2SK2(PqR2>`4@*mb!KLoV^-hftlros3;boY_DIRjcmgB3YP%=(?aTDxIs*0JZjrMUnqvZ~}bEmZ(6&F(?$w`MGyH(b2Ut4o$f<f{ZfLmfp=6K5ugblJQoUp3sDAf7dv?s4#j!N1I8)6~3w+@0FF3iGA3d<GYZMJMWMzTI?^Sl_?1J2kA%5A{ldm@cHY@%{X8_sD)sJnlxE#Y)VZ=Wlb^o^r%viF)3x+1Yf*-$6a|TfWNDUT6xhDL)suKrdOtChpZO0UcwJ1+eM-zynkgINJQWyzFhloNGjwIL5LTrO>Yyg(#atBrBLN4hPaZg=DR!kb&0dpTjje@cRcoC=TO(pyfz~H~o}bGk=k8k%MSU#+0#ZmIlNySWsb2lBKfz=1)Z6SGL>T;JL&!ruCPs)4nY^tecJJ6c)Bi`2=~SU_RHZu#+df<OX2ZMnUJ;|HY;5Vn+R!0VXie;RfRDsi5(3^P#u}@lo}P{&LvxDVA<HM()8pD$QXK6LTv+^Hi7aR~w2aD?pfNwpZoyIbyrzHvB>^X05UY^UnWV6tkm?-SNf?WY^c94Hza9EhJ`3rN=r6*7ka8W}OW>g=bkMvqBe*S|zByop{^v50XU+2bPNOQ0lDs-TnT8sOQAFfE&M||DsV99bcTf!4{V3WIOt*@>sj(1m444Gqqtbj%4RK(xkZwjW8&H1D>iFMugrFi6R>Eu0oZKsoGn`O0_qAaoR2qp3qZ^>(IoqK)~cyM3jqBd6M))F}1uri7T6-a0bScl{6-2>x0JZvN#~DVwjcmb#Ez*iVPxaox<;F;PYOu4Rd5gO1v~x#c>eW`>$mfS$aF+AcdvKpN*_qUv9H;5!DVkYvFRVTmtjyHp5kYQ0AIZt8S0se_^!!#PDyVX+c^7JUK+jv!b5ce@<mLI#suPf^yivjP<EfhYzZS%V~4OHr>)a551gq8DKAcO8@cbkVy3w=no`H_4bJ1=|D=%Q)4l^v1sv1zo8@U&jV!X1IHnr&G1Ne3Nl(#?6*$w(g-f?d$|D*(MlW=vGLmiEEeQ3+0s)x#_?2|<!`{n1VbV>47Si@Wxomd)UL_l;iZh{Tz(+tAh8n-#(~MF?QD2|#yw>=Sg1P^)mj3xnuHdAI8ZLH0@fLxu^iJHMe6AsH8THB2~q*=<M-1Ua>j=?Ku;DHNDL8wiwAn0yI`N?xqjsQvp($RQUyPvWTXL+R=cV!^DtNcqokVv0yr1&V!Kc7trk`3W~M#oo#H^;Y<K<#5>w9O_UpcjNb4w)q^T4l+_s^dQZZGv;~xQ_t^5_|raO}HcFPIRl1(FTH_xKPj;1Tqmn7C^%vK(RmU7W8>bBT!mCbF=0QlG|-wpFk9!zPZuYO%iEv~xkkb?3tqa;z_EoZVOUB-`x#d-Kn2hiwcGT!&<qphRxhDn|0*rEdW$buN3Zx;8Gnbu6?lq~bkM$)C~D^fF#53*-=beK<HLKc<uQBJ~{LE21qXKtc5D5eP>mY~jCPs0%Ikj(v1&w>(u1;IWRK^CAeg<%KjShe?gz5eVwtInd=kb{lwm%#oUxp1Ql&nT>lP42@Fd3V;W>F8ol%dl2y#RT`wo#V4W7L2<p*)dU7!TX-Gg=+)s56*b)`I7CZ4@vveq}uziKAIz06tJ@Beff+|a^`ItNa%UpFYEvSWVt>sGT7Vh8Bml~dA5#U>Afo_(IYU?e&?z-E8)V3FU%j>QXDyXo1#~gO}6mhx}|y0+jm(hQPdiqqpa-3Qn$>%Ico0BztDM@Z)jKTkXx^i{D$!+W1-$j&lzBbt1EDFy1oAhm<C>5bp4H!uof!-?Dm7Z*ev6`R&6TBx@^TL&3osXKbg(FWY!j%JXyrg*qd}QcrZ|lzTNX8nOo6WwuW$sB$08$N7bg;o!4?&Vmr%W)RiKJlhprtq2|8)w^wRF0Z+^PL>nJxM+rv9;g`t$2jLNTN={h&M}WGpW^0wWnZHiGhpjWdyV#0gRj7aDZt)uZ3r&*l*7WZ5Lf|CzzYPn2F2Q^yYndZ4hbkG2h8*Y(!IDx_$8S8zClguQ(5%`OMUEglrKS{Ry`!oyB|tU_jJko=Y{_a=BK%=&;Zxr$0I`^7nl~clFiO1iKj;SBw|AgiG=NAtC2Zj4_@56c?C5bf^j*t-#{UvsT4E((0O>2(8MZR}m#G{0DB#0EHQ&PRW&Wyf9`tW>nQ@b6$hgG(PW9x<N)AvU_2Mc0=z$;4xCUm8XhhJ-j2fHR8IcEd%Gk3I3qF);ZS{vd)$Z>ax(O@czL5U!;NWyp!I$iIoJS{y*CmZi?&TxT)Jj$YQi)qfDG?DM9@%NjXBGY9q`3aK3wP$yWE@r-SX9vB^R}}?dBGpn@Rc=Q4kDB-5uxwd$QnC>QBhaG94k0a1M2C88q2)`h2|7z=rV#_daFkv;Fo9RUY5TsBEP;*RmOLO%GPuv*|Qj~QA(nNen;n^C$VdXy=^p>9ApCS5;}?W666&?b>Vnrg)nJx=-lz2E?XQ5=z@NkF%eDP)bnLlUkPNvx4WMG+!5<Qs!3^{g$hOp@dA?`0dPBis)05D_*OOBsAhJs`Z)X&O}rqC!c85iAf=&B<v^;VMrY_&f&Z9Qoc9#c5YN3)sl~T<1RL>EUW&h|_Yb63BRiI?e#}jG6hLq2-rQ8rin-lC31uw`)#V9E&@Gf=v1RygVdq=1qrEhUjSd|>RhX;{7>=VKpoPeoWOkl1GoP$r%LIhSq!#csm3f@oL-bIqwyp>VWj`HkYX3aZV*~=;0~tud+(-bO|6KD+jaVn5shDh=yv|QHnmbu3R!Lz)MSxxJ8|9WlPawl-Z%_?DU6xUvhHFn834q1Q6leZSI{v@(<DmsJ(V1WgEK{$-XWMA5LgWl@lh}xEn#StTrG-9}6mz74w=%HU)qcEo{{~5OeW<jJY}uq1X=Nq-*1ZdJ6t`r^_3#<Q;3s5AN4a&y$3!KH2<&&ex!ypF#WLQ_`ROEp@5C{)M1>ELd3x<UQShLj(@K&~=2Nr*IX&k7!SM#37xltzYR<tw1^%V~(-N>WV~-o)t_=7Un){w8R2zhr)Y?qBQy=`YuJLs939%0+^9Z>W#1Ghf7Hdpk@mU`{_kmOlSW*DqIWYW6MzO@JxEF1+t-%7_s-!v18~!og9e*=eR5@Zu*g^fPydoPh-tVeTcGlv)4QR-+d_+=+*76!{L!R2FDQMJ~eQZR$EH?f|hcdp7?QYrEMPJrYJF8!C+RFsyluFPeVtDEyQ{jUGYf{bKb5!LVi;D!dw?RctNMqdKaMcX!u@NE6;b8!*Gw#Sa$`|u^t4#+D=9pMl%RZd_E`pnVGU9m8=Rz>}XR1_)YE7N$DoTe3^^7D1O5OK2+jxaxtF;erc-I``gZD64zzF^RXW#-L*!trbwKQ?AMde^8Kglt9q7gri?crZeML|)AcTr{~xZ4_yoKtAgKYm&sb7pYL=z=$$SNXcT0Cro}?X2G-u{HH#nh(mAr8~b14yZZ3?{Iv#$h8c6>wXq`zx*=PsAHwpaZPnRJJB_uz#JM1V^RIh%S_WeS_XN>o<HVV`}il)g$WAmGzao;s4l$7Mp?~fT}r@VEMNX<PcceJkg|Gw6%Z>wj9;3|<&E=Y#Pg~>q9`EumniFBnJ8|FLUmJl9wZ?(-=dPA3vMrgy=MG#EMvpU<=LcDRK}d+J&-d`vw(;HCmBW;xU&PWQ?b0^6|ql>-Bg?Y1+Ng}k;sCp76+nT>@(z2kDd?6J({bGZK5d?+D4%}GCq<qRo4K2)n$Yu2N(t5t^E4RTbXaF4w2zEIKPnQj!Uy`?FO{eivnl0|H19s+uy)3wrq(<+MoEZ3b1Bw6e2eWz8(*@e`bB^{yTdju+aci47MU8N_`EGzE_iCJipobc6EQ|K;7LzwlZO#WgxfgBa0>itxN93l2>oQrqYL=&S9*@K`F&>gTS}>bi=0WDpz^Lon41p)Ll=rQfy*@A3dLqEi&GrZNRV!aD^U=rI-W>61wTz)v2KN%{*5$louK(Gl-km^EjayU8^3AzjLuY{c|8-kfvWI*lZ{TuboI@0Eo5Mfr%Tr=3cxpw`k5#FU0N&6B5|v6DXu*>}R#*j!D^9L1p9hWfrRtnO`0UNoWBM=emuTOB({kuoY!PKm;okN?E8nLO8KSL^`(c698*(tnm_0+r_)<BPU$!E?LUe({Q@>&eYFpo~cVNgT@|31IVI3f7qAVdD_+K9ZnoiEakPNsTKx9K%v|Dggyzjf_--zZ<L038$A3CpXe~gY+uVw)E6bvcRj7AL-VmTJ`A@>QFWB<_36KlY55k*!w~v9^VkLXB`tU43vZQVIOtoOD`eyyXL>HWu6jOnp{1cSUP>m^SHr==WR5xYCm&T_)fAMxeBiMab$tVrJ5*8Zq>qlEIDiYSjAzLwD}#mtN5gahNXv_yD!hPDiK1hwSc#Gj-pE;|2Sj4=3!3;e4^##*W0MN9%va&c3j8e>3vZ=GhA%q-r-96vu0j4pJjA!H9S3nasC4UO={5vMfL9)xbb~5qQc1rjF@4}MRQ)jj8e*%3a%(kx85$!3oatq!6&^&4*x@A?Rpc)So(=sDjW6PV>I?QpPXPET>o>Nm34apN_b`3jOFbNA9aQW0LArR4x<)A^X$E<aUUx}E4l03wbbco;tk()Bwb?>O(fyDVnuA2!CKGw-(hyzbLX{#q`b`UQnCi_}Q)PSjyU`I3+*Qi@EL|+P%U6jF3KfmhC7P|=5sh)MHWFzLrh=#`o5sj<BDaESxDIURa}6U(7W_^-IMs3cG7D4<6KFX14HJLbN_&r^Y=-tDPD{)$z~1<Sm`#MTOX=?c<`t3)$JpYsC`t=HOYGz#NacP?E0LqH9zA@POr4$Gu~{VbGvB462h>$)$5h?$aj@b~0g+Xvb^h*;Zc-(<aTBg(Xcyhl!?Mn(FHo&xe%7^e*If$QwbaQd-6&$^9Gv_@e4Gh9hOtN=aMXB!P$@0bE#lnIi3LnMCg0CqQ#OQ$DkK@EClBvs<7r0B6|l$6VMpvpd5KE!TD%|?Wi3^NxE0+zo*TBVOcA}xe`OrA+)l8@*Q(?tIn;o`9-*>eMba2zp$QORGixU<7qr;c2l1aX<4{imHmht8lHn17=KcQyUj*X$vqc64Q_b3ry;T=TRWXQ5%NcAyzhQ~g(CsR@8%+1a5qdUUS-um2F2&ltSnZMEP^_aHs*h8>pi{u&zRb-7M(?sNYN`;}0Pv@}XTvdBy;r|yDe=1poUYJrjdw#*?=22A4Dbm|WTXRGI6Il6!Fn0h$!r-1y{q&NMM52t+S{5pI*DcsC$#T2*ro(iZHa8P_)F)q3Ct!kSEPOSAa^v9r2Ykn!)5UVcBnDk?W(<?G3|iw$SvMD%IW4MORq|wND<{4(Xl>gHn7KvIy|WnLw*|G=W|{%*0XCMUXi5-ygBfx#-vq-hosFFyo}Z?um(E^5{+Bvn63|RY3La~v`iT^-Y$DixS-_3+Wc9G_&BGKRHPet5a0Kxsh}nuN7?)inq9+SR7#N|dYPA%uo&MyR_zBc(j>+r7g7z<IJ;;DJiT!-QVS7$7hISQKHYH=xAcp%sPjS*=r@&51R90h^>T(^qfU3`Vf!Y0NjV?T5FiNl=x)M72<Px`jiT=epiuAtre+u01Rio*sMf@C3!MFwzf;#+a_{*)Yr|tY?>e*yNj695x-^D|@KKzAH_xBPxmwuLjOw*fDF5S*6hjdu!L4TkC7s*T?ZOzrCAM34L4LExgC(P^BHPU?yV42A7v&Ep2-Fwfvhi1dqhU!=v~yafL8Und*Eg3vLSC@Bau5kf{T6^?W<iZSJw=AOA<85%h$B0_;cOK09o^#h9@X}tusUh6EP8a-V5NTfRg3(%GkMv@9;D|vSeP5;7SxJ5RTSmt*9Qa46oBywXYihUdALzwsDnFrtE_lcf#9YuIg;wCvr!N{+;yRiUC3()iIhzsWIsCaj0}r=7VtJWM%#>*ZxRY&%MZU}duERfbVwK=Z-yJk!~#`xXnlj2v4GgnXmZ9f3)e0KDciUPkx1|Dy-Ugk$NOw>47}ERm;)M;<+OMaz!>bFfL!>1JcxU(24pnr!p8Ir^L7q@MQ<wMhcQp%EwzYtEOG>BaY&g(o;Kn6i8wxd_f$428d<?(&ZxA1h@5Jd@V|P`;1X0+ZPs#-RkWgJJA{=9JYti<re}D*tAz5>w7!D}4MEFAfJjw7>4SFdad%1&Q@uWD<9qFJkE!dstYQz@;P!TY&)*5BksVj0W~V?zoJDoF;Zh-S7S~l+uebUN?}JAZ3o;VV7m>)90y4bSV#WsLUz51B^XiwB27|kTfv8v(@yE$-IK<#sqP(-i*E<U_bGqC%2~L=(7#Bf0pTix<ivJwAp5z0Ft92{HbL>C?+v@odJ)Q6HCzkz_Qt_^^&K`2wl@Oc##&O!5-;PMk0!17IxJuK+-W^8pDugWyGHZ0z6!Omp+HnV|;D1L0QU%z1GY;q$=Qcl?l)U}eF-#c{hObwzcM1i(K}!*{SU_i{YlFwtlZ9Oli>)R#WR+cLkS2=4q%*ZUg$yK=<z8k-Mwk0oQ#PnQj;bNUK%}jj&L>W(bp%xX>i*8VKj7C{$nmo=lm(73j`&{fC>vu&bA0Sc<xtCkW7KD|P+zO}?Xmd5Ls{7o)_)AV+?!c=lujd_y|xsf*;l5HMWtA4*8o)NCM%N$779$AsyUu=SOF)q8(-CaZE{pC?iB}KQ02dJ&9#}o<s$?>k*sB9i|{x}!?qYO=Y$Nm#};oUD17lb(0)88ApFe=1sfXZjyFei6j6MJjGDKoc!L>xY%QS}8~UUteef17`$ZyHU1_LX>uO7F@~=!9Z?bX&vxlQ5g{JNe)|la5fsp2=CixH2uzBAK*!p6QPdM12#QaZJiXs%><z&UZp9M53_+kIRX*Rtr3*MEq|5UJWy$xDyH}fgjEmL4-J5X!UjJ`6&KzlQcIDXsfG6Ar;f+hVLae1jR4hR7`v+?e(e9bM6#zaktu)S8(?TFc{WSNbOw`=L?ynGePt|~hHhdBA7m=SmEUI|+MnU`?p9X)QLi7wQyU_RT7Map6T?xD6ZcoVwm-C>Lb8NEuaq?93?Na<;R)fv52?LBNOJ<dF#sE(FiO>SHs18!#A`~q+c7Wn$3LX{^i!wJtaACH1lL4`+OfPe|gA;a@X<tS}ekQbhMzsTLBT#ORFkv7JUJ#87Aj$x4KNYzU;{%>fnrXaE1rf^^nA4L=MT!A#4wAGI88Y?o~w@$rOgH9L?3qclS!eBbTvhyUlk=<pOg8#Z9EnqSd;n`>xK~!EA3z|;Gr?x5#t<=mf&hASaQ3hH5!8U6hdx6MU^pFz#k17xftYGb!ya-%#x0iRjp(PBTzNsf1P*Kyl*Mj;juf$zB3&=Q4U!+7oZjQk>%0{jGE$i$!Y_Ixt(Y~Ss#KQH#VB9P0oOO(VU*u@{^87HHeMn*?RVSscI2!-*!juM|E)uZ!?w`q}B7#F&6heE5=jq*HYk?TFAyY}~zuBoY;=@D7yOCN+&AfM+1^);{=J3P*h#2)(g5apVYJCCmef)Whlffc<9j!eaymQ~hIV7h{%g&Dgm>W~%wWH4&?<YfW#(>lc#4Fn6BuhJm-+j`F{3~70`dV{-%d~}R(crjUUcxn_TABneP+yCW&Tzz#<@1vqCy{0H<9|$vu(w)}U$&U)|L}GA?Lx0aTML&Jj0?7ZTd_0zI}0kPcu)N{|M14O4`s%5iOBSuU#GS3>V*kK?pf1ySf;8ZtM6}x7@#3%AkRmeW>R1k)v8;#%^GgQq-q@$f}i90zEu2buyoz9!uY}^FVDD)UT{R0hC`9Dwl?<PxG7J4`fFI+Z|kQ{$3@laER%)u#<&D1@)s&%Sw{#q9{94(%TDAmD0x5M?eZNIU$vuTn?e>i-y2uClk-tJc@*_TQ=9uUhJDze$Il#jAt#CIh&+RtbxbTwiP8f@#<s@q)xCkFJKF4}fb-PdV#W;sXgJ8}INzK}!aCDQ$9-G;<)}t`UzNXZUn(&2#Ax*7pKlMd6P&_OtjN7-Qv1hPrE~HWYtn#V=Bc=wGNlk^rX$0+(Y+}N@Q`j-spQE^X`Dk8|I<=SbY(%!TFHD}$5L7u&u4edX@07;aEaSP;Vk*ks%}iY)(=>uQX(98S<e##$r@IN+7nTb5<wlYT<h~Gd{Jy7Gr!jqVab)A00000<WdpRW275@00H^r0?d~l0(&povBYQl0ssI200dcD'
)
