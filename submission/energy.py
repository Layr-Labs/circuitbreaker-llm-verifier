"""Active -7089/20000 first-pass, gain join and starting-majority proof."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial, isqrt, prod

THRESHOLD = F(-7089,20000)

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
        ((4753, 3131, 3355, 3376, 3409, 2927, 4746, 6637),640,(4395, 3216, 3329, 3406, 3380, 3016, 4634, 6196),437),
        ((4767, 3108, 3375, 3353, 3314, 3362, 3413, 2921, 4745, 6637),809,(4397, 3193, 3339, 3387, 3289, 3396, 3374, 3008, 4625, 6182),706),
        ((4768, 3099, 3378, 3334, 3334, 3341, 3321, 3349, 3411, 2922, 4742, 6635),946,(4392, 3182, 3339, 3369, 3301, 3380, 3288, 3384, 3369, 3005, 4617, 6172),904),
        ((4768, 3094, 3376, 3329, 3338, 3324, 3343, 3330, 3323, 3343, 3407, 2924, 4740, 6635),1060,(4387, 3175, 3336, 3363, 3302, 3364, 3302, 3371, 3287, 3378, 3363, 3005, 4613, 6165),1056),
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
    ck(U5<=F(73670242301,2**39) and gain>0 and bonus>0 and 0<cost<=F(123066712217,10**12),'active energy bound')
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
    out-=y+2*penalty
    if out<=F(1,2):raise ValueError('initial strict majority after motif repair')
    def ceil(x):return up(x.numerator*10**12,x.denominator)
    return {'initial_lower_1e12':-ceil(-out), 'initial_before_motif_lower_1e12':-ceil(-before),'Y_upper_1e12':ceil(y),'cross_pair_upper_1e12':ceil(x),
            '0110_upper_1e12':ceil(penalty),'overlap_lower_1e12':-ceil(-tri),'11011_putback_lower_1e12':-ceil(-extra)}

PACKED_PREFIX_SECOND = (
    'Q8UWgJ${-ShXGBRGMe=~0i^r_6V&b)Ymug?9ljnc0$9Z(_u<XsesEEsM;ZHjATY(IHj#TrX$h*{(@fCAix|;xK4on@)os2?>Omhm2B3I)ymMc4i)Xe3ZqpvjwbH@H1EHby&|VG26Z`DH?v|D~<Fd4NaXk85HdK+qA|{J!xN!>QS@iry|SMQcmv1kn-g5Yjx&rZEsTaMLMVXl6)^`5$4XTXMs>~fzUC#D+Q~VL3A2zMeztI6M'
    '5$a8xGU!oym{isq7&p@XNp3nUCxo<IubWWOAA>+eP_h?qfPK#6<hETvbo{_I>oJFC@LS8#_D=2TSD~G+RBtA3x?suNaDm!RzUeWoyIJG%klpZ;y#jvBv_8w}fF?b<7-zd-b8poLJ7GxRfMpL)>1&5)Fu&Aw?O2Z-Tkf%0@^rRdvW3cm1Gj2!fv@U(x)?!H_RYR?E$CNI}^fT$ix<rnFyARnKBT%9E9lo5#);9!uelUQ$dh0w'
    'c0avypQeC0A?#yeGz;fTtt<pDcWgp}5m<bdK|*Ye(vn=)Zz}hPnM|CVg}j8(OhXNUzntd)Y~NdNYGKm3_Nf3buK3IxaxURb{-Qm2nC#439JsfdTs!;ockVfJ~Kml6(siNGLGNG6@7kECp)$p6>-gG~o7k`$WX9NJSu|+ze^lEt(n?GT12IppCangM(iril$jjBRvD8-q7t;)Obz?`I2E1`2y(sGadisvBb3RL7w;}B5P5azx'
    'Q28w@m`N^v?*NRx54t)<2lY{epO4F&J(<MoEJJ68U;dV?eQ_VC_2}wRjm`f9J79pXJP_5<+Jp-Rvb8bC4DA$die-#6Da9HJ=!(2njMVV(hZLo~35L4JwsVcT&Ut0^m^5LsK@GZ+U_xQlS3o5G$439wkOamnWsxn<6dEH$U_;^~QamIL%uyN7Z{E4sNaRm!}6%OMBArxF)o)t$lCW&8*>z%f!DbTt3JXx?2H|hYJdfwz+Ml6_'
    '+@1hCK}#(T2dx%7(}E9tS+0P+K~N%q%l3td&^BS?TK&1|vkB7;9d9^{tlD%|>aJXGSp8Q$moB3h8q%E+Zrz&LSTbOY1TNtgb}LXG{<64eSa9h{^z@Z8OOJXd52&&0rnYDD3-`3&BBp%X-XnY|}K7-QCN<pK%*38{XIB^1<@Mk0m@cFG|2895dmL|DI04N%mo_t0S*M22M%=PzaS|W>hST=QT7Mduna^*8IvKw&!~@LV11hLt'
    '}ddntP&OM>%2Ix=tlsv(lV62h%snptI;2IJj05rh4hzctpDXHag(zUm64Mm6z3DCI^Ikr^%6$l)F;_H*OI!i7493FI#ISKNjt)>>C$0+~1gO$jrB{W1Ac30YNYisAId8Y;zey(-6{`HhYeOOI-KHn90uhP*!ECWaV17>fcoEe-r=ql|ZfKfkOr&&ph!JGa?ymQr^i7Cg;rx*0g`gDza;X28<I<?(utB*Nf;{2UD+@B{%)tx9'
    ';Lbo#C1+A|L|*XeV4#VwSyVOdfrk+M}o0yyh>LBJtp+kZ@RJqyy#L{TC0!90wKj&@_v~Hir;xGbt+NGI^8H{KiO6DMCK1OVA2QZSL;83*2Re)%I7-wR!^*$gZ&8GR6U2meRfS<bG2DEluA??P3X}KBSWTLnw8xshJz{)G@7WZMRO}r4a-rlHFc4-zdR^{6!rmEQh!tbVP6xk-)nTq_;y9)~^(-aqy@{xM3$`7r+ZgT#Lbm%O'
    '<WGd?_@Dm~AqdP4*Ho1qYCPPmd~-1&n?fnH-|I9HE6IYyrQTz@#^n?PMN<eVv)KlI%2N^(9fp*F?%a(*N8J4&Q9YD=qDiw{4NjVj`=ZXRRCs)`sLezY_Pc0(JMVrxVuPU29*|Hj$UC&<^PP$&&bokNtkOLbNDQc5Hjw-w<^P`j){80F1=f^DBNYeG=oYphzEL1k_l)+nivWsJX1V9)Fspy~#1Qe838pzyvnoRdHzse^b;mYx'
    '6vo(eeoUx;rd&+j?cBg+fz6<8O2YYxMR}!jUg-td{`$8KB*PieZTTjbnqEW6oevOgv8JYsyyhPYptCp9)FPD`~%@xCxWN+_L#H3bXG75qRsz2+{Hba*gx#h%ZwU`9>D8$AwVMR_82-nQ@Mu)!cdvi4%*(bi{4z1mhN`p4ppYi~|UXv8^zsu>g!uGb{|Kfo2%O>1%koTFw8b8dCf&9B)G{Kwv}Mad`o3o%Q7fI`+Zo5CsM}&C'
    'z;JfpP-wsh!+@rtYiu%JCdD$8ozkx-^Bw5-&3*QC7#NN*k5r)*U!FFS?joAp}O)>lF>WrFrzchoEzdS36!sZf8;)>WT4?4j>;rw3!dL_b{5*(Oh$f1&M;tCWRZ*!+A(!S45UM18d(JAirZ!NImI#2*eT_Z~cgzP^#a35>pGP+zf{15aAthUB0(gwOd|nyNzRw-vtY2B!7V&x6x@P7@*W7LVmb#K8bO5_wE3xo>Pua55qCE*F'
    '2&Yo&yMfqQb!BEjG(X(&;Kpa2KuwpO|iP4!1A(YB%vSn0+{S@@|MwVEL`LTv}PNxos@%*#MJ0T0F~G-|&Ul4qa8yND0nZv)S7L7j__2pc$1Z6@CEk@CUAv{y+m3=c6dI(jV$0hFiOm>~=}SCq3>W&;M^lkf4fR0iocgrb}g;b*bM9(fxpdMhI~@DjumPC`YC%zD-9&@2ROxGfKZ#rs$D$2KIU!oHEUvb{5KDgf!5hI6$f!(G'
    'T!CpY568%|4SBw^b6s9M}I=YGBHY(|8Brbi{CnVaHUBp2=Mn4d90c@FwnQTtm1Dme(1F+TQk+0C2H{Q?AVz99B?sc0JB<mKr9`(cNL#js)U$jPL3G>MMGGG99}8M`Pi3Nr2u4KmMT;w>M2?L>%Z~wZNoL>eSichQgCocGV|PbN<)#^Ake|YyIu)Srd{AR{H)hhFCJ(=6T7L!O8u*bq{z4`g}GdH9egQpw1n1S2=K@6b!BykT'
    'zSH$3Zf&2^n&>q&S$88n{nJZ2Abyeb?}~P`LtK89)`yCg+*uH{VTeN;rV}wz$*-bnOSnl5kRR!rf&ACh!eG_sURl7*1RGKH>3$xk*wgVZKF6;>Og+B~u<)NYZC?Gxyi^ybY(Fg=`~i>{HBkq_0}vcf@vxd7@IGW`q2+(LvN>0KN>T?TmEV>Y~E9A>L`nq7t;co+|}X3j#plJs__w-_TZynU_qO&w29KvS(3KtfQd!yJpyNmf'
    '!q4M6FAPzTafbFq+M~H16g**cwOnKmM?#mcUxI;bJyYCek?ffsNV^<F&0x;OU|)V0sS)fC3#$4cG%#SH_4RxABRmHUa4uBdp<2Po%CGYrij^_HH&kFjaPTCRlVXaE@Q40h@;sn}6vKeHHBBSW!$;PmhU&CN(nObj8|J367C|`IspCHbgdI#bB#^apj;0L*6_>)p|!PPG+Bw-4c*vs%BmWEaH3@K}yfmykrm8p9;7`BE{J`F^'
    'aYFP916qtLP7cKKa9q>?CJ^5p`K=`)|6mh!8@wKYy(Vj0Xp2qL%GIb@)G%`!m|G7~O1>8kYp%3fr}r_h9$QI?*=wT>PeG6>yjQDt2BZ4V}f9oJhiEB9)Jp0p4B2#ZjvUz)mrAo*y|^P|I$fuJecyT=K-WrzPDRn87bsj6>2;EkPKLIC6dIOhTYz@6W+U!z4ib9pd-fkz~|~FzI*<iQs?9cFMOY|F{<a_td@(uWk7NT4d=#>w'
    'ZI^@QSYXuehDRGy&Se_M`oR$*Em$dM}<Pnrswpe8kKs&#Rb_ow67bb9vf>RdS9<Kv0+X>zSTX&WcseEVfoY95RErOe4`^UXN{Rhj_0!&BZo^gCyWJ`J9`r`GFDuKvlnzc&x{pwS4JIw7IFjIraV2L%Swh9I?AFJ$%PkWWb{6AADIZ0SP*-^W_G?^|=Oz3B_<01@usZG7$(nB`HwiF1m)&+l}~;cc8E>wFIOE<Q73sH<DnIwe'
    'AbRE05*w^Zy1C9|&Oz*2O8;h@i0mV>68ErB~gsg_49j4=c5PlL>O(3iU!s4d&G{08J+C>x`j49(~KI_1U^wB8-t=kl-^0Jc8hsCIER?_kl{1@d>Nl$kSdl$vBe8OugmqsT6)CUrn^$1H(sb@i7Ymlf(9((;abL-o=F@&i?o)qVp0fpjm`(L!dq)H_ydA+)@YTkQ#CJX*R6A`oN`#H9WnTf{oI*92RG=Qv;<$-pfx}IzvN<pI'
    'Fni<MM)@|1$>$h(5L!&yh_1B3xfz*`0%JsaqzcE;w!VP2o2Kh=;3_E7^w0yG`Xoa8=Jk#PO(-jPVTMxG#IT3fADBs7KXn)AEG}+4CYP;L!kWY>JMT6ZXL+_@hgAr&z?mM09E(@x&C!Cf&9L3}ho9k|RPQ;}2TI;^E^ib!YlR@Bx_*=DY)oNj7@_n?!?sG`Xg!vc!rwsS`LY5UoI~iG6*nb$;;k4|)9<?4ENDUl6sJ)n0#rge'
    'TfWbV5m94!aux<lw!LDyb@9O06~S{5jvm!ngDvtwucXz%beOQ+!jMlGH_Uo3Mi}A)^TOuEOv1U>{qyS?Z&}Ea2@ph-bwt;W;2-5*EU|KmS|&>9Zz~VK3bJmU%}rUUE$jML-k0Bk7Q}vfqQpy*FPdg9(bbwcaKTm$Uv<e})w#uy@zHS=12C=%1l=!5VV+eQsOw^npu{fsUPq(&dos)=n5=?bfRP)pz;EMsZCS*J>WNG+IL|!?'
    '-acLCeA+eKNxya%%X*CCx*)G75X9lfy|hL{ry0MRVp^<#SxUgvaRdy>~+*pn1Omj`s!rG>@6h1T;>AK>+T$RfX*>ITr4*8dd^j5Y^~X?qcr}^30zdd&76)=f?|TF04mE3WgS3tM842TA;kSV!{{eW}_4I9hlV0DhWw)Dzmm+?T~q&-IRoz$rkHacfQ=)d7QtGK@;O}w}o>haRJ*(z^?gfelo&PRIjvD3uxK7LG*q^5G1u~js'
    'tSVbF)9hel&dPaF97wx`gPo2e-ka%H$1@Zck7jw5=-JfP|;V9lYV`Z@706npR|mkSOxrsL}+=E3&e{sJ=f?7zM%WsyZd*w5@&Y1EnPuA9}%y%42Xk7g00>a6$&E`~i+Bntt!ry;p-z3DPYwo7TM4FE!NoOPHI?)_VgVl(KEuMgdl_%Gic^?y<rOHupivK-1r8_m&_fardZK{MUPmUwg>4iZrjOcSk&mwnFe8$braD5)S%;D('
    'tw;K#kXaVP3hwBFiIhtk<MDW8OXtV-Yty2BLOZLrRcJ<veY9QEjEnQghLQpL7A=)f^V;4!_w^GvhtNrS?)nud@P7bV|Z{--Eu1#d-vq(EnmZ2Q$k9^)Y+HK*EGyg_>H=UCaal)bSK`hp5L1%-Na1ch9N9-GC;3qDrJJ1${Oi-=XpE$jos$XHRD03(&@W0aaSJ`xO)|+J<pb6b`Wg*)RDeu35&AOOj)-n=N+7CvA?Qw_C+V3l'
    'vs#pQ<xCGf1G1z%lk9b)h7P8!kB)rA^N93!yOoB%b-gzvrHwXg*ofyPe0ixkwQIbVh9|(t!)Xsf-odw0TY)x<EyuH`LTH5HK)u;CJH{zvZk>OO--0S7&K5K4ntF*m5~-!k-uFmF_k`Zvc)Lat+qEgCWnL??I`q{p2%S-~>f_z6M{cq-`_}jjR4BvH*9I(h(9$f5!FR)AO29G=4JHe2F=~r?M1D&jd2^@zv7&I!=832Iz{fS>'
    'sk2)Pu^+VzbzFVwXF5;>It{vG$3IH~H?!mq6*R->Q!Ic4r}n84!0`v&pfRzf~MR+`ALiHJZ+XMq#rHI*~x3X0XkQO0IJ=!+(Yumtk}Llg3NtE$`QUhgZ@mfKYthvDeIvOItrFQtFBGQLw~~tL@pH(cq=7<(7pto)b#0AT6Og+QyCQNOL*8%<55v$_%eJq-)*lxd3!EN?j4=V+nhH>CAR6K_6Bl62)2|pl5xIG?pdAiKZV)wR'
    'vaPz&1K>!xRVbcr!>neZ7y}b<~U;#j3soOanFe8lZBEzle_w@hrxPN|PIl%)<b;xzrd3zaHXw&1?F(^jtAeXuYkGlOadm0f*(@4SI;!WQO6Q&Yv543vP~ql#Sel95Oy>;vUIDv46adtuqCi%ITv5+2wAq+nsLxcUq47_9%>49i!||@aD6dPhV(7=GFbsw=fa)n%ciw9OSsF#}LuZ8Z{J9YZ8dSA<rtNNfs0kkCjl$RVU?fgF'
    'JpBo@`Q6#kN~Ktv<Gy^eRW}YMVvsuc)-NB1)s{JmbU`vLwn)YUB-G)U}n_%$+B9B~Ik%g}2zKism%0nmJa4L$DDq9P!$QNE=%8J&hP21<}q9cjlOWfgu}h@a*>FV{*+AZV(7|=A8zv{=yP{CqX&=r+WJ6M?9R&47fG1h<nPxuSckbG_H9Qj?jH|_}bl@oU-HI+||r*RbO@J@Nm)&!T#;EmL^oftbtqCKd^RSPOD6*!@y-O6N'
    '8?yH+dwqIV{{<xK9AA<G!-tmj@|z8O8V6nQclWRRYx2UA>ykkjr74lrNxIhZ{*XyRnbSQh4-AlZ~i|_d90$F!(sXF``Tkp33~$-0@}GB96NL0F+9e5^SaRUp{IrOfMgJ#HM%U;kvJ-+%`WZefK&vOhA<OEdXpm0O3X9GkOoI62Wx+%$1{*Cod9M{;Pa8PUixF_@%O$(@|oM9?PcNctdGjEHG<Dt8<B)qYjmpZyB2mh^8fxNh'
    '0wKD@>^Hn6v{n@4;My&r7WmTyI#xv_Q^BK+7X@OoDwrg6={fW9v*Fo(bc?RKt<XX>wUohA8}zwnWI3A5W<mSIEWSrF_%D7$0yVST`;=_L#@1Vo(05^$6!Fw)8>%M#E#$IW1Ygb%*z_lUIYBIQ;GNTp;}ZH`zbGF(}W>Hr)_rht@1|@|{!OltJ?(Eeh|N<D6v{*w*D44+9V}QYM!n>ACmmLP|eO<`V%+SaAa^P7Ot+0ma}{A='
    'PvWHUnA9INr+!ADM9Tc|ah!a&mT*^a#5VlUWb1b0oZ&@NltLY_x27HbA)$7xxb#er_mx7{nlLZL-AlU<^rmI&d>cUGBU~>k?#{ASra4?25Dfht&^MllkzgFK}|r*jL|b#{PVqTsVnh|3xWOG(C{aMhYt;ddESzYVT%7CtG{ChaACRc=1RM@%2i>hct-bdkxhLEVSI;8m}^xpctJuY&vJbnmA1-G|eR7Lp%2eJDUIl-VF&bwh'
    'S4j2``cd-fIS>fr72?19S<c;`1WxojppU&kY_|;QNJX>CDHc68k#;LWwt0e%ke6d}4QRK8N8yjs2-S4W+IQ5KV4D31yG|RfLau>@o&H$`TdBB%$;IW%!)G;lYhhV{+epZ7w7pFGRM*jpBtEt?u58ce2#U3}@#xWXEP2swc_)CZ-weT#T!zk??pF_I$Y5FDQ(EaiJ*z=9nllP671A9FRSKk*pl}V8^z7O*E1WXP^_xhjU5#XA'
    'DN@MJ>!#>L99?2eNUaImQ#`j()R0uuS5N){Q(0?Q6$K6&9g!pzws{s5Jptd7?vGhK|bMNQI~NleQ@QE#$G4?|H$?G>U@{TY!uC%4JFwTON_4DtircR>CubB0Lx_n*^*>-_yYp4_BIp>zJ%1LtLwPn>33=G11PzjZ=tC!8m&|&QA@`qQ32O{WX*}S;0<R!^rD{(giGF6tl4Hv<&bUP{l(O<-<QU3obYY{!#PH2cq3(>dJNZmd'
    '}TP0TYha(Sqd{@CA0|iEtq=>p+Z)6u1l&&s>A%rS;$3oZnf|ai&<>?*X2eKzlHPqvK^bK2jV%V?hJLl@nRR>1sF4DkpW?U<N40?edJ%X44ZA%Q{EhuVKoOmBZo5kn3tchmRz(hV0oC-w<-Hm?7@gSDzj9nuzn^MUGXtzoL^gl|u-E+YMucVK#XJR0SdIGt+%STPU8~slB>ImUM@a1Eu1?B=6l)bnX~bUc|#07xG_Nb6XlzW>'
    'wxsMr8DoAL7{#c+6una?zw=ZAb)NB{j);$-oxl%y18@Ing69i>Vm^H+ZIh^ErW}tE%x;LkC9s(_a#WTd7t*nc4OA-|gn}XmR!GmY6&}Jo#xVnN-3!fR@@K8Ni{hY$rC<uMak2jPoJjcD{yBy11^&HQ5qHe3V9URB{Pg&_RNH{|lTB9I6q`u*Q9|*Ds}cc>ix1?;I7W*G>`I9HJ$0lW^Mtc8hC13}bRiy^?{lcLXu|b<$_34q'
    '6%6jPej&fcdHnj4(3~nc<PY+u1)sV-MGEgd^iqO_hqFda)yAQ#;7)SxX(p`rk&Il@&!RK^FRjmLn&r_jkb)6&Y)#0F#*>6GhlYkUE=iE9O^LfM^e(YV-_%J2N!nlO~cq;^EQA%}M|Hi_C93_k#Zn6~c(SxP<)Zwb%A+W}&^gx)Q|hx)(HB!?dF?IF}tMLv%sylO63q@@4U%p1l8AOUnXb&Bl@IKv$N;b@vp?;ihx+g-(;X;i'
    '5~(Is+G|TXa6_6-#^z{M@sn)9b~LX-w@^?-i!rEEIsz0YEB&&fUV?xEgCJbbnrey^coPKONAGB-`QJ`kpkuW|2SYawP4y+{ZH2?Rf6(xuI3R#3XNa0wflpFy=#V(5ewdB%f+b0dc-GU8Y@Q_+w4$gYaXc*YwgZ-)i<~D89J^NLbK9qhxu&P@r6!SjIPio`-rveXUxxlL2}=eciS?K)hBEzQi6O!o3M88CDmDq$I&s@~Fa};|'
    '>96w@#;i?_(Y5k0DWmCQuXpMZI18AwIFHu9&bsqem~AEHdb!Nrv%5bO@1ClAyb@%)~A0Lt>Wb;8$@+F~OHdY?B5(F30y4R=&Z_r8~pgMIM`z)#!_xh7jNF;=Y3`ad^iQ#3!0$X}^ed;OnEe3JJEA&^+{lGPIu5m(MMiy1Ky2_LJxHREBk&_TZ^!xxCvV*sv273{xMb>2@W@r_UWI<j`hwT$k}dx(~vrlWWM&Vw2!o(@VFMez'
    '<e^!r{#z+>?w4od;cnM^6)jQtWM?+(N~5M`<_jRb0@Zp?w|H0>}hN6}@>&h!7#-RtB`zG0y80a7!*zRP<I<M>SU0uyPGDJpoC-3fc`(H_yEj@tRc<R8&{LLT_=iSwk57vr$dgi~X0?31#u32h4KG$S(>8_=Os^aiHc-4&|?_IyMu?7(Tzq*rBGE6t&0aFZsTVtU)+Dv4iCK>mXc`I>YiD%|WNE6VX!Hutb7&?cvDYwCEzzo&'
    'X-Z%2}^RJYLx+=5R2nu!-34gk(@)JLF}D0j;Hx?Q2Kg3V4t1DjQ`#VLt_*Ylq8=fz(6Uk<N|^o=l#W@?=Wn(VuuD@HIlLdcAD8D(E0JI?e8#8MfBC9sd4PHx%uiBXnaCY0bw?&!zF^^QvS2pVU?jJ?SO|FjifD|BMXM_X-{I6)g2iMV%(Sh{7Pc24<#wtmZ89$)xem9)(sPqWEqk(EeHi%vYrgahrC40^J$~UTWZ$Pl5@Hj+'
    'RXI^SS7<j6|!-vJP2omXdW5&(C6EBr|h{epr_doK7NkTWr=L!6DpkrTl){{W_kRh{q}vm6Vz&t!JC5)Zlhp(RhD*d+-Ol<m(U1pk*69(WAEKd0L(~Fb$F;X^$()F+%$0OGdDI#wvq;L5pj?CDh(#$j~f>s|VPlu=Dkd9AovJTqD*!hCM693?(cWYCiB5mZ)V&nzsfxmj&K@E?M%~;E{J?B1Zy(7@+~n0ITMmJ?P5Uo1#txyW'
    '-$FeJr~S-;PU_5X3V(<tXpY6h(NRZ~FGukY{6g;l?L)(D@Se;C|$n5uz`N?5U^0n)LUUqdT}(x}!k~Mf*7DF1+T6;eNFC@P!Sg)m|A#uMy4m^X4)&YOe~%k&N<os$x95$3KOcW3b2GD{ZD*-{SJiD#lL2b4-Lv`x^{wDN67Jgve9Q-od}BB!*Plb0}_!8+S=IU!?!VLiGv{eK);IG~sUgs3l{ClioENNqW)`5-8)@y6$}llU'
    'r+p3t7((bd<ujdh8ZCR|cbQ1%h>!a7(?<U~TTdz;Gyf#bRuB0HB}$Y{_zmcSmfI!<jq8Fj~kxDDp*O-fV@vfW*e74wQkE^aw`X90_Ku<XRU`*LA@S;+1hbDOR^U*o>u#J=cmTTS+INzNBjAQ|glOR)E~0<LsKI;hRSxtFD4e$%ov$xeMBkza9~2BHH|?k&CtPMYg6Bh&MXZIb9LW$|`UGQTl%Ipl^7vZ&#ZBhv8_~RYv*tSV'
    'YzJ@q;#vuM3VfK6b!@;js;wOVOcVvZLF|xN7ULwxBQc-gv?bgAuQ%qzVZ{$<QtXKHgaq#qwZ7^uTQ5cj@1G;XSH*Tqh+CJs$==7rVZIfT=g2!kONhF?!ptI)_rZ?@uTP#3X7<yLoMbfQElE=QhXaK|vXc`WL)W|LEj2%&@P>;gxfbj4gDb;=v7Gt<hEt0UJ4{bBk|5$b3Th4w(qfKfyk^tTTNePSkf(4H)>RlZQi~U2hnh)2'
    'hae10(U}nZiIQ{$(2~uuKj@EFCX67YXLM=U$_;Mc14)P^jlXJl4>8%|c}8F`m?F8LoUwZY%+7Q?^4?s+2CWOW16xx>?5UmdhU8)9SQ7ex%I)`hM?9Po%R8xo<A3fmr9iqatW15ITVAmX`FQfoCs7bAl@R1>K<2+LRqlQtLr-<N*o*rz-h%zvRu7<7iiQS8ZI;cT~RIMr+`QZ+t1`RtZ^g#_7pO10i>VE7?|-{B)E2l(+-KNQ'
    '`p)16<W?ps!E$n0^y~bz!uw&|VER;Ar$LCF`_)_)0}sNZX>N(Ks&Eta{D#?pr?ZggnRDw&qK@bZSi6${U94KdLXL0(JztSCUu-L(q3i(Zu&H>>xG$gq&V%Nx(RrCKS+##LL~FM8J_lb(<D(jg$HSz)=&uGo46Pn9k;BPGwb;YJP>AP>GFTWWlEf$5Q3F%>7S9|L{I+D7K-iDI7fy5bP8Jeij8r{mxJ+5oh(oAN!JrmtRoNsI'
    'VtD=s?&%{*lc|@CuLc3;pS1SDgsQk#Z$oq!iU5D5qhN{mS@iM)w*J!QrJ3`IhMC$(_b?p(LjEwL7`z+{E^dL6>ez##JQ{i|KStg<pQ<4+8Zrpi*pL_UJza&a9b2PUX<#yA)Pc$Lu2QjN$1>LvUk^;=dXu;b>!4%@d%j4|yq>6M{j?^CSB!k*RFOwC(B@C`P{uMfd>>2A!}1MiChc?)LLhO&T!RWrpXMdrf>oO;E-ep-~?mHY'
    'U4|W`wlj&>F%JiRIHYBhnye7a~qIZre}Y4Wc=X@lR|NKC^kS8~Dc2$(3=z(Djj9{L0f<Il36{osG1%pd=k?-~TbxS^?(BQYCf5R|JWxNPz1uWFoReD2W)WnEMI~{LF^ihf$7qyWw1hvGJeN&9ngGU7pH&#<IApI86|5O`2;!F{2RA`>6va`dHl9b%zYTVS)<$*a7o>b&-HsPxn-cbnJ3gJR<CeIE=fl5_Q<EF=K4P#7(9&>f'
    'j%4nrohSoB$K}1Y9D#q@0lJhbh)k*wN_uNW(xBd0<t62=aM&-AM(OI+5vYBc<sz8S=W;vxRD^zZLO?e7rr7{v%un<az11$S1W``zRGZx|6bZjM+fgu(yB9xahP%vxCONm9AmF9iTqR8uFnrnDzuosQjs98MerRU`W~56qxMf?vya<16r~TJBn{;x+6!ps>6x8;O(T?e~GIfGb>j#U$_{d(uodNu5m;Yj=kBQ=Idt%Dhr-*rI'
    'R6Ml&gaNUZfm?W<0BI(@1U*_47NdOs#l+!7VHz;qi0b59(Dxsh`;oEQMR1;_U@wEe_g}&BXMM*&n1ZE_9^X4=g%$lf6DsH*RQQVzDwEia94-Hcu|i^%Q)7n8As&dO6qwn8F-yE(+J)ezoT|#M)qovl`_gz_(mhydp@C^nFCuvp4Y9OIqKb+C-r7GLC()^Ee^17;-=kv}cnLXqDMk#El*ER+JzO<8bQldFh)N@vaB75cj{bvg'
    ')w=0Wr5Vx7nm<-2eXZR2m0T=>y5YKt{|GMiIkWnlb?+xc3^OL`7PG%LmZ&64n|4l@Sp~4=IbqPk+EkZQeWV6#|m{Q;R9m$BGhQ3LZUX-b)kfHIly=wJLxsq~L{{+-eTX?Pp#~Joi@;++Jm!gfyf1p5KbLlV6yJ7v$m73&jR^hd^F1=*?97Bf@e<A^cTD=cCRrYe_x12av9De84swU)%O|csG^{n>YZeDkL7}!dqU!Pw!skqu'
    'Sn|y?|GDV%kCh6di9w7}<ELT~2kf`VNc5-Dm7_CC30NqZ=<@UOZLu?CPKgl2Jw6TuNCnU_S&Y&8aOLMcKzl`<^MAztdM#^ny-l0W*wvDo(6R^7rAt^lez@6K0P|ppTBp9bx%03B!W3JxWg~aeei0BSQkI$VlgZ0iETmzT)-Za~~<lfL^*1p>`osbu7{xiQO?W#{K3p{1ji}t8^wtQr*m&CT-WZ?=D9qIEq%#1AAPp`Fi8^1A'
    'YAbzHop6J->1j$y@PxIzd?Cqrlcz6iiX@6+&K?j$+~5B!QaY)$Q>Tj+0ac;{Rxh!?p#@NsSOp*>Bv8EMVH_&r6+sUwl6L6iOW9@9xRSyc5XTOASZ8y$uB+jM7DTwIZit3~*pHxfBQy+hqREY(up?6(#W*@)AD6ITE_Zw!V>DC@G+#zuTzuIsvPq8SgG2;?+F-9^wT97eTbz##W?@!|bwDDfga{^iFPY=JRaG$~J^2_B``8!('
    'iwlPVeK-!Hx$;@QJ4gmvPjSE%U_wJ-#B@t|k~sg}sRg=oIBnY594rt8D#fh&$}kEw<+NFs5f~{sCKQLb^z&2S;^j^r1FU+$@R6w>b5PaG)^qcBP0td<svOE;bU(`_>6*kD*6Ht&R~#Jo2wtEaF$JPUQyhM%St7rpONb4-GFr`eZF5D7<h_mSf!RCS7P+(vFrHEgjMwR+JTq&wl2xW6m2!?(FuBR^aEHzL$<LS;(six;{6>P|'
    'K#0vx&U21dRAHFw~=nZ~}UH(|5!79nR>^T#K)cWU|4ZXIhs)8gu_yyWoB|WfhUiA|x>V)0$}5D~OI6qikw@@K`;A!QRHqdW0>={Vg^$kE*MW>S8N~Kl9Mit9&MZ2Gct1+ur8}UhpP;4UMMaOBt@n(dNA)8}f|z1^0z}6-q<kkbAsgwi$oMe;k`HULb?v(Sc7^?b#Bz+2o{y$6R~i?vucNNB++STRTz_(<Fq(rqR?up9PZZIn'
    'vl`pJR|w>k_oG{$E5)PBC&iJ%awAkO9HWHf@AQCR2d5m@C~Ij#9~&QCz}=N{$;~qe2X$pPXdBtq1E5gxbm(__!lk5)JW$akPz_&)kP>)n3BpKWT>A)$ca}7rt!!_LQE_IR`cnV!`S5_@y8t!*JrDwrO)<N?)8Gh9!$@-&fa<oOvS;_PauIGGV0nV#cxb?fNC)O<a;qMR=-5Q+#${Crs*|$0xoU>7}c`<Atjc(|xat%w?VK`Y'
    'IT?Y8V6R%?Mi9JVvhLI<9-webvXbIc1B^bE0(6-F(axg!d9)`#5rZ-vSLODHkZnr4-ZrSU#YI*u*gKX(XX}D>d)CKUJSU3b46F<u$$G(Xw)3iitgs$H}uj$FgDA$t|g1?8JT{zZLZled%G8hnv+2PS&u$it2Hf&D8Tsy57@^16Kio<--w$Wq_}w>w+-t-K0LT)32e+2%wEl47YPucQCr^lq<Wow4^2uV01gvb1i1wJhB=v4x'
    'R+leeNOt>~|_KxI)W18~{Ekga{N4Of$X-6P+9GGp4{F+nOVt1^r&F)?=NA9BquJ|KCeXkN0BifcS`{QbiV1Tv%q#6(~}#)s`)w5#<)RdcMBjw|}!jREc67KHK-b!rr!7j^ZMQgB%?o+JoGtMv)%6l+yP)TTmu3Z0(gKxp!jiq|vI&)<ZN;K9Ki36t-P!{67XEWc6U*mzj<dwgER8Q@d@|q&T2TbP)Egtg7^6;CL!Rt)hHb=2'
    'JbBCCwf?c8C0;Yp8XR7?-gIGjvhP-UL=-J2r(c@C?wrk1G+afguknxe8zWo(p&6uVR1#oml~-f;0{n?Z2Y6o|H@4OFU?qzL?QV;dR`czNcoBoWUUGw<p(IX>jDHz+q!YJ?(7%ZF$TznN)DqE&x?_{LamQ#gx6^r(s-ch1x*tNQ>t<G1dCzopafGi6hLXc(>i7JfaRLIe>k&w#83ldk*XuNA26eB}+Z0`wlZ=Ux1Br`^Madwj'
    '3xOa=++t!8}sA8w;Pl47QD&MD4$ct{T=2;59t8lQzr-0IIaW<-iy9h95EV3b7ocfGilySE-4L5khSflLf3uV_rw*$d`_LNO}gIq@)@?-*s8_r+)MB0q{N@06y!=&_(mpfP>{ujIWBAw5y5p|1nPS@(&+wuf__+3Xg^XH6yk>7CKgO)w_UPi@l(Y@@qRU>-1Aeb%z0ApXvqLs`l5&faL7q*pDl?<=0!YKn$Ol%0Tgfz>Xw*mi'
    'v)lA9m91i$V&AwPUE@F7T>=t&Pw1>O#rWeXA9j?P?fw4u$KlhbPDEw;8zk4=H}sUOIs!J#7tv=XJ{sQ(rJDzkx$uK0XV43T{fDln29Lb3fVNT_drCw=(BAA4qo7g3~|Q1B`|qVPq#tq8+2JXZPOR^D{dwCn*Hii5RKSnWB}=wOi^PS;Rj}vPCuHL%{|>*(j^O%YJf+fh6gDcGvAQId3++;j|agm96jl6Nzhar)Dv-5XW^W5='
    'IBM%~>4c7a1z}$S(3_He5<L$0}p4&MmqH)gXv*_bnS{$iDTV*%2Y#1MtAo^g1&vV>ry+1bq+62gtf;@;NPtjW#wKmrSt8TVVW;8}yG3$WZF12Vwc+7SGF0B*8{p%t~T81b7+CtwjeJSq0iusVpSB5WRg{KQWo;Co7)cwh&`&mnmaUAM_J)v=8%Jr2anSv6Ta_gIED~i3dgDbUQx&;e7TiscoT<iix;2<KiJ(ANurVeyiNTdd'
    'VimmxVJFI*o6~Cstq~beJ4)k4<e7faz_pqXaqnK&HBddSH9l?ZsgQ-Mc^Sas7%bRKf)Lvk02cgxEIlFo=#!xYFwVC==qYA(;l1wQbC~5SKNHDM8}{f>N~ZcedHNIS|4)&r}y|VCk%%I5;mNvC5yvPgb;LlMgbHOBG8<_1P1XH|9F^dl#3jqi&?s@U^jhG&_(h<V{Ejx0VsG;DEeC%vF#1zPpoerXchfRo=q8Ay$j5qj@x3PC'
    'u2+?6}$bZQ)k^WR@_#;K2jJ%=TU4M%o1J>W?9vW(E!M=WnQPw4N5A|66Ub5hNK@E_D-j)&2iXSRM2c(l%X*Bnu_dS*$le|LQ&e4&>S^hZ58k+WTr3!oU7L?=jE~Yr*>Dl-R~Z2gOQ^($%$PbY|ITd!J`Nz(;6&XV_+`WoFJ%<{z<>xaQA+sqrdUTjz%t?mcTLghxO9TFy+?Wi{*b4$&o?{6E4q<>bDE-ypAV!jv<hR_FWK6Q'
    'mZyEuK!8qoctly4=soqg`^#g?XC7Xrc>WbR^@opk^J6rSl4y@<l;MS_dF5%=uHniD6&kMX9(yo$<O<c@1bi5Ql-#!mS0tPcGGYp3w$JnjWJb574RRlLQ^;+_Y5M1iM0eR%qwh<=QVK-H@=}@prd=X3A!nE?a4Kt>tJJ@osFJYnWV<(&nmH4_#v!c5#<2cMDzNRYRjQa6Tg*nE(g@FnQqFPH`}9W<*)RjkPYYmm;&?2t5_1a='
    'hi*NiCp0PLp8%i_G8H(nHp_2~ZunxrxRpZr1*4fmJn|;g(BrP`8uu&p@uasl0i*s`9sKzs4@_qCaTQDBYtFaHCp+a`Z<L3P#?QIOt)<z*Vn9r~ffod9VCpU%4Lh7PJsFd+ehGoRN~y(!$*u%T6x&RN?{le61bknIFG!{@5U3>Ln)6F=rwTUAYf!faDO&fQZrj=F%)}Qe*lUZFD@r8&P=8n>;IyJ*XR^m8!DeaB%1Fl6f}BBB'
    'i86IN`NmOR&Dy_&#0@CbmdYnYmoa@0dm)3_LG#d+|)%Rpc_<Oo}yGA(6`5VKO)u12V*O>~NvgM$|AK)KP8hIV;P1)a(^S4S0N4HCX<YVlnq9|2=f#V@rbn@<-r{bz@$~EDwuy(2>=KMI<ZsJh2De;rUS|N*0zBKY>1Cd>^;(HpvalV$k{oZ#BqANw1HJu9<h{%@z;#7SpLK?kDBxo9x)JPp6lA#u$dEzzol+8IPKBoZ5ig$n'
    'HSLSGLkIuU#=OC-Yy21xu4symAx9ep4EW2iWG>5Xs^kpV(N;m%qqHunK@YDua`7tD;|gy5i47IegwF*vg?IHp=pfFol8ih-dQ$*O(rVcdeER&4;}8&`cF%RA*f-W{B#EjHIl;3#}<(E;g(3HCk_6YA`|ikTwnon_K~AngV&N6k{fG*3Ul1mQK`5m=D8S_VU1G=Y+?M8HG-2+xitpA&U)&6r<_z)~64tng-OBj=051MYWVbuO'
    'IYNpehQ12L;f1`JI;VH%}dLYxLwp>`{Wr)q^40^Eb4-@Iz$kW-dgoWasL!Yf(Q94AYF3Dj{UlF-^G6`?H$WOie$`8N4o8Mpe-!;xM2@d>WiO6|G!zsCT=Y0&9wPx7V+5#@6KT*H8$2H)M5!thX`>rR?s#RALt-D@+{XRr;f8w6@|huYJ2(xa|XGytktPZGfoW4^s;unbI|`@c$gB5L(#B{quaF@qoHAQ6PEvkcxXPRS1xd1b'
    'G~?K|O&JLXVKwjmHv=|4yFJ?#Io)O;%4aB>y|mspE<&AUEde@cP1uzU!y=2;`}6CM_n`y}WjQf#k7NbNe4{cfrY42Lc8psZ;~%SQM5$Bd{sZagwB{plDdK1p|h<CFj_=wYuvQJ4Tl(m-@FHo>Z^IJj3e|zwHE@f>k11X1N@TW_#LV*l7_PS1N9{^Ap)x>@8|(91ONPXp2w5zV}I~78r)X7Ydf<R*KC0>*O!XF3v-haf{Phm%'
    'g&_&~YeUL5Q*U1uhD2H$U%_Z?s=xy@X?`osA!D@6POFtL83p%Yr&WAcI{AY{qTfMv~l=Leq-Sv)oue@8>nVO0=JU`bbTx^oySDH&qI9<&2u2@bddWFkw3%196}{K)S6_EYUqiTMGme56;GneZ>TDr;c-*WonO9YkH&n27!OyBwtce%Kpblx0}7%C#=bZ>}zE<#1CT;HZvoroIr?k$a-S-CpGh<(3ogpez>l)(fb!<JCuu<ug'
    'Sb#vs^ch>t~#I)P6W9oH2)W=6!kd-=^j>b7C?xv%J0{A0Dr8TY*Gl!tvO*5(=nyqP5L*;o<Erxbfz=$Q)?eY6p3w%kk`I9z&e7TXus9`*fxNESG_zMPA2}$t&zxo2o^J0yX2gy&Q*Dk_p>V{6M18)lsl-ZXHVle)`qE3qnHP{F-l0vJq`^!Eq`ON0VOo9?vIL{+rtR22#T}VD-?iK)c{V^#6Cdp3#ijeUga%cXW|;F1=0?Eo'
    'f|%St3X3vbOQF8EH@3msqf17D#(ozp4|{_jO)BhlyJle|R5wlz3>l!(vB`H&^4HlYOE+M-rQw?MyezY_*6?i+j_5R7k5)K|xgB#t}YoU!J$NwD^^MzJF|wDWeeH7E9tVt~jz?&zFqn>^l(Bs0Gie^H3}?7`KR%5V5?V2dn$V43|d5ZLb#;nHs4aU7IMylQb^;a0N8USV}Hns-_ix#Z))T<e2;&SF_y!cvEwzi(BDK#U8CZ4m'
    'Y7f4-7uxCl7?bn7GipOGGf6t?izSbBO12pW{w5<N8u?5)Gk^B{aJt2!P+At+eJ_+&c+hocu=VV6A%dv445!VHAMJRbS9x7ORk?Gh|{2j2jWh&-p{HHU?*h#Ed25RFCtM{5v(VV~K)6hA{+B1Y!BS>QD>%WaY64rOw6F`UfA_`h4O4XFqAdF>NQ|kwF>hAjdj3hvUvz34DYmQWK7{G53zOoh2=zc0LsS`YhW#Zcz7BK<JROZS'
    'e;x?k5TLr*|06=ew{(s;#2L1c+Wso%9o^m-k15OT;(c3&i`Vf-RCybY&Jy;)3s?ija10a`a7C$-oX@?IUy5kDJ3n?^ZC@qGGL&^KMGyp{w+Ds{dbx8;Cw8)V6$B5_v$D0J2y@v6silGT?)(L%c(J!^r5LKezWrt}3d9nja~$kDhmFFuT@Fbe+H4cm6rTca=jHXX<>`Zkc8xI{iDzTcHidvJsapq@?-&F_x0br~#6$Q^zd8>m'
    '2JC$Ei4b;{T~*7JD8VMKRLqQ<|Tc*NFN)E0XREXs;9~-+8jbmO6t4kENFK2f?F7t&&eb=iZ=6hjk#XQ}Jgz!V)HEjA2gJ?WP?})WT;n(~tULQc9y}#}6m9o#UZq%;vs$k9XrwR=?G;EPElV)7%f@Jhh9T_KZsOC0}{XZRV`0BQE@O6eSMZZmJf*W(}lg$vpEf8#u6&8dpHRwUVTwfYLPbHgP5T*rGs<UR(_4<(E0zV{L%$u`'
    '656>ctUxzQXz3Sd|_M#b+K76Mf_PYAt%$32Z@(M3`Yu<LDt*pWXGbGn^C*1&{lNnEeXoN<5OV<Mqsp4j>tWW$!mf-$(F}&7Y%0dal4By_ue3fT9KLS{EehOz$maM<cjI*Us*!97Twl9#QA6vFCTE;wIs0TL&9{r)!OwU#YlChzMx0B69>-?a_o%mZDErp{FTYZ4*qUq6g>zh7O~N6aob%zkk0)Tb0p$;-3IgXD#Ms*u^_3Ej'
    'xwCi7+Sr{l;Et^;qN)j{b<$du!*pGaY0=M@`(07hRNB>`pj^0zy;W?bTQFu@?Wg2))~xSXDWia0cU^Zn)Unpi@9-<jJ906msUMJXhX@wet3DZGZm8cb;AsTz4kuq}!Pkjo_&+QxoW3UMz^Df!N+o0yBLW*((<6Sp2@aAzL|xx}=R88jY?6d>)8QKj%h#_(=TSfR&#PTH!4LX#?L8nsVoIsS$gU-sq;*-E(tov=dOSh9Kujm@'
    'EUU!#vlT_R8K9h0Yfy_>0EuTMZol=Q48a_RJiuGtodAC*#-}W4y>!Ojh$e{ah4$?J2^}jS-_Lml{!kHRYJxtO9_=j0Z!UEgBoA3dU*jj?W_i>xvAEv(b_4#rmi~BF9iNb{C8Ed=zU7FA7?ukxmi;dGj)zZ$p`->S$pwuzm;W+`X{OKi%m@Wx7=JSo}ZL2QuBw#?7yI_#S;jT&(~m`{#fn(0h47w6(5X;6cg!5ExUOj>p{26!'
    '+F_FLsx?86|TV-n}sgCPNUiF=2j!M=oC3(tn`dtr3*PLrG#nKp)Z=`>1Q)i#W?pu7)SO0qpac@fbe<G`nS)73e^X|6}Qa@kRJxwSS$}pMY=xez(daQhVj%yy?_@9LU$|h=2jKYQYEc2@H4*;Ab{$1s^r&sJLuVbak+>@J5{&gW~j<H--3Ho{^bykoxA-YiE1ywxTjlNu>U&xkSVDMEPbXHB$O^!%{?|*8h#&7ptCCflZ+e)E'
    '7mck*xq>WxmoRy~mstqy9u#7lPv^&5h=vR2=o3)Qxr4OiHlq-tXc8Z1DWon$$RT1!vhh&%W??dqLw*Ed7~JkYKMlq|f^3SSV|CBHyQun!IkPUK@y7^B||d7~ZY`86~8>hzwGw+xn6xc2JOH?BmRGy7oo%nU(dA{o-md+6-GH*L_;HOV@--WFs76&eAg}I6U_XCb(hhP)BnS3%58t0GrV70G@%Gre;IMnN3IzGW!q>NFLmMH-'
    'i*e88jwiDj~O@xk`nlF{jpFkdiRkfw3ZtZ0(_H&XwK3)Mw&XYEZ7lw}pnQSu~UjxeSZ((i#IMECDd;#arB-uV#9}>d^TCLm322kHYfthZyZA#oo#M9OCNES4F60QfDkE^-t_P>w!OH&L6p`dF|<^N48i<I}4tA@dxi=FlY*%X}KKa#G}r+p0#DWc#Ah4esc>==JGZ`zq3rHLU!P=$%!x$Ujc#+lu?OwPU2-CY+md|Q42<W!m'
    'L%v?WsA{qb-Iciyh=<>0ea7Rc2=~k&0&9QhLg<vDg`9AN+;3n6t>rl%2iX>JFkZ>GZ1usDB|8pL{IYX#ydR-U+=g5@l$$$5^AG9~eqQ2u*;NO&2RwhloXRG_YiqCm@*mb?zdP>oV)_4*Ir1m$zza5D@$3RCZ(O9FIOI2{vh-r!+HPRib`+sp{-WfR;L^GHpNK)Og+_^#G4rxi&|4!zh^7bdH}uN&&{t?;ZhxbRK{Ap4L)oLu'
    'KIl3%9Ri8r?^qubQhPpX8ZNeQOMrjnTR~CIwQHgETU|_tLY;Kk{8DqKhWzcl9(9`D3skp(<Jo6VxS)qGaSY(xIY=rwe9&4ME6#$ltHpoL?jCCMT1bg|oA>oeWt33ci~(cDdlg;fqE%0~ULf)48S>%Ms{O=+B;4GNY~>3q1pRZK2K5`$rc$YM60apmEW!$LWE0xsKV&1NxAa+%N$tL3>d~$br}$tw}I({S*CHHvQVEf7C5Cl2'
    'g^L!A#ir<FpVwP6gkraDpRcd7+<+-(<&$Ff}^R!nRdx(Hlk%w$*-YK(hQYeyN62U&eyxR6#7yu+~n61R;E10uys&<+g%M)0j9<Q^CO)G&kjykW|A@ovGl!Fj=70>41yTKWjq{nA3T`Xi<T^Fw)k0sWd!{9Tn)0z(QkyyGFtGjK22y#=L!ImalLm3Jw0CHxU-0p|@ozFv<L_!gkP@j8693PKD-a`%cE9d9@X`WbuFUGZj3AWd'
    'ZvxMcp|lR_=Eex_WzePm?^a{`^?b{hz$(PgdVNv4={YaOz|>kz=PkN!&-t$00P=emckQXvmqYbK>zrZ6IRh1h+N({pP+Nx%*edH#EH#s3tw+R0{+Cni8OG7v)NDY8ifEw)h71&%$F*YH!s{+Q$<uNq^@X<X!=4^Xne1xL7_`)PYmZJ>Smm@h>9mg%|+g)FRAq<<M@S@sclwXvE>5@(Th`^aU5r1CGQy6}!pf7$Z-HV}<t=V~'
    '#F<f5B0P6sz6sniI$XgvppGkp;*hQ-Rs8k`BwF-YFc^w~BY>UH%34hPE$99VvD`wpx}|Zc$R0P#VM0k!|l`Q$)VmK4mjg7;A44U!s+W>w93bHNN^LdY7lb?OS*=4PYwob3W6{42R`}YCMZn6?8C-J!{;G^GTF3EIJQoja_hfu1h$I_@J)3gwJhu3HIK=5hltTKndemw+(RuW@)HpMJ*0``jcZel6aiIONsa^PqWZ=D1&eYrb'
    'xhPL*R5i2Y-aKu08P4Gh7P&Wdz|=*eOM=KH<|y3mID^f6gJUq6i(7A14EoPWnj60H0+XR&)#HUAEmnyv0?z2>}Dr()|0W<K>%J7s(Wr;A1lhhevX5j!b?=g+Lp@+-KmnqKhyQ)jVhpVRkX;WoA4s&9uD{)mX{$^4RV%Km_+ahz;+=*`-Mu_6V+4L^A=F#PT>%O}^rigy&!Q$H#ZOxK4fz&0^YAYb~^4$gM2&(v%js^9n#m%e'
    'pfGpUnry8)UL$kgFbOBe=_}DsY6e8U$@xsETKW1`E&#KJO4Vr+h`ds|yrxKIDWPvZ(@x=a_|#QpV^R5l}JJWUzGz#_%fZ4a=%nFjujoe$_};wrV90ebW}``?(p}Jlf&sy)828{USnU_cwe)vySdJSC(Yhx`=`FC5|Ah)qdex=`jGjrBjD7o+&)WaU`<XBCTg~c1ujhB2kdTP?kUjuU}%u#gN=#+4(MR{6g&%${1%PX=W0J+('
    '$0~ucT~CT0}~W*dD1%8%3&Q{v<X4{Ctd#C60V)28c3KMx=hDrI!Bh3CvYDOBan%Wl;sU9#1VJ(Y2+Ol)Qe39Ijo5(d=f{MHq|11(<VGK@K5F(4LWNxk$S8qQtA6HbATy0imB1L`U!0*`RLC4{xaURU+#+uw(>rQnXR=9tY&|l$wfZs&!u6xOm&BwU^&UwfQUf(d1W2Dq)+#f6wr-e|$PSyOIZknXOi3Pu37~;G-ZrWpf&3nR'
    '*Q;kcY|Qps$}|*XK}D!rS?kG{B~n5rZ@|>CfW`a%RmOCd#6^ht&)FRq*z~<vNbb&0p#P*7qh4cK$`eas-<<aCk_ed<+)I(%39}g3Q}D8qkb<Mn&_k3BYR?A+`=NG^e0}Q0&t>)$HO(E5+TBVV%69C^J02SrmTW5>2*vP_9Z_2T|Cn$Qau7*<INzU>|=-VBT5@w66kYpLcSTzR{Ib3u+yNcZmunU)j?~kHt)M-vgKvlnmIr*u'
    ';E<vjekp6y8u2SQ;TXpdyI75o=dcSQP91P~m9z21bfbIW3&m4hv`rWbp2m&o^S-13_JF9z~`@1&L&KVqOgyvtd=e$++1IRGnkOzGN;jzRm+<EaH%bps!4e3f34?E09ly-}Q+OcNRg9>v^^{vY_xsM-AP;gBAuB|DA~sLkh~fSJ$n;=OGpRS$VZ-e1={?s3`7-A)@UvvYmDhclDD~_%Ne;wwgO+rj*XQ#l*4zYDJIvzwMYIms'
    'Pnu)AnR0u~Q>8ttiTvAHm7i&eTqf{(S7DQ=cK*n=Y$Tw0f@v-#jZ$O`wjfVav#;hl@-F;jC#E`6sQ<XP_72PHZC8&V(kcSp`Dj?d~#z`p+f7F8V01<OfS4g4KX<ATY5%)$%-?90WO2X~#^cmTWe#4}ADcYI6&PiZ@d2(<rPX{BQDym(o0bOMwn8C3trd;TI`<za2K}vaR;iktSNb4RDfB8Q2<mLl=D%ptcSzt)ACo6sL;JIF'
    '$AVpVVML{X>i(++HT(oA(~YXIL^|^1B$G!BC-oHw4w1s6=BqyR=oyB`_IzTzqiB9bLua7<FEEo!}?6Pz^vY(9b~q4J<{t*7fD}mktn`kdeG~LL1r-Y=m-8=c1`UbP0;sW^!O=ck>sIq!ESZ@p72YyrDF5@hyw4*uZqj>BQS3pW7#I^sYXBYiBsoyc<rvaALElVP-SPZ8X;*iuxYIG$u<bq?-+(Z09xtjcMCwzRX25LSG1bFk'
    'V5h3q<JA-`XOc!cNUzZ0LUw_p{8>Pr2=Lt%P&<QnR~rd^wrOo^C^LDlEhwWExJK6uj64+gT*AIxPy+M;=dsz%SY-f=k!hWc@XWu@jk=iXVAJP902NOvMy$L>^cdML%>dh{syfJh4QL%KXPR@R)omMWm+HCVTlqtc60PRhJay!>75|1EqH;p6VY)!`F%Pj)o&^;8fH)nJ0EH)haS9E<_p^hs3nWfH)~XU#Gnl8BjxVfM$~O&+'
    'RSU{WaD5d;V3`H()(J?;y%HTzT)EpCN(QOC;3S!yP)(vRW*mGx*2YHJ=ioKeTlZwka>e9Z>aAo%Pw9h})<e4DZ0(!t3q{GT=gP;{YQ3MIFTY{xlp|M5z1#n>#X9MT@g<v{Edt9I#~cd}Av@=5ICUPw3Ji*GhSG6`-Qhq~B>e8}mS5)qsY3vfO6Mht}0a7QLq{Tt0Z0NI&4?ru^GGHp<obmNk+tJNbS$v{H%OLPEZYZ}IrhsF'
    'a%GD>;Dy!HPGr8}`;L=F*+;7P)#6dSNxwqa6o-B6Pa1g%fsGdlaKu^k!Tw!38B+Dj@R}Jq78$(B&cHS9?c{JnGQdmn}if2DS1seqr@EqjLZ9jwKR4ljopka0R~o7(Qm{Mwu2f-k=QYcA2?f3UwUv62jgun4M~v(#m_-aOZHPg495(_U;IdQCZ_p9?4wIf&|)7M=&Pe%IqVpQknk3EECXXXAk8}xSp;LrO<aRC?m3M<uK&r=g'
    '-C&{yY<5^kFkcL_Lu5mQxwr0xRJ>VGY|@<@U99IsjtOFVBZbR;Lz<Eedy})a0dEl+J(5CLkhHSx2-i!oE6*CVLIWQp-xvl-mAD@me{lavzGO%rli(%#IsDM*w;LhauYj>^f7#A%`}N$;4dq1B`p45@bO4`QyAz-R0Oy5vZdvFsqvE6a6w)#pfRvE44u*6{ky_QRYV5dJSmibf-KyQ_E(DBD(;m#3&>IpWq&=x&XPz$CppS4y'
    ';09^Xd5j+DPz;DTf5+h<w|&o*iv1_j{6Sm<#wlPUB_;*N3(;|MkmGf&>1l8?TSSCUFfWYIL8>%>RG{O&N0kb<;26HML5nmshavi14>qWmC{1Eb~@xMhltZFwO(1cxR?R6w?^g1q<`|nnc~u&H=5kwc!&H42JvwaRTvEr78rJuged0lB%Ea!O@2;fMS%Nseb`QpM_~#pZ@)qZiBSP#3xm(ZaVV4o)4sfSe=wSApPM$HwCs;nq'
    'sEw8fxb&(1^|ox|=tDf<i}vJ)UbPjMMs4Nm8QP<``kBcu{P=(g3Cq{fBUfsAeC(9j_SiL>%S}TJg}Jb}wb5r8nZknJ+=fC=Xa&QN0m8%@M4w{SM4HP<2O&M7f4<lUN7m$!!^|QkjA61MTOu_@w-4ZowJ{Up_KGlB5D-F&FTAGD4ET0RCV0XL;o=4X+hdZ}qJP);aFkkom!W#`oYh|NcWNT^pPBjvM#PfG3iq==?8+6LH{r@1'
    'cok74c$v_tC(8s+JyZxveq7c$}Cx!M^@mX_#?ym`eB4Fc$z1_t2q6c5DWJgZeC9oy>dcni=SJ@Z{8b0N1B6PfQ)oyEs7_Ot@<W=>lJ8&wfRb80F7lMP(YI*_XaCxv2_$do1lAj!|`t`?6{-)@yku2$KOpH2)uTcfeB|dUzG$y2nJuNqum4{v5(19|i?LnqP~IlJb|LeIz_eY1Xlf%4@fmotNofzcaQDvx?2wbX;duQBDEX8R'
    'lz`9_)`yC_38<c&8Re7e@t=x)aD;D!9|$w^d8#FBx9fJZ60O$}MLeJ%|{QQ~9g`X)TbekQdbe$;*gsgI>q}00EGy0*H(c%nbsjvBYQl0ssI200dcD'
)
