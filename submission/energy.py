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
    ck(U5<=F(73670242301,2**39) and gain>0 and bonus>0 and 0<cost<=F(123327920455,10**12),'active energy bound')
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
    'B4MdC#6aq)mscKCGYx10rnnhRzn;V<a4+o|q@1u0XIaFh(t`Gr3PONd>&I8RCX>Em*vTwVVzK7wZeyWOu%%hLA5MLNqnTVTrHQER'
    'KS?>?sUC_&l2>!Qea7nELOL*3l--%Pd+VAcK^Sx4t2>sp@hA-a_o<5ep-IC4CW{qXaN~3=EKIP>^e)``N%T)F*iEKwJ+=7Px+NJq'
    '<(ZXopj8Qju>-sstAyqt6BT*Uxk4~qjq3K2+7HhTN1er+uEc%O2gU0Zt-Dy_TD$8wv_s$Y&bEO^^rlevs-K@E2LBT*kcA+^HKKxQ'
    '^&|^}SlkHinDOoMEHht!wBXHRI!@ix4JXab+gN&Aler|AtYv}MfRpZwf(bej7eM|b{hW)8bhaG~L!qN&YrU`*AjgvH6RrEljhVSB'
    'Y%?*;80hR_HA^BDmw=9yn(_)%s~1;1B)Ou}6n?_dbDRK@6|8pnDxd_{Qb<3KZU;^qpeo$($l@yT59^FW*w%Dn(XC`kx@H=0L4(e^'
    '5t<VzMH(pd_$+EJ%M%*$xS!Qa4-yt>4sjElrZ!3OAY-e!LxM`lHB*?^O>E^LI=QC(Nk3GX84j6(20=hqes5X2JAg<nSp>_={&!LO'
    '+cjLs&asG-AI552AnrNxl!I(O0AKRFQx>0(ev|02N~LRqtSXjb=8XCbQN<*dosVQbO^sK2o_?33MPr_cy{?}C#57rTpaIemlH=XB'
    'ytg`!MMz{+xc3Xq4{{8JBPfjOAnUXX`Z4eRiadX9k5W@G+e+XR4LK7u5`45G%G7FD=iKpZEez?0D&Lx6=9%7J>3u*f0e!<xwA1KG'
    'K&c|KApn^&`o9f8^BfKY|8Sv(*6YTMi$_v<@MZ=@LHAZ(?err_%-RuWnV;~9HP>Km7u^IZkLJI!GI0S8Zj`rf3+y8y3^+WhwfK{q'
    'm)V?G6^F>bgxQLOrQ9I;eC}P|^G~M!QqxI1W%&U3@Z;Wi@lJ};BN=aUWV<C261^LvqJL(ubIdR=#%PF%ij*e$Mkyta7FbD8(Uspt'
    '%F4T+#jjKjDh0}$+6NRf9Bdl}PX^8kKJK;$L_&1F<lO6;7C6kKYb{jl@5_NHz%dxB$R&o)Xq8W3C;11uw=7HsF+zdZxOLh_v|H3c'
    '<qeANlf+xQYVrPW9*vCl6!2UU&b)tfFYSg1z2Zbwd#2onmjf&$11a)&??H4FTrb?MDG2LbeAq}OV1$&A;2$iM{f!igMxuhQ^`kSO'
    'z>UhJ5DNGs$Hx|F-cY+e_m)h=3s4Xv5Plh^tKeWrR3wWZTm9<vt%^7-sm<w5eAz%Y{8ghHhc}|HwSQr#K0~G8O0SfkQeYHsUQsgm'
    'Ufo>Su^JwZp6ACiz}p+tAv$Mp*jYG!4q@lBJVqOxu88AbQ&w<hhZN4(W}~J@y?M1dEK=HMuLH&?HMGzfIfw}*bX`(6)Z+_JPfg;G'
    '6KWoc_!H?sNg|6%nFaXfXdL2J7(W@JD;Ee)ku$uWptZ(?+zw2y&F=VFROp^3Iw|&?CnZMB+t)&Mv5I`xWg-lDamWu7V<Y{r?}?nt'
    '{0kB!z#)LTT#3RR(F}<<A;W>)u-4TkoNq6EO`rYX$bk&9g8v4>QHJY1f@oeF6T%ZBc)KOglB}UR#IZQuUlbb#yt(j+2@o4Os#FL)'
    '0ID>yI_rryrKuxdyPv-Nb}6w?!~pv4n;H&_2+_%6{<XF~v<)`RmlXvU#Op?Vw4?AB$jvM)wl)uI$G@<|_r8)N(Q|o{1okZKQ(O)a'
    'Pbfib#aX6p+OqZ+9Z?VdNs)w8>-ZTnF>7jQ31h7t5j8NqjzqvOks-m2qn&=X_E3FJ>UjxGbv_#ly|!MC{{t|6X&t@b`^ey3nKI|>'
    'lEYx^8S93fay#t{e%OjpjJ}Ttobu>O9S8;*Ao&)hlObY!6~$tKweSw0DhklUOOkEr=L&#4o%guu@we7`*(o2M=n&-eW5V<5y;p~R'
    'wS~J0m_vERKais3G7{{5;9&nmg3W)fB;~SCUBb3C%1~I1v~EE&z36&DXojYM8LDroq0Nzq1CEy?s<d?8)o(Uv)wHN?gw6k;I8rFc'
    '0Ttb}b3r((9>+eiVXIo_d-kkQCW>T+QdkUGBUN*J=T0FixrrqBjj}ghIzJ)GaB#&v+_^x~q@8zF{<&ykf}u=8Got$8E~=W;=o!87'
    'xU`5<=pp&=enEU!^$3KvrJ$pKnhPWj<XU@1$~!^BQ_DN|u`?Uxv0GR6S84AcMFUiahX?2N+$%BQVAF3&^>wf0t#{V*f-H)TE62)D'
    '7%tS_M(!D)YoV8lQh_&-%-+)M5OqfHGk__x>+om`Ao3n@uMC}2$s~qU^`;dtH(XyKil1a^f&SW=i-CK1XHz1ww<Nv8^kG>-jb%>1'
    'V+vX)o~^l@#ivF1&+B?hX8HxD;JOj<M0hp9l#0_(ElL8jD3Nrk(V+tF;cN>3lBmADx4qP*tQujxqX9k9{(-HrVjbO*EAHvmAYRgz'
    'DY3cT;stt1hWhwv-5nKJp=L{?5DFe0xuZqX7RAy8G_B{J3KG?t9ytiIm|=g1Eke|t#AXjMm5daUl2`)!Q$C@{P;99&%^L;2LC;W&'
    'A=@btO@pnV?!e=yIMnFcv;6^W+~+Gy5-C8L(`-gE&_2ubm2X$cC+1Su)CD#0A5|HuL1yx_{V9dXA?t`siA%b6KPA)-Ak*D=h4Kap'
    '*L>`w=)sO^cA9XCI`P2+3z`97yV=3v!oX9I2YDPF@%7z72`oJrLkg*+O~t`=#18>g)5N&W7opD7-Ob+SZ}doXj52Y*t@3gey+>0j'
    '74ryeTub&K9-wr6IJyEdcWQ#8v47KwDs?wQRxNRtC6zP5lVkkjw)Y#|d<hJf3@@jzv($rHJ$@yM6ZV90JQPrG<YS}6>wmk6b<&aG'
    'L#4k|z!rhvTPxT8c1Z3!!-FS7nv9e<>&+~`LW+M$bICOguO597)`-bi7X_1}hh<U!aAW6$vQ|wQK1Zi5*X^(HN4UAO+@&1pNq^Y?'
    'hqiEBec9h}^?M9EUbn7!3)WI`>N27tKrh4v?nu{)u)ane3~Hk7i34dWA7kU(Gim#JF|UHKJepFv0!>7bIb^(Ixf|9HmXvCS=c8EY'
    'b(qft+f=JrfsYcmd#A`VQF`=ZdK8s)U_t;{$~G7<RL1B6sL`3x*hv0@Pno<Fzi}BN?n5m|p{$O|ySdDFmL8Q6(6<8cK~*coN2u?e'
    '$E?p!BvDH8{}_S>B^k~Syrd0CYOKtJ2{K#ek4Q}@uZQkURIOos8B92A?)a=Cds3v_l&2!a!E+Nh$5W*~%0@#s@5Laa{P&Tb2g)&w'
    'G37CB@n7_f#cA%u0ji7wcvA2Sg`nC3XCiix(4tbBK8ScJ!mK*2#mmCj2=svN#Qp!Ss&J2`KkhJrbpnHJqn2?$m6tHqPsFp5_$7lW'
    ';{<6NKZ9vYJblG5A%_jGG{5%^$NkDXd^o8f4!ZgQwMtdQ4hXX>0$zVOQ?);_PJ<&Cnapfd!{33d{W}Y%Ci=>0>J6Er5W}<Ki9$JM'
    '1Di^{zm;#`OAS8(>@nZV(GM~8=Lh%A5Y0ezOA1`=%(~WH2m9UJW7lrrK?)~T9*y9+qnXT?V|<4Fb&-5a_6@~&YG~=5M1KaXErg6('
    'x$aa=P%T;k*pn5)cPO~|1lzvfII8i=?nRm>yX;;hhEnxHw)OzwzB#u5xQq{(ooVdbZ76_+hL;4+x|@!1Q;ygPYnx~;eepbTpxhhE'
    '1=w)>w7uH;4(PD*oIH#Yw8%&=O(F=H*MC*p>wT>AouNaf`CVG#6fL`rC#_v_v9Xb`LS@-{1dDx0u@7^G%4rD3ANhYocBP*-d|w+H'
    'R)G_S$p@cI)g=$~8+l0D*5-zkEtXpL1|nLnIwB{QI|}VhXgvx^cum@>SvnW$5E?_?05BGA*|coYYy;dzh)vfe9}n&IKm^(JpC$6O'
    '15!|S&g>WjR4ff~og|?Js=Wdb>?&4dsjoKz(^4)V!$Y}0`c4O~23F2}-(+LHgc;Ie&rEEF@9vqJJ;e6&&+wzlD#w-VZGbH<SRvSh'
    'ofiHF%3Xo3uE&`?wuJ?VLH|-!XgVx40c^$8+05;q-za^cG@~93jXzxPe;pjDpv|IgaavQJHSivOPAtAGMnXTWSli+2KK(n(Oab*Q'
    'kX#1tPI>)-kFt|Wlm<mWRa<HYpiN0b>R@Um3^zBk)o%;KS!}W4o12n31^JvOpx;~7NnXY;4|M8cZ<POK<AiSgEgpFpGHHN|gJv1;'
    '=87;#4zrsos__kd!<R`H(D-#c?p~aL49`~Q{7BZsjq?3RYTtAl{wJU!k;3vRBc+bkdj8(hGk?K7l0SBM01nxo8{B<)YClkpGXis0'
    't!>7J49?DdS}(zaiLg&{o_(CgPMb@TfuToIA|U6ioUfg=u$~j%>Vwtnn={G2%<+(_m5i-hBpQ8Yj}UUWVXK7$*_-3Z;V<6P2|6PO'
    'e~qwa#=5G!x>)T1$70=Ku_57QuCo+6Y4<>7EXJ6Ic2aKa#HjfG{CxqqUe|atWBNZBbE?UNogGl(R*m)%u{5Cj982gy{j-b6x`hXd'
    'CCnyq?ceaqAM;7vhK!-lk79r(Y%S&rpU`cyrN2tuEPl1JlE+U2&wpl{-@bLxq`fb4{<{RaBdZDeQ$IAfiFBEVkpK;Jd+oqob)kK0'
    'vk*cCzD*is5seWg8Q9rk_Bj+*onC_2p!d*`?U1#+4_Bf%lj?K?6m-s3t<2Htz{>|II;SQc{~H#W8>}aS)O%s{!?23l1w4JI$-TNl'
    '4)hJxdNlYx>Iwn0Y!SItd&~#=t<WA&4BZQ^f_%fI1xy*40?#U93~UrEx3=RS5z6XBX?6D^5h(I^$?7`JXZUtu$A0V%9~7<=oc}2~'
    '06))Q04pSw3UO_%GM)mTJUkg4=#~-<$Uxbh;epgMlG~TMg!42WOpo_B$988FcXc!jgrFHPLvojfT`omj;!4LXQT%f}L*z}vj{AJj'
    'AxTrJ7~ayD?^rc6BRr_PJ@vg2YZd#y6l+}Y1v@LtybkGqZo30A&wGsp4G$Fq=5;8p1XTL}b^)?O<?gzJbD1+#FMDgwIFA%1=N2ed'
    'GVCV1RLFn>m5DWfwy7;dzB<%xQV)qMN0hg2!d2sPbwdpGpj*Y9NrH5jd;Jf6#iZKXGO^yT|IgB|`-QSlnm?kO$RtWEprH`j4HIz|'
    'r+X)VxlOs5ZC?X0OEoYjE*EKActe5N>o+wgtJY9NVx?cPh&fgl(bBFbjM;UBH5-GAG<iX6)DM%>D}CCFuR}!4M;LLF%(yUCxSX}p'
    'TZ^CMuhU*}!hqYZ7#W~;%i^Ck#E~F}?u3GzMQC-W;9St3&;VH#k^Hp*pf?_s!=C~OWVhnebaPWCqxxQ|Y1u>Se(jvWs(+hi?PpMc'
    'Z#%v%uK}mU@~$cZVt#=_X~Bc-kdRiFD%j~tobYI$j|4tl;NE-Z@9PsEk)J00QnD%a78cec{q3Uq+6E|O(@GM<D~8&;n9uM{k??=F'
    'xWgDyqr=>_>&{Vk=ItWP8nbMaA5w4AWk&eJ;!(BIG)RA*yyKCmYJwA^3|@eqghr~_0!9~NIx3Hu1)V@Yj96>Pp#a@CRl=2a0y;V|'
    'Z$aFJ9={=!1_O8*ZsVbSLgA*0+0t^7=Z*@B=aGTmm14V{R{pm+W{;v)7kwUa+W~|M+`MLzL@_gp4mvAXKO)4@UCdCE1Idv-or<_Y'
    'pq><Hln+zLZv@fFShYn|<_%2}96%P;Hvx%;o9^xAE~?82H%_12H3yI7i;Q*5OmhySoQL4)9dIAKVzn5capVc*w@}?&J7n^sltaD~'
    'g(x-Jj67QH-PTJo@XWPaDJ3zh?V2l{*5BwPVnxoj-u6u89kfftWtSSkZw4|OlN7Mlhp-XO^FZ5CGrW2DTOXVb;`dII8hBY-Z<7uA'
    'DY5Mvd<R|nssF3y-tXhrYz>~oR$c2Fw8eca(MO1r9kGAk`;)+a*zNoKI`07Q-f$a`+wLOZM_IkS+S7|p5V_^!n(JSQ!;f+LdW>%C'
    '+q=&BBauz@4G7MPMTDO3-JC(;Y%VD8Q~{S!;np_j)m!ez&UJrhK@MR+4t(KTLqwKR64B>x$>ODCtjpQm9o_Ph^yPBHpaDxap=t)j'
    'B0lk(f1ePnlieVNlOfxi7(%xoq#E=6j(LPVJ#UvMt*NtLN&T{ZSR2%UeO|)>xQ@)8xS4Wj-}jRrvgOy^0$8OETIApi+k<$oo>A$1'
    'by<g;cHj8Sgi*=G>~I8>7pi1NT~n%veP=yXq2}?t@q_#V1FK^>=;mb40o`+zbb7S=U@V3iO~lUW2{{?duX<PHn|0~><plJCz}8!E'
    'yaatu5bO#3J$P#tK%$?JhNqaEEjS>a*W%gPOaL~NBmRfrbca2BD21|HSbQ<jr)ZB3@gqgY&ytTGmmxtE0R{}>p{<OAwqzGPM2ZuD'
    'HmYDtyP0y|V6XY3nT@u4H$({*&`xAec6$Jqeqe_zmj8nhn8U=c?o9VThi)+CrZr3jXE4<?z|vL#M~NRGs~og|bW*uv&<Lh?<>^Ae'
    'Gwnw=eAxjP79Cdla6gg?cXEH5gAVzli#lWN2e}ei5L{r+8yCpDFM<S?uK+o2Ml(2AgC>0=Fw96Q=&*L00CvkaxG%M2`Vaz$Uhs&d'
    'nZHtXvk2>;>TIo0O?0+WX&RZ>gG-jmEG?suAEj%&6?7+OU!qU#QcJ<ePzWUp<*^9pkjXniBZ~=l88{FnVk!izT+kR`xY-AH;PeM0'
    'lQ<y~?a9Gp)y55l2Ft{JkjbcbCZ>x4ujr{!t8qq$u)5NCW|yt_Z714M5sv`dbOll8uTJJdq{C7z(mhTY@juOk!q!sXE<D(vqNK;f'
    '0pAgAe8S70h?>%SLXU#>RD*t6=oMP0=yjR{+K%zT`rM>7d0Q5*mW^4t8B@ZY8x~y8n5DT2$l*kxAO<^qujvIoCkPgHed|G32|WPy'
    'QO<b5OJj5vQd(vql<4O0esbQ9st6s0v@AviKd{-a<JPF<5*jwCyZNawMlZtqjY@X1b-1~kV#ykm$?5mnz_O(ewt1%H8B`HNcYHG+'
    'M=csO(W~zCec78gPY_k_3%DC`c=pK)^O9Ur&+Sq2d3fc#-_B3WNF&3fTkd)atC&7k){lYK6yx<kh`Y!3olsCNkw7beZ7g%~YR&|o'
    '+%G`<{Av~pJ1JyIg`fie4)Z$Evluxvk5J6%xf`XJh~u>)t+u@Cc4S?LhHhtGA)4W6nhAYm>mLwT%e_IYcJ4@Rf0t{0UC1Nu0?=!|'
    '4QiQfV=R*%2Ajyq2wkX!XxSU}YPb1+rJ~0Xj)(l6ZB0Qh<di)Ft`)>+NhS5PBR7k;vKSryKO-s8&LpCLKF+PDQMmCZ-tf&YGGYQF'
    'uNve&{$WxJf}WaJ2BTFKI%jKoofAt@dq=Ln9&1~s@!TVz6VSwZ;^1REg5)<fSuxkaA}#nPBlUho2oYYHE;w?a!qx!VtVDJ*KaReI'
    '+DH(*JHH>(bxItjo%-oKM+pe^UGE8sT!iNY2tO>b42s$|3!|4>E6%FRFUs^gHsM`lehOrd>h#q3=V>ZyH?6OzV_8Y`W_M!Xw^_$+'
    'g8CfUq%{sWqw}w7SdOse@3Ja{=cQ5Qzs%Dv2Xjt@SOHZzyHPWNpMckKL8kRI+LjtCwAf}ny$0(+wkQ5zzFV8h6s#nBKmFk|`Q<5r'
    'kS326eL=l1+r{Z8a-)d;G!opw7}Apre#@#N61td+V(hWcrh%9-VKz<xg;3$db>KlH)W$ThJX-l%A~BtdXpXk57`cl|_cWGD?@XLO'
    'fgkX|DgDPO#QT+ANKwS@9MU*`U8Ll-r>KjS$1^&&kI(2J@+~4Ny!{6}zzBhWH8b1|8;$E(O`EWcUJ=SA8pA^_$V5HM1DyIc8bNPs'
    'tU{x5nUC$+$zc``MMdau;IkId?=s7$Xr=KBC!1ejAt+tXlc65xwfYEC5OS*+doWYfb%%VGfgQ*HCG((5g$-F|Z^|LpF7DvX>ZK&n'
    '22I}xHON^Vf{A2!^dYw$jh)qHh#2P+52Of%GCQ;@p77{9fuwicp!6Vm1&QFm4(ZpRYf&X(JhAsKrc$vo@9&*UM75XU6|*Yu7>W>c'
    '6csO$Jw#D4ih6!D!IjBeV8L9a719LeE$~Kcw=H^qmm^kmm-YI+o$@xwT7NPIvK`A(GMQ>6l;1mb&k0_?8$31p9NW-YRl1U1;+C_S'
    'G+HSsA5cQHb~N#hPS&<OO1D*^NgS8LeFKpS#VxM6taHmh9VXjMMh4Cl@ZG>=k`P1a{YZbj4$q;VsOVJf)-YD_nZSZ^adj1D`=NU_'
    'N9l?x!U{ge9|HpjBw5eH(l<WZkL$U?5}}<8NC&^2Qy@cVQtZx$@zBBQd-)ZMzByq&;hL6P;pN~_A$5~eqGGaG6@{|C!{9B5lPYHw'
    'PRgC0vCq%6CW%7kKF|n;)IhPI8Kg~;bz1Tu8eM=oy;A2J;P&RLiD(Y!J}TF0jvh*p?n!I^u9l@@nnSR~^FxyOFpyE8wfT_{aB{=r'
    '5t}9{S~gDIAK5c)x_o-3dnrc0QJhG!!xRZ85-9zC5|1Bv&DHd;QR)g&^_?fnoWI`upaC7e$mt32WI(A{;Tpw#qn!KMMD^Eb?xjPB'
    '2X%E>#L4=6*<~vI1=t(aANPr!<&&?C-xW%$(lFlXN*g2GXi7UIzlZdvjNCclaJOg)%MI9$OhJ6!_cV1LVRva~cM)tsQ!Qb?@u2I!'
    '3D!~_&id+??v|)ELmK6#1jCU2JuC|p{y2Rtk)yW2JFOa$2`2DE*2Uyy05-b5g=+j(MuXEc$7G;<T-89C)k<@T-)3<PnuJK*(qw%7'
    '#CS)y`jvkPM6w2oOdOp8-rgzB!qCyM#N&<Rp1N<5MUr&CRr&{&b#94B8KN4_SsI?;rl1qY4D*5Ew=w7~n2(JN33Fl*GG7u3AUc*5'
    '89zAXFb8z!%<yCJGidXx^g}~&JFcA)rZu+E6R!m4tvwVg38V~fafNF!=rwR;FGvDIQ7)@Iu6ts^)1pxfef6=x!6&wJ#iO(2W4#?V'
    'MYN@9kjYjSgRuo_)fNYACT6}mN=nX3JYUId$++V1{(<s7a93dBtdY*AoC4{$>({f@ROO8v|L_cEcrAE*EOo~|RcPdFT;w3#{g1B#'
    'q}pw@b2Osfg;P6ic~fn*60=PkvgFxPpk@@4LwR08WRO>tpL*0ecBeFonh7^5;zv+8x<P;VWOWWLi@E+zS7pV$T%+*oU`2l7?(kh)'
    'lPft>?nJxe7Xbu}@m!W4KR7B`F|&Ck&^U-xc<&TWid!W;0q&d1w)AI=LsyCLkzVk&hDL-Z7P;ZRn4BPR{ZWoGafe38?;512k>ZG~'
    'Kb{?&BhmHZX4Int+;QUo;Y1p$zo|e!EPmyJJ0Hwc2JBKW;R$W5y%uR*KhyV<@T2;dn1&ER>E4%P7N#zX>z~n?{~5HBRE+v$G}o<1'
    'gfI`1MUyw7Ar<loO-CP`L5cO%39wVC?jwt4l*~guss$PJH_l<tBx9JlQBII2xMf33=V@wt3@TExjAYLQSymv}z+870+}&d0__A5E'
    '>Y9~lwn?gc<I4e0zomrKu7t8_-c0(xx4Rdt63`nYK!*?SUvm36xs<__z0_=7(#u9vhI11EC|=rVdUVtZMsyLvl6?rB_fpcA`bL@_'
    'KfP`DQv>UpQD7jp0A*kj-NnL=BQW%8j&U^f(!~XlG2tsx4|-DgZmapwds-8S%+?40ivHJsN1K*|N;u_15RT5TmK*BlF2k#4F_R#d'
    'cDfkUI&TJ;(jAiC0=e}~RB^o^sJcJHF37p@cT}AEc%bwxhA;lEpD=q?j5o8B>31^uzL4KgmD(tjZY7si&;bK!Y>+NyUn#sP!+<SL'
    '02mEwWecx3a&nOpq>%kxsIr_F;}sN%Xe?D1Uo1eJjbIQ}R&Hnaotpx%>Lrct9xx8~o7!NzAS~nHn_5AleDc-pfVxSiEb<E4#8G0@'
    'hxP3)Qln)MvcHeEnAz^y2Vu7j*<_F8bp=l$rmGzFIE7K$mTB`wu1q6I?UYyVzIY$Q^<8^240&l3=a<f(pdjn-#iT2h;gaeyso<Ka'
    'Lbw)TS8<*+^I8sB<yB3t$+$2Mf}x1ayuTUM*yu}(i_iZKqI`!3zU-U4^T&Y|_G5=P*6FeX81IM`V!h9G-jC^c<;-Qok8-<q4)a;A'
    'C$`{5zK{Z$<91$C`QBt(3lt;@eaTwA3|74U+VzMgTNx{vpce-ZY?52dVqTHbxr7pZEOq2u3YavV>vA(~w<qpI=>*XT+S|3{^WP}O'
    'kV84z<CoLvm9z@yA}^_Uh98jDkS%_tuT_7E49viVb?#>F8*NYsQ>fZy)EOG>!>DOkxie0z-CriKQY!2VtiD)rZgtJhx*6j;<gO|('
    'B=Q!tw|cbYK9prpl|dV>y&pl)7sS6zohq;F@)dg?JmA8Vz*=ZCk-bf8r39xj>S^o@gj)$a+w6Wt-J5+mWm2r;!kDC0^fUK+??k(0'
    '@K6-2x3q+EHOqr4BtWez9_zmF1&Mvu&F#@5k)r#Oc@mIV?aIs$iY$Zr-1IPU0!=W?kYGf9R-o+LBRG*9cup2%%ij!w`cyaft*RC?'
    'CcI5?J*&5xnPB*H1~4luP;<rrA=n&4yH~Y4Mhf6J6w^3{tnCy1N0)wz{O`v+JOAMXi$+{lL<WQ`GxCpy!E%P0`>}0104Siacv)X;'
    'z(l!p^W5Oy=McRWuRMyw%}1ev1IZ7aI)RyFM<~^YM5gTLGn|Eq>elfyvDs>bOr9ovI1YDG@p=dep(?W#%W<qJbk?M!6VtcC)Es5z'
    'g{JK;X3O3^$ljzf|NAf$qNB%11M9T$JO1I~bnm_QIByfdky6c?SlZkO-E07cBKYQ5pLWZrNy8Q`Iwp@4^_UibW@4_=fYUcMH=!m<'
    'slntUCUfB8{mQ#J{?u>aaJSrXaxYraUggvsX*g(~uK}D%SKdlXp3v`l_;JDUG!+)CepQjIw|q*Q%l_NBt@2sIIEX?s&X9GiFq-<O'
    'M5&E3jXY%;IybG_ITId%6nojx2<%zvmL-KOk1U4L0MniJ_+`qzWjAWUHVu&|?Bsc%udP2|vmewmoC#k})H5*x6<}P~2ILc^>OvIr'
    '{@@+`8;mDGcBh|>XObgL+C~5oy2OAlN>L`$c#Y1Tjs6<NiR?_AP#jLrQojL*5O3zgIiyx~>Il>^+A$Ad={Q=Y9L_h8(W0)wi-2O5'
    'xh}WK*G)cUuBU1j&r91qWGX*r@0IO<#D3zdXYl>pI~6f8UebObYSc+`k-Zo6etc)^CbZ&tJ@3lU@``|4w}9LZ0{bg%)J{p~HFa2P'
    'gTwlv>1O{4NTNfLhzD9UF{ns|*|x{a|9LWp%d-C$tH>Fm*^xjv@w8$?=a}DQRsQx37l_GH^&pExxF1!*!_(hS*WmrkKLt4wcgxHC'
    'OW#d8Lh?#npbE)TH&liKNa7+Yh@35HnoBFQAdWmH9V-8QVz{*NOb#=lMv?!<`NV;Udcd#gH~t98!UkPcvw-py>a1=&9qHL4!+%cS'
    '3W*!D_J$u{fyNh55!!m_?jUPq3F(YnVw~RoDpO(zJq@tuj^6`ou-v-&Ed{!~TWF{JQ%m5^yk>2qID>o>dl0Q*wP)1!9fKX><@E1B'
    '=bozQE%<<SUrcNjw&@uB7&1fhROIZZ=Z4)w0@g>%6{9vN0Kz;0Pg#;OFT|4I2^teZ*&=G!^k@hH*zf`HCfJU>uwYtXh@x|69-oSD'
    'psFRZFb-o#29B@&z+3o3eh<aezALC%X1(yF1~5y~f}cky-K$MxV>dYO9IPE|hnefAiiK2NUy^twt$@Gq@v*sLUeaf<xpQFkyj+m}'
    'k0v*(C$+7Vs+?+7&?<myYpnv@6>Ar8>@eC-1W9s5m=?Ua;Lz2Ps-sfV(9B?&SQc7~*X%97e%GXomeX>^E6W1ZZ@NUd1}XTlNc(=Z'
    '(G^ysLi7vfDu74a-(Uod=NWVijG<46N#+h9HHm1Q>=GBd<_j^?_kIZzSvu&|x02&ryGQW@i+QXG@T_xt`js_dqzKa2)MiDTM9&n-'
    'Sqt6)L|cQE@;QY0{PH5*k5k0d-HC=0!{2d+g0>m>#)5;$zE12b9BJ5U>%GAhT*|A^*nCBw&J$%t_pnkD)z(x|b}rq8Y<4qDYlzbv'
    '6h$t+a6~UXfd*k_+1kZ%^8}Fl@R%oDI6W;0(pgIwQ<tlj-u?(9D17i;y!T{;JgX#M#^;d@@yeE_e|zi}#uh5cfy$jU*rIDxRN@Tw'
    'u&MDVkA_K}Qn@S_u<)k6;GaCzb;)xNY0yZ(XO@oH0VgZbT9NoCo4ZaJB5&x4Q;zxgD8EeLmwg{OtU^x;eU2C3@u>UDw7?Jz^_Xoc'
    'Y7hQN@Ov&aj<F~Kt5sIMk!*ZK+6t47b6Y~TqUVK_wFIK76};K3;fi;FSh@{WBh>jMms!lrMyWOZ8&Lb$V?x|!pLFmhlOZB=to5-N'
    '%{~7hx%dC*%LCvbWyRjR9Ir#*O^xZO@KVG2g?+>r>o=8mmqFNlmL5S5<EwoxoQ&f%v!b81UONq3LC2iYLwkqZdub<jLzFSbiaOY&'
    'O=1G}QL~7Y+>RR9WxzQ4FmZ?e&_$JCFITx!_bqVjAlVvUR537eGA-zLRxcFpA6J7pvZ0;sEYDLp1nTo!tXK##?1=lUX21b@a8QVH'
    '^+TBQ11M9asuq|noikv`Au!-jkfjn%GR|cHU!KN2mvQPM3Z~eF2@&71{&LGL^0|+!BS|RAT6zS6QHR^7dfwcp80c&V?=Vq+mXIV6'
    'T`R{6z!nSE<Cz1*m4uC`+KG^hJuftiakfHGB54)blK~^@&w26K2{dx7YM4LaDuxbp2w~_9u?d2S->_p=LdK}KvNxC(3@i2c=@>#c'
    '|M+zEI{&++{o~j3kBS`?5eq(AQ-yKvv^yKE4vNwA0i7=_iAR+s<&l_h>@`q$O^{HS8{>a$dFA8bBPi5fUnG$ktP$~fRf5{RhaHp8'
    'TDHEG&H=%bt4WF8dFE1&529{AhgluCT*$rvsddo%t&xh>O1;Db?g=9`L`hx3uL{DQAKjt}Q#wFQ665i{&Q;=W75_O!609`ZEcLnl'
    '>MSC1HZ->{v#^t}**Io=geLYR9dV&>ugHZrElJ8yL(lT>IzxMy(5-7@O?XBbG=fookKnv<Sz^K3o~MXVzp~EvtVqL|u0lb@Bq@#<'
    'Aybj733|5XhV2B*YltMNVBrR@?R`5K!-z!RMuQbm7*9y?Jqt7L0`aEOKuVHe3g(d}M$J}=hMaud=XQjWxA_GFDPn)4H#W6z_&mhG'
    'aEi`f0RJ+9A|R|*X#iOe^9KEcH_>h2AM;QU!$l4XqvU_CqFd;+?x++uYyA)=XBU-kUi<RsO^#|`#v;`AB#-%lf+y?~P=A4M4-7<;'
    'uX)flr7G*p$=MB!WV7atc9Xehx|jzDA1#Kgao=J*=EZdDBF!i9bM%_99p|_Fm2jDku+Z?FcYFGPSxCW<ZGOuk(RzT5bpKF7KiUox'
    '0Zd@l+pViEQ~Ig&75eN5z@17u%it4Qal!DH6-$srqM@#9j6W|9Z$t+WXoXWJ%gDa5ZY`hnof*|7o<NIEXse>Tneo05JEbkjZWw?3'
    'RF4h^FbEF8C+;ySn~Qs>@_KWc0b<#6A6h^dkK9aN$-v+qcJMEtk&UELW%(J>>i&e5GA>YQUr>h+eO0*G?0xAJH$J{(Q=L-@jlw0Y'
    'j-52;jthbJQWQ?@qg%wqdCAtHQfue-(lKck0Q=zj{WbvVdzc3DCxkHCTvh1IDbOJz;VSQ<RSvOs1Vf?vw6IhI?VBQHA5++HK3J(6'
    'js+-z0WNn1h+VL0?U*ae5^^zA#atVPbnrU3Gxhwah37VZip^BJEG()k1;Du8`$dY4;)GU_8s#bvzkHph6%!!+MeZ0!G^1V%_yYLB'
    'dJ@q*e{+3@5gTzkhD+i$&f)umVR%&+#i8G+{#64BNR-n$hKr2Bsc)GR_qP*D)`P=^9V&#wkWIxpCo=Pty|C}UXlxIh(ps~8P^80S'
    '&vNAFlFwO)-LLP@m1kijyLu2_A3D2$kTX#!XO%RpB`Si!No%yz6_{56K<-AU=O#R<^DvcxX?*cG8?hwE+kcZd^SK4PnbD;;u9Jg+'
    'lMKj3{ML=VX(6jaROw&ZH)`}o^U!nP-ITaJ$Ic%wFL)WE%whq#$60UKYGbw5X%M2s>2lGrHzqj^pG*(NQZ?^-(D#53s>(t@&gvEv'
    'r5<F?(RR?&UXsOhcB1-54jYy_;jr+@QGAJk-mZu&w=KlMDb7ZAVIEU|B(YSxj}=~`D&0a8R04?2OMEBqe%TIcI7=J8!dvwxJ7WhV'
    'jeaC4>%$wvj={~3KL<b3koA`q?Dr+l95^ydqS0jq2bSq?-xolZoUiE%rs@dUV72{CUGD=&vdMxhf?h00sbY$`$$v}nwON19>w7GI'
    'SbJ=B10GQn+}!)KCWl#WEqz*~mQ$sK*ON*yC8_dD>7B}(Mk#r_<=|p0Ss1ObszMnF5lgL}zXKPW%C)94uB^4xMpRaHJ<w*&Oy97Y'
    'U_j)ciiz@VTku~ORli|hJ8pw=NY6+mmU(*+xg$(<o*W&2aVNf68tHVs8Mx+l9mb`0yZY1*=5e0Iq?BciYdY0vF=IbzuT-}T3z3k7'
    'fRFm^bQWgCb`Ac5dNqYP%aYwL?~AW1&K~}Sd{^{~ZYF-|TyP&~-;1D1a{rFkxJF0SQ&RE3Ac3JKDPeAXM{>)$)OSy=P#GU6(NFv1'
    '7$lu&OBcy@`lpbs-5$1V|Ji^fyg8}3?B`I8D-4#geo%Ob(S4X}Uad9=nPfDIQ4)jK2LLf{lqA~lV0edL#LVP5X?xs0bi?-~lX&NL'
    'rL&j^a~mz%?VzLBfVZE$8Qh|BpStk{2TyQQPfnB7vBM#<QTrHesQVIm{(*Y9VS6vuvT@hOj;vs)fB$1miCS;dG*B9Zc%7>kV$~)m'
    'S+5q5-7%t~W?{~lF?-k#UN$b5S;+1J97I!Tz8Q!y;}|a6yNh5!7-?i>rMTi8+3f>YvKc<5WPC|BtdWBbY0^g*uY&OY(UaWs$hA&R'
    '!;bYT5^%uz(o02*0xgZ>+!<j4BUXzN{xnFh_HMvPC)H7g$T*Z9;=wa>RlZcm+NO!aSEv+?*)M3dgu4bZRZd0vEAP<^ljm{NSn|}p'
    'oLgi;Q2<?U2V&Am^xpVK_;dxh?ye0`>n~-No{UbVj^7f)m8maS;^BvfA$YlOtXuxDnRWT#cx=-UI@ANLWFH(7l2t%C8l0hOe>ZJP'
    '<!GMMnhdsPW5mlQ%Xmpb^dF|Byq{i9zpxnf0unxLoXtlXk)0ft+DHK`q|J@zfBG&c<xvX)9g52EmTVceHkV#gG%g=gch2lMG(09n'
    '=oSQ)k9ts&MY&<i8>n2=7G!w7%FDrck<+R1u3<}uYZPT|uW@b$_owxRhl?f1?a>b`G%uHVOGF2v%)04a(T3)58Kqqz9MCE{rl)XI'
    'B0p^clAaqzs?UpLXCRG^P=cgI3?KjIR|(i+%lO7)Jz|B^tL|VEy;@!DY5PIB@BgBGx3CuMerT=@OINgjbSJC3(gopoHnMKNdf}&B'
    'qXxWQRi1+-p55nf?Av*H=k6-{-N0AkYvQ5d(*amju0_bT`2MS}<&(+G?&<|}jpJ7cm37!f$+O4!;v};Zldvq{lV`C|zOmaOOV!|N'
    '6jJ#c8<=u5^uvEU^$BytU+Pf&pMN1zQORdFHcK0AhyBo9**osl$Wt%UG>Lu&6Zk+`?UaSnnv?q<`5NQ2M6DlvNr%ctL0QYSC8jmN'
    '0D<PO({xv?lqu4dsQ%{xkayZ$+X2)60rNg>bjuqx2W&n`3Kmzp<~}si!RvjEeRQoQ!bs73U8zGo!Ka8)<?FfkW78Ge2mKs+gmW7#'
    'oxV!rkOD8#(;W)Z8?E))HkZEoh|2b0&?l57_~d4U^H6q=Zt&J~^j_2*isAV0s%%wEPjxYHW6D+!d!Ey7kztE+1b-qpvIK62CGBx0'
    '^)+j&j*butU#H5z%~f~j%7d=tc_?Y&s<~#M_BiaEw!4hsWSi2yp$s;@R%+CEVCQZbJZJ_>)?02j5Sw2wnQ8ACIFo`T%ng$C-5W}f'
    'Uy~rK+*X}IU?C2Cl9lR?m8J$5W+2pLPILzZnc%hSFWoMqrt0=iGp6%jXkBalZ0?C3bHXoqmhjjLjJ{W9V}u(tlWMuYhl!w2>&6z>'
    '=&eTLwJw!Qml!eEHFCVuRKDYaP=`-16f8NAe(%@IPZ2LI@<sDeX%e9o@;w>u1e$hQ5=yG`G(^#=m$zUqa1`@OqHgOpBkbjPS|Y6F'
    '(yTyez7C{GbcPUe0L8qtDb1}3eEW5A`;xcqD49~>#(k$LGlaW4@nhOu2QRHT1^u?Tq^Bi=wq<4>U<+tUqWB=iB1}hO+D7|t>#@$s'
    '_YgONLRy@zs`j}s=nK_*Z23YSMg%%E^wGVP*CtwriOGem;>Vg({38`ECZM7c1rxwa(EAeh3Z{c{QW&*9k$b;h6p3#?!#4lCPg8nt'
    '-I!vN#Wv5L`I2U=6&Yk|3E7q^n(~8~@x2D5Z2e>HtMdrtw84#j&fd3^N5(ctHMB(g>|F;VY{(Zb6Q0)ium!?Yg!?au6=7<muq6+Q'
    'tj<1}4KeMoBeq=91$4X_vD+dpfb!&^Qpp?=_4A7BJoxU$2;1|~4Qj5Gc)l3WSI3uO>t5G&Y<Npc=N9yyc3(O7raH0fwk3^QKO=ZM'
    '@XHT{Cw<d1jYn?26t-QArk@bT&!mHuqQF!^p%l(0FhjxXCXkYs@8Sd&T<`3ucE7uN<e1ok6}mtit5v*U-sR=-Kg-uyXzTcjyb2iW'
    '9mtnun%P9=>w^P`U`4X%K14J7ZRn|}Zr0Rr9Fn*ChOTXid@f$1j3r(vi9>s!78Z8XuB|I&wpv#y^#mog$tp)Drf@>#ktZ$L$-f`*'
    'x@uF@TBA$0O3?t$^4Lo5xJ~#+Z4Z=2*iiObuqhlYL#7rzn31lEtCFJmMT1YP>Wz%cvOT)FlKHMI=jqTwGM2S{W{^$9L*_JR@0m5_'
    'tqcboi$Z1|Lg_dxK;teeO7gt2iycjGrYUnGc<<XpuQwkXbFai_h@%iWNx^1l!iu^3iCx&CU{X?^!R-$0I~0gE>_}!h-t-Th-2v1T'
    'MdLtlu-TYfq)4;Ig;c+Cr5VyfY>YcN<qLJjiRlS;YVoKjCeo)|3mM^ry#%}P59$BD$e?irIj@6}gB}*w?TmO{i^iY>t+Ov&cWS-V'
    'Fz@rY9}y0@)$ox0-8eM5)D$1~YmSGO?Q5_;*G$jtK{Mb-3?)m9u!L>?hSloR`+g*Kw-Zbb7<U;2bg<=Dw(AiIv-75A@*!kMXbzeF'
    'TylvAqhE6A8*dawWs3)6YeS=i*LKu<>WhEdpIEy>Ip*I)7-Jr?Xd&`KDCcnIi#CH#7rVO8Nhj5$HR|(MvzUyECzx;cyq0W=#@Z23'
    'P?hUDJv8pZ*)y;<bNACiOHd-H#rc23hej>+nHO^WFYE!1XrER}HZnd9XlhyMfNTHHiwlS;0{&)LsfC&KZ<@|64Ly+8>YS_n7vqpi'
    '=9)bm@VAA-=f0jH$zxXY(BLGhFGf^xE#YyjGLlS-dJ(F);y&bR%rpI=eJ#X~v#qw~;Czf;<8xlq9PNyvKpJdRS~-o@QswH0ZwOgi'
    'Z;0fY%`+xlMKXKko#Qbalm)e5l6jsBnU@+%$|@GoR1s>B$@ejL?y~Os%3ZqjZwC|}f11rPRXdCj1H!+%383EfY>cZ)6C8mwR=mLo'
    'G9K}?1URR{Bg$G~+{|aG!gFS}y<B-_RSxgJNW0ZHv)`6SotJ==ewsjVYv^EwHfs9wA_3U1p9OHR0HkhRv;F7Utwe8p^1kHs57TH9'
    'U-|b8HGo09*W;h4+3NDV`7Z_uQwk~ax8&u(-rfQ%?_eP~P$kLJ$2&j7`iGSa8Zl&PQ9GPJBl?3l&U?~NEf=V-s~5W3xxwM+nTuGR'
    '*m#_|3=`=zgBP8aeHpOy22K1T0U+6$sXx!=!4n>q=~aQ2OP1Pgz6k}6`6-#e(&nCw)tAZu+!(QoJ7v6bl{+xQEYRxd-_d;xt3;6W'
    'LYQkfTVGGhG*4{1a+pGkh7S4tyDo9$MnlqWr3ym-^<_+lof7qE(=oBX%jwji?f^)L$Eg8kek(F?tIC=d>|42UJB5oH(Rx|O06Sa#'
    '(>ur_qC|uw{T|97V~}N<?Dq}1>Anx+GZQ+i(68S<$NeFdo~LzAe0SiJ?DJ`mQrtR4@(aLHm}v>kDnBy!6>eO%0V1il)-|e6#F%~w'
    'YUzWh17&UjIMGK~$yHU8tsDo>8@J<*Wgfa{n!wcs`n+~}2@7O6t+XMC(?{#t_xDIrgv=rI%dFm~O}LXug5CNGMLa$Ng__TclI@<g'
    '4%+gC*HhV01Co6R7UGZduRk3U^Bm+`Ag(`+yo$2K%Eka(5r(RS!5E9#0<}yk6)Z|rYS#jNsK<6t{XbzoJ&q4fhTfLU5xM#*3lGP%'
    ';n-)4YWJRw5Wj$sc17bZG9aYMw~~o)_#72LfeZP%cCo7*Fy3?`$5({AjMXXM$btdZaH~XORk=N5KWbTOniUEhXXi?o$%+=vcvSuP'
    'c-9OoB=#oqMhFDJi=~*N;)z-^O`mLO2J;e1;;vM=a!Y&kmdaYjI&drY=1%&B2otp-tP=eufEk-lkCjLCaRBl(^o3p6et0zxEW|Cm'
    'QgYxQj;Ln8EvQtOrywIN!x&`Vp7RtngX#O=<D7s?&F36WH;=KppOmTxp1E;$wNpaqsNr&5t5QFUbYg(4m#hnk58XpkQn+yue}zZx'
    'lXJb%u?~6*?zcLuOg)@OTxUP3{G^_n#(=@B@mkD5rSyC8#k_4w0Z!g-)@V>LwaZa~rE{ExrbMNtDB;9p1IRhPdXP(Wlr{^7rmUjw'
    'IK9&X4=z?Wcn_1-0M5}N(y1hfySED%H_ZXm)AL8kbzc$1gW-yuNz(%4n5x;Xf}&r2l9BU3pP~`#`Y!uBL0vO(sdO_x%iRgO^okeX'
    ';PSI}GN@`WS>(;=afi~6I7B_lRk3u%<fHu)sWf>zwFv$n65d=&+gf5q_mhX1Tek@=VZ8R~I2$p!&Vf=UT7?}t?5xW0+U$}4m<*ie'
    '-raSKsvqwUZINjR#&!tNs9~BAq8WC`wq|$R6YL5uocpj3ciT#HNc>RNCa)e0frm5kb@Sa|+J%v{8}wa)FDA#oCr)s%U*G{6zskwx'
    'x&FGRGW}>`YK3d?KUSB;Ag`Z6%r99G*Ly`F04d1M<(P}UmM(}AcBK-u*gzXd>23mB#`W~W;f`A7Q4`Y#bbX=1u~))N6E6<ws3Zi;'
    'X9gERVhKJCTJJ#ktXjDkOG}^F=?%`|wIdOR{*Aw*--TcVw&X|kH{jL)9)m+M@etJ7Jx5Q#IEU#tu+r9{6>aN9Ip^V{V1E`wA+a4p'
    'AClJ+_gL|4d%+DD`lYA@iwGC`-JI$V%D7eX8Qlp8T^+ymZUU`41S_<5y;0B|Q(aAF_|8Q3o!hok_e-@CO9L26YiIJZjkE3xE#Es7'
    'O3j9wmrze)5i0uN*;vUJIWuz!6M?>lW281&l~sIw9(EnIs=(<Vm(IsKueyK`$liZ1h}IqEp9Nn=Ib>IIYEdL7rB++$dp=XMJf;BF'
    'Y{al(-n2YHeEa}oSB%Aoc`-N;_KQbycZ8^&N^amjHEI2YVkg$FkKLH^Y0<XjtYvkoaherR@G8FAxzTOJ%UGzsxy(#Y{-{wv+!So3'
    'FZR9T)mLhfK&ppDVB5eaypxCI^<`m0=a1%Wef6Z@a`0(a^1yXf`vw$v?W+G-@T_<uW>7sxprzZ&5Va!|iMvYbBMV5<YlGpdciC?8'
    'Lfao5Ha6|iUV2eI%EdUL)dSAvCqS?aFMYacc=;Hv3jmO3?_&JkP_v2y(I4Q3dPSMdlc>FDF}bM_{vbjK+nes;6giL94j}J-5j>Tt'
    'KOFU@h!Kn5$}D_UDWr@!kzJ8N(3o01VnQC>h4jv9RQ)DA^WIH<J@xDv5X_LQTnq6xecSL*ah}~gLYkiV=12W<x1yIAEmv=i9=`>Z'
    'Zx~cjCKm6W`=zrQtnXDlg7f%fft%VX0y4W@91=^yE~i|&p@~TC%O6ezXC+B`%+A8ocqxAVSbPhoBwpEhCMpYz<>NH{d<CK~k*xp;'
    '0B2eDTn|ylz$Cf&s^LrM6p@b^J*Hv3RO0}nrcHjLcPFY?C&z-E?CH*v!zw2{vh9%RSv!P`p_?eoGJqq&)5i%=6NyjK9Ys?K%$B2R'
    '6s_F$2=|!AKyi!0<8{Q{51Srqzs7H1L(Y3nkEk=Z-menrxi`#TtE%jBkI<Uq+ttyqxln#SISYyw;a}Dq<=M$_GGk|&s@(1=C1szt'
    'D2hOb>b8_kn9GF2u=})N%U5?ls4lY&>8(j}yEA*&S;%o`$2nyk-;}E!tPVsaTwTO=U7?ZhOWL6ekPx`mWp6)z8Qio{gW^PJ2T$9O'
    '@vcIFhpmtqo%OQmktKv%2}nB(dCE15-5SUy^PS(wr(ypt;pCDN$)A+o={;RYZVR2}palsn9czpA1vDZTC|zZL)Y`N)C;lSh=Seq!'
    'gdwC9la2WjeQJI_r18I?aoPl+y4n)R<pX@@pR?-4&)r>r{oUZD8bwg(zhr3znqA((6HY*Vcamap<jP#_>SA-!f3yvM?;-MVFN{1a'
    '|EV;aWtEF6uids4&s8iD^DCGw_|!#Mq{%_2to(6o`XxmdEo8*rHY=bxskqRQ7GU@1x|aBkSY;YJSc$q2*!y}HOV^jv%70ZKc8^cg'
    'puFX-=1f!B0FE3!@eoBK6w|yq|F0Cij=}FMY9DdIC@0w$hB3FingVU|gdeQgd+w7PRnboYDNELWm}DDjAc^HBjWEI4kdJV+=JmFe'
    '{=r2?m<k}E(>`ttxTux3ZOx;&Bo9K^QXHJ-lt2I^8VOIZ`7V1vPK%tLv9evl+#^(8A}???ue~l+MZ61sT!?$2pgGK4*%x7BnE3b|'
    'QdYbGGF10*MGF?8W=zo^sqDjizm=&(l9x2$nMa9IopaZRW+tP+nY=cRn$J4?ODJ0TuSC&$2}H4pOE75Gtk>15)MHOeu4lZNZkiem'
    '-Lqz|wIr3bh?IcFc0W9@Kz2ZDp1iyA1M~S<-K;3t_=c$UnxnrldH>f2=xVYsZhj+`RHF)k>?~gikLiqcGyyEFog%6@kVLg?e$BSH'
    '|D#H>OHp%AS|9;E?Mk1KbNy}d<Lffwem>kB<peO=cXI00<k8MKzaP=Dl)HdSQSB*cGly$F?N1SUG5sa%o%(JAw4AB^I<L>`bHHc2'
    'J9(Of_M988dAHQPaiFvc$P|DkAt%-i?*|Y>5GtZ;iFg7AVCVc!?u274b6{^)TAc!btCO$O{$eH+(w!ct*VcY#=uQ*2<)h{;uS`jT'
    'wwXvqNc1}$eTRvuAsDuolw=`KgR+8BP0e3|1Rc=VpV{V*(_j^B>@4M~D@g$7nONSovizn?kOdC>1(LU@`Kwfc#IGm<F%qGLubR`$'
    '8Z{~9?Luf_-`!$<*PjLxETf3vslb8;+ie<!p4*|8^@+Uid({Cp;WqBJ7io~5udzVE;cw46e<&U3o(W=44bm?r79ZIzvawL)Yf6JD'
    'ZHBJ2Tun$j#bd75-Q5W#4cJLNA{Wy5?fA&0QHLxfFv%JM0J>5w+8fn{mAk<+?a!vZ8vWj~5Z9I3Vd1#i5Z@1&2p#Z}^>&Ry&nRee'
    'c;Q9!a<lI9l$+p$3=r+OcPG1)RNL~y;@t&L6H{}RZ^Sywg2`jsO(Q>BxBVglD(P0OZh(7r1#lS4d{uu~%L(30R41<0R*0cq97C@t'
    '#BmBgPRwKdFlDZ6?3!Bh-zk1SwFzG0yKK+qO0a6RtAs!C268we_y<bOV0mXNGjcs$<<kIT@^AyGooM6g_!K{5Dv~t4LZW?CzqOh|'
    '49xl<)tEk24|h*1NC))pdV9)2mO?53kDxF-Aykn^L#(B6yKm=5S|-leJryRb%1B*v$y!GWsoCU)*<E(!#JU~go9m3WVB;R}fb<F!'
    '2w2&5KO5C9QW7RcDNv6m)({@tv~d3R;asye%M@4$Q8`yT4w@Q5^wBn<J}=v<#SlY#QmteNZ&9MHYu(J5zd6Uil!SLjy}NgR@<rrJ'
    '2|#;gmlCjUD2sIO7Qs%k!9s9N`Dee=B)ix<lWo#N`nG7k2=Op$e2(;l;7$^xe=XHt-kkWB<89{)A12=dy)npv2e}$2u?8@`80p%+'
    '%W>=V35w#+C0eiJO!Vl!u)?K0<*-I!(mLG=G>DjlR`}j;f&Rm{B0^0Z1IB(O$lQgH=bLM1HCW+M=s^4Y?$Xn54**4bmDM<{k$_!}'
    ')K?Sb^qDn1CnnTTkzVFefEu5}$bs)dVBcT7FIipz%c(zqNyI9t=5Rh@{KKPG-6aZr`(@mA;XNf3vnnPyq>lsV3Bg^yN2ZB5!4=2n'
    't;_@M_UBVLI{lcdHKo|7(#ucQYeX&?zN@RhBVTDY@@{p)VqrdF+LzLdba+kb)1y_<28UxrnzyhGJTLyT{qm6Sw2X;9z;K^=EiW{-'
    '%Cj4uwRQ$wCy9q_Ho0N!j^Z;KjGMKoIq(F1;r4X2VgywYm|s~jj6h(oqq-OwV5G%d`I`9G82F%Fxe(on+-jh^?Or8$m*P!`ce&|T'
    ')|>WuGh1;MQb0fX6;e!L>eBN(`Oy3ermr@bRbOelX*x6bi?S<Ay%WmeToh>$MwmyszeMh>JVpP=*-hJW3+wV6GF#bRgh%wTiphOI'
    '8rS@h$hP$&*0g*53GT;$I_|PHV`hTxR8y~)QQZJTw;DMg*v@cKH9C=!f<vo{cn5Lo_n;GD%<K*JV^#TrQn#jUD_#INolOuNsKzfp'
    'sRZnsHq@JG8If!d$*35o^&BF4ex*<9WIJM~{_iqEXXR`HQg35;&%+TuAm_?%;;W};@Cgx+v^j~qYNSw1PK3_F1a(QOtP9dmBpB;w'
    'c-oTn(@99V&X_4N4<nwS7gan^HM#R6!>^!GVx5$~9L{}P%ibcXJ83mkPr8FxBY)=k(&}tFH;V-XoDEa9;yPt=L(qoF9W~$ALQtuF'
    'Ll{q7SqcW1+SaX$CcxqJNOF%2e_<*o_M(WKENQLpcypK=f4YPk&`D&x?4MXhfZNx?yTdriYj5!{(kg0d=Dgnf*Yz=pVmP%giCa1?'
    '7Uns)ju>@21Xrrxj|xCf1zbD+J^y%RL(x6%L>0HIe?6OAw}~;BLOzwc-;B%}-X5IytrHNw<wnl&+=c%z0Hwym><bj-Ngk|o(N&jY'
    '3{GIdikwe)uQtjF4!t3J)@Ugd?GDZ79$G=~r(PrEYrwC^yY77IrR%7xnfll`(-#4o)b_1ZXU8YecM9@Bz6AC<&8dCE&YF)l>Noyw'
    'P!QR^GxPf_4{q7;=lE&zf$G6~V^45diAX}qN<PCDTuBw~Hs6Zh8tq>|2e#r(svmx^U0_6FkeJbZpvN%>S)<kk+DlmNcGt-~6*tR|'
    'IZuqYY3-=sXKv$ma6zmUGjO+u94Sge?e?PnwDU--raxF1W~`IaMu3<7jOqwNE8>CEug(2~W&nZ>X0>WR3y<o!wT@(l`rT<IDAA?0'
    '^T+&@I?}S0h|}}aX8}?9#9ev8{;?O4j>Yx7F1MbiO-xO(v^NdS=POXjm8sWwB*!_bizwZGuBYZFY#tAZgtu65ZWdPstW(TrWfE98'
    '#tI*z2E+SjFPjv1M;p8B`6(9(#yLcFMw%0!qS&pc<`Lc<Cl=F?f|g>m?^Yq{7e3_41o85V8E=NPW?r~6BM0xtojo}Ys=^6fyI-lO'
    'YAw*97*l|5zK~bj=!TuYt4|8Um#)aCQ<n}R4dj{>f7MV?0;L+#&a=F;7O>Iu^`AnSHu?sa6L33#fur$)MduKkuM@xM-N6;Smh+=K'
    ';to&?Q1_m2i0bI$Xh#AM>nPNzYPUOC&A}}HI%cocFg~<j)2|n&ax#a`J;DG-$lPi;W~_N9f|1%Cf-P#U<=$qN|GEU=o+{0;aKAhC'
    'p0C}}63E@#K!CUrLUU5W#-k}A+kx2|m%dVv-Q51HOGIBeLeIAjVFA?x^QyT(7=ms3be0a{TgFp-A^-e-C9}Cr;*rQuT#P|c4qfw5'
    'IRlDxgo_lXhK=5&NOlm-J8s-_*zcG5gVqtug9p}WoTqGVKktEXhREe@tG-}vqGDl(Hmh^*cNu#`Ja?Z9fV`!2sa`>@;|@Oha>a%1'
    'c&pCTJ^;h{TAp=$F@p(}f@-l9j0@-&Im&c}x2vpq5bY|3S$Kq-b2B1^pxCj<X>qIzkzG_tUkCjZ=C3XGnc`qN8DQ*8nAt4lXt8k_'
    ';{NM!UCvammZGqDH%WvvP>){vW;~c;Ks-5Y3#9s=lNJ+nTFOl*q{!(4(=bl2Lcw#wkZx8OsH23zbgM``g$(IO#W1dV^rhP4FGL_q'
    '<Jy?h8piW6FCsbot7|UG&Lh#Q%n>4=3dVv0OG5+zz&Gzfs@bxPixNM+oHG6H-C;UHG5MM?Nq0a88tTzRCN4=&{B$r`%WXagaqvWI'
    '9ofrzetNHdu<dn_;;qc1CRdQXS+PLvQcN2HXKA^xU5uYU?bz)Ww*^U(Bv*M>#7JP_=89pffSm9)>HwI9G(Zk7zw3DN%e*8TaXxOA'
    '`;9*ID&<2H`G!C?_Zp~_&jbSHIG(9ZxP~tQVR^b>mV8dhs`=@$#yI!Xa)jn%RvZtm#4{WH92zVfc!R@%Ng)tOn4Ck;t?a&|YU130'
    'bqqX*6MfN(K>V!KJAOdm;ZiL&ap(TF0>|9~u`qcK{ht8do?EM~=5SRVbo`0>)Zix$EaI9+09Q>j;}Lw$o~ZZ<nrn8Qo<<WslP{8)'
    'Xvd^mk3$sI*NQckpqJkH$1K<p0OQ{M6GY5J4A6|u6tG|j^d7h+Ka%^(?D{uU*W{z?i|l(a#u;8cyuI=2G>PiNU=};>E6eNgvu?~2'
    '09-bmk}!4ai=;BYu!X8BLJMI!DA_Bw1Hocu1z?L9)?;xRHUJN_9GOa-$LGzje3nJ!=JntRz>uGUKnAT`NARqk2b7xzehB?k&Grr*'
    'xhWa|e(b&IKf*ex^83!>bd?S;@K=6PXpo)=ie?H!6n8Jimnpb6I%z8xoJ=%0>PJj#j@NN^sdM<F!zj+ma#Rh{JtDnQe7FEy3HUnd'
    'U65g6epg;0x8xmX_`2<hTNZg82i<(`!L2cNI}2h!8lEbh@1qyT240NbCcgb>+GjA;LyR?&jk*>NjWs12H7Mk>g4ebCDZ{Vh<6>xi'
    'lj&BEoOn%sf8U>SAh;FS@eilI(uLz@<z`T>YCrM`o*hC0c(cyFCR1CKi*aAL*ibN0;!$dV1k^L(B@YUkKX?D~?JBnOH@984I_p@I'
    'o<YbmVOsAFL&n4K!lM4Lp5qRvQg$_g<hINEs$PB7okx$HZ9nuD>E+0)u;7@baq@-h>N)N8eCO~6ro!4PZ1rAD(|+K5iE0H^?(iEF'
    'G!WW_AV-$&_fg3oBRfw9%q)0S{*XO_U6%_#s=o4nW3)*eexI?aX*nFDbT}3EFO-1~*M3JMDH{naf}bUg0iiVLR1-Vq#P$_7T*^hb'
    '@BU}(x*<Dhp-$oUv&jb5TTYYtXQj0AveBjUQYJ~9-VTnJQyO=}L<>b(rd0lYLJq%4Clh}X@6VM<QNYnzlr_BAO5Ms!TjH5kPhHIb'
    '(C55Jc6c84C8EcCR@up|vpBXf|EdNTS4>l|*CJi7Tj3KaE0-yMTS<#9Ck{Xss}_6u7s>4E%B^t+2mlK4Kz-2TBRL3D2gnYNlGK*@'
    'n20rj^mA_k^|Yl$69{s&AZj$#ro~gGuTBvsKL`HA1eSc^;sO6UMq@@XgJ6(O^(fXPP0lEMZ!oyc*R8T$X{(Igt{Itl!%++IJ3n^l'
    '=7|@KnC+cUKZzKv?=*DUWh^3(Tu%$WAfZP|DGEQPy+rR{iA|ai&!FCX$^Vhlt1z4HJ^fG2Y3iFxx)e2J<X-N)&9^V#CIcOcW7W!{'
    'xLghVf8Po&$-$)yI1z{#b8R3nq3LvVo&`i)t*AekV(Z!x+D^=)VMPR#8{r};(r?fO`wDNAc@~RA_@Evcc9owR(U)dtH9}qKM6k#1'
    'qTsGV;yFIKH;ucmn{wQ=3D}khXml%iXmPv>+Fp0K+O~#4lLS)WAulIEs}8^KgRk?!I#?yp!xs)XPWs)d`mrRA-y_}ZZA~x7)7D@U'
    '{Yz&XZk~t7zWq4mxOwWp+4=(Hr`FV=Q`ETD#X+lxu8Hl-2eParzyBE{N=~M`1Im{dapL<Vv6sNU;47!c*E&sK+<Z$xm+CnbmnQ6i'
    '`3)H9zr&BeQboh8*0lT0lj~M*BdY)~hL?&LY;4J}?P-@?X3e{juxRG>nlLqdb8b#DpMUtBQ;2qMzoetN^J-UO$I|{UhpPfT%2lV>'
    'xHTk_cvtmeaJpAEFJQ4`27Fozf-!;($qub*5@HAlr{)1z3TXRvQWHa(fZ6N!F1GDRJ#I&JwXNU*Lcl3n7BP%MkIWY#y4=wYcDzI<'
    'MFXTBszPXx0zmPV%w{&?3!b39)JBXu6e+T_sN*>ramb!!tqrT{v^`K=udxR=l)KRjZYFKaDfCM2XnRP#Vwl3YVurn9ZKBffa(OK7'
    'RNY~D(5XKVNcg$*8=`S|-eQd4$QfyfJWKg(^0i~ct2qJrAeoMsF_ryQvHE2B^fWed9P8`J$Agu1Vt|87{nVv*v6yvR_?5K^lsK4S'
    'V67^nwpi(U9Q2<4pj+8uR0+X}ufKiKR66L@a0C#Xug2CBj989R6g9gxpgyHaj)JE{ZNC)TjtM)B@}JmEC*|g&Q)g3_O*zK}6n%an'
    '0V!$io+5}4e`j40*=OiaK30s$UM%`O9GPAPswRa@iPS2-GFWr@$O_@e{wO_U{aypI_}?u2YOirk8UjMQoGM0o-9pn#keYy_VGEs$'
    'r@gB+Y1H60@GK-fPMp*Q6}nWb7PumMjT9LQ5<<dZff&QRRZ>6E?CbeR+8(RxZKABwyJv*Sq+%W7;UX^aEcQb=fDNr*L`|GWCU&PJ'
    '`kqLiU(gb-el3jWTHwhXoK9sq+%A(Rjl$Z}Y;f=X0Kz*u_KDd(^9Ex`8HLsqMcRvy_`avHN;H^&+Svu`Pr!9$^*pUyLwT}9#dyM('
    'pa2a;v`)bhfo2<Ol)Z!nug!15Qkk6B9{NvfHer!vUckNrk8AQuCiKcB?0$%+$#^+d$C5M=N(+lBBAFRN_UccrI*`T+@O=Fg#6ocE'
    'pBt;rqF^^3UP6MVITjI;j~C=wnaM7fN|a0lN!xhtbh7ThXRz~W69iX<Jp)LBDmAg%9F+*CdbDgSn#7Q7M}#Nz2dhmc_HZ+OUrnIr'
    'If3paMxA4;JWt$4Bi+-u1Arv<pRjZ}es>b@_#3K8PF({(Va=bM5zI_I#*S1M@1)aN&!HY_I`FE^L@x`R+{}huVp}CE3JRN8sqBE*'
    '4F;u(tcq2*KyXzvYOrfYDU@37Jo7G?VFYf}it?F%Ka=;1Q+*>Wcxl5n{s~VPU{Z%Vs_j-K{h+0YsodYCRtui?6MY@!*>Jy~_y*Dg'
    'sS*9xYP-tSnXA|3c!{TN`=F5ydE|m@ra7<2PZIK^{lR44Vns9~Dwy%HQKJ{HB>TfEKW8b`ARzx-@a!I}fpS%Qz5Xq_Y#O>vYxUZc'
    '{6(NV&{l#M*!-+s;hU^b-LIF8-IVHTg>-E)I3ZjS57u-0ByjXGZQx4TG)Op^oEYrUc`eUMOecyUy6t+?R$-P?c{Nu=!($1{z0EcH'
    'W=ngE-v(*uIiyHyP8I=<FL$~t0jyBBT3^ppY9sw2TtqqDVr$=cC+~P#mRzGc#MveU`%mk-7>E#^IHZi+^kbApI`zHOC|bV{>%iG8'
    'SIrrQ^(BLl6l)X(e62!FZdO_AD*fcP)CEEIp{GBy;z0ZXQW?j5-8LK`KMh!2k`<Jc0zmGmTcmZjL38D;*4*47m|DX*R`K<^*`l2#'
    '0?kNR-=^7+bsaab2ewc)lt1nccf57kn>Js6O2sSoc)~Y&WlQHoT~Ql`6<y{zTXky;4FXdSZO{myyB53-))L!x{#w4-A-j_m5aQJ&'
    '042X<_H7giRDD9Tl5TAq;^&wn7g;gGEY;%vY@P3$g%8QUR9jymP*vLUp7PINZ-M>PCj#m$9X<Xxv|Q*n_mx4jA*Ytz-;->(;WHkU'
    ')K+w-!EED%hqAh@l&wnJx)Jqhr%@%pAYJ7$i3<AE#|8zDAU)o=WY&z<X8?JSa=k9vtp`)rs=V%)3>v<IEZHY6NQ0WpBg_q7or@SC'
    'efIYVSOt7`l(6TbHYbH}2v#clAAM{D0FrY?86qK98!av3#kCOZSGy;1G*%M4bBAR_U9)3`Zo5FT`TsZRNWG!TCJSxG-w43d8zPOI'
    'HCBibb`I;DM7e|k(t9JXm_P;oad6=dAZ871Hai6+LRl}_XE+Wt$jbqB%&x<VjcZEwRsx#IT&nVD1*-=d+OOOlg_m|+R>K+rmG!}n'
    'Axtol694kaxjZt{t^>}rHpk=`>DpqsVVzIsHU{_6GC?|Wh}p6!=mIv$EXh(ctJ1I6g@s7g6&G(x2+&RP<L5jwArR#}qxbmr(8tB9'
    'W|n}q!;bu_?e}`5zt_lI82-=swn1i6MwkZ~1{fMlAtRBSjVEb4<*KfFMkueOeT?A{uv24flfxqlw0gJ&bC5)w%RJr!lDHTNC*|EY'
    'n2D&9seNlqA!cu#o3t+wAI9$WXEM6GZrS<l#1Hu7$@T?eTbhD}oQI<9ZjWLf4lwJ!m`Ebr{6b7W<<bG<bvf@Sti(mVCvg{;6hen='
    'Pzl5?qk?Sqw~_bJ0jb7|bU^`T09D-3bsH0PPWWtkUfwildJU+uT~)a{@ydPoj(2-3nT`Nhu%)EuY;odO{^w&)V*O9ly7C|{-?iS`'
    'hEY`~!?jP~i*ZIqd~{{JlRRCUieAW}BebT+)WL$P%rh;IPR9?6uFRGDPqcH#>DdX#pc&nx<RBKFJaJI&^~}?pHtJCx&|a-!z`VT4'
    '@cm77;+gnCKvaKE@*e`I#L~@!m77w(WSS%~%;Pmrizp)qF#rGntp&rns5EYt00ENf0?d~l`Z_EJvBYQl0ssI200dcD'
)
