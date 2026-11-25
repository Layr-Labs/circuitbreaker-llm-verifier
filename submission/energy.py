"""Active -89/250 first-pass, gain join and starting-majority proof."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial, isqrt, prod

THRESHOLD = F(-89,250)

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
    ck(U5<=F(73670242301,2**39) and gain>0 and bonus>0 and 0<cost<=F(123000663514,10**12),'active energy bound')
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
    '8!2mGDGe$C8F{c723Ze-yLvUN)BvFB8V)2m+y>0#wMhRpt+bw>uB+L^C1t_E4xZsAj>!2q@iAgj2hP%`L!)k!xZ(2dJt8+w_$ZH!ng*!@1?Wd4f<C!U1+hACQM8NKs3>~VNLT5EXto1Zx&EuOjFEIqqQ^XgWN%hkO&PBSxOIx=TSecXw~CE4wS+W({P)j<Q#3P%IO52sDHWCT7%9?u@$gb89<3j<1xZUD!7ra!{XF^hy(bI6'
    'bt@4JRLIdwqEhG#XW|Ud-d772gDz00<^Y%8oUQ#N^aiaKAnesEIWG8Ltz=&23bFN^Np7E{62%Zzm02&v<`aV@k)na>E$tdqg^8LAJ7JuRzFd}By>}3ioB~(#@=@jqZD#V!uJJoRN6^@|LbZhjTzH)lPf7~3;>2Z6+Zk@@GlM2Kuhri(@OovV9sTZ;+};k&{n98sMqV;qzQ;a`QI<zppZ;ar*pCtDOWd}&tK=|?@|};fg>X9h'
    'dQxq6QWY?yh8oL_U?usV#54D;Nui#>EX}SUb3OuZ&SvA5>f=SCK~-kgHl?1Lw;BlIV)3w?hXu-K;!sbk7AMAc+Cq`brbR3lAnG#G*QlLx2S#e!zyN<=@VwDgGjL5B)Qd<jWhI4Vz)oaU2Swf@`bjr_+SwpKFLV9g1>Mk{f%%n4va%|UY(rLaq01~h;&XVyXl8nsy^7X%*~K+td2eC3lng{D0QiA<ls0)#-S@WqbnqBGSKSXN'
    'U`ovGnKzt<l7er-g8govy3TPrC$a4tx12#7=7hRX<)rC`s?}J`c8j{DFe&9F-wHDMF;7jF2*wc4JQ2PaKo-wFim{gE9YpE>togUefPqd04&i?%YQ^Aue&Q&L=p5kNL`H4Fl1CUq5HTs=2$9sfM)t0hE7ja2Vbv~3?E=z3W~PJ>E8x0db6L+>r3*nZntWLqHv}SWL*>~Q#t9jD7_vwy`m-QLQiQ0n=(RQC7wL+d+Zng-tTZVc'
    '+>t*^1lpwn9QVk^x6rqZAq6J+8@J=k?wiu)y6t82)NRD+=<Vu(x!K?N(9Kd2*+S1;K+whfp5@-;RgD}<Q6_yqQczn3TD@2i1Agvbe7Dt`h-{2bE!&X{u;mxD&J{-E6<h&f5_W00z*lCWFT&GM%Pc?T6fGNzv?*r+HVo8|!E=`8U`M_~W)ZD&S>@}xXijdUTLA*E^j&_EI<NIEUfR)~O^kHFLW<gP?s+gLRx3E>bB;(EuGVzc'
    'Rh(i(0GH+a14x7^BjDoiLjc3wJ|uGw{iKlUV*E3!6nl{f0`nvppsXDE%<SSLkCue1B2>PxVK!^oqK9W!#MkRyAF_F<cXv{+$#xF))Sf*_oG@A6Bzws5T_=v+nR5GZmRnUh@{rt}z~<s_*K-_2QG=K%g<A}`^ss44$3+mXJa!<Sl{O{Tak7+z_Rc19(UzS7i^UNyOJP3FSezgT&n7rRW`1oB4{j@_6<9GQMeI$}6g09~5`Y0k'
    '_5Vy3jd5v$A8dSjiL<8j^sD`qSkw<s-LNLTyQhXqnuT(PhWp}cAf<!<=3IpAq(~0qd|LgI_8F>gFzN$(q@Q%!4@v*zv{Cn~EBGMHpl6!Tca~=VW}-5`619J<L<}^B&7HfOI77>6rD~Tx?{ZdZ)|kYi9jC=UjPD)yusisc+O?mvjmwnXUT-=$G88>M^L5voM7H8T-o<hm2v7}$4PQ5#otAvdJ(aaAR<X;K)+Z8f0UM!3{DOxC'
    'w2Y3V1@t^=sh8E)mBZchx|@mkf^-00H~^H*tx1P1qlUb6!co(Mna~ia(p`X`IsfV40*HW_y!*;`CTdzjPH`5N!a?Ku4s-ecbj2YEH&5@pY!?;Y=vW-`a0Z-%55Fy%vuY^KFD2v1$VQ$^AwyW>E=0$Yd)yAHp^MT%OQIrQs4I4D3b4rox7+tSm>8d`+9h|_kA*Z9<ts2nAt;2Ei|#%H=n9(;<eG6%idC`SI`Z47P$?jF*1@%^'
    'T6lA@{89?w$>zqgIyX&5nzn4?ddIcLFlnbjYIZ!~|Dm!Ni0<>RVn{U7_|YWsp=XItIbDY*l$sSyjt0FBI}q)uy8gchBK2nTkve3L5%ez0<9xEBG+^Gq37mtK5IqE->_MCRhTuou6N7eC8?Gd*DvaV|R4o446BNLr&Tregugno`%BvXF*XtFby(!_x_d{QxerOa@EQwfygtv~lX$X_h!guk>mBZ8L>%!jt68*kug*e_lAH{r('
    've7GDmxKt{w9mQ<Toil!S_+3o_y`80sJn-59!dhYAot0kKiJs|<b0oetf+zIvGh4jE*T}-D#*sI!Qm3_XJ!PRB!hw<b4OMSLIyr>=_(-d;*26s#F9~-adR<73K*;IGK$_J5j&zqys#6>jmp`GRzsoJ9vvBmyT4!+%fZhkzLq#F-I&#oxjCTN(o$m89hn*z&D5-Ra1$aQg@2nBV27Ybd`uNSsz&JDx<yo&vt<{qfc}+0NfaZk'
    'Z4Hr`#3H<N$xV4b>*>$bSvYRG|M4P{ZU2gfga=IG-xiCR51WG4R{WBPBlFD2n8XT&i0C&cRb!Ost2WG*KhENoq_{LQNX6>sSD(D(E5=0Ui5cI8+_D*DpZc{AsX@S;R~!oW3MFv6^J85s<Sl!oAU(?$ZfPSFUPEB4APb4rlHr9|)K}+>rn7Q%l@t%$_Uxiw?fCN11H$1EZFDWgGL0ye2Ew#jMpwDA;%oq)KqzR?8};8Wl|M_{'
    'Wa$>OASRJVmt2z0)nY1rU?9e5xDKFI#uI{AxTjSA!dgJn%bQ}+g9il4z3z}-PTqKb<f(dvpu$sj^^=%siMVBGi{LL?k>+J-=U6u4SH0=eB<>tkZ~ViV!s*z`lb6O=tW5{AU9OsqET@FBO=_B5%;FT^ClLBkD@sg^V$b-P#q_?ww_3O%^Eek))}x|xP_4_Y<eNDZW9hMs`1d=ljVHTAP1psvApY#p8S=m+zSq!J(>H5MhIA!d'
    'd|0Z4!p)`qR@8pfjZ<uWFQGS0A6*lN?1#Bx+=yVLUeqxAsk@-hfd*xS?;uXxcSV-vVsoFblzI=XW-s#CFEY!Ubu84Xm06t`{uF4P{YLKwHn^2srujua^9lK&SsbVhLwPBcT4o8B^9M6X4<dU^`v?F#C)c*ELY*i^l<!Vl*?=IyQ5Qxnmdsj)KxO{>M_SQerP5QFJgIkb_A^#Js^#PCp3KaE(vXU;pV0#$64iRwGTezfadZf?'
    'UeQz}O-zSrXZ?_OyuVFb^Fsadd%w;Dlj?7X{jDJf8tg5S3_v@>j0Sc|Tx8s5{>=2S@6O7Y^Opfr{3iuGl<D9;n!0-6!^1Ci3qJ1xKaFQLCMSCuMYwTu{&OH+KH^I^Hk6Zxf|p0q>|^PUFFDT<WN_U9##Jk`Dj0t8r^ec1p^oRhSI{($a5o}1SInZ3j_o2$u;TiCAfoM4+)LcNO&J-(L+BPz-5qx|ci*EY5f!kni#uCR({WPy'
    'pt0FBmUguiVWHnyfEsUezes*jb~SU~>mPk#R70Yb^r&C){O&suZ;|F|`eBX;FsbD7@&#k2&BY=to6$kOO-K^}Oa_DVH2({`Dh$4a>~s8nVT#zH70><()W*@k8feLJ5L<dp-cK=D^wR++B!X8<jb^Iv#~08yhT@h0ko+@Y5b!<vsW%&m?r~c1GS#p0G{*~SoCMX-+`A*(757!-d}he*$ZzY&o{u$joRi6e38rV@SKIaFkqhAe'
    '%EPRM4RXBk)qRDs5KK75W-1_`h&>O%z~MJ%bn8OYh4OBT-JPDea_Z%5r08?+t*xrd^5`bB(PK{nMpWt4V#F+8%HG4aVHE+hLiLFQPxy%ug;h%iEkV{r>*scnZx}$k25Idf+(#^MQ;maw#Ec_5Q?M9!s0OuS_)9Rh=|ZsU{P2navsy>?>k{|*T1C}x22jtlS4Nu2AIFM~BTAz=8z4^R7)3<?silj#)Ft3rb<zMXG7dgK``PzN'
    'bj@Zwen|~zcg>m5=@C&{y!?RdU4mN!_YGxl>cl#4ws?Rh?y0FWyp(0ccRSz<e^_>jd|2b>U{np{gbFdQCJ)j<-?G*<#|E20NgsF|%d`Rlu1<+<3@MyE#R3;1#-_jLFkFJOyc{lOkK)%T9kTI?d@fQg>?%lzD92&%0#N<IFfS56bPj<jYRm#w>m#ZN#V4@~AE1VT;k+iRQ49Ep(W<~iQqmV1<Nrq&I0J2Rr(0LeLL&3ST&1VD'
    '?<IFm!0~qBcu%rw4r@z*m2a|f;A@|j1t;shVbiCn+tg{qvLs*}m%XlC^N`W7nZvS4g4t9nwW2aVr&rFQ$tH2F%-kb|YEr!i@h)*oR3G?kL$8@-Y~s6a!|BktyZ@oc+bSk|`-S0bA>q$+ypNysGllH=c}FAiP;P>LXTg4bj(EVmi^CF=KO?fT)7w?evF{RxL&bEL6LOL0AS`rWci-kxq?-W;OARe65CcxExryajOKcc%(l~OP'
    'Fu4Z;bmc>GxPh=Qr+2n=+#Gs}>R|1i^!H~AZJK<SWS`b8EM)fmn#6;Y7BE~svB>IycyBaU2U=YkDV%u-y-y%XibR@5<IL`45}0cL5CfvVYibUOK;W70)aFo$A9m@UC%DXNzDy!3iqGnej$-=$N7DEp6~ZAaZ)y_@O@cj{WKGS@>`%bm=S6kn{yb3G4|a!6D&2~;;@e3{Ysr>lFCs4gMPJtJKNVJ^Hh%M%B61ZPnEZDlP|mdg'
    'r?}HrkWjr0HW7AbhrB8XbCo&p>7%~vXZWX_?>fnf({LZI#RHfd!)aSQ(u1y#vcr#N%Gg6$E?22S^dl?BouYzo5M|&%#fqa1rs9EjDzw}l9bcVmxkX42DZZ=lgWfaWC6Jc`adIx4)bus21_chjIFki-mqQB#V$O~HD=i&=+vUq(V(T3LhnA6nIdUa^Y*FDGs2($nhbc<8S7}&8uUj%*TAougt0FATYd!s#Y&cmt(OT9$O(Q<;'
    '-{(#<xFAb27}lY3?-)oAU19##)G6%ffh6vF**K-+Yj{=f5-Q1n9IqGHayU@yOc~y+s~50BI5%x*v*$^Xb2R04)?_zfy1j{Z%wph7-tZs7ipMcc{RRgGJ{-6z*k1CMMa=Vp&4_c0KaMU)jH&_L2SPr>Aj-WsWV+hB1(r_Xvf0t%i)yq#YxC-JfFx&}0vN&>z!H}lkJLlZw874HI8{9->}CUc3-3z19n+>Z@wW|WE;Kv)#U>}u'
    'Qs-Uq;rXA+G=1Aej3yWYg587$Do+q>{yU$1r>P~jrs4wMARlIN2mPckgy2Qpmwrv2Yfz!1doy@<dkMkpZ8<fhtjOB|u-Zqr1IkbG-jCz))KnO~kvscS{c}Ixbun*WA%}!<;aBqBT**{icDO-qm-cn<lqE7r;FAi6_8<Mp>^;^vVch!f!j7Zq7fUbeR~ZymU7SG5z~eb>dScVL!cC|*L6?9BSK*Zh(hTZOW1|ei$FC<rUu44K'
    '!WJW|H4in1j2JEe)q0?3n5@b(7&kwoKpSGW<#Ju!T+&?Svj+waw8`j%`At4QDGwb2D4&PCl`K;4S<@&)xB7HLfc{RSw+XRjlU`{Enj3C@fKY?a`=i{pb*Rqz-7fU<NWFK4e1DWx=$PqTFKvw;Fxfe+5@4DJ;s0Zml;>1%gr^v`SBIX*uq97vWm{)77~#h)l<1v#i!HUW#Ap~x{cSp$@-rNFqMp2NJ6^;Ui`#n;ZmO30`0Zpm'
    'AXWCS&NVY~;x*Rdm&9zHb%5q2$fQ-You$imh(@2I9mOI`!IVBTHDPbn#-U{fBfR!S5I9qAfDe8$r)HKT4vQCtJAp8hgQ1UD?aI|+Jmh<0+X?Ff9W5MgQ;-C!(@%Uro_3L4z_k-Aq73QD2x9FVde+H<Q^-2e<<vKrKA3Y0znht)bA%rGsEl}|$@p&2H!tPr0#bU}g*ZUGolV;&dzKX;3w}#$2p<jA^p`y7+eSm8=%#SK_}7eG'
    'N$fG<?=UCAvPUz)uz$r6RwW}NJ$@kFvxpHf+Aq!Ciy9#l`KRw4iF&0+2)J$o?9%v(ZbDE8vP?z%5^>f6OsGtnq?0x~>%QT4w`NX>6;pv6*@$l0<ja2YgSbs3@`a%6E(+;7q2>!Ai+-yPllL7_1Hiqyrl|7Ls^3(L;7XcsEIzB}y<81aBUu|`DRpi-&Vbtu=_T@_AT;ORn{k}m5~^$wXBfyYJs~P7oU6kMqLJEXIDz*0ZotOa'
    'OUjKAGN!83bX*EKJm(&HNdLoK#Lm)F@@$f!m*+cPzYV1mpu!_?n%)6zCmzQ_Bqi5<L@?w(UFpQ+&2qVWJ^5nHgeGLukMUzeNO6{1L#Lw6Xl9+}-6$9J|2(zxl#I%wDJnmFTY1`3V8wD4a5PDStDIjeyDn0~l`y`>zkXMgc5T6d=UF#`#Iq}f{ZC>?GYXR`*AFfoOB<+k!fj^KC#Std|BTfBg%ND)_~&o2OXc-&G+R3;BR+BR'
    '^|i7p<F=ao0EOj<-&oU+0Lyo4SPdA}QtPq1uRJVHdm`rK&hA<zjZo`8UjOP^Rhl|jR}<9JaAJ+Q%f8dYL5Je>ka}zq=Tm{qMk^J4+x?~tzx7C!{>gd-N%x?sce}!e&c-ES&SlC+sqKE{KYKE(fHe+8VZlk!d>P&zJkhmRUz|fnv#nyvhL1=6lYf4424wyU{JKZ5Wq5};yrEC>87=W72|x7u)y{m-M?^DeF@jp>?{qS>xER&X'
    'mN2D=hFVBT?HOFW#GptqwSa{Ld=ld=p?%?xv(N`5@*8b?at##GOb9|!GEZ5tahXl<*4y&x9;wva1p#8XU?#4(UXJU$1r_Y#qIF9(*o($0T6*ug;-=qK#>N`{Hjt-wVRGIHz9Gs|C?=5cZI|m28#0N-g4J+sb}xHK`F3_IFz*dFIRJsbeDhkUZieh6$?O0)aIVdmD(^KBf!raF2^C6(7F+M%Ay0eB#c)M%&qgl;7mBt}P)%|;'
    '@(q8y4Soyjk^9g&O2t6w-?)0!X7iffH(qKq@<o_r+Rfmq_QQ&hk#Fwz&yTMKUTx_Q^lxJT27(A2mMBdG%tKTOxnz^(to<LPk9cj2Z6D;b702aS;r-`Kfw7K<v<H>s9QMd`^guux&4)aKz~$A85!e@WW(tnQx)jTNY3lYJh`L)ME&#8e>kk3pLC4WSlFGME{9BEKU07%NuFw3Yr5a$qbd)SwRi)E4ccLYB&_H>HyYYB9G}0vk'
    'Bjo6B*i72_(@g0}Y-0!U$G3Cfp2k6%6U&cw-o`@dU$5p76kowvZrM{DBu1^d#+)z<MwTi5MRdu?rRTXNr0$}?6>s9fRC(N-iq^7~$ym#K+Q>N0%Mo+7N*L?5{xJVTfy&lErtwUSTtu7VeE+WP_HE7g2xnWl8b&ym<GJ+uC#3jfMhU=N0B{76pE(gQgsr!K(4SuXAlzfXsAi(atKt^|7^UGM2T=xwokJKP<Bh8SeR&Hop`m2z'
    'C$ugxlO`-m&3@Pfr!<ti=mW~k+B5;CUhOzzSLtSwzDyqKq8Mc85#Y%iTk#V84G6SwEiLseiq28N<y(V~ekTTv^Xy*s=km-NRXvrNxfXz{<*!n1UVvLd;&9Fzxm5t9`~YC)d{!VyXc)Mb0F_=MI#s%J6zE~zE=ztWzL81OfkcFJYk_`sO-6=U0xv||qa}K5uQ&9$kX^ofZi7XLC&6aRRT4;^X+%4RC38t5bsToew3H{YHClG^'
    '>mpuO@S+FM@8?CNd~celB;Psr+*=QK3fc?CrW5|;=N`>+t3GJPpE7h>IQ?y(VQrZIoi9IsF7$FvpjMZt>GOz+T@|{L@}6j<_7vicWr;7;7j=*r(fjE13_xKPNV=<p_d|mX0VjDq@xJ9Iiv*9?8cR6lWeXe@5Hn)I1O23P*2H(OB97I+#^k9n&W|9fsXs_?C}rYkByW%*SqZg3-MiF=bYeHW-frUs(ZgUE5rEfrSO2eP^_?d)'
    'flT-!(*9=D++eb-UAG;)39=`D`bs9f-AOOuC$>HS)#NI&T`canKQOv62ujACk`A@C9Y>^GC-9?$Ds?7rYbn0*gDIf@h(eQjUA&!A)zp|P(T@omg)->bt83pv2n^DFltWhO`cS|ZzjQFK;tu>PN`ul>^O2)nwz5=L6W&BjDzth!bd*sg?l=ruBk`09c5!AbRd#f)Z$M_GLS5`=g`k2TYI&0XirVDV)npAl8=7cNB(4eh!{JFU'
    'Hk~GjwU-JfF@!}GN`uLXxWnKmINy|Y<mBnP7pxE)%p8GL9VHYatLGt^FiHWNXLW}Gs*@UnhqmLxRs_&(AK5(8e7<z=Or?BoV`Q}E!gQz}o8qs>x`M@k;#rdg9E4fh=o3OsS&<?PSktR~udBW%&NE1yo>Ix=^WPN3v<jbPul(%qxqvh@>!xG9?7@mRAjjeok8N`0oXol4Q<8sppTXWgD9Cln#C3Dw@Wh;&?Te~M&2Ot6H1}!j'
    '(`U)(6Z&<pIJw3VO-Jr|#`2N3QHGRDKE8kei5F28TQXcNXia3;J9I=OZh=eJ+(S?0L|nF_c}sclW~aU*GfVcVo}J4ZcBl<eh5xx@K-vmBS@Trhf75C#&42~YYgSM3fxXIbna?|9#GXa{UZ`B}Yu>3_dGS}X26cjR%7bu-)crIF2Ves%5A)_OjVs_SSeL94S~~(sn6N07IW^FgKfYgyjU+E`P18j<gJcvO7|3ZF>eTH`L*X#F'
    '?y&Pfua#^3p|w55@r4T{b3^?Xzbqz<x<L%Tv<ZQ(eu&|uk*X{ANctl$fK^j=L1UnRj8acjok430NAl|2e(=TBok<($fWGdJK}^FDtuvD2yvs&V?AHY<$SBDQ{C*tm97TrCaw0eN6;inbmm>~JU;+zY?7TjTm6go%5TxhtlUZx^Hdv<c=-%>dEerY30+_BFH<;y?2!#e*m$=h4`WLw*)Km}ay2QoKR1z8Xw4(j&upC?}3vssh'
    'Svviijk&qt1W{7k4B>Sgn~6auldwe+B1WXsJi(3&Ts3;$mrCk(>676V-dxCk0DE<X%zKP*TmnTh9Ji^5cQcCBr+Km@p!>epyG8x-&^04o5IRNrn8W=0xzR`ITKe!)Vvyry);M4VRBLQ~DYad_X3Gw&*%LcI1zzuL!1=}<h=qso2-pW+v(2|vTd7p~Vw-OIDFk>Z!NKK~w|IkovY{z*+94I4h7lMyq#ne;0JO=&MDm&te$2ia'
    'V^hpjDSajrIwrl_UBv%=B5;+0_ZpBa`Lu<HJ_9t;Kj6h`{2BhRXPxsRt{?sd{YJ2<D4Pjn1kz-r;{R{Kaitx>qV5oFH6~EI<ocQ>$-aKuR%06U>nGYQm`W+nGs+=6f&8xFx>H@2uMK|*IJ`?7cbOr!zMBZb2nEfACwbDIgVi$ru0obMqg5ML7F}f_9HdxnI&G|NtnW2dO6lXgNd61Gkh(*unP)6~VCf}7Y#*Ch=HLPzH*V*y'
    '43u@@#vV!KS=jSu5V~0?Q0GZrXfssy!g-J-0SY8xBV~YG7wp=GkasNGS&I<V){Tx<F9b#=eU65=!e;2s%?z*VCt{j~BfJ9aQg5k3hd4tR9%j$}t@=km+{&>{(7!|$cxk<vU{RY3MR8efK$wGn1iR{p^-k&yE9HFXRy0U<%*z;vqDh3)+=a!A6ea|EoizuR)soO(4^&IvQu@)yO$=W5Ye?ca730K7_g~?(;$<CIO!BQ~;Yz!m'
    'bE-pW0zsu0RFq9ntv__M1&Kl{HWS&(H7<j`n_Qggrq8v464?BBi+<ztu?Msg<?yj|gbpe0+uj2%GT}bF^ffzo(yi0g@i<Q{g%IsKs@$Q)uYAO7zj3&pddiMaMa<~&JO4t=x_A}9V^xbcbG=g2>xPcg!&;3{vMOCw_#65Zh=%@SR=!c9fm~QC<wSc<-|)>>yt<))dVJ3NOv$WirHwD?VkNPEWA!sk-jL=8{<wxNPJ8XT6gCWP'
    'G)RCPRF-r%3?nu_2sEjE_8foQ_J7ecyT!sCAh`gEq>k&TBHWN_zDvExJ_&)GS(#?Dj^VpQ*$Rp>QOZ&3U{VF8Efg;^K(GUC%+R=Oe&WUpb&+x`m4G=*D+d*Tr`4-h+)w)=-U1H3nHY$><!l}i)Nh`#tVLB98-Hn3q)8_fnPR$qSV`NrDe{<6^8+I$IWW6OmL39Hv{=M**-(`k5pl2Z%#{NviAY~+rLyJu;ekGg%d1r(>REnd'
    'o@Fu*3QOFPRa%$!rANYHR$m*5b*mXNEhT>4e?mtc4Es~y?9pu@nR)0nfKo-X)Se57x;P+1gj=EtDQR-_<{DLCIxI5ak>s<g$erk@-elSY^{gpx(oS9R`#ouH_6O{vx}~`@bcY3sREaQIm%Cs19*RIh-sYP(zy9&a*3f`UBW@XNAi%sraC1Mso6L|TSvjiI#L6dB+4i?!oNOq0v$sbjRMXkJe7}3@8@%J*(eV=iwkZUAL|yO6'
    'BrV(X-6=Aa!-JMXEQMx0=26|~hYcvYtIT&OXqwTyH|-Z5+Q67n?VX_*YrAC5>e4PRcW{jelH2-&FT5XoZwhzo0Qhl3$+P8Mi)dGhC@n*y`z*#xX0ROylra+7FAJJmGlujBXfFm*2IRw5IwF(1ziH75TMI)QD_2dd*@j|WJ0T|+$u8E4UkRMZKo3<a!p#zb*<^b%yCc%R#LO{nJLizU8Y4~U*`)O|bNz_zSfFf$y64`vi)wqm'
    'toFTv&N#eu^dp>-$!<@3tb{I;6+FExQutL&LmX#CInqpDweA2GbrZu@q{{8vTw(k8(v9?j?t9gVvXJW~8~ITHZ7j$^C6&e`AeYLNKw_x)$=eqa9p7QThisz_HI8sOJ^ksL<-fgwTU@@ttRo$5MNNameh8&cXF|R^j}2nowD1Jc<=?0nTsrJ*$PhnEc)G(_8hO%as5FBL8eQy5gTOGu+ki3b<bk5t`8*+cZy>}yEzmTlvbT=g'
    'hF>So1DNsoC@J0_$&&72%`4JFr%KYM9TdpG@&R-1Q~zDlz1|EXQXQQymm$H}Ds7BLLf&A9T<p8O@YofR>;#};YhjTkD?OBS6e+R0+1rNeb<wX?`R&b;Ean?`cOV_vfiGQvjVprbXeI|z2th$f1>q^{UU1@@oYv=9!K;9ql7ogZFl_(2MnTGjFp#E)zbbZL&Zcl_fVSY2q%oDyqKGC+9msDH^jNF4v-Fq|IAhIjvuoX7UEi(W'
    'iGS3xS}Opab1P5635cYF2-Xq(#Uk^jbU1u0(J(XU!|<f31cSdQi864;#hR}@;zZNT(9PLYRR;N`KD)e9@1kQZlBY}8;%1Ud{~XM0JD798^=_jxcw8NIuI5mu3m2&+vS{xdb(-Y^@xPc1s?7BA5xfQLFL^BPx6#YZx9|0m02Z%cWq(7zF9!YEU&Zp9;yS?d4Y&cAfobI>wzpbxlwa(&#=IKtmI1kogiYF7yts5>h3_cI0t>lg'
    'c*EEa3DO5)nRp9EeKsUva?q*d<b8KovAUm9Zkw$)kY>$TiNdx{mdse2>4YOG0(Wsjz*#pVIP~8NSbuO)Y%qP07^6p$s?=8R;EO3uUT4QIDUo+=@zZ^_!)2SU#W&*>-&)S%QturdY&%w$?l(R{Tj$>EK491yfPe-&9F*B~6T9$**nZ97>p<mAX?#SZirPAe&xpx)Z29~Md1P8r>8js!>5gkFQsYOAFSg$%=WqdlqS1GEJCN)='
    'gg8{&>GiGBzt^Tk;thh3Xa-VRy_)SG!KN72VYsgmh({ejE7g><M*Xd!lhcQuvusM}D|vnAmvkC%Gf>iw?Xi0B#`GuP2Yrbc=<pWI0s8@0Ji@>)IvsK4317?sTdrNVKPBR-_v2*<FXFWVi+QAJ?)1_>SR?(WA5x7TE0F~~|0)ajabj>3{y93p?;IDP=ON9!|0>H4jKI&K#-UdS%QNGRNUdo%p5I*SYD6pBhbN6$%)tJ4KRnnn'
    'rE4WpR`KV#dmwLDgu+hQsG4)%X8tbM$tr)#<N-3p{=VRvm~#!L+|3lvo+~%kg?=5l@ru)I7)<Hi+ATMyPtq|Aw%L8aJS~LM3Q*ZFhQ$K2sdyQn{v$^cG=3R41~A_17`kJnWqHwXgDXZR&R46{FG7koPiPN`aR4@&oC9Al4vF#z&Wp6C8_WF;4=7Q3V<f1++q8F-dva{T5l(OeW7C*+xV#cdOI?mWJrqXm7rA-s?Z4jZ(~u5W'
    '6+%++tGc#uuqPuqz8{$(>GgEWL8nk6;RTox_D}`krK(G}LBK;crI?oMh3Q1a*#ZLEON?jM@3qzF^GL5e!IADHM5SC$(acTKyyfGSbCX9A4pEPRS-xMGBpz-2mPJ1`doXQi4f4HFP+V5Hl;(C7xW3Q|PW~mn$pGn+MU|?(D1^&4;6>UTDTm=bkgM2bcp=(8&^d}h?zxHBj+RZ*==zrvw-pwY!Yw!Q#;{!p2vl(H)usg<o7^D6'
    'ZR}|1-8<I|&fAZ{avIVQ+&!r<3Po8}mxK)9<eHB{$7o=|Qc|uql8Bl5`;~}n&79`2C#|Qf0m_~Lt%`rE%n9xCr&IhxuV5Al^RHL`UJNQ?afl^~S)g8FK9Iu)0-A;VKGVEvaXK}ynM~(RkV=Aav-N<l6#jJZEfqS`rfvp3wcryWGpcf$V9qvUg+fMlckyS$q_(42PpC@{U6US)y>IxpuqT1z6V_@8QO)NwqV?FO@ctBnFtZb='
    'N`~+u3utl~(y-d}YMH%nGM5$XpoJ09!Yzn+a0E}Abt`|9quYYC>OXI5{;@ciQmE(ju_!OsYeOe~loJ8;K<je<`g@*neDiv0OclP1Ac$RR%dOHxzg143@&`Jwj||HBU)2TLGcb+V+Trt!@l#PPwy||p-<27Yp<2mdlmGf3c61$yRv`9=9BaF!v!8`CRd7k-@m<PxHT<c}NU6M4ibBGZ@*Igm>E>~Eji}=~5yi@PhTLRMAZK`M'
    'OFl{ArYmS7vbZQN3JLy?UptYV4imt`h*4S~9?3xFkJ#!2P_#m0G2nh(n6~BFcZo=nNBi%Vu9}UDT19PO>1<UJ0xPm|<#y+{XC+fMBdf@6p|06TI}~+|t9Nf+Reg>e%96hWEVh|I3=Kc3hu`BFglwB(AoG3CD6KRONs1DKZ=-(hOnnCgrgQL~LVgoevVzt|_YGJEd7q5H+A*MotX@mL{A1+fgcJVP{3Iz8rCh<lbbhw2J=beW'
    '*h8Z5=kr%n|1K3r$DQ#K7oEmOEMkWkM}0<zzUvqzb53Ln9(L2A_-Bf{ph0!qTzRUX$m&45YY&*vpsc+EqCmxpp2JH(3}Y3%7DLg`_v<#>uCxYNmj7oQ_7BQMpy|T-dApNpO^)WaDb80v;*?a@tj`b25%2)cHIO^8QtQwH-Ao2%VmQ6Y0#AWpn_@Pt^aKZ+nDzLODZ4^<@*qxBS+~doRq#7$%)vfEsLpL__bW1~^g>L#Dm`Dm'
    'rcrgL?WTWCI+~}y|Ap0&#poclgti1k>9`%qD0*l`_FMA^!J@Ds5Sg|#FWQ9HeU7)74FQp7yZC<n)d6vsnykjFTskl$stv*|bHU|~0o5}LW&Z;c<{+iqf_>x3s|wFHKSB-Bo;e-6FNMxjPU6B<pAzb#Q1cH|N;+RENG-gZ$t9v_%!M$cO-O{HDKns+jAKnSP9W)~r-*q0f7rU_Ja$Fi$b${LQjbJ+D!^f|Ap8Sek-7z&%DF1E'
    'MKd0o81PyCtgg;q1GGn9<{AC8M|H3ZKM>Vo{A-y><}kQbD7~hNQO|q@0{XQ!8cYTR_@3t{I>T1BbfS+T4afVrk?GMV2yGA+nN(xe)m(^1X-lD(lgIbd1jU?)S8-%0GLd!kI@YV61)-<b8aIn4uCyue0&&0kYr6q}fvhF{QH)S@nU-Frm&D=Rl99sr>LY_4iT=-EWL>JWi~V!zQQW}Z72{%A!6LXGJ?^#0ze+61G-80>()dNA'
    'T;sO)8jzY~Iz>B|&e$q#TBJ>*Slj%HE2h?En3RyNKS-%EX^MCp&qMqjd)kSBVk3rOrVB&nn>1;6@bc1wDx+TI-;$y14^(w(>;NZ5b*s8UeH1?pL=ZfVvlXG4e0m&&t!v_9nJT*aqJ`JM)0JGn+LW<z+qW?4OcX$Nxxg6%F8*oS9k60ShmivbJ=^luWEhjBYk{&lI-S$GgKR=ah0m${Ek6z7WAKbAIrY8p$wN1gE5Ct#Vj%T~'
    '6i0;$WXJ4VyxMd<?aXPwFj~iMy=yKM7O2{zHEz9m(TYMiRmjHFl}^P26{XbKoV%Gk>Wdpgt?}+V<v9;76Q>a+aSOf_yqge5<m@giF9OLwajnShf!36zO@C%N0PtP@S}NcCQ*}HZ9F^TfZEPd*=Nrt9-~^7=(<b;owOIQt0o+ez)z-Jj(PpT+XEn!&gQlQf+9o`qLFBwHqB<wj)?Yjvb3UQgk+&^-=ql2s=^(7m(Lg46JX!4p'
    'Bs~jmIcWJ-?}$MhpMd1>mDH_UuoW2KQ`cRI23;`sMan|7`MLQl0^?`ZLJ1X^uCBRwr~nd4`AL*69#bqHOE6Mc*<O5THjogh;WjIun@dMfpJ*c@LEjd!LbrN%(50f~yF+pASC|hMUVcR-d{pXbdFF@56&St6PYNbagUMT7t)1YR-CMa=B>4=bAlX$1tBY!+67|Hka!s}brc;EYt%|d^?Xj*q&8ZtW-YqX;0;h__C{M*J9sPaf'
    'S%Z?b=@A**{S4*FhX#x2LGpN^uJfS1B|xr$Z+*+-DD#MEx%l1hxh#fBc!~zn6-M*y9Bgoh_FYXQPkRo2aNZc?Z~2wcpmx4J^=93rF|VfAd%NB4k&1ItfyO%V5^;u-T%lmru`#3Mi^Wm1GgBmCl?irAw{W7k^j)h?(-CO*4K5N+NI8I4pg1>;oVRFsagP-<laZ_E7i$Ar#5bj??0O|MW^B~MRiG1pm{Q{%gMVBxg~|3u@0<|b'
    'wp$rxTgP#eFbO3VdtC}eOUAu$vi$if5r%hULmgG4<e>EI2$ufcu%7cs5_176*_KGyU8Vz73isM<as2WB90VgNYrB?S1BplDFTF!@L9?bx68Z4WGs6XKAUEYw4h%BSW=eJWkS9nT9|NZ}jUe!<XGiYGGasdMH{lBfJ;qX9&Axvse1)Vs`37DDL5?sWpnUmR_dOEuxI1PEYS)NcJM?3u;fy3{N?(cxs805$X~`R+KKd8_Xvofe'
    '<3qUX;zZF*HYecA*}5Nez|H^DEsdFt4q<3?C{-Mj<*@-;1q+>Y+x}+He<&Dy_vh4u4R1^7Z~dNEw95-ajPh6AfMFaaZuV{C6<{M~z)Nu$3;wuH&<p~AJS!?pH(6j|>mMl$5Z12<Tc#HZk|WcZ@Yk&S?HH@&;2I>45DWfyHtk7_V0+N?ADi4oYy8Y^&B)#0j$%+Yh$8n?LRNwg=QL?7H1%ow_{Y#E6-}FLA0W%#O+}7%+zPEu'
    'Bf!^dp6q;_(Q>?-C{wZ1WyViS-YS$|1u6*GeCAoN9>eYzpJRbE3WTTKMkNh1Yo|8NtF(p3i5<3TGkRR9AO!_^T`Q<mq)$<$<OPoB+}9=2@sgM&G+uWG=}d1gXxP-9;kHWS91S}O2SRr^@gqg{=36C913hKY%=^^Wpxp+xq}tGGh{qU*;$;<=4tTkDoczH)h6-t^pAPbo#a`Dn-V>?$>={WKRwiwzY#K{gU;Kml-&jmK{-FT7'
    'V+Tk~7{7EFK*~{T4w){Zkx={A!nR6#kwiVaTr`+Kyis^HQtPpp6)WM4|E{>V?##V-$4Vpdk3tUpa{vsVEZ4o>*8XM1|F0o_J@D%M6l0#T&F_j@cpGeLX%iSY(Bo_%#_ZhClx>bsLJy$=uUZn1kj;PRVV5_j{1x!cgGT$#XlwzLlttKt)IGd7JQvuw#_<oJ`TXLw_T^1}e|R=Xp(zh`wj4V90!Nb4fm>qg4ok5Y^<fZV`c}5$'
    '_9B07vmz$vEO<*qIEKDzVSw|q6yYaFoaQj1_*O>3Ga<QUK-v$w<#&!Jy2r@)Hp2G-i6*6#?uVdA#Q9FX%FN^qpQ?c#^l!cZ`i5p1C8S*+bOpM`6R!YzZBY1cPFA}i8ew5#ukHFD<?vWrmr;P<G;640WZMbt>44oL?)>sPue+E44UFYarMG#fh;t}VgHDb`LoKGPNg0LUAVaG;_Un^-Zn>k*wJBK?$+<EiBZh;x8Gmm%#4asr'
    '7faDxh@8_|yN+f5#5~)v{}uu~?2UD7wGEN5t>~&(gy7C^8tYLEJw)ZuqKwAXw98_3Bp}K&vmE%XzHAG0Eyf0~XrL`8#)V`FBogL?=P-9=l?R9#Q2Ne?GX*84ET2WT;ArvMKqM_6@4|rniFCL3riWmp!E@p76YV|p9u@UaDy#^q4=cH6=W#8LJE3oA?!Hk?+F7j66emv2fqyKd9HUAbUL)zwoIgjR=e@y;mWRQ0H2=P+@0+37'
    'Zoa#;CfI^3qcI0W$yw*EE;1%c0L{wAOsbLFy>z#%xHJc?y<;u&vFg))rEarYeKV_06sJ17!b)((Wqwzw_E!<)DS2qy(pt+GRU%p&WM6y4=`!7ZH|rp!(i<qon=>W1=d808=1iimb_Q&jqsQ9%G3K_~eP;Y1ONz1ZU86;LNZ=0NVf@x<@KwqpKxO|Nm=4`UlB7-oQkhN6%*Ns8?>l1(Z2i`S`EDm90_ChRMFg`qF4W2=7SCaM'
    'pJh?f$)B-xb>i~tcx2q5ke+MJ?gvg`t0EfNXrxc2`%ay?V<*HX)*lm|5VAq7>ah}@Sk_XDz`z!wQZ;`ruqq0m!k;zLacTG=?*#{$u}a9&?t!^Q@DY2NjH!0T#hq$&YRVD?w(xt}J4w{cy4W*FhG00*Xx8lGDG|hKPB(F@cmh{#a&Ne(CNgoAi>|hWBn4cpnjt#)XEqEUaf>qSKh}$T&p$o{kQqlc&NU&is=~lhcXo4+A~eqj'
    'M%`On<+tbYzFROhsdySzYaVAlL*<UDG&+cnq(Udc*`qEQ(q-G-Vz%>2`iX>=VJqJ;se`fkWZ{}B5vUSQgY8<X))(jyO;%mS;CR&Tz=~Ctuy?)W&!g&uJKpyNH2;Ic>bCa!HQr%YXXKl*Qj_F+H&G^Ys2W*oz)UIFlC8?g3R7B-1t}(e@~YaSYKQtXf@krzY0;ErD8fUiwO1ulgmRae_ne$yar1beYM|59zX%hF?w5L;)+*q1'
    '=lmRCQv&*3+Y~tsf{}X9&pN-c9GO{_ROuGkjT!IG<wZBLIxQ&TgG;kcgd-{=`78~bl|nNvmb_;xXpR^ghl=B3Y1+YK+oEk)xl|E~*p!Opd1lXA<Mp)Mw>O<rPa(7^O=fF4zkT{`j<Zn!a5hdvjn-*)PZ5eG_D83mWwq?R8_;pi|8lt;m0G+;*7EUk(8)XCqoTB-BOHHkI^N3jF88^^K#Ofs?}HJ?#vvAw5A365WO*HE^fHe_'
    '3YZSJI3F^ls7g{KUZ_v`1No8dSUWU+h1B|R&{-L>ps_p3?tKcYIsxPi$_0Se<_5m%&|A%@Yu`<;@SuIV=b<Zw_|3IJH=@%rkY&PTmMKc7`LB!O3+3+~dtEHBc@g8^xH}0@SqUCE*>>*wL?tYXQd2RdvcE3q2*zL!5cXbkkXGS;qs>G+JSxH!Oq{>5N7Q5jf@}5?*)#^ByOjJM^3X7jqpf5DEumCWdy}pczAh*CnW>_c)&#b6'
    'QKVHv=WiZYe)&7u#_2DpY1c@5(iCJJE~#l0tJLtT*?V#ukShvVm&%y-2x9nuEqu7AaNr6SGZN-7#YlLe*tHk4KN(%$TVoc4RN4H2%uslO%kQnd)!pW7%^*`ZQvaUlh}|1o;~VRDTzQg-`+t1;jH{%DFHD|LxdPE4DipJ9N!4oD+y#!}gme%Y>bjGLbh3t@JMq#`t2mNNJF<X_WjBm!jE00xo<Gp%D<`BP0fWOutJ-B5`gHKq'
    'qxy(1dn+3bLpm>?Nmb_0V~aOC>Js>?rXg{2h-`s0IZKeML)u2`>vNtyh4M4IsUxWdMS&vl1aL-x${kUC3mIXYF4e)zb+|x3GO?VjWRVClcHcw*Yx+zF6%p~|z)l_gik%jL%#hg~^PB<sbX5jql03O`=spq*XO>-l0YzVb3z^OC9$Rnx8+c}sk5XE2HC?G?F6(c?HY-J*+xR@Wx8MjJMg2Gal^0Rq(n2nEF{>T)VtFGdQ64`4'
    'RBO!GyccI*GV1y61f_cm-D{HR?JgsJtY^8$s%|tuiQvADLoGko1VC($DS?#SPNzbCNDTR-*hu{AJc1h9yKh9TVknSYAqh(%J=}rJ2~9pJzS_cQ++G%b!SH8?vfzyBlq>M*TT;S{S!H)m9JOk+U)vI^-FRvS2C&R;<ow-rWt(D&GU|N*n<fbc0Y1R~wbVvau-jp@Z2mIMMTGrC1m_v1fsl?e@t&MIvA`a5?46U`bvHWZd%jL5'
    ')JN`v$Vu|mvj$+py9$B1#c=pc51cpGQ#Wlx+LTB7xgxs|kWW*Wl^<cckTZETQqj**=UXd_cN;N$-s#C+t2-x$cOddgbt@dso|LH5<DyV8dCa<eijbDMDvT{&KOkCwH?n@buWM4~3iM2ND%V!m`c|+u1MsZQt`8uC62Nq^?8>}YANRz!0vWQ_xKuqZ#~I9DoZmkGTRSF@3$_NFF#S0+WsctrdfA)+khXmgiILA-d&9uGC3I*T'
    '=aVN#J~$Ta)dSX&6r1KOiNjT+1J_>0=4GCvWNmI-PLy(JUIsM2z>OhR_zJAPh+t0)nnYX2-4anp$r7iDXs`qK)S>ZKBaThv*!3t98@c3(n75s7uPMS)xvK7g&goBb6k^nj`O`5v%tsK$pi0*~fEwbSW^K*cN=${w#E$S=(smVZL%Eh^i2QbE#M!tkCBGe<(Miq<fO^T{`6gTXvuK4r%mGi_ykP+bm@J$^<auc9)r*u>*8^$D'
    'OO*Y_h>djekV`w#Ggr&xoFnw*b`a-9`J`Yf&2VpE$om68?I?~j1OrqBeu!k9oeM!t-Z<;v0P_t0q3ZB_y7{mX1h;Y`>Y7<uXs@2h*J`=E0-(T46{yO*uOjVjy{yDxpD<i<D!vot_3S)3;Sgs84GfLyR!zORRnl#oc<GdLN>;dEbNcW<1QHLl=VEh@Tijzob`MudQxrI8?K@;_=_ye6x&Z>W#yu=Kz~Rl+ts+rsN}eYfQ@jbp'
    'gtK&!!eDx*ih0I((fQWq?>#a}E)A!p#r6mxRq+dA$fBe9M7<^Tu&^fQqSzo3UFRd4TG6Zbp-~U%0mnKkK{DT9WD|e7`Id%tlU(>?H_J$`FD-r@zOA6?3G>8j8EhDCwC-=#QDQjc{e<>4jtcEj+1=Ll>rQL4$~G%W7f|as7$%tA4oiVjln!d;k-^`e*m*r$;+aZeg-zH}>oVgZ@my2{(&+vEjyzV6UxwMtvy;i~?DZix9gL{P'
    '<+XGd^Mq5q0)tY4_6GR0e+Fw4?gm=*d3EG&JdF1DuIJ^G?*V?89`rzHzcXN>6>lsk#tk}y13vE%hZ@lso*crSYS}cWMnRNd>s;mjaN37Gh;3jw3CNx9Z5OV>|McJTv*s9}j|f41V97Y`^102EG4{TxHmo=Mv>|cDu(hzLH`5JMkaz{nF#6IKq;!zx$rGc5Q@7){--g({37N3PU3isM-CO^r(e`Nr_?#OEme-p}e1ndPcy8_%'
    '|E!9L5|Ar~w<wrj-4HdB@;hiQ!rFJ2&8j<6X4~ro2z@}`D<B_kcpTf+4-4<Tc$4>^*BQE)D{J_RLyQXwM;0hGtb<>H5Ek7mN$!*~)_kf41I^W5ywGlpxwfTmqSQ%|Y$F1MK1S{M0xa+}IxByr+U4=uSGywCO#WaX)^dN{1njX?ygq`|UBj6SvgvOK!q1c9M+uC1KcUGpLLQ5oGLmQ8!zxGE)0`GarWRp8Su-&(`Snn)Yolju'
    'GxepK&rz-=F>2@~rTtDTbm$(k^_KEw(Ls5uY_|P6oy75oH;dXCCOd=SNdv`ykDj#@vwP7s?)^Z+QR>+k!gZjJ!h`{TH{uPR-?Qdhda<e=e$A~_XjG}tk=S&3=HGdl-~o~GlmV~rF(w%KpRlFoDkk}Ab4(54-&yH&VMq1BhrH$~Jbh%aB?EDGE><+7Po1=>d@SUJ0UnR_BqHNYBFFd}xFD&W(;PqRqVnl=(=<8Hhd$zs?6we>'
    '+Za;f%NC7W2sSzB3sBy%nXtY8-DGy9<pY?kmJlxUm9cGTtGaGz{>P;zh|<(p%^(6vo0`a@SMXy;5OG@)Dd&r9hTU6T=Tm#e$m9WQB0mUYZH;n^h)4%jAIkrQ1os3C&_%=+esnuiS#A^CB}YB%6zDCwl|KnC!0PaL3PF6TK(w#Y+Bf6vm`Hy58Mm}lf=DU62WWVOgK5DfCv^hBczg>btfPeP94>#sYj{E#)@_yQ7R6@=(9<rU'
    'ru0wDhG3~^3E3%Cy;q966dt%aebIR#<Ur;@0|dLr1#9<kKwJ{}4}N8}dFnoJYJf8%tec>>2$(@V_v(*rK=CrC_LDt|?`E~HG+q|m$3q0?Yg9R+l>_Hlh3R8RxO-*~5_DN+XW3~rZJeFH<VIrNakao``z#Yrt0G%U19-S#ea8Nphcu$DF-bPM*xU$0;qLne<lNU@cEqVV(TrqhJr^we<A9lIpf~Afq>_P;@{}wU$KZ@W#X{aJ'
    'I9w4Dgo=S(8txO@n_3DjWa;%U;q=cS`X4rjJsJ2`l1R}^P_r*7qUWzX5-)f<2~)cvtCqc*qfnGCgxw6%P;3Q2ZF4(p!;zjrIo`f?Yl#g(qrtB8z&z-+=UsPabjS{Q__p@2a#@jIQ&1Kt^8&8JoTojD$FOBW*hwl#wUsRD#atufCMdm!dgwX})8d)1>52~EuH>O-4i}Di;*U39P;-08%L*i`-YqeKx|A_%6XlRX*hDlUuj%Bj'
    't&}Q;e<iCl2Qcn9x<>NbQqd;dthujgWcc4-C}>#GvsO1_IHf^q$?Y=F)|I}`lzF;`S@@dR^9(*_vQ!PvI@vb=LwHpUybT<o6U*9+7i5la=ebr;d4s<^OV!Alv8J#$Oetuw%h?Hcc>U_6l9y%YumjDi>lV|r#Vg91^Cw#9K?F&cTnQ5T-#z)9qKx8RxWL-(JjK$a6lVeW?PwrwJAU!{w8X0B${hvMQIWi$%ANOIsBE;!OKTml'
    'exEa}{GC=y32GQW2_mC|yrb<XB7~)S{aXcWEwbUPpMX1byK2smK~P9MIDZf8eH1JM_)N-K(_rO7Lfl*BY2R5~n%*KD=tBH+*i2epD|@RDXV5R^fM-<Q_tb~O1HLlr@fopPLt=Kr0|XtPi8^k-x+3|B4KpK_axl=Oj(8u|cI`*Ht@ZC%`HM{EnV)<HS!msX?NI~-CHAFOQCLuX*zrIjJt)f?X~SXn0M9xrv=7Uof)pX~A^<TV'
    'AtnS)$%?wwmsz5$+xBQge&d@6nBim9BtzB{0Ayi!`y*Ti4PN86BXr?A!X+;_=Ja>%()>iJ7(o~`to*ZKZ;TRFUloiyt5ReSLAfiHk?SiX2t29VXdMq?Kp)v(uC92)xR0tq9o1!HeY4a1Lrxz57e<HbI`9Ow6>Vgrmm_|Eksz*Q39J0inNoxr+Iw>;^vjH5$WgZ1%6+>Df79PA&AyFI-RZWJv3Ep(2lB@9EIOPw3$ke`(yXHS'
    '^^~qEJR5UYDRM4cGLO>-Py0P=OEi1cm{(x<m(Dx**Ciuh)A;pPgAh%0FyfA>zzq(o;M5q#t$BEaYsWCl`9kS-uZt*=Z3&l&#W%RxrQoPX$AtgICrPN#=YV}%D9eLQA182z%3!X&pj!B#z^~{>q>Ah_U8g0)AZD)XRmFafuZTf}aSr&Ruhrb$b~b8W_8Fo~GTcg0Xn5!wyCCyZfE>-xwEuRJ@$GFcE`+?#`?vcZ%uWN;tFe$4'
    'jT<6YN%{+*H%%-jHAX8^XKxUcP&4D_Oh!AK(ocdhrcCph@CpgmcJ$vC$Jv}$$V`42890|o(+j>Iv?T1h)5_qRYp3EL9pR{5Am#KLl&|`8;?XYw^jUc_Sb!<~Y0&iKqSoo1@kgo5IwrGWP0n2*$ROHlgJ>%UgRh4^R1xe0#dZ~WwSj|g{l}pB+TySuuBR51mb`aS*3XaG@&ia-Qz2<KGCd`R2D#NR;X}1@-}fLwv1pLoC!7-d'
    'al9bu!Zt)m#&}mz&LM!xOS(xzUUuK@7&V-~UtkM&&MqyJA;yD1|B2D$!!pQHl|KzSt<l;@x%Vn^JC|RAhwJ0%w<3Y7Cd;<tpRUq|%ME>)sn*cC;r(LJd-(iPCmtu@BXR9I-e7^QUZT$4&T$zwE|CV)2C3@VV_f}EU!dn`P$t!%k{TMo3%hzE5a{FXJVdXjF6LDaT*Ji=1&B_$iK*z;`;`^#F^Hv!+o4t(%;r{-_pbk#?trso'
    'tLtQ*NWN@7${-Irvd6L`A(dF_`$XOOFOGY>XRtUDa_n;W2kP_P@N(*Fw!vI=*G<8XQe)++>8Eqnb0^^CAM0lA)Qsf2VcUJG3|ZB2qlN?+LfD+9V^EA=;zBHPP}hJ8D!x7UrgQ6~N@$Ck!2UCHZdS~H#V?^V_+W(JMJUr7kQm{ItuN*TQ-#@A^1ZqKQDOFer*q5?m`#6x^8nSOm#!w|RM}C@gte(zHC&S;8T-G{Nm9kmZzdlU'
    'BHP|kwE_GYO^y`3EEO3#p~k4^K)P?kN&3$==dV4lmJ$MVe!u&bDm?{i)+=tE8V2Bn2V+gAx!xz|yTfUJ?)bEdx8E`EEJkRiBlh(=W=PYt8S&*kBhOc#LZb~SjC+J^$E!mTxkynx5Zd$s0b|P6fOw%zV-fIBy(I3*jGsAJE?dzE9>+b#_$SPwV!+hIcueoVIArD`w96FGY-_lQ&k`;=!+pzeI*I1V*7a2$3N}xne~sPS_iLKw'
    'D=ZH_vWJPu4qR~sW^K``Zu+)&h?jXs&uEq|+q?@%GD70ZtxxINwA(sl5TIWS=m0BK&Z0OH4^lWMw<oI`-8GR!G0n}%FmRZo?3)_Q&*(C}n;Z2M$FThpL*FakZITcnv4C_QENal5`P`E>Fm$w9IFda~#j<djs@G#*CiXd({}FN_30FYrj7c?T@hs41;0Rp-KIn?ES{dhM)JtvPvvQx1+V*xoH~sB4!u6Qh#85o?W08|aMP4cQ'
    '<Mqg)#!w{zkPCM0FM1!icC+~HeX$C|MHpv=Vj&5#J5I6~7}VVQu_tIdmHF{wF0+K~3VgY=@Y#r<YHLW)8IAdCPc-JW^OKcGg`5zdd(|Wmz=VthKtY{FrsTg5Gq6>!jCw~)n*!#crZO3SWJAH_zlk8TTe9C+cgXVniPM1jHK^~{HjuF~91Q=CZNaaCzrN%ZP@fRcIbrs;K-!<m<z1%?D9*vGq#9*I1>(I%!)f^5jTSZPHw1{H'
    'ysy?O*a7|jNsUUfsqlU)cL6~|mLo1rVR$7ql2%o+V~>hu`h+$aEYs++Fg~GjmnmTwZJhnSTJ*=@@R4MN2I+8mCq(fW?Nc>m{xr=Yk|QL)OOdFK-OvKmr5gE%lM$Ya<;<b<;BF0ik$|4kzhy#}FYAX4#*HpDzY6z2b4e_Sr@;ln65AY$W{i>cKX}b8hwgHePWi!jp~=0Lfn+*qHj;EC4~e;2FWZ67Wv@5T=^DgSqu5dV0%%@3'
    'O(@CH4PHH(f|ETOcd!kaH8Vx(L@)b2z8;FV0W49i8zS~49Ty;=2+fQJyZ(GmsQAZlfTy*O>63*gqCJndY^Ln!-OgJ>H&vCWAgIj_l|8Z4oB9%}{ipuP$I{ylqS{aDeC8yY*2C!YgYpzL1E{oaqgfW(B*$&@tsrJZG$AkoCfwIXK!^I^dblXrwUUH^!v9KxR|`5@t5sFwe~l*gQRuKY>-gosU|&dxc{TS?N0+Y_J=>y4LGXSQ'
    'TUe$RF|HQ6C8pakKDQs|hUb-IG0SArhSZB>UnUpY`{^Xwd!CFYC5hi@&@>6MxS`DAa4<~6D@nyt2khdljFoKI0*vsLXOAZgSobl_X|&wB5!Jc3z|k@Yv6D2Hu8@mJ(iSXTW1nlt9ZoQXN+cxgbzahpQ1Jm@g?fGBl}1{SuwCrJPO?hw(41Eg;um=|I*<|eB(cYz$e$A@r|lp}3P1gi%EKX;|8A-#KO1UNYs62-NrF`USV9{Z'
    'F*GzAP3*<p&!PL}kH8#BGtV7@jETH4GIo<u=`!O*;uPI$$FR9ivg%~=gRIoBGDb`unM`)nu`@f*jSeW-;x*Uvx|ete*OdKUWq(DAM~Zyi;19cl@TFD&Vn()|oEUIY_=lHVaq0t&#QH>S2qDcAajD9P<zVSU5L8m5D6t^V6Kump?PI%(RDqq0?01zrSS3CAvm%94HUzA}hcWW?9t_xq3q`6chz3fC)w=XOsYEbh3`B5Dh^;{x'
    'g?F@<bmd{S`jiBJ-brV0B8sEFfEUrI6>t<5<w&x)$HM$mUz@@71It&8+Q^P8%gD6VRou`g@!03D0VARP99ZLf3x~aS?X4+<4cP>XRwiLJ7}q6v7bYOsZzMT)^!f5aB+#86NuuU;+`MUV8@&&t{6!BI<m04+=G*ret)r-?*>4D1q<Nfx2ZZ>siwy=xf#y}@8~4zqR?{Ag%1Dz>VV+2cUL?bE+~j`S{pKWv2?hV!O_`Jag5qtj'
    '=>DLrYFR3+wiY1o(wvjPBsWilPM<Cy?{V4nAZGB8_JWjl$*ps%+U&-GgrdqMGtO+8_d$6-L{2+})~KHyAdl3cBk}rNgt(?4|JAM-7u7l=^{T${t67)KFDh+N2;)vQ(K^5rHG`u(5r#Rg8|^#Fb1N}*E>=6Fqo{J8_?GG1OfTstjKwQ)NAR|<>kKk(!}n?W#r!8+;Cf8p2nepXjM=U1KGL?wdYtzHXaAXx5Et=@u>(IOh8)C4'
    'eo@OuMn~g-B9Xg%3Ownw>o~qcVtMComb!u~lb~x<M|xR+OYi;Q?=)mqHOEtjGx%<WJ#Tz=D~$VE;(Qy0J=!C)F816QV85%>gumV|fdAFgmgowSI1&)d$_f_<(6MHs?Ufx==fQ?P(g4WzD79K@=U1_RWOjNTNRDyY>X0d(9ZRg`GOQpS1G?hc5i~=zeW&A4$l*wcIIH|kJ4H2sVJ5qE^2&Q0<w>~;Sa`Oe=MiOxs%{S)tx?;D'
    'yP=5$uTnodFD6qRN)%KuuxbI6DC1ga!>AG@uz?7)bdyHDs4HcPwSne#$Z+oo61y|31$N@vs!U!5T<}9CD(EEL@MOI@pLb+Z>;G|*<pCPdiZ;lg_o~NbXMg*_{YwjEklv=m&LH<Cwg0#ka5!j)mZ#OJrO2y_%IOC=x1UOva4^8)`}#e5sHT>tDPonw6oErTeQyRPBcNetT`3ABEXD|;MOm+)5xM04G!Y(yWmBpsnQ8?a9Hu2p'
    'B+Jsh>AV<S484=nwZ@@~ObQo)^Frf;S0s4PaeqmN%7!2&nO5w@9V(k=we>>0=qHT4HrR^;3bD<Yi9-mI$~gR)M4-@o-pXKI@$f_P=Ai<16@mAr%i}b$F6=I?NGXc?FIn!)TelHj!#7K7M~RtNh32=T;gaB!pAUu?#8oBVBW8Ckrwx5N@UCBnMUWod$3u4hhB^w{4J;>Ies0(Q{bm&WXzpfOGZ{Tl;SZquaw)8UmV>hSh%ipU'
    'RYYb-AgIj$7LUUD`)#E`Yc3fBFbdjK-VLHAyF+C>@HPQ6Lm`-x7O{;>`XhRqCpXt{NZMsfO>*8XY$$bcp<+~k)}nTPY15r(*j%D+Fg*|VU@2XxFBUpqpAL$n)IA@Q+nBZuk_QShq-Oq}z?$ZhNo`sw%3lnBg8etXiI#$CFK?+IOyW{>20tza37*ZG_qD59!8Xcu@v&46qK{n-PP`CwM73>~YPxRY39uXH+~0ODNQ)qdM9;51'
    'cqQe1RxBc!?M2bboIdJd%yL_{M~e$r!6w<W8R!sL3>l~XS29roj&$Vl>ZBzI0`<XS`CfIV65#P^pRQ#!`Vd;nq+hkb3C1wTZXGT+4wtY(_0-YW4RmT^YD1w_cwfH#TR)rJXk^q!kjd}k+=IObHQKEm#M0mlX4&|P2v(1Qjs~*3=DIdygZM0maOYbMY58Cxu^>VVi=UJ}B%ZknfE^L6%gAQjgp>KGJ%R3o>q-r55Ibkauko$A'
    'F;Uy&;`@BE*du|%0B4pqlM1~^k3hlAzp;rn_X%BAVXzn_q=W+k-j_L!9TYgHOQn1!EI2~=ln7UW<zmNe2U=Q1;}N9o;A!AC`&YD5ZDWN1lJ)KNK@{-mZNKMFR6jTBm^oYrS<{~S5ml5(cJm@{aA~>cCA)fLrleXdNoD4bX8u%hUN|_k;4`MK4jRJk<@~(Eb_wpKID)a!^&YN3nkK+Rsu&!5R;m39E?CCDbrwo&{%C~@V~)Tx'
    'Q~eNnHj#@*4kFFnvC6$ewt$86&bll&o6k3sPT2yFkI$SeZqH~&{sO|vS#Cy#k{Y?Ki3r9*)I7>Rkh`JQ{g{<3jioK#HjaF`{=&sLnVDSw+q~`A$NSmyp68f~PodELf7a~33BvbO^aZ9XK(-~Kf~8<S(axH)Irt|%(dxG5;l+^q%GJcAD-t~X3SjL;b{R@*ijA7Z`$h;N2yib>&KiM;UU6A8Wx)@t$0t4xg6N#3<_J=9f^r!i'
    '_-n8g$Sz!-lmV3Rkg0ov;@gMK1L=HX{!>~NC+F#=i>oXc65v<j?cdM+V9^&2E!F_)K~ECxG*>PdHf00l(TsPU;4>_TcI>rX`%Od}s!P;^gF)L4;33*Ma#HGEq*%|Hrs}1;kG^lG#w996xa<|XCV}MY!kD`hKH`;nZ3=bR=?FAM`_GQ%?P2DJkU_0nkOU(=5mK3_&L+foK}9yIYh6=$shP)>2d~si)J|-&R1INDCJ`$y<I+r|'
    '(#(Qkm0h%`7C<QLeP(RMPO_I-Kx0`OfFQ<w#85xWn6@t~zVMC&?vpCHxsK;!E>jYb#YuYz)H2P?5rXz}6<J(6QtLkmVq1p?B%L`K1m;8XL??dc=CI<-oWku#@3sUuS!%1a%ngh1i0HRPE6cPz+e1@CFo4xIdx3w0Do5YycB$T9<{_K*DSExZUBj{EPxt@lthv#xe83i*r>7>Z<K~@&%St$gVOC+<D;Vl9t^#JD>vfQncU2W5'
    'lpt+75tPxNb&r7%b|5hPGSY+Yk<YGcOlq)LTw0rdIs5<k6PF}gZm>wwh?DH3u$NW41bTHSaJO|8&MNn@Ut;=#od~Dg(nw&iK|dO2b1*=}>25*l1P``l7Z_Tk30QdHMmb$K_}y>`^$qHu8vB@y^3!75aV+Ech;EPhZxr9k(AN=WukE{T>$bbX7?8=*l-|4Hx9EgOSb<{#lf*)%Qp4JrHw6PJwS{Lz#G201DyoUsu1-dOlt0nG'
    '!J2ePbS}re2h7V*ptB{hrCbbXKZM~(xF>?B#?(AH$NwNz*@=o_WCWOEaK{I_LN0P{0CZ^aO2jYxkI!)_A_(#?UI=bldP%8-o&oMeJmmnxgDdV42R7<n7E~e~)H(Nu?EDdVm%`*o((e6=>qn$Of`{fmCS~%zPmQx|SmCKO;^+E=_$+f>$?)O-FW%DwK^eW8AUR|gIxgedPp!s3HAFj)G<2O;9r79J2E`$#G#n64O)<vWTPSXW'
    '7AX-bWhalxyjdDqgz9fA=5^8K-ibT6dT4fN(R}v=!%(m7&*ws&zomAD1<P8vYaPVD!n1ZNh=C_M;Z>gM<ABx#OZ>G%#aeXJBe6knU;&={Z+An;FivD!j)kIZC;0@1uKq~Us8?W<<5A>>fh?V|2d#Amyoei(q=3B?P4e>1uC*3CyKCV{+rdUOko!}X33^S5jBZm>YY?(1Dv4a{BCL}5Vj+q>_wjO^bCVfD1^O7o`I$E@xwZga'
    'sLoR81Ap7fLs1{K-gS6f{#_uM)V-bD?A21(#CtR$tRss<nq<<_*rBZjO{zQTx*6P%0wqTAPEUxar)%v$Cr}76Yo!?=4q!qb?`z{EcYS!oe&UrOaJ?)RW3XvxmJ93z3k?I+KhPYvZq2LuOX}>Ut+)(wW3FcM9+&Jnv1>h@32rNu%Gn%tXP$2lOLmwPfV(u}S;_va0QF5u7G5Qk4mX|PZZcu&y)HPKxq7*a{3(~Hn4J<Q&c=ND'
    'DFr}9Ei?WtUvY;5rVtkWWzAhkLt=gBvIyLZM`C4$KSUV;VU<Bq{O>8#TK20eo6ieYPk2UMr5@cuGW>=+ZKvpHhT=UYrEygbrhleH=7?HAV%z;L=ofiT1o2az4_xnQVn{{ZWVubv6*R5>9uGd;ZbJ80!S)!&_XG8H?4pV8RNHbLJq>nMB>+Sja)`+hvoqK`)CPHvsi7`~F?7yPWm7v97u_J$XvLM&>@{Y#Dby&`<9K_(%`EFn'
    'kS=}&HY9y>{bmE{l0;(UfW@Vxnp)}efIqYtJ}y;~k!e6u*I%37?`wU2$gvfg=fF8!ggr{}GG#^ylq;Crr@P<cc>vHbM7m@H)L^b@Y4h@C(Uw-9=mB+wm&fzDdJG_ySjgT!ycn>!&tvOd#L&tV`EW1N`RYZc<Iz%xP&8K~iaY05EBsmXCIsL75f7XI;<%bpY~XCg*3KM9F{lo6GezDr(xVWDtmpx(?>$nv$%Pjr8#{^0ehadi'
    'DZuE})rEct-T9UCY<XToTUMHSyabdib?B<WNNMWwt9AciOeDXutOXmf`APX65{g@2Rl!vesUqktV2_qS-Eip=Ls+P>KA+2sv9Wb7wD^M_B!2>B9U^nf;)KXC?^oiWec{F>X+hc=ZZ8WYQQO(EQg`z|?J=|`fM{Lt(gl%nzEL$bL}QB^9=jmyjCFTOp3w?+j&ChT#$INVV2*az0Vf$>Z^9ztx$rVU;p!jdxNV>mkup$}_WMkc'
    '8S$Pob|BQ5rV8o96Jj<!^@F3aj7?H#6Mqz`p<~wkALT}BaO`28(9q*aAs07Fl0-2Q_URbMp*}a4VKKU6^kd>a$3~ry!2JdDq1h6Sq^fCzU*nOoh$oqE3g9y1h~2+BIlx>v>|fIEoy3*B0By~mGiZ0&<;UE~b@B~87+tcwM*njDcEkJliDi$aN$RFcL=JExU%Cz^S}Ik4T4d{G_y5>|Ulroo2uJ1WCnwBUsTZVP4Nc)$tBXSd'
    '8M{M7cJ+YmEoQe3C2}%hrlsSWrLfJY2U^Yb00000U$v31_GNz(00G+H0*;Ii%Qt>KvBYQl0ssI200dcD'
)
