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
    'e#<4uB#FnPyVlpr<p|Q3f^?F;%IKe_KXX@GHR-`WDEW2ZBxj)lKCOqcYw;|6JR!|GwcOk@XPA%}NVc@&pt+?cecK$b+0WFp`YWKfNq9XfV|Xz+V?GTV?vGj>)=qS-dXz6Oh+M7T^2c=fPD)b_Fs^i?ELRhl!GMYD?Te_kYn+3<mi)7U;2sjc6inZy?D@yrANLg=T*tBO=^u^U6Z{dnt?4H~(Kqtl8<%3sBT7P0w^_dMOx37}BoWYLDLtXT6@VPGD^2^h%U~H3vf^-Mt2vO0X0W?Fdd;Dh`fU2du=gua?lo37DJaR@$6V?YPR-KS%P^FLVdDKO&R$j>yT+}cu2&O~h|L_++qQ!&rQ7}+;^KPVzSaXGlxC&7v`XxhaDmRtQT-)<=sh>0TZ5sciZHjL?nhg0(7f&w50&=RIMF+0@kFRnQm(vImqS>&gfC`mHieE7A|biBi#BfN2CS-T$4-l7o_lN#C+PtB+B6-o&vrP4GWOtcc*{nzy}Xp^JMkkCsZnJ&>FxgFhM`gnLay&6+kUA>D>BGH?Cj}sy>+PcbepShS0?O=G9x7Mdy$6Qc*t}b5Dbg!FNV%yEFl~Tpp47U$rHC`dN>zwaKC#L?pEs?8ItKZ=8E;kjjvJBvv;GRd9ZB_?bf^aS0<TX>(!vv`&mYEt0*b(meI9TMq2lB3e%-2>>4I0AXN*FEbkGda==pbh0Xs+ETzXTN)HKIy9u+;xI6EA@z8>l8R!H&rz-uSIs`zVLPmZq+q^Dnt3Qi)qz;&D>mlZs6a-E>5PAX-&I@xmS6I!rYQ#`($FMk;`Q!SxoJSr=Hb{TF+}0Vp-;^tAi3TGm$W|s$2=Auq8ukjCKyz@m+F%b`J%-kRgc?dil_|&H5<R(p4B$3d-c5<ndGsqxD^Y{$v`{IldK9}5r(rQRkhO8jeC<vh?CrQaKdd=<yYa8wgG>7ZSL(9Z_(g$#?(>Z$-)H(V8A!k@{{-=v29Ir-vz{~%?DLA-pa|}7@gVU5+8<-E$M5auO-nanDlRJht(Ucq(q-z(IdDD3$e}s_loY*PK-_<Xn_uOKP`m0H!a>wdFk(wPrt4lO>{U~AG`QMh=f-FlZ#70l%2-)rp}DCnwrNM5?6&(ae5AO|LZfSbb?^^d8S{-sI#$=FJt(TFe<0j1Efb=Y<0Mih@MSWVFxtHnX~%l52pDs3{OO`Zrr1pVy?|^*@02-^--C%(vSHnaNDxbXtr^oXp=a`@<=H~dcJX=0OxS;QOkVu$VF!a#kwy&#3+f_vmN>Mh9RlSGt&!4hr+`AG3(?95%(Vj3uw~9l`9*R)tRs<XCkM8IQy8WS$<f}4L9v0@HBn&z)h9UX#(w2&G>T&K6Ss8S*?-kx#Uyjs7+I0Nl=JgKoofqFjc@;Xu(<2}Y{$EwwSv1H6Jr#;jtM|;3mFG$teK_RUzbTIyM2x2i+D`OqjTs(d4wA;`0F~!vTmC4J=DRvu2>g-mK!n5<^bzShk|feA3pAE3EnD^B{hr(k3g^tO&%JD^|Eb*C)xB_ksmV2g`I)(3&_-AHtko%dV%ct@Ig#FVuLJ!C3iIp!l@Aw@F|+PhRk=<%yIXAFaB&Hn2%i?w{^-d2}`@pU)W9ONLnzE(=ot!DZMswlCdBjCC+J;q@hpTj;dOZ6r;OaDsaPjGCm6F=g#a~Zw(387_N8jnRp~?laDY6zhZWUMEk%g-4M-y$P9;9rQ!^;0hr17{nUy|K7aei*#(r0Ljx4zKZvX+l<s=8Td~$gqWRi_g|e>zR2%tOhPn*Yh2j{d<x984uSQJ{Te;*43K4#aBBjb&I1_P`HZJ8(t8-UEr^s^?@P2}5yXGI#&PHGjb{wI7_J)Iu(9HNfg;D|6g9R$7u-;%HZJ_3IepKMvN&p(7v;T?g9Re9X=fr_<)BmefjE5Yh1)bJS{(Y*+vHD8~ueK-nD>I&Hi{h^KTax;X#z$*L?j<<wFQDDmM-S+BRj;!18z+iIOfA5(z>r3*0ejB#N>R7iy?cAI_?-o0qZ2dpaRca%iV9Jxp$9rs{KRkqdF3&G_k+ELWdn+gZ^{m&Ey)KrY~6ufiFhgm`*U3|?FPirNt)rf;b7f6o|S?S(1QOI8igNIpy!Sf95&0s3K}$g;ZP{*4&yy*Fx|2-ym1JK)W4bb;t5b+^QCtAwYf@~V*lsT>?pi-3hplUWJ7|uR5c>liJurTzvx3;Zs2|^bQ~HDhJP-;cnsRsZ^<g<{mpmCO-LrDB6RSe|BeC`>@`GajYVIKoX||B=wr%$u9a3##<H-;#L+;NPHJ>x5Eav3q9@D?d`0Ds0ymXLk(~$%f6Zt{UB3=`L9Pm+t8jH6b1iC0{s>BR$%C4K0t7ea{%Kkr#Ln`g23e_QEp3P_JXej#$IL*j@qKtwKbf5HbKU)I3V+lQnTaM<{*|4k3W-3UlHmmfr7Iu7`OsDc$;7yaM{j?g_WK12ztsLGpTi)cv1+fCvy|-LlEo1bnvVPPX=dqi1$;yD8!E*F0J<HIw{`RIE3CT+?#vQ`w&*#a7V}O3>#KUBP%lMO&)*h5mBy&3Qm^OAev~aFb`iDIkU?|gyl3{q;SK3stxL#==w`FZ)}>TTCq*!vsWPWhyknwsQZQAUsCQIpd|#ulsz18_vCe7>%caOrGaw9rzveu%FANeKWggSf1q6_rJP{clZ-0UE#d4<Flj62%1j_Pj_%Gm2mF>R-v3_dz$!UR!r}MAdR$aL8&Xr_=T<O!jkmHM{7+qcS43STlWfE;Cx<;qSCu-mT*<_kRC`CPwTopw8(%MK~<t>#@)m*@PXvJo6_{z_kv?Ve%>PLi~hU=!(dvW*+Ri9;*=6>MMivOx4+7RFl#Ctqqm6TNbR3s{^PO^ej!#!m-9R6bz;3p%-w>47Yue6X{qB_S262#@YQ1he=GEW|p-4+JkLuWa{@Yk32fa=Zz=j?K<8zLit*CwxSGzrsL#%1cmP-%RQ9+Xyk;S*dB#8C)%>{gRd9<zA}?j!}7R#DxVF5(vu8TlBg1qX0pxSqUb<C>KDM)y;|#hNfW#~FoY?J{WwQ2J+nX&@NtDt11ssd@UoaCTJ+bu#VS@rv9f)SHv)E&sPvvvBV`j>C*ktM=z$*X<_BAbou}Qa5=}XNZ72(T{pAQQEQV6t5pXnc4~|=^`;W2>Cgz7-xvTnBj%H8P&`O0KLKq)qIx_5Epo5Wc*>=FUhNiTdAR3^H`?p84UZq2;&bv^$Bdnn{~e3?D0%*h954Dk+BpiXe2ak?+_~^Z0Tsu%J0}_oTUpvxH=D8eV@gbv8{lceVNh`<#LmVfCrfAmpQ*)G9c0MHHUU}1Gf1{Z$!zTx&st!Qh+GEFFXLGYEY-{i-?xSRW{w?Fpdl}khv9m`TPOJEJftMz&7r%FqIVlj-Yf@D1ywPsyBL9XG4bHdZ$wF=Y?AX{ZZk=(MuV(%!i9j0f=<@>PCIew;-*9vj1ec<@ZEf0oq_fNzG8N`-GxLyv8CwI>s1R%Pta<<}BJq)~=YXy2U%AAIhbEtn*v!lk9zkB|918!UahZhR1gm_UFu=BIs$qD)ea<Z?ua3p90D{XOrIPljQIa(55Lj=RG)*MSmnX*I=ZR*G!R5=)PQPy4Y01C$k11uuj6HfnjC^vZNoe2R~3GkA@Hu*`^}M+(sC6-U?`sZDdvCJ}}~+&rg!x{@r8PlLfS{`He;+536#I!Fzulrrx`MfHIO0mnv*(7=!&#GQ}oQxX+=_YMQd55g>bU3$S3iZ=x;PI$8qeXAAf4w+`YvvV1Gtgk7G9uRV@SL-h0pz)uHit=6bfBm<{-K4lWRBR$7cOG@e&AwnRAbz+qCIx#4NNU2l9!(`USqHHkRhVeRSJ0<E7LIka(;TJ0{zDbKAp^iI9*$&4-dWP!{+{17v)q`}fW^D46Y*M0(`9U{PIm3libTE}&EdUiU)Kc+oz5x@!s0L=D*1&~PT~#`>)RG)LuuD|29U5}&n&BW=ZY5|vfVc1rbzL6vJ&kim{bfYL;^DUMrxorCD&MnX`9LWgZ^p{1@PTBeU$uoWLiXtB3`Sm~#X1M5H(ajVh59n0`qZ%=x4kuL2?v4E4RtNgazKuSdFnv{4_tDK(a})}AUyiz|FMPafJd|G0Zu})6{1C9cFT^BozNjQsV|lqSFhqjVHr3@WIT3<fn{CO=hFhJl~rJ4ADX@r#Pr5Iw4s{{a-yS*T>=q>sQk!uSZNBmbzjzOkMw2PIw;(iCEjZ)S@YdiFG@rIpTB35Nv?JqaPr8ri|w`&{j63n0WY8XK4S|oN|gfWC_$%(eJ{XROP8sEmo;s^iG}v_z)TCJmQ(eSReZ>;Yr5;}B}|I3`n#zUc0aj(z7MIjVQFwlst<d5!Tn|QZs5^R@-`-HThfj#$AB)BLg8s*S?hUQkWypOV|YZZP!OuKwYWENM1Qa1af9Z>gp02J{7vvTI8E=kyq|{%^d)c&l%~nGu6-gw2Ow2>4#QunsMOaCJduml51~tp3pF&oTk?|+d~RBjJ%N~+Pfjr_j%5ng4mQEpKOcp9OoiQi`BnF2M2J*EkKvvY7>z6*0}DkQ)+x1(K^#w(C-?pRboiXbdrR#Gn;XaCaL#8~Fw#Ai0v=7Wh=#WDRS#>Y-{KkwQ3xu#hdqmmjp{<-v6vX4_#T>aWNq@Lg;*lBC-V0J5Ssm!K_(`6Di^-oSS_A@rlvapythd@CrsA*NSnjt-1Rq3*dws^QpH*H{k;NCZ6cw@<0SL_P2@Tcc*v#M!;Lonw|8`ppDnF0%T`O-5ayFb)tLfw8Y+bUeU5nfGr)W3XdeS3D8_P&n1pe(33@|Wy-tnATjV!kh_XjE&Gs`Pc`>S0n;oZ9ixUmb&l1wpfOdPkKQfO?BQ>E|&UlgQV|mTRL=wgr`XY!#*`D<%n-29fIxgeN!|*s%;&wdJl0w4tm+5zkL%Z$0`1J-a-|OdoC$$XnY5O~%H59WU<=CKg=essDjmk^}M~SoY5sxw>2V7dv0`=TUa9L-8s~vGpYb@<(e7#`rhX!4!CW_CHiofw6r^Di5w9}w%zhKU7%s*D}(fC_tgnB9a;E8QBP+oKjH{26*v&YVGhDX|__zUZGGc)Ccan|U8&IO4uT<POnKt6Oa7k_?=i9Bj~o}o|$|CVH~{^Bkmft~@m{v`fEHH(QGZ)`<QydP~gs_CX$N~=#RT%=&J(8~nnd<@o5LD3fl34FEz+3K}coWX!iQ<-5}2*<_~K=E$8LeeF=%aH<g!yW0+O8;-^ls@CP>_OjK&M9Y()yU3Dp0J)gBc-ZNpKkQ~`f#+@8M)EF7<iDo?-Q+=)dSSA2*s+gKrg@K`x#a~Quh<LAXXKzKCxmN+DO}2AVSu0?Q0#5>56uLJQ%b|P4_l>w}3&WSb6kb=YLIqCy`=MWLq;{$O%QBBr`^LJ(p+@lIKLddjvb-PX4AD_*K5uT=KCrEwaIGa!M743qqSf%ZDJ+c(-zXKt~G+Q4-~jt}4yKQ)IDW*{1$u=Lw6@s+yi0h^Upk;NNGX3N`y=I-TWt==j><dgI91g%wkSxx5zf|M+<3fsrT0i8)FI(YnyC^RvKmtJAX@Os7Kd-~4ueswouvxwv0}H(wx3sryrhV*ujOPJJ3)DR6ot_phJ+mkq9W2Gp$71p~4~5%;EP?9TI-@+H~;nlEvzZqVO2A?3X#h5R{IL^BI(26vgY)l}|@P^Z3lDNl<^?I8O+uDB|xC^Emn{1O%_h7XkqCI%9>f8p#IJvD%OguWXNyw$E%kGlh>DuFGfyKjuCP`sWoXBpBF0Lml{ERdtTdk(djfyxfoT?AwCb0HsQ5f7H_4Sg@3F<;|=_kys|ia*+^-)1Qjc3J7m@Ga$%_lIrF<~W?-ccby)@g)oy;1Ip=M+Shc6~;Qy6HvRRlJENl1*!c(a*;}U0*@tIQw%&kb1J*=dKWQC7?Qt;U#p0kTt{>z_7@g&GM%oKc>f0u2u?!N{aPSkDSktf%FF#6`U3?UbK^;S0kI6k3{K7Rv}jmr5;dG3{9*5`)4Znckt#5{v8?R`%PS5Xi+l<SkY7eeO57NyjLXW&%ESiv6T7qNOVYl#BzNS8D;xuY%@#}YZR|1G7W=9f0yJhnjP0gx-rcRd!I=|rzHrIwmUW3vk7~ZZm{+vwZ>ksI*TobsaA*@>biqg)&<}-ChAn76N(8yw9j)O#L1^Rh4Ihe?GL<;g{dFL`vUq5QasREkuKoZ~i|C+wK}Q0}$FZ`pM}yWb)^9B|q|7sMkPT{MCX8{+E5!=-X%V??hMo=2s#p_XwX6?XD`)NBk8{mi=7PM>Ub?~{amN9W<}?7&%_*!Z?7V!?#f5@eFBGa((VD%oc4ds@Be-%t@Lhq*kU184SPuSf1awfWuhUMc^kerCEf_mw*8<@PS*D3^0yA`^7V@q25g7h1;6xpR^R}gq0EJ^=5BP-;IEtx;aC=8#FYkobpc53C#iAY>bd{YcnvZcNz3n>vNU3b@VV$Tx`zm}wHQ3<}B2}=K-01B7xL_K|PWJm{O0%+hF=fo$(heyxizSyhM*JdJn`}JO$c+%Nmc7CV;I#xDy(bBtGw<B-fesk1Q^Ks(y1<JtFYE^UssVvkz?KaN_DumZluoZC+!&x|?RoOVsT2~J`XsLI?4sxVG9vtfc_{9<1Rm5G{~<#-ntw)jac!CThU9v|fab)m#9YaDEeIg!4GHn(A54!{(sZgJByRv&n#_94)rLtqON{lk41$fPscnLoY}pP|m|N{dMw@2w;hD>#vD@SO17^}&PW(T})bqa?!lC9=3qnbiUsNb0VjrAw*qgujUL%+j^3u;@_&Y6V`On#kQ8D~~x)3lR<NHp-4?b^cUg}1XyGFevH>vlMu{J?avngcIK|^uB+Au&|s`!0ASieI$gCebAm`+Db$Ed%2tjK+l&t0dN%*}K%rnoIn%}Ursg-6<ABDi@lX+P-crj&1nqJaHELSMy){39^zP-XJZQl@@FHQNvDp|@jx`uI&Bhdfh99`$n;@<UtY|LAB)`_suFsepOC4w^R!!v8+r8=Hl&rcBQWZGC(-mYAvDp}<2Q>~hUg;6Nz*hT{>=i-p+sg&)X{V}@u`n=$0&O?JU$CST6<iXzjVmz$W6L#cf9-_sj*Dz&&LtSiaFes}WF-);AA%`7yJdm7NxyM7L4_6J<bm5D=W4B@l&gyI50(Jbu7nAr)%Bq6!HFlG#(E^;pGsXe#~J|aXa<^#mb<bE+Ae03j-+lRzIR@yy;j$x_52)+bjjwT)#EYDMMb5Q~Y$Y>lxtHM=!(TC3>L29La4id$8SPhU)p;X)kMUZ?3F>*RJ%-#n%y2X=~kzKPmTT|jMHVy2QQP={fiHyC)0!Z%**zj$X27%fkd_*Ghe?99kx|mVY5yh8lN)aDEjhjM)(aL`UFm=`GQ!z*Y#Z_k?FIb_ZjoHPVI`wa&Qygfe6Hug=cf#Si%m85dN}79Ujq6B$(-C2IM+%D|exC-4*nz6PLvaO+t%P$YF|)W2x(pFa!{(clHR<<@f!;qHJu8*hmU1*J^_z;^xfGeMW`9yy=NMJ53k&lW-)c_jci_6MlbDkU$zXo!;L3qbDp3nmu`<;f=SYxMCbS+WznW}O3})@*oi>mEK;UEakKDxBC-mQQjSZaY(9o$G9iRqV6YifL0-PE=FoLW&?LU<uhf{4fq2-&=-UQf%OUga2wiY4X_`%3B4gSi$I0SIBHywb(h3EqS;uGVGK?i}xoCp!G^NIyRqW>Ib$kp7f^lJA1LH_0SI^<tw!N!=i2ccq?{2Aa_z=?f(HHrK|j#SFW!w5z++9d$#`%kNWWlIft_Rl_C1_b&JI!A~UKj1*Lgbnd!)U9HP=S;EHS;Qu!E0oX(0wH-;A5b>P!sxD0c&7TbV^i!RF^5s}>V^3OK=U>utY~4?*JEu!a`R$R+@FrBh*TaC&yhMxcFf5UQ~Hw3A+8tscWJv^nz|+ENVGJ3s|aEURY+RE-NIvn8KhiT?C1Z$s5P;0X#HWMHSYw$)$v7F;z1bND(tXrn$q*5owGY1$0}Hdg!adO4%u}eb<j9<vDuI6OTjg{yl&`&?lGwd>NMp`9&(Z)Al>|w)fQM$7pg2{>GJKdwqBgBLjnimm7Zc8e8J!FG+uPUaR-X~$7lhOv?Vy!*g_JCt289@Bv0EfZrGF?oJEe)cA_-zqO;I>z1N2d57=OqA(KRl4~<8s$d^isb;vAKlc6Z$oJb)xZ{c9W?<}RLwC7{=(min_Gg%f9RCKD%0O2s}Jcz&@cqn5;cDLzaiK)D1Ma}nsHplT)Oe7p$lWaT8x2HPQVfxxY1;y=RC@#C^m-G~N$G<T2pI3DM>7Ck*xTk?;8Gw;m8c}d`cB-c=+Ywh^zrWLMT5yukTS1G6zzExm2Vpv^0R84s6)dq!a-2vX;8-IyU5&LO1JEakoZ>f(;h1(@gFp{2nkHO^K=nT=5goODPc5mK+3>bx4528B8725!F_8#YX2RQGbqBGfb}cE}GW%YCdJB}@8Ow`10IHA}Y2DB_GwuzVHS8YpBQ-hZ+u#Wjmnq-{xZUU@8wvV+gG{Ehm?rMq7P-Rk1T*!erc>!eckI|W`i1%$&hQbtydy=0NOUqydv4OI@S-!axxnVH(59AKPOHGgHxV_4%i+(QQ1QAxpUyqYw&E3c(L)uKlPGBc(Hkc-X>Q&BxJbdNOR0*v4|R(#CZFkQspDfpZsHifhO<y*Gde)&3=}+e)a`$lH97M+>arwTd0d97B+<z(A=atDl_9|_xfr2Q;!_%AlMtA=bRG9KbB)Jj_61Zn*egG}VVw(xoq-nC;}k}aIyFX*_Hdfdc??Bx*(8*5z4heUKpBpa2qLWPGX2jv9a`}(S<8i;z=hRtaeZzs_^AtgG}SgCm^+s3Bw4scn&Qrv$*^v(vW0!N+LV4FrP`Mg`aW5!KRIGvHixh|gJGeAh|6I~O$_IS><O_$4`ae={7ACy%?|DhodC6@-YjVFXA06MKuBtjj`+gT^xehkLfOny$<p+ngkn3F!W0<UciB#yOBbCt%VF&l`b6$4v9}DMhm4b_IGs*9x!k9wO(_9ld`TOXEOD<*N_ejU+M5Htk2K@!p2*rNAW^e^@hUL*ImBOvEyMMy(f<O53IcrG0X&K)xWI=TfqJB&?Nog7_tRnN1txggrT5NadV#n*CyG6gdmtY59X-A$(xbuIZ5R0Y$@G{m4jFo+iOgydK5FZb8Vb$BNB&%0y@ajR36*fwHan6&-?qNkAw%&{SSeaR#V)|yUjB7atNOt=Lp|X<G9OTxZjuE;4e6-XF)`x+$N6>0g3^|i(hslPvdv7Y>h)AqV1_aUNn=hXmtCDX0WZ&zmwW)|XyzP(%qS{tQj>KQ4=Qp=Y54c$`o0<L(B`D&$zn7lA4Nlv^dSy^)`+s<2}U;S`<xSC-2C@?jq@fQ{E_lX;sLuLyxEy6&f~9Pqt^YfDni3n5?-K78o^57*)u%q#7a`WapqoP(4=hs`<63TH?svA?lH@S+q`0m^IzW}R_iodX>)QyR3d|n9Bm75PVvEF=T|z3W6wiG=_{|5rZL4Qgz^{$rFa)0iOuZqG)qNCF*8P1(7#Yq41V~wec@Ri=UxWd3S)~jhY=Y<k3vEor^o`io13~#obw$X*x-EVHBTOMpI_qso-st&=*km2Z5l8dCmkfFvMk;Jg8B_GUhI9`i-|MVASuGU=ulbaTg)>eAi|{SU2^E}*>cM>yFTpBWE@`B1SzMaGi^$#C2J8oh-bQaLHK|tliI~1=JZ7(`OiJc%eb|!HdK&mAUKA+vD%8+=U2LYEsf#i8FovD%np9e?~+0CWH=)zIis)tS&7lq;tr*0vBPh2=-3fgL?&VPKA8%aQP`rO5Q;W!_S6JJ1(+rG$4SbQpp-bI(@UIoQ<ad{Ajqq4D37%8jA@Lnfn%QOg&{8{k&=Z9=@fX7Wo0~T=cHSbn$k6D4g8`_)$0@;4Y?>?{xvDCd)Ipr;uc|=gem)6x25^C@A@M-{83>3yNTuw#LhM?yf?GYu)bT(!Ju`gEcISqd(u2g{LQ#dEXdE5gjRq5&5BH(M2A;ie&O=wmbNll)S+PG;`qjy)+I99b~5VDAT)oHO~9`)T}Yt4K#XxH-b0OP1(jjFMv)}fCz*zA<ULtYI0Vm}17#|3yp6ZX?215O#{AD}Z}`%A7hag8uV&+am6`PV3&P|2!DIGR4D(cz!ckd4d5sDv#%7^ar`<X_+9^oI4E;Ly@uGL4h}UjCoNq)Sr@coiB0cM!(>Z4lCNyq9)X_7>{}rNE-gJ*I*TGIExS-kj<S)KW>Vzx(R{}~l(LZ0iZBC4Xvfah&4!eRhMPm$>xFm+rP8l$PfEoK7N5oe`)6!^o#THd%ZhHyFOUE@}z@RCOOBM!OU|ZIY_`m?fOw&Xzdk=j<y>6?l`%xaJ0DkM$5dQ{5GXX`>^`u^056GE;x5#f&Abg%9*RZLlPMevG>0tmLli{QR@w(lGnhNcXpOPV?WcGb1{b{cu?%9|W0U(P|OG?iO`*guK-ZekbWU_~l1EJw}afhXiW>*LbL8vkhO>0Tr-{+%0$KLx4PA)h`_`6*#6kPm(dfm*1+fTO|5H{aVhp7r^71OKeN2^&%7($9Sj2DzEea`7J-iQ9<btgfplsTD4)AR+%=4XY>1k~MWvXRB`RlKUjAxfpdbyU2z!0ZQNEr<tkGlhbE4#ozN5#fPAj^WBYs&iVKc(@0b+8#BoSgZ(vfG<C3)nyvB1(~5X&*9a1gx7aR<l{LeWi(t{9NMBySwIN+)8n6(8f}@4`TGb){b~74Gw-TRP)2yd9%|#*eB8H899QX#r=QrSp9_`>!U+@vVZ~xk5IQ^Tbcp$Sx7ye<xN>3Xd?>6dLPwBvPGg=52Ns4;mCAp%;<rBADXp`X3p1$4mna|ooOh5|*E2CV(RsC05px4W+@OC`FicsuM5E=_eDb<hNbEKy^czCw8?kk9b!P>AK(kw+6+~`=voQWv5pNpswAcr&7W`DDo`NG|x`?he{1$jrrj^_QqTws&=YHkNL*S=DJRkRbpL#Lljhf~da2rp8-S{v;2m$Gh)~USX`8IJn#46=-s1WwfFcsZSSp>tK9x`61``k;!<P*Fu(6kw&)7!X<fX~!Xs27|aEpaYzQ0V=OZogeboX1H{<E=Z1ziv9(_fq@gr#S6B6j%h{W2$@y!iTYq%z%4R9rySMi2NX8c+E{Np?HV!d!H#xOA|lhqLg~I{3eFiM>NswQ1|Rl{9+>S4CS@NEpHAOuas1BwH$Lu01w84+E_$QGHBzEBS%En+N_Rs2NVo4fu9?t3(!>dr1y>K{0WlFS0`^iof3q;g=KYva}E5T1GNkmm**9ml2KH1rIR8nKTN|_D|3tkz$`|PnrI$M!}>MfxX(O&SyEB(tm5h-UQ4$4FI_o?RNA$N0L_ZK#Who^EOOnfcd$#rGe*&W`%EZU0rcDse<zU52#UaE0(6+dXSn9CDU%=7M`AM3Yj!k_Y6HEArfSIkA|4yn!ls@rb%ubF^_8y}OtDpZh_qz7ItjP5g%r&_j-*I&6h`}Wv`Y!GiS7NL^&jgrSjK|B+E}oDbJ^ASuO4Mhtat$Yrz)5}_u`NG6ro8B8;p~m&gAqcT~X6i5Dl3xh4pdQ09q=(P^SH_XMN{|KwXIcRDzXB?%P%PPNdc$fk`p7gqOFaca_kRuLY`8H2!WrPOPJj@|T6a6dv(LBWi#|XZlM2m6Dw&Y@Hjh^(uWCwp+|2cflyx@*RD9SouN3j4LGpP8QwAOfz65#=K3cH`t&T(TD@ByDZR_qs^adJ~uk@+H^B`y|)s^V!OFk9;tWP<~#S0wFZi7cfZ$7Mz1x0et^11p}?~P6QltH?k=Zr_pOR(2{npt=F=HkFgrGa11B+D7Bk?3n_bmupvHIpnTpdd-dIH@N`JyBhTR-ehto;Ss$P@x^>p-*R!8YUXtnTtzBTU=Z8fQ-Ab9dL(8(u^j*=Xb75nc>XC6>LAufTbs~mA%5r{4H=}ser%rp(XbK1&W{l7KaS~=2ZtLuZU#6)o4dO`Q<I&rZaTW;vO_>=%Igb61j&iFa=J5yZWa3xspJ@LXdx~9U}sL^>tU0J)!$%FQhw+wLxU_{uv{X8)%$JI}zk2N$l!+mbyy2C;63vj7QRjLSV%9dxx6+X(`$C)I&s^$*vYtLZZ-JESRg(l!)INf!ZE*=W21V5NqQd3NTGL6)_P+80BO&|qJ-y6DSLgftJEH5I+Inp4T)I@*Ty~0@SqdByyBRb<wfi9yrdcrshgWjXo!PAV#Sn-odfrgM3Qc30BHU_pAZgyx~fOo(H3fG-1LB9Dj^Y&Lj!X0VsT&6KHh1Sd8cE#%13n_g#sm<o`UyH$~JFmC7%?90Pq?ae!pS(P}+^K~)jX$CW>Zjyt8^SgXFWerHB0y1?)42%83`%}@$ieH$pdVnie<E_h!SOfNvnCNV<Q$*)K9yNLO~kI<D1OEmQhP0j8v9=ADSmisJgBSSn}Y-oPkDt1t79r<H1J2hpzG8jkMd0tYc?*D!<w6wV)%ABMpSX3+gaFN3+WrHnOFY)tWl#B7@qJ8G*-Oon)~2UTx~-|uBYW#VjiO*Twn8R0({0S6;(zBb1Eh}so+JLj9SJJ`3f^@S5&6082&0zF^R#Vw~Z#g;)fS;SoL8rP+k$iz}f$hKmL_@ggG{qVjp&#3x)0Ym7fW-`o}g$|FaqeoYr-io;4H??EuIl>*GrAHeNFn&6T0~bFUkZ(GKoiUBQsL64wC6G9Gg;^erU4X|j%?zyv-5w0qaKY-P`(j>mR0hQk$=ORz3w(M`i<Ng3AydgKsirC(+Ol0-Yp4JqS}NQRRglL}Xo38>=`^d*D>yLBqYp+q*ueTWtd5oA*TRi(g5!#w0YFkJ1+20EWkW{X==SlSIZ)N@!$I+(EHpX!419wI~Pdrp0S9<d_yx<n&ispOSvi;h$590Os!0F*aDm^p;mj0{NtzkhE(_4cC@9dEmU4et<z6T|)YTGHOdc9(C5N6t#3WKJ)Kn%yx1y^e(2wcZdtTX+0fgl6cH^HU6He~i%~uyP_P4Zpi>fyxmj1A~xll~OoX|6^;OmW}28qmZJ0<*)OPCUJ5%Gy|KSL{Uq5vS_OLCly$WFe?u*{o7#oEr#H#XjKgnBk1Knjw(6X*Ij#7m=QexgSjcjM|qPM9i}n3L)OMyUB?nKz`gFv#9wf;g?;3xZ%OHc8s%M(@XFonJ|-*kboDMX<|M^~x(;Mc_@#L+BC?vK-I)XQ7>Emm%?1RxmAs7RtQUJ~6Mx*iVQgi)vBq6SSSNaWg=^94)oh8M4H*ke^*$+Uw-tXFmKer5qC7W(@lTSs<b{NZ1%FdP1$ROt|CaY53JGEk&Uu++J!4WCo=;*Vut^9oWgib@x-$2cGJDCn`-zyhdU9J>6Z3tY{V397h*20p=j&x<DKd+SjJEwz_?YmCu-TTa<&f^UdD0D%$wOwdqRf1+xCk9L1#pw#XbQ6mUyoO!dwP#!_|vxJOO+#4@R*FD=G-5@Nk|3`Z>8Oj#)aER6Ni7*s`_#f;%>MRDE+R?TQ(@btae8`!@QY@p~x(`GOs@@O|czB10Q)gxW|CJVdHIp@ce|sp@2Cu-djLc#_dlT5Doo5<(=%Ui6a$bRywfe(6dz77Aqs-(J_<ZZ*Ujs0-rO8z%9t^s)E+7FWL|G+?p&zIV(*p*z;k}M43k1{=9wprB?et(60Rm5ZaJglDY)}-FjXOr7eICutavF{AaRV&^1zbLI<fAdLFykoIz*<km95O6t+2k!K=Zb^uMxN*ROJDdhP1(4SzuUs=gle#O|*m9VAVo7ez!wOsgEtLBZ-q;mip<BF}s8t<R*p+!3=X5(qMxwsqI3DC(!b=AHq;C*rfXdYHeFa8Ap(9Fw-g#O!hMNbdCO_y>5+00)g%+p#yo#KSlyfSBnk_=}jd(V7rF-`%!xEzJ)_I)&pg-+^~qw!6@aNn!a{F==G2*K#Ttm9>dKfnit#N40%Ck;}HaUPPg6j<<_<EoC^ox#N5_UU?)Cs##T-tK_1-BFjZtM~v0s!d^QU^8xcNp1+X=MFD`^$d*FRN+dc!;zO&Jq-8K~-+rExJNEnc`j_C2J)nU1K+eYg7z21uDUzosR(f$#H=9^xlNoBYk{79R)LffERem<)IL}4(Jx~8D#UQGOQXUDgi*@K%RR|_^D7i+N#b0MlQA^23k&;BG<gVDutUUFwGUS0#ZR-A-uFiRDhTN)D_39ODr%c_wQ}1W7w|$Tz+st3z@^o0s_w}KRh?!8gIt%E*?U`dId1<xdb|}O2RfW$Gt)S6<(Vn#WQo^?Lf_)1=*Q4T)xVId?%6B@N9K2Di9KYkvz`ojw<6e0|W@rMzDTn-9BE%!WNRAJIm^=zzz)o67+3Oitin`sM`uSRYXhT##Og&|`e$aPqbpInP7NT(gFaDJZ##DYx@xmp9OC|A_+s*5%%UI_L2c5O4Z6T#FiH@e!c6342p+b<$B1S7>j{t9!iqB<+)V}SiZI<#)FU?IjcqAQ+;Ez*sPndlM1-v$fcts?(3OEMc4^WH6t-%4SY8Puv9<);z9<|}qIybSkMhMP=Km=&St}DC0uZjYf_PNXdPxOn>DyUd0KnQo=MuI`fPZ7#V(W{|jSi?u|_Jhp-(d?mBp{8&O1<_dHs}=z^U<88J0_(;p29p^ifd9CdDw&Rufs@w-EIk-_Thy_pFeix;zFXmcp@|^zDgLFEY?_TmPtGJzJW_$G#3G92i!aZ@izt6xZkt+a$S}(5<jlhBn(O59#mck+XUeWBTdjA)SW51QeKiiugxPPR+{dA5sIZk75GQhx<ji^sS1uik_H=Hvq1aEQLRGZ+nc`cl;keM9K6ayTIY=@pKnAxtKlEuZefp^J5XMCZIbpL(vAm0zL;(EUMbGQ+-i70a6VVXas76`!Ol}f|z}|vAr#6AaFCk{b_{SxN-f(;A<^bj|mXko;{JRWc*46cZ2DL?ca?ki~K!;CkG)=j@8K3C~K8_Kh6ch>xPm}mmGBiOVVPR$SdG?TNQIvQl%{v^P#~34Q@l*ha+A39tB%l)2R#@&*wRsASs0DZ{E@U8M1R_A{NJZDQ(xFQx7)QQ%U<5~UhGI1XG$5-}(uzNi&gJH&oO1NOnT?TddEAhe;jsrZ1z?@EyR4VkG(#Cgl1CPs;*5R+3b<1G@aX*0BHHys_Wr7N9DINpY?1gzd+Rq(qJ?N{hFBCZkX7eg-jz)(IJF!pSDu5<3l&)Etun5w)_^>j@?|a>rdgSKHSY1co#P?BG2otMA4&Xk44h)>KcCHEW^h>PFxw8b()ql!+iP7E14P$O*LL#R@m-Rmx&)>{+Xb|OStZ_J!Qnw6i#pI8ziCJ!ljfDdXTaYf*;<5Gz1`%eyl(-rB-0g<y*Q5V8=6ZVmGU)%R?@t}WV~iP=%;GFdr3(Iw#w-kaD4Mza^glloGGz<fp;>Wed~@SiSt{mgDV-qjv4eeM~P1s17sPh$#kx40{~WeW4^yCCyt4kf^%6#$uTR@an7^ldR0fYc#z_=?6-sjYlWzIrfgErVwQ$-LCaZ!2v;#p$#yvfP>Oduw!Koazlt#jtH|6K8@}oj+S_G#+~1~4VnV!&8EP_045!hs(3uIhwGepLZe-VtP`JA>PRhbHh27|YG0_pf0J^(ihHv@tA@z|F*&;O~tvSvIP4%MO?cT2+i)+^k1fCNuYVRB`F&%F4o7-@okv;$4it}YpF*SOeFfBjptN_Aov_T#XxKID(5FfVbDiND_JdH$upuOW(_-obXRWpK?bC=3&IcG9eOQLgc>WN-4)wA>d)I!jy@v(LF1N0HbFF;`AuW$T5sp%RO5FEzLcA7;aWD?B5i;60C->_W&oqRpSJD3>C0%FcQTN7t-u*kXfSt2bDy(WM6_MmTQ3!}C|ScUD0+nh&0(a>o|#@vaAIhm!!b+g~g_lUn$?3-V}nNb~7vIhHcY7$lLP$3alN9$_|jfwlq6(804MwP>k^8Q&y8e`@{gr+nDIU)u&aPtt#(q?#H=qKiLGIZt!gR~4)-xX@IxEHZqsa~?qnOV>9KEoNUlN0*|+O_$aL3ym7Y1PCaVKrBii~y<`u{kH;--<AE07JaJG{W{M5Z`XCck)Z+TGI$fqlM9ULH3vz%!gHeAzbzD4V(!R{XFBEVT2gs&zoP@Vo@HQo;=MujS|(|l#nT!HY8{m%YXlONZ*@`mE?>IZj{~rG7baWIq|fYX!W#w6_zP3jyBrem!a;ub(}Z!BRAv<*P5;;H%n9%y6ith@gcxdMOL2pWL*k;wIBYpkFotLEesh6`#;mZ*zdU;y#`(pTL8BXvP?=n$YYinyiGRR$zGe!mc%cwPsVEvQ-6l)<LwgSELYLyjPz640ZD-6f6ibhW2fmg*7Lhz<@r+EL+^0c*^f<{ne@N$1f9f+Wd8$z7&!wjG2NVi`t(UE21&etJMA?2=u)+^cF*%jYwnu9wQAZ1Dt|LW9mnIb$@n=!TnOyq>;hNrf;NF>agb=fG5lX-dzdn=1BTqt{$A=Z5BDu-U0=$c0$@Ogel%v9P;99@+4lnxr&@w85HH@d;QB~;<^yQwF1t_mG!J_BZ6M=+`nPwitafp4j`fd98jDjRfW_3F%5KS#rqOW7Ztm{Lk=`gbv0TXE3HO(0Y)fdB=p-{7720AH@-bCBeF1N_TpqctI74I?28sO$OM<BGnacd)M4VA4W#D%psPd6OL{eY_YDntJ%6`v)@A(XcNs18?YDxI*U=$k+WiB&bB}LP`_@fgx6el~Z$#=0mOtxs%yrIS=3Ktg2sV4+Q=tL~qiZFZ%mfW>Ngaz5l(7OB2Bp+&aCS=!boHBS<lW1=+GxBiTrP5dU6EYZ}%Q+9`jepLyKQ*67nk>|9dQRFVwW{lC5hLG<b(QmjUsl^CVaF}{K7T#tN%1tPyoc47xZo57x*zK9W?Ky6kiUFtw|;%cHcdnoO+EI1!^CG9ck~M<FjK5MP&e>w?=((l1Y{uO;Og~FJA%sRv^4&W1*Y-0y0z=OFIb95lkPoaGn$TgoxWxkyM2H%&RYm^0}BVhAJ*3Rop2$QC7l1k4al5kvIlS$GpAC`dl?5p?bH^GN86r@71fr*VzN=1r(euD9?t#|1%GYb!nDTRLxC`jAN2VqQ=J6)3MGtqOOmo2dHMM%TTP&RwEvH|t+6PKqoTGh()jCBAS^NT;!a0e>AeFB_DO%#Nk93k&^L43Hc#d`(!byKof$2^x?BkJRN12njQz8im+$(rOpVPw{zabc?-faScDyDR>d;*agB5eRQ`?m2K+M=M1+}opz8)x|aIo(v@Z%Yi&HFa-HSYb@HjUBeFr2-8(9pBUL94DY0(Qsdb$mIDawqO+02aqLx+Tx%-d4%9R5u5oUhqpRTA(GE_>y=Pm}@yin<%GhN=?)X8WMB&XU5F`FYRdm+I)8eZY;4}R)gjV!M68$#_TWu(i{rVwM2#;zNn}gl(p}WQFF|;1hZq?;6L6m5*h=VzsY3Mq87y4pKhl=rSUpaX+hk%slbKll4z{rK*4!dR5$#2PwtcAB*M}ttxNdOk6R>!jgP>NZlmA8sUUHotDc{S0j-+qiObS7=Mh0n;pJb|8{8zg%8>)@k32=^hB1&AY-epS(@VsWxZGbU@Kw^b1=;`uL)G87{YIk}gN32Jj>{H|*BiQhj{*;(w?jbbD37;MhDb#?zZ>FEu<kWbSgTjnd&ee0c>{WA`iGX_Z{|~ijTwjvhX|=r*7xY`LSKtpl?pKAw#z{>#kTB~E+?aAEyp`*s$}TOV1qeppFgo+HHX`Bkq5q^5~4Na_#VvDvHusQx|orJaAxP8L)*C(*Z5h7)M45bvrG|jcLORwKR$?h*#rb>7@!3|D(I~;vfBj}(83<3cP#=tF5Ykv7~)4=3B@<FpljaYXt@=v`EYHUEaynh_yzq3C~3?&K@;tMancS~hA{7(tkh6=G}~fM8=X6&0PoH5dRGR+R{kS@=DuxS%a)_EGqVBfpr3d=fKp2f$hxo*54<-a7E_dKkO?!B{Mx6p5qU439*`?{`=`%5&cqSlJqbyHp{3TU|FGEyTpy)U7}$qDg?2H1*M1PIKS6y3VFhLMwnc3nyay7Kpmj^h3g*g(O=O1@oV^<ZvlT*ZGR=)$_{RUAI_c?vf@Hab!T%rs2GQ#G0K9t2J7H8K-^=&>BIhAZn7+51Nyafir;bVXQcc2Au;$l8ZYnc4#od~f!107&o9Y+iElkkq`uOElE$FB27OTAs`78lDBDFtYS+!-p;YU=28Kt)|CBN`&U^VclEt(Y>j{jnlZkUX1r2a)@UNgA+itUt7b9P{w4hNHH3w^gx)sXjO?O$1CcLnKpEi<1An|c|-<_H8IBG1>Q*G#E13TTBBRzGa+c(&rWx7(>2b$UyH2cnE0;8+<o?$cp=)<+-a#R<vZzZSu${8x00JtPh1&htg!&7s<GWLuKQbc`(0+rOAOymSN%Cy6_^>oZIpD|B5M?~*b;ee{qDDKY+1)bHRP(lKh2-`}UJnvHL6JcNgFp-)dB;7?LZy7-vjnN4MllTb+38qEsUCy_!`xz5Bn+CSS4!>o+a)r_NKOinQ`90J%OSoYeU`bGbhNA$lfzOQW1d{#}e)m-6%%OM}x+AWYTlxv=dp2B1M+_zIpC>ONP$?d_AY@P5T)~CUr1DByaqJIT#0w+#A7Sc*f;Nqd7L}rGGU=5Jhz()`qb#^e>heZ%(Z^T+(#${A5FQ*AGnEdWt6l=jv0d3IXc;;i4(8o4HM)3vaJNK6(^d%Uxs|XK;GEAZ(bLf?mXXqEV3qg{R9UG$$MlT=*6~Nl^;PKmEnLPuWe=9_brLFA*Qqo1FxXC+%5zH<cOuQT$0#Q(bWtu=83z41-xIJd>3jTEXvbJYcP}C|b%2sVj24`%5WJbq7u^=k~Q4I~T7WElJPq6w=ACgj%NaKC-yyL*fc+2*Qq{3n`@iT9Vp|D<T0&=1xR|NMp5#zKkPX}BatEP863Ha0cQHc^A=R!>di~&OD!23lCNJdPJYQ;^(3Jz$heU$FPc1&%oA0<@kXe3nYe|Zvnm3gqA*#+LQh9AHqB<y608b~d^TPlwZx{Tm+9=OG;V|B|p>EBD0+TM3C#{+g?4(gTs(-hO2zsvCeI(VQpAy;okP+LJmh-3YLO914qLXEMRrh6bRig}}M7HWeD{h`c`lb!7wva5-wSA#<FZ$6rHxSIUC#A3a=dda+rf@qD~RAhq~Ih}iZ8%$CdQDM1txSfntPmNN4T&v`2GjnhG-&%=7-W{wng^az6`QI;Vy(&*;3Bz`$qlvhBu(z}`+Z5t;4i`2P&_XqwKuko})|!ncSXN2LVQs7^3kf$3Kl(P~c^X`_y~}>zTScI9`ek1G%^!Ptk$ubM<panU`1;%pGR+2qFl-VwLAA-+VFEUL+SA`z+BYAta<QV;uearw`Ly}gsvej+hMgGDiZ?;*=ukabK68Gfs;Jm|@9mwV^JETZXlB1<NPO15oIvYEEI7)BNCqf0I5rm7*dHS05>l@<<a2LqoeoM*NT|w6&dBBY`pTC0{psEcPX=?>({*S&d5uSEsgv|FO%Zqwr?%fkr8y25xMT8k`27Cc(Qo#1TVmn={~J<BLu$TxU7qg`I%(DHUuwHsZn@BCR}<?#ckX%+TMgX$%9cr1PUQ(Fqaep_fn#XGM7{M9v?l9F&Px`hkQ%${@!;HCZN+q%i{9FM-9<QOM_t>=83=;%iGcwB`wP^$Zk#kaja-x0k6Anj<CKejFBr`Qnfuo}f6<boOm*d$s0J|~Htz`l8y5+qrB1qUG$kWpc@XD{^NGhHrh*WktYagPAfY%vnCsPfq~LhP0^ww{y$w`+S1}h)@P28^D^Pg<hODh`VWMJdM<rAc4=i~(C8B*=qM7wp#<3YT)&6g{cai9y`*^?v00H`|FooO5ww5pb;M9;c&*3MiHp?hEX5jbtw8N)W8*E^LQ=@J)lz$4vF1i3Pe$JDkH&lI)Cn-S!=9jz*HGyg5We*9k#?FE5zn|vHA^Tn5WM3rof~`E0SksE%=v8gvQ^A|{YK2?U^mqp!;r^6T5?Z@9p73RQC!T8#Z>`tYJ40_F4c@tp9e7;N@!~wabm>iVQ1gD(>mj%dh|~1rhtdCH+*u&fufY<m&!B(Y?;?5G11!D#1jnofBSzD6{>5=wiE3&6CUY*`XihTwEf1=9V2S?ZZsh)$n)^GL?uZ(&%>!O>1zQ}fE;Zw;AJ&l4WZ1Jk483beJ6j90onTtO042@ZSMX9&99Auv5)3SGw9zaVx_XM75_(lk>6xk%qS_j<Cq_<m`)m$S>`Cf&H9;Zz8n`}tB*LW=g4Z{^TLe~lV}<dJgK_?1ePQV_z%COS52Xf{5_ND}%ges5-Wlqu-I2NqVo6#GIYGm@PG;kpCb3n8q=SAOTxnL#y>`832F^<>hzUNo1|SWO+f8Uz>Dg3-UXREaa{5NZ*yYy!;OwrMex9@{8R=qNnIqSHhV|ur5n1Gs3#jok36!6u2=i*z!@lmP-K;ozH6SKt%bUlAHU6rYh6jKoe|M|BzvcvHWUsVCLO(-PV0;~7f)qL_JS*}N)oxHpVZ1>;*b$kpnDi#=PD6{)`%nS+kJHyNqQ^K;Tn`-BHZf5Of8zJ`uqA*FUx@WnZ*2CQYDFi=NHSz3HnU8E{E!qJ5)SKL73hs2zW=({drb?p`GJG`lLg@Mgzlg&r?Yt9Mzt-4=vX8<u))znqnQv>y;QP^vc7}V=2RjVpq?C}-e3GR7G0v%*1*(L1bGM`td+rRD5nn)1Jr9a8F;{svXU&6x*)YnkAeAIyUZl}b$zy{%y$fZBkm8#GyEb4T_kFAg(6i-?)dV=7H0GC-OC9@S0%kcAs>`^P-AAu?gRiuY+I=${<ch%8`O0O*6(pCi^pTe%Bes}Q(#JEhBK8nBCu&?xaAlOlx!6vu~*-W#r->0PTmP#azt$&k$Y$5T!}VOrYr9NudKYjEnP2T4hMd@Xx8o;o9R3}Xd1wEv0fX^Q{QGxU85A;CD>Qd^K9;JB3uIA$4t5suFZ5PsEZ<LET0B~WGyQ)Nfq|O46Rj7u|sR4UwIY%JKm%qoMf~g0WBhi*_x+Kty8~rdx{yrk;&If=<S@MIFyi~46d|l6*?ITO{34#kX2Y<nJIMjauLiB@IxGKgQ90ovc}q%pt%fI>U<p8-r6KX@j!8JQ51Hsm#kb7Z4o3NWK5M`KNtr8)K3Ml&afLb>VquwUpu8~PX>hHb8%d|;XGzxqH_<KW<UwyQ8EsHK`T0&AVnI<(0SxO!&DwbEYFZ<wIs!Fx}>Kbn3sjaJl}4`#3(Fm*X-*d)mUQ;fRQNEYgB|<kD>)i)OfmYMbVWM5%0{ay_%=!7Vw0~kEE$0Ms^B?JtdDP*sDCX<=Bf+iIY?(c)j^FI0hNN)ePOb2if1rXVileMNiRcMc5A*v94G3Pv2K2LuyrH?x8gJ@w7qrju*j<TF_}s4vKX+t<mK>`GZjDvN+TkITvc~S_yknoESXy_HCX}A7i)LPwtAQCJI9Ez_Bg$s!zVtof>In#*bsGOE<=HX5%RvpVEwrcRWnEvOx)|35ntk;q8lGd2i`*u2Tl!R6GZohMm`+H#(q&#&1>dFOMUpB26-{{Ye83%3zTcA6|)*PHNTk>rYvVd7QY?>MH3$nG?6ez?Iw13VPgujPvO6>`Z?_$JNDbt_LoQeX%wK2eVa=ReGnU1iUSQ3Q$h|&X28o@L}s4IVSA5OaDWQVn~>i=%Q)u4EU)9X}(s0gC;G=jdd9Zg_=kZmxnnGgW3fpTieJn$V6&$JdEUSmPao{>?qy7rL<Jl*-V@Y>({0`S|JVBOQXK6it#8BpJ}I)w499G?|>AjP6?;7$il{+!6Jd44tebKQhMn%w)(G;wPZ^W#DSm^zhU|umjYUh|2okzi>-voU(BQ&S9}f#^iU~-8|H=u@1Q3%1@7aLAdpg#k%rV&`eT*1%d=F|6ne&f6Ag7BjSSBCX?LY|*)B9b7Hd=Q<E_`dQ`x=%g7~@3$JcILfX?ZLffCX3WdR<kp-Y1OZY6f;zR-c{QiyO@atI1St@R0<WE6Du@4l4}<`I!osCp*KY0)*@Mqo&S$V2VL<aT*XPm&v&14m~FbwIa7*RI>^CY{Y20`^!R-;uK#ZCdtdK{hoigeX!U!}PO-?D2A19qkwinVr`is*?)#srZVC5wEaw7Aqt+s__JrAFa2pv#m%`x{opy*N`kNgJga}@y=lE`a>}aM9!nLvcS?;La&l0r4oRiqe&GgL@1nq156I4Gk*<Ly)^ZLIZV3`AE@9-WQp8XQIoBNjR{K8{ecD|;{X}(5hF*rsc}+{&K%o$tmk3D+us=Av63!;@CBj7f7jxu>SNVN{A$jyenow2)Jg!@j;h|PJ(zzDf$h-RquO@dKZs38y4=#UBk1bkoDIVSL6YkG@y%TbB`S}=;ogT>z_0tDXCCP*$T`+X9cWz)@Q3PWDqPQCWjA^v0<4NP`24LnQDE6Mre8U`QY+-}c5Avf%PUyisoS7i<uRCR9eBM5fMGKJeHT~`mp9E`8L{oRjRv<4KB--YWLgboAvs&Z1<@<L{(G~j=N+cnF3D#}UaGyFK98;-@eevMkF^Uty&1Z!VGtkPTDys}Zz>(JRn&x;^<CiX=*H6upq?T1u8ypI<^TA)S;|07AV9&dt|gD`9FH%#a5ZCEn{;s$|DjSo43Gs$fOQIHD|ccO(Haa4a%bcpq%R<`9bt%IH()2ArVP(%2DrN4_CVs$!jufEkCE4F#fjcJGVt?v&n~O_8y+ip8@FAB1*R!<zxnWb_xu@it`120D<0_?8Ssxc>24N(nj$xM*?|mda&n3Cg8xx>0zi!6zzc)fcxj5(7^>>$WIHFR`VaW1A?$SUa~A}PALGbf+a<|`-3tkw5R0?unb*>FA6EN0%8vdtblCKa1n`5#F)V%NHeG1ec7OQ+l)~fy*<dFZ$%Wl8j9SGLh&>RFf~l)*q>)0e^67WgU}tqGMw^@@T(7%_+7E<)*Hc~)G%bE`JwL2#m#<?a8drNP^alra8%1+kenwQ6N@Ajn_T%jnpF#GD-tScTQZdLcB-L)l#}@3|NPH44kP>AGk0nE4`AA>NMLr1+0IvAvgQ9Y3Z%Zig?rzw)<!QbZlq3?d*JH1TR~tT)X#`j|Xl)O{rAt48VI6fnWiT<j;|x!Tp#E?ziGDWWS8<_#E`!Fuf5}G*z2pLsKuQPOAZf-c>FB1}4FLD+^4MpU^m2*i(D4F+io_e|!)b(aN;xL)1}@uRLTFP>N@LnL24Z@{v?6jHT(znYvRlSC+hFPM^(pYj15l_5?`0nWSb`<_OklTE{?=RpK#u>Zg5cA7TZ=qMkw)r!(=UnZXek;XSLjA35Yw8Ld<Bf~urt$L8$hsT^^?oT#$E@ApC(x1vFl$j6pG*jZT;K&-9rLKY_4q^@KkFQ!Wiqn0D+B%@hqTI&7Fw0=M^W~NyrdLtM}w?=p9&*mN?s~H((2=3xD(Mc-`*v32dST958D#_RJmDRNcEK{5|h&WrgSO8BnqH2Z~X<6Nh;kPvhw=lU;S4<aSaYoep}?X)+Iqm9&vgaJ&+f_YR1274X13p0eM0@*(SJivMS8r67@&N?C%M7BIW!=Q&xT01gGwfAvB0OtHB9Qe|rV)0pH(18HBu-tLbWIk0)MKF&hkhtP{y90C{T{6LMuPf-MaSn{LnsfcHM`_SZJ+>$V7e&5Z=c!m&CT*c4@_LhnYPVtpq<k-e-hz{|j%}XDO%-5zlNL=zX5xYj2`E+0*W7BC%WI^{m@(R;?S4pH8vfIu)g$CQ?wJBHM$zIie79u2<eff1>lM!@^OWA_dx_ewmbYvx%6#o+@M9N0S_!yxqLbG17k#(D(W-%zt`o#7(<X*pTjy@mGG+qN;M;Bw_Kp=jxI5oGv#dajnUJb!ZAviaSMhbyaMyN+8q+Y|#y-x!=GIbS65{N_6TqEzq!sW|u&W3~n*T!}wLSEi01W5*A0ek`0a!A@t8nk1pKe<rQH5@8rP>~Lo0EonrU(5zZp_g0b>joFExd~1ZBFHdm(H{h48M8=STO1o?rx2eg3zm6ziPy{~>CA&v0`N6vJ1*$64MLt&!rvkC&Lj#Y^)b?p2Kc!9Ppb0otZI0?L>7)Awvu)%Vk6nSR2rAmG<0uG5Y>Y$*IA82EUz3yD~g26jJ$S6J5IiRV7rkODZKq1d@zey+DHNsjOA|PlzXntE-&vq-Y#OrA5E557knK;i+t!9V(S7UbTF`$nlPp2J<yz{2Ls-Y91QEGl0#U7YK(PA`k6xI0bAxjbh$`bU@r69^C|K#8nu1_;)aAnoGvnsmd+dOQ4J7}1LC>bGT!ZhK2Jb70&|I*o^Ayx6X@|_8h=ih5C8i#i8P3Qab`kd-8|w(f)rg_ysq8&H!D7YAyhPr+YlM>Y__AVz(O7tW&Rv)uN>O&&cAyUSJ#`f%U$)*`+-mF7O=s1-+&U12tv=x&qw)=@+rxAKJND%#HaPMhoezeBaJ16`*`klwC)ki+!W=T?(ytj>wH73HojgDDnrL%Ryqa>Bexz1B#}W#-+&C3LOJVpVz4W!MXY5Voc@8=&ghpNa)4K0cbU&y{<F@0MBy;Bs6UO${_*Kny(B%eJ?*zo-#z>yU@k{%>-FDp%vGk&#qQ8egp8C?=VrxvArrdp{pP@bgYZ32s4a#bYS=Pa%YY(U#3(7sdL99q2{*>MsP2(8XycsATYYJCMj|}JwMt^^h-R4N8_DVSO{C{<TPDjR12Yvyze~BS>dLR5oe!was6;Vpex@f)7DMcsmO-ut?I(YC%#&r;j)}0*X{{kQeUc=g2LfhXYN81RUaVk(eG?^aaSqqr$~a@PP=c!Bvm!e`?;gD(Ft-}OOj5gVBy>~_D}%?iB;Nf(4t#EpvTYw+s9zVba~i<Z+UW1(x0zP;UiZ7Y3$6rpfIj6D-w5T$*{Bqcz7IS277>uR%-G!vM=Bz~k2K-^zKdDlZee_Scj>AIeb-`8iF}J38pnu@3O0v+8&8y;>Nhkm4y?3@6fR&Sdf~#g^{KMiIRyNxq(k$=sRgC{;t2C%Q)!93n+4Hg7r>-CAX%+e!ppV*>Te;>e8WAZe*O=!!va9G@X@5?);>_g$lf;$o^69E6@8UPkT~LYbChV=>FL6!a-KZkm(4^#!2OQgzBe2w320`PiO5!}RiboJlS0(s8mX|J+Kjx624z5~VTa~hh+lSV`h{|f`nyz*l2bW_<qBlLEDoaqc<;@tKyS_w{@h0_Kr~HFpwAtP=@<7)eT9eYJ7xn8X)H?!o+#ewO>A}X`FY^C$_`N+7)mU)q!}(=J=6e+TNOX}O=gYhoGcbKa7!C*1G#a8pZ_xhv_fPvdG7Z1cBiK^EJh_=IA@4`pUHbQz#WGC>I{EkikY`T$V)>5z-snd)RM0kj&BOqe&3RwwutV$@5TN3EF+tydOON-kuH78XWe?75LcuM#dma?{-O1@6xJ8_PHx-kdWyrnIFyYM!~^&0%tj}eSN-tgr+8+RSI{|d<T@a1l694{*Gip)O!X5Rmw^Nt?ms?*Q`R#8Os-sucmSVfBiD<*zfNw>Y_Z2a7rH)HR$#rba7nP0BEe0@0REmJARQFNa<MU@!IUCxN-jPh=7av<?m3I(`e}j9+gtHxZhNa73wFjh2$&PlY%HLkHfJ%$%<YDy2?p!o-+JS<;C3QAO&{j6w=l*zMBD?H*_UL<>y?VQz(m0-sjuR2U;z{DnT6kD)h3idB)R}^;?kR=J4%Y(y$E`J+5p4pprQ7161uaV3tu;5!O?QE5ZcvVTsblb<bIRY9KVfHq%t)~Vu{(^NWqsZE<VDSYEQ(5G><b~`-FzTY(k@#`Utr8fmnVThlM1a@v>KJw)<60dm*J7;&0S+J4{V1-?&GoXY+Jfsz`~WyUj9-x>=b7&jfK`B^!D`TaMOcK<)7f0DH_nf4imlag<d4&b5{l*yn_!$a!{>8&;!+OB)fd6&?pdiJfJ468kPve*Myg69aCnp!UOMkh;S8>6o4+*8RqYUg5Ra=O=55Q-_-^hr83Hiw3mD4hjFU)}^KpQ8)+wEB*rTv+B-xx%CSMaGC1PnsA?_1PzO$by|KSuXdR{*7O<&LYW0$9ljm8daLP;$OV~)0Vz!MMyt^NX?ICIrpH<z*v=h2Tj4r_Zv-u_-wOrgw7hem98G8dAJuSTs280Z5a@AhTg?chemMzknwBBc$aqK)5kb4kA;#Sz?3IA6z=SYe_{l3!a0U)z7vqaxIcTMJF0sup7LLk3E8^Cc>&-XzC)Kqw*b~n7hVgbhgvT?-$27;A*1e5liyhCA%w_+V@yqWU9|7X2h1;6TmaPRH9bmOA(W+bCP%0t@B!=+#Z$P<%B!IwRCKDH!NePV;Q_U=XP+wf|bJ0*y5@kA|%(@osivzH;$`W?Ru-o2hV^96OkMhJne_Z3_7ICDcHy&!ocKXqdBYL2gNlK(z$oPy{QWOo;rJb9jcZDIlcc7uZR7;nYRd1R8CX$=F^OC&uT<Gl&o1gN|o<T2YJm=`){eL9^W(01{G)_NbTVBoA7S0$iORmpu-xa4{9{ga8_cunbV76^|D~z#P=(a_7=H-5i?{t|2K-~EM5FT&J@QlIhR(v^@B`ioGHq1<5Pcq2H`f9IQq3CWv&3!QvdG2{^s1AK90PKS&eFG)Pr%Il)lIRq~@b56RUHpF~5f;znp)*m0PV9WFJ`9hMO6LREXZWx@`BNQ=Ylev_vV%R-k*E?jK$F&s51&gp4=9Z{iz9WavJ&(4mW?UZrI>F(OH+o;Q^eeWfs-*E){K`sC|~9{tFchs_JQtsA)Qwz`ZT72{M3W90xa<#h`Vmuu>p@VtAY8giED5Ch+lm(7OH)e(`lNw=SBTL6{!(m0C1)P9wx0R7*o(KBBqU~UK&?@#nCHex=h!-yuI5g*Uu6i?v>K)rV=6^lJm`4b~VFZEtC-^FG*lH=Ox9IM-dA4cP=)q(R~*3h7p4qkyC<p1lAc?#w71AiLX`E{lY^^v;tniwfY7;=*;be$iUFJAA3mt48}`9Sk8eEY?aSHzMSI1nba=Y3gHuwW+`^G?V9Xxtc~^e@0)3$<SQbds!9rOHsYV4CoCM@4%gGO$DBI-^yTnmA+~EeDdNiaD<d9iPC-m&ZNZ|5n1nk`j43&J6#ot(h7mE5jLDDV8mBoQD>&Dm8HsBAgSmunv;Sx{6J{K|29#aKjv^P$;83V87UM9Z)tH?}Vsf0mCwniH;6<2~J1=}|-_V$8@ccDfVq|MqR)hYt0({l%&0^TYpjdPWEvk~&F>^>iHMaC{-tJ=vQC~;Q+?*~a)hOGyzeQhXOk5o^j^y2j7MOy_1)h2%^%v!6W_3_<mMOhj{n{sV9Uvs`B`7-;XkuuiF2gbU-`<1$Yya<k8+Bb*8gI7!hhqtwF*(Hc`BcBdfQiXS*9UM}HalfQ>h*k|o*pHr({k7fxy*<F%#8l+LjSIdr7HZG95%EYRhi0JxZq#-V3>raqwc|kpYkWQ+-0F9LsqnxuD*~9*#(!4tu@7Ns%>I)Dlw5@Jc-a)43=kog9J;nT5RngLM9q+yUCuYeAI0w((G&rR}1Q?_kvQxe*ylzJ_RVA`&he*TgVWL$I+?+N4hCRU{q}DIA5W~%J=DVzOCCBSpi|!Cp-3poY~Epra>V>4#B^%|FcREF;(9sB)Aj(U4Y^%9@QT0n6XnJ&xLZ+Tp;q}-&c7t`9yW~FkaDzZ<rJ9w3M0%@<I@#wXzV0M=YCU1|{v!){bm|u{>8ZVU|Hz5=&hqs<xJY1Bvf9xvE`O7Yf!7ItK0ksL&!tYiPJHw*#F7Wo_R^o;(fpOjKeAAjiR>sEl(-GLLrb|17y@Nj%pm?Y?rlA)Bn(lT0+Y%s61q=jB1;y=_lH0d8VXGvcOlzigY#p2G<OH%l5NfW|n`(sK?4SNL2XKCv=I;}@dVCNL8tDm%MiILo{W|K@CP=AEE{wwf&TcjyZIRC8cN{c@ECTHsN)dwMb=7s=lu9Rkn!J2oncWu=NAPGhnGHfIVp6&`+7Nh98}tT-Z{WTr^P8H3C>0)N>Os~2L=az>`qm<!SJ@UFA&hkIW4dY|mvBMq%n9sfWq2wG_x)AB60pQZt9(GMRW!1b{Yxte69E{OFJry+dn@G+mKb9M5h{-O(w*qo|yDg;!m49MLR8GD8G#xn`Zop3HfwLq<;vLJAK`bak0IN~|xn-JWlA%)pY7&_F`Ij~YYyBmi^i-+vICFM+N@ZgL>a5B~;9KTFr{*5%KCyKhVw>cOHM?fLXQ*XAh<*KZl)Z&HX5*r2dMJRn~pDX+>hFx%X-8GD>bG+162H8jO=AJr}obu1(2oteDTe0SBtPb8*jY!HLjkQy5LY~~k=zmySAkZ?hyOm-8BSl7%Fw&)co^GE^R9AKTeISMbcq+PHqq<6a1M%H`)T>X72znh$uxRXncS70gczw5iBpenTw1=|~hYI4-%=cw`Tq6N?*}vlZ$-Vyap1;(lwCxBkQ*70J{39*eNI6&1P!Ra7S+&buYwRxk>*rE#7HM&bb=tD%tR3@vXP%miPU4p|ypHaTn4clnU$X`jCSOO%_u#M?rt(_Q-nF`Xt=p|3yR3X8jK;gePpwums6+)(f?oJE8a1)v9SmQk91L%>{`1KTTJtULVTf0>T>y!bSqy^6buB2V=|EI~sR=|?;<NQz1u(*nEq7GE^d<Nn%xdzsxfPLZ*1g=#^uQ&_vQ*Vf!0=#$gITX?J&NFEF{XX_UKQ!1rjbYDA{YC&RowZd`k0ak((0K_BKY~zua}-A&yu&P5<wrCR^d0DEHHVF=L#eWveedK9s}hS@w>7|WQ0QbQ5m)88bP!cct6FVCodBrxMalWXXjys6lNPim-e`b9+b_3K*^cQ?D3%#-CclD>qi9@o$v3CJ2vG0&Z_)F7GY(Bgk-?6;M=+vmLBT1lb$~K%+o{+9^$8tMvAC)v*I85l&ew>3H8o3!zN7~`g13D7p1V080k~E*Au~y)Q_Zmsf?yhUymio)*Irgae2GTsWx+9#lW4L!KwFQq&eD(b8+$MxCiK%xN?U4-|Jk&1w_ipe$cVWKp#qECjtyC2I2}Bf4<6Meim7-7+JBeSe{y?RWN(qH#n})3Xu40>2LLVt3WjN#efCPjj5cnvE|QXOjDiFPeXkh?c%Js^Q}PsuUIEv5&T#57H?$XJBUpsxZE^SKDY8_4oFNG?@Y^`wlJFXKWp3pA%M|TeN=m)$J}>R`W`4pFZg+WmH#olz^hs)Wg1a1ZQRLNX=27A<?ZA;*F*vSCwCr0s1St}ojYyA9%Xhr)nSi)$ZA}QRGHqqh_6BTEhO2t!+<6wxESV91tYA_qel-;tY;leKcWTJ)4R_&;0N40Z;SbBcZ&fuSueFNqWfN8!^6-W+|x;(Wp0CZ7&@v1<{QOXc|wt9@~6&57akDEr|8f7F{KTQg*j$M1wWU!sXf+sZ>r(b<e^4;rNh7c3tSV}9$LsWg--6R6?(J7|9C*2X2ipV8Jl!hzgDQ=#4|FGlPIT-o3$FMdsC5B15jnjVn}$zxAO{2ffU2biGZRi(ZKr`8+j(mojEckhTXKxA-$}YRJz-XNF{_dE_*5A>2|gHsN_D|{&-%B@+I_VZPCCvWT!)|s11m8^$R1sH^=^~{UZE+%b1Rm1a2Gb3win$OYXZ?Z)gCUyi1AaM2ArkEFIT=aa+L$R6OlRt<xY59%&7iSw|uOUcy$DrFoCQ00G470?d~lp!l2MvBYQl0ssI200dcD'
)
