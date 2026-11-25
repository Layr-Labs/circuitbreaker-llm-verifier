"""Active -699/2000 first-pass, gain join and starting-majority proof."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial, isqrt, prod

THRESHOLD = F(-699,2000)

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
        for k in range(1,n+1):
            for face in combinations(range(n),k):
                C=[[A[i][j] for j in face] for i in face];J,_=invdet(C)
                for sign in (-1,1):
                    v=[sign*b[i] for i in face];z=[sum(J[i][j]*v[j] for j in range(k)) for i in range(k)];count+=1
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
            for j in range(i+1,n):
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
    ck(U5<=F(73670242301,2**39) and gain>0 and bonus>0 and 0<cost<=F(123275637452,10**12),'active energy bound')
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
    'Jq}B_F>)>X#SAV7VF<-u{n$YX{pGGFsNsXY`V4IE2J_T%fOHxsAujmS!pF!s*a-)poo|l=9o&AAyzdVEJ(&ik&WfSJxFF|RwHYg-+Tf%$R)+l-M7g%=A4V`|yK_{P$$m7-G_sS>Q!xAOCO_L$Q0jxAVdFAl8yK*UWr_76gee4=IF*H9i)JezQMOlM1aT)}h9{5(TynoovAo*z&NCq_9REc5Mo{Fa+WNSG+~oaLcrT6Q|Eh3}'
    'X>Q78&{^&_Jf0wU2J17R;7A$_rv2tIhlpcD2RY?0e+i?{nxdCo3?SH;Xo*9K{EN@RBayp=0v<;X?q%DYsa;;Rwwu+=IAUZ37=tyD*6ZX>`g~pX{uT&(8&B=`PvRSJQ;u(f-@>~;RAxvd<~pC;O6?9oCHe;=4l+`@lw^IlhO8{%N!QMtfdMAuYAAT2IGL00c)@M<q54@=y&9}EUGu~t&<9u$7i*118nMj$wICf+I9I|e<!G>z'
    '-#<mUaEt)NH9jO-3+@R^<{owg_E!Tsd>hzz>9~t%u6x<W7}Gk2o>p!<yNXNNxb35vRshi1F-YSTr%gT`gSgSyxaM!u^`+gPwJehio2HP{x;D9LQaQ2HN}Ol^*e<y9)`or>L!dOXGkx68yNHiCkdb7TsKj^tOk*>#X4gtyMYsvz878Ky&Lx-)R?{}NJQLNyD773#8}upVQ9i=u+cD@krugAXx9V#`k`Wuhr5&@J^-=v+YDHw|'
    'k94r3q80razn7UL(FBpUR_tI1?@H?%`c@+Sd@z+K<-m$_)~BM>ejoSyM&CjL3_#^G%OTQgMO!f@c8B#K8>ZqMBB}>-%XQQ=e>ICl%LwFIp6?=(>VGbxUk9WTj92X&U5Yt=c@`#DrdwjQ_uZUUo={2TU-6f0xEu2ji@ak1dKARB&Gi_5BeqMoYOAVOaNLMx7R0D0+aTK9G#viT))URh_6LT_ODnS6lk#}sww(>;*^E_=KsNA&'
    '(PhuBQ<KqnG`aOp$8&Z9sOZIVz%vztrb-eA)&3KWazMh>iP=K9aX0KMTheGLY{|$!2Tut~>>*GyCd6tXX@`pH?CQr&=5nVd+wVI%lwoxm9M$q(xcaus8QGwEVvY1$Bk7|qCEb}ICi3e9nvpw@KFlzK$sle|`biD2QM3*|8VU@)e-h*jDR^1JVXDVe_l|M+G^&5B-tZ|;c*?#+`6m|6t#b?MPHuJfgiQfQ6fiosU`#h6@_9K*'
    '67ReXopaK?C!M0h2(TbOxrFpejZsd`rfnQr1&&4eD<kieENYyl=Jt0>4KbA%)$#BE*Toi|rSKt>9-EsjZNY68>mmrSz@Nm&ZG;{Gbs*KxKJqueSjIB;aT@WH36`uzrO@y=EP18Qw*BuuH<{Haa7k+fJ?{wTpNU>!KO4#c>bO^OK0q!8OJ~kOGrjsqyEx{GaKt)kJhMX|3kcq~WQjpMDVT%3M*ZVRH9?U_WX|OS{|POqedD~a'
    '1A_jPI87<{8Gy-Le}8f>44G6N!<(p2!JS61FB8keCKyTnF;uo3zI|I*`(VH{gOi|9Tz&I(u6{OGKST^{e!i5&Wk9oJ7jS)%;ceUEh(80B{mbQ!WwG+&n;2yXu%4R)?(e}!Ud*^T#JcYND&Ri%MFTc_FM2{Ej{6^Szj+PlKI*%&L!L^Aq)O$~Lhd)E%WUD4qexNL@!0!+9@(#43n#PYjNdF0YW9-?JOWeoZ3<6hhqFiG8h5iD'
    '-T}s^@bro?1@KG#&<C9Wh*NUzmi$%xqWRWfU-X>(_JD9>MC!|<TyQuqpm<opC@w&612!oV1&6?=`8sCU_5jZv>2PpXS(;fTT8MMN0M<tqt^01*G&T$e1xt>30%v2QkEVX<Uf;sLMX8H*Nh@(&{ncBO*E}wP1F8%lD`1fL{RVnuotl7>>gnvR8Yt1)({~><nYXe`Q4xsU#{u<2M`U|gXTKc@AwC~oRuJcI##8giC?UL>vZzJb'
    's(hbLtC=)`<D9xDXqrk^$=%_7nD};BM1!X_&k1GQhhIl$!Sg%|k&ekt*xg@LZRyLYCjy3+22wgY<fVbarjjzJY%Wg|+rb>NfuEkrwWIA@;#jZ)SLoP)0jkvvGlkd>2$J>{u>^ZtVQ^`gO|g7-=>mMd@nFI4r2lAcoZ+UxP)3)II<&2R>>(FS#`Bm{`#YZ54Q*r*HNGddHUQ3kKE#BXGxz*`Sw}OQw*?C1A~3=t=k9wk4eUp9'
    'RG8<@fV?$ARDwf}W<{urxwhPO#7)(iwNMLWVUzjLUcn}pUn~>fZ5E4l>h!i`0ZE4$e?SAQg0CsE|8ZrO2hJ-=;8L!mB_n3FxOJew)`%T*+cJwVh$HLW0VAU$pdtt3xwE|IzF0(we+xP<KqftEG#o?|g0+y^hRLyLhHS9JNU;vcNbzmhKzLAG76Q_uXEE$S`uO~dHkg_%P3g1WEYqRi$ME7gn1o~twgV!5O0y>W#XR1rXT=*?'
    'xTJCk+oLm3^sBSdf%e&xH#gN_`vTY_o2lv*{EJkJg#(8yO5dj(!MSXFtRP4GbzK!tU^@RELV@=ys$*8_L6#BjUr?WHVRrEywmuRbK`C12fk&i$oS<fBb+0LD`i8N$dt;-+2Wf3BHc+}@cHni09Hta%co5jIV3@f@RWMU9%9@EFS^~kE?UBN$F_$bXmy9<t@lp+<<Y_T?&x9Xg^L1>x>2|8Er)IQoisAyqoag&LlFZC12b8lz'
    '9>L_IK$6y^z^UNu<+J}$24BM>U}+N){vQ??PtfQl(yTbT8lRX|r;Q`MnklrWza(200R+M1-fNA%p5|3FPW<Y#k><2%=Ek4?k7>IK_9u<EE^JtDHd)bnCO+q!#zpKj0M!E&Rf^G$1V&?MRoVGSct3qy-=$8S_=n5>hf2JEid?_55`SbZ^cA?D8f===dTx@6Ri;djB%UBSVQu}T@NnBe$2*g4myGZK&KB}ivf{9|%M$LNfcC|Q'
    'pc6FY+3K^py4ZcJPeFKY>t&Mi_t0_YoT)?z1w&K$dKK{pI)8Wqo|q!0mo~Hc0wQbv{cWI)+rz&51|YV^a_{Q9b03Tyu3@IsZhB@lrgN50);*^wHIe?(z@PG)?j=ADxCG}+c0Qq2K!~Ye-xB!0ZZuc9puCm>6_)w~p476lOOlHjp+4@Kctx?>u217<g_f)B*@+5A&&|%%7al@46mR-msfxZ&PP9bNf?X&eQEd_s)_5&wJc-5r'
    'KpPcE@ndHZgGF29AO?g*^hbGLq+r*0I_frg*Lc4%pbXvXh+YAEgQ8|_EtSW?iG~0?j)KFNE>=^b&2oNboI?gZ69{2XcaziZzW@^V!5_@F0Qo!^=UyW-&DZa=Osd#{!!#4Hx3LOZ004yD$x#mm&`TqJvT7Z%w73?>%f^;8A+~70{T^IWID8_<I&Ate)w`Pfb*$zjSEUZ$CKW$Ll(T4}PJg59#Aa-lw%S(mW-G{-&zi=LJ}xwq'
    '4ei9sMxUh2MzR)H&<WbJ05nU1sgK6>kVqR;4NH0|TGpr92lyc86k981&>qO=N<p9YYf5#F;(=2+-h!tKN3>@Dy#c9kVA>fWw?zHWl;rFDon8(`qmFvHD$x3JaI)2u0}t1^owjEJ`;)N7$`iiEh1O{x;J&bfZU-a@mgTi!Jp3Qz#r6>ohvst%gOB{HeU;cvI<};EpMzDY;^#D43xNuJR^bM@M9fNhVUtG}eD+q?nCUG|UVmu-'
    '!VeRXpl4UoZi<w=@C*drmR_9CiXyjl(o4u@0U$kmt$LI64|Z1t6PS6g^6NlwA{EGC+Y~<g&yPJk#ckhY_A^$?i-dKTxSr1`HZehljoV_7tJUEeZ975_9m`^o40+(+!qDjIIJbt32_z{z0#(?_v%|}hvlJ@L^Tf2xXN(Di&f1df?2WSI*m}hiG;=ZYUXC6g{i=km|6ElQaqavMo3jQEA<r5Wtl0}5oHfZ`*cC<zt{-U9_vnT!'
    'aDEj8Y*`T^s0;h<Wr5Fe%;O?eNLJFEKOl%~LQW*B9d_MP;?qtS4HY}$fMnU9umq{u8+U$-AKjqk;>DA;q4t&o*x02vadmhF*Nx)kB0I+(N`4ATlLi6lUzFs9apf*v&_3p8{;*I(ASfJ&pWZ~k*5`=rr4mdsb?Y&lYW-vKS<$!qGE2vuWS>k-gj{n#{ZuR81ZI~VmB`<7%^#@D3$0juhKuF<G`-lLjxA5RnjIa@S$+m1%?j^C'
    'H4L;LesoFR{KpmP0y#j7s=#H~ypzH$K4;k<X#E6q8JBqgtxG;9p$G=vp<fDUGt02TA6jyqNmK_xG5!+*(j&BBQNXSpPSvshF}+gaoBSBd*Ya=OrYE&tX(X)4r;3Yi<zYbxEMq}!CU-fF@P^p>rz;@_eofoPXg@yAzU=4cJ57N<C3j|^8s5D4nzolKUer;9yfTKm%Tx7FbO%uoO^2BlLX{h78(LY=!L`zvFu<b)bg%^K);Eh<'
    'N138QdTmGNlOeQGK<}#AtgTZov2ekg#vJ`;;bkd($Y_IrFYAIMNzRb_H*MNqmMftX#Zkqs!f%NKMUl7Atiy~fC2U*&RnS?J*M%b>YtEyffUthl)0gSh>*!YJAaX#;ryN6rhuMOyHv$~ackm5ydUtX_7gOkF0%dSjH*hf4{>vJ&EU4+mJV$(nvQOpoD`HMz<d_JWG*Qrgd*j}4bWuE<e)w_ML*RWg@=zKSfJB*90vuhj(-+3o'
    'z{u~7UfpEXI$37K9i=bQmYCav>~483Z-nZxUK1HM-NuJ0%{`-9&!psWbJHc>c;=&h<o{uAPl{XG3PO`<e{i#G=}>d)mNoI&gXZp3ipvSWaMa|-MowiHagTO>H_pEH;R&@_*u8Qyfk=o+q3Ye@<;McI5l2%2Uf8bBos6%bi3|F6fj%dL7sPYIQLyxYnpKCL7`U-;Y@+=x+U}~16$g0NJ8db^YVvhrWK(`Rg*lG@l7jE*Yp*Yc'
    'LR!XDQkuvrl<A&M>(5h!v@s_qxwvgJ8Z=gMB6lx`_yWB$;=q)2IVYIy#3?hnpr?8+^iFu$*wI$LuqEFY<oq|x(T_*oyF8RDLu>`>2*K($t37;-WmJ7GZD?R|>VP7ndcMvw*oFHRdB~awhWMxR(^6YGhddT?p_OPy;BTCN((WFmK>+w1iCUa4Ml~1;@gyg<Mnw3DV!S{BZLZa60Q2&3PB<EKf>6~R;pdhCGK&TAn^<QjN5x{H'
    '*URnlOfU18ph4ydx%N72&eZg0*Dp_qz-oLMMDhGpPJ&w9TRsw2N&18UQ6Ar`II=T=^Da4~)m86Lz-H)hyY5wEBv8g|^r!s9fG*N=)&Y-hy368dyt~e@CZs)yKx^eJzV@B^!r_P!$etwI-|g?Q)4OaD8EXtlGoNB0^XYaNKV!x2(BF)49?75Z=hgRrzFTFh$iVsrfN?w!V@X(T?LhD~#zX)ftGYd7m;}qTPt+1zH0^)aa)1R$'
    '5%1*!q37Y1u%1K1KkBx2ms_`sEbKWsAA*O@NRH+L$DF>2K-53Ba_)A&ZG~_pE}Jd}fCE*Z`I2Z9TKf#F<1OLd%RLFNJ6GTK5W3|eYo_=4EcA2_Ldk26bVB`=Av__h=d<4Z!w$6LF^;s3^br^5Z@kdEyl`n967tw9m2Z>8ijUMRHo~8~aw`W%H%Bg2{2(tGYN<qYHq`Ir19Qn$L96d#_5eX@U#dDk-~y*vYt}LgGFlhrNSJR3'
    'xBCDrc}tFWNBauES~XX(3o0m>p)u4{Zk-INaL3E9`}+i3!SW|V@deO}bgl{x(U(X*1T)G_AvY4uYMZhIZbpvg7PKRSA|@405AKzPp*o^e<3^k&kb$+!GSBIbzK7<3Q(Uti93In?%j?Wz<A)XJN@v7lNv*oN5XY~grrYBuZ=@1+t#PukKJK}|q@4eLw+pnpdv#5fV^-YI7j+c?vEf~y6Gn3U8;vn4n;H~lP<r%u31VdbdwQ~M'
    'b^jR_5o4I#ObTXSauY}WTv#|!pdMQ12j)-nYH4LvV8h@zHv@7aFQJQNS^QRQ1UhRzFHB%v^$=KfJEc7{o6*I>x@BDc+h=#M62$!6Hv%6o_3c!z%1?6o>=$7>wxWu-7yAQn|JeaLobL%5kBe0xG&M46u>dVYcp3bTPTLXhD^$eTX4SR5Cf$dlBXX_iM6_#nwjPfi8laTLN*8FZD}A;=2}FMxySb*uv|}|EJTqJAC^taEnHz3c'
    'Q|1E5l*c}XJ1%rnI^VFM+!A(sipK%ddU%$yA)5AXXkw~a*k<kW_!>dujEb+Jw_g`P0k*l#dUOpsV~l#xbZYYz1K2>Y5RE@Zeohlh(N-vI*l`9UCt*vC{Z@|PiZK`~QBP7SH%Y$Y>=VZEb?udRynhH8J<B?YuqD>#U`9v-TzZnq{{$Csc`bWW0mKu#!X-Ak;G@Q3p&zAQ-ItsYS^@MW);O5~WQJ7f*2T8G3qpBbRTGdvX?%Bi'
    'CC=E%c0x5Q%!Lwt9uz%<^o}8ZPqGgdF;{N?I%H{H*y%Q2HQ;+#iH53hIBxGZjMOaMuc;@wRIYMKk46YC?kpfFm8iYZ3VP%YBQ8vMZ+$MH%~E=Y7glG7#skBBh`E2wfx(Sg=wmjusEK7dxCKZeYwiS+>viQSKayR-eP)9ueWJiOIRn1LF8o$?j7yOs>ON4#;&U2qz*~lh)#l-m%4~0@0!TL7+KpcYQRNL_)5@<PAlvG$wOMe5'
    'VyrtQc0oh~ah$s5#Z+nmj10~p5xf=JHB3!>-2C!h_hO_j#No9)qK7cng9(<Jdx6ku4F!ul6{RQ0e&7EgVRHn*rQAZMFKK6m`>CW<-m2l1OXs2<f#6%ch20Y|R}Q<H9V>QO3MaNnlhvK<!Eb!&Clg?ZNjm_H02;1`00%`DPX8k8tUx<xG>NKZvcI&A^3`0b6qH(g{Tt)WVeXME4FJ$HSXn8&FH8}mbt8pp<5#?_{1e~zD%TbH'
    ';}xR)k`oLtDa@jSTb@H~j9G>qi=>8icGnzK8OrTixd=A6*Wiu4GFMD<vWic<uTlv+QpWPiQ-sB_2jq$7lTzwxBnO;DveX#<g<})a>=d2ar*36Ud0K=!gW@mRes?A7$)NpDh){8WX%fS^NMop%KgiY$f*4)i_byEmA<*(7DHK<4^uL}afUn4_K{g`<K@1aA9C-^wf{^p0gni#8M)*mqEIlU^#cmRDfGblo^d<jg9kYX607tVF'
    '<OQ?ca~6GTQH~`iVxm{SH!5frkIl)V_!f6-oo)vS<W!pipk|p*Nu9}80MvThmeU%_4-GW*4e^&HY9R*8qSesq{eRV97A0SEkwa-WZYMcR@dtg&pv^33*{<P%VJD1Ym2!WrN2&FjM!vYtyhvV6`Jc`4tda@@pnccIm=|Wj`SoUb&v6@S*DNyfYot8KSi6&0k=P~W$5Y^|xqqv|%f2q8Kgqfa#MX^mt6@q?=Ds?Bv>`6lTQEs|'
    'f-%C%%Y(aoTvk`6ao)KXOg%I@9&eY~7cQ^dWSVL;?jTSxQYY())<_42&;NkI__5F#bY}W5-lJ8yDRRQCF4Jp7&6auM+LU3a0!?j-J-XKlaLpdDjb6UJ6<Bx+2J)<#RbiC<GA^A!u6_H9rYc>2%ohI7ya-r;sbxk~I{I@+_GvnTV*y^^{8=Yctzhh9g**&JM{m_59SccCk(cycr<gKR`vag*q{`o`RfMN;sVFQ}eKGa$FT(R)'
    't`<?iVpaMr8b6iuAE%HhaXWcA;o2^}ViGtlu#s^CYKy_v{c1jHG*rk|2F0Awp`ZFQzJa!w<y#Fu9hg$uEN^nr4uq#OJQvX6j5)F>VB4rwipIR}6FO5XoIwtx!Fy&yeI>rNjaP>^R5UUF4pZL_`7fH&m=k)$2!xavPv=p&vR(k=F4S+6%{JmyFE6!uvO)Q9$LkGa?BcA4j8<z96b<(}9o38s7woGlO}s7dDSUgf#B=4YSA4;s'
    '=(WuTrSR|l!GAM2?bTWAP^v18Ycm)|oxWVtj=9;JerFt(S^;erv@K4E@v7=(DRy^}bJ~N>+pK$CjB~~&MiFU)J)3d`cq$xM`M?JLsNP0%+CR&3v$|Np151F^5->nTVmKJg!Ln2vY#)dveWbcP1N-K1hAC-nobDY)9>O?mJ|B-4y|`u&=>HTXbKVmN54<2)Z!4pzhJQ0p0<WMuOO6Gby;TmG?mUJsr_ed`^)62gA*%(9+nPQM'
    '?jq1g=;`$!oeo%Ny8SAz(0YYNAd@x<08oF3V!e8Gmb-YE!~Z6Ca1#m9%R?Pp?=Dtw1Q3Bw_mkDmL=S*_!U+O;gE#D=MUvF*)sA^Wa=c7KOS1^0ZS7|8ub=>I7O1Pyf({?&4`aano}xBB3ue&%0IV~+LtC8}_~7cTWlNmCCqD#V36#|ZHu{@BU0yPFAVB#WAG4D|6`Qic?`sK6gpMkmB)bPS=pEm1_>8+1nH-<#wnl5-qI!;I'
    'On7^)I5+{!I|)F^0k39?)a*tSVb<%wYt2qR2&iZCdUiDSEHD&ug%pSie!g4eh){SPBIM7^+Ypc)B`fr%+iP>?0LdWZ8!=Dj?qt3a&fQk?2t`4*nPQGy?JxzBrC`u3Iy52w!nBV1!u(Jk&0jv@iyV(p%55BDHc7qngT#CSm9XbeWmUPeVD8sECik>hb4!81TwJ-*iz7P&GT-xk)0fA%{sd1W-w>-Uq!J9Lnzh?`?l?m7mQcS$'
    'GetT7Ss7<VMFPNbV&eC!;=jO?Qnsv(^>~{4)~G?|RrnfiogXO1c?TKA$sg+|p_eq-?+%i8%mz6JL4$zPGujF8I#BktT*LgL9Xz+9-r27E?tH+l-agkXgHm<vzmrN6M_U@kL`DRvj-&hdmqK5(JW=0VBL?O7Ybt=b;}pl_y~`z*M{VO+<P{dRw6ksAvFIUk3y}^@gK53*W%K^+@qh82j5Di)>~1L>;tb9!#=l30rUowN<?xN<'
    'Hv$z}8xnHA7G58`&QistCp9-c6h-HEf5t3jbOJh<K_VtX#(>DL3$2+@HpoZHU`e@lqn!y|dYfH4n=!hXVuW=Cva>C02xjnO6S6F!!{dO%(H)^L%+?Sd^rJy{iKh?~DP8VDhFKkRUgZECc>q;nnHTT15`?v_Lz{QMpsGMIe%0Ecd+e=sohTv)mmh3?nvbh#kT;0#4bj>j>07KBr_Zj3*slz;7e@J{*a%!D+sps4jJ;*WTk-_R'
    '#W!0evpPfZT-h3~cuQ5&tN!*R*CK-SX;|gS1*~Qbxk*Pai;OK>LI2L$uSu3@woQDrxMC`2XMdrTE+7$kpds+tbNPD!L;)sHHn&7)0{vfm3f+CbsBpT_qWM@RnFHx^)<mr8G}^Pw)fA~?>n{43_HLF+5XIm>B35>dLLFO8JUtLur6odnE3y&h^T}F?97D2XE`4?|kEChF>>;hQ330CD!0URWX;~<@e!6$oWe$D;);L-HyUi$8'
    'Sti(M)MBSD=gq<6Cx|$zgFxN|%8Za|1lj7UxRV`5^)`NJ_|qBzGKGc5I!|C9c0&MBE4mij*%M9sC8p4l(YB||+)mijZLGcMl4@@wR!gviZZ#Wv6ZA&Qg`@0*N`=RReJf%1;aZG6DdcQ29i%MOtv2X=3xG^{pPfhrpZ*5g2zyDacoZjgQ;1LI5lkFI#I>hczX+!HwUrThZMO_#dZ(9Vn9Lwq6eBDhv<C3A94)n36l!`5c2CnR'
    'RTk<Hr}<};l_;N$Jn80NeyZnzZ3({w_GR}8S!@Q`Q;I5r6>yuF92TUZ%5Kl-4h?h0km%UqfIW0ET0e-K<CVmCW|m1v+z0S57YGdkLadVwuR-$>rII$D;Rpr7Mi<#iF{Q&D=RdIOl}oJ7?3N<qjAFMRu5V*-2iQqUAS@~L4beT}-Xo_JYnXAr$|_Q{uxeEGQtZ%_*t6AEs{`xYR3cj6H_`_aKUdQ@-Ek91xq8~9OY*kgRGkK4'
    '_6#D^{o&Ai_&eYZkb9P9AiFPSil4j#VM(8jNqIkiC|?@H>JCkEv0c<Ldd~~2BHq#(2Su?-aQXSJYl69fsPj~=qhfyu79v2ejsajLBvf%rZ=-Z@ra8>AUGD#HJm1KE5E=-LFmvVjd&h0?`2^o`IB*AEqZxqjpKKU5G2(yc5moVnd?^kB-n6L2ES-Bc<H3h>iD}d&M1BAR<U3WPr8P#^v7L+=D7g)bGIYP%4-k*J;%T~o<8M$i'
    '6X|Ls6483CPT}y<Dp@UM#o&FVcK3Z$2+;aYKc>cLJmc;3Xx66+$3b_YLb_X&4cO}z*6Ql=3p?Mz#JhqHX(z}dY3%Vr(@~k-<Tttx{VaU59+vXS=<P_6Ry`Gr9bG&htJ61YSan4vjK}`Dmg`(F#@)pZlU<k#+@^?4d5rbW$~+rM6ul^_SfoEN06xq|h*WbG#Jf_R;YYC<I;F`jHSh{S<kLlCmi#*9=k*okP}2*0T3VgW+3hA='
    '(d%w`@7R;7$yX<Kk=sm9m+&1}e=RluY*u$sZ7APqd$@iQ!;DU3e#!%v7}rjJ2khf)DOqUTizYaJ&5;YfE=FOv<}uI`3=gwLNj8wqQxG12sLEhF!fnw+D;&zmNkbm|K$u@m1;j*-mWsSncgnyPGm)G|BKv3aqe@1MFUWI+SP~`Gg>mi`IdAiXObvlP>AQo_%>MT!YrNzwX@RO<OQvEjG`09+t!5bH4+)H51|qu-VNQ%Rh-fC0'
    'qQaUMmtN4Z&)1Y=1<<y_7#3^U*A8t2XQx9w_<Zb;VcxPj$~zL_*rXgW{bY-^mJa?Kx-ZJSi%1=4UeZMc1f}fGn9CU^;9Upc_Sac|A~+ma#|WE(ajIHL%}P}3Jzb4?aXvI<cnmi?^Qa6J$&DGM%L(;2?6ySaHRf2`$_BiFcoTNVns_V?TYEQ8Ari_E&dmw>FRJkhoBL62oo_k&G_(90`_|(d@HK@{=_WQi8Q`3@v?dx(a7%5z'
    'N|n-0cBEXe!;m$j*?OHn?uGw-h`ce%k~OI4VGPpP9~PU)MA-NX_~vjVv|%&JxJE|A+5|-6Nf+@DqAB5L-(TP$!F*#=H|Y+ps4c>Ly<-z%7%If^wD-Hj{=KMxxtie1NW_Zc;CMw)eFFXv{Ch@5bHrD$R`XZ5Y=(62(hL@N>y{qnkn!C*S#+Er4Hkm7X??|gGVXjv=uiPT+CwfLU2V~rwbupIPYmC;FOA&V7!<OqV^VXFTmkQW'
    '5E5)-5Lovs`sr1gsx!>mK)1PL_pA>6tK^U9B-wvF9z3Dgg89aad*-)3c10&uq5ydr;Sg$%#;6>8nE)*9>NKA<Q;avv=0K5sf=580khwPStr~#UNr-=vJI!14qO!|m@nxH>OCT(z_LC6Y^)bCta*Y=AYT&8qRpWSO=s$pX8YlS2>B^YJEx^hV{gDur<VEF}TaNPRx<l_~JwuBebTNZS%nJm{KnomE_j4fmDMlF`2c!V%jWcn;'
    '8TKnmYUOr;>j0yH0e`sYe`03GFYc^tkRN{b8z3f@f9p?PIGacJ{=}cH>BnKu=z|GY>J&w8#TQBWiEQ5d?aMadK(7@~-a`6J2sNA1RcVdTR0p1PnnPZ=Yd!tfDdHouDQJVw`EbAlQh&6SaVSPz<&SQtP^t(WIiFZl)bST}Q;kyrZs$-P$-nK?5(3c|PXT~F)QECmVMzWBT4Zid%+!!s-NzkGDE(6~!y{bT3Zy7yw#1E#gArka'
    'GpGk@>h|0`Nv1<*{TUyZHFwH*6ysg<tV&DLRbZII?sM@OzlFn;%_eppcp#MH|Hn5-94q)|u$|)2;>v(ldc=V25r^U#oPfc0bmg}46;fU8)Uj)5bsPt*P6omyeSGQ)a~#4ee|$ch&hTebn(aTHEp?A0BJn1Yqk$+U>&L9Qe<PS_^d0)&@D$4QY6sq*<_mrmo?D&XGF7pl%A#CF@AKy*Y3_&8Guh>4lJ%RtlL6ZxLZ=xf`8vUe'
    '@~!O1CL}JMI45^ThBs`e9GGUh+Pv|dWn9YHX6yALX7PF-A97NhUy`RZX9}A16QznTt!1&y1U|IoEp3`U+7@M@i25Oilt!6US+`OTUe@bTo}sTx682!{NG{At{Q98$G1CuO_KJeT9`=G<K!`F=m9%#|XKMQCb|130F_;d1n#xz}W%q?3S%Jtz)J9(SlHX+F*7VGzD$g3z)dH$>o*V_;nMI=69Sopc!iE`mFZXQT7|`|tdnQG8'
    'DwHP4lywPiW(Fj(mY_qkz<Gnq5YQ?vFJfOb{0ercw?*+53q`uInAPNx+AlRS&=mVvRb&KCu;XU3{Cez*AUm+%dz2apG9=RA751&&6PTA$X5-1v`sX1SCuK+{OBvu98`RnYIE!;aASVOw>6q)TNO7=fv*<UwOF(a>6u*QxvmCZ(HxVlAYu<^(G>?-}0!}d{8oj1EUyEcmV$898<2Qv2ls@YvrR;n_uk+^$f2KwzVMPBv?wZc2'
    'aU|jqm!#mjQ$=>TU&h~3IqV|xCO;b5E!Jak7>J_UCT3ETfMJv%UF88RsyEd*<Q<y@q;!fV%@T`qC{xAgfFY%q1`Z@2!Vq6<gnV7S1m%j|SkQMM{@^tEEP-gNAvFW6kR_}%l;iIWP7Cy8OQv`omI#k|0-hD*m|{HaDAE37A@O%EdN-BA5b}{6j*>0iN&mulmr3MCI_g7z^~O?i+P2k^w-~_Pm~^raM^FSXR(!WD9^4Ug0P=EK'
    'WuN)_RcO4ZKR=csox$#9G^sRlGm+2Se1+kjZk7sz>@w)@C!6QEI&Uu@Y*kev2WFdT_Ang~9<&RhVcE3`a#%}abjoLO@C`4ZUoK6_kFQ;a1ze4l!pfr>&c8!b)*O}+hyI&?_<t`r6f!Z+Bv(dK<;-~Bm&JSlmLz_gUMCfB`0y;=e=m3qvCvzeY9VF<o25`UTi3@+fV2Mt{4!q7p5#s+n2QTtRb>lkr?mptQv%mt@0-#7r%Gft'
    '8xY|#My)v_zg1?(Hb)BM8>~EMUyFMv>Ddta2VH2&o+vp@4iYwMcBrkocp_~@sL*65iCz213OEyc|JMuZczX3ox4a~CL6HkG9!F2n@Co>s6S=XxU5W!3)?R0hxFv*(mwjr7;~h*<KOl11&99jW(rB4eKOP3%DDs>vDjfdntg*y~($R0T%g~I>4f6lQ)Y2QGWmjktb1Lll;m}4!cpA7qAT*4dFFHoJJ>J(#G-e6_!*~HK-ZaE1'
    'dnMS40vto(ShYoQ-Ox<rGmE(6^lEp0I#OMZ=E9ALIRiPmSOBS({?v}{rjm#rBjE$5d->x*<3N}ja^T3Wbvb8)B7ctx&Q2BaaK1{nOQ+W-=I8DpGlRm9>8*`Eo?}nuUX!Yc>`((?VgH!h{k$1H!E89V8=UkFxumtxYE6wQ;Gj=mt_%HEe55NLcX%xDCXzgn#5K*>QL9~Pb4<<2;H6^?d?#-)xT?#KOu8FzlyB0BT2z4FgOBlR'
    '!mxKc;w|Z)GmHaI`!Wy~gsfVYiY8B53ITA=WVaka-MH#;{tue-r=t9#gijqOToxv4n&4p8B26%uz5udYW!Vamys2DpIdtn|m&^y``MDf^Qj4ddVYN8I>PhIPtyScFW>@iARTZyaD9b4xv#(7LDJ6175A_ul_o?FyC^&sYzSZ9-OTU+*Zd-I-Nm(l$_2a~wFxO85CQozv@8!hRabe(h38J{*^s(QL)(JsJcMS*9uSo%#(bC+T'
    'y026<j5KG>4jGn}=ei336H*ly$76{{)h&DuZ}ogRn&55w2rRH$sdM}>rQzHu8k1y5FwVjp)KGptt?t|e&{{ii$n%~pr}s(in+U&|E^#WMfUaf9ZK?x&_F%3xn7yZn^#{o>fNc;t7aid(sL2Z3+Ht{9#Y{AJ(S|aqZrSlf8gFJ3j=P7!@QjD7miNF3P<tNjN~qf6*!jEtY!USd6!SduAh^l@LSezsFCR{Fg#`O=IG|713BFmW'
    '&V0`V^VkX7CY$??sIH_zc_9GA!i5)MCn>@&bm%7ZTu1dUa*l$kOwj$}%nq2S_Cs)Whdsp=?OBv}3zdYAj2w_jT}cuvmus~SJYj9_GP9~1Cqs89LajV5G%igij;Bvcevq)FpyCBZ<B?c)M*pAWyFwRpmQolfQc(fjYrn^x^E)OqG?TA(mt+w3$Heatbu!2oDA!XxUqR4_lOc~mTL_)@bq7q|8vz%T%8MWUT`9x{BF%xQOQ)nQ'
    'w<Z{o^7(RqFbMR9OQ>S)NI@V`7f<znYBa=3H)w}sbgrCAEn^q=_|b7^_zi6NwRf=X{K_rLj?N6H@_(A+RaI7ve0c&P8i6?R*sX}P^nE;z2r`_J5DBxV_*Entzfpptt!14r(i#P%jH=0<LHsyu>cW&v1Q_IPxP%J+BE}O~C@%E(-@JTfU;C<60G&&`52LZNq9lKAioqZF3wE3F6Ora>3v?WsyRq@in1%)Os;d{naSv0=ryp)|'
    'B&^#rifEcrgdA3dR-#s5g*Q|$z}|Z*`skRClDIt!9^A-hmjEOWCGfeND{r2WiQ97$Iv^M6dYYa2D}xN+>hCkg!ZosQlo!MavvimIG*u(>4o;)o@LA(c3o*2k(|>2ybsvWErst+keP!t=EN;C<w3E-2yKM+Lgn5*^3%q339eg$!kF=5>i0dogBO;s`gh0~o=5TBi3~p9;=6W;j$``w<+|<X0AZ1Wq=bfmb+5Db{$36Z-iGjN<'
    ';;cFeTV+#5EH54J#>)AM66T~AmWruOxk$>@pzF|%9XfMi^M$6MRwOs5H*}dyYO+IQnTvxM@(aiY#%VD~G1S0Y3&+kloS=W~Fw?(i+mmVH$jLw48MIh6H1AJnTpZj*+goE(k*fj&PkKLDgxglwpf=<z9mPnx%eT`K)G)J%$=Tt;mHKUbJ|^Oc)8Te%-14`rt7Ld(dJT^#hHB=!qMtBTYfm-TEQavY8MNwrh;gN$<9YIx5p_&;'
    'ZE77PRp}#tf3BXnZXpg_*IRFRU#Zuxnleb6#80i;9uSPwS`nW4g<6OAtR#IUi*hhaH|Nv^$|XTFqT(#YsZX<Y)UqgxK;m=ifnSh=@jzH<y}i=<3%*@U7I1hs;}po=h6tGGV3Ay8RpkMj#!*rRmWgTj!MB!t{4<W)3Q?WZ&LlN!->P5$SjGy7oi<JmTnxKFqk=YpnnGan(Hdap4f)9&O|tFUYA0Rt`h>5~A%f7&Qw%ngCa@`P'
    'Ucv12dQrd%hTMYYA47@<EEUM;^U;=+vd=&tvo5O?2(Qp@6vZ5(hVfn(MiYMn==AkaG?tc`y8PkI<vhL>2WVHobt%!}h=4S*Z_!V>>K-p~6Li>6-RoRhkJjJ(xK%8Q(t|!hE{r)>x6Mgy<xH6O2GD9DZHuQzNl|+Tig7yxGm>K<*c*6<EIP0-ztc2k_Bd-Z-Dt>q@ok+sZ&wFxI4%FJRHw^q)X$+2EGjiD$A@+}r$Mwkt&Fi@'
    'shvAM3QVB;K>bg5P5Tx{Ux?QIl8KYSm=hriA8#9!IckVpI4O3r2;EMx=aVE^2U5v!p}vQN*$=Cmh?*xtO;A4rqE~T?AquH}-Uz0e-z5h)g!PB2KWZd*8Hnh%_C(MSj9a2_<F^-Kn{PE-djkH;CZKwUWp+cuD643DA!-6*6G#x4JsmNa0?H;_jfrYc+>oei3-1CcFWHyt$F-|w1~&IJvddMV?Gi7}v8h)-syoWpEFZLE@IJSK'
    '-yx+Pt0EpXfKOKUL*v4xS|IY7ni9^E|9B$!%LP+%UGIebr)^odrP_Q6rOpo1O-`|%s0KLpm&($BfMNW#L_^Z7`$*cQ8@hyDvdn&aOqIV{+G0ol%j_hf^dLLIh}ar2<xe;bL9QXl`PQf?v+2E+{d3@u1ylJt(;CY{x5vyH22*xNJ@1=u-y8ClTb*we;RK!VSBItY=0*s}EOLq(30Jo(smBAh;W?+6zeZkJiaN41hn0z4W^2{x'
    '?bl-AfW6;mng4v~mGeqQ{;VJ&nX%+;7Pf9uVCQ1@?%vOX=4+iIv}8I?v=FwXV?|I-=oLm+0wT{|XI(6e4MG;>_wTJ`7*hD^Lnp@(r?-Q5<5TLN=}`_GOb?6B@67RuStOpsQ2i%1J4&ushdD>ds{YfJY3?u?+o$5L_r_;T|L|^9G6v{Pu1a~wT;pL)%nh8V#f4-@6885EkGI-F3&SCVhy-px-MN@Hk@f9#u-+XLLv|)!#cXZL'
    'TR>pi?mFO+11k_yBFevQO=#pXq$*K|?Ne<n9nMDiipeimJNf3Vh^CivJDS=c?6{)7s>FSRGm<=7dCy()oG+6!-PxP9+^!(hL)iPZkQe-AkByKnn!XR%BPi|`=xkLjma=%<%-t(QM<)8xyvTxn(T#{hE)YU0NQrsR#-2*0o_M|E_`x|t!nh&zN(Ss^HXdVtDAx`i$}1l?Ue^_YOr<W(5=0A&Yk}sR^dee%bohiO(mo?0_G?%|'
    'v28MPM&0f_mL(pC<|7g=-@Rw;fA#w89QMF``vflF_sLWC8EHm$^mt&C<|VTc^;W=$k%D(@{tw$TL?S;xs;*SnN~KP46F1X$bGyby0qDr_mqzhlRV|tJWw_H^G2|bLb&HLAB~gnxV=Ney+vk)FyyGyI<Lu$V5Bf1GN4NgnV=dlWaC;Vz8xR5L2Cy6%K%p+lAZP;=*(fnR^Q*QH-=42-N-yT?$qR6st1^Qy-F(mDR=01%if0Ll'
    '77IQXJ!@!%jE1mU-}3d(;5cH;YW*&KK-5f_^c#P$SI|cnUY;G=z&|F}&ANRp^hL&riT+mq1h|zm!-1@t2v=DLQs2sVk5e4A?c<|_Ajd**rmm&{Sr%Xmk;9s`5gkh$wI88*&n^L`N%f_2UX)v{LKJ1Gf+Mk$8R|orh)24j+B|OLr0gaIaf$fGhLa%c3kSoIpX;3yZ>ZWYEB-bXS+ZRc=$ZTtuijJV2C>;-pxJ)>xuDFb%Z_ZV'
    'rxpV_xO4?ac0+Q<dH;X3%W4eK@JQc=p9F7>0>)qk<^#0lm};L81dD?GaQZfU)ES*@^Zb{xWtzl2fVtEmnXt>U#T9f+8BVIN4;&h@0r?LY&$^GJHGj)<X<1<+2oXXEyS}P3<COLfZQ?JI8oK3;Gy4}&mLf|1ww}g?NS}~3pV!0hEv(XNbA})CMTfHhFyfuE+(8Y;b??5f^8#ikMjXClycI$$jKmqUnqpMOPYBaifzM@gg5Kl*'
    'fY2pLI{1TdYG(9aH{g2gp$kV2Wtz$!8NsSZd6t!}Aokzgv85s^P-dAG%4OcGf))p4_e|T7^|;E;_Dh+~Gm?#PgtR-Zc^<o&;9N--(plrLoy6dsS|p#U<IrbnJZ@y@?>=dboe=;ZWsNNV{VzP(a%275EicmmdGrw#0t)lccRN<5F5{h9o8HkWl~^wz2gspWlaR*iTbcJo7O%^<<JA*zac`SP{fbKs|G!IYncRvao-LU0;=s|y'
    'v;ic~Os0p|#U?i`?7NZfsj|SxA4&3L#Zqe@P&v~sRLp<088Gro$bCqa=JHUt%T%P)AH#cH2YDxUeU8K=quoE3gnNf-J;-@T@h5^R(Mjc8QEWS?Sc5u?h|tVwQyWbyC0+V=!A98(HaK){&??q|U=;A(`7$cbt)%Z|jPKkAVo;RQmgW{taPdcdYZDlo##capzzg8*2-nTR3+Noxh1F>VBKGPEG8LgWPFbCLL{}JTXsr9b|NGmj'
    '02NBvA@I8uZdG-C(1l}iOA&rprdr7|JXdak_N=tr)arJ(PRXVLq%{2zGqev)o!}Cd>a4dg8|VGdsAFWA;0M0T9fD9)8-)S9O!$Ti1^7^JGP#1ho&i>bC4d_C%Nr2CtfO-Sdyxii(9zvX%?G@@y)Gmr4i`u&g%r!6ibJ6-u6m0tGzlq2uv-x~_^uE(Fylt3TA3$Zklh+yfUplD|D&L>sF94_b=#s7`;%f;FyVu`VZhB;2^>_S'
    '_@2snfRd?o0%~7nupEiX#4zUooFGm3?!LWVt>WiA2E@rMm&#{<jo{HL>~z&^?{OZB>oHZ?wc;8s+fk5bu8d)m;@bh~bK>A6(0Xz@WYv`tY6Hlx=yje?34IdrfDT@uY!!^y&)6XYrlfZtPhwPJ6-mxd&ExEk$qFy?lTCzUKBPk{sHuc9_@3gG)wWor&jmbt>kJ-BgXIuly<Wph1B`DWWP_+SJM{pUeU$^Y*0Nk*c-U6fchoL!'
    'JhpWyKF7wg3+$+_@vu&dDXUif%3?BU-y2vZ7w~t%l4ZW?g@cLsY62xhe~9lgo<eHg9}39l*|$&L!jWRZudp3MwN(q^y{@QS$KYov9J-)AgK&vpYH%__f*l_1oQ4FA$uFK1Pzys1PS{)ADuCNSQ&0d&W`)4q<GbBwCHFmGQ4H5gtKEWo7aa^&??&PLs-D&6OMoy?dK#Y2VA6fQSeZ_HK3wGWp>mg$ysm&i>Dm-7xbyKtkoBsi'
    'HF-b&Uz0*@eez}?GKDaXyS$QX^McwNF|Zo(rHpow&Ks#>6gIulMfI*_*txrvLDemuDR-%lNU4R8b6u#I*Hq*VEMDP}05jQOaMLSf6Vi>k)&a+VW-$}b=GkI~^i2R(?+)olO&<C34|7yJQ5ZVIvnnO_tS#=?l><SHLVffc6t#7h<ycLgOl#u`OT}}DPPnK7$RYFm_p@UI?NO*|;~^EB?r7fTsfACksy=qwyh#MXoJ`6SO6~OZ'
    ';3+J<%3P!aOGwE+{!?X;1QcnxXNKw@D0wn+3SX!^5^PNbLl#mXAZC-PV^Zw)iDG48U@3>PnnyFPb!}}{g)0+3xnNH9q=?)CO}?-)s=X28JSxq|fhQK0z4Z{K1D=jQJ80c#I}}T~8K4f1f<}gkZ1Tf~plG)$YrC+br$<JTgJYVUpI=69c(Wgp$f5(<PM$8L=fs@CT}YP6$n-a{ql$M$0%S3*1_8SLkj*&3TO#?LQrj5y>v|h;'
    'zo@R;j;17U?dU3~nSdYvo<j)2Q34heEN-F@AvA({qq*vRu!)}Qxv`5x7J1|gmF@)TN-RXwq_pr-j{Euq_jnF{%WV`V3%nDN3UzVtA`_@_(O?(bF@sa}B3&Bs&C_T~F_}yZc5f8rYHCINTBIWaF%6C9Hy9iU9PlH9jwKbthl;6!qEqg;;QjkxR^c#g!^`_LL&5O!x7(Jl7A=HcZS?eVuP=DPc>k}JDK02kIU_Yn7Fp390t&=h'
    'QtmW9t{ZYV%T6__R%V^K?*%mS|DnWgCZY#7y$2g4Bl9);H3Yr8DRrFSLH2$Q<Pqw0w-^{U+9bipi%2t}egj<;VmkRpKcWFMucMtgPqL&a_H9^9d`;%e=mO`5<L7Q{Gk=Xd&n+{52(kXx2ZJr~D|1KTqNe7#=s(hP_uU1MZyx6#(>}V$Yz5K}u9i{vny=JaF!xlA$){p+{fkdYDMlqfOQT5!3N{{mPjR87Ft!+Wt?Ah!@faqb'
    'Aq0o~oh|UrJ8EuM)gkkUACOi^>aCr{_QZkg&iThHJp`~WRf^t}w{J=t88E5tplvn~9Ie>?u!mOy+wJKx9X^pZuO2%E7a}MN`bs@=!iu3w>H-yBnV|*=fy@LyXpdKwqbs(y*TtS8lW0+7rnRoA{FsRGWyUkJA8x9)4wDYpql>|z;yqFL*__k7veg67s3`6)FzM2q8o*C0AF_-I{cBBKK3#6l48Q?0f8(R>0RO7+G49%D8h>!Z'
    '76Zr^*o91xY<X32&${MJU{6c5RO@MlLn@8DPKNqGP-B6ufs;8}y`g@oJ`12!dm-)Y&@SQHwVENQmg<y3uMw}-hCU+ri<R4;KLU7TO=W&7pY@4;d<>NK3tp}+2QR^k+n^$EC{UPaUT`EQJnWyoQ|Qz5V9(4yX}kP8?ds4xSmi}+U8sFNQvHoqGPfINr$FT(^vRv<lt!Qe#*BiI7?f7Gc4|Q0VL$7hmohjzadJjO5tsC8<Uvis'
    'gLsB$s7L7&+rvVKGEcaA(_lB~K-^aGXvW80CV}HAB&FFz8Hlc_0#rc=RG3_?s97UYoL`v#{sDCUVcqy#6e*A<hOQ61hd9ccCSjgJDRE)j%gR)P4w$G*5gQ5DGrq07o~i~Miwtm+${F{Ljmr*z*|CQtud7C_y|TyURc4#cP{Y_($;(m4F7__Ru9h7xl=m@B_5ClbT9^=UeSseC3VqXlyqRuDf*(b5TYZOUtMIu$PX}(fAbebG'
    'V}(|q=cRy#)38wlIwHdyw_lsl@=^Ntv2Gb#VIVO|7W^ih+1`BxNq)~=4nTNA<j57>K2)JzC^kP*=0CraPur>rjyqd^kpJ{~VONh{t1ImLznyegCtNJbzjmbgRusRkk>;So3WC~p-Cyfne=(9gW8djJ4toF=cSEaPvSsTEsHiSmWCMo{{G2E}srnxwO7oF`DTIHhwIY$7s>otYhtdmb*g@vxZ%~^_yPeuMe$iGWOiL?NJSBW?'
    'tbqXQ4Rl|UQE7gmRfVS!fv%8qYwRvzq)CS5-x6P>V)_T|YRP9_FsH0Vpp!){Ts*2C8s$cm%9Kp=Yqdd87%@Igwl69lMIk!4?xFkQ%}@92E<d}(IJRNdv-h!pL9K2(Q1wU)-V?>z&lrH11uGv`w3#9~!lS`;hs!m|kS7n>E<dgleZlD>Pfs&Z0+fyGIj_9Ru!Qqx1X82|C2?%;YwqOapc*G(H1TnSt1DHL2mW-1n1FM$W)w4n'
    'kc#9E&5UH<4K%#33oRAs2%CNgB1VUcneV<_Dw1Z^c`vq`?D0js{AgP#L|Teyt>U%?Wcb}y2>c$Ot9-euO>W|foA7%4e56KUa5;n8?AK8Cn6V|FOET@r)dmXrbZgDY5e(A^DeLnf;U@&CT}6fn#M`g7psf)TEp_z|_HY;6c8q)3*R);Ped_=Bgwk1K3USNl6dBwl#f~4i&P=viMbZNMH@>_){@3l$vSN|cX(nXpt<Y8+JQ_hC'
    'nyiKuB$Jy&fLF~p`e6wG2=rd042O7&tj%6+c3C}exMopx@aq+@(MYWDB0rd^_!94U5XC9X-Rm+>!q&B(lB666{}}<u^8GR;$GXJjnjGJXnt-qLV)A4>RrS-yZ7Y@-eBr0HngiBlO}kO|^1(2oB6m`y&M*vkJdi`e6c0}Uf}vWoscFB#GmFus;iB-k{q-^R(w(hxkyFMJZ*?EiMwjYq@;MCtruG{TUKWs>x@7a%^e^q)1ko~2'
    'q)A<uK-HK*ak)bRf3uXFM2>i%NA;^vsBD_O#mQ*0E7NjLb8fV(l~|?jc>w1nG5L^amNx2NG_NGKRB3j<xXPsR{UxK?d(p*p%4V@Gm79}%l49QGj_w_^E#3_PX}mAAmo~kB=*eWsj_w%NKoL&pw6dRLWhNk6k$*dIRbYfiU{SjnP&GsoO{`p_5^KA|6QjbXBiHb4HMK-90!JDUp9F3mj4Hw5?M|<QRXOw7(=y%~9=MQgE0c0i'
    'mi~)RLcGiqE+?E<{W&FBd|68sGJ}V+T5bw@<jM-Xk-7V=k4Kkhnf&72h>S<}j%(1)sG5j6K2!=+rwpz`6_NCtuiD&Gr8UL@Kb8XBErH^hZG^I0aGql<j|17A|L%TBw58Z{_F;av<sA5W8c@HSY_BqsD&SDd;&eE+H<t+v))cBo-|;~8rDzOy>C*|w1t;fhvCY?19hqGtf196C8_wRm&p@7Ez0w=+IWm|xCb}7S^+9~@Jgd?1'
    'LrhTTv}w}i@C|hlzvO|e{gt0aC{cl5D9I%`5yy&agSDk~xI<YH|7m5lHq+IEor6jSQz9IMfp%*}=ppX)j9`_mD8p^#k=h3T@7+UX3AG(J*VjvW=AwfJIZ9kJMW!ra>e!%%fp3HHK2ZkP&Pz_%RC9js$tvcrhC(CQqK_daFT^6#C5(a)vh5_9s}nMzFQpISOh|3|YRFTMc0+A!(Sr{RX!HuvEPx$i%#by2#~y^}TlP>^Cbq3k'
    'i54|8cmX@6^pN+^?puXwV`Uw3%ITkcuW4)?%D6{%Z(k2A&Mpwmg~AY;-%Q36dc^8Vug!o0lnU)2HkS$YeV*FkCW!aI1$v>5Eje6Ad4$(1<C{r0hd#2WcAhydsIdT|@JDw}S|55rJ?P<6uniUY6%Q_BF(%fj*8gDC&#jJoSbOex&rK{><ZGpUiv0)<O6HEcBL)4^Y@b)TclwhJ;>~@~;t|6$spV63hRP^c8>!6ODZ{$d0Gy7F'
    'K$&2l!OzB1XLs2k&1>C-z?qUbpn<vNX|60fR=Dw~7RYH@ONE#3$Sl|$K*Ic!tmPhh4|a-*t<2o}mn~QMPEVx&R6<$ic`Hb^sJ0!Ofu!j1Y}g)+Ua6?Cp-96)8-<8gE^ypcYjPJ07j1+kYtWYrV7~EE0;o8z;v?N4Yo-x1@aI<m>l^I9t!bC9nr*j(El8L{-#CC{%}EhEnI~LZYrE?|?i9aBv-Er$7xETOrStPA;~EtY6yL9x'
    '?^>oT%GXIay79YErDOv1vF1PvI7p-AJzT7lUcXo9W3hqng{tsz6F-(?P6bT%Ct2^lA4HPIUuz~ECOpxAy=W1TY_dY<ssQWg&y$Ix<rk-NP(P3PW{XeZW9ngmpkt?>Z9H=n<B@jbWd0vB0E!GcI?b?<-E=JLZ)=OO11>*3wpM}Nwdz>Au^}uj2N>oB;Ph?Eufz-Sy-UJ3L6>x4^_MYhrO{_Y4hFGMN|PLy!G*pHPsdyNE;jXg'
    'w?nyV0SdJvUEj%MT`feb=oU~6evfUaS|0hu;_Z*lA+QuV7idd?LD}Ov=z(F)Yw0fDCbSXMZtwGrY%+YSAfM1P8+wB~OVlE&ef+0opfYL_>GoOeZpdbi%hP<eMPQ+e_(Q9S{y9r$nqNox23_EH?V{jmyLg*`QDzmJT9J3P8HvS9cs^=_5M8TnUyQ#Vm1`#Q77Ig#pAeu|gJcr4`BJyg@7^eZdPINcbsTZvrXUV{rKV!qcd?+Q'
    'G_M|kHd;LP>?;qPLffH~z-Ryo!`dZ=G)*k-Q^{{F|KyV$<Jmz;cM+ez1_k&kV-a-jX1Fuf1<d8OI_!F!Au%>8&{04T3A<(FyE>QMm%eqf*buKca1w2NzR(#}d?PAQXxe=%wzLo2HUSA^@ApvfH9fG*G~%ggkX7*ZI;-@mac9s?m&4ML(8+F2xbVd_QnLb9$czpS?B)tWunQZ3=FC-5NxO^7^#knb7q4@X$mR9P|DDibjU^n6'
    'v0YW+Hu7sfporq~X&vydxJK4RhWU_u;6>GuDByZ^wlflhhZ`%Ro(LLsW&WEf7h3~qts}%%V5rMD7)@*^4lMkmAJ2O2<l71}Y0n2fLZMSL2Kex<r+A^GgcR0k;Bq3iBCK=BA9BwCABR#GXeJ!mDo4U-iFXA~f=)I&<mI|+>gpA30)I9WJX(qSu=LVJ9GEGf1uvw<^0zFZwUcu{oTZNqi-A?5ZL_eK7ZEmw9#-VcN{dG`mH*+E'
    'r>l9hXe7u$Vl;r_w$fh(GUqPr;zEccS-;%32iq}BDzrf^J%}THMu^ds0W6(L;iTy~pQJ!y!g1Ads#%uF#(}bX$B~10qLu^R4QPj>f!o&itqi2Py;C*t&y%HvNysx2m{wJV-@^<4T(cPB-<M8ekpsbdhfFeU2-N1YsSbKs5HeuFdi8Qp64{ex)E??UfyT%Qm4t^2j6~rzY!a&EZNIqT3Vpx~ilj>R*<`Y@KNDd116fjQ*lHhv'
    'uk+rer}#AX=o9jjas|pVci+Y(^cqAW;ckq;W5PzGpk(A}Hkhb_=X=fSWv|58V``O}ztLZveK)rFeN0ra$BvJiczgA!Q>)s)`*A3wSPt-C*}v17`~Rc=DE)U_a4EFSMAzX#KSS9Rpf65R$8^oRs`ow^Cj6=im5Hc@;OF2h;kvWj#>&LM#4FLk`iXND6!isT^4>#g-|s2zjIxq(C(38j^wUDMbAPsggo&rSWC8K=hrPyJ3ovI2'
    '_90Jq@YGUvQqKg}bWtL(-I9};^8)30iy1mpG5g$4<@zUDG|?KCJulqVd;n{+c3W2%8u%FUFpd-BPoI0;LE0WL@mumU(It&sc%k9(BHf05+=1O}O_QpNwGgc3=4$=Z6KEG({+m?h;t7_DOeMd=j>lZh7RA-JU_HP<%Di;LA^M}LGEQ?s4Z7u{`{FLRd^txO9!M@&)Oog(+v7<;$DaWXEIp$HWH5iPsq+_5JwA_l`Y@d*2~+7R'
    'r8qT=XW4y^=cLjy>IK_Axl}3E*N<<y><Y-El!!zD<+;+Untw7`DDOFQ&?&WkD&{Xd8IUU!H^(-U>Y{=|Gs9XPp*Zw9VCIEUn<Eg7OpbBA$P4I%c+?u~e~4bL3hoj-+?6*CH8=K8Rk0?v63K<D>Ebb{lwu?L1JQw8LKfWYeGlgH*w3@%6^GA*?>wDwjcgpJJMVYEh?=@q>~F*WX6H09K1g6p5!bh`VV+@oEMvx#bM?beLKteR'
    'WoeD2UonObX(9q+gnj@{h(5KXHkn2eR0!kSD1_E$V`3ASiJkr{3$jAZ@mx5lI^as-4Xsx_vtHkej^)_HH(&EkoC`va62)IjJV>jp5|(?-*Hav6$gA(&cs)?A-FF9@xJDTeqDl1J(Ml86by#7$93-CYJlEgA)VHDbqT3?)CQmh0KUg))sXg~Lh4+3_|DS;68F2s|Lt9z)kVKY*1WgH&a?Rc~??VNoRk@WY1egB=KFFlHbf~V`'
    'S)jnOM8Js(b5=OXBvTbv@sx%(-+*7xA2xmCU070g*boXq|DRDoe@EeqNZ<0yNj&2N^j7KpU6z(0nmKafL5?*ONM2d|gS0!tJ)%t3^&OkNr70Hr|D~)5oqjc@FWK#wOpkJ03MA{-i?$6VXw*t%79MO;__~F)ER#Mq=m2xX+BrAIy|>#^WEQCvb^rhXN1<x$6v(+I00G0R0*H(cC0=+DvBYQl0ssI200dcD'
)
