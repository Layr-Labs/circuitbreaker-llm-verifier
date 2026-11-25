"""Active -3529/10000 first-pass, gain join and starting-majority proof."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial, isqrt, prod

THRESHOLD = F(-3529,10000)

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
    ck(U5<=F(73670242301,2**39) and gain>0 and bonus>0 and 0<cost<=F(123132461967,10**12),'active energy bound')
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
    'd08X~Qti(x$}b2MrG>`pS)2311BnXoj129AGvp1fOed$*QxRgkvd;c?)SFfa65cAFdCD$}~jW*;`ks1Jxu9wV<3#niPPN$-GaI;%^*~sOICoT7Xp9GXiy`A5tUp1&Z8?{8GO<{!^rIuR&6si0?}Hl5-*iJh??ut#ewFVY`&97`@Q4?NhF2mf-7uQ3qGm<N8rAHy^Fkqrd;JZ#5;3j7JenNO+*4QdJGOcr4%9$sn&?&!<^~KB'
    'a1x0_bqlv)BW899@-7-fLZ)wpK){KL8QcKlq$fKP7ko<J5#=Y-&86j8wCOk-|9%yQt;oIYE6n6B4(Q3B$H7xD}_UrqhefEeMe)%t_|wfNc#+Pj{e@+{@Mal`^>+BRr<S7Q!cmn=F54onN3(N|7(bk}JUA_qtAzX>|qt{64eb%XiSren{aSu`ha+9MsocI1d++{J~JjyyU`Y#FmwLLrzyn$FB6gdt_S>H$fYb)@N>7kMZ)@+`'
    'R|QLNs8~n8<be-Hh%{<LiE>KXXFT7BB7$IE1gqCZm1wTb0V6cW&$@0rv(NooA%63~8sV{{18m90ao&Vww=C^}@6t^N0Fs^6;beubc&C1O{kh-xEsNF;SH|Qh0})oL977!j;}E%{j@mb1@WX`giCiXc%$zg$4Q<?4<MG#r&^E8&o+(cgHVm&|7X^>~@;GKfri)Ctp0fea#O@jQ(>?iv)a&zJ<RhdLsMP9bw!J;4vW|MM7jb{c'
    'bt|^aKmIaP$CSUrL=D=q<7iEPYfWw+Q0RK<QFeAOH}kNfbf&OE3JR$CF>*J9A2jae<lNBf?YWWePs2sik`n_CQ#hEGza9-JR?F05+O-9F-4h5m1%{{EM2*3o!i(+vVXBGOP3;&}5k<V99D>Dl0OH=y8<*FC=DZ^I+fOj^e@JN_{YUNe9_${y@#V-B`UnT>5b}hUF*osH4Um(`7E#8=Tng6Z$7uh=Gb-!e5_)@OZd?vEx)yVY'
    'j!bQuu3DSbDxoW*0GG+%nyq=X(ra^fHU0lS+@MHq+$}X9%W8TdQ^W*Yv$I<!hF_7KG4^(eNBNivi0oW8}#dmGh#Fo3A-)(CG^lz;}KIsQJt@faIvUlyqPil8gJ%I(8WHWj=JlUyzlf^o^w{U!;*@CM&w2I6uk5VY~=$7>CH-No@X9&qp$g##!mqK)-)sbHNuE(*hj-v@T+~FXOJj#R9at0uH3F<pl*NA-prB6i2A9)s@mO;J'
    'Lmo!xZqsu3_Pr_W0@+Vz3XL(rdvy*%HO`x{q*6JZC$J!wv$L^^>jEcAgVk0D!s`8F%21gj4VV<zAl*M<;yF($0~u8zXqb*F=)Wrd4?l>wGL_MZXHZOwwja?<4$DhG?{-6w;^$VWZ-+4|uCRwNq6-?KB1Jlgd}=esVwyJfG+#4n!jGZKNV^@8V=nD;I27mnr_EiEAoYQGAaGKb>tgg94)iBYwQ1n{bFr40QYTky!fnd+F?!Vg'
    'l0UK~$cHi+{ChR_)_UDh18P)BJi0p!i1n{>!HgVdFmHh!vUuy=;4%bWuX*Nw$27botDuzBZYqK{Zgehab3y$gw1=>+yKHyYn1_f;Qp?+TMtVNrFyXhVP6?EW8G3o&Gg?Tm2hBE2PDFjtm&MU}l;tqBk!YnaN<cKS*4Hcp!+P)9z4%Tk6cO|AH@tSEb!LbrQ*jue+uVjHPSBvZ9k&6v*`kH8aoX3wSt5j8RrzY@wN{C~GKz(-'
    'nf}-)bVp1yUt~C-FoNObEkX$0``u#}p9&$@J~32K+c?E@2rjhlyvF&6-kp={aG(O~xi4?%@4Xsl~o3c$YW8fHBb?)lOFkNu*O3TKGfLw+z*e>_9+Jh{`Pz2vrn!D5hI!aoH5@9_;6yDDz&ztDK*@s&*|jb}}vBV<)$^pEpH`17VOl66hVj@%&7P6?%V<gjF?2#-?qUYVJ!>m|B79iF7?Z3n4DOwu2bUSOiI<`(klmyy$_Bzc'
    '<%P`eNb*_0j3?5y$x)I_13`vRHXo3)!`9aRN8@pR$Ultw!KeUgBt_^ZrWY>%zH)*HyR8lE+0*5z121>0gh=2oaK$eOiA5wKc66*z8B`HBp${j$a4NQ@N*PtN46LD*T>1!&pl*RtJo!pZ!>Sm68z%Z8JrB)qubCZ%7U+QFYHS9M#8MVa-gzL65j{1xu|epx56U694Ri+@6&irc$^ZV`pm%X_pe|?VWv<^+fn-Y>DDPHD$rh5@'
    'BrQY-m-EA+iz0ghaP?T`k>ve9OxJi0w-ddA2Blpdzb>iS0Q;@ILUG*Yp5g4;7k6$KEdU7j`|Vv*HYmp4Yo;ZD{27n?UNyQS<#i2DSARx`RuvbACku3N-7y;w@&qX>}2X_fi0I>`}#`Pq0wva2BK{s?#g8-5W#GlVxb=bken;o|~|{4OX4lD~bX2DaHrlKu1P5$uD|P?Td_^nL(4F5_(gn$NWVbx=ihY6SbmJcOnPwNI9S>5I'
    '$b+2AeRzo>|P_@BT+x81id8%=U=gv_8}3-fk7;)(COT`<92gH>Y0}s}E&K5%1ONW2eAt<<3hl+o->!?jXzcFsIGuMF1j3OyKi#ges_eQih~Hf~moAA<RmP4O$;sM5$A0=~xqGtZd&0l{mTpU{^f&gs;{HBD<?;cLwDjA=@1U^#YM$0)Wbz#LAR7mad}5BwsclO!F|CvkV(Wq)-){EPX}FPYQdHl_^uYMh1(dpqxsd|NEYZ0?'
    'k|h59S@%!3FO>D}Q}l#b@-MyQ=;;Whk<RVi((569B7agBzLFvs-foF5&n$-qtX+&f)?1;sg|`%2j6RrET=+?N7{&^1URpHNQY5it<-fh0=`c-a7rk_9)*U!hQa))84x@p#qhgBz<#YnanuOq#@gK4ED5F2nw%=V^X)B$6(r_TJs-H2}ALEibO%aE>zn^1@n%rT&OnXGVNDj_|5yYLN)9QfGa2BDhHj?P}K0{AY(;^m(_f<(C'
    '(u$4iR&nFzJ-X{)HY$n6w<&n6%&lrS%x^KAQ=1Qkuv{WTJ!pT(@Bo7~ebS1jzyKMBNR?|DLYXJXyzaGn{y{$+E#N-BITyG<tRc&~&1e(T^>c_5U#KK{r=15lA072wOr+rvIZhV|4s&q52<ADl$TI;J!5QOC~jhuY)eZqx03`M@1$$-;tfj+o$sJD7^VnAjeX>kRmayF%AjA%jadAbS_l(!PJ$81@uFnSc3t{oRw2lr6|;fc!'
    'Kj(!{kOI_Jz~+@0K>}QD&?W%>B$_T{d2h9#y>{o+r8Oulh8b3=jjF&dH83SPV$d9oukA=T=mhbEA+L<#-hB+?38(=hq@`YQ4a6Y;E-5vaL%=sFqA3Qo1C4zJ2fr-zQE&?X!bXHMhRGeE7W$Ks4Y5O0l*;Fry%uJo{-b1_P65@Jvj8oGgm+nTp{3N?#{n!p>e7D>Id)533o-sx4lMqT$|55cTt;LSwsOxBFb!d(vs9B&sBYae'
    'IV%Ut3Rg=r&>e;~~oNSxeiD@@2%RFk+w!quLF9@M^4`EL+hN8d$Ii1HIb(_Gx6uKGAx9X52v!j?O%^KZSI&Q|;w|542W&=mqN;rG7)InEa$5Wx{X?W3#?WU{Ex&Q@dMU0k8~bQ2DYC@}Te<1n6G&67(`Yf_oGoVp-AV`MOCKF`Q}Qg9bN@&9S@Qqk}MjhO`PhP(s!qj@r7>#ROdYF_4VUPtgemgNc6!lYm5W{ORG<RLphh5-'
    'D;mQ)hcF*jAG9@+0>`wUXR;eB9=e`Ordq>+8%n^koz<#ao`^bB+uDMv@j0&+mSb2FL3X`(YQn@yHM$)~-()#rIJ}E7wK<UTYm`6>p2*sya4k4^@*r%HdVgImnWY-WC77`!qnGt?RQG44Tk$60#coBIa%tQ`+)@rpFa9?RaV%G@`SwQeHk_p^daQ7PUP%h<}cE?ZJrbXZ2qMcq^IVq&a3o&p{#^i7+3^H;<v~c+tj&`yY*fNW'
    '*(G0dpD$N2%=b&?GY&E?&WuG=J6KyO9_5ZLUW_WcvlVsjq}b_4DhYL~!Nww*lVoco$nLaEXDWZ9$^>wmJ%MLQo9%#-0jYsGx*XkVQpMgE0d4%3>;`NF$jwsbC8G_4oD2p)A<Kesu?nyqaJTywgX;bekSOK!F!#+!s8$AB~?1hN?4Axt`&!4mm0vK?zL7JM9KyHn>)5YxIEJ`s;AeJvEZk>38?*+j(i}nK6l=w%I5x<aRULu}'
    '>e9a>F&wMu6~O4IOK6*zOijyj<#W*KB$bG6unWSkKI{mFSLT;5g-1@r<;EkfEUu|Kk({<q3s8K@canHsqG&)Npg0hn`9>z&L9o)=`h?Gl7psM@J040032t=%AZHSg!;kf_O{~Y%A>$kmE~;K}db`%9Fmy7j0X-qFEdS1ThISTpchVHhv*t<_JUkh?9)@4RB<QR>+83Snv=7@x`|aYlTkC6kM1nu6^aR*S#pGVwEJ;_(WW>Bz'
    'r|Tu)$Fhr=isrLO%wfXN1CNtX$ylG(|h$o4?c|A5yTdg<j6$%)sR?U>y4F|9BISC~as5M(q7HO;M``11w4o-7}4OO6@1H@Gt;;)pzgm){t4K*JpYZKPh7l08vBoyGA3OOcmIV1&=(Xq1|Hk3SI#+=B73n1k&qJ`G3d~QmuP$89q+TCRjvbsz%t-$8fm<TV5`+{nXR99R9q;Q_<Oh6@6d&95jmhs&M6B-bzZ-F~`)O6?YAl+u'
    'P(y<)%Y{_wRx~K2V=f(1qFy^me?iC^WkLVx01$yQL_4=siZ7n?d*8TvD0)Y7CN<bBKBpk|svucgUs3npmnh`i5gs_cLwSXq60!V$vuXRFP}~#43N4eZZghylZ9^pK)G(1uwL|<@0qxsWMr``n!B^VIfIu0fn-$QpPXxgAo&#6B(x?Xt^Oh|Nn?ps};wMsHsk8VgT^5J|NzV(NkOPojywjzCmfk5bLiz-`R8Qgc*UiN*3e~vz'
    'wa_RBg;e8n#>8n~J`}p3zLWTnq4fMtlD5D2<_BBRyC8RbNuo%*Xh!Zu{M-fO+{BXetUu<YiJUdH|QgG#$j$F|wmn(~D<Rl3~RfU7fy_>O&)@^KX#UarQQ?iO1^SuV0<K%*z3RyGyL?wMTzYX)ca^aUj$NUylj?4Z0y}t{kz>vEO<dqM<vSBGUPISY$2Fayhsm+`+eSs%0_?yY?>>-k}82W@|pS_HGiW0m#z1XV$r`0ysFJ@f'
    'bV+T7O(EIS~vWy>jA~it$-mhTq<Y`|$^)wZX@fB*$9fqf4i^o}6R&;&^9XJq}a@q^?2UVD~d5LaL{reQ#m5x?EJ~NBSDny$`7lT$vWG?{%Dn#bQ7*5Swel58;(1flJy#g5|kX0C6tOzvqvew#Gr;9q~Gkq--YnW#hKI7ea&}yGNaIdRTGaVQ4%Qz<IR{LyOFFBgI?!b|DoZXT4|G;L*Q&heqqxD39bh@njy5K7ywN9eJEIH!'
    'S2N*ZB;v=P`AElZHXF-7nI%EA`!2J%A4RnDx4@1}f7ym<&W<fTDkCGVoHs`H+SI(WZ~5y(_#bu1t)vBUr`j%tguj5S^ROZPSUZMj?UXn$^jVAkQf;zem|?#{_<8`0mI3qOEOT{m+ED<L?z<`WkrPvtgNo3QI}5z;qtoL7~wYTUR|%dhUIjjru=?r~5W!0A<&KcT$U3{jJEq#!6JF)+FX=JXFKsjg0pJ3o@AWn*W=QBC86uMM'
    'c(9PhV8~GV)tJBt2Nq5X<85)e)h16(ax7Aib@?X1dJ9iYhgr?}{BBUT^PeB2o}OOWEu;Do~jE+q>rO^nL4R)P=-WxXEP=1O3l$_b(f}Q%Za-|M#bYtZG3XsmjeQ|B99hWbY#|9thy1qy{JkO1!-vDygFNy4a#@EeLlb{%fx<nYXgeDP3rnE*=VC?_l|aIzbXBXTOo&;iBdS%qAZBeNn-nRe%&v{8((GN81rEcEYzn6^ESk`$'
    'Gk<qf^#ikkg-!iP&UhdO_xAtD41(iYhZ8EZwsjV_!k5er`?qX1{>dDd5-ClRB@*%x~BP`&DEwG-$Wkl!XM}tNaV2M*~zXxJo)vfqS3v$nSsLfRXxLi0D0h67kNRzt&Jn`dd~}$EO+9WhC4rxR7oX=UM#KyMl2kG(TyyV-tysR+QZ}Q_OQIt=;)xeh|C`TR;*F*-kJ}V)QT=Q5!&6ZTG<Q7`_kJ$Talu?QNJ9Eqe3Lk(eW92g'
    '{cn)f-(E$znL6NP&#uh2W!K$7{MHg_21dWy-Q=Z&l7wMbEq+ztpQ1Udo?v>xFp*mVqLlRIPbWCKl1?EGZ|;j6na$Ql*GFejYSF>cM|lLh#1x;soLXw6qq3Q2!_y(MmdlCcPgEz=DuG(GdiP0Y9mC&>2%g^zU>_0w(BPL6I~2N|se<GAVQ+j`+_ygkir>4YT-#+QDw-_Pd7-u+3Drb+~ZM8HtOiMIHoWuUnOXpZ{64Sl*+c1x'
    'T988nho#g<OHRb`Wjf{(w17f~E@|`-R&HAY9;KBAyfoC$a(`f&ysRB#-^~{T=6Q;5@W)3p1l)0%~FGYUG3zp%2W?-ZoG@|Iei#wZpe?Ua#2!BZVSj2XS);38B93!}jnN!V#qrhe8n;*b@*o?;aR25VXk-qQRX=7C6JNh(UeUPQBgSx2*#T9Rkq-yU`fJ;K0;X84upsk^5kR8z{jQ^JQ0Ob;}FPoy)JkYi|$4S{c=G(_^UWh!'
    'zk^*m|(X-Wnn;n!EIAogtcW@72tisYXL$wY&SopPkz%Bo)mSyAf>&!nkx;v|v&|Wb7SGv{Uu|1cEiM8>vK<<UgOxnqWkR*11kF_DpP}xQT%0V<MEjyK==NE~aISvYz_Bhbs6>vQr{rA+N;r8x}{kJ#a2R>GLMQOOGY{;ZLML{Wue_x8}1ZExPUlybpe~Qj69Z&nH;iVvzNCAG%QP=1dr-=ep~QzUpMQntZ&vA}4HIZJ;S3)='
    'Ma0w#Ua^v@KQs8Q>ODuE(WTOo^z|N8lqC=4#y&2g~$$G|T|?fiLGG?Lle&t#HLchU5g1cC&yKr0@F$acCgS^(f^*|864_38mVjrNw{zCNsRLKX^JYw>`V~6YYWa1aL7?bs;ppyT#{zXNPU=8r3F;2NUW3_+QEHFkl%F!W=*-N&(uqNuyVi0zWDk??7{Qv>b)XW`uxF4;)5ss~<m6-O`n3n)Nyw!3qx>zm3LeSu?~*TD^echL'
    'oag9d{_2qY7uVyASA7QUW6T6;e!5NeHUs&7DBuqgh*1*jP6o(ltv8y_Cjy3HKyWEtnPlXL2DZn9*7Us&KxMPMFx82wN!*9=3KQPW102X{IrD55Bj;mV6m6Fj&}I{0mQ@DB^15bM?knaM{mFr+$i}BE;VDv-l%aAqy)c7W`jis>1?viU0v=L>{qb-9%{B?^K=`pyxQ%I0?h$UBU>Erj*Mg)9xK+omMiF^U7Ub{8_`bzFzLrvf'
    'xY{vW39;X~|XolZaddAA?iiQyJpxbS48Pu;jODMWt{mbZKqpJD7wwyu$wMJQ{{C3lhpca74`-^lt+~yL(kTIQ@ErbT=e3@DJWNu{G}F-&rS~EiV*JLh>&m`9|w=BD-1lpVGMLV*Uf8fcDey#Sa7JF9{3dE)XvGBZ(}X`4C?1>K^5YGfE?i88rSOZt(u^AdNEPC=LPxdSNB_vuPU;0=%wus0UxEz#@hXMMwNRv=aAE+iNbonX'
    '_kFZ?mvDUeM!}hx|yiEYDpdY}OeL0rX@(?Jf-4mFcW^@0a<*fT&j8Q#J`~J?rzCN^)`8+O<Kh)j)2F>SUjWfR+14bxHS9j``;ml`z!)qzRkL{bAMQ4rgyQ=G=>P;8X+cI<x7c%BBZ3A4;{|K<#&w%{=jrQc+74(EOX62OB;cUq)1L+_$yOVGLL?T|Ao!aojx*FEu^r@gnxI<%f8{8IJ%aDPv!$`F%`=o7wGcxi`uc+TiS$#i'
    '(PxyF%mO7$GRhk&FPcgz}qp+nS}oU7Bu}9)M*D36=sZuU;wI^*Fuuc^!Q1aAhojfy999+JU-;sNB%<#p;h+x)FcU`M>8KfOXs1)u3Y?j%yuuV&PIToBX)aG0>_sdAH@TIKOEAk>5o9C$SorJckh3tYyu#+XE^jQTnDtS?Fu^RLTyl3ok{TZJwx*LR}J*bLN$=2PncFI&@Y;BO&!5MSld7m{<F1ZJ0+DFxI<OhOLK>=Zz~Wp?'
    '`dQyNU0kvL;WB1%uuvX_wRs=gNP<x_&mC(a5sd8F>OYY1iNxVmkactg!px{%AvOgP`DOjOR}CU1wAJpIdylE3a_qMuEysI0&xH?z4(h)XgvieIolzGI*X{h+!)4#Ph&{FhqiadJW)CQHNaqwX+S?or}}(sbv<N+XW^WhJ;(IcMCaIQtC1rxRW6&C{u%wk`aNyYlJ|rQ{4lrIKixbQ!MtBaPVSCRYL&f1YpQId(Rks=d1+ET&'
    '}Dqc*eev_}0CX-khq4{1t8==pO4}U_D-h#(UBj=oA0z)$SI*-8wJm^>2kT%yGI&scFzv*f$?d!826txuVosCE(RKCg5H(5AQ@b@zcaHA>}MXX&SoZlO^P9tMR3e(3M63(^14BaB~4d-$NxRGumMy1yO|y#30jz8swJg*g1MVe!My3mY`*YQi1oPiCxd)!s30xUEp+@<YLxE%w3ke=<TFfOd)Dl*Q78y+K|j%5i{=G1B4gs&K'
    'fV&&T;FPn-3G1gOyOMpS+JnU1J7r%9W1WPqtV}HuW)ColNb5W!1ShKqSLAGdVF`oaz$e`;I0a`(>mOcLR9szmQ;T2{t`zkV+Simc~ONmR7mZf~vmWSG(XOr7_Fo7A2do4*N-4?ubxG$w@l2q!TXIek}!DIZiRMzIxGio3Ihiyx!PxSl91OsfMBxRWhO1!LFW?a94iX;NpAIWPRZPL9m|{A<HJ9;it~NbA;ZQt&Z-3S58BqCd'
    'DIMQURw(6ykO1&C$)mRpO>q3D2O2i{M(XSA~Tha*>QT5U$&ZReQ}=CUs(R0GJrdTnsxC*V{&Rr@6%~_J=z}i<td&;PbgtQEt=-;n$tM3vtMk5_GwPA<4#5ppH2JgAj7SfVZZw6@V%TUD&{sR<5R!8SP6rnIp|JNaAF`D22UB@|`<P8fId;-x)A>sr2z34^3ljqRnP$h&5c9h}4^#vVob4ACrT|zj!vlW?qUg;{-^e!G#S2t;'
    '-}qK8pnk`-AK)vB7s0zhgn*uNn><%(lPz^#zncWhqGY22!Pay8RnJNaLv=d4>>(4fx<tl#{PRq`t>q?Rm9((fiO3gr2S-(5@6+wnEn$n>hW?#%>vQ(4Ru46f6cWS+=#oAFN*TcI}YXt9;JHzz^W7Hlo)m{|J0NP;6WMM72RW%YqLDL)%*cjYB~{g#L~#Inf=2f$vk4A6%O5)xqvS_E28;25h{o=+T(alY|vD=S9YmK&nPlJ)'
    't$QHt?XJl7~7%`rbz6`y~3kXYsn*HB!9pu7;&}!U~yZ8j+@>n>K!FYAbN50w8#Dy^UjS@{NF4H2AYEG<Ybz(~#>1M4MYIe6s-WMc%s6y6<%f^pVqAkvdHu;q;^^YyeGMXIUF8Ck~LNr}Vy^ueRvg2zX7idBjbNG%OTMs!qF=7rY_{60p(cIl|r(Q!+voXD(QmymX$ESV~9K{~d&_YSH0sNmO}_AjFjw;yW<}AYQP1PTc$zcC'
    'p?uS%aJ_Dg5`yJOsPB$0eCZMn(#xw}k~YYL%uGf$Ak$d4O<jLYM6da?;=ozS7yZyXfO%9jBxS^mmT~at<wmqy~heiiu>TTG-8707S^n5)mtk>$*?p4%OazltqO&5lB_e=SZ$U+$kR1iY~+F7B!WK6&V)A6^*bi{_zEWy8F?;{wAPeAV8~$$X;uK+$ER=aQrc9J^orWT;&UrK|p2=q}Fy06W5vA4S<NvDV~mAj2yOt>SWcw@&'
    'iO3CL-it<w)Z54C>kUhoz;;-$-S9GGWJ!r~gL*?@_Hu@u=4XMzM}8G~LE+X*SYg30AwjzVBNcYYQ<zPJvm2Y)N2nen_9IwKwRS0v~GGR%pSA?=t2<<qVE?=M|d8ETEXi_zl>il8{NND?5g5)68u>s%5ZPOOFiG=<GTbP2yRC3Gt0O^DpGi>31C{{}`;;v`CXrZ?ac8OD@P-DPh8f6GfSv+!x5{qErp4uM1@)u<p3C-v($J9>'
    'pIQ!c&%ui(WNN_itXCDEsK(^^a9qqrz4?Jx_S1^B;;FH|YP@no$0<-TJ*L=knD!d9qOE6oFU<4Ma@=DyGs%res1IZTakbpr-1_6hZYYEKqmw?o}xOvyX3Xu42MylYFOjQo>+enyrAa*9#+RkUAKkk7%&Qm6zR`DKVPppSAQsD}5c{&qmJm+B=5BRWf$~0}HKO&RsfRY<7jzKeWt2GO(ELM&3F6i2R)dzYW!!LHS{96i})1%)'
    'U#{u#WZ(0xSW_b7d6Wo|734Y%529Ob$D)PaYm);>8CQ)__Z>qfJY`sGKkKdm!%*+IEgMpz>K(N@fJsm8^mXz6IU=AK(4UUJE9@25^7vjhF8yRrbDOo!gpdLdQKFt(4br%$LRH_nDYXSIh|%bR!PdiJtUrs~5l!qGCd4r$l3`FSq5^|1SmCR|Zc<e@AqO-9`MZ!Gi;G!<z|P)e@Z6Y{wPxrDtHV4NVRxv|2G#OW*}q1gQ9;(O'
    'b#_68|W`k<DKfEowrC+o@|&>ue>Q2q4pVSl|-=fxOcmfk<qQDdXG)G9sP?00J`wXl)zkh5;#r=3~_CMsALZcP-oB#RoS8PY$}H9Pz>`&z$1wb2o1PeK|wHA87kO33}nb?zvxZ*zdo!oS?Mcw1rHd0LboWfvhe-VYJ!Wj9g{FMDV1&ExPXsfHg!lJX2B3DMf`clot`}-W4kfw{z#Rk&$6Mg@8D#1ekB2Bq={n*bh+$d3p){CX'
    '6k{h=ibm;GYB-STh37M0)EKo)aeBU$&@D#({wR?Rh^W1zjwnqqYAJ>~E~cWZg<yw^4*NUb3Vi2SaT2`XXj&RMN#EtLL2?qnlEyvq3M9EQ!~;Z!ao6wLVm)17wcv9a)}HB=4W3(gI}gCh%oPdNGIy062l`rgNkR?dsgKXy}?Xx9_PIqo>y1qDNQfmIMfVc?-Di+lT6@?ObptFYKd8$M3<1set50`Mueei;;(Lh`H6D{;+0D1H'
    'MtlQ6-(9ajc^T`puTSa;b&J`>q5{kI9>&*BnN^e2zx(eH|thk`HjCETvk2i*M}&z16w&Tzx{#aD;6_9d>R%3nV}Zg8qKYXq9r?h3EqLpGFuspsY%FUJj1n{pE?VEaMf>Vu|$_Y70@@8gqu`#9m>O#tF8Xk9y_crBBvT>y|zs*b%J`^6D~7SZ)`Lwd!yR3d$d1a~nfUCM9~ZilXyJ)8b^(P<wV|F27^iR7IxfO?cw(^sGZ>6w'
    'v+oqcm!a3H*a97AFgrqT*@?nRq3z@=Ii#!y9K?IZu2v(yIQY&XpxU1C7c$hYHG=cEf|XPoCQ*Lj9TAeP!7i5V)fk0{Aj`0ltiAAnIJ;mSH>An_<kWH^}h~?+CJ#zhNrhM!$8{%7sw9D`Gw0Y=;g-djyi%r?I~a%ikatfX4r#2&M<M)Mno}>NxL+_G(quZ*Yxrx~2@`w7D=&7bKMnh)&jq!{NC5j16^7de|Nco4TbZg(-?E?W'
    '8BdFEaOo-R=vJ3PmSndqsM%cvL564bofnhTcA|7P<YsXaVM&^h)RbA-Hpw0E+1bTbsFCCx3G+cabhkrWnSMDJjA8Nx0`=Hv3JDATCJlbuI!#XE)I%%jSpkUT}cr9?69u&{?2P*Wuoa4O%kQqbzhl+y;eYMv?3rs{AtZ*KMNDx9HS+F<ANc1ZE}2*Bh4Epd)b6nJruO+C8;*1aRel58|bEdr3ZYPv_khSc?vF0aP@xhX$IYM;'
    'lBw6*__@ECJG2iMD#BP$^<?T`(-l)4BQuL}!yc?F~?l3Omkp)N}__zY)^cMUr}QBxlW6s)z@KhWw6v92E$Y+)<91N2}Wc8k|BB8N4TIEJlJ1Hpg?p4)s><g$v3L@mL(`CGS}fR7l1Rm`m*+fu4O+2XXM(+21MUJ+xQ2_Nm$F(Ld2X801v)7<Z9F&xO{n(ZwiP5*s9Bj4E9oIg+@!SFGmRk9hpbbh2fVJ=A%)zYqTS8Krb*&G'
    'TX(iys`Gg0i}^vc$Kv$rZXH_t>8lx1xStcA1ps^X5VVeLwaz!c3?quLqK@7P9F^IV-5C*LR{a$&6CYa(+2w4Vm#q;<gV|+HiKGn}tPSE%Lt~Gv%~JN4;-Ubv4cKrvg^UNK`r^Qa_;Ro0x}+$)n0}2^CvHbQT91@MIsTj$6JV(H6)}eQEsl)+0eLA={s`gE_#;D{Pw8qeRUL+it&hA79)E<9!}_)OVh&4hw?H%P6c3K;pWGHE'
    'UrvO9dVJdJ!cST7A5l|5`{w-Q}Uw6D3G{H_{!qx(}oGe+)?qu%@)hm%=Oo`<iplxo2yh3MWr2KhK-1+uWm-DF-$(ucpr>_qq!u)P+#KCC=WtQts|@voxihxeiOO|JDgE;Zda-Sa3b7vmgOVV#CJ4<x88Gb-M0$bh&hU^RWI%&KHVR#c)d)AkumeOF78Wav)bNuW`=HbeT3rh<9eJs-)*{>+M#sYicF0XvWDB3>cUw^fRFdA@'
    '}IJ7F<XkZEc=T(gZ@)oD)t%7L_h$7Q5{J2w2ZUZTr1O4XVli$iJ>}Fw|TfFY;^PsCH}=&M3kgH)ZvADL?VBrB4ku41@cs$VwY9EXc;1<euM_a=g!43n8wiXyOemQ2p|ZE49~5W$KE<+E&?~NGnBoc3W%5iNIr=4sZ%``6>BcUOxC2Fm-j=hk=Jiu<a{iMKE|Xm?Sv@^tf@X?XfTNb6P8V*~+5;Qa0||a?~IK$fx6dJ8=qs8w'
    'SM$n2SWCGVP~P_t^*bBQ<ZPvB*+i&uRg%;;c#4q}X!WxG@OwiC>1W6i>q@`*85TC!?bph;(2B@@xeusscW^<~buBd<aoMeQ?RcZ@y#>Qev^s)h2ONrB4=Le60{3dHI6TxtLYX+7=gknP>?J-g}#c8WYpXy1Lk~ACYVeNw0c{aD(TchbdA;#O>~=dB3_RXs>}FF1_+bV#c99U71iDW1-MGxClGxnFwyt)S@eX81H2Ym9vP3Oy'
    'hnSFGA2wNv9_ghU_0v;2x!`HW;lIaEP>|WR_-F^VmFFR+P0-3mt<B+{<Gb@vu$&`&G`57zp$|GzMa}Txp-rzf*Yl&QzP?3FtZ%2+fD^SWI?<x1<9GFJcv(Rcc9hLP9d_tmY99eBAeTIwhCP#4*||)-A(%L*#zedSYfYjgkyr(O(;$c*<*r@ioWOHnoz(aIn08?avDEQ6K>znt~de@*anbG1WyYTQ}i-_i7F&^8dB*A|Ycd={'
    'PFz!o-d;TWob!=_p?S)n_ft4ey%vMDE;{=bj3B>L!y?moFX#x<T_X2jpU(I9w9XECUsy-(VWBanS$YlaTj%8fRnPJqu)Y#a^c6(fCLwLB*W;3mD;yl$(SBh~%b)f7<-<-tTDO*17q$T%ijV)IjAP$;5SCFl@HAWAnqV7gKJq#LJym=5%W86*g}3XRrVu))u3QkRk3n7(_w=6bBF9)av5Y7;17nVB6r}$1@R61GINbaw|bW4Q'
    '7zwgc)5t*##h$+EVM1Q|Q>u8HRN#Xv^l%2fkbpCV`;!7R@S?cJmDH8k<Jho(oqDIQXZQ+{0{g==H~ABTa5+SMk-TqDa<J#na`H`;pKae8vQsm`Kll1W189$v%vnHL*Oz*834l!4SAz4WO8iOFc~8J;NRe9Q>fm|Iwkx5uf!CRyv=1=+{umZA_u$vY_1@oIPoduoeQ0tM?y)FtOeQx7LJvoD3b}Ed`Yj$w^Dx=fB`XuiFITvx'
    '8u6pRxUX9^1a}o3YwF9XTK_mtflUy$CQM>GrrkjJ}E~UMI!#sByN(p)*hLWF`&8+iU3u=q+kz$9H{k^bt5y_!q7X$4<(~j1dfbZVq4i{2~hkw^6Md;l4Rgeubtn14zfmij3zG#(T=j0fRZn$>I>hr8Tmh7ETm=aw#9W9%yVk&>Z>T#R!6O)=P|RCb($H6r7xpyv?YBs_Y%X!Df4lc8uAc%pU`Qc0=HMUU-k+&KRvSp_vHYqU'
    '^D>&G01s3xM@BM5cL51y|;^siDkp-{*NU%av#kP`@=Q<tgCk%PIoYkk^QGO9o<=jGa%$N{ZGK#vG+EubTP>t_(iZM6@|CBS0O!a6QAst-AOo|FD|6c|~$&Ow@rTiwQhzrmCp}JfruhwNIs-OZGoZo@exDT<}y&VC|8q-|DSY<Ys%Pyjh~R^)zwes^5J8uue30f1E(%?ZwR*$0iHY^IQ8X%i#B88j->R`+Z7BresC_pA!O$E^'
    'bF51L4bGrRhXg9Ob4Fua}6%^};MA9O+*;aNhc7Wrt^jPt(DnUAj}G?Wr;!N};307O$F!^W)^~A#VnfT)hFjpw|o$tA!f6QAt|WSXK8s`gu_Wq8McAm(9Dt@!k#8I69;Cv>R@2kXD4*{Uc7yhJt<k!mwR$tiK#ipHp*x!%TIlZlm=i_fWu@GZKNo$n{X{Jbw(4@BjU0BkjllnCkI6q&5S6G+6np%4SWCNd{YbT^KH*eR-kvzt'
    '+m8q7lTA3OTzZ1#E&b&|FDm&~@t7Z(bLI*;N`au_akFAS%U&6H6c|oiOkWW6WB<VF-f$k=(D$*}?&@jm7do#X|baethoaUXxh=jlDN=HTzc{FSj8EtlgC4cAfBi?)T@`0OZO4P!MiYNAsfdDHX}RkolYC(GE+%nfjN}ei#C{{oK4MWZBWQMWQAI*U$nOEdULD{ZZQdr=O6nhYp9Nl)-6TlNp0bW3KF4_%ykW-%?i$_yfzj@I'
    '_xyjwh=S#}6Gm*rLNo<}m7EWn3;2jwgGj)}PDLNnWp*V51FR1~f-4*@ZDXMflu=M5Y|1Jh5)GtccpH1k>Hs(N7?(1)S!36mq=agmJSDV1Kn7+U-K78DwV5UL>ZBb`&~|`ghpukSRM!YK+W>$SW?gh~n{?9c%%7aQR#dGJKwG1_-vtsT3G;pRRZySGEn*JO+Ypt!s3}GajK=kjpW3sE)#3-i)>fa82Es!D*ebTA-nF4Fi~Z;s'
    'e6;FpJ%wl&J)8N)&=LiDI6_;k(!D;AAt?I>iSODiim;W=s#FqLwn7=h{t@pdSwbnJg8d-6qJupj$CXUq@{BMtai%Q}QQ|@45fMnL!UITB_XjBj<^EBC~zyf<_l;C#5m5opUD8>a#inJBng%4t-pt7mhYmI{`Kqe<-L{rYyY6O^rh~?J$5%%3~>v4%48OWLZ!a%?YM%7V=#1h4AT@AhI&xS^wdhB#tGof-QTBy^yejICbru05'
    'UzD$apV{gSjz-7r1WyLb#;az!%z^C+R>j;n9<#=xsf<em0+RA=Zy}X^NX>Gkfq;vo;sb6A56;nzInv4afte;MppzFdQ-a>C?~0uA6ZT?!7f+K|s*O#jG~=-GE!h7bt2mP(Xqo^efHq2_>g%+|7K!SK5xw;{I48e_HAO?GN6K=V29pW+^kLTHqmJ*IRl0bE`*BHI7~MDH~glP_U52<q<F{R9CrH`ke)bBQ;e-JQ@1?=+2UTgL'
    'Xv&bL{@`{cBUG8AoIegMs9+{pZHjhAgr4(0)Y<r&SOAH>gC4Pz8~RB4!xs4OhT|J0n74>^6NWLlB<NQVUTWM_22WHzd<bMbq2Cbj)%~v^%Rag}Ne&Ic$zhvdfpx_fJL#B?qa(NZ)Or*BU12vtFxz<7*C?_|Dl{>=!V6k-@+pBlo|L0f|(ggO;!sZqq`Ld1#KHKX<0HGey|&<xl_Cb(W*WPA5BuN%DMDENx2zV|6hn@IgXcZi'
    'K#|xILkkV9a~sgF!G|7Axj3B%%{b;g2sIe)fp5EEgl8u^%k{rxSv^E7_oi$w4h`ge?WBkVmE1-P3CNpEk6W@5Ce?(FrT?{ssN_GKn(RaC;g8n!U<jNO`X`pt7$Fu2rcFea?511`*I|#Ex`E%AdU>Dvk0rnv)ovc@9W7BZn&^A!7@t!xSE1O;w8Fb;6u<S**rr2_5WT=S2;WtPU=*)2p$ql%?4AXwpC_2owtGdamuW1bM+?Kn'
    'C5Z)}fUW4s}w;8BZUCK}DQw@NnVmps}3+PvrRwBk5<JU81I%GdGPJ_?|F3UN`m*^f!rxv*IZjAq5}SuvQgZFeVAYj`<e|HH>oh4%gj1V~lZcZfAUc&4H8R6AmqT0`$ndeCm9})Ii7GH=rLCP|!&_V_UvGvI1Lxi7GPCz)`~otZsNL9v>)s*dDE`s?uPA#+~BsbWw=d7A@<=wc(`;6yYwx-(!{in55Cma%z?MANsrDE*qs%73'
    'Ibr5WWD26~&vq4eo)5TjH4Mlz;}-E6!0m{X*3H{>Z!+#YRh7fn&NC6_qvoI1;Q5#ZvpbMDH+Dw`NjY<!woIYx6wGWxl!o1qj+tyzzw9hp04KMM_nJtKb*}Vng|6XUWS?0pIwk*C5r72<%c(S~P{O@9$4!QX936x|Py85D&r>>fE`ETx|+Tw~*h(Wq<^SuRv^5wgcHcTH<i+x=Rf>8|jRb=?PNa5spalIMiqAI_^G%FWJuKuq'
    '93w+l{n(Y?pXZcZpMMSa<s1D-$Mc)93cvwc)oUFb6-X_gj>a*FImD;VxB6ayPi7gFnlHYS3#N4>*t!T?OlY*XRP7fvLosd4=BqWjAJymoeO>pGAZ@I&YkFuIiKZVWbB|>nugCAUASwE@0R{jRTl#eCyR*bQe){6-LuzgxuN9swqe9Bf0yBDBSsenr!@@adjSd$A6nJ_a!l7p5=xa&U#1{jN5774X6#unMfsW=Bac1a!d=z6P'
    'oBrYlFMKx#)_tNYi3UpmWw^zmnl}ZE=qH83>MhLY4k1Ge)|L;rssEH1g<tJd!FX@jq^t0Pi3jc*Q1#E84q407WF8k5vfFDI2rF;R>xQfg<95<OZ^*g$n>bH0i?_#u4I1LP1>BTH0!H;4<_?{K7moytrOPj))%|%3jF|&jj&xCIN+ZP<fWCh=NP7nluW_kE?X9rls7OzeYV{G_j(clK{Wi$O?M-1WYkvf=IG;@0I~g+4Vxs0Z'
    '-nNL0F5D*5*s`VW43hw^YPx{;E^>_i#oehw}N@IaIeT@4Q~$(>7-M^(8#>hj>52E8K8}5}}m#ieD@iS(rwWXS#k$i%ua>+$@{?-M4kxU2?DXH%2!r`LN!eiN*W|WhfKsW#r1N4C}XA+HGExyp}}{FGY}e`al+%W?tv!sF52&S%8z`1vgw4)K+!o!XZ=4OQTR~&?}IN+->AET4Sc@Y*R)Ql@>QWkRi8@bue9ZJ*#r68e)pQ*K'
    '%c#fRa-*8(>a?)E6>~e*6c12lf2C-V&>UMweW>a;z1JYC0guVls9On-NOfNS~Pg(y~)5XF<tEqRb)2P~b;_%rG<f;5R9qB^r0+p3cSq?y=N_&eb4-tW`7}>Oj_{hycxRIoM`B+C;v|kIy3~J&;McgJv`~Igx$+50>gcI6naIlyM7?EvD*JHDBcbKhp|omJ_UgJUne3G)<0mB09A`_>ldvH6D-<k+_~uqUR3lLAP#i3Y$HNeB'
    'T;g+NJRTiBFF9x>|31hL<kH!_hp%ZmbaB>jcagt1=^^8_noG*!I+E8i=b#s!8(HX=c&}J8-|zT@$N>WDkqVQnIm%yd1jwzwra1`nxvwvXv6n4%VUH^|Mu!Kl`ejm=a?#2kv#CaW`@O01D8-5}NEd*f+mICUPT*$kbgwij`LkdWZuw*m^8Au`upTUeLs_%Y7khvJ0+JggD*IvmBcr*?pP04Yvq`wJcO?Pm}a4oS(7~SeiJ-19'
    'sgL2;_!=by-%M3%_womQv;W4NQyHbjC+V<%=vwb|0S|rvZABr_FHv-daf$2<|y9Qu(G!f*6$kK&*&`e0AKPaXy;XXn+Y$jN8w6FysOkE0=Xd)R!J>^ZAD$M~yky+fn~Fx)D>&ST<^5Vv8V1zM$-3dILI8|FtFW8LW2&TS#7fX{xm>G0wgOPrXM6$t#*8vi3y?@XfagmPzYnPt$EC(NqT=_u}`DVxvWEO+bG^Wpf;liFY<bHt'
    '`{@EX`Tv2#~HC2WiCuX+FD|zQi+ms=Z2t`3VOL-Emn7i6X8!89FA#b-iq<;0#A8MrYi9s(~z~_RFSJ;Asr&tc(3>9_|De=-RV%Q&FWMSymD?OqP)BrV(A?#WFRRiE3<eiu*GvIaj}3VW#DgVHt6{ajyJZM3cHZ`dIIhJMi4pV8-gm<3db9i-q*w*zkJ)VwA9@$;+1*d&t^MVgqxqg`U?9vv)06O`8YhqZQrF4fqKCjY0?IYU'
    '{h(bjGiM%2y1b<Uey6)jYCSw=^)id437rOhMfC8PF_lVo;ibOVu5TALnM!Faow<J8~y|15+`h${mZ6%<wWN;-lh^_Vyhp2<jzgn{Q~di%nkp0o&CzFU!`(zQ9CkA8#q-9_-sbn$9)(yQ}?X4}E43PBi7qlgqjeJkQFNT$s`-l(8(qEv{`yNMxPA113jdoGvF5YNi9s%;F#<_nhU#^s|ej>$ZQC+3f(evd!2cegjN^q0o_L0-'
    'v=ARwlWv%%@-MOspIU(4LP#oXEq89ec>cu;X%Qa)QOZkj5`*JjYLu1>E4=(1;L~EA!E@XT<gz2jvg7d$ptLA5VuEMGfqlF$Z#=LcwQMEX0FU>hQ=v7eeVrBS&w(Lw3+=O{HXNJf_(VtvfCtRRTMb-w$(}5-OML(_vNY#;!I$U=8z`=g+Pv>DI9cxyi*ZycB(0D^5JuOfM-LPu^a~xk$DryDMSLQcV0D(o91%Hvf|72SStfvz'
    'ZZ<j0okL<6b9M&1)kLD~89NSoxIJ796yf(L8bF>~Idd^Qb1@ih_9<sdJww1p?ququSzA?6eI{kEeNz$u>>0nKD$ri&pCXfsxyM(cv2~KL(F(7blj@>30{2&b2cwn|#Qdq5S*z$?T-;tBOS6L!-yu>Hf{T>q0m>Q`(4Dv$-?pXq$d29O$c)jsFm|e+d`LGNpcqlpOTy53JIf+_zv(q+-^i5DwK9NQ&*h=ezte3;_F|Pvhd4#j'
    '4iwBG5gm$onLkinl6SKmFH>t<81dcwA#tHx<D|0~BBtk-5nQ`h5sc+oT?NmdCxrv{0ri=dIyxt2dQHNUs?L65Sh3*9+u7V5BeZwdsE7&FM~ElUg?}k=luI$f$Q(zcCh0v`iq_N<+rUx+`YYm!fzA_{to-Sc09UQVr_7343L=^9GQ1yz((FB{ffvN$D^=oc!?&@md4~j70HZ6@iKR%f*pya_Rz?U{MU?F)sw6#kBp)+IPxbl='
    '!dCQH4&63+%U7JJsmF^18%mM3y527@e@muhBnwC<1ro6@<i2Pd#_mYfaOQ1o)RaYq>E$!_QQp3T4V95Q6+c_&rJ5JL$Xu7WhNeu|$-qL~`N!tY{|^UYtiGsP)jr>^sQ&rBXB&wK)UMc_(LP1i;S-sx?^gFehB-oMECYrb4tP+V_Ji?7;g{s*NV>4|r2$(B6={8|xZQf2l~=K%ld3R$^lM!A8O||9<_0+wX(P9T2f@GXd4Nfn'
    'H>cn4Us)A<U~=Vk)}ih!pxcRp=cjZ|0qb0vgS`Wnko_!Lw#B;FN(<qQYGnn;6I#4gTTY(f)Oqh|3QDT)R`l%KY{S=E3>tTJOA<n7|E{cSND_ue^NAw5mSL4YPwp3a?s?>gbE%co|-bitmkux<GZ<IJ@$$Bj>X1oqrQbnU+6_(UU2bG5wtUMaMg*V<YH@>0PM>p_DVWV==zcO{9`@oPH(KFG+2tj@tuX^+tuk9aA)9o5|Bka%'
    'xtrC&Pe+@@`eDz}Jf*;IMIkua5jXi;qD7)?<c63y6*(`DEsWQds4}0NrjV(j{jVIvs^o88PF{;e-Uja0A;LCQzTTZ&=7)VCH{`9B?IZ&+C&tqBJ8f4dsO2UXz6z4BoOVhq|*5Tc@iMT-oc)Kve#ud0mW(p?dsV;CO$m+d4|G3lrgMt5d)|LgBsoGojU2yzn9*TGm<15yx~9pT#60<ED!)2F;GkE)~zBnBpCbf*r`r<Da=j)1'
    '11f9BRtmf(3zj%SC?<%HYZ#`ho|m=yW$4o{Rc+)egHl_s&I4Bh!;<@i8k%1p9kUK!4+!qYFWoCqr<*JibmBcaZ{=ZWbYd^aER`K{7-G9CiqhfK0?L)0aP#devW1yfwZZLnqOol~YT{hh?dWKJf95?TCl$QQxZpJbU_z!=KMD1iX^@w_`5-F;)Y+nKm?<EarYdBd4;0hnkvw81-`8u?8`t#E#AV5&&zNqzal`Q0ADkoi(KO{5'
    '(5quy8ECD{`Ym66ijp{>zA55TlJNEW{mI_<~HY{$ybvEoJ6~osV^)D0F|S=8fcH@nc_@UF;rB<!IwELeQO<7GhCL8^KZ;+Ajf1t7FmnL=>zeZ5V@a48EVM?9&g&p{&p8^<k-`cPUF1j?8poGV0j^rBLNbR0%#oa0ls8<*4uIH8!czm}x0bhWa=~XAV!>d(iH#Yn+w<Gy|O&5-C$c7VX_#!~1(84`mYfiS;CZZT5=-qsVKN4<'
    '=PV!;@?m_7^ta!bJrmRFA}r=SfR*&O}a#fK@aSV+!EF$?QDo@8tViLeaQJnnlsV>bE!5J!fh?j)FT^Vw9$0y4V;=<L?~$+g$t@6xiqtU(qYPz5aXfTMdK%kg|j~Qwlg$3TiqRp-kA_>gBU*;=ky>-$XQa5G#I?9^I)*r4j@fN7YvA*TYri3YO$k8T0%z@z*wZYn&|}_W{FRzvEnHRAa$I=x*i$YP&hdFsAQ<09QLl=N++2&h'
    'BcfjnB3sk4;ancRl_y3LcL1#8kQYm4#q!ZHcJ<@Pol};2`=gHP~t9bgeKixVk(GmdeP-Y+(<mR;SCHl$#!CBMQZF%qJ%un+?Y3KNSu0l(Nk^+S72`CUUBSuPookwyL`NZA@v&(AvDXRua#!5;;VDEgK&=6J6j>wGZTLB1M^;U+hTi(tLwoYL)#vb3V})hh6A`<;o%_JU!L^fFNNnvEQWWntxd@%j;A)PaV?=je9H|{9HGiyj'
    'aTMttP2~6@QvC$_>g=#=ED}x;3L~gu71d2FRd&!;O(0g`pC4cm$DM1%F7dH0tOUEYmw-u-26D5OL#NSxbi}X2Y?z${Akm0(iv5_D0iRBl}<3T|}%M`uE`ycq@I0G5HuM<FPFdeFwWq#TyuRSoBkTDKmL3$yY_2D_E#8BB(|?EWC}^MMnv6;jOg~J<-B7p@vsO|0`VY3Rwm}@JUbhW93HE3Lry@sWWA)gZVs@=r}u$Ba$s)KA'
    '&b%pr!GmB_TNef|x=;3V64#SLsOJ0BO|b&f$Xveznv;FobUMYiBtuyMxEq*pRPi({h<xOmE&pH+chhBizHeVB)ap+D9D-GJPwVU}5O{$^4X%aJ&@I!I6EB|3!c+<14;dGy!i$j$4WBygkoQ2@?If8ZA{rgeVq(%W5p|D+dU=zwk;hSck{;r@K5RV=BUh7t9s(xQ!&Hw$eK<Xf7YoN>E-mTsL96kd{CrmgQ%g!?e>7M<>Zuy@'
    'tEI+PF+*cBe@iP&#4q2-;3{ibuuf#x><hyi^=!%h{W0R3<sCm7<7rU1(EdOLv?>&`W(~7#wpBnuEQ?oi5XM1)a&aqjYqY-r2!h)LR%Zl)A(fN#*inI{5}1g_|DrW-7dfR9~<*q#sb5;$<Y;%1h(-C7b>!u#iEI66ah0E29#hNNR%9-pWga&^~6?Kf>3rO6)F_te%5}b=IDC-y}63Ey)qM)n(Qd-|-buR!KGWz3<%cy`D3f+O'
    'Xtjh2)c2;&=*A-eIY}6g_1$)KEjHntnu0wd-p|G#4B{A^T>^HuhHfjt{SkhH1zp7g<?^4CMvCfs6@CW)>7=G<912A7~}g>Q8E~RO!e@Y0}I<e@#l>_SR!Di$rgxD}cpspHa?>^ICB|rF;-sEQf%z=aVnR&5b912mVdw5oX;#cYR;5#Y_P9Itgcq7?K5tqZtid7c~!d?$G38Y_Pws6cU&n?Asbf+0f<sR<iUyly5Xt0VDK^`a'
    '4Sc%n}}vyr!;*!6x&{1SfKXA7-}Kp45Y2^<$DcSY*xhhyb=05eOjwnV*5MJaPMN(I#tVTt%KT7#(w-8=|#8278Abj*~KE#*3`A7bnO!J@)l1a{%TCTmi=>KuL|4>LhX5F$Aqd+9$JJiO-h<=BU2=e~ADu#!#0`r|o(Jr{gXg4?33m-A21qN`*b=v^)>4h8G~+--CQFt`z8U6y`%0AA8@LDZ2<5I(c^6`rt0L#~?xDO^-j_!#'
    'uL2oX49g7VKSN6nJx(eJ#UO6X9@$5KhRJz=k#DCPYq&eVh?x<i3YmaA#W0B0c&<2|^vFBkLlNMwt<Rk0XD}ce+3&LccRY3|^5MvwV}T<xXJE{d6L*r+aG%ks?Z)pKzE)HR_te@EAl6VRn+_ZQXr<MKb}LfU-8}J-It{qFM?%{&S_??2Aaca}0fh*%P3iaJ_)$f3pc<1QZ;^4v2S#69TS&CtWN;*!vT=%f5?)8G|#0M(H9qjx'
    'r+l9vKRjKTwn7RSBPmV?ytw(Y=5<h!@OtpY2N3A|Y8hbO&hT0puY)YbG2gRj!Gbyq5`Kr-Db|B)KeWB;LPbLV=wTLxLFN9n|TvL2p#$WF=rO^`vjxwVu0|U&mX*Tdy=^bOqwD)2_-gHNso5cmN9pR1rxDX9;zqPb>Mt76Qml?onI-R?%5(Gy}0G9A5`Iz7LZcE(gY1PBGeZAT6r96E-p1HRH&S!l*Hafkw?<*0@S|bI)+t7k'
    'tx4;WU32ekG=U-e6vQzr(eht;;IYkZIjudoEUI&K#7Si%;SqyRt`&oaaIQuB%E?+x(5~--DQ#4TFT4Kp>YhE2NFLhAnDNq-V-^RZT`e2~->KA~dTNBb0F~k8g9mflcz8v8q?;l(liHjsw+>(nl~3;-_<D+TIdMYRyRTlm@MC8);U{e-wJ)M{^r_@UZ!QyI_oP{jj}q^z%V32D5NSFGvYa0}N(yrm}3^mWpC8SRPkOkVb<m7v'
    '5MB6u0BE&aYXdd~m}T1UKdPyvo16{q`EHD+@?EjfhrB@XQKb$p!2KaYPVnYw@WnQR=c+0#&|J{@RvhiXXV3PsXmao5-xHA*iPftG^U4<anCqXy)hE#JRmezYs7;QX0aoU;S8|#mzv^TNb_RRdCRpzU2G$kcjHyfJr@88nSFtG7XEdg#cf~<W%&=O2sO(SPfOxic2W*L6xA-Jg~QFG>Ic0oj6#08t;h&x)wIWw{s@X+At{PR$'
    '3=t@als;kIR~w%s(R1YaKt{Q!+vzpHLL>k-C!7#X(cKF-yX*%9tN5BkA<oJ=3Z5v%r<M^o$=K0=Q8+9>LnXbrj?QrW=l`VJw*(LAtr+BQ@y`A{4bR>fK&hHyUQh7%x*0DA6AoJ#*G^sVZHoz?W!4t@r;r+S$K%EuA$;Keoz13j<gG-4g}~`RL)k-AL#5>3d|$Fl?WnUr^v)EBZ%h<&kCrg9Om)^zJLd`c1xp$~zMh63TkyCG'
    '{+7UJO#5oF1)Lx0q5JtM}`tK2gHxn>mXe99F41zbH)ROS~g150)IzUwX`68+t{Ro}gaRvlOLd{_F`uCxsY}#<+Qw@mm#8;hkAhNxF6;iGFK$Zw%{&RW3e*BW)rR&tIgpupk7U7mu@AcxR`Cub~sXv6AsO7v_?X(S7Y&RppKXOr_L*FwgDoxa-2_)}{ZA#s$Fyt*Wj2Bx({fR{`F)IryZ6)H#Rm#edr^#uzVDiA};#qr?=?e3'
    '5`eM6rl13O<*nhjaGn-)4rHQ^hH)4~wao^gb;yH1+l~`(wa*a_A2`db#9&R6N~KUX^;CPO3ohD33TPT_J0Sil3YSi{h_D(P8iMpUE4%ejq@MucW=j4+f)nTbJvoY&eoX%CoxYgK1O(CpuyNBy%3&MU-icHd0NPE;?I4zJAPExr&#i&`Pn>MajzLCEkhm-e=Ie6_kHwv?tFaC}3!frA2m7m`#3;ICA<TC09^yt=(WSJY^f@Yw'
    'u}e9Vj@p_Q}&cY{uWI=E<I4`M!!fMJcow0w^X;B+j*nfQ-)7cEE>*VwUVrHz)&KjuKew6#E-nE1Hew&+kgBe#Dt0D%I2XD_%Tbdvfva69Mfl!mdAEEGr1uR}3E5T@O*S7O*s9MI161so1d(2L<qPJ_3=GU`&9MQDX<zC0h&h000003fpURTEb)J00FI~0*H(cZjQE`vBYQl0ssI200dcD'
)
