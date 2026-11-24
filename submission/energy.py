"""Active -17/50 first-pass, gain join and starting-majority proof."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial, isqrt, prod

THRESHOLD = -F(17, 50)

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
    ck(U5<=F(73670242301,2**39) and gain>0 and bonus>0 and 0<cost<=F(123755621548,10**12),'active energy bound')
    ck(f4+f5+f6==62 and n4+n5+n6==218,'face accounting')
    return {'energy_upper_1e12': -(-cost.numerator*10**12//cost.denominator), 'first_pass_upper_1e12': -(-c1.numerator*10**12//c1.denominator), 'gain_lower_1e12':gain.numerator*10**12//gain.denominator, 'motif_gain_lower_1e12': bonus.numerator*10**12//bonus.denominator, 'motif_orthant_lower_1e12':motif['orthant_lower']*10**12//motif['scale'], 'U5_upper_1e12':-(-U5.numerator*10**12//U5.denominator)}

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
def initial(rows,M,G,first_pass_upper):
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
    out-=F(21,50)*first_pass_upper
    if out<=F(1,2):raise ValueError('initial strict majority after motif repair')
    def ceil(x):return up(x.numerator*10**12,x.denominator)
    return {'initial_lower_1e12':-ceil(-out), 'initial_before_motif_lower_1e12':-ceil(-before),'Y_upper_1e12':ceil(y),'cross_pair_upper_1e12':ceil(x),
            '0110_upper_1e12':ceil(penalty),'overlap_lower_1e12':-ceil(-tri),'11011_putback_lower_1e12':-ceil(-extra)}

PACKED_PREFIX_SECOND = (
    'jfZ1!FoVf12yOZsjGp&()=C{wBQOFQ-g31cOh@9NKnLZz9mo%n(J4&h;q0bQX8d$W2_sNa{LLc$x^NVO3u=fi4~M^B{Xth%nnyfx'
    'z>S!nm(ME)kY*Ls^jY`|;Mnnfwu`7}lNF0tB1p~5Z&uB!!&q^*vUaWdyalqf!$kIS>|k%3VL339Z&cr1BthN@nw)}D(-C|wzhC_m'
    'TKU9UJ(=|7OUHKs^zw^4uzB5nAssi6=BoYur4g%<A<kij4w9iJ7fYkiQVZS7xc`obh>;-DMta<pHo~THOVXY;YO;+K0N~&!H({IM'
    ';$j4&0>CE%ud(4NC;11t&@laufwa&gCx{)Vl~KbvvSm#xZ(>#V3mf~a%XBOqT&L3``Qc*0n`80#5s9<sQZfpe8Hoy06%!|ZrMQd-'
    ')>dtiKw#wc{^(L#im!fHG)C4qsudAq%+yE^)VUy7_t=lF1q0Tj!s#Q@{H42Fu;5{Ll{m;u<$VLSY1=LRxMS&Xx)gJ$7Gs_LD;(5;'
    '_W~o{dX~bv)+vF1HR&E4y@h;Whl^^{<;|564*imDtd41p0fI`}>hi>`6vn#WH%KJ!=475*Lw%{89-a@XjWIWQ>2ahrI;i!%fnXsj'
    '$~Z#q%S||vL=1CW9|aH7n;juxWvi{$e|?`PT{e3*m5E&+&fq&Y(-06?2GiR|v$w>O)htCs;9@K)-5fCab&#cT37fWBXS^A<+$mfa'
    'CiD(CMC-@TbiUM%?1+f!8B~BTY3qS~!U(E`|CUo~nzV8IJa6#9UQ-4CY1CD?fWh@@hpZONKhBC-&Q~OiEEVU~?l_IBNEF()EyzU?'
    'dquxBHq9Jo#})GIsq5XgIz(Dy<uby`+1#@MJQFp!k%hXl-EB_19_AHFnN{-T`^Ew$Z0g1djw9N3GD2_oCP+{MVOY+U8B)EIb+*e}'
    '*j5`f>Oq_#<wvP=6@Y9Rk!uOt?sL0E-)1ROZv!(2(TapGb{9;X`f;T{Cvrv*JnbmG!)`v1y4VY5qa>Jo;$R}7-Fc5<W3Q`u%p45B'
    '2|qA!j<8T;`4e#0JO-8YF7^$4znH5Numga15|v&CumfO`yI53#7}~WeaUw6>HTl#%h^66&>vwCvk!Swf9K~&`u6AEZ8FhPqD3&2!'
    'iCvpQYoeA51M($cK`q^|P+tWvO?uLTfmPI(6cCyt;wS*F+=qC|qAkCZx4=?|Aq#P*4;*(kzqXScn06iq9{CW~?7l5m(5YV;E|J>B'
    'Zz#QJ8zYxWv%q((%>@81s`XvekFM=tKZ-3p6niG4CxI|^TdDE4TRmL&Rg)v<`al6I_A)OHY(a*$-DK!N3q0b3RvO66xQphEb`TP('
    '*5<zmz8MHJMLq6r;{vCG_0JzML4<{Whi&za2-3r6-?x5KkC(0%TA&_Nk<ErvnC10hOM4m1PRHMw$sx;*jw#gbk$)Q|ny@aN3IDyd'
    'AP&35O`Q*saezRUOE%`KZR5~v=L_^m!J&2%(j;^9j7UwY1VoOrPFtbXrTOy|*!^(Q;pIAaj#!`gSGd!5Ul%+YuswWHVm((g{igdg'
    '<E_bkh8Pydee{m7hj)B;TLHiHm()a5Zo>LmQ=6>HQFH^$erW_9gA8!_wrsh_nuE4t_G4C8nk_R4c$^m!WEd#I-xxPc%(NmxfT4N&'
    'kNXNiQu;$0vk*YGau0V)>;*7x#Ny9UP*|CQtWBQ>?x#QcLv+UIT)6}CJ<A<*?Kb2MDE4BG@+7ZwMomW+mb@R^YR4q`f}`4{UIkD;'
    'T_A7;L7qsh9=T@X+p3{PdeS>i>wi5XaO~sHRp+=VtSxe>v9oz{UA`d^BBVhQmdfy#Q*OZHQs|&6-pYAj$vT)YU4F|V)=$_(wZ+Dl'
    'dB_;@vbbzd31|il9~{%DazTE^-~qJ<=H08%$rY=k{@r=_;a)}m%?24HP&p+I0cF1Bx;Hs<_i+vsgABKzNEZ_;tt<Z?DO$_PJS_$l'
    'rjH6%9!EV-{;Kgozv$S{!>Xh?UQB4v$TN|uA$n23Hed51<yf(e(b1l=%uA2tAmWx)i-VM1I8JlKv_i5suh9@v*Y7XS-3n(V@_>?W'
    'w(G|il)iyjxQoU619X{^Y-b#wOZO98`bz7AmIcf!H#v_M+aj6fu}ls1n}fM(Sonw>>){mvkD|Ojp@K%WDG$Y8@w|f|!`k?z{b>k#'
    'S`@kAx{>a8!bQyPqxx)60{3_%(@_&^iI%=CrjYysx~|^z;>Qhms-va}BEsC%o415`PIjy@SVgj?3nx~b7h7Gr@{JR%iT?gJVe5ff'
    'kphaDaUSQ5)g)$kG7OBWOvm9$ci%l<NRQVbWjm8#lH9m5yc)=leAF`-JHFH3*A<_=RBIRNrW{W2H`52MLyvFtmz>~td14!z&Y0f('
    'nSZi(ZQCU$$zfI}E6lGGbj4EIHtgIEPpFAdj6ZWDsa-?pv~Ml}oyPge>u>{%A{F#|-bE!3=VQLSf~yO(uj^*+F+tFDJLzH!uStR)'
    'HIOqU*n62UVFH~>bLtz)7z6#SVpSCSd3<;G?P3w2g2T3|zFn>u9w3p8-7#q=0ZGK8NpKlC;<+>=-p2+&2S+elLENtSmMz%!NhA5@'
    'p?09V@9PK31Uy6Mn2$6A%1go;dAh~Gb9$&n8?K<{u9H%SgDI{v|NnH-g(XiIKG%SMu>*~}DaM<e_5gm8_AF^Z(PKCXi*G~Yk?zUM'
    'nz&{mxc{TOKf!*pir5-u)4J+ha==fIF{0q)BZ~uWuKKh*@o#6A)D~PI!6OhPrgDo}PQ~G1X#{RyDNCrnLdon+jLCCwATjLeNE`&8'
    'h#02k3kYP!FY~WD9Srfz^&%X6d4*}5ZU^*_?O>$WtR`In``EAVCET>fRY#;%c$5w^puBB+VQL9j#a>UkdVnYPo6gV--L@t1BUh>2'
    'ENgxJ<+;aUdC->=5nXRK0k2|Ay1?L$arL<=5{S2Z!cHvjXSvnYA4m8D8|v;>Z<D=%dDlBsX(4k0cx5m*LG_n2mR5K?aUQP%%wsVz'
    'lbePbB0qIFS#AFs3fITA1f3~dLI8rSdr81y^GR;`JG~T36gq8SHU8@1<+-V<<U?GIdFMo+Y#vjWoB{FZpQX3nP3?2Ra$Wp=OrQz0'
    'XY{O61wQ)K-(!^W<2{F5K!v>e$kkh0afC(d4LkTW;bq5Tt;?Sz7xOZ$S3|}rwmQtoa340cAqRU%mnutI0GuphOMKqA#699e+M<9<'
    'd*~F|IOOW^`T_amRdA4YvG7QOIVV6)d&dtuChr&vm#zP^qTh|)-SgNLe@;JxVd=T~c<BM;%oX;a+zk_7cfCrYs_7;PM@DMT`kPmX'
    'pO*W(?1hN9cY3M;5umM3jcvX*4u|9;i%pE{pzDfVQs5bZ>o{mI=;ZsG7%phOSHDtynwffz!DYwb+(95owYx{cfw>HfrFy7Dw%d=w'
    'ZmY~5;j!MVNaKNZ09>6#K*wP~IigNUe!k8tsU8~o6UzNbW#JVyE(SknM#P*7*O*4vL5Lz+4Szn1THGa2)0#2M(}_wE%$&<>7=C@0'
    'W6Y*(N*X?P)2xgzdG`_}G4``ant1j^zO_Dng|d=zs72&&Ymkc%!DrO|1KMt|i{bI@st)r+`u@n+Pj5W>R)38O2)p>Vm_YLrcLwqo'
    'QbP^vHJ}6=j?};*QSzhhOYNzHDpusU!%PC~$i?81cU%+dF&)H|;2h9*Hk-}iUo|kkvM~yZQ){wXXk2A%XTM)O?<QyiOUU^6v`*<7'
    'w;vPXBaGaDYA<a>RHj#GB8U*J@OdG(E_Z3e`Mf!xsJTQ>6F8mh<n~nxKcnI_U%Cs_^AibYRDHZcP;$myY0{wo6(n3QAkR?qDWbuD'
    '#ny5HmjhC2iQU@;Wg*F@Ad|-3+`W3;)j8sCp6YfggrfQhv!0O;0aK>Q;;Q9I6GpLvCs*O*8zhxmaxld(>~^olm`NMZEo=IcTg{{S'
    'p}o!c=Liz*c;9%hXvJpSP64&@K%S2RBZr%R2sg%2=+j`b?yP)SdirsM-oi*l))!i+*e8^DPZ>IatUN&JKjPpldoXO265Q&7iaT#`'
    '-jnI~$Y2!2FyMC$<<id$&6PPty>Zu<u;_R|>ibWM19mKx@_?DNW*Dw)7Qw7u^o|i;NbYhq5K8y81~oZ@0jIGWMK^qn7c(WX<(FQz'
    '>u^f1`2yJrGqQmVSV7G+j%apeQR9Lp?D`bH%vt^Y={{&byx|P8qdGyS)y2IVhOo}E2W;(1Zn}CUZU50UrQ0|-Fy#PdzYqRFe(1-Y'
    'kFCEnc1F6~xyNqKYic4`Ls9}Z2BFx7W`0PTe=EzG5N6Sp8_0mjVDC0<`%2T_WE07h`PX=nZrWCEqiK>x6yn&tTP&$1FEA{25zcac'
    '(PSh$<h<`qRYw%_NlDdRbTmg8EP5g}fg<o5TZG8Y(p9=A)bh=15oPb*-bg@>BO9Z2W)kaAI#2-=>Lr0KLl1C;-f^H)!`5h_ko0S2'
    'pWa4D1`->?iVl$}@=l8FW@7WLj^Ipa(dY&l7HUN4*ZWU5jzq@`TK@4|{QRX_F7oG=MYh0UPeV2De5nk}S~kElF?it%S^wxu87m*b'
    'F>SbB^+W%Ei!1o&4x_C<47Ia}X<!S`%_yeFdm&K6$)0+BD(iwkk1@L^lz-V@$g1j`<ktWnm&H+{(}t3X8oh+K-whxi_>wRDAzbNu'
    'SL>C}^!s2{#8z>aX+H7}lsEe|DnrC;L3O4&MOTEZ=kl$3rr=r6ryz2rM_@NiC4|b5>(S<cR+3+1ldz(3RI`zDgC5nGAg%@SB_(w;'
    'TOv{y>4m$fT{AO84JGwC^2$acJ{nCY_^4Jt{Cvt>x-p`yyz5DlzcPm*9m*6-;vmcWAFk>Kf7yNM`kD}3KL_09C$m72sdzunt1sz8'
    '<6OD6DB(FDQ@+b00B`8;hJyL61{Ccbi{{)oNzsr+MEby9-dCm9h2+3i4G>8TG{Q<$#xSW|A{5TiKru7}B^eYUKb943)rqm&E^N;o'
    'aU_9wlWO%ii%K1Vt8(>wutTfkSja69Fx;`1g1Syyr~$m(%q1)9_<p{wFGDW83{N-*XYoCJ?i4+n5D3+zYJa^;y@bQRA;~zZyNm`$'
    'wlYofs4vLv8rtNEn`=@_Oav%M{Okw`c}4(&YF1PwCLKc2HWV$9j>j!;Ym0H(G;0ezMXdw2UunadIpox$Ig#9TmS~|vD+D{bf?FIy'
    'IT{%Kd-i|!+KI*vcm56R8yJY|u;zl9V1OzZqZZYxBEKj9e+u+nD@w9o^;@N3Q<th-#xV?5k~~?GH2n6~(1RaYxV=fVd<Luq9`a?x'
    'w!8s^x37X&>u>n2*J(<rJsdW}Qsha>5kOSuDG0x{sqeD=cFOu!dm*GKvyR(N_+<NDF0Xfg=72ylAhU0^Jjq3RglI~_YG#MUBI3Y#'
    'kd!G5Kmqi4ZE?)Kni5ef-6c{{`<zKOY&V_x4gJKAG%olE2{EFS+c4#CHP3-OeIm(_><K$_TmMLog0irxaDqNS>bjuDTIlX2F5JOn'
    'd13LAs?6<;HO{J7{$@P>4gxG|q3d{&1)^8AfmwpV9cp@Kviv}>wsQznSMPI3v87qPX>4IHS}jB@vk4iGsw~|Q-&h2?u-azhs0t5o'
    'twczyy%gI$S&2PN<gr#<##84(O~RNj$*|f6&VF7n5W_-MNY%^jZLtfB5UAPsE&i6cH3K9DOk$4D2$e@X0&rNCs2!_&L<|*Qn7mPt'
    '4It4eoBxW=(waq^q)eCy0J@3$ZVDdI^dB(KTn{|yDUe6Wf#Wa17O&hCw%{OJRH;x2&ts=(A;TJBBv7vgU34zNs8p2pvb&&a7KI{G'
    'TKcETL<nJW?-D~wGPFGwggb&Vz+6Aa#NGN`@dOYisG%ro?X~(OQ}g4;0`=E>_FVj>vueOq2Iwt+rUVokJ&mkww6MvSbkI0<Y~)oz'
    'K`(6?ET}d|rKs#S0e-d9Wvw`St0q+1qevnBidP%(0_t|qpW+<irUT*N6_HJbNKl?!&x`LAyDk?swKU^|uYB?-gtrT=<46HJx#Qy4'
    '-p46+!8`Y-UZ-)RTc?fU^VwQo3l^mV4go(OVu>MRvoe_AxA$Ap<K4cOeGAeA=Of@k)Oc**#RK~M37~lOjLO?Y5;xXL2qMEKU8RlX'
    'xDd;;61TUtsiLuc5HGF-Ti4dc{witx5KO86Nbbh{J17z*k+Ej|-a*qn0}*(Rfw(5O6$WSPBKEnRIkH<HaGV2r)8cCU*bj&GzLcey'
    'TQ%Rv_$N)5c}B$(&Q``Y0bpbPNvm)lDS>2<*Ds6Djc^qUU_bw7o7+c>dvUeMQ%(OGeYjaJPTrNeP0b5*?=Mhqf0vt=a^^Vk=Lzp>'
    'Bc6;7AjL==vA!HQW@uIvs^mDVBb|${tiwfq;XMR+R+i@)`n0a7!_o=iJjortx>1MCabtnaCD+r{Ac!~p8BT4BKJ%=7Z!^(Kuq&Cy'
    '<x+%R6u%uU(-ZgO*oX+rzCW}273S%@Tr>=|I7aChHfKh@EZzY{NSS36i;*FL9k9;g*f`^yvjx^5fyv#8MjMvzW=m6U%K%}I1=6u@'
    '&+~k#&gnpCSm}QAQT(!&Cah>@PCHix_TS&U8ahp%_5FX=!av6C4lpwt7^2RE0c&P~8^7qeSvEYgutTGMTf4uq9V<pugUY?KQgJLT'
    '?ki0gHnw&X9kzr4D_!Y_`AGKfoA(51yr(A@_W-CZIcKqdpYi7Da$2NxBnn709-*iaXZO$A9X%X_lY07GTuor^zYSyy-?Ovwg|O|i'
    'C7yzZ1%O|K9!R7Qme@xcPz(pJSkARwel&G+3vh;IwOX~xb$?UMK@>dZ9PI#`f=&?0>~(|2$D3RV+_nLf{@Y=r;sY?w^Ng+q?H{?4'
    '7vG+1XayW$4?Q*Q{dy!=Iq@eq=#Dw#Bu*Hs;D1jbJP0J;<nm>Rm{cf_x6i9&0@5oL?&7hTZmH}bB=SG@OK%(JqS6E{=O8#N*7!)N'
    '{a_V(z_nEVku(*##&|4^`g|gNQdX2Jt?v(~MdZnl*g4M`m8|Y2-Rm7iQiA7!eSBk|k06M6UfyST0wU7l?ROJNb+Y8`Fj<IDaQ{fS'
    '7$`fJsn9&S)7!3EGQn`QE)Z4I>$Jc8IHWN+a0HI!Qxu%JyC0COzIfqaX>O}q$z|5yx&W*Q*LhKuxX!Qw*~^{s=~8XLW1+U6`$H_V'
    '-Ek=r&f1WtX9IKTk=H|6!lPA5o~;3JZ3an~8eBa#n6$eQ88R+&xL7VIdnS)Z8(FH=;WzYlJbusL(cGgQ07CkU*4>ECVY8^c2x-#1'
    'D`Xs}Q5|MSQUJqOfTQu7s&}|t1xP8WvMZAD0Obpxz#8*=>;Wp$rkf>PwtoY!Xz<^te}xOB;aW-u3P`@3b=_0@3VnWFZ?;AC&!o+E'
    '2rxxZX&1gyQo=fv6{SwBe*Cas3~M}K>%8j&SIA-idgWxNlr!hZ4chwWdb+8n*oQ3fA|5H^#&ej{t&k{{is0=NV^5oZZ$FGM*hW0Q'
    'rj<N@-bdji^?Z%^g+)Xxs3hM7Y#rjb$j-v>o*+sML}+4MvVnECt4eiHJKj#nb<CDs(nUexe^v6V-t)GurYLU!JbTH3SxZ6mH@kfq'
    'Ah9;tgCev9gmbRM?uv=*AS}3ZJCPlLc9ay3_(Ti{xo(`vJQ<rO^$Q<fqnR!n@W8w*d!3kODq(ad1>jr8{#0R}6YFBHOk0%whXOWi'
    ';OK^kjv!XT?M%py*=d&{Q7=wqKnqT9^kQF9nyQP{%&OosIyf|8Y7F=+i@X)JAmwLlcdRc$bsWkal0L#Yb;$g~nb=9v=Fmezp-%Ud'
    '+m`ep5|W4@k0<(Q;0)T$idx($h4TO<9=`pJ8bS6E>h5=2X!>XFWyc9HUkhlYC~Jg}P#I~9@WGQ+xP2X;4&k|4{1@a!8QSt^O3^{L'
    '^7-zTBb@EjfzWPkq?3q-XDM1wUfei)l$1zn)0YW>SCK3Y%Rg7%ZP-?(IijZmNzvAv$0i(WrRg90V1Jy==hA>2<DTH<@^H7vri{|='
    'W@7;PXbLgGEI7WI6r?rVo4em=DM6Fs6Rd*X?{g9}$ClNv-Z}7;e5gV6kGrbkX!#Yt34+9K!PlU1c>n;QiR5+}&qgPXa~@Dr8hl5O'
    ';XhW}MP;yk1So@}+YwFmH6b3y{Nf*Ykg#633#YTb0dBMj@SI=X>kLVt;;X&mur$nFn72%LqsB*n%HM>nGTtyEa<|C>Y?eN?LeG#Y'
    '<PEjvKMPdN%YTg@_-#%@KuZ>+fy2aCpB;x@VO=UrUdk~#YeqKe2b+ZA%zo;J-HQI#O`N{z;9EJ4f#Sq1@y8y#&u;$`Eom2AD~GFK'
    'JJ(R?nITOPAB;UI0X1ZbDoHl?!|CD0lCCqp4|cPAa7yb)ItFJrLqT>#MNy_ZKVybNTz7KDj?gphg~I2gAUtcA+nquE(E>9r)G+*E'
    'd?Z?+>AtF6UzANM7nC8tcYgE4)Xbl#kG@!Z;D3IQQ9(a3MlF2M=Hkk42%548z=`s^c>_%Rx}XENJS!_b%e%Fx+jDHhF6_hRXiZ}J'
    '6`XG>)i((Uf|0Cf7jYv`&|Y88D5`8Nx^W?R?ZkfrdHv=4r&_DjT7=ey%th+ae97Y+7NfmOSt!l4RNYbp^g>R2!iuHiM6cKX??sw>'
    '$?n!q;ZPe9U<lvsrTW@x&e$vnqd*43I(0!fR!owYxZOsrLvk!MmIw8c0-+`6QJ&yJjX&{i<-7w7<m;Tqp-a4N<Tv{6{4<e(GGHnr'
    'HLF|47fkF@5oPDYBjXi1H_%f7MxJwPAFqI&LSQoG5FyfF2yk+1+@O(DTszgnSDwE7a%`_0nVo20B?e)5-gM=##)A1?W{VcD_Q%zc'
    'XBG5J%6<~#%>@^cjp2PHUuv*$H`nXDv0OCb?qRr@@y;(=D9-tY48f0vGnA#v2IKF>tWZR<+=AOJhTI06r5vsF;n6<u-<acmgS@Iz'
    'r=L0}1gd{p^^1OeEK);ph#*s$^MUZ^YJcWGE`Y&eVzT_P9M~||&xwlimz53a-xwGVR8|gyA?K3@$<8ypEw3=wTZ2h0hVNhR92piI'
    'ls8?hmZmgqiXf|hh}R|j@d?mv1zh>Z?otQ*XipUEpIZuBP06@b6T+g23HH?%c5;%`K{Ji{UR7=m5fiJX#%Ia$9`9A3?gFn!^eA>u'
    '>?*n<I>Ns(pEbgqZz~4eRvuP9h<BsEUKS2C;oC8y{BJIHiDdh`UIc}ThD3CS+lV=v#?PuCSQ!De8be7ILK$)FH9vR8<?(?Pbyd<B'
    'V{(R3jWdK?eUHvUNkNvNC8WF1K>5oB%KOxK4n4nJAP!pKWwGYdF&7b})LC*W4h3JM=>)ke#a6KXhdW619hoovm<)zV585J{o1S;#'
    'W9&-vEWt4oN|i|UkcRf!($UzTRz&~O*`3lM)ZXgsNSPaU5g<Ujm2K(t%};xR_Un3iS&`!?IZrOzuK*^cgI}$(0GJ?aAW>zh@tm=}'
    'YL`b(M2Kk^oR+pi6+F9pFtqDY=4QXhHm6)S!k8`ap*~+lqO~i-@iEmj4&6b8l4DA1$3Vxz9ArLmae?7pS<}+106i7ash>+JI?Zi#'
    '1+kc56na`-*Nk5PdO1(N4NVdR!I4?U7d7Vu&xR#urO~9HYXXL9>?pbM$SloE8LCDaih8z&JCY(#u;bx-pnr*NR<<cgrCXBDgA+o;'
    'R-4CrYMqK~E31=IyqU$uwI|?ombO^|bwy}mQBA<PcUe)8F>CES_uN|YXs1NBFq?f45lUl<WwE6ep>yN5F{EM^)~>`I-;8q)1Pp__'
    'J4O?Nj+H>|y3&_r)H@;t<j>{^Q%I<+mi<K$UsbR~YGP>E2ysxuIdCR{R*XdBYcSaCsC<#Y)<%!VrfeLRox}p(c#vNyT<NdN88JU~'
    'V}ora7F_H!s8v0p7HDM^FqWiP6HegUaYffZFpiAC!oa!RkWzd_^UqrgrCXuvs}j<(aoOWk|K+$DvQb9c2T1lF9Pz2VO8F`M3m~uw'
    '0}8$WDV7SO6z^_;Wyj5OV&1VOx()Jb;X!vs3~_9ro^#j0`%qF|XL$~nH(Q$H;!K2)pnUUNq~F?t>uc5VCpkw|eq=a@MKfwb=$?F<'
    'O@v*sG*w@d9leeN9p_4wrzxx4DM(_=LjgDd#Bt2AyAalrR|wyb4Y^B1MDE!&2Q8f1MauO8gJXfC#RJe~8t@GJtZo<C`<rPrkWbVd'
    '07W(Flf7)$Pe@-qWxtMHA!z1af&FgZa^XHF;vl75isHC%tJSeJ`SsBzFG${BV}9VT*D5eHK*lc>a<f2e0;<!oSY+6S&%$46h3FR2'
    '(IOJp^PCU={@X#OGjHxf6po}`uGYgMoWMjxP&0nwG7r~pyU~2sv;p92zi*1k5?~Ymt9LQ9$~h<p8<swaB8NOMG;8R24wAuJn&F2j'
    'bGK+vbx)#>5c6BaY1TG8qbLPR1n)7|D&t+pgP1}iJ{H+XSERv{>Lf7I$n^u?yiIfyEnw{oFZbhnozgF?*9Qk|@JiAxdQ$zI(~vLq'
    'S(^N((2nh-q@G}a+!!C1RgI0+`;X71l)Y?zzHq=rx02%-jH$R`?<(uOaec)do#{?GHkh`SJXSA6fmZx+Y3?qX--8iPc;qW9vNzPO'
    'K+6O>Yc;w8A%>x(UKw29LmB62lGXJHjbarycejE+>ZUZRs_SFjqy@T(r~h`#F%`S>ljVQJfjatk-JJZra2sGi<Qj>nFmPJ`QCT86'
    'TuO4uFz^W*d-H<UlEtGl*BspjAo3*|TPIDw1L)GxKys$xg_Jh17}P{r?icR9XOp~l`1^Ipud>nVK+DJ+tynL7<x5@Y&Z&4Cv_a_*'
    'KGYgcL980jd-dgg2wNgA7$Bpkn#u7RgY9Co3wGbNM7a%^*tUY_)j1-!vTcm_)_R{ZtJ`mVmlIGnCCA37siaZRTXRdyh@kQ%ad}Jg'
    '31XB;;Of3-U5(*Qbr5xMQW?LmKxq4?D$L^4acwJ!OvGo95kp<XI1zPCj;wWXA!^zWVvszJ6RWYq#Ra?0wQ{cNaY>vuR5OIiB((7>'
    '<sw7pl<hp1IR#Hsm7b1Ogo?m~nZz#&5H=THH51$dO99dt-~vg8<BLs1zq{X)GU&&|@1{3Gzf=40;ePY^`)$dMPs(nkAF?@Oa5<X5'
    '7fY#Uh#fX%$2eEe9wGQQ>RN!srMDXAMo4W5n};!4&pXWp57V}oYW$8*oaTNM*JEiH1Kv-I$KzE!I5?(z@zhVvnqX5_EJk6pJEZ}j'
    '{W*#R%!3HF^~`8>60L^}4#~g1`ht|PpZ_=uH74jn%|?kQ?(d{0-15)ee7D>6s!H}pm+g*b?G6ypq)N*g+c(a0Wfv<q95e%`v@m;#'
    'u%fPl5LXZky$>&8BGe=Qaa1c9v345jK+MCw_BQ+Hn6-^P1V~2pms(q4w**NSI!}Zz;s7Ut#4e`}iHzQkrq8EcJda%A)Y)WB+upX='
    'vwhD@vX{{~!u%;;+*($(l`^R+EhYiE;hYfdWW==1YKX9L>V*T0aYuM;m}EhX<t3Fr0@Ugo!|kKL_mMG8sNcrwwB)OYf9xC%<!Dkx'
    '_E#cXREnz&5U9$Kd`yUT2h&*>_Oi$^fa7EyBBfT7&ZOWER8++^IDwcBvi?QakIT?0y4~AA01gxL3760MoDfWwwem}~icz?GMM!B8'
    '%2arN_&@`zcF;xTcSB37+QV1bFXYKjbjYvNn&}>~mRN$0X%U?1k<LlfxD4k5tr6Uw>YVmbVp!Dn^5NJOf6^Yu4)Aj6gSriBgzwC;'
    '1i^6&B=h<sCx-$n<~Fziy5PAJhey}f490v*Sy%DR2cRpx6%lX)^R47QM6@{*Y;Wg_N;A4QOgc<~Ai-~0;ZY#0vtO~Kw(4%4IOVzd'
    'y9;90t0$REehD>m1k`pMPD4y0Qwdw?CeLsB@wudKCg!69m3~G-%}aYZ6?@a5mZ8JsDH$AS@nHc7ru&PBKI{i`A!(QL7#`dw-{>E$'
    '5by$088D<ehf5GEQrPAFxtoFOA)=x9M8kbi87$r>N^lwwq<+yw6WZMwi8JLnQ7jx(sUK|2Z}q5sY`$G_0vGTJuW$4W%+EpP%WG#('
    'X^(z~S~pE!H0sggC)UV&Fpe_C!6mF}!4eT`y<NYji}I=AgSRr>wC)^Bz~h>Bh0#p<68|k6$Sa(7xMB@`O~DO*PU%=C`B%LwCyyIC'
    'C8Ekxmo#~@N5N^`TO#cYXXhYSpFMrMy3U4Z_tXf9>1e!9V70wYPuR1_Hx(hFk0TZz-G|@=60S&|FYEAQ?WxWuF!HtyB7=30&Qzyi'
    'am|2B!OYh%PvhJKf1-cehq9G10Jf+=08=x;37)iu#243BCj-wCfqc|#U`azNahjRr6>BFR8*YP1k*Ve8{_BjQYas;JDt@ql`9t9V'
    'ruwOC{a;5&X+1RR+UJ(e=O5C%hTbORupoF(yymaPlPf2}KYqVdoAnQ$YEJ}~K^D!V))5k*VR@1A1Qi)-kCvHMwb#UrxCaBfDI#SQ'
    'jEqQdte`_!2E?G<Q7H^J&bk|^sh=YlIO{Jb>(r(x>Z4AODCv9-+>9umW56Lh*&gWjx3@$}#4VTMaUE{zabfCl4~Rk6=gb<4s&0l*'
    'vEB*DS78ny5R&z|m(BH0A^(sVAm(n&Yb&Jzi|VW&?dZ`O=U-&^zwDojIyeAP;~v6v&55Z=TiBswO%|=~y^UsOLBUrP(iK`)nUc-h'
    'ZN?iF(`lU1QzhNsmU<fy7w1bQwzBTaEE7GRTilID{jsi96vO(RYNrRc>~#(N8uKEj8DO4<_IWJ!iza7p;k|FOtY$!$xaEEP!l5kw'
    'E=(>_RYfKHD6&DKFXr&s*I=6S+<E=kbic57x3fy3f-t<_Ly+?U17;hSjA|)tzsdz9tEa^Glup|=Vg)V$=zG*YkbT4-Gx3nPnq%K9'
    'i^`1C$_*p0&Omd9OkN<buQ1Zh%94ra7;CEMdPz%Di#}A6ualMl+R@@z8;`j99yKS;+O<Ogleu2SuVEhuPb|v7S5VQ@8eNc_t!YQx'
    'CQ3{t4cf&&Ds5U7TpmF3hh=|81?8R|BzPtLV&dPmJD=~MZTZh&!`{ikVFB-cN1;)2%#smfo!V`~yx#J7SKv51_Jg*Y2j4)kb61-1'
    'F;CiLRA13YRK#*wFFPYz!<zi66E3RMB?bW?j-5yyn{Vm|IK?Z-^I#CG1ookXvscEFWx<`viCew_c>c>!L=;8XyT^K*u_Xn93^0ou'
    '4|Clhg2F69U}xz<1T(Xs?pL(n;;^asVm4Mn2)(rQEjkgsHuxC5lT?-WI*(=_vjM7X6?k;i(>slp-nfUgaLKe6C;0fTanJ1Ep*BZX'
    'a~?pl)jja$uI`wCOCcxgrd0u~ye-)o;{U@f=>0JU>!@lqLt$bpw+^HSYp5)8;J)BngF2q8L^-!POA7|T&m3BfcsVqbAP0da6nq2s'
    'G}YYRhc`U<)12WQ(<|k{voXeZc**7846J^ImP05VMc$K=&}1*h!r~lHp(4ktAV^8(#kFF(;$N!>218sZHo(4v%--6Nth(ZoSBE9s'
    'bW|#E*oj6N`()@K<QSNjdz}cD6yvw1K3^B-&*o9oNWE-ESwUjR=v-}f%g-8u8dYY~=8vYRbS2$7eubx|7|aPty#J+#Dt;`|1<j|s'
    '`S&m=`Qp>FZ&yb);Fv6IS7d`u@Fbx&p9UurfrH=HpD7hv0lx4`<a$eemg{5^=Zxp*t^(37ufKnk-)9|3{}Re?9^0j+lVu~v)^84{'
    'Lm*`<but2)-){*)ZAlDBmp{I7=iV+N^q1$|C29EJvk+|aHwhmdN+Jh`FZ|V{GT5*4w{-Tb9+690V6Vr3r|m3e{<s94mxkXXtVuiy'
    '^ho%=HJu2OFm9s{6!{7me1_XSkFEh|17zNQE4+Gt;T#D%%}FB{v7HaH-C?m_-<PRZXL=h=HcZrlt*i1y6>c&TT;K?q6q-tn5kD1_'
    'q4tpI5~p4|rzL&g&xst~%?L73Nc(2Vi{1`j{r@3^er-FMVYduMnd)>4*3L0;ND(T7)okUlrF_!Vz;HXB{y2ZXGcULu&6VfIk$DQH'
    '#}EJr;?f;5=y5sPu{v6Zgg;!M^`A_)%v+f%%q{gSjr0#kN#eac%y|QrEJ`^KTe3t+Fb<7<Xth87CQ@YF*ccet1j5qFM)m`7rZAb7'
    '2y_XXV+Yopxn>UShi$W9l8hHav0Y;of^zza5Ou~+%6E1+)Q1IyX4`jme>Qc)I{hnno+qO9#{ms&2(*RGe2k55L&f3rm2(e4U)^B>'
    ')oV84!(En+gEO&5Uu|A3)RSX>hLj)h%W?SAk2zR~xvq>SXE>@$d#MHsl=7HpxlkT=C5M6L`7FukjHg7V7-v(vj(g}(nyc3BEwk9q'
    'VC{8-C7X8!`<7%C0LtKOGq2xC2ZTQ#OBo+i5`@v+Obc7CezPK!M`Hfh)*gT46|V-)*@Q$^pCn<8Q3VDNbtlmj#5pnE`D7u6GdY09'
    'Vh3Fqta-!2(Ox4aL^({A@~iutxY|09wR7dAlTkK!c&I{=FR_YqI~SDnw~cgeTHAk_r&Y$i9VlLJB8D8*5)PLi`(@kNeL-5Zf{g(J'
    'M;5?SY1G6$PqT$n=wsWaAmHmS=~XIl!_6XgD<th^hmvh3WKhRmAyC*B`ay5B&s3m&JGgGX;jsDGV&%ZOu(X1Jfy`x!-sQ2%`e!7b'
    'wjy^pI*3|{X-<2_fJ<GF<E2E-_Z(9kDK+rT^CT=YU$?V6%r4o#dd{f7!B`!*p2Fa1ptiofG@ijyFXD`awz!1iqdZh?N#Imna9Vd&'
    'l`LB<GvigSF5GDSLpq0tf}XN0pvJZZST$z-S~;{uBR`9}CG-WuWI5ARMhL&pXKSF0QIT)A9Awihkxt(4R9$@Q=Qs+?CN38(2u{wl'
    'mXqqwPs!bw)24eTSNGoRSY%6IxdcX_b6NwGB%?2a2DsawZ#*i4w^NrlpyC~JOnr2+V@CouZvM^glbxDO4Q$hdfb?grVVgcc9V)42'
    'DN<yDZLk!z%dOaW8;0Bl*=zSr0gHp@QI&X^0dSMLT{)keho=5RAOy2mSIuW!44x;5xrasKV=ThPvo+hld6+R(V^3i)+m+bYCxi^3'
    'gX3!Jnm>ey5R;X`{ayx4r30a=0er*Pr`YYed_9?O^xY9>(T~zU^2SUE0p&~a`vd_KG>Otd%0qLwg0Ky~qB3WJ1KR~)!5n~1$3yFM'
    'jwvx6Ef_dN8%}sDh<?DIc+L3!GmN}?U0=(%e!WN)`5ol?wG|olLbe)xp=I`~-d6=?$%JS5st9*z+ZB9agXD8;L)jQ_Yg-y#ebiCI'
    'X_`Wt58d{b<-uPdAq6Fe!|C^H@O0yzzF5NbEuc!dCNar+6eQX)MK1}l;4|wjm1<Z*+q!Zr8M63k{~?e6NA3byiq`dOttSG>3~Ypb'
    'f?n6Sf)b=)MZ%2FSJPZuSX6wsJ@j7xgZEJ1c-w7wM4@;hV_Ya)W1yPU)ZDDC%Ny9+t9-3*%O^zxI%PdES$8{#wftUwK|aLr_Y|(@'
    'n1jAytv4YBpNw=HRsy0>2bniXx}F1Ux})G1t^3fw;f+BIM0(U*q*rVKGt<j8LWQy8E3VcqIPQ_=eT#r~*-AjlUoR(u6%@b07sF9$'
    '4q{e6?e9|(K2*y*P^!Uv8bjrx5ry1`EhLU@xIOCxyUN;_IH6$ig?TA|r^EK^y%J2jaxw+}mbOtvqNfBSd`4E9DsG&y@a&QFCF*zp'
    'dhD=-jEgxQ9fk(ZvphvSlvYx498Fv*|1@UhHJq1=dwsM83EVkz`cO|>{SOYE1lblO8F|aft5f8Ow|F+?SY!FA2lwbYwfzY4NC-sY'
    'gM!pKjVhb?0H6<R+K1#P03sRZ%k-4)c0NO?K%|fL#)?FMrN$KPV8c;zD4Ss(0mWZ*_2;FqWLgcylTb{690e%XmNN?Oy2})y3$)`I'
    'c#~JoAa%$H5wX}4G>7bVT2^;7{FmHU=Aj9&PTC8r@yO}Hti*D;mwpHZhy>d}b4PqDre_z(88j*OtJUHlg%<N&s&)Brf3v`}T07QS'
    'z*+&9Tc83PFPA(@bw8smon|rnDwoclrdA5OL;AXDw$I5HnICfNPD*3P(rVoqo7^#XzprWVAJt>8Y2ERu(O~7)>-$4TpnNR{jH@)F'
    '^_#u^zFBFY-f6){IU+Z5=MfzgIILBOlh|`d?K#a!Z1JFE$fWTMOpjf+kLIfQi^Jn-Xcm6)AlN52ujW3Q@GFfg4yWA4(pe4LW;T<;'
    'f~fg=Y}*5hnkMglF{4cD3SrluSEZ6eqP8T&t(B!ydXARL!kgmwsK^>tVN==B<8wPjg8!y8sWJeZs_i>={I_9mf8zNYZtXSF->6KA'
    '?kkMRp)YF;74r;RDIlyh#0ZgE2r7cm-cP=wEB0McSEbwVX34GS?xqra9)dJKFeKzstqY!l-ZbK8SIT)$cSmPF3WC8c+7k67+0FK3'
    'ln|;TJ)TlZp*)NfKa$J=hLIKTF+E^ol$7<m&vjdk5U$D&b2EcrkU4b#hoP@e?R=$W4*GcD4&A$}xWN4>@tm}VOJe`ZX6Ii!6y4+2'
    'rzSXZ{5^YeZ!!H%Pv=0)ex&Qtvhc=zq`jbxSFs=^t}vor+q3vY0v&+CVzrS?m4%H6+u*$Sbp<<~6A{;I)fJO2y6V}J`g0QZ8Yl2N'
    'SuFS6I2NN)sg#5tW_N}km;6u-)U$Fo$7`S8xmDF<hr5RwkFFVo>NPev`q{ss0^#D9#tdPUm!d>uF9rm+5V=INqyo*YG7Pk*9*jM0'
    'Gye1=k9>dOvhgLTR)B!V?mt6_It{i{D9+=%9rx|MoFbv{+4f!aHW`APxH_dBiZ*5D%b1a;5@Hl=VfMXshN<4L&;>FfFfV=s$W3`e'
    'X3Ce|A@YQ2PE2>W(ccPMx7ux*g80m;OK$XGsv;oF(F@m|eNgp1gX}Ga-ldd81?&hHT5Inae&TzUJLOs=zj>8Z7(J*XV{qK{5?G@c'
    'x}NO$mAKGR*nB{W9A%*JzSZoBr-tJ=ilBX^=f~VQsmm_tL4Bt+bNX-<@#rPu&Rp&!*Mq9F=;()nou3ULwmCr{#zA~6`(r`qPCX(W'
    'pO9#|r_N)QUE^9A0Knq0G9Z-g*wT{3sh3<R7q8S1Iam7H_=B<@LCvTXzGppUZa=zGN#x8!uQlGiBy-JiVZAGIO=g?P1gXWW*ZV$8'
    '=^rfA95y{wh0#4@xW(P;m_={OQ#Lu@lB*c+67m-M(S_J)?O3<QnHPOtM`YMfhdgeol-Af01w(Lm2o{)@;v;T=qQrw<<tmfWU!A%E'
    '+w_VBVpRx{o;Uj&kbpyYZ<glfk#Eh0a)&y~)4h<)vMFH`g*mj2-6Dwm{Y{C<m6Yc~GpeT<5MyLFpJ{6zWSSFBm)48v!F#fRBAi$V'
    '5oWpY0%z!gO{MS^*vlqTr5eAXQ8YlYmcK3}=c0A2gvKKBuxKd@h!X`az6<!ttjw_Wh@Pnp;})iGlDpfn41m9K<AWFPNj+<FeV`jV'
    'E~#{4CKoLq&CnDTL`Zn|V4LtKB*G#pStVWb=Jr%QNb~&I!Xh_GQHbygr(X75pj3H;$a+ayS*JKp`M}LW*>UD3ICSzB@N=<5KAdgn'
    'CdN;Pu-{8i&q55`G%hNER71uzVaFp@5^o4oi*$a`Mn+Z6;+=SX;!tQ}o$@$0yfWFFp}7@ZYnoNi;+#dH%@?j0A70^aWoiPfLti3D'
    '{)!L!OC$Z<(#OQRYI~~`XWS(U0qTsfI;?srb$Yy2vsYJdH~G~{-Td`zZ>+N;&w|m8^yf<XcO=opoQBXv&2lc^k=I3tcd!B(eJ;KO'
    '`3WPzUnSU3IKfa1U!q?()Ge=nU*UKfCdzx5)m`r+)hOU-86#u;MN`JEsq^A*=U;AIm}QezCIg%d9N^{%f!aiqz?j`!pzzPR!e%v$'
    'e7)?pH}s__;6Ddn=lY4JYRiI08tfhZu8_fpWkizu?E=q?*w#r+WcGY76F%|*nh^@k?#-rkwK+t`Mb{Fr>ppbd&P<U6PnbCZM+dlq'
    'hyUGKHfc(rN?<BFI9qKplFOfyhDOY>>mfYm<5Pg&-xGHdBV&5hO+6O#1Lnib|M??#UgK;FZ@{$Gs-Nq|$XVpzj$@O`B!Y4L!bCOm'
    'O3wyyb_`CYrSVE^vXKN%)R=J^6zV0|$yA1k8zb-ys|z)f36arAc9XEO{_&bxHn?t{Pfi8_HpzG0Sb;;2Xt=Y&wrCfTC+clOT>kL)'
    'EBKtA<~(yua3fM2Qi^d{v>Jn|o#JJl#T^!>JXkN(K~ew~0AE!w>>99!@13+Svr_Tq(X6CCXvTbG3KE5iFIC2CfyW@pPl8wsYjw)J'
    '*y6G=%(JH#?L%r$r#ieHeq@$DZ=^8C)4gO&oRbF^hdo$LVP;ws?x;YzB-WXOuwJp0<VYV_lx@L+Y;mm4SzBdd=|~lLsK$X|f#sLh'
    'S=T@cjlk+c={Ud1?89XH<4!e)VU=_|E%zY6=qr6{eu(@nBTCnfQeFVrSRYO~5s1a5P;d^m5dVaXuQQ|u^v!AiA<yKCviwJvp|&uk'
    '1Z;gKPNeVpZXh4bZhiN>i|9*vqx~@2e9t)AcT%e6$=eQ%yu7p%;V5fnZ{;xUkn?$Y`C>LDyZ8X{lsT2Z(B~T+j0alCFgeYHZ}jF}'
    'P}ha;4<bmmFwMyKn^kZr<6?1hFjoJM;ZXw}IdVf)V~dk_bJnQzu7vd7*i!fcis1hxMk)v6*siTNOrSt4jNQM!=R{eg04AKJGd#tC'
    'VL{-?s=}rB;6nrjWamVf`8cWJF^WlPaU#c__80pF)U_tenJn|^o5ZjK_a1ir7&BwZ(#U1q{fNI$3udq8GtOt2rVqZJfffU4tXB%-'
    'THoex(wyC=YhL>ZyE{I4DX_{csQ>L9TF;Whe^8{PawC$D#x#qUI_fqXKP!IEF{Mt|4{d=sAv%0um^Up6l#$S?G~z_suwx_a^m;g?'
    'u|nAV|Ita;8FpjOE{}oCs#%+yP!+3@Pl{O@v7)K39oDvIPPQetk6=ZYn%Q!JLcx*m`d9-HjRdaz><eGH3{D_vY?&y%Mc0E(=d)p='
    'uxOU5-LuylT<c}M!7hg0Ob%kE7drb3bN0c^6n@k3bvh=Y^!)f!RmMo`{Hxh>`2ylnv6+4mF2XInFqPt%0>(!eP`+tG1Lir?Wt}xh'
    'pzy=I7_{(U<!Q3}PQZRXA#9OL$WLCFllbh%F6r7;ycyWqC~T-;w+ANqFEm}0Hh}fa?(wTdPcKv`m;st)-s@?G+cbMB!L&Uu{Y!;)'
    '1=UQ_dh@qV!g_ImU5ppA>+5@engR}j@xY&-GYf!Tr_00R8onTs#~fqATA~>`y<L<WYR_hRUd~zlTNd|<l->6($tsqc2Qh|UvIZBl'
    '6%Wh_MMqgw_kSuADklNP=<=+<CYEaf2Ewpt0C5@FTJBLo>znH&P5h_ju0U#Q@BS0QU**f=XK#g{j4@6S%|H!8@$U9@5(eN6R}d=r'
    '@rU(qbXODO`SX8PHb-f(H()>$_lZK=$l!#agsQ>#G7v+zGgl)h2?~UE?w0^GUuR&@hV(DrCaB?1pv48eZpw{x-WWq6iV=K7DWvTK'
    'U2z^3zfL?y_1W9zpOI`VnW-hT!bt2nN=JghiQvRUvzvvMKT}sW;BheQL3TREYgV>4X{7pK$b%cwA8W8xL8yajk!iaee<R6&f^L0m'
    '8H)d`KLva5?_uUAirgum))DcR+<HG8buux3+ttPvjs2)ch)D@Yysmcr%avGR%`2omzZUEL)V)sIZWIv1cf4l=mePrZ;<9X3m3`i;'
    '2X|DHym2xEK_qUFk<<C-Mv2y2FtsCvd338a$DtEYk%_rnV)%V>K1ee?U;h}m`lWo(R38)eSNerLZB8K@+@D9cI`(5eX!+`rpfb;d'
    '!P*1ta{zmKHlZ2?_1$Z_uJ9SEy^T?WdIe%>&KaI&nZV|wsK;+`;mL&s4>*Y||6ZEgJl!6!%^H+Ld;YAPG92pr-O%O-+fs|NCc}zS'
    '43EgGD!09%%|NCB1s;=+BL-lZeSa}t&_qQJHkp=*?mIe4Q0ub=DQl!|MODa+-z~^&0zY<32LAA{yo9ZZ%D}Uo=p7K6-EAUoQUD3C'
    'G+z|3UQNgJ{G>S{Qjta#YhANnXVP0M)Oo9J(2=H%wayW!Bh-wt66zjVR39^0%nJoIxruNm9b^)Hjzor)lKnx$1#$+Ozb!k4lx=Sy'
    'SLVMYL<SC`r2(_3t^PpHEm~P^aYfdlQi5T<nLERFalK#}Siet&Ng2D7%Pr|lt`Y!|7jY6)Q_pr-{4MfPkaGaIQgV&U{`&Sp<sZ!n'
    'miVPn_6VGs9J~u5NJFWrM-MJyfS(pz{7Yig_eW}_S3giY9`v@wlzrCzVazg}U+&+|P;8_pQC0(EV(ER9wS0Egk;lk&Y1JEgio6`*'
    'P!itsN}D`LhlaIzw9a8n?m#(bu<!qR^3$cUaDi$j<d9C8;t5t=e{YjOwh4CRM-?eWeFJrjZI2~mKnop4iXdeADueHiX$Vl%Hyl)e'
    'ShmgGDc%}#S$EcZGJc4R{qO6iu$^_Yk(iaU-iy3;{PBRf>p{EB=ClYWB6%q1H>wLROd~3<^BN1Y!)vQPTdcxVlA9h@)BC_dLx%e`'
    'wOke}uBHNXJSD@LFB@jAtHgYF9k=Ph)d{xoH_P_szw(k(CBk4C2Ut5AkmVF3@VSlzu3wg0F#U9s$V)o9e{NtlL;eT_lA}qS_J8B;'
    'ypEeRgwSCIL)CxF&gANa<ib?in{ztz5ear=0*<q>zKTk^hJ4ReOUVxgJ81d*mSOE~Z2<(j?gdKtPr!qQAiqeG&RfHF=88l~#}l`~'
    'BKicE8V734$11WbiGe!EgM0WqEktg}&J#WgODnR{;9O?ff_KcsL5?po-AWpIVI;;iyor_asmdq7&;8otPDY+oMPkqX%t8=;c!!!j'
    '?Ve}y6r5ym(3%e-dpR7;DrU4MVZMkN**MSXBD#Tzp!1kQ#Yr_oiDF=J5H;_!7>mR(CETG}2&~@y$+qHzxaA7H>XV1WSNju{FcuBA'
    '{+EpOZXT}9@M5Z60Y-dF0BXI&vW8~8DU%Up{4UG!pT`wqvl{T<zr@@Q^<jpoRCOei>yjj8>Q<^%kwn77Zm>yQ1_SFF-WSx3U;RA3'
    'zVPDqQ~q{EWf&x_!FIHuy0DNAsU#)yEkUL1_<Me5y4DOe^f4Dv74CY9u2+J?8W0K(Y7qpeINc1nph8hANEsC)A3Hqcu1A^&+Q}m8'
    'URo~%m@(_+Xk^Ub&zWleP8pYl9cdIR@zinS?&QUFQH4ERokxp%2Q^H>QRPi4GNeb@$3L7Eh?sNyJB-NtzQ;kH>D%*vKd9tXocAf+'
    '^e*<kaH_XTsBiHrTWMN6d{(u|PC}F9X7N$@iEoMRScv1=!(B8^t7i1pgk2COC6?Vs#)0P<n8CAkOKU*dKci*x<W|2@Smu#^I2yi6'
    '<3Gl0d{N<jxjquDcD^s8Um>ar6Gz0&o$On<IGJoD-f-RLo4G@^MLMnL&9<ZQT})vFoGUAcW~v%0?khPs(WX1+xgov~km-P8D;BOX'
    'w&WT<AqJ}2w}#OvXwW`piRCi>bFxEdiVKk}bw2*7Kumo4K1w@aJ(gu2Mdxp<tI8Bzl<iBVq@8MgQLtL;WgZ<PP1Wc3Zu&;ruV2z`'
    'UrN3x&&7{tRG1CgP932NlsZq@-xa{Dlc{>as@90MOx!$32u*+=s4v}|bNF&Qa3Q`h%4O_1@qpO?!h*MzLsaUR00FqM0_=?&#H!L3'
    'vBYQl0ssI200dcD'
)
