"""Directed four-ray and adjacent-pair gain certificates."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt, comb, prod, lcm, factorial
from itertools import product
import json
from ast import literal_eval
from base64 import b85decode
from zlib import decompress as inflate
def unpack_literal(text):
    return literal_eval(inflate(b85decode(text)).decode("ascii"))
S=2**72
def need(c,m):
    if not c:raise ValueError(m)
def up(p,q):return -((-p)//q)
def iv(x):
    x=F(x);return (x.numerator*S//x.denominator,up(x.numerator*S,x.denominator))
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[1],-x[0]
def sub(x,y):return x[0]-y[1],x[1]-y[0]
def scalar(x,n):
    return (x[0]*n,x[1]*n) if n>=0 else neg(scalar(x,-n))
def divi(x,n):return x[0]//n,up(x[1],n)
def mul(x,y):
    a=[x[0]*y[0],x[0]*y[1],x[1]*y[0],x[1]*y[1]]
    return min(a)//S,up(max(a),S)
def plus(xs):
    l=h=0
    for x,y in xs:l+=x;h+=y
    return l,h
def sqrtf(x):
    need(x>0,'sqrt input');p=isqrt(x.numerator*S*S//x.denominator);return F(p,S),F(p+1,S)
def machin():
    def at(x):
        l=sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(16));return l,l+x**33/33
    x,y=at(F(1,5));t,u=at(F(1,239));return 16*x-4*u,16*y-4*t
PI=machin();mu=(iv(sqrtf(2/PI[1])[0])[0],iv(sqrtf(2/PI[0])[1])[1])
one=iv(1);zero=iv(0);mu2=mul(mu,mu)
raw=[one,mu]
for j in range(2,7):raw.append(scalar(raw[j-2],j-1))
pows=[one]
for k in range(6):pows.append(mul(pows[-1],neg(mu)))
mom=[plus(scalar(mul(raw[j],pows[k-j]),comb(k,j)) for j in range(k+1)) for k in range(7)]
v2=mom[2]
T=F(-291,800);a=F(28,57);t=F(41,57)
@lru_cache(None)
def r(d):return F(1) if d==0 else a*t**(d-1)+t**d/2+(F(d-1,2)*a*a*t**(d-2) if d>1 else 0)
@lru_cache(None)
def h(d):return F(1) if d==0 else a*t**(d-1)
def cov(u,z):
    i,x=u;j,y=z
    if x*y==0:return r(abs(x+y))
    if i==j:return r(abs(x-y))
    return h(abs(x))*r(abs(y))+h(abs(y))*r(abs(x))-h(abs(x))*h(abs(y))
@lru_cache(None)
def precision(coords):
    return inverse([[cov(x,y) for y in coords] for x in coords])
def inverse(C):
    n=len(C)
    if not n:return [],F(1)
    D=lcm(*(f.denominator for row in C for f in row))
    B=[[f.numerator*(D//f.denominator) for f in row] for row in C]
    M=[row+[int(i==j) for j in range(n)] for i,row in enumerate(B)];q=1
    for k in range(n):
        p=M[k][k];need(p>0,'positive covariance pivot')
        for j,R in enumerate(M):
            if k==j:continue
            b=R[k];M[j]=[(p*z-b*y)//q for z,y in zip(R,M[k])];M[j][k]=0
        q=p
    det=M[0][0];J=[r[n:] for r in M]
    need(all(sum(x*y for x,y in zip(B[i],col))==det*(i==j) for j,col in enumerate(zip(*J)) for i in range(n)),'adjugate identity')
    return [[F(D*b,det) for b in rr] for rr in J],F(det,D**n)
def coords(pat):
    return ((0,0),)+tuple((i//2,(j+1)*(-1)**(i%2)) for i,ss in enumerate(pat) for j in range(len(ss)))
def controls():
    def msg(bits):
        q=1
        for b in reversed(bits):q=int(q!=-1) if b=='1' else -int(q!=1)
        return q
    words=[''.join(w) for w in product('01',repeat=6)]
    def prefix(f):
        st={w for w in words if f(msg(w))}
        for j in range(5,-1,-1):
            for w in map(''.join,product('01',repeat=j)):
                if w+'1' in st and w+'0' in st:st.remove(w+'1');st.remove(w+'0');st.add(w)
        need(all(sum(w.startswith(x) for x in st)==int(f(msg(w))) for w in words),'ray cylinder partition')
        return sorted(st)
    m,h,z=[prefix(f) for f in (lambda x:x==-1,lambda x:x<=0,lambda x:x==0)]
    su={}
    for p,q in list(product(m,h))+list(product(z,m)):
        for c in m:
            k=(*sorted((p,q)),c);su[k]=su.get(k,0)+1
    need(all((msg(w)!=-1 or (w[0]=='0' and ('0' in w[1:3]))) and (msg(w)>0 or '0' in w[:2]) for w in words),'overlap cover')
    return su
def lower(pattern,sc,z,root_sign=1):
    coord=coords(pattern);K,det=precision(coord);sign=[root_sign]+[int(b)*2-1 for w in pattern for b in w]
    return lower_signed(K,det,sign,sc,z)
def propose_triangle(K, sign, scales):
    from math import sqrt
    n=len(scales)
    mean=sqrt(2/float(sum(PI)/2))
    variance=1-mean*mean
    precision=[[float(K[i][j])*sign[i]*sign[j] for j in range(n)] for i in range(n)]
    baseline=[scales[i]/10000 for i in range(n)]
    transform=[[baseline[i] if i==j else 0.0 for j in range(n)] for i in range(n)]
    product=[[precision[i][j]*baseline[j] for j in range(n)] for i in range(n)]
    sums=[sum(row) for row in product]
    offset=[0.0]*n
    offset_trace=[float(T)*sum(a*b for a,b in zip(row,sign)) for row in precision]
    for _ in range(32):
        for i,row in enumerate(transform):
            d=precision[i][i]
            for j in range(i+1):
                old=row[j]
                bias=variance*product[i][j]+mean*mean*sums[i]+mean*offset_trace[i]-d*old
                new=2/(sqrt(bias*bias+4*d)+bias) if i==j else max(0,-bias/d)
                change=new-old
                row[j]=new
                for k in range(n):
                    delta=precision[k][i]*change
                    product[k][j]+=delta
                    sums[k]+=delta
            old=offset[i]
            bias=offset_trace[i]+mean*sums[i]-d*old
            offset[i]=max(0,-bias/d)
            for k in range(n):
                offset_trace[k]+=precision[k][i]*(offset[i]-old)
    return [max(0,round(10000*x)) for x in offset], [
        [max(int(i==j),round(10000*x)) for j,x in enumerate(row)]
        for i,row in enumerate(transform)
    ]

def tilted_terms(K,sign,offset,transform):
    n=len(sign)
    need(len(offset)==len(transform)==n and min(offset)>=0,
         'nonnegative affine shift')
    need(all(len(row)==n and row[i]>0 and min(row)>=0 and not any(row[i+1:])
             for i,row in enumerate(transform)), 'positive lower triangular map')
    z=[T+F(a*b,10000) for a,b in zip(offset,sign)]
    columns=[[(i,F(sign[i]*transform[i][j],10000)) for i in range(n) if transform[i][j]]
             for j in range(n)]
    KU=[[sum(K[i][k]*v for k,v in column) for column in columns] for i in range(n)]
    Kz=[sum(a*b for a,b in zip(row,z)) for row in K]
    B=[[iv((int(i==j)-sum(v*KU[k][j] for k,v in column))/2) for j in range(n)]
       for i,column in enumerate(columns)]
    ell=[iv(-sum(Kz[k]*v for k,v in column)) for column in columns]
    return B,ell,iv(-sum(a*b for a,b in zip(Kz,z))/2),prod(F(transform[i][i],10000) for i in range(n))

def lower_signed(K,det,sign,sc,z,transform=None):
    n=len(K);v=[F(i,10000) for i in sc];s=F(z,10000)
    need(len(v)==n and min(v)>0,'witness scales')
    V=[x*j for x,j in zip(v,sign)];row=[sum(k) for k in K];D=sum(row)
    B=[[iv((int(i==j)-V[i]*V[j]*K[i][j])/2) for j in range(n)] for i in range(n)]
    ell=[iv(-T*V[i]*row[i]) for i in range(n)]
    constant=iv(-T*T*D/2)
    jacobian=prod(v)
    if transform is not None:
        B,ell,constant,jacobian=tilted_terms(K,sign,*transform)
    di=[B[i][i] for i in range(n)];rs=[plus(r) for r in B]
    mean=plus([mul(v2,plus(di)),mul(mu2,plus(rs)),mul(mu,plus(ell)),constant])
    bs=[add(mul(scalar(mu,2),rr),l) for rr,l in zip(rs,ell)]
    q=[(neg(mul(d,v2)),b,d) for d,b in zip(di,bs)]
    hz=[plus(mul(x,mom[j+1]) for j,x in enumerate(p)) for p in q]
    hzz=[plus(mul(x,mom[j+2]) for j,x in enumerate(p)) for p in q]
    var=plus(mul(mul(p[i],p[j]),mom[i+j]) for p in q for i in range(3) for j in range(3))
    third=plus(mul(mul(mul(p[i],p[j]),p[k]),mom[i+j+k]) for p in q for i in range(3) for j in range(3) for k in range(3))
    ed=[[scalar(x,2) for x in row] for row in B];v22=mul(v2,v2);v23=mul(v22,v2);cm32=mul(mom[3],mom[3])
    for i in range(n):
        for j in range(i):
            e=ed[i][j];e2=mul(e,e);var=add(var,mul(e2,v22))
            third=plus((third,scalar(mul(e,mul(hz[i],hz[j])),6),scalar(mul(e2,mul(v2,add(hzz[i],hzz[j]))),3),mul(mul(e2,e),cm32),
                      scalar(mul(mul(e,v23),plus(mul(ed[i][k],ed[j][k]) for k in range(j))),6)))
    def estimate(s):
        second = add(var, iv(s*s))
        second_lo = max(F(0), F(second[0], S))
        second_hi = F(second[1], S)
        remainder_lo = second_lo**2 / (6 * (4 + sqrtf(second_hi)[1]))
        # Jensen on the exact integral with |Y| replaced by sqrt(E[Y**2]).
        if second_lo:
            root_units = isqrt(second_lo.numerator * S**2 // second_lo.denominator)
            denominator = factorial(15) * S**15
            numerator = sum(
                (-root_units)**degree * (factorial(15) // factorial(degree)) * S**(15-degree)
                for degree in range(4, 16)
            )
            root = F(root_units, S)
            square = root*root
            strength = F(0)
            if root:
                exp_lower = 1-root+square/2-root*square/6+F(numerator,denominator)
                derivative = F(1,2)- (1-(1+root)*exp_lower)/square
                curvature = root*max(derivative, root/(3*(1+root)))/8
                lam = root**3/(240*curvature)
                anchors = (F(1,4),F(1,3),F(1,2),F(2,3),F(3,4),F(4,5),F(5,6),F(9,10))
                strength = curvature*min(F(1),max(F(0),*(
                    1-4*lam*at**5-(5*lam*at**4-2)**2/4 for at in anchors
                )))
            signed_gain = strength*max(F(0), second_lo-2*root*s+square)
            remainder_lo = max(remainder_lo, F(numerator, denominator)+signed_gain)
        fac = plus([iv(1-s+s*s/2-s**3/6), divi(mul(iv(1-s),var),2),
                    divi(third,6), iv(remainder_lo)])
        zl=F(mean[0],S)+s
        elo=sum(zl**j/factorial(j) for j in range(16))
        if transform is not None and (elo<=0 or s and fac[0]<=0):return 0
        need((not s or fac[0]>0) and elo>0,'positive cubic minorant')
        return iv(jacobian*elo*F(max(fac[0], S) if not s else fac[0],S)/(2**n*sqrtf(det)[1]))[0]

    if transform is None:
        return estimate(s)
    variance=max(0,sum(var)/(2*S))
    skew=sum(third)/(2*S)
    lo,hi=-10.0,10.0
    for _ in range(24):
        at=(lo+hi)/2
        if at**3+3*variance*at<skew:lo=at
        else:hi=at
    preferred=F(round(10000*(lo+hi)/2),10000)
    return max(estimate(t) for t in {s,F(0),*(preferred+F(k,8) for k in range(-2,3))})

def lower_motif(K, det, sign, scales, center):
    order = sorted(range(len(sign)), key=sign.__getitem__)
    precision = [[K[i][j] for j in order] for i in order]
    signs = [sign[i] for i in order]
    widths = [scales[i] for i in order]
    chart = propose_triangle(precision, signs, widths)
    return lower_signed(precision, det, signs, widths, center, chart)


def marked_gram(n, edges):
    adjacent=[[] for _ in range(n)]
    for u,v,tag in edges:
        adjacent[u].append((v,tag))
        adjacent[v].append((u,tag))
    rows=[]
    for root in range(n):
        row=[None]*n
        stack=[(root,-1,-1,F(1),F(0))]
        while stack:
            u,parent,prior,product,total=stack.pop()
            row[u]=F(1) if parent==-1 else product*(a+t/2+total/2)
            for v,tag in adjacent[u]:
                if v==parent:
                    continue
                if prior==-1:
                    stack.append((v,u,tag,product,total))
                else:
                    weight,term=(t,a*a/t) if tag==prior else (a,t)
                    stack.append((v,u,tag,product*weight,total+term))
        need(all(x is not None for x in row),'marked tree')
        rows.append(row)
    return rows

def edge_motif_gain():
    witnesses = unpack_literal('c-oy-Nsb#q4EsvoJOMgDX_wv-<p1A@sP;_k@tRy3)r*UyM5_Gj<Kwi}Y2QwS%dz~sTBjfT_Mx$S-?%LHjht=YTzTKjLibImR<9}XiH7W(QuLZKynp;t%>$b}!uh-(AYc7^P6v=y0~z<`%ANo+fRZ*S#hRDg&zYp3OB4Yvg6szF0`pRr7+3FA3j&N<RyA>_6~Ls%;blNh#@b?@C&-nr5tp={qFt(MI2t$tbzBtUri#5SlrI8YU`?Sl-}RLl6eSOdwv28o`amNuZoV}FpeUs;?$|r6VSH&<vV6oQuu_|YI!eQoIWn8ebHq!XRN$PB+U69rwQ5D*M=U_1gh!K$Aa&lUO!F0dl+WhL^=Te@X;U9-s4dh`t+!!uR6I}u2OkURs|Y?ZVc<DBtcs?f06;BqdtzDos`BC43WS5lv)I!7@DEx>iaw8{gyG+21Q$z8dBf0Fn<5k2Q`qyve}z4RMCC1o@38bbqCRO;Eku9v!-9!GE=c4QlXbArWn*}p5<#_zac359D<?qIv*I5512<a%S9=sd1~Z!L5R;4!!)$9LrcG=NUoTu6dw9z_QVsCbdj-EO!a2U-6H&|vD?|r&bBdWEgbpk0t*XZqk%{bh3>9kV@z_j^``A|m-o=&_dCcGiHoalBZh20wP;?tliDAMUAPsGn0-|uyXUltZ#OTL!N}Yazy(%BJTQ5jl$LukQcue9NSSAg{2d(cU?t1T<g9J?%a2g(`)7ewFH{_!C8`U^WgmHHDD6yIJ=W&X+DTnpd-#Du1zrfE>r7}spBJv%a-Z;CgFiCqXPdx$ON3cvgk@}gtXonQy$Bc@1@aBKOF>QtwgM+S>Ut{+79|K2%Is')
    def incoming(prefix):
        outcomes = set()
        for exterior in (-1, 0, 1):
            trit = exterior
            for pin in reversed(prefix):
                trit = step(trit, int(pin))
            outcomes.add(trit)
        return outcomes

    def step(trit, pin):
        value = trit + 2*pin - 1
        return (value > 0) - (value < 0)

    def line_bits(left, right):
        return int(left + step(right, 1) > 0) + int(right + step(left, 1) > 0)

    def matching_bit(main, perpendicular):
        return int((main == 2) + (perpendicular == 2) >
                   (main == 0) + (perpendicular == 0))

    for i, (prefixes, _, _) in enumerate(witnesses):
        need(len(prefixes) == 6 and
             all(set(prefix) <= {'0', '1'} for prefix in prefixes),
             'singleton prefix shape')
        for earlier, _, _ in witnesses[:i]:
            need(any(not (u.startswith(v) or v.startswith(u))
                     for u, v in zip(prefixes, earlier)),
                 'disjoint singleton cylinders')
        for x,y,l1,l2,r1,r2 in product(*(incoming(p) for p in prefixes)):
            left = line_bits(x, step(y, 1))
            right = line_bits(step(x, 1), y)
            left_perp, right_perp = line_bits(l1, l2), line_bits(r1, r2)
            need(step(x, 1) + step(y, 1) > 0 and left == 1 and right == 2 and
                 matching_bit(left, left_perp) == matching_bit(right, right_perp) == 0,
                 'singleton positive edge and degrees 1, 2')

    def lower(pattern, scales, center):
        x, y, l1, l2, r1, r2 = pattern
        word = x[::-1] + '11' + y
        edges = [(i,i+1,0) for i in range(len(word)-1)]
        signs = [2*int(bit)-1 for bit in word]
        for endpoint,perp in enumerate(((l1,l2),(r1,r2)),len(x)):
            for prefix in perp:
                previous=endpoint
                for bit in prefix:
                    k=len(signs)
                    signs.append(2*int(bit)-1)
                    edges.append((previous,k,endpoint-len(x)+1))
                    previous=k
        K,det=inverse(marked_gram(len(signs),edges))
        return lower_motif(K,det,signs,scales,center)
    oriented=sum(lower(*witness) for witness in witnesses)
    return {'oriented_lower': oriented, 'gain_lower': 2*oriented, 'scale': S, 'cases':len(witnesses)}
def negative_exp_upper(x):
    need(x >= 0, "nonnegative exponential argument")
    p,q=x.numerator,x.denominator
    den_power=1
    polynomial=1
    coefficient=1
    for k in range(12,0,-1):
        den_power*=q
        coefficient*=k
        polynomial=polynomial*p+coefficient*den_power
    return F(coefficient*den_power,polynomial)

def radial_overlap_upper(C, gamma):
    variance = C - T*T
    need(variance>0 and 0<gamma<1,"radial variance and tilt")
    limit = F(5)
    panels = 192
    total = 0
    for i in range(panels):
        left = limit * i / panels
        right = limit * (i + 1) / panels
        upper = negative_exp_upper(left + left*left/(2*variance))/(1-gamma*right/8)
        total += iv(upper * (right**8 - left**8) / factorial(8))[1]
    tail_start = limit*limit / (2*variance)
    tail_poly = 1 + tail_start + tail_start**2/2 + tail_start**3/6
    tail = negative_exp_upper((1-gamma)*limit+tail_start) * variance**4 * tail_poly / 105
    return negative_exp_upper(C/2) * (F(total, S) + tail)

def overlaps():
    tot=0
    for a,b,c,d in product((2,3),(1,2),(2,3),(2,3)):
        pts=((0,0),(0,1),(0,a),(0,-b),(1,1),(1,c),(1,-1),(1,-d))
        K,det=precision(pts);v=[T*(1 if i==0 else -1)*sum(k) for i,k in enumerate(K)]
        C=T*T*sum(map(sum,K));need(C>0 and min(v)>0,'positive shifted precision')
        v[0]-=T
        val = radial_overlap_upper(C,-T/v[0])/((2*PI[0])**4*sqrtf(det)[0]*prod(v))
        tot+=2*iv(val)[1]
    return tot
ROWS=unpack_literal('c-n=SOOE8Y4n?=(S&tHsNQ#p7(rf7b5Bnb4j<B;{<Dr14$_$$#c|WA;ufP9k?a-fg^rKtXM;mc?L>#eKJ8bUIhtDJCh<duWdftu^rg!E#+;bnJwROaZcKC{QjIrlod+5cTrWX~vT?hU2ZgKS96xvr;Abrf({qUhL1}E5xal{(q&al_3#P3~EJ^Rp;zIDhPhgmbFt&Af^KeDaOncw2FZM~?nEXO8eB^%D~o9nR2$Ku4bzR8|v{w$|C%FOijF+B5_PPVhF3#Q~vFY9<q1ywfLn|*s#<tR;+9Hz`ixH6yn2zNhx^J7kWQJF2ZKGeKaoqzxHiEsuP1ZkaPQm?k;Reh6`ubxo8&Tj(Dg(tnd>*=5NuT&Ik!YBGvJvwZ<3b5IQy?0up{+v;oekW|_RnI4tCGg#7x-~tpCcUGJYP1_M##DQ|VD#q7vw1H_-b5IRnaX5r)icKec!9WwAGxDLMSdjm9QxPH2e2nTHmAxV$W-s$8z=}P^@RLYiFpH_Dohj+Jn^N~86+nBPy;<x_2^`>w>|*Mms;sv7<to2qzIn<qS%9Eqev$HS(-9QYU=8Ox9xgM49cFW1l(+^5#oHm%6U*l@;H<L!V^9&6~PEKq^T;AOw4H@Kh-3Qs#!q9Fi3ND`NsQVaEjmE#S9xJsnW4&g%)}>WgLQ#1Bw*-oevOfMZ*mZFW|2fBgI1zUwz)lM{jg$Bx@;Cv_}8K647)go{A6~^{v?_-iJ*;4zb1VSXSZ>0ww`zbwY86;!;d;KbsPCQ_`cacuMvvd74TV@G4n)8Lpf}a^*Q)#^I=5R-hl#II00{1KJI=>nY!}WO&sB#g4?JHUft1TWOD^o31Kk3Q1yc1-NpwE8rvADQ3JZgtgpO;kV!$DO<5YJ98;?U#JI1p!)39=V2J%_`v)XwhF|>!LW)H^(LlZxga)j>boLBN8GO`5|KtFPk+ku>S9mmC<z!|r=%2n>35=|`n2gr8u`e&%*qZD{~C~Ke<OZ##DjaOAbET0TSY|UN;$W5W0^&mjYAVhtdUxHq>oE-Uk_2NxIJW*Cqo_YryYCKL$SPkikT+JZ3cSJ{*>=!B?BuFz4Qzi*rH=Vxq#A&i^OjV<7hlunF_2gkL+(-A(jFP><OL*X+kK<TVdQbMhB!pGAKx-*)JCMukD^Z{84W^50oJ~#wTH+2s}8+Y&Yz`0ySuZ1rE4r8=nY(_x)z%Wjze$A>Rv=XGc@T16MWpEh?e8URt3adPJOr$vB;3Xrp-Q@tj<@f;z4-#Bf0&H)x&UFku20BUCn5E~Mi0Pbf>%=yPMX<skT)Z|YxE9kkku2kal1+O&~Oq=9Bh9UDI#9J2K~r)(vviA1+w{&`1V6}?8VpT;(&2}PH&xuWns7R@wvJuIzMv&uN);{&fbdN)QYmz&m7gFCZ_WYk(RPj3p2FJCnVwPJH<k5BCW>C2NaMVIP32sRkaiqX|%wP8IFjH!q^nDa!`7~9TS<>sW+`%8FdQ}F)N9BR^KnRcy@j4c+cQ+aUmT-6M3jM}{hjNSHjh}(6PVAtjW#m703#SY}F;UJ4dpS0L^Z9=uGZ0wi)NpnYnL_ckRll~ZxIC1+AY)*Jl11mZ-Yc+x>Y1Br=F5%!>BO5c(@O~TRopkNAHMQA0a<a>S&c?)cH3FmUtknh>t08*!i#GlDR`JW)l`_3H67zGlQ1m2v$2Ku}Zx=H<%VDEAX?a?f)=zSeU#54tRI=c%>UUs!H7LNasQ9%##wob?3_k~S?N<-KdNKQUEJiV=UMryiqBB_JPDMY?g%|?~uxH<+7{#y{j{(@M8zXpwDxMp}t)m|<Fv*=`teP40qK0sbOYW08jCnKQQrop4K!h6Gm?WFU%W57Bx|T)53^<aBt4UQoXT@yAA$DP7pPuEp(0#%oXKK`Ff=DR-B<Z_lmsyt0D_-Cs0l~+FT?sIGZrRL>2}C|`C$8PWtC6UjO8?6?Z)@W(uRbUD8I?4Rn^JAPBijcQe;H?W4W?Ns=a<Pn&Em2Py>JPnSCVYUI8O9@*#(PGt$a0lR}Hhd{*OwY3yZEkG(XSf#MZ9eA+Zj=K-pSr+0t<51rwdbvB~QaRx3B8>3>+%eUhw6B&cgg?FuZytO>iL_l}>onrp1AIz%MY=}xf|sg?(@xV*>xMYfOhlhycCKCqVjS+_iDw!n}gTOunt0nPGrw+OBw*Of!*$t&zng<f=%mFK(^rl;3=NC4W&tjc-J8p!jrr<6pwSCsicfA}K}&Ht^T2UfoGhKh2K?78p;zF7aAuxy^|T1`f0=f0L#+aJ8}X=YYkk@|B|y#Nc$;l!)%C#z;%i7{THJr}k1nJoAHaeCY-yglzX5kCW=>cez)Rlr0Cm94>@T>s=orDeu!IY-^7PAVo=0afY>SRE;<{GqU4$MPQ&ojd&LJ~*y^DsQq^S8KBN_8%%xe|_3VBk!aC`WHFfcoq')
PACKED_EXTENSIONS='c-n<jZEou@47*DHbp!H4-?q<E^#0d)BD0~*n(d*0YbdrzQIzi=_nXLV-)<`Saog)--)!IR>aP3d5p&;CJ@+kUcErZ$&z|0Kw;a7ndVj54`&OocNJZb0VGsPt=e}Kn=?h$}PU7W#JzWb4y;qw~IHt8K!pRxVf>fJp$J2UWnsS2HzbXx+wTgu#xE`@OZOQCNj5$C<iZq-wP`f~xzN(}Fk`^FuHs7lj)vt>6RJEZ-9}qAdk{T=ReaDAKAL@LGGpNoL5%&QCE~Md$Gz>IkT?sCi&=C+Vr-HLXI%5k~b-g)Ik$DnYAiic1+6lehpyv<4^(p$B;D+Qm!M|uHbr{JPWn4_oZhmzc;zt;0tbkU-sKzLmbLel--?G05ZV3+yYA@_t(gRuveRMVnJqeBTfivka-f$(oT6DNv^!iP21+vUzZa`l&&in2h@F&$(A)|RTdAum?hvl#Q^-Gs^a;#8AvQiruu^yHINg;cQ!g#vCX$sq6>@^vmZRQs?)4wOX@EF*r<H3Nm2@0hUG7|l_2{}C`DD(@%Ce+hhG-J|?2E#@s;#%W>raKYgV=?4V1bS<+!g+Yw2<ty=XPEh*Cvw8~J)-)3J~9p>cCE;DQwVM&<XI0yE?H=b8FUf6X*tjo3(iRKC8X#3d7wo|<a!Jr{{Wh$7%K'
def gain_extensions():
    def message(bits):
        x=1
        for b in reversed(bits):x=int(x!=-1) if b=='1' else -int(x!=1)
        return x
    def regions(prefix):
        if len(prefix)==8:
            return {prefix:(message(prefix[:6]),message(prefix))}
        a,b=regions(prefix+'0'),regions(prefix+'1')
        if len(a)==len(b)==1:
            x,y=next(iter(a.values())),next(iter(b.values()))
            if x==y:return {prefix:x}
        return a|b
    parts=regions('')
    total=0
    seen=set()
    for term,multiplicity,scales,center in unpack_literal(PACKED_EXTENSIONS):
        p,q,r=term
        key=(*sorted((p,q)),r)
        need(key not in seen and multiplicity==1+(p!=q),'fresh ray classes')
        seen.add(key)
        (a,A),(b,B),(c,C)=(parts[x] for x in term)
        need(A+B<0 and C<0 and (a+b>=0 or c>=0),'fresh ray cylinder')
        total+=multiplicity*lower(tuple(term),scales,center)
    from math import sqrt
    gaussian_mean=sqrt(2/float(sum(PI)/2))
    for p,q,r in product(parts,repeat=3):
        if p>q or (p,q,r) in seen:
            continue
        (a,A),(b,B),(c,C)=(parts[x] for x in (p,q,r))
        if not (A+B<0 and C<0 and (a+b>=0 or c>=0)):
            continue
        term=(p,q,r)
        K,det=precision(coords(term))
        sign=[1]+[2*int(bit)-1 for ray in term for bit in ray]
        count=len(sign)
        diagonal=[float(K[i][i]) for i in range(count)]
        off_diagonal=[[gaussian_mean**2*float(K[i][j])*sign[i]*sign[j]
                       if i!=j else 0.0 for j in range(count)] for i in range(count)]
        threshold=[float(T*sum(K[i]))*gaussian_mean*sign[i] for i in range(count)]
        proposal=[1/sqrt(x) for x in diagonal]
        for _ in range(7):
            for i in range(count):
                slope=threshold[i]+sum(a*b for a,b in zip(off_diagonal[i],proposal))
                proposal[i]=2/(sqrt(slope*slope+4*diagonal[i])+slope)
        scales=[max(1,round(10000*x)) for x in proposal]
        total+=(1+(p!=q))*lower_signed(K,det,sign,scales,0)
    return total

def compute():
    budget=controls();seen={};one_side=0
    for line in ROWS.strip().splitlines():
        p,q,c,k,z,*v=line.split();key=(*sorted((p,q)),c);need(key not in seen and int(k)==budget.get(key),'disjoint cylinder multiplicity')
        seen[key]=int(k);one_side+=int(k)*lower((p,q,c),list(map(int,v)),int(z))
    need(seen==budget,'all cylinder classes represented')
    one_side+=gain_extensions();overlap=overlaps();gain=2*one_side-overlap
    need(gain*10**6>1287*S,'shifted ray-gain threshold missed')
    return {'method':'cubic halfnormal minorants with a remainder; depth-six and depth-eight cylinders','threshold':'-291/800','scale':S,'one_ray_lower':one_side,'overlap_upper':overlap,'gain_lower':gain,'cylinders':sum(seen.values()),'classes':len(seen)}
if __name__=='__main__':print(json.dumps(compute(),sort_keys=True))
def unpack_adjacent(payload):
    from lzma import decompress
    data = decompress(b85decode(payload))
    while data:
        h=data[0].bit_length()-1
        header=3+2*h
        need(h in (3,4) and len(data)>=header and all(data[:header]),'repair prefix header')
        words=[bin(x)[3:] for x in data[:header]]
        n=sum(map(len,words))
        start=header+2*n
        end=start+n*(n+1)//2
        need(len(data)>=end+6,'packed repair length')
        yield words[0],words[1:],data[header:header+n],data[header+n:start],data[start:end],data[end:end+6]
        data = data[end+6:]

BIN_CUTS = (F(1,2),F(3,4),F(1),F(5,4),F(3,2),F(2))

@lru_cache(None)
def halfnormal_bins(cut):
    terms = [(-cut*cut/2)**j/factorial(j) for j in range(19)]
    integral = [cut*v/(2*j+1) for j,v in enumerate(terms)]
    decay = iv(sum(terms[:-1]))[0], iv(sum(terms))[1]
    p = mul(mu, (iv(sum(integral[:-1]))[0], iv(sum(integral))[1]))
    q = sub(one,p)
    first = mul(mu,sub(sub(one,decay),p))
    second = neg(mul(mu,mul(iv(cut),decay)))
    need(p[0]>0 and q[0]>0,'halfnormal partition')
    reciprocals = [(S*S//w[1],up(S*S,w[0])) for w in (p,q)]
    return (p,q), reciprocals, first, second

def conditioned_jensen(columns, KU, totals, z, elow, codes):
    chosen = [divmod(code,len(BIN_CUTS)) for code in codes if code!=255]
    need(len(chosen)<=6 and len({j for j,_ in chosen})==len(chosen)
         and all(j<len(columns) for j,_ in chosen),'conditioning coordinates')
    shifts, moments, probabilities, pairs = [], [], [], []
    for rank,(j,cut_code) in enumerate(chosen):
        weights, reciprocal, first, second = halfnormal_bins(BIN_CUTS[cut_code])
        def sides(x):
            return mul(x,reciprocal[0]),mul(neg(x),reciprocal[1])
        diagonal = sum(v*KU[i][j] for i,v in columns[j])
        cross = sum(v*totals[i] for i,v in columns[j])-diagonal
        linear = sum(a*row[j] for a,row in zip(z,KU))
        affine = sub(iv(-linear),mul(mu,iv(cross)))
        shifts.append(sides(add(mul(iv((1-diagonal)/2),second),mul(affine,first))))
        moments.append(sides(first))
        probabilities.append(weights)
        for earlier,(k,_) in enumerate(chosen[:rank]):
            pairs.append((earlier,rank,iv(-sum(v*KU[i][j] for i,v in columns[k]))))
    result = zero
    for sides_chosen in product((0,1),repeat=len(chosen)):
        probability = one
        exponent = zero
        for rank,b in enumerate(sides_chosen):
            probability = mul(probability,probabilities[rank][b])
            exponent = add(exponent,shifts[rank][b])
        for i,j,c in pairs:
            exponent = add(exponent,mul(c,mul(moments[i][sides_chosen[i]],moments[j][sides_chosen[j]])))
        result = add(result,mul(probability,iv(elow(F(exponent[0],S)))))
    return F(max(S,result[0]),S)

def adjacent_floor():
    L=6
    def elow(x):
        k=0
        while abs(x)>F(1,2):x/=2;k+=1
        y=sum(x**j/factorial(j) for j in range(16));v=y.numerator*S//y.denominator
        need(v>0,'positive exp lower')
        for _ in range(k):v=v*v//S
        return F(v,S)
    @lru_cache(None)
    def source(h):
        edges=[(i,i+1,0) for i in range(h-1)]
        anchors=[0,h-1]+[j for j in range(h) for _ in range(2)]
        for ray,anchor in enumerate(anchors):
            tag=0 if ray<2 else anchor+1
            for node in range(h+L*ray,h+L*(ray+1)):
                edges.append((anchor,node,tag))
                anchor=node
        return marked_gram(h+len(anchors)*L,edges)
    def sg(x):return (x>0)-(x<0)
    def out(x,b):return sg(x+2*b-1)
    @lru_cache(None)
    def bounds(p):
        a,b=-1,1
        for x in p[::-1]:a,b=out(a,int(x)),out(b,int(x))
        return a,b
    def cbit(x,y,b):return int(x+out(y,b)>0),int(y+out(x,b)>0)
    def mbit(b,c,d):return b*int((c==2)+(d==2)>(c==0)+(d==0))
    def literal(core,p):
        N=0
        b=list(map(int,core))
        for x in product(*(range(a,z+1) for a,z in map(bounds,p))):
            left,right=[x[0]],[x[1]]
            for b1,b2 in zip(b,b[::-1]):
                left.append(out(left[-1],b1))
                right.append(out(right[-1],b2))
            e=[cbit(a,c,d) for a,c,d in zip(left,right[-2::-1],b)]
            need(all(u[1]==v[0]==1 for u,v in zip(e,e[1:])),'positive repair path')
            d=[]
            for j,z in enumerate(b):
                c=sum(e[j])
                c2=sum(cbit(x[2+2*j],x[3+2*j],z))
                d.append(c+mbit(z,c,c2))
            if len(b)==3:
                need(max(d[0],d[2])+d[1]<=4,'nonnegative adjacent gains')
            else:
                need(d[1:3]==[3,2] and d[0]<=1 and d[3]<=2,'strict middle repair')
            N+=1
        return N
    def coords(core,p):
        h=len(core)
        return list(range(h))+[h+i*L+j for i,w in enumerate(p) for j in range(len(w))], [2*int(b)-1 for b in core+''.join(p)]
    def conflict(a,b):
        c,p=a;d,q=b
        return c!=d or any(not(u.startswith(v) or v.startswith(u)) for u,v in zip(p,q))
    def lower_shift(sign,mat,shift,T,K,det,codes):
        n=len(K);L0=[[F(x,10000) for x in row] for row in mat];b=[F(x,10000) for x in shift]
        need(len(b)==n and min(b)>=0 and len(L0)==n and all(len(row)==n for row in L0),'shape')
        need(all(L0[i][i]>0 and all(x>=0 for x in L0[i]) and not any(L0[i][i+1:]) for i in range(n)),'triangular nonnegative')
        U=[[sign[i]*x for x in L0[i]] for i in range(n)];z=[T+sign[i]*b[i] for i in range(n)]
        columns=[[(i,row[j]) for i,row in enumerate(U) if row[j]] for j in range(n)]
        KU=[[sum(K[i][k]*v for k,v in column) for column in columns] for i in range(n)]
        tr=sum(v*KU[i][j] for j,column in enumerate(columns) for i,v in column)
        totals=[sum(row) for row in KU]
        q=sum(sum(row)*value for row,value in zip(U,totals))
        lin=sum(a*b for a,b in zip(z,totals))
        D=sum(z[i]*sum(a*b for a,b in zip(row,z)) for i,row in enumerate(K));p0,p1=PI;ml=sqrtf(2/p1)[0];mh=sqrtf(2/p0)[1]
        E=F(n,2)-tr/2+min((tr-q)/p0,(tr-q)/p1)+min(-lin*ml,-lin*mh)-D/2
        factor=conditioned_jensen(columns,KU,totals,z,elow,codes)
        return prod(L0[i][i] for i in range(n))*elow(E)*factor/(2**n*sqrtf(det)[1])
    w=list(unpack_adjacent('{Wp48S^xk9=GL@E0stWa8~^|S5YJf5;51ARi(LQ*flY9qok~<kThWBsq9||~I-QepDA&OvhK$vt;T|GKwOV-}8wX7pxgM$6ya|;buQp?^0wcHi-z&@TkI~=y>vhs9+aMfTame6$Ph{a#V)_~^E>Hra4QtB^430?22i;o5zfI07xQ}LN)h3$A{x2OhGTU<{Br6kTDjS~k6|2m9_|G0$OqPJC?B?~hn-u&w3g}U4_u2XZlp{x6wS*0!sLh`7eY1cczokAr!~JAig_T#Fia&kx4M-acS`Q7m61|xt-ta9drDCeOXW!x>ftLe?hlXm3cT&9YRjc);JV1sj$d|lC&JIxFWX+g8p-QL%{NVn>p>x8BK}+Y;W%~wReZK&7O8!(&Gl(ii0Zq8h4!EX{Y=rIPg6fLwYLjm7&{?dC8Q@X17#+`gXhN)E>?u|y%mZ0%TzAWSE}ExdPfhY!)XW`!=H}1XZ}dthDC8O1&4EzCi<-hAz}^)+^{RnO7gxCa;^9^L)Ayth3Qo)4rhB-Y3QGiMZlwgcd7GiSBBa8KhQL=LNzWs0Ym73xBGLn?N5ad_D&}DO)>sx|G7GwLniIKu3&kV(sTseq4o*eQ4j0n41h0jtj$Gpgbs<ytFhq+2ovg4-#qj<vWu;G<4E%amXPS;_X)0TW)SpeAoq4f>?qz%0=)1tX{0sAgm7UJL4^$a$ZQakPvz6;ar*gF=Kg6hd&6#1yr+mZ2@CQ=4dlzJyP=q+SY6Lcg@_r0Q+g!667%3aJfiJ6DiXo11>(@LF&WHf64p<3LqQOzj1-P-^O=Kccjw;!#=us0#4x&*yeV1@@FQc&)1*Gr1)m#Ha3cx*PY9e&#_-Ms9Yv9Y914xChm$jiQug;g8DwNg%RBV7TwR;&xJ{P6MuVst3z7Zi#xiz=IEn>;LnYzZLvU|#hlhRVgeR3m+fC-3ycnFD3Cf!pvQ>xo7-Zr49p4xC4oF>Qzu~b?Y0;+&oryTS|Hlxg>eyt~T3^M?hC98~OLZp&?H#?3)Xzr0HT}lC<$u%(UlzIh8vf^ONp-wEYk@KmZkM70a`E&UzNV}0W8AYQ)31`5;E@6287da%{9|r$GcphaFKjj(uAV2D#AF}>}raHy%G!{Y2aljMrql@30=uK$NtKddZ=kNxM#5OUOjKin?Q0vACmJk(x`P__raN>ER!Q;>BI-Uqis<ucBs^--BD!?_6UqTwc-HE>sI!%ybLaKqASeynuT{a33D2|;=^)f*+V0{}U=i1qS%h6*NmnT98CU1v2^^{<)#FVTM8`mW~Q||h!FXT3OS{;jNiCz)Dgc?le>f2)60fb;*Lv)oD{@#U9`FYA&yvPM;CO;q)LopiogAdyJM?3!$4VCep6Z>p8?$@```BMib4=<*$!quHn5%K{mt#c~<R3AE3S_{z_<0QBdVA7i0A-@kQpaO(>!y3=7L*Ib75d)gsC+!=097<EN;c$*?_$SIihF+vy1t4j*-5Lu;%(?g%4UtTrNEvto7KAzAnBvG40u8?w9fsQs!hd@+HJ~2hJ#3&fv_)Q)aH~Ac))CuwR~DjdWSrGEzcuiVF_^f}g`<$2^_i@Ym<rH2dRN}~GAwIGbwV)};xVzL(pYfJ_<8RK5R)hfUncO7ahvb!Bh+e#PK`(RBx(W0K?Irja^0vC`{s6j;84q<G>LTLiG>d;&7lSzTBr~mBF?KoTCQkjCLj_@g^_$1FyU=;1<1^t_Q^9<?GLpsl(%6nLpuOX2znO&;?y(394B{(e?s$7&Z*kf!h=_QFll4InfynI3Q!7W`WOQ#!8a?Ff$IBk&N#dVi9A#&O34CU0$4%ID$f@tsM4iP;%I)kIs$MN0KBRRxPIuJ9OGE-I-b(2nb7{x@@wt2$bi)}Q<Xz>0UYrSUEaW#mxl7vgLkm^cvqcqYG3&jn-<Oi4AV{3SIKEflJRXc9n3jNF4>xy3Y3nan_qO{LoHCh&ydcXXUHV4CYt+mbII_XFZ5L9&5ksoz&5M0rb5UDa?k?VpoA3^f*%8!hE}MtK#{y)P(4d$H^Lo#*gVdMlk9!wRn(*}7aQ&+@y#MW(8!()JdX8juIuhc?9&5tZ<G#s(?c+x?acP0r$Pi~Z`45W@_O9%P+x`Bvn;eMPidcZuosfM5p$6mbFyRyogwBu45?vK8*^#OlZn?R{9HfTj$LA)CyB_B?YbLotB^KDT06h5{motNMY9+PlHP`bG1rqozhqbyI>fuXvZtw4Mw81%J>@dsq`#euev7suNJi^i^2uB^<i|FHS7cD?@2n_8JUnS2Bt}-YZIY?m$tG4pI(^07po@>z-uJ|DqkR{Cvq~W;%~J2cq8QH8wM1d*BaBLGSq|3Pdu3C;PU0ExTs;xhHT}JwtQr?G{+hVEfvv%1oa$y<Ykijiu`QV)xEWkSz<1AeojtBdW?t_U*bNSAC^ZY0PY!g@UDNcy3VANfw(13QrAwUw5zdHxC8|K^+XN)yrK<7i5#ewa6q$IC=o)(b#>tG;tk|+u{wt#?*5If&8}q&p*F91<*I5{c{<Igkn1@L$#DiP_u_lN<5)Lrsk8=N`#kHat+99)`W}T8hpYnxVKkSIeU=JhWWp{yvSi`PFh}BuOoffi{bar-)8J|}zyPbYK1Nq%tLDVA-?}mywgI8<DgEri<W$l40`1}~<hx_QmBF8jXFVc1EfP~gB$=fK*y84qCiApZe{eqtfr(zhbqz|9!(G(9G8A%{iF0`sR6ikGS1j86Vt3cG9k4Pc^<?_-l<CI?agiz+3Cd7$WySo{0z{P7dGV-;0C4V_xqo~s;(eNZ>M70xO1(mz2qyp{4QSGKm`)xH(jwf0{C;*{C7Q!{%?E*LEm*a0)R?!^-dr?^cfFrRS+(aSJt}8KWCDDi&Fa3@k3Dfr8mDiG*_G!6jwR;ff{|=K+jd5Ox_51NXk{}zerG^^x$UkpT$XZd(Xmu9BIAdjSluIq?nCC8iCSoM+Ipo@Vz6m3vd~hO}WS$)U8{QBq@eWEj=Es8+ouQgM1+ziwGj(*mWY3e{K^Em}H1T3AbiIlyy70;TT|EGx_Z4h(YTM<nQ+rD`#o4taP+W7&55L9X$kbOXtsH=aeEhf+ByC>yjAM`SDON@Sl4MXYuT&LeBE!Ode_0j^oT(k7Ls=}XKBKCbV)M37o;{l-b)aFJYrrd3p9>r56_Y-PDMpYI$w<zR(#l?a(xIJw{U@aK2;+)=95R3)lf(X5(qCV@V(J{K3xo4m_1D$_Sx<`_Q`M(uzh=`IOC~Q@ag?m6;w7W`Jz)^E2nPuejkEofDt>=v-`54jV(tF$vrBtZ_Xl{PlCm6gXaoC+!&l<&$pj)tp38}wRlALe;hQ`SN_{_606S8DU4hp!MlDobQGjVcP+7U;&n^~~9ug3oHEuMC9-7G1!_A1gyWB$EjQWwbjxS55&PDu(oRNj@DRVnnR^K9vSZ$l3+}sx7%#-3l89n_bAMkmTa>K7*b8PLS=^69tkms~HWK9(-B11iLKF+>G`XxK`R{0H!=zgmVaE#UjL=OD~<9%7u&V})SM&Ej0Zk!GIb`OLjAEfHLpMMz6J|_Qe@jV+`?nfPt8n;AJwestea9DPy&L3_6RzxPf#utUFVBot{9Ch)Plb4>eUSQCCXtuc}(EHxcf9WSMLJ-JzT{dAdB&zH4%MIAsv80Jio=77dnf@bRBZ6%Tp}*03#*mULK;f$I)ttnus}7``tgkIW1GCk^%g%~w0dqRRzrwB==&-WzDShGKMytP#UH&EY2OW}7h#T4d1^EHRjsXW^`mT4*4vyjR$QeB%S9&2112dr0uyDHKG^iRbcQ~!4etZrfD$AtdEbIj%1-{i#n`Z0VBy&i$U#K#%-cpm6PIXmHU@8#5@^R)^)-iXtbLk*gsGIplF20&oBiz?sq0&CW1f7%3VIqJ_Z%CWQ*qXdKBJ`;~B;e?%hF$#%Dn0cRSH=6#u1y7===(>)8;bm(46wp&(6O96n@t{43%g|}7YiXimxjK2VK$&%V?7r_?#W*!IYXST#Orn~oifWI@2$bZohpVm3*Pwrc<kR2&3fVVMg^Zp1g#15g9WAzp}-isH1>;QPMT^6;*YqI#-Fw^HI@2H<Z{}+Rx62!sIMaaS8FwQ!0m0E4l|kfg0H-EeXK4|Vm!S5=#eco1f6yX2<Hs%a*c<l2z67<Ffz4anGewV;Q%9_OB*(z5PeP%i>H<WJtr4FHV#1e5vxfkm*ZH2qnKfcR8iN8tobH%@DIxY<8Dp>ePkTuo|GB(LnQYBCKhbKRm`RvRaZnhAN9lVzK6qLnkxLba4{As%=yD_42Z47*&Wui``;)q$3>w|1J7(p<J=72YdU%)*_gGp60{89%!~I3*6h%Zz#;|299;y4wpJF45Be+;8v%iw(nxxmPP>%&2gIN#v<M*7!PAAbIbeP#U9FPk{c6_BHfbD@@>RAxlz1|z;Xmp+RUvK=+}YzQkh2VWT+zc9E+m8NNfnWpmRLo^HnCbFUW>_|EP(LKd_A$o=h63d>5v`V42T{}<PL#9?A)J8k1g_CRMC_ac_i~9L#c*fZPHE-2Ewctd|cHnT-9-kRVKw^RnbShCQ&E8cniUNu-VlQ9VKXFCR^a>;j^<h*Oy(N?)V-laKO$S^fGfk6b(ON>l#{}2u}A0S%!bbCvt*a`s?+Wb-U59?qRW}xLSGrtTn52hG@#P<_Mv#%qisGv;_4+G>z<p<6XhsG~c@0HtBvu>(iU01M;Oy>Lo{N=g=wl%1@Od7s$7jDx0VFd@m9qXO~&SJhSa#X`ts_)z1M%>J8ad!I2@1y;ld8lw#Z*%ap{4kk~>?w`>VL#^>i7l(V8_Q+Bw8l!Lu;7>RYybzO4=N*sKN4g`nHcBNG)rZyQXFY4pCZfB}@$tYY{a%uzi<|X_kocz+8i4euY?XC8Z^)avhnsx*jLN%e?02gQ%U%KxpAMS`Z4^aBD3q}1@hwBQXQGUY$cuT;Z+{JCv<pdiQS%WLg(&WYITmsh5>TLzayFqy#Pu;5DeeNjl`82w!WXbsY4Tc#a46#AiQnl0}^Be0@DTpUDa3AqE+z!dL0Flp5<tR%<0%u=)l~O)8i2o}z?AZy84Gc~sZhZH#(zA7!C;%&8O}GF613rkw>rcA}00E~T&1e7sn^_;{vBYQl0ssI200dcD'))
    need(len(w)==95,'witness count')
    N=s=0
    for i,(core,p,perm,offset,triangle,codes) in enumerate(w):
        need(len(p)==2+2*len(core) and all(set(x)<=set('01') and len(x)<=L for x in [core]+p),'prefix shape')
        C=source(len(core))
        N+=literal(core,p)
        need(all(conflict((core,p),(row[0],row[1])) for row in w[:i]),'disjointness')
        ii,ss=coords(core,p);n=len(ii);need(sorted(perm)==list(range(n)),'permutation')
        ii=[ii[j] for j in perm];sp=[ss[j] for j in perm];cp=[[C[x][y] for y in ii] for x in ii];m=[[0]*n for _ in range(n)]
        for j in range(n):
            for k in range(j+1):
                m[j][k] = 250*triangle[j*(j+1)//2+k]
            need(m[j][j]>0,'factor diagonal support')
        # The exact lower_shift check validates the rounded proposal.
        from math import pi, sqrt
        K, det = inverse(cp)
        P = [[float(v)*sp[j]*sp[k] for k,v in enumerate(row)] for j,row in enumerate(K)]
        alpha = sqrt(2/pi)
        scale = 1-alpha*alpha
        L1 = [[v/10000 for v in row] for row in m]
        b1 = [v/40 for v in offset]
        mean = [v+alpha*sum(row) for v,row in zip(b1,L1)]
        g = [sum(v*u for v,u in zip(row,mean))+sp[j]*float(T*sum(K[j])) for j,row in enumerate(P)]
        for _ in range(32):
            for j,row in enumerate(P):
                diag = row[j]
                for k in range(j+1):
                    if not m[j][k]:
                        continue
                    v = L1[j][k]
                    t = alpha*g[j]+scale*sum(r*column[k] for r,column in zip(row,L1))
                    if k == j:
                        z = t-diag*v
                        u = 2/(sqrt(z*z+4*diag)+z)
                    else:
                        u = max(0,v-t/diag)
                    L1[j][k] = u
                    for l in range(n):
                        g[l] += alpha*(u-v)*P[l][j]
                u = max(0,b1[j]-g[j]/diag)
                for l in range(n):
                    g[l] += (u-b1[j])*P[l][j]
                b1[j] = u
        m = [[round(10000*v) for v in row] for row in L1]
        b = [round(10000*v) for v in b1]
        v=lower_shift(sp,m,b,T,K,det,codes);s+=(len(core)-2)*(v.numerator*10**15//v.denominator)
    need(N==7766 and s>=35925365600,'adjacent repair certificate')
    return F(s,10**15)


def semantic(part, cap_numerator):
    from bisect import bisect_left, bisect_right
    from collections import defaultdict
    from fractions import Fraction as F
    from functools import lru_cache
    from itertools import combinations, product
    from math import lcm


    S = 1 << 58


    def ceil(x, y):
        return -((-x)//y)


    Q = 1 << 80
    X_DEN = 1 << 51


    def ln_interval(num, den):
        """Bounds for Q*ln(num/den), where 1 <= num/den <= 2."""
        a, b = num-den, num+den
        aa, bb = a*a, b*b
        ax, bx = a, b
        lower = 0
        upper = 0
        for j in range(1,24,2):
            q, r = divmod(2*Q*ax,j*bx)
            lower += q
            upper += q+bool(r)
            ax *= aa
            bx *= bb
        return lower, upper+ceil(2*Q*ax*bb,25*bx*(bb-aa))


    ln2_lo, ln2_hi = ln_interval(2,1)
    inv_ln2_lo = (Q*Q)//ln2_hi
    inv_ln2_hi = ceil(Q*Q,ln2_lo)


    @lru_cache(maxsize=None)
    def phi_units(num):
        if num == 0:
            return 0,0,0,0
        need(0 < num <= X_DEN, 'entropy argument')
        shift = num.bit_length()-1
        ln_lo, ln_hi = ln_interval(num,1<<shift)
        e = (shift-51)*Q
        log_lo = e + (Q*ln_lo)//ln2_hi
        log_hi = e + ceil(Q*ln_hi,ln2_lo)
        return (-S*num*log_hi//(Q*X_DEN),
                ceil(-S*num*log_lo,Q*X_DEN),
                S*(-log_hi-inv_ln2_hi)//Q,
                ceil(S*(-log_lo-inv_ln2_lo),Q))


    need(part in (0,1), 'entropy phase')
    PART_ZERO = {(0,0),(0,2),(1,0),(1,2),(2,0)}
    split_pairs = {(d,s):{} for d in (1,2,4) for s in (0,1)}


    @lru_cache(maxsize=None)
    def binary_units(num):
        a,b,c,d = phi_units(num)
        p,q,r,s = phi_units(X_DEN-num)
        return a+p,b+q,c-s,d-r


    @lru_cache(maxsize=None)
    def entropy_support(masses, paired):
        units = binary_units if paired else phi_units
        lo,hi=min(masses),max(masses)
        if lo==hi:
            lower,upper,*_=units(2*lo)
            size=len(masses)
            return (lower,)*size,(upper,)*size
        exact_values=tuple(units(2*x)[0] for x in masses)
        origin=max(1,min(X_DEN-1,(lo+hi+2*sum(masses)//len(masses))//2))
        _,upper,neg,pos=units(origin)
        tangent=[]
        for x in masses:
            dx=2*x-origin
            tangent.append(upper+ceil(dx*(pos if dx>=0 else neg),1<<51))
        return exact_values,tuple(tangent)


    def clip_curve(entries):
        cuts, intercepts, slopes = [], [0], [sum(weight for _, weight in entries)]
        for cap, weight in sorted(entries):
            cuts.append(cap)
            intercepts.append(intercepts[-1] + weight * cap)
            slopes.append(slopes[-1] - weight)
        return tuple(cuts), tuple(intercepts), tuple(slopes)

    def grouped_clips(terms):
        groups = defaultdict(lambda: ([], []))
        for index, weight, cap in terms:
            groups[index][weight < 0].append((cap, weight))
        return [(index, *map(clip_curve, entries)) for index, entries in groups.items()]

    @lru_cache(maxsize=12000)
    def clip_value(masses, positive, negative):
        whole = 1 << 50
        center_twice = min(masses) + max(masses)
        folded_center = min(center_twice, 2 * whole - center_twice)
        cuts, intercepts, slopes = negative
        j = bisect_left(cuts, (folded_center + 1)//2)
        def entry(mass):
            r = min(mass, whole-mass)
            a = bisect_right(positive[0], r)
            oriented = mass if center_twice <= whole else whole-mass
            return positive[1][a] + positive[2][a]*r + intercepts[j] + slopes[j]*oriented
        return tuple(entry(mass) for mass in masses)

    def prove_simplex(linear, constant, terms, biases=(), symmetric=False, region=-1):
        dim = len(linear)
        grouped = defaultdict(F)
        for key,weight,*cap in biases:
            grouped[min(key,tuple(4-x for x in key)),cap[0] if cap else 32] += weight
        biases = [(key,weight,cap) for (key,cap),weight in grouped.items() if weight]
        common = lcm(constant.denominator, *(x.denominator for x in linear),
                     *(x.denominator for x in terms.values()), *(x.denominator for _,x,_ in biases))
        c = S*int(constant*common)
        line = [int(x*common) for x in linear]
        keys = tuple(dict.fromkeys((*terms, *(k for k,_,_ in biases))))
        indices = {k: dim+1+i for i,k in enumerate(keys)}
        paired = []
        for key in tuple(terms):
            other = tuple(4-x for x in key)
            if key<other and terms[key]*terms.get(other,0)>0:
                u,v = terms[key],terms[other]
                weight = min(u,v) if u>0 else max(u,v)
                terms[key] -= weight
                terms[other] -= weight
                paired.append((indices[key],int(weight*common),True))
        terms = paired+[(indices[k],int(v*common),False) for k,v in terms.items() if v]
        biases = grouped_clips((indices[k],int(v*common),cap<<44) for k,v,cap in biases)
        n = 1 << 48
        vertices = tuple(tuple(n*int(i==j) for i in range(dim))+(n*line[j],)
                         +tuple(int(n*k[j]) for k in keys) for j in range(dim))
        if symmetric:
            mid = tuple((x+y)//2 for x,y in zip(vertices[1],vertices[2]))
            vertices = (vertices[0],mid,vertices[2],vertices[3])
        pairs = tuple(combinations(range(dim),2))
        splits = split_pairs.setdefault((dim,symmetric,region),{})
        pending = [(vertices,int(symmetric),1)]
        if region>=0:
            def midpoint(i,j):
                return tuple((a+b)//2 for a,b in zip(vertices[i],vertices[j],strict=True))
            if region>=5:
                if region==5:
                    base=(midpoint(0,1),midpoint(0,2),midpoint(1,2))
                else:
                    corner=region-6
                    base=(vertices[corner],*(midpoint(corner,j) for j in range(3) if j!=corner))
                pending=[((*base,vertices[3]),1,1)]
            elif region:
                corner=region-1
                part=(vertices[corner],*(midpoint(corner,j) for j in range(4) if j!=corner))
                pending=[(part,1,1)]
            else:
                left,right=midpoint(0,1),midpoint(2,3)
                cycle=[midpoint(i,j) for i,j in ((0,2),(0,3),(1,3),(1,2))]
                pending=[((left,right,cycle[i],cycle[(i+1)%4]),1,4+i) for i in range(4)]
        checked = 0
        while pending:
            vertices, depth, shape = pending.pop()
            checked += 1
            bounds = [c+(v[dim]<<10) for v in vertices]
            if dim==4:
                a0,a1,a2,a3 = vertices
            for index,coefficient,paired in terms:
                masses = ((a0[index],a1[index],a2[index],a3[index]) if dim==4
                          else tuple(v[index] for v in vertices))
                values,tangent = entropy_support(masses,paired)
                if coefficient > 0:
                    bounds = [y+coefficient*x for x,y in zip(values,bounds)]
                else:
                    bounds = [y+coefficient*x for x,y in zip(tangent,bounds)]
            for index, positive, negative in biases:
                masses = tuple(v[index] for v in vertices)
                values = clip_value(masses, positive, negative)
                bounds = [y+(x<<8) for x,y in zip(values,bounds)]
            if min(bounds)>0:
                continue
            if depth>=46 or not pairs:
                raise ValueError('unproved entropy simplex')
            split = splits.get(shape)
            if split is None:
                def distance(pair):
                    p,q = pair
                    d = [vertices[p][k]-vertices[q][k] for k in range(dim)]
                    return sum(x*x for x in d)+(d[1]+d[3])**2+(d[2]+d[3])**2+d[3]**2 if dim==4 else sum(x*x for x in d)
                split = max(pairs,key=distance)
                splits[shape] = split
            a,b = split
            midpoint = tuple((x+y)//2 for x,y in zip(vertices[a],vertices[b]))
            for branch,j in enumerate((a,b)):
                part = list(vertices)
                part[j] = midpoint
                pending.append((part,depth+1,2*shape+branch))
        return checked


    def add(terms, key, weight):
        terms[tuple(key)] += weight


    def subset(mask, dim=4):
        return tuple(4*int(i in mask) for i in range(dim))


    def hb(terms, key, weight):
        add(terms,key,weight)
        add(terms,(4-i for i in key),weight)


    def case(eta, scores, b, o, j, t, old_mask, branch, marked, sign, region):
        beta,_,theta_new,_,kappa,*tail = scores
        xis, gammas, xi_rank0, xi_old = tail[:3],tail[3:6],tail[6],tail[7]
        deep_counts = (tail[8],tail[46])
        deep_biases = (tail[9],tail[47])
        marked_weight,marked_bias,adaptive = scores[53:-1]
        reserve = tail[-1]
        out_horizon = o-2 if o>=4 else int(o==2)
        if o>=4:
            o = 2
        in_horizon = min(out_horizon+1,3) if b==1 and j==1 else 0
        a = (F(0),F(1),F(1),F(1,2))[o]
        kind = (0,1+int(j==1),3)[b]
        cz = 1-eta[o]
        constant = -beta*j-reserve-xi_rank0*(t-int(b==0 and o==3))
        constant -= sum(xi*(b*int(kind==i)-int(o==i)) for i,xi in enumerate(xis,1))
        constant -= xi_old*(old_mask-int(b==0 and o==1))
        if in_horizon>1:
            constant -= deep_counts[in_horizon-2]
        if out_horizon>1:
            constant += deep_counts[out_horizon-2]
        eligible = b==2 and o==0 and t==2
        constant -= marked_weight*(2*eligible-marked)
        terms = defaultdict(F)
        biases = []
        if j == 0:
            cu,cv = eta[kind] if b else F(1),eta[kind] if b==2 else F(1)
            constant += cz+cu+cv
            line = [F(0)]*4
            az,au,av = subset([3]),subset([1,3]),subset([2,3])
            L0,L = (F(1),az) if b==0 else (F(0),au) if b==1 else (F(0),(0,0,0,0))
            wu,wv,wz = 1+int(branch==0)-cu,1+int(branch==0)-cv,-cz-int(branch==1)
            for i in range(4):
                add(terms,subset([i]), -1-int(branch==0)+int(branch==1))
            if region>=0:
                fraction = adaptive/2 if marked else kappa
                for i in range(4):
                    add(terms,subset([i]),-fraction)
                if region>5:
                    corner=region-6
                    for mask,weight in (([corner],1),([i for i in range(3) if i!=corner],1),([0,1,2],-1)):
                        add(terms,subset(mask),fraction*weight)
                elif region==5:
                    constant += fraction
                    line[3] -= fraction
                elif region:
                    hb(terms,subset([region-1]),fraction)
                else:
                    constant += fraction
            if b == 0:
                if region<0 or marked:
                    base_query = kappa-adaptive*marked
                    weight = -base_query*(1-a)/2
                    wu += weight
                    wv += weight
                    for i in (1,2):
                        for mask in ([3],[i],[k for k in range(3) if i!=k]):
                            add(terms,subset(mask),-base_query*a/2)
            else:
                for mu,corner in ((au,1),(av,2))[:b]:
                    for mask,weight in ((az,1),(subset([corner]),1),(mu,-1)):
                        add(terms,mask,-kappa*a*weight/b)
                    biases.append((mu,-gammas[kind-1]))
            for mask,weight in ((au,wu),(av,wv),(az,wz)):
                hb(terms,mask,weight)
            if b==1:
                constant -= theta_new/2
                line[0] += theta_new/2
                line[1] += 3*theta_new/2
                biases.append(((4,0,2,2),theta_new))
                if o==2:
                    chain=theta_new*(1-F(1,2**out_horizon))
                    constant += chain/2
                    line[1] -= chain
                    line[2] += chain/2
                    line[3] += chain/2
                    biases.append(((2,2,4,0),chain))
            if eligible:
                biases.extend(((au,-marked_bias),(av,-marked_bias)))
                line[1] -= adaptive/2
                line[2] -= adaptive/2
                line[3] -= adaptive
            if b==2:
                constant += theta_new
                line[3] -= 3*theta_new/2
                biases.append(((4,0,0,2),theta_new))
                if o==2:
                    chain=theta_new*(1-F(1,2**out_horizon))
                    constant += chain
                    line[3] += chain
        elif j == 1:
            w = eta[kind] if b else F(1)
            az = (2,0)
            L0,L = (F(0),(4,0)) if b else (F(1),az)
            line = [F(0),F(0)]
            constant += cz+w
            h_weight = -w
            if b:
                line[0] += -kappa*a+theta_new/2
                constant -= theta_new/2
                biases.append(((4,0),-gammas[kind-1]))
                if o==2:
                    chain=theta_new*(1-F(1,2**out_horizon))
                    constant += chain/2
                    line[0] -= chain/2
            elif o in (1,2):
                constant -= kappa
                h_weight -= kappa
                add(terms,az,kappa)
                add(terms,(0,4),kappa)
                add(terms,(2,4),-kappa)
            elif o:
                base_query = kappa-adaptive*marked
                fraction=adaptive*marked/2
                constant -= base_query + fraction
                h_weight -= base_query*a/2 + fraction
                for field,weight in ((az,1),((0,4),1),((2,4),-1)):
                    add(terms,field,weight*fraction)
            else:
                h_weight -= kappa
            hb(terms,az,-cz)
            hb(terms,(4,0),h_weight)
        else:
            az = (1,)
            L0,L = F(5,4),(0,)
            line = [F(0)]
            cost = F(5,4) if o in (1,2) else 1+a/2
            constant += cz-kappa*cost+adaptive*marked*(cost-F(5,8))
            hb(terms,az,-cz)
        if marked:
            biases.append((az,marked_bias))
        if o:
            biases.append((az,gammas[o-1]))
        if in_horizon>1:
            biases.append(((4,0),-deep_biases[in_horizon-2]))
        if out_horizon>1:
            biases.append((az,deep_biases[out_horizon-2]))
        for k,cap in enumerate((6,9,12,15,18,21,24,27,30)):
            credits=tail[10+4*k:14+4*k]
            for field in ((au,av)[:b] if b and j==0 else ((4,0),)*bool(b)):
                biases.append((field,-credits[kind-1],cap))
                if in_horizon>1:
                    biases.append((field,-credits[3],cap))
            if o:
                biases.append((az,credits[o-1],cap))
            if out_horizon>1:
                biases.append((az,credits[3],cap))
        coeff = -theta_new+(theta_new/2 if o==2 else 0)
        if out_horizon>1:
            coeff += theta_new*(F(1,2)-F(1,2**out_horizon))
        constant += theta_new*j+coeff*L0
        for i,weight in enumerate(L):
            line[i] += coeff*F(weight,4)
        return prove_simplex(line,constant,terms,biases=biases, symmetric=j==0 and b!=1 and sign!=1 and region<0,region=region)


    def prove_combined(eta,scores):
        need(len(eta)==4 and eta[0]==0 and all(-2<=x<=2 for x in eta), 'signed carrier coefficients')
        need(len(scores)==57 and min(scores[:5])>=0 and scores[-1]>0 and scores[1]==scores[3]==0 and 0<=scores[-2]<=scores[4], 'dual coefficients')
        total = 0
        for b,o,j,t,old_mask in product(range(3),range(6),range(3),range(3),range(2)):
            if not (j==0 or (j==1 and b<2) or (j==2 and b==0)) or (b!=2 and t) or (old_mask and not (b==1 and j==0)):
                continue
            if ((b,o) in PART_ZERO) != (part == 0):
                continue
            for branch in range(2 if j==0 else 1):
                for marked in range(2 if b==0 and o==3 else 1):
                    for sign in (0,):
                        regions=(-1,)
                        if b==0 and j==0:
                            if o==0:
                                regions=range(5)
                            elif o in (1,2,4,5) or marked:
                                regions=range(5,9)
                        for region in regions:
                            total += case(eta,scores,b,o,j,t,old_mask,branch,marked,sign,region)
        return total

    d = 100000000
    eta = [F(x,d) for x in (0, -200000000, 136569570, -8292733)]
    scores = [F(x,d) for x in (6630284, 0, 3270702, 0, 6499716, -287503283, 115098647, -117052289, -117500623, 98888890, 97583661, 844290, 0, -98032084, -300000000, 300000000, -280635214, 300000000, 300000000, 300000000, -300000000, 300000000, 300000000, 300000000, 8626228, -173197660, 300000000, 127966210, 93620091, 300000000, -129498628, 126894232, -65443153, -9138697, 4779753, 2751500, -34215222, 428161, 3321222, 154351170, -8620577, 15417435, -33633551, -93299849, -119674545, 18180322, 28727498, 300000000, -56283254, -96121600, 300000000, -98849760, -298773487, -7594272, 14127800, 6499716, 10)]

    checks = prove_combined(eta,scores)
    saving = sum(scores[i] for i in (0,3,4,56))
    cap = F(cap_numerator,10**20)
    coefficient = 1+1/(2*cap)+saving/3
    return {'checks': checks, 'local_reserve': str(scores[-1]), 'saving': str(saving),
            'semantic_increment': str(saving/3), 'coefficient': str(coefficient)}
