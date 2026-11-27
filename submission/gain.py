"""Regenerate the fixed T=-293/800 depth-six and depth-eight four-ray gain from covariance and optimality-free rational witnesses.
No external libraries, frozen probability tables, or random samples."""
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
T=F(-293,800);a=F(28,57);t=F(41,57)
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
def lower_signed(K,det,sign,sc,z):
    n=len(K);v=[F(i,10000) for i in sc];s=F(z,10000)
    need(len(v)==n and min(v)>0,'witness scales')
    V=[x*j for x,j in zip(v,sign)];row=[sum(k) for k in K];D=sum(row)
    B=[[iv((int(i==j)-V[i]*V[j]*K[i][j])/2) for j in range(n)] for i in range(n)]
    ell=[iv(-T*V[i]*row[i]) for i in range(n)]
    di=[B[i][i] for i in range(n)];rs=[plus(r) for r in B]
    mean=plus([mul(v2,plus(di)),mul(mu2,plus(rs)),mul(mu,plus(ell)),iv(-T*T*D/2)])
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
    ivs=iv(s)
    second = add(var, iv(s*s))
    second_lo = max(F(0), F(second[0], S))
    second_hi = F(second[1], S)
    remainder_lo = second_lo**2 / (6 * (4 + sqrtf(second_hi)[1]))
    fac = plus([iv(1-s+s*s/2-s**3/6), divi(mul(iv(1-s),var),2),
                divi(third,6), iv(remainder_lo)])
    zl=F(mean[0],S)+s
    elo=sum(zl**j/factorial(j) for j in range(16))
    need(fac[0]>0 and elo>0,'positive cubic minorant')
    return iv(prod(v)*elo*F(fac[0],S)/(2**n*sqrtf(det)[1]))[0]
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
        adjacent = [[] for _ in signs]
        for x,y,tag in edges:
            adjacent[x].append((y,tag))
            adjacent[y].append((x,tag))
        def path_covariance(x,y):
            stack = [(x,-1,())]
            while stack:
                node, parent, path = stack.pop()
                if node != y:
                    stack.extend((neighbor,node,path+(tag,))
                                 for neighbor,tag in adjacent[node] if neighbor!=parent)
                    continue
                if not path:
                    return F(1)
                weights = [t if left==right else a for left,right in zip(path,path[1:])]
                return (a+t/2)*prod(weights) + sum(
                    prod(weights[:i]+weights[i+1:])*(a*a if w==t else a*t)/2
                    for i,w in enumerate(weights))
            raise ValueError('motif tree is disconnected')
        count=len(signs)
        K,det=inverse([[path_covariance(i,j) for j in range(count)] for i in range(count)])
        return lower_signed(K,det,signs,scales,center)
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

def radial_overlap_upper(C):
    variance = C - T*T
    need(variance > 0, "radial variance")
    limit = F(5)
    panels = 192
    total = 0
    for i in range(panels):
        left = limit * i / panels
        right = limit * (i + 1) / panels
        upper = negative_exp_upper(left + left*left / (2*variance))
        total += iv(upper * (right**8 - left**8) / factorial(8))[1]
    tail_start = limit*limit / (2*variance)
    gamma = 1 + tail_start + tail_start**2/2 + tail_start**3/6
    tail = negative_exp_upper(limit + tail_start) * variance**4 * gamma / 105
    return negative_exp_upper(C/2) * (F(total, S) + tail)

def overlaps():
    tot=0
    for a,b,c,d in product((2,3),(1,2),(2,3),(2,3)):
        pts=((0,0),(0,1),(0,a),(0,-b),(1,1),(1,c),(1,-1),(1,-d))
        K,det=precision(pts);v=[T*(1 if i==0 else -1)*sum(k) for i,k in enumerate(K)]
        C=T*T*sum(map(sum,K));need(C>0 and min(v)>0,'positive shifted precision')
        val = radial_overlap_upper(C) / ((2*PI[0])**4 * sqrtf(det)[0] * prod(v))
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
    return total

def compute():
    budget=controls();seen={};one_side=0
    for line in ROWS.strip().splitlines():
        p,q,c,k,z,*v=line.split();key=(*sorted((p,q)),c);need(key not in seen and int(k)==budget.get(key),'disjoint cylinder multiplicity')
        seen[key]=int(k);one_side+=int(k)*lower((p,q,c),list(map(int,v)),int(z))
    need(seen==budget,'all cylinder classes represented')
    one_side+=gain_extensions();overlap=overlaps();gain=2*one_side-overlap
    need(gain*10**6>1287*S,'shifted ray-gain threshold missed')
    return {'method':'cubic halfnormal minorants with a remainder; depth-six and depth-eight cylinders','threshold':'-293/800','scale':S,'one_ray_lower':one_side,'overlap_upper':overlap,'gain_lower':gain,'cylinders':sum(seen.values()),'classes':len(seen)}
if __name__=='__main__':print(json.dumps(compute(),sort_keys=True))
def unpack_adjacent(payload, quantum=125):
    from base64 import b85decode
    from lzma import decompress
    data = decompress(b85decode(payload))
    witnesses = []
    while data:
        need(len(data) >= 9 and 0 not in data[:9], 'packed repair prefixes')
        words = [bin(x)[3:] for x in data[:9]]
        n = sum(map(len, words))
        start = 9 + 2 * n
        end = start + n * (n + 1) // 2
        need(end <= len(data), 'packed repair length')
        entries = [(i, j, quantum * value) for i in range(n) for j in range(i + 1)
                   if (value := data[start + i * (i + 1) // 2 + j])]
        witnesses.append((words[0], words[1:], list(data[9:9 + n]),
                          [quantum * x for x in data[9 + n:start]], entries))
        data = data[end:]
    return witnesses

def adjacent_floor():
    from lzma import decompress,FORMAT_RAW,FILTER_LZMA2
    A=F(28,57);t=F(41,57);delta=t-A;L=6
    def lin_add(*pairs):
        r={}
        for c,v in pairs:
            for k,x in v.items():r[k]=r.get(k,F(0))+c*x
        return {k:x for k,x in r.items() if x}
    def elow(x):
        k=0
        while abs(x)>F(1,2):x/=2;k+=1
        y=sum(x**j/factorial(j) for j in range(16));v=y.numerator*S//y.denominator
        need(v>0,'positive exp lower')
        for _ in range(k):v=v*v//S
        return F(v,S)
    def source(L):
        var=[]
        def fresh(q):
            v={len(var):F(1)};var.append(F(q));return v
        R,U,V=fresh(1),fresh(1),fresh(F(1,2));core={0:R};branch={0:(lin_add((A/2,R),(-A/2,U)),lin_add((F(97,114),R),(-F(41,114),U)))}
        def step(D,P):
            W=fresh(F(1,2));X=lin_add((1,P),(-delta,W));Z=lin_add((A,D),(-A*t,W));E=lin_add((t,D),(A,P),(A*A-t*t-A*t,W))
            return X,(Z,E),(lin_add((t,D),(2*A*A,W)),lin_add((t,P),(A,D),(2*A*A,W)))
        def ray(D,P):
            v=[]
            for _ in range(L):X,_,(D,P)=step(D,P);v.append(X)
            return v
        rays=[]
        for s in (-1,1):
            D=lin_add((A/2,R),(A/2,U),(s*A,V));P=lin_add((F(97,114),R),(F(41,114),U),(s*A,V))
            core[s],branch[s],(D,P)=step(D,P);rays.append(ray(D,P))
        for s in (-1,0,1):
            D,P=branch[s];V=fresh(F(1,2))
            for k in (-1,1):rays.append(ray(lin_add((1,D),(k*A,V)),lin_add((1,P),(k*A,V))))
        vec=[core[i] for i in (-1,0,1)]+[x for ray in rays for x in ray]
        C=[[sum(x*y.get(k,0)*var[k] for k,x in v.items()) for y in vec] for v in vec]
        need(all(C[i][i]==1 for i in range(len(C))),'unit variance')
        return C
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
        N=0;b0=list(map(int,core))
        for x in product(*(range(a,z+1) for a,z in map(bounds,p))):
            l=[x[0],out(x[0],b0[0]),out(out(x[0],b0[0]),b0[1])]
            r=[out(out(x[1],b0[2]),b0[1]),out(x[1],b0[2]),x[1]]
            e=[cbit(a,b,d) for a,b,d in zip(l,r,b0)]
            need(e[0][1]==e[1][0]==e[1][1]==e[2][0]==1,'positive adjacent C edges')
            d=[]
            for j in range(3):
                c=sum(e[j]);pp=sum(cbit(x[2+2*j],x[3+2*j],b0[j]));d.append(c+mbit(b0[j],c,pp))
            need(d[0]+d[1]<=4 and d[1]+d[2]<=4,'nonnegative individual gains')
            N+=1
        return N
    def coords(core,p):
        return list(range(3))+[3+i*L+j for i,w in enumerate(p) for j in range(len(w))], [2*int(b)-1 for b in core+''.join(p)]
    def conflict(a,b):
        c,p=a;d,q=b
        return c!=d or any(not(u.startswith(v) or v.startswith(u)) for u,v in zip(p,q))
    def lower_shift(C,sign,mat,shift,T):
        K,det=inverse(C);n=len(C);L0=[[F(x,10000) for x in row] for row in mat];b=[F(x,10000) for x in shift]
        need(len(b)==n and min(b)>=0 and len(L0)==n and all(len(row)==n for row in L0),'shape')
        need(all(L0[i][i]>0 and all(x>=0 for x in L0[i]) and not any(L0[i][i+1:]) for i in range(n)),'triangular nonnegative')
        U=[[sign[i]*x for x in L0[i]] for i in range(n)];z=[T+sign[i]*b[i] for i in range(n)]
        KU=[[sum(K[i][k]*U[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        tr=sum(U[i][j]*KU[i][j] for i in range(n) for j in range(n));q=sum(sum(U[i])*sum(KU[i]) for i in range(n))
        lin=sum(z[i]*sum(KU[i]) for i in range(n));D=sum(z[i]*K[i][j]*z[j] for i in range(n) for j in range(n));p0,p1=PI;ml=sqrtf(2/p1)[0];mh=sqrtf(2/p0)[1]
        E=F(n,2)-tr/2+min((tr-q)/p0,(tr-q)/p1)+min(-lin*ml,-lin*mh)-D/2
        return prod(L0[i][i] for i in range(n))*elow(E)/(2**n*sqrtf(det)[1])
    w=json.loads(decompress(b85decode(';4CEwL<9gW;ea56AqHmbjlSxff?8IH`BNVvh%(6g%AkGu)tXlipq00JH>H+wC-M6F0vBD#NyOr*v9QB?;O?PbjRAD&=2@aXKwN4e@)67TFeg0qQ?6O5z50P4Qffl?otUi1rp|clt&?#KF7sqFu5nfe2Z`5#4jLxnahk;aiL6rhsU*qDh$GLZhI2EtMm#E$VOe_A-}TDsiFS||iM+tnrjxb6W?5?Pa+t2TBj8B_NXg@Vm7U<V1+~W!Q^j2z0La^wf;-VI!$x;E!Yj?JkKf_T%Cw#>DZtL9%gX#qzRD6))2XJ4@u^dn<`qQuZ8aSB8W`3Hj(n9C4T5SPcc;^HpTn#L*g+FqIBAL&>x?p;^ztnSS}#vAU`2S^LB881^@0iS7<2*Q(7>gX!=8W?-C0i{tu_?(UF3Q_DKd=G0q;rZE3u5%aXKb5q2&Nq12FB<Ql?avlC88g7cNi+<Tj9fG!+3BY<g$E+&pXquH6(9cB&iEn<6bX2R1U`CL-C{dvK)UmUd`hQp&J4VW*k1l#N;!=zi6>$eyANgAyyvM;>qKt&A{*y(g=Mp+A*~K7PKL%<oSzD?3+%UI8*z>92~ko2<c9ORch;c4BmHCo!oJ_ARdeY;j;oAj>wY?<Q$5EHZ>RNT+Ib@i#}Ina5%GlQf*>G6}Eff?&adteMp8y#Ttt>BASbOYqX;{wGW54?TPp3aTDCsy=FbC}^Op<&^UKf@-fk6WKyFS9#p<XOtyT2vf+z_#hNdK@5Z9lN?c^M$R22#3gCAE%Q}yl{_Sh;5m44k6N)}*fF+`iiBD^@WT)gUW%wu%@E(dd0s1_r0}0HKr79QeM_P?`J}3blyd*?^&fKDaUZU9^LE^64OCO}j?0)U>j<F`Lo#5NA>@Eg3dsW_(sCOTVbeIL*t)**qP>^8>c0fz^OanG%BR67_^~)WVT()kG<x^rKjDRXC4TsvZ7i(pI_pU=ykowB{s?;DW?02la4FXD{Qf?crpdDr$VSQHH!Fz9mU+ycDOWjXs`r1qHvF8v7k}6gz}|C)o^C=@EpEF+#D84GKuB-xSbBc*j8pIn30@x1S0={x$Jg~b^QSzqa0rKHvTWZA7Be-XB~7o4M94&4)+Dj?n}@89UO!XtO-Fu??OU~*0?6qmq)}vtrUd0-v>i44r3o$3SuVkKV_T@2Lh=(<1hxqU&c#IrWF|gA?e@sQqwV&Aa!L=V;*J*8aNb6_LEZk7veZgE_@v0)tK{*_CWkwVeByJ%54&SjXs&4ry8Wgu%*Dt+Y!%jGgzRRVK8Q@ndFqkw?&!S6{Kgdrcj)!kG!O=bQ&VDm8mFY*`}J3kKoyRMe6b5?l%Y=Gu#KS3BFYuI4+$BBNhCy$V?!@vH~3K$fw9sQn69j-CV4R>W}Zmb6UzJIRI_UsZDTej_OyPuj*mlq<jH6x<1uBm6f-y9MmMR#v;Y&Utc65=kR~k>FGPL~W&-CZf_+6dzRKowm~-nCxJOif;3;A>!FaWU;mrhtuCEkOMa;y`CGv74uP6&~0+yR!HyfxcmJ`XUEFN~i2cVh(_aX1XkqBcYN>ve3A?x^1NN*#v!te~TacO&_qQmxOb@8392xfWyH1}W*rey@_ZGr;(b(wA0bE6wd3sMlfoz{)7?}d<TEy9-{*rv6?MW4Q-nZZ!>v^F+;JEa--agOfL;~hXBsuQguBp0+Ha!A~iivx|W4>O5>ZrB4u+NqPn?2PJL4&y3?1^z{umAUqsGMQJP7bZquU40|c1Tl}E$*jTl%`fFQ(MmmyrrN_mBfjS57w$5u7xK@@y@P19+V|T#XP$)~KulSUVE7Ryuf<q3l_xoA*P{Lk*1zh_EG97P0A@BE#Bf052C?e?SLlJ8h9lw(g0B{$FkO^=G;=dy=Sf9krVp$2*G1!as3`>P+~at3!wljNF_h-wb}_1vAf<4LzF>xT+o&+_F(h&HhQ<2Gi|k!tD=!i!<=Wis!3=A(-8Y4-HlAEez;?NF#wpy=iZ<8!wubSmk`@H{1XKS7p=+N{vPaV(RAtSkzPb3wNov0dTbC_8Rb>HJnfL*FL(eEZ+9rE9Q(n<tHvfA<<FyU$Rhj0m6gSLwmSn#CPKSPkuMTgjhegvSy_+|g1IDYNp$Kxv8z3T32<L#KunL5IFHM|Zr|DXHSo#CnQ*lqW7wO2qJvqgGv!vGYxs@n()iMpl2oASX=ZXL7W;wdklTrTV+#6HyLdc+w+52unK8zmETcYb&JRqb>5hBaG0_i7x!V=;ARNiQCiZ}ymn-cYGGbhN4XdbVUm5!Ga&_jU_7Bko+bdB$*IdNUvich2-U`Sp+uPSDYX>J|GrsT4i{0MJjH-jbd)+IlirLCa;r)WW7zB-}>LS~rU7)JEH1zA;G)}(%k0Iw6nNhRUIllcJV8FLPq@=B44^yEqemxa|RQ(t)~5kIKlv^;_8k^MYv;Gzc>tKdhLNuyhx{ka8%MpKFX*z^*vyeWTNa2VAFKxQzk$c;i8ki?uB9U%$Ctq53BC<~cOWl9}8>`uz35;Z2iB_g02v$aykO_Q(K8yC}_<M;c~;$4HWW+T>&IzZg1-*1kq_i;)%u^MyDe|y_+sh|QYM+Z4(jo*cO86t6vP363o3<Bac1Dw=T^SiFd>s5@RwQfOi?M)rz*glb@xY^#8E?z)U?qQ0!ma&JneA)zFyB!0QyLp2J^$mhj?2VkQ|3wFSDWCvR=WGA~'),format=FORMAT_RAW,filters=[{"id":FILTER_LZMA2,"dict_size":1<<20}]))
    w += unpack_adjacent('{Wp48S^xk9=GL@E0stWa8~^|S5YJf5;2B{Em|Xw|fk_n@Hymh}<V5|hP4P0A(i!4_(Nm=uw{@FVR6daNE=O{}nmnEO!yZZjleQI$`fN#h?Zw?7ib$*yf+)4nT1c@YHQsf~P3NwFbs`e6kY9b7XM>tTnD58W=8+2z&W~&N0PZ_L7MnK&>FsOZ7Gcz2GM?9?8!S_Bs2h*&<7f6Md!?c(6uCbSV0$IFiJDi57cLwUCP&bOl(nS#JyQM<#7?AVRUr5_S4dtQ5L8J6<1dz+o>()+v+iXrbxTNw<FkmwfG$RA`{Asc<Zy`^ovcJFpXlbIoj{E;@jPVs4G?E4*~DY4(9#!@$mILXJfIm9L5g(X&R@<TbX|&vL@@X$kD@x}F*cU~5NeYkW#f<uHjLEVuL<2P%gF07h3~;0K49$&mWb17+fU7s&DT>$xpH01FCyQ^zc-~(`aYG?FX`9$%A`XjnT@ipIaYecf6E)o<0_q#Ctbr1Yf~{p(T?u*aX05A8$(*(2MR)|UXVryHPT(_+}BgQQQocG;lhy%e)eK5d3R$j$BQPqXj${n|C0H(=*D$H$c5uV3u=c{?l-*g?>OAhw_ny5<PPLm)$>U^jB=c_!C?@o5Qocg^ORc7t!|W|4n7ND(NsT6HCGsY^}jD6N8b0VSa55gzb0GCm2Q-mK7CCy{o3_mdRH2KtnW79`h>gSmHB*KL9aI+9Ef{y{sq07kAc*7-K+n)EPTu5>&-UkAUTzw)X5sDE%I77heGDn0qw5$%dboi5@N>CF?%|`Pf>w5rC}}B(D$dIRB=2#ULGP(?*=&tE3)Vsgw)Sy>raQB!flaa$~BK3w&kZ43qH6P1VVG5ysii_8;Mc+zLZ`B9O*|TErC~+l}e_d68}Agkknj1ocpi!G%@nd*Uu1F20p2!nNIlV&VO+5Z}4qbHVARtr0q0q2w85^4@8%^oOE1F?m|l=F`f3=X#o~ID-_d@azFaTt2RT;LHb~nbWo$cF=;y(%>1fgXR&NRR`{^L#Pidm8k4DUDzy7IpX)n0F*9syy*(;{{5~?4kMxtU0ZTr@PkseNSNc~wsUhnoq%0(E<~w+#brSHuv`a<%1J2V38IR&;7DU4ifJrN9dbW~bDNR`Z?qjGILp!|msb;HPu2EWrrl*=>Svyvju;W%lvBm5#&Pu`E#q~`(1h5y_V~}*(gn1e*<3<umNrGRaVhpkZ143Lq(d^GmI;tZM%waZ2#R0iSFE7kjq80^p5b(ld51o6iCH<qa2*gpx*dfUu*OK>?u^X=jTMjqiI&3KZ5~UEiw{pa-`|YV%7THq38ip!0EzH#nD#<H>R*J?Q*{V0?{7ki%!|g@xQoS%R)_hT;yAE|&U7eSc#yhaDlZEZZY%?RDo2D%qNhgj3LF)+o+2J<W*I*!d;##l-&9nlR&t>t3N|FOP?YCu=#km#TMVMwtrc`?S7F(ud+8_EWE|ehCUTWN4#8nCs0YboMdqStCYbu>PFUo6L&k(yl*y(1&8X95~(H{-twWVk^p$I&Anfzgh;}{Yke0YHXebO8s%1d6Z(7Iz&9t|x^_n5nNm@(d%RK_TYt4RKqTxiJrd~<KSy3>24yC*KXp%9=?PQp<jEZ+;i>Z`4P$)hW0mv9%RC9%c<><JI=GgXcC6_c)x5m4QwV`N`5v1{P%J6B>YC#&8N!yoI>zEk7;7nOOOlKVc9Kk7wc?rIj|%@f0-Y$;1ODr$tdiBt&8?~3AyQfa?K@`NR%8inb(C&ybmZbE&YI^#}`$C}QxMk8Q*+6_#FR1-17aCxijG2WvUD5M6`ZsLI!VVV3c6{<Z^4ol&GYV;F%i?O+I?3ElyLJs{|V8pCEGG%)^cjd$$kV*D+b&41F(G9`qDyGlq^TDp?S}j^lua#B6Vuz4M{+rmkzt9w!Ol5O&<NoaQu$^08?!`!%8|tWvhuM^0^A1u3>y3?_TV7xIJYGlL+rPf6{WOzETzb&Rqj0ug0f~7`j%}er{RIdLdpmSL<Ej08H}Q$qHw33E<AAHij603Zi6Q<Qs>u1NDwf_HV%8=p^YE15;N!<>!m?Fa|A3(UImhcM0w`A{GQIXQVJnVD17_8s9|_OgcL&NXIMV!3-W+6!F)No<Ov(E;p;)ec)fk}LRTto~x`zR`^gCsF|NEYxZ0DFZq(ataBE;aCV-|FnoHz|Ll1YmMFO($SWtNHfi^^Brd>brCI+F5Q6CpDu9831a^dt;<!HWA=n@o5&w~^+1+(K8ps0G?rd5j=h5sbS8F0bLDN4gu=i3+Q(b>G8n>}4K3ZHS0!<VeN?i1Lo9AY4}B$tNUzRhwz2hoQt8u1MowIUv)>MB+*fc@~FwF?!hw7lfe$t%&$BLm?LH1+=42^f83s%8jgbE}t9FUNfUmGjxjIPaaKt5hO7Este3mq@_3&y*Tg_UZsPz|K;L6EuM*9u}4ofPLU|Y&gYa9=2NaByMe5+-BN?2d|^L${dR5y23F=f;382zeKk|d!blvyvDH!+842!d+r=@R-@||Zxj23j)8m^V6c97;R5?~(?7nK2P_&`#tH1NlFLh(~Cq?_MHiXuHy3an*RfXuMtL(Jn`H~%%@6x{JoS!$ui4in(`-#WE+%y|pHh6b%3J-mu1Ad^i(rQ<x0hPAvda<+G9wm~$WL&b77UlX@EZ4i)YG~>D8=J1wEkx4#zX34?R}=lTFs@k*4n-ZGGoQAW(Kwa4)JP9OU=n%S`|%5VBl4@n#d@IoR$F(-_X$ir;nwUC#cba>_xx4KIp5>oq55pJeWyF(H7mILUFPk{MY=?Mfn(pZh_UWuEr41euKB{Fp{HzZcuZhkfd`jZYg;Ei>;@sPnY;i1^fo<+B2#sl00FcS;xYgL!Gjk>vBYQl0ssI200dcD')
    need(len(w)==79,'witness count')
    C=source(L);N=0;s=0
    for i,(core,p,perm,b,entries) in enumerate(w):
        need(len(core)==3 and len(p)==8 and all(set(x)<=set('01') and len(x)<=L for x in [core]+p),'prefix shape')
        N+=literal(core,p)
        need(all(conflict((core,p),(row[0],row[1])) for row in w[:i]),'disjointness')
        ii,ss=coords(core,p);n=len(ii);need(sorted(perm)==list(range(n)),'permutation')
        ii=[ii[j] for j in perm];sp=[ss[j] for j in perm];cp=[[C[x][y] for y in ii] for x in ii];m=[[0]*n for _ in range(n)]
        for a,k,v in entries:
            need(0<=k<=a<n and type(v) is int and v>0 and m[a][k]==0,'sparse entry');m[a][k]=v
        v=lower_shift(cp,sp,m,b,T);s+=v.numerator*10**15//v.denominator
    need(N==5976 and s>=26772897894,'adjacent repair certificate')
    return F(s,10**15)

FIRST_CORRECTIONS_TAIL='p)&`rv_P;K@(ToeH9?(s$hCiKxHq)d#SL}Rg;G>P;dZa88!9iPS=9Osep&n|woXwfvF$)EXm#i$=-Os>$a6mypN<=DG%QV6e`=BRYPcLY*k(ok2=0vyeSNpE1jK`n!=qmV*)3lXPuk2Z=$oa`@9WgErXf}LrfX=fLw(MqrtA8^|Czh5hVKZ-+TulX8#5#>0x)cJzKO;L5i0s)k*F9()*6@C*TmSo=0HMY7SHHrY0n;`?nOaj#7k_2p4Ecbu&adt-0Z#~u1#{>9u~hql@ZYux*a!|2T8sT?me*OY~v@jKm0z201{{@z5kSY4=KUvU<y(zhib-pYC5(T^obcY>*jAJ^Vot)jEw%^kKc=d%r=dRgysL>U$cTJeax)3%JYuC(<4dKW%>hryS;u>bId_***@Uers?Clc!WCL=G}e=w;wM2BVCvBbrW&x5GpU0(%vs0b(4eRMWW=np|LgE1|OU)-cYW;wuP)}?Yw5|YMnXoZ2samh<7Dj2n`UO)y<Bc_SOs&Fyr``4|{<)`l_|SRbkg$3l9y2s&Be3^1^{zPu6Q`DUdTDSFUJQWWKQe26wZ{X@{zz?VSdj7r;Y^j2_p24`+B<FeQ~08#Sih0z7_w{{Cw@m<4g<Be{GNl{GIzrE3s@m97ksc&*9IW_UriL##gQ$@3l7p730P$|8WDu2ge6C5#TBU3)ad%u*+u%UyQkVt?8z$)e;CRl5dScJI5CoNL{~=YvU>^ceaHLm1?iD1O8_-T?WXU^;_Y85(R)+2T;}rR}x()i{w$LQ%CyJYsoXsV>RGW(!9x`bCxoSE+a(nu9B~q$<r6y>#E%5YE?%^~Y1Jkr!7@;Nom0oi9xtjW7bg`4zdgz=R=w9<ql09=9=Oe1&FHT8}~TWSw0>K#R|1k>X|eMH%gw(?3zMCZ<aJ7aPfto_mj_pqtSpC~x?Z0@vK*8hKO&i}Gx0RYoc7@SSSh*%E*-8lD=TOxFZ>*}^Ne>U~&@sJ9L1NaM7PK>&#w3KV)!^6j^YdAVt5UrZT(C3km44L`xuxc1bL#&PtdJ-l*;i+aN79+$1HE_7o>V$UaETMxL;9p_v2QXJRB=cK`<?xV)9<E{Gi{?NDPw}9yMPh<hm%f|-Ih9SvJc`$J8UZA{SZ4-g2FuhpHd5|SB_Bc@-H{ZgFsiqriiCrjlheLLdZwrk#)Mr%CsPo-KGBJD#h>R1A19fjVfw(~SO@tT$P*ay-Z40RhO+e6^>_J|J?)RH46S%bv8<goPypqoR*@!jCpe4^hN*Wx3<spKIpff8(6*Cs<6ZD@oTL8|DzHenxuux6WNK}PhjmPJJ7k^T!kCeE?%RG^edXP*p>#76g;g<6wxpbd-2LGLO57>t_v`e*gXxY^-n^!}0qe|~CF>}CTmZoc&AO+Y5ifhrNbveAz+)bM3#L6|q16v|Z3p*R%1BcQ6lXdpSUN2E#q-wm(W6_f$#NzBw=09jt#h{YV)An-3<nhA;7diSQ;J>f%=8fjlr2;OI%~fDXG?^4ujL%0#-(evoI`Sz-Om7#u?-S1w;6l~WRuH^7C~&w0XKmt!=z`jhRE+5;xV<5(Qje-5WNaKU9#`ul2&fSzj%F2(5B_DZ<qsk80iLi7d+XPb>I}_q7;7#SK3E!^;YbWZRzmTgDFrtrF0#Tv{EB!dNOT5GhL6#-W#oE#cZ08u0b3I{36&1#fmM!0!>+mNLwxxFDmK7^2o1ERJplP1LNi#?I(!uXeXdd}SPxM*2|3elSuRWyawnQHVZH)kd=<;oEHUWtG4CDA@EJU919du=Ax%Oo9}^_hrrA$F^kdi`8N>`n2IJ<@Dp7d|u;I$<zJiob9iniUlkbPtt9m8ic~BmoM!ap2A-Nid3zEymc|(Q_$gSWHy><p@{kVe?;@pdyPD7$ibQ@X@afke_mBdo0lycX2f{WIy!Z3zdzR#~gS-xpx-kV12N=N%L+4&plTaPY-8$qC#{kyuHaVg|Jadp+|Z6K0B#f`5*Y}O(`W8NQc#!^@;V8}Ish<Vuy^UZq938i<=2~OF{oq%;a59g!C%8b!V&cLa1wv&=U&cXZ6-32K-fe{S2IXgNq0aTj(Gj^Bmqvb_W81t4*aKDBsd9EcX+3C1AqxKFRP~g+*QFW?@p5A7zNFa4S8Nrj^QbyOXY4+E?FpV4431g6o42LqRbB1@v^>jbS5z^tBi7)61z80o1p>@)#s4%ru`<CtRYA4sSsBgM0E&Wy0wJ+SON;50&E2W<-FfCvi@7_{qifZyUBDx`8&gmeR{bV;!Zp?7z@zHnBT1nFz0NFlafl2?DfZ!E3gjXe=<p94YZ1#Qkkn&-p)sr*`mq7<sIedW3h9YhO+oa+JB-fpAC)W+U+b-yO%D~sv{uB(W1S$0(Of-qW6mpuF&KQR@ah^oHLJBKJ-(#?BMR%P!$VKbSUXWM%f1fXe^q%f2Y~z4BCOd!u2+!F;Aye3=c!K>Hs<@P12A!!pCh%(ml33Is@}X$B*&1^uweT@DHXZ3jO>{IBE@VkIlF)bjnjUfD2q*cj?wlly`I(j-^G}e$$zO5dtggw}(XT8;Zs{);MX7WmHOUo7^eYF4U|~)FS@4~xx*I7l$&F=caJb0}fL%2-3GEN3<f>llG^-NzJg{c3bdPftZ+RoB8vA%Z+UqddxiuO?lrq62X-R0*-m=ByBkuYuqcu5^IFbqqnLA<bpW=%F<C-ZzN$F!xy{sFCM1)HsIZ#Fpq*t6CrF!~A&(Ie~(RN~i<02msZ)juWipC7a`JW+w+lFLQ3zQl!>!a9-)D0%}3q8LxJ8NC7{jpI=wtGY<Hcc^|*@9>*=l=fb=Ir=iO;OCD+`0#Y-H@363PrhnceBM4=lITp!w#3a>xYCKpN8Xj3xJK5$KJ4B|MN^tePyWFuHhyVM4)AGvP<4YE+06^GE<Dau!tV)YGLWC-_z35&_rbeG5DTsKKW{d{AniCW6oTu_VTLR1%iFYHU4K_s^*Nd6B$jwbN7b8(-)coOUNdc5<3;&kK^GL5P+3Oh3SNQEa*-DsKbcj{i)se(Q3~*w;heiHycwOTP(#$bXd+fp-$!jUxcJ$V(ljn6mn{lQ-~Rz+=2LPD0co>?ZiTD;C|DWAVHZt;ifo+-{0VHLJHM$HKT?NX{!ke=01h5kS~{$BkP+-<j$3#P;O7w)^#N<q1+en)@#%D2a?XB@{W;X-z<Zg1=<7`{*v+McULeq`qZM)FQ1h8eMHko?hbgP)OgzOS+3Lm8sdzVp`6kv-Id4+x8>DH)VTCpJ+EKjx8*D=|M-cBm?b^U#bmb0OktqwA53Ye7+XvoHVaiRnR&6XBo%^fJZOlk%5BHt8u9ke*a9H?P}(Td*m`(zMDueah^E8k10)vu`Rm8y5MpS>n||0Teqc0vL&y=h>m!%*|J2P?Y1`(`#=)v?U?uWz4bb-6m7KGuq`w=wr&1p0pZGTid+x^DWQn;QVadARD1X=FI9#Quk3&Bp#+bX41-@+I-buLl3P|Q`8HvSU3t>px6atq$CJJn*7h0#jqF?=p{n@#HURXXNp*o3#=$Z?InXa^Z$9(KOIa1#v=8+&vbwTCX#AjF~(3Wy1=WlB((uHs~Z+=9yFc@y;fbce-<*`vck4dEY>z#K&<eh9(p*PzyOvpX@wM%6-@zDki*(Z%ZwAHG$EyEE&lWnCPoH^G%O!=FMuCzlv$Y_c$1xjyy1ahM!uWQxl-7|IQ#$B65(2oGq#${KSPCu5^=xMUuQCPFFTT>4aEQBfB8LUp++c2bXWD>K~5Ue)zNsaO7%$YL^djtP*AdUg+#A3VqSF$no!;=aM;hK}%{-ZWnjxrwJje+T9W+ax)2Lojm>FQvLH6i*BzBI2C>(}z$OS$syVL&JqD;z2}pDo5sEYeR9Wr|Vk;9w>HxAF3w-IjN7I>M3zU5<+_1Q3YEkF7l7wPRhE)d0q_0k;u7a9;i(y#|!hu3WW9i;u_3FMBA3HF0LNLe1=#oufP}fbSHCldEoNKYr5SilZr5rZ>Pz<aWxdledRW7Yj(w<oin;q(NQx$Em5-GM;*UZ}V=;e<?hn3C&UPVyZzP&webrBz~9>>E~}%5gu*8>a%@R+<Dd1BSsE_zOLel<70Xm|5uUis^xfwPH0oPZK2}*hF$EX-~o-#s_|RM;JY%A9Np8J?sYYpGct<Weo5{`x-M}S<aM%DUlY6ho_1ipkQW%+@mjpU-F~-m`vxmK>*pP0yBpGYs`+)6x7mqJG_a>95E$zU)MqUR0)WlG5HN*|k_XC6fLnW+49K5WsgsK{L<gu1y{Yi&V|f0-7Pl7+w3LdEEky(=;EM4grr}nY$75ICuyU3nDA@64b(+l9{0J<uIB~q(9B>;DeK6*VR|)?;zAA@luGKp9p$T%vE1I;O&&Fx&Evj_VMV&2Gerq|$t=IPM&LN<**j2>+NAqtRWwm-CkN<yadwr5I+zT3{uWlk%Tn58O#i>fV80~XJp7VpF_t|Qcpy-Dptn(PuEdB{A!4eQmjO<7o{wOwjIzi)n`j$g2($f>rK7u2xZ<DFTplM>1Mybz-c5-SMLNmT&N<o+7Bbm+vvfeE((+R0rt9o(Df)3PXF#ee?;`k<jb9%W3ndz4$@B(nrlOIGwl5?)0_A;gw=JV(7sJ`>cU;$m+)xvXsLYVG-zQX0)&)ycSD2OyIM%ju!Gt>xaXt%`W?c<JEyCz1+FVhzs*sue4maRqmx9J>)d|l`>C=LGLxo|1~n6r!**gq>0lL_M`U%l6ybdVyagii&m_wvb+_E0uU-zkEmc2*kOxJr5p{-1Z^H~M`Skd^Ib<AogG|G8?9G%?O6<kyPb&=FF*q)u23-8D5I%(}0_zc&o(<_j^&tbwyXCAYYBtS5tyBbX?*dg$puTy~1g+Ru)a<Ah(wAmMP>Cr#}#nQjl$@lRarFWxTe^5$ai7u;F?eALt}Nol$;-i_z&pjm7T*6jNGRSf2gT`|Z-a3F0N!hSo?TpD`r@g>!PI~RVm>F5Yy7D)GxZ(1X8elR@HTK}Q3Nv$l<1RvY)DKc95y!GQAm3{{Y!YzuFWU-?4jPC*57P7*A<()aG&*tzx0jWVzHV*y}OOJZJeHhY4wl*GM`pT6Xl-J%nkB?_CrKn{*GIbJ0`2e4BfWvg3bP}m?9h_e%4XEo6l!5Jj(Hqg77WF7YU71jtqaZ%^?s7jrSr?01*Ss9E@s}W#=4To=Ju4<?++wI^%+f3nP+#ewWRtKrm%y*I^vKg>K8Os-Yr-%pAs9i#m(B)g1%1W_OSru4RHk6vLvvx6$$!TT{m&4L@p1~tU9L*M??oWTiqGT5fv0MOzk!<N;{u~`{;JQ_RA(|<9Zd>Cbg#2eEEVD$wcx4m)2Hf!#0uBP2kwE1H1U=L<->ctzTOaz(M0?1SDY{s7z%x1!=6nKTc5!XpzZ|}WbLC_<=W2lQ;+Fb$=~&&bdw_zi9XZLy)^c+bne!IvkrDs_-ugDcSFl`HnrmV_WQ^~F#_9W?yF6<aU15xqjZ6nNg&WlB%slu^D+X$iG)&e_(Z6vYoh3m@U!A-2ZjB13Md7&GsRDD*UHlsFYLW7IKgl$&c^|OEW}GufGRR1`pS7-w`f*~jx$-FZ@E>Fo62o)*iCK~R4J$l9e;H6h*ZSa+`qQ?d$0omFJM^M5$!clDAn061ptQU?{4eDJJ<E1pFYUq3=yiYmk%DB))z3)?mlortb_&8XQ}Mp>(O3XpEI&L9!~<L$ezy$Kn^@B-9m-(u!}Mho7fa2VSukcU|!Q^0An~}O9Gy0(U1g{bYXa$CEpMz`)HUvz}O+nj`h}09l)+txEG4hy^4T4AYcv3P|3OUTlWp!7c0n~BGOvdbbGok$Q$uQn>Xb!8;_VVuk@a!AV?OV{f-Kn-nQFns1e)?Cd9{0@Lsquman^VA;+zrvtBZukYZv3vV{>8S$~Ltg7v-19oTVql@n3_Jqs`fbFfM=y>LQQ%ZH0D6Ro!#*+~OL4;M|xeoZ|jn{LYl$10Qw&f)eU%Mk@-dWz^e%NeqAi@?@*>}%!jL!kHqJ&{Mjqe-1MBn~O3s5&D3ZPry#MAF|h7d~z)$sR%ZYti*9KCzpYH1sdVx7#5Fr7g5no#&W_&L~-nDl54h7<R(*O03A@?VCh&Qh1W-ow=Coq5Q_drwsmU6zVw}0*78YC1739cL3~}m+hj4Ti>2nEG9qzZj`O%23_?_v8cAX?4*0+1{(I6&h%A|@LJpdqD;c8Y<A)G|L><pWHG~5AojP{DU1Fni<feh=VZU!wi}{<{QKg+VUN_j7eWsECQQSTP0i8=Y2*2%Uf`enfD_S!JuLPJ%PE&Dgsfmn+eYsy>U5`c=4vqJ4d$v5#}wL3v$^XsN`_Lzj(b3FEWYgTM;W-K0&_rZk%;nYyUF&WMfb5DwJt(~`%#LIbtgs$xip_-BnqV?)#<0lz*UpYs$X#(M<^wuPEi01F0!$Qilx8tctgPlj5CWuUSk}4FNv&ehmy!YfSe#<45(VPSO!r2c1Wny%M%CdI1key;Eh*r01jhLw`TMX<%>VznH~TAN-G6P$3+(6ip(?}jYZ8h%%r#Yq<na)SQwY(|7kqzpVm4g&mVPE7xA|;sMu$o8}C#7=bghu2FSy8m`3M9^-xabg0!vt@y@S0--a`mBxM#an<G`8X48ehAP3-5jhj>EXM3kdK&+xb1wiR?IRs;rER&+%xs>wX8!klwM{~IF68YH&y~|_mb@K`+9VzX@8##?{HDJEnjs``MskqXD7ffJXwO?tH9p7MP*Nz;fz<VpN)pl0W`Y@X@Y1Y<+#j3omW7JgKaCbrOFfch(30M)viMwRzZ1Cw!aC{A*beLm(IONz)Ier-v5HdlqM!M51@ahhAiOagI1&<La=!+GScn?Kma-I7>+@AQCmPZw>$%#<gHS8KzVnL|U4X0S9?kNy3#=Pz~@?FkmXS>yoxYwX;$)w9E0aan=zMZ@w&JXj)%@h8`&Ko}b7$VFL45GS^V<Lkt_!ef8D(Lk1I%8r>i$}qbMX^+tjSQMSN<!7&!=ysfaM$#&yCUVCWTFCpUN2{yuPqJrHiMA6i?g4o4_BR@f^lwzQ@o@^V?I)wKpjv34Zb{B@1oK5?O1Ehh>2P(%o~WO8w0<dr#nkzqeGNodn1~w7(0>N+HS={HE?^v*iD^vjxKRtS8H>1%k|Ax5k}mLvZtl4(>K$iZnqCJD3fY>0}|kS^@Yzf`zsclY#i3P+G0M`yajOR_?hvOfiSQB$BL)5*A-7E6=|SEnEIQk;<|CK6L&;nhy;aMVF~9)#=0=Lf&ZV<P@qG38%^Z)DCJAXC;bUo_v|q`d>L=4k%6O`TdYCJ=K5J`>65?AaXdXA087dWNu=$>k!JF!N0$3_wjhF&Yv5%A0>jy_20yvfw6ImbeJ9^eI2WA-aa!JKV_y1>#-<dEVR;3Pt)>Vk_e3=>j(fSkE=1D~QTY>H#2-J3;QC9Fy#!l9a#~_5Y@!w%UNyEL6rcizP7ZP45ymT?LBoh5IhqF;pf|u-&tmo=K?gtSiGgdk-28Y@v-U;`eVT;TbZ0&d1LN3X<xi&;6kP(;m}2lfaUK*ib>eEyGBsW4cx1cZyt(#r-}NBZtZ}GBX#Mz2!|j-`P2ENgRC{b(y4W36@=Tukip&$?pNQp_!cODna%J{(EcITe@D&%BhPhepaQHIM^LM6ab$KgBz6cj&YYA&>#ch{RS=IW6G|J0M#n+5CNEdzq2?<3n!*djb!XiaAbCW-pb1X;_(~A`UZHrM%f52QW-~31KD1|ecy@z4QwrQc<hlwm{9jYLy`pl)&c6}Y8Qby9C_m~mpai;#Y(b-XS5y=9(K-~N5?6C{=j)(lFWqU~^#hlWu2P8-@rm#}j1+UGW`f=K-`W$+V>N&P8D)NSLA8R^!FR}u7mie@<;+y(YRQFBoPPbY}Vle6u2l(^4wOHsu*G}SX*Mgtdt#s%yEF?P+VTDR}E2~>aP~%DlXPT|C!5CRR6Nq|(p&vT!+qP^5TO6Dg9T_n@f>%b{245jXKUX7~$_;LE&3)F$+o7oka5`irQH=VRl-`9g3tv6Slk#$cI`mI~?NQ7KR)KWb_v(a#2;O>w?h7wm2%Bm#6B+vt0Va5ItSK>V0Hp0aA$3>T&qaP%X(_mSLY=S4(=Psr;H%(Jvpx*FL>?Oh2-ZbQ6KewzVWEHvN~G8dyKS}S7JzRB@!hbfL)~aq+7K+N7!}sLX4Z`yV%zgZ7bb=AP<alFOx9R(Xet*^NwTtGR1H{v4LWNA<+n>pek$J`^<{mtrUz!O@U&>!HbjU+D8Gt@oj#`W^aZc$EIY-iJU}o++qI^>Alm+{25wmIK0Rsu-z-7$2Spm=#vs28?#k9}WB_9c*GSedrrLGO8AKuo5jMu}`%#$bb)?R_U|#iNJ|q6m4O*36UFB7nw+AhYzC33zYh1+OJT~g5JKgU$L_Ve>sG$}6{-PbTb-`p;ngGM}O^7u*V?%X4*SfnPdG9?&-kwMq8*}9AyOe_}_@cuu`obXLdr2?rHeH%ygL#zUZTG2cKs9@eFLTKiXac;JV?aS2TC(D$y;%Arn(nSdC_i)mIz_94?g>l)Sb0XfDVWxK8}`bW6lcc!Y~uh$Ob_VdRjug%l;u4#nfgQLEiC0eqQNG>pxtcJ&M1HjOqStk9i50;21uQ}tMy%RU-X`?NbA-X(*hg^u8iG;Cseawe>8N&*}jF_$c@3qpF(N*p)+92uKIO%V@Jif92ErmJryMa)wHs4vag-wI?&hmqfB(zX`*0}vfsieYa1k%Y&1li%B6l!1?>bW`6%hS^8>(ZBDs^Lr!iZ98?`t(-~{Zvizwl{_@)!&zj$;Ryzd`0pqmq52X+FHAhEh&0laAS#SrXW8_=)w8eCVcAA7FnnQAge^niTO{t2W%t_d2>$jhv4XvZI5&f9!k*LpvYxFW;Nsd1s7X)6zvO(sOHEqBx3GVjSjJ#c!7-Uwz{1~;t|)~5G5o^*|iNW2x_2V?ktI)N^DjIpLtqxTK=L5!yvA}tj7Oiz>t%n|+_8Dg&u!BKu+I1y)!3B3voa(z-eEJFC}cqNo|CdeTjd=kUFj^w~cNOhO5H43<mC=fuWrFi7N*Srm0*%{vOFR~D^6pFd0DNO`NVmQE!$A9QwmBr}#`U}-+%Ip#;+RkFs*_O{{B$-I0`ZZo>zQSTU-hBoekauOH#({*^b-f>U_0F^Lrgrd3uIW`;^<==AzNId*M3*W5!uPqAjRh(ce21i3+2(NW7;=c`xjCm!aHY5wKVEvsg>vO4OMlKStP9U1ND*c3@Q@z%RZZs|$$#2XQt0k1(zaj>+;672JQ1W&6-CZeL#W$^H<0G`&YP&52!CUTHr@LIX`TR-wgqFI0rI$%j>aSt8St4(FS2<$as&Z#NhyU}BKtAO-bhNZU$rTh%K&{1BCc%xaTlqB4Ez<14c8d^2x<K64?rC{;!+CI8zKG;5Jr>b0)2_;{91Y*wr~N?JCfp3gRPEnPEI@yrhn$PSe~;#Qd{SXdGk=AiQ%@y>Y3YO{1WPSwQto=1N^Sne<Dae{s`)`WOT)0f#9p#w&&%!;h1KrAVrWM&0N$Bpq1b7al%cBayCS|jphH${CP#2G0kEC_Gs-jf@7EMXFib3F0_gd$Bza>9jQ8>$xN(7V&~>-f7xNY*#>7dU?{MtXIsDVEu1b;;Z)JuAj%|>!2}B5P+cdCMS%C-;#c_lZAAxIP9+GY?dCb*M!&jdA3^Rg(>>!XaN<tnoF@_=7m=!Fv;{GFg4fGjK|V>j0?$6=%vWjEESmixy&5+I&cOWA`2KbT#<tMZhg~I-*+fW0!1`Fn6|xi`mEv(IE`bBtcOcD54_qUMb5oWcIdGaDj~pfbyC!kOY$493DCfaIA*Z)cVVvL%L%U8E`FE5D3+Vfx{tmpfbF(WIrP@~F<khe}3k3R%6Xp*>JJR}R>sT*AqjSo2aC<lDfX?9&mMVo$>61-sxLLFnz<3A7f(UD?2+()?zCZoj#KRBiy8rj%#j?(G+od|#*F`nSdw=Fq2>e1JR;70$2mhWwr=!H3XICqG(X3PtEiv|7LW?jQ8^d<V3Ze>P{AG)CHPS3l;B)lS9<y^@Eo0DPma|@E@)Q+l3W%1<30@gR=>s1;^0i4<)vnJgK<fbm>cz4EdzAX;^L7%)j!N<@lbK>xW-2Qy1qYELc?>3LEH^*X*ALCjQ?PLEb<5iQS<`2@R1)E37P9@T%wOe;M57Upp8S9@r%oYsqk?N<G$2(WxA^Q!0|N}r<SC;BJ}P2!ep=@Wu4Af|ItN+Tn@pwo-6$!XJT~7xe`;d-MXh{r0dnFWkPIaw#o4hpzp5h~LC;?ktzd5Mp5V@v4>E^!$J|-0D4*71DA~c0Ha20LAA)AUuU$+wSnm0~6oq0IeO}aYgPDZ6|I799s^M71yUv8_03Xfj&*vXGX<fWOk2*I%;IWE0hVSiygqm|XwB(l&HK6raVi%dkx+>YRV$I0S@1pa=twvGr=#ov-6K@?1mL&4|inH}VROTWja!b*DX|O~ugnw<4FKFLowjQ2qY44GxJ~v4?7no@5gFVHabhc7S`8#UyXc@^4lnug9WCpAUT5xI)>&RAXCQ%Hqc!=}X5Wq&IAZ!LF-IuqnWhn`SY<v-Kn5c0Qjs7!Z3O2@&^$w@|Zb)iV8YPj?&C@5YfDA$gyM&R>z@m<|fnpJ8t&Fsf<A(_9xL?tG?$bRnOvO37;sb?Ddut?c9sejUq)8So`JNj1D!pI?gn2*h-I`J*oifbadc}*PTFr@-=zUmEbCmj?I|oTS`z4kFV8sOe&IBUT9$(u6Y#g<2RPV6RGJ;!<N#2q+Lna+i-aqqu4#jP63Py#AJks*JU1<^2X<?*w8iRi6!jiBS$O&3UM04nGl=rrzNe?i0Th=Y9mTG$+K?Avw2+P%v8*(9}`NWOmeQJ8M1$bOki}f??gM#jBwU!;vjJ#zcXU@IqYQxr)o)rPSDqbILh~A>=#HkfqBh!CEP6V=UP|M(h&0gQEgMI1@6csTs3NXI{We`6nQ*dw}FqcNP3_q%b02+;eenIBmt%hK*GY|bkUG4K`GUt$RX6#M$F;~_p6LYM@`*z{&_h-P?TPYsT(&_A7W-D|WBkJBb=t)zdKgm|rn&ut^E!ldN#_Xr{GWhsggNGV|k&A{RzSC{k6KN_ve@Ahjpq!XaW8V$6lb&@8P|<4g{UrLc;PP=-aAy^ztl|<zt6ktQOK=^9uLXC?w+~T(#SiQPxZ1q3yw(mTgh()`3#${XApF!NT2%OJY`{TVK6`i#)_3=%rZ>P1-#76W%*HtBY)5{#_!23|*3LS`d$cm@@o~xO4ef?}i1xwEi(55jKbW;A+FGG0e1O?^<-8CTKT90bQnCAVJH%CVT{)}5aEFLOP7AUfM%D9CF=46IykJYilGr2OC8k#v@FY^fk`Yuu6e6&|q^gSL157t*st~{&3cK2Ez%g8Fa`!-E4iUD~wny>!5}nW7gf(0JvR~_kDiNu~Mv=$uKwfJ?5T~tr>vZQo+QFJGv|bN+Bq$`wOf8-If`YYzI>nxEzsRJ-sgs7Z-5!=4xS<l{q5E_22abq`oU?1q9f~G|w{CP<ovA#q3n7A?iM0m?zN28yp#K6G6B)GHh;FJ1@H`q0iP;>On-;RcEO<<wPBhANO22ZjfbssKjJ;hcWrMYQo2H?j$|mZU@HO=(PsCc=IiZ@?o<0z=)Ce%#rG!<k;l31W5Q>ULS`LNdw*Yf_4cM6P;k&ApKF+tbPm3ORT|wOT$b@k1g$)ozXc39GUDM*g*%Ri98eXlM*?JR*uP)K)1LQ{=RbBfSqziOtXUH$?CDvs8_Rn&pphz{<xfu3iPx@&v_admj`;}*89<|P%zIaVK1kMO&NwME0xSNWvhU>H4B}Da1<Z7e7-6-Z2am}f0s2c|chFVgf;rL<dpx<nSD&y=U1$$sKuKh4VZkr*OlnNwfA<aP6IMZd?e<s|>7z_E@Spc4nGh8wzlR;@gBby&$6Q@~ofXLV0QgCC|ELmPDYH|WWnU4>PKnbQyxop3q7>vuYoLTEnrXug7(x+(40fcb{A5rjNLtjAmf%QPcmcM5JfXmiz)+2J*q|;y60%LkOF0&^o%2%s?gQfz}>JrXMaa2BXk91iy=B9UURlqe0H-Q@(=7BM3O-fUchj+RCdZjBUaSbSyXcOd9_sZooe?k{Mx~92AX$GIWwmOI<1lL(;8Lqo;>gs}f;}6!J21@A>m=lC~&r}jxh9ZF_f3Npytn0xe80T89UBjLUNP~Zfb>3m5xihZbwnV-N>zi$X=%99&nTv<*F4O$Sa9U0>;3{rXSX2~_@j3U=-VLfJ!jrQemls<^f9(#?oEo;nbkFrd46r4yQJ6s$E;`L%^lNmlccNPZH{OW#-Ws<u<)qn^abUU97}$6N$pB3Pug(=?p1Atc*a1fKb428t0LB?m!YVWDM+&C@jOS*Rh*2$0XpjTnMbWK`A>!g=9xTThWu_y!if}pmtlOF=mATkwimaT9msYU{YKs%1Gd97jqsJO&ujfmd(roO>l$dhiB(PY<x_u2On1r5nOF$5Jw<re4?jx>dA*)Akz(RNq0lu;592sQFR1<BtB{t4C;R+5xIK{|xP(@UibG<UIUwLnu6AXF$lHoJfb`fenKiWRnUF#)Wk<s=wd-ioRtDI+6Ii?y1SED@Wnx@HqpSwAedP1fd(XgzMDA55e*LSF^p~2h1QP>x=$ch*Bq{3GVnFAM2Xfs-Lb7E)|sUQL}&vw+AF7#zP@QEbsOEO|7`($hEUs>3JH5KuwV=43uh>>)m6{<A2+*iovzzf=nzy*Vj>>R8Ce+O8<uc1(jXBhrYk>bgqKW=ia-+_nRyRw>o8_Ic+>&p7g@t<r|IsF(KN_o!_i2%qk?i=P+w|Wo6*dOJ9jP_aSj@^?NsJrXlOL1cUcq!vhZZv%IB@suT`d|Z5G0ODK*X{6zj+`e2>Q2qsXdGD?99LBM5v2ARINv^X|D|pq_iTy@1WaYTT@qEjo{QQOrzL*|mr;UljHrA2)n+s7axYcN3(^z3gqOHx%dquBG@6!fonhcIlFZ$FHt@_TNuuF(cVLhuC;Ch<H{9P~eK6?8H)6w5W(bpzuANBfubf5cqn|aCe8fPx<Yi}B7$c*|lKD~5w-Q1Kq)9yZBYlX+koh>1IawtV&o>POr4SF?*GpdE7cI`<6uXGDFY|URCl_OkJJ<byBIPG8=>4O8ixrc@o!+a(M!{osw=zgcDhp8DGDR6AG6&(24MIjPDpPIY7_B>Eyt!Kyc1o5rU_Aq*1jY^^7=Z>i`UKN6P=6|jSH+I@q#SmN7ry$zmuTtLQix+;m(*f7pLEp87bhagqm|ymsOi*!`r!Q!hcFAdC1j(c#-L1hv6h=!Hlyz<&eIt*)0yz#(i-@-!oK_m;)uc3-2ZaQ9Xr{(q*wy3yl7_#trP7-jYz+J98o7lJo?0F+DyyyC`sLo6&;QC;zlc3L|k!)rvfuwPW56rgcfBzbOMFh9DdQNyU3)i^Xef!mHnZ}8!bG7dvu{ZsDNA~X$~R9B@%cY%h0Is&NAKe(YTjnP}ZaLQ=o9;6Xi1ALA+52f5B={<DOs>vZX{N!O|(Y^VywmO~~j)OrTGu>W*OE>ZR9XNNIK)_XlZZDI#PRB$Y?}T#KKjFrXmpm)<-n2TAQtlsiT$iMVbV2I{^oj%cp>j~kyO<l#W73a$N^jHwTE1je^xbi*MZVaA|7fQKo0H{AS-d+(&H%CXts7sE_i8AF4e2Ez~A#9z1JGeQ(Np=i|e`_|x8B1prI?CBH+pn7LSz+q{D>6=ZE`pIb%0Pq@aN@Z2tG3RBD<@55Z245O+)=K_~9iQc1m4-!y-on;i-={r=GD(>8`4|W0?D8Qxqux^kfAIZ_q}}Xp`-m)yne@JZA1uU#9?x;)AMmhEWoOl+-pbxFg7kumscfO45etT%G_>riM5D`3s$r?I@wp}h&v^$crIVRsr`<~1IQz-csXRX=oG;E1Wq@Y#c8%3j*3Gl2@R_eL42hyS6;grG6yPmHI6inl2`={WU)iZNbmS-O=@b+&?YUU68R5J>EPaxB7E>+{I6*NW7PwD@?9fg!i2Rac!ojr)2Z`)pPl;mpY7;=({w=~Nb6t9c5{uB!A23hpXo?Q+X{>WDM4%>7>b;$&eP_kQ{|&w)A6TbjhSvxNeyk?AsnqbP-NIP=y`BVxwiH(?C?cRmsmUglg8(KPf@^ga>^WxOYLr<VRP|Ez9Q4ZctTx{F2`6j9H`nZ0BD+a>u5hK@myuNQE{Ss$E&=IyXhG2A=i#_db9fHBWWYSbpS4fAfN6J4!|<nR6O`NYOw#Hqyo9aWtHZvC*HJCL^{pvN*vV6Re!_~N(>mJf$GCWnL~gTXL1BN<N0iD~69V<&$$-$Ii1e&Yu~5VGz#CA^RiWFua=?gxQ4P}Ed84st53v3E<o<=c*Uj=4409wi=djyfRR%{>1d7`0-S#L5)*JI<EHmS`VZByfFivi3+@6CR6(we}qy~#rL2VDm{p6D^H9p1o*BK3&{f$sz?-81EC{lHYj+%F}pbIP($nZD@81>+wl#Eq(qLc=uAcAtFA}eF>n9cftU78On$HrzPOJu{5CA~U3Es6@^-_Yh*0IAk;)n@fB$42K|a2=<Y>*NX3=ruzjdN<Eh<tby%Z!)O(E$B+A0Alz)!a%2H)X%!Em+s_295W=KeMJKz80r;a$hyM@$}`4%YLO(<!K0bRT|~En<8?bBK$4Qa*s|EJ-BFGDMGbJrusTjTVZDnnO<}2lGc$gUeI9oCN1XDW@jbRj!BQIOAg5cB9@L}Ayv~{tkEm&%fscPtJyPJ-0^w|CxdX<>D$3we_&kTsrF7(;-~NSjlqsN!?#Ky_($8(yV8g8i>HOfkPTjzg<>^x?Z@aMCD74|}s-kMG-plwMQw^k(+MUZqi*rVSF6F>p3d-3H|H-Eh{yl|pPueGPNtt}|2>6WfMxG;|A8R_+uO3y^xAQ?$5)4EvuXi3(2)WhJ>5?5;ecxd2ziHD{dYVR7ttuDJzEi-Ck4`}$d|Zvds5AsYJtfg~YY<JSxpv?7|6{K=hFKNwg$A@cPyiC3H)6_{?DzAE!Hp&sBx@|ji&tsHMu=)^cfWgGIKbUTu5dvp@*Q=MT*y2bB5;GjSk_!57r*>)DDt`!2`@IxeA#^SiCG<<b|GKsDv1W?)h*5ed&b>5X_FZQd9?ZNy(|D9Ys0`cGQAiBz5sX<G1sOIVY`}m_g&=|FhTVeGptOJGLI7LMBLOBN}pU)${|KFBhcOEU6ICh1(dJaqg!{NuC>+6(~vbXCy@N*aR+pJ>p_iy$7D#^PEqzy^lt36o5YSqLY9ui_v1D7EF#0R<08_wg{yii4!nRh$TnC~v@Y0ew!!Zz>@YBeuJAvG2qwW6v%1Un0m^Z=rG>=tO<8k`yjiBo+Y35S3J8QFn{PLm%QUD#A=;-A**Kss1@V4!pt59wi=!L*MtxCw2uy|r?HJ1(JCY2;Y$7-LRwV9Z9yq%mX3ohdj*4E7u}jmbzK?fTt;d(;_BT%ocZ#AKsvZwsRzF(7OkJP7)3uXoNCY>7aN9xE=AIJCE3%du2aZfFGZ1H}P|LR`avr|eGvnKbVSAW_)Vp$MkE7{`1(lc?DW{i$OkXDWN}HLM*$MW=8lbHkGBFQF#1T%VSMztZosENm_Rx?f^8LI$qJUIw^mALEDAKDT7zuVg(}`Fa($*qnpbS)U8eBKrGLc6H-G2`l0Au#}8G%!}_aAbqfQKdRMxh1uJfj`+NF&NA0B&tZpVZTo&nh2#-P*jPE80!GK4#Xel?=h$^4?QQl=le=+3#PPR|DI_f8O+F0F6rPF#@dypoJKi)q9~iuqsD7eqtX7N<pJ_wH2~jWrkc;Og2uUJe%_pZmLP3m*#T-j|`V}6dgJ>MRiG#_hWWO6}~^a8MqW?{?$X4cV<&NQ@Af<Z@&T3+ANz|gm4Efgv6*Ph{2f-(*>wrPXZ~vd9}Akm<tQ3->(^(7QWzxKl0Szu$O@rzzCBWxFRqwalu|xyqu=wt?q6@Z#D;Hz<U-?Nnb8NRaKE8H=#Da_>N{ql+1pkLaZ>^>L6f;tnkvsx@8Pr5tZbR*k$x+*M?v>EQV3_>5QH~I`tDFPAM{|4<s2%cM>D%WPcQNmUUbsFl<&eNmk^pga-qtccloO+pFX1+XO@VC;#66WGPRL;MJl}D!<zt-In4TcZ1qJggE7eOMy-ez{K8CZE&HWrIQ-!g_7rU3Z_8`510kE21uE0h54Bfj+UBmt@zbAfW~;ttBW&4Z*_<{;!%l5WD_N5)YE#Rk2y0*6AV>8JwZWm`yHPi;j<Drn2{;#TILF!#;TRFqPt6;w9^u}rDx#B3twiuSO0Ss-v(cGDN(1^U7IcpSto&6IsJ7(s>mg(dx2mN1L28Rit79v4riOSIAro=d#+BStpfG=SjFEg`02+grdEw++A{!aHGXIhxvW-uz))QlfMpr`(}@=V=vlXN4Bu+lYq79g4aA_OAwPlIIFNdXX5@}yO(fGP00CyW2ES5^y&U_Fj2hBOW8ZkTT3U71s^gFW#2r9(Wb<$YW$k?zSt#+BOOq0BMA=4K_hD%u!<8i>-pN?}@98c6sd}w59zZbr9JF115(UMWg~PE})hXHs;1(g1XdRy1ANSkZ_K5e@VjyrTp*g)~!@mkqJ$*n4n-mbo3-FmywuCsJbE$viqXY0~yL+(!^=jjpY#M=p@ZVHqhMzo6q5Wrxp;mu@zW6Yb=Zh}(_rT5?TkS@VFS2kivUT7g!6+k~MTCfBv~*)U_S5k{lRj;@qgPASoPG`$Egc7$6gFr4&giwQsEZLX15}8((sj~=XRkaT6mkXQi~)t1mNWyiVl@0`?YbEO2z2Kw{AIOajcGgUDa*ua-?9y~&i;5BQ7b}L$VrU9dA?5$d13d#2RnmB5<wM0AJOf`*tqd327k93p=Mbv`TgW~tPl$n;%3khH9dS(ZHlH7$uw~|q_nj?wJxRJ=u!;6{2TsImUpb6i~x#ckmpQYz!OYF3!uO2|Np|h!|HCP!nP_1#*OXyKjD9|V%{FEGbJhHd(V`IASucxEld?2kE&Dd(bzPnEiHfDvT!<5q6qEUhp$^_9_wtEmgL96H-5&yelB?7GIyT<_D5lY_b0RjzS1$LQaaq{#4Y6J*sp)YL2Bzwr%01}68C@%#!o83P7A*FJ^qCyszXQ>quHZ570Q8)+TyE{4vk~b@e`C<HRnAj?U!1M3ckgr9K=0+^SEm$i=nR?(yZcxF;V@{;o9&1j+(0yYtSO-AEZPz>0m#;k6j-2;Bw~N411!xi&2CPu>qc_Zai=QRiei=1Ifp&f6)T)+eyiNfHI@c&@KcCCQ=bh7*@(A4%`@1Fg4Bn+VBLiCxyUE;XIw;*v#B}1#MdaAn`pEH9AL;;GDC%T}?T_J)<IEE#-4<TfY!(=@Cy#W(?J~czfqIMGa<TuhiJXCvt>8pJn&U+0C8Q*y!E(s^vgUh3~gz2xG4PORm2>=q8Bp2St`}o_k5Vas!KIgjCsTj`~x{P)jx8eek^&nfQ5Z9$bp8YfW@3y)K6#grn#@Gl;V<kB{|OGBNk%SY3dT&%RmrhR5I3=bdQruW0eHgbZs?<M!m;yV1$@!A3d_39eFl%LL4>evb12C+ndFV&!-f>awnv<j^pr-rM`&eOd|l5Q>1fO&<v3>i*$ps?rU-j~o5@6B#Hz>o>w{&5fycb{Lbp_xjlU4bks|!I;tgM3^sDM{Q{LPYe#v&~F1_&6v2S#x1R7b~8~^s<rX*h_*Fc{?BeOFisjgeYD<V8@GgT{g32v%a(*k08z+olAULs_ErWivH-A8nR-BQtG?*h1l-18+WdgdF7bN`g;T@NJ0+RC@chJ{kFxWIm@%{gttnO(9#$Wjz#5cR1DQxt9L;+b4&$&#yKPSeXBoDgJ)I`4w5zcmMsUS%c!J}9h}Mp`w0XS4DO*akWf#CWd1jrAa$oOND)X9)h-<5H#^)HRy~hjSsfjM9^R8Ot=U2gNn}lAEIepn!+ssCf2((Cb2XbCq5RGg@VGG;?W{xT*Ro_jV#L|-c#JR#xq;Z_3o)M^kr?KOyS`?x+SD?&0e&@c#TVcRZY+#v-fHo=UagdEr+Yf={BudmqsOI{yZ4~#Q*mZ=cqVWX(W6ENjxJ?^9Oic%b(QPe=)nCfR;qDsflEOqU=c;euWcx!-!<5MY@EDH=k9A5-G0xd{!CW6g(;r)Dyc1U!)03XH^z=hcmotfCII!^BGZK?Z)Bu>b<I&*8+MlHyb|1~Yf$O`L^VEYKGy6-Z2ARz)1=|WFy*rAkfK#y0wepqe$eQ639y{j<Rh><8R;c)CsoCLt$7C&9ppvF=7fGSO{B||Y^I+l!(w5)cq0*^6*PU%f;YAkLv|Q4O2}Ibd!wah;_N954?-VGSGUDVgYHkd7%r|0i>g}3#U&3-mRkbfbGJV_UomcfxI!>~y`w}ReO}+*8Y4)Ip$o1_FM8Sd8!Qp+K;}O~l#`eNYRi}ZS&TIE}-nA?<1;2)&O9_2#zDPfZnq8z+bGPEV&ow)V@W!j(?TL70&P#69QX+!BAhfK_A6<~a7bL>VHkq<7Npk1Y-Noge^w|GbycK-j2+Py{25HKA(__ZYn4fASzWz30!IYELGJfCTOtLK7E8Mn+nt0wPl7t}$#A15H*3IB$gyXn4v`{GkigO|Se)O!oU5ZkO9sUQ4;&8B>j~02e((4O}kuT_RhKd7Oh$>6P)6;xL5bK^<Fht)p)B>DM&g%;b$|t4p#o*)2`FMvuE;&zEO=mG!(n<s5?o4^>GtX^r_$=?>mjNAQ2HEhhM#SV!92T~?fp7)q7S}Pw`2!uy7E4g;QuE1H6`Reod>~Z=Uu2{dagCj5mMq9843UIg<4<25Y6}Flb`-Izs!P3MyRq@Et^96CIM47rlUPR8%l`k(<UPY%YLF%gfMk=JPW?Z{&y6|OPW*w4@|3)vsYGePvCr3s{sse*N_ly_nM(o4hQZymtWQP$RK)tdw#?dZ4rc!9&o?228VztfEzR-6&j$gKL{#>HJy*>UR4ZHL0Kfo&mv+>`8$Csl`3H7Z(GOaPx54V{vVr{+2+zG0xMQ#D(re{QX~!eW=9UBczkIyM*4e6n;Wg+Y89;vKds^_BtfP|}W>IT{sH{|<y0Y>N5%tJTsOI5hcM?gwoAt;*Ah|&|tH`tsh}K`>M#3ugSR@B2C9x?J%SHGz14D>zO+ZkBFMMukle(6@{O&!}H&vK@*Rof3aHWT3Gi4DNL9<Z20gfS^f&k1nICcvx_9S8mAC#Q()x2?cZX0khD>W_NXiO1I2dy3A*>?DK7p5-F4*}M~D@c_VII|2JtP}s|lWgETc!PN(Mcnk>v2wL%@&JX!ui9gOX8gE^rrv7f&gRrGq8PwYs$?}>oo>jnJ#lh+Iz#Sr%k6Fm(;@*i`NW>z|4AjoF0wC`^BX^oQtjiVM&3#Hxp93g1zLPBJ9~i4U3YkaF|ObzGcjy+27x?5P@bAbjb<+YFs*rkeC&^J<LzR2Q+JL5TD@R0jwFeI*_qg`n5XqRd|Cc*f}p&a<to!&t-?zn0YMJE4urfEr<Jhgt*%c@3PH*Ck4eIb|B-U~!`-OHi(10mssIzY{M>b770l?MAe~w{tK^TXzZGiM_)qZ=@(@2)+C;a;Sg1D$u)?PCavC;5a_)@=sjOq)e}+SZ)_6>nmH`u{1^ZEImVb>$SnFpYBuq26|CLox=|ZJZn^;Kiv=MGp^J?H{gNfNygNZ`92+G>h)KM|rUyU*^!tfY&@H+d8lXnp0d0^4o+ZXmJrZa?yjaLa?9`1zS4s>bXI}l_VSBOgy{2dH`THHGNS2r(gO}N+`#;Q4nHmrlf7xa;`SrFS`BF&bhUH<@x>+r0RO6bu+N#N2KpL>uT_N0Ghx7oNPYBzh(LDBg_@Mn{F*Z(s+zh2Tx6aI;3*kMv_2<)U^$FNf%bj@X@Qe_R7k=I_O4P#VA>m0o?78tRX4=Z8Tc(UQL5<zG#AgW<3Q50GiQU31AKmkf-=!N;RlVM2Fe+c}i@5<0-zEv(aP=%5~aJs)7-ClS8t*$B^CvTr*hb+w2??evM42355y`M>a(--r6!0s&D<01S<@(!>mb!BoMTF(YF&d|-WYP-{MHXItFsr1WwWkbrRYAw2va9XcTriesT0e1t=yBgpoTtK^7EQB1HRBSfmRvq(akGKs(;XSqxM>i1wQSk!K^d_|l<sW^u@O<8iz*B^NjZ~A|g;!fIM%%?nTbsumlBS_qd-7Jna$>7g_1bkq0@r{AVFuSk;8vU=A<q79uH39<-hipXM45>;b^#k|q?+DxjX1(<l3-$@SGo03Nc9kfl>GW#sTiGIql<djq8`3G29eSX{3qJEWwQr1Tr0YnbzdEkjiz9qrq!E~{JDY<?qW&LN#Icx`!}o>o0;jZKRH|y0X;WicOnXl6p_c4oe~#K&zj-g`pbXRA!h3u1_5FtmESgtz5i)lthpFAy0(jmeuz}*EGI6DmA+_PL9_V%cm(hRueGBzJc1WXR)f->vjF2{Fviwsdl(<<MpJgheX-X?=EUawVED1%IWuAo05_*TZRP7}ej7W+r&mkzrm<VdBRkSWO7QjPA#(yR6!DOrmXu!p2TX5bhqVBE*2W$bfRIinW?(t)PifA6-a?gt@pd7fZ_Kgs&THWd9DReZf2hW|4X4xeYoAUuo%aM6jMDEF$pOvE)1PNV@0x=;KESnv>JuVGybo4oB1!Q`9bRJ5v0j4B9--L3BQZcl@xOW8GbCVuBo7!7wR#R|>f>4L_oK6zMx=R+-$68e&IGehpFn4gFEBUmj~4RdaJJgDl6G!nH`lMp(BOzmI5e}TI<)9)A5m9asx(AK7RX2x+ahS!xDqgAVg{h97X@FQ1vXwtAb%rwa(rVdz2cP@QeGXGnpb2Tz`lw|Z7xg{q&vEM;ypdF$mU_UOk78Fu$+5eQ_Hm2!ZXy5giLnB*#q_AdhKzge4=dkp)pP*Ku!*54sUEw{sJnZ#T{9KDPR4O_%ZNPoHWCa7H}*XaAwoOwu`T-w|UHf8-P0!0L{}!86E!-oj1dYN{@eulevq_Q0i65kG&^Fc!LSIYzI=!^tA(v5-Ngz={g8*$CC)fNYqIVfEnE-v0+l+6!At>$zL`EY^V`*Aq(DEN7)cmVt&HjpDw9e4bLE*99uyBBI@67_@|qX$ssYl68*&Jx`<B=C-IZaU;obB*AmG^_TbfNL$1hj)(jMI#c9A*-suW6No2M%PMuICACV7*le(}iacn<iPpI5(6)kQ}X042kU}NQ=kCrgVhaF&t_NkuEIB10P@BAvb&##e-mEzS=g&`|PEBm!c9G_B*=2k9D`8+zPWgr%MZuG_vu4c3CxA2A^hRn{OL0H@lG$4Af*znvQQHj4ODU>J1J3EK-8pW1eXR=ezh5uP(NK^hTfAl@GFpn~)6AuTfqpgO9lBj3s|BTr|+Cq&1;gPSj7KVAxVV<y-VffrAYF60e2UD#D$EgiN2sLrj^oF+iY((dgQBT2agSpStG@nGt0RXOv`oRqo!~8~mss1`+80@-MAG+zi$^3~3p}&bRzb2Ty{x8FShH!pT`M-P_?7Vq$;$w1nVp~!O?CYH|T}E%vTh$8G2_+o}D=Fnt6VSmvB9n7kILHIRu8{X$*l6eA38~7Mka2uNSp^)F(5$d-$hv!m9Q}a)$+;k0$$*M##*#@i@Ht+eB=J0gl=~PE@|tSsLv;U125F~*25PYG_tgV+F$XIGrm{`+)%K(N$pSoaQPTOzj1G{eWQ-ADM0<~9%JTv}HP>U>U_-F<9-3kpPL5gvl>OM^_UMyIPPs*53$bvHLN-#eKkF_iw?<l*)8oJK6qu9qaYayP&^52N5__>`R-yg|<b8QtCf}S25J&LPqk^yhthgtU)=*{m;y0u}<Pfc_DZL~?0j195SuB3iVf;#S-+mAn_d2m0xKkdR{+($mTuY;w?vK(8b@<#<x3}-@X3erAcP+D)<vJa|Q8aZ_ZXb~iVw~<P?S;v1a`l2)s_nPT<zQVFx5Bqzs0q%DYflT4^s)6K9F8mUpP?qo00n3b3mC7blMbwQAe6#($=^rkG)CB+M!8B>Ya37N(O`v@;=A0SL>v|@?SUDYI0YgnDMib13wB>@qSEGkvXR*{u$ju?+Fuo6RTko>dUga$Jz8HcP;$ZH%3&s!=HDr+wPXI%`3_!YxO-TdPf3-~QBr9iX$ed8k2{id{TI6+t)~2n;$SE+CKcH;%Y?Zr$E{cDOLtE*B9ZVIXd#K;#G!Ve-FJ6%&+)WJIB#PTI>WKmwR?OBBbE7szYEPeR4pF1m_h2T+zu-gv&x~pwntUVuL}hcnuw)O*k!euU;%V1^`c!VDoGxtP!>pRJKL7T8tzM3rjDuTHfyL7Zg_<YeoDVq;%C1t1q*@^y_bpB4Z;E-j(`VMptN;05{cyu{7j(msd_quYqkh8%L7m=Wr`LE0=)eKgde0b{(Z(O4grvZxNXK03Q)PyTHkEvd?y#{X_5n^OEYW1Q=&q<{V~Qec?5Zr$!WS7uyN5A2C+5NPB=6itZvnW_w36wQ4)b#wdCVlKJh0{Kul?TZ{GVpUA5M@@-A;l6V^86EVv3$plb5H?QFy>jsI%Q!ac)-OQc=gVmV23UVt|!f}pyS-X2kog$ZcStxm4V52`KD>3dw&<%cQGVk<4fCRaTc61vkI#J6T*7)nD9KHvOHja_?$m|9ExrBGOA8E1X%0hWl8imD7I-^F?d7XHdbV1+{j4iY){VFCS6*98O4h0g#YfCO}rs~UQ~IQd(;6<|(Kz6)A>{=;)1u0-NyKZY`vwbtc&m=1}$izlk}ag^vzJ*$Xp$qt&NpYnKhTd_-izPl2*Upj7j;iG`c8U}9*$wzsV+iK!zUE3bRE{vNKU|G)1|E&iM?v({(95lBa>&u46MnZ;DOk@!y8Qg{4F21BkkuKG~_3rf8eBY7P;#9TQ*ThIOSJJI$!VV1?-&IvkFF2ogPeJEz==?lGGr12xpyMFLq#TIrp?AsGWv*6o+Gs1rHkXIW6S<$BW9&)mTfbIvtmA*D@`w~IvNzTqUi--&17Fi|p$=n2&P?H$K<H)4W#UK%ynxO09>_+rck9Y;kf<-jd3=w@bl)zbldXqsaqn?3BHE@o22>JA3Z42<137EJEPEA{uSP>U7wnMw13Z+xZ$*uH&@6W6NB8@@d_guGc@4=MrsTJ5XTXGCb2ss!C8?$NHNk|jN|jM=7OP+xup@(ng{Ls2W!5D&s+a(CGioMrl`Zph4+tm#!g^tCR*{xFpJtocaA=PgewQ4kw5*I(Gj|T($_O+GfLw_Q0ZRR|#b*U#RhEx@v+))_*DZ0TX?RUsnnXxClwISv1aBa?_N&jg$FF<4z060tZY3fyKt4tF;_V5-ros(GIHfk|1U7_^)yKj{fBv47Se2%@HHuNyD;x`3m+OI2fZ=VNyeNwQ=SIur>J-RzD(PGQib!nybYATl&I-$RXC|DD*BJ>VJ<+T9JD3wZB5jI)(`o3-FTXm*A)jG|3EC_KZ90S`?NN#FEI?BG%0=2c^!WHdj8TLflF8A1A*HpR*8EScBDh~Vx_e*Sge&xK3K-y2VSx^E*8e^$#$DSsI2XfXCfexHUX;hA<}$T*LnxJQV?y+sC!2m|VI2ejBOSBLRM};qNELHjB$qfKrY2LnSy+#qVABIF*1YP!@FsG=x(t}&Z^dc2r%(?G7j}F*&c)PS3<AI|&0cB*kV=yFoPIqM))9L@r^1Yrkd%#rUUVy`8w`rQhnb{;JS8IQt$Ig}km1O0G)|l;yUB#_rYiOPwkZb2WZ%$efeTSFkIr`aQ?=X$8}%dUmAO1s+y9$Yx2s3S93?(thXcz)fH%u%f@9V0U9<xqV<VHGFi~NWD;2Jz>q1JH9j*6TiRIO{y|TpE8(KU@--wD|u-cNG_A7L39+Y>G-&AljtrelnN4u`Fi%77?h~hSSAbi&B$5p<*Lr7A)x)#7>^^8I9^Zh;PEHMpi!_4i9O^WYdPRU|-qH}B!dUrx?4yL?YqTKOQz3z<&1fEHYO!b77ys8?c4XHb27C&o%WBbb&dV3q6>UthjTddqKNQ(;SMvFb|b09B$6kQ<Ia@uBw3||V2jeo`ZGDPL}bHK(Chz3vD*D>jmcc6tiomg}^eT;Oy#3^9#e$sty^;5{F9KS(8r>C(B)Tmq>+^+X3wS*jwom|BBQDw7y>YF5;c2+e&^~2(zdM-pTG13NpnT)R_k1{c1?tV9U0>SGX%@4%ZM_V;rqIay{T#8*{NxRd)W=_WjFTeB*oz{MkfQyc@2l6)cSaxDP6LFZQJ^H9Mi*KipMi;~0&$2mM<(J>arS9p8pSnAjA2~WLShM9hwKI>y4ZDn`josOF|9_7=Mv7;;Y(WUGn6^A*dU2>u^DmB=vNv8{nUl%H1LF0c{{Us7|0P^H&Fj6?n-9kTC#Kzt1J%Ael>L+rETO7MC^SRe;7rgjqObg!v$sdeAg}kfyHmYl@KBf6X@^NvM4`%g``UMxkj+2E4Ck?AK+>O~{_NzSly5oqYf(2XJUf??))8d_;)XkWU^C#vtt1=M#Ww;&{I&}MnVK@RAYV2))W6N7C`8KCWpk`5^Ip%)@dGXs<ed5i**86G5UT=c$8i2i-??#zk0<@4>cgveb6x3k53UQF6+%IMvBNq?D5_8dsD+=YrCfS|6<+LLWU6Ff-X>?ZVvt^{k$_TbjckO{l?mpf)re2}&2X;Rf*sL~<H&#0z|7P@6FrL6t~jNJZ5RVp2^C4+7fc{m?4fYW90msBtgx2c`-_6KzR<e37u!fE(uw>byV8p)#Q4xp5<W7JO6W`nuw62R$r894z%!Kd;dqECl<t_aa(E_{u}b$4xF2aO_K-Q*T2n-^HPWnB(Q7~40_>r|io_c{#8u$r)t>5@;MNwvJ<#4Zqkqei;2%+(ABC3bgBnuS!fx0-9Mgn@=ws;^EB1heEj<wUqCZGx>0Hd2qf;DHV|n|r{fA6B8sS;|Vf4l23C;1NGxIt8eCtF%%NOn3v0d5#<*qZg{Z6$Wz!-!s98~a*>7gzzg~b|I$tnj(YwhJzJ^u~VpX%=;OzYz$k}3wkA7Sy!r}g(%nQq8uKCzF`f!0y8G;60>hBQlp!Mh_1UALt+H6Q3-z{*U9zO##nXHNBq(dI>4k7x^~;3Yzwl7qYayb}1H09@v+6AIv9m{XmW2MAbAA>20PAdDxB*YutgbUSAfc~Mf`UIyh;@zgAD<Qy~<G<N2xffaZ#1C<H*ttxw?VVKP&x_Q0*5j?~hDDe*8-<}Ogs=Ty6eg7>)T7!)6joUdaE@GkTxmk5_{JP53B;yXtBgSvvOzT!i8U*~(w;X@|=tIIaRaECzX0K1q_?wcdqr9~Fsj@n4j&Qurj1!kM+y;?U8d7}tq=y8g+D{_Na*8q4G?v(6`3^6CZOgR~6Z>J9(K&Hm+!g<viHHB_mH;j6Sci*hLxqDNDl(*}4#lBYv>iITWl4nz>sd_<JWo3<q|`uTU5b-bu040f^&P0eTO#~Cc-P3}b{)nf0v2j#63+e2kV%5{Vvgdx$mjdg2OXdErKzlk^^2p7_V;nb3I-AeyElX;-4~NyL8U$m$<gRQVGnZbDcOJ7Z-sMw?X=45|ARfi<PmHWjyzJ)&7mJ8Wl<*2o19c6y6x{@d8TE7`KGK!Ztp4_QR)U(0N|<y<24Jj4<ZYMBoDK4XrMeyrapXuJ1RA5(%4Q{Q|IEU3_WLiYJV)Jet6Oqc5EwKm$xKk&6kkQaPg){ps%V%{5y<=yfgzkkl#Bf^pZj)z6!v(%r_aw1K8N|AJ=@rGWrf$Vw0Cf4vVqj^>|}B$wsB)tfDE-fQ}`s!0GzS{Hz;+)L-eo$2HyRyp34(2Ex2j;=7^qZl44)cVH({S?^rM7bj*Z5Mldj*cm+*Qy?o}WJQsn(rK=k6^U;PG9caYwx<s&H~)xP1(1{b1zHDY5BQGm-U1#Po|0}aNtqO`+b9A3IHxx#6PRWUsGZ?7J^%'
