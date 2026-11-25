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
    return {'method':'cubic halfnormal minorants with a remainder; depth-six and depth-eight cylinders','threshold':'-291/800','scale':S,'one_ray_lower':one_side,'overlap_upper':overlap,'gain_lower':gain,'cylinders':sum(seen.values()),'classes':len(seen)}
if __name__=='__main__':print(json.dumps(compute(),sort_keys=True))
def unpack_adjacent(payload):
    from lzma import decompress
    data = decompress(b85decode(payload))
    while data:
        need(len(data)>=9 and all(data[:9]), 'packed repair prefixes')
        words = [bin(t)[3:] for t in data[:9]]
        n = sum(map(len,words))
        start = 9+2*n
        count = n*(n+1)//2
        end = start+count
        need(len(data)>=end,'packed repair length')
        yield words[0],words[1:],list(data[9:9+n]),list(data[9+n:start]),data[start:end]
        data = data[end:]

@lru_cache(None)
def halfnormal_bins():
    terms = [F(-1,2)**j/factorial(j) for j in range(15)]
    integral = [v/(2*j+1) for j,v in enumerate(terms)]
    decay = iv(sum(terms[:-1]))[0], iv(sum(terms))[1]
    p = mul(mu, (iv(sum(integral[:-1]))[0], iv(sum(integral))[1]))
    q = sub(one,p)
    first = mul(mu,sub(sub(one,decay),p))
    second = neg(mul(mu,decay))
    need(p[0]>0 and q[0]>0,'halfnormal partition')
    reciprocals = [(S*S//w[1],up(S*S,w[0])) for w in (p,q)]
    return (p,q), reciprocals, first, second

def conditioned_jensen(columns, KU, totals, z, elow):
    weights, reciprocal, first, second = halfnormal_bins()
    mu0 = mu[0]/S
    estimate1, estimate2 = first[0]/S, second[0]/S
    def priority(j):
        diagonal = sum(float(v)*float(KU[i][j]) for i,v in columns[j])
        cross = sum(float(v)*float(totals[i]) for i,v in columns[j])-diagonal
        linear = sum(float(a)*float(row[j]) for a,row in zip(z,KU))
        return abs((1-diagonal)/2*estimate2-(linear+mu0*cross)*estimate1)
    chosen = sorted(range(len(columns)),key=priority,reverse=True)[:6]
    def sides(x):
        return mul(x,reciprocal[0]),mul(neg(x),reciprocal[1])
    mean_shift = sides(first)
    shifts, pairs = [], []
    for rank,j in enumerate(chosen):
        diagonal = sum(v*KU[i][j] for i,v in columns[j])
        cross = sum(v*totals[i] for i,v in columns[j])-diagonal
        linear = sum(a*row[j] for a,row in zip(z,KU))
        affine = sub(iv(-linear),mul(mu,iv(cross)))
        shifts.append(sides(add(mul(iv((1-diagonal)/2),second),mul(affine,first))))
        for earlier,k in enumerate(chosen[:rank]):
            pairs.append((earlier,rank,iv(-sum(v*KU[i][j] for i,v in columns[k]))))
    result = zero
    for sides_chosen in product((0,1),repeat=len(chosen)):
        probability = one
        exponent = zero
        for fields,b in zip(shifts,sides_chosen):
            probability = mul(probability,weights[b])
            exponent = add(exponent,fields[b])
        for i,j,c in pairs:
            exponent = add(exponent,mul(c,mul(mean_shift[sides_chosen[i]],mean_shift[sides_chosen[j]])))
        result = add(result,mul(probability,iv(elow(F(exponent[0],S)))))
    return F(max(S,result[0]),S)

def adjacent_floor():
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
    def lower_shift(sign,mat,shift,T,K,det):
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
        factor=conditioned_jensen(columns,KU,totals,z,elow)
        return prod(L0[i][i] for i in range(n))*elow(E)*factor/(2**n*sqrtf(det)[1])
    w=list(unpack_adjacent('{Wp48S^xk9=GL@E0stWa8~^|S5YJf5;3we<>RkW_flY9qok~<kThWBsq9||~I-QepDA&OvhK$vt;T|GKwOV-}8wX7pxgM$6yauY?wM0=Yi9Mk<Cz$#wG|a}Xg*S;iT7k6Y=oQ5?2FLR~yE%yFtp9k`9>wooSnR`KC=N2IN<&Zv9Qwg8Ex(?EiQBPyDh{S3!P;BAKIu!i{tSORT5oxwU00;v8ek!-nDY)P@I({Zhh7#@x83671<%x)Lsy65KNZPPzHg(#!+d)Cp9H&mRv9;<ihHm^OK0~%+i$Oi?7<WwFwhp&Q;8L=ob^8K)C*{F+0HcN7Nx=<N^X}+79kX5Y8f~YGVwCa68v(ZtE^1$CL<J@It|UUeg@vxO&P)B93uhFj5sZ^b7sZ?;}MU+=X^umK{AKPFbu^UzKi{$bnWn@x1VyK{n@|&FPg$r`{D_lpo@R`P>xq-O%{%E&3(d$1Dx<+Q~=%*khPoa4%BEcLXW}Onh{W_g>>a6^<l9T+s>i@N8wC>&gGa@WIaPi0xwoN5&cQ`4n^RFv$UTxaJ`^QS4&6urPZ@tmdrBMB8bYYM~O&<kP0aE`NtDF)*FjxOQ{0aBC^B;+X0D8PkHk;emCe;JUKge#nv#B)}$h@4wLiM8qY4lj~ty#PiN_~-RpBFs!n;H;wgD{!UcrxHg8ozu1l^Bw=EUxCmpTk9=<S2=%V85XZc*he~W&miB~uAN>r;gzEsB0=G3~k&2tJEQONWdp3!Vs-8oi0pOwHI=3v1^l+}_ws$PA%_b+zNaV{arH9uqODVit`09ViLKS@k9>MwNnK0l#;U3*nFp!FV;Lfno@N6@cr6+J)ETN8HTrVCjbPR_yeOaPSoIrPZTj!f0M4NjgctuPXcXwf@Ea?8Ha)H{kDWP#BNl%Yxlzn+thE0y*bKI02x|IUB1Pw;G6-3sygq)~M;Yp^K~C65#4;z=M?pJ@-!O`vJXuvw3aX8!xn2l%wT=W?zDi*$ChPw}i`#FBg03o_jJc;)p+2>?N;R9o=+exRk@!lEXZyX8s{<oSiTs}^qw6dF0-9Nv5p_Y&jEh@WGbz`2e|z{QFHGLz#f>Kf9~^Eo=?n~aotrE(;Ct9T`ZQLNr0F7|>Efy{SHY2eCqj!N6lFb_@Lg*_5SW)VeBGtQK?5n4$Lfz{ltK?DA-m24S$b<(#zls;>rNaIz#&k<T#MHb3Zq=7=06I4e}yJY6p?19`C@OAXij+CL6TRY?MD^H{XYK?B_g1NgD#J0JHF^%I%Us;eazaeH2xq~N;KI!1UNRl+hseFRT;w|VQ@a6xYMa(DYa-ReL%D|z<zL^XMAq7jHu0VuhD<;$f@7lbZsdFYENhhN#kQH1WF)y_veL!m<q%hWp0nItPRfL7_E?W>UA#2qq`LRwh47dGgo`B2~;6&FLrqsa$oc%)+#kXT-iO-l<hv)|Iov-E|PNF>v164i*J-ABT=%h3ALNw#@>`Uf99s#nqTUAd_=uqQ#(tC7b8jWnl%*dk#yPoN@>2I?6mgzOk5o!_CwrR6$L%a+J6=)Oc{D>u}v{uOk#=wrgp<k~<ut)*C+2bA_uXm1N#G17y+ls_+Xbq|W7ormNWFiDdjM=4AemHhvs6A4J7KNDLfjA=dL`Q+~BuED-&)&YZRTC8x2JJQo-e&P4zbmmW8q28@=Sc(*tD=LppwCr2KAxpxVg0uPs9J9kR_UOM+z4}%sPrt>;V&(8w5F`1MrcF<mct##+v!ZRl=6r9<^8~;G918nBbea%Cn~-EvmAG4-&kRjdKA=S%GFl%#ZacB6BD$d=F2%ezq{^^c$pMi{4oI@1Hmq>r*v(L<Cg50svp=nLLT<Bxg&Q_P0DeG#k{KifqWC#y<guwfH)v7n!I(RmLDuICi+s_Q<NypWZ;N855!>Sbxt4uvx}NhI>%OwC;c^<yFo2LeH4~_b6Q+UFP3Q<*D`loEj<TSWa=$g_^j0rU(UI(YpnGPFO#sG`J;PSV+S8yG}4>NU3m@k-*r7ZL;^dckUY2_iYp2iWfK<OpH|6(DR^L0Zlg^m@Yn>UnP?I}_#Z>;ZIbdu`MWprpXC+vxr16@_qb+ml+?Cq=T+o3?e3any<M5ZIX6&vHEhko26&>FddfY_dMf+S*a<;GX5;RYib$}xjYHMP&MqbTO^@njFE>_@Qa7mJk7YS>Td9N?m4N;qSL(pum|55Zkt&V*Y0TFe<L+#$Bgb@@rw5*6Sz{uxXna>8KT#&~OlGshH5Yg(=2K^h8c*goa;+*7Rew<Qh?lmH?s|G`{Z;Nk{G%q~(4V5LKPez)nc61xr}+7vRTH}@jgw)6>l<op3nELcCBj*z3$O&XZY9ghI_43qjwZiBXMPEf2oR}JyZzjYynHiXk?AXQ1KSw-V-?rb%RV1H=I4DbA%^m4B7DiuVN0Z18>e&}s!AN&6zN-H!nlm{0&Lj5qHFaNT{?hi_@_5Jvv1jrAu&R>%xJ!Hyd!;(9mdLd1woz<XJAp{yP>bW5G8OO8D#oddy`C{vGfNo8h4ga#vyIBFS?jni@JlmoJwL;e{Ho3{Z+cNH3h4oAQn!zEHd24>8-qGvSspzbrM`R(Zpr-*lLPl;3d@@lFo%7h%Df$N|3!xBFNd|b4Ace*6IrRMpqjL${hm%Sa6@b04{XEXU9UX0n0h;SVvQ<er&2Ej~2-j-Wz3(#e|+UFS0Ls#AEb4GB5V3y>jpsw%bZ>z9nc`{>kOuMz9$aP2gDSK^$$;We$zPTcu4u$5INEn&bBoSiKv{SKwjxIvCyrw!(5q9eDhu5D0?pd``EYni8X3RKR60gng;>PWz*mKb85EjYz5v1N(uvEhXajM<*JNN^Zq_cl@L9uw`*Jx6Er#V(P?)pHv}yVLidP0MwY);;yIM2-+@11^6h_2rOSf6VD3LT~(&qU#NL9MC)b&2I3g*1UuYKLzIhvkx~nqBHSG(kvmOQ(5MU9p3{ASjBTl%aJhE!99jUT!)XOnv=F9qy%w)PkmcOl1oa@R5Y?nS6>o}QP(Y)8#7qWNV32TD%$COcoLT5o*x~PxDPJtY%xv@U1XXf9Bv@qY`Gy^V*8)dwL_SZ~J+~fsIhb?B;J%R=<c?%_doh>!Vq3Q+iQ<pSaasiURRaBN5E}j9z}|IgVdKeM{>!e3h7j`8ZcZ|&*h6)W1my(kp6wK8+>keJpkf2_YY9FW^2Od%rrnZtA0-Y#YS_X*UJjT46UR>*lxY|Y@O|dWYmF>$STXZssJ&REnK6om<ftb__6R5JlBfPlmXm8ftb?(IZ!<Ng@;oGJJ&!0{d-5n#brOaFpxK#P=?>&_j>YPn3hr;wQD29@pSoU<q|Zt=5wS17`u*!5pOa)lKKS<7llab1QGUW4I3&oe(a8E@+%Ytzb;Hy}lX*>jfFuuTrw~!aoYDO(@I{JVer!jRdbYVmx{UjgC$NV5g)XDe)b7biu?kk5B2DW$mF6M^et(WUE=Qqdqwac+DCOiVJlH0X<$UPjSex_m2P>xw)t@5+ffp4nGA=hjvpaB_JC=u**J1Oy9*7Dxb%7l0b=sBj+OPeDV|swgcLn)L8C>3*+iXv!6lZ&w<0y?x{UU$va9H~uhi4Ky00000C!pI&V#nPF00D*;;!gkobONm7vBYQl0ssI200dcD'))
    need(len(w)==79,'witness count')
    C=source(L);N=0;s=0
    for i,(core,p,perm,offset,triangle) in enumerate(w):
        need(len(core)==3 and len(p)==8 and all(set(x)<=set('01') and len(x)<=L for x in [core]+p),'prefix shape')
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
        v=lower_shift(sp,m,b,T,K,det);s+=v.numerator*10**15//v.denominator
    need(N==5976 and s>=31265426500,'adjacent repair certificate')
    return F(s,10**15)
