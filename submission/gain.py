"""Regenerate the fixed T=-143/400 depth-six four-ray gain from covariance and optimality-free rational witnesses.
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
T=F(-143,400);a=F(28,57);t=F(41,57)
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
def compute():
    budget=controls();seen={};one_side=0
    for line in ROWS.strip().splitlines():
        p,q,c,k,z,*v=line.split();key=(*sorted((p,q)),c);need(key not in seen and int(k)==budget.get(key),'disjoint cylinder multiplicity')
        seen[key]=int(k);one_side+=int(k)*lower((p,q,c),list(map(int,v)),int(z))
    need(seen==budget,'all cylinder classes represented')
    overlap=overlaps();gain=2*one_side-overlap
    need(gain*10**6>1287*S,'shifted ray-gain threshold missed')
    return {'method':'cubic halfnormal minorants with a remainder; depth-six cylinders','threshold':'-143/400','scale':S,'one_ray_lower':one_side,'overlap_upper':overlap,'gain_lower':gain,'cylinders':sum(seen.values()),'classes':len(seen)}
if __name__=='__main__':print(json.dumps(compute(),sort_keys=True))
def adjacent_floor():
    import base64,zlib
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
    w=json.loads(zlib.decompress(base64.b64decode("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",validate=True)))
    need(len(w)==30,'witness count')
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
    need(N==2520 and s>=11918000000,'adjacent repair certificate')
    return F(s,10**15)
