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
    rows=[]
    for a,b,c,d in product((2,3),(1,2),(2,3),(2,3)):
        pts=((0,0),(0,1),(0,a),(0,-b),(1,1),(1,c),(1,-1),(1,-d))
        rows.append(pts)
        K,det=precision(pts);v=[T*(1 if i==0 else -1)*sum(k) for i,k in enumerate(K)]
        C=T*T*sum(map(sum,K));need(C>0 and min(v)>0,'positive shifted precision')
        v[0]-=T
        val = radial_overlap_upper(C,-T/v[0])/((2*PI[0])**4*sqrtf(det)[0]*prod(v))
        tot+=2*iv(val)[1]
    edges=bytes.fromhex('efdf7fce6e5d4c3bab9b2a198a08bf')
    part=list(range(16))
    credit=0
    for edge in edges:
        i,j=edge>>4,edge&15
        x,y=part[i],part[j]
        need(x!=y,'overlap tree acyclic')
        part=[x if label==y else label for label in part]
        pts=tuple(sorted(set(rows[i])|set(rows[j]),key=lambda p:(p!=(0,0),p)))
        n=len(pts)
        need(n==9,'intersection pin count')
        K,det=precision(pts)
        credit+=2*lower_motif(K,det,[1]+[-1]*(n-1),[5000]*n,0)
    need(len(set(part))==1 and 0<=credit<tot,'overlap tree credit')
    return tot-credit
from lzma import decompress as unzip
ROWS=literal_eval(unzip(b85decode('{Wp48S^xk9=GL@E0stWa761SMbT8$j;1V|mhFt&?mkGlRFTeYt4FEz|Z#$U(4{=E6hCv~g`Nb4qT(pco2;r3KP~6TD2#?5ZUGR%2)bLW73a|6ne9wZd%PB7mYjzT+{N(tP15^}co4BOU9AX`#Pd-DP8i4Cm$`Eqpibj-5@O#cTW->;pX*|3wTUHCReZZqEu-Lc?OJJ!^A<06|1re|E%WUpb0hl85RkR}>0G-(;^5ERpf3Fka;H~mqJwTYuU=gwwmKbhLwxOnS3(vhjdLSQV@{&DNoOjf42cVW3%hwGxb~_C$ATm-RAUd5)wC~EKOeZ#$Oy4hoBz2&B?jiV7zZA^=5!1^Nyf^)(rg(o47$vs#yu=#+W@K{m#6SAunw075#tK48XHTFd-McPQp*rSFT*ICRJ_JF3LZzgg9++ZVobQDkc4#QaeA?*jnX?}|!#nRjTPh$GVP|7*__FnOy!Wg&=NP6ZjUMoc%<sT&5%BMFWrGS%2j^NG_#@JR7i95;JdPr8lJParUo#b8^bA1SDWt~7E`rha-Eny4#w-vUKs^d7-C{0D5C!%s;D4iw7)sI^l;hnju^IwF^7)*d=EM&?X?&;h&@cx}b*1hRltetl-_PMV?!<XiO(0H7i7WLcSxZ*NY=EyJI+WgRx%W3QM%rrI60n$sDDD$lr4GC81_3YzatD_D+n>_;XByeTquN&X{JXPluSa<-V>>P>!k49A{Jm_4K*}EOSLwows6`c)t&{RF3s1C#B}G7~lg@pVi07FTd{wD?Dgz~g>K(-py$4Rb|3A)l7u(Sm3-L^k_(KtX#=g%qcpZ*jBVx(wmdwwRJ%^AgEsXOBFyX7HHRJYn@Im3@ZF_i#%?m4q*2Kjgvv(uSNX#t$)c@VVDnz(}7_e3HdKT(MjI;bU-J3<ytXmtZ;Kq*CdNh@vA_@PuTj;ic;|}GD1wVX6{$SIi4*LIK&u%>D?|JH7gJ8B~Nfns+(f_Qq2m;%`cVeL1o*7LoUhS<w7IYMKWRhz0Q^K7>t+McYAz315l7D4R)ym5ahbpKPZP+eSc_3t}%t1M(qT*qg-e*fe-Mzt3IF7HRb(dy)#KC0gm}~GU13c`dLuahRBJ@K)Vdy9^2I&mQx;zwAep!8&cCw1EN_Z*8%W=U=xoqkm*C`AdUf#*wpU7UU)x>xY>^h##<0TkWr#cddY7S3~w~g<AzXvBdqpvgQyDZ25rN1+rbp}tl(MK47Y(j{<Ftu1+UUUW^pfr&;ej0Y1(IKI88@Ow4La0{>B}cRRH|_velgtvh&V68g8hgE~FemyqF+laLm}nVdi<g9-{BAs)nltRu&F8x+{rV2mCe>H%?&hp%6{xLc%`%dLCr(Z1oa#1NA4oL8$)rhw{m|_Bv9WSbAi{SnIgwJ#+0Oe9oBvl&5ufeJ+HYd`VC_H9Rue`15#g(7-#5yMkd@--1#p5!5pfT->#0A2AG@`qpv5eSn@t|GWzGjeFt{r|eJafWyER{oan*)H0W~xMX-bJ{|D&gt4`<xLV20jSX7xsNrIFa{lsYiDWbc%ta+@OC1`1s8@sO@K(2eHXeQbbv>==-rLOzgejCxqGG_=p^T2~K!3}U?`{w5v-Qj3=;5@p}_|5Dgk9;dJnLf$q|Hc9N;yx<VmL04A@6d3TICAk^Z&UDfz;ze>g<RpQV)kfEfzbn8s5S28uX6eS7&{uNrZIR@z1XDl&fYzbM>@NP33QO>-(mUQ=IRYn2_Ery0SgpgpKEfKhPG|5R7Rd&`0kY~ft^#|u8{YIsokxIWq%(o;bHV-A5kA|YW;xqkCIrIKZApIcGoUx?@ojy;NnAT%tA3`p1d4vK47fAhG>4Mo1XH2v0N@eT4Mz$B<`S^$M8MdSR3OAxu;CZt`yYI&VCqlpTWD9nUfPMrYQQrrmhMkggg8us5r3oa0lL5{bdC*sW|fD&9TFyE0yu-xzI+aRy_a(q>b5>`k76bT>NY=Qt0oVsKkMGD0=l5Vy~9=%dkFUw8oKsfFvV?NdY&P_j-fPuRV5%wWp4sU+6PR*|9@2dzn#Gvck1k7xOffRK=l~n3Agi?Q~Pvr0}f;tm5r|otV>EQ&!N>PnPB;*Dx@FPM`L=<WUGHJ-Xx9)On;!V0T^z2E4$k#8LoWm&9U_!eMmKM#No%GafSFYPKrAK00000wov8{LZ5hV00E*6xFi4oz$lZ5vBYQl0ssI200dcD')).decode('ascii'))
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
def centered_moments(B, ell):
    """Directed variance and third centered moment of the tilted Gaussian log density."""
    n = len(ell)
    di=[B[i][i] for i in range(n)];rs=[plus(r) for r in B]
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
    return var, third

def cubic_factor(var, third, s):
    second = add(var, iv(s*s))
    second_lo = max(F(0), F(second[0], S))
    second_hi = F(second[1], S)
    remainder_lo = second_lo**2 / (6 * (4 + sqrtf(second_hi)[1]))
    # Jensen on the exact integral with |Y| replaced by sqrt(E[Y**2]).
    if second_lo:
        b = F(isqrt(second_lo.numerator * 2**24 // second_lo.denominator), 2**12)
        r = sum((-b)**j / factorial(j) for j in range(4,16))
        if b:
            kappa = b*b / (24*(1+b))
            lam = b**3 / (240*kappa)
            def tangent(j):
                a = F(j,8)
                u = 5*lam*a**4
                return u*(1-4*a/5)-u*u/4
            K = kappa*min(1,max(tangent(j) for j in range(9)))
            r += K*max(0,second_lo - 2*b*s + b*b)
        remainder_lo = max(remainder_lo, r)
    return plus([iv(1-s+s*s/2-s**3/6), divi(mul(iv(1-s),var),2),
                divi(third,6), iv(remainder_lo)])

def unpack_adjacent(payload, quantum=125, wide=False):
    from lzma import decompress
    from struct import unpack
    data=decompress(b85decode(payload))
    rows=[]
    while data:
        need(len(data)>=9 and 0 not in data[:9],'repair prefixes')
        words=[bin(x)[3:] for x in data[:9]]
        n=sum(map(len,words))
        start=9+n
        count=n+n*(n+1)//2
        end=start+(2 if wide else 1)*count
        need(end<=len(data),'repair length')
        body=unpack('<'+'H'*count,data[start:end]) if wide else data[start:end]
        entries=[(i,j,quantum*x) for i in range(n) for j in range(i+1)
                 if (x:=body[n+i*(i+1)//2+j])]
        rows.append((words[0],words[1:],list(data[9:start]),[quantum*x for x in body[:n]],entries))
        data=data[end:]
    return rows


def adjacent_floor(first=0, last=184, *, bridge=False):
    L=6
    def elow(x):
        k=0
        while abs(x)>F(1,2):x/=2;k+=1
        y=sum(x**j/factorial(j) for j in range(16));v=y.numerator*S//y.denominator
        need(v>0,'positive exp lower')
        for _ in range(k):v=v*v//S
        return F(v,S)
    def source(L):
        edges=[(0,1,0),(1,2,0)]
        for ray,anchor in enumerate((0,2,0,0,1,1,2,2)):
            tag=0 if ray<2 else anchor+1
            for node in range(3+L*ray,3+L*(ray+1)):
                edges.append((anchor,node,tag))
                anchor=node
        return marked_gram(3+8*L,edges)
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
        pins=list(map(int,core));n=len(pins)
        for x in product(*(range(a,z+1) for a,z in map(bounds,p))):
            l=[x[0]];r=[x[1]]
            for u,v in zip(pins[:-1],pins[:0:-1]):
                l.append(out(l[-1],u));r.append(out(r[-1],v))
            e=[cbit(a,b,d) for a,b,d in zip(l,reversed(r),pins)]
            need(all(a[1] and b[0] for a,b in zip(e,e[1:])),'positive adjacent C edges')
            d=[]
            for j in range(n):
                c=sum(e[j]);pp=sum(cbit(x[2+2*j],x[3+2*j],pins[j]))
                d.append(c+mbit(pins[j],c,pp))
            g=[4-a-b for a,b in zip(d,d[1:])]
            need((g[0]>=0 and g[-1]>=0 and set(g[1:-1])=={-1}) if n>3 else min(g)>=0,'repair gains')
            N+=1
        return N
    if bridge:
        rows = """1011 1 00 0 00 . . . 11 . 00
1101 00 1 . 00 . 11 . . 0 00
1111 00 00 . . 0 00 11 . 00 .
1111 00 00 00 . 11 . 0 00 . .
0111 11 00 . . 0 00 11 . 00 .
1110 00 11 00 . 11 . 0 00 . .
10101 1 1 00 0 . . 11 . . . 0 00
11101 00 1 . . 00 0 11 . . . 00 0
10111 1 00 0 00 . . 11 . 00 0 . .
01101 11 1 . . 0 00 11 . . . 0 00
10110 1 11 0 00 . . 11 . 0 00 . .
110101 00 1 00 . 11 . . . 11 . . . 0 00
101011 1 00 0 00 . . 11 . . . 11 . 00 .
111011 00 00 . . 00 0 11 . . . 11 . 00 .
110111 00 00 00 . 11 . . . 11 . 00 0 . .
011011 11 00 . . 0 00 11 . . . 11 . 00 .
110110 00 11 00 . 11 . . . 11 . 0 00 . ."""
        seen = []
        mass = 0
        for line in rows.splitlines():
            core,*rays = line.split()
            rays = ['' if w=='.' else w for w in rays]
            k = len(core)
            need(len(rays)==2+2*k,'run ray count')
            for old_core,old_rays in seen:
                need(core!=old_core or any(not (x.startswith(y) or y.startswith(x))
                     for x,y in zip(rays,old_rays)), 'disjoint same-core runs')
            seen.append((core,rays))
            literal(core,rays)
            edges = [(i,i+1,0) for i in range(k-1)]
            parents = [0,k-1]+[j for j in range(k) for _ in range(2)]
            for i,(anchor,word) in enumerate(zip(parents,rays)):
                tag = 0 if i<2 else anchor+1
                for _ in word:
                    edges.append((anchor,k,tag))
                    anchor=k
                    k+=1
            signs = [2*int(b)-1 for b in core+''.join(rays)]
            K,det = inverse(marked_gram(k,edges))
            mass += lower_motif(K,det,signs,[7000]*k,0)
        return F(mass,S)
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
        B = [[iv((int(i == j) - sum(v * KU[k][j] for k, v in columns[i])) / 2)
              for j in range(n)] for i in range(n)]
        ell = [iv(-sum(z[k] * KU[k][i] for k in range(n))) for i in range(n)]
        var, third = centered_moments(B, ell)
        lower = max([elow(E)]+[elow(E+c)*F(cubic_factor(var,third,c)[0],S)
                               for c in map(F, ('-1','-.75','-.5','-.25','0','.25'))])
        return prod(L0[i][i] for i in range(n))*lower/(2**n*sqrtf(det)[1])
    w=unpack_adjacent('{Wp48S^xk9=GL@E0stWa8~^|S5YJf5;2LoTc3l7mflY9qok~<kThWBsq9||~I-sSuj<YxQ9Azw$O^{3cmC1bbs!Jr};bS5HDqV|+;J$Zxg=F{NQop#{x;(*M=IPe#Yfd;6>NoK1NxVLn$_g~j&C=I<(WZv4S85q^7v9iax?hqE-=OymxN;ef`Q?dP;RFo2d3u!`dZdviU+0YtzNm<~!tkN)HD5*NTo_%{;LR-rk3!$7H!hzeob2*7sA7F28~6_8b@+=mSbK(#+y@7q*Wx~YkFX-h5=?aI4xp4+`NEAH<J{f!vuEO1noj7<A3yiIoe}R)a5N5m_5Kc?M@%2l@QPPkT=8lTqz}8f@1h*8B#Y>R8uWPZ`Mz;%A-_hj$J0uJq)4!;?_5tf<o*7GW7YN6ZONcsik{D}JG11%2<4PI<saC+apqDYy(%r8CN`lbqk#hj_coxQolZYaR1M*}2APJsX{w;};@1RuvRXu~BXwf>l{WCK-)zHS8;!dk*h+q%(!oC@=u?4nB0d^leXxsyA1i)Fw^IMBK*L5HJ+%sGaF+(=$33Vn4WR5i{R@kg!v0voF1V^*49&n1h1Zg{SGl(kkj-whVUw2hlW-(8XnSnV6Ev;u(ftgRs+Ijpdu<ud)urFmSZhJE!a1w(%&7gv*%8B|PYo(TGFAp0jzUZ@5koM`O)I+A0&sdYc31!zP)S|npAG>_#`E7owZ9_F4+B)tz~}5+e^T{M<YMYy3U@;z8=nBnkz@wXfoRLjZ7-RyUw{*NZboSJYqM;}cr(j`iQn^8^1&xS*FZ1P2Eebc9_oLG_2Eh);zXtr3^?4}3Qa%5S-WMu<>xCs4!T-QF_<u*%2GVGT^%8OT=9(0H?{ISw?_B)c(B%?L}XsTNeHnptj>A6GE5Y0@Ak&-DcU9@2f7Xz1#-Ab$t;|<Hu?inR&*=shUj3hduIb8@&?hUr4dOfQ?=PkK9tI@7v@oR7_Auz=8>g3y;!wKlmQ!|$|qbyilJHD)tQ6N@I0<`Qq~QSmU`5?y?rS&u2a7rzz1_f0MR)6ni)Ir%}tKGK5cgUo>B;vRh;k5b{`<O!xi8P>US$=pQ7YZ1t2m~3VnW2i=%lrpM>8Ba^hO3mw@K)d#7>j8ba{ZE*SMzCXUEx=$O2dr1eAB92XzzUT(|dix4U%BlUrH5Bw23saVECFbuNuS^|N@Wn_`MtoL4K4~zmc*Z>z1@@iU<!ra}1244*0%#_@jbvoc8K?UE-uSHl5KL3N_HG3TX5j*pepSa%Xk;I0K2dYPUf!<{vx$T&5XSA)Rj^$YQnY2sZ@@t%Rc%i@`4{sENYz~c2fvLdx_1d93_T}}{(`;Km_^TWVH!K0meUw+n!q?6N9_pNPnX6H`rS7WPew9K&28+LYW+#?c_6+tRZlQisiSlIIa6RlV|4TfKm@If`eE#ZJr}x(a;7q(Bvpq+PgNC+hoxkN#>GYmg=Ark~Giv377e`(J>Lnh)3V7jeWkR<pZ5(yrZW2lB226hUZVtXzz!=;z-JTSOa)7f89QqVXGsk>+jgtH+=`~P;XY30qe5lmVa!By+$e^F;a;0NJN6eA<qCoi>R}^>7x?>}Z8An#^f;vNENdSFk{TQ?Y(UJuw!`1Ewh6`R|f|;~ZO>|B=hip>i{j^9-P8O8Wi0p>9W@_BE4&F=ztSv;R_d91ABu{ll;1}k};K3fM%|SAtKh&3!_^FX1L2?f^Po95Zq1HsUKpMFb@i-fIZ-Ys+R#uA^#=WVHQk8#Vj~4j~@2-&(5&z6hHQJ|5$T(EqF%rMF4YvM}KvOAqWA#Y}3&U=QK}Tgm)r7>tTcQ08S}#+&x6})$dV07u#8n<mn&lCPL1lJK4JH|;6oS<TYh)T^6|0VUh*+B>BTjl6sfV09%i+q&^<SBnr)3gE9^&R|IOlq%+ny+{f4+fkX>@;TRS9p}qUBbBdeZ9b0%XsUdKk}bZ2z`;fFuoI6(*j|Ug(THSHo?i0guQQBmA}LN4pvi7Z<bP%E$hjut7rvBKp*Pz2LO}k!X%r4vxpa{L4ZFj`AiBK&Kh6l9t*2ib%EmxM94yZQDTYxnRygKel&b5;#=k?Ot#EnLT{Elt2+}ZrYwNJsIBw?dJ`{!8UvoQ1$YmY;b(s$H?tYf4rY+761SMzJ%~Rgfpcp00EK>@-zSdL^Nh(vBYQl0ssI200dcD',1,True)
    w += unpack_adjacent('{Wp48S^xk9=GL@E0stWa8~^|S5YJf5;2B{Em|Xw|fk_n@Hymh}<V5|hP4P0A(i!4_(Nm=uw{@FVR6daNE=O{}nmnEO!yZZjleQI$`fN#h?Zw?7ib$*yf+)4nT1c@YHQsf~P3NwFbs`e6kY9b7XM>tTnD58W=8+2z&W~&N0PZ_L7MnK&>FsOZ7Gcz2GM?9?8!S_Bs2h*&<7f6Md!?c(6uCbSV0$IFiJDi57cLwUCP&bOl(nS#JyQM<#7?AVRUr5_S4dtQ5L8J6<1dz+o>()+v+iXrbxTNw<FkmwfG$RA`{Asc<Zy`^ovcJFpXlbIoj{E;@jPVs4G?E4*~DY4(9#!@$mILXJfIm9L5g(X&R@<TbX|&vL@@X$kD@x}F*cU~5NeYkW#f<uHjLEVuL<2P%gF07h3~;0K49$&mWb17+fU7s&DT>$xpH01FCyQ^zc-~(`aYG?FX`9$%A`XjnT@ipIaYecf6E)o<0_q#Ctbr1Yf~{p(T?u*aX05A8$(*(2MR)|UXVryHPT(_+}BgQQQocG;lhy%e)eK5d3R$j$BQPqXj${n|C0H(=*D$H$c5uV3u=c{?l-*g?>OAhw_ny5<PPLm)$>U^jB=c_!C?@o5Qocg^ORc7t!|W|4n7ND(NsT6HCGsY^}jD6N8b0VSa55gzb0GCm2Q-mK7CCy{o3_mdRH2KtnW79`h>gSmHB*KL9aI+9Ef{y{sq07kAc*7-K+n)EPTu5>&-UkAUTzw)X5sDE%I77heGDn0qw5$%dboi5@N>CF?%|`Pf>w5rC}}B(D$dIRB=2#ULGP(?*=&tE3)Vsgw)Sy>raQB!flaa$~BK3w&kZ43qH6P1VVG5ysii_8;Mc+zLZ`B9O*|TErC~+l}e_d68}Agkknj1ocpi!G%@nd*Uu1F20p2!nNIlV&VO+5Z}4qbHVARtr0q0q2w85^4@8%^oOE1F?m|l=F`f3=X#o~ID-_d@azFaTt2RT;LHb~nbWo$cF=;y(%>1fgXR&NRR`{^L#Pidm8k4DUDzy7IpX)n0F*9syy*(;{{5~?4kMxtU0ZTr@PkseNSNc~wsUhnoq%0(E<~w+#brSHuv`a<%1J2V38IR&;7DU4ifJrN9dbW~bDNR`Z?qjGILp!|msb;HPu2EWrrl*=>Svyvju;W%lvBm5#&Pu`E#q~`(1h5y_V~}*(gn1e*<3<umNrGRaVhpkZ143Lq(d^GmI;tZM%waZ2#R0iSFE7kjq80^p5b(ld51o6iCH<qa2*gpx*dfUu*OK>?u^X=jTMjqiI&3KZ5~UEiw{pa-`|YV%7THq38ip!0EzH#nD#<H>R*J?Q*{V0?{7ki%!|g@xQoS%R)_hT;yAE|&U7eSc#yhaDlZEZZY%?RDo2D%qNhgj3LF)+o+2J<W*I*!d;##l-&9nlR&t>t3N|FOP?YCu=#km#TMVMwtrc`?S7F(ud+8_EWE|ehCUTWN4#8nCs0YboMdqStCYbu>PFUo6L&k(yl*y(1&8X95~(H{-twWVk^p$I&Anfzgh;}{Yke0YHXebO8s%1d6Z(7Iz&9t|x^_n5nNm@(d%RK_TYt4RKqTxiJrd~<KSy3>24yC*KXp%9=?PQp<jEZ+;i>Z`4P$)hW0mv9%RC9%c<><JI=GgXcC6_c)x5m4QwV`N`5v1{P%J6B>YC#&8N!yoI>zEk7;7nOOOlKVc9Kk7wc?rIj|%@f0-Y$;1ODr$tdiBt&8?~3AyQfa?K@`NR%8inb(C&ybmZbE&YI^#}`$C}QxMk8Q*+6_#FR1-17aCxijG2WvUD5M6`ZsLI!VVV3c6{<Z^4ol&GYV;F%i?O+I?3ElyLJs{|V8pCEGG%)^cjd$$kV*D+b&41F(G9`qDyGlq^TDp?S}j^lua#B6Vuz4M{+rmkzt9w!Ol5O&<NoaQu$^08?!`!%8|tWvhuM^0^A1u3>y3?_TV7xIJYGlL+rPf6{WOzETzb&Rqj0ug0f~7`j%}er{RIdLdpmSL<Ej08H}Q$qHw33E<AAHij603Zi6Q<Qs>u1NDwf_HV%8=p^YE15;N!<>!m?Fa|A3(UImhcM0w`A{GQIXQVJnVD17_8s9|_OgcL&NXIMV!3-W+6!F)No<Ov(E;p;)ec)fk}LRTto~x`zR`^gCsF|NEYxZ0DFZq(ataBE;aCV-|FnoHz|Ll1YmMFO($SWtNHfi^^Brd>brCI+F5Q6CpDu9831a^dt;<!HWA=n@o5&w~^+1+(K8ps0G?rd5j=h5sbS8F0bLDN4gu=i3+Q(b>G8n>}4K3ZHS0!<VeN?i1Lo9AY4}B$tNUzRhwz2hoQt8u1MowIUv)>MB+*fc@~FwF?!hw7lfe$t%&$BLm?LH1+=42^f83s%8jgbE}t9FUNfUmGjxjIPaaKt5hO7Este3mq@_3&y*Tg_UZsPz|K;L6EuM*9u}4ofPLU|Y&gYa9=2NaByMe5+-BN?2d|^L${dR5y23F=f;382zeKk|d!blvyvDH!+842!d+r=@R-@||Zxj23j)8m^V6c97;R5?~(?7nK2P_&`#tH1NlFLh(~Cq?_MHiXuHy3an*RfXuMtL(Jn`H~%%@6x{JoS!$ui4in(`-#WE+%y|pHh6b%3J-mu1Ad^i(rQ<x0hPAvda<+G9wm~$WL&b77UlX@EZ4i)YG~>D8=J1wEkx4#zX34?R}=lTFs@k*4n-ZGGoQAW(Kwa4)JP9OU=n%S`|%5VBl4@n#d@IoR$F(-_X$ir;nwUC#cba>_xx4KIp5>oq55pJeWyF(H7mILUFPk{MY=?Mfn(pZh_UWuEr41euKB{Fp{HzZcuZhkfd`jZYg;Ei>;@sPnY;i1^fo<+B2#sl00FcS;xYgL!Gjk>vBYQl0ssI200dcD')
    raw = inflate(b85decode('c-nPVNs`+z5EIGfzG<qM(|p5E@&Dfvpd}|6Njdn^Bzr*vkf(?fXYdPNTF)tN#?~@rIWa~|6<<u)Tua^NIGt3nIPlTrYB%YjB-f^MR>WGdzHgz%oQg9|kOUbVaz`({<TQ{w)HX#BHg2(I95r=gx02gjM;!Nq|J)bccc>}e7=)!FjxdU%Z6n@FV~9u2nhm)PN;~JyyPT)>8ka#X<2D}?!X}8Oz~y>gj<w)@POXm>P2`mQXX?d*nDUCr;_x?20POkzJ8RKlV(#ZO?UtS6SXmasW2Brq?t3i4aLOh>Z+GX~(r*0tgLqfS7atfxBtFEHOT|Q;8$$d=40|QCI4|vtzl3qXlr^3T;<#(L(VzJdZ;C;;0V^qO#5^~baL|+jAs(HMiDM)?xq*>xFMxta@o@oy6`tY*$m7i=USi}D4T@ofxD<R&6ENsxNL9J3p32~|-Xb4-Bsl-UXUbcia6>bY1qmNy=vJb2n<=0wZa4SXB-_8m<^9Q@iYLSKQJ0)*xBE3RUaf2~fRU9%cO@Teux-_0V)lQ@wy3mIIZQMM=G#(NFiEy6`h%LJrXKs3HZKhgiAs~}eQtL>@Bsv>i~E8iTtGvn^_T{>K_#LiDzP`Dq7qFVsuyF{Z38Mh@JHsID63I;54BfgZAfx9N-tn(k~?ZFHJ)pYUB@HiuWKOQ*yrZMalmGgXUTJ(8V@xYJYcHhL1ssuz}9*0ZD4QRHLT$vyN#o>01I#CFt-JByRN}Sw(}qwXV;T1l}J?07LEc~yMrTs+;h=x`e{!-HywaLTt)+3mSTM&p8C6$my@XJ=MJ^`W{1l1rT+S_TOP!QpMlr|#7-o(K@E-D<;R2RuKx`t?URVwKq6n0n!4C}m^xCy7%;)<XjnnZW$$r#lP2bH;XHvuJ)Ja8aXUHfu7Omgec>B6akPJqq|xD@j_`(<ME%$T^hr0VdGd3~g!_Cllot>Nz(bYYgWcE%)Nc9c;fWBk;4a>$h@r99qqgN2!1>{;-;7%Ns$f5kO(!<G4x9I9JG+o76rF{BVNF-8#`-=IIji3oXad<YEmQ^H>8+-1uBNivJR{YjAD^d%)c5~R^h0(puXD5;5l`o8oa8Wx{Yxk5RBL-HhpxE}Klr7#te;E413Wd~R=*P0O}yptYdm-02)(wag|GLC*5`{nY2^jbFM09c&MOry>if%xqJMg-yK3h*VuwfK'))
    extra = 0
    while raw:
        core,*prefixes = [bin(mask)[3:] for mask in raw[:9]]
        n = len(core)+sum(map(len,prefixes))
        need(len(core)==3 and len(prefixes)==8 and len(raw)>=9+n,'repair row packing')
        w.append((core,prefixes,list(raw[9:9+n]),[0]*n,[(j,j,5000) for j in range(n)]))
        raw=raw[9+n:]
        extra+=1
    need(extra==105 and len(w)==184 and 0<=first<last<=184,'witness counts')
    C=source(L);N=0;s=0
    for i,(core,p,perm,b,entries) in enumerate(w[first:last],first):
        need(len(core)==3 and len(p)==8 and all(set(x)<=set('01') and len(x)<=L for x in [core]+p),'prefix shape')
        N+=literal(core,p)
        need(all(conflict((core,p),(row[0],row[1])) for row in w[:i]),'disjointness')
        ii,ss=coords(core,p);n=len(ii);need(sorted(perm)==list(range(n)),'permutation')
        ii=[ii[j] for j in perm];sp=[ss[j] for j in perm];cp=[[C[x][y] for y in ii] for x in ii];m=[[0]*n for _ in range(n)]
        for a,k,v in entries:
            need(0<=k<=a<n and type(v) is int and v>0 and m[a][k]==0,'sparse entry');m[a][k]=v
        # Additional charts optimize the full lower triangle.
        from math import pi, sqrt
        K, det = inverse(cp)
        P = [[float(v)*sp[j]*sp[k] for k,v in enumerate(row)] for j,row in enumerate(K)]
        alpha = sqrt(2/pi)
        scale = 1-alpha*alpha
        L1 = [[v/10000 for v in row] for row in m]
        b1 = [v/10000 for v in b]
        mean = [v+alpha*sum(row) for v,row in zip(b1,L1)]
        g = [sum(v*u for v,u in zip(row,mean))+sp[j]*float(T*sum(K[j])) for j,row in enumerate(P)]
        for _ in range(32):
            for j,row in enumerate(P):
                diag = row[j]
                for k in range(j+1):
                    if i<79 and not m[j][k]:
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
    return F(s,10**15), N


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
    eta = [F(x,d) for x in (0, -200000000, 134069570, -10553063)]
    scores = [F(x,d) for x in (6552915, 0, 3213780, 0, 6609393, -287519204, 113217634, -118955334, -118939821, 101388890, 97369999, 784647, -103509, -98225777, -301273352, 302500000, -283135214, 297500000, 302500000, 297500000, -302500000, 302500000, 297500000, 300544506, 10994053, -174552605, 299577087, 130466210, 96120091, 302500000, -126998628, 124394232, -65747138, -6638697, 3602754, 4212825, -32034545, 1855374, 4565334, 155012462, -6120577, 12917435, -33139204, -95799849, -121987140, 20680322, 26227498, 302500000, -56855004, -95497910, 302500000, -99029221, -300068184, -7037821, 11627800, 6609393, 3000)]

    checks = prove_combined(eta,scores)
    saving = sum(scores[i] for i in (0,3,4,56))
    cap = F(cap_numerator,10**20)
    coefficient = 1+1/(2*cap)+saving/3
    return {'checks': checks, 'local_reserve': str(scores[-1]), 'saving': str(saving),
            'semantic_increment': str(saving/3), 'coefficient': str(coefficient)}
