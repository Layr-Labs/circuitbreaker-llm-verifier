"""Regenerate the fixed T=-7/20 depth-six four-ray gain from covariance and optimality-free rational witnesses.
No external libraries, frozen probability tables, or random samples."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt, comb, prod, lcm, factorial
from itertools import product
import json
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
T=F(-7,20);a=F(28,57);t=F(41,57)
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
    witnesses = [
        ((('00',''),('00','0')), (4706,2957,1625,642,5574,3190,5186,2473,4176,1876), -4567),
        ((('00',''),('01','00')), (4603,2867,1495,520,4268,3152,5144,1705,5228,2132,3791), -7715),
        ((('1','00'),('00','0')), (3529,1900,551,696,6107,4214,2058,3824,2577,4287,1960), -4810),
        ((('1','00'),('01','00')), (3521,1901,521,540,4574,4283,2018,3783,1749,5219,2200,3891), -5051),
        ((('01','00'),('00','0')), (2027,924,214,834,6649,853,6045,964,2089,2848,4629,2188), -4749),
        ((('01','00'),('01','00')), (1919,869,186,611,4979,770,6087,867,1913,1979,4303,2445,4276), -4807),
    ]
    words=('00','01','10','11')
    for p,q,r,u in product(words,repeat=4):
        event=(p=='00' or q=='00') and r[0]==u[0]=='0' and (r=='00' or u=='00')
        count=sum(all(w.startswith(prefix)
                      for w,prefix in zip((p,q,r,u),left+right))
                  for (left,right),_,_ in witnesses)
        need(count==int(event),'edge motif partitions')
    def lower(pattern, scales, center):
        edges = [(i,i+1,0) for i in range(4)]
        signs = [-1,-1,1,1,1]
        for endpoint,perp in enumerate(pattern,2):
            for prefix in perp:
                previous=endpoint
                for bit in prefix:
                    k=len(signs)
                    signs.append(2*int(bit)-1)
                    edges.append((previous,k,endpoint-1))
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
    limit = F(4)
    panels = 96
    total = 0
    for i in range(panels):
        left = limit * i / panels
        right = limit * (i + 1) / panels
        upper = negative_exp_upper(left + left*left / (2*C))
        total += iv(upper * (right**8 - left**8) / factorial(8))[1]
    tail_start = limit*limit / (2*C)
    gamma = 1 + tail_start + tail_start**2/2 + tail_start**3/6
    tail = negative_exp_upper(limit + tail_start) * C**4 * gamma / 105
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
ROWS='''
00 0 00 1 1211 1054 3444 4980 2794 2737 4754
00 100 00 2 542 1278 3369 5008 4540 3848 5597 2957 4962
011 00 00 1 -1124 1104 2381 5481 7491 3524 5100 2845 4855
00 0 0100 1 540 1157 3615 5139 2975 2202 3480 4451 6089
0100 0 00 1 609 1144 2659 3368 4431 6003 2873 2885 4888
011 00 0100 1 -831 1224 2488 5366 7391 3712 5281 2288 3426 4475 6102
00 100 0100 2 294 1427 3543 5179 4333 3943 5674 2359 3384 4497 6115
100 0100 00 2 254 1401 4476 3906 5670 2602 3361 4455 6037 3121 5108
011 0100 00 1 -864 1207 2417 5397 7388 2714 3318 4468 6037 3014 5007
00 10100 00 2 82 1322 3420 5087 4883 2809 3218 4571 6145 3036 5031
01011 00 00 1 -115 1067 2596 3758 2840 4953 6941 3440 4983 2750 4769
0100 0 0100 1 277 1267 2775 3306 4466 6027 3068 2308 3418 4480 6104
00 0 010100 1 479 1161 3637 5158 2998 2179 3821 3174 3074 4648 6189
010100 0 00 1 509 1151 2625 3696 3171 3093 4624 6163 2906 2911 4911
0100 100 0100 2 112 1584 2717 3295 4493 6066 4269 4006 5753 2475 3322 4527 6131
011 0100 0100 1 -673 1355 2529 5277 7283 2840 3251 4508 6067 2410 3358 4507 6119
10100 00 0100 2 46 1481 4669 2870 3193 4584 6153 3597 5264 2413 3355 4510 6122
10100 0100 00 2 -34 1454 4820 2844 3196 4588 6160 2632 3333 4476 6058 3208 5184
011 00 010100 1 -723 1230 2501 5353 7380 3735 5303 2265 3762 3190 3070 4649 6188
00 100 010100 2 306 1433 3564 5199 4309 3954 5683 2334 3716 3204 3066 4650 6188
100 010100 00 2 240 1412 4442 3927 5691 2570 3689 3185 3084 4632 6170 3152 5135
011 010100 00 1 -738 1217 2435 5371 7361 2681 3641 3195 3083 4631 6167 3045 5033
01011 00 0100 1 -162 1179 2712 3687 2862 4939 6933 3622 5159 2222 3464 4458 6093
01011 0100 00 1 -144 1163 2633 3711 2860 4935 6925 2663 3359 4438 6010 2910 4915
010101 00 00 1 1164 1067 2581 3714 3131 3426 3000 4832 3424 4965 2741 4761
010100 0 0100 1 257 1276 2740 3628 3194 3086 4628 6164 3105 2325 3408 4484 6107
0100 0 010100 1 276 1272 2789 3299 4470 6030 3092 2284 3753 3193 3069 4649 6188
0100 10100 0100 2 -58 1649 2750 3267 4516 6087 4605 2907 3169 4602 6168 2534 3292 4541 6139
100 010100 0100 2 124 1598 4236 4028 5775 2683 3618 3210 3076 4636 6172 2496 3311 4532 6134
100 0100 010100 2 146 1591 4245 4018 5763 2731 3288 4498 6069 2449 3649 3225 3061 4652 6187
011 010100 0100 1 -576 1367 2548 5251 7256 2805 3568 3221 3074 4635 6169 2432 3347 4512 6121
011 0100 010100 1 -586 1362 2543 5264 7271 2855 3243 4513 6070 2386 3688 3212 3064 4650 6187
00 10100 010100 2 85 1487 3619 5285 4643 2877 3190 4586 6153 2388 3685 3214 3064 4651 6187
010100 10100 00 2 -2 1465 2600 3659 3198 3078 4636 6172 4785 2857 3189 4592 6163 3239 5211
010101 00 0100 1 669 1179 2696 3644 3154 3417 3003 4832 3605 5139 2216 3468 4456 6092
01011 0100 0100 1 -198 1301 2755 3637 2884 4920 6916 2786 3291 4478 6039 2339 3397 4489 6109
010101 0100 00 1 725 1163 2618 3668 3153 3415 3004 4829 2652 3366 4433 6006 2900 4906
01011 00 010100 1 -124 1184 2727 3678 2864 4938 6932 3645 5180 2200 3803 3179 3073 4648 6189
01011 010100 00 1 -116 1172 2653 3696 2866 4931 6921 2630 3685 3176 3091 4626 6164 2940 4941
010100 0 010100 1 261 1281 2754 3621 3197 3085 4628 6164 3130 2301 3743 3196 3068 4649 6188
10100 010100 0100 2 -23 1663 4571 2920 3162 4606 6171 2716 3587 3224 3070 4640 6174 2556 3281 4546 6141
10100 0100 010100 2 -11 1657 4580 2915 3166 4604 6169 2764 3260 4520 6090 2508 3616 3234 3058 4652 6187
010101 0100 0100 1 426 1300 2738 3595 3178 3405 3007 4828 2774 3298 4473 6035 2331 3400 4487 6109
100 010100 010100 2 158 1605 4212 4039 5784 2697 3609 3213 3075 4637 6172 2470 3637 3228 3060 4652 6187
011 010100 010100 1 -500 1374 2562 5238 7244 2821 3560 3224 3073 4635 6169 2408 3676 3215 3064 4650 6187
010101 00 010100 1 608 1184 2710 3636 3157 3416 3003 4832 3627 5160 2194 3807 3178 3073 4648 6189
01011 010100 0100 1 -162 1312 2776 3622 2890 4916 6912 2752 3611 3202 3082 4630 6165 2359 3385 4494 6112
01011 0100 010100 1 -155 1307 2770 3629 2886 4918 6915 2801 3283 4482 6042 2316 3730 3200 3067 4649 6188
010101 010100 00 1 634 1172 2637 3653 3159 3412 3005 4828 2620 3692 3173 3092 4625 6163 2929 4932
10100 010100 010100 2 19 1671 4546 2928 3159 4608 6172 2730 3579 3227 3069 4640 6175 2530 3605 3238 3057 4652 6187
010101 010100 0100 1 396 1312 2759 3580 3185 3403 3008 4827 2741 3619 3199 3083 4629 6165 2352 3389 4493 6111
010101 0100 010100 1 411 1307 2753 3587 3181 3404 3007 4828 2789 3290 4478 6038 2308 3734 3199 3068 4649 6188
01011 010100 010100 1 -126 1319 2791 3613 2892 4914 6911 2767 3603 3205 3081 4630 6166 2336 3717 3203 3067 4650 6187
010101 010100 010100 1 391 1318 2774 3571 3188 3402 3008 4827 2755 3610 3202 3082 4630 6165 2329 3721 3202 3067 4650 6187
'''
def compute():
    budget=controls();seen={};one_side=0
    for line in ROWS.strip().splitlines():
        p,q,c,k,z,*v=line.split();key=(*sorted((p,q)),c);need(key not in seen and int(k)==budget.get(key),'disjoint cylinder multiplicity')
        seen[key]=int(k);one_side+=int(k)*lower((p,q,c),list(map(int,v)),int(z))
    need(seen==budget,'all cylinder classes represented')
    overlap=overlaps();gain=2*one_side-overlap
    need(gain*10**6>1287*S,'shifted ray-gain threshold missed')
    return {'method':'cubic halfnormal minorants with a remainder; depth-six cylinders','threshold':'-7/20','scale':S,'one_ray_lower':one_side,'overlap_upper':overlap,'gain_lower':gain,'cylinders':sum(seen.values()),'classes':len(seen)}
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
    need(N==2520 and s>=12300000000,'adjacent repair certificate')
    return F(s,10**15)
