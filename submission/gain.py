"""Regenerate the fixed T=-657/2000 depth-six four-ray gain from covariance and optimality-free rational witnesses.
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
T=F(-657,2000);a=F(28,57);t=F(41,57)
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
    C=[[cov(x,y) for y in coords] for x in coords];n=len(C)
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
    # Direct trit recursion from +1 boundary and exact finite cylinder coverage.
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
    # If minus it has negative first pin and at least one of second or third negative.
    # If nonpositive at least one of first two negative.
    need(all((msg(w)!=-1 or (w[0]=='0' and ('0' in w[1:3]))) and (msg(w)>0 or '0' in w[:2]) for w in words),'overlap cover')
    return su
def lower(pattern,sc,z,root_sign=1):
    coord=coords(pattern);K,det=precision(coord);n=len(K)
    sign=[root_sign]+[int(b)*2-1 for w in pattern for b in w];v=[F(i,10000) for i in sc];s=F(z,10000)
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
    fac=plus([iv(1-s+s*s/2-s**3/6),divi(mul(iv(1-s),var),2),divi(third,6)])
    zl=F(mean[0],S)+s
    # exp is increasing, odd Taylor is a lower polynomial for all real arguments.
    elo=sum(zl**j/factorial(j) for j in range(16))
    need(fac[0]>0 and elo>0,'positive cubic minorant')
    return iv(prod(v)*elo*F(fac[0],S)/(2**n*sqrtf(det)[1]))[0]
def overlaps():
    tot=0
    for a,b,c,d in product((2,3),(1,2),(2,3),(2,3)):
        pts=((0,0),(0,1),(0,a),(0,-b),(1,1),(1,c),(1,-1),(1,-d))
        K,det=precision(pts);v=[T*(1 if i==0 else -1)*sum(k) for i,k in enumerate(K)]
        C=T*T*sum(map(sum,K));need(C>0 and min(v)>0,'positive shifted precision')
        # Integral with radial variable r=sum v_i y_i: discard exp(-r^2/(2C)) on [0,1],
        # use e^-r<=3/8 on [1,infinity), and e^(-C/2)<=1/(1+C/2).
        J=F(1,factorial(8))+F(3,8)*C**4/105
        val=J/((1+C/2)*(2*PI[0])**4*sqrtf(det)[0]*prod(v))
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
    return {'method':'degree-three tilted halfnormal minorants; depth-six cylinders','threshold':'-657/2000','scale':S,'one_ray_lower':one_side,'overlap_upper':overlap,'gain_lower':gain,'cylinders':sum(seen.values()),'classes':len(seen)}
if __name__=='__main__':print(json.dumps(compute(),sort_keys=True))
