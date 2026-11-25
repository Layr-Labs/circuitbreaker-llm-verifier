"""Regenerate the fixed T=-7089/20000 depth-six four-ray gain from covariance and optimality-free rational witnesses.
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
T=F(-7089,20000);a=F(28,57);t=F(41,57)
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
    radius = sqrtf(second_lo)[0] if second_lo else F(0)
    exlo = sum((-radius/16)**j/factorial(j) for j in range(12))
    need(exlo>0,'remainder exponential positivity')
    for _ in range(4):
        exlo=F(iv(exlo*exlo)[0],S)
    remainder_lo = max(F(0),exlo-1+radius-radius**2/2+radius**3/6)
    fac = plus([iv(1-s+s*s/2-s**3/6), divi(mul(iv(1-s),var),2),
                divi(third,6), iv(remainder_lo)])
    zl=F(mean[0],S)+s
    # exp is increasing, odd Taylor is a lower polynomial for all real arguments.
    elo=sum(zl**j/factorial(j) for j in range(16))
    need(fac[0]>0 and elo>0,'positive cubic minorant')
    return iv(prod(v)*elo*F(fac[0],S)/(2**n*sqrtf(det)[1]))[0]
MOTIF_ROWS = """
00 1 00 - 00 0 0 -4572 4703,2958,1622,643,5615,3188,5178,2476,4178,1878
00 1 00 - 01 00 0 -7743 4599,2866,1494,521,4310,3148,5135,1706,5242,2135,3793
00 1 1 00 00 0 0 -4533 3527,1901,552,698,6150,4245,2059,3823,2581,4288,1963
00 1 1 00 01 00 0 -4719 3519,1901,523,541,4612,4313,2020,3784,1750,5235,2201,3890
00 1 01 00 00 0 0 -4237 2021,922,213,841,6704,850,6094,962,2083,2856,4634,2193
00 1 01 00 01 00 0 -4205 1913,866,185,615,5023,767,6138,864,1907,1982,4312,2449,4278
00 1 00 - 100 00 2 -9045 4616,2886,1445,535,4010,3131,5120,3170,3553,5341,2073,3730
0100 1 00 - 00 0 0 -4414 6174,4528,3480,2356,1743,657,5417,3359,5330,2531,4243,1925
00 1 0100 1 00 0 1 -4564 3706,2027,599,665,6015,1856,3437,4434,6096,3777,2550,4253,1937
0100 1 00 - 100 00 2 -8706 6169,4524,3486,2328,1529,551,3836,3300,5268,3076,3610,5399,2131,3807
00 1 00 - 0100 100 2 -8576 4700,2967,1357,570,3546,3123,5110,1853,3332,4465,6121,2935,3671,5481
0100 1 00 - 01 00 0 -7547 6167,4519,3500,2311,1583,535,4126,3315,5282,1749,5132,2194,3870
00 1 00 - 0100 01 0 -7540 4685,2950,1403,559,3826,3138,5124,1904,3327,4477,6127,1780,4930
00 1 1 00 100 00 2 -5268 3545,1923,512,548,4265,4346,2012,3776,3083,3619,5394,2125,3812
00 1 00 - 10100 00 2 -8811 4669,2938,1400,561,3793,3125,5113,3225,2719,3238,4689,6303,2152,3851
0100 1 1 00 00 0 0 -4792 5988,4313,3623,1707,587,687,5924,3790,2179,3984,2602,4315,1981
00 1 00 - 011 0100 0 -6884 4676,2941,1373,544,3605,3130,5116,1740,5291,7584,1876,3300,4476,6127
010100 1 00 - 00 0 0 -4269 6357,4769,3102,3191,3839,2329,1755,658,5375,3393,5361,2537,4251,1930
00 1 0100 1 100 00 3 -5759 3736,2060,554,526,4207,1827,3455,4419,6085,3875,3103,3587,5365,2090,3764
00 1 0100 1 01 00 1 -5107 3707,2036,566,519,4548,1833,3448,4423,6088,3842,1725,5259,2165,3843
0100 1 0100 1 00 0 1 -4613 6064,4405,3456,1828,643,651,5769,1974,3290,4517,6162,3327,2567,4275,1951
0100 1 1 00 100 00 2 -5919 5975,4303,3647,1735,545,547,4122,3904,2139,3947,3041,3645,5421,2151,3846
0100 1 00 - 0100 100 2 -8055 6201,4560,3443,2391,1432,588,3374,3297,5263,1902,3275,4501,6153,2837,3735,5549
0100 1 1 00 01 00 0 -5341 5973,4298,3653,1716,556,539,4457,3872,2145,3952,1769,5178,2226,3923
0100 1 00 - 0100 01 0 -7150 6200,4556,3457,2376,1482,576,3643,3311,5275,1955,3273,4511,6157,1828,4813
00 1 1 00 0100 100 2 -5036 3561,1940,482,586,3772,4341,1956,3701,1922,3265,4464,6099,2809,3778,5566
0100 1 00 - 10100 00 2 -8378 6184,4542,3462,2366,1480,577,3620,3297,5264,3126,2757,3219,4704,6317,2214,3932
00 1 1 00 0100 01 0 -4705 3540,1921,492,583,4095,4319,1967,3713,1984,3278,4469,6100,1851,4866
00 1 010100 1 00 0 1 -4542 3755,2063,611,658,5987,1865,3735,3160,3099,4734,6314,3688,2542,4245,1930
010100 1 00 - 100 00 2 -8511 6374,4785,3108,3184,3838,2303,1536,552,3799,3337,5299,3060,3618,5407,2139,3817
00 1 1 00 10100 00 2 -5274 3560,1937,498,577,4050,4360,1987,3742,3111,2777,3208,4647,6248,2223,3961
"""
def edge_motif_gain():
    def context(word,q):
        for b in reversed(word):
            z=q+2*int(b)-1
            q=int(z>0)-int(z<0)
        return q
    def coordinates(L,R,pairs):
        points=[(i-len(L),0) for i in range(len(L)+2+len(R))]
        signs=[2*int(b)-1 for b in L[::-1]+'11'+R]
        for endpoint,pair in enumerate(pairs):
            for side,word in enumerate(pair):
                points.extend((endpoint,(i+1)*(-1)**side) for i in range(len(word)))
                signs.extend(2*int(b)-1 for b in word)
        return points,signs
    def covariance(p,q):
        if p[0]==q[0]:d,k=abs(p[1]-q[1]),0
        else:d,k=abs(p[1])+abs(p[0]-q[0])+abs(q[1]),int(p[1]!=0)+int(q[1]!=0)
        if not d:return F(1)
        s=d-1-k
        need(s>=0,'marked path lengths')
        return a**k*t**s*(a+t/2+s*a*a/(2*t)+k*t/2)
    assignments=[];oriented=0;count=0
    for row in MOTIF_ROWS.strip().splitlines():
        L,R,A,B,C,D,mask,z,sc=row.split()
        L,R,A,B,C,D=[w.replace('-','') for w in (L,R,A,B,C,D)]
        mask,z=int(mask),int(z)
        need(0<=mask<=3 and len(L)<=6 and len(R)<=3 and max(len(A),len(B))<=6
             and max(len(C),len(D))<=5,'finite motif radius')
        need(context(L,1)==-1 and context(R,-1)>=0 and
             min(context(A,1),context(B,1))==-1 and
             max(context(C,1),context(D,1))<=0 and
             min(context(C,1),context(D,1))==-1,'motif trit implication')
        pairs=((A,B),(C,D));points,signs=coordinates(L,R,pairs)
        K,det=inverse([[covariance(p,q) for q in points] for p in points])
        value=lower_signed(K,det,signs,[int(x) for x in sc.split(',')],z)
        multiplicity=0
        for first in range(2 if mask&1 else 1):
            for second in range(2 if mask&2 else 1):
                reflected=(pairs[0][::-1] if first else pairs[0],
                           pairs[1][::-1] if second else pairs[1])
                pp,ss=coordinates(L,R,reflected)
                assignments.append(dict(zip(pp,ss)));multiplicity+=1
        oriented+=multiplicity*value;count+=1
    for i,p in enumerate(assignments):
        for q in assignments[:i]:
            need(any(q.get(k,v)==-v for k,v in p.items()),'disjoint motif cylinders')
    return {'oriented_lower':oriented,'gain_lower':2*oriented,'scale':S,
            'cases':count,'cylinders':len(assignments)}

def negative_exp_upper(x):
    need(x >= 0, "nonnegative exponential argument")
    # Reciprocal of the degree-12 positive Taylor sum, using integer
    # Horner evaluation before the single final rational division.
    x/=16
    p,q=x.numerator,x.denominator
    den_power=1
    polynomial=1
    coefficient=1
    for k in range(12,0,-1):
        den_power*=q
        coefficient*=k
        polynomial=polynomial*p+coefficient*den_power
    out=F(coefficient*den_power,polynomial)
    for _ in range(4):out=F(iv(out*out)[1],S)
    return out


def radial_overlap_upper(C):
    beta=C-T*T
    need(beta>0,'positive overlap cone bound')
    limit = F(4)
    panels = 192
    total = 0
    for i in range(panels):
        left = limit * i / panels
        right = limit * (i + 1) / panels
        upper = negative_exp_upper(left + left*left / (2*beta))
        mass=(right**8-left**8)/factorial(8)
        first=((right**9-left**9)/9-left*(right**8-left**8)/8)/(factorial(7)*(right-left))
        end=negative_exp_upper((1+left/beta)*(right-left))
        total += iv(upper*(mass-(1-end)*first))[1]
    tail_start = limit*limit / (2*beta)
    gamma = 1 + tail_start + tail_start**2/2 + tail_start**3/6
    tail = negative_exp_upper(limit + tail_start) * beta**4 * gamma / 105
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
00 0 00 1 574 1136 3489 5144 2877 2787 4909
00 100 00 2 118 1354 3404 5151 4469 3961 5754 2986 5132
011 00 00 1 468 1185 2326 5244 7095 3605 5296 2938 5029
00 0 0100 1 260 1229 3664 5295 3043 2285 3403 4498 6164
0100 0 00 1 355 1219 2720 3312 4474 6084 2949 2937 5042
011 00 0100 1 101 1303 2442 5203 7089 3796 5476 2386 3363 4545 6203
00 100 0100 2 66 1493 3579 5318 4307 4043 5825 2418 3332 4550 6204
100 0100 00 2 41 1469 4431 4004 5817 2648 3329 4513 6141 3146 5272
011 0100 00 1 182 1289 2375 5217 7071 2786 3268 4536 6145 3112 5190
00 10100 00 2 -67 1394 3453 5219 4819 2887 3205 4612 6233 3065 5196
01011 00 00 1 739 1134 2649 3816 2788 4836 6706 3495 5127 2804 4885
0100 0 0100 1 118 1335 2829 3269 4522 6128 3131 2386 3357 4536 6194
00 0 010100 1 553 1227 3679 5297 3059 2246 3840 3251 3042 4638 6211
010100 0 00 1 591 1221 2667 3721 3248 3059 4616 6187 2976 2958 5052
0100 100 0100 2 8 1643 2757 3273 4553 6171 4258 4097 5896 2527 3281 4584 6228
011 0100 0100 1 -10 1435 2494 5154 7038 2909 3212 4583 6184 2508 3302 4581 6227
10100 00 0100 2 -29 1545 4657 2937 3180 4628 6243 3633 5397 2466 3312 4573 6224
10100 0100 00 2 -95 1518 4778 2904 3187 4635 6254 2669 3310 4545 6174 3231 5346
011 00 010100 1 286 1304 2457 5186 7076 3815 5487 2351 3766 3263 3054 4674 6254
00 100 010100 2 352 1493 3596 5327 4274 4042 5817 2381 3749 3279 3043 4659 6235
100 010100 00 2 287 1474 4391 4014 5823 2602 3732 3258 3062 4648 6226 3174 5287
011 010100 00 1 340 1295 2396 5193 7049 2739 3650 3270 3064 4656 6233 3140 5211
01011 00 0100 1 402 1239 2760 3746 2824 4863 6753 3680 5302 2297 3405 4507 6172
01011 0100 00 1 461 1228 2686 3770 2821 4853 6739 2722 3315 4486 6097 2968 5039
010101 00 00 1 1313 1133 2617 3711 3172 3430 3075 4753 3450 5088 2775 4875
010100 0 0100 1 336 1340 2778 3661 3268 3064 4643 6220 3165 2399 3346 4531 6183
0100 0 010100 1 377 1336 2838 3259 4514 6115 3151 2348 3774 3268 3046 4658 6235
0100 10100 0100 2 -46 1708 2783 3248 4581 6198 4607 2963 3160 4651 6263 2580 3256 4605 6244
100 010100 0100 2 223 1654 4217 4111 5906 2712 3659 3277 3060 4663 6241 2545 3269 4580 6218
100 0100 010100 2 257 1647 4225 4100 5892 2768 3262 4546 6159 2491 3680 3295 3044 4673 6251
011 010100 0100 1 171 1445 2516 5127 7011 2861 3582 3293 3061 4674 6252 2526 3290 4579 6220
011 0100 010100 1 185 1439 2510 5136 7022 2921 3203 4579 6177 2472 3691 3283 3052 4685 6265
00 10100 010100 2 226 1547 3652 5410 4625 2939 3174 4619 6228 2430 3715 3286 3047 4675 6253
010100 10100 00 2 136 1525 2626 3702 3268 3063 4668 6249 4739 2913 3177 4629 6244 3260 5366
010101 00 0100 1 892 1238 2731 3660 3202 3436 3082 4783 3639 5265 2285 3387 4472 6126
01011 0100 0100 1 268 1358 2802 3693 2854 4867 6768 2840 3259 4535 6139 2412 3346 4546 6200
010101 0100 00 1 968 1226 2657 3681 3199 3431 3083 4776 2701 3302 4446 6048 2940 5028
01011 00 010100 1 561 1242 2772 3731 2828 4858 6750 3702 5317 2265 3817 3244 3052 4656 6232
01011 010100 00 1 586 1234 2702 3747 2830 4848 6737 2678 3702 3242 3068 4635 6210 2998 5064
010100 0 010100 1 542 1342 2789 3646 3264 3061 4636 6209 3188 2362 3758 3265 3045 4653 6227
10100 010100 0100 2 160 1720 4567 2973 3151 4647 6255 2739 3627 3289 3057 4676 6255 2599 3244 4602 6235
10100 0100 010100 2 182 1713 4575 2967 3154 4643 6251 2795 3238 4575 6188 2545 3644 3302 3046 4685 6263
010101 0100 0100 1 687 1356 2775 3618 3227 3429 3086 4794 2821 3254 4503 6100 2399 3335 4516 6163
100 010100 010100 2 430 1659 4186 4116 5904 2724 3644 3275 3057 4655 6231 2510 3662 3293 3042 4668 6242
011 010100 010100 1 341 1449 2532 5110 6998 2874 3570 3292 3058 4668 6244 2491 3676 3283 3050 4680 6257
010101 00 010100 1 1036 1240 2744 3648 3200 3431 3076 4783 3662 5280 2252 3811 3239 3037 4618 6185
01011 010100 0100 1 412 1368 2821 3671 2861 4859 6762 2796 3633 3265 3066 4656 6232 2431 3335 4546 6195
01011 0100 010100 1 435 1363 2815 3678 2857 4860 6762 2854 3250 4533 6133 2379 3743 3264 3051 4670 6248
010101 010100 00 1 1069 1233 2675 3662 3202 3426 3078 4777 2656 3701 3234 3053 4597 6162 2971 5052
10100 010100 010100 2 356 1726 4537 2978 3145 4641 6246 2752 3613 3287 3055 4669 6245 2565 3627 3301 3044 4679 6254
010101 010100 0100 1 810 1366 2795 3597 3230 3423 3082 4792 2775 3636 3259 3056 4625 6194 2418 3325 4518 6161
010101 0100 010100 1 839 1360 2789 3604 3226 3424 3081 4792 2834 3246 4502 6096 2366 3741 3260 3041 4639 6209
01011 010100 010100 1 559 1373 2835 3656 2864 4854 6759 2810 3620 3264 3064 4652 6226 2399 3727 3264 3050 4667 6243
010101 010100 010100 1 941 1371 2809 3585 3230 3419 3078 4791 2790 3623 3258 3055 4623 6190 2386 3725 3260 3040 4638 6207
'''
def compute():
    budget=controls();seen={};one_side=0
    for line in ROWS.strip().splitlines():
        p,q,c,k,z,*v=line.split();key=(*sorted((p,q)),c);need(key not in seen and int(k)==budget.get(key),'disjoint cylinder multiplicity')
        seen[key]=int(k);one_side+=int(k)*lower((p,q,c),list(map(int,v)),int(z))
    need(seen==budget,'all cylinder classes represented')
    overlap=overlaps();gain=2*one_side-overlap
    need(gain*10**6>1287*S,'shifted ray-gain threshold missed')
    return {'method':'cubic halfnormal minorants with a remainder; depth-six cylinders','threshold':'-7089/20000','scale':S,'one_ray_lower':one_side,'overlap_upper':overlap,'gain_lower':gain,'cylinders':sum(seen.values()),'classes':len(seen)}
if __name__=='__main__':print(json.dumps(compute(),sort_keys=True))
