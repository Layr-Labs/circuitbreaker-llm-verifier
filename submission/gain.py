"""Regenerate the fixed T=-29/80 depth-six four-ray gain from covariance and optimality-free rational witnesses.
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
T=F(-29,80);a=F(28,57);t=F(41,57)
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
    remainder_lo = second_lo**2 / (6 * (4 + sqrtf(second_hi)[1]))
    fac = plus([iv(1-s+s*s/2-s**3/6), divi(mul(iv(1-s),var),2),
                divi(third,6), iv(remainder_lo)])
    zl=F(mean[0],S)+s
    # exp is increasing, odd Taylor is a lower polynomial for all real arguments.
    elo=sum(zl**j/factorial(j) for j in range(16))
    need(fac[0]>0 and elo>0,'positive cubic minorant')
    return iv(prod(v)*elo*F(fac[0],S)/(2**n*sqrtf(det)[1]))[0]
def edge_motif_gain():
    # Each pair lists ordered prefixes on the two perpendicular rays.
    # The six cylinder patterns are disjoint within one main orientation.
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
        # Main pins 00111. The off-line paths leave the first and
        # second positive anchors along their separate paired incidences.
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
    # The two orientations of the main line are disjoint.
    return {'oriented_lower': oriented, 'gain_lower': 2*oriented, 'scale': S, 'cases':len(witnesses)}

def negative_exp_upper(x):
    need(x >= 0, "nonnegative exponential argument")
    # Reciprocal of the degree-12 positive Taylor sum, using integer
    # Horner evaluation before the single final rational division.
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
    return {'method':'cubic halfnormal minorants with a remainder; depth-six cylinders','threshold':'-29/80','scale':S,'one_ray_lower':one_side,'overlap_upper':overlap,'gain_lower':gain,'cylinders':sum(seen.values()),'classes':len(seen)}
if __name__=='__main__':print(json.dumps(compute(),sort_keys=True))

# Untrusted signed differences of two knot planes. replay.py verifies all comparisons.
KNOT_B85 = {
    21: (
        b'{Wp48S^xk9=GL@E0stWa761SMbT8$j;A00OxLp8xKYSoc`U1P1$C4l?gf0c#*}1&bbSXh$Md^^>=@KWqUJAm3oyfK<LW+&ywgwf0'
        b'#m>A5Gp6TuB83k0A%!$E`lIDPiOHi706TGB31<SJ?O}T^J>kawE=qkc!*zH8FwTa4@`e@s>!q)ls{Ech3C3d{6jtLeY^yrm1>^WN'
        b'gwbD#0R;Yx{p78=|Kg7F=geElwT<BPgs6Us@F_a}je5q|IT+51pxt-{V2w6>H|mQEFH*MH!Gst{d2oL(;*8)dkFg2i8K;qj53yp+'
        b'ewQPCw2D&!89xS)Mfrk`J^sgO(=)y1?9dj&<H-J)5i;Mdnaj!%OH(#Z#cNUp>C<T!J)0cVk-(P%@k|Kknzh5f-a>J)e+9a13Yqz}'
        b'=@Yl{T1?&98{KaN4?NfqJgMd6>l-MWf{GL&;S_%TU5#fO_c}<5I_ON*7L{%;JWY}3KT{p!o&2GCI~AALia}Gr<!u<zzx5ehr%Wl5'
        b'6M8@T!oQ)(ZF|ozunn|ZPKQZIcL+R~t@-RF7kxZEt>JwwGXyCd6uouF{d-Qh)j$sK0@$75{Zb~0xDh48QtiII&Q9uexyLyNln{R^'
        b'<@yQY+zV`=Y!0MuK-l15R!xEtgavNv1Imw5Ano3v!Q`^!jqFjjSu*#2;PFi0Yrd_~DSc0;dR#gPCIGoiX>-=GuwQt(^bjkY1<*4v'
        b'EKAA#a`CJ(^E|3R5F(KwWLW=|f`dkeXk`Ag4GQWU^$2j=co*fSIB{J=y5;H*h{uQ>Fh6HEs*@0Ud(oe2*H5tKl;-Af+SQ}@>OHIs'
        b'l-nE}$TpLU#N$ubLe*h^{_10ii)z?aJy3cB(n$1sPfZIf2yCJcBOnco*(hEG4B@brO6^5psA%0F&1glb!K0I3VOO%MCs%r+COmrI'
        b'O(b=<OKgedpY}7HP|~r37VJ49Qza+S;e2Im$it~$Gvwl=n{4-8#3_`Ctvfd>h09}z@cOWlp&9)-MW`Z(oQPj<Rh|ITa9%wnW&<c7'
        b'gadt~x5tQKZ2gS3#5dq~cZ?&wW5>gMkhtOtR)m~yb(@n4e_Uu_H(;h+HP!J<H*H-Ih?LuQxbCbC<NeWc+dtCWc$Bg^M-Q<N5bK&N'
        b'4U*23bUHE_d=xP;8R3r-GL(H9x`!#B?t~)qDBAnaUsSN>Kr*MoT}rJQVjdB<$NrPIe@SxiD}o*ONP^BRhs<pd^`I3U`PGZS4Il*_'
        b'qo!ewEhFEPtQ@QM9oTp0{crpxvPDC_0ZfQ4Dbt|gzqz!rTbH162=6Eh3smz(C7gxRS!W^W6zd>*HJV3wuM2{*?vpk^nL|dK(fWRt'
        b'%T7>nef&PHhlx$#%~B!%pzoa(;Jp?EXP5!!Ds92YZ&lBRWfcAI-+*;1+0yLA7#a3-CKsk_Cz=G4XTTs>7G=X{*%TkPa=bOO&#xOI'
        b'@-gl7tHFJ^ar&nqWsiBxuugY{$<_=)0Ba7DwOve^f+VfU0fNNM#6F!^yP}c5B_4Gh07N*Z5Zc@t5kzzJPT>2bzlwaanEdQYGa&Kk'
        b'yv!lJJjMdd9+7i|cNJS=${Fe67(I-JtJ#@yQp|v0ZQ!N@Et9kx)6_`vc_`*f#u+M%xU#e~Gq#EWWh)0N&Z+BTu;p-Ob-0@k85;87'
        b'3Qy4nn=p^XqVd5kS#H3XY~vL3)p3SNgKFMLh6%?HK)YF(+Vpt4#e@YlCfg7)&M}o{z}1Vy*jC0_W?ybGZ=%wOL9sT^)Y_?ID^Nr}'
        b'svi*{Xja^l!g?Q!z@Gw@)tW2KDVA8(d{+2_V2vj{)`O+XzwqhM;^#o5fTtO;N0;2@eC?O{oC$Os2HYnE*AvxnJiMUwpLO*Cx{<h;'
        b'm#~@w<W(rWmPh2c(L(r_X*M93o?NlPNnvMinX<lW-Rpk8A2UkODr{h#pCT17YT?0CoX)rI2@pBBS`R>svHkw*+|Z=%|7hlQth3}g'
        b'*AY5w#%8}y(+-FP9$p`U*&~3q_>_;!6OZoKxjGR88w!Mq<0Us=$w-ezlG4U<Tqe_p?PBAA&NvR=(}lypr(MOdZ=gTEBNla7z}NEg'
        b'WCy2nQYj_!H*Z5J7~R)}8EkPf=%8}NclARi!?Lw)nv8Os-JU7elPf)BK<1dNxWCl?6mzzJoEr5p_@kuWsTXu@0BCiAAm9sWQ-16E'
        b'l@zpAm*jDa)nViGiJi6Yaxq`0WDQ;<jsZ4aDB1moea-Pl0I=n2$=ii#xw4B4tt~R}*1xN~&ExkI6YtU~`7w#GZ(q;vrP59uYJA@E'
        b'P?jE1ef07rN-pnwaSsCW02ww4vl;l|O;9LA5G6K;#Rs)+lFzB#l`0QpaC?6vmRTP?caIRl`)IVMQm117wvnFE?3X;IcwAs?#vg$t'
        b'Wa*EjCHuG*$lQyb%GqmBG6|b9cV+@08sU;Qv_IBABLMJ0&41zeADqWehETjSmhg!uq>NgPjGJYY{{nK<g(&FeD5dyD@~<b-d@|t9'
        b'253bgIAUTRgr_kuF~diVlXI;&aq143$vZQS_HL>1k71a|o$WCg{LKTfX^!77qxKw(%OeHi`^;bb#UrTSeAd`U5b}k9Py<8*F1cdy'
        b'S9I1KV7|O+(hMU4we7Yk(-E61B&s?{+g{9>a#qKkE>?8AN|MhT8i<v}z*je+^*43qt?oKs89V9Mt46}*U&}~5zd%kf4o=wEwGz5b'
        b'!~7wNYc^*|tk(^P;txnn07?`(fs6V*zfoFK9QZj}W(&fFHuUIuX|)XwZtxLx=nbtXVQeah)xHfAwv~u}+Bs}!vycXkg_!XGW$)Xo'
        b'lF&|kL=Qo2f2=Dp;(pDBKkD2_;0WZ<Z+drmXVj*o8;DlhOV5}o+sPWn#&+y!aF&eVh{(h7V?@u_k{D#C<f4Bp3?w7wa%?6zOU+qk'
        b'Qzb^ye<2+yqtiJLeK2>XwAq(yMk{uA!parzu>6~s?TQ7+wtR}9G-h<5Zfq!AvqM~tNRLRf^$9~&He914N8JT>0GabYQ#G-v$XV;A'
        b'{+fHqoO(VIOPNtV`u;A#i~WY)vLZfa$_XrUVCk1qT^1~D+kKubs{^PaAcmF>o)ru+BhB&;PWPsu5LyH75{WPEgoF=TiisHKZS5y+'
        b')KdDu-~F)c>DHNU5Tr0)wm;LMP`@e_zJa}*3ulrZSByP7fi$^XWx2D_@5o;hJYRglojJrOMKy3j_m<!B9klyQ!sfFKK7}t>7R)ob'
        b'o1$<hl8z}ZM!z)_CD^2PePVa~x9LBEuqU61SJ+w>-|gkkr!J%AW%q6eDD?e7cT<Hjks+{%YcnblZ=Fcg1;@2sYY=rM<{rgaF~QZ7'
        b'ry`TUpV5d00{WoGlDfOVoou^@`uD<jAOu7<|GsTQL^u;`z47MKf&9Y?=pvca8!hndlq06Lt3Kn@cT#pMJwNuJ*-@~^Mc%FI2!MI)'
        b'mTb-zNr5p(I7K_WG3PWSY%tNj+Rn@OZ3M@H!ZQr@k$Ac*AB1&xKUK5RFv@s{B69g}#&W*p7j)nI<J}85E?SA9JAjSpH}*szv0)Y7'
        b'=nkHTZt0R{L^?`71Y^uIE(Lt6_{mkV@9D5fB{w6{1ZQwXM-T<97}>7tUro>AURMeFX3E!n>*k86xN_WxV(K+X=!jpHbZ0y&F#BAu'
        b'Cp#gAtW&9PRL%H&nCkFwN*GfSu{lgFZnkuiCvr;bgMep@C0eR(bl&Q>+$OLbc@jwP<@kXvVrOb&`dT!>t1?c#VE|dSASlgaU|8)c'
        b'Hp9`~88ok>os0@}bhUbzVE~@ykgT#O*6_ppEx3K(>CsOsbGaoc)D13;PpB41YS!w;Wzm{B<w`Qja8jH;%+{iZT)pB{az-(+B|d!p'
        b'Rxlaop^!x-1iONtoy3c0NK<jo#_`q;pC?H~K=}W(_JlAgU2+8Dw~8UpUn;zSZ>H33rB-=ZzqJIefsYX;+2%+??cE9MpbIh>5W(Jq'
        b'&kS{c6Lku%xZocgUr9fAg}O!KiryF|cI1@rS?~UbQc0Ns{iJ6;I3I(J18ej4AyW0%D4m|*<gVkwwXT&{K~{2`)r04uvSC!`i#MbT'
        b'g1FgevP0@rc@L~U+=D77FpftWqzvqbB8LNFnXlt>fI8d?DNHpJ@XE*^G*|PeclMFYi0aVh7%g*xE-#0sKKjaxj;q0jCpfM0Bi#wj'
        b'<LU|>xxzGG5QYYHWR~XxQA%bz6GEZ(wPG2@Znt~hSk$@bN7VZ<2Heec^+zrSo5Bm3Vd{B>&oOf?sRm7%?1#6YMP8=_0?l^PYzeIq'
        b'ENJU6v5(f224NP9zQ;5?wCO$yCC8Yq@Jg=?O?Du0uNt834srD#w{BZZQASbZ=C6p$hf2z1y{1UEzorodgfJL*8+1Ci-|{qj_vfBk'
        b'zR0UnQ`-2J0U%OP?yGJNu!uBLb2tAouF1RhGGS6MUL>?#P|{e?8dBY%yEx8l<3(nXc=0gz$<w=YkRzUsY@8yt*MH5&HmW!4Zl<ez'
        b'OXa@7|1FSW>#M^{7SLPd^#nX_@2&H?ZAkGVd1~LSl%F~_o}PQ3+<=&Nw&nO{E6c=q0#dygWusp6aT)|pWS5t!>&6DTAI5s<g%zE4'
        b')Eg-?v`EO4v=l-UYaYHFUO!;lgI_EB-XXL4wl%iqL8gEFiF>xLSlqWI#>707OwhnRNcFJywn$D-JoDyZ^Z`UwP97qq+zws~w%u+{'
        b'sxQMFiuhseU;cN(uSLNP5xS`@$^t}GCXjPbf|}~ti6a_;td9DepEcmFx{h?OG?23G`2{bccL7_3f3=sYY@>b<*i7Vr!(S>Ck#>Xo'
        b's(6vcIBg8-vF1V{K)L-XFGUR-v@m7*lff$zvUQHY2*vrTivc9>(YLwebpc}%1?@sH7{8X6+cH0S8ZxhCxDrP^nYkahXCW5AY81L='
        b'MI31%#pXBNQcjmCy|}>#N{&I_yTR<*BGmGs)gY)k=N6F2AEHai#Jdp@2%;EisbojmQ>DYBD;@ENY5H$@+h+euiM+JN3$eF8?O};G'
        b'bm%IfOB@t`8~gnwqu}B+;SDmX`}^su@re7IYoEK*+4y0smFL+m#m&D3uX|`SMNl>ToZ>i`_Ry4s7GK_8D^#(^GOR!aoHgwKa?UCf'
        b'2cfsya;i^uu#y$@Ejz$G2hQc=5!9`&MqjcBIe}`_f5hmbStX?frr$C%_x*e#V&De_Dbd<ObQA!cj8}mFh~1&*gXaHAbI23s-@zuc'
        b'd#=wY`(D`=##o((TiHCPMT&(sW`z@p5|@&?p1I?!Y{=GlnEfXjZDqCeCg%u_B?|J<#rqiE%>+0_kh(H;MB(Y~W9$4xlVrKa&CA}{'
        b'u<H(I{AIn#*yN#@s?0rJ;ec%14=VI_?;b4_q4J7c1S`B8zGM)i2Z<s*A`d}94qZ6~h{)w-wg#A^A$>fJb~<O!s;luzXE`9=AL7O2'
        b';|JemH?lf!;crS@7ee;%ZwI?TgU!-ca?!m%@GiD4-0uMxY{owf(L_>*>w@%8F+&Kh12G%KVJUVd*lM#Gksx5QIj{@G?PKv+Od_40'
        b'OQ;upw1CW*&QPXp4p@qhwd%9H^#<^#Fq%){EUt~xgn^Wye%i`sz7l$}Q7yAM86Ga=L!x+~r@z`wh$HSlyJD24Gea=2qOe<(aUSkk'
        b'r$y%h_hh;mobBGp_v!xCbz|K;W^39UKVW!&8Y1cr&8cqoFqUciQ6;_@SF`jkfbTfMBCc!-W8)G6wRj$`V80I8D&^Q#5c%OS-F?5E'
        b'l8z4Ck~XyG$_q~aOs`hfHcOwx!W<RUG=W+6O@kg$IG4nhA+;3Js_--*W2n1F39QR^!5}5b#-~+MV=$On-pWeoz<)3_Y4BAY;T6B|'
        b'dUuEH0OYf=E`eku&bZ#*ZskOgpc_~~?Hf->%QgG;4mkqDY8(|QF6s5lac}A?I><~g`M_B6<{-~nG@EYiN%_Z!kjhe7+NyQEq{qW{'
        b'_4NJSS3l6SjNHdCCqJ#F&o(gReQ)XnE7q?pqv^?>o71uEgn5sPgRA8{ZzmUo(KrdPRP;;(2>I5cok<$jzDQGty8kC8Il<<*F+cmX'
        b'+p#bbdBjMXD_fm!(nkmP;53jU9;*_SWR&{^>&)#3!p0?KtxDr<Owl74P^*mxgmnGASf5thYH~{2nZh?Kezlo;TLn*;j%Pl8EY6;v'
        b'm_*w|i{V&c8Eq{Ytrgs0(Z-j4$l#aPfCEa?Zie?pZz>%eOPw_!donRTzP~6I6=Jx}6T)<hD&s%%rCe!YWwd+hO$xmhMjhe683gB;'
        b'%7q`fG7ez3x4kDuw<^m+QX$m)pOlgskh2h@X`SWN_km(X&%=MQ*dFi=mY2i06mr+%Ij*K3%z)L`hpl&u%aJ?vD@anD{~QxD`%1y8'
        b'J-w`qs#dBtdrtD7fyOIJAZZ3Wt&Eqz>)>Sqix_raH3C<P9+PjTY9ji2L<b-<qqRPUt?GjKrorkz6jY@_hVgmN-M<X;@N5pf#Wh)>'
        b'IRajN_uQ;_U><Co>UOH(dUO17IFS(wEjdWo8riN4BQrz3E{#maJf07+eNIOW;Z{!!Y#wy{5{voia^Eaa+Y1BCA&GLi!Km<Me{mw@'
        b'K^cwyOsz~$AI`bH)!(eZa~s#Pvz@9m$l*L}w*4wjZ|ko{Hn92sqDuVFP&^IfK`Ri6CBpHwIqr>@-V4WN8(H@c@sg4jO6nPTB8NOX'
        b'i3vl~d~!(5v(4O;?8`-e6tBW_L6Zb~XEy?8X1I-#TmE-02AyG5|5VNVq~|JPsZYTzWOqHBp*`fZ3*sed(}Vu18%zk3JaAJH>R}xu'
        b'2pF_06mK+db~9#lbovRwl%55sQpnioH4{M7oi?Nf9HXcsbnZ?83-_-&Djx_4<hcGhpDU%`Fl4olzGns;SqDqLce2`rqgQi;VTksN'
        b'932_}Iev^?!)12|Tui_EoJ@G)RB4US@A=A0%iR`E6m(3TpCTRqYV(-MnsQu&u%*7PoNj%0O-}KE-}I3UVdtmZmrPuFC)9=YrIOxZ'
        b'2}HZRaAZ`Xw%`hqFDDWK3qJj0;Sj=ze&oY2Nsa&BpmVY~rUVs655^HEZf7t^!|!Y=i3#o{@*>s$;^~9tQ@9J3wE2J}MWM_<Jh|_p'
        b's(700UFc8{uw^W0WDG`7!3yV^GpnsFI2+6SZJFf;g!Qo2fG)NGij;m|8NVNcQ@=63zc+0xi_Q$38t~+M@R-DhH@}4m*AYF17?p4W'
        b'GlD#<;QsrtW#>}+IzmY>3iG7@0TpKVW_1Lr7is_vArOIe%S^SC34$({Bv(^q#341e8%Q?YnTOG2%%j6a`15`JcjYtC803!T>)g^b'
        b'^cGK)-WP-3Oj;6rBCPKlczyE~T68qg;^rv&HMzdbEQz_jzD8!^Z|b}iXOr{*Hbw{5GEx6^^zi_pqo#y|l{>qFa?+`QlY5i;Apwho'
        b'N=fpzy@&8-FAB!2VFx=0($f9m%IfHI^8;-M$l4PEBsCNx(W%%v+54BdkNRdPJ|$KyR4T6xQ92G)sd3TE(niDEY2~92oJUihi4AS1'
        b'9#Z)~qW{;H9l<}sRt%*pyq(-4M|2%>yLeo=;b~`A-ErF??O8-uuKBGl?2a<upLst85+cq9+(b@{`0xsgZ*>U5;pNRZ#Wi!m?6e7g'
        b'!>id6x^AXPs@rj~lA3F-j?(y~`#*9yyISGbJ|YD7(}5x$q7I`9Ee#jhRUoFS0=6~?qTsy0da!~3uI7kQ)cx(`Tn>k;o~`D(H$7$q'
        b'`LX)=WGMa6p&Jj}fQlkpT}~u`R0)&)#T9LqiL*?nO_XHo9L_9FY6yf-&7Oh!jr3q{?9QW`CC-hRc~w|%iOvm7R8(LZM}#v@2S+iy'
        b'IyYHvTHy<$7|5bREo01It}728_l5J)kF9p=jFqHXavbUglAQ<K*BYBsQ!IA9_9M2ZR0!lamv8KvwJFn%nu^O<U01Lhzq8vad$<X5'
        b'0t&Zk44b;SPC~ox$9VTuAT+II_9&WKK~k%ac}B&}s@?bTW(X-b+3DBdO1{dHKG2r$0M_87LNyN7tW!7iBk8a{US}UHz?-gRHHo;5'
        b'JY+fY%Km^!reW}kxbx(^NflgOC!AO*NSFJI(*MqGC;Sz1xsP!d{ldl%KX=Aj@|it0&FZ+vQ2wO=G|vlWhPS{~+ASCl$M6KCxrCJ4'
        b'K@N|v!?69kIP9D6DO?+}#NDT6;qeUI6y{qD{b*N68m%7FBlv~-ABGu*Y`tW?>6w${=r`CcG@&N<UmxMZ3nPt(nAQ`3;&-;MuolcQ'
        b';rW?;P8TP?CE5rY2pDc*BXl5g7cD6M9vL4|9CO%ce6c?7{JBNC?KN$|7Ae%4O0(RJxdTdq(JLYsu+u!b^RT;RyG#p$gxADzG4rTD'
        b'{aETED~k1h<>5NpRXmf#A?HYdU}!45x<dvQ+Jil$%_h8f+121(1+}Z?kAQ{M6&o*e?s!_LPL41kP|fuuSOBIi6xL5;Z6!K~t=tG+'
        b'XGOu($#y4)q23vCb7X(0;^!CY?$r*)!mTW{_7RF2*2ULr!(h6t#q1_9pZBoe0CaXpo&=i@N?v$BzT8wsC~6lfZu3A>k(%`CMWXNB'
        b'UHrAz5+6UlM@7E8Va5pp&2#h*sM#+W<bW0%cTqF}cf)cWOPCR#+)IKXb`!Ay247^#<1wNAi<9s_dZ&Dhtqm#L>lSR0c(ykOi{3<h'
        b'&6`h&F6kXWN%Y}4-N!ebA<XFH4zMI}>z!0#Fg9yCFu%2)r67<a4i=tXtq9v|<7WEgi*RgE=t9%BeP19N+7Q}45MGUoQiH2i_6h6e'
        b'KMthc$cW0z2EcUEKzO}8c)JEWWoT;Y7|hT}*XWVNBArK$(l*`tsUK};NF7O^!7{sq_m1I+e!IkS`G8QVdy{1VrGE`oB7$cgv6gQ^'
        b'q()oshxn&4_4u)Nu3c@D+o#DoyR$}cthWaF>ex1M?^1y8(_%5y>2i9GQwrz`3Uez;%Z`e;$OY$1Md>Q1GZ@?xy_00xT!K}qOi)@>'
        b'd9&sJZyL(BEB-sw`rT=aKX4^GqFyY^ZWCw;vZkp{g+kt)n&m~#M97Q{piztle|C&__+EWGdC&PHT(#UyghVkA{_kW<kc2-O^wagM'
        b'R4WGOO!3eu{V}^pC10p953`jIPOeJLKd3%tI(*xjeKvQ`uP}t=Qt~}s_y&CFTcsr~N3lqHib^{MKJXW+{TNYt7-3mEy9bC=r?-{O'
        b'Pm~9OBX5JPHnfvg)dT3JEPaXvKcTDjC*kV}n%bFbE(IkrMA6gKn9s!n((n+vVbN<{6_foA3aL~OWbet)FuiJIrN}L{k*sIlL|{GA'
        b'2-{^$b;Dq4tYuE^xGVu$VxlSa29hVEQI@KG#v4VX*vz^u^1>-m;gb%h{dI4(Z7n~Ck88NVQyoDT@C3ycrkbp&03v^wI=>v?PYvXO'
        b'K*jFfU1(M;Rn7ial}%7!dR}i#v{1;{bkNMj*c4=e7EPD}=O)r!c^92)9-}1lE4dI%&1}Z!bS1YoGmo8s<e3UL1M-*PNIrN9(mwRU'
        b'{TcBj*lsb2BgMw5pR9=C++Fv*(dpVg>MPVaUpvo~A00aj$}uJWb)IZBs#(A=+{%EX*bIK&xS0eXg9tS6t99dkmbBc8$+O9l_)I8('
        b'4bRE(5c5HmB)AdSygmYrvIjc4g19xs(J1(*jWlwyIF}be=~Gy$B49$+v4|%1J;qbNhOnHOv%vNwPYWKO=S2IZLD+5=Un=k7xoN6P'
        b'MX3gM_$R@}D5hc&um0AKq8Ahl6CfJkfM=`*P0>sCtADg}JWC#Y->9#Qg^_MSI17OK)kK51^m+@BZf|hOx5hO~iW`g}%5P`{-aNOp'
        b'_kt{n+nxcapveq>Rh}lBqX(0D2#h)5Xu6fgtH0$=`9hut@6b}GKuCej1$52k#q-dNmE`F{0&=QT#a?xjlI|y~3QP8oyxS-41(9U7'
        b'3#)?+{$g7dtNS8lj6^u%#qE0xxM_CDW&yB~IYC5QEi3No7!G5*m?!^7D@Q02_j@0g7+>ouGq7~Os3wk&?l=|_XQ!vp8{)furB1R_'
        b'+LNCs;I!vF*V8;4iQVZeQ&x#R&H|iF!io2g^zd4$TN8ClqFCdtIPwR`3p%F1#)$lwIkVr#jpRPC5t*sRf@sl02X9{AdTQjZ<p_N^'
        b'-U;{(MT5_lUB(*6r8!u^1Qcn_wR4ywuygn9QEG=Lg8^<wLE9;+J0LS1i+>kd`Q+>oX)T{OV+ihe-b-?`+kq52_@3S?9Wgx^pO}mx'
        b'6}y--)O-1IEf{!HOU;J=_$|11tzoOKvZ&w7;bLBl#(@;HxGu|7_PnpqMizDf4gleW%6vnS+@H9jK56IIIhMis@1;`FV17l`UX|op'
        b'M^PJa2WhVqJc!Y%0KbsRu)qW{3Z|VFc&32<U8o2<_T$hHcMbW9sy`u-e{v+oRo)8#^TkGHnl9UlZ~{0byQ;^VHiC$W_&ws4E%g2Y'
        b'^cF488Ox1hfHnH4ml4M-OE)R^YTL8BQOEKs9;Ag_)H$Cc^rti!T+!o#w|pGUai|C^Gmw-*27urP-(wNk;{s{IwSo&rEK923h><{;'
        b'yw6kl=SkRIokXUw-H&d2nQIk&@$nt5CG<WV>2K{>zJK|xKl-V4VMMbY=nbgtnJ&CjQ#+P7C+r>HMIF2tur`}DQ=t#BK}b34tdFC1'
        b'^SKj&0{*7msrxCqG*+vtr@R>r8-bC)#-ZRTKQ*>b{H7(<hRd;v`nI|j3{Agc?5+@Mk?X{i%&fuqzNKLEI_Vie;s`V};rmRr<X;P%'
        b'{RdzM4G2S0#fK4u{tNU7s?*sV3;`d1)Ns$>^@a-W<yZR>1nQ`})GTp5UR>BQNhwA1tm0r-U!$qEo7$D_usIr8(-_qR#cvKrG(7Vm'
        b'ha;pwy$|$?T?XsUv^1F+k(T83I_xTE^l$!FqL0s5ToPebC10Z{=jBk&t&eM7g^GDE(2Xsqnc(V@ceD*|!Gs7&F%7DAR^^XXx7F7k'
        b'_@$Pv^w^P9%gjX+`ldl8bkZesN4xtiUzad1Ku&gZl+Mv|u<<4(1kgCqi$jN+r%gDRl8hL6BYM`_VB}8K9!L!rOF*U1G#A=n;U09j'
        b'$<avdF`0E#OcE}m>1>cng3{*hLht-O2=`B^Rf~49MgZu%-u6WC85X}wE6o)-S8jl6{Xz`X)IKqmd{d3B;ue~Lsl4L@{|Cv7K%4t+'
        b'k<;%nUKBeOK6U1n+JR11ju%3s+eyL}c@gDG4WGy>XMRnZyJ}*kGhIUd*l;}4;>Dmh*7po=+)PB4GDz%bI*8vkRc?F?RSDSXm|gY3'
        b'lU|w0t{1uA)8b5j2OxT_qkT4c?rAR*rYfbsr3366_!2@(2?Q92Eg|%#FV!WD0P9^J>%k)-yY3;qkShbjS&~LhE>d{Q%K6+AS&PnH'
        b'g0j>o`zE=P57L8ocp>NLY$_DTBFFCQm_4C~q#+r+*l{e}2PH0LqYrHuz`C83$WX;}?J@xU>s9~H{C(;r-fkv@;%B0)yl6n!-G~#1'
        b'A0h4?NOn|Bn{^wEN$=hL<v+V+m)?M6axI|KUoSAb^1XMVAM!A~&#t~fQ)OiyLy~G@&^y@?{nutlI(}lH4DvE%T4%WMWHg(Wu7{;7'
        b'8(uQA-;1IFpGe=^YV$7dBcd`EpvN!%<uQDm_A;OhZ+Cq3h_2WG2B%MGh0Drn00GoRh{gc`x)?ovvBYQl0ssI200dcD'
    ),
    48: (
        b'{Wp48S^xk9=GL@E0stWa761SMbT8$j;9~|NI9&jGKYSoc`U1P1$C4l?gfoURIp_X!Xq|${fM~kersD3mb?qcQ?w1;*CTLZDBk{Tj'
        b'g;VMjYo~`l8p~dK%v5hsWj>Gwn>t!t-F+pWz?mO4`-`ymn3jos>8n8JI)$!h!e>CUAb?KFx5X|<8l!V1nVA>h<zDeEJvH3XCCpJV'
        b'-7Lb?0PlnbCnU^eytKmKQXV{vxhr+#tokB&S`<U2hznB|$yr>gLJoMUX>skddKX!%oCs~C=%;*#P4uj<n9+KE{3<(d(hD$Ri{aFx'
        b'hs!W*iEkXGwopNsS@wV|?1=Md-xl5NS!<{*|4+aa?itasPB?dz842?1E>KA*zv#@bJk7YJPeHL<v>7dYGS|aTv4l&~pRu-=1rk(H'
        b'xiM&&JSj3bK5-VES{qhbLoq|M!++^jkbk7d_^R)7V3jedgAqa6VsND{fZK6}D%;1pHpJI)vI#0)EJYI5kxgho#<SU@e1bF6QA2gF'
        b'fOQ3AjbhS)hR+oJ+jL^h1@42yo@ZPD<nWu=FwP`yR|s-yEAQmIb!$QM$`FCg5QezV*9b?a{$@FUah&f{p6b_qRb*j3^x^g)2K`5E'
        b'v6MbCFN2Cl@aA9!)C`KW7Q{rOxLW-ErH-Mv9)&>l*lrW;#-|LT9OZS~c%Tq*|D?0+5`(k7Mh+Vy$rHxs$%N?1m|$KdOnG1cOl!|F'
        b'6}HrSna@@~XDMS7Y0-_8|D7v=d|%?8$@@pBi-IcUXLM=fYOfA2XcvUEsGuZHUK;%wm0St(4?0^<8+z)?p!TNXj0Cy&Ysn>_V{AHE'
        b'!Dpb&bSxz~pXZa?4ug-mk1h<4bc4u6?PKb$+S6+8Ey<F^znDn~QGSY*uuYy1I=AUpsiZkfHS^u5Tcgd+epbg)@5+`k41c1uIo;4L'
        b'fNMuM*V=^HBMM!$uu$G|&2ZJl5Q9=1C<5ffV60&vV&A@!N>LX!x<plM667ii6YFPt0hHrs>rfFr9MSol5N!c6Nhum(E-IiFU4>c*'
        b'S;3^Ui1y71^<m(=pI)OJE!*~b<u)JDbZyji>}+eRc{}Ujy@P?+zc7?|^B<VEf@GXaGu5QoMxKK|=(vPi&_3}^ZQ^IbQQ2_^Y-|<>'
        b'4a?ZO(<)t>LdB>6ux62YjKnsdUbdO2Q<MBnA0Kq_m}}h_eoaGVqdrKGjr!wfT<@XLCDtkfJ^g5chnvFOTQaHl10kdF&EK{$%Db_|'
        b'N2cCTeMupI-|x3xjVDrYs<i&Cl57CNc})B`_Mph8h;Ex7*TOL&#}L!G2kR*(o`EMQyml<tB-hQ;2t@U##mmNdrtH^5P8Gz)EFrib'
        b'nu|eKgOP|rrEL32qogam<`DVZ-qVR8{JtM{GEO;Msv@{vmdogbpK6<_Q9R3mtS6uwQfAs+>;jcS#w+$^wM@F(LF6j#+>LwHmJ*`4'
        b';UfU6#m{!iUJjmjGv)Njfx8Y{?t{Pr`jYFn&y&-*AYRBocv#_EE8a<#`XMXpNGPziYa&4wZyK0LzESKXh*xoSx~^5o#x;#Sy>>Yd'
        b'8$eA0E14@#2ABzBH6=xfP_W0G6HNByzpt;a(69^CYP!db%c5p2J}p_(mw*zZGgaO^Y+Vrvl(-lHp>i%cr)wOnO9le-Y}xWs?b6!e'
        b'8f<Nnk0X}?j90toweG;+PMa`oo-~Z(pT-x<vte%oSTYreliH9^mJ&GC!}IOI61zQKY@W^x`DoreGn{|7J}k;beOZ<PGk#5}={FR0'
        b'39dCc_4{l`G;*B$*L7O^u1`7zOH<8ao=CfaR6h4uGZr54Bc}L|hA5WCTp{6~Qr2%CczI=7>O=oSZZe<}n<#h^K&SC@U-5Ox>2Q{%'
        b'-5;Ve01t*h3dX%my3Lish$>x-(3k^o7><8QH)LQ8NZzytX!xEpY4aqpy}F$!M`b#nVeA=NA|{eOHn6eOW-54aJPFv8KhL;)6+)=F'
        b'xA&zq-+K=(59OEDWozUKEzZS?*c8GKulzEtKhv#g!%9j8$(qc~(3$`FnvJ!?JHWIS1T)Kx#qY^PYJz2i2%K}@^4j`C_>PZdG9Z5Q'
        b'X-WW@SZTclxUl&R>tvF({Vci?E!^|EI)EI<Xxzm>8!vPcQNF(U;F_^3c=IEQ+4hAPj#>KwU`e>iWC3K^YalMn0j@7MobY85i1Xnt'
        b'sMlw;CyLneiOwn{S`7Em5>^q;B49n~-SLd68%X^qne5o1dvKS$-gx+Mndxf;MlIJ$1|bL0pdP`u@m^ZNrQO5v)J7ZGbWF#JYJS4N'
        b'NY%$MnH5e8iJ>~`wvyv{bYw7s?%XI**wqvtNpSMV`Yx0AZs~-yK^izXsysDCw`tDwPonmFtmc?}ROwkB9efmKK2cHw2wVTu`CNFR'
        b'1*&^@Fwrp}+4e`{@v8TdLPu1%(Ti}ZtS#=AP=6MwY|Auh0cOCHK}{ZJ85>GFG`7F&bbz}?QIS_6X+Yj}6lme{!eX-!UT;dLBbYl*'
        b'GiI+*kfdD1MThOsF`MDMyI`Os;=CIEEKyetA8M&@-{V{?4JQaw_p1f{jjNdYYanexma;sD{j!EltHg(CoYMoyx)w;cI4Jo{litOd'
        b'{eZHuwRgd&eit?^r(@8FS&bgYrmIpi3X{eSqM}-p7pQnysWHeXW6;3+tkX%lm#KoPtm9Z70aR3=))0`kLm#3(u-E{xu&-4~gA>Ed'
        b'?ohK#xDm7Z%cIm;e<14rsYZUpE>(RFK7g6Kvy|W>g&>vLBq#7Kmo_>rUpc^FGoT6_yTDqug%6)~E=q>~Pu`;T9LJ3@S*@Rml_lQ*'
        b'5x%U`(jRI~?r5zj;lK~!P|<hMY&F+iBj6LH7m7V5O%6L?>5p<b&5A+kPo|)x-b$eUvyLnKDl_xoeWBYZ629V#shm`PIwSovTV8ml'
        b'qF6vkWpZet%A$Iu=w561bW@S4qgG?ezyR*yVkam{*+VlQ!xc^i_u9`)y5R_=GsIDgL1j|WZ$%RC%KQIxELN0Abp2-B$ISH%>aC#G'
        b'+*<X){fBnAZO~LES=-1BvM(aTaZ0tRso39rC*LVY0SgW2IwN9<oHg7!xCE<iyTkEZW3ukXQVFMgP)wX%<ZkxVlUGx$+Sexe_4?;r'
        b'-qa>S6*LUeuUxb;@wXr*T}WVCZZp1G3e%F%O@*WbBbQZSa5j*^wY(m{48i)_a&;^IM-j;%Cis$l(NUNoxGB^0_3l%H@L+mB1p!$r'
        b'Knp>{QxVztfCJ|QmX~A1RN5%4>%bdZlF(uvwuMuBKgd)4FE}60-E0qn!E-)6^aZ`zIAyI0#Ti{(x^zWRO`!1G+$k>XE&lw-O_g=e'
        b'Lb8zos5j=<^-y-tnXh({g3}lyq%yo;z*8NCXa`Fw%3#d&5L3sNT;-y<M7u^P%1@W8v2va_f`YO%x~+b}a~Jw*oC{QwahN%!>aQL_'
        b')u&Sb#~E;>ZM^=B5PO#dBoi}3D^v=7!$Q(n=8R^3i4aN_Zf;L+`0_}?z=JHIp5SrcPaNk^i=hRhA>}(3=J2)tZ)7UGj!5fOiEHcQ'
        b'Yud-8)z0Y04b|(<a{}~Y^;KPGWOAZ29R7;iDtmK2j>;YJM^RoszZW@~@g2+BakXZr_80TS8b2sa+Z9v9^sh&$&)mX<^+h_7eSQFX'
        b'Vfh*q<lJlS%5nrtjg;Z8?f>J!`Z9_rSm$3~wiSqAA!)voh);SkeHnWWYs>`ihX>pnbFa8cabC+0b{48CvT*EyF=6orm3ID$`ee26'
        b'*LJ|(Fm|x}1?g)RXtpq^>i%8ZrM()?88O_LEJRa$i5xI0i|Rd~R>1&>0GW1p{WAuLcbFJzVf9bGx>9Foe(_UTlFHZ6s4$-j`GYDh'
        b'29$3oCO~}jgBJx~5NItaO_X*<Qqe!^l?t4%k*@?SuN8ZrE(VcUUDz{nFN~1|Il1htW^A;{YK%1yRsrD^L0MwLnZS^M?5@U%Jg`1%'
        b'whV$%GmmkPVmzaVGh-ERGW~?x4Nj=<%cHG^JjxiJ%Q?Xl#L9|af-+4zb{bm8vay<#J|g0m2SolU5~Ra|iB;uQf`x(|7<Mwzoh&yq'
        b'PPCa;t|bB|z0HfRItbQO@(nR$bPDJH6vAWT!F&PrZczco;<y!7120fzZtj{A2n#94)%y)^*U}h-wpEULVa(O<%FBAY@z5Ona$X=h'
        b'q5e%)fv87k44a36JM4-HRAe=q>@7EJ3n4IPc|l_8!Ev!BfGl+Elz=1|!SWzK=~E4>cR}lkcu)TdEaAHgQi!eMQI_e+6glBYaIGek'
        b'Oo12iA3GV$Vj`*0jRrtgz<lwt2SEDtzzhwParlqmewWQSs26{k#)E=b)gGs7+wc|w6J8#eIeV2^Jzn%Z@)13VlL|nvG))bJd(MU('
        b'h7P11uJGl!%|$Cwx3DS>haKly0JzhOr1S9&Wa~K4^E7=cxTZrM8$s4+O%<Dsr;4;ICYT|}IZ;wTR6j<&*w-<Soy%sPc$p=ve<piB'
        b'g5iDH$Ax@tmC5Bxj_50RWS5B5&fo)_Fw=f;#k7N|83Ha-F#N`xu>)z58#Ch-$i-FM7O$*%*%@PBE2iCak+-F;%w1UNsGYcZ@ch(G'
        b'y0Daeb3W=o*6kfcuHNJk(i(-_qx=$QL%UId7I$qv4YsVBEm{K16epQxo-7T`s39k#l$MrO5bh!Hm`uF;Qs6p8oHdxUl3yV9guysI'
        b'1yLaO$yNs&)U@4}d_1k{S#&$n{tND;U2CGypfX5oYo7FkpKHpiUpn$9%_<V`BSMPpS$6EJ`+(MuWo?L(c_7Kp?ftnV$&9e`DF^)U'
        b'-*(FHZ_z9<Lkh-~zw32euUnP#A-p-WT6pzQPk(~2WHt42oj2TC+EBt_oEc=_K)bp3oQYWV#0!or+goiT!x^YWn)3%01p8DS;v_>b'
        b'e#OtmyL+X~^DK;}$kr=<m|tPP5WVk-|B(?Qhm6HpzKn+}(}9&NnZyi3Qz{O@IHkMQdndTt*VFn%W9wxi)>1(Q08nWmT^3Z)eFDPJ'
        b'VF#fY(7zsr=qXg0a-MbgCNru4-~B`B0S=RkMo%mE7XN&@w_%_W892BTma05PntAXds6LxC0Z~|I(rF{(3s2?!F*vT4#RVaQw6Pme'
        b'QycQLBg*M0y%#a>DcSWv%;cnY84;408MzUDclN$_v#eOY7oEI*Z^2zNySgl|ZLw`}&%%u~(3CV`dNm<%<gub$_Ofefovpp5te?Gh'
        b'Z#&%5{3hSM636uC0X<t$8nrB*Y{ZcF_*Rcsm)NQFoQua99q>&9s^+tsWI>KjRe0J1RGI6vI%R8Lyi8`lml7BT_Q391Py5jXPn*`M'
        b'sIh)6n1Qxey(|UWNo|=#pTmA<r|25m)d5i=ymuuunN^YBQn$Cb+JnjNED9}2@?C2)Ezf=3`4;iSf#~7xqP3ri$)7dJpXexp+MM8P'
        b'@ov3t5r?y2Y1`BJ@g3-O+H+dzS6Vd6XM)q&X1e8TICjvm51QZ3t_Cs>I31Oe;(SIXk?8Uxy}9WLVKWI^Q_-Epxaa%mRl~Y^{fpnd'
        b'c?F)fnMtV95e(MvY7a#Fq6re)i>(P@tm7Te@S#?WhHq3Xmk@<90i5gR(PD}R^p|Hp5<W6J$NSh<9g3{LYs~g%Y`xbD-|O(E!P;Yi'
        b'G~`km21$Q5E>}OW7MQ9v`MAAtrbm$CUtj{l!e3k9>G!}@9fw!EL}&mURsD?f?iYYUh_ydT((<>#Y&wL6KRqICT<SG83T1UKQc>mz'
        b'B<Tlfa+)3;C<A#!s6}4O)o7)1`Qm==hA#K~u?CC59|UO_<&lG)C$9P}^74{qAd?!Ee;2xjC#)dldS>_*At3LRpZ5IWx<#DpWcwvp'
        b'!7C#N8dcP>*72}pbRgLq_MN7aRSdf#W#SeItdt_|Z({>&thXsbm9;s7xPVLlq6q;Lw;pDjcKX#bSXiTTiq1Kx3_&=>XVwK1nm^X*'
        b'7vUR*<>=5guLA_|pC{6?%5MBzEwLx)kmjW*J^Ik#_2fPU#DT#v9Yf-4W028?BvRMT_qxFx?=Ku)BB0o6wqJ9y_E&$oiktJ!Rq`#B'
        b'=od@oG@1iPkq^?$m+#gJ=zDCpJp!|yRJSHw^W#Q}sU-IHs7Az`6zSX2Xj>$k`}5{hMdW$`;Jw?@dylmaLZ&gTB*2*EX~&X51z&@o'
        b'ZO;;{u~WM$iuZ=braxc=YUWkj`<`d>achGQGEC>;5I)WCdHgAfQ>aLelF6DBUXE?wT9sMKDda}X#6!Zy9kiH?;IJvpsJ0<-r+A0|'
        b'M~lRW65_u4^W8p2nGLiof!wi5&`P;(a(saNzIckU7)?@NLR@>J=FZxFpi;aIP1o|RKZ?3#$~sr4X&OY3uZcpvX-jv7QiR>lyV8@N'
        b'ZDqr4EQ>C(JK1PgqCq@kddbrmtn>CnYxF+|^)R#zc_kw(Wj#{p41g@A06$ZLlbZ62L;12DkuK6df?Bz#m)aK(Kp5~UM^=z3(*ryl'
        b'RE2_dl3ZRlAYf+jZCa*NtoxC)1_;05c{5SP+TczU#>8+oJif7%RvzM>M9r!$nBb75kq%pj4_1u3tG1QVs$Wl0`iC%DLwij@D>c9X'
        b'*FQW(GnCq!ZfNf>{&P*km`!N$f3p2%f8`#e@}C@PaISrmZWQ#URjK2!%0L}8yis@+dJ8b<U=^rj8Tb5NPOzD+Cx_1ReilQshP<$|'
        b's*Z(If<BKCWi1e;a8>_ny!}iHSMw$ccuEH#88$P@?Mo?OrLhjfzWF4VcR~z~>Rbj$B}}7yytZ{!Xkt^)H0ljQY)laGzx?rUbg(L0'
        b'SgPHT>%~5q9&f6_xKs}CbFY(bxyNR4riJF?3;QB_wf#GoQQ;aYhggcbl?0s)VDSr#qYm>_G}`F=6EvfYH~Y?X?1eijSGb>iyp+y8'
        b'Qf34i>HXnt9l)1Se<*}I!IJ2G8^4o~qTk`KYKbh2WlE>j-u0CT*KH;ih2kl5@&>vd!4`qr;n41jU_+1o9472GtsI`4a$wDQXDtn='
        b'&|NCF1!p%w`&aNHW?DE|VBGW7{7Hu@sG0{jllZ`2kxr5c^0SNos=1JLC3)MjoHW4Prw{E%7#o7%20xx?y=YnNREFim$+^RH7C{?K'
        b's}`6u{Wn5hp~>P)H^j^5<=~1KHnu`H3s`Rw`R3)HgClX2(J^_Qh{Y9~)(NdJEu^Sb@N+ho2_>r`iarnSmS2<RygP3$r|?1EwaR=U'
        b'@piTxQ%tUqXhDD?Em%lMoarG1mzUZ%N5}wBQvLuY36yUuaJwpfVb+#BL$Q9Tm{8|!LUQA>nu=_CHl-!{ysIw+GRM>s;S;vK9C(ao'
        b'&8klhXXHoYR$HSo#g}e8M=JsuTJYhB`DtXclpIXsE%HC=5$S?prD_1|77*D<NT2n(s{VmlzV2D7!%a;vx&xc{H#z?#tP(H0u0Y61'
        b'S{As4Qtxr<2vZ0M#}6qbD#O{@rt4X{$Y5YkBw5xNg>Cuut5Cdz0o7bAO7?vt5i%*l_;G!E2j!14WJ-I+JHG~3XeLK2iOPivUl^*n'
        b'F2cM>J$1E!hDzZmI4L}t3}=LzcMX>a!07RcSeVzt>wap=bUv4r1LI{mK{f7f^{rkwGI86dKCZuHRn!1%XPE-!(BnMcgJuBa`ApUd'
        b'<Pm$psGo8W-(Bb-n~%jZQ&kEFQt=MMj6cw*OrR1qA^aJ&+nXwQNj9J{E^c^n`ZSJ0TlMVKUCtHhl5)uq8!xYRO<_nC4qTZKM8RPZ'
        b'j(xmUpw6KO<qwqVOsow3W#ahRGvWn-juiK--xHUWC!HK`!idK#?sCTubMS|-LT?9DGzYTXy%Q9^i~w#Lrd9hFQTckvlkonous>v6'
        b'3OiR{sODiocIasZXK)Hr{G+aF-YecL1v96U)PvaNg$;t?PEb<<1S@lp@#LF+3m=V@%yG`ipeA4bx1c{~rc@a`9N6OZ8C}5e&0ZLk'
        b'wH=e9&WOPmRPM)^+AtaHo>0ufS8(hnU57kBr|DPklSqZQnfWuis5L}wuF{vV9;s{trs0e4TMKWiFVdZzAxve;Jg_Wx@%d=b10dB('
        b'#GU1|-p}Xk0DJe;_{ySEjrL=i<bQxbqdHW5uWb6F0^<l>jDKL*)K#plkh$ZT<obRg4L{2U8Sk)Z#e?Kjw1B0K+|6d){S#K5ODXv5'
        b'HGxK{&F+0ND*+fMs_ChEJ|pEoD`y8aZpVXjno;Ue`VQ%jEbt1+qw?!obh8$s0Na-V(E%B0*H4tJG(F~RJ7I2EdH`HSg-o-8LK>gj'
        b'b|oVf2$IZ1QQa9<8C+5dtJBzSGnWYcoHbF5X*Mqrswc-#*@5YyS`d1>zSHu2!%+Z;_d!_qOa0>JBG{uQ=CyG!`})FCrooGG+-A#u'
        b'#%zCmtBJ3=Vn%<k!dm{Z(5`>BS{?|>RkO>vu85$Kf`~dpxnq46vykdPtoO4`?|G^P53({l*qrf2zFiv?G6iM^Pvyke*T!)Q<1$dr'
        b'jp81}I9%!Owi`{i$Te|l>Wfd-oYs#OT$)w=k_Gx1QGdY8nOxCvoNyP>+x=(`xi1{B0nO-O!NRm+Otsex3}lzPNcKq=`2v!~_#2xg'
        b'D*=wBd;Wz+@1?i>U|Z<T0zu!=m(m>gb_?88MnwG#8J%-U9E?0Qb2M-nb#3K-soZ?4I(8-3K@#5GtA)HRhcjK>8^^Lq5tG|#!;64Y'
        b'Xm?Tcs;x)rfhd2)7Rn5|XxBhZM3v70^(x?JhnB(4`OzAE9w(v}&!P^5h%LcYig51<CO?lTHQO*!)RmUgH-Cw!v4G*AKxbm)dipMk'
        b'a!Q-`>@#+muAhbLm<yfi(uk{x%j8vx{<W|76{k-^r^R~h6g9vNpG`VBF!$(jY-&{7R(aj9gC}<56<alc*k;9UN-@3-ZJKqeeD-OO'
        b'B)x!$;_1hvBFAY--sBgj!OAJ;s7gN;@$ejd@(9Z0wg~>}k_T^;#e$Usm+L{AH!Z|+tE|hJQtUgff3uKpnYE;z<^CPjzfKO=(Yuk#'
        b'IWRJH{>nBKRZjgB7g1qqYC7I`99R}ES|uG{T60I`EBp4^0Mi|wTZ-b?O=wr$FAS0VpLN-ddnD1gUT-wOvIX`zZi+?d3heS>oNqDA'
        b'd9d+H>n^jmT(4C0lMo^YHqaU_sC{V`R_=-IU{iu5f!nbfy1G|nlks*EYAQ^*7Y3Ea-&dgEXtI=02-a$!EB)!GF?~Su9a=cxwYF)S'
        b'HBsnYy{$dkcZcoHFXGcSyH5JvcOMRJIrFfC&a`=hkw;cqX$sX=69>2!^@v~Y?qhA!cn%Hb{Ib_Q7dTN8FjKAx{mV}v?IOYe7<8w8'
        b'?v7da%q*lyR}syKBxWYQ(+^`k+%&aBurTz}z^2i+cCsbi1^MfM!B8pjMpWJb5Fglhwvu?ovXyd(KF6dI4#Kmuu{2KVMuAzMEIfgJ'
        b'rM2Z0)I5vV^gO4kcEsPWx2)q#;!-WS-LRg&<8T%$pKk2Z_m5zr^>Pvp7uT4l?$~-@{xnVaL8*<#<=xvY6eGYDy=}WG6Q^2ori`>4'
        b'%_vw1P*QftW)2Df%AmfFtdO*Ug`a1Y0)bLpWG*D#r-SS7aUQ!(P@YrthIv??Cr2)zHWJ(i?uIizVp@}QLzZB%d8XhqSOtkV8b<5?'
        b'$Y$i`%E^h!$31P<AHt#|#Pjo9^OgmDew4nqCHq+(G&ehMdJQG?fycLl_7cT#1&_sEbsHWW0B1St!w;b@cb5iV;iu@T8_mGBNrgrB'
        b'pRg}?daJ8@E4c`UoGmko|M2lEWw%YXk<}T!CT(fJ+P6$wVqA2~kzd<!WGvK6XHMCXsF?gL?Wl8oFDN?JXr448{60DpZ(fa3^Qa|8'
        b'ai2b%k5-^U-N?*A76;P9jS5(0H%@++cF!;O|68q7X8T5w1g<;2|1Z_J<I)JF&>*jrS3O&ygt9N<f8$a>BuIHb&1JiSpiyrV7OQz_'
        b'ph6QwFtr<Y8$WyJ?;7^A01Vj4zd$W)LkcU?ttHx93F8&W`$baDQs`fKDsBZkM?f}S=1MUSjP%ZPpxbb7s5@#RX}r1Vt8Rj)z4e;a'
        b'dHSG@Th3V$<tX}%YpH%5yR3332fQX-``TQ<5UcQ?DW%MZGmD}&9{N&D#eyuKt2Yx(j@@$k#c0r<0sC=XF3}JjAG!i(mE6e(O(DUY'
        b'wi3i}-pQ9YN5$bP&Ad<u)CgMNm5G`8KF^RMZVWEw-As!iS7arVnN3w#n^7aiQb&-2Keu-K)#)8-WI>-484&`Vm+=I*&N<L#HW}|H'
        b'9&yE608MT74i0Co%&({>-luX<$JYT0HSjj4vFiBDFNyVFM)kXvTJV3Wceq=Q^*(m2e-3^R2cuBPI@5^nl!0^?`|O-m8FtRhpqOHq'
        b'Hh3#uv;da3`m++y2-kv62{e>hi$wl&O*5prem#Ek_sDk0x$R!tdbQ>w@^mCF3;;eNmii=!>L)WE)*Ft@*e2)wsRX@#NPd!LvXs79'
        b'9Elid#p0HJ8fS3%X9dH8$<H>8aRC0RgdvNPM#H$4hZUT36aR;zMROO*8v@t2$+uwXzwliJ`OPVGmCKjoojDV?+0}>O=pcn<;Q`x6'
        b'O_gd!3_We3XQX~f^jf&i_OP{MF3qR$6^Ubcn^H%qniB{mze%aEviI$+$O|{Dtxb3PV<{GokY*U;harNtSrB_pQBurh0S=u+*HFB)'
        b'*$P}u?BNLOl$**+<1X+<csq&hy3V4R_UQ%??8%#Q!2>yUH$}L0aq$byr)@HeLl3j?ri5Me_5(rvi2&TJ(ihlcSp<+>SHWxc^FW+r'
        b'PS~FGkHpWG_m}`nDR#02{4~J8oub!mw32ATb4J7Yq*3YLg$m9fjgXC-wAs7E#{o;!ZU0VXv#sT|fC%N~_~46HQ2fQhq9RSwWr@WI'
        b'k6JWo>#z}<gJp$MTCGlJe^W<X{9ek+Cyg$1+dxXvF)s^FFA=Rlzi4>UX3&XWhfR7Stt&&D2&7I!%CKPI)ce~Kc`D@M<dC13194B-'
        b'<;aDpg`rQ%w_2{x_~FJ`5Gzw!!abfVWQJ7Sv4d1j5Gcz=`X?P);;wf+HO836{)l9|hj;K(8G_Mmf7yPJJy^V@fSci|5}L9|+4o!4'
        b'F^Aw`^iYwnmnm)<F)r0xRkJL2Ku$M>ye|oB2FIP{N_uc)JLsU~SXOFIzHyZ{f`lj+{z#@^yW>he@MDi4s|7pG5u-heSzpCzs<w=U'
        b'OW9w+_~bza`9JNzMLA1nU6d4YPyE(ZvbjR2Aukij^^2U|rhb!K1^BH>g`AZ@R+gSt5x;BWVm{8%%1L;NMLVJGbry+pt!_QDkwes3'
        b'pGC-!Wn&d0g+!#6ZLYScE8#UIKoaZ@mevA;v3`Sja`zz201E6`?u-p%>WFcyd8sb!@vCgjtiRe<-F}4m+c=(Ox+vv8W7&uV8K1Cw'
        b'G3w@&p*64#-n!I5?Z=oaUrJ=X1-dCf#ukc;GA~&FmV^1b`qJF|Gy;(TQtp$MfT+<c(;_4`T6pA0Q?C7xIGwPYB6R=&SZ^mHmH3k#'
        b'00GoQhsFT_Tp>OivBYQl0ssI200dcD'
    ),
    75: (
        b'{Wp48S^xk9=GL@E0stWa761SMbT8$j;A00O2webrKYSoc`U1P1$C4l?gfoUQ*J$^wV>}Yzh1L<j2Eb6bFXo+D$RmSKR@a}oRS;c!'
        b'lCeGOm6Yd>np0ep+)s^u9TO=Li~C$xFS57$5ts@#hz6KhiRG)KF>~6x8k$Q!tZ)dKbc^Q>Z2{oRh!9e=Ojb40ckI-%?Q9Uw6Q*3$'
        b'2m{C>duf<nHM(kupxs#Uwx`M1h$qew0h3RtJ}&z`83&k8I`2Wa@qCYL&U&EzC`o?25QjU#9uHSK&?Um(C>CWNwee8v0SIxAeL%K<'
        b'o5VLDWI_jnSi?RIz;6?G?gQ^}YU>hQM7-f$sS_MRh063uiw`L!#Ynpak>62_u;)w+rf7rkA}!-VwJY8_FZ{JpYaGb+L4N{$J-%Wv'
        b'LTj%YK9&W)12q-jI}po@1iINB*%xB;D{WV1ZtPe2c!{RQOEe;IC-~)YV`enQ3R3M^rUwRtLgB1e5<-m4_AD;I{vXU};_VF{F#`zs'
        b'K6z>{W!sp+Uzt#;=s@+xAc5xT)I^D!4h4o>`IF~g3k-Rs5ygOsh?l!%-~|0S^vEliwn5wk0xYwSo)d5H$V-=h!IH}-mOg%O&f;rb'
        b'j6oLF6lX<jq|;LpwgDPVpsy`Daj=k98QByw2U37*WCRh{nYZcuWRWR7nbFwh16XUGvZytIqK_vKF}sqIbwYkHa;}!ryH(3iRy(j!'
        b'82SLC(36HK`D~tF0di8;DiVI<3x3UQQEj{s!;`FV(o&O0XdO!^htQjE;6}VL_<pl@rcm~JK!p!61b^0D1h48*+6bGdNTM#Etm>{*'
        b'7Ys2M?yOy+WN=Kr#Uaq5dD3u89B!09`4o6C{en{j!PKfxR&n&xcp{~-kU>qRdos~InoH0R)50CBdXP?D|6bH3Y6o%G)6s67p|_NJ'
        b';~<17<6<gb`r;Lb3vJPzihIGBRa{%P#GyU0SjD>XmB=CMAZU6&neMNPIdN!gY>hw`#mDuQGI<Vyzm4~}Df=c|Hp$!~m%R?GOXrg7'
        b'8dOupvK)Rt`6EANHcNkcCP-dvTB>J?za{cIrjK=L47WM;N*>$pObD8tN(fb^YH_LOkF{xb;X)Iy_5dE<9d}gzBD`+i|KfOO4-6;-'
        b'2U)#Q;`w7)X_8GT)Z_txLg4V);JQDe^ShFfb{(&JBHA4@1!W&-YEAC?_1o-fSf;^Hn}FL-0F9QqZ;*JPQk?>$mkIJ6QXNSyCR#>v'
        b'!?ubVUbP^~YNC`7bw$FIS?AC$29oexgL0owb*bCZgXT-Xv!o~-vZmQ!c6_Y-733~}zD@k0%=1@~WkN8uOq#!63``1}ro?blFAJ(Z'
        b'{BLd^#8IQld7<hE=fOC{R-HhOOtGY&eVFC(QGMonkPj(5Ga-eQVn%WNn~W(GtQFQaPOt?1<iv@MQ}#O4q`z#&1@dc0HIuX6Y=T;P'
        b'dS#Rv3SGI~D_#ds)v>IsE7w%Pr+)?<xSsjmOe0BS<UE82jmY&pudSU(koXWiKwB_1u1?6hnt$#nceM({+mdKdj{}8#!}8a`Cf~EP'
        b'{d4!>J97-$Z)?3{YPp4580NgEgJpsecPuwR+lx+TOS|ohpCS<%1z(gD66Ma|Xwb9w9!Lk>I2RGR4%aOcsB|(z`ThweM;sxQ-$A-3'
        b'{h~`3Xy$;$O;&_r7Yf1z=FJR1UG)hSl&YBSIQbNE-z!=C@&+_X%w{uJoj}EpbsSXsS<RhCxsLR0HSssBY&Psgq0Lj(O~Gu1`NvdC'
        b'7sD4Lqh0y9;=MoT8Z~YQpcPvA!Pn6wTe=iar>1NsZ$4QST|pP8B4$A7A2v7W&xEpIUk*PEp(#2uzOYzS|L)jvg>AP|7mdp(3t=}3'
        b'GowT5V{za5xrQDOt?Z-jG$U(K^00`*;?<l~BA!5vl^r%#WsDdmvJY#F4Mg)6CF+82L&$Eq_fC9{ma`#@2%25&(gaHIX;)AC3S$%%'
        b'#Ct%}`U(3gNEtA{vhh6lV|5NHw8;XFqqw<CNLNdS)>*kc$eL(0L6KE`{wB=#(l*XfTNz0!$9rJ}?Z{0$#{3S@c!9t-Vm?L9V7BaH'
        b'VZO7>q}u!<L;03qKkDjAf4SunA2e0d%Yezf3T|8^9}2?&Wp`9S1`3run0-9Wf{g3dPp-5=k3Qdr<}Td*!4XQ+&N**0|NVr$Z73L!'
        b'#cRHc5L*XscYHpKz(TtOEMk%&Fl|kfURViEpfWe!2|}=27><P!t(cATzVrpdhZ&-_B?b^AXU0v|&D)R}28kl)wlXg!e2wZ0e0*(v'
        b'Y}9WI!<V^pb&03R?6}iNxCW&Cl`6@!!A!8&(<nKP9M}-`6P22RQs+a1L7+)wwnL;MILjXaO3{g6SbtR`+p2J8na<P6x#qcJvvHPU'
        b'eRCQsnj^3etV3ZC0<oLtk4mbbuQ1FtI;2jRm*$QIH%}eQz^(OKyE%z@hLix5o*ttCFLC%q=h{3K6f`~ESyr$*<Z1qDQH`{g_?mAE'
        b'*^eKG@}Y6y5o{v93+0Ar-Vpl?(0J5?M`k*W$5iK8mz`f5akO*4d0CBw)5iXdQDu-dY_&EZ=GUVOgC4jIBc~#uE3krUF#m50Si3hf'
        b'CZpVpzRS&G<rIW58JMrzka$zw51jzy)t4wjUWpJ?qiMp5gcE>Qx3<-KsANP-Yl3|mb;jl7nAt3{j*2~GDvuG(n6VsDEY}ORW}H$a'
        b'+R&wjE(g?J!ZUoL%oVj*_b||LDj9JKDJpU*(eYm%DIa`Q>h;PCq0>)}CJes?1oE^rI<1MAIhBK9dOS}H&b`?T;Z_=PEzWT(a_esr'
        b'!6r<w(CbGT%KnbXLT_0)Djs<`<5=gPuqjWv^j{tz-D9dljZ&Qm3HoAik^4Qob}WHfaIQO^5mRzonM)NT3`LQelf-B##jSX$*6ms>'
        b'4o<kN<|1;wQOExG;g;NZ;9*+pa5MR!$JdsFxip5#zh4udxBh2U&b<ianD>V=`QgvrtD)ahB%ZTbl@I@6xS94*kC8LhBk|{63+ai*'
        b'hzD%`9kjDy2MT?8STEwuGdTdOYxX)7;x;!!-ogG&qJs}n354z}XS%pK;UaTItqDO4ONQ}gpglk5dHk|dAF<+GWQO+v*;~6yITO8m'
        b'qvu)s(hT2ewM3|ERNh_DYSPytr{rFVIg<<xmFqF1F#es2xz$#bSzJf=Jq>@4+ZgXAXE#ie!r9*0+hYh$XRSHZpiV`e{hXtl3dq}='
        b'%6+eK7=l|Kyw@oq89$f{8aON9%$(cJFJ9AjCGp0YBS3Qq>t6hB0I?#ZG*WR#8vha7u!6Cq4#~;)L|Qk>AQ={c8Wv~?dec4$er<f>'
        b'XBEfDqtX)<C2O%X0L-d4BUV`Eyv=QG>9~ebK2)Vlyw3qjLB*6@{6VrBc#Tkj3WQx&gGv&w{jK1gJa+xGE(yfTH8Mx_J38?wEv54{'
        b'*zKVPmEmBnDLz__p;NRLH0+<J%OikvM_9YL5yN>3#(JKAHgjW6W%?LEyY^A_6vjm;juF>JdG;=~GB7s&B~O>&eYK{DS6NVk_5yNP'
        b'%cY#;`siwpt<CV=<a{y4=v_!XAL_iEt9O(cL%JepI@IX!D*Y=s6j=s<_=vFOufb$at-YCLQtE_g*UOi|pRR?c)?nzwil!)rTwnm#'
        b'Ns2pNbuwDc76WBSg-P;WgJdO4jKc7lryT+bIJ<5}PZc*F<PhY9h|nPoiM5O24nojmNoEb2y=_L<F$C03{5{v|8Dy^YM=6V$o+VDN'
        b'|9X;YSfJhwBjh`SMvS7InmevCGp7y6LH|O2%AYe_QYbIECl8Cc5RBJ!X;#o3=Sag2;yJqVvKlO+9w221cPn1JgsEk_n!)C?xB3o%'
        b'OvCR&`bJBmAaYGtHDD@nBkM;kpS-o&p1f>Sr3Iy^rJ4|H3vCI;N7Wt`V9p93(`5_o#wf8qOwg+ab;P)fWtm9i%&JCB*6f&;0*Uo_'
        b'=)MW))ae8PggR-&2O3b9-Dt8r%MB>6B#})wYjJCieWiZciY)PLXOZWGSvR{ui`g&Hb{83tg}wk|bkDisr%5d#xRfCa@eYWc;KEg7'
        b'%E4MV$TC>6&$>L4bS%=aJ?nQJDK@2jqh76brMY-pHM10XRD!AWT=oqaP*2UY*P#@PR<U&hLg9AXe$bPLQg<6=J~(NpxzNce2MJN('
        b'oG6r-aas{-_cR=d6!1o#>yW0E3T-{&iSHu?y}$S?G`CvF&O?9Q*24B#BR_#MSbK}~ZRua`(W3HYoVNsZ2Si-Zjk#XVW=rFBOUhk#'
        b'Tv4JubR>i6exbP3lZ<SUgfX*@H0?>Gq|V1|D;M`wHQ2eyV>+uT$h;2A-Ueewr43|(KfHjKl<k2qzO|AuiwvhH#E~xc$%{~=IZnJd'
        b'i9B}gw+pChFSEKRT=75u#^@!))}n+LtB!s`rSnW6v-?)2dDUM(e#PT=n-TX#G7PnmWJTtcJUh$%%-y+_jv{vP`2T#~BdaFe>nYeQ'
        b'GHgP7>l<}=$U|*!_0ZJ}wiR`f8Fp@Vpb`nUMM^s)DZ0AE(bh^2P-m<N3TNxqSvdR4Hw4WqD1kl7T*<_Um}lid@JN0=vs^;bf>)|z'
        b'*;!#hdKg(|@pR6+9et9@ARVhODTh(mBSB4M5@eCdWJ#kUjQ?E=YW2WlfV^5t0|1?7v*Ft<VFQ>5!#(cb)#ck!lNOty<ykP6tDga1'
        b'EJ5N*eFxFyFt#QBj*I-??x)^aLN}e`xH>q0VshZsyb%=9B+(Ojsb?Z9T}^Y^fUqBey1%_haZkuVSSJV_ceffS=-tX%JMqmpU=?cu'
        b'gFp*%{etU<$X+sh7zBn87ka45dvUk)XX7)LbHc~&KTu>jg~=N&R`c5uJ0%<HT`I!P7d^+*A{QF?*#B2ZXmFW%;x9sW1O3p6AEtt`'
        b'=BU-{Q_g-xdlV7#gF*4n$NNwcbzq9=rtQxozi^l^Uk65@voYj31bS5txBPdOdk1y#{9oCaT44gxhC={`w2`k+2i`yf7=Xxdw$}u5'
        b'NJ6duZM$B|fBU)WfW>%FsQZjgWs~9<?RC<;2i{*+Xn22A1v|iK2(9wdywpUP@nIlqcTuYa?K{p$^+)b!XvkpGO*Pm<V#S9A*>)a|'
        b'cIVP4vY!a^<X8ToM^YOcT|UVdy4LhCyJ)SlCAqFPO$hn^I#lWx64u3g;adXF`)tJ$cGJ+7(gZX-2R7qLFvSnSQt)?^U>#Vi0pyGZ'
        b'X(mVQHdKzn{~pYdz48nJkENkrg<>CVBVG~_)C&}Rkq7?mn_meTnm`kD)Dvk4U77-5JDFX01i^=tTF7vz-_^4pSe#VAn5VL)4Z3DW'
        b'd9QGd2YKFCkqqX4Q40v#dhJzEhf>rnUmEO3AqM7KZeyXNo0G`fhEq5HAr5-Kp->__Y^a@zwMKo8b^y>#gz>xZ`yWl^Pw%<0Yqh&7'
        b'$Hjl8iIqT!Wq4HryBivv>+#VCUQD);syK%5<E|&|fHk1mX|VQfhyLx6<r`l~OE#*auuaL7lvx*@S#jF?bA^PY*K8&ap)Unh_+>sM'
        b'Z0rrZ!v!$MZPz$|;*X^kd1U9qA5v<`o40tH*0H0~eQLLh2l_f>f;2xWAjwJloRBhRP@#CEJ{%n2D!$*s-DSGvOkQORM9f}ghok2T'
        b'Q{1a;b<c@1Xuj_6^YN8R*)PC)HEFOPZ3I!Sf^Ey|`3n1T-RZGwGS_Hi|1huHOzq0;js1p6Qg{Qiv8uM(>j=l$`(Y>aTvmROWo`We'
        b'aWrjMBgc795a{v(EqHMPK;xw&eu%&`koO*hr(dZGFqJCe6Tz63(2}AEPT=anhH6v;wO*Gc4?Rub*O^EA$$%6*MKUnMoC>G7>Za@*'
        b'3h=Xz@*Mqa0L5b2qIT@@u2YzqBKKq=P|;K1mjK~`b_L(avTz{HU(Fky`1uC{^jeMIMm`F#O9<RJyz<|kgmSEyHN{SF|4`l}$0-75'
        b'&WnH%a>hG!c<$Vc3<#Z_NtbfT0y|%DTv$yxNtRb5PTE7Q=m`!{#~N`rywC4`tw(H3*lAHeQz+%W3pAat$yy0GVF-q6S)94I=vQ=&'
        b')$>^e-|5rSmrt^E6LM>Z)r`9)e^10ZlMQk~Mfl4CX6d02t3BiQ8t>~i<jm@x)Q}|Oi0uap0F?9QoY|VTL|kxb3jkHT!{w5#o=EX>'
        b'z%#b5fH94lOk+tSXar8dsSjb!Mz>9QCwd2((Z{F5ZtbhuS^#|%O9qzTTW}=a75*ykw=V=u;Fbz&aAB}@2yb4U7ja?6ZZL5rwGXB|'
        b'x^|Oj`5H=(hzR$n+a(@dh92LU77+{I*R}y9;q~hB>{nLyM%tESevuB~>_tHu%rQL5B=%d|T(<_<bejBb1NBhIko~sKxcJj_q;k+Q'
        b'ye2aTnT4MAL!%k@lb-zP%a0H|5qZ(nQ4O-#%qa5A#%(Y!&Pv>sOg8G=Rp}0cO=b4loO(n1k7x07buXW;JY)PMEo9F^68<uf0MQwa'
        b'Ylwn=F~IR?-I~$~gN`Mm*{st{lgLo;$R$xj{7Hf~#a?qNwn5}-x<KjHP1`oo;elwjZ`7zz0p$>J*6F8n{GFD#myZIEuu)Q!^q>z?'
        b'xa~D^*8=*(!SF%(q}dGJ1Hujnwh)d7uzfRM-6tG#ShqACfw{xmFY#D=D48Ak43B#lbCPiac@fypk?6mCVj->s0HJlLww2A_knUwA'
        b'6Q$od1=l2gXeYv=`t6bXOiMads4TOT+^YdIq4hSB%k0<okpX24=&ihti;5c0aZ8%8tJ@iN?F}|s**`STOdgnf%a!mOU!i)O1HaX6'
        b'=cFli?y3dX=BY*rC5O51T33GJCiE+Uzc*_%R8c(!@Jcr0hzIo)^6M33A$o}NMkUev(;r_Yn@$)_$(?dLB2<Y2YQ5JIRIkA0!K01v'
        b'GzxJw+=ZNQ={S&CY-jR43oV`MxT(7LA)3+>MMsRyGzSxe%0!2wM<k?EHZtn%MMAFnEruI0tj@@us&a_zB5<2x2?&TgdHJqz0d@FJ'
        b'wu`g;Go0$WN_QW+-r6x&Tt03cvC${mX;ToVLy%Lep=5qe85V1*Xb(%?Hd9l6oy{e2TXCD1jo5`f&G;Mo{bmYayn#H@$g_LH@K*Dv'
        b'paau^?a>xnj_7?A46Rf~Ct1cr_@MB4WiO+?U4CU(ZOD)LcS6UvukH2O_>9zXS;XpPr3T-iE-)>u;QJ96n{vNTF<Z<0`X}MHsFp)U'
        b'eJio|EJi1{a@gYC#^5!&jy#RDn+4Q$>B^LQPUDFhNrwns+?3A&L6CA9x#O}^ettUFv!fo5KCp3LB5Xk3&(sC+v4`1p--<W+yC;%0'
        b'8GN~!WP4R~e&n^2pgd8t44G0k{UV<-UdsOg+Y&b1DRn*vEX8(EvaA*lBKh&z+O7xz;G`}aRoT`s@-V|Bk+Q(900&co^peMID7TH8'
        b'T|bj}&RW;56CJG8wkwu~i*`9_Jr#>Z%P4roE$_!ohwSpOPF7X&F}@zrZsFEC7GYIV!`DnmHVOoW9??9vuyn=q%r$zFom}hpose}_'
        b'VN1+>+M8h&XYr#Mqr4gJQ`h?b9zu@7L0*>Gd*zH01(}YO?LMI(K4NUVrzVChf4?2cRoLtc+ghs4&9IdsYaiRR^~G7CfUKD*qx0eW'
        b'9q&IpVJz=1?J$)Eu~+X{$}%N*@ubjEB_L%Kr%V+!RhGw?EkAiICR5ewnknNj)&g7{>MAD6AOQ&M_7CqK-PY?<moLoW#RcMG40!E_'
        b'zKO@rVi~GxpOFf8<(%7!CIJ}`t;FOfy}gRYm+)lG@>8k@X@ReAFLi0t$ti*b*Cj*Ef;KZ3^7;o3K7ywNK^m^JRZ4XB|HWd;BHFZt'
        b'LOAfL1mNGADxEf`tvsxMfox@VE_?}@uux1?G0L)?oPofcTfgrXmklsdV`(ZDriGF9f^s3KentD#-4S1O%v|tCH?9X&11B)K4?=HQ'
        b'2%x|-E8$o7$P1c+i=>VR-=i^ylG;=VsUV8S3P(8l#Y>!%$MJ*ee8ickW*>Y@!WFZ7?cz>-R*6i=Wgt58a@a@M>esftPUa>}Z0(AK'
        b'4O!PTbNb5v$%lw96OcHpjAyXrBEIjPq54Qg0>X^gh{2%aeU?isF9%`#mObAjLeNA8Twq!&RK__Qw64IfOcs7s>>8(y3k9(nS};OB'
        b'+^rf**@i{TYE$OtxLi(Rb^tvIQ{q?Z*YeUe9d^Xpz#BPkPP;AvR?J`Yc>pg5$x+*!wpdz0Pwr4t7a}7M<5M?Sar|rp#e$~m!ZU*%'
        b'9|JdQWl4u>6DLg-x);D}Xwk=bj5!Dl-n@o=k1o6!K*zO2Gj}SrgLJR7Cva9ZZcFqM<B=!`P1_ocG&ODe<CDZMo)wG92?k-POlr4S'
        b'oA#yvr8X;B9MW=gFuP8xco2KGRHAW>fT|JrQxvolXGe7xU}#Ngj1>{Vvhuq}Aa=}sw;oSV6P!CQ!FrU);26-Th^hhHigDT0?vVmj'
        b'&jPc#4+=v&u2})i3N9*AVPTDtVutmAGvT7gWsaZi6Amh~-g?tP{;0!u7U=1TMe^goH#x9OzA`Mts1kjw;XK+L;GWcz0-T^&oWNX4'
        b'vyFo8SVxiJ-J|(NfnU+nO}OQEp#Bsuf=u{5$C~+fOG2L;C23fQC`)IA1Sj{(G;P15l&dU0t<sm&9jJ9w`g+&rWto#>`Zy7sa<rcR'
        b'CfWi9IxCr13fTzr$>t-IX^hI6$t43os+oPdrwq}VKJg6h*paXwZT*#cCDF!&)DB4uq!sf<TJjo<HToltazqh~AsUNxWj9ItnL3Tg'
        b'%Jddz;C@HRsZx7ujzGB{vR}e_;ql#o_0E}&;~%aR)hH;<Kl;l!zQE;GE;q#a-<rT99}38|GAj{#3P8sA4D+e~T*y|CtiT-DsW{Ir'
        b'+|mfNr~F?u&_MxCgGP36t1kH1vD%J1rbiwF<jIk~qP`PpGmT-#L-jm^+UTP}_x&puMFTr0@w*qEoE31VfPU+Z|GQ(gjf*LkO5w~s'
        b'Jg)<s4JW$vyugo&WgDQ3X>qY`^-!e#`5tL*!p*D9JUh~Ip<TBd$3W_v-(x)y(s@P`L|!5znLo_~C_bjYBhuZ(bohWVX7nA1@KWx2'
        b'PIS{6o`ZhF*0Acj3>*L_XY#4)P{uxWbSd+U?*KUJOdP<T>v?E&qFJ%w(ttas3uoIP_T^poiSE}ukZfqalR4mc+xu;J%WkJ3oP8Wx'
        b'lYhtyFg=A;UG`-ojPkUO+hsb|wERW!T=lGl9H&c~Pdg{a@uw0+$E;`?hWDmj!WAzGV(^qhv3j^~nju@Hq(Z17quwI8Q-}!o+j?XN'
        b'z;e<+Y6sb12@mr;Cae79^h}UeVb)#Z>0q9z*0#v#kE{@0z4JSCR}H^iz-9^W=#}c@P?V;6)IEbo4p>k<{yvBt^Go2abZ3fn1oJEG'
        b'XVkDCY8K`y&GnN#a`!4c?Y2;1Z%A)YoH-vDH5`U6LurTn(BNV6)DB-&9u%6>>-pDhnbLVd^<p*!YTHoFroh<^@$e`dK4B~96iEaV'
        b'Dxc~;GT=A`OAPyimMp$^*gukQkI!R%sVhQc8;Lv|@BMvENs`C2IzO-%u>XQGR@%)S59BZY6r$@>6W<%ZL@$pt5<qeXs3yM3QoRja'
        b'V-?r-est+y49{bYbCB4G06F%pV3-lkg=6v59cJ5;{C*!nWmnlX-u0)p$r>Da_H@BtN~!V9uPe}5nzwBRUesM9`2YPQn9vgbzVH|5'
        b'x+N@qeh$02Hs#p1_rqfIXmJf}WaqzqWhMP@QXp!CatTduLr&4g%6&=27orVf>h!!x-|h0~{xO47q62jJQDN%QWZtnwr~fc>--hRU'
        b'`Ia58ft5`#L)d!{)(0d!t|sMNKA<&_tB&@gb3keQu_b+*daHJDRexDneN$nzVJBf{_EC$KdE=_9H7cJbfD|Q-C00u{S~3x5Sd<2j'
        b'RFXXcHS`a2m^KgFOyZVFVY11uTj&Qi?02CfTBU58pZN_N9(2Bm=C^{u56x8`6ffWG2rk1gS8S;V4ZBKTsG|L-<K+f=Zb!9BBurEs'
        b'LfRgsa{34}1D8gcM_p{eR)Q!-7id<x*1a;jNQ}WT?*1X;O^#1$?DT>x7a&yMwReKSp{+N<lC6R+J!|$u{xqOH4&!GxD;?hMKcc?%'
        b't(qrA-1cCGM}&bmkn)gckp!&+-oB0`Fz7E&zJ%vBh_Gu9hEN%PFoxY4D2Z`4sOa3}VUw=$BM*1p`Xglvf7BB~t)Bf~)~+*8?_njp'
        b'PO1>c9p6D&G62J==<ydh+BXo1*FY$%pLLCl(_$@f!nG?V2n45+A*5OQ2-yMwJ{P+ekm901c35SLmYK@+)C()S8m&c$JXmC#u)*U8'
        b'Y!TYsx4m6A=(R}01)b8Nh65b_CL$Pfz-5P&<kXyUO?e=7nw$;iRjo*z9xF1}U+O2zJL;MUCXZko!XOQosLSjXPWM98Ef;gp=`+#Z'
        b'U!N<0Es*4PAueM1!w1QS`u|FBJP_7&!62xenL?#Nshi))N@CklO(9ik8cu99@>8JSf?w3P#wprj`Zox+E$;B?ylV8D<ui3vLr<0Z'
        b'eeI^{r~V-=IZP$8SrDbW724NKA+gX`V<iqOjB?};es-QPi&WR%jwu%M)xXld!It3K&zh_E+I|RKUGJ|qJb7w0uq)FN&kQmmn!J0%'
        b'2qkBTDdx*5ISAi*IUo!?ij1e4(Q&`me_Ux7hQqrw(~|i8qx1S`IxRf$4F(sS{e67Lt?2^37a<ZyCys^Z52IvJ3h`Cvdb+-3ODA}C'
        b'T7NDB7xS&+g1KT>rE)lKr70h%7oBJnqM1wo2H~#LW!Z}j<l*U6^Bb1pn{7XRYETay9a_~6v8_Lyro9(CzO#-M9L^2-3v}|GRQ8s1'
        b'#dbvAv$+aNAx1$X#NN3NX)AGGe`yN=%zhZiMrLgjS<RwLoZTEM%BB9-F2`ylyE%4wY$=&AAu=o1Qr=a$HILK)eUf#MgG-<;@~yQX'
        b'mN;{70pe>|dz$CYm)^3>N3sYHuc-3;=4H-ui6A-}1dkkZLdzM=<A!0QZ3tYU^w@GSk~@XdS~VI6A+M;1|8eyM?7!ymTDE~Nhz5<`'
        b'GFYz=Fpij~NxfD!?#V_y)iuuM)<+WSZf9m}g^!Ox`-sGx14vP`92Rg%;j5xtZrg=VE*|W+odVdf^xN}KchnaTsq`(Cm~E3JEP44a'
        b'<GK3m^oL=eLn5V=gGZWZ8jab)L87!~&|pnGjeYS#7)rv8VjdDXh?{JIf>7%ICQJ?QWJ(HI7y#VyRRO3bi_B<sELr^wL%Gs?Qy>x~'
        b'IZwB!AWOc`6$!pJ)`&BfAHS(y{(MQKlEEhil3MIm?Sjx#le7Mv>|R5zwt(Q!=QZcOIZ?!`8l#azrLu1;rpcuBzcF8@`E_QU{0NyW'
        b'A$nb(zf9VMpl^3-98Hjq%}E8`Qjz%pw9f^efT2Q}00E>#h{gc`%RGNbvBYQl0ssI200dcD'
    ),
}
