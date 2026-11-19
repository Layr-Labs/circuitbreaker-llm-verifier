"""Distance-seven and outer threshold arithmetic."""
from types import SimpleNamespace


"""Standalone exact D7 contraction with regenerated size-zero/two arm inputs."""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import factorial,isqrt,lcm
from itertools import permutations,product
import json

U,N,SPATIAL=10**18,36,80
STEPS=tuple(range(1,24))
FAC7=factorial(7)
C_HI,Q_HI=F(8579,50000),F(70711,100000)
SQRT2_LO,SQRT2_HI=F(141421,100000),F(70711,50000)
ARMS_DIGEST="8e900c99d144cf8f46deb165bf2844b6d47735fd5806b27632014ce8df9b6ae9"


def digest(value):
    return sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def exp_lower(x):
    assert 0<=x<=64
    n,d=x.numerator,64*x.denominator
    value=0
    for j in range(7,-1,-1):
        value=value*(-n)+(FAC7//factorial(j))*d**(7-j)
    value=max(0,U*value//(FAC7*d**7))
    for _ in range(6):value=value*value//U
    return value


CENTRAL_DEN=lcm(*[2**j*factorial(j)*(2*j+1)for j in range(24)])
CENTRAL_COEFF=[(-1)**j*(CENTRAL_DEN//(2**j*factorial(j)*(2*j+1)))*1000**(2*(23-j))for j in range(24)]
CENTRAL_DIV=10000*CENTRAL_DEN*1000**47


@lru_cache(None)
def central_strip(index):
    t=min(index,3300);assert t>=0;square=t*t;value=0
    for c in reversed(CENTRAL_COEFF):value=value*square+c
    return max(0,7978*U*t*value//CENTRAL_DIV)


def strip(index,s,sd,divisor=1):
    value=index*index*1000**2*(s.denominator-s.numerator)*sd.denominator**2//((SPATIAL*divisor)**2*s.numerator*sd.numerator**2)
    return central_strip(isqrt(value))


def scalar_controls():
    def atan(x,last):
        return sum(((-1)**j*x**(2*j+1)/(2*j+1)for j in range(last+1)),F(0))
    pi_lo=16*atan(F(1,5),15)-4*atan(F(1,239),14)
    pi_hi=16*atan(F(1,5),14)-4*atan(F(1,239),15)
    assert F(31415,10000)<pi_lo<pi_hi<F(3927,1250)
    assert SQRT2_LO**2<2<SQRT2_HI**2 and Q_HI**2>F(1,2)
    assert C_HI==3-2*SQRT2_LO
    assert pi_hi*F(27638,10000)**2<24
    assert pi_hi*F(13819,10000)**2<6
    assert pi_hi*F(7978,10000)**2<2
    assert 0<exp_lower(F(1))<3*U//8
    assert central_strip(0)==0 and 0<central_strip(1000)<central_strip(2000)<U
    return pi_hi


def arm_inputs():
    """Recompute all 23 pairs of length-64 arrays, without a receipt input."""
    cells,central=48,64
    def fd(x):return 13819*exp_lower(F(3,2)*x*x)//10000
    @lru_cache(None)
    def fa(index):return 27638*exp_lower(6*F(index,1000)**2)//10000
    def ceil_index(q):return (1000*q.numerator+q.denominator-1)//q.denominator
    roots=[[2*fd(F(x+1+2*(a+1),SPATIAL))//SPATIAL for a in range(cells)]for x in range(central)]
    first=[[fa(ceil_index((F(b+1)+SQRT2_HI*(a+1))/SPATIAL))//SPATIAL for b in range(cells)]for a in range(cells)]
    lower_m=[[int(2*SQRT2_LO*b+a)for b in range(cells)]for a in range(cells)]
    zmax=max(central-1,max(map(max,lower_m)))
    dz=[0]+[fd(F(z,SPATIAL))//SPATIAL for z in range(1,zmax+1)]
    ell=[F(0)]+[C_HI/2*sum(((j+1)*Q_HI**j for j in range(d)),F(0))for d in range(1,4)]
    rows=[]
    for step in STEPS:
        s=F(step,24)**2
        root_strip=[strip(a,s,ell[1])for a in range(cells)]
        child_strip=[strip(b,s,ell[2])for b in range(cells)]
        terminal_strip=[strip(z,s,ell[3],2)for z in range(zmax+1)]
        base_strip=[strip(z,s,C_HI)for z in range(zmax+1)]
        terminal=[sum(dz[z]*terminal_strip[m-z]for z in range(1,m+1))//U for m in range(zmax+1)]
        baseline=[sum(dz[z]*base_strip[x-z]for z in range(1,x+1))//U for x in range(central)]
        inner=[sum((first[a][b]*child_strip[b]//U)*terminal[lower_m[a][b]]for b in range(cells))//U for a in range(cells)]
        arms=[sum((roots[x][a]*root_strip[a]//U)*inner[a]for a in range(cells))//U for x in range(central)]
        assert all(a<=b for a,b in zip(baseline,baseline[1:]))
        assert all(a>=b for a,b in zip(arms,arms[1:]))
        assert all(0<=v<=U for v in baseline+arms)
        rows.append({'step':step,'s':[s.numerator,s.denominator],'baseline_scaled':baseline,'size2_arm_scaled':arms})
    assert digest(rows)==ARMS_DIGEST
    return rows


EXP_COEFF=[FAC7//factorial(j)*64000**(7-j)for j in range(8)]
EXP_DEN=FAC7*64000**7


@lru_cache(None)
def exp_index(index):
    assert 0<=index<=64000
    value=0
    for c in reversed(EXP_COEFF):value=value*(-index)+c
    value=max(0,U*value//EXP_DEN)
    for _ in range(6):value=value*value//U
    return value


def margin(a,b,c=0):return (282842*b+100000*(a+c))//100000


def geometry():
    initial=[]
    for i in range(N):
        row=[];I=2*i+1
        for j in range(N):
            J=2*j+1;n=1800000*I*I+1950000*J*J+2545596*I*J
            row.append((exp_index((n+1250000+2559999)//2560000),margin(j,i)))
        initial.append(row)
    inner=[];last=[];maximum=0
    for b in range(N):
        B=2*b+1
        for c in range(N):
            C=2*c+1;middle=[];final=[]
            for a in range(N):
                A=2*a+1;m=margin(a,b,c);maximum=max(maximum,m)
                n=1950000*C*C+2545596*B*C+600000*A*C
                middle.append((exp_index((n+650000+2559999)//2560000),m))
                n=1800000*C*C+2545596*B*C+600000*A*C
                final.append((exp_index((n+600000+2559999)//2560000),m))
            inner.append(middle);last.append(final)
    return initial,inner,last,maximum


from math import comb
from operator import mul

PARAMETERS = ((7, F(1997, 50), F(333, 2000), F(1129, 1000), F(1, 2500000), 532), (8, F(281, 5), F(327, 2000), F(553, 500), F(19, 250000000), 1040), (9, F(783, 10), F(101, 625), F(1091, 1000), F(7, 500000000), 1968), (10, F(546, 5), F(803, 5000), F(537, 500), F(53, 20000000000), 3620))
PATHS = PARAMETERS[1:]
PANELS = (2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 17, 20, 23)

def translate(polynomial, center):
    return tuple(sum(a * comb(i, j) * center**(i-j)
                     for i, a in enumerate(polynomial) if i >= j)
                 for j in range(len(polynomial)))

def longer_controls():
    pi_upper = scalar_controls()
    p8, p9, p10 = (24, -104, 70, -15, 1), (112, -200, 96, -17, 1), (-48, 320, -340, 126, -19, 1)
    assert translate(p8, 8) == (88, 184, 94, 17, 1)
    assert translate(p9, 8) == (48, 120, 72, 15, 1)
    assert translate(p10, 8) == (208, 640, 508, 158, 21, 1)
    assert 4 * 94 > 17**2 and 4 * 72 > 15**2
    assert translate((508, 158, 21, 1), -8) == (76, 14, -3, 1) and 4 * 14 > 3**2
    center = F(9, 2)
    negative_side_bound = 3 * sum(abs(a)*center**i for i, a in enumerate(translate(p9, center)))
    assert negative_side_bound < 1920*2*SQRT2_LO
    assert 384 * F(11, 48) == 88 and 1536 * F(13, 96) == 208
    assert 768*F(1, 8) == 48*2
    assert 736 * SQRT2_LO > PATHS[0][-1]
    assert 48 + 120 * 16 == PATHS[1][-1]
    assert 2560 * SQRT2_LO > PATHS[2][-1]
    scale = 10**12
    sqrt10_lower = F(isqrt(10 * scale**2), scale)
    assert sqrt10_lower**2 < 10
    theta_upper = 38 - 12 * sqrt10_lower
    determinants = (F(263, 288), 1-F(11, 48)**2, F(31, 32), 1-F(13, 96)**2)
    correlation_upper = (5*SQRT2_HI/24, F(11, 48), SQRT2_HI/8, F(13, 96))
    assert 7*8**3-65*8**2+144*8-44 == PARAMETERS[0][-1]
    assert 280*2*SQRT2_LO-708 > 0
    controls = []
    for params, endpoints_det, w in zip(PARAMETERS, determinants, correlation_upper, strict=True):
        distance, normalizer, density, mean, _, multiplier = params
        middle = distance-1
        middle_det = 9 * 12**middle * endpoints_det
        assert normalizer**2 * (2*pi_upper)**middle < 4 * middle_det
        assert (2*pi_upper*density)**2 * endpoints_det < 1
        assert 1 + F(2*(middle//2), middle)*theta_upper + 2*w/(middle*(1+w)) < mean
        controls.append({"distance": distance, "endpoint_determinant": str(endpoints_det),
                         "conditional_middle_precision_determinant": str(middle_det),
                         "orientation_counted_normalizer": str(normalizer), "endpoint_pair_density": str(density),
                         "middle_variance_mean": str(mean), "spectral_multiplier": multiplier})
    return controls

def incidence_time_moment(left, right, modes, mean):
    a, b = 1-left, 1-right
    degree = (modes+4)//2
    den = 1+(mean-1)*right
    if modes % 2:
        lo = F(isqrt(a.numerator*U**2//a.denominator), U)
        hi = F(isqrt(b.numerator*U**2//b.denominator)+1, U)
        moment = (a**degree*lo-b**degree*hi)/F(2*degree+1, 2)
        reciprocal = 1/den**modes
        determinant = F(isqrt(reciprocal.numerator*U**2//reciprocal.denominator), U)
    else:
        moment = (a**degree-b**degree)/degree
        determinant = 1/den**(modes//2)
    tangent = 1+F(modes, 4)*(mean-1)*(right-left)/den
    assert moment > 0
    return moment*determinant*tangent

def advance(state, weights):
    return [[sum(map(mul, column, coefficients))//U for coefficients in row]
            for column, row in zip(zip(*state, strict=True), weights, strict=True)]

def panel_masses(s, arm, shape):
    initial, interior, last, maximum = shape
    sides = [b+(s.denominator-s.numerator)*a//s.denominator
             for b, a in zip(arm["baseline_scaled"][:N], arm["size2_arm_scaled"][:N], strict=True)]
    strips = [strip(m, s, C_HI) for m in range(maximum+1)]
    inner, end = (
        [[[ex*strips[m]//U//SPATIAL for ex, m in row] for row in factors[i:i+N]]
         for i in range(0, N*N, N)]
        for factors in (interior, last)
    )
    state = [[ex*sides[i]//U*strips[m]//U//SPATIAL**2 for ex, m in row]
             for i, row in enumerate(initial)]
    for _ in range(3):
        state = advance(state, inner)
    tail = [[strips[margin(i, j)] for j in range(N)] for i in range(N)]
    for distance in range(7, 11):
        final = advance(state, end)
        mass = sum(sum(v*y//U*z//U for v, y, z in zip(row, sides, t, strict=True))
                   for row, t in zip(final, tail, strict=True))
        yield mass
        if distance < 10:
            state = advance(state, inner)

def certificate():
    controls = longer_controls()
    arms, shape = arm_inputs(), geometry()
    totals = [F(0)]*len(PARAMETERS)
    previous = F(0)
    for step in PANELS:
        s = F(step, 24)**2
        for i, (p, mass) in enumerate(zip(PARAMETERS, panel_masses(s, arms[step-1], shape), strict=True)):
            distance, normalizer, density, mean, _, _ = p
            totals[i] += density*normalizer*F(mass, U)*incidence_time_moment(previous, s, distance-1, mean)
        previous = s
    rows = []
    for (distance, _, _, _, floor, multiplier), total in zip(PARAMETERS, totals, strict=True):
        assert total > floor
        rows.append({"distance": distance, "exact_raw_price_sum": str(total),
                     "usable_raw_price_floor": str(floor), "K_increment": str(floor*multiplier)})
    return {"scope": "exact D7–D10 shared Price contraction", "middle_cells": N, "time_steps": list(PANELS),
            "path_controls": controls, "rows": rows,
            "global_K_increment": str(sum(F(row["K_increment"]) for row in rows))}


d7 = SimpleNamespace(certificate=certificate, scalar_controls=scalar_controls, C_HI=C_HI, Q_HI=Q_HI, SQRT2_HI=SQRT2_HI, SQRT2_LO=SQRT2_LO, U=U, exp_index=exp_index, exp_lower=exp_lower, margin=margin, strip=strip)


STEP, T = F(1, 100), F(3, 100)
LO, HI, ROOT3 = SQRT2_LO, SQRT2_HI, F(173205, 100000)
PHI = F(797, 2000)


def moving_controls():
    def atan(x,last):return sum(((-1)**j*x**(2*j+1)/(2*j+1)for j in range(last+1)),F(0))
    pi_lo=16*atan(F(1,5),15)-4*atan(F(1,239),14)
    pi_hi=16*atan(F(1,5),14)-4*atan(F(1,239),15)
    assert F(31415,10000)<pi_lo<pi_hi<F(3927,1250)
    assert LO**2<2<HI**2 and ROOT3**2<3
    assert F(229,50)**2*pi_hi<66 and F(69,50)**2*pi_hi<6
    assert F(27638,10000)**2*pi_hi<24 and F(7978,10000)**2*pi_hi<2
    assert F(81027,10000)*pi_hi<18*LO
    endpoint,child,closing=F(1,125),F(13,5000),F(13,2500)
    assert 0<2-14*HI/11 and T*(2-14*LO/11)<endpoint
    assert 0<3-2*HI and T*(3-2*LO)/2<child
    assert T*(3-2*LO)<closing
    assert max(T*(8*HI/11-1),T/(22*LO),endpoint)<STEP
    assert (1+LO)*STEP>child and (1+2*LO)*STEP>closing
    bias=HI/2-F(2,3);c=3-2*LO
    assert 0<LO/2-F(2,3) and bias*T<STEP
    assert T*(F(1,3)-LO/6)<STEP
    assert T*(bias+c)<F(13,2000)<endpoint
    assert F(627,250)**2>2*pi_hi and (1-T*T/2)/F(627,250)>PHI
    assert F(33,200)<F(1,6)<F(17,100)
    assert F(54,100)**2<F(11,36)<F(56,100)**2
    assert (1-(F(17,100)*T/F(54,100))**2/2)/(F(56,100)*F(627,250))>F(71,100)
    assert F(73,250)<1-HI/2<1-LO/2<F(147,500)
    assert F(707,1000)<LO/2<HI/2<F(177,250)
    assert (1-(F(147,500)*T)**2/(2*F(707,1000)**2))/(F(177,250)*F(627,250))>F(14,25)
    return pi_hi


def closing_strip(xlow,ylow):
    numerator=282842*ylow+100000*xlow-52000
    assert numerator>0
    return central_strip(173205*numerator//1000000000)


BASELINE_FLOORS = ('241/5000', '419/100000')
BASELINE_EXACT = ('2411227882761611/50000000000000000', '2100371733848019/500000000000000000')
BASELINE_FLOORS = tuple(map(F, BASELINE_FLOORS))
BASELINE_EXACT = tuple(map(F, BASELINE_EXACT))
def strip_box(left,right,offset,squared=False):
    a,b=left-offset,right-offset
    assert 0<=a<b<=F(3,5)
    terms=((2,F(1)),(4,F(-1)),(6,F(1,4)))if squared else((1,F(1)),(3,F(-1,2)))
    moments=[sum(c*(b**(k+j+1)-a**(k+j+1))/(k+j+1)for k,c in terms)for j in range(3)]
    m0,m1,m2=moments;assert m0>0
    return F(69,50)**(2 if squared else 1)*m0,m1/m0+offset,m2/m0+2*offset*m1/m0+offset**2

def baseline_controls():
    def atan(x,last):return sum(((-1)**j*x**(2*j+1)/(2*j+1)for j in range(last+1)),F(0))
    pi_lo=16*atan(F(1,5),15)-4*atan(F(1,239),14)
    pi_hi=16*atan(F(1,5),14)-4*atan(F(1,239),15)
    assert F(31415,10000)<pi_lo<pi_hi<F(3927,1250)
    assert LO**2<2<HI**2 and F(69,50)**2*pi_hi<6
    assert F(627,250)**2>2*pi_hi
    assert (1-F(3,100)**2/2)/F(627,250)>F(797,2000)
    normalizers=F(229,50),F(81027,10000),F(128,10),F(19)
    norm2,norm3,norm4,norm5=normalizers
    assert norm2**2*pi_hi<66
    assert norm3*pi_hi<18*LO
    assert norm4**2*(2*pi_hi)**3<4*10260
    assert norm5**2*(2*pi_hi)**4<4*7*12**4
    return pi_hi,normalizers

def weighted_distances(normalizers):
    budgets=[];left=F(1,100)
    for distance,normalizer,usable in zip(range(2,4),normalizers[:2],BASELINE_FLOORS,strict=True):
        single=distance==2;n=distance-1
        step,offset=(F(49,1000),F(1,125))if single else(F(1,25),F(13,2000))
        sides=[];mids=[]
        for i in range(12):
            a,b=left+i*step,left+(i+1)*step
            sides.append(strip_box(a,b,offset,squared=single))
            mids.append((step,(a+b)/2,(a*a+a*b+b*b)/3))
        coordinates=(sides,)if single else(sides,)+(mids,)*(n-2)+(sides,)
        total=0
        for indices in product(range(12),repeat=n):
            if indices>indices[::-1]:continue
            boxes=[coordinate[i]for coordinate,i in zip(coordinates,indices)]
            if single:
                exponent=F(33,2)*boxes[0][2];volume=boxes[0][0]
            else:
                exponent=18*(boxes[0][2]+boxes[-1][2])+F(39,2)*sum(row[2]for row in boxes[1:-1])
                exponent+=18*HI*sum(a[1]*b[1]for a,b in zip(boxes,boxes[1:]))
                exponent+=6*sum(a[1]*b[1]for a,b in zip(boxes,boxes[2:]))
                volume=boxes[0][0]*boxes[-1][0]*step**(n-2)
            index=(exponent.numerator*1000+exponent.denominator-1)//exponent.denominator
            weight=normalizer*volume*exp_index(index)
            total+=(weight.numerator//weight.denominator)*(1 if indices==indices[::-1]else 2)
        exact=F(total,U)
        assert exact==BASELINE_EXACT[distance-2] and exact>usable
        budgets.append({'distance':distance,'exact_probability_lower_sum':str(exact),'usable_probability':str(usable),
                        'cells':12**n,'cell_width':str(step),'offset':str(offset)})
    return budgets


CELLS = 60
SHORT_FLOORS = tuple(F(i, 10**7) for i in (552889, 313397, 114197, 75179, 142689, 32000))
LONG_FLOORS = tuple(F(i, 10**9) for i in (4439266, 874902, 166833))
LONG_PARAMETERS = (
    (F(64, 5), F(95, 144), F(52, 125), F(49, 100)),
    (F(96, 5), F(7, 9), F(66, 125), F(9, 20)),
    (F(28), F(55, 64), F(78, 125), F(429, 1000)),
)


def ceil(numerator, denominator):
    return (numerator+denominator-1)//denominator


def child_weight(previous, current, following):
    # Input indices locate upper cell endpoints, with zero denoting an absent
    # predecessor. These are the doubled mean and twelve times the variance of
    # 50000*y + 70711*x + 25000*p + 13000 in cell-index coordinates.
    mean2 = 26000+70711*(2*current-1)+50000*(2*following-1)
    variance12 = 70711**2+50000**2
    if previous:
        mean2 += 25000*(2*previous-1)
        variance12 += 25000**2
    index = ceil(3*mean2*mean2+variance12, 50_000_000_000)
    return 27638*exp_index(index)//1_000_000


def root_weight(middle_upper, first):
    # Use the worst middle endpoint, averaging only the arm coordinate.
    mean2 = 10*middle_upper+20*first-2
    return 69*exp_index(ceil(3*mean2*mean2+100, 2000))//2500


def dot(left, right):
    return sum(a*b for a, b in zip(left, right, strict=True))//U


def moving_arm_inputs():
    endpoints = range(2, CELLS+2)
    closing = [[closing_strip(a-1, b-1) for b in endpoints] for a in endpoints]
    first = [[child_weight(0, a, b) for b in endpoints] for a in endpoints]
    two = [dot(row, strip) for row, strip in zip(first, closing, strict=True)]
    continuation = [[[child_weight(a, b, c) for c in endpoints] for b in endpoints]
                    for a in endpoints]
    tail = closing
    for _ in range(2):
        tail = [[dot(weights, tail[j]) for j, weights in enumerate(matrix)]
                for matrix in continuation]
    four = [dot(row, strip) for row, strip in zip(first, tail, strict=True)]
    two_sided, baseline, four_sided = [], [], []
    for x in range(2, CELLS+1):
        root = [root_weight(x, a) for a in endpoints]
        two_sided.append(dot(root, two))
        baseline.append(central_strip(173205*(5*(x-1)-4)//50000))
        four_sided.append(dot(root, four))
    return two_sided, baseline, four_sided


def short_paths(two, baseline, four):
    sums = [0]*6
    earlier = [a+b for a, b in zip(two, baseline, strict=True)]
    for i, x in enumerate(range(2, CELLS+1)):
        density2 = 229*exp_index(ceil(11*(3*x*x-3*x+1), 20))//5000
        sums[0] += density2*2*two[i]*baseline[i]//U**2
        sums[1] += density2*two[i]**2//U**2
        sums[4] += density2*(2*earlier[i]*four[i]+four[i]**2)//U**2
        for j in range(i+1):
            y = j+2
            squares = 3*(x*x+y*y-x-y)+2
            product = (2*x-1)*(2*y-1)
            exponent = ceil(18*(200000*squares+3*70711*product), 6_000_000)
            density3 = 81027*exp_index(exponent)//100_000_000
            mult = 1 if i == j else 2
            sums[2] += mult*(density3*(two[i]*baseline[j]+baseline[i]*two[j])//U**2)
            sums[3] += mult*(density3*two[i]*two[j]//U**2)
            extra = earlier[i]*four[j]+four[i]*earlier[j]+four[i]*four[j]
            sums[5] += mult*(density3*extra//U**2)
    return tuple(F(n, U) for n in sums)


def longer_paths(two, baseline, four, *, last_distance=6):
    # Same central moment factors as the earlier long path contraction.
    start = [a+b for a, b in zip(two, four, strict=True)]
    end = [a+2*b for a, b in zip(start, baseline, strict=True)]
    centers = tuple(2*i+3 for i in range(len(start)))

    def weight(previous, current, following, terminal=False):
        diagonal = 900000 if terminal else 975000
        numerator = (diagonal*(3*following*following+1)
                     +54*70711*current*following+900000*previous*following)
        return exp_index(ceil(numerator, 6_000_000))//100

    first = [a*(exp_index(ceil(9*x*x+3, 20))//100)//U
             for a, x in zip(start, centers, strict=True)]
    initial = [[a*weight(0, x, y)//U for y in centers] for a, x in zip(first, centers, strict=True)]
    terminal = [[dot((weight(x, y, z, True) for z in centers), end) for y in centers]
                for x in centers]
    continuation = [[[weight(x, y, z) for z in centers] for y in centers] for x in centers]
    sums = []
    for distance, (norm, _, _, _) in enumerate(LONG_PARAMETERS[:last_distance-3], 4):
        total = sum(a*b for row, tail in zip(initial, terminal, strict=True)
                    for a, b in zip(row, tail, strict=True))//U
        sums.append(norm*F(total, U))
        if distance < last_distance:
            terminal = [[dot(weights, terminal[j]) for j, weights in enumerate(matrix)]
                        for matrix in continuation]
    return tuple(sums)


def certificate(*, last_distance=6):
    assert last_distance in (4, 5, 6)
    pi_hi = moving_controls()
    arms = moving_arm_inputs()
    short, long = short_paths(*arms), longer_paths(*arms, last_distance=last_distance)
    long_floors = LONG_FLOORS[:last_distance-3]
    long_parameters = LONG_PARAMETERS[:last_distance-3]
    assert all(exact > floor > 0 for exact, floor in zip(short, SHORT_FLOORS, strict=True))
    assert all(exact > floor > 0 for exact, floor in zip(long, long_floors, strict=True))
    short_rate = (6*F(33, 200)*F(71, 100)*(sum(SHORT_FLOORS[:2])+SHORT_FLOORS[4])
                  +12*F(73, 250)*F(14, 25)*(sum(SHORT_FLOORS[2:4])+SHORT_FLOORS[5]))
    long_rate = F(0)
    for distance, ((norm, variance, drift, density), floor) in enumerate(
            zip(long_parameters, long_floors, strict=True), 4):
        assert norm**2*(2*pi_hi)**(distance-1) < 36*12**(distance-1)*variance
        long_rate += 3*2**(distance-1)*drift*density*floor
    return {"scope": f"replacement moving-cover arm probabilities at distances two through {last_distance}",
            "threshold_abs_max": str(T), "common_normal_density_floor": str(PHI),
            "short_exact": list(map(str, short)), "short_floors": list(map(str, SHORT_FLOORS)),
            "long_exact": list(map(str, long)), "long_floors": list(map(str, long_floors)),
            "arm_rate_without_normal_density": str(short_rate+long_rate)}


associated_arm_certificate = certificate
EXTENT = F(9, 400)
REQUIRED_CREDIT = F(1935851777, 10**14)


def certificate():
    _, norms = baseline_controls()
    baselines = weighted_distances(norms)
    arm = associated_arm_certificate(last_distance=6)
    assert [row["distance"] for row in baselines] == [2, 3]
    assert EXTENT <= F(arm["threshold_abs_max"])
    assert F(arm["common_normal_density_floor"]) == PHI
    rates = [multiplicity*drift*density*F(row["usable_probability"])
             for row, (multiplicity, drift, density) in zip(
                 baselines, ((6, F(33, 200), F(71, 100)), (12, F(73, 250), F(14, 25))),
                 strict=True)]
    arm_rate = F(arm["arm_rate_without_normal_density"])
    rate = sum(rates)+arm_rate
    saving = PHI*rate*EXTENT**2/2
    assert saving > REQUIRED_CREDIT
    return {"scope": "outer with only no-arm D2-D3 and associated D2-D6 arms",
            "threshold_abs_max_of_central_band": str(EXTENT), "baselines": baselines,
            "baseline_rates_without_normal_density": list(map(str, rates)),
            "arms": arm, "rate_without_normal_density": str(rate),
            "outer_graph_credit": str(saving),
            "slack_above_uniform_saving": str(saving-REQUIRED_CREDIT)}


outer = SimpleNamespace(certificate=certificate)
