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
    for params, endpoints_det, w in zip(PARAMETERS, determinants, correlation_upper, strict=True):
        distance, normalizer, density, mean, _, multiplier = params
        middle = distance-1
        middle_det = 9 * 12**middle * endpoints_det
        assert normalizer**2 * (2*pi_upper)**middle < 4 * middle_det
        assert (2*pi_upper*density)**2 * endpoints_det < 1
        assert 1 + F(2*(middle//2), middle)*theta_upper + 2*w/(middle*(1+w)) < mean


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
    longer_controls()
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
            "rows": rows,
            "global_K_increment": str(sum(F(row["K_increment"]) for row in rows))}


d7 = SimpleNamespace(certificate=certificate, scalar_controls=scalar_controls, C_HI=C_HI, Q_HI=Q_HI, SQRT2_HI=SQRT2_HI, SQRT2_LO=SQRT2_LO, U=U, exp_index=exp_index, exp_lower=exp_lower, margin=margin, strip=strip)


STEP, T = F(1, 100), F(9, 400)
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
    endpoint,child,closing=F(49,10000),F(1,500),F(39,10000)
    assert 0<2-14*HI/11 and T*(2-14*LO/11)<endpoint
    assert 0<3-2*HI and T*(3-2*LO)/2<child
    assert T*(3-2*LO)<closing
    assert max(T*(8*HI/11-1),T/(22*LO),endpoint)<STEP
    assert (1+LO)*STEP>child and (1+2*LO)*STEP>closing
    bias=HI/2-F(2,3);c=3-2*LO
    assert 0<LO/2-F(2,3) and bias*T<STEP
    assert T*(F(1,3)-LO/6)<STEP
    assert T*(bias+c)<endpoint
    assert F(627,250)**2>2*pi_hi and (1-T*T/2)/F(627,250)>PHI
    assert F(33,200)<F(1,6)<F(17,100)
    assert F(54,100)**2<F(11,36)<F(56,100)**2
    assert (1-(F(17,100)*T/F(54,100))**2/2)/(F(56,100)*F(627,250))>F(71,100)
    assert F(73,250)<1-HI/2<1-LO/2<F(147,500)
    assert F(707,1000)<LO/2<HI/2<F(177,250)
    assert (1-(F(147,500)*T)**2/(2*F(707,1000)**2))/(F(177,250)*F(627,250))>F(14,25)
    return pi_hi


def closing_strip(xlow,ylow):
    numerator=282842*ylow+100000*xlow-39000
    assert numerator>0
    return central_strip(173205*numerator//1000000000)


BASELINE_FLOORS = ('241/5000', '419/100000', '89/100000', '89/500000')
BASELINE_FLOORS = tuple(map(F, BASELINE_FLOORS))

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
    paired_limit = F(23, 50)
    assert 1-F(3, 2)*paired_limit**2 > 3*paired_limit**2
    for distance,normalizer,usable in zip(range(2,6),normalizers[:4],BASELINE_FLOORS,strict=True):
        single=distance==2;n=distance-1
        cells=12
        step = F(49,1000) if single else F(1,25)
        offset = F((46,49,49,48)[distance-2],10000)
        sides=[];mids=[]
        for i in range(cells):
            a,b=left+i*step,left+(i+1)*step
            # Pair the reflected distance-two signs where the squared strip is convex.
            paired = single and b+offset <= paired_limit
            sides.append(strip_box(a,b,F(0) if paired else offset,squared=single))
            mids.append((step,(a+b)/2,(a*a+a*b+b*b)/3))
        coordinates=(sides,)if single else(sides,)+(mids,)*(n-2)+(sides,)
        total=0
        for indices in product(range(cells),repeat=n):
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
        floor = F((exact*10**8).__floor__()-1,10**8)
        assert exact>floor>usable
        usable = floor
        budgets.append({'distance':distance,'exact_probability_lower_sum':str(exact),'usable_probability':str(usable),
                        'cells':cells**n,'cell_width':str(step),'offset':str(offset)})
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
    # 50000*y + 70711*x + 25000*p in cell-index coordinates.
    mean2 = 70711*(2*current-1)+50000*(2*following-1)
    variance12 = 70711**2+50000**2
    if previous:
        mean2 += 25000*(2*previous-1)
        variance12 += 25000**2
    index = ceil(3*mean2*mean2+variance12, 50_000_000_000)
    return 27638*exp_index(index)//1_000_000


def root_weight(middle_upper, first):
    # Use the worst middle endpoint, averaging only the arm coordinate.
    mean2 = 100*middle_upper+200*first-100
    return 69*exp_index(ceil(3*mean2*mean2+10000, 200000))//2500


def dot(left, right):
    return sum(map(mul, left, right))//U


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
    for _ in range(4):
        for _ in range(2):
            tail = [[dot(weights, tail[j]) for j, weights in enumerate(matrix)]
                    for matrix in continuation]
        for i, (row, strip) in enumerate(zip(first, tail, strict=True)):
            four[i] += dot(row, strip)
    two_sided, baseline, four_sided = [], [], []
    for x in range(2, CELLS+1):
        root = [root_weight(x, a) for a in endpoints]
        two_sided.append(dot(root, two))
        baseline.append(central_strip(173205*(100*(x-1)-49)//1000000))
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
    penalty = 2*(F(3,2)*F(49,10000)**2+11*6*F(1,500)**2)
    paired_factor = F(99939,100000)
    assert 0 < paired_factor < 1-penalty
    short = tuple(paired_factor*value for value in short_paths(*arms))
    long = tuple(paired_factor*value for value in longer_paths(*arms, last_distance=last_distance))
    long_floors = LONG_FLOORS[:last_distance-3]
    long_parameters = LONG_PARAMETERS[:last_distance-3]
    assert all(exact > floor > 0 for exact, floor in zip(short, SHORT_FLOORS, strict=True))
    assert all(exact > floor > 0 for exact, floor in zip(long, long_floors, strict=True))
    short_floors = tuple(F((value*10**9).__floor__()-1,10**9) for value in short)
    long_floors = tuple(F((value*10**9).__floor__()-1,10**9) for value in long)
    assert all(a>b for a,b in zip(short_floors, SHORT_FLOORS, strict=True))
    assert all(a>b for a,b in zip(long_floors, LONG_FLOORS[:last_distance-3], strict=True))
    for distance, ((norm, variance, drift, density), floor) in enumerate(
            zip(long_parameters, long_floors, strict=True), 4):
        assert norm**2*(2*pi_hi)**(distance-1) < 36*12**(distance-1)*variance
    return {
        "threshold_abs_max": str(T), "common_normal_density_floor": str(PHI),
        "short_floors": list(map(str, short_floors)), "long_floors": list(map(str, long_floors)),
    }


associated_arm_certificate = certificate
EXTENT = F(9, 400)
REQUIRED_CREDIT = F(28189565029999,1000000000000000000)

ADJACENT_UNIT = 10**18
ADJACENT_WIDTH = F(1,100)
ADJACENT_EXTENT = EXTENT


def packed_chords(prefix, strips, convolve):
    n = len(strips)-1
    lower = convolve(prefix, strips[:n], n)
    upper = convolve(prefix, strips[2:]+[0], n)
    return [(3*a+b, a+3*b) for a,b in zip(lower, upper, strict=True)]


def chord_cell(i, m0, m1, j0, j1, width_squared):
    dm, dj = m1-m0, j1-j0
    assert dm >= 0 and dj >= 0
    c0 = m0*j0*j0
    c1 = dm*j0*j0+2*m0*j0*dj
    c2 = 2*dm*j0*dj+m0*dj*dj
    c3 = dm*dj*dj
    moment = 60*c0+30*c1+20*c2+15*c3
    tilt = (30*i+15)*c1+(30*i+16)*c2+(27*i+15)*c3
    value = 6*width_squared*moment-33*tilt
    assert value >= 0
    return value//(360*width_squared*U**2)


def chord_joint_rows(densities, middle, old, high, width, gain):
    divisor = 4*U*width
    old_row = joint_row = 0
    for i, (density, (a,b), (c,d)) in enumerate(zip(densities, old, high, strict=True)):
        old0, old1 = a//divisor, b//divisor
        j0 = min(a*gain//(divisor*U), c//divisor)
        j1 = min(b*gain//(divisor*U), d//divisor)
        assert j0 >= old0 and j1 >= old1
        args = (i, middle[i], middle[i+1])
        old_row += density*chord_cell(*args, old0, old1, width**2)//(U*width)
        joint_row += density*chord_cell(*args, j0, j1, width**2)//(U*width)
    return old_row, joint_row


def chord_clean_row(densities, middle, chords, width):
    divisor = 4*U*width
    return sum(
        density*chord_cell(i, middle[i], middle[i+1], a//divisor, b//divisor, width**2)//(U*width)
        for i, (density, (a,b)) in enumerate(zip(densities, chords, strict=True))
    )


@lru_cache(maxsize=None)
def half_density(prefactor, rate, mean_square):
    base = 1-rate*mean_square/64
    assert 0 <= base <= 1
    return (ADJACENT_UNIT*prefactor.numerator*base.numerator**64
            // (prefactor.denominator*base.denominator**64))


HALF_CDF_TERMS = tuple(F(-3, 2)**i/(factorial(i)*(2*i+1)) for i in range(8))
HALF_CDF_DENOMINATOR = lcm(*(term.denominator for term in HALF_CDF_TERMS))
HALF_CDF_COEFFICIENTS = tuple(int(term*HALF_CDF_DENOMINATOR) for term in HALF_CDF_TERMS)


@lru_cache(maxsize=None)
def half_cdf(value):
    assert 0 < value <= F(5, 4)
    n, d = value.numerator, value.denominator
    n2, d2 = n*n, d*d
    polynomial = 0
    denominator_power = 1
    for coefficient in reversed(HALF_CDF_COEFFICIENTS):
        polynomial = polynomial*n2+coefficient*denominator_power
        denominator_power *= d2
    assert polynomial > 0
    return ADJACENT_UNIT*11*n*polynomial//(8*HALF_CDF_DENOMINATOR*d**15)


def adjacent_certificate(extent=ADJACENT_EXTENT):
    assert 0 < extent <= ADJACENT_EXTENT
    scale = extent/ADJACENT_EXTENT
    closing_margin = scale*F(1,125)
    pi_upper = F(22, 7)
    assert pi_upper*F(11, 8)**2 < 6 and pi_upper*F(11, 4)**2 < 24
    one_side = 0
    width = ADJACENT_WIDTH
    radical = F(2829, 1000)
    assert radical > 2*HI
    a_edges = [scale*F(1,500)]+[j*width for j in range(1,61)]
    k_edges = [scale*F(47,10000)]+[j*width for j in range(1,81)]
    cells = 0
    for a, a_upper in zip(a_edges, a_edges[1:]):
        wa = a_upper-a
        amean = (a+a_upper)/2
        avar = wa**2/12
        neighbor = 0
        for k, k_upper in zip(k_edges, k_edges[1:]):
            wk = k_upper-k
            kmean = (k+k_upper)/2
            contrast = (2*kmean+radical*amean)**2+wk**2/3+radical**2*avar
            child_sum = F(2828, 1000)*k+a-2*closing_margin
            if child_sum <= 0:
                continue
            density = half_density(F(11, 8), F(3, 2), contrast)
            neighbor += (2*wk*density*half_cdf(min(child_sum, F(5, 4)))/ADJACENT_UNIT).__floor__()
            cells += 1
        amplitude = half_density(F(11, 4), F(6), amean**2+avar)
        one_side += (wa*amplitude*neighbor/ADJACENT_UNIT).__floor__()

    probability = F(one_side, ADJACENT_UNIT)**2
    usable = F((probability*10**8).__floor__()-1,10**8)
    assert probability > usable

    drift_lo, drift_hi = F(56, 1000), F(58, 1000)
    sd_lo, sd_hi = F(32, 100), F(35, 100)
    density_lower = F(11, 10)
    assert (1-drift_hi)**2 < F(8, 9) < (1-drift_lo)**2
    assert sd_lo**2 < F(1, 9) < sd_hi**2
    exponent = (drift_hi*extent)**2/(2*sd_lo**2)
    assert exponent < 1
    assert 2*pi_upper*(density_lower*sd_hi)**2 < (1-exponent)**2
    def correlation(distance, radical):
        coefficient = (1+F(distance,3))/2**((distance+1)//2)
        return coefficient*radical if distance%2 else coefficient

    for depth, margin in enumerate((a_edges[0], k_edges[0], closing_margin), 1):
        lower = (correlation(depth,LO)+correlation(depth+1,LO))/(1+correlation(1,HI))
        upper = (correlation(depth,HI)+correlation(depth+1,HI))/(1+correlation(1,LO))
        assert max(abs(1-lower),abs(1-upper))*extent < margin
    return {"scope": "exact adjacent outer probability and scalar rate",
            "pivotal_cells": cells, "exact_single_orientation": str(probability),
            "usable_single_orientation": str(usable),
            "usable_probability": str(2*usable),
            "threshold_abs_max": str(extent),
            "rate_without_normal_density": str(3*drift_lo*density_lower*2*usable)}


def long_baseline(distance, left, offset, normalizer, floor):
    step = F(1, 40)
    cells = range(12)
    sides, middle = [], []
    for i in cells:
        a, b = left+i*step, left+(i+1)*step
        sides.append(strip_box(a, b, offset))
        middle.append((step, (a+b)/2, (a*a+a*b+b*b)/3))

    def exponential(q):
        return exp_index(ceil(1000*q.numerator, q.denominator))

    adjacent = 18*HI
    states = [[(a[0]*b[0]*exponential(18*a[2]+F(39, 2)*b[2]+adjacent*a[1]*b[1])).__floor__()
               for b in middle] for a in sides]
    coordinates = (sides,)+(middle,)*(distance-3)+(sides,)
    for index in range(2, distance-1):
        A, B, C = coordinates[index-2:index+1]
        diagonal = 18 if index == distance-2 else F(39, 2)
        kernel = [[[(c[0]*exponential(diagonal*c[2]+adjacent*b[1]*c[1]+6*a[1]*c[1])).__floor__()
                    for c in C] for b in B] for a in A]
        states = [[sum(states[i][j]*kernel[i][j][k] for i in cells)//U
                   for k in cells] for j in cells]
    exact = normalizer*F(sum(map(sum, states)), U)
    assert exact > floor
    return {"distance": distance, "exact_probability_lower_sum": str(exact),
            "usable_probability": str(floor), "cells": 12**(distance-1),
            "cell_width": str(step), "offset": str(offset)}


def additional_baselines(pi_hi):
    T = EXTENT
    assert T*max(8*HI/11-1,F(1, 33),F(4, 33),1-20*LO/33,F(5, 33)) < F(1, 100)
    assert 0 < 2-14*HI/11 and T*(2-14*LO/11) < F(61, 10000)
    for a, b, den in ((-470,339,526),(-63,46,263),(-155,109,263),
                      (-48,10,263),(-430,243,526),(-135,61,263)):
        assert max(abs(a+b*LO), abs(a+b*HI))*T/den < F(3, 500)
    assert 0 < (1108-713*HI)/526 and T*(1108-713*LO)/526 < F(43, 10000)
    for w_lo, w_hi, variance, drift_lo, drift_hi, sd_lo, sd_hi, density in (
        (F(3, 8), F(3, 8), F(55, 64), F(78, 125), F(627, 1000),
         F(927, 1000), F(929, 1000), F(429, 1000)),
        (5*LO/24, 5*HI/24, F(263, 288), F(141, 200), F(353, 500),
         F(191, 200), F(239, 250), F(417, 1000)),
    ):
        assert drift_lo < 1-w_hi <= 1-w_lo < drift_hi
        assert sd_lo**2 < variance < sd_hi**2
        assert 2*pi_hi < F(627, 250)**2
        assert (1-(drift_hi*T)**2/(2*sd_lo**2))/(sd_hi*F(627, 250)) > density
    norm6, norm7 = F(28), F(399, 10)
    assert norm6**2*(2*pi_hi)**5 < 4*1924560
    assert norm7**2*(2*pi_hi)**6 < 4*24541056
    return [long_baseline(6, F(1, 100), F(61, 10000), norm6, F(272, 10**7)),
            long_baseline(7, F(3, 500), F(43, 10000), norm7, F(59, 10**7))]


def actual_base_deficit(eta, threshold):
    UNIT = U
    SCALE = 1000000000000000000000000
    LEAVES = ('000000', '0000010', '0000011', '00001', '00010', '00011', '0010000', '0010001', '0010010', '0010011', '0010100', '0010101', '0010110', '0010111', '001100', '0011010', '0011011', '001110', '0011110', '0011111', '01000', '010010', '0100110', '0100111', '010100', '010101', '010110', '010111', '0110000', '0110001', '0110010', '0110011', '0110100', '0110101', '0110110', '0110111', '011100', '0111010', '0111011', '011110', '011111', '10000', '10001', '10010', '100110', '100111', '101000', '101001', '101010', '101011', '10110', '101110', '101111', '110000', '110001', '110010', '110011', '110100', '110101', '110110', '110111', '111000', '111001', '111010', '111011', '111100', '111101', '111110', '111111')
    def sqrt_bounds(x):
        assert x >= 0
        n = isqrt((x*SCALE*SCALE).__floor__())
        lo = F(n, SCALE)
        return lo, lo if lo*lo == x else F(n+1, SCALE)

    def box_controls(left, right):
        wl, pl, zl = left
        wh, ph, zh = right
        def roots(value):
            lo, hi = sqrt_bounds(value)
            assert 0 < lo and lo*lo <= value <= hi*hi
            return lo, hi
        rlo, rhi = roots(F(2))
        tl, th = 4*rlo/11, 4*rhi/11
        C_lo = 1+tl/2
        a_lo, _ = roots(F(1, 2)+F(3, 2)*(1-wh)*pl)
        _, a_hi = roots(F(1, 2)+F(3, 2)*(1-wl)*ph)
        s_lo, _ = roots(F(1, 2)+F(3, 2)*wl)
        _, s_hi = roots(F(1, 2)+F(3, 2)*wh)
        v_lo, _ = roots(F(1, 4)+wl*(F(3, 4)+3*zl))
        _, v_hi = roots(F(1, 4)+wh*(F(3, 4)+3*zh))
        diagonal = ((1+tl+tl*v_lo-s_hi)/a_hi,
                    (1+rlo+tl/2+C_lo*v_lo-a_hi)/s_hi,
                    (C_lo*s_lo+tl*a_lo-F(8, 11))/v_hi)
        diagonal = tuple(F((value*10**12).__floor__(), 10**12) for value in diagonal)
        a, b, c = diagonal
        assert min(diagonal) > 0
        second = a*b-1
        determinant = a*b*c-a*F(129, 121)-b*F(32, 121)-c+F(32, 121)
        determinant += (2-a)*(tl if 2-a >= 0 else th)
        assert second > 0 and determinant > F(709, 10**6)
        return determinant


    def concavity_cover():
        leaves = set(LEAVES)
        assert len(leaves) == len(LEAVES) == 69 and max(map(len, leaves)) == 7
        prefixes = {path[:n] for path in leaves for n in range(len(path))}
        assert not leaves.intersection(prefixes)
        records = []
        def visit(path, left, right):
            if path in leaves:
                records.append(box_controls(left, right))
                return
            assert path in prefixes
            axis = max(range(3), key=lambda i: right[i]-left[i])
            midpoint = (left[axis]+right[axis])/2
            low_right, high_left = list(right), list(left)
            low_right[axis] = high_left[axis] = midpoint
            visit(path+"0", left, tuple(low_right))
            visit(path+"1", tuple(high_left), right)
        visit("", (F(0),)*3, (F(1),)*3)
        assert len(records) == len(LEAVES)
        return len(records), min(records)

    PAID_JOINT_COVARIANCE = F(153,50)
    def common_regression_controls():
        t_lo, t_hi = 4*SQRT2_LO/11, 4*SQRT2_HI/11
        spread = (F(94372, 100000)-F(9428, 10000))/(2*(1-F(83843, 100000)))
        a_lo, a_hi = F(2051,4000)-spread, F(103, 200)+spread
        error = F(9,2000)
        assert max(t_hi-a_lo, a_hi-t_lo) < error
        k3 = F(561, 50)
        assert k3 == F(561, 50)
        opposite = k3+F(1, 2)
        assert SQRT2_HI*F(109, 25) < opposite
        correction = error*(F(3, 2)+opposite)
        assert 3+correction < PAID_JOINT_COVARIANCE

    LIMIT = F(9,10000)
    CONTRAST_SD, DEFECT_SD = F(29,100), F(27,1250)
    PAIR_CORR = F(191,1000)
    VAR_LO, VAR_HI = F(279,10000), F(104,3125)
    ENDPOINT_MEAN, MIDDLE_MEAN = F(403,10000), F(443,10000)
    DEFECT_SHIFT, SIBLING_MEAN, SIBLING_SHIFT = F(443,100000), F(67,500), F(109,50000)
    DRIVER_LEFT = F(1, 25)


    def source_controls(eta, threshold, pi_hi):
        assert 0 < eta <= LIMIT and 0 <= threshold <= F(21, 1000)
        lam_lo, lam_hi = 2*SQRT2_LO, 2*SQRT2_HI
        rho_lo, rho_hi = F(8332, 10000), F(83843, 100000)
        edge_lo, edge_hi = F(9428, 10000), F(94372, 100000)
        mean_lo, mean_hi = F(2051,2000), F(103,100)
        assert edge_lo < lam_lo/3 and lam_hi/3+eta < edge_hi
        assert rho_lo < F(5, 6) and F(5, 6)+2*lam_hi*eta < rho_hi
        assert mean_lo < 2*edge_lo/(1+rho_hi) and 2*edge_hi/(1+rho_lo) < mean_hi
        assert VAR_LO < 1-(2*edge_hi)**2/(2*(1+rho_lo))-(edge_hi-edge_lo)**2/(2*(1-rho_hi))
        assert 1-(2*edge_lo)**2/(2*(1+rho_hi)) < VAR_HI
        assert mean_lo/2-(edge_hi-edge_lo)/(2*(1-rho_hi)) > F(1, 2)
        assert mean_hi/2+(edge_hi-edge_lo)/(2*(1-rho_hi)) < F(52, 100)
        assert all(w+d*d*r*eta < 1 for d, w, r in (
            (1, lam_hi/3, F(1)), (2, F(5, 6), SQRT2_HI),
            (3, SQRT2_HI/2, F(2)), (4, F(7, 12), 2*SQRT2_HI),
        ))

        root10_hi = F(79057, 25000)
        assert root10_hi**2 > 10
        sharp_two_step = F(367, 200)
        outside_cov = F(177, 40)*eta
        s_max = 3*eta
        assert max(F(104, 100)-SQRT2_LO/2, SQRT2_HI-1)*s_max/VAR_LO < ENDPOINT_MEAN
        assert max(mean_hi*SQRT2_HI-1, 1-mean_lo*SQRT2_LO/2)*s_max/VAR_LO < MIDDLE_MEAN
        assert 3*s_max/(1+rho_lo) < DEFECT_SHIFT
        assert (sharp_two_step*eta+F(52, 100)*outside_cov)/VAR_LO < SIBLING_MEAN
        assert outside_cov/(1+rho_lo) < SIBLING_SHIFT
        assert 0 < mean_lo+2-lam_hi and mean_hi+2-lam_lo < F(202, 1000)
        assert lam_lo*mean_lo > F(29, 10) and lam_hi*mean_hi < F(31, 10)
        assert DEFECT_SD**2 > (3-lam_lo)*s_max and CONTRAST_SD**2 > F(1, 12)
        assert 4*DEFECT_SD < F(1, 10) and F(2, 2**11) < F(1, 1000)
        assert (mean_hi-1)*threshold < DRIVER_LEFT
        assert (lam_lo-MIDDLE_MEAN)*DRIVER_LEFT-(F(1, 10)+DEFECT_SHIFT)*threshold > F(1, 10)

        raw_cov = SQRT2_HI*s_max/4
        pair_projection = (F(3, 2)/(1+rho_lo)+F(1, 2)/(1-rho_hi))*s_max*outside_cov
        driver_projection = ENDPOINT_MEAN*SIBLING_MEAN*VAR_HI
        assert raw_cov+pair_projection+driver_projection < PAIR_CORR*CONTRAST_SD*DEFECT_SD

        phi = F(797, 2000)
        assert F(627, 250)**2 > 2*pi_hi and (1-threshold**2/2)/F(627, 250) > phi
        assert 1-rho_hi > F(4039,25000)
        assert F(54, 100)**2 < 1-rho_hi**2 and 1-rho_lo**2 < F(56, 100)**2
        assert F(711, 1000)**2*2*pi_hi*F(56, 100)**2 < 1
        exponent = ((1-rho_lo)*threshold)**2/(2*F(54, 100)**2)
        assert F(711, 1000)*(1-exponent) > F(71, 100)
        return phi


    def mean_square(left, right):
        return (left**2+left*right+right**2)/3


    @lru_cache(maxsize=1)
    def defect_cells():
        step = F(1, 40)
        normalizer = F(3989, 5000)
        pi_hi = scalar_controls()
        assert normalizer**2*pi_hi < 2
        return tuple((right, step*normalizer*F(exp_lower(
            mean_square(right-step, right)/(2*(1-PAIR_CORR**2))), UNIT))
            for right in (step*i for i in range(1, 141)))


    def pair_probability(margin):
        total = F(0)
        for defect, mass in defect_cells():
            radius = margin-DEFECT_SD*defect/2
            if radius <= 0:
                break
            strip = central_strip((1000*radius/CONTRAST_SD).__floor__())
            total += mass*F(strip, UNIT)
        return total


    def directed_sqrt(value):
        return F(isqrt((value*UNIT**2).__floor__()), UNIT)


    def certificate(eta, threshold):
        pi_hi = scalar_controls()
        phi = source_controls(eta, threshold, pi_hi)
        densities = [(variance, directed_sqrt(2/(pi_hi*variance))) for variance in (VAR_LO, VAR_HI)]
        slope = F(1,2)-F(99,1000)
        shift = F(101, 1000)+DEFECT_SHIFT/2+SIBLING_SHIFT
        step = F(1, 400)
        total = F(0)
        for i in range(224):
            left, right = DRIVER_LEFT+i*step, DRIVER_LEFT+(i+1)*step
            probability = pair_probability(slope*left-shift*threshold)
            square = mean_square(left, right)
            amplitude = min(normalizer*F(exp_lower(square/(2*variance)), UNIT)
                            for variance, normalizer in densities)
            total += step*amplitude*F(999, 1000)*probability**2
        floor = F((total*10**6).__ceil__()-1, 10**6)
        assert total > floor > 0
        deficit = phi*F(4039,25000)*F(71, 100)*floor/2
        return {"scope": "actual-base D2 deficit from two symmetric convex endpoint events",
                "marking_cutoff": str(eta), "threshold_abs_max": str(threshold),
                "absolute_pair_correlation_bound": str(PAIR_CORR),
                "exact_probability_lower_sum": str(total), "usable_probability": str(floor),
                "deficit_coefficient": str(deficit), "band_weight": str(deficit*threshold**2)}


    assert eta == F(9,10000)
    common_regression_controls()
    assert F(99,1000)*VAR_LO > PAID_JOINT_COVARIANCE*eta
    control = concavity_cover()
    result = certificate(eta, threshold)
    result["hessian_box_count"] = control[0]
    result["minimum_determinant_lower"] = str(control[1])
    return result


def outer_time_moment(threshold, used, ratio):
    assert 0 < used <= threshold < 1 and 0 < ratio < 1
    moment = used**2/2-(F(1, 8)+ratio/24)*used**4+ratio*used**6/72
    moment += (threshold-used-(threshold**3-used**3)/6)*(used-ratio*used**3/6)
    assert moment > 0
    return moment


def short_adjacent_credit(full, short, adjacent_pair, pi_hi, threshold):
    assert adjacent_pair["distance"] == 1
    w = F(adjacent_pair["cover_correlation_upper"])
    ratio = (1-w)/(1+w)
    root = F(isqrt((ratio*10**24).__floor__()),10**12)
    assert 0 < root and root**2 <= ratio
    used = F(short["threshold_abs_max"])
    assert 0 < used <= F(full["threshold_abs_max"]) and used <= threshold
    gain = F(short["usable_probability"])-F(full["usable_probability"])
    assert gain > 0
    value = 3*gain*root*outer_time_moment(threshold, used, ratio)/(2*pi_hi)
    return F((value*U).__floor__(),U)


def integrated_pair_credits(baselines, arm, adjacent, pi_hi, threshold):
    assert 0 < threshold < 1
    used = min(EXTENT, threshold)
    probabilities = [F(adjacent["usable_probability"])]
    probabilities.extend(F(row["usable_probability"]) for row in baselines)
    b2, bb2, b3, bb3, c2, c3 = map(F, arm["short_floors"])
    probabilities[1] += b2+bb2+c2
    probabilities[2] += b3+bb3+c3
    for distance, floor in enumerate(arm["long_floors"], 4):
        probabilities[distance-1] += F(floor)

    credits = []
    for distance, probability in enumerate(probabilities, 1):
        critical_upper = (1+F(distance, 3))/2**((distance+1)//2)
        if distance % 2:
            critical_upper *= HI
        w = critical_upper+F(1, 10**9)
        assert 0 < w < 1 and 0 < EXTENT < 1
        ratio = (1-w)/(1+w)
        root = F(isqrt((ratio*10**24).__floor__()), 10**12)
        assert 0 < root <= 1 and root**2 <= ratio
        # Integrate the product of the two positive lower Taylor factors
        # over 0 <= s <= t <= EXTENT.
        moment = outer_time_moment(threshold, used, ratio)
        value = 3*2**(distance-1)*probability*root*moment/(2*pi_hi)
        credit = F((value*U).__floor__(), U)
        assert credit > 0
        credits.append({"distance": distance, "negative_mixed_response_floor": str(probability),
                        "cover_correlation_upper": str(w), "credit": str(credit)})
    return credits


def certificate():
    pi_hi, norms = baseline_controls()
    s_max = F(3, 100)
    variance = F(95, 144)
    assert 1 - F(7, 12) > F(52, 125)
    assert 2*pi_hi*variance*F(49, 100)**2 < (1-s_max**2/(2*variance))**2
    assert 0 < 14*LO/19-1 < 14*HI/19-1 < F(1, 15)
    assert s_max*max(F(2, 19), 1-12*LO/19) < F(1, 100)
    assert s_max*(2-24*LO/19) < F(13, 2000)
    variance = F(7, 9)
    assert 1-HI/3 > F(66, 125)
    assert 2*pi_hi*variance*F(9, 20)**2 < (1-s_max**2/(2*variance))**2
    assert s_max*max((17*HI-23)/28,(4*HI-5)/14,(3-LO)/14,(19-11*LO)/28) < F(1, 100)
    assert s_max*(61-39*LO)/28 < F(13,2000)
    offsets = (2-14*LO/11, F(7,3)-3*LO/2, 2-24*LO/19, (61-39*LO)/28)
    assert all(0 < EXTENT*shift < F(bound,10000)
               for shift,bound in zip(offsets,(46,49,49,48),strict=True))
    baselines = weighted_distances(norms)+additional_baselines(pi_hi)
    adjacent = adjacent_certificate()
    assert EXTENT <= F(adjacent["threshold_abs_max"])
    arm = associated_arm_certificate(last_distance=6)
    assert [row["distance"] for row in baselines] == [2, 3, 4, 5, 6, 7]
    assert EXTENT <= F(arm["threshold_abs_max"])
    assert F(arm["common_normal_density_floor"]) == PHI
    credits = integrated_pair_credits(baselines, arm, adjacent, pi_hi, F(49,2000))
    diamond = actual_base_deficit(F(9,10000), F(101,5000))
    shorter = adjacent_certificate(F(1,100))
    extra = short_adjacent_credit(adjacent, shorter, credits[0], pi_hi, F(49,2000))
    saving = sum(F(row["credit"]) for row in credits)+extra
    assert saving > REQUIRED_CREDIT
    return {
        'threshold_abs_max_of_central_band': str(F(49,2000)),
        'actual_base_deficit': diamond,
        'outer_graph_credit': str(saving),
        'slack_above_uniform_saving': str(saving-REQUIRED_CREDIT),
    }


outer = SimpleNamespace(certificate=certificate)


SCORE_UNIT = 1000000000
PROBABILITY_UNIT = 1000000000000
MARGIN_STEP_END = 500000

def midpoint_scales(eta):
    root_lo, root_hi = F(141421356,10**8), F(141421357,10**8)
    assert root_lo**2 < 2 < root_hi**2
    assert 0 < eta <= F(21,20000)
    variance = 3*(3-2*root_lo)*eta
    sigma = F(isqrt((variance*10**12).__floor__())+1,10**6)
    half = F(289,1000)
    assert sigma**2 > variance and half**2 > F(1,12)
    midpoint_sd = (1+root_hi/2)*sigma
    correlation = (5*root_hi/8)*(3*eta)/(midpoint_sd*half)
    assert 0 < correlation < 1
    return midpoint_sd,half,1-correlation**2


def positive_defect_bins(sigma, step):
    normal = F(79788,10**5)
    pi_upper = F(355,113)
    atan5 = sum((F((-1)**j,(2*j+1)*5**(2*j+1)) for j in range(9)),F(0))
    assert 16*atan5-4*(F(1,239)-F(1,3*239**3)) < pi_upper
    assert normal**2*pi_upper < 2
    scale = PROBABILITY_UNIT
    bins = []
    for lo in range(0,120000000,step):
        hi = lo+step
        exponent = F(lo*lo+lo*hi+hi*hi,6*SCORE_UNIT**2)/sigma**2
        assert exponent < 4096
        value = (scale*(1-exponent/4096)).__floor__()
        for _ in range(12):
            value = value*value//scale
        mass = (F(step,SCORE_UNIT)*normal*value/sigma).__floor__()
        bins.append((hi,mass))
    return tuple(bins)


def midpoint_strip_table(eta, strip):
    sigma,half,factor = midpoint_scales(eta)
    bins = positive_defect_bins(sigma,1000000)
    counts = []
    for margin in range(0,300000001,MARGIN_STEP_END):
        total = 0
        for hi,mass in bins:
            available = margin-hi
            if available <= 0:
                break
            index = available*10000*half.denominator//(SCORE_UNIT*half.numerator)
            total += mass*strip[min(30000,index)]
        counts.append(total*factor.numerator//(10**6*factor.denominator))
    assert counts == sorted(counts) and 0 < counts[-1] <= 10**15
    return tuple(counts)


ROOT_TICKS = (10000000, 16000000, 22000000, 28000000, 34000000, 40000000, 46000000)
MID_STEP = 1000000
MID_MAX_INDEX = 90
CORE_STEP = 1000000
CORE_SCALE = 1000000000000
SCORE_TICKS = 1000000000
CORR_PANELS = 32

def adaptive_core_table(eta, strip):
    lo2,hi2 = F(141421356,10**8), F(141421357,10**8)
    assert lo2**2 < 2 < hi2**2
    z,c = 3*eta,3-2*lo2
    assert 0 < 2*z < c
    vx,vy = c*z-z*z,c*z-z*z/2
    sx,sy = (F(isqrt((v*10**12).__floor__())+1,10**6) for v in (vx,vy))
    assert sx**2 > vx and sy**2 > vy
    bins = positive_defect_bins(sx,CORE_STEP)[:ROOT_TICKS[-1]//CORE_STEP]
    assert len(bins)*CORE_STEP == ROOT_TICKS[-1]
    table = [[CORE_SCALE]*(MID_MAX_INDEX+1) for _ in ROOT_TICKS]
    edges = [radius//CORE_STEP for radius in ROOT_TICKS]
    assert all(i*CORE_STEP == radius for i,radius in zip(edges,ROOT_TICKS,strict=True))

    for panel in range(CORR_PANELS):
        left = F(2*panel-CORR_PANELS,CORR_PANELS)
        right = F(2*panel+2-CORR_PANELS,CORR_PANELS)
        assert left*right >= 0
        nearest = min(abs(left),abs(right))
        sd_square = sy**2*(1-nearest**2)
        sd_ticks = isqrt((sd_square*SCORE_TICKS**2).__floor__())+1
        multiplier = max(abs(alpha+rho*sy/sx) for alpha in (2*lo2/3,2*hi2/3)
                         for rho in (left,right))
        mean_caps = [((hi*multiplier).__ceil__(),weight) for hi,weight in bins]
        for j in range(MID_MAX_INDEX+1):
            width = 2*j*MID_STEP
            partials = [0]
            total = 0
            for mean,weight in mean_caps:
                if width > mean:
                    inner = min(30000,10000*(width-mean)//sd_ticks)
                    outer = min(30000,10000*(width+mean)//sd_ticks)
                    total += weight*(strip[inner]+strip[outer])
                partials.append(total)
            for row,i in zip(table,edges,strict=True):
                row[j] = min(row[j],partials[i]//(2*10**9))

    assert all(row == sorted(row) for row in table)
    assert all(all(a <= b for a,b in zip(lo,hi,strict=True)) for lo,hi in zip(table,table[1:]))
    return table


def medium_payment_slacks(rho, shared, payment, anchor, positive10, epsilon, t, z02, k):
    """Check the medium and run budgets in units of the Price allowance."""
    sigma = F(2021,1250)
    a0,kappa = sigma-F(3,2),F(59599,10**6)
    c,d = a0/4,(sigma-1)/8
    g = F(2,3)*(shared-c-rho/4)
    B = c+g/2
    nu,nu2 = (F(1,4)+rho)*a0,shared-F(1,8)
    released = anchor-max(F(0),3*payment-2*B)
    assert k == F(33,100) and 0 <= t and 0 <= z02
    assert F(2,5) <= rho <= F(1,2)
    assert F(1,8)+rho/4 <= shared <= F(1,4)+rho/4
    return {
        "type01_pair": g/4-t,
        "type02_one_leaf": g/2-t-z02,
        "type02_two_leaves": kappa/2-2*t-z02,
        "short_high10": min(nu2,nu+F(1,8))-payment-positive10-t,
        "long_high10": 2*B+nu-positive10-3*payment-t,
        "medium_anchor_N2": 2*d+B-c-2*payment-t,
        "medium_anchor_N4": 3*B+k/8-4*payment-t,
        "closed_N2": d-payment-t,
        "closed_N4": B+min(F(1,2),2*shared)*a0-2*payment-t,
        "closed_N6": 2*B+k/8-3*payment-t,
        "long_high_anchor": released-t,
        "single_short": (5-3*sigma)/4+2*anchor-epsilon-payment-t,
        "type20_short_and_high10": kappa+2*anchor-epsilon-payment-positive10-t,
        "short_long": kappa+anchor-payment+released-epsilon-2*t,
        "double_short_ordinary": kappa+2*anchor-2*payment-epsilon-2*t,
    }

