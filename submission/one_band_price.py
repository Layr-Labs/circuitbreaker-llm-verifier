"""Central Price payments and the population inequalities."""
from types import SimpleNamespace
from one_band_outer import d7, central_strip


"""Fixed coarser adjacent panels with exact positive-power integration."""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import factorial,isqrt,lcm
from operator import mul
import json

U=10**18
LO,HI,Q=F(141421,100000),F(70711,50000),F(70711,100000)
C,ADJACENT_CAP=F(8579,50000),F(27,10000)
LANE_BYTES=18
GLOBAL_EXTRA_STEPS=(2,3,4,6,7,8,10,12,14,16,18,20,22,23)
CLEAN_EXTRA_STEPS=(2,4,6,7,8,10,12,14,16,18,20,22,23)


def digest(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()


panel_exp_lower = lru_cache(None)(d7.exp_lower)


panel_strip = central_strip


@lru_cache(None)
def amplitude_density(index):return 27638*panel_exp_lower(6*F(index,1000)**2)//10000


def conditional_factor(s,d):
    return 1-s/2**d if d else(1-s if s<=F(1,2)else 1/(4*s))


def adjacent_deviation(s,d,clean=False):
    if clean:
        variance=3*C*ADJACENT_CAP-9*ADJACENT_CAP*ADJACENT_CAP*s/2**d
        conditioned=C*C*conditional_factor(s,d)
        assert 0<variance<=conditioned
        variance=min(variance,conditioned)
        upper=F(isqrt(variance.numerator*10**18//variance.denominator)+1,10**9)
        assert upper**2>variance
        return upper
    factor=conditional_factor(s,d)
    upper=F(isqrt(factor.numerator*10**18//factor.denominator)+1,10**9)
    assert upper**2>factor
    return C*upper


def impulses(s,depth,clean=False):
    deviations=[adjacent_deviation(s,d,clean and d<=2)for d in range(depth)]
    return [F(0)]+[sum((j+1)*Q**j*deviations[d-1-j]for j in range(d))/2 for d in range(1,depth+1)]


def strips(count,s,sd,spatial,divisor=1):
    numerator=1000**2*(s.denominator-s.numerator)*sd.denominator**2
    denominator=(spatial*divisor)**2*s.numerator*sd.numerator**2
    return [panel_strip(isqrt(i*i*numerator//denominator))for i in range(count)]


def pack(values):
    assert len(values)<=680 and all(0<=v<=2*U for v in values)
    return int.from_bytes(b''.join(v.to_bytes(LANE_BYTES,'little')for v in values),'little')


def convolution(prefix_packed,values,length):
    assert len(values)<=680 and 0<=length<=len(values)
    assert all(0<=v<=U for v in values)
    product=prefix_packed*pack(values)
    raw=product.to_bytes(2*len(values)*LANE_BYTES,'little')
    return [int.from_bytes(raw[i*LANE_BYTES:(i+1)*LANE_BYTES],'little')for i in range(length)]


def controls():
    def atan(x,last):return sum(((-1)**j*x**(2*j+1)/(2*j+1)for j in range(last+1)),F(0))
    pi_hi=16*atan(F(1,5),14)-4*atan(F(1,239),15)
    assert pi_hi<F(31416,10000)
    assert LO**2<2<HI**2 and Q**2>F(1,2) and C==3-2*LO
    assert (3-F(4289,25000))**2>8 and 12*ADJACENT_CAP<F(4289,25000)
    assert pi_hi*F(27638,10000)**2<24 and pi_hi*F(13819,10000)**2<6
    assert pi_hi*F(7978,10000)**2<2 and pi_hi*F(954,1000)<3
    assert 680*2*U**2<1<<(8*LANE_BYTES)


def mean_square(k):
    # Twelve times the mean square on [k-1, k].
    return 12*k*(k-1)+4


def mean_product(k, j):
    # Twelve times the mean product on two independent unit cells.
    return 3*(2*k-1)*(2*j-1)


def adjacent_density_index(x, k):
    # Upper integer index for 1000 times the RMS shifted coordinate.
    a, b = 200000, 282843
    numerator = a*a*mean_square(k)+2*a*b*mean_product(k,x)+b*b*mean_square(x)
    denominator = 12*20000**2
    root = isqrt(numerator//denominator)
    return root+(root*root*denominator<numerator)


GLOBAL_BASELINE_STEPS = (*range(3, 79, 3), 79)
CLEAN_BASELINE_STEPS = (*range(6, 80, 6), 79)


def baseline_coefficients():
    spatial,nx,nk,nl=200,130,165,680
    az=[27638*panel_exp_lower(6*F(mean_square(i+1),12*spatial**2))//10000 for i in range(nx)]
    dz=[0]+[13819*panel_exp_lower(F(3,2)*F(mean_square(q),12*spatial**2))//10000 for q in range(1,nl)]
    @lru_cache(None)
    def fd(index):return 13819*panel_exp_lower(F(3,2)*F(index,1000)**2)//10000
    cells=[]
    for i in range(nx):
        row=[]
        for j in range(nk):
            index=adjacent_density_index(i+1,j+1)
            lower=(282842*j+100000*i)//100000
            coefficient=2*(az[i]*fd(index)//U)//spatial**2
            row.append((lower,coefficient))
        cells.append(row)
    return pack(dz),cells


def baseline(coeff,clean=False):
    dz,cells=coeff;previous=total=F(0);rows=[]
    for step in (CLEAN_BASELINE_STEPS if clean else GLOBAL_BASELINE_STEPS):
        s=F(step,80)**2;stds=impulses(s,3,clean)[1:]
        sx=strips(len(cells),s,stds[0],200);sk=strips(len(cells[0]),s,stds[1],200);sz=strips(680,s,stds[2],200,2)
        conv=[v//(U*200)for v in convolution(dz,sz,680)]
        side=0
        for i,row in enumerate(cells):
            value=sum(coefficient*sk[j]*conv[lower] for j,(lower,coefficient)in enumerate(row))//(U*U)
            side+=sx[i]*value
        side//=U
        credit=F(side,U)**2*((1-previous)**4-(1-s)**4)/4
        total+=credit
        rows.append({'s':str(s),'side_1e18':side,'stds':list(map(str,stds)),'time_credit_without_density':str(credit)})
        previous=s
    return F(954,1000)*total,rows


@lru_cache(maxsize=1)
def arm_transition():
    assert SQRT2_HI == HI and SQRT2_LO == LO
    size, spatial = 48, 80
    first = [[amplitude_density((50000*(k+1)+70711*(j+1)+3999)//4000)//spatial
              for k in range(size)] for j in range(size)]
    continuation = [[tuple(amplitude_density((50000*(k+1)+70711*(j+1)+25000*(i+1)+3999)//4000)//spatial
                           for i in range(size)) for k in range(size)] for j in range(size)]
    lower = [[int(2*LO*k+j) for k in range(size)] for j in range(size)]
    return first, continuation, lower


def extra_coefficients():
    size, spatial = 48, 80
    first, continuation, lower = arm_transition()
    root = [amplitude_density((25*(j+1)+1)//2)//spatial for j in range(size)]
    maximum = max(map(max, lower))
    dz = [0]+[13819*panel_exp_lower(F(3,2)*F(q,spatial)**2)//(10000*spatial)for q in range(1,maximum+1)]
    return root, first, continuation, lower, pack(dz), maximum


@lru_cache(maxsize=1)
def continuation_hankel():
    # For fixed j, continuation[j][k][i] depends on i+2*k.
    result = []
    for table in arm_transition()[1]:
        sequence = list(table[0])
        for row in table[1:]:
            sequence.extend(row[-2:])
        assert len(sequence) == 142 and all(0 <= w <= U for w in sequence)
        assert all(table[k][i] == sequence[i+2*k] for k in range(48) for i in range(48))
        result.append(pack(sequence))
    assert 48*2*U**2 < 1 << (8*LANE_BYTES)
    return tuple(result)


def hankel_dots(column, kernel):
    assert 0 <= min(column) <= max(column) <= 2*U
    product = pack(tuple(reversed(column)))*kernel
    raw = product.to_bytes(190*LANE_BYTES, 'little')
    return [int.from_bytes(raw[(47+2*k)*LANE_BYTES:(48+2*k)*LANE_BYTES], 'little')
            for k in range(48)]


def one_time(s,coeff,clean=False):
    root,first,continuation,lower,dz,maximum=coeff
    assert continuation is arm_transition()[1]
    hankel = continuation_hankel()
    ell=impulses(s,9,clean)
    selected={d:strips(48,s,ell[d],80)for d in range(1,9)}
    initial=[a*b//U for a,b in zip(root,selected[1])]
    state=[[initial[j]*first[j][k]//U*selected[2][k]//U for k in range(48)]for j in range(48)]
    probabilities={}
    for n in range(2,9):
        if n%2==0:
            last_strip=strips(maximum+1,s,ell[n+1],80,2)
            terminal=[v//U for v in convolution(dz,last_strip,maximum+1)]
            probabilities[n]=sum(state[j][k]*terminal[lower[j][k]]for j in range(48)for k in range(48))//U
        if n<8:
            state = [
                [value//U*strip//U
                 for value, strip in zip(hankel_dots(column, hankel[j]), selected[n+1], strict=True)]
                for j, column in enumerate(zip(*state, strict=True))
            ]
    b2=(1-s)*F(probabilities[2],U)
    added=sum(((1-s)**(n//2)*F(probabilities[n],U)for n in(4,6,8)),F(0))
    return probabilities,(1-s)*(2*b2*added+added*added),ell


def power_integral(left,right,power):
    return ((1-left)**(power+1)-(1-right)**(power+1))/(power+1)


def extra(coeff,clean=False):
    previous=total=F(0);rows=[]
    steps=CLEAN_EXTRA_STEPS if clean else GLOBAL_EXTRA_STEPS
    for step in steps:
        s=F(step,24)**2;probabilities,_,ell=one_time(s,coeff,clean)
        p={n:F(v,U)for n,v in probabilities.items()}
        integral=sum((2*p[2]*p[n]*power_integral(previous,s,2+n//2)for n in(4,6,8)),F(0))
        integral+=sum((p[n]*p[k]*power_integral(previous,s,1+(n+k)//2)for n in(4,6,8)for k in(4,6,8)),F(0))
        credit=F(954,1000)*integral;total+=credit
        rows.append({'s':str(s),'left_time':str(previous),'p_n_scaled':{str(n):v for n,v in probabilities.items()},
                     'impulse_deviations':list(map(str,ell[1:])),'integrated_time_credit':str(credit)})
        previous=s
    return total,rows


def certificate(progress=None):
    controls();bc=baseline_coefficients();ec=extra_coefficients()
    global_base,global_base_rows=baseline(bc)
    clean_base,clean_base_rows=baseline(bc,True)
    if progress:progress('cell-mean adjacent baselines exceed retained floors')
    global_extra,global_extra_rows=extra(ec)
    clean_extra,clean_extra_rows=extra(ec,True)
    if progress:progress('14-panel global and13-panel clean EXTRA reconstructed')
    assert global_extra==F('11599092875786601015788255843706960198011833614086425982882761/81423382595151116076318720000000000000000000000000000000000000000')
    assert clean_extra==F('10077934980568020478820506524069345124850497626881918982704299853/15633289458269014286653194240000000000000000000000000000000000000000')
    values=(global_base,global_extra,clean_base,clean_extra)
    floors=(F(39059609,10**10),F(1424539,10**10),F(88992826,10**10),F(6374849,10**10))
    assert all(v>f>0 for v,f in zip(values,floors))
    kbase=floors[2]-floors[0];kextra=floors[3]-floors[1]
    transfer=F(149647,10**9)-floors[1]
    assert transfer==F(71931,10**10) and kbase+kextra==F(54883527,10**10)
    return {'scope':'27/14 cell-mean adjacent baseline panels; identical retained floors and Price split',
            'external_numerical_inputs':[],'clean_cap':str(ADJACENT_CAP),'clean_window_tree_edges':29,
            'clean_residual_depths':[0,1,2],'larger_extra_global_residual_depths':[3,4,5,6,7,8],
            'components':[{'name':name,'exact_lower_sum':str(v),'usable_floor':str(f),'rows_sha256':digest(rows)}for name,v,f,rows in zip(('global_baseline','global_extra','clean_baseline','clean_extra'),values,floors,(global_base_rows,global_extra_rows,clean_base_rows,clean_extra_rows))],
            'global_K0_payment':str(floors[0]+floors[1]),'local_baseline_increment':str(kbase),'local_extra_increment':str(kextra),
            'local_adjacent_raw_cap':str(kbase+kextra),'global_to_clean_price_transfer':str(transfer),
            'global_baseline_time_steps':list(GLOBAL_BASELINE_STEPS),
            'clean_baseline_time_steps':list(CLEAN_BASELINE_STEPS),
            'global_EXTRA_time_steps':list(GLOBAL_EXTRA_STEPS),'clean_EXTRA_time_steps':list(CLEAN_EXTRA_STEPS),
            'replacement_scope':'global EXTRA floor decreases by transfer; the clean replacement increases by the same transfer'}


adjacent = SimpleNamespace(certificate=certificate, pack=pack, convolution=convolution)


def divided_convolution(prefix, values):
    # The omitted zero term is explicit in the packed prefix. Both ranges
    # use the same 10**18 scale as the existing adjacent convolution.
    assert prefix[0] == 0 and len(values) <= len(prefix)
    packed = adjacent.pack(prefix[:len(values)])
    return [value//10**18 for value in adjacent.convolution(packed, values, len(values))]


"""Exact standalone joint-endpoint baseline Price floor at distance two."""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import factorial,isqrt,lcm
import json

U,CELLS,SPATIAL,TIME=10**18,160,200,80
C_UPPER=F(8579,50000)
MULTIPLIER=F(141421,25000)
ROWS_SHA="1c07445b62a6b307979dc67dad458f8ef7aa7ededa2408f3cb7ad4c711211be4"


@lru_cache(None)


@lru_cache(None)


def baseline_controls():
    def atan(x,last):return sum(((-1)**j*x**(2*j+1)/(2*j+1)for j in range(last+1)),F(0))
    pi_lo=16*atan(F(1,5),15)-4*atan(F(1,239),14)
    pi_hi=16*atan(F(1,5),14)-4*atan(F(1,239),15)
    assert F(31415,10000)<pi_lo<pi_hi<F(31416,10000)
    assert F(45834,10000)**2*pi_hi<66
    assert F(13819,10000)**2*pi_hi<6
    assert F(7978,10000)**2*pi_hi<2
    assert pi_hi**2*11*F(2878,10000)**2<9
    assert 2*F(141421,100000)+C_UPPER==3
    assert F(141421,100000)**2<2 and MULTIPLIER**2<32


def variance_factor(s,distance):
    if distance:return 1-s/2**distance
    return 1-s if s<=F(1,2)else 1/(4*s)


def upper_root_integer(factor):
    value=isqrt(factor.numerator*10**18//factor.denominator)+1
    assert value*value*factor.denominator>factor.numerator*10**18
    return value


def strip_values(coef_num,coef_den,s):
    numerator=(1000*coef_num)**2*(s.denominator-s.numerator)
    denominator=(SPATIAL*coef_den)**2*s.numerator
    return [panel_strip(isqrt(i*i*numerator//denominator))for i in range(CELLS)]


def certificate(progress=None):
    baseline_controls()
    xd=[45834*panel_exp_lower(F(33*(i+1)**2,2*SPATIAL**2))//10000 for i in range(CELLS)]
    zd=[13819*panel_exp_lower(F(3*(j+1)**2,2*SPATIAL**2))//10000 for j in range(CELLS)]
    previous=old_total=new_total=F(0);rows=[]
    for step in range(1,TIME):
        s=F(step,TIME)**2;sn,sd=s.numerator,s.denominator
        h0=variance_factor(s,0)
        pm=upper_root_integer(variance_factor(s,1));po=upper_root_integer(h0);ph=upper_root_integer(min(h0,F(1,2)))
        middle=strip_values(141421*10**9,8579*pm,s)
        old_strip=strip_values(50000*10**9,8579*po,s)
        high_strip=strip_values(50000*10**9,8579*ph,s)
        gain=isqrt(4*sd*U**2//(4*sd-sn));assert gain>=U
        old_row=new_row=low_choices=high_choices=0
        for i,density in enumerate(xd):
            old_conv=sum(zd[j]*old_strip[i-j-1]for j in range(i))//(U*SPATIAL)
            high_conv=sum(zd[j]*high_strip[i-j-1]for j in range(i))//(U*SPATIAL)
            low_conv=old_conv*gain//U;joint_conv=min(low_conv,high_conv)
            assert joint_conv>=old_conv
            low_choices+=low_conv<high_conv;high_choices+=high_conv<low_conv
            weight=density*middle[i]//U
            old_row+=(weight*(old_conv*old_conv//U)//U)//SPATIAL
            new_row+=(weight*(joint_conv*joint_conv//U)//U)//SPATIAL
        if 2*sn>=sd:assert old_row==new_row
        root=isqrt(1000*(sd-sn)*U**2//(1000*sd+313*sn))
        determinant=F(sd-sn,sd)*F(root,U)
        factor=F(2878,10000)*(s-previous)*determinant
        old_total+=factor*F(old_row,U);new_total+=factor*F(new_row,U)
        rows.append({'s':str(s),'conditioning_row_1e18':old_row,'joint_row_1e18':new_row,
                     'low_case_strict_min_cells':low_choices,'high_case_strict_min_cells':high_choices})
        previous=s
    rows_sha=sha256(json.dumps(rows,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    assert rows_sha==ROWS_SHA
    old_floor=F((old_total*10**10).__floor__()-1,10**10)
    new_floor=F((new_total*10**10).__floor__()-1,10**10)
    assert old_floor==F(29238757,10**10) and new_floor==F(30239457,10**10)
    assert old_total==F(18712805569436220172821816543260857063389929,6400000000000000000000000000000000000000000000)
    assert new_total==F(619304114257511778362972946836495446012165997,204800000000000000000000000000000000000000000000)
    assert new_total>new_floor>old_floor
    if progress:progress('all79 D2 joint baseline rows match the original exact sums')
    return {'scope':'global D2 no-arm baseline Price only; other K0 terms and capped/local groups separate',
            'external_numerical_inputs':[],'spatial_denominator':SPATIAL,'cells_per_axis':CELLS,'time_denominator_root':TIME,
            'rows_sha256':rows_sha,'conditioning_exact_lower_sum':str(old_total),'joint_exact_lower_sum':str(new_total),
            'conditioning_usable_floor':str(old_floor),'joint_usable_floor':str(new_floor),
            'spectral_multiplier_lower':str(MULTIPLIER),'full_global_price_payment':str(MULTIPLIER*new_floor),
            'replacement_increment':str(MULTIPLIER*(new_floor-old_floor)),
            'replacement_scope':'new floor replaces its conditioning-only baseline; no arm class is included'}


baseline2 = SimpleNamespace(certificate=certificate, CELLS=CELLS, SPATIAL=SPATIAL, TIME=TIME, strip_values=strip_values)


"""Original D3-D4 and coarser positive-power D5-D6 Price payments."""

from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations
import json
from math import isqrt

C_HI = d7.C_HI
Q_HI = d7.Q_HI
SQRT2_HI = d7.SQRT2_HI
SQRT2_LO = d7.SQRT2_LO
U = d7.U
exp_index = d7.exp_index
exp_lower = d7.exp_lower
margin = d7.margin
scalar_controls = d7.scalar_controls
strip = d7.strip


TIMES = tuple(F(i, 24)**2 for i in range(1, 24))
GLOBAL_FLOORS = (F(1734471, 5000000000), F(561159, 10**10), F(99777, 10**10), F(2157, 1250000000))
LOCAL_FLOORS = (F(5711963, 10**10), F(189007, 2000000000), F(31747, 2000000000), F(27707, 10**10))
GLOBAL_MEANS = (F(1467, 1000), F(1281, 1000), F(1213, 1000), F(144, 125))
LOCAL_MEANS = (F(1111, 1000), F(227, 200), F(1193, 1000), F(1117, 1000))


def deviation(s, distance, cap=None, high=False):
    factor = 1-s/2**distance if distance else (1-s if s <= F(1, 2) else 1/(4*s))
    if cap is not None:
        variance = 3*C_HI*cap-9*cap**2*s/2**distance
        assert 0 < variance <= C_HI**2*factor
        return sqrt_upper(variance)
    if high and not distance:
        factor = min(factor, F(5, 9))
    return C_HI*sqrt_upper(factor)


def sqrt_upper(q):
    scale = 10**9
    bound = F(isqrt(q.numerator*scale**2//q.denominator)+1, scale)
    assert bound**2 > q
    return bound


def sqrt_lower(q):
    return F(isqrt(q.numerator*U**2//q.denominator), U)


@lru_cache(maxsize=1)
def arm_geometry():
    width, arms, central = 80, 48, 64

    def fd(x):
        return 13819*exp_lower(F(3, 2)*x*x)//10000

    roots = [[2*fd(F(x+2*a+3, width))//width for a in range(arms)]
             for x in range(central)]
    first, _, lower = arm_transition()
    largest = max(central-1, max(map(max, lower)))
    dz = [0]+[fd(F(z, width))//width for z in range(1, largest+1)]
    return roots, first, lower, dz


def arm_arrays(incoming, times=TIMES):
    for s in times:
        yield arm_panel(s, tuple(incoming(s, k) for k in range(3)))


@lru_cache(maxsize=None)
def arm_panel(s, deviations):
    roots, first, lower, dz = arm_geometry()
    arms, central = 48, 64
    largest = len(dz)-1
    amplitude = [sum((i+1)*Q_HI**i*deviations[d-1-i] for i in range(d))/2
                 for d in (1, 2, 3)]
    root_strip = [strip(i, s, amplitude[0]) for i in range(arms)]
    child_strip = [strip(i, s, amplitude[1]) for i in range(arms)]
    terminal_strip = [strip(i, s, amplitude[2], 2) for i in range(largest+1)]
    baseline_strip = [strip(i, s, deviations[0]) for i in range(central)]
    terminal = divided_convolution(dz, terminal_strip)
    baseline = divided_convolution(dz, baseline_strip)
    inner = [sum((first[a][b]*child_strip[b]//U)*terminal[lower[a][b]]
                 for b in range(arms))//U for a in range(arms)]
    size_two = [sum((row[a]*root_strip[a]//U)*inner[a] for a in range(arms))//U
                for row in roots]
    assert 0 <= min(baseline+size_two) <= max(baseline+size_two) <= U
    return baseline, size_two


def sides(s, baseline, size_two):
    return [s.denominator*b+(s.denominator-s.numerator)*a
            for b, a in zip(baseline, size_two, strict=True)]


def strip_table(s, sd):
    return [strip(i, s, sd) for i in range(margin(63, 63, 63)+1)]


def normalizer_check():
    pi = scalar_controls()
    assert F(81027, 10000)*pi < 18*SQRT2_LO
    assert F(225, 1000)*2*pi < SQRT2_LO
    numerators = (F(12862, 1000), F(19301, 1000), F(280379, 10000))
    for size, expected, norm in zip((3, 4, 5), (10260, 145152, 1924560), numerators, strict=True):
        total = 0
        for p in permutations(range(size)):
            rational, radical = 1, 0
            for i, j in enumerate(p):
                gap = abs(i-j)
                if not gap:
                    diagonal = 36 if i in (0, size-1) else 39
                    rational *= diagonal
                    radical *= diagonal
                elif gap == 1:
                    rational, radical = -36*radical, -18*rational
                elif gap == 2:
                    rational *= 6
                    radical *= 6
                else:
                    rational = radical = 0
            assert radical == 0
            inversions = sum(p[i] > p[j] for i in range(size) for j in range(i+1, size))
            total += (-1)**inversions*rational
        assert total == expected and norm**2*(2*pi)**size < 4*total
    assert F(195, 1000)**2*pi**2*95 < 36
    assert F(1804, 10000)**2*4*pi**2*7 < 9
    assert F(1716, 10000)**2*pi**2*55 < 16


def mean_check():
    root10_lo, root10_hi = F(3162277, 10**6), F(79057, 25000)
    assert root10_lo**2 < 10 < root10_hi**2
    lower = (root10_lo+3)/(2*SQRT2_HI)
    upper = (root10_hi+3)/(2*SQRT2_LO)
    theta, w5 = 38-12*root10_lo, SQRT2_HI/3
    assert 0 < theta < 1
    global_bounds = (SQRT2_HI+theta, 1+F(2, 3)*theta+F(14, 57),
                     1+theta+w5/(2*(1+w5)), 1+F(4, 5)*theta+F(6, 55))
    assert all(a < b for a,b in zip(global_bounds, GLOBAL_MEANS, strict=True))
    cap6, cap7 = F(3, 500), F(7, 1000)
    assert 12*cap7 < 3-2*SQRT2_HI
    assert 6*lower > F(468, 95)*SQRT2_HI and 12*lower > F(1038, 95)*SQRT2_HI
    assert 6*lower > F(24, 7)*SQRT2_HI and 12*lower > F(66, 7)*SQRT2_HI
    assert 6*upper < 12*SQRT2_LO and 12*lower > 18*SQRT2_HI
    local_bounds = (
        1+24*cap6*(upper-SQRT2_LO),
        1+cap6*(24*upper-F(2008, 95)*SQRT2_LO),
        1+cap7*(24*upper-F(123, 7)*SQRT2_LO),
        1+F(2, 5)*(4*cap7*(12*upper-18*SQRT2_LO)+F(3, 11)),
    )
    assert all(a < b for a,b in zip(local_bounds, LOCAL_MEANS, strict=True))


@lru_cache(maxsize=1)
def three_geometry():
    width = 80
    cells = range(40)
    return [(i, j, margin(j, i), margin(i, j),
                  81027*exp_lower(18*(F(mean_square(i+1)+mean_square(j+1), 12*width**2)
                                +SQRT2_HI*F(mean_product(i+1,j+1), 12*width**2)))//(10000*width**2))
            for i in cells for j in cells]


def three_paths(global_arms, high_arms, local_arms):
    geometry = three_geometry()
    global_total = local_total = previous = F(0)
    for s, old, high, local in zip(TIMES, global_arms, high_arms, local_arms, strict=True):
        st = strip_table(s, deviation(s, 1))
        lt = strip_table(s, deviation(s, 1, F(3, 500)))
        b, h, l = (sides(s, *row) for row in (old, high, local))
        ll = lh = hl = hh = clean = 0
        for i, j, a, c, weight in geometry:
            w = (weight*st[a]//U)*st[c]//U
            ll += w*b[i]*b[j]
            lh += w*b[i]*h[j]
            hl += w*h[i]*b[j]
            hh += w*h[i]*h[j]
            wl = (weight*lt[a]//U)*lt[c]//U
            clean += wl*l[i]*l[j]
        gain = sqrt_lower(1/(1-s/3))
        pair = min(gain*gain*ll, gain*lh, gain*hl, F(hh))
        common = F(225, 1000)*(s-previous)*(1-s)/s.denominator**2/U**3
        r = s/(1-s)
        global_total += common*pair*sqrt_lower(1/(1+GLOBAL_MEANS[0]*r)**2)
        local_total += common*clean/(1+LOCAL_MEANS[0]*r)
        previous = s
    return global_total, local_total


def four_paths(global_arms, local_arms):
    extent = 40
    geometry = [[[] for _ in range(i,extent)] for i in range(extent)]
    for i in range(extent):
        x = i+1
        for j in range(extent):
            y = j+1
            for k in range(i, extent):
                z = k+1
                exponent = (1800000*(mean_square(x)+mean_square(z))+1950000*mean_square(y)
                            +2545596*(mean_product(x,y)+mean_product(z,y))
                            +600000*mean_product(x,z)+7679999)//7680000
                if exponent <= 16000:
                    weight = 12862*exp_index(exponent)//(1000*80**3)
                    geometry[i][k-i].append((margin(j, i), margin(i, j, k), margin(j, k),
                                     weight*(1 if i == k else 2)))
    def replay(arms, steps, cap, mean):
        total = previous = F(0)
        for step in steps:
            s = TIMES[step-1]
            middle = strip_table(s, deviation(s, 2, cap))
            near = middle if cap is None else strip_table(s, deviation(s, 1, cap))
            base, extra = arms[step-1]
            bb = ab = aa = 0
            for i, row in enumerate(geometry):
                for k, cells in enumerate(row, i):
                    w = sum(((weight*near[a]//U)*middle[b]//U)*near[c]//U
                            for a,b,c,weight in cells)
                    bb += w*base[i]*base[k]
                    ab += w*(base[i]*extra[k]+extra[i]*base[k])
                    aa += w*extra[i]*extra[k]
            total += integrate_endpoint_panel(previous,s,(bb,ab,aa),mean,F(39,200),True,2)
            previous = s
        return total
    return (replay(global_arms,(3,4,5,6,7,8,10,12,14,17,23),None,GLOBAL_MEANS[1]),
            replay(local_arms,(4,6,8,10,12,15,18,23),F(3,500),LOCAL_MEANS[1]))

def long_geometry():
    extent = 32
    denominator = 12*640000
    endpoints, five, six = [], [], []
    for j in range(extent):
        y = j+1
        for k in range(extent):
            z = k+1
            center = 1950000*(mean_square(y)+mean_square(z))+2545596*mean_product(y,z)
            weight = 19301*exp_index((center+denominator-1)//denominator)//(1000*80**2)
            five.append((extent*j+k, extent*k+j, weight))
            row = []
            for i in range(extent):
                x = i+1
                numerator = (1800000*mean_square(x)+2545596*mean_product(x,y)
                             +600000*mean_product(x,z))
                exponent = (numerator+denominator-1)//denominator
                row.append((exp_index(exponent), margin(j,i), margin(i,j,k)))
            endpoints.append(row)
            for i in range(j, extent):
                x = i+1
                numerator = (1950000*(mean_square(y)+mean_square(z)+mean_square(x))
                             +2545596*(mean_product(z,y)+mean_product(z,x))
                             +600000*mean_product(y,x))
                exponent = (numerator+denominator-1)//denominator
                weight = 280379*exp_index(exponent)//(10000*80**3)
                six.append((extent*j+k, extent*i+k, margin(j,k,i), weight*(1 if j==i else 2)))
    return endpoints, five, six


def contract_arm_coefficients(geometry, arrays, near, middle):
    baseline, size_two = arrays
    bvalues, avalues = [], []
    for row in geometry:
        bsum = asum = 0
        for i, (exp, a, b) in enumerate(row):
            weight = (exp*near[a]//U)*middle[b]//U
            bsum += weight*baseline[i]//U
            asum += weight*size_two[i]//U
        bvalues.append(bsum//80)
        avalues.append(asum//80)
    return bvalues, avalues

def endpoint_coefficients_five(geometry, values):
    b, a = values
    return (sum(w*b[x]*b[y] for x,y,w in geometry),
            sum(w*(b[x]*a[y]+a[x]*b[y]) for x,y,w in geometry),
            sum(w*a[x]*a[y] for x,y,w in geometry))

def endpoint_coefficients_six(geometry, values, strips):
    b, a = values
    weighted = [(x,y,w*strips[m]//U) for x,y,m,w in geometry]
    return (sum(w*b[x]*b[y] for x,y,w in weighted),
            sum(w*(b[x]*a[y]+a[x]*b[y]) for x,y,w in weighted),
            sum(w*a[x]*a[y] for x,y,w in weighted))

def power_moment(left, right, arm_power, half=False, order=3):
    power = order+arm_power
    if not half:
        return ((1-left)**(power+1)-(1-right)**(power+1))/(power+1)
    a,b = 1-left,1-right
    low = F(isqrt(a.numerator*U**2//a.denominator),U)
    high = F(isqrt(b.numerator*U**2//b.denominator)+1,U)
    value = F(2,2*power+3)*(a**(power+1)*low-b**(power+1)*high)
    assert value > 0
    return value

def integrate_endpoint_panel(left,right,coefficients,mean,density,half=False,order=3):
    a = mean-1
    assert a > 0 and left < right
    denominator = 1+a*right
    q = F(2*order-1,2) if half else F(order-1)
    det = sqrt_lower(1/denominator**(2*order-1)) if half else 1/denominator**(order-1)
    integral = sum((F(coefficients[j],U**3)*power_moment(left,right,j,half,order)
                    for j in range(3)),F(0))
    tangent = 1+q*a*(right-left)/(2*denominator)
    return density*det*integral*tangent

def five_and_six_paths(global_arms,local_arms):
    endpoints,five,six = long_geometry()
    global_steps = (2,3,4,6,8,10,12,14,16,18,20,22,23)
    clean_steps = (2,4,6,8,10,12,14,16,18,20,22,23)
    g5=l5=g6=l6=previous_global=previous_clean=F(0)
    cap=F(7,1000)
    for step in global_steps:
        s=TIMES[step-1]
        std5=strip_table(s,deviation(s,2));std6=strip_table(s,deviation(s,3))
        g5s=contract_arm_coefficients(endpoints,global_arms[step-1],std5,std5)
        g6s=contract_arm_coefficients(endpoints,global_arms[step-1],std6,std6)
        a5=endpoint_coefficients_five(five,g5s)
        a6=endpoint_coefficients_six(six,g6s,std6)
        g5+=integrate_endpoint_panel(previous_global,s,a5,GLOBAL_MEANS[2],F(1804,10000))
        g6+=integrate_endpoint_panel(previous_global,s,a6,GLOBAL_MEANS[3],F(429,2500),True)
        previous_global=s
        if step not in clean_steps:
            continue
        near,mid,center=(strip_table(s,deviation(s,d,cap))for d in(1,2,3))
        shared=contract_arm_coefficients(endpoints,local_arms[step-1],near,mid)
        b5=endpoint_coefficients_five(five,shared)
        b6=endpoint_coefficients_six(six,shared,center)
        l5+=integrate_endpoint_panel(previous_clean,s,b5,LOCAL_MEANS[2],F(1804,10000))
        l6+=integrate_endpoint_panel(previous_clean,s,b6,LOCAL_MEANS[3],F(429,2500),True)
        previous_clean=s
    return (g5,l5),(g6,l6)


def certificate(progress=None):
    normalizer_check()
    mean_check()
    default = tuple(arm_arrays(deviation))
    high = tuple(arm_arrays(lambda s, d: deviation(s, d, high=True)))
    def clean(cap):
        return tuple(arm_arrays(lambda s, d: deviation(s, d, cap if d == 0 else None)))
    local6, local7 = clean(F(3, 500)), clean(F(7, 1000))
    results = (three_paths(default, high, local6), four_paths(default, local6),
               *five_and_six_paths(default, local7))
    rows = []
    for d, (global_sum, local_sum), g, l in zip(range(3, 7), results, GLOBAL_FLOORS, LOCAL_FLOORS, strict=True):
        assert global_sum > g and local_sum > l > g, d
        row = {"distance": d, "global_exact_lower_sum": str(global_sum), "global_floor": str(g),
               "local_exact_lower_sum": str(local_sum), "local_floor": str(l), "local_increment": str(l-g)}
        rows.append(row)
        if progress is not None:
            progress(row)
    return {"scope": "global and clean-window D3-D6 Price arithmetic", "rows": rows}


long_paths = SimpleNamespace(certificate=certificate, TIMES=TIMES, arm_arrays=arm_arrays, deviation=deviation, sqrt_lower=sqrt_lower, sqrt_upper=sqrt_upper, arm_geometry=arm_geometry, strip_table=strip_table, three_geometry=three_geometry, GLOBAL_MEANS=GLOBAL_MEANS, LOCAL_MEANS=LOCAL_MEANS)


"""Recompute the clean D2 baseline and the three retained size-two arm payments."""


d2_baseline = baseline2


CAP = F(3, 500)
INCIDENT = F(1, 50)
PIVOT = F(2878, 10000)
MEAN = F(1313, 1000)
CLEAN_MEAN = F(20981, 20000)
BASELINE_GLOBAL_FLOOR = F(30239457, 10**10)
BASELINE_LOCAL_FLOOR = F(42902073, 10**10)
ARM_FLOORS = (F(17745649, 10**10), F(27670235, 10**10), F(1257887, 400000000))


def arms2_mean_check():
    root10 = F(3162277, 10**6)
    assert root10**2 < 10
    assert F(571, 11)-16*root10 < MEAN
    # The sharp clean-window trace bound has one interior mode.
    root5_upper, root2_lower = F(559017, 250000), F(1414213, 1000000)
    assert root5_upper**2 > 5 and root2_lower**2 < 2
    assert 1+48*CAP*(4*root5_upper-5*root2_lower)/11 < CLEAN_MEAN


def clean_baseline():
    n, width, time = d2_baseline.CELLS, d2_baseline.SPATIAL, d2_baseline.TIME
    xd = [45834*exp_lower(F(33*(i+1)**2, 2*width**2))//10000 for i in range(n)]
    zd = [13819*exp_lower(F(3*(j+1)**2, 2*width**2))//10000 for j in range(n)]
    total = previous = F(0)
    for step in range(1, time):
        s = F(step, time)**2
        middle = 2*SQRT2_LO/deviation(s, 1, CAP)
        root = 1/deviation(s, 0, CAP)
        middle_strip = d2_baseline.strip_values(middle.numerator, middle.denominator, s)
        root_strip = d2_baseline.strip_values(root.numerator, root.denominator, s)
        row = 0
        for i, density in enumerate(xd):
            convolution = sum(zd[j]*root_strip[i-j-1] for j in range(i))//(U*width)
            weight = density*middle_strip[i]//U
            row += (weight*(convolution**2//U)//U)//width
        determinant = (1-s)*sqrt_lower(1/(1+CLEAN_MEAN*s/(1-s)))
        total += PIVOT*(s-previous)*determinant*F(row, U)
        previous = s
    return total


def clean_deviation(s, distance):
    return deviation(s, distance, CAP if distance == 0 else None)


def stratified_deviation(s, distance):
    if distance == 1:
        variance = C_HI*INCIDENT-s*INCIDENT**2/2
        assert 0 < variance < C_HI**2*(1-s/2)
        return sqrt_upper(variance)
    return clean_deviation(s, distance)


def values(s, row):
    baseline, arms = row
    return [(F(b, U), (1-s)*F(a, U)) for b,a in zip(baseline, arms, strict=True)]


def joint_arm_mass(old, high, gain):
    b, a = old
    bh, ah = high
    return min(gain**2*(2*b*a+a*a), gain*(b*ah+a*bh+a*ah), 2*bh*ah+ah*ah)


def clean_arm_mass(side):
    b, a = side
    return 2*b*a+a*a


def two_density():
    return [45834*exp_lower(F(33, 2)*F(x+1, 80)**2)//(10000*80) for x in range(64)]


def size_two_payments():
    density = two_density()
    high = arm_arrays(lambda s,d: deviation(s,d,high=True))
    baseline = arm_arrays(deviation)
    local = arm_arrays(clean_deviation)
    extra = arm_arrays(stratified_deviation)
    global_total = local_total = stratified_total = previous = F(0)
    for s, old, h, l, st in zip(TIMES, baseline, high, local, extra, strict=True):
        old, h, l, st = (values(s, x) for x in (old, h, l, st))
        gain = sqrt_lower(1/(1-s/3))
        global_middle = deviation(s, 1)/(2*SQRT2_LO)
        local_middle = deviation(s, 1, CAP)/(2*SQRT2_LO)
        global_mass = local_mass = stratified_mass = F(0)
        for x,w in enumerate(density):
            global_mass += F(w*strip(x,s,global_middle), U**2)*joint_arm_mass(old[x],h[x],gain)
            local_weight = F(w*strip(x,s,local_middle), U**2)
            local_mass += local_weight*clean_arm_mass(l[x])
            stratified_mass += local_weight*clean_arm_mass(st[x])
        assert stratified_mass >= local_mass
        common = PIVOT*(s-previous)*(1-s)
        r = s/(1-s)
        global_total += common*sqrt_lower(1/(1+MEAN*r))*global_mass
        factor = common*sqrt_lower(1/(1+CLEAN_MEAN*r))
        local_total += factor*local_mass
        stratified_total += factor*stratified_mass
        previous = s
    return global_total, local_total, stratified_total


def arms2_certificate():
    arms2_mean_check()
    baseline = clean_baseline()
    assert baseline > BASELINE_LOCAL_FLOOR > BASELINE_GLOBAL_FLOOR
    exact = size_two_payments()
    assert all(a > b for a,b in zip(exact, ARM_FLOORS, strict=True))
    g, l, s = ARM_FLOORS
    return {"scope": "D2 clean no-arm baseline and size-two closed-arm Price arithmetic",
            "baseline_local_exact_lower_sum": str(baseline), "baseline_local_floor": str(BASELINE_LOCAL_FLOOR),
            "baseline_local_increment": str(BASELINE_LOCAL_FLOOR-BASELINE_GLOBAL_FLOOR),
            "arm_exact_lower_sums": list(map(str, exact)), "arm_floors": list(map(str, ARM_FLOORS)),
            "clean_arm_increment": str(l-g), "stratified_arm_increment": str(s-l)}


arms2 = SimpleNamespace(certificate=arms2_certificate, CAP=CAP, CLEAN_MEAN=CLEAN_MEAN, MEAN=MEAN, PIVOT=PIVOT, clean_deviation=clean_deviation, two_density=two_density)


"""Recompute the global and clean-window size-four arm payments at D2 and D3."""




STEPS = (1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20, 22, 23)
PANEL_TIMES = tuple(F(j, 24)**2 for j in STEPS)
FLOORS = ((F(37143, 625000000), F(104039, 10000000000)), (F(130757, 1250000000), F(45859, 2500000000)))

def panel_arm_arrays(incoming):
    return arm_arrays(incoming, PANEL_TIMES)


@lru_cache(maxsize=1)
def continuation():
    # Adjacent contracts over the previous key. These arms sum over the next.
    _, kernel, _ = arm_transition()
    return [tuple(zip(*rows, strict=True)) for rows in kernel]


def prepend(suffix, next_strip, kernel):
    result = [[0]*48 for _ in range(48)]
    for current, rows in enumerate(kernel):
        weights = [v*t//U for v,t in zip(suffix[current], next_strip, strict=True)]
        for previous, row in enumerate(rows):
            result[previous][current] = sum(map(mul, row, weights))//U
    return result


def size_four_arms(incoming):
    roots, first, lower, dz = arm_geometry()
    kernel = continuation()
    for s in PANEL_TIMES:
        sd = [incoming(s, d) for d in range(5)]
        ell = [sum((j+1)*Q_HI**j*sd[d-1-j] for j in range(d))/2 for d in range(1, 6)]
        strips = [[strip(k, s, e) for k in range(48)] for e in ell[:4]]
        terminal_strip = [strip(k, s, ell[4], 2) for k in range(len(dz))]
        last = divided_convolution(dz, terminal_strip)
        suffix = [[last[m] for m in row] for row in lower]
        suffix = prepend(prepend(suffix, strips[3], kernel), strips[2], kernel)
        inner = [sum((a*b//U)*c for a,b,c in zip(row, strips[1], tail, strict=True))//U
                 for row,tail in zip(first, suffix, strict=True)]
        yield [sum((a*b//U)*c for a,b,c in zip(row, strips[0], inner, strict=True))//U for row in roots]


def size_four_moment(left, right, twice_power):
    a, b = 1-left, 1-right
    if twice_power % 2 == 0:
        power = twice_power//2+1
        return (a**power-b**power)/power
    low = F(isqrt(a.numerator*U**2//a.denominator), U)
    high = F(isqrt(b.numerator*U**2//b.denominator)+1, U)
    result = F(2, twice_power+2)*(a**(twice_power//2+1)*low
                                - b**(twice_power//2+1)*high)
    assert result > 0
    return result

def integrate_panel(left, right, coefficients, mean, density, distance):
    q, a = F(distance-1, 2), mean-1
    assert a > 0 and left < right
    denominator = 1+a*right
    determinant = sqrt_lower(1/denominator) if distance == 2 else 1/denominator
    tangent = 1+q*a*(right-left)/(2*denominator)
    integral = sum((F(value, U**3)*size_four_moment(left, right, distance+5+2*j)
                    for j, value in enumerate(coefficients)), F(0))
    return density*determinant*tangent*integral

def integrate(incoming, cap, mean2, mean3):
    density2, geometry3 = two_density(), three_geometry()
    total2 = total3 = previous = F(0)
    for s, (baseline, arm2), arm4 in zip(
            PANEL_TIMES, panel_arm_arrays(incoming), size_four_arms(incoming), strict=True):
        middle = deviation(s, 1, cap)
        strips2 = [strip(x, s, middle/(2*SQRT2_LO)) for x in range(64)]
        strips3 = strip_table(s, middle)
        c2, c3 = [0, 0, 0], [0, 0, 0]
        for w, t, B, A, C in zip(density2, strips2, baseline, arm2, arm4, strict=True):
            w = w*t//U
            c2[0] += w*2*B*C
            c2[1] += w*2*A*C
            c2[2] += w*C*C
        for i, j, x, y, w in geometry3:
            w = (w*strips3[x]//U)*strips3[y]//U
            c3[0] += w*(baseline[i]*arm4[j]+arm4[i]*baseline[j])
            c3[1] += w*(arm2[i]*arm4[j]+arm4[i]*arm2[j])
            c3[2] += w*arm4[i]*arm4[j]
        total2 += integrate_panel(previous, s, c2, mean2, PIVOT, 2)
        total3 += integrate_panel(previous, s, c3, mean3, F(225, 1000), 3)
        previous = s
    return total2, total3

def size_four_certificate():
    global_sums = integrate(deviation, None, MEAN, GLOBAL_MEANS[0])
    local_sums = integrate(clean_deviation, CAP, CLEAN_MEAN, LOCAL_MEANS[0])
    rows = []
    for i, d in enumerate((2, 3)):
        g, l = global_sums[i], local_sums[i]
        gf, lf = FLOORS[0][i], FLOORS[1][i]
        assert g > gf > 0 and l > lf > gf
        rows.append({"distance": d, "global_exact_lower_sum": str(g), "global_floor": str(gf),
                     "local_exact_lower_sum": str(l), "local_floor": str(lf), "local_increment": str(lf-gf)})
    return {"scope": "D2 and D3 size-four arm payments on sixteen positive-power panels", "rows": rows}


size_four = SimpleNamespace(certificate=size_four_certificate)


"""Exact one-band composition of the portable numerical replays."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt

POOL_SEAM = F('28909005122956699438697/36057402177020999792110000')
UPPER = F('2534356112695870650792437/3135474478647834971934850000')
PATH_UPPER = F('1633094581/12500000000')

CRITICAL_KAPPA = F(59,1000)
SCALE = 10000000000000000
LOW, HIGH = F('1/2'), F('2')

def bilateral033():
    S = 10**12
    def root(n, d=1):
        return isqrt(n*S*S//d)
    def up(n, d):
        return (n+d-1)//d
    l2, l3, l6 = (root(n) for n in (2, 3, 6))
    u2 = l2+1
    cl = 9*l6+12*l3-15*u2-20*S
    cu = 9*(l6+1)+12*(l3+1)-15*l2-20*S
    assert 0 < cl < cu < 2*S
    def phi(n, d, upper=False):
        z = root(n, d)
        return z+1-l2*(n+d)//(3*d) if upper else max(0, z-up(u2*(n+d), 3*d))
    def gap(n, d):
        return max(0, root(5*d-2*n, 2*d)+root(5*n*n-2*n*d, 2*d*d)
                   -up(cu*(root(n, d)+1), S)
                   -up(u2*(5*S-2*cl)*(n+d), 6*S*d))
    grid = [gap(50+j, 100) for j in range(151)]
    phi_grid = [phi(i, 5000) for i in range(2500, 10001)]
    def B(n, d):
        m = (200*d-100*n)//d
        if d == 5000:
            z = min(grid[j]+phi_grid[10000-n-50*j] for j in range(m+1))
        else:
            z = min(grid[j]+phi((250-j)*d-100*n, 100*d) for j in range(m+1))
        return max(0, min(z, gap(5*d-2*n, 2*d))-S//40000)
    margins = []
    for i in range(2750, 5000):
        z = max(0, B(i, 5000)-S//20000)
        z += i*max(0, B(5000, i+1)-up(5000*S, 4*i*(i+1)))//5000
        margins.append(z-up(33*phi(i+1, 5000, True), 100))
    minimum = min(margins)
    assert minimum > 3*S//10**5
    def shape(x):
        return (F(3, 2)-x)**2/(F(5, 2)-x)
    def product(a, b):
        h = 0 if a <= F(3, 2) <= b else min(shape(a), shape(b))
        return (F(root(a.numerator, a.denominator), S)+F(l2, 3*S)*(1+a))*h
    values = []
    for i in range(200, 220):
        a, b = F(i, 400), F(i+1, 400)
        values.append(F(16, 15)*(product(a, b)+product(1/b, 1/a))*S/(2*u2))
    analytic = min(values)
    assert analytic > F(40637, 100000) > F(33, 100)
    return dict(factor='33/100', cells=2250, direct_margin_floor=str(F(minimum, S)),
                cell=str(F(2750+margins.index(minimum), 5000)), analytic_min=str(analytic))

@lru_cache(maxsize=None)
def sqrt_bounds(x):
    lower = F(isqrt((x*SCALE**2).__floor__()), SCALE)
    upper = lower+F(1, SCALE)
    assert lower**2 <= x < upper**2
    return lower, upper

def parabola(t):
    return (t-LOW)*(HIGH-t)

@lru_cache(maxsize=None)
def phi_bounds(lo, hi):
    assert LOW <= lo <= hi <= HIGH
    best = max(lo, min(hi, F(5, 4)))
    root2_lo, root2_hi = sqrt_bounds(F(2))
    den_lo = sqrt_bounds(lo)[0]+root2_lo*(1+lo)/3
    den_hi = sqrt_bounds(hi)[1]+root2_hi*(1+hi)/3
    return F(2, 9)*min(parabola(lo), parabola(hi))/den_hi, F(2, 9)*parabola(best)/den_lo

@lru_cache(maxsize=None)
def star_residual_at(t):
    low2, high2 = sqrt_bounds(F(2))
    low3, high3 = sqrt_bounds(F(3))
    low6, high6 = sqrt_bounds(F(6))
    low = 9*low6+12*low3-15*high2-20+CRITICAL_KAPPA
    high = 9*high6+12*high3-15*low2-20+CRITICAL_KAPPA
    assert 0 < low <= high < F(5,2)
    assert high < 8*F(1,2)*low2/2
    return (sqrt_bounds(F(5,2)-t)[0]+sqrt_bounds(F(5,2)*t*t-t)[0]
            -high*sqrt_bounds(t)[1]+high2*(low-F(5,2))*(1+t)/3)


@lru_cache(maxsize=None)
def residual_lower(lo, hi):
    t = max(lo, min(hi, F(1)))
    coefficient = F(16, 15)*(t-1)**2/t-CRITICAL_KAPPA
    low, high = phi_bounds(lo, hi)
    quadratic = coefficient*(low if coefficient >= 0 else high)
    chord = min(star_residual_at(lo),star_residual_at(hi))-(hi-lo)**2/4
    return max(quadratic,chord)

def cover():
    remaining = [(LOW, HIGH, LOW, HIGH, 0)]
    discarded = accepted = maximum_depth = 0
    covered_area = F(0)
    while remaining:
        pl, ph, ql, qh, depth = remaining.pop()
        rl, rh = max(LOW, 3-ph-qh), min(HIGH, 3-pl-ql)
        if rl > rh:
            discarded += 1
            continue
        lower = residual_lower(pl, ph)+residual_lower(ql, qh)+phi_bounds(rl, rh)[0]
        if lower >= 0:
            accepted += 1
            maximum_depth = max(depth, maximum_depth)
            covered_area += (ph-pl)*(qh-ql)
            continue
        assert depth < 36, (pl, ph, ql, qh, lower)
        if ph-pl >= qh-ql:
            mid = (pl+ph)/2
            remaining.extend(((pl, mid, ql, qh, depth+1), (mid, ph, ql, qh, depth+1)))
        else:
            mid = (ql+qh)/2
            remaining.extend(((pl, ph, ql, mid, depth+1), (pl, ph, mid, qh, depth+1)))
    assert covered_area >= F(9, 8)
    return {"accepted_rectangles": accepted, "outside_rectangles": discarded,
            "maximum_depth": maximum_depth, "accepted_rectangle_area": str(covered_area)}

def population(q, u, d, p, w, A, beta, *, reverse=F(0),
               type20_payment=F(7, 10), pair_factor=F(8, 5), high_aware=F(0),
               endpoint_adjustments=(), radial_j_debit=None, joint_j_credit=F(0), v11=None):
    assert 0 <= high_aware <= u
    c = (pair_factor-F(3, 2))*A/4
    h = (pair_factor-1)*A/4
    if v11 is None:
        v11 = ((F(3, 2)+beta/2)*A+2*c-92*p-4*u+14*high_aware-F(232, 3)*w-42*reverse)/26
    v9 = p-u-v11
    shared = 13*v9+17*v11+12*u-high_aware+F(32, 3)*w+6*reverse
    assert v9 > 0 and v11 > 0
    assert (F(1, 8)+beta/4)*A < shared < (F(1, 4)+beta/4)*A
    g = F(2, 3)*(shared-c-beta*A/4)
    assert g > 0
    B = c+g/2
    if radial_j_debit is None:
        radial_j_debit = A*(10*beta-4)/8
    assert 0 <= joint_j_credit <= 2*min(high_aware, reverse)
    overlap = {
        "U2": 12*v9+16*v11+12*u+6*reverse-beta*A/2,
        "U3": 30*v9+40*v11+24*u+12*reverse-3*beta*A/4,
        "J": 4*v9+6*v11+2*high_aware+joint_j_credit-radial_j_debit,
    }
    assert all(value >= 0 for value in overlap.values())
    coefficients = {
        (0, 0): F(0), (0, 1): A*(F(1, 2)+beta/4)+B, (0, 2): A/2-g, (0, 3): 3*A/4-3*B,
        (1, 0): 3*A/4+h+shared, (1, 1): 3*A/4+h-B, (1, 2): A-h-2*B,
        (2, 0): type20_payment*A, (2, 1): 5*A/4-2*h-B, (3, 0): 3*A/2-3*h,
    }
    for kind, adjustment in endpoint_adjustments:
        coefficients[kind] += adjustment
    c7, w9, w11 = (0, 18, 12, 6), (0, 35, 22, 9), (0, 45, 28, 11)
    cW, cR = (0, F(88, 3), F(56, 3), 8), (0, 36, 24, 12)
    cH, c_reverse = (0, 14, 10, 6), (0, 16, 10, 4)
    slacks = {}
    for (a, b), value in coefficients.items():
        degree = a+b
        terminal_debit = cR[degree]+8*(3-degree)*(a > 0)
        alternate_debit = w9[degree]-2*(3-b)*(a == 0 and b > 0)
        required = (q*c7[a]+v9*w9[degree]+v11*w11[degree]+w*cW[degree]
                    +u*terminal_debit+high_aware*(alternate_debit-terminal_debit)
                    +d*cH[a]+reverse*c_reverse[degree])
        slacks[f"{a},{b}"] = value-required
    assert all(value >= 0 for value in slacks.values())
    credit = 6*q+9*v9+11*v11+12*u-3*high_aware+8*w+6*d+4*reverse
    return {"q": str(q), "u": str(u), "d": str(d), "v9": str(v9), "v11": str(v11),
            "shared": str(shared), "credit": str(credit),
            "type_slacks": {key: str(value) for key, value in slacks.items()},
            "overlap_slacks": {key: str(value) for key, value in overlap.items()}}

WEIGHTS = (139357900, 89691594, 36451375, 100, 100, 1446536023, 337910715)
SAVING = F('1981610720/100000000000000')
GAMMA = F('500000000/3827670601')


ETA = F(323, 400000)
EXTENT = F(9, 400)
SIGMA, ALPHA, LAM, SHARE = F(2021, 1250), F(21, 2), F(141421, 50000), F(13, 2)
MARGIN = F(1, 10**12)
H = (F(27, 10000), F(3, 500), F(7, 1000))
WEAK_WEIGHT = F(12719, 10**12)


def resource_table(adjacent, baseline2, arms2, long_paths, size_four):
    assert [row["distance"] for row in long_paths["rows"]] == [3, 4, 5, 6]
    assert [row["distance"] for row in size_four["rows"]] == [2, 3]
    paths = long_paths["rows"]
    four2, four3 = size_four["rows"]
    global_values = (
        F(adjacent["global_K0_payment"]),
        2*LAM*F(baseline2["joint_usable_floor"]),
        2*LAM*F(arms2["arm_floors"][0]),
        *(mult*F(row["global_floor"]) for mult,row in zip((19, 18*LAM, 120, 92*LAM), paths, strict=True)),
        2*LAM*F(four2["global_floor"]),
        19*F(four3["global_floor"]),
    )
    # Keep the previously spent rounding down of the original pair credit.
    k0 = sum(global_values)-F(20789, 25000000000000000000000)
    transfer = F(adjacent["global_to_clean_price_transfer"])
    assert transfer == F(71931, 10**10)
    assert k0+transfer == F(1070703902702099979211, 25000000000000000000000)
    caps = (
        F(adjacent["local_adjacent_raw_cap"]),
        F(arms2["baseline_local_increment"])+F(arms2["clean_arm_increment"])+F(four2["local_increment"]),
        F(paths[0]["local_increment"])+F(four3["local_increment"]),
        *(F(row["local_increment"]) for row in paths[1:]),
        F(arms2["stratified_arm_increment"]),
    )
    assert all(c > 0 for c in caps)
    assert caps[0]-transfer == F(13702899, 2500000000)
    common = (1, 2*LAM, 19, 18*LAM, 120, 92*LAM, 2*LAM)
    dirty = ((29,29,29), (0,18*LAM,18*LAM), (0,327,327),
             (0,522*LAM,522*LAM), (0,0,6040), (0,0,7788*LAM), (0,18*LAM,42*LAM))
    names = ("adjacent", "distance_two_nine_edge", "distance_three", "distance_four",
             "distance_five", "distance_six", "stratified_distance_two")
    return k0, tuple(zip(names, caps, common, dirty, strict=True))


def clean_window_budget(k0, groups):
    def resources():
        k_full = k0+sum(c*m for _,c,m,_ in groups)
        return k0, k_full, H, None, groups

    def class_collateral(eta, threshold, k):
        K0, _, H, _, _ = resources()
        assert 0 <= k <= 4
        return K0*max(F(0), SIGMA*threshold-k*H[0]-(4-k)*eta)/2

    def isolated_first_class_bonus(eta, adjacent_spend):
        eta, adjacent_spend = F(eta), F(adjacent_spend)
        K0, _, H, _, groups = resources()
        I1 = K0*(SIGMA*H[0]-4*eta)/2
        assert 0 < eta < H[0] and I1 >= 0 and 0 <= adjacent_spend <= groups[0][1]
        return min(I1, SHARE*eta*adjacent_spend)

    def dyad_payment(eta, margin=MARGIN):
        K0, _, H, _, _ = resources()
        assert margin > 0
        return max(F(0), K0*max(F(0), F(63, 50)*H[0]-F(41, 10)*F(eta))/2-margin)

    def extended_source(eta):
        eta = F(eta)
        assert POOL_SEAM <= eta <= UPPER
        K0, Kfull, H, _, groups = resources()
        retained = [cap for _, cap, _, _ in groups]
        S = ALPHA*eta*retained[0]
        T = eta*(LAM*(6*retained[1]+189*retained[3])+F(231, 2)*retained[2])
        t5, t6 = 2220*eta*retained[4], 2886*LAM*eta*retained[5]
        I = tuple(K0*(SIGMA*h-4*eta)/2 for h in H)
        J = {j: tuple(class_collateral(eta, H[j-1], k) for k in range(5)) for j in (2, 3)}
        dyad = {j: K0*((SIGMA-1)*H[j-1]-3*eta)/2 for j in (2, 3)}
        first_bonus = isolated_first_class_bonus(eta, retained[0])
        baseline = tuple(K0*(h-eta)-eta*sum(dirty[j]*x for x, (_, _, _, dirty)
                         in zip(retained, groups, strict=True))+min(S, I[j]) for j, h in enumerate(H))
        pool = baseline[0]+S-I[0]-MARGIN
        P = pool+S
        t = {j: MARGIN-baseline[j-1]-T for j in (2, 3)}
        v = max((t[2]+t[3])/2, t[3]-t5)
        gamma = dyad_payment(eta)
        branch = {"I1_positive": I[0], "S_above_I1": S-I[0], "I2_above_S": I[1]-S,
                  "I3_above_I2": I[2]-I[1], "isolated_overlap_above_I1": F(13, 2)*eta*retained[0]-I[0],
                  "top_above_middle_bonus": t[3]-t[2], "middle_bonus_positive": t[2],
                  "mixed_maximum_is_top_minus_D5": t[3]-t5-(t[2]+t[3])/2,
                  "usable_class_one_dyad_payment_positive": gamma}
        for j in (2, 3):
            for k in range(5):
                raw = SIGMA*H[j-1]-k*H[0]-(4-k)*eta
                branch[f"class_{j}_degree_{k}_clamp_sign"] = -raw if (j, k) == (2, 4) else raw
        assert min(branch.values()) > 0 and v == t[3]-t5 and first_bonus == I[0]
        assert pool <= 0
        v -= 6*LAM*eta*retained[6]
        guards = {"isolated_first_class": baseline[0]+first_bonus-MARGIN,
                  "pure_low_two_edge": 2*pool+gamma,
                  "pure_low_at_least_three": 3*pool+S, "pure_low_slope": P,
                  "middle_singleton": I[1]-S-T-t[2], "top_singleton": I[2]-S-T-t[3],
                  "pure_two_middle_dyad": dyad[2]/2-t[2],
                  "pure_top_middle_dyad": dyad[3]+18*LAM*eta*retained[6]-t[2]-t[3],
                  "pure_two_top_dyad": t5+t6+dyad[3]/2-t[3],
                  "pure_component_at_least_three": (S+T)/3-max(t.values()),
                  "individual_bonus_below_S": S-max(t.values()), "mixed_component_below_S": S-v,
                  "l1_q0_r1": P/2-v, "l1_q0_r2": (pool+3*S)/4-v,
                  "l2_nonsingleton": 2*(S-v)+2*pool}
        for j in (2, 3):
            singleton = J[j][1]-t[j]
            guards[f"l1_q1_r0_class_{j}"] = pool+singleton-T
            guards[f"l1_q2_r0_class_{j}"] = P/2+singleton-T
            guards[f"l1_q1_r1_class_{j}"] = pool+2*S-T+singleton-2*v
            for k in (1, 2):
                guards[f"l2_singleton_class_{j}_degree_{k}"] = S+J[j][k]-T-t[j]+2*pool
        strata = []
        for lower_edges, max_degree, beta in ((3, 3, pool/2+S/6), (5, 4, pool/2+3*S/10)):
            assert beta > 0
            values = {j: [S+J[j][k]+k*beta-T-t[j] for k in range(1, max_degree+1)] for j in (2, 3)}
            assert all(min(slacks) > 0 for slacks in values.values())
            strata.append({"lower_edge_stratum": str(lower_edges) if lower_edges < 4 else "at_least_5",
                           "beta": str(beta), "class_two_slacks": list(map(str, values[2])),
                           "class_three_slacks": list(map(str, values[3]))})
        assert min(guards.values()) >= 0
        coupled = (baseline[0]+first_bonus, baseline[1]+T+t[2], baseline[2]+T+t[3])
        assert min(coupled) >= MARGIN and coupled[1] == coupled[2] == MARGIN
        K = K0+sum(common*x for x, (_, _, common, _) in zip(retained, groups, strict=True))
        assert K == Kfull
        return {"scope": "exact full-capacity negative-pool source continuation; no event replay",
                "eta": str(eta), "scale": str(eta/F(19, 40000)), "reserve_floor": str(MARGIN),
                "regime": "negative_pool_class_one_dyads_paid", "capped_K": str(K), "eta_times_K": str(eta*K),
                "retained_spends": {name: str(x) for x, (name, _, _, _) in zip(retained, groups, strict=True)},
                "withdrawals": {name: "0" for name, _, _, _ in groups},
                "adjacent_overlap_saving": str(S), "isolated_collateral": list(map(str, I)),
                "isolated_first_class_bonus": str(first_bonus),
                "class_specific_degree_collateral": {str(j): list(map(str, J[j])) for j in (2, 3)},
                "pure_dyad_collateral": {str(j): str(dyad[j]) for j in (2, 3)},
                "usable_pure_class_one_dyad_payment": str(gamma), "mixed_pool": str(pool),
                "baseline_class_reserves": list(map(str, baseline)), "extra_high_class_reserve": str(T),
                "D5_top_overlap": str(t5), "D6_pure_top_dyad_overlap": str(t6),
                "extra_middle_class_reserve": str(t[2]), "extra_top_class_reserve": str(t[3]),
                "mixed_nonsingleton_effective_bonus": str(v),
                "coupling_guards": {name: str(value) for name, value in guards.items()}, "degree_strata": strata,
                "branch_guards": {name: str(value) for name, value in branch.items()},
                "exceptional_edge_reserves": list(map(str, coupled)), "reserves": list(map(str, coupled)),
                "capacity_upper_bound": str(Kfull),
                "D6_overlap_scope": "only whole E1 components consisting of two top edges"}

    def continuation_point(eta):
        point = extended_source(F(eta))
        S, pool = F(point["adjacent_overlap_saving"]), F(point["mixed_pool"])
        T, v = F(point["extra_high_class_reserve"]), F(point["mixed_nonsingleton_effective_bonus"])
        t = {2: F(point["extra_middle_class_reserve"]), 3: F(point["extra_top_class_reserve"])}
        J = {j: list(map(F, point["class_specific_degree_collateral"][str(j)])) for j in (2, 3)}
        beta = pool/2+S/4
        mu = {j: [S+J[j][k]+k*beta-T-t[j] for k in range(1, 5)] for j in (2, 3)}
        loss = max(F(0), -mu[2][3], -mu[3][3])
        compensation = min(4*beta, 2*(S-v)+beta, *(mu[j][k] for j in (2, 3) for k in range(3)))
        assert beta > 0 and compensation-loss > 0
        return {**point, "class_one_dyad_path_ratio": "63/100",
                "special_l4_control": {"beta": str(beta), "maximum_negative_degree4_cost": str(loss),
                                       "minimum_compensation": str(compensation), "slack": str(compensation-loss),
                                       "degree4_singletons_at_most": 1,
                                       "reason": "two degree4 singletons force a whole K4 with equal source excesses"}}

    return continuation_point(ETA)


TARGET_ALLOWANCE = F('600620827548562093285153/10000000000000000000000000000')
def stronger_resource_table(adjacent, baseline2, arms2, long_paths, size_four):
    k0, old_groups = resource_table(adjacent, baseline2, arms2, long_paths, size_four)
    groups = list(old_groups)
    rows = (*size_four["rows"], *long_paths["rows"][2:])
    floors = ((602786, 1060433), (106313, 187286), (100300, 164138), (17370, 28764))
    for row, (g, l) in zip(rows, floors, strict=True):
        slot = row["distance"]-1
        name, cap, common, dirty = groups[slot]
        new_global, new_local = F(g, 10**10), F(l, 10**10)
        old_global, old_local = F(row["global_floor"]), F(row["local_floor"])
        assert old_global < new_global < F(row["global_exact_lower_sum"])
        assert old_local < new_local < F(row["local_exact_lower_sum"])
        k0 += common*(new_global-old_global)
        groups[slot] = name, cap+(new_local-new_global)-(old_local-old_global), common, dirty
    return k0, groups

def cell_mean_resource_table(adjacent, baseline2, arms2, long_paths, size_four):
    global_credit, groups = stronger_resource_table(adjacent, baseline2, arms2, long_paths, size_four)
    # Floor increases at scale 10^-10. Local rows replace global rows.
    global_base, _, clean_base, _ = adjacent["components"]
    cells = [
        (0, global_base["exact_lower_sum"], clean_base["exact_lower_sum"],
         F(global_base["usable_floor"]), F(clean_base["usable_floor"]), 1869306, 4361809),
    ]
    for row, dg, dl in zip(long_paths["rows"][:2], (399518, 110407), (643497, 182424), strict=True):
        cells.append((row["distance"]-1, row["global_exact_lower_sum"], row["local_exact_lower_sum"],
                      F(row["global_floor"]), F(row["local_floor"]), dg, dl))
    for row, g, l, dg, dl in zip(
            (*size_four["rows"], *long_paths["rows"][2:]),
            (602786,106313,100300,17370), (1060433,187286,164138,28764),
            (1,10362,21310,4944), (1,17968,34662,3105), strict=True):
        cells.append((row["distance"]-1, row["global_exact_lower_sum"], row["local_exact_lower_sum"],
                      F(g,10**10), F(l,10**10), dg, dl))
    for slot, exact_global, exact_local, old_global, old_local, dg, dl in cells:
        increase_global, increase_local = F(dg,10**10), F(dl,10**10)
        assert old_global+increase_global < F(exact_global)
        assert old_local+increase_local < F(exact_local)
        name, cap, common, dirty = groups[slot]
        global_credit += common*increase_global
        groups[slot] = (name, cap+increase_local-increase_global, common, dirty)
    return global_credit, groups


def final_allocation(adjacent, baseline2, arms2, long_paths, size_four, d7, base, signed, nine, outer):
    bilateral, star_cover = bilateral033(), cover()
    k = F(bilateral["factor"])
    assert k == F(33, 100)
    assert star_cover["accepted_rectangles"] == 334 and F(117,2000) < CRITICAL_KAPPA
    k0, groups = cell_mean_resource_table(adjacent, baseline2, arms2, long_paths, size_four)
    budget = clean_window_budget(k0, groups)
    assert F(budget["eta"]) == ETA
    A = F(budget["eta_times_K"])+ETA*F(d7["global_K_increment"])
    assert A == TARGET_ALLOWANCE

    q, u, v9, v11, gap, X, E = (F(value, 10**14) for value in WEIGHTS)
    p, high_aware = u+v9+v11, u-gap
    assert F(signed["eta"]) == ETA and F(signed["threshold_abs_max"]) == EXTENT
    assert sorted(row["side"] for row in signed["records"]) == [-1, 1]
    assert 0 < q < min(F(row["usable_floor"]) for row in signed["records"])
    assert len(nine) == 2 and {row["phase"] for row in nine} == {"strict", "relaxed"}
    assert all(F(row["eta"]) == ETA and F(row["threshold_abs_max"]) == EXTENT for row in nine)
    caps = {row["phase"]: F(row["used_floor"]) for row in nine}
    assert 0 < high_aware < u < caps["relaxed"] and p < caps["strict"]
    assert F(base["weak"]["usable_floor"]) >= WEAK_WEIGHT

    rho, omega = 2*X/A, F(1, 4)
    a0 = SIGMA-F(3, 2)
    c, d = a0/4, omega*(SIGMA-1)/2
    shared = 13*v9+17*v11+11*u+gap+F(32, 3)*WEAK_WEIGHT
    s, epsilon = shared/A, E/A
    assert F(2, 5) < rho < F(1, 2)
    assert F(1, 8)+rho/4 < s < F(1, 4)+rho/4 and s >= omega/2
    assert F(1444, 111375) <= epsilon < 2*rho*F(117,2000)
    assert epsilon <= min(rho*a0, (2-SIGMA)/2)
    g = F(2, 3)*(s-c-rho/4)
    B = c+g/2
    old = 1/(1+16*(omega/2-c)/k)
    joint = (B+k/16)/(omega/2+k/16)
    assert g >= 0 and 0<old<joint<1
    share = old+F(908475033730, 10**12)*(joint-old)
    assert min((F(1, 2)+rho-omega)*a0, s-omega/2) > d*share
    pay = A*(d*share+(k/8)*(1-share))
    assert E <= (5-3*SIGMA)*A/4+pay
    beta = F(55903822, 10**14)
    nu = min((F(1, 2)+rho-omega)*a0, s-omega/2)
    assert 0 < beta < min(A*(s+(F(1, 2)+rho)*a0)/2, A*(nu-d*share), pay)
    adjustment = {(1, 0): beta, (1, 1): A*d*share, (2, 0): E-2*pay,
                  (1, 2): -pay, (2, 1): -2*pay, (3, 0): -3*pay}
    marking = population(q, u, F(0), p, WEAK_WEIGHT, A, rho, type20_payment=SIGMA/2,
                         pair_factor=SIGMA, high_aware=high_aware, endpoint_adjustments=adjustment.items(),
                         radial_j_debit=X-A/5, v11=v11)
    assert F(marking["shared"]) == shared and F(marking["credit"]) == F(1981610721, 10**14)
    assert F(outer["threshold_abs_max_of_central_band"]) == EXTENT
    assert F(base["path"]["gamma_upper"]) < PATH_UPPER
    saving = min(F(marking["credit"]), F(outer["outer_graph_credit"]))
    assert saving > SAVING and PATH_UPPER-saving < GAMMA
    return {"scope": "one band through .0225 with the joint medium pool",
            "bilateral_medium": bilateral, "paired_star_cover": star_cover,
            "global_K0": str(k0), "raw_price_groups": {name: str(cap) for name, cap, *_ in groups},
            "price_guard_slacks": budget["coupling_guards"], "price_allowance": str(A),
            "marking": marking, "outer_credit": outer["outer_graph_credit"],
            "circuit_coefficient": str(1+1/(2*GAMMA))}


allocation = SimpleNamespace(final_allocation=final_allocation)


