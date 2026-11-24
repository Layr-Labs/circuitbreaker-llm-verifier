"""Directed two-arm union bounds and batched log-concave interpolation."""
from array import array
from bisect import bisect_right
from functools import lru_cache
from math import gcd, isqrt
from operator import mul
from struct import Struct


def arm_ranges(branches, side, parent_cells, scale):
    """Enclose the sum of individual log-bilinear lower interpolants on each cell."""
    result = []
    for i in range(parent_cells):
        for j in range(side - 1):
            first = i * side + j
            corners = (first, first + 1, first + side, first + side + 1)
            minimum = sum(min(branch[k] for k in corners) for branch in branches)
            maximum = max(sum(branch[k] for branch in branches) for k in corners)
            assert 0 <= minimum <= maximum <= scale
            result.append((minimum, maximum))
    return result


def either_arm(scale, y, odd_y, range_y, z, odd_z, range_z):
    """A cell-average floor after factoring out the common corner means."""
    lo_y, hi_y = range_y
    lo_z, hi_z = range_z
    ordered_y = scale * y + odd_y * z
    ordered_z = scale * z + odd_z * y
    affine_y = (scale - lo_z) * y + (scale - hi_y) * z + hi_y * lo_z
    affine_z = (scale - hi_z) * y + (scale - lo_y) * z + hi_z * lo_y
    return max(ordered_y, ordered_z, affine_y, affine_z) // scale



@lru_cache(maxsize=8)
def companion_piece_index(grid):
    widths = [(end-start)//4 for start,end in zip(grid,grid[1:])]
    assert all(width*4 == end-start for width,start,end in zip(widths,grid,grid[1:]))
    stride = gcd(*widths)
    pieces = []
    for j, width in enumerate(widths):
        for k in range(4):
            piece = j, k, grid[j]+k*width, width
            pieces.extend([piece] * (width//stride))
    return stride,pieces

def row_lookup(family, parent_index, grid, scale):
    """Evaluate directed lower tangents for all itineraries in one row."""
    side = len(grid)
    start = parent_index * side
    rows = [table[start:start+side] if start+side <= len(table) else (0,)*side for table in family]
    columns = list(zip(*rows))
    stride, pieces = companion_piece_index(tuple(grid))
    size = len(family)
    zeros = (0,) * size
    layout = Struct(f"<{size}Q")
    subdivisions = 4
    bits = 16
    unit = 1 << bits
    mask = int.from_bytes(layout.pack(*([(1 << (64-bits))-1] * size)), "little")
    assert max(map(max, columns)) <= scale and 3*scale*unit < 1 << 64

    def pack(values):
        return int.from_bytes(layout.pack(*values), "little")

    @lru_cache(maxsize=None)
    def node(j, k):
        if k == 0:
            return columns[j]
        if k == subdivisions:
            return columns[j+1]
        half = k & -k
        return _edge_means(node(j, k-half), node(j, k+half))

    @lru_cache(maxsize=None)
    def tangent(j, k):
        low, up, down = [0]*size, [0]*size, [0]*size
        for i, (a, b) in enumerate(zip(node(j, k), node(j, k+1))):
            if not a or not b:
                continue
            if a < b:
                low[i] = a
                up[i] = 2*a*(b-a)//(b+a)
            else:
                low[i] = b
                down[i] = 2*b*(a-b)//(a+b)
        return pack(low) << bits, pack(up), pack(down)

    def evaluate(lo, hi):
        j, k, start, width = pieces[lo//stride]
        offset = lo-start
        if lo == hi and offset == 0:
            return node(j, k)
        end = offset + hi-lo
        if end > width:
            # A log-concave itinerary stays above its smaller endpoint.
            return tuple(map(min, evaluate(lo, lo), evaluate(hi, hi)))
        base, increasing, decreasing = tangent(j, k)
        scaled = base + increasing*(unit*offset//width) + decreasing*(unit*(width-end)//width)
        values = (scaled >> bits) & mask
        if not values:
            return zeros
        return layout.unpack(values.to_bytes(layout.size, "little"))

    def lookup(lo, hi):
        return evaluate(lo, hi) if 0 <= lo <= hi < grid[-1] else zeros
    return lookup



def scalar_row_lookup(table, parent, grid):
    side = len(grid)
    values = table[parent*side:(parent+1)*side]
    stride, pieces = companion_piece_index(tuple(grid))
    unit = 1 << 16

    @lru_cache(maxsize=None)
    def node(j, k):
        if k == 0:
            return values[j]
        if k == 4:
            return values[j+1]
        half = k & -k
        return isqrt(node(j, k-half) * node(j, k+half))

    @lru_cache(maxsize=None)
    def tangent(j, k):
        a, b = node(j, k), node(j, k+1)
        if not a or not b:
            return 0, 0, False
        low, high = (a, b) if a < b else (b, a)
        return low, 2*low*(high-low)//(high+low), a < b

    @lru_cache(maxsize=4096)
    def evaluate(lo, hi):
        j, k, start, width = pieces[lo//stride]
        offset = lo-start
        if lo == hi and offset == 0:
            return node(j, k)
        end = offset + hi-lo
        if end > width:
            return min(evaluate(lo, lo), evaluate(hi, hi))
        low, slope, increasing = tangent(j, k)
        position = offset if increasing else width-end
        return low + (slope*(unit*position//width) >> 16)

    def lookup(lo, hi):
        return evaluate(lo, hi) if 0 <= lo <= hi < grid[-1] else 0

    return lookup

def _edge_means(left, right):
    assert len(left) == len(right)
    if left == right:
        return left
    return tuple(map(isqrt, map(mul, left, right)))


def square_cell_mean(plane, i, j):
    """Downward geometric mean of the four corners of a square cell."""
    return isqrt(isqrt(plane[i][j] * plane[i][j+1])
                 * isqrt(plane[i+1][j] * plane[i+1][j+1]))


def companion_mean(previous, current, i, j, cache):
    """Sum the corner floors for the disjoint companion guards."""
    corners = (
        previous[i][j], previous[i][j+1],
        previous[i+1][j], previous[i+1][j+1],
        current[i][j], current[i][j+1],
        current[i+1][j], current[i+1][j+1],
    )
    if None in corners:
        return 0
    key = tuple(map(id,corners))
    result = cache.get(key)
    if result is None:
        result = sum(
            isqrt(isqrt(isqrt(a*b)*isqrt(c*d)) * isqrt(isqrt(e*f)*isqrt(g*h)))
            for a, b, c, d, e, f, g, h in zip(*corners, strict=True)
            if a and b and c and d and e and f and g and h
        )
        cache[key] = result
    return result


def mixed_arm_floors(grid, root_limit, center_cells, near_tables, far_tables,
                     half_width, density, scale, coordinate_scale):
    """Floor the mixed negative-arm chords on the target rectangles."""
    panels = 4
    cside = len(grid) - 1
    mixed_c2 = [0] * (center_cells*cside)
    mixed_support = set()
    near_next = row_lookup(near_tables, 0, grid, scale)
    far_next = row_lookup(far_tables, 0, grid, scale)
    for j in range(root_limit - 1):
        v0, v1 = grid[j], grid[j+1]
        near_prev, far_prev = near_next, far_next
        near_next = row_lookup(near_tables, j+1, grid, scale)
        far_next = row_lookup(far_tables, j+1, grid, scale)
        near, far = (near_prev, near_next), (far_prev, far_next)
        previous_faces = None
        for i, x in enumerate(grid[:center_cells+1]):
            mu = half_width(x, v0), half_width(x, v1)
            if mu[1][0] <= 0:
                continue
            current_nodes = []
            for kt in range(panels+1):
                row_nodes = []
                for k, (mu_lo, mu_hi) in enumerate(mu):
                    a, b = panels+kt, panels-kt
                    hi_lo = (mu_lo*a)//panels
                    hi_hi = -(-(mu_hi*a)//panels)
                    lo_lo = (mu_lo*b)//panels
                    lo_hi = -(-(mu_hi*b)//panels)
                    weight = density(-(-(mu_hi*kt)//panels))*mu_lo
                    row_nodes.append((weight, near[k](hi_lo,hi_hi), far[k](lo_lo,lo_hi),
                                      far[k](hi_lo,hi_hi), near[k](lo_lo,lo_hi)))
                current_nodes.append(row_nodes)
            current_faces = []
            for lower, upper in zip(current_nodes, current_nodes[1:]):
                a,b = lower
                c,d = upper
                weight = isqrt(isqrt(a[0]*b[0])*isqrt(c[0]*d[0]))
                strands = [tuple(isqrt(isqrt(p*q)*isqrt(r*s)) if p and q and r and s else 0
                            for p,q,r,s in zip(a[h],b[h],c[h],d[h]))
                           for h in range(1,5)]
                current_faces.append((weight, *strands))
            if previous_faces is not None:
                total = 0
                for low, high in zip(previous_faces, current_faces):
                    weight = isqrt(low[0]*high[0])
                    if not weight:
                        continue
                    p = [sum(map(isqrt, map(mul, a,b))) for a,b in zip(low[1:],high[1:])]
                    total += weight*(p[0]*p[1]+p[2]*p[3])
                idx=(i-1)*cside+j
                mixed_c2[idx] = total // (scale*scale*coordinate_scale*panels)
                if mixed_c2[idx]:
                    mixed_support.update(((i-1,j),(i,j),(i-1,j+1),(i,j+1)))
            previous_faces = current_faces

    return mixed_c2, mixed_support


def ceil_div(a, b):
    return -(-a // b)


def positive_ray_floors(branch):
    """Directed floors on mixed and double positive-ray pair counts."""
    n, panels = 48, 16
    S, U = branch.SCALE, branch.UNIT
    RLO, RHI = branch.ROOT_LO, branch.ROOT_HI
    EXP_DEN = branch.EXP_DEN
    coord = n * U
    coord2 = coord ** 2
    knots = [i*U//2 for i in range(13)]
    knots.extend(i*U for i in range(7, 25))
    knots.extend(i*U for i in range(26, 45, 2))
    knots.extend(i*U for i in range(48, 3*n//2+1, 4))
    assert knots[-1] == 3*coord//2
    side = len(knots)
    cside = side-1
    max_h = 5*coord//4
    m_lo = [RLO*y//U for y in knots]
    m_hi = [ceil_div(RHI*y,U) for y in knots]

    @lru_cache(maxsize=32768)
    def density(h):
        q=ceil_div(6*h*h*EXP_DEN,coord2)
        return branch.CG_NUM*branch.exp_floor(q)//S

    half_mesh = 1024
    prefix = [0]
    accumulated = 0
    for i in range(1,2*half_mesh+1):
        q=ceil_div(2*EXP_DEN*(3*i*i-3*i+1),half_mesh**2)
        accumulated += branch.CG_NUM*branch.exp_floor(q)
        prefix.append(accumulated//(S*half_mesh))

    def half_mass(b):
        if b <= 0:
            return 0
        if b >= 2*coord:
            return prefix[-1]
        j,rem=divmod(b*half_mesh,coord)
        return prefix[j]+((prefix[j+1]-prefix[j])*rem//coord)

    U0 = array('I',(half_mass(l-x//2) for x in knots for l in m_lo))
    B0 = array('I',(half_mass(l+x//2) for x in knots for l in m_lo))

    def table(upper, lower, geometry):
        result = array('I',[0])*(side*side)
        for j in range(side):
            hi_lookup = None if upper is None else scalar_row_lookup(upper,j,knots)
            lo_lookup = None if lower is None else scalar_row_lookup(lower,j,knots)
            for i,x in enumerate(knots):
                plus = x//2 if geometry in ('opposite_parent','alternating') else -x//2
                lo,hi = m_lo[j]+plus,m_hi[j]+plus
                if geometry in ('both_same','opposite_parent'):
                    first,last=0,min(lo,max_h)
                elif geometry in ('one_opposite','alternating'):
                    first,last=max(abs(lo),abs(hi)),max_h
                else:
                    first,last=0,min(-hi,max_h)
                if last <= first:
                    continue
                previous_h=first
                previous_value=0
                total=0
                for k in range(panels+1):
                    h=first+(last-first)*k//panels
                    if geometry=='two_opposite':
                        a,b=-hi+h,-lo+h
                        c,d=-hi-h,-lo-h
                    else:
                        a,b=lo+h,hi+h
                        if geometry in ('one_opposite','alternating'):
                            c,d=h-hi,h-lo
                        else:
                            c,d=lo-h,hi-h
                    large=S if hi_lookup is None else hi_lookup(a,b)
                    small=S if lo_lookup is None else lo_lookup(c,d)
                    value=density(h)*large*small
                    total+=(h-previous_h)*isqrt(previous_value*value)
                    previous_h,previous_value=h,value
                result[i*side+j]=total//(coord*S*S)
        return result

    G0 = table(U0,None,'both_same')
    safe = table(G0,G0,'opposite_parent')
    b1 = table(None,B0,'alternating')
    b2 = table(None,b1,'alternating')
    outside_even = table(None,b1,'one_opposite')
    tails = [table(f,safe,'one_opposite') for f in (U0,outside_even)]
    ray = [table(f,safe,'one_opposite') for f in tails]
    stable = [table(a,b,'two_opposite') for a in (B0,b2) for b in (B0,b2)]

    def cell_floors(values, x_count, y_count):
        return [[isqrt(isqrt(values[i*side+j]*values[i*side+j+1])*
                       isqrt(values[(i+1)*side+j]*values[(i+1)*side+j+1]))
                 for j in range(y_count)] for i in range(x_count)]

    ix = len(knots)-1
    iy = len(knots)-1
    def combined(bank):
        floors = [cell_floors(table,ix,iy) for table in bank]
        return [[sum(table[i][j] for table in floors) for j in range(iy)] for i in range(ix)]
    rays = combined(ray)
    stable = combined(stable)
    joint = [0,0]
    previous = None
    for i,x in enumerate(knots[:ix+1]):
        x_lo = 2*RLO*x//U
        current = [[0] * (iy + 1) for _ in range(iy + 1)]
        x_term = 33*x*x
        for j,y in enumerate(knots[:iy+1]):
            if not y:
                continue
            y_term = x_term + 12*y*y
            for k in range(j, iy+1):
                z = knots[k]
                magnitude = y + z
                if x_lo <= magnitude:
                    break
                exponent = y_term + 12*z*magnitude - (24*RLO*x*magnitude)//U
                q = ceil_div(exponent*EXP_DEN, 2*coord2)
                weight = branch.C3_NUM*branch.exp_floor(q)//S
                current[j][k] = current[k][j] = weight
        if previous is not None:
            dx = (knots[i]-knots[i-1]) // (U//2)
            ray_row, stable_row = rays[i-1], stable[i-1]
            for j in range(iy):
                dy = (knots[j+1]-knots[j]) // (U//2)
                ray_y, stable_y = ray_row[j], stable_row[j]
                if not (ray_y or stable_y):
                    continue
                for k in range(j, iy):
                    ray_z, stable_z = ray_row[k], stable_row[k]
                    if not (ray_y*ray_z or ray_y*stable_z or ray_z*stable_y):
                        continue
                    w = isqrt(square_cell_mean(previous,j,k)*square_cell_mean(current,j,k))
                    if not w:
                        continue
                    dz = (knots[k+1]-knots[k]) // (U//2)
                    volume = 3*dx*dy*dz
                    pair_multiplicity = 1 if j == k else 2
                    joint[0] += pair_multiplicity*volume*(w*ray_y*ray_z//S)
                    combined = (w*ray_y*stable_z)//S
                    if j != k:
                        combined += (w*ray_z*stable_y)//S
                    joint[1] += 2*volume*combined
        previous=current
    return joint, (2*n)**3*S*S




SCALE = 100_000_000       # 1e8
UNIT = 100_000_000        # 1e8 represents 1/n
ROOT_LO = 141_421_356     # < sqrt(2) * 1e8
ROOT_HI = 141_421_357     # > sqrt(2) * 1e8
INV_R2_LO = 70_710_678    # < (1/sqrt(2)) * 1e8
CG_NUM = 276_395_319      # < 2*sqrt(6/pi) * 1e8
C2_NUM = 47_746_482       # < 3/(2*pi) * 1e8
C3_NUM = 65_984_522       # < sqrt(108)/(2*pi)^(3/2) * 1e8
HALF_CDF_MESH = 2048
EXP_DEN = 1_000_000
EXP_CUT = 10_500_000
WEIGHT_CUT = 200_000


def curvature_floor(value, log_gain_num, log_gain_den):
    """Scale a floor by 1 + delta, a lower bound on exp(delta) for delta >= 0."""
    return value + value*log_gain_num//log_gain_den


def lmlo(a, b, _isqrt=isqrt):
    if not a or not b:
        return 0
    if a == b:
        return a
    m = _isqrt(a * b)
    return (_isqrt(a * m) + _isqrt(b * m)) >> 1


def geo_interp(a, b, num, den, _isqrt=isqrt):
    if a == b or num <= 0:
        return a
    if num >= den:
        return b
    if not a or not b:
        return 0
    lo, hi = a, b
    mn = a if a < b else b
    diff = (b - a) if b > a else (a - b)
    if diff * 48 > mn:
        num2 = num << 1
        if num2 >= den:
            lo = _isqrt(lo * hi)
            num = num2 - den
        else:
            hi = _isqrt(lo * hi)
            num = num2
        if diff * 8 > mn:
            num2 = num << 1
            if num2 >= den:
                lo = _isqrt(lo * hi)
                num = num2 - den
            else:
                hi = _isqrt(lo * hi)
                num = num2
        if not num:
            return lo
        if not lo or not hi:
            return 0
    return (den * lo * hi) // ((den - num) * hi + num * lo)


@lru_cache(maxsize=None)
def exp_floor(i):
    if i <= 0:
        return SCALE
    if i > EXP_CUT:
        return 0
    hi, lo = divmod(i, 4096)
    return (_small_exp_floor(hi*4096) * _small_exp_floor(lo)) // SCALE


@lru_cache(maxsize=None)
def _small_exp_floor(i):
    if i <= 0:
        return SCALE
    if i > EXP_CUT:
        return 0
    h = 256 * EXP_DEN
    u = -(-(i * SCALE) // h)
    l = (i * SCALE) // h
    v = SCALE - u + (l * i) // (2 * h) - (-(-(-(-(u * i) // (2 * h)) * i) // (3 * h)))
    if v <= 0:
        return 0
    for _ in range(8):
        v = (v * v) // SCALE
    return v


def build_graded_grid(n):
    assert n in (96, 104, 112, 320)
    U = UNIT
    pts = [0, U // 4, U // 2, U, (3 * U) // 2] if n == 104 else [0, U // 8, U // 4, U // 2, U, (3 * U) // 2]
    if n == 96:
        limits = (24, 48, 96, 144)
        steps = (1, 2, 4, 8)
    elif n == 112:
        limits = (28, 56, 112, 168)
        steps = (1, 2, 4, 8)
    elif n == 104:
        limits = (24, 48, 100, 156)
        steps = (1, 2, 4, 8)
    else:
        limits = (32, 92, 156, 220, 316, 480)
        steps = (1, 2, 4, 8, 16, 32)
    start = 1
    for end, step in zip(limits, steps, strict=True):
        pts.extend(k*U for k in range(start + step, end + 1, step))
        if pts[-1] < end*U:
            pts.append(end*U)
        start = end
    assert pts[-1] == 3 * n * U // 2
    return pts


def cubic_density_lookup(v_div):
    """Cache exact outward-rounded half-normal lattice densities."""
    # A coordinate contributes only its outward-rounded density index.
    # The final retained index is the largest t with 3*t*t <= 800*EXP_CUT.
    density_index_max = isqrt(800 * EXP_CUT // 3)
    density_table = [(CG_NUM * exp_floor((3*t*t + 799) // 800)) // UNIT
                     for t in range(density_index_max + 1)]
    density_support = density_index_max * v_div

    def density(v):
        magnitude = abs(v)
        if magnitude > density_support:
            return 0
        return density_table[(magnitude + v_div - 1) // v_div]

    return density


def packed_scaling(lane_count):
    layout = Struct('<' + 'I12x' * lane_count)
    word_mask = int.from_bytes(layout.pack(*([2**32 - 1] * lane_count)), 'little')
    reciprocal = ceil_div(1 << 96, SCALE)
    assert 3 * SCALE**3 < (1 << 96)

    def pack_row(values):
        return int.from_bytes(layout.pack(*values), 'little')

    return pack_row, layout, word_mask, reciprocal


def degree_two_weight_planes(grid, x_limit, y_limit, divisor, root_upper, prefactor):
    """Use symmetry and repeated horizontal means in the root density."""
    quantized = [ceil_div(point, divisor) for point in grid]
    pairs = [[(y+z, 12*(y*y+y*z+z*z)) for z in quantized[:y_limit+1]]
             for y in quantized[:y_limit+1]]
    denominator = 3200*UNIT
    for x in quantized[:x_limit+1]:
        base = 33*x*x
        cross = 24*root_upper*x
        plane = [[0]*(y_limit+1) for _ in range(y_limit+1)]
        for j in range(y_limit+1):
            for k in range(j, y_limit+1):
                linear, quadratic = pairs[j][k]
                exponent = -(-(UNIT*(base+quadratic)+cross*linear)//denominator)
                weight = prefactor*exp_floor(exponent)//UNIT
                plane[j][k] = plane[k][j] = weight
        horizontal = [[isqrt(row[k]*row[k+1]) for k in range(y_limit)] for row in plane]
        yield [[0]*j + [isqrt(horizontal[j][k]*horizontal[j+1][k])
                       for k in range(j,y_limit)] for j in range(y_limit)]

# Row-major claims, checked entry by entry against the negative recurrence.
ENDPOINT_CELL_FLOORS = [
13811474,14168391,14529473,15253422,16342891,17443622,19005812,21350249,23748781,26196514,28667791,31149997,33652003,36160773,38646195,41106176,43544416,45922777,48265261,50567939,52804750,54972113,57070527,59114877,61074316,62957398,64782967,66524150,68195414,69801785,71327638,72794018,74765022,77310001,79636594,81739549,83635405,85366951,86937457,88345994,89621842,90785580,91827447,92760731,93609194,94369562,95236614,96334081,97202575,97878938,98401744,98802903,99104000,99325649,99489239,99609431,99694899,99754668,99792937,99822910,99752189,99781284,99792978,99877419,0,0,0,
14015653,14284821,14647659,15374845,16465611,17569058,19136566,21485324,23887546,26337964,28811463,31295145,33797660,36305978,38790705,41249617,43685978,46062547,48402142,50701357,52934476,55098074,57193367,59232909,61188030,63068103,64889526,66626583,68294330,69896511,71417811,72880822,74846512,77383910,79703782,81801351,83691668,85416759,86982957,88388219,89660097,90820355,91859002,92789209,93634714,94392927,95256569,96349671,97214836,97888533,98408854,98807925,99108057,99329009,99491323,99610827,99696628,99756233,99793695,99823341,99752321,99781020,99792280,99876919,0,0,0,
14187766,14453522,14820924,15554241,16648189,17756199,19332147,21687517,24095144,26549700,29026370,31512180,34015583,36523131,39006893,41464209,43897658,46271359,48606448,50900712,53128279,55286171,57376843,59409289,61358139,63233801,65048037,66778415,68441204,70037270,71552249,73010015,74967311,77493719,79803919,81892764,83775049,85491372,87051470,88451037,89716920,90872257,91905934,92831522,93672287,94426662,95285452,96372495,97233082,97902899,98419392,98814494,99112679,99332978,99494563,99613011,99698124,99757191,99794521,99823879,99752756,99780384,99791579,99877154,0,0,0,
14521551,14787812,15163916,15910610,17012426,18130032,19722515,22091277,24509963,26972880,29455611,31945284,34450025,36955991,39437733,41891667,44304257,46685173,49011944,51297326,53512835,55659183,57741658,59760319,61695730,63562367,65363187,67080177,68732688,70316552,71819388,73265885,75207302,77713036,80002701,82074236,83940319,85638983,87187382,88576214,89828723,90972741,91997518,92914184,93746228,94493623,95341395,96416275,97268326,97930122,98440444,98829513,99123682,99341662,99501463,99618095,99700579,99758712,99796950,99825337,99752980,99780492,99792023,99877897,0,0,0,
14994718,15259709,15642170,16402124,17512802,18640485,20252012,22638457,25070754,27543799,30033512,32527459,35033376,37536678,40014259,42462672,44866829,47237089,49552633,51806265,54009560,56156977,58228102,60228850,62145638,63999280,65782611,67481966,69119811,70686313,72173562,73605647,75526095,78004396,80266732,82315988,84160294,85835968,87367210,88742515,89978356,91105937,92118875,93023863,93844163,94582291,95416133,96475985,97315616,97966514,98468328,98849633,99139299,99353473,99509660,99624065,99704845,99761632,99799086,99826731,99754185,99781041,99793016,99879093,0,0,0,
15468901,15735909,16124833,16897883,18017961,19156996,20787436,23189419,25634597,28116771,30612383,33109592,35616420,38117122,40589234,43031284,45427563,47787045,50090882,52332594,54521047,56652488,58712435,60695222,62592146,64431832,66197600,67880313,69502591,71051106,72523758,73941349,75840751,78292412,80527720,82554565,84377959,86032193,87543887,88905236,90127269,91239193,92239718,93133959,93942228,94668576,95491355,96537039,97362723,98002585,98495315,98869288,99153935,99364741,99517778,99629348,99708589,99764290,99800221,99827655,99755948,99781919,99793517,99879381,0,0,0,
16156369,16430123,16830805,17625523,18765409,19922243,21581037,24004357,26465630,28961366,31463007,33971671,36478711,38975311,41421538,43850507,46255677,48598719,50884416,53107290,55273514,57381794,59396015,61353098,63249147,65065797,66804827,68464399,70063237,71586429,73038026,74432788,76301896,78714480,80910295,82902951,84697536,86322443,87804044,89141761,90342449,91432002,92414100,93292227,94084022,94794669,95600396,96623778,97430008,98054150,98533604,98898579,99175642,99381226,99529979,99637636,99714432,99767891,99802561,99829039,99757333,99782715,99794296,0,0,0,0,
17144769,17423675,17834587,18650685,19810376,20985373,22674682,25126794,27608298,30117444,32626950,35137927,37641237,40129367,42565198,44978026,47362917,49682932,51943877,54138927,56274694,58352710,60333968,62255479,64122401,65878777,67581724,69238375,70774309,72263391,73717867,75082214,76910408,79270252,81415391,83362303,85121750,86709938,88148752,89453942,90626579,91687924,92644810,93499844,94270791,94962441,95743376,96736426,97517506,98122407,98585573,98937546,99204444,99402717,99545567,99649128,99723114,99772674,99805982,99831532,99759421,99783628,99800893,0,0,0,0,
18155066,18439556,18860117,19695777,20873590,22067439,23785139,26260498,28759023,31278661,33792961,36304173,38801309,41277820,43701846,46096417,48458304,50755129,52990725,55155646,57261164,59312686,61260415,63143801,64980262,66704449,68372869,69996752,71502543,72960430,74383277,75684561,77470496,79812804,81908500,83809702,85532880,87085873,88485521,89757760,90903386,91937737,92868983,93702376,94451988,95122686,95880782,96845699,97602859,98189245,98637468,98975872,99232564,99423016,99560326,99659780,99730695,99777603,99809634,99833810,99761600,99784340,99810529,0,0,0,0,
19187441,19477374,19906798,20759246,21955671,23165686,24908493,27403920,29916357,32443036,34959271,37453416,39939546,42409406,44830762,47206108,49542667,51790507,53995892,56157321,58233857,60244278,62167835,64018740,65805378,67513443,69149473,70739665,72214294,73642754,75034405,76308144,78053963,80343152,82390443,84246475,85927061,87445469,88815039,90054917,91172762,92179658,93086169,93899278,94627849,95279131,96015182,96952252,97686512,98254159,98687121,99013187,99259880,99443016,99574728,99669702,99737439,99782663,99813692,99835025,99762539,99785053,99819675,0,0,0,0,
20240837,20535806,20973196,21841392,23054029,24278409,26043004,28554889,31078501,33602802,36125293,38614294,41090694,43546107,45949932,48304450,50615653,52840174,55016172,57144373,59193458,61166451,63054036,64878847,66632355,68305428,69911192,71448455,72912082,74311083,75671856,76918566,78624791,80861289,82860626,84672122,86310549,87794595,89135599,90344146,91434219,92415422,93297869,94089886,94799066,95433030,96146577,97055752,97766789,98316060,98734443,99049399,99286442,99462216,99588523,99679573,99744146,99786306,99817192,99837373,99764236,99785753,99828170,0,0,0,0,
21313644,21613187,22057705,22939912,24167098,25404167,27187377,29708649,32236303,34770603,37290277,39771983,42235916,44675361,47058376,49390320,51677199,53875664,56020776,58117151,60138545,62075600,63927213,65723950,67443286,69081607,70657582,72163104,73595269,74965499,76273044,77515522,79183928,81369184,83320967,85089360,86689603,88138255,89448477,90626751,91689303,92645760,93504607,94275439,94965010,95580808,96272790,97154958,97843723,98375447,98781779,99085785,99311663,99479901,99601765,99689663,99751133,99790239,99820273,99840136,99767497,99787253,99836010,0,0,0,0,
22404143,22707806,23158671,24053221,25293206,26539401,28336022,30872095,33404536,35937888,38443342,40924321,43374091,45796146,48157696,50464825,52708826,54895464,57009673,59076015,61063717,62965397,64786526,66553133,68237805,69842439,71389285,72862174,74263314,75606773,76884546,78099134,79731298,81841380,83772114,85499994,87063678,88474579,89752676,90902861,91938485,92870043,93706093,94456646,95127602,95727998,96397696,97251021,97917767,98433050,98829041,99121107,99335766,99497344,99614629,99698078,99757109,99794574,99822409,99841370,99769297,99788136,0,0,0,0,0,
23510618,23817900,24274312,25179658,26430694,27685714,29495372,32039387,34573560,37103665,39600886,42069819,44504207,46905429,49245337,51527625,53744754,55901064,57984526,60021837,61971558,63836546,65631482,67365397,69014887,70588450,72106425,73547166,74917772,76235108,77483379,78669958,80265968,82327727,84212629,85905278,87434068,88804781,90050265,91171921,92180431,93088227,93903342,94633632,95285844,95868724,96516713,97344173,97989608,98487708,98869444,99151902,99359330,99514109,99626714,99706300,99762637,99798586,99825442,99843200,99770238,99788526,0,0,0,0,0,
24631242,24941657,25402880,26317571,27577664,28839665,30659753,33208434,35742018,38267783,40754792,43209067,45626185,48003829,50297901,52551855,54766780,56891796,58945234,60930896,62844578,64693604,66461746,68160499,69775039,71321055,72807861,74215998,75558672,76849182,78068695,79228721,80788229,82802356,84643191,86279700,87767858,89127601,90341423,91435346,92418486,93303290,94095935,94805402,95438220,96003606,96631643,97433982,98058836,98540015,98906551,99180673,99382162,99530384,99638471,99715333,99769105,99802926,99828971,99845623,99771765,99789079,0,0,0,0,0,
25764233,26077310,26542623,27465051,28729676,29999823,31827746,34378011,36908966,39426972,41901518,44340056,46737968,49092222,51363019,53590270,55772916,57867268,59891777,61844179,63722605,65536895,67275985,68939546,70520641,72039293,73493284,74869393,76186572,77447866,78638854,79775830,81299424,83267534,85063709,86657566,88111823,89442178,90603801,91669745,92649995,93511285,94281665,94971997,95588044,96137393,96745150,97522346,98127080,98591550,98944162,99209036,99403650,99546203,99649653,99723025,99774264,99806499,99831201,99847402,99773886,99789452,0,0,0,0,0,
26907817,27223109,27691782,28620313,29890414,31165030,32998416,35547344,38073627,40580011,43039959,45461973,47838414,50168910,52415940,54614223,56764071,58828591,60824650,62740212,64581939,66366469,68073754,69702273,71251722,72741524,74163546,75509865,76802471,78034264,79196582,80311551,81800361,83722070,85462361,87028285,88448858,89749055,90882443,91921822,92875292,93712645,94463325,95134064,95732290,96264918,96853747,97608639,98193637,98641103,98979441,99235323,99423856,99560782,99660067,99729982,99778993,99810204,99833465,99848699,99775348,99790172,0,0,0,0,0,
28060403,28377330,28848618,29781936,31056051,32333285,34169908,36715399,39234896,41718451,44171656,46575005,48928378,51219765,53455482,55623399,57740620,59775519,61742344,63622056,65427962,67180691,68829460,70421644,71968063,73428103,74818115,76136226,77402935,78606489,79742189,80835287,82290571,84165009,85866611,87395789,88780025,90032745,91155177,92167302,93094760,93893029,94624410,95291704,95872676,96390603,96960052,97690976,98256287,98688086,99015151,99262157,99443960,99575347,99670258,99737233,99783486,99813864,99835988,99849520,99775566,99790394,0,0,0,0,0,
29220081,29538305,30011540,30948403,32225063,33503160,35340844,37881160,40390677,42859574,45294621,47658137,49983308,52273117,54480116,56617180,58702670,60708063,62640264,64484209,66259384,67979781,69596870,71155672,72670634,74100191,75458004,76749020,77987680,79163334,80275761,81345249,82767149,84598005,86264107,87758879,89104720,90326491,91421209,92406814,93308812,94086815,94797875,95445538,96013070,96515613,97062933,97768850,98316093,98734208,99050077,99287844,99462951,99589398,99680355,99744512,99787962,99816834,99838089,99851102,99776851,0,0,0,0,0,0,
30385222,30704480,31179079,32118310,33396057,34673854,36510667,39043468,41539990,43992262,46408037,48750100,51051415,53313554,55489892,57595743,59650124,61627219,63522350,65329035,67075778,68763812,70348614,71874928,73358359,74757857,76085143,77350984,78561212,79708234,80798727,81842655,83231405,85021022,86651948,88113762,89421456,90612332,91681209,92640425,93503002,94275413,94966750,95595146,96145608,96632341,97162698,97847121,98376541,98780424,99084075,99312252,99481243,99602771,99689941,99751327,99792167,99819984,99840282,99853075,99778046,0,0,0,0,0,0,
31554549,31874341,32349702,33290176,34567695,35843850,37678037,40194492,42669720,45116819,47511794,49831050,52107035,54339872,56485169,58559745,60582949,62500459,64356050,66160623,67877527,69521979,71084170,72579602,74030984,75372547,76669359,77938202,79121815,80241409,81309133,82328959,83685076,85433952,87025705,88454285,89729724,90890951,91934290,92867193,93706908,94458773,95130449,95740309,96272898,96745800,97259869,97923942,98435837,98824185,99116762,99336337,99498764,99615015,99698290,99757136,99796095,99823141,99841951,99853661,99779133,0,0,0,0,0,0,
32726458,33046491,33522044,34462611,35738632,37011994,38841787,41345776,43804358,46232290,48593623,50900204,53148888,55350556,57465310,59509305,61501457,63384351,65204241,66977313,68664492,70275404,71803821,73270639,74690797,76003006,77271608,78510547,79668163,80762920,81807019,82803224,84128670,85835922,87384634,88781747,90030369,91162015,92179853,93087846,93905402,94637511,95290180,95871236,96397904,96856896,97354249,97996732,98492032,98866703,99148813,99359804,99514933,99626625,99706726,99763053,99799817,99825906,99843518,99854499,99780920,0,0,0,0,0,0,
33899563,34219441,34694689,35634393,36907658,38177145,40000898,42489947,44930022,47337840,49675816,51956809,54176829,56346347,58430509,60444159,62406695,64254153,66036317,67778668,69436433,71013729,72508845,73947224,75337321,76620732,77861344,79071226,80172475,81241870,82293019,83266015,84561417,86228306,87736677,89102194,90324427,91426753,92419206,93303459,94099255,94811949,95444674,96008359,96518400,96955160,97436819,98066943,98545934,98909783,99181466,99383047,99531197,99638051,99715761,99769372,99803182,99828513,99845181,99855906,99782084,0,0,0,0,0,0,
35072576,35391884,35866323,36804248,38073640,39338219,41154047,43626153,46046514,48432840,50745990,52999760,55190298,57327727,59381141,61364052,63273940,65106899,66854927,68565679,70191506,71734768,73198898,74609068,75952696,77224466,78439302,79620877,80697464,81743545,82768767,83718865,84984104,86613968,88085648,89416730,90610752,91683721,92650757,93512915,94287233,94980858,95595007,96143493,96638670,97061001,97525774,98134683,98597257,98951547,99213351,99404948,99547123,99649550,99723073,99774433,99806657,99830918,99846787,99856682,99782845,0,0,0,0,0,0,
36244285,36562671,37035820,37971043,39235504,40494050,42299826,44753221,47153198,49516215,51803480,54028379,56188293,58294307,60317010,62269417,64142686,65939648,67658649,69337800,70930068,72439331,73874983,75258172,76570854,77813489,79004130,80158555,81210552,82233153,83233692,84161202,85397072,86993419,88430421,89724042,90888122,91933764,92877993,93717928,94468880,95145028,95742026,96273283,96747190,97164098,97612656,98200758,98647518,98988063,99240727,99425225,99561599,99660717,99730300,99779104,99809974,99833552,99848325,99857629,99783880,0,0,0,0,0,0,
37413260,37730514,38201924,39133656,40392189,41643863,43437851,45870912,48249532,50587201,52847526,55043091,57171529,59245902,61238004,63162150,64999401,66758142,68447292,70094696,71651954,73128386,74538164,75893976,77176389,78390279,79556225,80664321,81712385,82711190,83687655,84564267,85768207,87360537,88764115,90023565,91158583,92178263,93098625,93917316,94646509,95303913,95883349,96397964,96858357,97262933,97695378,98263583,98695408,99022567,99265612,99444371,99575774,99670752,99738142,99784257,99813512,99835885,99849892,99858321,99784737,0,0,0,0,0,0,
38578309,38894262,39363589,40291023,41542660,42782189,44557893,46978785,49334713,51645784,53877371,56043065,58140455,60182890,62143934,64035041,65836772,67562259,69220825,70835758,72358148,73802979,75186885,76514325,77767808,78955310,80096827,81179409,82202969,83179359,84132483,84989998,86163257,87713577,89084483,90315574,91422870,92416543,93310323,94108248,94818601,95457394,96021412,96521728,96968167,97358892,97775509,98324599,98740907,99057101,99292142,99464019,99589779,99680083,99744781,99789365,99817175,99838048,99851225,99858752,99785396,0,0,0,0,0,0,
39738674,40052954,40519816,41442182,42686071,43916167,45677377,48075966,50407548,52690927,54892706,57027975,59094515,61104845,63035227,64888483,66653924,68351246,69979627,71561256,73050453,74465249,75822163,77119678,78344554,79508335,80625878,81682621,82682083,83637336,84567253,85404688,86548536,88057798,89396325,90600052,91679416,92647701,93517064,94293914,94985285,95606179,96154805,96638488,97072406,97452671,97854042,98384824,98786061,99089018,99316367,99483009,99602566,99689370,99751055,99793441,99820230,99840170,99852456,99859706,99786006,0,0,0,0,0,0,
40893257,41205628,41669568,42586430,43821346,45040902,46786860,49161857,51467549,53722462,55894407,57998193,60033222,62011597,63884173,65696455,67458154,69125549,70724138,72247747,73702505,75113667,76444279,77712002,78908560,80049039,81142731,82173919,83150104,84085373,84971573,85809314,86926093,88396396,89701429,90877327,91929502,92872355,93718481,94475000,95148008,95750126,96282903,96752541,97173247,97542650,97929518,98438672,98826465,99119508,99338824,99500085,99614876,99698422,99757388,99797609,99823183,99842228,99853647,99860242,99786605,0,0,0,0,0,0,
42040947,42351234,42811990,43722313,44946925,46155602,47885746,50236241,52514984,54739676,56880896,58953478,60957094,62903164,64743019,66520149,68247905,69884937,71453190,72944577,74367668,75747127,77051725,78290741,79460504,80577453,81646255,82652520,83607747,84523723,85387677,86205093,87297325,88730757,90001632,91147229,92173101,93091512,93914847,94650875,95305902,95890470,96410489,96867321,97272024,97629251,98003741,98495129,98867913,99150497,99361802,99516379,99626438,99706735,99763971,99802277,99826212,99843980,99854859,99860740,99787694,0,0,0,0,0,0,
43180985,43489076,43946321,44848907,46062517,47260119,48973309,51297634,53548864,55742365,57852048,59893846,61866024,63780409,65584207,67325418,69022653,70629710,72166978,73626109,75018884,76368274,77645538,78856513,80000798,81094411,82138298,83120094,84055673,84951535,85793747,86591841,87663189,89062017,90296149,91410599,92409866,93303380,94104709,94822172,95459050,96026163,96531551,96974949,97366555,97713253,98074633,98549125,98908971,99181456,99383802,99531751,99637895,99714692,99769498,99805660,99828653,99845592,99855908,99861303,99788520,0,0,0,0,0,0,
44313089,44618044,45071290,45965915,47168052,48353691,50049024,52345974,54568631,56715945,58808861,60819179,62759505,64624768,66411784,68117071,69782729,71360066,72866316,74293201,75656250,76976931,78225980,79408560,80528458,81599222,82618288,83577307,84493422,85368416,86189781,86969742,88016689,89358607,90581334,91667096,92641365,93512201,94291804,94988789,95607629,96157277,96647106,97077893,97457325,97794395,98142311,98600935,98949478,99210851,99403352,99546504,99649385,99722734,99774279,99808476,99831479,99847399,99856789,99861822,99789188,0,0,0,0,0,0,
45953058,46253801,46701036,47583712,48768296,49935458,51603752,53860503,56039483,58140465,60186131,62150309,64045026,65859400,67594081,69253296,70872159,72405674,73865401,75245840,76568057,77845282,79055440,80197419,81281429,82318206,83301249,84229745,85117794,85960845,86752476,87506900,88510639,89805026,90984916,92029826,92969108,93807022,94554909,95222308,95815510,96341101,96809283,97222439,97585078,97905931,98235348,98672106,99004081,99250334,99432079,99566545,99663636,99732081,99780505,99811523,99833308,99847957,99856714,99861302,99788944,0,0,0,0,0,0,
48156028,48451056,48889094,49753183,50911028,52050619,53678752,55876347,57992841,60029822,62010241,63910588,65739230,67482022,69149425,70748336,72304933,73778131,75174737,76496025,77764928,78985528,80116737,81204764,82268052,83258798,84195763,85084036,85932182,86733975,87488493,88212276,89165706,90394669,91516000,92506834,93398532,94191105,94897491,95528041,96088038,96585662,97024099,97410640,97752308,98051843,98356811,98765361,99074861,99304356,99471453,99594293,99682331,99744813,99789836,99817129,99837246,99850462,99858177,99861891,99789601,0,0,0,0,0,0,
50312794,50601118,51029026,51872734,53001918,54112384,55698054,57832629,59884961,61857969,63772873,65582908,67335290,69041016,70646665,72184641,73679944,75091343,76426592,77693458,78911032,80067164,81157158,82196526,83210557,84156163,85049932,85897815,86705474,87469132,88190802,88884211,89792049,90958318,92010420,92948588,93787491,94554401,95222164,95817108,96344406,96814231,97226674,97587676,97907163,98187770,98469973,98851209,99139406,99352574,99507276,99619912,99700557,99757625,99799218,99823853,99841569,99852702,99859086,99862147,99790350,0,0,0,0,0,0,
52418503,52699689,53116915,53938851,55037749,56117740,57658803,59728081,61715285,63624049,65475395,67219727,68900823,70540386,72085725,73564648,74997956,76315406,77588092,78839658,80008216,81112061,82151080,83145253,84111481,84986023,85835697,86673600,87443278,88170106,88862044,89520155,90380167,91492351,92491933,93382064,94175204,94898664,95529545,96089584,96585029,97025404,97413253,97753363,98052740,98314893,98575412,98930075,99197023,99395378,99539500,99643258,99717584,99769985,99806905,99829184,99845618,99855354,99860213,99862659,99791448,0,0,0,0,0,0,
54469698,54743263,55149163,55948439,57015890,58064179,59559167,61562077,63483463,65328687,67098122,68788908,70404228,71981264,73467200,74887714,76236598,77519618,78738354,79937274,81055482,82107393,83099865,84052356,84971941,85806817,86616471,87412930,88145999,88841220,89500345,90121286,90939774,91998384,92948737,93793850,94543534,95224377,95819760,96347611,96812369,97225360,97589734,97909149,98189714,98433848,98674299,99004206,99250687,99434572,99569125,99664821,99732865,99780046,99812231,99832941,99848678,99857324,99861640,99863370,99799548,0,0,0,0,0,0,
56463658,56729353,57123390,57899103,58934180,59949986,61398278,63334403,65189940,66971148,68671275,70294551,71849598,73364460,74792311,76154924,77444912,78671678,79839868,80986913,82053241,83056716,84007260,84919122,85780061,86589876,87361687,88117998,88816539,89479312,90103005,90693185,91474275,92478538,93380172,94181026,94889664,95521706,96082652,96579108,97028565,97416308,97757136,98055368,98317483,98545564,98767359,99074424,99302935,99472261,99596535,99684299,99745929,99789208,99818109,99837269,99850941,99858466,99861970,99863070,99815547,0,0,0,0,0,0,
58398676,58655860,59037713,59789211,60791263,61774306,63175349,65044175,66835503,68550638,70182855,71742860,73236872,74691392,76062832,77368399,78600495,79773840,80893525,81988480,83005039,83963735,84873915,85746210,86565824,87336284,88073305,88791180,89432760,90057923,90674840,91238880,91979741,92933242,93786765,94545060,95215445,95813874,96342268,96808247,97231239,97594967,97913256,98192414,98436673,98649778,98854663,99141118,99352889,99507589,99621276,99701907,99758580,99798000,99822753,99840838,99853214,99859700,99862258,99862786,99832470,0,0,0,0,0,0,
60272220,60521124,60890635,61617675,62586284,63536328,64890110,66692615,68417975,70065645,71633649,73133295,74567866,75964180,77280539,78528703,79705581,80829581,81900748,82944126,83913661,84829788,85701544,86533393,87311579,88047227,88751049,89436434,90046592,90636771,91222596,91757711,92459445,93364147,94170942,94890230,95525121,96089590,96587294,97024420,97412029,97752408,98051156,98320041,98547414,98746034,98934352,99201358,99398054,99540098,99644097,99718193,99770610,99806503,99827701,99844349,99854348,99859791,99862709,99862608,99846535,0,0,0,0,0,0,
62083423,62324327,62681444,63383902,64318949,65236070,66542855,68267188,69913393,71519608,73027261,74467196,75844193,77183580,78444407,79636468,80762075,81838350,82862126,83856990,84780261,85655110,86489880,87280599,88021698,88724981,89397685,90033172,90625517,91187960,91744043,92250642,92914639,93772107,94535443,95215296,95815947,96348141,96817774,97228763,97591438,97910127,98189569,98439872,98652025,98836368,99007287,99254377,99437449,99569748,99665436,99733431,99780876,99813393,99833213,99848105,99855950,99860426,99863545,99862826,99857931,0,0,0,0,0,0,
63832843,64065220,64409824,65087627,65989463,66874039,68135051,69793228,71373362,72916990,74364286,75746093,77067466,78350327,79555443,80694533,81772447,82802493,83780725,84729532,85607545,86441754,87239012,87990936,88698874,89371332,90014456,90612784,91174415,91713426,92238522,92698232,93324590,94158584,94880733,95521673,96088625,96590024,97033046,97420357,97759888,98057631,98318591,98551924,98749379,98920660,99075906,99303733,99473850,99597188,99684795,99746856,99790189,99819472,99837664,99850929,99858136,99861610,99863667,99862534,99866774,0,0,0,0,0,0,
65519914,65743999,66076160,66729362,67598327,68451979,69666087,71258582,72775920,74257498,75646332,76971970,78238703,79466178,80617824,81707899,82739322,83724812,84659482,85539481,86374183,87191167,87952733,88668881,89343676,89987717,90597685,91163431,91699698,92211615,92693247,93146221,93737936,94525015,95206603,95811291,96346283,96817076,97234351,97598553,97916513,98195659,98438848,98651390,98839548,98999623,99141429,99350872,99507718,99622136,99702549,99758845,99798271,99824370,99841207,99852993,99859432,99862120,99863089,99861759,99873105,0,0,0,0,0,0,
67145018,67360928,67680849,68310081,69147642,69969130,71133964,72663809,74121764,75542860,76874943,78146547,79360206,80532934,81633209,82675250,83661562,84604853,85497593,86334441,87130558,87906771,88633343,89314706,89960095,90573192,91147798,91686875,92198946,92685046,93141195,93569457,94128366,94870860,95514245,96084678,96588786,97031569,97422645,97764274,98062255,98324689,98551818,98749366,98923471,99072062,99201779,99396291,99540235,99644648,99718770,99770524,99806099,99828738,99844422,99854938,99859957,99862140,99862818,99861089,99877015,0,0,0,0,0,0,
68709119,68917031,69225607,69831424,70635991,71423515,72542228,74012904,75412297,76774603,78051363,79270336,80432817,81553385,82574546,83563915,84541729,85445595,86296321,87092306,87852451,88590825,89281790,89929841,90545775,91125568,91669588,92184616,92672555,93134826,93564705,93969181,94497972,95198121,95805067,96342442,96816822,97234373,97600221,97920002,98199251,98444822,98657015,98840165,99001535,99138824,99258772,99439508,99570312,99665701,99733714,99780392,99812947,99832946,99847647,99857087,99861170,99862465,99861383,99860446,99878650,0,0,0,0,0,0,
70216544,70414805,70710618,71291766,72063126,72819447,73895015,75307147,76648958,77954062,79177860,80345498,81458188,82508938,83503660,84448629,85383841,86246743,87055844,87814397,88540342,89241867,89899092,90516529,91099288,91648372,92166889,92656164,93121516,93559233,93964527,94348422,94849131,95493708,96079588,96584784,97029573,97422517,97765637,98064939,98326958,98556435,98754010,98923953,99071800,99199538,99308285,99477846,99597704,99685074,99747169,99789919,99819733,99837293,99850537,99858533,99862328,99863059,99861175,99866056,99879465,0,0,0,0,0,0,
72272808,72459581,72737952,73285855,74014069,74727678,75741303,77071514,78333981,79560267,80710411,81807709,82850722,83832313,84761644,85646174,86519693,87323600,88077686,88787996,89467061,90117950,90707568,91276362,91838854,92349219,92830478,93284728,93714823,94116704,94492119,94850252,95310065,95901269,96436824,96898726,97306598,97663609,97974746,98247568,98485895,98693744,98872375,99025512,99159536,99273582,99364592,99515867,99624692,99702778,99758058,99796593,99822182,99837057,99848793,99855453,99858121,99858573,99857417,99874378,99880131,0,0,0,0,0,0,
74937961,75111381,75368286,75873249,76543479,77199509,78130646,79353075,80511920,81634662,82689019,83693079,84642552,85535579,86382827,87189436,87981299,88709567,89394121,90040500,90660136,91238518,91776115,92288074,92790924,93250650,93683362,94091455,94474379,94833492,95171448,95492265,95898915,96425298,96889428,97303224,97663012,97975619,98247717,98486081,98694020,98873745,99027462,99160500,99275980,99374209,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
77390743,77550220,77786382,78250026,78865211,79466889,80320615,81442694,82504739,83530581,84493601,85408857,86272268,87085247,87858183,88592044,89308052,89966847,90588019,91175859,91734956,92253642,92737828,93200346,93649548,94061795,94447775,94811247,95153324,95474941,95776853,96061033,96421628,96888998,97298774,97655571,97967332,98248095,98486112,98693277,98873228,99027895,99160706,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
79641070,79787217,80003948,80429101,80992765,81543562,82325034,83353370,84324851,85239787,86135688,86968620,87754780,88495052,89198657,89864517,90509863,91107118,91668545,92197622,92698558,93164114,93603113,94020798,94410694,94773531,95115138,95455618,95760414,96045165,96313457,96565400,96883704,97295698,97657390,97971757,98244393,98487418,98693873,98873214,99027512,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
81702356,81835975,82034164,82423223,82938756,83442461,84156847,85096621,85983009,86816400,87630235,88387038,89101564,89774908,90413980,91016976,91581299,92113312,92614183,93114600,93563930,93983613,94380114,94755338,95101417,95423918,95730160,96030598,96300336,96551868,96788913,97010719,97289648,97653390,97971749,98247034,98485130,98691904,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
83587060,83709142,83890359,84245585,84715725,85175578,85827591,86684074,87491698,88250927,88989405,89675550,90323095,90933971,91514687,92058759,92564803,93041401,93490740,93938996,94341008,94718422,95073971,95407218,95715246,96003254,96276470,96540025,96765990,96987842,97209595,97404434,97648207,97968157,98246363,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
85308828,85419560,85585048,85908446,86337014,86755978,87349202,88128148,88862651,89553305,90222458,90843329,91428178,91981644,92505841,92994024,93446818,93875282,94279776,94678398,95037400,95374093,95690306,95986010,96260737,96517934,96760823,96984756,97191936,97387129,97580329,97751120,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
86877753,86978379,87128912,87423415,87813411,88194505,88732806,89439562,90106309,90733527,91338330,91898298,92425664,92923536,93393336,93831223,94238578,94623674,94986046,95339989,95658741,95957450,96237862,96500349,96744199,96972687,97187340,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
88304384,88396483,88532783,88800589,89154750,89500475,89988161,90627979,91232596,91800220,92346988,92829588,93297778,93768259,94189298,94582407,94948004,95293166,95615586,95928842,96210189,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
89600589,89684385,89807594,90049949,90371233,90684030,91124638,91702954,92249792,92763258,93241506,93686213,94105149,94525972,94902476,95253539,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
90774980,90850803,90962639,91181497,91471730,91754315,92151203,92673452,93165835,93627309,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
91837026,91905310,92006278,92203685,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
]
