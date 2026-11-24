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
    assert n in (96, 104, 320)
    U = UNIT
    pts = [0, U // 4, U // 2, U, (3 * U) // 2] if n == 104 else [0, U // 8, U // 4, U // 2, U, (3 * U) // 2]
    if n == 96:
        limits = (24, 48, 96, 144)
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
