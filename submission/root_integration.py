"""Directed two-arm union bounds and batched log-concave interpolation."""
from array import array
from bisect import bisect_right
from functools import lru_cache
from math import isqrt
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


def row_lookup(family, parent_index, grid, scale):
    """Evaluate directed lower tangents for all itineraries in one row."""
    side = len(grid)
    start = parent_index * side
    rows = [table[start:start+side] if start+side <= len(table) else (0,)*side for table in family]
    columns = list(zip(*rows))
    widths = [b-a for a, b in zip(grid, grid[1:])]
    size = len(family)
    zeros = (0,) * size
    layout = Struct(f"<{size}Q")
    subdivisions = 8
    bits = 16
    unit = 1 << bits
    mask = int.from_bytes(layout.pack(*([(1 << (64-bits))-1] * size)), "little")
    assert max(map(max, columns)) <= scale and 3*scale*unit < 1 << 64

    assert all(width % subdivisions == 0 for width in widths)

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
        low, up, down = [], [], []
        for a, b in zip(node(j, k), node(j, k+1), strict=True):
            if a < b:
                low.append(a)
                up.append(2*a*(b-a)//(b+a) if a else 0)
                down.append(0)
            else:
                low.append(b)
                up.append(0)
                down.append(2*b*(a-b)//(a+b) if b and a != b else 0)
        return pack(low) << bits, pack(up), pack(down)

    def evaluate(lo, hi):
        j = bisect_right(grid, lo) - 1
        if j == side - 1:
            return columns[j]
        width = widths[j] // subdivisions
        k, offset = divmod(lo-grid[j], width)
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


def _edge_means(left, right):
    assert len(left) == len(right)
    return tuple(map(isqrt, map(mul, left, right)))


def square_cell_mean(plane, i, j):
    """Downward geometric mean of the four corners of a square cell."""
    return isqrt(isqrt(plane[i][j] * plane[i][j+1])
                 * isqrt(plane[i+1][j] * plane[i+1][j+1]))


def companion_mean(previous, current, i, j):
    """Sum the corner floors for the disjoint companion guards."""
    corners = (
        previous[i][j], previous[i][j+1],
        previous[i+1][j], previous[i+1][j+1],
        current[i][j], current[i][j+1],
        current[i+1][j], current[i+1][j+1],
    )
    return sum(
        isqrt(isqrt(isqrt(a*b)*isqrt(c*d)) * isqrt(isqrt(e*f)*isqrt(g*h)))
        for a, b, c, d, e, f, g, h in zip(*corners, strict=True)
        if all((a, b, c, d, e, f, g, h))
    )


def mixed_arm_floors(grid, root_limit, center_cells, near_tables, far_tables,
                     half_width, density, scale, coordinate_scale):
    """Floor the mixed negative-arm chords on the target rectangles."""
    panels = 8
    cside = len(grid) - 1
    mixed_c2 = [0] * (center_cells*cside)
    mixed_support = set()
    for j in range(root_limit - 1):
        v0, v1 = grid[j], grid[j+1]
        near = [row_lookup(near_tables, k, grid, scale) for k in (j,j+1)]
        far = [row_lookup(far_tables, k, grid, scale) for k in (j,j+1)]
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
    n, panels = 64, 12
    S, U = branch.SCALE, branch.UNIT
    RLO, RHI = branch.ROOT_LO, branch.ROOT_HI
    EXP_DEN = branch.EXP_DEN
    coord = n * U
    coord2 = coord ** 2
    knots = [i*U//2 for i in range(21)]
    knots.extend(i*U for i in range(11, 33))
    knots.extend(i*U for i in range(34, 61, 2))
    knots.extend(i*U for i in range(64, 3*n//2+1, 4))
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
            hi_lookup = None if upper is None else row_lookup((upper,),j,knots,S)
            lo_lookup = None if lower is None else row_lookup((lower,),j,knots,S)
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
                    large=S if hi_lookup is None else hi_lookup(a,b)[0]
                    small=S if lo_lookup is None else lo_lookup(c,d)[0]
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
        x_hi = ceil_div(2*RHI*x,U)
        current = []
        for y in knots[:iy+1]:
            row=[]
            for z in knots[:iy+1]:
                if y <= 0 or z <= 0 or x_lo <= y+z:
                    row.append(0)
                    continue
                a=max(abs(3*y-x_lo),abs(3*y-x_hi))
                b=max(abs(y+2*z-x_lo),abs(y+2*z-x_hi))
                q=ceil_div((x*x+a*a+3*b*b)*EXP_DEN,2*coord2)
                row.append(branch.C3_NUM*branch.exp_floor(q)//S)
            current.append(row)
        if previous is not None:
            dx=knots[i]-knots[i-1]
            for j in range(iy):
                arm=rays[i-1][j]
                if not arm:
                    continue
                dy=knots[j+1]-knots[j]
                for k in range(iy):
                    dz=knots[k+1]-knots[k]
                    w=isqrt(square_cell_mean(previous,j,k)*square_cell_mean(current,j,k))
                    if not w:
                        continue
                    vol=3*(dx//(U//2))*(dy//(U//2))*(dz//(U//2))
                    joint[0]+=vol*(w*arm*rays[i-1][k]//S)
                    joint[1]+=2*vol*(w*arm*stable[i-1][k]//S)
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
    pts = [0, U // 8, U // 4, U // 2, U, (3 * U) // 2]
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
