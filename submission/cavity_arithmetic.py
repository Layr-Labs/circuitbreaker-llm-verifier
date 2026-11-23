"""Directed integer quadrature for the guarded cubic-tree probabilities."""
from bisect import bisect_right
from functools import lru_cache
from math import isqrt

SPACE = 2**24
SCALE = 10**10
PREC = 2**54
SQ2_LO = isqrt(2 * PREC**2)
SQ2_HI = SQ2_LO + 1


def ceil_div(num, den):
    return -((-num) // den)


def root2(lo, hi):
    return (lo * (SQ2_LO if lo >= 0 else SQ2_HI) // PREC,
            ceil_div(hi * (SQ2_HI if hi >= 0 else SQ2_LO), PREC))


# The endpoint mesh is in units of 1/512 throughout.
GRID = tuple(i * SPACE // 512 for i in
             [0, 1, 2, 4, 8, 12]
             + list(range(16, 193, 8))
             + list(range(208, 385, 16))
             + list(range(416, 641, 32))
             + list(range(704, 1025, 64)))
LAST = len(GRID) - 1
DISTANCES = tuple(r - l for l, r in zip(GRID, GRID[1:]))
assert len(GRID) == 55 and GRID[-1] == 2 * SPACE


@lru_cache(maxsize=None)
def _exp(ceil_millionth):
    # A cubic alternating Taylor truncation followed by downward squaring.
    denominator = 256 * 10**6
    k = ceil_millionth
    value = (6 * denominator**3 - 6 * k * denominator**2
             + 3 * k**2 * denominator - k**3)
    value = value * SCALE // (6 * denominator**3)
    for _ in range(8):
        value = value * value // SCALE
    return value


def exp_lower(num, den):
    k = ceil_div(num * 10**6, den)
    if k > 14 * 10**6:
        return 0
    return _exp(k)


G_LO = 27639531957  # downward integer at SCALE=10^10


@lru_cache(maxsize=None)
def _radial(step):
    return G_LO * exp_lower(6 * step**2, 40000**2) // SCALE


def radial(lo, hi):
    step = ceil_div(max(abs(lo), abs(hi)) * 40000, SPACE)
    return _radial(step)


MESH = 2048
MASS = [0]
for i in range(1, 2 * MESH + 1):
    value = exp_lower(2 * (3*i*i - 3*i + 1), MESH*MESH)
    MASS.append(MASS[-1] + G_LO * value // (SCALE * MESH))


def half_mass(lower):
    return MASS[max(0, min(2*MESH, lower * MESH // SPACE))]


def gm4(a, b, c, d):
    return isqrt(isqrt(a * b) * isqrt(c * d))


def gm8(a, b, c, d, e, f, g, h):
    return isqrt(gm4(a, b, c, d) * gm4(e, f, g, h))


def geometric_interpolation(a, b, delta, width):
    if delta == 0:
        return a
    if delta == width:
        return b
    if not a or not b:
        return 0
    if a > b:
        a, b = b, a
        delta = width - delta
    for ratio in (48, 8):
        if (b - a) * ratio <= a:
            break
        mid = isqrt(a * b)
        delta *= 2
        if delta <= width:
            b = mid
        else:
            a = mid
            delta -= width
    return width * a * b // ((width - delta) * b + delta * a)


def row_at(row, value):
    j = bisect_right(GRID, value) - 1
    if j == LAST or GRID[j] == value:
        return row[j]
    return geometric_interpolation(row[j], row[j+1], value-GRID[j], DISTANCES[j])


def interp_row(matrix, i, lo, hi):
    if lo < 0 or hi > GRID[-1]:
        return 0
    row = matrix[i]
    j = bisect_right(GRID, lo) - 1
    if j < LAST and hi <= GRID[j+1]:
        value = lo if row[j] <= row[j+1] else hi
        return geometric_interpolation(row[j], row[j+1], value-GRID[j], DISTANCES[j])
    return min(row_at(row, lo), row_at(row, hi))


def interval_mesh(start, end):
    if start >= end:
        return ()
    i = bisect_right(GRID, start)
    j = bisect_right(GRID, end - 1)
    return (start, *GRID[i:j], end)


def integral(nodes, values):
    return sum((b-a)*isqrt(x*y) for a,b,x,y in
               zip(nodes, nodes[1:], values, values[1:])) // SPACE


def matrix(function):
    return [[function(i,j) for j in range(LAST+1)] for i in range(LAST+1)]


RGRID = tuple(root2(x,x) for x in GRID)


def k_mid(i, j):
    lo, hi = RGRID[j]
    return lo + GRID[i]//2, hi + GRID[i]//2


def siblings(i, j):
    lo, hi = RGRID[j]
    return lo - GRID[i]//2, hi - GRID[i]//2


def k_weight(a, b, d, guard=None):
    mlo, mhi = root2(b,b)
    mlo += a//2
    mhi += a//2
    f = radial(d+mlo,d+mhi)
    if guard is None or not f:
        return f
    if guard == "U":
        u, _ = root2(a+d,a+d)
        return f * half_mass(7*b//2 + u) // SCALE
    name, table = guard
    lo, hi = 2*mlo+d, 2*mhi+d
    return f * interp_row(table, name, lo, hi) // SCALE


def k_steps(guards):
    result = []
    for guard in guards:
        steps = []
        for i,a in enumerate(GRID):
            row = []
            for j,b in enumerate(GRID):
                descriptor = guard if isinstance(guard, str) or guard is None else (j,guard)
                values = [k_weight(a,b,d,descriptor) for d in GRID]
                row.append(tuple((r, w*isqrt(values[r]*values[r+1])) for r,w in enumerate(DISTANCES)
                                 if values[r] and values[r+1]))
            steps.append(row)
        result.append(steps)
    return result


def apply_k(function, steps):
    return matrix(lambda i,j: sum(weight * isqrt(function[j][r]*function[j][r+1])
                                  for r,weight in steps[i][j]) // (SPACE*SCALE))


def side_values(i, j, one_child):
    mlo, mhi = siblings(i,j)
    lower = max(0, -(2 if one_child else 1)*mlo)
    nodes = interval_mesh(lower,2*SPACE)
    weights = tuple(radial(d+mlo,d+mhi) for d in nodes)
    return nodes,weights


def apply_s(functions, one_child, side_guard=None):
    output = [matrix(lambda i,j: 0) for _ in functions]
    for i,p in enumerate(GRID):
        for j,u in enumerate(GRID):
            nodes,weights = side_values(i,j,one_child)
            if not any(weights):
                continue
            mlo,mhi = siblings(i,j)
            if side_guard == "U":
                weights = tuple(w*half_mass(7*u//2 + root2(d-p,d-p)[0])//SCALE
                                for d,w in zip(nodes,weights))
            elif side_guard is not None:
                weights = tuple(w*interp_row(side_guard,j,2*mlo+d,2*mhi+d)//SCALE
                                for d,w in zip(nodes,weights))
            segments = tuple((r,b-a,weights[r],weights[r+1]) for r,(a,b) in
                             enumerate(zip(nodes,nodes[1:])) if weights[r] and weights[r+1])
            first = bisect_right(GRID,nodes[0])
            for result,f in zip(output,functions):
                # Apart from the clipped lower endpoint, side_values uses
                # the consecutive grid nodes from first through LAST.
                values = [row_at(f[j],nodes[0]), *f[j][first:]]
                result[i][j] = sum(width*isqrt((left*values[r]//SCALE)*
                                               (right*values[r+1]//SCALE))
                                   for r,width,left,right in segments)//SPACE
    return tuple(output)


def apply_d(branches, pairs):
    outputs = [matrix(lambda i,j: 0) for _ in pairs]
    mask = sorted({k for pair in pairs for k in pair})
    for i in range(LAST+1):
        for j in range(LAST+1):
            mlo,mhi=siblings(i,j)
            nodes = interval_mesh(0, min(-mhi, 2*SPACE))
            if not nodes:
                continue
            weights = [radial(h,h) for h in nodes]
            plus = {k:[interp_row(branches[k],j,-mhi+h,-mlo+h) for h in nodes] for k in mask}
            minus = {k:[interp_row(branches[k],j,-mhi-h,-mlo-h) for h in nodes] for k in mask}
            for output,(a,b) in zip(outputs,pairs):
                values = [w*x*y//(SCALE*SCALE) for w,x,y in zip(weights,plus[a],minus[b])]
                output[i][j] = integral(nodes,values)
    return tuple(outputs)


def side_pair_values(tables, i, j, h, upper_child):
    mlo,mhi = k_mid(i,j)
    sign = 1 if upper_child else -1
    zlo,zhi = mlo+sign*h,mhi+sign*h
    alo,ahi = root2(GRID[i]//2+sign*h,GRID[i]//2+sign*h)
    u = half_mass(3*GRID[j]//2+alo)
    return (u,*(interp_row(table,j,zlo,zhi) for table in tables))


POINT_PAIRS = ((0,0),(5,0),(0,3),(1,0),(0,1),(6,0),(0,4),(2,0),(0,2),
               (8,0),(7,0),(3,0),(4,0),(9,0),(10,0),(1,3),(5,1),(1,1),(2,3),(5,2))
POINT_UPPER = (0,1,2,5,6)
OUT = tuple(range(5))
V = tuple(range(5,11))
ALL_PAIRS = tuple((i,j) for i in OUT+V for j in OUT+V if i in OUT or j in OUT)
H_CUTOFF = 13*SPACE//10


def point_paths(guards):
    result = [matrix(lambda i,j: 0) for _ in POINT_PAIRS]
    upper = tuple(guards[k-1] for k in POINT_UPPER[1:])
    for i in range(LAST+1):
        for j in range(LAST+1):
            start,end=0,min(H_CUTOFF,k_mid(i,j)[0])
            nodes = interval_mesh(start,end)
            if not nodes:
                continue
            lower_values = [side_pair_values(guards,i,j,h,False) for h in nodes]
            upper_values = [side_pair_values(upper,i,j,h,True) for h in nodes]
            weights = [radial(h,h) for h in nodes]
            for output,(a,b) in zip(result,POINT_PAIRS):
                values = [d*lo[a]*hi[b]//(SCALE*SCALE)
                          for d,lo,hi in zip(weights,lower_values,upper_values)]
                output[i][j] = integral(nodes,values)
    return tuple(result)


def endpoint_cells(guards):
    panels = 16
    # Each grid node belongs to up to four cells. Reuse its directed
    # temporal means before taking the four-corner means of each cell.
    @lru_cache(maxsize=3 * (LAST + 1))
    def time_edges(i, j):
        cap = min(H_CUTOFF, k_mid(i,j)[0]) // panels * panels
        edges = []
        previous = None
        for t in range(panels + 1):
            h = cap*t // panels
            density = radial(h,h)
            lower = side_pair_values(guards,i,j,h,False)
            upper = lower if t == 0 else side_pair_values(guards,i,j,h,True)
            if previous is not None:
                pd, pl, pu = previous
                edges.append((isqrt(pd*density),
                              tuple(isqrt(a*b) for a,b in zip(pl,lower)),
                              tuple(isqrt(a*b) for a,b in zip(pu,upper))))
            previous = density,lower,upper
        return cap, edges

    output = [[0]*LAST for _ in range(LAST)]
    for i in range(LAST):
        for j in range(LAST):
            corners = [(i,j),(i,j+1),(i+1,j),(i+1,j+1)]
            caps, edges = zip(*(time_edges(a,b) for a,b in corners))
            if not all(caps):
                continue
            jacobian = gm4(*caps)
            total = 0
            for ends in zip(*edges):
                densities, lower, upper = zip(*ends)
                density = gm4(*densities)
                left = [gm4(*values) for values in zip(*lower)]
                right = [gm4(*values) for values in zip(*upper)]
                mass = (sum(left[k] for k in OUT)*sum(right)
                        + sum(left[k] for k in V)*sum(right[k] for k in OUT))
                total += density*mass
            output[i][j] = jacobian*total//(panels*SPACE*SCALE*SCALE)
    return output


def cell_k(guards):
    # Adjacent cells use the same directed kernel edge at a shared corner.
    @lru_cache(maxsize=3 * (LAST+1) * len(guards))
    def corner_edges(i, j, guard_index):
        guard = guards[guard_index]
        if guard is not None and not isinstance(guard,str):
            guard = j, guard
        weights = [k_weight(GRID[i],GRID[j],d,guard) for d in GRID]
        return tuple(isqrt(a*b) for a,b in zip(weights,weights[1:]))

    kernel = []
    for i in range(LAST):
        row = []
        for j in range(LAST):
            corners = [(i,j),(i,j+1),(i+1,j),(i+1,j+1)]
            values = [0]*LAST
            for guard_index in range(len(guards)):
                edges = [corner_edges(a,b,guard_index) for a,b in corners]
                for d,corners_at_d in enumerate(zip(*edges)):
                    values[d] += DISTANCES[d]*gm4(*corners_at_d)
            row.append(tuple((d,w) for d,w in enumerate(values) if w))
        kernel.append(row)
    return kernel


def apply_cell_k(values,kernel):
    return [[sum(w*values[j][k] for k,w in kernel[i][j])//(SPACE*SCALE)
             for j in range(LAST)] for i in range(LAST)]


def branch_tables(progress=lambda _:None):
    k,ku = k_steps((None,"U"))
    branches = [matrix(lambda i,j:half_mass(k_mid(i,j)[0]))]
    for _ in range(6):
        branches.append(apply_k(branches[-1],k))
    even = tuple(branches[t] for t in (0,2,4,6))
    odd = tuple(branches[t] for t in (1,3,5))
    A0 = apply_s(even,False)
    A = A0+apply_d(branches,tuple((i,j) for i in (1,3) for j in (0,2,4)))
    W1,W3,W5 = apply_s(odd,True)
    W2, = apply_d(branches,((1,1),))
    v0,v2,v4 = apply_s(even[:3],True,"U")
    v01,v21 = apply_s(even[:2],True,W1)
    v03, = apply_s(even[:1],True,W3)
    guards = (W1,W3,W5,W2,v0,v2,v4,v01,v21,v03)
    progress("branch and guard tables")
    p = point_paths(guards)
    hop1 = tuple(apply_k(apply_k(p[t],k),ku) for t in (0,1,2,3,5,7,9,15,16))
    hop2 = tuple(apply_k(apply_k(hop1[t],k),ku) for t in (0,1,3))
    hop3 = tuple(apply_k(apply_k(m,k),ku) for m in hop2[:2])
    P = apply_s(p+hop1+hop2+hop3,True)
    O_guard = tuple(apply_k(f,ku) for f in even)
    kw1,kw3 = k_steps((W1,W3))
    O_guard += tuple(apply_k(f,kw1) for f in even[:3])+(apply_k(even[0],kw3),)
    progress("point itineraries")
    plain = cell_k((None,))
    guarded = cell_k(("U",W1,W3))
    cell = endpoint_cells(guards)
    E = [list(row) for row in cell]
    for _ in range(4):
        cell = apply_cell_k(apply_cell_k(cell,plain),guarded)
        E = [[a+b for a,b in zip(base,growth)] for base,growth in zip(E,cell)]
    progress("cell itineraries")
    return A,A0,P,even,odd,O_guard,E
