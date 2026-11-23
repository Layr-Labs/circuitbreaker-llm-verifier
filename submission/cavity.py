"""Integer lower bounds for four finite signed cavity events on the cubic tree."""
from functools import lru_cache
from math import isqrt

SCALE = 4_000_000
UNIT = 100_000
ROOT_LO = 141_421
ROOT_HI = 141_422
HALF_CDF_MESH = 1024


def ceil_div(a, b):
    return -(-a // b)


def harmonic(a, b, t):
    if not t or a == b:
        return a
    if t == UNIT:
        return b
    if not a or not b:
        return 0
    w = 256*t // UNIT if b > a else ceil_div(256*t, UNIT)
    return 256*a*b // ((256-w)*b + w*a)


@lru_cache(maxsize=None)
def exp_floor(i):
    if i > 400_000:
        return 0
    h = 640_000
    u, l = ceil_div(i*SCALE, h), i*SCALE // h
    v = SCALE - u + l*i // (2*h) - ceil_div(ceil_div(u*i, 2*h)*i, 3*h)
    for _ in range(6):
        v = v*v // SCALE
    return v


def product(*values):
    acc = SCALE
    for v in values:
        acc = acc*v // SCALE
    return acc


def cavity(n=32):
    """Return four integer numerators and their common denominator."""
    m = 2*n
    side = m+1
    size = side*side
    bound = m*UNIT
    mesh = HALF_CDF_MESH

    @lru_cache(maxsize=None)
    def density(v):
        t = ceil_div(abs(v), 10*n)
        return 27_638*exp_floor(ceil_div(6*t*t, 10_000)) // 10_000

    mass = [0]
    for i in range(1, 2*mesh+1):
        e = exp_floor(ceil_div(20_000*(3*i*i-3*i+1), mesh*mesh))
        mass.append(mass[-1] + 27_638*e // (10_000*mesh))

    def half_mass(r):
        return mass[min(2*mesh, mesh*r//(n*UNIT))] if r > 0 else 0

    def blend(mat, x, y):
        i, f = divmod(x, UNIT)
        j, h = divmod(y, UNIT)
        k = i*side+j
        if not mat[k]:
            return 0
        a = harmonic(mat[k], mat[k+1], h) if h else mat[k]
        if not f:
            return a
        b = harmonic(mat[k+side], mat[k+side+1], h) if h else mat[k+side]
        return harmonic(a, b, f)

    def lower(mat, xl, xh, yl, yh):
        if xl < 0 or yl < 0 or xh >= bound or yh >= bound:
            return 0
        a = blend(mat, xl, yl)
        if not a:
            return 0
        if xh != xl:
            a = min(a, blend(mat, xh, yl))
        if yh != yl:
            a = min(a, blend(mat, xl, yh))
            if xh != xl:
                a = min(a, blend(mat, xh, yh))
        return a

    points = range(side)
    # Reuse the positive initial portion of the density stencil. Dropping
    # any later nonnegative contributions preserves the lower bound.
    kernels = []
    for i in points:
        for j in points:
            base = ROOT_HI*j+i*UNIT//2
            kernels_row = []
            for k in points:
                value = density(base+k*UNIT)
                if not value:
                    break
                kernels_row.append(value)
            kernels.append(kernels_row)
    branches = [[half_mass(ROOT_LO*j + i*UNIT//2) for i in points for j in points]]
    for _ in range(4):
        old = branches[-1]
        mat = []
        for i in points:
            for j in points:
                kernel = kernels[i*side+j]
                if not kernel:
                    mat.append(0)
                    continue
                row = j*side
                previous = kernel[0]*old[row]//SCALE
                total = 0
                for k in range(1,len(kernel)):
                    current = kernel[k]*old[row+k]//SCALE
                    total += isqrt(previous*current)
                    previous = current
                mat.append(total//n)
        branches.append(mat)

    def same(mat, strict):
        out = [0]*size
        for i in points:
            for j in points:
                ml, mh = ROOT_LO*j-i*UNIT//2, ROOT_HI*j-i*UNIT//2
                c = max(0, -(2 if strict else 1)*ml)
                if c >= bound:
                    continue
                k, rem = divmod(c, UNIT)
                row = j*side
                previous = harmonic(mat[row+k], mat[row+k+1], rem)*density(max(c+mh, 0))//SCALE
                total = 0
                for k in range(k+1, side):
                    edge = k*UNIT
                    weight = density(edge+mh)
                    if not weight:
                        break
                    current = weight*mat[row+k]//SCALE
                    total += isqrt(previous*current)*(edge-c)
                    previous, c = current, edge
                out[i*side+j] = total//(UNIT*n)
        return out

    even = branches[::2]
    odd = branches[1::2]
    same_even = [same(mat, False) for mat in even]

    closed = [0]*size
    private = [0]*size
    for i in points:
        for j in points:
            a = 3*UNIT*j//2+70_710*i
            previous = total = 0
            for k in points:
                current = density(k*UNIT)*half_mass(a-ROOT_HI*k)//SCALE*half_mass(a+ROOT_LO*k)//SCALE
                total += isqrt(previous*current)
                previous = current
            closed[i*side+j] = total//n

            ml, mh = ROOT_LO*j-i*UNIT//2, ROOT_HI*j-i*UNIT//2
            c = max(0, -2*ml)
            if c >= bound:
                continue
            level = 7*UNIT*j//2-ROOT_HI*i
            previous = density(c+mh)*half_mass(level+ROOT_LO*c//UNIT)//SCALE
            total = 0
            for k in range(c//UNIT+1, side):
                edge = k*UNIT
                current = density(edge+mh)*half_mass(level+ROOT_LO*k)//SCALE
                total += isqrt(previous*current)*(edge-c)
                previous, c = current, edge
            private[i*side+j] = total//(UNIT*n)

    # Disjoint extra classes for a protected continuation.
    weak = same(branches[1], True)
    rare = [0]*size
    def tail_mass(r):
        start = min(2*mesh, max(0, ceil_div(mesh*r, UNIT*n)))
        return mass[-1]-mass[start]
    for i in points:
        for j in points:
            ml, mh = ROOT_LO*j-i*UNIT//2, ROOT_HI*j-i*UNIT//2
            c = max(0, -2*ml)
            if c >= bound:
                continue
            row = j*side
            def term(edge):
                ix = min(m-1, edge//UNIT)
                u = harmonic(branches[1][row+ix], branches[1][row+ix+1], edge-ix*UNIT)
                radius = 7*j*UNIT//2-ROOT_LO*i+ceil_div(ROOT_HI*edge, UNIT)
                return density(edge+mh)*u//SCALE*tail_mass(radius)//SCALE
            previous = term(c)
            total = 0
            for k in range(c//UNIT+1, side):
                edge = k*UNIT
                value = term(edge)
                total += isqrt(previous*value)*(edge-c)
                previous, c = value, edge
            rare[i*side+j] = total//(UNIT*n)
    tails = [[0]*size for _ in range(7)]
    for i in points:
        for j in points:
            ml, mh = ROOT_LO*j+i*UNIT//2, ROOT_HI*j+i*UNIT//2
            upper = min(bound, ml)
            if not upper:
                continue
            previous, totals, left = [0]*7, [0]*7, 0
            r = 3*UNIT*j//2+70_710*i
            for k in range(ceil_div(upper, UNIT)+1):
                edge = min(k*UNIT, upper)
                density_here = density(edge)
                up = half_mass(r+ROOT_LO*edge//UNIT)
                down = half_mass(r-ceil_div(ROOT_HI*edge, UNIT))
                def side_bounds(mat, e):
                    return lower(mat, j*UNIT, j*UNIT, ml+e, mh+e)
                vp, vm = side_bounds(private, edge), side_bounds(private, -edge)
                rp, rm = side_bounds(rare, edge), side_bounds(rare, -edge)
                wm = side_bounds(weak, -edge)
                pairs = ((up, vm), (down, vp), (wm, vp), (rm, up),
                         (rp, vm), (rp, down), (rm, rp))
                for t, (a,b) in enumerate(pairs):
                    current = density_here*a//SCALE*b//SCALE
                    totals[t] += isqrt(previous[t]*current)*(edge-left)
                    previous[t] = current
                left = edge
            for mat, total in zip(tails, totals):
                mat[i*side+j] = total//(UNIT*n)

    def both(f, g):
        out = [0]*size
        for i in points:
            for j in points:
                s0, s1 = i*UNIT//2-ROOT_HI*j, i*UNIT//2-ROOT_LO*j
                if s0 <= 0:
                    continue
                previous = total = left = 0
                for k in range(ceil_div(s0, UNIT)+1):
                    edge = min(k*UNIT, s0)
                    a = lower(f, j*UNIT, j*UNIT, s0+edge, s1+edge)
                    b = lower(g, j*UNIT, j*UNIT, s0-edge, s1-edge)
                    current = density(edge)*a//SCALE*b//SCALE
                    total += isqrt(previous*current)*(edge-left)
                    previous, left = current, edge
                out[i*side+j] = total//(UNIT*n)
        return out

    def continue_arm(mat, collateral):
        out = [0]*size
        for i in points:
            for j in points:
                row = j*side
                previous = total = 0
                for k, weight in enumerate(kernels[i*side+j]):
                    value = weight*mat[row+k]//SCALE
                    if collateral:
                        value = value*half_mass(7*j*UNIT//2+ROOT_LO*(i+k))//SCALE
                    total += isqrt(previous*value)
                    previous = value
                out[i*side+j] = total//n
        return out

    ends = [closed, *tails[:2]]
    protected_odd = [continue_arm(mat, True) for mat in even]
    two = [continue_arm(continue_arm(mat, False), True) for mat in ends]
    four = [continue_arm(continue_arm(mat, False), True) for mat in two]
    guarded = ends + tails[2:] + two + four
    protected = [same(mat, True) for mat in guarded]
    mixed = protected + [both(a,b) for a in even for b in even]
    for mat in guarded:
        mixed.extend((both(mat, branches[1]), both(branches[1], mat)))
    for mat in ends:
        mixed.extend((both(mat, branches[3]), both(branches[3], mat)))
    # For each group take the geometric mean of each separate itinerary
    # before adding them. The four products count disjoint event classes.
    positive = same_even + [both(a, b) for a in odd for b in even]
    groups = [[positive, positive, even, mixed, mixed, guarded],
              [even, odd, even, same_even, protected, guarded, protected_odd, protected_odd, guarded]]
    counts = [0]*4
    for kind, group in enumerate(groups, 1):
        # Use a single cached evaluation for identical local classes.
        index = {id(mats): j for j, mats in enumerate(group)}
        shared = [index[id(mats)] for mats in group]
        offsets = [1]
        for mats in group:
            offsets.append(offsets[-1]+len(mats))
        dim = offsets[-1]
        width = (n+1)*dim
        for i in range((n if kind==1 else n//2)+1):
            @lru_cache(maxsize=None)
            def values(group_id, xl, xh, yl, yh):
                mats = group[group_id]
                if xl < 0 or yl < 0 or xh >= bound or yh >= bound:
                    return (0,)*len(mats)
                # A zero lower-left value forces the harmonic interpolant to zero.
                corner = (xl//UNIT)*side + yl//UNIT
                return tuple(lower(mat, xl, xh, yl, yh) if mat[corner] else 0 for mat in mats)

            plane = []
            for j in range(n+1):
                for k in range(n+1):
                    if kind == 1:
                        y, zz = j*UNIT, k*UNIT
                        xl, xh = 35_355*i*(j+k)//n, ceil_div(35_355*i*(j+k), n)
                        bl = max(0, y+zz-ceil_div(2*ROOT_HI*xh, UNIT))
                        bh = y+zz-2*ROOT_LO*xl//UNIT
                        jac = SCALE*35_355*(j+k)//(UNIT*n)
                        coords = [(y,y),(zz,zz),(bl,bh),(y,y),(zz,zz),(bl,bh)]
                    else:
                        xl = xh = i*UNIT
                        y, zz = -j*UNIT, -k*UNIT
                        bl, bh = 2*ROOT_LO*i-y-zz, 2*ROOT_HI*i-y-zz
                        jac = SCALE
                        coords = [(-y,-y),(-y,-y),(-zz,-zz),(bl,bh),(bl,bh),
                                  (-y,-y),(-zz,-zz),(-y,-y),(-zz,-zz)]
                    yy = ceil_div(max(abs(y-ceil_div(94_281*xh,UNIT)),abs(y-94_280*xl//UNIT)),10*n)
                    zbound = ceil_div(max(abs(2*zz+y-ceil_div(2*ROOT_HI*xh,UNIT)),abs(2*zz+y-2*ROOT_LO*xl//UNIT)),10*n)
                    xbound = ceil_div(xh, 10*n)
                    exponent = ceil_div(xbound*xbound+9*yy*yy+3*zbound*zbound,20_000)
                    weight = 6598*exp_floor(exponent)//10_000*jac//SCALE
                    if not weight:
                        plane.extend((0,)*dim)
                        continue
                    plane.append(weight)
                    for g, (lo, hi) in enumerate(coords):
                        plane.extend(values(shared[g], xl, xh, lo, hi))
            values.cache_clear()
            # Share the horizontal edge and square means between adjacent boxes.
            stride = n*dim
            lines = []
            for j in range(n+1):
                start = j*width
                left = plane[start:start+stride]
                right = plane[start+dim:start+width]
                lines.extend([isqrt(a*b) if a and b else 0 for a,b in zip(left,right)])
            squares = [isqrt(a*b) if a and b else 0 for a,b in zip(lines,lines[stride:])]
            if i:
                cube = [isqrt(a*b) if a and b else 0 for a,b in zip(old,squares)]
                for base in range(0,len(cube),dim):
                    w = cube[base]
                    if not w:
                        continue
                    terms = [sum(cube[base+first:base+last]) for first,last in zip(offsets,offsets[1:])]
                    if kind == 1:
                        counts[0] += 3*product(w, *terms[:3])
                        counts[1] += 3*product(w, *terms[3:])
                    else:
                        a, odd_y, b, regular, guarded_b, sy, wz, wy, sz = terms
                        counts[2] += 3*product(w, regular, a+product(odd_y,b))
                        counts[3] += 3*product(w, guarded_b, product(a,b)+product(sy,wz)+product(wy,sz))
            old = squares
    return counts, SCALE*n**3


if __name__ == "__main__":
    import sys
    print(cavity(int(sys.argv[1]) if len(sys.argv)>1 else 32))
