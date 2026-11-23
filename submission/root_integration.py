"""Lower sums for the four local cubic saving events."""

from bisect import bisect_right
from fractions import Fraction
from functools import lru_cache
from math import isqrt
from struct import Struct

from cavity_arithmetic import (DISTANCES, G_LO, GRID, LAST, SCALE, SPACE, exp_lower,
                               gm4, radial, root2)


TIME_PANELS = 24
PAIR_AT_ZERO = SCALE * 47_746_482 // 10**8
TRIPLE_AT_ZERO = SCALE * 65_984_522 // 10**8
SPACE_SQUARED = SPACE * SPACE


def _verify_density_constants():
    """Bound pi above by a short rational Machin enclosure."""

    def atan_bounds(inverse, pairs):
        partial = sum(Fraction((-1)**k, (2*k + 1) * inverse**(2*k + 1))
                      for k in range(2 * pairs))
        return partial, partial + Fraction(1, (4*pairs + 1) * inverse**(4*pairs + 1))

    _, upper_fifth = atan_bounds(5, 6)
    lower_239, _ = atan_bounds(239, 2)
    pi_upper = 16 * upper_fifth - 4 * lower_239
    assert Fraction(G_LO, SCALE)**2 * pi_upper < 24
    assert Fraction(PAIR_AT_ZERO, SCALE) * 2 * pi_upper < 3
    assert Fraction(TRIPLE_AT_ZERO, SCALE)**2 * (2*pi_upper)**3 < 108


def _last_node(numerator, denominator):
    coordinate, remainder = divmod(SPACE * numerator, denominator)
    assert remainder == 0
    return GRID.index(coordinate)


def _multiply_interval(interval, multiplier, denominator):
    low, high = interval
    return low * multiplier // denominator, (high * multiplier + denominator - 1) // denominator


def _row_lookup(family, parent_index):
    """Evaluate directed lower tangents for all itineraries in one row."""
    columns = list(zip(*(table[parent_index] for table in family)))
    size = len(family)
    zeros = (0,) * size
    layout = Struct(f"<{size}Q")
    subdivisions = 16
    bits = 16
    unit = 1 << bits
    mask = int.from_bytes(layout.pack(*([(1 << (64-bits))-1] * size)), "little")
    assert max(map(max, columns)) <= SCALE and 3*SCALE*unit < 1 << 64

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
            small, large = min(a, b), max(a, b)
            slope = 2*small*(large-small)//(large+small) if small else 0
            low.append(small)
            up.append(slope if a < b else 0)
            down.append(slope if a > b else 0)
        return pack(low) << bits, pack(up), pack(down)

    def evaluate(lo, hi):
        j = bisect_right(GRID, lo) - 1
        if j == LAST:
            return columns[j]
        width = DISTANCES[j] // subdivisions
        k, offset = divmod(lo-GRID[j], width)
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

    def lookup(interval):
        lo, hi = interval
        return evaluate(lo, hi) if 0 <= lo <= hi <= GRID[-1] else zeros
    return lookup


def _edge_means(left, right):
    return tuple(isqrt(a*b) for a, b in zip(left, right, strict=True))


def _sum_corner_means(first, second, third, fourth):
    return sum(gm4(a, b, c, d) for a, b, c, d in zip(first, second, third, fourth, strict=True)
               if a and b and c and d)


def _surface_means(family, n_rows, n_cols):
    return [
        [
            sum(gm4(table[i][j], table[i+1][j], table[i][j+1], table[i+1][j+1])
                for table in family)
            for j in range(n_cols)
        ]
        for i in range(n_rows)
    ]


def _success_bands(family, n_rows, n_cols):
    """Bounds on the sum of the log-interpolated arm probabilities in each cell."""
    bands = []
    for i in range(n_rows):
        row = []
        for j in range(n_cols):
            corners = [(table[i][j], table[i+1][j], table[i][j+1], table[i+1][j+1])
                       for table in family]
            lower = sum(map(min, corners))
            upper = max(map(sum, zip(*corners)))
            assert 0 <= lower <= upper <= SCALE
            row.append((lower, upper))
        bands.append(row)
    return bands


def _kind_one_slice(i, roots, b_panels, positive, switched, negative_middle):
    """One x slice of two-endpoint geometric means along the time axis."""
    x = GRID[i]
    rx_lo, rx_hi = roots[i]
    slices = []
    negative_small = (*switched, negative_middle)
    zero = ((0,) * len(positive),) * 2 + ((0,) * len(switched), (0,) * len(negative_small))
    pos_at = _row_lookup(positive, i)
    neg_at = _row_lookup(negative_small, i)
    for b in GRID[:b_panels + 1]:
        half_sum = (rx_lo + b//2, rx_hi + (b+1)//2)
        pair_exp = 9 * (x*x + b*b) + 12 * rx_hi * b
        pair_density = PAIR_AT_ZERO * exp_lower(pair_exp, 2 * SPACE_SQUARED)
        endpoint = None
        time_edges = []
        for t in range(TIME_PANELS + 1):
            y = _multiply_interval(half_sum, TIME_PANELS + t, TIME_PANELS)
            z = _multiply_interval(half_sum, TIME_PANELS - t, TIME_PANELS)
            h = _multiply_interval(half_sum, t, TIME_PANELS)
            density = pair_density * radial(*h) * half_sum[0] // (SCALE*SCALE*SPACE)
            values = (
                pos_at(y), pos_at(z), neg_at(y)[:-1], neg_at(z),
            ) if density else zero
            if endpoint is not None:
                prior_density, prior_values = endpoint
                time_edges.append((isqrt(density * prior_density),
                                   tuple(_edge_means(a, b) for a, b in zip(values, prior_values, strict=True))))
            endpoint = density, values
        slices.append(time_edges)
    return slices


def _kind_one(A, P, negative_middle, even, extended, roots):
    x_panels = _last_node(23, 32)
    b_panels = _last_node(5, 4)
    pos = neg = 0
    left = _kind_one_slice(0, roots, b_panels, A, P, negative_middle)
    scale_cube = SCALE**3
    for i in range(x_panels):
        right = _kind_one_slice(i+1, roots, b_panels, A, P, negative_middle)
        dx = GRID[i+1] - GRID[i]
        for j in range(b_panels):
            mass_pos, mass_neg = even[i][j], extended[i][j]
            if not (mass_pos or mass_neg):
                continue
            db = GRID[j+1] - GRID[j]
            for corners in zip(left[j], left[j+1], right[j], right[j+1], strict=True):
                density = gm4(*(point[0] for point in corners))
                if not density:
                    continue
                first, second, third, fourth = (point[1] for point in corners)
                if mass_pos:
                    upper = _sum_corner_means(first[0], second[0], third[0], fourth[0])
                    lower = _sum_corner_means(first[1], second[1], third[1], fourth[1])
                    pos += dx * db * (density * mass_pos * upper * lower // scale_cube)
                if mass_neg:
                    upper = _sum_corner_means(first[2], second[2], third[2], fourth[2])
                    lower = _sum_corner_means(first[3], second[3], third[3], fourth[3])
                    neg += dx * db * (density * mass_neg * upper * lower // scale_cube)
        left = right

    denominator = SCALE * SPACE_SQUARED * TIME_PANELS
    return Fraction(3 * pos, denominator), Fraction(3 * neg, denominator)


def _kind_two_slice(i, roots, yz_panels, positive, switched):
    """One x slice of density and itinerary edge means along the z axis."""
    x = GRID[i]
    rx_lo, rx_hi = roots[i]
    slices = []
    zero = ((0,) * len(positive), (0,) * len(switched))
    pos_at = _row_lookup(positive, i)
    neg_at = _row_lookup(switched, i)
    for j, y in enumerate(GRID[:yz_panels + 1]):
        endpoint = None
        z_edges = []
        for z in GRID[:min(j+1, yz_panels) + 1]:
            companion = (2*rx_lo + y + z, 2*rx_hi + y + z)
            exponent = 33*x*x + 12*(y*y + y*z + z*z) + 24*rx_hi*(y+z)
            density = TRIPLE_AT_ZERO * exp_lower(exponent, 2 * SPACE_SQUARED) // SCALE
            values = (
                pos_at(companion), neg_at(companion),
            ) if density else zero
            if endpoint is not None:
                prior_density, prior_values = endpoint
                z_edges.append((isqrt(density * prior_density),
                                tuple(_edge_means(a, b) for a, b in zip(values, prior_values, strict=True))))
            endpoint = density, values
        slices.append(z_edges)
    return slices


def _kind_two(A0, P, even, odd, bands, guard, extended, roots):
    x_panels = _last_node(17, 32)
    yz_panels = _last_node(17, 16)
    pos = neg = 0
    left = _kind_two_slice(0, roots, yz_panels, A0, P)
    scale_squared = SCALE**2
    for i in range(x_panels):
        right = _kind_two_slice(i+1, roots, yz_panels, A0, P)
        dx = GRID[i+1] - GRID[i]
        for j in range(yz_panels):
            dy = GRID[j+1] - GRID[j]
            for k in range(j+1):
                corners = (left[j][k], left[j+1][k], right[j][k], right[j+1][k])
                density = gm4(*(point[0] for point in corners))
                if not density:
                    continue
                dz = GRID[k+1] - GRID[k]
                first, second, third, fourth = (point[1] for point in corners)
                positive = _sum_corner_means(first[0], second[0], third[0], fourth[0])
                switched = _sum_corner_means(first[1], second[1], third[1], fourth[1])
                even_y, even_z = even[i][j], even[i][k]
                low_y, high_y = bands[i][j]
                low_z, high_z = bands[i][k]
                # Each expression bounds the same union integral on this cell.
                pos_bracket = max(
                    even_y + odd[i][j] * even_z // SCALE,
                    even_z + odd[i][k] * even_y // SCALE,
                    ((SCALE-low_z)*even_y + (SCALE-high_y)*even_z + high_y*low_z) // SCALE,
                    ((SCALE-low_y)*even_z + (SCALE-high_z)*even_y + high_z*low_y) // SCALE,
                )
                neg_bracket = (even_y * even_z + extended[i][j] * guard[i][k]
                               + guard[i][j] * extended[i][k]) // SCALE
                vol = dx * dy * dz * (1 if j == k else 2)
                pos += vol * (density * positive * pos_bracket // scale_squared)
                neg += vol * (density * switched * neg_bracket // scale_squared)
        left = right

    denominator = SCALE * SPACE**3
    return Fraction(3 * pos, denominator), Fraction(3 * neg, denominator)


def integrate_roots(A, A0, P, negative_middle, B_even, B_odd, O_guard, E_cells):
    """Return exact lower fractions for (pos1, neg1, pos2, neg2).

    Point tables contain one matrix per itinerary. E_cells instead contains
    lower sums of four-corner means for the extended E itineraries.
    """
    _verify_density_constants()
    roots = tuple(root2(coordinate, coordinate) for coordinate in GRID)
    n_rows = _last_node(23, 32)
    n_cols = _last_node(5, 4)
    even = _surface_means(B_even, n_rows, n_cols)
    odd = _surface_means(B_odd, n_rows, n_cols)
    guard = _surface_means(O_guard, n_rows, n_cols)
    bands = _success_bands(B_even, n_rows, n_cols)
    return (*_kind_one(A, P, negative_middle, even, E_cells, roots),
            *_kind_two(A0, P, even, odd, bands, guard, E_cells, roots))
