"""Lower sums for the four local cubic saving events."""

from bisect import bisect_right
from fractions import Fraction
from math import isqrt

from cavity_arithmetic import (DISTANCES, G_LO, GRID, LAST, SCALE, SPACE, exp_lower,
                               geometric_interpolation, gm4, interp_row, radial, root2)


TIME_PANELS = 16
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


def _table_values(family, parent_index, interval):
    if interval[0] < 0 or interval[1] > 2 * SPACE:
        return (0,) * len(family)
    lo, hi = interval
    j = bisect_right(GRID, lo) - 1
    if j == LAST or lo == hi == GRID[j]:
        return tuple(table[parent_index][j] for table in family)
    if hi <= GRID[j+1]:
        l_delta, h_delta = lo - GRID[j], hi - GRID[j]
        width = DISTANCES[j]
        values = []
        for table in family:
            a, b = table[parent_index][j:j+2]
            values.append(geometric_interpolation(a, b, l_delta if a <= b else h_delta, width))
        return tuple(values)
    return tuple(interp_row(table, parent_index, lo, hi) for table in family)


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


def _kind_one_slice(i, roots, b_panels, positive, switched):
    """One x slice of two-endpoint geometric means along the time axis."""
    x = GRID[i]
    rx_lo, rx_hi = roots[i]
    slices = []
    zero = ((0,) * len(positive),) * 2 + ((0,) * len(switched),) * 2
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
                _table_values(positive, i, y),
                _table_values(positive, i, z),
                _table_values(switched, i, y),
                _table_values(switched, i, z),
            ) if density else zero
            if endpoint is not None:
                prior_density, prior_values = endpoint
                time_edges.append((isqrt(density * prior_density),
                                   tuple(_edge_means(a, b) for a, b in zip(values, prior_values, strict=True))))
            endpoint = density, values
        slices.append(time_edges)
    return slices


def _kind_one(A, P, even, extended, roots):
    x_panels = _last_node(23, 32)
    b_panels = _last_node(5, 4)
    pos = neg = 0
    left = _kind_one_slice(0, roots, b_panels, A, P)
    scale_cube = SCALE**3
    for i in range(x_panels):
        right = _kind_one_slice(i+1, roots, b_panels, A, P)
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
    scaled_rho_x = (2*rx_hi + 2)//3
    slices = []
    zero = ((0,) * len(positive), (0,) * len(switched))
    for y in GRID[:yz_panels + 1]:
        factor_xy = x*x + 9 * (y + scaled_rho_x)**2
        endpoint = None
        z_edges = []
        for z in GRID[:yz_panels + 1]:
            companion = (2*rx_lo + y + z, 2*rx_hi + y + z)
            exponent = factor_xy + 3 * (2*z + y + 2*rx_hi)**2
            density = TRIPLE_AT_ZERO * exp_lower(exponent, 2 * SPACE_SQUARED) // SCALE
            values = (
                _table_values(positive, i, companion),
                _table_values(switched, i, companion),
            ) if density else zero
            if endpoint is not None:
                prior_density, prior_values = endpoint
                z_edges.append((isqrt(density * prior_density),
                                tuple(_edge_means(a, b) for a, b in zip(values, prior_values, strict=True))))
            endpoint = density, values
        slices.append(z_edges)
    return slices


def _kind_two(A0, P, even, odd, guard, extended, roots):
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
            for k, corners in enumerate(zip(left[j], left[j+1], right[j], right[j+1], strict=True)):
                density = gm4(*(point[0] for point in corners))
                if not density:
                    continue
                dz = GRID[k+1] - GRID[k]
                first, second, third, fourth = (point[1] for point in corners)
                positive = _sum_corner_means(first[0], second[0], third[0], fourth[0])
                switched = _sum_corner_means(first[1], second[1], third[1], fourth[1])
                even_y, even_z = even[i][j], even[i][k]
                pos_bracket = even_y + odd[i][j] * even_z // SCALE
                neg_bracket = (even_y * even_z + extended[i][j] * guard[i][k]
                               + guard[i][j] * extended[i][k]) // SCALE
                vol = dx * dy * dz
                pos += vol * (density * positive * pos_bracket // scale_squared)
                neg += vol * (density * switched * neg_bracket // scale_squared)
        left = right

    denominator = SCALE * SPACE**3
    return Fraction(3 * pos, denominator), Fraction(3 * neg, denominator)


def integrate_roots(A, A0, P, B_even, B_odd, O_guard, E_cells):
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
    return (*_kind_one(A, P, even, E_cells, roots),
            *_kind_two(A0, P, even, odd, guard, E_cells, roots))
