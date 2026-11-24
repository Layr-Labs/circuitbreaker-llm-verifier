"""Directed two-arm union bounds and batched log-concave interpolation."""
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
