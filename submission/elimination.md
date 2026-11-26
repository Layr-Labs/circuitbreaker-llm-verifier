# Gate lower bound for an affine disperser in full B2

This file gives a self-contained proof of the parameter-method gate-elimination
lower bound of Mishra and Santhanam (ECCC TR26-260, October 19, 2026,
Sections 2–4, Theorem 4.1) in the full Li–Yang B2 model (ECCC TR21-023,
Sections 1.1, 2.1, 2.4, 3.3). An *affine disperser of dimension `d`* on `n`
bits is a Boolean function `f: F_2^n -> F_2` that is nonconstant on every
affine subspace of `F_2^n` of dimension `d` (and hence on every affine
subspace of dimension at least `d`).

**Lemma (Full-B2 Gate Elimination).** Let `f: F_2^n -> F_2` be an affine
disperser of dimension `d >= 0`. Then every full-B2 circuit `C` computing `f`
has at least `(7/2)*n - 7*d` gates.

## 1. Computational representation, invariants, and potential

Start with an arbitrary acyclic deterministic full-B2 circuit `C` computing
`f` on `n` input bits `x_1, ..., x_n`, with unrestricted sharing, depth, and
fan-out, arbitrary binary truth tables (all 16 operations, including constant
and unary gates), possible repeated pins, and constant nodes. Each gate of `C`
is charged once; inputs, constant nodes, wires, and output selection cost
zero. Allowing the final output to be an arbitrary literal (`w` or `w + 1`) of
a node only enlarges the class of circuits, so proving the bound for this
enlarged model establishes it for standard full-B2 circuits.

Over `F_2`, every nondegenerate binary truth table on `(u, v)` is either:

- **XOR-type (affine):** `u + v + c` for `c in {0, 1}`, or
- **AND-type (quadratic):** `(u + b)(v + c) + e` for `b, c, e in {0, 1}`.

The remaining 6 tables are constant (`0`, `1`) or unary literals (`u`, `u+1`,
`v`, `v+1`). Negations on internal wires can always be absorbed into the truth
tables of descendant binary gates without changing the gate count.

Throughout the proof we maintain:

1. A set of `M` independent Boolean parameters `u = (u_1, ..., u_M)` and a set
   `P` of `q` pairwise disjoint unordered pairs `{p, h}` of parameters. Put
   `a = M - q`.
2. An injective map `Phi: F_2^M -> F_2^n` of the form

       Phi(u) = b + L*u + sum_({p, h} in P) c_(p,h) * p * h,

   where `b, c_(p,h) in F_2^n` (possibly zero) and `L` is an `n`-by-`M` matrix
   over `F_2`.
3. A computation consisting of:
   - An *affine system* `A * v = N * u + c` over `F_2`, where `v` is a vector
     of `r_sys >= 0` auxiliary variables, `A` is an invertible `r_sys`-by-`r_sys`
     matrix over `F_2`, and each row (equation) contains at most three
     distinct nonconstant variable terms (parameters and/or auxiliaries).
   - An acyclic sequence of `r_gate >= 0` binary B2 gates whose inputs come
     from constants, parameters, auxiliaries, or earlier gates (no feedback
     into the affine system).
   - An *output reference* that is either a constant (`0` or `1`) or a literal
     (`w` or `w + 1`) of a parameter, auxiliary, or gate.
   The semantic invariant is that for every `u in F_2^M`, the unique solution
   `v = A^(-1)(N*u + c)` and gate evaluations yield output equal to `f(Phi(u))`.

Let `s = r_sys + r_gate` be the total number of equations plus acyclic gates.
For any parameter `z`, its *degree* `deg(z)` is the number of distinct objects
(rows of the affine system and acyclic gates) in which `z` appears; repeated
pins within one gate count once, and the output reference is not counted in
`deg(z)`. Similarly, for a gate `G`, `deg(G)` is the number of distinct
descendant gates reading `G`.

Define the parameter weight of an unpaired parameter `z` by

    weight(z) = 0  if deg(z) = 0,
    weight(z) = 6  if deg(z) = 1,
    weight(z) = 7  if deg(z) >= 2.

Paired parameters carry no individual weight. Call an unpaired parameter with
`deg(z) = 1` *isolated unpaired*. Define the nonnegative integer potential

    H = 2*s + 3*q + sum_(unpaired z) weight(z).

Initially `M = a = n`, `q = 0`, `Phi(u) = u` (injective), the affine system is
empty (`r_sys = 0`), `s = r_gate = |C|`, and every parameter has weight at
most `7`, so `H_0 <= 2*|C| + 7*n`. Always `H >= 0` and `M, s >= 0`.

**Affine-subspace invariant.** Fixing one chosen member of every pair
`{p, h} in P` to `0` eliminates all quadratic terms of `Phi`, yielding an
affine map on `M - q = a` variables. As a restriction of the injective map
`Phi`, this affine map is injective, so its image is an affine subspace of
`F_2^n` of dimension `a`. Because `f` is a dimension-`d` affine disperser,
`f(Phi(u))` is nonconstant whenever `a >= d`.

**Legal reparameterizations and restrictions.** The following four operations
preserve the displayed form and injectivity of `Phi`:

- *(Affine reparameterization)* For an unpaired parameter `x` and an affine
  function `ell` of the remaining parameters, substitute `x = w + ell` with a
  new unpaired parameter `w`. Because `x` is unpaired, `x` does not appear in
  any quadratic term of `Phi`, so `Phi` keeps its form; the substitution is a
  bijection on `F_2^M`, preserving injectivity and `(M, q, a)`.
- *(Quadratic reparameterization)* For three distinct unpaired parameters
  `x, y, z` and constants `b, c, e in F_2`, substitute `z = w + (x+b)(y+c) + e`
  with a new unpaired parameter `w` and add `{x, y}` to `P`. Because `z` is
  unpaired, no existing quadratic term uses `z`; because `x, y` were unpaired,
  `{x, y}` is disjoint from all existing pairs in `P`. The map is a bijection
  on `F_2^M`, preserving injectivity with `Delta M = 0`, `Delta q = 1`,
  `Delta a = -1`.
- *(Unpaired restriction)* Fix an unpaired parameter `x` to a constant or to
  `y + c` for another unpaired parameter `y`. Because `x` is unpaired, `Phi`
  keeps its form and remains injective on `F_2^(M-1)`, with `Delta M = -1`,
  `Delta q = 0`, `Delta a = -1`.
- *(Paired restriction)* For `{p, h} in P`, fix `p = c` or `p = h + c` and
  remove `p` and `{p, h}` from `P`. Over `F_2`, `h * c` and `h*(h+c) = (1+c)*h`
  are linear in `h`, so `Phi` keeps its form and remains injective on
  `F_2^(M-1)`, with `Delta M = -1`, `Delta q = -1`, `Delta a = 0`.

No reparameterization introduces uncharged gates or equations; each is used
only inside the exact transformations in Section 2.

**Terminal-output exclusions while `a > d`.** While `a > d`:

1. The output cannot be a constant or a literal of a parameter or auxiliary.
   Indeed, every parameter and auxiliary evaluates to an affine function
   `alpha_0 + sum_i alpha_i u_i` of `u`. A constant output contradicts the
   affine-subspace invariant (`a > d`). If some `alpha_i = 1`, then if `u_i`
   is paired with `h`, first apply the paired restriction `h = 0` (`Delta a = 0`)
   so `u_i` becomes unpaired with coefficient still `1`; then impose the
   affine restriction `u_i = alpha_0 + sum_(j != i) alpha_j u_j` (`Delta a = -1`),
   leaving an injective map with `a - 1 >= d` on which the output is constant
   `0`, contradicting the affine-subspace invariant.
2. The output cannot be a literal of an AND-type gate `G = (x+b)(y+c) + e` on
   two unpaired parameters `x, y`: applying the unpaired restriction `x = b`
   gives `a - 1 >= d` and a constant output, again contradicting the
   affine-subspace invariant.

## 2. Affine moves, legal deletions, and pivots

**Moving eligible affine gates.** Whenever an acyclic gate `G` has XOR-type
table `u_1 + u_2 + c` and both inputs are constants, parameters, or
auxiliaries (no gate ancestor), add a new auxiliary variable `v_G` and the
equation `v_G + u_1 + u_2 = c` (at most three nonconstant terms) to the affine
system, reroute all uses of `G` (including the output reference if applicable)
to `v_G`, and remove `G` from the acyclic part. The new auxiliary matrix is
block lower triangular with `1` in the new diagonal entry, hence invertible;
`s`, `M`, `q`, `a`, and all parameter degrees are unchanged. Move all eligible
affine gates in topological order.

**Four legal deletions** preserve `(M, q, a)`, decrease `s` by `1`, and satisfy
`Delta H <= -1`:

- **(D1) Degenerate or linearly dependent gate.** Remove any gate `G` whose
  output is constant or a single-parent literal `U + c`. This covers constant
  and unary truth tables, gates with a constant input or repeated pins, and
  AND-type gates whose two inputs are affine values (parameters/auxiliaries)
  with linearly dependent nonconstant parts (so one equals the other up to a
  constant, reducing the table to one variable). Reroute all uses of `G`
  (including the output reference) to that constant or `U + c`, absorbing `+c`
  into descendant gate tables (or the output literal). Deleting `G` gives
  `-2` in `H`. If rerouting to a constant, no weight rises; if rerouting to a
  parent `U`, `U` already had `deg(U) >= 1`, so `U` can gain weight only if
  `U` was isolated unpaired, and by at most `7 - 6 = 1`. Thus `Delta H <= -1`.
- **(D2) Unused gate.** Remove any acyclic gate `G` with `deg(G) = 0` that is
  not the output reference. Then `Delta s = -1`, no degree rises, and
  `Delta H <= -2`.
- **(D3) Short equation.** Any row of the affine system with at most two
  nonconstant terms must contain at least one auxiliary (since `A` is
  invertible). Write it as `v = U + c` where `v` is an auxiliary and `U` is
  either `0` or one other variable (parameter or auxiliary). Substitute
  `v = U + c` into every other row (cancelling pairs mod 2), every gate (absorbing
  `+c` into its truth table), and the output reference if it uses `v`; then
  delete the row and `v`. Adding the deleted row to every other row using `v`
  is an elementary row operation leaving column `v` with a single `1`, so
  cofactor expansion proves the remaining auxiliary matrix invertible. Each
  remaining row still has at most three nonconstant terms. Only `U` can gain
  occurrences; since `U` already appeared in the deleted row, its weight can
  rise by at most `1` (from `6` to `7`), against `-2` from deleting the row.
  Thus `Delta H <= -1`.
- **(D4) Isolated-parameter affine pivot (while `a > d`).** Suppose an
  isolated unpaired parameter `x` (`deg(x) = 1`) occurs only in a row `E` of
  the affine system. Then column `N_x` is the standard basis vector `e_E`.
  Because `A` is invertible, `A^(-1) N_x != 0`, so some auxiliary `Q` has
  coefficient `1` on `x` in `v = A^(-1)(N*u + c)`; write `Q = x + ell(other
  parameters)`. Apply the affine reparameterization `x = w + ell` to `Phi`, so
  `Q = w` identically. Swap the roles of auxiliary `Q` and parameter `x` in the
  system (moving column `Q` of `A` to the parameter side as `w` and column
  `N_x` to the auxiliary side as `x`). No equation terms change. Multiplying
  the new auxiliary matrix `A_new` on the left by `A^(-1)` yields the identity
  matrix with column `Q` replaced by `A^(-1) N_x`, whose determinant is the
  `Q`-th entry of `A^(-1) N_x`, namely `1`; thus `A_new` is invertible.
  Because `a > d`, neither `x` nor `Q` is the output reference (Section 1).
  Now auxiliary `x` appears only in row `E` and in no gate; delete `x` and
  row `E`. Cofactor expansion along column `x` preserves invertibility.
  Old parameter `x` (weight `6`) is removed, row `E` is removed (`-2`), new
  parameter `w` has weight at most `7`, and no other degree rises. Thus
  `Delta s = -1`, `Delta H <= -2 - 6 + 7 = -1`.

**Zero-cost affine pivot.** Under the same algebra as D4, if an auxiliary `Q`
depends on an unpaired parameter `x` of `deg(x) >= 2` (so `weight(x) = 7`),
swapping `Q` and `x` via `x = w + ell` preserves `(M, q, a, s)`, preserves
invertibility of `A`, and replaces weight `7` by `weight(w) <= 7`, so
`Delta H <= 0`. No equations or gates are added.

**Quadratic pivot.** Suppose an AND-type gate `G = (x+b)(y+c) + e` on distinct
unpaired parameters `x, y` has unique child `D = G + z + h` (XOR-type), where
`z` is an isolated unpaired parameter (`deg(z) = 1`) distinct from `x, y`. By
Section 1 (while `a > d`), `z` is not the output reference. Apply the
quadratic reparameterization `z = w + (x+b)(y+c) + e + h` to `Phi` and add
`{x, y}` to `P`. Then `D = w` identically; reroute all uses of `D` (including
the output reference if `D` is the output) to `w`, and delete `D` and the now
unused `G`. Because `z` occurred nowhere else, no other equation or gate is
affected. The changes satisfy

    Delta s = -2,   Delta M = 0,   Delta q = 1,   Delta a = -1,
    Delta H <= -4 + 3 + (7 - 6) - weight(x) - weight(y).

When `deg(x), deg(y) >= 2`, `weight(x) = weight(y) = 7`, so `Delta H <= -14`.

## 3. Preliminary reductions (in priority order)

Iterate while `a > d`. Before checking cases 1–6, move all eligible affine
gates into the affine system. After any case executes, restart at the top.
Every reduction below preserves the invariants, satisfies `Delta a in {0, -1, -2}`
and `Delta H <= 14 * Delta a`, never increases `M` or `s`, and strictly
decreases `M + s`.

1. **Legal deletions (D1–D4).** Execute any available D1–D4 deletion. Then
   `Delta a = 0`, `Delta s = -1`, `Delta H <= -1`.
   *Consequence when absent:* Every gate has two distinct nonconstant parents
   and a nondegenerate table; AND-type gates on affine parents have
   independent linear parts; every non-output gate has `deg >= 1`; every row
   has exactly three nonconstant terms; and no isolated unpaired parameter
   appears in any row of the affine system.

2. **AND gate on parameters with an isolated unpaired input.** Suppose a gate
   `G` has an isolated unpaired parameter `z` and another parameter `w` as
   parents. (Since `G` is not in the affine system, `G` is AND-type; by case 1,
   `w != z`.) Fix `w` to the forcing value making `G` constant, delete `G` by
   D1, and reroute the constant (which raises no weights; if `G` were the
   output, the constant output at `a_new >= d` would contradict Section 1).
   Now `z` becomes unused (`-6`).
   - If `w` is unpaired, `deg(w) >= 1` gives `weight(w) >= 6`: `Delta a = -1`,
     `Delta M = -1`, `Delta s = -1`, `Delta H <= -6 - 6 - 2 = -14`.
   - If `w` is paired, use the paired restriction to remove `{w, h}` (`-3`);
     its partner `h` becomes unpaired of weight at most `7`: `Delta a = 0`,
     `Delta M = -1`, `Delta s = -1`, `Delta H <= -3 - 6 - 2 + 7 = -4`.
   *Consequence when cases 1–2 are absent:* No object (row or gate) contains
   both an isolated unpaired parameter and any other parameter. (Rows are
   excluded by case 1; XOR gates on parameters were moved into the system;
   AND gates on parameters are excluded by case 2.)

3. **Paired parameter of degree 0 or `>= 2`.**
   - If paired `p` has `deg(p) = 0`, apply the paired restriction fixing its
     partner `h = 0` and removing `{p, h}`. Surviving `p` still has weight `0`:
     `Delta a = 0`, `Delta M = -1`, `Delta H = -3`.
   - If paired `p` has `deg(p) >= 2`, let `C_1, C_2` be two distinct children
     of `p`. Fix `p = 0` and remove `{p, h}` (`-3`, partner `h` gains weight
     at most `7`). Eliminate any rows among `C_1, C_2` by D3 and then any
     gates among `C_1, C_2` by D1 (in topological order). Each row among
     `C_1, C_2` has at most two terms after `p = 0` and cannot vanish because
     `A` stays invertible; each gate has a fixed pin `p = 0`. Thus both `C_1`
     and `C_2` are deleted (`Delta s <= -2`, `-4` in `H`). By cases 1–2, no
     isolated unpaired parameter co-occurs with `p` in `C_1` or `C_2`, so no
     parameter other than `h` can gain weight. Thus `Delta a = 0`,
     `Delta M = -1`, `Delta s <= -2`, `Delta H <= -3 - 4 + 7 = 0`.
   *Consequence when cases 1–3 are absent:* Every paired parameter has degree
   exactly `1`.

4. **Partners sharing their unique child.** Suppose `{p, h} in P` share their
   unique child `B`. Impose `p = h + c` (a paired restriction, `Delta a = 0`)
   choosing `c` so that `B` becomes constant (any `c` if `B` is a row
   `v + p + h = d` or XOR gate; `c = alpha + beta + 1` if `B` is
   `(p+alpha)(h+beta) + gamma`). Delete `B` by D3 (if a row `v = d + c`) or
   D1 (if a gate) and reroute the constant. Surviving `h` is now unused
   (weight `0`). Thus `Delta a = 0`, `Delta M = -1`, `Delta s = -1`,
   `Delta H = -3 - 2 = -5`.

5. **Paired parameter in a row.** Suppose paired `p` appears in a row `E_p`
   and its partner `h` has unique child `B_h`. By case 4, `B_h != E_p`, and by
   `deg(p) = 1`, `p` does not appear in `B_h`. Apply the paired restriction
   fixing `h = 0`, removing `{p, h}` (`-3`), and delete `B_h` by D3 (if a row)
   or D1 (if a gate), giving `-2`. By cases 1–2, no isolated unpaired
   parameter co-occurs with `h` in `B_h`, so no parameter other than `p` gains
   weight. Because `p` did not occur in `B_h` and `A` stays invertible, `E_p`
   remains in the system as the unique child of now-unpaired `p` (weight `6`,
   so `+6`). Immediately execute D4 on `(p, E_p)` (`Delta s = -1`,
   `Delta H <= -1`). Combined: `Delta a = 0`, `Delta M = -1`, `Delta s = -2`,
   `Delta H <= -3 - 2 + 6 - 1 = 0`.
   *Consequence when cases 1–5 are absent:* No paired parameter occurs in the
   affine system; every paired parameter occurs in a single acyclic gate.

6. **AND gate on parameters with a paired input.** Suppose an acyclic gate `G`
   reads paired `p` (partner `h`, unique child gate `B_h != G`) and another
   parameter `w`. By case 1, `w != p`; by case 4, `w != h` (since `G != B_h`
   and `deg(h) = 1`). Fix `h = 0` (removing `{p, h}`, `-3`) and delete `B_h`
   by D1 (`-2`). Next fix `w` to force `G` constant and delete `G` by D1
   (`-2`), rerouting the constant. Now `p` is unused (weight `0`). By cases
   1–2, `B_h` had no isolated unpaired co-parent, so no other weight rises.
   - If `w` is unpaired, case 2 ensures `deg(w) >= 2` (`weight(w) = 7`), so
     `Delta a = -1`, `Delta M = -2`, `Delta s = -2`,
     `Delta H <= -3 - 7 - 4 = -14`.
   - If `w` is paired (necessarily in a pair distinct from `{p, h}`), remove
     that pair (`-3`); its partner gains weight at most `7`, so `Delta a = 0`,
     `Delta M = -2`, `Delta s = -2`, `Delta H <= -6 - 4 + 7 = -3`.

Unused parameters stay in `M` until explicitly restricted; syntactic usage is
never conflated with affine dimension.

## 4. Preparing a minimal AND gate

Assume `a > d` and none of cases 1–6 applies. By Section 1, the output is not
a constant, parameter, or auxiliary, so the acyclic part is nonempty. Choose a
first (minimal) gate `G` in topological order. Because all eligible affine
gates were moved, `G` is AND-type and both parents of `G` are affine values
(parameters or auxiliaries) with independent linear parts (case 1).

1. **Pivot auxiliary parents to parameters.** If the first parent `Q_1` of `G`
   is an auxiliary, its linear part is nonzero, so `Q_1` depends on some
   parameter `x_1` appearing in the affine system (`N_(x_1) != 0`). Because
   cases 1 and 3–5 are absent, `x_1` is unpaired with `deg(x_1) >= 2`. Perform
   the zero-cost affine pivot of Section 2 swapping `Q_1` and `x_1` (`Delta H
   <= 0`), move eligible affine gates, and recheck cases 1–6 (restarting if
   any applies). If the second parent `Q_2` of `G` is still an auxiliary,
   linear independence of `G`'s two parents implies `Q_2` depends on some
   parameter `x_2` distinct from `G`'s first parameter parent; again `x_2` is
   unpaired with `deg(x_2) >= 2`. Pivot `Q_2` with `x_2` (`Delta H <= 0`),
   move eligible affine gates, and recheck cases 1–6. After at most two pivots,
   either a preliminary reduction executed or `G` reads two distinct unpaired
   parameters `x, y` with `deg(x), deg(y) >= 2` (by cases 1, 2, 6).
2. **Eliminate children of `G` that also read `x` or `y`.** While `G` has a
   child gate `J` whose other input is `x` (respectively `y`), replace
   `J = h(x, G(x, y))` by the single B2 gate `h'(x, y) = h(x, G(x, y))` with
   inputs rerouted to `(x, y)`. This keeps `s, M, q, a` unchanged; only `y`
   (respectively `x`) gains an occurrence, and its weight is already `7`, so
   `Delta H = 0`. Move any newly eligible affine gates (which cannot affect
   `G`'s parents `x, y` or add children to `G`) and recheck cases 1–6. If a
   preliminary case executes, restart the main loop. Otherwise `G` remains an
   AND gate on `x, y` whose degree `deg(G)` strictly decreased by `1`, so this
   step terminates in at most `deg(G) <= s` iterations.

If preparation completes without triggering cases 1–6, then `G` is an AND gate
on distinct unpaired `x, y` with `deg(x) >= deg(y) >= 2` (after relabelling),
no child of `G` reads `x` or `y`, and by Section 1 and case 1, `G` is not the
output and has `h = deg(G) >= 1`. Put `u = deg(x) >= v = deg(y) >= 2`.

## 5. Exhaustive remaining cases

7. **Shared second child at `u = v = 2`.** Suppose `u = v = 2` and `x, y` share
   their second child `B != G`. Because `G` is AND-type, `G(x, y)` is constant
   on three of the four assignments `(x, y) in F_2^2`. The second function `B`
   (either a binary gate on `(x, y)` or a row `v + x + y = c`, whose linear
   part `x + y` is a binary function of `(x, y)`) takes the same value on at
   least two of those three points. Any two distinct points in `F_2^2` form an
   affine line defined by fixing `x`, fixing `y`, or setting `x = y + c`.
   Apply that unpaired restriction (`Delta a = -1`, `Delta M = -1`). Both `G`
   and `B` become constant (or `B` becomes `v = const`, eliminated by D3);
   delete both (`Delta s = -2`, `-4`) and reroute constants. The surviving
   parameter is now unused (`-7 - 7 = -14`), and constant rerouting raises no
   weights. Thus `Delta a = -1`, `Delta H <= -7 - 7 - 4 = -18 < -14`.

Henceforth, whenever `u = v = 2`, the other children `B_x` of `x` and `B_y` of
`y` are distinct.

8. **Forcing `x` to make `G` constant (all cases except `u = v = 2, h = h_1 = 1`).**
   Let `h_1` be the number of children `D_1, ..., D_h` of `G` whose other
   parent is an isolated unpaired parameter (`0 <= h_1 <= h`). Suppose
   `(u, h, h_1) != (2, 1, 1)`. Fix `x` to the forcing value `b` making `G`
   constant `e` (`Delta a = -1`, `Delta M = -1`, `x` loses weight `7`).
   Let `G, B_1, ..., B_(u-1)` be the `u` children of `x`. Because the circuit
   is acyclic and no child of `G` reads `x`, the `u + h` objects
   `G, B_1, ..., B_(u-1), D_1, ..., D_h` are pairwise distinct.
   Delete them in three stages:
   - Eliminate all rows among `B_1, ..., B_(u-1)` by D3 one by one. Each has
     at most two terms after `x = b` and stays non-vanishing because `A` stays
     invertible. By case 1, no isolated unpaired parameter occurs in any row,
     so no isolated unpaired parameter gains weight here.
   - Delete `G` by D1 and reroute all uses of `G` to constant `e`.
   - Delete all gates among `B_1, ..., B_(u-1), D_1, ..., D_h` in topological
     order by D1 (each has a fixed pin `x = b` or `G = e`). By cases 1–2, no
     gate `B_i` originally read an isolated unpaired parameter, and only the
     `h_1` counted children of `G` originally read an isolated unpaired
     parameter. By induction along topological order, every surviving parent
     created during rerouting that is an isolated unpaired parameter must be
     one of those at most `h_1` co-parents. Each can gain weight at most
     `7 - 6 = 1`, for a total weight increase of at most `h_1`.
   Thus `Delta s <= -(u + h)`, giving

       Delta H <= -7 - 2*(u + h) + h_1.

   - If `u >= 3`: since `h_1 <= h`, `Delta H <= -7 - 6 - h <= -14`.
   - If `u = 2`: then `v = 2`, `B_x != B_y`, and no child of `G` reads `y`.
     Consequently `y` does not appear in `B_x` or any `D_j`, so `y` gains no
     new occurrence during the deletions and loses `G`. Thus `deg(y)` drops
     from `2` to at most `1`, lowering `weight(y)` by at least `1` (`7 -> <= 6`),
     so

         Delta H <= -7 - 1 - 2*(2 + h) + h_1 = -12 - 2*h + h_1.

     Since `(h, h_1) != (1, 1)` and `h >= 1`, either `h >= 2` (so `-2*h + h_1
     <= -h <= -2`) or `(h, h_1) = (1, 0)` (so `-2*h + h_1 = -2`). In both
     cases `Delta H <= -14`.

In the sole remaining configuration, `u = v = 2` and `h = h_1 = 1`: `G` has a
unique child gate `D` whose other parent is an isolated unpaired parameter `z`
(distinct from `x, y`), and `B_x != B_y`. By case 1, `D` is nondegenerate,
hence either AND-type (case 9) or XOR-type (case 10).

9. **`u = v = 2, h = h_1 = 1` with `D` AND-type.** Write `D = (G+c)(z+b) + e`.
   Because `G` is nondegenerate, choose `(x_0, y_0) in F_2^2` with
   `G(x_0, y_0) = c`, and apply the two unpaired restrictions `x = x_0`,
   `y = y_0` (`Delta a = -2`, `Delta M = -2`, `-7 - 7 = -14`). The four
   objects `G, D, B_x, B_y` are pairwise distinct. Eliminate any rows among
   `B_x, B_y` by D3; delete `G` (constant `c`) and `D` (constant `e`) by D1;
   and delete any gates among `B_x, B_y` by D1. Rerouting constants `G, D`
   raises no weights; `z` loses its unique child `D` and becomes unused (`-6`);
   by cases 1–2, neither `B_x` nor `B_y` had an isolated unpaired co-parent,
   so no other weight rises. Thus `Delta s <= -4` and

       Delta a = -2,   Delta H <= -7 - 7 - 6 - 8 = -28 = 14 * Delta a.

   (If `D` was the output, the output becomes constant `e`; if `a - 2 >= d`
   this contradicts Section 1, while `a - 2 = d - 1 < d` terminates the loop
   with `H_T >= 0`.)

10. **`u = v = 2, h = h_1 = 1` with `D` XOR-type.** Execute the quadratic
    pivot of Section 2 on `G, D, x, y, z`. Because `weight(x) = weight(y) = 7`,

        Delta a = -1,   Delta s = -2,   Delta M = 0,
        Delta H <= -4 + 3 + (7 - 6) - 7 - 7 = -14.

Cases 1–10 are exhaustive whenever `a > d`.

## 6. Termination and final gate count

Each iteration while `a > d` performs a finite preparation (topological affine
moves, at most two zero-cost pivots, and at most `s` child reroutings) and
then executes one of cases 1–10, which preserves all invariants, never
increases `M` or `s`, and strictly decreases the nonnegative integer `M + s`.
Therefore the process terminates after at most `n + |C|` iterations in a valid
state with `a_T <= d` and `H_T >= 0`.

Every iteration satisfies `Delta H <= 14 * Delta a` (equivalently
`H_t - H_(t+1) >= 14 * (a_t - a_(t+1))`). Summing over all iterations gives

    2*|C| + 7*n >= H_0 >= H_T + 14*(n - a_T) >= 14*(n - d).

Dividing by `2` yields `|C| >= (7/2)*n - 7*d`. (If `n <= d` initially, the
bound `(7/2)*n - 7*d <= -(7/2)*n <= 0` holds trivially because `|C| >= 0`.)
This proves the Lemma.
