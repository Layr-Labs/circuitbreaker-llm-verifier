# The Ben-Sasson–Kopparty disperser bound for the certificate-prefix family

This file proves the disperser lemma used in `proof.md`, following the
construction and proof strategy of Ben-Sasson–Kopparty (SIAM J. Comput.
41(4), 2012; author preprint 2010-03-12, Sections 3, 6, 7), restated and
re-derived here for characteristic `2` and prime field degree `k`, which is
the case the family of `proof.md` uses. Prime degree is what makes the
subfield side conditions in the general argument vacuous; this is noted at
each use.

## D.1 Setup and basic facts

Let `k` be prime, `K = F_{2^k}`, and identify `K` with
`F_2[Z]/(Q)` for a monic irreducible `Q` of degree `k` when a basis is
needed. For a polynomial `h(X_1, ..., X_r)` over `K`, write `h mod (X_i^{2^k}
- X_i)` for the unique polynomial of degree `< 2^k` in each variable
obtained by replacing `X_i^{2^k}` by `X_i` repeatedly.

**Fact 1 (values and reduced forms).** Reduction does not change values on
`K^r`, and a function `K^r -> K` has a *unique* representing polynomial of
degree `< 2^k` in each variable (a nonzero such polynomial cannot vanish
on all of `K^r`: for `r = 1` it has more roots than its degree; induct).
Hence a reduced polynomial represents a nonconstant function iff it has a
monomial of positive degree with nonzero coefficient.

**Fact 2 (trace form).** Every `F_2`-linear map `pi : K -> F_2` is
`pi(z) = Tr(mu z)` for a unique `mu in K`, where
`Tr(z) = sum_{rho=0}^{k-1} z^{2^rho}`; `pi` is nontrivial iff `mu != 0`.
In particular this holds for the constant-coefficient functional of
`proof.md`. Consequently, for a subspace `Y <= K`, writing
`Y^\perp := {x in K : Tr(xy) = 0 for all y in Y}` for the common kernel
of the functionals `z |-> Tr(mu z)` with `mu in Y`, the map
`x |-> (y |-> Tr(xy))` from `K` to the `F_2`-dual of `Y` is surjective
(every `F_2`-functional on `Y` extends to `K`, where it has this form),
so `dim Y^\perp = k - dim Y` and `Y^{\perp\perp} = Y`; orthogonal
complementation is thus a bijection between the subspaces of `K` of
dimension `d` and those of dimension `k - d`, for every `d`.

A polynomial `P(X) = sum_j a_j X^{2^j}` (`a_j in K`) is *linearized*; it
computes an `F_2`-linear map `K -> K`. The `F_2`-*degree* of a monomial
`prod X_i^{beta_i}` is `sum_i wt_2(beta_i)`, where `wt_2` is the binary
digit sum; linearized polynomials are exactly the polynomials of
`F_2`-degree `<= 1` with zero constant term. Evaluation gives a bijection
between linearized polynomials `sum_{j<k} a_j X^{2^j}` and `F_2`-linear
maps `K -> K` (there are `2^{k^2}` of each, and a nonzero one has at most
`2^{k-1}` roots, so evaluation is injective).

For an `F_2`-linear subspace `Z <= K` of dimension `m`, the
*kernel-subspace polynomial* of `Z` is
`P_Z(X) = prod_{z in Z} (X - z)`.

**Claim (subspace polynomials).** `P_Z` is a monic linearized polynomial
of degree `2^m` whose root set is exactly `Z`; as a linear map it has
`ker P_Z = Z` and `dim im(P_Z) = k - m`. It is the unique monic
linearized polynomial of degree `2^m` vanishing on `Z`.

*Proof.* Induction on `m`. For `m = 0`, `P_{ {0} } = X`. Otherwise write
`Z = Z' \oplus {0, z_0}` with `dim Z' = m - 1` and `z_0 notin Z'`.
Splitting the product over the cosets `Z'` and `z_0 + Z'` gives
`P_Z(X) = P_{Z'}(X) \cdot P_{Z'}(X + z_0)` (in characteristic `2`,
`X - z = X + z`), and by induction `P_{Z'}` is additive, so
`P_{Z'}(X + z_0) = P_{Z'}(X) + P_{Z'}(z_0)`; with
`c := P_{Z'}(z_0)` this gives
`P_Z(X) = P_{Z'}(X)^2 + c \cdot P_{Z'}(X)`.
Squaring a linearized polynomial gives a linearized polynomial
(`(sum_j b_j X^{2^j})^2 = sum_j b_j^2 X^{2^{j+1}}`), as do scalar
multiples and sums, so `P_Z` is linearized; the square has degree `2^m`
and leading coefficient `1` while `c \cdot P_{Z'}` has degree at most
`2^{m-1}`, so `P_Z` is monic of degree exactly `2^m`. Its root set is
the union of the factors' root sets, namely `Z`. For uniqueness, if `L`
is monic linearized of degree `2^m` and vanishes on `Z`, then
`L - P_Z` has degree `< 2^m` (the leading terms cancel) and `2^m` roots,
so `L = P_Z` by Fact 1. ∎

**Lemma D0 (image-subspace polynomials).** For every `F_2`-linear
subspace `Y <= K` of dimension `d`, `1 <= d <= k`, there is a unique
monic linearized polynomial of degree exactly `2^{k-d}` whose image is
exactly `Y`.

*Proof.* For `Z <= K` of dimension `k - d`, let `Phi(Z) = im(P_Z)`,
which has dimension `d` by the Claim. We show `Phi` is injective. Let
`Z, Z'` have dimension `k - d` with `Phi(Z) = Phi(Z') =: Y'`, and put
`A = P_{Y'}`. Then `A \circ P_Z = X^{2^k} - X`: the composition is monic
of degree `2^d \cdot 2^{k-d} = 2^k` (leading coefficient
`1 \cdot 1^{2^d} = 1`), and it vanishes at every `x in K` because
`P_Z(x) in Y' = ker A`; `X^{2^k} - X` is likewise monic of degree `2^k`
and vanishes on all of `K`; the difference has degree `< 2^k` and `2^k`
roots, hence is zero by Fact 1. The same argument gives
`A \circ P_{Z'} = X^{2^k} - X`. Composition with a linearized polynomial
distributes over subtraction of the inner polynomial: writing
`A = sum_j b_j X^{2^j}`,
`A \circ (U - V) = sum_j b_j (U - V)^{2^j} = sum_j b_j (U^{2^j} -
V^{2^j}) = A \circ U - A \circ V`, so `A \circ (P_Z - P_{Z'}) = 0`. But
for `B := P_Z - P_{Z'} != 0`, the composition `A \circ B` is nonzero:
the term `B^{2^d}` has degree `2^d \cdot deg B` and every other term
`b_j B^{2^j}`, `j < d`, has strictly smaller degree, so
`deg (A \circ B) = 2^d \cdot deg B`. Hence `P_Z = P_{Z'}` and
`Z = ker P_Z = ker P_{Z'} = Z'`: `Phi` is injective.

By Fact 2, orthogonal complementation is a bijection between the
subspaces of dimension `d` and those of dimension `k - d`, so the two
Grassmannians have the same finite cardinality. An injective map
between finite sets of the same cardinality is surjective, so every
`d`-dimensional `Y` equals `im(P_Z)` for some `Z`, and `P_Z` is the
required polynomial. For uniqueness, if `L` is monic linearized of
degree `2^{k-d}` with `im L = Y`, the Claim applied to `ker L`
(dimension `k - d`, since `im L` has dimension `d`) gives
`L = P_{ker L}`, and injectivity of `Phi` gives `ker L = Z`, so
`L = P_Z`. ∎

For an `F_2`-linear subspace `Y <= K` of dimension `d`,
`1 <= d <= k`, the *image-subspace polynomial* `Q_Y` is the polynomial
of Lemma D0: the unique monic linearized polynomial with `Q_Y(K) = Y`;
it has degree `2^{k-d}`. (This is Lemma 3.7 of BK; the proof above is
self-contained and cites no external result.) Write
`Q_Y(X) = sum_{j=0}^{k-1} a_j X^{2^j}` (reduced form; `a_j = 0` for
`j > k - d` and `a_{k-d} = 1`). For completeness, set
`Q_{ {0} } := X^{2^k} - X`, whose associated function has image `{0}`;
this case is used only in the parametrisation of Theorem D6 and never
in Lemma D1 or Theorem D3.

**Lemma D1 (structural lemma, prime degree).** With `Q_Y` as above:
(1) every `d` cyclically consecutive coefficients `a_j` contain a nonzero
one; in particular `a_0 != 0` and `a_{k-d} != 0` (so the zero coefficients
above `a_{k-d}` form a run of length exactly `d - 1`, namely
`a_{k-d+1} = ... = a_{k-1} = 0`).
(2) every zero run of the cyclic sequence `(a_j)` other than that top run
has length at most `d - 2` (for `d >= 2`).

*Proof.* This is BK Lemma 3.10 specialised. We use the Frobenius shift
of coefficients (BK Claim 3.9): the reduced form of `Q_Y(X)^{2^s}` has
coefficient `(a_{j-s})^{2^s}` at position `j` (indices mod `k`); since
`a^{2^s} = 0` iff `a = 0`, only the zero pattern matters below. For (1):
suppose `a_j = 0` for all `j` in a cyclic window `{r, ..., r + d - 1}`.
With `s = k - r - d`, the coefficient of `Q'(X) := Q_Y(X)^{2^s}` at
position `j` is `(a_{j + r + d})^{2^s}`; for `j in {k-d, ..., k-1}` the
index `j + r + d` (mod `k`) lies in the window, so those coefficients
vanish and `Q'` has degree at most `2^{k-d-1}`. But as a function
`Q'(z) = Q_Y(z)^{2^s}`, so `Q'` vanishes on the kernel of `Q_Y`, which
has `2^{k-d}` elements (the map `Q_Y` has `d`-dimensional image). A
nonzero polynomial of degree at most `2^{k-d-1}` has at most `2^{k-d-1}`
roots, and `Q'` is not the zero function (its image is
`Y^{2^s} != {0}`), so by Fact 1 `Q'` is a nonzero polynomial —
contradiction. For (2): suppose `a_j = 0` on a cyclic window
`{w, ..., w + d - 2}` of length `d - 1` with `w != k - d + 1`. (If
`d = k` the only length-`(k-1)` zero run avoiding the nonzero position
`0 = k - d` is the top run, so the hypothesis with `w != k - d + 1`
cannot occur; assume `d <= k - 1`.) Put
`s = k - (w + d) + 1 mod k`; then `s in {1, ..., k - 1}` (as
`s = 0` iff `w = k - d + 1`). The coefficient of
`Q'(X) := Q_Y(X)^{2^s}` at position `j` is `(a_{j - s})^{2^s}`, and
`j - s = j + w + d - 1 - k`, so the index at position `j` is
`j + w + d - 1` (mod `k`). For `j in {k-d+1, ..., k-1}` this index
ranges over `{w, ..., w + d - 2}`, the window, so those coefficients
vanish and `Q'` has degree at most `2^{k-d}`; at position `j = k - d`
the index is `w - 1` (mod `k`), so that coefficient is
`c := (a_{w - 1})^{2^s}` (the coefficient just *before* the window).
Hence `Q' - c Q_Y` has degree at most `2^{k-d-1}` (the
position-`k-d` coefficients agree, as `a_{k-d} = 1`) and vanishes on
the `2^{k-d}` roots of `Q_Y` (where `Q'(z) = Q_Y(z)^{2^s} = 0`); a
nonzero univariate polynomial has at most as many roots as its degree,
so `Q' - c Q_Y` is the zero polynomial:
`Q_Y(z)^{2^s} = c Q_Y(z)` for all `z`. If `c = 0` then `Q' = 0` as a
polynomial, contradicting that `Q'` has image `Y^{2^s} != {0}`; so
`c != 0`, and
with `Y = Q_Y(K)` we get `y^{2^s} = c y` for all `y in Y`. By BK
Claim 3.8 (`a^{2^s} = c a` implies `a in b F_{2^{gcd(s,k)}}` for a
suitable `b`), `Y ⊆ b * F_{2^{gcd(s,k)}}`. Since `k` is prime and
`1 <= s <= k - 1`, `gcd(s, k) = 1`, so `Y ⊆ b * F_2`, a space of
dimension `1 < d` — contradiction. A zero run of length `>= d` would
contain a length-`(d-1)` window to which the same applies unless it is
the top run, whose length is exactly `d - 1` by part (1) (`a_{k-d} != 0`
and `a_0 != 0` bound it on both sides). ∎

For a linearized `Q(X) = sum_j a_j X^{2^j}` define its *associated matrix*
`M_Q` by `(M_Q)_{s,j} = (a_{j-s})^{2^s}` (indices mod `k`), and let `M_Q^0`
be its zero-one indicator, `(M_Q^0)_{s,j} = [a_{j-s} != 0]`. By the
Frobenius shift, `(M_Q)_{s,j}` is the coefficient of `X^{2^j}` in
`Q(X)^{2^s} mod (X^{2^k} - X)`.

**Lemma D2 (coefficient formula; BK Lemma 3.13).** Let `Q_1, ..., Q_t` be
linearized, `E_i, B_i ⊆ {0, ..., k-1}` with `|E_i| = |B_i| > 0`,
`alpha_i = sum_{j in E_i} 2^j`, `beta_i = sum_{j in B_i} 2^j`, `mu != 0`.
The coefficient of `M = prod_i X_i^{beta_i}` in
`R = Tr( mu * prod_i Q_i(X_i)^{alpha_i} ) mod (X_i^{2^k} - X_i)` is

```
c_M = sum_{rho=0}^{k-1} mu^{2^rho} * prod_i Perm( M_{Q_i}[E_i, B_i - rho] )^{2^rho},
```

where `Perm` is the permanent and `B_i - rho` is the cyclic shift.
Consequently, if for `rho = 0` every `M_{Q_i}^0[E_i, B_i]` is (after
ordering columns) lower triangular with nonzero diagonal, while for every
`rho in {1, ..., k-1}` some `M_{Q_i}^0[E_i, B_i - rho]` has an all-zero
column, then `c_M != 0` (the `rho = 0` summand is a product of nonzero
diagonal products; every other summand contains a zero permanent).

*Proof.* Expand the trace as a sum of Frobenius powers; the coefficient
factors over `i`. For one factor, the coefficient of `X_i^{beta_i}` in
`Q_i(X_i)^{alpha_i 2^rho} = prod_{j in E_i} Q_i(X_i)^{2^{j+rho}}` is, by the
matrix definition and the formula for the coefficient of a product of
linearized polynomials, the permanent of the minor with rows `E_i + rho`
and columns `B_i`; reindexing rows by the cyclic shift gives the stated
form. ∎

## D.2 The coefficient theorem (BK Theorem 6.1, linear prime-degree case)

**Theorem D3.** Let `Y_1, ..., Y_t <= K` be linear subspaces of dimensions
`d_i > 1`, let `Q_i` be the image-subspace polynomial of `Y_i`, let
`1 <= e_i < d_i - 1` integers and `alpha_i = sum_{j=0}^{e_i} 2^j`
(so `(Q_i)^{alpha_i}` is the product of the `e_i + 1` Frobenius powers
`Q_i^{2^j}`, `0 <= j <= e_i`). If

```
sum_i (d_i - (e_i + 1)) > k - max_i d_i + 1,
```

then for every `mu != 0` the polynomial
`R = Tr( mu prod_i Q_i(X_i)^{alpha_i} ) mod (X_i^{2^k} - X_i)` has a
monomial `prod_i X_i^{beta_i}` with `wt_2(beta_i) = e_i + 1` and nonzero
coefficient.

*Proof.* Reorder so `d_1 = max d_i`. Put `E_i = {0, ..., e_i}` (the row
sets of Lemma D2). Define offsets `l_1 = 0`,
`l_{i+1} = l_i + d_i - (e_i + 1)` (integers, not mod `k`), and let
`t_0` be minimal with `sum_{i <= t_0} (d_i - (e_i + 1)) >= k - d_1 + 1`
(it exists by the hypothesis of the theorem). For `i <= t_0` one has
`l_i <= k - d_1 <= k - d_i` (the partial sums before `t_0` are
`< k - d_1 + 1`), so `l_i` is a genuine index in `{0, ..., k - d_i}`,
not a residue class wrap. Writing `Q_i(X) = sum_j a_{i,j} X^{2^j}`, for
`i <= t_0` choose `B_i ⊆ {0, ..., k-1}`, `|B_i| = e_i + 1`, as follows.
If `a_{i, l_i} != 0`, let
`B_i = {l_i} ∪ {k - d_i + 1, ..., k - d_i + e_i}` (disjoint, since
`l_i <= k - d_i`). If `a_{i, l_i} = 0`, let `j_i` be the next index after
`l_i` with `a_{i, j_i} != 0` and `j_i'` the previous index before `l_i`
with `a_{i, j_i'} != 0`; indices here are genuine, not cyclic:
`j_i <= k - d_i` because `a_{i, k-d_i} != 0` (Lemma D1(1)) blocks the
zero run containing `l_i` from above, and `j_i' >= 0` because
`a_{i, 0} != 0`. The open interval `(j_i', j_i)` is a zero run containing
`l_i`, of length `q_i = j_i - j_i' - 1 >= 1`; with
`s_i = min(q_i, e_i)` let
`B_i = {j_i} ∪ {j_i' + 1, ..., j_i' + s_i} ∪ {k - d_i + s_i + 1, ..., k - d_i + e_i}`,
again with `e_i + 1` distinct elements (the middle block ends at
`j_i' + s_i <= j_i - 1 <= k - d_i - 1`, below the top block).
For `i > t_0` no covering will be needed from `B_i`, and we take the
canonical choice `B_i = {k - d_i, ..., k - d_i + e_i}`.
Set `beta_i = sum_{j in B_i} 2^j`, so `wt_2(beta_i) = e_i + 1` in all
cases. We verify the two hypotheses of Lemma D2 for
`M = prod X_i^{beta_i}`.

*Claim 1: each `M_{Q_i}^0[E_i, B_i]` is lower triangular with nonzero
diagonal, for a suitable column order.* The `(s, b)` entry is
`[a_{i, b-s} != 0]`, `s in E_i = {0, ..., e_i}`. For `i > t_0` (canonical
`B_i`), order the columns increasingly: the diagonal `(s, k - d_i + s)`
gives `a_{i, k-d_i} != 0` by Lemma D1(1), and an above-diagonal entry
(s < j) gives index `k - d_i + (j - s)` in `{k - d_i + 1, ..., k - 1}`
(using `j <= e_i < d_i`), a top-run zero by Lemma D1(1). Now let
`i <= t_0` and order the columns of `B_i` as listed. Diagonal entries:
`(0, l_i)` resp. `(0, j_i)` give `a_{i,l_i}` resp. `a_{i, j_i}`, nonzero
by construction; `(s, k - d_i + s)` gives `a_{i, k-d_i} != 0` by
Lemma D1(1); in the second construction, `(s, j_i' + s)` for
`1 <= s <= s_i` gives `a_{i, j_i'} != 0`. Above-diagonal entries: a
column `k - d_i + j` (`1 <= j <= e_i`) against a row `s < j` gives index
`k - d_i + (j - s)` in `{k - d_i + 1, ..., k - 1}`, a top-run zero by
Lemma D1(1); a column `j_i' + j` (`1 <= j <= s_i`) against a row `s < j`
gives index in the zero run `(j_i', j_i)`, since
`j_i' + (j - s) <= j_i' + s_i <= j_i' + q_i = j_i - 1`; a column
`k - d_i + j` with `j > s_i` against a row `s <= s_i` likewise lands in
the top run. Hence the permanent of each `rho = 0` minor is the product
of its diagonal, nonzero.

*Claim 2: for every shift `rho in {1, ..., k-1}` some
`M_{Q_i}^0[E_i, B_i - rho]` has an all-zero column.* The column of
`M_{Q_i}^0[E_i, B_i - rho]` indexed by `b in B_i` is
`([a_{i, b - rho - s} != 0])_{s=0..e_i}`; it is all zero if the cyclic
interval `[b - rho - e_i, b - rho]` lies inside the top zero run
`{k - d_i + 1, ..., k - 1}`, which happens whenever
`rho in [b + 1, b + d_i - (e_i + 1)]` (cyclic intervals in `Z_k`): then
`b - rho in [k - d_i + e_i + 1, k - 1]` and the whole window of length
`e_i + 1` stays in the run. So it suffices to cover
`{1, ..., k - 1}` by these intervals over all `i <= t_0` and
`b in B_i` (the partial sums give `l_{t_0 + 1} <= k - 1`, since the sum
before the last term is `< k - d_1 + 1` and each term is at most
`d_1 - 2`). For `i = 1`:
`l_1 = 0` and `a_{i,0} != 0` by Lemma D1(1), so `0 in B_1`, contributing
`[1, d_1 - (e_1 + 1)] = [1, l_2]`, and the high elements of `B_1`
contribute `[k - d_1 + 2, k - 1]`. For `1 < i <= t_0` we show
`∪_{b in B_i} [b + 1, b + d_i - (e_i + 1)] ⊇ [l_i + 1, l_{i+1}]` as integer
intervals (note `l_i <= k - d_1 <= k - d_i`). If `a_{i, l_i} != 0` then
`l_i in B_i` and its interval is exactly `[l_i + 1, l_{i+1}]`. If
`a_{i, l_i} = 0` and `q_i < e_i`, then
`{j_i} ∪ {j_i' + 1, ..., j_i' + q_i} = {j_i' + 1, ..., j_i}` contains
`l_i`, done. If `q_i >= e_i`, Lemma D1(2) bounds `q_i <= d_i - 2`: the
run `(j_i', j_i)` is not the top run (it contains `l_i <= k - d_i`, while
the top run lies in `{k - d_i + 1, ..., k - 1}` and `j_i' >= l_i - q_i`
with `a_{i,j_i'} != 0` excludes `j_i'` being inside the top run unless
`j_i' = k - d_i`, which would force `l_i > k - d_i`, contradiction), so
the run has length at most `d_i - 2` by D1(2). The intervals from
`b in {j_i' + 1, ..., j_i' + e_i} ⊆ B_i` then form the contiguous interval
`[j_i' + 2, j_i' + d_i - 1]`, which contains `j_i = j_i' + q_i + 1 <= j_i'
+ d_i - 1` and starts at `j_i' + 2 <= l_i + 1` (as `j_i' < l_i`); together
with `j_i`'s own interval `[j_i + 1, j_i + d_i - (e_i + 1)]`, whose right
end is at least `l_i + d_i - (e_i + 1) = l_{i+1}`, the union covers
`[l_i + 1, l_{i+1}]` (contiguity: the first interval ends at
`j_i' + d_i - 1 >= j_i - 1`, and the second starts at `j_i + 1`; the
possible single point `j_i` is covered because `j_i <= j_i' + d_i - 1`
lies in the first interval). Chaining `i = 1, ..., t_0` covers
`[1, l_{t_0 + 1}]` with `l_{t_0 + 1} >= k - d_1 + 1`, and the `i = 1` high
intervals cover `[k - d_1 + 2, k - 1]`; together, all of
`{1, ..., k - 1}`. ∎ (Theorem D3)

## D.3 Mixing and the main disperser theorem

**Lemma D4 (Bourgain decomposition).** Let `A ⊆ K^r` be `F_2`-affine,
`A = gamma + V`. Then there are linear subspaces `Y_i <= K` and `F_2`-linear
maps `sigma_{ij} : Y_j -> K` (`j < i`) with
`A = { (gamma_i + y_i + sum_{j<i} sigma_{ij}(y_j))_i : y_i in Y_i }` and
`dim A = sum_i dim Y_i`.
*Proof.* Induction on `r` via a linear section of the projection of `V`
to the first block: `V` is the graph over `Y_1 = pi_1(V)` of a linear map
into the remaining blocks plus `V ∩ (0 x K^{r-1})`; recurse. ∎

**Lemma D5 (triangular mixing does not change the key coefficient;
BK Lemma 7.2).** Let `P_i(X_i)` be linearized, `P_{ij}(X_j)` linearized
(`j < i`), `gamma_i in K`, `R_i = P_i(X_i) + sum_{j<i} P_{ij}(X_j) +
gamma_i`, `I_0 ⊆ [r]`, `|I_0| = t`, `e_i >= 1`, `alpha_i = sum_{j=0}^{e_i}
2^j`. For any `(beta_i)_{i in I_0}` with `wt_2(beta_i) = e_i + 1`, the
coefficient of `prod_{i in I_0} X_i^{beta_i}` in
`Tr(mu prod_{i in I_0} R_i^{alpha_i}) mod (X_i^{2^k} - X_i)` equals its
coefficient in `Tr(mu prod_{i in I_0} P_i(X_i)^{alpha_i}) mod (...)`.
*Proof.* Expand the trace into Frobenius twists `rho`. Fix `rho` and write
`P̂_i = sum_{j<i} P_{ij}(X_j) + gamma_i` (no `X_i`). Expand
`(P_i + P̂_i)^{alpha_i 2^rho} = prod_{l=0}^{e_i} (P_i + P̂_i)^{2^{l+rho}}`;
choosing from each factor either a `P_i`-power or a `P̂_i`-power writes it
as `(P_i)^{alpha_i 2^rho} + G_i`, where every monomial of `G_i` has
`F_2`-degree at most `e_i` in `X_i` (each `P̂_i^{2^{l+rho}}` factor
contributes `F_2`-degree `0` in `X_i`, the others at most `1`). Now argue
by descending induction on `I_0 = {i_1 < ... < i_t}`. For the largest
index `i_t`: factors indexed by `i < i_t` involve only the variables
`X_1, ..., X_i`, hence not `X_{i_t}`; so a monomial of `F_2`-degree
`e_{i_t} + 1` in `X_{i_t}` in the full product forces the
`(P_{i_t})^{alpha_{i_t} 2^rho}` choice at block `i_t` (the `G_{i_t}`
choice contributes at most `e_{i_t}`), and the coefficient of the target
monomial factors as the coefficient of `X_{i_t}^{beta_{i_t}}` in
`(P_{i_t})^{alpha_{i_t} 2^rho}` times the coefficient of the remaining
monomial in the product over `i < i_t`. Repeating the argument for
`i_{t-1}, ..., i_1` shows only the all-`P` choice contributes. Summing
over `rho` gives the claim. ∎

**Theorem D6 (BK Theorem 7.3, restated).** Let `1 <= t < r`, let `k` be
prime with `k t >= r(r + 1)`, and for `i in [r]` let `e_i = r + 1 - i`,
`alpha_i = sum_{j=0}^{e_i} 2^j` (so the exponent of the power map is
`alpha_i = 2^{e_i + 1} - 1 = 2^{r + 2 - i} - 1`). For `mu != 0` let
`F(x) = Tr( mu * Sym_r^t( x_1^{alpha_1}, ..., x_r^{alpha_r} ) )` on
`K^r`. Then `F` is nonconstant on every `F_2`-affine subspace
`A ⊆ K^r` with `dim A > kr/t + kt + r(r + 1)`.

*Proof.* Let `A` be such a subspace, decomposed as in Lemma D4 with
`d_i = dim Y_i`, `sum d_i = dim A`. Let `Q_i` be the image-subspace
polynomial of `Y_i`, and `Q_{ij}` the linearized polynomial with
`Q_{ij}(X) = sigma_{ij}(Q_j(X))` (it exists: `sigma_{ij} o Q_j` is an
`F_2`-linear map `K -> K`, and evaluation is a bijection). Put
`R_i(X_1, ..., X_r) = Q_i(X_i) + sum_{j<i} Q_{ij}(X_j) + gamma_i`; the map
`H : X |-> (R_1, ..., R_r)` has image exactly `A` (as `X_j` varies,
`Q_j(X_j)` ranges over all of `Y_j` independently). By Fact 1 it suffices
to show the reduced polynomial `h = F o H mod (X_i^{2^k} - X_i)` has a
positive-degree monomial with nonzero coefficient. Expanding the
symmetric sum,
`h = sum_{I ⊆ [r], |I| = t} Tr( mu prod_{i in I} R_i^{alpha_i} ) mod (...)`.

Let `S = {i : d_i > r + 1}`. Then `|S| >= t`: otherwise, since each
`d_i <= k`,
`dim A = sum d_i <= (t - 1) k + r (r + 1) < kr/t + kt + r(r + 1)`,
because `kr/t + kt - (t - 1) k = kr/t + k > 0` — contradicting the
hypothesis on `dim A`. Also
`sum_{i in S} (d_i - (r + 1)) >= sum_{i=1}^r (d_i - (r + 1)) = dim A -
r(r + 1) > kr/t + kt`, since the terms with `i not in S` are `<= 0`.
Choose `I_0 ⊆ S` with `|I_0| = t` maximising `sum_{i in I_0} d_i`; the
maximum is at least the average over `t`-subsets of `S`, so with
`|S| <= r`,
`sum_{i in I_0} (d_i - (r + 1)) >= (t / |S|) (kr/t + kt) >= (t/r)(kr/t +
kt) = k + k t^2 / r >= k + t (r + 1)`, where the last step uses the
hypothesis `kt >= r(r + 1)`, i.e. `kt/r >= r + 1`. For `i in I_0`,
`e_i + 1 = r + 2 - i <= r + 1`, hence
`sum_{i in I_0} (d_i - (e_i + 1)) >= sum_{i in I_0} (d_i - (r + 1)) >=
k + t(r + 1)`, which is `> k - max_{i in I_0} d_i + 1` because
`t(r + 1) + max_{i in I_0} d_i - 1 >= (r + 1) + (r + 2) - 1 = 2r + 2
>= 6 > 0` (using `t >= 1` and `max_{i in I_0} d_i >= r + 2`), and
`e_i <= r < d_i - 1` (as `d_i > r + 1`). Theorem D3 applied to
`(Y_i)_{i in I_0}` therefore gives a monomial
`M = prod_{i in I_0} X_i^{beta_i}`, `wt_2(beta_i) = e_i + 1`, with
nonzero coefficient in `Tr(mu prod_{i in I_0} Q_i(X_i)^{alpha_i})`,
hence by Lemma D5 in the `I_0`-summand of `h`.

No other summand contains `M` with nonzero coefficient. Indeed, fix
`J != I_0`, `|J| = t`. In the `J`-summand, expand the trace; for some
twist `rho`, `M` would have to appear in
`prod_{j in J} prod_{l=0}^{e_j} R_j^{(rho, l)}`, where each
`R_j^{(rho, l)} = R_j^{2^{l + rho}} mod (X^{2^k} - X)` is a sum of
monomials of `F_2`-degree at most `1`, each a constant or a power of a
single variable `X_s` with `s <= j`. If `M` appears in this product with
nonzero coefficient, some choice of one monomial per factor `(j, l)`
multiplies to `M`; set `lambda_j(l) = s` if the chosen monomial involves
`X_s`, else `lambda_j(l) = 0` (constant chosen). The `F_2`-degree in
`X_i` of any such product monomial is at most the number of factors
choosing `X_i` (each contributes `F_2`-degree at most `1` in `X_i`, and
reduction mod `X_i^{2^k} - X_i` never increases `F_2`-degree), while
`M` has `F_2`-degree `e_i + 1 = r + 2 - i` in `X_i` for `i in I_0`; so
necessarily `sum_j |lambda_j^{-1}(i)| >= r + 2 - i` for all `i in I_0`
(BK Claim 7.4's hypothesis (27)). Also
`sum_{i in I_0} |lambda_j^{-1}(i)| <= e_j + 1 = r + 2 - j` for each `j`.
Summing gives `sum_{i in I_0} i >= sum_{j in J} j`. More generally, for
every `s >= 0`, counting with weights `lambda_j(l)^s` and using
`lambda_j(l) <= j` yields (BK Claim 7.4)
`sum_{i in I_0} i^{s+1} - sum_{j in J} j^{s+1} >= (r + 2) (sum_{i in I_0}
i^s - sum_{j in J} j^s) >= 0` by induction on `s`. Hence all power sums of
`I_0` dominate those of `J`. The largest element of `I_0` is then at
least that of `J`; it cannot be smaller (a variable `X_i` with
`i > max J` is chosen by no factor, contradicting the count for
`i = max I_0 in I_0`), so they are equal, say `i_1 = j_1`, and every
factor of block `j_1` must choose `X_{i_1}` (its `r + 2 - j_1` factors
exactly meet the required count, and no other block can choose
`X_{i_1}`). Removing `i_1, j_1` and repeating the argument for the
second-largest elements, and so on, gives `I_0 = J` — contradiction.
Hence `M` appears in `h` with the nonzero coefficient from the `I_0`
summand alone, `h` is a nonconstant polynomial, and by Fact 1 `F` is
nonconstant on `A`. ∎

**Corollary D7 (the form used in `proof.md`).** Let `k` be prime,
`r >= 2`, `1 <= t < r`, `kt >= r(r + 1)`. Let `phi : F_2^m -> K^r` be an
injective `F_2`-linear map (block embedding with zero padding), let
`pi : K -> F_2` be the constant-coefficient functional, and
`F(z) = pi( Sym_r^t( y_1, ..., y_r ) )` with
`y_i = phi_i(z)^{2^{i+1} - 1}` (ascending exponents, as in BK Theorem 2.2
equation (4)). Then `F` is nonconstant on every affine subspace of
`F_2^m` of dimension `> kr/t + kt + r(r + 1)`, hence of dimension at
least `d := floor(kr/t) + kt + r(r + 1) + 1`.

*Proof.* By Fact 2, `pi(z) = Tr(mu z)` with `mu != 0`. Reversing the block
order is an `F_2`-linear permutation of `K^r` mapping affine subspaces to
affine subspaces of the same dimension, and it turns the ascending
exponents `2^{i+1} - 1` into the descending exponents `2^{r + 2 - i} - 1`
of Theorem D6, without changing `Sym_r^t`. So `F`, viewed on `K^r`, is
nonconstant on every affine subspace of dimension
`> kr/t + kt + r(r + 1)` by Theorem D6. An affine subspace
`A ⊆ F_2^m` maps under the injective linear `phi` to an affine subspace
of `K^r` of the same dimension on which the `K^r`-function is
nonconstant; hence `F` is nonconstant on `A`. The integer form follows
since `floor(kr/t) + kt + r(r + 1) + 1 > kr/t + kt + r(r + 1)`. ∎

Provenance: construction and overall proof strategy from Ben-Sasson–
Kopparty, *Affine dispersers from subspace polynomials*, SIAM J. Comput.
41(4):880–914, 2012 (preprint 2010-03-12, SHA-256
922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf):
Theorem 2.2 and equation (4); Lemma 3.7 (proved in full above as the
Claim and Lemma D0, without reliance on the citation); Lemma 3.10;
Lemma 3.13; Theorem 6.1 with
Claims 6.2–6.3 and Proposition 6.4; Lemma 7.1; Lemma 7.2; Theorem 7.3 with
Claim 7.4. The restriction to prime degree `k` (BK's construction also
takes the field degree prime), the vacuity of the subfield hypotheses by
Lemma D1(2), and the exact threshold `kr/t + kt + r(r + 1)` (in place of
the rounded `6 m^{4/5}` corollary) are as derived above.
