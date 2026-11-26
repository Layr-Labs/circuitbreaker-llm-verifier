# Self-contained algebraic affine disperser proof

This file proves the finite-field affine disperser lemma used in `proof.md`,
specialized to characteristic two and prime field degree `k`. All subspace
dimensions are over `F_2`. The construction follows Ben-Sasson and Kopparty
(BK, 2010, Theorem 7.3), with every algebraic and combinatorial step proved
directly here so no external disperser theorem is invoked as a black box.

**Lemma (Payload Affine Disperser).** Let `K` be the finite field of size
`q = 2^k` with `k` prime. Suppose integers `r, t` satisfy `1 <= t < r` and
`k*t >= r*(r+1)`. For any nonzero `F_2`-linear functional `pi: K -> F_2`, set
`w_i = r + 2 - i` for `1 <= i <= r` and define

    f(a_1, ..., a_r) = pi( sum_(|I| = t) product_(i in I) a_i^(2^(w_i) - 1) ).

Then `f` is nonconstant on every affine subspace `A` of `K^r` of
`F_2`-dimension strictly greater than `k*r/t + k*t + r*(r+1)`.

## 1. Finite-field and polynomial preliminaries

Every function `K^r -> K` has a unique polynomial representative in
`K[X_1, ..., X_r]` of degree `< q` in each variable, obtained by reducing
modulo the ideal `(X_1^q - X_1, ..., X_r^q - X_r)`. Uniqueness follows by
induction on `r` from the univariate root bound (a nonzero polynomial of
degree `< q` has fewer than `q` roots in `K`). In particular, if the reduced
representative of a function has a nonzero coefficient at a positive-degree
monomial, the function is nonconstant on `K^r`.

A linearized polynomial over `K` has the form `sum_(j=0)^(k-1) a_j X^(2^j)`
with `a_j in K`. Because Frobenius `z -> z^2` is `F_2`-linear, every
linearized polynomial defines an `F_2`-linear map `K -> K`. The space of
linearized polynomials and the space `End_(F_2)(K)` both have `F_2`-dimension
`k^2`, and the evaluation map is injective by the univariate root bound
(`2^(k-1) < q`). Hence linearized polynomials correspond bijectively to
`F_2`-linear maps `K -> K`.

The trace `Tr(z) = sum_(h=0)^(k-1) z^(2^h)` maps `K` into `F_2` because
`Tr(z)^2 = Tr(z)`. It is not identically zero on `K` because its reduced
polynomial has degree `2^(k-1) < q`. For any nonzero `mu in K`, `z -> Tr(mu*z)`
is therefore a nonzero `F_2`-linear functional on `K`. The map
`mu -> (z -> Tr(mu*z))` is an injective `F_2`-linear map between two
`k`-dimensional `F_2`-vector spaces, hence bijective. Thus our nonzero
functional satisfies `pi(z) = Tr(mu*z)` for a unique nonzero `mu in K`.

For an integer `b in {0, ..., q-1}`, let `wt(b)` be the Hamming weight (number
of 1s) of the binary expansion of `b`. Multiplying `h` factors
`X^(2^(v_1)), ..., X^(2^(v_h))` and reducing modulo `X^q - X` adds the
residues `2^(v_m mod k)` modulo `2^k - 1` (with `X^(2^k) = X`). Each binary
carry replaces `2^j + 2^j` by `2^((j+1) mod k)` and strictly decreases the
count of power-of-two summands by 1. Since a wraparound carry never creates
new summands, the reduced exponent `b` always satisfies `wt(b) <= h`.
Furthermore, when `h < k`, equality `wt(b) = h` occurs if and only if no
carry occurs at all—that is, the residues `v_1, ..., v_h (mod k)` are
pairwise distinct and `b = sum_(m=1)^h 2^(v_m mod k)`. Applying a Frobenius
power `z -> z^(2^rho)` cyclically shifts binary digits modulo `k` and
preserves `wt`. In particular, binary weight is subadditive under polynomial
multiplication and reduction modulo `X^q - X`.

## 2. Subspace image polynomials and zero-run uniqueness

For an `F_2`-linear subspace `V <= K` of dimension `d`, define
`P_V(X) = product_(v in V) (X - v)`. We show by induction on `d` that `P_V` is
a monic linearized polynomial of degree `2^d`: for `V = {0}`, `P_{0}(X) = X`;
for `V = U + F_2*b` with `b not in U`,

    P_V(X) = P_U(X) * P_U(X + b) = P_U(X)^2 + P_U(b) * P_U(X),

which is monic and linearized of degree `2^(dim(U)+1)`.

Assume now `2 <= d <= k`. Since `P_V` is `F_2`-linear with kernel `V`, its
image `V' = P_V(K)` has dimension `k - d`. Set `Q = P_(V')`, a monic
linearized polynomial of degree `2^(k-d)` whose `2^(k-d)` distinct roots are
all in `K`. The composition `Q(P_V(X))` is monic of degree `2^k = q` and
vanishes on all of `K`, so `Q(P_V(X)) = X^q - X`. We also have
`P_V(Q(X)) = X^q - X`: let `D(X) = P_V(Q(X)) - (X^q - X)`. Because `P_V` is
additive and has coefficients in `K = F_q` (so `c^q = c`),

    D(P_V(X)) = P_V(Q(P_V(X))) - (P_V(X)^q - P_V(X))
              = P_V(X^q - X) - (P_V(X^q) - P_V(X)) = 0.

Since `P_V(X)` is nonconstant, `D(P_V(X)) = 0` forces `D(X) = 0` as a
polynomial (if `D` were a nonzero constant `c` then `D(P_V(X)) = c != 0`; if
`deg(D) >= 1` then `deg(D(P_V)) = deg(D)*deg(P_V) >= 1`). Hence
`P_V(Q(X)) = X^q - X` vanishes on `K`, so `Q(K) <= ker(P_V) = V`. Because
`dim(ker(Q)) = k - d`, we have `dim(Q(K)) = d = dim(V)`, hence `Q(K) = V`.

Write `Q(X) = sum_(j=0)^(k-1) a_j X^(2^j)` with indices extended periodically
modulo `k`. Since `Q` is monic of degree `2^(k-d)`, we have `a_(k-d) = 1` and
`a_j = 0` for `j in {k-d+1, ..., k-1}`, a cyclic zero run of length `d - 1`.

*No cyclic zero run of length `d`.* Suppose `d` consecutive coefficients
modulo `k` were zero. Choosing `rho in {0, ..., k-1}` to shift that run to
`{k-d, ..., k-1}`, the reduced polynomial `Q^(2^rho) mod (X^q - X)` would be
nonzero of degree at most `2^(k-d-1)` while vanishing on all `2^(k-d)` roots
of `Q` in `K`, contradicting the root bound. In particular, `a_0 != 0` and
`a_(k-d) = 1 != 0`.

*Unique cyclic zero run of length `d - 1` (prime `k`, `d >= 2`).* Suppose
there were another cyclic zero run of length `d - 1`. Shifting it to
`{k-d+1, ..., k-1}` by Frobenius power `2^h` with `0 < h < k` produces a
reduced polynomial `Q^(h)(X) = Q(X)^(2^h) mod (X^q - X)` of degree at most
`2^(k-d)` that vanishes on every root of `Q` in `K`. Since `Q` has `2^(k-d)`
distinct roots in `K` and `Q^(h)` is nonzero, `Q^(h)(X) = c * Q(X)` for some
nonzero `c in K` (subtracting the appropriate multiple of `Q` yields a
polynomial of degree `< 2^(k-d)` with `2^(k-d)` roots). Evaluating at any
`x in K` and setting `z = Q(x) in V` gives `z^(2^h) = Q^(h)(x) = c*z`. Fix a
nonzero `z_0 in V`; then every `z in V` satisfies `(z/z_0)^(2^h) = z/z_0` as
well as `(z/z_0)^(2^k) = z/z_0`. Because `k` is prime and `0 < h < k`,
`gcd(h, k) = 1`, so Bezout gives integers `u, v` with `u*h + v*k = 1`.
Iterating the invertible Frobenius automorphism shows `(z/z_0)^2 = z/z_0`,
so `z/z_0 in F_2` for all `z in V`. Thus `dim(V) <= 1`, contradicting `d >= 2`.
Consequently, every cyclic zero run of `Q` other than `{k-d+1, ..., k-1}` has
length at most `d - 2`.

## 3. Nonzero coefficient for independent subspaces

Let `V_1, ..., V_t <= K` be linear subspaces with dimensions `d_i >= 2` and
image polynomials `Q_i(X) = sum_(j=0)^(k-1) a_(i,j) X^(2^j)`. Let integers
`e_i` satisfy `1 <= e_i < d_i - 1`, and put `h_i = e_i + 1` and
`delta_i = d_i - h_i = d_i - e_i - 1 >= 1`. Note `h_i < d_i <= k`. Assume

    sum_(i=1)^t delta_i >= k - max_i d_i + 1.

We prove that the reduced polynomial of

    Tr( mu * product_(i=1)^t Q_i(X_i)^(2^(h_i) - 1) )

has a nonzero coefficient at some monomial `product_(i=1)^t X_i^(b_i)` with
`wt(b_i) = h_i` for every `i`.

For each `i`, define the `k`-by-`k` matrix `T_i` (indices modulo `k`) by

    T_i[u, v] = a_(i, v - u)^(2^u),

which is the coefficient of `X^(2^v)` in `Q_i(X)^(2^u) mod (X^q - X)`. Let
`A_i = (0, 1, ..., e_i)`. Given an ordered tuple `B_i` of `h_i` distinct
residues modulo `k`, set `b_i = sum_(v in B_i) 2^v`, so `wt(b_i) = h_i < k`.
Expanding `Q_i(X_i)^(2^(h_i)-1) = product_(u=0)^(e_i) Q_i(X_i)^(2^u)` and
applying `Tr(z) = sum_(rho=0)^(k-1) z^(2^rho)`: in the `rho`-th trace
summand, raising `X_i^(2^v)` to `2^rho` shifts exponent bit `v` to `v + rho`.
By Section 1, since `h_i < k`, a product of `h_i` terms produces `X_i^(b_i)`
only when there are no carries, meaning the chosen columns form a bijection
from `A_i` to `B_i - rho (mod k)`. Therefore the coefficient of
`product_i X_i^(b_i)` in the reduced polynomial is

    sum_(rho=0)^(k-1) mu^(2^rho)
        * product_(i=1)^t perm(T_i[A_i, B_i - rho])^(2^rho).          (1)

We choose `B_1, ..., B_t` so that the `rho = 0` summand in (1) is nonzero and
every `rho in {1, ..., k-1}` makes at least one minor have an all-zero column.

Reorder the `t` indices so `d_1 = max_i d_i`. Let `J <= t` be the smallest
index with `sum_(i=1)^J delta_i >= k - d_1 + 1` (which exists by hypothesis).
Set `ell_i = sum_(j < i) delta_j` for `1 <= i <= J`, and `ell_i = 0` for
`i > J`. By minimality of `J`, for every `1 <= i <= J` we have
`ell_i <= k - d_1 <= k - d_i` (with `ell_1 = 0`). Moreover,
`sum_(i=1)^J delta_i = ell_J + delta_J <= (k - d_1) + (d_J - e_J - 1) <= k - 2`
because `d_J <= d_1` and `e_J >= 1`.

For each `i`, abbreviate `d = d_i`, `e = e_i`, `ell = ell_i`, `a_s = a_(i,s)`.
Notice `k - d + e <= k - 2 < k - 1` because `e <= d - 2`.

- **Case 1 (`a_ell != 0`).** Set the ordered tuple

      B_i = (ell, k - d + 1, k - d + 2, ..., k - d + e).

  Since `0 <= ell <= k - d`, these `e + 1` entries are strictly increasing in
  `{0, ..., k - 2}`, hence distinct modulo `k`. For `0 <= u, m <= e`, the
  entry `T_i[u, B_i[m]]` is `a_(B_i[m] - u)^(2^u)`. On the diagonal (`u = m`),
  the indices are `ell` (for `m = 0`) and `k - d` (for `1 <= m <= e`), all
  nonzero. Above the diagonal (`0 <= u < m <= e`), we have
  `B_i[m] - u = k - d + (m - u) in [k - d + 1, k - 2]`, lying in the final
  zero run. Thus `T_i[A_i, B_i]` is lower triangular with nonzero diagonal,
  so `perm(T_i[A_i, B_i]) != 0`.

- **Case 2 (`a_ell = 0`).** Since `a_0 != 0` and `a_(k-d) = 1 != 0`, we have
  `0 < ell < k - d`. Choose the largest `l < ell` and smallest `j > ell` with
  `a_l != 0` and `a_j != 0`, so `0 <= l < ell < j <= k - d` and `a_s = 0` for
  all `l < s < j`. Let `z = j - l - 1` be the zero-gap length; by Section 2,
  `1 <= z <= d - 2`. Put `v = min(z, e)` and set

      B_i = (j, l + 1, ..., l + v, k - d + v + 1, ..., k - d + e),

  where the last block is empty if `v = e`. These `e + 1` entries are distinct
  in `{1, ..., k - 2}`: `l + 1 < ... < l + v <= l + z = j - 1 < j <= k - d`,
  and when `v < e` every entry `k - d + m` (`m >= v + 1`) satisfies
  `k - d < k - d + m <= k - 2`. On the diagonal (`u = m`), `B_i[m] - m` equals
  `j` (`m = 0`), `l` (`1 <= m <= v`), and `k - d` (`v + 1 <= m <= e`), all
  nonzero. Above the diagonal (`0 <= u < m <= e`): if `1 <= m <= v`, then
  `B_i[m] - u = l + (m - u) in [l + 1, l + v] <= [l + 1, j - 1]`, giving 0;
  if `v + 1 <= m <= e`, then `B_i[m] - u = k - d + (m - u) in [k - d + 1, k - 2]`,
  again giving 0. Thus `T_i[A_i, B_i]` is lower triangular with nonzero
  diagonal, so `perm(T_i[A_i, B_i]) != 0`.

Next, for any column `b in B_i`, column `b - rho` of `T_i[A_i, B_i - rho]` is
all-zero whenever

    rho in [b + 1, b + delta_i]   (mod k).                            (2)

Indeed, writing `rho = b + s (mod k)` with `1 <= s <= delta_i`, for every
`0 <= u <= e_i` we have `b - rho - u = -s - u = k - (s + u) (mod k)`, and
`1 <= s + u <= delta_i + e_i = d_i - 1`, so `k - (s + u) in [k - d_i + 1, k - 1]`,
which is the final zero run of `Q_i`.

We now verify that the intervals (2) for `1 <= i <= J` cover `{1, ..., k - 1}`:

- For `i = 1`, `ell_1 = 0` and `a_(1,0) != 0`, so Case 1 applies:
  `B_1 = (0, k - d_1 + 1, ..., k - d_1 + e_1)`. Column `0` covers
  `[1, delta_1]`. The remaining columns `k - d_1 + 1, ..., k - d_1 + e_1`
  cover consecutive overlapping intervals whose union is
  `[k - d_1 + 2, k - d_1 + e_1 + delta_1] = [k - d_1 + 2, k - 1]`.
- For each `2 <= i <= J`, we show `B_i` covers `[ell_i + 1, ell_i + delta_i]`:
  - If `a_(i, ell_i) != 0`, then `ell_i in B_i`, directly covering
    `[ell_i + 1, ell_i + delta_i]`.
  - If `a_(i, ell_i) = 0` and `z < e_i`, then `v = z`, so `B_i` contains
    `l + 1, ..., l + z = j - 1`. Since `l < ell_i < j`, we have `ell_i in B_i`,
    again covering `[ell_i + 1, ell_i + delta_i]`.
  - If `a_(i, ell_i) = 0` and `z >= e_i`, then `v = e_i`, so `B_i` contains
    `l + 1, ..., l + e_i` and `j`. The first `e_i` columns cover
    `[l + 2, l + e_i + delta_i] = [l + 2, l + d_i - 1]`, and column `j` covers
    `[j + 1, j + delta_i]`. Because `z <= d_i - 2`, we have
    `j = l + z + 1 <= l + d_i - 1`, so there is no gap between the two
    intervals and their union is `[l + 2, j + delta_i]`. Since `l < ell_i < j`,
    this contains `[ell_i + 1, ell_i + delta_i]`.

Since `ell_i = sum_(j < i) delta_j`, the union of `[ell_i + 1, ell_i + delta_i]`
for `1 <= i <= J` is `[1, sum_(i=1)^J delta_i]`. Because
`sum_(i=1)^J delta_i >= k - d_1 + 1`, combining this with the tail interval
`[k - d_1 + 2, k - 1]` from `i = 1` covers every `rho in {1, ..., k - 1}`.
(And `j + delta_i <= (k - d_i) + (d_i - e_i - 1) = k - e_i - 1 <= k - 2`, so
no interval wraps around `k`.)

Consequently, in (1) every `rho != 0` term vanishes, leaving only the nonzero
`rho = 0` term `mu * product_(i=1)^t perm(T_i[A_i, B_i]) != 0`. This proves
the independent-subspace assertion.

## 4. Block-triangular parameterization of an arbitrary affine subspace

Let `A = gamma + W <= K^r` be an affine subspace of `F_2`-dimension `D`, with
`gamma in A` and `W` an `F_2`-linear subspace of `K^r`. Filter `W` by
`W_0 = W >= W_1 >= ... >= W_r = {0}`, where `W_i` consists of vectors in `W`
whose first `i` coordinates in `K` are zero. For each `i = 1, ..., r`, choose
a complementary subspace `U_i` with `W_(i-1) = U_i + W_i` (direct sum). By
definition of `W_i`, projection onto coordinate `i` maps `U_i` injectively
onto a subspace `V_i <= K` of dimension `d_i = dim(U_i)`, and expresses each
vector of `U_i` uniquely as `(0, ..., 0, z_i, sigma_(i+1,i)(z_i), ..., sigma_(r,i)(z_i))`
for `z_i in V_i` and `F_2`-linear maps `sigma_(j,i): V_i -> K` (`j > i`).
Since `W = U_1 + ... + U_r` (direct sum), we have `D = sum_(i=1)^r d_i` and

    A = { (gamma_i + z_i + sum_(j < i) sigma_(i,j)(z_j))_(i=1)^r : z_i in V_i }.   (3)

For each `i`, let `Q_i(X_i)` be the linearized image polynomial mapping `K`
onto `V_i` from Section 2 (taking `Q_i = 0` if `d_i = 0`, and any nonzero
linearized map onto `V_i` if `d_i = 1`). For each `j < i`, the composition
`sigma_(i,j) o Q_j: K -> K` is `F_2`-linear, so by Section 1 it is represented
by a linearized polynomial `Q_(i,j)(X_j)`. Define `R: K^r -> K^r` by

    R_i(X) = gamma_i + Q_i(X_i) + sum_(j < i) Q_(i,j)(X_j).

Then `R(K^r) = A`. To show `f` is nonconstant on `A`, it suffices to exhibit a
nonzero coefficient at a positive-degree monomial in the reduced polynomial of
`f(R(X))`.

## 5. Selecting `t` high-dimension blocks

Assume `D > k*r/t + k*t + r*(r+1)`. Let `S = {i in {1, ..., r} : d_i > r + 1}`.
Since at most `r` indices lie outside `S` and each has `d_i <= r + 1`,

    sum_(i in S) d_i >= D - r*(r+1) > k*r/t + k*t.

Since every `d_i <= k`, we have `|S| >= (sum_(i in S) d_i)/k > r/t + t > t`,
so `|S| >= t`. Let `I_0` be a `t`-element subset of `S` with largest values of
`d_i`. Because `|S| <= r`, averaging over `S` gives

    sum_(i in I_0) d_i >= (t / |S|) * sum_(i in S) d_i
                       >= (t / r) * sum_(i in S) d_i
                        > k + k*t^2/r.

Subtracting `t*(r+1)` and using the hypothesis `k*t >= r*(r+1)` yields

    sum_(i in I_0) (d_i - (r + 1)) > k + (t/r) * (k*t - r*(r+1)) >= k.   (4)

For every `i in I_0`, we have `1 <= i <= r`, so `e_i = r + 1 - i >= 1` and
`w_i = e_i + 1 = r + 2 - i <= r + 1 < d_i`. Thus `d_i >= r + 2 >= 3 > 1` and
`1 <= e_i < d_i - 1`. Furthermore, `delta_i = d_i - w_i >= d_i - (r + 1)`, so

    sum_(i in I_0) delta_i > k >= k - max_(i in I_0) d_i + 1.

Section 3 therefore produces a monomial `M = product_(i in I_0) X_i^(b_i)`
with `wt(b_i) = w_i` (note `w_i >= 2`, so `M` has positive degree) whose
coefficient in the reduced polynomial of
`Tr(mu * product_(i in I_0) Q_i(X_i)^(2^(w_i) - 1))` is nonzero.

## 6. Triangular mixing preserves the coefficient of `M`

Consider the `I_0`-summand of `f(R(X))`:

    Tr( mu * product_(i in I_0) R_i(X)^(2^(w_i) - 1) )
      = sum_(rho=0)^(k-1) mu^(2^rho) * product_(i in I_0)
          product_(u=0)^(w_i - 1) R_i(X)^(2^(u + rho)).

Fix `rho` and expand each factor `R_i(X)^(2^(u+rho))` as a sum of a constant
plus univariate linearized terms in `X_j` for `j <= i`. We claim that a choice
of one term from each factor can reduce modulo `(X_j^q - X_j)` to `M` only if
every factor of `R_i` selects a term from `Q_i(X_i)` (never a constant or an
`X_j` with `j < i`).

Prove this by descending induction over `i in I_0`. For the largest `i in I_0`,
no factor with smaller index `i' < i` contains `X_i`. Since `wt(b_i) = w_i`
and binary weight is subadditive (Section 1), all `w_i` factors from `R_i`
must select an `X_i`-term (from `Q_i(X_i)`). Having fixed those factors, none
of them contributes to any smaller variable, so the same argument applies to
the next-largest index in `I_0`, and so on down to the smallest. Therefore
each `rho`-summand has the exact same coefficient of `M` when `R_i(X)` is
replaced by `Q_i(X_i)`, and hence the `I_0`-summand of `f(R(X))` has the
exact same nonzero coefficient at `M`.

## 7. No other symmetric summand contributes to `M`

Let `J <= {1, ..., r}` be any `t`-element subset, and suppose some term
choice in the expansion of
`Tr(mu * product_(j in J) R_j(X)^(2^(w_j) - 1))` reduces to `M`. For each
`j in J` and `i in I_0`, let `n_(j,i) >= 0` be the number of the `w_j`
linearized factors of `R_j` that select a term in variable `X_i`. Because
`R_j` mentions only variables `X_i` with `i <= j`, and by binary-weight
subadditivity,

    sum_(i in I_0) n_(j,i) <= w_j   for all j in J,
    sum_(j in J)   n_(j,i) >= w_i   for all i in I_0,                 (5)
    n_(j,i) > 0  implies  i <= j.

Since `w_i = r + 2 - i` is strictly decreasing in `i`, every pair with
`n_(j,i) > 0` satisfies `i <= j` and therefore `w_i >= w_j`, with strict
inequality whenever `i < j`. Summing `n_(j,i)/w_i` and `n_(j,i)/w_j` gives

    t = sum_(i in I_0) w_i / w_i
     <= sum_(i in I_0, j in J) n_(j,i) / w_i
     <= sum_(i in I_0, j in J) n_(j,i) / w_j
     <= sum_(j in J) w_j / w_j = t.                                   (6)

Every inequality in (6) must be an equality. Hence no positive `n_(j,i)` can
have `i < j`, so `n_(j,i) > 0` forces `i = j`. For every `i in I_0`,
`sum_(j in J) n_(j,i) >= w_i > 0`, so `i in J`. Thus `I_0 <= J`; because
`|I_0| = |J| = t`, we conclude `J = I_0`. No other summand contributes to `M`.

Combining Sections 5–7, the coefficient of the positive-degree monomial `M` in
the reduced polynomial of `f(R(X))` is nonzero. Thus `f` is nonconstant on
`A`, proving the Lemma.
