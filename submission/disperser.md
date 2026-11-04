# D. A self-contained affine disperser of sublinear dimension

This file proves Theorem D, used by `proof.md` Section 4. It is a complete
reproduction, with full proofs, of the subspace-polynomial argument of
Ben-Sasson and Kopparty ("Affine dispersers from subspace polynomials",
author preprint dated 2010-03-12, Sections 3, 6 and 7; their Theorem 7.3).
Nothing below is cited without proof except the elementary facts of
`fields.md` (F1-F3). Small index typos of the preprint are corrected and
two steps (Lemmas D6.2 and D6.4) are reorganised; the mathematics is theirs.

Throughout, `k >= 3` is a **prime**, `K` is a field with `2^k` elements
(`fields.md` F3), `Tr(a) = sum_{j=0}^{k-1} a^(2^j)` and `sigma(a) = a^2`.
By F3: `a^(2^k) = a` for all `a` in `K`; `sigma` is an `F_2`-linear
automorphism; the only `a` with `a^(2^s) = a` for some `0 < s < k` are
`a in {0,1}`; and `Tr` maps `K` into `{0,1}` and is `F_2`-linear.
"Subspace" means `F_2`-linear subspace and "affine subspace" means a coset
of one. `wt(e)` is the number of ones in the binary expansion of an integer
`e >= 0`. Indices of coefficient sequences `a_0, ..., a_{k-1}` are read
modulo `k`.

## Theorem D

Let `r >= 2` and `1 <= t <= r` be integers with `r + 2 <= k`. For
`1 <= i <= r` put `e_i = r + 1 - i` and `alpha_i = 2^(e_i+1) - 1 =
sum_{l=0}^{e_i} 2^l`. Define `F : K^r -> {0,1}` by

    F(x_1,...,x_r) = Tr( sum_{I subset [r], |I| = t} prod_{i in I} x_i^(alpha_i) ).

If `A` is an affine subspace of `K^r` (viewed as `F_2^(kr)`) with
`dim A > kr/t + kt + r(r+1)`, then `F` is not constant on `A`.

## D1. Reduced polynomials and weights

Let `X_1,...,X_u` be variables. A polynomial over `K` is **reduced** if each
exponent of each variable lies in `[0, 2^k - 1]`. For an exponent `a >= 0`
let `red(a) = 0` if `a = 0`, and otherwise let `red(a)` be the unique integer
in `[1, 2^k - 1]` congruent to `a` modulo `2^k - 1`. For a polynomial `P`,
`red(P)` replaces every exponent `a` of every variable by `red(a)` and
collects terms. Since `x^a = x^(red(a))` for every `x in K` (because
`x^(2^k - 1) = 1` for `x != 0`, and `0^a = 0` for `a >= 1`), `red(P)` and `P`
define the same function on `K^u`. Also `red(P+Q) = red(P) + red(Q)` and
`red(PQ) = red(red(P) red(Q))`, because exponent addition is compatible with
reduction modulo `2^k - 1` and positivity is preserved.

**D1.1 (uniqueness).** If a reduced `P` vanishes on all of `K^u` then `P = 0`.
*Proof.* Induction on `u` (`u = 0` trivial). Write `P = sum_a P_a X_u^a`
with `a < 2^k` and `P_a` reduced in `X_1..X_{u-1}`. For each fixed
`x' in K^(u-1)`, the univariate polynomial `sum_a P_a(x') Y^a` has degree
`< 2^k` and `2^k` roots, so it is zero (F1); thus each `P_a` vanishes
everywhere and is zero by induction. QED.
*Consequence.* If a reduced `P` has a nonzero coefficient at some monomial
other than `1`, then `x -> P(x)` is not constant on `K^u`: otherwise
`P - c` would be a nonzero reduced polynomial vanishing everywhere.

For a monomial `c * prod X_i^(a_i)` (`c != 0`) its `X_i`-weight is
`wt(a_i)`; `deg_i(P)` is the largest `X_i`-weight of a monomial of `P`
with nonzero coefficient (`-infinity` for `P = 0`).

**D1.2.** For `a, b in [0, 2^k - 1]`, `wt(red(a+b)) <= wt(a) + wt(b)`. Hence
`deg_i(red(PQ)) <= deg_i(P) + deg_i(Q)` for reduced `P, Q`.
*Proof.* Binary addition gives `wt(a+b) <= wt(a) + wt(b)` (each carry turns
two ones into one). If `a + b < 2^k` we are done. Otherwise
`a + b = 2^k + c'` with `0 <= c' <= 2^k - 2`, `red(a+b) = c' + 1`, and
`wt(c'+1) <= wt(c') + 1 = wt(a+b)`. The second statement follows by
multiplying monomials. QED.

**D1.3.** Let `b_1,...,b_s in {0,...,k-1}` (repetitions allowed) and
`c = red(2^(b_1) + ... + 2^(b_s))`. Then `wt(c) <= s`, with equality only if
the `b_l` are pairwise distinct, in which case `c = sum_l 2^(b_l)`.
*Proof.* `wt(c) <= s` follows from D1.2 by induction. If the `b_l` are
distinct the sum is below `2^k` and has weight `s`. If `b_1 = b_2 = b`, then
`2^b + 2^b = 2^(b+1)`, which is congruent modulo `2^k - 1` to `2^(b')`
with `b' = b+1 mod k`; replacing the pair by `b'` keeps the residue modulo
`2^k - 1` and the positivity of the sum, hence keeps `c`, and leaves `s - 1`
terms, so `wt(c) <= s - 1`. QED.

## D2. Linearized and subspace polynomials

A **linearized** polynomial is `L(X) = sum_{j=0}^{k-1} a_j X^(2^j)` with
`a_j in K`. Since `(y+z)^2 = y^2 + z^2`, it satisfies `L(X+Y) = L(X) + L(Y)`
as polynomials and defines an `F_2`-linear map `K -> K`.

**D2.1.** Every `F_2`-linear map `T : K -> K` equals (as a function) exactly
one linearized polynomial. *Proof.* Linearized polynomials form an
`F_2`-space of dimension `k * k`, as do `F_2`-linear maps `K -> K`. The map
from the first space to the second is `F_2`-linear and injective, because a
nonzero linearized polynomial has degree at most `2^(k-1) < 2^k` and so does
not vanish on all of `K` (F1). Equal finite dimensions give bijectivity. QED.

**D2.2 (Frobenius twist).** If `L = sum_j a_j X^(2^j)` and `s >= 0`, then
`red(L^(2^s)) = sum_{j=0}^{k-1} (a_{j-s})^(2^s) X^(2^j)`.
*Proof.* Squaring is additive on polynomials in characteristic 2, so
`L^(2^s) = sum_j a_j^(2^s) X^(2^(j+s))`, and `red(2^(j+s)) = 2^((j+s) mod k)`.
QED.

**D2.3 (kernel-subspace polynomial).** For a subspace `V` of `K` of dimension
`d <= k-1`, `P_V(X) = prod_{v in V} (X - v)` is a monic linearized polynomial
of degree `2^d` whose roots in `K` are exactly the elements of `V`.
*Proof.* Induction on `d`; `P_{0} = X`. If `V = V' + F_2 b` with `b` not in
`V'`, then `P_V(X) = P_{V'}(X) P_{V'}(X - b) = P_{V'}(X)(P_{V'}(X) - P_{V'}(b))
= P_{V'}(X)^2 - P_{V'}(b) P_{V'}(X)`, which is linearized (square of a
linearized polynomial is linearized, D2.2 without reduction since the degree
stays below `2^k`) and monic of degree `2^d`. Its roots are `V` by definition.
QED.

**D2.4 (image-subspace polynomial).** Let `L` be a subspace of `K` of
dimension `d`, `1 <= d <= k`. There is a monic linearized polynomial `Q_L`
of degree `2^(k-d)` with `Q_L(K) = L` and with exactly `2^(k-d)` roots in `K`.
Writing `Q_L = sum_{j=0}^{k-1} a_j X^(2^j)`: `a_{k-d} = 1` and `a_j = 0` for
`k-d < j <= k-1`.
*Proof.* If `d = k` take `Q_L = X`. Otherwise `P_L` (D2.3) is an `F_2`-linear
map on `K` with kernel `L`, so `L' = P_L(K)` is a subspace of dimension
`k - d <= k - 1`; put `Q_L = P_{L'}` (D2.3), monic linearized of degree
`2^(k-d)` with root set `L'`. The polynomial `Q_L(P_L(X))` is monic of
degree `2^k` and vanishes on `K`, so `Q_L(P_L(X)) - (X^(2^k) - X)` has degree
`< 2^k` and `2^k` roots, hence `Q_L(P_L(X)) = X^(2^k) - X`. Let
`G(Y) = P_L(Q_L(Y)) - Y^(2^k) + Y`. Then
`G(P_L(X)) = P_L(X^(2^k) - X) - P_L(X)^(2^k) + P_L(X)
= P_L(X^(2^k)) - P_L(X)^(2^k) = sum_j (c_j - c_j^(2^k)) X^(2^(j+k)) = 0`,
where `c_j` are the coefficients of `P_L` and `c^(2^k) = c`. As `P_L` is
nonconstant, `G(P_L(X)) = 0` forces `G = 0`. Hence `P_L(Q_L(y)) = 0`, i.e.
`Q_L(y) in L`, for all `y in K`. The linear map `Q_L` has kernel `L'` of
size `2^(k-d)`, so its image has size `2^d = |L|`, and `Q_L(K) = L`. QED.

## D3. Structural lemma (Ben-Sasson-Kopparty Lemma 3.10, prime `k`)

Let `L` be a subspace of `K` of dimension `d >= 2` and
`Q_L = sum_j a_j X^(2^j)` its image-subspace polynomial (D2.4). Then:

(a) every `d` cyclically consecutive indices `{c, c+1, ..., c+d-1} mod k`
contain an index `j` with `a_j != 0`; in particular `a_0 != 0`;

(b) every `d-1` cyclically consecutive indices `{c, ..., c+d-2} mod k` with
`c != k-d+1 (mod k)` contain an index `j` with `a_j != 0`.

*Proof.* `Q_L` has exactly `2^(k-d)` roots in `K` (D2.4); call this root set
`Z`. For `s >= 0` put `Q' = red(Q_L^(2^s)) = sum_j a_{j-s}^(2^s) X^(2^j)`
(D2.2). `Q'` is nonzero (Frobenius is injective) and vanishes on `Z`
because `Q'(y) = Q_L(y)^(2^s)`.

(a) Suppose `a_j = 0` for `j in {c,...,c+d-1}`. Choose `s = k-d-c mod k`.
For `j in [k-d, k-1]` we have `j - s in {c,...,c+d-1}`, so the coefficient
of `X^(2^j)` in `Q'` is zero. If `d = k` this says `Q' = 0`, a contradiction.
If `d < k` it gives `deg Q' <= 2^(k-d-1) < 2^(k-d) = |Z|`, so the nonzero
`Q'` cannot vanish on `Z` (F1), a contradiction. For `a_0`: positions `k-d+1,...,k-1` are zero (D2.4), so by (a) applied
to `{k-d+1,...,k-1,0}` we get `a_0 != 0`.

(b) Suppose `a_j = 0` for `j in {c,...,c+d-2}` with `c != k-d+1`. Choose
`s = k-d+1-c mod k`; then `0 < s < k`. For `j in [k-d+1, k-1]`, `j - s` lies
in the zero run, so `deg Q' <= 2^(k-d)`. Let `c'` be the coefficient of
`X^(2^(k-d))` in `Q'`. Then `Q' - c' Q_L` has degree `< 2^(k-d)` and vanishes
on `Z`, so `Q' = c' Q_L` as polynomials, and `c' != 0` because `Q' != 0`.
Hence `y^(2^s) = c' y` for every `y in L = Q_L(K)`. Take a nonzero
`y_0 in L`; for any nonzero `y in L`, `(y/y_0)^(2^s) = y/y_0`, so by F3
`y/y_0 in {0,1}`, i.e. `y = y_0`. Thus `|L| <= 2`, contradicting `d >= 2`.
QED.

## D4. Coefficients of products (Ben-Sasson-Kopparty Lemma 3.13)

Let `Q_1(X_1),...,Q_u(X_u)` be linearized with coefficients
`Q_i = sum_j a_{i,j} X_i^(2^j)`. Put `M_i[rho, b] = (a_{i, b - rho})^(2^rho)`
for `rho, b in Z_k` (well defined because `a^(2^k) = a`); by D2.2 this is the coefficient of `X_i^(2^b)` in
`red(Q_i^(2^rho))`. Let `e_i >= 0` with `e_i + 1 <= k`, `E_i = {0,...,e_i}`,
`alpha_i = sum_{l in E_i} 2^l`, let `B_i` be subsets of `Z_k` with
`|B_i| = e_i + 1`, and `beta_i = sum_{b in B_i} 2^b`.

**Lemma D4.** The coefficient of `prod_i X_i^(beta_i)` in
`red( Tr( prod_{i=1}^u Q_i(X_i)^(alpha_i) ) )` equals

    sum_{w=0}^{k-1} prod_{i=1}^u Perm( M_i[E_i + w, B_i] ),

where `M_i[E_i + w, B_i]` has rows `l + w (mod k)`, `l in E_i`, and columns
`B_i`, and `Perm` is the permanent.

*Proof.* `Tr(P) = sum_w P^(2^w)`, and
`(prod_i Q_i^(alpha_i))^(2^w) = prod_i prod_{l in E_i} Q_i^(2^(l+w))`.
After reduction each factor is `sum_b M_i[l+w, b] X_i^(2^b)`. Expanding,
a term selects one `b_{i,l}` per factor and yields
`prod_i X_i^(red(sum_l 2^(b_{i,l})))`. By D1.3 the `X_i`-exponent has weight
`e_i + 1` only if the `b_{i,l}`, `l in E_i`, are distinct, and then it equals
`beta_i` exactly when `{b_{i,l}} = B_i` (binary expansions are unique). So
the contributions to our monomial come from the bijections `E_i -> B_i`, one
per `i`, and their sum is the product of the permanents. Summing over `w`
gives the formula. QED.

Note `M_i[rho + w, b] = (M_i[rho, b - w])^(2^w)`, so the entry of
`M_i[E_i + w, B_i]` in row `l + w` and column `b` is zero iff
`a_{i, b - w - l} = 0`.

**Criterion D4'.** Suppose (i) for every `i`, the columns of `M_i[E_i, B_i]`
can be ordered so that the matrix is lower triangular with nonzero diagonal;
and (ii) for every `w in {1,...,k-1}` some `i` has a column `b in B_i` with
`a_{i, b - w - l} = 0` for all `l in E_i`. Then the coefficient in Lemma D4
is nonzero. *Proof.* For `w >= 1` the factor of that `i` has a zero column,
so its permanent vanishes. For `w = 0`, a lower-triangular matrix has
permanent equal to the product of its diagonal entries (every other
bijection uses an entry above the diagonal), which is nonzero, and the
product over `i` of nonzero elements of the field `K` is nonzero. QED.

## D5. Product theorem (Ben-Sasson-Kopparty Theorem 6.1)

**Theorem D5.** Let `Y_1,...,Y_u` be subspaces of `K` with dimensions
`d_i >= 3`, let `1 <= e_i <= d_i - 2`, `alpha_i = 2^(e_i + 1) - 1`, and assume

    sum_{i=1}^u (d_i - e_i - 1) >= k - max_i d_i + 2.

Let `Q_i` be the image-subspace polynomial of `Y_i`. Then
`g'(X) = red( Tr( prod_i Q_i(X_i)^(alpha_i) ) )` has a nonzero coefficient at
some monomial `prod_i X_i^(beta_i)` with `wt(beta_i) = e_i + 1` for all `i`.

*Proof.* Relabel so that `d_1 = max_i d_i`. Write `a_{i,j}` for the
coefficients of `Q_i`, `delta_i = d_i - e_i - 1 >= 1` and
`Zr_i = {k-d_i+1, ..., k-1}`; by D2.4 `a_{i,j} = 0` for `j in Zr_i`,
`a_{i,k-d_i} = 1`, and by D3(a) `a_{i,0} != 0`. Put `ell_1 = 0` and
`ell_{i+1} = ell_i + delta_i` (ordinary integers). Since
`ell_{u+1} >= k - d_1 + 2`, let `t'` be the least `i` with
`ell_{i+1} >= k - d_1 + 1`. For `i <= t'`, `ell_i <= k - d_1 <= k - d_i`.

*Choice of `B_i`.* For `i > t'`, and for `i <= t'` with `a_{i,ell_i} != 0`,
let (using `ell := 0` when `i > t'`)

    B_i = {ell} U {k-d_i+1, ..., k-d_i+e_i}.

For `i <= t'` with `a_{i,ell_i} = 0` (so `0 < ell_i < k - d_i`), let
`jh_i = max{ j < ell_i : a_{i,j} != 0 }` (exists as `a_{i,0} != 0`),
`j_i = min{ j > ell_i : a_{i,j} != 0 }` (exists, `<= k-d_i`),
`q_i = j_i - jh_i - 1 >= 1`, `s_i = min(q_i, e_i)` and

    B_i = {j_i} U {jh_i+1, ..., jh_i+s_i} U {k-d_i+s_i+1, ..., k-d_i+e_i}.

In both cases the listed elements are distinct and lie in `[0, k-1]`:
`ell <= k-d_i`, `jh_i + s_i < j_i <= k - d_i`, and `k-d_i+e_i <= k-2`. So
`|B_i| = e_i + 1`.

*Condition (i) of D4'.* Row `rho in E_i`, column `b` has entry zero iff
`a_{i,b-rho} = 0`. First form: order columns `ell, k-d_i+1, ..., k-d_i+e_i`.
For column `k-d_i+sig` (`1 <= sig <= e_i`) and row `rho < sig`, the index
`k-d_i+sig-rho` lies in `[k-d_i+1, k-1]` (as `sig - rho <= e_i <= d_i - 2`),
so the entry is zero; the diagonal entries are `a_{i,ell}` and `a_{i,k-d_i}`,
nonzero. Second form: order columns `j_i, jh_i+1, ..., jh_i+s_i,
k-d_i+s_i+1, ..., k-d_i+e_i`. Entries above the diagonal have index in
`[jh_i+1, jh_i+s_i]` (inside the zero run `jh_i+1..j_i-1`, since `s_i <= q_i`)
or in `Zr_i`; diagonal entries are `a_{i,j_i}`, `a_{i,jh_i}`, `a_{i,k-d_i}`,
all nonzero. So (i) holds.

*Condition (ii) of D4'.* For `b in B_i` and an integer
`w in [b+1, b+delta_i]`, write `w = b + tau`, `1 <= tau <= delta_i`. For
`l in E_i`, `b - w - l = -(tau + l)` with `1 <= tau + l <= delta_i + e_i =
d_i - 1`, so `b - w - l mod k` lies in `Zr_i` and `a_{i,b-w-l} = 0`. Hence it
suffices to show that `[1, k-1]` is contained in
`U = union_i union_{b in B_i} [b+1, b+delta_i]`.

Block 1: `ell_1 = 0` and `a_{1,0} != 0`, so `B_1 = {0, k-d_1+1, ..., k-d_1+e_1}`
and `U` contains `[1, delta_1] = [ell_1+1, ell_2]` and
`union_{sig=1}^{e_1} [k-d_1+sig+1, k-d_1+sig+delta_1] = [k-d_1+2, k-1]`
(as `e_1 + delta_1 = d_1 - 1`).

Blocks `1 < i <= t'`: we show `[ell_i + 1, ell_{i+1}] subset U`.
If `a_{i,ell_i} != 0` then `ell_i in B_i` and `[ell_i+1, ell_i+delta_i]` is
one of the intervals. If `a_{i,ell_i} = 0` and `q_i < e_i`, then
`B_i` contains `[jh_i+1, jh_i+q_i] U {j_i} = [jh_i+1, j_i]`, which contains
`ell_i`; done as before. If `a_{i,ell_i} = 0` and `q_i >= e_i`, then the zero
run `jh_i+1, ..., j_i-1` starts at `jh_i + 1 <= ell_i <= k - d_1 <= k - d_i`,
which differs from `k - d_i + 1` modulo `k` (both lie in `[1, k-1]`). By
D3(b), applied to `Y_i` (dimension `d_i >= 2`), this run has length at most
`d_i - 2`, i.e. `j_i <= jh_i + d_i - 1`. The intervals for
`b = jh_i+1, ..., jh_i+e_i` cover `[jh_i+2, jh_i+e_i+delta_i] =
[jh_i+2, jh_i+d_i-1]`, which contains `[jh_i+2, j_i]`; the interval for
`b = j_i` is `[j_i+1, j_i+delta_i]`. Together they cover
`[jh_i+2, j_i+delta_i]`, which contains `[ell_i+1, ell_i+delta_i]` because
`jh_i < ell_i < j_i`.

Hence `U` contains `[1, ell_{t'+1}]`, with `ell_{t'+1} >= k-d_1+1`, and
`[k-d_1+2, k-1]`; so `[1, k-1] subset U`. By Criterion D4' the coefficient
of `prod_i X_i^(beta_i)`, `beta_i = sum_{b in B_i} 2^b`, is nonzero, and
`wt(beta_i) = |B_i| = e_i + 1`. QED.

## D6. Proof of Theorem D

**D6.1 (echelon decomposition).** Let `A` be an affine subspace of `K^r`
and `gamma in A`. There are subspaces `Y_1,...,Y_r` of `K` and `F_2`-linear
maps `sigma_{ij} : Y_j -> K` (`j < i`) such that
`(y_1,...,y_r) -> (gamma_i + y_i + sum_{j<i} sigma_{ij}(y_j))_{i=1..r}`
is a bijection from `Y_1 x ... x Y_r` onto `A`; in particular
`dim A = sum_i dim Y_i`.
*Proof.* Let `V = A - gamma`, `V_i = {v in V : v_1 = ... = v_{i-1} = 0}`,
`Y_i = {v_i : v in V_i}` and choose an `F_2`-linear section
`s_i : Y_i -> V_i` with `(s_i(y))_i = y`. Put `sigma_{ij}(y) = (s_j(y))_i`.
The map `(y_j) -> sum_j s_j(y_j)` has `i`-th block
`y_i + sum_{j<i} sigma_{ij}(y_j)`, because `s_j(y_j)` vanishes in blocks
before `j`. It is injective (if the sum is 0, block 1 gives `y_1 = 0`, then
block 2 gives `y_2 = 0`, and so on) and surjective onto `V` (given `v in V`,
take `y_1 = v_1`, then `v - s_1(y_1) in V_2`, and repeat). QED.

*Setup.* Let `A` be as in Theorem D, `d = dim A`, and take the decomposition
D6.1 with `d_i = dim Y_i`, `sum_i d_i = d`. Let `Q_i` be the image-subspace
polynomial of `Y_i` (D2.4; if `Y_i = 0` put `Q_i = 0`). For `j < i` let
`Q_{ij}` be the linearized polynomial with `Q_{ij}(x) = sigma_{ij}(Q_j(x))`
for `x in K` (D2.1; zero if `Y_j = 0`). Put
`R_i(X_1,...,X_i) = Q_i(X_i) + sum_{j<i} Q_{ij}(X_j) + gamma_i` and
`H = (R_1,...,R_r) : K^r -> K^r`. By D6.1 and `Q_j(K) = Y_j`, `H(K^r) = A`.
Let

    h = red( F(R_1, ..., R_r) ) = sum_{|J| = t} red( Tr( prod_{j in J} R_j^(alpha_j) ) ).

Then `h(K^r) = F(A)`. By D1.1 it suffices to find a monomial of positive
degree with nonzero coefficient in `h`.

*Choice of `I_0`.* Let `S = {i : d_i >= r + 2}`. Since `d_i <= k` and
`sum_{i not in S} d_i <= r(r+1)`, we get `sum_{i in S} d_i > kt`, so
`|S| > t`. Moreover `sum_{i in S} (d_i - r - 1) >= sum_{i=1}^r (d_i - r - 1)
= d - r(r+1) > kr/t + kt`. Let `I_0` consist of `t` indices of `S` with the
largest values of `d_i - r - 1` (all positive). The average over `I_0` is at
least the average over `S`, and `|S| <= r`, so

    sum_{i in I_0} (d_i - r - 1) >= (t/|S|) * (kr/t + kt) >= (t/r)(kr/t) = k.

*The monomial.* For `i in I_0`, `e_i = r + 1 - i` satisfies
`1 <= e_i <= r <= d_i - 2`, and `d_i >= r + 2 >= 4`. Also
`sum_{i in I_0} (d_i - e_i - 1) >= sum_{i in I_0} (d_i - r - 1) >= k
>= k - max d_i + 2`. Theorem D5 for the family `(Y_i)_{i in I_0}` gives
`M = prod_{i in I_0} X_i^(beta_i)`, `wt(beta_i) = e_i + 1`, with nonzero
coefficient in `g' = red(Tr(prod_{i in I_0} Q_i(X_i)^(alpha_i)))`. `M` has
positive degree.

**D6.2.** The coefficient of `M` in `g = red(Tr(prod_{i in I_0}
R_i^(alpha_i)))` equals its coefficient in `g'`.
*Proof.* Put `P_i = R_i - Q_i(X_i)`, a polynomial in `X_1..X_{i-1}` only.
For each `w`, `(prod_{i in I_0} R_i^(alpha_i))^(2^w) = prod_{i in I_0}
prod_{l=0}^{e_i} (Q_i(X_i)^(2^(l+w)) + P_i^(2^(l+w)))`. Expand this product
into terms indexed by sets `L_i subset E_i` (the factors where the `P_i`
part is chosen). The terms with all `L_i` empty sum to
`(prod_{i in I_0} Q_i(X_i)^(alpha_i))^(2^w)`, i.e. to the `w`-th summand of
`g'`. Consider a term with some `L_i` nonempty and let `i*` be the largest
such `i`. In this term the variable `X_{i*}` occurs only in the factors of
block `i*`: blocks `i > i*` of `I_0` contribute only powers of `Q_i(X_i)`,
and blocks `i < i*` only involve `X_1..X_i`. In block `i*` at most
`e_{i*} + 1 - |L_{i*}| <= e_{i*}` factors involve `X_{i*}`, each of
`X_{i*}`-weight at most 1 after reduction (D2.2). By D1.2 every monomial of
the reduced term has `X_{i*}`-weight at most `e_{i*} < wt(beta_{i*})`, so it
is not `M`. Summing over `w` gives the claim. QED.

**D6.3 (combinatorial lemma; Ben-Sasson-Kopparty Claim 7.4).** Let `I_0` and
`J` be `t`-subsets of `[r]`. For each `j in J` let
`lambda_j : {0,...,e_j} -> {0,1,...,j}`, and for `i >= 1` let `n_i` be the
number of pairs `(j,l)` with `lambda_j(l) = i`. If `n_i >= r + 2 - i` for all
`i in I_0`, then `J = I_0`.
*Proof.* Induction on `t`; `t = 0` is trivial. For `p >= 0` let
`g_p(0) = 0` and `g_p(y) = y^p` for `y >= 1`, and
`Delta_p = sum_{i in I_0} i^p - sum_{j in J} j^p`; `Delta_0 = 0`. Since
`g_p(lambda_j(l)) <= j^p` and the domain of `lambda_j` has `r + 2 - j`
elements,

    sum_{j in J} (r+2-j) j^p >= sum_{j,l} g_p(lambda_j(l)) >= sum_{i in I_0} n_i i^p
                             >= sum_{i in I_0} (r+2-i) i^p,

i.e. `Delta_{p+1} >= (r+2) Delta_p`. Hence `Delta_p >= 0` for all `p`. If
`max J > max I_0` then `sum_J j^p >= (max J)^p > t (max I_0)^p >=
sum_{I_0} i^p` for large `p`, a contradiction; so `max J <= max I_0 =: c`.
As `n_c >= r + 2 - c >= 1`, some `lambda_j` takes the value `c`, forcing
`j >= c`; so `c = max J`. Only `j = c` can take the value `c`, and its domain
has exactly `r + 2 - c` elements, so `lambda_c` is constantly `c`. Deleting
`c` from `I_0` and `J` (and `lambda_c`) leaves the counts `n_i`,
`i in I_0 - {c}`, unchanged, and induction gives `J - {c} = I_0 - {c}`. QED.

**D6.4.** For `J != I_0`, the coefficient of `M` in
`red(Tr(prod_{j in J} R_j^(alpha_j)))` is zero.
*Proof.* For each `w`, `(prod_{j in J} R_j^(alpha_j))^(2^w) =
prod_{j in J} prod_{l=0}^{e_j} R_j^(2^(l+w))`, and `red(R_j^(2^s))` is a sum
of a constant and monomials `c X_i^(2^b)` with `i <= j` (D2.2). Expanding,
each term chooses for every `(j,l)` either the constant (`lambda_j(l) = 0`)
or a monomial in `X_{lambda_j(l)}`, with `lambda_j(l) <= j`. By D1.2 a
monomial of the reduced term has `X_i`-weight at most `n_i`. If `M` occurred
in some term then `n_i >= wt(beta_i) = e_i + 1 = r + 2 - i` for all
`i in I_0`, and D6.3 would give `J = I_0`. QED.

*Conclusion.* By D6.2, D6.4 and the choice of `M`, the coefficient of `M` in
`h` equals its (nonzero) coefficient in `g'`. `M` has positive degree, so by
D1.1 `h` is not constant on `K^r`, i.e. `F` is not constant on `A`. This
proves Theorem D. QED.
