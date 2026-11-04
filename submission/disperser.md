# The algebraic disperser: a complete proof of the needed BK statement

**Disperser lemma (D).** Let `k` be a prime, let `K = F_{2^k}` be the field
with `2^k` elements, and let `1 <= t < r` with

    k t >= r (r + 1).                                          (D1)

For each `i in [r]` put `e_i = r + 1 - i` and `w_i = e_i + 1 = r + 2 - i`.
Let `pi : K -> F_2` be a nonzero `F_2`-linear map. Define

    f(a_1,...,a_r) = pi ( sum_{|I| = t} prod_{i in I} a_i^{2^{w_i} - 1} ) .  (D2)

Then `f` is nonconstant on every affine subspace of `K^r` of dimension
strictly greater than

    k r / t + k t + r (r + 1).                                 (D3)

This is the binary instance of Ben-Sasson and Kopparty, *Affine Dispersers from
Subspace Polynomials*, author preprint of 12 March 2010, Theorem 7.3, together
with the reduction from their Theorem 2.2 (their equation (4)). It is
reproduced here in full, specialised to `F_2` and to prime field degree `k`,
because `rules.md` requires unfamiliar or uncertain external premises to be
supported by argument reproduced and validated in the supplied evidence.

Everything in Sections 1-7 below is proved here. Section 8 records the
correspondence with the source, and Section 9 records the scope.

Throughout, `q = 2^k` is the size of `K` and all dimensions are over `F_2`.

## 1. Finite-field and polynomial preliminaries

**(1.1) Unique reduced representatives.** Every function `K^r -> K` has a
unique polynomial representative in which the degree in each variable is
less than `q`. Indeed, reduce modulo `X_i^q - X_i`; the reduced polynomial is
the unique one of degree `< q` in each variable that agrees with the function
on `K^r`. Uniqueness follows by induction on the number of variables using the
univariate root bound: if two such polynomials differ, take the difference; it
is not identically zero, and viewing it as a polynomial in `X_1` over
`K[X_2,...,X_r]` it has degree `< q` in `X_1`, so it has at most `q-1` roots
for each fixed value of the other variables, whereas it vanishes on all `q^{r-1}`
such values; by induction it is zero.

Consequently, a nonzero coefficient of *positive degree* in the reduced
representative proves that the function is nonconstant.

**(1.2) Linearized polynomials.** A **linearized** polynomial is one of the
form `sum_{j=0}^{k-1} a_j X^{2^j}` with `a_j in K`. These are exactly the
`F_2`-linear maps `K -> K`: each is `F_2`-linear by the Frobenius identity
`(a+b)^{2} = a^2 + b^2`. The space of such coefficient vectors has
`F_2`-dimension `k`, and the space of `F_2`-linear maps `K -> K` also has
dimension `k`. The map from coefficients to functions is injective by the root
bound (a nonzero linearized polynomial of degree `< q` has at most `q-1` roots,
and a linearized polynomial vanishing on all of `K` is zero), hence it is
surjective: every `F_2`-linear map `K -> K` has a unique linearized
representation.

**(1.3) Trace and the dual.** The absolute trace `Tr : K -> F_2` is
`Tr(z) = sum_{h=0}^{k-1} z^{2^h}`. It is not the zero function: the
corresponding polynomial is nonzero and has degree `2^{k-1} < q`, so it cannot
vanish on all of `K` (root bound). For `mu != 0` the functional
`z |-> Tr(mu z)` is nonzero, since multiplication by `mu` is a bijection of
`K`. The map `mu |-> (z |-> Tr(mu z))` is injective between two `k`-dimensional
`F_2`-spaces, hence bijective. Therefore

    every nonzero F_2-linear pi : K -> F_2 has the form
    pi(z) = Tr(mu z) for a unique nonzero mu in K.               (1)

**(1.4) Binary degree.** For `b` in `{0,...,q-1}` let `wt(b)` be the number of
`1`s in the binary expansion of `b`. Multiplying and reducing cannot increase
the total number of power-of-two factors: replace two equal powers by one
doubled power, with cyclic wrap at `k` because `X^{2^k} = X` in the reduced
range. Each such carry reduces the number of factors by one. Consequently an
exponent formed from `h` power-of-two factors has `wt <= h`; if `h < k` and
the final weight is exactly `h`, there were no carries, i.e. the chosen
residues modulo `k` were all distinct. Taking a Frobenius power cyclically
shifts residues, so it preserves the weight. This also proves the
subadditivity of binary degree, `wt(b + c mod (q-1)) <= wt(b) + wt(c)`, valid
after reduction modulo `X^q - X`.

## 2. Image polynomials and their zero patterns

**(2.1) Image-subspace polynomials.** For a linear subspace `V` of `K` of
dimension `d` define

    P_V(X) = prod_{v in V} (X - v).

This is a monic linearized polynomial of degree `2^d`. For an induction proof,
start with `P_{0}(X) = X`; adjoining a vector `b` to `U` gives

    P_{U+<b>}(X) = P_U(X) P_U(X + b)
                 = P_U(X)^2 + P_U(b) P_U(X) ,

which is linearized. Its kernel on `K` is exactly `V`, so its image `V'` has
dimension `k - d`.

Let `Q = P_{V'}`. Then `Q(P_V(X))` is monic of degree `q` and vanishes on `K`,
hence equals `X^q - X`. Also `P_V(Q(X)) = X^q - X`: composing their
difference with `P_V` and using `P_V(X^q - X) = P_V(X)^q - P_V(X)` shows the
difference is zero as a polynomial, since a nonzero polynomial stays nonzero
after a nonconstant polynomial substitution. It follows that `Q(K)` is
contained in `V`; since the kernel of `Q` has dimension `k - d`, the image of
`Q` has dimension `d`, so `Q(K) = V`.

Thus for every `d`-dimensional linear subspace `V` of `K` there is a monic
**image polynomial** `Q_V` of degree `2^{k-d}` with `2^{k-d}` distinct roots
in `K`, mapping `K` onto `V`. Write

    Q_V = sum_{j=0}^{k-1} a_j X^{2^j} ,                            (2)

with coefficients read cyclically modulo `k`.

**(2.2) The zero-run structure.** There is no cyclic run of `d` zero
coefficients in (2). Otherwise a Frobenius power of `Q_V`, reduced modulo
`X^q - X`, moves that run to the final `d` positions; the resulting nonzero
polynomial would have degree at most `2^{k-d-1}` while still having all
`2^{k-d}` roots of `Q_V`, contradicting the root bound. In particular
`a_0 != 0` and `a_{k-d} != 0`; the latter is also monicity. The interval
`k-d+1, ..., k-1` is therefore a zero run of length `d-1`.

**(2.3) Prime degree excludes other long runs.** Suppose `k` is prime, `d > 1`,
and there is a cyclic run of `d-1` zero coefficients other than the final one.
Shift it to the final positions by a Frobenius power `2^h` with `0 < h < k`.
The shifted polynomial has degree at most `2^{k-d}`, has the same roots as
`Q_V`, and so is `c Q_V` for some `c != 0` (subtract leading terms and use the
root bound). Hence every `z in V` satisfies `z^{2^h} = c z`.

Choose `z_0 != 0` in `V`. Then `z / z_0` is fixed by `Frob^h` and by `Frob^k`.
Since `gcd(h,k) = 1` (as `k` is prime and `0 < h < k`), Bezout gives integers
`u,v` with `u h + v k = 1`, and `Frob` is an invertible automorphism, so
`z / z_0` is fixed by `Frob` itself: `z / z_0 in F_2`. That forces
`dim(V) <= 1`, contradicting `d > 1`.

Therefore:

> **(Zero-run lemma.)** If `k` is prime and `d > 1`, then every cyclic run of
> zero coefficients of `Q_V` other than the final one has length at most
> `d - 2`.

This is exactly the prime-degree subfield exclusion that the source obtains
from its structural lemma; we prove it directly, so that no subfield-containment
hypothesis has to be carried through the rest of the argument.

## 3. The coefficient assertion for independent subspaces

**(3.1) Statement.** Let `V_1, ..., V_t` be linear subspaces of `K` of
dimensions `d_i > 1`. Let `1 <= e_i < d_i - 1`, put `h_i = e_i + 1` and
`delta_i = d_i - h_i > 0`, and suppose

    sum_i delta_i > k - max_i d_i + 1 .                            (3)

Then the reduced polynomial

    Tr( mu prod_i Q_{V_i}(X_i)^{2^{h_i} - 1} )                  (4)

has a nonzero coefficient at a monomial `prod_i X_i^{b_i}` with
`wt(b_i) = h_i` for every `i`.

Because `h_i < d_i <= k`, we have `h_i < k`, and the no-carry observation of
Section 1.4 applies.

**(3.2) The associated matrix of a linearized polynomial.** For
`Q = sum_{j=0}^{k-1} a_j X^{2^j}` define the `k x k` matrix `M_Q` over `K` by

    M_Q[u,v] = a_{v-u}^{2^u},        indices modulo k.             (5)

This is the source's associated matrix (its Definition 3.11) with `p = 2`:
its `(u,v)` entry is `(a_{v-u})^{p^u}`. Let `M_Q^0` be its **zero-one
indicator**: `M_Q^0[u,v] = 1` if `a_{v-u} != 0` and `0` otherwise. Write
`M_i = M_{Q_i}` and `M_i^0 = M_{Q_i}^0`.

**(3.3) The coefficient formula.** Let `A_i = {0,...,e_i}` and let `B_i` be a
set of `h_i` column indices with `b_i = sum_{v in B_i} 2^v`. Then the
coefficient of `prod_i X_i^{b_i}` in (4) is

    sum_{rho=0}^{k-1} mu^{2^rho} prod_i Perm( M_i[A_i, B_i - rho] )^{2^rho} ,  (6)

where `Perm` denotes the permanent of the indicated square submatrix and
`B_i - rho` is computed modulo `k`.

*Justification.* The exponent `2^{h_i} - 1` expands into `h_i` Frobenius
factors `X^{2^{v}}` for `v in B_i` (in the no-carry regime). A monomial of
binary weight `h_i` in `X_i` must use each of the `h_i` chosen powers exactly
once, since any repeated power would carry and lower the binary degree. The
sum over the bijections from rows to chosen columns is exactly the permanent.
Monomials involving constants or carries have smaller binary degree and so
cannot contribute to `X_i^{b_i}`. This is the source's Lemma 3.13 specialised
to `p = 2`; its equation (8) reads

    c_M = sum_{rho=0}^{k-1} mu^{2^rho} prod_i Perm( M_i[A_i, B_i - rho] )^{2^rho} ,

which is (6).

**(3.4) Two conditions that make (6) nonzero.** Following the source's
Remark 3.14, it suffices to produce column sets `B_i` with

  * **(T1)** `M_i^0[A_i, B_i]` is, up to reordering of rows and columns,
    triangular with every diagonal entry nonzero; and

  * **(T2)** for every `rho in {1,...,k-1}` there is an `i` such that
    `M_i^0[A_i, B_i - rho]` contains an all-zero column.

Indeed, (T1) makes the `rho = 0` summand of (6) a product of permanents of
triangular matrices with nonzero diagonal, hence nonzero; and (T2) makes every
other summand zero, since one matrix in the product then has a zero
permanent.

**(3.5) Choosing the column sets.** Reorder the sources so that
`d_1 = max_i d_i`. Following the source's construction, set

    ell_1 = 0 ,    ell_{i+1} = ell_i + d_i - (e_i + 1)   (mod k),   (7)

so that `ell_i = sum_{j<i} ( d_j - (e_j+1) )  (mod k)`. Let `t_0` be the
least index with `sum_{i <= t_0} ( d_i - (e_i+1) ) >= k - d_1 + 1`; it
exists by (3). Then

    sum_{i <= t_0} ( d_i - (e_i+1) ) < k ,                        (8)

because the sum just before `t_0` is `< k - d_1 + 1` and the last summand is
`d_{t_0} - (e_{t_0}+1) <= d_1 - 1`, so the total is `< k`. In particular
there is no wraparound in (7) for `i <= t_0`, and `ell_i <= k - d_1 <= k - d_i`
for `i <= t_0`.

For each source abbreviate `d = d_i`, `e = e_i`, `ell = ell_i` and
`a_j = a_{i,j}`.

*Case `a_ell != 0`.* Take

    B_i = { ell } u { k - d + 1 , k - d + 2 , ... , k - d + e } .    (9)

*Case `a_ell = 0`.* Let `j_i` be the least index `> ell` with `a_{j_i} != 0`
and `jb_i` the greatest index `< ell` with `a_{jb_i} != 0`. Let
`q_i = j_i - jb_i - 1` be the length of the zero-coefficient interval between
them, and `s_i = min(q_i, e)`. Take

    B_i = { j_i } u { jb_i + 1 , ... , jb_i + s_i }
              u { k - d + s_i + 1 , ... , k - d + e } .           (10)

(The last block is empty if `s_i = e`.)

Note that in this case the zero gap containing `ell` is *not* the final one,
because `ell <= k - d_1 <= k - d` by (8) and the final gap is
`k-d+1,...,k-1`. By the zero-run lemma of Section 2.3 that gap has length at
most `d-2`, so `q_i <= d-2`; in particular `q_i >= 1` (there is a nonzero
coefficient strictly between `jb_i` and `j_i`, since `a_{k-d} != 0` and
`k-d >= ell`).

**(3.6) Verification of (T1).** We check that `M_i^0[A_i, B_i]` is lower
triangular with nonzero diagonal.

*Case `a_ell != 0`.* The entries above the diagonal of `M_i^0[A_i, B_i]` are
`M_i^0[u,v] = a_{v-u}^0` with `u < v`. Since `v` ranges over
`{k-d+1,...,k-d+e}` for the columns beyond the first, we have
`v - u >= k - d + 1 - u >= k - d + 1 - e + 1`; and since `e < d`,
`v - u <= k - 1`. Thus `v - u` lies in the final zero gap
`{k-d+1,...,k-1}` whenever `v != ell`, so those entries are zero. The
diagonal entries are: the top-left one is `a_ell^0 = 1`, and the remaining
ones are `a_{k-d}^0 = 1`, which is nonzero because `a_{k-d} != 0` (Section
2.2). Hence `M_i^0[A_i,B_i]` is lower triangular with nonzero diagonal.

*Case `a_ell = 0`.* Write the uppermost-left `(s_i + 1) x (s_i + 1)`
submatrix. Its first column is indexed by `j_i` and the rest by
`jb_i + 1, ..., jb_i + s_i`; the first row is indexed by `0` and the rest by
`1,...,s_i`. Its entries above the diagonal are `a_{jb_i + t}^0` for
`t = 1,...,s_i`, which are zero by the definition of `s_i` (the coefficients
between `jb_i` and `j_i` are zero, and `s_i <= q_i`). Its diagonal entries are
`a_{j_i}^0 = 1` and `a_{jb_i}^0 = 1`, both nonzero by definition. The last
`e - s_i` columns, if they exist, are the same as in case (9) and are again
handled by the final zero gap. So `M_i^0[A_i,B_i]` is lower triangular with
nonzero diagonal.

This proves (T1).

**(3.7) Verification of (T2).** We follow the source's Proposition 6.4 and
Claim 6.3. For a column `k` of `B_i`, the entries of the column `k - rho` of
`M_i^0[A_i, {k} - rho]` are

    a_{k-rho}^0 , a_{k-rho-1}^0 , ... , a_{k-rho-e}^0 .

If

    [ k - (rho + e) , k - rho ]  (mod k)  is contained in  [ k_i - d_i + 1 , k_i - 1 ] ,  (11)

with `d_i` the dimension of `V_i`, then that column is all-zero. Hence:

> **(Zero-column rule.)** Whenever `k in B_i` and

>       rho in [ k + 1 , k + d_i - (e_i + 1) ]   (mod k) ,          (12)

> the matrix `M_i^0[A_i, B_i - rho]` contains an all-zero column.

(For (11) to hold we need the `e_i + 1` indices `k-rho, ..., k-rho-e_i` to
lie in the final zero gap `k-d_i+1,...,k_i-1`, i.e. `rho` in
`[k+1, k+d_i-e_i-1]`, which is (12).)

We now show that the intervals (12) cover every `rho in {1,...,k-1}`.

*Source 1.* Here `ell_1 = 0`, so by (9) `B_1 = {0} u [k-d_1+1, k-d_1+e_1]`.
The column `0` covers `rho in [1, d_1 - (e_1+1)]`, and the columns
`k-d_1+1, ..., k-d_1+e_1` together cover `rho in [k-d_1+2, k-1]`. So source
`1` covers

    [1, d_1 - (e_1+1)]  u  [k-d_1+2, k-1] .                     (13)

*Sources `1 < i <= t_0`.* We claim that

    union_{k in B_i} [ k+1 , k + d_i - (e_i+1) ]  contains
    [ ell_i + 1 , ell_{i+1} ] ,                                  (14)

where `ell_{i+1} = ell_i + d_i - (e_i+1)` by (7).

If `a_{ell_i} != 0`, then `ell_i in B_i` by (9), and the column `ell_i`
covers `[ell_i+1, ell_i + d_i - (e_i+1)] = [ell_i+1, ell_{i+1}]` directly.

If `a_{ell_i} = 0` and `q_i < e_i`, then (10) gives
`B_i = {j_i} u [jb_i+1, jb_i+q_i] u [k-d_i+q_i+1, k-d_i+e_i]`, and
`j_i = jb_i + q_i + 1`. Hence `ell_i in [jb_i+1, jb_i+q_i] subset B_i`, and
the column `ell_i` covers `[ell_i+1, ell_{i+1}]` as before.

If `a_{ell_i} = 0` and `q_i >= e_i`, then (10) gives
`B_i = {j_i} u [jb_i+1, jb_i+e_i]` with `j_i = jb_i + q_i + 1`. Here we use
the prime-degree zero-run lemma of Section 2.3: the gap between `jb_i` and
`j_i` has length `q_i`, and since `ell_i` lies strictly inside it we have
`q_i >= 1`; the zero-run lemma bounds the length of any non-final zero run by
`d_i - 2`, so `q_i <= d_i - 2`, i.e.

        j_i <= jb_i + d_i - 1 .                                    (15)

Moreover `jb_i < ell_i` and `ell_i <= k - d_1 <= k - d_i`, so `jb_i + 1 !=
k - d_i + 1`; this is the point at which the prime-degree hypothesis is used,
and it is exactly the source's use of its structural lemma. The columns
`jb_i + 1, ..., jb_i + e_i` of `B_i` cover, by (12), the union

    [ jb_i + 2 , jb_i + e_i + d_i - (e_i+1) ] = [ jb_i + 2 , jb_i + d_i - 1 ] ,

which contains `j_i` by (15). Since `jb_i < ell_i < j_i`, we obtain

    union_{k in B_i} [ k+1 , k + d_i - (e_i+1) ]
        contains [ jb_i + 2 , j_i + d_i - (e_i+1) ]
        contains [ ell_i + 1 , ell_{i+1} ] ,

the last inclusion because `ell_i + 1 > jb_i + 1` and
`ell_{i+1} = ell_i + d_i - (e_i+1) <= j_i + d_i - (e_i+1)`. This proves (14).

By (7), the intervals `[ell_i + 1, ell_{i+1}]` for `i = 1,...,t_0` are
contiguous and cover `[1, ell_{t_0+1}]`, and by (8)
`ell_{t_0+1} = sum_{i<=t_0} (d_i - (e_i+1)) < k`. Together with (13), which
covers the tail `[k-d_1+2, k-1]`, the coverage is complete: every
`rho in {1,...,k-1}` lies in one of the intervals (12), so (T2) holds. This
proves (3.1). QED

## 4. Decomposing an arbitrary affine subspace into blocks

Let `A` be an affine subspace of `K^r`, and let `gamma in A`. Gaussian
elimination on an `F_2`-basis of the direction space of `A`, respecting the
order of blocks, gives linear subspaces `V_i` of `K` and linear maps
`sigma_ij : V_j -> K` for `j < i` such that

    A = { ( gamma_i + z_i + sum_{j<i} sigma_ij(z_j) )_{i=1}^r : z_i in V_i } ,  (10)

    dim(A) = sum_i dim(V_i) .                                    (11)

For completeness: filter the direction space of `A` by successively requiring
its first `0, 1, ..., r` blocks to vanish; each successive quotient is split by
choosing an `F_2`-basis; the leading block images form the `V_i`, and the later
coordinates of those basis vectors give the `sigma_ji`. Their direct sum is
the original direction space, which proves (10) and (11). This is ordinary
linear Gaussian elimination over `F_2`; it is not an assumption of direct-sum
additivity of circuit complexity.

Let `Q_i = Q_{V_i}` be the image polynomial of `V_i`, with `Q_i = 0` when
`V_i = {0}`. The compositions `sigma_ij o Q_j` are linear maps on `K`, so by
(1.2) they have linearized polynomial representations `Q_ij`. Put

    R_i(X_1,...,X_r) = gamma_i + Q_i(X_i) + sum_{j<i} Q_ij(X_j) .    (12)

The map `R : K^r -> K^r` defined by `R_i` has image exactly `A`.

## 5. Selecting `t` high-dimension blocks

Suppose `dim(A) = D > k r / t + k t + r (r+1)`. Put

    S = { i : d_i > r + 1 },      d_i = dim(V_i).                  (13)

Since each `d_i <= k`, we have `D <= sum_{i in S} d_i + r (r+1)`, so
`D > k t + r (r+1)` implies `|S| >= t`. Moreover

    sum_{i in S} d_i >= D - r(r+1) > k r / t + k t .              (14)

Choose `I_0 subset S` with `|I_0| = t` and `d_i` as large as possible
(breaking ties by index). Since `|S| <= r`, the sum of the `t` largest `d_i` in
`S` is at least `t/r` times the sum over `S`, so

    sum_{i in I_0} ( d_i - (r+1) ) > k + k t^2 / r - t (r+1) >= k ,  (15)

the last inequality being exactly (D1) rearranged: `k t >= r(r+1)` implies
`k t^2 / r >= t (r+1)`.

For each selected `i`, `e_i = r + 1 - i` satisfies `1 <= e_i` and
`w_i = e_i + 1 = r + 2 - i <= r + 1 < d_i`. Thus `e_i < d_i - 1` and
`d_i > 1`, so the hypotheses of the coefficient assertion (3.1) hold, with

    delta_i = d_i - w_i = d_i - e_i - 1 ,                         (16)

    sum_{i in I_0} delta_i > k >= k - max_{i in I_0} d_i + 1 .

Applying (3.1) with the sources indexed by `I_0` and the parameters
`e_i` yields a monomial

    M = prod_{i in I_0} X_i^{b_i} ,      wt(b_i) = w_i ,          (17)

with a nonzero coefficient in

    Tr( mu prod_{i in I_0} Q_i(X_i)^{2^{w_i} - 1} ) .            (18)

## 6. Triangular mixing does not alter that coefficient

We must compare the coefficient of `M` in the polynomial for `f o R` with the
coefficient in (18).

Let `g` be the polynomial representing `f o R`, i.e. the reduced polynomial of

    sum_{|I|=t} Tr( mu prod_{i in I} R_i^{alpha_i} ) ,             (19)

where `alpha_i = sum_{j=0}^{e_i} p^j` in the source's notation; in the binary
case the exponent attached to block `i` is `2^{e_i} - 1` with the
corresponding Frobenius expansion, which is what `w_i = e_i + 1` records.

**Claim.** For the monomial `M` of (17), the coefficient of `M` in the
`I_0`-summand of (19) equals the coefficient of `M` in (18).

*Proof.* Expand each exponent into its `w_i` Frobenius factors and expand the
trace. Consider the selected indices in decreasing order. No lower-index
factor mentions `X_i`. To achieve binary degree `w_i` in `X_i`, all `w_i`
factors coming from `R_i` must contribute an `X_i` term; hence each must
contribute a term from `Q_i(X_i)`, not a constant and not an earlier
variable, since a constant contributes degree `0` and an earlier variable
contributes degree `0` in `X_i`. Carries cannot increase binary degree, by
Section 1.4, so a term contributing degree `w_i` in `X_i` must use exactly
`w_i` distinct powers and no carries. Having used those factors, repeat with
the next-largest selected index. Descending induction forces the diagonal
contribution in every block. Hence for each trace summand the coefficient of
`M` is the same as in (18), and therefore so is their sum. QED

## 7. No other symmetric summand cancels the monomial

It remains to show that the coefficient of `M` in the whole of (19) is the
nonzero coefficient of (18); other summands indexed by `J != I_0` with
`|J| = t` must not contribute to `M`.

In a trace expansion of the `J`-summand, index `j` supplies `w_j = r + 2 - j`
linearized factors, each mentioning only variables of index `<= j` or a
constant. Suppose a product choice could reduce to `M`. Let `n_ji` count the
factors from `j` choosing the variable `i in I_0`. Then necessarily

    sum_{i in I_0} n_ji <= w_j ,                                  (20)

    sum_{j in J} n_ji >= w_i   for every i in I_0 ,                (21)

    n_ji > 0  implies  i <= j .                                  (22)

The first inequality is immediate from the supply of factors; the second is
binary-degree subadditivity, valid after reduction; the third holds because
factor `j` mentions only variables of index `<= j`.

Now, since `w_i` strictly decreases with `i`, the implication `i <= j` gives
`w_i >= w_j`. Summing (21) over `i in I_0` and comparing with (20):

    t <= sum_{i in I_0, j in J} n_ji / w_i
       <= sum_{i in I_0, j in J} n_ji / w_j
       <= t .                                                     (23)

The first inequality uses the demands (21): `sum_j n_ji >= w_i` gives
`n_ji / w_i >= n_ji / w_j` because `w_j <= w_i`. The last uses the supplies
(20): `sum_i n_ji <= w_j`. All three quantities are `t`, so equality holds
throughout.

If some `n_ji > 0` had `i < j`, then `w_i > w_j` and that particular term
would be strict in the middle inequality, contradicting equality. Hence every
used factor has `i = j`. Each target index `i in I_0` has positive demand by
(21), so `I_0` is contained in `J`; as `|I_0| = |J| = t`, we get `J = I_0`.

Consequently no summand other than the `I_0` one contributes to `M`. The
coefficient of `M` in `f o R` is therefore the nonzero coefficient of (18).

`M` has positive degree, so by (1.1) the reduced polynomial representing
`f o R` is nonconstant. Since `R(K^r) = A`, `f` is nonconstant on `A`. This
proves the disperser lemma (D). QED

## 8. Correspondence with the cited source

The source is Eli Ben-Sasson and Swastik Kopparty, *Affine Dispersers from
Subspace Polynomials*, author preprint of 12 March 2010,
`https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf`, SHA256
`922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf`.

* Their **Theorem 7.3** (printed p. 33) states: for `t < r` integers, `n`
  prime with `n >= r(r+1)/t`, exponents `e_i = r + 1 - i` and
  `alpha_i = sum_{j=0}^{e_i} p^j`, and `mu != 0`, the function

      f(x_1,...,x_r) = Tr( mu sum_{|I|=t} prod_{i in I} x_i^{alpha_i} )

  is such that `|f(A)| > 1` for every `F_p`-affine `A` with
  `dim(A) > n r / t + n t + r (r+1)`. Our (D) is exactly this with `p = 2`,
  `n = k`, and the binary exponent `2^{w_i} - 1` in place of the `p`-ary
  `alpha_i`; the two coincide because over `F_2` the exponent
  `alpha_i = sum_{j=0}^{e_i} 2^j` is `2^{e_i+1} - 1 = 2^{w_i} - 1`.
* Their hypothesis `n >= r(r+1)/t` is exactly our (D1), `k t >= r(r+1)`.
* Their **Lemma 3.7** (existence of an image-subspace polynomial) is
  reproduced in Section 2.1.
* Their **Lemma 3.10** (the main structural lemma on subspace polynomials,
  in particular the exclusion of a run of `d-1` zero coefficients away from
  the final gap) is reproduced, in the binary prime-degree specialisation
  needed here, in Sections 2.2-2.3 (zero-run lemma).
* Their **Lemma 3.13** (the coefficient formula as a sum of permanents) is
  reproduced as equation (6).
* Their **Theorem 6.1** (the coefficient assertion for independent affine
  sources) is reproduced as (3.1), with the construction of the column sets
  `B_i` in (7)-(8) and the zero-column coverage argument in (9).
* Their **Lemma 7.2** (triangular mixing preserves the relevant coefficient)
  is reproduced as the Claim of Section 6.
* Their **Claim 7.4** (uniqueness of the contributing summand) is reproduced in
  Section 7, in the weighted-capacity form (23).
* Their **Theorem 2.2** and its equation (4) give the choice of parameters
  `k ~ m^{3/5}`, `r ~ m^{2/5}`, `t ~ m^{1/5}` and the disperser dimension
  `6 m^{4/5}`; `proof.md` applies (D) with concrete parameters instead of
  those asymptotic ones, so we do not need their dimension bookkeeping.

The only material we take from the source without re-deriving it is the
existence and uniqueness of the reduced polynomial representative (their
Proposition 3.1), the block decomposition by Gaussian elimination (their use of
Lemma 7.1), and the identification of the image polynomial with a linearized
polynomial; all three are elementary and are restated above in Sections 1-2
and 4.

## 9. Scope and what is not claimed

* (D) is a statement about a *fixed* field `K` of size `2^k` with `k` prime,
  and about the specific exponent pattern `2^{w_i}-1` with `w_i = r+2-i`
  decreasing in `i`. It is the source's binary instance, not a new
  disperser.
* (D) requires `pi != 0`; in `proof.md` the map `pi` is the constant
  coefficient functional of `K` relative to the polynomial basis, which is
  nonzero because `pi(1) = 1`.
* The proof of (D) uses the prime-degree hypothesis `k` in exactly one place,
  the zero-run lemma of Section 2.3, which is what excludes the proper-subfield
  obstruction. Without it, the source's corresponding hypothesis (the
  non-containment of the subspaces in an affine shift of a proper subfield)
  would be needed instead.
* No superpolynomial, superlinear or endpoint statement is made here.
* This is a written proof, not a formal proof-assistant certificate.
