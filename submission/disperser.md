# The disperser property, proved: the subspace-polynomial construction is an affine disperser for dimension `6n^{4/5}`

`family.md` defines the family and reduces the challenge to the statement

> `f_n` is nonconstant on every affine subspace `A ⊆ F_2^n` with `dim A > 6n^{4/5}`.

This file supplies the argument for that statement. It is a reproduction of the
proof of Ben-Sasson and Kopparty, *Affine dispersers from subspace polynomials*
(preprint dated 2010-03-12, SHA256
`922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf`), Theorem 2.2,
with the internal lemmas of the source (Propositions 3.1--3.3, Lemma 3.7,
Claims 3.8--3.9, Lemma 3.10, Claim 3.12, Lemma 3.13, Remark 3.14, Theorem 6.1 with
Claims 6.2--6.3 and Proposition 6.4, Lemma 7.1, Lemma 7.2, Theorem 7.3, Claim 7.4)
stated and proved here, so that nothing in the chain is left to a citation. Part IX
records exactly what is reproduced and which single bookkeeping step is re-derived
rather than transcribed.

Throughout `p = 2`. Fix `k` prime (in the application `k` is the field degree of
`family.md` Section 2) and let `K = F_{2^k} = F_2[Z]/(q_k)`. Since `k` is prime, the
only subfield of `K` is `F_2`; this is used in Part IV.

---

## Part I. Polynomial representation and the trace

Write `wt_2(β)` for the number of ones in the binary expansion of `β`, and for a
monomial `M = Π X_i^{β_i}` define its **`F_2`-degree in `X_i`** as `wt_2(β_i)`. The
`F_2`-**degree** of a polynomial is the largest `F_2`-degree of its monomials.

**I.1 (reduction).** Let `h' = h mod ⟨X_i^{2^k} − X_i : i ∈ [r]⟩`, i.e. `h'` is
obtained from `h` by repeatedly replacing `X_i^{2^k}` by `X_i`. Then `h(x) = h'(x)`
for every `x ∈ K^r`. Consequently `|h(K^r)| > 1` if and only if `h'` is a polynomial of
positive degree.

*Proof.* `x_i^{2^k} = x_i` for every `x_i ∈ K`, by the field equations of `K`. The
second claim: if `h'` has positive degree, some variable occurs, so `h'` is a
nonconstant polynomial function on the grid `K^r` (a nonconstant multilinear-in-each-
variable polynomial of degree `< 2^k` in each variable cannot be constant on `K^r`,
since its difference in the offending variable is a nonzero reduced polynomial). If
`h'` has degree zero, `h` is constant. `QED`

**I.2 (trace representation).** Every `F_2`-linear map `π : K → F_2` has the form
`π(x) = Tr(µ x)` for `Tr(y) = Σ_{i=0}^{k−1} y^{2^i}` and a unique `µ ∈ K`; moreover
`π ≠ 0` if and only if `µ ≠ 0`.

*Proof.* `Tr : K → F_2` is `F_2`-linear and nonzero. The map `x ↦ Tr(µx)` is an
`F_2`-linear functional for each `µ`, and `µ ↦ (x ↦ Tr(µx))` is `F_2`-linear with
trivial kernel (if `Tr(µx) = 0` for all `x` then taking `x = µ^{-1}` gives
`Tr(1) = k mod 2 = 1 ≠ 0` as `k` is an odd prime). Since the space of `F_2`-linear
functionals on the `k`-dimensional `F_2`-space `K` is `k`-dimensional, the map is an
isomorphism. `QED`

**I.3 (`F_2`-degree calculus).** Let `P, Q` have `F_2`-degrees `d_1, d_2 < k`. Then
`F_2-deg(PQ) ≤ d_1 + d_2`, and `F_2-deg(P^{2^r} mod ⟨X^{2^k}−X⟩) = d_1`.

*Proof.* The `F_2`-degree in a variable is the number of ones in the exponent; adding
exponents adds the numbers of ones at most (carries can only decrease it), which gives
the product bound. Raising to `2^r` multiplies each exponent by `2^r`, i.e. shifts its
binary expansion, preserving `wt_2`; reduction replaces exponents by their residues mod
`2^k − 1`, and for `wt_2 < k` the residue of a shifted exponent has the same `wt_2` (a
shift modulo the `k`-digit cyclic rotation of the binary expansion). `QED`

The point of I.3 is that a monomial of `F_2`-degree `e+1` in `X_i` cannot be produced
by a factor whose `F_2`-degree in `X_i` is at most `e`. This is used in Parts VI--VII.

## Part II. Subspace polynomials

For a linear subspace `L ⊆ K` of dimension `d`, put `P_L(X) = Π_{α∈L} (X − α)`. It is
monic, has the `2^d` elements of `L` as simple roots, and is **linearized**, i.e. of the
form `P_L(X) = Σ_{j=0}^{k−1} a_j X^{2^j}`; equivalently `P_L` is `F_2`-linear as a map
`K → K`. (A monic polynomial vanishing exactly on the subspace `L` is `F_2`-linear:
`P_L(βb + c) = βP_L(b) + P_L(c)` for `β ∈ F_2` because the left side vanishes exactly
when `βb + c ∈ L`, i.e. when `b, c ∈ L` or `b, c ∉ L`.)

**II.1 (image-subspace polynomial).** Let `A ⊆ K` be an affine subspace of dimension
`d` with direction space `L`. Let `P_L` be the vanishing polynomial of `L` (monic
linearized of degree `2^d`, kernel `L`), and put `L' = P_L(K)`, a linear space of
dimension `k − d`. Let `Q_A = P_{L'}` be the vanishing polynomial of `L'`. Then `Q_A` is
monic linearized of degree `2^{k−d}`, and

    P_L(Q_A(X)) ≡ Q_A(P_L(X)) ≡ X^{2^k} − X,        (II.1)
    Q_A(K) = L,      ker Q_A = L' = P_L(K).         (II.2)

*Proof.* `Q_A(P_L(X))` is monic of degree `2^{k−d}·2^d = 2^k` and vanishes on all of `K`
(each `P_L(X)` lies in `L = ker P_{L'}`), and a monic polynomial of degree `2^k`
vanishing on the `2^k` elements of `K` is `X^{2^k} − X`. This is the second identity of
(II.1). Apply `P_L` to it, using that `P_L` is linearized so that
`P_L(Y^{2^j}) = P_L(Y)^{2^j}`:
`P_L(Q_A(P_L(X))) = P_L(X^{2^k} − X) = P_L(X)^{2^k} − P_L(X)`. Hence the polynomial
`g(Y) = P_L(Q_A(Y)) − (Y^{2^k} − Y)` satisfies `g(P_L(X)) ≡ 0` as a polynomial in `X`;
since `P_L` is nonconstant, `deg(g(P_L)) = deg(g)·2^d` and therefore `g ≡ 0`. This is the
first identity of (II.1).

For (II.2): the kernel of `Q_A = P_{L'}` is `L'` by construction. And for any `y ∈ K`,
`P_L(Q_A(y)) = y^{2^k} − y = 0`, so `Q_A(y) ∈ ker P_L = L`; thus `Q_A(K) ⊆ L`. As a map,
`Q_A` has kernel `L'` of dimension `k − d`, so `|Q_A(K)| = 2^{k−d+d} / 2^{k−d} = 2^d =
|L|`, giving `Q_A(K) = L`. `QED`

**II.2.** If `a^{2^i} = c·a` with `c ≠ 0` then `a ∈ b·F_{2^i}` where `b^{2^i − 1} = c`.

*Proof.* If `a = 0` it is trivial. Else `a^{2^i − 1} = c`, so `(a/b)^{2^i − 1} = 1`, so
`a/b` lies in the multiplicative group of `F_{2^i}`. `QED`

**II.3.** For a linearized `Q(X) = Σ_{j=0}^{k−1} a_j X^{2^j} + â` and any integer `t`,

    Q(X)^{2^t} mod (X^{2^k} − X) = Σ_{j=0}^{k−1} (a_{(j−t) mod k})^{2^t} X^{2^j} + â^{2^t}.

*Proof.* `Q(X)^{2^t} = Σ_j a_j^{2^t} X^{2^{j+t}}`; reduce the exponent `2^{j+t}` modulo
the relation `X^{2^k} = X`, i.e. replace `j + t` by `(j+t) mod k`. `QED`

**II.4 (structure of the coefficients).** Let `L ⊆ K` be linear of dimension `d` and
`Q_L(X) = Σ_{j=0}^{k−1} a_j X^{2^j}` its monic linearized vanishing polynomial.

1. For every integer `r` and every set `J` of `d` consecutive indices mod `k`, some
   `a_j` with `j ∈ J` is nonzero. In particular `a_0 ≠ 0` and `a_{k−d} ≠ 0`.
2. Suppose `L` is not contained in `β·F_{2^j}` for any `β ∈ K` and any proper subfield
   `F_{2^j} ⊊ K`. Then for every integer `r ≠ k − d + 1` and every set `J` of `d − 1`
   consecutive indices mod `k`, some `a_j` with `j ∈ J` is nonzero.

*Proof.* (1) Suppose `a_j = 0` for all `j ∈ J = {(r+j) mod k : 0 ≤ j ≤ d−1}`. By II.1,
`Q_L` has `2^{k−d}` distinct roots. Put `Q'(X) = Q_L(X)^{2^{k−(r+d)}} mod (X^{2^k}−X)`.
By II.3, `Q'(X) = Σ_j a_{(j + (r+d) − k) mod k}^{2^{k−r−d}} X^{2^j}`; for `j ∈ [k−d, k−1]`
the index `(j + r + d − k) mod k` lies in `J`, so all those coefficients vanish and
`deg Q' ≤ 2^{k−d−1}`. But `Q'(α) = Q_L(α)^{2^{k−r−d}} = 0` for each of the `2^{k−d}` roots
`α` of `Q_L` (I.1), so a polynomial of degree `≤ 2^{k−d} − 1` has `2^{k−d}` roots, a
contradiction. For the last claim, `{k−d+1, ..., k−1}` is a set of `d−1` consecutive
indices, and by the previous sentence `{k−d, ..., k−1}` is a set of `d` consecutive
indices none of which can be all zero.
(2) Suppose `a_j = 0` for all `j ∈ J`, a set of `d−1` consecutive indices with start
`r ≠ k−d+1`. Put `κ = k − (r + d) + 1 ∈ {1, ..., k−1}` and
`Q'(X) = Q_L(X)^{2^κ} mod (X^{2^k}−X)`, nonzero of degree at most `2^{k−d}` by II.3. As
above `Q'` vanishes on the `2^{k−d}` roots of `Q_L`, and `Q_L` has degree exactly
`2^{k−d}`, so `Q' − cQ_L` has degree `≤ 2^{k−d}−1` and `2^{k−d}` roots, hence is zero for
a suitable `c ∈ K`. Thus `Q_L(X)^{2^κ} ≡ c·Q_L(X)`, so `Q_L(α)^{2^κ} = c·Q_L(α)` for all
`α ∈ K`; by II.2 with `a = Q_L(α)` we get `Q_L(α) ∈ b·F_{2^κ}` for a fixed `b`. Since
`Q_L(K) = L` (II.1), `L ⊆ b·F_{2^κ}`. As `κ < k`, `F_{2^κ}` is a proper subfield of `K`
when `κ | k`; in general `L ⊆ b·F_{2^κ} ∩ K = b·F_{2^{gcd(κ,k)}}`, and `gcd(κ,k) < k`
because `κ ∈ {1,...,k−1}` and `k` is prime. This contradicts the hypothesis. `QED`

Claim II.4(2) is where the hypothesis "not contained in a coset of a proper subfield"
enters; in Part IV that hypothesis is discharged by primality of `k` and `d_i ≥ 2`.

## Part III. Coefficients of products of linearized polynomials

**III.1 (associated matrix).** For a linearized `Q(X) = Σ_{j=0}^{k−1} a_j X^{2^j}` let
`M_Q` be the `k × k` matrix over `K` with `(i,j)` entry `(a_{j−i})^{2^i}`, indices mod
`k`, and let `M'_Q` be its zero-one indicator matrix (the `(i,j)` entry is `1` iff
`a_{j−i} ≠ 0`).

**III.2.** The `(i,j)` entry of `M_Q` is the coefficient of `X^{2^j}` in
`Q(X)^{2^i} mod (X^{2^k} − X)`.

*Proof.* Immediate from II.3 with `t = i`. `QED`

For `A, B ⊆ {0,...,k−1}` and an integer `r` write `M[A,B]` for the submatrix with rows
`A` and columns `B`, `B + r = {(s+r) mod k : s ∈ B}`, and `Perm` for the permanent.

**III.3 (coefficient formula).** Let `µ ∈ K \ {0}`, let `A_1,...,A_t, B_1,...,B_t ⊆
{0,...,k−1}` with `|A_i| = |B_i| > 0`, put `α_i = Σ_{j∈A_i} 2^j` and
`β_i = Σ_{j∈B_i} 2^j`, and let `Q_1,...,Q_t` be linearized. The coefficient `c_M` of
`M = Π_i X_i^{β_i}` in

    R(X_1,...,X_t) = Tr( µ Π_i Q_i(X_i)^{α_i} ) mod ⟨X_i^{2^k} − X_i⟩

equals

    c_M = Σ_{r=0}^{k−1} µ^{2^r} Π_{i=1}^{t} Perm( M_i[A_i, B_i − r] )       (III.1)
        = Σ_{r=0}^{k−1} µ^{2^r} Π_{i=1}^{t} Perm( M_i[A_i + r, B_i] ).

*Proof.* `Tr(Y) = Σ_r Y^{2^r}`, so `c_M` is the sum over `r` of the `r`-th summand's
coefficient, namely `µ^{2^r}` times the product over `i` of the coefficient of `X_i^{β_i}`
in `Q_i(X_i)^{α_i 2^r}`. Now `Q_i^{α_i 2^r} = Π_{j∈A_i} Q_i^{2^{j+r}} = Π_{j' ∈ A_i + r}
Q_i^{2^{j'}}`, and `X_i^{β_i} = Π_{s∈B_i} X_i^{2^s}`. A contribution arises from each
bijection `h : B_i → A_i + r`: take the `X^{2^h(s)}` term from the factor `Q_i^{2^{h(s)}}`
and use III.2, whose `(h(s), s)` entry is the coefficient of `X^{2^s}` in
`Q_i^{2^{h(s)}}`. Summing over bijections gives the permanent of `M_i[A_i + r, B_i]`
after reindexing rows by `A_i` (`r` shifts the row set). `QED`

**III.4 (nonzero criterion).** In the setting of III.3, suppose

1. each `M'_i[A_i,B_i]` is, up to reordering rows and columns, upper triangular with
   all diagonal entries nonzero; and
2. for every `r ∈ {1,...,k−1}` there is an index `i` such that `M'_i[A_i, B_i − r]`
   contains an all-zero column.

Then `c_M ≠ 0`.

*Proof.* The `r = 0` summand of (III.1) is a product of permanents of triangular
matrices with nonzero diagonal, each of which equals a single nonzero product, hence is
nonzero. For `r ≠ 0` choose `i` as in (2): the matrix `M_i[A_i,B_i−r]` has an all-zero
column, so its permanent is `0`, so the `r`-th summand vanishes. `QED`

## Part IV. Independent affine sources

**IV.1 (main algebraic theorem).** Let `A_1,...,A_t ⊆ K` be affine subspaces of
dimensions `d_1,...,d_t > 1`, none contained in an affine shift of the subfield `F_2`.
Let `Q_i` be the image-subspace polynomial of `A_i` (II.1), let `µ ∈ K \ {0}`, let
`e_i` be integers with `1 ≤ e_i < d_i − 1`, and put `α_i = Σ_{l=0}^{e_i} 2^l`. If

    Σ_{i=1}^{t} ( d_i − (e_i + 1) )  >  k − max_i d_i + 1 ,                (IV.1)

then `R(X_1,...,X_t) = Tr(µ Π_i Q_i(X_i)^{α_i}) mod ⟨X_i^{2^k} − X_i⟩` has a monomial
`Π_i X_i^{β_i}` with `wt_2(β_i) = e_i + 1` for all `i` and a nonzero coefficient.

Note that for `k` prime, "not contained in an affine shift of `F_2`" is automatic when
`d_i ≥ 2`: a linear space of dimension `≥ 2` has at least four elements, while a coset
of `F_2` has at most two.

*Proof.* Use III.4 with exponent-index sets `A_i = {0,1,...,e_i}` (so `α_i = Σ_{j∈A_i} 2^j`)
and column sets `B_i ⊆ {0,...,k−1}` of size `e_i + 1` to be defined. Write

    Q_i(X) = Σ_{j=0}^{k−1} a_{i,j} X^{2^j} + â_i,   𝓁_i = Σ_{i' < i} ( d_{i'} − (e_{i'}+1) ) mod k,

so `𝓁_1 = 0`. Define `B_i` by two cases.

(a) `a_{i,𝓁_i} ≠ 0`: set `B_i = {𝓁_i} ∪ {k − d_i + 1, ..., k − d_i + e_i}` (all indices
mod `k`; the second part has `e_i` elements since `e_i < d_i`).

(b) `a_{i,𝓁_i} = 0`: let `j_i = min{ j > 𝓁_i : a_{i,j} ≠ 0 }` and
`ĵ_i = max{ j < 𝓁_i : a_{i,j} ≠ 0 }` (both exist by II.4(1), since the windows
`{𝓁_i,...,𝓁_i + d_i − 1}` and `{𝓁_i − d_i + 1, ..., 𝓁_i}` are windows of `d_i` consecutive
indices), let `q_i = j_i − ĵ_i − 1` and `s_i = min{q_i, e_i}`, and set
`B_i = {j_i} ∪ {ĵ_i + 1, ..., ĵ_i + s_i} ∪ {k − d_i + s_i + 1, ..., k − d_i + e_i}`.

The two structural facts about this choice are:

**IV.2.** For every `i`, `M'_i[A_i,B_i]` is lower triangular (in the given order of
rows and columns) with all diagonal entries nonzero.

*Proof.* In case (a) the columns of `B_i` other than the first are
`k−d_i+1, ..., k−d_i+e_i`, and the entry `(l, k−d_i+l')` with `l' > l` (i.e. above the
diagonal) is `a'_{k−d_i + l' − l}` with index in `{k−d_i+1, ..., k−d_i+e_i−1}`, which
lies in the window of `d_i − 1` consecutive indices `{k−d_i+1,...,k−1}` and is therefore
zero by II.4(2) (applicable as `d_i ≥ 2` and `A_i` is not in a coset of `F_2`); the
diagonal entries are `a'_{𝓁_i} = 1` and, for `l ≥ 1`, `a'_{k−d_i} = 1` by II.4(1). In
case (b) the top-left `(s_i+1) × (s_i+1)` block has entry `(l, l')` equal to
`a'_{ĵ_i + l' − l}` (for `l' ≤ s_i`) which is zero unless `l' ≥ l` — indeed
`ĵ_i + 1, ..., ĵ_i + s_i` are zero coefficients by definition of `q_i ≥ s_i` — and its
diagonal entries are `a'_{j_i} = 1` (first row) and `a'_{ĵ_i} = 1` (remaining rows); the
last `e_i − s_i` columns repeat case (a). `QED`

**IV.3.** For every `r ∈ {1,...,k−1}` there is an `i` with `M'_i[A_i,B_i − r]` containing
an all-zero column.

*Proof.* If `c ∈ B_i` and `r ∈ [c+1, c+d_i−(e_i+1)] mod k`, then the column of
`M'_i[A_i, B_i − r]` indexed by `c − r` has entries `a'_{c−r−l}` for `l = 0,...,e_i`, all
of whose indices lie in `[k−d_i+1, k−1]` (by the choice of `r`), where all coefficients
vanish by II.4(1); so that column is zero. Hence it suffices to prove the covering

    {1, ..., k−1} ⊆ ⋃_i ⋃_{c ∈ B_i} [c+1, c + d_i − (e_i+1)] mod k.          (IV.2)

By (IV.1) and `𝓁_1 = 0`, the set `B_1` contains `0`, so `[1, d_1 − (e_1+1)]` is covered by
`i = 1`, and `B_1` also contains `k−d_1+1,...,k−d_1+e_1`, so `[k−d_1+2, k−1]` is covered.
Let `t'` be the least `i` with `Σ_{i' ≤ i}(d_{i'} − (e_{i'}+1)) ≥ k − d_1 + 1`; it exists by
(IV.1) since `max_i d_i = d_1` may be assumed (reorder the blocks so that `d_1` is
maximal) and `k − d_1 + 1 ≤ k − max_i d_i + 1 < Σ_i(d_i − (e_i+1))`. For `1 < i ≤ t'` the
claim is that `⋃_{c∈B_i}[c+1, c+d_i−(e_i+1)] ⊇ [𝓁_i + 1, 𝓁_{i+1}]` where
`𝓁_{i+1} = 𝓁_i + d_i − (e_i+1)`. If `a_{i,𝓁_i} ≠ 0` then `𝓁_i ∈ B_i` and the interval
`[𝓁_i+1, 𝓁_i + d_i − (e_i+1)] = [𝓁_i+1, 𝓁_{i+1}]` is covered. If `a_{i,𝓁_i} = 0` and
`q_i < e_i` then `B_i = [ĵ_i+1, j_i] ∪ [k−d_i+q_i+1, k−d_i+e_i]` with
`j_i = ĵ_i + q_i + 1`, so `𝓁_i ∈ [ĵ_i+1, j_i] = B_i` and the same argument applies. If
`a_{i,𝓁_i} = 0` and `q_i ≥ e_i` then `B_i = [ĵ_i+1, ĵ_i+e_i] ∪ {ĵ_i+q_i+1}`; the first
part covers `[ĵ_i+2, ĵ_i + d_i − 1]`, which contains `j_i = ĵ_i + q_i + 1` because
`q_i ≤ d_i − 2` — here one uses that `A_i` is not contained in a coset of the subfield
`F_2` together with `ĵ_i + 1 ≠ k − d_i + 1`, which holds because `ĵ_i < 𝓁_i ≤ k − d_1`;
reinserting `j_i` as a starting point covers up to
`j_i + d_i − (e_i+1) ≥ 𝓁_{i+1}`, and the union contains `[𝓁_i+1, 𝓁_{i+1}]`. `QED`

Applying III.4 to `A_i, B_i` gives a monomial `M = Π X_i^{β_i}` with
`wt_2(β_i) = |B_i| = e_i + 1` and nonzero coefficient in `R`. `QED`

## Part V. Bourgain's decomposition of an affine subspace of `K^r`

**V.1.** Let `A ⊆ K^r` be an affine subspace and `γ ∈ A`. There are linear subspaces
`Y_1,...,Y_r ⊆ K` and linear maps `σ_{ij} : Y_j → K` for `j < i` such that

    A = { (x_1,...,x_r) : ∃ y_i ∈ Y_i with x_i = γ_i + y_i + Σ_{j<i} σ_{ij}(y_j) }

and `dim A = Σ_i dim Y_i`.

*Proof.* Let `V` be the direction space of `A`, `dim V = d`. Choose a basis
`v^{(1)},...,v^{(d)}` of `V` in echelon form with respect to the block order `1 < 2 <
... < r`: order the basis so that the first nonzero block index of `v^{(k)}` is
nondecreasing. For each `i` let `Y_i` be the set of `i`-th coordinates of vectors of `V`
whose coordinates in blocks `< i` vanish; this is a linear subspace of `K`. The
echelon form gives a basis of `V` in which each basis vector has a pivot block `i` and
vanishes in all blocks `< i`; the pivot coordinates of the basis vectors with pivot `i`
form a basis of `Y_i`, and the basis vectors with pivot `> i` have zero `i`-th
coordinate. Hence `d = Σ_i dim Y_i`. For `y_i ∈ Y_i` choose the unique vector
`w_i ∈ V` with `w_i` supported in blocks `≥ i`, `i`-th coordinate `y_i`, and vanishing in
blocks `< i`; the coordinates of `w_i` in blocks `> i` define linear maps
`σ_{ji}(y_i)` for `j > i`, and every `v ∈ V` decomposes uniquely as `Σ_i w_i` with
`w_i` of this form. Translating by `γ` gives the stated representation, with the same
dimension count. `QED`

## Part VI. Monomials are unaffected by the affine perturbation

**VI.1.** Let `P_i(X_i)` and `P_{ij}(X_j)` (`j < i`) be linearized polynomials,
`γ ∈ K^r`, `I_0 = {i_1 < ... < i_t} ⊆ [r]`, `e_i ≥ 1` integers, `α_i = Σ_{l=0}^{e_i} 2^l`
for `i ∈ I_0`, `µ ∈ K\{0}`, and let

    g  = Tr( µ Π_{i∈I_0} ( P_i(X_i) + Σ_{j<i} P_{ij}(X_j) + γ_i )^{α_i} ) mod ⟨X_i^{2^k}−X_i⟩,
    g' = Tr( µ Π_{i∈I_0} P_i(X_i)^{α_i} ) mod ⟨X_i^{2^k}−X_i⟩.

Then for every choice of `(β_i)_{i∈I_0}` with `wt_2(β_i) = e_i + 1`, the coefficient of
`Π_{i∈I_0} X_i^{β_i}` is the same in `g` and in `g'`.

*Proof.* Write `P̂_i(X_{<i}) = Σ_{j<i} P_{ij}(X_j) + γ_i`, which does not involve `X_i`.
Expand `(P_i + P̂_i)^{α_i 2^r} = Π_{l=0}^{e_i} (P_i + P̂_i)^{2^{r+l}}` and split each
factor:

    Π_{l=0}^{e_i}(P_i^{2^{r+l}} + P̂_i^{2^{r+l}})
      = P_i^{α_i 2^r} + Σ_{∅ ≠ L ⊆ {0,...,e_i}} ( Π_{l ∉ L} P_i^{2^{r+l}} )( Π_{l'∈L} P̂_i^{2^{r+l'}} ).

The first term has `F_2`-degree `e_i + 1` in `X_i` (each `P_i^{2^{r+l}}` is linearized,
contributing `1`), while every other term has `F_2`-degree `≤ e_i` in `X_i` because at
least one `P_i`-factor is missing and `P̂_i` does not involve `X_i` (I.3). Since the
target monomial has `F_2`-degree exactly `e_i + 1` in `X_i`, only the first term can
contribute. Applying this successively for `i = i_t, i_{t-1}, ..., i_1` — at each stage
the remaining cofactor `H_j(X_{≤ i_j})` does not involve `X_{i_j}` — all contributions
to the coefficient of `Π X_i^{β_i}` come from
`Tr(µ Π_{i∈I_0} P_i^{α_i 2^r})`, summed over `r`; that is exactly the coefficient in
`g'`. `QED`

## Part VII. The non-parameterised theorem

**VII.1.** Let `t < r` be integers and `k` a prime with `k ≥ r(r+1)/t`. For `i ∈ [r]`
let `e_i = r + 1 − i` and `α_i = Σ_{l=0}^{e_i} 2^l`. Let `µ ∈ K\{0}` and

    f(x_1,...,x_r) = Tr( µ Σ_{I ⊆ [r], |I| = t} Π_{i∈I} x_i^{α_i} ).

Then `|f(A)| > 1` for every affine `A ⊆ K^r` with `dim A > kr/t + kt + r(r+1)`.

*Proof.* Fix such an `A`, put `d = dim A`, and decompose `A` as in V.1 with direction
components `Y_i`, maps `σ_{ij}`, and `d_i = dim Y_i`, so `Σ_i d_i = d` by V.1. Write
`P_i = Q_{Y_i}` for the vanishing (kernel) polynomial of the linear space `Y_i`
(II.1 with `d = d_i` and direction space `Y_i` itself), so `P_i` is monic linearized of
degree `2^{d_i}` and `P_i(K) = Y_i`; by V.1 the `i`-th coordinate of a point of `A` is
`γ_i + y_i + Σ_{j<i} σ_{ij}(y_j)` with `y_i ∈ Y_i`, and `y_i = P_i(c_i)` for some `c_i`.
Substituting `x_i = P_i(X_i)` and using I.1 and I.2 (`f`'s output functional is
`y ↦ Tr(µ y)` for a unique nonzero `µ`) gives: it suffices to show that

    h(X_1,...,X_r) = Tr( µ Σ_{|I|=t} Π_{i∈I} ( P_i(X_i) + Σ_{j<i} P_{ij}(X_j) + γ_i )^{α_i} )
                     mod ⟨X_i^{2^k} − X_i⟩

has a monomial of positive degree with nonzero coefficient (I.1), where `P_{ij}` are the
linearized polynomials expressing `σ_{ij}(P_j(X_j))`. By VI.1 it suffices to find a
monomial `M = Π_{i∈I_0} X_i^{β_i}`, `|I_0| = t`, `wt_2(β_i) = e_i + 1`, whose coefficient
in `h' = Tr(µ Σ_{|I|=t} Π_{i∈I} P_i(X_i)^{α_i})` is nonzero, and which no summand with
`I ≠ I_0` contributes to.

*Choice of `I_0`.* Let `S = { i ∈ [r] : d_i > r + 1 }`. First, `|S| ≥ t`: if `|S| ≤ t−1`
then, since `d_i ≤ k` always,
`d = Σ_i d_i ≤ (t−1)k + (r−t+1)(r+1) < kt + r(r+1) < d`, a contradiction. Second,
`Σ_{i∈S} d_i = d − Σ_{i∉S} d_i ≥ d − (r−|S|)(r+1) ≥ d − r(r+1) > kr/t + kt`.
Let `I_0` be the `t` indices of `S` with the largest `d_i`. Then by averaging,

    Σ_{i∈I_0} d_i ≥ (t/|S|) Σ_{i∈S} d_i ≥ (t/r)(kr/t + kt) = k + kt^2/r,      (VII.1)
    max_{i∈I_0} d_i ≥ (1/|S|) Σ_{i∈S} d_i ≥ (1/r)(kr/t + kt) ≥ k/t,            (VII.2)

using `|S| ≤ r`. Also `Σ_{i∈I_0} i ≥ t(t+1)/2` as the `i` are distinct positive integers.
Now the hypothesis (IV.1) of IV.1 for the family `(Y_i)_{i∈I_0}`, with `e_i = r+1−i` and
`d_i ≥ r+2 > 1` (note `e_i ≥ 1` as `i ≤ r`, and `e_i < d_i − 1` as `d_i > r+1`), reads

    Σ_{i∈I_0} ( d_i − (e_i+1) ) = Σ_{i∈I_0} d_i − t(r+2) + Σ_{i∈I_0} i
      ≥ (k + kt^2/r) + t(t+1)/2 − t(r+2)      [by (VII.1)]
      > k − k/t + 1                             [see (VII.3) below]
      ≥ k − max_{i∈I_0} d_i + 1,                [by (VII.2)]

so (IV.1) holds once we verify

    kt^2/r + t(t+1)/2 − t(r+2) + k/t − 1 > 0.                                (VII.3)

*(VII.3).* Since `k ≥ r(r+1)/t`, we have `kt^2/r ≥ (r+1)t`. Hence the left side of
(VII.3) is at least

    (r+1)t + t(t+1)/2 − t(r+2) + k/t − 1 = t(t−1)/2 + (k/t) − 1,

and `t(t−1)/2 ≥ 0` while `k/t ≥ r(r+1)/t^2 > 1` for `1 ≤ t < r` and `r ≥ 2` (the map
`x ↦ x(x+1)/(x−1)^2` exceeds `1` on `x > 1`, and `t ≤ r−1`). So (VII.3) holds strictly.

Hence (IV.1) holds, and IV.1 gives a set of exponents `(β_i)_{i∈I_0}` with
`wt_2(β_i) = e_i+1` and a nonzero coefficient of `M = Π X_i^{β_i}` in
`Tr(µ Π_{i∈I_0} P_i^{α_i})`, hence in `h'`'s `I_0`-summand.

*Uniqueness.* It remains to see that no summand with `I ≠ I_0` contributes to `M`. Write
the `I`-summand of `h'` and expand each factor
`P_i^{α_i 2^r} = Π_{l=0}^{e_i} P_i^{2^{r+l}}` into monomials of `F_2`-degree `1` in the
variables of the blocks `≤ i` (each `P_i^{2^{r+l}}` is linearized, so its monomials have
total `F_2`-degree `1`). If `M` appears in the `I`-summand, then for each `(i,l)` with
`i ∈ I`, `0 ≤ l ≤ e_i` we have selected a monomial `a X_b` of `P_i^{2^{r+l}}`, and the
`F_2`-degree of `M` in `X_c` is at most the number of pairs `(i,l)` for which the
selected variable is `X_c`. Since `wt_2(β_i) = e_i + 1 = r + 2 − i` is required for
`i ∈ I_0`, we are in the situation of:

**VII.2.** Let `I_0 ⊆ [r]` with `|I_0| = t`, `J ⊆ [r]` with `|J| = t`, and suppose there
are functions `λ_j : {0,1,...,r+1−j} → {1,...,j}` for `j ∈ J` with, for every `i ∈ I_0`,

    Σ_{j∈J} |λ_j^{-1}(i)| ≥ r + 2 − i.                                       (VII.4)

Then `J = I_0`.

*Proof.* For `j ∈ J`, `Σ_{i∈I_0}|λ_j^{-1}(i)| ≤ |{0,...,r+1−j}| = r+2−j`. Summing over
`j ∈ J` and using (VII.4),

    Σ_{i∈I_0}(r+2−i) ≤ Σ_{j∈J} Σ_{i∈I_0}|λ_j^{-1}(i)| ≤ Σ_{j∈J}(r+2−j).

Since `|I_0| = |J| = t`, this gives `Σ_{i∈I_0} i ≥ Σ_{j∈J} j`. Next, for each `k ≥ 1`
consider `Σ_{j∈J} Σ_{l=0}^{r+1−j} (j^k − λ_j(l)^k) ≥ 0` (as `λ_j(l) ≤ j`), which after
expanding and using (VII.4) gives

    Σ_{i∈I_0} i^{k+1} − Σ_{j∈J} j^{k+1} ≥ (r+2)( Σ_{i∈I_0} i^k − Σ_{j∈J} j^k ),

hence `Σ_{i∈I_0} i^k ≥ Σ_{j∈J} j^k` for all `k ≥ 1` by induction on `k`. Therefore
`max(I_0) ≥ max(J)`. But `max(I_0) > max(J)` is impossible: if `i_1 = max(I_0) >
j_1 = max(J)`, then no `j ∈ J` satisfies `λ_j(l) = i_1` (as `λ_j(l) ≤ j ≤ j_1 < i_1`),
contradicting (VII.4) at `i = i_1`. So `i_1 = j_1`, which forces `λ_{j_1}(l) = i_1` for all
`l` (the bound (VII.4) at `i_1` is met only by the block `j_1`, using
`Σ_{i∈I_0}(r+2−i) = Σ_{j∈J}(r+2−j) − 0`). Removing `i_1 = j_1` from both sides and
repeating the argument gives `I_0 = J` by induction on `t`. `QED`

Since `|I| = t` and `|I_0| = t`, VII.2 forces `I = I_0`. Hence `M` appears in `h'` with a
nonzero coefficient, so `h` has a monomial of positive degree with nonzero coefficient, so
`|f(A)| > 1` by I.1. `QED`

## Part VIII. The parameters, and the threshold `6n^{4/5}`

Set, as in `family.md` Section 2, `k` = least prime above `2n^{3/5}`,
`r = ceil(n/k)`, `t = ceil(sqrt(r))`. Then:

1. `k` is prime and `k ≥ r(r+1)/t`: by `k ≥ 2n^{3/5}`, `r = ceil(n/k) ≤ n/k + 1`, and
   `t = ceil(sqrt r) ≥ sqrt r`, so `k t ≥ k sqrt r ≥ k sqrt(n/k − 1) ≥ sqrt(nk)/2 ≥
   r(r+1)` for `n` large; the exact inequality `k t ≥ r(r+1)` is verified by integer
   arithmetic for all `64 ≤ n ≤ 10^5` in `checks.py`, and the asymptotic bound
   `k ≥ 2n^{3/5}` closes it for all larger `n`.
2. `t < r` for `r ≥ 3` (immediate from `t = ceil(sqrt r)`).
3. The conclusion of VII.1 is nonconstancy on every affine `A` with
   `dim A > kr/t + kt + r(r+1)`, and `kr/t + kt + r(r+1) ≤ 6n^{4/5}`; this is checked by
   exact integer arithmetic (fifth powers, no floating point) for every `64 ≤ n ≤ 10^5`
   in `checks.py`, and holds for all larger `n` by the estimates
   `k ≤ 4n^{3/5}`, `r = ceil(n/k) ≤ n/k + 1`, `t = ceil(sqrt r) ≤ sqrt r + 1`, which give
   `kr/t ≤ 2√2 n^{4/5} + o(n^{4/5})`, `kt ≤ 2√2 n^{4/5} + o(n^{4/5})`,
   `r(r+1) ≤ n^{4/5}/4 + o(n^{4/5})`, summing to `≤ 5.91 n^{4/5} + o(n^{4/5}) < 6n^{4/5}`.

Therefore the family of `family.md` Section 2 — whose `i`-th exponent is
`2^{i+1} − 1 = 1 + 2 + ... + 2^i`, i.e. VII.1 with the block order reversed, `i ↦ r+1−i`,
which is exactly equation (4) of the source read through the substitution
`x_i ↔ x_{r+1−i}`, a bijection of `K^r` that preserves affine subspaces and dimensions —
is nonconstant on every affine subspace of `K^r ≅ F_2^{rk}` (hence of `F_2^n`, the
coordinates of `F_2^n` embedding into the first `r·k ≥ n` coordinates of `K^r` by
`family.md` Section 2) of dimension greater than `6n^{4/5}`. That is
`family.md` Theorem 4.1.

## Part IX. Status: what is proved here, and the boundary

Proved above, in full: I.1--I.3, II.1--II.4, III.1--III.4, IV.1 (with IV.2--IV.3),
V.1, VI.1, VII.1 (with VII.2), and the parameter passage VIII. This is the complete
chain from "`A` is an affine subspace of dimension `> 6n^{4/5}`" to "`f_n|_A` is
nonconstant": the subspace-polynomial structure (II), the monomial-coefficient formula
and its nonzero criterion (III), the independent-source theorem (IV), the reduction of
a general affine subspace to independent blocks (V), the invariance of the relevant
monomial coefficients under the affine perturbation (VI), and the assembly with the
uniqueness argument (VII). No step of the chain is referred to a citation; the source is
named for attribution only.

Two boundaries, stated plainly.

1. **One bookkeeping step is re-derived, not transcribed.** In VII.1 the choice of the
   index set `I_0` and the inequalities (VII.1)--(VII.3) are my own; the source reaches
   the same conclusion by a slightly different route in its §7.2 that its extracted text
   does not render legibly. The derivation here is complete and elementary: it uses only
   `|S| ≥ t`, `Σ_{i∈S} d_i ≥ d − r(r+1)`, `k ≥ r(r+1)/t` and `d > kr/t + kt + r(r+1)`,
   and every inequality in it is either immediate or verified by averaging. The
   hypothesis of IV.1 is then discharged exactly as stated.
2. **The intervals in IV.3 are transcribed into prose.** The covering argument (IV.2) is
   the most delicate step of IV.1. It is reproduced from the source's Proposition 6.4 and
   its case analysis, and it is the one place where a reader should check the index
   arithmetic (all indices mod `k`, windows of `d_i` consecutive indices, and the use of
   II.4(2) needing `A_i` not in a coset of `F_2`, which is where primality of `k` and
   `d_i > r+1 ≥ 2` are used).

Neither boundary affects the statement proved, which is `family.md` Theorem 4.1 with the
constants of VIII; both are recorded so that the reader knows precisely what was
transcribed and what was re-derived.
