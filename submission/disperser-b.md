# Disperser, part B: independent sources and sublinear dimension (Ben-Sasson–Kopparty, §6–7)

This file reproduces, with proofs, §6 and §7 of Ben-Sasson–Kopparty,
*Affine Dispersers from Subspace Polynomials*: the independent-sources
disperser (Theorem 6.1) and the proof of Theorem 2.2 from it. It uses
`disperser-a.md` (their §3) freely. Statements keep BK’s numbering.

## 6. Disperser for independent affine sources

Throughout, p = 2 and F = F_{2^n}; "affine" and "dimension" are over F_2.
From §3: for an affine space A ⊆ F with image-subspace polynomial Q_A, the
coefficient support of Q_A is governed by Lemma 3.10; M_i and M_i′ are the
associated matrix of Definition 3.11 and its zero-one indicator;
coefficients of trace polynomials are computed by Lemma 3.13 and isolated
by the strategy of Remark 3.14.

**Theorem 6.1 (BK Theorem 6.1).** Let A_1, …, A_t ⊆ F be affine spaces of
dimensions d_1, …, d_t > 1, none of them contained in an affine shift of a
proper subfield of F. Let Q_i(X_i) ∈ F[X_i] be the image-subspace polynomial
of A_i, let µ ∈ F \ {0}, and let e_1, …, e_t satisfy 1 ≤ e_i < d_i − 1. Put

    α_i = Σ_{j=0}^{e_i} 2^j .

If

    Σ_{i=1}^t ( d_i − (e_i + 1) ) > n − max_{1≤i≤t} d_i + 1 ,

then the polynomial

    R = Tr( µ · ∏_{i=1}^t (Q_i(X_i))^{α_i} )  mod ⟨ X_i^{2^n} − X_i ⟩_{i∈[t]}

has a monomial ∏_{i=1}^t X_i^{β_i} with wt_2(β_i) = e_i + 1 whose coefficient
is nonzero. In particular |R(F^t)| > 1, by Proposition 3.1.

(Theorem 2.3 is the case e_1 = … = e_t = 1, µ = 1, via Propositions 3.1
and 3.3: then α_i = 1 + 2 and the conclusion is that
π(∏_{i=1}^t x_i^{1+2}) is nonconstant on A_1 × ⋯ × A_t whenever
Σ_i (d_i − 2) > n.)

*Proof.* WLOG d_1 = max_i d_i. In Lemma 3.13 take the row sets
𝒜_i = {0, …, e_i}, so α_i = Σ_{j∈𝒜_i} 2^j as required there. For subsets
B_i ⊆ {0, …, n−1} with |B_i| = e_i + 1, write β_i = Σ_{k∈B_i} 2^k (so
wt_2(β_i) = e_i + 1) and M = ∏_{i=1}^t X_i^{β_i}. By the second expression
in Lemma 3.13's formula (8), applied to the linearized parts of the Q_i as
justified in the next paragraph, the coefficient c_M of M in R is

    c_M = Σ_{r=0}^{n−1} ( µ · ∏_{i=1}^t Perm( M_i[𝒜_i, B_i − r] ) )^{2^r} ,    (14)

where B_i − r = { (k − r) mod n : k ∈ B_i } and M_i is the associated
matrix of the linearized part of Q_i.

*The constant terms.* Lemma 3.13 is stated for linearized polynomials,
while the Q_i here are affine linearized: Q_i = Q_i^{lin} + â_i. But c_M
is the same for the Q_i and for their linearized parts: expanding
(Q_i)^{α_i} = ∏_{j=0}^{e_i} (Q_i(X_i))^{2^j} with
(Q_i(X_i))^{2^j} = (Q_i^{lin}(X_i))^{2^j} + â_i^{2^j} (Freshman's dream),
a monomial of F_2-degree e_i + 1 in X_i — the X_i-part of M — arises only
by choosing the linearized term in every one of the e_i + 1 factors:
choosing a constant term â_i^{2^j} in some factor contributes F_2-degree 0
in X_i, the remaining e_i factors contribute at most e_i in total, and
reduction modulo X_i^{2^n} − X_i never increases F_2-degree. Hence
Lemma 3.13 applies verbatim with M_i as defined above. Following
Remark 3.14 we choose the sets B_i so that the r = 0 summand is nonzero —
each M_i′[𝒜_i, B_i] lower triangular with nonzero diagonal, hence each
permanent in the r = 0 product nonzero — while for every r ∈ {1, …, n−1}
at least one factor Perm(M_i[𝒜_i, B_i − r]) vanishes because
M_i′[𝒜_i, B_i − r] has an all-zero column. It therefore suffices to
construct B_i satisfying Claims 6.2 and 6.3 below.

Write Q_i(X_i) = Σ_{j=0}^{n−1} a_{i,j} X_i^{2^j} + â_i, and a′_{i,j} for the
zero-one indicator of a_{i,j}. Two support facts from §3 are used
repeatedly: a_{i,0} ≠ 0 and a_{i,n−d_i} ≠ 0, while a_{i,j} = 0 for
n − d_i + 1 ≤ j ≤ n − 1 (Lemma 3.10(1) and the degree bound
deg Q_i = 2^{n−d_i}).

**Choice of the sets B_i.** Define offsets ℓ_1 = 0 and, for 1 ≤ i < t,

    ℓ_{i+1} = ( ℓ_i + d_i − (e_i + 1) ) mod n ,
    i.e.  ℓ_i = ( Σ_{i′<i} ( d_{i′} − (e_{i′} + 1) ) ) mod n .

Let t_0 be the minimal index in {1, …, t} with

    S_{t_0+1} := Σ_{i′≤t_0} ( d_{i′} − (e_{i′} + 1) ) ≥ n − d_1 + 1 ;

such a t_0 exists by the hypothesis of the theorem, since
S_{t+1} = Σ_{i=1}^t (d_i − (e_i+1)) > n − d_1. For every i ≤ t_0,
minimality of t_0 gives S_i ≤ n − d_1 (for i = 1 this is S_1 = 0), hence
no mod-n reduction occurs in ℓ_i for i ≤ t_0: ℓ_i = S_i as a plain
integer, and

    ℓ_i ≤ n − d_1 ≤ n − d_i    (i ≤ t_0) .

The sets B_i are used substantively only for i ≤ t_0 (Claim 6.3), plus a
default choice for i > t_0 that needs nothing beyond Claim 6.2. Fix
i ≤ t_0. There are two cases, according to the coefficient a_{i,ℓ_i}.

*Case a_{i,ℓ_i} ≠ 0.* Set

    B_i = { ℓ_i } ∪ { n − d_i + 1, …, n − d_i + e_i } .

This has 1 + e_i elements: the top indices are e_i distinct residues,
none equal to ℓ_i since ℓ_i ≤ n − d_i.

*Case a_{i,ℓ_i} = 0.* Since a_{i,0} ≠ 0 we have ℓ_i ≥ 1, and since
a_{i,n−d_i} ≠ 0 and ℓ_i ≤ n − d_i we must have ℓ_i < n − d_i. Let

    j_i = min { j > ℓ_i : a_{i,j} ≠ 0 } ,   j̄_i = max { j < ℓ_i : a_{i,j} ≠ 0 } ,

which therefore exist and satisfy j̄_i < ℓ_i < j_i ≤ n − d_i. Let

    q_i = j_i − j̄_i − 1 ,   s_i = min { q_i, e_i } ;

thus a_{i,j} = 0 for j̄_i < j < j_i (a zero-gap of length q_i around ℓ_i),
and a_{i,j̄_i}, a_{i,j_i} are nonzero. Set

    B_i = { j_i } ∪ { j̄_i + 1, …, j̄_i + s_i } ∪ { n − d_i + s_i + 1, …, n − d_i + e_i } .

The three parts have 1, s_i and e_i − s_i elements respectively, totalling
e_i + 1; the last part is empty if s_i = e_i. They are pairwise disjoint:
j̄_i + s_i ≤ j̄_i + q_i = j_i − 1 < j_i ≤ n − d_i < n − d_i + s_i + 1.

*Remaining indices i > t_0 (if any).* Set
B_i = {0} ∪ { n − d_i + 1, …, n − d_i + e_i }; this is the first case's set
with ℓ_i replaced by 0, and a_{i,0} ≠ 0, so Claim 6.2 applies verbatim. No
covering property is needed from these i.

**Claim 6.2 (BK Claim 6.2).** For every i = 1, …, t the matrix
M_i′[𝒜_i, B_i] is, after reordering columns, lower triangular with nonzero
diagonal entries.

*Proof.* The matrix M_i′[𝒜_i, B_i] is the (e_i + 1) × (e_i + 1) minor with
rows 0, …, e_i and columns B_i; its (u, k)-entry is a′_{i, (k−u) mod n}.
Since e_i ≤ d_i − 2 < n, all indices below are read as the displayed
integers in {0, …, n−1}. Fix i and abbreviate a′_j = a′_{i,j}, d = d_i,
e = e_i.

*Case a′_{ℓ_i} = 1* (including every i > t_0, with 0 in place of ℓ_i).
Order the columns as (ℓ_i, n−d+1, …, n−d+e). The diagonal entries, at
(u, u), have coefficient indices: ℓ_i − 0 = ℓ_i, with a′_{ℓ_i} = 1; and, for
the column n−d+m at position m (1 ≤ m ≤ e), (n−d+m) − m = n−d, with
a′_{n−d} = 1 by Lemma 3.10(1). The entries above the diagonal are: in row 0,
the indices n−d+1, …, n−d+e; in row u ≥ 1 and column n−d+m with m > u, the
index n−d+(m−u) with 1 ≤ m−u ≤ e−1. All of these lie in
{n−d+1, …, n−d+e} ⊆ {n−d+1, …, n−1} (as e < d ≤ n), where a′ vanishes by
the degree bound. Hence the matrix is lower triangular with nonzero
diagonal.

*Case a′_{ℓ_i} = 0.* Order the columns as
(j_i, j̄_i + 1, …, j̄_i + s_i, n − d + s_i + 1, …, n − d + e). The first
s_i + 1 columns form the block

        ⎛ a′_{j_i}     a′_{j̄_i+1}   ⋯   a′_{j̄_i+s_i}  ⎞
        ⎜ a′_{j_i−1}   a′_{j̄_i}     ⋯   a′_{j̄_i+s_i−1}⎟
        ⎜     ⋮           ⋮         ⋱        ⋮        ⎟
        ⎝ a′_{j_i−s_i} a′_{j̄_i+1−s_i} ⋯   a′_{j̄_i}    ⎠

The (u, 0)-entry has index j_i − u; the (u, m)-entry (1 ≤ m ≤ s_i) has
index (j̄_i + m) − u. Above the diagonal: row 0 has indices
j̄_i + 1, …, j̄_i + s_i, all in the zero-gap (j̄_i, j_i), since
j̄_i + s_i ≤ j̄_i + q_i = j_i − 1; row u ≥ 1, column m > u has index
j̄_i + (m − u) with 1 ≤ m − u ≤ s_i − 1, likewise in the gap. So the
block's above-diagonal entries vanish. Its diagonal entries are
a′_{j_i} = 1 (at (0,0)) and a′_{j̄_i} = 1 (at (m, m), 1 ≤ m ≤ s_i). The
remaining columns n − d + s_i + 1, …, n − d + e (if any) sit at
positions s_i + 1, …, e; the column n − d + s_i + 1 + m′ at position
s_i + 1 + m′ has diagonal index
(n − d + s_i + 1 + m′) − (s_i + 1 + m′) = n − d, hence diagonal entry
a′_{n−d} = 1, and above-diagonal indices in
{n−d+1, …, n−d+e} ⊆ {n−d+1, …, n−1}, where a′ vanishes. Block-column
entries in rows below the block are below the diagonal, hence
unrestricted. So M_i′[𝒜_i, B_i] is lower triangular with nonzero
diagonal in this case too. ∎

**Proposition 6.4 (BK Proposition 6.4).** Let k ∈ B_i and

    r ∈ [k + 1, k + d_i − (e_i + 1)] mod n .

Then M_i′[𝒜_i, B_i − r] contains an all-zero column.

*Proof.* The column of M_i′[𝒜_i, B_i − r] indexed by k is the column vector

    ( a′_{i,k−r}, a′_{i,k−(r+1)}, …, a′_{i,k−(r+e_i)} )^T ,

with indices mod n: the (u, k−r)-entry of M_i′ is a′_{i,(k−r)−u}. Hence if

    [k − (r + e_i), k − r] mod n ⊆ [n − d_i + 1, n − 1] ,    (15)

every entry of this column is a′_{i,j} with j in the top zero window, and
the column is all-zero. Now r ∈ [k+1, k + d_i − (e_i+1)] means, for integer
representatives, k − r ∈ [−(d_i − e_i − 1), −1] = [1 − d_i + e_i, −1].
Therefore [k − (r + e_i), k − r] ⊆ [1 − d_i, −1] as ordinary intervals,
and reducing mod n sends [1 − d_i, −1] onto [n − d_i + 1, n − 1]. So (15)
holds. ∎

**Claim 6.3 (BK Claim 6.3).** For every r ∈ {1, …, n−1} there exists
i ∈ {1, …, t} such that M_i′[𝒜_i, B_i − r] contains an all-zero column.

*Proof.* By Proposition 6.4 it suffices to show the covering

    [1, n−1] ⊆ ∪_{i=1}^t ∪_{k∈B_i} [k + 1, k + d_i − (e_i + 1)] mod n .    (16)

(All intervals in this proof are intervals of residues mod n; for
i ≤ t_0 the offsets ℓ_i = S_i are plain integers < n, and we write
[ℓ_i + 1, ℓ_{i+1}] for the cyclic interval from ℓ_i + 1 to
ℓ_{i+1} = (S_{i+1}) mod n.)

For i = 1: ℓ_1 = 0 and a_{1,0} ≠ 0, so
B_1 = {0} ∪ [n − d_1 + 1, n − d_1 + e_1]. Proposition 6.4 with k = 0 gives
the interval [1, d_1 − (e_1 + 1)] = [ℓ_1 + 1, ℓ_2], and with k ranging over
the top part of B_1 the intervals [k+1, k + d_1 − (e_1 + 1)] chain
(consecutive starts, common length d_1 − e_1 − 1 ≥ 1) into
[n − d_1 + 2, (n − d_1 + e_1) + d_1 − (e_1 + 1)] = [n − d_1 + 2, n − 1].
Hence M_1′[𝒜_1, B_1 − r] has a zero column for every

    r ∈ [1, d_1 − (e_1 + 1)] ∪ [n − d_1 + 2, n − 1]
      = [ℓ_1 + 1, ℓ_2] ∪ [n − d_1 + 2, n − 1] .    (17)

We show that for 1 < i ≤ t_0,

    ∪_{k∈B_i} [k + 1, k + d_i − (e_i + 1)] ⊇ [ℓ_i + 1, ℓ_{i+1}] .    (18)

Together with (17) this proves (16): if S_{t_0+1} < n, then
ℓ_{t_0+1} = S_{t_0+1} ≥ n − d_1 + 1, so the chained intervals
[ℓ_1 + 1, ℓ_{t_0+1}] already reach n − d_1 + 1 and (17) covers the rest;
if S_{t_0+1} ≥ n, the chained intervals wrap past n − 1 and cover
[1, n−1] by themselves. Three cases establish (18).

*Case a_{i,ℓ_i} ≠ 0.* Then ℓ_i ∈ B_i, and Proposition 6.4 applied to
k = ℓ_i gives the interval [ℓ_i + 1, ℓ_i + d_i − (e_i + 1)], which is
[ℓ_i + 1, ℓ_{i+1}] (cyclically, since ℓ_{i+1} = (ℓ_i + d_i − (e_i+1)) mod n
for i ≤ t_0).

*Case a_{i,ℓ_i} = 0 and q_i < e_i.* Here s_i = q_i and, substituting
j_i = j̄_i + q_i + 1 and reordering,

    B_i = { j_i } ∪ { j̄_i + 1, …, j̄_i + q_i } ∪ { top part }
        = [j̄_i + 1, j_i] ∪ { n − d_i + q_i + 1, …, n − d_i + e_i } ,

because j̄_i + q_i + 1 = j_i makes the first two parts consecutive. Since
j̄_i < ℓ_i < j_i, we have ℓ_i ∈ [j̄_i + 1, j_i − 1] ⊆ B_i, and the previous
case's argument applies.

*Case a_{i,ℓ_i} = 0 and q_i ≥ e_i.* Here s_i = e_i, the top part of B_i is
empty, and

    B_i = { j_i } ∪ { j̄_i + 1, …, j̄_i + e_i } ,   j_i = j̄_i + q_i + 1 .

This is the step that uses the subfield hypothesis. Since i ≤ t_0 we have
ℓ_i ≤ n − d_1 ≤ n − d_i, and j̄_i < ℓ_i, so j̄_i + 1 ≤ n − d_i; in
particular j̄_i + 1 ≠ n − d_i + 1. Apply Lemma 3.10(2) to the coefficient
sequence of Q_i (i.e. to the image subspace of the direction of A_i;
the hypothesis consumed is the direction form — that subspace is not
contained in a constant multiple of a proper subfield. In the
application of Theorem 6.1 in §7 the distinction is vacuous: there n
is prime, the only proper subfield is F_2, and the linear spaces Y_i
there satisfy |Y_i| ≥ 2^{r+2}, so they are contained neither in an
affine shift of F_2 nor in a constant multiple of F_2, both of which
have 2 elements): the zero run
j̄_i + 1, …, j_i − 1 of length q_i starts at j̄_i + 1 ≠ n − d_i + 1, so it
cannot contain d_i − 1 consecutive indices; that is, q_i ≤ d_i − 2, or
equivalently

    j_i ≤ j̄_i + d_i − 1 .

Now take the union over the block part of B_i (all of B_i except j_i):
the intervals [k+1, k + d_i − (e_i + 1)] for k = j̄_i + 1, …, j̄_i + e_i
have consecutive starts and common length d_i − e_i − 1 ≥ 1, so they chain
into

    [j̄_i + 2, (j̄_i + e_i) + d_i − (e_i + 1)] = [j̄_i + 2, j̄_i + d_i − 1] ,

which contains j_i by the inequality just proved, and satisfies
j̄_i + 2 ≤ ℓ_i + 1 since j̄_i < ℓ_i. Reinserting k = j_i, whose interval
[j_i + 1, j_i + d_i − (e_i + 1)] starts at j_i + 1 ≤ (j̄_i + d_i − 1) + 1,
the union extends without a gap to j_i + d_i − (e_i + 1) ≥
ℓ_i + d_i − (e_i + 1). Hence

    ∪_{k∈B_i} [k + 1, k + d_i − (e_i + 1)] ⊇ [ℓ_i + 1, ℓ_i + d_i − (e_i + 1)]
                                            = [ℓ_i + 1, ℓ_{i+1}] ,

proving (18) in the last case as well. ∎

*Conclusion of the proof of Theorem 6.1.* Claims 6.2 and 6.3 are exactly
the two hypotheses of Remark 3.14 (lower triangular in place of upper
triangular is immaterial, by reversing the order of rows and columns):
the r = 0 summand of (14) is a product of permanents of triangular
matrices with nonzero diagonal, hence nonzero, and every summand with
r ∈ {1, …, n−1} contains a factor Perm(M_i[𝒜_i, B_i − r]) = 0 coming from
an all-zero column of M_i′[𝒜_i, B_i − r]. Therefore c_M ≠ 0: the monomial
M = ∏_{i=1}^t X_i^{β_i}, with β_i = Σ_{k∈B_i} 2^k and wt_2(β_i) = |B_i| =
e_i + 1, has nonzero coefficient in R. By Proposition 3.1, R is not
constant on F^t, i.e. |R(F^t)| > 1. ∎

## 7. Affine spaces of sublinear dimension: Lemmas 7.1, 7.2, Theorem 7.3, and Theorem 2.2

We now prove Theorem 2.2. Throughout this section p = 2; n denotes the
*block length* and is prime, F = F_{2^n}, r is the number of blocks, t < r,
and Tr(Y) = Σ_{s=0}^{n−1} Y^{2^s} is the trace F → F_2. We use from §3:
Proposition 3.1 (reduction modulo ⟨X_i^{2^n} − X_i⟩ preserves values on
F^r, and a reduced polynomial taking more than one value has a
positive-degree monomial), Proposition 3.3 (every F_2-linear π: F → F_2 is
π(x) = Tr(µx) for a unique µ, and π is nontrivial iff µ ≠ 0), Lemma 3.4
(F_2-linear maps F → F are given by linearized polynomials), Lemma 3.6 and
Lemma 3.7 (image-subspace polynomials Q_L: monic affine linearized,
Q_L(F) = L, linearized when L is linear); and from §6, Theorem 6.1.
All polynomials in X_1, …, X_r are reduced modulo ⟨X_i^{2^n} − X_i⟩; this
reduction does not increase F_2-degrees (reducing an exponent replaces
2^n by 1 and binary carries only decrease digit sums), and a linearized
polynomial remains a sum of monomials a·X_i^{2^j}, of F_2-degree ≤ 1 in
its variable, after reduction and after any Frobenius power.

**Lemma 7.1 (Bourgain's decomposition).** Let A ⊆ (F_{2^n})^r be an
F_2-affine subspace and γ ∈ A. Then there exist linear subspaces
Y_1, …, Y_r ⊆ F_{2^n} and F_2-linear maps σ_{ij}: Y_j → F_{2^n} (j < i)
such that

    A = { (x_1,…,x_r) : x_i = γ_i + y_i + Σ_{j<i} σ_{ij}(y_j),  y_i ∈ Y_i },

and dim A = Σ_{i∈[r]} dim Y_i.

*Proof.* Write A = γ + V with V linear. For 1 ≤ i ≤ r let
U_i = {v ∈ V : v_1 = … = v_{i−1} = 0} (U_1 = V) and
Y_i = {v_i : v ∈ U_i} ⊆ F_{2^n}, a linear subspace; the coordinate map
ε_i: U_i → Y_i, v ↦ v_i, is linear and surjective. Choose an F_2-linear
section s_i: Y_i → U_i with ε_i ∘ s_i = id, and for j < i set
σ_{ij}(y) = (s_j(y))_i ∈ F_{2^n}. Every v ∈ V is uniquely
v = Σ_i s_i(y_i) with y_i ∈ Y_i: put y_1 = v_1 ∈ Y_1, so
v − s_1(y_1) ∈ U_2; put y_2 = (v − s_1(y_1))_2 ∈ Y_2, so
v − s_1(y_1) − s_2(y_2) ∈ U_3; and continue. Uniqueness follows by reading
the coordinates in order: if Σ_i s_i(y_i) = 0, the first coordinate gives
y_1 = 0 (for i ≥ 2, s_i(y_i) ∈ U_i has first coordinate 0, and
(s_1(y_1))_1 = y_1), then inductively every y_i = 0. Hence
(y_1, …, y_r) ↦ Σ_i s_i(y_i) is an isomorphism ⊕_i Y_i → V and
dim V = Σ_i dim Y_i. The i-th coordinate of v = Σ_j s_j(y_j) is
v_i = y_i + Σ_{j<i} σ_{ij}(y_j): indeed (s_i(y_i))_i = y_i,
(s_j(y_j))_i = σ_{ij}(y_j) for j < i, and (s_j(y_j))_i = 0 for j > i
because s_j(y_j) ∈ U_j vanishes in all coordinates < j. Adding γ gives
the stated form of A. □ (This is the row-echelon form of a basis of V.)

**Lemma 7.2 (coefficient invariance under couplings and shifts).** For
each i ∈ [r] let P_i(X_i) ∈ F_{2^n}[X_i] be a linearized polynomial; for
each j < i let P_{ij}(X_j) ∈ F_{2^n}[X_j] be a linearized polynomial; let
γ = (γ_1, …, γ_r) ∈ (F_{2^n})^r. Let I_0 = {i_1 < i_2 < … < i_t} ⊆ [r],
let e_i ≥ 1 be integers for i ∈ I_0, α_i = Σ_{k=0}^{e_i} 2^k, and
µ ∈ F_{2^n} \ {0}. Let

    g  = Tr( µ · Π_{i∈I_0} ( P_i(X_i) + Σ_{j<i} P_{ij}(X_j) + γ_i )^{α_i} )
             mod ⟨X_i^{2^n} − X_i⟩,
    g′ = Tr( µ · Π_{i∈I_0} P_i(X_i)^{α_i} )  mod ⟨X_i^{2^n} − X_i⟩.

Then for every tuple (β_i)_{i∈I_0} with wt_2(β_i) = e_i + 1, the
coefficients of the monomial M = Π_{i∈I_0} X_i^{β_i} in g and in g′ are
equal.

*Proof.* Let X_{<i} = (X_1, …, X_{i−1}) and
P̂_i(X_{<i}) = Σ_{j<i} P_{ij}(X_j) + γ_i, which does not involve X_i.
Expanding the trace, Tr(z) = Σ_{s=0}^{n−1} z^{2^s},

    g = Σ_{s=0}^{n−1} ( µ · Π_{i∈I_0} (P_i(X_i) + P̂_i(X_{<i}))^{α_i} )^{2^s}
            mod ⟨X_i^{2^n} − X_i⟩                                        (19)
      = Σ_{s=0}^{n−1} µ^{2^s} · Π_{i∈I_0} Π_{l=0}^{e_i}
            ( P_i(X_i) + P̂_i(X_{<i}) )^{2^{s+l}}  mod ⟨X_i^{2^n} − X_i⟩, (20)

where (20) uses (z^{α_i})^{2^s} = Π_{l=0}^{e_i} z^{2^{s+l}}, since
α_i = Σ_{l=0}^{e_i} 2^l. Fix s and expand one factor of (20). By the
Freshman's dream in characteristic 2,

    Π_{l=0}^{e_i} ( P_i(X_i) + P̂_i(X_{<i}) )^{2^{s+l}}
        = Π_{l=0}^{e_i} ( P_i(X_i)^{2^{s+l}} + P̂_i(X_{<i})^{2^{s+l}} ),  (21)–(22)

and choosing the P̂-term exactly at the positions of a subset
L ⊆ {0, …, e_i} gives

    = ( Π_{l=0}^{e_i} P_i(X_i)^{2^{s+l}} )
        + Σ_{L≠∅} Π_{l∉L} P_i(X_i)^{2^{s+l}} · Π_{l′∈L} P̂_i(X_{<i})^{2^{s+l′}}.  (23)

Each P_i(X_i)^{2^{s+l}} is a linearized polynomial in X_i, hence a sum of
monomials of F_2-degree 1 in X_i; each P̂_i(X_{<i})^{2^{s+l′}} does not
involve X_i at all. So the first term of (23) has F_2-degree at most
e_i + 1 in X_i (a product of e_i + 1 such factors), while every term with
L ≠ ∅ has F_2-degree at most e_i + 1 − |L| ≤ e_i in X_i: only its
e_i + 1 − |L| factors P_i(X_i)^{2^{s+l}} involve X_i, each to
F_2-degree ≤ 1. Summarize (23) as

    (P_i(X_i))^{Σ_{l=0}^{e_i} 2^{s+l}} + G_i(X_{≤i}),

where G_i has F_2-degree at most e_i in X_i (and involves only
X_1, …, X_i, since P̂_i involves only X_{<i}).

Now fix M, with wt_2(β_i) = e_i + 1 for i ∈ I_0, and look at the s-th
summand of (20), which is µ^{2^s} times the product over i ∈ I_0 of the
brackets above. Peel off the largest index i_t:

    Π_{i∈I_0} [ (P_i)^{Σ_l 2^{s+l}} + G_i ]
        = [ (P_{i_t}(X_{i_t}))^{Σ_l 2^{s+l}} + G_{i_t}(X_{≤i_t}) ]
              · H_t(X_{≤ i_t − 1}),

where H_t, the product over i ∈ I_0 with i < i_t, does not mention
X_{i_t}. M has F_2-degree e_{i_t} + 1 in X_{i_t}; H_t has F_2-degree 0
in X_{i_t} and G_{i_t} has F_2-degree at most e_{i_t} in X_{i_t}, so no
monomial of G_{i_t} · H_t contributes to the coefficient of M: only the
pure term (P_{i_t}(X_{i_t}))^{Σ_l 2^{s+l}} can. The same peeling
applies at i_{t−1}: writing
H_t = [ (P_{i_{t−1}})^{Σ_l 2^{s+l}} + G_{i_{t−1}}(X_{≤i_{t−1}}) ]
· H_{t−1}(X_{≤ i_{t−2}}), no other factor mentions X_{i_{t−1}} —
H_{t−1} uses only X_{≤ i_{t−2}}, the brackets for i < i_{t−1} only
X_{≤i}, and the surviving pure term is a polynomial in X_{i_t} alone —
while M has F_2-degree e_{i_{t−1}} + 1 in X_{i_{t−1}} and G_{i_{t−1}}
at most e_{i_{t−1}}, so again only the pure term of the i_{t−1}-bracket
contributes. Continuing down through i_{t−2}, …, i_1, the
coefficient of M in the s-th summand of (20) equals its coefficient in

    µ^{2^s} · Π_{i∈I_0} (P_i(X_i))^{Σ_{l=0}^{e_i} 2^{s+l}}
        = ( µ · Π_{i∈I_0} P_i(X_i)^{α_i} )^{2^s}.

Summing over s = 0, …, n−1 gives exactly the coefficient of M in g′. □

**Theorem 7.3 (affine disperser, non-parameterized).** Let t < r be
integers and let n be prime with n ≥ r(r+1)/t. For each i ∈ [r] let
e_i = r + 1 − i and α_i = Σ_{k=0}^{e_i} 2^k. Let µ ∈ F_{2^n} \ {0} and

    f(x_1, …, x_r) = Tr( µ · Σ_{I⊆[r], |I|=t} Π_{i∈I} x_i^{α_i} ).

Then every affine A ⊆ (F_{2^n})^r with

    dim A  >  nr/t + nt + r(r+1)

satisfies |f(A)| > 1.

*Proof.* (1) *Parameterization of A.* Fix γ ∈ A and decompose A by
Lemma 7.1, with linear Y_i and linear σ_{ij}: Y_j → F_{2^n}. Let
Q_i ∈ F_{2^n}[X_i] be the image-subspace polynomial of Y_i (Lemma 3.7);
Q_i is linearized because Y_i is linear (Lemmas 3.6, 3.7). For j < i,
extend σ_{ij} arbitrarily to an F_2-linear map on F_{2^n} and let
Q_{ij} ∈ F_{2^n}[X_j] be the linearized polynomial of the F_2-linear map
x ↦ σ_{ij}(Q_j(x)) (Lemma 3.4); then Q_{ij}(x) = σ_{ij}(Q_j(x)) for all
x, independently of the extension, since Q_j(x) ∈ Y_j, the domain of
σ_{ij}. Set

    R_i(X_1, …, X_r) = Q_i(X_i) + Σ_{j<i} Q_{ij}(X_j) + γ_i,
    H = (R_1, …, R_r).

As X_j ranges over F_{2^n}, Q_j(X_j) ranges over all of Y_j (Lemma 3.7),
and the couplings track the same y_j, so H((F_{2^n})^r) = A. Let
h = f ∘ H mod ⟨X_i^{2^n} − X_i⟩, i.e.

    h = Σ_{I⊆[r], |I|=t} Tr( µ · Π_{i∈I} R_i^{α_i} )
            mod ⟨X_i^{2^n} − X_i⟩.                                   (24)–(25)

By Proposition 3.1, h((F_{2^n})^r) = f(A); it suffices to find a
positive-degree monomial with nonzero coefficient in h.

(2) *High-entropy blocks.* Let d_i = dim Y_i, so
Σ_i d_i = dim A > nr/t + nt + r(r+1), and let
S = {i : d_i > r+1}. Then |S| ≥ t: otherwise
Σ_i d_i ≤ |S|·n + r(r+1) ≤ (t−1)n + r(r+1) < nt + r(r+1) + nr/t
(each d_i ≤ n, and each of the at most r indices outside S has
d_i ≤ r+1). Also
Σ_{i∈S} d_i ≥ Σ_{i∈[r]} (d_i − (r+1)) = dim A − r(r+1) > nr/t + nt
(the terms dropped, for i ∉ S, are ≤ 0). Averaging over the t-subsets
I_0 of S, some I_0 satisfies

    Σ_{i∈I_0} (d_i − (r+1)) ≥ (t/|S|) · Σ_{i∈S} d_i − (r+1)t
        ≥ (t/r)(nr/t + nt) − (r+1)t = n + nt²/r − (r+1)t ≥ n,        (26)

the last step by the hypothesis n ≥ r(r+1)/t, i.e. nt/r ≥ r+1.

(3) *The I_0-summand.* Let g be the I_0-term of (24)–(25) and
g′ = Tr(µ · Π_{i∈I_0} Q_i(X_i)^{α_i}) mod ⟨X_i^{2^n} − X_i⟩. Apply
Theorem 6.1 to the spaces Y_i, i ∈ I_0, with these Q_i, this µ, and
e_i = r + 1 − i. Its hypotheses hold:

- d_i > r + 1 ≥ 3, in particular d_i > 1 (since i ∈ S and t ≥ 1);
- 1 ≤ e_i < d_i − 1: e_i = r + 1 − i ≥ 1 as i ≤ r, and
  e_i + 1 ≤ r + 1 < d_i gives e_i ≤ d_i − 2;
- Σ_{i∈I_0} (d_i − (e_i + 1)) ≥ Σ_{i∈I_0} (d_i − (r+1)) ≥ n by (26)
  (using e_i + 1 ≤ r + 1), and n > n − max_{i∈I_0} d_i + 1 because
  max_{i∈I_0} d_i > r + 1 ≥ 3 — so the required sum inequality holds;
- no Y_i is contained in an affine shift of a proper subfield of
  F_{2^n}: this is **where primality of n is used**. The subfields of
  F_{2^n} are F_{2^d} for d | n; n prime leaves only F_2 and F_{2^n}.
  An affine shift of F_2 has exactly 2 elements, while
  |Y_i| = 2^{d_i} ≥ 2^{r+2}, so containment is impossible.

Theorem 6.1 therefore gives a monomial
M = Π_{i∈I_0} X_i^{β_i} with wt_2(β_i) = e_i + 1 and nonzero coefficient
in g′. By Lemma 7.2, applied with P_i = Q_i (linearized), P_{ij} = Q_{ij}
and the shifts γ_i, the coefficient of M in g equals that in g′; in
particular it is nonzero.

(4) *Uniqueness.* We show no other summand of (24)–(25) contributes M.
Each R_j is a sum of monomials a·X_i^{2^b} with 1 ≤ i ≤ j, plus possibly
a constant term — a polynomial of F_2-degree at most 1 involving only
X_1, …, X_j. Let J ⊆ [r], |J| = t, J ≠ I_0, and suppose M had nonzero
coefficient in the J-summand Tr(µ · Π_{j∈J} R_j^{α_j}). Expanding the
trace, M must occur in Π_{j∈J} Π_{l=0}^{e_j} R_j^{2^{l+w}} for some
w ∈ {0, …, n−1}. Each factor R_j^{2^{l+w}} is again a sum of monomials
of F_2-degree ≤ 1 in X_1, …, X_j plus a constant. Choose, for each
(j, l), one monomial of the factor whose product is M; if the chosen
monomial is in the variable X_i (so i ≤ j), set λ_j(l) = i; if the
constant term was chosen at (j, l), set λ_j(l) := j — this only enlarges
the counts |λ_j^{−1}(i)| for i ∈ I_0 in what follows, preserves
λ_j(l) ≤ j, and keeps Σ_{i∈I_0} |λ_j^{−1}(i)| ≤ r + 2 − j, so Claim 7.4
applies to the extended functions.
The F_2-degree of M in X_i is at most the number of (j, l) with
λ_j(l) = i — each such factor contributes one power X_i^{2^b}, and
binary carries can only decrease the digit sum of the total exponent —
but for i ∈ I_0 that degree equals e_i + 1 = r + 2 − i. Hence

    Σ_{j∈J} |λ_j^{−1}(i)| ≥ r + 2 − i      for every i ∈ I_0.        (27)

**Claim 7.4.** Let I_0 ⊆ [r], |I_0| = t. Suppose J ⊆ [r], |J| = t, and
there exist functions λ_j: {0, …, r+1−j} → {1, …, j} for j ∈ J
satisfying (27). Then J = I_0.

*Proof of Claim 7.4.* For each j ∈ J,
Σ_{i∈I_0} |λ_j^{−1}(i)| ≤ |{0, …, r+1−j}| = r + 2 − j. Summing (27) over
i ∈ I_0 and then over j gives, using |I_0| = |J|,
Σ_{i∈I_0} i ≥ Σ_{j∈J} j. For general k ≥ 0, the expression
Σ_{j∈J} Σ_{l=0}^{r+1−j} (j^k − λ_j(l)^k) is ≥ 0 since λ_j(l) ≤ j;
expanding it and lower-bounding
Σ_{j} Σ_l λ_j(l)^k = Σ_i i^k·|{(j,l) : λ_j(l) = i}| ≥
Σ_{i∈I_0} (r+2−i)·i^k by (27) yields
Σ_{i∈I_0} i^{k+1} − Σ_{j∈J} j^{k+1} ≥ (r+2)·( Σ_{i∈I_0} i^k − Σ_{j∈J} j^k ),
which is ≥ 0 by induction on k (the case k = 0 is equality of
cardinalities). Hence for all k,

    Σ_{i∈I_0} i^k ≥ Σ_{j∈J} j^k.                                     (28)

Let i_1 = max I_0 and j_1 = max J. By (28) with k → ∞ (if j_1 > i_1 then
Σ_{j∈J} j^k ≥ j_1^k > t·i_1^k ≥ Σ_{i∈I_0} i^k for large k), i_1 ≥ j_1.
But i_1 ≤ j_1: by (27) some λ_j(l) equals i_1 (the count is
≥ r + 2 − i_1 ≥ 1), and λ_j(l) ≤ j ≤ j_1 forces j_1 ≥ i_1. So
i_1 = j_1, and (27) at i = i_1, j = j_1 then forces
|λ_{j_1}^{−1}(i_1)| ≥ r + 2 − j_1, the full number of l's; hence
λ_{j_1}(l) = i_1 for every l. Deleting i_1 from I_0 and j_1 from J
leaves (28) valid for the deleted sets, so the second-largest elements
satisfy i_2 ≥ j_2 by the same domination; and i_2 ≤ j_2 because no
remaining j can map to i_2 except one ≥ i_2, while j_1 is exhausted by
i_1 ≠ i_2. Thus i_2 = j_2. Inducting on the s-th largest elements gives
I_0 = J. □ (Claim 7.4)

So M, which has positive degree since wt_2(β_i) = e_i + 1 ≥ 2, appears
in h with the nonzero coefficient computed in step (3), contributed by
the I_0-summand alone. By Proposition 3.1, |h((F_{2^n})^r)| > 1, hence
|f(A)| > 1. □ (Theorem 7.3)

**Deduction of Theorem 2.2 from Theorem 7.3.** Fix the input length m
and set n = least prime > 2m^{3/5}, r = ⌈m/n⌉, t = ⌈√r⌉, as in
Theorem 2.2. By Bertrand's postulate there is a prime in
(2m^{3/5}, 4m^{3/5}], so

    2m^{3/5} < n ≤ 4m^{3/5}.

Hence m^{2/5}/4 ≤ m/n ≤ r ≤ m/n + 1 < m^{2/5}/2 + 1; in particular
r → ∞ with m, and t = ⌈√r⌉ < r once r ≥ 3. We verify the two hypotheses
of Theorem 7.3 and then its conclusion's threshold.

*Block-length hypothesis.* Since t ≥ √r,

    r(r+1)/t ≤ (r+1)·√r = (r+1)^{3/2} ≤ (m^{2/5}/2 + 2)^{3/2}
        = (1/(2√2)) · m^{3/5} · (1 + 4m^{−2/5})^{3/2}
        = (1/(2√2)) · m^{3/5} (1 + o(1)) < 2m^{3/5} < n

for all large m (the leading constant 1/(2√2) ≈ 0.354 is well below 2).
So n ≥ r(r+1)/t holds.

*Threshold.* We estimate nr/t + nt + r(r+1) *jointly*, using
r ≤ m/n + 1 rather than maximizing n and r independently. With t ≥ √r,

    nr/t ≤ n√r = (n²·(m/n + 1))^{1/2} = (nm + n²)^{1/2}
        ≤ (4m^{8/5} + 16m^{6/5})^{1/2} = 2m^{4/5}(1 + o(1)),
    nt ≤ n(√r + 1) ≤ 2m^{4/5}(1 + o(1)) + 4m^{3/5}
        = 2m^{4/5}(1 + o(1)),
    r(r+1) ≤ (m^{2/5}/2 + 1)(m^{2/5}/2 + 2) = (1/4)m^{4/5}(1 + o(1)).

Adding, nr/t + nt + r(r+1) ≤ (2 + 2 + 1/4)·m^{4/5}(1 + o(1))
= 4.25·m^{4/5}(1 + o(1)) < 6m^{4/5} for all large m. Thus every affine
A ⊆ F_2^m with dim A > 6m^{4/5} satisfies the dimension hypothesis of
Theorem 7.3 (after the identification of the ambient space below).

*Identification of the two f's.* By Proposition 3.3, the nontrivial
F_2-linear π of Theorem 2.2 is π(z) = Tr(µz) for a unique µ ≠ 0,
matching the form used in Theorem 7.3. The exponents differ only in
indexing order: equation (4) gives block i the exponent
Σ_{j=0}^{i} 2^j, while Theorem 7.3 gives block i the exponent
Σ_{k=0}^{e_i} 2^k with e_i = r + 1 − i, i.e. Σ_{j=0}^{r+1−i} 2^j. Under
the renaming of blocks i ↦ r + 1 − i, the two symmetric sums coincide:
f of Theorem 7.3 evaluated at (x_1, …, x_r) equals f of Theorem 2.2
evaluated at (x_r, …, x_1). Coordinate permutation preserves affinity
and dimension, so the conclusion of Theorem 7.3 transfers: f of
Theorem 2.2 is nonconstant on every affine subspace of dimension
> 6m^{4/5}, for all sufficiently large m. ∎ (Theorem 2.2)

*Application to the family of §1 (padded, all lengths).* Finally rename
(m, n) → (N, k) and apply this to the family (f_N) of §1 of the main
proof. There k is the least prime > 2N^{3/5}, r = ⌈N/k⌉, t = ⌈√r⌉, and
φ: F_2^N → K^r = (F_{2^k})^r is the pad-and-split map of §1 (append
rk − N zeros, split into r blocks of k bits, read each block in the
polynomial basis of K). This φ is injective and F_2-linear, and it is
admissible in Theorem 2.2 applied at length m := N: the theorem is
quantified over all such φ, so the short last block (rk − N < k padded
zeros) needs no divisibility assumption, and all parameters k, r, t are
computed from the true length N. The exponents agree
(Σ_{j=0}^{i} 2^j = 2^{i+1} − 1) and π, the constant coefficient, is a
nontrivial F_2-linear functional. Consequently, for every affine
A ⊆ F_2^N with dim A > 6N^{4/5}, the image φ(A) ⊆ (F_{2^k})^r is affine
of the *same* dimension (φ injective linear), and

    |f_N(A)| = |f(φ(A))| > 1

by Theorem 2.2. The estimates above are asymptotic, so this is proved
for all sufficiently large N; for N ≤ 7776 = 6^5 one has
6N^{4/5} ≥ N, no affine subspace of F_2^N has dimension > 6N^{4/5},
and the assertion is vacuous. Hence f_N is an affine disperser for
dimension d(N) = ⌊6N^{4/5}⌋ + 1 at every length where the property has
content, in particular for all sufficiently large N — which is all that
§2 of the main proof uses. This completes the proof of the disperser
property with no external premise beyond §3, §6 and the present §7. □
