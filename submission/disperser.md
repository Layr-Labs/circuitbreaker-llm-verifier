# Disperser: subspace polynomials, independent sources, Theorem 2.2 (BK §§3, 6, 7)

The §3, §6 and §7 material of Ben-Sasson–Kopparty, *Affine Dispersers
from Subspace Polynomials*, with proofs, specialized throughout to
p = 2 and prime block length. Statements keep BK's numbering.

**Main Lemma (consumed by `proof.md` §2).** Let N ≥ 1, k the least prime
> 2N^{3/5}, r = ⌈N/k⌉, t = ⌈√r⌉, and f_N the function of `proof.md` §1
over K = F_{2^k}: zero-pad to rk bits, split into blocks y_1, …, y_r ∈ K,
output the constant coefficient of Sym_r^t(y_1, …, y_r) with y_i raised
to 2^{i+1} − 1. Then f_N is nonconstant on every affine A ⊆ F_2^N with
dim A > 6N^{4/5}: an affine disperser for dimension
d(N) = ⌊6N^{4/5}⌋ + 1. For N ≤ 7776 = 6^5 this is vacuous, as
6N^{4/5} ≥ N. This is Theorem 2.2, proved at the end of §7.

## 3. Subspace polynomials

**Conventions.** p = 2; F = F_{2^n}, n prime (n = k in the application);
F̄ an algebraic closure; coefficient indices mod n.
Q(X) = Σ_{j=0}^{n-1} a_j X^{2^j} is *linearized*; with a constant term â
allowed, *affine linearized* — exactly the maps F → F of F_2-degree at
most 1. For affine L ⊆ F of dimension d, P_L (Definition 3.5) is its
kernel-subspace polynomial, Q_L (Lemma 3.7) its image-subspace
polynomial. wt_2(m) = sum of the binary digits of m, the *F_2-degree
weight*.

**Proposition 3.1.** For h ∈ F[X_1, …, X_r], let h′ = h mod
⟨X_i^{2^n} − X_i⟩ (unique representative of degree < 2^n in each
variable). Then h = h′ on F^r, and |h(F^r)| > 1 iff h′ has positive
degree.

*Proof.* a^{2^n} = a on F, so the replacements X_i^{2^n} ↦ X_i preserve
values. If h′ has positive degree then h′ − c ≠ 0 for every c, and a
nonzero polynomial of degree < 2^n in each variable cannot vanish on
F^r (r = 1: root bound; induct on r). Constant h′ gives constant h. □

**Proposition 3.3.** With Tr(Y) = Σ_{i=0}^{n-1} Y^{2^i}: F → F_2, every
F_2-linear π: F → F_2 is π(x) = Tr(µx) for a unique µ ∈ F; π is trivial
iff µ = 0.

*Proof.* (x, y) ↦ Tr(xy) is nondegenerate: Tr is a nonzero polynomial of
degree 2^{n-1} < 2^n, so Tr(z) ≠ 0 for some z, and for µ ≠ 0 the point
x = z/µ gives Tr(µx) ≠ 0. Hence µ ↦ Tr(µ·) is injective F_2-linear
between n-dimensional spaces, an isomorphism. □

**Lemma 3.4.** For an F_2-linear isomorphism φ: F_2^n → F, affine
transformations T: F_2^n → F_2^n correspond bijectively to affine
linearized polynomials P_T with P_T(φ(b)) = φ(T(b)).

*Proof.* Transport by φ to affine maps F → F. Affine linearized
polynomials induce F_2-affine maps, and both sides form F_2-spaces of
dimension n(n + 1) (coefficients a_0, …, a_{n-1}, â; resp. an F_2-matrix
plus a translation: n^2 + n). Evaluation is F_2-linear and injective:
if Σ a_j X^{2^j} + â vanishes on F then â = 0 (at 0) and the linearized
part, of degree ≤ 2^{n-1} < 2^n, has < 2^n roots, forcing all a_j = 0.
Equal finite dimensions give the bijection. □

**Definition 3.5 (kernel-subspace polynomial).** For affine L ⊆ F of
dimension d: P_L(X) = ∏_{α∈L} (X − α).

**Lemma 3.6.** P_L is monic affine linearized of degree 2^d, and P_L is
linearized iff L is linear.

*Proof.* Let V be linear of dimension d; induction on d. P_{ {0} } = X.
If W = V ⊕ ⟨b⟩, b ∉ V, then P_V(X + b) = P_V(X) + P_V(b) (P_V
linearized), so P_V(X)^2 + P_V(b)P_V(X) is monic linearized of degree
2^{d+1} vanishing exactly on V ∪ (V + b) = W — 2^{d+1} roots — hence
equals P_W; in particular ker P_V = V on F. For affine L = γ + V,
P_V(X) + P_V(γ) is monic affine linearized of degree 2^d with root set
{x : P_V(x) = P_V(γ)} = γ + V = L, hence equals P_L, with linearized
part P_V. It is linearized iff P_V(γ) = 0 iff γ ∈ ker P_V = V iff
0 ∈ L iff L is linear. □

**Lemma 3.7 (image-subspace polynomial).** If L ⊆ F is affine of
dimension d, there is a monic affine linearized Q_L(X) of degree
2^{n−d}, the image-subspace polynomial of L, with Q_L(F) = L and

    P_L(Q_L(X)) = X^{2^n} − X.                                        (6a)

If L is linear, Q_L is linearized (not merely affine linearized),
ker Q_L = P_L(F), Q_L has 2^{n−d} roots, and also

    Q_L(P_L(X)) = X^{2^n} − X.                                        (6b)

*Proof.* Let L be linear. L_0 := P_L(F) is linear of dimension n − d
(P_L linearized, kernel L). Set Q_L = P_{L_0} (Lemma 3.6): monic
linearized of degree 2^{n−d}. Q_L ∘ P_L is monic of degree 2^n and
vanishes on F (P_L(x) ∈ L_0, the root set of Q_L), so
Q_L ∘ P_L = ∏_{α∈F}(X−α) = X^{2^n} − X. With
g(Y) := P_L(Q_L(Y)) − (Y^{2^n} − Y) we get
g(P_L(X)) = P_L(X^{2^n} − X) − (P_L(X)^{2^n} − P_L(X)) = 0, since P_L
linearized gives P_L(X^{2^n} − X) = P_L(X^{2^n}) − P_L(X) and
P_L(X^{2^n}) = P_L(X)^{2^n} (coefficients lie in F, fixed by the
2^n-Frobenius). Since deg(g ∘ P_L) = (deg g)·2^d, g ≡ 0: (6a). Hence
image(Q_L) ⊆ ker P_L = L, and ker Q_L ⊇ L_0 with both of dimension
n − d, so |image Q_L| = 2^d = |L|, Q_L(F) = L, ker Q_L = L_0,
2^{n−d} roots, and (6b) is the identity just proved.
For affine L = γ + V, V its direction, set Q_L(X) = Q_V(X) + γ: monic
affine linearized of degree 2^{n−d}, Q_L(F) = V + γ = L, linearized part
Q_V. By Lemma 3.6's proof, P_L = P_V + P_V(γ) (characteristic 2), and
P_V linearized gives P_L(Q_L(Y)) = P_V(Q_V(Y) + γ) + P_V(γ)
= P_V(Q_V(Y)) = Y^{2^n} − Y: (6a). For affine L, (6b) holds only up to
the constant γ: Q_L(P_L(X)) = X^{2^n} − X + γ, by the same computation
and Q_V(P_V(γ)) = γ^{2^n} − γ = 0. All applications use (6a), the image
identity Q_L(F) = L, and — for linear L — linearizedness of Q_L. □

**Claim 3.8.** If a^{2^k} = ca with a, c ∈ F, k ≥ 1, and b ∈ F̄ is a
(2^k − 1)-th root of c, then a ∈ b·F_{2^k}.

*Proof.* a = 0 is clear. Else a^{2^k − 1} = c, so
(a/b)^{2^k − 1} = 1, i.e. (a/b)^{2^k} = a/b. □

**Claim 3.9.** For Q(X) = Σ_{j=0}^{n-1} a_j X^{2^j} + â and any integer s,

    (Q(X))^{2^s} mod (X^{2^n} − X) = Σ_j a_{j−s}^{2^s} X^{2^j} + â^{2^s},

indices mod n.

*Proof.* (Q(X))^{2^s} = Σ_j a_j^{2^s} X^{2^{j+s}} + â^{2^s}; reduction
replaces j + s by (j + s) mod n since X^{2^n} ≡ X; reindex. □

**Lemma 3.10 (main structural lemma).** Let L ⊆ F be linear of dimension
d, Q_L(X) = Σ_{j=0}^{n-1} a_j X^{2^j} its image-subspace polynomial.

1. Every d consecutive indices J = {(r + j) mod n : j = 0, …, d − 1}
   contain some j with a_j ≠ 0. In particular a_0 ≠ 0 and a_{n−d} ≠ 0.
2. If L is not contained in any constant multiple of a proper subfield
   of F (L ⊄ β·F_{2^k} for all β ∈ F, k < n), then every d − 1
   consecutive indices J = {(r + j) mod n : j = 0, …, d − 2} with
   r ≠ n − d + 1 contain some j with a_j ≠ 0.

*Proof.* If d = n then Q_L = X and both parts are immediate. Assume
d < n; then a_j = 0 for n − d + 1 ≤ j ≤ n − 1 (degree 2^{n−d}) and
a_{n−d} = 1 (monic).
Part 1: suppose a_j = 0 on J. Put s ≡ −(r + d) (mod n) and
Q′(X) = (Q_L(X))^{2^s} mod (X^{2^n} − X). By Claim 3.9 its
X^{2^j}-coefficient is a_{j+r+d}^{2^s}; for j ∈ [n − d, n − 1] the index
(j + r + d) mod n runs over J, so those coefficients vanish and
deg Q′ ≤ 2^{n−d−1}. Q′ ≠ 0 since a_{n−d} = 1 occurs among its
coefficients. But Q′ vanishes on ker Q_L
(Q′(α) = Q_L(α)^{2^s}), which has 2^{n−d} elements (Lemma 3.7), while a
nonzero polynomial of degree ≤ 2^{n−d−1} cannot have 2^{n−d} roots —
contradiction. "In particular": the windows {n−d+1, …, n−1, 0} and
{n−d, …, n−1} otherwise contain only coefficients vanishing by degree.
Part 2: suppose a_j = 0 on J = {r, …, r + d − 2}, r ≢ n − d + 1 (mod n).
Put k ≡ n − (r + d) + 1 (mod n), so k ∈ {1, …, n − 1} (k = 0 iff
r ≡ n − d + 1), and Q′(X) = (Q_L(X))^{2^k} mod (X^{2^n} − X). By
Claim 3.9 its X^{2^j}-coefficient is a_{j−k}^{2^k}; for
j ∈ [n − d + 1, n − 1], (j − k) mod n runs over J, so those coefficients
vanish: deg Q′ ≤ 2^{n−d}, and Q′ ≠ 0 as before. Choose c ∈ F with
deg(Q′ − c·Q_L) ≤ 2^{n−d−1}, matching leading coefficients (if
deg Q′ < 2^{n−d} take c = 0; the root count then gives Q′ ≡ 0, already
a contradiction). Q′ and Q_L both vanish on ker Q_L (2^{n−d} elements),
so Q′ = c·Q_L: y^{2^k} = cy for all y in the image L. Since d ≥ 1,
c ≠ 0, and Claim 3.8 gives L ⊆ b·F_{2^k}. For β ∈ L \ {0} and α ∈ L,
α/β ∈ F_{2^k} ∩ F = F_{2^{gcd(k,n)}} = F_2 (n prime), so
L ⊆ β·F_2, a constant multiple of a proper subfield — contradiction. □

**Remark (affine form used in §6).** §6 applies Lemma 3.10 to the
polynomials Q_i of *affine* spaces A_i, Q_i(X) = Σ_j a_{i,j} X^{2^j} +
â_i. The a_{i,j} are the coefficients of the linearized part of Q_i,
which by Lemma 3.7 is the image-subspace polynomial of the direction of
A_i — a linear space — so Lemma 3.10 applies to them: a_{i,j} = 0 for
j ∈ [n − d_i + 1, n − 1], a_{i,0} ≠ 0, a_{i,n−d_i} ≠ 0, and part 2's
hypothesis is consumed in the direction form (the direction lies in no
constant multiple of a proper subfield). This direction form is not
implied by Theorem 6.1's affine-shift hypothesis for arbitrary A_i; in
the §7 application n is prime and both forms are vacuous there.

**Definition 3.11 (associated matrix).** For Q(X) = Σ_{j=0}^{n-1} a_j
X^{2^j}, M_Q is the n × n matrix with (i, j)-entry, indices mod n,

    (M_Q)_{i,j} = (a_{j−i})^{2^i}

(row 0 is (a_0, …, a_{n−1}); each later row is the Frobenius square and
cyclic shift of the previous). The indicator M_Q′ has entry a′_{j−i},
where a′_m = 1 if a_m ≠ 0 and 0 otherwise. M[A, B] is the minor on rows
A, columns B (|A| = |B| ≠ 0); B − r = {s − r mod n : s ∈ B}.

**Claim 3.12.** (M_Q)_{i,j} is the coefficient of X^{2^j} in
(Q(X))^{2^i} mod (X^{2^n} − X).

*Proof.* Claim 3.9 with s = i. □

**Lemma 3.13 (permanent formula).** Let µ ≠ 0, A_i, B_i ⊆
{0, …, n − 1} with |A_i| = |B_i| > 0, α_i = Σ_{j∈A_i} 2^j,
β_i = Σ_{k∈B_i} 2^k, and Q_i linearized with associated matrices M_i.
The coefficient c_M of M = ∏_{i=1}^t X_i^{β_i} in

    R = Tr( µ · ∏_{i=1}^t (Q_i(X_i))^{α_i} ) mod ⟨X_i^{2^n} − X_i⟩_{i∈[t]}   (7)

is

    c_M = Σ_{r=0}^{n-1} µ^{2^r} · ∏_{i=1}^t Perm( M_i[A_i + r, B_i] )
        = Σ_{r=0}^{n-1} µ^{2^r} · ∏_{i=1}^t Perm( M_i[A_i, B_i − r] )^{2^r}.   (8)

*Proof.* Tr(Z) = Σ_r Z^{2^r} expands R; fix i and drop it from the
notation. ((Q^α))^{2^r} = ∏_{j∈A+r} (Q(X))^{2^j}; the monomial X^β
arises by choosing, for each k ∈ B, a factor index h(k) ∈ A + r and
taking X^{2^k} from (Q(X))^{2^{h(k)}}. Since |A| = |B| these choices
are exactly the bijections h: B → A + r, and by Claim 3.12 the chosen
terms contribute (M)_{h(k),k}; summing over h gives Perm(M[A + r, B]):
the first expression. For the second, (M)_{j+r,k} = ((M)_{j,k−r})^{2^r},
so M[A + r, B] is M[A, B − r] under entrywise 2^r-power and Perm
commutes with Frobenius. □

**Remark 3.14 (usage criterion).** In Lemma 3.13's notation, suppose
(1) each M_i′[A_i, B_i] is, up to reordering rows and columns,
triangular with all diagonal entries nonzero, and (2) for every
r ∈ {1, …, n − 1} some M_i′[A_i, B_i − r] has an all-zero column. Then
c_M ≠ 0: by (1) and (8) the r = 0 summand is µ·∏_i Perm(M_i[A_i, B_i])
≠ 0 (a triangular permanent is the product of its diagonal), and by (2)
every r ≠ 0 summand has a zero-column permanent factor, hence vanishes.
□

## 6. Disperser for independent affine sources

F = F_{2^n}; "affine" and "dimension" are over F_2. Support of Q_A is
governed by Lemma 3.10 via the affine-form remark; M_i, M_i′ are as in
Definition 3.11; trace coefficients by Lemma 3.13 and Remark 3.14.

**Theorem 6.1 (BK Theorem 6.1).** Let A_1, …, A_t ⊆ F be affine spaces
of dimensions d_1, …, d_t > 1, none contained in an affine shift of a
proper subfield of F. Let Q_i(X_i) be the image-subspace polynomial of
A_i, µ ∈ F \ {0}, and 1 ≤ e_i < d_i − 1. Put α_i = Σ_{j=0}^{e_i} 2^j. If

    Σ_{i=1}^t ( d_i − (e_i + 1) ) > n − max_{1≤i≤t} d_i + 1 ,        (†)

then R = Tr( µ · ∏_{i=1}^t (Q_i(X_i))^{α_i} ) mod
⟨ X_i^{2^n} − X_i ⟩_{i∈[t]} has a monomial ∏_{i=1}^t X_i^{β_i} with
wt_2(β_i) = e_i + 1 whose coefficient is nonzero. In particular
|R(F^t)| > 1, by Proposition 3.1. (Theorem 2.3 is the case
e_1 = … = e_t = 1, µ = 1: π(∏_{i=1}^t x_i^{1+2}) is nonconstant on
A_1 × ⋯ × A_t whenever Σ_i (d_i − 2) > n, via Propositions 3.1, 3.3.)

*Proof.* WLOG d_1 = max_i d_i. In Lemma 3.13 take row sets
𝒜_i = {0, …, e_i}. For B_i ⊆ {0, …, n−1} with |B_i| = e_i + 1, write
β_i = Σ_{k∈B_i} 2^k and M = ∏_{i=1}^t X_i^{β_i}. By (8), applied to the
linearized parts of the Q_i as justified next, c_M is

    c_M = Σ_{r=0}^{n−1} ( µ · ∏_{i=1}^t Perm( M_i[𝒜_i, B_i − r] ) )^{2^r} ,    (14)

with M_i the associated matrix of the linearized part of Q_i.

*Constant terms.* Lemma 3.13 needs linearized inputs, while
Q_i = Q_i^{lin} + â_i. But c_M is the same for both: expanding
(Q_i)^{α_i} = ∏_{j=0}^{e_i} (Q_i(X_i))^{2^j} with
(Q_i(X_i))^{2^j} = (Q_i^{lin}(X_i))^{2^j} + â_i^{2^j}, a monomial of
F_2-degree e_i + 1 in X_i arises only by choosing the linearized term
in all e_i + 1 factors — a constant choice contributes F_2-degree 0,
the remaining factors at most e_i in total, and reduction mod
X_i^{2^n} − X_i never increases F_2-degree. So Lemma 3.13 applies with
M_i as above. Per Remark 3.14, choose B_i so that each
M_i′[𝒜_i, B_i] is lower triangular with nonzero diagonal (r = 0 summand
≠ 0) while every r ∈ {1, …, n−1} has some M_i′[𝒜_i, B_i − r] with an
all-zero column: Claims 6.2 and 6.3. Support facts used throughout:
a_{i,0} ≠ 0, a_{i,n−d_i} ≠ 0, a_{i,j} = 0 for j ∈ [n − d_i + 1, n − 1].

**Choice of B_i.** Offsets: ℓ_1 = 0,
ℓ_{i+1} = (ℓ_i + d_i − (e_i + 1)) mod n, i.e.
ℓ_i = (Σ_{i′<i} (d_{i′} − (e_{i′} + 1))) mod n. Let t_0 be minimal with
S_{t_0+1} := Σ_{i′≤t_0} (d_{i′} − (e_{i′} + 1)) ≥ n − d_1 + 1; it exists
by (†), since S_{t+1} > n − d_1. For i ≤ t_0, minimality gives
S_i ≤ n − d_1 (S_1 = 0), so no wrap occurs: ℓ_i = S_i as a plain
integer and ℓ_i ≤ n − d_1 ≤ n − d_i. Fix i ≤ t_0.
Case a_{i,ℓ_i} ≠ 0: B_i = {ℓ_i} ∪ {n − d_i + 1, …, n − d_i + e_i},
which has 1 + e_i elements, ℓ_i distinct from the top part as
ℓ_i ≤ n − d_i.
Case a_{i,ℓ_i} = 0: then ℓ_i ≥ 1 (a_{i,0} ≠ 0) and ℓ_i < n − d_i
(a_{i,n−d_i} ≠ 0), so j_i = min{j > ℓ_i : a_{i,j} ≠ 0} and
j̄_i = max{j < ℓ_i : a_{i,j} ≠ 0} exist with j̄_i < ℓ_i < j_i ≤ n − d_i.
Put q_i = j_i − j̄_i − 1 (the zero-gap length around ℓ_i) and
s_i = min{q_i, e_i}; set
B_i = {j_i} ∪ {j̄_i + 1, …, j̄_i + s_i} ∪ {n − d_i + s_i + 1, …, n − d_i + e_i}.
The parts have 1, s_i, e_i − s_i elements, are pairwise disjoint
(j̄_i + s_i ≤ j̄_i + q_i = j_i − 1 < j_i ≤ n − d_i < n − d_i + s_i + 1),
and the last is empty if s_i = e_i.
For i > t_0: B_i = {0} ∪ {n − d_i + 1, …, n − d_i + e_i} (the first
case at ℓ = 0); only Claim 6.2 is needed from these i.

**Claim 6.2 (BK Claim 6.2).** Every M_i′[𝒜_i, B_i] is, after reordering
columns, lower triangular with nonzero diagonal.

*Proof.* The (u, k)-entry is a′_{i,(k−u) mod n}, rows u = 0, …, e_i;
since e_i ≤ d_i − 2 < n, indices below are plain integers. Fix i; write
a′_j, d, e for a′_{i,j}, d_i, e_i.
Case a′_{ℓ_i} = 1 (including i > t_0, with 0 for ℓ_i): order columns
(ℓ_i, n−d+1, …, n−d+e). Diagonal indices are ℓ_i (a′ = 1) and, for
column n−d+m at position m, (n−d+m) − m = n−d (a′_{n−d} = 1 by
Lemma 3.10(1)). Above-diagonal indices lie in {n−d+1, …, n−d+e} (row 0;
resp. n−d+(m−u) with 1 ≤ m−u ≤ e−1 for row u ≥ 1, column m > u), a
subset of [n−d+1, n−1] where a′ = 0 by the degree bound.
Case a′_{ℓ_i} = 0: order columns (j_i, j̄_i + 1, …, j̄_i + s_i,
n − d + s_i + 1, …, n − d + e). In the first block, the diagonal entries
are a′_{j_i} = 1 at (0,0) and a′_{j̄_i} = 1 at (m, m), 1 ≤ m ≤ s_i;
above-diagonal indices are j̄_i + 1, …, j̄_i + s_i (row 0) and
j̄_i + (m − u), 1 ≤ m − u ≤ s_i − 1 (row u ≥ 1, column m > u), all in
the zero-gap (j̄_i, j_i) since j̄_i + s_i ≤ j_i − 1. The remaining
columns have diagonal index n − d (a′ = 1) and above-diagonal indices
in {n−d+1, …, n−d+e}, where a′ = 0; entries below the first block are
below the diagonal. ∎

**Proposition 6.4 (BK Proposition 6.4).** If k ∈ B_i and
r ∈ [k + 1, k + d_i − (e_i + 1)] mod n, then M_i′[𝒜_i, B_i − r] has an
all-zero column.

*Proof.* The k-column is (a′_{i,k−r}, …, a′_{i,k−(r+e_i)})^T, indices
mod n; it is all-zero if

    [k − (r + e_i), k − r] mod n ⊆ [n − d_i + 1, n − 1] .    (15)

The range of r gives k − r ∈ [1 − d_i + e_i, −1] for integer
representatives, so [k − (r + e_i), k − r] ⊆ [1 − d_i, −1], which mod n
is exactly [n − d_i + 1, n − 1]. ∎

**Claim 6.3 (BK Claim 6.3).** For every r ∈ {1, …, n−1} some
M_i′[𝒜_i, B_i − r] has an all-zero column.

*Proof.* By Proposition 6.4 it suffices to show

    [1, n−1] ⊆ ∪_{i=1}^t ∪_{k∈B_i} [k + 1, k + d_i − (e_i + 1)] mod n .    (16)

Intervals are mod n; for i ≤ t_0, ℓ_i = S_i are plain integers and
[ℓ_i + 1, ℓ_{i+1}] denotes the cyclic interval ending at
ℓ_{i+1} = S_{i+1} mod n. For i = 1: ℓ_1 = 0, a_{1,0} ≠ 0,
B_1 = {0} ∪ [n − d_1 + 1, n − d_1 + e_1]; k = 0 gives
[1, d_1 − (e_1 + 1)] = [ℓ_1 + 1, ℓ_2], and k over the top part chains
(consecutive starts, common length d_1 − e_1 − 1 ≥ 1) into
[n − d_1 + 2, n − 1]. Hence a zero column for every

    r ∈ [1, d_1 − (e_1 + 1)] ∪ [n − d_1 + 2, n − 1]
      = [ℓ_1 + 1, ℓ_2] ∪ [n − d_1 + 2, n − 1] .    (17)

We show for 1 < i ≤ t_0 that

    ∪_{k∈B_i} [k + 1, k + d_i − (e_i + 1)] ⊇ [ℓ_i + 1, ℓ_{i+1}] .    (18)

With (17) this gives (16): if S_{t_0+1} < n the chained intervals reach
ℓ_{t_0+1} = S_{t_0+1} ≥ n − d_1 + 1 and (17) covers the rest; if
S_{t_0+1} ≥ n they wrap and cover [1, n−1] alone. Three cases for (18).
Case a_{i,ℓ_i} ≠ 0: ℓ_i ∈ B_i and k = ℓ_i gives
[ℓ_i + 1, ℓ_i + d_i − (e_i + 1)] = [ℓ_i + 1, ℓ_{i+1}].
Case a_{i,ℓ_i} = 0, q_i < e_i: s_i = q_i and
B_i = [j̄_i + 1, j_i] ∪ {n − d_i + q_i + 1, …, n − d_i + e_i}, the first
two parts consecutive since j̄_i + q_i + 1 = j_i; then
ℓ_i ∈ [j̄_i + 1, j_i − 1] ⊆ B_i and the previous case applies.
Case a_{i,ℓ_i} = 0, q_i ≥ e_i: s_i = e_i, the top part is empty,
B_i = {j_i} ∪ {j̄_i + 1, …, j̄_i + e_i}. *This step uses the subfield
hypothesis.* Since j̄_i < ℓ_i ≤ n − d_1 ≤ n − d_i, we have
j̄_i + 1 ≤ n − d_i, in particular j̄_i + 1 ≠ n − d_i + 1. Apply
Lemma 3.10(2) to Q_i's coefficients — i.e. to the direction of A_i; the
hypothesis consumed is the direction form. In the §7 application the
distinction is vacuous: n is prime, the only proper subfield is F_2,
and |Y_i| = 2^{d_i} ≥ 2^{r+2}, while an affine shift of F_2 and a
constant multiple of F_2 each have 2 elements. The zero run
j̄_i + 1, …, j_i − 1, of length q_i, does not start at n − d_i + 1, so
it contains no d_i − 1 consecutive indices: q_i ≤ d_i − 2, i.e.
j_i ≤ j̄_i + d_i − 1. The intervals for k = j̄_i + 1, …, j̄_i + e_i chain
into [j̄_i + 2, j̄_i + d_i − 1] ∋ j_i, with j̄_i + 2 ≤ ℓ_i + 1; adjoining
k = j_i extends the union gaplessly to j_i + d_i − (e_i + 1) ≥
ℓ_i + d_i − (e_i + 1). The union therefore contains
[ℓ_i + 1, ℓ_{i+1}]. ∎

*Conclusion of Theorem 6.1.* Claims 6.2 and 6.3 are Remark 3.14's two
hypotheses (lower for upper triangular is immaterial, reversing rows
and columns): the r = 0 summand of (14) is a product of permanents of
triangular matrices with nonzero diagonal, hence nonzero, and every
r ≠ 0 summand contains a vanishing permanent factor. So c_M ≠ 0 for
M = ∏_{i=1}^t X_i^{β_i} with wt_2(β_i) = |B_i| = e_i + 1, and
|R(F^t)| > 1 by Proposition 3.1. ∎

## 7. Affine spaces of sublinear dimension: Theorem 7.3 and Theorem 2.2

Now n is prime, F = F_{2^n}, r blocks, t < r. Used from §3:
Propositions 3.1 (reduction preserves values; two values ⇒ a
positive-degree monomial) and 3.3 (π = Tr(µ·)), Lemmas 3.4 (linear maps
are linearized polynomials), 3.6 and 3.7; from §6, Theorem 6.1.
Reduction mod ⟨X_i^{2^n} − X_i⟩ does not increase F_2-degrees, and
linearized polynomials stay sums of monomials a·X_i^{2^j} (F_2-degree
≤ 1 in their variable) under reduction and Frobenius powers.

**Lemma 7.1 (Bourgain's decomposition).** Let A ⊆ (F_{2^n})^r be an
F_2-affine subspace and γ ∈ A. Then there exist linear subspaces
Y_1, …, Y_r ⊆ F_{2^n} and F_2-linear maps σ_{ij}: Y_j → F_{2^n} (j < i)
such that

    A = { (x_1,…,x_r) : x_i = γ_i + y_i + Σ_{j<i} σ_{ij}(y_j),  y_i ∈ Y_i },

and dim A = Σ_{i∈[r]} dim Y_i.

*Proof.* Write A = γ + V, V linear. Let U_i = {v ∈ V : v_1 = … =
v_{i−1} = 0} (U_1 = V) and Y_i = {v_i : v ∈ U_i}; the coordinate map
ε_i: U_i → Y_i is linear and surjective. Choose linear sections
s_i: Y_i → U_i with ε_i ∘ s_i = id, and set σ_{ij}(y) = (s_j(y))_i for
j < i. Every v ∈ V is uniquely v = Σ_i s_i(y_i): put y_1 = v_1, so
v − s_1(y_1) ∈ U_2; put y_2 = (v − s_1(y_1))_2, so the remainder lies in
U_3; continue. Uniqueness follows by reading coordinates in order (the
first coordinate of Σ_i s_i(y_i) is y_1, since s_i(y_i) ∈ U_i has first
coordinate 0 for i ≥ 2). Hence ⊕_i Y_i ≅ V and dim V = Σ_i dim Y_i. The
i-th coordinate of Σ_j s_j(y_j) is y_i + Σ_{j<i} σ_{ij}(y_j):
(s_i(y_i))_i = y_i, (s_j(y_j))_i = σ_{ij}(y_j) for j < i, and
(s_j(y_j))_i = 0 for j > i since s_j(y_j) ∈ U_j vanishes in coordinates
< j. Adding γ gives A. □ (This is the row-echelon form of a basis of V.)

**Lemma 7.2 (coefficient invariance under couplings and shifts).** For
i ∈ [r] let P_i(X_i) ∈ F_{2^n}[X_i] be linearized; for j < i let
P_{ij}(X_j) ∈ F_{2^n}[X_j] be linearized; let γ ∈ (F_{2^n})^r,
I_0 = {i_1 < … < i_t} ⊆ [r], e_i ≥ 1 for i ∈ I_0,
α_i = Σ_{k=0}^{e_i} 2^k, and µ ∈ F_{2^n} \ {0}. Let

    g  = Tr( µ · Π_{i∈I_0} ( P_i(X_i) + Σ_{j<i} P_{ij}(X_j) + γ_i )^{α_i} )
             mod ⟨X_i^{2^n} − X_i⟩,
    g′ = Tr( µ · Π_{i∈I_0} P_i(X_i)^{α_i} )  mod ⟨X_i^{2^n} − X_i⟩.

Then for every (β_i)_{i∈I_0} with wt_2(β_i) = e_i + 1, the coefficients
of M = Π_{i∈I_0} X_i^{β_i} in g and in g′ are equal.

*Proof.* Put P̂_i = Σ_{j<i} P_{ij}(X_j) + γ_i, which does not involve
X_i. Expanding the trace, Tr(z) = Σ_{s=0}^{n−1} z^{2^s}, and using
(z^{α_i})^{2^s} = Π_{l=0}^{e_i} z^{2^{s+l}} (since α_i = Σ_{l=0}^{e_i}
2^l),

    g = Σ_{s=0}^{n−1} µ^{2^s} · Π_{i∈I_0} Π_{l=0}^{e_i}
            ( P_i(X_i) + P̂_i )^{2^{s+l}}  mod ⟨X_i^{2^n} − X_i⟩,   (19)–(20)

by the Freshman's dream in characteristic 2 each factor splits as
P_i(X_i)^{2^{s+l}} + P̂_i^{2^{s+l}}. Choosing the P̂-term exactly at the
positions of L ⊆ {0, …, e_i} gives, per i ∈ I_0,

    Π_{l=0}^{e_i} ( P_i(X_i) + P̂_i )^{2^{s+l}}
        = ( P_i(X_i) )^{Σ_{l=0}^{e_i} 2^{s+l}} + G_i(X_{≤i}),   (21)–(23)

where every L ≠ ∅ term has F_2-degree at most e_i + 1 − |L| ≤ e_i in
X_i (only its e_i + 1 − |L| factors P_i(X_i)^{2^{s+l}} involve X_i,
each of F_2-degree ≤ 1, while P̂-terms do not involve X_i at all); so
G_i has F_2-degree ≤ e_i in X_i and involves only X_1, …, X_i. Now fix
M and the summand s, and peel off the largest index i_t: the product of
the remaining brackets, H_t, does not mention X_{i_t}. M has F_2-degree
e_{i_t} + 1 in X_{i_t} while G_{i_t}·H_t has at most e_{i_t} there, so
only the pure term (P_{i_t}(X_{i_t}))^{Σ_l 2^{s+l}} contributes to M's
coefficient. Peeling at i_{t−1} is identical: no other factor mentions
X_{i_{t−1}} (H_{t−1} uses only X_{≤ i_{t−2}}, earlier brackets only
X_{≤i}, the surviving pure term only X_{i_t}), and G_{i_{t−1}} has
F_2-degree ≤ e_{i_{t−1}} there. Continuing down to i_1, the
coefficient of M in the s-th summand equals its coefficient in
(µ · Π_{i∈I_0} P_i(X_i)^{α_i})^{2^s}. Summing over s gives the
coefficient in g′. □

**Theorem 7.3 (affine disperser, non-parameterized).** Let t < r and
let n be prime with n ≥ r(r+1)/t. For i ∈ [r] let e_i = r + 1 − i and
α_i = Σ_{k=0}^{e_i} 2^k. Let µ ∈ F_{2^n} \ {0} and

    f(x_1, …, x_r) = Tr( µ · Σ_{I⊆[r], |I|=t} Π_{i∈I} x_i^{α_i} ).

Then every affine A ⊆ (F_{2^n})^r with dim A > nr/t + nt + r(r+1)
satisfies |f(A)| > 1.

*Proof.* (1) *Parameterization.* Fix γ ∈ A and decompose A by
Lemma 7.1, with linear Y_i and σ_{ij}. Let Q_i ∈ F_{2^n}[X_i] be the
image-subspace polynomial of Y_i (Lemma 3.7); Q_i is linearized since
Y_i is linear. For j < i, extend σ_{ij} arbitrarily to an F_2-linear
map on F_{2^n} and let Q_{ij} ∈ F_{2^n}[X_j] be the linearized
polynomial of x ↦ σ_{ij}(Q_j(x)) (Lemma 3.4); Q_{ij}(x) =
σ_{ij}(Q_j(x)) for all x, independently of the extension, since
Q_j(x) ∈ Y_j. Set R_i = Q_i(X_i) + Σ_{j<i} Q_{ij}(X_j) + γ_i and
H = (R_1, …, R_r). As X_j ranges over F_{2^n}, Q_j(X_j) ranges over
Y_j (Lemma 3.7) and the couplings track the same y_j, so H(F^r) = A.
With h = f ∘ H reduced,

    h = Σ_{I⊆[r], |I|=t} Tr( µ · Π_{i∈I} R_i^{α_i} )
            mod ⟨X_i^{2^n} − X_i⟩.                                   (24)–(25)

By Proposition 3.1, h(F^r) = f(A); find a positive-degree monomial
with nonzero coefficient in h.
(2) *High-entropy blocks.* Let d_i = dim Y_i, so Σ_i d_i = dim A >
nr/t + nt + r(r+1), and S = {i : d_i > r+1}. Then |S| ≥ t: otherwise
Σ_i d_i ≤ |S|·n + r(r+1) ≤ (t−1)n + r(r+1) < nt + r(r+1) + nr/t (each
d_i ≤ n; each of the ≤ r indices outside S has d_i ≤ r+1). Also
Σ_{i∈S} d_i ≥ Σ_{i∈[r]} (d_i − (r+1)) = dim A − r(r+1) > nr/t + nt
(terms dropped for i ∉ S are ≤ 0). Averaging over the t-subsets I_0 of
S, some I_0 satisfies

    Σ_{i∈I_0} (d_i − (r+1)) ≥ (t/|S|) · Σ_{i∈S} d_i − (r+1)t
        ≥ (t/r)(nr/t + nt) − (r+1)t = n + nt²/r − (r+1)t ≥ n,        (26)

the last step by n ≥ r(r+1)/t, i.e. nt/r ≥ r+1.
(3) *The I_0-summand.* Let g be the I_0-term of (24)–(25) and
g′ = Tr(µ · Π_{i∈I_0} Q_i(X_i)^{α_i}) mod ⟨X_i^{2^n} − X_i⟩. Apply
Theorem 6.1 to the spaces Y_i, i ∈ I_0, with these Q_i, this µ, and
e_i = r + 1 − i. Its hypotheses hold: d_i > r + 1 ≥ 3, in particular
d_i > 1; 1 ≤ e_i < d_i − 1, since e_i = r + 1 − i ≥ 1 (i ≤ r) and
e_i + 1 ≤ r + 1 < d_i; Σ_{i∈I_0} (d_i − (e_i + 1)) ≥
Σ_{i∈I_0} (d_i − (r+1)) ≥ n by (26) (using e_i + 1 ≤ r + 1), and
n > n − max_{i∈I_0} d_i + 1 since max d_i > r + 1 ≥ 3; and no Y_i is
contained in an affine shift of a proper subfield — **where primality
of n is used**: the subfields of F_{2^n} are F_{2^d} for d | n, so n
prime leaves only F_2, an affine shift of which has exactly 2 elements,
while |Y_i| = 2^{d_i} ≥ 2^{r+2}. Theorem 6.1 gives a monomial
M = Π_{i∈I_0} X_i^{β_i} with wt_2(β_i) = e_i + 1 and nonzero
coefficient in g′. By Lemma 7.2 (P_i = Q_i, P_{ij} = Q_{ij}, shifts
γ_i), the coefficient of M in g equals that in g′; in particular it is
nonzero.
(4) *Uniqueness.* Each R_j is a sum of monomials a·X_i^{2^b} with
i ≤ j, plus possibly a constant: F_2-degree ≤ 1 in variables X_{≤j}.
Let J ⊆ [r], |J| = t, J ≠ I_0, and suppose M had nonzero coefficient
in the J-summand; expanding the trace, M must occur in
Π_{j∈J} Π_{l=0}^{e_j} R_j^{2^{l+w}} for some w ∈ {0, …, n−1}. Choose one
monomial from each factor (j, l) whose product is M; set λ_j(l) = i if
the monomial is in X_i (so i ≤ j), and λ_j(l) := j if the constant
term was chosen — this only enlarges the counts |λ_j^{−1}(i)| for
i ∈ I_0 below, preserves λ_j(l) ≤ j, and keeps
Σ_{i∈I_0} |λ_j^{−1}(i)| ≤ r + 2 − j, so Claim 7.4 applies to the
extended functions. The F_2-degree of M in X_i is at most the number of
(j, l) with λ_j(l) = i (binary carries only decrease digit sums), and
for i ∈ I_0 it equals e_i + 1 = r + 2 − i. Hence

    Σ_{j∈J} |λ_j^{−1}(i)| ≥ r + 2 − i      for every i ∈ I_0.        (27)

**Claim 7.4.** Let I_0 ⊆ [r], |I_0| = t. If J ⊆ [r], |J| = t, and
λ_j: {0, …, r+1−j} → {1, …, j} for j ∈ J satisfy (27), then J = I_0.

*Proof.* For each j ∈ J, Σ_{i∈I_0} |λ_j^{−1}(i)| ≤ |{0, …, r+1−j}| =
r + 2 − j. For general k ≥ 0, Σ_{j∈J} Σ_l (j^k − λ_j(l)^k) ≥ 0 since
λ_j(l) ≤ j; expanding it and using (27) to lower-bound
Σ_j Σ_l λ_j(l)^k = Σ_i i^k·|{(j,l) : λ_j(l) = i}| ≥
Σ_{i∈I_0} (r+2−i)·i^k yields
Σ_{i∈I_0} i^{k+1} − Σ_{j∈J} j^{k+1} ≥ (r+2)·( Σ_{i∈I_0} i^k −
Σ_{j∈J} j^k ), which is ≥ 0 by induction on k (k = 0: equal
cardinalities). Hence for all k,

    Σ_{i∈I_0} i^k ≥ Σ_{j∈J} j^k.                                     (28)

Let i_1 = max I_0, j_1 = max J. By (28) with k → ∞, i_1 ≥ j_1 (if
j_1 > i_1 then Σ_{j∈J} j^k ≥ j_1^k > t·i_1^k ≥ Σ_{i∈I_0} i^k for large
k). But i_1 ≤ j_1: by (27) some λ_j(l) equals i_1 (the count is
≥ r + 2 − i_1 ≥ 1), and λ_j(l) ≤ j ≤ j_1. So i_1 = j_1, and (27) at
i = i_1, j = j_1 forces |λ_{j_1}^{−1}(i_1)| ≥ r + 2 − j_1, the full
number of l's: λ_{j_1}(l) = i_1 for all l. Deleting i_1 from I_0 and
j_1 from J leaves (28) valid for the deleted sets, so the second-largest
elements satisfy i_2 ≥ j_2 by the same domination; and i_2 ≤ j_2 since
no remaining j can map to i_2 except one ≥ i_2, while j_1 is exhausted
by i_1 ≠ i_2. Thus i_2 = j_2; inducting on the s-th largest elements
gives I_0 = J. □ (Claim 7.4)

So M, of positive degree since wt_2(β_i) = e_i + 1 ≥ 2, appears in h
with the nonzero coefficient from step (3), contributed by the
I_0-summand alone. By Proposition 3.1, |h(F^r)| > 1, hence
|f(A)| > 1. □ (Theorem 7.3)

**Deduction of Theorem 2.2 from Theorem 7.3.** Fix the input length m
and set n = least prime > 2m^{3/5}, r = ⌈m/n⌉, t = ⌈√r⌉. By Bertrand's
postulate there is a prime in (2m^{3/5}, 4m^{3/5}], so

    2m^{3/5} < n ≤ 4m^{3/5}.

Hence m^{2/5}/4 ≤ m/n ≤ r ≤ m/n + 1 < m^{2/5}/2 + 1; in particular
r → ∞ with m, and t = ⌈√r⌉ < r once r ≥ 3.
*Block-length hypothesis.* Since t ≥ √r,

    r(r+1)/t ≤ (r+1)·√r = (r+1)^{3/2} ≤ (m^{2/5}/2 + 2)^{3/2}
        = (1/(2√2)) · m^{3/5} · (1 + 4m^{−2/5})^{3/2}
        = (1/(2√2)) · m^{3/5} (1 + o(1)) < 2m^{3/5} < n

for all large m (the leading constant 1/(2√2) ≈ 0.354 is well below 2).
So n ≥ r(r+1)/t holds.
*Threshold.* Estimate nr/t + nt + r(r+1) *jointly*, using r ≤ m/n + 1
rather than maximizing n and r independently. With t ≥ √r,

    nr/t ≤ n√r = (n²·(m/n + 1))^{1/2} = (nm + n²)^{1/2}
        ≤ (4m^{8/5} + 16m^{6/5})^{1/2} = 2m^{4/5}(1 + o(1)),
    nt ≤ n(√r + 1) ≤ 2m^{4/5}(1 + o(1)) + 4m^{3/5}
        = 2m^{4/5}(1 + o(1)),
    r(r+1) ≤ (m^{2/5}/2 + 1)(m^{2/5}/2 + 2) = (1/4)m^{4/5}(1 + o(1)).

Adding, nr/t + nt + r(r+1) ≤ (2 + 2 + 1/4)·m^{4/5}(1 + o(1))
= 4.25·m^{4/5}(1 + o(1)) < 6m^{4/5} for all large m. Thus every affine
A ⊆ F_2^m with dim A > 6m^{4/5} satisfies the dimension hypothesis of
Theorem 7.3, after the identification below.
*Identification of the two f's.* By Proposition 3.3, the nontrivial
F_2-linear π of Theorem 2.2 is π(z) = Tr(µz) for a unique µ ≠ 0,
matching Theorem 7.3. The exponents differ only in indexing order:
equation (4) of Theorem 2.2 gives block i the exponent
Σ_{j=0}^{i} 2^j, while Theorem 7.3 gives block i the exponent
Σ_{k=0}^{e_i} 2^k with e_i = r + 1 − i, i.e. Σ_{j=0}^{r+1−i} 2^j. Under
the renaming of blocks i ↦ r + 1 − i, the two symmetric sums coincide:
f of Theorem 7.3 at (x_1, …, x_r) equals f of Theorem 2.2 at
(x_r, …, x_1). Coordinate permutation preserves affinity and dimension,
so Theorem 7.3's conclusion transfers: f of Theorem 2.2 is nonconstant
on every affine subspace of dimension > 6m^{4/5}, for all sufficiently
large m. ∎ (Theorem 2.2)

*Application to the family of §1 (the Main Lemma).* Rename (m, n) →
(N, k): k is the least prime > 2N^{3/5}, r = ⌈N/k⌉, t = ⌈√r⌉, and
φ: F_2^N → K^r = (F_{2^k})^r is the pad-and-split map of `proof.md` §1
(append rk − N zeros, split into r blocks of k bits, read each block in
the polynomial basis of K). φ is injective and F_2-linear, and it is
admissible in Theorem 2.2 at length m := N: the theorem is quantified
over all such φ, so the short last block (rk − N < k padded zeros)
needs no divisibility assumption, and all parameters k, r, t are
computed from the true length N. The exponents agree
(Σ_{j=0}^{i} 2^j = 2^{i+1} − 1) and π, the constant coefficient, is a
nontrivial F_2-linear functional. Hence for every affine A ⊆ F_2^N with
dim A > 6N^{4/5}, φ(A) ⊆ (F_{2^k})^r is affine of the *same* dimension
(φ injective linear), and |f_N(A)| = |f(φ(A))| > 1 by Theorem 2.2. The
estimates are asymptotic, so this holds for all sufficiently large N;
for N ≤ 7776 = 6^5 one has 6N^{4/5} ≥ N, no affine subspace of F_2^N
has dimension > 6N^{4/5}, and the assertion is vacuous. Hence f_N is an
affine disperser for dimension d(N) = ⌊6N^{4/5}⌋ + 1 at every length
where the property has content, in particular for all sufficiently
large N — which is all that `proof.md` §2 uses. □
