# Disperser: subspace polynomials, independent sources, Theorem 2.2 (BK §§3, 6, 7)

BK §§3, 6, 7; p = 2, prime block length; BK numbering kept.

**Main Lemma (consumed by `proof.md` §2).** N >= 1; k the least prime
> 2N^{3/5}; r = ⌈N/k⌉; t = ⌈√r⌉; f_N the function of `proof.md` §1
over K = F_{2^k}. Then f_N is nonconstant on every affine
A ⊆ F_2^N with dim A > 6N^{4/5}: an affine disperser for
d(N) = ⌊6N^{4/5}⌋ + 1. Vacuous for N <= 7776 = 6^5, as 6N^{4/5} >= N.
This is Theorem 2.2, proved at the end of §7.

## 3. Subspace polynomials

**Conventions.** p = 2; F = F_{2^n}, n prime (n = k in the
application); F̄ an algebraic closure; indices mod n.
Q(X) = Σ_j a_j X^{2^j} is *linearized*; with constant term â,
*affine linearized*: the F_2-degree <= 1 maps F → F.
P_L (3.5), Q_L (3.7): kernel-/image-subspace polys. of affine
L ⊆ F, dim d over F_2. wt_2(m) := binary digit sum.

**Proposition 3.1.** h ∈ F[X_1, …, X_r], h' := h mod
⟨X_i^{2^n} - X_i⟩ (unique representative of degree < 2^n in each
variable). Then h = h' on F^r, and |h(F^r)| > 1 iff h' has
positive degree. *Proof.* a^{2^n} = a on F: values agree. h' - c != 0 has
degree < 2^n per variable, so cannot vanish on F^r (root
bound; induct on r). □

**Proposition 3.3.** Tr(Y) := Σ_{i=0}^{n-1} Y^{2^i}: F → F_2.
Every F_2-linear π: F → F_2 is π(x) = Tr(µx) for a unique µ ∈ F;
π trivial iff µ = 0. *Proof.* µ -> Tr(µ·) is F_2-linear, and Tr(µ·) != 0 for
µ != 0 (degree 2^{n-1} < 2^n, root bound): injective
between n-dimensional spaces, hence onto. □

**Lemma 3.4.** φ: F_2^n → F an F_2-linear isomorphism => affine
maps T: F_2^n → F_2^n correspond bijectively to affine linearized
polynomials P_T with P_T(φ(b)) = φ(T(b)). *Proof.* Transport by φ. Both sides have F_2-dimension
n(n + 1): (a_0, …, a_{n-1}, â) vs. matrix + translation.
Evaluation is injective: ≡ 0 on F forces â = 0 (at 0) and
all a_j = 0 (degree <= 2^{n-1}, root bound). Hence bijective. □

**Definition 3.5 (kernel-subspace polynomial).** For affine L ⊆ F of
dimension d: P_L(X) = ∏_{α∈L} (X - α).

**Lemma 3.6.** P_L is monic affine linearized of degree 2^d,
and linearized iff L is linear. *Proof.* Induction on d, V linear
first. P_{{0}} = X. W = V ⊕ ⟨b⟩, b ∉ V: linearizedness gives
P_V(X + b) = P_V(X) + P_V(b), so P_V(X)^2 + P_V(b)P_V(X) is monic
linearized of degree 2^{d+1} vanishing exactly on V ∪ (V + b) = W
(2^{d+1} roots), hence = P_W; in particular ker P_V = V. Affine
L = γ + V: P_V(X) + P_V(γ) is monic affine linearized of degree
2^d with root set {x: P_V(x) = P_V(γ)} = γ + V = L, hence = P_L,
linearized part P_V. Linearized iff P_V(γ) = 0 iff γ ∈ ker P_V =
V iff 0 ∈ L iff L linear. □

**Lemma 3.7 (image-subspace polynomial).** L ⊆ F affine, dim d =>
there is a monic affine linearized Q_L of degree 2^{n-d} (the
image-subspace polynomial of L) with Q_L(F) = L and
P_L(Q_L(X)) = X^{2^n} - X. (6a)
L linear => Q_L linearized, ker Q_L = P_L(F), Q_L has 2^{n-d}
roots, and Q_L(P_L(X)) = X^{2^n} - X. (6b)
*Proof.* L linear: L_0 := P_L(F) is linear of dim n - d (P_L
linearized, ker = L). Q_L := P_{L_0} (3.6): monic linearized,
deg 2^{n-d}; Q_L ∘ P_L monic deg 2^n, vanishing on F
(P_L(x) ∈ L_0 = root set of Q_L), hence = Π_{α∈F}(X - α) =
X^{2^n} - X. g(Y) := P_L(Q_L(Y)) - (Y^{2^n} - Y) satisfies
g(P_L(X)) = P_L(X^{2^n} - X) - (P_L(X)^{2^n} - P_L(X)) = 0
(P_L linearized; P_L(X^{2^n}) = P_L(X)^{2^n}, coeffs in F);
deg(g ∘ P_L) = (deg g)2^d forces g ≡ 0: (6a). So im Q_L ⊆ ker
P_L = L and |im Q_L| = 2^d = |L| (ker Q_L ⊇ L_0, dims equal):
Q_L(F) = L, ker Q_L = L_0, and (6b) is the identity proved.
Affine L = γ + V: Q_L := Q_V + γ - monic affine linearized,
deg 2^{n-d}, Q_L(F) = V + γ = L, linearized part Q_V. 3.6 gives
P_L = P_V + P_V(γ); P_V linearized gives P_L(Q_L(Y)) =
P_V(Q_V(Y) + γ) + P_V(γ) = P_V(Q_V(Y)) = Y^{2^n} - Y: (6a).
(6b) holds only up to γ: Q_L(P_L(X)) = X^{2^n} - X + γ
(Q_V(P_V(γ)) = γ^{2^n} - γ = 0). Applications use (6a),
Q_L(F) = L, and linearizedness only for linear L. □

**Claim 3.8.** a^{2^k} = ca with a, c ∈ F, k >= 1, and b ∈ F̄ a
(2^k - 1)-th root of c => a ∈ bF_{2^k}. *Proof.* a != 0: (a/b)^{2^k-1} = a^{2^k-1}/c = 1 =>
(a/b)^{2^k} = a/b ∈ F_{2^k}. □

**Claim 3.9.** For Q(X) = Σ_{j=0}^{n-1} a_j X^{2^j} + â and any
integer s: (Q(X))^{2^s} mod (X^{2^n} - X) = Σ_j a_{j-s}^{2^s}
X^{2^j} + â^{2^s}, indices mod n. *Proof.* Frobenius-expand, use X^{2^n} ≡ X, reindex. □

**Lemma 3.10.** L ⊆ F linear, dim d, Q_L(X) = Σ_j a_j X^{2^j}
its image-subspace polynomial. (1) Every d consecutive indices
J = {(r + j) mod n: j = 0, …, d - 1} contain j with a_j != 0;
in particular a_0 != 0, a_{n-d} != 0. (2) If L ⊄ βF_{2^k} for
all β ∈ F, k < n, then every d - 1 consecutive indices
J = {(r + j) mod n: j = 0, …, d - 2} with r != n - d + 1 contain
j with a_j != 0. *Proof.* d = n: Q_L = X. d < n: a_j = 0 on
[n - d + 1, n - 1] (degree), a_{n-d} = 1 (monic).
Q' := (Q_L)^{2^s} mod (X^{2^n} - X) vanishes on ker Q_L
(Q'(α) = Q_L(α)^{2^s}), 2^{n-d} points (3.7), and its top coeffs
are window values of the a_j (3.9). (1) s ≡ -(r + d): coeffs at
j ∈ [n - d, n - 1] are a_{j+r+d}^{2^s}, indices over J; J zero =>
deg Q' <= 2^{n-d-1} < 2^{n-d} roots: contradiction (windows
{n-d+1, …, n-1, 0}, {n-d, …, n-1} give a_0, a_{n-d}).
(2) k := n - (r + d) + 1 mod n ∈ {1, …, n - 1} (k = 0 iff
r ≡ n - d + 1), s := k: J zero => deg Q' <= 2^{n-d}; matching
leading coeffs (deg(Q' - cQ_L) <= 2^{n-d-1}, both vanishing on
ker Q_L) gives Q' = cQ_L: y^{2^k} = cy on L. c != 0; 3.8:
L ⊆ bF_{2^k}; α/β ∈ F_{2^k} ∩ F = F_2 (n prime) =>
L ⊆ βF_2: contradiction. □

**Remark (affine form used in §6).** By 3.7, the linearized
part of A_i's Q_i is the image poly. of A_i's direction
(linear); §6 applies 3.10 to it, part (2) in direction form
(direction in no constant multiple of a proper subfield) —
not implied by 6.1's affine-shift hypothesis in general;
vacuous in §7 (n prime).

**Definition 3.11 (associated matrix).** For Q(X) = Σ_j a_j
X^{2^j}, M_Q is the n × n matrix with (M_Q)_{i,j} = (a_{j-i})^{2^i},
indices mod n (row 0 is (a_0, …, a_{n-1}); each later row is the
Frobenius square and cyclic shift of the previous). M_Q' has
entries a'_{j-i}, a'_m = 1 iff a_m != 0. M[A, B] := the minor on
rows A, columns B (|A| = |B| != 0); B - r := {s - r mod n: s ∈ B}.

**Claim 3.12.** (M_Q)_{i,j} is the coefficient of X^{2^j} in
(Q(X))^{2^i} mod (X^{2^n} - X). *Proof.* 3.9 with s = i. □

**Lemma 3.13 (permanent formula).** µ != 0, A_i, B_i ⊆
{0, …, n - 1} with |A_i| = |B_i| > 0, α_i = Σ_{j∈A_i} 2^j,
β_i = Σ_{k∈B_i} 2^k, Q_i linearized with matrices M_i. The
coefficient c_M of M = ∏_{i=1}^t X_i^{β_i} in
R = Tr( µ ∏_{i=1}^t (Q_i(X_i))^{α_i} ) mod
⟨X_i^{2^n} - X_i⟩_{i∈[t]} (7)
is c_M = Σ_{r=0}^{n-1} µ^{2^r} ∏_{i=1}^t Perm( M_i[A_i + r, B_i] )
= Σ_{r=0}^{n-1} µ^{2^r} ∏_{i=1}^t Perm( M_i[A_i, B_i - r] )^{2^r}. (8)
*Proof.* Expand Tr(Z) = Σ_r Z^{2^r}; fix i. X^β arises by
choosing h(k) ∈ A + r for each k ∈ B and the X^{2^k}-term of
(Q)^{2^{h(k)}}: coefficient (M)_{h(k),k} (3.12); |A| = |B| => the
choices are the bijections B → A + r, summing to
Perm(M[A + r, B]). Second form: (M)_{j+r,k} = ((M)_{j,k-r})^{2^r}
and Perm commutes with Frobenius. □

**Remark 3.14 (usage criterion).** In 3.13's notation, if (1)
each M_i'[A_i, B_i] is, up to row/column reordering, triangular
with nonzero diagonal, and (2) for every r ∈ {1, …, n - 1} some
M_i'[A_i, B_i - r] has an all-zero column, then c_M != 0: by (1)
the r = 0 summand of (8) is µ∏_i Perm(M_i[A_i, B_i]) != 0
(triangular permanent = diagonal product); by (2) every r != 0
summand has a zero-column permanent factor, hence is 0. □
## 6. Disperser for independent affine sources

**Theorem 6.1 (BK Theorem 6.1).** Let A_1, …, A_t ⊆ F be affine of
dimensions d_1, …, d_t > 1, none in an affine shift of a proper
subfield of F. Q_i(X_i) := image-subspace poly. of A_i,
µ ∈ F \ {0}, 1 <= e_i < d_i - 1, α_i = Σ_{j=0}^{e_i} 2^j. If
Σ_{i=1}^t ( d_i - (e_i + 1) ) > n - max_{1<=i<=t} d_i + 1, (†)
then R (as in (7), for these Q_i and α_i) has a monomial
∏_i X_i^{β_i} with
wt_2(β_i) = e_i + 1 and nonzero coefficient; hence |R(F^t)| > 1
by 3.1.
*Proof.* Wlog d_1 = max d_i. Apply (8) with row sets
𝒜_i = {0, …, e_i} to the linearized parts (next paragraph):
for |B_i| = e_i + 1, β_i = Σ_{k∈B_i} 2^k,
c_M = Σ_{r=0}^{n-1} ( µ ∏_{i=1}^t Perm( M_i[𝒜_i, B_i - r] ) )^{2^r}, (14)
M_i the matrix of that part. *Constant terms:* 3.13 needs
linearized inputs, but c_M is unchanged by â_i in
Q_i = Q_i^{lin} + â_i: in (Q_i)^{α_i} = Π_{j=0}^{e_i}
(Q_i^{lin}(X_i)^{2^j} + â_i^{2^j}), an F_2-deg-(e_i + 1) monomial
in X_i needs the linearized term in all e_i + 1 factors (a constant
factor gives F_2-deg 0, the rest <= e_i in total; reduction never
raises F_2-deg). By 3.14, need B_i with M_i'[𝒜_i, B_i] lower
triangular, nonzero diagonal (6.2), and every r ∈ {1, …, n-1}
killed by a zero column in some M_i'[𝒜_i, B_i - r] (6.3).

**Choice of B_i.** ℓ_1 = 0, ℓ_{i+1} = (ℓ_i + d_i - (e_i + 1))
mod n: ℓ_i = Σ_{i'<i} (d_{i'} - (e_{i'} + 1)) mod n. t_0 minimal
with S_{t_0+1} := Σ_{i'<=t_0} (d_{i'} - (e_{i'} + 1)) >= n - d_1 + 1
(exists by (†): S_{t+1} > n - d_1). i <= t_0: S_i <= n - d_1
(S_1 = 0), no wrap: ℓ_i = S_i <= n - d_1 <= n - d_i, plain integer.
Fix i <= t_0. a_{i,ℓ_i} != 0: B_i = {ℓ_i} ∪ {n - d_i + 1, …,
n - d_i + e_i} (1 + e_i elements; disjoint as ℓ_i <= n - d_i).
a_{i,ℓ_i} = 0: ℓ_i >= 1 (a_{i,0} != 0), ℓ_i < n - d_i
(a_{i,n-d_i} != 0); j_i := min{j > ℓ_i: a_{i,j} != 0},
j̄_i := max{j < ℓ_i: a_{i,j} != 0}: j̄_i < ℓ_i < j_i <= n - d_i.
q_i := j_i - j̄_i - 1 (zero-gap at ℓ_i), s_i := min{q_i, e_i}:
B_i = {j_i} ∪ {j̄_i + 1, …, j̄_i + s_i} ∪ {n - d_i + s_i + 1, …,
n - d_i + e_i} - sizes 1, s_i, e_i - s_i, disjoint
(j̄_i + s_i <= j_i - 1 < j_i <= n - d_i < n - d_i + s_i + 1); last
part empty if s_i = e_i. i > t_0: B_i = {0} ∪ {n - d_i + 1, …,
n - d_i + e_i} (first case at ℓ = 0); only 6.2 needed.

**Claim 6.2.** Every M_i'[𝒜_i, B_i] is, after reordering columns,
lower triangular with nonzero diagonal. *Proof.* (u, k)-entry a'_{i,(k-u) mod n}, u <= e_i <= d_i - 2 < n:
plain integers; fix i, write a'_j, d, e. Case a'_{ℓ_i} = 1
(incl. i > t_0, ℓ = 0): columns (ℓ_i, n-d+1, …, n-d+e); diagonal
indices ℓ_i, n-d (3.10(1)); above-diagonal indices in the zero
zone {n-d+1, …, n-d+e}. Case a'_{ℓ_i} = 0: columns (j_i,
j̄_i + 1, …, j̄_i + s_i, n - d + s_i + 1, …, n - d + e);
first-block diagonal j_i, j̄_i; above-diagonal indices in the
zero-gap (j̄_i, j_i), as j̄_i + s_i <= j_i - 1; rest as case 1. ∎

**Proposition 6.4 (BK Proposition 6.4).** If k ∈ B_i and
r ∈ [k + 1, k + d_i - (e_i + 1)] mod n, then M_i'[𝒜_i, B_i - r]
has an all-zero column. *Proof.* The k-column (a'_{i,k-r}, …, a'_{i,k-(r+e_i)})^T mod n
is all-zero if [k - (r + e_i), k - r] mod n ⊆ [n - d_i + 1, n - 1]. (15)
The r-range puts that window inside [1 - d_i, -1] ≡
[n - d_i + 1, n - 1] (mod n). ∎

**Claim 6.3.** For every r ∈ {1, …, n-1} some M_i'[𝒜_i, B_i - r]
has an all-zero column. *Proof.* By 6.4:
[1, n-1] ⊆ ∪_i ∪_{k∈B_i} [k + 1, k + d_i - (e_i + 1)] mod n. (16)
i <= t_0: ℓ_i = S_i plain; [ℓ_i + 1, ℓ_{i+1}] ends at S_{i+1}
mod n. i = 1: B_1 = {0} ∪ [n - d_1 + 1, n - d_1 + e_1]; k = 0
covers [1, d_1 - (e_1 + 1)] = [ℓ_1 + 1, ℓ_2]; top k's chain into
[n - d_1 + 2, n - 1]:
r ∈ [1, d_1 - (e_1 + 1)] ∪ [n - d_1 + 2, n - 1] covered. (17)
1 < i <= t_0: ∪_{k∈B_i} [k + 1, k + d_i - (e_i + 1)] ⊇
[ℓ_i + 1, ℓ_{i+1}]. (18)
Then (16): S_{t_0+1} < n => chain reaches ℓ_{t_0+1} >= n - d_1 + 1,
(17) the rest; >= n => the chain wraps over [1, n-1]. Cases for
(18). a_{i,ℓ_i} != 0: k = ℓ_i ∈ B_i gives it. a_{i,ℓ_i} = 0,
q_i < e_i: s_i = q_i, B_i = [j̄_i + 1, j_i] ∪ {n - d_i + q_i + 1,
…, n - d_i + e_i} (j̄_i + q_i + 1 = j_i); ℓ_i ∈ [j̄_i + 1,
j_i - 1] ⊆ B_i: same. a_{i,ℓ_i} = 0, q_i >= e_i: s_i = e_i,
B_i = {j_i} ∪ {j̄_i + 1, …, j̄_i + e_i}. *Subfield step:*
j̄_i < ℓ_i <= n - d_1 <= n - d_i => j̄_i + 1 <= n - d_i !=
n - d_i + 1; 3.10(2) in direction form (direction of A_i) -
vacuous in §7 (n prime; |Y_i| >= 2^{r+2} vs 2-element F_2
shifts/multiples). Zero run j̄_i + 1, …, j_i - 1 (length q_i),
not starting at n - d_i + 1, has no d_i - 1 consecutive indices:
q_i <= d_i - 2, j_i <= j̄_i + d_i - 1. k = j̄_i + 1, …, j̄_i + e_i
chain into [j̄_i + 2, j̄_i + d_i - 1] ∋ j_i (j̄_i + 2 <= ℓ_i + 1);
k = j_i extends to j_i + d_i - (e_i + 1) >= ℓ_i + d_i - (e_i + 1). ∎

*Conclusion of Theorem 6.1.* 6.2, 6.3 give 3.14's (1), (2)
(lower/upper immaterial): (14)'s r = 0 summand is a product of
nonzero triangular permanents; every r != 0 summand has a
vanishing factor. So c_M != 0, wt_2(β_i) = e_i + 1, and
|R(F^t)| > 1 by 3.1. ∎

## 7. Affine spaces of sublinear dimension: Theorem 7.3 and Theorem 2.2

Now n prime, r blocks, t < r; linearized polys. stay sums of
aX_i^{2^j} under reduction/Frobenius.

**Lemma 7.1 (Bourgain).** A ⊆ F^r F_2-affine, γ ∈ A => there are
linear Y_1, …, Y_r ⊆ F and F_2-linear σ_{ij}: Y_j → F (j < i) with
A = { (x_1, …, x_r): x_i = γ_i + y_i + Σ_{j<i} σ_{ij}(y_j),
y_i ∈ Y_i }, dim A = Σ_i dim Y_i.
*Proof.* A = γ + V. U_i := {v ∈ V: v_1 = … = v_{i-1} = 0},
Y_i := {v_i: v ∈ U_i}, s_i a linear section of U_i → Y_i,
σ_{ij}(y) := (s_j(y))_i. Greedy: killing coordinates in
order gives v = Σ_i s_i(y_i), unique by coordinate order;
hence ⊕ Y_i ≅ V, with coordinate i equal to y_i +
Σ_{j<i} σ_{ij}(y_j). Add γ. □

**Lemma 7.2.** P_i(X_i) ∈ F[X_i] linearized (i ∈ [r]); P_{ij}(X_j)
linearized (j < i); γ ∈ F^r; I_0 = {i_1 < … < i_t} ⊆ [r];
e_i >= 1 (i ∈ I_0); α_i = Σ_{k=0}^{e_i} 2^k; µ != 0. Let
g := Tr( µ Π_{i∈I_0} ( P_i(X_i) + Σ_{j<i} P_{ij}(X_j) + γ_i )^{α_i} ),
g' := Tr( µ Π_{i∈I_0} P_i(X_i)^{α_i} ), both mod
⟨X_i^{2^n} - X_i⟩. Then for all (β_i) with wt_2(β_i) = e_i + 1,
the coefficient of M = Π_{i∈I_0} X_i^{β_i} in g equals that in g'.
*Proof.* P̂_i := Σ_{j<i} P_{ij}(X_j) + γ_i (no X_i). Expanding
Tr, (z^{α_i})^{2^s} = Π_{l=0}^{e_i} z^{2^{s+l}}:
g = Σ_{s=0}^{n-1} µ^{2^s} Π_{i∈I_0} Π_{l=0}^{e_i} ( P_i(X_i) +
P̂_i )^{2^{s+l}} reduced. (19)–(20)
Freshman split per factor; choosing P̂ exactly at L ⊆ {0, …, e_i}:
Π_{l=0}^{e_i} ( P_i(X_i) + P̂_i )^{2^{s+l}} = ( P_i(X_i) )^{Σ_l
2^{s+l}} + G_i(X_{<=i}), (21)–(23)
L != ∅ terms of F_2-deg <= e_i + 1 - |L| <= e_i in X_i. Peel
i_t, …, i_1 in decreasing order: at i_t the other factors omit
X_{i_t} and G_{i_t}(rest) has F_2-deg <= e_{i_t} < e_{i_t} + 1
(M's degree), so only the pure term contributes; likewise down.
Sum over s. □

**Theorem 7.3.** t < r, n prime, n >= r(r+1)/t. i ∈ [r]:
e_i = r + 1 - i, α_i = Σ_{k=0}^{e_i} 2^k, µ != 0,
f(x_1, …, x_r) = Tr( µ Σ_{I⊆[r], |I|=t} Π_{i∈I} x_i^{α_i} ).
Then every affine A ⊆ F^r with dim A > nr/t + nt + r(r+1) has
|f(A)| > 1. *Proof.* (1) γ ∈ A; decompose by 7.1. Q_i := image poly. of
Y_i (3.7); Q_{ij} := linearized poly. of x -> σ_{ij}(Q_j(x))
(3.4; equality since Q_j(x) ∈ Y_j).
R_i := Q_i(X_i) + Σ_{j<i} Q_{ij}(X_j) + γ_i; H := (R_1, …, R_r)
has H(F^r) = A. h := f ∘ H reduced:
h = Σ_{|I|=t} Tr( µ Π_{i∈I} R_i^{α_i} ) mod ⟨X_i^{2^n} - X_i⟩. (24)–(25)
(2) d_i := dim Y_i; Σ d_i > nr/t + nt + r(r+1).
S := {i: d_i > r + 1}: |S| >= t, else Σ d_i <= |S|n + r(r+1) <=
(t-1)n + r(r+1) < nt + r(r+1) + nr/t. Σ_{i∈S} d_i >= Σ_{[r]}
(d_i - (r+1)) > nr/t + nt. Averaging t-subsets I_0 ⊆ S:
Σ_{i∈I_0} (d_i - (r+1)) >= (t/|S|)Σ_{i∈S} d_i - (r+1)t >=
(t/r)(nr/t + nt) - (r+1)t = n + nt²/r - (r+1)t >= n, (26)
by n >= r(r+1)/t. (3) g := I_0-term of (24)–(25);
g' := 7.2's g' at P_i = Q_i. 6.1 to Y_i
(i ∈ I_0), e_i = r + 1 - i: d_i > r + 1 >= 3 > 1;
1 <= e_i < d_i - 1 (e_i + 1 <= r + 1 < d_i);
Σ_{I_0} (d_i - (e_i + 1)) >= Σ_{I_0} (d_i - (r+1)) >= n (26) > n -
max d_i + 1; subfields: n prime => only F_2 (2-element shifts),
|Y_i| = 2^{d_i} >= 2^{r+2}. 6.1 gives M = Π_{i∈I_0} X_i^{β_i},
wt_2(β_i) = e_i + 1, coeff != 0 in g'; 7.2: same coeff in g.
(4) R_j is a sum of aX_i^{2^b} (i <= j) and a constant:
F_2-deg <= 1 in X_{<=j}. For J != I_0 contributing to M, pick in
Π_{j∈J} Π_{l=0}^{e_j} R_j^{2^{l+w}} (some w) one monomial per
(j, l) with product M; λ_j(l) := its variable index (<= j), or j
if constant (counts below only grow; Σ_{i∈I_0} |λ_j^{-1}(i)| <=
r + 2 - j, so 7.4 applies). F_2-deg of M in X_i <= #{(j, l):
λ_j(l) = i}, and equals e_i + 1 = r + 2 - i for i ∈ I_0:
Σ_{j∈J} |λ_j^{-1}(i)| >= r + 2 - i ∀i ∈ I_0. (27)

**Claim 7.4.** Let I_0 ⊆ [r], |I_0| = t. If J ⊆ [r], |J| = t, and
λ_j: {0, …, r+1-j} → {1, …, j} for j ∈ J satisfy (27), then
J = I_0. *Proof.* Per j: Σ_{i∈I_0} |λ_j^{-1}(i)| <= r + 2 - j
(domain size). For k >= 0, Σ_j Σ_l (j^k - λ_j(l)^k) >= 0
(λ_j(l) <= j); expanding, and bounding via (27)
Σ_j Σ_l λ_j(l)^k = Σ_i i^k|{(j, l): λ_j(l) = i}| >=
Σ_{i∈I_0} (r+2-i)i^k, gives Σ_{i∈I_0} i^{k+1} - Σ_{j∈J} j^{k+1} >=
(r+2)( Σ_{i∈I_0} i^k - Σ_{j∈J} j^k ) >= 0 by induction on k
(k = 0: equal cardinalities). Hence
Σ_{i∈I_0} i^k >= Σ_{j∈J} j^k for all k. (28)
i_1 := max I_0, j_1 := max J. k → ∞ in (28) gives i_1 >= j_1 (else
j_1^k > ti_1^k dominates); (27) gives some λ_j(l) = i_1, j <= j_1,
so i_1 <= j_1: i_1 = j_1. (27) at i_1 then forces
|λ_{j_1}^{-1}(i_1)| >= r + 2 - j_1 = domain size: λ_{j_1}(l) = i_1
∀l. Delete i_1, j_1: (28) persists for the deleted sets; repeat
for the next-largest (i_2 <= j_2 since j_1 is exhausted by
i_1 != i_2): i_2 = j_2. Induction on rank: I_0 = J. □

So M occurs in h only from I_0, with nonzero coefficient
(positive degree: wt_2(β_i) >= 2); 3.1 gives |f(A)| > 1. □

**Theorem 2.2 (deduction from 7.3).** Fix m; n := least prime
> 2m^{3/5}, r := ⌈m/n⌉, t := ⌈√r⌉. Bertrand: a prime in
(2m^{3/5}, 4m^{3/5}]: 2m^{3/5} < n <= 4m^{3/5}; hence
m^{2/5}/4 <= m/n <= r <= m/n + 1 < m^{2/5}/2 + 1;
t = ⌈√r⌉ < r for r >= 3. *Block length:* t >= √r:
r(r+1)/t <= (r+1)√r = (r+1)^{3/2} <= (m^{2/5}/2 + 2)^{3/2} =
(1/(2√2)) m^{3/5} (1 + 4m^{-2/5})^{3/2} = (1/(2√2)) m^{3/5}
(1 + o(1)) < 2m^{3/5} < n for large m (1/(2√2) ≈ 0.354 < 2):
n >= r(r+1)/t. *Threshold, jointly via r <= m/n + 1:*
nr/t <= n√r = (nm + n²)^{1/2} <= (4m^{8/5} + 16m^{6/5})^{1/2} =
2m^{4/5}(1 + o(1)); nt <= n(√r + 1) <= 2m^{4/5}(1 + o(1)) +
4m^{3/5}; r(r+1) <= (m^{2/5}/2 + 1)(m^{2/5}/2 + 2) =
(1/4)m^{4/5}(1 + o(1)). Sum: 4.25 m^{4/5} (1 + o(1)) < 6m^{4/5}
for large m: every A ⊆ F_2^m, dim A > 6m^{4/5}, meets 7.3's
hypothesis. *Identification:* by 3.3, 2.2's π is Tr(µ·), µ != 0
unique - 7.3's form. Exponents: 2.2's (4) gives block i
Σ_{j=0}^{i} 2^j; 7.3 gives Σ_{k=0}^{e_i} 2^k, e_i = r + 1 - i,
i.e. Σ_{j=0}^{r+1-i} 2^j. Renaming blocks i -> r + 1 - i:
f_{7.3}(x_1, …, x_r) = f_{2.2}(x_r, …, x_1); coordinate
permutation preserves affinity and dimension. So 2.2's f is
nonconstant on every affine subspace of dim > 6m^{4/5}, all
large m. ∎

*Application (Main Lemma).* Rename (m, n) → (N, k). φ: F_2^N →
K^r = (F_{2^k})^r is `proof.md` §1's pad-and-split map (append
rk - N zeros; k-bit blocks; polynomial basis of K): injective
linear. 2.2 at m := N quantifies over all injective linear φ:
φ admissible, parameters k, r, t from the true length N, short
last block (rk - N < k zeros) needing no divisibility. Exponents
agree (Σ_{j=0}^{i} 2^j = 2^{i+1} - 1); π (constant coefficient)
is nontrivial F_2-linear. For affine A ⊆ F_2^N, dim A >
6N^{4/5}: φ(A) affine, same dim, |f_N(A)| = |f(φ(A))| > 1 by
2.2, all large N. N <= 7776 = 6^5: 6N^{4/5} >= N, vacuous. Hence
f_N disperses for d(N) = ⌊6N^{4/5}⌋ + 1 at all large N - what
`proof.md` §2 uses. □
