# Disperser, part A: subspace polynomials (Ben-Sasson–Kopparty, §3)

Reproduces, with proofs, the §3 material of Ben-Sasson–Kopparty,
*Affine Dispersers from Subspace Polynomials*, on which the disperser
argument of `proof.md` §2 rests; statements keep BK's numbering.

**Conventions for this file.** p = 2 throughout; n is the block length,
F = F_{2^n}, F̄ an algebraic closure of F, coefficient indices mod n.
Q(X) = Σ_{j=0}^{n-1} a_j X^{2^j} is *linearized*, and *affine linearized*
with a constant term â allowed: Q(X) = Σ_j a_j X^{2^j} + â. Linearized
polynomials induce F_2-linear maps F → F; affine linearized ones are
exactly those of F_2-degree at most 1. For affine L ⊆ F of dimension d,
P_L (Definition 3.5) is its kernel-subspace polynomial and Q_L
(Lemma 3.7) its image-subspace polynomial. The *F_2-degree weight* of m
is wt_2(m) = sum of the binary digits of m.

**Proposition 3.1.** Let h(X_1, …, X_r) ∈ F[X_1, …, X_r] and let
h′ = h mod ⟨X_i^{2^n} − X_i⟩_{i∈[r]} (the unique representative of degree
< 2^n in each variable). Then h(x) = h′(x) for all x ∈ F^r. Consequently
|h(F^r)| > 1 iff h′ has positive degree.

*Proof.* a^{2^n} = a on F, so X_i^{2^n} ↦ X_i preserves values on F^r,
and h′ is obtained from h by such replacements. Constant h′ gives constant
h on F^r. Nonconstant h′ has h′ − c a nonzero polynomial for every c, and
a nonzero polynomial of degree < 2^n in each variable cannot vanish on
all of F^r (r = 1: < 2^n roots; induct on r): so h′, hence h, takes more
than one value on F^r. □

**Proposition 3.3.** Let Tr(Y) = Σ_{i=0}^{n-1} Y^{2^i} be the trace from F
to F_2. For every F_2-linear map π: F → F_2 there is a unique µ ∈ F with
π(x) = Tr(µx) for all x ∈ F. Furthermore π is trivial iff µ = 0.

*Proof.* (x, y) ↦ Tr(xy) is nondegenerate: Tr, a nonzero polynomial of
degree 2^{n-1} < 2^n, has Tr(z) ≠ 0 for some z; for µ ≠ 0, x = z/µ gives
Tr(µx) ≠ 0. Hence µ ↦ Tr(µ·): F → Hom_{F_2}(F, F_2) is injective
F_2-linear between spaces of dimension n — an isomorphism (existence
and uniqueness); π ≡ 0 iff µ = 0 by the same nondegeneracy. □

**Lemma 3.4.** Let φ: F_2^n → F be an F_2-linear isomorphism. There is a
one-to-one correspondence between affine transformations T: F_2^n → F_2^n
and affine linearized polynomials in F[X]: for every affine T there is a
unique affine linearized polynomial P_T with P_T(φ(b)) = φ(T(b)) for all
b ∈ F_2^n.

*Proof.* Transport by φ to affine maps S: F → F. An affine linearized
polynomial induces an F_2-affine map F → F. Both form F_2-spaces of
dimension n(n + 1): the polynomials via their coefficients
a_0, …, a_{n-1}, â ∈ F (n each); the maps since n^2 + n = n(n + 1) (an
F_2-matrix plus a translation). Evaluation, F_2-linear in the
coefficients, is injective: if Σ a_j X^{2^j} + â vanishes on F, the
linearized part is constant as a function, so â = 0 (at 0) and it
vanishes on F; a nonzero linearized polynomial has degree
≤ 2^{n-1} < 2^n, hence < 2^n roots: all a_j = 0. Injective, equal
finite dimensions: an isomorphism. □

**Definition 3.5 (kernel-subspace polynomial).** For an affine subspace
L ⊆ F of dimension d, the kernel-subspace polynomial of L is

    P_L(X) = ∏_{α∈L} (X − α).

**Lemma 3.6.** If L ⊆ F is an affine subspace of dimension d, then P_L is
a monic affine linearized polynomial of degree 2^d. Furthermore, P_L is
linearized iff L is a linear space.

*Proof* (BK state this without proof; the Moore determinant). Let V be
linear of dimension d, basis β_1, …, β_d; Δ(u_1, …, u_d) =
det(u_j^{2^{i-1}})_{1≤i,j≤d} is the Moore determinant.
Induction on d: Δ as a function of u_d is a nonzero linearized polynomial
vanishing exactly on span(β_1, …, β_{d-1}) — it has those 2^{d-1} roots,
degree 2^{d-1} in u_d, and leading coefficient Δ(β_1, …, β_{d-1}) ≠ 0 —
hence Δ(β_1, …, β_d) ≠ 0. Set

    P(X) = Δ(X, β_1, …, β_d) / Δ(β_1, …, β_d),

the (d + 1) × (d + 1) determinant with columns X, β_1, …, β_d and rows
2^0, …, 2^d. First-column expansion gives
P(X) = Σ_{i=0}^d c_i X^{2^i} with X^{2^d}-coefficient Δ(β)/Δ(β) = 1:
P is monic linearized of degree 2^d. For x ∈ V the first column is an
F_2-combination of the β_j-columns (Frobenius is F_2-linear), so P(x) = 0;
monic of degree 2^d = |V| and vanishing on V, P = ∏_{α∈V}(X − α) = P_V.
For affine L = γ + V, P_V(X) − P_V(γ) is monic affine
linearized of degree 2^d with root set {x : P_V(x) = P_V(γ)} = γ + ker P_V
= γ + V = L, hence equals P_L; its linearized part is P_V. Finally, P_L is
linearized iff its constant term P_V(γ) is 0, iff γ ∈ ker P_V = V, iff
0 ∈ L, iff L is linear. □

**Lemma 3.7 (image-subspace polynomial).** If L ⊆ F is an affine subspace
of dimension d, there exists a monic affine linearized polynomial Q_L(X)
of degree 2^{n−d}, the image-subspace polynomial of L, such that
Q_L(F) = L and

    P_L(Q_L(X)) = X^{2^n} − X.                                        (6a)

If L is linear, then Q_L is linearized (not merely affine linearized),
ker Q_L = P_L(F), Q_L has 2^{n−d} roots, and also

    Q_L(P_L(X)) = X^{2^n} − X.                                        (6b)

*Proof.* Suppose first that L is linear. Let L_0 = P_L(F); as P_L is
linearized with kernel L, L_0 is a linear space of dimension n − d. Set
Q_L = P_{L_0}, monic linearized of degree 2^{n−d} by Lemma 3.6. The
composition Q_L ∘ P_L is monic of degree 2^n and vanishes on F (since
P_L(x) ∈ L_0, the root set of Q_L), so Q_L ∘ P_L = ∏_{α∈F}(X − α) =
X^{2^n} − X. Let g(Y) = P_L(Q_L(Y)) − (Y^{2^n} − Y). Then
g(P_L(X)) = P_L(X^{2^n} − X) − (P_L(X)^{2^n} − P_L(X)) = 0, because P_L,
being linearized, satisfies P_L(X^{2^n} − X) = P_L(X^{2^n}) − P_L(X) and
P_L(X^{2^n}) = P_L(X)^{2^n} (the coefficients lie in F, hence are fixed by
the 2^n-Frobenius). Since deg(g ∘ P_L) = (deg g)·2^d, g ≡ 0, proving
(6a); in particular image(Q_L) ⊆ ker P_L = L, and ker Q_L ⊇ P_L(F) = L_0
with both of dimension n − d, so |image Q_L| = 2^d = |L| and
Q_L(F) = L; also ker Q_L = L_0, Q_L has 2^{n−d} roots, and (6b) is the
identity Q_L ∘ P_L = X^{2^n} − X just proved.
For general affine L = γ + V with V its direction, let Q_V be as above
and set Q_L(X) = Q_V(X) + γ: monic affine linearized of degree 2^{n−d}
with Q_L(F) = V + γ = L and linearized part Q_V. By the proof of
Lemma 3.6, P_L = P_V + P_V(γ) (characteristic 2), and since P_V is
linearized, P_L(Q_L(Y)) = P_V(Q_V(Y) + γ) + P_V(γ) = P_V(Q_V(Y)) =
Y^{2^n} − Y, giving (6a). (For affine L the identity (6b) holds only up
to the constant γ: Q_L(P_L(X)) = X^{2^n} − X + γ, by the same computation
and Q_V(P_V(γ)) = γ^{2^n} − γ = 0. All applications below use (6a), the
image identity Q_L(F) = L, and — for linear L — the linearizedness of
Q_L.) □

**Claim 3.8.** Let k ≥ 1 and suppose a, c ∈ F satisfy a^{2^k} = ca. Let b
be any (2^k − 1)-th root of c in F̄. Then a ∈ b·F_{2^k}.

*Proof.* If a = 0, clear. Else a^{2^k − 1} = c gives
(a/b)^{2^k − 1} = c/c = 1, i.e. (a/b)^{2^k} = a/b: a/b ∈ F_{2^k}. □

**Claim 3.9.** For an affine linearized polynomial
Q(X) = Σ_{j=0}^{n-1} a_j X^{2^j} + â and any integer s,

    (Q(X))^{2^s} mod (X^{2^n} − X) = Σ_{j=0}^{n-1} a_{j−s}^{2^s} X^{2^j} + â^{2^s},

indices mod n.

*Proof.* F_2-linearity of Z ↦ Z^{2^s}:
(Q(X))^{2^s} = Σ_j a_j^{2^s} X^{2^{j+s}} + â^{2^s}; reducing
X^{2^{j+s}} mod (X^{2^n} − X) replaces the index j + s by (j + s) mod n,
since X^{2^n} ≡ X; reindex j ↦ j − s. □

**Lemma 3.10 (main structural lemma).** Let L be a d-dimensional linear
subspace of F and let Q_L(X) = Σ_{j=0}^{n-1} a_j X^{2^j} be its
image-subspace polynomial.

1. For any integer r, the set J = {(r + j) mod n : j = 0, …, d − 1} of d
   consecutive indices in Z_n contains some j with a_j ≠ 0. In particular
   a_0 ≠ 0 and a_{n−d} ≠ 0.
2. Suppose L is not contained in any constant multiple of a proper
   subfield of F, i.e. L ⊄ β·F_{2^k} for any β ∈ F and any k < n. Then
   for any integer r ≠ n − d + 1, the set J = {(r + j) mod n :
   j = 0, …, d − 2} of d − 1 consecutive indices contains some j with
   a_j ≠ 0.

*Proof.* If d = n then L = F, Q_L = X, a_0 = 1: both parts immediate.
Assume d < n; note a_j = 0 for n − d + 1 ≤ j ≤ n − 1 (degree 2^{n−d})
and a_{n−d} = 1 (monic).

Part 1: suppose a_j = 0 for all j ∈ J. Set s ≡ −(r + d) (mod n),
s ∈ {0, …, n − 1}, and Q′(X) = (Q_L(X))^{2^s} mod (X^{2^n} − X). By
Claim 3.9 its X^{2^j}-coefficient is a_{j+r+d}^{2^s} (indices mod n); as
j runs over [n − d, n − 1], (j + r + d) mod n runs over J, so its
coefficients at indices ≥ n − d vanish: deg Q′ ≤ 2^{n−d−1}. Q′ ≠ 0, as
a_{n−d} = 1 occurs among its coefficients. By Proposition 3.1, Q′
vanishes on ker Q_L (Q′(α) = Q_L(α)^{2^s}), which has 2^{n−d} elements
(Lemma 3.7, L linear): a nonzero polynomial of degree ≤ 2^{n−d−1}
cannot have 2^{n−d} roots — contradiction. "In particular": apart
from a_0 resp.
a_{n−d}, the windows {n − d + 1, …, n − 1, 0} and {n − d, …, n − 1}
contain only coefficients vanishing by degree, so a_0 ≠ 0 and
a_{n−d} ≠ 0.

Part 2: suppose a_j = 0 for all j ∈ J = {r, …, r + d − 2},
r ≢ n − d + 1 (mod n). Set k ≡ n − (r + d) + 1 (mod n), so
k ∈ {1, …, n − 1} (k = 0 iff r ≡ n − d + 1), and
Q′(X) = (Q_L(X))^{2^k} mod (X^{2^n} − X). By Claim 3.9 its
X^{2^j}-coefficient is a_{j−k}^{2^k}; for j ∈ [n − d + 1, n − 1],
(j − k) mod n runs over J, so those coefficients vanish:
deg Q′ ≤ 2^{n−d}; Q′ ≠ 0 as before. By Proposition 3.1, Q′ vanishes on
ker Q_L. Choose c ∈ F with deg(Q′ − c·Q_L) ≤ 2^{n−d−1} (match the
leading coefficients; if deg Q′ < 2^{n−d}, take c = 0 — the root count
below then gives Q′ ≡ 0, already a contradiction). Q′ and Q_L both
vanish on ker Q_L, which has 2^{n−d} elements, so Q′ − c·Q_L ≡ 0:
Q′ = c·Q_L. Thus Q_L(α)^{2^k} = c·Q_L(α) for all α ∈ F — i.e.
y^{2^k} = cy for all y in the image L. Since d ≥ 1, some y ≠ 0
exists: c ≠ 0. By Claim 3.8, L ⊆ b·F_{2^k} for a (2^k − 1)-th
root b of c in F̄. Fix β ∈ L \ {0}: for any α ∈ L, α/β ∈ F_{2^k} (a
ratio of two elements of b·F_{2^k}) and α/β ∈ F, so
β^{−1}·L ⊆ F_{2^k} ∩ F = F_{2^{gcd(k,n)}}. As 1 ≤ k ≤ n − 1,
gcd(k, n) < n: F_{2^{gcd(k,n)}} is a proper subfield of F and
L ⊆ β·F_{2^{gcd(k,n)}} — contradicting the hypothesis. □

**Remark (affine form used in §6).** In §6, Lemma 3.10 is applied to the
image-subspace polynomials Q_i of *affine* spaces A_i of dimension d_i,
written Q_i(X) = Σ_{j=0}^{n-1} a_{i,j} X^{2^j} + â_i as in §6. The
coefficients a_{i,j} are those of the linearized part of Q_i, which by
Lemma 3.7 is the image-subspace polynomial of the direction of A_i —
a linear space — so Lemma 3.10 applies to them. Conclusions:
a_{i,j} = 0 for j ∈ [n − d_i + 1, n − 1] (since deg Q_i = 2^{n−d_i});
a_{i,0} ≠ 0 and a_{i,n−d_i} ≠ 0; and part 2's hypothesis in the
direction form: the direction of A_i is not contained in any constant
multiple of a proper subfield of F. (This is the form Claim 6.3 of
`disperser-b.md` consumes. It is not implied by Theorem 6.1's
affine-shift hypothesis for arbitrary A_i, but in the §7 application
n is prime and both forms are vacuous for the spaces at hand.)

**Definition 3.11 (associated matrix).** For a linearized polynomial
Q(X) = Σ_{j=0}^{n-1} a_j X^{2^j} over F, its associated matrix M_Q is the
n × n matrix over F whose (i, j)-entry, rows and columns indexed by
{0, …, n − 1} with index arithmetic mod n, is

    (M_Q)_{i,j} = (a_{j−i})^{2^i} :

row 0 is (a_0, a_1, …, a_{n−1}), row 1 is ((a_{n−1})^2, (a_0)^2, …,
(a_{n−2})^2), and so on. The zero-one indicator matrix M_Q′ has (i, j)-entry
a′_{j−i}, where a′_m = 0 if a_m = 0 and a′_m = 1 otherwise; thus the
(i, j)-entry of M_Q′ indicates whether the (i, j)-entry of M_Q is nonzero.
For nonempty A, B ⊆ {0, …, n − 1} with |A| = |B|, M[A, B] denotes the
|A| × |A| minor on rows A and columns B, and B − r = {s − r mod n : s ∈ B}.

**Claim 3.12.** The (i, j)-entry of M_Q is the coefficient of X^{2^j} in
(Q(X))^{2^i} mod (X^{2^n} − X).

*Proof.* Claim 3.9 with s = i: the coefficient is
a_{j−i}^{2^i} = (M_Q)_{i,j}, by definition. □

**Lemma 3.13 (permanent formula).** Let µ ∈ F \ {0}. Let
A_1, …, A_t, B_1, …, B_t ⊆ {0, …, n − 1} satisfy |A_i| = |B_i| > 0, and
set α_i = Σ_{j∈A_i} 2^j, β_i = Σ_{k∈B_i} 2^k. Let Q_1, …, Q_t be
linearized polynomials with associated matrices M_1, …, M_t. Then the
coefficient c_M of the monomial M = ∏_{i=1}^t X_i^{β_i} in

    R(X_1, …, X_t) = Tr( µ · ∏_{i=1}^t (Q_i(X_i))^{α_i} )
                     mod ⟨X_i^{2^n} − X_i⟩_{i∈[t]}                     (7)

is given by

    c_M = Σ_{r=0}^{n-1} µ^{2^r} · ∏_{i=1}^t Perm( M_i[A_i + r, B_i] )
        = Σ_{r=0}^{n-1} µ^{2^r} · ∏_{i=1}^t Perm( M_i[A_i, B_i − r] )^{2^r}.   (8)

*Proof.* Since Tr(Z) = Σ_{r=0}^{n-1} Z^{2^r},

    Tr( µ · ∏_i (Q_i(X_i))^{α_i} ) = Σ_{r=0}^{n-1} µ^{2^r} · ∏_i ((Q_i(X_i))^{α_i})^{2^r}.

Let c_M^{(r)} be M's coefficient in the r-th summand, so
c_M = Σ_r c_M^{(r)} with c_M^{(r)} = µ^{2^r} times, over i, the
X_i^{β_i}-coefficient of ((Q_i(X_i))^{α_i})^{2^r}. Fix i, dropped from
the notation. Since α = Σ_{j∈A} 2^j,

    ((Q(X))^{α})^{2^r} = ∏_{j∈A} (Q(X))^{2^{j+r}} = ∏_{j∈A+r} (Q(X))^{2^j}.

With β = Σ_{k∈B} 2^k: X^β arises by choosing h(k) ∈ A + r for each
k ∈ B and taking X^{2^k} from the factor (Q(X))^{2^{h(k)}}; as
|A| = |B|, the contributing choices are exactly the bijections
h: B → A + r, and by Claim 3.12 the term from factor h(k) has
coefficient (M)_{h(k),k} = a_{k−h(k)}^{2^{h(k)}}. Summing over h gives
the permanent:

    Σ_{h: B → A+r} ∏_{k∈B} (M)_{h(k),k} = Perm( M[A + r, B] ).

So c_M^{(r)} = µ^{2^r} · ∏_i Perm( M_i[A_i + r, B_i] ): (8), first
expression. For the second:
(M_i)_{j+r,k} = (a_{i,k−j−r})^{2^{j+r}} = ((M_i)_{j,k−r})^{2^r}, so
M_i[A_i + r, B_i] is M_i[A_i, B_i − r] under entrywise 2^r-power; the
permanent, a polynomial in the entries, commutes with this Frobenius:
Perm( M_i[A_i + r, B_i] ) = Perm( M_i[A_i, B_i − r] )^{2^r}. □

**Remark 3.14 (usage criterion).** Keep the notation of Lemma 3.13.
Suppose:

1. each M_i′[A_i, B_i] is, up to reordering of rows and columns,
   triangular with every diagonal entry nonzero; and
2. for every r ∈ {1, …, n − 1} there is some i with M_i′[A_i, B_i − r]
   containing an all-zero column.

Then c_M ≠ 0. Indeed, by condition 1 and (8), the r = 0 summand
µ · ∏_i Perm( M_i[A_i, B_i] ) is nonzero: the permanent of a triangular
matrix is the product of its diagonal entries, nonzero by the indicator
condition, and µ ≠ 0. By condition 2, every summand with r ≠ 0 has a
factor the permanent of a matrix with an all-zero column, hence
vanishes. Hence c_M is the r = 0 summand, nonzero. □
