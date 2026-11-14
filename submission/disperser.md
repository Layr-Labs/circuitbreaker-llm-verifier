# The affine disperser: complete proof for p = 2

This file is supplied evidence for the joint theorem declared in `bound.json`.
It proves that the family defined in `field.md` is an affine disperser for
dimension d(n) = floor((77/25)*n^(4/5)) + 1 at every input length
n >= N0 := 10^20.
The construction is the Ben-Sasson-Kopparty construction (SIAM J. Comput.
41(4):880-914, 2012), specialized to p = 2; every proof is reproduced in
full, and no external premise is used beyond elementary finite-field theory
proved here or in `field.md`.

Throughout: n is the actual input length; k = k(n) is the odd prime and
K = K_n = F_(2^k) the field produced by `field.md` (for k prime the only
proper subfield of K is F_2, since F_(2^a) contains F_(2^b) iff b divides a);
r = r(n) = ceil(n/k); t = t(n) = ceil(sqrt(r)). Recall from `field.md` the
injective F_2-linear map phi = (phi_1,...,phi_r): F_2^n -> K^r (zero-pad to
r*k bits, split into r blocks) and the nonzero F_2-linear map pi: K -> F_2
returning the coefficient of Z^0. The disperser is

    f_n(x) = pi( Sym^t_r(y_1, ..., y_r) ),  y_i = phi_i(x)^(2^(i+1)-1),

where Sym^t_r = sum over t-element subsets S of [r] of prod_{i in S} y_i.

An affine subspace of F_2^n is b + V for a linear subspace V; its dimension
is dim V. g: F_2^n -> F_2 is an affine disperser for dimension D if g is
nonconstant on every affine subspace of dimension at least D.

**Theorem A.** For every integer n >= 10^20, f_n is an affine disperser for
dimension d(n) = floor((77/25)*n^(4/5)) + 1.

## 1. Polynomial facts over K

**Proposition 1.1.** Let h in K[X_1,...,X_s] and h_red the unique polynomial
congruent to h modulo the ideal (X_i^(2^k) - X_i : i) with degree less than
2^k in each variable (obtained by repeatedly replacing X_i^(2^k) by X_i; the
representative with all exponents below 2^k is unique). Then
h(x) = h_red(x) for all x in K^s, and the monomials X_1^(b_1)...X_s^(b_s),
0 <= b_i < 2^k, are linearly independent as functions on K^s. Hence if
h_red has a positive-degree monomial with nonzero coefficient, h is
nonconstant on K^s.

*Proof.* The polynomial X_i^(2^k) - X_i vanishes at every point of K (all
2^k elements of K are roots of X^(2^k) - X), so each replacement preserves
all values on K^s. Independence: by induction on s. For s = 1: a linear
combination F = sum_b c_b X^b of degree < 2^k vanishing on K is a
polynomial with 2^k roots and degree < 2^k, hence the zero polynomial and
all c_b = 0. For s >= 2: write F = sum_{b_1} X_1^(b_1) G_{b_1}(X_2..X_s)
with G_{b_1} the coefficient polynomials; for each fixed (a_2..a_s) in
K^(s-1), the univariate polynomial
sum_{b_1} G_{b_1}(a_2..a_s) X_1^(b_1) has degree < 2^k and vanishes on K,
so it is the zero polynomial, whence G_{b_1}(a_2..a_s) = 0 for every b_1
and every point of K^(s-1); by induction each G_{b_1} is the zero
polynomial. Nonconstancy: if h_red has a positive-degree monomial with nonzero
coefficient c* X^b* and h were constant on K^s with value c, then the
polynomial h_red - c vanishes on K^s, so by the independence just proved
all its coefficients vanish — contradicting the coefficient c* of the
positive-degree monomial. Hence h, which agrees with h_red on K^s, is
nonconstant on K^s. qed

For an integer i >= 0 let wt(i) be the number of ones in its binary
expansion. Uniqueness of binary representation gives: a sum of distinct
binary powers determines its summands, and wt(sum of distinct powers) = the
number of summands.

**Lemma 1.3 (weight inequality under reduction).** Let
E = sum_{a in S} 2^(h_a) be a sum of |S| >= 1 binary powers (h_a >= 0,
repetitions allowed), and let q := E mod (2^k - 1) with 0 <= q < 2^k - 1
(so q = 0 iff 2^k - 1 divides E). Then wt(q) <= |S|. Moreover, if all
h_a < k, then: (1) if the h_a are pairwise distinct and |S| <= k - 1, then
q = E and wt(q) = |S|; (2) if some two h_a coincide, then wt(q) < |S|;
(3) if |S| <= k - 1, then q != 0.

*Proof.* Write C_j := #{a : h_a = j (mod k)} for 0 <= j < k, so
E = sum_j C_j 2^j modulo 2^k - 1 and sum_j C_j = |S|. Apply
the carry rule "two units at position j become one unit at position j+1"
(indices cyclic, using 2 * 2^(k-1) = 2^k = 1), which preserves the
represented residue and strictly decreases the total count sum_j C_j. Since
the count is a positive integer, only finitely many carries occur, and the
process stops exactly when every count is 0 or 1. Let the terminal counts be
c_j in {0,1}: the terminal configuration represents q (a {0,1}-valued
configuration with total count c represents sum c_j 2^j, an integer in
[0, 2^k - 1] whose residue is q; since 0 <= q < 2^k - 1, either
sum c_j 2^j = q, or sum c_j 2^j = 2^k - 1 (all ones) and q = 0), and its
total count is at most |S|. If q > 0, the representing {0,1}-configuration
is unique (two distinct subsets of {0,...,k-1} with sums in [0, 2^k - 1]
that are congruent mod 2^k - 1 are either equal or one has sum 0 and the
other 2^k - 1), so its count is wt(q) and wt(q) <= |S|. If q = 0, then
wt(q) = 0 <= |S|. Claim (1): distinct h_a < k give all C_j in {0,1}
already, so no carry occurs and the terminal configuration is the original
one; E is the sum of |S| <= k - 1 distinct powers 2^j with j < k, so
1 <= E <= 2^1 + ... + 2^(k-1) = 2^k - 2 < 2^k - 1, E is its own residue q,
and wt(q) = wt(E) = |S|. Claim
(2): if two h_a coincide, some C_j >= 2, so at least one carry occurs, the
terminal count is at most |S| - 1, and wt(q) <= terminal count <= |S| - 1
(when q > 0 the terminal count is exactly wt(q); when q = 0, wt(q) = 0).
Claim (3): if q = 0, the terminal configuration is the all-ones state
(the empty configuration is unreachable from a positive count, and a
{0,1}-configuration represents 0 only if it is empty or all-ones), whose
count is k, so k <= terminal count <= |S| <= k - 1, a contradiction.
This proves the lemma. qed

**Proposition 1.2.** (i) The functions X -> X^(2^j), 0 <= j <= k-1, together
with 1, are linearly independent over K as functions on K, and every F_2-affine
map L: K -> K equals sum_j a_j X^(2^j) + c for unique a_j, c in K (the
displayed family has F_2-dimension k(k+1) — k coefficients a_j and c in K —
which equals the dimension k^2 + k of the space of F_2-affine maps K -> K;
by independence the evaluation map is injective, hence onto; restricting to
c = 0, the linearized family has dimension k^2 = dim of Hom_{F_2}(K, K), so
every F_2-linear map is a unique linearized polynomial). A polynomial
sum_j a_j X^(2^j) (respectively + c) is called linearized (respectively
affine linearized); it is additive up to the constant: P(u+v) = P(u) +
P(v) + P(0) over any extension field (freshman's dream).
(ii) Tr(a) := sum_{j=0}^{k-1} a^(2^j) is a nonzero F_2-linear functional on
K, and pi(a) = Tr(mu*a) for a unique mu in K; mu != 0 iff pi is nonzero.

*Proof.* (i) A nontrivial combination sum_j a_j X^(2^j) + c is a nonzero
polynomial of degree at most 2^(k-1) < 2^k, which cannot vanish at all 2^k
points of K; this gives independence. The maps of the displayed form make up
an F_2-linear space of dimension k(k+1) (k coefficients a_j in K, plus c in
K), equal to dim_F2 of the space of F_2-affine maps K -> K (k^2 linear maps
plus a k-dimensional constant term); by independence the evaluation map is
injective, hence onto. Additivity is freshman's dream. (ii) Tr is additive and F_2-linear; it is a nonzero
polynomial of degree 2^(k-1) < 2^k, hence a nonzero function. The map
mu |-> (a |-> Tr(mu*a)) from K to the F_2-dual of K is linear and injective
(if mu != 0 and Tr(mu*a) = 0 for all a, then Tr vanishes on the image of the
bijection a |-> mu*a, i.e. on all of K, a contradiction); dimensions agree
(both k), so it is a bijection. qed

## 2. Subspace polynomials

**Lemma 2.1.** Let L = b + V be an affine subspace of K, dim V = D. Then
P_L(X) := prod_{alpha in L}(X - alpha) is a monic affine linearized
polynomial of degree 2^D. P_L is linearized iff 0 in L iff L is linear.
Moreover the root set of P_L in K is exactly L.

*Proof.* Induction on D for V linear. D = 0: P = X. D >= 1: write
V = V' + F_2*v, v in V\V', Q := P_{V'} (linearized, monic, degree
2^(D-1) by induction), c := Q(v). Then Q(X)*(Q(X)+c) is monic of degree 2^D;
its zero set is {X : Q(X) = 0} union {X : Q(X) = Q(v)} = V' union {X :
Q(X+v) = 0} (Q linearized) = V' union (v+V') = V; hence P_V = Q^2 + cQ,
which is linearized since Q^2 = sum_j a_j^2 X^(2^(j+1)) is. For L = b + V:
P_L(X) = P_V(X - b) = P_V(X) + P_V(b), affine linearized of degree 2^D. The
constant term P_V(b) is the product of all roots up to sign, nonzero iff
b not in V. The root set statement holds by definition. qed

**Lemma 2.2.** Let L = b + V, dim V = D >= 1, P = P_L.

1. P(K) is an affine subspace L_0 of dimension k - D.
2. Q_L := P_{L_0} is a monic affine linearized polynomial of degree 2^(k-D)
   whose root set in K is exactly L_0 (an affine subspace of dimension
   k - D, so Q_L is the kernel-subspace polynomial of its root set), and
   Q_L(P_L(X)) = X^(2^k) - X as polynomials.
3. If moreover L is linear, then P_L is linearized, P_L(Q_L(Y)) =
   Y^(2^k) - Y as polynomials, and Q_L(K) = L.

*Proof.* (1) Write P = P̂ + c with P̂ = P - P(0) linearized. For any fixed
x_0 in L and all x in K: P(x) = 0 iff P̂(x) = c = P̂(x_0) iff
x in x_0 + ker(P̂). So L is a coset of ker(P̂), whence |ker P̂| = 2^D and
|P̂(K)| = 2^k / 2^D; thus P(K) = c + P̂(K) has F_2-dimension k - D.

(2) Q_L(P_L(X)) is monic of degree 2^(k-D) * 2^D = 2^k and vanishes at
every X in K, because P_L(X) lies in L_0, which is the root set of
Q_L = P_{L_0} by definition. The unique monic polynomial of degree 2^k
vanishing on all of K is X^(2^k) - X (any polynomial vanishing on K is
divisible by X^(2^k) - X), so Q_L(P_L(X)) = X^(2^k) - X. This identity
also shows that every element of L_0 = P_L(K) is a root of Q_L, so the
root set of Q_L contains L_0 and, having at most deg Q_L = 2^(k-D)
elements while |L_0| = 2^(k-D), equals L_0.

(3) Now let L be linear. Then P_L(0) = 0 (0 is a root), so c = 0, P_L is
linearized, and P_L(X^(2^k)) = (P_L(X))^(2^k) (coefficients a of P_L lie
in K, a^(2^k) = a). Let g(Y) := P_L(Q_L(Y)) - Y^(2^k) + Y. Using the part
(2) identity and linearity:

    g(P_L(X)) = P_L(Q_L(P_L(X))) - (P_L(X))^(2^k) + P_L(X)
              = P_L(X^(2^k) - X) - (P_L(X))^(2^k) + P_L(X)
              = P_L(X^(2^k)) - P_L(X) - (P_L(X))^(2^k) + P_L(X) = 0.

If g != 0, then g(P_L(X)) is a nonzero polynomial in X (a nonzero
constant when deg g = 0, and of degree deg(g)*2^D >= 1 otherwise), while
the computation above shows it is identically 0 — a contradiction; so
g = 0 and P_L(Q_L(Y)) = Y^(2^k) - Y.
For y in K: P_L(Q_L(y)) = 0, i.e. Q_L(K) is contained in the root set of
P_L, which is L. Writing Q_L = Q̂ + c', the solutions of Q_L(y) = 0 form a
coset of ker(Q̂) with at most 2^(k-D) elements, so |ker Q̂| <= 2^(k-D) and
|Q_L(K)| = 2^k/|ker Q̂| >= 2^D = |L|; hence Q_L(K) = L. qed

## 3. The structural lemma

Fix an algebraic closure F̄ of F_2. The nonzero elements of F_(2^u) are
exactly the roots of X^(2^u - 1) - 1; hence for a, b, c in F̄ with a, c != 0,
b^(2^u - 1) = c and a^(2^u) = c*a together give (a/b)^(2^u - 1) = 1, i.e.
a/b in F_(2^u)*, i.e. a in b*F_(2^u). (Such b exists in F̄ since X^(2^u-1)-1
splits over F̄.) Also F_(2^u) intersect K = F_(2^gcd(u,k)).

**Claim 3.1.** For Q(X) = sum_j a_j X^(2^j) + a-hat in K[X] and u >= 0,
Q(X)^(2^u) mod (X^(2^k) - X) = sum_j a_j^(2^u) X^(2^((j+u) mod k))
+ a-hat^(2^u).

*Proof.* Freshman's dream u times, then X^(2^(j+u)) = X^(2^((j+u) mod k))
modulo the ideal, by induction on j + u - k using X^(2^k) = X. qed

**Lemma 3.2 (structural lemma).** Let A be an affine subspace of K of
dimension D >= 1, Q_A(X) = sum_{j=0}^{k-1} a_j X^(2^j) + a-hat its
image-subspace polynomial (monic affine linearized of degree 2^(k-D); so
a_j = 0 for k - D < j <= k - 1, and by Lemma 2.2(2) Q_A has exactly
2^(k-D) roots in K, namely L_0 = P_A(K)).

1. Among any D consecutive cyclic indices J = {(r+j) mod k : 0 <= j <= D-1}
   some a_j is nonzero. Consequently a_(k-D) != 0 (the window starting at
   k - D ends in the D - 1 zero positions) and a_0 != 0 (the window starting
   at k - D + 1 covers {k-D+1,...,k-1,0}).
2. Suppose A is a linear subspace (b = 0) of dimension at least 2 that is
   not contained in a constant multiple of a proper subfield of K (for k
   prime: any linear subspace of dimension >= 2, since constant multiples
   of F_2 have at most 2 elements). If r != (k - D + 1) mod k and
   a_j = 0 for all j in J = {(r+j) mod k : 0 <= j <= D-2}, then the
   direction space of A is contained in a constant multiple of F_2, so
   dim A <= 1 — a contradiction; hence some a_j with j in J is nonzero.

*Proof.* (1) Assume all a_j, j in J, vanish; set u := (k - r - D) mod k. By
Claim 3.1 the coefficient of X^(2^j) in Q'(X) := Q_A(X)^(2^u)
mod (X^(2^k)-X) is a_((j-u) mod k)^(2^u). For j in {k-D,...,k-1} the index
(j - u) mod k = (j + r + D) mod k, and as j runs through the D consecutive
integers k-D..k-1 the residues form exactly J. Hence all coefficients of Q'
at positions k-D..k-1 vanish and deg Q' <= 2^(k-D-1). Q' is nonzero (the
coefficient positions (j-u) mod k permute {0..k-1}, so some nonzero a_j
appears). Every root alpha of Q_A satisfies Q'(alpha) = Q_A(alpha)^(2^u) = 0
(evaluation and reduction agree on K, as in Proposition 1.1). A nonzero
polynomial of degree 2^(k-D-1) cannot have 2^(k-D) roots. Contradiction.

(2) Assume all a_j, j in J, vanish; set u := (k - r - D + 1) mod k. Note
u != 0 because r != k - D + 1 mod k. As in (1), with the D - 1 positions
j in {k-D+1,...,k-1} mapping to J, all coefficients of Q'(X) :=
Q_A(X)^(2^u) mod (X^(2^k)-X) above position k - D vanish; Q' is nonzero, so
it has degree exactly 2^(k-D) (degree < 2^(k-D) would again contradict the
2^(k-D) roots) and, having the same roots as Q_A with monic Q_A of equal
degree, Q' = c*Q_A where c in K* is the top coefficient of Q'. Evaluating on
K: Q_A(X)^(2^u) = c*Q_A(X) for all X in K. Since A is linear, Q_A(K) = A by
Lemma 2.2(3). Pick b in F̄ with b^(2^u-1) = c; for every nonzero v = Q_A(X)
in A, (v/b)^(2^u - 1) = 1, so v = b*e with e in F_(2^u)*. Fix
v_1 != v_2 in A with v_1 != 0 (possible: dim A >= 2 gives |A| >= 4 and a
nonzero value exists). For any w in A, w and v_2 are b*e_w, b*e_2 with
e_w, e_2 in F_(2^u) union {0}, so (w - v_2)/(v_1 - v_2) =
(e_w - e_2)/(e_1 - e_2) lies in F_(2^u); the left side also lies in K
(all values in K and v_1 != v_2), hence in
F_(2^u) intersect K = F_(2^gcd(u,k)). Since 1 <= u <= k-1 and k is prime,
gcd(u,k) = 1, so every difference w - v_2 lies in (v_1 - v_2)*F_2: the
direction space of A is contained in a constant multiple of F_2 and
dim A <= 1, contradicting dim A >= 2. qed

## 4. The permanent formula

Throughout Sections 4-6, "mod ideal" abbreviates reduction modulo the ideal
(X_i^(2^k) - X_i : i) of Proposition 1.1, always followed by taking the
unique representative of degree < 2^k in each variable.

For linearized Q(X) = sum_j a_j X^(2^j) define the associated matrix
M(i,j) = (a_(j-i mod k))^(2^i) in K (rows/columns indexed 0..k-1) and the
0-1 indicator M'(i,j) = [M(i,j) != 0]. M[A,B] is the submatrix with row set
A and column set B; Perm is the permanent; B + s = {(b+s) mod k}.

**Claim 4.1.** M(i,j) is the coefficient of X^(2^j) in Q(X)^(2^i)
mod (X^(2^k)-X). *Proof.* Claim 3.1 with u = i. qed

**Lemma 4.2.** Let mu in K*; A_1..A_s, B_1..B_s subset {0..k-1} with
|A_i| = |B_i| in [1, k-1]; alpha_i = sum_{a in A_i} 2^a, beta_i = sum_{b in B_i} 2^b;
Q_i linearized with associated matrices M_i. The coefficient of
M := prod_i X_i^(beta_i) in R := Tr(mu prod_i Q_i(X_i)^(alpha_i))
mod (X_i^(2^k)-X_i)_i equals

    sum_{w=0}^{k-1} mu^(2^w) * prod_i Perm(M_i[A_i, B_i - w])^(2^w).

*Proof.* Expand Tr: R = sum_w mu^(2^w) prod_i Q_i^(alpha_i 2^w) mod ideal,
and Q_i^(alpha_i 2^w) = prod_{a in A_i} Q_i^(2^(a+w)) = prod_{a in A_i+w}
(sum_j M_i(a,j) X_i^(2^j)) (freshman's dream, Claim 4.1). Choosing for each
a in A_i + w an exponent h(a) in {0,...,k-1} (all exponents of the already
reduced factors Q_i^(2^u) are binary powers 2^j with j < k, by Claim 4.1),
the resulting X_i-exponent before reduction is
E_i = sum_{a in A_i+w} 2^(h(a)), a sum of |A_i + w| = |B_i| <= k - 1 binary
powers with exponents below k. A selection contributes to the coefficient
of X_i^(beta_i) only if beta_i = E_i mod (2^k - 1) (for E >= 1, reduction
sends X^E to X^(E mod (2^k - 1)) with the convention that residue 0 gives
exponent 2^k - 1; here beta_i <= 2^k - 2 because |B_i| = e_i + 1 <= k - 1,
so the residue-0 convention never applies to the target). Now
wt(beta_i) = |B_i|
(beta_i is a sum of |B_i| distinct binary powers), and Lemma 1.3 applied to
E_i, q = beta_i, |S| = |B_i| gives wt(beta_i) <= |B_i| with equality only
in the distinct case; so equality must hold, Lemma 1.3(2) rules out
repeated exponents among the h(a), and then Lemma 1.3(1) (with
|S| = e_i + 1 <= k - 1) gives E_i = beta_i exactly. By uniqueness
of binary representation, h is a bijection from A_i + w onto B_i, and the
product of chosen entries is
prod_{a in A_i + w} M_i(a, h(a)) = Perm(M_i[A_i + w, B_i]). The identity
Perm(M[A+w, B]) = Perm(M[A, B-w])^(2^w) follows since the (a,b)-entry of the
left minor is a_(b-a-w mod k)^(2^(a+w)) while the (a,b+w)-entry of
M[A, B-w], raised to the 2^w, gives the same value; the two minors differ by
the column relabeling b <-> b + w. Summing over w proves the formula. qed

**Remark 4.3.** If (a) each M'_i[A_i,B_i] is triangular (after reordering
rows and columns together) with all diagonal entries 1, and (b) for each
w in {1,...,k-1} some M'_i[A_i, B_i - w] has an all-zero column, then the
coefficient of M is nonzero: the w = 0 summand is a product of nonzero
permanents (a triangular 0-1 matrix with unit diagonal has permanent 1),
while a zero column forces permanent 0. qed

## 5. Independent sources

**Theorem 5.1.** Let A_1,...,A_s (s >= 1) be LINEAR subspaces of K of
dimensions D_1 >= ... >= D_s >= 2 (so ordered). Let
Q_i = sum_j a_(i,j) X^(2^j) + a-hat_i be the image-subspace polynomial of
A_i, mu in K*, and integers 1 <= e_i <= D_i - 2; alpha_i := 2^(e_i+1) - 1.
If

    sum_i (D_i - (e_i + 1)) > k - D_1 + 1,

then R(X) := Tr(mu prod_i Q_i(X_i)^(alpha_i)) mod (X_i^(2^k)-X_i)_i has a
monomial prod_i X_i^(beta_i) with wt(beta_i) = e_i + 1 for every i and
nonzero coefficient.

*Proof.* Dropping constants: replacing Q_i by its linearized part
Q°_i = sum_j a_(i,j) X^(2^j) changes no coefficient of any monomial whose
F_2-degree in every variable X_i is e_i + 1. Fix a trace power 2^w and
expand (Q°_i + a-hat_i)^(alpha_i 2^w) =
prod_{l=0}^{e_i}(Q°_i^(2^(l+w)) + a-hat_i^(2^(l+w))): a selection choosing
the constant for at least one l produces, before reduction, a product in
which X_i occurs as a sum of at most e_i binary powers (one per remaining
factor), so its X_i-exponent E_i is a sum of at most e_i powers and, by
Lemma 1.3, its reduced weight is at most e_i (and by Lemma 1.3(3) the
q = 0 corner cannot arise, as e_i + 1 <= k - 1 summands cannot represent
0 modulo 2^k - 1). Our target
monomials have X_i-exponent of weight exactly e_i + 1 in every variable,
and before reduction each variable X_i appears only in the i-th factor; so
by Lemma 1.3 (weight inequality), only the all-Q° selection can contribute
to them. Hence we may work with Q°_i and its coefficients a_(i,j).

*Index sets.* Let l_1 := 0 and l_{i+1} := (l_i + D_i - (e_i+1)) mod k.
Note D_i - (e_i+1) >= 1. Define B_i (an e_i + 1 element subset of
{0,...,k-1}):

- If a_(i,l_i) != 0:  B_i := {l_i} union {k-D_i+1, ..., k-D_i+e_i}.
- If a_(i,l_i) = 0: let jb_i := max{j < l_i : a_(i,j) != 0} (exists since
  a_(i,0) != 0, in the non-wrapping sense: for the indices i used below we
  will have l_i >= 1) and j_i := min{j > l_i : a_(i,j) != 0} (exists since
  a_(i,k-D_i) != 0; for the indices used below, no wrap occurs). Let
  q_i := j_i - jb_i - 1 (the length of the maximal run of zero coefficients
  jb_i+1, ..., j_i-1 containing l_i), s_i := min(q_i, e_i), and

      B_i := {j_i} union {jb_i+1, ..., jb_i+s_i} union
             {k-D_i+s_i+1, ..., k-D_i+e_i}.

The displayed sets are pairwise disjoint (run positions have zero
coefficients, j_i a nonzero one, tail positions lie above k-D_i while the run
lies strictly below j_i <= k-D_i... the last containment is proved where it
is used below), so |B_i| = 1 + s_i + (e_i - s_i) = e_i + 1.

**Claim 5.2.** For every i, M'_i[{0..e_i}, B_i] is lower triangular with all
diagonal entries 1 (columns ordered as displayed).

*Proof.* The entry at row a in {0..e_i}, column c, is
[a_(i,(c-a) mod k) != 0].

Case a_(i,l_i) != 0, columns c_0 = l_i, c_b = k-D_i+b (b = 1..e_i).
Entries above the diagonal: (a, c_b) with a < b has index
(k - D_i + b - a) mod k with 1 <= b - a <= e_i <= D_i - 2, an integer in
{k-D_i+1, ..., k-1} where all coefficients vanish. Diagonal: (0, c_0) = 1 by
assumption; (b, c_b) = [a_(i,k-D_i) != 0] = 1.

Case a_(i,l_i) = 0, columns c_0 = j_i, c_b = jb_i+b (b = 1..s_i), then
c_b = k-D_i+b (b = s_i+1..e_i). (0, c_0) = 1. (a, c_0) for 1 <= a <= s_i:
index (j_i - a) mod k lies in the zero run {jb_i+1, ..., j_i-1}: 0.
(a, c_b), 1 <= b <= s_i, b > a (above diagonal): index jb_i + (b-a) with
1 <= b-a <= s_i <= q_i lies in the run: 0. Diagonal (b, c_b) = [a_(i,jb_i)] =
1. Tail columns: (a, c_b) with s_i < b <= e_i, a < b: index k-D_i+(b-a),
1 <= b-a <= e_i - s_i <= e_i <= D_i-2, again in the vanishing range: 0;
diagonal (b, c_b) = [a_(i,k-D_i)] = 1. Entries below the diagonal are
unconstrained. qed

**Claim 5.3.** (i) If c in B_i and 1 <= v <= D_i - (e_i+1), then the column
c - w with w = (c + v) mod k of M'_i[{0..e_i}, B_i - w] is all-zero.
(ii) For every w in {1,...,k-1} there is i such that M'_i[{0..e_i}, B_i - w]
has an all-zero column.

*Proof.* (i) The column c - w has entries a_(i,(c-w-a) mod k), a = 0..e_i.
With w = c + v: (c - w - a) mod k = (v + a) mod k (sign), more precisely
c - w - a = -(v + a) with 1 <= v + a <= (D_i - e_i - 1) + e_i = D_i - 1, so
the index is k - (v+a) in {k-D_i+1, ..., k-1}: all zero.

(ii) Let S(i) := sum_{i' <= i}(D_{i'} - (e_{i'}+1)) and let t0 be minimal
with S(t0) >= k - D_1 + 1; it exists since S(s) > k - D_1 + 1. Each summand
lies in [1, D_1 - 2], so S(t0) < (k - D_1 + 1) + D_1 - 1 = k. For
1 <= i <= t0, l_i = S(i-1) (no reduction), and for 2 <= i <= t0 in
particular l_i >= S(1) = D_1 - (e_1+1) >= 1 and
l_i <= S(t0-1) <= k - D_1 <= k - D_i (using D_i <= D_1). Finally
l_{t0+1} = S(t0) lies in [k - D_1 + 1, k - 1] (used only for the final
concatenation below).

The case a_(i,0) = 0 never occurs at i = 1 (a_(1,0) != 0 by Lemma 3.2(1)),
so B_1 is of the first form: B_1 = {0} union {k-D_1+1, ..., k-D_1+e_1}.
For 2 <= i <= t0, jb_i exists as an ordinary maximum (l_i >= 1 and
a_(i,0) != 0); and if l_i >= k - D_i then, being <= k - D_1 <= k - D_i,
equality l_i = k - D_i holds and a_(i,l_i) = a_(i,k-D_i) != 0 — first case;
hence in the second case l_i < k - D_i and j_i exists as an ordinary minimum
with j_i <= k - D_i, so the run {jb_i+1..j_i-1} and j_i lie below or at
k - D_i, while the tail {k-D_i+s_i+1, ..., k-D_i+e_i} lies strictly above
k-D_i and below k: the three parts of B_i are disjoint as promised.
For indices i > t0 (whose B_i never enter the coverage below), if the
ordinary extrema jb_i, j_i do not exist or a_(i,l_i) != 0 fails in a way
that leaves the second case undefined, set instead

    B_i := {0} union {k-D_i+1, ..., k-D_i+e_i},

the first-case shape at position 0, which is legitimate because a_(i,0) != 0
and a_(i,k-D_i) != 0 always; Claim 5.2's argument for the first case applies
verbatim with l_i replaced by 0. This makes B_i well defined for every
i in [s], as required by the w = 0 summand of Lemma 4.2.

Coverage. By (i), a column c in B_i covers all w in the (ordinary) interval
[c+1, c+D_i-(e_i+1)]:

- i = 1: column 0 covers [1, D_1-(e_1+1)] = [l_1+1, l_2]; each tail column
  k-D_1+b (b = 1..e_1) covers [k-D_1+b+1, k+b-e_1-1], and the union over b
  is [k-D_1+2, k-1]: the left end is at b = 1, the right end k-1 at
  b = e_1, and consecutive intervals meet or overlap since
  D_1 - e_1 - 1 >= 1.
- 2 <= i <= t0, first case (a_(i,l_i) != 0): column l_i covers
  [l_i+1, l_i+D_i-(e_i+1)] = [l_i+1, l_{i+1}].
- second case, q_i < e_i (s_i = q_i): {j_i} union {jb_i+1..jb_i+q_i} =
  {jb_i+1, ..., j_i}, which contains l_i; so l_i in B_i and the previous
  bullet applies.
- second case, q_i >= e_i (s_i = e_i): first, q_i <= D_i - 2. For if
  q_i >= D_i - 1, the D_i - 1 consecutive indices
  {jb_i+1, ..., jb_i+D_i-1} would all lie in the run and vanish; this
  window starts at r' = jb_i+1 with r' <= l_i <= k - D_i < k - D_i + 1,
  so r' != k - D_i + 1 and Lemma 3.2(2) applies (A_i is linear of
  dimension D_i >= 2), giving a nonzero coefficient inside the window —
  contradiction. So j_i <= jb_i + D_i - 1. The run columns jb_i+b (b=1..e_i)
  cover
  [jb_i+2, jb_i+D_i-1] (consecutive intervals [jb_i+b+1, jb_i+b+D_i-e_i-1],
  b = 1..e_i, each of length D_i-e_i-1 >= 1, with overlaps), which contains
  [jb_i+2, j_i]; the column j_i covers [j_i+1, j_i+D_i-e_i-1]. Together they
  cover [jb_i+2, j_i+D_i-e_i-1], which contains [l_i+1, l_{i+1}] because
  jb_i+2 <= l_i+1 (jb_i < l_i) and j_i+D_i-e_i-1 >= l_i+D_i-e_i-1 =
  l_{i+1} (j_i > l_i).

All intervals being ordinary (their right ends are at most
l_{t0+1} <= k-1 resp. k-1), the columns of B_1,...,B_(t0) cover

    [1, l_2] union [l_2+1, l_3] union ... union [l_(t0)+1, l_(t0+1)]
        = [1, S(t0)],   together with   [k-D_1+2, k-1],

and since S(t0) >= k - D_1 + 1 this is all of {1, ..., k-1}. qed

By Claims 5.2, 5.3 and Remark 4.3 (with A_i = {0,...,e_i}, |B_i| = e_i+1),
the monomial prod_i X_i^(beta_i), beta_i := sum_{b in B_i} 2^b, has nonzero
coefficient in R. This proves Theorem 5.1. qed

## 6. Proof of Theorem A

**Lemma 6.1 (decomposition).** Let A be an affine subspace of K^r and
gamma in A. There are F_2-linear subspaces Y_1..Y_r of K and F_2-linear
maps sigma_ij: Y_j -> K (j < i) with

    A = {(x_1..x_r) : x_i = gamma_i + y_i + sum_{j<i} sigma_ij(y_j), y_i in Y_i}

and dim A = sum_i dim Y_i.

*Proof.* Let V (dim D) be the direction of A with basis v_1..v_D; form the
(kr) x D matrix M over F_2 with these columns and column-reduce it to
echelon form M~ = MT (T invertible; column space = V unchanged): column c
has pivot row j_c, entries above the pivot zero, j_1 < ... < j_D, pivot
entries 1. For y in F_2^D let y_i in K be the i-th block of
sum_{c : i(c)=i} y_c (column c), where i(c) is the block containing row j_c;
Y_i := the set of such y_i. Row j of block i of M~ y: pivot rows of block i
carry the coordinates of y_i; a non-pivot row j of block i receives, by
echelon structure (zeros above pivots), contributions only from columns with
j_c <= j, hence with i(c) <= i; columns of block i are absorbed in y_i, and
a column c with i(c) = i' < i contributes a fixed F_2-linear function of the
coordinate y_c of y_{i'}; call it sigma_{i,i'}(y_{i'}). So M~ y = x with
x_i = y_i + sum_{i'<i} sigma_{i,i'}(y_{i'}), and A = gamma + V gives the
displayed form. The map y |-> (y_1..y_r) is F_2-linear, injective (within each block,
the pivot coordinates of y_i can be read off the vector y_i by downward
induction on the pivot rows, since the echelon form is unitriangular on
the pivot rows of each block) and surjective onto prod Y_i, so
D = sum dim Y_i. qed

**Lemma 6.2 (perturbation invariance).** Let I_0 = {i_1 < ... < i_s}
subset [r]. For i in I_0 let P_i(X_i) be linearized and P̂_i(X_1..X_{i-1}) :=
sum_{j<i} P_ij(X_j) + gamma_i with P_ij linearized, gamma_i in K; e_i >= 1,
alpha_i = 2^(e_i+1)-1, mu in K*. Then for beta_i with wt(beta_i) = e_i+1,
the coefficient of prod_{i in I_0} X_i^(beta_i) in
g := Tr(mu prod_{i in I_0}(P_i(X_i)+P̂_i(X_{<i}))^(alpha_i)) mod ideal
equals the coefficient in g' := Tr(mu prod_{i in I_0} P_i(X_i)^(alpha_i))
mod ideal.

*Proof.* Fix a trace power 2^w and expand
(P_i + P̂_i)^(alpha_i 2^w) = prod_{l=0}^{e_i}(P_i^(2^(l+w)) + P̂_i^(2^(l+w))).
A selection choosing P̂_i^(2^(l+w)) for some l uses the constant-or-
other-variables factor in that slot; the all-P selection is
P_i^(alpha_i 2^w). Process i = i_s first: the variable X_{i_s} occurs in no
other factor (P-parts involve only X_i, i < i_s; P̂-parts only X_{<i},
i < i_s), so in a selection with a P̂-choice at i_s, the X_{i_s}-exponent
before reduction is a sum of at most e_s binary powers 2^h (reducing each
factor P_i^(2^u) or P̂_i-part modulo the ideal first — which does not change
the fully reduced product — turns its every term into c*X_i^(2^h) or a
term c*X_j^(2^h) in a lower-indexed variable X_j (j < i), with h < k, by
X^(2^k) = X; one such power per remaining factor); by the weight
inequality (Lemma 1.3), its reduced weight is at most e_s (and the q = 0
corner is excluded by Lemma 1.3(3)), while
wt(beta_{i_s}) = e_s + 1, so such selections cannot contribute to the
coefficient of our monomial. Having fixed the pure part at i_s (which
involves X_{i_s} only), the same reasoning applies to i_{s-1}, ..., i_1
(the X_{i_v}-exponent can come only from the i_v factor, whose P̂-choices
cap the pre-reduction exponent at e_v powers). Only the all-pure selections
contribute, on every trace power, giving identical contributions in g and
g'. qed

**Theorem 6.3.** Let 2 <= t < r with k >= r(r+1)/t; e_i := r+1-i and
alpha_i := 2^(e_i+1)-1 for i in [r]; mu in K*;

    h(x_1..x_r) := Tr( mu * sum_{I subset [r], |I|=t} prod_{i in I} x_i^(alpha_i) ).

Then h is nonconstant on every affine subspace A of K^r with
dim A > rk/t + kt + r(r+1).

*Proof.* Let D = dim A, gamma in A, and decompose as in Lemma 6.1 with
spaces Y_i, dims D_i. Let Q_i be the image-subspace polynomial of Y_i (set
Q_i := 0 if D_i = 0) and Q_i = Q°_i + a-hat_i its split into linearized part
and constant. For j < i the map sigma_ij ∘ Q°_j: K -> K is F_2-linear, hence
(by Proposition 1.2(i)) equal to a unique linearized Q_ij. Set

    R_i(X_1..X_i) := Q_i(X_i) + sum_{j<i} Q_ij(X_j) + gamma_i.

By Lemma 6.1, R := (R_1..R_r): K^r -> K^r maps onto A, and h ∘ R is
represented (Proposition 1.1) by

    H := Tr( mu * sum_{|I|=t} prod_{i in I} R_i^(alpha_i) ) mod ideal.

Let S := {i : D_i > r+1}. If |S| <= t-1 then
D <= |S| k + (r-|S|)(r+1) <= (t-1)k + (r-t+1)(r+1), and
(t-1)k + (r-t+1)(r+1) = tk + r(r+1) - (k + (t-1)(r+1)) < tk + r(r+1)
< rk/t + kt + r(r+1) < D,
a contradiction; hence |S| >= t. Also every i
outside S has D_i <= r+1, so sum_{i in S} D_i >= D - r(r+1) > rk/t + kt.

Pick I_0 subset S with |I_0| = t consisting of t indices with the largest
D_i; the sum of the t largest among values in S (|S| <= r) is at least
t/r of the total: sum_{i in I_0} D_i >= (t/r)(rk/t + kt) = k + kt^2/r.
Each i in I_0 has D_i >= r+2, so e_i = r+1-i satisfies 1 <= e_i <= r <=
D_i - 2; Y_i is a linear space of dimension >= r+2 >= 4, as required by
Theorem 5.1 (linear of dimension >= 2). Relabel I_0 = {i'_1 > ...} so
that dimensions are non-increasing (keeping (D, e) pairs attached) and apply
Theorem 5.1 to these t spaces with mu: its numerical hypothesis holds since

    sum_{i in I_0}(D_i - (e_i+1)) = sum_{i in I_0}(D_i - (r+1)) +
    sum_{i in I_0}(i-1) >= (k + kt^2/r) - t(r+1) + 0 >= k > k - max D_i + 1,

using kt^2/r >= t(r+1), i.e. k >= r(r+1)/t, and max D_i >= 2. So there is a
monomial M* = prod_{i in I_0} X_i^(beta_i), wt(beta_i) = e_i+1, with nonzero
coefficient in g'(X) := Tr(mu prod_{i in I_0} Q°_i(X_i)^(alpha_i)) mod ideal
(we applied 5.1 to the linearized parts; by the opening paragraph of its
proof, the same coefficient is obtained from the affine Q_i, and by Lemma 6.2
— with P_i := Q°_i, P̂_i := the remaining terms of R_i, whose (D_i, e_i)
structure is exactly as required — the coefficient of M* in the summand
g := Tr(mu prod_{i in I_0} R_i^(alpha_i)) of H is nonzero).

Uniqueness of the summand. Let I != I_0, |I| = t, fix a trace power 2^w and
expand the inner product prod_{i in I} prod_{l=0}^{e_i} R_i^(2^(l+w)):
reducing each factor R_i^(2^u) modulo the ideal first (this does not change
the fully reduced product, and X^(2^k) = X brings every exponent below
2^k), each factor is a sum of terms c*X_j^(2^v) with v < k, j <= i, and a
constant. Selecting one term per factor: each nonconstant selected term is
c*X_j^(2^v) for exactly one variable X_j, contributing exactly one binary
power to that variable's exponent, so the X_j-exponent before reduction is
a sum of 2^v (v < k), one per selected nonconstant term reading X_j. By the
weight
inequality (Lemma 1.3), the reduced weight of the X_j-exponent is at most
the number of nonconstant selected terms reading X_j (and the q = 0 corner
cannot produce beta_j != 0 by Lemma 1.3(3)). If the reduced product
contains M*, then for every
j in I_0 the number lambda(j) := sum_{i in I} lambda(i,j) of selected terms
reading X_j is at least wt(beta_j) = e_j + 1 = r+2-j, where lambda(i,j) :=
number of selected terms from factor i reading X_j; note lambda(i,j) = 0 for
j > i and sum_{j} lambda(i,j) <= e_i + 1 = r+2-i.

*Claim 6.4.* If I, I_0 are t-element subsets of [r] and lambda(i,j) >= 0
integers satisfy lambda(i,j) = 0 for j > i, sum_j lambda(i,j) <= r+2-i, and
sum_{i in I} lambda(i,j) >= r+2-j for all j in I_0, then I = I_0.

*Proof.* Sum the demands over I_0 and the capacities over all j:
sum_{I_0}(r+2-j) <= sum_I (r+2-i), so with |I_0| = |I| = t,
(*) sum_{i in I} i <= sum_{j in I_0} j. For every p >= 1, weighting the
demands by j^(p-1) and using j <= i on the support of lambda(i, j):

    sum_{j in I_0}(r+2-j) j^(p-1) <= sum_{i in I} sum_j lambda(i,j) i^(p-1)
    <= sum_{i in I}(r+2-i) i^(p-1).                     (**)

Expanding (**) gives sum_{I_0} j^p - sum_I i^p >= (r+2)( sum_{I_0} j^(p-1)
- sum_I i^(p-1) ); at p = 1 the right side is 0, so sum_{I_0} j >=
sum_I i. Inductively, if sum_{I_0} j^q >= sum_I i^q for all q < p, then
the expanded (**) at p gives sum_{I_0} j^p - sum_I i^p >= (r+2)(nonnegative)
>= 0. So sum_{I_0} j^p >= sum_I i^p for all p >= 1.

Next, max I_0 <= max I: if j := max I_0 > max I, the demand at j is
r+2-j > 0 but lambda(i, j) = 0 for every supplier i in I (all satisfy
i <= max I < j), so the demand at j is unmet. Conversely max I_0 >= max I:
if max I_0 < max I, then for large p the left power sum is at most
t (max I_0)^p < (max I)^p <= sum_I i^p, contradicting the display. Hence
max I_0 = max I =: m_1, and m_1 lies in both sets. The only supplier able
to serve X_(m_1) is i = m_1 (suppliers satisfy i >= m_1), so its full
capacity r+2-m_1 goes to m_1 and lambda(m_1, j) = 0 for j < m_1. Delete
m_1 from both sets: the demands at the remaining j in I_0 must be met by
I \ {m_1} with unchanged capacities, and the power-sum inequalities
persist for the remaining sets (subtract m_1^p from both sides), so the
same two steps apply again. Inducting over the sorted elements, I_0 = I.
qed

By Claim 6.4 (whose hypotheses are exactly the counting relations derived
above, applied once per trace power 2^w), no summand I != I_0 produces M*,
and the I_0-summand produces it with nonzero coefficient. So H has a
positive-degree monomial with nonzero coefficient, and h is nonconstant on A
by Proposition 1.1. qed

**Proof of Theorem A.** Let n >= 10^20 and u := n^(1/5) >= 10^4. Set
c := k/u^3, so by `field.md` (Bertrand at A = ceil(u^3 + 2u^2)):

    1 + 2/u < c < 2 + 4/u + 2/u^3.                          (P)

In particular 1 < c <= 2.001 (using 4/u + 2/u^3 <= 0.0005 < 0.001 at
u >= 10^4), r = ceil(n/k) = ceil(u^2/c), and t = ceil(sqrt(r)).

*Parameter conditions.* From r = ceil(u^2/c): u^2/c <= r <= u^2/c + 1.

(i) t >= 2: r >= u^2/c >= u^2/2.001 >= 10^8/2.001 > 4*10^7, so
    t = ceil(sqrt(r)) > sqrt(4*10^7) > 6000 >= 2.

(ii) t < r: t <= sqrt(r) + 1 <= sqrt(u^2/c + 1) + 1 <= sqrt(u^2 + 1) + 1
    <= u + 2 (c >= 1, sqrt(a+1) <= sqrt(a)+1 <= u+1), and
    u + 2 < u^2/2.001 <= r for u >= 5 (5 + 2 = 7 < 25/2.001 = 12.49).

(iii) k >= r(r+1)/t: since t >= sqrt(r),
     r(r+1)/t <= (r+1)*sqrt(r) <= (u^2/c + 2) * (u + 1),
 using sqrt(r) <= sqrt(u^2/c + 1) <= u + 1. Expanding:
 (u^2/c + 2)(u + 1) = u^3/c + u^2/c + 2u + 2. It suffices that
     c*u^3 >= u^3/c + u^2/c + 2u + 2,
 which, after multiplying by c > 0, reads
     c^2*u^3 >= u^3 + u^2 + (2u + 2)*c.                    (*)
 Now c > 1 + 2/u, so c^2 > (1 + 2/u)^2 = 1 + 4/u + 4/u^2, whence
 c^2*u^3 > u^3 + 4u^2 + 4u; and (2u+2)*c <= (2u+2)*2.001 =
 4.002u + 4.002 (c <= 2.001), so the right side of (*) is at most
 u^3 + u^2 + 4.002u + 4.002. Since
 4u^2 + 4u >= u^2 + 4.002u + 4.002  iff  3u^2 - 0.002u - 4.002 >= 0,
 which holds for every u >= 2 (3*4 - 0.004 - 4.002 > 0, and 3u^2 - 0.002u
 is increasing), condition (*) holds, hence k >= r(r+1)/t.

*Dimension bound.* We bound rk/t + kt + r(r+1) separately, in terms of c.

    rk/t <= (n+k)/sqrt(n/k)            (r <= n/k+1, t >= sqrt(n/k))
          = (u^5 + c*u^3)*sqrt(c)/u = sqrt(c)*u^4 + c^(3/2)*u^2,
    kt   <= k*(sqrt(n/k) + 2)          (t <= sqrt(r)+1 <= sqrt(n/k)+2)
          = c*u^3*(u/sqrt(c) + 2) = sqrt(c)*u^4 + 2c*u^3,
    r(r+1) <= (u^2/c + 1)(u^2/c + 2) = u^4/c^2 + 3u^2/c + 2.

Summing:

    rk/t + kt + r(r+1) <= h(c)*u^4 + 2c*u^3 + (c^(3/2) + 3/c)*u^2 + 2,
    h(c) := 2*sqrt(c) + 1/c^2.

*The maximum of h on the interval (1, 2.001].*
h'(c) = 1/sqrt(c) - 2/c^3 = (c^(5/2) - 2)/c^3, so h is strictly decreasing
on (1, 2^(2/5)) and strictly increasing on (2^(2/5), 2.001]. Hence h attains
its maximum on the closed interval [1, 2.001] at an endpoint:
h(1) = 2 + 1 = 3, and h(2.001) = 2*sqrt(2.001) + 1/2.001^2 < 2*1.4146 +
1/4 = 3.0792 (1.4146^2 = 2.00109... > 2.001 gives sqrt(2.001) < 1.4146
strictly, and 2.001 > 2 gives 1/2.001^2 < 1/4 strictly). So
h(c) < 3.0792 on c in (1, 2.001].

*Collecting.* With c <= 2.001: 2c <= 4.002; c^(3/2) = c*sqrt(c) <=
2.001*1.4146 < 2.8312 (sqrt(c) < 1.4146 from (F)); 3/c <= 3 (c >= 1).
Therefore

    rk/t + kt + r(r+1) <= 3.0792*u^4 + 4.002*u^3 + 5.8312*u^2 + 2.

For u >= 10^4 the lower-order part satisfies, exactly,

    4.002/u + 5.8312/u^2 + 2/u^4 <= 4.002/10^4 + 5.8312/10^8 + 2/10^16
                                 < 0.000401,

so the sum is at most (3.0792 + 0.000401)u^4 = 3.079601u^4 < 3.08u^4, i.e.

    rk/t + kt + r(r+1) < (77/25) * n^(4/5).

Let A be an affine subspace of F_2^n with
dim A >= d(n) = floor((77/25)*u^4) + 1,
so dim A > (77/25)*u^4 > rk/t + kt + r(r+1). Let w(z)_i := z_(r+1-i) be the
block-reversal bijection of K^r, and let mu_0 in K* be such that
pi(a) = Tr(mu_0 a) (Proposition 1.2(ii)). For x in F_2^n, z = phi(x):

    f_n(x) = Tr( mu_0 * sum_{|I|=t} prod_{i in I} z_i^(2^(i+1)-1) ) = h(w(z)),

because h(w(z)) = Tr(mu_0 sum_{|I|=t} prod_i w(z)_i^(2^(r+2-i)-1)) and
w(z)_i = z_(r+1-i), so each z_j carries the exponent 2^(r+2-i)-1 with
i = r+1-j, i.e. 2^(j+1)-1. The map w composed with phi is injective
F_2-linear, so w(phi(A)) is an affine subspace of K^r of dimension
dim A > rk/t + kt + r(r+1), with 2 <= t < r and k >= r(r+1)/t
(conditions (i)-(iii) above). Theorem 6.3 shows h is nonconstant on
w(phi(A)); hence f_n is nonconstant on A. This proves Theorem A. qed
