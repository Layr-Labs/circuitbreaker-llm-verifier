# The disperser lemma for the payload family

This file re-proves, from first principles and in our own notation, the
prime-degree binary case of the affine-disperser construction of
Ben-Sasson and Kopparty (*Affine Dispersers from Subspace Polynomials*,
author version 2010-03-12, Lemmas 3.7, 3.10, 3.13, Theorem 6.1, Lemma 7.2,
Theorem 7.3, Claim 7.4). Every step below is proved; nothing is imported
as a black box beyond elementary finite-field algebra and linear algebra.

**Lemma D.** Let k be a prime, K the field with 2^k elements, and suppose
integers t, r satisfy 1 <= t < r and k*t >= r*(r+1). Set w_i = r+2-i.
For any nonzero F_2-linear map pi : K -> F_2, the function

    f(a_1,...,a_r) = pi( sum over subsets I of {1..r}, |I|=t,
                          of  product over i in I of  a_i^(2^(w_i)-1) )

is nonconstant on every affine subspace of K^r of F_2-dimension strictly
greater than  kr/t + k*t + r*(r+1).

Throughout, q = 2^k, and "dimension" is over F_2.

## 1. Polynomial representatives and linear maps

Every function K^r -> K is represented by a polynomial, unique once each
variable's degree is below q: reduce modulo X_i^q - X_i. Uniqueness
proves by induction on the number of variables: a reduced polynomial that
vanishes on all of K^r, read as a univariate in its last variable with
coefficients in the earlier variables, has degree < q and q^m >= q roots
for every fixed assignment of the others, so each coefficient polynomial
vanishes identically, and induction applies. Consequently, if the reduced
representative of a function has a term of positive degree, the function
is not constant: the reduced representative of a constant is the constant
itself.

A *linearized* polynomial is L(X) = sum_{j=0}^{k-1} a_j X^(2^j). Each such
L is F_2-linear because Frobenius z -> z^2 is a field automorphism of
characteristic two. Conversely every F_2-linear map K -> K arises this
way: coefficients and maps both form F_2-spaces of dimension k^2, and a
nonzero linearized polynomial has degree at most 2^(k-1) < q, hence
cannot vanish at all q points of K (root counting), so the map from
coefficients to functions is injective, therefore surjective.

The trace Tr(z) = sum_{h=0}^{k-1} z^(2^h) lands in F_2 and is not the
zero function, having degree 2^(k-1) < q. For mu != 0 the functional
z -> Tr(mu z) is nonzero (multiply by mu^(-1)), and mu -> Tr(mu .) is an
injective linear map between two k-dimensional spaces, hence bijective.
So every nonzero linear pi : K -> F_2 equals Tr(mu .) for exactly one
mu != 0; we fix this mu for the rest of the file.

For 0 <= b <= q-1 let wt(b) count the ones in the binary expansion of b.
If a monomial X^b arises as a product of h powers-of-two factors
X^(2^(c_1)) ... X^(2^(c_h)), possibly after reductions X^(2^k) -> X, then
wt(b) <= h: reduction X^(2^k) = X identifies exponent digits cyclically
modulo k, and each merge of two equal digits 2^j + 2^j = 2^(j+1), as well
as each cyclic wrap, weakly decreases the number of digits. Moreover
wt(b) = h only if no merge ever happened, i.e. only if the residues
c_1,...,c_h are pairwise distinct modulo k. This degree subadditivity
survives further reductions, since those are themselves merges. We use it
constantly: a monomial of binary weight h in a variable can only be
produced by h distinct Frobenius factors in that variable.

## 2. Subspace polynomials and their zero gaps

For a linear subspace V <= K of dimension d >= 0 define

    P_V(X) = product over v in V of (X - v) = product over v in V of (X+v).

P_V is monic of degree 2^d. It is linearized: for d = 0 this is X, and

    P_(U+<b>)(X) = P_U(X) * P_U(X+b) = P_U(X)^2 + P_U(b) * P_U(X),

using P_U(X+b) = P_U(X) + P_U(b) for linearized P_U. So P_V is a linear
map with kernel exactly V: its roots are the points of V, and they are
distinct because P_V divides X^q - X (every element of K is a root of
X^q - X, which is squarefree since its derivative is 1).

Let V' be the image of P_V, of dimension k - d, and let Q = P_(V'), the
*image polynomial* of V, monic linearized of degree 2^(k-d).

*Claim 1.* Q(P_V(X)) = X^q - X. The left side is monic of degree
2^k = q and vanishes on all of K (P_V maps K into V'), so it equals the
unique monic degree-q polynomial with all q field elements as roots.

*Claim 2.* P_V(Q(X)) = X^q - X. Let D = P_V(Q(X)) - (X^q - X). Using
Claim 1, D(P_V(X)) = P_V(Q(P_V(X))) - P_V(X)^q + P_V(X)
= P_V(X^q - X) - P_V(X)^q + P_V(X) = 0, because P_V(X)^q = P_V(X)
pointwise on K. A nonzero polynomial composed with the nonconstant
polynomial P_V is nonzero (degrees multiply in the domain F_2[X]), so
D = 0.

*Claim 3.* Q maps K onto V. For z in K, Claim 2 gives
P_V(Q(z)) = z^q - z = 0, so Q(z) lies in ker P_V = V. Q is linear with
2^(k-d) distinct roots (Claim 2 and squarefreeness), so its kernel has
dimension k-d and its image has dimension d; containment in V forces
equality.

Write Q(X) = sum_{j=0}^{k-1} a_j X^(2^j); a_j = 0 for j > k-d because the
degree is 2^(k-d), and a_(k-d) = 1 by monicity. So the final interval of
coefficient indices k-d+1, ..., k-1 is a cyclic *zero run* of length d-1.

*Claim 4 (no run of length d).* Q has no cyclic run of d consecutive zero
coefficients. A cyclic run can be moved to the final d positions by
replacing Q with a Frobenius power Q^(2^h) reduced modulo X^q - X (this
cyclically shifts coefficient indices by h). That power is nonzero, has
degree at most 2^(k-d-1), and has the same 2^(k-d) roots as Q (its value
at z is Q(z)^(2^h)), contradicting root counting. In particular
a_0 != 0: together with the final interval, a_0 = 0 would create a run of
length d.

*Claim 5 (prime k: uniqueness of the long run).* Let k be prime and
d >= 2. The final interval k-d+1,...,k-1 is the ONLY cyclic zero run of
length d-1. Suppose another exists; shift it to the final positions by
some Frobenius power 2^h with 0 < h < k, giving a reduced polynomial R
with the same 2^(k-d) roots as Q. Its coefficient at X^(2^(k-d)) is
c = a_((k-d-h) mod k)^(2^h), and its degree is at most 2^(k-d). If c = 0
then R is nonzero of degree below 2^(k-d) with 2^(k-d) roots, which is
impossible; so c != 0 and c^(-1) R is monic of degree exactly 2^(k-d).
Then c^(-1) R - Q has degree below 2^(k-d) and at least 2^(k-d) roots,
so R = c Q. Evaluating at w in K, R(w) = Q(w)^(2^h), so
Q(w)^(2^h) = c Q(w) for all w; every z in Q(K) = V satisfies
z^(2^h) = c z. Fix z_0 != 0 in V (d >= 2). Every quotient z/z_0 with
z in V satisfies (z/z_0)^(2^h) = cz/(cz_0) = z/z_0, so it is fixed by
Frobenius^h and by Frobenius^k = identity; since k is prime
and 0 < h < k, gcd(h,k) = 1 and Bezout expresses Frobenius as a composition
of powers of Frobenius^h and Frobenius^k, so z/z_0 is fixed by Frobenius
itself, landing in F_2. Then |V| <= 2, contradicting d >= 2.

Consequently, for prime k and d >= 2, every cyclic zero run of Q other
than the final interval has length at most d-2.

## 3. A surviving coefficient for independent subspaces

Let V_1,...,V_t <= K be linear subspaces with dimensions d_i >= 2, and
let Q_i = sum_j a_(i,j) X^(2^j) be their image polynomials. Fix exponents
e_i with 1 <= e_i <= d_i - 2 and write h_i = e_i + 1,
delta_i = d_i - h_i >= 1. Assume

    sum_i delta_i > k - max_i d_i + 1.                       (H)

*Claim 6.* Under (H), the polynomial

    Tr( mu * product_i Q_i(X_i)^(2^(h_i)-1) )

has a nonzero coefficient at some monomial product_i X_i^(b_i) with
wt(b_i) = h_i for every i.

For a k-by-k matrix over K indexed by Z/k define
T_i[u,v] = a_(i,v-u)^(2^u); this is the coefficient of X^(2^v) in
Q_i(X)^(2^u) after cyclic reduction. Fix row set A_i = {0,...,e_i} (h_i
rows). Choosing a column set B_i of size h_i and expanding
Q_i^(2^(h_i)-1) = product_{u=0}^{h_i-1} Q_i^(2^u), the coefficient of
X_i^(b_i), b_i = sum_{v in B_i} 2^v, in one trace summand with Frobenius
shift rho is

    sum over bijections sigma : A_i -> B_i - rho
        of  product_{u in A_i} a_(i, sigma(u)-u)^(2^u)
    = perm( T_i[A_i, B_i - rho] ).                            (1)

Monomials of weight h_i can only be assembled from h distinct
power-of-two factors hitting h distinct columns (Section 1); anything
else (repeated columns, merges) lands on a monomial of strictly smaller
weight and does not interfere. In characteristic two the permanent equals
the determinant, though we only need that it is nonzero when a triangular
argument applies.

*Construction of B_i.* Reindex so that d_1 is maximal. Let J be the first
index with sum_{i<=J} delta_i >= k - d_1 + 1; it exists by (H). For
i <= J put ell_i = sum_{j<i} delta_j; then ell_i <= k - d_1 <= k - d_i
for i <= J, and the crossing sum satisfies
sum_{i<=J} delta_i <= (k - d_1) + delta_J <= k - 2. For i > J put
ell_i = 0 (these sources are not needed for coverage).

Fix one source i <= J and abbreviate d = d_i, e = e_i, ell = ell_i,
a_j = a_(i,j), delta = delta_i = d - e - 1. Recall 1 <= e <= d-2, so
k-d+1 <= k-e-1 = k-d+delta, and note ell <= k-d.

* If a_ell != 0, take the ordered column list
  B = (ell, k-d+1, k-d+2, ..., k-d+e).
* If a_ell = 0, let l be the largest index < ell with a_l != 0 (it
  exists: a_0 != 0 and ell >= 1) and j the smallest index > ell with
  a_j != 0 (it exists: a_(k-d) = 1 and ell <= k-d). Then l < ell < j <= k-d
  and positions l+1,...,j-1 are zero, a cyclic run shorter than the final
  one, hence of length z = j-l-1 <= d-2 by Claim 5. Put v = min(z, e) and
  take B = (j, l+1, ..., l+v, k-d+v+1, ..., k-d+e).

Both lists have 1 + e = h_i entries, all distinct: in the first case
ell <= k-d < k-d+1; in the second l+v < l+z+1 = j <= k-d < k-d+v+1.

*The rho = 0 minor is lower triangular with nonzero diagonal.* Order rows
u = 0,...,e. In the first case pair row 0 with column ell and row
s in {1..e} with column k-d+s. The diagonal entries are
a_ell (row 0) and a_(k-d) = 1 (rows 1..e). Entries above the diagonal
are T[0, k-d+s] = a_(k-d+s) with s >= 1, and T[s, k-d+s'] with s' > s,
i.e. a_(k-d+s'-s) with s'-s >= 1: all indices lie in the final zero
interval k-d+1,...,k-1.

In the second case pair row 0 with column j, rows 1..v with columns
l+1,...,l+v, and rows v+1..e with columns k-d+v+1,...,k-d+e. The diagonal
entries are a_j, then a_l^(2^s) = Frobenius powers of a_l (nonzero), then
a_(k-d) = 1. Entries above the diagonal fall in the gap l+1,...,j-1
(rows 0..v-1 against columns strictly right of their diagonal up to l+v
or j), or in the final interval (any row against a final column strictly
right of its diagonal): for row s in {1..v} against final column
k-d+v+s'', the coefficient index is k-d+v+s''-s >= k-d+1; for rows
v+1..e against later final columns the index is k-d+s'-s with s' > s.
So every above-diagonal entry has its coefficient index in a zero run.

Hence perm(T_i[A_i, B_i]) is a product of nonzero diagonal entries
(block-lower-triangular permanent), in particular nonzero.

*Zero columns for shifts.* A chosen column b of B_i makes the column
b - rho of T_i[., .] vanish on all rows u = 0..e exactly when every index
b-rho-u with 0 <= u <= e lies in the final interval k-d+1,...,k-1, i.e.
when rho lies in the cyclic interval [b+1, b+delta] (both ends mod k):
the interval has exactly delta = d-1-e members since
(b - (k-1)) ... (b - (k-d+1+e)) is that cyclic interval.             (2)

*Coverage of all nonzero shifts.* Source 1 has ell_1 = 0. Its columns
include 0, which by (2) covers rho in [1, delta_1], and the final columns
k-d_1+1,...,k-d_1+e_1: the final column b_s = k-d_1+s (1 <= s <= e_1)
covers [k-d_1+s+1, k-d_1+s+delta_1] by (2). Consecutive intervals overlap
or abut (delta_1 >= 1), and since e_1+delta_1 = d_1-1 their union is
[k-d_1+2, k-1].

Each middle source 1 < i <= J covers the whole interval
[ell_i+1, ell_i+delta_i]:
if a_ell != 0, the column ell itself is chosen and covers exactly that
interval. If a_ell = 0: when z < e the chosen columns include every
position l+1,...,l+z, in particular ell. When z >= e the columns
l+1,...,l+e cover [l+2, l+e+delta] = [l+2, l+d-1] (using
e + delta = d - 1), which reaches j = l+z+1 <= l+d-1, and the column j
extends the coverage through j+delta; since
ell + delta <= (j-1) + delta <= (l+d-2) + delta < l+d-1+delta, the union
[l+2, l+d-1+delta] contains [ell+1, ell+delta].

Because sum_{i<=J} delta_i >= k - d_1 + 1, the middle sources cover
[delta_1 + 1, k - d_1 + 1], which meets both the head [1, delta_1] and
the tail [k-d_1+2, k-1] of source 1. Together all rho in 1..k-1 are
covered. Sources after J were not used and cannot remove coverage.

*Conclusion of Claim 6.* For every rho != 0 some factor
perm(T_i[A_i, B_i - rho]) in trace-summand rho has a zero column, hence
vanishes, so the coefficient of product_i X_i^(b_i) equals the rho = 0
summand, mu * product_i perm(T_i[A_i, B_i]), a product of nonzero field
elements. This proves Claim 6.

## 4. Block decomposition of an affine subspace

Let A <= K^r be an affine subspace of dimension D. Filter the direction
space U by the number of leading zero blocks: U^(i) = U intersected with
(0,...,0, K^(r-i)). Each quotient U^(i)/U^(i+1) embeds into K by reading
the (i+1)-st block; choose an F_2-basis of each quotient. The span V_i of
the embedded basis in block i+1, together with the tail coordinates of
the chosen basis vectors (which define linear maps sigma_ij : V_j -> K
for j < i), parametrizes A: there is gamma in A with

    A = { ( gamma_i + z_i + sum_{j<i} sigma_ij(z_j) )_(i=1..r) : z_i in V_i },
    D = sum_i dim V_i.                                        (3)

This is Gaussian elimination respecting the block order; the chosen
vectors form a basis of U because the filtration quotients compose to it.

Let Q_i be the image polynomial of V_i (Section 2; Q_i = 0 if V_i = 0),
and let Q_ij be the linearized polynomial of the linear map
sigma_ij composed with Q_j (Section 1 represents all linear maps). Define

    R_i(X) = gamma_i + Q_i(X_i) + sum_{j<i} Q_ij(X_j).

The map R = (R_1,...,R_r) : K^r -> K^r has image exactly A by (3), since
Q_i maps K onto V_i (Claim 3).

## 5. Choosing t high-dimensional blocks

Assume D > kr/t + kt + r(r+1). Let S = { i : dim V_i > r+1 }.

* |S| >= t: the indices outside S contribute at most r(r+1) dimensions in
  total, so sum_(i in S) d_i >= D - r(r+1) > kt; since each d_i <= k and
  |S| <= r, we get |S| * k >= sum_(i in S) d_i > kt, so |S| > t and
  |S| >= t + 1.
* sum_{i in S} d_i >= D - r(r+1) > kr/t + kt.

Choose I0 consisting of t elements of S with the largest dimensions
(breaking ties by index). The chosen sum is at least t/r of the total
over S because |S| <= r. For i in I0:

    d_i - (r+1) >= d_i - (r+1),   and summing,
    sum_(i in I0) (d_i - (r+1)) > (t/r)(kr/t + kt) - t(r+1)
                                = k + kt^2/r - t(r+1) >= k,    (4)

the last step exactly by the guard kt >= r(r+1).

Set e_i = r+1-i for i in I0. Then e_i >= 1 (i <= r), and

    d_i - 1 - e_i = d_i - r - 2 + i >= (r+2) - r - 2 + i = i >= 1,

using d_i > r+1, so 1 <= e_i <= d_i - 2; also h_i = e_i + 1 = w_i =
r+2-i <= r+1 < d_i <= k, so h_i < k. With delta_i = d_i - w_i:

    delta_i = d_i - (r+2-i) = (d_i - (r+1)) + (i-1) >= d_i - (r+1),

so sum_(i in I0) delta_i >= k by (4), and the hypothesis (H) of Claim 6
holds: k >= k - max d_i + 1 because max d_i >= 2. Claim 6 supplies a
monomial M = product_(i in I0) X_i^(b_i), wt(b_i) = w_i, with nonzero
coefficient in Tr( mu * product_(i in I0) Q_i(X_i)^(2^(w_i)-1) ).

## 6. Triangular substitution preserves the coefficient

We claim M has the same coefficient in the substituted polynomial
Tr( mu * sum over |I| = t of product_(i in I) R_i(X)^(2^(w_i)-1) )
restricted to the summand I = I0; Section 7 removes the other summands.

Expand R_i^(2^(w_i)-1) = product_{u=0}^{w_i-1} R_i^(2^u) and the trace
into its Frobenius summands; fix one summand (an overall 2^rho shift).
Process the selected indices in DECREASING order. Let i* be the largest
element of I0. Since R_j mentions X_{j'} only for j' <= j (Section 4),
the variable X_i* is mentioned among the factors of the I0 summand only
by the factors R_i*^(2^u). The monomial M requires binary weight
w_i* = wt(b_i*) in X_i*; by subadditivity (Section 1) at most the w_i*
factors R_i*^(2^u) can contribute to it, and weight exactly w_i* forces
every one of them to contribute a single distinct power of X_i* with no
merges. Each such contribution is a term a X^(2^c) of Q_i*(X_i*) itself
(the constant gamma_i* contributes degree 0, and the Q_i* j' terms for
j' < i* do not mention X_i*). Therefore the X_i*-degree-(b_i*) part of
the summand is exactly the corresponding part of
Q_i*(X_i*)^(2^(w_i*)-1), multiplied by an honest polynomial in the
variables X_i, i < i*. Strip it and repeat with the next-largest element
of I0. Descending induction shows the coefficient of M in the
rho-summand of the I0 summand equals its coefficient in the rho-summand
of Tr( mu product Q_i^(2^(w_i)-1) ). Summing over rho, the coefficients
agree.

## 7. No other symmetric summand contributes to M

Consider a summand J != I0, |J| = t, in the symmetric sum, still inside a
fixed trace summand rho. A factor R_j^(2^u) mentions only variables
X_1,...,X_j (Section 4: R_j reads X_{j'} for j' <= j). Let n_ji count
the factors of the J summand whose chosen term mentions X_i, for
i in I0. Producing M needs binary weight w_i in each X_i, so

    sum_(j in J) n_ji >= w_i    for every i in I0,               (5a)
    sum_(i in I0) n_ji <= w_j   for every j in J,                (5b)
    n_ji > 0 implies i <= j.                                     (5c)

Since w_i = r+2-i decreases strictly in i, (5c) gives
n_ji > 0 => w_i >= w_j. Now

    t = sum_(i in I0) w_i / w_i
      <= sum_(i in I0, j in J) n_ji / w_i          by (5a)
      <= sum_(i in I0, j in J) n_ji / w_j          since w_i >= w_j when n_ji>0
      <= sum_(j in J) w_j / w_j = t                by (5b).        (6)

Equality throughout forces n_ji = 0 whenever i < j (strict w_i > w_j
would make the middle inequality strict). So every used factor has
i = j. Each i in I0 has positive demand w_i > 0 in (5a), so some j in J
has n_ji > 0; combined with n_ji > 0 => i = j, we get i in J.
Hence I0 is contained in J, and |I0| = |J| = t gives J = I0 —
contradiction. So no other summand contributes to the coefficient of M.

## 8. Assembling Lemma D

The coefficient of M in f(R(X)) is now computed: by Section 7 only the
I0 summand contributes, by Section 6 its contribution equals the nonzero
coefficient from Claim 6, and by Section 1 a positive-degree term in the
reduced representative means f composed with R is a nonconstant
polynomial function on K^r. Since R maps K^r onto A (Section 4), f is
nonconstant on A. This holds for every affine subspace A of dimension
D > kr/t + kt + r(r+1), proving Lemma D.

All hypotheses used were: k prime (Claim 5), 1 <= t < r and
kt >= r(r+1) (Section 5), and pi != 0 (Section 1). The lemma is applied
in proof.md with pi = the constant-coordinate functional, whose mu is
nonzero because pi(1) = 1.
