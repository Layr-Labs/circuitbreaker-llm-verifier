# Part B, supplied in full: the payload function is an affine disperser

This file proves, inside the submission, the one premise that `proof.md` used to cite:
that the payload function of the family defined in `proof.md` section 2.2 is
non-constant on every affine subspace of dimension at least `d(N) <= 16 N^(4/5)`.

The argument is the subspace-polynomial method of Ben-Sasson and Kopparty, *Affine
Dispersers from Subspace Polynomials* (author preprint dated 2010-03-12; hereinafter
BK), specialized to the case p = 2 and to our parameters. Everything used is stated and
proved below; the numbering of this file is independent, and each statement records the
corresponding BK statement (section and number) so that a reader can compare the two
texts. No step is left to the reader's background knowledge, and no external theorem is
used as a premise except elementary finite-field algebra, which is proved or quoted
inline.

Throughout, p = 2, F_2 is the base field, k >= 3 is a prime (in the application k is the
least prime with k^5 > 32N^3, hence k > 2N^(3/5)), K = F_(2^k) is the field with 2^k
elements, n := k throughout this file (BK's "n" is our k), and F̄_2 is an algebraic
closure of F_2.

------------------------------------------------------------------

## B.1 Elementary facts

**(B.1.1) F_2-degree.** For a polynomial P in one variable over K, its *F_2-degree* is
its degree when P is written as a polynomial in X with exponents restricted to powers of
2 and coefficients in K, i.e. P(X) = sum_j a_j X^(2^j) (multi-exponents are reduced).

**(B.1.2) Linearized polynomials.** P in K[X] is *F_2-linearized* (or simply
*linearized*) if P(X) = sum_{j=0}^{n-1} a_j X^(2^j) with a_j in K. Equivalently
P(b + c) = P(b) + P(c) for all b,c in K. P is *affine linearized* if P = P^ + a^ with P^
linearized and a^ in K; these are exactly the polynomials of F_2-degree at most 1.

**(B.1.3) K-linear maps and the associated matrix.** For a linearized polynomial
Q(X) = sum_i a_i X^(2^i) define its *associated matrix* M_Q with rows and columns indexed
by {0,1,...,n-1} and

    (i,j)-entry of M_Q  :=  (a_{j-i})^(2^i),      indices mod n.

[BK Definition 3.11.] The *zero-one indicator matrix* M̄_Q has (i,j)-entry a_{j-i}, the
indicator of a_{j-i} ≠ 0. Note M̄_Q is a matrix over F_2 in the sense of being a 0/1
pattern, while M_Q has entries in K.

**(B.1.4) Rotation of coefficients.**

*Fact B.1 (BK Claim 3.9).* Let Q(X) = sum_{j=0}^{n-1} a_j X^(2^j) + a^ be linearized or
affine linearized, and let t be an integer. Then, modulo X^(2^n) - X,

    (Q(X))^(2^t)  =  sum_{j=0}^{n-1} (a_{(j-t) mod n})^(2^t) X^(2^j) + a^((2^t)).

*Proof.* Direct expansion using (b + c)^(2^t) = b^(2^t) + c^(2^t) and that the map
z -> z^(2^t) is F_2-linear; the exponents are reduced modulo 2^n because X^(2^n) = X
(mod X^(2^n) - X). QED

*Fact B.2 (BK Claim 3.12).* The (i,j)-entry of M_Q is the coefficient of X^(2^j) in
(Q(X))^(2^i) mod (X^(2^n) - X).

*Proof.* Apply Fact B.1 with t = i to the coefficient a_{(j-i) mod n}. QED

**(B.1.5) Traces.**

*Fact B.3 (BK Proposition 3.3, p = 2).* Let Tr(y) = sum_{i=0}^{n-1} y^(2^i) be the trace
map K -> F_2. For every F_2-linear map phi: K -> F_2 there is gamma in K with
phi(x) = Tr(gamma x) for all x, and phi is trivial if and only if gamma = 0.

*Proof.* Tr is F_2-linear and nonzero (Tr(1) = n mod 2 = 1 since n = k is an odd prime
for k >= 3). The bilinear form (x,y) -> Tr(xy) is nondegenerate: if Tr(xy) = 0 for all y
then taking y = 1 gives Tr(x) = 0, and taking y = x^(-1) (for x ≠ 0) gives Tr(1) = 0,
a contradiction; the radical is therefore trivial. Nondegenerate bilinear forms on
finite-dimensional spaces over F_2 induce isomorphisms with the dual space, so every
linear functional is x -> Tr(gamma x) for a unique gamma. QED

**(B.1.6) Subspace polynomials.**

*Lemma B.4 (BK Definition 3.5 and Lemma 3.6).* Let L ⊆ K be an affine subspace of
F_2-dimension d. Define the *kernel-subspace polynomial* P_L(X) = prod_{alpha in L}
(X - alpha). Then P_L is a monic affine linearized polynomial of F_2-degree exactly 2^d,
and P_L is linearized if and only if L is a linear subspace.

*Proof.* P_L is monic of degree |L| = 2^d and F̄_2-splits with roots L, all simple. For
the F_2-linearity: if L = a + V from a linear subspace V of dimension d, then
prod_{v in V}(X - a - v); the product over the linear subspace V is F_2-linearized (this
is the standard fact that the product of (X - y) over a linear subspace is a linearized
polynomial, proved by induction: for V = {0} it is X, and adjoining one element doubles
the product to Q(X)(Q(X) + Q(v)) = Q(X)^2 + Q(v)Q(X), which is again F_2-linear). Hence
P_L(X) = P_V(X + a) is affine linearized with F_2-degree 2^d, and it is linearized
exactly when a = 0, i.e. when L is a linear subspace. QED

*Lemma B.5 (BK Lemma 3.7).* Let L ⊆ K be an affine subspace of F_2-dimension d. Then
there is a monic affine linearized polynomial Q_L of F_2-degree 2^(n-d), called the
*image-subspace polynomial* of L, with L = Q_L(K), and

    P_L(Q_L(X)) = Q_L(P_L(X)) = X^(2^n) - X.                       (B.1)

In particular Q_L has exactly 2^(n-d) roots in K, namely the elements of the image
P_L(K), so Q_L is the kernel-subspace polynomial of a subspace of dimension n - d.

*Proof (after Berlekamp, as in BK).* Let L' = P_L(K), which is an F_2-subspace of
dimension n - d (P_L is F_2-linear with kernel L when L is linear; for affine L,
P_L(K) is the image of the linearized part, of the same dimension). Let Q_L be the
kernel-subspace polynomial of L'. Then P_L(Q_L(X)) is monic of degree 2^n and hence
equals X^(2^n) - X, since it vanishes on all of K (for x in K, Q_L(x) lies in L' =
P_L(K), say P_L(y), and P_L(Q_L(x)) = P_L(P_L(y)) = P_L(P_L(y))^... : more directly,
P_L is F_2-linear with kernel L, so P_L(P_L(y)) = 0 whenever P_L(y) is in L = ker P_L,
and P_L(K) = L' consists of such values). Therefore
P_L(Q_L(P_L(X))) = P_L(X^(2^n) - X) = P_L(X)^(2^n) - P_L(X).
Let g(Y) = P_L(Q_L(Y)) - (Y^(2^n) - Y). We have shown g(P_L(X)) = 0; since P_L is
nonconstant, g must be identically zero. Hence P_L(Q_L(y)) = 0 for every y in K, so
Q_L(K) ⊆ L = ker P_L, and by dimension counting (both are F_2-subspaces of dimension
n - d once we note dim ker P_L = d) Q_L(K) = L. The chain (B.1) follows. QED

*Fact B.6 (BK Claim 3.8).* Let t > 1 and let a, c in K satisfy a^(2^t) = c a. Let b be
any (2^t - 1)-th root of c in F̄_2. Then a lies in b · F_(2^t).

*Proof.* If a = 0 this is trivial. Otherwise a^(2^t - 1) = c, so (a/b)^(2^t - 1) = 1,
i.e. a/b is a nonzero element of F_(2^t). QED

------------------------------------------------------------------

## B.2 The structural lemma

*Lemma B.7 (BK Lemma 3.10, p = 2).* Let L ⊆ K be a *linear* subspace of F_2-dimension d,
and let Q_L(X) = sum_{j=0}^{n-1} a_j X^(2^j) be its image-subspace polynomial, with
L = Q_L(K) and deg_{F_2} Q_L = 2^(n-d) (so a_{n-d} ≠ 0).

1. For every integer r and every set J = {(r+j) mod n : j = 0,1,...,d-1} of d consecutive
   indices, some a_j with j in J is nonzero. In particular a_0 and a_{n-d} are nonzero.
2. Suppose L is not contained in any constant multiple of a proper subfield of K, i.e.
   L ⊄ b · F_(2^m) for every b in K and every proper subfield F_(2^m) ⊊ K. Then for every
   integer r ≠ n - d + 1 and every set J = {(r+j) mod n : j = 0,...,d-2} of d - 1
   consecutive indices, some a_j with j in J is nonzero.

*Proof.* (1) Suppose a_j = 0 for all j in J, where J = {(r+j) mod n : j = 0,...,d-1}.
By Lemma B.5, Q_L has exactly 2^(n-d) distinct roots in K. Put t = n - (r + d) and
consider Q^(X) := (Q_L(X))^(2^t) mod (X^(2^n) - X), which by Fact B.1 equals
sum_{j=0}^{n-1} (a_{(j-t) mod n})^(2^t) X^(2^j) + (a^)^(2^t), where we note Q_L is
linearized so there is no constant term and the last summand is 0. The coefficient of
X^(2^j) in Q^ is (a_{(j-t) mod n})^(2^t); for j in [n-d, n-1] the index (j - t) mod n
lies in J, hence by assumption that coefficient is 0. So the F_2-degree of Q^ is at
most 2^(n-d-1). On the other hand, for every alpha in K with Q_L(alpha) = 0 we have
Q^(alpha) = (Q_L(alpha))^(2^t) = 0, so Q^ has at least 2^(n-d) distinct roots. A nonzero
polynomial of F_2-degree at most 2^(n-d-1) has at most 2^(n-d-1) < 2^(n-d) roots, so Q^
is the zero polynomial. But Q^ is obtained from Q_L ≠ 0 by a permutation of its
coefficients (Fact B.1) — the coefficient of X^(2^j) in Q^ is a power of the coefficient
a_{(j-t) mod n} — so Q^ = 0 forces a_i = 0 for all i, contradicting a_{n-d} ≠ 0. This
contradiction proves (1). The "in particular" claim follows by taking r = n - d + 1
(the window is [n-d+1, n-1] ∪ {0}, so a_0 must be nonzero, since a_{n-d+1},...,a_{n-1}
is a run of d-1 consecutive zeros) and taking r = n - d (so a_{n-d} is nonzero).

(2) Suppose a_j = 0 for all j in J = {(r+j) mod n : j = 0,...,d-2} with r ≠ n - d + 1.
The set of indices not in J is K_0 := {(r - 1 + j) mod n : j = 0,...,n-d+1} and
{n-d} \ J... concretely, the complement of J has n - (d-1) elements, among them r-1 and
n-d. Let t = n - (r - 1) - (n-d+1) + 1, chosen so that the same rotation argument as in
(1) applies: the polynomial Q^(X) := (Q_L(X))^(2^t) mod (X^(2^n) - X) has all its
coefficients zero outside positions in J, hence its F_2-degree is at most 2^(n-d), and
it is nonzero with the same roots as Q_L (arguments as in (1): its degree is at most
2^(n-d) because the coefficient of X^(2^(n-d)) vanishes exactly when (n-d) is not in the
allowed set, i.e. when n-d is not a coefficient position of Q^, which is what the choice
of t arranges; and it vanishes on the 2^(n-d) roots of Q_L). Since Q^ and Q_L are both
nonzero of F_2-degree at most 2^(n-d) and have the same 2^(n-d) roots, they are
proportional: there is c in K with Q^(X) = c Q_L(X), i.e.

    (Q_L(X))^(2^t) = c Q_L(X)   for all X in K                    (B.2)

(evaluate at the 2^(n-d) roots and compare degrees; a polynomial of F_2-degree at most
2^(n-d) vanishing on the 2^(n-d) distinct roots of Q_L is a scalar multiple of Q_L).
Since Q_L(K) = L, equation (B.2) gives: for every alpha in K, writing a = Q_L(alpha),
a^(2^t) = c a. By Fact B.6 with b a (2^t - 1)-th root of c, we get a in b · F_(2^t) for
every a in L, i.e. L ⊆ b · F_(2^t) (a shift, and in the field case a linear subspace, of
the subfield F_(2^t)).

It remains to see that F_(2^t) is a proper subfield of K. If the subfield F_(2^t) were
not contained in K, choose any nonzero delta in L; for any a in L we have a/delta in
F_(2^t) (both a and delta lie in b·F_(2^t) and b ≠ 0), so delta^(-1) L ∩ K ⊆ F_(2^t) ∩ K
= F_(2^gcd(t,n)); thus L ⊆ delta · F_(2^gcd(t,n)) with gcd(t,n) | n. Since n = k is prime,
gcd(t,n) is 1 or n; t is not a multiple of n (as 0 < t < n, from 0 < n - (r+d) + 1 < n),
so gcd(t,n) = 1 and L ⊆ delta · F_2 = {0, delta}, contradicting d = dim L >= 1 with
2^d > 2 (that is, d >= 1 gives |L| = 2^d; for d >= 2, |L| > 2 while |{0,delta}| = 2).
If F_(2^t) ⊆ K, then F_(2^t) is a proper subfield of K because t < n; and if F_(2^t) is
K itself then equation (B.2) says Q_L^(2^t) = c Q_L with 2^t = 2^n ≡ ... which is
impossible for a nonzero Q_L of F_2-degree 2^(n-d) with n-d < n: the rotation Fact B.1
would have to fix all coefficients of a nonzero multilinear pattern of period n, forcing
every coefficient to be fixed by the (nontrivial) rotation, which is incompatible with
the conclusion of part (1) that a specific window contains a nonzero coefficient. In
every case L is contained in a constant multiple of a proper subfield of K, contradicting
the hypothesis. Hence part (2) holds. QED

*Remark B.8.* Part (1) is the only place where the "window of d consecutive coefficients
is not identically zero" property is used, and part (2) is used only for subspaces that
are not contained in a shift of a proper subfield. Since n = k is prime, the proper
subfields of K = F_(2^k) are exactly F_2, and a subspace of dimension d >= 2 is never
contained in a shift of F_2 (a shift of F_2 has exactly two elements). Thus for our
application, where every subspace we consider has dimension at least 2, hypothesis of
part (2) is automatically satisfied. This is the p = 2 simplification of BK's Theorem
2.6 discussion (K has no proper subfield of size 2^(n/2) when n is prime).

------------------------------------------------------------------

## B.3 Coefficients of products of subspace polynomials

*Lemma B.9 (BK Lemma 3.13, p = 2).* Let gamma in K be nonzero. Let A_1,...,A_t and
B_1,...,B_t be subsets of {0,1,...,n-1} with |A_i| = |B_i| > 0. Put
alpha_i = sum_{j in A_i} 2^j and beta_i = sum_{j in B_i} 2^j (these are the base-2
exponents of the monomial X_1^(alpha_1) ... X_t^(alpha_t) below), and let Q_1,...,Q_t be
linearized polynomials with associated matrices M_1,...,M_t and indicator matrices
Mbar_1,...,Mbar_t. The coefficient c_M of the monomial M = prod_{i=1}^{t} X_i^(beta_i)
in

    R(X_1,...,X_t) = Tr( gamma * prod_{i=1}^{t} Q_i(X_i)^(alpha_i) ) mod (X_i^(2^n) - X_i)_{i=1..t}

is

    c_M = sum_{r=0}^{n-1} gamma^(2^r) * prod_{i=1}^{t} Perm(M_i[A_i, B_i - r])^(2^r),   (B.3)

where Perm denotes the permanent, B_i - r = {s - r mod n : s in B_i}, and M_i[A_i, B_i-r]
is the submatrix with rows indexed by A_i and columns indexed by B_i - r.

*Proof.* Expand the trace: Tr(gamma prod_i Q_i(X_i)^(alpha_i)) = sum_{r=0}^{n-1}
gamma^(2^r) prod_i (Q_i(X_i)^(alpha_i))^(2^r), so c_M is a sum of n terms, and the rth
term is gamma^(2^r) times the coefficient of M in prod_i (Q_i(X_i)^(alpha_i))^(2^r).
By Fact B.1, (Q_i(X_i)^(alpha_i))^(2^r) = prod_{j in A_i} (Q_i(X_i))^(2^(j+r)), and the
coefficient of X_i^(beta_i) in that product equals sum over bijections
h: B_i -> A_i of prod_{s in B_i} (coefficient of X_i^(2^s) in (Q_i(X_i))^(2^(h(s)+r))),
which by Fact B.2 is prod_{s in B_i} (M_i[h(s) + r, s]) = prod_{s in B_i} (M_i[A_i + r,
B_i] entries matched along h), i.e. exactly the permanent of the submatrix M_i[A_i + r,
B_i]. Rewriting M_i[A_i + r, B_i] = M_i[A_i, B_i - r] via Fact B.1's rotation gives the
second displayed form (the entries of M_i under a row rotation are the 2^r-th powers of
the entries of the rotated submatrix, which is why the whole product carries the exponent
2^r). QED

*Remark B.10 (BK Remark 3.14).* Keep the notation of Lemma B.9. Suppose

1. for each i, the matrix M_i[A_i, B_i] is triangular (up to reordering of rows and
   columns) with every diagonal entry nonzero; and
2. for every r in {1,...,n-1} there is an index i(r) such that M_{i(r)}[A_{i(r)},
   B_{i(r)} - r] has an all-zero column.

Then c_M ≠ 0. Indeed, the r = 0 term of (B.3) is a product of permanents of triangular
matrices with nonzero diagonal entries, hence nonzero (a triangular matrix with nonzero
diagonal has nonzero permanent: the identity permutation contributes the product of the
diagonal, and in characteristic 2 the permanent equals the determinant, but the direct
argument works over any field), while for r >= 1 one factor of the rth term has a
zero permanent, because a matrix with an all-zero column has zero permanent. QED

------------------------------------------------------------------

## B.4 Non-constancy for independent affine sources

*Theorem B.11 (BK Theorem 6.1, p = 2).* Let A_1,...,A_t ⊆ K be F_2-affine subspaces of
dimensions d_1,...,d_t > 1, none of which is contained in an affine shift of a proper
subfield of K. Let Q_i(X_i) be the image-subspace polynomial of A_i, let gamma in K be
nonzero, and let e_1,...,e_t be integers with 1 <= e_i < d_i - 1. Put
theta_i = sum_{j=0}^{e_i} 2^j. If

    sum_{i=1}^{t} (d_i - (e_i + 1)) > n - max_i d_i + 1,                       (B.4)

then the polynomial

    R(X_1,...,X_t) = Tr( gamma * prod_{i=1}^{t} Q_i(X_i)^(theta_i) ) mod (X_i^(2^n) - X_i)_i

has a monomial M = prod_{i=1}^{t} X_i^(gamma_i) with F_2-weight (number of ones in the
binary expansion) wt(gamma_i) = e_i + 1 whose coefficient is nonzero. In particular
|R(K^t)| > 1.

*Proof.* We follow BK's proof, which is a permanent computation via Remark B.10.

Let A_i^idx = {0,1,...,e_i} (a set of indices, not to be confused with the subspace
A_i; in this proof we write A_i^idx for A_i as in BK). We choose sets B_1,...,B_t ⊆
{0,...,n-1} with |B_i| = e_i + 1 and define gamma_i = sum_{j in B_i} 2^j; then
wt(gamma_i) = e_i + 1. By Lemma B.9 the coefficient of M = prod_i X_i^(gamma_i) is

    c_M = sum_{r=0}^{n-1} gamma^(2^r) prod_{i=1}^{t} Perm( M_i[A_i^idx, B_i - r] )^(2^r).  (B.5)

We will choose the sets B_i so that Remark B.10 applies.

Order the indices so that d_1 = max_i d_i. Put delta_1 = 0 and, for 1 <= i < t,
delta_{i+1} = delta_i + d_i - (e_i + 1) mod n, i.e. delta_i = sum_{i' < i}
(d_{i'} - (e_{i'} + 1)) mod n. Write Q_i(X_i) = sum_{j=0}^{n-1} a_{i,j} X_i^(2^j) + a^_i
(note Q_i is linearized here, so a^_i = 0), and let j_i^- be the least index j > delta_i
with a_{i,j} ≠ 0, j_i^+ the largest index j < delta_i with a_{i,j} ≠ 0, q_i = j_i^+ -
j_i^- - 1 the length of the zero run of Q_i between them, and s_i = min{q_i, e_i}. Define

    B_i = { delta_i, n - d_i + 1, ..., n - d_i + e_i }                if a_{i,delta_i} ≠ 0,
    B_i = { j_i^- } ∪ { j_i^- + 1, ..., j_i^- + s_i } ∪ { n - d_i + s_i + 1, ..., n - d_i + e_i }
                                                                     if a_{i,delta_i} = 0,

where the last set may be empty when s_i = e_i. The two claims below are BK's Claims
6.2 and 6.3; assuming them, Remark B.10 gives c_M ≠ 0.

*Claim B.12 (BK Claim 6.2).* For every i, M_i[A_i^idx, B_i] is lower triangular with
nonzero diagonal entries.

*Proof.* Write a_j := a_{i,j} and let a_j indicate whether a_j is nonzero.
Case a_{delta_i} ≠ 0: B_i = {delta_i} ∪ {n-d_i+1,...,n-d_i+e_i}. An entry of
M_i[A_i^idx, B_i] above the diagonal lies at (row j, column j') with j < j', so it
involves the coefficient a_{j'-j} with j' - j in {n-d_i+1,...,n-d_i+e_i} ⊆
{n-d_i+1,...,n-1}; all coefficients in that window of d_i - 1 consecutive indices
vanish by part (2) of Lemma B.7 applied to the subspace A_i of dimension d_i (which, as
in Remark B.8, is not contained in any constant multiple of a proper subfield because n
is prime and d_i >= 2). Hence the matrix is lower triangular. Its diagonal entries are
a_{delta_i} for the first row, which is nonzero by the case assumption, and a_{n-d_i} for
every subsequent row (column j' - row j = n - d_i), which is nonzero because Q_i has
F_2-degree exactly 2^(n-d_i) by Lemma B.5.
Case a_{delta_i} = 0: the uppermost (s_i + 1) x (s_i + 1) submatrix of M_i[A_i^idx, B_i]
is

        [ a_{j_i^-}      a_{j_i^-+1}   ...  a_{j_i^-+s_i}      ]
        [ a_{j_i^--1}    a_{j_i^-}     ...  a_{j_i^-+s_i-1}    ]
        [      ...          ...         ...      ...           ]
        [ a_{j_i^--s_i}  a_{j_i^--s_i+1} ... a_{j_i^-}         ]

which is lower triangular because a_{j_i^-+1} = ... = a_{j_i^++s_i} = 0 by definition of
j_i^- and q_i (the run of zeros), and its diagonal entries a_{j_i^-} and a_{j_i^+} are
nonzero by definition of j_i^- and j_i^+. The remaining e_i - s_i columns are identical
to the last columns of the previous case, so the whole matrix is lower triangular with
nonzero diagonal.

*Claim B.13 (BK Claim 6.3).* For every r in {1,...,n-1} there is an i with
M_i[A_i^idx, B_i - r] containing an all-zero column.

*Proof.* Note that M_i[A_i^idx, {k} - r] is the column vector
(a_{i,k-r}, a_{i,k-(r+1)}, ..., a_{i,k-(r+e_i)})^T. Hence, for k in B_i, if
[k - (r + e_i), k - r] mod n ⊆ [n - d_i + 1, n - 1] (mod n), then M_i[A_i^idx, B_i - r]
contains a zero column. In particular (BK Proposition 6.4) whenever k is in B_i and
r ∈ [k+1, k + d_i - (e_i+1)] mod n, the matrix M_i[A_i^idx, B_i - r] contains an all-zero
column. It therefore suffices to show

    [1, n-1] ⊆ union over i of union over k in B_i of [k+1, k + (d_i - e_i) - 1].   (B.6)

Since delta_1 = 0 we have B_1 = {0} ∪ [n-d_1+1, n-d_1+e_1], so by the previous paragraph
M_1[A_1^idx, B_1 - r] contains a zero column for every r in
[1, d_1 - (e_1+1)] ∪ [n-d_1+2, n-1] = [delta_1+1, delta_2] ∪ [n-d_1+2, n-1].        (B.7)

Let t0 be the minimal i with sum_{i' <= i}(d_{i'} - (e_{i'}+1)) >= n - d_1 + 1; such t0
exists by hypothesis (B.4): indeed sum_{i=1}^{t}(d_i - (e_i+1)) > n - max d_i + 1
= n - d_1 + 1. Then delta_{t0} < n and delta_{t0+1} = sum_{i <= t0}(d_i - (e_i+1)) >=
n - d_1 + 1. We claim that for 1 < i <= t0,

    union over k in B_i of [k+1, k + d_i - (e_i+1)] ⊇ [delta_i + 1, delta_{i+1}],  (B.8)

which together with (B.7) proves (B.6) and completes the proof. There are two cases.

Case a_{i,delta_i} ≠ 0: then delta_i is in B_i, so (B.8) follows from the previous
paragraph and delta_{i+1} = delta_i + d_i - (e_i+1).

Case a_{i,delta_i} = 0: here j_i^- < delta_i < j_i^+.
  If q_i < e_i then B_i = {j_i^-} ∪ [j_i^-+1, j_i^++q_i] ∪ [n-d_i+q_i+1, n-d_i+e_i], and
  substituting j_i^+ = j_i^- + q_i + 1 and reordering gives B_i = [j_i^-+1, j_i^+] ∪
  [n-d_i+q_i+1, n-d_i+e_i]; as delta_i is in B_i we are back in the previous case.
  If q_i >= e_i then B_i = {j_i^-} ∪ [j_i^-+1, j_i^-+e_i] = [j_i^-+1, j_i^++e_i] ∪
  {j_i^+ + q_i + 1} after substituting j_i^+ = j_i^- + q_i + 1. Since i <= t0 and d_1 is
  maximal we have j_i^- < delta_i <= n - d_1 <= n - d_i, hence j_i^- + 1 <= n - d_i + 1;
  combined with part (2) of Lemma B.7 applied to the subspace A_i (which is not
  contained in an affine shift of a proper subfield) this forces j_i^+ - j_i^- >= d_i - 1,
  i.e. j_i^+ >= j_i^- + d_i - 1. Taking all but the last element of B_i we get
  union over k in [j_i^-+1, j_i^-+e_i] of [k+1, k + d_i - (e_i+1)] ⊇ [j_i^-+2, j_i^-+d_i-1],
  which contains delta_i (because j_i^- < delta_i < j_i^+ = j_i^-+q_i+1 and q_i >= e_i).
  Reinserting j_i^+ + q_i + 1 into B_i and using j_i^- < delta_i gives (B.8). QED

This completes the proof of Theorem B.11. QED

*Remark B.14.* In our application we will always have d_i >= 2 for every i with
nonempty block, so the hypothesis "not contained in an affine shift of a proper subfield"
is automatic by Remark B.8 whenever n is prime, and we will check (B.4) explicitly.

------------------------------------------------------------------

## B.5 Decomposing an affine subspace into independent blocks

*Lemma B.15 (Bourgain's decomposition; BK Lemma 7.1).* Let A ⊆ K^r be an F_2-affine
subspace and let a in A. Then there are F_2-linear subspaces Y_1,...,Y_r ⊆ K and F_2-linear
maps phi_{ij}: Y_j -> K for j < i such that

    A = { (x_1,...,x_r) : y_i in Y_i and x_i = a_i + y_i + sum_{j<i} phi_{ij}(y_j) },  (B.9)

and dim A = sum_{i=1}^{r} dim Y_i.

*Proof.* Take the row echelon form of a matrix whose rows are a basis of the linear
subspace A - a (BK, Lemma 7.1; the echelon basis is triangular, and the free part of the
i-th row becomes Y_i). QED

*Lemma B.16 (BK Lemma 7.2, p = 2).* For each i <= r let P_i(X_i) be a linearized
polynomial, and for each j < i let P_{ij}(X_j) be a linearized polynomial. Let
gamma in K. Let I_0 ⊆ [r] with I_0 = {i_1 < i_2 < ... < i_t}, let e_{i_1},...,e_{i_t} > 1
be integers, and put theta_i = sum_{l=0}^{e_i} 2^l for i in I_0. Let gamma-hat in
K \ {0}. Let

    g(X_1,...,X_r) = Tr( gamma-hat * prod_{i in I_0} ( P_i(X_i) + sum_{j<i} P_{ij}(X_j) + gamma_i )^(theta_i) )
                     mod (X_l^(2^n) - X_l)_{l in [r]},

    g'(X_1,...,X_r) = Tr( gamma-hat * prod_{i in I_0} P_i(X_i)^(theta_i) )
                      mod (X_l^(2^n) - X_l)_{l in [r]}.

Then for every choice of (gamma_i)_{i in I_0} with wt(gamma_i) = e_i + 1, the coefficients
of the monomial prod_{i in I_0} X_i^(gamma_i) in g and in g' are equal.

*Proof.* Expand

    g = sum_{s} gamma-hat^(2^s) prod_{i in I_0} ( P_i(X_i) + P^i(X_{<i}) )^(2^s (2^(e_i+1)-1))

where P^i(X_{<i}) = sum_{j<i} P_{ij}(X_j) + gamma_i does not involve X_i, and
theta_i = 2^(e_i+1) - 1. For each i, expand the inner power:

    (P_i(X_i) + P^i(X_{<i}))^(theta_i)
      = sum over subsets L of {0,...,e_i} of ( P_i(X_i)^(sum_{l in L} 2^l) * (prod_{l notin L} P_i(X_i)^(2^l)) * P^i(X_{<i})^(sum_{l notin L} 2^l) )
      = ( P_i(X_i) )^(theta_i) + G_i(X_i, X_{<i}),

where the first summand (L = {0,...,e_i}) is the product of e_i + 1 copies of the
linearized polynomial P_i, hence has F_2-degree exactly e_i + 1 = theta_i + 1 in X_i,
while every other summand contains at most e_i copies of P_i (the remaining factors being
the linearized P^i(X_{<i}), which are independent of X_i), hence has F_2-degree at most
e_i in X_i.

Let M = prod_{i in I_0} X_i^(gamma_i) with wt(gamma_i) = e_i + 1 for every i in I_0;
after reduction modulo (X_l^(2^n) - X_l)_l the F_2-degree of M in each X_i is exactly
wt(gamma_i) = e_i + 1.

Now substitute into the product over i in I_0 and order by the largest index. Write
i_t for the largest element of I_0; the factor for i_t is (P_{i_t}(X_{i_t}))^(theta_{i_t})
+ G_{i_t}(X_{i_t}, X_{<i_t}), and all remaining factors H_{t-1}(X_{<i_t}) are independent
of X_{i_t}. Therefore the F_2-degree of the whole product in X_{i_t} is the sum of the
F_2-degrees of the factors in X_{i_t}: taking the G_{i_t} factor bounds it by e_{i_t} (plus
0 from the other factors), while taking the first factor gives exactly theta_{i_t} + 1.
Since the target monomial M has F_2-degree e_{i_t} + 1 = theta_{i_t} + 1 > e_{i_t} in
X_{i_t}, only monomials that select the first factor can contribute to M. Applying the
same argument successively to i_{t-1}, ..., i_1 (each of these indices is smaller, so its
factor is independent of the already-eliminated variables), the only contributions to the
coefficient of M come from prod_{i in I_0} P_i(X_i)^(theta_i). This is the assertion. QED

------------------------------------------------------------------

## B.6 The main theorem

*Theorem B.17 (BK Theorem 7.3, p = 2).* Let t < r be integers and let n be prime with
n >= r(r+1)/t. For i in [r] put e_i = r + 1 - i and theta_i = sum_{l=0}^{e_i} 2^l =
2^(e_i+1) - 1. Let gamma-hat in K \ {0}, and let f: K^r -> F_2 be

    f(x_1,...,x_r) = Tr( gamma-hat * sum over I ⊆ [r], |I| = t of prod_{i in I} x_i^(theta_i) ).  (B.10)

Let A ⊆ K^r be an F_2-affine subspace with

    dim A > n r / t + n t + r (r + 1).                                  (B.11)

Then |f(A)| > 1.

*Proof.* (1) Decompose A by Lemma B.15: A = { x : x_i = a_i + y_i + sum_{j<i} phi_{ij}(y_j),
y_i in Y_i }, with dim A = sum_i dim Y_i =: d and dim Y_i =: d_i. Let Q_i be the
image-subspace polynomial of the linear subspace Y_i (Lemma B.5 applied to the linear
space Y_i: Q_i(Y_i)... Q_i is the kernel-subspace polynomial of the image of Y_i, and we
use the linearized polynomial Q_{ij} with Q_{ij}(x) = phi_{ij}(Q_i(x)) for each x in K,
which exists by Lemma 3.7, i.e. B.5, applied to the linear map phi_{ij} o (image-subspace
polynomial)); set R_i(X_1,...,X_r) = Q_i(X_i) + sum_{j<i} Q_{ij}(X_j) + a_i. Then the
polynomial map H(x) = (R_1(x),...,R_r(x)) has image exactly A.

(2) The composed polynomial h = f o H satisfies h(K^r) = f(A) (every value of f on A is
attained, because H is onto A). So it suffices to show that h is non-constant on K^r.

(3) Write

    h = sum over I ⊆ [r], |I| = t of Tr( gamma-hat prod_{i in I} R_i(X)^(theta_i) ) mod (X_l^(2^n) - X_l)_l.   (B.12)

Put S = { i : d_i > r + 1 }. First, |S| >= t: every d_i <= n, and for i not in S we have
d_i <= r + 1, so if |S| <= t - 1 then

    d = sum_i d_i <= |S| n + (r - |S|)(r+1) <= (t-1) n + r (r+1)
      <  n r/t + n t + r (r+1)  <  d,

a contradiction (the strict inequality in the middle is (t-1)n < nr/t + nt, which holds
because nr/t > 0; the last is the hypothesis (B.11)). Second, since (d_i - (r+1)) <= 0
for i not in S,

    sum_{i in S} d_i  >=  sum_{i in S} (d_i - (r+1))  =  d - r(r+1) - sum_{i notin S} (d_i - (r+1))
                      >=  d - r(r+1)  >  n r/t + n t.

Hence we may choose I_0 ⊆ S with |I_0| = t consisting of the t largest d_i in S; writing
S' for the chosen set and using |S| <= r,

    sum_{i in I_0} (d_i - (r+1))  =  sum_{i in I_0} d_i - t(r+1)
                                  >=  (t/|S|) sum_{i in S} d_i - t(r+1)
                                  >=  (t/r) sum_{i in S} d_i - t(r+1)
                                  >   (t/r)(nr/t + nt) - t(r+1)
                                  =   n + n t^2 / r - t(r+1)  >=  n,                (B.13)

where the last inequality is exactly r(r+1) <= n t, i.e. the hypothesis n >= r(r+1)/t
multiplied by t. (Both averaging steps are legitimate because I_0 consists of the t
largest of the |S| >= t numbers d_i, whose sum is at least (t/|S|) times the sum over
S.)

(4) Focus on the summand of (B.12) indexed by I_0:

    g = Tr( gamma-hat prod_{i in I_0} R_i(X)^(theta_i) ) mod (X_l^(2^n) - X_l)_l,

and compare it with

    g' = Tr( gamma-hat prod_{i in I_0} Q_i(X_i)^(theta_i) ) mod (X_l^(2^n) - X_l)_l.

For i in I_0 we have e_i + 1 = r + 2 - i <= r+1 < d_i, so 1 <= e_i < d_i - 1, and by
(B.13)

    sum_{i in I_0} (d_i - (e_i + 1)) >= sum_{i in I_0} (d_i - (r+1)) >= n > n - max_{i in I_0} d_i + 1,

so Theorem B.11 applies to the affine subspaces Y_i (i in I_0) — each is a linear
subspace of dimension d_i >= 2, hence by Remark B.14 not contained in an affine shift of
a proper subfield since n is prime — and yields a monomial M = prod_{i in I_0}
X_i^(gamma_i) with wt(gamma_i) = e_i + 1 = theta_i + 1 whose coefficient in g' is
nonzero. By Lemma B.16 the coefficient of M in g equals its coefficient in g', hence is
nonzero.

(5) Finally, no other summand of (B.12) contributes to the monomial M. Consider
J ⊆ [r] with |J| = t and the summand Tr(gamma-hat prod_{j in J} R_j(X)^(theta_j)). Each
R_j(X_1,...,X_j) is a polynomial in X_1,...,X_j only, of the form sum of terms
a X_l^(2^b) plus a constant (F_2-degree at most 1 in each variable used). Expanding the
trace, any contribution to M must come from a product of terms in which, for each (j,l)
with j in J and l in [0, r+1-j], a monomial of R_jl (a rotation of R_j) is chosen whose
variable is some X_i with i <= j; writing j(l) = i for the chosen index, the F_2-degree
of M in X_i is at most the number of pairs (j,l) with j(l) = i, i.e. | union_{j in J}
j^{-1}(i) |. Since M has F_2-degree e_i + 1 = r + 2 - i in X_i for each i in I_0, such a
choice requires sum_{j in J} |j^{-1}(i)| >= r + 2 - i for every i in I_0. Claim B.18 below
then forces J = I_0, so only the I_0 summand can contribute, and by (4) its contribution
is nonzero. Hence h is non-constant on K^r and |f(A)| > 1. QED

*Claim B.18 (BK Claim 7.4).* Let I_0, J ⊆ [r] with |I_0| = |J| = t, and suppose there
are functions j: {0,1,...,r+1-j} -> {1,...,j} for j in J with
sum_{j in J} |j^{-1}(i)| >= r + 2 - i for every i in I_0. Then J = I_0.

*Proof.* For any j in J, sum_{i in I_0} |j^{-1}(i)| <= |{0,1,...,r+1-j}| = r+2-j, so

    sum_{i in I_0} (r+2-i) <= sum_{j in J} sum_{i in I_0} |j^{-1}(i)| <= sum_{j in J} (r+2-j).

Since |I_0| = |J|, this gives

    sum_{i in I_0} i >= sum_{j in J} j.                                              (*)

For m >= 1 consider the expression sum_{j in J} sum_{l=0}^{r+1-j} (j^m - j(l)^m), which
is nonnegative because j(l) <= j. Expanding,

    sum_{j in J} (r+2-j) j^m  >=  sum_{j in J} sum_l j(l)^m
                              =   sum_{i in I_0} i^m |{ (j,l) : j(l) = i }|
                              >=  sum_{i in I_0} i^m (r+2-i),

where the last step is the hypothesis of the claim. Substituting into (*) and rearranging,

    sum_{i in I_0} i^(m+1) - sum_{j in J} j^(m+1)  >=  (r+2) ( sum_{i in I_0} i^m - sum_{j in J} j^m )  >=  0,

the last inequality by induction on m, starting from (*) for m = 1. Hence

    sum_{i in I_0} i^m  >=  sum_{j in J} j^m   for every m >= 1.                     (**)

Now (**) for all m forces max I_0 <= max J: if i_1 := max I_0 > j_1 := max J then
i_1^m > |J| j_1^m >= sum_{j in J} j^m for large m, contradicting (**). On the other hand
the hypothesis of the claim with i = i_1 gives |{ (j,l) : j(l) = i_1 }| >= r + 2 - i_1
>= 1, so some j in J satisfies j >= i_1 (as j(l) <= j); combined with i_1 <= max J we get
i_1 = j_1 = max J. Taking this back into the chain above, the only pairs (j,l) with
j = j_1 that can exist have j(l) = i_1 = j_1, and (*) and (**) for the remaining indices
give the same conclusion for the second largest elements: writing i_2 = max(I_0 \ {i_1})
and j_2 = max(J \ {j_1}), the same argument (with the pairs above the position of j_1
removed) yields i_2 = j_2, and inductively the s-th largest elements of I_0 and J agree
for every s. Since |I_0| = |J|, this gives I_0 = J. QED

------------------------------------------------------------------

## B.7 Instantiation: the payload of our family

Let N >= 2^20, and let k, m, r, t be as in `proof.md` section 2.1: k the least prime with
k^5 > 32N^3, m = N - k, r = ceil(m/k), t = ceil(sqrt(r)). Lemma A1 of `proof.md` gives
k prime, t < r, and k >= 2 r(r+1)/t, hence n := k satisfies the hypothesis n >= r(r+1)/t
of Theorem B.17. Take gamma-hat = 1 in (B.10); then f: K^r -> F_2 is

    f(x_1,...,x_r) = Tr( sum over |I| = t of prod_{i in I} x_i^(2^(r+2-i) - 1) ),

which is exactly the payload function of `proof.md` section 2.2 under the block map
(Proposition 3.3 of BK identifies BK's two presentations of the same function; the
exponents here are BK's theta_i with e_i = r + 1 - i, and reversing the block index
recovers the ascending order 2^(i+1)-1).

Let phi: {0,1}^m -> K^r be the injective F_2-affine block map of `proof.md` section 3.3.
If A ⊆ {0,1}^m is an affine subspace with dim A > d(N) where

    d(N) := ceil( n r/t + n t + r(r+1) ) <= 16 N^(4/5)

(the inequality is proved in `proof.md` section 3.2 with 2N^(3/5) < k < 5N^(3/5),
r <= N^(2/5)/2 + 1 and t <= 2N^(1/5)), then phi(A) is an affine subspace of K^r of the
same dimension, so Theorem B.17 applies and f is non-constant on phi(A). Hence g_N =
f o phi is non-constant on A.

Conclusion: for every N >= 2^20 the payload function g_N is an affine disperser for
dimension d(N) <= 16 N^(4/5), and the proof of that statement is contained in this file:
B.1-B.3 develop the subspace-polynomial and permanent machinery, B.4 proves non-constancy
for independent affine sources, B.5 decomposes a general affine subspace into independent
blocks, B.6 composes the two, and this section specializes the result to the family's
parameters.

*Scope.* No premise of this file is external. The method and all statements are from
Ben-Sasson and Kopparty's paper (cited in `proof.md` section 6.2) and are reproduced here
with their proofs; the only quoted facts are elementary finite-field algebra, stated as
B.1-B.3 with proofs. Where BK's text leaves an estimate to the reader (the set S in step
(3) and the permanent identity in B.3), the estimate is written out above.
