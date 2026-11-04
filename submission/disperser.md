# The explicit affine disperser: construction, quoted theorem, and what is proved here

This file supplies the family's key property. Section 1 fixes the notion used by
rules.md, states the construction, and states the quoted theorem. Sections 2-4
reproduce the parts of the source argument that are short and self-contained, and
Section 5 is an explicit ledger of what is quoted rather than re-proved, so that a
reviewer can see exactly where the boundary of the supplied evidence lies.

Source: E. Ben-Sasson and S. Kopparty, *Affine dispersers from subspace
polynomials*, SIAM J. Comput. 41(4):880-914, 2012; author preprint dated
12 March 2010 (file `affine-disperser.pdf`), Theorem 2.2 with equation (4),
with proofs in Sections 3, 6 and 7. Below, "BK" refers to that preprint and all
section, lemma and theorem numbers are BK's.

## 1. Definition, construction, quoted theorem

**Definition 1.1 (affine disperser, BK Definition 2.1; rules.md).** A function
f : F_2^m -> F_2 is an *affine disperser for dimension d* if it is non-constant on
every affine subspace of F_2^m of dimension d. Equivalently (and this is the form
used in `framework.md`), for every affine subspace A with dim A >= d there are
x, y in A with f(x) != f(y): a subspace of dimension greater than d contains one of
dimension d.

**Construction (BK, Theorem 2.2 and equation (4)).** Let m be the input length and
let k be the least prime strictly larger than 2 m^(3/5). Put r = ceil(m/k) and
t = ceil(sqrt(r)). Let F = F_{2^k}, so that F is a field with 2^k elements and F_2
is its prime field. Let phi : F_2^m -> F^r, y -> (phi_1(y), ..., phi_r(y)) be an
injective F_2-linear map (split the m bits into r consecutive blocks of k bits,
padding the last block with zeros: this is injective and linear). Let
pi : F -> F_2 be a nontrivial F_2-linear map (the constant-coefficient
functional). Define

    g(y) = pi( Sym_r^t( phi_1(y)^(1+2), phi_2(y)^(1+2+2^2), ..., phi_r(y)^(1+2+...+2^r) ) ),  (BK-4)

where Sym_r^t is the t-th elementary symmetric polynomial in r variables. Since
1 + 2 + ... + 2^i = 2^(i+1) - 1, this is the function used in `proof.md` with
y_i = phi_i(y)^(2^(i+1)-1).

**Quoted theorem (BK Theorem 2.2, case p = 2).** With the parameters above, for
every affine subspace A of F_2^m with dim A > 6 m^(4/5), one has |g(A)| > 1; that
is, g is an affine disperser for dimension greater than 6 m^(4/5).

g is computable in time polynomial in m: the symmetric polynomial is computed in
characteristic two by the division-free recurrence of `proof.md` Section 2.3 in
O(r t) field operations, and powers by binary exponentiation (BK make the same
remark, via the Newton-Girard identities).

## 2. The subspace-polynomial toolkit

The proof represents affine subspaces as images of linearised polynomials. This
section records the definitions and the one lemma whose proof is short enough to
reproduce; the others are quoted in Section 5.

A polynomial P in F_{2^k}[X] is *F_2-linearised* if P(X) = sum_i a_i X^(2^i) with
a_i in F_{2^k}; this is equivalent to P(beta b + gamma c) = beta P(b) + gamma P(c)
for all b, c in F_{2^k} and beta, gamma in F_2. It is *affine linearised* if
P(X) = P_0(X) + a with P_0 linearised.

**Lemma 2.1 (BK Lemma 3.4).** For every F_2-linear isomorphism varphi : F_2^k -> F_{2^k}
and every affine transformation T : F_{2^k} -> F_{2^k} there is a unique affine
linearised polynomial P_T with P_T(varphi(b)) = T(varphi(b)) for all b in F_2^k.

**Lemma 2.2 (BK Lemma 3.6).** If L is an affine subspace of F_{2^k} of dimension d
then its *kernel-subspace polynomial* P_L(X) = prod_{alpha in L} (X - alpha) is a
monic affine linearised polynomial of degree 2^d, and it is linearised exactly when
L is a linear space.

**Lemma 2.3 (BK Lemma 3.7).** If L is an affine subspace of F_{2^k} of dimension d
then there is a monic affine linearised polynomial Q_L, the *image-subspace
polynomial* of L, of degree 2^(k-d), with L = Q_L(F_{2^k}) = {Q_L(c) : c in F_{2^k}}.
Moreover P_L(Q_L(X)) = Q_L(P_L(X)) = X^(2^k) - X.

*Proof (BK's, reproduced).* Let L' = P_L(F_{2^k}) be the image of P_L and define
Q_L = P_{L'}, the kernel-subspace polynomial of L'. Then Q_L(P_L(X)) is a monic
polynomial of degree 2^k that vanishes on all of F_{2^k}, hence equals
X^(2^k) - X. Therefore P_L(Q_L(P_L(X))) = P_L(X^(2^k) - X) = P_L(X)^(2^k) - P_L(X).
So the polynomial G(Y) = P_L(Q_L(Y)) - (Y^(2^k) - Y) satisfies G(P_L(X)) = 0; since
deg(G(P_L(X))) = deg G * deg P_L, this forces G = 0. Hence P_L(Q_L(Y)) = 0 for every
Y, that is, the image of Q_L is contained in L, and by dimension counting the image
equals L. []

Two consequences are used constantly: an affine subspace of dimension d is exactly
the image of a monic affine linearised polynomial of degree 2^(k-d), and the
coefficients of these polynomials are highly structured (Section 3.2 of BK).

## 3. From the non-parameterised theorem to the quoted one

BK prove a general theorem and then specialise. The specialisation is short and is
reproduced here, because it is what fixes the constant 6 in the dimension bound.

**Quoted theorem (BK Theorem 7.3, non-parameterised version).** Let t < r be
integers, let k be prime with k >= r(r+1)/t, put e_i = r + 1 - i and
alpha_i = 1 + 2 + ... + 2^(e_i), and let mu be a nonzero element of F_{2^k}. Define

    g_7(x_1,...,x_r) = Tr( mu * sum over t-subsets I of {1..r} of prod_{i in I} x_i^(alpha_i) ).

Then for every F_2-affine subspace A of F_{2^k}^r with dim A > kr/t + kt + r(r+1)
one has |g_7(A)| > 1. (In BK the three terms of the threshold are written with k in
place of the abstract prime n; the middle and last terms are k*t and r(r+1).)

**Lemma 3.1 (BK's proof that Theorem 7.3 implies Theorem 2.2).** With
k = least prime > 2 m^(3/5), r = ceil(m/k), t = ceil(sqrt(r)):

  (a) the hypotheses of Theorem 7.3 hold, i.e. k is prime and k >= r(r+1)/t;
  (b) kr/t + kt + r(r+1) < 6 m^(4/5) for all sufficiently large m;
  (c) hence Theorem 7.3 with dim A > 6 m^(4/5) applies, and Proposition 3.3 of BK
      identifies g_7 with the function (BK-4) up to renaming the variables.

*Proof.* (a) For m >= 4: k < 4m^(3/5) by Bertrand's postulate, so r = ceil(m/k) >=
m^(2/5)/4 and r(r+1)/t <= sqrt(r)(r+1) <= m^(3/5) < k, as BK observe. (Only m >= 512
is used below.)

(b) Write m for the input length. Each of the three terms is estimated separately.
Because t = ceil(sqrt(r)) >= sqrt(r) and r >= m/k, we have t >= sqrt(m/k), so
kr/t <= k sqrt(r) <= k sqrt(m/k + 1) = sqrt(k m + k^2) <= sqrt(k m) + k
<= sqrt(4 m^(8/5)) + 4 m^(3/5) = 2 m^(4/5) + 4 m^(3/5), using k < 4 m^(3/5).
Next, r = ceil(m/k) <= m/k + 1 < m^(2/5)/2 + 1 and therefore
t <= sqrt(r) + 1 <= m^(1/5)/sqrt(2) + 2, so kt <= 2 sqrt(2) m^(4/5) + 8 m^(3/5).
Finally r(r+1) <= (m^(2/5)/2 + 1)(m^(2/5)/2 + 2) <= m^(4/5)/4 + 3 m^(2/5)/2 + 2.
Adding, the threshold is at most (2 + 2 sqrt(2) + 1/4) m^(4/5) + O(m^(3/5)), and
2 + 2 sqrt(2) + 1/4 < 2 + 17/6 + 1/4 < 6 (since 2 sqrt(2) < 17/6). To make the range
explicit, bound the lower-order terms crudely: kr/t <= sqrt(km) + k <= 2m^(4/5) + 4m^(3/5),
kt <= sqrt(km) + 2k <= 2m^(4/5) + 8m^(3/5) and r(r+1) <= m^(4/5)/4 + (3/2)m^(2/5) + 2, so
the threshold is at most 4.25 m^(4/5) + 12 m^(3/5) + 1.5 m^(2/5) + 2. The last three terms
are below 1.75 m^(4/5) exactly when 1.75 u^2 - 12 u - 1.5 > 0 with u = m^(1/5), which holds
for u >= 7 by monotonicity (its value at 7 is 1/4), that is for m >= 16807. Hence the
threshold is below 6 m^(4/5) for every m >= 16807, and `check_cases.py` verifies the same
integer inequality exactly for every m with 512 <= m < 16807. BK's own estimates are
sharper, and the constant 6 is theirs.

(c) is immediate from Theorem 7.3 and (b). []

## 4. The skeleton of BK's proof of Theorem 7.3

The argument is reproduced here at the level of its structure, so that a reader can
see what the quoted lemmas are used for; the lemmas themselves are listed in
Section 5.

1. *Decomposition.* By BK Lemma 7.1 (Bourgain's decomposition), for any affine
   subspace A of F_{2^k}^r and any gamma in A there are linear spaces Y_1..Y_r and
   linear maps sigma_ij : Y_j -> F_{2^k} with
   A = { (x_1..x_r) : exists y_i in Y_i with x_i = gamma_i + y_i + sum_{j<i} sigma_ij(y_j) }
   and dim A = sum_i dim Y_i. This is just the echelon form of a basis of the
   underlying linear subspace, and we reproduce nothing else about it.
2. *A polynomial map onto A.* Let Q_i be the image-subspace polynomial of Y_i
   (Lemma 2.3 above) and let Q_ij be the linearised polynomial with
   Q_ij(x) = sigma_ij(Q_j(x)) (BK Lemma 3.7). Then R_i = Q_i(X_i) + sum_{j<i} Q_ij(X_j) + gamma_i
   defines a map H with H(F_{2^k}^r) = A exactly.
3. *Reduction to a monomial.* h = g_7 o H is a polynomial in the X_i modulo
   X_i^(2^k) - X_i, and h(F_{2^k}^r) = g_7(A) by Proposition 3.1 of BK. So it
   suffices to exhibit a positive-degree monomial with a nonzero coefficient in h.
4. *High-entropy blocks.* With d_i = dim Y_i and d = dim A, let S be the set of
   blocks with d_i > r+1; then |S| >= t and sum_{i in S} d_i >= d - r(r+1), so
   averaging over the t-subsets of S some I_0 carries
   sum_{i in I_0} d_i >= (t/|S|) sum_{i in S} d_i >= (t/r)(d - r(r+1)) > k + k t^2/r,
   which is the hypothesis of BK Theorem 6.1 in the form used there. On that subset,
   BK's Theorem 6.1 (disperser for independent affine sources, algebraic form)
   produces a monomial in the I_0-summand of h whose coefficient is nonzero.
5. *Excluding the other summands.* BK Lemma 7.2 shows that for any other t-subset I
   and any monomial with the p-adic digit weights used in step 4, the coefficients
   of the corresponding monomial agree with those of a product of linearised
   polynomials in disjoint blocks, so the I_0-monomial cannot be produced by any
   other summand. Hence it survives in h with a nonzero coefficient.

## 5. Ledger: proved here versus quoted

Proved or reproduced in this file: Definition 1.1; the construction (BK-4); the
computation of the elementary symmetric polynomial; BK Lemma 3.7 (with BK's proof);
the specialisation Lemma 3.1 (i.e. BK's own reduction of Theorem 2.2 to Theorem
7.3, with explicit estimates); the structure of the proof of Theorem 7.3
(Section 4).

Quoted by number, not re-proved here:

* BK Theorem 2.2 (the statement used as the premise) and BK Theorem 7.3.
* BK Lemma 3.4, Lemma 3.6, Proposition 3.1, Proposition 3.3 (subspace polynomials
  and their coefficients).
* BK Lemma 3.10 (structural lemma on coefficients of image-subspace polynomials),
  Lemma 3.13 (the coefficient formula, equation (14)), Remark 3.14.
* BK Theorem 6.1 (independent affine sources, algebraic form), Claim 6.2, Claim
  6.3, Proposition 6.4, Lemma 7.1, Lemma 7.2.

The residual premise of this submission is therefore exactly BK Theorem 2.2: the
construction is exhibited, its computation is polynomial, its parameter choices are
matched to BK's, and the reduction that gives the constant 6 is reproduced, but the
internal p-adic weight argument of BK Sections 3 and 6 (about twelve pages of
algebra) is quoted, not copied. Anything that depends on that argument inherits the
status of a published theorem rather than of evidence reproduced in this
submission; the submission note says so explicitly.
