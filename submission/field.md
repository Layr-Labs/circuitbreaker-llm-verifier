# Field construction, uniform evaluator, and NP membership

This file is part of the supplied evidence for the joint theorem declared in
`bound.json`. It defines, with full proofs, one fixed deterministic
polynomial-time procedure that constructs the finite fields used by the family,
and one fixed polynomial-time evaluator computing the family at every input
length. No external premise is used beyond elementary finite-field arithmetic,
which is proved here.

Throughout, n >= 1 is the actual input length. Write

    k(n) := the least prime strictly greater than 2*n^(3/5),

and let d(n) := floor((4251/1000)*n^(4/5)) + 1 (the disperser dimension
proved in `disperser.md`; the constant 4251/1000 < 4.251 is established
there). By Bertrand's postulate (for every
integer N > 1 there is a prime p with N < p < 2N), applied to
N = ceil(2*n^(3/5)), we have k(n) < 2*ceil(2*n^(3/5)) <= 4*n^(3/5) + 2. Since
2*n^(3/5) >= 2, k(n) >= 3, so k(n) is an odd prime. The quantities

    r(n) := ceil(n / k(n)),   t(n) := ceil(sqrt(r(n)))

are used by the disperser construction in `disperser.md`.

## 1. Deterministic factoring over F_2 (Berlekamp)

**Lemma 1.1 (irreducibility test and factoring).** There is a deterministic
algorithm that, given a monic squarefree h in F_2[X] of degree D >= 1,
outputs the factorization of h into monic irreducible polynomials, in time
polynomial in D. (Only squarefree inputs are needed below: the two
polynomials ever factored are Phi_k and f^(i-1)(X^k), both proved squarefree
in Lemma 1.2.)

*Proof.* The algorithm and its correctness follow.

*Squarefree input.* h is monic squarefree of degree D: gcd(h, h') = 1, so
h has no repeated roots and is the product of its distinct monic
irreducible factors.

*Berlekamp's subalgebra.* For a polynomial g with g = sum_i g_i X^i
(g_i in F_2), freshman's dream gives g^2 = sum_{g_i = 1} X^(2i). Hence the
map S(g) := g^2 mod h is F_2-linear on the D-dimensional algebra F_2[X]/(h).
Its matrix in the basis 1, X, ..., X^(D-1) has as its j-th column the
coefficient vector of S(X^j) = X^(2j) mod h, computed by one polynomial
squaring and reduction per column, O(D^2) field operations each, O(D^3) in
total. Let B := { g : g^2 = g mod h } = ker(S - I), computed by
Gaussian elimination over F_2 (deterministic, polynomial). Let
h = f_1 * ... * f_u be the factorization into monic irreducibles (u
unknown). By the Chinese remainder
theorem, F_2[X]/(h) is isomorphic to the product of the fields
F_2[X]/(f_j); an element g satisfies g^2 = g if and only if each component
g_j in the field F_2[X]/(f_j) satisfies g_j^2 = g_j, i.e. g_j is a root of
Z^2 - Z, whose only roots in a field of characteristic 2 are 0 and 1. Hence
B has exactly 2^u elements, so dim B = u, and u = 1 (h irreducible) holds if
and only if B = {0, 1}.

*Splitting.* If B contains a nonconstant g, then g is congruent to 0 or 1
modulo each f_j and to both values for at least one index each (g not congruent
to a constant mod h). Therefore h divides g*(g+1), and

    h = gcd(h, g) * gcd(h, g + 1)

with both factors monic of degree at least 1 and less than D: every f_j divides
exactly one of g, g+1; at least one f_j divides g (otherwise g = 0 mod h) and at
least one f_j divides g+1 (otherwise g = 1 mod h). Both gcd's are computed by
the Euclidean algorithm. Recursing on the two proper monic squarefree factors
until every factor has B = {0,1} yields all monic irreducible factors; the
recursion depth is at most D and each level costs O(D^3). This proves the
lemma. qed

**Lemma 1.2 (cyclotomic factors).** Let k >= 3 be an odd prime and let
m := ord_k(2), the multiplicative order of 2 modulo k (the least j >= 1 with
2^j = 1 mod k; it exists since gcd(2, k) = 1 and divides k - 1). Let
Phi_k(X) := (X^k - 1)/(X - 1) in F_2[X], of degree k - 1. Then:

1. Every irreducible factor of Phi_k over F_2 has degree exactly m, and its
   roots (in the algebraic closure of F_2) are primitive k-th roots of unity.
2. For i >= 1 let f^(i) be any monic irreducible factor of f^(i-1)(X^k), where
   f^(1) is any monic irreducible factor of Phi_k. Then every root of f^(i)
   is a primitive k^i-th root of unity, and f^(i) is well defined whenever
   f^(i-1)(X^k) is a nonconstant squarefree polynomial.
3. f^(i-1)(X^k) is nonconstant and squarefree for every i >= 2.

*Proof.* (1) X^k - 1 is squarefree since gcd(X^k - 1, k*X^(k-1)) = 1 in
characteristic 2 (k odd). Its k roots form the group of k-th roots of unity;
the roots of Phi_k are those of order exactly k (the primitive ones), k - 1 in
number, since the only root of X - 1 is 1. Let zeta be a primitive k-th root and
F its minimal polynomial over F_2, of degree e. The field F_2(zeta) is the
finite field F_{2^e}, which contains zeta if and only if zeta^(2^e) = zeta,
i.e. zeta^(2^e - 1) = 1, i.e. k divides 2^e - 1, i.e. m divides e.
Conversely k divides 2^m - 1, so zeta lies in F_{2^m}, whence
F_2(zeta) = F_{2^e} is contained in F_{2^m}, i.e. e divides m. With m
dividing e, e = m, proving the degree claim; F's roots are the Frobenius
conjugates zeta, zeta^2, zeta^(2^2), ..., all primitive k-th roots of
unity, so all roots of F are such.

(2) Let zeta be a root of f^(i-1)(X^k), i.e. zeta^k = w for some root w of
f^(i-1). By induction every root of f^(i-1) is a primitive k^(i-1)-th root
of unity, so w has order exactly k^(i-1). Write the order of zeta as
k^a * v with gcd(v, k) = 1 and a >= 0. The order of zeta^k is
k^(a-1) * v when a >= 1, and is v itself (coprime to
k) when a = 0; matching order exactly k^(i-1) (a power of the prime k)
forces v = 1 and a = i >= 1. So zeta has order exactly k^i. (Only this
forward statement is claimed and used: f^(i) is a nonconstant polynomial
whose roots are primitive k^i-th roots of unity; f^(i) need not contain
all of them in its Frobenius orbit.) This proves the roots claim.

(3) f^(i-1) has only simple roots (it is irreducible over F_2, hence
squarefree), all nonzero (they are roots of unity). The roots of f^(i-1)(X^k)
are then also simple: for a root zeta with zeta^k = w a simple root of
f^(i-1), the derivative of f^(i-1)(X^k) at zeta is
f^(i-1)'(zeta^k) * k * zeta^(k-1) = f^(i-1)'(w) * zeta^(k-1) != 0
because k = 1 in F_2, zeta != 0, and w is a simple root of f^(i-1). So
f^(i-1)(X^k) is squarefree, and it is nonconstant by (2). qed

## 2. Constructing an irreducible polynomial of odd prime degree k

Fix an odd prime k. All steps below are deterministic and polynomial in k; any
free choice (a kernel element, a factor, an index) is resolved by taking the
lexicographically first option, so the procedure is a single fixed algorithm.

**Step A.** Compute m := ord_k(2) by repeated doubling modulo k.

**Step B.** Set f^(1) := the lexicographically first monic irreducible factor
of Phi_k (Lemma 1.1, Lemma 1.2(1)); it has degree m. Compute the largest s with
k^s dividing 2^m - 1 by trial division (s >= 1 since k | 2^m - 1, and
s <= m). For i = 2, ..., s set f^(i) := the lexicographically first monic
irreducible factor of f^(i-1)(X^k) (Lemmas 1.1, 1.2). Write m_s for the degree
of f := f^(s) and set

    K := F_2[X]/(f),  alpha := the class of X in K.

By Lemma 1.2(2), alpha has multiplicative order exactly k^s. Note m_s divides m
(the degree of any f^(i) is ord_{k^i}(2), which divides ord_{k^s}(2) = m_s for
i <= s, and in particular every degree divides m because k^i | 2^m - 1 for
i <= s); also m_s <= m.

**Lemma 2.1 (alpha is a k-th nonresidue).** alpha is not a k-th power in the
multiplicative group K*.

*Proof.* K has 2^(m_s) elements, and alpha in K* has order k^s, so
k^s | 2^(m_s) - 1. Since m_s | m, the integer 2^(m_s) - 1 divides
2^m - 1 = ell * k^s with gcd(ell, k) = 1, so the k-adic valuation of
2^(m_s) - 1 is exactly s (it is at least s and at most that of 2^m - 1). If
alpha = beta^k for some beta in K*, then the order of beta^k, which is
order(beta)/gcd(order(beta), k), equals k^s, forcing order(beta) = k^(s+1)
(writing order(beta) = k^a*v with gcd(v,k)=1, beta^k has order k^(a-1)*v,
so v = 1 and a = s+1). Then k^(s+1) divides |K*| = 2^(m_s) - 1, contradicting
the valuation bound. qed

**Lemma 2.2 (binomial irreducibility, prime exponent).** Let E be any field, k
an odd prime, a in E*, and suppose a is not a k-th power in E*. Then
X^k - a is irreducible in E[X].

*Proof.* Suppose X^k - a = g * h with g monic irreducible of degree e,
1 <= e <= k - 1, and let beta be a root of g, so beta^k = a and
[E(beta):E] = e. Let beta = beta_1, ..., beta_e be the roots of g (the
conjugates of beta over E) and c := beta_1 * ... * beta_e the constant term
of g up to sign; c lies in E and

    c^k = (beta_1 * ... * beta_e)^k = (beta^k)_1 * ... * (beta^k)_e = a^e,

because each conjugate of beta^k = a is a itself. Since k is prime and e < k,
gcd(e, k) = 1, so there are integers u > 0, v > 0 with u*e - v*k = 1
(Bezout, adjusted by adding multiples of e*k). Then

    a = a^(u*e - v*k) = (a^e)^u * (a^v)^(-k) = (c^k)^u * ((a^v)^(-1))^k
      = (c^u / a^v)^k,

exhibiting a as a k-th power of the element c^u / a^v of E*, a
contradiction. qed

**Step C.** Define the field extension

    E := K[Y]/(Y^k - alpha),   beta := the class of Y,

which is a field by Lemmas 2.1 and 2.2, of size 2^(m_s * k). Since k divides
m_s * k, E contains a unique subfield H of size 2^k, and Gal(E/H) is cyclic of
order m_s generated by the Frobenius x -> x^(2^k). Define

    gamma := T(beta),  T(x) := sum_{j=0}^{m_s - 1} x^(2^(k*j)),

the trace from E to H. So gamma lies in H.

**Lemma 2.3 (gamma has degree k over F_2).** gamma is an element of H not
lying in K; consequently gamma has degree exactly k over F_2, and its minimal
polynomial g(X) = prod_{i=0}^{k-1} (X - gamma^(2^i)) is a monic irreducible
polynomial of degree k over F_2.

*Proof.* For 0 <= j <= m_s - 1 write 2^(k*j) = nu_j * k + rho_j with
0 <= rho_j < k. Since 2^k = 2 mod k (Fermat), rho_j = 2^j mod k. Because
j < m_s <= m = ord_k(2), the residues rho_0, ..., rho_(m_s - 1) are distinct
and nonzero. Using beta^k = alpha,

    beta^(2^(k*j)) = beta^(nu_j*k + rho_j) = alpha^(nu_j) * beta^(rho_j),

so gamma = P(beta) where P(X) := sum_{j=0}^{m_s-1} alpha^(nu_j) * X^(rho_j)
is a nonzero polynomial over K of degree at most k - 1 (nonzero because its
monomials X^(rho_j) are distinct with nonzero coefficients alpha^(nu_j)).
If gamma were in K, then beta would be a root of the nonzero polynomial
P(X) - gamma in K[X] of degree at most k - 1, contradicting
[K(beta):K] = k. Hence gamma not in K, so in particular gamma is not in
F_2. Now gamma lies in H, a field of 2^k elements with k prime; the subfields
of H are F_{2^c} for c dividing k, i.e. F_2 and H only. So gamma has degree
exactly k over F_2. The elements gamma, gamma^2, ..., gamma^(2^(k-1)) are the
k distinct Frobenius conjugates of gamma, so g(X) as defined is its minimal
polynomial: monic, irreducible, of degree k, with all coefficients fixed by
Frobenius and hence in F_2. qed

**Step D (output).** Compute gamma in E (m_s exponentiations by repeated
squaring, each O(m_s * k) squarings in E), then its k conjugates gamma^(2^i),
then expand g(X) = prod_i (X - gamma^(2^i)) inside E[X]; by Lemma 2.3 the
coefficients of g lie in F_2 (each is fixed by the Frobenius automorphism
x -> x^2 of E). Output the coefficient list of g. Every arithmetic step is a
polynomial number of operations in fields of size at most 2^(m_s*k),
2^(m_s*k) <= 2^(k^2), each operation costing poly(k^2) bit operations.

**Proposition 2.4 (fixed deterministic irreducible polynomial).** There is a
single fixed deterministic algorithm IRRED which, on input an odd prime k,
outputs a monic irreducible polynomial q_k in F_2[X] of degree exactly k, in
time polynomial in k. In particular, on input k = k(n) this yields the field
K_n := F_2[Z]/(q_k) of 2^k elements, with the fixed polynomial basis
1, Z, ..., Z^(k-1).

*Proof.* Steps A-D above, with all free choices resolved lexicographically.
Deterministicity, polynomial running time, and correctness (monic, irreducible,
degree k) follow from Lemmas 1.1, 1.2, 2.1, 2.2, 2.3. The primality of k and
the value k(n) itself are computed by trial division, polynomial in n. qed

## 3. The family, its evaluator, and NP membership

**Definition of the family (one fixed total language).** Define f_0 := 0 on the
empty input. For n >= 1:

- k := k(n), r := r(n), t := t(n) as above (k odd prime, r = ceil(n/k),
  t = ceil(sqrt(r)) >= 1; also t <= r because r >= 1);
- run IRRED(k) to obtain q_k and the field K_n = F_2[Z]/(q_k) with basis
  1, Z, ..., Z^(k-1);
- pad the input x with r*k - n zeros and split it into r consecutive blocks
  B_1, ..., B_r of k bits each; let phi_i(x) in K_n be the field element whose
  coordinate vector in the basis is (first bit of B_i = coefficient of Z^0,
  ..., last bit = coefficient of Z^(k-1)); the map
  phi = (phi_1, ..., phi_r): F_2^n -> (K_n)^r is F_2-linear and injective
  (the padding is fixed and the blocks partition all r*k >= n bits; a
  restriction of a bijection F_2^(r*k) -> (K_n)^r);
- compute y_i := phi_i(x)^(2^(i+1) - 1) in K_n for 1 <= i <= r (binary
  powering: the exponent has i + 1 bits, so O(i) multiplications);
- compute the elementary symmetric polynomial
  Sym^t_r(y_1, ..., y_r) = sum over t-element subsets S of {1..r} of
  prod_{i in S} y_i, by the division-free recurrence: initialize E_0 = 1 and
  E_j = 0 for 1 <= j <= t; after processing y_i (for i = 1, ..., r), update
  E_j := E_j + y_i * E_(j-1) for j = min(i, t) down to 1. The invariant
  E_j = Sym^j_i(y_1, ..., y_i) follows by induction on i, splitting the
  t-subsets of {1..i} into those containing i and those not;
- output pi(E_t), where pi: K_n -> F_2 returns the coefficient of Z^0
  (a fixed nontrivial F_2-linear map).

Set f_n(x) := the above value, and define the language
L := { x in {0,1}* : f_{|x|}(x) = 1 }.

**Proposition 3.1 (polynomial-time evaluator; P and NP membership).** There is
a single deterministic algorithm A which, on input (n, x) with |x| = n, outputs
f_n(x), in time polynomial in n. Consequently L is in P. Moreover L is in NP
via the fixed verifier V(x, w) := A(|x|, x) with witness bound W(|x|) = 0 and
empty witness w: for every x, x in L iff there exists w, |w| <= W(|x|), with
V(x, w) = 1; both directions are immediate from the definition. The witness
bits are never inputs to the circuits considered in `proof.md`.

*Proof.* The steps above are a fixed deterministic procedure: the prime search
is polynomial in n: candidates are the odd integers from 3 up to
4*n^(3/5) + 2, at most 2*n^(3/5) + 2 of them, each tested by trial division
by all integers up to its square root, i.e. at most 2*n^(3/10) divisions per
candidate;
IRRED(k) is polynomial in k = O(n^(3/5)) by Proposition 2.4; each field
operation in K_n is a polynomial multiplication or division modulo q_k,
O(k^2) bit operations (or better); the powering y_i uses O(r) = O(n^(2/5))
multiplications per block, O(r^2) = O(n^(4/5)) in total across blocks, each
O(k^2) = O(n^(6/5)) bit operations; the symmetric-polynomial recurrence uses
O(r*t) = O(n^(3/5)) field operations. All bounds are polynomial in n. The
evaluator is total on every input length n >= 0 (f_0 = 0 by definition), so
A computes the characteristic function of L in polynomial time, L is in P,
and the empty-witness verifier above establishes both directions of the NP
condition. qed

**Remark (no advice).** Every quantity above (k, r, t, q_k, the basis, the
padding) is a deterministic polynomial-time function of n alone. There is no
per-length advice, no truth-table oracle, and no nonconstructive choice: the
same algorithm and the same language L work at all lengths.

**Remark (why k prime).** `disperser.md` needs the field K_n = F_(2^k) to have
no proper subfield other than F_2; this holds exactly when k is prime, and is
used in the disperser proof (the hypothesis "no affine space is contained in an
affine shift of a proper subfield" then holds automatically for spaces of
dimension at least 2). The primality of k(n) also makes Step B's cyclotomic
argument (Lemma 1.2) apply with exponent k.
