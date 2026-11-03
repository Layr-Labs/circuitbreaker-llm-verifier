# (7/2, 1): a full B2 lower bound of 3.5n - o(n) for one fixed uniform family in P

## 0. What is declared

This submission declares the finite joint pair

    (c, p) = (7/2, 1)

for ONE fixed total Boolean family (f_N) that is uniformly computable in
polynomial time, hence in NP with the empty witness, and proves

    C_B2(f_N) >= (7/2) N - 67 N^(4/5)      for every integer N >= 2^10,

in the full Li-Yang B2 model at actual input length N. Consequently

    liminf_{N -> infinity} C_B2(f_N) / N >= 7/2,

which is the declared finite projection. The exponent is exactly 1 and the
coefficient is exactly 7/2 = 3.5, the leading coefficient of the best published
B2 lower bound for an explicit function (Mishra-Santhanam, ECCC TR26-260,
Theorem 1.1: 3.5n - O(log n); Theorem 4.1: 3.5n - 7d for every affine disperser
of dimension d).

No new circuit theorem is claimed. This is a complete, self-contained exposition
of the current-record class theorem (reproduced in `framework.md`, every case and
every arithmetic step included) composed with the explicit family of Ben-Sasson
and Kopparty (2012) whose affine-disperser property is quoted in `disperser.md`.
The contribution of the submission is the assembled end-to-end proof at actual
input length, the explicit field construction that removes the need for a
deterministic irreducible-polynomial algorithm, and the exact remainder
bookkeeping with an explicit N_epsilon.

## 1. The family

Fix the following total function f_N : {0,1}^N -> {0,1} for every integer N >= 1.
Everything is deterministic; there is no advice, no oracle, no search over
exponentially many candidates, and no dependence on unproved existence.

**Step 1 (field degree).** For an integer m >= 1 put

    k(m) = the least prime strictly greater than 2 * m^(3/5).

This is well defined because there are infinitely many primes. By Bertrand's
postulate, k(m) < 4 * m^(3/5); hence the scan for k(m) (test 2,3,5,... by trial
division up to the square root) runs in time polynomial in m.

**Step 2 (payload length).** Put N(m) = m + k(m). Both terms are nondecreasing in
m and m is strictly increasing, so N(m) is strictly increasing. For the given N let

    m(N) = the largest m >= 1 with N(m) <= N.

If no such m exists (that is, N < N(1) = 1 + k(1)), set f_N identically 0.

**Step 3 (the input split).** Write m = m(N), k = k(m). Read the input
x in {0,1}^N as the concatenation of

* a *prefix* q = (q_0,...,q_{k-1}) in {0,1}^k, the first k bits;
* a *payload* y in {0,1}^m, the next m bits;
* a *tail* in {0,1}^{N-k-m}, the remaining bits.

The tail is ignored by f_N: replacing those bits never changes the value. This is
what makes f_N total at every integer length; the tail plays no role in the lower
bound, where it is fixed (Section 3).

**Step 4 (the field).** Let Q(Z) = Z^k + q_{k-1} Z^{k-1} + ... + q_1 Z + q_0 in
F_2[Z]. If Q is not irreducible over F_2, set f_N(x) = 0. Otherwise
F = F_2[Z]/(Q) is a field with 2^k elements, and we use the polynomial basis
1, Z, ..., Z^{k-1}, in which the *constant coefficient* of an element is its first
basis coordinate.

**Step 5 (blocks).** With r = ceil(m/k), the payload y is split into r consecutive
blocks of k bits each, the last block padded on the right by zeros. The i-th block
defines an element phi_i(y) in F, with its first bit as the constant coefficient.
The map y -> (phi_1(y),...,phi_r(y)) is F_2-linear and injective.

**Step 6 (the value).** Put

    y_i = phi_i(y)^(2^(i+1) - 1)  in F,   for i = 1..r,

and let Sym_r^t be the t-th elementary symmetric polynomial with t = ceil(sqrt(r)),

    Sym_r^t(y_1,...,y_r) = sum over t-element subsets I of {1..r} of prod_{i in I} y_i.

Finally

    f_N(x) = the constant coefficient of Sym_r^t(y_1,...,y_r)  in F,

that is, f_N(x) = pi(Sym_r^t(y_1,...,y_r)) where pi : F -> F_2 is the nontrivial
F_2-linear map "take the constant coefficient".

For a payload of length m with field degree k = k(m), r = ceil(m/k) and
t = ceil(sqrt(r)), this is exactly the function of equation (4) of Ben-Sasson and
Kopparty (2012) with p = 2, since 1 + 2 + 4 + ... + 2^i = 2^(i+1) - 1.

**Remark 1.1.** The prefix construction is deliberate. The published family fixes
a deterministic irreducible polynomial of degree k by an external algorithm. Here
the candidate polynomial is carried by the input and tested, which
* (i) needs only an irreducibility *test*, proved below in Section 2.2 from
   scratch, instead of a deterministic *construction* of an irreducible
   polynomial;
* (ii) keeps the family uniformly polynomial-time computable at every length;
* (iii) costs only the k = O(m^(3/5)) prefix bits, which are charged to the actual
   input length by the restriction argument of Section 3.

## 2. The family is total and computable in polynomial time

### 2.1 The schedule

m(N) is computable in polynomial time by the following deterministic procedure on
input N: for m = 1,2,... compute k(m) by trial division over the integers above
2 m^(3/5) (each trial division costs O(m^(3/5)) bit operations at most, and at most
4 m^(3/5) candidates are tried), test whether m + k(m) <= N, and stop at the first m
with m + k(m) > N, returning m - 1. Because N(m) is strictly increasing, the loop
stops at m = m(N) + 1, so it performs m(N) + 1 = O(N) iterations, each polynomial
in N. Hence m(N), k(m(N)) and r, t are computable in polynomial time.

### 2.2 Irreducibility testing

**Lemma 2.1.** Let k >= 2 be prime and let Q in F_2[Z] be monic of degree k. Then Q
is irreducible over F_2 if and only if

  (i) Z^(2^k) = Z  in  F_2[Z]/(Q),  and
  (ii) gcd(Z^2 + Z, Q) = 1.

*Proof.* If Q is irreducible then F_2[Z]/(Q) is the field F_{2^k}, in which every
element a satisfies a^(2^k) = a, giving (i); and Q has degree k >= 2 so it has no
root in F_2, giving (ii).

Conversely assume (i) and (ii). Let Q_1 be an irreducible factor of Q of degree e,
and let R = F_2[Z]/(Q_1) = F_{2^e}. Reducing (i) modulo Q_1 gives a^(2^k) = a for
a = Z mod Q_1. The fixed points of the map b -> b^(2^k) on F_{2^e} are exactly the
subfield F_{2^g} with g = gcd(k, e): indeed b^(2^k) = b iff b lies in the splitting
field of Z^(2^k) - Z inside F_{2^e}, which is the unique subfield of order
2^gcd(k,e). Hence a lies in F_{2^g}, so its minimal polynomial has degree dividing
g, but that degree is e, so e = g and therefore e divides k. As k is prime, e = 1
or e = k. If e = 1 then Q_1 is Z or Z + 1, contradicting (ii). Hence every
irreducible factor of Q has degree k, and since deg Q = k, Q is irreducible. []

Both conditions are testable in time polynomial in k: compute Z^(2^k) mod Q by k
successive squarings of polynomials of degree < k (each squaring and reduction is
O(k^2) bit operations), compare with Z, and run the Euclidean algorithm on
Z^2 + Z and Q. This is Rabin's test specialised to prime degree.

### 2.3 Field arithmetic, powers and the symmetric polynomial

For irreducible Q the ring F = F_2[Z]/(Q) is a field with 2^k elements, and
addition, multiplication and reduction modulo Q of k-bit words cost O(k^2) bit
operations. A power a^e with e < 2^(r+1) is computed by binary exponentiation in
O(r) multiplications, that is O(r k^2) bit operations.

The elementary symmetric polynomial is computed in characteristic two without any
division, by the recurrence obtained by splitting the t-subsets of {1..i} according
to whether they contain i:

    E_0 = 1,  E_j = 0 for 1 <= j <= t;
    for i = 1..r:
        for j = min(i,t) down to 1:
            E_j = E_j + y_i * E_(j-1);

the descending order of j ensures that y_i is used exactly once per update.
Correctness is by induction on i, and the cost is O(r t) field operations, that is
O(r t k^2) bit operations.

**Proposition 2.2 (uniformity).** The map x -> f_N(x) is computable by a single
deterministic algorithm running in time polynomial in N, for every N. Hence the
language L = {x : f_N(x) = 1, N = |x|} is in P, and therefore in NP: the verifier
V(x,w) accepts exactly when w is empty and the algorithm outputs 1, with witness
bound W(n) = 0, and both implications x in L <=> (exists w, |w| <= 0, V(x,w) = 1)
hold trivially.

*Proof.* Section 2.1 computes the schedule and k, r, t in polynomial time; Lemma
2.1 tests irreducibility in polynomial time (and the answer 0 is returned
immediately if the test fails); Section 2.3 computes the powers in O(r^2 k^2) and
the symmetric polynomial in O(r t k^2) bit operations. With m = O(N),
k = O(N^(3/5)), r = O(N^(2/5)), t = O(N^(1/5)), every stage is polynomial in N. []

**Remark 2.3 (all lengths).** f_N is defined for every integer N >= 1, including
all small lengths, by the schedule of Section 1. The tail is ignored, so each
length has a genuine total Boolean function; nothing is left undefined, and no
per-length advice is used, since the same algorithm computes every f_N.

## 3. The lower bound

### 3.1 A slice on which the family is the published one

Fix N >= 2^10, and let m = m(N), k = k(m) as in Section 1. Choose once and for all
one irreducible polynomial Q_0 in F_2[Z] of degree k; such a polynomial exists
because the number of monic irreducible polynomials of degree k over F_2 is
(2^k - 2)/k >= 1 (each element of F_{2^k} outside F_2 has a minimal polynomial of
degree exactly k, since k is prime, and each such polynomial has k roots).

Let S be the set of inputs x in {0,1}^N whose prefix equals the coefficient vector
of Q_0 and whose tail is identically zero. Then S is an affine subspace of F_2^N
of dimension m: it is the image of the injective affine map that sends the payload
y in F_2^m to (Q_0-coefficients, y, 0...0). On S,

    f_N(x) = pi( Sym_r^t( phi_1(y)^(2^2-1), ..., phi_r(y)^(2^(r+1)-1) ) ),

which is the function of `disperser.md`, equation (BK-4), for the field F = F_2[Z]/(Q_0)
of order 2^k, with r = ceil(m/k), t = ceil(sqrt(r)) and the injective F_2-linear
block map phi = (phi_1,...,phi_r).

### 3.2 The disperser property of the slice

**Quoted theorem (Ben-Sasson-Kopparty 2012, Theorem 2.2 with p = 2).** Let m be a
positive integer, let k be the least prime larger than 2 m^(3/5), let
r = ceil(m/k), t = ceil(sqrt(r)), let phi : F_2^m -> F_{2^k}^r be an injective
F_2-linear map and let pi : F_{2^k} -> F_2 be a nontrivial F_2-linear map. Then

    g(y) = pi( Sym_r^t( phi_1(y)^(1+2), phi_2(y)^(1+2+4), ..., phi_r(y)^(1+2+...+2^r) ) )

is an affine disperser for dimension greater than 6 m^(4/5): for every affine
subspace A of F_2^m with dim A > 6 m^(4/5) one has |g(A)| > 1.

*Hypotheses check.* k is prime and k is the least prime larger than 2 m^(3/5)
(Section 1, Step 1, with the payload length m); r = ceil(m/k), t = ceil(sqrt(r));
the block map of Step 5 is F_2-linear and injective, with zero padding of the last
block only; pi is the constant-coefficient functional, which is F_2-linear and not
identically zero. All hypotheses are met exactly, and the exponents agree because
1 + 2 + ... + 2^i = 2^(i+1) - 1.

BK state Theorem 2.2 for every positive integer m, and that is the form quoted above.
`disperser.md` reproduces the construction, the reduction of Theorem 2.2 to the
non-parameterised Theorem 7.3 of the same work (that reduction is carried out there
for every m >= 512, which is all this submission uses, since m = m(N) >= N/2 >= 512
whenever N >= 2^10), and the statements of the three internal results on which that
reduction rests; it also states precisely which parts are quoted rather than
re-proved.

### 3.3 Gate elimination applied to the slice

Let C be any B2 circuit computing f_N on {0,1}^N, of size |C|. Substituting the
fixed prefix bits and the fixed zero tail into C simplifies gates (constant or
unary gates are deleted, with wires and negations free) and yields a B2 circuit C'
computing the m-bit function g of Section 3.1, with |C'| <= |C|. By the quoted
theorem, g is an affine disperser for dimension greater than 6 m^(4/5), hence for
dimension

    d = floor( 6 m^(4/5) ) + 1,

since any subspace of dimension d has dimension greater than 6 m^(4/5). Theorem A
of `framework.md` applied to g and C' gives

    |C| >= |C'| >= (7/2) m - 7 d >= (7/2) m - 42 m^(4/5) - 7.

This holds for every circuit C computing f_N, so

    C_B2(f_N) >= (7/2) m(N) - 42 m(N)^(4/5) - 7.                              (3.1)

### 3.4 The remainder, with explicit constants

**Lemma 3.1.** For every integer N >= 2^10, with m = m(N) as in Section 1,

    m >= N - 7 N^(3/5) - 1     and     m >= N/2.

*Proof.* By Bertrand's postulate k(m) < 4 m^(3/5), so N(m) = m + k(m) <= m + 4 m^(3/5).
By maximality of m, N < N(m+1) = (m+1) + k(m+1) <= (m+1) + 4 (m+1)^(3/5) <= m + 1 + 7 m^(3/5)
for m >= 1 (using (m+1)^(3/5) <= 2^(3/5) m^(3/5) <= (7/4) m^(3/5)). Hence
N - m <= 7 m^(3/5) + 1 <= 7 N^(3/5) + 1, which is the first claim; and for N >= 2^10
one has 7 N^(3/5) + 1 <= N/2, which is the second. []

Substituting Lemma 3.1 into (3.1), and using m <= N and m^(4/5) <= N^(4/5),

    C_B2(f_N) >= (7/2)(N - 7 N^(3/5) - 1) - 42 N^(4/5) - 7
              = (7/2) N - (49/2) N^(3/5) - 42 N^(4/5) - 21/2
              >= (7/2) N - (49/2 + 42) N^(4/5) - 21/2
              >= (7/2) N - 67 N^(4/5)                                          (3.2)

for every N >= 2^10. Here the third line uses N^(3/5) <= N^(4/5) for N >= 1, the
first coefficient is 24.5 + 42 = 66.5, and the additive constant is absorbed by the
slack 0.5 N^(4/5) >= 10.5, which holds for every N with N^(4/5) >= 21, i.e. every
N >= 2^10. The coefficient 42 in (3.1) is 7*6, where 6 is the constant in the quoted
dimension bound; any other explicit constant C in d = C m^(4/5) would only change
the remainder, not the liminf.

### 3.5 The liminf

**Proposition 3.2.** liminf as N -> infinity of C_B2(f_N)/N >= 7/2. More precisely,
for every epsilon > 0 and every integer

    N >= N_epsilon := max( 2^10, ceil( (67/epsilon)^5 ) )

one has C_B2(f_N) >= (7/2 - epsilon) N.

*Proof.* For N >= (67/epsilon)^5 we have 67 N^(4/5) <= epsilon N, so (3.2) gives
C_B2(f_N) >= (7/2 - epsilon) N. []

Hence the full guarantee and the declared finite projection both hold for ONE fixed
eligible family, uniformly in the actual input length N, with all quantifiers
explicit: the family, its polynomial-time evaluator and its NP verifier do not
depend on epsilon, and coverage is over every integer length N >= N_epsilon.

## 4. Conformance with the model and the rules

* **Model.** Deterministic acyclic circuits, fan-in two, one output, all sixteen
  binary operations available, unrestricted depth and fan-out, inputs x only (the
  NP witness is never a circuit input), every retained gate counted once, inputs,
  wires, output selection, constants and negations free, exact agreement on every
  input, non-uniform circuits allowed. Negations are absorbed into the reading
  gate's operation, which is the normalisation rules.md requires; a retained
  negated projection would cost one gate, and the argument never relies on such a
  gate being free, because every deletion in `framework.md` is justified by the
  computed function becoming constant or a literal of an existing value.
* **Eligibility.** One fixed language L in P, one fixed verifier V and witness
  bound 0, both directions established in Proposition 2.2. Per-length advice,
  truth-table selection and non-constructive existence are not used; the field is
  built from bits of the input and an irreducibility test that is proved in Lemma
  2.1.
* **Lengths.** Every integer N >= 1 has its own total function (Section 1 and
  Remark 2.3); the bound covers every integer N >= N_epsilon, so no sparse
  subsequence, no limsup and no selected lengths are involved. Blocks, field
  encodings, zero padding of the last block and the k prefix bits are all charged
  against the actual input length: the prefix and tail are *fixed* in the
  restriction of Section 3.1, and the O(N^(3/5)) loss they cause is retained
  explicitly in (3.2).
* **Same family for both coordinates.** The coefficient 7/2 and the exponent 1 come
  from the single theorem (3.2) about the single family f_N; the projection
  (7/2, 1) is the exact leading behaviour of that one expression, and
  `expressions.py` certifies that the finite projection is implied by the declared
  polynomial_log guarantee (c = 7/2, p = 1, q = 0, r = 0).
* **Reduced rationals.** 7/2 is in lowest terms with positive denominator; the
  exponent is 1/1.
* **Evidence.** `framework.md` contains the complete class theorem (Theorem A),
  `disperser.md` the family's disperser construction and the quoted inputs, this
  file the family, eligibility, assembly, remainder and quantifiers, and
  `check_cases.py` an exact-rational audit of every case inequality, of the final
  summation, and of the constants in Lemma 3.1 and (3.2).

## 5. What is used, and what is not claimed

Used, and nothing else:
1. The gate elimination theorem of `framework.md`, proved in full here; it is a
   reproduction of Mishra-Santhanam (ECCC TR26-260), Theorem 4.1.
2. The explicit affine disperser of Ben-Sasson-Kopparty (2012), Theorem 2.2, in the
   form quoted in Section 3.2 and reproduced (construction and reduction, with the
   internal lemmas quoted by number) in `disperser.md`.
3. Bertrand's postulate, for the size of the least prime above 2 m^(3/5).

Not claimed:
* No improvement over the published 3.5n - O(log n) bound is claimed, and no claim
  is made about robustness of the constant 7/2 under a different potential: within
  the three-parameter family of potentials described in `framework.md` Section 10
  the constant 7/2 is optimal, and a stronger coefficient would need a different
  elimination step.
* Nothing is claimed about exponents other than 1, about lane or global endpoints,
  or about any function other than the family defined in Section 1.
* `check_cases.py` verifies arithmetic, not mathematics; a passing run does not
  discharge any proof obligation.
