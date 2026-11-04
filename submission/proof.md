# One total uniform family and its coefficient-four assembly

This file gives the complete family, eligibility, parameter and all-length
assembly. Its new circuit conclusion uses exactly the class lemma stated in
Section 5 and proved separately in `elimination.md`. A missing or unproved
class lemma is a gap in the main theorem; neither a prior review nor the
public coefficient-four idea note establishes that lemma.

## 1. Provenance and the precise theorem

The certificate-prefix family and the reproduced BK disperser argument are
reused from mooselumph's public submission
`f4f2c464-3045-4c32-a6e4-32478c05f343`, candidate commit
`8db0d3171ce7d19a48503ae34b5945331da55535`. They also appear in the
repository's `tests/fixtures/candidate-3_5` fixture. This reuse concerns
mathematical text and construction, not any previous validity judgment.
The public idea note of submission
`5e9fcf7f-8338-4729-9faa-a18994661af6` suggested the coefficient-four
elimination target. It supplies an idea rather than an
accepted theorem premise; the necessary argument is reproduced and must
stand on its own in `elimination.md`.

For the single language below, with N the entire actual input bit length,
the resulting theorem is

    C_B2(F_N) >= 4N - O(N^(4/5))

for every sufficiently large integer N. Thus
liminf_(N->infinity) C_B2(F_N)/N >= 4. The coefficient and exponent refer
to this same family. The whole family is not asserted to be an affine
disperser: the property is established after a legitimate prefix restriction.

## 2. One exact evaluator at every length

Every integer choice below uses ordinary nonnegative integer arithmetic.
On a binary input x of length N, the deterministic evaluator A does this:

1. If N=0, return 0.
2. Find the least positive integer B with B^5>32*N^3. Find the least prime
   integer k>=B by testing successive integers using deterministic trial
   division. Here primality means k>=2 with no divisor from 2 through
   floor(sqrt(k)).
3. If k>=N, return 0. Otherwise set m=N-k, r=ceil(m/k), and let t be the
   least nonnegative integer with t^2>=r. If t>=r or k*t<r*(r+1), return 0.
4. Read the first k input bits b_0,...,b_(k-1) and form the monic polynomial

       Q(Z)=Z^k + sum_(j=0)^(k-1) b_j Z^j  over F_2.

   Check both Z^(2^k)=Z modulo Q and gcd(Q,Z^2+Z)=1. If either fails,
   return 0. Compute the first test by k repeated squarings modulo Q,
   starting at Z modulo Q; no integer of bit length 2^k is constructed.
5. Work in K=F_2[Z]/(Q). Append exactly rk-m zero bits to the payload,
   divide into r consecutive blocks of k bits, and regard block i as
   a_i=sum_(j=0)^(k-1) a_(i,j) Z^j. Thus the first bit in every block
   is its constant coefficient. The prefix is not part of any field block.
6. Put w_i=r+2-i and y_i=a_i^(2^(w_i)-1), for i=1,...,r. Return the
   constant coefficient of

       S_t(y_1,...,y_r)=sum_(I subset {1,...,r}, |I|=t) product_(i in I) y_i.

The exponents descend strictly with i. They are exactly those in
`disperser.md`. The output functional pi is the constant coefficient of
the unique representative of degree below k; it is F_2-linear and pi(1)=1.
In particular it is nonzero. We do not replace this functional with an
inconsistent trace or reverse the exponent order.

Step 6 has a polynomial-time implementation. Initialize E_0=1 and all
E_j=0 for 1<=j<=t. For i=1,...,r, first compute y_i by binary powering;
then update E_j <- E_j+y_i*E_(j-1) for j=min(i,t),...,1 in descending
order. Partitioning subsets according to whether they contain i proves
that after processing i, E_j is the order-j elementary symmetric
polynomial in the first i values. Hence E_t is exactly the specified sum.
No division is used. Each exponent 2^(w_i)-1 has w_i one-bits, so binary
powering uses O(w_i) quotient-ring multiplications.

All failed guards return the specified constant, and the valid branch is
also specified on every input. There is no hidden choice of a modulus,
truth table, advice string, witness or hard prefix in this evaluator.

## 3. Exact irreducibility test, existence, and uniform eligibility

We justify the special prime-degree test, including its converse. In an
algebraic closure of F_2 the polynomial H(Z)=Z^(2^k)-Z has derivative 1,
so it has exactly 2^k distinct roots. Its roots form a field: Frobenius
shows closure under addition and multiplication, and if a is a nonzero
root, (a^(-1))^(2^k)=a^(-1). Call this field K_0. The automorphism
sigma(a)=a^2 has sigma^k=id on K_0. An orbit has length dividing k, since
its period divides any exponent returning that element to itself.

For an orbit O, the monic polynomial product_(a in O)(Z-a) is fixed by
squaring all its coefficients. These coefficients are in K_0 and satisfy
c^2=c, so belong to F_2. This polynomial is irreducible over F_2: the roots
of any nonconstant F_2 factor would be a nonempty subset of O invariant
under sigma, and an orbit has no proper such subset. Conversely any monic
irreducible factor of H must consist of one orbit, by this construction
and irreducibility. Therefore every irreducible factor of H has degree
1 or k, because k is prime. Its degree-one roots are exactly 0 and 1.
There are 2^k>2 roots, so at least one orbit has length k. Thus an
irreducible monic degree-k Q exists for every prime k used by A.

The first test says Q divides the squarefree polynomial H. Consequently Q
is a product of distinct irreducible factors of degrees 1 or k. The gcd
test excludes exactly the two degree-one factors Z and Z+1. Since deg Q=k,
a passing Q is one irreducible degree-k factor. Conversely, every monic
irreducible degree-k Q divides H and has neither degree-one factor, so
passes. This proves the test in both directions. In particular every
accepted quotient is a field, and a passing prefix exists at every length
whose length guards pass. The existence proof makes no algorithmic search
for that prefix part of A.

All length comparisons are exact. Binary search finds B by fifth-power
comparisons with 32*N^3. Integer division computes r, and binary search
with square comparisons computes t. For N>=1, B>=3. Bertrand's theorem
in the form "for every integer B>1 there is a prime p with B<p<2B"
therefore applies, giving

    2*N^(3/5) < B <= k < 2B <= 4*N^(3/5)+2.                (1)

(The last inequality can be made strict, but is not needed.) Thus the
prime scan and its trial divisions take polynomial time in N. The
standard theorem's hypotheses are checked here; no prime-gap estimate
stronger than Bertrand is needed.

In the branch k<N, the degree-k quotient representatives use k bits.
Schoolbook multiplication followed by reduction modulo the monic Q takes
O(k^2) bit operations. Polynomial Euclid and k modular squarings take
polynomial time even on rejected prefixes. The powers in Step 6 use
O(r^2) multiplications in total, and the symmetric recurrence O(rt).
Here k,r,t<=N and all scalar integer operands have polynomial bit length.
Together with the prime search, a coarse fixed O((N+1)^6) bit bound
suffices for A on every string. The early guard covers the other branch.

Let L={x:A(x)=1}, and F_N be its characteristic function on {0,1}^N.
Take the one fixed witness bound W(N)=0 and verifier

    V(x,w)=1 iff w is empty and A(x)=1.

If x is in L, the empty witness is permitted and accepted. If a permitted
witness is accepted then A(x)=1, so x is in L. Both directions hold with
one polynomial-time verifier. In fact L is in P. Circuit inputs are only
the N bits of x; witness bits are never circuit inputs.

## 4. Every sufficiently large integer length is hard

Inequality (1) holds for every positive integer N. It gives
k=Theta(N^(3/5)) and k/N->0. Hence k<N eventually, and m=N-k=Theta(N).
Then the exact ceiling definitions give

    r=Theta(N^(2/5)),     t=Theta(N^(1/5)),
    r < (1/2)N^(2/5)+1.

In particular r->infinity, so t<r eventually. To check the other guard,
use t>=sqrt(r) and the displayed upper bound for r:

    r(r+1)/t <= sqrt(r)(r+1)
      <= sqrt((1/2)N^(2/5)+1)*((1/2)N^(2/5)+2).

The last quantity divided by N^(3/5) tends to 1/(2*sqrt(2)), whereas
k/N^(3/5)>2. It follows that k*t>=r(r+1) for every sufficiently large
integer N. These are uniform inequalities on all integers, not estimates
along a prime or sparse subsequence. Finite guarded lengths require no
special circuit claim.

For each such N fix any accepted irreducible prefix Q; Section 3 proves
that one exists. With these k original input bits fixed, the remaining
function G_Q has exactly m=N-k Boolean inputs. Put

    d_N = floor(kr/t) + kt + r(r+1) + 1.                  (2)

The proved lemma in `disperser.md` applies because k is prime, 1<=t<r,
and kt>=r(r+1). It says the field function in Step 6 is nonconstant on
every affine subspace of K^r with dimension greater than
kr/t+kt+r(r+1), for the exact nonzero functional pi already specified.

Appending the fixed zero bits defines an injective F_2-linear map from
F_2^m into K^r. The image of each affine input subspace has exactly its
original dimension. By (2), G_Q is therefore nonconstant on every affine
subspace of dimension d_N (and on all larger ones). This is the dimension-d_N
affine-disperser property used below. The rk-m<k padding bits are fixed
zeros and do not contribute to m or N. They lose no affine dimension;
we apply the gate lemma at m, never at the artificially enlarged rk.

Finally, kr<=m+k=N gives kr=O(N). From the uniform estimates for k,r,t,

    floor(kr/t)=O(N^(4/5)),   kt=O(N^(4/5)),
    r(r+1)=O(N^(4/5)),       d_N=O(N^(4/5)).             (3)

In particular d_N/m->0, so 1<=d_N<=m eventually. No vacuous
high-dimensional disperser condition is being used to obtain a positive
leading coefficient.

## 5. Full-B2 model, gate lemma, and the actual-length bound

The sole class lower bound needed from `elimination.md` is:

    For 1<=d<=m, if g:{0,1}^m->{0,1} is nonconstant on every
    dimension-d affine subspace, then every exact deterministic
    full-B2 circuit C for g has |C| >= 4m-9(d+1).           (4)

Here full B2 means the pinned general model with all sixteen binary truth
tables, acyclic circuits, one output, arbitrary sharing and fan-out, two
pins per gate including repeated pins, free input/wire/output selection,
and the gate counting and constant-substitution normalization in rules.md.
Unary, constant and negated-projection retained gates are covered by (4),
as are terminal outputs. The complete proof of (4), rather than any
citation to an unproved four-coefficient theorem, is required evidence in
`elimination.md`. No formula, monotone or direct-sum transfer is assumed.

Take any B2 circuit for F_N. Fix its k prefix inputs to the coefficients
of the chosen Q. It now computes G_Q exactly on all m payload inputs.
Substitution adds no gates: it only attaches fixed constants to former
input pins, and the no-larger normalized circuit is within the model of
(4). Constants introduced here are dealt with by (4)'s normalization;
no free computational basis is being appended. Therefore (2)--(4) give

    C_B2(F_N) >= 4m-9(d_N+1)
              = 4N - [4k+9(d_N+1)].                     (5)

All k certificate bits are thus charged at their actual input length.
Restriction does not need an efficient rule selecting Q: an exact circuit
for F_N must work on every prefix. The evaluator itself reads and tests
its input prefix and was fixed before any asymptotic quantifier.

By (1) and (3), R(N)=4k+9(d_N+1)=O(N^(4/5))=o(N) eventually.
Extend R arbitrarily over the finitely many failed-guard lengths. For any
real epsilon>0 choose N_epsilon large enough that all the preceding
conditions hold and R(N)/N<=epsilon for every integer N>=N_epsilon.
Then (5) yields

    C_B2(F_N) >= (4-epsilon)N for every integer N>=N_epsilon.

This proves precisely the claimed liminf and the coefficient-four,
exponent-one polynomial_log guarantee with both logarithmic exponents zero.
In the equivalent relative-epsilon version rendered by `expressions.py`,
replace epsilon here by 4*epsilon. There is one fixed family throughout.
No finite experiment, explicit finite threshold or stronger exponent is
asserted.

## 6. Source scope and supplied evidence

`disperser.md` reproduces the needed algebraic argument specialized to
binary fields of prime degree, including coefficient isolation and the
triangular mixing argument. Its source attribution is Ben-Sasson and
Kopparty, *Affine Dispersers from Subspace Polynomials*, author preprint
2010-03-12, Lemmas 3.7, 3.10, 3.13, Theorem 6.1, Lemma 7.2, Theorem 7.3
and Claim 7.4:
https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
(SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf,
as recorded in the reused package). The reproduced argument, not retrieval
of that link, supplies its mathematical premises.

The predecessor's parameter-method background is Mishra--Santhanam,
ECCC TR26-260, 2026-10-19, Theorem 4.1 and Sections 2--4. Its earlier
7/2 conclusion does not imply (4). The proof of the improved gate lemma
is separately supplied in `elimination.md`. Li--Yang, ECCC TR21-023,
2021-02-20, Sections 1.1, 2.1, 2.4 (footnotes 3--4) and 3.3 identify the
model pinned by rules.md; no Li--Yang LP bound is used here.

The package is proof-only: `experiments` and `checks` are empty. There is
no declared computation whose output is a mathematical premise. The two
main claims are identical canonical renderings of `bound.json`, and their
identical sorted evidence lists contain all supplied proof and bound files.
