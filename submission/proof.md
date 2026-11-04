# A uniform certificate-prefix family: coefficient 4, exponent 1

## Claim and provenance

For the single total Boolean language defined below, with **N denoting the
entire input length**, we prove

    C_B2(F_N) >= 4N - R(N),       R(N)/N -> 0.

The joint pair is (4,1). The declared guarantee is the explicit eventual
bound proved in section 5a below, stronger than the leading liminf alone. This applies existing mathematics; it does not
claim discovery of a new circuit lower bound. The gate argument is the ternary-pair refinement supplied in production
submission 237016f9-ece3-4a05-9933-dfcf6056e787, reproduced in elimination.md.
Its fair affine systems and role exchanges derive from Mishra--Santhanam,
ECCC TR26-260, October 19, 2026, Sections 2-4. The algebraic construction is based on Ben-Sasson
and Kopparty (BK), *Affine Dispersers from Subspace Polynomials*, author
version March 12, 2010, Sections 3, 6 and 7. The necessary arguments are
supplied in our notation in `elimination.md` and `disperser.md`.

Our input-prefix wrapper removes the deterministic irreducible-polynomial
construction prerequisite of the repository's BK/Shoup seed. It explicitly
pays for all prefix bits. The **whole** function is not an affine disperser:
only a fixed-prefix restriction is asserted to be one. No conjecture,
restricted circuit basis, direct-sum assumption or exponent gain is used.
We assert no exponent above 1, endpoint, or optimality over all families.

## 1. One evaluator on every bit string

The first bit of each field block is its constant coefficient. On x of length N:

1. If N=0, output 0.
2. Let B be the smallest positive integer with B^5 > 4*N^3. Let k be the
   smallest prime at least B. Use deterministic trial division to find k.
3. If k>=N, output 0. Otherwise put m=N-k, r=ceil(m/k), t=ceil(sqrt(r)).
   If t>=r or k*t<r*(r+1), output 0.
4. Read the first k bits b_0,...,b_(k-1) as

       Q(Z)=Z^k + sum_(j=0)^(k-1) b_j Z^j    in F_2[Z].

   Check Z^(2^k)=Z modulo Q and gcd(Q,Z^2+Z)=1. If either fails, output 0.
5. Use K=F_2[Z]/(Q). Append rk-m zero bits to the remaining m payload bits
   and interpret the r blocks of k bits as a_1,...,a_r in K. This is an
   injective F_2-linear map of the payload.
6. For i=1,...,r set e_i=r+1-i and y_i=a_i^(2^(e_i+1)-1). Output the
   constant coefficient of the elementary symmetric polynomial of order t
   in y_1,...,y_r, in its representative of degree less than k.

For explicit evaluation, initialize E_0=1 and E_1=...=E_t=0. For each i,
update E_j <- E_j+y_i*E_(j-1), in decreasing order j=min(i,t),...,1.
Splitting subsets by whether they contain i proves this division-free
recurrence. Binary powering takes O(e_i+1) multiplications per y_i.

This specifies every choice and every input length. Q is part of x, not
advice or a value supplied by an unspecified construction algorithm.

## 2. Exact prefix test and existence of hard prefixes

In an algebraic closure of F_2, Z^(2^k)-Z has 2^k distinct roots because
its derivative is 1. The roots are closed under addition, multiplication
and nonzero inversion, hence form a field of size 2^k. Frobenius z->z^2
permutes them and its k-th iterate is the identity. A Frobenius orbit of
length d gives an irreducible degree-d polynomial over F_2: its product
of linear factors is Frobenius-fixed, while a proper factor over F_2 would
require a proper invariant subset of the orbit. Conversely, an irreducible
factor's roots form one orbit. Thus factor degrees divide k.

Since k is prime, the orbits have size 1 or k. Only 0 and 1 are fixed,
and 2^k>2, so an irreducible degree-k polynomial exists for every such k.
No algorithm selecting it is needed in the uniform evaluator.

The first prefix test says Q divides the squarefree polynomial Z^(2^k)-Z;
the second excludes the factors Z and Z+1. Remaining irreducible factors
all have degree k. Since Q has degree k, it must itself be irreducible.
Conversely, every irreducible degree-k Q passes both tests. The test is
therefore exact. Compute its first part by k squarings modulo Q and its
second by polynomial Euclid, both in polynomial bit time.

## 3. Uniform P and NP with empty witnesses

Bertrand's prime theorem gives k<2B, and so

    beta*N^(3/5) < k < 2*beta*N^(3/5)+2,  beta=4^(1/5).

The prime scan and trial division are polynomial in N. B is computed by
integer binary search and fifth-power comparisons. All subsequent integer,
polynomial and field operations have polynomial bit cost. Field multiplication
can use schoolbook O(k^2) operations; the powers take O(r^2) multiplications
in total and the symmetric recurrence O(rt). A coarse O((N+1)^6) bound
suffices for the whole algorithm, including rejected inputs.

Let A be this fixed evaluator, L={x:A(x)=1}, W(N)=0, and let V(x,w)
accept exactly when w is empty and A(x)=1. If x is in L the empty witness
is accepted; if a permitted witness is accepted then A(x)=1 and x is in L.
This proves both directions with one uniform polynomial-time verifier.
Witnesses are not inputs to the target circuits.

## 4. All sufficiently large integer lengths pass the guards

We have k=Theta(N^(3/5)), m=N-k, r=Theta(N^(2/5)), t=Theta(N^(1/5)).
Consequently k<N and t<r eventually. Also r<N/k<beta^(-1)*N^(2/5) and

    r(r+1)/t <= sqrt(r)(r+1)
              <= sqrt(beta^(-1)*N^(2/5))*(beta^(-1)*N^(2/5)+1).

Dividing the last expression by N^(3/5) gives limit beta^(-3/2)<1,
whereas k/N^(3/5)>beta>1. Thus k*t>=r(r+1) eventually. These estimates hold
for every sufficiently large integer N, not just selected lengths.

## 5. Hard restriction and its precise loss

For each sufficiently large N fix any irreducible Q of degree k and fix
the first k bits of x to its coefficients. The remaining function G_Q
has m inputs. Put

    d_N = floor(r(r+3)/2+max(k(t-1),kr/t))+1.

Section 9 of `disperser.md` proves that, for prime field degree k, r>t>=1 and
k*t>=r(r+1), the function

    (a_i) -> pi(sum_(|I|=t) product_(i in I) a_i^(2^(r+2-i)-1))

is nonconstant on every affine subspace of K^r of dimension greater than
r(r+3)/2+max(k(t-1),kr/t), for any nonzero F_2-linear pi. Our pi is the constant
coefficient, with pi(1)=1. The payload embedding is injective and linear,
so G_Q is an affine disperser of dimension d_N on its m actual input bits.
Zero padding preserves dimension and does not add circuit inputs.

Since kr=O(N), r=Theta(N^(2/5)) and t=Theta(N^(1/5)),

    d_N=O(N^(4/5)),          k=O(N^(3/5)).

We use this conservative sublinear disperser dimension, not later
logarithmic-dimension constructions. Either suffices for the leading pair.

Fixing prefix inputs of any circuit for F_N yields a circuit for G_Q with
no additional gates. `elimination.md` proves the full-B2 bound

    |C| >= 4m-9d

for a dimension-d affine disperser: elimination.md assumes nonconstancy
on dimensions strictly greater than an integer h and proves 4m-9(h+1).
Set h=d-1, since nonconstancy on every dimension-d affine subspace also
holds on every larger one. This gives exactly 4m-9d. For the eventual
application d>=1. The proof includes arbitrary sharing,
repeated pins, constants, degenerate gates and terminal-output cases.
Apply it to the restricted circuit (or its no-larger normalized form):

    C_B2(F_N) >= 4(N-k)-9d_N
               = 4N-(4k+9d_N).

Set R(N)=4k+9d_N eventually and extend it on the finitely many guarded
lengths. Then R(N)=O(N^(4/5))=o(N). The choice of Q in this lower-bound
restriction is not advice to the evaluator: every circuit must work for
all prefixes, while the uniform algorithm simply reads and tests its input.

For each epsilon>0, take N_epsilon large enough that the guards pass and
R(N)/N<=epsilon for every N>=N_epsilon. This gives exactly
C_B2(F_N)>=(4-epsilon)N for all those integers. Both coordinates in
bound.json belong to this one theorem and one language.

## 5a. A single explicit eventual loss using prime density

We use the standard prime number theorem in its precise classical form:
if pi(x) is the number of primes at most the positive real x, then
pi(x)/(x/ln(x)) tends to 1 as x tends to infinity. This is an
unconditional standard theorem, not a conjecture about prime gaps.
No effective error term or computable numerical threshold is assumed.

Its needed consequence is that the least prime p at least the integer B
satisfies p/B -> 1. For any fixed real eta>0, the theorem gives

    pi((1+eta)B)-pi(B) ~ eta*B/ln(B) > 0

for every sufficiently large B. Thus there is a prime between B and
(1+eta)B, and 1<=p/B<=1+eta eventually. Letting eta be arbitrary proves
the stated ratio limit over all integers B. This derivation also proves
the required all-length scope; no selected subsequence is used.

Our integer threshold B^5>4*N^3 makes B the least integer strictly
above beta*N^(3/5), where beta=4^(1/5)=2^(2/5)>1. Thus the prime-density
consequence and elementary rounding give

    k/N^(3/5)->beta,
    r/N^(2/5)->1/beta,
    t/N^(1/5)->1/sqrt(beta),
    kr/N->1.

Here r=ceil(N/k)-1 and t=ceil(sqrt(r)); additive rounding errors divided
by the diverging scales tend to zero. Each of k(t-1)/N^(4/5) and
(kr/t)/N^(4/5) tends to sqrt(beta)=2^(1/5). Their maximum has the same
limit. Also H_r/N^(4/5)->1/(2*beta^2)=2^(1/5)/4, and floor plus one
changes the scaled limit by zero. Therefore

    d_N/N^(4/5)->(5/4)*2^(1/5),
    (4k+9d_N)/N^(4/5)->(45/4)*2^(1/5).

The final limit is strictly less than the declared rational 12923/1000.
The exact verification is

    12923^5-2*11250^5=20011730845880843>0.

Since positive fifth powers are strictly increasing, this implies
2^(1/5)<12923/11250. Multiplying by 45/4 gives
(45/4)*2^(1/5)<12923/1000. By convergence with this strict fixed slack,
there is one threshold N0 such that at every integer N>=N0,

    4k+9d_N <= (12923/1000)*N^(4/5).

Increase N0 to include the all-length guard threshold from section 4.
The previously proved circuit restriction theorem now gives

    C_B2(F_N)>=4N-(12923/1000)*N^(4/5)

at every such integer length. This is a fixed explicit loss, not an
unknown Big-O constant or a liminf-only guarantee. No numerical threshold
or effective prime-gap estimate is claimed. Its finite projection [4,1]
follows by dividing the loss by N. Prime density affects only the analysis;
the evaluator still uses elementary trial division and Bertrand's bound
for polynomial running time.

## 5b. Provenance and refinement scope

The original seven-file certificate-prefix construction and algebraic
proof come from production f4f2c464-3045-4c32-a6e4-32478c05f343 and its
matching protected fixture. Positive-excess selection was developed in
ab6c806f-a0bd-41c9-b00a-bdd84324bb3f. The exact block-weight refinement
is supplied in disperser.md section 9. The ternary-pair coefficient-four
gate argument is reused verbatim from production
237016f9-ece3-4a05-9933-dfcf6056e787. No prior judgment is a premise.

This candidate uses the fixed degree threshold beta*N^(3/5), implemented
by B^5>4*N^3, and standard prime density to prove the single explicit loss
above. It retains the complete construction, all-length uniformity, input
accounting and circuit-model proof. The theorem concerns this one new
total family, not a family chosen in response to a coefficient or exponent
quantifier. The rational loss constant is verified by exact integer arithmetic.

## 6. Evidence and limits

`disperser.md` contains the coefficient separation and block-mixing proof.
`elimination.md` contains the complete ternary-pair gate argument. `audit.py` supplies
supplementary finite checks of Boolean identities, rational estimates and
the prefix test; these do not replace either universal proof. No claimed
review, checker success or prior acceptance is a mathematical premise.

Sources, for attribution and version identification:

* Mishra--Santhanam, ECCC TR26-260, 2026-10-19, Theorem 4.1 and Sections 2-4,
  https://eccc.weizmann.ac.il/report/2026/260/download/
  SHA256 38e9ff3cea4ac0399723756f00c2425036c438887ab5ae4abb10d09fa971a3bc.
* Ben-Sasson--Kopparty, author preprint 2010-03-12, Lemmas 3.7, 3.10, 3.13,
  Theorem 6.1, Lemma 7.2, Theorem 7.3 and Claim 7.4,
  https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
  SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf.
* Li--Yang, ECCC TR21-023, 2021-02-20, Sections 1.1, 2.1, 2.4 and 3.3,
  for the protected full-B2 conventions. Their LP and lower-bound case
  analysis are not invoked by this submission.
