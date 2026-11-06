# A uniform certificate-prefix family: coefficient 3, exponent 1

## Claim and provenance

For the single total Boolean language F defined below, with **N the entire
input length**, we prove

    C_B2(F_N) >= 3N - R(N),       R(N)/N -> 0.

The declared joint pair is (3,1). This applies existing mathematics; it
claims no new circuit lower bound. The algebraic construction is the
affine-disperser family of Ben-Sasson and Kopparty (BK), *Affine
Dispersers from Subspace Polynomials*, author version 2010-03-12,
Sections 3, 6, 7, re-proved in our notation in `disperser.md`. The gate
argument in `elimination.md` is a self-contained flat-potential parameter
argument in the tradition of classical gate elimination (Schnorr 1974;
Wegener 1987, Ch. 5); it is deliberately simpler and weaker than the full
parameter method of Mishra--Santhanam (ECCC TR26-260, 2026-10-19), whose
public application by this challenge's organizer reaches 7/2 and strictly
dominates this submission.

Our input-prefix wrapper removes any need to construct an irreducible
polynomial by an unspecified algorithm: the polynomial is part of the
input, tested explicitly, and all prefix bits are paid for. The WHOLE
function is not claimed to be an affine disperser; only fixed-prefix
restrictions are. No conjecture, restricted basis, direct-sum assumption
or exponent above 1 is used.

## 1. One evaluator on every bit string

The first bit of each field block is its constant coefficient. On input x
of length N:

1. If N=0, output 0.
2. Let B be the smallest positive integer with B^5 > 64*N^3. Let k be the
   smallest prime at least B, found by deterministic trial division.
3. If k>=N, output 0. Otherwise put m=N-k, r=ceil(m/k), t=ceil(sqrt(r)).
   If t>=r or k*t<r*(r+1), output 0.
4. Read the first k bits b_0,...,b_(k-1) as the monic polynomial

       Q(Z) = Z^k + sum_(j=0)^(k-1) b_j Z^j    in F_2[Z].

   Test Z^(2^k) = Z (mod Q) and gcd(Q, Z^2+Z)=1. If either fails,
   output 0.
5. Otherwise Q is irreducible of degree k (Section 2); put K=F_2[Z]/(Q).
   Append r*k-m zero bits to the remaining m payload bits and read the r
   blocks of k bits as a_1,...,a_r in K. This padding is an injective
   F_2-linear map of the payload into K^r.
6. For i=1,...,r put w_i=r+2-i and y_i=a_i^(2^(w_i)-1). Output the
   constant coefficient of the elementary symmetric polynomial sigma_t
   of y_1,...,y_r, as its canonical remainder of degree < k modulo Q.

For explicit evaluation, initialize E_0=1, E_1=...=E_t=0 and, for each i,
update E_j <- E_j + y_i*E_(j-1) in decreasing order j=min(i,t),...,1;
splitting the subsets counted by sigma_t according to whether they
contain i proves this recurrence without division. Each y_i is obtained
by binary powering with O(w_i) multiplications in K.

Every choice at every length is fixed by this rule; Q is a piece of the
input, not advice and not the output of an unspecified construction.

## 2. The prefix test is exact; hard prefixes exist

In an algebraic closure of F_2, Z^(2^k)-Z has 2^k distinct roots
(derivative 1). The root set is closed under addition, multiplication
and nonzero inversion, hence is a field of 2^k elements, on which
Frobenius z->z^2 acts with k-th iterate the identity. The roots of an
irreducible factor over F_2 form exactly one Frobenius orbit, and a
Frobenius orbit of length d yields an irreducible degree-d factor; so
the factor degrees of any divisor of Z^(2^k)-Z divide k.

Since k is prime, orbits have size 1 or k; only 0 and 1 are fixed, and
2^k>2, so at least one orbit of size k exists: an irreducible degree-k
polynomial exists for every prime k. No selection algorithm is needed by
the evaluator.

The first test says Q divides the squarefree Z^(2^k)-Z; the second
excludes the factors Z and Z+1. Every other irreducible factor has
degree k, and Q itself has degree k, so Q is irreducible. Conversely
every irreducible degree-k polynomial passes both tests. So the test is
exact. It runs in polynomial time: k squarings modulo Q for the first
part, the Euclidean algorithm for the second.

## 3. Uniform P and NP with empty witnesses

We use Bertrand's postulate (for every integer n>1 some prime p satisfies
n<p<2n; a standard textbook theorem, e.g. Hardy--Wright Thm. 418):
applied to B (when B is itself prime, k=B), it gives B<=k<2B, so with
B<=2^(6/5)*N^(3/5)+1,

    2^(6/5)*N^(3/5) < B <= k < 2^(11/5)*N^(3/5)+2 = O(N^(3/5)).

B is found by binary search on fifth powers; the prime scan and trial
division are polynomial; schoolbook arithmetic in K costs O(k^2) bit
operations per multiplication; each y_i costs O(w_i)<=O(r) such
multiplications, and the symmetric recurrence costs O(r*t). Rejected
inputs halt earlier. A coarse O((N+1)^6) bound covers the whole
algorithm.

Let A be this one fixed evaluator, L={x : A(x)=1}, W(N)=0, and
V(x,w)=1 iff w is empty and A(x)=1. If x is in L the empty witness is
accepted; conversely an accepted witness forces A(x)=1. Both directions
hold for the single uniform polynomial-time verifier A; witness bits are
never inputs to the target circuits.

## 4. The guards hold at every sufficiently large length

Write m=N-k. From Section 3, k=Theta(N^(3/5)), hence r=ceil(m/k) satisfies
r <= N/k < 2^(-6/5)*N^(2/5), and for large N, k<=N/2 so m>=N/2 and
r >= m/(2k) > 2^(-16/5)*N^(2/5); so r=Theta(N^(2/5)) and
t=ceil(sqrt(r))=Theta(N^(1/5)); in particular t<r eventually. For the
second guard, using t>=sqrt(r) and r+1<=2r,

    r(r+1)/(k*t) <= sqrt(r)*(r+1)/k <= 2*r^(3/2)/k
                  <= 2*2^(-9/5)*N^(3/5) / (2^(6/5)*N^(3/5))
                  = 2^(-2) = 1/4 < 1,

for every N large enough that the two upper bounds on r and sqrt(r)
apply; an explicit threshold is computed in `audit.py`. So k*t>r*(r+1)
for every sufficiently large integer N. These estimates hold for EVERY
large integer, not for a subsequence.

## 5. The payload restriction is an affine disperser

Fix any N passing the guards and any irreducible Q of degree k; hardwire
the first k input bits to the coefficients of Q. The remaining function
G_Q has the m payload bits as inputs, maps them linearly and injectively
into K^r by zero padding, and outputs

    G_Q(a_1,...,a_r) = pi( sum_(|I|=t) product_(i in I) a_i^(2^(w_i)-1) ),

pi = the constant-coordinate functional (nonzero, with pi(1)=1). This is
exactly the function of Lemma D of `disperser.md`. That lemma, with the
guards t<r and k*t>=r*(r+1), proves G_Q nonconstant on every affine
subspace of K^r of dimension strictly greater than kr/t+kt+r*(r+1).
Zero padding is linear and injective, so it carries affine subspaces of
F_2^m of a given dimension to affine subspaces of K^r of the SAME
dimension (an injective linear map between F_2-spaces preserves affine
dimensions). Therefore G_Q is an affine disperser of dimension

    d_N = floor(k*r/t) + k*t + r*(r+1) + 1

on its m actual input bits: nonconstant on every affine subspace of
dimension at least d_N.

With r=Theta(N^(2/5)) and t=Theta(N^(1/5)):

    kr/t = O(N*N^(-1/5)) = O(N^(4/5)),   kt = O(N^(4/5)),
    r(r+1) = O(N^(4/5)),   so d_N = O(N^(4/5)),   k = O(N^(3/5)).

## 6. Assembly of the lower bound

Let C be any full-B2 circuit computing F_N. Hardwiring the first k input
bits of C to the coefficients of Q replaces those inputs by constants —
a free operation in the full-B2 model that changes no gate — and yields
a circuit of size at most |C| computing G_Q on the m payload bits. By
Lemma E of `elimination.md`, since G_Q is an affine disperser of
dimension d_N,

    |C| >= 3m - 6*d_N = 3(N-k) - 6*d_N = 3N - R(N),

    R(N) = 3k + 6*d_N = O(N^(4/5)).

Extend R arbitrarily on the finitely many guarded lengths. Then
R(N)/N -> 0. Given epsilon>0, choose N_epsilon larger than the guard
threshold of Section 4 and with R(N)<=epsilon*N for all N>=N_epsilon
(possible since R(N)/N->0). Then

    C_B2(F_N) >= (3-epsilon)*N   for EVERY integer N >= N_epsilon,

and liminf C_B2(F_N)/N >= 3. The declared pair (3,1) and the finite
projection in `bound.json` are two renderings of this one theorem and
this one language.

The restriction fixing Q is a lower-bound device, not advice: every
circuit must compute F_N on ALL prefixes, while the uniform evaluator
simply reads and tests whatever prefix it is given.

## 7. Evidence, checks, and honest scope

* `disperser.md` proves Lemma D from first principles (prime-degree zero
  gaps, the permanent argument, block decomposition, coefficient
  isolation).
* `elimination.md` proves Lemma E from first principles (states, moves,
  the flat potential H = 2s+6*(used parameters), two eliminating cases,
  termination).
* `family.py` implements the evaluator of Section 1 verbatim.
* `audit.py` performs finite supplementary checks: exactness of the
  prefix test against brute-force factoring for small prime k; agreement
  of the fast evaluator with a direct subset-sum implementation of
  sigma_t; the guard threshold; the R(N)/N decay; the prime-k zero-run
  uniqueness of `disperser.md` Claim 5 on random subspaces; and
  nonconstancy of the payload on random affine subspaces at the smallest
  non-vacuous parameters (k=37, r=10, t=3, dimension 345 of 370). None
  of these replaces the general proofs; they guard against slips.
* No claimed review, checker success, or prior acceptance is a premise
  of anything above.

This submission is strictly dominated by the organizer's public baseline
at (7/2,1), which remains under review at the time of writing. We claim
only: one fixed uniform NP language, an unconditional 3N-o(N) full-B2
circuit lower bound for it, proved self-containedly.

Sources, for attribution and version identification:

* Ben-Sasson--Kopparty, author preprint 2010-03-12, Lemmas 3.7, 3.10,
  3.13, Theorem 6.1, Lemma 7.2, Theorem 7.3, Claim 7.4,
  https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
  SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf.
* Mishra--Santhanam, ECCC TR26-260, 2026-10-19 (context of the full
  parameter method; not invoked by our proofs),
  https://eccc.weizmann.ac.il/report/2026/260/download/
  SHA256 38e9ff3cea4ac0399723756f00c2425036c438887ab5ae4abb10d09fa971a3bc.
* Li--Yang, ECCC TR21-023, 2021-02-20, Sections 1.1, 2.1, 2.4, 3.3, for
  the protected full-B2 conventions (their LP analysis is not used).
* Schnorr's gate elimination and its textbook treatment in Wegener,
  *The Complexity of Boolean Functions*, 1987, Chapter 5, for the
  classical tradition our `elimination.md` belongs to.
* Bertrand's postulate as stated in Hardy--Wright, *An Introduction to
  the Theory of Numbers*, Thm. 418.
