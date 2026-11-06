# A uniform certificate-prefix family: coefficient 4, exponent 1

## Claim and provenance

For the single total Boolean language defined below, with **N denoting the
entire input length**, we prove

    C_B2(F_N) >= 4N - R(N),       R(N)/N -> 0.

The joint pair is (4,1). Coefficient 4 and the evaluator use existing
mathematics; the contribution here is the sharper supplied block-selection
and dimension estimate, and its resulting certified remainder. The gate argument is the
parameter method of Mishra and Santhanam, ECCC TR26-260, revision 1, November 3, 2026,
Sections 2-4. The algebraic construction is based on Ben-Sasson
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
2. Let B be the smallest positive integer with B^5 > 32*N^3. Let k be the
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

    2*N^(3/5) < k < 4*N^(3/5)+2.

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
Consequently k<N and t<r eventually. Also r<N/k+1<(1/2)N^(2/5)+1 and

    r(r+1)/t <= sqrt(r)(r+1)
              <= sqrt((1/2)N^(2/5)+1)*((1/2)N^(2/5)+2).

Dividing the last expression by N^(3/5) gives limit 1/(2*sqrt(2)),
whereas k/N^(3/5)>2. Thus k*t>=r(r+1) eventually. These estimates hold
for every sufficiently large integer N, not just selected lengths.

## 5. Hard restriction and its precise loss

For each sufficiently large N fix any irreducible Q of degree k and fix
the first k bits of x to its coefficients. The remaining function G_Q
has m inputs. Put

    d_N = kt + r(r+3)/2 + 1.

`disperser.md` proves that, for prime field degree k, r>t>=1 and
k*t>=r(r+1), the function

    (a_i) -> pi(sum_(|I|=t) product_(i in I) a_i^(2^(r+2-i)-1))

is nonconstant on every affine subspace of K^r of dimension greater than
max(kt,kr/t)+r(r+3)/2, for any nonzero F_2-linear pi. Our pi is the constant
coefficient, with pi(1)=1. The payload embedding is injective and linear,
and kt>=kr/t because t^2>=r. Thus G_Q is an affine disperser of dimension d_N on its m actual input bits.
Zero padding preserves dimension and does not add circuit inputs.

Since kr=O(N), r=Theta(N^(2/5)) and t=Theta(N^(1/5)),

    d_N=O(N^(4/5)),          k=O(N^(3/5)).

We use this conservative sublinear disperser dimension, not later
logarithmic-dimension constructions. Either suffices for the leading pair.

Fixing prefix inputs of any circuit for F_N yields a circuit for G_Q with
no additional gates. `elimination.md` proves the full-B2 bound

    |C| >= 4m-8d

for a dimension-d affine disperser. It includes arbitrary sharing,
repeated pins, constants, degenerate gates and terminal-output cases.
Apply it to the restricted circuit (or its no-larger normalized form):

    C_B2(F_N) >= 4(N-k)-8d_N
               = 4N-(4k+8d_N).

Set R(N)=4k+8d_N eventually and extend it on the finitely many guarded
lengths. Then R(N)=O(N^(4/5))=o(N). The choice of Q in this lower-bound
restriction is not advice to the evaluator: every circuit must work for
all prefixes, while the uniform algorithm simply reads and tests its input.

For each epsilon>0, take N_epsilon large enough that the guards pass and
R(N)/N<=epsilon for every N>=N_epsilon. This gives exactly
C_B2(F_N)>=(4-epsilon)N for all those integers. Both coordinates in
bound.json belong to this one theorem and one language.

## 6. Evidence and limits

`disperser.md` contains the coefficient separation and block-mixing proof.
`elimination.md` contains the complete gate argument. `audit.py` supplies
supplementary finite checks of Boolean identities, rational estimates and
the prefix test; these do not replace either universal proof. No claimed
review, checker success or prior acceptance is a mathematical premise.

Sources, for attribution and version identification:

* Mishra--Santhanam, ECCC TR26-260 revision 1, 2026-11-03, Sections 2-4,
  https://eccc.weizmann.ac.il/report/2026/260/download/
  The revision is identified by its date and title; the older version hash does not identify this revision.
* Ben-Sasson--Kopparty, author preprint 2010-03-12, Lemmas 3.7, 3.10, 3.13,
  Theorem 6.1, Lemma 7.2, Theorem 7.3 and Claim 7.4,
  https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
  SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf.
* Li--Yang, ECCC TR21-023, 2021-02-20, Sections 1.1, 2.1, 2.4 and 3.3,
  for the protected full-B2 conventions. Their LP and lower-bound case
  analysis are not invoked by this submission.

## 7. Reused submission material

The evaluator, certificate-prefix family, algebraic disperser proof and bounded
audit are adapted from submission db5f9c6e-607f-424a-89f8-5a49d29d07df
(candidate commit 02e21edbaa931b85f664acd2da6873155dc57716).
The full gate argument in elimination.md, including its 4n-8d finite refinement,
is reproduced from submission 28feddfa-dfc6-412a-a507-a525b15567fe
(candidate commit a28162b3c8200880a6cabcdb07572205030bc491).
All premises are supplied here; neither source's review outcome is a premise.
The resulting fixed family is exactly the source certificate-prefix family,
but the proved coefficient is 4 rather than 7/2. This package does not claim
a larger coefficient than the other coefficient-4 candidates or a smaller
remainder than their polylogarithmic remainders. Its purpose is a compact,
complete algebraic route to the coefficient-4 theorem.

## 8. Sharper certified remainder for the same family

The new block-selection proof in disperser.md Section 5 selects by positive
deficit d_i-w_i. This replaces the former sum of coarse thresholds by
max(kt,kr/t)+H, where H=sum_i w_i=r(r+3)/2. Since t=ceil(sqrt(r)),
kt>=kr/t, so our integer d_new=kt+H+1 strictly exceeds that threshold.
The evaluator and its guards are unchanged. Define

    A=kr/t, B=kt, R=r(r+1), H=r(r+3)/2,
    d_old=floor(A)+B+R+1, d_new=B+H+1.

The former sufficient threshold is also proved by the supplied argument:
it is no smaller than the new threshold eventually, as shown below. No
prior acceptance is assumed for either assertion. Both are O(N^(4/5)).
As r tends to infinity, A/B=r/t^2 tends to 1. Moreover H=R/2+r and
r/B=r/(kt)<=sqrt(r)/k tends to zero, using r=O(N^(2/5)) and
k=Theta(N^(3/5)). Floors and constants are o(B). Therefore

    2*d_new-d_old = B-A+2r+O(1) = o(B),
    d_old>=B, and hence d_new/d_old -> 1/2.

This proof does not require convergence of the prime ratio k/N^(3/5).
In particular d_new<d_old for every sufficiently large integer length.
The exact certified circuit losses are L_new=4k+8d_new and
L_old=4k+8d_old. Since k/B=1/t tends to zero,

    L_new/L_old -> 1/2 < 1.

This compares two proved lower-bound remainders for the SAME fixed family,
not actual minimum circuit complexities. It improves the earlier compact
algebraic package's certified remainder while retaining coefficient 4 and
exponent 1. It does not claim to improve the stronger polylogarithmic
remainders supplied for other families.
