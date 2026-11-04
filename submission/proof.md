# Coefficient 7/2 with a sharper finite bound for the same certificate-prefix family

## Claim and provenance

For the single total Boolean language defined below, with **N denoting the
entire input length**, we prove

    C_B2(F_N) >= (7/2)N - R(N),       R(N)/N -> 0.

The joint pair is (7/2,1). This applies existing mathematics; it does not
claim discovery of a new circuit lower bound. The gate argument is the
parameter method of Mishra and Santhanam, ECCC TR26-260, October 19, 2026,
Theorem 4.1, pp. 11-28. The algebraic construction is based on Ben-Sasson
and Kopparty (BK), *Affine Dispersers from Subspace Polynomials*, author
version March 12, 2010, Sections 3, 6 and 7. The necessary arguments are
supplied in our notation in `elimination.md` and `disperser.md`.

The family, evaluator, gate proof, and most of the algebraic proof are reused
with credit from CircuitBreaker submission f4f2c464-3045-4c32-a6e4-32478c05f343,
immutable commit 8db0d3171ce7d19a48503ae34b5945331da55535 (the certificate-prefix
package). We do not claim those arguments or the coefficient as new. That
package's input-prefix wrapper removes a deterministic irreducible-polynomial
construction premise and explicitly counts its prefix bits.

Our contribution is the block-selection argument in disperser.md section 5.
For exactly the same language, it improves the sufficient disperser dimension
from floor(kr/t)+kt+r(r+1)+1 to

    d_N=floor(r(r+3)/2+max(kr/t,(t-1)(k-2)))+1.

We prove the resulting explicit bound at every integer N>=2^20:

    C_B2(F_N) >= (7/2)(N-k)-7d_N
               >= (7/2)N-17N^(4/5)-18N^(3/5)-6N^(2/5)-7.

This improves the finite lower bound on the same family; it does not increase
the leading coefficient, change the exponent, or improve the v1 comparator
metric, which remains [7/2,1]. The full function need not be an affine
disperser: the proof uses a fixed-prefix restriction. No conjecture or
direct-sum assumption is involved.

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

We give an explicit threshold, N>=2^20. From section 3,

    2N^(3/5)<k<4N^(3/5)+2<=5N^(3/5).

Since N^(2/5)>=256, k<N and m=N-k is positive. The ceiling defining r gives

    r >= N/k-1 > N^(2/5)/5-1 >= 256/5-1 > 3,
    r = ceil(N/k-1) < N/k < N^(2/5)/2 <= N^(2/5).

For any integer r>=3, t=ceil(sqrt(r))<r. Moreover

    r(r+1)/t <= sqrt(r)(r+1)
              <= N^(1/5)(N^(2/5)+1)
              <= 2N^(3/5) < k.

Thus both guards pass for every N>=2^20. These statements require no
information about prime gaps beyond Bertrand's theorem. In addition,
r=Theta(N^(2/5)), t=Theta(N^(1/5)), and k=Theta(N^(3/5)).

## 5. Hard restriction and its precise loss

For each sufficiently large N fix any irreducible Q of degree k and fix
the first k bits of x to its coefficients. The remaining function G_Q
has m inputs. Put

    T_N = r(r+3)/2 + max(kr/t,(t-1)(k-2)),
    d_N = floor(T_N)+1.

`disperser.md` proves that, for prime field degree k and r>t>=1, the function

    (a_i) -> pi(sum_(|I|=t) product_(i in I) a_i^(2^(r+2-i)-1))

is nonconstant on every affine subspace of K^r of dimension greater than
T_N, for any nonzero F_2-linear pi. Our pi is the constant
coefficient, with pi(1)=1. The payload embedding is injective and linear,
so G_Q is an affine disperser of dimension d_N on its m actual input bits.
Zero padding preserves dimension and does not add circuit inputs.

Here is an explicit remainder estimate. Since t>=sqrt(r) and
t-1<sqrt(r), both kr/t and (t-1)(k-2) are at most k*sqrt(r).
Also rk<N: the defining ceiling gives rk<m+k=N. Therefore

    k*sqrt(r)=sqrt(k*rk)<sqrt(kN)<=sqrt(5)*N^(4/5),
    r(r+3)/2 < (1/8)N^(4/5)+(3/4)N^(2/5).

It follows that

    d_N <= (sqrt(5)+1/8)N^(4/5)+(3/4)N^(2/5)+1.             (8)

In particular d_N=O(N^(4/5))=o(N), with every constant explicit.
The dimension assertion is meaningful even if d_N exceeds the payload
dimension: in that case it is vacuous and the gate lower bound is
nonpositive, which causes no difficulty.

Fixing prefix inputs of any circuit for F_N yields a circuit for G_Q with
no additional gates. `elimination.md` proves the full-B2 bound

    |C| >= (7/2)m-7d

for a dimension-d affine disperser. It includes arbitrary sharing,
repeated pins, constants, degenerate gates and terminal-output cases.
Apply it to the restricted circuit (or its no-larger normalized form):

    C_B2(F_N) >= (7/2)(N-k)-7d_N
               = (7/2)N-((7/2)k+7d_N).

Set R(N)=(7/2)k+7d_N for N>=2^20. By (8),

    R(N) <= (7sqrt(5)+7/8)N^(4/5)
            +(35/2)N^(3/5)+(21/4)N^(2/5)+7
         <= 17N^(4/5)+18N^(3/5)+6N^(2/5)+7.              (9)

For the last inequality, sqrt(5)<9/4, so 7sqrt(5)+7/8<17;
also 35/2<18 and 21/4<6. All powers are of the entire actual
input length N. Define R(N) arbitrarily large on the finitely many smaller
lengths if a statement at those lengths is desired. Equation (9) gives
R(N)=O(N^(4/5))=o(N). Choosing Q in the lower-bound restriction is not
advice to the evaluator: every circuit works for every prefix, whereas
the uniform evaluator reads and checks its own input prefix.

The improvement over the credited package is strict in its dimension
threshold: W_r=r(r+3)/2<r(r+1) for r>1, and each member of our maximum
is strictly less than kr/t+kt. Hence T_N<kr/t+kt+r(r+1); indeed the
integer term r(r-1)/2 already supplies a positive gap. The unchanged
circuit lemma consequently proves a strictly stronger finite lower bound
for these same parameters. This is a lower-order improvement, not a
stronger expression under the protected comparator.

For each epsilon>0, take N_epsilon large enough that the guards pass and
R(N)/N<=epsilon for every N>=N_epsilon. This gives exactly
C_B2(F_N)>=(7/2-epsilon)N for all those integers. Both coordinates in
bound.json belong to this one theorem and one language.

## 6. Evidence and limits

`disperser.md` contains the coefficient separation and block-mixing proof.
`elimination.md` contains the complete gate argument. `family.py` is the
unchanged explicit evaluator; its polynomial-time correctness is established
in this proof. This is a proof-only submission with no declared experiment.
The credited package's supplementary audit.py is omitted because none of
its finite computations is used as a premise. No alleged prior acceptance,
review, finite test, or checker success is a mathematical premise.

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
