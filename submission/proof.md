# Coefficient 7/2 with a fully explicit remainder

## Claim and provenance

For the single total Boolean language defined below, with **N denoting the
entire input length** and s=N^(1/5), we prove that for EVERY integer
N>=16807=7^5

    C_B2(F_N) >= (7/2)N - E(N),
    E(N) = (10927/300)N^(4/5) + 42N^(3/5) + (7/2)N^(2/5)
           + (119/24)N^(1/5) + 35/2.

This is the version-2 `eventual_lower_bound` in bound.json: leading term
(7/2)N, explicit remainder with the five monomials of E, threshold 16807.
Since E(N)/N tends to zero it implies liminf C_B2(F_N)/N >= 7/2, which is
the finite projection (7/2,1). The remainder, not only the liminf, is the
declared theorem, and section 5 proves it with explicit constants.

This applies existing mathematics; it claims no new circuit lower bound.
The gate argument is the parameter method of Mishra and Santhanam, ECCC
TR26-260, October 19, 2026, Theorem 4.1. The algebraic construction is
Ben-Sasson and Kopparty (BK), *Affine Dispersers from Subspace Polynomials*,
author version March 12, 2010, Sections 3, 6 and 7. The necessary arguments
are supplied in `elimination.md` and `disperser.md`. The family, those two
files, family.py and audit.py are unchanged from Yukon submission
f4f2c464-3045-4c32-a6e4-32478c05f343, which stated only an unspecified
O(N^(4/5)) loss. The new content is the explicit, all-length parameter
analysis of sections 4 and 5 and the exact check in remainder.py.

The input-prefix wrapper avoids any deterministic irreducible-polynomial
construction. It explicitly pays for all prefix bits. The **whole** function
is not an affine disperser: only a fixed-prefix restriction is asserted to
be one. No conjecture, restricted circuit basis, direct-sum assumption or
exponent gain is used. We assert no exponent above 1 and no endpoint.

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

## 4. Explicit parameter bounds for every N>=16807

Fix an integer N>=16807 and put s=N^(1/5), so s>=7 and N=s^5.
All statements below hold for every such integer N, not for a subsequence.

(P1) Since (2s^3)^5=32N^3, the least positive integer B with B^5>32N^3
is B=floor(2s^3)+1. Hence 2s^3<B<=2s^3+1, and B>=2.

(P2) k is the least prime with k>=B, so k>=B>2s^3. By Bertrand's postulate
(Chebyshev's theorem: for every integer B>=2 some prime p has B<p<2B),
k<2B, so k<=2B-1<=4s^3+1.

(P3) k<N, because 4s^3+1<s^5 for s>=7. Thus m=N-k>=1 and r=ceil(m/k)
satisfies m/k<=r<m/k+1. Multiplying by k gives

    kr < m+k = N,              r < N/k < s^5/(2s^3) = s^2/2.

For the lower bound, r>=m/k=N/k-1>=s^5/(4s^3+1)-1. We claim this is at
least (9/40)s^2. Clearing the positive denominator 40(4s^3+1), the claim
is 4s^5-160s^3-9s^2-40>=0. For s>=7 we have s^2-40>=9, so
4s^5-160s^3=4s^3(s^2-40)>=36s^3>=9s^2+40. Hence

    r >= (9/40)s^2 >= 441/40 > 11,     so r>=12.

(P4) t=ceil(sqrt(r)) satisfies sqrt(r)<=t<sqrt(r)+1, so t>=1. As r>=3,
sqrt(r)+1<=r, so t<r. Since (17/24)^2=289/576>1/2,

    t < sqrt(s^2/2)+1 = s/sqrt(2)+1 < (17/24)s+1.

(P5) The guard k*t>=r(r+1) holds: k*t>=k*sqrt(r)>2s^3*sqrt(r), while

    r(r+1) = sqrt(r)*sqrt(r)*(r+1) < sqrt(r)*(s/sqrt(2))*(s^2/2+1)
           = sqrt(r)*(s^3/(2*sqrt(2)) + s/sqrt(2)) <= sqrt(r)*2s^3,

because s/sqrt(2)<=s^3 and 1/(2*sqrt(2))+1<2.

Thus every integer N>=16807 passes all guards of section 1: N>0, k<N,
t<r and k*t>=r(r+1). (Exact computation with family.py shows that the
guards also pass for every 112<=N<16807; this finite fact is not used.)

## 5. Hard restriction and the explicit loss

For N>=16807 fix any irreducible Q of degree k (it exists by section 2)
and fix the first k input bits to its coefficients. The remaining function
G_Q has m=N-k inputs. Put

    d_N = floor(kr/t) + kt + r(r+1) + 1.

`disperser.md` proves that, for prime field degree k, r>t>=1 and
k*t>=r(r+1), the function

    (a_i) -> pi(sum_(|I|=t) product_(i in I) a_i^(2^(r+2-i)-1))

is nonconstant on every affine subspace of K^r of dimension greater than
D0=kr/t+kt+r(r+1). Section 4 verified these hypotheses. Our pi is the
constant coefficient, a nonzero functional with pi(1)=1, and the exponents
are those of section 1 step 6, since 2^(e_i+1)-1=2^(r+2-i)-1. The zero-padded
payload embedding is injective and F_2-linear, so it maps an affine
subspace of F_2^m of dimension d_N to one of K^r of the same dimension,
and d_N=floor(D0)+1>D0, because kt and r(r+1) are integers and so
floor(D0)=floor(kr/t)+kt+r(r+1). Hence G_Q is an affine disperser of
dimension d_N on its m actual input bits. Padding adds no circuit inputs.

Fixing prefix inputs of any circuit for F_N to constants yields a circuit
for G_Q with no additional gates. The fixed prefix inputs become constant
sources; under the pinned conventions constants cost zero, and gates reading
them are removed by deletion D1 of elimination.md (fixed pins), at no cost.
`elimination.md` proves the full-B2 bound |C|>=(7/2)n-7d for every n-bit
affine disperser of dimension d, for every n and d, including arbitrary
sharing, repeated pins, constants, degenerate gates and terminal-output
cases. Applied to the restricted circuit,

    C_B2(F_N) >= (7/2)(N-k) - 7d_N = (7/2)N - R(N),
    R(N) = (7/2)k + 7d_N.

**Explicit estimate of d_N.** Using (P3) and (P4):

* floor(kr/t) <= kr/t < N/t <= N/sqrt(r) <= s^5/(s*sqrt(9/40))
  = (2*sqrt(10)/3)s^4 < (53/25)s^4, since (2*sqrt(10)/3)^2=40/9<2809/625.
* kt < (4s^3+1)((17/24)s+1) = (17/6)s^4 + 4s^3 + (17/24)s + 1.
* r(r+1) < (s^2/2)(s^2/2+1) = (1/4)s^4 + (1/2)s^2.

Adding these and 1, with 53/25+17/6+1/4=(636+850+75)/300=1561/300,

    d_N < (1561/300)s^4 + 4s^3 + (1/2)s^2 + (17/24)s + 2.

**Explicit loss.** By (P2), (7/2)k <= 14s^3 + 7/2. Therefore

    R(N) < (10927/300)s^4 + 42s^3 + (7/2)s^2 + (119/24)s + 35/2 = E(N),

using 7*1561=10927, 7*4+14=42, 7/2, 7*17/24=119/24 and 7*2+7/2=35/2.
Combining, for every integer N>=16807,

    C_B2(F_N) >= (7/2)N - R(N) > (7/2)N - E(N).

The right side is negative for N below roughly 2*10^5, so the bound is
vacuous (but still true) there; it is informative for larger N.

This is exactly the declared eventual lower bound with leading term
(7/2)N^1 (zero log and loglog exponents) and the explicit remainder
monomials (10927/300)N^(4/5), 42N^(3/5), (7/2)N^(2/5), (119/24)N^(1/5),
(35/2)N^0, strictly descending in power and each below the leading power.

**Finite projection.** For N>=1 each monomial of E is at most its
coefficient times N^(4/5), so E(N)/N <= S*N^(-1/5) with coefficient sum
S=10927/300+42+7/2+119/24+35/2<105. For each epsilon>0 let N_epsilon be
the larger of 16807 and ceil((105/epsilon)^5). Then
C_B2(F_N)>=(7/2-epsilon)N for EVERY integer N>=N_epsilon. Hence liminf
over all integer N of C_B2(F_N)/N >= 7/2, the projection (7/2,1). Both
coordinates and the remainder belong to this one theorem and this one
language. The choice of Q in the lower-bound restriction is not advice to
the evaluator: every circuit must work for all prefixes, while the uniform
algorithm simply reads and tests its input.

## 6. Evidence and limits

`disperser.md` contains the coefficient separation and block-mixing proof.
`elimination.md` contains the complete gate argument. `audit.py` supplies
supplementary finite checks of Boolean identities, rational estimates and
the prefix test. `remainder.py` is the single declared experiment: it runs
audit.py, recomputes k,r,t,d_N with family.py, and checks R(N)<E(N) exactly
with rational arithmetic for every N in [16807,100000] and for a geometric
sequence of lengths up to about 1.09e15; it also confirms that every guard
failure is below 112. These finite checks
do not replace the universal proofs above. No claimed review, checker
success or prior acceptance is a mathematical premise.

Sources, for attribution and version identification:

* Mishra--Santhanam, ECCC TR26-260, 2026-10-19, Theorem 4.1 and Sections 2-4,
  https://eccc.weizmann.ac.il/report/2026/260/download/
  SHA256 38e9ff3cea4ac0399723756f00c2425036c438887ab5ae4abb10d09fa971a3bc
  (as recorded by submission f4f2c464).
* Ben-Sasson--Kopparty, author preprint 2010-03-12, Lemmas 3.7, 3.10, 3.13,
  Theorem 6.1, Lemma 7.2, Theorem 7.3 and Claim 7.4,
  https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
  SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf.
* Li--Yang, ECCC TR21-023, 2021-02-20, Sections 1.1, 2.1, 2.4 and 3.3,
  for the protected full-B2 conventions. Their LP and lower-bound case
  analysis are not invoked by this submission.
* Bertrand's postulate (Chebyshev, 1852): for every integer n>=2 there is
  a prime p with n<p<2n.
