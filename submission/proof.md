# An explicit all-length bound for the certificate-prefix family

## Claim and provenance

For the single total Boolean language defined below, with **N denoting the
entire input length**, we prove for every integer N>=1

    C_B2(F_N) >= (7/2)N - 35*N^(4/5).

In particular, for every epsilon>0 the sufficient threshold

    N_epsilon = max(1, ceil((35/epsilon)^5))

ensures C_B2(F_N)>=(7/2-epsilon)N for every integer N>=N_epsilon.
The evaluator, including its behavior at N=0, is unchanged from the
predecessor identified in Section 7. The new contribution is an explicit
remainder and all-length threshold, not an increased leading score.

The joint pair is (7/2,1). This applies existing mathematics; it does not
claim discovery of a new circuit lower bound. The gate argument is the
parameter method of Mishra and Santhanam, ECCC TR26-260, October 19, 2026,
Theorem 4.1, pp. 11-28. The algebraic construction is based on Ben-Sasson
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

## 4. An explicit threshold for every guard

For N>=10^5 put a=401/100. Section 3 and N^(3/5)>=1000 give

    2*N^(3/5) < k < 4*N^(3/5)+2 < a*N^(3/5).

Thus k/N<a/N^(2/5)<=a/100<1. Since

    r=ceil(N/k)-1,
    N/k-1 <= r < N/k,

we have r>N^(2/5)/a-1>=100/a-1>23, hence r>=24. For any
integer r>=3, ceil(sqrt(r))<=r-1, so t<r and t>=1.

For the remaining guard, r<N/k and k^5>32*N^3 imply

    k^2 > 32*N^3/k^3 > 32*r^3,
    k*t >= k*sqrt(r) > sqrt(32)*r^2 > 2*r^2 >= r*(r+1).

Therefore every guard passes at every integer N>=10^5. This calculation
uses only the exact choices of B,k,r,t, including the ceiling in r.
It makes no restriction to prime or otherwise selected input lengths.

## 5. Hard restriction and its precise loss

For each integer N>=10^5 fix any irreducible Q of degree k and fix
the first k bits of x to its coefficients. The remaining function G_Q
has m inputs. Put

    d_N = floor(kr/t) + kt + r(r+1) + 1.

`disperser.md` proves that, for prime field degree k, r>t>=1 and
k*t>=r(r+1), the function

    (a_i) -> pi(sum_(|I|=t) product_(i in I) a_i^(2^(r+2-i)-1))

is nonconstant on every affine subspace of K^r of dimension greater than
kr/t+kt+r(r+1), for any nonzero F_2-linear pi. Our pi is the constant
coefficient, with pi(1)=1. The payload embedding is injective and linear,
so G_Q is an affine disperser of dimension d_N on its m actual input bits.
Zero padding preserves dimension and does not add circuit inputs.

Since kr=O(N), r=Theta(N^(2/5)) and t=Theta(N^(1/5)),

    d_N=O(N^(4/5)),          k=O(N^(3/5)).

We use this conservative sublinear disperser dimension, not later
logarithmic-dimension constructions. Either suffices for the leading pair.

Fixing prefix inputs of any circuit for F_N yields a circuit for G_Q with
no additional gates. `elimination.md` proves the full-B2 bound

    |C| >= (7/2)m-7d

for a dimension-d affine disperser. It includes arbitrary sharing,
repeated pins, constants, degenerate gates and terminal-output cases.
Apply it to the restricted circuit (or its no-larger normalized form):

    C_B2(F_N) >= (7/2)(N-k)-7d_N
               = (7/2)N-((7/2)k+7d_N).

The choice of Q in this lower-bound restriction is not advice to the
evaluator: every circuit must work for all prefixes, while the uniform
algorithm simply reads and tests its input.

## 6. Quantitative remainder for every actual length

Let R(N)=(7/2)k+7d_N for N>=10^5. Since sqrt(r)<=t<=sqrt(r)+1,

    d_N <= 2*k*sqrt(r)+k+r^2+r+1,
    R(N) <= 14*k*sqrt(r)+(21/2)*k+7*r^2+7*r+7.

Use r<N/k, k<a*N^(3/5), and k>2*N^(3/5), all proved above:

    R(N) < (14*sqrt(a)+7/4)*N^(4/5)
           +(21*a/2)*N^(3/5)+(7/2)*N^(2/5)+7.

Here a=401/100 and sqrt(a)<801/400, since
(801/400)^2-a=1/160000>0. Also N^(1/5)>=10. Division by N^(4/5)
therefore gives the entirely rational upper estimate

    R(N)/N^(4/5)
      < 5957/200 + 8421/2000 + 7/200 + 7/10000
      = 42539/1250
      < 35.

Together with Section 5 this proves the displayed all-length bound for
N>=10^5. For 1<=N<10^5, N^(1/5)<10 and

    (7/2)*N-35*N^(4/5)
      = N^(4/5)*((7/2)*N^(1/5)-35) < 0.

Nonnegativity of circuit size proves the same bound at these remaining
lengths, whether a guard passes or fails. The evaluator at N=0 is already
defined, but the quantitative inequality is only claimed for N>=1.

Finally N>=max(1,ceil((35/epsilon)^5)) implies
35*N^(-1/5)<=epsilon. This proves the explicit threshold in the claim
and both original asymptotic quantifiers for one fixed family. The two
native coordinates remain exactly (7/2,1). Making the remainder explicit
does not establish any larger coefficient or exponent.

## 7. Predecessor, declared reference and evidence limits

The material predecessor is Yukon submission
`db5f9c6e-607f-424a-89f8-5a49d29d07df`. Its original declared joint
claim is the same total family in P/NP with
C_B2(F_N)>=(7/2)N-O(N^(4/5)) and liminf C_B2(F_N)/N>=7/2.
The present family.py, disperser.md and elimination.md retain that package's
bytes. This submission makes its quantitative remainder explicit, replaces
the eventual guard argument by Section 4, and extends the finite audit.
The predecessor is declared both in prior_results and in the separate
platform contribution-dependency declaration.

An authenticated accepted version of the predecessor's original declared
claim supports the asymptotic component for this identical family. Its
incidental proof prose and old checker outputs are not additional accepted
claims and do not establish the new constant 35. The full arguments needed
here are supplied independently in the current evidence, including Sections
4 and 6 and the two auxiliary proof files. The result is unconditional and
does not assume that the predecessor is available or human-approved. A
platform unavailability notice supplies no mathematical premise or negation.

`disperser.md` contains the coefficient separation and block-mixing proof.
`elimination.md` contains the complete gate argument. `audit.py` supplies
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
