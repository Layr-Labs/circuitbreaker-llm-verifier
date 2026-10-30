# A uniform certificate-prefix family: coefficient 7/2, exponent 1

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

Set R(N)=(7/2)k+7d_N. Section 5a proves the explicit estimate
R(N)<=38*N^(4/5) and that all guards pass, for EVERY integer N>=112.
Hence R(N)=O(N^(4/5))=o(N), with an explicit constant and threshold. The choice of Q in this lower-bound
restriction is not advice to the evaluator: every circuit must work for
all prefixes, while the uniform algorithm simply reads and tests its input.

For each epsilon>0, take N_epsilon large enough that the guards pass and
R(N)/N<=epsilon for every N>=N_epsilon. This gives exactly
C_B2(F_N)>=(7/2-epsilon)N for all those integers. Both coordinates in
bound.json belong to this one theorem and one language.

## 5a. Explicit remainder and threshold

**Claim.** For every integer N>=112 all guards of Section 1 pass (k<N, t<r,
k*t>=r(r+1)) and R(N)=(7/2)k+7d_N<=38*N^(4/5). Consequently

    C_B2(F_N) >= (7/2)N - 38*N^(4/5)    for every integer N>=112.

*Range 112<=N<10^7.* `explicit_check.py` uses the evaluator's exact integer
parameter computation (B, k, m, r, t as in Section 1; primality by trial
division; integer ceilings and isqrt). B, hence k, and r=ceil((N-k)/k) are
nondecreasing step functions of N, and t, the guards t<r and k*t>=r(r+1),
and d_N depend only on (k,r). On a maximal run of N with fixed (k,r), the
remaining conditions k<N and R<=38*N^(4/5) only become easier as N grows.
The script therefore checks the first N of every such run exactly; the runs
are generated consecutively (next run starts at k*r+k+1 or where B changes)
and tile the interval with no gap. It checks k<N, t<r, k*t>=r(r+1) and the
integer inequality (2R)^5<=76^5*N^4, equivalent to R<=38*N^(4/5) since 2R is
a nonnegative integer. It reports zero failures over 32,435 runs; the
largest ratio R/N^(4/5) is 25.78, at N=112. This covers every integer of the
declared interval, not a sample.

*Range N>=10^7.* Put a=N^(1/5)>=25.1. All inequalities below are of the form
(polynomial of lower degree)<=c*a^j with the left side divided by a^j
decreasing in a, so checking them at a=25 proves them for all a>=25.

1. B is least with B^5>32N^3=(2a^3)^5, so 2a^3<B<=2a^3+1.
2. Bertrand's postulate gives a prime in [B,2B), so 2a^3<k<2B<=4a^3+2<=4.01a^3
   (2<=0.01a^3). In particular k<a^5=N, so the first guard passes.
3. m=N-k and r=ceil(m/k). From m/k<N/(2a^3)=a^2/2 we get r<=a^2/2+1; from
   m/k>(a^5-4.01a^3)/(4.01a^3) we get r>a^2/4.01-1>=(a/2.01)^2 (true at a=25:
   a^2(1/4.01-1/4.0401)>=1). Hence sqrt(r)>=a/2.01 and r>=155.
4. t=ceil(sqrt(r)) satisfies sqrt(r)<=t<sqrt(r)+1<r, so t<r passes.
5. k*t>2a^3*sqrt(r)>=0.995*a^4, while r(r+1)<(a^2/2+1)(a^2/2+2)
   =a^4/4+1.5a^2+2<=0.2526*a^4 (1.5a^2+2<=0.0026a^4). So k*t>=r(r+1) passes.
6. d_N<=kr/t+kt+r(r+1)+1. Since r<m/k+1, kr<m+k=N, so kr/t<=N/sqrt(r)
   <=a^5/(a/2.01)=2.01a^4. Next sqrt(r)<=sqrt(a^2/2+1)<=a/sqrt(2)+0.03, so
   t<=0.7072a+1.03 and kt<=4.01a^3(0.7072a+1.03)<=2.836a^4+4.131a^3<=3.002a^4 (4.131a^3<=0.166a^4).
   With item 5 and 1<=0.0001a^4: d_N<=(2.01+3.002+0.2526+0.0001)a^4<=5.265a^4.
7. R=(7/2)k+7d_N<=14.04a^3+36.86a^4<=(14.04/25+36.86)a^4<=37.43a^4<=38*N^(4/5).

The two ranges cover every integer N>=112. For N<112 the language is still
defined (guarded outputs) and no bound is claimed. The scored liminf
statement follows because 38*N^(4/5)/N tends to zero.

## 6. Evidence and limits

`disperser.md` contains the coefficient separation and block-mixing proof.
`elimination.md` contains the complete gate argument. `explicit_check.py` is the complete enumeration used in Section 5a.
`audit.py` supplies
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
