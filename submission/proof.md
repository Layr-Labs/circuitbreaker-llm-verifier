# A uniform certificate-prefix family: coefficient 4, exponent 1

## Claim and provenance

For the single total Boolean language defined below, with **N denoting the
entire input length**, we prove

    C_B2(F_N) >= 4N - R(N),       R(N)/N -> 0.

The joint pair is (4,1). Coefficient 4 and the evaluator use existing
mathematics; the current contribution is the exact adjacent-order selection in Section 13.
It improves the supplied certified remainder on infinitely many lengths and
weakly improves it at every sufficiently large length. The gate argument is the
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
2. Let B be the smallest positive integer with B^5 > 4*N^3. Let k be
   the smallest prime at least B, found by deterministic trial division.
3. If k>=N, output 0. Otherwise put m=N-k, r=ceil(m/k),
   t_minus=floor(sqrt(r)), t_plus=ceil(sqrt(r)). If t_plus>=r, output 0.
   For each t in {t_minus,t_plus}, compute H,T_t,U,d_t and B_t as in
   Sections 5 and 12. Choose the lexicographically least (8d_t-4B_t-1,t).
   This specifies t uniquely by exact integer arithmetic; put d_N=d_t.
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

Bertrand's prime theorem gives k<2B, so

    4^(1/5)*N^(3/5) < k < 2*4^(1/5)*N^(3/5)+2.

Thus k=Theta(N^(3/5)) and is prime. The runtime and eligibility argument
needs only this elementary size bound, not an effective prime-gap rate.
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

We have k=Theta(N^(3/5)), m=N-k, r=Theta(N^(2/5)),
t=Theta(N^(1/5)). Consequently k<N and t<r for every sufficiently
large integer N. The former guard kt>=r(r+1) is removed: the sharpened
deficit-selection lemma does not need it. These estimates cover all large
integer lengths, not a sparse subsequence.

## 5. Hard restriction and its precise loss

For each sufficiently large N fix any irreducible Q of degree k and fix
the first k bits of x to its coefficients. The remaining function G_Q
has m inputs. Put

    H=r(r+3)/2, T=sum_(j=2)^t max(k-j,0), U=sum_(j=1)^r max(k-j,0),
    d_N = H+max(T,floor(U/(t+1)))+1.

`disperser.md` proves that, for prime field degree k and r>t>=1, the function

    (a_i) -> pi(sum_(|I|=t) product_(i in I) a_i^(2^(r+2-i)-1))

is nonconstant on every affine subspace of K^r of dimension greater than
H+max(T,U/(t+1)), for any nonzero F_2-linear pi. Our pi is the constant
coefficient, with pi(1)=1. The payload embedding is injective and linear,
so G_Q is an affine disperser of dimension d_N on its m actual input bits.
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
The current evaluator retains the certificate-prefix design, uses the direct
threshold B^5>4N^3, and removes the obsolete guard. The proved coefficient
is 4 rather than the original 7/2. This package does not claim
a larger coefficient than the other coefficient-4 candidates or a smaller
remainder than their polylogarithmic remainders. The immediate own predecessor 440a0962-3f37-4edb-8097-86c3b209b5c1
supplies the refined positive-deficit lemma and remainder comparison. Its
publication or review outcome is not a mathematical premise. Its complete
threshold argument is reproduced here. The preceding two-stage degree rule and its resource comparison are
credited to own predecessor dce46c96-6868-4a18-875f-290e1e06b980 and
are now an explicit reference, not the current evaluator. The full-hypothesis
threshold and gate improvement in Section 10 are carried from c42988f5;
the finite four-gate ledger is carried from
3b8a3eb4-444e-4250-a94f-b8e5578f7772. The new contribution
is the adjacent-order selection and its exact certificate comparison in Section 13.
The direct-degree and degree-sensitive terminal arguments are retained material.

## 8. Threshold-halving comparison at the current parameters

The general supplied positive-deficit proof applies to every 1<=t<r.
For the current k,r,t define

    A=kr/t, B=kt, R=r(r+1), H=r(r+3)/2,
    d_old=floor(A)+B+R+1,
    d_base=H+max(B,floor(A))+1.

Both are sufficient dimensions for the same current family: the exact
threshold in Section 5 obeys T<=kt and U/(t+1)<kr/t, so d_N<=d_base,
and d_old>=d_base. Since t=sqrt(r)+O(1), A/B->1. Also r/B->0 and
k/B=1/t->0. Consequently

    2d_base-d_old=|A-B|+2r+O(1)=o(B),
    d_base/d_old -> 1/2.

Section 10 proves d_base-d_N=O(k+r)=o(B), so the exact sharper dimension
has the same half-loss limit. The finite degree term in Section 12 is O(r)
and does not affect this comparison. Thus the complete current loss
L=4k+8d_N-4B_deg-1 satisfies

    L/(4k+8d_old) -> 1/2.

This compares supplied sufficient remainders for the SAME current family,
not actual circuit complexities. It does not compare this family to later
polylogarithmic-dimension constructions. The threshold argument and its
original half-loss comparison are credited to the supplied predecessors.

## 9. Direct degree threshold and leading certified loss

We use the classical unconditional prime number theorem in its standard
form pi(x)~x/ln(x), where pi counts primes <=x. This familiar theorem is
not re-proved here; its exact statement is also recorded in NIST DLMF
Equation 27.2.3, https://dlmf.nist.gov/27.2.E3 (Hadamard/de la Vallee
Poussin). The application needs only the following elementary consequence.
For every fixed eta>0,

    pi((1+eta)x)-pi(x) ~ eta*x/ln(x) > 0

eventually. Hence the least prime >=x is x*(1+o(1)), for all real x
tending to infinity. No special sequence of lengths is selected. Applying
this with B as in Section 1 gives

    k/N^(3/5) -> a_star=4^(1/5)=2^(2/5).

For any degree sequence k/N^(3/5)->a in (0,infinity), the same rounding
calculation used below gives r=N/k+O(1), t=sqrt(r)+O(1), and

    d_base=kt+r(r+3)/2+O(k)
          =sqrt(kN)+N^2/(2k^2)+o(N^(4/5)).

Section 10 shows d_N differs from d_base by O(k+r)=o(N^(4/5)). Thus
for the stronger certified loss L=4k+8d_N-4B_deg-1 from Section 12,
where B_deg=t(t+3)/2=O(r),

    L/N^(4/5) -> 8*F(a),  F(a)=sqrt(a)+1/(2a^2).

For a>0, F'(a)=(a^(5/2)-2)/(2a^3). Hence its unique global minimum
is at a=a_star, with F(a_star)=(5/4)*2^(1/5). For the actual current
family this proves

    L/N^(4/5) -> 10*2^(1/5),
    C_B2(F_N) >= 4N-(10*2^(1/5)+o(1))*N^(4/5).

The optimization concerns this supplied certificate expression among
degree sequences with a positive finite scaling limit. It is not an
optimality claim for true circuit complexities or other constructions.

### 9.1 Explicit preceding-family comparison

The reference uses k0 as the least prime at least the least positive B0
with B0^5>32N^3, then k_pre as the least prime at least ceil(12k0/25).
Otherwise it has the same family algorithm, guard, positive-deficit proof,
sharp dimension and finite ledger. All these arguments are supplied here
and apply to the reference because its degree is prime and Theta(N^(3/5)).
Thus it independently has certified loss L_pre=4k_pre+8d_N(k_pre)-4.
The same prime-number-theorem consequence, applied twice, yields

    k0/N^(3/5)->2, k_pre/N^(3/5)->24/25,
    L_pre/N^(4/5)->8*F(24/25),
    L/L_pre -> ((5/4)*2^(1/5))/(sqrt(24/25)+625/1152).

This limit is strictly below 19/20. Here is a rational certificate without
rounding the irrational quantities: 2^(1/5)<23/20 since 23^5>2*20^5,
and sqrt(24/25)>979/1000 since 25*979^2<24*1000^2. Therefore

    (5/4)*2^(1/5) < 23/16
        < (19/20)*(979/1000+625/1152)
        < (19/20)*F(24/25).

Consequently limsup L/L_pre<=19/20<1, comparing complete certified
remainders including every prefix bit. Neither the functions, fields nor
the actual minimum circuit sizes are equated. The reference is credited
material, not a presumed accepted claim.

## 10. Exact capacity threshold with either adjacent order

The supplied selection proof uses sum delta_i+max d_i>k+1. For arbitrary
1<=t<r its exact dimension is the Section 5 formula. The clipped capacity
sums can be computed without enumeration: with u=min(t,k-1), v=min(r,k-1),

    T_t=(u-1)k-u(u+1)/2+1,
    U=vk-v(v+1)/2.

Here u>=1, and its u=1 case gives the empty sum zero. Eventually k>r+1
and k>t, and t=sqrt(r)+O(1). Hence

    T_t=kt+O(k+r),
    U/(t+1)=kt+O(k+r),
    d_N=H+kt+O(k+r).

For the second estimate, r=t^2+O(t), so kr/(t+1)=kt+O(k), and
r(r+1)/(2(t+1))=O(r^(3/2))=O(k) for the actual direct-degree rule.
For any fixed positive finite degree-scaling limit in Section 9 the same
last comparison holds, with a constant depending on that limit.
Also d_base=H+kt+O(k), so d_base-d_N=O(k+r), nonnegative by Section 8.
This proves the asymptotic comparisons used above.

Crucially, the maximum in d_N need not be attained by T_t when the lower
adjacent order is selected. No all-length assertion of that stronger form
is used here. Section 13 proves T dominance only on its specified infinite
subsequence, where its exact saving calculation needs it.

## 11. Retained finite ledger

Section 3 of elimination.md proves 4m-8d+4 for a dimension-d affine
disperser. Applying it to the fixed-prefix restriction gives

    C_B2(F_N)>=4(N-k)-8d_N+4.

This retained result is credited to own ledger provider
3b8a3eb4-444e-4250-a94f-b8e5578f7772, including its cleanup-(C)
clarification. The stronger degree-sensitive bound below is also retained.
No previous platform judgment is a mathematical premise.

## 12. Use the high-degree coefficient already isolated

Section 9 of disperser.md proves more than nonconstancy: the fixed-prefix
payload G_Q has algebraic degree at least B_deg=t(t+3)/2 on every
d_N-dimensional affine subspace. The isolated reduced field monomial
has total binary weight at least this value; the supplied conversion to
algebraic normal form and affine-composition argument establish the claim.
The whole prefix function is not asserted to have this property.

Eventually B_deg=Theta(N^(2/5)), d_N=Theta(N^(4/5)), and m=N-k is
Theta(N), so 1<=B_deg<=d_N<m. Apply Section 5 of elimination.md to
the same hard prefix restriction, obtaining

    C_B2(F_N)>=4(N-k)-8d_N+4B_deg+1
              =4N-(4k+8d_N-4B_deg-1).

Relative to the preceding same-family finite ledger bound with +4, this
is a certified improvement of 4B_deg-3 gates, tending to infinity. The
new complete loss differs by O(r) from that in Section 11, so the leading
constant and reference ratio in Section 9 and the threshold-halving limit
in Section 8 remain valid. No leading coefficient or exponent gain is
claimed. The direct degree rule and leading-loss proof are inherited from
811df8a3-8fde-441f-aee7-1b18dfcadce0. The high-degree consequence and
terminal gate bound are inherited from the degree-gate package first
uploaded as bef2a7bf-62bb-4961-8fe0-860b62ca879a and retried unchanged as
6fd4cdbe-2d59-4478-9d8c-81812124629f. The former had an operational
provider-server failure, not a mathematical verdict. Neither outcome is
assumed. The current evaluator changes the symmetric order as follows.


## 13. Select the adjacent order by the complete certificate

Let F_up denote the fully specified reference using the same k, prefix
validation, guards and symmetric recurrence, but always t_plus=ceil(sqrt(r)).
This is exactly the preceding degree-gate evaluator. Define d_up and B_up
by the same general formulas evaluated at t_plus, and its certified loss
L_up=4k+8d_up-4B_up-1. The current evaluator chooses the smaller value of
8d_t-4B_t-1 between t_minus and t_plus, breaking ties toward the smaller t.
All values are integers of polynomial bit length and use the clipped sums
in Section 10. Selection is polynomial time and preserves every eligibility
argument, with t=sqrt(r)+O(1) at every sufficiently large actual length.

The general algebraic and gate arguments give the stated loss for either
order. Both options eventually obey 1<=B_t<=d_t<m and 1<=t<r. Since the
reference option is included, the complete current certificate satisfies

    L(N)<=L_up(N) for every sufficiently large integer N.

This compares certified losses for two DIFFERENT fixed languages, not
their actual minimum circuit sizes. It asserts no strict saving at every
length. Its strict improvement is nevertheless unbounded, as we now prove.

The integer sequence r(N)=ceil(N/k(N))-1 tends to infinity, and k(N) is
nondecreasing. Thus r(N+1)<=r(N)+1: if k stays fixed this follows from
ceil((N+1)/k)<=ceil(N/k)+1; if it increases, replacing it by the old k
only increases that upper bound. An unbounded integer sequence whose
upward steps are at most one attains every sufficiently large integer
(take its first crossing of that integer). In particular infinitely many
actual lengths have r=u^2+1 with u tending to infinity.

Along those lengths, Section 9 gives

    k/r^(3/2) -> (2^(2/5))^(5/2)=2, hence k/u^3->2.

Eventually k>r+1, and the two orders are u and u+1. With
U=(u^2+1)k-(u^2+1)(u^2+2)/2 and
T_t=(t-1)k-t(t+1)/2+1, direct division gives

    (T_u-U/(u+1))/k -> 1/4,
    (T_(u+1)-U/(u+2))/k -> 9/4.

For example (u^2+1)/(u+1)=u-1+2/(u+1), whereas the quadratic-capacity
term divided by (u+1)k tends to 1/4; all t(t+1)/(2k) terms vanish.
The corresponding division by u+2 uses u-2+5/(u+2). Both limits are
positive, so both dimension maxima are eventually attained by T_t.
The floors then disappear exactly. Therefore

    d_up-d_minus=k-u-1,
    B_up-B_minus=u+2,
    L_up-L_minus=8k-12u-16.

This is positive eventually, so the current rule selects t_minus, and

    (L_up-L)/k -> 8

on this infinite set of actual lengths. Since k/N^(3/5)->2^(2/5), the
strict saving there is (8*2^(2/5)+o(1))*N^(3/5). The weak comparison
holds at every large length. These are differences of complete proved
certificates, not differences of actual circuit complexities. The native
projection remains (4,1), and Section 9's leading constant is unchanged.
