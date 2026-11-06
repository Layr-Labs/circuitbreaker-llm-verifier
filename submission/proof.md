# A uniform certificate-prefix family: coefficient 4, exponent 1

## Claim and provenance

For the single total Boolean language defined below, with **N denoting the
entire input length**, we prove

    C_B2(F_N) >= 4N - R(N),       R(N)/N -> 0.

The joint pair is (4,1). Coefficient 4 and the evaluator use existing
mathematics; the current contribution is the two-prime selection in Section 14.
It improves the supplied adjacent-order certificate on infinitely many actual
lengths and weakly improves it at every sufficiently large length. The gate argument is the
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
2. Let B be the smallest positive integer with B^5 > 4*N^3. Find
   k_plus, the smallest prime >=B, and, if it exists, k_minus, the
   largest prime <B, by deterministic trial division.
3. A degree-k branch is eligible when k<N and, with m=N-k and
   r=ceil(m/k), t_plus=ceil(sqrt(r))<r. If the k_plus branch is not
   eligible, output 0. For each eligible branch k in {k_minus,k_plus},
   consider t_minus=floor(sqrt(r)) and t_plus. Compute d_t and B_t by
   Sections 5 and 12, and L(k,t)=4k+8d_t-4B_t-1. Select the
   lexicographically least (L(k,t),k,t). This uniquely specifies k,m,r,t
   and d_N=d_t by exact integer arithmetic. Absent k_minus is ignored.
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

Bertrand's prime theorem gives k_plus<2B and k_minus>B/2 when the latter
exists (the bounded B cases are checked directly). Thus the selected k obeys

    (4^(1/5)/2)*N^(3/5) < k < 2*4^(1/5)*N^(3/5)+2.

Thus k=Theta(N^(3/5)) and is prime. The runtime and eligibility argument
needs only this elementary size bound, not an effective prime-gap rate.
Both prime scans and trial division are polynomial in N. B is computed by
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
is the two-prime selection and its exact certificate comparison in Section 14.
The adjacent-order rule and its complete comparison are retained from own
5629bae8-dba0-4c3b-8569-722a079fe560, with the relevant arguments supplied
below. No source publication or judgment is a mathematical premise.

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
tending to infinity. No special sequence of lengths is selected. Similarly pi(x)-pi((1-eta)x)>0 eventually for fixed eta in (0,1),
so the largest prime <x is x*(1+o(1)) as x tends to infinity. Strict
endpoints do not matter: the interval ((1-eta)x,x) contains primes
by the same asymptotic, including when x itself is prime. Both branch
degrees therefore satisfy k_plus/B->1 and k_minus/B->1. Selecting
either branch at each length preserves this limit. With B as in Section 1,
this gives

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
assumed. The retained adjacent-order reference is described next; the current two-prime
selection is in Section 14.


## 13. Retained adjacent-order comparison at the upper prime

For this section ONLY, let k=k_plus and r=ceil((N-k_plus)/k_plus). These
are reference parameters, not necessarily the degree and block count
selected by the current evaluator. Define F_up using this upper prime,
the same prefix validation and guards, and always t_plus=ceil(sqrt(r)).
Define F_adj using the same upper prime and the smaller certificate among
t_minus and t_plus, with ties toward smaller t. These are respectively
the supplied degree-gate and adjacent-order predecessors. Each is a fixed
total P language and has its certificate by the general arguments above.
Let their complete losses be L_up and L_adj. Inclusion gives L_adj<=L_up
at every sufficiently large actual length.

The integer reference sequence r(N)=ceil(N/k_plus(N))-1 tends to infinity.
Since k_plus(N) is nondecreasing, r(N+1)<=r(N)+1. An unbounded integer
sequence with upward steps at most one attains every sufficiently large
integer by its first-crossing property. Thus infinitely many actual lengths
have this reference r=u^2+1 tending to infinity.

Along those lengths Section 9 gives k_plus/u^3->2. Eventually k_plus>r+1,
and the two orders are u and u+1. With

    U=(u^2+1)k_plus-(u^2+1)(u^2+2)/2,
    T_t=(t-1)k_plus-t(t+1)/2+1,

exact division gives

    (T_u-U/(u+1))/k_plus ->1/4,
    (T_(u+1)-U/(u+2))/k_plus ->9/4.

For example (u^2+1)/(u+1)=u-1+2/(u+1), whereas the capacity correction
divided by (u+1)k_plus tends to1/4. The analogous division by u+2 uses
u-2+5/(u+2). Both limits are positive, so both dimension maxima are
attained by T_t eventually on this set, and the floors disappear exactly.
Therefore

    d_up-d_minus=k_plus-u-1,
    B_up-B_minus=u+2,
    L_up-L_adj=8k_plus-12u-16 ~8k_plus.

All these exact equalities concern the two upper-prime REFERENCES. The
current loss L in Section 14 is at most L_adj, so its saving from L_up
is AT LEAST the displayed value, not asserted equal to it. In particular
the infinite N^(3/5)-order improvement over F_up remains available.
Comparisons here and below are between complete proved certificates for
DIFFERENT fixed languages, not actual minimum circuit complexities.

## 14. Compare both neighboring prime degrees

The exact evaluator includes the upper-prime adjacent-order branch just
defined. Thus at every sufficiently large integer length its complete
certificate satisfies

    L<=L_adj<=L_up.

No strict saving at every length is asserted. We now prove an additional
unbounded improvement over L_adj at infinitely many actual lengths,
using only the standard PNT already stated and elementary counting.
No uniform short-interval prime theorem is assumed.

### 14.1 Infinitely many short gaps with unchanged block count

For each sufficiently large odd prime q let p be its next prime and put

    N_q=ceil((q^5/4)^(1/3)),       c=4^(-1/3).

The equality 4N_q^3=q^5 is impossible for odd q by its power of two.
The rounding error is O(1), so

    q^5<4N_q^3<(q+1)^5

for all sufficiently large q: the left quantity above q^5 differs by
O(q^(10/3)), while (q+1)^5-q^5 is asymptotic to5q^4. Consequently
B(N_q)=q+1, k_minus=q, k_plus=p, and N_q/q=c*q^(2/3)+O(1/q).
These are actual integer input lengths, not an added promise on inputs.

Fix a real X tending to infinity and consider primes q in [X,2X].
The PNT gives (1+o(1))*X/log(X) such primes. Their consecutive gaps p-q
sum to at most3X: the sum telescopes and the next prime after the last
one is below4X by Bertrand. Therefore at most X/(4log(X)) of these gaps
exceed12log(X). At least a positive constant times X/log(X) primes have

    Delta=p-q<=12log(X).

Put epsilon_X=100log(X)/X^(1/3), which eventually is less than1/4.
For every such short-gap q, the estimates above give

    0<=N_q/q-N_q/p=N_q*Delta/(pq)<=epsilon_X,
    |N_q/q-c*q^(2/3)|<=1/q<epsilon_X

for all sufficiently large X, uniformly for q in [X,2X]. The constant100
is more than12*c*2^(5/3), so the first bound follows directly from
N_q<=c*(2X)^(5/3)+1 and p,q>=X.

Only o(X/log(X)) INTEGER q in [X,2X] have

    dist(c*q^(2/3), Z)<=2epsilon_X.

Here is a direct count. The increasing function f(x)=c*x^(2/3) has
derivative at least a positive constant times X^(-1/3) on [X,2X].
For each integer j, the preimage of [j-2epsilon_X,j+2epsilon_X]
therefore has length O(epsilon_X*X^(1/3))=O(log(X)). At most O(X^(2/3))
integers j meet f([X,2X]) plus its endpoints. Each preimage interval
contains at most its length plus one integers. Their total is
O(X^(2/3)log(X))=o(X/log(X)), as claimed. This count includes all bad
primes without invoking any equidistribution theorem for primes.

There are thus short-gap primes q outside this exceptional set for every
sufficiently large X. For each of them, N_q/q and N_q/p lie in the same
open interval between consecutive integers: f(q) is more than2epsilon_X
from either endpoint, its rounding displacement is below epsilon_X, and
the displacement between the two ratios is at most epsilon_X. Hence

    ceil(N_q/q)=ceil(N_q/p),
    r_minus=r_plus=R.

Taking such q in unbounded dyadic ranges proves infinitely many distinct
actual lengths with this common block count.

### 14.2 The exact complete loss decreases

At one of these lengths let t be the upper-prime reference's selected
adjacent order. It is also an available lower-prime order because R is
the same. Eventually q>R+1 and t tends to infinity. The unclipped sums
satisfy

    T_t(p)-T_t(q)=(t-1)Delta,
    U(p)-U(q)=R*Delta.

Both adjacent orders obey R>=(t-1)^2, so R/(t+1)>=t-3. For real a>=b,
floor(a)-floor(b)>=floor(a-b). Therefore

    floor(U(p)/(t+1))-floor(U(q)/(t+1)) >=(t-3)Delta.

Each entry of the maximum defining d_t falls by at least (t-3)Delta
when p is replaced by q. The same holds for the maximum and hence

    d_t(p)-d_t(q)>=(t-3)Delta.

The degree B_t=t(t+3)/2 is unchanged. Including every prefix bit gives

    L_adj-L(q,t)>=4Delta+8(t-3)Delta=(8t-20)Delta.

The current minimum is at most L(q,t), so the same lower bound holds
for L_adj-L. For large odd primes Delta>=2, and t=Theta(N_q^(1/5)),
because q~4^(1/5)*N_q^(3/5) and R~N_q/q. Thus

    L_adj-L>=(8t-20)Delta=Omega(N_q^(1/5))>0

on infinitely many actual lengths. At all large lengths weak dominance
holds by inclusion. No asymptotic equality in terms of the gap, no strict
all-length saving, and no comparison of actual circuit sizes is claimed.
All leading limits and the native pair (4,1) remain as proved above.
