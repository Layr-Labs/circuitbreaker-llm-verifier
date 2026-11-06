# Explicit11.54 remainder at every integer input length

## 1. One fixed total language at every actual length

At integer N>=0 enumerate every prime2<=k<=N. Let r0=ceil(N/k) and
consider r=r0 and r=r0-1, retaining only r>=3,k>r+1. Put
m=min(N,rk) and p=N-m. Then r=ceil(m/k),0<=p<=k, and all eligible
k are odd primes at least5. Enumerate every integer2<=t<r and put

    H=r(r+3)/2, T=(t-1)k-t(t+1)/2+1, B=t(t+3)/2.

For each positive count s=t,...,r and A=t-1,...,k-s-2 define

    X=min(floor(A/(t-1)),floor((k-t-A)/2),
          floor(((t-1)k-t(t-1)/2-tA)/(t-1))).

Skip X<1; otherwise set C=k-A,J=min(s,C-2X), and

    e_(s,A)=A+(J-t+1)X+floor((C-J-1)^2/4)-floor((C-s-1)^2/4).

For1<=s<t put

    a_s=min(k-s-1,floor((2sk-s(s+1))/(2(s+1)))),
    M_0=max(0,T-r+t-1).

Let E0 be the maximum of M_0, every a_s with a_s>=s, and all the listed
e_(s,A). Empty ranges add nothing. Let d=H+E0+1. Require B<=d<m.
Select the least key(4p+8d-4B-3,k,t,p), together with the explicit zero
option(4N,0,0,0). Let L_N be its first component. The zero option returns0
on the whole length; in particular N=0 is defined. Other choices ignore
the first p actual bits and use the remaining m bits as the payload.

The threshold calculation in family.py uses the equivalent one-A scan
proved in disperser.md Section13. Its small-positive branches are exactly
as above. The selector skips only pairs with lower dimension H+M_0+1>=m
or a lower loss key already at least its incumbent; both follow from
E0>=M_0. These exact skips do not alter the minimum or tie rule.

For a nonzero choice construct Q_k by family.py's irreducible_prime(k),
whose complete deterministic algorithm and correctness are proved in
field.md Sections1--3. Use K=F2[Z]/Q_k in its polynomial basis. Pad the
m payload bits by rk-m zeros, obtaining blocks a_i in K, first bit the
constant coordinate. Put y_i=a_i^(2^(r+2-i)-1), and output the constant
coordinate of the order-t elementary symmetric polynomial in the y_i.
Its division-free recurrence starts E_0=1, all other E_j=0, and updates
E_j<-E_j+y_i E_(j-1) in descending j. Splitting subsets by whether they
contain i proves this recurrence. The field is a deterministic function
of the selected degree, independent of every target input bit.

The two r choices are an exact optimization over all ignored prefixes
0<=p<=k. For any such prefix, m=N-p lies in[N-k,N], so its ceiling
block count is r0 or r0-1. At fixed k,r,t, the threshold and degree do
not depend on m within((r-1)k,rk]. Taking m=min(N,rk) maximizes m, preserves
admissibility and minimizes p and the loss4p+8d-4B-3. Thus no other prefix
length can improve the selected key. Every pair with r>=3 has N>=m>2k,
so k<=N does not omit any eligible choice.

## 2. Uniform P, elementary field facts and NP

In an algebraic closure of F2, the polynomial Z^(2^u)-Z has exactly2^u
distinct roots, since its derivative is1. Those roots are closed under
addition, multiplication and nonzero inversion, and form a field.
Frobenius z->z^2 permutes them. A Frobenius orbit gives a monic irreducible
polynomial over F2 by taking the product of its linear factors: the
coefficients are fixed by Frobenius, and a proper binary factor would
require a proper invariant subset of the orbit. Conversely the roots of
an irreducible polynomial form one orbit. These elementary facts apply
to every degree u and underlie the supplied constructor.

field.md gives every additional algorithmic argument: binary fixed-space
factor extraction with a Chinese-remainder proof; bounded-degree
cyclotomic root lifting and a nonresidue; the prime Kummer extension;
a concrete relative trace whose positive basis coefficients prove it has
degree k; and its binary minimal polynomial. No factoring, nonresidue,
or irreducible-construction oracle is left as an external premise.
Every factor matrix has fewer than k^2 rows and columns. A generous
O((k+1)^12) bit bound covers the constructor, whose represented extension
has fewer than k^2 binary coordinates. It does not enumerate that field.

Parameter selection scans polynomially many prime/order choices with
trial primality and exact integer loops. k,r,t<=N+1, and all threshold
integers have O(log(N+1)) bits. The represented field constructor is
performed only for the selected k, followed by schoolbook field arithmetic
and the polynomial-time symmetric recurrence. A coarse O((N+1)^12) bound
covers the complete evaluator, including all exceptional lengths. No
practical runtime at every large N is asserted. Cache entries memoize
pure functions and are unnecessary to this cold-start polynomial bound.

Let A be this single deterministic evaluator and L={x:A(x)=1}. The NP
verifier accepts exactly empty witnesses w when A(x)=1. The empty witness
exists exactly for members of L, establishing both directions uniformly.
Witnesses are never target-circuit inputs. No field coefficients are
advice or existential witnesses: Q_k is constructed by the stated code.

For retained finite audits, the old paid-prefix irreducibility criterion
is exact at prime degree k: Z^(2^k)=Z modQ makes Q divide a squarefree
polynomial whose irreducible factors have degrees dividing k, and
 gcd(Q,Z^2+Z)=1 removes both linear factors. All remaining factors have
degree k, hence a monic degree-k Q passing both tests is irreducible.
Conversely irreducible Q passes. This test is not needed to choose the
new field; it independently checks finite constructor outputs.

## 3. Prime-degree affine payload degree

For every selected PRIME k, disperser.md proves the required field and
integer threshold statements. Section12 permits t selected blocks of
nonnegative deficits delta_i=d_i-w_i, at least one positive, provided
sum_positive delta_i+max_positive d_i>k+1. Positive sources kill every
unwanted trace shift; zero-deficit sources supply nonzero minors by the
linearized root bound. Triangular mixing, symmetric noncancellation and
Section9 then prove Boolean degree>=B=t(t+3)/2 on each such affine space.

The integer optimization in Sections11--12 is separate: H+E0 is exactly
the largest dimension of an integer vector0<=d_i<=k failing the enlarged
selection predicate. This combinatorial maximum holds for the stated
integers; its FIELD consequence requires prime k. Fewer than t eligible
indices give excess T-r+t-1. With s<t positive indices but enough
eligible indices, all positives can be padded with zeros, giving precisely
a_s via individual upper bounds and their integer sum. With s>=t, replacing
selected zeros by unused positives gives the positive-only top/tail
formula. Every branch has an attaining failed vector. Section13's count
monotonicity computes the identical maximum in fewer steps.

Consequently every affine payload dimension at least d=H+E0+1 has degree
at least B. The guard k>r+1 keeps all exponent weights below k. No
composite-degree coefficient assertion is made; the rejected extra claim
of provider bf5b6832 is not reused. Predicate optimality does not assert
optimal affine dispersers, circuit sizes or limits of other methods.

## 4. Full B2 and the finite actual-length inequality

Fix the p ignored actual input bits in any circuit for F_N. Constant
normalization cannot add gates. The remaining circuit computes the
payload over the deterministically constructed field K. The m actual
payload bits embed by zero padding into rk bits; that map is affine and
injective, preserving every affine dimension. Therefore the payload has
degree>=B on every dimension-d affine space.

The supplied elimination.md Section6 proves that this hypothesis with
1<=B<=d gives C_B2(g)>=4m-8d+min(5B-1,4B+3). Its complete charged
representation covers all sixteen binary truth tables, constants, unary
gates, repeated pins, arbitrary sharing and fanout, pair promises,
closure and terminal cases. No free affine basis or formula restriction
is assumed. Since t>=2 gives B>=5, its offset is4B+3. Hence

    C_B2(F_N)>=4(N-p)-8d+4B+3=4N-L_N.

For the zero option L_N=4N and the same assertion is the zero bound.
All ignored bits are counted in actual N. No input bits specify the field.

## 5. Every sufficiently large integer length

The new search contains every old zero-deficit certificate: its old p=k
and m=N-k belong to the allowed prefix range, and the two-r reduction
only increases m within the same block count and decreases p. The new
threshold also satisfies d<=d_old, where the original sufficient dimension
is H+max(T,floor((rk-r(r+1)/2)/(t+1)))+1. Sections10--12 prove this
inclusion using explicit failed vectors and sufficient subset bounds,
not any source verdict.

Take the least prime k at least the least integer b with b^5>32N^3,
and t=ceil(sqrt(r)) at the old payload m=N-k. Bertrand gives
k=Theta(N^(3/5)),r=Theta(N^(2/5)),t=Theta(N^(1/5)) at every large integer
N. Eventually all guards hold and B<=H<d<=d_old=O(N^(4/5))<N-k.
Thus the new selected loss is at most O(N^(4/5)). It is nonnegative,
since p>=0,B<=H and d>=H+1, including the nonnegative zero option.
Therefore L_N=o(N) over ALL integer lengths. For each epsilon>0 there
is N_epsilon such that C_B2(F_N)>=(4-epsilon)N for every integer
N>=N_epsilon, with one fixed uniform language and verifier throughout.

## 6. Complete reference comparisons

The original, weighted, positive-exact and zero-deficit references use the
same old prime/order loops and unconstrained k-bit prefix, differing only
in their explicitly supplied thresholds. Write their losses L_old,
L_weighted,L_exact,L0, with4n on empty searches. The compressed reference
has L_comp(N)=min(4N,min_(s=0,...,5)(L0(N+s)-4s)), with ties(loss,k,t,5-s).
Its k-5-bit affine decoder and existence proof are retained in family.py
and field.md Section4; the ignored bits make its actual prefix p=k-s.
Every such prefix has0<=p<=k, so the new optimum satisfies L_N<=L_comp(N).
It also dominates each earlier reference whenever present. These are
numeric certified-loss comparisons, not actual optimal circuit sizes.

Complete exact searches give

    N      original  weighted  positive-exact  zero-deficit  compressed  new
    500      1461      1437         1429           1413        1393     1257
    1000     2533      2493         2485           2485        2465     2265
    2000     4513      4449         4433           4417        4397     4061
    10000   16809     16601        16577          16521       16501    15365.

Finite experiments check these and the declared bounded algebraic
predicates. Universal constructor correctness, circuit coverage and
asymptotics rest on the written proofs. No strict improvement at every
length or new leading coefficient is claimed.

## 7. Provenance and claim scope

The field payload and full-B2 degree-sensitive theorem originate in
8eef1799-c27f-46f3-8e61-3c4e0d28d145, with the supplied complete arguments
retained. The compact proof and audit patterns build on2309de11-70f3-4d9d-99be-be46178f7983.
The marked-subset reference is retained from8d0d93b5-9da6-405e-98f3-aef47d9b9c09,
positive-only integer optimization from673f938b-d06e-4423-9514-ed84afe6b13f,
and prime-only zero-deficit extension frombf5b6832-24d7-4ceb-9d64-abd59ca97414.
That provider was rejected on an overbroad extra composite-degree field
claim; this package explicitly does not reuse that assertion. The separate
integer maximum and prime-only proof material remain supplied and credited.
The compressed-prefix reference and its audit are retained from
705a962d-524f-424b-8fe1-28e2f3f1cb37. None of their judgments supplies a premise.

The deterministic field construction is a prime-degree specialization of
Shoup1990, with the needed algorithm and elementary proof supplied in
field.md. It is not claimed as a new irreducible-construction algorithm.
The field constructor, its proof, ignored-prefix optimizer and Section13
count-loop reduction are retained from37b2928a-14a5-4e21-a621-12723a6848ff.
The explicit-prime and finite-cover contributions are credited in Section9. Literature for the retained construction and gate argument
remains BK and Mishra--Santhanam ECCC TR26-260 revision1. Unfamiliar
arguments needed for the theorem are reproduced; ordinary field algebra,
root counting, Euclid, Gaussian elimination and Bertrand are used as stated.


## 8. An elementary explicit prime interval, with proof

For every real x>=2^22 there is a prime p with x<p<=121x/100.
Here is a complete proof, using elementary factorial accounting rather
than an external explicit prime-gap theorem. The classical Chebyshev
factorial combination is used; no novelty of that combination is claimed.

Write L(z)=log(floor(z)!) for z>=0 (natural logarithms), and
psi(z)=sum_(p^j<=z) log p, theta(z)=sum_(p<=z)log p. Put

 T(z)=L(z)-L(z/2)-L(z/3)-L(z/5)+L(z/30),
 A=(14log2+9log3+5log5)/30.

The elementary factorial estimate

 |L(z)-(z log z-z)|<=1+log^+(z)                       (8.1)

holds also0<=z<1 with0log0=0. For z>=1, integrate the increasing
function log t between1 and n=floor z to bound log(n!) between
nlog n-n+1 and nlog n-n+1+log n. The change in zlog z-z between
n and z is between0 and log z. These bounds imply(8.1).
For0<=z<1, L=0 and |zlog z-z|<=1. The coefficients of zlog z
and z cancel in T, so for z>=1,

 |T(z)-Az|<=5(1+log z).                              (8.2)

Legendre's factorial valuation formula, obtained by counting the
multiples of each prime power, gives

 T(z)=sum_(p^j<=z) log p * M(z/p^j),
 M(v)=floor v-floor(v/2)-floor(v/3)-floor(v/5)+floor(v/30).

M is30-periodic: its increment is30-15-10-6+1=0. It is constant
on each interval between consecutive integers, takes only0 and1,
and equals1 throughout1<=v<6. Checking its30 values proves these
finite facts; the declared computation repeats that exact floor check.
It follows that

 psi(z)-psi(z/6)<=T(z)<=psi(z).

Iterating the first inequality until the argument is below1 and using
(8.2), the geometric series, and at most1+log z nonzero terms gives

 Az-5(1+log z)<=psi(z)
       <=(6/5)Az+5(1+log z)^2,                       (8.3)

for z>=1. The term count bound uses log6>1. The geometric sum of
arguments is at most6z/5. No estimate for prime counts is assumed.

We also need theta(z)<4z for z>=1. For completeness, the elementary
primorial induction gives product_(p<=n)p<=4^n. For n=2m the
primes m<p<=2m divide binomial(2m,m)<=4^m; multiply by the induction
hypothesis at m. For n=2m+1, the primes m+1<p<=2m+1 divide
binomial(2m+1,m)<=4^m; multiply by the induction hypothesis at m+1.
The latter binomial bound holds at m=1 and its ratio from m to m+1
is(2m+3)(2m+2)/((m+2)(m+1))<4. Base cases n<=2 are immediate.
Thus theta(z)<=floor(z)log4<4z. Consequently,

 0<=psi(z)-theta(z)
       <=4sqrt(z)+4(log2 z)z^(1/3).                  (8.4)

Indeed separate the exponent j=2, and bound all j>=3 by z^(1/3);
there are at most log2 z exponents, since p>=2.

Let y=121x/100. Equations(8.3)--(8.4) give

 theta(y)-theta(x)>=Ax/100-E(x),
 E(x)=5(1+log y)+5(1+log x)^2
                       +4sqrt(y)+4(log2 y)y^(1/3).

A>9/10. One fully rational way to see this is to use the positive
series log q=2sum_(j>=0) ((q-1)/(q+1))^(2j+1)/(2j+1): its first
three terms for q=2 and3 and first five terms for q=5 imply
(14log2+9log3+5log5)/30>9/10. The declared computation checks the
rational partial sums. This series follows by integrating the geometric
series for2/(1-t^2) from0 to(q-1)/(q+1); only lower bounds are used.

At x0=2^22, y0=(121/100)2^22, we have

 log x0<16, log y0<16, log2 y0<23,
 sqrt(y0)<2253, y0^(1/3)<172,
 E(x0)<85+1445+9012+15824=26366<(9/1000)x0.

For the logarithm comparisons, log2<7/10 follows from the same series
and its geometric tail bound; thus log x0=22log2<15.4. Also
log(121/100)<21/100 by integrating1/t<1 over[1,121/100], so log y0<15.61.
The other bounds follow by squaring, cubing and y0<2^23.
Each summand of E(x)/x decreases for x>=x0: direct differentiation
reduces this to log x>1 and log y>3/2. Both hold since x0>9>e^2;
for example e<3 follows from its factorial series. Therefore
E(x)<(9/1000)x<Ax/100 for all x>=x0. We obtain theta(121x/100)>theta(x),
which is possible only if a prime lies in(x,121x/100]. This proves the
stated prime interval for every real x in the range.


## 9. Every-length remainder577/50=11.54

For EVERY integer N>=0 the same fixed evaluator satisfies

 L_N<=(577/50)N^(4/5),
 C_B2(F_N)>=4N-(577/50)N^(4/5).                       (9.1)

For N>=2^39 choose h=ceil((12023/10000)N^(3/5)), and the least prime k>=h.
Since h>=2^22, Section8 implies k<=121h/100. With
u=k/N^(3/5),

 12023/10000<=u<=1454783/1000000+(121/100)N^(-3/5)<1454784/1000000.

Here N^(3/5)>=2^(117/5)>2^22 and121/(100*2^22)<1/1000000.
Take r=ceil(N/k),t=ceil(sqrt(r)),p=0. All guards follow explicitly:
1<u<2, k<N, r>N^(2/5)/2>2048, and with x=N^(1/5)>=64,

 k> (9/10)x^3,
 r+1< (10/9)x^2+2< (9/10)x^3.

Thus k>r+1, k is an odd prime at least5, and2<=t<r.

The ORIGINAL sufficient dimension is

 d_ref=H+max(T,floor(U/(t+1)))+1,
 H=r(r+3)/2, T=(t-1)k-t(t+1)/2+1,
 U=rk-r(r+1)/2.

Since t=ceil(sqrt(r)), T<=k sqrt(r) and U/(t+1)<=k sqrt(r).
The exact zero-deficit threshold used by the selector is no larger
than this retained sufficient threshold, so

 d_new<=d_ref<=H+k sqrt(r)+1.

Using r<=N/k+1 and sqrt(1+z)<=1+z/2 gives

 H<=N^2/(2k^2)+(5/2)N/k+2,
 k sqrt(r)<=sqrt(Nk)+k^(3/2)/(2sqrt(N)),
 d_ref/N^(4/5)<=f(u)+g(u)delta+3delta^2,
 f(u)=1/(2u^2)+sqrt(u),
 g(u)=5/(2u)+u^(3/2)/2, delta=N^(-2/5).              (9.2)

Both f and g are convex on the relevant interval: f''(u)=
u^(-4)(3-u^(5/2)/4)>0 for u<2, and g''(u)=5/u^3+3/(8sqrt u)>0.
Their endpoint bounds are fully rational. Put a=12023/10000,
u_hi=1454784/1000000, v=1454783/1000000, and

 sa=54824721/50000000, su=120614427/100000000.

Squaring gives sa^2>a and su^2>u_hi. Direct rational arithmetic gives

 1/(2a^2)+sa<1803/1250,
 1/(2v^2)+su<1803/1250,
 25/12+a*sa/2<3,
 50/29+u_hi*su/2<3.

These upper-bound f and g at their two endpoints: a>6/5, u_hi>v>29/20.
The declared computation repeats all comparisons. Since delta<1/32768
for N>=2^39, convexity yields

 d_ref/N^(4/5)<1803/1250+3/32768+3/1073741824<577/400.

Thus d_ref<N and B<=H<d_new<=d_ref. The full-payload candidate is
admissible in the selector and its loss8d_new-4B-3 is below
(577/50)N^(4/5). This proves the infinite tail of(9.1).

For0<=N<=199 the zero option suffices:200^5*199<=577^5<200^5*200,
with N0 separately equality. The remaining200,...,2^39-1 has a
complete finite generated certificate supplied in audit.py. Its73969
rows are produced by exact integer root, bounded sieve and finite
parameter trials; no success of that search over an infinite range is
assumed. A separate checker verifies primality, all guards, original
sufficient dimensions, left>(r-1)k, right<=rk+k, d_ref<min(left,rk),
and consecutive coverage, plus

 (50cost(v))^5<=577^5*v^4

at BOTH endpoints of each row. Here

 cost(x)=8d_ref-4B-3+4max(0,x-rk).

The original sums are checked by exact arithmetic-series identities:
H=(r+1)(r+2)/2-1, T=sum_(j=2)^t(k-j), U=r(2k-r-1)/2.
These equal the defining sums and use different forms from generation;
all earlier bounded literal dimension comparisons remain executed.
Primality is independently tested by division by every prime through
the square root of the largest row degree, supplied by a separate small
sieve. Every composite has a prime divisor at most its square root,
so this is equivalent to full trial division; the generator's large
sieve is not reused by the checker.

For every interior integer N in a row, the actual selector includes
its k,r,t with m=min(N,rk),p=N-m in[0,k]. Since
B<=H<d_new<=d_ref<m, L_N<=cost(N). The function
(577/50)x^(4/5)-cost(x) is concave for x>0: the first term is
concave and negative max(0,x-rk) is concave. Hence nonnegative
endpoints imply nonnegative values throughout. First left200,
adjacency, and final right549755813887=2^39-1 prove complete finite
coverage. Together with the zero range and analytic tail this proves
(9.1) at EVERY integer N. No circuit upper bound is asserted.

The supplied complete prime lemma is classical Chebyshev accounting
with explicit constants. The finite-cover/concavity method builds on
d10eb869; the generated-certificate pattern builds on5ec6176f and
2dc9a80c; the supplied elementary-prime-proof route builds on39c9e977.
Earlier weaker explicit claims/certificates are omitted from this
package to keep its declared computation focused on the stronger bound.
All foundational field, disperser, elimination and comparison proofs
remain supplied; none of the predecessor judgments is a premise.
