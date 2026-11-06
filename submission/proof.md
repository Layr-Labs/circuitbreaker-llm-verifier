# Explicit11.5213 remainder at every integer input length

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

For every real x>=2^47 there is a prime p with x<p<=(11669/10000)x.
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

We strengthen BOTH estimates with a complete elementary simultaneous argument.
Let l=51072/50999, u=59595/50999, L=lA and U=uA. For every real z>=0,

 Lz-32sqrt(z)<=psi(z)<=Uz+80sqrt(z).                  (8.5)

The same30-value table shows M=1 on[7,10) and M=0 on[24,29).
Together with M=1 on[1,6) and0<=M<=1, this gives

 T(z)>=psi(z)-psi(z/6)+psi(z/7)-psi(z/10),
 T(z)<=psi(z)-psi(z/24)+psi(z/29).

The constants obey the exact rational identities

 l=1+l/24-u/29,  u=1+u/6-l/7+u/10.                  (8.6)

Using(8.2) and the appropriate upper/lower induction hypotheses yields

 psi(z)>=Lz-[5(1+log z)+(32/sqrt24+80/sqrt29)sqrt(z)],
 psi(z)<=Uz+[5(1+log z)+(80/sqrt6+32/sqrt7+80/sqrt10)sqrt(z)]

for z>=1. Here1+log z<=(3/2)sqrt z: differentiation shows that its
ratio to sqrt z is maximized at z=e with value2/sqrt e<sqrt2<3/2.
Squaring proves

 1/sqrt24<21/100, 1/sqrt29<19/100,
 1/sqrt6<41/100, 1/sqrt7<38/100, 1/sqrt10<32/100.

Thus the first error coefficient is below7.5+32*.21+80*.19=29.42<32;
the second is below7.5+80*.41+32*.38+80*.32=78.06<80. These prove
both assertions inductively. For0<=z<1, psi=0 and the inequalities
hold because0<L<32 and U>0. The crude A<5 follows from log q<q-1,
and l<2. Every recurrence argument is at most z/6, so its tree is finite.
This establishes(8.5) for every real argument. No prime-gap theorem is used.

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

Let q=11669/10000 and x>=x0=2^47, and put y=qx. Since
U/L=u/l=19865/17024, we obtain

 theta(y)-theta(x)>=L(q-19865/17024)x-E(x),
 E(x)=80sqrt(x)+36sqrt(y)+4(log2 y)y^(1/3).

Here q-19865/17024=191/10640000. Also L>A>9/10. To verify the latter
bound on A without any numerical logarithm oracle, use the positive
series log v=2sum_(j>=0)((v-1)/(v+1))^(2j+1)/(2j+1). Its first
three terms for v=2 and3, and first five terms for v=5, give
(14log2+9log3+5log5)/30>9/10. The series follows by integrating
2/(1-t^2); the declared computation checks these rational partial sums.

At x0 the following integer/rational bounds hold:

 sqrt(x0)<11863284, sqrt(qx0)<12815092,
 (qx0)^(1/3)<54763, log2(qx0)<48,
 E(x0)<80*11863284+36*12815092+192*54763=1420920528,
 1420920528<(1719/106400000)2^47.

All power comparisons and the final cross-multiplication are repeated by
the computation. Since q<2, the logarithm bound follows from qx0<2^48.
Each term of E(x)/x decreases for x>=x0: the square-root terms are
constant multiples of x^(-1/2), and the last is a positive constant
times log(qx)x^(-2/3), decreasing when log(qx)>3/2. That condition
holds since x0>9>e^2, with e<3 from its factorial series. Hence

 E(x)<(1719/106400000)x<L(q-19865/17024)x

for every x>=x0. Therefore theta(qx)>theta(x), proving a prime in
(x,qx] for EVERY real x>=2^47.

## 9. Universal remainder115213/10000=11.5213

For EVERY integer N>=0, the SAME fixed evaluator satisfies

 L_N<=(115213/10000)N^(4/5),
 C_B2(F_N)>=4N-(115213/10000)N^(4/5).                 (9.1)

We first give common sufficient-dimension estimates. For any selected
prime k with1<k/N^(3/5)<2, put r=ceil(N/k),t=ceil(sqrt r),p=0.
For N>=2^30, all guards hold: k<N,r>N^(2/5)/2>2048, and with
x=N^(1/5)>=64,

 k>(9/10)x^3,
 r+1<(10/9)x^2+2<(9/10)x^3.

Thus k>r+1, k is odd and at least5, and2<=t<r. The original sufficient
dimension is d_ref=H+max(T,floor(U0/(t+1)))+1, where

 H=r(r+3)/2,T=(t-1)k-t(t+1)/2+1,U0=rk-r(r+1)/2.

Since t=ceil(sqrt r), both T and U0/(t+1) are at most k sqrt r.
The exact zero-deficit selector threshold is no larger, so

 d_new<=d_ref<=H+k sqrt r+1.

Using r<=N/k+1 and sqrt(1+z)<=1+z/2 gives

 H<=N^2/(2k^2)+(5/2)N/k+2,
 k sqrt r<=sqrt(Nk)+k^(3/2)/(2sqrt N),
 d_ref/N^(4/5)<=f(v)+g(v)delta+3delta^2,              (9.2)
 f(v)=1/(2v^2)+sqrt v,g(v)=5/(2v)+v^(3/2)/2,
 v=k/N^(3/5),delta=N^(-2/5).

Both f and g are convex on(1,2): f''(v)=v^(-4)(3-v^(5/2)/4)>0,
and g''(v)=5/v^3+3/(8sqrt v)>0. Thus endpoint upper bounds suffice.

For N>=2^79 take a=30583/25000=1.22332,h=ceil(aN^(3/5)), and the
least prime k>=h. Since h>=2^47, Section8 gives k<=qh with
q=11669/10000. Hence

 a<=v<=qa+qN^(-3/5)<1427492109/1000000000=:b.

Here b=qa+10^(-9), and q/2^47<10^(-9). Put
sa=11060379741/10000000000,sb=11947770123/10000000000.
Squaring gives sa^2>a,sb^2>b. Exact rational comparisons give

 1/(2a^2)+sa<28803/20000,
 1/(2(qa)^2)+sb<28803/20000,
 25/12+a*sa/2<3,25/14+b*sb/2<3.

These bound f at both endpoints by1.44015 and g by3, since a>6/5
and b>qa>7/5. As delta<2^(-31), (9.2) yields

 d_ref/N^(4/5)<28803/20000+3/2^31+3/2^62<115213/80000.

Thus d_ref<N,B<=H<d_new<=d_ref, the candidate is admissible, and
8d_new-4B-3<(115213/10000)N^(4/5). This proves the infinite tail.

For the finite but large middle range2^30<=N<2^79, use instead
h=ceil((33/25)N^(3/5)). It lies in[2^18,2^48]: indeed
N^(3/5)<2^(237/5)=2^47*2^(2/5), with2^(2/5)<3/2 by fifth powers;
therefore h<=(99/50)2^47+1<2^48. Section10's finite prime-chain
certificate proves there is a prime at least h and at most(101/100)h
for EVERY real argument h in that full interval. The least prime k>=h
therefore gives

 33/25<=v<(133321/100000)=:b0.

The rounding allowance follows from101/(100*2^18)<1/100000.
Take s0=1148913/1000000 and s1=115465/100000; their squares exceed
the two endpoints. Exact endpoint comparisons establish

 1/(2(33/25)^2)+s0<359/250,
 1/(2(13332/10000)^2)+s1<359/250,
 25/12+(33/25)*s0/2<3,
 25/13+b0*s1/2<3.

For the upper endpoint reciprocals we used b0>13332/10000>13/10.
Thus f<1.436 and g<3 throughout this interval. Since delta<=1/4096,

 d_ref/N^(4/5)<359/250+3/4096+3/2^24<115213/80000.

The same guard/admissibility argument proves(9.1) on this ENTIRE middle
range, not just at sampled N. All displayed rational comparisons are
reproduced by the declared computation.

For0<=N<=198 the zero option suffices:40000^5*198<=115213^5<40000^5*199,
with N0 separately equality. The remaining199,...,2^30-1 has a complete
finite generated certificate in audit.py:6484 adjacent rows, first left199,
last right1073741823. A separate checker verifies primality, all guards,
original sufficient dimension, left>(r-1)k,right<=rk+k,d_ref<min(left,rk),
consecutive coverage and, at both endpoints v,

 (10000cost(v))^5<=115213^5*v^4,
 cost(x)=8d_ref-4B-3+4max(0,x-rk).

It computes the original sums by exact arithmetic-series identities
H=(r+1)(r+2)/2-1,T=sum_(j=2)^t(k-j),U0=r(2k-r-1)/2. The older
bounded literal dimension checks are also retained. Every row includes
an actual selector candidate with m=min(N,rk),p=N-m in[0,k] and
B<=H<d_new<=d_ref<m. Thus L_N<=cost(N). Concavity of
(115213/10000)x^(4/5)-cost(x) transfers endpoint inequalities to
EVERY interior integer. Together with the zero, chain-based middle and
analytic infinite ranges this proves(9.1) at every integer N.

## 10. A finite prime chain with independently checked primality witnesses

The supplied computation generates a finite list p0<...<pM with

 2^18<=p0<=(101/100)2^18,
 p_(i+1)<=(101/100)p_i,
 pM>=2^48,

and verifies that every entry is prime. Therefore EVERY real h in
[2^18,2^48] has a prime in[h,(101/100)h]: use p0 when h<=p0, and
the first chain entry at least h otherwise. Its predecessor is below h,
so its ratio bound implies the result. This is a finite coverage proof,
not an extrapolation of prime gaps beyond the last entry.

The chain contains2596 entries:1045 small primes checked by independent
trial division, and1551 entries with Proth witnesses. A small-prime entry
is below2^32 and is checked by division by every prime through its square
root using a separate small sieve. Every composite has such a prime factor.
A Proth entry supplies integers(m,c,b) with

 p=c2^m+1, m>=1,0<c<2^m,c odd,
 1<b<p, b^((p-1)/2)=-1 modulo p.

Here is the complete primality proof for this criterion. For any prime
ell dividing p, b is invertible modulo ell. Its order divides p-1=c2^m,
since squaring the displayed congruence gives b^(p-1)=1. Its order does
not divide c2^(m-1), because -1 is not1 modulo odd ell. Thus its order
has2-part exactly2^m. Fermat's theorem implies2^m divides ell-1,
so ell>=2^m+1>sqrt p, using c<2^m. A composite integer p has a prime
factor at most sqrt p, a contradiction. Hence p is prime. This is the
classical Proth criterion, with the necessary argument reproduced.

The generator uses integer candidate trials, Jacobi symbols to choose
potential witnesses, and modular exponentiation. Its Jacobi heuristic
is not a premise: the checker verifies the actual power congruence for
every emitted witness, together with the form, inequalities and coverage.
If finite generation did not succeed, the declared check would fail.
Only this fixed finite range is asserted, and the actual computation
reproduces and checks the complete chain from supplied source, without
external data, model credentials, network or a prime table. The chain
and interval generator are evidence, not part of the target evaluator.

The finite-cover/concavity method builds on d10eb869, generated certificates
on5ec6176f and2dc9a80c, elementary prime proofs on39c9e977, balanced
certificate code on049df349, and the strengthened one-sided psi recurrence
on7ae3664d. The coupled leading constants are classical Sylvester1881 values,
reproduced in Gantumur, arXiv:2512.02466v2, Section6, Equation86
(https://arxiv.org/abs/2512.02466v2). The needed argument is fully
supplied above; that citation is provenance, not a missing premise.
The additions here are the explicit square-root errors and all-real
prime threshold, finite Proth-witness chain coverage, and application
to the universal11.5213 circuit remainder. All needed
proofs are supplied; no provider judgment or promotion is a premise.
