# Explicit11.4891 remainder at every integer input length

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


## 8. A verified finite floor scheme and an explicit prime interval

We prove that EVERY real x>=2^220 has a prime in(x,(109/100)x].
All unfamiliar ingredients, including the finite step certificates, are
supplied and checked; no external explicit prime-gap theorem is a premise.

Let F(z)=log(floor(z)!) for z>=0, psi(z)=sum_(p^j<=z)log p and
theta(z)=sum_(p<=z)log p, with natural logarithms. Use the classical
Sylvester scheme

 Dplus=(1,6,10,14,105),
 Dminus=(2,3,5,7,11,13,385,1001),
 V(z)=sum_(d in Dplus)F(z/d)-sum_(d in Dminus)F(z/d),
 M(v)=sum_(d in Dplus)floor(v/d)-sum_(d in Dminus)floor(v/d).

This is Sylvester1892's scheme nu8, reproduced in Gantumur,
arXiv:2512.02466v2, Equation90. That citation identifies the classical
scheme; every needed argument and finite certificate is supplied below.
Neither the scheme nor general Sylvester refinement is claimed as invented
here. The selected finite step certificates are not claimed optimal.

The reciprocal sums of Dplus and Dminus are equal, so M is30030-periodic:
every divisor divides30030 and its increment is30030 times that zero
reciprocal difference. M is constant on each[j,j+1). Its30030 integer
values have minimum-1 and maximum4, checked exactly in the computation.

Define A=sum_(d in Dminus)(log d)/d-sum_(d in Dplus)(log d)/d. Factoring
the divisors gives the exact expression

 A=(17/105)log2+(11/70)log3+(43/462)log5
       +(281/4290)log7+(43/455)log11+(6/77)log13.      (8.1)

All coefficients are positive. The positive series
log v=2sum_(i>=0)((v-1)/(v+1))^(2i+1)/(2i+1), obtained by integrating
2/(1-t^2), with20 terms for each of the six primes proves A>49/50.
The computation verifies the exact rational truncations and factorization.
Also A<5 follows from log v<v-1 and the rational coefficients in(8.1).

The elementary factorial estimate

 |F(z)-(z log z-z)|<=1+log^+(z)                       (8.2)

holds, with0log0=0. For z>=1 put n=floor z and integrate the increasing
logarithm to bound log(n!) between nlog n-n+1 and nlog n-n+1+log n.
The change in zlog z-z from n to z is between0 and log z. For0<=z<1,
F=0 and |zlog z-z|<=1. This proves(8.2). Cancellation of the zlog z
and z terms by the reciprocal identity yields, for z>=1,

 |V(z)-Az|<=13(1+log z).                             (8.3)

Counting multiples of every prime power in each factorial gives

 V(z)=sum_(p^j<=z)(log p) M(z/p^j).                  (8.4)

Two explicit integer step lists LOWER_STEP and UPPER_STEP are supplied
in audit.py. They have59 and57 nonzero coefficients, respectively; both
start with(1,1). Write their coefficients as P_j and Q_j, and set
P(v)=sum_j P_j 1_(v>=j), Q(v)=sum_j Q_j 1_(v>=j). The supplied finite
check proves

 Q(v)<=M(v)<=P(v) for EVERY real v>=0.               (8.5)

Indeed it checks all30030 integer cells, including0; the functions are
constant on their unit cells. Beyond30030 all list jumps are finished,
P is4 and Q is-1, and periodic M remains in[-1,4]. This proves the
universal extension. The finite arrays and complete verification code
are in the eight-file archive, not an external data attachment.

Since log p>=0, (8.4)--(8.5) imply

 psi(z)+sum_(j>=2)Q_j psi(z/j)<=V(z)
                <=psi(z)+sum_(j>=2)P_j psi(z/j).     (8.6)

Let pP=sum_(j>=2,P_j>0)P_j/j, nP=sum_(j>=2,P_j<0)P_j/j, and define
pQ,nQ similarly. Put c=96863/100000,d=105557/100000,L=cA,U=dA. Direct
rational arithmetic on the supplied lists verifies

 (1+nP)c+pP*d<=1,
 pQ*c+(1+nQ)d>=1.                                   (8.7)

We now prove for EVERY real z>=0 that

 Lz-400z^(9/10)<=psi(z)<=Uz+400z^(9/10).             (8.8)

For each list and j>=2 let r_j=floor((j^9)^(1/10)). Its positive integer
value satisfies r_j^10<=j^9, so j^(-9/10)<=1/r_j. The exact finite
check establishes for BOTH lists

 sum_(j>=2)|P_j|/r_j<189/200,
 sum_(j>=2)|Q_j|/r_j<189/200.                        (8.9)

Thus each corresponding sum with j^(-9/10) is below189/200.
For z>=1, 1+log z<=(3/2)sqrt z: differentiation gives maximum ratio
2/sqrt e<sqrt2<3/2, since e>2. Also sqrt z<=z^(9/10). Apply the
induction hypotheses in(8.6), taking the correct upper or lower bound
for each coefficient sign. Equations(8.3),(8.7) bound the lower linear
term by at least Lz and the upper by at most Uz. Either error is below

 [13*(3/2)+400*(189/200)]z^(9/10)
                  =397.5z^(9/10)<400z^(9/10).

This proves both induction steps. For0<=z<1, psi=0 and the bounds hold
because0<L<5<400,U>0 and z^(9/10)>=z. All nontrivial list indices
are at least15, so the recurrence tree has finite depth, with arguments
shrinking by at least15 at every level. This completes the simultaneous
proof of(8.8); no convergence or uncomputed infinite iteration is assumed.

For completeness, theta(w)<4w for w>=1 follows by a primorial induction.
The primes m<p<=2m divide binomial(2m,m)<=4^m. The primes m+1<p<=2m+1
divide binomial(2m+1,m)<=4^m. The latter bound holds at m=1 and its
successive ratio is(2m+3)(2m+2)/((m+2)(m+1))<4. Induction yields
product_(p<=n)p<=4^n for every integer n, with n<=2 as base cases.
Therefore theta(w)<=floor(w)log4<4w. Separating exponent2 and bounding
all higher exponents by w^(1/3) gives

 0<=psi(w)-theta(w)<=4sqrt w+4(log2 w)w^(1/3).        (8.10)

Take q=109/100,y=qx. Then U/L=d/c=105557/96863 and
q-d/c=2367/9686300. Since L>9/10 by(8.1), equations(8.8)--(8.10) give

 theta(y)-theta(x)>=L(q-d/c)x-E(x),
 E(x)=400(1+q)x^(9/10)+4sqrt(qx)+4(log2(qx))(qx)^(1/3).

Here q^(9/10)<q was used. At x0=2^220, qx0<2^221, so

 E(x0)<400(1+q)2^198+4*2^111+884*2^74
                   <(21303/96863000)2^220.

The last strict inequality is an exact rational check in the computation;
21303/96863000=(9/10)(q-d/c). Each term of E(x)/x decreases for
x>=x0: the first two are constant multiples of x^(-1/10) and x^(-1/2),
and the last decreases when log(qx)>3/2. This holds because x0>9>e^2;
e<3 follows from its factorial series. Hence theta(qx)>theta(x) for
every x>=x0, proving the stated prime interval on its entire real range.

### 8.1 Two rational kernels with simultaneous error control

We improve the prime interval: EVERY real x>=2^270 has a prime
in(x,(10387/10000)x]. The complete rational data are RATIONAL_SCHEMES in
audit.py. The lower certificate has174kernel coefficients and436step
jumps with common denominator36036000000; the upper has162kernel coefficients
and430step jumps with denominator630000000. A pair(d,v) in its coefficient
list means c_d=v/denominator; a pair(j,v) in its step list means
b_j=v/denominator. Both lists begin(1,denominator).

For each certificate define M(v)=sum c_d floor(v/d),
V(z)=sum c_d F(z/d), A=-sum c_d(log d)/d,
and B(v)=sum b_j 1_(v>=j). Every d divides360360 and sum c_d/d=0.
The exact finite checker verifies both identities. Hence M is360360-periodic.
It verifies ALL cells0,...,360359 and the final jump at360360, showing
M<=B for the lower certificate and B<=M for the upper. The lower
constant tail is at least max M, and the upper tail at most min M.
The checker computes every integer cell exactly by a difference sieve. For
integer j>=1, floor(j/d)-floor((j-1)/d) is1 precisely when d divides j
and0 otherwise. Adding each numerator at all multiples of its divisor, then
prefix summing, therefore gives exactly the numerator of M(j), with M(0)=0.
This checks every cell and omits no point or predicate.
Cell constancy and periodicity therefore extend these inequalities to
EVERY real v>=0. Both M and B vanish on[0,1). The lists have distinct
increasing indices; every nonleading index j is at least2.

The reproduced factorial estimate and prime-power valuation argument give

 |V(z)-Az|<=K(1+log z), K=sum|c_d|, z>=1,
 V(z)=sum_(p^h<=z)(log p)M(z/p^h).

Write P_j for the lower certificate's b_j, Q_j for the upper's.
Their verified envelopes imply

 psi(z)+sum_(j>=2)Q_jpsi(z/j)<=V_upper(z),
 V_lower(z)<=psi(z)+sum_(j>=2)P_jpsi(z/j).            (8.11)

For each list let p=sum_(j>=2,b_j>0)b_j/j and
n=sum_(j>=2,b_j<0)b_j/j. Put

 L=981285/1000000,U=1019145/1000000.

The exact checker proves

 (1+nP)L+pP*U < A_lower,
 pQ*L+(1+nQ)U > A_upper.                            (8.12)

Here A_lower and A_upper refer to the two different kernels. To check
these logarithmic statements, factor every d into2,3,5,7,11,13 and
combine rational coefficients of their logs. For each prime p put
u=(p-1)/(p+1). Its log is enclosed by the first60terms of
2sum_(i>=0)u^(2i+1)/(2i+1), with positive remaining error at most
2u^121/(121(1-u^2)). This follows from the geometric sum after replacing
each remaining denominator by121. The endpoints are reversed for
negative rational coefficients. Thus(8.12) consists of finite exact
rational inequalities, not floating-point conclusions.

For BOTH lists the checker further proves

 (3/2)K+2000sum_(j>=2)|b_j|/floor((j^9)^(1/10))<2000. (8.13)

Each integer root is independently bounded by consecutive tenth powers.
Together with1+log z<=(3/2)z^(9/10), this proves for EVERY real z>=0

 Lz-2000z^(9/10)<=psi(z)<=Uz+2000z^(9/10).           (8.14)

For completeness this is a finite simultaneous induction, not an assumed
fixed-point limit. At0<=z<1, psi(z)=0 and both bounds hold because
0<L<U<2 and2000z^(9/10)>=2z. At z>=1 all arguments z/j in(8.11)
are at most z/2, so repeated substitutions reach that base interval in
finitely many levels. For the lower bound subtract the sum involving P:
use upper estimates for P_j>0 and lower estimates for P_j<0. Its leading
term is at least(A_lower-pP*U-nP*L)z>=Lz by(8.12). For the upper
bound subtract the sum involving Q, using lower estimates when Q_j>0
and upper estimates when Q_j<0. Its leading term is at most
(A_upper-pQ*L-nQ*U)z<=Uz. The absolute error in either case is at most
the left side of(8.13) times z^(9/10), strictly below2000z^(9/10).
This proves the induction and handles every real argument, with all
coefficient signs and finite recursion depths explicit.

Put q=10387/10000. Exact arithmetic verifies Lq-U>0 and

 2000(1+q)2^243+4*2^136+1084*2^91 < (Lq-U)2^270.

At x0=2^270, qx0<2^271, so this bounds the complete higher-prime-power
and error contribution from(8.10). The corresponding error divided by x
decreases thereafter, as proved above. Thus theta(qx)>theta(x) for
EVERY real x>=2^270, proving the improved prime interval. The rational
kernels were discovered by a heuristic linear program, but every literal
coefficient, envelope and arithmetic obligation used here is checked
exactly. No optimizer, numerical optimum or prior verdict is a premise.

## 9. The11.4891 remainder at EVERY integer input length

For the SAME fixed evaluator and EVERY integer N>=0,

 L_N<=(114891/10000)N^(4/5),
 C_B2(F_N)>=4N-(114891/10000)N^(4/5).                        (9.1)

We first give common sufficient-dimension estimates. For a prime k with
1<k/N^(3/5)<2, take r=ceil(N/k),t=ceil(sqrt r),p=0. For N>=2^30,
k<N,r>N^(2/5)/2>=2048, and with x=N^(1/5)>=64,

 k>(9/10)x^3,
 r+1<(10/9)x^2+2<(9/10)x^3.

Thus k>r+1, k is odd and at least5, and2<=t<r. Put
H=r(r+3)/2,T=(t-1)k-t(t+1)/2+1,U0=rk-r(r+1)/2. The original
sufficient dimension d_ref=H+max(T,floor(U0/(t+1)))+1 satisfies

 d_new<=d_ref<=H+k sqrt r+1,
 H<=N^2/(2k^2)+(5/2)N/k+2,
 k sqrt r<=sqrt(Nk)+k^(3/2)/(2sqrt N).

Here d_new is the actual exact zero-deficit threshold. The first bound
uses t=ceil(sqrt r); the next two use r<=N/k+1 and sqrt(1+z)<=1+z/2.
Writing v=k/N^(3/5),delta=N^(-2/5), we get

 d_ref/N^(4/5)<=f(v)+g(v)delta+3delta^2,              (9.2)
 f(v)=1/(2v^2)+sqrt v,g(v)=5/(2v)+v^(3/2)/2.

Both are convex on(1,2): f''(v)=v^(-4)(3-v^(5/2)/4)>0 and
g''(v)=5/v^3+3/(8sqrt v)>0. Hence endpoint bounds suffice.

For N>=2^450, put a=129481/100000,h=ceil(aN^(3/5)), and let k be the
least prime at least h. Since h>=2^270, Section8.1 gives k<=qh,
q=10387/10000. Thus

 a<=v<1344920147/1000000000=:b,
 b=qa+1/1000000,

using q/2^270<1/1000000. Define
sa=(floor(10^10 sqrt(a))+1)/10^10,sb=(floor(10^10 sqrt(b))+1)/10^10. Squaring and exact
rational arithmetic verify

 sa^2>a,sb^2>b,
 1/(2a^2)+sa<1436133/1000000,
 1/(2(qa)^2)+sb<1436133/1000000,
 25/12+a*sa/2<3,25/13+b*sb/2<3.

The reciprocal g estimates use a>6/5 and b>qa>13/10. Therefore
f<1.436133,g<3 throughout[a,b]. Since delta<=2^(-180), (9.2) yields

 d_ref/N^(4/5)<1436133/1000000+3/2^180+3/2^360<114891/80000.

Thus d_ref<N,B<=H<d_new<=d_ref and the candidate is admissible.
Its loss8d_new-4B-3 is below(114891/10000)N^(4/5). This proves the infinite tail.

For2^35<=N<2^45, use h=ceil(a1*N^(3/5)), a1=65327/50000,
and the least prime k>=h. Here h belongs to[2^21,2^28], since
N^(3/5)>=2^21 and h<=a1*2^27+1<2^28. Section10 also supplies a
separate short chain proving a prime in[h,(51/50)h] for EVERY real h
in that interval. Consequently

 a1<=v<b1=(51/50)a1+1/1000000,

because(51/50)/2^21<1/1000000. Put
s1=(floor(10^10 sqrt(a1))+1)/10^10,
s2=(floor(10^10 sqrt(b1))+1)/10^10. Exact rational checking gives

 s1^2>a1,s2^2>b1,
 1/(2a1^2)+s1<287189/200000,
 1/(2((51/50)a1)^2)+s2<287189/200000,
 25/12+a1*s1/2<3,25/13+b1*s2/2<3.

Thus convexity gives f<1.435945,g<3. With delta<=1/16384,

 d_ref/N^(4/5)<287189/200000+3/16384+3/2^28<114891/80000.

The common guards, admissibility and selector comparison prove(9.1)
throughout this entire range. This uses a complete finite prime-chain
cover, not a sample of values of N.

For2^45<=N<2^450, instead set h=ceil((6501/5000)N^(3/5)). It lies in
[2^18,2^271], because N^(3/5)<2^270 and
h<=(6501/5000)2^270+1<2^271. Section10's finite prime chain proves
for EVERY real h in that interval a prime between h and(103/100)h.
Thus the least prime k>=h gives

 a0=6501/5000<=v<83701/62500=:b0.

The rounding allowance uses(103/100)/2^18<1/100000. Set
s0=(floor(10^10 sqrt(a0))+1)/10^10,s1=(floor(10^10 sqrt(b0))+1)/10^10. Exact comparisons give

 s0^2>a0,s1^2>b0,
 1/(2a0^2)+s0<35901/25000,
 1/(2((103/100)a0)^2)+s1<35901/25000,
 25/12+a0*s0/2<3,25/13+b0*s1/2<3.

Consequently f<1.43604,g<3. Since delta<=1/262144,

 d_ref/N^(4/5)<35901/25000+3/262144+3/2^36<114891/80000.

The same admissibility and selector comparison proves(9.1) throughout
this ENTIRE middle range. Every displayed rational bound is repeated
by the declared computation.

For0<=N<=195 the zero option suffices:40000^5*195<=114891^5<40000^5*196,
with N0 separately equality. The remaining196,...,2^35-1 has a complete
generated25041-row finite certificate in audit.py. A separate checker
verifies primality, guards, original sufficient dimensions,
left>(r-1)k,right<=rk+k,d_ref<min(left,rk), consecutive coverage,
and at both endpoints v,

 (10000cost(v))^5<=114891^5*v^4,
 cost(x)=8d_ref-4B-3+4max(0,x-rk).

The defining sums are computed by exact arithmetic-series identities
H=(r+1)(r+2)/2-1,T=sum_(j=2)^t(k-j),U0=r(2k-r-1)/2; older bounded
literal threshold checks are retained. Each interior N has the actual
selector candidate with m=min(N,rk),p=N-m in[0,k] and
B<=H<d_new<=d_ref<m, so L_N<=cost(N). Concavity of
(114891/10000)x^(4/5)-cost(x) transfers endpoints to every interior integer.
First left196, adjacency, and last right34359738367 exhaust this finite
range. Together with the zero, chain-based middle and analytic tail,
this proves(9.1) at EVERY integer N. No circuit upper bound is asserted.

## 10. A finite prime chain with independently checked primality witnesses

The supplied computation generates a finite list p0<...<pM with

 2^18<=p0<=(103/100)2^18,
 p_(i+1)<=(103/100)p_i,
 pM>=2^271,

and verifies that every entry is prime. Therefore EVERY real h in
[2^18,2^271] has a prime in[h,(103/100)h]: use p0 when h<=p0, and
the first chain entry at least h otherwise. Its predecessor is below h,
so its ratio bound implies the result. This is a finite coverage proof,
not an extrapolation of prime gaps beyond the last entry.

The chain contains6285 entries:421 small primes checked by independent
trial division, and5864 entries with Proth witnesses. A small-prime entry
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

A separate short list q0<...<qT covers[2^21,2^28] with ratio51/50.
It has259 primes, first2097169 and last269604541. The generator uses
trial division; a separate sieve of primes through the square root of
the largest entry independently checks every proposed prime by complete
division. Every composite has a prime divisor at most its square root.
The checker also verifies2^21<=q0<=(51/50)2^21,
q_(i+1)<=(51/50)q_i, and qT>=2^28. The same first-entry/predecessor
argument therefore proves the short-chain coverage for EVERY real h
in its entire interval. Its canonical row digest is
b1d8c4aa61b0187db2baaf5f189a1ec4e8cff4733bb1f69a876262bce7cf28c6.
It supplies only this bounded interval; no extrapolation is asserted.

The finite-cover/concavity method builds on d10eb869; generated certificates
on 5ec6176f and 2dc9a80c; the elementary factorial/prime-proof route on
original 39c9e977, also delivered byte-identically as published 14c284dd;
balanced certificate code on 049df349; finite-tree errors on 7ae3664d;
and coupled estimates/Proth-chain coverage on 6927dad0.

The retained nu8 step lists and count-oracle cache come from 23856b84.
The nu8 rational parameters and proof framework come from 7099a2f9.
The rational-envelope/logarithm/error checker builds on 240f549a, with
the direct simultaneous induction developed in 8131ff33. Source
af5ec9e9-e61e-4b69-b397-1f6db657991c supplied the preceding period60060
certificates and envelope/length proof framework. The current period360360
lists replace those earlier lists; the language remains identical.

This package supplies the current exact rational kernels, full-period
difference-sieve checks, simultaneous error bounds and short-chain bridge
used for the11.4891 remainder. Adaptive discovery constraints are not proof
premises: every period cell and tail is checked. No selection optimality,
historical invention of Sylvester's method, prior verdict or promotion is
assumed. All necessary arguments and finite generators are supplied.
