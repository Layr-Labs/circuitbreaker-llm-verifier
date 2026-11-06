# Nearest-rounding certificate and a positive second correction

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
Appendix A.1--A.3. Use K=F2[Z]/Q_k in its polynomial basis. Pad the
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

Appendix A gives every additional algorithmic argument: binary fixed-space
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
and Appendix A.4; the ignored bits make its actual prefix p=k-s.
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
Appendix A. It is not claimed as a new irreducible-construction algorithm.
The field constructor, its proof, ignored-prefix optimizer and Section13
count-loop reduction are retained from37b2928a-14a5-4e21-a621-12723a6848ff.
The explicit-prime and finite-cover contributions are credited in Section10. Literature for the retained construction and gate argument
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
For real w>=2, theta(w)<=floor(w)log4<4w. Separating exponent2 and bounding
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

We improve the prime interval: EVERY real x>=2^240 has a prime
in(x,(1028118367/1000000000)x]. The complete rational data are RATIONAL_SCHEMES in
certificates.py, loaded and checked by audit.py. The lower certificate has256kernel coefficients and682step
jumps with common denominator1531530000000; the upper has259kernel coefficients
and1056step jumps with denominator1531530000000. A pair(d,v) in its coefficient
list means c_d=v/denominator; a pair(j,v) in its step list means
b_j=v/denominator. Both lists begin(1,denominator).

For each certificate define M(v)=sum c_d floor(v/d),
V(z)=sum c_d F(z/d), A=-sum c_d(log d)/d,
and B(v)=sum b_j 1_(v>=j). Each scheme supplies its own period D:
D=3063060 for the lower list and D=3063060 for the upper list.
Every d in that list divides its D, and sum c_d/d=0. The exact checker
verifies these identities, making each M periodic with its own D.
It verifies ALL cells0,...,D-1 and the final jump at D for EACH list, showing
M<=B for the lower certificate and B<=M for the upper. The lower
constant tail is at least max M, and the upper tail at most min M.
The checker computes EVERY integer cell exactly using blocks of at most65536.
For j>=1, floor(j/d)-floor((j-1)/d) is1 precisely when d divides j and0
otherwise. Within block[left,right), the first multiple is
max(d,d*ceil(left/d)); adding each numerator at all such multiples gives
exactly these differences. The max with d excludes j=0, where M(0)=0.
The accumulated M and B values carry across every block boundary. The
blocks partition0,...,D-1, so induction on blocks and then cells proves
that every computed M(j) is exactly sum v_d floor(j/d), in numerator units.
No cell or predicate is sampled or omitted. At each positive block start,
M(left-1) is independently compared with that direct floor sum; M(right-1)
is compared after EVERY block. There are93direct boundary checks per list.
Extrema include every cell; the final jump and constant tail are checked
as before. This uses bounded memory without weakening any envelope check.
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

 L=986117/1000000,U=1013811/1000000.

The exact checker proves

 (1+nP)L+pP*U < A_lower,
 pQ*L+(1+nQ)U > A_upper.                            (8.12)

Here A_lower and A_upper refer to the two different kernels. To check
these logarithmic statements, factor every d into2,3,5,7,11,13,17 and
combine rational coefficients of their logs. For each prime p put
u=(p-1)/(p+1). Its log is enclosed by the first60terms of
2sum_(i>=0)u^(2i+1)/(2i+1), with positive remaining error at most
2u^121/(121(1-u^2)). This follows from the geometric sum after replacing
each remaining denominator by121. The endpoints are reversed for
negative rational coefficients. Thus(8.12) consists of finite exact
rational inequalities, not floating-point conclusions.

For BOTH lists the checker further proves

 (3/2)K+18000sum_(j>=2)|b_j|/floor((j^7)^(1/8))<18000. (8.13)

Each integer root is independently bounded by consecutive eighth powers.
For z>=1 the function(3/2)z^(7/8)-1-log z is positive at1 and
has derivative[(21/16)z^(7/8)-1]/z>0. Thus
1+log z<=(3/2)z^(7/8). These bounds prove for EVERY real z>=0

 Lz-18000z^(7/8)<=psi(z)<=Uz+18000z^(7/8).           (8.14)

For completeness this is a finite simultaneous induction, not an assumed
fixed-point limit. At0<=z<1, psi(z)=0 and both bounds hold because
0<L<U<2 and18000z^(7/8)>=2z. At z>=1 all arguments z/j in(8.11)
are at most z/2, so repeated substitutions reach that base interval in
finitely many levels. For the lower bound subtract the sum involving P:
use upper estimates for P_j>0 and lower estimates for P_j<0. Its leading
term is at least(A_lower-pP*U-nP*L)z>=Lz by(8.12). For the upper
bound subtract the sum involving Q, using lower estimates when Q_j>0
and upper estimates when Q_j<0. Its leading term is at most
(A_upper-pQ*L-nQ*U)z<=Uz. The absolute error in either case is at most
the left side of(8.13) times z^(7/8), strictly below18000z^(7/8).
This proves the induction and handles every real argument, with all
coefficient signs and finite recursion depths explicit.

Put q=1028118367/1000000000. Exact arithmetic verifies Lq-U>0 and

 18000(1+q)2^210+4*2^121+964*2^81 < (Lq-U)2^240.

At x0=2^240, qx0<2^241, so this bounds the complete higher-prime-power
and error contribution from(8.10). The first normalized error is proportional to x^(-1/8); the other
normalized terms decrease by the derivative calculation following(8.10). Thus theta(qx)>theta(x) for
EVERY real x>=2^240, proving the improved prime interval. The rational
kernels were discovered by a heuristic linear program, but every literal
coefficient, envelope and arithmetic obligation used here is checked
exactly. No optimizer, numerical optimum or prior verdict is a premise.

## 9. Complete quantitative bound with an adaptive cut

For EVERY integer N>=2^20 this SAME family satisfies

 C_B2(F_N)>=4N-(2872021921/250000000)N^(4/5)
       +(6552747/1000000)N^(3/5)-(196601707/10000000)N^(2/5)
       +(733109/200000)N^(1/5)-37.                         (9.1)

We give the complete parameter, adjacent-choice and normalization proof.
All constants below are exact rationals; radicals are analysis devices,
not advice or operations missing from the fixed integer selector.

### 9.1. Parameters at every integer length

Write x=N^(1/5)>=16, q=1028118367/1000000000, q0=51/50, and put

 a=1301401524908/1000000000000, b=q*a+1/1000000000000,
 a0=130157859/100000000, b0=q0*a0+1/4000=6639300809/5000000000.

For2^20<=N<2^400 choose the least prime k>=ceil(a0*x^3).
The integer argument is in[2^12,2^241]: a0*16^3>2^12 and
ceil(a0*2^240)<2^241. Section10's checked chain gives k<=q0*ceil(a0*x^3).
Since q0/16^3<1/4000, v=k/x^3 lies in[a0,b0].
For N>=2^400 choose k>=ceil(a*x^3) least prime. Section8.1 gives
k<=q*ceil(a*x^3), so v lies in[a,b] since q/2^240<1/1000000000000.
Exact arithmetic gives a<a0<b0<b and13/10<a<b<67/50.

For z among a,b,a0,b0 define uz=(floor(10^14 sqrt(z))+1)/10^14.
Consecutive rational squares bracket each sqrt(z); in particular
ua=114078986886631/10^14. With f(v)=1/(2v^2)+sqrt(v), exact comparisons give

 1/(2a^2)+ua<2872021921/2000000000,
 1/(2(qa)^2)+ub<2872021921/2000000000,
 1/(2a0^2)+ua0<14360083/10000000,
 1/(2(q0*a0)^2)+ub0<14360083/10000000.

The reciprocal denominators in the right endpoints are smaller than
b,b0 respectively, hence these are upper bounds there. On(1,2),
f''(v)=v^(-4)(3-v^(5/2)/4)>0. Convexity therefore proves
8f(v)<2872021921/250000000 on BOTH entire prime-selection intervals.
These endpoint comparisons and all chain coverage checks execute exactly.

Put r=ceil(N/k), s=sqrt(r), and w=k/s^3. Since r>=x^2/v and
r<=x^2/v+1,

 s>13, k<N, r+1<=x^2/v+2< (13/10)x^3<=k,
 wl=19169/10000<w<wh=20709/10000.

For the last assertion use the exact squared comparisons
wl^2<a^5/(1+b/256)^3 and b^5<wh^2. The r+1 inequality follows from
10/13+2/256<(13/10)*16; also256/(67/50)>13^2.
Thus k is odd prime, k>s^3, and both t=j=floor(s) and t=j+1 satisfy
2<=t<r. Let

 H=r(r+3)/2, U=rk-r(r+1)/2>0,
 T_t=(t-1)k-t(t+1)/2+1,
 E_ref(t)=max(T_t,floor(U/(t+1))), d_ref(t)=H+E_ref(t)+1.

The supplied integer-selection proof in disperser.md gives
 d_new(t)<=d_ref(t). Both choices also have B=t(t+3)/2>=r/2.

### 9.2. Adaptive quotient comparison

Define

 c=1/8+1/(4w)-1/(32w^2), h0=191/1000,
 eta=1/2+1/(4w)-c/s, beta=1-eta, beta_min=36958/100000,
 S=ks-k/2-s^3/4+c*k/s-r/2.

On the whole enclosure,0<c<1/4,1/2<eta<3/4,
and0<beta_min<1/2-1/(4w)<beta<2/5. To verify these guards,
c'(w)<0 for w>1,0<c(wh)<c(wl)<1/4, s>13>wh, and
1/2-1/(4wh)+(1/4)/13<2/5. All are exact rational comparisons.

For a denominator s+beta, exact division gives

 U/(s+beta)=ks-k*beta+k*beta^2/s-s^3/2+beta*s^2/2-R,
 R=[(1+beta^2)/2+w*beta^3]*s^2/(s+beta).

The difference between the expression before R and S is exactly rF,
where F=w*((1/2+c/s)^2-2c)+1/2-1/(16w). The chosen c makes
w*(1/4-2c)+1/2-1/(16w)=0, so rF=wc*s+wc^2.
Put bstar=1/2-1/(4w) and A(w)=(1+bstar^2)/2+w*bstar^3.
Then R>A(w)(s-beta) and direct expansion gives

 A(w)-wc=3/16+1/(64w^2)>h0,
 A(w)=w/8+7/16-1/(32w)+1/(64w^2)<7/10.

The last expression increases for w>1 and is below7/10 at21/10.
Also wc^2<21/160 and A(w)*beta+wc^2<7/25+21/160<1.
Consequently R-rF>h0*s-1, proving the BASE denominator bound

 U/(s+beta)<S-h0*s+1.                                    (9.2)

### 9.3. Move the cut to exploit both adjacent branches

Put z=339/1000, delta=z/(w*s^2), eta1=eta-delta and beta1=beta+delta.
Then1/2<eta1<3/4 and beta1<1/2: explicitly

 1/(4wh)-(1/4)/13-z/(wl*13^2)>0,
 2/5+z/(wl*13^2)<1/2.

Changing the denominator lowers the quotient by

 D=U*delta/[(s+beta)(s+beta1)]
  =z*[s^3-s^2/(2w)-1/(2w)]/[(s+beta)(s+beta1)].

The numerator without z exceeds(s-2)(s+1/2)^2, because their difference is
s^2*(1-1/(2w))+(7/4)*s+1/2-1/(2w)>0. The denominator is below
(s+1/2)^2, so D>z(s-2)>zs-1. Equation(9.2) therefore yields

 U/(s+beta1)<S-(53/100)s+2,                              (9.3)

since h0+z=53/100.

Let e=s-j in[0,1). If e<=eta1 choose t=j. Its quotient has
j+1=s+1-e>=s+beta1, so(9.3) applies. Its triangular term is at most
k(s-1). The gap of S-(53/100)s+2 above this term is positive:

 s^3*(w/2-1/4)+(wc-1/2)*s^2-(53/100)s+2
 >(7/10)s^3-(1/2)s^2-(53/100)s+2>0 for s>13.

If e>=eta1 choose t=j+1. Its exact triangular expression is
T(e)=k(s-e)-(s+1-e)(s+2-e)/2+1.
Throughout[eta1,eta], its derivative -k+s+3/2-e lies in(-k,0).
Thus T(eta1)<T(eta)+k*delta=T(eta)+zs. Exact substitution gives
T(eta)=S-(1/2+beta)s-beta(1+beta)/2+1.
Using beta>beta_min and T decreasing, we obtain
T(e)<S-(53058/100000)s+1<S-(53/100)s+2.
For its quotient, j+2>=s+1 implies

 U/(j+2)<=ks-k+k/(s+1)-s^3/2+s^2.

Subtract this from S-(53/100)s+2. The gap is
k/2-k/(s+1)+ck/s+s^3/4-(3/2)s^2-(53/100)s+2>0,
since s>13 and k>0. Both branches prove that the SAME finite selector
contains a t with

 E_ref(t)<=S-(53/100)s+2.                               (9.4)

### 9.4. Admissibility and normalized gate count

Admissibility need not be assumed. Formula(9.4) is below ks+2, because
c<1/4,s>13 and all other nonconstant corrections are negative.
The elementary ceiling estimates give

 H<=N^2/(2k^2)+(5/2)N/k+2,
 ks<=sqrt(Nk)+k^(3/2)/(2sqrt(N)).

Therefore d_ref<=H+ks+3 and d_ref/N^(4/5)<2: f(v)<1.437,
5/(2v)+v^(3/2)/2<3 on[a,b], x>=16, and
1.437+3/256+5/16^4<2. Since2N^(4/5)<N, we have d_new<=d_ref<N.
Also B<=H<d_new. These are admissible p=0 candidates of the actual selector.

Expand ck/s=k/(8s)+r/4-s^5/(32k). Since r>=N/k, (9.4) and the
same ceiling estimates imply

 d_new<=f(v)N^(4/5)-K1(v)N^(3/5)+g8(v)N^(2/5)
                 -(53/100)*N^(1/5)/ub+5,
 K1(v)=v/2+1/(4v^(3/2)),
 g8(v)=9/(4v)+(5/8)v^(3/2)-1/(32v^(7/2)).

Here -s^5/(32k)<=-N^(5/2)/(32k^(7/2)); retaining this sign is essential.
K1 increases for v>1 and its exact root bracket gives
8K1(v)>8*(a/2+1/(4a*ua))>6552747/1000000.
Further,

 g8''(v)=(9/2)v^(-3)+(15/32)v^(-1/2)-(63/128)v^(-11/2)>0 (v>1),
 9/(4a)+(5/8)*a*ua-1/(32a^3*ua)<G8=264436809/100000000,
 9/(4b)+(5/8)*b*ub-1/(32b^3*ub)<G8,
 8G8-2/b<196601707/10000000,
 8*(53/100)/ub>733109/200000.

The root upper bounds give upper bounds even in the negative reciprocal
term. Convexity controls every v in[a,b]. Finally use4B>=2N/k and
Section4's gate theorem4N-8d_new+4B+3. This proves(9.1); its constant
is -8*5+3=-37, and every N in the asserted domain has been covered.

### 9.5. All lengths and retained consequences

For EVERY integer N>=0 the family also satisfies
C_B2(F_N)>=4N-(11488091/1000000)N^(4/5).
At N>=2^20 this follows from(9.1): with x>=16, its positive correction
K*x^3-D*x^2+P*x-37 is positive, where K,D,P are its rational constants.
For0<=N<=195 the zero option suffices, by
4000000^5*195<=11488091^5<4000000^5*196, and N=0 separately.
The remaining range196,...,2^28-1 is covered by the supplied3784-row
finite certificate, overlapping the analytic domain. The exact checker
verifies primality and all parameter guards, original sufficient dimension,
left>(r-1)k,right<=rk+k,d_ref<min(left,rk), adjacency, and at both endpoints u

 (1000000*cost(u))^5<=11488091^5*u^4,
 cost(u)=8d_ref-4B-3+4max(0,u-rk).

For each interior integer N, m=min(N,rk),p=N-m in[0,k] is an eligible
candidate, d_new<=d_ref<m and its certified loss is at most cost(N).
Concavity of(11488091/1000000)u^(4/5)-cost(u) transfers the endpoint
inequalities to the interior. First left196, adjacency and final right
268435455 exhaust the range. No finite sample is extrapolated.

The asymptotic bound is therefore4N-O(N^(4/5)) on EVERY sufficiently
large integer length, with one fixed uniform family and empty-witness NP
verifier. At N>=2^27, x>40 and K-D/40-37/40^3>6, so the retained
+6N^(3/5) corollary follows. At N>=2^65, x>=8192 and
K-D/8192-37/2^39>655029/100000 gives the retained late corollary.
All inequalities are exact; the positive N^(1/5) term may be discarded here.

### 9.6. Attribution of this consolidated argument

The original field/disperser/gate and prime/finite-cover framework is the
complete f663633a-c1f8-4922-8d4b-267eb7eb0a06 package, whose provenance
is retained in Sections7 and10. Enclosing guards and constants were developed
in e7a21b64-d265-44ea-8811-c0c344e251b8; exact division remainder in
6685cf68-24e8-485e-bbc7-d1aa64d9fb71; adaptive cancellation in
2d1a29e3-93be-4d08-b266-740c6c25c85e. The moved cut was delivered in f4aa34ff-acca-4351-9b39-de344cd74587.
The present extension tightens the existing prime ratio by its exact error
slack, then recomputes all root endpoints. The old all-length coefficient
is retained for the small finite cover. No target or kernel changes.
This chapter consolidates those arguments rather than relying on their
verdicts. Historical provenance elsewhere describes its source's development.
The target evaluator, full circuit model and leading projection remain unchanged.

## 10. A finite prime chain with independently checked primality witnesses

The supplied computation generates a finite list p0<...<pM with

 2^12<=p0<=(51/50)2^12,
 p_(i+1)<=(51/50)p_i,
 pM>=2^241,

and verifies that every entry is prime. Therefore EVERY real h in
[2^12,2^241] has a prime in[h,(51/50)h]: use p0 when h<=p0, and
the first chain entry at least h otherwise. Its predecessor is below h,
so its ratio bound implies the result. This is a finite coverage proof,
not an extrapolation of prime gaps beyond the last entry.

The chain contains8222 entries:698 small primes checked by independent
trial division, and7524 entries with Proth witnesses. A small-prime entry
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

The finite-cover/concavity method builds on d10eb869; generated arithmetic
certificates on 5ec6176f and 2dc9a80c; the elementary factorial/prime route
on original39c9e977, also delivered byte-identically as14c284dd; balanced
certificate generation on049df349; finite-tree errors on7ae3664d;
and coupled estimates/Proth-chain coverage on6927dad0.
The retained nu8 lists and count cache come from23856b84, with rational
parameter refinement from7099a2f9. Rational-envelope verification builds
on240f549a, direct simultaneous induction on8131ff33, period refinement
onaf5ec9e9 and exact difference checking on9cff9e1c.
The nearest-square-root saving and positive3N^(3/5) correction are
retained from ad149b6e-cc4a-4a29-b09a-fc9d3dce9371.
The credited720720 kernels introduced the17/20 error estimate, retained
by the prior1441440 refinement. The preceding mixed-period proof instead
used7/8 and9000. The intervening1021020 pair used9/10 and6000;
the present3063060 lists use7/8 and18000.
The original field/algebraic/circuit sources are credited in Section7.
All data and proof premises are supplied; no optimizer, verdict or
unverified attribution decision is a mathematical premise.

## 11. A supplied analytic prime theorem and the limiting error coefficient

This section proves for the SAME fixed family

 for every eta>0 there is N_eta such that, for ALL integers N>=N_eta,
 C_B2(F_N)>=4N-(10*2^(1/5)+eta)*N^(4/5).                (11.1)

In particular C_B2(F_N)>=4N-(11487/1000)N^(4/5) for every sufficiently
large integer N. No numerical threshold for this last assertion is claimed.
This strengthens the secondary asymptotic error coefficient, not the native
leading coefficient4. The value10*2^(1/5) minimizes the particular f-envelope
used here; no optimality assertion about actual circuits or other proofs follows.

We supply the needed prime distribution argument rather than treating an
external theorem as evidence. Its analytic route is classical: D.J.Newman's
Tauberian argument, as presented by D.Zagier, American Mathematical Monthly
104(8),1997,705--708, DOI10.2307/2975232, especially the analytic lemma and
the theta-asymptotic steps. Source PDF:
https://people.mpim-bonn.mpg.de/zagier/files/doi/10.2307/2975232/fulltext.pdf
SHA256 a2cdbbcc38b9f21cd68903bcd0262be2c211677cbda7ae903b9e124bd0d12561.
The proof below is in our notation; no external analytic-number-theory lemma
is left unproved. Ordinary Cauchy/residue theory, locally uniform convergence
of holomorphic functions, and dominated integration are used with explicit
hypotheses. These standard complex-analysis facts do not concern circuit size.

### 11.1. Bounded Laplace transforms: the contour argument

Let u:[0,infinity)->C be bounded and locally integrable, |u|<=M, and
G(z)=integral_0^infinity u(t)e^(-zt)dt for Re(z)>0. Suppose G extends
holomorphically to an open neighborhood of the closed right half-plane.
Then integral_0^T u(t)dt tends to G(0) as T tends to infinity.

Here are the details. G_T(z)=integral_0^T u(t)e^(-zt)dt is entire by
uniform differentiation on compact sets. Fix any R>0. Compactness of the
imaginary segment[-iR,iR] gives a thin strip to its left where G is analytic.
Take a rectilinear path J from iR to-iR in that strip, with Re(z)<0 except
at its endpoints, staying to the left of0. Together with the right semicircle
Cplus of radiusR it bounds a region where G is holomorphic and containing0.
Cauchy's formula applied separately to G_T on the full circle and to G on
Cplus+J gives

 G_T(0)-G(0)=(1/(2*pi*i))*(
  integral_Cplus (G_T-G)(z)e^(Tz)(1+z^2/R^2)dz/z
 +integral_Cminus G_T(z)e^(Tz)(1+z^2/R^2)dz/z
 -integral_J G(z)e^(Tz)(1+z^2/R^2)dz/z),

where Cminus is the left semicircle. All contours are positively oriented;
the only kernel pole inside is0, with residue equal to the displayed value.
On Cplus, |(G_T-G)e^(Tz)|<=M/Re(z), by the tail integral. On Cminus,
|G_Te^(Tz)|<=M/|Re(z)|, by integrating e^((T-t)Re(z)) for0<=t<=T.
For |z|=R,

 |1+z^2/R^2|=2*|Re(z)|/R.

Thus each semicircle integrand has absolute value at most2M/R^2 away
from its endpoints, and contributes at most M/R after division by2pi.
These bounds are integrable up to the endpoints; their endpoint values do
not affect the integrals. On the fixed path J the kernel times G is bounded
and Re(z)<0 in its interior. Dominated convergence makes its integral tend
to0. Consequently limsup_T |G_T(0)-G(0)|<=2M/R. Since R was arbitrary,
the asserted convergence follows. No uniform-width extension for all R is
assumed: the strip and J may depend on the fixed R.

### 11.2. Zeta and the zero-free boundary

For Re(s)>1 set Z(s)=sum_(n>=1)n^(-s) and
P(s)=sum_(p prime)(log p)p^(-s), using real logarithms. These series and
their needed derivatives converge locally uniformly by comparison with
sum(log n)^j*n^(-sigma), sigma>1. Unique factorization and absolute
convergence give the Euler product and its logarithm

 Z(s)=product_p(1-p^(-s))^(-1),
 log Z(s)=sum_(p,m>=1)p^(-ms)/m.

The last series converges absolutely, so its exponential is nonzero and
agrees with Z. Differentiation gives

 -Z'(s)/Z(s)=P(s)+B(s),
 B(s)=sum_(p,m>=2)(log p)p^(-ms).

B is holomorphic for Re(s)>1/2: on any compact subset, its absolute sum
is bounded by a fixed multiple of sum_(n>=2)(log n)n^(-2sigma), sigma>1/2.
Termwise identities follow from local uniform convergence.

Integration of floor(x)=sum_(n>=1)1_(x>=n) gives for Re(s)>1

 Z(s)=s/(s-1)-s*integral_1^infinity {x}x^(-s-1)dx.

Since0<={x}<1, the integral and all its local complex derivatives converge
uniformly on compact sets with Re(s)>0. This continues Z meromorphically
there, with its only possible pole at1, where its residue is1.

Z has no zero at1+it for real t!=0. For sigma>1 the Euler logarithm gives

 log[Z(sigma)^3*|Z(sigma+it)|^4*|Z(sigma+2it)|]
 =sum_(p,m>=1) p^(-m*sigma)/m *
       [3+4cos(mt*log p)+cos(2mt*log p)]>=0,

because the bracket is2(1+cos(mt*log p))^2. If Z(1+it)=0, analyticity
would give |Z(sigma+it)|=O(sigma-1). The factor Z(sigma)^3 is
O((sigma-1)^(-3)), while Z(sigma+2it) stays bounded since2t!=0.
The product would tend to0, contradicting its lower bound1.
At s=1, write Z(s)=1/(s-1)+holomorphic; direct logarithmic differentiation
shows -Z'/Z-1/(s-1) is holomorphic nearby. Subtracting B now proves

 P(s)-1/(s-1) extends holomorphically near every point Re(s)>=1.

Away from1 use the zero-free property and continuity to take a zero-free
neighborhood of each boundary point. Local continuations agree with the
original series on Re(s)>1 and hence with each other on overlap. This is
exactly the neighborhood condition needed below, not a zero-free strip of
fixed width or a hypothesis about the Riemann zeros.

### 11.3. Theta asymptotics and primes in every fixed relative interval

Let theta(y)=sum_(p<=y)log p as in Section8, where theta(y)<4y was proved.
Define u(t)=e^(-t)theta(e^t)-1 for t>=0. It is bounded and locally integrable:
on every finite interval theta(e^t) has finitely many jumps. Absolute
integration and summation give, for Re(z)>0,

 integral_0^infinity u(t)e^(-zt)dt=P(1+z)/(1+z)-1/z.

The right side extends holomorphically to a neighborhood of Re(z)>=0.
At z=0 its two simple poles cancel by Section11.2; at other boundary points
there is no pole. Section11.1 therefore proves convergence as X tends to
infinity of

 integral_1^X [theta(y)-y]dy/y^2.

This convergence and monotonicity imply theta(X)/X tends to1. Indeed,
if theta(X)>=(1+epsilon)X at arbitrarily large X, then over
[X,(1+epsilon/2)X] its integrand has integral at least
 epsilon^2/[4(1+epsilon/2)^2]>0. Convergence requires integrals over such
intervals with both endpoints tending to infinity to tend to0, a contradiction.
If theta(X)<=(1-epsilon)X for some fixed0<epsilon<1 at arbitrarily large X,
integrate over[(1-epsilon/2)X,X]. Monotonicity makes the integral at most
-epsilon^2/4, again a contradiction. Thus the limit is1 over all real X.

For EVERY fixed real q>1 it follows that theta(qX)-theta(X)>0 for every
sufficiently large real X: choose epsilon<(q-1)/(q+1) and use
(1-epsilon)qX-(1+epsilon)X>0. A positive difference means that a prime lies
in(X,qX]. This proves the only new prime-distribution fact used below,
with no finite search or external retrieval as a premise.

### 11.4. Apply the prime theorem to the unchanged selector

The entire Chapter9 adjacent-choice and normalization argument applies to
ANY prime k whose v=k/N^(3/5) lies in its enclosing interval[a,b]; the
particular prime-selection routes there were only existence proofs.
Let vstar=2^(2/5). It lies strictly inside[a,b], and

 f'(v)=-v^(-3)+1/(2sqrt(v))=0 iff v=vstar,
 8f(vstar)=10*2^(1/5)=:Astar.

The derivative changes from negative to positive there, so this is the
minimum of this particular envelope for v>0.
Fix eta>0. By continuity, choose fixed rationals a1>0 and q1>1 so close
to vstar and1 respectively that[a1,q1*a1] lies inside(a,b) and
8f(v)<Astar+eta/4 on that interval. For all sufficiently large integers N,
the same strict inequalities hold on[a1,q1*a1+q1/N^(3/5)].
Set h=ceil(a1*N^(3/5)). Section11.3 supplies a prime k in[h,q1*h] at
every sufficiently large N. Thus its v lies in that smaller interval and
satisfies all Chapter9 guards. This is an analysis choice among parameters
already enumerated by the SAME fixed polynomial-time selector; neither the
family nor its NP verifier depends on eta,a1 or q1.

The selected certified loss L_N is no greater than this admissible witness
candidate loss. Section4 proves C_B2(F_N)>=4N-L_N for the actual selected
function; no transfer between different functions is being assumed. Hence
the Chapter9 candidate-loss estimate yields

 C_B2(F_N)>=4N-8f(v)N^(4/5)+K*N^(3/5)-D*N^(2/5)+P*N^(1/5)-37,

with its fixed positive constants K,D,P. Discarding the positive terms,
the loss beyond(Astar+eta/4)N^(4/5) is at most D*N^(2/5)+37, which is
less than(3eta/4)N^(4/5) for every sufficiently large N. This proves(11.1)
with all its quantifiers on one family. Thresholds may depend on eta.
The equivalent statement is

 liminf_(integer N->infinity) [C_B2(F_N)-4N]/N^(4/5)>=-10*2^(1/5).

Finally, (11487/10000)^5>2 is an exact rational inequality. It implies
Astar<11487/1000, so a fixed positive eta in(11.1) proves the stated
eventual rational11.487 error bound. No numeric or symbolic projection
other than the existing leading[4,1] is requested for this secondary theorem.

## Appendix A. Deterministic prime-degree fields over F2

This is a prime-degree specialization of the deterministic construction
route of Victor Shoup, *New Algorithms for Finding Irreducible Polynomials
over Finite Fields*, Mathematics of Computation54(189),1990,435--447,
https://www.shoup.net/papers/detirred.pdf, Sections2--3. The supplied
algorithm and argument below prove the particular statement used here;
no appeal to an external factoring or irreducible-construction oracle is
needed. All operations used by the constructor are specified below.

### A.1. Extract an irreducible factor using binary linear algebra

Polynomials are represented by binary integers, with bit j the coefficient
of X^j. Euclidean division, multiplication and gcd are exact finite bit
operations. On a squarefree monic polynomial f of degree D, form the
D-dimensional quotient algebra A=F2[X]/f. The map b->b^2 is F2-linear.
Its matrix has column j the coefficients of X^(2j) modf. Gaussian
elimination computes the kernel V of that matrix minus the identity,
with pivots and free columns in increasing order and the usual free-bit
basis. The factor algorithm returns f if dimV=1. Otherwise it selects
the first basis vector b other than0 or1, computes g=gcd(f,b), divides
h=f/g, and repeats on the numerically smaller of g,h. All choices are
fixed and every repeat strictly reduces degree.

For correctness, write the mathematical factorization f=product f_i
into distinct monic irreducibles. Existence follows by repeated division
in the Euclidean domain; this does not assume a factoring algorithm.
The Chinese remainder map A->product F2[X]/f_i is an isomorphism:
Euclid gives inverses of the pairwise coprime factors and hence explicit
idempotents for each component. Each component is a field, whose only
solutions of z^2=z are0 and1. Thus V corresponds exactly to one binary
choice per irreducible factor, and its dimension equals their number.
Dimension1 means f is irreducible. Otherwise a nonconstant b in V has
some component0 and some1. Therefore its gcd g includes exactly the
zero components, while h includes the others. Both have positive degree
less than D, are squarefree, and divide f. Repetition terminates at an
irreducible divisor. A nonconstant basis vector exists whenever dimV>1,
since the constants form a one-dimensional subspace. The construction
requires no enumeration of quotient-field elements.

Each matrix has D rows of D bits. Its entries and Gaussian elimination
are computable in polynomial time. There are at most D descents. A
coarse O((D+1)^5) bit-operation bound covers this elementary version;
no fast factorization algorithm is required.

### A.2. A cyclotomic field and a nonresidue

Let k>=3 be odd prime. Compute m, the least positive integer with
2^m=1 modk, by repeated multiplication modulo k. Fermat gives
m|k-1, so2<=m<=k-1 and gcd(m,k)=1. Set e=v_k(2^m-1), using repeated
exact integer division. In particular1<=e<=m. Start with the polynomial
Phi_k=1+X+...+X^(k-1), and extract an irreducible factor f_1 by Appendix A.1.
For i=2,...,e form f_(i-1)(X^k) and extract an irreducible factor f_i.
Every f_i has degree m, as shown next. At the end put f=f_e,
K=F2[A]/f and a=A modf.

We use the elementary cyclicity of the multiplicative group of a finite
field. For completeness, for every prime-power ell^u dividing its order
Q-1, some element has order divisible by ell^u. Otherwise every nonzero
element would be a root of X^((Q-1)/ell)-1, contradicting the root bound.
Taking the appropriate powers isolates such an ell-primary element.
The product of these elements for the different primes has order Q-1:
commuting elements of coprime orders have product order equal to their
product, as follows by raising to each order. Lagrange's theorem supplies
the divisibilities used here. Thus a finite field contains all d roots
of unity when d divides Q-1, and primitive ones exist.

The roots of Phi_k are precisely the primitive k-th roots of unity. Their
Frobenius orbits have size ord_k(2)=m, so every irreducible factor of Phi_k
has degree m. This polynomial is squarefree since X^k-1 has derivative
X^(k-1) in characteristic2 and does not vanish at its nonzero roots.
Inductively a root beta of f_(i-1)(X^k) has beta^k primitive of order
k^(i-1), so beta has order exactly k^i. For i<=e, k^i divides2^m-1, so
all these roots are in the field of2^m elements. Their orbit length is
ord_(k^i)(2), which divides m and is divisible by ord_k(2)=m, hence is m.
The composition is squarefree: its roots are nonzero, k is odd, and the
original irreducible polynomial has simple roots. Thus every factor
extracted by Appendix A.1 has the required degree and primitive-root order.
In particular a has order k^e in K. It is not a k-th power in K: a k-th
root would have order k^(e+1), impossible since that does not divide
|K*|=2^m-1. These statements apply also to e=1, with no lifting loop.

All factored polynomials have degree at most mk<k^2. There are e<=k
factor extractions. Neither an integer factorization of2^m-1 nor an
exponentially large polynomial is used; only division by the given k
and polynomials of degree less than k^2 occur.

### A.3. A prime extension and descent by trace

The polynomial Y^k-a is irreducible over K. To prove this directly, take
a root beta in an algebraic closure. K contains the k-th roots of unity.
Put q=2^m. Since beta^k=a and a^q=a, beta^q/beta is a k-th root of unity.
It cannot be1, as that would put beta in K (the roots of Z^q-Z are exactly
K), contradicting the nonresidue property. As k is prime it has order k.
Frobenius^m fixes that ratio. The conjugates beta,beta^q,... therefore
form an orbit of length exactly k. Their minimal polynomial over K has
degree k and divides Y^k-a, proving irreducibility.

We may consequently calculate in E=K[Y]/(Y^k-a), representing elements
as k coefficients in K. Its dimension over F2 is mk, which is less than
k^2. Let beta be the residue of Y. Define the explicit element

    gamma=sum_(i=0)^(m-1) beta^(2^(ki)).

The automorphism sigma:z->z^(2^k) permutes this sum cyclically, because
sigma^m is the identity on E. Thus gamma^(2^k)=gamma. Also gamma is not
in K. Indeed write2^(ki)=k u_i+r_i with1<=r_i<k. The residues r_i are
distinct as i runs0,...,m-1: Fermat gives2^k=2 modk, so they are exactly
the first m powers of2 modulo k. Therefore in the K-basis1,beta,...,
beta^(k-1),

    gamma=sum_i a^(u_i) beta^(r_i)

has distinct nonzero positive-degree terms, with no constant term.
It cannot lie in K. In particular gamma is not in F2. Since its binary
Frobenius orbit has size dividing prime k, that size is exactly k.
Hence its minimal polynomial over F2 is

    Q_k(Z)=product_(j=0)^(k-1)(Z+gamma^(2^j)).

The product is computable using represented E arithmetic. Frobenius
permutes the factors, so each coefficient is fixed by squaring and is
therefore0 or1. It is monic, degree k and irreducible, because gamma has
orbit length k. This defines one completely deterministic Q_k for every
odd prime k, with no field-polynomial input, advice, random choice,
number-theory conjecture or external algorithm.

The represented arithmetic never enumerates E. A K element has m bits,
and an E element has k such coefficients. Reduce Y^k to a in products.
For squaring, characteristic2 kills all cross terms and a term of degree
2i>=k needs just one such reduction. Compute successive sigma powers by
k squarings, repeated m times, rather than writing enormous exponents.
Compute the final product with ordinary polynomial multiplication over E.
These exact procedures use polynomially many bit operations. Combining
Appendix A.1--A.3, the generous bound O((k+1)^12) covers the implementation;
its largest factor matrix has fewer than k^2 rows and columns. Only the
selected degree's Q_k need be constructed by the language evaluator.


### A.4. Retained compressed-prefix reference class

Let k>=5 be odd prime, K a field of2^k elements by the elementary finite-field facts in proof.md Section2, and
Tr(z)=sum_(i=0)^(k-1)z^(2^i). Tr takes values in F_2 and Tr(1)=1,
so its kernel has M=2^(k-1) elements. There exists nonzero alpha with
Tr(alpha)=Tr(alpha^-1)=0. Suppose otherwise. The monic polynomial

    T(X)=Tr(X)/X=sum_(i=0)^(k-1)X^(2^i-1)

has precisely the M-1 distinct nonzero trace-zero elements as roots.
All would then have inverse trace1. Thus T divides the polynomial

    U(X)=X^M(Tr(1/X)+1)=X^M+sum_(i=0)^(k-1)X^(M-2^i).

U is monic degree M, while T has degree M-1, so the quotient is X+c.
The coefficient of X^(M-1) forces c=1. But U has coefficient1 at
X^(M-2), while (X+1)T has coefficient0 there. Its exponents are powers
of2 or powers of2 minus1; M-2 is neither when M>=8. This contradiction
proves the claim using polynomial division and the elementary root bound.

Neither alpha nor alpha^-1 is1 since Tr(1)=1. Replace alpha by alpha^-1
if necessary to arrange also Tr(1/(alpha+1))=0. This is always possible:

    Tr(1/(alpha^-1+1))=Tr(alpha/(alpha+1))
                       =1+Tr(1/(alpha+1)).

The first two zero traces are preserved under inversion. Since k is
prime and alpha is not in F_2, its Frobenius orbit has size k. Its monic
irreducible polynomial Q has constant1, coefficient of Z^(k-1) equal
Tr(alpha)=0, and coefficient of Z equal Tr(alpha^-1)=0. These follow by
expanding the product over the conjugates; the product of conjugates is
nonzero in F_2, hence1. Also Q(1)=1 since Q has no root1. The logarithmic
derivative identity, obtained by differentiating that product, gives

    Q'(1)/Q(1)=sum_(i=0)^(k-1)1/(1+alpha^(2^i))
              =Tr(1/(1+alpha))=0.

Consequently Q+1 is divisible by Z^2 (its constant and linear coefficients
vanish) and by (Z+1)^2 (its value and derivative at1 vanish). These
coprime factors give

    Q=1+Z^2(Z+1)^2 R=1+(Z^4+Z^2)R.

R is monic degree k-4. Its coefficient at degree k-5 is0 because Q's
coefficient at degree k-1 is0. All remaining coefficients are exactly
k-5 free bits, as in the decoder. Hence every eligible k has at least
one encoded irreducible polynomial. Multiplication by Z^4+Z^2 is
injective in F_2[Z], so different bit prefixes encode different Q. The
prefix test from proof.md Section2 remains exact for every decoded monic Q.
This proof constructs no field polynomial by advice or an external
algorithm: it proves existence for restricting circuits; evaluation just
decodes its actual input and checks it in polynomial time.

The complete published refined-threshold carrier c6e8ef0a-3f72-4b51-ab27-7edc2fe2ef58 supplies retained foundations and earlier quantitative arguments.
The original720720 source d6babd50-b62e-40ea-bfaa-b8334d0d9f95 and1441440 source8436d193-a7be-4aa4-b7fc-146c0ed108aa supply actual discovery step indices.
The mixed sourceea82aa89-47fe-4043-a109-65fc99a2bc4d supplies additional actual seed indices and the prior discovery bounds. Original field/algebraic/circuit sources remain credited in Sections7/10; no verdict or pending graph decision is a premise.

The immediate complete source80251948-bc9e-41de-9360-2d4c14f48540 supplies
the retained family/proof framework. The present3063060 data, biased
adjacent-choice estimate, tighter middle chain and continuous correction
are independently supplied refinements. Small-block checking preserves
every cell; larger-wheel screening and generation strides are only hints.

The shifted-choice correction builds directly on own complete source a0f5f28f-6d5f-4481-be91-2472f5c7770c, whose selected literal data and evaluator are unchanged. Its publication-dependent attribution edge will be added when available.

The earlier-length extension builds directly on own complete shifted-choice source3c7430ac-7ee5-4960-bb42-7178b539e682; its direct attribution edge is deferred until publication.

This denser-chain extension builds directly on complete own source5777a2bc-d991-4c66-bd60-24b925b57ce4; its attribution edge is deferred until publication.
