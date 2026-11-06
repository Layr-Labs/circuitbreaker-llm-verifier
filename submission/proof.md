# Explicit coefficient-four bound at every actual integer length

For the fixed deterministic-field evaluator from
37b2928a-14a5-4e21-a621-12723a6848ff, the new Sections8--9 prove
C_B2(F_N)>=4N-(25/2)N^(4/5) at EVERY integer actual length N>=0.
The proof combines an explicit Bertrand-based analytic tail with a
220-interval finite certificate. The evaluator is unchanged and the
bounded audit now checks that certificate. All construction and full-B2
arguments needed for the theorem are supplied below. The native pair is
(4,1), and no actual circuit upper bound is asserted.

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
The new contribution is the explicit every-length bound in Sections8--9.
Its proof uses the original sufficient dimension and a supplied finite
certificate, without a sharp-asymptotic-constant theorem as a premise. Literature for the retained construction and gate argument
remains BK and Mishra--Santhanam ECCC TR26-260 revision1. Unfamiliar
arguments needed for the theorem are reproduced; ordinary field algebra,
root counting, Euclid, Gaussian elimination and Bertrand are used as stated.


## 8. A fully explicit all-length bound without PNT

For the same unchanged evaluator we also prove the concrete guarantee

    N>=2^18 implies L_N<(25/2)N^(4/5),
    C_B2(F_N)>=4N-(25/2)N^(4/5).                         (8.1)

Here N is every integer at least262,144. This section uses only
Bertrand's postulate for its prime comparison. The result gives an explicit threshold; it
does not assert that this threshold is the smallest possible one.

Let b be the least integer with b^5>=N^3, and choose the least prime
k>=b. Bertrand gives b<=k<2b<=2N^(3/5)+2. Put

    u=k/N^(3/5), r=ceil(N/k), t=ceil(sqrt(r)), p=0,m=N.

For N>=2^18,

    1<=u<2+2N^(-3/5)<=2+2/1024<1001/500.                 (8.2)

These are eligible base parameters. First k<N since N^(2/5)>144
and k/N<(1001/500)N^(-2/5). Also r>=N/k>(500/1001)144>3.
We have r+1<=N^(2/5)+2<N^(3/5)<=k: set x=N^(1/5)>12 and use
x^2+2<x^3. Therefore k>r+1, and k is an odd prime at least5.
For integer r>=3,2<=ceil(sqrt(r))<r, so the order guard holds.

Let d be the exact dimension for this candidate. The retained original
sufficient threshold, with its proved dominance, gives

    d<=H+max(T,floor(U/(t+1)))+1<=H+k sqrt(r)+1,

since t-1<sqrt(r), T<k sqrt(r), and U/(t+1)<=rk/(t+1)<k sqrt(r).
Using r<=N/k+1 and sqrt(1+z)<=1+z/2 for z>=0,

    H<=N^2/(2k^2)+(5/2)N/k+2,
    k sqrt(r)<=sqrt(Nk+k^2)
                 <=sqrt(Nk)+k^(3/2)/(2sqrt(N)).

Consequently

    d/N^(4/5)<=f(u)+g(u)N^(-2/5)+3N^(-4/5),
    f(u)=1/(2u^2)+sqrt(u),
    g(u)=5/(2u)+u^(3/2)/2.                              (8.3)

On the interval[1,1001/500], f is convex, since

    f''(u)=u^(-4)(3-u^(5/2)/4)>0.

For example u<9/4 and (9/4)^5<144 prove the strict inequality.
Its maximum is therefore at an endpoint. At1 it equals3/2;
at1001/500 it is strictly less than

    1/8+283/200=77/50,

because1001/500>2 and sqrt(1001/500)<283/200. The function g is also convex, since g''(u)=5/u^3+3/(8sqrt(u))>0.
Its endpoint values satisfy g(1)=3 and

    g(1001/500)<5/4+(1001/500)(283/200)/2<3.

Thus g(u)<=3 throughout the interval. For N>=2^18,
N^(2/5)>=2^(36/5)>144: raising to the fifth power reduces this last
comparison to2^16>9^5. Hence N^(-2/5)<1/144 and N^(-4/5)<1/144^2.
Substitution in(8.3) gives

    d/N^(4/5)<77/50+3/144+3/144^2
              =77/50+1/48+1/6912<25/16.

Indeed25/16-77/50-1/48=1/600>1/6912. In particular d<N, because
N^(1/5)>12>25/16. Also B<=H<d, so the candidate is admissible in
the selector's full-payload r0 branch. Its loss is8d-4B-3, strictly
less than(25/2)N^(4/5). The global minimum is no larger, and the finite
circuit inequality from Section4 proves(8.1).

All ceilings and prime choices in this comparison have exact integer
definitions and are contained in the fixed selector. The elementary
real inequalities above supply a proof for the entire infinite range;
no finite numerical scan is used to infer an untested interval. Field
and coefficient consequences continue to apply only to the selected
PRIME degree. This explicit
threshold does not change the evaluator or asserts an actual circuit upper bound.


## 9. Finite interval certificate extends the bound to every N

For every integer N>=0, the same evaluator satisfies

    L_N<=(25/2)N^(4/5),
    C_B2(F_N)>=4N-(25/2)N^(4/5).                        (9.1)

For0<=N<=298 the zero option gives L_N<=4N, which is at most the
stated loss: for positive N raise the inequality to the fifth power
and use8^5*N<=8^5*298<=25^5. The N=0 case is equality. Section8 covers
all N>=2^18. It remains to certify the finite range299,...,2^18-1.

The literal FINITE_LOSS_COVER list in audit.py supplies220 integer rows
(left,right,k,r,t), in increasing consecutive intervals, beginning at299
and ending at262143. For each row the declared audit verifies primality
of k, r>=3,k>r+1,2<=t<r, and puts

    H=r(r+3)/2, B=t(t+3)/2,
    d_ref=H+max(T,floor(U/(t+1)))+1,
    A=8d_ref-4B-3, cost(N)=A+4max(0,N-rk).

It checks left>(r-1)k, right<=rk+k, d_ref<min(left,rk), B<=d_ref,
and verifies at BOTH endpoints v in{left,right} the exact integer test

    (2cost(v))^5<=25^5*v^4.                             (9.2)

All costs are nonnegative. Thus(9.2) is equivalent to
cost(v)<=(25/2)v^(4/5), without floating-point arithmetic.

For every N in the row interval, m=min(N,rk) has ceil(m/k)=r,
p=N-m lies in[0,k], and d_ref<m. The selector's two-count reduction
therefore includes this parameter choice. Its actual threshold obeys
d_new<=d_ref, with B<=H<d_new, so it is admissible and

    L_N<=4p+8d_new-4B-3<=cost(N).

This comparison is numeric; it does not substitute a row's function
for the fixed selected evaluator. The finite selected-family circuit
bound of Section4 still uses that evaluator's own L_N.

On the positive real interval[left,right], the function

    g(x)=(25/2)x^(4/5)-A-4max(0,x-rk)

is concave: its power term is concave, and minus the maximum of two
affine functions is concave. Nonnegative endpoint values therefore imply
nonnegative values throughout that interval. The two integer tests prove
the required inequality for every integer in the row, not just endpoints.
The audit verifies that each row begins exactly one above its predecessor's
right endpoint, so the220 rows cover the entire remaining finite range.
Together with the zero-option range and Section8 this proves(9.1).

The finite certificate does not claim optimal intervals, optimal reference
parameters, or that numerical samples prove an untested infinite range.
Every finite endpoint check is reproduced by the declared credential-free
experiment, while concavity, exact candidate inclusion and the analytic
tail argument provide the stated mathematical coverage. The field degree
of every row is PRIME, with the same supplied construction and full-B2
proof; no rejected composite-degree claim is reused.
