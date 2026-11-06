# Sharp certificate-loss constant for the deterministic-field family

For the unchanged fixed evaluator from37b2928a-14a5-4e21-a621-12723a6848ff,
Section8 proves L_N/N^(4/5)->10*2^(1/5) over all integer lengths. The
same family therefore has full-B2 complexity at least
4N-(10*2^(1/5)+epsilon)N^(4/5) at every sufficiently large N for each
epsilon>0. Sharpness concerns this selector's certificate arithmetic,
not actual circuit upper bounds. Its native leading pair remains(4,1).
The complete field construction, finite selection and circuit proofs are
retained below; the evaluator and bounded audit are unchanged.

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

The deterministic field constructor, its complete proof, prefix optimizer
and Section13 count-loop reduction are retained from
37b2928a-14a5-4e21-a621-12723a6848ff. Its construction specializes Shoup1990;
field.md supplies the algorithm and proof rather than an external oracle.
The current contribution is Section8's sharp asymptotic certificate-loss
constant for this unchanged evaluator. It does not claim a new general
irreducible-construction algorithm or actual circuit-size optimality. Literature for the retained construction and gate argument
remains BK and Mishra--Santhanam ECCC TR26-260 revision1. Unfamiliar
arguments needed for the theorem are reproduced; ordinary field algebra,
root counting, Euclid, Gaussian elimination and Bertrand are used as stated.

## 8. Sharp asymptotic constant of this selector's certified loss

Let C_star=10*2^(1/5). For the EXACT fixed selector and threshold of
Section1, not for optimal actual circuits, we prove

    lim_(N->infinity, N integer) L_N/N^(4/5)=C_star.         (8.1)

Consequently for every real epsilon>0, at all sufficiently large integer
N the same language satisfies

    C_B2(F_N)>=4N-(C_star+epsilon)N^(4/5).                  (8.2)

This does not assert an upper bound on actual circuit complexity or
optimality of this remainder over other families or proof methods.
The universal quantitative improvement is (8.2); the equality in (8.1)
concerns only the declared computable certificate loss.

We use the familiar prime number theorem in its standard form
pi(x)~x/ln(x), with the natural logarithm, as x tends to infinity over
positive reals. For every fixed
eta>0 it implies pi((1+eta)x)-pi(x)~eta*x/ln(x)>0. Hence the least prime
at least x has ratio to x tending to1. This is the only additional
external standard theorem used in this section. All unfamiliar field
and circuit components remain supplied elsewhere in the package.

### 8.1 Elementary bounds on the exact integer maximum

For every guarded k,r,t, the proofs in disperser.md give

    E0<=max(T,floor(U/(t+1))), U=rk-r(r+1)/2,
    E0>=T-r+t-1=(t-1)(k-t/2)-r, E0>=0.                   (8.3)

The first follows from inclusion in the original sufficient threshold;
the second is Section12's fewer-than-t-eligible attaining vector.
There is also the useful lower bound

    E0>=r*floor((k-r)/(t+1)).                            (8.4)

To prove it, put x=floor((k-r)/(t+1)) and assign every weight w_j=j+1
the dimension d_j=w_j+x. The guard k>r+1 gives k-r>=2 and t>=2, so
r+1+x<=k. If x=0 the vector fails because it has no positive deficits.
If x>0, every t-subset has positive deficit sum tx and maximum dimension
at most r+1+x, whence

    tx+(r+1+x)<=(t+1)x+r+1<=k+1.

It therefore fails the enlarged predicate too. Its total dimension is
H+rx. The exact maximum theorem of Section12 proves (8.4). This is an
integer selection argument and makes no composite-degree field assertion.

### 8.2 All selected parameters lie on the relevant scales

Section5 already gives a constant M with L_N<=M N^(4/5) at every large
integer N. Such lengths use a nonzero option since4N is larger. Write
k,r,t,p for their selected parameters. Since B<=H and d=H+E0+1,

    L_N=4p+8H+8E0-4B+5>=4p+4H+8E0.

Thus p=O(N^(4/5)), r=O(N^(2/5)) and E0=O(N^(4/5)), uniformly for all
these selected choices. As t>=2 and t<r<k, (8.3) implies

    E0>=tk/4-r, hence tk<=4(E0+r)=O(N^(4/5)).

In particular k=O(N^(4/5))=o(N). The actual payload has m=N-p, and
0<=rk-m<k, so rk=N+O(N^(4/5))~N. Combined with the upper bound on r,
this gives k=Omega(N^(3/5)), and then t=O(N^(1/5)). Also r/k tends to0.
Using floor(z)>=z-1 in (8.4) gives

    E0+r>=r(k-r)/(t+1).

Here r(k-r)~N. Since E0+r=O(N^(4/5)), it follows that
t=Omega(N^(1/5)). The bound on tk now implies k=O(N^(3/5)), and
rk~N gives r=Omega(N^(2/5)). We have proved simultaneously

    k=Theta(N^(3/5)), r=Theta(N^(2/5)),
    t=Theta(N^(1/5)).                                    (8.5)

Finally the selector's exact prefix restriction gives p<=k, strengthening
p=O(N^(3/5)) and rk=N+O(N^(3/5)). No assumption about an optimum's
shape was made; these estimates follow from the actual global minimum.

### 8.3 Uniform asymptotic expansion and lower constant

On the scales (8.5), with p<=k and 0<=rk-m<k, elementary algebra gives

    H=N^2/(2k^2)+O(N^(2/5)),
    T=tk+O(N^(3/5)),
    U/(t+1)=N/t+O(N^(3/5)),
    r*floor((k-r)/(t+1))=N/t+O(N^(3/5)).

For the first, r=N/k+O(1). For the remaining estimates use
rk=N+O(k), N/(t(t+1))=O(N^(3/5)), r^2/(t+1)=O(N^(3/5)),
and k+t^2+r=O(N^(3/5)); rounding changes a floor by less than one.
Bounds (8.3)--(8.4) therefore sandwich the exact maximum as

    E0=max(tk,N/t)+O(N^(3/5)).

The implicit constants are uniform for the selected parameters, because
(8.5) gives fixed positive upper and lower scale bounds. Since
B=O(N^(2/5)) and p=O(N^(3/5)), the selected loss satisfies

    L_N=4N^2/k^2+8max(tk,N/t)+O(N^(3/5)).                 (8.6)

Put a=k/N^(3/5), b=t/N^(1/5). After division by N^(4/5),

    L_N/N^(4/5)=4/a^2+8max(ab,1/b)+O(N^(-1/5)).

For all positive a,b, max(ab,1/b)>=sqrt(a). Apply arithmetic-geometric
mean to the five positive terms4/a^2 and four copies of2sqrt(a):

    4/a^2+8sqrt(a)>=5*(64)^(1/5)=10*2^(1/5)=C_star.

Equality is possible at a=2^(2/5), b=2^(-1/5). Equation(8.6) proves
liminf L_N/N^(4/5)>=C_star, with no assertion about circuit upper bounds.

### 8.4 Matching prime candidate at every sufficiently large N

Let b_N be the least integer with b_N^5>=4N^3 and let k_N be the least
prime at least b_N. The prime number theorem consequence stated above
gives k_N/N^(3/5)->2^(2/5) along ALL integer N, not a subsequence.
Take the actual full payload m=N, ignored prefix p=0,
r_N=ceil(N/k_N), and t_N=ceil(sqrt(r_N)), using the exact integer ceiling.
These are among the selector's r0 choices. Then

    r_N/N^(2/5)->2^(-2/5),
    t_N/N^(1/5)->2^(-1/5).

Eventually r>=3,k>r+1,2<=t<r. The original sufficient threshold gives
d=O(N^(4/5))<m and B<=H<d, so these are admissible choices at EVERY
sufficiently large integer N. The same expansion above applies to these
explicit scaled parameters, using (8.3)--(8.4), and gives candidate loss
(C_star+o(1))N^(4/5). The global selected loss is no greater. Thus
limsup L_N/N^(4/5)<=C_star. Together with Section8.3 this proves(8.1),
and the already proved finite circuit inequality gives(8.2).

The evaluator has not changed: the prime number theorem is used only to
prove the existence of these comparison parameters inside its complete
finite search. Its runtime and uniformity remain the constructive bounds
of Sections1--2 and field.md. The sharp constant concerns this particular
certificate selector, whose arithmetic and loss are specified exactly.
