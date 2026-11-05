# sum-product-and-ordinary-condensation.md

Component references resolve through proof.md Section 0.

# Explicit finite approximate rings and four-coordinate expansion

Status: complete original quantitative proofs of the finite lemmas below,
including E_p(1/10000) from cz-conditional-entropy.md for every prime
p >= 2^20000. These proofs do not use CZ Lemma 7.3, Bourgain's filling
lemma, or an unspecified absolute constant. The surrounding condenser arguments
require their separate audits; this file claims no external review.

Attribution: the strategy is the Katz--Tao subset construction and the
self-improving approximate-ring construction of Terence Tao, *The
sum-product phenomenon in arbitrary rings*, Contributions to Discrete
Mathematics 4(2), 2009, pp. 59--82, especially Lemma 1.4, Proposition 3.1,
and Theorem 1.2. Official source:
https://cdm.ucalgary.ca/article/view/61994
The sumset tools below are Ruzsa's inequalities and Petridis's minimal-
growth argument. The explicit powers, normalization, proofs and application
below are independently supplied here. We work only in commutative rings
with identity. No inclusion AB+AC subset A(B+C) is used: that inclusion is
false in general. Scalar distributions and every product of sets below
are distinguished explicitly.

Notation: nA denotes an n-fold sumset, with 0A={0}; AA={ab:a,b in A};
aA is multiplication by the single element a. A unit is an invertible
ring element. All sets used in cardinality estimates are finite and
nonempty unless indicated otherwise. Products in F_p^d are coordinatewise.

## 1. Fully quantified additive tools

### 1.1 Ruzsa triangle and covering

For finite nonempty subsets X,Y,Z of an abelian group,

  |X-Z| |Y| <= |X-Y| |Y-Z|.                              (1)

For each w in X-Z choose one representation w=x_w-z_w. The map
(w,y) -> (x_w-y,y-z_w) is injective: summing the two coordinates recovers
w and the chosen x_w then recovers y. This proves (1).

If |X+Y| <= M|Y|, there is T subset X with |T|<=M and

  X subset Y-Y+T.                                        (2)

Choose a maximal collection of disjoint translates t+Y, t in X. Their
number is at most M. Every other x+Y intersects one of them, which gives
x in t+Y-Y. The same conclusion follows from |X-Y|<=M|Y| by using -Y.

### 1.2 Explicit iterated sumset estimate

If |U+V|<=K|U| for finite nonempty U,V and K>=1, then for all nonnegative
integers m,n with m+n>0,

  |mV-nV| <= K^(m+n) |U|.                                (3)

Here is a proof rather than an invocation of Plunnecke's theorem. Choose
nonempty X subset U minimizing k=|X+V|/|X|. Thus k<=K and |Y+V|>=k|Y|
for every Y subset X, including Y empty. For any finite set C, order its
elements and add the translates X+c successively. For the next c, let
Y={x in X:x+c belongs to the previous union X+C_old}. Then Y+V+c is
contained in the previous union X+V+C_old. Consequently the newly added
part of X+V+c has size at most |X+V|-|Y+V|<=k(|X|-|Y|). The latter
parenthesis is exactly the newly added size of X+c. Summing proves

  |X+V+C| <= k |X+C|.

Induction gives |X+jV|<=k^j|X| for every j>=0. Apply (1) to mV,-X,nV
to obtain |mV-nV|<=k^(m+n)|X|<=K^(m+n)|U|. This includes m=0 or n=0.

In particular, |A+A|<=K|A| implies

  |mA-nA| <= K^(m+n)|A|.                                 (4)

If P contains a set C of size N and |P-P|<=L N, then N<=|P|<=L N,
and (3) with U=P,V=-P yields

  |mP-nP| <= L^(m+n+1)N.                                 (5)

The hypotheses of (5) include C subset P. The bound |P|<=|P-P| follows
by fixing one element of P. Interchanging m and n in (3) gives the stated
signs. The harmless extra power in (5) avoids any unrecorded constants.

## 2. Explicit Katz--Tao subset lemma

Let C be a set of units in a commutative ring, M=|C|, K>=2, and suppose

  |C+C|, |CC| <= K M,
  |(C-C) outside the units| < K^(-25) M.                  (6)

Then there exists B subset C with

  |B| >= M/(2K),
  |BB-BB| <= K^45 |B|.                                   (7)

Proof. Write r(z)=#{(c,c') in C^2:cc'=z}. Cauchy--Schwarz gives
sum_z r(z)^2 >= M^4/|CC| >= M^3/K. Because multiplication by each unit
is injective, the same energy equals sum_{a,b in C}|aC intersect bC|.
Each intersection has size at most M. It follows that at least M^2/(2K)
ordered pairs satisfy |aC intersect bC|>=M/(2K). Fix b in C for which

  B={a in C:|aC intersect bC|>=M/(2K)}

has size at least M/(2K).

For a in B put H=aC intersect bC. The sets aC+H and H+bC have size
at most KM, since H subset aC and H subset bC. Applying (1) with middle
set -H gives

  |aC-bC| <= 2K^3 M.

By (2), for each a in B there is X_a with |X_a|<=2K^3<=K^4 such that

  aC subset b(C-C)+X_a.                                  (8)

For a,a' in B, first multiply the inclusion for a' by the single a,
then use commutativity and the inclusion for a twice. This gives

  aa'C subset b^2(2C-2C)+T_(a,a'),
  T_(a,a')=bX_a-bX_a+aX_(a'),
  |T_(a,a')|<=K^12.                                      (9)

For every d in BB-BB choose one representation d=aa'-a''a'''. Subtract
the corresponding inclusions (9). Since dC subset aa'C-a''a'''C, we get

  dC subset b^2(4C-4C)+T_d,  |T_d|<=K^24.                (10)

No injectivity of multiplication by d is needed here. Assign every
c in C to one translate in (10) containing dc. One translate receives
at least M/K^24 elements c. Fix c0 among them. For each such c,

  f=c-c0 in C-C,   df in b^2(8C-8C).

These f are distinct. By (6), more than
M/K^24-M/K^25 >= M/K^25 of them are units (K>=2). For each such f there
is a unique g in 8C-8C with df=b^2g because b is a unit. Conversely,
a pair (f,g) with f a unit determines at most one d. By (4), there are
at most

  |C-C| |8C-8C| <= K^2 M * K^16 M = K^18 M^2

possible pairs. Counting them for all d gives

  |BB-BB| <= K^43 M <= 2K^44 |B| <= K^45 |B|,

as asserted.

## 3. Explicit homogeneous approximate-ring theorem

Let B be a finite nonempty subset of a commutative ring with identity,
containing a unit a. Let N=|B| and L>=2. If

  |BB-BB| <= L N,
  |(B-B) outside the units| < L^(-36) N,                  (11)

then there is a finite unital subring S such that

  B subset aS,   |S| <= L^13 N.                          (12)

Proof. Normalize C=a^(-1)B. This set has size N and contains 1.
Commutativity gives CC-CC=a^(-2)(BB-BB); multiplication by a unit
also preserves which differences are units. Thus the assumptions (11)
hold for C. Put P=CC. Then C subset P, |P-P|<=LN, and (5) applies.
Also C-C subset P-P, so |C-C|<=LN.

Define

  S={x in R: |xC+C|<=L^36 N}.                            (13)

We prove a self-improvement:

  x in S implies |xC+C|<=L^10 N.                        (14)

The map (c1,c2)->xc1+c2 maps C^2 into at most L^36 N values. Its
collision count is at least N^3/L^36. Each collision gives
x(c1-c1')=c2'-c2. For any fixed difference d=c1-c1', there are at most
N^2 such collisions: there are at most N representations of d and at
most N representations of its determined image xd. Therefore at least
N/L^36 distinct d in C-C satisfy xd in C-C. By (11), one such d is a
unit. Fix this d. Because d in C-C and xd in C-C,

  dC subset P-P,   xdC subset P-P.

It follows that P-dC subset 2P-P. By (5), |P-dC|<=L^4N. Since |dC|=N,
covering (2) supplies T of size at most L^4 with

  P subset dC-dC+T.

Multiply this inclusion by the single x to get

  xP subset xdC-xdC+xT subset 2P-2P+xT.

Since C subset P,

  xC+C subset 3P-2P+xT.

By (5), its size is at most L^6N * L^4 = L^10N, proving (14).
This argument uses only inclusions with explicitly fixed multipliers.

Next, S is finite and obeys the size bound in (12). Repeating the collision
argument using (14) gives at least N/L^10 differences d in C-C with
xd in C-C for each x in S. Discarding nonunits leaves more than
N/L^10-N/L^36 >= N/L^11 choices. Each pair (d,xd) with d a unit uniquely
determines x. There are at most |C-C|^2<=L^2N^2 such pairs. Hence

  |S| <= L^13 N.                                        (15)

For each x in S, (14) and (2) give a set T_x with

  xC subset C-C+T_x,   |T_x|<=L^10.                      (16)

If x,y in S, then

  (x+y)C+C subset 3C-2C+T_x+T_y.

Because C subset P, (5) bounds |3C-2C| by L^6N. Thus this sumset has
size at most L^26N, so x+y in S. Also 0 in S. Every finite additive
semigroup containing zero is a group: for any x, a repetition in the
sequence 0,x,2x,... gives a positive t with tx=0, and (t-1)x is -x.
Thus S is an additive subgroup.

For multiplication, apply (16) for y, multiply by x, and apply (16) for x:

  xyC+C subset xC-xC+xT_y+C
           subset 3C-2C+T_x-T_x+xT_y.

The right side has size at most L^6N * L^30 = L^36N, hence xy in S.
Finally, for c in C, cC+C subset P+P, whose size is at most L^3N by (5).
Thus C subset S, in particular 1 in S. We have proved that S is a finite
unital subring with C subset S and (15). Multiplying by a gives (12).

## 4. Coordinate-fiber bound on bad differences

Let C subset F_p^d be nonempty, |C+C|<=K|C|, and let pi_j denote the
j-th coordinate projection. Put

  Z_j=(C-C) intersect ker(pi_j).

Choose one representative c from C for each value of pi_j(C). Their
translates c+Z_j are disjoint, so

  |pi_j(C)| |Z_j| <= |C+Z_j| <= |2C-C| <= K^3|C|.         (17)

A vector in F_p^d is a nonunit exactly when at least one coordinate is
zero. Taking a union bound in (17) proves

  |(C-C) outside the units|
      <= K^3 |C| sum_(j=1)^d 1/|pi_j(C)|.                (18)

This bounds the number of distinct bad differences, not just the number
of pairs yielding them. That distinction is needed in Sections 2 and 3.

## 5. Numerical four-coordinate sum-product theorem

Let p be a prime with p>=2^20000, set gamma=1/10000, and let A subset
F_p^4 be nonempty. Assume:

(a) |A|<p^(7/3-1/16);
(b) at least half of A consists of units;
(c) every D subset A with |D|>=p^(-1/16)|A| has both coordinate-pair
    projections pi_12(D), pi_34(D) of size at least p^(5/4);
(d) no such D is contained in either kind of plane

    {(x,y,c1*x,c2*y):x,y in F_p},
    {(x,y,c2*y,c1*x):x,y in F_p},

    for any c1,c2 in F_p.

Then

  max(|A+A|,|AA|) > p^gamma |A|.                          (19)

In particular, E_p(gamma) of cz-conditional-entropy.md holds with the
explicit gamma and prime threshold above.

Proof. Suppose instead both sets in (19) have size at most p^gamma|A|.
Set C=A intersect (F_p^*)^4, M=|C|, N=|A|, and K=2p^gamma. Then
M>=N/2 and |C+C|,|CC|<=KM. As p^(1/16)>=2, C qualifies for (c).
Every individual coordinate projection has size at least p^(1/4):
each paired projection has size at least p^(5/4), while its other
coordinate has at most p values. Equation (18) gives

  Z:=|(C-C) outside the units| <=4K^3 M p^(-1/4).          (20)

The following strict numerical inequality will handle all bad-difference
requirements at once:

  8K^1624 p^(-1/4)
    =2^1627 p^(-219/2500)
    <=2^(1627-1752)=2^(-125)<1.                          (21)

In particular 4K^28 p^(-1/4)<1, so (20) meets (6). Section 2 yields
B subset C with |B|>=M/(2K) and |BB-BB|<=K^45|B|. Since M<=2K|B|,
(20) and (21) imply

  |(B-B) outside the units|
    <=8K^4 |B| p^(-1/4) < K^(-1620)|B|.

Apply Section 3 with L=K^45 to obtain a unit a in B and a finite unital
subring S of F_p^4 such that

  B subset aS,   |S|<=K^585|B|.                          (22)

The subset B still qualifies for (c)--(d), since

  |B|>=N/(4K)=N/(8p^gamma)>=p^(-1/16)N.

For the final inequality, p^(1/16-gamma)>=2^1248>8 by the prime threshold.
Moreover 2<=p^(1/20000) gives K<=p^(3/20000), hence K^585<p^(1/8).
Equations (a) and (22) give

  |S|<p^(7/3-1/16+1/8)=p^(115/48)<p^3.                  (23)

Every additive subgroup of F_p^4 is an F_p-vector subspace. Thus |S| is
p to an integer power, and (23) implies dim(S)<=2. Scaling by the unit a
preserves every paired projection cardinality. By (c), both paired
projections of S have size at least p^(5/4), so their dimensions are at
least two. Consequently dim(S)=2 and both projections are bijections
onto F_p^2.

It remains to classify this elementary ring situation. Write S as the
graph of the unital F_p-linear ring isomorphism phi:F_p^2 -> F_p^2 given
by its two projections. The images under phi of (1,0) and (0,1) are
nonzero complementary idempotents. An idempotent in F_p^2 has each
coordinate in {0,1}. The only nonzero complementary pair whose members
are distinct from (1,1) is (1,0),(0,1), in either order. F_p-linearity
therefore makes phi either the identity or the coordinate swap.

In the identity case, writing a=(a1,a2,a3,a4),

  aS={(x,y,(a3/a1)x,(a4/a2)y):x,y in F_p}.

In the swap case,

  aS={(x,y,(a3/a2)y,(a4/a1)x):x,y in F_p}.

All denominators are nonzero because a is a unit. These are exactly the
two forbidden kinds of plane in (d), but B subset aS and B qualifies for
(d). This contradiction proves (19).

## 7. Integration and scope

The input used by the CZ conditional-entropy module is precisely the
four-coordinate premise in Section 5. Thus its choice
eta<=1/100000, 1000eta<=gamma, p^eta>=2^64 may now use

  gamma=1/10000,  eta=1/10000000.

This forces p>=2^640000000, much stronger than p>=2^20000, so there is
no additional large-prime caveat beyond that module's stated guard.
Choosing a smaller positive eta remains valid. These are deliberately
conservative numerical constants, not useful practical parameters.

The present file supplies the missing four-coordinate set-expansion premise.
It does not by itself prove an extractor,
condenser, advice-breaker, or nonmalleable code: distributional reductions
and parameter plumbing in the respective companion modules still need
to be checked. No finite numerical experiment is substituted for a
uniform theorem over primes.


# Ordinary low-rate somewhere condensation: finite scalar construction

Date: 2026-11-05 UTC. Status: **research proof, pending independent audit**.
This note gives a finite ordinary condenser, including scalar sum-product and
entropy growth. It uses the complete elementary graph/BSG and minimal-growth
proofs in cz-conditional-entropy.md Section 3, not that file's unproved 4D
expansion premise. Every finite padding, amplification, composition, selector,
mixture, full-bit output, and evaluator step is specified below. This is not a
completed circuit candidate, and the new proof needs independent checking.
It does not change or submit any frozen M=1000 file. No external/model review, installation, paid call,
credentials, or submission was used. Public papers were read only.

All logarithms are binary, TV is half-L1, and H_inf is min-entropy. Parameters
called constants are fixed independently of the outer input length. Statements
with real entropy thresholds mean literal atom bounds, so integer ceilings are
not silently omitted.

## 1. Useful result and precise boundary

Fix rational 0<beta<tau<1. One fixed multiple-source map at input rate beta/2
suffices for ALL amplification rounds. If that map has the quantitative
rectangle bound R below, choose its arity ell>=2 and rational gamma with

    0 < gamma <= min(beta/100, (1-tau)/100, 1/100),
    a = gamma/(16 ell),    g = a/4,
    T = ceil((tau-beta)/g),
    q = ceil(r/ell^T),     D = (ell+1)^T.                    (1)

For all r passing the explicit finite guards below, a deterministic D-row map
on r bits, with q bits per row, maps every beta*r source within

    e_B = T 2^(-gamma*q/4)                                  (2)

of a selector-somewhere source of entropy at least tau*q. For any 0<nu<tau,
it is within

    e_C = T 2^(-gamma*q/4) + D 2^(-nu*q)                    (3)

of a convex mixture whose selected fixed row has entropy (tau-nu)q on every
component. For ceil(63q/64) use tau=255/256,nu=1/256: rate127/128 exceeds
that threshold for q>=128.

Section 5 reduces R to a scalar entropy-growth assertion for AB+C in prime
fields. Sections 6--7 prove it with explicit constants using an elementary
prime-field pivot and the shared finite BSG proof. The construction takes
2^O(r) bit operations. Thus r=O(log outer N) is sufficient for polynomial
outer-N time; polynomial time in r itself is not claimed or needed here.

The gain over the earlier conditional outline is that no separate unknown
padding theorem, variable-rate extractor family, composition constant, prime
finding algorithm, full-bit output conversion, or small-input evaluator is
left in the dependency list. The scalar assertion G is proved below; it was
an isolated premise during the initial derivation. The separate four-dimensional expansion needed by the
CZ advice condenser is NOT proved here.

## 2. The finite rectangle premise R

For the fixed beta,tau, let sigma=(1+tau)/2. Assume integers ell>=2,m_R>=1, a
rational gamma obeying (1), and specified maps

    E_m : ({0,1}^m)^ell -> {0,1}^m,    m>=m_R,

such that for EVERY collection of sets A_1,...,A_ell of size at least
2^(beta*m/2), and EVERY B subset of {0,1}^m of size at most 2^(sigma*m),

 |{x in product A_i : E_m(x) in B}| / product |A_i|
       <= 2^(1-gamma*m).                                   (R)

Only this one-sided rectangle bound is used. In particular no strongness,
retained seed, conditional independence between output rows, or uniformity of
E_m on arbitrary correlated blocks is assumed. The marginal input blocks of
the condenser need NOT be independent. Independent uniform A_i variables
appear only in a counting contradiction.

A many-independent-source extractor producing m bits with error 2^(-gamma*m)
for entropy beta*m/2 implies R: the target set has uniform mass at most
2^(-(1-sigma)m), and gamma<=1-sigma. Conversely, R is all the elementary
condenser argument needs; claiming a full extractor is unnecessary here.

Let j0 be the least integer such that 2^j0>=ell+1. Define

    m_* = max(m_R, ceil(4/a), ceil(2(j0+1)/a), ceil(4/gamma)). (4)

All constants in (1)--(4) are explicit rational/integer expressions in the
scalar-growth constants, whose explicit choices are supplied in Section 7.

## 3. One finite boosting step, with arbitrary source probabilities

Suppose a source X on n bits has entropy rho*n, where

    beta <= rho < tau,    m=ceil(n/ell)>=m_*.

Append ell*m-n zero bits and split into ell ordered m-bit blocks. The padded
source has entropy at least u*ell*m, where u=rho-a/4. Indeed its rate is at
least rho-1/m>=rho-a/4, since ell*m-n<ell and rho<1. Output those ell blocks
and E_m of them, in that order; call this ell+1-row matrix M.

For each row i define its heavy set H_i by the threshold

    P(M_i=y) >= 2^(-(u+a)m).

Every H_i has at most 2^((u+a)m) elements. Put epsilon=2^(-gamma*m/4). We
first show that the event W on which all ell+1 rows are heavy has mass below
epsilon. Suppose its mass p is at least epsilon. Because every full padded
input atom is at most 2^(-u*ell*m), the number of input tuples in W is at least

    |W| >= p 2^(u*ell*m) >= 2^((u*ell-gamma/4)m).             (5)

Here W can be viewed as a subset of product H_i for i<=ell, whether or not
every tuple in that product lies in the original padded support. Therefore
for each i<=ell,

 |H_i| >= |W| / product_(j!=i,j<=ell)|H_j|
         >= 2^((u-gamma/4-(ell-1)a)m)
         >= 2^(beta*m/2).                                   (6)

For the last inequality, u>=beta-a/4 and
 gamma/4+(ell-1)a+a/4 <= gamma/4+ell*a+a/4 <= gamma;
then gamma<=beta/100 gives considerably more than beta/2. Also

    u+a = rho+3a/4 < tau+3a/4 <= sigma.                     (7)

Apply R to H_1,...,H_ell and B=H_(ell+1). Its upper bound on the fraction
mapped into H_(ell+1) is 2^(1-gamma*m). On the other hand, (5) and the upper
bounds on |H_i| make this fraction at least

    2^(-(gamma/4+ell*a)m) = 2^(-5gamma*m/16).                (8)

These are incompatible: m>=2/gamma implies
(11/16)gamma*m>=11/8>1. Hence P(W)<epsilon.

Condition on W not occurring. This changes the input distribution in TV by
P(W)<epsilon. Choose the first nonheavy row I. For any row value y, its
selected-row probability in this conditioned law is at most

    (ell+1) 2^(-(u+a)m)/(1-epsilon)
       <= 2^(j0+1) 2^(-(u+a)m).                             (9)

The second inequality uses gamma*m>=4, hence epsilon<=1/2. Since
m>=2(j0+1)/a, equation (9) gives selected entropy at least

    (u+a/2)m = (rho+a/4)m = (rho+g)m.                      (10)

Thus the boosting step has gain g and error epsilon. Notice that no reduction
to flat sources or integer-valued 2^(rho*n) support size is needed. For reference
the full guard is

    m_* = max(m_R, ceil(4/a), ceil(2(j0+1)/a), ceil(4/gamma)). (11)

## 4. Composition, mixing, and exact evaluator

### 4.1 Selector composition does not condition away entropy

Suppose a matrix law is e-close to a selector-somewhere k-source. This means
there is a nearby matrix Mbar and a jointly distributed selector I such that
Z=Mbar_I has H_inf(Z)>=k. It does not assert H_inf(Z|I=i)>=k. Apply a new
condenser rowwise. On the nearby matrix, the selected input Z has the promised
unconditional entropy, so the new condenser supplies a new selector J and a
nearby selected output with the next entropy threshold. A joint coupling can
alter just the chosen row to realize that output; all other rows can be kept.
The error adds. The composite selector is (I,J), and there is no union bound
over the number of rows. This argument iterates using product row indices.

For completeness, whenever a selected row Z is e-close to a target Z*, use
maximal coupling of Z and Z* and sample Z* by its conditional law given Z.
Change only the selected matrix row. Matrix disagreement probability is at
most e. Conversely a coupling of two nearby matrices lets one lift a selector
from either matrix to the common probability space without additional error.
These elementary couplings justify each use of the preceding paragraph.

### 4.2 Fixed schedule

For j=0,...,T, put rho_j=beta+j*g and n_j=ceil(r/ell^j). At round j<T,
rho_j<tau by the definition of T. Apply Section 3 to every current row,
appending zeros separately and using E_(n_(j+1)). The identity

    ceil(ceil(r/ell^j)/ell)=ceil(r/ell^(j+1))

follows by comparing both sides to the least integer v with ell^(j+1)*v>=r.
All row lengths are therefore exactly n_j; the final q=n_T and D in (1) are
literal, not asymptotic. The gains yield rho_T>=tau. Since every n_j>=q,
error is at most the sum of T copies of 2^(-gamma*q/4), proving (2).

It suffices to require

    r >= ell^T m_* .                                       (12)

This slightly wasteful integer guard ensures q>=m_* and hence all intermediate
guards. Furthermore r/ell^T<=q<r/ell^T+1, so q=Theta(r) with explicit constants.
Row order is lexicographic in {1,...,ell+1}^T, where 1,...,ell are the padded
input blocks and ell+1 is the E output. No selector, heavy set, coupling, or
mixture label is computed by the evaluator.

### 4.3 Selector to convex mixture

Let Z=Mbar_I have entropy tau*q. Whenever P(I=i)>=2^(-nu*q),

    P(Z=z | I=i) <= 2^(-tau*q)/P(I=i)
                  <= 2^(-(tau-nu)q).

The total mass of the other indices is at most D*2^(-nu*q). On those indices
replace the selected row by an independent uniform q-bit row, leaving a
well-defined matrix law. Its mixture indexed by I has a fixed high-entropy
row on EVERY component. This modification and the initial coupling prove (3).
If the bound exceeds one it can be capped at one; for a useful guarantee take
q large enough that it is below one. No assumption that synthetic rows remain
in the true evaluator's support is made.

To expose one exponential rate with no hidden polynomial prefactor, additionally
require

 q >= max( (8/gamma) ceil(log(2T)), (2/nu) ceil(log(2D)) ).    (13)

Then (3) is at most

    2^(-min(gamma/8,nu/2)*q)
      <= 2^(-c*r),    c=min(gamma/8,nu/2)/ell^T>0.           (14)

Ceilings can be applied to the right sides of (13). For beta>=tau an identity
row already has the selector target; no boosting is needed. For P3 alone,
beta>=4/5 likewise makes the identity sufficient.

### 4.4 Small lengths and running time

For r not passing (12)--(13), define the evaluator to output D all-zero q-bit
rows, with q as in (1). This is a total deterministic function but makes NO
entropy guarantee at those lengths. It is not a finite fallback for the full
circuit theorem, whose outer schedule still has to ensure the guarded regime.
Once E_m is evaluated in time 2^(C*m) for fixed C, the constant number of tree
nodes gives total time 2^O(r). The zero-padding and row enumeration cost O(Dq).
Thus no advice table or source-dependent search is hidden in this wrapper.

## 5. Reducing R to one quantitative scalar entropy-growth premise

The following finite interface is inspired by BIW Lemma 3.1. It is first used
as a premise to keep the construction separate from its proof. Sections 6--7
then prove it with explicit rational constants and a finite threshold.

**G(beta).** There are specified rational theta,kappa>0 and integer p0>=2 such
that for every prime p>=p0, every real h with

    (beta/4) log p <= h <= (9/10)log p,

and every three independent F_p sources A,B,C of min-entropy at least h,
AB+C is within p^(-kappa) of a distribution having min-entropy at least

    min((1+theta)h, (9/10)log p).                            (G)

This statement includes quantitative smoothing, not just support growth. It
only needs prime fields. Section 7 supplies theta=1/21500,
kappa=beta/172000, and p0=2^ceil(44720000/beta).
The earlier weaker values theta=1/1000000, kappa=beta/8000000 and
p0=2^ceil(512000000/beta) remain valid: lower output entropy, larger
error and a larger prime threshold only weaken this contract. Inner
breaker constructions may retain those values unchanged.

### 5.1 Canonical prime and arithmetic in exponential input time

For m>=16 choose the least prime p in [2^m,2^(m+1)]. A deterministic sieve up
to 2^(m+1) finds it using 2^O(m) bit operations and space. Arithmetic in F_p
uses integer reduction modulo p; canonical representatives are 0,...,p-1.
Embed an m-bit string as its binary integer, injectively.

Here is a short existence proof sufficient for these lengths, avoiding an
unproved prime-finding oracle. Let P(v) be the product of primes at most v.
Induction gives P(v)<4^v: for odd v=2n+1, primes in (n+1,2n+1] divide the
binomial coefficient (2n+1 choose n), which is less than 4^n because the two
central terms are equal and their sum is less than the whole binomial sum;
multiply by P(n+1)<4^(n+1). For even v>2 use P(v)=P(v-1), and check v<=3.
If no prime lies in (n,2n], prime valuations in (2n choose n) imply

 4^n/(2n+1) <= (2n choose n)
              <= (2n)^(sqrt(2n)) 4^(2n/3).                 (15)

Indeed each prime p<=sqrt(2n) contributes a power at most 2n, since each term
of its factorial-valuation sum is 0 or 1. Primes above sqrt(2n) have exponent
at most one. Those in (2n/3,n] have exponent zero; the remaining ones have
product at most P(floor(2n/3))<4^(2n/3). For n=2^m, m>=16, (15) would imply

 (2/3)2^m <= log(2n+1)+(m+1)2^((m+1)/2)
            < (m+2)+(m+1)2^((m+1)/2),

which is false at m=16 and remains false by induction: each of the two terms
on the right grows by a factor less than 2 when m increases by one, while the
left doubles. This proves the needed prime interval for all guarded m.

This deliberately uses the permitted 2^O(m) regime. It does not import AKS,
a prime-density algorithm, or an unproved uniform polynomial-in-m constructor.

### 5.2 A finite entropy-growth tree

Let J be the least nonnegative integer satisfying

    (1+theta)^J (beta/4) >= 9/10.                            (16)

It can be computed by rational multiplication; no floating-point logarithm is
needed. Define F_0(x)=x and F_(j+1)=F_j(x_left)F_j(x_mid)+F_j(x_right), on
three disjoint equal-sized groups of inputs, in that order, in F_p. Thus F_J
has 3^J input leaves. Each m-bit leaf has entropy at least beta*m/2, which is
at least (beta/4)log p because log p<=m+1<=2m. Repeatedly using G shows that
F_J is within

    e_J = ((3^J-1)/2) p^(-kappa)                            (17)

of an F_p source of entropy (9/10)log p. At each level independently replace
the three disjoint subtrees by their high-entropy approximations. Product TV
is at most the sum of their errors; a deterministic function contracts TV.
One more application of G adds p^(-kappa). The recurrence e_(j+1)<=3e_j+
p^(-kappa), e_0=0, proves (17) and preserves independence across disjoint
output groups. No assertion of conditional strongness is used.

### 5.3 Elementary full-length extraction from eight independent high-rate groups

If A,B are independent F_p sources of point masses at most 1/K_A,1/K_B, then
for every nontrivial additive character chi_t, Cauchy--Schwarz and character
orthogonality give

 |E chi_t(AB)|^2
 <= (sum_a P_A(a)^2) (sum_a |sum_b P_B(b)chi_t(ab)|^2)
 =  (sum_a P_A(a)^2) p(sum_b P_B(b)^2)
 <= p/(K_A K_B).                                           (18)

The equality follows by expanding the square: sum_a chi_t(a(b-b')) is p
when b=b' and zero otherwise; the latter is the finite geometric-series
identity. For four independent sources A_1,...,A_4 with entropy
(9/10)log p, put W=A_1 A_2+A_3 A_4. Equation (18) and independence show each
nontrivial Fourier bias of W has magnitude at most p^(-4/5). Applying the
same character orthogonality to P_W-1/p, followed by Cauchy--Schwarz, yields

    TV(W,U_(F_p)) <= (1/2)sqrt(p-1)p^(-4/5)
                  <= (1/2)p^(-3/10).                       (19)

Make two independent W values using eight disjoint F_J trees. Let
b_1=floor(m/2), b_2=ceil(m/2). Output the low b_1 bits of the first canonical
integer W and the low b_2 bits of the second, concatenated. For uniform
V in {0,...,p-1}, reduction modulo M=2^b has TV from uniform on M values
exactly s(M-s)/(Mp), where s is the remainder of p modulo M. This is at most
M/(4p). Consequently the two reductions together add at most 2^(-m/2).
Thus this fixed map on

    ell = 8*3^J                                             (20)

independent m-bit inputs outputs exactly m bits and has error at most

 8e_J + p^(-3/10) + 2^(-m/2)
    <= C 2^(-c0*m),
 c0=min(kappa,3/10,1/2),    C=4(3^J-1)+2.                  (21)

There is no lost seed to subtract from the output. All input groups used here
are independent because they are disjoint groups of the independent E inputs.
This does not claim the blocks of the final single-source condenser are
independent; only the R test distribution is a product distribution.

### 5.4 Concrete reduction of constants

Set gamma=min(c0/2,beta/100,(1-tau)/100,1/100). Let c_C be the least integer
with 2^c_C>=C, and c_p be the least integer with 2^c_p>=p0. Take

    m_R = max(16,c_p,ceil(2*c_C/c0)).                        (22)

Then (21) is at most 2^(-gamma*m) for every m>=m_R. Its uniform comparison
implies R as explained in Section 2. Equations (16),(20),(22),(11)--(14) now
specify every wrapper parameter, including total row count, row length, source
threshold, error, asymptotic coefficient, and finite guard, in terms of G's
three explicit constants from Section 7. Every tree choice and arithmetic operation is
fixed. All constants may be very large; no claim about practical feasibility
or final offline proof-size fit follows.

### 5.5 An optional smaller entropy cap

Suppose G_rho has G form with its own constants and both9/10 caps replaced
by rho. Fix integer k>=2 and rho>1/2+1/(2k). Choose J by
(1+theta)^J*(beta/4)>=rho. Tree induction gives(17) at entropy rho*log p. For W=sum_(i=1)^k A_(2i-1)A_(2i), (18) and
independence bound each Fourier bias by p^(k(1/2-rho)), hence
TV(W,U)<=p^(-d)/2, d=k(rho-1/2)-1/2>0.
Use two independent W values from4k disjoint F_J trees. Their low-bit
map outputs m bits with error 4k*e_J+p^(-d)+2^(-m/2).
Thus ell=4k*3^J, C=2k*(3^J-1)+2 and c0=min(kappa,d,1/2)
replace (20)--(21); (22) and Sections2--4 apply with this ell and the
same T,M,D0,gamma formulas. Inner calls keep k=2,cap9/10.
The outer G_rho is proved in Section7.3.

## 6. An elementary finite scalar sum-product theorem

The argument uses Garaev's prime-field pivot (R58, Lemma 3), specialized to
one set. We give its proof and a deliberately weak explicit bound. It does not
invoke Garaev's full unequal-size theorem or an unspecified implied constant.
The sumset tool used here is proved in cz-conditional-entropy.md Section 3.1:
if finite nonempty X,C in an abelian group satisfy |X+C|<=K|X|, there is a
nonempty X' subset of X with

    |X'+jC|<=K^j |X'|, hence |jC|<=K^j|X|, for every j>=1.  (23)

Indeed choose X' minimizing |Y+C|/|Y| over nonempty Y subset of X. Its
minimal-growth translate argument applies to every additional set, so induction
gives the displayed bound. This is the exact elementary lemma, not a call to
an unproved asymptotic Plunnecke estimate.

We also use the elementary difference triangle inequality

    |U-V||Z| <= |U+Z||V+Z|.                                (24)

For each d in U-V fix a representation u_d-v_d. The map
(d,z) -> (u_d+z,v_d+z) is injective, proving (24).

### 6.1 Prime-field pivot, including its counting proof

Let B be a nonzero subset of F_p with b=|B|>=2. Put Q=|BB|.
Multiplicative energy and Cauchy--Schwarz give

    sum_(x,y in B) |xB intersect yB| >= b^4/Q.

Choose x0,y0 in B so their intersection has size at least b^2/Q. Let
C subset of B be defined by x0*C=x0*B intersect y0*B, and s=x0/y0. Then
|C|>=b^2/Q and s*C subset of B. Commutativity and b>=2 imply Q<b^2,
so C has at least two elements.

Let R=(B-B)/(C-C), using nonzero denominators. If R is a proper subset of
F_p, then R+s is not R. Otherwise invariance under the nonzero translation s
would imply invariance under all its multiples, which fill the additive prime
field. Choose r in R with r+s outside R, and write
r=(x1-x2)/(y1-y2), with x1,x2 in B and distinct y1,y2 in C. The map

    (b',c) -> b'+(r+s)c,    (b',c) in B x C,

is injective: a collision with different c would put r+s in R (R=-R), while
equal c forces equal b'. Its image lies inside B+B+rB because sC subset B.
After multiplying by the nonzero y1-y2, we obtain a set of size at least b|C|
inside

    (x1-x2)B + (y1-y2)B + (y1-y2)B.                        (25)

If R=F_p, average over t in F_p the collision count of b'+tc on B x C.
The average is at most b|C|+(b|C|)^2/p: equal c contribute exactly b|C|,
and each ordered pair of different c has at most one t for each pair of b'
values. Thus some t has image size at least

    b|C|/(1+b|C|/p) >= (1/2)min(b|C|,p),

by Cauchy--Schwarz on representation counts. Represent t=(x1-x2)/(y1-y2)
using R, with x1,x2 in B and distinct y1,y2 in C. Scaling gives a subset of

    (x1-x2)B + (y1-y2)B                                   (26)

with that size. In either case four coefficients from B produce a sum of four
or six signed dilates of B whose size is at least

    (1/2)min(b^3/Q,p).                                     (27)

Zero differences cause no difficulty; the chosen denominator is nonzero.

### 6.2 Popular overlaps and an explicit bound

Let A subset of F_p^*, N=|A|>=2, and set

    K=max(|A+A|,|AA|)/N >=1.

Multiplicative energy gives sum_(a,b in A)|aA intersect bA|>=N^3/K.
Choose b0 in A with sum_(b in A)|bA intersect b0A|>=N^2/K, and let

    B={b in A : |bA intersect b0A| >= N/(2K)}.

The complementary contribution is at most N^2/(2K), while any one summand is
at most N. Consequently |B|>=N/(2K). For each b in B put I=bA intersect b0A
and C_b=bA union b0A. Since I lies in both dilates,

    |I+C_b| <= 2|A+A| <=2KN.

Apply (23) with j=2 and base I; this gives

    |bA+b0A| <= |2C_b| <= (2KN)^2/|I| <=8K^3 N.

Equation (24) gives |bA-b0A|<=2K^3N, so both signs are bounded by 8K^3N.
If |B|>=2, apply Section 6.1 to B. Its product set has size at most KN.
The resulting four- or six-term set therefore has size at least

    min(N^2/(16K^4),p/2).                                  (28)

Each term is one of x1*A, -x2*A, y1*A, -y2*A, with repeated terms allowed.
Let P be their union. The preceding signed bounds imply

    |b0A+P|<=32K^3N,    |6P|<=(32K^3)^6N=2^30 K^18 N.     (29)

The latter follows from (23) with base b0A. A four-term sum has cardinality
at most |6P| as well: add two fixed elements of nonempty P to embed a translate
of it into 6P. Hence (28)--(29) imply at least one of

    K^22>=N/2^34,      K^18>=p/(2^31 N).

Since K>=1, in all cases

    K^22 >= 2^(-34) min(N,p/N).                            (30)

If |B|<2 then K>N/4. For N>=2 this already implies K^22>=N/2^34, so (30)
also holds. Finally, when N<=p^(9/10), we have p/N>=N^(1/9). Thus

 max(|A+A|,|AA|) >= (1/4) N^(1+1/198).                     (31)

The weaker constant 1/4 uses 2^(-34/22)>1/4. This very coarse scalar theorem
has enough room for the entropy-growth losses below. All set arguments above
are finite and explicit. Prime order was used only in the translation step of
the pivot, which is why arbitrary extension fields cannot be substituted.

## 7. Finite scalar entropy growth from the shared BSG proof

Set

    eta=1/21500, a0=1/512,
    theta=eta,   kappa=beta*eta/8,
    P0=ceil(2080/(beta*eta)),   p0=2^P0.                    (32)

We prove G with these constants. Let p>=p0 be prime, and let independent
A,B,C have entropy at least h in the range stated in G. The prime guard gives h>=520/eta>2. Put
M=floor(2^h)-1. Then M>=2^(h-1), M<=p^(9/10), and

    h>=520/eta,     M^eta>=2^64, M^(eta/256)>4.                           (33)

The last bound uses eta*log M>=520-eta>512. Consequently
16<=M^(eta/128), 2^26<=M^(eta/16), and 2<=M^(eta/256).
Put c_s=1/128,c_p=1/16 and

  L1=4+a0+c_s, P1=23+8a0+c_p,
  L2=10+3a0+3c_s, P2=44+15a0+7c_s+c_p.

A distribution with atom cap 2^(-h)<=1/(M+1) is a convex mixture
of uniform laws on (M+1)-point sets. The extreme points of
0<=p_x<=1/(M+1), sum p_x=1 have every coordinate at its bounds:
two intermediate coordinates could be perturbed oppositely; if only one
remained, the sum and rational cap force it to an endpoint as well.
Apply this independently to the three sources. The zero-removal and
re-flattening step below reduces each flat triple to nonzero M-point
supports. Mixing preserves the final cap and total variation bound.

### 7.1 A finite collision-growth lemma at scale M

The elementary energy lemma in cz-conditional-entropy.md Section 3.3, written
at scale M instead of p^2, says: for sets U,V in an abelian group, with
sizes in [M^(1-d),M], if cp(U_U+U_V)>M^(-(1+e)), there are U',V' with

 |U'|,|V'|>=M^(1-2d-e)/16,
 |U'+V'|<=2^26 M^(1+7d+8e).                                (34)

This is a literal substitution in that proof's general finite formula, not
an extension of its hypotheses. Its multiplication version applies in F_p^*.
The full graph, popular-sum, and three-edge-walk proof is in that section and
contains no unproved expansion premise.

For every S subset of F_p^* with M^(1-11eta)<=|S|<=M, equation (31) and
M^(eta/256)>4 give

 max(|S+S|,|SS|)/|S|
    >= (1/4) M^((1-11eta)/198)
    >= M^((1-11eta)/198-eta/256)
    > M^(217eta/2).                                         (35)

The last strict inequality has exact margin
(1-11eta)/198-eta/256-217eta/2=223/181632000>0.
The prefactor1/4 is absorbed by M^(eta/256)>4.

We now show that if A_1,A_2,C_0 are nonzero M-point sets, the distribution U_A1 U_A2+U_C0 is M^(-eta)-close to a law
with collision probability at most M^(-(1+(2+a0)eta)). The following is the
scalar version of the finite peeling argument; its full parameters are given
so it does not import the CZ four-dimensional hypothesis.

Call D plus-good if cp(U_D+U_V)<=M^(-(1+(2+a0)eta)) for EVERY nonzero V whose
size is exactly M; define times-good using multiplication. Unions of
disjoint plus-good sets are plus-good because the squared L2 norm is convex.
Starting with R=A_1, repeatedly remove plus-good pieces until R is times-good
or |R|<=M^(-eta)|A_1|. If neither stopping condition holds, |R|>M^(1-eta).
A witness to its failure of times-goodness and (34), with d=eta,e=(2+a0)eta,
supply D subset R and B' with

 |D|,|B'|>=M^(1-L1*eta),    |DB'|<=M^(1+P1*eta).             (36)

The constants in (34) are charged using c_s,c_p above. If D were not
plus-good, (34) again, with d=L1*eta,e=(2+a0)eta, supplies D' subset D and C' with

 |D'|,|C'|>=M^(1-L2*eta),   |D'+C'|<=M^(1+P2*eta).           (37)

The self-sum inequality in the shared Section 3.1 gives

 |D'+D'|<=M^(1+(2P2+L2)*eta),
 |D'D'|<=|DD|<=|DB'|^2/|B'|<=M^(1+(2P1+L1)*eta).

Dividing their sum by |D'| bounds it by
M^(13869eta/128)+M^(7705eta/128)
 <=M^(27739eta/256)<M^(217eta/2).
Here L2=5135/512<11, 2P2+2L2=13869/128 and
2P1+L1+L2=7705/128. This contradicts (35), since
|D'|>=M^(1-L2*eta)>=M^(1-11eta).
Therefore D is plus-good and can be
removed. The analysis-only peeling terminates after finitely many points.

For its plus-good union, fix a2 in A_2 and write
U_D*a2+U_C0=a2*(U_D+a2^(-1)U_C0); scaling is bijective and preserves the
allowed test-set sizes. For a times-good remainder, fix c0 and use translation.
Both components have the stated collision bound; convexity preserves it after
mixing. If a remainder of relative mass at most M^(-eta) was discarded, replace
it by the good mixture. That mixture is nonempty because M^(-eta)<1. This
proves the claimed TV and collision bounds.

### 7.2 Smoothing, zero handling, and constants of G

A flat (M+1)-point source loses at most 1/(M+1) in TV when 0 is
removed. Each resulting support has M or M+1 nonzero points.
A uniform (M+1)-point law is the average of its uniform M-point
subsets (omit one point uniformly), so all three nonzero laws are
independent mixtures of exact M-point laws. Apply the preceding result
componentwise. Removing zero in all three sources costs at most3/M.

If a probability vector P has cp(P)<=M^(-(1+(2+a0)eta)), the total mass on points
with P(x)>M^(-(1+(1+a0)eta)) is at most

    cp(P) / M^(-(1+(1+a0)eta)) <= M^(-eta).

Move just the excess above this atom cap to coordinates below it. There is
enough total capacity because M^(1+(1+a0)eta)<=p: M<=p^.9 and
.9(1+(1+a0)eta)<1. This produces a law of entropy (1+(1+a0)eta)log M at TV cost at
most M^(-eta). Its entropy is at least (1+eta)h: a0*eta*h>=65/64,
while 1+(1+a0)eta<65/64.
Consequently AB+C is within

    3/M + M^(-eta) + M^(-eta) <=5M^(-eta)
      <=5*2^eta*p^(-beta*eta/4) <=p^(-beta*eta/8)            (38)

of a law of entropy at least (1+theta)h, hence at least the smaller target
in G. The last inequality follows from log p>=2080/(beta*eta), whose spare
exponent beta*eta*log p/8 is at least260, larger than log(5)+eta.
Averaging over the independent flat decompositions preserves the same cap
and error. This proves G with (32).

### 7.3 Restricted outer scalar growth

Put k_f=1280,eta_f=5/2326,sigma=1/1120 and
rho_f=1/2+1/(2k_f)+eta_f/(8k_f)=11918429/23818240. In(8), if
M^(1-d)<=|U|<=M,|V|=M,L=M^(1+e), then alpha=1/(2M^e) gives

 |U'|>=M^(1-d-e)/16, |V'|>=M^(1-e)/16,
 |U'+V'|<=2^26 M^(1+d+8e).                                 (39)

Likewise multiplicatively. Use eta=eta_f,a0=c_s=c_p=2^(-20).
Use remainder cutoff M^(-sigma*eta) and collision threshold M^(-(1+(1+sigma+a0)eta)). Thus first
left deficit d=sigma*eta and collision exponent e=(1+sigma+a0)eta.
Both witnesses have size M. Formula(39) twice gives left deficits
L1*eta,L2*eta, right deficit E*eta, growth P1*eta,P2*eta, where

 L1=1+2sigma+a0+c_s, L2=L1+1+sigma+a0+c_s,
 E=1+sigma+a0+c_s, P1=sigma+8(1+sigma+a0)+c_p,
 P2=L1+8(1+sigma+a0)+c_p.

The self-sum inequalities give normalized expansion at most
M^((2P2+E+L2)eta)+M^((2P1+E+L2)eta)<=M^((771490807/36700160)eta),
charging2 by c_s*eta. For N<=p^rho_f, (30) gives expansion
at least(1/4)N^s_f,s_f=1081801/23836858. For
M^(1-L2*eta)<=|S|<=M the lower bound exceeds this upper bound:

 s_f*(1-L2*eta_f)-c_s*eta_f-(771490807/36700160)eta_f
   =4201250227/25435289810108416>0.

Take P0=ceil(2^28/(beta*sigma*eta_f)),h<=rho_f*log p.
Then h>=2^26/(sigma*eta_f). Hence c_s*eta_f*log M>26, absorbing16,2^26,4 and2 by their stated charges. Repeat normalization
and peeling; scaling/mixing preserve it.
Smoothing now discards at most M^(-sigma*eta): divide the collision
bound by the unchanged atom cap M^(-(1+(1+a0)eta)). Capacity and entropy
floor are unchanged, since rho_f*(1+(1+a0)eta_f)<1 and
 a0*eta_f*h>=65/64. Zero removal costs3/M, so total error is
5M^(-sigma*eta_f)<=p^(-beta*sigma*eta_f/8), by the new prime guard.
Thus G_rho_f has theta=eta_f,kappa=beta*sigma*eta_f/8,p0=2^P0.
Fourier decay eta_f/8>=kappa gives c0=kappa in Section5.5.
Section7.2 and inner calls stay unchanged.

The ordinary primitive uses only the finite proofs here and shared Section3. 

## 8. Sources, scope, and proof-size accounting

Primary sources actually inspected:

- Barak, Kindler, Shaltiel, Sudakov, Wigderson, *Simulating Independence: New
  Constructions of Condensers, Ramsey Graphs, Dispersers, and Extractors*, full
  version ECCC TR10-037, Theorem 3.6, Section 5, Remark 4.9:
  https://eccc.weizmann.ac.il/report/2010/037/download
- Barak, Impagliazzo, Wigderson, *Extracting Randomness Using Few Independent
  Sources*, Sections 3.1--3.4 and Appendix A:
  https://www.math.ias.edu/~avi/PUBLICATIONS/MYPAPERS/BIW04/biw.pdf
- M. Z. Garaev, *A Quantified Version of Bourgain's Sum-Product Estimate in F_p
  for Subsets of Incomparable Sizes*, Electronic Journal of Combinatorics
  15 (2008), R58, especially Lemma 3 and its proof:
  https://www.combinatorics.org/ojs/index.php/eljc/article/download/v15i1r58/pdf/

The heavy-set/counting architecture is attributed to BKSSW. The AB+C growth
tree and additive/multiplicative peeling architecture are attributed to BIW.
The prime-field pivot is attributed to Garaev. Only the proved pivot is used.
The finite BSG and minimal-growth tools are proved in the companion module,
with attribution there to the classical BSG/Petridis methods. The coarse
constants, finite schedule, padding, smoothing, prime sieve, and Fourier ending
here are independently written. No claim of a new extractor theorem or
ownership of the cited constructions is made; this is original proof support.

Partitions and flat labels are proof devices. The evaluator uses fixed
arithmetic, padding, AB+C trees, field-product sums and bit truncation;
source-linear extraction is not replaced by seed enumeration.

