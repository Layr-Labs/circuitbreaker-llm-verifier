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

### 1.3 A sharper form of (5)

If C subset P, |C|=N and |P-P|<=LN, then, for k=m+n>=1,

    |mP-nP| <= L^(m+n) N.                                (5s)

Indeed put Q=|P|>=N and k0=|P-P|/Q<=LN/Q. Applying (3) with
U=P,V=-P, and interchanging the signs if necessary, gives
|mP-nP|<=k0^k Q<=(LN)^k/Q^(k-1)<=L^k N. All sets are nonempty;
k0>=1 since a translate of P lies in P-P. This improves (5) without
changing its hypotheses. The older bound (5) remains a valid corollary.

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

## 2a. The same subset lemma with its constants retained

Under the exact hypotheses (6), the subset B in Section 2 can be chosen with

    |B|>=M/(2K),   |BB-BB|<=256 K^37 |B|.                (7s)

Here is the full change to that proof's count. Its energy argument gives
B={a:|aC intersect bC|>=M/(2K)} for a fixed unit b in C. The Ruzsa triangle
bound is |aC-bC|<=2K^3 M. Covering therefore gives the same inclusion (8)
with |X_a|<=2K^3, keeping the factor 2 instead of replacing it by K.
For a,a' in B, the explicit set bX_a-bX_a+aX_a' has size at most8K^9.
Subtracting two such inclusions, for each represented d in BB-BB,

    dC subset b^2(4C-4C)+T_d,   |T_d|<=64K^18.

Assign each c in C to one translate containing dc. One translate receives
at least M/(64K^18) distinct c. Fix c0 in that translate. The differences
f=c-c0 are distinct elements of C-C and satisfy df in b^2(8C-8C).
The bad differences number less than M/K^25<=M/(128K^18), because
K>=2 implies K^7>=128. Thus strictly more than M/(128K^18) of these f
are units. Each such f gives one g in8C-8C with df=b^2g; a unit f and g
uniquely determine d. By (4), at mostK^18M^2 pairs(f,g) are available.
Counting for all d yields |BB-BB|<=128K^36M<=256K^37|B|.
Every inclusion uses a fixed multiplier, not an invalid distributive
set inclusion. Since256<=K^8, the old K^45 bound is retained as well.

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

## 3c. A finite ring from unit actions

Let B be a nonempty set of UNITS,|B|=N,L>=2,and suppose
 |BB-BB|<=LN, |(B-B) outside units|<L^(-7)N.
Then there is a finite unital subring T with
 B subset aT, |T|<=L^7N, for any fixeda in B.             (12a)
Normalize C=a^(-1)B. Then1 in C,everyc in C is a unit,P=CC contains C,
|P-P|<=LN,and normalization preserves nonunitdifferences. Define
S={x:|xC+C|<=L^7N}. Collisioncounting gives at leastN^3/L^7 collisions
xc1+c2=xc1'+c2'. For fixed differenced=c1-c1',at mostN^2 representations
are possible. Hence at leastN/L^7 distinctd in C-C satisfyxd in C-C.
The strict bad-difference bound leaves a UNITd. ThusdC andxdC are both
subsets ofP-P. This fixedunit acts bijectively,and distributivity gives
 |xC+C|=|d(xC+C)|=|xdC+dC|<=|2P-2P|<=L^4N.              (14a)
This uses no distributivity identity for products of sets.

Repeating collisioncounting with(14a) leaves more than
N/L^4-N/L^7>=N/L^5 unitchoices,as1-L^(-3)>=1/L forL>=2.
Each pair(d,xd)withd a unit determinesx,andthere are at most
|C-C|^2<=L^2N^2 such pairs. Selectingonepair already provesS finite;
doublecounting then gives|S|<=L^7N. Covering(2) gives,for eachx in S,
 xC subsetC-C+T_x, |T_x|<=L^4.

We need two actions,not arbitrary closure properties of S.
First,the unitd chosen above proves
 |(x+1)C+C|<=|xdC+2dC|<=|3P-3P|<=L^6N<=L^7N.
ThusS+1 subsetS. Second,for eachfixedc in C,
 cxC+C subsetcC-cC+C+cT_x subset2P-P+cT_x,
whose size is at mostL^3N L^4=L^7N by(5s). ThereforecS subsetS.
Becausec is a unit andS is finite,both these actions restrict to
permutations ofS;their inverseactions preserveS as well.

For any monomialm in elements ofC(includingemptyproduct1),multiplication
M_m and its inverse preserveS. Conjugating translationby1 gives
 M_m(x->x+1)M_m^(-1)=(x->x+m).
ConsequentlyS is invariant under translations by all such monomials,
their negatives and every finite integerlinearcombination. Since0 in S,
the additivegroupT generated by these monomials is contained in S.
It contains1 andC,and is closed under multiplication because products
of monomials are monomials and finite sums distribute. ThusT is a finite
unital subring,|T|<=|S|<=L^7N,andB=aC subsetaT,proving(12a).
Every inverse multiplier used above exists because all ofC consists of units.

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

## 5d. A stronger rate from unit actions

Under exactly Section5's hypotheses(a)--(d),forp>=2^20000000,
 max(|A+A|,|AA|)>p^(1/1054)|A|.                           (19u)
Put gamma=1/1054 and assume failure. WriteN=|A|,C=A intersectunits,
M=|C|,K=2p^gamma. ThenM>=N/2,|C+C|,|CC|<=KM. Sincep^(1/16)>2,
C qualifies for(c),so each individualcoordinateprojection has at least
p^(1/4)values. Equation(18) bounds its bad differences by
 Z<=4K^3M p^(-1/4).
The decisive strict inequality is
 2^59 K^263 p^(-1/4)=2^322 p^(-1/2108)
 <=2^(322-9487)<1.                                      (21u)
Here1/4-263/1054=1/2108 and2108*9487=19998596<20000000.
This also implies4K^28p^(-1/4)<1,so the hypotheses(6)hold.
Section2a givesB subsetC,|B|>=M/(2K),|BB-BB|<=L|B|,L=256K^37.
Its bad differences number at most8K^4|B|p^(-1/4),sinceM<=2K|B|.
BecauseL^7=2^56K^259,(21u)makes this strictly less thanL^(-7)|B|.
Every element ofB is a unit,sinceB subsetC. Section3c therefore
provides a unita in B andfiniteunitalringT with
 B subsetaT, |T|<=L^7|B|.

The subsetB qualifies for(c)--(d):|B|>=N/(8p^gamma)>=p^(-1/16)N,
sincep^(1/16-gamma)>8 at the displayedprime threshold. Also
 L=2^45p^(37/1054)<=p^(37/1054+45/20000000)<p^(1/28).
ThusL^7<p^(1/4),and|T|<p^(7/3-1/16+1/4)=p^(121/48)<p^3.
Every additive subgroup ofF_p^4 is anF_p-vector subspace,so dimT<=2.
Unit scaling preserves both pairedprojectioncardinalities. By(c),each
projection ofT has size at leastp^(5/4),hence dimensionat leasttwo.
ConsequentlydimT=2 and both projections bijectontoF_p^2.
WriteT as the graph of the resulting unitalF_p-linear ringisomorphism
phi:F_p^2->F_p^2. The images of(1,0),(0,1)are nonzero complementary
idempotents,so are that pair in eitherorder. F_p-linearity makesphi
identity orcoordinateswap. Writinga=(a1,a2,a3,a4),allai nonzero,
aT is therefore one of
 {(x,y,(a3/a1)x,(a4/a2)y):x,y inF_p},
 {(x,y,(a3/a2)y,(a4/a1)x):x,y inF_p}.
These are precisely the forbiddenplane types in(d). SinceB qualifies
andB subsetaT,this is a contradiction,proving(19u)with allsubsetquantifiers.

The selectedfamily useseta=1/707235,671eta<gamma=1/1054,
because671*1054=707234<707235,andtheta=eta/11=1/7779585.
The CZ guardeta*u>=64 givesu>=45263040>20000000,so itsprimep>2^u
meets this theorem's threshold. All other CZ/uniform-source guards retain
their complete supplied proofs. No prior verdict or numericalexperiment
is a premise of this theorem.

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


# Scalar sum-product and one support-growth step

From the original ordinary-condenser module only the scalar sum-product
estimate (30) of Section 6 is used here, by the new Lemmas RG and S1' of
Section 5; its somewhere condenser, the entropy-growth statement G and G's
proof are not part of this package. All
logarithms are binary, TV is half-L1, and H_inf is min-entropy. Statements
with real entropy thresholds mean literal atom bounds, so integer ceilings are
not silently omitted.

## 5. One support-growth step

This section replaces the original reduction to scalar entropy growth. It
uses only (30) of Section 6.2 (for A subset F_p^*, N=|A|>=2 and
K=max(|A+A|,|AA|)/N: K^22>=2^(-34)min(N,p/N)); the sumset inequality (6)
and the energy lemma (8) of restricted-nonmalleability-and-breaking.md
Sections 3.1 and 3.3, cited below as (R6) and (R8) (not this file's (6) and
(8)), both stated in any abelian group and used here in (F_p,+) and
(F_p^*,x); and prime-fields Sections 1-2, recalled in Section 5.2. U_D is
uniform on D and cp is collision probability.

Constants: eta_RG=1/996, xi=1/5977, d=eta_RG, e=2eta_RG, r=(1-2xi)/(1+2xi)=5975/5979.

    r(1-d-2e)-22(3d+21e)=2015/5955084>0,   eta_RG/2-3xi=1/11906184>0.       (RG0)

Exactly: r(1-5eta_RG)-990eta_RG=[5975*991/5979-990]/996=2015/(5979*996), and
1/1992-3/5977=(5977-5976)/(1992*5977).

### 5.1 Lemma RG

**Lemma RG.** Let p be prime and M an integer with log M>=4871000 and
p>M^(2/(1+2xi)). Let A1,A2,C0 be subsets of F_p^* of size exactly M and
X1,X2,X3 independent and uniform on them. Then for every W subset F_p,

    P[X1*X2+X3 in W] <= M^(-d)+(|W|M^(-(1+e)))^(1/2).                 (RG)

Proof. Call D subset F_p^* plus-good if cp(U_D+U_V)<=M^(-(1+e)) for every
V subset F_p^* with |V|=M, and times-good if cp(U_D*U_V)<=M^(-(1+e)) for
every such V. A disjoint union of plus-good sets is plus-good: the cp of a
mixture is at most the largest cp, by convexity of the squared L2 norm.
Start with R=A1 and repeat: stop if R is times-good or |R|<=M^(1-d);
otherwise find a nonempty plus-good D subset R (below) and move D out of R.
Each round removes a point, so the procedure stops; it is a proof device,
never evaluated.

Finding D. R is not times-good: some V with |V|=M has cp(U_R*U_V)>1/L0,
L0=M^(1+e). Apply (R8) in F_p^*, with its L equal to L0, m=|R| in
(M^(1-d),M], n=M and alpha=M/(2L0)=1/(2M^e): there are D subset R and B'
subset V with

    |D|>=|R|/(16M^e)>=M^(1-d-e)/16,   |B'|>=M^(1-e)/16,
    |DB'|<=2^18(2L0)^3/(alpha^5|R|M)<=2^26 M^(1+d+8e).

Suppose D is not plus-good: some W' with |W'|=M has cp(U_D+U_W')>1/L0.
Apply (R8) in (F_p,+) with m=|D|, n=M and the same alpha: there are D'
subset D and C' subset W' with |D'|>=|D|/(16M^e), |C'|>=M^(1-e)/16 and
|D'+C'|<=2^26 M^(2+8e)/|D|. By (R6), additively and then multiplicatively
(using D' subset D),

    |D'+D'|/|D'| <= |D'+C'|^2/(|C'||D'|) <= 2^60 M^(3+18e)/|D|^3
                 <= 2^72 M^(3d+21e),
    |D'D'|/|D'| <= |DD|/|D'| <= |DB'|^2/(|B'||D'|)
                <= 2^60 M^(1+2d+18e)/|D| <= 2^64 M^(3d+19e),

so K(D')<=2^73 M^(3d+21e) (2^72+2^64<2^73). Also N'=|D'|>=M^(1-d-2e)/256
>=2 (as (1-5eta_RG)log M>=9), D' subset F_p^*, N'<=M and p/N'>=p/M>M^r (from
p>M^(2/(1+2xi))); as r<=1, min(N',p/N')>=M^(r(1-d-2e))/256. By (30),
2^(-42)M^(r(1-d-2e))<=K(D')^22<2^1606 M^(22(3d+21e)), that is,
[r(1-d-2e)-22(3d+21e)]log M<1648. By (RG0) the bracket is 2015/5955084, and
2015*4871000/5955084>1648 (the exact threshold is
1648*5955084/2015=4870460.76...). Contradiction: D is plus-good.

Distribution. Let G be the plus-good union and R the final remainder. X1 is
uniform on G with probability |G|/M and on R with probability |R|/M. Given
X1 in G and X2=a2: X1*a2+X3=a2*(X1+a2^(-1)X3), where a2^(-1)C0 is a subset
of F_p^* of size M, so this law has cp<=M^(-(1+e)) (scaling is a
bijection), and P[in W]<=(|W|cp)^(1/2) by Cauchy-Schwarz; average over a2.
Given X1 in R with R times-good, fix X3=c0: cp(U_R*U_A2)<=M^(-(1+e)) gives
the same bound; average over c0. If instead |R|<=M^(1-d), the R branch has
probability at most M^(-d). This proves (RG).

### 5.2 The one-step map and its support lemma

For an integer m>=1 let p_m be the least prime above 2^m; prime-fields
Sections 1-2 (restricted-nonmalleability-and-breaking.md, first part) prove
2^m<p_m<2^(m+1) and find p_m by trial division in O(2^(3m/2+2)(m+1)^2) bit
operations. Read an m-bit string as a binary integer in [0,2^m), hence as an
element of F_(p_m) (injective, since 2^m<p_m). Define

    E_m(x1,x2,x3) = (x1*x2+x3 mod p_m) mod 2^m,                        (S2)

the low m bits of the canonical representative in {0,...,p_m-1}. For alpha
in {0,1}^(3m) write alpha=(alpha1,alpha2,alpha3) with m-bit blocks; its four
ROWS are alpha1, alpha2, alpha3 and E_m(alpha).

**Lemma S1' (support growth).** Let m>=300000000. For every Q subset
{0,1}^(3m) with |Q|>=2^(3m/2-15), the union U of the four row images of Q
(a subset of {0,1}^m) satisfies |U|>2^((1/2+xi)m).

Proof. Suppose N=|U|<=2^((1/2+xi)m). As Q is a subset of U^3,
N>=2^(m/2-5), so M:=|U\{0}|>=N/2>=2^(m/2-6) and log M>=m/2-6>=149999994
>=4871000. Also M<=N<p_m^(1/2+xi) (as 2^m<p_m), so p_m>M^(2/(1+2xi)). Let
X1,X2,X3 be independent and uniform on U, viewed in F_(p_m). Lower bound: on
Q the fourth row lies in U, so P[E_m(X) in U]>=|Q|/N^3>=2^(-15-3xi*m).
Upper bound: E_m(X) lies in U exactly when Y=X1X2+X3 lies in
U*={y<p_m: y mod 2^m in U}, and |U*|<=2N since p_m<2^(m+1). P[some Xi=0]
<=3/N, and conditioned on no zero the Xi are independent and uniform on
U\{0}, of size M. Lemma RG with |W|=|U*|<=2N<=4M gives

    P[E_m(X) in U] <= 3/N+M^(-eta_RG)+2M^(-eta_RG) <= 4M^(-eta_RG),

since 3/N<=3/M<=M^(-eta_RG) (M^(1-eta_RG)>=3). So 4M^(-eta_RG)>=2^(-15-3xi*m), that
is, eta_RG*log M<=17+3xi*m, and with log M>=m/2-6: (eta_RG/2-3xi)m<=17+6eta_RG. But
(eta_RG/2-3xi)m=m/11906184>=300000000/11906184>25>18. Contradiction.

Only set counting and independent uniform variables on U appear; no
independence of the rows of an actual source is used. sampler-and-
disperser.md Sections 3-4 use Lemma S1' at m=P3*B0*L>=2990*309*1024.

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
is used once, by Lemma RG in Section 5. All set arguments above
are finite and explicit. Prime order was used only in the translation step of
the pivot, which is why arbitrary extension fields cannot be substituted.

## 8. Sources

Barak, Impagliazzo, Wigderson, *Extracting Randomness Using Few Independent
Sources*, Sections 3.1--3.4 and Appendix A:
https://www.math.ias.edu/~avi/PUBLICATIONS/MYPAPERS/BIW04/biw.pdf
(the AB+C growth step and the additive/multiplicative peeling architecture);
M. Z. Garaev, *A Quantified Version of Bourgain's Sum-Product Estimate in F_p
for Subsets of Incomparable Sizes*, Electronic Journal of Combinatorics 15
(2008), R58, Lemma 3 (the prime-field pivot, of which only the pivot is used).
The finite BSG and minimal-growth tools are proved in the companion module,
with attribution there. The exact-size witness form of (R8) used in Lemma RG
(|V|=M exactly, alpha=1/(2M^e)) and the use of (30) near rate 1/2 follow
mpjunior92's submission f7447a6f-a996-404e-bb7f-df63d9488837 (its
sum-product Section 7.3 and (39)), kept in its revisions
04e91ff4-9abd-474d-810a-1e48d1c4db6c and 7d6289bb-6624-4b01-ae0c-cb7dd033f0e0.
The coarse constants, Lemma RG and Lemma S1' here are independently written;
no new extractor theorem is claimed. All peeling partitions are proof devices;
the evaluator uses only fixed integer arithmetic, AB+C and bit truncation.

