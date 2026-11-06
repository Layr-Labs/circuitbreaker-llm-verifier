# restricted-nonmalleability-and-breaking.md

Component references resolve through proof.md Section 0.

# Canonical prime selection and exact field-to-bit encodings

This original finite proof removes a setup/encoding issue, not the sum-product
or condenser theorem. The prime-existence argument is the central-binomial
method of Erdos for the Bertrand--Chebyshev theorem, independently written
with a finite threshold and a short explicit small-case certificate.
Primary source: P. Erdos, *Beweis eines Satzes von Tschebyschef*, Acta Litt.
Sci. Szeged 5, 194--198 (1932):
https://combinatorica.hu/~p_erdos/1932-01.pdf .

All logarithms are base two. Integers are represented in binary. The actual
outer input length is N; the local block length u will be O(log N).

## 1. A completely finite prime interval proof

For real x>=1, let P(x) be the product of all primes at most x. We first
prove P(m)<4^m for every positive integer m, by strong induction.
For m=1,2 this is direct. If m=2r>=4 then m is composite and
P(m)=P(m-1)<4^(m-1)<4^m. If m=2r+1>=3, each prime in (r+1,2r+1]
divides binomial(2r+1,r). The two equal middle coefficients in that odd
binomial row show binomial(2r+1,r)<=2^(2r)=4^r. Thus

  P(2r+1)<=P(r+1)*binomial(2r+1,r)<4^(r+1)*4^r=4^(2r+1).

Consequently P(x)<4^x for every real x>=1, using its integer floor.

We prove that for every integer n>=2 there is a prime p with n<p<2n.
Assume no such prime exists. For n>=512, consider B=binomial(2n,n).
The maximum binomial coefficient is B. The 2n-1 interior coefficients
sum to at most (2n-1)B, and the two endpoints sum to 2<=B. Hence

  B>=4^n/(2n).                                               (1)

For any prime p, its valuation in B is

  sum_(j>=1) [floor(2n/p^j)-2floor(n/p^j)].

Each term is zero or one. For p<=sqrt(2n), the entire prime-power
contribution is therefore at most 2n; there are at most sqrt(2n) such
primes. For p>sqrt(2n), its exponent is at most one. A prime in
(2n/3,n] has exponent zero, since floor(2n/p)=2 and floor(n/p)=1.
There is by assumption no prime in (n,2n), and 2n is composite. Thus

  B <= (2n)^sqrt(2n) * P(2n/3)
    < (2n)^sqrt(2n) * 4^(2n/3).                             (2)

Combining (1)--(2) gives 4^(n/3)<(2n)^(sqrt(2n)+1).
This is impossible for n>=512. To verify the comparison explicitly, put
x=sqrt(2n)>=32. The reverse strict inequality is equivalent to

  x/log_2 x > 6*(1+1/x).

At x=32 it reads 32/5>99/16, equivalently 512>495. The left side is
strictly increasing for x>=32 (differentiate x/ln x); the right side is
strictly decreasing. This proves the contradiction for every n>=512.

For 2<=n<512 use the following increasing prime chain:

  2, 3, 5, 7, 13, 23, 43, 83, 163, 317, 631.

Each term after the first is strictly smaller than twice its predecessor.
All entries are prime: trial division by primes at most sqrt(631)<26,
namely 2,3,5,7,11,13,17,19,23, verifies the nontrivial entries. Given n in
this range, take the first chain entry p greater than n; its predecessor
q is at most n, so p<2q<=2n. This completes the proof for all n>=2.

No asymptotic prime distribution theorem or probabilistic primality test
is needed.

## 2. The exact prime algorithm and bit cost

For integer u>=1, set p(u) to the least prime strictly greater than 2^u.
Enumerate successive integers, and for each q test division by integers
2<=d<=floor(sqrt(q)). All comparisons may be performed as d*d<=q without
floating point. A candidate is prime exactly when none of those integers
divides it: every composite has a factor at most its square root.

Section 1 proves that this search terminates with

  2^u<p(u)<2^(u+1),    ceil(log p(u))=u+1.                    (3)

There are fewer than 2^u candidate integers, each with at most
2^((u+1)/2) trial divisors and O((u+1)^2) bit cost per division. A safe
bound is O(2^(3u/2+2)*(u+1)^2). If u<=c log N for one fixed c, this is
polynomial in the actual outer length N. The claim is NOT polynomial in u.
A primality table, uncomputed parameter or randomized setup is unnecessary.

Use the standard representatives 0,...,p-1 for F_p. Addition,
multiplication and inversion by Euclid have polynomial bit cost in u;
all encodings below are fixed leading-zero binary encodings of those
representatives. The actual nonlinear polynomial map and its degree must
be separately specified by its construction proof.

## 3. Source-rate loss from the prime interval

Embed an arbitrary u-bit source X into F_p by its integer value. This is
injective because p>2^u, so its min-entropy is unchanged. Suppose

  H_infinity(X)>=(1-gamma)*u,
  0<=gamma<delta<1.

It meets a theorem requiring entropy at least (1-delta)log p whenever

  (delta-gamma)*u >= 1-delta.                               (4)

Indeed log p<u+1 and the right side of (4) rearranges to
(1-gamma)u>=(1-delta)(u+1). A convenient conservative choice is
gamma<=delta/2 and u>=2/delta. The error p^-c is at most 2^(-cu),
whereas lower bounds involving p^c must use p>2^u in the correct direction.
Source partitions and conditioning may change the entropy rate first;
condition (4) must use the POST-partition entropy, not the original rate.

The bit length of one output field element is u+1, not u. A pair in F_p^2
therefore uses exactly

  nField=2*(u+1)                                              (5)

bits under the fixed injective encoding. This encoding preserves
min-entropy and TV distance exactly on the image. Any polynomial output
law is regarded as a distribution on the full nField-bit cube supported
on that image. No density or full-support property is inferred.



# Finite restricted CZ conditional entropy: one isolated deep premise

Status: complete finite reduction with the explicit expansion input E below.
The selected ring proof supplies E at gamma256/20481,p>=2^3200000.
The six-source evaluator is a restricted CZ variant; its fixed-point-free
condition is on the first four coordinates,as needed for appended advice.
Encoding rigidity,graph/energy estimates,fiber partition,smoothing and all
finite parameter inequalities are supplied here.

Attribution: the encoding and original eight-input polynomial map are due to Eshan
Chattopadhyay and David Zuckerman, *Non-Malleable Codes Against Constant
Split-State Tampering*, ECCC TR14-102, Sections 6--7:
https://eccc.weizmann.ac.il/report/2014/102/download/ . Their expansion theorem
extends Jean Bourgain's work over the product of two prime fields. The graph
argument below is a self-contained quantitative Balog--Szemeredi--Gowers
argument; the minimal-growth argument is Petridis's short Plunnecke method.
The constants below are conservative original bookkeeping, not quoted CZ
constants. All logarithms are base two, TV is half L1, and all products of
vectors are coordinatewise. U_A means uniform on a finite nonempty set A.

## 1. The exact remaining external premise

Fix a prime p>3 and a positive real gamma. E_p(gamma) is this statement:

For every A subset of F_p^4 with |A|<p^(7/3-1/16), assume
(i) at least half of A has all four coordinates nonzero; and
(ii) every D subset of A with |D|>=p^(-1/16)|A| has both coordinate-pair
projections of size at least p^(5/4), and is contained in neither of the
following kinds of two-dimensional subspace, for any c1,c2 in F_p:

  {(x,y,c1*x,c2*y): x,y in F_p},
  {(x,y,c2*y,c1*x): x,y in F_p}.

Then |A+A|+|A*A|>p^gamma |A|.

For the selectedgamma256/20481 and everyprimep>=2^3200000,
sum-product-and-ordinary-condensation.md Section5g proves this exact input
statement,with its full subset/projection quantifiers. The local theorem
below is stated conditionally to expose that input; the final selected
construction supplies it through the explicit ring proof. No external
Bourgain/CZ theorem or unknown finite constant is used as a premise.

## 2. Main conditional theorem and evaluator

Assume E_p(gamma), and choose a positive rational eta satisfying

  eta <= 1/7000, (92+2^-20)*eta <= gamma, eta*logp >= 64*2^20.       (1)

Let X1,...,X6 be independent F_p-valued sources with

  H_inf(Xi) >= log ceil(p^(1-eta)).                          (2)

Let fi be arbitrary maps on the respective supports into F_p. Assume at
least one of f1,...,f4 is fixed-point-free on its honest support. In particular, the
case in which every coordinate has distinct appended advice qualifies.
Define q(x)=x^4+x^2+x, e(x)=(x,q(x)), and

  E6(x1,...,x6) = (e(x1)+e(x2))*(e(x3)+e(x4))
                 +e(x5)+e(x6).               (3)

Put W=E6(X1,...,X6), V=E6(f1(X1),...,f6(X6)). Their joint law is within

  zeta = min(1,45*p^(-eta))                                 (4)

of a law (W*,V*) on F_p^2 x F_p^2 such that, for EVERY v with positive
V*-probability,

  H_inf(W* | V*=v) >= kappa = (19/6)*eta*log p.                   (5)

This is the joint smoothing needed for (N1) of condenser-integration.md.
No generic fixed-point/same? theorem, seeded nonmalleable primitive, or
coordinate-wise conditional claim about the original unsmoothed law is used.

Evaluator (3) takes a constant number of field operations. All proof
partitions and mixture labels are analysis objects and are never computed
by the evaluator. Section 8 gives the exact prime, bit-encoding and selected
weak-AGHP matching interface.

## 3. Two elementary additive-combinatorial tools

### 3.1 Small cross-sumset gives small self-sumset

In any abelian group, finite nonempty A,B satisfy

  |A+A| <= |A+B|^2/|B|.                                    (6)

Proof. Choose a nonempty X subset of B minimizing K=|X+A|/|X|.
Thus K<=|B+A|/|B|, and |Y+A|>=K|Y| for every Y subset of X,
including the empty set. For any finite C, add its translates X+c one at
a time. If Y subset of X consists of points whose translate by the next c
was already covered by X+C_old, then Y+A+c was already covered by
X+A+C_old. The newly added part of X+A+c has size at most
|X+A|-|Y+A|<=K(|X|-|Y|), namely K times the newly added part of X+c.
Summing proves |X+A+C|<=K|X+C|. Taking C=A gives
|X+A+A|<=K^2|X|. A translate of A+A lies inside this set, so
|A+A|<=K^2|X|<=|A+B|^2/|B|. This proves (6).
The same proof works multiplicatively in any abelian group.

### 3.2 A dense graph produces many three-edge paths

Let a bipartite graph on sets A,B of sizes m,n have at least alpha*m*n
edges, where 0<alpha<=1. There are A' subset of A, B' subset of B with

  |A'|>=alpha*m/8,  |B'|>=alpha*n/8,

such that every (a,b) in A' x B' is joined by at least

  alpha^5*m*n/2^18                                          (7)

three-edge walks a,b1,a1,b. Vertices may repeat.

Proof. Repeatedly remove a left vertex of degree below alpha*n/4 or a
right vertex of degree below alpha*m/4. In total this removes fewer than
alpha*m*n/2 edges. Call the remaining graph H, retaining the original
part-size parameters m,n. Set beta=alpha^3/2048. Call an ordered pair of
left vertices bad if it has fewer than beta*n common neighbors in H.
The sum, over right vertices b, of the number of bad ordered pairs in
N_H(b) is at most beta*m^2*n.

Choose b with probability proportional to its H-degree. Its degree is
below alpha*m/4 with probability at most 1/2. The expectation of
(number of bad pairs in N_H(b))/deg_H(b) is at most 2*beta*m/alpha.
Markov's inequality shows the ratio exceeds 8*beta*m/alpha with
probability at most 1/4. Thus some b0 has degree at least alpha*m/4
and ratio at most 8*beta*m/alpha. Put A0=N_H(b0). Delete from A0 every
a having more than 16*beta*m/alpha bad partners in A0; at most half of
A0 is deleted. Call the result A', of size at least alpha*m/8.

Every remaining left vertex has H-degree at least alpha*n/4. Therefore
the set B' of right vertices having at least alpha*|A'|/8 neighbors in
A' has size at least alpha*n/8: otherwise its incident edges, even
counting all |A'| possible neighbors per vertex, cannot supply the
alpha*n*|A'|/4 total edges after the other vertices contribute at most
alpha*n*|A'|/8.

For a in A', b in B', at least alpha^2*m/64 vertices a1 in A' are
neighbors of b. At most 16*beta*m/alpha=alpha^2*m/128 of these are bad
partners of a. Each of the remaining alpha^2*m/128 choices has at least
beta*n common neighbors b1 with a. Multiplication gives (7).

### 3.3 Explicit energy-to-small-sumset lemma

Let A,B be finite nonempty sets in an abelian group, m=|A|, n=|B|.
If cp(U_A+U_B)>1/L, there exist nonempty A',B' such that, with
alpha=max(m,n)/(2L),

  |A'|>=alpha*m/8, |B'|>=alpha*n/8,
  |A'+B'| <= 2^18*(2L)^3/(alpha^5*m*n).                    (8)

Indeed write r(s)=#{(a,b):a+b=s}. Let S contain sums with
r(s)>=mn/(2L). The contribution of other sums to sum r(s)^2 is at most
m^2*n^2/(2L). Since r(s)<=min(m,n), the graph of edges whose sum lies
in S has more than alpha*m*n edges. Also |S|<=2L. Apply 3.2.
Each three-edge walk yields a+b=s1-s2+s3 with si in S. For a fixed
(a,b), the triple determines b1=s1-a and a1=s2-b1, so different walks
have different triples. Distinct sums a+b use disjoint triples.
Counting gives (8). The assumption itself implies alpha<=1/2 because
cp(U_A+U_B)<=1/max(m,n).

In particular, if |A|,|B| lie in [p^(2-d),p^2] and
cp(U_A+U_B)>p^(-(2+e)), then

  |A'|,|B'| >= p^(2-2d-e)/16,
  |A'+B'| <= 2^26*p^(2+7d+8e).                            (9)

The same lemma applies to multiplication inside G=(F_p^*)^4. The proofs
in this section are complete and have no sum-product assumption.

## 4. Encoding rigidity supplies the expansion hypotheses

Suppose f,g have fibers of size at most L=p^(20eta), and f has no fixed
point on its honest domain. Let A be any subset of

  {(e(x)+e(y),e(f(x))+e(g(y))): x,y in F_p},                 (10)

where only the given domains are used. Every S subset of A with
|S|>=p^(2-1000eta) satisfies the hypotheses of E_p(gamma).
Consequently it expands by the factor p^gamma under sum/product.

Here and below p^eta>=2^64 implies L>16 and absorbs all displayed small
constant factors. The detailed verification is as follows.

First, any prescribed value of e(x)+e(y) has at most four preimages.
The first coordinate gives y=s-x; substituting in q(x)+q(s-x) gives
a degree-four polynomial with leading coefficient 2, nonzero since
p>3. A fixed second coordinate-pair value in (10) has at most 4L^2
preimages. Thus for any D subset of A,

  |projection_12(D)| >= |D|/4,
  |projection_34(D)| >= |D|/(4L^2).                       (11)

The number of distinct points of A with a zero coordinate is at most
5p(1+L)<=10pL. For the four coordinates, fixing x leaves at most
1,4,L,4L possible y values respectively. Since
10pL/|S|<=10p^(-1+1020eta)<1/2, S has the required nonzero fraction.

For the first forbidden plane, fix y=a and write g(a)=b. Its equations
imply f(x)=c1*x+d, d=c1*a-b, and

  q(c1*x+d)+q(b)-c2*(q(x)+q(a))=0.                         (12)

If c1=0 there are at most L possible x. If c1!=0 and (12) is a nonzero
polynomial, there are at most four. If it is identically zero, its cubic
coefficient forces d=0; its fourth, second and first coefficients give
c1^4=c1^2=c1=c2, hence c1=1. The first equation would then force
f(x)=x, so no honest x is permitted. This bounds the plane intersection
by pL. For the swapped plane the first equation gives
f(x)=c2*(q(x)+q(a))-b. If c2=0 there are at most L solutions. Otherwise
substitution into the second equation gives a degree-sixteen polynomial
with nonzero leading coefficient c2^4, so there are at most 16<=L
solutions. This intersection also has size at most pL.

For |D|>=p^(-1/16)|S| we have |D|>=p^(31/16-1000eta)>pL, so neither
plane can contain D. From (11), its smaller projection has size at least
p^(31/16-1041eta)>=p^(5/4). Finally |S|<=p^2<p^(7/3-1/16).
Ateta<=1/7000,the explicit margins are1-1021eta>=5979/7000>0,
31/16-1020eta-1>=2217/2800>0,and31/16-1041eta-5/4>=7543/14000>0.
Thus10p^(-1+1020eta)<=10p^-eta<1/2 andallstrict/weak
projection/plane comparisons above hold. This proves the
claim. It uses E_p(gamma) only for the FINAL expansion conclusion, not
for the preceding polynomial rigidity and counting statements.

## 5. Finite distributional energy growth

Let A1,A2,C be subsets of G=(F_p^*)^4, all of sizes in [p^(2-7eta),p^2].
Assume every S subset of A1 of size at least p^(2-1000eta) satisfies
|S+S|+|S*S|>p^gamma |S|. Then every W subset of F_p^4 with |W|<=T,
T=4p^(2+(16/3)eta), has P[U_A1*U_A2+U_C in W]<=2p^(-eta).

Call D plus-good if P[U_D+U_C0 in W]<=p^(-eta) for every C0 subset of G of
size in [p^(2-7eta),p^2] and every W with |W|<=T; define times-good
similarly. A disjoint union of plus-good sets is plus-good: for fixed C0,W
its hitting probability is the size-weighted average of the pieces'.

Start with R=A1 and peel off disjoint plus-good subsets, until R is
times-good or |R|<=p^(-eta)|A1|. If R is still large, |R|>=p^(2-8eta). If
it is not times-good, a witnessing B and W give P[U_R*U_B in W]>alpha,
alpha=p^(-eta): the graph on R x B with edges {(x,b):xb in W} has more than
alpha|R||B| edges. Section 3.2 gives D subset R, B' subset B with
|D|>=alpha|R|/8, |B'|>=alpha|B|/8, each pair in D x B' joined by at least
alpha^5|R||B|/2^18 walks x,b1,x1,b. In G a walk gives xb=s1*s2^(-1)*s3 with
s1=xb1, s2=x1b1, s3=x1b in W; for fixed (x,b) the triple fixes b1 and x1,
and the triple fixes the product. Counting triples,

  |D|>=p^(2-9eta)/8, |B'|>=p^(2-8eta)/8, |D*B'|<=2^24 p^(2+36eta).      (13)

Here 2^18 T^3/(alpha^5|R||B|)<=2^24 p^(6+21eta-4+15eta). We prove that D is
plus-good. Otherwise a witnessing C0 and W, with the same counting
additively (x+c=s1-s2+s3), give D' subset D, C' subset C0 with

  |D'|>=p^(2-10eta)/64, |C'|>=p^(2-8eta)/8, |D'+C'|<=2^27 p^(2+37eta).  (14)

Using (6) in the additive and multiplicative groups gives

  |D'+D'| <= |D'+C'|^2/|C'| <= 2^57 p^(2+82eta),
  |D'*D'| <= |D*D| <= |D*B'|^2/|B'| <= 2^51 p^(2+80eta).

After division by |D'|, their sum is at most
(2^63+2^57)p^(92eta)<2^64 p^(92eta)<=p^((92+2^-20)eta)<=p^gamma,
sinceeta*logp>=64*2^20. This explicitly pays the fixed factor2^64.
But |D'|>=p^(2-10eta)/64>=p^(2-1000eta), contradicting its expansion.
Thus D is plus-good. The peeling removes at least one point each time
and therefore terminates; no algorithmic runtime is asserted for this
analysis-only decomposition.

Let A+ be the union of peeled pieces and Ax the remainder. On the A+
component fix a2 in A2: P[U_A+*a2+U_C in W]=P[U_A+ + a2^(-1)*U_C in
a2^(-1)W]<=p^(-eta), as multiplication by a2 is a bijection, a2^(-1)C is an
allowed C0 and |a2^(-1)W|=|W|. If Ax is times-good, fix c in C:
P[U_Ax*U_A2 in W-c]<=p^(-eta). Otherwise |Ax|<=p^(-eta)|A1| and the Ax
component has probability at most p^(-eta). Averaging proves the claim.

## 6. Small fibers for the restricted six-source map

Suppose X1,...,X6 are independent flats withsupportsizes>=p^(1-2eta),
everyfi hasfibers<=p^(20eta),andoneoff1,...,f4 isfixed-point-free.
Then(W,V)forE6 iswithin36p^-eta ofa law withjointmin-entropy
 (2+(19/3)eta)logp.                                         (15)
Swappingwithineitherfirstpairorswappingthetwofactors preservesthemap,
so thefixed-point-freecoordinate can be putfirst. No symmetrybetween
theadditivepairandmultiplicativepairs isasserted.

Foreachpair(1,2),(3,4),(5,6),partitionthejointimage
(e(Xi)+e(Xj),e(fi(Xi))+e(fj(Xj))) bymultiplicityr=1,...,4.
Honestpairpreimages are<=4bySection4,so eachpart isuniform.
Discardpartswithimagesize<p^(2-6eta). Theirprobabilityperpair is
atmost(1+2+3+4)p^(2-6eta)/p^(2-4eta)=10p^-2eta.
Renormalizeallthreeindependentmixtures,TVcost<=30p^-2eta.
Eachsurvivingimage hasatmost10p^(1+20eta) zero-coordinatepoints;
theirrelativefraction<=10p^(-1+26eta)<=p^-eta. Deletethem at
totalcost<=3p^-eta. EveryremainingunitsetAi has sizein
[p^(2-7eta),p^2],using1-p^-eta>=p^-eta.

Thejointoutput is nowU_A1*U_A2+U_A3 DIRECTLY. Hereditaryexpansion
forA1followsfromSection4andthefirstcoordinate's fixed-point-freeness.
Section5 applies with C=A3: every target of at mostT=4p^(2+(16/3)eta)
points has probability<=eps_hit=2p^-eta. Every surviving mixture component
obeys this bound,so the mixture does too. No extra pair or multiplier is used.
Since eta<=1/7000 and p^eta>=2^64, 2<=T=4p^(2+(16/3)eta)<=p^4.
Choose W0 as floor(T) heaviest points. Since floor(T)>=T/2,
every point outsideW0 has mass<=2eps_hit/T=p^(-(2+(19/3)eta))=:tcap.
Clipping to tcap removes mass<=P(W0)<=eps_hit. Ambient p^4 has
at least1/tcap points;redistribute the removed mass belowtcap.
This costs at mosteps_hit inTV. Totalcost<=30p^-2eta+5p^-eta
<=35p^-eta<36p^-eta,proving the same supplied joint minentropy.
This proves(15)directlyforthesix-source map,notviaanunsupported
extensionoftheoriginaleight-source theorem.

## 7. Arbitrary fibers, arbitrary sources, and conditional smoothing

First take independent flat Xi with support Si of size at least
ceil(p^(1-eta)). Partition each Si into

  Si,0 = {x in Si: |Si intersect fi^(-1)(fi(x))|<=p^(20eta)},
  Si,1 = Si \ Si,0.                                       (16)

These parts are unions of COMPLETE fibers. The argument of fi^(-1) is
fi(x), not x; this corrects the displayed fiber-partition typo in CZ.
Discard any part with size below p^(-eta)|Si|. There is at most one
such part per coordinate because p^(-eta)<1/2. The total probability
of a branch using any discarded part is at most 6*p^(-eta).
Every retained product branch has independent flat coordinates of
support size at least p^(1-2eta).

On the all-small-fiber branch, Section 6 gives a joint law Q with
largest atom p^(-(2+(19/3)eta)), at TV cost at most 36*p^(-eta).
The Q-mass of values v having Q_V(v)<p^(-(2+(19/6)eta)) is at most p^(-(19/6)*eta),
since V has at most p^2 values. On every other v,
max_w Q(w|v)<=p^(-(19/6)*eta). Replace the conditional W law on bad v by
uniform on F_p^2, keeping Q_V. This costs at most p^(-(19/6)*eta), and gives
(5) on every v. Total branch cost is at most 37*p^(-eta).

Ifaretainedbranch hasalarge-fibercoordinatej,fixfj(Xj)andallother
honestcoordinates. V isfixed,andXj isuniformonacompletefiberofsize
>p^(20eta). Ifj<=4,thefirst-factor/pairsymmetries putitfirst.
Thehonestfirstcoordinate hascoefficientX3+X4inX1,so isinjective
exceptateventX3+X4=0,ofprobability<=1/|S4|<=p^(-1+2eta)<=p^-eta.
ReplaceWtherebyuniformonF_p^2,keepingV. Ifj=5or6,thehonestfirst
coordinate is(X1+X2)(X3+X4)+X5+X6,withcoefficientONEinXj,
so isinjectivewithoutanyexception. Thus everyfinebranch hasconditional
Wmin-entropy>=20eta logp,which implies(5),atworstTVcostp^-eta.
Foreachv,mixingthesefinebranchespreservestheconditionalmaximum-atom
bound. Thefirst-fourfixed-point-freehypothesis isretainedoneverybranch;
nofixed-point-freehypothesisoncoordinates5or6isneeded.

On each discarded product branch replace W by uniform independently of
its actual V. The error is at most its probability. Convex mixtures of
joint laws satisfying max_w P(W=w|V=v)<=p^(-(19/6)*eta) still satisfy that
inequality: multiply by P(V=v), sum the inequalities, then divide.
Adding the discarded-branch cost to the worst retained-branch cost
givesatmost(6+37)p^-eta<=45p^-eta,withoutpayingforall2^6branches separately.

Finally, a distribution with largest atom at most 1/M, M an integer,
is a convex mixture of uniform distributions on M-point sets. For
completeness, in the finite polytope 0<=p_x<=1/M and sum p_x=1, an
extreme point cannot have two coordinates strictly between 0 and 1/M,
since those coordinates can be perturbed in opposite directions. Nor
can it have exactly one, because the sum is an integer multiple of
1/M. Thus every extreme point has exactly M nonzero atoms, all 1/M.
Apply this with M=ceil(p^(1-eta)) independently to each Xi. Every
component preserves the fixed-point-free hypothesis on its support.
The proved uniform error and conditional maximum-atom bound survive
mixing, establishing the theorem under (2).

## 8. Prime choice, exact encoding and the selected weak-AGHP match

The canonical least prime above2^u obeys2^u<p<2^(u+1),by the supplied
finite prime proof. Trialdivision costsO(2^(3u/2+2)(u+1)^2),polynomial
in actualN wheneveru=O(logN). Encode the pairW in2(u+1)bits.
For the integer flat-source threshold ceil(p^(1-eta)),the safe guard
 k>=(1-eta)(u+1)+1                                    (17)
implies(2),since logceil(p^(1-eta))<=(1-eta)logp+1.
When k>=(1-rho)u,it suffices that(eta-rho)u>=2-eta; inparticular
rho<=eta/2 andeta*u>=4 suffice. Distinct fixed appendedadvices ensure
first-four fixed-point-freeness on actual honest supports; no global
extension offsupport is required. The extra ring guardu>=3200000 is
separate from the energy constant-absorption guardeta*logp>=64*2^20,which every selected call supplies viaeta*u>=64*2^20. This stronger guard also impliesp^eta>=2^64.

The selected conditional entropy is nowkappa=(19/6)eta logp. We do NOT
assert the former eight-input bilinear matching bound at this weaker
entropy target. The selectedseven-input weak-AGHP match is proved
explicitly in affine-and-standard-wrappers.md Section3.4,using only
the firstsix CZ inputs plusone independentrawr-bit seed andthe supplied
AGHP-specific W2 collision proof. Its public deficit/field/seed/prime
andoutput guards,whole seed retention,partnerkernel andcost2zeta are
all stated there. The earlierunused ordinarybilinear specialization
has been removed;the generic bilinear primitive itself remains proved.

# Conditional tools and the selected antichain multicopy breaker

Alternating extraction and the local A/Z/B operations are existing methods,
adapted from Cohen's local correlation breakers and the whole-tuple MC-LA
proof of publishedc78568c1. All required proofs are supplied below; no external
theorem or source verdict is a premise. The selected schedule and known-A
ledger are separately proved in8.2. The old one-hot schedule is not used.

All logarithms are binary and TV is half L1. WriteD(V|T)=TV((V,T),(U,T)),
where U is independent uniform on V's length. Define
 Havg_inf(X|T)=-log(sum_t Pr[T=t] max_x Pr[X=x|T=t]).

## 1. Primitive interfaces

Fix positive integers n,d,h,a,s,s0,q and0<eps<1 with2^-q<=eps,s<=h.
Supply these deterministic strong seeded extractors,all error at mosteps:
  I : {0,1}^n x {0,1}^s0 -> {0,1}^h, entropy threshold KI;
  F : {0,1}^n x {0,1}^a  -> {0,1}^h, entropy threshold KF;
  G : {0,1}^d x {0,1}^s  -> {0,1}^a, entropy threshold KG;
  J : {0,1}^h x {0,1}^a  -> {0,1}^s, entropy threshold KJ.

Strongness retains the complete seed. All pads retain ignored suffixes.
The actual selected dimensions and entropy guards are stated in8.2.

## 2. Local alternating operations

On a current h-bit stateW and d-bit original seedY,compute
 T=prefix_s W,A=G(Y,T),Z=J(W,A),B=G(Y,Z).
The selected round usesC=A orC=B according to the fixed public label code,
thenWnew=F(X,C). Initialization isI(X,prefix_s0Y). Section8.2 gives the
complete fixed schedule and proof for the selected multicopy contract.

## 3. Elementary conditional tools, with proofs

Let L=(X,X') and R=(Y,Y'). A rectangular history is a sequence of messages,
each a deterministic function of one whole side and previous messages.
Initially the sides are independent. On conditioning on a left message,
the left factor in the conditional product law changes and the right factor
does not; the symmetric statement holds for a right message. Induction
proves that L and R are independent conditional on every history value.
This works for arbitrary correlated tuples inside each side.

If the total left-message length is bL, then

  Havg_inf(X|H)>=H_inf(X)-bL.                               (5)

For a left revelation T with at most 2^b values, summing the largest joint
atom over its fibers gives

  sum_t Pr[t] max_x Pr[x|t] <= 2^b max_x Pr[x].

This proof also holds conditional on the preceding history. A right
revelation does not change the conditional X distribution once the preceding
history is fixed, by conditional independence. Thus only left lengths
contribute in the induction. Likewise right lengths bR imply
Havg_inf(Y|H)>=H_inf(Y)-bR. A message already determined by earlier messages
has zero additional cost.

Average entropy at least K+q implies conditional entropy at least K except
on history mass at most 2^-q<=eps: apply Markov to the largest conditional
atom. These bounds apply also to functions of a side instead of X or Y.

Here is the robust extraction rule used repeatedly. Suppose conditionally
on a rectangular history H that P is a p-bit function of the right side,
D(P|H)<=eta, and X is on the left with Havg_inf(X|H)>=K+q. If E is a strong
(K,eps) extractor, then

  D(E(X,P) | H,R) <= eta+2eps.                              (6)

For fixed H=t, write g_t(p)=TV(E(X|H=t,p),U), which lies in [0,1]. Conditional
independence implies that the left side of (6) averages g_t(P) even when the
ENTIRE right side is retained. Replacing P's marginal by uniform changes
this expectation by at most TV(P|t,U). On the good entropy histories the
uniform-seed expectation is at most eps; the other histories cost at most
eps altogether. Average to obtain (6). This proof requires no pointwise
claim about a fixed seed. The symmetric rule retains the entire left side
when extracting from the right.

If V is right-sided and Q is left-sided conditional on H, then

  D(V|H,Q)=D(V|H).                                         (7)

This follows by averaging conditional TV over H: V and Q are independent on
each branch. This allows an already nearly uniform seed to survive extra
revelations on the other side. Finally, if joint laws P,Q on (V,T) have TV
at most eta, the triangle inequality gives

  D_P(V|T) <= D_Q(V|T)+2eta.                                (8)

The second eta accounts for the possibly different T marginals.

## 8. Selected antichain multicopy instance

The selected NB inner map uses the five-round theorem8.2 below and the
explicit rounded-rational PXT instance of wrappers4.2. The earlier local
A,Z,B,F operations and robust lemmas are fully supplied in Sections1--3.

## 8.2 Five-round antichain with two copies per other label

The local A,Z,B,F operations are those of Section2. The selected
schedule below is different from the earlier one-hot schedule: its complete
proof is supplied here. The whole-tuple look-ahead method is adapted from
published c78568c1-1aa4-4da9-8556-877fb4218f6f; the antichain schedule,
separated-set invariant and the bounds below are this construction's changes. The five-round precursor is our submitted2f83b733-9c02-4d12-96cb-eb0f18326b01; its theorem is fully supplied here,not an accepted-result premise. This later variant additionally exploits NB's two-copies-per-other-label property through an explicit known-A look-ahead proof.
No source verdict is a mathematical premise.

Fix t>=1,c=t+1. Assign the eight fixed labels1,...,8 respectively the codes
 11000,10100,10010,01100,01010,00110,10001,01001.           (AC)
These are eight distinct weight-two binary words. For any two different
labels alpha,beta, some coordinate has AC(alpha)=1,AC(beta)=0: otherwise
the two sets of ones would be nested, hence equal since both have size2.

Let I,F,G,J have the dimensions and strong error eps of Section1, with
q>=lg(1/eps). Require
 h>=c*s+KJ+q, h>=s, s0<=d,
 k>=max(KI,(5c+t)*h+KF+q),
 d>=c*s0+(5c+2)*a+KG+q.                                  (AC1)
The actual evaluator starts W=I(X,prefix_s0 Y). For r=1,...,5 it computes
 T=prefix_s W,A=G(Y,T),Z=J(W,A),B=G(Y,Z),
 C=B if AC(alpha)_r=1 and C=A otherwise; then W=F(X,C).
It returns W,or any requested prefix. Seed pads retain all ignored suffixes.
This is a fixed five-round algorithm; it never computes a proof history or
a separated set. Invalid public length/label guards return the zero string.

For any source tuple L=(X0,...,Xt) independent of any seed tuple
R=(Y0,...,Yt), Hinf(X0)>=k,Y0 uniform on d bits,k<=n,and fixed advice
labels with alpha0 different from every other label, the honest h-bit output
is76eps-close to uniform jointly with ALL t other outputs and ALL t+1
original input seed strings. Correlations inside each tuple are arbitrary;
each nonhonest label occurs at most twice among the tampered labels. Thus t<=14; repetitions within that bound are allowed. The guarantee is uniform over those labels.

Here is a full proof. The conditional-product, average-entropy and robust
extraction proofs (5)--(8) work for the two WHOLE tuples L,R. Only the honest
X0,Y0 need entropy. Reveal all c initial seed prefixes, costing c*s0 right
bits; call the history H0. Strong I gives D(W0|H0,R)<=eps: X is independent
of the entire seed tuple, and the honest used prefix is uniform. Other right
coordinates append through their conditional kernel given that prefix.
Given H0 every W_i is left-sided and the two whole tuples remain independent.

Maintain a REAL rectangular history H_r containing the old W tuples and
actual selected C tuples for rounds1,...,r. Its left cost is at most c*r*h;
its right cost is at most c*s0+c*r*a. This product assertion is exact:
given old H,first reveal the old W tuple (a left function); once it is fixed,
EVERY actual C_i is a right function, even when it is a B_i. Revealing those
c right values preserves the product. Every next W_i=F(X_i,C_i) is left-sided
on the new history. Auxiliary unselected A values need not be maintained.
Indeed its unnormalized conditional density factors as
 P_L(L|H)1{Wtuple(L)=w} P_R(R|H)1{Ctuple(R,w)=c}.
This proves product directly for the kept history after forgetting auxiliaries.

Let S_r be the tampered indices which at some round<=r chose A while the
honest copy chose B. It is an analysis set determined only by the fixed labels.
The invariant,with S_0 empty,is
 D(W0 | H_r,(W_i:i in S_r), ENTIRE original R)<=delta_r.   (AC2)
In particular it never asserts independence from unseparated W_i.

At a new round first reveal ONLY current W_i for i in S_old,left-only,
then their actual chosen C_i,right-only; call the result H1. Those C_i are
right functions after their W_i are fixed. The previous invariant retains
both W_S and whole R,so appending these right functions cannot increase
its distance. Rectangle implies that the honest W0 law is unchanged by
the right revelation. Thus D(W0|H1,R)<=delta and the two sides are still
conditionally independent. Do NOT reveal unseparated W_i yet: they can
equal the honest state.

For reference the whole-tuple local look-ahead statement is
 D(B0|H1,W0,...,Wt,A0,...,At,Z0)<=2delta+6eps,             (AC-LA)
provided Havg_inf(Y0|H1)>=KG+(c-b)*a+q and h>=c*s+KJ+q, where b counts old separated copies currently choosing A.
Its complete proof is as follows. Replace W0 by exact uniform conditional
H1,preserving the conditional kernel of the other W_i given W0 and the
independent right law. This changes the LOCAL law by at most delta and
makes no new assertion about real X entropy. Reveal the c prefixes T_i.
Strong G with the honest uniform prefix gives A0 error<=2eps,retaining the
whole left tuple. The auxiliary honest W0 has average entropy at least
h-c*s>=KJ+q after those prefixes. Robust J therefore gives Z0 error<=4eps
retaining the whole right tuple. For each old separated copy currently choosing A, its W_i and C_i=A_i
were fixed in H1. These b A_i values are already known. The auxiliary
replacement keeps each such W_i fixed through the conditional other-state
kernel and leaves the independent right law unchanged. At honest-state values outside its original support,extend that kernel with the same fixed H1 coordinates; this is an arbitrary auxiliary law,not a new real-source assertion. On its support,
G(Y_i,prefix W_i) is therefore still exactly that fixed A_i. This does
not impose a new relation involving the replaced honest W0. With T fixed
every remaining A_i is right-sided; reveal only these c-b NEW A_i values,
losing at most(c-b)*a right bits. All A values are now fixed. Strong G in the symmetric
form of (6) gives B0 error<=6eps retaining the whole left tuple,including
all W_i and Z0. The real and auxiliary local laws differ by at most delta;
(8) transfers this conditional-uniform statement at cost2delta. This proves
(AC-LA). In the REAL law Haux=(H1,Wtuple,Aall) is rectangular; Z0 is determined
there and B0 is right-sided. No real-source entropy assertion was replaced
by an auxiliary one.

If the honest choice is A,apply robust G BEFORE revealing unseparated W_i.
By (AC2) its honest prefix is delta-close and independent of the right side
conditional H1. The Y entropy guard gives A0 distance<=delta+2eps jointly
with the WHOLE left tuple. Now reveal the full old W tuple; the previously
revealed W_S is a subset of it,not an additional c*h charge. Each next
separated output F(X_i,C_i),i in S_old,is left-sided because its C_i was
fixed in H1. Reveal those outputs temporarily. A0's strong distance survives
these left revelations. The X entropy guard and robust F give the next honest
state distance<=delta+4eps jointly with next separated outputs and WHOLE R.
Only AFTER that strong conclusion append the remaining actual C_i,right
functions of R and the now fixed old W tuple. Keep the actual history
(Hold,Wtuple,Ctuple); drop temporary analysis information. This is the real
rectangular history described above. It proves (AC2) for the same S with
 delta_new<=delta+4eps.                                   (AC-A)

If the honest choice is B,use (AC-LA) without first fixing unseparated W_i.
Let S_new=S_old union the indices currently choosing A. For old separated
copies their C_i was fixed in H1; for newly separated copies C_i=A_i is
fixed in Haux. Their next outputs are therefore LEFT functions. Temporarily
reveal them. Conditional product/property(7) preserves B0's distance.
The X entropy guard and robust F then give the next honest state distance
<=2delta+8eps jointly with ALL next S_new outputs and WHOLE R.
After this strong conclusion append EVERY actual selected C_i,including
unseparated B_i and the honest B0. They are right functions of R once the
old W tuple is fixed. Keep only (Hold,Wtuple,Ctuple),dropping auxiliary
unselected A_i. Since R is retained,the added C values are functions of
already retained variables; adding them and then forgetting auxiliary
information cannot increase TV. The kept history is rectangular directly
by the real left-W/right-C message order,not by an auxiliary-law assertion.
All next W_i are left-sided there. This proves (AC2) with
 delta_new<=2delta+8eps.                                  (AC-B)

For clarity,forgetting an auxiliary A does NOT discard an original seed.
The entire original R is retained in every final uniformity statement.
Entropy of Y in an intermediate G call is measured relative to the kept
rectangular history,not conditional on revealing its entire source R.
No statement claims uniformity after revealing unseparated W_i.

Here are every entropy charge used above. After r-1 rounds the maintained
left cost is c(r-1)h. The pre-revealed W_S is part of the full old W tuple,
so after that full reveal the left cost is<=c*r*h; at most t next separated
outputs add t*h. This is at most(5c+t)h. Right maintained cost is
c*s0+c(r-1)a. Let S_A,S_B partition old separated copies by their current
choice. Pre-revealing C_S costs |S|a; the look-ahead proof only reveals
c-|S_A| fresh A values, so the total at a B call is
 c*s0+[cr+|S_B|]a.                                        (AC-R)
At the first honest B, S is empty. At the second B in coordinate r,
old-separated currently-B labels contain r but not the first honest-B
coordinate. With at most two copies per other label,
|S_B|<=2(deg(r)-1),where deg=(4,4,3,3,2) are the code-coordinate degrees.
For r=2,3,4,5 the bounds cr+|S_B| are respectively
2c+6,3c+4,4c+4,5c+2,all <=5c+2 for c>=2. Coordinate1 cannot be a secondB.
An ordinary A call uses no auxiliary A tuple; it costs at most
c*s0+[c(r-1)+t]a<=c*s0+(4c+t)a<c*s0+(5c+2)a since t=c-1.
Thus every honest G call is legal under the new right guard in(AC1).
The auxiliary J call separately uses h-c*s>=KJ+q. The source ledger
remains(5c+t)h; no left information was saved or omitted.
Only honest extraction calls incur an entropy tail,not every tampered copy.

Each code has two B choices and three A choices. Starting delta0<=eps,
(AC-A),(AC-B) give final error<=76eps: moving an A before a B increases
the upper recurrence,so the worst order is three A then two B,namely
 eps ->13eps ->34eps ->76eps. At the end every tampered index lies in S
by the weight-two antichain property. Dropping the history and taking output
prefixes proves the claimed joint-output/all-original-seed guarantee.
For t14,c15,the selected source charge is89h and right charge77a. NB4.6step7 contains exactly two copies of each nonhonest label:the other seven honest rows and seven primed rows. Conditioning its proof histories changes no label,so it meets this narrower contract on every branch.
