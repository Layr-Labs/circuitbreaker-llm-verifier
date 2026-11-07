# restricted-nonmalleability-and-breaking.md

Component references resolve through proof.md Section 0.

# Canonical prime selection and exact field-to-bit encodings

The finite prime certificate below uses Erdos's central-binomial proof of
Bertrand--Chebyshev.
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

  B>=4^n/(2n). (1)

For any prime p, its valuation in B is

  sum_(j>=1) [floor(2n/p^j)-2floor(n/p^j)].

Each term is zero or one. For p<=sqrt(2n), the entire prime-power
contribution is therefore at most 2n; there are at most sqrt(2n) such
primes. For p>sqrt(2n), its exponent is at most one. A prime in
(2n/3,n] has exponent zero, since floor(2n/p)=2 and floor(n/p)=1.
There is by assumption no prime in (n,2n), and 2n is composite. Thus

  B <= (2n)^sqrt(2n) * P(2n/3)
    < (2n)^sqrt(2n) * 4^(2n/3). (2)

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

## 2. The exact prime algorithm and bit cost

For integer u>=1, set p(u) to the least prime strictly greater than 2^u.
Enumerate successive integers, and for each q test division by integers
2<=d<=floor(sqrt(q)). All comparisons may be performed as d*d<=q without
floating point. A candidate is prime exactly when none of those integers
divides it: every composite has a factor at most its square root.

Section 1 proves that this search terminates with

  2^u<p(u)<2^(u+1), ceil(log p(u))=u+1. (3)

There are fewer than 2^u candidate integers, each with at most
2^((u+1)/2) trial divisors and O((u+1)^2) bit cost per division. A safe
bound is O(2^(3u/2+2)*(u+1)^2). If u<=c log N for one fixed c, this is
polynomial in the actual outer length N (not in u).

Use the standard representatives 0,...,p-1 for F_p. Addition,
multiplication and inversion by Euclid have polynomial bit cost in u;
all encodings below are fixed leading-zero binary encodings of those
representatives.

## 3. Source-rate loss from the prime interval

Embed an arbitrary u-bit source X into F_p by its integer value. This is
injective because p>2^u, so its min-entropy is unchanged. Since p>2^u, an error p^-c is at most 2^(-cu). Suppose

  H_infinity(X)>=(1-gamma)*u,
  0<=gamma<delta<1.

It meets a theorem requiring entropy at least (1-delta)log p whenever

  (delta-gamma)*u >= 1-delta. (4)

Indeed log p<u+1 and the right side of (4) rearranges to
(1-gamma)u>=(1-delta)(u+1). Condition (4) uses the POST-partition entropy.

The bit length of one output field element is u+1, not u. A pair in F_p^2
therefore uses exactly

  nField=2*(u+1) (5)

bits under the fixed injective encoding. This encoding preserves
min-entropy and TV distance exactly on the image. Any polynomial output
law is regarded as a distribution on the full nField-bit cube supported
on that image.



# Finite restricted CZ conditional entropy: one isolated deep premise

Attribution: the encoding and original eight-input polynomial map are due to Eshan
Chattopadhyay and David Zuckerman, *Non-Malleable Codes Against Constant
Split-State Tampering*, ECCC TR14-102, Sections 6--7:
https://eccc.weizmann.ac.il/report/2014/102/download/ . Their expansion theorem
extends Jean Bourgain's work over the product of two prime fields. The graph
argument below is a self-contained quantitative Balog--Szemeredi--Gowers
argument; the minimal-growth argument is Petridis's short Plunnecke method.
All logarithms are base two, TV is half L1, and all products of
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
sampler-and-disperser.md Section SP.5 proves this exact input
statement,with its full subset/projection quantifiers.

## 2. Main conditional theorem and evaluator

Assume E_p(gamma), and choose a positive rational eta satisfying

  eta <= 1/6000, (80+7*2^-20)*eta <= gamma, eta*logp >= 64*2^20. (1)

Let X1,...,X6 be independent F_p-valued sources with

  H_inf(Xi) >= log ceil(p^(1-eta)). (2)

Let fi be arbitrary maps on the respective supports into F_p. Assume at
least one of f1,...,f4 is fixed-point-free on its honest support. In particular, the
case in which every coordinate has distinct appended advice qualifies.
Define q(x)=x^4+x^2+x, e(x)=(x,q(x)), and

  E6(x1,...,x6) = (e(x1)+e(x2))*(e(x3)+e(x4))
                 +e(x5)+e(x6). (3)

Put W=E6(X1,...,X6), V=E6(f1(X1),...,f6(X6)). Their joint law is within

  zeta = min(1,45*p^(-eta)) (4)

of a law (W*,V*) on F_p^2 x F_p^2 such that, for EVERY v with positive
V*-probability,

  H_inf(W* | V*=v) >= kappa = (19/6)*eta*log p. (5)

Evaluator (3) takes a constant number of field operations. Section 8 gives the exact prime, bit-encoding and selected
weak-AGHP matching interface.

## 3. Two elementary additive-combinatorial tools

### 3.1 Small cross-sumset gives small self-sumset

In any abelian group, finite nonempty A,B satisfy

  |A+A| <= |A+B|^2/|B|. (6)

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

  |A'|>=alpha*m/8, |B'|>=alpha*n/8,

such that every (a,b) in A' x B' is joined by at least

  alpha^5*m*n/2^18 (7)

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
  |A'+B'| <= 2^18*(2L)^3/(alpha^5*m*n). (8)

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
  |A'+B'| <= 2^26*p^(2+7d+8e). (9)

The same lemma applies to multiplication inside G=(F_p^*)^4. The proofs
in this section are complete and have no sum-product assumption.

## 4. Encoding rigidity supplies the expansion hypotheses

Suppose f,g have fibers of size at most L=p^(20eta), and f has no fixed
point on its honest domain. Let A be any subset of

  {(e(x)+e(y),e(f(x))+e(g(y))): x,y in F_p}, (10)

where only the given domains are used. Every S subset of A with
|S|>=p^(2-1000eta) satisfies the hypotheses of E_p(gamma).
Consequently it expands by the factor p^gamma under sum/product.

Here and below p^eta>=2^64 implies L>16 and absorbs all displayed small
constant factors.

First, any prescribed value of e(x)+e(y) has at most four preimages.
The first coordinate gives y=s-x; substituting in q(x)+q(s-x) gives
a degree-four polynomial with leading coefficient 2, nonzero since
p>3. A fixed second coordinate-pair value in (10) has at most 4L^2
preimages. Thus for any D subset of A,

  |projection_12(D)| >= |D|/4,
  |projection_34(D)| >= |D|/(4L^2). (11)

The number of distinct points of A with a zero coordinate is at most
5p(1+L)<=10pL. For the four coordinates, fixing x leaves at most
1,4,L,4L possible y values respectively. Since
10pL/|S|<=10p^(-1+1020eta)<1/2, S has the required nonzero fraction.

For the first forbidden plane, fix y=a and write g(a)=b. Its equations
imply f(x)=c1*x+d, d=c1*a-b, and

  q(c1*x+d)+q(b)-c2*(q(x)+q(a))=0. (12)

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
Ateta<=1/6000,the explicit margins are1-1021eta>=4979/6000>0,
31/16-1020eta-1>=307/400>0,and31/16-1041eta-5/4>=257/500>0.
Thus10p^(-1+1020eta)<=10p^-eta<1/2 andallstrict/weak
projection/plane comparisons above hold. This proves the
claim. It uses E_p(gamma) only for the FINAL expansion conclusion, not
for the preceding polynomial rigidity and counting statements.

## 5. Finite distributional energy growth

Putnu=2^-20,a_set=5+nu,tau=16/3,alpha=p^-eta,T=4p^(2+tau*eta).
LetA1,A2,C beunitsets inG=(F_p^*)^4 withsizes in
[p^(2-a_set*eta),p^2]. AssumeEVERYS subsetA1 withsize>=p^(2-1000eta)
has|S+S|+|SS|>p^gamma|S|. ThenEVERYW subsetF_p^4,|W|<=T,satisfies
 Pr[U_A1*U_A2+U_C inW]<=2alpha.
CallD plus-good ifPr[U_D+U_C0 inW]<=alpha forEVERYunitC0 inthat
sizerange andEVERY|W|<=T;define times-good analogously. Disjoint
unionsofplus-goodsets areplus-good bytheirsize-weightedmixture.
Peelplus-good subsetsfromR=A1 untilR istimes-good or|R|<=alpha|A1|.
Ifneitherholds,|R|>=p^(2-(a_set+1)eta),andawitnessB,W hasproduct
edge-density>alpha. Section3.2 andthethree-edgewalkidentity
xb=s1*s2^(-1)*s3 giveD subsetR,B' subsetB with
 |D|>=p^(2-(a_set+2)eta)/8, |B'|>=p^(2-(a_set+1)eta)/8,
 |DB'|<=2^18*T^3/(alpha^5|R||B|)
       <=2^24*p^(2+(3tau+2a_set+6)eta). (13)
Forfixedx,b,thetriple fixesb1=x^(-1)s1 andx1=s2*b1^(-1);
distinctproducts usedisjointtriples. Everyedgelabelisaunit,evenifW
containsnonunits. IfD isnotplus-good,itsadditivewitnessC0,W and
theidentityx+c=s1-s2+s3 giveD' subsetD,C' subsetC0 with
 |D'|>=p^(2-(a_set+3)eta)/64, |C'|>=p^(2-(a_set+1)eta)/8,
 |D'+C'|<=2^27*p^(2+(3tau+2a_set+7)eta). (14)
Thefactor8 fromthelowerboundon|D| ischargedinthedenominator.
Apply(6)additivelytoD',C',andmultiplicativelytoD,B':
 |D'+D'|<=2^57*p^(2+(6tau+5a_set+15)eta),
 |D'D'|<=|DD|<=2^51*p^(2+(6tau+5a_set+13)eta).
OnlyD'D' subsetDD isused;noset-distributivityassertionoccurs.
DividingbytheLOWERboundon|D'| andusingp>=1,theirrelativesum is
 <=(2^63+2^57)p^((6tau+6a_set+18)eta)
 <2^64*p^((80+6nu)eta)<=p^((80+7nu)eta)<=p^gamma.
Here nu*eta*logp>=64 pays2^64 EXACTLY. Also
(1000-a_set-3)eta*logp=(992-nu)eta*logp>6,so
|D'|>=p^(2-1000eta),contradictingitshereditaryexpansion.
ThusD isplus-good. Eachpeelremovesanonemptyset,sofinitepeeling
terminates;thisisanalysis,notanevaluatoralgorithm.
LetA+ betheunionofpeeledsetsandAx theremainder;omit emptyzero-weight
components. Forfixeda2 inA2,
multiplicationbyitsunitinverse gives
 Pr[U_A+*a2+U_C inW]=Pr[U_A+ +a2^(-1)U_C ina2^(-1)W]<=alpha;
thepartnerandtargetsizesarepreserved. IfAx istimes-good,fixc and
useW-c. Otherwiseitsmixtureweightis<=alpha. Averagingprovesthe
stated2alpha boundwithallpartnerandtargetquantifiers.

## 6. Small fibers for the restricted six-source map

SupposeindependentflatXi havesupportsizes>=p^(1-2eta),eachfi has
fibers<=p^(20eta),andonefirst-fourmapisfixed-point-free. Swapwithin
afirstpairorbetweenthetwofactors toputthatcoordinatefirst;nosymmetry
withtheadditivepairisused. Foreachof(1,2),(3,4),(5,6),partitionits
jointpairimagebymultiplicityr=1,2,3,4. Thehonestquarticpreimagebound
makesthisexact:animagepartA_r hasweight r|A_r|/M andisuniform,
whereM>=p^(2-4eta)istheinputsupportproductsize.
Discardpartswith|A_r|<p^(2-5eta). Theirpairmass isatmost
(1+2+3+4)p^(2-5eta)/p^(2-4eta)=10p^-eta. Conditioningoneachpair's
retainedcomponents keepspairsindependentandcostsexactlyitsdiscarded
massinTV. ProductTVcostis<=30p^-eta;thereisnodivisionbykeptmass.
Eachretainedimagehas<=10p^(1+20eta)zero-coordinatepoints,afraction
<=10p^(-1+25eta)<=p^-eta:1-27eta>0andp^eta>=2^64 pay10.
Replaceeachcomponentbyuniformonitsunitsubset,keepingitscurrent
mixtureweight. Thiscosts<=3p^-eta acrosspairs. Withnu=2^-20,
 (1-p^-eta)p^(2-5eta)>=p^(2-(5+nu)eta),
because1-p^-eta>1/2 andp^(-nu*eta)<=2^-64<1/2 bytheexistingguard.
Uppercardinalityremainsp^2. Conditionedoncomponentlabels,thejoint
outputisEXACTLYU_A1*U_A2+U_A3. Section4supplieshereditaryexpansion
forA1;Section5giveshittingerror eps_hit=2p^-eta onEVERYtargetofsize
<=T=4p^(2+(16/3)eta). Mixingpreservesthebound.
Ateta<=1/6000,2<=T<=p^4and2+(19/3)eta<=4;theguardp^eta>=2^64
pays4. Thefloor(T)heaviestpointsW0 havemass<=eps_hit. Everyoutside
point hasmass<=2eps_hit/T=p^(-(2+(19/3)eta))=:tcap.
Clippingtotcapremovesmass<=eps_hit;ambientp^4hasenoughcapacityto
redistributeitbelowtcap. TheTVcostistheremovedmass. Totalcostis
30p^-eta+3p^-eta+2p^-eta=35p^-eta<36p^-eta. Thus(W,V)iswithin
36p^-eta ofajointlawofmin-entropy(2+(19/3)eta)logp. (15)

## 7. Arbitrary fibers, arbitrary sources, and conditional smoothing

First take independent flat Xi with support Si of size at least
ceil(p^(1-eta)). Partition each Si into

  Si,0 = {x in Si: |Si intersect fi^(-1)(fi(x))|<=p^(20eta)},
  Si,1 = Si \ Si,0. (16)

These parts are unions of COMPLETE fibers. The argument of fi^(-1) is
fi(x), not x.
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
finite prime proof, in time polynomial in actualN (Section2). Encode the pairW in2(u+1)bits.
For the integer flat-source threshold ceil(p^(1-eta)),the safe guard
 k>=(1-eta)(u+1)+1 (17)
implies(2),since logceil(p^(1-eta))<=(1-eta)logp+1.
When k>=(1-rho)u,it suffices that(eta-rho)u>=2-eta; inparticular
rho<=eta/2 andeta*u>=4 suffice. Distinct fixed appendedadvices ensure
first-four fixed-point-freeness on actual honest supports; no global
extension offsupport is required. The extra ring guardu>=3200000 is
separate from the energy constant-absorption guardeta*logp>=64*2^20,which every selected call supplies viaeta*u>=64*2^20. This stronger guard also impliesp^eta>=2^64.

The selected seven-input weak-AGHP match (first six CZ inputs plus one
raw r-bit seed) is proved in affine-and-standard-wrappers.md Section3.4.

# Conditional tools and the selected antichain multicopy breaker

The local A/Z/B operations adapt Cohen's local correlation breakers and the
whole-tuple MC-LA proof of publishedc78568c1; all proofs are supplied below.

All logarithms are binary and TV is half L1. WriteD(V|T)=TV((V,T),(U,T)),
where U is independent uniform on V's length. Define
 Havg_inf(X|T)=-log(sum_t Pr[T=t] max_x Pr[X=x|T=t]).

## 1. Primitive interfaces

Fix positive integers n,d,h,a,s,s0,q and0<eps<1 with2^-q<=eps,s<=h.
Supply these deterministic strong seeded extractors,all error at mosteps:
  I : {0,1}^n x {0,1}^s0 -> {0,1}^h, entropy threshold KI;
  F : {0,1}^n x {0,1}^a -> {0,1}^h, entropy threshold KF;
  G : {0,1}^d x {0,1}^s -> {0,1}^a, entropy threshold KG;
  J : {0,1}^h x {0,1}^a -> {0,1}^s, entropy threshold KJ.

Strongness retains the complete seed. All pads retain ignored suffixes.

## 2. Local alternating operations

On a current h-bit stateW and d-bit original seedY,compute
 T=prefix_s W,A=G(Y,T),Z=J(W,A),B=G(Y,Z).
The selected round usesC=A orC=B according to the fixed public label code,
thenWnew=F(X,C). Initialization isI(X,prefix_s0Y). Section8.2 proves the schedule.

## 3. Elementary conditional tools, with proofs

Let L=(X,X') and R=(Y,Y'). A rectangular history is a sequence of messages,
each a deterministic function of one whole side and previous messages.
Initially the sides are independent. On conditioning on a left message,
the left factor in the conditional product law changes and the right factor
does not; the symmetric statement holds for a right message. Induction
proves that L and R are independent conditional on every history value.
This works for arbitrary correlated tuples inside each side.

If the total left-message length is bL, then

  Havg_inf(X|H)>=H_inf(X)-bL. (5)

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

  D(E(X,P) | H,R) <= eta+2eps. (6)

For fixed H=t, write g_t(p)=TV(E(X|H=t,p),U), which lies in [0,1]. Conditional
independence implies that the left side of (6) averages g_t(P) even when the
ENTIRE right side is retained. Replacing P's marginal by uniform changes
this expectation by at most TV(P|t,U). On the good entropy histories the
uniform-seed expectation is at most eps; the other histories cost at most
eps altogether. Average to obtain (6). The symmetric rule retains the entire left side
when extracting from the right.

If V is right-sided and Q is left-sided conditional on H, then

  D(V|H,Q)=D(V|H). (7)

This follows by averaging conditional TV over H: V and Q are independent on
each branch. Finally, if joint laws P,Q on (V,T) have TV
at most eta, the triangle inequality gives

  D_P(V|T) <= D_Q(V|T)+2eta. (8)

## 8. Selected antichain multicopy instance

NB's inner map uses theorem8.2 with the PXG instance of wrappers4.2.

## 8.2 Selective-transcript five-round multicopy theorem

Equation numbers in this subsection are local.

### 8.2.1. Statement

Fix the eight public label codes, in label order 1 through 8:
11000, 10100, 10010, 01100, 01010, 00110, 10001, 01001.
The bit 1 means select B; the bit 0 means select A.

Let t>=1, c=t+1. The honest label differs from every other label, and each
nonhonest label has multiplicity at most 2. Necessarily t<=14.
All labels, dimensions, and maps below are fixed, not selected from private data.

Let n,d,h,a,s,s0,q be positive integers, 0<eps<1, 2^(-q)<=eps, s<=h,
s0<=d. Supply deterministic strong seeded extractors of error at most eps:

 I: n source bits, s0 seed bits, h output bits, threshold KI;
 F: n source bits, a seed bits, h output bits, threshold KF;
 G: d source bits, s seed bits, a output bits, threshold KG;
 J: h source bits, a seed bits, s output bits, threshold KJ.

"Strong" retains the complete seed of the padded map. Every primitive has its
own legal dimension/field/threshold hypotheses. Require

 h >= 11s+KJ+q,
 k >= max(KI, (2c+t)h+3cs+KF+q), k<=n,
 d >= c*s0+5c*a+KG+q. (1)

For every source tuple L=(X0,...,Xt) on n-bit strings and seed tuple
R=(Y0,...,Yt) on d-bit strings such that the TWO WHOLE TUPLES are independent,
Hinf(X0)>=k, and Y0 is uniform, run each copy with its assigned label:

 W_i = I(X_i, Prefix_s0(Y_i)).
 For r=1,...,5:
   T_i=Prefix_s(W_i);
   A_i=G(Y_i,T_i);
   Z_i=J(W_i,A_i);
   B_i=G(Y_i,Z_i);
   C_i=B_i for code bit 1, C_i=A_i for code bit 0;
   W_i=F(X_i,C_i).

Then the final honest h-bit output satisfies

 TV((W_0,W_1,...,W_t,R), (U_h,W_1,...,W_t,R)) <= 76eps. (2)

Here U_h is independent of the displayed original target marginal.
Taking any requested output prefix preserves the guarantee. All ignored suffixes of all c original seeds remain in R.

There is no entropy requirement on X_i or Y_i for i>0, no within-side
independence, and no functional relation to X0 or Y0. In particular sources
or seeds can coincide.

### 8.2.2. Conditional tools, including retained-side hypotheses

Write D(V|T)=TV((V,T),(U,T)), with the ideal law retaining the actual T
marginal. All variables have finite alphabets; zero-probability histories
can be omitted. Define

 2^(-Havg(X|H)) = sum_h P(h) max_x P(X=x|H=h).

A rectangular history is a sequence of messages, each a deterministic
function of one WHOLE side and previous messages. If L and R begin
independent, their conditional density after any fixed transcript factors:
a left message restricts only the left factor and a right message only
the right factor. Thus L and R are independent on every history branch.
The same argument permits independent private coins assigned to their
respective side.

A left message with at most 2^b possible values increases
sum_h P(h) max_x P(X0=x|h) by at most 2^b. Indeed, for fixed previous h,

 sum_m max_x P(X0=x,M=m|h) <= 2^b max_x P(X0=x|h).

This holds even if the message is a function of other correlated coordinates
of L. Right messages leave the conditional X0 law unchanged. Induction
therefore proves Havg(X0|H)>=Hinf(X0)-bL, where bL is total left-message
length. Symmetrically Havg(Y0|H)>=d-bR. The induction permits adaptive
message functions; no independence inside a side is used. Messages already
determined by the history cost zero.

Havg(X|H)>=K+q implies Hinf(X|H=h)>=K outside history mass at most
2^(-q)<=eps, by Markov's inequality on the conditional largest atom.

Robust extraction rule. Suppose L,R are independent given rectangular H,
S is a right function given H with D(S|H)<=eta, X is a left coordinate,
Havg(X|H)>=K+q, and E is strong (K,eps). Then

 D(E(X,S)|H,R) <= eta+2eps. (3)

For fixed H=h and seed s, let
g_h(s)=TV(E(X|h,s),U). Conditional independence means that retaining the
ENTIRE right side gives the average of g_h(S), since conditional source
laws depend on R only through h. Replacing S's marginal by uniform changes
this expectation by at most TV(S|h,U). On good-entropy h the uniform-seed
average is at most eps; the exceptional h have total mass at most eps.
Averaging proves (3).
Exchanging the sides gives the symmetric rule, retaining all of L.

Two additional facts will be used:
(a) If V is right-sided and Q is left-sided given rectangular H, then
D(V|H,Q)=D(V|H). Revealing Q preserves product and the right conditional law.
(b) If two joint laws on (V,T) are within eta, then
D_P(V|T)<=D_Q(V|T)+2eta. One eta compares the joint laws, the second compares
their T marginals after adjoining independent uniform V.

Projection of retained side information contracts D. Projection of a
rectangular transcript need NOT remain rectangular. Every compressed real
history below is separately proved rectangular by its message order.

### 8.2.3. The kept real history and separated-set invariant

Reveal all c initial s0-bit seed prefixes as H0. This is right-only.
Strong I with the honest uniform prefix gives

 D(W0|H0,R)<=eps. (4)

For completeness, conditional output laws of W0 depend on R only through
the honest prefix because X0 is independent of the WHOLE R. Appending the
rest of R using its conditional kernel given that prefix therefore keeps
the strong bound. All current W_i are left
functions given H0.

In each round define the actual selected summary Q_i to be:
 - Prefix_s(W_i) for an A-choice copy;
 - the complete h-bit W_i for a B-choice copy.

The kept real history H_r consists of H0 and, for each completed round,
all c selected summaries followed by all c actual selected C_i values.
Given old H, all current W_i are left functions. Thus Q is a left message.
After Q is fixed, every selected C_i is a right function:
an A-choice needs only its fixed prefix; a B-choice has its full state
fixed, so A_i,Z_i,B_i are functions of R. Revealing selected C is right-only.
Consequently H_r is genuinely rectangular. Each next W_i=F(X_i,C_i) is
left-sided on H_r. Unselected A,Z,B values are evaluated by the algorithm
but are not revealed in this history.

The right cost after r rounds is at most c*s0+c*r*a.
Each copy's code has two B's and three A's. At every prefix of the rounds,
the total left summary cost is therefore at most

 2c*h+3c*s. (5)

Let S_r be the tampered indices that have chosen A in at least one
completed round in which the honest copy chose B. This set depends only
on labels. Maintain

 D(W0 | H_r, W_(S_r), R) <= delta_r. (6)

Initially S_0 is empty and delta_0<=eps by (4).

Before a round, let H=H_(r-1), S=S_(r-1). Temporarily reveal W_S on the
left, then current actual selected C_S on the right; call this H1.
C_S is right-sided after W_S is fixed, regardless of the choices of those
copies. Equation (6), which already retained W_S and all R, implies

 D(W0|H1,R)<=delta. (7)

The real whole sides remain independent given H1. Its right cost is at
most c*s0+c*(r-1)*a+|S|a. Whole W_S is TEMPORARY, not kept in H_r.

### 8.2.4. Selective local look-ahead, including the auxiliary law

Suppose the honest choice is B. Let
Jset={0} union {tampered indices currently choosing A}, and let

 N={0} union {tampered indices outside S currently choosing A},
 nu=|N|.

All old separated A-choice copies have their W_i and C_i=A_i fixed in H1.
Assume, for this local step,

 Havg(Y0|H1)>=KG+nu*a+q, h>=nu*s+KJ+q. (8)

We prove the REAL-law bound

 D(B0|H1, W0,...,Wt, (A_i:i in Jset), Z0)
       <= 2delta+6eps. (9)

Given H1, the whole current-state tuple and R are independent. Replace
the honest STATE W0 by exact uniform conditional on H1, keeping H1's
marginal and the conditional kernel
Law((W_i:i>0)|W0,H1), and keeping Law(R|H1) independently.
For honest-state values outside original support, extend this kernel
arbitrarily while keeping every old separated coordinate equal to its
fixed H1 value. This extension always exists: S omits coordinate 0.
The resulting local joint law differs from the real local law by exactly
D(W0|H1), at most delta by (7) and conditional product.

This law contains current states and seeds. It need not admit the old
state/source equations, and it makes NO assertion about real X0 entropy.
Because W_S stays fixed and the conditional R law is unchanged, every
already fixed old C_S equation remains true. In particular old separated
A_i values still cost zero. An independent replacement of all states,
or an extension that changes fixed separated coordinates, would be invalid.

In the auxiliary law W0 is uniform given H1, and the state tuple is
independent of R given H1. Compute the nu new prefixes T_i, i in N,
and write M for this prefix vector. Old separated prefixes are known.
The first G extraction is justified BEFORE conditioning on M: its honest
T0 is exactly uniform given H1. For fixed H1 the conditional A0 law given
the entire state tuple depends on it only through T0. On good Y-entropy
branches, strong G thus gives error eps even retaining the whole state
tuple; bad branches cost eps. Therefore

 D_aux(A0|H1, whole state tuple)<=2eps,
 hence D_aux(A0|H1,M)<=2eps. (10)

On H2=(H1,M), the state tuple and R remain independent. Uniform W0
given H1 and the nu*s-bit left message M imply
Havg_aux(W0|H2)>=h-nu*s>=KJ+q. A0 is now right-sided.
Robust J gives

 D_aux(Z0|H2,R)<=4eps. (11)

Reveal only the nu new A_i for i in N. They are right functions after
M is fixed; let H3=(H2,A_N). All other A_i in Jset were already known
in H1. Since (11) retains R, appending A_N costs no extra TV.
On H3, Z0=J(W0,A0) is left-sided, and its seed distance is at most 4eps.
The right entropy loss in this step is at most nu*a. Prefix M was left-only,
so (8) gives Havg_aux(Y0|H3)>=KG+q. Product still holds. Symmetric robust
G gives

 D_aux(B0|H3, whole state tuple)<=6eps. (12)

The retained information here determines Z0, and contains precisely all
A_i in Jset; M is already a function of the full state tuple.
Thus (12) is the auxiliary version of (9). Deterministic evaluation of
all relevant local A/Z/B values contracts the delta distance between
local laws. Fact (b) transfers the conditional-uniform statement, including
its changed target marginal, at cost 2delta. This proves (9).

In the REAL law the information (H1,whole state tuple,A_Jset) has its own
rectangular order: reveal states on the left, then selected A values on
the right. Z0 is determined there. This fact is separate from the
auxiliary-law argument.

### 8.2.5. The universal 11-prefix bound

View each weight-two code as a two-element subset of {1,...,5}.
The coordinate degrees among the eight labels are (4,4,3,3,2).

At the first honest B, S is empty. Its coordinate can be only 1, 2, or 3,
whose degree is at least 3. There are at most 8-3=5 other labels omitting
that coordinate, each with at most two copies. Hence nu<=1+2*5=11.

At the second honest B, let the first and second coordinates be i<j.
Every unseparated copy contains i. To choose A now it must omit j.
Among labels containing i, only the honest label contains j as well:
a weight-two set containing both is exactly {i,j}. Thus at most deg(i)-1
other labels can be unseparated and choose A. Consequently
nu<=1+2*(deg(i)-1)<=7.

These statements hold for missing labels and arbitrary multiplicities
0,1,2. They prove the local h guard in (8) from h>=11s+KJ+q.

### 8.2.6. Projection to a real F history

Define Snew=S in an honest A round, and
Snew=S union {tampered indices currently choosing A} in an honest B round.
Let P be the REAL history obtained from old kept H by revealing:
 1. all current selected summaries Q;
 2. actual selected C_i only for i in Snew.

This is independently rectangular: Q is left-only, and each of these C_i
is right-only once its selected summary is fixed. It contains every next
separated output's seed, but does not append the honest selected seed yet.

If the honest choice is B, P is a deterministic projection of the retained
information in (9). Q is a function of the full W tuple. For old S, C_S
was fixed in H1; newly separated copies choose A and their C_i=A_i is
among A_Jset. Therefore (9) implies

 D(B0|P)<=2delta+6eps. (13)

Because honest Q0 is its full W0, B0 is right-sided on the REAL P.

If the honest choice is A, apply symmetric robust G at H1 BEFORE revealing
unseparated states. By (7) its honest seed Prefix_s(W0) has distance at
most delta. The Y entropy bound proved below gives

 D(A0|H1,L)<=delta+2eps. (14)

All current summaries are functions of retained L and old H. Here Snew=S,
so P is a deterministic projection of (H1,L). Thus
D(A0|P)<=delta+2eps. The honest summary contains its prefix, so A0 is
right-sided on the real P.

In either case, each next separated output
W_i^+=F(X_i,C_i), i in Snew, is left-sided on P, since its C_i is fixed.
Reveal these at most t h-bit outputs on the left. Product survives and
fact (a) preserves the honest selected seed's distance. Call the result
Pplus. Once the real source entropy bound below is supplied, robust F
therefore yields, retaining all R and all next separated outputs,

 honest A: delta_new <= delta+4eps;
 honest B: delta_new <= 2delta+8eps. (15)

Now append all remaining actual C_i, including the honest selected seed.
They are functions of retained R and selected summaries in P. This adds
no TV cost after the seed-retaining F conclusion. Drop temporary
information and keep exactly H_r and W^+_(Snew), still retaining all R.
This proves invariant (6) for the next round.

The independence at F comes from the REAL message order defining Pplus.

### 8.2.7. Every entropy ledger and error recurrence

For each real F call, H contains past selected summaries, P adds only
current selected summaries, and Pplus adds at most t next outputs.
By (5) its total left cost from original X0 is at most

 2c*h+3c*s+t*h = (2c+t)h+3cs. (16)

Its right messages do not change the conditional left law. Equation (1)
therefore gives Havg(X0|Pplus)>=KF+q. Initialization needed only KI.

For an honest B in round r, H1 costs at most
c*s0+c*(r-1)a+|S|a right bits. The local proof adds only nu new A values.
N and S are disjoint index sets among the c copies, so |S|+nu<=c. Hence

 c*s0+[c*(r-1)+|S|+nu]a <= c*s0+c*r*a <= c*s0+5c*a. (17)

This proves (8) for every local G call. The first local G has no greater
right cost; prefixes cost no right entropy. The auxiliary law keeps the
same H1 marginal and right conditional law, so this entropy statement
is valid there as well.

For an honest A call the right cost is at most
c*s0+[c*(r-1)+|S|]a <= c*s0+[c*(r-1)+t]a < c*s0+5c*a.
Thus its G call in (14) is legal. No unselected A_i of an old separated
B-choice copy is revealed. A newly separated A-choice value is charged
once, and an already known old separated A-choice value is free.

The h guard is Section8.2.5. The robust tail at each honest extraction
has mass at most eps; no tail is charged per tampered copy.

Starting delta0<=eps, (15) has A:delta+4eps, B:2delta+8eps.
An adjacent A followed by B gives 4eps more than B followed by A, so
among all three-A/two-B orderings the maximum has AAA BB:
eps ->13eps ->34eps ->76eps.
For the actual ordered code list the final coefficients are
40,44,48,52,56,64,52,60, all at most 76.

Every different weight-two label has an honest-1/other-0 coordinate.
Therefore every tampered copy is in S_5. Invariant (6) retains all other
outputs and ALL c original seed strings at the end. Dropping H_5 proves
(2). Output truncation contracts TV.

At t=14,c=15, (17), (16), and the local h guard are respectively:
75a; 44h+45s; 11s+KJ+q.

