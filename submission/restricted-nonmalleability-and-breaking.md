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

Status: a proved finite reduction CONDITIONAL ON the explicit four-dimensional
sum-product premise E below. It does not supply a numerical expansion constant
or an unconditional self-contained CZ primitive. The selectedsix-source evaluator is a restrictedvariantderivedhere from
theCZ eight-source method; its fixed-point-freecondition is explicitly
on thefirstfourcoordinates,sufficient for distinctappendedadvice. The
encoding rigidity, graph/energy estimates, fiber partition, smoothing, and
finite parameter ledger are independently proved here. No claim of external
review, tested final candidate, or submission is made.

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

E_p(gamma) IS NOT PROVED HERE. The restated CZ Theorem 6.2 in Section 7
asserts the qualitative existence of gamma>0 for these fixed thresholds.
A fully offline self-contained proof still must prove a numerical rational
gamma>0 and a numerical prime threshold p0 such that E_p(gamma) holds for
every prime p>=p0. Merely declaring such constants or searching until a finite
sample of primes passes does not establish this assertion.

The source's proof uses two further Bourgain inputs: its Lemma 7.4, filling
F_p^2 by a bounded iterated sum/product of a set of size p^(1+tau0), and its
Lemma 7.3, extracting a subset with controlled bounded sum/product growth
from simultaneous small sumset and product set. Lemma 7.3 is not proved in
CZ. The elementary BSG proof below does NOT replace either of those inputs.
The source's initial Theorem 6.2 statement and its Section 7 restatement also
have different displayed quantification of tau0,tau1; this note uses the
precise fixed-threshold premise above instead of inferring finite constants.

A direct visual check of CZ PDF page 28 (printed page 27) finds an additional
quantitative obstruction: Lemma 7.3 literally prints tau1=p^(3*k^2)*tau,
and the following proof asks for p^(3*k^2)*tau<tau1. With fixed positive
tau and k, this cannot hold uniformly over unbounded primes. This is not
just a text-extraction error. It does not refute the qualitative theorem,
but the displayed formula cannot calibrate gamma; an intended corrected
Bourgain estimate must be recovered and proved, not silently substituted.

## 2. Main conditional theorem and evaluator

Assume E_p(gamma), and choose a positive rational eta satisfying

  eta <= 1/8000,   103*eta < gamma,   p^eta >= 2^64.       (1)

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

  H_inf(W* | V*=v) >= kappa = 4*eta*log p.                   (5)

This is the joint smoothing needed for (N1) of condenser-integration.md.
No generic fixed-point/same? theorem, seeded nonmalleable primitive, or
coordinate-wise conditional claim about the original unsmoothed law is used.

Evaluator (3) takes a constant number of field operations. All proof
partitions and mixture labels are analysis objects and are never computed
by the evaluator. Section 8 gives the exact prime, bit-encoding and final
ordinary bilinear-extraction match.

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
Ateta<=1/8000,the explicit margins are1-1021eta>=6979/8000>0,
31/16-1020eta-1>=81/100>0,and31/16-1041eta-5/4>=4459/8000>0.
Thus10p^(-1+1020eta)<=10p^-eta<1/2 andallstrict/weak
projection/plane comparisons above hold. This proves the
claim. It uses E_p(gamma) only for the FINAL expansion conclusion, not
for the preceding polynomial rigidity and counting statements.

## 5. Finite distributional energy growth

Let A1,A2,C be subsets of G=(F_p^*)^4, all of sizes in [p^(2-7eta),p^2].
Assume every S subset of A1 of size at least p^(2-1000eta) satisfies
|S+S|+|S*S|>p^gamma |S|. Then every W subset of F_p^4 with |W|<=T,
T=4p^(2+7eta), has P[U_A1*U_A2+U_C in W]<=2p^(-eta).

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

  |D|>=p^(2-9eta)/8, |B'|>=p^(2-8eta)/8, |D*B'|<=2^24 p^(2+41eta).      (13)

Here 2^18 T^3/(alpha^5|R||B|)<=2^24 p^(6+26eta-4+15eta). We prove that D is
plus-good. Otherwise a witnessing C0 and W, with the same counting
additively (x+c=s1-s2+s3), give D' subset D, C' subset C0 with

  |D'|>=p^(2-10eta)/64, |C'|>=p^(2-8eta)/8, |D'+C'|<=2^27 p^(2+42eta).  (14)

Using (6) in the additive and multiplicative groups gives

  |D'+D'| <= |D'+C'|^2/|C'| <= 2^57 p^(2+92eta),
  |D'*D'| <= |D*D| <= |D*B'|^2/|B'| <= 2^51 p^(2+90eta).

After division by |D'|, their sum is at most
(2^63+2^57)p^(102eta)<2^64 p^(102eta)<=p^(103eta)<=p^gamma (p^eta>=2^64).
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
 (2+8eta)logp.                                         (15)
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
Section5 applies with C=A3: every target of at mostT=4p^(2+7eta)
points has probability<=eps_hit=2p^-eta. Every surviving mixture component
obeys this bound,so the mixture does too. No extra pair or multiplier is used.
Since eta<=1/8000 and p^eta>=2^64, 2<=T=4p^(2+7eta)<=p^4.
Choose W0 as floor(T) heaviest points. Since floor(T)>=T/2,
every point outsideW0 has mass<=2eps_hit/T=p^(-(2+8eta))=:tcap.
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
largest atom p^(-(2+8eta)), at TV cost at most 36*p^(-eta).
The Q-mass of values v having Q_V(v)<p^(-(2+4eta)) is at most p^(-4eta),
since V has at most p^2 values. On every other v,
max_w Q(w|v)<=p^(-4eta). Replace the conditional W law on bad v by
uniform on F_p^2, keeping Q_V. This costs at most p^(-4eta), and gives
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
joint laws satisfying max_w P(W=w|V=v)<=p^(-4eta) still satisfy that
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

## 8. Prime choice, exact encodings, and final ordinary extraction

Use prime-fields-and-encoding.md for the independently proved canonical
prime p(u), the least prime above 2^u, satisfying 2^u<p<2^(u+1). Its
trial-division algorithm takes O(2^(3u/2+2)*(u+1)^2) bit operations, hence
polynomial time in an outer length N when u=O(log N). No polynomial-time
claim in u is made. This module adds only a constant number of F_p
operations. Encode every field element in exactly u+1 bits and every
pair in nField=2(u+1) bits.

Because this theorem uses integer flat-source size ceil(p^(1-eta)), a
convenient extra-safe input guard for a u-bit source of entropy k is

  k >= (1-eta)*(u+1)+1.                                   (17)

Indeed log ceil(p^(1-eta)) <= (1-eta)log p+1 < the right side.
For k>=(1-rho)u, it suffices that

  (eta-rho)u >= 2-eta.                                    (18)

For example rho<=eta/2 and eta*u>=4 suffice. Distinct appended advice
is incorporated into the u-bit string before the injective embedding.
AtleastoneoftheFIRSTFOURcoordinate maps mustbefixed-point-free
onitshonestsupport;actualdistinctappendedadvice makesALLsixmaps
fixed-point-free. Sucha map can
be extended arbitrarily off support. No global no-fixed-point extension
is required by the theorem above.

Let inputs7and8 be independent u-bit sources of entropy at least
(1-rho)u, independent of thefirstSIX. Pad their concatenation with
two zero bits to nField. Their total entropy is unchanged. By (5),
kappa>=4eta*u. The previously proved bilinear-extraction.md primitive,
strong retaining this padded concatenation, gives h output bits with

  eIP <= (1/2)*sqrt((2^h-1)*2^(nField-kappa-2(1-rho)u)).    (19)

The exact reduction in condenser-integration.md Section 6.2 gives
restricted nonmalleability error at most 2*zeta+eIP.
For rho<=eta/2 and h=floor(eta*u)>=1, (19) implies

  eIP <= 2^(-eta*u),
  eNM <= 91*2^(-eta*u).                                   (20)

To check: kappa+2(1-rho)u-nField >= (4eta-2rho)u-2
>=3eta*u-2, so the exponent inside the square root in (19) is
at most -2eta*u+2; its factor 1/2 cancels the resulting factor 2.
The output obeys h<=nField. Equation (20), including the coefficient 91,
is conditionalonlyonE_p(gamma),(1),thestatedsourceguardandthe
FIRSTFOURfixed-point-freecondition. WusesthefirstSIXpolynomialinputs,
andtheindependenthashinginputisconcatenatedfrom7and8.

In the advice condenser, u=r+a and k=r-2b. Guard (17) is exactly

  eta*(r+a) >= 2b+a+2-eta.                                 (21)

One conservative regime is a<=b<=eta*r/1000, eta*r>=4.
Then the entropy deficit rho=(a+2b)/(r+a)<=3eta/1000<eta/2.
The output h=floor(eta*(r+a)) has positive linear size in r; subsequent
advice-condenser error amplification and its B=(ell+9)b ledger must
still use (20), not the original CZ's unspecified exponential constant.

The enormous finite guard p^eta>=2^64 is deliberate. It is met by
eta*u>=64. Any final family may handle smaller lengths separately only
with its own explicitly proved fallback. This note supplies no unknown
truth-table search as a claimed fallback.



# A fully expanded constant-advice strong correlation breaker

Status: complete reduction to the finite seeded extractor proved in
seeded-linear.md, with an explicit numerical instantiation of that primitive.
Section 8.1 keeps the noncircular one-tampering form; Section 8.2 proves the
direct 14-copy form that NB's shift-free inner call uses. It uses neither Li 2023's new breaker nor a black-box flip-flop or merger. It is
not a completed outer condenser, circuit candidate, or submission certificate.

Attribution: alternating extraction and look-ahead extraction are existing
techniques. The two-step A/Z/B operation and the idea of successively giving
rows the later output are from Gil Cohen, *Local Correlation Breakers and
Applications to Three-Source Extractors and Mergers*, ECCC TR15-038 (2015),
Sections 3--5 and Appendix A:
https://eccc.weizmann.ac.il/report/2015/038/download/ .
The target strong contract is Li, CCC 2019, Definition 40 and Theorem 41:
https://doi.org/10.4230/LIPIcs.CCC.2019.28 .
We replace the optimized construction there by a deliberately inefficient
constant-advice leadership schedule. The algorithm, finite ledgers, arbitrary-
copy adaptation, and proof below are independently written. This is an
exposition/specialization, not a claim to invent correlation breaking.

All logarithms are binary, all string lengths are nonnegative integers, and
TV is half L1. Write D(V|T)=TV((V,T),(U,T)), with U independent uniform of
V's length. Average conditional min-entropy means

  Havg_inf(X|T)=-log(sum_t Pr[T=t] max_x Pr[X=x|T=t]).

## 1. Primitive-parameter theorem

Fix advice length v>=1 and a set of M<=2^v advice strings (M>=2), identified
with 1,...,M in lexicographic order; advice below means an element of this
set (an advice string outside the set uses c=A in every round of Section 2).
Fix positive integers n,d,h,a,s,s0,q and 0<eps<1 such
that 2^-q<=eps and s<=h. Supply these deterministic strong seeded extractors,
all of error at most eps:

  I : {0,1}^n x {0,1}^s0 -> {0,1}^h, entropy threshold KI;
  F : {0,1}^n x {0,1}^a  -> {0,1}^h, entropy threshold KF;
  G : {0,1}^d x {0,1}^s  -> {0,1}^a, entropy threshold KG;
  J : {0,1}^h x {0,1}^a  -> {0,1}^s, entropy threshold KJ.

Strong means the entire seed is retained. Assume s0<=d and

  h >= 2s+KJ+q,                                            (1)
  k >= max(KI, (2M+1)h+KF+q),                              (2)
  d >= 2s0+2Ma+KG+q.                                      (3')

There is a deterministic map

  C : {0,1}^n x {0,1}^d x {0,1}^v -> {0,1}^h

such that for every arbitrary correlated pair (X,X') independent of every
arbitrary correlated pair (Y,Y'), with H_inf(X)>=k and Y exactly uniform,
and for every distinct advice alpha,alpha' in the fixed set,

  D(C(X,Y,alpha) | C(X',Y',alpha'),Y,Y') <= (8M+2)eps.        (4)

Neither X' nor Y' needs entropy; neither must be a deterministic function
of its honest partner. Output may be truncated to any m<=h. Without the
unequal-advice requirement, the honest output alone is still close to
uniform with error at most (8M+2)eps, retaining Y,Y'.

The separate conditions k<=n and legal primitive dimensions are required.
Section 8.1 specifies actual primitives for the one-tampering form; the
primitives of the selected 14-copy form are fixed at its instantiation,
affine-and-standard-wrappers.md (F# slots), which checks every condition.

## 2. Algorithm

On input (x,y,alpha), initialize

  w = I(x,Prefix(y,s0)).

For r=1,...,M do:

  A = G(y,Prefix(w,s));
  Z = J(w,A);
  B = G(y,Z);
  c = B if index(alpha)=r, and c=A otherwise;
  w = F(x,c).

Return w. All calls use the original x and y where indicated. There is no
source deletion, randomized setup, choice based on a proof transcript, or
search for an extractor truth table. The same algorithm is evaluated on
(x',y',alpha'). Primed variables below always mean that actual evaluation.

The only advice dependence is one choice of B among M choices. Exponential
dependence on advice length is intentional: v is a fixed outer constant.

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

## 4. Two-step look-ahead lemma, including arbitrary copies

Suppose a rectangular history H has been fixed in the analysis. Let W,W'
be h-bit functions of the left side and H, and assume

  D(W|H)<=delta,
  Havg_inf(Y|H)>=KG+2a+q.

Compute A,Z,B from (W,Y) and A',Z',B' from (W',Y') by Section 2. Under (1),
the history

  Hla=(H,W,W',A,A',Z)

is rectangular and

  D(B|Hla)<=2delta+6eps.                                   (9)

The additional left cost is at most 2h and right cost at most 2a. Z costs
nothing once W,A are fixed, and the prefixes are determined by W,W'.

Proof first for W exactly uniform conditional on H. Write
T=(Prefix(W,s),Prefix(W',s)). Reveal T on the left. Since W's prefix was
uniform before its revelation, strong extraction G and entropy smoothing
give

  D(A | H,W,W') <= 2eps.

The entire W-pair can be retained here by the symmetric form of (6): the
output's conditional distribution depends on that pair only through W's
prefix. In particular D(A|H,T)<=2eps.

Conditional on H,T, the W-pair and the right pair remain independent, and
W has average entropy at least h-2s>=KJ+q. Apply (6) to J(W,A), obtaining

  D(Z | H,T,R) <= 4eps.

With T fixed, A,A' are right-side functions, so this implies
D(Z|H,T,A,A')<=4eps. The history H2=(H,T,A,A') is rectangular. Its right
leakage is 2a, hence Havg_inf(Y|H2)>=KG+q. Z is now left-sided. Extracting
B=G(Y,Z) and using the symmetric form of (6) gives

  D(B | H2,L) <= 6eps.

This retains W,W',Z, and yields (9) with delta=0.

For the general case, construct an auxiliary law on the LOCAL tuple
(H,W,W',Y,Y'). Keep H's marginal and (Y,Y')|H unchanged. Replace W|H by
uniform, keeping W'|W,H on its support and defining it arbitrarily elsewhere.
Maintain conditional independence of the W-pair and right pair. This law
is at TV distance exactly D(W|H)<=delta from the real local tuple. The
conditional Y entropy is unchanged, so the exact-uniform proof applies.
Apply the deterministic A/Z/B algorithm and then (8), obtaining (9).
This auxiliary law is used only to prove (9); it NEVER replaces the real
source X, whose entropy is subsequently bounded in the real distribution.
There is no need to lift an absent W value to the original left side.

In the real distribution Hla is rectangular: reveal both W's on the left,
then A,A' on the right, then the already determined Z. B is right-sided
given Hla. This verifies the independence needed later, without relying
on the auxiliary law for any real entropy assertion.

## 5. Invariants and entropy ledger

Let j=index(alpha), j'=index(alpha') with j!=j'. The proof selects histories
depending on these fixed indices; the evaluator never sees them.

Initially reveal S=Prefix(Y,s0) and S'=Prefix(Y',s0), setting H0=(S,S').
Then W0=I(X,S) and W0'=I(X',S') are left-sided given H0, and strongness of I
implies

  D(W0 | H0)<=eps.                                         (10)

It also retains the full right pair because X is independent of that pair.
Initial left cost is zero and right cost at most 2s0.

After r rounds the history Hr has both Wr,Wr' left-sided, and left/right
sides independent. Its left cost is at most 2rh. Its right cost is at most
2s0+2ra before round j and 2s0+(2r+1)a from round j on (6.1 and 6.3 add 2a
per round, 6.2 adds 3a).

Before round j, maintain marginal closeness D(Wr|Hr)<=delta_r. At and after
round j, maintain joint closeness

  D(Wr | Hr,Wr')<=eta_r.                                  (11)

At any round the local extra left cost needed to prove extraction is at
most 3h, including a temporary fixing of the tampered next output. Before
round r this gives a maximum (2r+1)h<=(2M+1)h. Thus (2) supplies KF+q
residual source entropy in EVERY replenishment call in the proof.

The G calls see at most 2s0+2Ma right cost: in round r<j the call at
H_(r-1) (2s0+2(r-1)a); in round j Section 4 needs KG+2a+q at H_(j-1)
(2s0+2(j-1)a); in round r>j the call at H1=(H_(r-1),W',C')
(2s0+(2r-1)a+a, the fixing of one a-bit tampered seed included). So (3')
supplies KG+q at every G call, by (5). Consequently
all applications of (6) and (9) below meet their entropy hypotheses.

## 6. Marginal maintenance, gaining independence, preserving it

### 6.1 Before the good advice leads

Suppose r<j, with current H,W,W' and D(W|H)<=delta. The honest choice is
C=A. Reveal W,W' on the left. Robust strong G extraction, before that
revelation and retaining the whole left side, gives

  D(A | H,W,W') <= delta+2eps.                             (12)

The source Y had the required entropy by Section 5. Once W,W' are fixed,
both C=A and the actual tampered choice C' (A' OR B') are right-sided:
even J(W',A') is a function of the right side when W' is fixed.

Use (6) on F(X,A), with base history (H,W,W'). The X entropy requirement
holds after the 2h-bit left revelation. It follows that

  D(Wnew | H,W,W',R) <= delta+4eps.

Set Hnew=(H,W,W',C,C'). This history is rectangular, costs at most 2h on
the left and 2a on the right, and both new W values are left-sided. The
last distance retains Hnew, proving the next marginal invariant.
The actual C' may be correlated with A in any way. It is safely retained
as part of R, rather than wrongly assumed independent of A.

### 6.2 The good advice leads

Now r=j, so C=B whereas C'=A' because j'!=j. Apply Section 4 to obtain
Hla and D(B|Hla)<=2delta+6eps. With A' fixed, the actual tampered next
output Wnew'=F(X',A') is left-sided. Reveal it temporarily. Since B is
right-sided, (7) gives

  D(B | Hla,Wnew') <= 2delta+6eps.

The real X source retains KF+q average entropy after this extra h-bit
revelation. Apply (6) to F(X,B), retaining the ENTIRE right side:

  D(Wnew | Hla,Wnew',R) <= 2delta+8eps.                     (13)

Take Hnew=(Hla,B), but DO NOT include Wnew' in the maintained history.
It is retained separately in (11). This history's extra left cost is
2h, and extra right cost is 3a. Both new W values are left-sided given
Hnew. Project (13) to (Hnew,Wnew') to establish (11).

### 6.3 After the good advice leads

Suppose r>j and D(W|H,W')<=eta. The honest choice is always A; the tampered
choice may be either A' or B'. Reveal W' on the left. This makes the
actual C' a right-side function, so reveal C' on the right. Let

  H1=(H,W',C').

Conditioning on this right revelation does not change W's conditional law
once H,W' are fixed. Therefore

  D(W|H1)=D(W|H,W')<=eta.

The right entropy loses at most a bits. Apply robust G extraction using
the prefix of W, retaining the whole left side, to get

  D(A | H1,L) <= eta+2eps.

Both W and Wnew'=F(X',C') are left-sided given H1. Reveal them, obtaining
H2=(H1,W,Wnew'). Equation (12)'s counterpart therefore remains
D(A|H2)<=eta+2eps. Now A is right-sided because W is fixed. The total new
left cost, including W',W,Wnew', is at most 3h. Apply (6) to F(X,A):

  D(Wnew | H2,R) <= eta+4eps.

Maintain Hnew=(H,W,W',C',A), again omitting Wnew' from the history while
retaining it in (11). This costs 2h on the left and 2a on the right. The
claimed next joint invariant follows by projection.

This is the key preservation step. It does not assert that A remains
uniform after revealing an arbitrary correlated C'. Instead it fixes W'
first, then C' on the other side, and only THEN extracts A using the still
uniform honest prefix. A strong extractor retains all of the left side,
including the actual next tampered output. The ordering avoids any
unproved nonmalleable-merger claim.

## 7. Error and full joint strongness

From (10) and Section 6.1, immediately before round j,

  delta <= [1+4(j-1)]eps = (4j-3)eps.

Section 6.2 gives eta_j<=(8j+2)eps. Each of the M-j remaining rounds adds
at most 4eps. Thus

  eta_M <= (4M+4j+2)eps <= (8M+2)eps.

On every final history branch, W,W' are left-sided and R=(Y,Y') is
independent of the left tuple. Hence, exactly,

  D(W | HM,W',Y,Y')=D(W | HM,W').

Discard HM. The retained W' is the actual complete tampered evaluation,
and Y,Y' are the complete original input seeds. This proves (4), rather
than a weaker guarantee retaining only a derived seed or the tampered seed.

For a marginal-only guarantee, choose any auxiliary advice different from
the honest advice, run the already proved joint theorem with that advice
and the same seed pair, and discard its tampered output. This proves the
claimed marginal bound regardless of any advice associated with an output
that is not being retained. No extra history or entropy cost is needed.

## 8. Retained single-copy instance and selected multicopy theorem

### 8.1 CA# at the eight fixed labels (retained, unselected)

Fix v4 and advice set encode(1),...,encode(8),M8 as in NB.
For n,m,e_c>=1 and integer w>=max(2^20,e_c+2lg(4n)+14), use the
following integer slots, writing z=m only in this table:
 a=46w+ceil((31485/1048576)z)
 s=45w+ceil((5549/524288)z)
 s0=45w+ceil((33907/1048576)z)
 hc=220w+ceil((1)z)
 d=957w+ceil((617285/1048576)z)
 k=10961w+ceil((14008165/524288)z)

Here hc is Section1's h; q=e_c,eps=2^-e_c. Use I=PXR[6,287/100](inputn,outputhc,seeds0),
F=PXR[6,317/100](inputn,outputhc,seeda),G=PXR[2,147/25](input d,output a,seed s),
J=PXR[2,29/5](input hc,output s,seed a), with output truncated to m.
Each is the exact-degree strong source-linear primitive of Section7e in
linear-and-bilinear-primitives.md. Pad seeds to the slots and retain suffixes.
If k>n output0^m; otherwise run Section2's algorithm. It has error
66*2^-e_c retaining both seeds for distinct advice in the fixed set.
Every exact source threshold and padding guard is checked.

Proof. The PX constants are
 I: c=387/11968,pi=185632852377/4902092800,Lam=558845013849409/42761175875209,g=115293835544519911404377/17514977638485606400
 F: c=417/13888,pi=223995010587/5688524800,Lam=1014741853230169/104413920565969,g=217328960369425748802737/42767941863820902400
 G: c=43/122,pi=92040865327/3198156800,Lam=21609/14884,g=6163161169473/97543782400
 J: c=17/48,pi=1198172971/41943040,Lam=841/576,g=47620034669/754974720
Each row is (slot, required w coefficient, required z coefficient, ceil cost):
 a: [261571714587/5688524800]w+[417/13888]z+[417/13888].
 a: [1866640171/41943040]w+[94333/25165824]z+[17/48].
 s: [143892948527/3198156800]w+[1353855/127926272]z+[43/122].
 s0: [220506196377/4902092800]w+[387/11968]z+[387/11968].
 hc: [165926932589/754974720]w+[11059157/301989888]z+[1993/576].
 hc: [0]w+[1]z+[0].
 hc: [45]w+[5549/524288]z+[1].
 d: [1530266219093/1599078400]w+[9187666781/15607005184]z+[289521/14884].
 d: [45]w+[33907/1048576]z+[1].
 k: [468764248760256733609137/42767941863820902400]w+[2789778502851642/104413920565969]z+[2789778502851642/104413920565969].
 k: [165652477432517855212377/17514977638485606400]w+[558845013849409/42761175875209]z+[558845013849409/42761175875209].
 k: [957]w+[617285/1048576]z+[1].
For every row B_slot>=the printed z coefficient and
(A_slot-the printed w coefficient)*2^20>=the printed ceil cost, by
the displayed exact rational comparisons. Since w>=2^20 the complete
inequality follows for every integer z>=1: each component ceiling costs
at most one. These rows and constants are the full finite certificate,
not a premise imported from a parameter-search computation.

The rows respectively prove both a seeds, s, s0, h>=2s+KJ+q,
h>=m,h>=s,the right-source(3'),s0<=d,the left-source(2),KI<=k,d<=k.
They also imply KJ<=hc,KG<=d,KF<=k. Thus all primitive inputs are<=n
when k<=n, and their unit condition holds: heights<=6 give
e_c+h+6+2lg(4u)<=e_c+2lg(4n)+12<=w. The exact PXR internal padding follows
from its exact threshold guard. The full primitive-parameter theorem and
correlated-pair proof above now apply, with no change to the invariants.
This is a total finite algorithm. For its selected outer use all input,
output,error and base/hash degrees areO(logactualN); canonical exact-field
enumerations cost2^O(logN), and larger extensions use explicit towers.
This actual-length argument is not a generic polynomial bound in arbitrary
unrelated component dimensions.


## 8.2 Direct multicopy constant-advice breaker

This section extends the proved algorithm of Sections 1--7, without an
additive-shift interface. Fix t>=1 and M>=2, put c=t+1, and use the SAME
one-hot M-round algorithm of Section 2 on each input. Supply I,F,G,J as in
Section 1, each of error eps, with q>=lg(1/eps), and require

 h>=c*s+KJ+q, s<=h, s0<=d,
 k>=max(KI,(c*M+t)*h+KF+q),
 d>=c*s0+c*M*a+KG+q.                                    (MC)

If k<=n, for every source tuple (X0,...,Xt) independent of a seed tuple
(Y0,...,Yt), with Hinf(X0)>=k, Y0 uniform on d bits, and the honest advice
different from every other advice in the M-element set, the honest h-bit
output is (8M+2)eps-close to uniform jointly with ALL t other outputs and
ALL t+1 original seeds. The copies inside either tuple are arbitrary;
other advice labels may repeat. A prefix of the output has the same guarantee.

Here is a complete extension of the proof. Rectangular histories, (5)--(8),
and their proofs apply verbatim to two whole tuples instead of two pairs:
conditional independence is still between exactly two sides, and only the
honest X0,Y0 need entropy. Initially reveal all c seed prefixes, costing
c*s0 right bits. Strong I gives D(W0|H0)<=eps. After r rounds retain the
whole old W tuple and all actual chosen C seeds, so the maintained left cost
is c*r*h. The right cost before the distinguished round is c*s0+c*r*a;
from that round onward it is c*s0+(c*r+1)*a. All current W values are
left-sided. Let j be the honest advice index; at round j all other calls
choose A, since their labels differ from the honest label.

The needed multicopy look-ahead statement is
 D(B0|H,W0,...,Wt,A0,...,At,Z0)<=2delta+6eps              (MC-LA)
when D(W0|H)<=delta and Havg(Y0|H)>=KG+c*a+q. To prove it first make
W0 exactly uniform conditional on H, preserving the conditional kernel of
(W1,...,Wt) given W0 and the independent right-tuple law. This costs delta.
Reveal the c prefixes T; robust G, retaining the whole left tuple, gives
D(A0|H,W0,...,Wt)<=2eps. The honest W0 has average entropy at least
h-c*s>=KJ+q after T. Robust J therefore gives D(Z0|H,T,right)<=4eps.
With T fixed all A values are right-sided. Reveal them, losing c*a right
bits; strong G and the symmetric form of (6) give
D(B0|H,T,A0,...,At,left)<=6eps. This retains every W and Z0 and proves the
uniform case. The original and auxiliary LOCAL laws differ by delta;
(8) gives (MC-LA). In the real law its history is rectangular: reveal the
W tuple, then the A tuple; Z0 is determined. Its extra costs are c*h on
the left and c*a on the right. The auxiliary law changes no real-source
entropy assertion, exactly as in Section 4.

For r<j the honest choice is A0. Reveal the W tuple, apply robust G before
that revelation retaining the whole left tuple, then robust F. This gives
D(Wnew0|Hnew)<=delta+4eps, retaining the right tuple. Every tampered choice
is right-sided once its W is fixed, whether it is A or B. Maintain the
W tuple and the c chosen C values: costs c*h and c*a.

At r=j use (MC-LA). All tampered choices are A_i and are fixed in that
history, so all t next tampered W outputs are left-sided. Temporarily
reveal them; (7) preserves B0's distance. Robust F gives distance
2delta+8eps from uniform jointly with these t outputs and the whole right
tuple. Maintain the look-ahead history and B0, omitting the next tampered
outputs from the history but retaining them as side information. Maintained
costs are c*h and (c+1)*a; temporary left cost is (c+t)*h.

For r>j the honest choice is A0. Given the joint invariant, first reveal
all t other W values on the left and their actual C choices on the right.
The honest W remains equally close to uniform: a right revelation cannot
change its conditional law given the preceding left history. Robust G
then extracts A0 retaining the whole left side. Reveal W0 and all t next
tampered W outputs (which are left-sided because their C choices are fixed).
Robust F gives the next joint invariant with only 4eps extra error.
Maintain the old W tuple, tampered C tuple and A0, omitting the new tampered
outputs. Maintained costs are c*h and c*a; the temporary left cost is
(2t+1)*h=(c+t)*h.

At round r the largest left leakage is
[c*(r-1)+c+t]*h=(c*r+t)*h<=(c*M+t)*h. A G call sees at most
c*s0+c*M*a right bits: in the distinguished look-ahead it sees c*r*a;
afterward the previous excess a plus the t tampered C seeds gives
[c*(r-1)+1+t]*a=c*r*a. Thus (MC) pays EVERY real-source smoothing guard.
The uniform local look-ahead source pays h-c*s>=KJ+q separately. There is
one entropy tail per honest extraction, not one per tampered copy.
Before j, delta<=(4j-3)eps; at j the joint error is <=(8j+2)eps;
afterward it increases by 4eps per round. The final error is at most
(4M+4j+2)eps<=(8M+2)eps. On each final history the full output tuple is
left-sided and the entire original seed tuple is independent of it.
Appending that tuple therefore preserves distance; discard the history.
This proves the claimed joint retention, with no restrictions on copy
correlations and no retrieval of an external multicopy theorem.
