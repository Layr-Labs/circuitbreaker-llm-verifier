# restricted-nonmalleability-and-breaking.md

Component references resolve through proof.md Section 0.

<!-- COMPONENT prime-fields-and-encoding.md SELECTED SHA256 a0f9d1338b22df92cf8787359b89e61a583f2161d8dfb9a947664ce9e450c621 -->

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



<!-- COMPONENT cz-conditional-entropy.md SELECTED SHA256 7728410551b00af990f06320a68e7bc1b99eafa3bcd75128583372d60aaf503d -->

# Finite restricted CZ conditional entropy: one isolated deep premise

Status: a proved finite reduction CONDITIONAL ON the explicit four-dimensional
sum-product premise E below. It does not supply a numerical expansion constant
or an unconditional self-contained CZ primitive. The polynomial evaluator,
encoding rigidity, graph/energy estimates, fiber partition, smoothing, and
finite parameter ledger are independently proved here. No claim of external
review, tested final candidate, or submission is made.

Attribution: the encoding and eight-input polynomial map are due to Eshan
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

  eta <= 1/100000,   1000*eta <= gamma,   p^eta >= 2^64.       (1)

Let X1,...,X8 be independent F_p-valued sources with

  H_inf(Xi) >= log ceil(p^(1-eta)).                          (2)

Let fi be arbitrary maps on the respective supports into F_p. Assume at
least one fi is fixed-point-free on its honest support. In particular, the
case in which every coordinate has distinct appended advice qualifies.
Define q(x)=x^4+x^2+x, e(x)=(x,q(x)), and

  E(x1,...,x8) = (e(x1)+e(x2))*(e(x3)+e(x4))
                +(e(x5)+e(x6))*(e(x7)+e(x8)).               (3)

Put W=E(X1,...,X8), V=E(f1(X1),...,f8(X8)). Their joint law is within

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
All numerical comparisons follow from eta<=1/100000. This proves the
claim. It uses E_p(gamma) only for the FINAL expansion conclusion, not
for the preceding polynomial rigidity and counting statements.

## 5. Finite distributional energy growth

Let A1,A2,C be subsets of G=(F_p^*)^4, all of sizes in
[p^(2-7eta),p^2]. Assume every S subset of A1 of size at least
p^(2-1000eta) satisfies |S+S|+|S*S|>p^gamma |S|.
Then U_A1*U_A2+U_C is p^(-eta)-close to a law with collision probability
at most p^(-(2+10eta)).

Call a set D plus-good if cp(U_D+U_C0)<=p^(-(2+10eta)) for every
C0 subset of G of size in [p^(2-7eta),p^2]. Define times-good similarly.
Any disjoint union of plus-good sets is plus-good, since a uniform law
on the union is a convex mixture and the squared L2 norm is convex.

Start with R=A1 and peel off disjoint plus-good subsets, until R is
times-good or |R|<=p^(-eta)|A1|. If R is still large, its size is at
least p^(2-8eta). If it is not times-good, choose a witnessing B of size
at least p^(2-7eta). Formula (9), multiplicatively with d=8eta,e=10eta,
gives D subset of R and B' subset of B with

  |D|,|B'| >= p^(2-27eta),  |D*B'| <= p^(2+137eta).          (13)

Here 16 and 2^26 in (9) are each at most p^eta.
We prove that D is plus-good. Otherwise choose a witnessing C0, and
apply (9) additively with d=27eta,e=10eta, obtaining D' subset of D,
C' subset of C0 with

  |D'|,|C'| >= p^(2-65eta), |D'+C'| <= p^(2+270eta).         (14)

Using (6) in the additive and multiplicative groups gives

  |D'+D'| <= p^(2+605eta),
  |D'*D'| <= |D*D| <= |D*B'|^2/|B'| <= p^(2+301eta).

After division by |D'|, their sum is at most
p^(670eta)+p^(366eta)<=p^(671eta)<p^gamma. But
|D'|>=p^(2-65eta)>=p^(2-1000eta), contradicting its expansion.
Thus D is plus-good. The peeling removes at least one point each time
and therefore terminates; no algorithmic runtime is asserted for this
analysis-only decomposition.

Let A+ be the union of peeled pieces. If a times-good remainder exists,
call it Ax. Otherwise discard the small remainder, whose relative mass
is at most p^(-eta). On the A+ component fix a2 in A2. Its output law is
U_A+*a2+U_C = a2*(U_A+ + a2^(-1)*U_C), so it has the required collision
bound: multiplication by a2 is a bijection and a2^(-1)C is an allowed
set in the plus-good definition. On the Ax component fix c in C;
times-goodness gives the same bound, unchanged by translation. Convexity
of collision probability gives the bound after mixing. Replace discarded
mass by this good mixture. Its TV cost is at most p^(-eta), proving the
claim.

## 6. Small fibers: the finite first-eight joint-entropy estimate

Suppose Xi are independent flat sources on sets of sizes at least
p^(1-2eta), every fi has fibers of size at most p^(20eta) on its support,
and at least one fi is fixed-point-free there. Then (W,V) for (3) is
within 36*p^(-eta) of a law with joint min-entropy at least

  (2+8eta)*log p.                                          (15)

The expression (3) is invariant under swapping variables within a pair,
swapping its two factors in either product, and swapping the two product
terms. These permutations move any chosen fixed-point-free coordinate
to coordinate 1 without changing the function.

For each of pairs (1,2),(3,4),(5,6), partition the image of
(e(Xi)+e(Xj),e(fi(Xi))+e(fj(Xj))) by its multiplicity r, with 1<=r<=4
as proved in Section 4. Conditional on such a part, the image is uniform.
Discard parts with image size below p^(2-6eta). Their total mass per
pair is at most (1+2+3+4)*p^(2-6eta)/p^(2-4eta)
=10*p^(-2eta). At least one part survives, since this bound is below one.
Renormalize the surviving mixture, at that TV cost. Independence allows
doing this on all three pairs, for total cost at most 30*p^(-2eta).

For each surviving flat image A, at most 10p^(1+20eta) points have a zero
coordinate. Its relative bad mass is at most 10p^(-1+26eta)<=p^(-eta).
Delete them and renormalize; three such deletions cost at most 3*p^(-eta).
Each resulting set has size at least p^(2-7eta), since
1-p^(-eta)>=p^(-eta). All have size at most p^2.

For the original pair (7,8), the probability that its four-coordinate
image has any zero coordinate is at most
5*(1+p^(20eta))/|S8|<=10*p^(-1+22eta)<=p^(-eta).
Discard that event and fix its surviving value alpha in G. The desired
four-coordinate joint output is now A1*A2+alpha*A3. Section 4 verifies
the hereditary expansion condition for A1; multiplication by alpha
preserves the size of A3. Section 5 gives a further p^(-eta) TV cost
to a collision-probability bound p^(-(2+10eta)). All mixture components
obey that bound, so their mixture does too. The total cost so far is
at most 30*p^(-2eta)+5*p^(-eta)<=35*p^(-eta).

For any probability law P with cp(P)<=c, clipping each atom to a cap t
removes mass at most c/t: the removed mass is bounded by
sum_{P(z)>t} P(z)<=sum P(z)^2/t. If the ambient space has at least 1/t
points, redistribute that removed mass among unfilled capacities to
obtain a probability law of largest atom at most t. This changes TV by
at most c/t. Take t=p^(-(2+8eta)); the ambient p^4 space has sufficient
capacity, and c/t<=p^(-2eta). This proves (15) and the stated 36 bound.

## 7. Arbitrary fibers, arbitrary sources, and conditional smoothing

First take independent flat Xi with support Si of size at least
ceil(p^(1-eta)). Partition each Si into

  Si,0 = {x in Si: |Si intersect fi^(-1)(fi(x))|<=p^(20eta)},
  Si,1 = Si \ Si,0.                                       (16)

These parts are unions of COMPLETE fibers. The argument of fi^(-1) is
fi(x), not x; this corrects the displayed fiber-partition typo in CZ.
Discard any part with size below p^(-eta)|Si|. There is at most one
such part per coordinate because p^(-eta)<1/2. The total probability
of a branch using any discarded part is at most 8*p^(-eta).
Every retained product branch has independent flat coordinates of
support size at least p^(1-2eta).

On the all-small-fiber branch, Section 6 gives a joint law Q with
largest atom p^(-(2+8eta)), at TV cost at most 36*p^(-eta).
The Q-mass of values v having Q_V(v)<p^(-(2+4eta)) is at most p^(-4eta),
since V has at most p^2 values. On every other v,
max_w Q(w|v)<=p^(-4eta). Replace the conditional W law on bad v by
uniform on F_p^2, keeping Q_V. This costs at most p^(-4eta), and gives
(5) on every v. Total branch cost is at most 37*p^(-eta).

If a retained branch has a large-fiber coordinate j, use the expression's
symmetries to put it first. Fix f1(X1) and X2,...,X8. The tampered output
V is then fixed. The conditional X1 law is uniform on a full fiber of
size greater than p^(20eta), because (16) retained complete fibers.
The first coordinate of the honest output is an affine function of X1
with coefficient X3+X4. On X3+X4!=0 this function is injective, giving
conditional W min-entropy at least 20eta*log p. The exceptional
probability is at most 1/|S4|<=p^(-1+2eta)<=p^(-eta), by independence.
Replace W on these exceptional fixed-data branches by uniform on F_p^2,
without changing V. Every fine branch now has conditional entropy at
least 20eta*log p. For each v, their mixture has the same maximum-atom
bound. Thus this large-fiber case gives (5), with TV cost p^(-eta).

On each discarded product branch replace W by uniform independently of
its actual V. The error is at most its probability. Convex mixtures of
joint laws satisfying max_w P(W=w|V=v)<=p^(-4eta) still satisfy that
inequality: multiply by P(V=v), sum the inequalities, then divide.
Adding the discarded-branch cost to the worst retained-branch cost
proves 45*p^(-eta), without paying for all 2^8 branches separately.

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
A coordinate map need only be fixed-point-free on that support; it can
be extended arbitrarily off support. No global no-fixed-point extension
is required by the theorem above.

Let inputs 9 and 10 be independent u-bit sources of entropy at least
(1-rho)u, independent of the first eight. Pad their concatenation with
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
is conditional only on E_p(gamma), (1), and the stated source guard.

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



<!-- COMPONENT constant-advice-breaker.md SELECTED SHA256 ba965970498b192b125e6f49cbc293fbd9ee4fd6dcd8d4517952601c3161940f -->

# A fully expanded constant-advice strong correlation breaker

Status: complete reduction to the finite seeded extractor proved in
seeded-linear.md, with an explicit numerical instantiation of that primitive.
This supplies an OLD, noncircular one-tampering standard breaker for the
constant-advice use inside the generalized affine-breaker transformation. It
uses neither Li 2023's new breaker nor a black-box flip-flop or merger. It is
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

Fix advice length v>=1 and M=2^v. Identify advice strings with 1,...,M in
lexicographic order. Fix positive integers n,d,h,a,s,s0,q and 0<eps<1 such
that 2^-q<=eps and s<=h. Supply these deterministic strong seeded extractors,
all of error at most eps:

  I : {0,1}^n x {0,1}^s0 -> {0,1}^h, entropy threshold KI;
  F : {0,1}^n x {0,1}^a  -> {0,1}^h, entropy threshold KF;
  G : {0,1}^d x {0,1}^s  -> {0,1}^a, entropy threshold KG;
  J : {0,1}^h x {0,1}^a  -> {0,1}^s, entropy threshold KJ.

Strong means the entire seed is retained. Assume s0<=d and

  h >= 2s+KJ+q,                                            (1)
  k >= max(KI, (2M+1)h+KF+q),                              (2)
  d >= 2s0+(3M+2)a+KG+q.                                  (3)

There is a deterministic map

  C : {0,1}^n x {0,1}^d x {0,1}^v -> {0,1}^h

such that for every arbitrary correlated pair (X,X') independent of every
arbitrary correlated pair (Y,Y'), with H_inf(X)>=k and Y exactly uniform,
and for every distinct advice alpha,alpha',

  D(C(X,Y,alpha) | C(X',Y',alpha'),Y,Y') <= (8M+2)eps.        (4)

Neither X' nor Y' needs entropy; neither must be a deterministic function
of its honest partner. Output may be truncated to any m<=h. Without the
unequal-advice requirement, the honest output alone is still close to
uniform with error at most (8M+2)eps, retaining Y,Y'.

The separate conditions k<=n and legal primitive dimensions are required.
Section 8 specifies actual primitives and satisfies every condition.

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
2s0+3ra. The stronger bound 2s0+(2r+1)a is also valid after the good
leadership round, but is unnecessary.

Before round j, maintain marginal closeness D(Wr|Hr)<=delta_r. At and after
round j, maintain joint closeness

  D(Wr | Hr,Wr')<=eta_r.                                  (11)

At any round the local extra left cost needed to prove extraction is at
most 3h, including a temporary fixing of the tampered next output. Before
round r this gives a maximum (2r+1)h<=(2M+1)h. Thus (2) supplies KF+q
residual source entropy in EVERY replenishment call in the proof.

The largest extra right cost needed for local G calls is 2a, beyond the
existing history. Condition (3) supplies KG+q after that leakage. It also
covers a post-leadership fixing of one a-bit tampered seed. Consequently
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

## 8. Explicit finite instantiation from seeded-linear.md

Here all numerical constants are deliberately wasteful. No existential
extractor constant remains. Fix integers n,m,e,v>=1 and M=2^v. We aim for
primitive error eps=2^-e and final error (8M+2)2^-e.

Use Section 1 of seeded-linear.md with internal height 11 and B0=2048.
For any source length u>=1 define

  b=e+15;
  rho_u=b+2 ceil(log(4u))+2;
  S_u=least power of 2 at least 45 rho_u;
  G_u=4096(4b+S_u+5);
  P_u=18S_u+18b+18.

Its exact construction for target z bits, with leaf output ceil(z/2048)
and final truncation, has entropy threshold and seed length bounded by

  K(u,z)<=2z+G_u,   D(u,z)<=6z/2048+P_u.                   (14)

These follow directly from the theorem's displayed bounds, with the
ceiling allowances included. In particular 128P_u<=G_u.

Choose w as the LEAST power of two satisfying

  w>=e,   128w>=m,
  w>=G_U, where U=max(n,128w,(3M+8)w).                      (15)

Define

  h=128w, a=s=s0=w,
  d=(3M+8)w,
  K=(2M+4)h.                                              (16)

The finite guard is K<=n. A source with entropy at least K suffices.
Instantiate I=F with the primitive for (input n, output h, error 2^-e);
G with (input d, output w, error 2^-e); J with (input h, output w,
error 2^-e). Each primitive uses its exact canonical algorithm, and ignores
an unused suffix to accept a seed of exactly w bits. Its strong guarantee
retains that suffix as well.

Monotonicity of G_u and (15) give G_n,G_d,G_h<=w and
P_n,P_d,P_h<=w/128. Therefore

  KI=KF<=257w;
  KG,KJ<=3w;
  D(n,h)<=49w/128<w;
  D(d,w),D(h,w)<=11w/1024<w.

Thus every padded seed fits. Each exact primitive's threshold is at most
its source length: for I,F use K<=n, for G use d>3w, and for J use h=128w.
Equation (1) holds since 2s+KJ+e<=6w<h. Equation (2) holds because

  (2M+1)h+KF+e <= (2M+1)h+258w <= (2M+4)h=K.

Equation (3) holds because

  2s0+(3M+2)a+KG+e <= (3M+8)w=d.

Run Section 2 and return the first m bits. The resulting map has all-seeds-
retained error at most (8M+2)2^-e and needs source entropy K. To obtain
error at most 2^-E, choose e=E+ceil(log(8M+2)).

The selector (15) halts by a literal dyadic search. For a finite bound,
S_u<90(e+21+2 ceil(log u)), whence

  G_u <= 2^23(e+ceil(log u)+1).

Set c=max(128,3M+8), and
F0=e+ceil(log n)+ceil(log c)+m+1>=2. Let wstar be the least power of two
at least 2^30 F0. Then wstar<2^31 F0 and
ceil(log max(n,c wstar))<=ceil(log n)+ceil(log c)+31+ceil(log F0).
Since ceil(log F0)<=F0 and 31<=16F0, the bound above gives
G_max(n,c wstar)<=2^28 F0<wstar. The other two guards hold too. Therefore
w<=wstar<2^31 F0. This proves an explicit finite stopping bound, and in
particular w=O_v(m+e+log n), rather than a noncomputable choice of constants.

For parameter choices failing K<=n, define the map to be the all-zero m-bit
string. This makes the family total while claiming extraction only under
its explicit guard. No search over distributions is performed.

The finite primitive is polynomial-time in its dimensions and error
exponent at fixed internal height 11. All dimensions in (16) and the
number of calls (1+4M) are bounded by explicit polynomials in n,m,e,M.
The dyadic selector uses O(log(F0)+log(M)+1) trials. Thus, for fixed v and
outer-logarithmic m,e, the evaluator is uniform polynomial bit-time in the
outer input length, with no unbounded-advice efficiency claim.

