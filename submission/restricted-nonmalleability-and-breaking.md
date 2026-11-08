# restricted-nonmalleability-and-breaking.md
Component references resolve through the map in sampler-and-disperser.md.
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
# Finite quadratic-field conditional entropy
Attribution: the encoding and original eight-input polynomial map are due to Eshan
Chattopadhyay and David Zuckerman, *Non-Malleable Codes Against Constant
Split-State Tampering*, ECCC TR14-102, Sections 6--7:
https://eccc.weizmann.ac.il/report/2014/102/download/ . Their expansion theorem
extends Jean Bourgain's work over the product of two prime fields. The graph
argument below is a self-contained quantitative Balog--Szemeredi--Gowers
argument; the minimal-growth argument is Petridis's short Plunnecke method.
The selected variant uses quadratic-field multiplication; its changed
expansion and encoding rigidity are proved explicitly below.
All logarithms are base two and TV is half L1. The nonlinear product below
is in Kp, while the joint ring Kp x Kp uses coordinatewise Kp-products.
U_A means uniform on a finite nonempty set A.
## 1. Canonical quadratic field and exact expansion input
For prime p>3 let D be the least positive nonsquare in F_p. There is one:
nonzero squaring has exactly two preimages per square, so only(p-1)/2
nonzero elements are squares. Testing all t in F_p for each candidate D
is a deterministic O(p^2 poly(logp)) algorithm using polynomial space.
Define Kp=F_p[omega]/(omega^2-D), represented by pairs(a,b), with
 (a,b)*(c,e)=(ac+Dbe,ae+bc), conjugation(a,b)=(a,-b).
For(a,b)!=0 its norm a^2-Db^2 is nonzero, since a zero norm with b!=0
would make D a square. Thus(a,-b)/(a^2-Db^2) is its inverse and Kp is a
field of size p^2. All operations and coordinate encodings are canonical.
Any unital F_p-field automorphism fixes F_p and sends omega to a root of
T^2-D, hence to omega or -omega. It is identity or conjugation.
This setup costs2^O(u) when p<2^(u+1), and is polynomial in actual outer
N at every selected u=O(logN); it is not an oracle or nonuniform table.
The fixed delta_R,q_R,gamma are defined in sampler SP.4; eta below may vary subject to (1). E_p(gamma) means: for A subset Kp x Kp,p^(2-h_R*eta0)<=|A|<=p^2, eta0=gamma/D_E, its unit
part C0 has at least half its size, both field projections>=p^q_R, and field fibers<=4p^((5/2)eta0);
every D0 subset A of size>=p^-delta_R*|A| has both projections>p and lies
in neither graph(z,c*z),(z,c*conj(z)). Then |A-A|+|C0/C0|>p^gamma*|A|.
For the selected gamma=6078385309659225236809123636178660293392119168499327/46544625331781931598897533495840739357496317520642048 and p>=2^(2^90), sampler SP.5 supplies a
complete proof of this exact statement, not an external premise.
## 2. Mixed-source conditional theorem and evaluator
Assume E_p(gamma), and choose a positive rational eta satisfying
 eta<=1/200,(4159/154+5*2^-80)*eta<=gamma,eta*logp>=64*2^80. (1)
Let X1,X2 in F_p and Z34,Z56 in F_p^2 be four independent sources with
 H_inf(Xi)>=log ceil(p^(1-eta/1232)), i=1,2,
 H_inf(Z34),H_inf(Z56)>=log ceil(p^(2-eta/616)). (2)
The coordinates inside either Z may be correlated. Let f1,f2 be arbitrary
maps on the scalar supports and g34,g56 arbitrary maps on the joint
supports into F_p^2. At least one of f1,f2 is fixed-point-free on its
honest support. Define q(x)=x^4+x^2+x,e(x)=x+q(x)*omega and h(x,y)=e(x)+e(y).
The evaluator on six scalar coordinates is
 W=h(X1,X2)*h(Z34)+h(Z56),
 V=h(f1(X1),f2(X2))*h(g34(Z34))+h(g56(Z56)). (3)
Their joint law is within zeta=min(1,45p^-eta) of a law(W*,V*) with
 H_inf(W*|V*=v)>=kappa=(5/4)eta*logp (5)
for EVERY v of positive V*-probability. The scalar inputs and two joint
inputs are the four independent groups; no independence within a joint
input is required. Evaluator(3) uses a constant number of field operations.
Smoothing the scalar inputs leaves support exponent1-(1233/1232)eta;
retained joint supports have exponents2-(617/616)eta and2-eta/616.
The three pair images have honest multiplicities at most4. The cutoff
p^(2-(1849/616)eta) discards at most10p^-eta from each image. Only the
first two images need unit deletion; the additive image is arbitrary.
Put a_set=1849/616+nu,tau=5/4,nu=2^-80. Section5's energy exponent is
4tau+4a_set+10=4159/154+4nu. The remaining nu pays fixed constants.
The joint atom cap in Section6 is2+(9/4)eta; Section7 converts it to(5).
Wrappers(U0') supplies eta*u>=4s+11sig+8, hence
kappa>=5s+(55/4)sig+10, paying its four-group conditioning loss.
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
### 3.2 Dense graphs give small self-differences
Let a bipartite graph on finite nonempty A,B have at least alpha*m*n
edges, m=|A|,n=|B|,0<alpha<=1. If every edge(a,b) has label a+b in T,
there is nonempty D subset A with
 |D|>=alpha*m/4, |D-D|<=2^14*|T|^4/(alpha^5*m*n^2). (7)
Call an ordered pair(x,z) in A^2 bad when it has fewer than alpha^2*n/64
common neighbors. Summing bad ordered pairs inside N(b) over b gives
less than alpha^2*m^2*n/64. Vertices b of degree below alpha*m/2
account for at most half the edges. The other vertices have degree sum
at least alpha*m*n/2 and squared-degree sum at least alpha^2*m^2*n/4,
by Cauchy. Thus some b0 of degree at least alpha*m/2 has at most
|N(b0)|^2/16 bad ordered pairs in N(b0). Put A0=N(b0).
Delete vertices with more than |A0|/4 bad partners in A0. At most
|A0|/4 are deleted, leaving D of size>=3|A0|/4>=alpha*m/4.
For any x,y in D, at least |A0|/2 vertices z in A0 are good with both.
Each gives at least(alpha^2*n/64)^2 walks x,b1,z,b2,y. There are thus
at least alpha^5*m*n^2/2^14 such four-edge walks, allowing repetitions.
Their labels satisfy x-y=t1-t2+t3-t4. For fixed x,y the labels recover
b1=t1-x,z=t2-b1,b2=t3-z, so different walks have different label tuples.
Tuples for distinct differences are disjoint; counting T^4 proves(7).
The same proof in a multiplicative abelian group gives
 |D|>=alpha*m/4, |D/D|<=2^14*|T|^4/(alpha^5*m*n^2). (8)
Here a/b=c/d iff ad=bc, and four-path labels give x/y=t1*t2^-1*t3*t4^-1.
When A,B are units of a commutative ring, every edge product is a unit;
only these labels are inverted, even if the ambient T contains nonunits.
## 4. Quadratic-field encoding rigidity
Suppose f,g have fibers<=L=p^((5/4)eta), with f fixed-point-free on its
honest support. Let A subset of
 {(e(x)+e(y),e(f(x))+e(g(y))):x,y in their F_p domains} subset Kp x Kp.
Every S subset A with |S|>=p^(2-h_R*eta) satisfies E_p(gamma)'s hypotheses.
Here eta<=eta0, so |S|>=p^(2-h_R*eta0). First |S|<=p^2. Any prescribed e(x)+e(y) has at most four preimages:
y=s-x from the real coordinate, then q(x)+q(s-x) has degree4 with
nonzero leading coefficient2. Hence for D0 subset A its honest projection
has size>=|D0|/4 and its tampered projection>=|D0|/(4L^2).
A zero honest Kp coordinate gives y=-x, at most p pairs. A zero tampered
coordinate gives g(y)=-f(x), at most pL pairs. Thus A has at most2pL<=10pL
nonunits. At eta<=1/200, 10pL/|S|<=10p^(-1+(h_R+5/4)eta)<=10p^-eta<1/2.
For a forbidden graph, fix y and write g(y)=b0,c=A0+B0*omega.
On the graph e(f(x))+e(b0)=c*(e(x)+e(y)), the real equation forces
 f(x)=A0*x+D*B0*q(x)+t.
If B0!=0, substitution into the imaginary equation gives a degree16
polynomial with nonzero leading coefficient(D*B0)^4, so at most16 roots.
If B0=0,A0=0, f is constant on solutions, at most L points. If B0=0,A0!=0,
the polynomial q(A0*x+t)-A0*q(x)+constant either has at most4 roots or is
identically zero. In the latter case its cubic coefficient4A0^3*t forces
t=0; quadratic and quartic coefficients give A0^2=A0,A0^4=A0, hence A0=1.
Then every solution has f(x)=x, forbidden. Each y therefore contributes<=L.
For the conjugate graph the real equation is f(x)=A0*x-D*B0*q(x)+t.
The B0!=0 case again has nonzero degree16. With B0=0,A0!=0, the polynomial
is q(A0*x+t)+A0*q(x)+constant. An identity would first force t=0, then
its linear coefficient2A0 is nonzero, impossible. The A0=0 case is again
a constant fiber. Both graph intersections thus have size<=pL; L>16 is
paid by p^eta>=2^64. No assumption that f is a polynomial is made: only
solutions of the real equation are substituted into the imaginary one.
For D0 subset S with |D0|>=p^-delta_R*|S|,
 |D0|>=p^(2-delta_R-h_R*eta)>pL.
The bounds delta_R<1/7,h_R<7,eta<=1/200 give
 1-(h_R+9/4)eta>=763/800>0,
 1-delta_R-(h_R+5/4)eta>=4569/5600>0,
proving the unit and graph exclusions. Since p^eta>=2^64>8, the unit
part C0 of S has min field projection>=|S|/(8L^2)
 >=p^(2-(h_R+7/2)eta)>=p^q_R.
Every qualifying D0 has min projection>=|D0|/(4L^2)
 >=p^(2-delta_R-(h_R+7/2)eta)>=p^(q_R-delta_R)>p.
The preimage bounds also give field fiber caps4 and4L^2<=4p^((5/2)eta0), inherited by C0. Both projection inequalities use eta<=gamma/D_E from(1), the identities
q_R=2-(h_R+7/2)*eta0 and q_R-delta_R=14gamma-eta0>1, eta0=gamma/D_E.
Thus SP.4 applies to every such S. Its expansion conclusion alone uses
E_p(gamma); the encoding, unit, graph and projection bounds were proved first.
## 5. Finite distributional energy growth
Put nu=2^-80,a_set=1849/616+nu,tau=5/4,alpha=p^-eta,T=4p^(2+tau*eta).
Let A1,A2 be unit sets in G=(Kp^*)^2 and C any subset of Kp x Kp, all with sizes in
[p^(2-a_set*eta),p^2]. Assume EVERY S subset A1 with |S|>=p^(2-h_R*eta)
has |S-S|+|S/S|>p^gamma|S|. Then EVERY W subset(Kp x Kp),|W|<=T,satisfies
 Pr[U_A1*U_A2+U_C in W]<=2alpha.
Call D plus-good if Pr[U_D+U_C0 in W]<=alpha for EVERY arbitrary C0 in that
size range and EVERY |W|<=T; define times-good using UNIT multiplicative partners B. Size-weighted
mixtures show disjoint unions of plus-good sets are plus-good.
Peel plus-good subsets from R=A1 until R is times-good or |R|<=alpha|A1|.
If neither holds, |R|>=p^(2-(a_set+1)eta), and a witness B,W has product
edge-density>alpha. The multiplicative four-path lemma gives D subset R:
 |D|>=p^(2-(a_set+2)eta)/4,
 |D/D|<=2^22*p^(2+(4tau+3a_set+6)eta). (13)
Indeed |W|^4<=2^8*p^(8+4tau*eta), and alpha^5|R||B|^2 has exponent
6-(3a_set+6)eta. All product edge labels are units.
If D is not plus-good, an additive witness C0,W and the same lemma give
D' subset D with
 |D'|>=p^(2-(a_set+3)eta)/16,
 |D'-D'|<=2^24*p^(2+(4tau+3a_set+7)eta). (14)
The factor4 from the lower bound on |D| enters this denominator.
Also D'/D' subset D/D; no sumset or product-set distributivity is used.
Divide both bounds by the LOWER bound on |D'|. Their relative sum is
 <=(2^28+2^26)*p^((4tau+4a_set+10)eta)
 <2^29*p^((4159/154+4nu)eta)<=p^((4159/154+5nu)eta)<=p^gamma.
Here nu*eta*logp>=64 pays2^29. Moreover
(h_R-a_set-3)eta*logp=nu*eta*logp>=64>4, so
|D'|>=p^(2-h_R*eta), contradicting hereditary expansion.
Thus D is plus-good, and the peeling always removes a nonempty set
until it terminates. Unit multiplication preserves plus-goodness: rescale
the arbitrary additive partner and target by the inverse unit; its size is preserved and no unit assumption on C0 is used. For a plus-good piece fix the
A2 coordinate; for the times-good remainder fix C and use W-c. Their
contributions are at most alpha. A leftover of size<=alpha|A1| has mixture
weight<=alpha. Averaging gives2alpha, with all partner/target quantifiers.
## 6. Small fibers for the mixed-source map
Suppose the four independent inputs are flat, with scalar support sizes
>=p^(1-(1233/1232)eta), joint Z34 support>=p^(2-(617/616)eta), and
joint Z56 support>=p^(2-eta/616). Assume f1,f2 have fibers<=L=p^((5/4)eta),
g34 has fibers<=4L on its support, and g56 is arbitrary. One of f1,f2
is fixed-point-free; exchange these two coordinates if necessary.
For each of the three pair images in Kp x Kp, the honest coordinate
h(x,y) has at most four preimages by Section4's quartic calculation.
Partition each image by exact multiplicity r=1,2,3,4. A component A_r has
weight r|A_r|/M and is uniform, where M is its input support size.
The three lower bounds on M have exponents2-(1233/616)eta,
2-(617/616)eta and2-eta/616 respectively. Discard image parts of size
<p^(2-(1849/616)eta). Their masses are at most10p^-eta,10p^(-2eta),
10p^(-3eta), hence total TV cost<=30p^-eta. Pair images stay independent
on their own component labels. There is no division by retained mass.
For the first image, zero honest coordinate occurs at at most p input
pairs, and zero tampered coordinate at at most pL. For the second,
zero honest coordinate occurs at at most p inputs; at most p possible
output pairs of g34 have h=0, with at most4L preimages each. Thus its
bad-input count is at most p+4pL<=5pL. Both images have at most10pL
nonunit points. Each retained component's nonunit fraction is at most
10p^(-1+(1849/616+5/4)eta)<=10p^(-1+7eta)<=p^-eta,
using1-9eta>0 and p^eta>=2^64. Replace the first two components by
uniform laws on their unit subsets, keeping the component weights.
This costs<=2p^-eta (also<=3p^-eta). Their sizes are at least
 (1-p^-eta)p^(2-(1849/616)eta)>=p^(2-(1849/616+nu)eta),
since p^(-nu*eta)<=2^-64<1/2<1-p^-eta. The third image needs no unit
restriction and already meets the same lower size bound. All sizes<=p^2.
The first image has Section4's scalar-pair rigidity. Thus Section5, with
arbitrary additive partner C, applies to the exact component law
U_A1*U_A2+U_C and gives mass<=eps_hit=2p^-eta on every target of size
<=T=4p^(2+(5/4)eta). Mixtures preserve this bound.
The guards imply2<=T<=p^4 and2+(9/4)eta<=4. The floor(T) heaviest
points have mass<=eps_hit. Outside them each atom is at most
2eps_hit/T=p^(-(2+(9/4)eta)). Clipping to this cap removes at most
eps_hit; ambient size p^4 permits redistribution below the cap. Its TV
cost is the removed mass. Total cost<=30p^-eta+3p^-eta+2p^-eta<36p^-eta.
The resulting joint law has min-entropy(2+(9/4)eta)logp. (15)
## 7. Arbitrary fibers and arbitrary mixed sources
First take independent flat inputs satisfying(2), with scalar supports
S1,S2 and joint supports S34,S56. Partition each scalar Si into its
fibers of size<=p^((5/4)eta) and those of larger size, measuring fibers
inside Si. Partition S34 using the threshold4p^((5/4)eta) for g34.
Each part is a union of COMPLETE fibers on the original support.
Do not partition g56. Discard parts of size<p^-eta times their original
support; at most one part per partitioned group is discarded, total
cost<=3p^-eta. Every retained product branch keeps the four groups
independent and flat. Its support sizes meet Section6's three bounds.
On the all-small branch, Section6 gives Q with largest joint atom
p^(-(2+(9/4)eta)) at cost<=36p^-eta. Q-mass on v with
Q_V(v)<p^(-(2+eta)) is at most p^-eta, since V has at most p^2 values.
On other v, max_w Q(w|v)<=p^(-(5/4)eta). Replace the honest conditional
law on bad v by uniform Kp, keeping Q_V. Its entropy is sufficient
because(5/4)eta<2. This costs<=p^-eta and proves(5) at cost<=37p^-eta.
If a scalar fiber is large, fix its partner value and all other honest
inputs. V is fixed and the honest W is affine in e(Xj), with coefficient
h(Z34). If nonzero this coefficient is invertible and e is injective,
so the complete fiber of size>p^((5/4)eta) gives(5). A zero coefficient
lies among at most p pairs of S34 and has probability
<=p/|S34|<=p^(-1+(617/616)eta)<=p^(-1+2eta)<=p^-eta.
On that event replace W by uniform Kp keeping V. Exchange X1,X2 to
handle the second scalar. No exchange with a joint group is needed.
If g34 has a large fiber, fix g34(Z34),X1,X2,Z56. V is fixed. Off the
zero event h(X1,X2)=0, W is an invertible affine function of h(Z34),
whose fibers have size<=4. A complete g34 fiber of size>4p^((5/4)eta)
therefore gives atom cap<p^(-(5/4)eta), with no constant entropy loss.
The zero event has probability<=1/|S2|<=p^(-1+2eta)<=p^-eta; replace W
there by uniform Kp. Fine-branch mixtures preserve the conditional cap.
Choose a large scalar first if one exists, otherwise the large joint
branch. Every retained branch has cost<=37p^-eta. The first-pair
fixed-point-free hypothesis persists on every restricted scalar support.
On discarded branches replace W by uniform independently of the actual
V. The error is their mass. Conditional caps survive mixtures: multiply
the cap inequality by each component's V marginal, sum, and divide.
Thus total error<=(3+37)p^-eta<=45p^-eta, without summing over all branches.
For general sources, a distribution of atom cap1/M is a convex mixture
of uniform M-point sets: an extreme point of the finite polytope
0<=p_x<=1/M,sum p_x=1 cannot have two fractional coordinates (perturb
one against the other), nor exactly one (integrality of M*sum p_x).
Use M=ceil(p^(1-eta/1232)) for each scalar and M=ceil(p^(2-eta/616))
for each joint group. Product mixtures preserve independence and the
scalar fixed-point-free property. The uniform error and conditional cap
are preserved by mixing, establishing the theorem for all inputs in(2).
## 8. Prime choice and exact mixed-source encoding
The canonical least prime above2^u satisfies2^u<p<2^(u+1) by the supplied
finite prime proof, with polynomial cost in actual N. Encode W in
Rrow=2(u+1) bits. Each r-bit scalar is mapped injectively by appending
the fixed a-bit advice, u=r+a. Split each 2r-bit joint block into two
r-bit halves and append the same advice to each half. This is an
injective map into F_p^2 and preserves all joint entropy.
For a scalar support, k>=(1-eta/1232)(u+1)+1 guarantees(2), since
logceil(p^(1-eta/1232))<=(1-eta/1232)logp+1. (17)
For a joint support, k2>=(2-eta/616)(u+1)+1 gives its bound in(2).
When raw block entropies are r-sig and2r-sig, the scalar guard suffices
for both:2r-sig>=2(r-sig)>=(2-eta/616)(u+1)+2.
Different appended advices make the two scalar partner maps fixed-point-
free on actual honest supports. The joint partner maps may be arbitrary.
All functions may be extended arbitrarily off their honest supports.
The separate guard u>=2^90 pays ring expansion. Also eta*u>=64*2^80
implies eta*logp>=64*2^80 and p^eta>=2^64, paying energy constants.
The four-group conditional-entropy match and five returned rows are
proved in affine-and-standard-wrappers.md Sections3.4--3.5.
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
## 8.2 Selective-transcript four-round multicopy theorem
Equation numbers in this subsection are local.
### 8.2.1. Statement
Fix the six public label codes, in label order 1 through 6:
1100, 1010, 1001, 0110, 0101, 0011.
The bit 1 means select B; the bit 0 means select A.
Let t>=1, c=t+1. The honest label differs from every other label, and each
nonhonest label has multiplicity at most 2. Necessarily t<=10.
All labels, dimensions, and maps below are fixed, not selected from private data.
Let n,d,h,a,s,s0,q be positive integers, 0<eps<1, 2^(-q)<=eps, s<=h,
s0<=d. Supply deterministic strong seeded extractors of error at most eps:
 I: n source bits, s0 seed bits, h output bits, threshold KI;
 F: n source bits, a seed bits, h output bits, threshold KF;
 G: d source bits, s seed bits, a output bits, threshold KG;
 J: h source bits, a seed bits, s output bits, threshold KJ.
"Strong" retains the complete seed of the padded map. Every primitive has its
own legal dimension/field/threshold hypotheses. Require
 h >= 7s+KJ+q,
 k >= max(KI, (2c+t)h+2cs+KF+q), k<=n,
 d >= c*s0+4c*a+KG+q. (1)
For every source tuple L=(X0,...,Xt) on n-bit strings and seed tuple
R=(Y0,...,Yt) on d-bit strings such that the TWO WHOLE TUPLES are independent,
Hinf(X0)>=k, and Y0 is uniform, run each copy with its assigned label:
 W_i = I(X_i, Prefix_s0(Y_i)).
 For r=1,...,4:
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
Each copy's code has two B's and two A's. At every prefix of the rounds,
the total left summary cost is therefore at most
 2c*h+2c*s. (5)
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
### 8.2.5. The universal seven-prefix bound
View each weight-two code as a two-element subset of {1,...,4}.
All four coordinate degrees are3. At the first honest B, S is empty.
At most6-3=3 other labels omit its coordinate, each with at most two
copies, so nu<=1+2*3=7. At the second honest B, the unseparated copies
contain the first coordinate and must omit the second to choose A.
Only the honest label contains both. Hence at most3-1=2 other labels
qualify and nu<=1+2*2=5<=7. Missing labels or copies only decrease these
bounds. Thus h>=7s+KJ+q supplies the local guard in(8).
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
 2c*h+2c*s+t*h = (2c+t)h+2cs. (16)
Its right messages do not change the conditional left law. Equation (1)
therefore gives Havg(X0|Pplus)>=KF+q. Initialization needed only KI.
For an honest B in round r, H1 costs at most
c*s0+c*(r-1)a+|S|a right bits. The local proof adds only nu new A values.
N and S are disjoint index sets among the c copies, so |S|+nu<=c. Hence
 c*s0+[c*(r-1)+|S|+nu]a <= c*s0+c*r*a <= c*s0+4c*a. (17)
This proves (8) for every local G call. The first local G has no greater
right cost; prefixes cost no right entropy. The auxiliary law keeps the
same H1 marginal and right conditional law, so this entropy statement
is valid there as well.
For an honest A call the right cost is at most
c*s0+[c*(r-1)+|S|]a <= c*s0+[c*(r-1)+t]a < c*s0+4c*a.
Thus its G call in (14) is legal. No unselected A_i of an old separated
B-choice copy is revealed. A newly separated A-choice value is charged
once, and an already known old separated A-choice value is free.
The h guard is Section8.2.5. The robust tail at each honest extraction
has mass at most eps; no tail is charged per tampered copy.
Starting delta0<=eps, (15) has A:delta+4eps, B:2delta+8eps.
An adjacent A followed by B costs4eps more than B followed by A.
Thus the two-A/two-B worst ordering is AABB:
eps ->5eps ->9eps ->26eps ->60eps<=76eps.
The six actual code coefficients are36,40,44,48,52,60.
Every different weight-two label has an honest-1/other-0 coordinate.
Thus every tampered copy is in S_4. Invariant(6) retains all other
outputs and all c original seed strings. Dropping H_4 proves(2), with
the conservative76eps bound. Output truncation contracts TV.
At t=10,c=11, (17),(16) and the local h guard are respectively
44a;32h+22s;7s+KJ+q.
## 8.3 Five-label corollary with five prefixes
Use the first five codes from8.2:1100,1010,1001,0110,0101, with at most
two copies of each nonhonest label. Put t=8,c=9; missing copies are allowed.
The entire proof of8.2 applies with guards
 h>=5s+KJ+q,
 k>=max(KI,26h+18s+KF+q),k<=n,
 d>=9s0+36a+KG+q.
Indeed the coordinate degrees are(3,3,2,2). The first honest B can only
occur at1 or2, both of degree3, so nu<=1+2*(5-3)=5. At the second B,
its first coordinate has degree3 and at most two other labels contain
it while omitting the second, so nu<=5 again. These are the only local
prefix counts changed from8.2. The transcript, auxiliary law, rectangularity,
whole-tuple independence and retained-seed proof is otherwise identical.
Substituting c=9,t=8 into its two-A and four-round ledgers gives26/18/36.
The actual codes have errors36,40,44,48,52 times eps; in particular the
conservative76eps bound retains all9 original seeds and8 other outputs.
