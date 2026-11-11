# restricted-nonmalleability-and-breaking.md
Component references resolve through the map in sampler-and-disperser.md.
# Canonical prime selection and exact field-to-bit encodings
The finite prime certificate below uses Erdos's central-binomial proof of
Bertrand--Chebyshev.
Primary source: P. Erdos,*Beweis eines Satzes von Tschebyschef*,Acta Litt.
Sci. Szeged 5,194--198 (1932):
https://combinatorica.hu/~p_erdos/1932-01.pdf .
All logarithms are base two. Integers are represented in binary. The actual
outer input length is N;the local block length u will be O(log N).
## 1. A completely finite prime interval proof
For real x>=1,let P(x) be the product of all primes at most x. We first
prove P(m)<4^m for every positive integer m,by strong induction.
For m=1,2 this is direct. If m=2r>=4 then m is composite and
P(m)=P(m-1)<4^(m-1)<4^m. If m=2r+1>=3,each prime in (r+1,2r+1]
divides binomial(2r+1,r). The two equal middle coefficients in that odd
binomial row show binomial(2r+1,r)<=2^(2r)=4^r. Thus
P(2r+1)<=P(r+1)*binomial(2r+1,r)<4^(r+1)*4^r=4^(2r+1).
Consequently P(x)<4^x for every real x>=1,using its integer floor.
We prove that for every integer n>=2 there is a prime p with n<p<2n.
Assume no such prime exists. For n>=512,consider B=binomial(2n,n).
The maximum binomial coefficient is B. The 2n-1 interior coefficients
sum to at most (2n-1)B,and the two endpoints sum to 2<=B. Hence
B>=4^n/(2n). (1)
For any prime p,its valuation in B is
sum_(j>=1) [floor(2n/p^j)-2floor(n/p^j)].
Each term is zero or one. For p<=sqrt(2n),the entire prime-power
contribution is thus at most 2n;there are at most sqrt(2n) such
primes. For p>sqrt(2n),its exponent is at most one. A prime in
(2n/3,n] has exponent zero,since floor(2n/p)=2 and floor(n/p)=1.
There is by assumption no prime in (n,2n),and 2n is composite. Thus
B <= (2n)^sqrt(2n) * P(2n/3)
< (2n)^sqrt(2n) * 4^(2n/3). (2)
Combining (1)--(2) gives 4^(n/3)<(2n)^(sqrt(2n)+1).
This is impossible for n>=512. To verify the comparison explicitly,put
x=sqrt(2n)>=32. The reverse strict inequality is equivalent to
x/log_2 x > 6*(1+1/x).
At x=32 it reads 32/5>99/16,equivalently 512>495. The left side is
strictly increasing for x>=32 (differentiate x/ln x);the right side is
strictly decreasing. This proves the contradiction for every n>=512.
For 2<=n<512 use the following increasing prime chain:
2,3,5,7,13,23,43,83,163,317,631.
Each term after the first is strictly smaller than twice its predecessor.
All entries are prime: trial division by primes at most sqrt(631)<26,
namely 2,3,5,7,11,13,17,19,23,verifies the nontrivial entries. Given n in
this range,take the first chain entry p greater than n;its predecessor
q is at most n,so p<2q<=2n. This completes the proof for all n>=2.
## 2. The exact prime algorithm and bit cost
For integer u>=1,set p(u) to the least prime strictly greater than 2^u.
Enumerate successive integers,and for each q test division by integers
2<=d<=floor(sqrt(q)). All comparisons may be performed as d*d<=q without
floating point. A candidate is prime exactly when none of those integers
divides it: every composite has a factor at most its square root.
§1 proves that this search terminates with
2^u<p(u)<2^(u+1),ceil(log p(u))=u+1. (3)
There are fewer than 2^u candidate integers,each with at most
2^((u+1)/2) trial divisors and O((u+1)^2) bit cost per division. A safe
bound is O(2^(3u/2+2)*(u+1)^2). If u<=c log N for one fixed c,this is
polynomial in the actual outer length N (not in u).
Use the standard representatives 0,...,p-1 for F_p. Addition,
multiplication and inversion by Euclid have polynomial bit cost in u;
all encodings below are fixed leading-zero binary encodings of those
representatives.
## 3. Source-rate loss from the prime interval
Embed an arbitrary u-bit source X into F_p by its integer value. This is
injective because p>2^u,so its min-entropy is unchanged. Since p>2^u,an error p^-c is at most 2^(-cu). Suppose
H_infinity(X)>=(1-gamma)*u,
0<=gamma<delta<1.
It meets a theorem requiring entropy at least (1-delta)log p whenever
(delta-gamma)*u >= 1-delta. (4)
Indeed log p<u+1 and the right side of (4) rearranges to
(1-gamma)u>=(1-delta)(u+1). Condition (4) uses the POST-partition entropy.
The bit length of one output field element is u+1,not u. A pair in F_p^2
thus uses exactly
nField=2*(u+1) (5)
bits under the fixed injective encoding. This encoding preserves
min-entropy and TV distance exactly on the image. Any polynomial output
law is regarded as a distribution on the full nField-bit cube supported
on that image.
# Finite quadratic-field conditional entropy
Attribution: the encoding and original eight-input polynomial map are due to Eshan
Chattopadhyay and David Zuckerman,*Non-Malleable Codes Against Constant
Split-State Tampering*,ECCC TR14-102,§§6--7:
https://eccc.weizmann.ac.il/report/2014/102/download/ . Their expansion theorem
extends Jean Bourgain's work over the product of two prime fields. The graph
argument below is a self-contained quantitative Balog--Szemeredi--Gowers
argument;the minimal-growth argument is Petridis's short Plunnecke method.
The selected variant uses quadratic-field multiplication;its changed
expansion and encoding rigidity are proved explicitly below.
All logarithms are base two and TV is half L1. The nonlinear product below
is in Kp,while the joint ring Kp x Kp uses coordinatewise Kp-products.
U_A means uniform on a finite nonempty set A.
## 1. Canonical quadratic field and exact expansion input
For prime p>3 let D be the least positive nonsquare in F_p. There is one:
nonzero squaring has exactly two preimages per square,so only(p-1)/2
nonzero elements are squares. Testing all t in F_p for each candidate D
is a deterministic O(p^2 poly(logp)) algorithm using polynomial space.
Define Kp=F_p[omega]/(omega^2-D),represented by pairs(a,b),with
(a,b)*(c,e)=(ac+Dbe,ae+bc),conjugation(a,b)=(a,-b).
For(a,b)!=0 its norm a^2-Db^2 is nonzero,since a zero norm with b!=0
would make D a square. Thus(a,-b)/(a^2-Db^2) is its inverse and Kp is a
field of size p^2. All operations and coordinate encodings are canonical.
Any unital F_p-field automorphism fixes F_p and sends omega to a root of
T^2-D,hence to omega or -omega. It is identity or conjugation.
This setup costs2^O(u) when p<2^(u+1),and is polynomial in actual outer
N at every selected u=O(logN);it is not an oracle or nonuniform table.
The fixed delta_R,q_R,gamma are defined in sampler SP.4;eta below may vary subject to (1). E_p(gamma) means: for A subset Kp x Kp,p^(2-h_R*eta0)<=|A|<=p^2,eta0=gamma/D_E,its unit
part C0 has at least half its size,both field projections>=p^q_R,and field fibers<=4p^((5/2)eta0);
every D0 subset A of size>=p^-delta_R*|A| has both projections>p and lies
in neither graph(z,c*z),(z,c*conj(z)). Then |A-A|+|C0/C0|>p^gamma*|A|.
For the selected gamma=((2-nu)/(14+(h_R+5/2)/D_E)) and p>=2^(2^90),sampler SP.5 supplies a
complete proof of this exact statement,not an external premise.
## 2. Mixed-source conditional theorem and evaluator
Assume E_p(gamma),and choose a positive rational eta satisfying
eta<=1/120,(16089/937+5*2^-80)*eta<=gamma,eta*logp>=64*2^80. (1)
Let X1,X2 in F_p and Z34,Z56 in F_p^2 be four independent sources with
H_inf(Xi)>=log ceil(p^(1-eta/(937/20))),i=1,2,
H_inf(Z34)>=log ceil(p^(2-(1017/1874)eta)),
H_inf(Z56)>=log ceil(p^(2-2eta/(937/20))). (2)
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
inputs are the four independent groups;no independence within a joint
input is required. Evaluator(3) uses a constant number of field operations.
Partition complete fibers without discarding small branches. Section6 weights
image losses against the ORIGINAL source supports;the first/third image cutoffs
are p^(2-(977/937)eta),and the second is p^(2-(2891/1874)eta).
Each costs at most10p^-eta globally against its original support.
Only the first two images need unit deletion;the additive image is arbitrary.
Put a_set=c_set=977/937+nu,b_set=2891/1874+nu,tau=5/4,nu=2^-80.
The maximum relative energy exponent is4tau+4a_set+8=16089/937+4nu.
The remaining nu pays fixed constants. Section7 combines the weighted
small branch with all complete large fibers,without support-size conditioning.
The joint atom cap in Section6 is2+(9/4)eta;Section7 converts it to(5).
Wrappers(U0') supplies eta*u>=4s+11sig+8,hence
kappa>=5s+(55/4)sig+10,paying its four-group conditioning loss.
## 3. An elementary additive-combinatorial tool
### 3.2 Dense graphs give small self-differences
Let a bipartite graph on finite nonempty A,B have at least alpha*m*n
edges,m=|A|,n=|B|,0<alpha<=1. If every edge(a,b) has label a+b in T,
there is nonempty D subset A with
|D|>=alpha*m/4,|D-D|<=2^14*|T|^4/(alpha^5*m*n^2). (7)
Call an ordered pair(x,z) in A^2 bad when it has fewer than alpha^2*n/64
common neighbors. Summing bad ordered pairs inside N(b) over b gives
less than alpha^2*m^2*n/64. Vertices b of degree below alpha*m/2
account for at most half the edges. The other vertices have degree sum
at least alpha*m*n/2 and squared-degree sum at least alpha^2*m^2*n/4,
by Cauchy. Thus some b0 of degree at least alpha*m/2 has at most
|N(b0)|^2/16 bad ordered pairs in N(b0). Put A0=N(b0).
Delete vertices with more than |A0|/4 bad partners in A0. At most
|A0|/4 are deleted,leaving D of size>=3|A0|/4>=alpha*m/4.
For any x,y in D,at least |A0|/2 vertices z in A0 are good with both.
Each gives at least(alpha^2*n/64)^2 walks x,b1,z,b2,y. There are thus
at least alpha^5*m*n^2/2^14 such four-edge walks,allowing repetitions.
Their labels satisfy x-y=t1-t2+t3-t4. For fixed x,y the labels recover
b1=t1-x,z=t2-b1,b2=t3-z,so different walks have different label tuples.
Tuples for distinct differences are disjoint;counting T^4 proves(7).
The same proof in a multiplicative abelian group gives
|D|>=alpha*m/4,|D/D|<=2^14*|T|^4/(alpha^5*m*n^2). (8)
Here a/b=c/d iff ad=bc,and four-path labels give x/y=t1*t2^-1*t3*t4^-1.
When A,B are units of a commutative ring,every edge product is a unit;
only these labels are inverted,even if the ambient T contains nonunits.
## 4. Quadratic-field encoding rigidity
Suppose f,g have fibers<=L=p^((5/4)eta),with f fixed-point-free on its
honest support. Let A subset of
{(e(x)+e(y),e(f(x))+e(g(y))):x,y in their F_p domains} subset Kp x Kp.
Every S subset A with |S|>=p^(2-h_R*eta) satisfies E_p(gamma)'s hypotheses.
Here eta<=eta0,so |S|>=p^(2-h_R*eta0). First |S|<=p^2. Any prescribed e(x)+e(y) has at most four preimages:
y=s-x from the real coordinate,then q(x)+q(s-x) has degree4 with
nonzero leading coefficient2. Hence for D0 subset A its honest projection
has size>=|D0|/4 and its tampered projection>=|D0|/(4L^2).
A zero honest Kp coordinate gives y=-x,at most p pairs. A zero tampered
coordinate gives g(y)=-f(x),at most pL pairs. Thus A has at most2pL<=10pL
nonunits. At eta<=1/120,10pL/|S|<=10p^(-1+(h_R+5/4)eta)<=10p^-eta<1/2.
For a forbidden graph,fix y and write g(y)=b0,c=A0+B0*omega.
On the graph e(f(x))+e(b0)=c*(e(x)+e(y)),the real equation forces
f(x)=A0*x+D*B0*q(x)+t.
If B0!=0,substitution into the imaginary equation gives a degree16
polynomial with nonzero leading coefficient(D*B0)^4,so at most16 roots.
If B0=0,A0=0,f is constant on solutions,at most L points. If B0=0,A0!=0,
the polynomial q(A0*x+t)-A0*q(x)+constant either has at most4 roots or is
identically zero. In the latter case its cubic coefficient4A0^3*t forces
t=0;quadratic and quartic coefficients give A0^2=A0,A0^4=A0,hence A0=1.
Then every solution has f(x)=x,forbidden. Each y thus contributes<=L.
For the conjugate graph the real equation is f(x)=A0*x-D*B0*q(x)+t.
The B0!=0 case again has nonzero degree16. With B0=0,A0!=0,the polynomial
is q(A0*x+t)+A0*q(x)+constant. An identity would first force t=0,then
its linear coefficient2A0 is nonzero,impossible. The A0=0 case is again
a constant fiber. Both graph intersections thus have size<=pL;L>16 is
paid by p^eta>=2^64. No assumption that f is a polynomial is made: only
solutions of the real equation are substituted into the imaginary one.
For D0 subset S with |D0|>=p^-delta_R*|S|,
|D0|>=p^(2-delta_R-h_R*eta)>pL.
The bounds delta_R<1/7,h_R<7,eta<=1/120 give
1-(h_R+9/4)eta>=443/480>0,
1-delta_R-(h_R+5/4)eta>=883/1120>0,
proving the unit and graph exclusions. Since p^eta>=2^64>8,the unit
part C0 of S has min field projection>=|S|/(8L^2)
>=p^(2-(h_R+7/2)eta)>=p^q_R.
Every qualifying D0 has min projection>=|D0|/(4L^2)
>=p^(2-delta_R-(h_R+7/2)eta)>=p^(q_R-delta_R)>p.
The preimage bounds also give field fiber caps4 and4L^2<=4p^((5/2)eta0),inherited by C0. Both projection inequalities use eta<=gamma/D_E from(1),the identities
q_R=2-(h_R+7/2)*eta0 and q_R-delta_R=13gamma-eta0>1,eta0=gamma/D_E.
Thus SP.4 applies to every such S. Its expansion conclusion alone uses
E_p(gamma);the encoding,unit,graph and projection bounds were proved first.
## 5. Finite distributional energy growth
Put nu=2^-80,a_set=c_set=977/937+nu,b_set=2891/1874+nu,tau=5/4,
alpha=p^-eta,T=4p^(2+tau*eta). Let A1,A2 be unit sets in(Kp^*)^2,
and C arbitrary in Kp x Kp,with respective size ranges
[p^(2-a_set*eta),p^2],[p^(2-b_set*eta),p^2],[p^(2-c_set*eta),p^2].
Assume EVERY S subset A1 of size>=p^(2-h_R*eta) has
|S-S|+|S/S|>p^gamma|S|. For EVERY target W of size<=T we prove
Pr[U_A1*U_A2+U_C in W]<=2alpha.
Call nonempty D subset A1 plus-good if Pr[U_D+U_C0 in W]<=alpha for
EVERY arbitrary C0 in the C size range and EVERY target of size<=T.
For nonempty R subset A1 put w_R=|R|/|A1|. Call R times-good if
Pr[U_R*U_B in W]<=alpha/w_R for EVERY unit B in the A2 size range
and EVERY target of size<=T. Empty R and w_R<=alpha are times-good.
Peel plus-good subsets from R=A1 until R is times-good. If it is not,
rho=alpha/w_R lies in(0,1),and a multiplicative witness has density>rho.
The four-path lemma gives nonempty D subset R with
|D|>=rho|R|/4=alpha|A1|/4>=p^(2-(a_set+1)eta)/4,
|D/D|<=2^14*T^4/(rho^5*|R|*|B|^2)
<=2^22*p^(2+(4tau+a_set+2b_set+5)eta). (13)
Indeed rho^5|R|=alpha^5|A1|/w_R^4>=alpha^5|A1|. All inverted edge labels are units.
If D is not plus-good,its additive witness has density>alpha. The same
lemma produces Dprime subset D with
|Dprime|>=alpha|D|/4>=p^(2-(a_set+2)eta)/16,
|Dprime-Dprime|<=2^24*p^(2+(4tau+a_set+2c_set+6)eta). (14)
The factor4 enters through the lower bound on |D|. Also
Dprime/Dprime subset D/D. Divide by the LOWER bound on |Dprime|.
The relative quotient and difference exponents are respectively
4tau+2a_set+2b_set+7=16089/937+4nu and
4tau+2a_set+2c_set+8=16089/937+4nu. Their relative sum is
<=(2^26+2^28)*p^((16089/937+4nu)eta)
<2^29*p^((16089/937+4nu)eta)<=p^((16089/937+5nu)eta)<=p^gamma.
Here nu*eta*logp>=64 pays2^29. Further h_R-a_set-2=nu,so that
same margin pays16 and |Dprime|>=p^(2-h_R*eta). This contradicts the
hereditary hypothesis. Thus D is plus-good and can be removed.
Each removal is nonempty,so the finite process terminates.
For removed D fix b in A2. As b is a unit,scaling C,W by b^-1
preserves their cardinalities and C's class;plus-goodness bounds its
conditional hit probability byalpha. The disjoint original weights
|D|/|A1| sum to<=1,so all removed pieces contribute<=alpha.
For remaining R fix c in C and use W-c. Times-goodness bounds its
conditional hit probability byalpha/w_R;original weight w_R gives
contribution<=alpha (zero if R is empty). This proves2alpha.
## 6. Weighted small-fiber smoothing
Take independent original flat supports S1,S2,S34,S56 satisfying(2),of
sizes M1,M2,M34,M56. Let T1 subset S1,T2 subset S2,T34 subset S34 be
nonempty sets on which f1,f2 have fibers<=L=p^((5/4)eta),and g34 has
fibers<=4L. Keep S56 whole and its partner arbitrary. One scalar map
is fixed-point-free. Write N1=|T1||T2|,N2=|T34| and
w=N1*N2/(M1*M2*M34)<=1.
For the normalized law P_T on these four independent uniform supports,
we construct Q_T with the conditional cap(5) and
w*TV(P_T,Q_T)<=36p^-eta. (15)
No lower bound on |T1|,|T2| or |T34| is required.
Each of the three pair images has honest multiplicity at most4 by
Section4. Partition each by its exact multiplicity r=1,2,3,4. An image
part A_r is uniform with weight r|A_r|/N,where N=N1,N2,M56 respectively.
Put B1=B3=p^(2-(977/937)eta),B2=p^(2-(2891/1874)eta).
For image j discard parts with |A_r|<Bj in this analysis.
For the first image their weighted mass is at most
w*sum_r r|A_r|/N1<=10B1/(M1*M2)<=10p^-eta.
For the second it is <=10B2/M34<=10p^-eta,and for the third
<=10B3/M56<=10p^-eta. These inequalities use the ORIGINAL bounds
M1*M2,M56>=p^(2-2eta/(937/20)) and M34>=p^(2-(1017/1874)eta);
each cutoff exponent exceeds its original deficit by1. Other weights are<=1.
The union of discarded component events thus has weighted mass<=30p^-eta.
On each discarded event replace the FULL output pair(W,V) by uniform
Kp x Kp,retaining its component weight. This is permitted: the theorem
does not require Q_T's partner marginal to equal the real marginal.
Uniform Kp x Kp has joint atom cap p^(-4),which is sufficient below.
On the retained product components,the first image has at most2pL
nonunit points. The second has at most p+4pL<=5pL: h=0 permits at
most p input pairs or output pairs,and g34 has fibers<=4L. Thus each
has at most10pL nonunits. Their nonunit fractions are bounded by
10p^(-1+(2891/1874+5/4)eta)<=10p^(-1+7eta)<=p^-eta,
using1-9eta>0 and p^eta>=2^64. Replace their uniform laws by uniform
laws on their unit subsets,preserving the component weights. Weighted
TV cost is <=2p^-eta,also <=3p^-eta. Their respective sizes are at least(1-p^-eta)B1>=p^(2-a_set*eta)
and(1-p^-eta)B2>=p^(2-b_set*eta),since p^(-nu*eta)<=2^-64<1/2.
The third image already has size>=B3>=p^(2-c_set*eta) and needs no units.
The first component has Section4's scalar-pair rigidity on T1,T2.
Section5 thus applies to every retained product component,with
arbitrary additive partner,giving mass<=2p^-eta on every target of
size<=T=4p^(2+(5/4)eta). Each normalized component can be clipped to
joint atom cap p^(-(2+(9/4)eta)) at cost<=2p^-eta: outside its floor(T)
heaviest atoms the cap follows from2eps_hit/T,and excess removed from
those atoms is<=eps_hit=2p^-eta. Ambient size p^4 permits redistribution.
This weighted cost is also <=2p^-eta. The uniform replacements on
discarded components have the same cap,since2+(9/4)eta<=4.
Mixtures preserve a joint atom cap,so the resulting normalized Q has
that cap and w*TV(P_T,Q)<=35p^-eta. Q-mass on v with
Q_V(v)<p^(-(2+eta)) is <=p^-eta,since V has at most p^2 values.
On the other v,Q(W=w|V=v)<=p^(-(5/4)eta). Replace W by uniform Kp
on the bad v,keeping Q_V. Its entropy suffices because(5/4)eta<2.
This costs at most p^-eta in the normalized law and at most p^-eta
after multiplication by w. This proves(15),even if every component
was discarded. If w=0 there is no small branch to construct.
## 7. Complete fibers without branch discard
Start with the original independent flat supports in(2). Partition each
scalar support into complete fibers of size<=L and those of size>L;
partition S34 at the threshold4L for g34. Do not partition S56.
Do not discard ANY partition part. Empty parts simply have zero weight.
The resulting at most eight product branches are independent and flat.
On the all-small branch use Section6: its probability is exactly w,
and its weighted error is at most36p^-eta.
On each remaining branch choose the first large scalar fiber,if any,
and otherwise a large g34 fiber. For a large scalar Xj,fix its partner
value and all other honest variables. V is fixed and W is affine in
e(Xj),with coefficient h(Z34). Off coefficient zero,invertibility and
injectivity of e give atom cap<1/L on the COMPLETE original fiber.
No lower bound on the partition-part probability was used.
For a large g34 fiber fix g34(Z34),X1,X2,Z56. Off h(X1,X2)=0,W is
an invertible affine function of h(Z34),whose fibers have size<=4.
A complete g34 fiber has size>4L,so again the cap is<1/L.
On a zero-coefficient event replace W by uniform Kp,keeping the actual
V. Bound these errors with their original weights,not separately by
inverse branch probabilities. Across scalar-large branches their sum
is <=Pr[h(Z34)=0]<=p/M34,since at most p original pairs have h=0.
Across g34-large branches their sum is <=Pr[h(X1,X2)=0]
<=1/min(M1,M2). Branches are disjoint by the chosen priority;summing
their conditional event probabilities times their weights cancels all
partition normalizers. Each of these two original probabilities is
<=p^-eta by(2) and eta<=1/120. Total weighted large-branch error is
thus<=2p^-eta. Fine-branch conditional caps survive mixtures:
Q(w,v)<=p^(-(5/4)eta)Q_V(v) can be summed before dividing by Q_V(v).
Mix all small and large ideal laws with their original branch weights.
The same conditional cap holds,and total error is at most
(36+2)p^-eta<=45p^-eta. This proves(5) without branch support cutoffs.
For general sources,a distribution with atom cap1/M is a convex mixture
of uniform M-point sets: an extreme point of0<=p_x<=1/M,sum p_x=1
cannot have two fractional entries (perturb one against the other),or
exactly one (integrality of M*sum p_x). Use M=ceil(p^(1-eta/(937/20))) for
each scalar M=ceil(p^(2-(1017/1874)eta)) for Z34,and M=ceil(p^(2-2eta/(937/20)))
for Z56. Independent
product mixtures preserve all four-group independence and the scalar
fixed-point-free property. Uniform error and conditional caps survive
mixing,proving the theorem for every source in(2).
## 8. Prime choice and exact mixed-source encoding
The canonical least prime above2^u satisfies2^u<p<2^(u+1) by the supplied
finite prime proof,with polynomial cost in actual N. Encode W in
Rrow=2(u+1) bits. Each r-bit scalar is mapped injectively by appending
the fixed a-bit advice,u=r+a. Split the third raw block into r,r-m bits,m=floor(eta*r/2),
and zero-pad its second part to r;split the fourth into r,r and append
the same advice to all four coordinates. This is an
injective map into F_p^2 and preserves all joint entropy.
For a scalar support,k>=(1-eta/(937/20))(u+1)+1 guarantees(2),since
logceil(p^(1-eta/(937/20)))<=(1-eta/(937/20))logp+1. (17)
For a joint support,k2>=(2-2eta/(937/20))(u+1)+1 gives its bound in(2).
The full fourth block meets this bound since2r-sig>=2(r-sig).
The third needs k2>=(2-(1017/1874)eta)(u+1)+1;the scalar guard gives
2r-sig>=(2-2eta/(937/20))(u+1)+2,and m<=eta*r/2<=eta*(u+1)/2
pays its shorter input. These injections preserve all source entropy.
Different appended advices make the two scalar partner maps fixed-point-
free on actual honest supports. The joint partner maps may be arbitrary.
All functions may be extended arbitrarily off their honest supports.
The separate guard u>=2^90 pays ring expansion. Also eta*u>=64*2^80
implies eta*logp>=64*2^80 and p^eta>=2^64,paying energy constants.
The four-group conditional-entropy match and five returned rows are
proved in affine-and-standard-wrappers.md Sections3.4--3.5.
# Conditional tools and the selected antichain multicopy breaker
The local A/Z/B operations adapt Cohen's local correlation breakers and the
whole-tuple MC-LA proof of publishedc78568c1;all proofs are supplied below.
All logarithms are binary and TV is half L1. WriteD(V|T)=TV((V,T),(U,T)),
where U is independent uniform on V's length. Define
Havg_inf(X|T)=-log(sum_t Pr[T=t] max_x Pr[X=x|T=t]).
## 8. Selected antichain multicopy instance
The ordinary-G reference proof is8.2;selected NB uses the two-sided-affine variant8.6.
## 8.2 Selective-transcript four-round multicopy theorem
Equation numbers in this subsection are local.
### 8.2.1. Statement
Fix the five label codes,in order:
1100,1010,1001,0110,0101.
The bit 1 means select B;the bit 0 means select A.
Let t>=1,c=t+1. The honest label differs from every other label,and each
nonhonest label has multiplicity at most 2. Necessarily t<=8.
All labels,dimensions,and maps below are fixed,not selected from private data.
Let n,d,h,a,s,s0,q be positive integers,0<eps<1,2^(-q)<=eps,s<=h,
s0<=d. Supply deterministic source-linear seeded maps of error at most eps. I/F/J are strong on uniform affine sources;G is strong on arbitrary sources:
I: n source bits,s0 seed bits,h output bits,threshold KI;
F: n source bits,a seed bits,h output bits,threshold KF;
G: d source bits,s seed bits,a output bits,threshold KG;
J: h source bits,a seed bits,s output bits,threshold KJ.
"Strong" retains the complete seed of the padded map. Every primitive has its
own legal dimension/field/threshold hypotheses. Require
h >= 5s+KJ+q,
k >= max(KI,(2c+t)h+2cs+KF+q),k<=n,
d >= c*s0+4c*a+KG+q. (1)
For every tuple L=(X0,...,Xt) of fixed affine images of one uniform affine variable Z,and seed tuple R=(Y0,...,Yt) on d-bit strings such that the TWO WHOLE TUPLES are independent,
Hinf(X0)>=k,and Y0 is uniform,run each copy with its assigned label:
W_i = I(X_i,Prefix_s0(Y_i)).
For r=1,...,4:
T_i=Prefix_s(W_i);
A_i=G(Y_i,T_i);
Z_i=J(W_i,A_i);
B_i=G(Y_i,Z_i);
C_i=B_i for code bit 1,C_i=A_i for code bit 0;
W_i=F(X_i,C_i).
Then the final honest h-bit output satisfies
TV((W_0,W_1,...,W_t,R),(U_h,W_1,...,W_t,R)) <= 76eps. (2)
Here U_h is independent of the displayed original target marginal.
Taking any requested output prefix preserves the guarantee. All ignored suffixes of all c original seeds remain in R.
There is no entropy requirement on X_i or Y_i for i>0,no within-side
independence,and no functional relation to X0 or Y0. In particular sources
or seeds can coincide.
### 8.2.2. Conditional tools,including retained-side hypotheses
Write D(V|T)=TV((V,T),(U,T)),with the ideal law retaining the actual T
marginal. All variables have finite alphabets;zero-probability histories
can be omitted. Define
2^(-Havg(X|H)) = sum_h P(h) max_x P(X=x|H=h).
A rectangular history alternates functions of one whole side and its past.
Its conditional density factors: each message restricts only its own factor.
Independent private coins can be included on their respective side. For each
old history h,a b-bit left message obeys
sum_m max_x P(X0=x,M=m|h)<=2^b max_x P(X0=x|h).
Right messages leave X0's conditional law unchanged. Induction gives
Havg(X0|H)>=Hinf(X0)-bL,and symmetrically Havg(Y0|H)>=d-bR.
Messages may depend on correlated coordinates or previous messages;a value
already determined costs zero. Markov on the largest conditional atom gives
Hinf(X|H=h)>=K outside mass2^-q<=eps whenHavg(X|H)>=K+q.
Robust extraction rule. Suppose L,R are independent given rectangular H,
S is a right function given H with D(S|H)<=eta,X is a left coordinate,
Havg(X|H)>=K+q,and E is strong (K,eps). Then
D(E(X,S)|H,R) <= eta+2eps. (3)
For fixed H=h and seed s,let
g_h(s)=TV(E(X|h,s),U). Conditional independence means that retaining the
ENTIRE right side gives the average of g_h(S),since conditional source
laws depend on R only through h. Replacing S's marginal by uniform changes
this expectation by at most TV(S|h,U). On good-entropy h the uniform-seed
average is at most eps;the exceptional h have total mass at most eps.
Averaging proves (3).
Exchanging the sides gives the symmetric rule,retaining all of L.
Two additional facts will be used:
(a) If V is right-sided and Q is left-sided given rectangular H,then
D(V|H,Q)=D(V|H). Revealing Q preserves product and the right conditional law.
(b) If two joint laws on (V,T) are within eta,then
D_P(V|T)<=D_Q(V|T)+2eta. One eta compares the joint laws,the second compares
their T marginals after adjoining independent uniform V.
Projection of retained side information contracts D. Projection of a
rectangular transcript need NOT remain rectangular. Every compressed real
history below is separately proved rectangular by its message order.
### 8.2.3. The kept real history and separated-set invariant
Reveal all c initial s0-bit seed prefixes as H0. This is right-only.
Strong I with the honest uniform prefix gives
D(W0|H0,R)<=eps. (4)
For completeness,conditional output laws of W0 depend on R only through
the honest prefix because X0 is independent of the WHOLE R. Appending the
rest of R using its conditional kernel given that prefix thus keeps
the strong bound. All current W_i are left
functions given H0.
In each round define the actual selected summary Q_i to be:
- Prefix_s(W_i) for an A-choice copy;
- the complete h-bit W_i for a B-choice copy.
The kept real history H_r consists of H0 and,for each completed round,
all c selected summaries followed by all c actual selected C_i values.
Given old H,all current W_i are left functions. Thus Q is a left message.
After Q is fixed,every selected C_i is a right function:
an A-choice needs only its fixed prefix;a B-choice has its full state
fixed,so A_i,Z_i,B_i are functions of R. Revealing selected C is right-only.
Consequently H_r is genuinely rectangular. Each next W_i=F(X_i,C_i) is
left-sided on H_r. Unselected A,Z,B values are evaluated by the algorithm
but are not revealed in this history.
Affine real-history invariant. Given every real kept history,the latent Z is uniform on an affine fiber. Initially this is true;I is source-linear at its fixed seed prefix. Each old W_i=F(X_i,C_i) is affine in Z because its previous selected C_i is fixed in H_r. Selected summaries,prefixes,temporary W_S and next separated outputs are therefore affine messages once their seeds are fixed. Right-only messages leave the conditional left law unchanged. Each b-bit left message loses at most b rank on every nonempty fiber. Thus every REAL I/F call sees a uniform affine marginal of the required rank. The same reasoning applies to the shortened Fz. The auxiliary affine extension below justifies J;G remains ordinary.

The right cost after r rounds is at most c*s0+c*r*a.
Each copy's code has two B's and two A's. At every prefix of the rounds,
the total left summary cost is thus at most
2c*h+2c*s. (5)
Let S_r be the tampered indices that have chosen A in at least one
completed round in which the honest copy chose B. This set depends only
on labels. Maintain
D(W0 | H_r,W_(S_r),R) <= delta_r. (6)
Initially S_0 is empty and delta_0<=eps by (4).
Before a round,let H=H_(r-1),S=S_(r-1). Temporarily reveal W_S on the
left,then current actual selected C_S on the right;call this H1.
C_S is right-sided after W_S is fixed,regardless of the choices of those
copies. Equation (6),which already retained W_S and all R,implies
D(W0|H1,R)<=delta. (7)
The real whole sides remain independent given H1. Its right cost is at
most c*s0+c*(r-1)*a+|S|a. Whole W_S is TEMPORARY,not kept in H_r.
### 8.2.4. Selective local look-ahead,including the auxiliary law
Suppose the honest choice is B. Let
Jset={0} union {tampered indices currently choosing A},and let
N={0} union {tampered indices outside S currently choosing A},
nu=|N|.
All old separated A-choice copies have their W_i and C_i=A_i fixed in H1.
Assume,for this local step,
Havg(Y0|H1)>=KG+nu*a+q,h>=nu*s+KJ+q. (8)
We prove the REAL-law bound
D(B0|H1,W0,...,Wt,(A_i:i in Jset),Z0)
<= 2delta+6eps. (9)
Given H1,the whole current-state tuple and R are independent. Replace
the honest STATE W0 by exact uniform given H1,keeping H1's
marginal and the conditional kernel
Law((W_i:i>0)|W0,H1),and keeping Law(R|H1) independently.
Use this affine off-support extension on each nonempty H1 fiber. The real state tuple is uniform on an affine space A with origin(x0,v0),where x denotes W0 and v the other states. Let P be its projected direction on x and K={v:(0,v) lies in direction(A)}. Choose a linear section s:P->V so direction(A)={(x,s(x)+k):x in P,k in K};extend s linearly to the whole h-bit x space. All W_S coordinates of K and s vanish since W_S is fixed;choose the extension with those coordinates zero. Set v=v0+s(x-x0)+U_K with independent uniform U_K and x uniform on the whole h-bit space. On the original x support this is exactly the real conditional kernel;elsewhere it keeps W_S fixed. Thus the new WHOLE state tuple is uniform affine and its joint distance from the real law is exactly the marginal replacement cost<=delta. Keep the right law independently. Old C_S equations survive. No original X equations off support are required.
In the auxiliary law W0 is uniform given H1,and the state tuple is
independent of R given H1. Compute the nu new prefixes T_i,i in N,
and write M for this prefix vector. Old separated prefixes are known.
The first G extraction is justified BEFORE conditioning on M: its honest
T0 is exactly uniform given H1. For fixed H1 the conditional A0 law given
the entire state tuple depends on it only through T0. On good Y-entropy
branches,strong G thus gives error eps even retaining the whole state
tuple;bad branches cost eps. Thus
D_aux(A0|H1,whole state tuple)<=2eps,
hence D_aux(A0|H1,M)<=2eps. (10)
On H2=(H1,M),the state tuple and R remain independent. Uniform W0
given H1 and the nu*s-bit left message M imply
rank_aux(W0|H2)>=h-nu*s>=KJ+q on every nonempty fiber,because M is affine in the auxiliary state tuple. A0 is right-sided.
The seed-marginal lift of strong affine J (error<=eps plus the2eps seed distance in(10),bounded conservatively by4eps) gives
D_aux(Z0|H2,R)<=4eps. (11)
Reveal only the nu new A_i for i in N. They are right functions after
M is fixed;let H3=(H2,A_N). All other A_i in Jset were already known
in H1. Since (11) retains R,appending A_N costs no extra TV.
On H3,Z0=J(W0,A0) is left-sided,and its seed distance is at most 4eps.
The right entropy loss in this step is at most nu*a. Prefix M was left-only,
so (8) gives Havg_aux(Y0|H3)>=KG+q. Product still holds. Symmetric robust
G gives
D_aux(B0|H3,whole state tuple)<=6eps. (12)
The retained information here determines Z0,and contains precisely all
A_i in Jset;M is already a function of the full state tuple.
Thus (12) is the auxiliary version of (9). Deterministic evaluation of
all relevant local A/Z/B values contracts the delta distance between
local laws. Fact (b) transfers the conditional-uniform statement,including
its changed target marginal,at cost 2delta. This proves (9).
In the REAL law the information (H1,whole state tuple,A_Jset) has its own
rectangular order: reveal states on the left,then selected A values on
the right. Z0 is determined there. This fact is separate from the
auxiliary-law argument.
### 8.2.5. Five-prefix bound
Code coordinate degrees are(3,3,2,2). The first honest B occurs only at1 or2,both degree3. At most two other labels omit it,so nu<=1+2*2=5. At the second honest B,unseparated A-choice labels contain the first B coordinate but omit the second;at most two qualify. Thus nu<=5 again. Missing copies only reduce these counts;h>=5s+KJ+q pays(8).
### 8.2.6. Projection to a real F history
Define Snew=S in an honest A round,and
Snew=S union {tampered indices currently choosing A} in an honest B round.
Let P be the REAL history obtained from old kept H by revealing:
1. all current selected summaries Q;
2. actual selected C_i only for i in Snew.
This is independently rectangular: Q is left-only,and each of these C_i
is right-only once its selected summary is fixed. It contains every next
separated output's seed,but does not append the honest selected seed yet.
If the honest choice is B,P is a deterministic projection of the retained
information in (9). Q is a function of the full W tuple. For old S,C_S
was fixed in H1;newly separated copies choose A and their C_i=A_i is
among A_Jset. Thus (9) implies
D(B0|P)<=2delta+6eps. (13)
Because honest Q0 is its full W0,B0 is right-sided on the REAL P.
If the honest choice is A,apply symmetric robust G at H1 BEFORE revealing
unseparated states. By (7) its honest seed Prefix_s(W0) has distance at
most delta. The Y entropy bound proved below gives
D(A0|H1,L)<=delta+2eps. (14)
All current summaries are functions of retained L and old H. Here Snew=S,
so P is a deterministic projection of (H1,L). Thus
D(A0|P)<=delta+2eps. The honest summary contains its prefix,so A0 is
right-sided on the real P.
In either case,each next separated output
W_i^+=F(X_i,C_i),i in Snew,is left-sided on P,since its C_i is fixed.
Reveal these at most t h-bit outputs on the left. Product survives and
fact (a) preserves the honest selected seed's distance. Call the result
Pplus. Once the real source entropy bound below is supplied,the affine F lift
thus yields,retaining all R and all next separated outputs,
honest A: delta_new <= delta+4eps;
honest B: delta_new <= 2delta+8eps. (15)
Now append all remaining actual C_i,including the honest selected seed.
They are functions of retained R and selected summaries in P. This adds
no TV cost after the seed-retaining F conclusion. Drop temporary
information and keep exactly H_r and W^+_(Snew),still retaining all R.
This proves invariant (6) for the next round.
The independence at F comes from the REAL message order defining Pplus.
### 8.2.7. Every entropy ledger and error recurrence
For each real F call,H contains past selected summaries,P adds only
current selected summaries,and Pplus adds at most t next outputs.
By (5) its total affine left rank loss from original X0 is at most
2c*h+2c*s+t*h = (2c+t)h+2cs. (16)
Its right messages do not change the conditional left law. Equation (1)
thus gives rank(X0|Pplus)>=KF+q on every nonempty real affine fiber. The affine F lift is the usual seed-marginal comparison on each such fiber: replacing the honest seed by uniform costs its conditional TV distance;strong affine extraction costs at most eps;remaining seeds are appended by their source-independent right kernel. No nonlinear-fiber claim or source-entropy tail is used. Conservative bounds in(15) remain valid. Initialization needed only KI.
For an honest B in round r,H1 costs at most
c*s0+c*(r-1)a+|S|a right bits. The local proof adds only nu new A values.
N and S are disjoint index sets among the c copies,so |S|+nu<=c. Hence
c*s0+[c*(r-1)+|S|+nu]a <= c*s0+c*r*a <= c*s0+4c*a. (17)
This proves (8) for every local G call. The first local G has no greater
right cost;prefixes cost no right entropy. The auxiliary law keeps the
same H1 marginal and right conditional law,so this entropy statement
is valid there as well.
For an honest A call the right cost is at most
c*s0+[c*(r-1)+|S|]a <= c*s0+[c*(r-1)+t]a < c*s0+4c*a.
Thus its G call in (14) is legal. No unselected A_i of an old separated
B-choice copy is revealed. A newly separated A-choice value is charged
once,and an already known old separated A-choice value is free.
The h guard is Section8.2.5. The robust tail at each honest extraction
has mass at most eps;no tail is charged per tampered copy.
Starting delta0<=eps,(15) has A:delta+4eps,B:2delta+8eps.
An adjacent A followed by B costs4eps more than B followed by A.
Thus the two-A/two-B worst ordering is AABB:
eps ->5eps ->9eps ->26eps ->60eps<=76eps.
The five actual code coefficients are36,40,44,48,52.
Every different weight-two label has an honest-1/other-0 coordinate.
Thus every tampered copy is in S_4. Invariant(6) retains all other
outputs and all c original seed strings. Dropping H_4 proves(2),with
the conservative76eps bound. Output truncation contracts TV.
At t=8,c=9,the right and left ledgers and local guard are36a;26h+18s;5s+KJ+q.
## 8.5 Affine fibers and generated auxiliary records
Let Z be uniform on an affine space and X a fixed affine image of rank k. Conditioning an r-bit affine map H of Z on a nonempty value leaves a uniform affine fiber,with rank(X|H)>=k-r. Indeed its direction is ker(H) inside the old direction V,and dim(T_X V/T_X ker(H))<=dim(V/ker(H))<=r. A linear image of a uniform affine fiber is uniform because all its nonempty point fibers have equal size. If an auxiliary is sampled by a kernel depending on Z only through fixed H,conditioning that auxiliary leaves the same Z fiber. Apply this argument inductively to adaptive affine messages whose coefficients are fixed by the old history. Arbitrary nonlinear leakage does not meet these hypotheses.

## 8.6 Selected two-sided-affine variant
Use wrappers4.2's actual four-round algorithm and five codes. States have h1 bits before their FIRST B,h2 before their SECOND B,s before A;early F_l requests the next public width and Fz outputs z. I_h1/I_s have seeds s0/s0s with s0s<=s0. All selected maps are strong affine and source-linear. Both complete input tuples are affine images of independent uniform affine latents,honest X has rank>=k,honest Y is uniform. Suppose az<=a,sz<=s,s<=h2<=h1,z<=h2,
 h1>=9s+KJ+q,h2>=max(5s+KJ+q,5sz+KJz+q),
 k<=n,k>=max(KI_h1,KI_s) and every weighted threshold below,
 d>=max(4s0+5s0s+40a+KG+q,4s0+5s0s+42a+11az+KGz+q).
Then error76eps retains the complete9seeds and8other outputs. Labels and their honest index are fixed in analysis;the evaluator receives no such index.
### Real histories and initialization
For each fixed honest code choose the following analysis history. In a round where its bit is A,keep the old Q_i:whole W_i for B copies,whole s-bit state for A copies. Reveal A_i for B copies before their C_i as in the affine order below. In a round where its bit is B,keep instead T_i for EVERY copy and Z_i for B copies. The real order is T_all left,then A_all right,then Z_B left,then selected C_all right. Given T,each A_i is affine right;given A_i,each Z_i is affine left;given Z_B,each B-choice C_i is affine right. A-choice C_i=A_i is already fixed. At every step coefficients are fixed by earlier history. Hence both latents stay independent uniform affine on every nonempty fiber,and actual bits bound rank loss. Current W_i of either width are affine left with old selected C fixed. Old history and full R determine every old right revelation.
The initial prefix record costs <=4s0+5s0s. Two codes start B,three start A:with two copies of each other label,honest B gives3s0+6s0s;honest A gives4s0+5s0s. Missing copies decrease costs. Each selected I map uses its own uniform prefix and paid KI_h1/KI_s. Strong retention extends the remaining R kernel. Thus the initial separated-state invariant is delta0<=eps,as in8.2.3.
### Local B comparison before projection
Let S be old separated copies. Form H1 by temporarily revealing their full W_S,then A_S_B and selected C_S. This is affine in rectangular order,and costs |S|+|S_B| current right seed blocks. Since H_r and full R were retained,8.2.3's invariant gives D(W0|H1,R)<=delta. Use8.2.4's affine section extension on the WHOLE current-state tuple,keeping W_S fixed and the independent right law unchanged. Honest W0 becomes full uniform honest state width with joint replacement cost<=delta;other states may have different widths. Old C_S equations survive.
Now let N={0} union ALL tampered copies outside S,nu=|N|<=9. At the SECOND honest B,nu<=5:every other B-choice label was already separated at the honest code's first B coordinate,otherwise its two B coordinates would equal the honest code. The unseparated A labels contain the first B coordinate,whose degree is at most3,so at most two such other labels remain. Put M=(T_i:i in N),of nu*s bits (nu*sz short). Old S prefixes are already known. Before conditioning on M,strong affine G at the uniform T0 gives the same bound(10),retaining the whole state tuple. On each M fiber uniform honest W0 loses at most nu*s rank. The h1/h2 guard pays strong J on the actual width;the seed-marginal lift proves(11),error<=4eps retaining full R. Reveal A_N on the right. Since R was retained,this costs no further TV. Their right loss is at most nu*a. The displayed d guard pays KG+q after these affine right revelations. The second strong G,with Z0 within4eps of uniform,gives(12),error<=6eps retaining the whole state tuple and all A values. Other Z_B are determined by that tuple and A_all. Transfer to the REAL law costs2delta,including the changed target marginal. Thus real C0 has distance<=2delta+6eps from uniform even retaining H1,M,A_N,whole state tuple and all Z_B.
### Projection and honest F
Project that real conclusion to P consisting of old history,T_all,A_all,Z_B and selected C_S. All these are determined by the retained local record. The conditional kernel of the full raw left latent given the whole current-state tuple and that local record is independent of the right latent:the record had an explicit rectangular affine order. Append that raw latent by this kernel before projecting. Consequently C0 is within2delta+6eps of uniform even given P and raw X. This does NOT claim that projection alone implies rectangularity. The separate REAL construction of P is T_all left,A_all right,Z_B left,C_S right,so both independent affine latent fibers remain. Temporary full W_S is not kept in P and incurs no real source-rank charge.
For Snew=S union current A-choice copies,the selected C_Snew is fixed in P. Reveal their next F outputs on the left. They are affine functions of raw X at fixed seeds,of the actual next widths. They may be appended to the preceding C0 comparison since raw X was retained. The weighted k bounds pay the honest source threshold after this real left record. The strong affine F seed lift then gives error<=2delta+8eps retaining FULL R and the separated next outputs. Only afterward append the other selected C_i;their Z_B or A_i is already fixed,so they are affine right. This establishes the next compressed real history and separated-state invariant.
When the honest bit is A,use8.2.6's A proof with the corresponding G/Gz and F/Fz roles with whole-state Q_i:the honest near-uniform prefix is used BEFORE further state revelation. Project its comparison to the actual Q,C_S,new separated outputs,then apply F retaining full R. Append remaining B-choice A_i,C_i only AFTER this F conclusion. Its error is delta+4eps. The previous compressed histories are independent affine fibers,which is all this step needs;they need not have exposed whole earlier states. No unneeded tampered A_i is conditioned on before an honest A extraction.
### Finite counts and closure
Only real kept messages charge ranks. For honest A,each copy contributes its actual current h1/h2/s width;for honest B,each B copy contributes2s,each A copy s. Add actual next widths of Snew. With two copies of each other label,the table gives(A,B,C,t) for rounds1--4,t=1/2/3 meaning next F output h1/h2/s,and t=4 finalFz:
1100:(4,0,14,2),(2,4,28,3),(2,8,33,3),(2,8,34,4).
1010:(4,0,14,3),(6,4,17,2),(4,6,33,3),(4,6,34,4).
0110:(4,0,5,1),(6,2,19,2),(4,4,37,3),(4,4,38,4).
0101:(4,0,5,1),(6,2,19,3),(6,6,24,2),(6,4,34,4).
0011:(4,0,5,3),(8,2,8,1),(8,4,24,2),(8,2,34,4).
For t<=3 require k>=A*h1+B*h2+C*s+KF_(l_t)+q. For t=4 require k>=A*h1+B*h2+C*s+8z+KFz+q. Positive counts are nondecreasing in other-label multiplicities,so missing copies reduce them. These20rows pay every actual honest F invocation;J_l uses the actual current B ambient with the same gamma and output s. The short Jz ambient is h2 because a final B is the second B. No prefix-of-extractor-output assertion is used.
For right costs,before round r with c<=9 copies put bprev=number of completed B choices over all copies. Prior cost is c*(r-1)+bprev early seed blocks. In an honest A round current local cost is |S|+|S_B|. In an honest B round it is |S|+|S_B|+nu=c+|S_B|. Direct five-code substitution gives early totals by rounds1,2,3: (9,25,40),(9,22,40),(0,22,40),(0,22,35),(0,13,37). Thus40a suffices. Final(prior a,current az) pairs are(41,12),(41,12),(41,12),(42,11),(42,11),all<=42a+11az since az<=a. The initial prefix bound is paid separately. These rank guards are used before extraction.
The selected codes are distinct and each has two B coordinates. For every other code,one honest B coordinate is A there;this directed antichain property gives separation by round4. The second-B prefix bound was proved above,without assuming the first B occurs before round3. The recurrence is A:delta+4eps,B:2delta+8eps,starting delta0<=eps;each two-B code costs<=76eps. Every other label is separated by round4. Dropping surplus real history retains all required original seeds and outputs,proving the theorem.
