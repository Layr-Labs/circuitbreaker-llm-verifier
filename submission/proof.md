# A quantitative eventual bound for the norm-one language

For a vertex and three incident edge labels, count the labels that disagree
with the anchor at that vertex. For respectively zero, one, two and three
disagreements, charge \(0,1,1+\delta,\lambda\), with
\(\lambda=47091/100000\) and \(\delta=2102599/1000000000\).
Let \(U=13064756648/10^{11}\) bound the expected crossing-cover density
and \(X=8113009/50000000\) bound the crossing-edge density.
The replay proves lower bounds \(s_1=10154921783199/3276800000000000\) and
\(s_2=380936142031507263/287729254400000000000\) for disjoint populations of
saved vertices of crossing degree one and two on a single cover.
A further local adjustment at initially isolated anchors gives a
charge rebate \(J=406564145190519541121747715463/8092385280000000000000000000000000\), proved in [family.md](family.md).
The penalized anchor objective and the unit-factor interpolator below give
a uniform simple cubic cutwidth constant
\[
 a_0=(1-\delta)U+\delta X-(1-\lambda-2\delta)s_1-(1-\lambda)s_2-J
   =\frac{148364259508455484820287437447791}{1156055040000000000000000000000000}
.
\]
The quantitative graph estimates in [family.md](family.md) sharpen the
uniform width loss to \(O(n/(\log_2 n)^2)\) and transfer it at the
exact coefficient \(1+1/(2a_0)\).

The four component lower bounds, in positive and negative order within each crossing degree, are
\[
\frac{6132802254699}{3276800000000000},\frac{8044239057}{6553600000000},\frac{647713169760441}{838860800000000000},\frac{39692631200919}{71932313600000000}.
\]

This submission refines the negative quadrature of hrushikesh’s [provisional 4.895796 proof](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/a49ccfaf5eb17002ea7edfb90c9fe76bcc9e2c82/submission). The graph argument, norm transfer and event definitions below follow that work. The positive replay retains the published mesh.

The two-phase graph conversion, norm-one transfer and radial certificate
are inherited through xlib's [4.89223 submission](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/3e58bdc08d3f5f0e216606c77776f1cc6a0815f4/submission).
The expanded endpoint, far-arm and mixed negative-neighbor guards come
from the [4.89253 submission](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/51bd01f184d39d2105a2764df51bf57792b300e9/submission).
The curvature correction, expanded root events and integer integration
are adapted from xlib's [PR 895](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/cec7eb07dec8daa984c7138f119365364449c414/submission).
The directed integration uses the full endpoint bank for the compressed
branch and a point subset for individual histories. These histories extend the covered
negative same-sign neighbor on the same crossing cover. We directly
recompute both root pairs. The sign-specific grids come from
[PR 902](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/6cd43769f675208d27ba8ec380a7eb5180843121/submission).
The guarded positive rays and mixed-arm center calculation come from
[PR 905](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/40f0a2292a512fb3cf23df9010c66741acc1bb19/submission). The incumbent enlarged the ray aperture. The exact
positive and negative runs together verify all eleven floors.

## Cover descent and a physical anchor order

Let \(G(b)\) be the minimum total charge among edge bits at anchor word
\(b\). The interpolation proof below verifies joint submodularity in
edge and anchor bits. Complementing all signs preserves the objective.
Changing one anchor changes \(G(b)\) by at most \(\lambda\), since the
paired star charges differ by either \(\lambda\) or \(\delta\).

For a degree-\(d\) finite cover the minimum for repeated anchors is
\(d\) times the base minimum. Lifting a base solution gives one
inequality. To descend, average the edge bits of an upstairs solution
in each base edge fiber. Threshold all averages at one uniform point
of \([0,1]\). For a submodular function \(F\) on finitely many bits,
the expected value under such a common threshold of coordinates
\(y_i\in[0,1]\) is convex in \(y\). Indeed for every permutation
\(\sigma\), use the corresponding nested prefixes \(S_i\) and affine
function \(F(\varnothing)+\sum_i y_{\sigma(i)}[F(S_i)-F(S_{i-1})]\).
At a vertex of the Boolean cube this affine function is at most \(F\),
by telescoping over only the present coordinates and using diminishing
returns. A permutation sorting \(y\) in descending order attains
equality in expected value, since every threshold set is one of its
prefixes almost surely. The extension is thus the maximum of these
affine functions. Apply convexity to the four coordinates of each
star, matching the incident edge fibers by the covering maps. Its
expected threshold charge is at most its average upstairs charge.
The anchors are already Boolean, so the induced base optimum is at
most the normalized upstairs optimum.

Finite covers of arbitrarily large girth exist for each connected
graph. Contract a spanning tree to identify its fundamental group
with a free group on the other edges. A short nonbacktracking closed
walk gives a nonidentity reduced word. Separate each such word in a
finite permutation group by making its successive positions distinct
points. Realize its letters by partial permutations taking each
position to the next, using the inverse step for an inverse letter.
Reducedness prevents contradictions of domains or images. Extend
each partial permutation arbitrarily. The word moves its first
point. Use the product over all words up to the required length and
its regular action as voltages. The resulting cover has no closed
lift of one of those walks.

On the cubic tree take unit rows proportional to
\(\alpha^{d(v,z)}\), \(0<\alpha<1/\sqrt2\). There are
\(3\cdot2^{i-1}\) vertices at distance \(i>0\). For center distance
\(j>0\), partition column vertices by projection onto the center
path. The two endpoints contribute transverse generating series
\(1/(1-2\alpha^2)\) and the \(j-1\) internal vertices contribute
\((1-\alpha^2)/(1-2\alpha^2)\). The squared norm before normalization
is \((1+\alpha^2)/(1-2\alpha^2)\). The resulting correlations are
\(\alpha^j[(j+1)-(j-1)\alpha^2]/(1+\alpha^2)\), tending to
\(r_j=(1+j/3)2^{-j/2}\).

Truncate at radius \(s\) and renormalize. At fixed \(\alpha\) the
geometric squared tails make the row error uniformly small. In a
high-girth cover every local truncated row is a tree row. Take
Gaussian scalar products against independent normals at the cover
columns. For a physical row use the norms of the cover entries in
each base column fiber, with independent normals at the base
columns. Cauchy–Schwarz shows that the repeated physical keys have
entrywise larger off-diagonal covariance than the upper keys. Add
independent Gaussian tie noise separately in the physical graph
and cover, with the same subsequent normalization. The lifted
physical noise still preserves domination.

Increasing off-diagonal covariances decreases the expected value of
a submodular bit function at a common threshold. Smooth the
thresholds in the multilinear extension of that function.
Differentiating along the covariance interpolation and integrating
the Gaussian density by parts multiplies each covariance change by
the expected mixed second derivative. This derivative has the sign
of the nonpositive mixed discrete difference, since the smoothing
derivatives have the same sign. Approximate singular matrices by
positive definite matrices, and take the sharp threshold, using
continuity of the one-dimensional marginals.

For these centered Gaussian laws the expected cost is largest at
threshold zero. Its derivative at a common threshold \(t\) sums
the density at each \(X_i=t\) times the corresponding mean first
discrete difference. Conditional on \(X_i=t\), the other tests are
\(X_j-r_{ij}X_i<(1-r_{ij})t\), with residual law independent of
\(t\). The conditional first difference is decreasing, by
submodularity, and odd in \(t\), by complement symmetry. This
gives a nonpositive derivative at \(t>0\), and the singular case
follows by approximation.

Choose \(\alpha\), truncation radius and noise before the size of
the graph, so that all finite neighborhoods used by the probability
bounds have the limiting law within the desired tolerance. The
cover identity then implies for every physical threshold an
expected cost at most \((a_0+\epsilon)n\).
At fixed truncation radius each independent coordinate changes
only a bounded number of key bits. The resampling variance of
\(G(b_t)\) and of the key count in a fixed interval is \(O(n)\).
Take finitely many standard-normal quantile bins of mass at most
\(\epsilon\), including the tails. Chebyshev and a union bound
give one realization with all boundary costs within \(\epsilon n\)
of expectation and at most \(2\epsilon n\) keys in each bin.
The Lipschitz bound extends this to every prefix of the order by
those distinct keys. Apply the two-phase interpolation. Letting
the fixed tolerance tend to zero gives the uniform width
coefficient \(a_0\).

## One cover and its saving events

The limiting wave satisfies \(\sum_{u\sim v}X_u=2\sqrt2X_v\):
the difference has variance zero by \(r_1,r_2\). Consequently
three opposite neighbors are impossible. In each high-girth
approximation remove vertices of crossing degree three and the
bounded-radius witness centers near them. The expected discarded
fraction tends to zero. Each removed vertex costs at most a
bounded number of star charges and bounded-length walks.

The remaining crossing graph consists of paths and even cycles.
Cover an even-length path by its smaller parity class, and an
odd-length path or a cycle by its positive class. Ignore isolated
vertices. The size is minimum. In an odd path list positive-negative
pairs from the positive endpoint. Switch the first consecutive
pairs whose positive members pass the following test. A side
positive is an outside guard if it has no opposite child or
all its crossing arms close at even length. It is an unsaved
guard when it has one closing odd arm and a same-positive sibling
direction with an outside guard. Switch the initial endpoint if
both of its sides are outside, or one is outside and the other
unsaved. Switch an internal positive if its one side is outside.
All switched positives and their same-positive neighbors were
ineligible for an old saving. The outside guards stay out of
the cover since their components were even paths or isolated.
Thus old saved positives and all their needed neighbors survive.
Switching keeps the cover size. A proper nonempty switched prefix
has one double-covered edge at its end.

Give singly covered crossing edges the uncovered sign and a
double-covered boundary edge the negative sign. Each vertex in the chosen
cover \(K\) has either one or two crossing disagreements. In fact at a
double-covered edge the negative endpoint loses this disagreement. It has
another crossing edge to the positive vertex just switched out of the
cover. That edge gives it a disagreement. Each crossing edge has exactly
one disagreeing endpoint, so there are precisely \(|E_{\rm cross}|-|K|\)
vertices with two crossing disagreements. The remaining edges start
with their common endpoint sign.

Choose the witnessed sets \(S_1,S_2\subseteq K\) of crossing degrees one
and two, whose crossing edges all disagree with their own signs and whose
same-sign neighbors all belong to \(K\). The proofs below certify the
members of these two sets. Flip the noncrossing edges incident to them.
Every flipped edge has both endpoints in \(K\). All witnessed vertices
then become unanimous disagreements and cost \(\lambda\). Before the
flips, the vertices of \(S_2\) were among the two-disagreement stars.
After the flips, an additional two-disagreement star outside the
witnessed sets needs at least one flipped edge. There are at most
\(2|S_1|+|S_2|\) such edges. Thus the final number of two-disagreement
stars is at most \(|E_{\rm cross}|-|K|+2|S_1|\). Each other member of
\(K\) costs at most the split base charge one, apart from the extra
penalty just counted. Stars outside \(K\) still match their anchors.
The total charge is at most
\[
 (1-\delta)|K|+\delta|E_{\rm cross}|
 -(1-\lambda-2\delta)|S_1|-(1-\lambda)|S_2|.
\]
This bound includes the negligible deletion debit in the high-girth
approximation, where the finite pattern probabilities converge. Both
coefficients of the savings are positive. The bound for the crossing
cover uses the reference enclosure below. The limiting density of
crossing edges is \(3\arccos(2\sqrt2/3)/(2\pi)\), since a rotationally
invariant Gaussian pair represented by two unit directions of angle
\(\theta\) has opposite signs on angular measure \(2\theta\) out of
\(2\pi\). Here \(\arccos(2\sqrt2/3)=\arctan(1/(2\sqrt2))\).
For an outward rational bound, the replay takes
\(u=353553391/10^9>1/(2\sqrt2)\), upper-bounds \(\arctan u\) by the
nine-term alternating series through \(u^{17}/17\), and divides by
the lower endpoint of the certified \(\pi\) interval. This gives the
bound \(X\) above. Apply the additional isolated-center adjustment in [family.md](family.md)
after this labeling. It reduces the expectation of the displayed upper
bound by at least \(J\), yielding \(a_0\) on the limiting local model.

## Conditional branches and root densities

Write \(r=\sqrt2\), \(g(h)=2\sqrt{6/\pi}\exp(-6h^2)\) on
\(h\ge0\), and \(H(u)=\int_0^{\min(2,u)}g(h)\,dh\) for \(u>0\),
zero otherwise. At a vertex with parent value \(a\) and its own
value \(b\), the children have values \(rb-a/2\pm\eta\),
where \(\eta\) has variance \(1/12\). The sum follows from the
eigenvalue identity. The half difference has zero covariance with
every score outside its two subtrees by equality of distances.
Distinct forward contrasts are also uncorrelated: in the nested
case the deeper contrast is orthogonal to both coordinates of
the earlier contrast. Gaussianity makes all these innovations
independent. The root star therefore gives independent conditional
branch tests in the three forward subtrees.

For an opposite parent/current pair of magnitudes \(a,b\), put
\(B_0(a,b)=H(rb+a/2)\) and \(B_{j+1}=K B_j\), where
\[
 Kf(a,b)=\int_0^2 g(d+rb+a/2) f(b,d)\,dd .
\]
These certify \(j\) extra crossing edges beyond the current vertex,
followed by closure. The next opposite child is unique because the
midpoint in the current sign convention is positive. Use the two
closing parities \(I=(0,2,4,6,8)\), \(O=(1,3,5,7)\). A subscript
\(I\) or \(O\) in a root formula denotes the sum of those components.
For a same-sign entry \(a,b\), let \(m=rb-a/2\). The operator
\[
 Sf(a,b)=\int_{\max(0,-m)}^2 g(d+m)f(b,d)\,dd
\]
chooses as opposite child the more negative forward value in
positive-current coordinates. Starting instead at \(\max(0,-2m)\)
gives \(S_1\), where the sibling is nonnegative. When both are
opposite use
\[
 D(f,k)(a,b)=\int_0^{\min(2,-m)}
   g(h)f(b,-m+h)k(b,-m-h)\,dh
\]
for \(m<0\), zero otherwise. The family \(A_0=\{SB_j:j\in I\}\) certifies that a positive vertex
belongs to the original cover. The more negative child has odd
terminal distance from the vertex. In an even path the vertex then
belongs to the smaller parity class, and an odd path uses its
positives. Include all twenty \(D(B_i,B_j)\), \(i\in O,j\in I\), to
form \(A\). These give the more negative child the complementary
parity and require the other child to have odd terminal distance.
They are disjoint from \(A_0\) and from one another.

### Which odd paths switch

For the outside guards use the disjoint family
\[
\mathcal W=\{U=H(m)\}\ \cup\ \{W_j=S_1 B_j:j=1,3,5,7\}
   \ \cup\ \{W_{ij}=D(B_i,B_j):i,j=1,3\}.
\]
An outside vertex has zero, one or two opposite children. Every
terminal arm in these guards has even length. This puts the vertex
outside the minimum cover of an even path or leaves it isolated.
The notation \(W_{ij}\) denotes a two-child guard.
For a vertex with just one opposite child of magnitude \(d\) and
the other child in outside class \(o\), define
\[
V_{i,o}(p,u)=\int_{\max(0,-2m)}^2
       g(d+m) B_i(u,d)\,o(u,2m+d)\,dd,\qquad m=ru-p/2 .
\]
The family \(\mathcal V\) takes \(i=0,2,4,6\) for \(o=U\),
\(i=0,2,4\) for \(o=W_1\), \(i=0,2\) for \(o=W_3\), and \(i=0\)
for \(o=W_{11}\). Its unique crossing arm has odd terminal length
and its same-positive child is outside the cover. Consequently the
vertex cannot save before switching. The ten classes and nine
outside guards separate by number, lengths, ordered sides or the
additional outside guard.

At a positive endpoint entered from an opposite parent of magnitudes
\((a,b)\), put \(M=rb+a/2\). Each ordered pair \((f,k)\) in
\[
\mathcal T=(\mathcal W\times\mathcal W)\ \sqcup\
(\mathcal V\times\mathcal W)\ \sqcup\ (\mathcal W\times\mathcal V)
\]
certifies its admission. Use the integral
\[
T_{f,k}(a,b)=\int_0^{\min(13/10,M)}
        g(h) f(b,M-h) k(b,M+h)\,dh .
\]
The order is by child magnitude and ties have zero measure. For an
internal positive, multiply the \(K\) integrand by \(o(b,2rb+a+d)\)
and denote the operator by \(K'_o\). Use
\(\mathcal W_{\rm int}=\{U,W_1,W_3,W_5,W_{11}\}\).
The near-arm family \(E\) consists of the 261 endpoint tests and
all one through five guarded two-step continuations of those tests,
with a nonzero continuation given by \(K'_oK\) applied to the shorter
history. Thus its terminal distance, measured from the negative
neighbor that enters the arm, is odd and at most eleven. Every
positive along the path from the positive endpoint is admitted.
Lengths, endpoint guards and successive internal guards distinguish
the histories.

For a far arm of even terminal distance define
\[
 O_g=\{K'_U B_j:j\in I\}\ \cup\
 \{K'_{W_1}B_j:j=0,2,4,6\}\ \cup\
 \{K'_{W_3}B_j:j=0,2\}\ \cup\ \{K'_{W_{11}}B_0\}.
\]
These twelve cases guard the next positive across the negative root.

For point evaluation use a subfamily \(E_{\rm pt}\) of 74 histories
in \(E\). Start with the first 40 ordered endpoint pairs of the
predecessor. Within the original eleven guard names, enumerate the
ordered pairs in lexicographic index order, keeping those with at
least one outside guard and absent from the first 40. Append the
first four of these. The resulting point bank has 44 endpoints.
The numerical IDs \(0,\ldots,13\) used in that list
denote, in order,
\[
 U,W_1,W_3,W_5,W_{11},V_{0,U},V_{2,U},V_{4,U},V_{0,W_1},
 V_{2,W_1},V_{0,W_3},W_{13},W_{31},V_{0,W_{11}} .
\]
Every pair has an outside member. The remaining histories have one
guarded hop (16 using \(U\), three \(W_1\), two \(W_3\)), two hops
(six with outer guard \(U\), one with \(W_1\)) or three hops (two
with outer guard \(U\)). The listed indices choose among existing
shorter histories. Each differing sequence certifies a disjoint
event even if a downward numerical floor happens to coincide.

### Covered negative neighbors and the four root integrals

Write \(P_S=\{S_1e:e\in E_{\rm pt}\}\). Its same-negative vertex
has only one opposite child. The odd path from the positive
endpoint switches through this negative endpoint, which belongs
to the resulting cover. Also use the following two-child cases
at a same-negative vertex.

* In \(D(B_i,B_j)\), \(i,j\in I\), both crossing arms have odd
  length. This is a negative vertex in the smaller class of an
  even path. It stays covered.
* In either ordered orientation of \(D(e,B_k)\), \(k=1,3,5,7\),
  the near \(e\) arm finishes at a positive endpoint and has odd
  length. The opposite far arm has even length. All positives on
  the near arm pass the admission tests. The last of their switched
  pairs includes the negative vertex under consideration. It is
  covered even when it is the negative endpoint of the double-covered
  boundary edge. No saving at this particular neighbor is needed.

For these mixed cases the near-arm list starts with the first 12
point endpoints, the first four one-hop \(U\) histories and the first
two two-hop \(U\) histories, in this order, indexed from zero.
Use all 18 entries for \(k=1\), entries \(0,1,2,3,4,5,12,13\)
for \(k=3\), the first seven for \(k=5\) and the first four for
\(k=7\). This gives 74 ordered mixed cases alongside the 25
odd/odd cases. Their union is \(P_D\). The number of children
separates \(P_D\) from \(P_S\), then ordered lengths, parities and
near-arm guards distinguish its cases. Put \(P=P_S\sqcup P_D\).
For any disjoint family in a probability formula use its sum.

For crossing degree one, fix the sign, root magnitude \(x>0\)
and the opposite magnitude \(b>0\). Their joint density is
\(f_2(x,-b)=(3/2\pi)\exp[-9(x^2+b^2)/2-6rxb]\).
Set \(s=b+2rx,\ h=ts/2,\ y=s/2+h,\ z=s/2-h\), with \(0\le t\le1\).
The other two magnitudes are \(y,z\) and the remaining density and
Jacobian are \(g(h)s/2\). At \(y\) the forward midpoint is
\(ry-x/2\ge 3x/2\), so both-opposite cases are impossible there.
For positive and negative roots respectively multiply by
\[
 A_0(x,y) A(x,z) B_I(x,b),\qquad
 P_S(x,y) P(x,z) E(x,b).
\]
In the first case the root had an old saving. In the second case
the unique opposite arm switches completely through the negative
endpoint, with both its same-sign neighbors now covered. Integrate
on the sign-specific box in the finite-arithmetic section and multiply by three
for the placement of the unique opposite neighbor.

For crossing degree two denote the root magnitude by \(x\), the
opposite magnitudes by \(y,z\) and the same-sign magnitude by
\(b=2rx+y+z\). From covariance \(2r/3\), the first neighbor given
\(x\) has variance \(1/9\); the second neighbor uses the independent
forward contrast. Thus the joint density for one root sign is
\[
 f_3(x,-y,-z)=\frac{\sqrt{108}}{(2\pi)^{3/2}}
 \exp[-x^2/2-9(y+2rx/3)^2/2-3(2z+y+2rx)^2/2].
\]
The forward midpoint at \(b\) is nonnegative. With
\(u=B_I(x,y)\), \(v=B_I(x,z)\), the positive and negative
multipliers are
\[
 A_0(x,b)(u+v-uv),\qquad
 P_S(x,b)\,[uv+E(x,y)O_g(x,z)+O_g(x,y)E(x,z)].
\]
For the positive root either odd-distance terminal arm puts it
in the old cover, giving an old saving with the covered same-sign
neighbor. Independence of its two arm tests gives the conditional
union \(u+v-uv\). For the negative root, two odd distances put it
in the smaller class of an even path. A mixed near and far pair
switches through the negative and the guarded next positive.
This additional guard places the negative strictly before any
double-covered edge, so both crossing labels disagree with its
anchor. The three negative cases are disjoint by ordered parity.
Integrate on the sign-specific box below, with factor three
for placement of the same-sign neighbor. The four root functions
count disjoint classes by root sign and crossing degree in the
same final cover. All their tests inspect a bounded neighborhood.

## Directed lower integration

Every individual itinerary is jointly log-concave in its entering
magnitudes. It integrates products of Gaussian densities on domains
defined by affine constraints on magnitudes and the ordered absolute
innovations. This includes all path truncations and the aperture
for \(T\). Marginalization preserves log-concavity: transport the
normalized one-dimensional sections by their increasing quantile
map \(L\), use joint log-concavity at the interpolated points, and
bound the Jacobian \((1-\theta)+\theta L'\ge(L')^\theta\). This yields
the geometric interpolation of the two original masses as a lower
bound. Iteration and approximation cover additional coordinates
and sections of mass zero.

In a one-dimensional panel a log-concave integrand lies above
the exponential chord of its two positive endpoint lower bounds.
Its average is at least their geometric mean. The replay for
individual point tables uses the mean of the two chord values
at one-quarter and three-quarters, by Jensen on the two half
intervals. A zero corner may give a zero lower contribution.
Apply the geometric-mean bound successively in all coordinates
of a root cell. The first root formulas are separately log-concave
for each itinerary in \(x,b,t\). All arguments and the nonnegative
Jacobian are affine in each individual variable with the others
fixed. In the second root formula the arguments are jointly affine.
Compute the eight-corner mean *per itinerary* before summing.
Independent factors in a single mean multiply.

For a near-arm endpoint with parameter corners \((a_i,b_i)\),
one may choose a different cap \(C_i\le\min(13/10,M_i)\) at each
corner and substitute \(h=C_i t\). For a left \(U\) guard its
argument is \(3b_i/2+a_i/r-rh\). In this case use the smaller cap
inside its positive support; other left guards keep the full
inward-rounded aperture. Four-term Hölder gives
\[
 \prod_i\Big(\sum_j q_{ij}\Big)^{1/4}
       \ \ge\ \sum_j\prod_i q_{ij}^{1/4}\qquad(q_{ij}\ge0).
\]
Apply it to the common panels at the four corners of each endpoint
integral, including each corner's own Jacobian. Within each panel,
multiply the four-corner means of the two guards and the density factor.
Adding the admissible ordered guard products constructs a floor
for the sum of individual four-corner means for the 261 endpoints.
For a hop \(K\) or \(K'_o\) use the same inequality on the common
integration grid for \(d\). Each corner value of the inner function
at \((b,d)\) occurs twice among the eight kernel corners. Its mean
is the four-corner branch mean on that input cell. Thus a positive
kernel propagates floors for sums of itinerary means, first by
\(K\), then by the sum of the guarded \(K'_o\). Repeat five times
and add the contributions of the six different lengths. The
point subfamily is integrated individually, allowing its members
to be used in another \(S_1\) or \(D\) test.


For the positive two-arm union we also need bounds throughout each cell.
At a given \((x,y)\) cell define \(Y\) as the sum, over its five
admitted odd-distance itineraries, of the log-bilinear interpolants
of their four certified corner probabilities. If a corner is zero
use zero in the cell interior for that term. Separate log-concavity
gives \(0\le Y\le u\). Let \(l_y\) be the sum of the five individual
corner minima and \(h_y\) the largest of the four corner sums.
Each interpolant is bounded by the weighted arithmetic average at
its corners, with the same weights for all five. Thus
\(l_y\le Y\le h_y\le1\). Define \(Z,l_z,h_z\) on \((x,z)\) in the
same way. Independence and monotonicity give the lower integrand
\(Y+Z-YZ\). Since \((h_y-Y)(Z-l_z)\ge0\), it dominates
\[
 (1-l_z)Y+(1-h_y)Z+h_y l_z .
\]
Exchanging \(y,z\) gives a second bound. Each coefficient is nonnegative.
For each term multiply by the common positive density and the
same-neighbor itinerary, then apply the corner-mean inequality.
The repeated \(Y\) values on the two \(z\) faces reduce to the four-corner
mean on \((x,y)\), and symmetrically for \(Z\). Two further valid,
disjoint descriptions within each orientation count success of the
first arm or certified even distance in the first arm and success of
the other. The replay takes the largest of the four lower bounds for
each cell. Their maximum still bounds the integral of the same union
there. Only disjoint spatial cells are added.

### Recovering the Gaussian curvature between corners

For a single nonnegative itinerary write its integrand as a known
Gaussian \(g_*\) times separately log-concave factors. In a positive
cell, Jensen bounds its mean below by the exponential of its mean
logarithm. The mean log of each remaining factor is at least the
uniform average of its corner logs, by successive concavity in each
coordinate. For \(g_*\) retain the extra amount
\[
 \delta=\mathbb E_{\rm cell}\log g_*-\mathbb E_{\rm corners}\log g_* .
\]
Thus the corner geometric-mean floor can be multiplied by
\(\exp(\delta)\ge1+\delta\). A term whose corner floor is zero
may be discarded. Applied separately to the nonnegative terms
in either positive-arm affine bound, the factor is common to all
terms, because all have the same root density. The Jacobian in
degree one is a further log-concave factor.

For a transition panel of width \(d\), the Gaussian log density
is a constant minus \(6(v+s)^2\) or \(6(s-v)^2\), where \(s\)
is fixed. At relative coordinate \(t\), its excess above the chord
of its endpoint logs is \(6d^2 t(1-t)\). The average of \(t(1-t)\)
is \(1/6\) on both the full interval and each of its two halves.
The ordinary geometric-mean panel and the two quarter-point
panel rule therefore both permit a multiplier \(1+d^2\).
When the two quarter points are used for \(K,S,S_1\) and their
guarded versions, compute the quarter-point chord values of the
kernel and the input branch separately, then multiply corresponding
values. This is the chord of the product, with all square roots
rounded down. On a partial first panel the branch first uses its
directed interpolation at the inward-rounded integration limit.
Point endpoint and \(D\) integrals use the same curvature correction
on their individual integrands. Each compressed hop has one Gaussian
factor. Its plain kernel supplies this correction once, including
when multiplied by a side guard. For a compressed endpoint the
temporal width is \(C_i/18\) at parameter corner \(i\); use its own
factor at that corner before applying Hölder. Successive independent
innovations in a nested itinerary each contribute a separate
Gaussian density.

For the degree-two root write \(d_x,d_y,d_z\) for its cell widths.
The mixed monomials in its quadratic log density have the same
uniform cell and corner averages. The average of \(x^2\) at the
ends of an interval exceeds its uniform mean by \(d_x^2/6\).
The root correction is therefore
\[
 \delta_2=(33d_x^2+12d_y^2+12d_z^2)/12.
\]
For degree one let \(d_x,d_b\) be the widths and take \(n_t\) panels in \(t\).
On \(t\in[k/n_t,(k+1)/n_t]\), let \(\bar s\) be the cell midpoint
of \(s=b+2rx\), and write \(V=8d_x^2+d_b^2\). The pair density
contributes \(9(d_x^2+d_b^2)/12\). The remaining Gaussian has
logarithm \({\rm const}-3t^2s^2/2\). The corner mean of \(s^2\)
is \(\bar s^2+V/4\), and its cell mean is \(\bar s^2+V/12\).
The cell mean of \(t^2\) is \((k^2+k+1/3)/n_t^2\), and its corner
mean exceeds this by \(1/(6n_t^2)\). Since \(t\) is independent
of \(x,b\) under both averaging rules, multiplying the two means
and taking their difference gives the total
\[
 \delta_1=
 \frac{36 n_t^2(d_x^2+d_b^2)+12\bar s^2+
             [12k(k+1)+7]V}{48 n_t^2}.
\]
The squared midpoint uses a nonnegative rational lower bound on
\(2\bar s\). The replay scales each surviving base floor by
\(1+\delta\), with a downward integer division. The Gaussian
factors in the root star are distinct from the forward innovations.


### Finite arithmetic

The finite replay restricts companion-variable panels to
\([0,3/2]\) on the shared grid, using zero for values outside its
tabulated domain. These restrictions bound nonnegative integrals from
below. The separate half-normal scalar \(H\) keeps its aperture two.
The exponent cutoff is \(21/2\), above which the density floor is
zero. Dropping low-weight root cells is likewise a lower bound.
Point endpoint quarter-panel values factor the chord of the density
and first guard from that of the second guard. Multiplying their
two directed quarter values gives the corresponding lower product
chord, by additivity of logarithms before rounding. Kernel indices
and fixed-point coefficients are stored in integer arrays. Computing
one guarded kernel row before accumulating its branch integrals
imposes the same outward divisions as tabulating all rows first.

The declared initial `{"part":"negative"}` checks both negative base
roots and all five center floors. Request `{"part":"positive"}` to
check the radial enclosure, both positive base roots and both ray floors. The two outputs together validate all
numerical premises in the width formula. Each checks the interpolation,
crossing-edge ceiling and the same exact rational transfer.

Positive mesh: 103 knots at scale 320, first positive knot 1/2560; degree-one upper coordinates (9/16, 63/80, 1); degree-two upper coordinates (41/80, 87/80, 87/80). Use 40 temporal panels for the degree-one root. The positive replay uses integer transition vectors, releases its
plain table after building the branches, and avoids the guarded
negative tables. Contiguous panel rows with equal directed width
share that positive integer width outside their exact dot product.
Positive and negative runs cache repeated exact integer square roots
on pairs of old endpoint columns. A zero companion corner contributes
zero in every guarded itinerary. Nonzero repeated companion corners
retain their original seven nested integer square roots. These
operations preserve the individual rounding and the certified lower
contribution. Each companion piece can also be located by an exact
integer lattice index. The symmetric degree-two density reuses equal
reflected weights and horizontal geometric means. The scalar ray
chord evaluates the same 16-bit directed tangent without packing
other itineraries. Exact integer dot products preserve the radial
matrix contraction term by term.
Negative mesh: 68 knots at scale 112, first positive knot 1/896; degree-one upper coordinates (15/28, 11/14, 1); degree-two upper coordinates (15/28, 15/14, 15/14). Use 20 degree-one temporal panels and 28 compressed-endpoint panels. The integer knots are exactly build_graded_grid() in root_integration.py. Both grids end at 3/2.
Table values in both replays use scale \(1/10^8\).
Each grid lookup is checked at the lattice boundaries by binary search.
In the degree-two density twice the negative log density, up to its additive constant, is
\[
 33x^2+12(y^2+yz+z^2)+24r x(y+z).
\]
All terms are increasing in the nonnegative magnitudes, so raising
coordinates and \(r\) gives an upper exponent and a lower density.
Both exact root integrals are symmetric in \(y,z\). The calculation
takes a diagonal cell once and a representative of an off-diagonal
pair twice, after taking the best union bound within that cell.

The half-normal and root-density prefactors are rounded downward
and certified by the outward Machin interval for \(\pi\). The two
integer bounds for \(r=\sqrt2\) and the lower one for \(1/r\) are
checked by squaring. Gaussian magnitudes are rounded upward to
\(1/40000\), the exponential argument upward to \(10^{-6}\).
The cubic Taylor polynomial for \(\exp(-x/256)\) is a lower
bound for nonnegative \(x\). Only positive lower values are squared,
eight times with downward divisions. Above the exponent cutoff
use zero. For an integer exponent numerator \(i=4096q+r\), multiply the
separate downward floors at \(4096q\) and \(r\), then divide down by the
common scale. The exponential product identity keeps the result downward.

Tabulate \(H\) on a separate \(1/2048\) mesh. On interval \(i\)
the exact mean of \(6h^2\) is
\(2(3i^2-3i+1)/2048^2\). Jensen for the exponential of the negative
argument gives a lower bound for the cell average. Add the integer
density numerators *before* dividing to form each cumulative node.
Between consecutive nodes, use linear interpolation of the two
lower cumulative values. The half-normal density decreases on the
positive half-line, so \(H\) is concave there; the chord is a lower
bound throughout. For \(t\ge2\) use the truncated \(H(2)\).

At a non-grid point a log-concave branch is bounded below by the
geometric interpolation of the two certified node values. The
implementation bisects a chord up to twice by downward integer
square roots, then uses a weighted harmonic mean, no larger than
the remaining geometric interpolation. On a single segment use
the smaller of the directed endpoint bounds; across segments the
true log-concave branch is at least the smaller of its two boundary
values. An out-of-grid lookup contributes zero. In the endpoint
and side integrations we may discard an initial or final piece
to put the sampled endpoint strictly inside the positive support
of \(H\). Integration limits are rounded inward, and all omitted
cells and dropped terms are nonnegative. The root density-weight cutoff is
\(200000/10^8\). Dropping the entire lower contribution of such a cell or
endpoint remains valid. Products, means and sums
for each type use integer lower arithmetic. Both hosted inputs return their exact root pair and their new
ray or center counts, as well as the same rational transfer.


For the root and guarded endpoint lookups at a fixed parent coordinate, bisect each
companion interval into four equal pieces, using downward geometric
midpoints. On a piece let \(m\le h\) be its two certified endpoints,
and \(q\) the fractional distance from the smaller endpoint. If \(m>0\),
the true value is at least \(m(h/m)^q\), which is at least
\[
 m+q\,2m(h-m)/(h+m).
\]
Indeed \(\log z\ge2(z-1)/(z+1)\) for \(z\ge1\), by differentiation,
and \(\exp v\ge1+v\). If \(m=0\) use zero. Round the slope down,
and round the minimum \(q\) on the argument interval down to a
multiple of \(2^{-16}\). When an interval crosses a piece boundary,
use the smaller of its two evaluated endpoint bounds by log-concavity.
These calculations for several itineraries use distinct 64-bit digits
of one integer. Each digit before rescaling is at most
\(3\cdot10^82^{16}<2^{64}\), so no addition or multiplication by
the common fraction numerator carries into the adjacent digit.
After shifting right sixteen places a mask discards the adjacent
digit's low bits. The resulting values are individual lower bounds.

## The reference crossing-path enclosure

On a three-regular tree the normalized critical wave has covariance
\(r_j=(1+j/3)2^{-j/2}\) across distance \(j\). Its path recurrence is
\(r_{j+1}=\sqrt2\,r_j-r_{j-1}/2\), with \(r_1=2\sqrt2/3\).
Thus \(X_{i+1}-\sqrt2 X_i+X_{i-1}/2\) is uncorrelated with
the preceding path values, including \(X_i\). In the joint Gaussian
law these innovations are independent of the past and have
variance \(1/12\). Negate alternating coordinates and put the
magnitudes in the positive orthant. The first-pair and transition
densities for one alternating pattern are
\[
 f_2(x,y)=c_2 e^{-9(x^2+y^2)/2-6\sqrt2 xy},\qquad
 f(z\mid x,y)=c_1 e^{-6(z+\sqrt2 y+x/2)^2},
 \qquad c_2=\frac3{2\pi},\quad c_1=\sqrt{\frac6\pi}.
\]
A directed \(j\)-edge path has two alternating patterns of identical
probability \(p_j\). Counting one half of the directed
nonbacktracking paths gives \(3\cdot 2^{j-1}p_j\) per vertex.
An \(\ell\)-edge crossing path contributes \((\ell-j+1)_+\)
to this half-count and has minimum crossing cover
\(\lceil\ell/2\rceil\). The alternating sum through an odd \(J\)
equals this value for \(\ell<J\), and equals
\(\ell-(J-1)/2\ge\lceil\ell/2\rceil\) for \(\ell\ge J\).
An even cycle contributes its length to each half-count, counting a
nonbacktracking walk from each vertex in either direction and halving.
The odd alternating sum equals its length, which bounds its cover
size of half its length.
After the vanishing crossing-degree-three correction described
above, it therefore suffices to bound
\[
 \Gamma_{25}=\sum_{j=1}^{25}(-1)^{j+1}3\cdot2^{j-1}p_j .
\]

### Tails outside the aperture

The aperture is \([0,L]^m\) with \(L=5/2\) and \(m=j+1\).
One more positive step has conditional probability at most \(1/2\)
because its conditional mean is nonpositive. For two steps, the
integrand consisting of the first positive transition density
times the second positive continuation probability decreases in
each of the two initial magnitudes. Its supremum occurs at zero
initial magnitudes. There, after dividing out the innovation
scale, the two independent standard Gaussian innovations occupy a
wedge of angle \(\arctan(1/\sqrt2)<\pi/5\). The strict inequality
follows from \(\cos(2\pi/5)=(\sqrt5-1)/4<1/3\).
Each arm of length \(k\) consequently costs at most
\(H(k)^{-1}\), where \(H(k)=10^{\lfloor k/2\rfloor}2^{k\bmod2}\).

Anchor the tail in coordinate \(i\), indexed from zero, at
the adjacent pair \((a,a+1)\), \(a=\min(i,m-2)\).
Reversibility and the second-order path recurrence separate the
outward positive arms conditional on this pair. The pair density
with one specified coordinate above \(L\) has integral at most
\[
 c_2\int_L^\infty e^{-9u^2/2}\frac{du}{6\sqrt2 u}
 \le \frac{c_2e^{-9L^2/2}}{54\sqrt2 L^2}.
\]
The first inequality integrates out the other coordinate after
dropping its quadratic penalty. The second uses \(1/u\le u/L^2\).
The pair density is symmetric, so this bound also applies when the
tail coordinate is the second in its anchor pair. Summing over \(i\)
with multipliers \(H(a)^{-1}H(m-2-a)^{-1}\) bounds the full omitted
part of one positive orthant.

### Tensor Gauss error on the aperture

The density on the positive orthant is
\(f_m(0)\exp(-x^\top Qx/2)\), where \(f_m(0)=c_2c_1^{m-2}\).
The precision matrix has bandwidth two with nonnegative entries,
and its diagonal is \((9,9)\) at \(m=2\), \((12,33,12)\) at
\(m=3\), and \((12,36,39,\ldots,39,36,12)\) at \(m\ge4\).
These entries follow by expanding the pair and transition
exponents.

Use a \(q=48\) point Gauss–Legendre rule on \([0,L]\).
Completing the square in derivative coordinate \(i\), its
\(2q\)-th derivative has absolute value at most
\[
 f_m(0) Q_{ii}^q\frac{(2q)!}{2^q q!}
\]
times the Gaussian of the remaining Schur complement. The
coefficient is the \(2q\)-th moment of a centered Gaussian with
variance \(Q_{ii}\), obtained by differentiating its Fourier
representation. At least \(\max(0,m-9)\) interior coordinates
are more than two indices from \(i\). Their diagonal \(39\)
and their nonnegative couplings to all remaining coordinates
are unchanged by the Schur subtraction. After keeping their
individual \(e^{-39x_k^2/2}\) factors, the other diagonal block
remains positive semidefinite and the discarded cross terms
are nonnegative on the aperture. Both the one-dimensional
integral and the checked Gaussian rule sum for
\(e^{-39x^2/2}\) are below \(1/4\). For the integral use
\(\sqrt{\pi/78}<1/4\).

The one-dimensional rule error is bounded by the sup norm of
the \(2q\)-th derivative times
\[
 \frac{L^{2q+1}(q!)^4}{(2q+1)((2q)!)^3}.
\]
Indeed, use Hermite interpolation at the \(q\) quadrature nodes.
The monic shifted Legendre polynomial has squared norm
\(L^{2q+1}/[(2q+1)\binom{2q}{q}^2]\), and the Hermite remainder
divides it by \((2q)!\). Telescope the \(m\) coordinates, using
positivity of all weights and preserving the inactive Gaussian
factors above in each term. The resulting tensor error is
\[
 \frac{f_m(0)L^{m+2q}(q!)^3\sum_i Q_{ii}^q}
      {(2q+1)2^q((2q)!)^2(4L)^{\max(0,m-9)}} .
\]

### Directed replay

The standalone arithmetic uses 100-digit Decimal contexts
with directed addition, multiplication and division.
For the two elementary functions, take the previous and
next representable decimals around the standard library's
correctly rounded square root and exponential. Machin's identity
\(\pi=16\arctan(1/5)-4\arctan(1/239)\) supplies its enclosure.
For each arctangent the even 80-term alternating sum is a
lower endpoint, and its first omitted term
\(1/(161a^{161})\) gives the upper endpoint for input \(1/a\).
Machin's identity follows by tangent addition with both sides
in the first quadrant after division by four.

The numerical initializers only propose centers for the nodes.
The program certifies a bracket of radius \(10^{-72}\) around
each center. For a bracket endpoint \(a/b\) the signs are exact
via the integer recurrence for \(F_k=k!b^kP_k(a/b)\):
\[
 F_0=1,\quad F_1=a,\quad
 F_k=(2k-1)aF_{k-1}-(k-1)^2b^2F_{k-2}.
\]
There are 48 disjoint intervals inside \((-1,1)\), each with
opposite signs of \(P_{48}\), which identifies all the roots.
The derivative identity
\((x^2-1)P'_q(x)=q(xP_q(x)-P_{q-1}(x))\) gives an
enclosure at each exact center. Since \(\|P_q\|_\infty\le1\)
on \([-1,1]\), applying Markov's inequality twice bounds
\(\|P''_q\|_\infty\) by \(q^4\). Multiply this bound by the
actual outward bracket radius to enclose the derivative at the
root. The mapped node is \((1+x)L/2\) and its weight is
\(L/[(1-x^2)P'_q(x)^2]\), all enclosed outward.

The tensor contraction works on pair states. Its step factors are
\[
 A_{xy}=e^{-3x^2/2-6\sqrt2 xy},\quad
 B_{xz}=e^{-6xz},\quad
 C_{yz}=w_z c_1 e^{-12y^2-6z^2-12\sqrt2 yz}.
\]
Each of \(A,B\) is at most one. The script proves
\(\sum_z w_zc_1e^{-6z^2}<1\), so every \(C\) row sum is below
one and the exact transition contracts total nonnegative mass.
The script also bounds the initial pair mass by one.
At scale \(S=2^{224}\), let \(D_A,D_B,D_C\) be the largest
widths of the fixed-point enclosures for each corresponding
table. Switching to the lower table endpoints loses at most
\((D_A+D_B+qD_C)/S\) per step on a state of total mass at most
one. This follows by changing one factor at a time, using
the other exact factors as upper bounds. The downward shifts
after multiplying by \(A\) and after completing each new pair
entry each lose at most \(q^2/S\) in total. Contraction propagates
the old error without increase. The initial enclosure debit is
2308 grid units and each of the 24 steps adds 4755, giving
116428 grid units. Both counts are recomputed from the outward
tables in the replay.

For each \(j\), subtract the tensor error from the lower
fixed-point endpoint. Add the tensor error and aperture tail
to the upper fixed-point endpoint. Accumulate with the
appropriate alternating endpoint for \(\Gamma_{25}\). The
program proves
\[
 \Gamma_{25}
 \le 0.130647566474459387249570528851780
 < U=\frac{1633094581}{12500000000}.
\]
The rounded displayed number is an upper bound on the
program's full directed upper endpoint.

## Two-phase interpolation for a cubic edge order

For a finite simple cubic graph \(G\), give every vertex \(v\) an anchor
bit \(b_v\) and every edge a bit. A vertex is split if its three incident
edge bits are unequal. Let the charge be \(\phi(k)\) at \(k\) disagreements
with its anchor, where
\[
 \phi(0)=0,\qquad \phi(1)=1,\qquad \phi(2)=1+\delta,\qquad \phi(3)=\lambda .
\]
The following argument uses \(3/7\le\lambda<1/2\) and
\(0\le\delta\le\lambda\). Set
\[
 A_* = 3/4-\lambda/2,\qquad D_*= A_*-\lambda,\qquad
 v_* = \frac{\lambda\delta}{1+\delta}.
\]
In addition assume
\[
 (D_*+3v_*/2)^2-4(A_*+v_*/2)v_*\le 0.                 \tag{P}
\]

**Unit-factor interpolation theorem.** If the minimum total charge is at
most \(M\) at every anchor word induced by the prefixes of one order on an
\(n\)-vertex graph, then \(\operatorname{cw}(G)\le M+o(n)\), with uniform
error in the graph and order. We use the two-phase grouping and retain
the endpoint penalty and an extra degree-one slack in its analysis.
The charge convention was used in the 4.8833 submission. [Charge source](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/92f3238c3d1a039578dd9fc5d91401728c48b18f/submission) · [Two-phase grouping source](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/237f72b0e798b7f65329cd0551d526e8f1e41e6c/submission/proof.md#L9-L62)

## The changing edges form one connected component

For a fixed anchor, the successive mismatch increments
\(1,\delta,\lambda-1-\delta\) decrease. With the anchor fixed at one,
increments as edge bits increase are the negatives of the increments
in reverse order, and again decrease. For a mixed pair consisting of
one edge and the anchor, let \(j\) be the number of ones in the other
two edges. The edge increment with anchor zero is \(\Delta_{j+1}\),
and with anchor one it is \(-\Delta_{3-j}\), with \(\Delta\) the three
mismatch increments. The mixed submodular inequality reduces to
\(\lambda-\delta\ge0\) for \(j=0,2\) and \(2\delta\ge0\) for \(j=1\).
All mixed discrete differences are therefore nonpositive. The total
charge is jointly submodular as a sum of these local charges. Applying
its lattice inequality to minimizing edge assignments at two arbitrary
anchor words, then bounding the meet and join charges by their minima,
proves submodularity of the minimized charge \(G(b)\).

For fixed anchors, minimizers are closed under meet, so choose the least minimizer. If \(b\le b'\), applying the joint lattice inequality to two endpoint optima \((b,x),(b',x')\) shows that their meet and join are optima for the earlier and later anchors. Leastness gives \(x\le x'\). The all-zero and all-one anchors give the unique all-zero and all-one edge states.

In one anchor step, let \(D\) consist of the changed edges and their endpoints. A component of \(D\) that misses the changed anchor has identical anchors and identical fixed incident bits at both endpoints of the step. Exchanging only that component's edges preserves all costs elsewhere. Optimality in both states forces the two local costs to be equal. Returning those edges to zero at the later endpoint would then give a strictly smaller optimum, contrary to leastness. Thus a nonempty \(D\) is connected and contains the changed anchor.

We work on \(D\), since splits outside it stay constant in this step. Put \(M' = M\) minus the number of these fixed splits. The local soft endpoint costs \(M_0,M_1\) are both at most \(M'\). Each local optimum is optimal on this edge set at its own anchors. Comparing the two edge states under either anchor and using the one-anchor change gives \(|M_1-M_0|\le\lambda\). If \(|V(D)|\le L\), flipping all its edges in any order adds at most \(L\) splits above the frozen outside count, so only the connected case \(|V(D)|>L\) needs the two phases.

All changed bits move from zero to one. Let \(P\) count vertices split at both endpoints, \(A\) those split only at the destination and \(B\) those split only at the source. A degree-two vertex of \(D\) has a single fixed bit, making it type \(A\) for fixed zero and type \(B\) for fixed one. A degree-one vertex may also be of type \(P\). The degree-three vertices of \(D\) move from unanimity zero to unanimity one. Partition these vertices into counts \(z_0,z_1\) according to whether they pay a mismatch at the source or destination. Omit the changed anchor so that \(z_0+z_1\) is their total up to one. Let \(a\) count type \(A\) vertices paying a source mismatch and \(b\) count type \(B\) vertices paying a destination mismatch, with the same omission. Let \(d_0,d_1\) be the counts of two-disagreement split stars at the
source and destination, among the vertices of \(D\). Then, up to a
uniform \(O(1)\) from the changed anchor,
\[
\begin{aligned}
 M_0&=P+B+\lambda(z_0+a)+\delta d_0,\\
 M_1&=P+A+\lambda(z_1+b)+\delta d_1,\\
 B-A&=\lambda(z_1-z_0+b-a)+\delta(d_1-d_0)+O(1).
\end{aligned}                                                   \tag{1}
\]

We may take \(B\ge A\) by reversing the step and complementing its bits and anchors if needed. The resulting order is reversed back at the end. Connectivity concerns the original changed-edge set and survives this normalization.

## Grouping changed edges at bounded cost

To order an edge set that changes from zero to one, give the edges group labels and assign each group one independent uniform priority. At time \(p\), the expected number of split vertices starts from the linear interpolation of the endpoint split counts. The extra contribution at a degree-two changing vertex whose edges belong to two groups is \(p(1-p)\). At a degree-three changing vertex with \(k\) distinct groups for \(k=2,3\), it is \(1-p^k-(1-p)^k=k p(1-p)\). One group gives no extra cost. A vertex with at most one changing edge also needs no extra budget. These formulas follow by inspecting the fixed bit at a degree-two vertex and the two unanimous outcomes at a degree-three vertex. We count one unit for each paid degree-two vertex and at most three for each paid degree-three vertex.

We use the following grouping construction on a subcubic changed-edge graph \(E\). Fix integers \(L,J\) with \(2\le L\le J\), and use one group for each connected component of size at most \(L\). For a larger component \(D_0\), write \(\beta(D_0)=|E(D_0)|-|V(D_0)|+1\). Let \(K\) range over the components induced by the vertices of degree three in \(D_0\). There is a grouping of all these larger components with groups of at most \(3J\) edges and split units bounded by

\[
 \sum_{D_0}\left(\beta(D_0)-\sum_{K\subseteq D_0}\beta(K)
             +\frac32\sum_{K\subseteq D_0,\ |K|>L}\beta(K)\right)
       +O(nL/J).                                               \tag{2}
\]

To see this, first eliminate cycles outside the induced components \(K\). Every cycle outside them contains a vertex of degree two. Cut that vertex in two, keeping one of its incidences at each copy. The alternate path on the cycle keeps the modified graph connected and the extra vertex lowers cycle rank by one. Continue until every edge outside the induced components is a bridge. The exact number of such cuts is \(\beta(D_0)-\sum_K\beta(K)\). Recombining each such vertex later can cost at most one unit.

Keep an induced component with at most \(L\) vertices intact. In each larger \(K\), delete a degree-three vertex whenever its removal leaves the current graph connected. Each deletion lowers cycle rank by two, and two deleted vertices cannot be adjacent since the later one would already have lost an incidence. The remainder has only cycles or bridges as its blocks. Otherwise a two-connected noncycle block has a vertex with three incidences in that block, and its deletion leaves the graph connected. Restore the three incidences of each deleted vertex by three separate pendant copies. They do not reconnect distinct blocks. This connected replacement has the same boundary as \(K\), since the removed vertices had all three edges inside it. Recombination costs at most three units per removed vertex, for at most \(3\beta(K)/2\) units in \(K\).

The expanded graph has at most \(3n\) vertices. It admits a tree decomposition with bags of size at most \(\max(L,3)\), using one bag for a small \(K\), and endpoint pairs for bridges and bags of three along cycles for the larger replacements. For a cycle, fix one vertex and take bags consisting of that vertex and each consecutive pair along the remaining path. Glue the decompositions at common vertices and along the bridges. In such a decomposition, assign each vertex to one containing bag. A weighted centroid bag separates the graph into pieces with at most half the vertices, since outside that bag the connected bag sets of surviving vertices cannot cross between centroid branches. Remove centroid bags recursively until all pieces have at most \(J\) vertices. At every dyadic size range above \(J\) the processed pieces are disjoint, so there are \(O(n/J)\) separators. Assign one group to all edges touching each surviving piece, with edges between deleted vertices grouped individually. The only additional split units come from separator vertices and are \(O(nL/J)\), which proves (2).

## Retire the components that can move before the large pieces

In the connected \(D\), let \(K_L\) be the union of the induced degree-three components with size greater than \(L\). Let \(z_L=|V(K_L)|\) and let \(d_L\) be the number of their boundary edges in \(D\). Delete their vertices and all incident edges, and discard isolated remaining vertices. Call the resulting graph \(D_0\). Its degree-three vertices are precisely the original degree-three vertices outside \(K_L\), because none of those is adjacent to \(K_L\).

In a component \(C_0\) of \(D_0\), write \(m\) for its degree-three count, \(s_A,s_B\) for the numbers of type \(A,B\) vertices missing one edge to \(K_L\), and \(l\) for its original degree-one count in \(D\). Set \(s=s_A+s_B\). The leaves of \(C_0\) number \(s+l\), so \(\beta(C_0)=(m-s-l)/2+1\). Let \(B_0\) count its type \(B\) vertices that have no missing edge. Use the grouping above on this component, with integer intrinsic split charge \(u(C_0)\) before the separator allowance. All its induced trivalent components are small. Thus \(0\le u(C_0)\le\beta(C_0)\), and \(u(C_0)=0\) when the component itself is small.

Let \(\mathcal F\) consist of all components with \(A(C_0)+u(C_0)\le B_0(C_0)\). The first phase flips the edges of \(\mathcal F\) by the groups above while every other changing edge remains zero. For a component in \(\mathcal F\), the \(B_0\) vertices have interpolated charge \(1-p\), the type \(A\) vertices have interpolated charge \(p\), and the \(s_B\) border vertices remain split. The net change from its source split count is bounded by \(p[A(C_0)+(1-p)u(C_0)-B_0(C_0)]\) plus its separator allowance. The sum is nonpositive before the allowance. The total expected number of splits in this phase is at most \(P+B+O(nL/J)\).

For the second phase, retain the flipped bits in \(\mathcal F\) and group the remaining changing graph \(D_1\) anew. At its two ends a type \(A\) vertex of \(\mathcal F\) is already split and a \(B_0\) vertex is already unanimous. A type \(B\) border vertex of \(\mathcal F\) is still split at the start and becomes unanimous at the second endpoint. All its original trivalents are already unanimous. Thus, if \(D_F=B_0(\mathcal F)-A(\mathcal F)\ge0\), the second-phase expectation is at most

\[
 E(p)=P+pA+(1-p)B-(1-p)D_F+p(1-p)S_1+o(n),             \tag{3}
\]

where \(S_1\) is the split-unit bound for the new groups. The error in this and the following estimates will be \(O(n/L+nL/J+1)\).

Write \(m_F\) for the trivalents in the retired components and \(m_R\) for those in the other components of \(D_0\). Write \(l_1\) for the leaf count of \(D_1\). Every component of \(D_1\) has more than \(L\) vertices. It contains either a large \(K\) or a large component of \(D_0\), since a small unretired component must attach to a large \(K\) in the connected \(D\). The large induced trivalent components of \(D_1\) are exactly the large \(K\) retained from \(D\). Summing the cycle-rank identity by degrees gives

\[
\sum_{D_1\text{ components}}\beta\le (z_L+m_R-l_1)/2+O(n/L),\qquad
\sum_{K\text{ large}}\beta(K)=(z_L-d_L)/2+O(n/L).
\]

In (2), discard the nonpositive contribution of small induced components. We obtain

\[
 S_1\le\tfrac34(z_0+z_1)-\tfrac14 H+o(n),
 \qquad H=3m_F+m_R+d_L+2l_1.                              \tag{4}
\]

We also have \(H+2D_F\ge 2(B-A)\). Here are the contributions to this last inequality. In a retired component with counts \(m,s_A,s_B,l\) as before, its share of \(H\) is \(3m+3(s_A+s_B)\), since its border vertices become leaves in \(D_1\). Adding its contribution \(2(B_0-A)\) to \(2D_F\) pays \(2(B-A)\) there. An unretired component satisfies \(B_0-A<u(C_0)\le(m-s-l)/2+1\). Integrality gives \(2(B_0-A)\le m-s-l\), so \(2(B-A)\le m+s-l\). This is bounded by its share \(m+s+2l\) of \(H\). The remaining nontrivalent vertices touch \(K_L\) on all their changed edges, since none occurs in a nonempty component of \(D_0\). They contribute three units per original leaf and two per original degree-two vertex through \(d_L+2l_1\). Each such share pays its own contribution to \(2(B-A)\). This proves

\[
 H\ge0,\qquad H+2D_F\ge 2(B-A).                            \tag{5}
\]

## The extra slack from leaves and split penalties

Put \(J_*=H+2D_F-2(B-A)\), nonnegative by (5), and let \(B_\ell\)
be the number of degree-one vertices of type \(B\) in \(D\).
A retired component's original type \(B\) leaves belong to its
\(B_0\), as their changed edge lies in that component. In total these
leaves number at most \(D_F+A(\mathcal F)\). In an unretired component
the proof of (5) supplies slack at least \(3l\), with \(l\) its original
leaf count. An original leaf outside these components contributes
three units through \(d_L+2l_1\) and requires at most two toward
\(2(B-A)\). The nonretired original type \(B\) leaves therefore number
at most \(J_*\). Consequently \(B_\ell\le D_F+A+J_*\).

Every type \(B\) vertex of changing degree two has a fixed one bit. It
has two changing zero bits at the source. Unless it pays a destination
mismatch counted by \(b\), its anchor is one and its source has two
disagreements. This gives \(d_0\ge B-b-B_\ell-O(1)\).
Use (1) and the leaf estimate to conclude
\[
 (1+\delta)d_0\ge
 \lambda(z_1-z_0+b-a)+\delta d_1-D_F-J_*-b-O(1).       \tag{6}
\]

## Pay both phases from the endpoint costs

The first phase fits \(M'+o(n)\) since it stays at or below \(P+B\)
up to sublinear error. For the second phase take \(p\ge1/2\) first.
The mean of \(M_0,M_1\) pays (3). Its split term exceeds the
interpolated split term by \((p-1/2)(B-A)\ge0\).
The trivalent mismatch allowance pays the extra
\(3p(1-p)(z_0+z_1)/4\) in (4), because \(3/16\le\lambda/2\).
The other allowances are nonnegative.

For \(0<p\le1/2\) use the interpolation \(I=(1-p)M_0+pM_1\) and
\[
 \theta=\min\{1,(4\lambda-3p)/(2\lambda p)\}.
\]
Equations (4) and (5) yield
\[
 S_1\le \tfrac34(z_0+z_1)-\tfrac{\theta}{2}(B-A)
             +\tfrac{\theta}{2} D_F-\tfrac{\theta}{4}J_*+o(n)
\]
by using weight \(\theta\) of the exact identity
\(H=2(B-A)-2D_F+J_*\) and weight \(1-\theta\) of \(H\ge0\).
Substitute (1) in this inequality and use (3). The coefficients of
\(z_0,z_1\) in the bound for \(E(p)-I\) are respectively
\[
 (1-p)[p(3/4+\theta\lambda/2)-\lambda],\qquad
 p[(1-p)(3/4-\theta\lambda/2)-\lambda].
\]
The remaining coefficients are
\[
\begin{array}{c|r}
 D_F & -(1-p)(1-\theta p/2)\\
 a   & -\lambda(1-p)(1-\theta p/2)\\
 b   & -\lambda p(1+\theta(1-p)/2)\\
 d_0 & -\delta(1-p)(1-\theta p/2)\\
 d_1 & -\delta p(1+\theta(1-p)/2)\\
 J_* & -\theta p(1-p)/4 .
\end{array}
\]
All these remaining coefficients are nonpositive. The choice of
\(\theta\) makes the \(z_0\) coefficient nonpositive. For \(\theta<1\)
the \(z_1\) coefficient is
\(p(3/2-3p/2-\lambda/p)\le p(3/2-\sqrt{6\lambda})\le0\).
For \(\theta=1\) it is \(p(D_*-A_*p)\), so only
\(0<p<D_*/A_*\) needs a further allowance. In that range
\(p<1/4\), since \(\lambda\ge3/7\) gives \(D_*/A_*\le1/5\).

For this last range, put
\[
 h=\min\left\{\frac{\delta(1-p/2)}{1+\delta},\frac p4\right\}
\]
and pay \(h(1-p)(1+\delta)d_0\) from the \(d_0\) allowance. Replacing
its negative by the upper bound from (6) subtracts
\(h(1-p)\lambda\) from the \(z_1\) coefficient.
If \(h=p/4\), the new coefficient is nonpositive since
\((1-p)(A_*-\lambda/4)\le A_*-\lambda/4\le\lambda\).
For \(h=\delta(1-p/2)/(1+\delta)\), the same conclusion follows
from (P), which states that the quadratic
\[
 (A_*+v_*/2)p^2-(D_*+3v_*/2)p+v_*
\]
has nonpositive discriminant and positive leading coefficient.
This quadratic is precisely \(h(1-p)\lambda-p(D_*-A_*p)\).

Here the \(J_*\) coefficient stays nonpositive because \(h\le p/4\).
For \(D_F\) and \(a\) use \(h\le1-p/2\). For \(b\), its increase
is \(h(1-p)(1-\lambda)\le p(1-\lambda)/4\le p\lambda\).
The \(d_1\) coefficient only decreases. The new \(z_0\) coefficient
is at most
\((1-p)[p(3/4+\lambda/2)+(p/4)\lambda-\lambda]\le0\),
using \(p<1/4\) and \(3/7\le\lambda<1/2\). The unused \(d_0\)
coefficient stays nonpositive. These checks give \(E(p)\le I+o(n)\)
throughout the second phase. At \(p=0\) the first-phase endpoint
has the required bound. Restoring the frozen outside splits at
unit cost gives \(M+o(n)\).

## Realize the whole order simultaneously

Use \(L=\max(2,\lceil n^{1/16}\rceil)\) and \(J=L^2\). At the \(L+1\) equally spaced mesh times in each phase, the split counts have variance \(O(nJ)\). Indeed resampling a group's priority can change a count by at most twice its edge count, and the sum of squared group sizes is \(O(nJ)\). The number of edges flipped in each mesh interval also has variance \(O(nJ)\) and expectation \(O(n/L)\). A union bound and Chebyshev give one positive-probability realization for both phases in which all these counts differ from their means by at most \(n^{3/4}\), since the failure probability is \(O(LJ/\sqrt n)=o(1)\). Between mesh times, the number of additional splits is at most twice the number of edges flipped in the interval. Thus every edge prefix in both phases obeys the bound with a uniform sublinear error. Within a group the edges may be ordered arbitrarily, since the whole group belongs to a single mesh interval almost surely. Fix one successful choice of the two phase orders, for instance the first passing choice under any fixed ordering of the finite permutations of groups. The end of the first phase equals the start of the second as an actual edge state.

Choose one such order for each step of the fixed nested least-minimizer chain and concatenate. Each original edge flips exactly once. Order the original vertices by the index of their median incident edge, with arbitrary tie breaking. For a vertex prefix, let \(j\) be its largest median index. A crossing edge with index below \(j\) is the unique lowest incident edge of its excluded endpoint and charges that endpoint as split just before \(j\). A crossing edge above \(j\) charges its included endpoint as split just after \(j\), using its unique highest edge. These two sets of endpoints are disjoint, and the two split sets differ only at the endpoints of edge \(j\). The edge \(j\) itself contributes at most one. The vertex prefix therefore incurs at most three more crossing edges than the edge-prefix split bound. Together with the empty prefix this proves the bound on the resulting single order.

For \(\lambda=47091/100000\) and \(\delta=2102599/1000000000\), the left side
of (P) is \(-105331260465423/100420961892255480100000\). The replay checks
this rational inequality and the exact width and projected coefficient above.
