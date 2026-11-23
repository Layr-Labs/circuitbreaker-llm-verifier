# An eventual 4.89223N gate lower bound for the norm-one language

The exact experiment gives four disjoint saving events with total
probability at least \(p=5827595462269497/1342177280000000000\).
The crossing-cover expectation is at most
\(U=13064756648/10^{11}\). At \(\lambda=49257/100000\), the two-phase
conversion uses \(R=4411483561/4410931512\). The uniform simple cubic
cutwidth constant satisfies
\[
 a_0=R[U-(1-\lambda)p]
   =\frac{76051814369729864129395244548969}{592025205904244736000000000000000}
   <\frac{50000}{389223}=\frac{1}{2(4.89223-1)}.
\]
The companion family proof transfers any fixed constant strictly
between these width values to the asserted eventual gate count.

The conversion builds on [PR 891](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/237f72b0e798b7f65329cd0551d526e8f1e41e6c/submission/proof.md).
The path switches use the guards of [PR 890](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/92f3238c3d1a039578dd9fc5d91401728c48b18f/submission/proof.md).
We use the graph proof, norm-one transfer and directed certificate of the
[4.89153 submission](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/4321d1d8c25ed11c15bdbfad2d367ecc8efd17b5/submission).
The larger saving comes from certifying success on either of two opposite
arms, including an already-covered negative neighbor, extending the
finite branch cases and integrating their log-concave probabilities
more accurately. The proof and computation follow below.

## Cover descent and a physical anchor order

A star costs one when its incident edge signs split, zero when
all match its anchor, and \(\lambda\) for unanimous disagreement.
Let \(G(b)\) denote the minimum over edge signs at anchor word \(b\).
The interpolation proof below establishes its submodularity.
Complementing all signs preserves the cost. Changing one anchor
changes \(G(b)\) by at most \(\lambda\).

For a degree-\(d\) finite cover the minimum for repeated anchors is
\(d\) times the base minimum. A lift shows one inequality. For the
other, express each charge as the minimum over an auxiliary vertex
bit of nonunanimity of that bit and the three incident edges plus
\(\lambda\) times the auxiliary-anchor mismatch. Average the bits of
an optimum upstairs in their respective fibers. Threshold all these
averages at one uniform point in \([0,1]\). The expectation of each
nonunanimity charge is the range of the four averages, at most the
average of their fiber ranges. The expected anchor charge equals its
fiber average. Thus the induced base optimum is at most the
normalized optimum upstairs.

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
expected cost at most \([U-(1-\lambda)p+\epsilon]n\).
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
double-covered boundary edge the negative sign. Flip any equal-sign
edge needed by a witnessed covered vertex whose equal-sign
neighbors are all covered. An uncovered vertex's incident signs
agree with its anchor. A witnessed covered vertex pays \(\lambda\)
when all its crossing edge signs disagree with its anchor. This
can include the switched negatives except the last one at a
double-covered boundary. The old cover cardinality minus
\((1-\lambda)\) per witnessed vertex therefore upper bounds the
new cost, with the negligible deletion debit already specified.

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
for \(m<0\), zero otherwise. The family \(A_0\) consists of \(SB_i\)
for \(i\in I\). Add \(D(B_i,B_j)\) for \(i=1,3\) and \(j=0,2,4\)
to get \(A\). These added cases have the more negative arm of the
other parity, so they are disjoint from \(A_0\) and each other.
A positive vertex with an odd terminal arm belongs to the old
cover: in an even path it is in the smaller class, and in an odd
path positives are chosen.

The chosen outside family is \(U=H(m)\), the three
\(W_j=S_1(B_j)\) for \(j=1,3,5\), and \(D(B_1,B_1)\).
They are distinguished by opposite degree and closing lengths.
The chosen unsaved guards are the six pairs
\((i,o)=(0,U),(2,U),(4,U),(0,W_1),(2,W_1),(0,W_3)\).
Use \(S_1(B_i)\) with the additional multiplicative guard
\(o(b,2m+d)\) on the positive sibling. The analytic \(U\)
formula is valid at every nonnegative companion magnitude,
including outside the tabulated square.

At an opposite entry to the positive end of an odd path, put
\(M=rb+a/2\). For the allowed ordered outside/unsaved pairs
defined above, integrate \(f(b,M-h)k(b,M+h)\) against \(g(h)\)
on \(0\le h\le\min(2,M)\), and call this \(T_{f,k}\).
At an internal hop across a positive vertex, multiply the
\(K\) integrand by the outside guard \(o(b,2rb+a+d)\) and call
the result \(K'_o\). Use \(o=U,W_1,W_3\).
The extended odd-arm family \(E\) uses all endpoint pairs
and the first four iterates of the two-step hop \(K'_oK\).
All lengths and side guard choices distinguish the itineraries.
For a same-negative neighbor, its unique opposite child can be
the start of a completely switched odd path. A subfamily of \(E\)
using 30 distinct ordered endpoint pairs and 12, 3 and 2 continuations
at one, two and three \(U\)-guarded hops is indexed in the certificate.
Indices 0 through 4 denote the outside guards in the order just given.
Indices 5 through 10 denote the six unsaved guards in their displayed order.
Each selected pair has an outside member. Applying \(S_1\) to these
47 distinct itineraries gives the family \(P\).
For an even far arm use
\[
 O_g=\{K'_UB_j:j\in I\}\cup\{K'_{W_1}B_j:j=0,2,4\}
          \cup\{K'_{W_3}B_0\}.
\]
The extra first-positive guard allows a further pair to switch.
The negative root in the mixed near-odd/far-even case therefore
does not lie on the double-covered boundary.

For crossing degree one fix the root sign and magnitude \(x>0\).
With the opposite neighbor of magnitude \(b>0\) their joint
density is
\[
 \frac{3}{2\pi}\exp[-9(x^2+b^2)/2-6rbx].
\]
The two same-sign magnitudes have sum \(s=b+2rx\). With
\(h=ts/2\), \(0\le t\le1\), these magnitudes are \(s/2+h\)
and \(s/2-h\). The extra density in \(t\) is \((s/2)g(h)\).
The positive multiplier is
\(A(x,s/2+h)A(x,s/2-h)B_I(x,b)\).
For the negative root there is an additional way its smaller same-sign
neighbor can stay covered. If that vertex has two opposite children,
and each has no further opposite child, their crossing component is
the two-edge path with the negative vertex at its center. The old
cover contains that center. The switching rule only changes odd paths,
so the center stays covered. Its conditional witness has probability
\(T=D(B_0,B_0)\) under the stated truncations. Every member of \(P\)
requires exactly one opposite child, making \(T\) disjoint from \(P\).
The midpoint of the forward children of the larger companion is
\(r(s/2+h)-x/2>0\), so this two-child witness is needed only on
the smaller companion. The negative multiplier we integrate is
\(P(x,s/2+h)[P(x,s/2-h)+T(x,s/2-h)]E(x,b)\).
Integrate on \([0,23/32]\times[0,5/4]\times[0,1]\),
with a factor three for the location of the opposite neighbor.

For crossing degree two fix a root sign, its magnitude \(x\),
and two opposite magnitudes \(y,z\). The first opposite neighbor
conditional on \(x\) has mean \(2rx/3\) in the root sign
convention and variance \(1/9\); the other uses the independent
contrast above. The resulting joint density of the three
magnitudes is
\[
 \frac{\sqrt{108}}{(2\pi)^{3/2}}\exp[-x^2/2
       -(9/2)(y+(2r/3)x)^2-(3/2)(2z+y+2rx)^2].
\]
The same-sign magnitude is \(b=2rx+y+z\). For a positive root it
suffices for either opposite arm to certify an odd terminal distance.
If the crossing path is even, both root-to-end distances are then odd
and the root is in the smaller parity class. If the path is odd, its
positive vertices form the old cover. Requiring a covered same-sign
neighbor gives an old saving that survives the switches as shown above.
Conditioned on the root star the two arm tests are independent. With
\(u=B_I(x,y)\) and \(v=B_I(x,z)\) the positive multiplier is
\(A_0(x,b)(u+v-uv)\). The negative multiplier is
\[
 P(x,b)[B_I(x,y)B_I(x,z)+E(x,y)O_g(x,z)+O_g(x,y)E(x,z)].
\]
In the first negative case both arms are short certified odd
terminal arms, so the negative is in the smaller parity class
of the even path. In the mixed cases the odd near arm and the
guarded even far arm switch past the negative. The brackets are
disjoint by arm parity. Integrate on \([0,17/32]\times[0,17/16]^2\)
with a factor three for the location of the same-sign neighbor.
The two root signs and the two degrees distinguish the four
events, each of which certifies a vertex saving in the same
simultaneously switched cover. Interchanging the two named opposite
neighbors preserves both integrals, including the rectangular domain.
Expanding the negative log density above, apart from its normalization,
gives \([33x^2+12(y^2+yz+z^2)+24rx(y+z)]/2\).
The replay integrates the diagonal grid cells once and one representative
of each off-diagonal pair twice.

## Directed lower integration

A fixed itinerary has jointly log-concave probability as a
function of the conditioned scores. Its innovation density is
Gaussian and the specified signs, ordered magnitudes and
truncations are affine constraints. Marginalization preserves
log-concavity. One proof for a one-dimensional integration uses
the increasing quantile transport \(L\) between two normalized
sections. At the intermediate points \((1-\theta)t+\theta L(t)\)
apply log-concavity. The Jacobian \(1-\theta+\theta L'(t)\)
is at least \(L'(t)^\theta\), giving the geometric interpolation
of the two original masses as a lower bound. Approximation
extends this calculation to zero sections.

A log-concave integrand on an interval integrates to at least
the width times the geometric mean of its two endpoint values.
Apply log-concavity along the interval and the arithmetic-geometric
inequality to the exponential interpolant. Apply this successively
to all three coordinates in a root cell. In the first root
formula the arguments and the nonnegative Jacobian are affine
in each coordinate separately when the other two are fixed.
Separate log-concavity suffices. Take each eight-corner mean
before summing the disjoint itineraries. Independent factors
in the same mean can be processed independently and multiplied.

Two univariate refinements improve the branch tables. For positive
lower endpoints \(a,b\), the integral of their geometric interpolant
on the unit interval is \(\sqrt{ab}\sinh(t)/t\), where
\(t=|\log(b/a)|/2\) uses the natural logarithm. At \(a=b\) use its
continuous value. Since \(\tanh(t)=|b-a|/(a+b)\le t\) and
\(\sinh(t)/t\ge1+t^2/6\), a directed lower average is
\[
 \sqrt{ab}\left[1+\frac{1}{6}\left(\frac{a-b}{a+b}\right)^2\right].
\]
The average is taken on the endpoint products for a single itinerary
in \(D,S,S_1\) and the point endpoint integrals. A zero endpoint gives
the safe lower average zero. For a point hop \(K\) or \(K'_o\), split
each integration panel in half. Compute lower values of the geometric
chord at the two half-panel midpoints, separately for the kernel and
the branch. The two chord products are values of the product chord at
one-quarter and three-quarters. Jensen on each half bounds its integral
below by its width times that value. All square roots and divisions in
these averages are rounded down.

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

For the extended \(E\) integrals, each of the four parameter
corners has its own cap \(C_s\le M_s\), with \(h=C_su\).
Use sixteen equal \(u\) panels. Four-term Hölder bounds the
four-corner mean of each integral below by the sum of the means
of its panel contributions. The cap contributes its Jacobian.
On one panel the density and the two guard products factor in
the mean. The sum over disjoint admitted endpoint pairs splits
as outside times all plus unsaved times outside. For a hop
integral the old argument \(f(b,d)\) appears at each of its four
corners twice in the eight-corner kernel evaluation. Its mean
is exactly the four-corner input mean, giving a nonnegative
kernel on these cell means. Apply this first to \(K\) and then
to the sum of \(K'_o\) applied on the disjoint guard branches.
Repeat four times and add the disjoint lengths.

The implementation uses coordinate scale \(2^{24}\) and value
scale \(10^{10}\). Integer square root supplies an outward
\(\sqrt2\) bracket. Exact Machin bounds on \(\pi\) verify all
density constants. Increase each exponent before using the
alternating cubic lower truncation at \(1/256\) of that exponent
and eight downward squares. Values past the cutoff give zero.
The cumulative \(H\) uses full \(1/2048\) cells and Jensen's
inequality in the square variable with its exact cell mean. Between
cell boundaries interpolate the two lower cumulative masses linearly.
The density decreases on the positive half-line, so its cumulative is
concave there and this interpolant stays below it. The degree-one
root integral uses 24 equal panels in \(t\).
For other arguments in a grid interval, geometric interpolation
uses downward square-root bisection and then a weighted harmonic
mean, which is below the remaining weighted geometric mean.
A monotone interpolant gives the smaller endpoint across one
segment. Across separate segments log-concavity of the true
function bounds the whole interval by the smaller of the two
certified endpoint values. Drop any wholly or partially
out-of-grid lookup. The analytic \(U\) and \(H\) branches remain
valid outside the grid. Every square root, interval-domain cut
and mass division has the safe direction.

For the many root lookups at a fixed parent coordinate, bisect each
companion interval into sixteen equal pieces using downward geometric
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
\(3\cdot10^{10}2^{16}<2^{64}\), so no addition or multiplication by
the common fraction numerator carries into the adjacent digit.
After shifting right sixteen places a mask discards the adjacent
digit's low bits. The resulting values are individual lower bounds.
In the point endpoint and side integrations, skipping panels outside
the positive support of either endpoint table removes only zero lower
contributions.

The replay returns each individual rational lower integral and
uses their sum in the combination above, independently of
the crossing-cover computation.

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

The aperture is \([0,L]^m\) with \(L=3\) and \(m=j+1\).
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

Use a \(q=60\) point Gauss–Legendre rule on \([0,L]\).
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
There are 60 disjoint intervals inside \((-1,1)\), each with
opposite signs of \(P_{60}\), which identifies all the roots.
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
3606 grid units and each of the 24 steps adds 7503, giving
183678 grid units. Both counts are recomputed from the outward
tables in the replay.

For each \(j\), subtract the tensor error from the lower
fixed-point endpoint. Add the tensor error and aperture tail
to the upper fixed-point endpoint. Accumulate with the
appropriate alternating endpoint for \(\Gamma_{25}\). The
program proves
\[
 \Gamma_{25}
 \le 0.130647566471217912420492299889353
 < U=\frac{1633094581}{12500000000}.
\]
The rounded displayed number is an upper bound on the
program's full directed upper endpoint.

## Two-phase interpolation for a cubic edge order

For a finite simple cubic graph \(G\), give each vertex \(v\) an anchor \(b_v\in\{0,1\}\) and each edge \(e\) a bit \(x_e\in\{0,1\}\). A vertex is split when its incident edge bits are unequal. Charge a split vertex one, a unanimous vertex with bit different from its anchor \(\lambda\), and every other vertex zero. Let \(Q_b(x)\) be the sum of these charges and let \(g(b)=\min_x Q_b(x)\). A vertex order determines a sequence of anchors by turning the vertices of each prefix to one.

**Interpolation theorem.** Fix \(3/8\le\lambda\le1\). If the maximum of \(g(b)\) along some anchor order of an \(n\)-vertex graph is at most \(M\), then

\[
\operatorname{cw}(G)\le M\left(1+\frac{9(1-2\lambda)_+^2}{16\lambda(3-2\lambda)}\right)+o(n).
\]

Here cutwidth is the smallest maximum number of crossing edges at a vertex prefix. The error is uniform in the graph and anchor order. The proof constructs one edge order and then one vertex order. We reconstruct the two-phase converter from the pinned #891 proof and give the connectedness argument for its transitions. [Source](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/237f72b0e798b7f65329cd0551d526e8f1e41e6c/submission/proof.md#L9-L62).

## The changing edges form one connected component

Introduce a Boolean variable \(s\) at each vertex. Its charge is the minimum over \(s\) of the indicator that \(s\) and the three incident edge bits are not all equal, plus \(\lambda|s-b_v|\). Indeed a split among the three edges always costs one, and a unanimous triple admits the cheaper charge \(\lambda\) for a mismatch. The indicator is submodular since its values by Hamming weight are \(0,1,1,1,0\) and its successive increments decrease. The mismatch cost is also submodular. Minimizing a submodular function over auxiliary Boolean variables preserves submodularity, by applying its lattice inequality to two minimizing completions. Thus the total is jointly submodular in the anchors and edge bits.

For fixed anchors, minimizers are closed under meet, so choose the least minimizer. If \(b\le b'\), applying the joint lattice inequality to two endpoint optima \((b,x),(b',x')\) shows that their meet and join are optima for the earlier and later anchors. Leastness gives \(x\le x'\). The all-zero and all-one anchors give the unique all-zero and all-one edge states.

In one anchor step, let \(D\) consist of the changed edges and their endpoints. A component of \(D\) that misses the changed anchor has identical anchors and identical fixed incident bits at both endpoints of the step. Exchanging only that component's edges preserves all costs elsewhere. Optimality in both states forces the two local costs to be equal. Returning those edges to zero at the later endpoint would then give a strictly smaller optimum, contrary to leastness. Thus a nonempty \(D\) is connected and contains the changed anchor.

We work on \(D\), since splits outside it stay constant in this step. Put \(M' = M\) minus the number of these fixed splits. The local soft endpoint costs \(M_0,M_1\) are both at most \(M'\). Each local optimum is optimal on this edge set at its own anchors. Comparing the two edge states under either anchor and using the one-anchor change gives \(|M_1-M_0|\le\lambda\). If \(|V(D)|\le L\), flipping all its edges in any order adds at most \(L\) splits above the frozen outside count, so only the connected case \(|V(D)|>L\) needs the two phases.

All changed bits move from zero to one. Let \(P\) count vertices split at both endpoints, \(A\) those split only at the destination and \(B\) those split only at the source. A degree-two vertex of \(D\) has a single fixed bit, making it type \(A\) for fixed zero and type \(B\) for fixed one. A degree-one vertex may also be of type \(P\). The degree-three vertices of \(D\) move from unanimity zero to unanimity one. Partition these vertices into counts \(z_0,z_1\) according to whether they pay a mismatch at the source or destination. Omit the changed anchor so that \(z_0+z_1\) is their total up to one. Let \(a\) count type \(A\) vertices paying a source mismatch and \(b\) count type \(B\) vertices paying a destination mismatch, with the same omission. Then, up to a uniform \(O(1)\) from the changed anchor,

\[
M_0=P+B+\lambda(z_0+a),\qquad M_1=P+A+\lambda(z_1+b),
\qquad B-A=\lambda(z_1-z_0+b-a)+O(1).               \tag{1}
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

## Pay both phases from the endpoint costs

The first phase already fits \(M'+o(n)\). In the second phase for \(p\ge1/2\), the mean of \(M_0,M_1\) pays (3). Its split term exceeds the interpolated split term by \((p-1/2)(B-A)\ge0\). The mismatch allowance for the trivalent vertices pays the possible extra \(3p(1-p)(z_0+z_1)/4\), since \(3p(1-p)/4\le3/16\le\lambda/2\). The other terms in (3) and (4) can only lower the requirement. All statements allow the uniform residual error.

For \(0<p\le1/2\) use the interpolation \((1-p)M_0+pM_1\). Combine weight \(\theta\) of the second inequality in (5) with weight \(1-\theta\) of its first one and use (1) in (4). This gives

\[
 S_1\le (3/4-\theta\lambda/2)z_1+(3/4+\theta\lambda/2)z_0
        +\theta\lambda(a-b)/2+\theta D_F/2+o(n).
\]

Subtract the endpoint interpolation from (3). The coefficients of \(D_F,a,b\) are respectively

\[
 (1-p)(\theta p/2-1),\quad -\lambda(1-p)(1-\theta p/2),\quad
 -\lambda p(1+\theta(1-p)/2),
\]

all nonpositive. Set \(\theta=\min\{1,(4\lambda-3p)/(2\lambda p)\}\), which is nonnegative in the stated range. The coefficient of \(z_0\) is then nonpositive. If \(\theta<1\), the coefficient of \(z_1\) is \(p(3/2-3p/2-\lambda/p)\), which is at most \(p(3/2-\sqrt{6\lambda})\le0\). If \(\theta=1\), that coefficient is at most

\[
\max_{p\ge0}p\left[(1-p)(3/4-\lambda/2)-\lambda\right]
 = \frac{(3/4-3\lambda/2)_+^2}{4(3/4-\lambda/2)}.
\]

Since \(\lambda z_1\le M_1+O(1)\le M'+O(1)\), the extra multiplier is at most \(9(1-2\lambda)_+^2/[16\lambda(3-2\lambda)]\). At \(p=0\) the end of Phase 0 supplies the same bound. Restoring the frozen splits outside \(D\) at unit cost gives the asserted multiplier on \(M\).

## Realize the whole order simultaneously

Use \(L=\max(2,\lceil n^{1/16}\rceil)\) and \(J=L^2\). At the \(L+1\) equally spaced mesh times in each phase, the split counts have variance \(O(nJ)\). Indeed resampling a group's priority can change a count by at most twice its edge count, and the sum of squared group sizes is \(O(nJ)\). The number of edges flipped in each mesh interval also has variance \(O(nJ)\) and expectation \(O(n/L)\). A union bound and Chebyshev give one positive-probability realization for both phases in which all these counts differ from their means by at most \(n^{3/4}\), since the failure probability is \(O(LJ/\sqrt n)=o(1)\). Between mesh times, the number of additional splits is at most twice the number of edges flipped in the interval. Thus every edge prefix in both phases obeys the bound with a uniform sublinear error. Within a group the edges may be ordered arbitrarily, since the whole group belongs to a single mesh interval almost surely. Fix one successful choice of the two phase orders, for instance the first passing choice under any fixed ordering of the finite permutations of groups. The end of the first phase equals the start of the second as an actual edge state.

Choose one such order for each step of the fixed nested least-minimizer chain and concatenate. Each original edge flips exactly once. Order the original vertices by the index of their median incident edge, with arbitrary tie breaking. For a vertex prefix, let \(j\) be its largest median index. A crossing edge with index below \(j\) is the unique lowest incident edge of its excluded endpoint and charges that endpoint as split just before \(j\). A crossing edge above \(j\) charges its included endpoint as split just after \(j\), using its unique highest edge. These two sets of endpoints are disjoint, and the two split sets differ only at the endpoints of edge \(j\). The edge \(j\) itself contributes at most one. The vertex prefix therefore incurs at most three more crossing edges than the edge-prefix split bound. Together with the empty prefix this proves the bound on the resulting single order.

For \(\lambda=49257/100000\), direct substitution in the theorem gives
the multiplier \(4411483561/4410931512\).
