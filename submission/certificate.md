# Directed numerical evidence for the norm-one circuit bound

The notation, graph construction and circuit theorem are in [the proof](proof.md). The branch integration, center adjustments, radial enclosure and covariance estimates here are reproduced from [the predecessor's proof and family supplement](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/a49ccfaf5eb17002ea7edfb90c9fe76bcc9e2c82/submission).

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
Apply it to the 18 panels at the four corners of each endpoint
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

For the positive base replay, use 103 vertices on the scale-320 grid.
The lattice scale is \(1/(320\cdot10^8)\), with first positive
node \(1/2560\). The degree-one positive box ends at
\(9/16,63/80,1\) and the degree-two positive box ends at
\(41/80,87/80,87/80\). Use 40 panels in the degree-one
\(t\) coordinate. The positive replay uses integer transition vectors, releases its
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
The scale-104 grid has 60 vertices and first positive node \(1/416\).
For the degree-one and degree-two negative base replays the upper
coordinates are \((7/13,10/13,1)\) and \((7/13,27/26,27/26)\). It uses 16 panels
in the degree-one \(t\) coordinate and 18 panels in the compressed endpoint
temporal coordinate. Both graded grids end at \(3/2\).
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

### Positive rays between a center and guarded negative endpoints

A positive center of crossing degree zero can also have the
following same-positive arm. Its neighbor \(p\) and the next
vertex \(q\) each have exactly one crossing neighbor, a negative
endpoint of a single-edge crossing component. The other
same-positive neighbor of \(q\) is outside the original cover
by either zero crossing degree or a single crossing arm of
even length two. Before switching, both \(p,q\) are in the
positive cover and each has an outside same-positive neighbor,
the center or the tail respectively. They are unsaved guards
for one another. Thus both single-edge crossing paths switch,
and \(p,q\) become uncovered while their negative endpoints
enter the cover.

Impose this additional safety guard at each of the two negative
endpoints. Each same-negative neighbor has crossing degree
zero and has at least one other same-negative neighbor of
crossing degree zero. The endpoint therefore has charge one
before the ray adjustments. Its two same-sign neighbors are
outside the cover, so neither the base savings nor any
secondary target flip touches it. It is ineligible as a
degree-one target for either neighbor and neither neighbor
has all three neighbors covered. In particular the negative
center adjustments also leave it untouched.

Select a positive center when at least two incident arms are
available, counting both ordinary stable degree-two targets
and these rays. Use the stable target adjustment as before.
For a ray assign the negative label to the edge from the center
to \(p\), to \(p q\), and to the two crossing edges of \(p,q\).
The negative endpoints fall from charge one to zero. If a ray
is selected only in that direction then \(p\) pays \(\lambda\)
and \(q\) pays \(1+\delta\). Apart from the center this saves
\(1-\lambda-\delta\). The same pair \(p,q\) can also be selected
from the other end. Now both pay \(\lambda\), and the combined
decrease \(2-2\lambda\) pays both arm credits. Distinct ray
pairs are vertex-disjoint because each participating positive
vertex has exactly one same-positive neighbor of crossing
degree one on a single-edge odd crossing path. Its other
same-positive neighbor is the center or the guarded outside
tail. Negative endpoints are unique by crossing degree one.
None of the ray vertices is a stable target.

Each ordinary degree-two arm saves \(1+\delta-\lambda\).
If exactly two arms are used, and \(k\) of them are stable
degree-two targets, the center pays \(1+\delta\). Together
the decrease is at least \(\beta_k\). If all three arms are
used, the center pays \(\lambda\), and the same three-pair
calculation above shows that all three pair credits fit.
Stable pairs keep their existing credit. Add the counts
\(\sigma_0,\sigma_1\) of pairs of two rays or of a ray and
one stable degree-two target with weights \(\beta_0,\beta_1\).

To certify these counts, let \(H\) be the half-normal mass
function used above. For a same-sign entry \((a,b)\), write
\(m=rb-a/2\) and \(U_0(a,b)=H(\max(0,m))\). To witness a
neighbor of the negative endpoint, use
\[
 G_0(a,b)=\int_0^m g(h)U_0(b,m+h)\,dh
\]
for \(m>0\), zero otherwise. This certifies both children
same-sign and the larger child also of crossing degree zero.
For an opposite entry set \(M=rb+a/2\). The negative endpoint
with both guarded same-negative neighbors is certified by
\[
 B_*(a,b)=\int_0^M g(h)G_0(b,M+h)G_0(b,M-h)\,dh.
\]
Use either \(U_0\) or the distinct outside-even-path guard
\(S_1 B_1\) for the tail. If its table is \(T\), the guard at
\(q\) is
\[
 Q_T(a,b)=\int_{|m|}^\infty g(h)T(b,h+m)B_*(b,h-m)\,dh.
\]
The guard at \(p\) has the same integral with \(Q_T\) in
place of \(T\). Each fixed tail history gives a log-concave
event by the Gaussian branch argument and the linear sign,
ordering and path predicates. Stable arms use all four
\(D(B_i,B_j)\) with \(i,j\in\{0,2\}\) for this additional
calculation.

The replay uses a separate graded scale-48 lattice with half-unit
steps up to six units, unit steps to 24, steps of two to 44, and
steps of four from 48 to 72. It uses 16 panels in each inward half-normal
aperture and estimates the joint star integral through the last
lattice vertices at \(3/2\) for the center and both selected neighbors. At table nodes the child magnitudes
are bracketed using the two rational square-root bounds.
The upper contrast aperture is \(5/4\). All allowed ranges are
rounded inward. For the one-opposite integrals, the lower
endpoint is the larger absolute midpoint bracket. For two
same-sign children it uses the lower positive midpoint as an
upper integration limit, and for two opposite children it uses
the negative of the upper midpoint. Geometric endpoint means
bound each panel integral below. The mass table uses a 1024-mesh
cumulative lower Riemann-Jensen sum and concave interpolation.
The same symmetric star polynomial allows the ray replay to pair reflected cells, with both mixed-arm assignments counted. The lower table for each fixed history is used separately for
the four-corner mean of each target rectangle. Multiply the two
neighbor factors and the joint star density over boxes whose
third neighbor is positive. Sum with factor three for the
unordered position pair and factor two for the two roles in
a mixed pair. These are lower expectations on the pair counts,
including multiplicity when all three arms qualify.

### Directed calculation for the centers

Put \(r=\sqrt2\), \(m=rv-x/2\). The stable degree-two target uses
\(D(B_i,B_j)(x,v)\) for \(i,j\in\{0,2,4,6,8\}\), separately and
then summed. For negative targets add both orientations of all
\(74\cdot12\) near/far combinations from the banks above.
Both degree-two target classes have support \(v<x/(2r)\).

For a mixed target write \(\mu=-m>0\) and set \(h=t\mu\)
in the two-opposite-child integral. The half-normal contrast
has \(0\le t\le1\), since all used \(\mu\) are below the
contrast aperture. The two child magnitudes are
\(\mu(1+t),\mu(1-t)\) and the density factor is \(\mu g(t\mu)\).
Use either the near history at the larger child and the far
history at the smaller, or the reverse. The ordered parity
classes are disjoint. Each factor is log-concave in \((x,v)\)
at fixed \(t\), and log-concave in \(t\) at fixed \((x,v)\),
since all composed magnitude arguments are affine under each
restriction. The factor \(\mu g(t\mu)\) has these same properties
by concavity of \(\log\mu\) and the negative square in its log.
Repeated one-coordinate interpolation therefore gives the
three-coordinate geometric corner lower chord on each box.

Use four uniform \(t\) panels. On each box bracket
\(\mu=x/2-rv\), require positive lower bounds, and bracket
both child magnitudes at each endpoint of a panel. The
integer replay bounds the product density and Jacobian
by the density at the upper bracket times the lower
Jacobian. Form directed four-corner means on the two
\(x\)-faces and combine them geometrically. Expand the two
guard sums in each orientation before taking the corner
means, so the two sums of means can be multiplied.
Jensen over \(t\) gives a sum of log-bilinear lower chords
on the target rectangle as in the one-crossing case below.
Add these floors to the stable degree-two target floors
for negative centers. The positive terms use only the
stable degree-two bank.
For a negative degree-one target the corresponding pair of guards
\(E_i\) and \(P_j\) contributes
\[
 \int_{d\ge\max(0,-2m)} g(d+m)E_i(v,d)P_j(v,2m+d)\,{\rm d}d
\]
for \(x,v>0\). The two child subtrees are conditionally independent.
The guard \(E_i\) starts at the positive vertex of magnitude \(d\)
and certifies the switches on its odd crossing path. Each integrand
is log-concave in \((x,v,d)\) on the linear domain by the branch
argument in the graph proof.

For this integral use rectangles in \((x,v)\) and slabs in \(d\)
from the scale-96 grid. Bisect a \(d\) slab when its width
is at least \(4/96\), using a downward geometric mean for the
interposed \(E_i\) value. Discard a box unless all its corners have
\(d\ge0\) and \(2m+d\) in the table aperture. The rational brackets
on \(r\) certify these linear conditions on the full box. In
particular \(d+m\ge0\), so the upper bracket in the decreasing
half-normal density gives a lower bound. For \(P_j\) bracket
its affine sibling argument at each corner and use the directed
row lookup.

For a fixed \(i,j\) and one slab, the product of the three factors
is bounded below by its trilinear geometric chord through
the eight downward corner values times the known Gaussian
curvature factor. At any fixed \((x,v)\), Jensen in the uniform
\(d\) coordinate gives a lower bound on its integral. It is the
slab width times the geometric mean of the values at the two ends
of \(d\), as functions of \((x,v)\), times the exponential of the
mean log-curvature in \(d\). This is a log-bilinear chord on the
outer rectangle with an explicit additional log-curvature term
in \(x,v\). For the density log \(-6(d+rv-x/2)^2+\mathrm{const}\),
the interpolation error separates into the three independent
coordinate variances. After averaging the log contribution over
the uniform coordinates, the correction is
\((\Delta_d)^2+2(\Delta_v)^2+(\Delta_x)^2/4\),
where the \(\Delta\)'s are physical cell widths.

For the final root cube, extend each such chord by ignoring the
unused neighbor coordinate. Expand products of the nonnegative
sums over slabs and guard choices, and apply Jensen to each
product over the uniform root cube. The mean logs add. Thus scalar
floors for all factors can be multiplied when they are their
means of the logarithmic bounds, exponentiated. Within a fixed
slab, the eight-corner means of the density, each \(E_i\) and
each \(P_j\) factor separately. Therefore the code multiplies
the sum of \(E_i\) means by the sum of \(P_j\) means, the density
mean and the slab width. It takes integer floors in the downward
direction and uses \(1+u\le e^u\) separately for the \(d\)
curvature and outer curvature. Summing the slab floors includes
all ordered pairs of the 74 guards without per-pair point tables.

For one center sign and two of its neighbor positions, let their
magnitudes be \(y,z\). The third magnitude is \(2rx-y-z\) when all
signs match. On the positive magnitude region the joint Gaussian
density for \((x,y,z)\) is
\[
 \frac{\sqrt{108}}{(2\pi)^{3/2}}\,
 \exp[-(x^2+(3y-2rx)^2+3(y+2z-2rx)^2)/2].
\]
This density is symmetric in the two neighbor magnitudes. Its exponent
numerator expands to \(33x^2+12(y^2+yz+z^2)-24rx(y+z)\), so use the
lower bound on \(r\) to get an upper exponent. Integrate symmetric
pair terms over one triangular half and give off-diagonal cells
multiplicity two. Directed shell terms use their actual selected order.
All corner floors retain the support \(2rx>y+z\), \(y,z>0\).
The two target subtrees are independent conditional on the star. The
selected neighbor magnitudes are integrated up to \(9/13\). The
center magnitude is integrated up to \(25/26\) for degree-one target
factors and \(29/26\) for degree-two target factors. Outside its
integration aperture each factor contributes zero.
Each is integrated by the directed corner mean. At a fixed cell
its two-variable mean is the eight-corner mean with repeated
coordinates. For a degree-one target add the individual means,
preserving the disjoint guard choices. Multiply the Gaussian
cell mean floor and the two target mean floors, then the cell volume.
For the Gaussian root density, the mean log-curvature between it
and its eight-corner chord is
\([33\Delta_x^2+12(\Delta_y^2+\Delta_z^2)]/12\).
This follows from the squared affine terms and the three uniform
coordinate variances, and uses a downward \(1+u\) factor.
Both assignments of the degrees contribute for a mixed pair.

The three unordered position pairs are counted with multiplicity
when all three neighbors qualify, so multiply each one-pair
integral by three. For positive degree-two pairs use just the
stable even-path bank, and for negative ones use that bank and
the mixed-arm extension. The other degree combinations use only
negative centers. This gives lower expectations
\(\tau_0,\tau_1,\tau_2\) of the three pair counts.

For the additional negative-center event, use the same grid and
density, with the target at position \(y\). Enforce \(y\le z\)
and \(y\le2rx-y-z\) at the corners, using the lower bound on
\(2rx\). These linear conditions make \(y\) the smallest
of the three magnitudes. The second and third magnitudes \(z\)
and \(2rx-y-z\) each use a covered-neighbor guard. Use the
74 point histories \(P_j=S_1E_j\) and the 25 stable \(D(B_i,B_j)\)
for \(i,j\in\{0,2,4,6,8\}\). These are disjoint by path histories and
crossing degree. For the histories with indices
\(0,1,2,3,4,5,6,7,9,10,15,16,44,45,46,47\), evaluate throughout
the center aperture. For the other point histories use the
available smaller parent rows and zero beyond. The corner means
for the second guard repeat the unused coordinate. For the third guard the replay brackets
its affine magnitude and forms the eight-corner mean before
adding the separate guards. Multiply by the density and the
two target-degree alternatives at \(y\), and sum over the cells.
Ties have probability zero. The factor of three counts the
distinct choices of the smallest neighbor, and only the negative
center sign is used. The resulting lower probability is \(\tau_3\).

For the additional-target count replace the second covered-arm
factor by the sum of the two target-degree factors. Keep the
third covered-arm guard and the minimum test for \(y\).
Count both choices of the distinguished second neighbor by
using a factor of six in place of three. This counts each
additional eligible neighbor at most once after the minimum
position is fixed. A third neighbor that is also eligible can
contribute in the other ordered role. This gives a lower
expectation \(\tau_4\) on the additional-target count independently
of the three-covered floor \(\tau_3\); both use events that satisfy
the combinatorial three-covered adjustment eligibility.

The resulting expectation debit is at least
\(J=\sum_{k=0}^2\beta_k\tau_k+\beta_0\tau_3+\alpha\tau_4+
\beta_0\sigma_0+\beta_1\sigma_1\),
with exact rational value as checked by the replay. The finite
witnesses discard centers near crossing-degree-three defects as
in the graph proof, with the exclusion radius enlarged to contain
the new guard histories.

### Coupling finite neighborhoods with the limiting wave

Fix a radius-\(D\) tree ball. Its limiting covariance has kernel
spanned by the local eigenvalue relations
\(\sum_{u\sim w}X_u-2\sqrt2 X_w=0\) at the internal vertices \(w\).
To verify this, generate the limiting scores from the root score,
the root's two-dimensional neighbor residual, and one independent
forward contrast at each other internal vertex. The root residual
has conditional covariance \(1/9\) on the diagonal and \(-1/18\)
off the diagonal, with rank two. The other contrasts have variance
\(1/12\) and are independent by the branch calculation. The input
innovations are linear functions of the scores. If there are \(b\)
internal vertices, the rank is \(1+2+(b-1)=b+2\), which is the count
of leaves. The \(b\) displayed relations are independent. In any
nonzero combination, a deepest internal vertex with nonzero
coefficient contributes that coefficient at one of its children,
where no other nonzero relation can contribute.

Let \(V\) be the orthogonal complement of that kernel. The covariance
\(A_V\) of the limiting wave on \(V\) is positive definite on this
fixed finite-dimensional space. Write \(X_R\) for the approximate
Gaussian vector on the ball, \(Y_R\) for its projection onto \(V\),
and \(B_V\) for the covariance of \(Y_R\). The entrywise covariance
estimate gives \(\|A_V-B_V\|_{\rm op}=O_D(R^{-2})\).
The eigenvalue residual variance bound gives
\(\mathbb E\|X_R-Y_R\|^2=O_D(R^{-4})\). Indeed the fixed kernel
basis consists of the internal star residual vectors, each with
variance \(O(R^{-4})\), and any orthonormal kernel basis is made
of fixed linear combinations of them.

For large \(R\) both covariances on \(V\) have a fixed positive
lower eigenvalue bound. Their square roots then differ by
\(O_D(R^{-2})\) in operator norm. For completeness, if
\(S=\sqrt{A_V}\) and \(T=\sqrt{B_V}\), solve
\(S(S-T)+(S-T)T=A_V-B_V\) by integrating
\(e^{-uS}(A_V-B_V)e^{-uT}\) over \(u\ge0\). The positive
lower eigenvalue bounds give the asserted norm estimate.
With \(Z=B_V^{-1/2}Y_R\), the vector \(X=A_V^{1/2}Z\)
has the limiting law. It is coupled with \(X_R\) so that
\(\mathbb E\|X-X_R\|^2=O_D(R^{-4})\), including the separate
orthogonal kernel error.

The combinatorial witnesses use finite signs, an ordering of two
child magnitudes within fixed signs, and bounded sign paths and
stars. Their numerical
integrals use inward apertures within these events. Each
combinatorial decision is a Boolean combination of finitely many
centered linear sign tests involving one score or two distinct
scores. These forms have positive limiting variance, because
\(|r_j|<1\) at every positive distance. For each form, standardize
its limiting and approximate values. Their mean-square difference
is still \(O_D(R^{-4})\), by the triangle inequality and the
fixed positive variance bound. If their correlation is \(\rho\),
this squared difference is \(2(1-\rho)\). Their sign
disagreement probability is \(\arccos(\rho)/\pi=O_D(R^{-2})\)
by the two-dimensional Gaussian angle formula. A finite union
bound gives the same error for every local event.

This controls the bounded alternating path counts, edge sign
counts, base witnesses and isolated-center witnesses. A star
with three opposite neighbors has zero limiting probability,
so the expected fraction of defects and of excluded neighborhoods
is \(O_D(R^{-2})\). The eligible-target count, the additional
covered-arm predicates and the smallest-target decision are all
controlled by those finite sign tests. The simultaneous
adjustments apply equally to the approximate wave.
At zero threshold the expected cover charge per corresponding
set of \(v\) anchors is at most \((a_0+O_D(R^{-2}))v\).
Submodular comparison and the threshold-maximum argument give
this bound at every physical threshold.

### Exact batched vertical scaling

For each old point-history row let \(S=10^8\), and let \(v\) be
one directed branch probability ordinate and \(d\) the directed
half-normal density. Here \(0\le v\le S\) and \(0\le d\le3S\).
With \(R=\lceil 2^{96}/S\rceil\),
\(0\le dvR/2^{96}-dv/S<3S^2/2^{96}<1/S\). Thus
\(\lfloor dvR/2^{96}\rfloor=\lfloor dv/S\rfloor\), since each
noninteger quotient has distance at least \(1/S\) to the next
integer. The product before shifting fits within 128 bits. In the
packed row each lane starts at a multiple of 128 bits. Taking the
low 32 bits of each lane after the right shift by 96 gives precisely
the independent directed quotient. Lane packing and reusing the
first denominator division across equal root tables preserve each
integer floor.
