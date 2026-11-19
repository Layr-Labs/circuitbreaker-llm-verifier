# A joint medium reserve for a marked Gaussian ordering

We prove a limiting cutwidth coefficient strictly below
`gamma=500000000/3827670601` for every simple cubic graph.
Section 3 of `family.md` gives one fixed total uniform NP language
with `C_B2(L_N)>=(4827670601/1000000000)N` at every sufficiently large
actual integer length N, with zero subtracted remainder.

Covariance comparison pays for dirty windows. Guarded ordering adds
savings near zero on the same field.

## 1. Crossing matchings and Gaussian comparison

For vertex bits `b` on a cubic graph, write `F(b)` for the minimum cover
of edges joining opposite bits. It equals their maximum matching size
by bipartiteness. For edge bits `x` put

```
Q(b,x)=sum_v NAE(b_v,(x_e : e incident to v)),
F(b)=min_x Q(b,x).
```

NAE indicates both bits occur. Nonzero-cost vertices cover all crossing
edges; conversely a crossing cover lets outside vertices consistently
fix their incident edge bits. NAE has nonpositive mixed differences.
Taking meet/join of joint minimizers proves F submodular. F is
complement-invariant and changing one bit costs at most one.

For centered unit Gaussian fields `X,Y` with covariance matrices `A,C`,
put `Z_s=sqrt(1-s)X+sqrt(s)Y` with the fields independent. At threshold
zero, Gaussian integration by parts gives the Price identity

```
(d/ds) E F(1[Z_s<0]) = sum_(i<j) (C_ij-A_ij) phi_(ij,s)(0,0)
                     * E[Delta_ij F | Z_(s,i)=Z_(s,j)=0].            (1)
```

Here phi is the bivariate density and Delta the four-corner difference
with other coordinates thresholded. Smooth the multilinear extension's
thresholds: the Gaussian Hessian identity pairs the off-diagonal terms;
the two boundary derivatives give (1), and diagonal coefficients vanish.
Positive-definite approximations extend by marginal continuity and
bounded F. The proof applies at any common threshold. Entrywise covariance
domination therefore decreases the expectation.

For a single field of covariance `C` let `f(t)=E F(1[X<t])` and
`W_j=X_j-C_ij X_i` under conditioning at `i`. Set `k_j=1-C_ij` for
that `i` and put `g_i(t)=E[D_i F(1[X_-i<t]) | X_i=t]`. The density
of the conditioned boundary gives `f'(t)=phi(t) sum_i g_i(t)`.
Symmetry gives `g_i(-t)=-g_i(t)`. On interpolating from `-t` to `t`
the conditioned tests are `W_j<k_j s`, so differentiating their
boundary integrals gives

```
2g_i(t)=integral_(-t)^t sum_(j!=i) k_j f_(W_j)(k_j s)
           E[Delta_ij F | X_i=X_j=s, other bits tested at s] ds.     (2)
```

Condition first on X_i and then on W_j=k_j*s; f_(W_j) is its density.
Thus f is even and decreases for t>=0. Independent coordinate noise
and marginal limits include singular fields.

## 2. From nested cuts to a layout

A vertex ordering with every prefix crossing cover <=M yields width
<=M+6 by the following edge interpolation.
[Cutwidth proof, Lemmas 2.2 and 3.1–3.2](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/737818f9911e00e9d4909c644c771e0724fb6c37/submission/cutwidth.md#L279).

Clamp same-anchor edge bits. Submodularity gives monotone least
minimizers along prefixes. Flipping one anchor costs at most M+1.
In the monotone interpolation let A gain and B lose nonzero list cost,
P retain it. Resetting any subset of A in the destination optimum
gives Hall's inequality: its changing B-neighborhood is at least its
size. Process P–P,B–P edges and three-edge A stars first. Remaining
A blocks have at most two B neighbors, hence form edges or half-edges
on B. Hall implies tree components with at most one half-edge, or
unicyclic components with attached trees. Peel leaves, then finish the
half-edge or cycle. At most one block is unfinished between blocks,
two within a block, so split lists cost <=M+3. The edge-position bag
adds at most two endpoints. In median-incident-edge vertex order,
a crossing edge before the median charges the outside vertex's lowest
edge, one after it the inside vertex's highest; both enter the bag.
At most one edge at the median remains. Thus width <=M+6.

## 3. Radial keys on the physical graph

Let the base be simple connected cubic. Assign free generators outside
a spanning tree. A forbidden nonbacktracking closed walk gives a nonempty
reduced word. For length j, prescribe partial permutations along points
0,...,j; reducedness makes these injections, extendible to permutations.
The product of the finitely many separating groups, acting regularly on
fibers, gives any required finite girth. On a degree-d cover repeated
bits have objective dF: lifting a base cover gives <=; averaging an
upstairs cover gives a fractional base cover whose inequalities,
summed over a maximum matching, give >=.

Use the nonnegative unit row `f_v(z)=sqrt(c_a) a^dist(v,z)` on the
universal cubic tree, with `7/10<=a<1/sqrt(2)` and
`c_a=(1-2a^2)/(1+a^2)`. The sphere sizes `3*2^(r-1)` give the norm.
At distance `j`, grouping by the nearest point on the joining path gives
the covariance and its critical limit

```
r_j(a)=a^j [(j+1)-(j-1)a^2]/(1+a^2),
w_j=(1+j/3)2^(-j/2).
```

Take fiber norms to obtain unit physical rows b_v. Fiberwise
Cauchy–Schwarz dominates the independent-cover covariance. Define
`D_e=<b_u,b_v>-r_1(a)>=0` from complete rows and keep these marking
numbers fixed. For the actual field truncate at R, renormalize before
fiber norms, scale by sqrt(1-epsilon^2), and add independent physical
noise epsilon. Geometric squared-mass tails and contraction under fiber
norms give uniform covariance convergence. Choose a, then R,epsilon,
then cover girth for every window. Apply the same truncation/noise
upstairs, repeating physical noise along fibers; domination includes
repeated images.

Corrected keys inherit original keys within distance three, each used
at most 22 times. For fixed R each independent input affects O_R(1)
corrected keys. Bounded differences gives variance O_R(n) for a fixed
threshold's cover and for original-key interval counts. Chebyshev controls
a fixed mesh simultaneously. Between nodes, changes cost <=22 times
the interval count, of mean <=n*length/sqrt(2pi). Tail counts handle
unbounded thresholds. Choose mesh/tail errors before n. Independent
noise separates original keys almost surely; inherited ties have size
<=22, so partial tie blocks add <=22. Section 2 now turns uniform
expected bounds into limiting cutwidth, with all source/mesh choices
independent of the graph and its order.

## 4. Reference expectation at zero

For odd J let `N_j` be half the directed nonbacktracking j-walk count
of the crossing graph. The sum `sum_(j<=J)(-1)^(j+1)N_j` bounds
the cover on paths, with values `ceil(r/2)` for r<J edges and
`r-(J-1)/2` otherwise, and on cycles, where it equals their length.
Deleting b degree-three vertices changes that bound by `O_J(b)`,
since each deletion affects at most `(j+1)3*2^(j-1)` directed walks.
At the critical covariance, the zero-variance star relation
`sum_(u~v) X_u=2sqrt(2)X_v` makes a three-edge crossing star
impossible. Its probability tends uniformly to zero on the
high-girth reference. There are `3*2^(j-2)` unoriented j-paths per
vertex. The resulting limiting expected cover per vertex is at most
`Gamma_J=sum_(j<=J)(-1)^(j+1)3*2^(j-2)p_j`, with p_j the
probability of alternating signs on that tree path.

We use the same quadrature argument as the published path certificate,
with `J=25,N=60,L=3,S=2^224` and an independent mandatory replay.
[Path certificate, Sections 2–8](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/737818f9911e00e9d4909c644c771e0724fb6c37/submission/certificate.md).

The AR(2) recurrence has coefficients `sqrt(2),-1/2` and innovation
variance `1/12`. With alternating signs, its positive pair
density is `c2 exp[-9(x^2+y^2)/2-6sqrt(2)xy]` and extension
density is `c1 exp[-6(z+sqrt(2)y+x/2)^2]`, where
`c2=3/(2pi),c1=sqrt(6/pi)`. Reversal permits any anchor pair.
For `m=j+1>=2`, an outlier at i above L is anchored at
`h=min(i,m-2)`. Integrating the other anchor coordinate and
bounding each one-step/two-step completion by `1/2,1/10`
gives a single-orthant tail at most
`c2 exp(-9L^2/2)/(54sqrt(2)L^2) * H(h)H(m-2-h)`,
with `H(r)=10^-floor(r/2)2^-(r mod 2)`. Sum over i.

The 2Nth derivative bound in coordinate i is
`f_m(0)q_i^N(2N)!/(2^N N!)` by the Gaussian Fourier moment,
with diagonal precision q_i. The Gauss monic-square remainder
coefficient is `L^(2N+1)(N!)^4/[(2N+1)((2N)!)^3]`.
Telescoping the tensor approximation multiplies by `L^(m-1)`
and sums over i. For each term divide by `(4L)^max(0,m-9)`
using the unchanged interior factors `exp(-39x^2/2)`, whose
integral and quadrature mass are at most `1/4`.
The checker isolates Legendre roots in disjoint sign brackets,
encloses weights via `(x^2-1)P_N'=N(xP_N-P_(N-1))` and
`||P_N''||<=N^4`, and rounds radicals and exponentials outward.

Factor the positive weighted transition as `A_ij B_ik C_jk`
by separating the `xy,xz,yz` products. The checker verifies
`A,B<=1` and `sum_k C_jk<=1`. If the coefficient enclosure
widths are `D_A,D_B,D_C` at scale S, the per-step additional
lost mass is at most `(D_A+D_B+ND_C+2N^2)/S`, including two
sets of downward state roundings. Accumulate it from the anchor
tensor using contraction. Doubling orthant bounds and applying
the correct endpoint to each alternating term verifies
`Gamma_25<U=0.13064756648`.

## 5. Graph margin

Put `T=9/400`, `U=.13064756648` and
`C=1981610721/10^14`. Sections 8–11 prove a reduction C
per vertex for `|t|<=T` and a larger reduction at every outer
threshold. The exact strict margin is

```
gamma-(U-C)=22330492665321/382767060100000000000000 > 0.
```

Reserve part of this slack for the finite profile, truncation, jitter,
concentration and layout errors before choosing graph size. Concatenate
component layouts; components below the common size threshold have
bounded width. This proves the simple cubic hypothesis in `family.md`.

## 6. Safe advances and strict matching savings

Apply Gaussian comparison to the original field, then account for
advances at each threshold. Let B be its below-threshold set. All
named vertices below are physically distinct and listed stars complete;
lightness proves this faithfulness when counting events.

For the seven-key advance, name the root `o` with neighbors `a,b,u`.
Use complete stars at `a`, with further neighbors `a1,a2`, at `u`,
with further neighbors `w,w'`, and at `w`, with further neighbors
`w1,w2`. Put `a,b,a1,a2,w,w1,w2` below the threshold and `o,u,w'`
above. Delete `o` from the crossing graph. Then `a` is isolated, so
restoring the original `o` adds one to the maximum matching. After
advancing `o` into `B`, the only new crossing edge is `ou`. It can
be replaced by `uw` in any matching using it, since `w` has no
other crossing neighbor. The maximum matching has decreased by one.

For the nine-key advance take stars on the tree with edges

```
oa, ob, ou, aa1, aa2, bw, bw', uz, uz', ww1, ww2, zz1, zz2.
```

Put `a,a1,a2,w,w1,w2,z,z1,z2` below threshold and `o,b,u,w',z'`
above. After deleting `o`, the external crossing graph has the
isolated vertex `a` and closed pairs `b–w` and `u–z`. Originally
`oa` adds one to a maximum matching. After the advance the new
edges `ob,ou` can be covered by choosing `b,u` in minimum covers
of their closed pairs. Thus this advance also saves one.

The positive support alone is strict: starting from the guarded
configuration, turn other outside named bits inside. Submodularity
keeps the root marginal <=-1, its minimum possible value.
Guards only distinguish counted events. A seven-key guarded
success has two inside neighbors of its root, and determines `u,w`.
Among twelve names at the root, at most two can succeed, weighted
`1/2` each. In a nine-key success the inside neighbor is unique and
each of the other two neighbors has a unique inside child. Its
twelve guarded names are disjoint. Its successes are also disjoint
from the seven-key guarded successes at the same root.

A weak support consists of `a,a1,a2` where `a` neighbors the root and
`a1,a2` are its other neighbors. When below threshold, `a` is isolated
after deleting the root. Restoring the original root adds one, and
restoring it inside adds at most one. To count the weak-pair event,
choose a center `m` outside whose neighbors are all outside. Qualify
an outside neighbor `r` with both other neighbors `a,b` inside,
where the lower-key neighbor `a` is a private leaf. Break key ties by a
fixed vertex order. The weak proposal
uses `H=max(Z_a,Z_a1,Z_a2,Z_b)` if `Z_m+(3/2)H >= (5/2) Z_r`, with
all keys taken before the pass. This filter preserves the weak
marginal. Qualify only roots whose proposal passes the filter.

Each such root has unique original outside neighbor `m`, so their
moving sets `R_m` are disjoint. The private leaves are also distinct.
Originally `m` is isolated, and removing its `k` qualifying private
pairs reduces the matching by `k`. After the moves the roots are
leaves at `m`, adding one when `k>0`. The group improves by
`(k-1)_+` and each individual marginal is zero. After previous
disjoint group moves the same group marginal can only decrease, by
telescoping its singleton marginals and using submodularity. A
strict root is outside these groups by original neutrality. Twelve
canonical weak-pair names at a center count `binomial(k,2)` successes
when all qualifying names are allowed, or fewer under lightness.
Weight each by `2/3`, so their total is bounded by the group saving.

Lower each root to the minimum of its original key and admitted
support maxima, including filtered weak proposals. Graph/original-key
masks are fixed across thresholds. A moved root has new_key<t<=old_key,
so its certificate applies. Process weak groups, strict roots, then
other weak moves; submodularity preserves savings and safety.
Minima/maxima select original keys within radius three.

## 7. Radial excess and conditional defects

Small excess ensures faithful patches; neighboring excesses pay marking.
Use complete rows before truncation/jitter, keeping all marks fixed
throughout comparison.

Write `P_v(p)=b_v(p)^2` for the squared radial mass on fibre `p` of
the row at `v`. With `q=a^2` and `B` the base adjacency matrix, its
exact matrix is

```
P=(1-3beta)(I-beta B)^(-1), beta=q/(1+2q^2).
```

Indeed the counts of nonbacktracking walks start with `I,B,B^2-3I`
and obey `N_j=B N_(j-1)-2N_(j-2)` for `j>=3`. The normalized power
series `c_a sum_j q^j N_j` gives the formula. Hence for each `v,p`,

```
sum_(u~v) P_u(p)=beta^(-1) P_v(p)-(beta^(-1)-3)1[p=v].             (3)
```

The ratio `h=P_u(p)/P_v(p)` on an edge lies in `[q,1/q]`. On the
universal tree each point lies on one side of the edge and its two
squared row entries have ratio `q` or `1/q`. The fibre ratio is the
weighted mean of these ratios. Set

```
phi_q(h)=sqrt(h)-sqrt(q)(1+h)/(1+q),
d_(v,u)(p)=P_v(p) phi_q(h),   D_(v,u)=sum_p d_(v,u)(p).
```

The summands are nonnegative and symmetric in `v,u`, and the
last equality uses `r_1(a)=2a/(1+a^2)`. The ratios and these
summands include repeated physical projections.

On a fibre the excess at neighboring lifts is half the squared
difference of the unit row directions times the product of their
fibre norms. Along a tree path of length `j`, telescope the directions
and apply Cauchy–Schwarz. The quotient of the two end norms' product by
any successive pair's product is at most `a^(-(j-1))`. Thus, whenever
every edge of the lifted path has excess at most `H`,

```
w_j-o(1) <= Cov(Y_v,Y_u) <= w_j+j^2 (sqrt(2))^(j-1) H+o(1).          (4)
```

Here and below `o(1)` denotes a uniformly small source allowance
chosen before graph size. For the ordering patches take `H=eta` or
`eta/2` as prescribed. The relaxed nine-key patch allows `eta` at its
four terminal child edges and uses `eta/2` elsewhere. The maximum
distance is six and `w_j+j^2(sqrt(2))^(j-1) eta<1` for `1<=j<=6` at
the fixed `eta=323/400000`. Thus no two distinct patch names project
to the same physical vertex. The same argument with the actual
per-path upper bounds can be used when a cleaner inner window is
needed within a larger Price event.

Let `Lambda=2sqrt(2)`, `c=3-Lambda` and `R_v=sum_(u~v)Y_u-Lambda Y_v`.
For the complete rows, applying the triangle inequality to the sum of
the three lifted neighboring vectors gives the coordinatewise lower
bound `(1/a+2a)b_v-e_v` before subtracting `Lambda b_v`. The correction
`e_v` is supported at the image of the root and has norm at most
`(1/a-a)sqrt(c_a)`. For the upper bound apply Cauchy–Schwarz across the
three fibre norms and (3), giving
`sum_(u~v)b_u <= sqrt(3(a^(-2)+2a^2)) b_v`. The row of `R_v` therefore
differs by norm `o(1)` from a nonnegative row `r` with `r<=(c+o(1))b_v`.
Put `S=<r,b_v>` so `0<=S<=c+o(1)` and `S=sum_(u~v)D_(v,u)+o(1)`.
For any fixed lifted distance `j` to `x`, the row ratios give

```
||r||^2 <= cS+o(1),
(sqrt(2))^(-j) S-o(1) <= Cov(R_v,Y_x)
                         <= (sqrt(2))^j S+o(1).                    (5)
```

For a clean star with bound `H` on its three edges, use `S<=3H+o(1)`.
For a pivotal score at zero in the independent interpolation of
Section 1, conditioning subtracts at least `s 2^(-j)S^2-o(1)` from
the actual defect variance at distance `j` from that endpoint.
Maximizing over `S<=c` gives `c^2 h_j(s)` where `h_0=1-s` for
`s<=1/2` and `1/(4s)` thereafter, and `h_j=1-s/2^j` for `j>=1`.
For the clean caps of the Price calculation the parabola is
increasing up to `3H`, yielding `3cH-9sH^2/2^j`. For the stratified
first outward star replace `3H` by `.006+2*.007=1/50`. Additional
pivots and a centered positive quadratic tilt decrease the covariance.
The resulting standard deviations can be used in centered strips.
In a branch at depth `d` the propagated midpoint error has
coefficients `(i+1)2^(-i/2)/2`, `0<=i<d`, against defects at depths
`d-1-i`. They are the impulse response of the repeated root
`1/sqrt(2)` in the tree midpoint recurrence. Apply the triangle
inequality for conditional standard deviations to that sum.

For the energy allocation we need a paired inequality at an edge.
For one fibre at one endpoint put `x` equal to the fraction of mass
in the branch towards the other endpoint. With `b=1/a` set
`f(x)=sqrt(a^2+(b^2-a^2)x)-a-(b-a)x`. The excess integrand is
`P_v(p)f(x)`. The two away fractions sum to `1-x` minus the root
atom's fraction. Concavity bounds their total integrand from below by
`P_v(p) f(1-x)` minus a uniform Lipschitz bound times the root atom.
At the other endpoint the mass ratio is `h=a^2+(b^2-a^2)x` and the
reverse fraction is `x'=a^2(1-x)/h`, including the first endpoint's
atom. For `y=sqrt(h)+1/sqrt(h)` and `z=a^2+b^2`, the sum of the two
lower integrands divided by `P_v(p)sqrt(h)` is `R(y)-zy/(a+b)` with

```
R(y)=sqrt(z^2-(u-1)^2),   u=sqrt((z+1)^2-zy^2).
```

Here `2<=y<=a+b` and `1<=u<=z-1`. The second derivative of
`(u-1)^2` is `2z((z+1)^2-u^3)/u^3>0`, so `R` is concave. Its
values at the ends are `2sqrt(z-1)` and `z`. The chord and
`f(x)/sqrt(h)=1-y/(a+b)` show that the four off-edge excesses are
at least `C_a D_e-o(1)`, where

```
C_a=(2(a+b)sqrt(z-1)-2z)/(a+b-2).
```

The two root-atom losses are bounded by a fixed constant times
`c_a`, uniformly over the base. The limiting `C_a` is
`9sqrt(6)+12sqrt(3)-15sqrt(2)-20 > 2021/1250`. The inequality
follows using the lower square-root bounds `2.44948974,1.73205080`
and upper bound `1.41421357` for the respective radicals. We may
use `sigma=2021/1250` with arbitrarily small source debit.

Finally, at a vertex the weighted triangle inequality is
`phi(t)<=r phi(r)+s phi(s)` for `t,r,s` in `[1/2,2]` of sum three,
where `phi=phi_(1/2)`. Here is a scalar verification. Fix `t`; an
interior stationary point of the right side has `r=s` or
`sqrt(r)+sqrt(s)=9/(4sqrt(2))`. At the boundary set `r=1/2`,
`s=5/2-t`. In the rationalized `phi`, the numerators at `s,t`
agree. If `D(x)=sqrt(x)+sqrt(2)(1+x)/3`, then
`sD(t)-D(s)=sqrt(s)(sqrt(st)-1)+sqrt(2)(st-1)/3>=0`.
At `r=s` a minimum requires `r<=81/128` by the second derivative
of `x phi(x)`. Put `r=x^2/2`, `u=x-1` so `0<=u<=1/8`.
The difference after squaring the positives
`8+3x^3-4x^2-x^4` and `3sqrt(6-2x^2)` is
`u^2(15-6u-5u^2+8u^3+3u^4+2u^5+u^6)>=0`.
For the other stationary points put `L=r+s`, in `[81/64,41/32]`.
The difference multiplied by `sqrt(2)` is `Q_0(L)-sqrt(6-2L)`, with
`Q_0(L)=-729/256+27L/8+(2/3)(4-2L-L^2+(81/32-L)^2/2)`.
Here `Q_0'(L)=17/48-2L/3>=-1/2`, so the difference is increasing
and its value at `81/64` is the verified balanced case. Continuity
covers the boundary zeros.

By (3), the expectation with weights `P_v(p)` of the absolute
error in the sum of the three fibre ratios from three is at most
`2(beta^(-1)-3)`. Clamp the ratios to `[1/2,2]` and adjust them in
that box to sum three. Their weighted total displacement is `o(1)`.
The bounded derivatives of `phi_q` and `h phi_q(h)` on the common
compact interval show that `D_i<=sum_(j!=i) sum_p h_j d_j(p)+o(1)`.
In particular `D_i<=2 sum_(j!=i)D_j+o(1)`. When one low edge is
charged from both ends, the two weights of its symmetric summand
are `h` and `1/h` and their sum is at most `q+1/q`. This improves
the budget available from edges between unmarked vertices.

## 8. The paired weak move on the actual graph

Use the weak group from Section 6. Write `m` for its center and
`r,s,h` for its neighbors. At the two moving roots use child pairs
`(a,a'),(b,b')`, with ending pairs at `a,b`. This inspects the
eleven edges of the five stars at `m,r,s,a,b`. All must have
`D_e<=theta=eta/2`. The path bound (4), with their diameter six,
gives distinct physical images. For each star let
`D_v = sum_(u~v) Y_u-lambda Y_v`, here a Gaussian variable.
By (5) its variance is at most `3(3-lambda)theta+o(1)`, and its
covariance with a named key at distance j lies between
`-o(1)` and `3 theta sqrt(2)^j+o(1)`.

Use the driver vector

```text
O=Y_m                       R=2lambda Y_m-3(Y_r+Y_s)
H=sqrt(3)(Y_r-Y_s)           G_r=sqrt(3)(Y_a-Y_a')
G_s=sqrt(3)(Y_b-Y_b')        J_a=sqrt(3)(Y_a1-Y_a2)
J_b=sqrt(3)(Y_b1-Y_b2).
```

On the ideal tree the covariance is the identity. The two
center contrasts and the later sibling contrasts are orthogonal
by the sibling permutations and have variance one using
`w_2=5/6`. On the actual graph use `K=I+E`. For each pair of
driver rows with coefficients u,v, sort the products
`u_i v_j` by sign. The pairwise covariance errors are nonnegative.
For a sharper bound, use the weighted Cauchy–Schwarz inequality
on the telescope of unit directions in one fiber. Its coefficient
is the sum of `b_0 b_j/(b_(i-1)b_i)` over the j edges, where b
denotes the fiber row norms. Let `r=1/a`; the successive differences
of `log_r b_i` lie in [-1,1]. The coefficient is convex in those
differences, so its maximum is attained on the endpoints. Swapping
an up-step followed by a down-step increases both affected ratios.
Thus it suffices to use a run of down-steps followed by up-steps.
The coefficient is `r^(j-1)` times the sum of the two geometric
series in `r^(-2)` with their respective run lengths. It is largest
for lengths differing by at most one. At the critical parameter
this gives the upper error
`theta j sqrt(2)^(j-1)(4-2^(1-floor(j/2))-2^(1-ceil(j/2)))`.
Sum the positive products and separately the negatives to bound both
ends of the error interval. The larger absolute value gives
`P_ij`. The replay checks `max_i sum_j P_ij<1`. The geometric
matrix series gives `|K^(-1)|<=(I-P)^(-1)=W` entrywise.
Hadamard's determinant bound on K gives the checked density prefactor
`1/636`. The conditional defect means similarly use the absolute
majorants for `Cov(D_v,drivers)` multiplied by W, and projection
can only decrease each variance.

Use the coordinates `Y_m=t+z,Y_r=t+xz,Y_s=t+yz`, and divide the
four later contrasts by z. The Jacobian is `6sqrt(3) z^6`.
Use `.22<=z<=.8`, six intervals of width `.03` on `x,y` and six
of width `.15` on the negative first-contrast magnitudes,
starting at `.05`. For the two ending contrasts use the checked
radii up to `.9`. Center the five conditional defects and restrict
every one to `.075z`. Their joint probability is at least the
product of the five strip floors by Gaussian correlation for
centered symmetric convex sets.

For one moving root write `x` also for the upper box endpoint,
`u=(3-lambda_lo)T/.22`, and `mu,nu` for the defect mean
majorants divided by z at its star and selected child's star.
The children have midpoint increment at most
`M z=(lambda_hi*x-1+u+mu+.075)z/2`. The larger child's sign
uses `M+g_hi/(2sqrt(3)_lo)<0`, for first contrast magnitude in
`[g_lo,g_hi]` after division by z. The lower child is selected
and has increment at most `A z=(M-g_lo/(2sqrt(3)_hi))z`.
Its ending contrast cap is
`sqrt(3)_lo(-lambda_lo*A-u-nu-.075)`; round strictly below it.
The affine filter follows from the lower pair midpoint, which
is no larger than the maximum H in that filter. Its sufficient
inequality after division by z is
`1-(3/4)(1+mu+.075)+(3lambda_lo/4-5/2)x>3u/4`.
The remaining center sign is checked from the star equation at m.
These tests certify every admitted box with the actual conditional
means and with no hypothesis on residual independence.

The check uses the seven driver amplitude bounds on each box
and W to upper-bound the exponent by `a z^2`. Integrate
`z^6` on 116 consecutive slabs. Take the right endpoint for
the decreasing exponential and the left endpoint for the increasing
defect strips. The mesh has disjoint interiors. The exponential
and strip routines use odd Taylor lower sums and downward scaling.
The resulting exact sum strictly exceeds `12719/10^12`. Uniform
finite-source errors fit the strict gaps. The driver covariance
remains nonsingular and threshold signs have no atoms, so the
same finite-profile choice works for all weak namings.

## 9. Reductions at the outer thresholds

Apply (2) to the ideal cover. A distance-d ordered pair contributes
its conditional mixed difference times the drift `1-w_d` and
the density of `N(0,1-w_d^2)` at `(1-w_d)s`. The cover has
`3*2^(d-1)` ordered distance-d pairs at each vertex for each
retained distance. Every omitted pair has nonpositive difference.

For endpoints a distance d>=2 apart, use alternating middle
colors. Require every middle third neighbor to match its middle
color. At each endpoint either require both outer roots opposite
the adjacent middle, or remove a same-color root whose outside
crossing component is a closed even path with all its boundary
colors specified. Its perfect matching contributes a fixed
constant. Let A and B be the two separate matching increases when
attaching the surviving exterior edges at endpoints in the
adjacent middle colors, and C the joint increase. For odd d
the additions lie on opposite sides of the crossing graph.
If A=B=1, the symmetric difference of two exterior optimum
matchings has only even alternating paths and cycles. The two
exposed vertices lie in opposite sides so cannot be the ends of
one such path. Choose components independently to expose both,
giving C=2. Thus C>=A+B. For even d the additions lie on the same
side. If C=2, deleting either added vertex shows A=B=1. If A=B=0,
an augmenting path after adding both can use only one new vertex
as an endpoint, so C=0. Consequently C<=A+B. In flip bits x,y
relative to the middle colors, the path term is
`floor((d-1+x+y)/2)`. Combining its difference with C-A-B and
converting back to physical bits gives `Delta<=-1` in both parities.

Under the critical tree conditioning `Y_0=Y_d=s`, write J for the
centered middle block. The middle increments over threshold are
`J_i+b_i s` with `b_i=(w_i+w_(d-i))/(1+w_d)-1`.
The off-path third-neighbor increments are
`lambda J_i-J_(i-1)-J_(i+1)+(lambda b_i-b_(i-1)-b_(i+1)-(3-lambda))s`.
The endpoints of J and b vanish. The two pairs of exterior siblings
have midpoint increments
`-(J_1+(b_1+3-lambda)s)/2` and its reflected counterpart.
Their contrasts are independent of J and one another by the two
subtree swaps. The checker uses upper absolute offsets on `|s|<=.03`.
For retained distances two through six the magnitudes start at
`.01` and the opposite-child strips use positive distances from
the shifted midpoint.

For d>=3, the middle precision has diagonal 36 at its ends and
39 in the interior, off-diagonal entries `-18sqrt(2)` at distance
one and 6 at distance two. This follows from the critical path
recurrence `Y_j=sqrt(2)Y_(j-1)-Y_(j-2)/2+epsilon_j`, with
innovation variance `1/12` and initial pair variance determinant
`1/9`. The conditional precision determinant is
`9*12^(d-1)*(1-w_d^2)`. Alternating signs make the exponent on the
positive magnitudes a sum of positive diagonal, adjacent and
distance-two products. The normalizing factor includes both signs.
For d=2 the exponent is `33x^2/2` on the one middle magnitude.
Integrate the two side-strip lower polynomials (a square when
they use one coordinate) and use their weighted first and second
moments to upper-bound the mean exponent on each box. Jensen
therefore gives a lower bound for the box integral. The exact
contractions use the joint middle density.

For extra arms, put the same-colored root and its outward path in
alternation. On a closed path of size two its other child and the
ending children have the required noncrossing signs. For size four
use two more alternating vertices and the corresponding complete
stars. The conditional successive branch contrasts are independent
by subtree swaps. For parent amplitude a and previous amplitude p,
the continuation density at b is bounded below by the directed
normal density at `b+sqrt(2)*a+p/2+.0026`. The terminal strip
uses `2sqrt(2)*b+a-.0052`. The two root half-differences use
the same bound on the endpoint midpoint shift, paid by `.008`
in their root bounds. The amplitude cells have width `.01`.
Each terminal strip uses its lower endpoint.
An arm on either side or on both gives disjoint sign cases, separated
by the root sign and the number of vertices in the full exterior
crossing component. At distance two the one-arm contribution is paid once, separately
from the baseline.

Keep the terminal strip at its lower corner and the fixed middle
magnitude in the arm root factor at its upper corner. The remaining
product of density factors has exponent `-sum Q_j`. Under the
uniform measure on the box, Jensen gives
`E exp(-sum Q_j) >= exp(-E sum Q_j) = prod exp(-E Q_j)`.
The artificial product measure is used only to compute the cell
moments. The factors may share arm amplitudes.

For upper cell index i the uniform mean is `(2i-1)/200` and
variance is `1/120000`. The mean of each squared affine form is
its squared mean plus the sum of coefficient-squared variances.
The checker uses these means and rounds the positive cross terms
with the upper radical. For a continuation with preceding, current
and following upper indices p,a,b, take the exponent index as the
ceiling of `(3M^2+V)/50000000000`, where
`M=26000+25000(2p-1)+70711(2a-1)+50000(2b-1)` and
`V=25000^2+70711^2+50000^2`. Omit the preceding terms at the
first transition. The root index is the ceiling of
`(3(10x+20a-2)^2+100)/2000` with upper middle index x.
This gives an integral lower bound for every arm box. Contract
positive factors with downward integer rounding. At D2 and D3
also average the middle density, with fixed arm weights on each
cell. The longer middle contractions already use uniform moments.
For D4–D6 the end mass is `A_1(2B_2+A_2)` where A is the
sum of size-two and size-four arms and B is the no-arm strip.
Reflection accounts for an arm only at the other endpoint.

Retain no-arm D2–D4 floors `.0482,.00419,.00089`. Family §13 checks D4. The three arm
floors (single size two, double size two, at least one size four)
are `.0552889,.0313397,.0142689` for D2 and
`.0114197,.0075179,.0032000` for D3. The D4–D6 positive-arm floors are
`.004439266,.000874902,.000166833`, with strict replay slack.
At distances 2–6 the drift floors are `.165,.292,.416,.528,.624`
and the residual-density floors are `.71,.56,.49,.45,.429`. For
`v=1-w_d^2` bound the density by `(1-(.03)^2/(2v))/sqrt(2pi v)`
from below. D4–D6 use `v=95/144,7/9,55/64`. Every retained distance
and arm class contributes once; all omitted terms pay zero.

The middle blocks are nonsingular and all strict sign boundaries
have probability zero. Thus the finitely many strict floors
hold uniformly for sufficiently accurate finite source choices.
The transfer does not depend on graph size. Equation (2) gives
`-f'(t)/n>=B*t*phi(t)` on `0<=t<=T`, with B the sum of these
rates including the directed pair multiplicities.
Using `phi(t)>.3985` gives the exact outer credit
`794046328543260351/40000000000000000000000>C` at `T=9/400`.
Evenness, monotonicity, same-threshold covariance comparison and safe
moves transfer this gain to every larger absolute threshold.

### Price transport

Set Z=sqrt(1-s)X+sqrt(s)Y, r=s/(1-s), and condition both pivots zero.
Ideal whitened middle/sibling drivers g remain independent standard,
independent of X_end,Y. Apply the same linear maps to Y to define H.
Drivers are g+sqrt(r)H(Y), with
Section 7's impulse errors. Conditional Y is centered. Fix sizes/labels
and both reflections; pairing g,-g uses cosh>=1, leaving
exp(-r||H||^2/2), of mass det(I+rV)^(-1/2), V=Cov(H|Z_end=0).
The centered tilt decreases error variances; Gaussian correlation
pays strip products. Only ideal coordinates factor; Y may repeat.
At zero the exact arm root/continuation densities are 2fD(x+2a) and
2phi_(1/sqrt(12))(b+sqrt(2)a+p/2), fD=density of |N(0,1/3)|.
Use Section 9's zero-offset closing margins and Section 7 impulses.
Pivot density is >=1/[2pi sqrt(1-w_d^2)].

For d>=3 let Q be full path precision, CE the endpoint covariance,
and S the conditional middle covariance with regression L. For
Bx=x_mid-L*x_end, M=B^T S^-1 B=Q-embed(CE^-1) is PSD and
tr(M Cov(X))=d-1.
Its adjacent entries are -12sqrt(2) at ends,
-18sqrt(2) inside; second neighbours have entry6; endpoint subtraction
adds w_d/(1-w_d^2). The source bound is
D2<=C*(D_left+D_right), C*=(sqrt(10)+3)/(2sqrt(2)). For fibre ratios
u,v>=1/2, u+v<=5/2 put A=sqrt(2u),B=sqrt(2v),S=A+B,P=AB,
z=(sqrt(10)+3)/2, psi=sqrt(uv)-(u+v+.5)/3. Then
6[C*(phi(u)+phi(v))-psi] equals
z(3S-4-S^2)+(2z-5)P+S^2+1. Min P gives
(z-1)(S-2)(3-S) on [2,3], or 3(S-3)[z-(S+3)/2] on [3,sqrt(10)],
both >=0. Use D2<=1/6 with alternating cap weights
0,theta,0,theta,..., theta=38-12sqrt(10). Adjacent terms cancel, giving
the actual trace <=(d-1)+2floor((d-1)/2)theta+2w_d/(1+w_d).
Schur/AM-GM/Hadamard bound the inverse determinant by
(1-s)^(1+(n+m)/2)/(1+mu_d*r)^((d-1)/2) for arm sizes n,m.
Use mu3,...,mu7=1.467,1.281,1.213,1.152,1.129.
For d2 H0=sqrt(33)[Y_w-(4sqrt(2)/11)(Y_u+Y_v)].
Var(H0)=1+(192/11)D2-24sqrt(2)(D_left+D_right)
<=571/11-16sqrt(10)<1.313. Clean means:

```
1+48H2(4sqrt(5)-5sqrt(2))/11 <1.04905,
1+24H2(C*-sqrt(2)) <1.111,
1+H2(24C*-2008sqrt(2)/95) <1.135,
1+H3(24C*-123sqrt(2)/7) <1.193,
1+(2/5)[4H3(12C*-18sqrt(2))+3/11] <1.117.
```

Joint root cases use Var(H)+3Var(R)/(2c^2)<=3/2, c=3-2sqrt(2),
H=sqrt(3)(Y_a-Y_b) at an endpoint.
Indeed for branch square roots A,B,W>=1/sqrt(2), sum squares3,
f=A+B+W-2sqrt(2) in [0,c], one has 2(A-B)^2+f^2/c^2<=1.
For f<=2c/5 use
(A-B)^2<=1/2-sqrt(2)f; otherwise use
(A-B)^2<=6-2(2sqrt(2)+f)^2/3 and convexity, with endpoint values
1 and (160-96sqrt(2))/25<1. Integrate. Choose h=Var(H) before tilt:
h<=p gives gain g=[1-(1-p)s]^(-1/2); otherwise use root variance
c^2 min(h0(s),1-2p/3). Baseline2 uses p=3/4 and min(gA,B)^2,
A,B the low/high side convolutions. Arms2 uses p=2/3 and the minimum
of LL,LH,HH EXTRA polynomials. D3 uses p=2/3 and minimum of four
complete middle integrals. Never add cases.

### D7–D10 and interfaces

Distance adjacency on a high-girth cubic cover obeys
`P1=t,P2=t^2-3,Pd=t*P_(d-1)-2P_(d-2)`. Put L=2sqrt(2).
For P7=t^7-13t^5+48t^3-44t its tangent residual at L is
(t-L)^2R(t), R=t^5+2Lt^4+11t^3+6Lt^2+56t+64L.
Here R'=5(t^2+4Lt/5)^2+(37/5)(t+30L/37)^2+632/37>0
and R(-3)=280L-708>0.
For z=t^2-8 the next polynomials are

```
P8=88+184z+94z^2+17z^3+z^4,
P9=t(48+120z+72z^2+15z^3+z^4),
P10=208+640z+508z^2+158z^3+21z^4+z^5.
```

Their tangent residuals are respectively

```
184(t-L)^2+z^2[(z+17/2)^2+87/4],
(t-L)^2[120(t+2L)+t(t+L)^2(z^2+15z+72)],
640(t-L)^2+z^2[76+t^2((t^2-3/2)^2+47/4)].
```

The odd expression is nonnegative for t>=0. For -3<=t<=0,
expanding P9/t at t^2=9/2 gives coefficients
`271/16,-17/4,-12,1,1`. Thus `|P9|<=9363/4<1920L`,
whereas its tangent `1968t-1920L<=-1920L`.
Let Sd be half the sum of excesses over directed distance-d pairs;
S1=sum_e D_e. All four residuals are nonnegative on [-3,3]. Pairing their
positive-semidefinite matrices with the physical row Gram matrix,
using unit diagonal and Pd(L)=3*2^(d-1)w_d, gives
`Sd>=multiplier*S1-o(n)` with multipliers 532,1040,1968,3620.
The last three are bounded by derivatives 368L,1968,1280L.
Finite-reference replacement has uniform o(n) error.

Use the alternating-path witness above, two endpoint contrasts and
optional size-two arms. All sign/arm classes and all distances are
disjoint. The common banded precision and trace formula give:

```
d   precision determinant   both-sign normalizer   pair density   mu
7   24541056                39.94                  .1665          1.129
8   305550144               56.2                   .1635          1.106
9   3748902912              78.3                   .1616          1.091
10  45586457856             109.2                  .1606          1.074
```

Endpoint determinants are 263/288,2183/2304,31/32,9047/9216.
The checker derives every entry from `9*12^(d-1)*(1-w_d^2)`,
the pi enclosure and the trace bound, without marginal factorization.

Retain 36 cells of width 1/80 per middle coordinate. The two-position
state starts with two Gaussian factors and endpoint-side arrays
`L0+(1-s)L2`. An interior insertion pays its new square, adjacent
and distance-two products and its strip; the final insertion uses
endpoint diagonal 36. After the eligible number of interior steps,
apply the final factor to a copy, close the side and endpoint strip,
then extend the interior state to the next distance. For each new
state divide the nonnegative sum downward at scale 10^18. All prior
factor floors survive. At time-panel right endpoint r retain the
spatial factors and integrate
`(1-s)^(1+(d-1)/2)/(1+(mu-1)s)^((d-1)/2)` over [l,r].
The convex denominator tangent at r and the decreasing weight's
average `r-s>=(r-l)/2` give a lower factor
`1+(d-1)(mu-1)(r-l)/[4(1+(mu-1)r)]`.
Integer powers integrate exactly; half powers use directed square roots.
Use endpoints `(j/24)^2`, j=2,3,4,5,6,7,8,10,12,14,17,20,23.
The four raw sums exceed floors `4e-7,76e-9,14e-9,2.65e-9`.
Multiplication by the spectral constants adds exactly
`65797/200000000` to K. No dirty-window debit is introduced.

certificate(): price.{adjacent,baseline2,arms2,long_paths,size_four}
and the D7–D10 contraction. Baseline2 excludes arms; arms2 uses only
2L0*a+a^2, a=(1-s)L2; size_four requires at least one4-arm.
Long paths include baseline/size2 once. Adjacent's overlap has curvature two.
Time nodes are (j/80)^2 for the adjacent and D2 baselines and
(j/24)^2 for the other paths. Section 9 of family.md gives the
cell-average densities for adjacent baseline and D3–D6.
D4 uses 11 global and 8 clean panels. An extra 2k arm
contributes (1-s)^k. For D4–D6 retain the strips at the right
endpoint r of [l,r] and integrate powers of 1-s, including the
factor (1-s)^(1+nu), nu=(d-1)/2. The tangent to the convex
(1+(mu-1)s)^(-nu) at r applies since each power weights
r-s by average at least (r-l)/2. Omitted tails contribute
nonnegative mass. Uniform source errors, times<1,
nonsingular ideal blocks, atomless boundaries and strict floor margins
fix one profile/truncation/jitter before graph size.

## 10. Signed and nine-key events on the actual graph

Put `eta=323/400000`, `theta=eta/2`, `T=9/400`,
`lambda=2sqrt(2)` and `k=(0,1,109/25,561/50,586/25)`.
Use the physical unit-variance field and defects
`A_v=sum_(u~v)Y_u-lambda Y_v`. Section 7 bounds their covariances by
incident budgets times `sqrt(2)^dist`, up to uniform errors.

Use the source lemmas: d-path excess is at most `k_d` times its cutoff,
and the sibling-contrast bounds below hold. Calibration is separate.

**Signed seven.** Use complete faithful stars:
`o:{a,b,u}`, `a:{o,a1,a2}`, `u:{o,w,w'}`, `w:{u,w1,w2}`.
All nine edges in these stars have excess at most eta. The desired
event has `a,b,a1,a2,w,w1,w2<t` and `o,u,w'>t`.
Condition on

```
O=Y_o=t+x, R=2Y_u-Y_a-Y_b, Q=sqrt(3)(Y_w'-Y_w).
```

Set `g=Cov((R,Q),O)`, `K=Cov((R,Q)|O)` and integrate
precision coordinates defined by `(R,Q)-gO=K(y,z)`, `x,y,z>0`.
The density is `(2pi)^(-3/2)sqrt(det K)` times
`exp(-[(t+x)^2+(y,z)K(y,z)^T]/2)`.
The source bounds give

```
|g_R|<=2eta, |g_Q|<=sqrt(3)*3.67eta,
K_RR>=1-8k2*eta-(2eta)^2,
K_QQ>=1-6k2*eta-(sqrt(3)*3.67eta)^2,
|K_RQ|<=sqrt(3)*11.05eta+(2eta)*(sqrt(3)*3.67eta),
K_RR<=1+2k2*eta, K_QQ<=1.
```

Thus `sqrt(det K)/(2pi)^(3/2)>1547/25000`. Its quadratic is bounded
above by `a*y^2+b*z^2+2c*yz`, where
`(a,b,c)=(1007041401,1000000001,15463245)/10^9`.

Mean upper rows / `10^5` have columns x,y,z,h.
Rows 1--4 mean score minus t; rows 5--6 mean t minus score.

```
a,b     94363 -15810   1583  315
a1,a2   83686 -21997   2500  917
w       83686  47438 -27827  917
w1,w2   71618  50820 -38984 1611
u      -94280 -32476    625  315
w'     -83333 -45193 -27827  917
```

Here h is a bias calibrated at `T0=11/200`. For t in
`[-T,0]` use `ceil(h*T/T0)/10^5` on the first four rows and zero
on the last two; reverse for t in `[0,T]`. Indeed the actual constant
is `(Cov(Y_v,O)-1)t`, reversed for guards, and covariance is at most one.
The ideal y coefficients are
`(-1/6,-lambda/12,lambda/6,1/2,-1/3,-lambda/6)` and ideal z
coefficients are `(0,0,-1/(2sqrt(3)),-lambda/(4sqrt(3)),0,
-1/(2sqrt(3)))`. Add the raw signed errors

```
y: (2k2,17.8,2,2k2,2k2,2k3)*eta,
z: sqrt(3)*(7.85,14.8,k2,10.53,1,k2)*eta,
```

and root-regression costs `c_v*|g_R|,c_v*|g_Q|`, where
`c_v=w_d+k_d*eta` at depths `(1,2,2,3,1,2)`. Round upward.
The first four x entries use this covariance upper bound; the last
two use the negative lower bound. Precision-coordinate means need no
inverse-K term.

Subtract the conditional means. Four defects and three half-differences
remain jointly centered Gaussian. Reconstruction, for defect radii
`r_o,r_u,r_a,r_w` and half radii `d1,d2,d3`, score losses at most

```
r_o/3+d1; 1.415*(r_o/3+d1)+r_a/2+d2;
.472*r_o+.5*r_u; .5*r_o+.708*r_u+.5*r_w+d3;
r_o/3; .472*r_o+.5*r_u.
```

Fix `r_o=.028,r_u=.0365`. Their joint strip is
`830353*926583/10^12=769390973799/10^12`. The third-radius
ceiling is `.009333334` and middle guard `.031466`.
The private recipes retain their full defect radii. Each ending closing
radius includes the extra `.5*(.028-.0235)=.00225`; its defect mass
is unchanged.
Integrate the centered normal with variance upper
`3*(3-2sqrt(2))*eta`, rounding SD upward at `10^-6`, to obtain them.
Half-difference variances are at most `1/12`, so standardized strip
width `3.46d` is conservative.

For the conditional centered pair `((Y_a-Y_b)/2,(Y_a1-Y_a2)/2)`,
the absolute density has lower bound `(84/11)exp(-(x^2+y^2)/(2*.0788))`
by Section 8 of `family.md`.
Use `|D1|<=u, |D2|<=c-(c-v)|D1|/u`, with
`u<=first-r_o/3`, `c<=second-r_a/2-1.415*r_o/3` and
`v<=c-1.415u`. Integration of the linear lower exponential
pays this trapezoid. Gaussian correlation combines it with
defect and ending strips. Each box uses one private recipe.

Use worst affine score corners and Jensen's cell-average exponential
bound. For root density the worse endpoint is t=0 or the signed ±T.
Transverse widths `(dy,dz)` are `(.0075,.03)` negative and
`(.0225,.04)` positive, inside `[0,3.6]x[0,6.8]`.
The root partitions in units of `.001` have endpoints

```
negative: 0,2,4,6,8,10,12,14,16,19,22,25,28,31,34,37,40,43,46,49,
          52,56,60,64,68,72,77,82,87,93,99,105,111,117,123,129,135,141
positive: 0,4,8,12,16,20,24,29,34,39,44,49,54,60,66,73,81,89,97,105
```
Each box uses one core/private/ending recipe, evaluated at its actual
worst margins. A cached choice may select a recipe, but never supplies
its mass. Separate density tables use the actual widths. Directed
initial division and twelve downward squarings at scale `10^18`
bound each exponential; strip integrals retain every floor.
The `1239545,167778` admitted boxes have disjoint interiors and each
signed sum exceeds `139358000/10^14`. Omitted boxes pay zero.

**Strict/relaxed nine.** Use faithful edges
`oa,ob,ou,aa1,aa2,bw,bw',uz,uz',ww1,ww2,zz1,zz2`.
Require `a,a1,a2,w,w1,w2,z,z1,z2<t` and `o,b,u,w',z'>t`.
Strict names have all thirteen edges at cutoff theta; relaxed names
allow eta only on four terminal child edges. Both are faithful.
Use raw coordinates:

```
O=t+x, R=Y_b+Y_u-2Y_a=y,
Q_b=sqrt(3)(Y_w'-Y_w)=z, Q_u=sqrt(3)(Y_z'-Y_z)=h.
```

For `g=Cov((R,Q_b,Q_u),O)` and its conditional covariance K, put
`gamma=(2theta,1.733*3.67theta,1.733*3.67theta)`.
The comparison matrix M has diagonal
`(1-8k2*theta-gamma0^2,1-6k2*theta-gamma1^2,
1-6k2*theta-gamma2^2)` and negative off-diagonals
`(1.733*24.55theta+gamma0*gamma1,
1.733*24.55theta+gamma0*gamma2,35.3theta+gamma1*gamma2)`.
Source bounds imply `K_ii>=M_ii`, `|K_ij|<=-M_ij` and `|g|<=gamma`.
Diagonal dominance/Neumann series give `|K^-1|<=W=M^-1` entrywise.
With `s=x+T`, `U=(y,z,h)+gamma*s`, the density is at
least `(63/2500)exp[-(s^2+U^T W U)/2]`. Hadamard gives the prefactor
from `det K<=1+2k2*theta`; no extra Jacobian occurs.

For any scalar L with `|Cov(L,O)|<=a` and raw driver covariance
bounds b, its conditional mean has absolute value at most
`a*s+(b+a*gamma)^T W U`. This constructs all ten emitted rows:
six defects at `o,a,b,u,w,z`, the root half `H0=(Y_b-Y_u)/2`,
and halves `Ha,Hw,Hz` of the three terminal sibling pairs.
For a defect at v use `p_j=1.415^dist(v,j)`,
`a=B*p_o`, `b=B*(max(p_b+p_u,2p_a),1.733max(p_w,p_w'),
1.733max(p_z,p_z'))`, with B=3theta, except B=5theta at w,z
in relaxed naming. For H0 use
`a=theta/2`, `b=theta*(k2,1.733*(1+k3)/2,1.733*(1+k3)/2)`.
Ha's two Q bounds are `1.733*14.8theta`; the opposite-ending
Q bound is `1.733*26theta` strict/`1.733*28.9theta` relaxed.
For an ending leaf v put p_d=k_d strict, `(d+1)k_d/d` relaxed.
Its remaining raw bounds are `a=theta*p_dist(v,o)/2`,
`b_R=theta*(p_dist(v,b)+p_dist(v,u)+2p_dist(v,a))/2` and
`b_Qb=1.733theta*(p_dist(v,w)+p_dist(v,w'))/2` (similarly Q_u),
except for the sharper opposite-ending bounds above.
For private v=a1 keep p_d=k_d and the two sharper 14.8 bounds.

Subtract these mean bounds in the algebraic star reconstructions.
Centered score losses, in order a; b,u; a1,a2; w,z; terminal pairs,
are `r/3; r/3+d0; .972r+da; .972r+1.415d0;
1.208r+e/2+1.5d0+(dw or dz)`.
Four core defects have SD at most `.01442`; the two ending defects
have that SD in strict naming and at most `.01862` in relaxed naming.
The four half variances are at most `1/12`.
For either centered conditional pair `(H0,Hw)` or `(H0,Hz)`,
the off-diagonal is at most
`(p2+p4)*theta/4+p3*theta^2/4+26*113*theta^2/.95<.0044`;
here p_d=k_d strict and p_d=(d+1)k_d/d relaxed.
Since `1/12+.0044<.09`, the same covariance enlargement and convex
trapezoid argument applies, now with slope 1.5. Keep the other ending
strip and all six defect strips; Gaussian correlation pays their
product with the selected pair probability. No residual independence
or simultaneous use of both alternative pairs is assumed.

Calibration uses rational arithmetic and a Machin pi bound for SDs,
radius masses and radicals. Each disjoint
cell pays one witness. Envelope and witness reflect under z↔h; factor
two needs no exchangeability of the actual law. On `[0,.2]x[0,3]x[0,4]^2`, the same supplied choices are
checked independently in both phases at `T=9/400`. Retain cells
whose zero-based ending indices sum to 22 through 43. Each phase
checks 58493 nonzero codes. Section 11 of family.md gives the density
moments. The strict and relaxed floors are
`1276/10^9` and `926/10^9`. These exceed `u+v9+v11`
and `u`, respectively, with the weights in Section 11.
Omitted cells pay zero; selection carries no probability assertion.

Nonsingular regression blocks, atomless sign boundaries, covariance
compactness and strict coefficient/sum gaps permit common profile,
truncation and positive jitter before graph size, uniformly over
eligible names and |t|<=T. With P9 strict and R9 relaxed eligible
names, P9<=R9. The nine-key saving is at least the strict floor
times P9 plus the relaxed floor times R9-P9.

The scalar source bounds are proved and replayed as described in
Section 4 of `family.md`.

## 11. Resource accounting for the central interval

The checker replays the 2250 bilateral cells proved in `family.md` §7.

Write eta=323/400000, D'=min(D,eta), H=(.0027,.006,.007),
lambda=141421/50000 and sigma=2021/1250. Tests use complete-source D
before truncation/jitter. At zero the global coefficient and increments are

```
K0=1118787644545099979211/25000000000000000000000
group             raw x         common a    dirty b1,b2,b3
adjacent          .0057376030    1           29,29,29
D2 nine-edge      .0023044849     2lambda     0,18lambda,18lambda
D3                .0002575579    19          0,327,327
D4                .0000455893    18lambda    0,522lambda,522lambda
D5                .0000077190    120         0,0,6040
D6                .0000009555    92lambda    0,0,7788lambda
stratified D2     .000377694     2lambda     0,18lambda,42lambda
K=K0+sum a*x=1851282581032699979211/25000000000000000000000.
```

Here 71931/10^10 moved from global adjacent credit into its clean
increment. Each row gives a*x*sum D', less eta*b_j*x per exceptional
edge in H1<D<=H2, H2<D<=H3 or D>H3. Projection cannot increase the
tree-name dirty bounds. The path-source sums give the common multipliers.
Local increments replace global floors; larger arms exclude smaller
paid arms. Section 10 of family.md lists the updated global and local floors.
The stratified increment has separate cleanliness tests.

Let Ej={e:D_e>Hj}; components connect incident edges. The component
calculation absorbs dirty losses. Put m=10^-12; x7 is stratified D2,

```
S=(21/2)eta*x_adj,
T=eta[lambda(6x2+189x4)+(231/2)x3],
t5=2220eta*x5, t6=2886lambda*eta*x6,
Ij=K0(sigma*Hj-4eta)/2,
Jjk=(K0/2)[sigma*Hj-kH1-(4-k)eta]_+,
Rj=K0(Hj-eta)-eta*sum_g b_jg*xg+min(S,Ij),
p0=R1+S-I1-m, P=p0+S,
tj=m-Rj-T (j=2,3), v0=max((t2+t3)/2,t3-t5),
v=v0-6lambda*eta*x7,
Jj*=K0[(sigma-1)Hj-3eta]/2.
```

### The 63/100 dyad source law

Put k=sqrt(2)/3, phi(t)=sqrt(t)-k(1+t), and
A(t)=sqrt(t(5t/2-1))-k(7t/2-1). Rationalization gives A/phi=d/a,
where d=sqrt(t)+k(1+t), a=sqrt(t(5t/2-1))+k(7t/2-1). Thus

```
d-ca=sqrt(t)(1-c sqrt(5t/2-1))+k[(1+c)+(1-7c/2)t].
```

Splitting at t=1, using k<1/2 and square-root bounds 5/4 below one,
25/16 through 11/8, 27/20 through 9/8, and 7/4 through 13/8 gives
A>=.63phi through 11/8, A>=.73phi through 9/8, A>=.53phi through
13/8, and A>=.5phi everywhere. The first two upper-piece margins
in d-ca are 7/3200 and 77/16000; for .53 both terms are nonnegative.
For z=5/2-t the rationalized numerators agree. Hence phi(z)>=phi(t)
for t>=5/4, and phi(z)>=1.3phi(t) for t>=13/8, using
t>=1.69z and 1+t>=1.3(1+z). Consequently
A(t)+.1phi(5/2-t)>=.63phi(t) for t>=11/8.

For h,r,s in [1/2,2] with h+r+s=3, if h,r<=11/8 use the first bound. Otherwise h>11/8,
r<9/8 and concavity gives phi(r)+phi(s)>=phi(5/2-h); combine the
last bound with A(r)>=.73phi(r). Thus, including endpoint zeroes by continuity,
A(h)+A(r)+.1phi(s)>=.63[phi(h)+phi(r)].

For a--b--c let h,r,s be the three ratios at b. At a, the sum
of its two outward ratios is `q^-1+2q-1/h` minus the atom correction
`delta_q*1[p=a]/P_a(p)`, where `q=a_profile^2,delta_q=q^-1+2q-3`.
Shift by q and use concavity and `phi_q(q)=0` to get the outward
payment `P_b(p) A_q(h)` with `A_q(h)=h phi_q(q+q^-1-1/h)`.
Clamping and the uniform Lipschitz bound charge at most `O(delta_q)`
after summing at each endpoint. Section 7's ratio correction to sum
three costs `O(1/2-q+delta_q)`. Integrate the displayed scalar
inequality to get `sum_O D+.1D_g>=.63(D_ab+D_bc)-o(1)`,
with any physical closing edge counted twice in O.

### Derivation of the dirty losses and their simultaneous payment

Copy complete physical rows/D to a cover of girth >twice the window
radius, retaining all tree names and repeated images. Each bad edge
has cover-degree many lifts; count upstairs and divide. Lift incident
edges together for union bounds: a projected new-site difference is
contained in the image of the upstairs difference. Thus short cycles
and aliases preserve every dirty bound.

Write L=2sqrt(2), R_v=sum_(u~v)b_u-Lb_v and S_v=sum_(e~v)D_e.
Uniformly, ||R_v^-||=o(1), <R_v,b_v>=S_v+o(1), and the two-step star
sum is >=L*S_v-o(1). The outward row sums, across w--u, are

```
U1=L*b_u-b_w+R_u,
U2=6b_u-L*b_w+L*R_u+sum_(z~u,z!=w)R_z.
```

Expand U1/U2 inner products between different branches, subtract the
ideal-row expansion, and collect a residual only when all its required
centres remain. This gives the following positive terms and losses.
Here a site carrying q*S_v loses at most 3q*eta; lambda<L replaces L.

```
group  retained terms                      dirty sites      loss/eta
adj    D_e                                 29 edges         29
D2     lambda*S_w                          6 centres        18lambda
D3     9D_ab; 5S_a                         13 edges;14 sites 327
D4     6lambda*S_w;3lambda*S_u             14;30 sites      522lambda
D5     44D_ab;24S_u;14S_z                  29 edges;30;62   6040
D6     24lambda*S_w;14lambda*S_u;8lambda*S_z 30;62;126        7788lambda
```

For example, D3's four-pair sum is
9D_ab+<R_a,Lb_b-b_a>+<R_b,Lb_a-b_b>+<R_a,R_b>+o(1).
Collecting at a gives <R_a,5b_a+L R_a> >=5S_a-o(1).
For D4 the primary term is L^2*sum D2-2L*S_w>=6L*S_w;
the collected R_u bracket is 3L*b_u+6R_u+L*sum_(w~u)R_w.
D5's primary term is 44D_ab; its collected R_u and R_z brackets are
24b_u+6L R_u and 14b_z+4L R_z+6sum_(u~z)R_u.
D6's primary term is 36sum D2-12L*S_w>=24L*S_w. Its L*R_u
bracket is 14b_u+4L R_u+6sum_(w~u)R_w; its R_z bracket is
8L*b_z+20R_z+4L*sum_(u~z)R_u+6sum_(u~z,w~u,w!=z)R_w.
All omitted directional brackets are nonnegative near criticality by
b_u>=alpha*b_w; residual products have vanishing negative error.
Counting each S_v twice over edge incidences yields common multipliers
2lambda,19,18lambda,120,92lambda, respectively.

The edge ball of line-graph radius r has 1+4+...+2^(r+1) edges.
Endpoint half-balls through depth r have 2(1+2+...+2^r) vertices.
Consequently D3 uses 13 middle edges and their 14 endpoints; D4 uses
14 centres and their 30 neighbours; D5 uses 29 middle edges, their
30 endpoints and 62 neighbours; D6 uses 30 centres and one/two further
layers 62,126. D2's nine-edge test fails at 2+4=6 middle vertices.
Stratified D2 adds eight outer positions: its debit is 18lambda above
H2 plus 24lambda above H3, hence (0,18lambda,42lambda). Longer arms
keep global residual bounds outside these windows; their event support
does not enlarge the local test.

For incident exceptional edges the respective numbers of NEW sites are

```
adj:8; D2:2; D3:4,4; D4:4,8; D5:8,8,16; D6:8,16,32.
```

Order a connected component's edges along a spanning tree. Subtracting
these union bounds from the independent dirty debits saves, per additional
edge, 21*x_adj, 12lambda*x2, 231*x3, 378lambda*x4, 4440*x5,
5772lambda*x6, all times eta. Thus adjacent gives 2S(k-1) on E1
components; D2--D4 give 2T(k-1) on E2 components; D5 gives
2t5(k-1) on E3 components. Only whole two-top E1 components spend D6's `2t6` saving.
Other D6 union savings are discarded.

For stratified D2 let `C_r(e)` be the vertices within r of an endpoint
of e. Its dirty centres lie in `C_1(e)` for middle bad edges and
`C_2(e)` for top edges, of sizes at most 6 and 14. Two incident edges
have three endpoints and at most five outgoing edges; their one- and
two-neighborhoods have sizes at most 8 and 18. Thus two middle, two
top, or mixed edges save at least 4,10,6 centres, respectively.
For a mixed pair `C_1(middle)` lies in `C_2(top)`. Each centre costs
at most `3lambda*eta*x7`. Use the mixed refund `18lambda*eta*x7`
only on whole top/middle E1 dyads. In each nonsingleton E2 block of a
mixed E1 component select one incident pair and use the uniform refund
`12lambda*eta*x7`, giving the stated v. Selected pairs have disjoint
edge sets; pure dyads and mixed blocks belong to different components.
Partition remaining edges into singletons. Summing within-part union
bounds overestimates all dirty centres, even if parts overlap. These
refunds spend only stratified D2 debits, independently of D2–D6 refunds
and outside radial energy.

The same outside energy pays all strata. Select isolated E2 edges, isolated class-one E1 edges, and
one maximum-D edge from each pure higher dyad (the top edge in a mixed
dyad). These form a matching. The paired law sum_off D>=sigma*D-o(1)
then gives I_j for isolated E1 edges, J_jk for isolated E2 edges with
k lower E1 neighbours, and J_j* per pure higher dyad. For J_jk retain
the actual lower excesses and transfer (K0/2)(D_f-H1) from each lower
neighbour: there are at most two accesses to its unused K0(D_f-H1).
Outside E1, each edge is accessed at most once at either endpoint,
so its discarded energy pays all half-charges together.

Pure class-one dyads instead use the stated 63/100 source inequality
on their four outer incidences and 1/10 of the centre's third edge.
Subtracting the caps gives [(63/50)H1-(41/10)eta]_+. Their accessed
vertices lie in separate whole E1 components; including them preserves
the same multiplicity-two outside charge. This gives `Gamma=(K0/2)[(63/50)H1-(41/10)eta]_+-m`.

The isolated first-class upgrade also survives these added charges.
An isolated edge unshared by all selected outside accesses receives
full weight, hence an extra I1. A shared edge is joined by an outside
edge to an exceptional edge in another component. Their 29-pivot balls
intersect in 13 (distance two). Add shared isolated edges to anchored
components, or build an unanchored component of p>=2 isolated edges:
each adds <=16 pivots, saving >=13p/2 beyond all within-E1 savings.
Thus the uniform extra is min(I1,13eta*x_adj/2)=I1. Full outside
charges and these extra union savings are used on disjoint classes.

It remains an exact scalar component reduction. The checker verifies
I1<S<I2<I3, p0<=0, 0<=t2<=t3<=S, and the displayed v0,v.
After leaving m per edge, pure class-one components have Pl-2S:
size two uses 2p0+Gamma>=0; sizes >=3 use 3p0+S>=0 and P>=0.
An isolated first-class edge uses R1+I1>=m. Higher singletons use
Ij-S-T>=tj; higher dyads use J2*/2>=t2,
J3*+18lambda*eta*x7>=t2+t3 and
t5+t6+J3*/2>=t3. Larger pure higher components leave
(S+T)(h-2)-sum tj, nonnegative when max(tj)<=(S+T)/3.

In a mixed component let l count lower edges, q_jk higher singleton
blocks of degree k, and r nonsingleton higher blocks. Its lower bound is

```
F=sum (S+J_jk-T-tj)q_jk+2(S-v)r+Pl-2S.
```

For a nonsingleton block, if it has a class-two edge its remaining
amount is >=2S-t2-t3; if all its edges are top it is >=2S-2t3+2t5.
These give `2(S-v0)`; the stratified refund raises it to `2(S-v)`. Connectivity gives q+r>=1 and
sum k*q_jk+r<=2l. For l=1 enumerate (q,r)=(1,0),(2,0),(0,1),(1,1),(0,2):
with L_j=J_j1-tj the tests are

```
p0+L_j-T>=0, P/2+L_j-T>=0, P/2-v>=0,
p0+2S-T+L_j-2v>=0, (p0+3S)/4-v>=0.
```

For l=2 each coefficient is >=-2p0>=0, cancelling Pl-2S=2p0.
For l=3 use beta=p0/2+S/6 and k<=3; for l>=5 use
beta=p0/2+3S/10 and k<=4. The checked beta>=0,
S+J_jk+k*beta-T-tj>=0 and 2(S-v)+beta>=0 give
F>=(P-2beta)l-2S>=0, the latter remainder being `(2S/5)l-2S`
for l>=5.

For l=4 put beta=p0/2+S/4, mu_jk=S+J_jk+k*beta-T-tj,
L4=max(0,-mu_24,-mu_34), and
C4=min(4beta,2(S-v)+beta,mu_jk for j=2,3 and k<=3).
The checker verifies C4>=L4 and beta>=0. At most one degree-four
singleton E2 block can occur: two consume all eight lower incidences,
forcing all four links between their two high edges. With the high
edges this is a whole K4, whose complete resolvent has diagonal
`(1-2beta_source)/(1+beta_source)` and off-diagonal
`beta_source/(1+beta_source)`, beta_source=q/(1+2q^2). All six
adjacent covariances, hence
excesses, are equal, contradicting the different E1/E2 classes.
This is exact for every finite profile on the simple kernel.
Now write

```
F=sum mu_jk*q_jk+[2(S-v)+beta]r
                  +beta[8-sum k*q_jk-r].
```

If there is one degree-four singleton its loss is <=L4. Another block
pays >=C4; if none exists, the degree slack is four and pays 4beta.
If there is no degree-four singleton all terms are nonnegative.
This exhausts components, including cycles.

Boundedly many bounded-norm rows give error <=C0*delta*n, with C0
independent of graph/cover degree. Choose profile/truncation/jitter
before graph size. The m reserve pays bounded-multiplicity exceptions;
the final strict slack pays the remainder. We obtain K*sum D'-o(n),
using each global/local increment, outside half-charge and union refund once.

Section 9's disjoint D7–D10 cylinders add 65797/200000000 to
the capped coefficient. The actual crossing-cover expectation is maximal at zero
by Section 1, so this payment holds throughout the central interval.
The marked budget is

```
A=eta(K+65797/200000000)=600620827548562093285153/10^28.
```

Now let n_ab count vertices with a high and b medium incident edges:
high is D>eta, medium eta/2<D<=eta. Let sS count low marked--marked
edges, Ui count unmarked vertices with i marked neighbors, and J count
U1--U1 edges. Write rho=2X/A, c=(sigma-3/2)A/4,
h=(sigma-1)A/4, g=2(z-c-rho*A/4)/3, B=c+g/2.
The budget `(K+65797/200000000)*sum D' + 2z*sS + X*U2 + (3X/2)*U3 +
(X-A/5)*J` pays the sum of the following vertex coefficients, up to
the already incurred `o(n)`.

```
01 A(1/2+rho/4)+B    02 A/2-g          03 3A/4-3B
10 3A/4+h+z         11 3A/4+h-B        12 A-h-2B
20 sigma*A/2        21 5A/4-2h-B       30 3A/2-3h.
```

Section 6 of `family.md` derives the table and joint medium/type11
reserve. Normalize `cbar=c/A`, `gbar=g/A`, `Bbar=B/A` and put

```
omega=1/4, k=33/100, d=omega*(sigma-1)/2,
alpha0=1/[1+16(omega/2-cbar)/k],
alpha1=(Bbar+k/16)/(omega/2+k/16),
alpha=alpha0+(908475033730/10^12)(alpha1-alpha0),
pay=A[d*alpha+(k/8)(1-alpha)].
```

Full internal low retention gives the following additional endpoint
payments, including Section 12 of family.md with b10=55903822/10^14:

```
10:b10, 11:A*d*alpha, 20:E-2pay,
12:-pay, 21:-2pay, 30:-3pay.
```

The exact checks include `0<alpha0<alpha1<1`, `g>=0`,
`2/5<rho<1/2`, `1/8+rho/4<z/A<1/4+rho/4`, `z/A>=omega/2`,
and

```
min((1/2+rho-omega)(sigma-3/2),z/A-omega/2)>d*alpha,
1444/111375<=E/A<2rho*.0585,
E/A<=min(rho*(sigma-3/2),(2-sigma)/2,(5-3sigma)/4+pay/A).
```

The first inequality removes the high 10--11 fee. The direct type20
bonus uses its own mixed-star reserve and the unused type11 debit
at type20 centers with only one H1 high neighbor. Section 12
bounds b10 using reserves left after these payments.

Let L7 and LW count the eligible seven-key and weak-pair names
with weights 1/2 and 2/3. Write L9=L11=P9 and LR=LA=R9.
The incidence proof in Section 5 of family.md gives the following
lower bounds. For j=a+b define

```
c7=(0,18,12,6), c9=(0,35,22,9), c11=(0,45,28,11),
cW=(0,88/3,56/3,8),
cR(a,b)=(0,36,24,12)_j+8(3-j)1[a>0],
cA(a,b)=c9_j-2(3-b)1[a=0<b].
```

The lower counts are respectively (6,9,11,12,9,8)*n minus sum c*n_ab
(c7 uses a; c9,c11,cW use j), plus (sS,U2,U3,J) coefficients

```
(0,0,0,0), (26,12,30,4), (34,16,40,6),
(24,12,24,0), (22,12,24,2), (64/3,0,0,0).
```

Negative affine bounds remain valid. Seven/weak weights 1/2,2/3 and
guard/name disjointness are included. Take

```
(q,u,v9,v11,gap,X,E)
=(139357900,89691594,36451375,100,100,1446536023,337910715)/10^14,
e=u-gap, w=12719/10^12,
z=13v9+17v11+12u-e+32w/3.
```

Pay qL7+v9L9+v11L11+(u-e)LR+eLA+wLW. Require signed>=q,
strict>=u+v9+v11, relaxed>=u. L9,L11 are alternative counts sharing
one probability budget; the relaxed counts have total weight u.
Section 6's root-sign/neutral-root disjointness makes savings additive.

Its sS coefficient is exactly 2z. The three overlap checks verify that
its positive U2,U3,J coefficients cover X,3X/2,X-A/5. The ten marked
checks verify that each adjusted source coefficient covers

```
q*c7_a+v9*c9_j+v11*c11_j+(u-e)cR(a,b)+e*cA(a,b)+w*cW_j.
```

Adding these inequalities cancels every debit and allowance once.
Price plus forward saving is at least Cn-o(n), where

```
C=6q+9v9+11v11+12u-3e+8w=1981610721/10^14.
```

For |t|<=9/400 bound the base expectation by its value at zero,
then apply forward savings at that t. Fix source errors below the
marked-check margins and Section 5's graph slack before graph order.
The central reduction equals C; Section 5's outer/reference bounds give
the strict graph claim.

