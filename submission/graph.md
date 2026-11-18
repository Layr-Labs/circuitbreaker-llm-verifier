# Gaussian crossing matchings and their conversion to physical cutwidth

This file gives the graph and finite path-integration arguments used by
`proof.md`. A cut counts physical edges with multiplicity; loops never
cross. Throughout the probabilistic argument the base is connected, simple
and cubic. The final section extends the result to every cubic multigraph.
All random rows used on the base have fixed finite support before its order
tends to infinity. Covers provide comparison fields, not a replacement
ordering source on the base.

## Crossing objective and Gaussian identities

Given vertex bits b, let F(b) be the minimum vertex cover number of the
graph of edges whose endpoint bits differ. This crossing graph is
bipartite. Its maximum matching has size F(b): starting with a maximum
matching, follow alternating paths from unmatched vertices on the zero
side; the unreachable zero vertices and reachable one vertices form a
cover of size exactly the matching. Absence of an augmenting path proves
that every crossing edge is covered.

For edge bits x define

    Q(b,x)=sum_v NAE(b_v,(x_e:e incident with v)),
    C(x)=sum_v NAE((x_e:e incident with v)).

NAE is zero on empty or constant lists and one otherwise. Then
`F(b)=min_x Q(b,x)`: zero-cost vertices force all their edges to their
own color, so the other vertices cover the crossing edges. Conversely a
crossing cover permits the requirements of all uncovered vertices to be
satisfied simultaneously. NAE is submodular: a mixed difference is -1
when all other arguments agree and zero otherwise (for two arguments it
is -2). Sums and identified arguments preserve this property. Minimizing
over x also preserves it, by applying the meet/join inequality to two
minimizers. Thus F is symmetric and submodular. Changing one vertex bit
changes F by at most one, since it changes just one summand of Q.

If X,Y are centered Gaussian vectors with unit diagonal and
`Cov(Y_i,Y_j)>=Cov(X_i,X_j)` off the diagonal, then

    E F(1[Y<t]) <= E F(1[X<t]) <= E F(1[X<0]).       (G1)

Here and below vectors may have repeated coordinates. For the first
inequality, extend F multilinearly to the cube and convolve each threshold
with a small independent Gaussian. The off-diagonal second derivatives
are nonpositive averages of mixed differences. With independent X,Y,
Gaussian integration by parts on `sqrt(1-s)X+sqrt(s)Y` gives derivative
one half the covariance difference contracted against that Hessian.
The diagonal terms cancel, and all remaining terms are nonpositive.
Pass to the sharp threshold by bounded convergence of each continuous
unit Gaussian marginal.

For the second inequality first add positive independent coordinate noise.
Moving the orthant boundaries gives

    f'(t)=phi(t) sum_i E[D_iF(1[X_-i<t]) | X_i=t].

For t>=0, the conditional vector `U_j=X_j-t` has means
`m_j=(Cov(X_i,X_j)-1)t<=0`. Its reflection `-U+2m` has the same law
and is at most -U coordinatewise. The reflected threshold bits dominate
the complementary bits. D_iF decreases in every other bit, while symmetry
gives `D_iF(1-b)=-D_iF(b)`. Hence each displayed conditional mean is
nonpositive. f is even, proving (G1); remove the extra noise by a limit.

Two quantitative versions will be used. For `Z_s=sqrt(1-s)X+sqrt(s)Y`,
with C,A the covariance matrices of Y,X respectively,

    d/ds E F(1[Z_s<0])
      =sum_(i<j)(C_ij-A_ij) phi_(ij,s)(0,0)
                   E[Delta_ij F | Z_(s,i)=Z_(s,j)=0].   (G2)

This follows from the same smoothing argument: integration of the two
threshold derivatives leaves the bivariate density and the four-corner
mixed difference. The two ordered Hessian terms cancel the factor 1/2.
Apply first on compact subintervals of s<1, where reference coordinate
noise gives a positive covariance matrix; bounded convergence passes to
the endpoint. All terms discarded from (G2) are nonpositive.

For one fixed positive definite field with covariance C, put
`W_j=X_j-C_ij X_i`, `k_j=1-C_ij`. The preceding threshold derivative can
equivalently be written

    f'(t)=phi(t) sum_i g_i(t),
    2g_i(t)=integral_(-t)^t sum_(j!=i) k_j f_(W_j)(k_j s)
              E[Delta_ij F | X_i=X_j=s] ds.          (G3)

The other bits in the last expectation are tested at s. To derive this,
condition on X_i=t, use symmetry to subtract the marginal at -t, and
move the remaining conditional thresholds from -k_j t to k_j t one
coordinate at a time. Differentiation contributes k_j and a mixed
difference. Conditioning `W_j=k_j s` is precisely the original exterior
law given `X_i=X_j=s` after adding the deterministic means `C_iv s`.

## A deterministic layout lemma

For every vertex order of a loopless cubic graph,

    cw(H) <= max_(prefix b) F(b)+6.                 (G4)

The proof below allows the optimizer of Q to change on an arbitrarily
large set at a single vertex step. Clamp an edge by
`b_u AND b_v<=x_uv<=b_u OR b_v`. This does not change min Q: changing a
same-color edge to its forced color cannot destroy an all-equal local
list. The minimizers in a clamped cube are closed under meet and join.
Choose the least one x(b). For b<=c, meet/join of the least minimizers
are feasible at b,c respectively; joint submodularity and optimality
force both to remain minimizers. Thus `x(b)<=x(c)`.

We need an interpolation fact. Suppose x<=y are clamped at b and y is
optimal. Among vertices incident to changing edges classify the endpoint
Q-costs as A:(0,1), B:(1,0), P:(1,1). The case (0,0) is impossible.
Vertices in P have a permanent disagreement among fixed arguments.
Changing edges join opposite anchors, so no changing edge joins two A
or two B vertices. If q is the fixed contribution elsewhere, then
`Q(b,x)=q+|P|+|B|`.

The A--B changing-edge graph has a matching covering A. Indeed, resetting
to zero in y all changing edges incident to a set S of A saves |S|,
and can cost at most its number of B-neighbors. Optimality gives Hall's
inequality. Hall's conclusion follows from an alternating search in a
maximum matching: an uncovered A either has an augmenting path or yields
a deficient neighbor set.

Flip first the P--P and B--P edges. Next process each A of changing
degree three as one consecutive block. It has no fixed edge and is
unsplit before and after its block; at most one extra split occurs
during it. All remaining A have changing degree at most two. Regard an
A with two different B-neighbors as an edge of a multigraph on B; one
with a single B-neighbor gives a half-edge, also when its two changing
edges are parallel. Its possible A--P edge stays in its block. Hall
ensures that a component has at most as many edge/half-edge objects as
B-vertices. Therefore it is a tree with at most one half-edge, or a
unicyclic graph without half-edges. Peel tree leaves away from a
half-edge; finish the half-edge last. Peel attached trees of a unicyclic
component and then traverse its cycle consecutively. At completed blocks,
processed A-blocks minus finished B-vertices is at most one: each peeled
edge finishes a B; only the first cycle edge can fail to do so. Inside
a block the difference is at most two. Split cost is consequently at
most `q+|P|+|B|+2=Q(b,x)+2` throughout.

At one vertex-prefix step b->c, first clamp x(b) at c. This flips at
most three edges and costs at most `F(b)+1` in Q_c. Then interpolate
to x(c) by the preceding fact. The resulting monotone edge order has
split cost at most `M+3`, where `M=max F(prefix)`; its endpoints are the
all-zero and all-one fields. Number edges by this order. The bag
`B_i={v:first(v)<=i<=last(v)}` has size at most M+5, since only the two
ends of edge i can be additional to the split vertices. Order vertices
by their median incident edge positions. At a vertex cut with last
median theta, a crossing edge below theta charges its outside endpoint's
unique lowest edge, and one above theta charges its inside endpoint's
unique highest edge. No vertex receives two charges. Both kinds lie
in B_theta; the single edge at theta costs at most one more. This proves
(G4), including parallel edges.

## Covers, radial rows and concentration

Every finite connected graph has finite covers of any prescribed girth.
Here is a construction. Choose a spanning tree and label the remaining
oriented edges by free generators. A closed nonbacktracking walk gives
a nonempty reduced word: an inverse cancellation would enclose a closed
tree walk or an immediate backtrack. A reduced word of length k acts
nontrivially on a finite set: on points 0,...,k prescribe each letter's
permutation to carry the corresponding successive points, using the
inverse prescription for inverse letters. Reducedness prevents conflicting
domains or images, so complete the partial injections to permutations.
The word carries 0 to k. Take a product of these representations for the
finitely many words up to the girth target, and let Q be the finite
generated permutation group. Use fibers Q, trivial edge action on the
spanning tree and right multiplication by each label on other edges.
No forbidden word fixes a fiber point in the regular action. Thus the
cover has the desired girth and is finite; a connected component suffices.

In a degree-d cover, repeated bits satisfy `F_up=d F_base`. The upper
bound pulls back a cover. A cover upstairs gives fractional vertex
weights downstairs by division by d. For every crossing edge their sum
is at least one. Summing on a maximum bipartite matching proves the
opposite inequality.

On the universal cubic tree use complete nonnegative unit rows

    f_v(z)=sqrt(c_alpha) alpha^dist(v,z),
    c_alpha=(1-2alpha^2)/(1+alpha^2),
    7/10<=alpha<1/sqrt(2).

There are `3*2^(r-1)` points on a positive radius-r sphere. Their squared
sum gives the normalization above. Counting off-path cones gives

    r_d(alpha)=alpha^d[(d+1)-(d-1)alpha^2]/(1+alpha^2),
    r_d(alpha)->w_d=(1+d/3)2^(-d/2).

Push the squared rows through the tree covering map and take square roots
of the fiber masses to obtain unit base rows b_v. Define adjacent excess
`D_uv=<b_u,b_v>-r_1(alpha)>=0` and `E=sum_e D_e`.
Fiberwise Cauchy--Schwarz proves both this nonnegativity and covariance
domination from the independent cover field to the repeated base field.

The actual keys use radius-R truncated, renormalized rows pushed to the
base and then positive independent vertex jitter:

    Y_v=sqrt(1-epsilon^2)<b_v^R,g>+epsilon xi_v.

All g coordinates and xi coordinates are independent standard normals.
The full-alpha omitted squared mass is a geometric tail, uniform in
every graph and root. Taking norms on each fiber contracts row distance.
Thus fixed R sufficiently large and epsilon sufficiently small make
every pair covariance uniformly close to that of the complete rows.
On a sufficiently high-girth cover each fixed neighborhood used in the
proof is a genuine tree window. Repetition includes the base jitter:
equal images add its covariance and preserve domination. This also
handles singular repeated-key endpoint laws in (G1)-(G2).

All local rules in this submission inherit an original key within
distance three. A vertex key is inherited at most
`1+3+6+12=22` times. Once R is fixed, their threshold contributions depend
only on independent noises in bounded neighborhoods, so variances of
sums of bounded local terms are O_R(n). The original F itself has
variance O_R(n) by bounded differences in the independent Gaussian input
coordinates and its one-vertex Lipschitz bound: resampling one input
can alter only O_R(1) keys. The same is true after the local rule.

For finitely many threshold mesh points, Chebyshev gives simultaneous
o(n) deviations with positive probability. Within a mesh interval, the
change is bounded by 22 times the number of original keys in that interval;
their expectation is at most n times the interval length divided by
sqrt(2pi), and their count has variance O_R(n). Gaussian tail counts
handle thresholds outside a fixed bounded interval. Choose that cutoff,
then a fine fixed mesh, before n. Resolve groups of equal inherited keys
arbitrarily; their size is at most 22. Apply (G4). Consequently any
uniform all-threshold expected bound `(gamma+arbitrarily_small_error)n`
for these corrected keys gives `cw(H)<=(gamma+o(1))n` on every sufficiently
large connected simple cubic base, uniformly in its choice.

## A finite path upper bound for the reference objective

Let N_j(D) be half the number of directed nonbacktracking j-walks in a
subcubic graph D. For odd J put

    A_J(D)=sum_(j=1)^J (-1)^(j+1)N_j(D),
    C_J=1+3J 2^(J-1).

If B is the set of degree-three vertices, then
`F(D)<=A_J(D)+C_J|B|`. Indeed, after deleting B, paths and cycles have
cover number at most their odd alternating path-count sum: on an r-edge
path it is ceil(r/2) when r<J, and r-(J-1)/2 otherwise; on a cycle it
is its order. Restoring B costs |B|. A specified vertex can occur in at
most `(j+1)3*2^(j-1)` directed j-walks, so all changes of signed counts
are covered by the stated C_J.

The limiting reference field has covariance w_distance. At a cubic star,
`sum_neighbors X=2sqrt(2)X_root` almost surely, by expanding its variance
using `w_1=2sqrt(2)/3,w_2=5/6`. A degree-three crossing star therefore
has zero probability at threshold zero. Fixed finite profile windows
converge to this law. For any fixed odd J the reference expectation is
thus at most `(Gamma_J+arbitrarily_small_error)n`, where

    Gamma_J=sum_(j=1)^J (-1)^(j+1)3*2^(j-2)p_j

and p_j is the probability of alternating signs on a j-edge tree path.
The finite certificate below uses J=47 and proves

    Gamma_47 < U=26129513284633533/200000000000000000.   (G5)

On a tree path, covariance recurrence gives the stationary Gaussian
autoregression `X_i=sqrt(2)X_(i-1)-X_(i-2)/2+epsilon_i/sqrt(12)`.
After alternating signs the positive-state densities are

    f2(x,y)=c2 exp[-9(x^2+y^2)/2-6sqrt(2)xy],
    k(z|x,y)=c1 exp[-6(z+sqrt(2)y+x/2)^2],
    c2=3/(2pi), c1=sqrt(6/pi).

Thus `p_j=2 integral_(positive orthant) f2 product k`.
The identity `f2(x,y)k(z|x,y)=f2(y,z)k(x|z,y)` permits anchoring the
product at any adjacent pair. A single positive extension has mass at
most 1/2. Two extensions have mass at most
`atan(1/sqrt(2))/(2pi)<1/10`, by decreasing their nonnegative starting
values to zero and integrating the wedge `U>0,V>sqrt(2)U`.
The last strict inequality follows from
`cos(2atan(1/sqrt(2)))=1/3>cos(2pi/5)`.
Write `H(k)=10^-floor(k/2) 2^-(k mod 2)`.

Truncate the n-dimensional integral, n=j+1, to [0,L]^n. A tail on either
coordinate of the anchored pair is at most

    Tpair=c2 exp(-9L^2/2)/(54sqrt(2)L^2),

by integrating first `exp(-6sqrt(2)xy)` in the non-outlying coordinate,
then using `integral_L^infty exp(-9x^2/2)dx<=exp(-9L^2/2)/(9L)`.
For an outlier at i choose anchor start a=min(i,n-2). Its other two arms
cost H(a)H(n-2-a). Hence the single-orthant omitted mass is at most
`Tpair sum_(i=0)^(n-1) H(a)H(n-2-a)`.

Use N-node Gauss--Legendre quadrature on [0,L]. The precision diagonals
of the n-coordinate density are `(9,9)` for n=2, `(12,33,12)` for n=3,
and `(12,36,39,...,39,36,12)` for n>=4. Completing the square in coordinate
i bounds its amplitude by `f_n(0)=c2*c1^(n-2)`. The Gaussian characteristic
function and its 2N-th moment give derivative bound
`f_n(0) q_i^N (2N)!/(2^N N!)`. For completeness, integrating the Gaussian
by parts gives its characteristic function `exp(-q_i t^2/2)` and moment
recurrence `M_(2k)=q_i(2k-1)M_(2k-2)`; taking the absolute value under
the Fourier integral proves this derivative bound.

The one-dimensional Gauss error is
`L^(2N+1)(N!)^4/[(2N+1)((2N)!)^3]` times the 2N-th derivative bound.
This coefficient follows by integrating the squared monic orthogonal
polynomial in the Hermite-interpolation remainder: the interpolation
error has one repeated zero at every node, and positivity leaves
`h^(2N)(xi)/(2N)!` times that squared polynomial. The Legendre norm and
leading coefficient yield the displayed coefficient after rescaling.
Telescoping the tensor product gives single-orthant error

    E_n= f_n(0)L^(n+2N)(N!)^3 sum_i q_i^N
                     /[(2N+1)2^N((2N)!)^2].

This can be divided by `(4L)^max(0,n-9)`: for each varying coordinate i,
at least n-9 internal coordinates j remain at distance more than two.
Their independent quadratic factors `exp(-39x_j^2/2)` survive outside
the one-coordinate derivative. Both their integral and their positive
quadrature mass are below 1/4 (the latter is checked). Pulling out these
factors in every mixed tensor product proves the refinement.

The program uses N=93,L=3, 100-digit outward Decimal intervals and exact
fixed-point tensor contraction with denominator S=2^224. Legendre roots
are isolated in disjoint rational sign-changing intervals. At x=a/b the
integer recurrence `t_0=1,t_1=a,t_k=(2k-1)a t_(k-1)-(k-1)^2 b^2 t_(k-2)`
has the sign of P_k(x). These N brackets isolate all roots. Derivatives
and weights are enclosed by `(x^2-1)P_N'=N(xP_N-P_(N-1))` and
`w=2/[(1-x^2)(P_N')^2]`. The bound `||P_N''||<=N^4` controls the whole
bracket: `||P_k||<=1` follows from the Legendre integral representation,
and telescoping `P_(k+1)'-P_(k-1)'=(2k+1)P_k` bounds the first and second
derivatives by N^2,N^4.

Only O(N^2) factors are needed. Set

    A_ij=exp(-3x_i^2/2-6sqrt(2)x_i x_j),
    B_ik=exp(-6x_i x_k),
    C_jk=c1 w_k exp(-12x_j^2-6x_k^2-12sqrt(2)x_j x_k).

Starting from `M_ij=w_i w_j f2(x_i,x_j)`, update
`M'_jk=C_jk sum_i (M_ij A_ij)B_ik`. The implementation stores one lower
array and an explicit upper bound for its lost L1 mass. It checks A,B<=1,
`max_j sum_k C_jk<=1`, and initial mass<=1. The C-row check follows from
`c1 sum_k w_k exp(-6x_k^2)<1`; the remaining exponential factors are at
most one. Thus this positive tensor update has L1 operator norm at most
one. Let DA,DB,DC bound the maximum factor-enclosure width in grid units.
At a step, rounding M*A down loses at most N^2/S; replacing A and B by
their lower factors loses at most DA/S and DB/S since state mass and
the C-row mass are at most one. Replacing C loses at most N*DC/S, since
the summed intermediate mass over its output index is at most N.
Rounding the next array down loses at most N^2/S. Therefore add
`(DA+DB+N*DC+2N^2)/S` to the L1 deficit per step, starting with the sum
of the initial entry widths. Adding this accumulated bound to the
computed lower mass encloses the exact tensor quadrature. Every product
and deficit is an exact integer calculation; it is a directed enclosure,
not an estimated floating rounding error. Pi is enclosed by Machin's arctangent identity;
Decimal square roots and exponentials are expanded outward one representable
step after correctly rounded evaluation. Final p_j bounds double the
quadrature error and include the one-sided tail above. Negative terms
of Gamma_J reverse interval endpoints. This is the finite predicate
checked by `certify.py`, rather than a floating approximation to (G5).

## Uniform extension to physical cubic multigraphs

Suppose the simple-cubic conclusion with constant gamma has been proved
as above. For each epsilon>0 it gives a constant C_epsilon such that
all connected simple cubic graphs K satisfy
`cw(K)<=(gamma+epsilon)|K|+C_epsilon`, absorbing finitely many small K.

For a connected bridgeless loopless subcubic block, contract an edge uv
when either endpoint has degree at most two or there are parallel uv
edges; delete all r>=1 internal edges. The new degree is
`deg(u)+deg(v)-2r<=3`, and edges-minus-vertices decreases by r-1.
If no such move remains, contract a triangle. Before that move the graph
is simple cubic; its contraction is cubic and decreases edges-minus-
vertices by one. Alternate paths show bridgelessness is preserved.
Termination leaves one vertex or a simple triangle-free cubic kernel K
with `|K|<=2(nu(B)-1)`, where `nu=M-N+1` is cycle rank. Such a kernel is
also 2-connected: a cut vertex in a bridgeless cubic graph would have
at least two edges into each of two components.

Record contractions as binary/ternary rooted trees of pieces, each with
original boundary at most three. Expand each kernel vertex depth first,
visiting children in decreasing size. Before a missing leaf the earlier
pieces include at most two per level; whenever one exists the chosen
child has at most half its parent's leaves. At most `2ceil(log2 N)`
pieces are therefore exposed, contributing at most `6ceil(log2 N)`
crossing edges. Expansion adds only this amount to kernel width.

Delete bridges of any connected loopless subcubic graph. Their blocks
form a tree, and cycle ranks add. A block with attached components has
width at most its width plus the maximum attached width plus three,
by listing each attached component just before its attachment vertex.
Choose a centroid block, so every attached component has at most half
the vertices (move into a component larger than half until none remains).
Induction on the resulting recursion gives, for every epsilon>0,

    cw(G)<=(2gamma+2epsilon)nu(G)
                         +D_epsilon ceil(log2(|G|+1))^2.

Indeed each block costs at most `(2gamma+2epsilon)nu(B)` plus
`6ceil(log2 N)+C_epsilon`; the depth is at most `floor(log2 N)+1`.
For a cubic multigraph on D vertices, remove its loops. Its cycle rank
is at most D/2+1 and its cuts are unchanged. Thus uniformly over connected
physical cubic multigraphs `cw<=gamma D+o(D)`. Let e(d) be the largest
positive excess on connected d-vertex cubic multigraphs, or zero if there
are none, and set `h(D)=max_(d<=D)e(d)`. The finite maximum exists because
degree bounds all multiplicities. The uniform estimate implies
`h(D)/D->0`: split the maximum at a fixed large d_0, bounding its large-d
part by epsilon D and its finite small-d part by a constant. For disconnected
graphs concatenate component orders, taking their maximum width. This
proves one nonnegative nondecreasing physical remainder h, as required
by `family.md`, with the same gamma for all physical cubic graphs.
