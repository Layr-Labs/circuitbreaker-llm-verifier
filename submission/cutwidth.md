# Cubic cutwidth at most 0.1341 times the number of vertices, plus o(n)

## 0. Statement, provenance and the finite numerical premise

An undirected multigraph may have parallel edges and loops. A loop contributes
two to degree and never crosses a vertex cut. A graph is cubic when every
degree is three. For an ordering of its vertices, its width is the maximum
number of edges with one endpoint in a proper prefix and the other outside.
The cutwidth cw(G) is the minimum width of an ordering. Put

    gamma = 1341/10000 = 0.1341.

**Theorem 1.** There is a nonnegative nondecreasing function h on the positive
integers, with h(D)/D tending to zero, such that every cubic multigraph K on D
vertices satisfies cw(K)<=gamma D+h(D). Disconnected graphs, loops and parallel
edges are included. This is hypothesis H(gamma,h) in `circuit.md`.

The proof first treats simple bridgeless cubic graphs using a perfect matching
and its complementary cycles. Gaussian keys are assigned to matching edges;
the edges on bounded cycle blocks are inserted jointly to minimize the number
of split vertices. The minimum is a symmetric submodular function of the fixed
edge sides and equals a vertex-cover number on the transitions of a binary
word. Edge-rooted Hellinger lifting compares this quantity on every graph with
the cubic tree. Short cycles are handled by periodic unrolling, with no penalty
per short cycle. The final sections extend the result to all cubic multigraphs.

**Credit.** The route from cubic cutwidth to circuit lower bounds, Gaussian
edge orderings, and the kernel compression used below are due to McGuire,
Schlesinger and Williams. The Hellinger pushforward comparison for vertex
profiles appears in Arthur Ostapenko's CircuitBreaker submission
5fef109bf7f1b664295244a0dee6df066ed2232d. Here it is applied to edge-rooted
profiles. Joint insertion along complementary cycles, the transition-cover
bound and the periodic-unrolling argument give the improvement in this file.
All nonstandard graph and comparison arguments needed here are proved below.
Petersen's classical theorem is used in precisely this form: every finite
bridgeless cubic graph has a perfect matching.

The sole numerical premise is the following finite Gaussian inequality, proved
by the interval certificate in `certificate.md` and its declared computation.
Let (X_0,...,X_5) be a centered Gaussian vector with covariance C_ij=c_|i-j|,
where

    c_0=1,
    c_d=(d+1+2sqrt(2))/(4*2^(d/2))             for 1<=d<=5.

For 1<=j<=5, let p_j be the probability that all j consecutive pairs
(X_0,X_1),...,(X_(j-1),X_j) have opposite signs. Define

    g_5=p_1-p_2+p_3-p_4+p_5.

**Numerical fact C.** g_5<=gamma.

No empirical graph experiment or floating-point estimate is used as a
mathematical premise. The rest of this file proves that fact C implies
Theorem 1, including the uniformity of its lower-order term.

## 1. Symmetric submodular Gaussian comparison

For a function F on {0,1}^q, its mixed difference at distinct coordinates i,j
is F(s^11)-F(s^10)-F(s^01)+F(s^00), the other coordinates being fixed. It is
submodular when all mixed differences are nonpositive, equivalently
F(a)+F(b)>=F(a meet b)+F(a join b) for all a,b. It is symmetric when
F(1-a)=F(a). Write s_t(x)_i=1[x_i<t].

**Lemma 1.1.** Let X,Y be centered Gaussian vectors with unit coordinate
variances, and let F be symmetric submodular.

(a) If Cov(Y_i,Y_j)>=Cov(X_i,X_j) for every i!=j, then
E F(s_t(Y))<=E F(s_t(X)) for every t.

(b) For every t, E F(s_t(X))<=E F(s_0(X)).

**Proof of (a).** The multilinear extension is

    Ftilde(p)=sum_b F(b) product_i p_i^b_i(1-p_i)^(1-b_i).

For i!=j, its mixed derivative is the expectation of the mixed difference of
F at independent Bernoulli coordinates with parameters p, hence is
nonpositive. Put sigma_delta(u)=Phi(-u/delta), where Phi is the standard
normal distribution function, and define

    f_delta(z)=Ftilde(sigma_delta(z_1-t),...,sigma_delta(z_q-t)).

This function is smooth and bounded with bounded first and second
derivatives. Its off-diagonal mixed derivatives are nonpositive: they are the
corresponding derivative of Ftilde times two nonnegative normal densities
divided by delta^2.

Take X,Y independent and Z_r=sqrt(1-r)X+sqrt(r)Y. For any such smooth f,

    E f(Y)-E f(X)
      = (1/2) integral_0^1 sum_(i,j)
          [Cov(Y_i,Y_j)-Cov(X_i,X_j)] E partial_i partial_j f(Z_r) dr.  (1)

To justify (1), write jointly Gaussian vectors as linear combinations of
independent standard normals. One-dimensional integration by parts gives
E[U_i H(V)]=sum_j Cov(U_i,V_j) E partial_j H(V). Differentiate E f(Z_r)
for 0<r<1 and apply this identity separately to X and Y. The resulting
derivative is the integrand in (1). It is bounded, and expectations extend
continuously to the endpoints, so integration proves the identity. Equal
variances make the diagonal summands zero. Each other summand is nonpositive
by covariance domination and submodularity. Apply this to f_delta and let
delta decrease to zero. The coordinates have continuous unit normal
marginals, so none equals t with positive probability; bounded convergence
proves (a). Joint covariance matrices may be singular.

**Proof of (b).** First replace X by sqrt(1-eps^2)X+eps Z with independent
standard coordinate noise Z; it is enough to prove the claim for this
positive definite covariance and then let eps decrease to zero. The
expectation is differentiable, with derivative

    phi(t) sum_i E[D_i F(s_t(X_-i)) | X_i=t],                (2)

where D_iF(b)=F(1,b)-F(0,b). Formula (2) follows by integrating the joint
Gaussian density over the orthants with common moving boundary t, one
boundary coordinate at a time. Density tails justify differentiation.

For t>=0 condition on X_i=t and put U_j=X_j-t for j!=i. Its conditional
mean is m_j=(C_ij-1)t<=0. The vector U'=-U+2m has the same conditional
distribution as U and satisfies U'<=-U coordinatewise. Its bits 1[U'<0]
therefore dominate 1-1[U<0] almost surely. By submodularity D_iF is
nonincreasing in the other bits, and symmetry gives D_iF(1-b)=-D_iF(b).
Consequently the conditional expectation of D_iF is nonpositive. Thus (2)
is nonpositive for t>=0. Symmetry makes the expectation even in t, proving
the maximum occurs at zero. The noise limit follows by bounded convergence.
QED.

**Lemma 1.2 (partial minimization).** If G(b,x) is submodular then
F(b)=min_x G(b,x) is submodular. If G is invariant under complementing all
bits, so is F. For fixed b the set of minimizers in x has a least element in
the coordinatewise order, and this least minimizer is nondecreasing in b.

**Proof.** Choose minimizers x,y for b,c. Submodularity gives

    F(b)+F(c)>=G(b meet c,x meet y)+G(b join c,x join y)
             >=F(b meet c)+F(b join c).

Complementation bijects the feasible sets. For a fixed b, meet and join of
two minimizers are minimizers, again by submodularity and optimality. The
meet of the finite nonempty set of minimizers is its least element. If b<=c
and x,y are least minimizers, the same submodular inequality and optimality
force x meet y to minimize at b. Leastness implies x<=x meet y, hence x<=y.
QED.

## 2. The cover and edge-rooted Gaussian keys

For now H is a finite connected simple cubic graph with a perfect matching
M. Its other edges form disjoint cycles. Let pi:T->H be its universal cover.
Here is the construction and the needed facts. Fix a base vertex v. A vertex
of T is a nonbacktracking walk from v, including the empty walk; two are
adjacent when one extends the other by one edge. Each nonempty walk has one
parent and two extensions, and the empty walk has three extensions. Thus T
is the cubic tree. The map pi sends a walk to its endpoint and each directed
tree edge to the corresponding directed graph edge. At every tree vertex
the three incident directed edges map bijectively to the three at its image.
Hence every graph walk lifts uniquely from any specified preimage of its
initial vertex. A nonbacktracking graph walk lifts to a simple tree path.

If tree vertices a,b have the same image, there is an automorphism theta of
T preserving pi and carrying a to b. For each x, take the tree path from a
to x, project it to a nonbacktracking graph walk, and lift that walk from b;
define theta(x) to be its endpoint. Adjacency and pi are preserved. The
analogous construction from b to a inverts it. Thus any two lifts of an
oriented edge are carried to each other by such an automorphism, and the
automorphism permutes every vertex fibre of pi. These are the only cover
facts used below. The preimage of M is a perfect matching of T. A factor
cycle lifts to a doubly infinite nonmatching path, with a distinct matching
edge leaving the path at each vertex.

For an integer R>=1 put A_R=1+3R/2 and define tree vertex profiles

    a_u(z)=A_R^(-1/2) 2^(-dist_T(u,z)/2) 1[dist_T(u,z)<=R].

They are nonnegative unit vectors: the radius-j sphere has 3*2^(j-1)
vertices for j>=1, contributing 3/2 to squared norm before normalization.
Write r_d(R)=<a_u,a_v> when dist_T(u,v)=d. In particular r_0(R)=1. For a
tree edge e=uv use the nonnegative unit edge profile

    K_e(z)=(a_u(z)+a_v(z))/sqrt(2+2r_1(R)).

It is symmetric in u,v and equivariant under tree automorphisms. For a graph
edge e choose a lift e_tilde and push its squared profile to H:

    P_e(x)=sum_(z:pi(z)=x) K_e_tilde(z)^2.

The sum is finite and independent of the lift, because the automorphisms
above permute each fibre and carry one profile to the other. It sums to 1
and is supported on graph vertices at distance at most R from an endpoint
of e. For 0<eps<1 let g_x,xi_e be independent standard normals on graph
vertices and graph edges, and set

    Z_e=sqrt(1-eps^2) sum_x sqrt(P_e(x))g_x+eps xi_e.          (3)

Every Z_e is standard normal. Distinct graph edges have distinct scores
almost surely, since their difference has variance at least 2eps^2. On T
use independent vertex and edge noises to define

    Z^T_e=sqrt(1-eps^2) sum_z K_e(z)g'_z+eps xi'_e.           (4)

Each defining sum has finite support.

**Lemma 2.1 (edge lifting).** For any two distinct tree edges e,f,

    Cov(Z_pi(e),Z_pi(f))>=Cov(Z^T_e,Z^T_f).

**Proof.** In each fibre, Cauchy--Schwarz yields

    sqrt(P_pi(e)(x)P_pi(f)(x))
      >=sum_(z:pi(z)=x) K_e(z)K_f(z).

Sum over x and multiply by 1-eps^2. If pi(e)=pi(f), the graph covariance
also includes eps^2, whereas the distinct tree edges' extra noises are
independent. This only strengthens the inequality. Both variances are 1.
The argument applies simultaneously to every pair in any finite tree-edge
family, even if several edges project to the same graph edge. QED.

The square in P_e is taken after summing the two TREE vertex profiles. A sum
of independently normalized graph vertex profiles is not being used.

## 3. The stationary tree process

**Lemma 3.1.** For every fixed d>=0,

    r_d(R)->w_d=(1+d/3)2^(-d/2) as R tends to infinity.

**Proof.** For d=0 the assertion is exact. For d>=1 take the tree path
u=q_0,...,q_d=v. Vertices attached to q_i at depth j>=1 contribute the
unnormalized product 2^(-d/2-j), when both endpoint distances are at most R.
There are 2^j such vertices at each endpoint i=0,d and 2^(j-1) at each
interior i. Each permitted j range differs from 1,...,R by O(d) terms. The
two endpoint branches contribute (2R+O(d))2^(-d/2), the d-1 interior
branches contribute ((d-1)R/2+O(d^2))2^(-d/2), and path vertices contribute
O(d)2^(-d/2). Divide by A_R=(3/2)R+1. For fixed d the limit is the stated
(d+3)/3 times 2^(-d/2). QED.

Along a lifted factor path let e_i be the matching edge at its i-th vertex.
These edges are distinct. For |i-j|=d>=1 their four endpoint distances are
d,d+1,d+1,d+2. The stationary Gaussian process (Z^T_e_i) therefore has
unit variances and covariances

    c_d^(R,eps)=(1-eps^2)
      [r_d(R)+2r_(d+1)(R)+r_(d+2)(R)]/[2+2r_1(R)].           (5)

As R tends to infinity and eps to zero, this tends to

    c_d=[w_d+2w_(d+1)+w_(d+2)]/[2+2w_1]
        =(d+1+2sqrt(2))/(4*2^(d/2)), d>=1;     c_0=1.        (6)

The last identity follows by substituting w_d and multiplying out. The
d>=1 formula must not be used at d=0: the key at one position has variance
1. Every finite matrix (c_|i-j|) is positive semidefinite, being a limit of
the finite-radius covariance matrices. Thus these formulas define the
finite Gaussian vectors required in fact C and in the limits below.

Let p_j^(R,eps) be the probability of j consecutive sign transitions in
the process (5) at zero, for 1<=j<=5, and let

    g_5^(R,eps)=p_1^(R,eps)-p_2^(R,eps)+p_3^(R,eps)
                         -p_4^(R,eps)+p_5^(R,eps).

Finite-dimensional Gaussian convergence gives p_j^(R,eps)->p_j and
g_5^(R,eps)->g_5. Indeed the event boundary is contained in a union of
coordinate-zero hyperplanes, each of probability zero since its coordinate
variance is 1. This proves convergence even if a limiting joint covariance
were singular. No quantitative convergence rate is needed in this theorem.

## 4. Optimized insertion on a path or cycle

For three bits define h(a,b,c)=1[the three bits are not all equal]. It is
complement symmetric and has every mixed difference equal to -1, so is
submodular. Sums of these penalties are submodular; identifying coordinates
preserves the meet/join inequality and hence submodularity.

Given a binary word b=(b_1,...,b_k), define

    F_k(b)=min_(x_0,...,x_k) sum_(i=1)^k h(b_i,x_(i-1),x_i).

The analogous cyclic minimum C_k(b) imposes x_k=x_0. In the path, x_0,x_k
are free end bits. By Lemma 1.2 both F_k and C_k are symmetric submodular.

**Lemma 4.1 (transition covers).** F_k(b) is the minimum vertex-cover size
of the path consisting of the edges i--(i+1) with b_i!=b_(i+1). C_k(b) is
the minimum vertex-cover size of the corresponding transition subgraph of
the cycle.

**Proof.** An unpenalized vertex forces its incident path or cycle bits to
its label b_i. Consequently endpoints of a transition edge cannot both be
unpenalized: penalized vertices cover every transition edge. Conversely,
given a vertex cover, assign every edge incident to a vertex outside the
cover that vertex's label. These assignments agree whenever both endpoints
are outside the cover, because their connecting edge is not a transition.
Fill remaining bits arbitrarily. At most the cover vertices are penalized.
This proves both equalities. QED.

A run of r consecutive transition edges is a path with minimum cover
ceil(r/2), proved by an alternating cover and a matching of ceil(r/2)
disjoint edges. Put t_i=1[b_i!=b_(i+1)] for a path. Products extending beyond
its last edge are zero, and set

    A_5(t)=sum_i [t_i-t_i t_(i+1)+t_i t_(i+1)t_(i+2)
                      -t_i...t_(i+3)+t_i...t_(i+4)].

On a maximal run of r ones this sum is ceil(r/2) for r<=5 and r-2 for
r>=5. The second quantity is at least ceil(r/2). Therefore F_k<=A_5
pointwise. In any stationary binary process let p_j be the probability of
j consecutive transitions. For k>=5,

    E F_k<=E A_5
      =(k-1)p_1-(k-2)p_2+(k-3)p_3-(k-4)p_4+(k-5)p_5
      =k g_5+B,
    B=-p_1+2p_2-3p_3+4p_4-5p_5<=p_2+p_4<=2.               (7)

The first inequality for B uses p_1>=p_2, p_3>=p_4, p_5>=0. Thus

    E F_k<=k g_5+2.                                        (8)

Equations (7)-(8) also hold for the finite-radius process (5). Comparison
is always applied to the submodular F_k; the truncation A_5 itself need not
be submodular.

## 5. Short cycles have no density penalty

For a cyclic word b of length k, let b^m be the word repeated m times and
read on a path of mk vertices. Then

    m C_k(b)-1<=F_(mk)(b^m)<=m C_k(b).                      (9)

To prove this, first regard b^m as a ring. If some edge is not a transition,
its transition components are path components repeated m times, with total
cover number m C_k(b). If every edge is a transition, the cyclic binary
word alternates, so k is even; the original and repeated ring both have
cover number half their vertex count. Finally remove the single wrap edge
to get the path. Deleting one edge changes minimum cover number by at most
one. This proves (9), also for odd k.

Fix a factor cycle of the graph H, of length k. Unroll it m times along a
lifted nonmatching tree path. Its successive matching edges are distinct
in T even if their projections repeat. The projected Gaussian matching
labels at threshold t are the same graph word repeated m times. Lemmas
2.1 and 1.1 applied to the symmetric submodular F_(mk), followed by (8),
give for mk>=5

    m E C_k-1
      <=E F_(mk)(graph repeated labels at t)
      <=E F_(mk)(tree labels at zero)
      <=mk g_5^(R,eps)+2.

Divide by m and let m tend to infinity. For every threshold t, every graph
factor cycle, and every length k,

    E C_k<=k g_5^(R,eps).                                  (10)

There is no girth assumption. Unrolling is only an expectation argument;
the final construction operates on bounded graph blocks.

## 6. Bounded blocks and the vertex ordering

Fix an integer L>=5. Keep each factor cycle of length less than L as one
block, with all its cycle edges free. For each other cycle of length k,
write k=qL+r with 0<=r<L and divide its vertices into q contiguous paths:
one of length L+r and q-1 of length L. Designate the cut edges as fixed and
all other cycle edges free. When q=1, designate one closure edge fixed, so
the block is still a path. Matching edges are always fixed. Every block
has fewer than 2L vertices, and there are at most n/L path blocks. A fixed
edge has endpoints in at most two blocks.

Use the Gaussian scores (3) on fixed edges. At threshold t their bits are
1[Z_e<t]. Independently within each block, minimize the sum of h over its
vertices by choosing the free-edge bits. Choose the least minimizing free
vector. Lemma 1.2 says these vectors are nondecreasing as t increases.
When all boundary bits are zero or all are one, the unique zero-cost
extension is respectively all zero or all one: each block vertex has a
fixed matching edge, and zero cost forces agreement with it. Thus every
free edge changes exactly once, at one of its block's fixed key values.
Give it that value as its entry time. Fixed edges use their own key as
entry time. Order all edges by entry time, breaking ties arbitrarily.

At a threshold t the number of vertices with at least one edge strictly
below t and at least one edge at or above t equals the sum of block minima;
write this number as Cost(t). The strict comparison convention agrees with
the least minimizer before a key changes; the free edge flips immediately
after its chosen key value.

**Lemma 6.1 (edge bags to cutwidth).** If the edges of a loopless cubic
graph are in an order with positions 1,...,m, define

    B_i={v: first(v)<=i<=last(v)},

where first and last are the smallest and largest of its three incident
edge positions. Ordering vertices by the position of their median edge,
with fixed tie-breaking, gives cutwidth at most max_i |B_i|+1.

**Proof.** Fix a nonempty proper vertex prefix and let theta be the median
position of its last vertex. For a crossing edge xy with x in the prefix:
if its position is less than theta, charge it to the outside vertex y.
The median at y is at least theta, so this is its unique lowest edge and
y lies in B_theta. If the edge position is greater than theta, charge it
to x. The median at x is at most theta, so this is its unique highest
edge and x lies in B_theta. No vertex is charged twice: prefix vertices
only receive high-edge charges and outside vertices only low-edge charges.
At most one edge has position exactly theta. The cut size is therefore at
most |B_theta|+1. QED.

Fixed keys are pairwise distinct almost surely. Any free edge whose entry
time is a given fixed key belongs to one of that key's at most two blocks.
These blocks contain fewer than 4L vertices. Including the fixed edge
itself adds at most two endpoints. For a bag at position i with entry time
tau, a vertex in the bag either is split at threshold tau or has an
incident edge entering at tau. Hence

    |B_i|<=Cost(tau)+4L+2,
    cw(H)<=sup_t Cost(t)+4L+3<=sup_t Cost(t)+12L+10.          (11)

The last deliberately looser bound will absorb all finite tie costs.

## 7. Expectations and concentration at all thresholds

For a whole-cycle block use (10). For a path block of k vertices, first
allow its two fixed cut edges to have arbitrary sides. This reduces its
optimal cost to F_k on its matching labels. Restoring the prescribed
sides changes the penalties only at its two endpoints, increasing cost
by at most two. Its matching edges lift to consecutive off-path edges of
the universal-cover factor path. By Lemmas 2.1 and 1.1 and (8), its
expected cost at every t is at most k g_5^(R,eps)+4. Summing over blocks,

    E Cost(t)<=[g_5^(R,eps)+4/L] n for every real t.          (12)

We now prove that a single realization achieves this leading bound
simultaneously at all thresholds. All parameters L,R,eps below are fixed
before the graph size.

A block term lies in [0,2L] and depends only on the independent noises of
fixed scores with an endpoint in that block. Vertex noises for such a
score lie at distance at most R+1 from the block. Choose one representative
vertex per block. Since its block has diameter at most 2L, the dependence
neighborhood has radius at most R+2L+1 around this representative. Two
block terms are independent when their representatives are farther apart
than 2R+4L+2. A radius-s ball in a cubic graph has fewer than 3*2^s vertices.
There are at most n representatives, and covariance in absolute value of
two variables in [0,2L] is at most L^2. Thus

    Var Cost(t)<=3 L^2 2^(2R+4L+2)n=:K_1 n.                (13)

For an interval I let W(I) count fixed edges whose scores lie in I. There
are at most 3n/2 such edges, each standard normal, so

    E W([a,b))<=(3n/2)(b-a)/sqrt(2pi),
    E W((-infinity,-T)), E W([T,infinity))
                  <=(3n/2)Phi(-T).                        (14)

There is a constant K_2 depending only on R with Var W(I)<=K_2 n,
uniformly over intervals. For an explicit loose bound assign each fixed
edge to one endpoint; at most three edges are assigned to a vertex.
Let Y_v be their interval indicator sum, in [0,3]. It depends on vertex
noises in radius R+1 about v and on the independent noises of its assigned
edges. Distinct Y_v are independent when their vertices are farther than
2R+2 apart. Hence

    Var W(I)<= (27/4)2^(2R+2)n=:K_2 n.                     (15)

Changing the side of one fixed edge changes the total optimized cost by
at most two. Use the old free-edge assignments as feasible assignments
after the change: only its endpoint penalties can change, by at most one
each; then reverse the comparison. Consequently, for s<t,

    |Cost(t)-Cost(s)|<=2 W([s,t)).                          (16)

The cost is zero at the all-zero and all-one fixed assignments, so the
costs in either tail are bounded by twice the corresponding tail count.

For precision, let delta>0. Choose T so that 3Phi(-T)<delta/8 and a finite
mesh -T=t_0<...<t_J=T with every gap at most delta sqrt(2pi)/24. Consider
the finitely many random variables Cost(t_j), the J interval counts
W([t_j,t_(j+1))), and the two tail counts. By (13),(15), Chebyshev's
inequality, and the union bound, for all sufficiently large n there is a
realization in which each variable differs from its expectation by less
than delta n/16. Intersect this event with pairwise distinctness of fixed
keys, which has probability one. Equations (12),(14),(16) then imply

    sup_t Cost(t)<= [g_5^(R,eps)+4/L+delta]n.                (17)

Indeed in an interior mesh interval its cost is at most the cost at the
left endpoint plus twice that interval count. The expectation contribution
of twice the count is at most delta n/8 and the deviations contribute
3delta n/16, less than delta n in total. Either tail has cost less than
delta n/4. The interior bound is nonnegative, since it bounds an expected
nonnegative cost, so these tail bounds are covered after adding delta n.
Using a larger sufficiently-large-n threshold if needed does not change
any preceding choices.

**Theorem 7.1.** For every eta>0 there is N_eta such that every connected simple
cubic graph with a perfect matching on n>=N_eta vertices satisfies

    cw(H)<=(gamma+eta)n.

**Proof.** First choose L>=5 with 4/L<eta/4. By Section 3 and fact C,
choose R and 0<eps<1 such that g_5^(R,eps)<=gamma+eta/4. Apply (17) with
delta=eta/4. Finally choose n sufficiently large for its concentration
statement and for (12L+10)/n<eta/4. Equation (11) completes the proof.
All choices depend only on eta, not on the graph. If eta is large the same
argument applies, or it follows from a smaller positive eta. QED.

## 8. Kernel compression

A kernel is a finite simple triangle-free 2-connected cubic graph. It is
bridgeless: a bridge would have a component containing at least one other
vertex on either side, since its endpoints have degree three in a simple
graph; removing an endpoint would disconnect it. Petersen's theorem
therefore supplies a perfect matching. Theorem 7.1, with finitely many
smaller sizes covered by the trivial bound cw<=3n/2, implies: for every
eta>0 there is C_eta>=0 such that all kernels satisfy

    cw(H)<=(gamma+eta)|V(H)|+C_eta.                         (18)

For a connected loopless multigraph B with N vertices and M edges define
its cycle rank mu(B)=M-N+1.

**Lemma 8.1.** If B is connected, bridgeless, loopless and subcubic on N
vertices, either cw(B)<=6ceil(log2 N), or there is a kernel H on
h<=2(mu(B)-1) vertices such that any width-q ordering of H expands to a
width at most q+6ceil(log2 N) ordering of B.

**Proof.** A one-vertex B has width zero. Otherwise repeatedly apply the
following rules, preferring the first whenever possible:

1. If adjacent p,q are joined by r>=1 edges, and either one has degree at
most two or r>=2, contract {p,q} and delete its r internal edges.
2. If rule 1 does not apply and three vertices form a triangle, contract
the triangle and delete its three internal edges.

Rule 1 produces degree deg(p)+deg(q)-2r<=3 and decreases M-N by r-1>=0.
Before rule 2 all degrees are three and the graph is simple. Each triangle
vertex has one outgoing edge, so the contraction has degree three and
decreases M-N by one. These contractions reduce vertex count and preserve
connectedness, looplessness and bridgelessness. For the last assertion an
alternate path avoiding any surviving edge projects to a walk avoiding it.

At termination, either one vertex remains, or H is cubic, simple and
triangle-free. Then h/2=M_H-N_H<=M_B-N_B=mu(B)-1. It is 2-connected: if
deleting v left at least two components, each would send at least two
edges to v, since a unique such edge would be a bridge. This would give
degree at least four.

Record every contraction by a rooted tree with original vertices as
leaves and two or three children at each internal node. A subtree's leaf
set is a piece. At its creation at most three original edges left it,
because the contracted vertex had degree at most three. Its original
boundary never changes.

List a piece's leaves depth first, visiting children in nonincreasing
order of their leaf counts. In a proper prefix let l be the first missing
leaf. Along the root-to-l path the prefix is the union of earlier siblings'
pieces. There are at most two earlier siblings at each level. Whenever one
exists, the chosen child has at most half its parent's leaves. This occurs
at most ceil(log2 N) times. Thus the prefix is a union of at most
2ceil(log2 N) pieces and at most 6ceil(log2 N) edges leave it.

If one contracted vertex remains, this proves the first alternative.
Otherwise expand the vertices of H in its specified order using these
leaf orders. A prefix consists of whole pieces forming a prefix of H and
part of one more piece. A crossing edge either crosses that H-prefix or
leaves the partial piece. This gives the asserted additive width. QED.

## 9. Assembling bridge blocks

The blocks of a connected loopless subcubic graph G are the components
after deleting all bridges; singleton blocks are allowed. Fix a block B
and let F_1,...,F_m be the components of G-V(B).

**Lemma 9.1.** Every F_j joins B by exactly one bridge e_j, at a vertex
v_j in B, and

    mu(G)=mu(B)+sum_j mu(F_j).

If B has width w and the F_j have widths w_j, then G has width at most
w+max_j w_j+3, with an empty maximum interpreted as zero.

**Proof.** An edge between vertices of B belongs to B: if it were a
bridge, a path within B would complete a cycle through it. Every edge
leaving B is a bridge. There is at least one from each F_j by connectedness;
two would give a cycle through either one using paths within B and F_j.
Thus it is unique. Every edge lies in B, in an F_j, or is an e_j. Counting
vertices and edges proves the mu identity. B itself is bridgeless, because
a cycle of G uses no bridge and therefore stays in a single block.

Walk along the ordering of B, and immediately before a vertex v list the
orderings of all F_j attached there. In a proper prefix let u be its first
unlisted B-vertex. At most w internal B-edges cross. At most one F_j at u
is partially listed, contributing at most max w_j crossing edges. All
other F_j are complete or untouched. A bridge crosses only when its
B-endpoint is u, giving at most three more edges. If all B vertices have
been listed, so has the graph. QED.

**Lemma 9.2 (centroid block).** Some block B has every component of
G-V(B) of size at most |V(G)|/2.

**Proof.** Start with a block. If a component F after its removal is
larger than half, it is unique. Move across its unique joining bridge to
the adjacent block B' contained in F. The component of G-V(B') containing
the old block is V(G)\V(F), of size less than half. A new component larger
than half must therefore be a proper subset of F. Repetition strictly
decreases the size of the unique large component, so terminates. QED.

**Proposition 9.3.** For every eta>0 there is D_eta such that every connected
loopless subcubic graph G on N vertices satisfies

    cw(G)<=(2gamma+2eta)mu(G)+D_eta ceil(log2(N+1))^2.        (19)

**Proof.** Put c=2gamma+2eta and K_N=6ceil(log2 N)+C_eta+3, where C_eta
comes from (18). Lemma 8.1 shows that any bridgeless block B of size at
most N has width at most c mu(B)+K_N-3: either it already has logarithmic
width, or expand a kernel of h<=2(mu(B)-1) vertices using (18).

Induct on N'<=N to prove every connected subcubic G' of size N' has width
at most c mu(G')+(floor(log2 N')+1)K_N. Choose its centroid block B. Each
hanging component has size at most N'/2, so its induction overhead is at
most floor(log2 N')K_N. Lemma 9.1 bounds the new width by

    c mu(B)+K_N-3+c max_j mu(F_j)+floor(log2 N')K_N+3
       <=c mu(G')+(floor(log2 N')+1)K_N.

The one-vertex and no-hanging-component cases follow from the block bound.
For L_N=ceil(log2(N+1))=floor(log2 N)+1, the overhead is at most
L_N(6L_N+C_eta+3)<=(C_eta+9)L_N^2. Take D_eta=C_eta+9. QED.

## 10. Completion for all cubic multigraphs

First let K be connected cubic on D vertices, allowing loops. Deleting its
ell loops leaves a connected loopless subcubic graph with the same cuts
and the same cutwidth. It has 3D/2-ell edges, so its cycle rank is
D/2+1-ell<=D/2+1. By (19), for every eta>0,

    cw(K)<=(gamma+eta)D+2gamma+2eta
                          +D_eta ceil(log2(D+1))^2.        (20)

Thus the positive excess over gamma D is o(D), uniformly over connected
cubic multigraphs.

For each positive integer d let e(d) be the maximum of
max(0,cw(K)-gamma d) over connected cubic multigraphs on d vertices, or zero
if none exists. This is a maximum over a finite set up to relabeling:
degree three bounds all edge and loop multiplicities. Define

    h(D)=max_(1<=d<=D) e(d).

It is nonnegative and nondecreasing. To see h(D)/D->0, fix epsilon>0.
Equation (20), first choosing eta small and then d large, gives
e(d)<=epsilon d for all d>=d_0. Consequently

    h(D)<=max(max_(d<d_0)e(d),epsilon D).

Its ratio has limsup at most epsilon; let epsilon decrease to zero.

For a disconnected cubic multigraph concatenate optimal vertex orderings
of the components. The resulting width is the maximum component width.
A component on d<=D vertices has width at most gamma d+h(d)<=gamma D+h(D).
This proves Theorem 1 with the stated common h. No explicit rate for h is
asserted or needed by the circuit transfer. QED.
