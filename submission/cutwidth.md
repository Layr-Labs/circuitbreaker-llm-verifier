# Cubic cutwidth at most 0.13065 times the number of vertices, plus o(n)

## 0. Statement, provenance, and the finite numerical premise

A multigraph may have loops and parallel edges. A loop contributes two to
degree and never crosses a vertex cut. A graph is cubic when every degree
is three. The width of a vertex ordering is the maximum number of edges
with exactly one endpoint in a proper prefix; cutwidth, denoted cw, is the
minimum width of an ordering. Put

    gamma = 2613/20000 = 0.13065.

**Theorem 1.** There is a nonnegative nondecreasing function h on the
positive integers, with h(D)/D tending to zero, such that every cubic
multigraph K on D vertices satisfies

    cw(K) <= gamma D + h(D).

Disconnected graphs, loops, and parallel edges are included. This is the
hypothesis H(gamma,h) used by the circuit transfer in `circuit.md`.

The new construction optimizes a vertex cover of the entire crossing
graph at each threshold of a Gaussian vertex ordering. A deterministic
interpolation theorem converts these covers into an actual cutwidth
ordering with an absolute additive loss, even when the optimal edge
assignment changes on a large set at one threshold. Exact lifting of
crossing-cover numbers through finite covers permits Gaussian comparison
with graphs of arbitrarily large girth. On the limiting cubic-tree wave,
crossing components have maximum degree two; a finite odd path-count sum
bounds their cover cost. All constructions on the graph itself use fixed
finite-radius Gaussian profiles. No girth assumption or infinite-volume
ordering assumption is made about the input graph.

**Credit.** The circuit transfer, Gaussian ordering approach, and kernel
compression used here build on McGuire, Schlesinger, and Williams. The
Hellinger pushforward comparison for vertex profiles appears in Arthur
Ostapenko's CircuitBreaker submission
5fef109bf7f1b664295244a0dee6df066ed2232d. The global crossing-cover
optimization, exact cover-lifting argument, Hall-matching interpolation,
and finite alternating-path certificate provide the improvement here.
The submodularity of the matching cut function and monotone canonical
crossing covers are established in the maximum matching-width literature,
including Jeong, Sæther, and Telle, *Maximum matching width: new
characterizations and a fast algorithm for dominating set*,
https://arxiv.org/abs/1507.02384. The forms used here, and all graph and
Gaussian comparison arguments needed below, are proved in this file.

The sole numerical premise is finite. Let (X_0,...,X_13) be a centered
Gaussian vector with unit variances and covariance

    Cov(X_i,X_j) = w_|i-j|,
    w_d = (1+d/3) 2^(-d/2),  d >= 0.

For 1<=j<=13, let p_j be the probability that all j consecutive pairs
(X_0,X_1),...,(X_(j-1),X_j) have opposite signs. Define

    Gamma_13 = sum_(j=1)^13 (-1)^(j+1) 3*2^(j-2) p_j.

**Numerical fact C.** Gamma_13 <= 0.130648091104 < gamma.

The accompanying `certificate.md` and its declared interval computation
prove this fact, including quadrature, tail, root-isolation, and arithmetic
errors. The covariance is positive semidefinite by the finite-profile
limit proved below; the numerical certificate also gives an explicit
stationary Gaussian realization. No empirical graph experiment, Monte
Carlo estimate, or infinite alternating-series assumption is a premise.
The rest of this file proves that fact C implies Theorem 1, including a
single uniform sublinear function h for all the stated graphs.

## 1. Symmetric submodular Gaussian comparison

For a function F on {0,1}^q, its mixed difference at distinct coordinates
i,j is F(s^11)-F(s^10)-F(s^01)+F(s^00), with the other coordinates fixed.
It is submodular when these differences are nonpositive, equivalently
F(a)+F(b)>=F(a meet b)+F(a join b). It is symmetric when F(1-a)=F(a).
Write s_t(x)_i=1[x_i<t].

**Lemma 1.1.** Let X,Y be centered Gaussian vectors with unit coordinate
variances, and let F be symmetric submodular.

(a) If Cov(Y_i,Y_j)>=Cov(X_i,X_j) for all i!=j, then
E F(s_t(Y))<=E F(s_t(X)) for every t.

(b) For every t, E F(s_t(X))<=E F(s_0(X)).

**Proof of (a).** The multilinear extension is

    Ftilde(p)=sum_b F(b) product_i p_i^b_i(1-p_i)^(1-b_i).

For i!=j its mixed derivative is an average of mixed differences of F,
so is nonpositive. For delta>0 let sigma_delta(u)=Phi(-u/delta) and set

    f_delta(z)=Ftilde(sigma_delta(z_1-t),...,sigma_delta(z_q-t)).

This is smooth and bounded, with bounded first and second derivatives;
its off-diagonal mixed derivatives are nonpositive. Take X,Y independent
and put Z_r=sqrt(1-r)X+sqrt(r)Y. Gaussian integration by parts gives

    E f(Y)-E f(X)
      = (1/2) integral_0^1 sum_(i,j)
          [Cov(Y_i,Y_j)-Cov(X_i,X_j)] E partial_i partial_j f(Z_r) dr.  (1)

For completeness, express the Gaussian vectors as linear combinations
of independent standard normals. One-dimensional integration by parts
gives E[U_i H(V)]=sum_j Cov(U_i,V_j) E partial_j H(V). Differentiating
E f(Z_r) for 0<r<1 and applying this identity separately to X and Y
gives the integrand in (1). Bounded derivatives justify integration and
extension to the endpoints. Equal variances cancel the diagonal terms;
every other summand is nonpositive. Apply (1) to f_delta and let delta
decrease to zero. Each coordinate has a continuous unit normal marginal,
so bounded convergence proves (a), including singular covariance matrices.

**Proof of (b).** First replace X by sqrt(1-eps^2)X+eps Z with independent
standard coordinate noise Z. It suffices to prove the statement for the
resulting positive definite covariance and then take eps to zero. Writing
D_iF(b)=F(1,b)-F(0,b), differentiation of Gaussian orthant integrals gives

    d/dt E F(s_t(X))
      = phi(t) sum_i E[D_i F(s_t(X_-i)) | X_i=t].            (2)

One obtains (2) by moving one orthant boundary at a time; Gaussian tails
justify differentiation. For t>=0, condition on X_i=t and put U_j=X_j-t
for j!=i. Its conditional mean is m_j=(C_ij-1)t<=0. The vector
U'=-U+2m has the same conditional distribution as U and satisfies U'<=-U
coordinatewise. Thus the bits 1[U'<0] dominate 1-1[U<0] almost surely.
Submodularity makes D_iF nonincreasing in the other bits, and symmetry
gives D_iF(1-b)=-D_iF(b). Its conditional expectation in (2) is therefore
nonpositive. The expectation is nonincreasing for t>=0 and, by symmetry,
even in t. Finally pass to the noise limit by bounded convergence. QED.

**Lemma 1.2.** Partial minimization preserves submodularity: if G(b,x) is
submodular, then F(b)=min_x G(b,x) is submodular. Complement symmetry is
also preserved when all bits of G are complemented together. For fixed b,
the minimizers in x are closed under meet and join and have a least element.

**Proof.** If x,y minimize at b,c, then

    F(b)+F(c)>=G(b meet c,x meet y)+G(b join c,x join y)
             >=F(b meet c)+F(b join c).

Complementation bijects feasible assignments. At a fixed b, the same
inequality and optimality show that meet and join of minimizers are
minimizers. Their finite nonempty set therefore has a least element. QED.

## 2. Crossing covers and monotone optimal assignments

For a graph H and vertex bits b, let D_b be the spanning graph whose edges
join opposite bits, and put F_H(b)=tau(D_b), its minimum vertex-cover
number. The graph D_b is bipartite, with the bits as its two classes. For
edge bits x define

    Q_b(x)=sum_v NAE(b_v,(x_e:e incident with v)),
    C(x)=sum_v NAE((x_e:e incident with v)),

where NAE is one unless all its arguments agree. An empty or singleton
list has NAE zero. On cubic graphs, C(x) counts the vertices split by the
edge cut x.

**Lemma 2.1.** F_H(b)=min_x Q_b(x). It is symmetric submodular, and
changing one coordinate of b changes F_H by at most one.

**Proof.** A vertex with Q term zero forces all incident edge bits to its
own bit. Thus the penalized vertices cover every opposite-bit edge.
Conversely, given a cover, assign every edge incident with an uncovered
vertex that vertex's bit. Two such requirements agree, since an edge
between uncovered vertices cannot cross. Fill the other edges arbitrarily.
At most the cover vertices are penalized, proving the identity.

NAE of any arity is symmetric submodular. At two specified coordinates
its mixed difference is -1 when all remaining coordinates are zero or all
are one, and zero otherwise; for arity two it is -2. Sums and coordinate
identifications preserve submodularity. Apply Lemma 1.2. Changing one b_v
changes Q for any fixed x at only v, by at most one. Evaluate the old
optimal assignment at the new bits and reverse the argument. QED.

Call x clamped at b when, for every edge uv,

    b_u AND b_v <= x_uv <= b_u OR b_v.                       (3)

Thus every same-bit edge has that bit. Setting a same-bit edge to its
endpoint bit cannot increase Q: at either endpoint it cannot destroy an
all-equal list. Hence the clamped and unrestricted minima both equal F_H(b).

**Lemma 2.2.** The least clamped minimizer x(b) exists and is nondecreasing
in b.

**Proof.** At fixed b the clamped cube is closed under meet and join, so
the proof of Lemma 1.2 supplies its least minimizer. Suppose b<=c, and let
x,y be the least clamped minimizers at these anchors. The lower and upper
bounds in (3) are nondecreasing in b. Therefore x meet y is feasible at b
and x join y at c. Joint submodularity of Q and endpoint optimality force
both of these assignments to be optimal at the corresponding anchors.
Leastness gives x<=x meet y, hence x<=y. Only comparable anchors are used;
no claim about a joint feasible sublattice for incomparable anchors is
needed. QED.

## 3. A deterministic conversion from crossing covers to cutwidth

The following interpolation is the step allowing global optimization.
In this section H is loopless with maximum degree at most three.

**Lemma 3.1.** Fix b. Suppose x<=y are clamped at b and y is a global
minimizer of Q_b in the clamped cube. The edges on which x and y differ
can be changed one at a time from zero to one so that every intermediate
assignment z satisfies

    C(z) <= Q_b(x)+2.                                      (4)

**Proof.** Let U be the differing edges. Each joins opposite anchors.
Classify vertices incident with U by their Q terms at x and y: A has
values (0,1), B has (1,0), and P has (1,1). The case (0,0) is impossible
because its incident U bits change from all zero to all one while its
anchor stays fixed. Vertices in A have anchor zero and every fixed edge
zero; vertices in B have anchor one and every fixed edge one. A P vertex
has a disagreeing pair among its anchor and fixed edges, so its Q term
is permanently one. No U edge joins two A vertices or two B vertices.
If q is the total Q contribution outside the vertices incident with U,

    Q_b(x)=q+|P|+|B|.                                      (5)

The U edges between A and B have a matching covering A. Indeed, take
S subset A and start from y. Set to zero all U edges incident with S.
This remains clamped, reduces the Q terms at S by |S|, and can increase
terms only at B neighbors of S, by at most |N_B(S)|. P terms stay one,
and there are no A-A edges. Optimality of y implies |N_B(S)|>=|S|.
Hall's condition gives the matching. One may prove the version of Hall
used here by taking a maximum matching: an alternating search from an
unmatched A vertex would either augment it or exhibit a set S with fewer
neighbors than vertices, contradicting the condition.

First flip all U edges inside P and between B and P. A vertices are
unaffected, so the split cost is at most q+|P|+|B|. Next, for each A
vertex incident with three U edges, flip those three consecutively.
Such a vertex has no fixed incident edge. It has actual split cost zero
before and after its block, and at most one during it. Completed and
unstarted such A vertices cost zero, all B vertices cost at most |B|,
and P at most |P|. Thus this phase costs at most q+|P|+|B|+1.

Let A' be the remaining A vertices. Each has U-degree at most two and
at least one B neighbor, by the matching. Each A' vertex with two
distinct B neighbors represents an ordinary edge in a multigraph on B;
each with one B neighbor represents a half-edge there. This includes
the case of two parallel U edges to that same B vertex: both are one
A' block and give one half-edge. An A'-P edge, if present, will be
flipped within its A' block.

The matching restricted to A' implies that in every component of this
multigraph the number of ordinary edges plus half-edges is at most the
number of B vertices. Every component is consequently a tree with at
most one half-edge, or a unicyclic graph with no half-edge; isolated B
vertices are allowed. Process one component at a time. In a tree peel
ordinary leaf edges, taking leaves away from the half-edge if present,
and process that half-edge last. In a unicyclic component peel attached
trees and then process cycle edges consecutively around the cycle.

Call B finished when all its U edges have been flipped. Every B-P edge
was flipped in the first phase. At completed A' blocks, the number of
processed blocks minus the number of finished B vertices is at most one.
Each peeled tree edge finishes at least one B vertex; a last half-edge
finishes its vertex; only the first edge of a cycle can fail to finish
a vertex, and the last finishes two. Completed components have a
nonpositive final difference. B vertices finished in earlier phases
only improve the bound. Within a block, counting that block as started
can raise the difference by at most one, so it is at most two.

Unstarted A' vertices have split cost zero; started ones cost at most
one. Finished B vertices cost zero, and unfinished ones at most one.
The previous A vertices of U-degree three cost zero. Therefore

    C(z) <= q+|P|+(number of started A' blocks)+(number of unfinished B)
         <= q+|P|+|B|+2 = Q_b(x)+2.

All U edges have now been flipped, proving (4). QED.

**Lemma 3.2.** For any vertex ordering of a loopless cubic graph H, let M
be the maximum of F_H over its prefix bits. Then

    cw(H) <= M+6.                                         (6)

**Proof.** Consider consecutive prefix vectors b<=c differing at one
vertex, and their least clamped minimizers x<=y from Lemma 2.2. Clamp x
at c by flipping to one the edges that have just become same-anchor-one
edges. There are at most three, all incident with the changed vertex.
Call the result z; then x<=z<=y and z is clamped at c. Changing the
anchor affects Q at one vertex by at most one. Each clamping edge flip
is individually nonincreasing for Q_c. Hence throughout this first phase

    C <= Q_c <= Q_c(x) <= Q_b(x)+1 = F_H(b)+1.

Apply Lemma 3.1 to z<=y at c. Since Q_c(z)<=F_H(b)+1, the remaining
flips have C<=F_H(b)+3. The initial and final all-equal vertex vectors
have respectively all-zero and all-one clamped edge fields. Concatenating
these monotone refinements therefore orders every edge exactly once,
with every edge-prefix split cost at most M+3.

Number the edges in that order. Put

    B_i={v: first(v)<=i<=last(v)},

where first and last are the extreme incident edge positions. Relative
to the split cost immediately after edge i, B_i can add only its at most
two endpoints. Thus |B_i|<=M+5. Order vertices by their median incident
edge positions, breaking ties arbitrarily. For a proper prefix let theta
be the median position of its last vertex. For a crossing edge uv with u
inside the prefix, an edge position below theta is the unique lowest
incident edge at v, and can be charged to v in B_theta. A position above
theta is the unique highest incident edge at u and can be charged to u
in B_theta. No vertex is charged twice. At most one edge has position
exactly theta, so the cut has size at most |B_theta|+1<=M+6. QED.

By the bipartite matching-cover equality, (6) also says that cubic
cutwidth is at most linear maximum matching width plus six. The proof
uses no bound on the number of edges changing at a single anchor step.

## 4. Finite covers and Hellinger Gaussian profiles

For now H is a connected finite simple cubic graph on n vertices.
Its universal cover is the cubic tree T: choose a graph vertex v, take
nonbacktracking walks from v as tree vertices (including the empty walk),
and join a walk to its one-edge extensions. Map each walk to its endpoint.
The empty walk has three extensions, and every other walk has a parent
and two extensions. The map is a bijection on incident directed edges
at each tree vertex, so every graph walk lifts uniquely from a chosen
initial preimage.

Any two tree vertices over the same graph vertex are related by a
fibre-preserving tree automorphism: project the path from the first
vertex to any other tree vertex and lift it from the second. The reverse
construction is its inverse. This proves the lift independence of the
profiles defined below. The same tree covers every connected finite
cover of H, compatibly with its projection to H.

**Lemma 4.1 (finite covers of large girth).** For every integer g, H has
a finite regular cover H'->H of girth greater than g.

**Proof.** Choose a spanning tree and orient the remaining edges; assign
one free generator to each oriented non-tree edge and its inverse to the
reverse orientation. A closed nonbacktracking walk records a nonempty
reduced word by deleting its tree edges. The word is nonempty because a
tree has no nonbacktracking closed walk. Consecutive inverse letters
would be inverse traversals of one non-tree edge separated by a closed
walk in the spanning tree. A nonempty such tree walk has an immediate
backtrack, and an empty one makes the two traversals immediate inverses.
Both contradict nonbacktracking. Thus a walk of length at most g records
a nontrivial reduced word of length at most g.

Each nontrivial reduced word s_1...s_k has a finite permutation
representation in which it is not the identity. On points 0,...,k,
prescribe for each generator a the maps i-1->i when s_i=a, and i->i-1
when s_i=a^(-1). These prescriptions are partial injections: a conflicting
domain or image could only come from adjacent inverse letters. Complete
each partial injection to a permutation. The word sends 0 to k and is
not the identity. Take the product of these representations over all
nontrivial reduced words of lengths at most g, and let Q be the finite
group generated by the resulting generator images. Every such word is
nonidentity in Q.

The cover has vertices (u,q) with u in V(H), q in Q. A spanning-tree
edge joins equal q coordinates; an oriented non-tree edge carrying a
joins (u,q) to (v,qa). This is a finite cover of degree |Q|. A closed
nonbacktracking walk of length at most g upstairs would project to one
downstairs and multiply q by its nonidentity recorded word. The regular
action of Q has no fixed point for a nonidentity element, a contradiction.
The covering graph has no such walk, so has girth greater than g. QED.

**Lemma 4.2 (exact cover lifting).** For any finite degree-m cover
pi:H'->H and repeated bits b'=b composed with pi,

    F_H'(b')=m F_H(b).                                     (7)

**Proof.** Pulling back a minimum crossing cover gives the upper bound.
For any crossing cover C' upstairs put

    y_v=|C' intersect pi^(-1)(v)|/m.

There are exactly m lifted copies of each crossing edge uv, so
y_u+y_v>=1. Thus y is a fractional vertex cover downstairs. A maximum
matching in a bipartite graph has the same size as a minimum cover:
given a maximum matching, mark the vertices reachable by alternating
paths from unmatched left vertices; the unmarked left vertices together
with the marked right vertices form a cover with one vertex per matched
edge. An uncovered edge would extend the search, and a reachable unmatched
right vertex would give an augmenting path. This proves the equality.
Summing y over the disjoint endpoints of such a maximum matching shows

    F_H(b) <= sum_v y_v = |C'|/m.

Apply this to a minimum cover upstairs for the lower bound in (7). QED.

For an integer R>=1 put A_R=1+3R/2 and use the nonnegative unit profiles
on T

    a_u(z)=A_R^(-1/2) 2^(-dist_T(u,z)/2) 1[dist_T(u,z)<=R].

They have squared norm one: every nonzero radius-j sphere has
3*2^(j-1) vertices and contributes 3/2 before normalization. Write
r_d(R)=<a_u,a_v> at tree distance d. Push squared profiles to a finite
cubic graph H by

    P_v^H(x)=sum_(z:pi(z)=x) a_u(z)^2,  pi(u)=v.

The fibre automorphisms make this independent of u. It is a probability
measure supported within graph distance R of v. Let (g_x)_x and (xi_v)_v
be independent standard normals and, for 0<eps<1, set

    Z_v^H=sqrt(1-eps^2) sum_x sqrt(P_v^H(x))g_x+eps xi_v.    (8)

Each key has variance one, and distinct keys are distinct almost surely:
their difference has variance at least 2eps^2. On T use independent
vertex noises in the corresponding finite sums, plus independent
coordinate noises of variance eps^2.

For a finite cover pi:H'->H the measures satisfy

    P_pi(v')^H(x)=sum_(y:pi(y)=x) P_v'^H'(y).

Fibrewise Cauchy--Schwarz, followed by summation, shows that the covariance
of two repeated base-graph keys Z_pi(u')^H,Z_pi(v')^H is at least the
covariance of the corresponding H' keys. For distinct images, the
independent coordinate noises contribute zero in both comparisons.
For equal images, repetition adds eps^2 on the base side whereas
distinct cover coordinates have independent noises, which only
strengthens the inequality. Variances are one.

Apply Lemma 1.1 to the symmetric submodular F_H', and then (7). For every
t and every degree-m cover,

    m E F_H(s_t(Z^H))
       =E F_H'(s_t(Z^H repeated))
       <=E F_H'(s_0(Z^H')).                                (9)

Repeated coordinates and singular covariances cause no problem for the
comparison. The cover is only an expectation device; all random inputs
and eventual orderings on H use (8) on H itself.

## 5. Finite path counts and the critical tree wave

For a graph D of maximum degree three, let N_j(D) be half the number of
directed nonbacktracking walks of j edges. For any fixed odd J put

    A_J(D)=sum_(j=1)^J (-1)^(j+1) N_j(D),
    C_J=1+sum_(j=1)^J (j+1)3*2^(j-2)=1+3J*2^(J-1).

**Lemma 5.1.** If B is the set of vertices of D having degree three, then

    tau(D) <= A_J(D)+C_J |B|.                              (10)

**Proof.** First suppose D has maximum degree two. On a path with r
edges, N_j=max(r-j+1,0). Its odd alternating sum equals ceil(r/2) for
r<J, and r-(J-1)/2 for r>=J; both are at least its minimum cover
ceil(r/2). On a cycle, N_j is its number of vertices for every j,
including walks that wrap around it. Its odd alternating sum equals
that number, at least the cover size. Isolated vertices contribute zero.
Summing over components proves tau(D)<=A_J(D) in this case.

In general D-B has maximum degree two, and adding B to a cover of D-B
covers D. A directed nonbacktracking j-walk containing a specified
vertex in a specified position has at most 3*2^(j-1) possibilities:
choose one incident edge in at most three ways and each subsequent
continuation in at most two; the same count applies when the specified
position is internal, extending in both directions. Consequently

    |N_j(D)-N_j(D-B)| <= (j+1)3*2^(j-2)|B|.

Using tau(D)<=|B|+tau(D-B)<=|B|+A_J(D-B), and bounding each signed
difference in absolute value, gives (10). QED.

**Lemma 5.2.** For each fixed d>=0,

    r_d(R) -> w_d=(1+d/3)2^(-d/2) as R tends to infinity.   (11)

**Proof.** The d=0 case is exact. For d>=1 use the path u=q_0,...,q_d=v.
A vertex attached to q_i at depth j>=1 contributes the unnormalized
product 2^(-d/2-j). At each endpoint there are 2^j vertices at that
depth, and at each internal path vertex 2^(j-1). The permitted depths,
requiring both profile distances at most R, differ from 1,...,R by
O(d) terms. The two endpoint contributions are
(2R+O(d))2^(-d/2); the d-1 internal contributions are
((d-1)R/2+O(d^2))2^(-d/2). Path vertices contribute O(d)2^(-d/2).
Divide by 1+3R/2 to obtain (11). QED.

Thus every finite collection of tree keys in (8) converges, as R tends
to infinity and eps to zero, to a centered unit-variance Gaussian vector
with correlations w_distance. These limiting matrices are positive
semidefinite, being limits of finite covariance matrices. Finite Gaussian
convergence applies to sign events because their boundaries are contained
in coordinate-zero hyperplanes, each of probability zero. No nonsingularity
of a limiting joint covariance is required.

For a root v and its three neighbors, the limiting Gaussian vector obeys

    sum_(u adjacent v) X_u = 2sqrt(2) X_v                  (12)

almost surely. Indeed w_1=2sqrt(2)/3 and w_2=5/6, so the variance of the
difference is 3+6w_2+8-12sqrt(2)w_1=0. If X_v is positive, all three
neighbors cannot be negative; if X_v is negative, they cannot all be
positive. Therefore the limiting probability of a degree-three crossing
star at threshold zero is zero.

Let beta_(R,eps) be that probability for the finite tree profiles. Let
p_j^(R,eps) be the probability of alternating signs on a fixed j-edge
tree path, and define

    Gamma_J^(R,eps)=sum_(j=1)^J (-1)^(j+1)3*2^(j-2)p_j^(R,eps).

By the preceding finite-dimensional convergence, for the fixed J=13,

    beta_(R,eps) -> 0,
    Gamma_13^(R,eps) -> Gamma_13.                          (13)

These are limits of finitely many probabilities only. The coefficient
C_13=159745 is finite and is fixed before taking either limit.

## 6. Expectations and concentration at all thresholds

Fix R,eps, and J=13. Choose, by Lemma 4.1, a degree-m cover H' of H
with girth greater than 4R+2J+10. Every neighborhood used by profiles
on a path of length at most J or a root star lifts injectively to the
cubic tree. To see this, its union lies in a ball of radius R+J+1;
a collision or an extra edge in this ball would produce a cycle of
length at most 2(R+J+1)+1, below the chosen girth. The pushed profiles
therefore have exactly the tree covariances needed for these events.

The graph H' has n'=mn vertices and exactly n'3*2^(j-1) directed
nonbacktracking j-walks. Its girth ensures these walks have distinct
vertices for j<=J. For its zero-threshold crossing graph D, linearity
of expectation gives

    E A_J(D)=n' Gamma_J^(R,eps),
    E |B|=n' beta_(R,eps).

Use (10), divide (9) by m, and cancel n'=mn. For every input graph H
and every real threshold t,

    E F_H(s_t(Z^H)) <= a_(R,eps) n,
    a_(R,eps)=Gamma_13^(R,eps)+159745 beta_(R,eps).          (14)

There is no loss depending on the cover degree or on short cycles of H.
The left side of (14) is nonnegative, so a_(R,eps)>=0. By (13), this
quantity tends to Gamma_13.

We next obtain one realization good at every threshold, without any
locality assumption on its optimal cover. Write F(t)=F_H(s_t(Z^H)).
A radius-R ball has at most B_R=3*2^R-2 vertices. Changing one independent
input g_x in (8) can change keys only in this ball around x. Each changed
anchor changes F by at most one, by Lemma 2.1; thus that input changes
F(t) by at most B_R, regardless of its new value. Changing one xi_v
changes it by at most one.

For independent inputs and a function whose change upon replacing input
i is at most c_i, its variance is at most sum_i c_i^2. Indeed, expose
the inputs in order and subtract successive conditional expectations.
The resulting martingale differences are orthogonal and bounded in
absolute value by c_i; sum their second moments. Applying this here gives

    Var F(t) <= n(B_R^2+1)=K_R n.                         (15)

For any interval I let W(I) count graph vertices whose keys lie in I.
The same input-replacement bounds apply to W, so Var W(I)<=K_R n.
Every key is standard normal, hence

    E W([a,b)) <= n(b-a)/sqrt(2pi),
    E W((-infinity,-T)), E W([T,infinity)) <= n Phi(-T).    (16)

Changing the anchors one at a time also gives, for s<t,

    |F(t)-F(s)| <= W([s,t)).                               (17)

The all-zero and all-one anchor vectors have cost zero, so F in either
tail is at most its corresponding tail count.

Fix delta>0. Choose T with Phi(-T)<delta/8 and a finite mesh
-T=t_0<...<t_k=T of gaps at most delta sqrt(2pi)/8. Apply Chebyshev's
inequality and the union bound to the finitely many F(t_i), interval
counts W([t_i,t_(i+1))), and two tail counts. Equations (15)-(16) imply
that, for all sufficiently large n depending only on R,delta, there is
a realization in which each deviates from its expectation by less than
delta n/8. Intersect with pairwise distinctness of keys, a probability-one
event because eps>0.

For a threshold in a mesh interval, (14),(16),(17) bound its cost by
(a_(R,eps)+3delta/8)n. Either tail has cost below delta n/4. Because
a_(R,eps)>=0, this proves the convenient uniform bound

    sup_t F(t) <= (a_(R,eps)+delta)n.                      (18)

All constants and the mesh were fixed before the input graph size.
Distinct keys turn all vertex prefixes of their increasing order into
threshold vectors; Lemma 3.2 then supplies an actual vertex ordering
with cutwidth at most the right side of (18) plus six.

## 7. The uniform finite-parameter conclusion

**Theorem 7.1.** For every eta>0 there is N_eta such that every connected
simple cubic graph H on n>=N_eta vertices satisfies

    cw(H) <= (gamma+eta)n.                                (19)

**Proof.** First fix J=13. By numerical fact C and (13), choose a finite
R and a positive eps<1 such that 0<=a_(R,eps)<=gamma+eta/3. Apply (18)
with delta=eta/3. Finally take n sufficiently large for its concentration
statement and for 6/n<=eta/3. Lemma 3.2 proves (19). Each parameter and
size threshold depends only on eta, not on H. No simultaneous infinite
radius and graph-size limit is being used. QED.

The following sections extend this conclusion uniformly to loops,
parallel edges, bridges, and disconnected graphs. They also construct
the required nondecreasing function h explicitly from finite maxima.

## 8. Kernel compression

A kernel is a finite simple triangle-free 2-connected cubic graph.
Theorem 7.1, with finitely many
smaller sizes covered by the trivial bound cw<=3n/2, implies: for every
eta>0 there is C_eta>=0 such that all kernels satisfy

    cw(H)<=(gamma+eta)|V(H)|+C_eta.                         (20)

For a connected loopless multigraph B with N vertices and M edges define
its cycle rank nu(B)=M-N+1.

**Lemma 8.1.** If B is connected, bridgeless, loopless and subcubic on N
vertices, either cw(B)<=6ceil(log2 N), or there is a kernel H on
h<=2(nu(B)-1) vertices such that any width-q ordering of H expands to a
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
triangle-free. Then h/2=M_H-N_H<=M_B-N_B=nu(B)-1. It is 2-connected: if
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

    nu(G)=nu(B)+sum_j nu(F_j).

If B has width w and the F_j have widths w_j, then G has width at most
w+max_j w_j+3, with an empty maximum interpreted as zero.

**Proof.** An edge between vertices of B belongs to B: if it were a
bridge, a path within B would complete a cycle through it. Every edge
leaving B is a bridge. There is at least one from each F_j by connectedness;
two would give a cycle through either one using paths within B and F_j.
Thus it is unique. Every edge lies in B, in an F_j, or is an e_j. Counting
vertices and edges proves the cycle-rank identity. B itself is bridgeless, because
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

    cw(G)<=(2gamma+2eta)nu(G)+D_eta ceil(log2(N+1))^2.        (21)

**Proof.** Put c=2gamma+2eta and K_N=6ceil(log2 N)+C_eta+3, where C_eta
comes from (20). Lemma 8.1 shows that any bridgeless block B of size at
most N has width at most c nu(B)+K_N-3: either it already has logarithmic
width, or expand a kernel of h<=2(nu(B)-1) vertices using (20).

Induct on N'<=N to prove every connected subcubic G' of size N' has width
at most c nu(G')+(floor(log2 N')+1)K_N. Choose its centroid block B. Each
hanging component has size at most N'/2, so its induction overhead is at
most floor(log2 N')K_N. Lemma 9.1 bounds the new width by

    c nu(B)+K_N-3+c max_j nu(F_j)+floor(log2 N')K_N+3
       <=c nu(G')+(floor(log2 N')+1)K_N.

The one-vertex and no-hanging-component cases follow from the block bound.
For L_N=ceil(log2(N+1))=floor(log2 N)+1, the overhead is at most
L_N(6L_N+C_eta+3)<=(C_eta+9)L_N^2. Take D_eta=C_eta+9. QED.

## 10. Completion for all cubic multigraphs

First let K be connected cubic on D vertices, allowing loops. Deleting its
ell loops leaves a connected loopless subcubic graph with the same cuts
and the same cutwidth. It has 3D/2-ell edges, so its cycle rank is
D/2+1-ell<=D/2+1. By (21), for every eta>0,

    cw(K)<=(gamma+eta)D+2gamma+2eta
                          +D_eta ceil(log2(D+1))^2.        (22)

Thus the positive excess over gamma D is o(D), uniformly over connected
cubic multigraphs.

For each positive integer d let e(d) be the maximum of
max(0,cw(K)-gamma d) over connected cubic multigraphs on d vertices, or zero
if none exists. This is a maximum over a finite set up to relabeling:
degree three bounds all edge and loop multiplicities. Define

    h(D)=max_(1<=d<=D) e(d).

It is nonnegative and nondecreasing. To see h(D)/D->0, fix epsilon>0.
Equation (22), first choosing eta small and then d large, gives
e(d)<=epsilon d for all d>=d_0. Consequently

    h(D)<=max(max_(d<d_0)e(d),epsilon D).

Its ratio has limsup at most epsilon; let epsilon decrease to zero.

For a disconnected cubic multigraph concatenate optimal vertex orderings
of the components. The resulting width is the maximum component width.
A component on d<=D vertices has width at most gamma d+h(d)<=gamma D+h(D).
This proves Theorem 1 with the stated common h. No explicit rate for h is
asserted or needed by the circuit transfer. QED.
