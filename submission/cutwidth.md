# Cubic cutwidth with a quantified structural improvement

Theorem 1 and Sections 0--10 establish the retained finite-path baseline.
Theorem 2 in Section 11 proves the stronger graph constant
130647566423167663341749381/1000000000000000000000000000 by an additional
uniform structural gain of 10^-18. Theorem 2 supplies the graph hypothesis
used in the current circuit theorem. The baseline symbol gamma below is
not silently replaced by the stronger constant.

## 0. Baseline statement, provenance, and finite numerical premise

A multigraph may have loops and parallel edges. A loop contributes two to
degree and never crosses a vertex cut. A graph is cubic when every degree
is three. The width of a vertex ordering is the maximum number of edges
with exactly one endpoint in a proper prefix; cutwidth, denoted cw, is the
minimum width of an ordering. Put

    gamma = 130647566423167664341749381/1000000000000000000000000000 = 0.130647566423167664341749381.

**Theorem 1.** There is a nonnegative nondecreasing function h on the
positive integers, with h(D)/D tending to zero, such that every cubic
multigraph K on D vertices satisfies

    cw(K) <= gamma D + h(D).

Disconnected graphs, loops, and parallel edges are included. This has the
form H(gamma,h) of `circuit.md`; the current circuit theorem uses the
stronger instance proved as Theorem 2 in Section 11.

The inherited construction optimizes a vertex cover of the entire crossing
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
and original finite alternating-path certificate are retained from
submission fcf21648-1b1f-463c-b5bc-c50f53f63436, source commit
81470c9513639af985c88699dbcd2e2efbef7463. The baseline graph argument and the localized 65-term certificate are
retained from submission ae372f38-09d4-452a-af05-4954cdd2f168, commit
2a8b03fe0be18fdf4918dc5fc34da51534f1aca1, with earlier sources credited
in proof.md. Its finite graph constant is called gamma_base in Section 11.
The new contribution is that section's buffered-curvature and majority-gain
argument, supported by the complete local Gaussian probability proofs in
certificate.md Sections 9 and 10. All inherited arguments remain supplied.
The submodularity of the matching cut function and monotone canonical
crossing covers are established in the maximum matching-width literature,
including Jeong, Sæther, and Telle, *Maximum matching width: new
characterizations and a fast algorithm for dominating set*,
https://arxiv.org/abs/1507.02384. The forms used here, and all graph and
Gaussian comparison arguments needed below, are proved in this file.

The baseline numerical premise is finite. Let (X_0,...,X_65) be a centered
Gaussian vector with unit variances and covariance

    Cov(X_i,X_j) = w_|i-j|,
    w_d = (1+d/3) 2^(-d/2),  d >= 0.

For 1<=j<=65, let p_j be the probability that all j consecutive pairs
(X_0,X_1),...,(X_(j-1),X_j) have opposite signs. Define

    Gamma_65 = sum_(j=1)^65 (-1)^(j+1) 3*2^(j-2) p_j.

**Numerical fact C.** Gamma_65 <= 0.1306475664231676643417493807 < gamma.

The accompanying `certificate.md` and its declared interval computation
prove this fact, including quadrature, tail, root-isolation, and arithmetic
errors. The covariance is positive semidefinite by the finite-profile
limit proved below; the numerical certificate also gives an explicit
stationary Gaussian realization. No empirical graph experiment, Monte
Carlo estimate, or infinite alternating-series assumption is a premise.
Sections 1--10 prove that fact C implies Theorem 1, including a single
uniform sublinear function h for all the stated graphs. Section 11 then
proves Theorem 2 using its additional explicit local probability bounds.

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

By the preceding finite-dimensional convergence, for the fixed J=65,

    beta_(R,eps) -> 0,
    Gamma_65^(R,eps) -> Gamma_65.                          (13)

These are limits of finitely many probabilities only. The coefficient
C_65=3597115094373362565121 is finite and is fixed before taking either limit.

## 6. Expectations and concentration at all thresholds

Fix R,eps, and J=65. Choose, by Lemma 4.1, a degree-m cover H' of H
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
    a_(R,eps)=Gamma_65^(R,eps)+3597115094373362565121 beta_(R,eps).          (14)

There is no loss depending on the cover degree or on short cycles of H.
The left side of (14) is nonnegative, so a_(R,eps)>=0. By (13), this
quantity tends to Gamma_65.

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

**Proof.** First fix J=65. By numerical fact C and (13), choose a finite
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

## 11. A quantified structural improvement

The following theorem strengthens Theorem 1. The symbol gamma in
Sections 0--10 retains its original value; write it here as gamma_base.
Put

    delta = 1/1000000000000000000,
    gamma_new = gamma_base-delta
              = 130647566423167663341749381/1000000000000000000000000000.

**Theorem 2.** There is a nonnegative nondecreasing function h_new on
the positive integers, with h_new(D)/D tending to zero, such that every
cubic multigraph K on D vertices satisfies

    cw(K) <= gamma_new D+h_new(D).

Thus this theorem supplies H(gamma_new,h_new). Its additional numerical
premises are the two finite-window probability bounds proved in
`certificate.md` Sections 9 and 10. All uses of those bounds, and the
graph argument giving the improvement, are supplied below.

### 11.1. Conventions and the two local events

Until the final extension let H be connected, simple and cubic, with n
vertices. Retain F_H, Q_b, C, the least clamped minimizers x(b), and the
finite-radius field Z^H from Sections 2--6. In this section always use

    b_t^+(v)=1[Z_v>t].

Lowering t gives increasing anchor prefixes in decreasing-key order.
Complement symmetry gives F_H(b_t^+)=F_H(s_t(Z)) almost surely at each
fixed t. Thus the comparison and threshold statements about F remain
valid. The gain count defined below is not assumed complement symmetric.

Fix the following constants before any radius, noise, cover or graph:

    q=3/100000,       theta=1/100,       zeta_c=1/700000,
    p=1/10^17,        zeta_g=1/35000,
    eta=99 zeta_c/100=99/70000000,
    tau=zeta_g/21=1/735000,   alpha=3 tau/26=1/6370000.

All critical-window covariances below are w_distance from Section 5.

For the **curvature event**, take adjacent roots u,v in the cubic tree.
Each has two other neighbors A, one designated positive and the other
negative. Each of the four A vertices has children B,C: prescribe B
opposite to A and C equal to A. Prescribe the two children of each B
equal to B. This is a window of 22 vertices: two roots and twenty
prescribed nonroots. Its diameter is at most seven. Call the prescribed
nonroot pattern E_uv. Certificate Section 9 proves: if a centered unit
Gaussian vector on this window has covariance within zeta_c entrywise
of the critical covariance, then its root pair is nonsingular and, for
every |t|<=theta,

    P(E_uv relative to t | X_u=X_v=t) > q/2.                (23)

Here positive means strictly above t and negative strictly below t.
The full covariance may be singular.

On a graph in which this window is embedded with all indicated neighbors
of u,v,A,B, this event forces

    Delta_uv F_H=F_H(11)-F_H(10)-F_H(01)+F_H(00)=-2.        (24)

Indeed each A-B edge crosses, while A-C and the other two edges of B
do not. The crossing components involving u,v,A,B are therefore isolated
from all outside anchors in each of the four root settings. Equal roots
give two two-edge paths and two isolated A-B edges, of total cover
number four. Opposite roots give one five-edge path and two isolated
A-B edges, of total cover number five. The rest of the crossing graph
is unchanged. The mixed difference is 4-5-5+4=-2.

For the **gain event**, orient a tree edge as (w,v), with center w.
The other neighbors of v are a,b and those of w are c,d. Include the
two further neighbors of each of a,b,c,d. This window W(w,v) has
fourteen vertices, thirteen edges, diameter five, and contains the tree
path between every pair of its vertices. Prescribe v,w,d and the two
children of d positive; prescribe a,b,c and their six children negative.
Call this pattern P_w. Certificate Section 10 proves that for every
centered unit Gaussian vector on this window whose covariance is within
zeta_g entrywise of critical,

    P(P_w relative to t) >= p/2  for every |t|<=theta.       (25)

Positive definiteness is unnecessary. For graph applications it will
suffice that these fourteen vertices have distinct images in H. All
neighbors of the six nonbuffer vertices are listed; extra edges among
the eight buffer vertices or from them to the outside are allowed.

### 11.2. Strict covariance and threshold comparison

For a positive definite unit Gaussian covariance C write
f_C(t)=E F_H(b_t^+(X)). If D is symmetric with zero diagonal, Gaussian
differentiation along a positive definite covariance path gives

    d/ds f_(C+sD)(t)
      =sum_(i<j) D_ij phi_(C+sD),ij(t,t)
           E[Delta_ij F_H | X_i=X_j=t].                    (26)

Here phi_ij is the bivariate marginal density; the other bits in the
conditional expectation are 1[X_k>t]. To justify (26), smooth each
indicator and apply the integration-by-parts identity (1). The diagonal
terms vanish; the two ordered terms for each distinct pair cancel the
factor 1/2. The two indicator derivatives converge to masses at t, and
the resulting mixed difference is Delta_ij F_H. Positive definiteness
on a compact subinterval supplies smooth conditional Gaussian densities;
finite sums, bounded mixed differences and Gaussian domination justify
the limit. Every mixed difference is nonpositive by Lemma 2.1.

Let H'->H be a degree-m cover of the girth used in Section 6, and put
n'=mn. Let X be its field (8), and let Y be Z^H repeated on the fibers.
Write c_R=(1-eps^2)r_1(R) and

    delta_e=Cov(Z_u^H,Z_v^H)-c_R >=0  (e=uv),
    Delta=sum_(e in E(H)) delta_e.

The comparison of Section 4 gives C_Y>=C_X off diagonal. All these
covariances are nonnegative, since every profile coefficient is
nonnegative; unit variances give

    0 <= (C_Y-C_X)_ij <= 1.

Require R>=10, initial critical covariance error at most zeta_c/100
on the curvature window, and

    1-c_R >= b_0:=(1-w_1)/2.

Finite-profile convergence permits these conditions. The girth ensures
that each initial window has exactly its finite-tree covariances. Along
C_s=(1-s)C_X+sC_Y, 0<=s<=eta, each window error is at most
zeta_c/100+eta=zeta_c. The covariance is positive definite for s<1
because X has independent positive coordinate noise. By (23), (24) and
submodularity, every adjacent pair has conditional expected mixed
difference at most -q when |t|<=theta.

For a unit Gaussian pair with nonnegative correlation c<1,

    phi_c(t,t)=exp(-t^2/(1+c))/(2 pi sqrt(1-c^2))
               >= d_0:=exp(-theta^2)/(2 pi),  |t|<=theta.

Drop the other, nonpositive, summands of (26) and integrate on [0,eta].
Monotonicity handles the rest of the interpolation. Each base edge has
m lifts, and exact cover lifting (7) identifies the expectation at Y.
The result, using f_X(t)<=f_X(0)<=a_(R,eps)n' from Sections 1 and 6, is

    E F_H(b_t^+(Z^H)) <= a_(R,eps)n-A Delta,
    A=eta q d_0,                         |t|<=theta.       (27)

The covariance at Y can be singular. The strict estimate was obtained
on a fixed initial interval; extension to s=1 uses the non-strict
comparison and Gaussian continuity, since each coordinate has variance
one and hence zero probability of equaling the fixed threshold.

An additional identity supplies loss at the remaining thresholds:

    f_C''(t)+t f_C'(t)
      =2 sum_(i<j)(1-C_ij) phi_C,ij(t,t)
           E[Delta_ij F_H | X_i=X_j=t].                    (28)

To derive it, add an independent common normal of variance a>=0 and
divide the vector by sqrt(1+a). Its covariance is
(C+a 11^T)/(1+a), whose derivative at zero is 1-C off diagonal.
Conditional on that common standard normal G, its expected value is
f_C(t sqrt(1+a)-sqrt(a)G). Differentiation at zero gives
(t f_C'(t)+f_C''(t))/2. Equation (26) gives the other derivative and
proves (28). Positive definiteness justifies both differentiations.

Apply (28) to X on H'. The 3n'/2 adjacent pairs each contribute at
most -2 b_0 d_0 q; every other contribution is nonpositive. Consequently

    f_X''(t)+t f_X'(t) <= -B n',
    B=3q b_0 d_0,                         |t|<=theta.

Symmetry gives f_X'(0)=0. Multiplication by exp(t^2/2), followed by
two integrations, yields

    f_X(t) <= f_X(0)-(B/2)exp(-theta^2/2)n' t^2,
                                                   0<=t<=theta.

The evenness and positive-threshold monotonicity from Lemma 1.1 extend
the loss at theta to all |t|>=theta. Comparison and cover lifting give

    E F_H(b_t^+(Z^H)) <= (a_(R,eps)-B_*)n, |t|>=theta,
    B_*=(3/2)q b_0 d_0 exp(-theta^2/2)theta^2.             (29)

The elementary bounds pi<22/7, exp(-x)>=1-x and sqrt(2)<99/70 give
d_0>159/1000, b_0>1/35 and exp(-theta^2/2)>999/1000. In particular
(also verified in Certificate Section 9),

    A > 47223/7000000000000000,
    B_* > 1429569/70000000000000000.                       (30)

### 11.3. Stability along paths and selection of gain windows

We prove the local covariance control needed to use (25) on H, which
can have arbitrarily short cycles. Let u_0,...,u_d be a tree path, with
R>=2d, and use the tree profiles a_i=a_(u_i) from Section 4. For the
projection pi to H set

    p_i(x)=sqrt(sum_(z:pi(z)=x) a_i(z)^2),
    D(i,j)=<p_i,p_j>-<a_i,a_j> >=0.

Then

    D(0,d) <= d 2^((d-1)/2) sum_(i=0)^(d-1)D(i,i+1)
                                                +2sqrt(2d/R). (31)

Here is a proof including truncation errors. Restrict every profile to
K=B_T(u_0,R-d), obtaining b_i. The ball B_T(u_i,R-2d) lies in K,
so at most 2d radial shells are removed and the removed squared norm
e_i is at most 3d/(1+3R/2)<=2d/R. On K all profiles are positive and
their pointwise ratios lie between 2^(-|i-j|/2) and 2^(|i-j|/2).
On each nonempty fiber write b_i=q_i e_i^*, where q_i is its norm and
e_i^* a unit vector. The truncated excess satisfies

    D_K(i,j)=(1/2)sum_x q_i(x)q_j(x)||e_i^*(x)-e_j^*(x)||^2.

The triangle and Cauchy--Schwarz inequalities give
||e_0^*-e_d^*||^2<=d sum_i||e_i^*-e_(i+1)^*||^2. The ratio bounds give
q_0q_d<=2^((d-1)/2)q_iq_(i+1). This proves (31) for D_K without
the error term. Splitting a fiber into K and its complement cannot
increase its summed affinities, by

    sqrt((a+c)(b+d)) >= sqrt(ab)+sqrt(cd).

The complement has nonnegative excess, so D_K(i,j)<=D(i,j).
Restoring a profile changes its collapsed norm vector by at most sqrt(e_i)
in Euclidean norm. Thus the collapsed endpoint affinity increases by at
most sqrt(e_0)+sqrt(e_d), while the restored uncollapsed inner product
is nonnegative. This gives (31).

For d<=5 the multiplier in (31) is at most 20 and its additive term
is less than 7/sqrt(R). Adjacent lifts have distinct base images because
H has no loops, so (1-eps^2)D(i,i+1)=delta_(pi(u_i u_(i+1))).
For a general pair, the repeated base-coordinate noise contributes
eps^2 if its images coincide and zero otherwise. The lower direction
follows directly from D>=0. Therefore, before assuming any injectivity,

    |Cov(Z_pi(u_0)^H,Z_pi(u_d)^H)-w_d|
       <=20 sum_(i=0)^(d-1)delta_(pi(u_i u_(i+1)))+E,
    E=7/sqrt(R)+max_(1<=d<=5)|r_d(R)-w_d|+2eps^2.          (32)

For every oriented edge (w,v) of H lift its fourteen-vertex window to
the universal cover, and let L(w,v) be the sum of delta_e over its
thirteen edges, counting lift multiplicities. There is the exact identity

    sum_(oriented (w,v)) L(w,v)=26 Delta.                  (33)

Indeed a fixed base edge is counted twice as root edge, eight times
as a first-level edge, and sixteen times as a second-level edge.
For the latter counts choose its endpoint nearer the root (two ways),
each successive nonbacktracking continuation (two ways per step),
and the root-edge orientation (two ways): 2*2*2=8 and
2*2*2*2=16. Paths with repeated images are still counted with their
lift multiplicities, so short cycles do not change the identity.

At each w choose a neighbor v(w) minimizing L(w,v), with deterministic
tie breaking, and set L_w=L(w,v(w)). Averaging the three choices gives

    sum_w L_w <= (26/3) Delta.

Require E<=tau. Call w good when L_w<=tau. The window is path-closed,
so (32) puts every pair covariance in a good window within
20tau+E<=21tau=zeta_g of critical. Its fourteen images are distinct:
otherwise the covariance of a coincident pair would be one, whereas
its critical covariance is w_d<=w_1 and
1-w_1>2/35>zeta_g. Thus the precise graph condition for the gain event
in Section 11.1 holds; no induced tree among buffer vertices is required.

If G is the set of good centers, the preceding sum bound gives

    |G| >= max(0,n-26 Delta/(3tau)).

Greedy selection gives an independent subset S of G with |S|>=|G|/4:
each selection removes the chosen vertex and at most three neighbors.
The orientations, branch labels, good set and independent subset are
fixed from the graph and its profile covariances, before sampling any
keys or choosing a threshold. This single S will be used throughout.

### 11.4. Gains survive the deterministic interpolation

For anchor bits b set

    N_S(b)=sum_(w in S)1[P_w(b)],    K_S(b)=F_H(b)-N_S(b).

First check the local gain at any clamped minimum x(b). On P_w the
complete crossing components meeting the roots are a-v-b and w-c.
The equal-sign child edges isolate them from every buffer vertex and
the outside. The path a-v-b has the unique minimum cover {v}. The
penalized vertices of any minimum Q_b form a minimum crossing cover,
so a,b are unpenalized and va,vb have bit zero. Clamping fixes vw and
wd to one, the two child edges of a,b,c to zero, and those of d to one.
These conclusions are unaffected by extra edges at the buffer vertices.

Let m(z) assign to each vertex the majority of its three incident edge
bits. The forced edges give majorities m_v=0, m_w=m_d=1 and
m_a=m_b=m_c=0. The bit of wc is immaterial: both w and c already
have two fixed incident bits. Hence updating w to the majority of its
neighbors saves exactly one crossing edge.

For vertex bits a let T_S(a) update every w in S to its neighbors'
majority and leave the others unchanged. Both m and T_S are monotone
maps. Since S is independent, no edge joins two updated vertices, so
if d_a(w) counts the neighbors opposite to a_w,

    |cut(T_S(a))|=|cut(a)|-sum_(w in S)max(0,2d_a(w)-3).

Also |cut(m(z))|<=C(z): every edge crossing the majority cut has
exactly one minority endpoint incidence, and each split cubic vertex
has exactly one minority incidence. At a minimum, C(x(b))<=Q_b(x(b))
=F_H(b), and each occurring pattern supplies the stated saving. Thus

    0 <= K_S(b) <= F_H(b).                                (34)

We claim that for every anchor ordering,

    cw(H) <= max_(prefix b) K_S(b)+49.                     (35)

For consecutive prefixes b<=c differing at u, Section 3 supplies a
monotone one-edge refinement from x(b) to x(c) with C(z)<=F_H(b)+3.
If u is outside a window whose P_w holds at b, that pattern holds also
at c. Every forced edge just identified has the same value at the two
endpoint minima, hence at every intermediate z between them, including
the initial clamping moves. The saving at w persists through the whole
refinement. Each selected window lies within distance three of its
center, so only centers at distance at most three from u can lose a
pattern. There are at most 1+3+6+12=22 such centers. Consequently

    |cut(T_S(m(z)))| <= F_H(b)+3-N_S(b)+22=K_S(b)+25.

Concatenate these edge refinements. They start at all zero and end at
all one. One edge flip changes at most two vertex majorities. Changing
one majority coordinate can affect T_S only there and at its at most
three neighbors, so at most eight transformed vertex bits change at
one step. The transformed sequence is monotone. Refine each such step
into individual vertex changes; each costs at most three cut edges,
so the width increases by at most 24 over the preceding endpoint.
Every vertex changes exactly once along the refined sequence. This is
an ordering, and proves (35).

Changing one anchor can change at most 22 pattern indicators and changes
F_H by at most one. Thus K_S is 23-Lipschitz under one anchor change.
All these statements also hold for empty S.

### 11.5. Continuous tradeoff and all-threshold concentration

By (25), linearity of expectation, and the size of S, for |t|<=theta,

    E N_S(b_t^+) >= (p/2)|S|
                 >= (p/8)max(0,n-26 Delta/(3tau)).

No independence among the local events is used. Combining this with
(27), the improvement below a_(R,eps)n is at least

    A Delta+(p/8)max(0,n-26 Delta/(3tau))
       >= min(A alpha,p/8)n.                              (36)

Indeed put x=Delta/(alpha n)>=0. The left side divided by n becomes
A alpha x+(p/8)(1-x)_+. On [0,1] this is a convex combination of
A alpha and p/8; for x>=1 it is at least A alpha. This includes the
case with no good centers. For |t|>=theta, use K_S<=F_H and (29).

The exact constants give

    A alpha > 47223/44590000000000000000000 > 1/10^18,
    p/8=1/(8*10^17) > 1/10^18,
    B_* > 1429569/70000000000000000 > 1/10^18.

For the first strict comparison, the integer inequality is 47223>44590.
Thus, whenever the finite profile parameters meet the local conditions
already specified, the same deterministic S satisfies

    E K_S(b_t^+(Z^H)) <= (a_(R,eps)-delta)n
                                      for every real t.  (37)

In particular a_(R,eps)-delta>=0 by (34). The expectation of K_S,
rather than its complement or its threshold monotonicity, is used here.

For completeness the concentration argument remains uniform although
the optimizing cover is global. An underlying g_x in (8) affects keys
at at most B_R=3*2^R-2 vertices; xi_v affects only its own key. The
23-Lipschitz bound and the variance argument of Section 6 therefore give

    Var K_S(b_t^+) <= 23^2 n(B_R^2+1).

For a key-count W(I) on an interval I, the same argument gives
Var W(I)<=n(B_R^2+1). Unit normal marginals give
E W((s,t])<=n(t-s)/sqrt(2pi), and either tail count beyond +/-T
has expectation n Phi(-T). Moreover

    |K_S(b_t^+)-K_S(b_s^+)| <=23 W((s,t])  (s<t).

In either tail, (34) and the one-anchor bound for F show that K_S is
at most the corresponding tail count. Fix any error u>0. Choose T
with tail expectation at most un/4 and a finite mesh of [-T,T] so
that 23 times each interval's expected count is at most un/4.
Chebyshev and a union bound over this fixed mesh imply, for all
sufficiently large n, a positive-probability event on which: every
mesh K_S differs from its expectation by at most un/4; 23 times each
interval count's deviation is at most un/4; and each tail deviation is
at most un/4. Hence (37) and interval interpolation give

    sup_t K_S(b_t^+) <= (a_(R,eps)-delta+u)n.              (38)

The tails cost at most un/2 and the coefficient a_(R,eps)-delta is
nonnegative. All variance constants and mesh sizes depend only on the
fixed radius and u, not on H or S. Intersect the event with distinctness
of all keys, which has probability one because eps>0. The decreasing-key
prefixes are then exactly the anchor vectors needed in (35).

### 11.6. The uniform limit and all cubic multigraphs

Keep J=65 fixed. For each r>0 choose one finite R>=10 and positive
eps<1 so that all of the following hold:

    max_(1<=d<=7)|(1-eps^2)r_d(R)-w_d| <= zeta_c/100,
    1-c_R>=b_0,       E<=tau,
    a_(R,eps) <= Gamma_65+r/3.

These conditions are simultaneous consequences of (11)--(14) and
E->0. The large constant C_65 remains fixed. For every input graph use
the finite cover from Section 6; the cover degree cancels from every
estimate, and no input graph is required to have large girth.
Apply (38) with u=r/3 and take n sufficiently large for its uniform
concentration statement and for 49/n<=r/3. Equation (35) gives

    cw(H) <= (Gamma_65-delta+r)n
           <= (gamma_new+r)n.                            (39)

The radius, positive noise, mesh and size threshold depend only on r.
In particular (39) preserves the entire delta as a uniform little-o
improvement; no fixed fraction of delta is spent on approximation errors.

We now reapply the graph operations of Sections 8--10 with gamma_new
in place of gamma. This parameter change uses (39), not the old
Numerical fact C with gamma replaced. Explicitly, for every r>0,
finitely many small sizes and (39) give a constant C'_r>=0 for which
every kernel satisfies

    cw(H)<=(gamma_new+r)|V(H)|+C'_r.

Lemma 8.1 is independent of the coefficient. Its bound on kernel size,
and the centroid and bridge argument of Proposition 9.3, consequently
give, for every connected loopless subcubic G on N vertices,

    cw(G)<=(2gamma_new+2r)nu(G)
                         +(C'_r+9)ceil(log2(N+1))^2.

Here gamma_new>0, so the cycle-rank inequalities have the same direction.
Deleting the loops of a connected cubic multigraph K on D vertices
preserves its cuts and leaves cycle rank at most D/2+1. Hence

    cw(K)<=(gamma_new+r)D+2gamma_new+2r
                         +(C'_r+9)ceil(log2(D+1))^2.       (40)

Define e_new(d) as the maximum of max(0,cw(K)-gamma_new d) over
connected cubic multigraphs of order d, or zero if there are none, and
put h_new(D)=max_(1<=d<=D)e_new(d). These are finite maxima, since
the degree bound limits all loop and edge multiplicities. The function
is nonnegative and nondecreasing. Given any s>0, first choose r<s/2
in (40), then d_0 sufficiently large; this gives e_new(d)<=s d for
all d>=d_0. Thus

    h_new(D)<=max(max_(d<d_0)e_new(d),sD),

and h_new(D)/D->0. Finally concatenate optimal component orderings of
an arbitrary cubic multigraph. Its width is the maximum component
width, bounded by gamma_new D+h_new(D) because gamma_new>0 and h_new
is nondecreasing. This proves Theorem 2. QED.
