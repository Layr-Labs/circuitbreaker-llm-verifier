# Cubic cutwidth below the limiting Gaussian crossing-cover coefficient

## 0. Statement and proof map

A cubic multigraph has degree three at every vertex, counting a loop twice. For a vertex ordering, count the edges with exactly one endpoint in each prefix. Its cutwidth is the minimum over orders of the maximum of those counts. Each crossing edge counts once.

**Theorem 1.** There exists a single nonnegative nondecreasing function h with h(n)/n tending to zero such that every cubic multigraph H on n vertices satisfies

    cw(H) <= gamma n + h(n),    gamma = 0.130647519423167665.

Let w_d=(1+d/3)2^(-d/2). For the centered unit-variance Gaussian wave on the cubic tree with covariance w_d at distance d, let p_j be the probability that consecutive vertices alternate signs along a fixed j-edge path. Define the finite sums

    Gamma_J = sum_(j=1)^J (-1)^(j+1) 3*2^(j-2) p_j.

Gamma is the limit of the alternating sums and the infimum of their odd truncations, as proved in certificate.md. The strengthened comparison proves that the worst-case asymptotic cubic cutwidth per vertex is at most Gamma minus 0.000000047. The displayed rational uses Gamma <= Gamma_47 < 0.130647566423167665. The finite computations check the series enclosure and the two independent parts of the layout improvement.

Sections 1–6 give galarnon42's matching-cover layout and finite-profile comparison, building on Arthur Ostapenko and McGuire, Schlesinger and Williams. Section 7 supplies the earlier retiming lemma and the elementary pivotal configuration from our first submission 4c133ffa-7452-42c6-9629-d7821ee890d4. Section 7.5 applies the quantitative proofs in certificate.md Sections 8 and 9. Section 7.6 gives the new two-pass ordering and the nine-edge count. Sections 8–10 transfer the result to arbitrary cubic multigraphs with a uniform remainder. Matching-cut submodularity is also treated by Jeong, Sæther and Telle, arXiv:1507.02384.

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

By the preceding finite-dimensional convergence, for each fixed odd J,

    beta_(R,eps) -> 0,
    Gamma_J^(R,eps) -> Gamma_J.                          (13)

These are limits of finitely many probabilities only. The coefficient
C_J=1+3J*2^(J-1) is finite and is fixed before taking either limit.

## 6. Expectations and concentration at all thresholds

Fix R,eps, and an odd J. Choose, by Lemma 4.1, a degree-m cover H' of H
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
    a_(R,eps)=Gamma_J^(R,eps)+C_J beta_(R,eps).          (14)

There is no loss depending on the cover degree or on short cycles of H.
The left side of (14) is nonnegative, so a_(R,eps)>=0. By (13), this
quantity tends to Gamma_J.

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

## 7. A uniform improvement of 10^(-14) by retiming vertex keys

For every fixed odd J for which the path bound of Sections 5–6 is used, the asymptotic cutwidth coefficient is at most

    Gamma_J - 10^(-14).                                             (S7)

The same saving therefore holds below the infimum of the odd bounds. The size thresholds and Gaussian profiles may depend on J. This section proves (S7) for connected simple cubic graphs, after which the kernel and bridge reductions apply as before. Write F(B) for the crossing-cover number of a vertex set B. By bipartite matching-cover equality, it is also the maximum size of a matching in the crossing graph. Changing one vertex changes F by at most one.

### 7.1 A local change in a vertex ordering

At a center o name the neighbors a,b,u. Let a1,a2 be the two neighbors of a beyond o, let w be a neighbor of u other than o, and let w1,w2 be the two neighbors of w beyond u. Require these named vertices to be distinct. Put

    T(o)={a,b,a1,a2,w,w1,w2}.

If o is outside B and T(o) is contained in B, then F(B+o)=F(B)-1. To prove this, delete o from the crossing graph, obtaining K. The vertex a is isolated there, so the old edge oa augments a maximum matching of K. After changing the bit at o, the only possible crossing edge at o is ou. If this edge crosses, u is outside B and w has just one crossing neighbor, u. Any matching using ou can replace that edge by uw, obtaining a matching of K of the same size. This proves the decrease. Complementing all the bits gives the corresponding decrease when o is inside B and T(o) is outside. The proof uses the specified full neighbor lists at o,a,w; extra edges between boundary vertices at distance three are allowed.

Given any set of eligible centers, select a subset S in which no selected center lies in another's trigger set. If there are g eligible centers, the conflict graph has at most 7g edges. Give its vertices independent continuous priorities and keep every vertex of lowest priority in its closed neighborhood. The expected size is sum 1/(1+degree), at least g/15 by Cauchy-Schwarz. Fix any such independent set of size at least g/15.

Let Z be the original real vertex keys and order by increasing key. At v in S replace the key by

    min(Z_v,max_{x in T(v)} Z_x)  when Z_v>=0,
    max(Z_v,min_{x in T(v)} Z_x)  when Z_v<0.                     (21)

Other keys keep their original values. At a threshold t, take B={v:Z_v<t}. Each vertex whose bit changes under (21) saves exactly one unit of F, by the preceding matching argument or its complement. The changes can be performed in any order since no changed vertex is in another's trigger set. Thus the new value equals the old value minus the number of changed vertices. Each new key is one of the original keys in its radius-three ball. The original keys are almost surely distinct owing to the small independent noise. At most 22 new keys can equal one original key, by the cubic ball-size bound. Resolving ties arbitrarily therefore adds at most 22 to the largest threshold F. The layout conversion in Lemma 3.2 adds at most six.

### 7.2 A finite-radius Gaussian construction and its comparison

Let f_x be the normalized profile on the universal cubic tree proportional to alpha^d(x,z), where 0.7<=alpha<1/sqrt(2). Its squared norm before normalization is (1+alpha^2)/(1-2alpha^2). Counting the off-path branches by their attachment vertex as in Lemma 5.2 gives its cosine at distance d,

    r_d(alpha)=alpha^d[(d+1)-(d-1)alpha^2]/(1+alpha^2).

Thus r_d(alpha) tends to w_d=(1+d/3)2^(-d/2) as alpha increases to 1/sqrt(2). Write b_x(p) for the norm of f_x on the fiber over the base vertex p. The pushed row b_x depends only on the projection of x. On an adjacent tree pair define its excess

    D_xy=<b_x,b_y>-<f_x,f_y>.

It is nonnegative by Cauchy-Schwarz on each fiber and is the same for all lifts of a given base edge. Define eta=10^(-4). Call a center good if all edges in its lifted radius-three tree ball have excess at most eta.

Use the infinite alpha-profile only to define these excesses and estimate covariances. On the actual graph truncate the profile at radius R, normalize it, push it to the base, and add independent tie noise as in (8). For fixed alpha the omitted squared mass is a geometric tail tending to zero uniformly in the root, and pushforward contracts row distance by the reverse triangle inequality on each fiber. Consequently these finite rows and their covariances approach the full profiles uniformly over all graphs and all individual pairs. On a finite cover of sufficiently large girth the local truncated profiles have exactly the corresponding tree covariances. The comparison, exact lifting and concentration arguments of Sections 4–6 apply to these finite rows unchanged. For each fixed J the alternating-path expectation of the high-girth cover at zero has upper bound Gamma_J+e, with e arbitrarily small: first take alpha close to its critical value, then choose R large and the positive tie noise small. This follows from convergence of the finitely many path and three-neighbor probabilities used in (10), including vanishing of the degree-three error. These choices can meet simultaneously all the fixed-patch accuracy conditions below, before the graph size is chosen.

For a good center o we need a more precise local Gaussian coupling. Let u_x(p) be the unit direction of f_x on fiber p, and let U_x map the pth coordinate vector to u_x(p) supported on that fiber. This is an isometry and U_x b_x=f_x. Then

    2 D_xy=sum_p b_x(p)b_y(p)||u_x(p)-u_y(p)||^2.

Adjacent positive alpha-profile entries have ratios in [alpha,1/alpha], so their fiber norms do as well. For a tree path o=x_0,...,x_k=x, k<=3, one has

    b_x(p)^2 <= alpha^(-(2(k-i)+1)) b_{x_(i-1)}(p)b_{x_i}(p).

Telescope the differences of the isometries along the path, applying each to b_x. Their norm bounds from the last two displays give

    ||f_x-U_o b_x|| <= sqrt(2 eta) sum_{j=0}^{k-1} alpha^(-(2j+1)/2).

The right side is zero at k=0 and respectively less than .01691,.04107,.07556 at k=1,2,3 (use alpha>=0.7). All lifted vertices of this ball have distinct base projections: otherwise their ideal rows would be within 2(.07556) in norm, whereas distinct rows at distances 1 to 6 have limiting norm distance at least sqrt(2-2w_1)>.33. The trigger in Section 7.1 is therefore available at each good center.

Using the same Gaussian white noise on the rows f_x and U_o b_x couples the ideal alpha-values and full base values, with difference standard deviations bounded by the last display. At a fixed finite patch the ideal alpha Gaussian vector also couples to the critical one with arbitrarily small difference standard deviations by continuity of the positive matrix square root. Truncation and tie noise give another arbitrarily small rowwise error. These Gaussian couplings can be joined over their common Gaussian vector by sampling conditionally independent residuals, and standard deviations of successive differences add as upper bounds. Thus for the actual base keys versus the critical wave we may take standard-deviation bounds

    .0001, .017, .0412, .0757 at respective distances 0,1,2,3.      (22)

The coupling is needed separately for each center to bound expectations.

### 7.3 A quantified change at good centers

For every 0<=t<=1.25, the critical wave has an event of probability at least 4.06725 times 10^(-13) on which its root lies above t by more than .004 and all seven trigger keys lie below t by more than .15 at distance one, .33 at distance two and .64 at distance three. To verify this set lambda=2 sqrt(2), rho=lambda/3. Target the following values:

    o=t+.02, a=t-.33, b=t-.20, u=lambda o-a-b, w=t-.45,
    a1=a2=(lambda a-o)/2,   w1=w2=(lambda w-u)/2.

The Gaussian root O is standard normal. With independent standard G1,G2 its three neighbors are rho O plus (G1, -G1/2+sqrt(3)G2/2, -G1/2-sqrt(3)G2/2)/3. A subsequent sibling pair with parent C and grandparent P is (lambda C-P)/2 plus and minus G/(2sqrt(3)), with a fresh independent standard G. These formulas follow from the wave equation: sibling differences have variance 2(1-w_2)=1/3 and zero covariance with earlier variables outside their two descendant branches, by equal tree distances. In exposing the tree a level at a time all these new contrasts are independent.

Center the six independent normals O,G1,G2,G_a,G_u,G_w at

    o, 3(a-rho o), (2/sqrt(3))[3(b-rho o)+(3/2)(a-rho o)], 0,
    2sqrt(3)[w-(lambda u-o)/2], 0.

Use interval half-widths .015,.06,.06,.30,.12,.40, respectively. The absolute centers of the first, second, third and fifth intervals are at most 1.27,1.047,1.363,4.398; the last three bounds follow from their affine expressions on [0,1.25]. The induced errors in a,b,u are below .0342,.0415,.0415, in an a-child below .1424, in w below .1008 and in a w-child below .2788. The target a-child is (.91421356...)t-.47669047..., and the target w-child equals t-.92968037.... These bounds give all claimed margins. Each normal interval has mass at least its length times its least density. The six masses exceed .005,.025,.017,.22,.000003,.29, respectively, whose product is 4.06725 times 10^(-13).

From (22) and the Gaussian tail bound, the probability of failure under coupling is at most

    2 exp(-.004^2/(2·.0001^2)) + 4 exp(-.15^2/(2·.017^2))
        + 6 exp(-.33^2/(2·.0412^2)) + 4 exp(-.64^2/(2·.0757^2))
      < 7.2 times 10^(-14).

Hence every selected good center changes under (21) with probability greater than 3.34 times 10^(-13) on this threshold range. Indeed the actual root is above t>=0 and all trigger keys are below t. Reflection gives the same guarantee on [-1.25,0] for the branch that increases the center's key.

If at least half the n base vertices are good, the independent selection has size at least n/30. At each threshold |t|<=1.25 the expected saving is therefore more than 1.11 times 10^(-14)n. For mesh thresholds the number of changed centers has variance O_R(n), since its indicators depend only on Gaussian inputs within a fixed radius. Between mesh points the changed-center count varies by at most 22 times the number of original keys in that window, since only centers within distance three can depend on a crossing key. The Section 6 finite-grid and concentration argument therefore proves one realization with throughout the central interval a loss arbitrarily close to this expected saving from the original upper bound. In the tails the new F is at most the original F, hence at most the corresponding vertex tail count. The standard normal tail at 1.25 is below .11 (integrate the lower polynomial 1-z+z^2/2-z^3/6 for exp(-z) with z=x^2/2 from 0 to 1.25, including the Gaussian normalizing factor). Every odd Gamma_J exceeds .113: a matching of size at least N_1-N_2 exists by removing one edge of each overlapping edge pair as needed. On the high-girth covers this lower bound on F and the incumbent upper bound give, after taking profile limits,

    Gamma_J >= (3/2)p_1-3p_2
             = (3/(2pi))[arccos(5/6)-arccos(2sqrt(2)/3)] > .113.

Here 2p_1=2p_2+Pr(the two outer path signs differ). Alternating cosine Taylor bounds give the first arccos above .580 and the second below .341. Choose all profile, mesh and concentration allowances to total less than 5 times 10^(-16)n, and finally absorb the constant tie and layout overhead. This proves (S7) for the good-center case.

### 7.4 A strictly decreasing covariance comparison at an edge

If fewer than half the vertices are good, at least n/88 distinct base edges have D_e>eta; each observed edge can be charged to at most 44 centers, by the radius-three ball bound at its two ends. Choose the finite profiles so on each such edge their base-field covariance exceeds its high-girth-cover value by at least .99 eta. This is possible by the uniform row approximations in Section 7.2.

Let X denote the high-girth-cover field and Y an independent copy of the base field repeated in the covering fibers. Write A,C for their covariances. Hellinger comparison gives C>=A entrywise, and both diagonals equal one. For Z=sqrt(1-s)X+sqrt(s)Y, differentiating the zero-threshold expectation of F on the cover gives a sum over unordered pairs i,j with terms

    (C_ij-A_ij) density_(Z_i,Z_j)(0,0) E[Delta_ij F | Z_i=Z_j=0].

This is Gaussian integration by parts as in Lemma 1.1, integrating the two derivative coordinates across zero; the independent tie noise makes the covariance positive definite for s<1. All terms are nonpositive. We exhibit a quantitative strict decrease on adjacent pairs.

Condition the ideal tree wave at the endpoints i,j of an edge to be zero. At i let its positive outer neighbor be a and its negative outer neighbor b; at j let its positive outer neighbor be c and its negative outer neighbor d. Require each of a,d to have one opposite-sign outward child and one like-sign outward child, with both outward children of the opposite child of its own sign. Require absolute margin at least m=.04 on each of a,b,c,d and on the specified children and grandchildren of a,d.

On these signs the mixed difference is at most -1, irrespective of the remaining crossing graph. Indeed the opposite child of each of a,d is a leaf in that graph before the two central bits are chosen, and its parent has no other outward crossing. The two private leaf edges can always be used in some maximum matching, by replacing a match to the center by the leaf edge, or adding the leaf edge when both ends are otherwise free. Delete the four private vertices and the two central vertices. Let M be the matching size of the remaining crossing graph, A_b and A_c the indicators that its vertex b, respectively c, can be exposed by a maximum matching. If both are one there is a maximum matching exposing both simultaneously. The symmetric difference of two witnessing maximum matchings has only even paths and alternating cycles, since an odd path would augment one matching. Endpoints of a path are in the same color class, whereas b,c are in opposite classes. Choose independently the witnessing matchings on the relevant components; vertices unmatched in both present no constraint.

After subtracting M+2, the matching sizes for the two equal central colors are A_b and A_c. One opposite-color assignment allows only the central edge besides the private matches, contributing one. The other contributes at least max(1,A_b+A_c) by the central edge or the exposed attachments. This proves the claimed mixed difference.

For clarity the pivotal probability has a short, exact finite lower bound. At either centered endpoint the magnitude x of its two outer scores has the half-normal law |N(0,1/12)|, independently on the two sides of the edge. Let f_A be that density and f_D,H_D the density and distribution function of |N(0,1/3)|. Using the sibling contrasts from Section 7.3, the amplitude k of a private parent's opposite child is (|D|-lambda x)/2 when positive. The other child has amplitude lambda x+k. Both children of the opposite child have its sign with margin at least m precisely when their absolute contrast is at most lambda k+x-2m. Thus the probability of the event is

    Q = [ integral_{x>=m} f_A(x)
          integral_{k>=m} 2 f_D(2k+lambda x) H_D(lambda k+x-2m) dk dx ]^2.

No prescribed neighbor index is required: exactly one neighbor at each end has the specified sign. The following integer-rational sum gives Q>0.018079200775. Put h=.01, A(z)=(11/4)(1-6z^2/64)^64, D(z)=(11/8)(1-(3/2)z^2/64)^64, and

    H(z)=(11/8) sum_{r=0}^7 (-3/2)^r min(z,5/4)^(2r+1)/(r! (2r+1)).

For ai=4,...,59 and ki=4,...,79 sum, with x=ai h and k=ki h, the terms

    2h^2 A(x+h) D(2(k+h)+2.829(x+h)) H(2.828k+x-2m).

The floor of 10^12 times the square of the sum is 18079200775. The exact-arithmetic checker for the pivotal probability reproduces this finite sum. The expressions A,D bound the respective densities below on the ranges used because exp(-z)>=(1-z/64)^64 for 0<=z<=64, and pi<22/7 bounds both normalizing factors. Integrating the odd Taylor lower polynomial of exp(-3z^2/2) yields H as a lower bound for the distribution, using monotonicity for the cap. Taking the upper box endpoints for decreasing densities and the lower ones for the increasing distribution, with 2.828<lambda<2.829, therefore gives the stated lower bound.

For sufficiently accurate finite parameters, the cover law conditional on X_i=X_j=0 has the same pivotal event with margins .039 and probability at least .0179. Conditional Gaussian convergence holds because the two-by-two endpoint block is invertible, and the ideal .04-margin event lies strictly inside the finite target inequalities. For 0<=s<=5.5 times 10^(-6) this event holds for Z given Z_i=Z_j=0 with probability at least .0165 after removing the margin condition. To see this, write X_v=L_v X_P+E_v for each of at most 28 outside-patch values, with P=(i,j). In the ideal wave the regression rows on one side at depths 1,2,3 are (sqrt(2),-1/2), (3/2,-1/sqrt(2)), (sqrt(2),-3/4), by the conditional child means; on the other side interchange the coordinates. For accurate finite profiles each corresponding row has l1 norm less than 2.5. The residual vector E is independent of X_P,Y. Conditionally on Z_P=0,

    Z_v=sqrt(1-s) E_v + sqrt(s)(Y_v-L_v Y_P).

The centered perturbation before multiplying by sqrt(s) has conditional variance at most its unconditional variance, at most (1+2.5)^2. The probability any perturbation exceeds .038 in absolute value is at most

    56 exp[-.038^2/(2·3.5^2·5.5·10^(-6))] < .00125.

As .039 sqrt(1-s)>.038, the stated pivotal probability follows. The adjacent zero-pair density is at least 1/(2pi), so integration of these terms from 0 to 5.5 times 10^(-6) yields a per-vertex saving of at least

    (1/88)(.99·10^(-4))(5.5·10^(-6))(.0165/(2pi)) > 1.6 times 10^(-14).

All other pairs and the rest of the interpolation have nonpositive contributions. Exact lifting divides the cover expectation by its degree, and Lemma 1.1(b) gives the same upper bound at every threshold on the base. Select the finite-profile errors and the original concentration and layout errors below the slack between this saving and 10^(-14). This proves (S7) also for the bad-center case.

### 7.5 Combining the two quantitative improvements

Use the excess D_e of Section 7.2 with the new fixed cutoff eta=0.000425. A good center has every edge of its lifted radius-three ball below or equal to this cutoff. Let g be their number on an n-vertex connected simple cubic graph. Take finite profiles as close to critical as required in both companion proofs. The number of retained good centers is at least g/15 by Section 7.1. For every |t|<=.67 the good-center calculation in certificate.md Section 8 gives a uniform trigger probability q with q/15>0.00000002. Section 7.1 shows that each trigger saves one unit of the crossing-cover number, with no double counting.

For a center that is not good, some edge in its ball has D_e>eta. One endpoint of that edge is within distance two of the center. A given edge can therefore be observed by at most 2+4+8=14 centers. There are at least (n-g)/14 such distinct heavy edges. The finite-profile covariance surplus on each is at least .99 eta. The integrated two-orientation witness estimate of certificate.md Section 9 gives per unit of surplus at least K=0.000672905379119700, up to arbitrarily small finite-profile losses. All other derivative terms have nonpositive sign. After exact division by the cover degree, this saves at least .99 K eta (n-g)/14 from the unretimed expectation. Lemma 1.1(b) transfers this zero-threshold upper bound to every threshold. The good-center retiming then saves the additional qg/15 at each central threshold. In particular

    .99 K eta /14 > 0.00000002,      q/15 > 0.00000002.

For completeness the outer thresholds already have slack. For adjacent centered normal scores of correlation r, the probability of opposite threshold bits is

    p(t,r) = (1/pi) integral from r to 1 of exp(-t^2/(1+s))/sqrt(1-s^2) ds
           <= exp(-t^2/2) arccos(r)/pi.

Indeed differentiate the bivariate normal lower-quadrant probability at (t,t) with respect to r. The density satisfies partial_r f=partial_x partial_y f, and the integrated derivative equals f(t,t). At r=1 the two threshold bits coincide. This yields the formula by summing the two crossing orientations. Count every crossing edge as an upper bound on the matching size. Hellinger comparison gives adjacent correlation at least the truncated tree value, arbitrarily close to 2sqrt(2)/3. Since arccos(2sqrt(2)/3)<.34, pi>3.14, and exp(-.67^2/2)<1-.67^2/2+(.67^2/2)^2/2, the limiting upper bound (3/2)p(.67,2sqrt(2)/3) is less than .1302. The cosine estimate follows from the Taylor upper bound cos(.34)<=1-.34^2/2+.34^4/24<2sqrt(2)/3. Every odd Gamma_J is at least Gamma >= Gamma_46 > .1306 by certificate.md. Thus the outer thresholds satisfy the same claimed saving. Retiming never increases their crossing-cover number.

All inequalities used for the reported saving have positive slack. Choose alpha near critical, then a finite truncation radius and positive tie noise, to meet their finitely many accuracy conditions simultaneously. At any fixed threshold the retimed value is the original F minus a sum of finite-range trigger indicators. Their variances are O_R(n), as in Section 7.3. On central mesh intervals use the at-most-22 local multiplicity bound there. On the bounded outer intervals use the original mesh bounds of Section 6. Beyond a sufficiently large fixed tail cutoff use the vertex tail counts. With concentration and a sufficiently fine fixed mesh these losses are arbitrarily small per vertex. Resolving the ties and converting the matching cover to a graph ordering add at most 22+6. For each fixed odd J this proves the asymptotic graph bound Gamma_J minus 0.00000002.

### 7.6 Two passes on nine-edge patches

Use the named neighbors and Hellinger excess of Sections 7.1–7.2 of this manuscript. Let η = 17/40000, and call an edge heavy when its excess exceeds η. A named tree patch comprises the four complete stars at o, a, u, w. Their union has nine tree edges. Every covariance comparison in Section 8 of the certificate uses a path inside this union, and every defect estimate there uses one of the four complete stars. Its bounds also show distinct projections for the named vertices when those nine edges are light: for distinct named tree vertices at distance d ≤ 6, the limiting correlation is bounded above by w_d + d²(√2)^(d−1)η < 1. Thus one light named patch suffices at a center. Its favorable one-orientation trigger probability exceeds q = 361900864/10^15, up to arbitrarily small uniform finite-profile error.

Let h be the number of heavy edges and g the number of vertices admitting at least one light named patch. Choose a root uniformly, a uniform ordering of its three incident edges, and a uniform choice for the child w of u. The projection of any of the nine tree edges is uniform on base edges. Indeed each root-to-edge route is a stationary nonbacktracking walk started with a uniform directed edge, whose continuation chooses either permitted edge uniformly. This marginal argument also applies to unlabeled sibling pairs ordered by a fair coin. There are 3n/2 base edges. The union bound gives n−g ≤ 6h even if some projected edges coincide. At each eligible center fix one light patch.

Let T(v) be its seven trigger vertices. Mark each eligible center with independent probability 1/5, independently of its Gaussian scores. In the first pass replace a marked center's score by the minimum of its score and the maximum score in T(v), reading all scores before the pass. In the second pass replace every eligible center's score by the maximum of its first-pass score and the minimum first-pass score in T(v), reading all scores after the first pass. At any threshold every change in the first pass adds an outside center whose entire trigger was already inside. All additions can be made successively and each saves one matching edge by Section 7.1. Removing inside centers in the second pass has the complementary property. Both passes can only decrease the crossing matching size at every threshold.

For 0 ≤ t ≤ .67 the favorable original event therefore saves one unit when the center is marked, with probability at least q/5. For −.67 ≤ t ≤ 0 reflect the original Gaussian field. A center in the favorable removal event survives to the second pass whenever none of the at most seven eligible vertices in T(v) was marked. Its conditional probability is at least (4/5)^7 > 1/5. Each successful removal counts once. This yields qg/5 in expected savings on the entire central interval, up to arbitrarily small finite-profile losses. In either tail use the already stronger outer-threshold bound in Section 7.5.

Independently, the preceding covariance comparison saves at least .99 Kηh at all thresholds before the two passes, where K = 672905379119700/10^18. Combining these on the same graph gives a saving per vertex at least min(q/5,.99 Kη/6), which exceeds 0.000000047. Both strict inequalities are rational. Each final score is an original score within radius six, obtained using only local scores and marks. At most 190 vertices can inherit any given score. At a fixed threshold the crossing matching size is 1-Lipschitz in each final bit. A change in any one independent input therefore affects a bounded number of bits for each fixed profile radius. The bounded-difference variance estimate and mesh argument of Section 6 apply, with interval counts bounded by 190 times the original counts. This provides one realization valid at every threshold. Ties cost at most 190 more matching edges, and Lemma 3.2 adds six. Taking the profile and graph-size limits for each odd path truncation gives the stated bound relative to Γ. Sections 8–10 of the graph proof extend it to all cubic multigraphs with a uniform sublinear remainder.


**Theorem 7.1.** For every positive epsilon there is N_epsilon such that every connected simple cubic graph H on n>=N_epsilon vertices satisfies cw(H)<=(gamma+epsilon)n. Apply the preceding proof with the fixed J=47 and Gamma_47 < 0.130647566423167665. All cutoffs are uniform in the input graph. The following reductions complete the proof of Theorem 1.

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
