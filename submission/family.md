# Norm-one language and transfer to full \(B_2\)

The graph constant from [proof.md](proof.md) is
\[
a_0=\frac{3696211249941850775474886159547}{28796518400000000000000000000000},\qquad c=1+1/(2a_0)=\frac{18094470449941850775474886159547}{3696211249941850775474886159547}.
\]
For connected simple cubic graphs of order \(v\), the graph proof and
the quantitative estimate below give
\[
\operatorname{cw}(H)\le a_0 v+O\bigl(v/(\log_2(v+8))^2\bigr).
\]
The constant in the error is independent of the graph. Cutwidth counts
crossing edges at the cuts of a vertex order, including multiplicity
in a multigraph.

The family is the same norm-one language as in the Circuit Breaker team's *A marked Gaussian ordering with paired signed and nine-key Q bounds*, companion `family.md`, Sections 1–3. It has many norm-one inputs but excludes large accepting rectangles on every partition of the input. A circuit layout supplies a read-once program for the same accepting set, so the rectangle restriction gives a lower bound on the layout width.

## The language at every length

For \(N<16\) reject every word. At all other lengths define
\[
t=2^{\lfloor\log_4 N\rfloor},\qquad
m=\max\{j\text{ positive odd}:(t+1)j\le N-t\},\qquad
\ell=N-t-(t+1)m.
\]
The set defining \(m\) is nonempty, since \(N\ge t^2\) and \(t\ge4\). The first \(t\) input bits are the coefficients in increasing degree order below the leading term of a monic \(R(W)\) of degree \(t\) over \(\mathbb F_2\). The next \(m\) bits encode a monic \(Q(Z)\) of degree \(m\) in the same way. Require both to be irreducible over \(\mathbb F_2\). Put
\[
F=\mathbb F_2[Z]/(Q),\qquad q=2^m,\qquad E=F[W]/(R).
\]
An irreducible factor of \(R\) over \(F\) of degree \(d\) would have a root in a field of size \(2^{md}\). That root generates a copy of the \(2^t\)-element field over \(\mathbb F_2\), since its minimal polynomial over \(\mathbb F_2\) is \(R\). The tower law gives \(t\mid md\), so \(t\mid d\) by coprimality of \(m\) and \(t\). Therefore \(R\) is irreducible over \(F\) and \(E\) is a field with \(q^t\) elements.

The next \(tm\) bits encode \(z=\sum_{v=0}^{t-1}\sum_{u=0}^{m-1}z_{u,v}Z^uW^v\), with \(u\) varying fastest. This is an additive bijection from the payload space to \(E\). The last \(\ell\) bits are ignored. With a valid header, accept exactly when
\[
z^{\,e}=1\quad\text{in }E,\qquad e=1+q+\cdots+q^{t-1}.
\]
Write \(S_N\) for the accepting set at length \(N\), \(L=\bigcup_N S_N\), and \(L_N\) for the characteristic function of \(S_N\).

For a monic polynomial \(P\) of degree \(d\), irreducibility over \(\mathbb F_2\) is equivalent to
\[
P\mid X^{2^d}-X,\qquad
\gcd(P,X^{2^j}-X)=1\quad\text{for every proper positive divisor }j\mid d .
\]
Indeed \(X^{2^j}-X\) has simple roots and contains exactly the monic irreducible factors of degrees dividing \(j\). A factor of smaller degree in the first divisibility test would consequently violate a gcd test. This finite-field factor identity also supplies a header at every required degree. If \(I_d\) is the number of monic irreducibles of degree \(d\), then \(\sum_{j\mid d}jI_j=2^d\). The contribution from \(j<d\) is at most \(\sum_{j=1}^{\lfloor d/2\rfloor}2^j<2^d\), so \(I_d>0\).

The tests use modular squaring and polynomial gcd. The exponent \(e<q^t\) has at most \(tm\) binary digits, and quotient-field arithmetic and binary exponentiation have polynomial cost in \(N\). These algorithms give a single deterministic polynomial-time evaluator on every word, including invalid headers and small lengths. Take witness bound \(W(N)=0\) and let the verifier accept exactly when the witness is empty and the evaluator accepts. Membership gives the empty accepting witness, and any accepted witness implies membership. Thus \(L\in\mathrm P\subseteq\mathrm{NP}\), with all field data supplied by the actual input.

## Norm equations restrict rectangles

Fix a valid header and zero tail. The function \(E^*\to E^*\), \(z\mapsto z^e\), is a group homomorphism. Its image solves \(x^{q-1}=1\), since \(e(q-1)=q^t-1\), and hence has at most \(q-1\) elements. Its kernel has at least \((q^t-1)/(q-1)=e\) elements. The degree-\(e\) polynomial \(X^e-1\) bounds the kernel size from above by \(e\), so the fixed header and tail supply exactly \(e\ge2^{(t-1)m}\) accepting words.

Over any field \(K\), the system
\[
\prod_{j=1}^t(x_j-a_{ij})=b_i\qquad(1\le i\le t)
\]
has at most \(t!\) solutions in \(K^t\) if the \(a_{ij}\) are distinct in each fixed column \(j\). This is Theorem 3.3 of Kollár, Rónyai and Szabó, *Norm-graphs and bipartite Turán numbers*, Combinatorica 16(3), 399–406 (1996). The following argument gives the needed unrestricted-field statement. [Published paper](https://doi.org/10.1007/BF01261323)

Over an algebraic closure \(k\), set \(X=(\mathbb P^1_k)^t\) and homogenize the equations as sections of \(\mathcal L=\mathcal O_X(1,\ldots,1)\).
\[
\prod_j(u_j-a_{ij}v_j)-b_i\prod_jv_j=0.
\]
If some \(v_j=0\), the right-hand product vanishes, and every row must have a vanishing factor from a finite coordinate. A finite coordinate can supply that factor for at most one row. With fewer than \(t\) finite coordinates there is no common zero. The common zero scheme \(Z\) lies in the affine chart with all \(v_j\ne0\). It is closed there and proper over \(k\) as a closed subscheme of \(X\), so it is finite. [Affine and proper implies finite](https://stacks.math.columbia.edu/tag/01WN)

At each point of \(Z\) the local ring of \(X\) is regular of dimension \(t\). The \(t\) local equations have zero-dimensional quotient and form a system of parameters, hence a regular sequence. Outside \(Z\) one equation is a unit locally. Thus the Koszul resolution of \(\mathcal O_Z\) has terms \(\mathcal L^{-i}\) with multiplicity \(\binom ti\). After twisting by \(\mathcal L^d\) with \(d\ge t\), the Euler characteristic of each \(\mathcal L^{d-i}\) is \((d-i+1)^t\), by projective-line cohomology in nonnegative degree and the product formula. The Euler characteristic of \(\mathcal L^d\) on \(Z\) is the length of \(Z\), because the line bundle is free of rank one on each local Artinian factor. Additivity gives
\[
\operatorname{length}(Z)
 =\sum_{i=0}^t(-1)^i\binom ti(d-i+1)^t=t!,
\]
the \(t\)-th backward difference of a monic degree-\(t\) polynomial. Each geometric point contributes at least one to the length, which bounds the number of \(K\)-points. [Regular sequences](https://stacks.math.columbia.edu/tag/02JN) · [Koszul resolution](https://stacks.math.columbia.edu/tag/062F) · [Projective cohomology](https://stacks.math.columbia.edu/tag/01XS) · [Product formula](https://stacks.math.columbia.edu/tag/0BEC)

For distinct \(d_1,\ldots,d_t\in E\), the equations \((y+d_i)^e=1\) have the displayed form with \(x_j=y^{q^{j-1}}\), \(a_{ij}=-d_i^{q^{j-1}}\), and \(b_i=1\). Each column is distinct because Frobenius powers permute \(E\), and \(x_1=y\) makes the substitution injective. There are at most \(t!\) possible \(y\). In a partition of the payload bits, a side assignment injectively specifies its summand of \(z\), by setting the other side's coordinates to zero. A one-rectangle with \(t\) assignments on the first side and \(t!+1\) on the second would contradict the bound on the translates.

Set \(h=t+m+\ell\) and \(K_N=2^h(t!+1)\) for \(N\ge16\). Suppose a partition of all \(N\) bits had a rectangle \(D_1\times D_2\subseteq S_N\) with both sides of size at least \(K_N\). Fix separately one header-and-tail pattern on each side, retaining at least \(t!+1\) assignments per side. The two patterns jointly form a valid header because the rectangle is accepting. Payload restriction is injective on each retained side, giving the forbidden payload rectangle. Thus \(S_N\) excludes rectangles with both side sizes at least \(K_N\) on every partition.

The parameters satisfy
\[
\frac{\sqrt N}{2}<t\le\sqrt N,\qquad
\sqrt N-4<m<2\sqrt N,\qquad
0\le\ell\le2t+1.
\]
The bounds on \(t\) follow from \(t^2\le N<4t^2\). For \(m\), use \((N-t)/(t+1)-2<m\le(N-t)/(t+1)\) with the bounds on \(t\). The tail estimate follows from maximality of the odd \(m\). Consequently, for base-two logarithms,
\[
N-\log_2|S_N|\le t+2m+\ell=O(\sqrt N),\qquad
\log_2 K_N=O(\sqrt N\log N).
\]
These estimates apply to every integer length \(N\ge16\).

## A read-once program must have many transitions

Use a layered nondeterministic program with a single start vertex and a single accept vertex. Each layer transition either reads one input coordinate, labelling transitions by its bit, or guesses using unlabelled edges. A read coordinate is used at just one layer. An input accepts if some start-to-accept path agrees with it. Count edges with their multiplicities.

Suppose an accepting set \(S\) excludes rectangles of side sizes at least \(K\ge2\) on every partition. Let \(r\) be the number of coordinates omitted by a program computing \(S\), with \(|S|\ge2^rK\). The set \(S'\) on the read coordinates has size \(|S|/2^r\). Fixing the omitted coordinates would extend any forbidden rectangle in \(S'\) to one for \(S\). At a program vertex \(v\), let \(A_v\) be the read prefixes reaching it and \(B_v\) the read suffixes reaching accept from it. Every prefix in \(A_v\) concatenates with every suffix in \(B_v\) to give an accepting assignment. In particular, \(|B_v|\ge K\) implies \(|A_v|\le K-1\).

Choose one accepting path per \(x\in S'\). The suffix set at the start has size \(|S'|\ge K\), and the suffix set at the accept vertex has size \(1<K\). Charge \(x\) to the first traversed edge \(u\to v\) for which \(|B_v|<K\). The prefix in \(A_u\), the suffix in \(B_v\), and the bit label if the edge reads a coordinate determine \(x\). Each edge receives at most \((K-1)^2\) charges, giving the program lower bound
\[
\frac{|S|}{2^r(K-1)^2}.
\]
The first-drop argument and circuit-to-cutwidth method follow McGuire, Schlesinger and Williams. [Pinned proof](https://github.com/rrwilliams0/circuit-lb-4.6155n/blob/a56500014ba159397903dea81d41ac836f6880d5/proof.md)

## A full \(B_2\) circuit gives a subcubic wire graph

A circuit is an acyclic shared DAG on the actual \(N\) input bits. Each gate has two ordered pins and one of the sixteen binary truth tables. All gates are counted, including constant-valued tables, unary tables, projections and negated projections. Fan-out and depth are unrestricted, and both pins may have the same supplier. If a free constant is used as a pin, absorb it into the truth table and feed the other supplier into both pins. For two constant pins, use the corresponding constant table with any actual raw input feeding both pins. These changes preserve the function and add zero gates.

Write \(b_N=\lceil\log_2 K_N\rceil\). A nonempty set satisfying the rectangle bound can ignore at most \(2b_N-1\) coordinates. Otherwise start from an accepting word, take two disjoint sets of \(b_N\) ignored coordinates, put them on opposite sides of a partition, and vary them independently to make an accepting rectangle of side sizes \(2^{b_N}\). Since \(b_N=o(N)\), the selected output cannot be a free constant or an input literal at large lengths. Prune nonancestors of the selected output gate and unused inputs. Let \(s'\le s\) be the remaining gate count and \(n_u\) the syntactically used input count. The omitted count \(r=N-n_u\) satisfies \(r\le2b_N-1=o(N)\).

For each supplier with \(j\ge1\) outgoing pin occurrences, distribute its value to those pins by a binary tree with \(j-1\) degree-three copy vertices. A copy vertex requires all three incident edge bits to agree. The circuit inputs have degree one in the resulting undirected graph \(G\), the selected gate has degree two, and every other vertex has degree three. \(G\) is connected because all retained circuit objects lead to the selected output, and loopless by acyclicity. Parallel edges retain their original pin identities. The \(n_u+s'-1\) suppliers send \(2s'\) original pin occurrences, so there are \(s'-n_u+1\) copy vertices. Writing \(v=|V(G)|\) and \(\nu(G)=|E(G)|-|V(G)|+1\) gives
\[
v=2s'+1,\qquad |E(G)|=3s'-n_u+1,\qquad \nu(G)=s'-n_u+1.
\]

For a vertex order of width \(w\), put one program layer at each cut, with states the assignments to the crossing edge bits. Adjacent states agree on persistent edges and satisfy the local rule at the vertex traversed. A copy vertex checks equality. An internal gate checks its table on the ordered pins and outgoing bit. At the selected gate enforce output one. An input layer reads its bit and enforces its incident edge value. There are at most four valid assignments to a vertex's incident bits, so consecutive layers have at most \(4\cdot2^w\) transitions. The start and accept states are the empty cuts.

A program path specifies one consistent value per edge throughout its interval, so topological evaluation recovers the actual old gate values. Conversely an accepting circuit evaluation gives a valid program path. The program uses at most \(4v2^w\) edges and reads each used input once. For large \(N\), the estimates for \(S_N\) and \(K_N\) give \(|S_N|\ge2^r K_N\) uniformly over the circuits. The transition bound therefore yields
\[
\operatorname{cw}(G)\ \ge\
\log_2|S_N|-r-2\log_2(K_N-1)-2-\log_2 v.
\]

## Compress subcubic wire graphs to simple cubic blocks

Write \(F(v)=v/(\log_2(v+8))^2\), an increasing function for
\(v\ge0\) since \(\ln(v+8)>2\). The graph premise gives, for every connected loopless
subcubic multigraph \(G\) of order \(v\),
\[
\operatorname{cw}(G)\le2a_0\,\nu(G)+O(F(v)).
\]

First take a connected bridgeless block \(B\) of order at most \(v\). Contract adjacent vertices and delete their \(j\) internal edges whenever either endpoint has degree at most two or \(j\ge2\). The new degree is the sum of the old degrees minus \(2j\), which is at most three. Every surviving edge still has an alternate walk avoiding it, obtained by projecting an old alternate path, so bridgelessness persists. The operation decreases \(|E|-|V|\) by \(j-1\). After finitely many contractions the result is a singleton or a connected simple cubic graph \(H\). In the cubic case,
\[
|V(H)|/2=|E(H)|-|V(H)|\le |E(B)|-|V(B)|=\nu(B)-1.
\]

Each final vertex represents a binary tree of mergers. Every piece in this tree has boundary at most three, as recorded when it was formed. Expand an order of the final vertices by listing the leaves in each merger tree with the larger child first. A proper prefix of the order within a tree is the union of earlier siblings along the path to the first unlisted leaf. Each earlier sibling occurs when the next child has at most half the current size. There are at most \(\lceil\log_2 v\rceil\) such pieces, with total boundary at most \(3\lceil\log_2 v\rceil\). At a cut inside a piece, the additional edges crossing the cut are bounded by its already listed part's boundary. Complete pieces follow the final order. This proves
\[
\operatorname{cw}(B)\le2a_0\,\nu(B)+O(F(|V(B)|)+\log_2(v+1)).
\]
For a contraction to one vertex the expansion begins with width zero.

In general, removing the bridges of \(G\) gives a tree of bridgeless blocks, with singletons allowed. Give each block weight equal to its number of vertices. There is a block \(B\) whose hanging components \(F\) in \(G-V(B)\) all have at most \(v/2\) vertices. To find one, move across any incident bridge whose far side has weight greater than half. The side just left has weight less than half, and any subsequent side of weight greater than half is a strictly smaller subtree. Each \(F\) joins \(B\) by one bridge. Place its recursive order immediately before its attachment vertex in the order of \(B\). A cut has at most one partially listed \(F\), and all crossing attachment bridges meet the next unlisted vertex of \(B\), giving at most three such bridges. Thus
\[
\operatorname{cw}(G)\le\operatorname{cw}(B)+\max_F\operatorname{cw}(F)+3
\]
with the maximum zero when there are no hanging components. Cycle rank is nonnegative and adds across bridges. On each
recursion path the linear costs sum to at most \(2a_0\nu(G)\). At
depth \(j\) the current order is at most \(v/2^j\). The sum of
\(F(v/2^j)\) for \(j\le(\log_2 v)/2\) is \(O(F(v))\), by bounding
each logarithm below by \((\log_2 v)/2\) and summing the geometric
numerators. At later depths the geometric numerator sum is
\(O(\sqrt v)\) and the denominator is at least one. The residual
bridge and binary-expansion costs are \(O(\log^2(v+1))\), which fit
\(O(F(v))\). This proves the subcubic bound.

## The bound at the actual input length

A circuit of size \(s\ge\lceil c N\rceil\) already meets the claim.
For smaller circuits \(v=2s'+1=O(N)\) uniformly. The two cutwidth
bounds give
\[
 2a_0(s-n_u+1)+O(F(v))
 \ \ge\ \log_2|S_N|-r-2\log_2(K_N-1)-2-\log_2 v .
\]
Here \(n_u=N-r\) with \(r=O(\sqrt N\log N)\), and the loss in
\(\log_2|S_N|\) counts every actual input bit, including both headers
and the tail. Since \(F\) is increasing, \(F(v)=O(F(N))\).
All other errors are also \(O(F(N))\). Rearranging with the exact
\(a_0>0\) yields
\[
 C_{B_2}(L_N)\ \ge\
 \left(1+\frac1{2a_0}\right)N
       - O\left(\frac{N}{(\log_2 N)^2}\right)
\]
at every sufficiently large actual integer length. The error denotes
one fixed positive constant after enlarging it to absorb the
different finite error terms.

## Additional reduction in the Gaussian charge

Use the fixed signed crossing cover and label rule of [proof.md](proof.md).
A center has crossing degree zero. Its incident edges initially agree
with its anchor, since it lies outside the cover. Write \(x>0\) for its
magnitude. Eligible neighbors have the same sign, lie in the cover,
and have all their crossing edge labels disagreeing with their anchors.

A degree-two neighbor is eligible if its two crossing arms each
end after one, three, five, seven or nine edges. Its crossing component is an even
path and it lies at odd distance from each endpoint. Its cover
membership is stable, both crossing labels disagree, and the
center is its only same-sign neighbor.

For a negative center also allow degree-one neighbors under this
certificate. The neighbor \(v\) is the negative endpoint of an odd
path whose unique crossing arm certifies some \(E_i\) in the point
bank. The guarded switch reaches \(v\), giving cover membership
and a crossing disagreement against the positive vertex that
leaves the cover. The other same-negative neighbor of \(v\) is
in the cover by \(P_j=S_1E_j\). The certificate selects the fixed
96 ordered index pairs \((i,j)\) listed in the numerical source.
Distinct choices give disjoint histories by the lengths, endpoint
guards and retained internal guards of \(E_{\rm pt}\).

For a center with exactly two eligible neighbors, assign the opposite of
the common anchor sign to its two edges to them. For each degree-one
target, do the same to its other same-sign edge. Call either assignment
a flip, even if the edge already has the assigned value. All targets
then cost \(\lambda\). A target of degree two originally cost
\(1+\delta\), with its only same-sign edge to the uncovered center.
A target of degree one originally cost at least one, because its
crossing edge disagreed and its edge to the center agreed. Each
additional same-sign edge at a degree-one target may increase a
nontarget charge by at most \(\delta\). Its other endpoint is already
in the cover and had at least one crossing disagreement. Unanimously
disagreeing stars already have their edge flipped.

These changes run simultaneously. A target has a unique same-sign
neighbor outside the cover, its center. Thus targets for distinct
centers are distinct. Secondary flips join two cover vertices.
Such a center receives just its two flips in the pair adjustment. All
flips on shared edges request the same value. Previously saved vertices remain saved.
The center's new charge is \(1+\delta\). If \(k\) of its two targets
have crossing degree two, the net charge decreases by at least
\[
 \beta_k=1-2\lambda+(2k-3)\delta,\qquad k=0,1,2.
\]
This calculation debits at most one secondary edge for each
degree-one target, even if a neighbor receives more than one flip.
The three values of \(\beta_k\) are positive.

### Centers with three covered neighbors

For a negative center whose smallest-magnitude neighbor is an
eligible target, use the three-edge adjustment if the other two
neighbors each have either a covered-arm certificate or their own
target certificate. Assign the opposite anchor sign to all three
center edges and each degree-one target's secondary edge.
At two eligible neighbors this saves at least \(1-\lambda\)
beyond the pair debit. The center changes from \(1+\delta\) to
\(\lambda\), and the last neighbor adds at most \(\delta\).
With one eligible neighbor the target originally paid at least one, the
center pays \(\lambda\), and the total debit at non-targets
is at most \(3\delta\). This gives at least \(\beta_0\).

For any center with three eligible targets, flip all its incident
edges and the secondary edges at degree-one targets. Write \(t\)
for the number of degree-two targets. The decrease is at least
\[
 3-4\lambda+(2t-3)\delta.
\]
The sum of the pair debits for all three unordered pairs is
\(3-6\lambda+(4t-9)\delta\), since each degree-two target appears
in two pairs. Their difference is at least \(2\lambda\).
Thus each pair can be counted even when a center has all three,
and the remaining decrease is enough for the following two debits.

For a certified three-covered negative center with smallest
neighbor eligible, count \(\beta_0\) once and \(\alpha=8963129/20000000\)
for each eligible neighbor other than the smallest. The rational
checks \(\beta_0+\alpha\le1-\lambda\) and
\(\beta_0+2\alpha\le2\lambda\) make these debits valid for
centers with two and three eligible neighbors, respectively.
For one eligible neighbor the debit is just \(\beta_0\).

Target uniqueness still holds. Every additional edge from a center
goes to a covered vertex, so it cannot modify the charge of a
different center. Its covered endpoint had a crossing disagreement,
and the same penalty debit applies. Shared requests on edges agree.
The approximate wave uses the actual smallest neighbor magnitude.
The unrounded witness predicates depend on finitely many signs
and pairwise comparisons.

### Directed calculation for the centers

Put \(r=\sqrt2\), \(m=rv-x/2\). The degree-two target uses
\(D(B_i,B_j)(x,v)\) for \(i,j\in\{0,2,4,6,8\}\), separately and then summed.
Its support has \(v<x/(2r)\). For a negative degree-one target,
a selected pair \((i,j)\) with arm \(E_i\) and sibling guard
\(Q=P_j\) has probability
\[
 \int_{d\ge \max(0,-2m)}
  g(d+m) E_i(v,d) Q(v,2m+d)\,{\rm d}d
\]
for \(x,v>0\). The two child subtrees are conditionally independent.
The guard \(E_i\) starts with the positive vertex of magnitude \(d\)
and certifies each switch along the odd crossing path. The integration
domain and the restriction on \((x,v)\) are convex linear constraints.
Each individual integrand is log-concave by the branch argument in
the graph proof, and marginalization preserves log-concavity.

At each grid pair \(x,v\), the replay uses a lower bound on \(m\)
to move the integration start inward. It uses the upper bound
on \(m\) in the decreasing Gaussian density and brackets the
sibling argument in \(Q\). It evaluates directed quarter chords
of the two factors \(gE_i\) and \(Q\), including the known Gaussian
curvature factor, with the same method as the point endpoints.
It stores a distinct table for each pair of guards.

For one center sign and two of its neighbor positions, let their
magnitudes be \(y,z\). The third magnitude is \(2rx-y-z\) when all
signs match. On the positive magnitude region the joint Gaussian
density for \((x,y,z)\) is
\[
 \frac{\sqrt{108}}{(2\pi)^{3/2}}\,
 \exp[-(x^2+(3y-2rx)^2+3(y+2z-2rx)^2)/2].
\]
This is the joint star law of the graph proof with the two neighbor
signs changed to match the center sign. We use corner density bounds
supported on \(2rx>y+z\), \(y,z>0\). The calculation brackets \(r\)
and takes the larger absolute endpoint for each squared affine term.
The two target subtrees are independent conditional on the star. The
selected neighbor magnitudes are integrated up to \(12/13\). The
center magnitude is integrated up to \(33/26\) for degree-one target
factors and \(35/26\) for degree-two target factors. Outside its
integration aperture each factor contributes zero.
Each is integrated by the directed corner mean. At a fixed cell
its two-variable mean is the eight-corner mean with repeated
coordinates. For a degree-one target add the individual means,
preserving the disjoint guard choices. Multiply the Gaussian
cell mean floor and the two target mean floors, then the cell volume.
Both assignments of the degrees contribute for a mixed pair.

The three unordered position pairs are counted with multiplicity
when all three neighbors qualify, so multiply each one-pair
integral by three. For degree-two pairs the two possible center
signs give equal contributions. The other degree combinations
use only negative centers. This gives lower expectations
\(\tau_0,\tau_1,\tau_2\) of the three pair counts.

For the additional negative-center event, use the same grid and
density, with the target at position \(y\). Enforce \(y\le z\)
and \(y\le2rx-y-z\) at the corners, using the lower bound on
\(2rx\). These linear conditions make \(y\) the smallest
of the three magnitudes. The second and third magnitudes \(z\)
and \(2rx-y-z\) each use a covered-neighbor guard. Use the
sixteen fixed point histories \(P_j=S_1E_j\) with indices
\(0,1,2,3,4,5,6,7,9,10,15,16,44,45,46,47\), and the 25 \(D(B_i,B_j)\)
for \(i,j\in\{0,2,4,6,8\}\). These are disjoint by path histories and
crossing degree. The corner means for the second guard repeat
the unused coordinate. For the third guard the replay brackets
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
\(J=\sum_{k=0}^2\beta_k\tau_k+\beta_0\tau_3+\alpha\tau_4\),
with exact rational value as checked by the replay. The finite
witnesses discard centers near crossing-degree-three defects as
in the graph proof, with the exclusion radius enlarged to contain
the new guard histories.

## Uniform rate of the graph estimate

Let \(v\) be the order of the physical cubic graph and
\(R=\lfloor(\log_2 v)/16\rfloor\). Instead of the geometric rows for
qualitative convergence, use the following radius-\(R\) rows on the
cubic tree. If a column is at distance \(d\) from its row center, its
unnormalized entry is \(2^{-d/2}f(d/R)\), where
\[
 f(t)=t^2(1-t)^2\quad(0\le t\le1), \qquad f(t)=0\quad\text{otherwise}.
\]
Normalize each row to Euclidean norm one. Its squared norm before
normalization is \(\Theta(R)\), since each positive-distance shell
has \(3\cdot2^{d-1}\) entries and \(f\) is bounded below on
\([1/4,3/4]\).

For a row \(h\) centered at a tree vertex, write \(Ah\) for its
adjacency sum over columns. At distance \(d\ge1\) the unnormalized
row residual \((A-2\sqrt2 I)h\) equals \(2^{-d/2}\sqrt2\) times the
second difference of \(f\) at step \(1/R\). This second difference
is \(O(R^{-2})\), because the extended \(f\) has a Lipschitz first
derivative. At distance zero the residual is
\(3\cdot2^{-1/2}f(1/R)=O(R^{-2})\). Summing the squared residual
over the \(O(R)\) nonzero shells and dividing by the squared
row norm gives
\[
 \|(A-2\sqrt2 I)h\|_2=O(R^{-2})
\]
for the normalized row. Summing rows whose centers neighbor a
given center is exactly the same adjacency sum. For each fixed
column the two sums have the same near and far distance counts.

For the covariance \(c_j\) between the normalized rows at distance
\(j\), the residual bound and Cauchy–Schwarz give
\(3c_1=2\sqrt2+O(R^{-2})\) and, for fixed \(j\ge1\),
\(c_{j-1}+2c_{j+1}=2\sqrt2 c_j+O(R^{-2})\).
With \(c_0=1\), induction at each fixed distance gives
\(c_j=r_j+O_j(R^{-2})\), where
\(r_j=(1+j/3)2^{-j/2}\) obeys the exact recurrence and initial values.
Add independent Gaussian tie noise of amplitude \(R^{-2}\) at each
center before variance normalization. Its effects on fixed covariance
entries are \(O(R^{-4})\), and each star residual still has
variance \(O(R^{-4})\).

Choose a finite cover of girth greater than \(4(R+D)+5\), where
\(D\ge2\) is a fixed radius containing every finite witness and
defect exclusion region. The row calculations and their inner
products on the needed neighborhoods are the tree ones. Cover
descent and entrywise covariance domination in the graph proof apply
to these rows and noise for every \(R\). They use locality and the
fiber comparison, independently of the specific radial weights.

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

### A single physical order with the same error

Each independent coordinate determining the physical keys can
affect at most \(O(2^R)\) keys. Projection of the cover fibers
preserves radius-\(R\) locality. The anchor Lipschitz bound and
Efron–Stein inequality give variance \(O(v2^{2R})\) for the
minimum charge at a fixed threshold and for the number of keys
in a fixed interval. Use \(O(R^2)\) standard-normal quantile bins
of mass at most \(R^{-2}\). With tolerance \(v/R^2\), the union
of the Chebyshev failure probabilities for endpoints and bin
counts is \(O(2^{2R}R^6/v)=o(1)\). There is a realization
with every endpoint minimum at most \(a_0v+O(v/R^2)\)
and \(O(v/R^2)\) keys in each bin. Its keys are distinct almost
surely by the tie noise. Lipschitzness then gives the same
bound at all prefixes of one order.

The explicit grouping error in the unit-factor interpolation is
\(O(v/L+vL/J_g+v^{3/4}+1)\), with \(L=\lceil v^{1/16}\rceil\)
and \(J_g=L^2\). The first terms cover cycle-rank and separator
errors and interpolation mesh increments. The third is the
deviation tolerance for the group priorities. Median-edge
ordering adds bounded cost. These errors are \(O(v^{15/16})\)
uniformly, giving the simple cubic cutwidth loss
\(O(v/(\log_2(v+8))^2)\) together with the physical-key error.

## Strict declared floor

The checked limiting coefficient is greater than 4.8954101.
Using any fixed intermediate coefficient strictly between
4.895410 and that limit absorbs the uniform, sublinear transfer
loss for all sufficiently large actual input lengths. Thus the
declared eventual lower bound is 4.895410 N with zero explicit
remainder. The language and verifier are the original fixed
objects above.
