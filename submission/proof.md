# An eventual coefficient of 4.8957965 for the norm-one language

For the total norm-one language defined below there is an integer \(N_0\) such that at every actual input length \(N\ge N_0\), every exact, single-output circuit over all sixteen binary truth tables has at least \((9791593/2000000)N\) gates. This includes arbitrary fan-out and repeated pins. A deterministic polynomial-time evaluator establishes uniform NP membership with the empty witness.

The construction, graph layout, probability witnesses and quantitative transfer are due to the predecessor and the earlier work credited below. The contribution here is the exact eventual floor extracted from their unrounded estimate. All ingredients used in the inference are supplied here. The finite integral formulas and their directed evaluation details appear in the numerical appendix. [Hrushikesh's source at revision a49ccfa](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/a49ccfaf5eb17002ea7edfb90c9fe76bcc9e2c82/submission) · [Numerical appendix](certificate.md)

The graph layout uses the exact constants
\[
a_0=\frac{415782528257175539331450637327469}{3239608320000000000000000000000000},\qquad c=1+1/(2a_0)=\frac{2035586688257175539331450637327469}{415782528257175539331450637327469}.
\]
For connected simple cubic graphs of order \(v\), the graph proof and its uniform-rate estimate below give
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

For a vertex and three incident edge labels, count the labels that disagree
with the anchor at that vertex. For respectively zero, one, two and three
disagreements, charge \(0,1,1+\delta,\lambda\), with
\(\lambda=47091/100000\) and \(\delta=2102599/1000000000\).
Let \(U=13064756648/10^{11}\) bound the expected crossing-cover density
and \(X=8113009/50000000\) bound the crossing-edge density.
The replay proves lower bounds \(s_1=1712234861992131/553779200000000000\) and
\(s_2=2436373455716264877/1842977177600000000000\) for disjoint populations of
saved vertices of crossing degree one and two on a single cover.
A further local adjustment at initially isolated anchors gives a
charge rebate \(J=78136407780990579173253934703/1619804160000000000000000000000000\), proved in the section on additional Gaussian charge reduction.
The penalized anchor objective and the unit-factor interpolator below give
a uniform simple cubic cutwidth constant
\[
 a_0=(1-\delta)U+\delta X-(1-\lambda-2\delta)s_1-(1-\lambda)s_2-J
   =\frac{415782528257175539331450637327469}{3239608320000000000000000000000000}
.
\]
The quantitative graph estimates below sharpen the
uniform width loss to \(O(n/(\log_2 n)^2)\) and transfer it at the
exact coefficient \(1+1/(2a_0)\).

The four component lower bounds, in positive and negative order within each crossing degree, are
\[
\frac{6132802254699}{3276800000000000},\frac{168947820237}{138444800000000},\frac{647713169760441}{838860800000000000},\frac{3958389147471}{7199129600000000}.
\]

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
[PR 905](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/40f0a2292a512fb3cf23df9010c66741acc1bb19/submission). The predecessor used here extends
its ray aperture and refine its center integration. The exact
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
bound \(X\) above. Apply the additional isolated-center adjustment in the section on additional Gaussian charge reduction
after this labeling. It reduces the expectation of the displayed upper
bound by at least \(J\), yielding \(a_0\) on the limiting local model.

## Certification of the numerical inequalities

For a marked path the cubic-tree wave gives the children of a vertex with parent and current values \(a,b\) as \(\sqrt2 b-a/2\pm\eta\), with \(\eta\) a fresh centered normal of variance \(1/12\). Conditional on the incident fields, different forward subtrees have independent contrasts. A directed alternating crossing path therefore has a one-dimensional conditional transition in each successive magnitude. Truncate the alternating sum of half-counts of nonbacktracking walks at 25 edges. An \(\ell\)-edge path contributes exactly \(\lceil\ell/2\rceil\) for \(\ell<25\) and at least that amount otherwise. Even cycles are also covered by this upper count. Gaussian tail bounds outside the aperture and the signed tensor-integration errors on it give an upper bound of \(0.130647566474459387249570528851780<U\). The appendix derives the signs, normalizations and every error term.

For the savings, fix the marked root sign, crossing degree and the finite branch history witnessing the switch. The outside guard and endpoint tests are simultaneous tests on the same root configuration. The positive and negative cases of the two crossing degrees give the four floors stated above. Conditional on the local Gaussian values the independent forward contrasts give products of the indicated branch weights. Histories are separated before summing wherever they are incompatible. A discarded admissible event costs no lower mass. On each rectangular panel the nonnegative itinerary integral is bounded below by its volume times the geometric mean of its corner bounds, using the log-concavity and marginalization argument in the appendix. The quadratic Gaussian factor admits the stated Jensen curvature correction. All floor divisions are directed down, and interval endpoints for a domain are moved inward. These bounds, the four-root counting formulas, the positive arm calculation and the radial upper calculation are detailed in the numerical appendix. The two declared replays check the disjoint parts with exact arithmetic. Together they furnish precisely the inequalities used in the charge estimate. [Integration formulas and enclosures](certificate.md)

## Additional reduction in the Gaussian charge

Use the fixed signed crossing cover and label rule above.
A center has crossing degree zero. Its incident edges initially agree
with its anchor, since it lies outside the cover. Write \(x>0\) for its
magnitude. Eligible neighbors have the same sign, lie in the cover,
and have all their crossing edge labels disagreeing with their anchors.

A degree-two neighbor is eligible if its two crossing arms each
end after one, three, five, seven or nine edges. Its crossing
component is an even path and it lies at odd distance from each
endpoint. Its cover membership is stable, both crossing labels
disagree, and the center is its only same-sign neighbor.
For a negative neighbor allow a second degree-two certificate.
One crossing arm is an odd-distance guarded near history
from the 74-entry point bank \(E_{\rm pt}\). The other is one of
the twelve even-distance far histories in \(O_g\), with its next
positive certified by an outside guard.
The guarded prefix switches through the negative vertex and
through the next positive on its far arm. The two adjacent
positives both leave the cover, so both crossing labels disagree
with the negative vertex. These mixed-arm certificates are
disjoint from the even-path ones by ordered terminal parity.

For a negative center also allow degree-one neighbors under this
certificate. The neighbor \(v\) is the negative endpoint of an odd
path whose unique crossing arm certifies an \(E_i\) from the full
74-history point bank. The guarded switch reaches \(v\), giving
cover membership and a crossing disagreement against the positive
vertex that leaves the cover. Its other same-negative neighbor
is covered by some \(P_j=S_1E_j\) from the same bank. Both choices
are made among disjoint histories by lengths, endpoint guards and
retained internal guards. They inspect conditionally independent
branches given the star at the degree-one target.

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
neighbor eligible, count \(\beta_0\) once and \(\alpha=889947797/2000000000\)
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

The two positive ray patterns specified in the numerical appendix end at guarded vertices in the cover. Each admitted flip has the covered endpoint and the charge reduction specified there. The mixed and same-type arm integrals count the admissible pairs with their multiplicities. Combining the five center floors \(\tau_i\) and the two ray floors \(\sigma_i\), the extra saving is at least
\[
 J=\sum_{k=0}^2\beta_k\tau_k+\beta_0\tau_3+\alpha\tau_4+
       \beta_0\sigma_0+\beta_1\sigma_1 .
\]
Here \(\beta_k=1-2\lambda+(2k-3)\delta\) and \(\alpha=889947797/2000000000\). Direct substitution gives the rational \(J\) above. The target-uniqueness and conflict rules limit the debit at a shared endpoint. The appendix proves the ray and center counts and evaluates their directed lower masses. [Center and ray calculations](certificate.md)

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

For completeness, at fixed neighborhood radius the local relations of the limiting wave span the covariance kernel. The independent root value, two root residuals and one contrast per other internal vertex account for its whole remaining rank. A deepest internal vertex proves independence of the relations through a child coordinate. Project onto the complementary support, whose limiting covariance has a positive smallest eigenvalue. The square roots of the restricted covariance matrices differ by \(O(R^{-2})\). The complementary residual has coordinate standard deviations \(O(R^{-2})\) by the star-residual estimate. Coupling with common normals then changes each relevant sign or strict pairwise comparison with probability \(O(R^{-2})\): for a centered Gaussian form of positive limiting variance this follows from the angle between its two normalized coefficient vectors. The finite union bound covers all local witness predicates and all defects. The same proof applies to the finite alternating radial upper count. These comparisons transfer the expected central charge bound with a graph-independent \(O(R^{-2})\) density loss. The covariance-kernel calculation and the coupling estimates are given in full in the numerical appendix. [Covariance estimate](certificate.md)

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

## Strict floor at every sufficiently large length

Put \(d=9791593/2000000\). The circuit transfer above gives \(C_{B_2}(L_N)\ge cN-KN/(\log_2 N)^2\) for one constant \(K>0\) and every integer \(N\) beyond a fixed \(N_1\). Exact subtraction gives
\[
 c-d =
 \frac{83309088867973844534353753831883}
      {831565056514351078662901274654938000000}>0 .
\]
Choose \(N_0\ge\max(2,N_1)\) so that \(\log_2 N_0\ge\sqrt{K/(c-d)}\). Then the loss is at most \((c-d)N\) for every \(N\ge N_0\). Hence \(C_{B_2}(L_N)\ge dN\) at every such length, as claimed. In particular its all-length liminf divided by \(N\) is at least \(d\).
