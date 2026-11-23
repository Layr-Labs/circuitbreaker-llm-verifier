# Norm-one language and transfer to full \(B_2\)

Put \(c=244649/50000=4.89298\). The graph theorem used here is a uniform bound \(\operatorname{cw}(H)\le a_0|V(H)|+o(|V(H)|)\) for connected simple cubic graphs, with
\[
a_0=\frac{30070896090723854017009185977}{234131155797934080000000000000}
      <\frac{25000}{194649}=\frac{1}{2(c-1)}.
\]
Its proof is supplied in [proof.md](proof.md). The transfer needs only \(\operatorname{cw}(H)\le a|V(H)|+A\) for fixed \(0<a<1/[2(c-1)]\) and \(A\ge0\). The uniform graph bound gives this premise with any fixed \(a\) strictly between \(a_0\) and the threshold, by increasing \(A\) to cover the finitely many smaller orders. Cutwidth counts crossing edges at the cuts of a vertex order, including multiplicity in a multigraph.

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

The graph premise gives, for every connected loopless subcubic multigraph \(G\) of order \(v\),
\[
\operatorname{cw}(G)\le2a\,\nu(G)+O_A(\log^2(v+1)).
\]
First take a connected bridgeless block \(B\) of order at most \(v\). Contract adjacent vertices and delete their \(j\) internal edges whenever either endpoint has degree at most two or \(j\ge2\). The new degree is the sum of the old degrees minus \(2j\), which is at most three. Every surviving edge still has an alternate walk avoiding it, obtained by projecting an old alternate path, so bridgelessness persists. The operation decreases \(|E|-|V|\) by \(j-1\). After finitely many contractions the result is a singleton or a connected simple cubic graph \(H\). In the cubic case,
\[
|V(H)|/2=|E(H)|-|V(H)|\le |E(B)|-|V(B)|=\nu(B)-1.
\]

Each final vertex represents a binary tree of mergers. Every piece in this tree has boundary at most three, as recorded when it was formed. Expand an order of the final vertices by listing the leaves in each merger tree with the larger child first. A proper prefix of the order within a tree is the union of earlier siblings along the path to the first unlisted leaf. Each earlier sibling occurs when the next child has at most half the current size. There are at most \(\lceil\log_2 v\rceil\) such pieces, with total boundary at most \(3\lceil\log_2 v\rceil\). At a cut inside a piece, the additional edges crossing the cut are bounded by its already listed part's boundary. Complete pieces follow the final order. This proves
\[
\operatorname{cw}(B)\le2a\,\nu(B)+A+3\lceil\log_2 v\rceil .
\]
For a contraction to one vertex the expansion begins with width zero.

In general, removing the bridges of \(G\) gives a tree of bridgeless blocks, with singletons allowed. Give each block weight equal to its number of vertices. There is a block \(B\) whose hanging components \(F\) in \(G-V(B)\) all have at most \(v/2\) vertices. To find one, move across any incident bridge whose far side has weight greater than half. The side just left has weight less than half, and any subsequent side of weight greater than half is a strictly smaller subtree. Each \(F\) joins \(B\) by one bridge. Place its recursive order immediately before its attachment vertex in the order of \(B\). A cut has at most one partially listed \(F\), and all crossing attachment bridges meet the next unlisted vertex of \(B\), giving at most three such bridges. Thus
\[
\operatorname{cw}(G)\le\operatorname{cw}(B)+\max_F\operatorname{cw}(F)+3
\]
with the maximum zero when there are no hanging components. Cycle rank is nonnegative and adds across bridges, so the costs \(2a\nu(B)\) along any recursion path sum to at most \(2a\nu(G)\). There are \(O(\log(v+1))\) levels and at most \(A+3\lceil\log_2 v\rceil+3\) additional width per level. This proves the subcubic bound.

## The bound at the actual input length

A circuit of size \(s\ge\lceil cN\rceil\) already meets the claim. For any smaller circuit, \(v=2s'+1=O(N)\), with a constant independent of the circuit and of \(N\). Combining the two cutwidth bounds gives
\[
\begin{aligned}
2a(s-n_u+1)
 &\ge \log_2|S_N|-r-2\log_2(K_N-1)-2-\log_2 v
                  -O_A(\log^2(v+1)),\\
s &\ge \left(1+\frac{1}{2a}\right)N-O_{a,A}(\sqrt N\log N).
\end{aligned}
\]
Here \(n_u=N-r\), \(r=O(\sqrt N\log N)\), and the length estimates count the full original input, including both headers and the tail. Fix \(a,A\) before \(N\) and the circuit. The gap \(1+1/(2a)-c>0\) absorbs the entire sublinear term for every sufficiently large integer \(N\). Thus this fixed total language satisfies
\[
 C_{B_2}(L_N)\ge \frac{244649}{50000}\,N
\]
at every sufficiently large actual input length, with exponent \(1\) and zero subtracted remainder.
