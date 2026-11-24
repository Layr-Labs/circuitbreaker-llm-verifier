# A fair local cubic cut and a full-basis circuit lower bound above five
Submission proof for mathematical review. Only the joint theorem in bound.json is scored. The graph proof is computer-assisted; graph-a.md and graph-b.md supply its dependencies and replay.py reconstructs the arithmetic witnesses.

## 1. Fixed family and eligibility: an elementary norm predicate

The previous submission invoked a general sumset extractor. The signed-cut argument below needs only a sufficiently dense accepting set with no large accepting rectangles. Use the classical norm-graph idea of Kollár, Rónyai and Szabó, *Norm-graphs and bipartite Turán numbers*, Combinatorica 16 (1996), 399–406, DOI 10.1007/bf01261323. Norm-one languages also precede this work in the CircuitBreaker competition (including Arthur Ostapenko's November 2026 entries). No novelty of that family is claimed. We reproduce all the required family properties, using a weaker elementary common-neighbor bound; no external norm-graph theorem, extractor, correlation breaker, sampler, or irreducible-polynomial algorithm is assumed. The original graph-law and numerical work is unchanged.

### 1.1 A fixed uniform evaluator at every actual input length

For n<64 define f_n=0, including n=0. For every integer n>=64 put

```
L=floor(log2 n), B=ceil(log2(L+1)), a=L-B, p=2^a,
ell=floor(n/p), m=p*ell, d=2^floor(a/2), k=m/d, Q=2^k.
```

Here a>=3, d>=2 and d divides m. From 2^L<=n<2^(L+1) we get 2^B<=ell<2^(B+1)<4(L+1). Further, 0<=n-m<p<=n/(L+1), and p>=n/(4(L+1)). Thus m/n tends to 1, d is of order sqrt(n/log n), and k and d*log2 d are O(sqrt(n log n))=o(n).

Construct a field F of size 2^m with an ordered binary basis by the following fixed deterministic procedure. Enumerate monic binary polynomials of degree ell in order of their coefficient encodings; choose the first having no monic divisor of degree 1 through floor(ell/2), by polynomial division. This is exactly the irreducibility test. Existence: in a splitting field the 2^ell distinct roots of X^(2^ell)-X (derivative -1) form a field by closure under addition, multiplication and inversion. An element of degree j<ell over F2 generates a proper subfield, with j dividing ell by the dimension formula and hence j<=ell/2. Such elements belong to the roots of X^(2^j)-X. Their union has cardinality at most sum_{j=1}^{floor(ell/2)} 2^j<2^ell, so an element of degree ell exists. The case ell=1 is immediate. Take the quotient by the selected irreducible polynomial with its ordered polynomial basis. Since ell=O(log n), exhaustive enumeration and division require n^{O(1)} bit operations.

Repeat the following extension exactly a times. In the current field K of size 2^r define T(z)=sum_{j=0}^{r-1} z^(2^j). This is binary-linear and takes values in F2 since T(z)^2=T(z), using z^(2^r)=z. Its defining polynomial is nonzero of degree 2^(r-1)<2^r, so it is not identically zero. Choose the first basis vector u with T(u)=1. The map z->z^2+z has kernel {0,1} and image of size 2^(r-1); its image has trace zero by telescoping, and is therefore precisely the trace-zero space. Hence X^2+X+u has no root and is irreducible. Adjoin a root beta using the quotient, retaining the old basis followed by beta times the old basis. The final binary dimension is ell*2^a=m. Addition is XOR; a multiplication at one extension level uses at most five multiplications and finitely many additions in the preceding level (reduce beta^2=beta+u). The depth a is O(log n); trace evaluations and a search of at most m basis vectors per level use polynomially many field operations. All field operations and the entire field construction are uniform polynomial time. No advice, random choices, oracle, or large field search is used.

Map the first m bits of x to z in F in the resulting ordered binary basis; ignore the remaining n-m bits. Define

```
e=(2^m-1)/(Q-1)=1+Q+...+Q^(d-1),
f_n(x)=1 iff z^e=1 in F.
```

The exponent is a positive integer with at most m bits, so repeated squaring evaluates the predicate in polynomial time. All small lengths, field choices, and rounding are fixed by this algorithm. The language {x:f_|x|(x)=1} belongs to P and hence NP with the single verifier V(x,w)=1 iff w is empty and f_|x|(x)=1, and witness bound W(n)=0. Both directions of the verifier equivalence follow from the definition.

### 1.2 Accepting density and essential input count

Write g for the restriction to the m active bits. By the finite group order, every nonzero field element is a root of X^(2^m-1)-1; there are exactly 2^m-1 such roots and they are distinct. Because e divides 2^m-1, the polynomial X^e-1 divides that split squarefree polynomial and has exactly e roots in F. Here e is odd since Q is even. The uniform accepting density delta=e/2^m is at least 2^(-k), since e>=Q^(d-1)=2^(m-k). Each of the m active bits is essential: otherwise the involution flipping that bit pairs all accepting inputs without fixed points, contradicting the odd count e. An arbitrary circuit for f_n can be restricted by fixing its tail bits to zero and normalized without increasing its gate count, as in §2. The restricted circuit for g therefore reads all b=m active bits. All asymptotic conclusions still use the original actual length n, with b/n tending to one.

### 1.3 Elementary common-neighbor lemma

Let K be any field and a_{ij} in K for i,j=1,...,d, pairwise distinct down each fixed column j. The system

```
product_{j=1}^d (z_j-a_{ij}) = 1   (i=1,...,d)
```

has at most d^d solutions in K^d. The classical sharper d! bound is not required. Here is a direct argument for the weaker bound.

In the polynomial ring on d pairs (X_j,Y_j) let L_{ij}=X_j-a_{ij}Y_j, W=product_j Y_j, and I be the ideal generated by H_i=product_j L_{ij}-W. All H_i have degree one in each pair. For each column j, the d forms Q_{ij}=product_{r != i} L_{rj} form a basis of homogeneous degree d-1 in that pair: evaluation at (a_{ij},1) is zero on every Q_{rj} except Q_{ij}, which is nonzero. The space has dimension d. For any degree E>=d-1, products of these forms with monomials of degree E-(d-1) span homogeneous degree E, because every degree-E monomial factors as monomials of those two degrees.

For any D>=d, every polynomial P homogeneous of degree D in each pair is congruent modulo I to W times a polynomial homogeneous of degree D-1 in each pair. Process columns j=1,...,d. At the start of column j, P is congruent to (product_{h<j}Y_h)R, where R has degree D-1 in processed pairs and D in the rest. Terms divisible by Y_j already have the next desired factor. Each remaining term contains X_j^D. Expand its coefficient in the other d-1 pairs using the spanning products above; the degrees there are at least D-1>=d-1. In each resulting summand one Q_{r,h} is selected in each other column h. Each excludes only its own row r, so a row i among the d rows is excluded by none of the d-1 selections. That summand contains L_{ih} for all h!=j. Write

```
X_j^D = L_{ij} X_j^(D-1) + a_{ij}Y_j X_j^(D-1).
```

The second term is divisible by Y_j. The first term has the full product of the L_{ih}, including column j; replacing that product by W modulo I supplies Y_j as well. The old outside factors Y_h are retained, and each pair degree is unchanged. This proves the induction, and after column d yields the claimed common factor W.

Dehomogenize by setting every Y_j=1. Let J be the ideal generated by the original affine equations with their right sides subtracted. Any polynomial of maximum individual degree D>=d can first be homogenized separately to degree D in all pairs. The preceding argument and dehomogenization give a representative modulo J of maximum individual degree D-1. Iterate. The quotient ring is therefore spanned over K by the d^d monomials with every individual exponent in {0,...,d-1}. For any finite collection of s distinct solutions, evaluation induces a surjective linear map from that quotient onto K^s: for a chosen point and each other point choose a differing coordinate, use the linear factor in that coordinate vanishing at the other point, normalize it to one at the chosen point, and multiply the factors. These polynomials give all unit vectors on the finite collection. Therefore s<=d^d. This applies to every finite subset of the solution set and proves the lemma, with no algebraic-geometry or multiplicity assumption.

### 1.4 All coordinate partitions have no large accepting rectangles

Let A,B be subsets of F such that (u+v)^e=1 for every u in A and v in B. If |A|>=d choose distinct u_1,...,u_d. For any v in B, repeated binary Frobenius gives

```
product_{j=0}^{d-1} (v^(Q^j)+u_i^(Q^j))=1   (i=1,...,d).
```

In each column j the u_i^(Q^j) are distinct, because Frobenius is injective on a field. Apply the lemma with z_j=v^(Q^j) and a_{ij}=-u_i^(Q^j). Since the coordinate with j=0 equals v, different v give different solution tuples. Thus |B|<=d^d. Consequently for K_0=d^d+1 it is impossible that both sides have size at least K_0. For any partition of the m binary coordinates, identify each side assignment with its element of F supported on that side of the fixed basis. Each map is injective, and their field sum is the full input. Every accepting monochromatic rectangle is therefore thin on at least one side (cardinality <K_0). This conclusion is uniform over all coordinate partitions, without any independence or balance requirement on the zero-inputs.

We have log2 K_0=O(d log d)=o(n) and log2(1/delta)<=k=o(n). These, together with b=m=n-o(n), are exactly the hypotheses of the one-sided signed-cut inequality in §4 and Evidence P2. The quantitative entropy loss is sublinear at every sufficiently large integer length; no sparse-length extrapolation is involved.

## 2. Full B2 normalization and physical graph
All 16 truth tables on two ordered pins, repeated pins, arbitrary acyclic topology and unrestricted fanout allowed. First restrict the ignored tail to zero; the constant normalization below never increases the gate count. Prune nodes not reaching the output. Constants/projections/unary gates are allowed as truth tables. If free constant sources arise under substitutions, evaluate constant-only gates, substitute constants into consumer truth tables, fill an unused pin by repeating its other input; no new counted gate is required. Raw variable or constant output cannot compute this hard family at large n. For the remaining pruned circuit of t gates and b raw inputs, a signal used f>=2 times is expanded into a path of f-1 copy vertices. Copies represent wires, not gates. There are 2t uses and t+b-1 nonoutput signals, hence t+1-b copies. The undirected wire multigraph is connected, loopless, subcubic with N=2t+1 vertices, M=3t+1-b edges and cycle rank mu=t-b+1. Every genuine internal gate and copy has degree three, raw inputs degree one, output degree two; parallel edges and repeated pins retain multiplicity.

Credit McGuire, Schlesinger and Williams (MSW), public proof.md §§2.3–2.4 in rrwilliams0/circuit-lb-4.6155n (observed blob 63de7d5bf662378bf5921c031f5d70a0cffde4ba), solely for the wire-graph representation and the block/kernel compression. Their Li-family counting argument and all-prefix bound are not premises of this theorem; the norm-one family and accepting-only transcript count are proved in §1, §4 and Evidence P2. In a bridgeless block contract adjacent vertices if one has degree <=2 or there is a parallel pair; otherwise contract any triangle. Each step reduces vertex count, preserves bridgelessness, cannot increase rank, and has a resulting vertex of degree <=3. The terminal graph is a singleton or simple triangle-free cubic. A cut vertex in a bridgeless cubic graph would have at least two edges into each of two components, impossible. Thus every nonsingleton kernel is 2-connected and has h<=2(mu(B)-1). The perfect-matching consequence of Tutte’s theorem follows directly: for every nonempty vertex set S, each odd component of H−S has odd edge boundary of size at least three by cubicity and bridgelessness. The sum of those boundaries is at most 3|S|. For S empty each component of a cubic graph has even order. This verifies Tutte’s odd-component condition. The binary/ternary contraction forest pieces each have <=3 original block boundary edges. Bridge blocks have tree incidence and additive rank. No MSW numerical Gaussian certificate is imported.

## 3. Exact graph input
Graph-a.md supplies the marked-tree score, first-pass cost and directed sign-moment/gain bounds. Graph-b.md supplies fair patches, eleven-pin deletion, the recalibration at T=-7/20, and directed finite probability proofs. The conditional-quantile transfer below folds final physical vertex bits, not the potentially nonsubmodular intermediate costs. Every finite simple cubic graph H with a perfect matching admits a law of vertex bits with fair marginals, uniformly bounded dependence and expected crossing count <=c|V(H)|, c=12405759/100000000. Fix finite shell, window and patch radii using strict margins independent of H. Efficiency of generating the cut is not needed. Numerical gate: first-pass c1<.125645569; gain >=851968/2^29; motif saving >=300311352/2^48; final relaxed majority <=33422386/67108864<.499. Exact arithmetic gives c1-gain-motif<c and starting majority >.5054287. The score-bearing check is `cone-check.py` via `replay.py:run_cone`, deriving `62405759/12405759` from Evidence B6. The frozen `fair_shift35_20261124/EXACT.stdout` instead records the earlier cost `12423221/100000000` and coefficient `62423221/12423221` (approximately 5.024721125); replaying that predecessor is a calibration/regression check, not certification of the stronger coefficient. Evidence B5–B6 explains the distinct certificate roles. Continuous activation yields exact fair marginals. This is a one-cut theorem, not a cutwidth theorem at c.

## 4. All-lengths circuit transfer

For the restricted normalized circuit from §1, b=m. Place W=2*1_raw+indegree-outdegree on each physical vertex; sum 2b and bounded magnitudes. On side A with q inputs, r outgoing and s incoming crossing arcs, w=r+s, the directional surpluses q-r and b-q-s have sum b-w and difference W(A)-b. Simulate actual topological circuit evaluation by sending each crossing arc bit and a final output bit. Complete transcripts define monochromatic rectangles by induction on local checks. Only accepting transcripts need to be thin. For fixed opposite-side assignment, at most 2^(r+1) transcripts arise by varying the first side. The uniform mass of accepting transcripts thin on that side is at most K_0*2^(r+1-q); symmetrically the mass thin on the second side is at most K_0*2^(s+1-(b-q)). Their union includes every accepting input, of uniform mass delta. Therefore

```
delta <= 2*K_0*(2^(r-q)+2^(s-(b-q))).
min(q-r,b-q-s) <= log2(4*K_0/delta).
w+|W(A)-b| >= b-2*log2(4*K_0/delta) = b-o(n).
```

The detailed block/forest lifting in Evidence P3 transfers the unchanged fair kernel laws to the actual restricted physical graph. For fixed cutoff D, pieces have O(D) mass, charged overhead O(N/D), and total kernel order <=2mu(G). Expected cut cost is at most 2c mu(G)+O(N/D), with signed charge sum b+o(N) with probability tending to one by bounded dependence and Chebyshev. Conditional expectation of the nonnegative cost on that event is at most its unconditional expectation divided by 1-o(1); some outcome attains both. Consider any sequence of original computing circuits with bounded original gate count divided by n; unbounded ratios are automatic. The restricted gate count t is no greater than the original count, and N=2t+1=O(n). Fix D, let n grow, then let D grow. The rectangle inequality gives n-o(n)<=2c(t-b+1)+o(n) in this iterated-limit sense, since b/n->1. If original gate count <(1+1/(2c)-epsilon)n infinitely often for fixed epsilon>0, the restricted counts along that subsequence contradict the inequality. Hence liminf over all integer actual lengths C_B2(f_n)/n >= 1+1/(2c)=62405759/12405759. No gate-value independence, direct sum, restricted basis, or finite-n exact endpoint assertion is used.

## 5. Frozen proofs
The following sections give the conditional-quantile transfer, signed cut inequality and weighted block/forest lifting in full. Only the one-cut argument in §§1–4 is scored. The auxiliary all-prefix graph estimates retained for historical context in graph-a.md are not used as premises; no all-threshold circuit comparison is asserted here.


---

## Evidence P1: conditional-quantile transfer

# Conditional-quantile transfer of nonlinear local scores from graph covers

**Conditional-quantile transfer theorem.** A bounded-radius equivariant scalar score depending arbitrarily on iid noise on a universal cover can be replaced by one scalar score on the base. All one-score laws are retained. For every common (or separate) threshold, every finite tuple of lifted scores, and every submodular cost of the resulting bits, the base expectation is no greater. This applies to nonlinear edge keys, not just Gaussian linear rows. The key is to compress the **final scalar score**, not to declare a nonlinear composition of Boolean costs submodular.

## 1. Precise theorem and construction

Let pi:T->G be the universal cover of a finite connected graph; bounded maximum degree is needed only for uniform layout asymptotics below. A finite or countable regular cover with transitive deck action on each fiber works as well. Lift any supplied deterministic base markings. Vertex noises Z_z on T are iid with a standard Borel law. The score indices may be vertices, edges, or a specified equivariant countable index family over finitely many base orbits. For each lifted index a there is a real-valued Borel score Y_a depending only on a finite set N_a of iid noises, equivariant under deck transformations (up to a common null set). Assume that N_a is likewise chosen equivariantly or enlarged to a bounded ball. Then there exist base scores X_i indexed by the deck orbits, Borel functions of independent standard uniforms U_x on base vertices, with:

1. X_[a] uses only U_x for x in pi(N_a) and has exactly the distribution of Y_a.
2. For arbitrary lifts a_1,...,a_k, thresholds t_1,...,t_k and any submodular function F on {0,1}^k,
```
E F(1[X_[a1] <= t1],...,1[X_[ak] <= tk])
 <= E F(1[Y_a1 <= t1],...,1[Y_ak <= tk]).
```
A uniform choice of strict inequalities instead is also valid. The *same* base scores work for all tuples, all thresholds and all costs. F need not be nonnegative, normalized, monotone or complement-symmetric. No covariance matrix, linearity or Gaussian interpolation is assumed. Repeated base indices need not be separated.

Fix any order of the finitely many vertices of G. At step x let the remaining outside variables be the uniforms already introduced and all as yet unprocessed lifted noises over base vertices other than x. They are mutually independent of the untouched fiber at x and of the new U_x. For every lifted row a replace its current scalar value by **the increasing conditional quantile**, evaluated at the same independent U_x, of its current law over the unprocessed fiber pi^-1(x), with these outside variables held fixed.

More explicitly, for a rational q integrate the indicator that the current value is <=q against the product of the original iid noise laws in this fiber; each row uses only finitely many of those variables. The resulting CDF values are Borel in its outside arguments. For 0<u<1 take the infimum of rationals q where the conditional CDF is >=u. Assign arbitrary real values on the null endpoints u=0,1. This standard quantile rule has the conditional law, including its atoms. If the row does not use the current fiber, the quantile is its unchanged value for almost all u. Processing never expands its support: none of the Borel operations just described can introduce a variable not already used by that row other than the new U_x for the replaced fiber. One can choose bounded-radius scores literally as orbit translates of a representative, or use Fubini at the finitely many steps and the countably many rows/relevant translate equalities to disregard versions on null sets.

The row indexed by a deck translate is still the corresponding translate in the unprocessed coordinates. In fact a deck transformation permutes the original independent coordinates *within each fiber*, preserves the law integrated at x, and fixes the previously introduced base uniforms pointwise. Thus the same conditional quantile recipe is equivariant in the remaining arguments. After all base vertices have been processed there are no remaining lifted coordinates: all rows in the same orbit are equal almost surely. Selecting one representative defines X. Order-independence or canonicity with respect to every automorphism of G is not asserted or needed.

## 2. Why a fiber replacement lowers every submodular bit cost

Condition on the outside variables at one step. The old threshold bits are arbitrarily coupled Bernoulli variables of fixed marginal probabilities p_i. The new ones have these same conditional marginal probabilities, but their events are initial intervals of a *single* uniform U_x, hence nested. With atoms and '<=' this holds up to endpoints of uniform measure zero; the analogous strict-threshold event has the corresponding strict CDF mass. The entire real-score quantile is used, so threshold-by-threshold compatibility is automatic.

Here is an elementary proof of the classical minimum-coupling fact. Relabel the coordinates so p_1>=...>=p_k and let S_j={1,...,j}. Put w_j=F(S_j)-F(S_{j-1}). For any A, adding its members in this order and using diminishing returns gives F(A)>=F(empty)+sum_{j in A}w_j. Therefore the expected cost in any coupling is >=F(empty)+sum_j p_j w_j. In the nested coupling only the prefix sets S_j occur, where these inequalities are equalities. Ties between marginal probabilities do not change the value. This is the usual Lovasz-extension/comonotonic coupling statement; it is not a new inequality. Taking the outer expectation proves a decrease at that step. Iteration gives the displayed theorem. Exact preservation of each single-score law follows by the conditional quantile construction at each step. The conclusion holds for all thresholds as a deterministic inequality between expectations; it does not assert that a random outcome simultaneously realizes infinitely many strict comparisons.

## 3. Specialization to the scored finite-radius fair cut

Apply the proved transfer to the final physical vertex scores (Bernoulli bits), not to a nonlinear intermediate optimization cost. Index the rows by physical vertices, using the two labeled incidence pairs per contracted anchor as in graph-b.md. Apply the inequality of §1 to the two lifted endpoint scores of each edge with the submodular XOR cost. The sum of expected crossing indicators cannot increase; no atomless marginal or threshold sweep is needed. All one-row marginals, in particular exact fairness, are preserved. Every row still depends only on the independent base marks in its projected finite-radius support. A fixed degree bound and fixed radius give a uniform bound on the number of overlapping supports, hence a dependency graph of bounded maximum degree. This is the only specialization needed in Evidence B5 and the weighted cut assembly in P3. It does not assert a cutwidth bound or a polynomial-time algorithm for conditional integrals.

## 4. Prior work, controls and assurance

The Gaussian linear common-edge-row fold is recovered transparently: each conditional fiber contribution is a normal with conditional mean the outside sum and standard deviation the l2 norm of the fiber coefficients. A common standard Gaussian quantile produces the norm compression. The nonlinear graph-cover theorem here uses conditional **quantiles**, not a covariance bound on the input of a Boolean circuit.

For classical context compare Francis Bach, *Submodular functions: from discrete to continuous domains*, arXiv:1511.00394, on continuous/distributional submodular extensions; Chernozhukov, Fernandez-Val and Galichon, *Improving point and interval estimates of monotone functions by rearrangement* (Biometrika, 2009), on sequential coordinate rearrangement. Neither the classical nested-Bernoulli coupling fact nor conditional quantile sampling is attributed here as new. Comparison to the exact graph-cover statement is source-relative; exhaustive mathematical novelty remains unverified.

The all-size statement is the written argument, not an inference from finite fixtures.


---

## Evidence P2: one-sided signed cut inequality

# A single signed cut using accepting transcripts only

Let a pruned acyclic Boolean circuit on b primary input bits and t binary gates have the wire graph of §2: each primary input has one outgoing arc, each genuine nonoutput gate indegree two and outdegree one, each copy indegree one and outdegree two, and the gate output indegree two. There are N=2t+1 vertices and cycle rank mu=t-b+1. A cut partitions vertices and counts crossing arcs with multiplicity. The target on the b input bits has accepting density delta>0 under uniform input and no all-ones rectangle under any coordinate partition with both side cardinalities at least K_0. Section 1 proves these hypotheses for the present target, b=m, with log2(K_0/delta)=o(n).

Put W(v)=2*1_input(v)+indeg(v)-outdeg(v) on each vertex. Thus input and ordinary gate weights are +1, copy weights -1, and output weight +2; total weight is 2b. For a vertex side A carrying q raw inputs and r outgoing, s incoming arcs, w=r+s, cancellation of internal arcs gives

```
(q-r)+(b-q-s)=b-w,
(q-r)-(b-q-s)=W(A)-b.
```

Two parties receive their respective raw input assignments. Evaluate the DAG in public topological order and send every crossing arc bit when available, followed by the output bit. For any fixed complete transcript, prescribing its values turns every party's local consistency checks into conditions on only its own assignment (induction in the fixed topological order). Thus each nonempty transcript class is a product P x Q, and its output bit is fixed.

Fix the second party's assignment. Varying the first party's q bits gives at most 2^(r+1) complete transcripts, because the first party's r outgoing bits and, if necessary, the final output bit determine the remaining messages together with the fixed second-party input. Each accepting transcript with |P|<K_0 contains fewer than K_0 such local assignments. Summing over all fixed second-party assignments gives uniform mass at most K_0*2^(r+1-q) of accepting transcripts thin on the first side. Symmetrically the accepting mass thin on the second side is at most K_0*2^(s+1-(b-q)). Every accepting transcript is thin on some side by the hypothesis; rejecting transcripts need not satisfy any rectangle restriction. Therefore

```
delta <= 2*K_0*(2^{-(q-r)}+2^{-(b-q-s)})
      <= 4*K_0*2^{-min(q-r,b-q-s)}.
```

It follows, using the two identities above, that every cut satisfies

```
w(A)+|W(A)-b| >= b-2*log2(4*K_0/delta).
```

The output bit is explicitly counted, and parallel arcs and repeated pins retain multiplicity. Only the existence of one low-cost charge-balanced cut is needed. This is the general one-sided form of the previously submitted signed-cut argument. Under the kernel law stated in §3, Evidence P3 constructs such a cut and proves the coefficient. No two-source extractor is used.


---

## Evidence P3: weighted block and forest assembly

# From one kernel cut to the actual graph

G is any connected loopless subcubic multigraph on N vertices. Use the bridge blocks and contraction forest constructed in §2 of this proof, following MSW §2.4. For each bridgeless block B (singletons allowed), every merge has 2 or 3 children and every piece has at most 3 B-edges leaving it; the final contracted graph is a singleton or a simple triangle-free 2-connected cubic kernel H_B with h_B≤2(μ(B)−1). The sum of μ(B) equals μ(G). Select an integer d≥2.

1. Root the bridge-block tree anywhere. Process postorder. A block receives the unfinalized vertices from each child; each such child remainder has <d vertices and its only connection to this block is its parent bridge. Bundle each retained child remainder with its actual attachment vertex in the block. Each such starting bundle has ≤1+3(d−1) original vertices. Finalize this job, consisting of the block plus all these remainders, if its mass is ≥d or it is the root. Otherwise return it unfinalized. Each nonroot finalization incurs at most one severed bridge to its parent, even when that parent is later finalized jointly with other remainders. Finalized jobs are disjoint; except possibly the root each has ≥d vertices, so ≤N/d such bridges suffice. Every nonsevered bridge inside a job belongs to a starting bundle. Every job has one anchor block; the unused child blocks within bundles are intentionally treated as monochrome here.
2. In the contraction forest of each anchor block, start with its vertex bundles from step 1, successively merging the still-attached original vertices of the children. If the remaining mass of a piece reaches 4d, detach this *residual bundle* as its own part and return zero attached mass to the parent. Starting bundles have mass <4d, so every detached residual bundle has <12d vertices. The sum of the detached sizes is ≤N: at most N/(4d) detachments. Charge at most 3 edges *outside the full corresponding contraction piece* in the anchor block on each detachment. Edges from a new residual bundle to an earlier detached descendant, although potentially numerous in total over all descendants, were already charged when the innermost previously detached endpoint was separated. More explicitly, for an edge with endpoints in differently colored bundles, if at least one endpoint was detached then take its chronologically earliest detached endpoint. If its other endpoint was inside the same full piece at that moment it was either detached earlier or belongs to that same residual bundle; both contradict the choice/separation. Thus the edge is among the at most three originally leaving that full piece. Attachments in step 1 follow their block vertex together and add no new uncharged edge. At termination each vertex of the kernel or singleton carries <4d original vertices; zero masses allowed. Any remaining edge between different final bundles projects to an edge of the contracted kernel (or it has already been charged). Thus total overhead in physical edges is ≤N/d+3N/(4d), before any cuts inside individual bundles.
3. **Nonnested bisection corollary.** Give G arbitrary deterministic weights in [−2,2], with total 2b. Color each detached small part by its own fair independent sign; color each final singleton independently the same way. On each nonsingleton kernel independently use the fair sign rule of graph-b.md, Evidence B5 (all vertices carrying their full bundles). The expected cut cost is ≤(c+η)·2μ(G)+C_η(N/d+1)+7N/(4d). The number of kernels used is at most the number of jobs, and their sizes sum to at most 2μ(G). Each independent piece or kernel spin governs at most 12d original vertices, so its total absolute weight ≤24d; the number of terms (including empty kernel carriers) is O(N), and the dependency graph has uniformly bounded degree D_η. The weighted side mass therefore has expectation b and variance O_{η,d}(N). Chebyshev gives imbalance o(N) with probability 1−o(1) for fixed d,η. The nonnegative cut cost conditional on this good event has mean at most its unconditioned expectation divided by 1−o(1); some outcome attains this bound. No variance of the cut cost is needed. For N=O(n), first let n→∞, then d→∞ at fixed η, then η→0. Evidence P2's one-sided law, with log2(4*K_0/delta)=o(n), forces b≤2c μ(G)+o(n) in the liminf sense. If t/n is unbounded along a subsequence the claimed lower bound there is automatic. This proves the circuit conclusion under the stated kernel hypothesis, which is required only for the compressed kernels.
