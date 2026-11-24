# Graph evidence A: marked-tree source, line optimization and directed gain

Technical graph evidence. The auxiliary all-prefix bounds below provide context; the scored proof uses only the specified one-cut threshold and its supplied directed witnesses. Read with proof.md and graph-b.md.


---

## Evidence A1: crossing-lines construction

# A crossing-lines Gaussian field and a uniform cubic cutwidth bound

**Auxiliary graph statement.** Every simple cubic graph H with a perfect matching, uniformly over the graph and the matching, has an original-vertex layout of cutwidth at most
```
(73765053306/549755813888 + o(1)) |V(H)| < (0.13418 + o(1)) |V(H)|.
```
In particular no girth or local weak limit assumption is imposed on the simple, triangle-free, 2-connected cubic kernels in the Williams reduction. The numerator comes from an outward upper enclosure, not a floating-point sample. There is also a genuine limitation of **this specified field and one-stage pinned matching policy**: the limiting mean density of splits on an ideal straight line at zero threshold is at least 72537121230/549755813888 > .130648. This is a lower bound for the method, not the cutwidth of an actual finite graph.

## A. Exact cycle objective; comparison on an actual finite graph

Write M for a fixed perfect matching, and C for its complementary 2-factor. Contract M to the loopless 4-regular multigraph Q. At each contracted vertex its four half-edges come in two pairs, one per original endpoint; continuation within a physical C-cycle means taking the other half-edge of the same pair. Each vertex on a C-cycle has one matching bit s prescribed by its incident M key. Minimize the sum of split indicators (NAE among three incident edge bits) over the C-edge bits. In an individual C-cycle the minimum is exactly the maximum matching size of the subgraph of disagreement edges i--i+1 with s_i different from s_{i+1}. Indeed, the split vertices cover those edges. Given a vertex cover, force each unsplit vertex outside it to have both C-edge bits equal to its own s. The forced bits are consistent; fill remaining edges arbitrarily. The disagreement graph consists of paths, or an alternating even cycle. Each has minimum vertex cover equal to maximum matching, by pairing consecutive vertices on each path or around the even cycle. The analogous cost F_L on L consecutive anchor positions with free line-end incidences is the maximum matching on the L-1 path edges. It is symmetric submodular in the L anchor bits, since it is a partial minimization of a sum of NAE terms (endpoints have one incident path edge; equivalently free dangling bits), using the lattice argument of the earlier general edge-folding theorem.

For every binary cycle word s of length k (including k<L), the deterministic upper bound
```
opt_cycle(s) <= (1/(L-1)) sum_{i=0}^{k-1} F_L(s_i,...,s_{i+L-1 mod k})      (L>=2)
```
holds with periodic indexing throughout. Repeat a cycle maximum matching on the two-way periodic line; restrict to each window's L-1 edges. Each restricted matching has size at most its window maximum. Each original selected edge is counted exactly L-1 times in the k windows, including multiplicities on windings. No independence of window costs is assumed.

The universal cover of Q is a 4-regular tree with the pairing of half-edges lifted from Q. Each pair of straight directions belongs to a bi-infinite line. Each anchor vertex belongs to exactly two lines. Any marked physical C walk lifts as a straight, nonbacktracking path of distinct anchors, even on short odd cycles or when the base contains parallel connectors: at each step the next labeled half-edge is distinct from the reverse. A bounded, invariant, unit Gaussian coefficient row on this marked tree defines a single unit physical Gaussian score per M edge by fiberwise Euclidean norm compression to Q (the prior universal folding theorem applies with anchor vertices indexing the iid noises). Every physical window may have repeated matching-edge keys. Choose its straight lift; its physical covariance dominates every entry of the chosen lifted covariance by fiberwise Cauchy--Schwarz. The same physical rows serve all windows. Since F_L is symmetric submodular, Gaussian comparison and the zero-threshold maximum from the prior theorem bound E F_L at **every** physical threshold by its tree zero-threshold mean. We do not compare probabilities of individual alternating sign patterns after folding.

## B. An invariant field with prescribed along-line correlations

Let h_0=1 and h_d=a t^(d-1) for d>=1, where t=41/57, a=28/57, so t^2+2a^2=1 and 2 sum_{d>=1} h_d^2=1. First truncate to d<=K, rescaling the positive h_1,...,h_K by their common factor that makes 2 sum h_d^2=1, and set h_d=0 for d>K. For a row rooted at u and coordinate x on the marked tree, segment the u-x path into maximal straight runs of lengths d_1,...,d_j. Give x coefficient product h_{d_i} when j<=J (empty product 1) and zero otherwise. Normalizing by sqrt(1+2J) produces a unit row: the squared mass of the first run across either line is 2, and every next run must change to the other line, providing factor 2 sum h_d^2=1. The support radius is at most JK. The definition uses neither orientations nor arbitrary choices at crossings.

Take roots distance d>0 on the same line, parametrized by 0,d. Project every tree vertex onto that line at an integer b. From the projection the remaining path takes zero or more alternating runs, the first off-line run having only one choice of line. For a fixed projection, the squared-products mass of such remainders with j runs is exactly 1 for every j>=0. A nonzero prefix to the projection consumes one run; for both rows to be nonzero, at least one of the roots differs from the projection. Thus the allowed remainder has 0 through J-1 runs. Both root prefixes contribute h_|b| h_|b-d|, with h_0=1. The exact off-diagonal inner product is therefore
```
J/(1+2J) sum_{b in Z} h_|b| h_|b-d|.
```
The diagonal is 1. Send J->infinity and then K->infinity for each fixed d, using convergence in l2 of the truncated sequences. In the resulting stationary Gaussian line law the correlations are r_0=1 and
```
r_d = a t^(d-1) + t^d/2 + (d-1) a^2 t^(d-2)/2  (d>=1; omit last term for d=1).
```
This follows by separating b<0, b=0,d, 0<b<d and b>d. Its finite covariance matrices are positive definite: they are Gram matrices of distinct integer translates of a nonzero absolutely summable sequence h on Z, divided by 2. A finite linear dependence would give a Laurent polynomial times the nonzero continuous Fourier series of h identically zero; the latter is nonzero on an interval, so the polynomial is zero. The finite truncations on the actual graph need not be positive definite.

## C. Simultaneous physical layout, finite radius, and ties

To turn all the threshold estimates into one physical order choose an integer S and split each C-cycle longer than S by pre-scoring ceil(k/S) predetermined edges into segments of length at most S vertices; leave shorter cycles intact. Score a cut edge by the matching key at a predetermined endpoint. There are at most 2n/S such edges. At any threshold forcing each such edge to its prescribed bit increases the previous global optimum by at most 2 per forced edge: start with an unconstrained optimum and change that edge, affecting at most its two endpoints. Hence the total expected split cost at every threshold is at most
```
n E_tree[F_L at zero]/(L-1) + 4n/S
```
for the actual bounded row (J,K). Here independence between blocks is neither asserted nor needed. All remaining unscored edge components are bounded paths or cycles of at most S vertices. Internal free-bit minimizers exist; their least choices increase monotonically with the exterior threshold bits, by jointly submodular NAE and the meet/join inequality. This applies to cycles as well as paths, with the 2-factor edges itself being optimized; the comparison in A is what justifies the bound for cyclic cores. Identical M and cut scores at two incidences are permitted.

Each free edge flips once at an exterior score. Each M key occurs at two physical vertices, supplies at most four cut-edge scores (two incident per endpoint), and is relevant at most to the two blocks of those endpoints plus two blocks for each incident cut edge. Thus O(S) physical edges flip at any one source-key event, counting multiplicity; this refers to the number a given M key can trigger, not the count of keys tied numerically. Use n equal-Gaussian-measure windows. The number of M keys in a window has mean <=1 and its indicators have bounded dependency degree for fixed J,K. Each block cost <=S depends on a ball of radius bounded in terms of J,K,S. Dependency-colored Hoeffding plus a union bound over windows and thresholds supplies simultaneous deviations O_{J,K,S}(sqrt(n log(2n))). At a grid cut, the number of open vertices is the optimized split count. At a partial tie/window the newly touched vertices are at most twice O(S) times the number of M keys in that same window; the bound therefore includes all actual edge-prefix cuts. Order physical vertices by middle incident edge as in the predecessor: earlier/later crossing edges inject into the bag of vertices with inclusive first/last interval containing the current median, with at most one extra crossing edge. This establishes a physical cutwidth bound, not just a balanced or mean cut.

For fixed L,S choose K then J sufficiently large to make the zero-threshold expectation for this L close to its infinite-law limit. It is a bounded cost continuous off finitely many coordinate hyperplanes, each of probability zero. For these fixed parameters let n grow so that the concentrated remainder divided by n vanishes. Then let L,S grow. This order of choices proves that the limiting upper coefficient is at most the infinite-line limiting maximum matching density. It asserts no small practical radius at moderate n.

## D. A uniform upper and an intrinsic obstruction from sign moments

Let P_j be the probability that j consecutive adjacent sign transitions all alternate in the limiting line law. The left-to-right greedy matching of disagreement edges on a path is maximum: when the first available edge is present it can replace any edge using its right endpoint. Write D_i for the indicator of disagreement at edge i and B_i for its greedy selection. Then B_i=D_i(1-B_{i-1}), initialized by B_0=0. Repeated substitution yields the Bonferroni-type bounds for all interior edges i>=5:
```
L4=P_1-P_2+P_3-P_4 <= E B_i <= U5=P_1-P_2+P_3-P_4+P_5 .
```
The remainder has the next alternating sign and involves products of bits times a 0/1 selection indicator. The finitely many leading edges cost O(1/L) on averaging. The limit exists, for example by subadditivity of the path maximum after cutting one central edge with additive discrepancy at most 1, or by the standard tiling argument with fixed long windows. Thus U5 is an upper bound valid for physical layouts by A-C. L4 bounds the ideal-tree path mean from below. It is not a lower bound on the actual folded finite-graph cost.

The certificate gives
```
scale = 549755813888 = 2^39
U5 in [73520857720,73765053306]/scale
L4 in [72537121230,72678012835]/scale.
```
For each j, the sign indicator of either alternating word equals the complement-symmetric sign expansion, so
```
P_j = 2^(-j) sum_{I even subset {0,...,j}} (-1)^(sum_{i in I}i) E(prod_{i in I} sign X_i).
```
The checker derives the coefficients of L4 and U5 directly by enumerating subsets, grouping by translation and reversal. Only two-, four- and six-spin moments remain. To evaluate an even moment on a positive definite correlation matrix A, scale every correlation in its first row/column (except the diagonal) by x from 0 to 1; this is a convex combination of positive definite matrices. The starting moment vanishes because the first centered sign is independent. Gaussian integration by parts (or differentiating a Gaussian density at its two moving signs using smooth sign approximation) gives
```
M(A) = (2/pi) sum_{j>0} A_0j integral_0^1
         M(C_j(x))/sqrt(1-x^2 A_0j^2) dx .
```
Here C_j(x) is the correlation matrix on the other coordinates conditional on coordinates 0 and j being zero. The product of sign derivatives contributes four times the bivariate density at zero, hence 2/pi. Gaussian mollification justifies the formula by dominated convergence since all matrices in each fixed integration segment are positive definite on the compact interval. For two signs M=2 asin(A_01)/pi. To form C, first remove coordinate j and then 0. For any remaining i,k the unnormalized conditional covariance is
```
A_ik - A_ij A_kj - [x^2/(1-x^2 A_0j^2)] (A_0i-A_0j A_ij)(A_0k-A_0j A_kj),
```
and normalization uses these positive diagonal variances. The checker uses rectangles of width 1/384 for each integral, with **whole-cell intervals**, including the outer cell as an interval when evaluating the inner integral. No quadrature accuracy assumption or statistical samples enter the certificate. See the checked source for integer-directed rounding, Machin pi, square root and alternating Taylor enclosures.

## Boundary of the advance

The ratio 41/57 selects this particular explicit line filter. No optimality over all line processes is claimed. On this filter and policy, infinitely extending the path window cannot lower the tree mean below the stated L4 barrier. To beat .130648 one must change at least the pinned-edge choices, the field, the policy, or exploit graph-specific folding slack. A finite-profile numerical optimizer, exterior-record novelty, formal proof and independent external review are not claimed.


---

## Evidence A2: two-pass construction and analytic bounds

# A matched-pair second pass beats 0.1306476 on every matched cubic graph

**Auxiliary computer-assisted real-arithmetic theorem.** For every finite simple cubic graph H with a perfect matching, uniformly as its order n grows, the original-vertex cutwidth satisfies cw(H) <= (5603815728133241/42949672960000000+o(1))n; the coefficient is 0.13047400228989406 to the displayed precision. In particular it clears 0.1306476; the strict slack therefore also allows fixed finite choices of every approximation parameter for the looser 0.1306476 coefficient. This applies to all connected simple triangle-free 2-connected cubic kernels in the Williams reduction, without a girth, expansion or local-limit hypothesis. This is a graph ordering theorem and only implies a circuit bound through the hypotheses of the chosen kernel-to-circuit transfer. No assertion of independent public priority, peer review, Lean verification, or leaderboard acceptance is made.

## 1. One bounded-radius local policy on the cover suffices

Fix a perfect matching M; its complement consists of cycles C. Contract each matched pair to an anchor in the loopless 4-regular multigraph Q with four distinct half-edge incidences, paired into two straight continuations (one original physical vertex per continuation). Equal-neighbor parallel edges are not identified as incidences. The marked universal cover is a 4-regular tree and its straight walks are double rays. We will construct a finite-radius equivariant atomless scalar key for each straight-line edge there.

Here is why this gives a finite H order, not just a cut expectation on a large-girth subsequence. For four exterior C-edge threshold bits around a matching pair use the **submodular** cost
```
F(e1,e2,e3,e4)=min_q {NAE(e1,e2,q)+NAE(e3,e4,q)}.
```
NAE is the binary not-all-equal indicator; submodularity survives joint partial minimization by meet and join. The already proved conditional-quantile fiber transfer (proof.md, Evidence P1) puts these *final nonlinear scalar edge keys* simultaneously on all physical edges, preserving each marginal, their bounded local dependence and the per-pair upper bound for expectation of F at each threshold. This does not apply submodularity to a composition on the initial matching bits; that different statement is false. A single lifted star handles every physical matched pair including coincident projected endpoints.

Choose the least minimizing matching bit at each physical pair. It is monotone in the exterior bits and switches only at one of the four exterior keys. Thus every edge has a scalar key; sorting gives an edge order. With n marginal quantile bins for exterior keys (their common law is atomless), concentration sums of split indicators at the grid thresholds and exterior-edge-bin counts: at fixed radius R each depends on independent base noise in a bounded ball, so dependency graphs have maximum degree O_R(1). Hoeffding after a dependency coloring (or bounded-differences on independent uniforms) and a union bound give one outcome with all fluctuations O_R(sqrt(n log(2n))). Each grid bin contains O_R(sqrt(n log(2n))) exterior edges and hence at most three times as many edges of either type counting all tied entries, since an exterior edge touches at most two pairs. Comparing a prefix within a bin with its preceding grid threshold changes the inclusive split bag by at most twice its bin edge count, so arbitrary internal ties cost only this error. The bag at edge position i consists of physical vertices whose interval from first through last incident edge contains i. Ordering physical vertices by their *median (second) incident edge position* and breaking ties arbitrarily converts the maximum bag size into a cutwidth bound: at any prefix let i be the median-position boundary. A crossing edge whose position is below i charges its endpoint on the later side (only its first edge can be below i there); one above i charges its earlier endpoint (only its last edge can be above i there). These endpoints lie in bag i and each receives at most one charge. The unique edge at position i can be counted separately. This handles a cut dividing a median tie as well. Consequently a per-physical-vertex cover expectation bound c for **every** threshold gives cw(H)<=cn+O_R(sqrt(n log(2n)))+1. No efficient evaluation of arbitrary conditional quantiles is claimed.

## 2. The limiting marked Gaussian law and its finite realizing scores

Reuse the finite crossing-lines construction (Evidence A1 above). Set t=41/57, a=28/57, s=sqrt(1-t^2)=sqrt(2)a, delta=t-a. Define h_0=1 and h_d=a t^(d-1) for d>=1. Truncate h to [1,K], rescaling its positive terms to have twice their squared sum equal to one. The score at an anchor is a linear row on iid standard normals at cover vertices: coefficient proportional to the product of the h-lengths of maximal straight runs on the unique path from the anchor to each vertex, supported up to J runs; divide by sqrt(1+2J). The empty path uses 1. Squared mass at each nonzero run depth is 2 (first line has two choices, later line is forced). For each finite collection of positions, let J and then K grow. Along one line the limiting correlations are r_0=1 and
```
r_d = a t^(d-1) + t^d/2 + (d-1)a^2 t^(d-2)/2, d>=1
```
(omit the last term at d=1). Equivalently r_d=(1/2) sum_i h_|i| h_|i-d|. For positions i,j away from the shared root on **different** lines the covariance is
```
h_|i| r_|j| + h_|j| r_|i| - h_|i| h_|j|.
```
Indeed project each noise vertex onto the union of the two infinite lines. If its projection k is on the i-line away from root, the two prefix factors are h_|i-k| and h_|j| h_|k|. Any nonempty remaining path must first leave that line; each number of further runs has squared weight mass one. The available remainder depths are J+O(1) independently of k (with the run-truncated coefficients fixed); dividing by 1+2J gives 1/2 in the limit. The j-line projections give the swapped contribution. The root has no direction leaving the union and contributes only O(1/J). Sum, using the along-line formula, then remove run truncation in l2. Thus all cross/colliding-source statistics used below come from **one** compatible limit of finite local score rows. Along lines the limiting sequence also has the absolutely summable moving-average representation X_i=sum h_|i-j| W_j/sqrt(2) with iid standard Gaussian W_j. This gives ergodicity by the iid shift factor.

## 3. True straight-line minimization, and why the finite factor attains its limit

At threshold T the bit on a pinned anchor is b_i=1[X_i>T]. Minimize total NAE splits over exterior edge bits on each line. A path or cycle has minimum equal to the maximum matching on its disagreement edges i--i+1 (equivalent vertex cover: unsplit anchor vertices force both exterior edges, consistently; the disagreement runs are paths or alternating even cycles). On the double ray the conditional dynamic-programming message through a bit b, from incoming difference m in {-1,0,1}, is
```
f(m,1)=1[m!=-1],      f(m,0)=-1[m!=1].
```
This follows from subtracting the two fixed-exiting-bit optimal path losses, with difference cost(0)-cost(1). Two equal successive bits erase any entering message. With incoming messages from both sides of an edge, that edge's least minimizing bit is one precisely when their sum is positive. These updates and threshold edge decisions are monotone in all inputs. A finite equivariant edge key uses the j nearest anchors on either side of that edge and boundary message 0 on each outer side, then the strict-positive sum rule. It is endpoint-reversal invariant. As T varies its monotonicity and homogeneous ends give a single switch that is an input key, hence atomless.

In the stationary limit, each message stabilizes if the ray ever sees an identical adjacent pair. The probability of failure within an even length 2k is bounded (uniformly in T) by the zero-threshold probability of length-2k alternation: each directed alternating orthant under common shift is logconcave in the shift by Gaussian logconcavity under marginalization (Prékopa); it is even by reversal combined with central reflection. This probability decreases to zero, since the stationary Gaussian line has positive density of identical adjacent signs at threshold zero by the ergodic theorem. Thus finite free-boundary keys approximate the exact infinite messages uniformly over thresholds, and the four-edge F costs differ in mean by o(1). The same argument with bounded path endpoints ignored shows the expectation of first-pass splits at a vertex equals the stationary matching density c1(T): a deterministic path of N anchors has minimizing edge bits from the two path-end messages (ties choose the coordinatewise least global minimizer), and boundary influence at a uniform interior vertex tends to zero. This gives c1(T)=lim_{N->infty} E F_N/N, where F_N is the symmetric submodular min-split path cost.

Each fixed finite window then has the limiting Gaussian law described in §2 as J,K increase. Bounded threshold-cost expectation converges at each T by Gaussian convergence (each threshold hyperplane has zero marginal probability). Convergence is uniform in T by a finite marginal-quantile grid: fixed dimension, and changing a threshold across a marginal bin changes cost only if some coordinate belongs to that bin. Thus for each epsilon fix one finite window and finite run/length cutoffs, the same for all H. Apply §1 and let epsilon decrease. At no point is a nonlocal argmin or a global common random tie coin pushed to the finite graph.

## 4. Exact analytic first-pass profile including tails

Let A_k(T) be the probability of k consecutive disagreement edges (k+1 alternating anchor bits). Left-greedy matching on a finite path has accepted-edge indicator equal to the alternating sum of consecutive disagreement-run indicators ending there; the odd truncations are upper bounds, regardless of unknown preceding bits. Averaging or sending path length to infinity gives c1(T)<=U7=sum_{k=1}^7 (-1)^(k-1) A_k(T). Also c1 is even and nonincreasing on [0,infty): each symmetric submodular F_N has that Gaussian threshold property, as follows (condition at changing coordinate X_i=T, reflect the others about their conditional means <=T, apply diminishing returns and complement symmetry). Pass to the limit. We do **not** claim the final two-pass cost is monotone in |T|.

At T=0 the existing exact integer Plackett sign-moment checker for U5 re-run at grids 2048,1024 gives U5(0)<=73670242301/2^39. For general T, A1-A2=S2(T)=P(X_0>T,X_2<=T), by the two-disagreement XOR identity and stationarity. For k=3,5 each alternating word has even length, and reflection/reversal and Gaussian logconcavity as above yield A_k(T)<=A_k(0). For any centrally symmetric union of words of length L, Gaussian change of density on translation by T*1 gives its mass at T at least exp(-D_L*T^2/2) times its mass at zero, where D_L=1' C_L^{-1} 1 for the line covariance. The opposite exponentials combine into a cosh >=1. In particular this lower-bounds A4(T). Writing kappa=sqrt((1-r_2)/(1+r_2)), direct bivariate differentiation gives
```
-S2'(T) = phi(T) [2 Phi(kappa T)-1]  (T>=0).
```
Integration and monotonic density imply S2(0)-S2(T)>=kappa phi(T)phi(kappa T)*T^2. A3(0)<.019 by the four-sign expansion, A4<=A3 and D5<1.34; the exact rational check, using exp(-x)>=1-x, yields U5(T)<=U5(0)-.046 T^2 for |T|<=.28.

The U7 improvement is A6-A7, an eight-position word (+,-,+,-,+,-,+,+) and its complement. Its zero-threshold mass has a short independent variational certificate. Put K=diag(sign) C_8^{-1} diag(sign), choose independent halfnormals Y_i of scales v_i=(457,308,331,334,332,303,474,651)/1000 in the positive orthant, and let E=E[Y' K Y]. Then the probability of both words is at least
```
(prod v_i)/(128 sqrt(det C_8)) * exp((8-E)/2).
```
This is Jensen for the Gaussian density divided by the product halfnormal density; E=sum_i Kii vi^2 +(2/pi) sum_{i!=j} Kij vi vj, since E|N(0,1)|=sqrt(2/pi). Using 1+(8-E)/2 in place of the exponential gives >.000346 by exact rational inequalities. D8<1.63 and central symmetry thus give A6(T)-A7(T)>.00032 uniformly for |T|<=.275. In summary,
```
c1(T) <= 73670242301/2^39 -.046 T^2-.00032   (|T|<=.275).
```
Beyond this range bound *first-pass cost* by its value at .275 using monotonicity of c1; second-pass cost never exceeds first-pass cost pointwise.

## 5. The actual pairwise gain, calculated on the tree

At one endpoint of a root anchor, let x,y be the two incoming messages from its line rays, and b the root pinned bit. The two adjacent first-pass edge bits are [x+f(y,b)>0] and [y+f(x,b)>0]. Flipping the root matching bit gains D=+1 if both bits oppose b, D=-1 if both equal b, D=0 if mixed. Reoptimizing each matching anchor independently therefore saves max(0,D1+D2) per matched pair, even after all other matching anchors change (pairs partition the physical vertices). For b=1, write Y for both x,y<=0 and at least one equal -1, and N for at least one equal -1. Then Y={D=+1}, N={D>=0}; both are decreasing in the incoming messages. The positive gain equals 1[Y1 N2]+1[Y2 N1]. The two lines through the matched pair are exchangeable, so dividing by two physical vertices per pair gives the probability of the single event b=1,Y1,N2 as a lower bound on expected gain per physical vertex. Gains for b=0 are nonnegative and are omitted entirely.

Here is a low-dimensional realization that makes the joint event certifiable. Generate independent standard normals R and U. R is the root score. For each of the two lines ell=+1,-1 take an independent V of variance 1/2. Its two rays initially have coordinates
```
D = a(R+ell U)/2 +/- a V,        E=((a+t)R+delta ell U)/2.
```
At each ray step generate a fresh N(0,1/2) innovation W, independently across rays/steps. The emitted score and next coordinates are
```
X = D+E-delta W,
D' = tD+2a^2 W,                  E' = tE+aD.
```
The root plus j scores on each of four rays have precisely the joint Gaussian law from §2. For transparency this can be derived from one-sided filters B=s sum_{k>=1}t^(k-1)W_k of standard Gaussian source coordinates along a ray. The four root B modes are iid standard normals; the root score is their half-sum. On one whole line write scores as convolution by h/sqrt(2). Its central source W_0 is the other line's two B modes summed over sqrt(2); symmetrically that other line's W_0 is this line's two modes summed over sqrt(2). The two within-line marginal laws are correct, the shared root agrees, and expansion gives the cross covariance in §2. For one current ray mode B, the next mode B' and current standard source satisfy B'=tB+sZ, W_current=sB-tZ with independent standard Z. Let L be the contribution of the opposite half-line and central source in units before the 1/sqrt(2) line normalization; L'=a W_current+tL. With D=aB and E=(L+tB/sqrt(2))/sqrt(2) this gives the displayed update for W=Z/sqrt(2). For the first ray site L=(a*(other-line pair sum)+t*(opposite same-line ray))/sqrt(2). Writing the given line's pair as (R+ell U)/2 +/- V and the other pair's sum as R-ell U gives the stated initialization. Independently of this derivation, the exact checker verifies every covariance on the root and four rays out to depth 8, more than the numerical proof uses.

For a lower bound use six sites of each ray and the largest boundary message +1. This dominates the actual incoming message by monotonicity. Conditional on D,E, if p_m and p_p denote probabilities of the messages -1 and +1, the update is the integral over X<=T of (1-child p_p) and over X>T of (1-child p_m), respectively. Conditional on R,U,V the rays are independent. For adjacent ray marginals m1,p1 and m2,p2, the probabilities of Y,N are
```
m1(1-p2)+m2(1-p1)-m1m2,         m1+m2-m1m2.
```
These are decreasing in stochastically larger messages. Integrate each conditional line outcome over its own V; conditional on R,U the two lines are now independent with opposite U. Finally integrate the product Y(R,U) N(R,-U) over R>T,U. This handles all cross-line dependence, not an independent-endpoint guess.

## 6. Directed finite sums and all-threshold check

The source cert_ray.cpp contains no floating arithmetic. Coordinates and all input grid widths are integers in units 1/40000. The finite 6-step convolution bounds apply throughout each closed threshold band [l,h]: use the more restrictive root event R>h, and use threshold l on all ray scores, making each ray message larger. For every parent exact state D,E, increasing either coordinate increases all future scores with the same innovations, because both transition coefficients are nonnegative. The integrator rounds D and E **up**, treats an above-table coordinate as adversarial incoming message +1, and clamps a below-table argument to its lower grid boundary (an upward change). In an innovation bin it rounds the *next* D using the upper W endpoint, and integrates the current bit exactly to interval precision at W=(D+E-l)/delta; E' is independent of W. Unknown tails are counted as adversarial where a child message is needed. The recurrence stores integer lower bounds for p_m and upper bounds for p_p in units 2^-29. It computes the latter from a lower bound on its complement: the full zero-bit probability plus the one-bit integral of child p_m. This avoids subtracting uncertain mass in the wrong direction.

For the first line at R,U query each ray at the upward D endpoint in its V interval; since its two coefficients have opposite sign, bound them at opposite ends. Use V symmetry, omit its tails and sum positive Gaussian bin masses times lower probabilities; the second line is analogous. After V integration both line success functions decrease with R and with their own contrast ell U, since initial D,E coefficients of these variables are positive. In each R,U cell take the higher R and the higher own contrast separately for each line, multiply the resulting two lower bounds, and omit tails. All products use explicit integer rounding. In the Y formula first decrease the true m_i to their lower bounds with the true p_i held fixed (partial derivative 1-p_other-m_other >=0), then increase p_i to its upper bound or 1-m_i_lower, whichever is smaller. This remains a valid lower bound. No independence is asserted across physical cycles or across private uniforms added by folding.

The Gaussian CDF routine uses scale 2^50, midpoint cells of width 1/32768 up to 8, and exact integer brackets for the density at each midpoint. At positive argument x the latter brackets exp(-x^2/q), q=1 or 2, by 11th/12th alternating Taylor sums at x^2/(256q)<=1/4 followed by eight squarings. The two normalizations are bracketed by consecutive 15-place decimals verified by rational Machin bounds. The density has first derivative bounded by 1 and second by 2 for both variances, so cumulative midpoint error is at most 8/(12*32768^2), and the partial cell error at most 1/32768^2. Reflection handles negative arguments; above 8 the omitted normal tail is <2^-29 by the standard Mills inequality. Bin masses use lower CDF differences (upper CDF at the left), and every dropped rounding residue is unfavorable. All p accumulators are nonnegative and <=2^58; Gaussian fixed-point intermediate products use signed 128-bit integers. The supplied range of parameters fits the stated integer types. Compiler OpenMP parallelism changes only the order of exact integer additions; it is optional.

The general band argument explains how all-prefix bounds could be certified. For the scored fixed threshold set l=h=-14000/40000 and use the directed root-above lower column in the accompanying transcript. For the scored threshold, replay.py verifies the embedded witness hashes and replays the relevant exact arithmetic; replay.py --full regenerates the directed transcripts.
