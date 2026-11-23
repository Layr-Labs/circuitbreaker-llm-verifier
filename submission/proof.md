# An eventual 4.8833N gate lower bound by split-penalized softened crossing covers

For the fixed total uniform NP language in "From cubic layouts to the norm-one language", unrestricted full-$B_2$ circuits require at least $(48833/10000)N = 4.8833\,N$ gates for every sufficiently large actual length $N$, with zero explicit remainder.

On the limiting radial Gaussian wave of the infinite cubic tree, the hard crossing-cover density has certified upper bound $U = 1633094581/12500000000 = 0.13064756648$ and the crossing-edge density is $W_1 = \frac{3}{2\pi}\arcsin(1/3) < 4056504299/25000000000 = 0.16226017196$. We introduce a four-level submodular vertex star cost $(0, 1, 1+\delta, \lambda)$ charging $1+\delta$ when two incident edges disagree with the vertex anchor and $\lambda$ when all three disagree. At $\lambda = 27/50$ and $\delta = 41/20000$, this $\delta$-penalty on degree-two disagreements eliminates both secondary tie-breaking and layout inflation, yielding exact unit layout factor $R = 1$: if every prefix of a vertex ordering of a simple cubic graph $H$ has generalized anchor cost at most $M$, then $\operatorname{cw}(H) \le M + o(|H|)$ uniformly. A simultaneous odd-path switching construction using a pairwise-disjoint $(\mathrm{OUT}, V_{\mathrm{odd}})$ partition and selective noncrossing flips yields crossing-degree-one and crossing-degree-two saving densities at least
$$
s_1 = \frac{1256465834907}{419430400000000}, \qquad s_2 = \frac{345262219467}{268435456000000},
$$
certified by integer-only log-concave and Hölder-compressed quadrature. Hence the limiting cubic width and implied full-$B_2$ gate coefficient are
$$
a_0 = (1-\delta)U + \delta W_1 - (1-\lambda-2\delta)s_1 - (1-\lambda)s_2 = \frac{135009406533170243909}{1048576000000000000000},
$$
$$
1 + \frac{1}{2a_0} = \frac{659297406533170243909}{135009406533170243909} > \frac{48833}{10000} = 4.8833.
$$
Choosing $a$ strictly between $a_0$ and $1/[2(48833/10000 - 1)]$ absorbs all uniform $o(n)$ errors, giving $\operatorname{cw}(H) \le a|V| + O(1)$ for every connected simple cubic graph $H$. The norm-one reduction then proves the circuit claim for every sufficiently large integer length $N$.

## Layout from split-penalized softened anchors

Let $H = (V, E)$ be a finite simple cubic graph on $n$ vertices and fix $0 < \delta \le \lambda \le 1$. Given anchor bits $b \in \{0,1\}^V$ and edge bits $x \in \{0,1\}^E$, let $k_v(b_v, x_{\partial v}) = \sum_{e \ni v} |x_e - b_v| \in \{0,1,2,3\}$ and define
$$
q_{\lambda,\delta}(b_v, x_{\partial v}) = c_{k_v}, \qquad (c_0, c_1, c_2, c_3) = (0,\; 1,\; 1+\delta,\; \lambda),
$$
$Q_{\lambda,\delta}(b, x) = \sum_{v \in V} q_{\lambda,\delta}(b_v, x_{\partial v})$, and $G_{\lambda,\delta}(b) = \min_{x \in \{0,1\}^E} Q_{\lambda,\delta}(b, x)$.

**Theorem 1 (Unit-factor layout interpolation).** *Let $0 < \delta \le \lambda \le 1$, $r = (1-\lambda)/(3-\lambda)$, $\alpha_d = \frac{3}{4} - \frac{\lambda(1+r)}{4}$, and $\alpha_s = \frac{3}{4} + \frac{\lambda(1+r)}{4}$. Suppose*
$$
1 - \frac{r}{4} - \frac{\delta}{1+\delta} \ge 0,\quad
\lambda - \frac{r(2-\lambda)}{2} - \frac{(2-r)(1-\lambda)}{4} \ge 0,\quad
\lambda - \frac{\alpha_s}{2} - \frac{\lambda\delta(1 - r/4)}{1+\delta} \ge 0,\quad
\frac{7\lambda - 3}{4} \ge 0,
$$
*and that $A = \alpha_d + \frac{\lambda\delta r}{2(1+\delta)}$, $B = \alpha_d - \lambda + \frac{\lambda\delta(1 + r/2)}{1+\delta}$, $C = \frac{\lambda\delta}{1+\delta}$ satisfy $B \le 0$ or $B^2 \le 4AC$. Then for every finite simple cubic graph $H$ on $n$ vertices and every nested single-vertex anchor ordering $0 = b^{(0)} < b^{(1)} < \dots < b^{(n)} = 1$ with $\max_t G_{\lambda,\delta}(b^{(t)}) \le M$,*
$$
\operatorname{cw}(H) \le M + o(n)
$$
*uniformly over all $H$ and anchor orderings.*

At $(\lambda, \delta) = (27/50, 41/20000)$, exact rational arithmetic in `certify.py` verifies $r = 23/123$, all four inequalities strictly, and $B^2 - 4AC < 0$.

**Joint submodularity and monotonicity.** Each $q_{\lambda,\delta}(b_v, x_1, x_2, x_3)$ is symmetric in $(x_1, x_2, x_3)$ and complement-invariant: $q_{\lambda,\delta}(1-b_v, 1-x_{\partial v}) = q_{\lambda,\delta}(b_v, x_{\partial v})$. All two-coordinate cross-differences on $\{0,1\}^4$ are nonpositive:
* For $(x_1, x_2)$ at $b_v = 0$ and $x_3 = m \in \{0,1\}$, $c_{m+2} - 2c_{m+1} + c_m$ is $\delta - 1 \le 0$ ($m = 0$) or $\lambda - 1 - 2\delta \le 0$ ($m = 1$); complement symmetry handles $b_v = 1$.
* For $(b_v, x_1)$ at $m = x_2 + x_3 \in \{0,1,2\}$, $c_{2-m} - c_{3-m} - c_{m+1} + c_m$ is $\delta - \lambda \le 0$ ($m \in \{0,2\}$) or $-2\delta \le 0$ ($m = 1$).

Hence $Q_{\lambda,\delta}(b, x)$ is jointly submodular in $(b, x)$. Partial minimization preserves submodularity (applying the lattice inequality to optimizers at $b, b'$) and complement invariance. Flipping one anchor $b_v$ sends $k_v \mapsto 3 - k_v$, changing $q_{\lambda,\delta}$ by at most $\max(|c_0 - c_3|, |c_1 - c_2|) = \lambda$, so $G_{\lambda,\delta}$ is $\lambda$-Lipschitz per bit. Minimizers of $Q_{\lambda,\delta}(b, \cdot)$ are closed under meet, giving a unique least minimizer $x^*(b)$ that is nondecreasing in $b$, with $x^*(0) = 0$ and $x^*(1) = 1$ since $\lambda > 0$.

**Endpoint states and absence of $B_g$–$B_g$ edges.** Fix one transition $b^{(t)} \to b^{(t+1)}$ with least minimizers $x^{(0)} \le x^{(1)}$, and let $D$ be the graph of changed edges ($0 \to 1$) without isolated vertices. Fix integers $L \ge 2$ and $J = L^2$. Give each component of $D$ on at most $L$ vertices one label (zero fluctuation); let $D_L \subseteq D$ be the union of components of size $> L$.

Vertices of $D$ are type $P, A, B$ if split at both endpoints, only at $x^{(1)}$, or only at $x^{(0)}$, respectively; remaining vertices are trivalent in $D$, flipping $(0,0,0) \to (1,1,1)$. Every $P$ vertex in $D$ has degree $1$ in $D$ (degree $2$ leaves one fixed edge, forcing unanimity at one endpoint). Omitting the single changing anchor (an $O(1)$ perturbation to all counts and costs), partition vertices of $D_L$ by unchanging anchor:
* Trivalents: $z_d$ (anchor $0$, costs $(0, \lambda)$) and $z_s$ (anchor $1$, costs $(\lambda, 0)$); $z = z_d + z_s + O(1)$.
* Type $B$: good $B_g$ (anchor $1$) and bad $B_b$ (anchor $0$), with degree subscript $1$ or $2$ in $D_L$. Endpoint costs: $B_{g,2}: (1+\delta, 0)$, $B_{g,1}: (1, 0)$, $B_{b,2}: (1, \lambda)$, $B_{b,1}: (1+\delta, \lambda)$.
* Type $A$: good $A_g$ (anchor $0$) and bad $A_b$ (anchor $1$). Endpoint costs: $A_{g,2}: (0, 1+\delta)$, $A_{g,1}: (0, 1)$, $A_{b,2}: (\lambda, 1)$, $A_{b,1}: (\lambda, 1+\delta)$.
* Type $P$ (degree $1$ in $D_L$): $P_{1,1}$ (anchor $1$, costs $(1+\delta, 1)$) and $P_{1,0}$ (anchor $0$, costs $(1, 1+\delta)$).

Set $D_0 = B_{g,2} + B_{b,1} + P_{1,1}$ and $D_1 = A_{g,2} + A_{b,1} + P_{1,0}$ for the $k=2$ vertices in $D_L$ at source and destination. Write $P^{\mathrm{all}}, A^{\mathrm{all}}, B^{\mathrm{all}}$ for full-graph counts and $P, A, B, \ell = A_1 + B_{g,1} + B_{b,1} + P_1 + O(1)$ for counts in $D_L$. Since every vertex outside $D_L$ pays $\ge 1$ when split and $\ge 0$ otherwise, $(1-p)Q_{\lambda,\delta}(b^{(t)}, x^{(0)}) + p Q_{\lambda,\delta}(b^{(t+1)}, x^{(1)}) \le M$ gives for $p \in [0,1]$:
$$
M \ge P^{\mathrm{all}} + p A^{\mathrm{all}} + (1-p) B^{\mathrm{all}} + \mathrm{Budget}(p),
\tag{1}
$$
$$
\mathrm{Budget}(p) := (1-p)\delta D_0 + p\delta D_1 + \lambda\big[p z_d + (1-p) z_s + p B_b + (1-p) A_b\big].
$$

No edge $e \in E(D_L)$ joins two vertices $u, v \in B_g$: at $x^{(0)}$, each $w \in B_g$ has anchor $1$ and $k_w \in \{1,2\}$ incident $0$-edges, so flipping $e$ from $0$ to $1$ reduces $k_w$ by $1$ at both ends, dropping cost by $1$ ($1 \to 0$) or $\delta$ ($1+\delta \to 1$) per end, for a strict total drop $\ge 2\delta > 0$, contradicting optimality of $x^{(0)}$. No secondary tie-breaking functional is needed. Symmetrically, no edge of $D_L$ joins two vertices of $A_g$.

Because $D_L$ is a union of components of $D$, switching all edges of $D_L$ at either endpoint cannot improve $Q_{\lambda,\delta}$; accounting for the single changing anchor gives $Q_0(D_L) = Q_1(D_L) + O(1)$:
$$
B - A + \delta(D_0 - D_1) = \lambda(z_d - z_s + B_b - A_b) + O(1).
\tag{2}
$$
Expanding $B = B_{g,2} + B_{g,1} + B_{b,2} + B_{b,1}$ in (2) gives
$$
(1+\delta) D_0 = A + \delta D_1 + \lambda(z_d - z_s) + (\lambda B_b - B_{b,2}) - \lambda A_b + (P_{1,1} - B_{g,1}) + O(1).
$$
Using $A, \delta D_1, P_{1,1} \ge 0$, $B_{b,2} \le B_b$, and $B_{g,1} \le \ell$ yields
$$
(1+\delta) D_0 \ge \lambda(z_d - z_s) - (1-\lambda) B_b - \lambda A_b - \ell - O(1).
\tag{3}
$$

Let $\delta_{\mathrm{triv}}$ count changing edges between trivalents and nontrivalents in $D_L$. Since $B_g$ has no internal edges in $D_L$, counting half-edges out of $B_g$ gives $2 B_g - B_{g,1} \le \delta_{\mathrm{triv}} + (2A - A_1) + (2 B_b - B_{b,1}) + P_1 + O(1)$. Adding $\ell = A_1 + B_{g,1} + B_{b,1} + P_1 + O(1)$ and using $B_g = B - B_b + O(1)$ gives $\delta_{\mathrm{triv}} + \ell \ge 2(B - A - 2B_b) - O(1)$. Substituting (2) yields
$$
\delta_{\mathrm{triv}} + \ell \ge 2\lambda(z_d - z_s) - 2\big[(2-\lambda) B_b + \lambda A_b\big] - 2\delta(D_0 - D_1) - O(1).
\tag{4}
$$

**Weighted cycle-rank fragmentation.** Assign each edge of $D$ a priority label given an independent $\mathrm{Unif}[0,1]$ priority; edges flip at their label's priority. Degree-one vertices contribute $P, pA, (1-p)B$. Degree-two vertices with two labels split with probability $2p - p^2$ ($A$) or $1 - p^2$ ($B$), adding $p(1-p)$. Trivalents with $k \in \{2,3\}$ labels split with probability $1 - p^k - (1-p)^k = k\,p(1-p)$. Thus expected splits at time $p$ are at most $P^{\mathrm{all}} + p A^{\mathrm{all}} + (1-p) B^{\mathrm{all}} + p(1-p)S$, where $S$ charges multi-label degree-two vertices $1$ and multi-label trivalents $k$.

For each component $C \subseteq D_L$, $\beta(C) = |E(C)| - |V(C)| + 1 = (z_C - \ell_C)/2 + 1$. Let $\mathcal{K}(C)$ be the connected components $K$ induced by trivalents of $C$. Every cycle of $C$ outside $\bigcup_K K$ meets a degree-two vertex; splitting it into two pendant copies preserves connectivity and lowers $\beta(C)$ by $1$ without altering any $K$. Repeating until every edge outside $\bigcup_K K$ is a bridge splits $\beta(C) - \sum_K \beta(K)$ degree-two vertices (cost $\le 1$ in $S$ each). Leave each $K$ with $|K| \le L$ intact. In each large $K$ ($|K| > L$), repeatedly delete an internal degree-3 vertex whose removal keeps $K$ connected; each such vertex has no external incidence, deleted vertices are independent, and each deletion drops $\beta(K)$ by $2$. The connected remainder is a cactus (any non-cycle nonbridge block has a degree-3 internal vertex with no external edge, contradicting termination). Reinserting each deleted vertex as three pendant copies preserves connectivity and boundary attachments, costing $\le 3$ units of $S$ per deleted vertex, or $\le \frac{3}{2}\beta(K)$ per large $K$.

The expanded graph has a tree decomposition of bag size $\le \max(L, 3)$: one bag per small $K$, block-cut triple bags for each large cactus (triangulating cycles $v_0, \dots, v_{q-1}$ by $\{v_0, v_i, v_{i+1}\}$), and pair bags for remaining bridges. Repeatedly deleting weighted tree centroid bags from connected pieces of size $> J$ halves piece sizes; pieces in each dyadic range $(2^j J, 2^{j+1}J]$ are disjoint, so over the expanded graph ($\le 3n$ vertices) at most $O(nL/J)$ separator vertices are removed. Give all edges incident to each surviving piece ($\le J$ vertices, $\le 3J$ edges) one label and edges with both ends removed individual labels; re-identifying separator vertices costs $\le 3$ in $S$ each. With at most $n/L$ components $C \subseteq D_L$ and $n/L$ large $K$'s, let $m$ count trivalents in small $K$'s ($|K| \le L$), $d(K)$ count outgoing changed edges of $K$, and $\delta_L = \sum_{|K|>L} d(K)$. Using $\beta(K) = (|K| - d(K))/2 + 1$, dropping $-\sum_{|K|\le L}\beta(K) \le 0$, and substituting $\sum_{|K|>L}|K| = z - m$:
$$
S \le \frac{z - \ell}{2} + \frac{1}{2}\sum_{|K|>L}\frac{|K| - d(K)}{2} + O\!\left(\frac{n}{L} + \frac{nL}{J}\right)
  = \frac{3}{4} z - \frac{1}{4} H + O\!\left(\frac{n}{L} + \frac{nL}{J}\right),
  \qquad H := m + \delta_L + 2\ell.
\tag{5}
$$

**Destination reset and unit-factor proof.** In $x^{(1)}$, resetting all edges incident with a large $K$ to $0$ saves $\lambda[Z_d(K) - Z_s(K) - e(K)]$ inside $K$ ($e(K) \in \{0,1\}$ for the omitted anchor). Any exterior vertex touching $t$ outgoing edges changes $k$ by at most $t$ unit steps; since $0 < \delta \le \lambda \le 1$, each unit step $k \to k \pm 1$ increases $c_k$ by at most $\max(1, \delta, 1+\delta-\lambda) = 1$, so exterior cost rises by at most $d(K)$. Optimality of $x^{(1)}$ gives $d(K) \ge \lambda[Z_d(K) - Z_s(K) - e(K)]$. Summing over large $K$ (small $K$'s have signed trivalent imbalance $\le m$) yields
$$
H \ge \lambda(z_d - z_s) + (1-\lambda)m + 2\ell - O(1).
\tag{6}
$$
Since $\delta_L \ge \delta_{\mathrm{triv}} - 3m$, (4) gives
$$
H \ge 2\lambda(z_d - z_s) - 2\big[(2-\lambda)B_b + \lambda A_b\big] - 2\delta(D_0 - D_1) + \ell - 2m - O(1).
\tag{7}
$$
Combining (6) with weight $1-r$ and (7) with weight $r = (1-\lambda)/(3-\lambda)$ eliminates $m$ and inserts into (5) to give
$$
S \le \alpha_d z_d + \alpha_s z_s + \frac{r}{2}\big[(2-\lambda) B_b + \lambda A_b\big] + \frac{r\delta}{2}(D_0 - D_1) - \frac{2-r}{4}\ell + O\!\left(\frac{n}{L} + \frac{nL}{J} + 1\right).
\tag{8}
$$

For $p \in [0, 1/2]$, set $\mu(p) := \min(\bar{\mu}(p),\, p\frac{2-r}{4})$ where $\bar{\mu}(p) := \frac{\delta(1 - pr/2)}{1+\delta}$. Subtract $(1-p)\mu(p)\big[(1+\delta)D_0 - \lambda(z_d - z_s) + (1-\lambda)B_b + \lambda A_b + \ell\big] \ge 0$ (by (3), up to $O(1)$) from $\mathrm{Budget}(p) - p(1-p)S$ and inspect every variable's coefficient:
* $D_1$: $p\delta + p(1-p)\frac{r\delta}{2} \ge 0$.
* $D_0$: $(1-p)[\delta(1 - pr/2) - \mu(p)(1+\delta)] \ge 0$ since $\mu(p) \le \bar{\mu}(p)$.
* $\ell$: $(1-p)[p\frac{2-r}{4} - \mu(p)] \ge 0$ since $\mu(p) \le p\frac{2-r}{4}$.
* $A_b$: $(1-p)\lambda[1 - \frac{pr}{2} - \mu(p)] \ge (1-p)\lambda[1 - \frac{r}{4} - \frac{\delta}{1+\delta}] \ge 0$.
* $B_b$: using $\mu(p) \le p\frac{2-r}{4}$, at least $p\big[\lambda - (1-p)\big(\frac{r(2-\lambda)}{2} + \frac{(2-r)(1-\lambda)}{4}\big)\big] \ge 0$.
* $z_s$: using $\mu(p) \le \bar{\mu}(p)$, at least $(1-p)[\lambda - p\alpha_s - \lambda\bar{\mu}(p)]$; the bracket is linear in $p$ with values $\frac{\lambda}{1+\delta} > 0$ at $p = 0$ and $\lambda - \frac{\alpha_s}{2} - \frac{\lambda\delta(1-r/4)}{1+\delta} \ge 0$ at $p = 1/2$.
* $z_d$: the coefficient is $p\lambda - p(1-p)\alpha_d + (1-p)\lambda\mu(p)$.
  - If $\mu(p) = p\frac{2-r}{4}$, since $\alpha_d - \lambda\frac{2-r}{4} = \frac{3}{4}(1-\lambda)$, this is $p[\lambda - (1-p)\frac{3}{4}(1-\lambda)] \ge p\frac{7\lambda-3}{4} \ge 0$.
  - If $\mu(p) = \bar{\mu}(p)$, expanding gives $Ap^2 - Bp + C$, which is nonnegative for all $p \ge 0$ since $A, C > 0$ and either $B \le 0$ or $B^2 \le 4AC$.

Thus $\mathrm{Budget}(p) \ge p(1-p)S - O(n/L + nL/J + 1)$ for all $p \in [0, 1/2]$. Reversing $0 \leftrightarrow 1$ and $p \leftrightarrow 1-p$ covers $p \in [1/2, 1]$. Hence expected splits at every $p \in [0,1]$ are at most $M + O(n/L + nL/J + 1)$.

**Concentration and vertex ordering.** Take $L = \max(2, \lceil n^{1/16}\rceil)$, $J = L^2$, and $L$ equal time intervals. Each label covers $\le 3J$ edges, affecting splits by $\le 6J$, with influence sum $\le 2|E(D)| \le 3n$, giving variance $O(nJ)$ for each mesh split count and interval edge count. Since $LJ/\sqrt{n} = o(1)$, Chebyshev and a union bound give a realization within $n^{3/4}$ of expectation at all $2L+1$ counts; within each interval splits rise by at most twice its edges. Concatenating across all anchor transitions gives an edge ordering of $H$ with $\le M + o(n)$ split vertices at every prefix. Ordering vertices by their median incident edge index: for any proper nonempty prefix with largest included median $j$, crossing edges $< j$ charge distinct excluded vertices splitting just before $j$, crossing edges $> j$ charge distinct included vertices splitting just after $j$, and these two split sets differ only at edge $j$'s endpoints while edge $j$ adds $\le 1$ crossing, yielding $\operatorname{cw}(H) \le M + o(n)$ uniformly.

## From the infinite cubic tree to every finite cubic graph

**Exact cover descent.** For a finite degree-$d$ cover $\pi : \widehat{H} \to H$ with lifted anchors $\hat{b} = b \circ \pi$, lifting gives $G_{\lambda,\delta}(\hat{b}) \le d\,G_{\lambda,\delta}(b)$. Conversely, for any upstairs edge state $\hat{x}$, let $\bar{y}_e = d^{-1} \sum_{\hat{e} \in \pi^{-1}(e)} \hat{x}_{\hat{e}} \in [0,1]$. At $b_v = 0$, thresholding $y = (y_1, y_2, y_3) \in [0,1]^3$ (descending order $y_{(1)} \ge y_{(2)} \ge y_{(3)}$) at $U \sim \mathrm{Unif}[0,1]$ gives the Lovász extension
$$
\hat{q}_0(y) := \mathbb{E}_U q_{\lambda,\delta}\big(0, \mathbf{1}[y > U]\big)
= (1-\delta)\max_i y_i + \delta \max_{i < j}(y_i + y_j) + (1+\delta-\lambda)\max_i(-y_i),
$$
whose coefficients on $y_{(1)}, y_{(2)}, y_{(3)}$ are $1 = c_1$, $\delta = c_2 - c_1$, and $\lambda - 1 - \delta = c_3 - c_2$. Since $1-\delta, \delta, 1+\delta-\lambda \ge 0$, $\hat{q}_0$ (and by symmetry $\hat{q}_1(y) = \hat{q}_0(1-y)$) is convex on $[0,1]^3$. Applying Jensen on each fiber star gives
$$
\mathbb{E}_U Q_{\lambda,\delta}\big(b, \mathbf{1}[\bar{y} > U]\big)
= \sum_{v \in V(H)} \hat{q}_{b_v}(\bar{y}_{\partial v})
\le \frac{1}{d} \sum_{\hat{v} \in V(\widehat{H})} q_{\lambda,\delta}(\hat{b}_{\hat{v}}, \hat{x}_{\partial \hat{v}})
= \frac{1}{d} Q_{\lambda,\delta}(\hat{b}, \hat{x}),
$$
so some threshold $U$ proves $G_{\lambda,\delta}(b) \le d^{-1} G_{\lambda,\delta}(\hat{b})$, establishing equality.

**Gaussian comparison and threshold monotonicity.** For submodular complement-invariant $G : \{0,1\}^m \to \mathbb{R}$ and $X \sim N(0, C)$ with unit diagonals, set $f_C(t) = \mathbb{E} G(\mathbf{1}[X_i < t])$.
1. If $C^{(1)} \ge C^{(0)}$ entrywise, then $f_{C^{(1)}}(t) \le f_{C^{(0)}}(t)$: along $C_\eta(s) = (1-\eta)((1-s)C^{(0)} + s C^{(1)}) + \eta I$ ($\eta \downarrow 0$), Gaussian integration by parts on smoothed multilinear extensions yields $\frac{d}{ds} f_{C_\eta(s)}(t) = (1-\eta)\sum_{i<j} (C^{(1)}_{ij} - C^{(0)}_{ij}) \phi_{ij}(t,t)\,\mathbb{E}[\Delta_{ij} G \mid X_i = X_j = t] \le 0$ by submodularity.
2. For all $t \in \mathbb{R}$, $f_C(t) \le f_C(0)$: in the positive-definite case $f'_C(t) = \phi(t)\sum_i g_i(t)$ with $g_i(t) = \mathbb{E}[G_{i=1} - G_{i=0} \mid X_i = t]$. Writing $W_j = X_j - C_{ij} X_i$ (independent of $X_i$), the conditional tests $W_j < (1-C_{ij})t$ increase in $t$, so submodularity makes $g_i$ nonincreasing; complement symmetry $g_i(-t) = -g_i(t)$ forces $g_i(t) \le 0$ for $t \ge 0$ and $\ge 0$ for $t \le 0$.

**Radial waves, high-girth covers, and prefix concentration.** On the cubic tree $T$, fix $0 < \alpha < 1/\sqrt{2}$ and unit rows $h_v(z) = \sqrt{(1-2\alpha^2)/(1+\alpha^2)}\,\alpha^{d(v,z)}$. Projecting onto a length-$j$ path yields inner product $\alpha^j[(j+1) - (j-1)\alpha^2]/(1+\alpha^2) \to r_j = (1 + j/3)2^{-j/2}$ as $\alpha \uparrow 1/\sqrt{2}$. Truncating at radius $R$ and renormalizing converges uniformly on every finite ball as $R \to \infty$.

For any finite connected simple cubic $H$, let $D$ exceed all witness radii below and take a finite cover $\pi : \widehat{H} \to H$ of girth $> 4R + 2D + 10$. (Contracting a spanning tree identifies $\pi_1(H)$ with the free group on nontree edges; for each nonbacktracking closed walk of length $\le 4R+2D+10$, extend partial bijections along its reduced word on a finite line to permutations, and take a connected component of the regular cover from the finite product group.) On $\widehat{H}$, define truncated radius-$R$ rows by graph distance to get $\widehat{X}_0$, mixed with tie noise $\widehat{X} = \sqrt{1-\eta^2}\,\widehat{X}_0 + \eta\widehat{Z}$. On $H$, collapse squared entries of each tree row fiberwise under the universal cover, take per-column square roots to get $X_0$, and mix with base tie noise as $X = \sqrt{1-\eta^2}\,X_0 + \eta Z$. By Cauchy–Schwarz per fiber and nonnegative noise correlations on repeated base vertices, $\operatorname{Cov}(X \circ \pi) \ge \operatorname{Cov}(\widehat{X})$ entrywise. Cover descent, Gaussian comparison, and threshold monotonicity give
$$
|H|^{-1} \mathbb{E} G_{\lambda,\delta,H}(\mathbf{1}[X < t])
\le |\widehat{H}|^{-1} \mathbb{E} G_{\lambda,\delta,\widehat{H}}(\mathbf{1}[\widehat{X} < 0])
\qquad \forall t \in \mathbb{R}.
$$
Taking $\alpha \uparrow 1/\sqrt{2}$, $R \to \infty$, $\eta \downarrow 0$, every radius-$D$ ball upstairs converges uniformly to the limiting Gaussian wave (covariance $r_j$), whose continuous marginals assign probability $0$ to sign boundaries.

Set $\rho = r_1 = 2\sqrt{2}/3$. On the limiting wave, $\sum_{u \sim v} X_u = 2\sqrt{2} X_v$ a.s. (since $\operatorname{Var}(\sum_{u \sim v} X_u - 2\sqrt{2} X_v) = 3 + 6r_2 - 12\sqrt{2} r_1 + 8 = 0$). Given $X_v = x$, one named neighbor has law $N(\rho x, 1/9)$. At a vertex of value $b$ reached from parent $a$, its two forward children have values $\sqrt{2} b - a/2 \pm \zeta$ with $\zeta \sim N(0, (1-r_2)/2) = N(0, 1/12)$ orthogonal to all coordinates outside the two forward subtrees; contrasts at distinct vertices are mutually uncorrelated, hence independent. Conditionally on the root and its three neighbors, all forward innovations in the three subtrees are independent.

If the tree bound gives upstairs expected cost $\le s_0$ per vertex, then for any $\epsilon > 0$, fixed $(\alpha, R, \eta)$ before $H$ give $\mathbb{E} G_{\lambda,\delta,H}(\mathbf{1}[X < t]) \le (s_0 + \epsilon)n$ for all $t$. Each base Gaussian affects $O(2^R)$ vertices, so $\lambda$-Lipschitzness gives variance $O_R(n)$ for each threshold cost and quantile-interval count. Partitioning $\mathbb{R}$ into standard-normal quantile intervals of mass $\le \epsilon$, Chebyshev and a union bound give (for large $n$) a realization with distinct keys, $\le 2\epsilon n$ keys per interval, and cost $\le (s_0 + 2\epsilon)n$ at every finite boundary. Since flipping each anchor changes $G_{\lambda,\delta}$ by $\le \lambda \le 1$, every key-order prefix has cost $\le (s_0 + 4\epsilon)n$, and Theorem 1 yields $\operatorname{cw}(H) \le (s_0 + 5\epsilon)n + o(n)$ uniformly.

## Certified upper bounds for the reference crossing cover and crossing density

Color each vertex by $\operatorname{sgn}(X_v)$, calling an edge crossing if its endpoint signs differ.

**Exact crossing-edge density $W_1$.** For each edge $(u, v)$ with correlation $\rho = 2\sqrt{2}/3$, $\mathbb{P}(X_u X_v < 0) = \frac{1}{\pi}\arccos(\rho) = \frac{1}{\pi}\arcsin(\sqrt{1-\rho^2}) = \frac{1}{\pi}\arcsin(1/3)$. With $3/2$ edges per vertex, $W_1 = \frac{3}{2\pi}\arcsin(1/3)$. In `reference.py`, the positive series $\arcsin(1/3) = \sum_{k \ge 0} \frac{(2k)!}{4^k (k!)^2 (2k+1)} 3^{-(2k+1)}$ has consecutive ratio $< 1/9$, so 36 terms plus tail bound $\frac{9}{8}$ times the next term certify $W_1 < 4056504299/25000000000 = 0.16226017196$.

**Alternating crossing-path bound for $U$.** On the limiting zero-threshold wave, a crossing-degree-three vertex contradicts $\sum_{u \sim v} X_u = 2\sqrt{2} X_v$ a.s. Deleting crossing-degree-three vertices on a high-girth cover changes the minimum crossing cover $M$ by $\le 1$ and $j$-walk counts by $\le (j+1)3\cdot 2^{j-1}$ per deletion ($o(1)$ per vertex in the limit). In the remaining degree-$\le 2$ crossing graph, let $W_j$ be half the number of directed nonbacktracking crossing $j$-walks. An $\ell$-edge crossing path has $M = \lceil \ell/2 \rceil$ and contributes $\max(0, \ell-j+1)$ to $W_j$; an even $\ell$-edge crossing cycle has $M = \ell/2$ and contributes $\ell$ to every $W_j$. Thus for odd $J$, $\Gamma_J = \sum_{j=1}^J (-1)^{j+1} W_j \ge M$. On the tree, $W_j$ equals $3\cdot 2^{j-1}$ times the orthant probability of one alternating pattern on $m = j+1$ vertices.

**Gaussian path density, tails, and product quadrature.** From $r_{j+1} = \sqrt{2} r_j - r_{j-1}/2$ ($j \ge 1$) and $(3/2)r_1 = \sqrt{2}$, innovations $X_{i+1} - \sqrt{2} X_i + X_{i-1}/2 \sim N(0, 1/12)$ are independent along a path. Negating alternating coordinates puts the pattern in $\mathbb{R}_{>0}^m$ with start density $f_2(x_1, x_2) = c_2 \exp(-\frac{9}{2}(x_1^2 + x_2^2) - 6\sqrt{2} x_1 x_2)$ ($c_2 = 3/(2\pi)$) and step density $c_1 \exp(-6(z + \sqrt{2} y + x/2)^2)$ ($c_1 = \sqrt{6/\pi}$). The joint density is $f_m(0)\exp(-\frac{1}{2} x^\top Q x)$ with $f_m(0) = c_2 c_1^{m-2}$ and nonnegative bandwidth-two precision $Q$ having diagonals $(9,9)$ ($m=2$), $(12,33,12)$ ($m=3$), and $(12, 36, 39, \dots, 39, 36, 12)$ ($m \ge 4$).

A one-step positive continuation has probability $\le 1/2$; a two-step continuation is pointwise maximized at $(0,0)$, giving $\arctan(1/\sqrt{2})/(2\pi) < 1/10$ (since $\cos(2\pi/5) = (\sqrt{5}-1)/4 < 1/3$). By reversibility, an $a$-step arm has probability $\le 1/H_{\mathrm{arm}}(a)$ with $H_{\mathrm{arm}}(a) = 10^{\lfloor a/2\rfloor} 2^{a \bmod 2}$. At cutoff $L = 11/5$, anchoring coordinate $i$ at pair $(h, h+1)$ ($h = \min(i, m-2)$) bounds the $x_i > L$ tail by $H_{\mathrm{arm}}(h)^{-1} H_{\mathrm{arm}}(m-2-h)^{-1} \cdot \frac{3}{2\pi} \frac{e^{-9L^2/2}}{54\sqrt{2}\,L^2}$.

With a $q$-point Gauss–Legendre product rule ($q = 46$) on $[0, L]^m$, Gaussian Fourier moments bound $|\partial_i^{2q} f_m(x)| \le Q_{ii}^q \frac{(2q)!}{2^q q!} f_m(0)$ times the Schur-complement Gaussian, retaining $e^{-39x_j^2/2}$ on at least $\max(0, m-9)$ interior coordinates outside radius $2$ of $i$. Both $\int_0^L e^{-39x^2/2}dx < \sqrt{\pi/78} < 1/4$ and the 1D rule sum (checked in `reference.py`) are $< 1/4$. Combining with the 1D error constant $L^{2q+1}(q!)^4/[(2q+1)((2q)!)^3]$ yields the tensor quadrature bound
$$
\frac{f_m(0)\, L^{m+2q}\, (q!)^3 \sum_{i=1}^m Q_{ii}^q}{(2q+1)\, 2^q\, ((2q)!)^2\, (4L)^{\max(0, m-9)}}.
$$
In `reference.py`, 36-digit outward decimal intervals (using only explicit `Q_LOW`/`Q_HIGH` directed $+, -, \times, /$ and unrounded sign-bit operations `copy_abs`/`copy_negate`, with zero dependence on the ambient decimal context) enclose $\pi$ via Machin's formula, evaluate $\exp(x)$ via elementary range reduction $y = |x|/2^k \le 2^{-10}$ (forming $y$ by directed division of $x.\mathrm{copy\_abs}()$ by $2^k$) with a 12-term Taylor enclosure and remainder $\le 2y^{12}/12!$ (plus $k$ squarings and inversion for $x < 0$), certify $\sqrt{\cdot}$ endpoints directly by directed squaring, bracket all $46$ roots of $P_{46}$ to radius $10^{-29}$ by exact scaled integer Legendre recurrences, bound $|P''_{46}| \le 46^4$ by Markov's inequality, and contract pair states using factored kernels $A(x,y) = e^{-3x^2/2 - 6\sqrt{2}xy}$, $B(x,z) = e^{-6xz}$, $C(y,z) = w_z c_1 e^{-12y^2 - 6z^2 - 12\sqrt{2}yz}$ at scale $2^{104}$ with per-step debit $(D_A + D_B + q D_C + 2q^2)/2^{104}$. At $J = 33$, `reference.py` certifies $\Gamma_{33} < 0.13064756645 < U = 1633094581/12500000000$.

## Disjoint path switching and the pointwise cost bound

First suppose $H'$ is a finite cubic graph whose crossing subgraph $E_{\mathrm{cross}}(H')$ has maximum degree at most $2$. Let $K^+$ be the canonical independent minimum crossing cover: on even-edge paths take the uniquely smaller color class; on odd-edge paths and even cycles take the positive class; leave crossing-isolated vertices out. Any positive vertex with an odd-length terminal crossing arm lies in $K^+$, and any vertex of either sign with two odd-length terminal arms lies in the unique minimum cover.

**Pairwise-disjoint guard classes.** At a positive vertex $w$ entered from a same-positive parent, define four mutually disjoint forward-subtree events:
* $U$: both forward children of $w$ are positive (crossing degree $0$, so $w \notin K^+$).
* $W_1$: one forward child is negative and closes after an odd number ($1, 3, 5, \dots$) of further crossing edges beyond that child ($B_1, B_3, B_5, \dots$); the other child is positive. Then $w$ is an endpoint of an even-length crossing path ($2, 4, 6, \dots$ edges), so $w \notin K^+$.
* $W_2$: both forward children are negative and each closes after an odd number ($1, 3, \dots$) of further crossing edges ($B_1, B_3, \dots$). Then each arm from $w$ has even length ($2, 4, \dots$), placing $w$ at even distance from both ends of an even-length crossing path, so $w \notin K^+$.
* $V_{\mathrm{odd}}$: one forward child is negative and closes after an even number ($0, 2, 4, \dots$) of further crossing edges ($B_0, B_2, B_4, \dots$, giving an odd-length arm from $w$), and the other forward child is positive and belongs to $\mathrm{OUT} := U \sqcup W_1 \sqcup W_2$ relative to $w$.

Thus $U, W_1, W_2, V_{\mathrm{odd}}$ are pairwise disjoint, every vertex in $\mathrm{OUT}$ lies outside $K^+$, and every vertex in $V_{\mathrm{odd}}$ has a same-positive child in $\mathrm{OUT} \subseteq V \setminus K^+$ and therefore cannot save in $K^+$.

**Simultaneous odd-path switches.** List each odd-edge crossing path from its positive endpoint as $(u_1^+, v_1^-, \dots, u_m^+, v_m^-)$. Admit $u_1^+$ if its two off-path positive children $(w_-, w_+)$ (ordered by increasing magnitude; ties have probability $0$) satisfy
$$
(w_-, w_+) \in (\mathrm{OUT} \times \mathrm{OUT}) \;\sqcup\; (V_{\mathrm{odd}} \times \mathrm{OUT}) \;\sqcup\; (\mathrm{OUT} \times V_{\mathrm{odd}}).
$$
These three regions are disjoint since $\mathrm{OUT} \cap V_{\mathrm{odd}} = \emptyset$; every admitted pair contains an uncovered child in $\mathrm{OUT}$ (so $u_1^+$ cannot save in $K^+$) and neither child can save in $K^+$. Admit an internal positive vertex $u_j^+$ ($j \ge 2$) if its unique off-path positive neighbor lies in $\mathrm{OUT}$ relative to $u_j^+$. Let $r \in \{0, \dots, m\}$ be largest such that $u_1^+, \dots, u_r^+$ are all admitted, and switch $\{u_1^+, \dots, u_r^+\}$ to $\{v_1^-, \dots, v_r^-\}$ simultaneously on all odd paths to get $K$. Then $|K| = |K^+|$, $K$ covers all crossing edges of $H'$, and if $1 \le r < m$ the single boundary edge $e^* = (v_r^-, u_{r+1}^+)$ is double-covered.

**Selective labeling and pointwise cost bound.** Set edge bits $x_e$ on $H'$ as follows:
1. On each singly-covered crossing edge $(w, z)$ ($w \in K, z \notin K$), set $x_e = b_z$. On each double-covered edge $e^* = (v_r^-, u_{r+1}^+)$, set $x_{e^*} = -1 = b_{v_r^-}$ (disagreeing with $u_{r+1}^+$, agreeing with $v_r^-$).
2. Let $c_{\mathrm{dis}}(v) \in \{0,1,2\}$ count crossing edges at $v$ disagreeing with $b_v$. Define the saving set $S \subseteq K$ as all $v \in K$ whose incident crossing edges all disagree with $b_v$ (excluding $v_r^-$ at any double-covered edge) and whose same-sign neighbors in $H'$ all lie in $K$. Partition $S = S_1 \sqcup S_2$ by crossing degree $c(v) \in \{1,2\}$.
3. Flip a noncrossing edge $e = (u, v)$ ($x_e = -b_u$) if and only if $u, v \in K$ and $\{u, v\} \cap S \ne \emptyset$; otherwise set $x_e = b_u$.

Every crossing edge of $H'$ disagrees with exactly one endpoint in $K$ (including $e^*$, which disagrees only with $u_{r+1}^+$), and every $v \in K$ has $c_{\mathrm{dis}}(v) \in \{1,2\}$ (at $v_r^-$, the edge $(u_r^+, v_r^-)$ disagrees with $v_r^-$). Partitioning $K = K_1 \sqcup K_2$ by $c_{\mathrm{dis}}(v) \in \{1,2\}$ gives $|K_1| + |K_2| = |K|$ and $|K_1| + 2|K_2| = |E_{\mathrm{cross}}(H')|$, hence $|K_2| = |E_{\mathrm{cross}}(H')| - |K|$.
* Every $v \notin K$ has $k_v = 0$ (cost $0$). Every $v \in S$ has $k_v = 3$ (cost $\lambda$).
* For $v \in K \setminus S$ with $f(v)$ incident flipped noncrossing edges: if $v \in K_2 \setminus S_2$, flipping $1$ noncrossing edge drops cost from $1+\delta$ to $\lambda$, so $v$ pays $\le 1+\delta$; if $v \in K_1 \setminus S_1$, $k_v \in \{1,2,3\}$ costs $1, 1+\delta, \lambda$, so $v$ pays $\le 1 + \delta\,\mathbf{1}[f(v) \ge 1]$.
* Every flipped noncrossing edge meeting $v \in K_1 \setminus S_1$ has its other endpoint in $S$; since vertices of $S_1$ have $2$ noncrossing edges and vertices of $S_2$ have $1$, $|\{v \in K_1 \setminus S_1 : f(v) \ge 1\}| \le 2|S_1| + |S_2|$.

Summing over $V(H')$ yields the exact pointwise bound
$$
\begin{aligned}
Q_{\lambda,\delta}(b, x)
&\le |K_1 \setminus S_1| + (1+\delta)|K_2 \setminus S_2| + \lambda(|S_1| + |S_2|) + \delta(2|S_1| + |S_2|) \\
&= (1-\delta)|K| + \delta|E_{\mathrm{cross}}(H')| - (1-\lambda-2\delta)|S_1| - (1-\lambda)|S_2|.
\end{aligned}
\tag{9}
$$
Every vertex saving in $K^+$ keeps its saving in $K$ (guards exclude it and its same-positive neighbors from being switched out), and every switched negative vertex strictly before a double-covered boundary (or at the negative end of a fully switched odd path) earns its saving whenever its same-negative neighbors lie in $K$.

On a general finite cubic graph $\widehat{H}$, let $Z_3$ be the set of crossing-degree-three vertices, set $H' = \widehat{H} \setminus Z_3$, build $K$ on $H'$ as above, restrict $S_1, S_2$ to vertices outside $N_{\widehat{H}}(Z_3)$ (so every vertex of $V_{\mathrm{int}} := V(\widehat{H}) \setminus (Z_3 \cup N_{\widehat{H}}(Z_3))$ has full degree $3$ in $H'$), define $x_e$ on $E(H')$ by steps 1–3, and set $x_e$ arbitrarily on edges incident to $Z_3$. Every $v \in V_{\mathrm{int}}$ obeys the exact star bound above, while each of the $\le 4|Z_3|$ vertices in $Z_3 \cup N_{\widehat{H}}(Z_3)$ pays at most $1+\delta$, yielding
$$
Q_{\lambda,\delta,\widehat{H}}(b, x) \le (1-\delta)|K| + \delta|E_{\mathrm{cross}}(H')| - (1-\lambda-2\delta)|S_1| - (1-\lambda)|S_2| + 4(1+\delta)|Z_3|.
$$

## Integral cavity operators and the four local events

Put $r = \sqrt{2}$, $g(h) = 2\sqrt{6/\pi}\,e^{-6h^2}$ ($h \ge 0$), and $H(t) = \int_0^{\min(2,t)} g(h)\,dh$ for $t > 0$ (zero for $t \le 0$). All branch functions take nonnegative magnitudes $(a, b) \in [0,2]^2$ of (parent, current).

**Alternating branches and cover events.** For opposite parent/current signs $(a, b)$, set $B_0(a, b) = H(rb + a/2)$ and $B_{j+1} = K B_j$ ($j = 0, \dots, 5$) where $(K f)(a, b) = \int_0^2 g(d + rb + a/2) f(b, d)\,dd$. Set $I = \{0, 2, 4, 6\}$, $O = \{1, 3, 5\}$. For same-sign $(p, u)$ with $m = ru - p/2$, set
$$
S(f)(p, u) = \int_{\max(0,-m)}^2 g(d+m) f(u, d)\,dd, \qquad
S_1(f)(p, u) = \int_{\max(0,-2m)}^2 g(d+m) f(u, d)\,dd,
$$
and $D(f, k)(p, u) = \int_0^{\min(2,-m)} g(h) f(u, -m+h) k(u, -m-h)\,dh$ for $m < 0$ (zero if $m \ge 0$). Put $A_0 = \{S(B_j) : j \in I\}$ and $A = A_0 \cup \{D(B_i, B_j) : i \in \{1,3\}, j \in \{0,2,4\}\}$ (pairwise disjoint, certifying membership in $K^+$ via an odd terminal arm).

**Disjoint guard and switched-arm families.** At a same-positive child $(p, u)$, $U$ has probability $H(m)$; set $\mathcal{W} = \{W_{1,1}, W_{1,3}, W_{1,5}, W_{2,11}\} := \{S_1(B_1), S_1(B_3), S_1(B_5), D(B_1, B_1)\}$ and $\mathrm{OUT}_{\mathrm{cert}} = \{U\} \cup \mathcal{W} \subseteq \mathrm{OUT}$. For $i \in \{0,2,4\}$ and $o \in \{U, W_{1,1}, W_{1,3}\}$, define
$$
V_{i,o}(p, u) = \int_{\max(0,-2m)}^2 g(d+m) B_i(u, d) o(u, 2m+d)\,dd,
$$
with $U(u, 2m+d) = H(7u/2 - rp + rd)$, and set $\mathcal{V} = \{V_{0,U}, V_{2,U}, V_{4,U}, V_{0,W_{1,1}}, V_{2,W_{1,1}}, V_{0,W_{1,3}}\} \subseteq V_{\mathrm{odd}}$.

At a positive vertex $b$ with negative parent $a$, put $M = rb + a/2$. For $(f, k) \in (\mathrm{OUT}_{\mathrm{cert}} \times \mathrm{OUT}_{\mathrm{cert}}) \sqcup (\mathcal{V} \times \mathrm{OUT}_{\mathrm{cert}}) \sqcup (\mathrm{OUT}_{\mathrm{cert}} \times \mathcal{V})$, define
$$
T_{f,k}(a, b) = \int_0^{\min(2, M)} g(h) f(b, M-h) k(b, M+h)\,dh,
$$
with $U(b, M \pm h) = H(\frac{3}{2}b + a/\sqrt{2} \pm rh)$. For an internal positive guard $o \in \{U, W_{1,1}, W_{1,3}\}$, define
$$
(K'_o f)(a, b) = \int_0^2 g(d + rb + a/2) o(b, 2rb + a + d) f(b, d)\,dd,
$$
with $U(b, 2rb + a + d) = H(\frac{7}{2}b + r(a+d))$. Let $\mathcal{T}$ be the full disjoint family of pairs $T_{f,k}$, and define up to four guarded two-step hops (arm lengths $1, 3, 5, 7, 9$):
$$
\mathcal{E} = \bigcup_{s=0}^{4} \Big(\{K'_o K : o \in \{U, W_{1,1}, W_{1,3}\}\}\Big)^s (\mathcal{T}), \qquad
\mathcal{O}_g = \{K'_U B_j : j \in I\} \cup \{K'_{W_{1,1}} B_j : j \in \{0,2,4\}\} \cup \{K'_{W_{1,3}} B_0\}.
$$
Let $\mathcal{E}_{\mathrm{pt}} \subset \mathcal{E}$ be the explicit 34-element subfamily in `cavity.py` (20 pairs in $\mathcal{T}$ plus 14 $U$-guarded one-, two-, and three-hop continuations) and set $P_1 = \{S_1(e) : e \in \mathcal{E}_{\mathrm{pt}}\}$. Every element of $P_1$ certifies a same-negative neighbor in $K$ via a completely switched odd path; every element of $\mathcal{O}_g$ guards the first positive on an even-length far arm, ensuring a crossing-degree-two negative root is not at a double-covered boundary.

**Four disjoint root events.** For one specified root sign and magnitude $x > 0$:
* **Kind 1 ($\mathrm{pos}_1, \mathrm{neg}_1 \subseteq S_1$).** Condition on the root $x > 0$ and unique opposite neighbor $b > 0$, with joint density $f_2(x, -b) = \frac{3}{2\pi}\exp(-\frac{9}{2}(x^2 + b^2) - 6\sqrt{2}xb)$. The two same-sign neighbors have sum $s = b + 2\sqrt{2}x$ and values $\frac{s}{2} \pm \eta$ ($\eta \sim N(0, 1/12)$); by symmetry $|\eta| = h \in [0, s/2]$ has density $g(h)$. Setting $h = \frac{t}{2}s$ ($t \in [0,1]$, Jacobian $\frac{s}{2}$) gives $y = \frac{1+t}{2}s$, $z = \frac{1-t}{2}s$. Multiplying by $3$ for neighbor placement and restricting to the graded-grid box $[0, 23/32]\times[0, 5/4]\times[0, 1]$ (grid cutoffs above $7/10$ and $6/5$) bounds $\mathbb{P}(\mathrm{pos}_1)$ and $\mathbb{P}(\mathrm{neg}_1)$ below by
$$
3 \iiint_{[0,23/32]\times[0,5/4]\times[0,1]} f_2(x, -b)\, \frac{s}{2}\, g\!\left(\frac{ts}{2}\right) F(x, b, y, z)\, dx\,db\,dt,
$$
with $F = A(x,y) A(x,z) B_I(x,b)$ ($\mathrm{pos}_1$) and $F = P_1(x,y) P_1(x,z) E(x,b)$ ($\mathrm{neg}_1$), where each family denotes the sum of its disjoint itinerary functions.
* **Kind 2 ($\mathrm{pos}_2, \mathrm{neg}_2 \subseteq S_2$).** Condition on $x > 0$ and its two opposite neighbors $y, z > 0$, with joint density $f_3(x, -y, -z) = \frac{\sqrt{108}}{(2\pi)^{3/2}}\exp(-\frac{x^2}{2} - \frac{9}{2}(y + \rho x)^2 - \frac{3}{2}(2z + y + 2\sqrt{2}x)^2)$ and same-sign neighbor $b = 2\sqrt{2}x + y + z$. Multiplying by $3$ and restricting to the graded-grid box $[0, 17/32]\times[0, 17/16]^2$ (grid cutoffs above $1/2$ and $1$) bounds $\mathbb{P}(\mathrm{pos}_2)$ and $\mathbb{P}(\mathrm{neg}_2)$ below by
$$
3 \iiint_{[0,17/32]\times[0,17/16]^2} f_3(x, -y, -z)\, F(x, b, y, z)\, dx\,dy\,dz,
$$
with $F = A_0(x,b)[B_I(x,y) + B_O(x,y) B_I(x,z)]$ ($\mathrm{pos}_2$) and $F = P_1(x,b)[B_I(x,y) B_I(x,z) + E(x,y) O_g(x,z) + O_g(x,y) E(x,z)]$ ($\mathrm{neg}_2$).

Branch factors use disjoint forward subtrees (conditionally independent given the root star), bracket terms are disjoint by arm parity, and the four events are disjoint by root sign and crossing degree, all within bounded radius. Passing to a high-girth cover, deleting crossing-degree-three vertices $Z_3$, discarding witness centers within the witness radius of $Z_3$ (automatically excluding $N_{\widehat{H}}(Z_3)$), and applying the general cubic form of (9) gives limiting expected cost per vertex at most $a_0 = (1-\delta)U + \delta W_1 - (1-\lambda-2\delta)s_1 - (1-\lambda)s_2$.

### Log-concavity, Hölder compression, and integer replay

Every individual elementary itinerary $f_\alpha : [0,2]^2 \to [0,\infty)$ is jointly log-concave: Gaussians are log-concave, affine restrictions preserve log-concavity, and by Prékopa's theorem marginalizing one variable preserves log-concavity (via the increasing quantile transport $L$ between normalized 1D marginals $f(x_0,t)/m_0$ and $f(x_1,t)/m_1$ with $L'(t) = \frac{f(x_0,t)m_1}{m_0 f(x_1, L(t))}$ and $(1-\theta) + \theta L' \ge (L')^\theta$). Thus 1D integrals are bounded below by length times endpoint geometric mean; 3D integrals of coordinatewise log-concave integrands are bounded below by volume times 8-corner geometric mean (in Kind 1, $s, y, z, h,$ and Jacobian $s/2$ are each affine separately in $x, b, t$; in Kind 2, all arguments are jointly affine); and the minimum of log-concave $f$ on $[y_{\mathrm{lo}}, y_{\mathrm{hi}}]$ is $\min(f(y_{\mathrm{lo}}), f(y_{\mathrm{hi}}))$.

**Hölder compression of 4-corner geometric means.** On a graded grid $0 = g_0 < g_1 < \dots < g_S = 2$, write $\mathrm{GM}_4(f_\alpha; i, j) = \big(f_\alpha(g_i,g_j) f_\alpha(g_i,g_{j+1}) f_\alpha(g_{i+1},g_j) f_\alpha(g_{i+1},g_{j+1})\big)^{1/4}$.
1. *Endpoint pairs ($T_{\mathrm{GM4}}$).* On cell $[g_i, g_{i+1}] \times [g_j, g_{j+1}]$ with corners $(a_s, b_s)$ ($s=0,1,2,3$) and $M_s = r b_s + a_s/2$, pick $0 < U_s \le \min(13/10, M_s)$ and substitute $h = u U_s$ ($u \in [0,1]$, Jacobian $U_s$). Subdividing $[0,1]$ into $\nu$ uniform steps and applying 1D log-concavity in $u$, the four-term Hölder inequality $\prod_{s=0}^3 (\sum_k c_{s,k})^{1/4} \ge \sum_k \prod_{s=0}^3 c_{s,k}^{1/4}$ summed over all disjoint pairs $(f, k) \in \mathcal{T}$ (factoring into sums over $\mathrm{OUT}_{\mathrm{cert}}$ and $\mathcal{V}$) yields a single table $T_{\mathrm{GM4}}(i, j) \le \sum_{(f,k) \in \mathcal{T}} \mathrm{GM}_4(T_{f,k}; i, j)$.
2. *Guarded hops ($K$ and $K'_o$).* For $(L f_\alpha)(a_s, b_s) = \int_0^2 \kappa(a_s, b_s; d) f_\alpha(b_s, d)\,dd$, 1D log-concavity per interval $[g_r, g_{r+1}]$ ($r = 0, \dots, S-1$) and four-term Hölder across the four corners (two at $b_s = g_j$, two at $b_s = g_{j+1}$, whose four-corner product of edge roots equals $\mathrm{GM}_4(f_\alpha; j, r)$) summed over $\alpha$ and guards $o \in \{U, W_{1,1}, W_{1,3}\}$ give
$$
\sum_{\alpha, o} \mathrm{GM}_4(K'_o f_\alpha; i, j)
\ge \sum_{r=0}^{S-1} \Delta g_r \left(\sum_o \mathrm{GM}_4(\sqrt{\kappa_o(\cdot; g_r)\kappa_o(\cdot; g_{r+1})}; i, j)\right) F_{\mathrm{GM4}}(j, r),
$$
and identically for $K$. Iterating four times and summing yields $E_{\mathrm{GM4}}(i, j) \le \sum_{e \in \mathcal{E}} \mathrm{GM}_4(e; i, j)$ covering all endpoint pairs and four-hop continuations.

**Integer arithmetic (`cavity.py`).** All computations use exact integers at `SCALE = 10^8` on the 55-point graded grid (`n = 64`, `UNIT = 10^8` for $1/n$, spacings exact multiples of `UNIT/8`: $\{0, \frac{1}{8n}, \frac{1}{4n}, \frac{1}{2n}, \frac{1}{n}, \frac{3}{2n}\}$ near zero, $1/n$ on $[2/n, 3/8]$, $2/n$ on $(3/8, 3/4]$, $4/n$ on $(3/4, 5/4]$, and $8/n$ on $(5/4, 2]$).
* Scalar bounds `ROOT_LO = 141421356 < sqrt(2)*1e8 < 141421357 = ROOT_HI`, `INV_R2_LO = 70710678 < 1e8/sqrt(2)`, `RHO_HI = 94280905 > (2*sqrt(2)/3)*1e8`, `CG_NUM = 276395319 < 2*sqrt(6/pi)*1e8`, `C2_NUM = 47746482 < 3/(2*pi)*1e8`, and `C3_NUM = 65984522 < sqrt(108)/(2*pi)^(3/2)*1e8` are verified programmatically in `certify.py` by exact rational inequalities against `reference.py`'s Machin enclosure of $\pi$.
* `exp_floor(i)` evaluates the cubic Taylor lower bound $1 - x + x^2/2 - x^3/6$ at $x = i/(256\cdot 10^6)$ rounded downward and squares $8$ times with floor division (zero for $i > 1.4\cdot 10^7$). `density(v)` rounds $|v|$ upward to multiples of $1/40000$ and $6v^2$ upward to $10^{-6}$. `half_mass(r)` sums lower bounds on $2048$-mesh cells via Jensen on the exact cell mean $\frac{2(3i^2-3i+1)}{m^2}$, truncating $r$ downward.
* `geo_interp(a, b, num, den)` bounds $a^{1-t}b^t$ ($t = \mathrm{num}/\mathrm{den} \in [0,1]$) below via up to two downward-rounded `isqrt` bisections (skipping bisections when $|a-b|/\min(a,b) \le 1/48$, or the second when $\le 1/8$) to $[\mathrm{lo}, \mathrm{hi}]$ followed by the weighted harmonic mean $\lfloor \frac{\mathrm{den}\cdot\mathrm{lo}\cdot\mathrm{hi}}{(\mathrm{den}-r)\mathrm{hi} + r\cdot\mathrm{lo}} \rfloor \le \mathrm{lo}^{1-r/\mathrm{den}}\mathrm{hi}^{r/\mathrm{den}}$ (weighted AM–GM on reciprocals, monotone in $\mathrm{lo}, \mathrm{hi}$). On $[y_{\mathrm{lo}}, y_{\mathrm{hi}}]$ within one cell $[g_j, g_{j+1}]$, `eval_loc` uses `y_lo` if $v_{j+1} \ge v_j$ and `y_hi` otherwise (monotonicity of $v_j^{1-t} v_{j+1}^t$); across a grid node it takes the minimum of both endpoints (log-concavity). Cell lookups use an O(1) table `idx_lut` asserted identical to `bisect_right`.
* In Kind 1 ($n_t = 16$) and Kind 2, 8-corner `isqrt` geometric means are computed per itinerary (or read from `E_gm4`, `BI_gm4`, `BO_gm4`, `Og_gm4`) before summing disjoint itineraries (conservatively skipping nonnegative cells with density weight $< 1000$), certifying `(772452651795, 484013183112) / 419430400000000` and `(10161195286116, 7101915687234) / 13421772800000000`.

## From cubic layouts to the norm-one language

Write $c = 48833/10000$ and $g = 1/[2(c-1)] = 5000/38833 > a_0$. Fix $a \in (a_0, g)$ and $C$ such that every connected simple cubic graph $H$ has $\operatorname{cw}(H) \le a|H| + C$. All logarithms below are base $2$.

**The language.** Reject lengths $N < 16$. For $N \ge 16$, set $t = 2^{\lfloor \log_4 N \rfloor}$ and let $m$ be the largest odd integer $\le (N-t)/(t+1)$. Read $t$ bits for the low coefficients of monic $R(W) \in \mathbf{F}_2[W]$ of degree $t$ and $m$ bits for monic $Q(Z) \in \mathbf{F}_2[Z]$ of degree $m$ in ascending degree. Reject $m < 1$ or reducible $R$ or $Q$ over $\mathbf{F}_2$. Set $q = 2^m$, $F = \mathbf{F}_2[Z]/(Q)$, and $E = F[W]/(R)$. An irreducible factor of $R$ over $F$ of degree $b$ gives a field of order $2^{mb}$ containing $\mathbf{F}_{2^t}$, so $t \mid mb$; $\gcd(t, m) = 1$ forces $b = t$, hence $E$ is a field of order $q^t$. Read the next $tm$ bits as $z \in E$ on basis $Z^i W^j$ ($i$ fastest), ignore the remaining $\ell = N - t - m - tm$ bits, and accept iff $z^e = 1$ for $e = 1 + q + \dots + q^{t-1}$.

By Frobenius, $X^{2^d} - X$ is the square-free product of all monic binary irreducibles of degree dividing $d$, giving a polynomial-time irreducibility test via modular squaring and gcd. Degree-$d$ irreducibles exist since proper-divisor roots total at most $\sum_{j=1}^{\lfloor d/2\rfloor} 2^j < 2^d$. Field operations and exponentiation by $e$ ($< tm$ bits) run in polynomial time, placing this total language in $\mathbf{P} \subseteq \mathbf{NP}$ with empty witness (bound $0$).

**Rectangles and translate bound.** On any fixed valid header and tail, $z \mapsto z^e$ maps $E^*$ into the roots of $X^{q-1} - 1$ with kernel of size $\le e$, hence exactly $(q^t-1)/(q-1) = e$ accepting payloads. For distinct $d_1, \dots, d_t \in E$, at most $t!$ elements $y \in E$ satisfy $(y + d_i)^e = 1$ for all $i$: writing $x_j = y^{q^{j-1}}$ and $a_{ij} = -d_i^{q^{j-1}}$, the equations $\prod_{j=1}^t (x_j - a_{ij}) = 1$ have distinct entries in every column and $x_1 = y$. Homogenizing over $\overline{\mathbf{F}}_2$ in $(\mathbf{P}^1)^t$ gives $\prod_j (u_j - a_{ij} v_j) = \prod_j v_j$. If any $v_j = 0$, fewer than $t$ finite coordinates can each annihilate at most one row's left product, leaving a nonzero row, so no boundary solutions exist; projecting to each $\mathbf{P}^1$ shows the zero set is finite. Under the Segre embedding $(\mathbf{P}^1)^t \hookrightarrow \mathbf{P}^{2^t-1}$, each equation is a hyperplane. The degree-$h$ Hilbert function is $(h+1)^t$, so the Segre variety has degree $t!$; successive hyperplane cuts never increase the sum of degrees of irreducible components by the hyperplane Bézout theorem, leaving at most $t!$ points.

Partial payload assignments embed injectively in $E$ and add to $z$, forbidding any payload one-rectangle with $t$ rows and $t!+1$ columns. Put $K = 2^{t+m+\ell}(t!+1)$. Any $K \times K$ one-rectangle on $\{0,1\}^N$ fixes a common valid header and tail on each side by pigeonholing, leaving $\ge t!+1$ payloads per side, impossible. For $S_N = L \cap \{0,1\}^N$, $t, m = \Theta(\sqrt{N})$ and $0 \le \ell < 2(t+1)$ give $\log K = O(\sqrt{N}\log N)$ and $|S_N| \ge e \ge 2^{(t-1)m}$, so $N - \log |S_N| = O(\sqrt{N})$.

**Full-$B_2$ circuit lower bound.** At most $2\lceil \log K \rceil - 1$ inputs can be ignored (else two disjoint sets of $\lceil \log K \rceil$ ignored bits yield a $K \times K$ one-rectangle). Take an exact full-$B_2$ circuit of size $s$ for $f_N$. For large $N$ the output is a gate; trim to output ancestors ($s' \le s$ gates), absorb constant pins/subgates into descendant truth tables (replacing a constant gate by a surviving feeder if needed), and drop newly unused inputs, leaving $b$ used inputs with $r = N - b \le 2\lceil \log K \rceil - 1 = o(N)$. Expand each fan-out of $h \ge 1$ pins into a binary tree of $h-1$ copy vertices (enforcing equality). The connected loopless multigraph $G$ has maximum degree $3$, $V = 2s' + 1$ vertices ($s' - b + 1$ copies), and cycle rank $\nu(G) = |E(G)| - V + 1 = s' - b + 1$.

A width-$w$ ordering of $G$ gives a layered nondeterministic read-once program whose states are edge-cut bit assignments and whose transitions satisfy each vertex's table (reading an input at primary inputs and requiring output $1$), with at most $4V 2^w$ transitions. By acyclicity, paths accept exactly the circuit's inputs. Fixing one accepting path per accepted assignment on the $b$ read bits (at least $|S_N|/2^r$ assignments), charge each to the first transition where accepting suffixes drop from $\ge K$ to $< K$ (existing eventually since $|S_N|/2^r \ge K$). Each state defines a one-rectangle on read bits (transferred by fixing the $r$ unread bits), so its source has $< K$ prefixes and at most $K^2$ assignments charge any transition, yielding
$$
w \ge \log |S_N| - r - 2\log K - \log(4V).
\tag{10}
$$

**Subcubic width-versus-rank bound.** Every connected loopless multigraph $G$ of maximum degree $3$ on $V$ vertices with cycle rank $\nu$ has $\operatorname{cw}(G) \le 2a\nu + O(\log^2(V+1))$.
1. *Bridgeless case.* Repeatedly contract any adjacent pair containing a degree-$\le 2$ vertex or parallel edges, deleting internal edges. Bridgelessness, connectivity, and max degree $3$ persist, and cycle rank never rises, ending at a singleton or simple cubic $H$ with $|H| = 2(\nu(H)-1) \le 2\nu$ and width $\le 2a\nu + C$. Expand each contracted node in-place via its binary contraction tree larger-child-first: at any cut inside a partly expanded node, completed siblings along the path to the first missing leaf each halve the remaining size, so at most $\lceil \log V \rceil$ completed siblings occur, each with $\le 3$ boundary edges, adding $\le 3\lceil \log V \rceil$.
2. *General case.* Delete bridges and pick a vertex-weighted centroid block $B_0$ so every hanging component has $\le V/2$ vertices. Order $B_0$ by step 1 and insert each hanging layout just before its attachment vertex in $B_0$: only one hanging layout is partial at any cut, and crossing attachment bridges meet the next attachment vertex ($\le 3$ edges). Cycle ranks sum over bridges, and depth $\le \log(V+1)$ adds $O(\log^2(V+1))$.

If $s < cN$, then $V = O(N)$, so (10) and the subcubic bound give $2a(s' - b + 1) + O(\log^2 N) \ge N - O(\sqrt{N}\log N)$, hence $s \ge s' \ge (1 + \frac{1}{2a})N - O(\sqrt{N}\log N)$, contradicting $s < cN$ for all sufficiently large integer lengths $N$ since $1 + \frac{1}{2a} > c$.
