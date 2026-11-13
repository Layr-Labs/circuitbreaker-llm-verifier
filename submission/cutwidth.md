# Cutwidth of cubic multigraphs is at most γD + 24D/log₂(D+2), with γ = (3/(2π))·arccos((1+2√2)/4) < 1/6

Status: one of the three proof files of this submission (with `circuit.md` and `krs.md`). Theorem 1 below is
the hypothesis H(γ,h) of `circuit.md` §0 for the explicit γ and h displayed there; `proof.md` §6 combines this
file with `circuit.md`.

**Credit.** The argument (Gaussian edge scores, median ordering, the constant γ, and the compression and
expansion steps that reduce to simple cubic graphs) is from S. McGuire, S. Schlesinger, R. Williams (MSW),
"Explicit functions requiring larger than 4.5n size circuits" (preprint, November 2026); the authors report a
Lean 4 formalization (`proof.md` §8). This text is an independent write-up and uses neither as a premise. It
adds an explicit remainder; the weight 2^{−d/2} with an explicit 1/R truncation loss; a tie rule that costs +1;
the extension of the reduction to cubic multigraphs with loops and parallel edges; and rational witnesses for
every numerical inequality.

**Labels.** [P] = proved here. The only outside facts are those of undergraduate analysis and probability
listed in §1.0. Notation is local to this file: N (§2) and n (§3, Theorem 2) are vertex counts, not the input
length N of `circuit.md`.

## 0. Statements

A **multigraph** K has a finite vertex set V and a finite edge set E. Each edge has an unordered pair of ends
{u, v}, and u = v is allowed (a *loop*). deg(v) is the number of non-loop edges at v plus twice the number of
loops at v. K is **cubic** if every degree is 3. For P ⊆ V, **cut(P)** is the number of edges with exactly one
end in P; loops are never counted. An **ordering** σ = (v₁, …, v_D) of V has cutwidth
cw(σ) := max_{1≤i≤D−1} cut({v₁,…,v_i}) (0 if D = 1), and cw(K) := min_σ cw(σ). This agrees with
`circuit.md` Definition 4.3 (there the empty and the full prefix are also allowed; their cuts are empty).

Throughout,

    γ := (3/(2π)) · arccos(x₀),   x₀ := (1 + 2√2)/4,   θ₀ := arccos x₀ ∈ (0, π/2).

**Theorem 1** [P]. Every connected cubic multigraph K on D vertices satisfies

    cw(K) ≤ γ·D + h(D),   h(D) := 24·D / log₂(D + 2).

Numerically 0.1403509686 < γ < 0.1403509687 [P, (W3)]. In particular γ < 14035097/10⁸ < 1/6, and
1/6 − γ > 0.026315.

h is nonnegative and nondecreasing on [0, ∞), and h(D)/D = 24/log₂(D+2) → 0 [P]. Nonnegativity is clear
(log₂(D+2) ≥ 1). Indeed d/dx [x/ln(x+2)] has the sign of g(x) = ln(x+2) − x/(x+2), and g(0) = ln 2 > 0,
g′(x) = x/(x+2)² ≥ 0. So **H(γ,h) holds** in the sense of `circuit.md` §0 (there h is used only at positive
integers, where it takes values in [0, ∞)). The bound needs no hypothesis on D.

**Theorem 2** [P] (simple cubic graphs, sharper remainder). Every simple cubic graph H on n ≥ 32 vertices satisfies

    cw(H) ≤ γ·n + 9.45·n/(log₂n − 4) + √(3·log₂n + 6) · n^{3/4} + 1.

**Plan.** §1 collects elementary facts. §2 reduces cubic multigraphs to simple cubic graphs at an additive cost of
3⌊log₂D⌋ + 3 (Corollary 2.3). §3 proves Theorem 2: random median ordering (§3.1–3.3), one-threshold bound
(§3.4), and all cut positions at once (§3.5). §4 proves Theorem 1, and §5 gives every numerical inequality with
rational witnesses.

## 1. Preliminaries

**1.0 Assumed.** Fundamental theorem of calculus. Differentiation under the integral sign when the integrand and
its t-derivative are continuous and the derivative has an integrable bound independent of t. Change of variables,
including polar coordinates in R². Probability on R^n with a density. Linearity of expectation. The union bound.
Independence of functions of disjoint sets of independent random variables. Chebyshev's inequality,
P[|X − EX| ≥ λ] ≤ Var X/λ², which is Markov's inequality applied to (X − EX)². For indicators X, Y,
|Cov(X,Y)| ≤ (Var X · Var Y)^{1/2} ≤ 1/4, by Cauchy–Schwarz.

**1.1 Two cut facts** [P].
(F1) cw(σ) is at most the number of non-loop edges, so cw(K) ≤ 3D/2 for cubic K.
(F2) If K′ is K with all loops deleted, cut_K(P) = cut_{K′}(P) for every P. Hence cw(K) = cw(K′).

**1.2 Gaussian facts.** Let V be a finite set and ω = (ω_z)_{z∈V} have independent standard normal coordinates.
Write φ(s) = e^{−s²/2}/√(2π) and Φ(s) = ∫_{−∞}^s φ.

**(G0)** [P]. If u₁, …, u_k ∈ R^V are orthonormal, then ⟨u₁,ω⟩, …, ⟨u_k,ω⟩ are independent standard normals.
*Proof.* Extend to an orthonormal basis and let U be the orthogonal matrix with these rows. Then g = Uω has density
(2π)^{−|V|/2} exp(−|Uᵀg|²/2)·|det U| = (2π)^{−|V|/2} exp(−|g|²/2). Now take marginals. ∎

**Lemma 1.1 (Sheppard)** [P]. Let y, y′ be unit vectors with r = ⟨y,y′⟩ and X = ⟨y,ω⟩, Y = ⟨y′,ω⟩. Then
P[XY < 0] = arccos(r)/π.
*Proof.* If y′ = y, then X = Y and the probability is 0 = arccos(1)/π. If y′ = −y, then XY = −X² < 0 a.s., and
the probability is 1 = arccos(−1)/π. Otherwise choose an orthonormal basis u₁, u₂ of span(y, y′). By (G0),
g = (⟨u₁,ω⟩, ⟨u₂,ω⟩) is standard normal in R², with X = ⟨ŷ, g⟩ and Y = ⟨ŷ′, g⟩ for unit vectors
ŷ = (cos α, sin α) and ŷ′ = (cos β, sin β) in R². We may take |α − β| = ψ := arccos r ∈ (0, π). In polar
coordinates g = ϱ(cos ϑ, sin ϑ) the density is (2π)^{−1} e^{−ϱ²/2} ϱ dϱ dϑ, so ϑ is uniform on [0, 2π). Now
X > 0 iff ϑ lies in the open half-circle centred at α, and Y > 0 iff ϑ lies in the one centred at β. Their symmetric
difference is two arcs of length ψ each. XY < 0 iff ϑ is in one of them and ϑ is not on the boundary, a null set.
Hence the probability is 2ψ/(2π). ∎

**Lemma 1.2 (zero is the worst threshold)** [P]. Let y, y′, r, X, Y be as in Lemma 1.1 with r < 1. Then for
every t ∈ R,

    P_t(r) := P[min(X,Y) < t < max(X,Y)] ≤ arccos(r)/π.

*Proof.* If r = −1, the bound is 1 and holds trivially. Let |r| < 1 and s := √(1 − r²). With u₁ = y and
u₂ = (y′ − r y)/s, (G0) gives X = g₁ and Y = r g₁ + s g₂ with g₁, g₂ independent standard normals. The pair
(X, Y) has density f(u, v) = (2πs)^{−1} exp(−(u² − 2ruv + v²)/(2s²)), which is symmetric in u and v. Each of X, Y
equals any fixed value with probability 0. So P_t = P[X < t] + P[Y < t] − 2P[X < t, Y < t] = 2Φ(t) − 2G(t), with

    G(t) = P[X < t, Y < t] = ∫_{−∞}^{t} φ(u) Φ((t − r u)/s) du.

The integrand and its t-derivative φ(u)φ((t−ru)/s)/s are continuous, and the latter is at most φ(u)/(s√(2π)).
Hence G′(t) = φ(t)Φ((t − rt)/s) + ∫_{−∞}^t f(u, t) du. By symmetry, ∫_{−∞}^t f(u,t) du = ∫_{−∞}^t f(t,u) du,
and since f(t, u) = φ(t)φ((u − rt)/s)/s this equals φ(t)Φ((t − rt)/s). So G′(t) = 2φ(t)Φ(βt), where
β := (1 − r)/s ≥ 0, and

    P_t′ = 2φ(t) − 2G′(t) = 2φ(t)·(1 − 2Φ(βt)).

This is ≥ 0 for t ≤ 0 and ≤ 0 for t ≥ 0, so P_t ≤ P₀. Finally P₀ = P[XY < 0] = arccos(r)/π by Lemma 1.1. ∎

**Lemma 1.3** [P]. Let y be a unit vector and Z = ⟨y, ω⟩, a standard normal by (G0).
(a) P[Z ∈ [a, b)] ≤ (b − a)/√(2π), and P[Z = t] = 0 for every t.
(b) If y′ ≠ y is a unit vector, then P[⟨y,ω⟩ = ⟨y′,ω⟩] = 0.
(c) For T > 0, P[Z < −T] = P[Z > T] ≤ φ(T)/T.
*Proof.* (a): φ ≤ 1/√(2π). (b): ⟨y − y′, ω⟩ = |y − y′|·(standard normal) by (G0), and |y − y′| > 0; apply (a).
(c): ∫_T^∞ φ(u) du ≤ ∫_T^∞ (u/T) φ(u) du = φ(T)/T. ∎

**Lemma 1.4 (arccos on [0,1])** [P]. Let f = arccos.
(a) f is decreasing on [−1, 1].
(b) For m ∈ [0, 1) and x ∈ [0, 1], f(x) ≤ f(m) + (m − x)/√(1 − m²).
(c) For c₁, c₂, c₃ ∈ [0, 1], f(c₁) + f(c₂) + f(c₃) ≤ 3 f((c₁ + c₂ + c₃)/3).
*Proof.* (b): On [0,1), f′(s) = −(1 − s²)^{−1/2} is nonincreasing. For m ≤ x < 1, f(x) − f(m) = ∫_m^x f′ ≤
f′(m)(x − m). For 0 ≤ x ≤ m, f(m) − f(x) = ∫_x^m f′ ≥ f′(m)(m − x). Both rearrange to (b). The case x = 1
follows by continuity. (c): Let m be the mean. If m = 1, all c_i = 1. Otherwise apply (b) to each c_i and add;
the terms (m − c_i) sum to 0. ∎

## 2. From cubic multigraphs to simple cubic graphs

In this section G is a **connected loopless** multigraph with maximum degree at most 3 and N ≥ 1 vertices. For
S ⊆ V(G), ∂S is the set of edges with exactly one end in S, so |∂S| = cut(S). A *partition with trees* is a
partition V(G) = B₁ ⊔ … ⊔ B_k into *blocks* together with, for each block, a rooted tree T_i whose leaves are the
vertices of B_i and in which every internal node has exactly two children. For a node ν, L(ν) is the set of leaves
below ν. The **quotient** Q has vertex set {1, …, k} and one edge {i, j} for each edge of G with one end in B_i and
the other in B_j, i ≠ j. Edges inside a block are dropped, so Q is loopless, and deg_Q(i) = |∂B_i|.

**Lemma 2.1 (compression)** [P; MSW's compression step, rewritten]. There is a
partition with trees such that
(a) |∂L(ν)| ≤ 3 for every node ν of every tree;
(b) if k ≥ 2, then Q is a simple cubic graph.

*Proof.* Start with singleton blocks and one-node trees. Then (a) holds because |∂{v}| = deg_G(v) ≤ 3, as G is
loopless, and Q = G is connected. Repeat the following **merge** while k ≥ 2 and it applies:

> pick blocks B ≠ B′ joined by at least one Q-edge such that deg_Q(B) ≤ 2, or deg_Q(B′) ≤ 2, or B and B′ are
> joined by at least two Q-edges; replace them by B ∪ B′, whose tree has a new root with children the roots of
> T_B and T_{B′}.

Let m ≥ 1 be the number of G-edges between B and B′. Then |∂(B ∪ B′)| = |∂B| + |∂B′| − 2m. This is
≤ 2 + 3 − 2 = 3 if one of the degrees is ≤ 2, and ≤ 3 + 3 − 4 = 2 if m ≥ 2. So (a) holds at the new root, and old
nodes are unchanged. The new quotient is the old one with the Q-edges between B and B′ contracted, so it stays
connected. Each merge lowers k, so the process stops.

Suppose it stops with k ≥ 2. Q is connected with at least 2 vertices, so every block has a Q-neighbour. If some
block had Q-degree ≤ 2, a merge would apply. So every Q-degree is ≥ 3, and by (a) at the roots, deg_Q(i) =
|∂B_i| ≤ 3, so every Q-degree is exactly 3. No two blocks are joined by two Q-edges, and Q is loopless. Hence Q is
simple and cubic. ∎

**Lemma 2.2 (expansion)** [P; MSW's expansion step, with the count made exact]. Take the partition of
Lemma 2.1 and any ordering π = (i₁, …, i_k) of {1, …, k}. Let w_Q be its cutwidth in Q, with w_Q := 0 if k = 1.
In each T_i fix the *big-first leaf order*: the order in which a depth-first traversal visits the leaves when, at
every internal node, it enters first the child with more leaves (ties broken arbitrarily). Let σ list B_{i₁}, then
B_{i₂}, …, each block in its big-first leaf order. Then

    cw(σ) ≤ w_Q + 3(⌊log₂ N⌋ + 1).

*Proof.* Let P be a proper nonempty prefix of σ. Write P = A ∪ S, where A = B_{i₁} ∪ … ∪ B_{i_{j−1}} and S is a
nonempty prefix of the leaf order of T := T_{i_j}.

*Case S = B_{i_j}.* Then ∂P is exactly the set of G-edges between the first j blocks and the others. These are the
Q-edges crossing the Q-prefix {i₁, …, i_j}, which is proper because P ≠ V. So |∂P| ≤ w_Q.

*Case S ≠ B_{i_j}.* Let λ be the last leaf of S, and let ν₀ (the root), ν₁, …, ν_q = λ be the path from the root
to λ. Let J be the set of p < q such that ν_{p+1} is the child of ν_p visited *second*, and let c_p be the other
child. In a depth-first traversal the leaves below a node are consecutive, and those below the first child come
before those below the second. A leaf μ ≠ λ precedes λ iff, at the deepest common ancestor ν_p of μ and λ, μ is
below the first child and λ below the second, that is, iff p ∈ J and μ ∈ L(c_p). Hence

    S = {λ} ∪ ⋃_{p∈J} L(c_p).

For p ∈ J the second child has at most half the leaves: |L(ν_{p+1})| ≤ |L(ν_p)|/2. So 1 = |L(λ)| ≤ N·2^{−|J|},
and |J| ≤ ⌊log₂N⌋.

Take e ∈ ∂P with end a ∈ P and end b ∉ P.
- If a ∈ S, then a lies in one of the at most ⌊log₂N⌋ + 1 sets Z ∈ {{λ}} ∪ {L(c_p) : p ∈ J}, and b ∉ P ⊇ Z, so
  e ∈ ∂Z. By Lemma 2.1(a) (λ and the c_p are nodes), at most 3(⌊log₂N⌋ + 1) edges arise this way.
- If a ∈ A, then j ≥ 2 and b lies in B_{i_j} ∪ … ∪ B_{i_k}. So e is a Q-edge crossing the proper Q-prefix
  {i₁, …, i_{j−1}}, and there are at most w_Q such edges. If j = 1 then A = ∅ and this case is empty.

Adding the two counts gives the bound. ∎

**Corollary 2.3** [P]. Let K be a connected cubic multigraph on D vertices. Then either cw(K) ≤ 3⌊log₂D⌋ + 3, or
there is a simple cubic graph H on n vertices, 4 ≤ n ≤ D, with

    cw(K) ≤ cw(H) + 3⌊log₂D⌋ + 3.

*Proof.* By (F2), cw(K) = cw(G) for G = K minus its loops. G is connected (a loop never connects anything),
loopless, has maximum degree ≤ 3, and has N = D vertices. Apply Lemma 2.1. If k = 1, Lemma 2.2 gives
cw(G) ≤ 3⌊log₂D⌋ + 3. If k ≥ 2, let H := Q, simple cubic on k ≤ D vertices (so k ≥ 4), and let π be an optimal
ordering of H; Lemma 2.2 gives the bound. ∎

## 3. Simple cubic graphs: the random median ordering

Fix a simple cubic graph H with vertex set V, |V| = n ≥ 4, graph distance d, and an integer R ≥ 1. H need not
be connected: d(x,z) := +∞ if x and z lie in different components (all weights below are then 0). Write
S_j(x) := {z : d(x,z) = j} and B_j(x) := {z : d(x,z) ≤ j}. Set

    r* := 2√2/3,   ρ_R := r*·(1 − 1/R) ∈ [0, r*),   κ(ρ) := (3/(2π))·arccos((1 + 3ρ)/4)   (0 ≤ ρ ≤ 1).

Note κ(r*) = γ, since (1 + 3r*)/4 = (1 + 2√2)/4 = x₀. Since ρ_R < r*, Lemma 1.4(a) gives κ(ρ_R) ≥ γ.

### 3.1 Field rows

For x ∈ V let a_x ∈ R^V be the row a_x(z) := 2^{−d(x,z)/2} if d(x,z) ≤ R, and 0 otherwise. Let â_x := a_x/|a_x|;
note |a_x| ≥ a_x(x) = 1.

**Lemma 3.1 (spheres)** [P]. |S₁(x)| ≤ 3, and |S_{j+1}(x)| ≤ 2|S_j(x)| for j ≥ 1. Hence |S_j(x)| ≤ 3·2^{j−1}
(j ≥ 1), |B_j(x)| ≤ 3·2^j − 2, and the sequence |S_j(x)|·2^{−j} is nonincreasing for j ≥ 1.
*Proof.* Every w ∈ S_{j+1}(x) has a neighbour in S_j(x), the previous vertex on a shortest path. Every u ∈ S_j(x)
with j ≥ 1 has a neighbour in S_{j−1}(x), so at most 2 of its 3 neighbours lie in S_{j+1}(x). Counting adjacent
pairs (u, w) ∈ S_j × S_{j+1} gives |S_{j+1}| ≤ 2|S_j|. Then |B_j| ≤ 1 + 3(1 + 2 + … + 2^{j−1}) = 3·2^j − 2. ∎

**Lemma 3.2 (adjacent rows)** [P; a bound of the same 1/R form is Lemma 3.4 of S. Gaspers, arXiv:2611.05146v1
(preprint, November 2026), not used here]. For every edge xy, ⟨â_x, â_y⟩ ≥ ρ_R. If moreover B_R(x) = B_R(y) = V
(for instance R ≥ diameter), then ⟨â_x, â_y⟩ ≥ r*.

*Proof.* Since x and y are adjacent, they lie in the same component, and |d(x,z) − d(y,z)| ≤ 1 for every z in
that component; for z outside it, d(x,z) = d(y,z) = +∞ and a_x(z) = a_y(z) = 0, so such z contribute to no sum
below. We claim that for z ∈ B_R(x) ∩ B_R(y),

    a_x(z) a_y(z) ≥ (√2/3)·(a_x(z)² + a_y(z)²).                                                     (3.1)

If d(x,z) = d(y,z) = j, the left side is 2^{−j} and the right side is (2√2/3)·2^{−j} < 2^{−j}. If
{d(x,z), d(y,z)} = {j, j+1}, the left side is 2^{−j−1/2} and the right side is (√2/3)·(3/2)·2^{−j} = 2^{−j−1/2},
with equality. Sum (3.1) over z ∈ B_R(x) ∩ B_R(y), the support of a_x·a_y:

    ⟨a_x, a_y⟩ ≥ (√2/3)·( |a_x|² − ℓ_x + |a_y|² − ℓ_y ),   ℓ_x := Σ_{z ∈ B_R(x)∖B_R(y)} a_x(z)².

If z ∈ B_R(x)∖B_R(y), then d(x,z) ≤ R < d(y,z) ≤ d(x,z) + 1, so d(x,z) = R. Hence ℓ_x ≤ |S_R(x)|·2^{−R}. By the
monotonicity in Lemma 3.1, R·|S_R(x)|2^{−R} ≤ Σ_{j=1}^R |S_j(x)|2^{−j} = |a_x|² − 1, so ℓ_x ≤ (|a_x|² − 1)/R.
The same holds for y. Therefore

    ⟨a_x, a_y⟩ ≥ (√2/3)(1 − 1/R)(|a_x|² + |a_y|²) ≥ (√2/3)(1 − 1/R)·2|a_x||a_y|,

which is the first claim. If B_R(x) = B_R(y) = V, then ℓ_x = ℓ_y = 0 and the factor (1 − 1/R) disappears. ∎

### 3.2 Edge vectors

For an edge e = xy put b_e := (â_x + â_y)/|â_x + â_y|. This is well defined because
|â_x + â_y|² = 2 + 2⟨â_x,â_y⟩ ≥ 2.

**Lemma 3.3 (edge vectors are distinct)** [P]. If e ≠ f, then b_e ≠ b_f.
*Proof.* Let e = xy, α = 1/|a_x| and β = 1/|a_y|. The coordinates of the nonnegative, nonzero vector α a_x + β a_y, a
positive multiple of b_e, are α + β 2^{−1/2} at x and α 2^{−1/2} + β at y. At every z ∉ {x, y}, both d(x,z) ≥ 1 and
d(y,z) ≥ 1, so the coordinate is at most (α + β)2^{−1/2}, which is strictly less than both values above. So {x, y}
is the unique 2-set whose coordinates in b_e all exceed every other coordinate. Uniqueness holds because two
distinct such sets would each contain a coordinate exceeding a coordinate of the other. Thus b_e determines e,
since H is simple. ∎

**Lemma 3.4 (star)** [P]. Let v have edges e₁, e₂, e₃ and write b_i := b_{e_i}. Then, for i ≠ j,

    ⟨b_i, b_j⟩ ≥ ρ_R ≥ 0,   and   ⟨b₁,b₂⟩ + ⟨b₁,b₃⟩ + ⟨b₂,b₃⟩ ≥ (3/4)(1 + 3ρ_R).

*Proof.* Let e_i = v x_i and λ_i := ⟨â_v, â_{x_i}⟩ ≥ ρ_R (Lemma 3.2). Then

    c_i := ⟨b_i, â_v⟩ = (1 + λ_i)/√(2 + 2λ_i) = √((1 + λ_i)/2) ≥ c := √((1 + ρ_R)/2) ∈ [0, 1].

*Pairs.* Write b_i = c_i â_v + w_i with w_i ⊥ â_v and |w_i|² = 1 − c_i². By Cauchy–Schwarz,
⟨b_i,b_j⟩ = c_i c_j + ⟨w_i,w_j⟩ ≥ c_i c_j − √((1−c_i²)(1−c_j²)) ≥ c² − (1 − c²) = 2c² − 1 = ρ_R.

*Sum.* With s := b₁ + b₂ + b₃, Cauchy–Schwarz gives |s|² ≥ ⟨s, â_v⟩² = (c₁ + c₂ + c₃)² ≥ 9c² = 9(1 + ρ_R)/2.
Also |s|² = 3 + 2Σ_{i<j}⟨b_i,b_j⟩. So Σ_{i<j}⟨b_i,b_j⟩ ≥ (9(1+ρ_R)/2 − 3)/2 = (3/4)(1 + 3ρ_R). ∎

### 3.3 Scores, the median ordering, and the charging lemma

Let ω be a standard Gaussian vector in R^V, and define the **score** of an edge e by k_e := ⟨b_e, ω⟩, a standard
normal by (G0). For a vertex v let ℓ_v ≤ m_v ≤ u_v be its three scores sorted in nondecreasing order (minimum,
median, maximum); these are defined for every ω. On the event

    Ω₀ := {the 3n/2 scores are pairwise distinct},   P(Ω₀) = 1

(by Lemma 3.3 and Lemma 1.3(b)), the inequalities are strict: ℓ_v < m_v < u_v. The **median ordering** σ_ω
lists the vertices by nondecreasing m_v, with ties broken by a fixed labelling of V. Call v **open at t** if
ℓ_v < t < u_v, and let

    O(t) := #{v ∈ V : ℓ_v < t < u_v}   (t ∈ R).

**Lemma 3.5 (charging; ties cost +1)** [P]. On Ω₀, cw(σ_ω) ≤ 1 + max_{t∈R} O(t).

*Proof.* Let P be a proper nonempty prefix of σ_ω, u its last vertex, and t := m_u. Let e = xy cross P with x ∈ P
and y ∉ P. Since σ_ω is sorted by medians, m_x ≤ t ≤ m_y. There are three cases.
- k_e > t. Since k_e is a score of x and k_e > t ≥ m_x, k_e is the maximum of x's three distinct scores. So
  ℓ_x < m_x ≤ t < k_e = u_x: x is open at t, and e is the **only** edge at x with score > t. Charge e to x.
- k_e < t. Symmetrically k_e = ℓ_y < t ≤ m_y < u_y: y is open at t, and e is the only edge at y with score < t.
  Charge e to y.
- k_e = t. On Ω₀ at most one edge has score exactly t.

The charging is injective. Two edges charged to the same v ∈ P are both its unique edge with score > t, so they
are equal. Two edges charged to the same v ∉ P are both its unique edge with score < t, so they are equal. A vertex
cannot be in P and outside P. Hence cut(P) ≤ O(t) + 1. O takes integer values in [0, n], so the maximum
exists. ∎

### 3.4 One threshold

**Lemma 3.6** [P]. For every vertex v and every t ∈ R, P[v open at t] ≤ κ(ρ_R). Hence E O(t) ≤ κ(ρ_R)·n.

*Proof.* Let X₁, X₂, X₃ be the scores of the edges e₁, e₂, e₃ at v, and r_{ij} := ⟨b_i, b_j⟩ = Cov(X_i, X_j).
By Lemma 3.4, r_{ij} ∈ [ρ_R, 1] ⊆ [0, 1]. By Lemma 3.3, r_{ij} < 1, since unit vectors with inner product 1 are
equal. Fix t. With probability 1, t differs from X₁, X₂, X₃ (Lemma 1.3(a)), and on that event

    1{min_i X_i < t < max_i X_i} = (1/2)·Σ_{i<j} 1{t strictly between X_i and X_j}.                (3.2)

Indeed, if t is below or above all three, both sides are 0. Otherwise t splits the three values into nonempty
groups of sizes 1 and 2, exactly the 2 pairs with one member in each group are separated, and the left side is 1.
Take expectations in (3.2) and apply Lemma 1.2 to each pair, then Lemma 1.4(c), Lemma 3.4 and Lemma 1.4(a):

    P[v open at t] ≤ (1/(2π)) Σ_{i<j} arccos r_{ij} ≤ (3/(2π)) arccos((r₁₂ + r₁₃ + r₂₃)/3)
                   ≤ (3/(2π)) arccos((1 + 3ρ_R)/4) = κ(ρ_R).  ∎

### 3.5 All cut positions at once

A bound on E O(t) for each fixed t does **not** bound E max_t O(t), and the cut positions of σ_ω are random. We
therefore use a union bound over finitely many **fixed** thresholds. Fix T > 0 and an integer M ≥ 2. Put
τ_j := −T + 2T(j − 1)/(M − 1) for 1 ≤ j ≤ M (so τ₁ = −T, τ_M = T, and the spacing is 2T/(M−1)), and τ₀ := −∞,
τ_{M+1} := +∞. For 0 ≤ j ≤ M define

    N_j := #{v : ℓ_v < τ_{j+1} and u_v > τ_j}.

**Lemma 3.7 (sandwich and means)** [P].
(i) For every ω, max_{t∈R} O(t) ≤ max_{0≤j≤M} N_j.
(ii) E N₀ ≤ 3nΦ(−T) and E N_M ≤ 3nΦ(−T).
(iii) E N_j ≤ κ(ρ_R)·n + 6T·n/((M − 1)√(2π)) for 1 ≤ j ≤ M − 1.

*Proof.* (i): Every t ∈ R lies in [τ_j, τ_{j+1}] for some 0 ≤ j ≤ M, with the convention (−∞, τ₁] and [τ_M, ∞)
at the ends. If v is open at t, then ℓ_v < t ≤ τ_{j+1} and u_v > t ≥ τ_j, so v is counted in N_j. Hence
O(t) ≤ N_j.
(ii): N₀ counts v with ℓ_v < −T, so E N₀ ≤ Σ_v Σ_{e∋v} P[k_e < −T] = 3nΦ(−T). N_M is symmetric.
(iii): If v is counted in N_j, then either ℓ_v < τ_j, in which case ℓ_v < τ_j < u_v and v is open at τ_j, or
τ_j ≤ ℓ_v < τ_{j+1}. So N_j ≤ O(τ_j) + #{v : some score at v lies in [τ_j, τ_{j+1})}. Lemma 3.6 and
Lemma 1.3(a) give E N_j ≤ κ(ρ_R)n + 3n·(2T/(M−1))/√(2π). ∎

**Lemma 3.8 (locality and variance)** [P]. For each j, N_j = Σ_v X_v, where X_v ∈ {0,1} is a function of
(ω_z)_{z ∈ B_{R+1}(v)} only, and Var N_j < 3n·4^R.

*Proof.* For an edge e = vx, b_e is supported on B_R(v) ∪ B_R(x) ⊆ B_{R+1}(v). So the three scores at v, and
hence X_v, depend only on ω restricted to B_{R+1}(v). If d(u,v) ≥ 2R + 3, then B_{R+1}(u) ∩ B_{R+1}(v) = ∅, so
X_u and X_v are independent and Cov(X_u, X_v) = 0. Otherwise |Cov| ≤ 1/4 (§1.0). For each v there are at most
|B_{2R+2}(v)| ≤ 3·2^{2R+2} − 2 vertices u with d(u,v) ≤ 2R+2 (Lemma 3.1). Therefore
Var N_j = Σ_{u,v} Cov(X_u, X_v) ≤ n(12·4^R − 2)/4 < 3n·4^R. ∎

**Proposition 3.9** [P]. For all integers R ≥ 1, M ≥ 2 and every real T > 0 with 3Φ(−T) ≤ κ(ρ_R),

    cw(H) ≤ 1 + n·[ κ(ρ_R) + 6T/((M − 1)√(2π)) ] + √(3(M + 1)·n·4^R).

*Proof.* Let λ := √(3(M+1)n4^R). By Chebyshev and Lemma 3.8, P[N_j ≥ E N_j + λ] ≤ Var N_j/λ² < 1/(M+1) for each
of the M + 1 indices j. By the union bound, the event "N_j < E N_j + λ for all j" has probability > 0, and so does
its intersection with Ω₀. Fix ω in that intersection. By Lemmas 3.5 and 3.7,

    cw(H) ≤ cw(σ_ω) ≤ 1 + max_t O(t) ≤ 1 + max_j N_j < 1 + max_j E N_j + λ
          ≤ 1 + n·max{ κ(ρ_R) + 6T/((M−1)√(2π)), 3Φ(−T) } + λ,

and the maximum is the first entry by hypothesis. ∎

### 3.6 Parameters: proof of Theorem 2

**Lemma 3.10** [P]. For R ≥ 1, κ(ρ_R) ≤ γ + S/R, where S := 3√2/(π√(7 − 4√2)) < 1.16527.
*Proof.* x(ρ) := (1 + 3ρ)/4 maps [0, r*] into [1/4, x₀] ⊂ [0, 1), and x₀ − x(ρ_R) = 3r*/(4R) = 1/(√2 R). Also
1 − x₀² = (16 − (9 + 4√2))/16 = (7 − 4√2)/16. By Lemma 1.4(b) with m = x₀,

    arccos x(ρ_R) ≤ θ₀ + (1/(√2R))·4/√(7 − 4√2).

Multiplying by 3/(2π) gives κ(ρ_R) ≤ γ + 3√2/(π√(7−4√2)R). The numerical bound on S is (W4) in §5. ∎

*Proof of Theorem 2.* Let n ≥ 32 and L := log₂n ≥ 5. Choose

    R := ⌊L/4⌋ ≥ 1,   M := 4R + 1,   T := 2.

The hypothesis of Proposition 3.9 holds: 3Φ(−2) < 0.0813 < γ ≤ κ(ρ_R) by (W5). Next,
6T/((M−1)√(2π)) = 3/(R√(2π)) < 1.19683/R by (W4). With Lemma 3.10, the bracket in Proposition 3.9 is at most
γ + (1.16527 + 1.19683)/R < γ + 2.3622/R. Since R > L/4 − 1 = (L − 4)/4 > 0, we get 2.3622/R < 9.4488/(L − 4) ≤
9.45/(L − 4). Finally 4R + 2 ≤ L + 2, and 4^R ≤ 2^{L/2} = n^{1/2}, so √(3(M+1)n4^R) ≤ √(3L + 6)·n^{3/4}.
Proposition 3.9 then gives Theorem 2. ∎

## 4. Proof of Theorem 1

Let K be a connected cubic multigraph on D vertices, and recall h(D) = 24D/log₂(D+2).

**(a) D < 2¹⁶.** By (F1), cw(K) ≤ 3D/2. Since D + 2 < 2¹⁷, (3/2 − γ)·log₂(D+2) < (3/2 − γ)·17 < 24 by (W6), so
3D/2 < γD + h(D).

**(b) D ≥ 2¹⁶.** Put L := log₂D ≥ 16. Corollary 2.3 leaves three cases.
- (b1) cw(K) ≤ 3L + 3.
- (b2) Kernel H with n < 64: cw(H) ≤ 3n/2 ≤ 93 (n ≤ 62 is even), so cw(K) ≤ 96 + 3L.
- (b3) Kernel H with 64 ≤ n ≤ D: by Theorem 2, cw(H) ≤ γn + ε(n), where
  ε(x) := 9.45x/(log₂x − 4) + √(3log₂x + 6)·x^{3/4} + 1.
  ε is nondecreasing on [64, ∞). The second and third terms clearly are. For the first, d/dx[x/(log₂x − 4)]
  has the sign of (log₂x − 4) − 1/ln 2, which is ≥ 2 − 1/ln 2 > 0 (W7). With γn ≤ γD this gives
  cw(K) ≤ γD + ε(D) + 3L + 3.

Since ε(D) ≥ 9.45·2¹⁶/12 > 93, all three cases give cw(K) ≤ γD + ε(D) + 3L + 3. It remains to show
ε(D) + 3L + 3 ≤ h(D). As D + 2 ≤ 2D, log₂(D+2) ≤ L + 1, so it suffices that

    Ψ(L) := (L + 1)·[ 9.45/(L − 4) + √(3L + 6)·2^{−L/4} + (3L + 4)·2^{−L} ] ≤ 24   for all real L ≥ 16,

using D^{3/4}/D = 2^{−L/4} and (1 + 3L + 3)/D = (3L+4)2^{−L}. Each of the three products is nonincreasing on
[16, ∞):
- (L+1)/(L−4) = 1 + 5/(L−4).
- The logarithmic derivative of (L+1)√(3L+6)2^{−L/4} is 1/(L+1) + 1/(2L+4) − (ln 2)/4 ≤ 1/17 + 1/36 − (ln 2)/4 < 0.
- The logarithmic derivative of (L+1)(3L+4)2^{−L} is 1/(L+1) + 3/(3L+4) − ln 2 ≤ 1/17 + 3/52 − ln 2 < 0.
(Both by (W7).) So Ψ(L) ≤ Ψ(16) = 17·(9.45/12 + √54/16 + 52/2¹⁶) < 21.21 < 24 (W8). ∎

## 5. Numerical inequalities, with rational witnesses

Each inequality below is a comparison of rational numbers after the stated rational bounds on π, √2, e and ln 2
are inserted; monotonicity of the inserted bounds is used in the obvious direction each time. The reader may
verify each by direct computation from the intermediate values given (decimals are rounded in the safe
direction).

- **(W1)** 3.14159265358 < π < 3.14159265359. Machin's formula π = 16 arctan(1/5) − 4 arctan(1/239): for
  α = arctan(1/5), the tangent doubling formula gives tan 2α = 5/12 and tan 4α = 120/119, and the subtraction
  formula gives tan(4α − arctan(1/239)) = (120·239 − 119)/(119·239 + 120) = 28561/28561 = 1; as
  0 < 4α − arctan(1/239) < 4/5 < π/2, it equals π/4. For 0 < x < 1 the partial sums
  s_n(x) := Σ_{k=0}^{n} (−1)^k x^{2k+1}/(2k+1) alternate around arctan x (the terms decrease), so
  arctan(1/5) ∈ (s_7(1/5), s_8(1/5)) ⊂ (0.1973955598498063, 0.1973955598498835) and
  arctan(1/239) ∈ (s_3(1/239), s_2(1/239)) ⊂ (0.0041840760020747, 0.0041840760020748). Hence
  π > 16·0.1973955598498063 − 4·0.0041840760020748 = 3.1415926535886016 and
  π < 16·0.1973955598498835 − 4·0.0041840760020747 = 3.1415926535898372. Also
  1.4142135623730² = 1.999999999999731161391129 < 2 < 2.00000000000001400410360361 = 1.4142135623731²;
  e > Σ_{k=0}^{7} 1/k! = 685/252 > 2.718; and 0.693 < ln 2 < 0.6932, from ln 2 = ln((1 + 1/3)/(1 − 1/3)) =
  2Σ_{k≥0} 3^{−(2k+1)}/(2k+1): the terms are positive and each is less than 1/9 of the previous one, so the
  partial sum up to k = 3 (> 0.69313) is a lower bound and adding 9/8 of the term k = 4 (< 0.0000128) gives an
  upper bound.
- **(W2)** Let a := 0.2939503814 and a′ := 0.2939503813. For 0 < s ≤ 1 the Taylor partial sums of cos s
  alternate around cos s with decreasing terms. Then cos a ≤ Σ_{k=0}^{8} (−1)^k a^{2k}/(2k)! < 0.9571067811654
  < 0.95710678118650 = (1 + 2·1.4142135623730)/4 < x₀, and cos a′ ≥ Σ_{k=0}^{7} (−1)^k a′^{2k}/(2k)!
  > 0.9571067811943 > 0.95710678118655 = (1 + 2·1.4142135623731)/4 > x₀. Since cos decreases on [0, π],
  a′ < θ₀ < a.
- **(W3)** γ = 3θ₀/(2π) ∈ (3a′/(2·3.14159265359), 3a/(2·3.14159265358)), and 3a′/(2·3.14159265359)
  > 0.1403509686229, 3a/(2·3.14159265358) < 0.1403509686712; so 0.1403509686 < γ < 0.1403509687.
- **(W4)** S² = 18/(π²(7 − 4√2)) < 18/(3.14159265358²·(7 − 4·1.4142135623731)), where 3.14159265358²
  > 9.869604401027 and 7 − 4·1.4142135623731 = 1.3431457505076 with product > 13.25631721043; so
  S² < 18/13.25631721043 < 1.3578431864 < 1.3578541729 = 1.16527². Also (3/√(2π))² = 9/(2π) < 9/6.28318530716
  < 1.4323944879 < 1.4324020489 = 1.19683². Hence S + 3/√(2π) < 1.16527 + 1.19683 = 2.3621 < 2.3622 and
  4·2.3622 = 9.4488 ≤ 9.45.
- **(W5)** 3Φ(−2) ≤ 3φ(2)/2 = 3e^{−2}/(2√(2π)) < 3/(2.718²·2·2.5) = 3/36.93762 < 0.081219 < 0.0813 < 0.1403
  < γ. Uses Lemma 1.3(c), (W1), (W3) and 2π > 6.25.
- **(W6)** (3/2 − γ)·17 < (3/2 − 0.1403509686)·17 = 23.1140335338 < 24, by (W3).
- **(W7)** 2 − 1/ln 2 > 2 − 1/0.693 > 2 − 1.4431 > 0; 1/17 + 1/36 = 53/612 < 0.0867 < 0.17325 = 0.693/4
  < (ln 2)/4; and 1/17 + 3/52 = 103/884 < 0.1166 < 0.693 < ln 2.
- **(W8)** √54 < 7.3485 since 7.3485² = 54.00045225; so Ψ(16) = 17·(9.45/12 + √54/16 + 52/2¹⁶)
  < 17·(0.7875 + 0.45928125 + 0.00079346) = 21.20877007 < 21.21.
