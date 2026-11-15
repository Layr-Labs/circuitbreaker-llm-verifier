# Complete proof of the joint theorem: C_B2(f_n) >= (114062441973/25000000000)*n - o(n)

This file is the master proof for the joint theorem declared in `bound.json`.
The theorem reproduces, with full supplied-source evidence, the wire-graph
route of McGuire, Schlesinger and Williams ("Explicit functions requiring
circuits larger than 4.6n", preprint, November 2026; public appendix
repository `github.com/rrwilliams0/circuit-lb-4.6155n`, commit
`858a9b6875a2cdab63f711cb2faf2e04326d22fb`, whose complete `proof.md` this
file transcribes and compresses for its Sections 1-3), instantiated with
the *edge-score* constant of that route:

    kappa_E := (3/(2 pi)) * arccos(mu0),   mu0 := (1 + 2 sqrt 2)/4,

    liminf_{n -> inf} s_n / n  >=  1 + 1/(2 kappa_E)  >=  4.56249767892.


Every lemma below is proved here in full, with two classes of exception,
each identified where used: (i) a short list of classical textbook facts
each used in one line and each carrying a standard name (Zorn's lemma,
the existence of a splitting field of a polynomial over a field
(obtained by quotienting out one irreducible factor at a time), Gauss's
theorem that a polynomial ring over a UFD is a UFD, the Chinese
remainder theorem for pairwise comaximal ideals, the unique
factorization of a one-variable polynomial ring over a field (the
Euclidean algorithm), Petersen's matching theorem, the Levy continuity
theorem and the portmanteau theorem); the Hilbert Nullstellensatz and
Zariski's lemma, by contrast, are PROVED in `krs.md` (Section 4.1);
(ii) the graph-theoretic reduction and
matched-insertion framework, which reproduce the preprint [MSW] of
McGuire, Schlesinger and Williams (credited, traced section by section;
every proof is re-proved here). The hard family and its counting core
are fully self-contained: the norm-one rectangle bound (Lemma HARD via
Sublemma NORM) reduces to the Kollar-Ronyai-Szabo theorem on systems of
product equations, which `krs.md` proves from first principles
(integral extensions, valuation rings via Zorn, Nakayama, localization,
and the tools just listed). All arithmetic constants are certified by
the supplied script `checks.py`, whose complete transcript
`checks-output.txt` is supplied. The source's stronger constant 4.6155...
(their Theorem M_tree: the *exact* value of the cone supremum M_123,
shown equal to a Gaussian functional M_tree at the 3-regular tree point)
is not used or claimed here; this file needs only the source's rule bound
(Lemma RB below), which caps the cone supremum by 2*kappa_E.

## 0. The model, the family, and the theorem

**Model.** A B_2-circuit on input variables x_1, ..., x_n is a directed
acyclic graph whose sources are the variables and whose other vertices,
the *gates*, each have two ordered predecessors (possibly equal) and
compute an arbitrary Boolean function of them. One gate is the output.
Fan-out is unrestricted and subcircuits may be shared. The *size* s is
the number of gates. This is exactly the pinned full B_2 model of
`rules.md`: two input pins per gate, all sixteen binary truth tables,
unrestricted depth and sharing, exact computation on every input; inputs,
wires and the output selection cost nothing; each retained gate,
including any trivial or degenerate gate, costs one.

**The hard function (a norm-one indicator family).** The family is a
total explicit family defined from finite fields; its rectangle
property (Lemma HARD below) rests on the Kollar-Ronyai-Szabo norm-graph
counting theorem, which is proved in full in the supporting file `krs.md`. (The choice of a
norm-one indicator as the hard function for this circuit route is due to
ArthurOstapenko's concurrent accepted entry 808da37; our layout,
parameters and proofs below are independent.)

Fix an input length N. Put

    t := 2^floor(log2(N)/3),   m := the largest odd m >= 1 with
    (t+1)m + t <= N,   lambda := N - m - t - mt,

and if no such m exists set f_N := 0 identically. Note t is a power of
2, m is odd, so gcd(m, t) = 1; t <= N^(1/3) and m <= (N-t)/(t+1) <= N/t
and lambda < 2(t+1). An input x in {0,1}^N is parsed as four
consecutive blocks: the *Q-block* (m bits: the coefficients of
X^(m-1), ..., X^0 of a monic Q in F_2[X] of degree exactly m), the
*R-block* (t bits: likewise a monic R in F_2[Y] of degree exactly t),
the *z-block* (mt bits z_ij, i < m, j < t, the coefficients of
z(X,Y) := sum z_ij X^i Y^j), and the *leftover* (lambda bits, ignored).
With

    e := (2^(mt) - 1)/(2^m - 1) = 1 + 2^m + 2^(2m) + ... + 2^(m(t-1)),

define f_N(x) := 1 iff (i) Q is irreducible over F_2, (ii) R is
irreducible over F_2, and (iii) z(X,Y)^e = 1 in F_2[X,Y]/(Q,R);
otherwise f_N(x) = 0. Condition (iii) is well-defined because (i)-(ii)
make the ideal (Q, R) maximal (Fact FIELD, `field.md`), so the quotient is a
field A_N with 2^(mt) elements in which z^e is computed. Throughout we
write H_{Q,R} := {w in A_N : w^e = 1} for the set of e-th roots of
unity; |H_{Q,R}| = e, which is odd, and e >= 2^(mt-m) (Fact UNIT,
`field.md`).

*The evaluator and NP membership.* The fixed evaluator parses the
blocks, tests Q and R for irreducibility by the standard deterministic
test (Fact IRRED, `field.md`: X^(2^d) = X mod h and gcd(h, X^(2^(d/p)) - X)
= 1 for each prime p | d), and, if both pass, computes z^e in A_N by
binary exponentiation, doing all arithmetic on bivariate polynomials
reduced modulo Q in X and modulo R in Y. Each step is polynomial in N,
the family is total and uniform at every length, so L := {x :
f_{|x|}(x) = 1} is in P; the NP verifier uses the empty witness
(bound 0). Witness bits never enter the circuits; no per-length advice
occurs; the language, the evaluator and the verifier are fixed once and
for all.

**Main Theorem.** Every B_2-circuit computing f_n has, for every eps > 0
and all large n, at least (1 + 1/(2 kappa_E) - eps) n gates; in
particular

    liminf_{n -> inf} s_n/n >= 1 + 1/(2 kappa_E) >= 4.56249767892
    = 114062441973/25000000000,

the declared coefficient. (The enclosure of kappa_E and the final exact
comparison are certified in `checks.py`, parts 2 and 3; the assembly of
the theorem is Section 3 below.)

*Proof strategy.* (1) A circuit with s gates reading b variables has a
wire graph of cycle rank s - b + 1 (Section 1.3), which compresses to a
cubic *kernel* on at most 2(s - b) vertices at logarithmic cutwidth cost
(Section 1.4). (2) An ordering of the wire graph of width w yields an
unambiguous read-once program with about 2^w edges (Lemma SIM), while
the norm-one family admits no such program with fewer than
2^{N - o(N)} edges (Lemmas HARD and COUNT); together, a kernel cutwidth
bound cw(H) <= (kappa + o(1)) h forces s >= (1 + 1/(2 kappa)) N - o(N)
(Theorem K2C). (3) Every kernel has cutwidth at most (M_123/2 + o(1)) h,
where M_123 is the supremum of an explicit Gaussian functional over a
123-type cone (Section 2.5, Theorem A); and the rule bound (Lemma RB)
caps M_123 by 2 kappa_E with kappa_E = (3/(2 pi)) arccos((1 + 2 sqrt
2)/4) (Section 2.4). Sections 0 and 1.1 supply the hard family and its
rectangle property, whose counting core (Theorem KRS) is proved in
`krs.md`; Section labels [MSW Sx] trace the graph-theoretic parts to
the source; every proof is included.

## 1. From circuits to kernels [MSW S2.1-2.4]

### 1.1 Rectangles, the rectangle bound, and the counting lemma

The five finite-field facts used here (LINPOLY, FROB, KER, IRRED, UNIT,
FIELD: root counting, Frobenius order, the irreducibility test and
existence, cyclicity and root-of-unity counts, and maximality of
composite ideals) are stated and proved in the supporting file
`field.md`, to keep every file of this package under its size limit.
**Rectangles.** Split the coordinates {1, ..., N} into disjoint sets U,
V. For P in {0,1}^U and Q in {0,1}^V the *rectangle* P x Q is the set
of inputs whose U-part lies in P and V-part in Q; its sides are P and
Q. It is *accepting* for f if f = 1 on all of it.

**Sublemma NORM (fixed-field rectangle bound).** Let F_0 be a field
with 2^(mt) elements (m, t as in the family; such fields exist by Facts
IRRED(b) and FIELD), let phi: w -> w^2 be its Frobenius, which by Fact
FROB has order mt, so sigma := phi^m : w -> w^(2^m) is an automorphism
of order exactly t. Let H_0 = {w in F_0 : w^e = 1} with
e = sum_{j=0}^{t-1} 2^(mj), and let A, B ⊆ F_0 satisfy
A + B ⊆ H_0. If |A| >= t then |B| <= t!.

*Proof.* Pick distinct a_1, ..., a_t in A. For each b in B and each i:
(a_i + b)^e = 1. Expand e = sum_{j=0}^{t-1} 2^(mj) and use the
freshman's-dream additivity of w -> w^(2^(mj)) (Fact FROB):

    (a_i + b)^e = prod_{j=0}^{t-1} (a_i + b)^(2^(mj))
                = prod_{j=0}^{t-1} (a_i^(2^(mj)) + b^(2^(mj))).

Introduce unknowns y_j := b^(2^(mj)) in F_0 and constants
c_ij := a_i^(2^(mj)): each b in B yields a solution (y_0, ..., y_(t-1))
of the system prod_j (y_j + c_ij) = 1, i = 1, ..., t. For each column
j the map w -> w^(2^(mj)) is a bijection (a power of sigma), so
i -> c_ij is injective: the columns are distinct in the sense of
Theorem KRS (`krs.md`). Theorem KRS bounds the number of solutions by
t!, and b -> (y_j)_j is injective since y_0 = b. So |B| <= t!. qed

**Lemma HARD.** For all N large enough that t = t(N) >= 2 and
3m + 2t + lambda + t·ceil(log2 t) + 2 <= N: (a) f_N has at least
2^(N - 2m - t) accepting inputs; (b) with K := 2^(m+t+lambda+1)·t!,
every accepting rectangle of f_N, under every split of the coordinates,
has a side of size < K; (c) f_N depends on each of its mt z-block
variables.

*Proof.* (a) The accepting inputs are exactly, for valid (Q, R), the
patterns with z(X,Y)^e = 1 in A_N; z -> its class in A_N is a bijection
from z-patterns to A_N (unique representative with X-degree < m,
Y-degree < t), so for each valid pair the count is |H_{Q,R}| = e >=
2^(mt-m) (Fact UNIT). With at least one valid pair (Fact IRRED(b)) and
2^lambda free leftover bits:

    #accepting >= 2^lambda · 2^(mt-m) = 2^(N - m - t - m) = 2^(N-2m-t).

(b) Let P x Q be an accepting rectangle under the split U, V. Let
Prefix denote the m + t coordinates of the Q-block and R-block. For an
assignment a of the coordinates in U ∩ Prefix, let P_a := {p in P :
p restricted to U ∩ Prefix is a}; the P_a (a in the image, at most
2^(m+t) of them) partition P. Likewise Q = disjoint union of Q_b over
assignments b of V ∩ Prefix, at most 2^(m+t) of them. Fix any a, b with
P_a, Q_b both nonempty. Every combined input from P_a x Q_b has the
same prefix blocks: the Q-block string is a's U-part completed by b's
V-part, and similarly the R-block; call these Q_[a,b], R_[a,b]. All
pairs are accepting, so Q_[a,b], R_[a,b] are valid and F_0 :=
F_2[X,Y]/(Q_[a,b], R_[a,b]) is a field with 2^(mt) elements (Fact
FIELD). The z-block of a combined input has U-coordinates from p and
V-coordinates from q; since z(X,Y) is F_2-linear in its coefficient
pattern, the class of z(p,q) in F_0 equals alpha_p + beta_q where
alpha_p (resp. beta_q) is the class of the pattern carried by p's
U-coordinates (resp. q's V-coordinates) alone. Acceptance of every pair
says (alpha_p + beta_q)^e = 1. So with A := {alpha_p : p in P_a} and
B := {beta_q : q in Q_b} we have A + B ⊆ H_0. Let lambda_U, lambda_V be
the numbers of leftover coordinates in U, V (lambda_U + lambda_V =
lambda). An element of P_a is determined by its alpha and its
U-leftover (two elements of P_a with the same alpha and the same U-leftover
have the same U z-pattern - the class of a pattern determines the
rectangular pattern itself, the classes of the monomials X^i Y^j,
i < m, j < t forming the F_2-basis of (a), so distinct patterns give
distinct classes - and the U z-pattern together with the U-prefix
pattern (fixed in P_a) and the U-leftover determines the element,
hence they are equal), so |P_a| <= |A| · 2^lambda_U; likewise |Q_b| <= |B| ·
2^lambda_V. Apply Sublemma NORM: either |A| >= t, and then |B| <= t!,
whence |Q_b| <= t!·2^lambda_V; or |A| <= t - 1, whence |P_a| <=
(t-1)·2^lambda_U. In both cases min(|P_a|, |Q_b|) <= t!·2^lambda (in
the first case min <= |Q_b| <= t!·2^(lambda_V) <= t!·2^lambda; in the
second min <= |P_a| <= (t-1)2^(lambda_U) <= t!·2^lambda as t <= t!).
This holds for every pair (a, b) of nonempty cells. Therefore: either
|P_{a*}| > t!·2^lambda for some a* — then every cell Q_b with
Q_b nonempty must have |Q_b| <= t!·2^lambda (otherwise the pair
(a*, b) contradicts the previous paragraph, as
min(|P_{a*}|, |Q_b|) >= min(t!·2^lambda + 1, t!·2^lambda + 1) >
t!·2^lambda), so |Q| <= 2^(m+t)·t!·2^lambda; or all cells are small,
|P_a| <= t!·2^lambda for all a, and |P| <= 2^(m+t)·t!·2^lambda. Either
way min(|P|, |Q|) <= 2^(m+t+lambda)·t! < K.

(c) Fix a valid pair (Q, R) (Fact IRRED(b)) and let c be a z-block
coordinate, v the z-pattern with a single 1 at c. The classes of the
monomials X^i Y^j (i < m, j < t) form an F_2-basis of A_N (A_N =
E[Y]/(R) with E = F_2[X]/(Q); the same basis underlies the pattern-to
-element bijection used in (a) and (b)), so vbar := the class of v is a
nonzero element. H_{Q,R} and its translate H_{Q,R} + vbar have equal,
odd cardinalities (Fact UNIT), so they differ (equality would partition
H_{Q,R} into pairs {w, w + vbar}, an even count): pick z in H_{Q,R}
with z + vbar not in H_{Q,R}. The two inputs (Q, R, z, 0) and
(Q, R, z XOR v, 0) differ only in coordinate c and have f = 1 and
f = 0. So every z-block variable is essential, and f_N depends on at
least mt = N - (m + t + lambda) variables. qed

### 1.2 Read-once programs and the simulation

A *read-once program* on x_1, ..., x_n has pairwise disjoint vertex sets
V_0, ..., V_l, V_0 a single *start* and V_l a single *accept* vertex,
and directed edges each from some V_{i-1} to V_i (the i-th *layer*). A
layer is a *query layer*, tagged with a variable x_j, whose edges carry
labels 0 or 1, or a *guess layer* with unlabelled edges. Every variable
tags exactly one query layer. Parallel edges count separately; the
*size* is the number of edges. A start-to-accept path is *consistent*
with an input a if each query edge tagged x_j carries the label a_j;
the program *accepts* a if some consistent path exists; it *computes* f
if it accepts exactly the inputs with f = 1; it is *unambiguous* if
every accepted input has exactly one consistent path. The *past* of V_i
is the set of variables tagging layers 1..i; the *future* is the rest.

**Lemma COUNT.** Let f have exactly M >= K > 1 accepting inputs and let
every accepting rectangle of f have a side of size < K. Then every
read-once program computing f has more than M/K^2 edges.

*Proof.* For a vertex v in V_i let P_v be the set of assignments to the
past of V_i that admit a consistent path from start to v, and Q_v the
set of assignments to the future that admit a consistent path from v to
accept. The two paths constrain disjoint variable sets and concatenate,
so P_v x Q_v is an accepting rectangle, whence min(|P_v|, |Q_v|) < K.
Fix one consistent path for each accepted input. Along it |P_v| goes
from 1 (< K) at the start to M (>= K) at the accept vertex (whose past
is everything and whose P is the set of accepted inputs). Let e = (u, v)
be the first edge of the path with |P_u| < K <= |P_v|, and charge the
input to e; by the above |Q_v| < K. An input charged to e is determined
by its assignment on the past of u (an element of P_u), its assignment
on the future of v (an element of Q_v), and, if e lies in a query layer,
the value of the queried variable, which is the label of e (in a guess
layer past and future already exhaust the variables). So fewer than K^2
inputs are charged to any one edge, and M < |E| K^2. qed

By Lemma HARD(a),(b) and Lemma COUNT, for all large N every read-once
program computing f_N has more than 2^{N-2m-t}/K^2 edges. Unwinding
K = 2^{m+t+lambda+1}·t! and t! <= t^t (so log2 t! <= t·ceil(log2 t)),
this is more than

    2^{N - 4m - 3t - 2 lambda - 2 t ceil(log2 t) - 2} =: 2^{N - L_2(N)}
                                                                        (1.1)

edges, with L_2(N) = o(N) since t <= N^(1/3), m <= N/t and
lambda < 2(t+1) (each term is O(N^(2/3)) or O(N^(1/3) log N)).

### 1.3 Wire graphs

A variable or gate is a *signal*; a signal is *read* if it has a
directed path to the output. *Pruning* deletes everything unread; the
pruned circuit computes the same function; assume it has at least one
gate. Route every signal with F >= 2 uses through F - 1 *copy vertices*,
each with one incoming and two outgoing wires. The *wire graph* G of a
pruned circuit with t gates and b read variables is the undirected
multigraph with one vertex per gate, read variable and copy vertex, and
one edge per wire. Each of the t + b - 1 non-output signals has at least
one use and there are 2t uses, so there are t + 1 - b copy vertices.
Counting edges into gates and copy vertices:

    N = 2t + 1 vertices,   M = 3t + 1 - b edges,   mu(G) = M - N + 1 = t - b + 1.   (1.2)

Gates have degree 3 except the output (degree 2), copy vertices degree
3, variables degree 1; G is connected (every vertex reaches the output)
and loopless (the circuit is acyclic); a gate reading one signal twice
gives a parallel pair. For an ordering u_1, ..., u_N of the vertices,
the *cut* F_i after u_i is the set of edges with exactly one endpoint
among u_1, ..., u_i; the *width* of the ordering is max_i |F_i|, and the
*cutwidth* cw(G) is the least width over all orderings.

**Lemma SIM.** Let G be the wire graph of a pruned circuit with at most
s gates on n variables, and let G have an ordering of width w. Then
there is an unambiguous read-once program on x_1, ..., x_n computing the
circuit's function with at most 4(2s + 1)2^w + 2n edges.

*Proof.* (Checks) Each edge carries one bit. A non-output gate computing
g requires z = g(a, b) for the bits a, b on its incoming wires (slot
order) and z on its outgoing wire; the output gate requires g(a, b) = 1;
a copy vertex requires its three bits to agree; a variable allows any
bit. For each input there is exactly one bit assignment to the edges
that puts each read variable's value on its edge and satisfies all
non-output and copy checks (evaluate the circuit in topological order);
it satisfies the output check iff the circuit accepts the input.

(The program) V_i := {(i, sigma) : sigma in {0,1}^{F_i}}; V_0 = {start}.
To process u_i let B be its edges to earlier vertices, O its edges to
later vertices, R = F_{i-1} \ B (so F_i = R u O). For every sigma in
{0,1}^{F_{i-1}} and omega in {0,1}^O such that sigma|_B u omega
satisfies the check of u_i, add an edge from (i-1, sigma) to (i,
sigma|_R u omega). If u_i is a gate or copy vertex this is a guess
layer. If u_i is the variable x_j it is a query layer tagged x_j, and
the edge carries the bit that sigma|_B u omega puts on u_i's single
edge. After layer N append one two-edge query layer (labels 0 and 1)
for each of the n - b unread variables, ending at accept.

(Correctness) A consistent start-to-accept path is exactly a bit
assignment to the edges of G satisfying every check and putting each
read variable's value on its edge: every edge is recorded when its
first endpoint is processed and checked when its second is. By the first
paragraph there is exactly one such path if the circuit accepts and
none otherwise; the program computes the function and is unambiguous.

(Size) A gate check has 4 satisfying triples, a copy check 2, the
output check at most 4; fixing bits only lowers these numbers. So each
vertex has at most 4 outgoing edges. Each V_i (i < N) has at most 2^w
vertices, there are N <= 2s + 1 such layers, and the final layers have
2 edges each: at most 4(2s+1)2^w + 2n edges. qed

### 1.4 Kernels and the reduction

Graphs here are finite loopless multigraphs; parallel edges count in
degrees. *Cubic*: every degree 3. *2-connected*: connected, at least
three vertices, and connected after deleting any one vertex. A *kernel*
is a simple, triangle-free, 2-connected cubic graph; h denotes its
number of vertices. A *bridge* is an edge whose deletion disconnects
the graph; an edge is a bridge iff it lies on no cycle (two parallel
edges form a cycle). *Contracting* a vertex set S replaces S by one
vertex, deletes the edges inside S and redirects the others.

**Lemma COMPRESS.** Let B be a connected bridgeless subcubic graph with
N_B vertices. Either B has an ordering of width at most
6 ceil(log2 N_B), or there is a kernel H on h <= 2(mu(B) - 1) vertices
such that, for every q, an ordering of H of width at most q yields an
ordering of B of width at most q + 6 ceil(log2 N_B).

*Proof.* If N_B = 1 the first alternative holds (width 0); let
N_B >= 2. Apply, while one applies, preferring the first: Rule 1: if
adjacent p != q, joined by r >= 1 edges, have deg p <= 2 or deg q <= 2
or r >= 2, contract {p, q}. Rule 2: if Rule 1 does not apply and a, b,
c are pairwise adjacent, contract {a, b, c}. Rule 1 creates a vertex of
degree deg p + deg q - 2r <= 3 and changes M - N by r - 1 >= 0. Under
Rule 2 every degree is 3 and the graph is simple (Rule 1 inapplicable),
so each of a, b, c has exactly one edge leaving the triangle, the new
vertex has degree 3, and M - N drops by 1. Each rule lowers N, so the
process stops; contraction keeps the graph connected, loopless and
bridgeless (a path avoiding an edge maps to a walk avoiding it). When it
stops, either one vertex remains, or the graph H is cubic (otherwise
Rule 1 would apply to a vertex of degree <= 2), simple (otherwise Rule
1 with r >= 2) and triangle-free (otherwise Rule 2); then h/2 = |E(H)|
- |V(H)| <= M_B - N_B = mu(B) - 1, and H is 2-connected: if deleting v
left components C_1, ..., C_m with m >= 2, each C_i sends at least two
edges to v (a single edge would be a bridge of H), so deg v >= 4.

(Pieces) Record the contractions as a rooted tree per current vertex
(leaves = original vertices, internal nodes with 2 or 3 children). The
leaf set of a subtree is a *piece*; a piece was a current vertex of
degree <= 3 when it was formed, and its set of leaving edges never
changes afterwards, so at most 3 edges of B leave a piece. (Expanding
one piece) List the leaves depth first, visiting children in
nonincreasing order of their leaf counts. A proper prefix L, with first
missing leaf l and root path x_0, ..., x_m = l, is the union over i < m
of the pieces of the children of x_i that precede x_{i+1}; if x_{i+1}
is the j-th child it contributes j - 1 <= 2 pieces, and j >= 2 forces
x_{i+1} to have at most half the leaves of x_i, so j >= 2 happens for
at most ceil(log2 N_B) indices; hence L is a union of at most
2 ceil(log2 N_B) pieces and at most 6 ceil(log2 N_B) edges leave it.
(Expanding an ordering) If one vertex remains, its leaf list is an
ordering of B of width <= 6 ceil(log2 N_B). Otherwise list the pieces
of the vertices of H in the given order, each expanded as above: a
prefix is a prefix of H's ordering (whole pieces) plus a proper prefix
of one more piece, so every cut has at most q + 6 ceil(log2 N_B) edges
(an edge leaving the prefix either crosses the corresponding cut of H
or leaves the partial piece). qed

The *blocks* of a connected subcubic G are the components of G with all
bridges deleted (a block may be a single vertex).

**Lemma ASSEMBLY.** Let G be connected subcubic, B a block, F_1, ...,
F_m (m >= 0) the components of G - V(B). (a) B is bridgeless; each F_j
is joined to B by exactly one edge e_j, ending at v_j in B; every edge
of G lies in B, in some F_j, or is some e_j; and mu(G) = mu(B) +
sum_j mu(F_j). (b) If B has an ordering of width w and each F_j one of
width w_j, then G has an ordering of width at most w + max_j w_j + 3
(the max of the empty set is 0).

*Proof.* (a) An edge with both endpoints in V(B) lies in B (otherwise
it is a bridge of G whose endpoints are joined by a path in B, closing
a cycle). An edge with exactly one endpoint in V(B) is a bridge of G.
If two distinct edges xy, x'y' joined F_j to B (x, x' in F_j), a path
x' -> x inside F_j and a path y -> y' inside B would close a cycle
through the bridge xy; so exactly one edge joins F_j to B, and no edge
joins F_i to F_j. Any edge f of B lies on a cycle of G, which uses no
bridge, hence lies in B; so f is not a bridge of B, and B is
bridgeless. Counting: each F_j is connected with |E(F_j)| =
|V(F_j)| - 1 (a tree, being joined to the rest by a single edge), so
M - N + 1 = (M_B - N_B + 1) + sum_j (M_j - N_j) + m = mu(B) +
sum_j mu(F_j). (b) Walk B's ordering; on reaching v, first list each
F_j with v_j = v (by its own ordering), then v itself. Let Q be a
prefix not containing all of B, u the first vertex of B not in Q. At
most w edges of B cross the cut. An F_j with v_j before u lies inside
Q, one with v_j after u lies outside. Of the F_j with v_j = u (listed
consecutively) at most one is split, in a prefix of its ordering: at
most max_j w_j of its edges cross. A bridge e_j crosses only if v_j =
u, and at most 3 edges at u are bridges (deg u <= 3). If Q contains
all of B the cut is empty. qed

**Lemma CENTROID.** Every connected subcubic graph with N vertices has
a *centroid* block: a block B such that every component of G - V(B)
has at most N/2 vertices.

*Proof.* Start at any block. If the current block B has a component F
of G - V(B) with more than N/2 vertices (at most one can), let xy be
the edge joining F to B with x in F, and move to the block B'
containing x. B' is bridgeless and connected, so B' lies in V(F). Every
path from B into F passes through x, so the component of G - V(B')
containing B is V(G) \ V(F), of fewer than N/2 vertices; hence a large
component of G - V(B'), if any, lies inside F \ V(B') and is smaller
than F. The walk stops, at a centroid. qed

**Theorem K2C (kernels to circuits).** (uses Lemma HARD) Let kappa > 0
and suppose that for every eta > 0 there is C_eta >= 0 such that every
kernel H on h vertices has cw(H) <= (kappa + eta) h + C_eta. Then (a)
for every eta > 0 there is D_eta such that every connected subcubic G
with N vertices has an ordering of width at most (2 kappa + 2 eta)
mu(G) + D_eta ceil(log2(N+1))^2; (b) liminf_{n -> inf} s_n/n >= 1 +
1/(2 kappa).

*Proof.* (a) Fix eta, put c := 2 kappa + 2 eta and K_N := 6 ceil(log2
N) + C_eta + 3. (Blocks) Let B be connected, bridgeless, subcubic with
N_B <= N. By Lemma COMPRESS either B has an ordering of width <=
6 ceil(log2 N_B), or its kernel H has h <= 2(mu(B) - 1) <= 2 mu(B) and
an ordering of width <= (kappa + eta) h + C_eta <= 2(kappa + eta) mu(B)
+ C_eta = c mu(B) + C_eta (as h <= 2 mu(B) and -2(kappa+eta) <= 0).
Either way B has an ordering of width <= c mu(B) + K_N - 3 (mu(B) >= 0).
(Recursion) Claim: every connected subcubic G' with N' <= N vertices has
an ordering of width at most c mu(G') + (floor(log2 N') + 1) K_N;
induction on N'. Take a centroid block B (Lemma CENTROID) and the
components F_j of G' - V(B). B is bridgeless (Lemma ASSEMBLY(a)), so it
has an ordering of width <= c mu(B) + K_N - 3. Each F_j has <= N'/2
vertices, so floor(log2 |F_j|) <= floor(log2 N') - 1 and by induction
F_j has an ordering of width <= c mu(F_j) + floor(log2 N') K_N. By
Lemma ASSEMBLY(b) and mu(G') = mu(B) + sum_j mu(F_j) >= mu(B) +
max_j mu(F_j), the graph G' has an ordering of width <= c mu(B) + K_N -
3 + c max_j mu(F_j) + floor(log2 N') K_N + 3 <= c mu(G') +
(floor(log2 N') + 1) K_N (m = 0 gives c mu(B) + K_N - 3; this includes
N' = 1). Finally L_N := ceil(log2(N+1)) = floor(log2 N) + 1 >=
ceil(log2 N), so the additive term is at most L_N (6 L_N + C_eta + 3)
<= D_eta L_N^2 with D_eta := C_eta + 9 (as 6 + 3/L_N <= 9).

(b) Fix eta, c := 2 kappa + 2 eta. Let a circuit with s = s_N gates
compute f_N (N the input length; write m, t, lambda, k_0 := m + t +
lambda for the family parameters at N). By Lemma HARD(c) it reads
b >= mt = N - k_0 variables (a circuit's function depends only on its
read variables, and f_N depends on each z-block variable), and it has
at least one gate. Its wire graph is connected, subcubic, with
N_g = 2 t_g + 1 <= 2s + 1 vertices and mu = t_g - b + 1 <= s - mt + 1 =
s - N + k_0 + 1 by (1.2) (t_g <= s). Put l := log2(s + N + 2) >= 1. By
(a) it has an ordering of width w <= c (s - N + k_0 + 1) +
D_eta ceil(log2(2s + 2))^2 <= c(s - N + k_0 + 1) + 9 D_eta l^2
(ceil(log2(2s+2)) <= l + 2 <= 3l). By Lemma SIM, f_N has a read-once
program with at most 4(2s+1)2^w + 2N <= 2^{w+3+l} edges: indeed
4(2s+1) + 2N = 8s + 4 + 2N <= 8(s + N + 2) (as 4 + 2N <= 8 + 8N),
and with l := log2(s + N + 2) taken as a real number,
8(s + N + 2) = 2^{l+3} exactly, so 4(2s+1)2^w + 2N <=
2^w·(4(2s+1) + 2N) <= 2^w·2^{l+3} = 2^{w+3+l}. Lemma COUNT applies for all large N: M_N >= 2^{N-2m-t} >= K,
because the hypothesis of Lemma HARD gives
N - 2m - t >= m + t + lambda + t ceil(log2 t) + 2 > log2 K
(recall K = 2^{m+t+lambda+1}·t! <= 2^{m+t+lambda+1+t ceil(log2 t)}).
With (1.1), for all large N, the program's edge count 2^{N - L_2(N)}
is below the bound 2^{w+3+l}, so N - L_2(N) < w + 3 + l <= c(s - N +
k_0 + 1) + 9 D_eta l^2 + 4 l^2 (3 + l <= 4 l^2 for l >= 1):

    N - L_2(N) < c (s - N + k_0 + 1) + (9 D_eta + 4) l^2.

Now l^2 <= 2 log2^2(s+2) + 2 log2^2(N+2), and for any eta > 0 there is
C'_eta with 2(9 D_eta + 4) log2^2(s+2) <= eta s + C'_eta for all s >= 0
(a fixed quadratic of log2 is eventually below any positive linear
function; C'_eta can be taken as the maximum of the difference over the
finitely many small s). Hence

    (c + eta) s_N > (1 + c) N - c (k_0 + 1) - L_2(N)
                    - 2(9 D_eta + 4) log2^2(N+2) - C'_eta,

and since k_0 + L_2(N) = 5m + 4t + 3 lambda + 2 t ceil(log2 t) + 2
= o(N) (t <= N^(1/3), m <= N/t, lambda < 2(t+1): each term is
O(N^(2/3)) or O(N^(1/3) log N)), dividing by N and letting N -> inf
gives liminf s_N/N >= (1 + c)/(c + eta) = (1 + 2 kappa + 2 eta)/(2
kappa + 3 eta). Let eta -> 0: liminf s_N/N >= 1 + 1/(2 kappa). qed


## 2. The matched-insertion bound and the rule bound [MSW S2.5-2.8, S3]

### 2.1 Gaussian sign patterns

For x in R^j and real theta let s_theta(x) in {0,1}^j have
(s_theta(x))_i = 1 exactly when x_i <= theta. F : {0,1}^j -> R is
*symmetric* if F(1 - s) = F(s) and *submodular* if every mixed
difference Delta_ij F(s) = F(s + e_i + e_j) - F(s + e_i) - F(s + e_j) +
F(s) (i != j, s_i = s_j = 0) is <= 0. A *correlation matrix* is a
positive semidefinite (PSD) matrix with unit diagonal (it may be
singular). For a j x j correlation matrix C put Psi_theta(F; C) :=
E F(s_theta(X)), X ~ N(0, C).

**Lemma G (sign patterns).** Let F : {0,1}^j -> R be submodular.
(a) If F is also symmetric, then Psi_theta(F; C) <= Psi_0(F; C) for
every theta. (b) If C_0, C_1 are correlation matrices with (C_1)_ij >=
(C_0)_ij for all i != j, then Psi_theta(F; C_1) <= Psi_theta(F; C_0)
for every theta. (c) With r_ij := C_ij: for j = 2,
Pr[X_1, X_2 differ in sign] = arccos(r_12)/pi; for j = 3,
Pr[X_1, X_2, X_3 not all of one sign] = (1/(2 pi)) sum_{i<j}
arccos(r_ij).

*Proof.* (Reduction to positive definite C.) Each X_i has variance 1,
so a.s. no coordinate equals a given real number. Let W ~ N(0, I_j) be
independent and X^eps = sqrt(1 - eps) X + sqrt(eps) W; its correlation
matrix (1 - eps) C + eps I_j is positive definite, and the entrywise
order in (b) is preserved. As eps -> 0, X^eps -> X pointwise and
s_theta is locally constant at X a.s., so E F(s_theta(X^eps)) ->
E F(s_theta(X)) by bounded convergence. It suffices to prove (a), (b)
for positive definite matrices.

(a) Evenness: -X has the law of X and s_{-theta}(-X) = 1 - s_theta(X),
so symmetry gives Psi_{-theta} = Psi_theta. Monotonicity on [0, inf):
write F = sum_A c_A prod_{i in A} s_i (multilinear expansion, valid on
the cube), so Psi_theta = sum_A c_A Pr[X_i <= theta for all i in A].
Differentiating each orthant probability,

    d Psi_theta/d theta = phi(theta) sum_a E[ d_a F(s_theta(Y)) | X_a = theta ],

d_a F(t) = F(1_a, t) - F(0_a, t), Y := X without coordinate a, phi the
standard normal density. Fix theta >= 0 and a. Given X_a = theta, Y is
Gaussian with mean m = theta (C_aj)_{j != a}; every coordinate
satisfies m_j = theta C_aj <= theta, since |C_aj| <= 1 and theta >= 0
(if C_aj <= 0 then m_j <= 0 <= theta; if 0 < C_aj <= 1 then m_j =
theta C_aj <= theta). Take the reflection Y' := 2m - Y, which has the
same conditional law. For every j: Y_j > theta implies Y'_j = 2 m_j -
Y_j < 2 m_j - theta <= theta, so s_theta(Y') >= 1 - s_theta(Y)
coordinatewise. Submodularity says d_a F is nonincreasing in each
coordinate; symmetry gives d_a F(1 - t) = -d_a F(t). Hence,
conditionally on X_a = theta,
E d_a F(s_theta(Y)) = E d_a F(s_theta(Y')) <= E d_a F(1 - s_theta(Y))
= -E d_a F(s_theta(Y)),
so each summand is <= 0 and d Psi_theta/d theta <= 0. Psi_theta is
even and nonincreasing on [0, inf), hence largest at theta = 0.

(b) Let D = C_1 - C_0 (D_ii = 0, D_ij >= 0), C_u = (1 - u) C_0 + u C_1
and p_u the density of N(0, C_u). The Fourier transform of p_u is
exp(-xi^T C_u xi / 2), whose u-derivative is
(1/2) sum_{i != j} D_ij (i xi_i)(i xi_j) exp(...); consequently
d p_u/du = (1/2) sum_{i != j} D_ij d_i d_j p_u (differentiation under
the integral is justified by dominated convergence: the eigenvalues of
C_u lie in [lambda_*, j] for some lambda_* > 0, which bounds p_u and
its first two derivatives by one integrable Gaussian envelope). Hence
d/du Psi_theta(F; C_u) = (1/2) sum_{i != j} D_ij I_ij with
I_ij = Int F(s_theta(x)) d_i d_j p_u(x) dx. Write x = (x_i, x_j, x');
for fixed x', F(s_theta(x)) is constant on each quadrant of the
(x_i, x_j)-plane around the point (theta, theta); the integral of
d_i d_j p_u over a quadrant is +p_u(theta, theta, x') for the quadrants
{x_i, x_j <= theta} and {x_i, x_j > theta}, and -p_u(theta, theta, x')
for the other two. So I_ij = Int p_u(theta, theta, x') Delta_ij
F(s_theta(x')) dx' <= 0 by submodularity; Psi_theta(F; C_u) is
nonincreasing in u. (Explicitly, d_i d_j p_u integrated against a
quadrant-constant function reduces, by two integrations in x_i and
x_j, to the second mixed difference of the function at theta, i.e. to
Delta_ij F <= 0.)

(c) For j = 2 and r = r_12 < 1, write (X_1, X_2) = (U, r U +
sqrt(1 - r^2) V) with U, V independent N(0, 1). The event {X_1 > 0,
X_2 > 0} is a sector of opening pi - arccos r of the
rotation-invariant plane, so it has probability (pi - arccos r)/(2 pi);
by X -> -X, the signs differ with probability 1 - 2(pi - arccos
r)/(2 pi) = arccos(r)/pi (r = 1: X_1 = X_2 a.s.). For j = 3: three
nonzero numbers have 3 pairs of equal sign if all three signs agree and
exactly 1 pair of equal sign otherwise; taking expectations of
1[X_i, X_j same sign] = 1 - arccos(r_ij)/pi over the three pairs,
1 + 2q = sum_{i<j} (1 - arccos(r_ij)/pi), q = Pr[all signs agree],
so 1 - q = (1/(2 pi)) sum_{i<j} arccos(r_ij). qed

### 2.2 The field, the keys, and the matched-edge cost

Fix an integer R >= 1 (the *field radius*). In a graph H give each
vertex x the row a_x in R^{V(H)} with a_x(z) = 2^{-d(x,z)/2} if
d(x, z) <= R, and a_x(z) = 0 otherwise (d = graph distance); the
entries of a_x are nonnegative and a_x(x) = 1. Put a-hat_x :=
a_x/|a_x|. Let omega = (omega_z) be independent standard normal
variables (the *noise*). For an edge e = xy the *key* is the unit
vector y_e := (a-hat_x + a-hat_y)/|a-hat_x + a-hat_y| (the denominator
is sqrt(2 + 2 <a-hat_x, a-hat_y>) >= sqrt 2), and the *score* of e is
k_e := <y_e, omega>. The scores of finitely many edges form a centred
Gaussian vector with Cov(k_e, k_f) = <y_e, y_f>; in particular every
score is N(0, 1).

*Matched edges.* Let M be a matching of a simple cubic graph H and uv
in M. Label N(u) = {v, a, b} and N(v) = {u, c, d}. The *boundary
edges* of uv are ua, ub, vc, vd, always in this order; they are four
distinct edges and none lies in M. At a threshold theta the *boundary
sides* are s_uv(theta) := (1[k_e < theta])_{e = ua, ub, vc, vd}. If uv
is put on side x in {0, 1}, the number of split vertices among u, v is

    g(x; s) := 1[x, s_1, s_2 not all equal] + 1[x, s_3, s_4 not all equal].

The cost function used throughout, with coordinates (s_1, s_2, s_3,
s_4) = (ua, ub, vc, vd) abbreviated (a, b, c, d), is

    F(s) := 1[s_a != s_b] + 1[s_c != s_d] + 1[s_a = s_b != s_c = s_d],    (2.3)

and Phi_F(C) := Psi_0(F; C) = E F(s_0(X)), X ~ N(0, C) (a 4 x 4
correlation matrix).

**Lemma MI (the optimal side).** For every s in {0,1}^4:
g(1; s) - g(0; s) = n00(s) - n11(s), where n11 (n00) counts the pairs
{s_1, s_2}, {s_3, s_4} whose two entries are both 1 (both 0); so
x*(s) := 1[n11(s) > n00(s)] is a least minimiser of g(.; s). Moreover
min_x g(x; s) = F(s); F is symmetric and submodular, with the mixed
differences Delta_12 F = -1 - 1[s_3 != s_d], Delta_34 F = -1 -
1[s_1 != s_2] and Delta_ij F = -1[s_k != s_l] for the cross pairs
ij in {13, 14, 23, 24} (kl the complementary pair).

*Proof.* For one pair (s_1, s_2), comparing x = 1 with x = 0 in
1[x, s_1, s_2 not all equal]: if s_1 = s_2 = 1 the difference is -1;
if s_1 = s_2 = 0 it is +1; if s_1 != s_2 it is 0. Adding the two pairs
gives n00 - n11. min_x g(x; s) = F(s): if both pairs are split, g = 2
for both x, and F = 1 + 1 + 0 = 2; if exactly one pair is split, g = 1
for x equal to the unsplit pair's value and 2 otherwise, and F = 0 + 1
+ 0 = 1; if all four coordinates are equal, g = 0 and F = 0; if s_1 =
s_2 != s_3 = s_4, g = 1 for every x and F = 0 + 0 + 1 = 1. F(1 - s) =
F(s) because every term only compares coordinates; F is invariant under
s_1 <-> s_2, s_3 <-> s_4 and (s_3, s_4) <-> (s_1, s_2). Mixed
differences, by direct evaluation over the free coordinates (all
settings of each free pair are machine-verified in `checks.py`,
part 1). Pair (1,2), (s_3, s_4) free: Delta_12 F = -1 - 1[s_3 != s_4]
(if s_3 != s_4 the four values are F(11) = 2, F(10) = 2, F(01) = 2,
F(00) = 1, sum -1; if s_3 = s_4 they are F(11) = 1, F(10) = 1,
F(01) = 1, F(00) = 0, sum -2). Pair (3,4) is symmetric. Cross pair
(1,3), (s_2, s_4) free: Delta_13 F = -1[s_2 != s_4] (the four values
are 2, 1, 1, 0 when s_2 = s_4 = 0, sum 0; and 1, 2, 0, 1 when
s_2 = 0, s_4 = 1, sum -1); the other cross pairs by the symmetries.
All mixed differences are <= 0, so F is submodular. qed

### 2.3 The 123 types and the cone

Order the six labelled vertices of a matched edge as (u, v, a, b, c, d)
and its five *internal edges* uv, ua, ub, vc, vd. A vector t =
(t_u, ..., t_d) in Z^6 with t_u = 0 is *generic* if (E) |t_x - t_y| <=
1 for the five internal edges xy; (Ru) min(t_v, t_a, t_b) = -1; (Rv)
min(t_u, t_c, t_d) = t_v - 1. The *special types* are the 4 vectors
(0, 1, 1, 1, t_c, t_d) with t_c, t_d in {1, 2} and the 4 vectors
(0, -1, t_a, t_b, 0, 0) with t_a, t_b in {0, 1}. T123 is the set of
generic and special types.

**Lemma TYPES.** Let uv be an edge of a simple cubic graph H, labelled
as above, and z any vertex. Then t(z) := (d(x, z) - d(u, z))_{x =
u,v,a,b,c,d} lies in T123, and |T123| = 123. Every type has |t_x| <= 2.

*Proof.* t_u = 0, and |d(x,z) - d(y,z)| <= d(x,y) = 1 on each internal
edge: (E). So t_v, t_a, t_b in {-1, 0, 1} and t_c, t_d in t_v +
{-1, 0, 1}; hence |t_x| <= 2. If z not in {u, v}: the second vertex of
a shortest u-z path is a neighbour of u at distance d(u,z) - 1, so
min over neighbours of u of (d(x,z) - d(u,z)) = -1, i.e. (Ru);
likewise at v: min over x in {u, c, d} of (d(x,z) - d(u,z)) =
(d(v,z) - 1) - d(u,z) = t_v - 1, i.e. (Rv). So t(z) is generic. z = u:
t = (0, 1, 1, 1, d(c,u), d(d,u)) with d(c,u), d(d,u) in {1, 2} (c != u
and c-v-u is a path of length 2); z = v: t = (0, -1, d(a,v) - 1,
d(b,v) - 1, 0, 0) with d(a,v), d(b,v) in {1, 2}: both are special
types. Count of generic types by t_v: t_v = -1: (Ru) automatic, 9 pairs
(t_a, t_b) in {-1,0,1}^2; (Rv) asks min(-1+0... min(t_u, t_c, t_d) =
min(0, t_c, t_d) = -2 with t_c, t_d in {-2,-1,0}: 5 pairs; 45. t_v =
0: 5 pairs (t_a, t_b) with min(t_v, t_a, t_b) = min(0, t_a, t_b) = -1,
and 5 pairs (t_c, t_d) with min(0, t_c, t_d) = -1 with t_c, t_d in
{-1,0,1}: 25. t_v = 1: 5 pairs (t_a, t_b) with min = -1; (Rv)
min(0, t_c, t_d) = 0 with t_c, t_d in {0,1,2} always holds: 9 pairs;
45. Total 115 generic types. No special type is generic (the z = u
types have min(t_v, t_a, t_b) = 1 != -1, the z = v types have
min(t_v, t_a, t_b) = -1 but min(t_u, t_c, t_d) = 0 != t_v - 1 = -2),
and the 8 special types are pairwise distinct. So |T123| = 123. (The
count and the axioms are machine-verified in `checks.py`, part 1.) qed

*The cone.* For t in T123 the *generator* is the vector g_t :=
(2^{-t_x/2})_{x = u..d} in R^6, with positive entries and g_t(u) = 1.
For a mass vector m in [0, inf)^{T123} \ {0} put

    G(m) := sum_{t in T123} m_t g_t g_t^T.

For a 6 x 6 PSD matrix G with positive diagonal put G-hat := corr G
(G-hat_xy = G_xy / sqrt(G_xx G_yy)), let B be the 4 x 6 incidence
matrix of the four boundary edges (row ua has ones in columns u and a,
etc.), W := B G-hat B^T, and C(G) := corr W. Put C(m) := C(G(m)) and

    f(m) := Phi_F(C(m)),      M_123 := sup { f(m) : m >= 0, m != 0 }.

f is invariant under scaling m (all of G, G-hat, W, C scale
homogeneously), so the supremum may be taken over the compact simplex
sum_t m_t = 1, where it is attained: G(m)_uu = sum m_t > 0; every
G_xx lies in [G_uu/4, 4 G_uu] because |t_x| <= 2 and g_t(u) = 1; every
W_ee >= 2 because G-hat >= 0 entrywise; so C(m) is well defined and
continuous on the simplex, and Phi_F is continuous on the compact set
of 4 x 4 correlation matrices (Lemma A.4, Section 2.5). Example: a
single generator (m concentrated at one type) gives G-hat with all
entries 1 on the support, and if all mass is on one type, G = m_t g g^T
has G-hat the all-ones matrix, W all entries 4, C the all-ones matrix,
and F(all ones) = 0: mixing types is what produces cost.

*Rows and keys.* Fix m. In R^{T123} with the standard inner product the
*unnormalised rows* are r_x := (sqrt(m_t) g_t(x))_t for x in {u, ...,
d}; their Gram matrix is G(m). Put a-hat_x := r_x/|r_x|, so that
G-hat_xy = <a-hat_x, a-hat_y> >= 0. For an internal edge e = xy put
y_e := (a-hat_x + a-hat_y)/|a-hat_x + a-hat_y| and, with a standard
Gaussian vector omega in R^{T123}, Z_e := <y_e, omega>. The *edge
cosines* are e_xy := G-hat_xy for the five internal edges.

**Lemma REAL.** (Z_ua, Z_ub, Z_vc, Z_vd) is a centred Gaussian vector
with covariance C(m); Z_uv is jointly Gaussian with it; every Z_e has
variance 1; and f(m) = E F(s_0(Z_ua, Z_ub, Z_vc, Z_vd)).

*Proof.* Each y_e is a unit vector, so Z_e ~ N(0, 1) and Cov(Z_e, Z_f)
= <y_e, y_f>. For boundary edges e = xy and f' = x'y',
<y_e, y_f'> = <a-hat_x + a-hat_y, a-hat_x' + a-hat_y'> /
(|a-hat_x + a-hat_y| |a-hat_x' + a-hat_y'|) = W_{ef'} /
sqrt(W_ee W_{f'f'}) = C(m)_{ef'}. Hence the four boundary scores are
centred Gaussian with covariance C(m) = corr of the key Gram block, and
E F(s_0 of the boundary keys) = Psi_0(F; C(m)) = f(m). qed

### 2.4 Edge cosines, the star lemma, and the rule bound

**Lemma EC (edge cosines in ratio form).** For each internal edge xy,
with r := G_yy/G_xx: e_xy >= (rho/2)(sqrt r + 1/sqrt r) >= rho, where
rho := 2 sqrt 2/3.

*Proof.* For each type, |t_x - t_y| <= 1, so g_t(y)/g_t(x) lies in
[q, Q] with q = 2^{-1/2}, Q = 2^{1/2}; qQ = 1 and q + Q = 3/sqrt 2.
For alpha, beta > 0 with beta/alpha in [q, Q]: (beta - q alpha)(Q
alpha - beta) >= 0, i.e. beta^2 + alpha^2 <= (3/sqrt 2) alpha beta.
With alpha = g_t(x), beta = g_t(y), weight m_t and summing over t:
G_yy + G_xx <= (3/sqrt 2) G_xy (all terms nonnegative). Divide by
sqrt(G_xx G_yy): sqrt r + 1/sqrt r <= (3/sqrt 2) e_xy. Since
sqrt 2/3 = rho/2, this is e_xy >= (rho/2)(sqrt r + 1/sqrt r) >=
(rho/2) * 2 = rho by AM-GM. qed

For e in [-1, 1]^3 put mu(e_1, e_2, e_3) := ((sum_i sqrt((1+e_i)/2))^2
- 3)/6; it is continuous and increasing in each argument on [-1,1]^3,
and mu(rho, rho, rho) = (9(1+rho)/2 - 3)/6 = (1 + 3 rho)/4 =: mu0 =
(1 + 2 sqrt 2)/4. The *key correlations* at u are k_a := <y_ua, y_uv>,
k_b := <y_ub, y_uv> and k_ab := <y_ua, y_ub> = C_ab; at v, k_c, k_d
and k_cd := C_cd likewise. Their mean at w in {u, v} is kbar_w, and
S_w is the sum of the three arccos values at w.

**Lemma STAR (the star lemma).** Let w be u or v, with edge cosines
e_1, e_2, e_3 (at u: e_uv, e_ua, e_ub; at v: e_vu, e_vc, e_vd) and key
correlations k_1, k_2, k_3 (at u: k_a, k_b, k_ab; at v: k_c, k_d,
k_cd). Then (a) each k_i lies in [rho, 1]; (b) kbar_w >= mu(e_1, e_2,
e_3) >= mu0; (c) S_w <= 3 arccos(kbar_w); (d) S_w <= 3 arccos(mu0) =
2 pi kappa_E.

*Proof.* Let y_1, y_2, y_3 be the keys of the three edges at w and
c_i := <y_i, a-hat_w> = (1 + e_i)/sqrt(2 + 2 e_i) = sqrt((1+e_i)/2) >=
kappa_0 := sqrt((1+rho)/2) (Lemma EC gives e_i >= rho). (a) Write y_i =
c_i a-hat_w + w_i with w_i orthogonal to a-hat_w, and psi_i := arccos
c_i, which lies in [0, arccos kappa_0] within [0, pi/2] (c_i >= kappa_0
> 0). Then <y_i, y_j> = c_i c_j + <w_i, w_j> >= c_i c_j - |w_i||w_j| =
cos psi_i cos psi_j - sin psi_i sin psi_j = cos(psi_i + psi_j) >=
cos(2 arccos kappa_0) = 2 kappa_0^2 - 1 = 1 + rho - 1 = rho (using
sin psi_i <= sin(arccos kappa_0) etc.: |w_i| = sin psi_i). And
<y_i, y_j> <= 1 by Cauchy-Schwarz. (b) |y_1 + y_2 + y_3| >= <y_1 +
y_2 + y_3, a-hat_w> = c_1 + c_2 + c_3 >= 0, and |y_1+y_2+y_3|^2 =
3 + 2 sum_{i<j} k_ij = 3 + 6 kbar_w, so kbar_w >= ((c_1+c_2+c_3)^2 -
3)/6 = mu(e_1, e_2, e_3) >= mu(rho, rho, rho) = mu0, by monotonicity
of mu in each argument and e_i >= rho. (c) arccos is concave on [0, 1]
(second derivative -x(1 - x^2)^{-3/2} <= 0) and the k_i lie in [rho,
1] within [0, 1]; Jensen: (arccos k_1 + arccos k_2 + arccos k_3)/3 <=
arccos((k_1+k_2+k_3)/3) = arccos(kbar_w). (d) arccos is decreasing and
kbar_w >= mu0, so by (c)
S_w <= 3 arccos(kbar_w) <= 3 arccos(mu0). 3 arccos(mu0) = 2 pi
kappa_E is the definition
of kappa_E with mu0 = (1 + 2 sqrt 2)/4, certified in `checks.py` part
2. qed

**Lemma RB (the rule bound).** For every m: f(m) <= RB(m) := (S_u +
S_v)/(2 pi) <= 2 kappa_E. Hence M_123 <= 2 kappa_E.

*Proof.* By Lemma MI, F(s) = min_x g(x; s) <= g(X; s) for every X. Take
X := 1[Z_uv <= 0] and s := s_0(Z): pointwise,
F(s_0(Z)) <= 1[Z_uv, Z_ua, Z_ub not all of one sign] + 1[Z_uv, Z_vc,
Z_vd not all of one sign]. (Z_uv, Z_ua, Z_ub) is centred Gaussian with
unit variances and correlations <Z_uv, Z_ua> = k_a, <Z_uv, Z_ub> =
k_b, <Z_ua, Z_ub> = k_ab (Lemma REAL applied to the five internal-edge
keys), so by Lemma G(c) the first term has expectation (arccos k_a +
arccos k_b + arccos k_ab)/(2 pi) = S_u/(2 pi); likewise the second is
S_v/(2 pi). Take expectations and use Lemma STAR(d):
f(m) = E F(s_0(Z)) <= (S_u + S_v)/(2 pi) <= 6 arccos(mu0)/(2 pi) =
(3/pi) arccos(mu0) = 2 kappa_E. Taking the supremum over m,
M_123 <= 2 kappa_E. qed


### 2.5 Matched insertion [MSW S3]

**Lemma A.2 (the cone with a shell).** Fix a simple cubic graph H, the
field of Section 2.2 with radius R >= 3, and an edge uv labelled as
above. Let G be the 6 x 6 Gram matrix G_xy = <a_x, a_y>, x, y in {u,
v, a, b, c, d}. Then G = G_in + G_sh, where

    G_in = sum_{t in T123} m_t g_t g_t^T,   m_t = sum { 2^{-d(u,z)} : t(z) = t, d(u,z) <= R - 2 },

with m_{t(u)} >= 1, and G_sh is PSD with tr G_sh <= (60/(R - 2)) tr
G_in.

*Proof.* G = sum_z r_z r_z^T with r_z = (a_x(z))_x. If d(u,z) <= R -
2, then d(x,z) <= d(x,u) + d(u,z) <= 2 + (R - 2) = R for each of the
six x (each x is within distance 2 of u), so a_x(z) = 2^{-d(x,z)/2} =
2^{-d(u,z)/2} 2^{-t_x(z)/2}, i.e. r_z = 2^{-d(u,z)/2} g_{t(z)} with
t(z) in T123 (Lemma TYPES). Summing over these z gives G_in, and z = u
(t(u) = (0,1,1,1,d(c,u),d(d,u)), a special type, contributing
2^{-d(u,u)} = 1) gives m_{t(u)} >= 1. G_sh := sum over d(u,z) >= R -
1 of r_z r_z^T is PSD; if d(u,z) >= R + 3 then all six d(x,z) >=
(R + 3) - 2 = R + 1 and r_z = 0, so only R - 1 <= d(u,z) <= R + 2
contributes. Let N_k := #{z : d(u,z) = k} and s_k := N_k 2^{-k}. For
k >= 1 every vertex at distance k + 1 has a neighbour at distance k,
and every vertex at distance k has a neighbour at distance k - 1, hence
at most 2 at distance k + 1: N_{k+1} <= 2 N_k and s_{k+1} <= s_k. For
z at distance k: d(x,z) >= k - d(x,u) with d(x,u) in {0,1,1,1,2,2}, so
|r_z|^2 <= sum_x 2^{-d(x,z)} <= 2^{-k}(1 + 3*2 + 2*4) = 15 * 2^{-k}.
Hence tr G_sh = sum_{R-1 <= d(u,z) <= R+2} |r_z|^2 <= 15 sum_{k =
R-1}^{R+2} s_k <= 60 s_{R-2} (s_k nonincreasing). And tr G_in = sum_x
G_in,xx >= sum_{d(u,z) <= R-2} a_u(z)^2 = sum_{k=0}^{R-2} s_k >=
(R - 2) s_{R-2}. Combining, tr G_sh <= 60 s_{R-2} <= (60/(R - 2)) tr
G_in. qed

**Lemma A.3 (the shell moves C by O(1/R)).** Let G_in = G(m) be a
nonzero cone point and S a PSD 6 x 6 matrix with tr S <= eps tr G_in,
where delta := 96 eps <= 1/12. Then G := G_in + S has positive
diagonal, W(G) has positive diagonal, and |C(G)_ef - C(G_in)_ef| <= 24
delta for all e, f.

*Proof.* g_t(x)^2 = 2^{-t_x} lies in [1/4, 4] and g_t(u) = 1, so for
every type g_t(x)^2 <= 4 g_t(u)^2; summing over t (m >= 0):
G_in,xx <= 4 G_in,uu and G_in,xx >= G_in,uu/4 for every x. Thus tr
G_in <= 24 G_in,uu and G_in,xx >= tr G_in/96 > 0 for every x. By
Cauchy-Schwarz and PSDity, |S_xy| <= sqrt(S_xx S_yy) <= tr S <=
delta * min(G_in,xx, G_in,yy), and 0 <= S_xx <= delta G_in,xx; so all
diagonals of G are positive and at least (1 - delta) G_in,xx.
(Normalized entries) G-hat(G)_xy = G_xy/sqrt(G_xx G_yy) = (G-hat(G_in)_xy
+ e1) phi where |e1| <= delta (from the entrywise bounds) and phi =
((1 + e2)(1 + e3))^{-1/2} with e2, e3 in [0, delta]; so phi lies in
[1/(1 + delta), 1] and |G-hat(G)_xy - G-hat(G_in)_xy| <= 2 delta.
(W entries) W(G_in)_ee = 2 + 2 G-hat(G_in)_xy >= 2 (G-hat >= 0
entrywise); each entry of W is a sum of 4 entries of G-hat whose
diagonal entries are exactly 1, so |W(G)_ef - W(G_in)_ef| <= 8 delta,
and W_ee <= 8 gives 8 delta <= 6 delta sqrt(W_ee W_ff) whenever
W_ee, W_ff in [4·... one checks 8 delta <= 12 delta <= 6 delta
sqrt(W_ee W_ff) as soon as sqrt(W_ee W_ff) >= 2, which holds since
W_ee, W_ff >= 2. (C entries) C(G)_ef = W(G)_ef / sqrt(W(G)_ee W(G)_ff)
= (C(G_in)_ef + e1') phi' with |e1'| <= 6 delta <= 1/2, e2', e3' in
[-6 delta, 6 delta] with |e2'|, |e3'| <= 1/2, so phi' lies in
[1 - 6 delta, 1 + 12 delta] (as sqrt(1 +- x) >= 1 +- x for |x| <=
1/2), and |C(G)_ef - C(G_in)_ef| <= 6 delta(1 + 12 delta) + 12 delta
<= 24 delta (using 72 delta^2 <= 6 delta for delta <= 1/12). qed

**Lemma A.4 (Phi_F is uniformly continuous).** On the compact set CC of
4 x 4 correlation matrices, Phi_F is continuous, hence uniformly
continuous.

*Proof.* (uses the classical Levy continuity theorem and portmanteau
theorem) If C_n -> C in CC, the characteristic functions
exp(-xi^T C_n xi / 2) converge pointwise to exp(-xi^T C xi / 2), so
N(0, C_n) converges to N(0, C) in distribution (Levy's continuity
theorem). The map x -> F(s_0(x)) is bounded and continuous off the
hyperplanes {x_i = 0}, which have probability 0 under N(0, C) (all
variances are 1); by the portmanteau theorem, E F(s_0(X_n)) ->
E F(s_0(X)). CC is compact (closed and bounded: PSD with unit diagonal
has entries in [-1, 1]), and a continuous function on a compact set is
uniformly continuous. qed

For R >= 3 put eta_R := sup { |Phi_F(C) - Phi_F(C')| : C, C' in CC,
max |C_ef - C'_ef| <= 24 * 96 * 60/(R - 2) }. By Lemmas A.3 and A.4,
eta_R -> 0 as R -> inf. No explicit modulus is needed: R is chosen
from the target slack alone.

**Lemma A.6 (the mixed order; deterministic).** Fix H, a perfect
matching M, and real scores k_e for all edges (ties allowed). The
*entry time* of an edge is tau_e := k_e if e not in M, and for uv in
M, tau_uv := inf { theta : x*(s_uv(theta)) = 1 }. List the edges as
e_1, ..., e_{3h/2} by increasing (tau_e, a fixed index of e). For a
vertex w, first(w) and last(w) are the least and largest positions of
its three edges; the *bags* are X_i := {w : first(w) <= i <=
last(w)}. The *cost* at theta is Cost(theta) := #{w : w has an edge
with tau < theta and an edge with tau >= theta}. For a grid theta_0 <
... < theta_J the *window counts* are W_j := #{e not in M : theta_j
<= k_e < theta_{j+1}} and the *tail counts* W_- := #{e not in M : k_e
< theta_0}, W_+ := #{e not in M : k_e >= theta_J}. Then:

(a) For uv in M, tau_uv is one of its four boundary scores, and
    1[tau_uv < theta] = x*(s_uv(theta)) for every theta.
(b) For every theta, Cost(theta) = sum_{uv in M} F(s_uv(theta)).
(c) The bags X_i form a path decomposition of H (every vertex lies in
    X_{first(w)}, every edge e_i = xy has both endpoints in X_i, and
    the bags containing any fixed w are consecutive).
(d) For every grid, every bag satisfies |X_i| <= 6 W_-, or |X_i| <=
    6 W_+, or |X_i| <= Cost(theta_j) + 6 W_j for some j.
(e) Ordering the vertices by the position of their median edge (ties
    by a fixed numbering) gives every cut at most max_i |X_i| + 1
    edges.

*Proof.* (a) Each coordinate of s_uv(theta) is nondecreasing in theta
and changes only at a boundary score, so n11(s_uv(theta)) -
n00(s_uv(theta)), and with it x* = 1[n11 > n00] (Lemma MI), is a
nondecreasing step function of theta, constant between consecutive
boundary scores (with the convention 1[k_e < theta]). For theta below
the least boundary score, s = 0000 and x* = 0; above the largest, s =
1111 and x* = 1. So {theta : x*(s_uv(theta)) = 1} is a half-line (k_l,
inf) for one boundary score k_l, tau_uv = k_l, and
1[tau_uv < theta] = x*(s_uv(theta)) for every theta. (b) A vertex is
counted in Cost(theta) iff its three edges are not all on one side of
theta. Every vertex is covered by M. For uv in M the side of uv is
1[tau_uv < theta] = x*(s), and the sides of the boundary edges are s =
s_uv(theta); the number of counted vertices among u, v is g(x*(s); s)
= F(s) (Lemma MI). (c) As stated. (d) Fix i and w in X_i, t_i :=
tau_{e_i}. Then w has an edge with tau <= t_i and one with tau >= t_i.
If the entry times of the edges in a set I of reals lie in an interval:
each edge f with tau_f in I is either unmatched with k_f in I, or
matched with tau_f = k_e in I for a boundary edge e of f; an unmatched
edge is a boundary edge of at most two matched edges; so #{f : tau_f
in I} <= 3 #{e unmatched : k_e in I}, and these edges have at most 6
#{...} endpoints. If t_i < theta_0: w is an endpoint of an edge with
tau < theta_0, so |X_i| <= 6 W_-. If t_i >= theta_J: likewise |X_i|
<= 6 W_+. If theta_j <= t_i < theta_{j+1}: either w has an edge with
tau < theta_j - then w is counted in Cost(theta_j) - or its edge with
tau <= t_i has tau in [theta_j, theta_{j+1}), and at most 6 W_j
vertices are of this kind. (e) Fix a nonempty proper prefix P of the
vertex ordering, its last vertex x_0 and theta := the median-edge
position of x_0. The three edges at a vertex have distinct positions
(H is loopless), so at least one is below and one above the median. A
crossing edge e = xy (x in P) with position < theta is the lowest edge
of y, and y lies in the bag at theta (first(y) < theta <= med(y) <
last(y)): charge it to y. With position > theta it is the highest edge
of x and x lies in the bag at theta (first(x) < med(x) = theta <
last(x)): charge it to x. At most one edge has position exactly theta.
No vertex is charged twice (a vertex is charged at most one edge of
each kind), so the cut has at most |X_theta| + 1 edges. qed

**Lemma A.7 (the expected cost at every threshold).** Let
R >= 2 + 12 * 96 * 60 = 69122. For every simple cubic graph H, every
perfect matching M and every real theta:

    E Cost(theta) <= (h/2)(M_123 + eta_R).

*Proof.* For uv in M, Z := (k_ua, k_ub, k_vc, k_vd) is a centred
Gaussian vector with covariance C(G), G the Gram matrix of Lemma A.2:
for boundary edges e = xy, f' = x'y', Cov(k_e, k_{f'}) = <y_e, y_f'>
= W(G)_{ef'}/sqrt(W(G)_ee W(G)_{f'f'}) = C(G)_{ef'} (Lemma REAL's
computation, now for the field Gram matrix), and each score is N(0, 1).
s_uv(theta) and s_theta(Z) differ only on a null event (each score has
a continuous law), so E F(s_uv(theta)) = Psi_theta(F; C(G)) <=
Phi_F(C(G)) by Lemma G(a) (F symmetric and submodular, Lemma MI). By
Lemma A.2, G = G(m) + G_sh with tr G_sh <= eps tr G(m), eps =
60/(R - 2) (so delta = 96 eps <= 96 * 60/69120 < 1/12); by Lemma A.3,
max |C(G)_ef - C(m)_ef| <= 24 * 96 * 60/(R - 2), so Phi_F(C(G)) <=
f(m) + eta_R <= M_123 + eta_R (by the definition of eta_R and
Lemma A.4). By Lemma A.6(b), E Cost(theta) = sum_{uv in M} E
F(s_uv(theta)), a sum of h/2 terms. qed

**Fact PETERSEN (cited).** [Diestel, Graph Theory, Corollary 2.2.2.]
Every bridgeless cubic (multiigraph-free, i.e. simple or with parallel
edges) graph has a perfect matching. (The classical theorem covers
bridgeless cubic multigraphs, hence simple cubic graphs.)

**Theorem A (matched insertion, perfect matchings).** For every xi in
(0, 1] there is h_0(xi) such that every simple cubic graph H on h >=
h_0(xi) vertices with a perfect matching has a path decomposition
whose bags all have size at most (M_123/2 + xi) h, and a vertex
ordering with every cut of size at most (M_123/2 + xi) h + 1.

*Proof.* (Constants, depending on xi only.) Choose an integer R >=
69122 with eta_R <= xi/8 (possible by Lemma A.4). Delta := xi/100;
choose T >= 1 with (3/2)(1 - Phi(T)) <= xi/100 (Phi the standard
normal distribution function). J := ceil(2T/Delta), theta_j := -T +
j Delta (j = 0, ..., J); there are L := 2J + 3 counts in all:
Cost(theta_j) (J + 1 of them), W_j (J), W_-, W_+. K := 3 * 2^{2R+6},
lambda(h) := xi h / 40, h_0 := ceil(1600 L K / xi^2) + 1. Fix H with h
>= h_0 and a perfect matching M; build the field with this R; draw
omega. (Means) E Cost(theta_j) <= (h/2)(M_123 + xi/8) = (M_123/2) h +
xi h/8 (Lemma A.7). Each unmatched score is N(0, 1), so the expected
number of them with score in a window of width Delta is at most Delta
times the maximum density 1/sqrt(2 pi) < 0.4, times the 3h/2 edges:
E W_j <= (3h/2) Delta (1/sqrt(2 pi)) < 1.5 * (xi/100) * 0.4 * h <
0.006 xi h. E W_-, E W_+ <= (3h/2)(1 - Phi(T)) <= xi h / 100. (Locality) The score of an edge xy
depends only on omega on the ball B(x, R + 1) (a_x uses vertices within
R of x, y_e uses rows of x, y). The term F(s_uv(theta)) of uv depends
only on omega on B(u, R + 2) (the rows of u, v, a, b, c, d are
supported within distance R of their vertices, all within distance 2
of u). (Variances) Each count X has Var X <= K h: Cost(theta) is a sum
of h/2 terms with values in [0, 2] (variance <= 1), centred at the
first endpoint of each matched edge; two terms are independent when
their centres are at distance > 2R + 4 (their supporting balls are
disjoint); a ball of radius 2R + 4 in a cubic graph has fewer than
3 * 2^{2R+4} vertices (each sphere at most doubles the previous), so
the dependency neighbourhood of a term has at most 3 * 2^{2R+4} other
centres and Var Cost(theta) <= h * 1 * 3 * 2^{2R+4} <= K h. For W_j
assign each unmatched edge to one endpoint: W_j = sum_x Y_x with Y_x
in [0, 3] (variance <= 9/4), dependency radius 2R + 2, so Var W_j <=
h * (9/4) * 3 * 2^{2R+2} <= K h; likewise W_-, W_+. (One good omega)
By Chebyshev, each count deviates from its mean by >= lambda(h) with
probability <= K h / lambda(h)^2 = 1600 K / (xi^2 h); by the union
bound over the L counts, the probability that some count deviates is
<= 1600 L K / (xi^2 h) < 1 for h >= h_0. Fix an omega with no
deviation. (Bags) By Lemma A.6(d) each bag satisfies |X_i| <= 6 W_-
< 6(xi/100 + xi/40) h = 0.21 xi h (or the same with W_+), or
|X_i| <= Cost(theta_j) + 6 W_j < (M_123/2) h + (xi/8 + xi/40) h +
6(0.006 xi + 0.025 xi) h = (M_123/2) h + xi h (0.125 + 0.025 + 0.186)
= (M_123/2 + 0.336 xi) h. In all cases |X_i| <= (M_123/2 + xi) h.
Lemma A.6(c) and (e) finish: the bags form a path decomposition, and
the median-edge ordering has every cut at most max |X_i| + 1. qed

**Corollary A1.** (i) For every xi in (0, 1] there is C_xi such that
every kernel H on h vertices has cw(H) <= (M_123/2 + xi) h + C_xi.
(ii) M_123 > 0. (iii) Every family of B_2-circuits computing the
norm-one family f_N of Section 0 satisfies
liminf_{N -> inf} s_N/N >= 1 + 1/M_123.

*Proof.* A kernel H is bridgeless: if uv were a bridge with components
A (containing u), B of H - uv, then u has two neighbours in A (H is
simple and cubic), so A \ {u} is nonempty; likewise B \ {v}; every
path from A \ {u} to B \ {v} uses uv, hence u, so deleting u
disconnects H - against 2-connectivity. By Fact PETERSEN, H has a
perfect matching. (i) For h >= h_0(xi), Theorem A gives cw(H) <=
(M_123/2 + xi) h + 1; for h < h_0 every cut has at most 3h/2 <=
3 h_0/2 edges, so C_xi := 3 h_0(xi)/2 + 1 works. (ii) Take m with
mass 1 on each of the generic types A = (0, 1, -1, 1, 2, 2) and B =
(0, 1, 1, -1, 2, 2). With masses (1, 1) the unnormalised rows are r_u
= (1, 1), r_a = (2^{1/2}, 2^{-1/2}), r_b = (2^{-1/2}, 2^{1/2}), all
nonzero. a-hat_u + a-hat_a and a-hat_u + a-hat_b are positive vectors
that are not parallel (their coordinates are (p, q) and (q, p) with p
= a-hat_u(1) + a-hat_a(1) > q = a-hat_u(1) + a-hat_a(2)... computing:
a-hat_u = (1,1)/sqrt 2, a-hat_a = (sqrt 2, 1/sqrt 2)/sqrt(5/2) =
(2/sqrt 5, 1/sqrt 5), a-hat_b = (1/sqrt 5, 2/sqrt 5); so a-hat_u +
a-hat_a = (1/sqrt 2 + 2/sqrt 5, 1/sqrt 2 + 1/sqrt 5) and a-hat_u +
a-hat_b = (1/sqrt 2 + 1/sqrt 5, 1/sqrt 2 + 2/sqrt 5); two positive
vectors of the form (p, q), (q, p) with p != q are not parallel), so
the keys y_ua, y_ub are distinct unit vectors, both with positive
coordinates, hence <y_ua, y_ub> < 1, i.e. C(m)_ab < 1 (equality would
force y_ua = y_ub). Since F(s) >= 1[s_1 != s_2], Lemma G(c) gives
f(m) >= Pr[Z_1, Z_2 differ in sign] = arccos(C(m)_ab)/pi > 0. (The
parallelism check and arccos positivity are machine-verified in
`checks.py`, part 1.) (iii) For every eta > 0, (i) with xi :=
min(eta, 1) is the hypothesis of Theorem K2C with kappa := M_123/2 >
0, and Theorem K2C(b) gives liminf s_n/n >= 1 + 1/(2 kappa) = 1 +
1/M_123. qed

## 3. Assembly and the constant

By Lemma RB, M_123 <= 2 kappa_E; by Corollary A1(ii), M_123 > 0; the
map t -> 1 + 1/t is decreasing on (0, inf), so Corollary A1(iii)
gives

    liminf_{n -> inf} s_n/n >= 1 + 1/M_123 >= 1 + 1/(2 kappa_E),
    kappa_E = (3/(2 pi)) arccos(mu0),  mu0 = (1 + 2 sqrt 2)/4.

`checks.py` part 2 certifies, in directed-rounding interval arithmetic
(mpmath's iv type at 80 bits of precision, cross-validated against an
independent enclosure by exact rational interval bisection of arccos
via its integral representation):

    mu0 in [0.95710678118654752, 0.95710678118654760],
    arccos(mu0) in [0.2939503813270101, 0.2939503813270102],
    kappa_E in [0.14035096863582366, 0.14035096863582374],
    1 + 1/(2 kappa_E) = 4.5624976789250185... >= 4.56249767892.

Part 3 verifies, in exact rational arithmetic, that the declared reduced
fraction 114062441973/25000000000 equals 4.56249767892 (decimal) and
that it does not exceed the certified lower endpoint of 1 + 1/(2
kappa_E). Part 1 machine-checks the finite tables used above: the
mixed-difference table of F (all 4 choose 2 pairs x 16 settings), the
Lemma MI minimisation table, the count |T123| = 123 with the type
axioms (E), (Ru), (Rv) re-derived by enumeration, the special-type count
8 and generic count 115, and the non-parallelism of the two vectors in
Corollary A1(ii). The complete transcript is supplied as
`checks-output.txt`. This completes the proof of the Main Theorem. qed

*Relation to the source and to the board.* This file proves the
reduction (their Sections 2.1-2.4) and the matched-insertion framework
(their Sections 2.5-2.8 and 3) of [MSW] in full, and instantiates them
with the source's rule bound, which caps the cone supremum by
2 kappa_E; this is the constant their Section 1.2 attributes to the
edge-score ordering ("4.5624977 n"). The source's sharper Theorem
M_tree (M_123 = M_tree = 0.2765836..., coefficient 4.6155), which
requires their computer-assisted Sections 4-9, and the frontier
entry's further analytic improvement (coefficient 4.67235) are not used
or claimed. The hard family (Section 0) is the norm-one indicator
family; the same choice of hard function for this route appears in
ArthurOstapenko's accepted entry 808da37 (who also prove their counting
core self-containedly); our family, layout, rectangle argument and
KRS file are written independently, and `krs.md` proves the
counting theorem from first principles. Apart from the named classical
facts listed in the preamble (each used in one line), every lemma of
the chain is proved above in full.
