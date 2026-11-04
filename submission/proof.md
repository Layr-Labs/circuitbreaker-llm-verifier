# Coefficient 1 at exponent 1 for the OR family: complete self-contained proof

This proof is fully self-contained: it defines one fixed total family, proves
NP eligibility from first principles, states the pinned B2 model exactly as in
`rules.md`, and proves the matching lower and upper bounds with all steps
shown. It cites no external theorem for any scored conclusion. (Remarks on the
repository's prior baseline appear only in Section 7 and are not part of the
scored joint theorem.)

## 1. The fixed total family

Let

    L = { x in {0,1}* : x contains at least one bit equal to 1 }.

For each n >= 0 define f_n : {0,1}^n -> {0,1} by

    f_0(empty) = 0,
    f_n(x_1,...,x_n) = x_1 OR x_2 OR ... OR x_n   for n >= 1.

Thus f_n is the characteristic function of L at length n. Every length is
defined, including the small cases n = 0 and n = 1. This is one fixed total
Boolean family; no per-length advice, oracles, or choices are used.

## 2. Uniform polynomial-time evaluation and NP membership

One deterministic algorithm A computes every f_n: on input x of length n,
scan the n bits from left to right and output 1 if any bit is 1, otherwise
output 0. This takes O(n) time, polynomial in the actual input length n.

Define the verifier V(x, w) = 1 iff w is the empty string and A(x) = 1, and
take the witness bound W(n) = 0, a polynomial. We establish both directions
of NP membership for the same fixed V and W at every length:

- If x in L, then A(x) = 1. Take w = empty: |w| = 0 <= W(|x|) and
  V(x, w) = 1.
- If x is not in L, then A(x) = 0. Any nonempty w is rejected by the
  empty-witness check, and for w = empty we have V(x, w) = 0 because
  A(x) = 0. Hence no w with |w| <= W(|x|) satisfies V(x, w) = 1.

Therefore x in L iff there exists w with |w| <= W(|x|) and V(x, w) = 1, so
L is in NP (in fact in P, via A). The witness is never an input to the
circuits whose size we bound below.

## 3. The pinned circuit model

We work in exactly the model stated in `rules.md` (Li-Yang's general B2
model): deterministic acyclic Boolean circuits with one output (output
selection), two input pins per gate, each gate computing one of the 16 binary
Boolean truth tables, unrestricted depth and fan-out, exact computation on
every input assignment. Circuits may be nonuniform. The circuit inputs are
x_1,...,x_n alone, not NP witnesses. Cost: every retained gate counts once,
including trivial or degenerate gates; input nodes, wires, and output
selection cost zero. Repeated pins (both pins of a gate fed by the same node)
and shared subcircuits are allowed. Zero-cost constant nodes, if the
formalism admits them, are handled explicitly in Section 4. Our argument uses
only these stated conventions.

For each n, C_B2(f_n) denotes the minimum number of gates over all such
circuits computing f_n exactly on {0,1}^n.

## 4. Lower bound: C_B2(f_n) >= n - 1 for every n >= 1

Fix n >= 1 and let C be any B2 circuit computing f_n exactly. For n >= 2,
f_n is non-constant and every input coordinate is essential (shown below),
so the output node of C must be a gate: it can be neither a constant node
nor an input node, since otherwise some essential input would have no
directed path to it. (The case n = 1 needs no argument here: gate counts are
nonnegative, so C_B2(f_1) >= 0 = 1 - 1 holds trivially, and Section 5 shows
equality.) Trim C to the output gate and its ancestors: keep precisely those
nodes (gates, input nodes, and any constant nodes) that have a directed path
to the output gate, and discard everything else. Trimming cannot increase the
gate count, and the trimmed circuit still computes f_n.

Every input coordinate is essential for f_n: f_n(0,...,0) = 0 while
f_n(e_i) = 1, where e_i has a single 1 in position i. Hence f_n is not
constant in any coordinate. If some input node x_i had no directed path to
the output gate, then the output value would be unchanged when x_i varies
while the other inputs are fixed, contradicting essentiality. Indeed, order
the nodes topologically: a node not reachable from x_i has all its
pin-sources also not reachable from x_i, so by induction its value is
unaffected by varying x_i; the output gate is such a node. (Equivalently,
after trimming, x_i would not even be retained.) So each of the n input
nodes has a directed path to the output gate.

Consequently, in the underlying undirected graph of the trimmed circuit
(direction forgotten), every retained node is connected to the output node:
the graph is connected. Let g be the number of retained gates (each costing
one, including degenerate ones) and k the number of retained zero-cost
constant nodes (k = 0 if the formalism has none). The graph has
N = n + g + k distinct nodes, so it has at least N - 1 distinct undirected
edges. Each gate supplies at most 2 distinct such edges through its two
input pins; repeated pins can only merge edges, never create new ones, so
they cannot increase this upper bound. Constant nodes supply no input-pin
edges. Hence

    2g >= n + g + k - 1,

i.e. g >= n - 1 + k >= n - 1.

Since C was arbitrary, every B2 circuit for f_n has at least n - 1 gates:

    C_B2(f_n) >= n - 1   for all n >= 1.

The small cases check out directly: f_1(x_1) = x_1 is computed with 0 gates
by outputting the input wire (output selection costs zero), so
C_B2(f_1) = 0 = 1 - 1; f_0 is the constant 0, computed with 0 gates, so
C_B2(f_0) = 0.

## 5. Upper bound and tightness: C_B2(f_n) = n - 1 for n >= 1

A full binary tree of 2-input OR gates over the n inputs computes f_n
exactly and uses precisely n - 1 OR gates (each OR is the truth table 0111,
one of the 16 allowed). For n = 1 the tree is empty (0 gates). Hence
C_B2(f_n) <= n - 1 for all n >= 1, and combined with Section 4,

    C_B2(f_n) = n - 1   for every n >= 1, and C_B2(f_0) = 0.

In particular the ratio C_B2(f_n)/n tends to exactly 1, so no coefficient
larger than 1 at exponent 1 can hold for this family.

## 6. The declared (1, 1) guarantee

Put g(n) = 1 * n^1 * (log2 n)^0 * (log2 log2 n)^0 = n, matching the exact
rational parameters in `bound.json` (coefficient 1/1, exponent 1/1, log
exponents 0). For n >= 1,

    C_B2(f_n) / g(n) = (n - 1) / n = 1 - 1/n  ->  1   as n -> infinity,

so liminf over all integer n of C_B2(f_n)/g(n) = 1 >= 1. (The ratio is taken
for n >= 1; the single n = 0 term, where g(0) = 0, does not affect the
liminf.)

Equivalently, for every epsilon > 0 set N_epsilon = max(2, ceil(1/epsilon)).
For every integer n >= N_epsilon we have C_B2(f_n) = n - 1 and

    (n - 1) - (1 - epsilon)*n = epsilon*n - 1 >= 0,

hence C_B2(f_n) >= (1 - epsilon)*n = (1 - epsilon)*g(n). The same fixed
family, verifier V, and witness bound W work for every epsilon and every
sufficiently large integer length n = |x|; there is no padding or encoding
loss since f_n is defined directly on n input bits, and witness bits are not
circuit inputs. The finitely many small lengths n < N_epsilon are irrelevant
to the asymptotic quantifiers. This is exactly the joint declaration in
`bound.json`: full guarantee polynomial_log with (c, p, q, r) = (1, 1, 0, 0)
and finite projection (1, 1).

## 7. Remarks (not part of the scored joint theorem)

The repository's starting baseline records a published (31/10, 1) claim for a
different (BK/Shoup) family. Its proof relies on external results (Li-Yang
ECCC TR21-023 Theorem 1.1, an affine-disperser preprint Theorem 2.2, and
Shoup's deterministic irreducible-polynomial Theorem 3.2) whose proofs are
not reproduced in this submission, so that claim is not asserted here. The
only scored theorem of this submission is the (1, 1) OR-family result proved
in Sections 1-6 above. The tightness identity C_B2(f_n) = n - 1 is recorded
as an unscored lemma in `submission.json`.
