# Coefficient 1 at exponent 1 for the total OR family: a complete, self-contained proof

This submission declares the pair (c, p) = (1, 1) for one fixed total family and proves it
unconditionally in the pinned model. It uses no external theorem, no citation and no
computation: every premise is discharged below. It does not claim a new research result; it
is the elementary connectivity bound for a function that depends on all of its inputs, stated
and proved exactly in the form the rules require (family, NP membership, model, bit
accounting, small lengths, quantifiers).

## 1. The fixed total family

Let L be the set of binary strings that contain at least one symbol 1. For every integer
n >= 0 define f_n : {0,1}^n -> {0,1} by

    f_n(x_1, ..., x_n) = x_1 OR x_2 OR ... OR x_n,

with f_0 = 0 (the empty string contains no 1, so it is not in L). The family is total: f_n is
defined on every input of every length. The input length is the actual bit length n = |x|.
There are no blocks, field encodings, paddings, per-length parameters or advice strings; one
rule defines every length, including the small exceptions n = 0 and n = 1 (f_1(x_1) = x_1).

## 2. Uniform NP membership (both directions)

Evaluator. The algorithm A reads x left to right and outputs 1 as soon as it reads a 1, and
outputs 0 if it reaches the end without reading a 1. It runs in time O(n) on every input and
computes f_n at every length, including n = 0 (it reads nothing and outputs 0).

Verifier and witness bound. Let W(n) = 0 and let V(x, w) = 1 iff w is the empty string and
A(x) = 1. V runs in polynomial time.

(=>) If x in L then x contains a 1, so A(x) = 1 and V(x, empty) = 1 with |empty| = 0 <= W(|x|).
(<=) If V(x, w) = 1 for some w with |w| <= W(|x|) = 0, then w is empty and A(x) = 1, so x
contains a 1, i.e. x in L.

Hence L is in P, therefore in NP, with the same fixed algorithms at every length. The witness
is empty and is never an input to the circuits below; the circuit inputs are the n bits of x.

## 3. The circuit model

We use the general B2 model fixed by the rules: a circuit is a finite directed acyclic graph
whose nodes are input nodes (one for each of the n input bits), constant nodes (labelled 0 or
1), and gates. Every gate has exactly two input pins; each pin is attached to one node (an
input node, a constant node, or a gate), and the two pins of a gate may be attached to the
same node. Each gate carries one of the sixteen binary Boolean truth tables. Depth and fan-out
are unrestricted and subgraphs may be shared. One node is designated as the output. The
circuit computes f_n when, for every x in {0,1}^n, the value at the output node equals f_n(x).
The cost C_B2(f_n) is the minimum, over all circuits computing f_n, of the number of gates;
every retained gate costs one, including degenerate gates (gates whose truth table ignores a
pin, computes a constant, a projection or a negated projection), while input nodes, constant
nodes, wires and output selection cost zero. Nothing in the argument below depends on the
particular truth tables, so no normalization of degenerate or constant gates is needed: we
count every gate that is present, which is exactly what the model charges.

## 4. The lower bound

**Theorem.** For every integer n >= 1, every B2 circuit computing f_n has at least n - 1
gates. Hence C_B2(f_n) >= n - 1 for every n >= 1.

*Proof.* Let n >= 1 and let C be a B2 circuit computing f_n with output node o.

Step 1 (every input is essential). For 1 <= i <= n, f_n(0, ..., 0) = 0 while f_n(e_i) = 1,
where e_i is the input with a single 1 in position i. So changing x_i alone changes the value
of f_n. If no directed path led from the input node of x_i to o, the value at o would be a
function of the other inputs only, because the value at every node is determined by the
values at the nodes attached to its pins, and by induction along the acyclic order no node
reachable to o would read x_i; this contradicts the displayed pair of inputs. Hence for every
i there is a directed path from input node i to o. In particular, for n >= 2 the output o is a
gate (an input node or a constant node as output would compute a function of at most one
variable or a constant, contradicting Step 1 for the other inputs; for n = 1 the output may
be the input node itself, and the claimed bound n - 1 = 0 holds trivially). From now on let
n >= 2.

Step 2 (trim). Let C' be the subcircuit consisting of all nodes from which o is reachable by
a directed path (including o), with all pins among them. C' is again a B2 circuit computing
f_n with the same output: deleting nodes that cannot reach o does not change any value on a
path to o. C' contains all n input nodes by Step 1, some number k >= 0 of constant nodes, and
some number g' of gates, with g' <= g, the number of gates of C. It suffices to prove
g' >= n - 1.

Step 3 (connectivity). Consider the undirected multigraph G whose vertices are the nodes of C'
and whose edges are the pins of the gates of C' (one edge from each pin to the node it is
attached to; two pins attached to the same node give two parallel edges, and a pin attached
to the gate itself is impossible since the graph is acyclic). Every vertex of C' is connected
to o in G, because every node of C' has a directed path to o and each step of such a path is
a pin of a gate of C', hence an edge of G. So G is connected. The number of vertices of G is
n + k + g' (n input nodes, k constant nodes, g' gates). The number of edges of G is at most
2 g', since every edge is one of the two pins of one of the g' gates (pins of deleted gates
are gone with them, and input and constant nodes have no pins).

Step 4 (count). A connected multigraph on V vertices has at least V - 1 edges: a spanning tree
of a connected graph has exactly V - 1 edges and is a subgraph of it. Therefore

    2 g'  >=  (number of edges of G)  >=  (n + k + g') - 1,

so g' >= n + k - 1 >= n - 1. Since g >= g', the circuit C has at least n - 1 gates. As C was
an arbitrary circuit computing f_n, C_B2(f_n) >= n - 1. This completes the proof. QED

Remarks on the model. Degenerate gates are counted like any other gate; they only help the
bound (each still supplies two pins, i.e. at most two edges, and costs one). Repeated pins
cannot increase the edge count above 2 g'. Constant nodes cost zero but are vertices of G,
which is why the bound is g' >= n + k - 1; dropping k >= 0 gives the stated n - 1. Shared
subgraphs and unrestricted fan-out are allowed and are already covered, since the counting
uses only the number of pins. No step uses a substitution, so no normalization of constants
arising under substitutions is needed.

## 5. Small lengths and the quantifiers

n = 0: f_0 = 0 is a constant; the theorem is stated for n >= 1 and the asymptotic statement
below does not involve n = 0. n = 1: f_1(x_1) = x_1 is computed by the circuit whose output is
the input node itself, with 0 gates, and the bound n - 1 = 0 holds. For every n >= 2 the
theorem gives C_B2(f_n) >= n - 1.

Asymptotic statement. With c = 1 and p = 1, g(n) = c * n^p = n. For every integer n >= 1,
C_B2(f_n) / n >= (n - 1) / n = 1 - 1/n. Hence for every epsilon > 0 and every integer
n >= N_epsilon := max(2, ceil(1/epsilon)), C_B2(f_n) >= (1 - epsilon) * n, which is exactly the
required form C_B2(f_n) >= (c - epsilon) * n^p for every integer n >= N_epsilon. Equivalently
liminf_{n -> infinity} C_B2(f_n) / n^1 >= 1. The remainder r(n) = 1 satisfies r(n)/n -> 0. The
bound holds at every integer length, not on a sparse subsequence, and the same family and the
same algorithms are used throughout.

## 6. Exactness of the coefficient (why no larger coefficient is claimed)

The binary OR tree of n leaves uses exactly n - 1 two-input OR gates and computes f_n, so
C_B2(f_n) <= n - 1 for n >= 1. Together with the theorem, C_B2(f_n) = n - 1 exactly, and
C_B2(f_n)/n -> 1. The declared pair (1, 1) is therefore valid, and no coefficient above 1 is
claimed for this family.

## 7. Summary of declared obligations

- Family: one total language L (strings containing a 1), f_n = OR_n, f_0 = 0; every length defined.
- NP membership: polynomial-time evaluator A, verifier V with empty witness, both directions shown.
- Model: Li-Yang general B2 as fixed by rules.md; every retained gate counted once, degenerate
  gates included; inputs, constants, wires and output selection free.
- Bit accounting: the circuit inputs are exactly the n = |x| bits; no encoding or padding.
- Bound: C_B2(f_n) >= n - 1 for every integer n >= 1 (connectivity count).
- Quantifiers: for every epsilon > 0 and every integer n >= max(2, ceil(1/epsilon)),
  C_B2(f_n) >= (1 - epsilon) n; liminf C_B2(f_n)/n >= 1.
- Declared pair: coefficient 1, exponent 1, no logarithmic factors; projection (1, 1).
