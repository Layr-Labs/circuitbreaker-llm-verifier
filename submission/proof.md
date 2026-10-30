# OR at every length: coefficient 1, exponent 1

## Fixed total family and uniform NP verifier

Let L be the set of finite binary strings containing at least one 1. For every
integer n >= 0, let f_n:{0,1}^n -> {0,1} be its characteristic function.
Thus f_0(empty)=0, f_1(x_1)=x_1, and f_n(x)=OR(x_1,...,x_n) for n>=2.
The parameter n is exactly the number of input bits, with no encoding or padding.

Fix W(n)=0. Define a single deterministic verifier V(x,w): reject if w is
nonempty; otherwise scan x and accept if and only if some bit equals 1.
This algorithm takes O(|x|+|w|+1) time and is the same at all lengths.
If x belongs to L, the empty witness has length at most W(|x|) and is accepted.
Conversely, any accepted witness is empty and V's scan has found a 1 in x,
so x belongs to L. This proves both directions of the required NP equivalence.
In fact L is in P. Witness bits are never inputs to the circuits below.

## Circuit model

Consider any deterministic finite acyclic circuit with one designated output,
inputs x_1,...,x_n, and gates with two ordered input pins computing any of the
16 binary Boolean functions. Depth and fan-out are unrestricted. Repeated pins
and shared subgraphs are permitted. Each retained gate costs one, including
constant-function gates and gates computing a projection or a negated projection.
Input nodes, wires and selection of the output have cost zero. The argument also
allows cost-free constant source nodes, so it covers the constants allowed during
Li-Yang substitutions. It makes no assumption that the gates are nondegenerate,
monotone, XOR-free, or arranged as a formula.

## Essential inputs reach the output

For every n>=1 and every i in {1,...,n}, changing the all-zero input only in
coordinate i changes f_n from 0 to 1. Thus all n coordinates are essential.

In an acyclic circuit, the function at a node depends only on variables labelling
its ancestor input nodes. To prove this, take a topological ordering. An input
node depends only on its label and a constant node on no variable. At a gate,
its output is a fixed binary function of its two predecessor values; hence it
depends only on the union of their ancestor variables. Induction in the ordering
proves the assertion, including constant and unary gate truth tables.
Consequently, an output computing f_n must have an ancestor input node with
each of the n variable labels.

## Counting gates, with arbitrary sharing

Fix n>=2. Remove every node that is not an ancestor of the designated output
(including the output itself among its ancestors). This cannot increase the
number of gates or change the computed function. Let g be the number of retained
gates, m the number of retained input nodes, and k the number of retained
cost-free constant source nodes. If duplicate input nodes carry the same label,
keep them in this count; essentiality still gives m>=n.

The output cannot be an input or a constant source node: those depend on at most
one coordinate, whereas f_n depends on n>=2 coordinates. It is therefore a gate.
Forget edge directions and replace any parallel edges by a single edge. The
result is a finite connected undirected graph: every retained node has a directed
path to the output before directions are forgotten. It has v=m+k+g vertices.
It has at most 2g distinct edges, because each edge originally enters a gate
through one of its two pins. Repeated pins only reduce the number of distinct
edges. Sharing and unbounded fan-out do not change this incoming-pin count.

Every connected graph on v>=1 vertices has at least v-1 edges. Here is a proof:
start with the output vertex alone. While some vertex is not yet reached,
connectivity supplies an edge from the reached set to a vertex outside it
(take the first such edge on a path to any unreached vertex). Add that vertex
and edge. After v-1 additions this constructs v-1 distinct edges of the graph.
Applying this fact gives

    m+k+g-1 <= 2g,  hence  g >= m+k-1 >= n-1.

The original circuit has at least g gates, so every such circuit computing f_n
has at least n-1 gates. This argument counts degenerate and terminal gates
directly; it needs no elimination step, no free negation gate, and no assumption
about disjoint descendants. In particular it covers the full B2 model in rules.md.

## Matching construction and all-length asymptotics

For n>=2, first compute x_1 OR x_2, then successively OR the result with x_3,
through x_n. Each OR is an allowed binary truth table; there are n-1 gates.
For n=1 select x_1 as output at zero cost. The empty-length value was defined
above; its circuit convention has no effect on a limit at infinity.
It follows that C_B2(f_n)=n-1 for every n>=1.

For every real epsilon>0 choose the integer N_epsilon=max(2,ceil(1/epsilon)).
For EVERY integer n>=N_epsilon,

    C_B2(f_n)/n = 1-1/n >= 1-epsilon.

Thus liminf over all integer lengths of C_B2(f_n)/n is 1. In bound.json the
full guarantee is g(n)=1*n^1*(log2 n)^0*(log2 log2 n)^0=n for sufficiently
large n, and its finite projection is exactly [1,1]. The same fixed family
and verifier establish the entire joint declaration in both source claims.

## Attribution and evidence scope

This is the elementary essential-input connectivity argument, also recorded in
the challenge's scientific-cases.md, section "OR: coefficient 1 (positive)"
at repository commit 2cd3a17f3f296b7f0ac248a63a9b54a0fabdf851, lines 8-20.
The complete argument is reproduced above; that file is not needed as an
external premise. This submission claims no novelty and relies on no finite
experiment, checker output, prior accepted submission, or unproved external result.
