# Exact gate complexity of the uniform parity family

## Scope and status

This is an elementary, self-contained theorem, not a claim of new research or a stronger bound than the published 31/10 seed. The claimed joint projection is (1,1). Every argument below concerns the one fixed family defined here. No prior accepted submission, numerical experiment, LP certificate, or external theorem is a premise.

## One uniform total language, including every length

For every integer n >= 0 and x in {0,1}^n, define f_n(x) as the XOR of all input bits. The empty XOR is zero. Let L contain exactly the binary strings of odd parity. A single deterministic algorithm initializes a bit b to zero, scans x from left to right, replaces b by b XOR x_i, and returns b. Induction on the number of scanned bits proves this computes f_n at every length, including zero. It takes O(n+1) time and constant work space on a sequential input model.

Set the witness bound W(n)=0. Define V(x,w) to reject every nonempty w and otherwise run the algorithm just given and accept if and only if its output is one. The empty string is the only witness of length at most W(n). If x is in L, that witness is accepted. Conversely, if an allowed witness is accepted, the algorithm returned one, so x is in L. Thus L is in NP with a fixed polynomial-time verifier and a fixed polynomial witness bound; indeed it is in P. There is no advice, length-dependent choice, padding, field encoding, or auxiliary input. The circuit input length is exactly n=|x|.

## Circuit model and syntactic reachability

Consider any deterministic acyclic circuit computing f_n in the full B2 model of rules.md. Each retained gate has two incoming pins and may use any of the 16 binary Boolean truth tables; depth and fan-out are unrestricted. Repeated pins and arbitrary shared descendants are allowed. Count each retained gate once. Inputs, wires, and output selection cost zero. We also allow free constant source nodes in the following argument, so constants introduced by normalization cannot invalidate the lower bound.

A variable is syntactically reachable from the selected output when there is a directed path from its input source to that output; selecting that input itself counts as a path of length zero. If an input variable is not reachable, changing only that variable cannot change the output. To prove this assertion, put the gates in a topological order (an order exists for a finite acyclic graph: repeatedly remove a vertex with no predecessors; if none existed, following predecessors in a finite graph would produce a directed cycle). Induct through this order in the subgraph of ancestors of the selected output. Each source in that subgraph is either a constant or an input other than the absent variable, and hence has unchanged value. Each gate applies its fixed truth table to two unchanged pin values and hence also has unchanged value. The output therefore has unchanged value. This proof uses no property of a particular basis and works equally for constant gates and truth tables that ignore one or both pins.

For n>=1, parity depends essentially on every input: the all-zero assignment has output zero, while changing only bit i to one yields output one. Thus every input variable must be syntactically reachable from the output of any correct circuit.

## Fan-in counting lemma with sharing and degeneracies

Suppose the selected output is a gate, and retain just its syntactic ancestors, including that gate and every source reached by following its incoming pins recursively. Let G be the number of retained gates in this ancestor subgraph, S the number of source vertices, and E the number of directed pin edges counted with multiplicity. In particular, two repeated pins give two edges. A source vertex is an input or a free constant; a retained constant-valued gate still counts as a gate, not as a free source. All pins of each retained gate have their predecessor present. Hence E=2G. If a normalized representation omits ignored pins or uses explicit unary nodes, the same argument gives E<=2G, which is sufficient.

Every vertex of this subgraph has a directed path to the selected output. Forget directions and multiplicities; the resulting simple undirected graph on G+S vertices is connected. Any connected finite undirected graph on v>=1 vertices has at least v-1 distinct edges: start with one visited vertex and, while some vertex is unvisited, connectivity provides an edge from the visited set to its complement; add that edge and its new endpoint. Exactly v-1 such steps visit all vertices, and all chosen edges are distinct. Restoring multiplicities cannot reduce the edge count. Therefore

    G + S - 1 <= E <= 2G,

so G >= S-1. If all n input variables are reachable, their n distinct input sources occur in the subgraph, so S>=n and G>=n-1. Extra constant sources only strengthen this counting inequality. If the representation permits duplicate source nodes carrying the same variable, at least n sources are still needed to represent n distinct variables. Repeated pins only increase edge multiplicity; sharing never invalidates connectivity or the bound on the number of pins.

The full circuit has at least G gates because deleting vertices outside the output's ancestor subgraph cannot increase gate count. Thus the original circuit has at least n-1 gates.

## Terminal cases and normalization

If the selected output is a free constant, it cannot compute parity for any n>=1. If it is an input wire, it depends on at most one variable and therefore cannot compute parity for n>=2. For n=1, a zero-gate circuit selecting x_1 is correct and the claimed lower bound is zero. If a terminal negation is represented as a retained negated-projection gate, it is counted once and is covered by the preceding gate argument. The model does not grant an uncharged terminal computational operation. Negations absorbed into descendant truth tables remain covered because all truth tables were allowed from the start.

No substitution or normalization is required in this proof. Nevertheless the counting argument explicitly permits free constant nodes and gates with fewer than two effective inputs, and the reachability argument covers all constant and unary functions. Consequently Li-Yang normalization does not introduce a basis restriction or an uncounted resource into the argument. We never infer that the number of gates is a potential, assume a direct sum, restrict fan-out, or turn the circuit into a formula.

Combining these cases, every full-B2 circuit computing f_n has at least n-1 gates for every n>=1.

## Matching upper bound and all-length asymptotics

For n=1, select the input x_1. For n>=2, use a chain of n-1 XOR gates: first x_1 XOR x_2, then XOR that output with x_3, continuing until x_n. XOR is one of the allowed binary truth tables. Induction on chain length proves its output is f_n. This gives C_B2(f_n)<=n-1. Together with the lower bound,

    C_B2(f_n)=n-1  for every integer n>=1.

The n=0 family member was defined as zero and is irrelevant to the asymptotic assertion. In particular, for positive n the ratio C_B2(f_n)/n equals 1-1/n, so its limit and liminf over all integer lengths are one. Explicitly, for any real epsilon>0 choose N_epsilon=max(1,ceil(1/epsilon)). Then for every integer n>=N_epsilon, 1/n<=epsilon and C_B2(f_n)>=(1-epsilon)n.

In bound.json the full guarantee is g(n)=1*n^1*(log2 n)^0*(log2 log2 n)^0=n for sufficiently large n. The declared finite projection is exactly coefficient 1 and exponent 1, both attached to this entire same theorem. It is not an endpoint, a sparse-length claim, an improved literature record, or a claim that NP is outside P/poly. No finite computation is needed to establish the universal assertion.
