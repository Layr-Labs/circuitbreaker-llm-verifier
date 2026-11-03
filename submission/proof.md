# OR family: tight coefficient 1 at exponent 1

## Fixed total family

Define the total Boolean family OR_n : {0,1}^n → {0,1} for all n ≥ 0. For n ≥ 1, OR_n(x) = 1 iff at least one input bit x_i is 1. For n = 0, f_0(empty) = 0. This defines a single language: the set of nonempty binary strings containing at least one 1 bit, together with the empty string mapped to 0.

## Uniform P and NP membership

A deterministic algorithm computes OR in linear time: scan the n input bits and output 1 as soon as a 1 is found, else 0. This is polynomial in n. The language L = { x ∈ {0,1}* : x contains a 1 } is therefore in P, hence in NP. Use V(x,w) = 1 iff w is empty and the linear-time evaluation A(x) = 1. The witness bound W(n) = 0 is polynomial. Both directions are immediate: if x ∈ L then A(x) = 1, so V(x, ε) = 1; conversely if w is empty and V(x,w) = 1 then A(x) = 1, so x ∈ L. Every integer n ≥ 0 is covered, including n = 0. The same verifier and witness bound work for all lengths. No per-length advice, truth-table oracle or nonconstructive existence is used.

## Full B2 circuit model

Li–Yang's general B2 model applies: deterministic acyclic circuits, one output, two input pins per gate, all 16 binary Boolean truth tables, unrestricted depth and fan-out, exact computation on every input, nonuniform circuits allowed. Inputs are x alone (witness bits are not circuit inputs). Count each retained gate once. Inputs, wires and output selection cost zero. Constant nodes during substitutions cost zero and are handled by Li–Yang normalization. Negations can be absorbed into descendant truth tables; a retained negated-projection gate costs one. No free uncharged computational basis is assumed. These conventions follow ECCC TR21-023, 2021-02-20, Sections 1.1, 2.1, 2.4 (footnotes 3–4) and 3.3.

## Lower bound argument

Fix n ≥ 2 and consider any B2 circuit computing OR_n with g gates. Trim the circuit to the ancestors of the output gate: discard any gates that do not affect the output. The resulting trimmed graph is connected (every gate is reachable from the output via the input pins) and contains all n essential inputs (otherwise the circuit would ignore some input bit, contradicting correctness of OR on that bit).

Let the trimmed graph have g′ ≤ g gates and k constant nodes (zero or one constant 0 node, possibly one constant 1 node). The graph, with direction forgotten, has:

- Vertices: n input bits + g′ gates + k constants
- Edges: each of the g′ gates has exactly 2 input pins (each connecting to a predecessor vertex), giving at most 2g′ edges. Repeated pins do not increase this upper bound.

A connected graph with V vertices has at least V − 1 distinct edges. Therefore:

    n + g′ + k − 1  ≤  2g′

Rearranging:

    n + k − 1  ≤  g′

Since k ≤ 2 (at most one constant 0 node and one constant 1 node), we have:

    n + 2 − 1  ≤  g′   ⇒   n + 1  ≤  g′

But this is too crude. A more careful argument: the trimmed circuit for n ≥ 2 must satisfy g′ ≥ n − 1.

**Refined connectivity argument:** In the trimmed graph, consider the n input nodes (which have no incoming edges, only outgoing edges to gates). The graph has n source nodes. Contract each maximal connected subgraph that contains no input node into a single node. The remaining graph has exactly n input nodes and some internal nodes (gates). The number of edges is at most 2g′, and connectivity among the n inputs and g′ gates, together with the requirement that each of the n inputs feeds into the circuit, forces g′ ≥ n − 1.

A standard graph-theoretic proof: each gate has indegree 2, and the output node must be able to reach every input. Consider the underlying undirected graph. It has at least n + g′ vertices and at most 2g′ edges. A connected graph on V vertices requires at least V − 1 edges. Hence:

    n + g′ − 1  ≤  2g′
    n − 1  ≤  g′

Thus g′ ≥ n − 1, and since g ≥ g′, we have g ≥ n − 1.

## Matching upper bound

A binary OR tree using n − 1 two-input OR gates computes OR_n exactly. Each internal gate combines two values, and the tree has n − 1 internal nodes. This is an explicit B2 circuit with exactly n − 1 gates, so C_B2(OR_n) ≤ n − 1.

## Exact complexity

From the lower bound g ≥ n − 1 and the upper bound C_B2(OR_n) ≤ n − 1, we obtain:

    C_B2(OR_n) = n − 1   for every n ≥ 2.

For n = 0, C_B2(f_0) = 0 by convention (constant output needs no gates). For n = 1, OR_1(x) = x can be computed by a single wire (no gates) or by a projection gate costing 0 after normalization, so C_B2(OR_1) = 0.

## Declared guarantee

We have C_B2(OR_n) / n = (n − 1)/n → 1 as n → ∞. Therefore:

    liminf_{n→∞} C_B2(OR_n) / n = 1

Equivalently, for any ε > 0, there exists N_ε such that for every integer n ≥ N_ε, C_B2(OR_n) ≥ (1 − ε)×n. The declared pair (1, 1) with polynomial guarantee c = 1, p = 1, q = 0, r = 0 is established.

## Smaller coefficients

Any rational coefficient 0 < c ≤ 1 at exponent p = 1 follows from the same OR family, because (n − 1)/n ≥ c for sufficiently large n. These are valid smaller finite pairs but do not represent new scientific progress beyond the maximal coefficient 1.

## Sources

- LY-full: ECCC TR21-023, front-page date 2021-02-20; Sections 1.1, 2.1, 2.4 (footnotes 3–4) and 3.3 for B2 model conventions. https://eccc.weizmann.ac.il/report/2021/023/download/; SHA256 4530cbaf38f0f9a09704e94f7d1f57109895b71b92d0752a285e07d1109a4062
