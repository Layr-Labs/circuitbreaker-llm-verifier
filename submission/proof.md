# OR positive control: exact B2 complexity and coefficient 1

This submission is a self-contained positive control for CircuitBreaker's
hosted review. It follows the OR coefficient-one case in the organizer's
`scientific-cases.md` at source revision
`bcd3c4eae718b65e5e53ce6506ffe81919c384e3`. It claims no improvement over the
published coefficient 31/10 in the starting submission.

## 1. One fixed total family and an NP verifier

Let L be the set of finite binary strings that contain at least one 1. For
every integer n >= 0 and every x in {0,1}^n, define

    f_n(x) = 1 if some input bit is 1, and 0 otherwise.

In particular, f_0(empty) = 0 and f_1(x_1) = x_1. No lengths are omitted.
The actual input length is n; there are no encodings, blocks, padding bits,
advice strings, field elements, or additional circuit inputs.

Fix the following deterministic verifier V for all lengths. Reject if the
witness w is nonempty. Otherwise scan x from left to right and accept if and
only if a 1 occurs. Set the polynomial witness bound p(n) = 0. This algorithm
has running time O(|x|+|w|+1), including reading or rejecting its witness.

If x belongs to L, the empty witness has length zero, is within p(|x|), and
V accepts it. Conversely, if a witness of length at most p(|x|) makes V
accept, it is empty and the scan found a 1 in x, so x belongs to L.
Thus both directions of NP membership hold, with one fixed verifier.
The same scan, with no witness, evaluates f_n in linear time, so L is in P.
The witness is never an input to a circuit computing f_n.

## 2. Model and essential input bits

We use exactly the deterministic acyclic, single-output B2 circuits in
`rules.md`: two input pins per gate, any of the 16 binary Boolean tables,
arbitrary depth and fan-out, repeated pins, and shared subgraphs. Inputs,
wires, and output selection are free. Each retained gate is counted once.
Free constant nodes, when present under the stated substitution convention,
are allowed in the lower-bound argument.

For each n >= 1, every one of the n input bits is essential for f_n: the
all-zero assignment yields 0, while setting only bit i to 1 yields 1.
In any exact circuit computing f_n, input i must therefore have a directed
path to the selected output node. Indeed, a node with no input-i ancestor
computes a function independent of input i, as induction in a topological
ordering shows. The same is true of a selected input or constant output.

This dependence argument uses no restriction on a gate's truth table, and
does not assume different inputs have disjoint downstream computations.

## 3. Gate lower bound, including degenerate circuits

Fix n >= 2 and any exact B2 circuit for f_n. Delete nodes and gates that
are not ancestors of its output, retaining the output itself. Write g for
the number of gates in the retained graph and k for its constant nodes.
All n input vertices remain, by essentiality. The original circuit has
at least g gates, so a lower bound for this retained graph suffices.

For n >= 2 the output cannot be an input or a constant: either choice
depends on at most one input, whereas f_n depends on all n. Thus the
retained graph has an output gate, along with its ancestors.

Forget edge directions and collapse parallel edges caused by repeated
pins. The resulting finite undirected graph is connected: every vertex
has a directed path to the output before directions are forgotten, and
collapsing parallel edges preserves connectivity. It has

    v = n + g + k

vertices. If e is its number of distinct undirected edges, connectedness
implies e >= v-1. For completeness, start with one vertex and repeatedly
add an edge from the visited set to an unvisited vertex. Such an edge
exists until all vertices are visited, by connectedness. The v-1 added
edges are distinct, proving the inequality.

Every edge comes from an input pin at some retained gate, since inputs
and constant nodes have no incoming pins. There are at most two such
pins per gate. Collapsing parallel edges can only reduce the edge count,
so e <= 2g. Therefore

    n + g + k - 1 <= e <= 2g,
    g >= n + k - 1 >= n - 1.

Every original circuit for f_n consequently has at least n-1 gates.

Here are the model boundary cases in this counting argument. A shared
descendant is one vertex and is counted once, even when fan-out is
unbounded. Repeated pins produce no extra distinct edge. A constant
truth-table gate, projection, negated projection, or other degenerate
gate remains a counted gate if retained; it still has at most two pins.
A free constant node contributes to k and is never counted as a gate.
Allowing those nodes only strengthens the displayed counting inequality.
No substitution or normalization is performed by this proof, so it
does not assume that constant gates or terminal negations can be removed
for free. Negation absorbed into a descendant truth table changes no
dependency edges needed for the argument. No additional basis, free
computational gate, or direct-sum additivity is introduced.

The remaining claimed exact-complexity case is n = 1. Selecting x_1
as output computes f_1 with zero gates; no circuit can have a negative
gate count. Thus C_B2(f_1) = 0 = n-1. We have already defined f_0
and its evaluator. The theorem makes no gate-complexity assertion
at n = 0, so no zero-input output convention is needed for the bound.

## 4. Matching upper bound

For n >= 2, take a binary OR tree with leaves x_1,...,x_n. Each internal
vertex is a two-input OR gate, one of the permitted B2 tables. A chain
is an explicit such tree: compute x_1 OR x_2, then OR its result with
x_3, and continue until x_n. There is one gate for each index 2,...,n,
giving exactly n-1 gates and a single output. This circuit computes
f_n on every input. For n = 1 use the zero-gate output selection above.

Together with Section 3 this proves

    C_B2(f_n) = n-1 for every integer n >= 1.

## 5. Exact coefficient and asymptotic quantifiers

The declared rational coefficient is 1/1, whose numerator and denominator
are canonical nonnegative/positive decimal integers and are relatively
prime. It has an exact finite nonzero binary64 representation.

For every integer n >= 1,

    C_B2(f_n)/n = 1 - 1/n.

Hence the limit, and therefore the liminf, equals 1. Explicitly, for
each real epsilon > 0, choose

    N_epsilon = max(1, ceil(1/epsilon)).

For every integer n >= N_epsilon, 1 <= epsilon*n and consequently
C_B2(f_n) = n-1 >= (1-epsilon)*n. The remainder is the constant 1,
which divided by n tends to zero. The family, evaluator, verifier,
witness bound, and circuit model are unchanged as epsilon varies.
This covers every sufficiently large integer length, not just a
subsequence. The upper bound also shows that this family's limiting
coefficient cannot exceed 1 and does not reach the superlinear endpoint.

## 6. Evidence scope

The proof above establishes the mathematical claim directly for all
circuits and all stated lengths. It does not depend on finite searches,
an LP certificate, the baseline's affine-disperser construction, or an
external circuit theorem. The only organizer sources used are
`rules.md` (model and score definition) and the section
“OR: coefficient 1 (positive)” in `scientific-cases.md`, both at the
source revision identified above. The argument includes that case's
essential-input graph proof and makes its small-length, NP-verifier,
sharing, constant-node, degenerate-gate, and epsilon cases explicit.
