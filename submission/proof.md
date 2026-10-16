# OR on every input length

## One uniform total family and NP membership

For every integer n >= 0 and x in {0,1}^n, define f_n(x) = 1 exactly when
some bit of x is 1. In particular f_0(empty) = 0 and f_1(x_1) = x_1.
Let L be the set of all finite binary strings containing a 1. This is a
single total language, with no advice, parameter choice, oracle, or
per-length construction. Its actual input length is n = |x|.

Use the fixed witness bound p(n) = 0. The fixed deterministic verifier
V(x,w) rejects a nonempty w, and otherwise scans x and accepts if and
only if it encounters a 1. It runs in O(|x| + |w| + 1) time.
If x is in L, the empty witness has length at most p(|x|) and is accepted.
Conversely, any accepted witness is empty and acceptance implies that x
contains a 1, hence x is in L. If x is empty or all zero, no witness is
accepted. This proves both directions of the NP definition for all lengths;
the same scan also shows L is in P. Witness bits are not circuit inputs.

## Exact gate count in the full B2 model

Allow deterministic acyclic one-output circuits with every binary Boolean
truth table, two pins per gate, arbitrary depth and fan-out, sharing,
repeated pins, and free input nodes, wires, output selection and constants.
Count each retained gate once, including unary, constant and degenerate
truth tables. Thus this argument also covers the free constants used
during Li-Yang substitutions; it does not require a restricted basis or
a gate-elimination normalization.

Fix an integer n >= 2 and any circuit computing f_n exactly on all inputs.
Each of its n input variables is essential: changing just x_i from 0 to 1
at the all-zero input changes f_n from 0 to 1. A variable with no directed
path to the output cannot affect that output, by induction along a
topological order. Therefore all n variables occur among output ancestors.
The output cannot be a constant or a single input, since either would
depend on at most one variable.

Delete all nodes outside the ancestors of the output. Let g be the number
of remaining gates and k the number of remaining free constant nodes.
Identify any repeated occurrences of an input terminal with the same
variable node. This does not change any gate or computed value. Forget
directions and merge parallel edges to obtain a simple undirected graph.
It is connected: each retained node has a path to the output. It has
n + g + k vertices. A finite connected graph with v vertices has at least
v - 1 distinct edges: starting at one vertex, each newly reached vertex
requires an edge, giving a spanning tree with v - 1 edges.

Every edge enters some gate before directions are forgotten. The two
pins at each gate supply at most 2g distinct edges in total, even when
pins coincide or descendants are shared. Consequently

    n + g + k - 1 <= number of edges <= 2g,
    g >= n + k - 1 >= n - 1.

The original circuit has at least g gates. This argument never restricts
gate truth tables, so it includes XOR, negation, projections, constants
and every degeneracy. A retained negated-projection gate is charged once;
absorbing a negation into a descendant's truth table changes no counting
step. Unary and constant gates still have at most two incoming edges.
No uncharged unary computation, monotonicity assumption, formula
assumption, depth bound, or direct-sum claim is used.

For the matching upper bound, set y_1 = x_1 and, for i = 2,...,n, use one
binary OR gate y_i = OR(y_(i-1),x_i); select y_n as output. This acyclic
B2 circuit has exactly n - 1 gates. At n = 1 select the input with no
gate. At n = 0 the defined function is the constant-zero output, covered
by the constant-output convention. Therefore, for every integer n >= 1,

    C_B2(f_n) = n - 1.

## All-length asymptotics

For every integer n >= 1, C_B2(f_n)/n = 1 - 1/n, so the limit and liminf
are both 1. Equivalently, given any epsilon > 0, take
N_epsilon = max(2, ceil(1/epsilon)). For every integer n >= N_epsilon,
C_B2(f_n) = n - 1 >= (1 - epsilon)n. The remainder r(n) = 1 satisfies
r(n)/n -> 0. This covers all sufficiently large integer lengths with one
fixed family and verifier, and needs no finite experiments or citations.
