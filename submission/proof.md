# OR_n: C_B2(OR_n) >= n-1 at every length n >= 2, hence the pair (1, 1)

This proof is fully self-contained. It cites no external theorem: every step is
proved below.

## 1. The family and its uniform NP membership

L is the set of binary strings that contain at least one 1. Thus f_n = OR_n:
{0,1}^n -> {0,1}, f_n(x) = 1 iff some bit of x is 1. Every length is defined: in
particular f_0(empty) = 0, since the empty string contains no 1, and f_1(x) = x_1.

Evaluator A: scan x once and output 1 iff a 1 is seen. This is one fixed
deterministic algorithm running in O(n) time at every length.

NP verifier: witness bound W(n) = 0 and V(x, w) = 1 iff w is the empty string and
A(x) = 1. V runs in linear time.
- If x is in L, then A(x) = 1, so w = empty (|w| = 0 <= W) gives V(x, w) = 1.
- If V(x, w) = 1 for some w with |w| <= W(|x|) = 0, then A(x) = 1, so x is in L.

Both directions hold, so L is in P and hence in NP. There is no advice, no
per-length choice and no nonconstructive step. Witnesses are not circuit inputs:
the circuits below read exactly the n bits of x.

## 2. Circuit model

We use the general B2 model fixed in rules.md. A circuit C is a finite directed
acyclic graph whose nodes are:
- the n input nodes x_1..x_n, each with no incoming pins;
- possibly constant nodes 0/1, with no incoming pins and costing zero;
- gates, each with exactly two input pins.

Each pin is a wire from an earlier node; both pins may come from the same node.
A gate computes any of the 16 functions {0,1}^2 -> {0,1}, including degenerate
and trivial ones. Fan-out and depth are unrestricted. One node is selected as the
output, and C must compute f_n exactly on all 2^n inputs. The size |C| is the
number of gates. Inputs, constants, wires and output selection cost zero. This
proof only needs the facts listed here: two pins per gate, acyclicity, and that
only gates are counted. So it applies verbatim to every circuit in the model,
including circuits with redundant, degenerate or unreachable gates, repeated pins,
constant nodes, or negations absorbed into truth tables. No normalization is used.

## 3. Lower bound: every circuit for OR_n has at least n-1 gates (n >= 2)

Fix n >= 2 and any B2 circuit C computing OR_n, with output node o.

Step 1 (every input is essential). For each j, let e_j be the input with a single 1
in position j. Then OR_n(0...0) = 0 and OR_n(e_j) = 1, so the output value depends
on x_j.

Step 2 (every input is an ancestor of the output). Let A be the set of nodes with a
directed path to o, including o itself. Suppose x_j is not in A. By induction in
topological order, the value of each node in A is a function of the values of input
nodes in A only. (An input or constant node has its own value. A gate in A has both
pin sources in A, because a path from the source to o runs through the gate, and
its value is a function of their values.) Hence the value at o would be the same on
0...0 and on e_j, which differ only in x_j. This contradicts Step 1. So all n input
nodes lie in A.

Step 3 (counting). Let h be the number of gates in A and k the number of constant
nodes in A, so A has n + h + k nodes. Form the undirected multigraph H on vertex set
A with one edge per pin of each gate in A, joining the gate to its pin's source,
which lies in A by Step 2. H has exactly 2h edges, counting parallel edges from
repeated pins separately. H is connected, because every node of A has a directed
path to o along pin edges, and that path lies in A. A connected multigraph with N
vertices contains a spanning tree with N - 1 distinct edges, so its edge count is
at least N - 1. Therefore

    2h >= (n + h + k) - 1  >=  n + h - 1,   so   h >= n - 1 + k >= n - 1.

Since every gate of A is a gate of C, |C| >= h >= n - 1. As C was arbitrary,

    C_B2(OR_n) >= n - 1   for every integer n >= 2.        (*)

(For n = 1, C_B2 = 0 because the output can be the input node. This does not affect
the asymptotic statement.)

Remark (tightness, not needed for the claim). A balanced tree of n - 1 binary OR
gates computes OR_n, so C_B2(OR_n) = n - 1 for n >= 2. The declared coefficient 1 is
therefore exactly the true constant, and no larger coefficient holds for this family.

## 4. The declared asymptotic guarantee

The declaration is the polynomial_log guarantee g(n) = 1 * n^1 * (log2 n)^0 *
(log2 log2 n)^0 = n, with finite projection (c, p) = (1, 1).

Given eps > 0, let N_eps = max(2, ceil(1/eps)). For every integer n >= N_eps, (*)
gives C_B2(OR_n) >= n - 1 = (1 - 1/n) n >= (1 - eps) n = (1 - eps) g(n). This holds
at every integer length, not a subsequence. Equivalently, C_B2(OR_n)/g(n) >= 1 - 1/n
for all n >= 2, so liminf over all integers n of C_B2(OR_n)/g(n) >= 1. With
c = 1 > 0 and p = 1 this is the declared full guarantee. Its finite projection
liminf C_B2(OR_n)/n^1 >= 1 follows immediately for the same fixed family, verifier
and witness bound, so both score coordinates come from this one theorem.

## 5. Scope

This is a deliberately modest, fully proved entry. It does not approach the
literature value 31/10 and claims no new circuit lower-bound technique. It matches
the "OR: coefficient 1 (positive)" case in scientific-cases.md, written out in full.
