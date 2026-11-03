# Exact full-B2 circuit size of the OR family

This is the known-answer example in `scientific-cases.md`, section “OR:
coefficient 1 (positive)”, at source revision
`2cd3a17f3f296b7f0ac248a63a9b54a0fabdf851`. It is not a new lower bound.
The following proof is self-contained and does not rely on a prior submission
or on an external theorem about circuit complexity.

## One fixed total family and its NP verifier

Let L be the set of finite binary strings containing at least one 1. For
every integer n >= 0 let f_n be its characteristic function on {0,1}^n:
f_n(x) is 1 exactly when some bit of x is 1. In particular f_0 on the empty
string is 0. A single deterministic algorithm scans its input and returns 1
if it sees a 1 and otherwise returns 0. Its running time is O(n+1).

Set the witness bound W(n)=0. The fixed verifier V(x,w) accepts exactly when
w is the empty string and the scan of x returns 1. If x is in L then the
empty witness has length at most W(|x|) and is accepted. Conversely, if a
witness of length at most W(|x|) is accepted, it is empty and the scan has
found a 1, so x is in L. Thus L is in NP; indeed it is in P. The witnesses
are not inputs to the nonuniform circuits below.

## Lower bound in the full B2 model

Fix any n >= 2 and any deterministic acyclic one-output B2 circuit that
computes f_n exactly on all n-bit inputs. Each gate has two input pins and
may compute any of the 16 binary Boolean functions. Depth and fan-out are
unrestricted, repeated pins and constant nodes are allowed, and every gate
costs one. Inputs, constant nodes, wires and output selection cost zero.

Discard all vertices and wires that are not ancestors of the selected output,
retaining the output itself. This cannot increase the number of gates and
does not change the computed function. For each i, the two inputs 0^n and
the string with only bit i set differ in output. Therefore every input x_i is
an ancestor of the output: if x_i were absent, changing it alone could not
change the output of the retained circuit. Since n >= 2, the output cannot
be an input or a constant, so it is a gate.

Forget edge directions and identify parallel wires between the same vertices.
Every retained vertex has a path to the output, so the resulting finite simple
undirected graph is connected. Let g be the number of retained gate vertices
and k the number of retained constant vertices; the n distinct input vertices
are also present. Thus the graph has n+g+k vertices and at least n+g+k-1
distinct edges: a connected finite graph contains a spanning tree on all its
vertices. Every edge comes from a wire ending at a gate input pin, and there
are exactly two pins per retained gate. Identifying parallel wires can only
decrease the number of distinct edges, so there are at most 2g edges. These
two inequalities imply

    n + g + k - 1 <= 2g,   hence   g >= n + k - 1 >= n - 1.

This argument counts shared subgraphs only once and permits arbitrary fan-out,
repeated pins, constants, unary or degenerate truth tables, and all binary gate
types. It does not assume a formula or monotonicity. The original circuit had
at least g gates, giving C_B2(f_n) >= n-1.

## Upper bound and joint asymptotic theorem

For n >= 2, start with x_1 and successively OR in x_2,...,x_n with one
two-input OR gate per step. This full-B2 circuit uses exactly n-1 gates and
computes f_n, proving the matching upper bound. For n=1 the output may select
the sole input at zero cost; for n=0 it may select a zero constant at zero
cost. In particular, for every n >= 1,

    C_B2(f_n) = n-1.

For any epsilon > 0, choose an integer N_epsilon >= 2 with
N_epsilon >= 1/epsilon. Then for every integer n >= N_epsilon,

    C_B2(f_n) = n-1 >= (1-epsilon)n.

Consequently liminf over all integer n tending to infinity of C_B2(f_n)/n
equals 1. The declared full guarantee is g(n)=1*n^1*(log2 n)^0*
(log2 log2 n)^0 = n for sufficiently large n, and its joint finite
projection is the pair [1,1]. The upper bound also shows that this proof
cannot support a coefficient greater than 1 or an exponent greater than 1
for this family. No stronger lower bound or new research result is claimed.
