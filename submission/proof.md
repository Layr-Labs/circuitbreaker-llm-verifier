# Odd parity: an exact full-B2 bound from the ancestor multigraph

## 1. Statement and limits

For the language L of bit strings with odd Hamming weight, let f_n be its
characteristic function on exactly n input bits; f_0=0. We prove
C_B2(f_n)=n-1 for every integer n>=1. The declaration uses the eventual
lower bound n-1 and joint transport projection [1,1]. This is a standard
elementary baseline theorem, not a new circuit-complexity result, not an
improvement over the challenge's literature seed, and not a claim to a
record, payment, or resolution of the open challenge. The proof below is
self-contained; no prior participant submission or prior judge decision is
a premise. No finite computation or participant executable is required.

## 2. One uniform total NP language, including small lengths

Define the same algorithm for every finite string x: set b=0; scan its n
actual input bits once; replace b by b XOR x_i at step i; return b.
This defines f_0=0 and a total function at all lengths. Define V(x,w) to
reject unless w is empty and otherwise return this scan result. Use the
fixed polynomial witness bound W(n)=0. If x belongs to L, the empty
witness makes V accept; if some witness of length at most W(n) makes V
accept, it is empty and x has odd weight. Thus both NP directions hold.
The evaluator and verifier are uniform and take O(n+1) bit operations.
There is no per-length advice, truth-table oracle, nonconstructive choice,
extra input padding, or counted witness bit.

## 3. Circuit model and an essential-variable observation

Use precisely the challenge's deterministic acyclic B2 model: one selected
output; each retained gate has two input pins and an arbitrary one of the
16 binary Boolean functions; inputs and wires cost zero; gates cost one;
fan-out and shared descendants are unrestricted. Repeated pins, constant
nodes, unary projections, unused gate arguments, and trivial truth tables
are allowed. We do not assume that a formula lower bound transfers to a
circuit, and we do not assume direct-sum additivity.

For n>=1, each variable x_i is essential to parity: fix all other bits,
flip x_i, and the output flips. Any exact circuit must therefore have a
structural directed path from input x_i to its selected output node.
Indeed, by induction in a topological ordering, any node's value depends
only on variable sources that have structural paths to it. Constant
sources have no variable dependence. If no path from x_i existed, the
output would be unchanged under every flip of x_i, a contradiction.
An output inversion does not change this observation; it only complements
the selected node's value. The observation is about necessary paths, not
an assumption that every structural path carries effective dependence.

## 4. Ancestor multigraph lemma, with sharing and degeneracies

Take an arbitrary exact circuit for n>=2. Its selected output cannot be
a variable or a constant (or an inversion of either): that would depend
on at most one variable, whereas f_n depends on all n variables. Thus its
selected output is a gate v. Retain all structural ancestors of v, v
itself, and all input/constant sources feeding those ancestors. Let g be
the number of gates in this retained subgraph, and c the number of
retained constant-source nodes. All n variable-source nodes occur by
Section 3. We may represent each original input x_i by one source node:
if a circuit drawing duplicates the same free input label, identifying
those source copies preserves the directed connections and computation.
Constants may be merged or kept separate; either convention only changes
c and does not weaken the inequality below.

Forget edge directions but preserve one edge for each input pin feeding
each retained gate. This yields a connected undirected multigraph: every
retained vertex has a directed path to v, hence an undirected path to v.
All pins feeding a retained gate come from retained ancestors/sources,
so the edge count is E=2g. Repeated pins give parallel edges; fan-out and
shared descendants give ordinary additional edges. The vertex count is
V=n+g+c. Acyclicity excludes a gate feeding itself or a directed cycle;
we do not assume the undirected graph is a tree.

Every finite connected undirected multigraph has E>=V-1. For a direct
proof, start with any one vertex; while not all vertices have been
reached, connectedness supplies an edge from the reached set to a new
vertex. Keep that edge and add the new vertex. After V-1 additions,
V-1 distinct edges have been kept. Additional or parallel edges cannot
reduce the edge count. Applying this fact gives

    2g >= n+g+c-1, hence g >= n+c-1 >= n-1.

The original circuit has at least g gates. This proves its gate count
is at least n-1, including arbitrary sharing and fan-out. Keeping a
constant/unary/trivial gate charges it just as required. If a constant
or unary output normalization removes a gate, its normalized exact
circuit still has the essential-variable property and the same argument
applies to that circuit. Negations absorbed in gate tables add no extra
basis operation or uncharged gate. A free terminal complement, if used,
leaves all ancestor vertices and edges unchanged. No degenerate or
terminal-output convention permits fewer than n-1 gates.

## 5. Matching upper construction and n=1

For n>=2 set y_2=x_1 XOR x_2 and y_k=y_(k-1) XOR x_k for k=3,...,n.
Each XOR is one allowed binary truth table; this chain uses exactly n-1
gates and computes parity by induction. For n=1 select x_1 as output
with zero gates, so C_B2(f_1)=0=n-1. The n=0 definition is only needed
for the total language and is outside the claimed exact n>=1 theorem.
Thus Section 4 and the construction prove the complete exact statement
for every claimed integer length, not merely infinitely many lengths.

## 6. Full eventual declaration and finite projection

Choose the single threshold N=1. For every integer n>=N, the exact
result implies C_B2(f_n)>=n-1. The leading monomial is n and the explicit
remainder is the constant monomial 1, strictly lower order. For every
real epsilon>0 choose an integer N_epsilon>=max(1,ceil(1/epsilon));
then for all n>=N_epsilon, (n-1)/n>=1-epsilon. In fact the exact ratio
is 1-1/n, so its all-integer-length limit and liminf are 1. The same
fixed language and evaluator supply both transport coordinates c=1 and
p=1. Neither a sparse subsequence, different family, asymptotic notation
alone, nor an accepted prior result is used to justify the explicit loss.

## 7. Provenance and execution scope

The dependency-on-all-inputs/connected-ancestor edge count and the XOR
chain are standard elementary observations and have been completely
reproduced here. This candidate was prepared from the protected model
and a fresh argument; it does not use the other account's MOD3 candidate,
proof, statement, cost files, or reviews. The parent task has historical
context about earlier work; no claim is made that that context was erased
or that elementary parity was discovered here for the first time.
No previous score, review, acceptance, or account identity is mathematical
evidence. The claim is analytic and includes zero experiments. Local
format/hash checks and advisory critique are separate from the official
cloud judgment, publication, human approval, and any prize eligibility.
