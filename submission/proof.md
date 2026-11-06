# Exact-weight class: strengthened finite lower bound

## Joint theorem and model
For all n>=0 let f_n(x)=1 iff the Hamming weight of x is floor(n/2),
including f_0(empty)=1. One fixed polynomial-time procedure counts the
input bits and compares the count to floor(n/2). A fixed verifier accepts
exactly the empty witness and the inputs satisfying this test. Therefore the
family is in NP in both directions, with witness bound zero, while the
circuit receives only x. We prove C_B2(f_n)>=2n-2 for every n>=3. The
joint score is coefficient 2 and exponent 1 with explicit constant loss 2.

The full B2 model has a single output, acyclic circuits, all 16 binary
Boolean tables, arbitrary sharing and depth, free wires, constant and input
nodes, and one unit of cost per retained gate. Repeated pins and degenerate
gates are allowed. The normalization proof below handles them.

## Normalization
A gate whose output is a constant, a parent literal or a complemented parent
literal can be removed while preserving the represented function whenever that
function depends on at least two distinct input variables. Constants and
uncomplemented literals are rewired directly. For a complemented literal,
rewire to its uncomplemented parent and complement the corresponding input
of every child truth table. This works with shared descendants and repeated pins:
transform every occurrence. If the removed gate is the output and its parent is
another gate, complement the parent's truth table, and complement the input
positions of all other children using that parent, preserving their values; the
parent itself becomes the output. If the output parent is a primary input, the
whole function would depend on at most one input, contrary to the hypothesis.
Constant outputs contradict the same hypothesis. A gate depending on at most
one distinct nonconstant parent has one of the four unary tables, hence is
removable by this rule. Each removal reduces the gate count by one. A surviving
previously degenerate gate remains degenerate: its new truth table is obtained
by complementing some inputs, and at most one effective parent remains. Thus
any k initially distinct degenerate gates can be removed in at least k steps;
if another removal already deletes one, it has already been counted.



## Four-gate base for the three-bit exact-one function
Write E(x,y,z) for exact-one on three bits. Every variable is essential.
Fixing any one variable to zero gives XOR of the other two, with two 1s in
its truth table; fixing it to one gives NOR of the other two, with one 1.
Also, fixing any two variables gives the three pairwise distinct residuals
identity, negation, and constant zero (the pair weights 0, 1, and 2).

Consider a circuit with at most three gates. Delete irrelevant gates and
normalize degenerate gates as above. A direct input cannot be the output.
If two inputs of the first gate each occur nowhere else, all restrictions of
those inputs affect the output through a single bit, giving at most two
residual functions, a contradiction. This also excludes circuits with at
most two gates: if the second gate uses the first, its only other pin must
be the third variable because the target depends on all three. Both parents
of the first gate then occur nowhere else. If the second gate does not use
the first, a single-output circuit cannot use both gates. Thus a
three-gate candidate must use all three gates.

Let the last gate be H. First suppose one input of H is a primary variable
x, and its other input is the output of a two-gate subcircuit G. If H uses
no two-gate subcircuit, at most two gates are relevant, already excluded.
For each fixed x=b, H is a unary function of G_b. Since both E_0=XOR(y,z)
and E_1=NOR(y,z) are nonconstant, that unary function is identity or
negation; G_0 must therefore be XOR or XNOR, and G_1 must be NOR or OR.

Write G as a binary gate of a first gate A and a primary input v; a direct
input, repeated A, or constant as the other pin cannot yield an essential
three-variable subcircuit. If v=x, A must depend only on y,z to make both
essential, and G_0 and G_1 are unary functions of the same A. Nonconstant
ones must be equal or complementary, but XOR/XNOR cannot equal NOR/OR.
If v is, say, z, then A must use y. If A does not use x, both G cofactors
are identical, again impossible. If A uses x and y, then for each b the
unary map A(b,y) must be nonconstant because G_b depends on y. It is hence
y or not-y. As z enters directly, the two truth tables G_0 and G_1 differ
only by a possible negation of the y input, so they have the same number of
ones. XOR and XNOR have two ones; NOR and OR have one or three. Contradiction.
This covers permutations of x,y,z and all 16 truth tables.

Now suppose H takes the outputs A and B of the first two gates. If B uses A,
its other pin must be a primary input v. The other two primary variables can
affect the output only through A, giving at most two residual functions when
that pair is fixed, contradicting the three residuals above. If B does not
use A, each of A and B is a gate on two primary inputs. All three variables
must be covered, so their input pairs overlap in exactly one variable x;
after renaming let A=A(x,y) and B=B(x,z). For each x=b, both E_b depend on
y and z, so A(b,y) and B(b,z) must both be nonconstant unary maps, hence
identity or negation. The two tables E_0 and E_1 are therefore obtained
from the same truth table H by independently permuting its rows and columns,
and have equal numbers of ones. But E_0 has two and E_1 has one, a
contradiction. The output cannot use constants, repeated pins or unary gates
in a normalized minimum circuit. Thus C_B2(F_3)>=4.


## Stronger class theorem
For integers n>=3 and 1<=k<=n-1, E(n,k) is the function which is one iff
the n-bit input has weight k. Every input is essential: choose k-1 ones
among the other n-1 inputs and flip the selected input.

At n=3 the only nontrivial exact-weight functions are E(3,1) and E(3,2).
The preceding section proves C_B2(E(3,1))>=4. E(3,2)(x) equals
E(3,1)(not x_1,not x_2,not x_3). Complementing every primary input does
not change B2 circuit complexity for a function depending on at least two
inputs: at each input pin of each consuming gate replace its truth table
by the table with that argument flipped, simultaneously for repeated pins.
If the output is a primary input or a constant this transformation needs
special handling, but neither is possible for these essential three-input
functions. Repeated fan-out merely applies the same replacement at each
consumer. The transformation is reversible and adds no gates. Thus the
same four-gate lower bound holds for E(3,2).

For n>=4, fix any two distinct input bits. The assignments 00,01,10,11
induce E(n-2,k), E(n-2,k-1), E(n-2,k-1), E(n-2,k-2), treating an
out-of-range weight as the constant-zero function. The three functions
with consecutive weights k,k-1,k-2 are distinct: at least two weights
lie in [0,n-2], because 1<=k<=n-1. Each in-range weight has an input
accepted only by its corresponding function, and the remaining function
is either a different in-range weight or the zero function.

Let S(n) be the minimum gate count over 1<=k<=n-1. Take a minimum circuit
for a minimizing k, trim unused gates and normalize it. Every input is
essential, so the first topological gate has two distinct primary inputs.
If each occurs only as an input of that gate, the other gates can observe
the pair only through its one output bit, allowing at most two distinct
residual functions. This contradicts the three distinct pair restrictions.
The output cannot be a primary input because every one of n>=4 variables
is essential. Hence one of the two parents also feeds a distinct gate.

Choose its fixed value b so that the restricted weight stays nontrivial:
take b=1 if k=n-1, and b=0 otherwise. The remaining function is
E(n-1,k-b) with 1<=k-b<=n-2, so it still depends on every input. Each
of the two distinct gates now has a constant pin; normalization removes
both without adding gates, including for repeated pins, sharing, and
terminal negations. Thus S(n)>=S(n-1)+2. Since S(3)>=4, induction yields
S(n)>=2n-2 for every n>=3 and all 1<=k<=n-1.

The middle-weight function is E(n,floor(n/2)), which is in that class for
every n>=3. This establishes the full explicit lower bound. For each
epsilon>0 and n>=max(3,ceil(2/epsilon)), it implies
C_B2(f_n)>=(2-epsilon)n, proving the all-integer-length liminf claim.

## Provenance and limits
This work combines and generalizes our exact-one and exact-middle-weight
proofs, using the former's four-gate base and the latter's choice of fixed
input in the induction. It also uses the MOD3 normalization template and
its public predecessor. These are declared attribution dependencies.
The complete proof is repeated here: no prior model judgment is assumed.
The improvement over our earlier exact-middle-weight theorem is in its
explicit constant loss, from 3 to 2; the asymptotic pair remains (2,1),
which is below the challenge frontier. No novelty or reward is guaranteed.
