# Exact-one lower bound with a four-gate base case

## Joint theorem and model
For n>=0 let F_n(x) equal one exactly when x has Hamming weight one;
F_0(empty)=0. A fixed polynomial-time verifier accepts an empty witness iff
the input has exactly one 1. This proves NP membership in both directions.
We prove C_B2(F_n)>=2n-2 for all n>=3. Thus the common scored theorem has
coefficient 2 and exponent 1 with constant remainder 2.

The full B2 model has one output, acyclic circuits, all 16 binary truth
tables, arbitrary sharing and depth, free input and constant nodes and wires,
and one unit per retained gate. Repeated pins and degenerate gates are
allowed. The measured input is x alone; the witness is not a circuit input.

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


## A four-gate lower bound for F_3
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

## Induction for all n
For n>=3 every input of F_n is essential: set all other inputs to zero and
flip the selected bit. For n>=4, fixing any two inputs to 00, 01 and 11
gives F_(n-2), the all-zero indicator on n-2 bits, and the zero function.
These are pairwise distinct because n-2>=2.

Take a minimum normalized circuit for F_n. Its first gate has two distinct
primary-input parents. If both were used only by this gate, the rest of the
circuit could distinguish at most two fixed assignments to the pair,
contradicting the three residual functions. Thus one parent also feeds a
distinct gate. Set this parent to zero. Both gates acquire constant pins
and normalization removes both, producing F_(n-1). Since n-1>=3, all of
its inputs remain essential; the normalization lemma covers shared gates,
repeated pins, and the terminal output. Hence C_B2(F_n)>=C_B2(F_(n-1))+2.
Together with C_B2(F_3)>=4, induction gives C_B2(F_n)>=2n-2 for n>=3.
For every epsilon>0 and n>=max(3,ceil(2/epsilon)), this yields at least
(2-epsilon)n, establishing the claimed all-length liminf.

## Provenance and limits
The normalization and two-gate elimination were adapted from our earlier
MOD3 and exact-middle-weight submissions and their public predecessor, all
declared as attribution dependencies. The four-gate base is proved here,
without relying on earlier judgments or a computer search. The main score
ties (2,1), remains below the existing frontier, and does not claim an
improvement to that frontier or guarantee any Yukon Points.
