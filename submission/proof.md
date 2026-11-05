# Exact-middle-weight: an explicit 2n-3 gate lower bound

## Family, model, and main theorem
For each n>=0 let f_n(x)=1 exactly when the Hamming weight of x equals
floor(n/2). At n=0 the empty string has weight zero and f_0(empty)=1.
A fixed deterministic algorithm counts the 1 bits, computes floor(n/2),
and compares the numbers in polynomial time. A fixed verifier accepts iff
its witness is empty and this algorithm accepts, proving both directions of
NP membership with witness bound zero. The circuit receives only x.

We prove C_B2(f_n)>=2n-3 for every integer n>=2. The full B2 model has
one output, acyclic circuits, arbitrary depth and sharing, all 16 two-input
Boolean gate tables, free primary inputs, wires and constant nodes, and one
unit for every retained gate. Repeated pins and degenerate gates are allowed.
The explicit v2 remainder is the constant 3; the pair (2,1) refers to this
same family and theorem, and no superlinear result is asserted.

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


## Class theorem: every nontrivial exact-weight function
For each n>=2 and 1<=k<=n-1, define E(n,k)(x)=1 iff |x|=k.
Every variable is essential: fix the other n-1 inputs to have weight k-1;
flipping the variable changes the total from k-1 to k.

For n>=3, fix two distinct variables. The assignments 00,01,10,11 yield,
on the remaining n-2 inputs, respectively E(n-2,k), E(n-2,k-1),
E(n-2,k-1), E(n-2,k-2). When a target weight is outside 0,...,n-2,
the corresponding function is identically zero. The three functions at
weights k,k-1,k-2 are pairwise distinct: at least two of these consecutive
integers lie in 0,...,n-2 because 1<=k<=n-1 and n-2>=1; each represented
weight has an assignment accepted only by its own exact-weight function.
At most one of the three can be outside the range.

Let S(n) be the minimum B2 gate count over k in 1,...,n-1. Choose an
optimal circuit for a minimizing E(n,k), n>=3. Normalize degenerate gates.
Every input is essential, so a first topological gate G has two distinct
primary input parents x,y. If both parents fed only G, all four fixed
assignments to x,y would affect the rest of the circuit through only the
one-bit value of G and would yield at most two residual functions. This
contradicts the three pairwise distinct residual functions above. Direct
output of x or y is impossible because E(n,k) depends on all n>=3 inputs.
Thus one of x,y also feeds a distinct gate H.

Choose the value b of this input so that the restricted target weight stays
nontrivial: set b=1 when k=n-1, and otherwise set b=0. Then the restricted
function is E(n-1,k-b), with 1<=k-b<=n-2. Both G and H acquire a constant
input and can be removed by normalization. The restricted function still
depends on all n-1>=2 inputs, so the normalization lemma applies including
the terminal-output case. At least two distinct gates are removed; any
pre-normalization deletions only strengthen this bound. Therefore
S(n)>=S(n-1)+2. At n=2, E(2,1) is XOR and depends on both inputs, hence
requires at least one gate. Induction proves S(n)>=2n-3 for all n>=2.

For every n>=2, floor(n/2) belongs to 1,...,n-1, so f_n is a member of
the class and the same bound holds for every n>=2. This proves the main
explicit-remainder theorem and the separate unscored class lemma. In
particular, for every epsilon>0 the bound also implies
C_B2(f_n)>=(2-epsilon)n for all n>=max(2,ceil(3/epsilon)).

## Provenance and limitations
This is a classical two-gate elimination technique, not a claim of a new
research result. The gate normalization and general proof template are adapted
from our MOD3 submission 43c5fb3c-7bae-424c-a0ad-3dcbc72c82a6 and the public
MOD3 note 466a2973-4bf2-4e79-be86-b9c6a87260d5; both are declared attribution
dependencies. The exact-weight class, its choice of restricted value and its
base and induction are established above and use no prior model judgment.
The main theorem is weaker than known 2.5n bounds for suitable symmetric
functions and the displayed 4n frontier; it makes no promotion or reward claim.
