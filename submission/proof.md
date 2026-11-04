# A self-contained 2n-3 gate lower bound for MOD3

## Statement and uniform family
For n>=0, define M_a^n(x)=1 iff the Hamming weight of x is congruent to a
modulo 3, for a in {0,1,2}. Our single family is f_n=M_0^n, including
f_0(empty)=1. A deterministic algorithm maintains the count modulo 3 while
reading x and outputs whether the count is zero. It runs in O(n) time.
A fixed verifier accepts iff the witness is empty and this algorithm outputs 1;
the witness bound is zero. Thus the yes-language is in P and NP, with both
implications of the witness condition immediate. The circuit inputs are x only.

We prove that for n>=2, every full B2 circuit for f_n has at least 2n-3
gates. The model allows arbitrary acyclic sharing and depth, all 16 two-input
Boolean tables, repeated pins, constants and unary gates, with every gate counted;
inputs, wires, constant nodes, and output selection cost zero.

## Removing degenerate gates
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

## Properties of residue functions
For each n>=2 and residue a, M_a^n depends on every variable. To see this,
fix the other n-1 variables to have weight t congruent to a-1 modulo 3;
when a=1 choose t=0, when a=2 choose t=1, and when a=0 choose t=2 if
n>=3. For n=2,a=0 instead choose t=0: weights 0 and 1 give distinct outputs.
All chosen t are achievable and changing the remaining bit flips the output.
Also, for m>=1 the three functions M_0^m,M_1^m,M_2^m are pairwise distinct:
inputs of weights 0 and 1 separate their three characteristic functions
(the third is zero on both).

Fix n>=3 and two distinct variables xi,xj. Setting these variables to 00,
01,10,11 respectively gives M_a^(n-2), M_(a-1)^(n-2),
M_(a-1)^(n-2), M_(a-2)^(n-2), with subscripts modulo 3. Since n-2>=1,
these include three distinct functions of the remaining variables.

## Gate-elimination induction
Let S(n)=min over a in {0,1,2} of the minimum number of B2 gates computing
M_a^n. For n>=3 take a minimum circuit for any M_a^n and remove degenerate
gates until all remaining gates depend on two distinct nonconstant parents.
Every variable is essential. Select a topologically first gate G; its parents
must be two distinct primary inputs xi,xj (constant nodes, repeated pins and
unary gates have been removed).

If each of xi and xj feeds only G, then the rest of the circuit can access
these two variables only through the single bit G(xi,xj). Fixing xi,xj
could then produce at most two distinct residual functions of the other inputs.
This contradicts the three distinct residual functions proved above. Direct
output of one of these inputs cannot occur because M_a^n depends on all n>=3
variables. Hence at least one of xi,xj feeds some distinct second gate H.
Fix that input to 0. Both G and H now have a constant pin, so both are
initially degenerate and at least two gates can be removed. The restricted
function is another residue function on n-1 variables, all of which are
essential because n-1>=2. The removal lemma applies throughout. Therefore
S(n)>=S(n-1)+2. In the remaining case in which a minimal circuit has already
undergone deletions, the resulting bound only becomes stronger.

For n=2 each residue function depends on both inputs and cannot be computed
with zero gates, so S(2)>=1. Induction gives S(n)>=2n-3 for every n>=2.
In particular C_B2(f_n)>=2n-3. For any epsilon>0, if n>=max(2,ceil(3/epsilon)),
then C_B2(f_n)>=(2-epsilon)n. Thus liminf over all integer lengths of
C_B2(f_n)/n is at least 2. The two score coordinates are the same theorem's
coefficient 2 and exponent 1; no upper bound or improvement over the 4n
frontier is asserted.

## Attribution and scope
The modulus-three 2n-O(1) lower bound is classical (Schnorr 1974). The author
consulted the public note of Yukon submission 466a2973-4bf2-4e79-be86-b9c6a87260d5
by kaankolcu, including its gate-elimination outline; that submission is declared
as a dependency. This proof provides the argument above in full and does not
rely on a previous judgment. It makes no originality or record claim.
