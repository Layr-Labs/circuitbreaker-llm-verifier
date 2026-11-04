# MOD5 lower bound and a modulo-m class theorem

## The jointly scored theorem
For n>=0, f_n(x) is 1 if and only if the Hamming weight of x is divisible by
5; in particular f_0(empty)=1. One deterministic algorithm scans x while
maintaining the count modulo 5 and accepts exactly residue zero; it runs in
O(n) time. A fixed verifier accepts iff its witness is empty and this algorithm
accepts. This proves both directions of NP membership with witness bound zero.
The circuit receives only the input x, not the witness. We prove below that
C_B2(f_n)>=2n-9 for every n>=5. Hence for each epsilon>0 this is at least
(2-epsilon)n when n>=max(5,ceil(9/epsilon)), so the all-length liminf is at
least 2. The two native scores are coefficient 2 and exponent 1 for this same
single family and theorem.

The full B2 model has acyclic circuits, one output, all 16 binary truth tables,
free wires, input and constant nodes, arbitrary fan-out and depth, and a charge
of one for every retained gate including unary and constant gates. The argument
applies to repeated pins and shared descendants.

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


## Supporting class theorem for every fixed modulus
Fix m>=3. For each residue r in {0,...,m-1}, let M^n_r(x) be one exactly when
the input weight equals r modulo m. We show C_B2(M^n_r)>=2n-2m+1 for n>=m.

First, if n>=m, every variable is essential in every M^n_r. Fix the other
n-1 bits to have weight t in {0,...,m-1} congruent to r-1 modulo m. This is
possible since n-1>=m-1; changing the remaining input from zero to one changes
the value from zero to one. Second, for every l>=m-1, the m functions M^l_r
are pairwise distinct: for each r an input of weight exactly r exists and is
accepted by M^l_r but by no other residue function.

For n>=m+1 fix two distinct variables. Their settings 00,01,10,11 give the
residual functions M^(n-2)_r, M^(n-2)_(r-1), M^(n-2)_(r-1),
M^(n-2)_(r-2), respectively, with residues modulo m. The three residue indices
are distinct because m>=3, and the functions are distinct because n-2>=m-1.

Let S_m(n) be the minimum gate count over all residues r. Consider a minimum
circuit for a minimizing residue at n>=m+1. Normalize its degenerate gates as
above. Every input is essential, and the first remaining topological gate has
two distinct input variables x,y. If each variable feeds only this gate, the
rest of the circuit depends on x,y only through its one output bit and hence
the four restrictions to x,y have at most two distinct residual functions,
contrary to the preceding paragraph. The selected output cannot be a bare input
because the function depends on all n>=4 inputs. Therefore one of x,y feeds a
second distinct gate. Fix that variable to 0: both gates now have constant pins
and are degenerate. Removing them yields a circuit with at least two fewer gates
for some M^(n-1)_s. Every variable remains essential since n-1>=m, so the
normalization lemma applies throughout. Initial deletions only strengthen the
bound. Thus S_m(n)>=S_m(n-1)+2 for n>=m+1.

For n=m every residue function depends on every input and so needs at least
one gate. Therefore S_m(m)>=1, and induction yields S_m(n)>=2(n-m)+1 for
every n>=m. Taking m=5,r=0 gives the jointly scored MOD5 bound 2n-9.
The argument proves only these lower bounds, not their optimality and not a
new 4n frontier result.

## Provenance and limitations
The gate-elimination technique and normalization were adapted from our prior
MOD3 submission 43c5fb3c-7bae-424c-a0ad-3dcbc72c82a6 and its predecessor
466a2973-4bf2-4e79-be86-b9c6a87260d5, which are explicitly declared as
attribution dependencies. The MOD3 2n-O(1) bound is classical (Schnorr 1974).
The class extension is a direct elementary adaptation; no novelty claim is made.
No previously obtained judgment is used as a premise. The proof is fully given
above, without experiments or reliance on external mathematical results.
This scored bound has the same (2,1) coordinates as our prior MOD3 submission
and is dominated by a (4,1) result; it does not guarantee Yukon Points.
