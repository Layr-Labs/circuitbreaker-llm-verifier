# MOD5 lower bound and a modulo-m class theorem

## The jointly scored theorem
For n>=0, f_n(x) is 1 if and only if the Hamming weight of x is divisible by
5; in particular f_0(empty)=1. One deterministic algorithm scans x while
maintaining the count modulo 5 and accepts exactly residue zero; it runs in
O(n) time. A fixed verifier accepts iff its witness is empty and this algorithm
accepts. This proves both directions of NP membership with witness bound zero.
The circuit receives only the input x, not the witness. We prove below that
C_B2(f_n)>=2n-5 for every n>=4. Hence for each epsilon>0 this is at least
(2-epsilon)n when n>=max(4,ceil(5/epsilon)), so the all-length liminf is at
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
the input weight equals r modulo m. We show C_B2(M^n_r)>=2n-m for n>=m-1.

First, if n>=m-1, every variable is essential in every M^n_r. The
other n-1 bits realize all weights from 0 through n-1. When n=m-1 this
includes exactly 0,...,m-2. If r-1 modulo m is among these values, fix
the other bits to that weight and flipping the selected bit changes rejection
to acceptance. Otherwise r=0: fix the other bits to weight 0, so the value
changes from acceptance to rejection. The same argument works for larger n.

Second, if l>=m-2, any three different residue functions M^l_a, M^l_b,
M^l_c are pairwise distinct. Each residue represented by a weight in
0,...,min(l,m-1) has an input it accepts and the other functions reject.
At most one residue is unrepresented, so at most one of the functions is
identically zero. Thus all three are distinct, including at l=m-2.

For n>=m fix two distinct input bits. The assignments 00,01,10,11 leave
M^(n-2)_r, M^(n-2)_(r-1), M^(n-2)_(r-1), M^(n-2)_(r-2),
respectively, with indices modulo m. The three indices are distinct as
m>=3 and the functions are distinct because n-2>=m-2.

Let S_m(n) be the minimum gate count over all residues r. At n>=m take
an optimal circuit and normalize degenerate gates as above. Every input is
essential, and its first topological gate has two distinct primary inputs
x,y. If each fed only that gate, the remainder of the circuit would see
x,y only through a single bit, so at most two residual functions could
result from fixing x,y. This contradicts the three distinct residuals.
The output cannot be a bare input or constant since n>=m>=3 and all inputs
are essential. Hence at least one of x,y feeds a second distinct gate. Fix
it to zero. Both original gates acquire a constant pin and at least two
gates can be eliminated by the normalization lemma. The restriction is
another residue function on n-1 variables, which remains essential on all
its variables since n-1>=m-1>=2. So S_m(n)>=S_m(n-1)+2 for n>=m.

For the base n=m-1, all inputs are essential. Trim a circuit to the
ancestors of its output. Its undirected underlying graph is connected and
contains n input nodes, g gates and k explicit constant nodes. It has at
least n+g+k-1 distinct edges, but at most 2g edges enter the two input pins
of its gates. Repeated pins cannot increase that count. Thus g>=n+k-1>=m-2.
A direct output of an input or constant is impossible because the function
has at least two essential variables. Therefore S_m(m-1)>=m-2. Induction
yields S_m(n)>=m-2+2(n-m+1)=2n-m for all n>=m-1. For m=5,r=0 this gives
C_B2(MOD5_n)>=2n-5 for n>=4, proving both score-source claims.

## Provenance and limitations
The gate-elimination technique and normalization were adapted from our prior
MOD3 submission 43c5fb3c-7bae-424c-a0ad-3dcbc72c82a6 and its predecessor
466a2973-4bf2-4e79-be86-b9c6a87260d5, and our modulo-m submission
0bc96713-fadd-4bee-9eb8-0f47713519fb, and finite-base improvement
6fb8035a-1838-407e-babf-3ae0da985c66, declared as
attribution dependencies. The MOD3 2n-O(1) bound is classical (Schnorr 1974).
The class extension, boundary-residue refinement and connectivity base are elementary; no novelty claim is made.
No previously obtained judgment is used as a premise. The proof is fully given
above, without experiments or reliance on external mathematical results.
This scored bound has the same (2,1) coordinates as our prior MOD3 submission
and is dominated by a (4,1) result; it does not guarantee Yukon Points.
