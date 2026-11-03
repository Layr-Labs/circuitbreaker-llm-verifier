# A self-contained 2n-4 bound for MOD-3 in full B2

## 1. One fixed total family and the theorem

For every n>=0 let f_n(x)=1 exactly when sum_i x_i=0 modulo 3.
In particular f_0(empty)=1. We prove C_B2(f_n)>=2n-4 for every integer
n>=3. The declared asymptotic guarantee is 2n-o(n), with finite projection
(2,1). This is a conservative elementary bound, not a record or an
improvement on the organizer's 31/10 baseline.

Let L be all finite binary strings with Hamming weight divisible by three.
One fixed deterministic algorithm evaluates L by scanning x, maintaining
its residue in {0,1,2}, initially zero. Its runtime is O(n+1). Set W(n)=0
and V(x,w)=1 iff w is empty and this scan accepts. If x is in L, the empty
witness is accepted; conversely acceptance of a permitted witness means
that the scan found residue zero, so x is in L. Thus L is in P and NP.
The circuit inputs are exactly the n bits of x: there are no witness bits,
padding, advice, field encodings, or selected subsequences of lengths.

## 2. Normalization with explicit transfer to the charged B2 model

For the lower-bound proof only, use an auxiliary class of signed circuits.
Its charged nodes are arbitrary binary Boolean tables. Constants, wire
complements, and an optional final output complement are free. This is
an enlargement of the competitor class, not an assertion that a retained
negated-projection gate is free in the original B2 model.

Every original B2 circuit with g gates maps to a signed representation
with at most g charged nodes. Process gates in topological order: replace
constant gates by constants, bypass unary projections, carry negated
projections as signs on wires, and absorb incoming signs in the truth
table of each remaining binary node. Repeated pins on the same underlying
signal, with either signs, restrict a binary table to a unary function
and are bypassed likewise. Keep any output complement explicitly. Nothing
is duplicated. A negated-projection gate charged in the original circuit
may disappear in this relaxation; ignoring that cost only helps the
competitor. Therefore a lower bound for EVERY signed representation is
a lower bound for every original full B2 circuit. Only this direction of
model transfer is used.

Trim nodes outside the output's ancestry. In a minimum signed circuit,
every retained node has two distinct nonconstant predecessor nodes and
a truth table depending on both arguments; otherwise it can be removed
by the preceding operations. Absorb incoming signs into its table. After
restrictions, repeat the same simplifications; each removes a charged node,
so normalization terminates and never increases the count. A terminal
unary, constant, or negated-projection node is also removable, with its
output represented by a constant or signed wire. Thus terminal effects
are explicitly covered in the auxiliary model and no uncharged gate
is smuggled into the original B2 conclusion.

We allow arbitrary depth, fan-out, and shared descendants. A distinct
successor below always means a distinct charged DAG node, not a pin.

## 3. A restriction-closed class used for induction

For m>=0, weights a_i in {1,2}, and t in {0,1,2}, define
F_(a,t)(x)=[sum_i a_i*x_i=t modulo 3]. The fixed scored family is the
special case a_i=1,t=0. The other members are auxiliary induction cases,
not alternative scored languages or per-length advice.

Lemma A (residue coverage). Any at least two nonzero weighted bits can
realize all three residues. For two weights a,b, the four sums are
{0,a,b,a+b}. If a=b, this contains {0,a,2a}, all three residues. If a!=b,
it already contains {0,1,2}. Fix all further bits to zero.

Lemma B (essential inputs). If m>=3, every variable of F_(a,t) is
essential. For variable i the other m-1 bits realize t by Lemma A.
Their sum t makes output one at x_i=0 and zero at x_i=1, since a_i!=0.

Lemma C (no isolated pair compression). Suppose m>=4. For any pair i,j,
its four assignments realize all three partial sums s=a_i*x_i+a_j*x_j.
The resulting cofactors on the remaining m-2 bits are
[sum_(k!=i,j) a_k*x_k=t-s], for s=0,1,2. They are pairwise distinct,
since the remaining bits realize every residue by Lemma A. If the pair's
only influence on the output passed through one Boolean signal h(x_i,x_j),
there could be at most two cofactor functions, one per value of h. This
contradicts the three distinct cofactors. No assumption about the
structure or sharing of descendants is involved.

## 4. Classification of the binary gate tables

There are two constant tables and four essentially unary tables. Of the
six tables with two ones, four are the unary projections or negated
projections and the other two are XOR and XNOR. The eight tables with one
or three ones have the form
((z xor alpha) AND (w xor beta)) xor gamma, for bits alpha,beta,gamma.
This follows by placing the unique one of an AND table at an arbitrary
input pair and optionally complementing the output. These counts exhaust
all 16 tables. Normalization therefore leaves only XOR/XNOR or these eight
AND-type tables. Each AND-type table has a canalizing value for EACH
argument: set its corresponding literal to zero, forcing output gamma.

## 5. Two distinct nodes can be eliminated per lost variable

Let m>=4 and choose a minimum normalized signed circuit for F_(a,t),
with g charged nodes. Choose a first retained gate G in topological order.
Its two distinct predecessors must be original input nodes x_i,x_j.
It cannot be the output: Lemma B says the function depends essentially
on all m inputs, but G sees only two. All incoming signs have been absorbed.

Case 1: G is AND-type. If both x_i and x_j had no charged successors other
than G, their only influence on the output would pass through G, violating
Lemma C. They cannot themselves be the output, by essentiality of all m
inputs. Hence one, say x_i, has another distinct charged successor H.
Set x_i to its canalizing value b for G. Gate G becomes constant and is
removed. Gate H loses its x_i predecessor and becomes a function of at
most one remaining signal; it is removed as a constant or signed
projection. If removing G already removes H, that still deletes the two
distinct original nodes G,H. They may share other predecessors or be
ancestrally related; neither event prevents these removals. No node is
added. The restricted function has m-1 inputs, their old nonzero weights,
and target t-a_i*b modulo 3. It remains in the induction class.

Case 2: G is XOR or XNOR, so G=x_i xor x_j xor e for some bit e.
If a_i=a_j, impose x_i=x_j. Then G=e and the weighted sum replaces the
pair by 2*a_j*x_j, a nonzero weight modulo 3; the target is unchanged.
If a_i!=a_j, impose x_i=1-x_j. Then G=1 xor e, and the pair's contribution
is a_i+(a_j-a_i)*x_j. The new weight a_j-a_i is nonzero modulo 3, and the
new target is t-a_i. In both cases exactly m-1 independent Boolean inputs
remain with nonzero weights. They range over ALL assignments, not a
promise on the new input space. Replace every reference to x_i with
the chosen (possibly complemented) reference to x_j. No cycle is created
and no charged node is duplicated or added.

Since G is a retained output ancestor but not the output, it has a charged
successor H. Its becoming constant leaves H with at most one nonconstant
predecessor, so H is removed, including if H is terminal. Thus G and H
are two distinct deleted nodes. If simultaneous input identification
simplifies H earlier, H is still a second deleted node. Other fan-out
and arbitrary shared subgraphs do not change this counting.

In either case there is an (m-1)-input member of the auxiliary class
computed by a signed circuit with at most g-2 charged nodes.

## 6. Induction, termination, and all lengths

Every member of the auxiliary class on m>=3 variables needs at least
2m-4 charged signed nodes. For m=3 there are three essential inputs by
Lemma B. Zero charged nodes can depend on at most one input; one binary
node can depend on at most two. Constants and output signs do not increase
essential arity. Thus at least two nodes are required, exactly 2m-4.

For m>=4, Section 5 produces a class member on m-1 variables with at most
g-2 nodes. The induction hypothesis gives g-2>=2(m-1)-4, or g>=2m-4.
The integer m decreases by one each step to the base m=3, so every integer
length is covered. Minimum circuits exist because finite circuits exist
and their nonnegative integer sizes have a minimum. No bound on depth
or fan-out, no direct-sum assumption, and no restriction to formulas,
monotone gates, or a smaller truth-table basis was imposed.

By the model transfer in Section 2 this holds for the original B2 gate
count, in particular C_B2(f_n)>=2n-4 for every integer n>=3. The total
family is already defined at n=0,1,2; no bound at those lengths is needed
for the asymptotic theorem.

## 7. Exact guarantee and quantifiers

For every epsilon>0 choose N=max(3,ceil(2/epsilon)). For every integer
n>=N the proved finite bound implies C_B2(f_n)>=2n-4>=(1-epsilon)*2n.
Consequently liminf C_B2(f_n)/(2n)>=1. The canonical polynomial-log
expression has coefficient 2, exponent 1, and both logarithmic exponents
zero. Equivalently, for every eta>0 and every integer
n>=max(3,ceil(4/eta)), C_B2(f_n)>=(2-eta)*n. Both score coordinates
refer to this same theorem for one fixed NP language.

## 8. Provenance and limits

The organizer's baseline was read before editing. It uses a different
family, affine disperser constructions, and the Li-Yang 31/10 analysis.
Those construction and case lemmas are not premises of this proof.
MOD-3 and gate elimination are classical circuit-complexity topics; no
novelty, priority, optimality, record improvement, or P-versus-NP separation
is claimed. This submission supplies every fact its theorem needs,
without relying on a previous submission or judgment. The optional LP
certificate is not used. There is no computational premise: small local
sanity checks do not replace the all-circuit argument above.
