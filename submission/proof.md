# A complete (7/3)n - 4 lower bound for MOD-3 in full B2

## 1. Fixed family, uniformity, and the precise theorem

For every integer n >= 0, define

    f_n(x_1,...,x_n) = [sum_i x_i = 0 modulo 3].

Thus f_0(empty)=1. The one scored language is the set of all finite binary
strings whose Hamming weight is divisible by three. We prove

    C_B2(f_n) >= (7/3)n - 4                  for every integer n >= 2.   (1)

In particular, the finite coefficient/exponent pair is (7/3,1), for the
same family and the actual input length n. This is an exposition of a
classical technique, not a new record or an improvement over the 31/10 seed.

One fixed algorithm starts with residue 0, reads the bits of x in order,
adds each bit modulo 3, and accepts exactly when the final residue is 0.
It takes O(n+1) time. Set W(n)=0. Define V(x,w) to accept if and only if w
is empty and the residue algorithm accepts x. Membership in the language
provides the empty accepted witness. Conversely, an accepted permitted
witness means precisely that x has weight divisible by three. This proves
both directions of the NP definition. In fact the family is in P.
The n circuit inputs are just x: there are no witnesses, padding, advice,
field encodings, or restrictions to a subsequence of lengths.

We use addition and multiplication modulo 2 only in Sections 5--6, where
they are explicitly distinguished from sums modulo 3.

## 2. An auxiliary circuit model and its direction of transfer

For the proof, enlarge the allowed class to signed circuits: constants,
wire signs (complements), and an optional output complement are free.
Each retained binary node is still a charged gate with any of the sixteen
binary Boolean tables. The graph is acyclic, with arbitrary fan-out and
one designated output. A wire, or a constant, may itself be the output.

Every circuit with g gates in the problem's full B2 model yields such a
signed circuit with at most g retained gates. Process the original DAG
in topological order. A constant node is replaced by its value. A unary
node is bypassed as a literal of its predecessor. If both pins now refer
to the same signal, with either choice of signs, the binary table reduces
to a unary function and can also be bypassed. Incoming signs are absorbed
into descendant tables, and an output sign is retained explicitly. Delete
nodes outside the output ancestry. These operations duplicate no gate.
In the original model a retained negated-projection node costs one; its
removal here only enlarges the competitors and is not a free-negation
claim about the original model. A lower bound for every signed circuit
therefore implies that same lower bound for the original charged model.

After this normalization, each gate depends on both its pins, its pins
are distinct nonconstant predecessor nodes, and it is an output ancestor.
There are exactly ten possible tables left:

* two parity tables u XOR v XOR e, called X-type;
* eight tables ((u XOR a) AND (v XOR b)) XOR c, called A-type.

Here a,b,c,e are bits. To check exhaustion, the two constant and four
essentially unary tables account for six of the sixteen tables. The
one-one and three-one tables are exactly the eight A-type tables. Of the
six two-one tables, four are unary and the other two are X-type.
Complementing either argument or the output preserves the X/A type of
each surviving gate. Each A-type gate has a value of either argument
that makes its output constant; we call such a value canalizing.

For a normalized signed circuit D, let X(D) and A(D) count its X-type
and A-type nodes, respectively. Write

    s(D)=X(D)+A(D),             mu(D)=3X(D)+2A(D).                     (2)

On restricting an input to a constant, bypass every resulting unary or
constant node and trim unused nodes again. These operations never add
a charged node. Every retained node has its old type, possibly with
changed signs, unless it becomes unary and is deleted. Thus mu cannot
increase, and deleting distinct original nodes subtracts at least the
sum of their old weights. A gate with a constant pin is deleted even if
its remaining value is a complement; terminal gates are covered by the
free output sign. Counting is by distinct DAG nodes, never by paths or
occurrences, so shared descendants and repeated pins cause no double
counting. If a node is already deleted earlier in a simplification, it
is still one of the distinct original nodes removed.

## 3. The auxiliary residue functions and a one-bit bottleneck

For t in {0,1,2}, put

    M_(n,t)(x)=[sum_i x_i = t modulo 3].

These functions support induction; the scored family remains t=0 at
every length. Fixing one variable to b gives M_(n-1,t-b), with the target
read modulo 3.

For n>=2 each M_(n,t) depends essentially on every input. Indeed, for
an input x_i, compare targets t and t-1 on the other n-1 bits. If n>=3,
those bits can realize all three residues, so the two target indicators
differ. If n=2, the other bit realizes 0 and 1; any two distinct targets
in {0,1,2} still give different indicators (the three are 1-x, x, and 0).
Consequently none of these n>=2 functions is constant or unary.

For n>=3, choose any two distinct inputs x_i,x_j. Their four assignments
have partial sums 0,1,1,2 modulo 3. The induced cofactors on the other
n-2 inputs are the three indicators of targets t,t-1,t-2. When n-2=1,
these indicators are 1-x, x, and 0 and are pairwise distinct. When n-2>=2,
the remaining bits realize all three residues, so again the indicators
are pairwise distinct.

It follows that the influence of this input pair on the output cannot
pass only through a single Boolean signal h(x_i,x_j). If it did, fixing
the pair would yield at most two cofactor functions, one for each value
of h. This observation makes no assumptions about the part of the DAG
after the signal, including its fan-out and all sharing there.

## 4. Weighted elimination, with all successor coincidences covered

We prove that every normalized signed circuit computing M_(n,t) obeys

    mu(D) >= 6n - 10                   (n>=2, t in {0,1,2}).          (3)

For n=2 all inputs are essential, so the circuit has at least one gate,
of weight at least 2. This is 6*2-10.

Let n>=3, and suppose the assertion is proved on n-1 variables. Take
any normalized signed circuit D for M_(n,t). It has a first gate G in
topological order, whose predecessors are two distinct original inputs,
say x and y. It cannot be the output, since all n>=3 inputs are essential.
If both x and y had G as their only gate successor, their whole influence
would pass through G, contradicting the three cofactors in Section 3.
Thus, after interchanging x and y if necessary, x feeds a further gate
P distinct from G. We consider every possible pair of gate types.

Case I: G and P are both X-type. Fix x to either constant. Both G and
P now have at most one nonconstant argument, so both are deleted. Their
old weights sum to 6. This remains true if P also reads G: if further
simplification occurs, it can only delete more nodes.

Case II: P is A-type. Set x to the canalizing value for its pin at P.
Then P is constant, irrespective of its other pin's value. Gate G is
also removed because its pin x is constant. The restricted output is
M_(n-1,t-b), which is nonconstant since n-1>=2. Hence P cannot be the
designated output. As an output ancestor, P has a gate successor H.
It is distinct from P by acyclicity, and from G because G precedes every
other gate and has only original inputs as pins. The constant value of
P makes H unary or constant, so H is removed as well. These are three
distinct gates G,P,H, each of old weight at least 2. Thus the weight loss
is at least 6. This covers P reading G, P sharing another ancestor with
G, and any number of additional successors of either node.

Case III: P is X-type and G is A-type. Set x to its canalizing value for
G. Then G is constant and P is removed because it reads the fixed input
x. The restricted output is again nonconstant. If G has a successor H
different from P, its constant pin makes H removable. G,P,H are distinct,
and their old weights sum to at least 6 (in fact at least 7).
The only remaining configuration is that the sole gate successor of G
is P. Such a successor must exist because G is a non-output ancestor.
In this configuration the two pins of P are x and G, since P is already
a successor of x. Both pins are now constants, so P also becomes
constant. It cannot be the output because the restricted function is
nonconstant. Take a successor H of P. By acyclicity H is different from
both G and P. It is removed on receiving the constant from P. Again
three distinct original nodes are deleted, with total old weight at
least 6. No assumption of independent or unshared successor subgraphs
has been made.

The cases are exhaustive: either both gates are X-type, P is A-type,
or only G is A-type. In each case one input is fixed and a circuit D'
for an (n-1)-variable residue function is left with

    mu(D') <= mu(D)-6.

Apply normalization if needed. By induction mu(D')>=6(n-1)-10, proving
(3). This induction terminates at n=2 and includes every integer n>=2.
We did not assume D was optimal or had bounded depth or bounded fan-out.

## 5. Why algebraic degree forces A-type gates

All algebra in this section is in the Boolean polynomial ring over F_2,
where x_i^2=x_i. Every Boolean function h has a unique multilinear
polynomial. Existence follows by summing, modulo 2, the products of
literals that indicate its accepting assignments and then expanding.
For uniqueness, a multilinear polynomial can be written P_0+x_n P_1.
If it vanishes on all Boolean inputs, evaluating at x_n=0 makes P_0
identically zero; evaluating at x_n=1 then makes P_1 identically zero.
Induction on the number of variables proves both polynomials are zero.
The largest monomial size is its algebraic degree. A constant has degree
at most zero. Substitution of affine forms cannot increase degree:
a degree-d monomial becomes a product of at most d degree-one forms,
and reduction with x_i^2=x_i can only decrease degree.

We prove directly that a circuit with a A-type gates has degree at most
a+1, regardless of the number of X-type gates, free constants, signs,
or shared descendants. Replace each A-type gate by one AND and affine
changes of its inputs and output. This uses exactly one multiplication
per such gate; every other operation is affine. We may therefore work
with arbitrary acyclic circuits of binary AND gates and affine operations
with free constants. The assertion is by induction on the number a of
AND gates, simultaneously over every finite number of inputs.

For a=0 the computed function is affine, of degree at most 1. If a>0,
consider the first AND gate in topological order. Its two input functions
are affine. If one is constant, this AND is affine as well and can be
removed without adding an AND; the induction hypothesis applies.
Otherwise choose a nonconstant affine input l(x), and an original input
x_j with coefficient 1 in l. Change coordinates to

    t=l(x),                z_i=x_i for i != j.

This is an affine bijection: x_j=t+lambda(z) for an affine form lambda.
In these coordinates the first AND has one pin t. Restrict t to b, for
b=0 and b=1. In either restriction that gate becomes affine and can be
removed. Replacing original inputs with their affine expressions and
performing this restriction adds no AND gate; it can be implemented by
affine operations preceding the old DAG, so no cycle is introduced.
Each restricted output h_b(z) consequently has a circuit with at most
a-1 AND gates. The induction hypothesis gives degree(h_b)<=a for both b.

For every Boolean (t,z), the unrestricted output is

    h(t,z) = h_0(z) + t*(h_0(z)+h_1(z)).                              (4)

The right side has degree at most a+1. Substituting the original affine
coordinates back cannot increase degree, so degree(h)<=a+1. This finishes
the induction, including arbitrary reuse of intermediate computations.
In particular, every signed circuit D for h satisfies

    A(D) >= degree(h)-1.                                             (5)

This proof is not the invalid assertion that ordinary arithmetic degree
increases by only one at each multiplication; the affine restriction
induction is what proves (5) for Boolean circuits with shared inputs.

## 6. The required degree of MOD-3 at every length

For any subset S of input positions, its monomial coefficient in h is

    a_S = sum_(T subseteq S) h(1_T)                     in F_2.       (6)

Here 1_T sets exactly positions in T to one and every other input to zero.
To verify (6), expand the unique polynomial and sum its values on those
assignments. A monomial on U contributes zero unless U is a subset of S.
For U subseteq S it contributes 2^(|S|-|U|) modulo 2, which is zero
unless U=S, when it is one. Thus precisely a_S remains.

For M_(n,t), (6) depends only on k=|S|. Let b_(k,t) denote the parity
of the number of subsets of a k-element set with size t modulo 3. For
k=1 the vector (b_(1,0),b_(1,1),b_(1,2)) is (1,1,0). Splitting subsets
according to whether they contain the last element gives

    b_(k+1,t) = b_(k,t) + b_(k,t-1)                    in F_2.       (7)

For every k>=1 the vector has two ones and its single zero is in position
-k modulo 3. Indeed this is true at k=1. If a three-entry vector has
its zero at z and its ones at z+1,z+2, recurrence (7) has its only zero
at z+2=z-1 modulo 3, proving the next case. This is a proof for every k,
not a finite computation or a periodicity assumption.

For n>=2, either b_(n,t)=1, giving the degree-n monomial, or
t=-n modulo 3. In the latter case t is different from -(n-1), so
b_(n-1,t)=1 and a degree-(n-1) monomial is present. Therefore

    degree(M_(n,t)) >= n-1            for all n>=2 and all t.        (8)

Combining (5) and (8), every normalized signed circuit for M_(n,t) has
A(D)>=n-2. Constants and output complements do not change the
positive-degree coefficients used here.

## 7. Converting the two inequalities to the one scored theorem

Take any original full B2 circuit for f_n, n>=2, with g gates, and let
D be its normalized signed representation from Section 2. Then s(D)<=g.
Both (3) and A(D)>=n-2 apply to this very same circuit D. Consequently

    3g >= 3s(D)
       = (3X(D)+2A(D)) + A(D)
       >= (6n-10) + (n-2)
       = 7n-12.

Dividing by three proves (1). The potential mu has not been equated to
gate count: the explicit identity 3s=mu+A supplies the needed conversion.
There is one output and no direct-sum assumption. All sixteen binary
tables, nonuniform DAGs, unrestricted sharing, and terminal-output cases
have been covered by the transfer and the deletion arguments above.

The canonical full expression is g_*(n)=(7/3)n, with polynomial exponent
1 and both logarithmic exponents 0. Given any real epsilon>0, let

    N_epsilon=max(2,ceil(12/(7*epsilon))).

For every integer n>=N_epsilon, (1) implies

    C_B2(f_n) >= (7/3)n-4 >= (1-epsilon)*(7/3)n.

Thus liminf over all integer lengths of C_B2(f_n)/g_*(n) is at least 1.
For the finite projection, given eta>0 use n>=max(2,ceil(4/eta)) to get
C_B2(f_n)>=(7/3-eta)n. The family and its verifier are fixed before any
of these quantifiers. Both structured score claims describe this exact
same joint theorem; no independent theorem supplies a coordinate.

## 8. Provenance and the extent of the claim

The weight 3 for parity gates and 2 for the other normalized gates, and
its combination with multiplicative complexity, are the method of
Arist Kojevnikov and Alexander S. Kulikov, "Circuit Complexity and
Multiplicative Complexity of Boolean Functions", CiE 2010, pp.239--245.
The consulted eight-page author version is at
https://www.lirmm.fr/~ashen/nafit/2010-kulikov.pdf
(SHA256 213372dab2d3e64f72204784b057de6e513907432628662d5dbb753f3e59bfbc).
Its Section 2.2 Lemma 1 and Section 3 Lemmas 2--3 explain the method.
Here all necessary facts are proved directly, including the degree
lemma and MOD-3 coefficient recurrence; no external theorem is a premise.

The earlier submission 46347322-9007-4a86-9e69-666592c7fad6 supplied the
same MOD-3 family, empty-witness verifier, signed-model transfer, and
three-cofactor viewpoint for the weaker bound 2n-4. The public note of
e3a45b32-6bc0-4471-a22e-e12e8d7aa084 was also consulted: it applies the
same published weighted method to MOD-3 and explains the 7/3 coefficient,
degree recurrence, and constant loss. Both are declared for attribution.
Their judgments and scores are not premises and no prior_results import
is used. The full argument required by this submission stands above.

The organizer's unchanged 31/10 baseline and its LP were read before
editing. Its disperser constructions and long case analysis are not used
in this proof. The optional LP checker is irrelevant to (1). There is
no computational premise or claimed discovery, optimality, superlinear
bound, or separation of complexity classes. This result strengthens our
previous coefficient 2 but remains below the organizer's literature seed.
