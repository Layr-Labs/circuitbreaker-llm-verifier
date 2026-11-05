# MOD3 in full B2: an explicit lower bound of 2n - 2

## 0. The claim and its provenance

For every integer n >= 0 and r in {0,1,2}, define

    M(n,r)(x) = 1 iff the Hamming weight of x is congruent to r modulo 3.

The ONE scored family is f_n = M(n,0), on exactly n input bits. We prove

    C_B2(f_n) >= 2n - 2 for every integer n >= 4.                 (T)

In fact the same inequality holds for each of the three auxiliary residues.
The auxiliaries serve an induction; they are not different families used to
supply the two score coordinates. The declaration's leading expression is
2n, its explicit remainder is the constant 2, and its finite projection is
(2,1). No superlinear lower bound, separation of complexity classes, exact
upper bound, leaderboard record, or historical novelty is claimed.

The restriction and two-gate deletion strategy is adapted from kaankolcu's
Yukon submission 466a2973-4bf2-4e79-be86-b9c6a87260d5, original source commit
8c58caa6d5bd7dd75ef67c792c11255658e1d870, file submission/proof.md, especially
Sections 3-7. That source proves 2n-3 for MOD3. Its proof was read while
preparing this submission and is declared as a dependency. The contribution
developed here is the explicit three-input obstruction in Section 5 and the
choice of a restriction value in Section 6, improving that finite estimate
by one gate. We do not claim that this refinement is new in the literature.
No previous evaluation or approval is a premise. All arguments needed for
(T), including normalization, are provided below.

## 1. A total uniform language in P and NP

Let L be the set of all finite binary strings whose number of ones is
divisible by 3. On the empty input f_0=1; on one bit f_1(x)=1-x. The same
definition applies at all other lengths. A single deterministic evaluator
starts a counter at zero, updates it modulo 3 for each input bit, and accepts
exactly when the final counter is zero. It takes O(n) time and constant
counter space. It has no advice, padding, certificate input, or parameter
chosen differently for each length.

Set W(n)=0 and let V(x,w) accept exactly when w is empty and the evaluator
accepts x. If x is in L, the empty witness is accepted. Conversely, if there
is a witness of length at most W(|x|) accepted by V, the evaluator accepts x,
so x is in L. This proves both directions of NP membership with one fixed
polynomial-time verifier. Circuit inputs below are x alone, of actual length
n; the witness introduces no input bits.

## 2. Model and normalization with every gate charged

We work with acyclic deterministic circuits having one output, all sixteen
binary Boolean truth tables, two ordered pins per gate, unrestricted fan-out
and depth, and exact computation on every assignment. Each gate is counted
once. Repeated pins, shared descendants, and degenerate truth tables are
permitted. Input nodes, wires and output selection do not count as gates.
We also allow any number of free constant source nodes. This is at least as
powerful as the challenge's full B2 model, so a lower bound proved here
transfers to that model: its minimum size cannot be smaller. We never allow
free gates implementing a new basis or free computational edges.

A gate is proper if its pins are distinct nonconstant nodes and its truth
table depends on both pin values. Otherwise it is degenerate. At a
degenerate gate v its value is a constant, a pin value, or the complement
of a pin value. This follows by substituting a constant pin, identifying
repeated pins, or deleting an argument ignored by the truth table. A unary
Boolean function has exactly these four possibilities.

Normalization lemma. If a circuit computes a function depending on at
least two variables, a degenerate gate can be deleted without adding any
gate. Moreover every other currently degenerate gate remains degenerate
until it too is removed. Consequently k distinct degenerate gates permit
at least k gate deletions, even when their descendants are shared.

Proof. If v computes a constant, replace every use of v by that constant
source. If it forwards a pin w, redirect every use to w. If it complements
w, redirect each use to w and complement the corresponding argument in the
receiving gate's truth table. Do this for every occurrence, including both
pins if necessary. Each receiving gate still has an allowed binary truth
table and its value is unchanged. The substituted node precedes v, so
acyclicity is preserved. No descendant is duplicated.

If v is the output, a constant value or a function of an input node alone
cannot equal the assumed function of at least two variables. In the
forwarding case designate w as output. In the complement case w must be
another gate. Complement w's truth table, designate w as output, and
after the redirections above, compensate at every surviving gate pin now
pointing to w by complementing that argument in its receiving truth table.
The value of w becomes the old value of v; all other
surviving gate values remain unchanged. This repairs a terminal negation
without introducing a free NOT gate.

A pre-existing constant pin remains constant or is supplemented by a new
constant pin. Equal pins are redirected identically and remain equal.
Complementing truth-table arguments or the whole output cannot create
dependence on an ignored argument. Thus every other degenerate gate stays
degenerate. Each operation removes exactly one gate, so the process
terminates. Removing gates outside the output's ancestor subgraph afterward
only decreases size and does not change retained pins. This proves the
lemma and yields a circuit in which every gate is proper and relevant. Its
output is a gate, since a constant or an input is insufficient. QED.

All lower bounds will be proved for such normalized circuits. Given any
original circuit, normalization supplies one with no more gates, so a
lower bound for the normalized circuit is also a lower bound for the
original one.

## 3. Residues and the one-bit bottleneck

Restriction identity. Fixing any variable to b in {0,1} turns M(n,r) into
M(n-1,r-b), where residues are reduced modulo 3 and the other variables
are relabeled in their original order. This follows directly by subtracting
the fixed bit from the Hamming weight.

Full dependence. For n>=2 every M(n,r) depends on every input. There is an
assignment of weight r because r<=2<=n. Flipping any chosen coordinate of
that assignment changes its weight by 1 or -1; the output changes from 1
to 0 because neither change is divisible by 3.

Three distinct pair restrictions. For n>=3, fixing two distinct inputs to
00, 01, and 11 leaves respectively M(n-2,r), M(n-2,r-1), and M(n-2,r-2).
For m=n-2>=1 the three residue functions on m bits are pairwise distinct:
the all-zero assignment singles out residue 0, and an assignment of weight
one singles out residue 1. This also distinguishes residue 2 from each of
the other two when it is the zero function at m=1.

Thus M(n,r) cannot have a representation

    F(g(x_i,x_j), x_rest)                                      (B)

for any Boolean g and F and any pair i!=j when n>=3. Once x_i,x_j are
fixed, (B) has at most two possible residual functions, F(0,.) and F(1,.),
contrary to the three just exhibited. This observation allows arbitrary
sharing inside F; it is an information bottleneck, not a formula bound.

## 4. Deleting two gates for either restriction value

For n>=3, take a normalized circuit C for M(n,r). The first gate G in a
topological ordering has two distinct primary inputs, say x_i and x_j:
it has no earlier gates and proper gates have no constant or repeated pins.
If both inputs feed no other gate, all their influence on the output passes
through the single value of G. This is the forbidden representation (B).

Therefore one of these input nodes, call it x, feeds at least two DISTINCT
gates. Fan-out here counts successor gates, not pin occurrences. Fix x to
any b in {0,1}. Both successors have a constant pin, so they are two
distinct degenerate gates. The restricted function is M(n-1,r-b); it
depends on all n-1>=2 remaining variables. The normalization lemma deletes
both successors, even if one is an ancestor of the other. It cannot
resurrect the second as a proper gate. Thus

    |C| >= C_B2(M(n-1,r-b)) + 2, for EACH b in {0,1}.           (D)

The same identified input works for both choices of b. Equation (D) is
not a direct-sum assertion and does not require disjoint descendant cones.

For completeness M(2,r) needs at least one gate, since it depends on two
variables. Applying (D) gives a preliminary bound of three gates for every
M(3,r). This will rule out circuits of fewer than three gates in the next
section; their impossibility is not assumed.

## 5. A three-input terminal obstruction

Write E1(x,y,z)=1 iff exactly one of its three inputs is 1, and E2=1 iff
exactly two are 1. These are M(3,1) and M(3,2), respectively. We show that
each requires at least FOUR gates in the full model above. No computer
enumeration is used.

A Boolean function of y,z is called affine here if it has the form
alpha*y XOR beta*z XOR gamma with alpha,beta,gamma in {0,1}. Composition
with unary Boolean functions, or complementing the output, preserves
affinity. XOR and XNOR are affine. AND and NOR are not: the XOR of their
values on 00,01,10,11 is 1, while for every affine function it is 0.

For any chosen variable x, the two restrictions of E1 are

    E1(0,y,z)=y XOR z;    E1(1,y,z)=NOR(y,z).

The restrictions of E2 are

    E2(0,y,z)=AND(y,z);   E2(1,y,z)=y XOR z.

In both cases one slice is parity and the other is nonaffine. Both
functions are symmetric in the variable names, so this statement holds
whichever variable is named x.

Suppose either function had a circuit with at most three gates. Normalize
and trim it. Section 4 gives a lower bound of three, so there are exactly
three proper relevant gates. Name them A,B,H in topological order. H is
the output. There are only the following possibilities, which exhaust the
graph without imposing a fan-out restriction.

Case I: B does not use A. Both A and B have primary-input pins. Relevance
forces both to feed H. Their two-element input sets together must contain
all three variables. Consequently they overlap in exactly one variable,
and after renaming the representation is

    H(A(x,y), B(x,z)).                                        (I)

Fix x to the parity-slice value b (b=0 for E1, b=1 for E2). Both A(b,y)
and B(b,z) must depend on their variables, or (I) would miss that variable.
Each is therefore a bijection on {0,1}, either identity or complement.
As y,z vary, H is tested on all four possible input pairs. Since (I) is
parity on this slice, its full truth table must be XOR or XNOR. Now on
the other slice A(1-b,y) and B(1-b,z) are unary functions, possibly
constant. XOR or XNOR of two unary functions is affine in y,z. This
contradicts the nonaffine other slice.

Case II: B uses A. A is a gate on two primary variables u,v. Its other
input at B is a primary variable w; there is no other earlier gate, and
properness forbids using A twice. H must use B for relevance. Its other
pin is either A or a primary variable t.

If w is one of u,v, the third variable z has not appeared in A or B.
Essential dependence on z forces the other input of H to be z. But
then the circuit is H(B(A(u,v),w),z), a representation of type (B) in
which the Boolean pair function is B(A(u,v),w). Section 3 forbids it.

It remains that w is the third variable z, so B=B(A(u,v),z). If H's
other input is A or z, all influence of u,v again passes through the
single bit A. Section 3 forbids both possibilities. The only remaining
case has H's other input equal to u or v. Renaming that input x and the
other member of the pair y gives

    H(B(A(x,y),z),x).                                         (II)

Fix x to the parity-slice value b. The unary map A(b,y) must be a
bijection, and the unary map s -> H(s,b) must also be a bijection;
otherwise the output cannot be parity in y,z. Inverting those identity
or complement maps shows that B's FULL truth table on its two inputs
is XOR or XNOR. On the other slice B(A(1-b,y),z) is therefore affine
in y,z, and composing it with the unary map s -> H(s,1-b) keeps it
affine (including the constant cases). This contradicts the other slice.

There is no third topology: if B is independent of A it is Case I;
otherwise it is Case II. Unused gates were trimmed; repeated or constant
pins and terminal unary gates were eliminated in Section 2; all three
variables are essential. Shared use of A in B and H was considered in
Case II, rather than silently excluded. We have proved

    C_B2(M(3,1)) >= 4 and C_B2(M(3,2)) >= 4.                   (E)

No equality or existence of a four-gate circuit is asserted.

## 6. All sufficiently large integer lengths

Take any residue r at n=4. Section 4 identifies an input with two distinct
successor gates. Choose the fixed value b so that r-b is 1 or 2 modulo 3:
for r=0 take b=1; for r=1 or r=2 take b=0. Equation (D) and (E) imply

    C_B2(M(4,r)) >= 2 + 4 = 6 = 2*4 - 2.

Suppose now that n>=5 and the stated bound holds for every residue at
length n-1. Apply (D), choosing for instance b=0. It gives, for each r,

    C_B2(M(n,r)) >= 2 + [2(n-1)-2] = 2n-2.

Induction proves this for every integer n>=4 and every residue. In
particular it proves (T) for the fixed scored family f_n=M(n,0).

This is an eventual bound with the explicit threshold N=4 and constant
remainder 2, not merely a liminf statement. Its declaration is exactly
L(n)=2*n^1*(log2 n)^0*(log2 log2 n)^0 and a single remainder monomial
2*n^0*(log2 n)^0*(log2 log2 n)^0. Thus it uses no rounding of an unsupported
expression. Dividing by n yields C_B2(f_n)/n >= 2-2/n. For any epsilon>0,
N_epsilon=max(4,ceil(2/epsilon)) suffices for C_B2(f_n)>=(2-epsilon)n at
EVERY integer n>=N_epsilon. This proves the matching finite projection
(2,1) for that same family. No sparse-length argument, parameter-dependent
family, encoding loss, external theorem, or computational witness is used.

## 7. Scope of this entry

The earlier attributed argument supplies a 2n-3 estimate. The explicit
terminal proof and the adaptive choice of b at n=4 above supply 2n-2.
This is a finite-estimate refinement, not an improved leading coefficient.
The entry is below the challenge's 31/10 literature seed and the currently
displayed coefficient-four candidates. It does not establish a frontier
advance, a prize entitlement, or a formal machine-checked proof. Its
mathematical content must be assessed from the complete proof above.

Source of the adapted method:
https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/8c58caa6d5bd7dd75ef67c792c11255658e1d870/submission/proof.md

This link documents provenance only; no retrieval of it is needed to check
the submitted argument. Neither the source's text nor any previous review
is being incorporated as an unproved premise.
