# C_B2(MOD3_n) >= 2n - 3 for every n >= 2: a self-contained proof

This file is the complete proof of the joint declaration in `bound.json`. It
uses no external theorem: the family, NP membership, circuit model, every
normalization step and the gate-elimination induction are proved here. The
argument is the classical Schnorr-style two-gates-per-variable elimination
(attribution only, not a premise), written out in full for the pinned
general B2 model. It claims no new
circuit lower bound beyond what is proved below and does not improve the
31/10 literature seed. Its coefficient is 2.

## 1. The family and NP membership

For n >= 0 and x in {0,1}^n let |x|_1 = x_1 + ... + x_n (an integer). Define

    f_n(x) = 1  iff  |x|_1 = 0 (mod 3).

In particular f_0(empty string) = 1 (the empty sum is 0). Every length,
including n = 0, 1, 2, is defined by this one formula; there are no exceptions.

Let L = {x in {0,1}* : |x|_1 = 0 (mod 3)}. One fixed algorithm A decides L at
every length: scan x once, keeping a counter c in {0,1,2}, replacing c by
(c+1) mod 3 on each 1; accept iff c = 0 at the end. A runs in O(n) time, so
L is in P. The NP verifier is V(x,w) = 1 iff w is the empty string and
A(x) accepts, with witness bound W(n) = 0.

- If x is in L, then w = empty satisfies |w| <= W(|x|) and V(x,w) = 1.
- If V(x,w) = 1 for some w, then A(x) accepts, so x is in L.

Both directions hold, so L is in NP with one fixed polynomial-time verifier
and one fixed polynomial witness bound. The characteristic function of L on
{0,1}^n is exactly f_n. The witness is not a circuit input; the circuits below
read only x, with n = |x| the actual input length. No block encoding or
padding is used, so there is no length loss.

## 2. Circuit model

We prove the lower bound in a model that contains every circuit allowed by
rules.md (Li-Yang general B2, ECCC TR21-023 Sections 2.1 and 2.4 conventions
as restated in rules.md). A *circuit* over variables y_1,...,y_m is a finite
directed acyclic graph whose nodes are

- input nodes, one for each variable (some may have no successors);
- any number of constant nodes labelled 0 or 1;
- gates; each gate G has an ordered pair of pins (u,v) pointing to nodes
  preceding it (u = v is allowed: repeated pins) and an arbitrary operation
  phi_G : {0,1}^2 -> {0,1} (any of the 16 truth tables).

One node (of any type) is selected as the output. Fan-out and depth are
unrestricted and subgraphs may be shared. The value of a gate is
phi_G(val(u), val(v)). The *size* |C| is the number of gates; input nodes,
constant nodes, wires and output selection cost zero, and every gate,
including degenerate or useless ones, counts one. For a Boolean function h
let s(h) be the minimum size of a circuit computing h exactly on every input.

Every circuit of the pinned B2 model is a circuit in this sense (if that model
has no constant nodes, it is the special case with none). Hence
C_B2(f_n) >= s(f_n), and it suffices to prove s(f_n) >= 2n - 3 for n >= 2.
Allowing free constant nodes only makes the lower bound stronger; this is
exactly why substitutions below may create constant nodes without cost.

Renaming variables does not change s, and adding unused input nodes does not
change the function computed or the size.

## 3. Normalization lemmas

Call a gate G with pins (u,v) and operation phi **degenerate** if

- (a) u = v, or
- (b) u or v is a constant node, or
- (c) phi does not depend on one of its two arguments.

**Lemma 3.1.** If G is degenerate, there is a node w in {u,v} (w is not G,
by acyclicity) such that, on every input, val(G) is a constant, or val(G) =
val(w), or val(G) = 1 - val(w).

*Proof.* In case (a), val(G) = phi(val(u), val(u)) = psi(val(u)) for the
unary function psi(a) = phi(a,a). In case (b), say u is the constant c; then
val(G) = phi(c, val(v)) = psi(val(v)) with psi(a) = phi(c,a); symmetrically
for v. In case (c), phi depends on at most one argument, so val(G) = psi of
that pin's value. Every unary psi is one of 0, 1, identity, negation. []

**Lemma 3.2 (negation transfer).** Let w be a gate of a circuit C. Replace
phi_w by its complement 1 - phi_w, and for every gate K having a pin on w,
replace phi_K by the operation obtained by complementing each argument that
is read from w (both arguments if both pins of K are on w). The resulting
circuit C' has the same size, val'(w) = 1 - val(w), and every other node has
the same value as in C on every input.

*Proof.* Induction in topological order. Nodes before w are unchanged; w is
complemented by construction. For K after w, each argument read from w is
complemented twice (once in val'(w), once inside phi'_K), and all other
arguments are unchanged by induction, so val'(K) = val(K). Each modified
operation is again one of the 16 binary truth tables. []

**Lemma 3.3 (removing a degenerate gate).** Let C compute h, and let G be a
degenerate gate, with w and the case as in Lemma 3.1. Then there is a circuit
of size |C| - 1 computing h, except possibly when G is the output and
val(G) = 1 - val(w) with w not a gate.

*Proof.* We remove G and keep the value of every remaining node.

- val(G) is a constant c: redirect every pin on G to a constant node c (add
  one, free, if needed); if G is the output, select that constant node.
- val(G) = val(w): redirect every pin on G to w; if G is the output, select w.
- val(G) = 1 - val(w) and G is not the output: for each gate K with a pin on
  G, redirect those pins to w and complement the corresponding arguments of
  phi_K (as in Lemma 3.2). Each new operation is a binary truth table and
  K's value is unchanged.
- val(G) = 1 - val(w), G is the output and w is a gate: first apply Lemma 3.2
  to w. Now val'(w) = 1 - val(w), and G's operation was adjusted so that G
  keeps its value; thus val(G) = val'(w), and the previous case "val(G) =
  val(w)" removes G.

Redirecting a pin from G to a predecessor w of G, or to a constant node,
cannot create a cycle because w precedes every successor of G. Deleting G
then leaves an acyclic circuit with one fewer gate computing h. []

## 4. The function classes

For m >= 1 and r in {0,1,2} let g_{m,r}(y_1..y_m) = 1 iff |y|_1 = r (mod 3),
and let F_m = {g_{m,r}, 1 - g_{m,r} : r in {0,1,2}}. Then f_n = g_{n,0} is
in F_n. Indices of r are taken mod 3 throughout.

**Lemma 4.1 (closure).** For h in F_m with m >= 2, any variable y_i and any
a in {0,1}, the restriction h|_{y_i=a}, viewed as a function of the other
m-1 variables, is in F_{m-1}: g_{m,r}|_{y_i=a} = g_{m-1,r-a}, and the
negations restrict to the negations.

*Proof.* The remaining sum equals |y|_1 - a. []

**Lemma 4.2 (dependence).** For k >= 2, every h in F_k depends on every one
of its k variables. In particular h is not constant, and h is not equal to a
variable, a negated variable, or a constant.

*Proof.* Negation preserves dependence, so take h = g_{k,r}. Fix y_i. The
other k-1 >= 1 variables can be given sum t = 0 or t = 1 (all zero, or a
single one). Choose t = 0 if r is 0 or 1, and t = 1 if r = 2. Then setting
y_i = 0, 1 gives total sum t, t+1, i.e. {0,1} or {1,2}; exactly one of these
is congruent to r mod 3, so flipping y_i changes h. Since h depends on
k >= 2 variables, it is not constant and not a (negated) single variable. []

**Lemma 4.3 (three distinct pair restrictions).** For m >= 3, h in F_m and
distinct variables y_i, y_j, the four restrictions h|_{y_i=a, y_j=b},
(a,b) in {0,1}^2, include at least three pairwise distinct functions of the
remaining m-2 >= 1 variables.

*Proof.* For h = g_{m,r} (negations behave identically) the restriction is
g_{m-2, r-a-b}. The pairs (0,0), (0,1), (1,1) give residues r, r-1, r-2,
which are the three distinct residues mod 3. It remains that g_{k,0}, g_{k,1},
g_{k,2} are pairwise distinct for k >= 1. If k >= 2, for each s in {0,1,2}
there is an assignment with sum exactly s; there g_{k,s} = 1 and
g_{k,s'} = 0 for s' != s. If k = 1, they are 1-y, y and the constant 0,
which are pairwise distinct. Complementing all three preserves distinctness. []

## 5. Structure of a minimum circuit

**Lemma 5.1.** Let m >= 3, h in F_m, and let C be a circuit computing h with
|C| = s(h). Then every gate of C is non-degenerate: its two pins are distinct
nodes, neither is a constant node, and its operation depends on both
arguments. Moreover s(h) >= 1 and the output is a gate.

*Proof.* By Lemma 4.2, h is neither constant nor a literal, so a circuit with
zero gates (whose output is an input or constant node) cannot compute h;
thus s(h) >= 1 and the output is a gate. Suppose some gate G is degenerate.
If Lemma 3.3 applies, it yields a circuit for h of size s(h) - 1, contradicting
minimality. The only exception is G being the output with val(G) = 1 - val(w)
for an input or constant node w; then h is a negated variable or a constant,
contradicting Lemma 4.2. []

**Lemma 5.2.** In such a minimum circuit there is a gate G whose two pins are
distinct input nodes y_i, y_j (i != j), and at least one of y_i, y_j has a
pin from some gate G' != G.

*Proof.* Take a gate G that is first in a topological order (one exists since
s(h) >= 1). Its pins point to non-gate nodes, which by Lemma 5.1 are distinct
and not constants; so they are input nodes y_i, y_j with i != j.

Suppose that G is the only gate with a pin on y_i and the only gate with a pin
on y_j. Then every node other than y_i, y_j computes a function of
val(G) and the variables other than y_i, y_j: by induction in topological
order, inputs other than y_i, y_j and constants trivially qualify, G
qualifies, and every other gate reads only nodes other than y_i, y_j (only G
reads them). The output is a gate other than y_i, y_j, so
h = H(val(G)(y_i,y_j), rest) for some function H. Then, for each fixed
assignment of the rest, h|_{y_i=a,y_j=b} depends on (a,b) only through the
single bit val(G)(a,b), so among the four restrictions there are at most two
distinct functions (those with val(G)(a,b) = 0 coincide, as do those with
val(G)(a,b) = 1). This contradicts Lemma 4.3. Hence one of y_i, y_j, say
y_i after renaming, is read by a gate G' != G. []

## 6. Elimination step

**Lemma 6.1.** For m >= 3 and h in F_m there is a variable y_i such that,
for a = 0 (indeed for either a), s(h|_{y_i=a}) <= s(h) - 2, where
h|_{y_i=a} is in F_{m-1}.

*Proof.* Let C be a minimum circuit for h, and take G, G', y_i, y_j from
Lemma 5.2. The gate G' != G has a pin on y_i; by Lemma 5.1 its other pin is a
node z != y_i. G is not the output: otherwise h would depend only on y_i, y_j,
contradicting Lemma 4.2 (h depends on all m >= 3 variables).

Substitute y_i = a: replace the input node y_i by a constant node a, keeping
all pins. The resulting circuit C_1 has s(h) gates and computes
h' = h|_{y_i=a}, which is in F_{m-1} by Lemma 4.1, with m - 1 >= 2.

Step 1. In C_1, G has a constant pin, so G is degenerate. G is not the output,
so Lemma 3.3 removes it: a circuit C_2 with s(h) - 1 gates computes h'. This
modification only redirects pins that pointed to G and changes operations of
gates that read G; the pin of G' on the constant node a is untouched.

Step 2. In C_2, G' still has a pin on the constant node a, so G' is
degenerate (its other pin may now be a different node or also a constant;
this does not matter). If G' is not the output, Lemma 3.3 removes it. If G'
is the output, by Lemma 3.1 h' is a constant, val(w) or 1 - val(w) for a pin
node w of G'. By Lemma 4.2 (applied with k = m - 1 >= 2), h' is not constant
and not a literal or constant node value, so w is a gate and Lemma 3.3
applies (its exceptional case does not occur). Either way we obtain a circuit
with s(h) - 2 gates computing h'.

G and G' are distinct gates of C, so exactly two gates were deleted. []

## 7. The theorem

**Theorem.** For every m >= 2 and h in F_m, s(h) >= 2m - 3. In particular,
for every integer n >= 2, C_B2(f_n) >= s(f_n) >= 2n - 3.

*Proof.* Induction on m. For m = 2, h depends on both variables (Lemma 4.2),
so it is not computed by a circuit with no gates (whose output is a single
input or constant node); thus s(h) >= 1 = 2*2 - 3. For m >= 3, Lemma 6.1
gives h' in F_{m-1} with s(h) >= s(h') + 2 >= 2(m-1) - 3 + 2 = 2m - 3.
Apply this to h = f_n = g_{n,0} in F_n and use C_B2 >= s from Section 2. []

## 8. The declared guarantee and projection

The declaration is a version-2 `eventual_lower_bound` with leading term
L(n) = 2*n^1*(log2 n)^0*(log2 log2 n)^0 = 2n and explicit remainder
R(n) = 3*n^0 = 3. The theorem gives C_B2(f_n) >= L(n) - R(n) = 2n - 3 for
EVERY integer n >= N with N = 2; there are no skipped or sparse lengths, and
the family, verifier and witness bound are fixed before all quantifiers. The
remainder coefficient 3 is the exact declared constant (not an unknown K).

The finite projection (c,p) = (2,1) follows: (2n - 3)/n = 2 - 3/n tends to 2,
so for every epsilon > 0 and every n >= max(2, 3/epsilon),
C_B2(f_n) >= (2 - epsilon)n; hence liminf C_B2(f_n)/n >= 2. The loss 3 is
o(n). This is exactly the projection implied by the protected algebra in
`expressions.py` for an eventual bound whose leading term is 2n. Both score
transport fields are coordinates of this one theorem: coefficient 2,
exponent 1.

## 9. Proofs of the unscored string claims, and a remark

1. *Restriction lemma (Lemma 6.1).* Proved in Section 6.
2. *Eligibility.* L is in P and in NP with V and W = 0 as in Section 1; both
   directions are proved there.

*Remark (not a claim).* f_0 = 1 and f_1(y) = 1 - y; f_1 is computed by a single
   negation gate, and s(f_1) >= 1 because f_1 is neither a constant nor an
   un-negated variable; the main theorem's bound 2n - 3 is negative or zero
   for n <= 1 and is not claimed there. For n = 2 the bound is 1, which holds
   since f_2 = NOR(y_1,y_2) is a single gate and depends on both variables.

## 10. Scope and limitations

- The coefficient 2 is below the 31/10 literature seed of this challenge and
  below other known explicit B2 lower bounds of the form 3n - o(n) or more.
  These are context only, not premises. This theorem is submitted because every
  step is self-contained and verifiable without relying on external case
  analyses.
- No upper bound is claimed; the lower bound is not asserted to be tight.
- No finite computation or checker is used, so nothing here relies on
  experimental coverage.
