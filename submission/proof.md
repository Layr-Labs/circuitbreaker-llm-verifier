# Coefficient 2, exponent 1 for MOD3 by gate elimination: complete proof

This submission proves a weaker, fully self-contained theorem than the 31/10
literature seed. It is the classical gate-elimination argument for MOD3-type
functions (Schnorr 1974; Wegener, *The Complexity of Boolean Functions*, 1987,
Chapter 5, which covers functions in the class Q_{2,3}). It claims no new circuit
lower bound. Every step is proved here, so no external theorem is used as a
premise. The citations are only for attribution.

## 0. Model and notation

Fix a finite set X of input variables. A circuit D over X is a finite directed
acyclic graph whose nodes are

* one input node for each variable in X,
* the constant nodes 0 and 1 (each may occur at most once; they cost zero),
* gates. Each gate G has an ordered pair of pins (u, v) pointing to nodes that
  come strictly earlier in a fixed topological order, and an arbitrary truth
  table phi_G : {0,1}^2 -> {0,1} (any of the 16).

Repeated pins (u = v), pins to constant nodes and unrestricted fan-out are
allowed. One node is designated as the output. On an assignment
x in {0,1}^X every node gets a value: inputs get their bit, constants get
their constant, and gate G with pins (u, v) gets phi_G(val u, val v). D computes
h : {0,1}^X -> {0,1} if val(output) = h(x) for every x. The size |D| is the
number of gates. Inputs, constant nodes, wires and the choice of output cost
nothing. C(h) is the minimum size of a circuit computing h. This is the
Li–Yang B2 model of rules.md (ECCC TR21-023, Sections 1.1, 2.1, 2.4): we
count every gate, including degenerate gates, in the circuit being bounded.
Degenerate gates are only *removed* by the size-reducing transformations in
Lemma 2, which can only help an upper bound and are used only inside a proof
by minimality. A "literal" is a function of the form b, x_i or not-x_i
(constant, projection, or negated projection).

The fan-out of an input node x_i is the number of distinct gates that have
at least one pin equal to x_i.

## 1. The family, uniformity and NP membership

Let L be the set of all finite binary strings (including the empty string)
whose number of 1-bits is divisible by 3. For each integer n >= 0 let
f_n : {0,1}^n -> {0,1}, f_n(x) = 1 iff x_1 + ... + x_n is congruent to 0
(mod 3). So f_0 = 1, and f_n is defined at **every** length with no exceptions
or advice.

One fixed algorithm A reads x, keeps a counter mod 3 in O(1) states, and
accepts iff the counter is 0 at the end. It runs in time O(|x|) at every
length, so L is in P. For NP, take witness bound W(n) = 0 and verifier
V(x, w) = 1 iff w is the empty string and A(x) accepts. Then V runs in
polynomial time. If x is in L, then w = empty satisfies |w| <= 0 and
V(x, w) = 1. Conversely, if V(x, w) = 1 for some w with |w| <= W(|x|), then A
accepts x, so x is in L. Hence L is in NP. The witness is never a circuit
input. The circuits for f_n read only the n bits of x, and n = |x| is the
actual input length. There are no blocks, padding or encodings.

For a variable set Y with |Y| = m and a residue a in Z_3, write

    M_Y^a(y) = 1  iff  sum_{y in Y} y  is congruent to a (mod 3).

C(M_Y^a) depends only on m and a (renaming variables preserves circuit size),
so we write M_m^a. Note f_n = M_n^0.

## 2. Properties of M_m^a

**Lemma 1.** Let |Y| = m and a in Z_3.

(P1) Restriction. For y_i in Y and c in {0,1}, M_Y^a restricted to y_i = c is
     M_{Y \ {y_i}}^{a-c}.

(P2) Pair test. If m >= 4 and y_i, y_j are distinct variables of Y, the four
     subfunctions M_Y^a|_{y_i=b_1, y_j=b_2} ((b_1,b_2) in {0,1}^2) include at
     least 3 pairwise distinct functions of the other m-2 variables.

(P3) Non-literal. If m >= 2, M_Y^a is neither constant nor a literal.

*Proof.* (P1) With y_i = c the sum of all variables is c plus the sum of the
others, so the sum is congruent to a iff the sum of the others is congruent
to a - c.

(P2) Let Y' = Y \ {y_i, y_j}, with |Y'| = m-2 >= 2. By (P1) applied twice, the
subfunction for (b_1, b_2) is M_{Y'}^{a-s} with s = b_1 + b_2. The pairs (0,0),
(1,0) and (1,1) give s = 0, 1, 2, i.e. the three residues a, a-1, a-2, which are
pairwise different in Z_3. Because m-2 >= 2, the sums of the variables of Y'
take every value 0, 1, 2. So for each residue e there is an assignment of Y'
with sum congruent to e. On that assignment M_{Y'}^e = 1, while M_{Y'}^{e'} = 0
for every e' != e. Hence the three functions M_{Y'}^a, M_{Y'}^{a-1} and
M_{Y'}^{a-2} are pairwise distinct.

(P3) Since m >= 2, some assignment has sum congruent to a, so the function
takes value 1 there. Flipping any one variable changes the sum by +1 or -1,
so the residue changes and the value becomes 0. Therefore the function is not
constant, and it depends essentially on every variable, i.e. on at least
m >= 2 variables. A literal depends on at most one variable. QED.

## 3. Removing a degenerate gate

Call a gate G of a circuit D *degenerate* if there is a node w that is a pin
of G, or a constant node, and a map l : {0,1} -> {0,1} such that
val(G) = l(val w) on every assignment, **as forced by the syntax**, namely in
one of these cases:

(a) a pin of G is a constant node b. Then val G = phi_G(b, val v) (or
    phi_G(val u, b)). Take w the other pin and l = phi_G(b, .) (or phi_G(., b)).
(b) both pins are the same node u. Take w = u and l(z) = phi_G(z, z).
(c) phi_G depends on at most one of its two arguments. Take w the pin of the
    argument it depends on (either pin if neither) and l the induced map.

In each case w precedes G in the topological order.

**Lemma 2.** Let D compute h, where h is neither constant nor a literal, and
let G be a degenerate gate of D with witness (w, l). Then there is a circuit
D' with |D'| = |D| - 1 that computes h. D' is obtained from D by deleting G,
redirecting some pins that pointed to G, and changing some truth tables.
All other gates and all pins that did not point to G keep their nodes. In
particular, any pin that pointed to a constant node in D still does in D'.

*Proof.* There are four cases, by l.

1. **l is the constant b.** h is not constant, so G is not the output.
   Redirect every pin that points to G to the constant node b (create it at
   cost zero if absent). Delete G. Every remaining node keeps its value.
2. **l is the identity.** Redirect every pin pointing to G to w. If G was the
   output, make w the output. Delete G. Values are preserved, and the output
   still computes h.
3. **l is negation, and G is not the output.** Redirect every pin pointing to
   G to w, and in each gate K that has such a pin, change phi_K by negating
   the corresponding argument(s): phi'_K(z_1,z_2) = phi_K(z_1 xor e_1, z_2 xor e_2),
   where e_t = 1 iff pin t of K was redirected. Then K sees not-val(w) = val(G)
   as before, so all values are preserved. Delete G.
4. **l is negation, and G is the output.** If w were an input or a constant
   node, h would be a literal or a constant, which is excluded. So w is a
   gate. Replace phi_w by its negation, so the new val(w) equals the old
   not-val(w) = old val(G) = h. For every gate K other than G that has a pin
   to w, negate the corresponding argument(s) of phi_K as in case 3, so K's
   value is unchanged. Redirect every pin that pointed to G to w, keeping its
   table (the new val(w) equals the old val(G)). Make w the output and delete G.
   All gates other than w keep their values, and the output computes h.

Redirecting a pin of a gate K (which comes after G) to w or a constant (which
comes before G) preserves acyclicity in the same topological order. Exactly
one gate is deleted. Pins that did not point to G are untouched. QED.

**Corollary 3.** If D is a minimum-size circuit for a function h that is
neither constant nor a literal, then no gate of D is degenerate. In
particular, every gate of D has two distinct pins, and neither pin is a
constant node.

*Proof.* Otherwise Lemma 2 gives a smaller circuit for h. QED.

## 4. The elimination step

**Lemma 4.** For every m >= 4 and a in Z_3 there is a' in Z_3 with
C(M_m^a) >= C(M_{m-1}^{a'}) + 2.

*Proof.* Let h = M_Y^a with |Y| = m >= 4. Let D be a minimum-size circuit for
h. By (P3) h is neither constant nor a literal. So the output is a gate, and
D has at least one gate. By Corollary 3 no gate of D is degenerate.

*Step 1 (a bottom gate).* Let G be the first gate in the topological order.
Its pins point to earlier nodes, which are not gates, so they are inputs or
constants. By Corollary 3 they are not constants and are distinct. So the pins
of G are two distinct input nodes y_i and y_j.

*Step 2 (one of them has fan-out at least 2).* Suppose both y_i and y_j have
fan-out 1, i.e. G is the only gate with a pin to either. Neither is the output,
because the output is a gate. Form D* from D by turning G into a new input node
z, removing its pins. In D* the nodes y_i and y_j have no successors and are not
the output. So D* computes a function F(z, y_rest) that does not depend on y_i
or y_j, where y_rest is the other m-2 variables. On every assignment, D gives G
the value phi_G(y_i, y_j), and all other node values agree with D* evaluated at
z = phi_G(y_i, y_j). Hence h(y) = F(phi_G(y_i, y_j), y_rest) identically. Then
for every (b_1, b_2) the subfunction h|_{y_i=b_1, y_j=b_2} equals
F(phi_G(b_1,b_2), .), which is one of only the two functions F(0, .) and
F(1, .). This contradicts (P2), which needs m >= 4 and gives at least 3
distinct subfunctions. Hence, by symmetry, we may assume y_i has fan-out at
least 2. Since G's two pins are distinct, G has exactly one pin to y_i. So
there is a gate H != G that also has a pin to y_i.

*Step 3 (substitute y_i = 0 and remove G and H).* Let D_0 be D with every pin
to y_i redirected to the constant node 0 and the input node y_i removed.
|D_0| = |D|. D_0 computes h_0 = h|_{y_i=0} = M_{Y\{y_i}}^{a} by (P1). That
function has m-1 >= 3 >= 2 variables, so by (P3) it is neither constant nor a
literal. In D_0 both G and H have a pin to the constant node 0, so both are
degenerate (case (a)).

Apply Lemma 2 to G in D_0. This gives D_1 with |D_1| = |D| - 1 computing h_0.
H is a different gate, so it survives in D_1. Its pin that pointed to y_i, and
then to the constant 0, did not point to G. By the last sentence of Lemma 2 it
still points to a constant node (H's other pin may have been redirected, and
its table changed, which does not matter). So H is degenerate in D_1 by case
(a). Apply Lemma 2 to H in D_1. This gives D_2 with |D_2| = |D| - 2 computing
h_0. Therefore

    C(M_{m-1}^a) = C(h_0) <= |D| - 2 = C(M_m^a) - 2.

This is the claim with a' = a. QED.

## 5. The theorem and the scored projection

**Theorem.** For every integer m >= 3 and every a in Z_3,
C(M_m^a) >= 2m - 6. In particular C_B2(f_n) >= 2n - 6 for every integer
n >= 3.

*Proof.* Induction on m. For m = 3, 2m - 6 = 0 <= C(M_3^a). For m >= 4,
Lemma 4 and the induction hypothesis give
C(M_m^a) >= C(M_{m-1}^{a'}) + 2 >= 2(m-1) - 6 + 2 = 2m - 6. Taking a = 0 and
m = n gives the bound for f_n = M_n^0. QED.

**Joint (2, 1) guarantee.** The declaration has c = 2, p = 1, q = r = 0, so
g(n) = 2n. For every integer n >= 3, C_B2(f_n)/g(n) >= (2n - 6)/(2n) = 1 - 3/n.
Given epsilon > 0, let N_epsilon = max(3, ceil(3/epsilon)). Then
C_B2(f_n) >= (1 - epsilon) g(n) for **every** integer n >= N_epsilon. So
liminf C_B2(f_n)/(2n) >= 1, equivalently liminf C_B2(f_n)/n >= 2. The remainder
is the constant 6, which is o(n). The finite projection (c, p) = (2, 1)
follows by exact algebra (equal to the guarantee). Both coordinates come from
this one theorem about this one fixed family. The bound covers every integer
length, not a sparse subsequence.

## 6. Scope and checks

* No restricted basis, formula, monotone or average-case transfer is used.
  The argument is directly in full B2, with all 16 gate types, repeated pins,
  free constants and unrestricted sharing. All constant, unary (identity or
  negation) and output cases of eliminated gates are handled in Lemma 2.
  Negated output handling (case 4) adds no free basis element: it rewrites an
  existing gate's truth table.
* Sanity check against upper bounds: MOD3 has B2 circuits of size O(n) (an
  adder chain of mod-3 counters), so a coefficient-2 *lower* bound is not
  contradicted by any known construction. We claim nothing about optimality.
* The declared coefficient 2 is below the 31/10 literature seed. This
  submission offers a fully verified, self-contained (2, 1) result. It does not
  claim progress over the seed.
