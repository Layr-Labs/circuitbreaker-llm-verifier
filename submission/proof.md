# Coefficient 2 for MOD_3: a self-contained gate-elimination proof

This file proves, without external premises, that the fixed family below has
C_B2(f_n) >= 2n-6 for every integer n >= 3, hence the joint pair (2,1). The
argument is the classical gate-elimination method (Schnorr, "Zwei lineare untere
Schranken fuer die Komplexitaet Boolescher Funktionen", Computing 13, 1974); it
is reproduced in full here so that nothing depends on the citation. No new
circuit lower bound is claimed: the coefficient 2 is weaker than the published
31/10 seed. Its purpose is a fully supplied, checkable proof.

## 1. The family and NP membership

Let L be the set of binary strings x (including the empty string) whose number
of ones is divisible by 3. For every integer n >= 0, f_n : {0,1}^n -> {0,1} is
the characteristic function of L at length n:

    f_n(x) = 1  iff  x_1 + ... + x_n == 0 (mod 3).

So f_0 = 1, f_1(x) = [x_1 = 0], and every length is defined by the same rule;
there are no exceptions, advice or per-length choices.

Evaluator A: scan x once, keeping a counter in {0,1,2} incremented mod 3 at each
1; accept iff the final counter is 0. This takes O(n) steps at every length.
Verifier V(x,w) = 1 iff w is the empty string and A(x) = 1; witness bound
W(n) = 0. If x in L, then w = empty has |w| <= W(|x|) and V(x,w) = 1. If some w
with |w| <= W(|x|) = 0 has V(x,w) = 1, then w is empty and A(x) = 1, so x in L.
Hence L in P, so L in NP. The witness is never a circuit input; the circuits
below have exactly the n input bits of x.

More generally, for m >= 0 and a residue s in Z_3 write

    M_{m,s}(y_1,...,y_m) = 1  iff  y_1 + ... + y_m == s (mod 3).

Thus f_n = M_{n,0}. Residues are always taken mod 3.

## 2. Circuit model

We use Li–Yang's general B2 model (ECCC TR21-023, 2021-02-20, Sections 1.1,
2.1, 2.4 with footnotes 3–4, and 3.3), restated here in full.
A circuit is a finite acyclic directed graph whose nodes are: the input nodes
x_1..x_m, constant nodes 0 and 1, and gates. Each gate G has two ordered pins;
each pin is connected to one earlier node (input, constant or gate), and G
carries an arbitrary truth table T_G : {0,1}^2 -> {0,1} (all 16 allowed).
Both pins may be connected to the same node; fan-out and depth are
unrestricted. One node is designated as the output. The circuit computes the
function v_out(x), where node values are defined in topological order. The
size |C| is the number of gates, every gate counted once, including trivial,
degenerate or negated-projection gates; inputs, constants, wires and the
output designation cost zero. C(h) = C_B2(h) is the least size of a circuit
computing h exactly on all 2^m inputs. Circuits are nonuniform. A minimum
circuit exists for every h since sizes are nonnegative integers and some
circuit (e.g. a DNF built from binary gates) computes h.

The lower bound is proved for a minimum circuit, so it holds for every circuit.
All simplifications below are applied to circuits we construct inside the
proof; they never add a gate, never add a basis element and never treat a
constant node or negation as free computation beyond what the model already
allows (constant nodes cost zero, and a truth table may absorb negation of a
pin because all 16 tables are allowed).

We say h depends essentially on y_i if some two inputs differing only in
y_i get different values.

## 3. Bypassing a degenerate gate

**Lemma 1 (bypass).** Let C compute h and let G be a gate of C. Suppose that
either (i) v_G is a constant b for all inputs, or (ii) v_G = u(v_w) for all
inputs, where w is a node connected to a pin of G and u : {0,1} -> {0,1} is
the identity or negation. Assume moreover that G is not the output, or that
G is the output in case (ii) with w a gate. Then there is a circuit computing h
whose gates are exactly the gates of C other than G (with possibly changed
truth tables and pin connections); in particular its size is |C| - 1.

*Proof.* Step 1: for every pin of every gate S connected to G: in case (i)
reconnect it to the constant node b; in case (ii) reconnect it to w and
replace T_S by the table obtained by applying u to that pin's argument. Since
w is a pin node of G, w precedes G, which precedes S, so the graph stays
acyclic. Every gate other than G sees the same pin values as before, so all
other node values are unchanged. Step 2: if G was the output (so case (ii)
and w is a gate): if u is the identity, designate w as output. If u is
negation, replace T_w by its complement 1 - T_w, then for every pin connected
to w (including pins reconnected in Step 1) compose the receiving table with
negation on that pin, and designate w as output. Each receiving gate again
sees exactly its previous value, and the new v_w equals the old v_G. Step 3:
delete G, which now has no outgoing connections and is not the output. The
output value is unchanged on every input. No gate was added. QED.

**Lemma 2 (normal form).** Let h depend essentially on at least two
variables and let C be a minimum circuit for h. Then the output is a gate, and
for every gate G of C: its two pins are connected to two distinct nodes,
neither is a constant node, and T_G(a,b) depends essentially on both a and b.

*Proof.* If the output were an input or constant node, h would depend on at
most one variable. Suppose some gate G violates the conclusion. If a pin of
G is at a constant node, or both pins are at the same node, or T_G depends on
at most one argument, then v_G is a constant or v_G = u(v_w) for a pin node w
of G and u in {identity, negation}. If G is not the output, Lemma 1 gives a
smaller circuit for h, contradicting minimality. If G is the output and w is
a gate (case (ii)), Lemma 1 again contradicts minimality. Otherwise G is the
output and v_G is a constant, or a literal of an input variable, or a function
of a constant node; then h depends on at most one variable, contrary to
hypothesis. QED.

## 4. Properties of M_{m,s}

**Lemma 3.**
(a) If m >= 2, each M_{m,s} takes value 1 somewhere, and for s != s' the
functions M_{m,s} and M_{m,s'} are different.
(b) If m >= 3, M_{m,s} depends essentially on every variable.
(c) Fixing y_i = c in M_{m,s} gives M_{m-1,s-c} on the other m-1 variables.
(d) If m >= 4 and i != j, the four restrictions of M_{m,s} obtained by
fixing (y_i,y_j) = (a,b) in {0,1}^2 include at least three distinct
functions of the remaining m-2 variables.

*Proof.* (a) With m >= 2 we can choose exactly 0, 1 or 2 ones, so every
residue is attained; the supports of M_{m,s} and M_{m,s'} are disjoint and
nonempty, so the functions differ. (b) Fix i. By (a) on the other m-1 >= 2
variables choose z with sum == s. Setting y_i = 0 and the rest to z gives 1;
setting y_i = 1 gives sum == s+1, value 0. (c) is immediate. (d) The
restriction at (a,b) is M_{m-2, s-a-b}. The assignments (0,0), (0,1), (1,1)
give residues s, s-1, s-2, which are pairwise distinct mod 3; since
m-2 >= 2, (a) shows these are three distinct functions. QED.

## 5. The elimination step

**Lemma 4.** For every m >= 4 and every residue s,

    C(M_{m,s}) >= 2 + min over s' of C(M_{m-1,s'}).

*Proof.* Let C be a minimum circuit for h = M_{m,s}, with |C| = Gamma. By
Lemma 3(b) h depends essentially on m >= 3 >= 2 variables, so Lemma 2 applies.
The output is a gate, so C has a gate. Let G be a gate that is first among the
gates in a topological order. Its pins are connected to nodes that are not
gates; by Lemma 2 they are not constant nodes and are distinct. Hence they are
two input variables x_i and x_j with i != j. Call the fan-out of a variable
the number of gate pins connected to it.

*G is not the output.* Otherwise h = T_G(x_i,x_j) depends on at most two
variables, contradicting Lemma 3(b) with m >= 3.

*x_i or x_j has fan-out at least 2.* Suppose both have fan-out exactly 1, so
the only pins at x_i and x_j are the two pins of G. Replace G by a fresh input
node z (keeping all connections out of G as connections out of z). The
resulting circuit computes some Phi(z, x_rest), where x_rest is the other m-2
variables, and since no other gate reads x_i or x_j and the output is not x_i
or x_j, h(x) = Phi(T_G(x_i,x_j), x_rest) for every input. If T_G(a,b) =
T_G(a',b') then the restrictions of h at (a,b) and at (a',b') are both
Phi(T_G(a,b), .), hence equal. As T_G takes at most two values, the four
restrictions comprise at most two distinct functions, contradicting
Lemma 3(d) (m >= 4).

By symmetry let x_i have fan-out at least 2. G has exactly one pin at x_i
(its pins are at distinct nodes), so some gate H != G has a pin at x_i.

*Substitution.* Reconnect every pin attached to x_i to the constant node 0
and remove x_i. The resulting circuit C_0 has the same Gamma gates and
computes, by Lemma 3(c), h' = M_{m-1,s} on the other m-1 variables. By
Lemma 3(b) with m-1 >= 3, h' depends essentially on at least three variables.

In C_0, G has a pin at a constant, so v_G = T_G(0, v_{x_j}) is a constant or
u(v_{x_j}) with u the identity or negation. G is not the output, so Lemma 1
yields C_1 computing h' with gates those of C_0 except G. Lemma 1 changes only
pins that were connected to G, so H still has a pin at the constant node 0
(H's other pin may have been reconnected from G to x_j or to a constant).
Hence in C_1, v_H is a constant or u(v_w) where w is the node at H's other pin.
If H is not the output, Lemma 1 removes H. If H is the output: if v_H is
constant, or w is an input or constant node, then h' would depend on at most
one variable, impossible; so w is a gate and Lemma 1 removes H. Either way we
obtain a circuit for M_{m-1,s} with Gamma - 2 gates. Therefore
min over s' of C(M_{m-1,s'}) <= C(M_{m-1,s}) <= Gamma - 2. QED.

The degenerate situations are all covered: repeated pins and constant pins
are excluded in C by Lemma 2 and created only by the substitution, where they
are handled by Lemma 1; H may read G (then it is reconnected and still has a
constant pin); H or G as output is treated explicitly; constant nodes cost
zero throughout; and no gate is ever added.

## 6. The main theorem

Let a_m = min over s of C(M_{m,s}). Lemma 4 applies to every s, so
a_m >= a_{m-1} + 2 for every m >= 4, and a_3 >= 0 since sizes are
nonnegative. By induction a_m >= 2(m-3) for every m >= 3.

**Theorem.** For every integer n >= 3, C_B2(f_n) = C_B2(M_{n,0}) >= a_n >= 2n - 6.

**Asymptotic form.** For n >= 3, C_B2(f_n)/n >= 2 - 6/n. Given epsilon > 0, let
N_epsilon = max(3, ceil(6/epsilon)); then for EVERY integer n >= N_epsilon,
C_B2(f_n) >= (2 - epsilon) n. In the declared form with g(n) = 2n: for every
epsilon > 0 and every integer n >= max(3, ceil(3/epsilon)),
C_B2(f_n) >= 2n - 6 >= (1 - epsilon) * g(n). Equivalently liminf over all integers n of
C_B2(f_n)/n >= 2, i.e. C_B2(f_n) >= 2n - o(n) with remainder r(n) = 6,
r(n)/n -> 0. This is exactly the declared guarantee g(n) = 2 * n^1 *
(log2 n)^0 * (log2 log2 n)^0 with finite projection (c,p) = (2,1). The bound
holds at every integer length n >= 3, not on a subsequence, and the family,
evaluator, verifier and witness bound are fixed before any quantifier.

The coefficient and exponent belong to this one theorem about this one family.

## 7. Scope and remarks

- The bound is measured against the actual input length n; there is no
  padding, blocking or encoding loss.
- Only a lower bound is claimed; nothing here asserts tightness.
- The two supplementary string claims in submission.json are the explicit
  inequality of Section 6 and the step inequality of Lemma 4, both proved above.
- The proof is complete as written; the Schnorr citation is attribution only.
