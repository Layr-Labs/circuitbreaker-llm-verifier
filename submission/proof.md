# Coefficient 2 for the threshold-two family: a self-contained gate-elimination proof

This submission proves, with every step supplied here, that the threshold-two
family needs at least 2n-3 gates in the full B2 model at every input length
n >= 2. The argument is the classical gate-elimination method, written out
completely for this one family. No external result is used; the following
reference is attribution only and no step depends on it: C. P. Schnorr, "Zwei
lineare untere Schranken fuer die Komplexitaet Boolescher Funktionen",
Computing 13 (1974), 155-171.

## 1. The family and its NP membership

Let L be the set of binary strings that contain at least two 1-bits. For every
integer n >= 0 let f_n : {0,1}^n -> {0,1} be its characteristic function:
f_n(x) = 1 exactly when x has at least two coordinates equal to 1. In
particular f_0 (on the empty string) and f_1 are the constant 0 function. We
also write T2^n for f_n and OR^m for the m-variable disjunction.

A fixed algorithm A scans x once, counting 1-bits and stopping at two; it runs
in time O(n) at every length and outputs f_n(x). Define the verifier
V(x, w) = 1 iff w is the empty string and A(x) = 1, with witness bound W(n) = 0.
If x is in L then V(x, empty) = 1. If V(x, w) = 1 then w is empty and A(x) = 1,
so x is in L. Hence L is in P, and therefore in NP, with one fixed verifier and
witness bound for all lengths. Witnesses are not circuit inputs; the circuits
below read only x.

## 2. Circuit model and conventions

We use the Li-Yang general B2 model exactly as fixed in the rules. A circuit on
input variables x_1..x_n is an acyclic directed graph whose nodes are the n
input nodes, any number of constant nodes (each labelled 0 or 1), and gates.
Each gate has an ordered pair of pins; each pin is attached to an earlier node
(an input, a constant node or a gate), and both pins may be attached to the
same node. Each gate carries one of the 16 truth tables g : {0,1}^2 -> {0,1}.
One node is selected as the output. Fan-out and depth are unrestricted. Every
gate is counted once, including trivial or degenerate gates; inputs, constant
nodes, wires and the output selection cost zero. The size |C| is the number of
gates, and C_B2(f) is the least size of a circuit that computes f on every
input. A node computes a function of the inputs in the evident way.

A function h is a *literal* if h = x_k or h = NOT x_k for one input variable.

This model allows zero-cost constant nodes anywhere, which is at least as
permissive as the rules' model (constant nodes during substitutions): every
circuit of the rules' model is a circuit here of the same size, so a lower
bound proved here transfers to the rules' model. By Lemma 1(1)-(2) below, a
minimum circuit for a function that is neither constant nor a literal reads no
constant node, so the two minima coincide for such functions.

C_B2 does not depend on the names of the input variables: renaming inputs maps
circuits to circuits of the same size. Since T2 is symmetric, a circuit for
T2 on any n-1 named variables bounds C_B2(T2^(n-1)).

**Absorbing a negation (re-attachment).** If a pin of a gate H is attached to a
node u and we re-attach it to a node v with value(v) = NOT value(u), replacing
H's truth table h(a,b) by h(NOT a, b) (first pin) or h(a, NOT b) (second pin)
keeps the value computed by H unchanged. If both pins of H are re-attached this
way, both arguments are complemented. No gate is added.

**Complementing a gate in place.** Replacing the truth table h of a gate H by
NOT h makes H compute NOT value(H). Absorbing a negation into every pin that
reads H (complementing that argument of the reading gate's truth table, as
above) keeps every reading gate's value unchanged. No gate is added or removed.

**Acyclicity.** Every re-attachment below moves a pin from a gate G to a node u
that feeds a pin of G (or to a new constant node). Since u precedes G, and G
precedes all of its readers, u precedes every gate whose pin is moved, so the
graph stays acyclic.

## 3. Normalization of minimum circuits

**Lemma 1.** Let f be a function on n inputs that is neither constant nor a
literal, and let C be a circuit of minimum size computing f. Then:

1. the output node of C is a gate;
2. no pin of any gate of C is attached to a constant node;
3. no gate of C has both pins attached to the same node;
4. the truth table of every gate of C depends on both of its arguments.

*Proof.* (1) If the output were an input or a constant node, f would be a
literal or constant.

For (2)-(4), call a gate G *degenerate* if one of (2)-(4) fails at G. Then the
value of G equals one of: a constant c; the value of one node u; or the
negation of the value of one node u. Indeed, if a pin of G is attached to a
constant node, G is a function of its other pin alone; if both pins are
attached to the same node u, G is a function of u alone; if its truth table
ignores an argument, G is a function of the other pin's node alone; and every
function of one bit is 0, 1, the identity or the negation.

Suppose G is degenerate and is not the output. Delete G and re-attach every pin
that was attached to G: to a new constant node c in the first case; to u in the
second; to u with the negation absorbed into the receiving gate in the third
(Section 2). Every gate keeps its value, so the new circuit computes f with one
gate fewer, contradicting minimality.

Suppose G is degenerate and is the output; the output gate may also be read by
other gates. If G's value is a constant, f is constant. If it equals the value
of a node u: re-attach every pin attached to G to u, select u as the output and
delete G; every gate keeps its value and one gate fewer computes f, a
contradiction. If it equals NOT value(u): when u is an input, f is a literal;
when u is a constant node, f is constant; when u is a gate H, complement H in
place, absorb a negation into every pin attached to H other than G's pin
(Section 2), then re-attach every pin that was attached to G to H (H now has
the old value of G, so these readers keep their values), select H as the output
and delete G. All remaining gates keep their values and the new output computes
NOT value(H_old) = f with one gate fewer, a contradiction. Every case is
impossible, so no gate is degenerate. QED

Consequently a minimum circuit for such f contains no constant, unary or
negated-projection gate at all, and no gate takes a constant pin.

**Lemma 2.** In a minimum circuit for f as in Lemma 1 (with n >= 1 inputs and
at least one gate), the first gate G in a topological order has its two pins
attached to two distinct input variables x_i and x_j.

*Proof.* No gate precedes G, so its pins are attached to inputs or constant
nodes. By Lemma 1(2) they are inputs, and by Lemma 1(3) they are distinct. QED

## 4. Three distinct restrictions

For n >= 3, distinct indices i, j and bits a, b, write f|(a,b) for the function
of the other m = n-2 >= 1 variables obtained from T2^n by fixing x_i = a and
x_j = b. Counting 1-bits gives f|(0,0) = T2^m, f|(0,1) = f|(1,0) = OR^m and
f|(1,1) = 1.

**Lemma 3.** For m >= 1 the functions T2^m, OR^m and the constant 1 are pairwise
distinct.

*Proof.* At the all-zero assignment T2^m = OR^m = 0, while the constant is 1, so
the constant differs from both. At the assignment with exactly one 1-bit (it
exists since m >= 1) T2^m = 0 and OR^m = 1. QED

**Lemma 4.** For n >= 2, T2^n is neither constant nor a literal, and it
depends on every variable.

*Proof.* T2^n(0..0) = 0 and T2^n(1..1) = 1. For any k pick l != k (n >= 2);
at the assignment with x_l = 1 and all other bits 0, flipping x_k changes the
value from 0 to 1. A literal depends on one variable only, while T2^n depends
on n >= 2 variables. QED

## 5. One elimination step

**Lemma 5.** Let n >= 3 and let C be a minimum circuit for T2^n. Then some input
variable is attached to pins of at least two distinct gates.

*Proof.* By Lemma 4 the function is neither constant nor a literal, so Lemmas
1 and 2 apply (a circuit without gates outputs an input or constant). Let G be
the first gate, with pins on distinct inputs x_i and x_j. By Lemma 1(3) no
gate has both pins on x_i. Suppose that x_i and x_j are each attached to only
one pin, namely their pin of G. Fix any assignment z to the other n-2
variables. By induction along the topological order, every node other than x_i
and x_j has a value of the form F_u(z, G(x_i,x_j)): this holds for G itself
(F_G(z,g) = g), for the other inputs and the constants, and for every other
gate because neither of its pins is attached to x_i or x_j. By Lemma 1(1) the
output is a gate, not x_i or x_j, so the output value is F(z, G(x_i, x_j)). Hence f|(a,b) = f|(a',b') whenever
G(a,b) = G(a',b'), so at most two of the four restrictions f|(a,b) are
distinct. This contradicts Lemma 3, which gives three distinct restrictions
for m = n-2 >= 1. So one of x_i, x_j, say x_i, is also attached to a pin of a
gate G' other than G. Both G and G' have a pin on x_i. QED

**Lemma 6.** Let n >= 3. Then C_B2(T2^n) >= C_B2(T2^(n-1)) + 2.

*Proof.* Take a minimum circuit C for T2^n and the variable x_i of Lemma 5,
attached to distinct gates G and G'. Substitute x_i = 0: re-attach every pin on
x_i to a constant node 0 and remove x_i from the inputs. The resulting circuit
computes T2^n restricted by x_i = 0, which is T2^(n-1) on the other n-1
variables (zero adds no 1-bit). It has |C| gates.

Now G has a constant pin, so its value is a constant, the value of one node u,
or the negation of the value of one node u (Lemma 1's case analysis). We remove
it with the operations from the proof of Lemma 1, which are function-preserving
transformations of any circuit and do not rely on minimality. If G is not the
output, delete it and re-attach every pin that was attached to G to a new
constant node, to u, or to u with the negation absorbed into the receiving gate
(Section 2). If G is the output, the circuit computes T2^(n-1) with
n-1 >= 2, which by Lemma 4 is neither constant nor a literal; so G's value is
not a constant, and u is neither an input nor a constant node (otherwise the
output would be a constant or a literal). Hence u is a gate. In the identity
case, re-attach every pin attached to G to u, select u as the output and delete
G. In the negation case, complement u in place, absorb a negation into every
pin attached to u other than G's pin, re-attach every pin that was attached to
G to u (which now carries G's old value), select u as the output and delete G.
The output gate may have readers; in every case they are re-attached before G
is deleted, so no pin is left dangling. Each such removal deletes exactly one
gate and keeps the computed function.
Re-attaching pins, complementing truth-table arguments and complementing a gate
in place never remove a constant pin from another gate, so G' still has its constant pin from x_i (or
two constant pins) after this step, and it is still degenerate. Remove G' in
the same way. The final circuit computes T2^(n-1) with at most |C| - 2 gates.
Thus C_B2(T2^(n-1)) <= C_B2(T2^n) - 2. QED

## 6. The theorem and the declared pair

**Theorem.** For every integer n >= 2, C_B2(f_n) >= 2n - 3.

*Proof.* Induction on n. For n = 2, f_2 = x_1 AND x_2 is neither constant nor a
literal (Lemma 4), so a circuit with no gate, whose output is an input or a
constant node, cannot compute it; hence C_B2(f_2) >= 1 = 2*2-3. For n >= 3,
Lemma 6 and the induction hypothesis give
C_B2(f_n) >= C_B2(f_(n-1)) + 2 >= 2(n-1) - 3 + 2 = 2n - 3. QED

**Asymptotic guarantee.** Let g(n) = 2n, the polynomial_log expression with
coefficient 2, exponent 1, log exponent 0 and loglog exponent 0. For every
integer n >= 3 the factors (log2 n)^0 and (log2 log2 n)^0 equal 1 (log2 log2 n
is positive there), so g(n) = 2n, and C_B2(f_n)/g(n) >= (2n-3)/(2n) =
1 - 3/(2n). For each epsilon > 0 put N_epsilon = max(3, ceil(3/(2*epsilon)));
then C_B2(f_n) >= (1-epsilon)g(n) for every integer n >= N_epsilon. Hence liminf over all integers n of
C_B2(f_n)/g(n) >= 1, and equivalently liminf C_B2(f_n)/n >= 2. The declared
finite projection (2, 1) is this same statement with c = 2, p = 1. No
sparse subsequence, limsup or finite experiment is used: the bound holds at
every integer length n >= 2, and the family, verifier and witness bound are
fixed before the quantifiers.

**Orientation only.** A circuit that scans x_1..x_n while maintaining
s_k = x_1 OR .. OR x_k and t_k = [at least two 1-bits among x_1..x_k] by
t_k = t_(k-1) OR (s_(k-1) AND x_k) and s_k = s_(k-1) OR x_k shows
C_B2(f_n) <= 3n. This upper bound plays no role in the proof.

## 7. Coverage of the model conventions

- Degenerate, constant and unary gates count once; Lemma 1 shows a minimum
  circuit has none, so the lower bound applies to every circuit, including
  circuits that contain such gates (they are only larger).
- Constant nodes created by substitution cost zero (Lemma 6) and are never
  charged; a gate that reads one is either deleted by an explicit
  function-preserving step or retained and counted. No constant is created as
  a gate.
- Negations are absorbed into descendant truth tables with every receiving
  pin compensated; at the output, a negated gate is complemented in place, and
  a negated input or constant output contradicts Lemma 4. No uncharged
  computational basis or free output negation is assumed.
- Repeated pins and shared subgraphs are allowed by the model; Lemma 1(3)
  shows a minimum circuit has no repeated pin, and Lemma 5 accounts for every
  pin attached to x_i and x_j through the factoring argument.
- Inputs are exactly x; the length n is the actual input bit length; there is
  no padding, block encoding or advice.
