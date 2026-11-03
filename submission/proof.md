# Coefficient 2 for the exactly-one family in the full B2 model: a self-contained gate-elimination proof

This file proves, with no external premises, that the exactly-one counting
family below satisfies C_B2(f_n) >= 2n - 3 for every integer n >= 2, hence
liminf C_B2(f_n)/n >= 2. The declared pair is (c, p) = (2, 1). The argument
is the classical gate-elimination method (Schnorr 1974; see also Wegener,
*The Complexity of Boolean Functions*, 1987, Ch. 5). Those works are named
only for provenance: every step used here is proved below, for this concrete
family, in the circuit model of rules.md. No new circuit result is claimed;
the bound is strictly weaker than the 31/10 literature seed.

## 1. The family, at every length

For a binary string x let |x|_1 be its number of ones. The scored family is

    f_n(x) = 1  if |x|_1 = 1,   and 0 otherwise,

defined for EVERY integer n >= 0 and x in {0,1}^n. In particular
f_0(empty string) = 0 (the empty string has zero ones), f_1(x) = x_1, and
f_2 = XOR. There are no exceptional lengths. The language is
L = { x in {0,1}* : |x|_1 = 1 }, and f_n is its characteristic function on
{0,1}^n. Auxiliary functions on fewer variables appear only inside the
induction and are defined here when first used.

## 2. Uniform polynomial time, hence NP

Evaluator A: scan x once, maintaining the saturating count
c = min(number of ones seen so far, 2) in the constant-size state space
{0, 1, 2}; on each input bit equal to 1 replace c by min(c + 1, 2), and at
the end output 1 iff c = 1 (correct, since |x|_1 = 1 iff the saturated count
is 1). Each step updates a constant-size state, so A performs O(n) bit
operations in total — linear time — uses no advice, and is the same single
algorithm at every length; by definition it outputs f_{|x|}(x).

NP verifier: V(x, w) = 1 iff w is the empty string and A(x) = 1, with witness
bound W(n) = 0. Both directions are immediate: if x in L then the empty
witness is accepted with |w| = 0 <= W(|x|); conversely V(x, w) = 1 forces
A(x) = 1, so x in L. V runs in polynomial time. So L is in P, hence in NP,
with one fixed verifier and witness bound. Witness bits are not inputs of the
target circuits; the circuits below read only x.

## 3. Circuit model and conventions

A circuit on variables x_1, ..., x_n is a finite directed acyclic graph whose
nodes are: the n input nodes; any number of constant nodes labelled 0 or 1;
and gates. Each gate v has an ordered pair of pins, each pin an earlier node
(pins may coincide), and an arbitrary truth table
g_v : {0,1}^2 -> {0,1} (all 16 are allowed). One node is designated the
output. Every node receives a value on each input x in the evident way, and
the circuit computes h if the output value equals h(x) for EVERY x in
{0,1}^n. The size |C| is the number of gates, each counted once whether or
not it is degenerate; inputs, constant nodes, wires and the output selection
cost zero. Fan-out, sharing and depth are unrestricted. C_B2(h) is the
minimum size of a circuit computing h; it is well defined because every h has
a circuit of binary gates (for example a DNF). Renaming variables preserves
sizes.

**Model transfer.** This model contains the full Li-Yang B2 model of
rules.md: every circuit of that model is a circuit here of the same size, and
the free constant nodes are merely permitted here, never required. Hence the
minimum here is at most the minimum there, and every lower bound proved here
holds in the rules.md model. Section 4 below plays the role of the Li-Yang
normalization (TR21-023 Sections 2.4 and 3.3): constant, unary, repeated-pin
and terminal-output cases are all handled explicitly, and each removal only
exhibits a strictly smaller circuit, so every gate of every circuit
considered stays charged and no free extra computational basis appears.

A function h depends on x_i if some two inputs differing only in coordinate
i have different h-values.

**Proper gates.** A gate v is *proper* if (i) its two pins are distinct
nodes, (ii) neither pin is a constant node, and (iii) g_v depends on both of
its arguments. A gate that is not proper is *degenerate*.

**Observation 3.1.** If v is degenerate then, on every input, value(v) =
u(value(w)) for one fixed node w that is a pin of v and one fixed u in
{0, 1, id, neg}.
*Proof.* If the two pins coincide, say both are w, then value(v) =
g_v(w, w), a unary function of w. If a pin is a constant node, substituting
that constant into g_v leaves a unary function of the other pin. If g_v does
not depend on one of its arguments, value(v) is a unary function of the
other pin. Every unary Boolean function is one of 0, 1, id, neg, so in every
case w may be chosen as a pin of v (when both pins are constants, take u to
be the resulting constant and w either pin). ∎

## 4. The removal lemma (normalization)

**Lemma 4.1 (removing one degenerate gate).** Let D be a circuit computing a
function h that depends on at least two variables, and let v be a degenerate
gate of D. Then there is a circuit D' with |D'| = |D| - 1 computing h, whose
gates are exactly the gates of D other than v (with possibly rewired pins
and modified truth tables), such that every gate that was degenerate in D is
still degenerate in D'.

*Proof.* Write value(v) = u(value(w)) as in Observation 3.1; the node w
precedes v and w is not v. Delete v and repair the circuit.

1. *Successor pins.* For each gate z and each pin of z equal to v: if u is
   the constant 0 or 1, point that pin to a constant node with value u; if
   u = id, point it to w; if u = neg, point it to w and replace g_z by the
   truth table obtained by negating that one argument. After this, z
   receives exactly the same argument value as before on every input, so all
   remaining node values are unchanged. Acyclicity is preserved because w
   precedes v, which precedes z.
2. *Output.* If v is not the output, nothing more is needed. If v is the
   output, then h = u(value(w)) is nonconstant, so u is id or neg and w is
   not a constant node; and w is not an input node either, because h would
   then depend on at most one variable, contradicting the hypothesis. So w
   is a gate. If u = id, designate w as the output. If u = neg, replace g_w
   by 1 - g_w, and for every gate z with a pin on w negate the argument of
   g_z corresponding to that pin (as in step 1); then designate w as the
   output. All nodes other than w keep their values, readers of w see their
   old argument values back (double negation), w now computes the old value
   of v, and the output value equals h as before.

So D' computes h with one gate fewer. Degeneracy of every other gate z is
preserved: in step 1 the only pins that change are pins pointing to v, which
are redirected to w or to a constant node — a gate with a constant pin keeps
a constant pin or gains one; a gate whose two pins were both v keeps two
equal pins; and negating an argument of a table, or a whole table, cannot
create dependence on an argument the table ignored. In step 2 only truth
tables change, by such negations. ∎

**Corollary 4.2 (normalization).** If h depends on at least two variables
and D computes h, then there is a circuit D* computing h with |D*| <= |D| in
which every gate is proper and every gate is an ancestor of the output or
the output itself. Moreover, if D contains k distinct degenerate gates, then
some circuit with at most |D| - k gates computes h.
*Proof.* The function h is fixed throughout, so Lemma 4.1's hypothesis holds
at every step. While some gate is degenerate, apply Lemma 4.1 to it. Each
application deletes exactly one gate, and every originally degenerate gate
stays degenerate until it is itself deleted, so at least k applications
occur; since the size strictly decreases, the process terminates with all
gates proper. Finally delete every gate that is not an ancestor of the
output: this changes no computed value, keeps all pins of kept gates (the
pins of an ancestor are ancestors), so cannot make a proper gate degenerate,
and does not increase the size. ∎

In D*: since h depends on at least two variables, the output is neither an
input node nor a constant node, hence the output is a gate; no gate has a
constant pin; any constant nodes left over are unused and free.

## 5. Facts about the exactly-one family

For m >= 0 let NOR_m : {0,1}^m -> {0,1} be 1 iff the input is the all-zero
string (so NOR_0(empty) = 1).

**Lemma 5.1 (full dependence).** For every n >= 2, f_n depends on every one
of its n variables.
*Proof.* Fix i and let e_i in {0,1}^n be the string with a single 1 in
position i. Then f_n(e_i) = 1 and f_n(0^n) = 0, and e_i and 0^n differ only
in coordinate i. ∎

**Lemma 5.2 (one-variable restriction).** For n >= 1 and every i: fixing
x_i = 0 in f_n yields f_{n-1} on the remaining variables (renamed in
order); fixing x_i = 1 yields NOR_{n-1}.
*Proof.* With x_i = 0 the output is 1 iff the remaining n-1 bits contain
exactly one 1. With x_i = 1 the count is already 1, so the output is 1 iff
all remaining bits are 0. ∎

**Lemma 5.3 (three distinct subfunctions).** Let n >= 3 and i != j, and let
m = n - 2 >= 1. Among the four restrictions of f_n obtained by fixing
(x_i, x_j) to the values in {0,1}^2 there are at least three pairwise
distinct functions of the remaining m variables, namely f_m, NOR_m and the
constant 0.
*Proof.* The restriction (0,0) leaves exactly-one on the remaining bits:
f_m. The restrictions (0,1) and (1,0) force all remaining bits to be 0:
NOR_m. The restriction (1,1) makes the count 2, so the output is constantly
0. They are pairwise distinct for m >= 1: at input e_1 (one 1) f_m = 1 and
NOR_m = 0; at input 0^m f_m = 0 and NOR_m = 1; the constant 0 differs from
f_m at e_1 and from NOR_m at 0^m. ∎

## 6. The elimination lemma

Let S(n) = C_B2(f_n).

**Lemma 6.1.** For every n >= 3, S(n) >= S(n-1) + 2.

*Proof.* Let C be any circuit of size s computing f_n. Since f_n depends on
n >= 3 >= 2 variables (Lemma 5.1), Corollary 4.2 yields a circuit C* of size
s* <= s, all of whose gates are proper and reach the output, with the output
a gate. For an input node, its *fan-out* is the number of distinct gates
having it as a pin. Choose a gate G that is minimal in a topological order of
the gates, i.e. no pin of G is a gate. Since G is proper, its two pins are
distinct, non-constant and not gates: two distinct input nodes x_i and x_j.

*Case 1: every gate reading x_i or x_j is G itself* (both fan-outs equal
1). Form D from C* by deleting the nodes x_i, x_j and the gate G, adding a
fresh input node y, pointing every former reader of G to y, and designating
y as the output of D if G was the output of C*. Since no
node other than G reads x_i or x_j, each remaining node of C* computes, on
input x, exactly the value that its counterpart in D computes on
(y, x_rest), where x_rest lists the other n - 2 variables in order and
y = g_G(x_i, x_j). Hence

    f_n(x) = F(g_G(x_i, x_j), x_rest)

for the function F that D computes. Fixing (x_i, x_j) to any of the four
constant pairs therefore produces F(0, ·) or F(1, ·), at most two distinct
functions of the remaining n - 2 variables. This contradicts Lemma 5.3,
which supplies at least three. So Case 1 is impossible.

*Case 2: some gate other than G reads x_i or x_j.* Say x_i has a reader
H != G. Replace the input node x_i by a constant node with value 0. The
circuit keeps its size s* and computes f_{n-1} on the remaining variables
(Lemma 5.2). In it, G and H are two distinct gates, and each now has a
constant pin, so both are degenerate. The computed function f_{n-1} depends
on all n - 1 >= 2 variables (Lemma 5.1), so Corollary 4.2 (with k = 2) gives
a circuit with at most s* - 2 <= s - 2 gates computing f_{n-1}. Hence
S(n-1) <= s - 2.

In both surviving cases S(n-1) <= s - 2; minimizing over all C computing f_n
gives S(n) >= S(n-1) + 2. ∎

All repeated-pin, constant-node, unary-table and negated-output situations,
including the ones created by the substitution in Case 2, are absorbed by
Lemma 4.1, which charges every gate and introduces no uncounted basis
element.

## 7. The bound at every length

**Base.** S(2) >= 1: f_2 = XOR depends on two variables, so its output node
can be neither an input node nor a constant node; it is a gate. (One XOR
gate computes f_2, so S(2) = 1 exactly.)

**Induction.** By Lemma 6.1, S(n) >= 2n - 3 for every integer n >= 2.

**Theorem.** For every integer n >= 2, C_B2(f_n) >= 2n - 3.

**Asymptotics.** For every integer n >= 2, C_B2(f_n)/n >= 2 - 3/n. Given
epsilon > 0, put N_epsilon = max(2, ceil(3/epsilon)); then
C_B2(f_n) >= (2 - epsilon) n for EVERY integer n >= N_epsilon. In the
multiplicative form of the declared guarantee, g(n) = 2n, so
N'_epsilon = max(2, ceil(3/(2 epsilon))) gives
C_B2(f_n) >= (1 - epsilon) g(n) for every integer n >= N'_epsilon, since
n >= 3/(2 epsilon) implies 3/n <= 2 epsilon and hence
2n - 3 >= 2(1 - epsilon)n. Thus
liminf_{n -> infinity} C_B2(f_n)/n >= 2. The loss is the constant 3, and 3/n
tends to 0, so no linear loss reduces the coefficient below 2. This holds at
every sufficiently large integer length — not a subsequence, no padding, no
special remainder form — for the one fixed family f_n with the fixed
evaluator and verifier of Section 2, measured against the actual input
length n = |x|.

This establishes the declared full guarantee
g(n) = 2 n^1 (log2 n)^0 (log2 log2 n)^0 with liminf C_B2(f_n)/g(n) >= 1 and
its finite projection (c, p) = (2, 1).

## 8. Scope and provenance

The method is classical gate elimination (Schnorr 1974; Wegener 1987, Ch. 5),
cited for provenance only: no step above uses any external premise, and the
whole argument can be checked line by line. No new circuit result is
claimed: the coefficient 2 is below the 31/10 literature seed (Li-Yang, ECCC
TR21-023) and also below stronger results submitted or pending elsewhere on
this challenge. As a side remark, the elimination template of Section 6
applies verbatim to any total uniform family {h_n} that (i) for every
n >= 2 depends on all of its variables, (ii) for every n >= 3 and every
pair i != j has at least three pairwise distinct two-variable restrictions
on (x_i, x_j), and (iii) for every n >= 3 and every i has some constant b
with h_n restricted by x_i := b equal to h_{n-1} on the remaining
variables; all three properties are used by the induction (Lemma 6.1 needs
(ii) in Case 1 and (iii) with b = 0 in Case 2). The exactly-one family
satisfies (iii) with b = 0, and its restriction triple
(f_m, NOR_m, 0) is immediate to verify.
