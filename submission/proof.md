# Coefficient 2 for MOD3 in the full B2 model: a self-contained gate-elimination proof

This file proves, with no external premises, that the counting family
MOD3 below satisfies C_B2(f_n) >= 2n - 3 for every integer n >= 2, hence
liminf C_B2(f_n)/n >= 2. The declared pair is (2, 1). The argument is the
classical gate-elimination method (Schnorr 1974; see also Wegener, *The
Complexity of Boolean Functions*, 1987, Ch. 5). Those works are named only
for provenance: every step used here is proved below, for this concrete
family, in the circuit model of rules.md. No new circuit result beyond this
classical bound is claimed; the bound is weaker than the 31/10 literature seed.

## 1. The family, at every length

For a binary string x let |x|_1 be its number of ones. For n >= 0 and a
residue a in Z_3 = {0,1,2} define

    M^n_a(x) = 1  if |x|_1 ≡ a (mod 3),   and 0 otherwise,   x in {0,1}^n.

The scored family is f_n = M^n_0 for EVERY integer n >= 0. In particular
f_0(empty string) = 1 (the empty count is 0), f_1(x) = 1 iff x = 0, and
f_2(x) = 1 iff x = 00. There are no exceptional lengths. The language is
L = { x in {0,1}* : |x|_1 ≡ 0 (mod 3) }, and f_n is its characteristic
function on {0,1}^n. The auxiliary functions M^n_1, M^n_2 are used only
inside the induction.

## 2. Uniform polynomial time, hence NP

Algorithm A: scan x once, keeping a counter c in Z_3 (start c = 0; on each
1 replace c by c+1 mod 3); output 1 iff c = 0 at the end. A runs in time
O(n) on input length n, uses no advice and is the same algorithm at every
length; it outputs f_{|x|}(x) by definition.

NP verifier: V(x, w) = 1 iff w is the empty string and A(x) = 1. Witness
bound W(n) = 0. Both directions: if x in L then w = empty gives
V(x, w) = A(x) = 1 with |w| = 0 <= W(|x|); conversely if V(x, w) = 1 for
some w with |w| <= W(|x|) then A(x) = 1, so x in L. V runs in polynomial
time. So L is in P, hence in NP, with one fixed verifier and witness bound.
Witnesses are not circuit inputs; the circuits below read only x.

## 3. Circuit model and conventions

We use the general B2 model of rules.md (Li–Yang, ECCC TR21-023, Sections
1.1, 2.1, 2.4, 3.3). A circuit on variables x_1, ..., x_n is a finite
directed acyclic graph whose nodes are: the n input nodes; any number of
constant nodes labelled 0 or 1; and gates. Each gate v has an ordered pair
of pins (p_1(v), p_2(v)), each pin being an earlier node (pins may coincide),
and an arbitrary truth table g_v : {0,1}^2 -> {0,1} (all 16 are allowed).
One node is designated the output. On input x every node gets a value in
the obvious way; the circuit computes h if the output value equals h(x) for
EVERY x in {0,1}^n. The size |C| is the number of gates, each counted once
whether or not it is degenerate; inputs, constant nodes, wires and output
selection cost zero. Fan-out, sharing and depth are unrestricted.
C_B2(h) is the minimum size of a circuit computing h; it is well defined
because every h has a circuit (e.g. a DNF built from binary gates). Renaming
the variables does not change sizes.

**Model transfer.** This model contains the Li–Yang B2 model of rules.md:
every circuit there is a circuit here with the same size (free constant
nodes are merely permitted here, never required). So the minimum here is at
most the minimum in rules.md, and every lower bound proved here holds in the
rules.md model. Section 4 below plays the role of the Li–Yang normalization
(Sections 2.4 and 3.3 of TR21-023): it covers constant, unary, repeated-pin
and terminal-output cases explicitly, and each removal only exhibits a
smaller circuit, so every gate of every circuit considered stays charged.

A function h depends on x_i if some two inputs differing only in x_i get
different h-values.

**Proper gates.** A gate v is *proper* if (i) its two pins are distinct
nodes, (ii) neither pin is a constant node, and (iii) g_v depends on both
of its arguments. A gate that is not proper is *degenerate*.

**Observation 3.1.** If v is degenerate then, on every input, value(v) =
u(value(w)) for one fixed node w that is a pin of v, or a constant node,
and one fixed u in {0, 1, id, neg}.
*Proof.* If the pins coincide (= w), value(v) = g_v(w, w), a unary function
of w. If a pin is a constant node κ, substitute κ into g_v; what remains is
a unary function of the other pin (if the other pin is also constant, take
w to be that constant node). If g_v does not depend on one argument, it is
a unary function of the other pin. Every unary function is 0, 1, id or neg. ∎

## 4. The removal lemma (normalization)

**Lemma 4.1 (removing one degenerate gate).** Let D be a circuit computing
a function h that depends on at least two variables, and let v be a
degenerate gate of D. Then there is a circuit D' with |D'| = |D| - 1
computing h, whose gates are exactly the gates of D other than v (possibly
with rewired pins and modified truth tables), such that every gate that was
degenerate in D is still degenerate in D'.

*Proof.* Write value(v) = u(value(w)) as in Observation 3.1; note w is not v
and w precedes v. Delete v and repair as follows.

1. *Successor pins.* For each gate z and each pin of z that equals v: if u
   is the constant κ, point the pin to a constant node κ; if u = id, point
   it to w; if u = neg, point it to w and replace g_z by the table obtained
   by negating that argument. Then z receives exactly the same argument
   value as before on every input, so all node values are unchanged.
   Acyclicity is kept because w precedes v, which precedes z.
2. *Output.* If v is not the output, nothing else is needed. If v is the
   output: u is not constant and w is not a constant node, since otherwise
   h would be constant. If u = id, make w the output. If u = neg, then w is not an input node or constant node,
   because then h would depend on at most one variable; so w is a gate.
   Replace g_w by 1 - g_w, and for every gate z with a pin on w, negate that
   argument in g_z (as in step 1). All nodes other than w keep their values,
   w now computes the old value of v, and w becomes the output.

So D' computes h and has one gate fewer. Degeneracy is preserved for every
other gate z: if z had a constant pin, that pin is untouched (the only
pins changed are those pointing to the gate v in step 1; step 2 changes
truth tables only) — in step 1 a pin may additionally become constant,
which keeps z degenerate; if z had two equal pins, both are
rewired identically, so they stay equal or both become constants (a
constant pin); if g_z ignored an argument, negating arguments or negating
the whole table does not create dependence on it. ∎

**Corollary 4.2 (normalization).** If h depends on at least two variables
and D computes h, then there is a circuit D* computing h with |D*| <= |D|
in which every gate is proper and every gate is an ancestor of the output
(or the output itself). Moreover, if D has k distinct degenerate gates, a
circuit computing h of size <= |D| - k exists.
*Proof.* While some gate is degenerate, apply Lemma 4.1 to it. Each
application removes exactly one gate, and every originally degenerate gate
stays degenerate until it is itself removed, so all k of them are removed
(at least k applications) and the process stops, since size decreases, with
every gate proper. Finally delete all gates that are not ancestors of the
output; this changes no output value, cannot make a gate degenerate (the
pins of kept gates are kept) and does not increase size. ∎

In D*: the output is a gate (h depends on two variables, so it is neither
an input node nor a constant node); no gate has a constant pin; constant
nodes, if any remain, are unused and free.

## 5. Facts about MOD3

**Lemma 5.1 (restriction).** Fixing x_i = b in M^n_a gives, on the other
n-1 variables (renamed in order), the function M^{n-1}_{a-b}.
*Proof.* The count of the remaining bits must be ≡ a - b. ∎

**Lemma 5.2 (full dependence).** For n >= 2 and every a, M^n_a depends on
every variable.
*Proof.* Since n >= 2, some t in {0, ..., n} has t ≡ a (mod 3). Take any x
with exactly t ones. Flipping x_i changes the count to t ± 1, which is not
≡ a, so the value changes from 1 to 0. ∎

**Lemma 5.3 (three subfunctions).** Let n >= 3, a in Z_3, i ≠ j. The four
restrictions of M^n_a obtained by fixing (x_i, x_j) in {0,1}^2 include at
least three pairwise distinct functions of the remaining m = n-2 >= 1
variables.
*Proof.* Fixing (x_i, x_j) = (b, c) yields M^m_{a-b-c}; the residues for
(0,0), (0,1), (1,1) are a, a-1, a-2, i.e. all three residues. The functions
M^m_0, M^m_1, M^m_2 are pairwise distinct for m >= 1: the all-zero input
(count 0) gives value 1 only to M^m_0, and an input with exactly one 1
(possible since m >= 1) gives value 1 only to M^m_1. ∎

## 6. The elimination lemma

Let S(n) = min over a in Z_3 of C_B2(M^n_a).

**Lemma 6.1.** For every n >= 3, S(n) >= S(n-1) + 2.

*Proof.* Fix a and a circuit C of size s computing M^n_a. By Lemma 5.2,
M^n_a depends on all n >= 3 variables, so Corollary 4.2 gives a circuit C*
of size s* <= s computing it, with all gates proper and the output a gate.
Choose a gate G that is minimal in a topological order, i.e. no pin of G is
a gate. Its pins are neither constant nodes nor equal (G is proper), so they
are two distinct input nodes x_i and x_j. For an input node, its *fan-out*
is the number of distinct gates having it as a pin.

*Case 1: x_i and x_j both have fan-out 1* (only G). Form D from C* by
deleting x_i, x_j and turning G into a fresh input node y. Since every edge
leaving x_i or x_j enters G, each remaining node of C* computes, on input x,
exactly the value its counterpart in D computes on (y, x_rest) with
y = g_G(x_i, x_j), where x_rest are the other n-2 variables. Hence
M^n_a(x) = F(g_G(x_i, x_j), x_rest) for the function F computed by D. Every
restriction of (x_i, x_j) to constants therefore yields F(0, ·) or F(1, ·):
at most two distinct functions. This contradicts Lemma 5.3. So Case 1 is
impossible.

*Case 2: one of them, say x_i, has fan-out >= 2.* Then besides G there is a
gate H ≠ G with a pin on x_i. Replace the input node x_i by a constant node
with value b = 0. The resulting circuit has the same
size s* and computes M^{n-1}_{a} on the remaining variables (Lemma 5.1).
In it, G and H are two distinct gates each having a constant pin, so both
are degenerate. The computed function depends on all n - 1 >= 2 variables
(Lemma 5.2), so Corollary 4.2 (with k = 2) gives a circuit for M^{n-1}_a of
size <= s* - 2 <= s - 2.

Hence C_B2(M^n_a) >= S(n-1) + 2 for every a, i.e. S(n) >= S(n-1) + 2. ∎

Note that all repeated-pin, constant-node, unary-gate and negated-output
situations, including ones created by the substitution, are absorbed by
Lemma 4.1, which charges every gate and adds no free basis element.

## 7. The bound at every length

**Base.** S(2) >= 1: M^2_a depends on two variables (Lemma 5.2), so no
circuit whose output is an input or constant node computes it; the output
is a gate.

**Induction.** By Lemma 6.1, S(n) >= 2n - 3 for every integer n >= 2.

**Theorem.** For every integer n >= 2,
C_B2(f_n) = C_B2(M^n_0) >= S(n) >= 2n - 3.

**Asymptotics.** For every integer n >= 2, C_B2(f_n)/n >= 2 - 3/n. Given
ε > 0, put N_ε = max(2, ceil(3/ε)); then C_B2(f_n) >= (2 - ε) n for EVERY
integer n >= N_ε. Thus liminf_{n→∞} C_B2(f_n)/n >= 2. The loss 3 is a
constant, so 3/n -> 0 and the coefficient is exactly 2, not reduced by a
linear loss. This holds at every sufficiently large integer length, not on
a subsequence, for the one fixed family f_n = M^n_0 with the fixed
evaluator A and verifier V of Section 2, measured against the actual input
length n with no padding or encoding.

This establishes the declared full guarantee
g(n) = 2 n^1 (log2 n)^0 (log2 log2 n)^0 with liminf C_B2(f_n)/g(n) >= 1,
and its finite projection (c, p) = (2, 1).

## 8. Scope

The coefficient 2 is below the 31/10 literature seed (Li–Yang). This
submission's value is a complete, citation-free proof that a reviewer can
check line by line. It makes no claim about upper bounds beyond what is
stated (MOD3 has linear-size circuits, so no superlinear claim is made).
