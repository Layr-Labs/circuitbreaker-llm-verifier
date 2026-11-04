# Coefficient 7/3 for MOD3 in the full B2 model: weighted gate elimination plus a degree bound

This file proves, with no external premises, that the counting family MOD3
below satisfies C_B2(f_n) >= (7/3)n - 4 for every integer n >= 2, hence
liminf C_B2(f_n)/n >= 7/3. The declared pair is (7/3, 1).

The argument combines two classical ideas, each proved in full below:
(i) gate elimination (Schnorr 1974; Wegener, *The Complexity of Boolean
Functions*, 1987, Ch. 5) run on a *weighted* measure that charges 3 for an
XOR-type gate and 2 for every other gate; (ii) the GF(2)-degree bound on the
number of non-affine gates (Schnorr's multiplicative-complexity bound). The
combination of a weighted measure with a multiplicative-complexity bound is
the idea of Kojevnikov and Kulikov (CiE 2010). These works are named only for
provenance: every step used here is proved below for this concrete family in
the circuit model of rules.md. No result beyond this classical 7n/3 - O(1)
bound is claimed; it is weaker than the 31/10 literature seed.

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
rules.md model. Section 4 plays the role of the Li–Yang normalization
(Sections 2.4 and 3.3 of TR21-023): it covers constant, unary, repeated-pin
and terminal-output cases explicitly, and each removal only exhibits a
circuit of smaller measure, so every gate of every circuit considered stays
charged.

A function h depends on x_i if some two inputs differing only in x_i get
different h-values. For an input node x_i, its *fan-out* is the number of
distinct gates having x_i as a pin.

**Table classes and weights.** Of the 16 tables, the two tables
(y_1, y_2) -> y_1 ⊕ y_2 and (y_1, y_2) -> 1 ⊕ y_1 ⊕ y_2 form the *XOR class*.
The *affine* tables are the 8 tables (y_1, y_2) -> b_1 y_1 ⊕ b_2 y_2 ⊕ b_0;
the other 8 tables are *non-affine*. Every XOR-class table is affine. The
weight of a gate v is w(v) = 3 if g_v is in the XOR class and w(v) = 2
otherwise. For a circuit C let X(C) be the number of XOR-class gates,
A(C) = |C| - X(C) the number of all other gates (including degenerate ones),
N(C) the number of non-affine gates (so N(C) <= A(C)), and

    μ(C) = 3 X(C) + 2 A(C) = Σ_v w(v).

Negating one argument of a table, or negating its output, maps the XOR class
onto itself; since these operations are involutions they also map the
complement of the XOR class onto itself. Hence such operations never change
a gate's weight.

**Proper gates.** A gate v is *proper* if (i) its two pins are distinct
nodes, (ii) neither pin is a constant node, and (iii) g_v depends on both
of its arguments. A gate that is not proper is *degenerate*.

**Observation 3.1.** If v is degenerate then, on every input, value(v) =
u(value(w)) for one fixed node w that is a pin of v, or a constant node,
and one fixed u in {0, 1, id, neg}. If a pin of v is a constant node κ,
then u is the unary function y -> g_v(κ, y) (or g_v(y, κ)) of the other pin.
*Proof.* If the pins coincide (= w), value(v) = g_v(w, w), a unary function
of w. If a pin is a constant node κ, substitute κ into g_v; what remains is
a unary function of the other pin (if the other pin is also constant, take
w to be that constant node). If g_v does not depend on one argument, it is
a unary function of the other pin. Every unary function is 0, 1, id or neg. ∎

**Observation 3.2 (proper tables).** A table depending on both arguments
is either in the XOR class or of *AND type*, i.e. of the form
((y_1 ⊕ α) ∧ (y_2 ⊕ β)) ⊕ γ. *Proof.* Exactly 10 of the 16 tables depend on
both arguments (16 minus 2 constants minus 4 tables y_1, ¬y_1, y_2, ¬y_2);
2 are in the XOR class and the 8 AND-type tables (α, β, γ in {0,1}) are
pairwise distinct and depend on both arguments. ∎
An AND-type table has a *controlling value* for each argument: putting
y_1 = α (resp. y_2 = β) makes it the constant γ.

## 4. The removal lemma (normalization)

**Lemma 4.1 (removing one degenerate gate).** Let D be a circuit computing
a function h that depends on at least two variables, and let v be a
degenerate gate of D. Then there is a circuit D' computing h, whose gates
are exactly the gates of D other than v (possibly with rewired pins and
truth tables modified only by negating arguments or outputs), such that
every gate that was degenerate in D is still degenerate in D'. In
particular |D'| = |D| - 1 and μ(D') = μ(D) - w(v) <= μ(D) - 2. Moreover,
if u (from Observation 3.1 for v) is a constant, then every gate that had a
pin on v has a constant-node pin in D', hence is degenerate in D'.

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
   h would be constant. If u = id, make w the output. If u = neg, then w is
   not an input node or constant node, because then h would depend on at
   most one variable; so w is a gate. Replace g_w by 1 - g_w, and for every
   gate z with a pin on w, negate that argument in g_z (as in step 1). All
   nodes other than w keep their values, w now computes the old value of v,
   and w becomes the output.

So D' computes h and has one gate fewer. Tables change only by negating an
argument or the output, so by Section 3 every remaining gate keeps its
weight and μ drops by exactly w(v) >= 2. Degeneracy is preserved for every
other gate z: if z had a constant pin, that pin is untouched (the only pins
changed are those pointing to v in step 1; step 2 changes truth tables
only) — in step 1 a pin may additionally become constant, which keeps z
degenerate; if z had two equal pins, both are rewired identically, so they
stay equal or both become constants; if g_z ignored an argument, negating
arguments or negating the whole table does not create dependence on it.
The last sentence of the lemma is step 1 with u constant. ∎

**Corollary 4.2 (normalization).** Let h depend on at least two variables
and let D compute h. If Z is a set of distinct gates of D that are all
degenerate in D, then some circuit D* computing h has
μ(D*) <= μ(D) - Σ_{v in Z} w(v). Also some circuit D° computing h with
μ(D°) <= μ(D) has every gate proper and every gate an ancestor of the
output (or the output itself).
*Proof.* Apply Lemma 4.1 to the members of Z one at a time; each stays
degenerate until it is itself removed, and each removal lowers μ by that
gate's weight. For D°, apply Lemma 4.1 while some gate is degenerate (μ
strictly decreases, so this stops), then delete all gates that are not
ancestors of the output; this changes no output value, keeps the pins of
kept gates, and does not increase μ. ∎

## 5. Facts about MOD3

**Lemma 5.1 (restriction).** Fixing x_i = b in M^n_a gives, on the other
n-1 variables (renamed in order), the function M^{n-1}_{a-b}.
*Proof.* The count of the remaining bits must be ≡ a - b. ∎

**Lemma 5.2 (full dependence).** For n >= 2 and every a, M^n_a depends on
every variable; in particular it is not constant.
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

**Lemma 5.4 (bottleneck).** Let n >= 3, let C be a circuit on n variables,
let x_i ≠ x_j be input nodes, and let B be a set of gates containing a gate
H such that: every pin of a gate in B is x_i, x_j or a gate of B; no gate
outside B has a pin on x_i, on x_j, or on a gate of B \ {H}; and the output
is not in B \ {H}. Then C does not compute M^n_a for any a.
*Proof.* Suppose C computes M^n_a. Every node outside B ∪ {x_i, x_j}
computes, on input x, a value determined by value(H) and the other n-2
variables x_rest (induction along the topological order: its pins are
outside B ∪ {x_i, x_j} or equal H). The output is H or such a node: it is
not in B \ {H} by hypothesis, and not x_i or x_j because M^n_a depends on
n >= 3 variables (Lemma 5.2). Hence M^n_a(x) = F(value(H), x_rest) for a
fixed F, so the four restrictions of (x_i, x_j) yield at most the two
functions F(0, ·), F(1, ·), contradicting Lemma 5.3. ∎

Lemma 5.4 is used below with B = {G} (x_i, x_j each of fan-out one, H = G)
and with B = {G, H} (G has pins x_i, x_j and its only successor is H, whose
pins lie in {x_i, x_j, G}, and no other gate reads x_i or x_j).

## 6. The weighted elimination lemma

Let W(n) = min over a in Z_3 and over circuits C computing M^n_a of μ(C).

**Lemma 6.1.** For every n >= 3, W(n) >= W(n-1) + 6.

*Proof.* Fix a and a circuit C computing M^n_a with μ(C) = W(n); we show
W(n-1) <= μ(C) - 6. M^n_a depends on all n >= 3 variables (Lemma 5.2), so by
Corollary 4.2 and minimality we may assume every gate of C is proper and an
ancestor of the output. The output is a gate depending on n >= 3 variables,
so no gate whose pins are both input nodes is the output; any non-output
gate, being an ancestor of the output, has at least one successor gate.

Throughout, "fix x_k = b" means: replace the input node x_k by a constant
node b. The resulting circuit has the same gates and the same μ, and on the
remaining n-1 >= 2 variables it computes M^{n-1}_{a-b} (Lemma 5.1), which
depends on at least two variables and is not constant (Lemma 5.2). So
Corollary 4.2 applies to it, and it suffices in each case to exhibit gates
whose successive removal lowers μ by at least 6; the result is a circuit for
a member of the family at length n-1, whence W(n-1) <= μ(C) - 6.

*Case 0: some input x_k has fan-out >= 3.* Fix x_k = 0. Three distinct
gates have a constant pin, so they are degenerate; each has weight >= 2.
Corollary 4.2 lowers μ by >= 6.

From now on every input has fan-out <= 2. Choose a gate G minimal in a
topological order. Its pins are not gates, not constant nodes and not
equal (G is proper), so they are two distinct input nodes x_i, x_j; and G is
not the output. If x_i and x_j both had fan-out 1, Lemma 5.4 with B = {G}
would give a contradiction. So after renaming, x_i has fan-out exactly 2:
its consumers are G and one other gate H ≠ G. Let L denote the other
consumer of x_j if x_j has fan-out 2.

*Case 1: G is of AND type* (Observation 3.2), with controlling values α
for its x_i-pin and β for its x_j-pin.
(1a) G has a successor K ≠ H. Fix x_i = α. G and H have a constant pin.
Remove G first: u is constant (Observation 3.1), so by Lemma 4.1 K gets a
constant pin; G, H remain/are degenerate as well. G, H, K are distinct
(K ≠ G since K has a gate pin; K ≠ H by choice). Three removals of weight
>= 2 each lower μ by >= 6.
(1b) The only successor of G is H. Since H has pins x_i and G, which are
distinct, H's pins are exactly {x_i, G}. If x_j has fan-out 1, then
B = {G, H} satisfies Lemma 5.4 (x_i is read only by G and H, x_j only by G,
G only by H), a contradiction. So L exists; L has pin x_j, and x_j is not a
pin of H, so L ≠ H, and L ≠ G. Fix x_j = β. G is constant, so removing it
gives its only successor H a constant pin; L has a constant pin. G, H, L
are distinct degenerate gates: μ drops by >= 6.

*Case 2: G is in the XOR class* (weight 3).
(2a) H is in the XOR class. Fix x_i = 0. G and H both have a constant pin;
both have weight 3: μ drops by 6.
(2b) H is of AND type, with controlling value α for its x_i-pin. Fix
x_i = α. Then H computes the constant of its table, so u is constant for H.
The restricted function is not constant, so H is not the output; H is an
ancestor of the output, so it has a successor gate M. Remove H first; by
Lemma 4.1 M gets a constant pin. G also has a constant pin. G, H, M are
distinct (M ≠ G because M has the gate pin H while G's pins are inputs;
M ≠ H by acyclicity). μ drops by >= 3 + 2 + 2 >= 6.

By Observation 3.2 the cases are exhaustive (G and H are proper). In every
case a circuit for some M^{n-1}_{a'} has μ <= μ(C) - 6. ∎

**Corollary 6.2.** W(2) >= 2 (M^2_a depends on two variables, so the output
is a gate, of weight >= 2), and by Lemma 6.1, W(n) >= 6n - 10 for every
integer n >= 2.

## 7. The degree bound on non-affine gates

Every h : {0,1}^m -> {0,1} has a unique representation as a multilinear
polynomial over GF(2), h(x) = ⊕_S c_S Π_{i in S} x_i (its algebraic normal
form); deg h is the largest |S| with c_S = 1 (deg 0 for constants). Uniqueness:
the 2^m monomial functions are linearly independent, since evaluating at
the indicator 1_T gives h(1_T) = ⊕_{S ⊆ T} c_S, which determines c_T by
induction on |T|. Inverting, c_S = ⊕_{T ⊆ S} h(1_T), because
⊕_{T ⊆ S} ⊕_{U ⊆ T} c_U = ⊕_U c_U · #{T : U ⊆ T ⊆ S} and the count
2^{|S|-|U|} is odd only for U = S. Reducing any polynomial with x^2 = x
does not raise degree, so deg(p ⊕ q) <= max(deg p, deg q) and
deg(p · q) <= deg p + deg q.

**Lemma 7.1 (hyperplane).** Let h : {0,1}^m -> {0,1}, let
ℓ(x) = x_k ⊕ λ(x_{-k}) ⊕ c_0 with λ a GF(2)-linear form in the other
variables, and for b in {0,1} let h_b(x_{-k}) = h(x) evaluated at
x_k = λ(x_{-k}) ⊕ c_0 ⊕ b. Then deg h <= 1 + max(deg h_0, deg h_1).
*Proof.* For every x, ℓ(x) = b exactly when x_k = λ ⊕ c_0 ⊕ b, so
h(x) = h_{ℓ(x)}(x_{-k}) = h_0(x_{-k}) ⊕ ℓ(x) · (h_0 ⊕ h_1)(x_{-k}). The right
side is a polynomial of degree <= max(deg h_0, 1 + max(deg h_0, deg h_1)). ∎

**Lemma 7.2 (non-affine gates).** If a circuit C (model of Section 3, any
number of variables) computes h, then deg h <= N(C) + 1.
*Proof.* Induction on N(C). If N(C) = 0, every node computes an affine
function of x (inputs and constants do; an affine table applied to affine
functions is affine), so deg h <= 1. If N(C) >= 1, let g be the first
non-affine gate in topological order. Every earlier node computes an affine
function, so the pins of g compute affine functions ℓ_1, ℓ_2.
If ℓ_1 is a constant function c, replace the table of g by the table
(y_1, y_2) -> g_g(c, y_2), a unary and hence affine function: all values
are unchanged and N drops by one, so deg h <= N(C) by induction.
Otherwise ℓ_1 = x_k ⊕ λ(x_{-k}) ⊕ c_0 for some variable x_k. For b in {0,1}
build C_b on the variables x_{-k}: add new XOR-class gates (or a constant
node, if λ = 0) computing λ(x_{-k}) ⊕ c_0 ⊕ b from the input nodes, placed
before all old gates, and point every former pin on x_k to that node.
Acyclicity holds and every old node now takes the value it had in C at
x_k = λ ⊕ c_0 ⊕ b, so C_b computes h_b of Lemma 7.1, and pin 1 of g now
computes the constant b. Replace g's table as in the constant case. The new
gates are affine and no other table changes, so N(C_b) <= N(C) - 1 and by
induction deg h_b <= N(C). Lemma 7.1 gives deg h <= N(C) + 1. ∎

**Lemma 7.3 (degree of MOD3).** For n >= 2 and every a, deg M^n_a >= n - 1.
*Proof.* For S of size k, c_S = ⊕_{T ⊆ S} M^n_a(1_T) = N_k(a) mod 2, where
N_k(r) = #{T ⊆ {1..k} : |T| ≡ r (mod 3)}. Claim: for k >= 1, N_k(r) is
even iff r ≡ -k (mod 3). For k = 1, (N_1(0), N_1(1), N_1(2)) = (1, 1, 0) and
-1 ≡ 2. Splitting by whether k+1 is in T gives
N_{k+1}(r) = N_k(r) + N_k(r-1). Exactly one residue makes N_k even, so
N_{k+1}(r) is even iff both N_k(r) and N_k(r-1) are odd, i.e. iff
r ≢ -k and r - 1 ≢ -k, i.e. iff r ≡ -k-1. The residues -n and -(n-1) differ,
so N_n(a) or N_{n-1}(a) is odd: the monomial x_1⋯x_n, or x_1⋯x_{n-1}
(with n - 1 >= 1), has coefficient 1. ∎

## 8. The bound at every length

**Theorem.** For every integer n >= 2, C_B2(f_n) >= (7/3)n - 4.
*Proof.* Let C be any circuit computing f_n = M^n_0. By Corollary 6.2,
3X(C) + 2A(C) = μ(C) >= W(n) >= 6n - 10. By Lemmas 7.2 and 7.3,
A(C) >= N(C) >= deg f_n - 1 >= n - 2. Adding,
3|C| = 3X(C) + 3A(C) >= (6n - 10) + (n - 2) = 7n - 12. ∎

**Asymptotics.** For every integer n >= 2, C_B2(f_n)/n >= 7/3 - 4/n. Given
ε > 0, put N_ε = max(2, ceil(4/ε)); then C_B2(f_n) >= (7/3 - ε) n for EVERY
integer n >= N_ε. Thus liminf_{n→∞} C_B2(f_n)/n >= 7/3. The loss 4 is a
constant, so the coefficient is exactly 7/3, not reduced by a linear loss.
This holds at every sufficiently large integer length, not on a
subsequence, for the one fixed family f_n = M^n_0 with the fixed evaluator A
and verifier V of Section 2, measured against the actual input length n
with no padding or encoding. By the model transfer of Section 3 it holds
in the rules.md full B2 model.

This establishes the declared full guarantee
g(n) = (7/3) n^1 (log2 n)^0 (log2 log2 n)^0 with liminf C_B2(f_n)/g(n) >= 1,
and its finite projection (c, p) = (7/3, 1).

## 9. Scope

The coefficient 7/3 is below the 31/10 literature seed (Li–Yang). This
submission's value is a complete, citation-free proof that a reviewer can
check line by line. It makes no claim about upper bounds (MOD3 has
linear-size circuits, so no superlinear claim is made).
