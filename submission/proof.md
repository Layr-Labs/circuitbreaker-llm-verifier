# Coefficient 2, exponent 1: a self-contained gate-elimination proof for MOD3

This submission proves, with every step written out, the classical lower bound
C_B2(f_n) >= 2n - 3 (Schnorr 1974) for one fixed total family in P. It applies
known mathematics and claims no new circuit lower bound. It is weaker than the
published 31/10 seed; its purpose is a fully supplied proof that relies on no
unreproduced external premise. Section numbers below refer to this file.

Notation: for a bit string x, |x|_1 is the number of ones in x. Indices of
residues are taken mod 3. For m >= 0 and a in {0,1,2} define

    M_a^m : {0,1}^m -> {0,1},   M_a^m(x) = 1  iff  |x|_1 = a (mod 3).

## 1. The fixed family and its NP membership

The language is L = { x in {0,1}* : |x|_1 = 0 (mod 3) }, with characteristic
functions f_n = M_0^n on {0,1}^n for every integer n >= 0. Every length is
covered with no exceptions: f_0 is the constant 1 on the empty string (it has
zero ones), f_1(x) = NOT x_1, and so on. Small lengths play no role in the
asymptotic claim but are defined.

Evaluator A: scan x once, keeping a counter s in {0,1,2}, initially 0; on
reading a 1 replace s by s+1 mod 3. Output 1 iff s = 0 at the end. By induction
on the number of bits read, s = (number of ones read) mod 3, so A(x) = 1 iff
x in L. A runs in O(n) time on a RAM and O(n log n) on a Turing machine; in
any case polynomial time. The same algorithm is used at every length, with no
advice and no length-dependent choice.

NP verifier: V(x, w) = 1 iff w is the empty string and A(x) = 1; witness
bound W(n) = 0. V is polynomial time. If x in L then w = empty has |w| <= 0 and
V(x,w) = 1. Conversely, if V(x,w) = 1 for some w with |w| <= W(|x|), then
A(x) = 1, so x in L. Hence L in P, so L in NP. Witnesses are not circuit inputs;
the circuits below read only the n bits of x.

## 2. Circuit model

A B2 circuit on variables x_1..x_m is a finite directed acyclic graph whose
source nodes are the variables; every other node (a gate) has exactly two
ordered input pins, each pointing to an earlier node (both pins may point to
the same node), and is labelled by one of the 16 functions {0,1}^2 -> {0,1}.
One node is designated the output. Fan-out and depth are unrestricted. The
size is the number of gates; every gate is counted, including degenerate
ones (constant tables, projections, negated projections, repeated pins).
Inputs, wires and the choice of output cost nothing. The circuit computes F if
the output's value equals F(x) on every x in {0,1}^m. C_B2(F) is the minimum
size. This is exactly the Li-Yang general B2 model in rules.md.

For the proof it is convenient to use a *generalized circuit*: the same
definition, except that the source nodes may additionally include constant
nodes labelled 0 and 1, which cost nothing. Let C*(F) be the minimum size of a
generalized circuit computing F. Every B2 circuit is a generalized circuit,
so

    C_B2(F) >= C*(F).                                          (2.1)

Hence any lower bound proved for C* holds in the pinned model. Allowing free
constants only weakens the circuits we must defeat; nothing uncharged is added
to the circuits being counted in C_B2.

A *literal* is a function of the form x_k or NOT x_k for a single variable.

## 3. Elimination and substitution lemmas

Call a gate G *degenerate* if one of: (i) one of its pins points to a constant
node; (ii) both pins point to the same node; (iii) its truth table, as a
function of its two pins, depends on at most one argument. In each case the
value of G on every input equals one of: a constant b; the value of a single
node u that is a pin of G; or the negation of the value of such a node u.
(In case (i) fix the constant argument; in case (ii) substitute equal
arguments; in case (iii) drop the irrelevant argument. In all cases G becomes
a function of one bit, and the four functions of one bit are 0, 1, identity,
negation.) Note u precedes G in the topological order.

**Lemma 3.1 (elimination).** Let D be a generalized circuit of size s
computing F, and let G be a degenerate gate of D. If F is neither a constant
nor a literal, there is a generalized circuit of size s - 1 computing F,
obtained by deleting G and changing only (a) pins that point to G, (b) truth
tables of gates, and (c) the output designation. In particular every pin that
pointed to a constant node in D still points to that constant node.

*Proof.* Three cases.

* G equals the constant b. Redirect every pin that points to G to the constant
  node b (create it if absent; it is free). If G is the output, designate the
  constant node b as the output. Delete G.
* G equals the value of node u. Redirect every pin that points to G to u; if
  G is the output, designate u as output. Delete G.
* G equals NOT u. For every gate H with a pin pointing to G, redirect that pin
  to u and replace H's truth table h by the table that negates that argument
  (if both pins of H pointed to G, negate both arguments). H's value on every
  input is unchanged. If G is not the output, delete G. If G is the output:
  if u is a constant node, designate the opposite constant node as output; if
  u is a gate, replace u's truth table by its negation, then for every gate
  whose pin points to u (now including the redirected former successors of G)
  negate that argument in its table, so all their values are unchanged, and
  designate u, which now computes NOT(old u) = old G = F, as output; delete G.
  If u is a variable, F = NOT x_k would be a literal, excluded by hypothesis.

In every case each redirected pin points to u or a constant, which precede
every former successor of G, so the graph stays acyclic. Exactly one gate is
removed and every remaining node computes the same function as before (or, for
the re-labelled u, its complement, compensated in all successors), so the
output still computes F. Pins pointing to constant nodes are never redirected,
because only pins pointing to G are redirected. QED.

**Corollary 3.2 (normal form).** If F is neither constant nor a literal and D
is a generalized circuit computing F of size C*(F), then D has no degenerate
gate. In particular every gate of D has two pins pointing to two distinct
non-constant nodes. Moreover the output of D is a gate (a source would compute a
constant or a variable).

*Proof.* A degenerate gate would, by Lemma 3.1, give a smaller circuit. QED.

**Lemma 3.3 (substitution).** If D is a generalized circuit of size s
computing F(x_1..x_m), and c in {0,1}, then replacing the source x_k by the
constant node c (redirecting all pins pointing to x_k, and the output if it is
x_k) gives a generalized circuit of size s on the remaining m - 1 variables
computing the restriction F|_{x_k=c}. Every gate that had a pin pointing to
x_k now has a pin pointing to a constant node, so it is degenerate.

*Proof.* Immediate: on each input of the remaining variables, every node
takes the value it took in D on the extended input with x_k = c. QED.

## 4. Facts about M_a^m

**F1.** For m >= 2 and every a, M_a^m depends on every one of its m
variables; in particular it is neither constant nor a literal.
*Proof.* M_a^m is symmetric, with value v(w) = [w = a mod 3] at weight w. The
weights 0,1,2 all occur for m >= 2 and exactly one of them is = a mod 3, so v
is not constant on {0,...,m}; hence there is w with 0 <= w < m and
v(w) != v(w+1). For any variable x_k take an input with x_k = 0 and weight w
(possible since w <= m - 1); flipping x_k gives weight w + 1 and changes the
value. A constant depends on no variable and a literal on exactly one, while
m >= 2. QED.

**F2.** Fixing one variable x_k = c turns M_a^m into M_{a-c}^{m-1} on the
remaining variables (indices mod 3), since the remaining weight must be
= a - c (mod 3). Fixing two variables to (b1, b2) gives M_{a-b1-b2}^{m-2}.

**F3.** For m' >= 1 the three functions M_0^{m'}, M_1^{m'}, M_2^{m'} are
pairwise distinct.
*Proof.* Each input has exactly one residue of its weight, so their supports
are pairwise disjoint; two of them can be equal only if both are identically
zero. M_0^{m'} is 1 at 0^{m'} and M_1^{m'} is 1 at 10^{m'-1}, so neither is
zero; hence M_0 != M_1, and M_2 differs from both (if M_2 = M_0, both would be
zero). QED.

## 5. Main theorem

**Theorem 5.1.** For every integer m >= 2 and every a in {0,1,2},
C*(M_a^m) >= 2m - 3. Consequently, by (2.1), C_B2(f_n) = C_B2(M_0^n) >= 2n - 3
for every integer n >= 2.

*Proof.* Induction on m, simultaneously for all a.

*Base m = 2.* 2m - 3 = 1. A generalized circuit with no gates designates a
source as output and so computes a constant or a variable, which M_a^2 is not
by F1. Hence at least one gate.

*Step m >= 3.* Assume the claim for m - 1 (and all residues). Let F = M_a^m
and let D be a generalized circuit computing F with s = C*(F) gates. By F1,
F is neither constant nor literal, so by Corollary 3.2 D has no degenerate
gate and its output is a gate; in particular s >= 1. Let G be the first gate
in a topological order of D. Its pins point to sources; by Corollary 3.2 they
are two distinct non-constant sources, i.e. two distinct variables x_i, x_j.

*Case B: neither x_i nor x_j is a pin of any gate other than G.* The output is
a gate, so it is neither x_i nor x_j. Fix any assignment rho to the other
m - 2 variables. For each node v other than x_i, x_j, by induction along the
topological order its value under (rho, x_i = b1, x_j = b2) depends on
(b1, b2) only through the value g(b1, b2) of G: v = G, or every pin of v
points to G or to a node other than x_i, x_j that precedes v. Hence the
restrictions F|_{x_i=b1, x_j=b2}, as functions of the other m - 2 variables,
coincide whenever g(b1, b2) coincides, and the four assignments yield at most
two distinct restrictions. But by F2 the assignments (0,0), (0,1), (1,0),
(1,1) yield M_a^{m-2}, M_{a-1}^{m-2}, M_{a-1}^{m-2}, M_{a-2}^{m-2}, and by F3
(m - 2 >= 1) the three functions M_a, M_{a-1}, M_{a-2} are pairwise distinct.
Contradiction; Case B cannot occur.

*Case A: some x_k in {x_i, x_j} is a pin of a gate H != G.* Apply Lemma 3.3
with x_k = 0. The result D' has s gates and computes F' = M_a^{m-1} (F2), on
m - 1 >= 2 variables, so F' is neither constant nor a literal (F1). In D' both
G and H have a pin pointing to the constant node 0, so both are degenerate.
Apply Lemma 3.1 to G; by its last sentence H still has its pin pointing to the
constant node 0 (that pin pointed to a constant, not to G), so H is still
degenerate in the new circuit (re-labelling H's truth table in the third case
of Lemma 3.1 does not change its pins). Apply Lemma 3.1 to H. The result is a
generalized circuit with s - 2 gates computing M_a^{m-1}. By the induction
hypothesis s - 2 >= 2(m - 1) - 3, i.e. s >= 2m - 3.

Since Case B is impossible, Case A holds and C*(M_a^m) = s >= 2m - 3. QED.

## 6. The declared joint theorem

For the single family f_n = M_0^n of Section 1, Theorem 5.1 gives, for EVERY
integer n >= 2,

    C_B2(f_n) >= 2n - 3 = 2n - r(n),  r(n) = 3,  r(n)/n -> 0.

Take g(n) = 2 * n^1 * (log2 n)^0 * (log2 log2 n)^0 = 2n, which matches the
declared guarantee c = 2, p = 1, q = 0, r = 0 (kind polynomial_log). For every
epsilon > 0 let N_epsilon = max(2, ceil(3/(2 epsilon))). For every integer
n >= N_epsilon, C_B2(f_n) >= 2n - 3 >= 2n - 2 epsilon n = (1 - epsilon) g(n).
Thus liminf over all integers n of C_B2(f_n)/g(n) >= 1, i.e. C_B2(f_n) >=
2n - o(n) at every sufficiently large length, not a sparse subsequence. The
declared finite projection (c, p) = (2, 1) is the guarantee itself with
q = r = 0, so liminf C_B2(f_n)/n >= 2 holds for the same family. Both
coordinates come from this one theorem about this one family.

No upper bound, optimality or tightness is claimed; the family and remainder
are fixed and stated, n is the actual input length, and there are no blocks,
padding or encodings to account for.

## 7. Scope, sources and relation to the seed

* Source of the argument: C. P. Schnorr, "Zwei lineare untere Schranken fuer
  die Komplexitaet Boolescher Funktionen", Computing 13 (1974), 155-171, which
  proves 2n - 3 for functions in the class Q_{2,3}^n, MOD3 among them. The
  argument is reproduced completely above, in the form presented in standard
  textbooks (e.g. Wegener, The Complexity of Boolean Functions, 1987, Ch. 5,
  Sec. 2; Jukna, Boolean Function Complexity, 2012, Sec. 16.1). No step of
  the proof depends on the citation.
* This is weaker than the published 31/10 seed (Li-Yang, ECCC TR21-023). It is
  submitted because every premise is supplied here, whereas the seed relies on
  a long external case analysis that this submission does not reproduce.
* No computation, LP certificate, conjecture, restricted basis, direct-sum
  assumption or prior submission is used.
