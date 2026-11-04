# Case analysis for Theorem 6.1 — Part B: Cases 6–8

Notation and standing setup as in `proof.md` §§3–5 and `cases-a.md`.
After Case 5, no ∧-type gate is directly fed by two variables. Let G be a
topologically minimal ∧-type gate fed by I_1 and I_2, each computed by a
fair cyclic xor-circuit (or a variable), so xor-reconstruction applies.

## Case 6: substituting a constant to a protected variable

Let x be a protected variable; by Cases 1–2 it is a 1-variable feeding an
⊕-type gate P. Substitute x ← d and normalize; at least P is eliminated by
Rule 3. The nominal outcome — one gate eliminated, x non-influential, x's
quadratic equation killed, no potential change — gives:

    Δμ = 1 + α_I + α_Q.

Case 6 rules out every failure of this nominal outcome, leaving the
Case-6 step defined below, reused as the first substitution of several
Case-8 leaves. Note 1 + α_I + α_Q = 12.4 < δ at the optimum:
alone it would not suffice.

**Case 6.1: a second gate Q is eliminated during the normalization —
realizes 2 − 2α_φ + α_Q.** Q is not eliminated by Rule 1: if P is the output
gate, Q becomes the output gate after P's removal, which Rule 1 cannot
remove; otherwise Q's out-degree never decreases when P is removed by Rule 3,
so Rule 1 cannot apply to Q either. Two gates are eliminated (P, Q),
ΔΦ ≤ 2 by the normalization lemma, Δi = −1 (x), Δq = −1 (x's quadratic
equation). Hence Δμ ≥ 2 − 2α_φ + α_Q + α_I, the slack term.

For the rest of Case 6, P is the ONLY gate eliminated and ΔΦ = 1; we show
this is avoidable. Let R be P's other input besides x. To produce the
potential increment, R is an ∧-type gate or a variable.

**Case 6.2.1: R is an ∧-type gate, troubled after normalization (vacuous).**
R receives no new inputs from Rule 3 and its feeding variables' out-degrees
do not decrease, so R was already troubled — no increment. If R were newly
troubled, its out-degree must decrease, forcing P to be the output gate
and R a 2-gate; but then R's other descendant cannot reach the output gate
P (else a cycle containing an ∧-type gate), contradicting Case 0.2. So
ΔΦ = 1 is impossible here.

**Case 6.2.2: R = u is a variable.** Its out-degree does not decrease
during Rule 3, so u can feed a troubled gate after normalization only as
a "2−-variable" in the paper's sense: a 1-variable feeding P whose
out-degree rises to 2 through P's Rule-3 elimination (P's descendant
becomes fed by u), so u then feeds two gates, one possibly troubled.
Subcases split on u's protection status before the substitution.

*6.2.2.1: u is an unprotected 1-variable.* Let y be the couple of x.
Substitute a constant to y (making x unprotected), then the affine
substitution x ← u, trivializing P; t = 2, Δi = −3 (x, y, u),
Δμ ≥ (3/2)α_I per substitution.

*6.2.2.2: u is a protected 1-variable coupled to x.* Apply the affine
substitution x ← u, trivializing P; t = 1, Δi = −2, Δμ ≥ 2α_I.

*6.2.2.3: u is a protected 1-variable with couple v ≠ x — realizes α_I/3.*
Let y be the couple of x. Substitute constants to y and v (making both x
and u unprotected), then the affine substitution x ← u, trivializing P;
t = 3, Δi = −4 (x, y, u, v). Hence Δμ ≥ (4/3)α_I per substitution,
attaining the α_I/3 term.

*6.2.2.4: u is a 2-variable.* Then P must be a 1-gate. Let A be the other
descendant of u and B the descendant of P; A ≠ B (else P is useless). As
in 6.2.1, A cannot be the newly troubled gate, so B is: B is an ∧-type
gate fed by P and a variable t. Substitute a constant to t, trivializing B
(then P can be removed by Rule 1), then a constant to the couple of x
(making x unprotected); t = 2, Δi = −3 (x, t, couple of x),
Δμ ≥ (3/2)α_I per substitution.

*Standing outcome of Case 6:* apart from the configurations above,
substituting a constant to a protected variable eliminates exactly one
gate, introduces no troubled gate, and gives Δμ = 1 + α_I + α_Q. This is
the **Case-6 step** used below.

## Case 7: I_1 is a 2-variable, I_2 a 1-gate — realizes 4 − 4α_φ

Let I_1 = x be a 2-variable and I_2 = Q a 1-gate, both feeding the minimal
∧-type gate G. Then x is unprotected (Case 1) and x does not feed Q (else
Q is useless). Substitute a constant to x, trivializing G; t = 1. Four
gates are eliminated (G, Q, two descendants of G and x), ΔΦ ≤ 4 itemized
as: Q by Rule 1 (contributes ≤ 2); G trivialized, neither of its inputs
becoming or feeding a troubled gate (contributes 0); the other two by
Rule 3 (≤ 1 each). With Δi = −1 (x):

    Δμ ≥ 4 − 4α_φ + α_I,

attaining the 4 − 4α_φ term.

## Case 8: both inputs of G are gates, or a 2-gate remains

By symmetry let Q := I_2 be a gate (the variable/variable and
2-variable/1-gate configurations are Cases 5 and 7), and P := I_1 a gate
or a variable, with out-degree ≤ that of Q if it is a gate. By minimality
of G, Q computes an affine function of its inputs.

### Case 8.1: Q depends on an unprotected variable

By Proposition 2.7 (path property of fair xor-circuits), a path exists
from some unprotected variable x to Q. Substitute an appropriate constant
to Q through x via xor-reconstruction (t = 1): x becomes non-influential
and G is trivialized.

**Case 8.1.1: Q is a 2+-gate — realizes 4 − 4α_φ.** The reconstruction
replaces the constant output b of Q, which feeds at least three gates:
Q is removed and one new gate Z is added; normalization then removes the
three descendants of b and at least one descendant of G (Rules 2 and 3) —
4 gates net. The reconstruction changes no potential (the empty packing
persists), and the removals give ΔΦ ≤ 4; Δi = −1 (x). Hence
Δμ ≥ 4 − 4α_φ + α_I.

**Case 8.1.2: Q is a 1-gate and P = t is a variable.** Perform
xor-reconstruction on Q; both the unprotected 1-variable t and x become
non-influential (Δi = −2), so Δμ ≥ 2α_I.

**Case 8.1.3: P and Q are both 1-gates — realizes 4 − 4α_φ.** Q computes an
affine function. Xor-reconstruction replaces Q's output by a constant b of
out-degree 2; normalization removes G, P, and the other two descendants of
b and G — four distinct gates. Distinctness: b's two descendants differ
from G and from each other (a coincidence closes a cycle through the cyclic
xor-circuit), and P is not among them — fed by b, P would, unwinding the
reconstruction, feed Q in the original circuit, excluded since P, Q are the
two inputs of the minimal gate G, neither feeding the other. Itemized: P removed by Rule 1 contributes
≤ 2; G contributes 0; the other two descendants contribute ≤ 2 in total;
so ΔΦ ≤ 4. With Δi = −1 (x): Δμ ≥ 4 − 4α_φ + α_I.

*Standing condition for the rest of Case 8:* Q — and P as well, if it is an
⊕-type gate of equal out-degree — depends only on protected variables.

### Case 8.2: Q depends only on protected variables

Let I be the set of protected variables Q depends on.

**Case 8.2.1: two coupled protected variables feed Q.** Suppose Q is fed
by two protected variables x_j, x_k coupled to each other; both are
1-variables (Case 1). Apply the affine substitution x_j ← x_k ⊕ c, making
Q's output constant; t = 1, Δi = −2 (x_j, x_k), so Δμ ≥ 2α_I.

### Case 8.2.2: Q is a 2+-gate

Fix x_j ∈ I and let x_k be its couple. First substitution: x_k ← d plus
normalization; by Case 6 exactly one gate (a descendant of x_k) is
eliminated, with Δμ = 1 + α_I + α_Q, and x_j becomes unprotected.

**Case 8.2.2.1: Q is not fed by x_k — realizes (5 − 2α_φ + α_Q)/2.**
After x_k ← d, Q remains a 2+-gate and now depends on the unprotected
variable x_j. Second substitution: xor-reconstruction with x_j,
trivializing G. The reconstruction degenerates three gates A, B, C and
trivializes G: 4 gates eliminated in this step, with ΔΦ ≤ 2 (A contributes
≤ 1, C contributes ≤ 1; that eliminating B and G is ΔΦ-free is assumed here
and discharged in Case 8.2.5), so the second step gives
Δμ ≥ 4 − 2α_φ + α_I. Averaging the two substitutions:

    Δμ ≥ (1/2)[ (1 + α_I + α_Q) + (4 − 2α_φ + α_I) ]
       = (1/2)(5 − 2α_φ + 2α_I + α_Q) per substitution,

attaining (5 − 2α_φ + α_Q)/2: the sub-δ Case-6 step (12.4) is compensated
by the stronger second step (13.2).

**Case 8.2.2.2: Q is fed by x_k (reduction).** After x_k ← d, the gate A
that was fed by x_k is now fed by Q's other input Q′; Q′ is a 2+-gate
depending on the unprotected x_j. Apply Case 8.2.2.1 with Q := Q′; the
same bound holds.

*Standing condition after 8.2.2:* Q is a 1-gate, and P is a variable or a
1-gate.

### Case 8.2.3: P = t is a variable

Let x_j ∈ I with couple x_k. Substitute a constant to x_k (making x_j
unprotected) and normalize — G is then fed by t and by Q′, where Q′
computes an affine function of x_j and may equal Q — then apply
xor-reconstruction on Q′ trivializing G; t = 2. Three influential inputs
die (Δi = −3: t and x_j become non-influential through the reconstruction,
x_k was substituted), so Δμ ≥ (3/2)α_I per substitution.

### Case 8.2.4: P is an ⊕-type 1-gate

If P depended on an unprotected variable, Case 8.1.3 would apply with P, Q
exchanged; so both P and Q depend only on protected variables.

**Case 8.2.4.1: some x_j ∈ I has its couple x_k not feeding the other gate
— main branch realizes (5 − 2α_φ + α_Q)/2.** By symmetry treat Q and x_j,
with x_k not feeding P. First substitution x_k ← d: a Case-6 step
(Δμ = 1 + α_I + α_Q), and x_j becomes unprotected. Second substitution:
xor-reconstruction with x_j trivializing G, as in Case 8.1.3 (4 gates
removed). *Main branch:* if the second substitution has ΔΦ ≤ 2, it gives
Δμ ≥ 4 − 2α_φ + α_I, and the two-substitution average is
(1/2)(5 − 2α_φ + 2α_I + α_Q) per substitution, attaining
(5 − 2α_φ + α_Q)/2.

For the rest of 8.2.4.1, assume ΔΦ ≥ 3 in the second substitution.
Eliminating G and B (the node preceding P on the x_j → G path) is assumed
ΔΦ-free — failures are discharged in Case 8.2.5. The descendant of G
removed by Rule 3 contributes ≤ 1. So both inputs of P must each
contribute +1, giving four sub-configurations: (a) P fed by a gate T —
impossible: T must be ⊕-type (Case 5) and cannot become troubled by
degeneracy, since both its feeding variables stay influential with
undiminished out-degree; (b) P fed by two protected variables — impossible:
by Case 1 they are 1-variables and neither feed nor become troubled gates;
(c) variables t_1, t_2 are passed to P by the degeneracy of B and of the
eliminated descendant of x_k — Case 8.2.4.1.2; (d) t_1 directly feeds P,
and t_2 is passed via B and via P′ (the descendant of x_k), making t_2 a
3-variable feeding an ∧-type gate T that is troubled after P's elimination
— Cases 8.2.4.1.1 and 8.2.4.1.3.

*8.2.4.1.1: t_1 (or t_2) is protected.* Add a third substitution t_1 ← c
(or t_2 ← c), eliminating at least three gates via Rules 2 and 3 while
killing a quadratic equation; this step alone gives
Δμ ≥ 3 − 3α_φ + α_Q + α_I. Together the three substitutions give
Δμ ≥ (1/3)(8 − 8α_φ + 3α_I + 2α_Q) per substitution (value ≈ 12.93 at the
optimum, above δ).

*8.2.4.1.2 (vacuous): t_1, t_2 both passed to P during degeneracy* (t_2
only via B, t_1 only via P′). Then t_2 feeds Q and B feeds P before
x_k ← d; B must be a gate (B = x_j is impossible, x_j being a protected
2-variable before the substitution); Q = B ⊕ t_2 ⊕ c_1 forces B to depend
on t_2; then P = x_k ⊕ t_1 ⊕ B ⊕ c_2 depends on the unprotected t_2,
contradicting the standing condition of Case 8.2. No bound is needed.

*8.2.4.1.3 (vacuous): t_1 directly feeds P before x_k ← d.* If t_2 also
feeds P directly or only via P′, then P depends on the unprotected t_1,
impossible by Case 8.1; if t_2 is passed by B, then as in 8.2.4.1.2 either
P or Q depends on the unprotected t_1, again impossible by Case 8.1.
*Standing consequence:* P depends only on variables whose couples feed Q,
and vice versa.

**Case 8.2.4.2: P fed by two non-coupled protected variables.** Let P be
fed by x_j, x_k, protected variables not coupled to each other; both are
1-variables (Case 1) feeding only P. Substitute appropriate constants to
both, trivializing P and G; Q and at least one descendant of G are also
eliminated: 4 gates (P, G, Q, a descendant of G), with ΔΦ ≤ 3 (0 for P and
G; 2 for Q, one per input; 1 for the descendant). Two quadratic equations
are killed (Δq = −2) and Δi = −2, so over t = 2 substitutions
Δμ ≥ (1/2)(4 − 3α_φ + 2α_I + 2α_Q) per substitution (value 13.1, above δ).

**Case 8.2.4.3 (vacuous).** By 8.2.4.2, P and Q are each fed by at most
one variable. Suppose P depends on a protected x_j whose couple x_k feeds
Q. Then Q's other input is a gate R (by 8.2.4.2), non-constant (Case 0.3).
There is no path from x_k to R (Q is fed by x_k and R), so R does not
depend on x_k; but Q depends on x_k and on some variable x_i that R
depends on, and x_i ≠ x_j (else Q depends on the coupled pair x_j, x_k,
against the assumption). Hence Q depends on two non-coupled protected
variables, so P is fed by two protected variables, contradicting the
assumption of 8.2.4.2 under which this case arises — impossible.

### Case 8.2.5: discharging the ΔΦ-freeness assumption of 8.2.2.1 and 8.2.4.1

Cases 8.2.2.1 and 8.2.4.1 assumed that eliminating G and B cannot increase
the potential; we now discharge this. Let x_j ∈ I with couple x_k; a path
x_j → Q exists (Proposition 2.7). Substitute x_k ← d and normalize: x_j
becomes unprotected, and by Case 6 only the descendant of x_k is
eliminated, with no new troubled gate. If x_k feeds Q, its other input Q′
is passed to feed G (Q′ is ⊕-type, identified with Q). Xor-reconstruction
with x_j trivializes G. Let B be the node preceding Q on the x_j → Q path.
Eliminate G first, then B, then a descendant of G; no other elimination in
this process can increase the potential. The cases split on whether
eliminating G (Case 8.2.5.1) or eliminating B (Case 8.2.5.2 and its
reductions) gives ΔΦ = +1.

**Case 8.2.5.1: eliminating G increases the potential by 1.** Then P must
be a variable x_m whose out-degree dropped from 3 to 2 through G's
elimination; after x_k ← d, x_m is a 3-variable.

*8.2.5.1.1 (vacuous): x_m was already a 3-variable before x_k ← d.* Then
x_m cannot feed G before the substitution (Case 3), so it was passed to G
via the degeneracy of P′ (the descendant of x_k) when x_k ← d. For the
potential to increase, x_m must feed an ∧-type gate T right before G's
elimination, hence right after x_k ← d. Feeding T before the substitution
is impossible (Case 3); feeding T via P′'s degeneracy would require P′ to
be a 2+-gate, but T becoming troubled needs x_m to be a 2-variable after
G's elimination, forcing P′ to be only a 1-gate — contradiction.

*8.2.5.1.2: x_m was a 1-variable before x_k ← d — realizes α_I/3.* Then
x_m and x_k feed the same ⊕-type gate P′. Work in the circuit *before* the
x_k ← d substitution (after it, x_k is a constant and x_m a 3-variable, so
the analogy is only available there): P′ is fed by protected x_k (whose
couple is x_j) and the 1-variable x_m — exactly the configuration of
Case 6.2.2 with (x, u) := (x_k, x_m). All three subcases of 6.2.2 apply:
x_m unprotected — the 6.2.2.1 maneuver (constant to the couple x_j, then
the affine substitution x_k ← x_m), t = 2, Δi = −3; x_m = x_j, i.e.
coupled to x_k — the 6.2.2.2 affine substitution x_k ← x_m, t = 1,
Δi = −2; x_m protected with a couple different from x_k — the 6.2.2.3
maneuver (constants to both couples, then x_k ← x_m), t = 3, Δi = −4 —
the read-once property of rdq-sources makes the four killed inputs
distinct. The worst case is the last, so Δμ ≥ (4/3)α_I per substitution,
attaining the α_I/3 term.

*8.2.5.1.3: x_m was a 2-variable originally, raised to 3 by x_k ← d.*
Then P′ is an ⊕-type 2-gate fed by the protected x_k and x_m; x_m either
feeds G directly or is passed to G through P′'s degeneracy.

— *8.2.5.1.3.1: x_m feeds G directly, and P′ feeds an ∧-type gate T before
x_k ← d.* Substitute constants to x_m and x_k, trivializing P′, T and G;
their descendants are degenerated by normalization; t = 2. Six gates are
eliminated via Rules 2 and 3, with ΔΦ ≤ 6 (at most one per removed gate);
one quadratic equation is killed (Δq = −1, x_k's); Δi = −2. Hence
Δμ ≥ (6 − 6α_φ + α_Q + 2α_I)/2 per substitution — the average of a
Case-3-shaped and a Case-1-shaped decrement (value 12.9, above δ).
Henceforth x_m is passed to G through P′'s degeneracy.

— *8.2.5.1.3.2 (reduction): Q is a 2-gate and P′ is a 2+-gate depending on
the unprotected x_m.* Apply Case 8.1 with P′ and Q exchanged; the bound is
whichever Case 8.1 provides for that configuration, in particular ≥ δ.

— *8.2.5.1.3.3: Q is a 3+-gate.* First substitution x_k ← d eliminates P′
(Q becomes an ⊕-type gate depending on the unprotected x_j): a Case-6
step, Δμ ≥ 1 + α_Q + α_I. Second substitution: xor-reconstruction removing
five gates by Rules 2 and 3 (three descendants of Q, the first gate on the
reversed x_j → Q path, a descendant of G), with ΔΦ ≤ 5, so
Δμ ≥ 5 − 5α_φ + α_I; t = 2 in total. The average
(1/2)(6 − 5α_φ + 2α_I + α_Q) ≥ (1/2)(6 − 6α_φ + 2α_I + α_Q) per
substitution (value ≥ 12.9) is above δ.

**Case 8.2.5.2: eliminating B increases the potential by 1.** B feeds Q
already before x_k ← d. After the reconstruction, B's other input besides
the constant b is a node D that becomes (or feeds) a troubled gate; D feeds
Q before the reconstruction, originally either directly, or via P′ (the
descendant of x_k), passed to B during normalization after x_k ← d. D cannot be an ∧-type gate (minimality of G), so D is a
variable x_ℓ.

*8.2.5.2.1: x_ℓ is unprotected.* Then Q does not depend on x_ℓ (Case 8.1).
Since Q is fed by B and x_ℓ (or by B, x_ℓ, x_k), Q = B ⊕ x_ℓ ⊕ c (or
Q = B ⊕ x_ℓ ⊕ x_k ⊕ c), so B depends on x_ℓ.

— *8.2.5.2.1.1 (vacuous): B degenerates to a 0-gate right before its
elimination.* A 0-gate's out-degree must have decreased, so B must be the
gate Z added by the xor-reconstruction, i.e. B = x_j before the
reconstruction. Then Q is fed by x_j and x_ℓ and depends on both, but
B = x_j cannot depend on x_ℓ, while the setup says B depends on x_ℓ —
contradiction. Henceforth B is a 1+-gate before elimination; x_ℓ is passed
to the new troubled gate T during B's elimination, and the other input of
T is a variable t.

— *8.2.5.2.1.2 (vacuous): B and t both feed T before x_k ← d.* Then T is a
minimal ∧-type gate fed by the ⊕-type gate B, and B depends on the
unprotected x_ℓ — impossible by Case 8.1.

— *8.2.5.2.1.3 (vacuous): B is passed to T when P′ degenerates, and t
feeds T originally.* Since x_k is protected, P′ is ⊕-type; T is a minimal
∧-type gate fed by P′, and P′ depends on the unprotected x_ℓ — impossible
by Case 8.1. Henceforth t is passed to T when P′ degenerates and B feeds T
originally.

— *8.2.5.2.1.4: t is a 1-variable — realizes α_I/3.* Then P′ (the
descendant of x_k) is an ⊕-type gate fed by two 1-variables: protected x_k
(whose couple is x_j) and t — the configuration of Case 6.2.2 with
(x, u) := (x_k, t), the analogy driven by the protected feeder x_k
(without it the 6.2.2 maneuvers, which substitute couples, would not
apply). The three subcases: t unprotected — the 6.2.2.1 maneuver (t = 2,
Δi = −3); t = x_j, coupled to x_k — the 6.2.2.2 affine substitution
x_k ← t (t = 1, Δi = −2); t protected with another couple — the 6.2.2.3
maneuver (t = 3, Δi = −4, the four killed inputs distinct by the read-once
property). The worst case gives Δμ ≥ (4/3)α_I per substitution, attaining
the α_I/3 term.

— *8.2.5.2.1.5 (vacuous): t is a 2-variable.* For T to become troubled,
P′ must then be a 1-gate; but B depends on the unprotected x_ℓ, which is
impossible after Case 8.1.

*8.2.5.2.2: x_ℓ is protected.* Then x_ℓ is a protected 2-variable feeding
an ∧-type gate (it either feeds a truly troubled gate or is the inner
variable of a pack). Use three substitutions. (1) x_k ← d: a Case-6 step,
Δμ ≥ 1 + α_I + α_Q. (2) Xor-reconstruction removing B, G and a descendant
of G: G's elimination is ΔΦ-free — P ≠ B (else B = P would let G be
eliminated by input substitution rather than reconstruction), and P cannot
make G troubled (if P is a variable it is an unprotected 1-variable or a
2+-variable; if P is a gate, a path P → Q would close a cycle
G → P → Q → G through an ∧-type gate once P and Q are degenerate) — so
ΔΦ ≤ 2 and Δμ ≥ 3 − 2α_φ + α_I. (3) x_ℓ ← c: three gates are removed by
Rules 2 and 3; the originally troubled gate fed by x_ℓ is eliminated first
with both inputs 2-variables (no increment), the other two removals add
≤ 2, and destroying the troubled gate (or breaking the pack) subtracts 1,
netting ΔΦ ≤ 1; Δq = −1, so Δμ ≥ 3 − α_φ + α_I + α_Q. In total,
Δμ ≥ (7 − 3α_φ + 3α_I + 2α_Q)/3 per substitution over t = 3 (value ≈ 12.93
at the optimum, above δ; the displayed blend
(1/3)(7 − 4α_φ + 3α_I + 2α_Q) ≈ 12.87 is weaker and also suffices).

**Case 8.2.5.3 (reduction).** B feeds B′ and B′ feeds Q originally, and B′
is the gate eliminated by x_k ← d. Case 8.2.5.2's analysis still applies:
B is fed by a variable x_ℓ right after the reconstruction and x_ℓ feeds Q
originally; if x_ℓ is unprotected, Q does not depend on x_ℓ, and
Q = B ⊕ x_k ⊕ x_ℓ ⊕ c forces B to depend on x_ℓ, with B not a 0-gate
before elimination (otherwise B = x_j); if x_ℓ is protected, substitute
x_ℓ ← c after the whole process, mirroring step (3) of 8.2.5.2.2. The
bounds of 8.2.5.2 are inherited.

**Case 8.2.5.4 (reduction).** B feeds Q′ and Q′ feeds Q; Q is eliminated
by the substitution x_k ← d and Q′ takes the place of Q. Case 8.2.5.2's
analysis applies verbatim with Q′ in the role of Q, and its bounds are
inherited.

This exhausts the configurations: every step of the elimination procedure
falls under one of Cases 0–8, and in each leaf the per-substitution
decrement of μ is at least

    δ = α_I + min{ α_I/3, 2 − 2α_φ + α_Q, 4 − 4α_φ, 3 + α_φ, 5 − α_Q,
                   (5 − 2α_φ + α_Q)/2 },

with the six quantities attained respectively by the α_I/3 leaves
(5.4.1.4.1.4, 5.4.1.4.2.5, 6.2.2.3, 8.2.5.1.2, 8.2.5.2.1.4), Cases 1 and 6.1,
the 4 − 4α_φ leaves (3, 5.4.2.1, 7, 8.1.1, 8.1.3), the base leaf of
Case 5, Case 5.4.1.4, and Cases 8.2.2.1 and 8.2.4.1; pair elimination
(Case 0.4) and normalization are never worse. This proves Theorem 6.1. ∎
