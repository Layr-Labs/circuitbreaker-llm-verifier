# Case analysis for Theorem 6.1 — Part A: Cases 0–5

Notation and standing setup are as in `proof.md` §§3–5. We are given C
computing f | R, f an affine disperser for dimension d, the dimension of R
above 2d + 2, and a packing P. In each step we perform t ∈ {1, 2, 3}
substitutions (t = 0 for the preliminary simplifications of Case 0) and then
normalize; Δμ, ΔΦ, Δi, Δq denote the changes caused (Δμ the *decrease* of
μ). The per-substitution decrement must be ≥ δ of Theorem 6.1. Cases are
ordered: each case assumes the configurations of all previous cases are
absent. Lemma 5.3 (normalization) is used throughout without mention: after
any substitution we normalize freely since β ≥ 1 − 2α_φ ≥ 0 for Rule 1 and
β ≥ 1 − α_φ for the other rules, and any potential increase during
normalization comes only from gates fed by the eliminated gate's inputs.

## Case 0: preliminary simplifications (only Δμ ≥ 0 is needed)

**0.1** If C is not normalized, normalize it. By Lemma 5.3, Δμ ≥ 0.

**0.2** Suppose a gate in the acyclic part of C does not reach the output
gate; take a topologically maximal such gate G. Then G is a 0-gate, and
removing it by Rule 1 costs ΔΦ ≤ 2, hence Δμ ≥ 1 − 2α_φ ≥ 0. Henceforth
every gate reaches the output.

**0.3** If a gate G outputs a fixed constant c under every assignment,
replace G by c: one gate is removed and the potential rises by at most 2
(one for each input of G), so Δμ ≥ 1 − 2α_φ ≥ 0. If a gate G fed by I_1, I_2
computes a function depending only on I_1, rewire it to (I_1, 0) — this
raises the potential by at most 1 — after which G is degenerate and is
removed by Rule 3, again with total ΔΦ ≤ 2 and Δμ ≥ 1 − 2α_φ ≥ 0. Henceforth
every gate's output depends on both its inputs.

**0.4** If the packing P is non-empty, apply troubled pair elimination
(Lemma 5.4): the decrement per substitution is ≥ δ_p = min{ (3/2)α_I,
3 − 3α_φ + α_I + α_Q } = 69/5 > δ at our parameters. Henceforth P = ∅ can be
assumed in the remaining cases: the packing is empty (troubled pairs may
still exist in C, unpacked), so every troubled gate is truly troubled with
respect to P, and Φ counts troubled gates directly.

## Case 1: a protected variable feeding two gates, or an ∧-type gate

Substitute an appropriate constant to this protected variable; t = 1. The
substitution kills one quadratic equation (Δq = −1) and one influential
input (Δi = −1), and at least two gates become degenerate and are removed
by Rules 2 and 3, while ΔΦ ≤ 2 (only gates fed by the substituted variable
can become troubled, at most one per eliminated gate along the two
affected descendants). Hence

    Δμ ≥ 2 − 2α_φ + α_I + α_Q,

realizing the term 2 − 2α_φ + α_Q (value 17/5 at the optimum; slack).
*Standing consequence:* every protected variable now feeds exactly one
gate, and that gate is ⊕-type.

## Case 2: a protected 0-variable

Let p be a protected variable feeding no gate, and q its couple in their
quadratic equation. Substitute a constant to q; t = 1. The quadratic
equation becomes affine, and both p and q become non-influential
(Δi = −2). Hence Δμ ≥ 2α_I, far above δ. *Standing consequence:* every
protected variable is a 1-variable feeding an ⊕-type gate.

## Case 3: a variable feeding an ∧-type gate, with total out-degree ≥ 4 — realizes 4 − 4α_φ

Suppose a variable x feeds an ∧-type gate G and the total out-degree of x
and G is at least 4. Substitute an appropriate constant to x, trivializing
G; t = 1. At least three descendants of x and G become degenerate; a gate
fed by both G and x (counted twice in the out-degree sum) is trivialized
and its descendants become degenerate; altogether at least four gates are
eliminated by Rules 2 and 3. The potential rises by at most one per
eliminated gate, ΔΦ ≤ 4, and Δi = −1 (x). Hence

    Δμ ≥ 4 − 4α_φ + α_I,

realizing the term 4 − 4α_φ. *Standing consequence:* every variable feeding
an ∧-type gate has out-degree 1 or 2, and every ∧-type gate fed by a
variable is a 1-gate or 2-gate subject to the later cases.

## Case 4: an ∧-type gate fed by two variables, one a 1-variable

Let G be ∧-type, fed by variables x and y with x a 1-variable. By Case 1,
x and y are unprotected. Substitute an appropriate constant to y,
trivializing G; t = 1. Both x and y become non-influential (Δi = −2), and
no potential increase is possible (the only gates whose inputs change are
fed by x or y; x's single descendant is G itself). Hence Δμ ≥ 2α_I.
*Standing consequence:* an ∧-type gate fed by two variables is fed by two
2-variables — i.e. it is a troubled gate (Case 5) — or by variables one of
which is handled in Case 8.

## Case 5: an ∧-type gate fed by two 2-variables (the troubled configuration)

Let G be ∧-type, fed by two 2-variables x and y. By Case 3, G is a 1-gate,
hence troubled. Let B (resp. C) be the other descendant of x (resp. y) and
D the descendant of G; B ≠ D and C ≠ D since C is normalized, and x, y are
unprotected free variables by Case 1. Substituting a constant to x
trivializes G and removes B, D, G by normalization. If no troubled gate is
introduced while eliminating these three gates, we are in the base leaf:

**Case 5, base leaf — realizes 3 + α_φ.** Gates G, B, D are eliminated
(3 gates); the eliminated G was troubled, so ΔΦ = −1; Δi = −1 (x); Δq = 0.
Hence, exactly,

    Δμ = 3 + α_φ + α_I,

realizing the term 3 + α_φ. For Cases 5.1–5.4 we therefore assume that at
least one troubled gate IS introduced by the elimination above.

**Case 5.1: B = C.** By truth-table enumeration over the finitely many
∧-type/⊕-type choices for G and B (the enumeration is finite: 16 two-input
functions for each gate, and only the induced behavior on the four
assignments to (x, y) matters), there is a constant substitution to x or y,
or an affine substitution x ← y ⊕ c, making both G and B output fixed
constants; t = 1. Both x and y become non-influential (Δi = −2), so
Δμ ≥ 2α_I.

**Case 5.2: D feeds both B and C.** Run the base plan (constant to x;
eliminate G, D, B in this order). By the normalization lemma
(Lemma 5.3), a troubled gate can be introduced only at a gate fed by the
other input I of D after the elimination; if I is a gate it cannot
become troubled (its out-degree does not decrease) and no troubled gate
is introduced, so we are in the base leaf. Assume then that a troubled
gate is introduced; then I is a variable u. The gates whose inputs
change are exactly the descendants of D and of B; C receives the
nonconstant affine function of u that D passes when it degenerates,
and each descendant of B receives the signal B passes, likewise a
nonconstant affine function of u (if B instead
passes a constant, its descendants are eliminated on constants already
in the first step, which only adds gates there). Hence *every* newly
introduced troubled gate is fed by a u-dependent signal after the base
plan. Let A be one of them. (Remark: the source paper's own operation
for this case — constants to the other input of the descendant below B
and to y — covers only the configuration where u ends up feeding
exactly two gates and the troubled gate sits below B; the operation
below is uniform and covers C itself being troubled, further
descendants of D, and B being fed by x and D in any position.)

First, B is not the output gate in this branch: by Case 0 the output
is reachable from C, but a path from C to B would have to enter B
through one of its inputs x or D, and D through one of its inputs G
or u — every such entry point is a variable, which no gate feeds; so
no path exists. Hence B has a descendant, B is eliminated in the base
plan by Rule 3 without the output exception, and after the base plan
the u-dependent signals feed C and — unless B passed a constant, in
which case B's descendants were already eliminated in the first step,
adding gates there — at least one descendant of B (plus possibly
further descendants of D).

Append a second substitution: a constant to u, at its controlling
value for the ∧-type gate A (Case 3's move on the pair (u, A); if u is
protected, Δq = −1, which only helps); t = 2 in total. We ledger the
combined operation end to end. Gates eliminated in the second step:
A itself, trivialized; every other gate fed by a u-dependent signal
(at least one), which degenerates on its now-constant input and is
eliminated (by Rule 3, or under (T2) if it is the output gate); and
one descendant of A, which degenerates on receiving A's constant
output — A is not the output gate, by (T1), and being useful it has a
descendant. If that descendant is itself one of the u-fed gates just
counted, take its descendant instead: it is not the output gate
either (both of its inputs become constant in this step, against
(T1)), so it has a descendant, which degenerates and is eliminated.
The two steps therefore eliminate at least 6 gates in total (at least
3 of them in the second step, or at least 4 in the first step in the
constant-passing case of the first paragraph). Influential inputs killed: Δi = −3 — x (its descendants G and
B are eliminated), u (substituted), and y (its descendants are G and
C; C is fed by a u-dependent signal after the base plan, so C is
trivialized or degenerated in the second step). Potential: G is destroyed; every
introduced troubled gate is fed by a u-dependent signal and is
destroyed in the second step; a gate of the final circuit can be newly
troubled only if it is fed by a variable passed during the second
step (Lemma 5.3), and the gates introduced through any one passed
variable pairwise share that variable, hence contribute at most one
pack among them (the packing observation of Rule 1 in `proof.md` §5).
In the worst case — the u-dependent signals feed exactly two gates —
the passed variables are two (the other input of the one degenerated
u-fed gate, and the other input of A's descendant), so
ΔΦ ≤ −1 + 2 = +1; each further gate fed by a u-dependent signal adds
one eliminated gate against at most one further passed variable, and
the constant-passing case has only one passed variable, so the
estimate only improves. Hence

    Δμ ≥ 6 + 3α_I − α_φ = 6 + 144/5 − 1/5 = 173/5 ≥ 2δ

over t = 2 substitutions. (If the base plan introduces a second
troubled gate as well, it is also fed by a u-dependent signal, so it
is trivialized or degenerated in the second step; the count is
unchanged and the bound stands.)

**Case 5.3: C feeds D and D feeds B (vacuous).** Run the base plan
(constant to x; eliminate G, D, B in this order). The only nodes that could
become troubled, or feed a gate that becomes troubled, are y and the
introduced constant c; neither can become or create a troubled gate. Hence
no troubled gate is introduced, contradicting the standing assumption of
Case 5; the configuration falls back to the base leaf,
Δμ = 3 + α_φ + α_I.

**Case 5.4: the remaining topology.** By 5.3, D is disconnected (in both
directions) from one of B, C; by symmetry from B. Write B's other input as
I_1 and D's other input as I_2 (D is fed by G and I_2). Substitute a
constant to x to trivialize G. A newly introduced troubled gate is either
I_1 or I_2, or is fed by I_1 or I_2 after normalization: I_1, I_2 themselves
cannot become troubled, since their out-degrees do not decrease unless
their own inputs are involved. So some new troubled gate is fed by I_1 or
I_2 after normalization. Case 5.4.1 treats a troubled gate fed by I_2;
Case 5.4.2 (in the continuation of this file) treats one fed by I_1.

### Case 5.4.1: a new troubled gate E fed by I_2

Let E be the new troubled gate, fed by I_2. First, I_2 must be a variable
z: if I_2 were a gate it could not feed a *newly* troubled gate in this
position. E is an ∧-type 1-gate, and since B and D are disconnected, D
itself must feed E before normalization.

**Case 5.4.1.1: z is a 2-variable.** Then D must be a 1-gate for E to be
troubled. The other input of E is either a variable t, or B — where B,
after x is substituted, degenerates and passes its other input t to E (B is
fed by x and t; x is unprotected by Case 1).

*5.4.1.1.1: the other input of E is a variable t.* Substitute a constant
to t (not to x), trivializing E; eliminate E, D, G and the descendant of E
in this order; t = 1. Four gates are eliminated. The potential rises by at
most 1 (only the other input of the descendant of E can introduce a
troubled gate), and Δi = −1 (t). Hence Δμ ≥ 4 − α_φ + α_I, stronger than
the 4 − 4α_φ shape.

*5.4.1.1.2: the other input of E is B, and B is ∧-type.* Substitute
constants to t and y, trivializing B and G; x becomes a 0-variable;
t = 2. Three influential inputs die (Δi = −3), so Δμ ≥ (3/2)α_I per
substitution.

*5.4.1.1.3: the other input of E is B, and B is ⊕-type.* Apply the affine
substitution x ← t ⊕ c, fixing B's output to a constant; replace B by that
constant. With c chosen appropriately E is also trivialized; D and G are
removed by Rule 1 and the descendant of E degenerates by Rule 3; t = 1.
Five gates are eliminated (B, E, D, G, descendant of E). Only t and the
other input of the descendant of E can introduce troubled gates (y and z
are 2-variables and cannot), so ΔΦ ≤ 2, and Δi = −1 (x). Hence
Δμ ≥ 5 − 2α_φ + α_I.

For the rest of Case 5.4.1, z is a 1-variable; then D must be a 2-gate for
E to be troubled.

**Case 5.4.1.2: D is an ∧-type 2-gate.** Substitute a constant to z,
trivializing D; eliminate D, G and the two descendants of D in order;
t = 1. Four gates are eliminated; eliminating G causes no potential
increment, and overall ΔΦ ≤ 3; Δi = −1 (z). Hence Δμ ≥ 4 − 3α_φ + α_I.

**Case 5.4.1.3: z is protected.** Substitute a constant to x (G is
trivialized; B and D become degenerate), then a constant to z (E is
trivialized; the descendant of E and the copy of z passed through D become
degenerate); t = 2. Six gates are removed in total (G, B, D, E, the
descendant of E, and one further gate degenerated by the two substitutions;
the paper states the count of 6 removed gates without itemizing them, and
it is this count that the accounting below uses), all via Rules 2 and 3,
each removal introducing at most one troubled gate, so ΔΦ ≤ 6. The quadratic equation involving the
protected z is removed (Δq = −1), and Δi = −2 (x, z). The paper's
accounting decomposes the two substitutions as a Case-1-shaped step and a
Case-3-shaped step:

    Δμ ≥ [ (2 − 2α_φ + α_I + α_Q) + (4 − 4α_φ + α_I) ] / 2 per substitution,

a blend of the slack term and the 4 − 4α_φ term (value 12.9 at the
optimum, above δ).

**Case 5.4.1.4: z unprotected, D an ⊕-type 2-gate — main line realizes 5 − α_Q.**
Assume first that eliminating the descendant of D (step (4) below) does
not introduce a troubled gate; the contrary is Cases 5.4.1.4.1/5.4.1.4.2
below. Apply ONE quadratic substitution

    z ← ((x ⊕ c_1) ∧ (y ⊕ c_2)) ⊕ c_3,

with constants chosen so that D's output becomes a constant b insensitive
to its inputs, which also trivializes E; t = 1. This adds one quadratic
equation (Δq = +1). Normalize in this order: (1) replace D by b
(ΔΦ = 0); (2) eliminate G by Rule 1 (ΔΦ = 0; x and y become 1-variables);
(3) eliminate E by Rule 2 (ΔΦ = 0: E's other input is B, fed by the
1-variable x and hence not troubled, or a variable t that becomes a
1-variable); (4) eliminate the descendant of D by Rule 3 (ΔΦ = 0, by the
standing assumption); (5) eliminate the descendant of E by Rule 3
(ΔΦ ≤ 1). Five gates are eliminated (D, G, E, descendant of D, descendant
of E). The originally troubled G is destroyed, contributing −1 to the
potential and cancelling step (5): net ΔΦ ≤ 0. Only z becomes
non-influential (Δi = −1). Hence

    Δμ ≥ 5 − α_Q + α_I,

realizing the term 5 − α_Q.

### The contrary of 5.4.1.4: eliminating the descendant of D introduces a troubled gate

Suppose step (4) above does introduce a troubled gate A. Then F, the
descendant of D eliminated in step (4), is directly fed by some variable u,
and A is directly fed by F and a 2+-variable v (once F is eliminated in
step (4), A is fed by the signal F passes — a function of u — and v; this
is the form used in 5.4.1.4.2). By the initial assumption of
Case 5, a variable t either directly feeds E, or is passed to E while B
degenerates. Case 5.4.1.4.1 is the latter (t feeds B, B feeds E);
Case 5.4.1.4.2 is the former (t feeds E directly).

**5.4.1.4.1: t feeds B and B feeds E.**

*5.4.1.4.1.1: B is a 1-gate.* Substitute constants to y and z so that G, D,
E are all trivialized; B's out-degree drops to 0 and B is removed by
Rule 1; t = 2. Three influential inputs die (Δi = −3: x becomes a
0-variable; y, z substituted), so Δμ ≥ (3/2)α_I per substitution.
Henceforth B is a 2-gate and t is a 1-variable.

*5.4.1.4.1.2: B is an ∧-type 2-gate.* Then t, a 1-variable, is unprotected
(by Case 1). Substitute a constant to x, trivializing B and making t
non-influential; t = 1, Δi = −2, so Δμ ≥ 2α_I.

*5.4.1.4.1.3: B is an ⊕-type 2-gate, t unprotected.* Substitute a constant
to y, trivializing G (x becomes a 1-variable), then the affine substitution
x ← t, trivializing B; t = 2, Δi = −3 (x, y, t), so Δμ ≥ (3/2)α_I per
substitution.

*5.4.1.4.1.4: B is an ⊕-type 2-gate, t protected — realizes α_I/3.*
The couple of t is neither x nor y (both unprotected). Substitute a
constant to the couple of t, making t unprotected, then apply 5.4.1.4.1.3;
t = 3, and four influential inputs die (Δi = −4: the couple of t, y, x, t).
Hence Δμ ≥ (4/3)α_I per substitution, i.e. per-substitution
Δμ − α_I = α_I/3. This leaf attains the α_I/3 term.

**5.4.1.4.2: t feeds E directly.** After the quadratic substitution of
5.4.1.4, A is fed by u and v; x and y are 1-variables; u, v, x, y are
pairwise distinct (so F ≠ B, C) and u, v ≠ z.

*5.4.1.4.2.1: u = t.* Then u is a 2-variable feeding E and F, hence
unprotected (Case 1). Substitute a constant to x, trivializing G; remove B
and D by Rule 3; then E and F are each fed by two 2-variables among z, u,
so substitute a constant to z or u, or apply the affine substitution
u ← z ⊕ c, making both E and F output input-independent functions; t = 2
in total, Δi = −3 (x, z, u), Δμ ≥ (3/2)α_I per substitution. Henceforth
u ≠ t.

*5.4.1.4.2.2: F is a 1-gate.* Substitute a constant to v, replacing the
troubled A by a constant; F becomes a 0-gate and is removed by Rule 1; a
descendant of A is removed by Rule 3; another descendant of v is removed
by Rule 3 (v has out-degree ≥ 2). Four gates are eliminated. The potential
rises by ≤ 2 in total: trivializing A contributes 0 (v is constant by then
and F is fed by the ⊕-type D), removing F contributes 0 (D is ⊕-type and
u becomes a 1-variable), and the two descendant eliminations contribute
≤ 1 each. With Δi = −1 (v): Δμ ≥ 4 − 2α_φ + α_I for this substitution.

*Pair accounting for 5.4.1.4.2.2 (t = 2).* Theorem 6.1 requires the two
substitutions together to drop μ by at least 2δ, so we account for the pair
jointly, listing each eliminated gate once: D, G, E, the descendant of E
(first substitution), and F, A, the descendant of A, the descendant of v
(second substitution) — 8 gates; Δq = +1 (the quadratic equation of
5.4.1.4). The potential events are: G, troubled, is destroyed (−1); A is
created by the first substitution's step (4) (+1) and destroyed by the
second (−1); the descendant-of-E elimination contributes ≤ +1; the two
Rule-3 removals of the second substitution contribute ≤ +1 each; the
removal of F contributes 0 (shown above). Net ΔΦ ≤ +2. For the influential
count, z dies at the first substitution and v at the second; moreover u
also dies: u feeds exactly F and A in this configuration — this is the
content of the paper's remark above that removing F makes u a 1-variable
(before that removal u feeds F and the still-standing A, and nothing else)
— so once A is replaced by a constant and F is removed, u is a 0-variable,
hence non-influential: Δi = −3 over the pair. The pair therefore drops μ by

    Δμ ≥ 8 − α_Q + 3α_I − 2α_φ   (over t = 2 substitutions),

which at the parameters of §7 is far above 2δ; in particular the pair
clears 2δ with slack to spare, and this leaf is not binding.

*5.4.1.4.2.3: F is an ∧-type 2-gate.* (F is forced to be exactly a 2-gate:
it is a 2+-gate originally, feeds none of G, D, E, and must be exactly 2
for A to be troubled; u is a 1-variable, u ≠ t.) Substitute constants to u
and t, trivializing E and F; D becomes a 0-gate and is removed by Rule 1,
making the unprotected z a 0-gate, hence non-influential; t = 2,
Δi = −3 (u, t, z), Δμ ≥ (3/2)α_I per substitution.

*5.4.1.4.2.4: F is an ⊕-type 2-gate, u unprotected.* Substitute a constant
to t, trivializing E (D becomes a 1-gate). Since
F = D ⊕ u ⊕ c_1 = G ⊕ z ⊕ u ⊕ c_2, the affine substitution z ← u ⊕ c makes
F = G ⊕ c ⊕ c_2; rewire all descendants of F to be fed by G directly,
removing D and F; z and u become 0-variables; t = 2, Δi = −3 (t, z, u),
Δμ ≥ (3/2)α_I per substitution.

*5.4.1.4.2.5: F is an ⊕-type 2-gate, u protected — realizes α_I/3.*
The couple of u is none of x, y, z, t, v (none of them is protected).
Substitute a constant to the couple of u, making u unprotected, then apply
5.4.1.4.2.4; t = 3, Δi = −4 (the couple of u, t, z, u). Hence
Δμ ≥ (4/3)α_I per substitution, attaining the α_I/3 term.

### Case 5.4.2: the new troubled gate E is fed by I_1

Now E does not receive a variable from D; instead E is fed by I_1 = t, the
input passed from B. The other input of E must be a variable z.

**Case 5.4.2.1: B is a 2+-gate — realizes 4 − 4α_φ.** Then B is ⊕-type
(by Case 3). Apply the affine substitution x ← t ⊕ c (x is unprotected, by
Case 1), fixing B's output regardless of its inputs; replace B by a
constant, choosing c so that E is also trivialized; rewire G to be fed by
t instead of x, disconnecting x — this introduces no troubled gates;
t = 1. Four gates are eliminated (B replaced by a constant; E and two
descendants removed by Rules 2 and 3). The total potential increase is
≤ 4: the other three removals contribute ≤ 1 each by the normalization
lemma accounting, and replacing B by a constant contributes at most 1, by
the following inspection of the local configuration. The only gates whose
inputs change in that step are the descendants of B and G (rewired from x
to t). A descendant of B gains a constant as one input, so it cannot become
troubled (a troubled gate is fed by two variables); it can at most
degenerate and pass its other input through. G was already the troubled
gate of the standing configuration and remains fed by two variables (t
replacing x), so it is not a *new* troubled gate. Hence a new troubled
gate can arise only downstream of a degeneration of a descendant of B,
mediated by the single variable t that those degenerations can pass; at
most one such gate (or one pack) results. With Δi = −1 (x):

    Δμ ≥ 4 − 4α_φ + α_I,

attaining the 4 − 4α_φ term.

**Case 5.4.2.2: B is a 1-gate.** Substitute constants to z and y,
trivializing both E and G; B becomes a 0-gate and is removed by Rule 1;
t = 2, Δi = −3 (z, y, and x, which loses both its descendants G and B),
so Δμ ≥ (3/2)α_I per substitution.

*Standing result after Case 5:* no ∧-type gate is directly fed by two
variables. From here on, let G be a topologically minimal ∧-type gate, fed
by I_1 and I_2, each computed by a fair cyclic xor-circuit; the
xor-reconstruction operation of `proof.md` §3 applies to I_1 and I_2.
