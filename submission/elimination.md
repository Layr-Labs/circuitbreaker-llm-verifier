# Case analysis for Theorem 6.1: elimination

Notation and setup as in `proof.md` §§3–5. Given C computing f | R, f an
affine disperser for dimension d, dim R above 2d + 2, and a packing P. Each
step performs t ∈ {1, 2, 3} substitutions (t = 0 in Case 0), then
normalizes; Δμ (the *decrease* of μ), ΔΦ, Δi, Δq are the resulting changes,
and the per-substitution decrement must be ≥ δ. Cases are in priority
order: each assumes all previous configurations absent. Lemma 5.3
(normalization) is used throughout without mention: normalization is free,
since β ≥ 1 − 2α_φ ≥ 0 for Rule 1 and β ≥ 1 − α_φ for the other rules, and
any potential increase during normalization comes only from gates fed by
the eliminated gate's inputs. Parameters: α_I = 48/5, α_Q = 9/5,
α_φ = 1/5, δ = 64/5.

## Case 0: preliminary simplifications (only Δμ ≥ 0 is needed)

**0.1** C not normalized: normalize. By Lemma 5.3, Δμ ≥ 0.

**0.2** A gate in the acyclic part does not reach the output; take a
topologically maximal such G. G is a 0-gate; Rule 1 removal costs
ΔΦ ≤ 2, hence Δμ ≥ 1 − 2α_φ ≥ 0. Henceforth every gate reaches the output.

**0.3** A gate G outputs a fixed constant c under every assignment:
replace G by c; one gate removed, potential rises ≤ 2 (one per input),
Δμ ≥ 1 − 2α_φ ≥ 0. A gate G fed by I_1, I_2 depending only on I_1: rewire
to (I_1, 0) — potential +≤ 1 — then G is degenerate and Rule 3 removes it,
again total ΔΦ ≤ 2, Δμ ≥ 1 − 2α_φ ≥ 0. Henceforth every gate depends on
both inputs.

**0.4** Packing P non-empty: troubled pair elimination (Lemma 5.4) gives
per-substitution decrement ≥ δ_p = min{ (3/2)α_I, 3 − 3α_φ + α_I + α_Q }
= 69/5 > δ. Henceforth P = ∅ (troubled pairs may still exist, unpacked);
every troubled gate is truly troubled w.r.t. P and Φ counts troubled gates
directly.

## Case 1: a protected variable feeding two gates, or an ∧-type gate

Substitute an appropriate constant to this protected variable; t = 1. One
quadratic equation killed (Δq = −1), one influential input killed
(Δi = −1), ≥ 2 gates degenerate and are removed by Rules 2 and 3,
ΔΦ ≤ 2 (only gates fed by the substituted variable can become troubled,
≤ 1 per eliminated gate along the two affected descendants). Hence

    Δμ ≥ 2 − 2α_φ + α_I + α_Q,

realizing 2 − 2α_φ + α_Q (17/5 at the optimum; slack). *Standing
consequence:* every protected variable feeds exactly one gate, ⊕-type.

## Case 2: a protected 0-variable

p protected, feeding no gate; q its couple. Substitute a constant to q;
t = 1. The equation becomes affine; p and q both non-influential
(Δi = −2). Δμ ≥ 2α_I, far above δ. *Standing consequence:* every protected
variable is a 1-variable feeding an ⊕-type gate.

## Case 3: a variable feeding an ∧-type gate, total out-degree ≥ 4 — realizes 4 − 4α_φ

x feeds an ∧-type gate G; total out-degree of x and G ≥ 4. Substitute an
appropriate constant to x, trivializing G; t = 1. ≥ 3 descendants of x and
G become degenerate; a gate fed by both G and x (counted twice in the
out-degree sum) is trivialized and its descendants degenerate; ≥ 4 gates
eliminated in all by Rules 2 and 3. ΔΦ ≤ 4 (≤ 1 per eliminated gate),
Δi = −1 (x). Hence

    Δμ ≥ 4 − 4α_φ + α_I.

*Standing consequence:* every variable feeding an ∧-type gate has
out-degree 1 or 2; every such ∧-type gate is a 1-gate or 2-gate.

## Case 4: an ∧-type gate fed by two variables, one a 1-variable

G ∧-type, fed by variables x, y with x a 1-variable; x, y unprotected
(Case 1). Substitute an appropriate constant to y, trivializing G; t = 1.
Both x, y non-influential (Δi = −2); no potential increase (only gates fed
by x or y change inputs; x's single descendant is G). Δμ ≥ 2α_I.
*Standing consequence:* an ∧-type gate fed by two variables is fed by two
2-variables — a troubled gate (Case 5) — or is handled in Case 8.

## Case 5: an ∧-type gate fed by two 2-variables (the troubled configuration)

G ∧-type, fed by 2-variables x, y. By Case 3, G is a 1-gate, hence
troubled. B (resp. C) the other descendant of x (resp. y), D the descendant
of G; B ≠ D, C ≠ D (C normalized); x, y unprotected free variables
(Case 1). Base plan: constant to x trivializes G; B, D, G removed by
normalization.

**Base leaf — realizes 3 + α_φ.** No troubled gate introduced: G, B, D
eliminated (3 gates); G was troubled, so ΔΦ = −1; Δi = −1 (x); Δq = 0.
Exactly,

    Δμ = 3 + α_φ + α_I.

Cases 5.1–5.4 assume a troubled gate IS introduced by this elimination.

**5.1: B = C.** Truth-table enumeration over the ∧-type/⊕-type choices
for G and B (16 two-input functions each; only the induced behavior on the
four assignments to (x, y) matters) gives a constant substitution to x or
y, or an affine substitution x ← y ⊕ c, making both G and B output fixed
constants; t = 1. Δi = −2 (x, y), so Δμ ≥ 2α_I.

**5.2: D feeds both B and C.** Run the base plan (constant to x; eliminate
G, D, B in this order). By Lemma 5.3, a troubled gate can be introduced
only at a gate fed by the other input I of D after the elimination; if I
is a gate it cannot become troubled (its out-degree does not decrease) and
none is introduced — the base leaf. Assume one is introduced; then I is a
variable u. The gates whose inputs change are exactly the descendants of
D and of B. When D degenerates, it passes a nonconstant affine function
s(u) of u, and D and B are then removed by Rule 3, which rewires their
descendants to be fed by u itself; hence every newly introduced troubled
gate is a descendant of D or of B — a descendant of C keeps a gate input,
so it cannot become troubled — and, apart from the constant cascade below
a trivialized B treated in 5.2.γ, it is fed by u itself after the base
plan. Let A be one of them. (Remark: the source paper's operation for this case —
constants to the other input of the descendant below B and to y — covers
only the configuration where the troubled gate sits below B; the operation
below is uniform apart from one residual branch, 5.2.β4, which is exactly
that configuration and uses the source's operation there.)

First, B is not the output gate in this branch: by Case 0 the output is
reachable from C, but a path from C to B would enter B through x or D, and
D through G or u — every such entry point is a variable, which no gate
feeds, so no path exists. The same entry-point argument shows C is not the output gate either: B reaches the output (the circuit is normalized, so every gate does), but a path from B to C would have to enter C through one of its inputs, y or D — through the variable y, which no gate feeds, or through D, which is entered only through G or u, and G only through the variables x and y; and D feeds B, so such a path would in any case close a cycle. Hence C, like B, has a descendant. Hence B has a descendant and is eliminated in
the base plan by Rule 3 without the output exception; after the base plan
the u-dependent signals feed C and — unless B passed a constant, in which
case B's descendants were already eliminated in the first step, adding
gates there — ≥ 1 descendant of B (plus possibly further descendants
of D).

If the constant output of G is controlling for D, D is trivialized too,
and the base plan is already a leaf: B, fed by two constants (the
substituted x-value and D's output), is trivialized, so a descendant of B
is eliminated on constants in the same step, and C, fed by D's constant
output, is eliminated as well, degenerating (passing y) or trivialized.
At least five gates are eliminated (G, D, B, C, and a descendant of B),
Δi = −1 (x; it is −2 if u feeds no gate besides D, for then u is dead
too), and ΔΦ ≤ +1 in the minimal configuration (G is destroyed; at most
one pack forms through the passed y and one through the other input of
B's eliminated descendant); every further pack forms through a variable
passed by a gate eliminated in the same step beyond the minimal count,
and every such further step contributes β ≥ 1 − 2α_φ ≥ 0 to Δμ by
Lemma 5.3, so the minimal configuration is the worst case. Hence
Δμ ≥ 5 + α_I − α_φ = 72/5 ≥ δ at t = 1. Assume henceforth that D
degenerates, passing s(u).

*5.2.γ: B is trivialized in the base plan.* Then a descendant W of B is
eliminated on B's constant output in the first step — at least four gates
are eliminated there (G, D, B, W) — and a troubled gate can be introduced
in the constant cascade below B, fed not by u but by a variable passed
through the cascade. Two cases. If a troubled gate A fed by u is also
introduced (below D, or in the cascade through a passer itself fed by u),
substitute a constant to u as in the sub-cases below — controlling for C
if C is ∧-type, else controlling for A: the second step eliminates A and
C, and a descendant of whichever of them is trivialized by the chosen
constant (if the constant is controlling for C but not for A, A
degenerates instead and the third gate is a descendant of C; A's pins are
variables, so none of these gates coincides with another) — at least
seven gates in total, with Δi = −2 and ΔΦ ≤ +2 in the minimal
configuration, so Δμ ≥ 7 + 2α_I − 2α_φ = 129/5 ≥ 2δ. If A = C, the second
step eliminates only C and a descendant of C, but y dies as well
(Δi = −3), giving Δμ ≥ 6 + 3α_I − 2α_φ = 172/5. Each further cascade gate
created in the first step is accompanied by a further gate eliminated
there (β ≥ 1 − 2α_φ ≥ 0 per extra step, by Lemma 5.3), so the minimal
configuration is the worst case in both γ sub-cases. If the only troubled
gates introduced are in the cascade below B, let Z be one of them, fed by
the passed variable J and a second variable (J is a 2-variable, since Z
is troubled); substitute J₀, controlling for Z: Z is trivialized, and a descendant of Z and the other
descendant of J are eliminated on constants. If the descendant of Z
just counted is the other descendant of J itself, it receives constants
on both pins (Z's output and J₀), so it is trivialized, and its own
descendant is eliminated on a constant instead; the count is unchanged.
At least seven gates in total, Δi = −2 (x, J), ΔΦ ≤ +2 in the minimal
configuration, and Δμ ≥ 7 + 2α_I − 2α_φ = 129/5 ≥ 2δ. Assume henceforth that B degenerates
in the base plan, passing s(u); then every newly introduced troubled gate
is fed by u.

The second substitution is a constant to u (t = 2 in total), but its
effect on y depends on the type of C: the descendants of y are exactly G
and C, and G is already eliminated, so y is killed exactly when C is
eliminated without passing y on. Split on the type of C.

*Sub-case 5.2.β: C is an ⊕-type gate.* With its signal input fixed, an
⊕-type C passes y on, so no choice of the constant to u trivializes C; y
survives the second substitution in every branch below, and the ledger is
closed at Δi = −2 by gates and by the potential. In 5.2.β0–β3 the constant
to u is chosen controlling for a troubled gate A (A is ∧-type, and s takes
both values), so A is trivialized and a descendant of A is eliminated on a
constant in the same step; 5.2.β4 is the one configuration where that
count is not enough.

*5.2.β0: u is protected.* The lost third death is paid for by Δq = −1:
the second step eliminates at least C and A (five gates in total), with
ΔΦ ≤ +1, so Δμ ≥ 5 + 2α_I + α_Q − α_φ = 129/5 ≥ 2δ. Assume henceforth
that u is unprotected.

*5.2.β1: two or more troubled gates are introduced.* Say k ≥ 2 are
introduced. The second step eliminates all k of them, a descendant of the
trivialized one, and C: at least k + 5 gates in total. The k troubled
gates are created in the first step and destroyed in the second, G is
destroyed, and the passed variables of the second step contribute at most
k − 1 packs in total (the trivialized gate passes a constant; each of the
other k − 1 troubled gates and the descendant of the trivialized one
passes at most one variable; no pack forms through y, a 1-variable after
the base plan), so ΔΦ ≤ k − 1. Hence
Δμ ≥ (k + 5) + 2α_I − (k − 1)α_φ, which increases with k (each further
troubled gate adds a gate and at most one pack, a net ≥ 1 − α_φ > 0), and
the worst case is k = 2: Δμ ≥ 7 + 2α_I − α_φ = 130/5 ≥ 2δ. Assume
henceforth that exactly one troubled gate A is introduced.

*5.2.β2: A is a descendant of D.* The second step eliminates A, a
descendant of A, C, and a descendant of B (the last is fed by s(u) after
the base plan; if B was trivialized in the first step, that descendant was
eliminated on B's constant there instead, and the count is unchanged). If
the descendant of B just counted is the descendant of A itself, it
receives constants on both inputs — A's, from the trivialized A, and
s(u₀) through B — so it is trivialized, and its own descendant is
eliminated on a constant instead; the count is unchanged. At least seven
gates in total, and ΔΦ ≤ +1 (A is created in the first step and destroyed
in the second; no pack forms through y, a 1-variable after the base plan).
Hence Δμ ≥ 7 + 2α_I − α_φ = 130/5 ≥ 2δ.

*5.2.β3: A is a descendant of B, and B has a descendant other than A.*
The second step eliminates A, a descendant of A, C, and another
descendant of B; a
coincidence between the two descendants is discharged as in 5.2.β2. At
least seven gates in total, ΔΦ ≤ +1, and Δμ ≥ 7 + 2α_I − α_φ = 130/5
≥ 2δ.

*5.2.β4: A is the only descendant of B (B is a 1-gate).* Here the uniform
second substitution guarantees only six gates (G, D, B, then C, A, and a
descendant of A), which at Δi = −2 is not enough; use the source paper's
operation for this configuration instead, and do not run the base plan
first. The other input of A is a variable v: v is a 2-variable, since A
becomes troubled, and v is unprotected, by Case 1; v is neither x nor y,
for a variable feeding G among three gates is a Case 3 configuration.
Substitute the constant v₀ controlling for A. Then A is trivialized,
whatever signal B passes; a descendant of A and v's other descendant V
each receive a constant and are eliminated (by Rule 3 under (T2) if one
of them is the output gate; neither is trivialized at the output, by
(T1)); and B, whose only descendant was A, is a 0-gate and is removed by
Rule 1. At least three gates are eliminated in this step (four, unless
A's descendant and V coincide), one influential input is killed (v), and
ΔΦ ≤ +2 in the minimal configuration (at most one pack through each of
the two passed signals). Now
substitute y₀, controlling for G: G is trivialized, and x dies with it,
since x's only descendants were G and B; y dies as well. C is eliminated —
trivialized on the constants it receives, or degenerating and passing
s(u) — and D, if still present, is eliminated with it (trivialized if
G's constant output is controlling for it, else degenerating and passing
s(u)). If D was already eliminated in the first step (the case v = u,
where V is D), C is trivialized on G's output directly, and a descendant
of C is eliminated on C's constant output instead; and if the constant
chosen for A was also controlling for D, then C was already eliminated in
the first step (fed by D's constant output and by y), and the second step
eliminates G and that descendant of C. In total at least six gates are
eliminated over the two steps and Δi = −3 (v, y, x). The potential grows
by at most +2 in the minimal configuration (one pack through each signal
passed in the first step, and one through the descendant of C in the
second, against G's destruction); every further pack forms through a
variable passed by a gate eliminated in the same step beyond the minimal
count, and every such further step contributes β ≥ 1 − 2α_φ ≥ 0 to Δμ by
Lemma 5.3 — with one exception: if D, eliminated in the second step, has
a further descendant beyond B and C, that descendant can become troubled
when D passes s(u), a pack accompanied by no further elimination, raising
the bound to ΔΦ ≤ +3 in that corner. Absorbing it, the floor is

    Δμ ≥ 6 + 3α_I − 3α_φ = 171/5 ≥ 2δ

over t = 2 substitutions.

Assume henceforth that C is ∧-type.

*Sub-case 5.2.α: C is an ∧-type gate.* Choose the constant for C, not for
A: s is a nonconstant affine function of u, hence takes both values, so
choose u₀ with s(u₀) the controlling value of C's signal input; substitute
u ← u₀. Then C is trivialized and eliminated, and y loses its last
descendant. The gate A, fed by the now-constant signal s(u₀), is eliminated
as well — trivialized if s(u₀) is controlling for A, degenerated, passing
its other input, if not. In the minimal configuration the second step
eliminates A, C, and a descendant of C (degenerating on C's constant output,
under (T2) if it is the output gate) and, in the trivialized shape, also a
descendant of A (degenerating on A's constant output; A is not the output
gate, by (T1), and being useful it has a descendant). If that descendant of A is
also the descendant of C, it receives constants from both gates, is
eliminated and counted once, and the trivialized shape has the
degenerated shape's count and ledger. Every further gate fed
by a u-dependent signal, and every further descendant of C, is eliminated
in the same step; charge each such extra elimination (+1) against the at
most one pack introduced by the variable it passes (−α_φ), a net
≥ 1 − α_φ > 0, so the minimal configuration is the worst case.
Influential inputs killed: Δi = −3 — x (first step), u (substituted),
y (here). Potential: G is destroyed, and A is created in the first step and
destroyed in the second; the passed variables of the second step are the
other input of C's descendant and — in the degenerated shape — A's
other input, or — in the trivialized shape — the other input of A's
descendant, each contributing ≤ 1 pack; worst case ΔΦ ≤ +1. Hence in the
trivialized shape (≥ 7 gates in total)

    Δμ ≥ 7 + 3α_I − α_φ = 178/5 ≥ 2δ,

and in the degenerated shape (≥ 6 gates in total)

    Δμ ≥ 6 + 3α_I − α_φ = 173/5 ≥ 2δ,

over t = 2 substitutions. If A = C, the two descendant counts of the
trivialized shape coincide, and the potential bound improves to ΔΦ ≤ 0
(G is destroyed; C is created in the first step and destroyed in the
second; the only passed variable of the second step is the other input
of C's descendant); that shape then gives Δμ ≥ 5 + 3α_I = 169/5, and
the bound stands. (If the base plan introduces a second troubled gate,
it too is fed by a u-dependent signal, so it is eliminated in the second
step; the counts are unchanged and the bounds stand.)

**5.3: C feeds D and D feeds B (vacuous).** Run the base plan (constant to
x; eliminate G, D, B in this order). The only nodes that could become
troubled, or feed a gate that becomes troubled, are y and the introduced
constant c; neither can become or create a troubled gate. Hence none is
introduced, contradicting the standing assumption of Case 5; the
configuration falls back to the base leaf, Δμ = 3 + α_φ + α_I.

**5.4: the remaining topology.** By 5.3, D is disconnected (both
directions) from one of B, C; by symmetry, from B. Write I_1 for B's other
input and I_2 for D's other input (D is fed by G and I_2). Substitute a
constant to x to trivialize G. A newly introduced troubled gate is either
I_1 or I_2, or is fed by I_1 or I_2 after normalization: I_1, I_2 cannot
themselves become troubled, their out-degrees not decreasing unless their
own inputs are involved. Case 5.4.1: one fed by I_2; Case 5.4.2: one fed
by I_1.

### Case 5.4.1: a new troubled gate E fed by I_2

I_2 must be a variable z: a gate I_2 could not feed a *newly* troubled gate
in this position. E is an ∧-type 1-gate; B, D disconnected forces D itself
to feed E before normalization.

**5.4.1.1: z is a 2-variable.** Then D is a 1-gate (else E not troubled).
E's other input is a variable t, or B — B, after x is substituted,
degenerates and passes its other input t to E (B fed by x and t; x
unprotected by Case 1).

*5.4.1.1.1: E's other input is a variable t.* Substitute a constant to t
(not to x), trivializing E; eliminate E, D, G and E's descendant in this
order; t = 1. ΔΦ ≤ 1 (only the other input of E's descendant can introduce
a troubled gate), Δi = −1 (t). Δμ ≥ 4 − α_φ + α_I, stronger than the
4 − 4α_φ shape.

*5.4.1.1.2: E's other input is B, B ∧-type.* Substitute constants to t and
y, trivializing B and G; x becomes a 0-variable; t = 2. Δi = −3, so
Δμ ≥ (3/2)α_I per substitution.

*5.4.1.1.3: E's other input is B, B ⊕-type.* Affine substitution
x ← t ⊕ c fixes B's output to a constant; replace B by it. With c chosen
so E is also trivialized, D and G are removed by Rule 1 and E's descendant
degenerates by Rule 3; t = 1. Five gates eliminated. Only t and the other
input of E's descendant can introduce troubled gates (y, z are 2-variables
and cannot), so ΔΦ ≤ 2; Δi = −1 (x). Δμ ≥ 5 − 2α_φ + α_I.

For the rest of 5.4.1, z is a 1-variable; then D is a 2-gate (else E not
troubled).

**5.4.1.2: D is an ∧-type 2-gate.** Substitute a constant to z,
trivializing D; eliminate D, G and D's two descendants in order; t = 1.
Eliminating G causes no potential increment; overall ΔΦ ≤ 3; Δi = −1 (z).
Δμ ≥ 4 − 3α_φ + α_I.

**5.4.1.3: z is protected.** Constant to x (G trivialized; B, D
degenerate), then constant to z (E trivialized; E's descendant and the
copy of z passed through D degenerate); t = 2. Six gates removed in total
(the paper states the count of 6 without itemizing; this count is what the
accounting uses), all by Rules 2 and 3, each removal introducing ≤ 1
troubled gate: ΔΦ ≤ 6. Δq = −1 (z's quadratic equation), Δi = −2 (x, z).
The paper's accounting decomposes the two substitutions as a Case-1-shaped
step and a Case-3-shaped step:

    Δμ ≥ [ (2 − 2α_φ + α_I + α_Q) + (4 − 4α_φ + α_I) ] / 2 per substitution,

a blend of the slack term and the 4 − 4α_φ term (12.9 at the optimum,
above δ).

**5.4.1.4: z unprotected, D an ⊕-type 2-gate — main line realizes
5 − α_Q.** Assume first that eliminating D's descendant (step (4) below)
introduces no troubled gate; the contrary is 5.4.1.4.1/5.4.1.4.2. Apply
ONE quadratic substitution

    z ← ((x ⊕ c_1) ∧ (y ⊕ c_2)) ⊕ c_3,

constants chosen so D's output becomes a constant b insensitive to its
inputs, also trivializing E; t = 1; Δq = +1. Normalize in this order:
(1) replace D by b (ΔΦ = 0); (2) eliminate G by Rule 1 (ΔΦ = 0; x, y
become 1-variables); (3) eliminate E by Rule 2 (ΔΦ = 0: E's other input is
B, fed by the 1-variable x and hence not troubled, or a variable t that
becomes a 1-variable); (4) eliminate D's descendant by Rule 3 (ΔΦ = 0, by
assumption); (5) eliminate E's descendant by Rule 3 (ΔΦ ≤ 1). Five gates
eliminated. Troubled G destroyed (−1) cancels step (5): net ΔΦ ≤ 0.
Δi = −1 (z only). Hence

    Δμ ≥ 5 − α_Q + α_I.

### The contrary of 5.4.1.4: eliminating D's descendant introduces a troubled gate

Step (4) introduces a troubled gate A. Then F, D's descendant eliminated
in step (4), is fed directly by some variable u, and A is fed directly by
F and a 2+-variable v (once F is eliminated in step (4), A is fed by the
signal F passes — a function of u — and v; this is the form used in
5.4.1.4.2). By the initial assumption of Case 5, a variable t either feeds
E directly, or is passed to E while B degenerates. 5.4.1.4.1 is the
latter; 5.4.1.4.2 the former.

The coincidence F = C (hence u = y) is discharged before the split.
(F ≠ B is already forced by the topology of 5.4: F is fed by D, while D
is disconnected from B in both directions.) Suppose F = C. Then step (4) of
the main line introduces no troubled gate after all, so the main line's
accounting applies verbatim and the coincidence never reaches the split:
when C is eliminated at step (4), D has already been replaced by the
constant b (step (1)) and G already removed (step (2)), so C is the only
descendant y has left, and every gate whose inputs change at step (4)
receives only affine functions of the 1-variable y; a gate fed by a
1-variable is not troubled (compare the opening of 5.4.1: a gate fed by the
1-variable z cannot be troubled). The degree change in the coincidence is
exactly this drop of y to a 1-variable, and no protection hypothesis is
consumed by it: y is unprotected by the standing hypothesis of Case 5, and
no subcase operation is invoked. Henceforth F ≠ C and u ≠ y, and the
pairwise distinctness assumed in 5.4.1.4.2 is in force.

**5.4.1.4.1: t feeds B and B feeds E.**

*5.4.1.4.1.1: B is a 1-gate.* Substitute constants to y and z so G, D, E
are all trivialized; B's out-degree drops to 0, B removed by Rule 1;
t = 2. Δi = −3 (x becomes a 0-variable; y, z substituted), so
Δμ ≥ (3/2)α_I per substitution. Henceforth B is a 2-gate and t a
1-variable.

*5.4.1.4.1.2: B is an ∧-type 2-gate.* Then t, a 1-variable, is unprotected
(Case 1). Substitute a constant to x, trivializing B and making t
non-influential; t = 1, Δi = −2, Δμ ≥ 2α_I.

*5.4.1.4.1.3: B is an ⊕-type 2-gate, t unprotected.* Constant to y
trivializes G (x becomes a 1-variable); then affine substitution x ← t
trivializes B; t = 2, Δi = −3 (x, y, t), Δμ ≥ (3/2)α_I per substitution.

*5.4.1.4.1.4: B is an ⊕-type 2-gate, t protected — realizes α_I/3.* The
couple of t is neither x nor y (both unprotected). Substitute a constant
to the couple of t, making t unprotected, then apply 5.4.1.4.1.3; t = 3,
Δi = −4 (couple of t, y, x, t). Δμ ≥ (4/3)α_I per substitution, i.e.
per-substitution Δμ − α_I = α_I/3.

**5.4.1.4.2: t feeds E directly.** After the quadratic substitution of
5.4.1.4, A is fed by u and v; x, y are 1-variables; u, v, x, y pairwise
distinct (so F ≠ B, C) and u, v ≠ z.

*5.4.1.4.2.1: u = t.* Then u is a 2-variable feeding E and F, hence
unprotected (Case 1). Substitute a constant to x, trivializing G; remove
B and D by Rule 3; then E and F are each fed by two 2-variables among
z, u, so substitute a constant to z or u, or apply the affine substitution u ← z ⊕ c,
making both E and F output input-independent functions; t = 2,
Δi = −3 (x, z, u), Δμ ≥ (3/2)α_I per substitution. Henceforth u ≠ t.

*5.4.1.4.2.2: F is a 1-gate.* Substitute a constant to v, replacing
troubled A by a constant; F becomes a 0-gate, removed by Rule 1; a
descendant of A is removed by Rule 3; another descendant of v is removed
by Rule 3 (v has out-degree ≥ 2). ΔΦ ≤ 2 in total: trivializing A
contributes 0 (v constant by then, F fed by the ⊕-type D), removing F
contributes 0 (D ⊕-type, u becomes a 1-variable), the two descendant
eliminations ≤ 1 each. Δi = −1 (v): Δμ ≥ 4 − 2α_φ + α_I for this
substitution.

*Pair accounting for 5.4.1.4.2.2 (t = 2).* Theorem 6.1 requires the two
substitutions together to drop μ by ≥ 2δ; account for the pair jointly,
each eliminated gate listed once: D, G, E, E's descendant (first
substitution), and F, A, A's descendant, v's descendant (second) — 8
gates; Δq = +1 (the quadratic equation of 5.4.1.4). Potential events: G,
troubled, destroyed (−1); A created by the first substitution's step (4)
(+1), destroyed by the second (−1); E's-descendant elimination ≤ +1; the
two Rule-3 removals of the second substitution ≤ +1 each; F's removal 0
(shown above). Net ΔΦ ≤ +2. Influential count: z dies at the first
substitution, v at the second; moreover u dies: u feeds exactly F and A in
this configuration — the content of the remark above that removing F makes
u a 1-variable (before that removal u feeds F and the still-standing A,
and nothing else) — so once A is replaced by a constant and F removed, u
is a 0-variable, hence non-influential: Δi = −3 over the pair. The pair
drops μ by

    Δμ ≥ 8 − α_Q + 3α_I − 2α_φ   (over t = 2 substitutions),

far above 2δ at the parameters of `proof.md` §7; this leaf is not binding.

*5.4.1.4.2.3: F is an ∧-type 2-gate.* (F is exactly a 2-gate: 2+-gate
originally, feeds none of G, D, E, and must be exactly 2 for A to be
troubled; u is a 1-variable, u ≠ t.) Substitute constants to u and t,
trivializing E and F; D becomes a 0-gate, removed by Rule 1, making the
unprotected z a 0-gate, hence non-influential; t = 2, Δi = −3 (u, t, z),
Δμ ≥ (3/2)α_I per substitution.

*5.4.1.4.2.4: F is an ⊕-type 2-gate, u unprotected.* Substitute a constant
to t, trivializing E (D becomes a 1-gate). Since
F = D ⊕ u ⊕ c_1 = G ⊕ z ⊕ u ⊕ c_2, the affine substitution z ← u ⊕ c
makes F = G ⊕ c ⊕ c_2; rewire all descendants of F to be fed by G
directly, removing D and F; z, u become 0-variables; t = 2,
Δi = −3 (t, z, u), Δμ ≥ (3/2)α_I per substitution.

*5.4.1.4.2.5: F is an ⊕-type 2-gate, u protected — realizes α_I/3.* The
couple of u is none of x, y, z, t, v (none protected). Substitute a
constant to the couple of u, making u unprotected, then apply 5.4.1.4.2.4;
t = 3, Δi = −4 (couple of u, t, z, u). Δμ ≥ (4/3)α_I per substitution.

### Case 5.4.2: the new troubled gate E is fed by I_1

E is fed not by a variable from D but by I_1 = t, the input passed from B;
its other input must be a variable z.

**5.4.2.1: B is a 2+-gate — realizes 4 − 4α_φ.** Then B is ⊕-type
(Case 3). Apply the affine substitution x ← t ⊕ c (x unprotected,
Case 1), fixing B's output regardless of its inputs; replace B by a
constant, c chosen so E is also trivialized; rewire G to be fed by t
instead of x, disconnecting x — this introduces no troubled gates; t = 1.
Four gates eliminated (B replaced by a constant; E and two descendants
removed by Rules 2 and 3). Total ΔΦ ≤ 4: the other three removals ≤ 1 each
by the normalization-lemma accounting, and replacing B by a constant ≤ 1,
by this inspection: the only gates whose inputs change are the descendants
of B and G (rewired from x to t). A descendant of B gains a constant input,
so cannot become troubled; it can at most degenerate and pass its other
input through. G was already the troubled gate of the standing
configuration and remains fed by two variables (t replacing x), so is not
a *new* troubled gate. Hence a new troubled gate can arise only downstream
of a degeneration of a descendant of B, mediated by the single variable t
those degenerations can pass; ≤ 1 such gate (or one pack) results.
Δi = −1 (x):

    Δμ ≥ 4 − 4α_φ + α_I.

**5.4.2.2: B is a 1-gate.** Substitute constants to z and y, trivializing
both E and G; B becomes a 0-gate, removed by Rule 1; t = 2,
Δi = −3 (z, y, and x, which loses both its descendants G and B), so
Δμ ≥ (3/2)α_I per substitution.

*Standing result after Case 5:* no ∧-type gate is directly fed by two
variables. Henceforth G is a topologically minimal ∧-type gate, fed by
I_1 and I_2, each computed by a fair cyclic xor-circuit (or a variable);
the xor-reconstruction operation of `proof.md` §3 applies to I_1, I_2.

## Case 6: substituting a constant to a protected variable

x protected; by Cases 1–2 a 1-variable feeding an ⊕-type gate P.
Substitute x ← d and normalize; at least P is eliminated by Rule 3. The
nominal outcome — one gate eliminated, x non-influential, x's quadratic
equation killed, no potential change:

    Δμ = 1 + α_I + α_Q.

Note 1 + α_I + α_Q = 12.4 < δ: alone it would not suffice. Case 6 rules
out every failure of the nominal outcome, leaving the Case-6 step below,
reused as the first substitution of several Case-8 leaves.

**6.1: a second gate Q is eliminated during normalization — realizes
2 − 2α_φ + α_Q.** Q is not eliminated by Rule 1: if P is the output gate,
Q becomes the output gate after P's removal, which Rule 1 cannot remove;
otherwise Q's out-degree never decreases when P is removed by Rule 3, so
Rule 1 cannot apply to Q either. Two gates eliminated (P, Q), ΔΦ ≤ 2 by
the normalization lemma, Δi = −1 (x), Δq = −1. Hence
Δμ ≥ 2 − 2α_φ + α_Q + α_I, the slack term.

For the rest of Case 6, P is the ONLY gate eliminated and ΔΦ = 1; we show
this avoidable. Let R be P's other input besides x; to produce the
increment, R is an ∧-type gate or a variable.

**6.2.1: R is an ∧-type gate, troubled after normalization (vacuous).** R
receives no new inputs from Rule 3 and its feeding variables' out-degrees
do not decrease, so R was already troubled — no increment. If R were newly
troubled, its out-degree must decrease, forcing P to be the output gate
and R a 2-gate; but then R's other descendant cannot reach the output gate
P (else a cycle containing an ∧-type gate), contradicting Case 0.2. So
ΔΦ = 1 is impossible here.

**6.2.2: R = u is a variable.** Its out-degree does not decrease during
Rule 3, so u can feed a troubled gate after normalization only as a
"2−-variable" in the paper's sense: a 1-variable feeding P whose
out-degree rises to 2 through P's Rule-3 elimination (P's descendant
becomes fed by u), so u then feeds two gates, one possibly troubled.
Subcases split on u's protection status before the substitution.

*6.2.2.1: u an unprotected 1-variable.* Let y be the couple of x.
Substitute a constant to y (making x unprotected), then the affine
substitution x ← u, trivializing P; t = 2, Δi = −3 (x, y, u),
Δμ ≥ (3/2)α_I per substitution.

*6.2.2.2: u a protected 1-variable coupled to x.* Affine substitution
x ← u, trivializing P; t = 1, Δi = −2, Δμ ≥ 2α_I.

*6.2.2.3: u a protected 1-variable with couple v ≠ x — realizes α_I/3.*
Let y be the couple of x. Substitute constants to y and v (making both x
and u unprotected), then x ← u, trivializing P; t = 3, Δi = −4 (x, y, u,
v). Δμ ≥ (4/3)α_I per substitution, attaining the α_I/3 term.

*6.2.2.4: u is a 2-variable.* Then P is a 1-gate. Let A be u's other
descendant, B the descendant of P; A ≠ B (else P useless). As in 6.2.1, A
cannot be the newly troubled gate, so B is: B is an ∧-type gate fed by P
and a variable t. Substitute a constant to t, trivializing B (then P can
be removed by Rule 1), then a constant to the couple of x (making x
unprotected); t = 2, Δi = −3 (x, t, couple of x), Δμ ≥ (3/2)α_I per
substitution.

*Standing outcome of Case 6:* apart from the configurations above,
substituting a constant to a protected variable eliminates exactly one
gate, introduces no troubled gate, and gives Δμ = 1 + α_I + α_Q: the
**Case-6 step** used below.

## Case 7: I_1 is a 2-variable, I_2 a 1-gate — realizes 4 − 4α_φ

I_1 = x a 2-variable and I_2 = Q a 1-gate, both feeding the minimal
∧-type gate G. x is unprotected (Case 1) and does not feed Q (else Q
useless). Substitute a constant to x, trivializing G; t = 1. Four gates
eliminated (G, Q, two descendants of G and x); ΔΦ ≤ 4 itemized: Q by
Rule 1 (≤ 2); G trivialized, neither input becoming or feeding a troubled
gate (0); the other two by Rule 3 (≤ 1 each). Δi = −1 (x):

    Δμ ≥ 4 − 4α_φ + α_I,

attaining the 4 − 4α_φ term.

## Case 8: both inputs of G are gates, or a 2-gate remains

By symmetry Q := I_2 is a gate (variable/variable and 2-variable/1-gate
are Cases 5 and 7) and P := I_1 a gate or a variable, of out-degree ≤
that of Q if a gate. By minimality of G, Q computes an affine function of
its inputs.

### Case 8.1: Q depends on an unprotected variable

By Proposition 2.7 (path property of fair xor-circuits), a path exists
from some unprotected variable x to Q. Substitute an appropriate constant
to Q through x via xor-reconstruction (t = 1): x becomes non-influential
and G is trivialized.

**8.1.1: Q is a 2+-gate — realizes 4 − 4α_φ.** The reconstruction replaces
the constant output b of Q, which feeds ≥ 3 gates: Q is removed and one
new gate Z is added; normalization then removes the three descendants of
b and ≥ 1 descendant of G (Rules 2 and 3) — 4 gates net. The
reconstruction changes no potential (the empty packing persists); the
removals give ΔΦ ≤ 4; Δi = −1 (x). Δμ ≥ 4 − 4α_φ + α_I.

**8.1.2: Q is a 1-gate and P = t is a variable.** Xor-reconstruction on Q;
both the unprotected 1-variable t and x become non-influential
(Δi = −2), so Δμ ≥ 2α_I.

**8.1.3: P and Q are both 1-gates — realizes 4 − 4α_φ.** Q computes an
affine function. Xor-reconstruction replaces Q's output by a constant b of
out-degree 2; normalization removes G, P, and the other two descendants of
b and G — four distinct gates. Distinctness: b's two descendants differ
from G and from each other (a coincidence closes a cycle through the
cyclic xor-circuit), and P is not among them — fed by b, P would,
unwinding the reconstruction, feed Q in the original circuit, excluded
since P, Q are the two inputs of the minimal gate G, neither feeding the
other. Itemized: P removed by Rule 1 ≤ 2; G contributes 0; the other two
descendants ≤ 2 in total; so ΔΦ ≤ 4. Δi = −1 (x): Δμ ≥ 4 − 4α_φ + α_I.

*Standing condition for the rest of Case 8:* Q — and P as well, if an
⊕-type gate of equal out-degree — depends only on protected variables.

### Case 8.2: Q depends only on protected variables

Let I be the set of protected variables Q depends on.

**8.2.1: two coupled protected variables feed Q.** x_j, x_k protected,
coupled to each other, both feeding Q; both are 1-variables (Case 1).
Apply the affine substitution x_j ← x_k ⊕ c, making Q's output constant;
t = 1, Δi = −2 (x_j, x_k), so Δμ ≥ 2α_I.

### Case 8.2.2: Q is a 2+-gate

Fix x_j ∈ I, x_k its couple. First substitution: x_k ← d plus
normalization; by Case 6 exactly one gate (a descendant of x_k) is
eliminated, Δμ = 1 + α_I + α_Q, and x_j becomes unprotected.

**8.2.2.1: Q is not fed by x_k — realizes (5 − 2α_φ + α_Q)/2.** After
x_k ← d, Q remains a 2+-gate and depends on the unprotected x_j. Second
substitution: xor-reconstruction with x_j, trivializing G. The
reconstruction degenerates three gates A, B, C and trivializes G: 4 gates
eliminated in this step, ΔΦ ≤ 2 (A ≤ 1, C ≤ 1; eliminating B and G is
ΔΦ-free — assumed here, discharged in 8.2.5), so the second step gives
Δμ ≥ 4 − 2α_φ + α_I. Averaging the two substitutions:

    Δμ ≥ (1/2)[ (1 + α_I + α_Q) + (4 − 2α_φ + α_I) ]
       = (1/2)(5 − 2α_φ + 2α_I + α_Q) per substitution,

attaining (5 − 2α_φ + α_Q)/2: the sub-δ Case-6 step (12.4) is compensated
by the stronger second step (13.2).

**8.2.2.2: Q is fed by x_k (reduction).** After x_k ← d, the gate A that
was fed by x_k is now fed by Q's other input Q′; Q′ is a 2+-gate depending
on the unprotected x_j. Apply 8.2.2.1 with Q := Q′; the same bound holds.

*Standing condition after 8.2.2:* Q is a 1-gate; P is a variable or a
1-gate.

### Case 8.2.3: P = t is a variable

Let x_j ∈ I with couple x_k. Substitute a constant to x_k (making x_j
unprotected) and normalize — G is then fed by t and by Q′, where Q′
computes an affine function of x_j and may equal Q — then apply
xor-reconstruction on Q′, trivializing G; t = 2. Δi = −3 (t and x_j
become non-influential through the reconstruction, x_k was substituted),
so Δμ ≥ (3/2)α_I per substitution.

### Case 8.2.4: P is an ⊕-type 1-gate

If P depended on an unprotected variable, 8.1.3 would apply with P, Q
exchanged; so both P and Q depend only on protected variables.

**8.2.4.1: some x_j ∈ I has its couple x_k not feeding the other gate —
main branch realizes (5 − 2α_φ + α_Q)/2.** By symmetry treat Q and x_j,
x_k not feeding P. First substitution x_k ← d: a Case-6 step
(Δμ = 1 + α_I + α_Q), x_j becomes unprotected. Second substitution:
xor-reconstruction with x_j trivializing G, as in 8.1.3 (4 gates removed).
*Main branch:* if the second substitution has ΔΦ ≤ 2, it gives
Δμ ≥ 4 − 2α_φ + α_I, and the two-substitution average is
(1/2)(5 − 2α_φ + 2α_I + α_Q) per substitution, attaining
(5 − 2α_φ + α_Q)/2.

For the rest of 8.2.4.1, assume ΔΦ ≥ 3 in the second substitution.
Eliminating G and B (the node preceding P on the x_j → G path) is assumed
ΔΦ-free — failures are discharged in 8.2.5. The descendant of G removed by
Rule 3 contributes ≤ 1. So both inputs of P must each contribute +1, giving
four sub-configurations: (a) P fed by a gate T — impossible: T must be
⊕-type (Case 5) and cannot become troubled by degeneracy, both its
feeding variables staying influential with undiminished out-degree;
(b) P fed by two protected variables — impossible: by Case 1 they are
1-variables and neither feed nor become troubled gates; (c) variables
t_1, t_2 are passed to P by the degeneracy of B and of the eliminated
descendant of x_k — 8.2.4.1.2; (d) t_1 directly feeds P, and t_2 is passed
via B and via P′ (the descendant of x_k), making t_2 a 3-variable feeding
an ∧-type gate T troubled after P's elimination — 8.2.4.1.1 and 8.2.4.1.3.

*8.2.4.1.1: t_1 (or t_2) is protected.* Add a third substitution
t_1 ← c (or t_2 ← c), eliminating ≥ 3 gates via Rules 2 and 3 while
killing a quadratic equation; this step alone gives
Δμ ≥ 3 − 3α_φ + α_Q + α_I. Together the three substitutions give
Δμ ≥ (1/3)(8 − 8α_φ + 3α_I + 2α_Q) per substitution (≈ 12.93 at the
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

**8.2.4.2: P fed by two non-coupled protected variables.** P is fed by
x_j, x_k, protected, not coupled to each other; both are 1-variables
(Case 1) feeding only P. Substitute appropriate constants to both,
trivializing P and G; Q and ≥ 1 descendant of G are also eliminated:
4 gates (P, G, Q, a descendant of G), ΔΦ ≤ 3 (0 for P and G; 2 for Q, one
per input; 1 for the descendant). Two quadratic equations killed
(Δq = −2), Δi = −2, so over t = 2 substitutions
Δμ ≥ (1/2)(4 − 3α_φ + 2α_I + 2α_Q) per substitution (13.1, above δ).

**8.2.4.3 (vacuous).** By 8.2.4.2, P and Q are each fed by ≤ 1 variable.
Suppose P depends on a protected x_j whose couple x_k feeds Q. Then Q's
other input is a gate R (by 8.2.4.2), non-constant (Case 0.3). There is no
path from x_k to R (Q is fed by x_k and R), so R does not depend on x_k;
but Q depends on x_k and on some variable x_i that R depends on, and
x_i ≠ x_j (else Q depends on the coupled pair x_j, x_k, against the
assumption). Hence Q depends on two non-coupled protected variables, so P
is fed by two protected variables, contradicting the assumption of 8.2.4.2
under which this case arises — impossible.

### Case 8.2.5: discharging the ΔΦ-freeness assumption of 8.2.2.1 and 8.2.4.1

8.2.2.1 and 8.2.4.1 assumed that eliminating G and B cannot increase the
potential; we now discharge this. Let x_j ∈ I with couple x_k; a path
x_j → Q exists (Proposition 2.7). Substitute x_k ← d and normalize: x_j
becomes unprotected, and by Case 6 only the descendant of x_k is
eliminated, with no new troubled gate. If x_k feeds Q, its other input Q′
is passed to feed G (Q′ is ⊕-type, identified with Q). Xor-reconstruction
with x_j trivializes G. Let B be the node preceding Q on the x_j → Q path.
Eliminate G first, then B, then a descendant of G; no other elimination in
this process can increase the potential. The cases split on whether
eliminating G (8.2.5.1) or eliminating B (8.2.5.2 and its reductions) gives
ΔΦ = +1.

**8.2.5.1: eliminating G increases the potential by 1.** Then P must be a
variable x_m whose out-degree dropped from 3 to 2 through G's elimination;
after x_k ← d, x_m is a 3-variable.

*8.2.5.1.1 (vacuous): x_m was already a 3-variable before x_k ← d.* Then
x_m cannot feed G before the substitution (Case 3), so it was passed to G
via the degeneracy of P′ (the descendant of x_k) when x_k ← d. For the
potential to increase, x_m must feed an ∧-type gate T right before G's
elimination, hence right after x_k ← d. Feeding T before the substitution
is impossible (Case 3); feeding T via P′'s degeneracy would require P′ to
be a 2+-gate, but T becoming troubled needs x_m to be a 2-variable after
G's elimination, forcing P′ to be only a 1-gate — contradiction.

*8.2.5.1.2: x_m was a 1-variable before x_k ← d — realizes α_I/3.* Then
x_m and x_k feed the same ⊕-type gate P′. Work in the circuit *before*
the x_k ← d substitution (after it, x_k is a constant and x_m a
3-variable, so the analogy is only available there): P′ is fed by
protected x_k (whose couple is x_j) and the 1-variable x_m — exactly the
configuration of Case 6.2.2 with (x, u) := (x_k, x_m). All three subcases
of 6.2.2 apply: x_m unprotected — the 6.2.2.1 maneuver (constant to the
couple x_j, then the affine substitution x_k ← x_m), t = 2, Δi = −3;
x_m = x_j, i.e. coupled to x_k — the 6.2.2.2 affine substitution
x_k ← x_m, t = 1, Δi = −2; x_m protected with a couple different from
x_k — the 6.2.2.3 maneuver (constants to both couples, then x_k ← x_m),
t = 3, Δi = −4 — the read-once property of rdq-sources makes the four
killed inputs distinct. The worst case is the last, so
Δμ ≥ (4/3)α_I per substitution, attaining the α_I/3 term.

*8.2.5.1.3: x_m was a 2-variable originally, raised to 3 by x_k ← d.*
Then P′ is an ⊕-type 2-gate fed by the protected x_k and x_m; x_m either
feeds G directly or is passed to G through P′'s degeneracy.

— *8.2.5.1.3.1: x_m feeds G directly, and P′ feeds an ∧-type gate T before
x_k ← d.* Substitute constants to x_m and x_k, trivializing P′, T and G;
their descendants are degenerated by normalization; t = 2. Six gates are
eliminated via Rules 2 and 3, ΔΦ ≤ 6 (≤ 1 per removed gate); one quadratic
equation is killed (Δq = −1, x_k's); Δi = −2. Hence
Δμ ≥ (6 − 6α_φ + α_Q + 2α_I)/2 per substitution — the average of a
Case-3-shaped and a Case-1-shaped decrement (12.9, above δ). Henceforth
x_m is passed to G through P′'s degeneracy.

— *8.2.5.1.3.2 (reduction): Q is a 2-gate and P′ is a 2+-gate depending on
the unprotected x_m.* Apply Case 8.1 with P′ and Q exchanged; the bound is
whichever Case 8.1 provides for that configuration, in particular ≥ δ.

— *8.2.5.1.3.3: Q is a 3+-gate.* First substitution x_k ← d eliminates P′
(Q becomes an ⊕-type gate depending on the unprotected x_j): a Case-6
step, Δμ ≥ 1 + α_Q + α_I. Second substitution: xor-reconstruction removing
five gates by Rules 2 and 3 (three descendants of Q, the first gate on the
reversed x_j → Q path, a descendant of G), ΔΦ ≤ 5, so
Δμ ≥ 5 − 5α_φ + α_I; t = 2 in total. The average
(1/2)(6 − 5α_φ + 2α_I + α_Q) ≥ (1/2)(6 − 6α_φ + 2α_I + α_Q) per
substitution (≥ 12.9) is above δ.

**8.2.5.2: eliminating B increases the potential by 1.** B feeds Q already
before x_k ← d. After the reconstruction, B's other input besides the
constant b is a node D that becomes (or feeds) a troubled gate; D feeds Q
before the reconstruction, originally either directly, or via P′ (the
descendant of x_k), passed to B during normalization after x_k ← d. D
cannot be an ∧-type gate (minimality of G), so D is a variable x_ℓ.

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

— *8.2.5.2.1.5 (vacuous): t is a 2-variable.* For T to become troubled, P′
must then be a 1-gate; but B depends on the unprotected x_ℓ, impossible
after Case 8.1.

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
Δμ ≥ (7 − 3α_φ + 3α_I + 2α_Q)/3 per substitution over t = 3 (≈ 12.93 at
the optimum, above δ; the displayed blend
(1/3)(7 − 4α_φ + 3α_I + 2α_Q) ≈ 12.87 is weaker and also suffices).

**8.2.5.3 (reduction).** B feeds B′ and B′ feeds Q originally, and B′ is
the gate eliminated by x_k ← d. 8.2.5.2's analysis still applies: B is fed
by a variable x_ℓ right after the reconstruction and x_ℓ feeds Q
originally; if x_ℓ is unprotected, Q does not depend on x_ℓ, and
Q = B ⊕ x_k ⊕ x_ℓ ⊕ c forces B to depend on x_ℓ, with B not a 0-gate
before elimination (otherwise B = x_j); if x_ℓ is protected, substitute
x_ℓ ← c after the whole process, mirroring step (3) of 8.2.5.2.2. The
bounds of 8.2.5.2 are inherited.

**8.2.5.4 (reduction).** B feeds Q′ and Q′ feeds Q; Q is eliminated by the
substitution x_k ← d and Q′ takes the place of Q. 8.2.5.2's analysis
applies verbatim with Q′ in the role of Q, and its bounds are inherited.

This exhausts the configurations: every step of the elimination procedure
falls under one of Cases 0–8, and in each leaf the per-substitution
decrement of μ is at least

    δ = α_I + min{ α_I/3, 2 − 2α_φ + α_Q, 4 − 4α_φ, 3 + α_φ, 5 − α_Q,
                   (5 − 2α_φ + α_Q)/2 },

the six quantities attained respectively by the α_I/3 leaves (5.4.1.4.1.4,
5.4.1.4.2.5, 6.2.2.3, 8.2.5.1.2, 8.2.5.2.1.4), Cases 1 and 6.1, the
4 − 4α_φ leaves (3, 5.4.2.1, 7, 8.1.1, 8.1.3), the base leaf of Case 5,
Case 5.4.1.4, and Cases 8.2.2.1 and 8.2.4.1; pair elimination (Case 0.4)
and normalization are never worse. This proves Theorem 6.1. ∎
