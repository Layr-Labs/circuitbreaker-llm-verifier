# Case analysis for Theorem 6.1: elimination

Notation as in `proof.md`. C computes f | R, f an affine disperser for d,
dim R > 2d + 2, packing P. Step: t in {1, 2, 3} subst. (t = 0 in Case 0),
then normalize; Dmu per subst. unless stated; each leaf needs Dmu >= δ.
Priority order; each case assumes earlier configs. absent. Lemma 5.3
silent throughout - the blanket justification for every per-leaf DPhi bound
not separately argued: norm. never hurts (β >= 1 - 2ap >= 0 for R1,
β >= 1 - ap otherwise) and raises Φ only via gates fed by the elim.
gate's inputs. aI = 48/5, aQ = 9/5, ap = 1/5, delta = 64/5. Signs: Dmu = decrease of μ; Di/Dq negative when an input/equation is killed; DPhi positive when potential grows. c->x: constant
subst. to x; R1-R5: `proof.md` §5's rules.

## Case 0: preliminaries (only Dmu >= 0 needed)

**0.1** C unnormalized: normalize. Dmu >= 0.
**0.2** A gate in the acyclic part not reaching the output; topologically maximal G is a
0-gate; R1: DPhi <= 2, Dmu >= 1 - 2ap >= 0. Now all gates reach the
output.
**0.3** G outputs a fixed constant c: replace by c; DPhi <= 2,
Dmu >= 1 - 2ap >= 0. G depending only on I_1: rewire to (I_1, 0)
(DPhi <= 1), then R3; again DPhi <= 2, Dmu >= 1 - 2ap >= 0. Now
every gate depends on both inputs.
**0.4** P != {}: Lemma 5.4 gives decrement >= δ_p = min{ (3/2)aI,
3 - 3ap + aI + aQ } = 69/5 > δ. Now P = {}; Φ counts TGs directly.

## Case 1: a protected variable feeds two gates, or an &-type gate

c->(that var); t = 1: Dq = -1, Di = -1, >= 2 gates degen. (R2, R3),
DPhi <= 2. Dmu >= 2 - 2ap + aI + aQ (realizes 2 - 2ap + aQ; 17/5 at the
optimum, slack). *Standing:* every prot. var feeds exactly one
xor-type gate.

## Case 2: a protected 0-variable

p prot., feeding no gate; q its couple. c->q; t = 1: equation affine,
Di = -2. Dmu >= 2aI >> δ. *Standing:* every prot. var is a
1-var feeding a xor-type gate.

## Case 3: a variable feeding an &-type gate, total out-degree >= 4 - realizes 4 - 4ap

x feeds &-gate G, out-deg(x) + out-deg(G) >= 4. c->x triv. G; t = 1:
>= 3 desc. of x and G degen.; a gate fed by both (counted twice)
is triv., its desc. degen. - >= 4 gates (R2, R3);
DPhi <= 4, Di = -1. Dmu >= 4 - 4ap + aI. *Standing:* vars feeding
&-gates have outdeg 1 or 2; such gates are 1- or 2-gates.

## Case 4: an &-type gate fed by two variables, one a 1-variable

G fed by x (1-var), y, both unprot. (Case 1). c->y triv. G;
t = 1: Di = -2, DPhi <= 0 (only gates fed by x or y change inputs; x's single desc. is G). Dmu >= 2aI. *Standing:* an &-gate fed by two
vars is fed by two 2-vars (Case 5) or falls to Case 8.

## Case 5: an &-type gate fed by two 2-variables (the troubled configuration)

G fed by 2-vars x, y; by Case 3 a 1-gate, hence troubled. B, C: other
desc. of x, y; D: desc. of G; B, C != D (normalized); x, y
unprot. (Case 1). Base plan: c->x triv. G; B, D, G removed.
**Base leaf (3 + ap).** No TG introduced: gates 3; DPhi = -1 (G troubled);
Di = -1 (x). Exactly Dmu = 3 + ap + aI. Cases 5.1-5.4 assume one IS
introduced.

**5.1: B = C.** Truth-table enumeration over the &/xor choices for G, B
(only the four (x, y) assignments matter) gives c->x, c->y, or x <- y + c
making G and B both output constants; t = 1, Di = -2. Dmu >= 2aI.


**5.2: D feeds both B and C.** Base plan (c->x; elim. G, D, B in order).
By Lemma 5.3 a TG arises only at a gate fed by D's other input I; I a
gate: outdeg undiminished, cannot become troubled - base leaf. So I = var
u. Only desc. of D, B change inputs; C receives the nonconstant affine
function of u D passes on degen.; each B-descendant B's signal, likewise
a nonconstant affine function of u (B constant-passing => its desc.
elim. on constants in step 1, adding gates there). Hence *every* new TG
is u-fed; let A be one. (Paper's operation - constants to the other input
of the desc. below B and to y - covers only u feeding exactly two gates,
TG below B; ours uniform: C itself troubled, further desc. of D, B fed
by x and D in any position.)
B not the output: output reachable from C (Case 0), but any C->B path
enters B via x or D, D via G or u - all entry points vars, which no gate
feeds. So B has a desc., elim. by R3 without the output exception; then
u-signals feed C and - unless B constant-passed (its desc. already elim.
in step 1, adding gates) - >= 1 desc. of B, plus possibly further desc.
of D.
Second subst.: constant to u at its controlling value for &-gate A
(Case 3's move on (u, A); u prot. -> Dq = -1, only helps); t = 2. Step 2
elim.: A (triv.); every other u-fed gate (>= 1), degen. on its
now-constant input (R3, or (T2) if the output); one desc. of A, degen.
on A's constant output - A not the output (T1), and useful, so has a
desc. If that desc. is itself a counted u-fed gate, take its desc.
instead: not the output either (both inputs constant this step, vs (T1)),
so has a desc., which degen. and is elim. Total >= 6 gates (>= 3 in
step 2; >= 4 in step 1 if constant-passing).
Di = -3: x (desc. G, B elim.), u (substituted), y (desc. G, C; C u-fed
after base plan, triv./degen. in step 2).
Potential: G destroyed; all introduced TGs u-fed, destroyed in step 2. A
final gate is newly troubled only if fed by a var passed in step 2
(Lemma 5.3); gates introduced via one passed var pairwise share it =>
<= 1 pack among them (R1 packing observation, `proof.md` §5). Worst case
(u-signals feed exactly two gates): passed vars are two (other input of
the degen. u-fed gate; other input of A's desc.), so DPhi <= -1 + 2 = +1;
each further u-fed gate adds an elim. gate against <= 1 further passed
var; constant-passing case: one passed var only - estimate only improves.

    Dmu >= 6 + 3*aI - ap = 6 + 144/5 - 1/5 = 173/5 >= 2*delta

over t = 2. (Second TG from base plan: likewise u-fed, triv./degen. in
step 2; count unchanged, bound stands.)

**5.3: C feeds D and D feeds B (vacuous).** Base plan: only y and the
introduced constant c could become troubled, or feed one; neither can.
None is introduced - against Case 5's assumption; base leaf applies,
Dmu = 3 + ap + aI.

**5.4: remaining topology.** By 5.3, D is disconnected (both directions)
from one of B, C; by symmetry from B. I_1: B's other input; I_2: D's other
input. c->x triv. G. A new TG is I_1 or I_2, or fed by one after
norm. (I_1, I_2 cannot become troubled themselves: outdegs
undiminished). 5.4.1: fed by I_2; 5.4.2: fed by I_1.

### Case 5.4.1: a new troubled gate E fed by I_2

I_2 is a var z (a gate could not feed a new TG here). E is an &-type
1-gate; B, D disconnected forces D to feed E before norm..

**5.4.1.1: z a 2-variable.** Then D is a 1-gate (else E not troubled). E's
other input: a var t, or B (after c->x, B degen. and passes its
other input t; B fed by x, t).
*5.4.1.1.1: other input t.* c->t triv. E; eliminate E, D, G, E's
desc. in order; t = 1: DPhi <= 1, Di = -1 (t). Dmu >= 4 - ap + aI.
*5.4.1.1.2: other input B, &-type.* c->t, c->y triv. B, G; x a
0-var; t = 2, Di = -3. Dmu >= (3/2)aI.
*5.4.1.1.3: other input B, xor-type.* x <- t + c fixes B's output constant,
c also triv. E; replace B by it; D, G by R1, E's desc.
degen. (R3); t = 1, gates 5: DPhi <= 2 (only t and the other input of E's desc. can introduce TGs; y, z are 2-vars and cannot), Di = -1 (x).
Dmu >= 5 - 2ap + aI.
Now z is a 1-var; then D is a 2-gate (else E not troubled).

**5.4.1.2: D an &-type 2-gate.** c->z triv. D; eliminate D, G, D's
two desc.; t = 1: DPhi <= 3 (eliminating G causes no potential increment), Di = -1 (z). Dmu >= 4 - 3ap + aI.

**5.4.1.3: z protected.** c->x (G triv.; B, D degen.), then c->z
(E triv.; E's desc. and the copy of z through D degen.);
t = 2, gates 6 (the paper's count, unitemized; R2, R3, each <= +1
troubled): DPhi <= 6, Dq = -1, Di = -2 (x, z). The paper's accounting
(Case-1 step plus Case-3 step): Dmu >= [ (2 - 2ap + aI + aQ) + (4 - 4ap + aI) ] / 2
(12.9 > δ).

**5.4.1.4: z unprot., D a xor-type 2-gate - main line realizes
5 - aQ.** Assume step (4) introduces no TG; the contrary is 5.4.1.4.1/.2.
ONE quadratic subst. z <- ((x + c_1) & (y + c_2)) + c_3, constants making
D's output a constant b insensitive to its inputs, also triv. E;
t = 1, Dq = +1. Normalize in order: (1) D replaced by b (DPhi = 0);
(2) G by R1 (DPhi = 0; x, y become 1-vars); (3) E by R2 (DPhi = 0:
E's other input is B, fed by 1-var x, or a var t becoming a
1-var); (4) D's desc. by R3 (DPhi = 0, by assumption); (5) E's
desc. by R3 (DPhi <= 1). Gates 5; G destroyed (-1) cancels (5): net
DPhi <= 0. Di = -1 (z). Dmu >= 5 - aQ + aI.

### The contrary of 5.4.1.4: step (4) introduces a troubled gate

Step (4) introduces troubled A. Then F, D's desc. elim. in (4),
is fed directly by a var u, and A directly by F and a 2+-var v
(once F is elim., A is fed by F's passed signal - a function of u -
and v). A var t feeds E directly (5.4.1.4.2) or is passed to E as B
degen. (5.4.1.4.1).

**5.4.1.4.1: t feeds B and B feeds E.**
*5.4.1.4.1.1: B a 1-gate.* c->y, c->z triv. G, D, E; B's outdeg
drops to 0 (R1); t = 2, Di = -3 (x a 0-var; y, z). Dmu >= (3/2)aI.
Now B is a 2-gate, t a 1-var.
*5.4.1.4.1.2: B an &-type 2-gate.* t unprot. (Case 1). c->x triv.
B, t non-influential; t = 1, Di = -2. Dmu >= 2aI.
*5.4.1.4.1.3: B a xor-type 2-gate, t unprotected.* c->y triv. G
(x a 1-var); then x <- t triv. B; t = 2, Di = -3. Dmu >= (3/2)aI.
*5.4.1.4.1.4: B a xor-type 2-gate, t protected (aI/3).* Couple of t is
neither x nor y. c->(couple of t) unprotects t; then 5.4.1.4.1.3; t = 3,
Di = -4. Dmu >= (4/3)aI, i.e. Dmu - aI = aI/3.


**5.4.1.4.2: t feeds E directly.** After 5.4.1.4's subst., A is fed by
u, v; x, y are 1-vars; u, v, x, y pairwise distinct (F != B, C),
u, v != z.
*5.4.1.4.2.1: u = t.* u a 2-var feeding E, F, hence unprot.
(Case 1). c->x triv. G; B, D by R3; then E, F are each fed by two
2-vars among z, u, so c->z, c->u, or u <- z + c makes both
input-independent; t = 2, Di = -3 (x, z, u). Dmu >= (3/2)aI. Now
u != t.
*5.4.1.4.2.2: F a 1-gate.* c->v replaces troubled A by a constant; F
becomes a 0-gate (R1); a desc. of A by R3; another desc. of v by
R3 (out-deg(v) >= 2). DPhi <= 2: triv. A 0 (v constant by then, F
fed by xor-type D), F's removal 0 (D xor-type, u becomes a 1-var),
the two desc. removals <= 1 each. Di = -1 (v): Dmu >= 4 - 2ap + aI for
this subst.
*Pair accounting (t = 2).* The pair must drop μ by >= 2δ; account jointly,
each gate once: D, G, E, E's desc. (subst. 1); F, A, A's desc., v's desc.
(subst. 2) - 8 gates; Dq = +1 (5.4.1.4's quadratic eq.). Potential: G,
troubled, destroyed (-1); A created at step (4) (+1), destroyed at
subst. 2 (-1); E's-desc. elim. <= +1; the two R3 removals of subst. 2
<= +1 each; F's removal 0 (above). Net DPhi <= +2. Influential: z dies at
subst. 1, v at subst. 2; moreover u dies - u feeds exactly F and A here:
removing F makes u a 1-var, and once A is constant and F removed, u is a
0-var, non-influential: Di = -3 over the pair. Total: Dmu >= 8 - aQ + 3aI
- 2ap over t = 2 - far above 2δ; not binding.
*5.4.1.4.2.3: F an &-type 2-gate* (exactly a 2-gate; u a 1-var,
u != t). c->u, c->t triv. E, F; D a 0-gate (R1), z a 0-gate; t = 2,
Di = -3 (u, t, z). Dmu >= (3/2)aI.
*5.4.1.4.2.4: F a xor-type 2-gate, u unprotected.* c->t triv. E
(D a 1-gate). Since F = D + u + c_1 = G + z + u + c_2, z <- u + c makes
F = G + c + c_2; rewire F's desc. to G, removing D, F; t = 2,
Di = -3 (t, z, u). Dmu >= (3/2)aI.
*5.4.1.4.2.5: F a xor-type 2-gate, u protected (aI/3).* Couple of u is
none of x, y, z, t, v. c->(couple of u) unprotects u; then 5.4.1.4.2.4;
t = 3, Di = -4. Dmu >= (4/3)aI.

### Case 5.4.2: the new troubled gate E is fed by I_1

E is fed by I_1 = t, B's passed input; its other input is a var z.
**5.4.2.1: B a 2+-gate (4 - 4ap).** B xor-type (Case 3). x <- t + c
(x unprot.) fixes B's output; replace B by a constant, c also
triv. E; rewire G to be fed by t instead of x; no TGs introduced;
t = 1, gates 4 (B replaced; E and two desc. by R2, R3). DPhi <= 4:
the three removals <= 1 each (Lemma 5.3); replacing B <= 1 by inspection -
only desc. of B and G change inputs; a desc. of B gains a
constant input, so can only degen. and pass its other input; G remains
fed by two vars (t for x), so is not *new*; a new TG arises only
downstream of such a degeneration, mediated by the single passed var
t: <= 1 (or one pack). Di = -1 (x). Dmu >= 4 - 4ap + aI.
**5.4.2.2: B a 1-gate.** c->z, c->y triv. E, G; B a 0-gate (R1);
t = 2, Di = -3 (z, y, x - losing both desc.). Dmu >= (3/2)aI.

*Standing after Case 5:* no &-gate is directly fed by two vars.
Now G is a topologically minimal &-gate fed by I_1, I_2, each a
fair cyclic xor-circuit or a var; xor-reconstruction (`proof.md` §3)
applies to I_1, I_2.

## Case 6: substituting a constant to a protected variable

x prot.; by Cases 1-2 a 1-var feeding a xor-type gate P. x <- d,
normalize; at least P elim. (R3). Nominal outcome - 1 gate, x
non-influential, x's quadratic equation killed, DPhi = 0:
Dmu = 1 + aI + aQ = 12.4 < δ, insufficient alone. Case 6 rules out every
failure of the nominal outcome.

**6.1: a second gate Q is eliminated (2 - 2ap + aQ).** Q is not elim. by
R1: if P is the output gate, Q becomes it after P's removal, which R1
cannot remove; otherwise Q's outdeg never decreases when P is removed by
R3. Gates 2 (P, Q); DPhi <= 2; Di = -1 (x), Dq = -1.
Dmu >= 2 - 2ap + aQ + aI, the slack term. Now P is the ONLY gate elim.
and DPhi = 1; we show this avoidable. R := P's other input - an &-gate
or a var.

**6.2.1: R an &-type gate, troubled after normalization (vacuous).** R
receives no new inputs from R3, its feeders' outdegs do not decrease, so R
was already troubled - no increment. If R were newly troubled, its outdeg
must decrease, forcing P the output gate and R a 2-gate; but then R's
other desc. cannot reach output gate P (else a cycle containing an
&-gate), contradicting Case 0.2. So DPhi = 1 impossible.

**6.2.2: R = u a variable.** Its outdeg does not decrease during R3, so u
can feed a TG after norm. only as a "2--var": a 1-var feeding P whose
outdeg rises to 2 through P's R3 elim. (P's desc. becomes fed by u); u
then feeds two gates, one possibly troubled. Split on u's protection
status.
*6.2.2.1: u an unprotected 1-variable.* y := couple of x. c->y (x becomes
unprot.), then x <- u triv. P; t = 2, Di = -3 (x, y, u). Dmu >= (3/2)aI.
*6.2.2.2: u a protected 1-variable coupled to x.* x <- u triv. P; t = 1,
Di = -2. Dmu >= 2aI.
*6.2.2.3: u a protected 1-variable, couple v != x (aI/3).* y := couple of
x. c->y, c->v (x, u both unprot.), then x <- u triv. P; t = 3, Di = -4
(x, y, u, v). Dmu >= (4/3)aI, attaining aI/3.
*6.2.2.4: u a 2-variable.* Then P is a 1-gate. A := u's other desc.,
B := P's desc.; A != B (else P useless). As in 6.2.1, A cannot be the new
TG; B is: B &-type, fed by P and a var t. c->t triv. B (P then falls by
R1), then c->(couple of x); t = 2, Di = -3 (x, t, couple of x).
Dmu >= (3/2)aI.
*Standing outcome of Case 6:* apart from the configs. above, a constant
to a prot. var elim. exactly one gate, introduces no TG, and gives
Dmu = 1 + aI + aQ: the **Case-6 step** below.

## Case 7: I_1 is a 2-variable, I_2 a 1-gate - realizes 4 - 4ap

I_1 = x a 2-var, I_2 = Q a 1-gate, both feeding minimal &-type G.
x unprot. (Case 1), not feeding Q (else Q useless). c->x triv.
G; t = 1. Gates 4 (G, Q, two desc. of G and x); DPhi <= 4 itemized:
Q by R1 (<= 2); G triv., neither input becoming or feeding a TG (0);
other two by R3 (<= 1 each). Di = -1 (x). Dmu >= 4 - 4ap + aI.

## Case 8: both inputs of G are gates, or a 2-gate remains

By symmetry Q := I_2 is a gate (var/var and 2-var/1-gate
are Cases 5 and 7); P := I_1 a gate or var, outdeg <= Q's if a
gate. By minimality of G, Q computes an affine function of its inputs.

### Case 8.1: Q depends on an unprotected variable

By Prop 2.7 a path runs from some unprot. var x to Q. Substitute
a constant to Q through x via xor-reconstruction (t = 1): x
non-influential, G triv.
**8.1.1: Q a 2+-gate (4 - 4ap).** Q's constant output b feeds >= 3 gates:
Q removed, one new gate Z added; norm. removes b's three
desc. and >= 1 desc. of G (R2, R3) - 4 gates net. Reconstruction
changes no potential; DPhi <= 4; Di = -1 (x). Dmu >= 4 - 4ap + aI.
**8.1.2: Q a 1-gate, P = t a variable.** Xor-reconstruction on Q;
Di = -2 (t, x). Dmu >= 2aI.
**8.1.3: P, Q both 1-gates (4 - 4ap).** Xor-reconstruction replaces Q's
output by a constant b of outdeg 2; norm. removes G, P, and
the other two desc. of b and G - 4 distinct gates: b's desc.
differ from G and each other (a coincidence closes a cycle through the
cyclic xor-circuit); P is not among them (fed by b, P would feed Q
originally; P, Q are G's two inputs, neither feeding the other). DPhi <= 4
(P by R1 <= 2; G 0; other two <= 2 total). Di = -1 (x).
Dmu >= 4 - 4ap + aI.
*Standing, rest of Case 8:* Q - and P too, if a xor-type gate of equal
outdeg - depends only on prot. vars.

### Case 8.2: Q depends only on protected variables

I := the prot. vars Q depends on.
**8.2.1: two coupled protected variables feed Q.** x_j, x_k prot.,
coupled, both 1-vars (Case 1). x_j <- x_k + c makes Q's output
constant; t = 1, Di = -2 (x_j, x_k). Dmu >= 2aI.

### Case 8.2.2: Q is a 2+-gate

Fix x_j in I, x_k its couple. Subst. 1: x_k <- d + norm.; by
Case 6 exactly one gate (a desc. of x_k) is elim.,
Dmu = 1 + aI + aQ; x_j becomes unprot.
**8.2.2.1: Q not fed by x_k ((5 - 2ap + aQ)/2).** Q remains a 2+-gate
depending on unprot. x_j. Subst. 2: xor-reconstruction with x_j
triv. G; it degen. three gates A, B, C and triv. G:
4 gates, DPhi <= 2 (A <= 1, C <= 1; elim. B and G DPhi-free -
assumed here, discharged in 8.2.5), so Dmu >= 4 - 2ap + aI. Average:
Dmu >= (1/2)[ (1 + aI + aQ) + (4 - 2ap + aI) ]
= (1/2)(5 - 2ap + 2aI + aQ), attaining (5 - 2ap + aQ)/2: the sub-δ
Case-6 step (12.4) is compensated by the stronger second step (13.2).
**8.2.2.2: Q fed by x_k (reduction).** After x_k <- d, gate A fed by x_k
is now fed by Q's other input Q'; Q' is a 2+-gate depending on
unprot. x_j. Apply 8.2.2.1 with Q := Q'; same bound.
*Standing after 8.2.2:* Q is a 1-gate; P is a var or a 1-gate.

### Case 8.2.3: P = t is a variable

x_j in I, couple x_k. c->x_k (x_j unprot.) + normalize - G then fed by
t and by Q', Q' an affine function of x_j, possibly Q itself - then
xor-reconstruction on Q' triv. G; t = 2. Di = -3 (t, x_j killed by
the reconstruction, x_k substituted). Dmu >= (3/2)aI.


### Case 8.2.4: P is a xor-type 1-gate

If P depended on an unprot. var, 8.1.3 applies with P, Q exchanged; so
both P, Q depend only on prot. vars.
**8.2.4.1: some x_j in I has couple x_k not feeding the other gate - main
branch ((5 - 2ap + aQ)/2).** By symmetry treat Q, x_k not feeding P.
Subst. 1: x_k <- d, a Case-6 step (Dmu = 1 + aI + aQ); x_j unprot.
Subst. 2: xor-reconstruction with x_j triv. G, as in 8.1.3 (4 gates
removed). *Main branch:* DPhi <= 2 in subst. 2 gives Dmu >= 4 - 2ap + aI;
average (1/2)(5 - 2ap + 2aI + aQ), attaining (5 - 2ap + aQ)/2. Rest:
DPhi >= 3 in subst. 2. Elim. G and B (node preceding P on the x_j -> G
path) assumed DPhi-free (failures: 8.2.5); G's R3 desc. <= 1. So both
inputs of P contribute +1 each - four sub-configs: (a) P fed by a gate T:
impossible - T xor-type (Case 5), cannot become troubled by degen., both
feeders staying influential, outdegs undiminished; (b) P fed by two prot.
vars: impossible - 1-vars (Case 1), feeding and becoming no TGs; (c) t_1,
t_2 passed to P by the degen. of B and of x_k's elim. desc. - 8.2.4.1.2;
(d) t_1 feeds P directly, t_2 passed via B and via P' (x_k's desc.): t_2
a 3-var feeding an &-gate T troubled after P's elim. - 8.2.4.1.1/.3.
*8.2.4.1.1: t_1 (or t_2) protected.* Subst. 3: t_1 <- c (or t_2 <- c):
>= 3 gates (R2, R3), one quadratic equation killed; this step alone
Dmu >= 3 - 3ap + aQ + aI. Together Dmu >= (1/3)(8 - 8ap + 3aI + 2aQ)
(~ 12.93 > δ).
*8.2.4.1.2 (vacuous): t_1, t_2 both passed to P in degeneracy* (t_2 only
via B, t_1 only via P'). Then t_2 feeds Q, B feeds P before x_k <- d; B a
gate (B = x_j impossible: x_j a prot. 2-var beforehand);
Q = B + t_2 + c_1 forces B to depend on t_2; then
P = x_k + t_1 + B + c_2 depends on unprot. t_2 - against Case 8.2's
standing condition.
*8.2.4.1.3 (vacuous): t_1 feeds P directly before x_k <- d.* If t_2 also
feeds P directly or only via P', P depends on unprot. t_1 (Case 8.1:
impossible); if t_2 is passed by B, then as in 8.2.4.1.2, P or Q depends
on unprot. t_1 - impossible by Case 8.1.
*Standing consequence:* P depends only on vars whose couples feed Q, and
vice versa.
**8.2.4.2: P fed by two non-coupled protected variables.** x_j, x_k
prot., not coupled, both 1-vars (Case 1) feeding only P. Constants to
both triv. P and G; Q and >= 1 desc. of G also elim.: 4 gates (P, G, Q, a
desc. of G), DPhi <= 3 (0 for P, G; 2 for Q, one per input; 1 for the
desc.). Dq = -2, Di = -2, t = 2: Dmu >= (1/2)(4 - 3ap + 2aI + 2aQ)
(13.1 > δ).
**8.2.4.3 (vacuous).** By 8.2.4.2, P, Q each fed by <= 1 var. Suppose P
depends on prot. x_j whose couple x_k feeds Q. Q's other input is a gate
R (by 8.2.4.2), non-constant (Case 0.3). No path x_k -> R (Q fed by x_k,
R), so R does not depend on x_k; but Q depends on x_k and on some x_i
that R depends on, x_i != x_j (else Q depends on coupled x_j, x_k,
against assumption). Hence Q depends on two non-coupled prot. vars, so P
is fed by two prot. vars - contradicting the 8.2.4.2 assumption this case
arises under. Impossible.

### Case 8.2.5: discharging the DPhi-freeness assumption of 8.2.2.1 and 8.2.4.1

8.2.2.1 and 8.2.4.1 assumed elim. G and B DPhi-free; discharged here.
x_j in I, couple x_k; path x_j -> Q exists (Prop 2.7). Subst. x_k <- d,
normalize: x_j unprot.; by Case 6 only x_k's desc. elim., no new TG.
If x_k feeds Q, its other input Q' is passed to feed G (Q' xor-type,
identified with Q). Xor-reconstruction with x_j triv. G. B := node
preceding Q on the x_j -> Q path. Elim. G, then B, then a desc. of G; no
other elim. here raises the potential. Split: elim. G raises it (8.2.5.1)
or elim. B does (8.2.5.2 and reductions).
**8.2.5.1: eliminating G raises the potential by 1.** Then P is a var
x_m, outdeg dropped 3 -> 2 through G's elim.; after x_k <- d, x_m is a
3-var.
*8.2.5.1.1 (vacuous): x_m already a 3-variable before x_k <- d.* x_m
cannot feed G beforehand (Case 3), so was passed to G via P''s degen.
(P' := x_k's desc.) at x_k <- d. The rise needs x_m feeding an &-gate T
just before G's elim., hence just after x_k <- d: feeding T beforehand is
impossible (Case 3); feeding T via P''s degen. needs P' a 2+-gate, but T
troubled needs x_m a 2-var after G's elim., forcing P' a 1-gate -
contradiction.
*8.2.5.1.2: x_m a 1-variable before x_k <- d (aI/3).* x_m, x_k feed the
same xor-type gate P'. Work *before* x_k <- d (after it x_k is constant,
x_m a 3-var): P' fed by prot. x_k (couple x_j) and 1-var x_m -
Case 6.2.2's config., (x, u) := (x_k, x_m). By 6.2.2: x_m unprot. -
.1 maneuver, t = 2, Di = -3; x_m = x_j - .2, x_k <- x_m, t = 1, Di = -2;
x_m prot., other couple - .3 maneuver, t = 3, Di = -4 (four killed
inputs distinct, read-once). Worst: Dmu >= (4/3)aI, attaining aI/3.
*8.2.5.1.3: x_m a 2-variable originally, raised to 3 by x_k <- d.* P' a
xor-type 2-gate fed by prot. x_k and x_m; x_m feeds G directly or is
passed to G through P''s degen.
- *8.2.5.1.3.1: x_m feeds G directly; P' feeds an &-gate T before
x_k <- d.* Constants to x_m, x_k triv. P', T, G; desc. degen.; t = 2:
six gates (R2, R3), DPhi <= 6 (<= 1 each), Dq = -1 (x_k's equation),
Di = -2. Dmu >= (6 - 6ap + aQ + 2aI)/2 (12.9 > δ). Now x_m passed
to G through P''s degen.
- *8.2.5.1.3.2 (reduction): Q a 2-gate, P' a 2+-gate depending on
unprot. x_m.* Case 8.1 with P', Q exchanged; bound applies (>= δ).
- *8.2.5.1.3.3: Q a 3+-gate.* Subst. 1: x_k <- d elim. P' (Q xor-type,
depending on unprot. x_j): Case-6 step, Dmu >= 1 + aQ + aI. Subst. 2:
xor-reconstruction removing five gates (R2, R3: three desc. of Q, first
gate on reversed x_j -> Q path, a desc. of G), DPhi <= 5, so
Dmu >= 5 - 5ap + aI; t = 2. Average (1/2)(6 - 5ap + 2aI + aQ) >=
(1/2)(6 - 6ap + 2aI + aQ) (>= 12.9 > δ).
**8.2.5.2: eliminating B raises the potential by 1.** B feeds Q already
before x_k <- d. B's other input besides constant b is a node D becoming
(or feeding) a TG; D feeds Q beforehand, directly or via P' (x_k's
desc.), passed to B in norm. after x_k <- d. D not &-type
(minimality of G): D is a var x_l.
*8.2.5.2.1: x_l unprotected.* Q then does not depend on x_l (Case 8.1);
Q fed by B, x_l (or B, x_l, x_k) gives Q = B + x_l + c (or + x_k + c), so
B depends on x_l.
- *8.2.5.2.1.1 (vacuous): B degenerates to a 0-gate just before
elimination.* Its outdeg fell, so B is the reconstruction's added gate Z,
i.e. B = x_j beforehand; then Q is fed by x_j, x_l and depends on both,
yet B = x_j cannot depend on x_l - contradiction. Now B a 1+-gate;
x_l passed to new TG T at B's elim.; T's other input a var t.
- *8.2.5.2.1.2 (vacuous): B and t both feed T before x_k <- d.* T a
minimal &-gate fed by xor-type B, B depending on unprot. x_l -
Case 8.1 excludes.
- *8.2.5.2.1.3 (vacuous): B passed to T as P' degenerates, t feeds T
originally.* x_k prot. => P' xor-type; T a minimal &-gate fed by P',
P' depending on unprot. x_l - Case 8.1 excludes. Now t passed
to T as P' degen., B feeds T originally.
- *8.2.5.2.1.4: t a 1-variable (aI/3).* P' xor-type, fed by two 1-vars:
prot. x_k (couple x_j) and t - Case 6.2.2's config., (x, u) := (x_k, t)
(the analogy needs x_k prot.). Subcases as 8.2.5.1.2: t unprot., .1
maneuver - t = 2, Di = -3; t = x_j, .2 - t = 1, Di = -2; t prot., other
couple, .3 maneuver - t = 3, Di = -4 (read-once). Worst: Dmu >= (4/3)aI,
attaining aI/3.
- *8.2.5.2.1.5 (vacuous): t a 2-variable.* T troubled forces P' a 1-gate;
but B depends on unprot. x_l, impossible after Case 8.1.
*8.2.5.2.2: x_l protected* - a prot. 2-var feeding an &-gate (a truly
troubled gate, or a pack's inner var). Three subst. (1) x_k <- d: Case-6
step, Dmu >= 1 + aI + aQ. (2) Xor-reconstruction removing B, G, a desc. of
G: G's elim. DPhi-free - P != B (else G falls by input subst., not
reconstruction); P cannot make G troubled (P a var: unprot. 1-var or
2+-var; P a gate: path P -> Q closes a cycle G -> P -> Q -> G through an
&-gate once P, Q degen.) - so DPhi <= 2, Dmu >= 3 - 2ap + aI.
(3) x_l <- c: three gates removed (R2, R3); the troubled gate fed by x_l
goes first, both inputs 2-vars (no increment); the other two add <= 2;
destroying it (or breaking the pack) subtracts 1: DPhi <= 1; Dq = -1, so
Dmu >= 3 - ap + aI + aQ. Total, t = 3: Dmu >= (7 - 3ap + 3aI + 2aQ)/3
(~ 12.93 > δ; weaker blend (1/3)(7 - 4ap + 3aI + 2aQ) ~ 12.87 also
suffices).
**8.2.5.3 (reduction).** B feeds B' and B' feeds Q originally; B' is the
gate elim. by x_k <- d. 8.2.5.2's analysis applies: B is fed by a var
x_l right after the reconstruction, x_l feeds Q originally; if x_l
unprot., Q does not depend on x_l, and Q = B + x_k + x_l + c forces B
to depend on x_l, B not a 0-gate before elim. (else B = x_j); if x_l
prot., substitute x_l <- c after the whole process, mirroring step (3)
of 8.2.5.2.2. Bounds inherited.
**8.2.5.4 (reduction).** B feeds Q', Q' feeds Q; Q elim. by x_k <- d and
Q' takes Q's place. 8.2.5.2's analysis applies verbatim with Q' in Q's
role; bounds inherited.

This exhausts the configurations: every step of the elimination procedure
falls under one of Cases 0-8, and in each leaf the per-subst. decrement
of μ is at least

    delta = aI + min{ aI/3, 2 - 2ap + aQ, 4 - 4ap, 3 + ap, 5 - aQ,
                      (5 - 2ap + aQ)/2 },

the six quantities attained respectively by the aI/3 leaves (5.4.1.4.1.4,
5.4.1.4.2.5, 6.2.2.3, 8.2.5.1.2, 8.2.5.2.1.4), Cases 1 and 6.1, the
4 - 4ap leaves (3, 5.4.2.1, 7, 8.1.1, 8.1.3), the base leaf of Case 5,
Case 5.4.1.4, and Cases 8.2.2.1 and 8.2.4.1; pair elimination (Case 0.4)
and normalization are never worse. This proves Theorem 6.1. ∎
