## 7. Gate elimination: the full case table

Elimination Lemma (Li-Yang Section 4, statement). Let C be a
normalized circuit computing f|R for an rdq-source R, where f is
the disperser of Section 4. Then either C can be simplified
without a substitution (measure strictly decreases), or there are
1 <= s <= 3 substitutions with measure decrement at least delta
per substitution, where delta = 64/5 at the published point, and
with the parameters alpha_I, alpha_Q, alpha_phi satisfying the
six inequalities of Section 8. The proof in Li-Yang Section 4 is
a case analysis on a topologically minimal AND gate: its inputs
are classified by whether they are variables or gates, by
out-degree, by protected/unprotected status, and by whether the
local structure is a troubled gate, a troubled pair, or one of
the bottleneck configurations in the paper's Figures. Each case
gives an explicit substitution (constant, affine, or quadratic,
or xor-reconstruction to a gate) and a count of gates removed,
influential inputs killed, quadratic equations killed, and
potential change, from which Delta mu >= delta is read off.

The case table proving that lemma, following Li-Yang Section 4
in its case order, is given in Sections 7.1-7.9 below. Conventions
used throughout, in addition to Sections 5-6: Delta mu =
mu_before - mu_after, so a removed gate contributes +1, a killed
influential input +alpha_I, a killed quadratic equation
+alpha_Q, an added quadratic equation -alpha_Q, and a potential
change Delta Phi contributes -alpha_phi * Delta Phi. A
k-variable / k-gate has out-degree k; a 2+-gate has out-degree
at least 2. The couple of a protected variable is the other
variable in its quadratic equation. Cases are taken in order:
in a later case, earlier cases are assumed not to apply.
At the published point alpha_I = 48/5, alpha_Q = 9/5,
alpha_phi = 1/5, delta = 64/5; every arithmetic line below is
evaluated at that point, and "per substitution" means the
total Delta mu divided by the number of substitutions in the
block (1 to 3), as in the lemma statement. Normalization
bounds from Section 6 are used without repeating them:
Rule 1 gives Delta mu >= 1 - 2*alpha_phi per removed gate,
Rules 2-5 give >= 1 - alpha_phi, and a newly troubled gate
can only be a removed gate's input or a gate fed by it.

### 7.1 Case 0: simplification, no substitution needed

0.1 If C is not normalized, normalize it: Delta mu >= 0 by
Section 6, since alpha_phi < 1/2. 0.2 If a gate in the acyclic
part does not reach the output, take a topologically maximal
such gate; it is a 0-gate, removed by Rule 1 with
Delta mu >= 1 - 2*alpha_phi >= 0. 0.3 If a gate computes a
constant for all inputs, replace it by that constant:
one gate removed, Delta Phi <= 2 (one per input), so
Delta mu >= 1 - 2*alpha_phi >= 0. If a gate depends only
on one input, rewire it to that input plus a constant 0,
make it degenerate, remove it by Rule 3, with the same
bound. 0.4 If the packing P is nonempty, troubled-pair
elimination (Section 6, Lemma 3.12) gives
Delta mu >= delta_p = min(3*alpha_I/2,
3 - 3*alpha_phi + alpha_I + alpha_Q) = min(72/5, 69/5)
= 69/5 >= delta per substitution. Henceforth C is
normalized, P is empty, and every gate depends on both
of its inputs.

### 7.2 Case 1: protected variable feeding two gates or an AND gate

Substitute a constant to the protected variable so that a
quadratic equation is killed and the variable becomes
non-influential, while at least two gates are removed by
Rules 2 and 3. Delta mu >= 2 - 2*alpha_phi + alpha_I +
alpha_Q = 2 - 2/5 + 48/5 + 9/5 = 65/5 = 13 >= delta.

### 7.3 Case 2: protected 0-variable

Let p be a protected 0-variable and q its couple. Substitute
a constant to q, turning the quadratic equation into an
affine one; both p and q become non-influential.
Delta mu >= 2*alpha_I = 96/5 >= delta. Henceforth every
protected variable is a 1-variable feeding an XOR-type gate.

### 7.4 Case 3: variable feeding an AND gate, total out-degree >= 4

Let x feed an AND-type gate G with out(x) + out(G) >= 4.
Substitute a constant to x trivializing G; at least three
descendants of x and G become degenerate, and if a gate is
fed by both G and x it is trivialized and its descendants
degenerate instead (no double counting loss). At least
four gates are removed by Rules 2 and 3, and x becomes
non-influential. Delta mu >= 4 - 4*alpha_phi + alpha_I
= 4 - 4/5 + 48/5 = 64/5 = delta exactly. Henceforth a
variable feeding an AND-type gate has out-degree 1 or 2.

### 7.5 Case 4: AND gate fed by two variables, one a 1-variable

Let G be fed by variables x, y, with x a 1-variable. By
Case 1, x and y are unprotected. Substitute a constant to
y trivializing G; both x and y become non-influential.
Delta mu >= 2*alpha_I >= delta.

### 7.6 Case 5: AND gate fed by two 2-variables (troubled gate)

Then G is a 1-gate, hence troubled; otherwise Case 3 would
apply (out(x) + out(G) >= 4 if out(G) >= 2). Let B, C be
the descendants of x, y respectively (the gates fed by x
and y other than G, in the paper's notation), and D the
descendant of G. By Case 1, x and y are unprotected free
variables, and by normalization B != D and C != D.
Substituting a constant to x (or y) trivializes G and
removes B (or C), D and G by normalization. If no troubled
gate is introduced, G itself is destroyed, so Phi drops
by 1, x becomes non-influential, and
Delta mu = 3 + alpha_phi + alpha_I = 3 + 1/5 + 48/5
= 64/5 = delta exactly. Assume henceforth that a troubled
gate is introduced by that substitution.
5.1 If B = C: a constant substitution to x or y, or an
affine substitution x <- y xor c, makes both G and B
constant (checked by enumerating the AND-type truth
tables); x and y become non-influential, so
Delta mu >= 2*alpha_I >= delta.
5.2 If D feeds both B and C: instead substitute constants
to z and y (z the other input of the new troubled gate E,
see below) so that E and G are trivialized, and remove B
by Rule 1. The argument for the substitution is: by the
normalization lemma, the only node that can become or
feed a troubled gate is the other input I of D; if I is
a gate its out-degree does not decrease, so it cannot
become troubled; hence I must be a variable u, B must
feed another AND-type gate E, and C != E. After the new
substitution, u is a 2-variable feeding exactly C and E
and B was a 1-gate originally. Three influential inputs
x, y, z are killed in two substitutions, so
Delta mu >= 3*alpha_I/2 = 72/5 >= delta per substitution
(plus the gates removed). Henceforth D does not feed
one of B, C.
5.3 If C feeds D and D feeds B: after substituting a
constant to x and eliminating G, D, B in order, the
normalization lemma says the only nodes that can become
or feed troubled gates are y and the other input of D,
and neither can do so in this configuration; hence no
troubled gate is introduced, contradicting the assumption
of this subcase. Henceforth D is disconnected (in both
directions) from one of B, C; by symmetry, from B.
5.4 Suppose B is fed by x and I_1, and D is fed by G and
I_2. Substitute a constant to x trivializing G. By the
normalization lemma, a newly introduced troubled gate is
either I_1 or I_2, or fed by I_1 or I_2 after normalization.
Since B and D are disconnected, out(I_1) and out(I_2)
do not decrease, so they cannot themselves be newly
troubled unless their inputs are involved; hence a new
troubled gate must be fed by I_1 or I_2.
5.4.1 The new troubled gate E is fed by I_2 = z, and D
fed E before normalization.
5.4.1.1 z is a 2-variable, D is a 1-gate, and the other
input of E is a variable t or B.
5.4.1.1.1 Other input t, a variable: substitute a constant
to t trivializing E, and eliminate E, D, G and a descendant
of E in order. Delta Phi <= 1 (only the other input of
that descendant can introduce potential), so
Delta mu >= 4 - alpha_phi + alpha_I = 4 - 1/5 + 48/5
= 67/5 >= delta.
5.4.1.1.2 Other input B, B AND-type: substitute constants
to t and y trivializing B and G, making x a 0-variable.
This block is two substitutions, not three: the initial
constant to x from 5.4 is not part of it — x becomes
non-influential here without a further substitution.
Three variables (t, y, x) become non-influential, so
Delta mu >= 3*alpha_I = 144/5 total, i.e. 72/5 per
substitution >= delta, before counting the four gates
E, D, G, B removed as well.
5.4.1.1.3 Other input B, B XOR-type: substitute
x <- t xor c fixing B to a constant (and, with c chosen
suitably, trivializing E), then remove D and G by Rule 1
and degenerate the descendant of E by Rule 3. Only t
and the other input of that descendant can introduce
troubled gates (y and z are 2-variables after the
substitution), so Delta mu >= alpha_I + 5 - 2*alpha_phi
= 48/5 + 5 - 2/5 = 71/5 >= delta.
5.4.1.2 z is a 1-variable, D is an AND-type gate:
substitute a constant to z trivializing D, and eliminate
D, G and the two descendants of D. Delta Phi <= 3, so
Delta mu >= 4 - 3*alpha_phi + alpha_I
= 4 - 3/5 + 48/5 = 65/5 = 13 >= delta.
5.4.1.3 z is protected: after substituting a constant to
x (trivializing G, degenerating B and D), substitute a
constant to z trivializing E, degenerating the descendant
of E and z. Six gates are removed, at most six troubled
gates are introduced (all removals by Rules 2 and 3),
two influential inputs are killed and one quadratic
equation is killed. Total Delta mu >= 6 - 6*alpha_phi
+ 2*alpha_I + alpha_Q = 6 - 6/5 + 96/5 + 9/5 = 129/5;
per substitution (two substitutions) this is 129/10
= 12.9 >= delta.
5.4.1.4 z is unprotected, D is an XOR-type gate: perform
the quadratic substitution
z <- ((x xor c_1) AND (y xor c_2)) xor c_3, with constants
chosen so that D's output is a constant b independent of
the inputs, trivializing E. Normalize in the order:
replace D by b (no potential increase); eliminate G by
Rule 1 (x and y become 1-variables, no increase);
eliminate E by Rule 2 (no increase, since its other input
is B or a variable t, and in the first case B is fed by
the 1-variable x, in the second t becomes a 1-variable);
eliminate the descendant of D by Rule 3; eliminate the
descendant of E by Rule 3 (Delta Phi <= 1). The originally
troubled gate G is destroyed (Phi -1), so net
Delta Phi <= 0. Five gates are removed, one quadratic
equation is added, and one influential input is killed:
z, which was a 1-variable before the substitution and is
now the quadratic variable, a 0-variable that is not
protected (it is the left-hand side of the new equation,
not a variable occurring on its right-hand side). x and y
are NOT killed here: they were 2-variables, become
1-variables, and become protected by occurring in the new
quadratic equation, so they remain influential on both
counts. Hence, with the added quadratic equation,
Delta mu >= 5 - alpha_Q + alpha_I = 5 - 9/5 + 48/5
= 64/5 = delta exactly. The remaining work in this
subcase is to rule out a potential increase from the
descendant-of-D step; that is done in 5.4.1.4.1-5.4.1.4.2
below, by showing that either that step is clean, or an
explicit alternative block gives at least delta per
substitution.
5.4.1.4.1 The variable t that feeds B (and B feeds E):
if B is a 1-gate, substitute constants to y and z
trivializing G, D and E, remove B by Rule 1, and make x
a 0-variable: three variables non-influential in two
substitutions, Delta mu >= 3*alpha_I/2 >= delta per
substitution. Otherwise B is a 2-gate, so t is a
1-variable. If B is AND-type, t is unprotected by
Case 1; substitute a constant to x trivializing B,
making t non-influential: Delta mu >= 2*alpha_I >= delta.
If B is XOR-type and t is unprotected, substitute a
constant to y trivializing G (x becomes a 1-variable),
then an affine substitution x <- t trivializing B:
three variables (x, y, t) non-influential in two
substitutions. If B is XOR-type and t is protected,
its couple cannot be x or y (unprotected), so first
substitute a constant to its couple to make t
unprotected, then apply the previous block: four
variables non-influential in three substitutions, and
with the gates removed this is at least delta per
substitution (4*alpha_I/3 = 64/5 before counting gates).
5.4.1.4.2 The variable t directly feeds E: after the
quadratic substitution, a troubled gate A introduced by
eliminating the descendant F of D is fed by u and v,
with u, v distinct from x, y (now 1-variables) and from
z (a 1-variable feeding D), so F is neither B nor C.
If u = t: u is a 2-variable feeding E and F, hence
unprotected by Case 1; substitute a constant to x
trivializing G, remove B and D by Rule 3, then E and F
are fed by exactly the two 2-variables z and u, so a
constant substitution to z or u, or an affine
substitution u <- z xor c, makes both outputs constant:
three variables (x, z, u) non-influential in two
substitutions. If F is a 1-gate: substitute a constant
to v replacing A by a constant, remove F by Rule 1 and
a descendant of A and of v by Rule 3; v is a 2+-variable,
so another descendant of v is also removed or made
trivial. Four gates are removed with Delta Phi <= 2,
so Delta mu >= 4 - 2*alpha_phi + alpha_I
= 4 - 2/5 + 48/5 = 66/5 >= delta. If F is a 2+-gate, then
to make A troubled F must be exactly a 2-gate and u a
1-variable initially. If F is AND-type and u != t:
substitute constants to u and t trivializing E and F,
D becomes a 0-gate, removing it by Rule 1 makes z a
0-variable: three variables non-influential in two
substitutions. If F is XOR-type and u is unprotected:
substitute a constant to t trivializing E (D becomes a
1-gate), then note F = D xor u xor c_1 = G xor z xor u
xor c_2, so an affine substitution z <- u xor c makes
F = G xor c xor c_2; rewire descendants of F to G,
removing D and F and making z and u 0-variables:
three variables non-influential in two substitutions.
If F is XOR-type and u is protected, its couple cannot
be x, y, z, t or v: x and y are protected as the
right-hand side of the new quadratic equation, but
read-once disjointness bars them as a couple, t is
barred by Case 1, and v is a 2+-variable, so
unprotected; so first substitute a constant to its
couple to make u unprotected, then apply the previous
block: four variables non-influential in three
substitutions.
5.4.2 The new troubled gate E is not fed by a variable
from D, hence is fed by I_1 = t from B, and its other
input is a variable z.
5.4.2.1 B is a 2+-gate, hence XOR-type by Case 3 (it
cannot be AND-type, or Case 3 would apply to x feeding
it, since x is unprotected by Case 1): substitute
x <- t xor c fixing B to a constant, replace B by a
constant, and trivialize E; E and two of its descendants
are removed by Rules 2 and 3, and replacing B by a
constant can introduce at most one potential increment
(caused by t), but the paper's bound for the whole block
is Delta Phi <= 4, so Delta mu >= 4 - 4*alpha_phi +
alpha_I = 64/5 = delta exactly.
5.4.2.2 B is a 1-gate: substitute constants to z and y
trivializing E and G, making B a 0-gate removed by
Rule 1; x, y, z become non-influential, so
Delta mu >= 3*alpha_I/2 >= delta per substitution
(two substitutions, three variables, plus gates).
This completes Case 5: every AND gate fed by two
variables is handled, so henceforth no AND-type gate is
directly fed by two variables.

### 7.7 Case 6: protected variable (preparatory block)

Let x be protected. By Cases 1 and 2, x is a 1-variable
feeding an XOR-type gate P. Try substituting a constant
x <- d and normalizing: at least P is removed by Rule 3.
If the result is not exactly one gate removed with no
new troubled gate, one of the following gives at least
delta per substitution, so we may assume that case for
later averaging in Case 8.
6.1 Another gate Q is eliminated second (P first): Q is
not eliminated by Rule 1 (if P is the output, Q becomes
the output; otherwise Q's out-degree never decreases
after P is removed), and the two eliminations give
Delta Phi <= 2, so Delta mu >= 2 - 2*alpha_phi +
alpha_Q + alpha_I = 2 - 2/5 + 9/5 + 48/5 = 65/5
>= delta.
6.2 P is the only gate eliminated and Delta Phi = 1:
let R be the other input of P apart from x.
6.2.1 R is an AND-type gate troubled after normalization:
it was originally troubled (Rule 3 gives it no new
inputs and does not decrease the out-degrees of its
feeding variables, and the only way it could become
troubled is if P is the output and R was originally a
2-gate, which would require a path from the other
descendant of R to P, impossible by Case 0.2 since that
would be a cycle containing an AND-type gate), so it
does not produce a potential increment, contradiction.
6.2.2 R = u, a variable, whose out-degree does not
decrease under Rule 3.
6.2.2.1 u unprotected 1-variable: let y be the couple
of x; substitute a constant to y making x unprotected,
then an affine substitution x <- u trivializing P:
three variables (x, y, u) non-influential in two
substitutions, Delta mu >= 3*alpha_I/2 >= delta
per substitution.
6.2.2.2 u protected 1-variable coupled with x: an affine
substitution x <- u trivializes P, both become
non-influential: Delta mu >= 2*alpha_I >= delta.
6.2.2.3 u protected 1-variable with couple v != x: let
y be the couple of x; substitute constants to y and v
making x and u unprotected, then an affine substitution
x <- u trivializing P: four variables non-influential
in three substitutions, Delta mu >= 4*alpha_I/3
= 64/5 = delta per substitution before counting gates.
6.2.2.4 u a 2-variable, P a 1-gate: let A be the other
descendant of u and B the descendant of P; A != B,
otherwise P is useless. As in 6.2.1, A cannot be the
newly troubled gate, so B is: B is an AND-type gate fed
by P and a variable t. Substitute a constant to t
trivializing B, remove P by Rule 1, and substitute a
constant to the couple of x making x unprotected:
three variables (x, t, couple of x) non-influential
in two substitutions, Delta mu >= 3*alpha_I/2 >= delta
per substitution.
Henceforth, substituting a constant to a protected
variable may be assumed to eliminate exactly one gate
with no new troubled gate, giving the building-block
decrement Delta mu = 1 + alpha_I + alpha_Q = 1 + 48/5
+ 9/5 = 62/5, which is used only in averages with a
following substitution in Case 8, never alone as a
delta block.

### 7.8 Case 7: minimal AND gate fed by a 2-variable and a 1-gate

Let I_1 = x be a 2-variable and I_2 = Q a 1-gate, with x
unprotected by Case 1 and x not feeding Q (otherwise Q
is useless). Substitute a constant to x trivializing G
and eliminate four gates in total: G, Q and two
descendants of G and x. Potential bounds: Q removed by
Rule 1 gives Delta Phi <= 2; G trivialized, with neither
input becoming nor feeding a troubled gate, gives
Delta Phi <= 0; the other two gates removed by Rule 3
give Delta Phi <= 1 each. Hence Delta Phi <= 4 and
Delta mu >= 4 - 4*alpha_phi + alpha_I = 64/5 = delta.

### 7.9 Case 8: minimal AND gate fed by gates (xor-reconstruction)

Let G be a topologically minimal AND-type gate, fed by
I_1 = P and I_2 = Q; G exists, otherwise the function
could be trivialized by a single affine substitution.
By minimality, P and Q compute affine functions, so
xor-reconstruction (Section 5, Proposition 2.8) can be
applied to them: a path from a variable to Q in a fair
cyclic xor-circuit is reversed, making Q an input and
the variable a gate, so a constant substitution to Q
is well defined. If P is a gate, assume out(P) <= out(Q),
swapping the names if needed.
8.1 Q depends on an unprotected variable x, with a path
x -> Q (Proposition 2.7). Substitute a constant to Q via
xor-reconstruction trivializing G.
8.1.1 Q a 2+-gate: after reconstruction Q is replaced by
a constant b of out-degree at least 3 (Q is removed and
a new gate Z is added, with no new troubled gate, since
the packing is empty by Case 0.4, so potential is
unchanged by the reconstruction itself). Normalizing
removes the three descendants of b and at least one
descendant of G by Rules 2 and 3, so Delta Phi <= 4 and
Delta mu >= 4 - 4*alpha_phi + alpha_I = delta.
8.1.2 Q a 1-gate and P = t a variable: by Cases 1 and 7,
t is an unprotected 1-variable. Xor-reconstruction makes
both t and x non-influential, so Delta mu >= 2*alpha_I
>= delta.
8.1.3 P and Q both 1-gates: after reconstruction Q is
replaced by a constant b of out-degree 2. Normalizing
removes G, P and the other two descendants of b and G
(four gates; P and the descendant of G do not coincide,
or there would be a cycle containing an AND-type gate,
and P and the descendant of b apart from G do not
coincide, since P's only descendant before reconstruction
was G). Potential bounds: P by Rule 1 gives <= 2, G
trivialized with both inputs eliminated gives <= 0,
the other two gates by Rules 2 or 3 give <= 2 in total,
so Delta Phi <= 4 and Delta mu >= 4 - 4*alpha_phi +
alpha_I = delta.
Henceforth, for any minimal AND gate as in Case 8 with
out(P) <= out(Q), Q depends only on protected variables,
and similarly for P if it is an XOR-type gate with the
same out-degree as Q.
8.2 Q depends only on protected variables, say in set I.
8.2.1 P or Q is fed by two coupled protected variables
x_j, x_k (both 1-variables by Cases 1 and 2): substitute
x_j <- x_k xor c trivializing Q and killing the quadratic
equation containing them; both variables become
non-influential, so Delta mu >= 2*alpha_I >= delta.
8.2.2 Q a 2+-gate: let x_j in I be protected, x_k its
couple. Substitute x_k <- d and normalize; by Case 6
exactly one gate is eliminated, with
Delta mu = 1 + alpha_I + alpha_Q = 62/5.
8.2.2.1 Q is not fed by x_k originally. Original-state facts,
all used below. (i) After the substitution x_k <- d and
normalization, exactly one gate is eliminated (Case 6 building
block), with Delta mu = 1 + alpha_I + alpha_Q = 62/5 for that
substitution; call the eliminated gate P', the descendant of
x_k. (ii) Q still feeds G and is still a 2+-gate: Q is not fed
by x_k, so Q is not P', and Q's descendants are unchanged by
step 1 except as in (iv) below. (iii) x_j is now unprotected:
the quadratic equation containing x_j and x_k is killed by
x_k <- d, and by read-once disjointness x_j occurs on the
right-hand side of no other quadratic equation. Q still depends
on x_j: substituting x_k <- d does not remove x_j from Q's
affine function. (iv) If A below is P' (A fed by x_k as well
as Q), then after step 1 Q's descendants are G together with
A's former descendants, rewired to Q. If A had a descendant
other than G, Q is still 2+ and the analysis below applies with
any such descendant as A. If A's only descendant is G, then A
is an input of G, so A = P (G's inputs are P and Q, and A != Q
since A is fed by Q); then step 1 eliminates both A and G (A's
removal rewires G to Q twice, and Rule 5 removes G), which is
Case 6.1, giving Delta mu >= 2 - 2*alpha_phi + alpha_I +
alpha_Q = 65/5 >= delta for step 1 alone, so no second
substitution is needed. If A was a 0-gate, then A is not the
output gate: if A were the output, then after step 1 Q is the
output, and the second substitution below would make the output
constant, contradicting Proposition 2.4 (two substitutions,
well within the duration); and a non-output 0-gate is removed
by Case 0.2 before Case 8. Hence A = P' with A a 0-gate cannot
occur, and Q is a 2+-gate after step 1 in all cases where a
second substitution is needed.
Definitions, in words (no figure). Let B be the preceding node
of Q on a path from x_j to Q after step 1 (such a path exists
by Proposition 2.7, since Q depends on x_j; take the path to
be shortest, so its nodes are distinct): B feeds Q, and B
is an xor/equiv gate, or B is x_j itself when the path has
length 1. If the original preceding node was P' (eliminated in
step 1), then B is P's other input, to which Q was rewired when
P' was removed. Let A be a descendant of Q other than G (exists
since Q is 2+). Let C be a descendant of G: G is not the output
gate, because the second substitution below trivializes G, and
a trivialized output gate would make f|R' constant,
contradicting Proposition 2.4; hence G has a descendant.
Distinctness. B != G (B feeds Q, G is fed by Q; equality would
be a cycle Q -> G -> Q containing the AND-type gate G,
impossible since cycles in a semicircuit are xor-only).
A != G by definition. C != G (C is a descendant of G).
C != B (C is fed by G, B feeds Q; equality would be a cycle
Q -> G -> B -> Q containing G, impossible). The remaining
possible coincidences are B = A and C = A. B = A is a genuine
xor-only coincidence (a cycle Q -> A -> Q of xor/equiv gates
is permitted in a fair cyclic xor-circuit), so it is not
excluded; it is handled by replacement-gate accounting in the
subcase below. C = A is handled separately below. Otherwise
G, A, B, C are four distinct gates.
Second substitution. Substitute a constant to Q via
xor-reconstruction with the now-unprotected x_j (Section 5:
the operation applies only for an unprotected variable), choosing
the constant b so that G is trivialized. The reconstruction
replaces Q by the constant b and adds the gate Z (x_j becomes a
gate), net gate count unchanged, potential unchanged (Section 5:
Z is xor-type, so no new troubled gate). After reconstruction,
b feeds G, A, and B (B is the preceding node on the reversed
path, now fed by b and its other input), so: G is trivialized,
removed by Rule 2; B is degenerate, removed by Rule 3; A is
fed by b and its other input, hence degenerate, removed by
Rule 3; C is fed by G and its other input, hence after G is
removed it is fed by a constant, removed by Rule 2 or 3.
Case split on potential. In this subcase assume the elimination
of B and of G does not increase Phi; the complementary cases
are 8.2.5.1 (G's elimination increases Phi) and 8.2.5.2 (B's
elimination increases Phi), proved there. Then Delta Phi <= 0
for G, <= 0 for B, <= 1 for A (Rule 3), <= 1 for C (Rule 2 or 3),
so Delta Phi <= 2 in total. Also x_j becomes non-influential:
the reconstruction makes x_j a linear variable in R', so it is
no longer free. Hence, in the distinct case,
Delta mu >= 4 - 2*alpha_phi + alpha_I
= 4 - 2/5 + 48/5 = 66/5 for this substitution.
Coincidence B = A (replacement-gate accounting). Suppose
B = A, and write M for this gate. M is an xor/equiv gate
feeding Q and fed by Q. Take the path from x_j to Q shortest,
as above, with nodes x_j = I_0 -> I_1 -> ... -> I_k = Q,
k >= 2 (if k = 1 then B = x_j is a variable, so B != A).
Then B = I_{k-1} = M, and Q feeds M only through M's side
input: M's inputs are I_{k-2} and T, and Q cannot be I_{k-2}
(the path nodes are distinct), so T = Q. Write U for Q's
other input, so Q's inputs are M and U. After
xor-reconstruction with Q <- b: M is fed by b and U, hence
degenerate, removed by Rule 3; and the gate on the reversed
path fed by M and Q — namely I_{k-2} if k >= 3, or the new
gate Z replacing x_j if k = 2 — is fed by M and b, hence
degenerate, removed by Rule 3. (If k >= 3,
I_{k-2} = g_{k-1}(M, Q) becomes g_{k-1}(M, b); if k = 2,
Z = g_1(M, Q) becomes g_1(M, b); in either case the gate
computes M or NOT M, up to a constant, since g_j is in
{oplus, equiv}.) G is trivialized and removed by Rule 2, and
C, fed by G, is removed by Rule 2 or 3. So four distinct
gates are removed: G, M, C, and the replacement gate
I_{k-2} (or Z). Distinctness: the replacement gate is an
xor/equiv gate in the cyclic part, distinct from M (path
nodes distinct, and Z is fed by M), from G (AND-type), and
from C (C is fed by G, hence downstream of G in the acyclic
part, while the replacement gate is in the cyclic part
feeding M -> Q -> G; no edge runs from the acyclic part back
into the cyclic part). In this subcase C != A automatically,
for the same reason: C is fed by G, while A = M is in the
cyclic part. Potential: G <= 0 by the case assumption; M <= 1
by Rule 3; the replacement gate <= 0 by the Section 6
refinement (it is removed before M, so its I_1 is the gate
M); C <= 1. Total Delta Phi <= 2, so Delta mu >=
4 - 2*alpha_phi + alpha_I = 66/5 for this substitution,
exactly as in the distinct case. The x_j credit is
unaffected: x_j is a linear variable after reconstruction
whether or not Z survives as a gate.
Shared-descendant case C = A. Let D be the gate fed by both Q
and G. D cannot be a 0-gate: a non-output 0-gate is removed by
Case 0.2, and if D were the output gate, then after the
substitution D is fed by the constants b and the trivialized G,
so f|R' would be constant, contradicting Proposition 2.4.
Hence D has a descendant E. The gates removed are G, B, D, E: B != D, since B is a node
of the cyclic xor-circuit while D is fed by the AND-type gate
G and so lies downstream of G in the acyclic part, and no edge
runs from the acyclic part back into the cyclic part.
E is distinct from G, B, D (E is fed by D; E = G would need
D = P or D = Q, both cycles; E = B would need a cycle
Q -> G -> D -> B -> Q). Their potential contributions: G <= 0
and B <= 0 by the case assumption; D is removed by Rule 3 with
I_1 = G, a gate, so by the Section 6 refinement no rewired
descendant can become troubled and Delta Phi = 0 for D; E is
removed by Rule 2 or 3 after G's removal, Delta Phi <= 1.
Total Delta Phi <= 1, four distinct gates, x_j killed, so
Delta mu >= 4 - alpha_phi + alpha_I = 4 - 1/5 + 48/5 = 67/5.
In each of the three cases (distinct, B = A, shared) the
second substitution gives at least 66/5.
Averaged with the previous substitution,
[(1 + alpha_I + alpha_Q) + (4 - 2*alpha_phi + alpha_I)]/2
= (5 - 2*alpha_phi + 2*alpha_I + alpha_Q)/2
= (5 - 2/5 + 96/5 + 9/5)/2 = 64/5 = delta per substitution
(using the smaller 66/5 bound; the shared case gives more).
8.2.2.2 Q is fed by x_k: after the substitution, the
other input Q' of Q is passed to feed G; Q' is also a
2+-gate depending on the unprotected x_j, so apply
8.2.2.1 with Q' in place of Q. Henceforth Q is a 1-gate,
so P is a variable or a 1-gate.
8.2.3 P = t a variable, unprotected 1-variable by Cases
1, 3 and 7: let x_j in I be protected with couple x_k.
Substitute a constant to x_k making x_j unprotected;
after normalization G is fed by t and a gate Q' (possibly
Q) depending on x_j. Substituting a constant to Q' via
xor-reconstruction trivializes G and makes t, x_j, x_k
non-influential: three variables in two substitutions,
so Delta mu >= 3*alpha_I/2 >= delta per substitution.
8.2.4 P an XOR-type 1-gate: if P depends on an
unprotected variable, apply 8.1.3 with P and Q exchanged,
so assume both P and Q depend only on protected variables.
8.2.4.1 P or Q depends on x_j with couple x_k not feeding
P or Q: substitute a constant to x_k and normalize; by
Case 6 exactly one gate is eliminated,
Delta mu = 1 + alpha_I + alpha_Q, and P (or Q) is still
an XOR-type 1-gate, now depending on the unprotected
x_j. Substituting a constant to Q (or P) via
xor-reconstruction trivializes G, giving
Delta mu >= 4 - 4*alpha_phi + alpha_I for that
substitution, and >= 4 - 2*alpha_phi + alpha_I if
Delta Phi <= 2 for it; in the latter case the average
with the first substitution is
(5 - 2*alpha_phi + 2*alpha_I + alpha_Q)/2 = delta,
as in 8.2.2.1. If Delta Phi >= 3 for that substitution,
then both inputs of P must increase potential by one.
The possibilities are: (a) P is fed by a gate T just
before its removal: since there is a path from T to G
after xor-reconstruction, T must be XOR-type and cannot
become troubled, impossible; (b) P is fed by two
variables even before the x_k substitution: both are
protected 1-variables by Case 1, which can neither feed
nor become troubled gates, impossible; (c) P is fed by
two gates, with variables t_1, t_2 passed to it by the
degeneracy of B (the gate on the reversed path from G
to x_j) and of the descendant of x_k; (d) P is fed by a
variable t_1 and an XOR-type gate passing a variable
t_2 to G via the degeneracy of B and of the descendant
of x_k. Cases (c) and (d) are handled in 8.2.4.1.1-3:
8.2.4.1.1 If t_1 or t_2 is protected, after normalization
it feeds a troubled gate or is an inner variable of a
pack; substituting a constant to it eliminates at least
three gates by Rules 2 or 3 while killing a quadratic
equation, giving Delta mu >= 3 - 3*alpha_phi + alpha_Q
+ alpha_I = 3 - 3/5 + 9/5 + 48/5 = 69/5. Together with
the two previous substitutions (62/5 and at least
64/5), the average over three substitutions is at
least delta (the sum is at least 62/5 + 64/5 + 69/5
= 195/5, divided by 3 is 13 >= delta).
8.2.4.1.2 If t_1 and t_2 are both passed to P during
degeneracy, with t_2 passed only via B and t_1 only via
the descendant of x_k: B must be a gate (if B = x_j it
would be a protected 2-variable, impossible by Case 1),
t_2 is unprotected, Q does not depend on it, and
Q = B xor t_2 xor c_1 while P = x_k xor t_1 xor B xor
c_2, so P depends on t_2, contradicting the assumption
that P depends on no unprotected variable.
8.2.4.1.3 If t_1 directly feeds P before the x_k
substitution: if t_2 directly feeds P or is passed only
via the descendant of x_k, then P depends on the
unprotected t_1, impossible by the Case 8.1 assumption.
Otherwise t_2 is passed by B (alone or with the
descendant of x_k); since B feeds P or Q directly in
the original circuit, either P depends on t_1, or B
depends on t_1 and hence Q depends on t_1, again
impossible. After this case, P depends only on
variables whose couples feed Q, and Q only on variables
whose couples feed P.
8.2.4.2 P is fed by two protected variables x_j, x_k,
not a couple (by 8.2.1): both are 1-variables by
Case 1. Substitute constants to both, trivializing P
and G and eliminating Q and at least one descendant of
G: four gates removed, Delta Phi <= 3 (P and G give 0,
Q gives at most 2, the descendant at most 1), two
quadratic equations and two influential inputs killed.
Total Delta mu >= 4 - 3*alpha_phi + 2*alpha_I + 2*alpha_Q
= 4 - 3/5 + 96/5 + 18/5 = 131/5; per substitution (two)
this is 131/10 >= delta.
8.2.4.3 P and Q are each fed by at most one variable:
by Case 0.3 neither computes a constant, so P depends
on at least one protected variable x_j whose couple
x_k feeds Q. The other input of Q must be a gate R
(by 8.2.4.2, it cannot be a second variable), R does
not compute a constant by Case 0.3, and R does not
depend on x_k (there is no path from x_k to R, since
x_k is a 1-variable by Case 1 and Q is a 1-gate feeding
an AND-type gate). Let x_i be a variable depended on
by R, also depended on by Q; x_i != x_j (no path from
x_j to R). Then Q depends on two protected variables
x_k and x_i, not a couple of each other, so P must be
fed by two protected variables, impossible by 8.2.4.2.
Hence this subcase cannot occur.
8.2.5 The deferred Delta Phi > 2 cases from 8.2.2 and
8.2.4. Substitute x_k <- d and normalize as in Case 6
(one gate, no new troubled gate), then substitute a
constant to Q via xor-reconstruction with x_j now
unprotected, eliminating G, B (the preceding node of Q
on the path from x_j to Q) and the descendant of G.
Consider a potential increase from eliminating G or B.
8.2.5.1 Elimination of G increases potential by 1: then
P must be a variable x_m whose out-degree decreased
from 3 to 2, so x_m is a 3-variable after the x_k
substitution.
8.2.5.1.1 x_m a 3-variable even before x_k <- d: by
Case 3 it cannot feed G before, so it is passed to G
via the degeneracy at x_k <- d, through a degenerate
gate P'. To increase potential, x_m must feed an
AND-type gate T just before G is eliminated; by the
xor-reconstruction picture it feeds T even before
xor-reconstruction, just after x_k <- d. Either it
feeds T before x_k <- d, impossible by Case 3 (a
3-variable feeding an AND-type gate), or it is passed
via P', in which case P' must be a 2+-gate, impossible
since to make T troubled while eliminating G, x_m
must be a 2-variable after elimination, so P' can only
be a 1-gate.
8.2.5.1.2 x_m a 1-variable before x_k <- d: then
x_m and x_k feed a common XOR-type gate P', and the
same argument as in Cases 6.2.2.1-6.2.2.3 (affine
substitutions through a common XOR gate) gives a block
killing at least three variables in at most three
substitutions, at least delta per substitution on
average.
8.2.5.1.3 x_m a 2-variable originally, out-degree
increased to 3 after x_k <- d: then P' must be an
XOR-type 2-gate fed by x_m and x_k, and x_m either
directly feeds G originally or is passed via P'.
8.2.5.1.3.1 x_m initially feeds G: P' must feed an
AND-type gate T even before x_k <- d to form a troubled
gate. Substitute constants to x_m and x_k trivializing
P', T and G, degenerating their descendants; six gates
are removed by Rules 2 and 3, with total
Delta mu >= 6 - 6*alpha_phi + alpha_Q + 2*alpha_I
= 129/5, i.e. 129/10 >= delta per substitution (two).
8.2.5.1.3.2 Q is a 2-gate: P' is a 2+-gate depending on
the unprotected x_m, so apply Case 8.1 with P' and Q
exchanged, giving at least delta.
8.2.5.1.3.3 Q is a 3+-gate: first substitute x_k <- d
and normalize (Case 6, Delta mu = 1 + alpha_I +
alpha_Q), so P' is eliminated and Q becomes an XOR-type
gate depending on the unprotected x_j. Then substitute
a constant to Q via xor-reconstruction trivializing G,
eliminating five gates by Rules 2 and 3 (three
descendants of Q, the first gate on the reversed path
from Q to x_j, and the descendant of G), giving
Delta mu >= 5 - 5*alpha_phi + alpha_I = 68/5 for that
substitution; averaged with the first, (62/5 + 68/5)/2
= 65/5 >= delta.
8.2.5.2 Elimination of B increases potential by 1: B
feeds Q even before x_k <- d. After xor-reconstruction,
the other input of B apart from the constant b must be
a node D that becomes or feeds a troubled gate; by the
xor-reconstruction picture D feeds Q before
xor-reconstruction, either directly or via the
descendant of x_k, and is passed to B at the
normalization after x_k <- d. By the normalization
lemma, D either becomes or feeds a troubled gate T
after normalization; it cannot be an AND-type gate by
minimality of G, so D is a variable, say x_l.
8.2.5.2.1 x_l unprotected: Q cannot depend on x_l after
Case 8.1, and Q = B xor x_l xor c (or with x_k xor c),
so B must depend on x_l.
8.2.5.2.1.1 B is a 0-gate just before elimination, so
its out-degree change decreases out(x_l): then B = Z
(the new gate from xor-reconstruction), otherwise it
feeds the next node on the path, so B = x_j just before
xor-reconstruction, Q is fed by x_j and x_l and depends
on them, contradicting B depending on x_l.
Henceforth B is a 1+-gate before elimination, so
out(x_l) does not decrease, and x_l must be passed to
the newly introduced troubled gate T during elimination
of B (if x_l fed T originally, T would already be
troubled). The other input of T apart from B is a
variable t.
8.2.5.2.1.2 Both B and t feed T even before x_k <- d:
T is a topologically minimal AND-type gate fed by an
XOR-type gate B depending on the unprotected x_l,
impossible after Case 8.1.
8.2.5.2.1.3 B is passed to T when the descendant of
x_k degenerates, t feeds T originally: since x_k can
only be XOR-type in that configuration, T is a
topologically minimal AND-type gate fed by an XOR-type
gate (the descendant of x_k) depending on the
unprotected x_l, impossible.
8.2.5.2.1.4 t is passed to T when the descendant of
x_k degenerates, B feeds T originally, and t is a
1-variable: the descendant of x_k is an XOR-type gate
fed by two 1-variables, and affine substitutions as in
Cases 6.2.2.1-6.2.2.3 give at least delta per
substitution on average over the alternative block's own
substitutions.
8.2.5.2.1.5 Same, with t a 2-variable: the descendant
of x_k must be a 1-gate to make T troubled, and B
depends on the unprotected x_l, impossible after
Case 8.1.
8.2.5.2.2 x_l protected: after normalization x_l either
feeds a truly troubled gate or is an inner variable of
a pack, so in either case x_l is a protected 2-variable
feeding an AND-type gate. Substitute a constant to x_l
to remove three gates. Total over three substitutions:
x_k <- d gives >= 1 + alpha_I + alpha_Q = 62/5;
xor-reconstruction removes B, G and the descendant of
G with Delta Phi <= 2, giving >= 3 - 2*alpha_phi +
alpha_I = 61/5; x_l <- c removes three gates with
Delta Phi <= 1 (one originally troubled gate eliminated
first gives 0, the other two give at most 2, and
destruction of the troubled gate or pack gives -1),
giving >= 3 - alpha_phi + alpha_I + alpha_Q = 71/5.
Sum = 194/5, average 194/15 >= delta.
8.2.5.3 B feeds B' and B' feeds Q originally, with B'
the gate eliminated in the x_k substitution. Configuration:
B' is the Case 6 gate, an XOR-type gate fed by B and x_k;
after step 1 (x_k <- d) B' is removed by Rule 3 and Q is
rewired to B, so after step 1 B feeds Q and B is the
preceding node of Q on the x_j-to-Q path. Original-state
facts: (i) x_k does not feed Q: otherwise Q, fed by B' and
x_k, becomes degenerate when B' is removed and x_k <- d,
so step 1 eliminates two gates, which is Case 6.1 with
Delta mu >= 65/5 >= delta for step 1 alone, and no deferred
analysis is needed. (ii) The path from x_j to Q after step
1 passes through B, and the xor-reconstruction in step 2
applies to this path exactly as in 8.2.5.2. (iii) The
equation is Q = B xor x_k xor x_l xor c for a constant c
(LY Case 8.2.5.3, fact 2; read from Figure 46), with x_k
constant after step 1, so every dependence or
out-degree claim in 8.2.5.2 transfers with the x_k term
absorbed into c. Fact 1: the potential increase from B's
elimination is fed by a variable x_l just after
xor-reconstruction, and x_l feeds Q in the original
circuit. Proof: by the normalization lemma, a newly
troubled gate or pack created by B's elimination is
created from B's other input D after reconstruction: D is
topologically before G (D feeds B, B feeds Q, Q feeds G),
so if D were a gate it would be XOR-type in the cyclic
part and could neither become a troubled gate (troubled
gates are AND-type) nor feed one (troubled gates are fed
only by variables); hence D is a variable x_l. And Q
depends on x_l in the original circuit: x_l feeds B, B
feeds B', B' feeds Q, so x_l is in Q's dependency cone.
Fact 2: if x_l is unprotected, then Q does not depend on
x_l (Case 8.2 hypothesis: Q depends only on protected
variables), while Q = B xor x_k xor x_l xor c forces B to
depend on x_l (otherwise Q would depend on x_l through
the explicit term, since x_k is constant after step 1).
Also B is not a 0-gate just before elimination: if it
were, then as in 8.2.5.2.1.1 B = Z, so B = x_j just before
xor-reconstruction, and Q would be fed by x_j and x_l
and depend on x_l, contradicting Fact 2. Fact 3: if x_l
is protected, substituting a constant to x_l after the
whole process gives sufficient decrement; see the
protected accounting below. Unprotected case: the
8.2.5.2.1 subcases transfer as follows. 2.1.1: B a 0-gate
just before elimination is excluded by Fact 2 as above.
2.1.2: a minimal AND gate T fed by B and t before step 1
is fed by an XOR-type gate B depending on the unprotected
variable x_l, so Case 8.1 would have applied to T, a
contradiction; the same holds if T is fed by B' instead
of B, since B' = B xor x_k xor c depends on x_l through
B. 2.1.3: if t is passed to T when B' degenerates, then
T was fed by B' before step 1, and T is a topologically
minimal AND-type gate fed by the XOR-type gate B'
depending on the unprotected variable x_l, impossible
by Case 8.1. 2.1.4: Identifying the configuration by the
counterfactual step-1 analysis: if t is passed to T when B'
degenerates and t is a 1-variable, then t was an input
of B'; B' is fed by B and x_k, and t cannot be x_k
(x_k is substituted away in step 1 and cannot be passed
to T), so t = B, and B is a variable, hence B = x_j,
the only variable on the path. x_j is a 1-variable by
Cases 1-2 and x_k is a 1-variable by Case 6, so B' is
an XOR-type gate fed by two 1-variables and the affine
substitutions of Cases 6.2.2.1-6.2.2.3 apply from the
original state, in place of the deferred sequence (steps 1-3
are not performed in this subcase), giving at least delta
per substitution within that block. 2.1.5: if t is a 2-variable, then t was an input
of B', so t = B (t cannot be x_k, substituted away in
step 1); B is a variable only if B = x_j, a 1-variable
by Cases 1-2, so t cannot be a 2-variable; and in any
case B depends on the unprotected variable x_l while
feeding the minimal AND gate T, contradicting Case 8.1. Protected case: x_l
is protected, and the third substitution x_l <- c gives
the same accounting as 8.2.5.2.2, with the distinctness
proved as follows. By Fact 1 the potential increase from
B's elimination is fed by x_l, so after step 2 x_l
either feeds a truly troubled gate T or is an inner
variable of a pack; in either case x_l is a protected
2-variable feeding an AND-type gate, exactly as in
8.2.5.2.2. The three gates removed by x_l <- c are: (a)
the troubled gate T fed by x_l (or a gate of the pack
containing x_l), removed first; (b) the other gate U
fed by x_l (x_l is a 2-variable after step 2), removed
by Rule 2 or 3 when fed by the constant c; (c) a third
gate W, a descendant of T if distinct from U, or a
descendant of U if T's descendant is U itself. The
three are distinct: T and U are distinct gates both fed
by x_l; W is a descendant of T or of U, hence distinct
from both (a descendant is never the gate itself, and
if W is a descendant of T equal to U, then W is a
descendant of U instead). W exists: T is not the output
gate (if T were the output, x_l <- 0 would make the
output constant, contradicting Proposition 2.4, since
only three substitutions have been made, well within
the duration), so T has a descendant; and if that
descendant is U, then U is not the output either (if U
were the output and fed by T and x_l, then x_l <- 0
with T = 0 would make the output constant, again
contradicting Proposition 2.4), so U has a descendant.
Also none of T, U, W coincides with B, G, or G's
descendant removed in step 2 (those are removed
already, and the third substitution counts only gates
present after step 2). The potential bound is as in
8.2.5.2.2: T's removal gives Delta Phi <= 0 (both its
inputs are 2-variables, which cannot feed troubled
gates after normalization), the other two give at most
2, and destroying the troubled gate or pack decreases
Phi by 1, so Delta Phi <= 1. Hence the three
substitutions give Delta mu >= 62/5 (step 1), >= 61/5
(step 2, Delta Phi <= 2), and >= 3 - alpha_phi +
alpha_I + alpha_Q = 71/5 (step 3, three gates, x_l and
its quadratic equation killed), for a total of 194/5
and an average of 194/15 >= delta, exactly as in
8.2.5.2.2. The potential bound Delta Phi <= 2 in step 2
uses the same case split as 8.2.5.2.2: G's elimination
does not increase Phi, and B's elimination increases
Phi by at most 2 (a newly troubled gate or pack is fed
by x_l, and x_l's other descendants, if any, are
accounted for in the third substitution above, not
double-counted here).
8.2.5.4 B feeds Q' and Q' feeds Q, with Q eliminated in
the x_k substitution and Q' taking Q's place.
Configuration (the 8.2.2.2 configuration; LY Figure 47):
Q is fed by x_k and Q'; Q' is a 2+-gate depending on x_j,
fed by B (the preceding node of Q' on the path from x_j
to Q') and one other node R. Step 1 (x_k <- d) makes Q
degenerate (fed by the constant d and Q'), so Q is removed
by Rule 3 and G is rewired to Q'; this is the Case 6
building block (exactly one gate eliminated,
Delta mu = 1 + alpha_I + alpha_Q = 62/5). If D = G (i.e.
Q' already feeds G, so P = Q'): then step 1 eliminates Q
and G, not Q alone. Q is removed by Rule 3 as above, and
G, rewired to Q' twice, has coincident inputs and is
removed by Rule 5 (which removes exactly such gates, with
no input fixed and no substitution, as in Section 6). Two
gates are removed, Delta Phi <= 1 + 1 by Rules 3 and 5,
and x_k and its quadratic equation are killed, so
Delta mu >= 2 - 2*alpha_phi + alpha_I + alpha_Q = 65/5
>= delta for step 1 alone (the same bound as Case 6.1),
and no deferred analysis is needed in this subcase.
Henceforth assume D != G. After step 1,
Q' feeds G and is still a 2+-gate: Q' fed Q and at least
one other gate D before step 1, and after Q's removal Q'
feeds G in Q's place plus D, so its out-degree is
unchanged. B still feeds Q', and B is the preceding node
of Q' on the x_j-to-Q' path (Q' still depends on x_j by
8.2.2.2, so the path exists by Proposition 2.7).
Transfer: the analysis of 8.2.5.2 applies with Q' in
place of Q, as Li-Yang state ('Similar to the previous
case, one may carefully check that the analysis in
Case 8.2.5.2 still works'), not the analysis of 8.2.5.3:
the gate eliminated in step 1 here is Q, not Q', and Q'
survives to take Q's place. Explicitly, Fact 1 holds with
Q' for Q (the potential increase from B's elimination is
fed by a variable x_l just after xor-reconstruction, and
x_l feeds Q' in the original circuit (D feeds Q' before
reconstruction, directly or via Q, the descendant of x_k;
via Q is impossible, since Q's inputs are Q' and x_k, and
x_l is live after step 1, while x_k is substituted away,
so x_l != x_k and x_l does not feed Q, hence directly);
Fact 2 holds with the
equation Q' = B xor x_l xor c, with no x_k term, since
x_l reaches Q' directly and not through Q, the descendant
of x_k. The subcases 2.1.1-2.1.3 transfer with Q' for Q,
giving contradictions as in 8.2.5.2 (2.1.3, if its
configuration arises); 2.1.4 and 2.1.5 are vacuous here:
2.1.4 requires the step-1 gate to be an XOR-type gate fed
by two 1-variables, but Q is fed by Q' (a gate) and x_k,
and 2.1.5 requires the descendant of x_k to be a 1-gate,
but Q is a 2+-gate (Case 8.2.2). The protected case gives
the same three-substitution accounting
(62/5 + 61/5 + 71/5, average 194/15 >= delta), since step 2
removes B, G, and the descendant of G with
Delta Phi <= 2, exactly as in 8.2.5.2.2. The 8.2.5.1 branch
cannot arise in this configuration at all, so its
one-line transfer in the previous wording is withdrawn.
Its premise is that G's elimination increases Phi by 1,
so P must be a variable x_m with out-degree 3 after
step 1; but P fed G before step 1, step 1 removes only Q
(fed by x_k and Q', under D != G as assumed above), so no
live variable is passed to any gate and no live variable's
out-degree changes, and Case 3
caps a variable feeding an AND-type gate at out-degree 2,
a contradiction. In particular the P' subcases are vacuous:
P' is a gate fed by x_m and x_k, but the only gate
eliminated in step 1 is Q, fed by Q' and x_k, so x_m
would have to feed Q, hence x_m = x_k, impossible for a
live variable after step 1; indeed P' cannot exist at all,
since x_k is a 1-variable feeding only Q. The 1.1 former
branch (x_m feeding T before step 1) is impossible by
Case 3, as in 8.2.5.1.1, and the 1.2 affine block and the
1.3.1 block are standalone alternatives giving at least
delta per substitution whenever their configurations
arise, though they do not arise here, since P' cannot
exist. Fan-outs after step 1
are: Q' feeds G and D (out-degree unchanged, still 2+);
B feeds Q'; G is fed by P and Q'. No gate is credited
twice: Q is eliminated in step 1 and counted there
(62/5), and step 2 counts only B, G, and the descendant
of G, all present after step 1.
This completes Case 8, and with it the case table:
every configuration of a normalized circuit with empty
packing and no constant or degenerate gate falls into
Case 0, 1, 2, 3, 4, 5, 6, 7 or 8, and each gives either
a simplification with Delta mu >= 0 or a block of
1 to 3 substitutions with Delta mu >= delta per
substitution, at the published point. Hence the
Elimination Lemma holds as stated.

From the Elimination Lemma, the measure lower bound follows
without further case work, and this derivation IS proved here.
Start with C_0 computing f, R_0 empty, P_0 any maximal packing,
so mu_0 = mu(C_0, P_0, empty) = |C_0| + alpha_I*i_0 +
alpha_phi*Phi(C_0,P_0) with q = 0 and i_0 <= n (not i_0 = n:
a variable of out-degree 0 that is unprotected is not
influential). Apply the lemma repeatedly, with the duration
invariant of Section 5: after k substitutions R has dimension
n - k, C computes f|R, and the lemma is applied only while
k < n - 2d - 5, each block adding at most 3 substitutions, so
the total is between n - 2d - 5 and n - 2d - 3 and never
exceeds n - 2d - 3, so f|R stays nonconstant by Proposition
2.4 with the terminal-case margin of Section 5. After
s >= n - 2d - 5 substitutions,
mu_s >= 0, and each substitution reduced mu by at least delta,
so mu_0 >= mu_s + delta*s >= delta*(n - 2d - 5). Combined
with mu_0 <= |C| + alpha_I*n + alpha_phi*(n/2 + 5d/2)
(i_0 <= n, Phi_0 <= n/2 + 5d/2 by Section 6, q_0 = 0), this
gives the inequality used in Section 8:

    delta*(n - 2d - 5) <= mu_0
        <= |C| + alpha_I*n + alpha_phi*(n/2 + 5d/2).

No step here claims mu can go below 0, and no single-gate or
constant-circuit boundary step is hidden in the -5: if the
circuit were trivialized before s substitutions, that would
contradict Proposition 2.4, since the total is at most n - 2d.

