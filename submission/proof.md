# Complete proof of the joint theorem: C_B2(f_n) >= (7/2)n - (539/25)n^(4/5) - 7

This file is the master proof for the joint theorem declared in
`bound.json`. The hard function is the family f_n of `field.md`, an affine
disperser for dimension d(n) = floor((77/25)*n^(4/5)) + 1 at every
n >= N0 := 10^20 (`disperser.md`, Theorem A), computable in polynomial time
(`field.md`, Proposition 3.1), hence a single total language in P with an
empty-witness NP verifier. The lower-bound argument is the gate-elimination proof of
Mishra-Santhanam (ECCC TR26-260, 19 October 2026, "An Explicit 3.5n - o(n)
Size Lower Bound against General Boolean Circuits"), reproduced below in full
with all cases and all potential computations; the theorem we prove is their
Theorem 4.1, and the final instantiation with the disperser above is new
text. No external premise is used beyond `field.md` and `disperser.md`.

## 0. The model and the statement

**Model.** A B_2-circuit is a directed acyclic graph in which every vertex
has indegree zero or two, the two incoming edges ("pins") of a gate are
ordered, indegree-zero vertices are inputs labelled X_1, ..., X_n (or
constants 0/1), and each gate is labelled by one of the sixteen binary
Boolean operations; one gate of outdegree zero is the output. Fan-out is
unrestricted and subcircuits may be shared. The size |C| is the number of
gates, each gate counted once, including trivial (constant) or degenerate
(unary) gates; inputs, wires, and the output designation cost nothing. This
is the Li-Yang general B_2 model (ECCC TR21-023, Sections 1.1, 2.1, 2.4
footnotes 3-4, 3.3), and it matches the conventions of TR26-260 Section 2.1:
negations on wires are absorbed for free into the operations of the gates
reading them; a retained gate that computes a negated projection is a unary
gate and is charged like any gate while present. Identify {0,1} with F_2.
Every operation h: F_2^2 -> F_2 has the unique multilinear form
h(t_1, t_2) = c + b_1 t_1 + b_2 t_2 + b_12 t_1 t_2; operations depending on
both inputs are XOR-type (t_1 + t_2 + c) or AND-type
((t_1 + al)(t_2 + be) + ga). An *affine* operation is b_1 t_1 + b_2 t_2 + c
with possibly vanishing coefficients.

**Main Theorem.** Let f: F_2^n -> F_2 be an affine disperser for dimension d
(nonconstant on every affine subspace of dimension >= d), and let C be a
B_2-circuit computing f. Then |C| >= (7/2)n - 7d.

**Corollary (the joint theorem).** For the family f_n of `field.md`, for
every integer n >= 10^20: C_B2(f_n) >= (7/2)n - (539/25)*n^(4/5) - 7;
consequently liminf_{n->inf} C_B2(f_n)/n >= 7/2.

*Proof of the corollary from the Main Theorem.* By `disperser.md` Theorem A,
f_n is an affine disperser for dimension
d(n) = floor((77/25)n^(4/5)) + 1 <= (77/25)n^(4/5) + 1 for every
n >= 10^20, so the Main Theorem gives, for every such n and every B_2-circuit
C for f_n,

    |C| >= (7/2)n - 7*((77/25)n^(4/5) + 1)
         =  (7/2)n - (539/25)*n^(4/5) - 7

using 7*77/25 = 539/25 exactly. This is the declared eventual lower bound
with the two explicit remainder monomials (539/25)*n^(4/5) and 7 at
N = 10^20, and dividing by n and letting n -> inf gives the liminf claim
since (539/25)n^(-1/5) -> 0 and 7/n -> 0. The
language L = {x : f_{|x|}(x) = 1} is in P by `field.md` Proposition 3.1, with
the fixed verifier V(x, w) := A(|x|, x) and witness bound 0. Both directions
of NP membership are immediate; witnesses never enter the circuits. qed

The rest of this file proves the Main Theorem.

## 1. States: parameters, the map Phi, fair semi-circuits

**Parameters and Phi.** We maintain an injective map
Phi: F_2^m -> F_2^n of the form

    Phi(u) = b + L*u + sum_{ {p,r} in P } c_{p,r} * p*r,                (1)

where b in F_2^n, L an n x m matrix over F_2, P a set of disjoint unordered
pairs of parameters, and c_{p,r} in F_2^n. Parameters in a pair are *paired*
(the other member is the *partner*); all others are *unpaired*. A pair
contributes one degree of freedom (the term p*r becomes affine once either
member is fixed or p is identified with a literal of r). Initially Phi =
identity (m = n, b = 0, L = I, P empty).

**Fair semi-circuits.** The current *circuit* is a fair semi-circuit reading
the parameters u = (u_1, ..., u_m):

(i) an *affine system* A v = N u + c of r_eq >= 0 equations over F_2 in the
    parameters and in auxiliary variables v = (v_1, ..., v_r_eq), where each
    equation contains at most three nonconstant terms (each row of (A | N)
    has at most three nonzero entries), and the auxiliary matrix A is
    invertible over F_2 (*fairness*: every parameter assignment extends
    uniquely to a solution, v = A^{-1} N u + A^{-1} c);
(ii) an *acyclic part*: an ordinary acyclic circuit whose gates, in a fixed
    topological order, are labelled by B_2 operations and read constants,
    parameters, auxiliary variables, or earlier gates; no wire goes back
    from the acyclic part to the affine system;
(iii) an *output*: a constant, or a literal of a parameter, an auxiliary
    variable, or a gate.

The *size* s of a fair semi-circuit is (number of equations) + (number of
gates in the acyclic part). Initially the affine system is empty and the
acyclic part is the given circuit C, so s = |C|; gates and equations are
charged one each (constants are free, as in the pinned conventions).

**Main invariant.** C(u) = f(Phi(u)) for every u in F_2^m. The procedure
below maintains it; in particular, whenever f∘Phi is nonconstant, the output
of the circuit cannot be a constant, since the computed function is exactly
f∘Phi. Since the original C computes f and Phi_0 = id, the invariant holds
initially.

**Moving affine gates into the system.** If a gate J of the acyclic part has
an affine operation h(t_1,t_2) = b_1 t_1 + b_2 t_2 + c and both inputs are
constants, parameters, or auxiliary variables, replace J by the equation
J + b_1 t_1 + b_2 t_2 = c and regard J as a new auxiliary variable. The new
equation has at most three nonconstant terms; the new auxiliary matrix is
block-triangular with the old A and a 1, hence invertible; the function
computed is unchanged; and s is unchanged (one gate replaced by one
equation). We use this freely as a preparatory step.

## 2. Parameter changes

Four kinds of change to the parameters, each replacing Phi by Phi∘tau; in
all four cases tau is injective F_2-affine(-plus-quadratic) of the same form
(1), so Phi∘tau is injective of the form (1) (a fresh parameter w is never
paired; quadratic terms created by substitution involve only unpaired
parameters, and pairs stay disjoint):

- **AR** (affine reparameterisation): for an unpaired parameter x and affine
  phi, set x <- w + phi(u_{-x}) with w fresh. tau is a bijection F_2^m ->
  F_2^m; m, q, a = m - q, and im Phi are unchanged.
- **QR** (quadratic reparameterisation): for distinct unpaired x, y, z and
  constants al, be, ga, set z <- w + (x+al)(y+be) + ga, w fresh, and add the
  pair {x, y} to P. m unchanged, q increases by 1, a decreases by 1; tau is
  a bijection.
- **UF** (fixing an unpaired parameter): for an unpaired parameter x and an
  affine map psi, impose x = psi(u_{-x}) (constants allowed). m decreases by
  1, q unchanged, a decreases by 1; im Phi shrinks (or stays).
- **PF** (fixing a paired parameter): for a pair {p, r} in P, fix p to a
  constant or impose p = r + c, then remove p and the pair {p, r}. m and q
  each decrease by 1, a is unchanged.

The algorithm applies UF only with a constant or x = y + c, and PF only with
a constant or p = r + c; Lemma 2.3 below additionally uses a general affine
UF, purely as a counterfactual. Fixing a parameter to a constant, or to a
literal of another parameter, deletes a column of N (updating c), or
replaces x-occurrences by (y + c)-occurrences; neither changes A, and the
number of nonconstant terms per equation does not increase. (AR) and (QR)
are implemented in Section 3.

**Proposition 2.1 (images contain affine subspaces).** The image of every
injective Phi of the form (1) contains an affine subspace of dimension
a = m - q. Hence, if f is an affine disperser for dimension d and a >= d,
then f∘Phi is nonconstant.

*Proof.* Fix one parameter of each pair in P to 0: q parameters become
constants, every quadratic term becomes affine in the remaining m - q
parameters, and the restricted map is injective affine F_2^(m-q) -> F_2^n,
so its image is an affine subspace of dimension m - q = a. If a >= d and
f∘Phi were constant, then f would be constant on this a-dimensional affine
subspace, hence on some d-dimensional affine subspace of it, contradicting
dispersal. qed

**Lemma 2.2 (weights).** Let wt(x) = 0, 6, 7 for an unpaired parameter x of
degree 0, 1, >= 2 respectively (degree = number of distinct children, where
the children of a value are the gates reading it and the equations
containing it with nonzero coefficient; a gate reading the same value on
both pins counts once). Then wt(x) <= 7 always, and W := 2s + 3q +
sum_{x unpaired} wt(x) satisfies W >= 0.

*Proof.* All terms are nonnegative. (Stated for reference; degree and
children are syntactic notions.) qed

**Lemma 2.3 (the output is not simple).** Let a > d. Then the output of the
current circuit is not a constant, not a nonconstant affine function of the
parameters, and not a literal of an AND-type gate whose two parents are
unpaired parameters.

*Proof.* A constant output contradicts Proposition 2.1 (a >= d implies
f∘Phi nonconstant, but the circuit computes f∘Phi). Let the output be
affine in the parameters and nonconstant, with coefficient 1 on some
parameter x. If x is paired, fix its partner to a constant (PF): a is
unchanged and x becomes unpaired; a > d still holds. Now substitute for x
the affine expression making the output constant (UF with
x <- (sum of the other affine terms) + (constant offset)); a decreases by 1,
so a - 1 >= d and the restricted image still contains an affine subspace of
dimension >= d (Proposition 2.1), on which f — and hence the computed output
— is nonconstant; but the output is now constant. Contradiction. Finally,
let the output be a literal of an AND-type gate (x+al)(y+be)+ga on unpaired
parameters x, y: fix y := be (UF); the gate, hence the output, becomes
constant, while a - 1 >= d, the same contradiction. qed

## 3. Reparameterisation in the circuit, rerouting, and deletions

**Rerouting.** When a gate or auxiliary variable becomes a constant or a
literal of another internal value, we replace all its occurrences by that
constant or literal (if it is the output, the output designation moves to
the constant or literal). Two facts are used later, each argued where it is
needed (Section 4.2 and Case 8):

(R1) For parameter weights: rerouting to a constant cannot increase any
     parameter's degree; when a gate is D1-deleted and rerouted to a
     literal of one of its inputs, only that input (and, for a D3-deleted
     equation, only the parameter or auxiliary variable it solves for) can
     gain children — so a parameter gains children only if it is such an
     input or term of a deleted value.

(R2) Readers of a deleted value V are rewired to constants or to literals
     of V's input values — so a gate can acquire new readers only while one
     of its own children (a gate reading it) is being deleted and rerouted
     to a literal of it. No blanket "children only shrink" claim is made;
     instead, Section 4.2 and Case 8 each prove the termination/accounting
     they need from the local order of operations.

**AR in the circuit (exchanging a parameter and an auxiliary variable).**
Let Q be an auxiliary variable whose solved expression in terms of the
parameters is Q = x + phi(u_{-x}) for an unpaired parameter x and affine phi
(equivalently, the Q-entry of the column A^{-1}N_x of x equals 1, where N_x
is the column of x in N). Then x occurs in some equation (otherwise its
column in N, hence in A^{-1}N, would be zero). Exchange the roles: substitute x <- w + phi(u_{-x})
in Phi, let Q act as the new parameter w and let x act as a new auxiliary
variable. Concretely, exchange the column of x in N with the column of Q in
A. Every equation keeps exactly the same terms (the relabeling swaps the two
symbols), and the new auxiliary matrix A' = A·B where B is the identity with
its Q-column replaced by t := A^{-1} N_x (N_x the column of x in N); since
the coefficient of x in Q's expression is t_Q = 1, det B = 1 and A' is
invertible. No gate or equation changes, so degrees of all parameters other
than x are unchanged, m, q, a, s are unchanged, and Delta W = wt'(w) - wt(x),
where the new parameter w has exactly the children of Q (deg'(w) = deg(Q) >=
1, since the column of Q in the invertible A is nonzero, so w occurs in an
equation). The algorithm applies AR only when either x is isolated (inside
D4, where the possible +1 is offset by an equation deletion) or deg(x) >= 2
and Q is an input of the minimal AND (so deg'(w) >= 2 and Delta W = 0).

**QR in the circuit.** QR is applied only in the bottleneck configuration:
distinct unpaired parameters x, y, z; an AND-type gate G = (x+al)(y+be)+ga
whose only child is an XOR-type gate D = G + z + c; z isolated with only
child D. After the substitution z <- w + G(x,y) + c, gate D computes
G + (w + G + c) + c = w, so D is rerouted to the literal w of the fresh
parameter w and deleted; G becomes unused (D was its only child) and is
deleted as in D2 below. No equation mentions z (z was isolated with child a
gate), so the affine system is untouched; x and y do not gain children (G is
deleted, D's readers now read w).

**Deletions.**

- **D1 (gate computing a constant or a literal).** A gate J of the acyclic
  part that computes a constant, or a literal of one of its inputs, is
  deleted and all its uses are rerouted to that constant or literal. This
  covers constant and unary gates, gates with a constant input, gates whose
  two inputs are the same value, and AND-type gates whose two inputs are
  parameters or auxiliary variables with dependent (proportional) linear
  parts: in each case J's operation composed with its inputs depends on at
  most one value or is constant, and the multilinear form of a function of
  at most one F_2-variable is constant or a literal. The affine system is
  unchanged; Delta m = Delta q = Delta a = 0; Delta s = -1. An unpaired
  parameter can increase its weight only by 6 -> 7, and only if it was
  isolated with unique child J and J is rerouted to a literal of it and
  thereby hands over at least two children; so Delta W <= -2 + 1 = -1, and
  Delta W <= -2 if J is rerouted to a constant.

- **D2 (unused gate).** A gate J with no children that is not the output is
  deleted. Delta m = Delta q = Delta a = 0, Delta s = -1, and the degrees of
  J's parents (if parameters) each drop by one; so no weight increases, and
  weights may only decrease (a 7 may drop to 6, a 6 to 0). Hence
  Delta W <= -2, with strict inequality when a parent's weight drops.

- **D3 (short equation).** An equation E with at most two nonconstant terms
  must contain an auxiliary variable (its row of A is nonzero by
  invertibility), so E reads v = U + c with v auxiliary and U a parameter,
  an auxiliary variable, or 0. Substitute v <- U + c everywhere, delete v
  and E. In the matrix picture, substituting into another equation adds
  E's row to it — an elementary row operation on A; afterwards the column
  of v has a single nonzero entry (in E's row), so deleting v's column and
  E's row multiplies det A by that single nonzero entry: A stays
  invertible. Term counts do not increase (a term replaced by a constant or
  a literal of one term). Delta s = -1, Delta m = Delta q = Delta a = 0,
  and Delta W <= -1 (only an isolated unpaired parameter whose unique
  child was E can gain, by 6 -> 7, when its rerouted equation-handover
  leaves it with >= 2 children); Delta W <= -2 if no weight increases.

- **D4 (unit parameter).** Let x be an unpaired, isolated parameter whose
  unique child is an equation E. The column N_x of x in N has a single
  nonzero entry, in E's row (x occurs in no other equation), so
  t := A^{-1} N_x != 0; choose an auxiliary variable Q with t_Q = 1. Then
  Q's solved expression is Q = x + phi(u_{-x}), and we apply AR to exchange
  x and Q (Delta W = wt'(w) - 6 <= 7 - 6 = 1, all else unchanged). After
  the exchange the new auxiliary variable x occurs only in E; delete x and
  E exactly as in D3 (single nonzero column entry; the values of all other
  auxiliary variables are unchanged, the remaining equations are the old
  ones minus E, and no gate read x). Delta s = -1, so Delta W <= 1 - 2 =
  -1.

## 4. The gate elimination algorithm

Write a := m - q. The algorithm runs while a > d. Each iteration applies
preparatory operations (moving affine gates into the system; AR; rerouting)
that leave m, s, a unchanged and do not increase W, and then performs one
*reduction*. Throughout, "Delta" refers to the whole iteration.

### 4.1 Preliminary reductions (Cases 1-6)

**Case 1: a deletion is available.** Repeatedly apply D1-D4 until none
applies. Each has Delta m = Delta q = Delta a = 0, Delta s = -1,
Delta W <= -1. If none applies: every gate is XOR- or AND-type with two
distinct nonconstant inputs (D1); every gate other than the output has a
child (D2); every equation has exactly three nonconstant terms (D3); no
unpaired isolated parameter occurs in the affine system (D4).

**Case 2: an AND-type gate on two parameters has an unpaired isolated
input.** Let G be AND-type with inputs z, w, both parameters, z unpaired and
isolated. Fix w (UF if unpaired, PF if paired) to the constant that makes G
constant, delete G by D1, rerouting it to that constant; z becomes unused
(its only child was G), losing weight 6. No other parameter gains: the
rerouting is to a constant.

- w unpaired: Delta m = -1, Delta s = -1, Delta q = 0, Delta a = -1; the
  weights lost are wt(z) = 6 and wt(w) >= 6; Delta W <= -6 - 6 - 2 = -14.
- w paired with r: PF removes w and the pair; r becomes unpaired with
  weight at most 7 (whatever its degree), and no other parameter is
  affected (the rerouting is to a constant): Delta m = -1, Delta s = -1,
  Delta q = -1, Delta a = 0; Delta W <= -3 (pair) - 6 (z) - 2 (G) + 7 (r)
  = -4 <= 14*Delta a = 0.

*Structural consequence of Cases 1-2.* Let z be unpaired and isolated with
unique child J. J is a gate (D4-inapplicable). If J's other input were a
parameter w, then: J XOR-type would have been moved into the affine system
(preparatory step), making z's child an equation — Case 1; J AND-type is
Case 2 (which applies to any AND-type gate on two parameters with an
unpaired isolated input, including the output gate); and if J is the output
and XOR-type reading parameter w, the output is affine in the parameters —
Lemma 2.3. Contradictions all around; so no gate or equation has both an
unpaired isolated parameter and a different parameter among its parents.

**Case 3: a paired parameter has degree 0 or >= 2.** Let p be paired with r.
If deg(p) = 0, apply PF fixing r to a constant: p becomes unpaired of
degree 0 (weight 0), Delta m = -1, Delta q = -1, Delta a = 0,
Delta W = -3. If deg(p) >= 2, apply PF fixing p to a constant and delete
two of p's children (D1-D4 apply: gates reading the constant p become
D1-deletable; equations containing p lose a term, leaving at most two, and
are D3-deletable): Delta m = -1, Delta s = -2, Delta q = -1, Delta a = 0;
no initially unpaired parameter gains weight (its child would need p and it
as parents, contradicting the structural consequence above), and r, newly
unpaired, has weight <= 7: Delta W <= -4 - 3 + 7 = 0.

If Case 3 is inapplicable, every paired parameter is isolated.

**Case 4: two partners have a common child.** Let p, r be partners with a
common child H. If H is an equation, it reads v + p + r = c (it has an
auxiliary variable by invertibility, and exactly three nonconstant terms);
if H is a gate, its two inputs are p and r (it is a common child of both,
and a gate has two pins). Impose p = r + c' by PF, choosing c' as follows:
any c' if H is an equation or an XOR-type gate; c' = al + be + 1 if H is
AND-type (p+al)(r+be)+ga (then (p+al)(r+be) = (r + c' + al)(r + be) equals
r(1 + c' + al + be) + const = const since c' + al + be = 1, using r*r = r).
H becomes constant, is deleted (D1 or D3) and rerouted to its constant; r,
whose only child was H (paired parameters are isolated), becomes unused.
Delta m = -1, Delta s = -1, Delta q = -1, Delta a = 0,
Delta W <= -3 - 2 = -5.

**Case 5: a paired parameter occurs in the affine system.** Let p be paired
(r partner) and occurring in the system; p's unique child H_p is an
equation, and r's unique child H_r != H_p is an equation or gate (Case 4
excluded a common child; r may or may not occur in the system). Fix r to a
constant (PF) and delete H_r (D1 if a gate reading the constant r; D3: the
equation loses the term r, leaving at most two nonconstant terms). H_r has
no unpaired isolated parameter as its other parent (structural consequence:
r is a parameter and a parent of H_r). So far Delta W <= -3 (pair) - 2
(H_r) + 6 (p becomes unpaired, isolated with child H_p: weight 6) = +1.
Then apply D4 to p and H_p: Delta W <= -1 more. In total Delta m = -1,
Delta s = -2, Delta q = -1, Delta a = 0, Delta W <= +1 - 1 = 0.

**Case 6: an AND-type gate on two parameters has a paired input.** Let G be
AND-type with inputs p, w, p paired with r. p is isolated with unique child
G; r is isolated with unique child H_r; r does not occur in the system
(Case 5), so H_r is a gate; H_r != G (Case 4); w != r (else p, r share the
child G: Case 4); w != p (D1: G would read the same value twice); w is not
unpaired isolated (Case 2). So w is unpaired of degree >= 2, or paired.
Fix r to a constant (PF; p becomes unpaired and unused since its only child
G will be deleted), then fix w to the constant making G constant (UF if
unpaired, PF if paired), delete G and H_r by D1 (both become constant
computers; H_r reads the constant r). No initially unpaired parameter gains
weight (G is rerouted to a constant; H_r's other parents are not unpaired
isolated parameters, by the structural consequence applied to H_r, whose
parent r is a parameter).

- w unpaired (weight 7 lost): Delta m = -2, Delta s = -2, Delta q = -1,
  Delta a = -1, Delta W <= -7 - 4 - 3 = -14.
- w paired with r': Delta m = -2, Delta s = -2, Delta q = -2, Delta a = 0;
  r' becomes unpaired with weight <= 7; p is unused (weight 0);
  Delta W <= -4 - 6 + 7 = -3 <= 0.

After Cases 1-6: every parameter occurring in the affine system is unpaired
of degree >= 2, and both inputs of any AND-type gate that reads two
parameters are unpaired of degree >= 2.

### 4.2 Restructuring the minimal AND

Assume a > d and no preliminary reduction applies. The acyclic part is
nonempty: otherwise the output is a constant, a literal of a parameter or
auxiliary variable — an affine function of the parameters — contradicting
Lemma 2.3. Let G be the first gate of the acyclic part in topological
order. Its inputs are constants, parameters, or auxiliary variables; if G's
operation were affine it would have been moved into the system
(preparatory step), so G is AND-type: a *minimal AND*. Its two inputs are
distinct nonconstant values with independent linear parts as affine
functions of the parameters (D1-inapplicability).

*Making the inputs parameters.* If an input is an auxiliary variable Q, its
solved expression is nonconstant (independence), so it depends on some
parameter x, whose column in A^{-1}N is nonzero at Q; x occurs in some
equation (else its column in N, and in A^{-1}N, is zero), so x is unpaired
of degree >= 2 (post-Case-6 structure). Apply AR exchanging Q and x:
Delta W = wt'(w) - wt(x) = 7 - 7 = 0, since the new parameter w has the
children of Q, which include an equation (Q's column in A is nonzero) and
the gate G, so deg'(w) >= 2. Re-run the preparatory moves and Cases 1-6; if
a case applies, the iteration ends with it. Otherwise repeat for the other
input if it is auxiliary: the two inputs of G were linearly independent
before the exchange, and AR is an invertible affine change of parameters,
so they remain independent, and the second auxiliary input depends on some
parameter other than the new w (otherwise the two inputs would be
dependent). After at most two exchanges both inputs are parameters x, y —
unpaired, of degree >= 2 (Cases 2, 6 inapplicable).

*Rerouting children that read x or y.* Let J be a gate child of G reading
also the parameter x (J cannot be an equation: no wire returns from the
acyclic part). J computes h(x, G(x, y)) for some operation h, a function of
the two parameters x, y; relabel J by that function of (x, y) and reroute
its two inputs to x and y. The computed function is unchanged; J stops
being a child of G; only y can gain occurrences, and wt(y) = 7 already
(deg(y) >= 2), so Delta W <= 0. (Symmetrically with x, y exchanged.) After
each rerouting, re-run preparatory moves (moving affine gates into the
system) and inspect Cases 1-6, ending the iteration if one applies. If no
case applies, no deletion has occurred since the previous rerouting: the
only operations between two reroutings are affine-gate moves (no rewiring),
case inspections, and the rerouting itself. The rerouting relabels J — a
child of G — so that J reads x and y instead of x and G: J's own readers
are untouched (they keep reading J), so the only parent-set change to any
value is that G loses the child J, and y possibly gains the child J (J now
reads y). Hence, within the rerouting loop and between its steps, the
number of children of G strictly decreases at each rerouting and never
increases (by (R2), gaining a reader would require deleting and rerouting
a gate that reads G; no deletions occur inside the loop). The loop
therefore terminates. At the end, no child of G reads x or y; G is not the
output (Lemma 2.3) and has k >= 1 children (D2-inapplicable). Swapping
names, assume deg(x) >= deg(y) >= 2.

### 4.3 Eliminating the minimal AND

**Case 7: deg(x) = deg(y) = 2 and x, y have a common second child B.** B
is an equation or a gate distinct from G. If an equation: exactly three
nonconstant terms including an auxiliary variable (invertibility), reading
v + x + y = c. If a gate: its two inputs are x and y. Let h := x + y in the
equation case, else the operation of B. As a function on F_2^2, the
AND-type G is constant on the three points off its exceptional corner; two
of these three points share an h-value, and they span a line
al*px + be*py = ga (al, be) != (0,0) on which both G and h are constant
(an affine line in F_2^2 has exactly two points). Restrict once by UF:
x := ga if be = 0; y := ga if al = 0; x := y + ga if al = be = 1. Now G is
constant; B is constant if it is a gate, and degenerates to v = c_0 if it
is an equation. Delete G (D1) and B (D1, or D3 rerouting v to c_0). All
occurrences of x and y were in G and B, so the parameter not fixed by the
UF becomes unused; both lose weight 7; no other weight changes (all
reroutings are to constants). Delta m = -1, Delta s = -2, Delta q = 0,
Delta a = -1, Delta W <= -7 - 7 - 4 = -18.

**Case 8: deg(x) >= 3, or k >= 2, or the unique child of G has no unpaired
isolated other parent.** Let k_1 = the number of children of G whose other
parent is an unpaired isolated parameter. Fix x to the constant making G
constant (UF). Every child of x is then deletable (D1 for gates, D3 for
equations — an equation with three terms containing the constant x drops
to at most two), and every child of G is a gate reading the constant G,
deletable by D1. The deg(x) + k deleted values are distinct: children of x
and children of G are disjoint, since no child of G reads x (Section 4.2).
Delete them all (starting with G, a child of x). The unpaired parameter x
loses weight 7; each deletion contributes -2. For the weight changes of
other parameters, only (R1) is used: a deletion reroutes the deleted
value's readers to constants or to literals of the deleted value's inputs,
so a parameter can gain children only if it is an input of a deleted
value; the weight cap is 7. Other unpaired parameters:
those of degree >= 2 cannot gain (7 is the cap); the unused stay unused;
an isolated one, z, can only gain (6 -> 7) if its unique child is deleted
and rerouted to a literal of z; its unique child is not a child of x
(structural consequence: it would have z and x as parents), so it is a
child of G, one of the k_1 with other parent z; so the total gain is at
most k_1. Therefore

    Delta m = -1, Delta s = -(deg(x) + k), Delta q = 0, Delta a = -1,
    Delta W <= -7 - 2(deg(x) + k) + k_1.

If deg(x) = 2 then deg(y) = 2 and, writing B_x, B_y for the second children
(B_x != B_y by Case 7-inapplicability; y is not a parent of B_x since y's
children are G and B_y), no child of G reads y, and G is rerouted to a
constant, so y loses the child G and gains none: y becomes isolated and its
weight drops from 7 to 6, one further unit. Hence

    deg(x) >= 3:  Delta W <= -7 - 6 - 2k + k <= -13 - k <= -14  (k >= 1);
    deg(x) = 2, k >= 2:  Delta W <= -12 - 2k + k_1 <= -12 - k <= -14;
    deg(x) = 2, k = 1:  Delta W <= -14 + k_1,  i.e. -14 if k_1 = 0.

So Case 8 satisfies Delta W <= 14 Delta a = -14 unless deg(x) = deg(y) = 2,
k = k_1 = 1: the *bottleneck configuration* — G has the unique child D,
whose other parent is an unpaired isolated parameter z. D is AND- or
XOR-type (D1-inapplicable). These two cases follow.

**Case 9: bottleneck, D AND-type.** D = (G + de)(z + eta) + th. Pick x_0,
y_0 with G(x_0, y_0) = de (G is AND-type, hence nonconstant on F_2^2).
Fix x := x_0 and y := y_0 (two UFs). Then G = de and D = th are constants;
delete G and D (D1); delete the second children B_x, B_y of x and y (D1 or
D3: x, y are now constants), which are four distinct values: B_x != B_y
(Case 7), and neither equals D (no child of G reads x or y; D reads G and
z, and z != x, y). The parameter z becomes unused (its child D is rerouted
to a constant); x and y lose weight 7 each; no unpaired parameter gains
(B_x, B_y do not have unpaired isolated other parents — x is a parent of
B_x, y of B_y, and the structural consequence applies; G, D reroute to
constants). Delta m = -2, Delta s = -4, Delta q = 0, Delta a = -2,
Delta W <= -7 - 7 - 6 - 8 = -28 = 14 Delta a.

**Case 10: bottleneck, D XOR-type.** D = G + z + c. This is exactly the QR
configuration of Section 3: apply QR with z <- w + G(x, y) + c, add the
pair {x, y}, reroute D to w and delete it, then delete the unused G (D2).
No equation mentions z; no other unpaired parameter gains weight (D's
readers now read w; G is deleted unrouted). Delta m = 0, Delta s = -2,
Delta q = +1, Delta a = -1; weights: x and y lose 7 each and gain pair
status (+3 for the pair), z's weight 6 is replaced by w's weight
(deg(w) = deg of D's former children, at most 7). Delta W <= -4 - 7 - 7
- 6 + 3 + 7 = -14 = 14 Delta a.

### 4.4 End of the algorithm and the count

**Lemma 4.1 (one iteration).** If a > d, one iteration transforms the state
with Delta a in {0, -1, -2}, Delta W <= 14 Delta a, and m + s strictly
decreasing.

*Proof.* The iteration performs preparatory operations (Delta m = Delta s =
Delta a = 0, Delta W <= 0: moving affine gates preserves everything; AR is
used only with Delta W = 0 here; a relabeling rerouting can only give the
child J to a parameter already of weight 7) and then one
reduction: Cases 1-6 if any applies, else the restructuring of Section 4.2
(each of whose preparatory sub-steps — AR with deg(x) >= 2 as in Section
3, affine-gate moves, and the relabeling reroutings, where only a
parameter already of weight 7 can gain a child — either ends the iteration
in a preliminary case or has Delta W <= 0 and unchanged m, s, a), else one
of Cases 7-10. The tables above give Delta a in {0,-1,-2} and
Delta W <= 14 Delta a in every case. Deletions and QR decrease s with m
unchanged; every other reduction decreases m with s not increased; so
m + s strictly decreases. (The procedure is finite, and it is easily
implemented in polynomial time: all tests are lookups of local structure,
degree counters, truth tables of at most 16 entries, and Gaussian
elimination over F_2 for the solved expressions.) qed

**Proof of the Main Theorem.** Initialize: Phi = id, P empty, empty affine
system, C the given circuit; m = n, q = 0, s = |C|, a = n. Initially every
parameter is unpaired with weight <= 7, so W_0 <= 2|C| + 7n. While a > d,
iterate Lemma 4.1. The algorithm halts with a_T <= d after
T <= n + |C| iterations (m + s decreases, both nonnegative). Summing:

    W_T - W_0 = sum Delta W <= 14 sum Delta a = 14 (a_T - n).

Since W_T >= 0 and a_T <= d: W_0 >= 14(n - d), so

    2|C| + 7n >= W_0 >= 14(n - d)  =>  |C| >= (7/2) n - 7 d.

This proves the Main Theorem. qed

## 5. Consistency with the pinned model conventions

For completeness we record why the bound applies to the original acyclic
B_2-circuit as defined in Section 0. Fair semi-circuits, parameters, and the
map Phi are devices of the proof; the initial state is the original circuit
with Phi the identity and an empty affine system, and the initial potential
is W_0 <= 2|C| + 7n where |C| counts every gate of the original circuit
exactly once, including trivial, unary, and negated-projection gates
(negations on wires being free, a retained unary gate is charged while
present, and is removed — with a charge of -1 to s — by D1 when it arises).
Constants introduced by substitutions are free (they are never gates), and
constant gates already present in C are charged and deleted by D1/D2.
Substitutions never turn a witness bit into a circuit input: the circuit
always reads only the m current parameters, and the invariant C(u) =
f(Phi(u)) pins the computed function to the disperser restricted to the
current image, whose affine dimension a = m - q is at least d throughout,
keeping the restricted function nonconstant (Proposition 2.1). The final
inequality |C| >= (7/2)n - 7d uses only the initial potential and the
invariant, so it applies to every B_2-circuit computing f, uniform or not.
