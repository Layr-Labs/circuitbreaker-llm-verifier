# The class theorem: every affine disperser needs 3.5n - 7d gates

This file proves the following statement, which is the mathematical core of the
submission. It is a complete reproduction, in self-contained notation, of the
gate elimination theorem for affine dispersers obtained by Mishra and Santhanam
(ECCC TR26-260, 19 October 2026, Theorem 4.1). No new theorem is claimed here;
what is supplied is a full proof with every case and every arithmetic step
written out, so that the argument can be checked without retrieving any source
other than the two quoted classical inputs listed at the end.

**Theorem A.** Let f : {0,1}^n -> {0,1} be an affine disperser for dimension d,
that is, f is non-constant on every affine subspace of F_2^n of dimension d.
Then every B2 circuit C computing f has size at least 3.5n - 7d.

The B2 model used is the one fixed by rules.md: acyclic circuits with fan-in two,
one output, any of the sixteen binary Boolean operations at each gate,
unrestricted fan-out and depth, each retained gate counted once, inputs, wires,
constants and negations free, agreement with f required on every input.

## 1. Notation and normalisation

Identify {0,1} with F_2, so that addition is XOR, multiplication is AND, and
negation is addition of 1. Every operation h : F_2^2 -> F_2 has a unique
multilinear form h(t_1,t_2) = c + b_1 t_1 + b_2 t_2 + b_12 t_1 t_2. A gate is
called *constant* if its operation ignores both inputs, *unary* if it depends on
exactly one input, *XOR-type* if h(t_1,t_2) = b_1 t_1 + b_2 t_2 + c with
b_1 = b_2 = 1, and *AND-type* if h(t_1,t_2) = (t_1 + alpha)(t_2 + beta) + gamma
with constants alpha, beta, gamma. Of the sixteen operations two are constant,
four are unary, two are XOR-type and eight are AND-type. Negations on wires are
absorbed into the operation of the gate that reads the wire, at no cost, and
constants are free; this is the Li-Yang normalisation required by rules.md.

A gate is *trivial* if it computes a constant or a literal of one of its
parents; it can be deleted while rerouting its uses, at no change of the computed
function and without adding gates, because wires cost nothing. Throughout, the
current circuit is kept normalised: every gate is XOR-type or AND-type and has two
distinct parents, neither of which is a constant.

## 2. Parameters

Gate elimination proceeds on a space of independent parameters rather than on the
input variables. A *state* consists of

* an injective map Phi : F_2^m -> F_2^n of the form
  Phi(u) = b + L u + sum over {p,r} in P of c_{p,r} p r,
  where b in F_2^n, L in F_2^{n x m}, P is a set of disjoint unordered pairs of
  parameters, and c_{p,r} in F_2^n;
* a circuit C, represented below as a fair semi-circuit, that reads the
  parameters u_1..u_m and satisfies the main invariant
  C(u) = f(Phi(u)) for every u in F_2^m.

A parameter in a pair of P is *paired*, its mate is its *partner*, and the
remaining parameters are *unpaired*. A pair contributes only one degree of
freedom, since the quadratic term p r becomes affine once p is fixed to a
constant or identified with a literal of r. We track the state by four numbers:
s = |C| (the current size), m = the number of parameters, q = |P|, and
a = m - q.

**Proposition 2.1 (non-triviality).** The image of every injective map of the
above form contains an affine subspace of dimension a = m - q. Consequently
f o Phi is non-constant whenever a >= d.

*Proof.* There are q pairs and q quadratic monomials, pairwise involving disjoint
parameters. Fix one parameter of each pair to an arbitrary constant. This fixes
exactly q parameters, every quadratic term becomes affine in the remaining m - q
parameters, and the resulting map is still injective (it is Phi composed with the
injective map that fixes those coordinates). Its image is therefore an affine
subspace of dimension m - q. If a >= d then f is non-constant on that subspace,
hence on the larger set im(Phi). []

**Lemma 2.2 (no simple output).** Let a > d. Then the output of C is neither a
constant, nor a non-constant affine function of the parameters, nor a literal of
an AND-type gate whose two parents are unpaired parameters.

*Proof.* A constant output contradicts Proposition 2.1. Suppose the output is a
non-constant affine function of the parameters and pick a parameter x with
nonzero coefficient. If x is paired, first fix its partner to a constant: this
decreases m and q by one each, so a is unchanged, and it makes x unpaired. Now substitute
for x the affine expression in the other parameters that makes the output
constant; this is a (UF) restriction of Section 4, it decreases a by one, and the
restricted image still contains an affine subspace of dimension at least a - 1 >= d,
on which f would then be constant, contradicting Proposition 2.1. If instead the
output is a literal of an AND-type gate on two unpaired parameters x, y, fix y to
the constant that makes that AND-type gate constant; the output becomes constant
while a decreases by only one, the same contradiction. []

## 3. Fair semi-circuits

A *fair semi-circuit* on the parameters u consists of

1. an affine system A v = N u + c in r auxiliary variables v_1..v_r, where
   A in F_2^{r x r} is invertible (fairness: a unique solution v for every u) and
   every row of (A|N) has at most three nonzero entries;
2. an acyclic part, an ordinary circuit whose gates are listed in topological
   order and read constants, parameters, auxiliary variables, or earlier gates;
3. an output, which is a constant or a literal of a parameter, an auxiliary
   variable, or a gate.

No variable of the affine system is designated as an output; this is what allows
exchanging the roles of an auxiliary variable and a parameter. The *size* of the
semi-circuit is s = r + (number of gates of the acyclic part). Initially the
affine system is empty and the acyclic part is the original circuit, so s = |C|.

**Moving a gate into the system.** If a gate J of the acyclic part has an affine
operation h(t_1,t_2) = b_1 t_1 + b_2 t_2 + c whose parents are constants,
parameters or auxiliary variables, replace J by the equation
J + b_1 t_1 + b_2 t_2 = c and declare J a new auxiliary variable. The new
equation has at most three nonconstant terms, the new auxiliary matrix is
[[A,0],[*,1]] and hence invertible, and the function computed is unchanged since
every use of J now reads the auxiliary variable J. One gate is replaced by one
equation, so s is unchanged. Applying this in topological order moves every gate
with an affine operation and no AND-type ancestor into the affine system.

**Degrees and rerouting.** The *degree* deg(t) of an internal value t is its
number of distinct children, where the children of a value are the equations that
contain it with a nonzero coefficient and the gates that read it (a gate reading
the same value twice counts once). A parameter is *unused* if its degree is zero
and *isolated* if its degree is one. When a gate or auxiliary variable becomes a
constant or a literal of another internal value, every occurrence of it is
replaced by that constant or literal (*rerouting*) and it is deleted. Rerouting
can increase the degree of a parameter that was a parent of the replaced value.

**The potential.** The weight of an unpaired parameter x is

    wt(x) = 0 if deg(x) = 0,   wt(x) = 6 if deg(x) = 1,   wt(x) = 7 if deg(x) >= 2,

and paired parameters have no individual weight. The potential is

    W = 2s + 3q + sum over unpaired x of wt(x),

which is always nonnegative. For any quantity X, write X' for its value after an
operation and Delta X = X' - X; a *step* may group several operations, and then
Delta refers to the beginning and the end of the whole step.

**Preservation of the representation.** Fixing a parameter to a constant deletes
its column from N and changes c; imposing x = y + c changes neither A nor the
number of nonconstant terms. An affine substitution for a parameter *inside the
system* could expand one term into many and is therefore not performed directly:
it is implemented by exchanging the roles of a parameter and an auxiliary
variable, which leaves every equation unchanged and keeps A invertible (the new
auxiliary matrix is A times the identity with one column replaced by t = A^{-1}N_x,
whose relevant entry is 1, so the determinant is preserved). Rerouting an
auxiliary variable v that solves v = U + c: replacing v by U + c in every equation
is an elementary row operation (determinant preserved), after which the column of
v has a single nonzero entry and deleting that row and column preserves
invertibility. Deleting an auxiliary variable whose only use is a single equation,
together with that equation, is the same determinant expansion. Hence every
operation below keeps the two invariants (at most three nonconstant terms per
equation, invertible auxiliary matrix).

## 4. The four changes of parameters

Let u - x denote all parameters other than x, and let w be a fresh parameter.

**(AR) affine reparameterisation.** For an unpaired parameter x and an affine map
phi, substitute x <- w + phi(u - x). This is a bijection of the parameter space,
so im(Phi) is unchanged, and m, q, a, s are unchanged. If the new parameter w
takes over the children of an auxiliary variable Q with deg(Q) >= 2, and deg(x) >= 2,
then Delta W = 0.

**(QR) quadratic reparameterisation.** For distinct unpaired parameters x, y, z and
constants alpha, beta, gamma, substitute z <- w + (x + alpha)(y + beta) + gamma
and add the pair {x,y} to P. This is again a bijection. It is applied only in the
following configuration: the circuit contains the AND-type gate
G = (x + alpha)(y + beta) + gamma, an XOR-type gate D = G + z + c, z is an
isolated parameter whose only child is D, and the only child of G is D. After the
substitution D computes w, so every use of D is rerouted to w and D is deleted;
G then has no child and is deleted as unused. Thus Delta m = 0, Delta q = +1,
Delta a = -1, Delta s = -2, and, since x and y lose their weights, the new pair
contributes +3, the two deletions contribute -4, and w replaces the isolated
parameter z of weight 6 by a parameter of weight at most 7,

    Delta W <= -4 + 3 + 7 - 6 - wt(x) - wt(y) = -wt(x) - wt(y) <= -14,

the last inequality holding because (QR) is used only when deg(x), deg(y) >= 2.

**(UF) fixing an unpaired parameter.** For an unpaired parameter x and an affine
map psi, impose x = psi(u - x); psi may be constant. Then Delta m = -1, Delta q = 0,
Delta a = -1, and no gate is added: each occurrence of x becomes a constant or a
literal of another parameter, and the gates and equations so simplified are
counted in s until they are deleted.

**(PF) fixing a paired parameter.** For a pair {p,r}, either fix p to a constant
or impose p = r + c, and remove p and the pair. Then Delta m = -1, Delta q = -1,
Delta a = 0.

## 5. Deletions

**D1 (gate depending on at most one input).** If a gate J computes a constant or a
literal of one of its parents, reroute every use of J and delete J. Then
Delta m = Delta q = Delta a = 0 and Delta s = -1. An unpaired parameter that was
a parent of J can increase its weight from 6 to 7 only if it was isolated and has
at least two children after the rerouting; hence Delta W <= -2 - 6 + 7 = -1, and
Delta W <= -2 if J is rerouted to a constant (then no weight can increase).

**D2 (unused gate).** A gate with no children that is not the output is deleted:
Delta s = -1 and Delta W <= -2.

**D3 (equation with at most two nonconstant terms).** Such an equation must
contain an auxiliary variable (otherwise its row in the invertible auxiliary
matrix would vanish), so it reads v = U + c. Substitute v = U + c everywhere and
delete v and the equation: Delta s = -1, Delta m = Delta q = Delta a = 0, and
Delta W <= -1, with Delta W <= -2 if no unpaired parameter's weight increases.

**D4 (unit parameter deletion).** Let x be an unpaired, isolated parameter whose
unique child is an equation E. Because the column of x in N has a single nonzero
entry, the vector A^{-1}N_x is nonzero; choose an auxiliary variable Q with
coordinate 1 in it, so Q = x + phi(u - x) and Q depends on x. Apply (AR) with this
Q, which exchanges the roles of x and Q without changing any equation, and then
delete x (now an auxiliary variable) and E. Then Delta m = Delta q = Delta a = 0,
Delta s = -1, the weight 6 of x is replaced by a weight of at most 7 for the new
parameter Q, no other unpaired parameter increases its weight, so
Delta W <= -2 - 6 + 7 = -1.

## 6. Preliminary reductions

Each iteration of the algorithm first moves every gate with an affine operation and
no AND-type ancestor into the affine system; this keeps m, q, a, s and never
increases W. Then the following cases are tested in order and the first applicable
one is performed. Consequently, when Case j is applied, all earlier cases are
known to be inapplicable.

**Case 1: a deletion is available.** Apply D1, D2, D3, D4 repeatedly until none
applies. Each such deletion has Delta m = Delta q = Delta a = 0, Delta s = -1 and
Delta W <= -1. If no deletion applies, then: every gate is XOR-type or AND-type
with two distinct non-constant parents, and two parents of an AND-type gate that
are both parameters or auxiliary variables have independent linear parts (if their
linear parts were dependent, the second parent would be the first plus a constant,
and the gate would be a unary or constant function of one value, which is exactly
what D1 deletes); every gate other than the output has a child; every equation has
exactly three nonconstant terms; no unpaired isolated parameter occurs in the
affine system.

**Case 2: an AND-type gate on parameters has an unpaired, isolated input.** Let G
be an AND-type gate with parents z and w, where z is an unpaired, isolated
parameter. Fix w to the constant that makes G constant, delete G by D1 and reroute
it to that constant; since G was the only child of z, the parameter z becomes
unused, and the degree of every other parameter is unchanged because G was
rerouted to a constant. If w is unpaired, fix w to a constant by (UF): both z and w
had weight at least 6, so Delta m = -1, Delta s = -1, Delta q = 0, Delta a = -1 and

    Delta W <= -6 - 6 - 2 = -14.

If w is paired, fix w by (PF) and let r be its partner, which becomes unpaired with
weight at most 7: Delta m = -1, Delta s = -1, Delta q = -1, Delta a = 0 and
Delta W <= -3 - 6 - 2 + 7 = -4.

*Consequence (used in Cases 3, 6, 8 and 9).* Assume Cases 1 and 2 are inapplicable.
Then an unpaired isolated parameter z is not the output (Lemma 2.2, which applies
because a > d), does not occur in the affine system (else D4 applies), and is not
a parent of a gate together with another parameter (else the gate is affine with
two parameter parents, hence already moved into the system, or Case 2 applies).
Consequently, when a gate or equation is deleted in the cases below, no unpaired
isolated parameter can be rerouted to a literal of itself, so no parameter's weight
rises from 6 to 7 unless a case explicitly counts such a rise; this is what makes
the term -2 per deletion exact in Cases 3, 6, 8 and 9. Consequently the unique child
of such a z is a gate G of the acyclic part (it is not an equation, and it exists
because z is not the output), and no gate and no equation has both an isolated
unpaired parameter and a different parameter as parents.

**Case 3: a paired parameter has degree zero or at least two.** Let p be paired
with partner r. If deg(p) = 0, fix r by (PF): the pair is removed and p becomes an
unpaired parameter of degree zero, so Delta m = -1, Delta q = -1, Delta a = 0,
Delta s = 0 and Delta W <= -3. If deg(p) >= 2, fix p by (PF) and delete two of its
children by D1-D4: Delta m = -1, Delta s = -2, Delta q = -1, Delta a = 0, and only
r, which becomes unpaired, can gain weight (at most 7), so
Delta W <= -3 - 4 + 7 = 0.

*Consequence.* If Case 3 is inapplicable, every paired parameter is isolated and
has exactly one child.

**Case 4: two partners have the same child.** Let p, r be partners with a common
child H. If H is an equation it reads v + p + r = c for an auxiliary variable v; if
H is a gate its two parents are p and r. Impose p = r + c' by (PF), choosing c' so
that H becomes constant (any c' works if H is an equation or XOR-type; for
H = (p + alpha)(r + beta) + gamma take c' = alpha + beta + 1). Delete H and reroute
it; r survives but becomes unused because H was its only child. Then Delta m = -1,
Delta s = -1, Delta q = -1, Delta a = 0 and Delta W <= -3 - 2 = -5.

**Case 5: a paired parameter occurs in the affine system.** Let p be paired with
partner r, and let p occur in the affine system. By Case 3, p has a unique child
H_p, which must be an equation; similarly r has a unique child H_r, and H_r != H_p
by Case 4. Fix r by (PF) and delete H_r (by D1 if a gate, by D3 if an equation).
By Cases 1 and 2, H_r has no unpaired isolated parameter as a parent, so no
initially unpaired parameter gains weight. Now p is an unpaired isolated parameter
whose unique child is H_p, so it has weight 6; removing the pair contributes -3,
deleting H_r contributes -2, and the new weight contributes +6, giving
Delta W <= +1 so far. Then apply D4 to p and H_p, which gives Delta W <= -1. Hence
Delta m = -1, Delta s = -2, Delta q = -1, Delta a = 0 and Delta W <= 0.

**Case 6: an AND-type gate on parameters has a paired input.** Let G be an
AND-type gate with parents p and w, where p is paired with partner r. By Case 3, G
is the unique child of p and r has a unique child H_r; by Case 5, r does not occur
in the affine system, so H_r is a gate; by Case 4, H_r != G; by Case 1, w != p; by
Case 2, w is not an unpaired isolated parameter, so w is either an unpaired
parameter of degree at least two or a paired parameter. Fix r by (PF), then fix w
to the constant that makes G constant (by (UF) if w is unpaired, by (PF) if w is
paired), and delete G and H_r by D1. After r is restricted, p becomes unpaired and
unused. If w was unpaired its degree was at least two, so it contributed +7 and
now contributes nothing; hence Delta m = -2, Delta s = -2, Delta q = -1,
Delta a = -1 and Delta W <= -7 - 3 - 4 = -14. If w was paired with partner r', both
pairs are removed, only r' can be newly unpaired with weight at most 7, and
Delta m = -2, Delta s = -2, Delta q = -2, Delta a = 0,
Delta W <= -6 - 4 + 7 = -3.

*Consequence.* If Cases 1-6 are inapplicable, every parameter occurring in the
affine system is unpaired with degree at least two, and every AND-type gate whose
parents are parameters has two unpaired parents of degree at least two.

## 7. Restructuring the minimal AND

Assume a > d and that no preliminary reduction applies. The acyclic part has at
least one gate, since otherwise the output would be an affine function of the
parameters, contradicting Lemma 2.2. Let G be the first gate of the acyclic part in
topological order, the *minimal AND*. All affine gates whose parents both lie in
the affine system have been moved into the system, so G is AND-type; by Case 1 its
two parents are parameters or auxiliary variables with independent linear parts.

**Parents to parameters.** Suppose a parent of G is an auxiliary variable Q. Then
Q is a non-constant affine function of the parameters and depends on some parameter
x; x occurs in some equation, because otherwise the column of x in A^{-1}N would
vanish. By the consequence of Cases 1-6, x is unpaired with degree at least two.
Apply (AR) exchanging Q and x. The new parameter w replacing Q has degree at least
two (it occurs in an equation and is a parent of G), and x has degree at least two,
so both weights are 7 and Delta m = Delta s = Delta q = Delta a = Delta W = 0. Then
move affine gates into the system again and retest the preliminary reductions; if
one applies, perform it and end the iteration. If both parents of G are auxiliary
variables, repeat the exchange for a parameter on which the second one depends
(it exists because the two parents were linearly independent and an invertible
affine change of parameters preserves that). After at most two applications of
(AR), both parents of G are parameters x and y, and by Cases 2 and 6 they are
unpaired with deg(x), deg(y) >= 2.

**Rerouting children of G.** Suppose a child J of G also has x as a parent. Since
J must be a gate of the acyclic part, its value is h(x, G(x,y)) for some
h in B2; relabel J by the operation (x,y) -> h(x, G(x,y)) and reroute its parents
to x and y. The function computed is unchanged, J is no longer a child of G, and
only y can gain occurrences; since deg(y) >= 2, its weight stays 7. Hence
Delta m = Delta s = Delta q = Delta a = 0 and Delta W <= 0. Each rerouting removes
a child of G and creates none, so the process terminates. Afterwards no child of G
reads x or y, G is not the output (Lemma 2.2). If G has no child at that point,
then it is an unused gate and D2 applies, which ends the iteration with Delta a = 0
and Delta W <= -2; otherwise G has at least one child, so k >= 1 below.

## 8. Eliminating the minimal AND

Swap x and y so that deg(x) >= deg(y) >= 2, and let k >= 1 be the number of
children of G.

**Case 7: deg(x) = deg(y) = 2 and x, y have a common second child B != G.**
If B is an equation, by Case 1 it has exactly three nonconstant terms and at least
one auxiliary variable, so it reads v + x + y = c; set h(x,y) = x + y. If B is a
gate, its parents are x and y; set h to its operation. As a function of (x,y), the
AND-type gate G takes the same value at three points of F_2^2, and at least two of
those points have equal h-values, so those two points form a line
alpha x + beta y = gamma on which both G and h are constant. Restrict to that line
by one (UF): fix x = gamma if beta = 0, fix y = gamma if alpha = 0, or impose
x = y + gamma if alpha = beta = 1. Then G is constant and B is constant (or, if an
equation, becomes v = c_0). Delete G and B by D1 or D3, rerouting them to their
constant values. The restriction removes one of x, y and the other becomes unused,
so both weights 7 are lost and no other unpaired parameter gains weight. Hence
Delta m = -1, Delta s = -2, Delta q = 0, Delta a = -1 and
Delta W <= -7 - 7 - 4 = -18.

From now on Case 7 is inapplicable; in particular, if deg(x) = deg(y) = 2 then the
second children B_x of x and B_y of y are distinct.

**Case 8: deg(x) >= 3, or k >= 2, or the unique child of G has no unpaired
isolated input.** Let k_1 be the number of children of G whose other parent is an
unpaired isolated parameter, so 0 <= k_1 <= k. Fix x to the constant that makes G
constant (a (UF)), delete all deg(x) children of x starting with G, and then all k
children of G. These deg(x) + k internal values are distinct because no child of G
has x as a parent. Each deletion lowers W by 2. The parameter x loses its weight 7.
An unpaired parameter of degree at least two already has the maximal weight 7 and
cannot gain. An initially unused parameter stays unused. An initially isolated
parameter z can gain weight only if its unique child is deleted and rerouted to a
literal of z; by the consequences of Cases 1 and 2 the unique child of z is not a
child of x, so it must be one of the k_1 children of G whose other parent is
unpaired and isolated. Hence at most k_1 parameters gain exactly 1. Therefore
Delta m = -1, Delta s = -(deg(x) + k), Delta q = 0, Delta a = -1 and

    Delta W <= -7 - 2(deg(x) + k) + k_1.

If deg(x) = 2, then also deg(y) = 2, and there is one more decrease: the children
of y are G and B_y with B_y != B_x, no child of G has y as a parent, so after G is
rerouted to a constant y loses G, gains nothing, and becomes isolated; its weight
drops from 7 to at most 6. Combining,

    Delta W <= -7 - 2 deg(x) - k                <= -14   if deg(x) >= 3,
    Delta W <= -12 - 2k + k_1                    <= -14   if deg(x) = 2 and k >= 2,
    Delta W <= -14 + k_1                                  if deg(x) = 2 and k = 1.

Since Delta a = -1, Case 8 satisfies Delta W <= 14 Delta a unless
deg(x) = deg(y) = 2 and k = k_1 = 1. That remaining situation is the *bottleneck
configuration*: G has a unique child D whose other parent is an unpaired isolated
parameter z. It is checked before Case 8 is applied, and then one of the two
following cases is used instead. Since Case 1 is inapplicable, D is AND-type or
XOR-type.

**Case 9: the bottleneck configuration with D of AND-type.** Write
D = (G + delta)(z + eta) + theta. Since G is non-constant, choose x_0, y_0 with
G(x_0,y_0) = delta. Fix x = x_0 and y = y_0 by two (UF) restrictions; then G
computes delta and D computes theta. Delete G and D by D1, rerouting them to those
constants, and delete B_x and B_y by D1 or D3, which is possible since x and y are
fixed. The four deleted values are distinct: B_x != B_y, and neither equals D since
no child of G has x or y as a parent. The only child of z was D, so z becomes
unused. No unpaired parameter gains weight: G and D are rerouted to constants, and
neither B_x nor B_y has an unpaired isolated parameter as a parent, because x is a
parent of B_x, y is a parent of B_y, and Cases 1 and 2 are inapplicable. Thus
x, y, z lose weights 7, 7, 6 and the four deletions contribute -8:

    Delta m = -2, Delta s = -4, Delta q = 0, Delta a = -2,
    Delta W <= -7 - 7 - 6 - 8 = -28 = 14 Delta a.

**Case 10: the bottleneck configuration with D of XOR-type.** Write D = G + z + c
and apply (QR) with z <- w + G(x,y) + c, adding the pair {x,y}: D computes w, so D
is rerouted to w and deleted, and G becomes unused and is deleted by D2. No
equation and no other gate is affected, because the only child of z was D. Two
deletions contribute -4, the new pair contributes +3, x and y lose 7 each, and the
new parameter w (whose children are those of D) has weight at most 7 and replaces
the isolated parameter z of weight 6. No other unpaired parameter gains. Hence
Delta m = 0, Delta s = -2, Delta q = +1, Delta a = -1 and

    Delta W <= -4 - 7 - 7 - 6 + 3 + 7 = -14 = 14 Delta a.

## 9. The iteration lemma and Theorem A

**Lemma 9.1.** Suppose a > d. In polynomial time the algorithm transforms the
current state into a new state with m + s strictly smaller, Delta a in {0,-1,-2},
and Delta W <= 14 Delta a.

*Proof.* If the current state admits a preliminary reduction, Cases 1-6 give
Delta a in {0,-1} and Delta W <= 14 Delta a (Case 2 with w unpaired and Case 6 with
w unpaired have Delta a = -1 and Delta W <= -14; the other subcases have Delta a = 0
and Delta W <= 0). Otherwise Section 7 produces, using only (AR), moves into the
system and reroutings, a state with the same m, s, a and no larger W, in which the
minimal AND has two unpaired parameter parents of degree at least two and no child
reading them. Then Case 7 (Delta a = -1, Delta W <= -18), Case 8 (Delta a = -1,
Delta W <= 14 Delta a), Case 9 (Delta a = -2, Delta W <= -28) or Case 10
(Delta a = -1, Delta W <= -14) applies. In every case Delta W <= 14 Delta a.
Finally, the deletions and (QR) leave m unchanged, every other reduction decreases
m, and no reduction increases s; hence m + s strictly decreases. Each iteration
performs at most s moves into the system, at most two exchanges, and at most one
rerouting per child of the minimal AND; Gaussian elimination over F_2 expresses
auxiliary variables as affine functions of the parameters and supports the tests of
linear dependence and the (AR) exchanges; all remaining tests involve the lists of
gates and equations and Boolean truth tables of constant size. Hence each iteration
runs in polynomial time. []

*Proof of Theorem A.* Initialise Phi as the identity on F_2^n, with no pairs, an
empty affine system and the acyclic part equal to C. Then m = a = n, q = 0 and
s = |C|, and since every input is an unpaired parameter of weight at most 7,

    W_0 <= 2|C| + 7n.

While a > d apply Lemma 9.1. Each iteration strictly decreases m + s, which starts
at n + |C| and stays nonnegative, so the process halts with a_T < d after finitely
many iterations. Summing Delta W <= 14 Delta a over the T iterations gives
W_T - W_0 <= 14(a_T - n), and since W_T >= 0 and a_T <= d,

    2|C| + 7n >= W_0 >= W_T + 14(n - a_T) >= 14(n - d).

Rearranging, |C| >= 3.5n - 7d. []

## 10. On the constant

The proof uses a potential of the form V = s + beta q + alpha I - sigma L, where I
counts used unpaired parameters and L counts unpaired isolated parameters, with
alpha = 7/2, beta = 3/2, sigma = 1/2; multiplied by 2 this is the potential W used
above. Within this three-parameter family no choice of weights improves the leading
coefficient 7/2: the binding estimates are Case 2 with w unpaired, the second
subcase of Case 6 with w unpaired, the subcase deg(x) >= 3 with k = 1 of Case 8,
and Cases 9 and 10, all of which attain equality 14 Delta a. A stronger coefficient
therefore requires a genuinely different elimination step and not a reweighting of
this potential; no such step is asserted in this submission.

## 11. Quoted inputs

Two classical statements are used, and nothing else.

1. **Bertrand's postulate** (Chebyshev; elementary proof by Erdős, 1932): for every
   real x >= 2 there is a prime p with x < p < 2x. Used only to show that the least
   prime k > 2m^{3/5} satisfies k < 4m^{3/5}, i.e. that the field degree is
   O(m^{3/5}) and the scan for k is polynomial in m.
2. **The existence of an explicit affine disperser for sublinear dimension**, in
   the precise form of Theorem 2.2 of Ben-Sasson and Kopparty, *Affine dispersers
   from subspace polynomials*, SIAM J. Comput. 41(4):880-914, 2012 (author preprint
   of 12 March 2010, Theorem 2.2 and equation (4); proof in Sections 3, 6 and 7 of
   that preprint). Its statement, its hypotheses as met by the family of
   `proof.md`, and a reproduction of its construction are in `disperser.md`.

No other external premise is used: the gate elimination theorem above is proved in
full, and the arithmetic in Section 8 and Section 9 is verified exactly by
`check_cases.py`.
