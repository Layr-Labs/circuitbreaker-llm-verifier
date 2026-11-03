# The gate-elimination theorem: `4n - 10d` gates for every dimension-`d` affine disperser

This file proves, from its own definitions, the lower-bound engine used by the
bound:

> **Theorem 10.1.** Let `f : F_2^n -> F_2` be an affine disperser for dimension
> `d` (Definition 1.1 of `family.md`), where `0 <= d < n`. Then every B2 circuit
> computing `f` has at least `4n - 10d` gates.

The proof is a self-contained exposition of the gate-elimination argument of
Mishra and Santhanam ("An Explicit `4n - O(log n)` Size Lower Bound against
General Boolean Circuits", ECCC TR26-260, revision 1, 3 November 2026;
SHA256 of the revision PDF `a123943ed05c9569385d3b9d98c695630a110559f33a4c990ec44e2a22b4d921`).
The attribution is recorded in Section 12; every step below is stated and
argued here, so the file can be audited without the source. The one idea that
is theirs and not classical is in Section 4: substitutions are performed on a
space of *parameters* that evolves by affine reparameterisation, which lets the
argument pair two inputs of a minimal AND-gate and *fix the value of that gate
without performing any substitution*. That is what removes the bottleneck case
analysis that kept earlier proofs at `3n`--`3.1n`.

Throughout, `n` is fixed once and for all. All statements are for this single
`n`; the passage to all large `n` happens only in `proof.md`, where the loss
`10 d(n)` with `d(n) = o(n)` is evaluated.

## 1. Conventions

`F_2 = {0,1}`, addition is XOR, multiplication is AND, negation is `+1`. All
linear algebra in this file is over `F_2`. The input variables are
`X = (X_1, ..., X_n)`.

The circuit model is the one fixed by the challenge rules (Li--Yang's B2 model):

- a circuit is a finite directed acyclic graph with `n` input nodes of in-degree
  zero and gates of in-degree exactly two; one gate (or a literal, or a
  constant) is the designated output;
- every gate is labelled by one of the sixteen functions `F_2^2 -> F_2`; all
  sixteen are available;
- depth and fan-out are unrestricted, sharing is unrestricted, and the gate may
  read the same value twice;
- the circuit must compute the target function exactly, on every input;
- `|C|` counts gates only: input nodes, wires and the output selection cost
  nothing, and every retained gate is counted once, degenerate or not;
- negations on wires are free and are absorbed into the operation of the gate
  that reads the wire, which is why the output may be a literal.

Classify the sixteen operations. Writing `h(t_1,t_2) = c + b_1 t_1 + b_2 t_2 +
b_{12} t_1 t_2`:

- two operations are **constant** (`h` ignores both inputs);
- four are **unary** (exactly one of `b_1, b_2` is nonzero);
- two are **XOR-type** (`b_{12} = 0` and `b_1 = b_2 = 1`): `t_1 + t_2` and
  `t_1 + t_2 + 1`;
- eight are **AND-type** (`b_{12} = 1`): `(t_1 + alpha)(t_2 + beta) + gamma` for
  `alpha, beta, gamma in F_2`.

These counts are exhaustive: `2 + 4 + 2 + 8 = 16`, and they are re-derived by
brute force over all sixteen truth tables in `checks.py`. An operation is
**affine** when `b_{12} = 0`, i.e. the eight constant, unary and XOR-type
operations.

Two elementary facts are used repeatedly.

**Fact 1.1.** An AND-type gate `(t_1 + alpha)(t_2 + beta) + gamma` is constant
if one of its inputs is constant, and is a literal of the other input if
`t_1 = t_2`.

**Fact 1.2.** If two affine functions `l_1, l_2` of the same variables have
linearly dependent linear parts, then `(l_1 + alpha)(l_2 + beta) + gamma` is
constant or unary (hence computable with no gate, after absorbing the
negation).

## 2. The analysis device: fair semi-circuits

Gate elimination repeatedly substitutes expressions for variables. To keep the
bookkeeping exact we use the representation of Find, Golovnev, Hirsch and
Kulikov, in the "parameter" form of the source.

**Definition 2.1 (fair semi-circuit).** Let `u = (u_1, ..., u_m)` be variables,
called *parameters*. A fair semi-circuit on `u` consists of

(i) an **affine system** of `q >= 0` equations in the parameters and in `q`
    further variables `v = (v_1, ..., v_q)`, written `A v = N u + c` with
    `A in F_2^{q x q}`, `N in F_2^{q x m}`, `c in F_2^q`, such that each row of
    `(A | N)` has at most **three** nonzero entries and `A` is invertible
    (*fairness*: for every assignment of `u` there is exactly one solution `v`);

(ii) an **acyclic part**: an ordinary circuit in the sense of Section 1 whose
     gates read constants, parameters, auxiliary variables `v`, or earlier
     gates, listed in a topological order; there are no wires from the acyclic
     part back into the affine system;

(iii) an **output**, which is a constant or a literal of a parameter, an
     auxiliary variable, or a gate. No variable of an equation is an output.

Parameters, auxiliary variables and gates are the *internal values*. `s` denotes
the **size**: the number `q` of equations plus the number of gates in the
acyclic part. For the original circuit `C_0` the affine system is empty and the
acyclic part is `C_0`, so `s = |C_0|`.

Fairness makes every auxiliary variable an affine function of the parameters,
`v = A^{-1} N u + A^{-1} c`, so a semi-circuit computes a well-defined Boolean
function of `u`, and every internal value is an affine function of `u` or a
function of parameters computed by the acyclic part. Cyclic dependencies among
auxiliary variables are allowed; only invertibility of `A` matters.

**Lemma 2.2 (moving an affine gate into the affine system).** Let `J` be a gate
of the acyclic part labelled by an affine operation, whose two inputs are
constants, parameters or auxiliary variables. Replacing `J` by the equation
`J + b_1 t_1 + b_2 t_2 = c`, regarding `J` as a new auxiliary variable, and
rewriting every use of the gate `J` as a use of the variable `J`, yields a fair
semi-circuit computing the same function of the parameters, with the same size,
and with every equation still having at most three nonconstant terms.

*Proof.* The new equation has at most three nonconstant terms by construction.
The new auxiliary matrix is `[[A, 0], [*, 1]]` in block form, hence invertible.
For each assignment of the parameters the old solution extends uniquely with `J`
equal to the value the gate computed, so the function is unchanged and no gate
is created or destroyed: one gate becomes one equation. `QED`

Applying Lemma 2.2 in topological order moves every affine gate with **no
AND-type ancestor** into the affine system, and leaves every affine gate that
does have an AND-type ancestor in the acyclic part. This operation is used in
Section 9 and changes neither `s`, nor `m`, nor `q`, nor the degree of any
parameter (Section 5), hence does not increase the potential.

## 3. Parameters, pairs, and the four operations

During the elimination we keep an **injective** affine map

    Phi : F_2^m -> F_2^n,   Phi(u) = b + L u,   b in F_2^n, L in F_2^{n x m},

whose variables `u` are the *current parameters*, together with a set `P` of
**disjoint pairs** of parameters. A parameter in a pair is *paired*, its partner
is the other member, the rest are *unpaired*. Initially `Phi` is the identity,
`m = n`, `P` is empty, and the parameters are the input variables `X`.

For a set `P` of disjoint pairs put

    Z(P) = { u in F_2^m : p * r = 0 for every {p,r} in P }.

So `Z(P)` is the union of the `3^{|P|}` coordinate subspaces obtained by forcing,
in each pair, one of the three lines `p = 0`, `r = 0`, `p = r + 1`.

The **state** is the triple (current semi-circuit `C`, map `Phi`, set `P`), and
it is tracked by four numbers: the size `s`, the number `m` of parameters, the
number `q = |P|` of pairs, and

    a = m - q.

A parameter may remain in `Phi` after `C` stops reading it; `m` counts it.

**The four operations.** Below, `u - x` denotes all parameters except `x`, and
`w` is a fresh parameter not occurring in `Phi`.

- **(AR) affine reparameterisation.** For an unpaired parameter `x` and an
  affine map `phi`, substitute `x <- w + phi(u - x)`. This is a bijection of the
  parameter space: `m, q, a, s` are unchanged, and `im Phi` is unchanged.
- **(UF) fixing an unpaired parameter.** For an unpaired `x` and an affine map
  `psi`, impose `x = psi(u - x)` (the constant map is allowed). This removes one
  parameter: `m` decreases by one, `q` is unchanged, so `a` decreases by one.
- **(PF) fixing a paired parameter.** For a pair `{p,r}`, impose `p = 0`, or
  `r = 0`, or `p = r + 1` — the three lines inside `p*r = 0` — and delete the
  substituted parameter together with the pair. Both `m` and `q` decrease by
  one, so `a` is unchanged.
- **(PC) adding a pair.** For two distinct unpaired parameters `x, y`, add
  `{x,y}` to `P`. `Phi` is unchanged, but `Z(P)` gains the constraint
  `x*y = 0`. `m` is unchanged and `q` increases by one, so `a` decreases by one.

Each operation replaces `Phi` by `Phi o tau`, where `tau` expresses the old
parameters in the new ones. In (AR) `tau` is a bijection; in (UF) and (PF)
`tau` fixes one parameter as a function of the others, hence is injective; in
(PC) `tau` is the identity. In every case `Phi o tau` is again injective and of
the form `b + L u`, and `tau` maps `Z(P')` into `Z(P)`. Since each operation
turns `C` into `C'` with `C'(u) = C(tau(u))` for all `u in Z(P')`, the
**fundamental invariant**

    C(u) = C_0( Phi(u) )   for every u in Z(P)                              (*)

holds at every state. `(*)` is the only thing used at the end: the current
semi-circuit agrees with the original circuit on the affine image of `Z(P)`.

**Implementing (AR) without breaking the three-term bound.** A general affine
substitution for a parameter into `A v = N u + c` could turn one term into many
and violate the three-nonconstant-term condition. So (AR) is implemented by
**exchanging the roles of a parameter and an auxiliary variable**, changing only
`Phi`. Suppose the auxiliary variable `Q` has the solved expression
`Q = x + phi(u - x)` in the parameters, where `x` is unpaired and `phi` affine.
Substitute `x <- w + phi(u - x)` in `Phi`, rename `Q` as the new parameter `w`,
and regard `x` as an auxiliary variable. No equation changes; concretely the
column of `x` in `N` and the column of `Q` in `A` are exchanged, so every
equation keeps exactly the same terms.

*Invertibility is preserved.* Let `j` be the column index of `Q`, let `N_x` be
the column of `x` in `N`, and let `t = A^{-1} N_x`. Because `Q` depends on `x`,
the `j`-th coordinate of `t` is `1`. The new auxiliary matrix `A'` is `A` with
its `j`-th column replaced by `N_x`, i.e. `A' = A B` where `B` is the identity
with its `j`-th column replaced by `t`. Then `det B = t_j = 1`, so `A'` is
invertible. `QED`

**Rerouting.** When a gate or an auxiliary variable `z` becomes a constant or a
literal of another internal value, replace every occurrence of `z` by that
constant or literal and delete `z`. If `z` is an auxiliary variable occurring in
an equation with at most two nonconstant terms, say `z = U + c_0`, then
substituting `z = U + c_0` elsewhere is the same as adding one row of
`(A|N|c)` to another: restricted to `A` these are elementary row operations
(determinant unchanged), and afterwards the column of `z` in `A` has a single
nonzero entry, so deleting its row and column leaves an auxiliary matrix whose
determinant is `det A` up to sign, by expansion along that column. Hence
rerouting an auxiliary variable, and deleting an auxiliary variable that occurs
in a single equation together with that equation, both preserve fairness of the
system; and replacing a term by a constant, or by a literal of a single other
term, never increases the number of nonconstant terms in an equation.

## 4. Degrees, weights and the potential

**Parents, children, degree.** The *parents* of a gate are its inputs (constants
or internal values); a gate is a *child* of each internal value it reads, and an
equation is a child of each internal value occurring in it with a nonzero
coefficient. A gate reading the same value twice is one child of that value. The
*degree* `deg(t)` of an internal value `t` is its number of distinct children.
A parameter is *unused* if `deg = 0` and *isolated* if `deg = 1`. Degree is
syntactic: it counts direct occurrences only, so a solved dependence of an
auxiliary variable on a parameter does not contribute to that parameter's degree.

**Weights and the potential.**

    wt(x) = 0 if deg(x) = 0,    wt(x) = 6 if deg(x) >= 1      (unpaired x)

    W = s + 4q + sum over unpaired parameters x of wt(x).

Paired parameters contribute only through the term `4q`; the potential is
nonnegative, and `W = s + 4q + 6 I` where `I` is the number of used unpaired
parameters. For any tracked quantity `Y` we write `dY = Y' - Y` for the change
across one step (possibly a group of operations).

**Lemma 4.1 (rerouting does not create uses).** Suppose an object (a gate or an
equation) is replaced everywhere by a constant, or by a literal of an internal
value that already occurs in that object, and then deleted. Then every parameter
that is unused before the replacement is unused after it, and no unpaired
parameter's weight increases.

*Proof.* Every child occurrence of a parameter `z` is either unchanged or
destroyed. New occurrences arise only by copying occurrences that were already
present in the deleted object; if the deleted object mentioned `z`, then
`deg(z) >= 1` before the operation, so `wt(z) = 6` before and after, which is
the maximum. If `z` was unused, it is mentioned by nothing deleted and gains no
new child, so it stays unused. `QED`

Consequently, in every accounting below, only three effects can raise `W`: a
pair being destroyed (each pair contributes `-4` to `dW` when removed, so
destroying a pair without recreating one *lowers* `W`), a parameter that becomes
unpaired (it can gain at most `+6`), and the `+4` of a newly created pair. Every
other effect is a decrease coming from `ds <= 0` and from removed pairs.

## 5. The subspace lemma

**Lemma 5.1.** Suppose `a > d` and the output of the current semi-circuit `C`
is either

(i) an affine function of the parameters, or
(ii) a literal of an AND-type gate whose two parents are affine functions of the
    parameters.

Then, given those affine functions, one can construct in polynomial time an
affine subspace `U subseteq Z(P)` with `dim U >= a - 1 >= d`, and a constant
`b`, such that `C(u) = b` for every `u in U`.

*Proof.* Write an affine function as `l(u) = sum_j lambda_j u_j + kappa`. Choose
one parameter from each pair of `P` and fix it to `0`; this gives an affine
subspace `U_0 subseteq Z(P)` of dimension `m - q = a`, represented by the point
`0` and the standard basis vectors `e_j` of the unfixed parameters. (Each pair
loses exactly one of its two parameters, so `3^{|P|}` becomes one line per pair
and the dimension drops from `m` to `m - q`.)

*Step 1 (linear case).* Let `l` be an affine function. Since `U_0` is a linear
subspace, `l` is constant on `U_0` iff `lambda_j = 0` for all basis indices `j`
of `U_0`, in which case the value is `kappa`. Otherwise pick `j_0` with
`lambda_{j_0} = 1`; then for `u in U_0`,
`l(u) = 0` iff `u_{j_0} = kappa + sum_{j != j_0} lambda_j u_j`, so
`{u in U_0 : l(u) = 0}` is the affine subspace with point `kappa e_{j_0}` and
basis `{ e_j + lambda_j e_{j_0} : j != j_0 }`, of dimension `a - 1`.

*Step 2 (AND case).* If the output is a literal of an AND-type gate with affine
parents, it equals `(l_1 + alpha)(l_2 + beta) + gamma`. If `l_1 + alpha` and
`l_2 + beta` are both constant on `U_0`, with values `c_1, c_2`, take `U = U_0`
and `b = c_1 c_2 + gamma`. Otherwise one of them, say `l_1 + alpha`, is
nonconstant on `U_0`; by Step 1 the set `U` of points of `U_0` where
`l_1 + alpha = 0` is an affine subspace of dimension `a - 1`, and the output
equals `gamma` there: take `b = gamma`.

In both cases `C` is constant equal to `b` on `U`, and `dim U >= a - 1 >= d`
because `a > d` and dimensions are integers. The construction compares at most
two coefficient vectors against the basis of `U_0`, so it is polynomial time.
`QED`

Because `Phi` is injective, `S = Phi(U)` then has `dim S = dim U >= d`, and by
`(*)` the original circuit satisfies `C_0(x) = b` for every `x in S`.

## 6. Deletions

Each deletion either removes a gate of the acyclic part or an equation of the
affine system, or removes an auxiliary variable together with its equation. In
all four cases `dm = dq = da = 0`, `ds <= -1`, and no unpaired parameter
increases its weight, hence `dW <= -1`.

**D1 (delete a gate depending on at most one input).** If a gate `J` of the
acyclic part computes a constant or a literal of one of its inputs — true for
every constant or unary gate, every gate with a constant input, every gate whose
two inputs coincide, and (Fact 1.2) every AND-type gate whose inputs are
parameters or auxiliary variables with dependent linear parts — reroute every
use of `J` to that constant or literal and delete `J`. We use D1 in exactly this
form, including the case where a gate *becomes* constant on all of `Z(P)` at
some state: it is then deleted by D1, rerouted to the constant it computes.
`ds = -1`, `dW <= -1` by Lemma 4.1.

**D2 (delete an unused gate).** A gate with no child that is not the designated
output can be deleted with the function unchanged. `ds = -1`, `dW <= -1`.

**D3 (delete an equation with at most two nonconstant terms).** Such an equation
must contain an auxiliary variable, else the corresponding row of the auxiliary
matrix would be zero, contradicting invertibility. Write it as `v = U + c` with
`v` auxiliary and `U` a parameter, an auxiliary variable, or `0`. Substitute and
delete `v` with the equation (Section 3, rerouting); `dm = dq = da = 0`,
`ds = -1`, `dW <= -1`.

**D4 (unit parameter deletion).** Let `x` be an unpaired, isolated parameter
whose unique child is an equation `E`. Since the column `N_x` is supported on the
row of `E`, `A^{-1}N_x != 0`; choose an auxiliary variable `Q` whose coordinate
in `A^{-1}N_x` is `1`, i.e. a `Q` that depends on `x`; then `Q = x + phi(u - x)`.
Apply (AR) with this `Q`, exchanging the roles of `x` and `Q`: no equation
changes, and `x` (now auxiliary) still occurs only in `E`. Delete `x` and `E`.
`dm = dq = da = 0` (one parameter replaces another), `ds = -1`; the old
parameter `x` had weight `6` and is replaced by a parameter of weight at most
`6`, `dW <= -1`.

**Affine reparameterisation (AR).** `dm = dq = da = ds = 0`, and
`dW = wt'(w) - wt(x) = 6 - 6 = 0`: `x` occurs in some equation (otherwise
`A^{-1}N_x = 0` and `Q` could not depend on it), and the new parameter `w = Q`
occurs in some equation (the column of `Q` in the invertible `A` is nonzero), so
both have degree at least one. Hence (AR) is free: it never increases `W` and
never changes `a`.

## 7. The iteration lemma

**Lemma 7.1.** Suppose `a > d`. In deterministic polynomial time the algorithm
below either constructs an affine subspace `S subseteq F_2^n` of dimension at
least `d` and a constant `b` with `C_0(x) = b` for all `x in S`, or transforms
the state so that `m' + s' < m + s`, `da in {0,-1}` and `dW <= 10 da`.

**The algorithm.** Call a state *terminal* when the output is a constant, or a literal of a
parameter or auxiliary variable, or a literal of an AND-type gate whose two parents are
constants, parameters or auxiliary variables — exactly the two forms of Lemma 5.1. The
algorithm is a loop, and the terminal test is re-run after **every** primitive operation,
not only at the start of an iteration:

0. *Terminal test.* If the state is terminal, apply Lemma 5.1 and return `S = Phi(U)`, `b`.
1. *Normalise.* Apply one of: moving an affine gate without an AND-type ancestor into the
   affine system (Lemma 2.2), or a deletion D1--D4, subject to guards G1--G3 below. Then
   return to step 0. When no operation applies, go to step 2.
2. *Preliminary reductions.* Apply the first applicable of Cases 1, 2, 3, then return to
   step 0.
3. *Minimal AND.* Otherwise find the first gate `G` of the acyclic part in topological
   order; turn its inputs into parameters by at most two (AR) applications, returning to
   step 0 after each; then apply the first applicable of Cases 4, 5.

**Guards, and why the control flow is sound.** The failure modes to exclude are that
normalisation exposes a terminal state but the algorithm continues, that a deletion removes
the designated output, and that the acyclic part is emptied so that step 3 has no gate to
find. The three guards are:

- **G1 (the output is never deleted or rerouted).** Deletions apply only to objects that are
  not the designated output. If the output gate becomes constant or a literal, that is
  precisely a terminal state, detected at step 0 *before* any deletion is applied, and the
  algorithm returns instead of deleting. So D1 is never applied to the output gate.
- **G2 (D4 never deletes the output).** D4 deletes an unpaired isolated parameter. Such a
  parameter is the output only if the output is a literal of a parameter, which is a
  terminal state; hence, whenever D4 fires, the state is not terminal, so the output is not
  a parameter and not an auxiliary variable, and no parameter deleted by D4 is the output.
  The same argument covers the auxiliary variable deleted by D3.
- **G3 (step 3 always has a gate).** If the acyclic part becomes empty, then the output
  gate itself has been moved into the affine system, so the output is a literal of a
  parameter or auxiliary variable — a terminal state, caught at step 0. So on reaching
  step 3 the acyclic part is nonempty, and its first gate in topological order exists. It is
  AND-type with independent inputs, because all affine gates with no AND-type ancestor have
  been moved into the affine system (step 1) and D1 does not apply.

The terminal test and the guards change none of the tracked quantities: no parameter is
fixed or paired, no pair is removed, no gate or equation is deleted and `Phi` is not
altered. Hence they contribute `dm = dq = da = ds = dW = 0`, and the accounting of
Cases 1--5 and of the deletions below is unaffected by where the test is performed.

**Case 1: a paired parameter has degree zero or at least two.** Let `p` be paired
with partner `r`.

- If `deg(p) = 0`: apply (PF) fixing `r` to `0`. The pair is removed, `p` becomes
  unpaired with degree `0`, so it contributes nothing. `dm = -1`, `dq = -1`,
  `da = 0`, `ds = 0`, and `dW <= -4`.
- If `deg(p) >= 2`: apply (PF) fixing `p` to `0`, then delete two children of
  `p` by D1 or D3. `dm = -1`, `dq = -1`, `da = 0`, `ds <= -2`. Only the newly
  unpaired `r` can gain weight, by at most `6`. Hence
  `dW <= -4 - 2 + 6 = 0`.

If Case 1 does not apply, every paired parameter is isolated.

**Case 2: two partners have the same child.** Let `p` be paired with `r`, both
isolated, with a common child `H`; by Case 1's absence, `H` is the unique child
of both. If `H` is a gate, its two inputs are `p` and `r`; the two-input function
`H` takes equal values at two of the three points `(0,0), (0,1), (1,0)` of the
line `p*r = 0` (pigeonhole on three points and two values), and the line through
those two points is one of `p = 0`, `r = 0`, `p = r + 1`; restrict to it with
(PF). Then `H` is constant on `Z(P')` and is deleted by D1. If `H` is an
equation, then it has three nonconstant terms (D3 does not apply) and at least
one auxiliary variable, so it reads `v + p + r = c`; impose `p = r + 1` with
(PF), which turns it into `v = c + 1`, and delete it by D3. In both cases the
surviving parameter of the pair becomes unused (its only child was `H`), and
only constants are rerouted. `dm = -1`, `dq = -1`, `ds = -1`, `da = 0`, so
`dW <= -4 - 1 = -5`.

**Case 3: a paired parameter occurs in the affine system.** Let `p` be paired
with `r`, with unique children `H_p` and `H_r` (Case 1 absent), both equations
(`p` occurs in the affine system), and `H_p != H_r` (Case 2 absent). Apply (PF)
fixing `r` to `0` and delete `H_r` by D1 or D3; `p` does not occur in `H_r`, so
this creates no new child of `p`. Now `p` is unpaired and isolated with unique
child `H_p`. If `H_p` has at most two nonconstant terms, delete it by D3. If it
has three, apply D4 to `p` and `H_p`. `dm = -1`, `dq = -1`, `ds <= -2`,
`da = 0`; the only possible increase is `p`'s weight, at most `6`, so
`dW <= -4 - 2 + 6 = 0`.

If Cases 1--3 do not apply, every paired parameter is isolated, occurs only in
the acyclic part, and has a child different from its partner's child.

**Making the minimal AND explicit.** Suppose Cases 1--3 do not apply. The acyclic
part has at least one gate (otherwise the output is an affine function of the
parameters and step 1 would have returned). Let `G` be its first gate in
topological order. All affine gates reading only parameters and auxiliary
variables have been moved into the affine system, so `G` is AND-type; and D1
does not apply, so both inputs of `G` are parameters or auxiliary variables with
linearly **independent** linear parts.

If an input of `G` is an auxiliary variable `Q`, it is nonconstant (its linear
part is independent of the other input's) and depends on some parameter `x`;
that `x` must occur in some equation, else the column of `x` in `N`, hence in
`A^{-1}N`, would vanish, contradicting the dependence. By the absence of Case 3,
every parameter occurring in the affine system is unpaired, so `x` is unpaired
and `deg(x) >= 1`. Apply (AR), exchanging `Q` and `x`: the new parameter `w`
replacing `Q` has degree at least two (it occurs in an equation, since the
column of `Q` in the invertible `A` is nonzero, and it is a parent of `G`), and
`x` has weight `6`, so `dW = 0` and `dm = ds = dq = da = 0`. Re-run steps 1--3
and, if none applies, continue. Repeating at most twice (the linear parts of the
two inputs stay independent under invertible affine changes of parameters) makes
both inputs parameters; call them `x` and `y`, and write
`G = (x + alpha)(y + beta) + gamma`.

If the output is a literal of `G`, step 1 applies (Lemma 5.1). Otherwise `G` is
not the output and, D2 being inapplicable, it has a child; every child of `G` is
a gate of the acyclic part, since no wire returns to the affine system.

**Case 4: both `x` and `y` unpaired.** First absorb the constants: if
`alpha = 1`, substitute `x <- w + 1` with a fresh `w` (an (AR) with `phi = 1`),
so that every occurrence of `x` becomes `w + 1` and the negation is absorbed into
the constant terms of the equations and the operations of the gates reading `x`;
`dm = ds = dq = da = dW = 0` and `w` has the same degree as `x`. Do the same for
`y` if `beta = 1`. Now `G = x*y + gamma`.

Add the pair `{x,y}` by (PC). Because every `u in Z(P')` satisfies `x*y = 0`,
`G` equals the constant `gamma` on `Z(P')`; reroute `G` to `gamma` and delete it
by D1. Its child (there is at least one) now has a constant input, so delete one
child by D1 as well. Before the step both `x` and `y` had weight `6` (they are
parents of `G`); now they are paired with no individual weight, and the new pair
contributes `+4`. No other unpaired parameter increases (Lemma 4.1; rerouting
goes to a constant or to a literal of another input). Hence

    dm = 0,   ds <= -2,   dq = +1,   da = -1,
    dW <= -2 + 4 - 6 - 6 = -10 = 10 da.

**Case 5: at least one of `x`, `y` is paired.** By Cases 1--3 every paired input
of `G` has `G` as its unique child, and its partner has a unique child, a gate of
the acyclic part different from `G`. Let `k in {1,2}` be the number of paired
inputs.

*Subcase `k = 1`.* Say `y` is paired with partner `r`, and `x` is unpaired. Apply
(PF) fixing `r` to `0`, which removes the pair; then (UF) fixing `x` to the
constant `alpha`, so that `G = gamma` on the new `Z(P')`. Delete `G` by D1 and
reroute it to `gamma`; `y` becomes unused. `H_r`, the unique child of `r`, now
has a constant input, so delete it by D1. Weight changes: `x` loses `6`; `y` had
no individual weight (paired) and has weight `0` afterwards (unused); `r` is
removed. No initially unpaired parameter gains by Lemma 4.1. Hence

    dm = -2,   ds <= -2,   dq = -1,   da = -1,   dW <= -2 - 4 - 6 = -12 <= 10 da.

*Subcase `k = 2`.* Let `r_x, r_y` be the partners of `x, y`; the pairs
`{x,r_x}` and `{y,r_y}` are distinct by the absence of Case 2. Apply (PF) twice,
fixing `r_x` and `r_y` to `0`, removing both pairs; then (UF) fixing `x` to
`alpha`, so that `G = gamma`. Delete `G` by D1 and reroute; `y` becomes unused.
Delete `H_{r_x}` (the unique child of `r_x`) by D1. Weight changes: `x` and `y`
were paired, hence had no individual weight, and afterwards `x` is fixed (not a
parameter) and `y` is unused, so both contribute `0`. Hence

    dm = -3,   ds <= -2,   dq = -2,   da = -1,   dW <= -2 - 8 = -10 = 10 da.

**Proof of Lemma 7.1.** If step 0 returns, the conclusion is Lemma 5.1 together
with `(*)` and injectivity of `Phi`. By G1--G3 the loop is well defined: every deletion
acts on an object that is not the output, the states in which the acyclic part is empty or
the output is affine are terminal and are caught at step 0, and each pass through steps 1--3
either returns or applies exactly one of the moves analysed below. Otherwise, the normalisation step performs
moves of affine gates into the affine system (`dm = ds = da = 0`, `dW <= 0`) and
deletions (`da = 0`, `dW <= -1`), and then exactly one of Cases 1--5 fires with
`da in {0,-1}` and `dW <= 10 da`; Case 4's preparatory (AR) steps and the
minimal-AND restructuring have `dW = 0` and `da = 0`. Summing the pieces gives
the stated bounds for the whole iteration. Termination of the iteration count:
`ds < 0` in Cases 4 and 5 and in every deletion, and `dm < 0` in Cases 1, 2, 3,
5 with no step increasing `s`; so `m + s` strictly decreases. Each iteration
costs at most `s` moves into the affine system, at most two (AR) exchanges,
Gaussian elimination over `F_2` to express auxiliary variables as affine
functions of the parameters (used for dependence tests and (AR)), and constant
size truth-table tests, all polynomial in `n`. `QED`

## 8. The main theorem

**Theorem 10.1 (restated as a bound).** Let `f : F_2^n -> F_2` be an affine
disperser for dimension `d` with `d < n`. Then every B2 circuit `C_0` on `n`
inputs computing `f` satisfies `|C_0| >= 4n - 10d`.

*Proof.* Suppose `|C_0| < 4n - 10d`. Start the algorithm of Lemma 7.1 from the
initial state: `Phi` the identity, `P` empty, the affine system empty, the
acyclic part `C_0`. Then `m = a = n`, `q = 0`, `s = |C_0|`, and

    W_0 = |C_0| + 6 * #{i : deg(X_i) >= 1} <= |C_0| + 6n < (4n - 10d) + 6n = 10(n-d).

Run Lemma 7.1 repeatedly, and let `a_t`, `W_t` be the values of `a` and `W` after
`t` returned iterations. Summing `dW <= 10 da` over the `t` iterations,

    0 <= W_t <= W_0 + 10 (a_t - n) < 10(n - d) + 10(a_t - n) = 10(a_t - d),

so `a_t > d` after every returned iteration: the hypothesis `a > d` of Lemma 7.1
is again satisfied, and the algorithm does not stop for `a <= d`. Since each
returned iteration strictly decreases `m + s`, which is a nonnegative integer
starting at `n + |C_0|`, after at most `n + |C_0|` iterations Lemma 7.1 must
return an affine subspace `S subseteq F_2^n` with `dim S >= d` and a constant `b`
such that `C_0(x) = b` for every `x in S`.

But `f` is an affine disperser for dimension `d`, so `f|_S` is nonconstant, and
there is `x in S` with `f(x) != b = C_0(x)`. So `C_0` does not compute `f`,
contradicting the hypothesis. Hence `|C_0| >= 4n - 10d`. `QED`

**Remark 10.2 (what the argument does and does not give).** The proof is a
*refutation* argument: it shows that a too-small circuit cannot be correct. The
source's algorithmic form additionally finds such an `x` in time
`poly(n, 2^d)` by evaluating `f` on a `d`-dimensional subspace; for the
`d = O(n^{4/5})` family used here that is exponential, and the constructive form
is not claimed. Only the lower bound is claimed. The coefficient `4` is exactly
the coefficient that this potential analysis yields with weights `(6,4)`; the
source's Remark 4.3 states that no other weights `alpha, beta` in
`V = s + beta q + alpha I` improve the leading coefficient, and that remark is
attributed, not used, here.

## 9. Accounting summary

| step | `dm` | `dq` | `da` | `ds` | bound on `dW` |
| --- | --- | --- | --- | --- | --- |
| D1, D2, D3 | 0 | 0 | 0 | -1 | `<= -1` |
| D4 | 0 | 0 | 0 | -1 | `<= -1` |
| (AR) | 0 | 0 | 0 | 0 | `= 0` |
| move a gate into the system | 0 | 0 | 0 | 0 | `<= 0` |
| Case 1, `deg(p)=0` | -1 | -1 | 0 | 0 | `<= -4` |
| Case 1, `deg(p)>=2` | -1 | -1 | 0 | `<= -2` | `<= 0` |
| Case 2 | -1 | -1 | 0 | -1 | `<= -5` |
| Case 3 | -1 | -1 | 0 | `<= -2` | `<= 0` |
| Case 4 | 0 | +1 | -1 | `<= -2` | `<= -10` |
| Case 5, `k=1` | -2 | -1 | -1 | `<= -2` | `<= -12` |
| Case 5, `k=2` | -3 | -2 | -1 | `<= -2` | `<= -10` |

Every row satisfies `dW <= 10 da`; `checks.py` recomputes each row from the
stated weight changes with exact integers and confirms the table, and confirms
that the only rows with `da = -1` are exactly those needing `dW <= -10`.

## 10. What is machine-checked, and what is not

`checks.py` verifies mechanically: the classification of the sixteen operations;
each row of the accounting table (as integer arithmetic on the stated
`ds, dq, dm` and weight changes); the identity `a = m - q` along each row; the
summation inequality of Theorem 10.1 on a grid of rational values of
`|C_0|, n, d`; and the parameter inequalities of Lemma 5.1 of `family.md` for a
large range of `n` by exact integer arithmetic. It does not, and cannot, verify
the case analysis: that each case is applicable when its predecessors are not,
that the geometric claims (Case 2's pigeonhole line, Lemma 5.1's dimension
count, the independence of the two inputs of the minimal AND) hold, or that the
invariant `(*)` is preserved. Those are the mathematical content of Sections 2--8
and are argued there in prose, as they must be.

## 11. Where each rule requirement is met

- *Fixed family, total, uniformly defined*: `family.md` Sections 2--3, all
  lengths including `n = 0`.
- *NP with fixed verifier and witness bound*: `family.md` Corollary 3.2.
- *Inputs are `x` alone, no witness bits in the circuit*: the family is a
  function of `x`; the witness is empty.
- *Full B2 model, no transfer from a restricted model*: Theorem 10.1 is stated
  and proved for B2 circuits directly; no monotone, formula, bounded-depth or
  XOR-free statement is used. The semi-circuits of Section 2 are an analysis
  device that *contains* the original circuit as its initial state and is never
  used as a stronger hypothesis.
- *Gate counting*: gates only, each retained gate once (Section 1); constants
  and negations are free, absorbed as the model allows.
- *All-length quantifiers and explicit `o(n)` loss*: `proof.md` Section 5.

## 12. Attribution

The gate-elimination argument, its five cases and its potential are due to
Mishra and Santhanam, ECCC TR26-260, revision 1 of 3 November 2026,
"An Explicit `4n - O(log n)` Size Lower Bound against General Boolean Circuits"
(SHA256 `a123943ed05c9569385d3b9d98c695630a110559f33a4c990ec44e2a22b4d921`);
it builds on Find--Golovnev--Hirsch--Kulikov (FOCS 2016) and Li--Yang
(STOC 2022). The exposition here is rewritten for this submission, follows the
source's structure, and is stated in the challenge's model conventions. No new
mathematical result is claimed in this file; the contribution of the submission
is the assembly in `proof.md`: applying this theorem to the subspace-polynomial
family, with the complete eligibility, model, bit-accounting and quantifier
argument, and with the arithmetic of every case re-verified.
