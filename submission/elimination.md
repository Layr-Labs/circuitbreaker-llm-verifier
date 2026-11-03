# Gate elimination for an affine disperser: every full-B2 circuit has 4n - 10d gates

**Class theorem (T).** Let `d` be a nonnegative integer and let
`f : F_2^n -> F_2` be nonconstant on every affine subspace of `F_2^n` of
dimension `d`. Then every circuit for `f` in the full binary basis `B2` has at
least

    4n - 10d

gates.

This file is a complete, self-contained proof of (T). It covers arbitrary
sharing, repeated pins, constant and unary gates, degenerate tables, a single
terminal output, and the full unrestricted-depth, unrestricted-fan-out model.
No finite experiment is a premise of any step.

**Attribution and scope.** The gate-elimination argument is the parameter method
of Mishra and Santhanam, *An Explicit 4n - O(log(n)) Size Lower Bound against
General Boolean Circuits*, ECCC TR26-260 **Revision 1** (3 November 2026), and
specifically its Theorem 4.1 and Lemma 4.2, together with Sections 2.2-2.6 and
3.1-3.2. The pinned source conventions of `rules.md` (Li-Yang ECCC TR21-023,
Sections 1.1, 2.1, 2.4 footnotes 3-4 and 3.3) are the model in which (T) is
stated. Section 8 below records the exact algebraic steps that are reproduced
here rather than cited, and Section 9 records what is *not* claimed.

Everything in Sections 1-7 is proved in this file. Citations name a version, a
theorem or a section, the relevant statement, and how its hypotheses are met;
nothing here rests on an unreproduced external premise beyond the model
conventions of `rules.md`, which the organizer's own rules fix.

## 1. Model, notation, and the class of functions

### 1.1 The circuit model

A `B2` circuit is an acyclic network of gates over `F_2`. Each gate has two
input pins and one output, and is labelled by one of the sixteen binary
Boolean functions `h : F_2^2 -> F_2`. Depth and fan-out are unrestricted. The
output is the value of one gate or one input variable. Inputs are the `n`
variables `X_1,...,X_n` and are not gates; wires and the output selection cost
nothing. Every retained gate is counted once, including constant and unary
gates.

We identify `F_2 = {0,1}` with addition `+` (XOR) and multiplication `*` (AND).
Every `h : F_2^2 -> F_2` has a unique multilinear expression

    h(t1,t2) = c + b1 t1 + b2 t2 + b12 t1 t2 ,      c,b1,b2,b12 in F_2 .

A gate is **constant** if `b1 = b2 = 0`; **unary** if exactly one of `b1,b2` is
`1`; otherwise it is **XOR-type** (`b12 = 0`, so `h = t1 + t2 + c`) or
**AND-type** (`b12 = 1`, so `h = (t1 + alpha)(t2 + beta) + gamma` for unique
`alpha,beta,gamma in F_2`). Of the sixteen tables, two are constant, four are
unary, two are XOR-type and eight are AND-type.

A gate is **affine** if its operation is `b1 t1 + b2 t2 + c`, i.e. `b12 = 0`
and at most one of `b1,b2` is `1`; such a gate is constant or unary or
XOR-type. A **literal** of a value `t` is `t` or `t + 1`. Negations on wires
are free, because they are absorbed into the operation of the gate that reads
the wire; a negated projection that is itself a retained gate costs one gate.

The affine-disperser hypothesis of (T) says: `f` is nonconstant on every affine
subspace of dimension exactly `d`. (An affine subspace of dimension `> d`
contains one of dimension `d`; so `f` is then nonconstant on it as well. We
use this only in the direction "nonconstant on a `d`-dimensional subspace".)

### 1.2 Internal values, parents, children, degree

The **internal values** of a circuit are its gates, its input variables, and
any auxiliary variables introduced later. The **parents** of a gate are its two
input values; the **parents** of an equation are the internal values appearing
in it with a nonzero coefficient. The **children** of an internal value `t` are
the equations that mention `t` and the gates that read `t`; a gate that reads
the same value twice contributes only one child. The **degree** `deg(t)` is the
number of distinct children of `t`. A value is **unused** if `deg(t) = 0` and
**isolated** if `deg(t) = 1`.

Degree is *syntactic*: it counts direct occurrences only.

**Rerouting.** When a gate or an auxiliary variable becomes a constant or a
literal of another value, every occurrence of it is replaced by that constant
or literal. Rerouting can increase the degree of a value that was a parent of
the replaced object. We use one elementary fact repeatedly:

> **(Rerouting lemma.)** In each deletion described in Section 3, whenever we
> reroute to a literal of a parameter, that parameter already occurs in the
> deleted object. Consequently a parameter that is unused before a deletion is
> unused after it.

*Proof.* A deleted gate or equation has parents among the internal values that
occur in it. Rerouting to a literal of `z` is possible only if `z` is a parent
of the deleted object, so `z` already had an occurrence there and was already
used. A constant rerouting adds no occurrence at all. No other value gains an
occurrence, since only the deleted object's occurrences are affected. QED

The same conclusion holds for the restrictions of Section 2: fixing a parameter
to a constant, or removing a parameter, replaces the removed parameter's
occurrences by constants, and any value that gains an occurrence was a parent
of the removed object's children and hence was already used.

## 2. States: fair semi-circuits, the map Phi, and pairs

Gate elimination is carried out on a slightly richer object than a circuit. A
**state** consists of

* parameters `u = (u_1,...,u_m)` over `F_2`, together with an injective affine
  map

      Phi : F_2^m -> F_2^n ,     Phi(u) = b + L u ,                  (1)

  with `b in F_2^n` and `L in F_2^{n x m}`. `Phi` records the affine changes
  made so far and is injective on the *whole* parameter space;
* a set `P` of **disjoint pairs** of parameters; a parameter in a pair is
  **paired** and the other member is its **partner**; all others are
  **unpaired**. Write `q = |P|` and

      Z(P) := { u in F_2^m : p r = 0 for every {p,r} in P }.            (2)

  So `Z(P)` is the set of parameter assignments in which at least one member of
  each pair is `0`. Each pair therefore contributes a **three-point** set
  `{00,01,10}` of the four corners of `F_2^2`;
* a **fair semi-circuit** `C` reading the parameters:

  * an **affine system** `A v = N u + c` over `F_2`, where `v = (v_1,...,v_r)`
    are **auxiliary variables**, `A in F_2^{r x r}`, `N in F_2^{r x m}`,
    `c in F_2^r`. Every row of `(A | N)` has at most three nonzero entries
    ("each equation has at most three nonconstant terms"), and `A` is
    invertible ("fairness": the system has a unique solution for every `u`).
    The auxiliary dependencies may be cyclic; invertibility is what gives
    uniqueness.
  * an **acyclic part**: gates in a topological order, each reading a constant,
    a parameter, an auxiliary variable, or an earlier gate. No wire goes back
    from the acyclic part into the affine system.
  * an **output**: a constant or a literal of a parameter, an auxiliary
    variable, or a gate.

The output is a *designation* rather than a gate, so it is not a child of the
value it designates; in particular a parameter designated by the output may
still be "isolated" in the sense of Section 1.2. We therefore require, as an
explicit part of the procedure, that **the output remain well defined at
every state**. Section 5, Step 1a gives the safeguard that enforces this.

Write `s` for the size of the fair semi-circuit, i.e. the number `r` of
equations plus the number of gates in the acyclic part.

**Invariant.** At every state,

    C(u) = C_0(Phi(u))   for all u in Z(P),                       (3)

where `C_0` is the original circuit. Initially `Phi` is the identity, `P` is
empty, the affine system is empty, the acyclic part is `C_0`, and
`m = a = n`, `q = 0`, `s = |C_0|`.

**Dimension.** Put `a = m - q`.

> **(Domain lemma.)** If `a >= d` then `f(Phi(u))` is nonconstant over
> `Z(P)`.

*Proof.* Fix, for each pair `{p,r} in P`, the member `p` to the value `0`, and
leave all unpaired parameters free. The resulting set of parameter assignments
is an affine subspace `U_0` of dimension `m - q = a` inside `Z(P)`. Since
`Phi` is affine and injective, `Phi(U_0)` is an affine subspace of `F_2^n` of
dimension `a`. If `a >= d` it contains an affine subspace `T` of dimension
exactly `d`, and `f` is nonconstant on `T` by the hypothesis of (T); hence `f`
is nonconstant on `Phi(U_0)`, hence on `Z(P)`. QED

### 2.1 The four operations

All four operations are implemented by replacing `Phi` with `Phi o tau` for an
injective map `tau` of parameter spaces; `tau` maps `Z(P')` into `Z(P)`, so
the invariant (3) is preserved.

* **(AR) Affine reparameterisation.** For an unpaired parameter `x` and an
  affine map `phi` of the other parameters, substitute `x <- w + phi(u - x)`
  with `w` a fresh parameter. This is a bijection, so `im Phi` is unchanged,
  and `m`, `q`, `a` are unchanged. It is implemented by exchanging the roles
  of a parameter and an auxiliary variable, as in Section 3.2; the dense
  affine map `phi` is *not* written into the computational representation.
* **(UF) Fixing an unpaired parameter.** For an unpaired parameter `x` and an
  affine map `psi` of the other parameters, impose `x = psi(u - x)`. `psi` may
  be constant. `m` decreases by one, `q` is unchanged, so `a` decreases by one.
* **(PF) Fixing a paired parameter.** For a pair `{p,r} in P`, substitute
  `p <- 0`, `r <- 0`, or `p <- r + 1`, and remove the substituted parameter
  together with the pair. These are exactly the three lines in `F_2^2`
  contained in `p r = 0`, so the restriction stays inside `Z(P)`. `m` and `q`
  each decrease by one, so `a` is unchanged.
* **(PC) Adding a pair.** For two distinct unpaired parameters `x,y`, add the
  pair `{x,y}` to `P`. `Phi` is unchanged; `Z(P)` gains the constraint
  `x y = 0`. `m` is unchanged and `q` increases by one, so `a` decreases by one.

In (UF) and (PF) the implementation substitutes only constants or literals into
the computational representation, so no gate and no longer equation is created.

### 2.2 The potential

Define the weight of an unpaired parameter by

    wt(x) = 0   if deg(x) = 0,
             6   if deg(x) >= 1,                                  (4)

and the potential by

    W = s + 4q + sum_{x unpaired} wt(x).                          (5)

Equivalently, with `I` the number of used unpaired parameters,

    W = s + 4q + 6 I .                                           (5')

Paired parameters contribute only through the `4q` term. `W >= 0` always.

Initially `q = 0`, every parameter is unpaired and has weight at most `6`, so

    W_0 <= |C_0| + 6n .                                          (6)

## 3. Deletions and affine reparameterisation

Each of the four deletions below has `Delta m = Delta q = Delta a = 0`,
`Delta s = -1`, and `Delta W <= -1`.

**(D1) Deleting a gate that depends on at most one input.** Suppose a gate
`J` in the acyclic part computes a constant or a literal of one of its inputs.
This covers every constant gate, every unary gate, every gate with a constant
input, every gate whose two inputs are identical, and every AND-type gate
whose two inputs are parameters or auxiliary variables with dependent linear
parts: two dependent affine forms are constant or equal up to a constant, so
the gate's operation depends on at most one effective input. We reroute every
use of `J` to that constant or literal and delete `J`. We also delete by (D1)
any gate that has become constant on all of `Z(P)`.

By the rerouting lemma, no unpaired parameter increases its weight, so
`Delta W <= -1`.

**(D2) Deleting an unused gate.** A gate with no children that is not the
designated output can be deleted. Its parents only lose occurrences, so again
no unpaired parameter increases its weight and `Delta W <= -1`.

**(D3) Deleting an equation with at most two nonconstant terms.** Let `E` be
such an equation. It must contain an auxiliary variable, since otherwise the
corresponding row of `A` would be zero, contradicting invertibility of `A`. So
`E` reads `v = U + c` for an auxiliary variable `v` and `U` a parameter, an
auxiliary variable, or `0`. Substitute `v = U + c` for all occurrences of `v`
and delete `v` together with `E`.

The two invariants are preserved. *Term counts:* replacing a term by a
constant, or by a literal of a single other term, cannot increase the number
of nonconstant terms of any equation. *Invertibility:* substituting `v` in
another equation is the same as adding the row of `(A | N | c)` for `E` to that
row, which is an elementary row operation and preserves `det A`. After the
substitution the column of `v` in `A` has a single nonzero entry, in the row
of `E`; expanding the determinant along that column shows that deleting that
row and column leaves a matrix with the same determinant up to sign, so the
remaining auxiliary matrix is invertible. No unpaired parameter increases its
weight, so `Delta W <= -1`.

**(D4) Unit parameter deletion.** Let `x` be an unpaired, isolated parameter
whose unique child is an equation `E`, and suppose the output does not
designate `x`. (The safeguard of Step 1a guarantees this: if the output were
a literal of the parameter `x`, the output check of Step 0 would already
succeed and the procedure would have stopped.) Note that the output also
cannot designate the auxiliary variable `Q` used below, for the same reason:
a literal of an auxiliary variable is one of the two forms of the Step 0
check. Since (AR) merely exchanges the names of `x` and `Q`, and the output
designates neither, the deletion below destroys neither the designated value
nor the output. Write the affine system as
`A v = N u + c`. Since `x` occurs only in `E`, the column `N_x` of `N` has a
single nonzero entry, so `A^{-1} N_x != 0`; choose an auxiliary variable `Q`
whose coordinate in `A^{-1} N_x` is `1`, i.e. an auxiliary variable that
depends on `x`. Then `Q = x + phi(u - x)` for an affine map `phi`. Apply (AR)
with this `Q`: the parameter `x` and the auxiliary variable `Q` exchange roles,
no equation changes, and the new auxiliary variable `x` still occurs only in
`E`. Now delete `x` and `E`.

*Invertibility:* the column of `x` in the new `A` has a single nonzero entry,
in the row of `E`, so the determinant expansion of (D3) applies directly and
the remaining auxiliary matrix is invertible. The removed parameter `x` had
weight `6` and is replaced by the new parameter `Q`, whose weight is at most
`6` (it occurs in an equation, because the column of `Q` in the invertible
matrix `A` is nonzero, and it is read by the gate that read `Q` before, so it
is used). No other unpaired parameter increases its weight. Hence
`Delta W <= -1`.

### 3.1 Moving affine gates into the affine system

Suppose a gate `J` of the acyclic part is labelled by an affine operation
`h(t1,t2) = b1 t1 + b2 t2 + c`, where each `t_i` is a constant, a parameter, or
an auxiliary variable. Replace `J` by the equation

    J + b1 t1 + b2 t2 = c                                        (7)

and regard `J` as a new auxiliary variable. The new equation has at most three
nonconstant terms, and the new auxiliary matrix is block triangular with an
extra `1` on the diagonal, hence invertible. The solution of the old system
extends uniquely to the new system with `J` equal to the value the gate
computed, so the computed function is unchanged. One gate became one equation:
`Delta s = 0`. No parameter gains an occurrence, so `Delta W <= 0`.

Applying this to the gates in topological order moves every affine gate whose
two inputs are already in the affine system into that system, while every
affine gate with an AND-type ancestor stays in the acyclic part. This is
finite, since each move removes a gate from the acyclic list.

### 3.2 The implementation of (AR)

Let `Q` be an auxiliary variable depending on an unpaired parameter `x`, so
`Q = x + phi(u - x)`. The parameter `x` must occur in some equation: otherwise
the column `N_x` of `N`, and hence of `A^{-1} N`, would be zero, contradicting
that `Q` depends on `x`.

Substitute `x <- w + phi(u - x)` in the map `Phi`, regard `Q` as the new
parameter `w`, and regard `x` as an auxiliary variable. **No equation or gate
of the circuit changes**; only the roles of `x` and `Q` are exchanged.

*Well-definedness and invertibility.* If `Q` is column `j` of `A` and `N_x`
is the column of `x` in `N`, the new auxiliary matrix `A'` is obtained from
`A` by replacing column `j` by `N_x`; that is, `A' = A B` where `B` is the
identity matrix with its `j`-th column replaced by the vector
`t = A^{-1} N_x`. Since `Q` depends on `x`, the coefficient of `x` in
`Q = x + phi(u-x)` is `t_j = 1`. Hence `det B = t_j = 1` and `A'` is
invertible. Every equation keeps exactly the same terms.

*Effect on the potential.* The new parameter `w` replaces the parameter `x`;
all pairs, equations and gates are unchanged. So `Delta m = Delta q =
Delta a = Delta s = 0`, and `Delta W = wt'(w) - wt(x)`. The parameter `x`
occurs in some equation of the old system, so `wt(x) = 6`. The new parameter
`w` occurs in some equation, since the column of `Q` in the invertible matrix
`A` is nonzero, so `wt(w) = 6`. Therefore `Delta W = 0`.

*Domain.* Since `x` was unpaired, the new domain is the same product of
three-point sets and unrestricted bits, with `w` in place of `x`; no constraint
on a pair is changed.

## 4. The terminal case: a constant output on a large affine subspace

We need one fact about the output of the current circuit.

> **(Lemma 4.1, terminal case.)** Suppose `a > d` and the output of the
> current circuit is either (i) an affine function of the parameters, or
> (ii) a literal of an AND-type gate whose two parents are affine functions of
> the parameters. Then one can construct, in polynomial time, an affine
> subspace `U subset Z(P)` of dimension at least `a - 1 >= d` and a bit `b`
> such that `C(u) = b` for all `u in U`.

*Proof.* Fix one member of each pair to `0` and leave the unpaired parameters
free. This gives an affine subspace `U_0 subset Z(P)` of dimension
`m - q = a`, represented by the point `0` and the `a` standard basis vectors
`e_j` of the free parameters.

Case (i): the output is `l(u) = sum_j lambda_j u_j + kappa`. If `lambda_j = 0`
for every basis vector `e_j` of `U_0`, then `l` is constant on `U_0` with
value `kappa`; take `U = U_0` and `b = kappa`. Otherwise choose a basis vector
`e_{j0}` with `lambda_{j0} = 1`. Then `{u in U_0 : l(u) = 0}` is the affine
subspace `u_{j0} = kappa + sum_{j != j0} lambda_j u_j`, of dimension `a - 1`,
represented by the point `kappa e_{j0}` and the basis `e_j + lambda_j e_{j0}`
for `j != j0`. Take `U` to be this subspace and `b = 0`.

Case (ii): the output equals `(l1 + alpha)(l2 + beta) + gamma`, where `l1,l2`
are affine in the parameters. If both `l1 + alpha` and `l2 + beta` are constant
on `U_0`, with values `c1, c2`, then the output is constant on `U_0` with
value `c1 c2 + gamma`; take `U = U_0` and `b = c1 c2 + gamma`. Otherwise one
of them, say `l1 + alpha`, is nonconstant on `U_0`, and the set
`{u in U_0 : l1(u) + alpha = 0}` is an affine subspace of dimension `a - 1` as
in case (i); on it the output equals `gamma`. Take that subspace for `U` and
`b = gamma`.

In all cases `C(u) = b` for `u in U` and `dim U >= a - 1 >= d`. The
construction compares the coefficients of at most two affine functions against
the basis of `U_0` and writes down a point and at most `a` basis vectors, so it
takes polynomial time. QED

The check for case (i) or case (ii) is syntactic and constant time: inspect the
output and the labels of the one or two values it reads. To decide whether a
value is affine in the parameters, express the auxiliary variables explicitly:
`v = A^{-1} N u + A^{-1} c`, by Gaussian elimination over `F_2`.

## 5. The one-iteration lemma

> **(Lemma 5.1.)** Suppose `a > d`. Then, in deterministic polynomial time,
> either

>  * the algorithm constructs an affine subspace `S subset F_2^n` of dimension
>    at least `d` and a bit `b` with `C_0(x) = b` for all `x in S`; or

>  * the state is replaced by a new state with `m' + s' < m + s`,
>    `Delta a in {0,-1}`, and `Delta W <= 10 Delta a`.

*Proof.* We describe one iteration.

**Step 0 (terminal check).** If the output has one of the two forms of
Lemma 4.1, apply that lemma to obtain `U subset Z(P)` with `dim U >= a - 1 >=
d` and `b` with `C(u) = b` on `U`. Since `Phi` is affine and injective,
`S = Phi(U)` is an affine subspace of `F_2^n` with `dim S >= dim U >= d`, and
by the invariant (3), `C_0(x) = b` for every `x in S`. This is the first
alternative.

**Step 1 (preparation).** Otherwise, repeatedly apply the deletions (D1)-(D4)
and move every eligible affine gate into the affine system (Section 3.1),
until none applies. Each deletion has `Delta a = 0` and `Delta W <= -1`; each
move has `Delta m = Delta q = Delta a = Delta s = 0` and `Delta W <= 0`.
Hence after this step `Delta a = 0` and `Delta W <= 0`.

Throughout Step 1 we have `a > d`: the entry hypothesis of Lemma 5.1 is
`a > d`, and every operation of Step 1 has `Delta a = 0`, so `a` is constant
on Step 1. Lemma 4.1 is therefore applicable at every point of Step 1.

**Step 1a (re-check, and the output safeguard).** The output must remain well
defined at every state, so after *each individual* operation of Step 1 we
re-run the syntactic output check of Step 0: is the output a constant or a
literal of a parameter or of an auxiliary variable, or is it a literal of an
AND-type gate whose two parents are constants, parameters or auxiliary
variables? If it is, we stop the preparation immediately and apply Lemma 4.1,
which returns the first alternative of Lemma 5.1. The check is also run once
more after the last operation of Step 1.

We record why this re-check is what makes the procedure safe. The only
operation of Section 3 that can make the output *undefined* is (D4), which
deletes a parameter `x` together with the equation `E` that is `x`'s unique
child: if the output were a literal of `x`, then `x` could not be deleted
without destroying the output. But if the output is a literal of a parameter
`x`, then the output check of Step 0 *already succeeds*, so the re-check of
Step 1a stops the procedure before (D4) can be applied to `x`. The same
holds for the auxiliary variable `Q` that (D4) uses in its (AR) step, since a
literal of an auxiliary variable is also one of the two forms of the Step 0
check; and since (AR) only exchanges the names of `x` and `Q`, the output,
which designates neither, survives the exchange and the subsequent deletion.
Conversely, if the output check succeeds at a re-check, the output is well
defined and Lemma 4.1 applies, so no further operation is performed.

It is worth being explicit about the order, because it is the order that
makes the argument work: the output check of Step 0 is performed *before* any
operation of Step 1 is attempted, and the re-check of Step 1a is performed
*after* each operation. So at the moment (D4) is applied, the output has
already been tested against both of the forms that would make it designate
`x` or `Q`.

We record, operation by operation, why the output survives Step 1.

* **(D1)** deletes a gate `J` that computes a constant or a literal of one of
  its inputs. The output designation of `J` is a use of `J`, and (D1)
  reroutes *every* use of `J` to that constant or literal. So if the output
  designated `J`, the output now designates the literal, which is still a
  well-defined value. The output therefore survives (D1).
* **(D2)** never deletes the output: its hypothesis is that the gate "is not
  the designated output".
* **(D3)** deletes an equation `E` and the auxiliary variable `v` that `E`
  solves for. If the output designated `v`, the output check of Step 0 would
  already have succeeded, since a literal of an auxiliary variable is one of
  the two forms. Otherwise the output does not designate `v`, and (D3)
  substitutes `v = U + c` everywhere, leaving the output designation
  untouched.
* **(AR)** exchanges the names of a parameter `x` and an auxiliary variable
  `Q`, without changing any equation and without changing the acyclic part.
  If the output designated `x` or `Q`, it now designates the other one,
  which is again a well-defined value, and the output check then succeeds,
  so Step 1a stops the procedure.
* **(D4)** deletes a parameter `x` and the equation `E` that is `x`'s unique
  child. Two separate observations apply. First, by Step 1a the output does
  not designate `x`, nor the auxiliary variable `Q` that (D4) uses, because
  both would have been caught by the output check. Second, even without
  Step 1a, (D4) could not apply to a parameter that the output designates:
  the output designates a *value*, and if it designated `x` then `x`'s
  children would include the value the output designates, so `x`'s unique
  child would not be the equation `E` that (D4) requires. (Recall from
  Section 1.2 that the output is a designation and not a child of the value
  it designates, so this argument uses the fact that the output designates a
  gate when the output check fails.)

Thus the output is well defined at every state of the procedure, and the
re-check adds at most one extra terminal outcome; the potential accounting
is unaffected, because the deletions and moves of Step 1 already have
`Delta a = 0` and `Delta W <= 0`, and stopping early only makes the
potential smaller.

We now check the following three cases in order, and perform the reduction of
the first applicable case; for each case we may therefore assume all earlier
cases are inapplicable.

**Case 1: a paired parameter has degree zero or at least two.** Let `p` be a
paired parameter with partner `r`.

*Subcase `deg(p) = 0`.* Use (PF) to fix `r` to `0`. The pair `{p,r}` is
removed and `p` becomes unpaired. Since `deg(p) = 0`, `p` is unused and
contributes nothing to `W`. Every initially unused parameter stays unused, so
no unpaired parameter increases its weight. We have `Delta m = -1` (`r` is
fixed), `Delta s = 0`, `Delta q = -1`, hence `Delta a = Delta m - Delta q = 0`,
and `Delta W <= -4`.

*Subcase `deg(p) >= 2`.* Use (PF) to fix `p` to `0`, and delete two of its
children using (D1) or (D3). No parameter except the newly unpaired `r` can
increase its weight after the restriction: every rerouting target that is a
parameter was a parent of a deleted object and was already used. The weight of
`r` can increase by at most `6`. Hence `Delta m = -1`, `Delta s = -2`,
`Delta q = -1`, `Delta a = 0`, and

        Delta W <= -4 - 2 + 6 = 0 .

If Case 1 is inapplicable, every paired parameter is isolated, i.e. has
exactly one child.

**Case 2: two partners have the same child.** Let `p` be a paired parameter
with partner `r` and common child `H`. By Case 1, `H` is the unique child of
both `p` and `r`.

If `H` is a gate, its two pins are `p` and `r`. Among the three points of
`Z` with `pr = 0`, namely `(0,0)`, `(0,1)`, `(1,0)`, the gate `H` takes the
same value at two of them, because its two inputs take only two values there.
The affine line through those two points is exactly one of `p = 0`, `r = 0`, or
`p + r = 1`. We restrict to this line using (PF), so that `H` becomes
constant, delete `H` by (D1) and reroute it to its constant value.

If `H` is an equation, then by (D3) it has exactly three nonconstant terms and
contains an auxiliary variable `v` (invertibility of `A`), so it reads
`v + p + r = c`. Imposing `p = r + 1` reduces it to `v = c + 1`, which we
delete by (D3), rerouting `v` to that constant.

In both cases the surviving parameter of the pair becomes unused, since `H` was
its only child; and we reroute only to constants, so no unpaired parameter
increases its weight. Hence `Delta m = -1`, `Delta s = -1`, `Delta q = -1`,
`Delta a = 0`, and

        Delta W <= -4 - 1 = -5 .

If Case 2 is inapplicable, every paired parameter is isolated and has a
different child from its partner.

**Case 3: a paired parameter occurs in the affine system.** Let `p` be a
paired parameter occurring in the affine system, with partner `r`. Since Case 1
is inapplicable, `p` has a unique child `H_p`, and `H_p` must be an equation
(since `p` occurs in the affine system). Similarly `r` has a unique child
`H_r`, and by Case 2, `H_r != H_p`.

Use (PF) to fix `r` to `0`, and delete `H_r`: by (D1) if it is a gate, by (D3)
if it is an equation. Deleting `H_r` creates no new child of `p`, since `p`
does not occur in `H_r`. So `p` is now an unpaired, isolated parameter whose
unique child is the equation `H_p`. If `H_p` has at most two nonconstant
terms, delete it by (D3). If `H_p` has exactly three nonconstant terms, apply
(D4) to the unpaired isolated parameter `p` and the equation `H_p`.

The two objects `H_r` and `H_p` are deleted, and by (D4) the new parameter
replacing `p` has weight at most `6` while `p` had weight `6`; no initially
unpaired parameter increases its weight. Hence `Delta m = -1`, `Delta s = -2`,
`Delta q = -1`, `Delta a = 0`, and

        Delta W <= -4 - 2 + 6 = 0 .

If Cases 1-3 are inapplicable, then every paired parameter is isolated,
occurs only in the acyclic part, and has a different child from its partner.

**Step 2 (restructuring the minimal AND).** Assume `a > d` and that none of
Steps 0-1 applies. The acyclic part is nonempty, since otherwise the output
would be an affine function of the parameters and Step 0 would have fired.
Let `G` be the first gate in the topological order of the acyclic part. By
Section 3.1, `G` is AND-type, since every affine gate whose two inputs lie in
the affine system has been moved into it. By (D1), the two inputs of `G` are
parameters or auxiliary variables with independent linear parts, and by (D2)
`G` has at least one child.

We now make both inputs of `G` parameters, using at most two applications of
(AR).

Suppose one input of `G` is an auxiliary variable `Q`. Since the affine system
has a unique solution, `Q` is an affine function of the parameters, and since
the two inputs of `G` have independent linear parts, `Q` is nonconstant and
depends on some parameter `x`. The parameter `x` occurs in some equation of
the affine system: otherwise the column `N_x` of `N` would be zero, so the
column of `x` in `A^{-1} N` would be zero, contradicting that `Q` depends on
`x`. By Cases 1-3, every parameter occurring in the affine system is unpaired.
So `x` is unpaired, with `deg(x) >= 1`.

Apply (AR) with this `Q` and `x` (Section 3.2): `Delta m = Delta s =
Delta q = Delta a = Delta W = 0`. The new parameter `w` that replaces `Q`
occurs in an equation (the column of `Q` in the invertible matrix `A` is
nonzero) and is a parent of `G`, so it is used; `x` was used as well. If the
other input of `G` is also an auxiliary variable, repeat: the two inputs of `G`
were linearly independent before the first (AR), and an invertible affine
change of parameters preserves that independence, so the second auxiliary
variable depends on a parameter different from `w`. After at most two
applications, both inputs of `G` are parameters; call them `x` and `y`. Write

        G = (x + alpha)(y + beta) + gamma ,   alpha,beta,gamma in F_2 .

If the output is a literal of `G`, Step 0 applies and we stop. Otherwise `G`
is not the output and, by (D2), it has at least one child. Every child of `G`
is a gate in the acyclic part, since no wire goes back from the acyclic part to
the affine system.

If a preliminary reduction became applicable during the restructuring, we
perform it and end the iteration; it is one of Cases 1-3, hence it satisfies
the required bounds.

**Step 3 (eliminating the minimal AND).** We distinguish the two cases on
whether the inputs `x,y` of `G` are paired.

**Case 4: both inputs `x,y` of `G` are unpaired.** We first complement `x` and
`y` if necessary, so that `G = x y + gamma`. Explicitly, if `alpha = 1`,
substitute `x <- w + 1` in `Phi` and replace every occurrence of `x` in the
circuit by `w + 1`; the negation is absorbed into the constant term of each
equation and into the operation of each gate that reads `x`. This changes no
gate count and no degree, so it changes no term of `W`. Do the same for `y` if
`beta = 1`.

Now use (PC) to add the pair `{x,y}` to `P`. Every assignment in the new `Z(P)`
satisfies `x y = 0`, so `G = gamma` on `Z(P)`. Use (D1) to reroute `G` to the
constant `gamma` and delete it. Since `G` has at least one child, and every
child of `G` now reads a constant, we also delete one child of `G` by (D1).

The parameters `x` and `y` were previously inputs of `G`, hence used, and each
had weight `6`. Now both are paired, so they lose these individual weights.
The new pair contributes `+4` to `W`. No other unpaired parameter increases its
weight: `G` is rerouted to a constant, and the deleted child of `G` is
rerouted either to a constant or to a literal of its other input, which was
already used. Altogether

        Delta m = 0 ,  Delta s <= -2 ,  Delta q = +1 ,  Delta a = -1 ,

    and      Delta W <= -2 + 4 - 6 - 6 = -10 .

**Case 5: at least one input of `G` is paired.** By Cases 1-3, each paired
input of `G` has `G` as its unique child, and the partner of each paired input
has a unique child, which must be a gate in the acyclic part different from
`G`. Let `k in {1,2}` be the number of paired inputs of `G`.

*Subcase `k = 1`.* Assume `x` is unpaired and `y` is paired; let `r` be the
partner of `y` and `H_r` the unique child of `r`. Use (PF) to fix `r` to `0`,
which removes the pair `{y,r}`. Then use (UF) to fix `x` to the constant
`alpha`, so that `G` computes the constant `gamma`. Delete `G` by (D1) and
reroute all its occurrences to `gamma`. As `G` was the only child of `y`, the
parameter `y` becomes unused. Delete `H_r` by (D1), since `r` is now the
constant `0`. No initially unpaired parameter increases its weight: `G` is
rerouted to a constant and `H_r` to a constant or to a literal of its other
input, which was already used. The parameter `x` loses its weight `6`; `y` had
no individual weight while paired and has weight `0` afterwards, being unused.
Hence `Delta m = -2`, `Delta s <= -2`, `Delta q = -1`, `Delta a = -1`, and

        Delta W <= -2 - 4 - 6 = -12 .

*Subcase `k = 2`.* Let `r_x`, `r_y` be the partners of `x`, `y`, and let
`H_{r_x}` be the unique child of `r_x`. By Case 2 the pairs `{x,r_x}` and
`{y,r_y}` are distinct. Use (PF) twice to fix `r_x` and `r_y` to `0`, removing
both pairs. Then use (UF) to fix the now unpaired parameter `x` to the constant
`alpha`, so that `G` computes the constant `gamma`. Delete `G` by (D1) and
reroute all its occurrences to `gamma`. As `G` was the only child of `y`, the
parameter `y` becomes unused. Delete `H_{r_x}` by (D1), since `r_x` is now the
constant `0`. No initially unpaired parameter increases its weight. The
parameters `x` and `y` had no individual weight while paired, and contribute
zero afterwards: `x` is fixed to a constant and `y` becomes unused. Hence
`Delta m = -3`, `Delta s <= -2`, `Delta q = -2`, `Delta a = -1`, and

        Delta W <= -2 - 8 = -10 .

Cases 4 and 5 exhaust the possibilities, because every input of `G` is either
unpaired or paired.

**Step 4 (the iteration bound).** In every reducing case we have
`Delta a in {0,-1}` and `Delta W <= 10 Delta a`: in Step 1 the deletions have
`Delta a = 0` and `Delta W <= -1 <= 0`; the moves have `Delta a = 0` and
`Delta W <= 0`; Cases 1, 2, 3 have `Delta a = 0` and `Delta W <= 0`, `-5`, `0`
respectively; Cases 4 and 5 have `Delta a = -1` and `Delta W <= -10`, `-12`,
`-10`. In each case `m + s` strictly decreases: the deletions and Case 4
decrease `s` and leave `m` unchanged, every other reduction decreases `m`, and
no reduction increases `s`. QED

## 6. Termination, the inequality, and the class theorem

Suppose `C_0` is a `B2` circuit computing `f`, and `|C_0| < 4n - 10d`. Run the
algorithm from Section 5.

Initially `m = a = n`, `q = 0`, `s = |C_0|`, and by (6), `W_0 <= |C_0| + 6n`,
so

    W_0 < (4n - 10d) + 6n = 10 (n - d).                            (8)

Suppose the first `t` iterations each returned a new state, with values `a_t`
and `W_t`. Summing `Delta W <= 10 Delta a` over those iterations gives

    W_t <= W_0 + 10 (a_t - n) < 10 (n - d) + 10 (a_t - n) = 10 (a_t - d).

Since `W_t >= 0`, we get `a_t > d`. Hence Step 0 never fires with `a <= d`,
and Lemma 5.1 applies again at the new state.

Each iteration that returns a new state strictly decreases `m + s`, a
nonnegative integer whose initial value is `n + |C_0|`. Therefore after at most
`n + |C_0|` iterations the algorithm must return the first alternative of
Lemma 5.1: an affine subspace `S subset F_2^n` of dimension at least `d` and a
bit `b` with `C_0(x) = b` for all `x in S`.

But `f` is nonconstant on every affine subspace of dimension `d`. Since
`dim S >= d`, `S` contains an affine subspace `T` of dimension exactly `d`;
`f` is nonconstant on `T`, so `f` is nonconstant on `S`. This contradicts
`C_0` being identically `b` on `S`.

Therefore `|C_0| >= 4n - 10d`. This proves the class theorem (T). QED

**Complexity of the algorithm.** `Phi` is stored as `b` and `L`, and `P` as a
list of pairs; all are updated in polynomial time after each operation. Each
iteration performs at most `s` moves of affine gates into the affine system, at
most two applications of (AR), and a constant number of deletions and
reductions, each of which inspects only lists of gates and equations and
Boolean truth tables of constant size. Auxiliary variables are expressed as
`v = A^{-1} N u + A^{-1} c` by Gaussian elimination over `F_2`, which is
polynomial. Since `m <= n` and `s <= |C_0|` throughout, each iteration takes
time polynomial in `n` and `|C_0|`.

**Where the gate count comes from.** Rows are introduced only one-for-one from
a charged gate (Section 3.1), so `s` never exceeds `|C_0|` plus the number of
gates moved in, minus the number deleted; and every row so introduced is
charged to a gate that existed. Affine reparameterisations change roles
without adding computation. Restrictions never add a gate. Deletions only
remove. Hence the inequality bounds the number of gates of the *original*
circuit: no potential value was ever mistaken for a gate count, and the bound
holds for arbitrary sharing, repeated pins, degenerate tables, constants,
unary tables and a single terminal output.

## 7. Summary of the accounting

| Step | `Delta s` | `Delta I` | `Delta q` | `Delta a` | `Delta W` upper bound |
|---|---:|---:|---:|---:|---:|
| (D1)-(D4) | -1 | <= 0 | 0 | 0 | -1 |
| move affine gate | 0 | <= 0 | 0 | 0 | 0 |
| (AR) | 0 | 0 | 0 | 0 | 0 |
| Case 1, `deg(p)=0` | 0 | 0 | -1 | 0 | -4 |
| Case 1, `deg(p)>=2` | -2 | <= +1 | -1 | 0 | 0 |
| Case 2 | -1 | <= 0 | -1 | 0 | -5 |
| Case 3 | -2 | <= +1 | -1 | 0 | 0 |
| Case 4 | <= -2 | -2 | +1 | -1 | -10 |
| Case 5, `k=1` | <= -2 | <= -1 | -1 | -1 | -12 |
| Case 5, `k=2` | <= -2 | <= 0 | -2 | -1 | -10 |

Here `Delta I` is the change in the number of used unpaired parameters. In
Case 4, `x` and `y` each lose a used-unpaired-parameter contribution of `6`
(the `Delta I = -2` row), so that row already includes their loss; in Cases 1
and 3 the newly unpaired partner may gain at most one used parameter
(`+6`).

## 8. Relation to the cited source, and what is reproduced here

The class theorem (T) is Theorem 4.1 of Mishra–Santhanam, ECCC TR26-260
Revision 1 (3 November 2026), *An Explicit 4n - O(log(n)) Size Lower Bound
against General Boolean Circuits*, at
`https://eccc.weizmann.ac.il/report/2026/260/revision/1/download/`, SHA256
`a123943ed05c9569385d3b9d98c695630a110559f33a4c990ec44e2a22b4d921`. That
source's Theorem 4.1 reads: if `f : F_2^n -> F_2` is an affine disperser for
dimension `d`, computable in polynomial time, then any `B2` circuit computing
`f` has size at least `4n - 10d`. Our (T) drops the hypothesis that `f` be
polynomial-time computable, because (T) is a statement about circuits only;
`proof.md` supplies the polynomial-time family that instantiates it.

Everything the source's proof uses is reproduced above rather than cited:
the fair semi-circuit formalism and its two invariants (Section 2), the four
parameter operations with their effect on `Phi` and on `Z(P)` (Section 2.1),
the weight function and the potential (Section 2.2), the four deletions with
their invertibility arguments (Section 3), the implementation of the affine
reparameterisation and its invertibility argument (Sections 3.1-3.2), the
terminal-output lemma with its explicit construction (Section 4), the
one-iteration lemma with all five cases and their inequalities (Section 5),
and the summation argument with termination (Section 6).

**One place where this file is more explicit than the source.** The source
describes its algorithm as: check the output; if the check fails, repeatedly
apply the deletions (D1)-(D4) and move eligible affine gates into the affine
system until none applies; then perform the first applicable case. It does
not spell out a safeguard for the possibility that a deletion in the
preparation phase makes the output *undefined* -- in particular, that (D4)
deletes a parameter together with the equation that is its unique child, which
would destroy the output if the output designated that parameter. The source
states in Section 3.1 that "in our setting, no variable of an equation is
designated as an output", which is the structural reason the output is a
designation rather than a gate, and the source's own Lemma 2.2 is invoked at
the *start* of each step. This file makes the safeguard explicit and
self-contained: Section 2 records that the output is a designation rather than
a child of the value it designates; Section 5 adds **Step 1a**, which re-runs
the syntactic output check after each individual operation of Step 1 and stops
the preparation (returning the terminal alternative of Lemma 5.1) as soon as
it succeeds; and Section 5 gives an operation-by-operation argument that the
output survives every operation, in particular that (D4) cannot apply to a
parameter the output designates. Step 1a is an addition of this file; it is
not taken from the source, and it does not change any bound. It was added in
response to the review of an earlier submission of this package, which
correctly identified the missing safeguard as a reduction defect.

The source's own statement that the output is not designated inside the
affine system is used here only for the fact that the output is a
*designation*; it is not needed to show that the output survives the
operations, which the Step 1a argument shows directly.

The pinned model conventions of `rules.md` come from Li and Yang, ECCC
TR21-023 (front-page date 2021-02-20), Sections 1.1, 2.1, 2.4 (footnotes 3-4)
and 3.3: deterministic acyclic circuits, one output, two input pins per gate,
all sixteen binary truth tables, unrestricted depth and fan-out, exact
computation on every input, nonuniform circuits permitted, inputs `x` alone
rather than NP witnesses, each retained gate counted once including trivial or
degenerate gates, inputs/wires/output selection free, repeated pins and shared
subgraphs allowed, and constant nodes during substitutions handled by
normalization. Our treatment matches these conventions point by point: the
sixteen tables are enumerated in Section 1.1; normalization of constant and
degenerate gates is (D1) and (D3); the terminal-output cases are handled in
Lemma 4.1 rather than by assuming a free extra output gate; and no
uncharged computational basis is introduced anywhere.

## 9. What is not claimed

* This file proves a lower bound of the form `c n - O(d)` for affine
  dispersers. It does not prove a superlinear bound, an endpoint at any
  exponent `p > 1`, or any statement about the class `P` versus `NP`.
* It does not claim a coefficient above `4`. Section 6 of `optimality.md`
  shows, by an explicit exact-arithmetic calculation, that within the family of
  potentials `V = s + beta q + alpha I` and the row estimates of Section 7 no
  choice of weights improves the leading coefficient `4`. That is an
  optimisation bound for these estimates, not a barrier to any other proof.
* It is a written proof, not a formal proof-assistant certificate.
* No checker output, reviewer verdict, or prior acceptance is used as a
  premise of any step.
