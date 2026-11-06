# Gate elimination in the pair domain: C(f) >= 4m - 9(d+1) for affine dispersers

This file proves the elimination theorem used in `proof.md`. The argument is a
gate-elimination in the style of Find–Golovnev–Hirsch–Kulikov and of
Mishra–Santhanam (ECCC TR26-260, Revision 1), re-derived here with the
potential `W = s + 5U + 3q`, for which every step satisfies `Delta W <= 9 *
Delta a`. The factor 9 (in place of 10) is what the bookkeeping below
delivers; the two binding cases are Case 4 and Case 5.2, each attaining -9.

## E.1 States, semi-circuits, invariants

Let `f : F_2^m -> F_2` and let `C_0` be a B_2 circuit computing `f`. We run a
process whose state is a triple `(C, Phi, P)`:

* `Phi : F_2^m' -> F_2^m`, `Phi(u) = b + L u`, an injective affine map; the
  coordinates `u = (u_1, ..., u_m')` are called *parameters*. (We reuse `m`
  for the current number of parameters; the original input length is the
  initial value.)
* `P` is a set of disjoint pairs of parameters. A parameter in a pair is
  *paired*; its mate is its *partner*; other parameters are *unpaired*.
  Write `q = |P|` and `a = m - q` (the *safe dimension*).
* `Z(P) = {u : u_p * u_r = 0 for every {p, r} in P}`: on `Z(P)` at least one
  member of each pair is `0`.
* `C` is a *fair semi-circuit*: an affine system `A v = N u + c` over `F_2`
  in auxiliary variables `v = (v_1, ..., v_r)`, followed by an acyclic list
  of B_2 gates whose inputs are constants, parameters, auxiliaries, or
  earlier gates. Every row of `(A | N)` has at most three nonzero entries
  (*nonconstant terms*), and `A` is invertible (*fairness*), so every
  parameter assignment extends uniquely to all *internal values*
  (parameters, auxiliaries, gates). The *size* is
  `s = (#equations) + (#gates)`. The output is a constant or a literal of an
  internal value.

Invariants, maintained by every operation below:

* (I1) `C(u) = C_0(Phi(u))` for all `u in Z(P)`.
* (I2) `Phi` is injective affine and every operation composes `Phi` with a
  map `tau` sending `Z(P_new)` into `Z(P_old)`.
* (I3) the affine system is fair with at most three nonconstant terms per
  equation.

Initially `Phi` is the identity, `P` is empty, the system is empty and the
acyclic part is `C_0`; thus `m = a =` (input length), `q = 0`, `s = |C_0|`.

The *degree* `deg(t)` of an internal value is its number of distinct
children (equations containing it, plus gates reading it; a gate reading `t`
on both pins is one child). A parameter is *unused* if `deg = 0` and
*isolated* if `deg = 1`. Let `U` = number of unpaired parameters with
`deg >= 1`. The potential is

```
W = s + 5U + 3q  >= 0,        W_initial <= |C_0| + 5m.
```

*Rerouting*: when an internal value becomes a constant or a literal of
another internal value, replace all its occurrences by that constant or
literal. In every deletion below, rerouting targets a value that already
occurs in the deleted gate or equation; hence a parameter unused before a
deletion stays unused after it, and no unused parameter acquires weight.

## E.2 Operations on parameters

Each operation replaces `Phi` by `Phi o tau`; injectivity and (I2) are
immediate in each case.

* **(AR) Affine reparameterisation.** For an unpaired `x` and affine `phi`,
  substitute `x <- w + phi(u_{-x})` with `w` fresh: `Delta m = Delta q =
  Delta a = 0`. Implemented by exchanging the roles of `x` and an auxiliary
  `Q` whose solved form is `Q = x + phi(u_{-x})`: the equations are
  unchanged as relations, the column of `x` in `N` and the column of `Q` in
  `A` are swapped, and the new auxiliary matrix is `A' = A B` where `B` is
  the identity with column `j(Q)` replaced by `A^{-1} N_x`; since the
  coefficient of `x` in `Q`'s solved form is `(A^{-1} N_x)_{j(Q)} = 1`,
  `det B = 1`, so `A'` is invertible and (I3) is preserved. Both `x` (it
  occurs in some equation, else no auxiliary could depend on it) and `w`
  (the column of `Q` in the invertible `A` is nonzero, so `w` occurs in
  some equation) are used before/after, hence `Delta W = 0`.
* **(UF) Fixing an unpaired parameter.** Impose `x = c` for a constant
  `c`, or `x = y + c` for another parameter `y` and constant `c` (the only
  forms ever used below; substituting either form into an equation
  replaces the term `x` by a constant or by the single term `y`, so the
  bound of three nonconstant terms per equation is preserved):
  `Delta m = -1`, `Delta q = 0`, `Delta a = -1`.
* **(PF) Fixing a paired parameter.** For `{p, r} in P`, substitute one of
  the three lines `p = 0`, `r = 0`, `p + r = 1` contained in `pr = 0`, and
  remove the substituted parameter and the pair:
  `Delta m = -1`, `Delta q = -1`, `Delta a = 0`.
* **(PC) Pair creation.** For distinct unpaired `x, y`, add `{x, y}` to
  `P`: `Delta m = 0`, `Delta q = +1`, `Delta a = -1`; `Phi` unchanged and
  `Z(P)` shrinks, preserving (I1).
* **Complementation.** Substitute `x <- w + 1` in `Phi` and absorb the `+1`
  into the constant terms of equations and the operations of gates reading
  `x`: degrees and weights are unchanged, `Delta W = 0`.
* **Moving an affine gate into the system.** If a gate `J` has an affine
  operation `h(t_1, t_2) = b_1 t_1 + b_2 t_2 + c` and both inputs are
  constants, parameters, or auxiliaries, replace `J` by the equation
  `J + b_1 t_1 + b_2 t_2 = c` with `J` as a new auxiliary. The new row has
  at most three nonconstant terms; the new `A` is block triangular with a
  `1` on the diagonal, hence invertible; `s` is unchanged (one gate becomes
  one equation) and degrees of parameters are unchanged, so `Delta W = 0`
  and `Delta a = 0`.

Fixing a parameter to a constant or to another parameter plus a constant
(the uses of (UF)/(PF) below) only deletes a column of `N` or substitutes
inside equations; it never increases the number of nonconstant terms and
never changes `A`'s invertibility.

## E.3 Deletions

Each deletion has `Delta m = Delta q = Delta a = 0` and `Delta s = -1`,
hence `Delta W <= -1` (weights cannot increase, by the rerouting remark).

* **(D1)** A gate computing a constant or a literal of one of its inputs:
  reroute it and delete it. The instances used are recognisable in
  polynomial time: every constant or unary gate, every gate with a
  constant input, every gate with two identical inputs, every AND-type
  gate whose two input *values* have dependent affine forms in the
  parameters (then one value is the other plus a constant, and the gate is
  a literal of one of its inputs — of the second if the first form is
  constant; the forms are computed from
  `v = A^{-1} N u + A^{-1} c`), and every gate that the steps below make
  constant on all of `Z(P)` by fixing an input or adding a pair (reroute
  to that constant; sound by (I1)). The accounting `Delta W <= -1` holds
  for any such deletion, recognisable or not.
* **(D2)** A gate with no children that is not the output: delete it.
* **(D3)** An equation with at most two nonconstant terms. It contains an
  auxiliary variable `v` (otherwise its row of `A` would be zero,
  contradicting fairness), so it reads `v = U + c` with `U` a parameter, an
  auxiliary, or `0`. Substitute `v = U + c` everywhere and delete `v` and
  the equation. No equation gains nonconstant terms (a term is replaced by
  a constant or a single term), and fairness is preserved: substituting in
  other equations is a sequence of elementary row operations on
  `(A | N | c)`, after which the column of `v` has a single nonzero entry,
  in the deleted row; expanding `det A` along that column shows the
  remaining auxiliary matrix is invertible.
* **(D4)** An unpaired, isolated parameter `x` whose unique child is an
  equation `E`, provided the output is not a literal of `x` (if it is, the
  output is an affine function of the parameters and Lemma E1 applies
  instead; in the Case 3 use below this proviso is automatic, since the
  terminal check has excluded an output that is a literal of the paired
  parameter `p`, and no step of Case 3 reroutes anything to `p`, as `p`
  does not occur in `H_r`). The column `N_x` has a single nonzero entry (row `E`), so
  `A^{-1} N_x != 0`; pick an auxiliary `Q` with `(A^{-1} N_x)_Q = 1`, i.e.
  `Q = x + phi(u_{-x})` in solved form. Apply (AR) to exchange `x` and
  `Q`. Now `x` is an auxiliary occurring only in `E` (its column in the new
  `A` is the old `N_x`) and read by no gate (it was isolated with child
  `E`). Delete `x` together with `E`: fairness of the remainder follows by
  expanding the determinant along the column of `x`, which has its single
  nonzero entry in row `E`. The removed `x` had weight `5`; the new
  parameter `Q` has weight at most `5`; `Delta W <= -1`.

## E.4 Terminal lemma

**Lemma E1.** Suppose `a > d` and the output of `C` is (i) an affine
function of the parameters, or (ii) a literal of an AND-type gate whose two
inputs are affine functions of the parameters. Then in polynomial time one
can find an affine subspace `U ⊆ Z(P)` of dimension at least `a - 1 >= d`
and `b in F_2` with `C(u) = b` on `U`; consequently `S = Phi(U)` is an
affine subspace of dimension at least `d` on which `C_0` is constantly `b`.

*Proof.* Fix one member of each pair to `0`; this gives `U_0 ⊆ Z(P)` of
dimension `m - q = a` with point `0` and the standard basis in the unfixed
coordinates. An affine function `l(u) = sum lambda_j u_j + kappa` is
constant on `U_0` iff `lambda_j = 0` on all basis coordinates; if not, its
zero set inside `U_0` is an affine subspace of dimension `a - 1` (solve for
one basis coordinate with `lambda_j = 1`). In case (i) apply this to the
output itself. In case (ii) the output is `(l_1 + alpha)(l_2 + beta) +
gamma`; if both `l_1 + alpha, l_2 + beta` are constant on `U_0` so is the
output; else the output equals `gamma` on the zero set of a nonconstant
one of them, an affine subspace of dimension `a - 1`. Finally `Phi` is
injective affine, so `dim Phi(U) = dim U >= a - 1 >= d`, and (I1) gives
`C_0 = b` on `S`. ∎

## E.5 The step lemma

**Lemma E2.** Suppose `a > d`. Then in polynomial time the process either
produces `(S, b)` as in Lemma E1, or reaches a new state with `m + s`
strictly smaller, `Delta a in {0, -1}` and `Delta W <= 9 * Delta a`.

*Proof.* One iteration begins with a preparation loop. Run the terminal
check: if the output is a constant, a literal of a parameter or
auxiliary, or a literal of an AND-type gate whose two parents are
constants, parameters, or auxiliaries (all affine functions of the
parameters, auxiliaries via `v = A^{-1} N u + A^{-1} c`), Lemma E1
applies and we stop with `(S, b)`. Otherwise exhaust (D1)–(D4) and move
every affine gate with no AND-type ancestor into the system, and run the
terminal check again: the moves can create the configurations of
Lemma E1 that were absent before them (an output that was an affine gate
becomes a literal of an auxiliary; an AND gate whose parents were affine
gates acquires auxiliary parents), so the check is repeated after the
moves, not only before them. Every pass that does not stop strictly
decreases `s` or the number of acyclic gates, so the loop ends with no
deletion or move applicable and the terminal check failing. Then perform
the first applicable of Cases 1–3; if none applies, restructure and
eliminate the minimal AND (Cases 4–5), re-running the preparation loop
and the check of Cases 1–3 after each exchange of the restructuring, as
detailed there. Preparatory operations (moves, (AR), deletions) have
`Delta a = 0` and `Delta W <= 0`; each case below satisfies
`Delta W <= 9 Delta a` and decreases `m + s`; hence so does the whole
iteration.

**Case 1: a paired parameter has degree `0` or `>= 2`.** Let `{p, r} in P`.
If `deg(p) = 0`, apply (PF) fixing `r = 0`. Then `Delta m = Delta q = -1`,
`Delta s = 0`, `p` becomes an unpaired *unused* parameter, no weight
increases: `Delta W <= 3 * (-1) = -3`, `Delta a = 0`.
If `deg(p) >= 2`, apply (PF) fixing `p = 0` and delete two distinct
children of `p`: each child is a gate with a constant input, deleted by
(D1), or an equation that drops to at most two nonconstant terms, deleted
by (D3). Only `r` can gain weight (it becomes unpaired; at most `+5`).
`Delta W <= -2 (deletions) - 3 (pair removed) + 5 = 0`, `Delta a = 0`.

**Case 2: two partners share a child.** Let `p, r` be partners with a
common child `H`; by Case 1 both are isolated, so `H` is the unique child
of each. If `H` is a gate, its inputs are `p, r`; `H` takes the same value
at two of the three points `(0,0), (1,0), (0,1)` of `pr = 0`; the line
through those two points is one of `p = 0`, `r = 0`, `p + r = 1`. Restrict
by (PF) to that line; `H` becomes constant and is deleted by (D1). If `H`
is an equation, (D3) being inapplicable it has exactly three nonconstant
terms, one of them an auxiliary (fairness), so `H` reads `v + p + r = c`;
impose `p = r + 1` by (PF), reducing `H` to `v = c + 1`, deleted by (D3).
In both subcases the surviving partner's only child is gone, so it is
unused; rerouting is to constants only. `Delta m = Delta q = -1`,
`Delta s = -1`: `Delta W <= -1 - 3 = -4`, `Delta a = 0`.

If Cases 1–2 do not apply, every paired parameter is isolated and has a
different child from its partner.

**Case 3: a paired parameter occurs in the affine system.** Let `p` be
paired, occurring in the system, with partner `r` and unique children
`H_p` (an equation) and `H_r != H_p`. Apply (PF) fixing `r = 0` and delete
`H_r` by (D1) (gate with constant input) or (D3) (equation dropping to two
terms); `H_r` does not contain `p`, so this creates no new child of `p`.
Now `p` is unpaired and isolated with unique child the equation `H_p`. If
`H_p` has at most two nonconstant terms, delete it by (D3); else it has
exactly three, and (D4) applies to `p` and `H_p`. Weight changes over
the whole case: fixing `r` makes `p` an unpaired used parameter (`+5`);
in the (D3) subcase `p` becomes unused again when `H_p` is deleted
(`-5`), and in the (D4) subcase `p` is replaced by a parameter of weight
at most `5`. In both subcases the net weight change is at most `+5`.
`Delta m = Delta q = -1`, `Delta s = -2`:
`Delta W <= -2 - 3 + 5 = 0`, `Delta a = 0`.

If Cases 1–3 do not apply, every paired parameter is isolated, occurs only
in the acyclic part, and has a child different from its partner's.

**Restructuring the minimal AND.** After the preparation loop, the
acyclic part is nonempty: if it were empty, the output would be a
constant or a literal of a parameter or auxiliary, and the terminal check
would have stopped us. Let `G` be its first gate in topological order.
Every affine gate with no AND-type ancestor was moved into the system,
and the first gate has no gate ancestors, so `G` is AND-type; by
(D1)-exhaustion its two inputs are parameters or auxiliaries whose
affine forms in the parameters have independent linear parts (in
particular both are nonconstant). Moreover `G` is not the output: the
output being a literal of `G` is a literal of an AND gate whose parents
are parameters or auxiliaries, which the terminal check (run after the
moves) would have caught; and `G` has a child, by (D2)-exhaustion.
If an input of `G` is an auxiliary `Q`, then `Q`'s form depends on some
parameter `x`; `x` occurs in some equation (else the column of `x` in
`A^{-1} N` would vanish, contradicting the dependence), and `x` is
unpaired (a paired parameter occurring in the system would trigger
Case 3). Apply (AR) to exchange `Q` and `x` (`Delta W = 0`). For the
second input, choose the parameter from the support of its form
*different from the parameter introduced by the first exchange*: such a
parameter exists, for if the second form depended only on that one
parameter, its linear part would be `0` or that parameter's linear part,
contradicting the independence of the two forms (independence is
preserved by the first exchange, an invertible reparameterisation). So
after at most two exchanges both inputs of `G` are parameters, say
`x, y`, and `G = (x + alpha)(y + beta) + gamma`. After each exchange we
re-run the preparation loop and re-check Cases 1–3; if the terminal
check or a case fires, the iteration ends with that accounting (the
exchanges contributing `Delta W = 0`, `Delta a = 0`).

**Case 4: both inputs of `G` are unpaired.** Complement `x` and `y` if
needed (`Delta W = 0`) so that `G = xy + gamma`. Add the pair `{x, y}` by
(PC): now `xy = 0` on `Z(P)`, so `G = gamma` there; delete `G` by (D1),
and delete one child of `G` by (D1) (it now has the constant input
`gamma`). Both `x, y` were used and unpaired (weight `5` each) and become
paired; the new pair contributes `+3`; no other weight increases (reroutes
are to a constant or to a literal of an input already present).
`Delta m = 0`, `Delta s <= -2`, `Delta q = +1`, `Delta a = -1`:
`Delta W <= -2 - 10 + 3 = -9 = 9 * Delta a`.

**Case 5: at least one input of `G` is paired.** Each paired input has
`G` as its unique child, and its partner has a unique child, a gate in the
acyclic part different from `G` (Cases 1–3 inapplicable).

*Subcase 5.1: exactly one input paired.* Say `x` unpaired, `y` paired with
partner `r` and `r`'s unique child `H_r`. Apply (PF) fixing `r = 0`
(pair removed), then (UF) fixing `x = alpha`, so `G` computes the constant
`gamma`. Delete `G` by (D1) and delete `H_r` by (D1) (constant input
`r`); if `G` has a child other than `H_r`, one such child may be deleted
by (D1) as well (constant input `gamma`), but the bound does not need it.
After `G` is deleted, `y` is unused. Weight changes: `x` (used, unpaired)
is removed, `-5`; the removed pair gives `-3`; no other weight increases.
`Delta m = -2`, `Delta s <= -2`, `Delta q = -1`, `Delta a = -1`:
`Delta W <= -2 - 5 - 3 = -10 <= 9 * Delta a`.

*Subcase 5.2: both inputs paired.* Let the partners be `r_x != r_y` (the
pairs are distinct, by Case 2) with unique children `H_{r_x}, H_{r_y}`.
Apply (PF) twice, fixing `r_x = 0` and `r_y = 0` (both pairs removed),
then (UF) fixing `x = alpha`, so `G` computes `gamma`. Delete `G` by
(D1), and delete the children `H_{r_x}` and `H_{r_y}` by (D1) (each has a
constant input). These are three distinct deletions: if
`H_{r_x} != H_{r_y}` there is nothing to check; if
`H_{r_x} = H_{r_y} = H`, then `H` is a child of both `r_x` and `r_y`, so
its two pins are exactly `r_x` and `r_y` (a gate has only two pins), in
particular `H` does not read `G`, and one child of `G` — which exists,
since `G` is not the output and (D2) is exhausted — is a third distinct
gate with the constant input `gamma`, deleted by (D1) in its place.
After `G` is deleted, `y` is unused; `x, y` carried no individual weight
(paired) and contribute none after. `Delta m = -3`, `Delta s <= -3`,
`Delta q = -2`, `Delta a = -1`: `Delta W <= -3 - 6 = -9 = 9 * Delta a`.

Every step is a list operation or a Gaussian elimination over `F_2`
(expressing auxiliaries as affine functions of the parameters), so each
iteration runs in polynomial time; deletions and Cases 1–3, 5 decrease
`s` or `m`, Case 4 decreases `s`, and no step increases `m + s`. ∎

## E.6 The elimination theorem

**Theorem E3.** Let `f : F_2^m -> F_2` be an affine disperser for
dimension `d` (nonconstant on every affine subspace of dimension at least
`d`; equivalently of dimension exactly `d`, since a larger subspace
contains one of dimension `d`). Then every B_2 circuit computing `f` has
size at least `4m - 9(d + 1)`.

*Proof.* Let `C_0` compute `f`. Run the process of Lemma E2. Suppose for
contradiction that the initial potential satisfies `W_0 <= 9(m - d - 1)`.
Summing `Delta W <= 9 Delta a` over any number of state-returning
iterations gives `0 <= W_t <= W_0 + 9(a_t - m) <= 9(a_t - d - 1)`, so
`a_t >= d + 1 > d` at every iteration: Lemma E2 stays applicable, but
`m + s` strictly decreases, so the process cannot return states forever;
it must exit through the `(S, b)` outcome of Lemma E1. Then `C_0` is
constantly `b` on the affine subspace `S` of dimension at least `d`,
while `f`, a disperser for dimension `d`, is nonconstant on `S` —
contradiction with `C_0` computing `f`. Hence `W_0 >= 9(m - d - 1) + 1`,
and since `W_0 <= |C_0| + 5m`,

```
|C_0| >= 9(m - d - 1) + 1 - 5m = 4m - 9(d + 1) + 1 >= 4m - 9(d + 1).  ∎
```

Provenance: the semi-circuit framework and the parameter-domain viewpoint
follow Find–Golovnev–Hirsch–Kulikov (FOCS 2016; Comput. Complexity 2023)
and Mishra–Santhanam (ECCC TR26-260, Revision 1, Sections 2–4), whose
published analysis uses the weights `(6, 4)` and charge `10`. The present
accounting with weights `(5, 3)` and charge `9`, including the deletion of
both partners' children in Subcase 5.2 and of one child of `G` in
Subcase 5.1, is the derivation given here; every case is verified above
directly from the definitions.
