# A coefficient-four gate bound from three-point pair domains

## Statement and credit

**Theorem.** Let `d` be a nonnegative integer, and let
`f : F_2^n -> F_2` be nonconstant on every affine subspace of dimension at
least `d`. Every deterministic acyclic, one-output circuit for `f` over
all sixteen binary Boolean truth tables, with arbitrary sharing, has at
least

    4n - 9(d+1)

gates. Input variables, wires, and output selection have zero cost; each
retained gate, including a trivial gate, has cost one. This theorem is a
class theorem; eligibility and an all-length family instantiating it must
be supplied separately.

The three-point-pair-domain idea and the potential `s+5I+3q` are taken
from the public, unjudged note titled "Coefficient 4: three-point pair
domains give 4n - 9(d+1) for every affine disperser", associated with
submission `5e9fcf7f-8338-4729-9faa-a18994661af6`. That note credits GPT 6
Astra, GPT 6.1 Sol, and Claude Opus 5.5 for the idea, accounting, and
writing. The proof below reconstructs and proves the claim independently;
neither its submission status nor its reported audits are premises.
It uses the fair affine-system/parameter framework of Mishra--Santhanam,
ECCC TR26-260 (19 October 2026), Sections 2--4. Everything needed from
that framework is proved below. Three-point fibers were also used in
weighted gate elimination by Golovnev--Kulikov, ECCC TR15-170 (ITCS 2016);
no theorem from that work is a premise here.

All arithmetic and affine maps below are over `F_2`.

## 1. Representation and the domain invariant

We work in a relaxation of the circuit model, starting from the given
circuit without increasing its size. There are `M` independent Boolean
parameters, a collection `P` of `q` disjoint pairs of these parameters,
and a domain `D` in parameter space. For every pair `{p,h}` exactly one
corner `(alpha,beta)` is excluded; the allowed values of that pair are
the other three points of `F_2^2`. Unpaired parameters range freely.
The domain is the product of these pair domains and the free bits.

There is an affine map

    Phi : F_2^M -> F_2^n

which is injective on the **entire** parameter cube, not just on `D`.
There are auxiliary variables `v` specified by an affine system

    A v = N u + c,

where `A` is square and invertible. Each equation has at most three
distinct nonconstant terms after cancellation. An acyclic collection of
binary Boolean gates follows the system. Gate inputs can be parameters,
auxiliaries, or earlier gates. No gate is an input to an equation.
The output is a constant or a literal of a parameter, auxiliary, or
gate. Permitting the final output literal for free is a relaxation of
the stated model and therefore strengthens the lower-bound theorem.

The representation is required to compute `f(Phi(u))` for `u in D`.
Its behavior outside `D` is unrestricted. Auxiliaries nevertheless have
unique affine values on the entire cube because `A` is invertible.
The representation's size `s` is its number of equations plus retained
gates. Initially `M=n`, `Phi` is the identity, `q=0`, there are no
equations, and `s` is the original circuit size.

A value's degree is the number of distinct equations in which it occurs
as a nonzero term, plus the number of distinct gates having it as a
parent; repeated pins count once. Equations have no distinguished
output term for this definition. The output reference is not counted.
Let `I` be the number of **used unpaired
parameters**, namely unpaired parameters of positive degree. Put

    a = M-q,                 W = s + 5I + 3q.

Always `W>=0`, and initially `W<=s_initial+5n`. We retain unused
parameters in `M` unless a specified restriction removes them; syntactic
nonuse does not by itself justify a dimension decrease or increase.

For a pair excluding `(alpha,beta)`, fixing `p=1+alpha` leaves `h`
free, and fixing `h=1+beta` leaves `p` free. We call these the safe
values. Fixing one member of every pair safely gives an affine subspace
of the parameter cube, contained in `D`, of dimension `M-q=a`.
Its image under `Phi` is an affine subspace of exactly that dimension.
In particular, a constant output is impossible whenever `a>=d`.

We perform reductions only while

    a >= d+2.                                      (guard)

Two useful terminal exclusions follow directly, without requiring that
the following affine slices themselves be implemented as reductions.
On a safe affine subspace of dimension `a`, any affine output is
constant on a subspace of dimension at least `a-1`. Thus a constant or
affine output, including a literal of a parameter or auxiliary, is
impossible under the guard. An AND-type operation on two affine values
also has a constant affine slice of codimension at most one: write it
`(L+b)(R+c)+e`. If `L` is nonconstant, set `L=b`; if `L` is constant,
the operation is either constant or affine in `R`, to which the previous
argument applies. Consequently a gate of this form with affine parents
cannot be the output under the guard. Each resulting slice has dimension
at least `a-1>=d`, and `Phi` preserves its dimension.

## 2. Boolean tables and legal simplifications

The sixteen two-input tables consist of two constants, four literals
of one input, two tables `U+V+e`, and eight tables

    (U+b)(V+c)+e.

These last eight are called AND-type. Their majority value is `e`,
and their sole exceptional input is `(1+b,1+c)`. Constant nodes introduced
by substitution are free. Negations are absorbed into descendant truth
tables or into affine equation constants; a terminal negation is covered
by the relaxed output-literal convention. No retained negation or
projection gate is silently omitted from the size count.

Here are four simplifications. Each deletes at least one gate or
equation, keeps `M,q,a` unchanged, and does not increase `I`.

**D1: a degenerate gate.** A gate which is a constant or a literal of
one actual parent is rerouted to that value and deleted. This includes
constant or unary tables, fixed pins, and repeated pins. It also applies
when both parents are affine values with dependent linear parts: they
are constant or equal up to a constant, so composing the two-input table
gives a constant or a literal of a nonconstant actual parent. A rerouted
parameter was already used in the deleted gate, so it cannot newly
enter `I`.

**D2: a dead gate.** Delete a gate of degree zero unless it supplies the
output. This cannot increase parameter usage.

**D3: a short equation.** An equation with at most two distinct terms
contains an auxiliary, since otherwise its row in `A` would be zero.
Solve it as `v=U+c` with `U` either one other value or zero. Substitute
in all other equations, gates and the output, and delete this equation
and `v`. Substitution replaces a single term by at most one term, so
equation lengths do not increase. Ordinary row elimination followed by
expansion along the `v` column proves that the remaining auxiliary
matrix is invertible. If `U` is a parameter, it was already used in the
deleted equation, so `I` does not increase.

**D4: an isolated unpaired parameter in an equation.** Suppose an
unpaired parameter `x` has just one child, an equation `E`. Some
auxiliary `Q` has coefficient one on `x` in its affine solution, because
the column `N_x` is nonzero and `A^{-1}N_x` is nonzero. Exchange the
roles of `x` and `Q`, leaving the equations' terms unchanged. Taking
`w=Q` as the new unpaired parameter, write

    Q=x+ell(other parameters),   x=w+ell(other parameters),

and make this affine substitution in `Phi`. The coordinate change is
a bijection on the full cube. It preserves the product domain: all
paired coordinates stay fixed by the change, and for every choice of
them the old free coordinate `x` is exchanged bijectively with the new
free coordinate `w`. Hence `Phi` stays affine and injective.

The new auxiliary matrix replaces the column `A_Q` by `N_x`. Its
determinant is nonzero because the `Q` coordinate of `A^{-1}N_x` is
one. The old parameter `x`, now auxiliary, occurs only in `E`; it is
not the output, since that would be affine under the guard. Delete it
and `E`. Expanding the determinant on its single-occurrence column
proves invertibility of the remaining matrix. No gate or equation has
been added. Before deletion a used unpaired parameter was exchanged
for another used unpaired parameter: `Q` occurs in some equation
because the original column `A_Q` is nonzero. Thus `I` cannot increase.

All four simplifications give `Delta W<=-1` and `Delta a=0`.
The deletion arguments count original distinct objects. In particular,
if propagation already deletes a second designated object, it counts
as deleted; if it survives, a short equation stays short under D3 and
a fixed-pin gate stays degenerate until it is deleted. These observations
will justify simultaneous deletion counts below.

An affine gate whose parents are parameters or auxiliaries can instead
be moved into the system: add its defining equation and make its output
a new auxiliary. The new matrix is block triangular with a final
diagonal entry one; the equation has at most three terms. This removes
one gate and adds one equation, so `s` does not change, and `I` cannot
increase. Repeat in topological order whenever possible. This is a
transfer of one charged object to another, not an additional free
computation.

The role exchange used in D4 is also legal without the final deletion,
whenever an unpaired `x` has coefficient one in the solution for an
auxiliary `Q`. It preserves `s,M,q,a,I`, for the same reasons. We call
this an affine role exchange. The possibly dense affine expression is
substituted only in `Phi`; the computational equations retain their
original terms. No dense affine gate is inserted for free.

## 3. Clean up the paired parameters

At the start of each iteration, move eligible affine gates and perform
any available D1--D4 simplification. A simplification is a successful
iteration, after which we restart. Otherwise try the following pair
cases in order, restarting after any successful case.

**P0: an unused member.** If paired `p` is unused, fix its partner `h`
safely, and remove `h` and this pair. The surviving free `p` is unused.
Thus `Delta M=Delta q=-1`, `Delta a=0`, `Delta I<=0`, and
`Delta W<=-3`.

**P2: a member of degree at least two.** Fix such a member `p` safely,
remove it and the pair, and delete at least two distinct children of
`p`. A gate child has a fixed pin and is deleted by D1; an equation
child becomes short and is deleted by D3. The other member becomes an
unpaired parameter, so it is the only possible new member of `I`.
Consequently

    Delta s<=-2,  Delta I<=1,  Delta q=-1,  Delta a=0,
    Delta W<=-2+5-3=0.

When neither P0 nor P2 applies, every paired member has exactly one
child.

**PC: partners have the same child.** Suppose partners `p,h` have the
same unique child `B`. If `B` is a gate, its Boolean table has the same
value on two of the three allowed points, by the pigeonhole principle.
Those two points form an affine line contained in the allowed domain.
Restrict the pair to that line. The line is expressible by fixing one
coordinate or identifying the coordinates up to a constant; remove one
parameter and the pair. Now `B` is constant and can be deleted. If `B`
is an equation, its other term is an auxiliary and it has form
`v+p+h=c`. Among the three allowed points two have the same `p+h`;
on their line `v` is constant, and D3 deletes the equation. The
surviving parameter is unused. Therefore

    Delta s<=-1,  Delta I<=0,  Delta q=-1,  Delta a=0,
    Delta W<=-4.

The restricted `Phi` is still affine and injective. This includes the
case that the deleted gate supplies the output: if the resulting
constant output conflicts with the domain invariant, the supposed
configuration simply cannot occur under the guard.

**PE: a paired member occurs in an equation.** Let `p` have sole child
`E`, an equation, and let its partner `h` have sole child `B`. By PC,
`E` and `B` are distinct. Fix `h` safely and remove it and the pair.
Delete `B` using D1 or D3. It contained no `p`, since `p` had only
child `E`; none of these substitutions creates new occurrences of `p`.
If `E` survives, `p` is an unpaired parameter isolated in `E`, so D4
deletes `E`. If `E` has already disappeared, `p` is unused and it is
already the second distinct deletion. The only possible increase in
`I` is one, for `p` or its replacement in D4. Thus

    Delta s<=-2,  Delta I<=1,  Delta q=-1,  Delta a=0,
    Delta W<=0.

If none of these pair cases applies, every paired parameter has exactly
one child, that child is a gate, and partners have distinct children.
In particular no paired parameter occurs in any equation. All auxiliary
values therefore depend only on unpaired parameters.

## 4. Expose a first AND-type gate

Assume none of the preceding reductions applies under the guard. The
acyclic gate part is nonempty, since otherwise the output would be
affine. Its first gate `G` has only parameter or auxiliary parents.
It is neither degenerate nor affine, by D1 and the topological moves,
and hence is AND-type. Its parents' affine linear parts are independent,
or D1 would apply.

If one parent `Q` is auxiliary, its affine form depends on an unpaired
parameter `x`: it is nonconstant, and paired parameters occur in no
equation. Use an affine role exchange to make `Q` an unpaired parameter.
If the other parent is still auxiliary, linear independence guarantees
that its form involves some unpaired parameter other than the first
new parent; exchange using that parameter. Thus at most two exchanges
make both parents of `G` parameters. Paired coordinates and all their
occurrences are unchanged. Newly exposed degeneracies or eligible affine
gates can be simplified before continuing; otherwise these exchanges
have zero cost and the same `G` has two distinct parameter parents.

The terminal exclusion in Section 1 shows that `G` is not the output.
D2 consequently gives at least one gate child of `G`. Gates cannot
have equation children. No preparation which reroutes children of `G`
away from its input parameters is needed.

Exactly three cases remain, according to the pair status of its parents.

**GP: both parents are unpaired.** Write

    G=(x+b)(y+c)+e.

Restrict `(x,y)` to the three-point majority fiber, excluding
`(1+b,1+c)`, and record this as a new pair. No parameter or coordinate
map is changed. On the new domain `G=e`, so delete `G` and each of
its original children by constant substitution and D1. There is at
least one child, and it is distinct from `G`. The previously used
unpaired parameters `x,y` cease to contribute to `I`. Ordinary
rerouting cannot create another used parameter, since any parameter
rerouted from a deleted gate was already used there. Hence

    Delta s<=-2,  Delta I<=-2,  Delta q=1,  Delta a=-1,
    Delta W<=-2-10+3=-9.

The new computation need only be correct on the new domain. In
particular `G` is not asserted constant on the full parameter cube.
The coordinate map remains affine and injective on that full cube.

**AP1: exactly one parent is paired.** Let the parents be paired `p`
and unpaired `w`, and let `h` be the partner of `p`. Let `H` be the
unique child of `h`; by the pair cleanup `H` is a gate distinct from
`G`. Fix `w` to the value forcing the AND-type `G` to its majority
constant. Also fix `h` safely in its pair. Remove `w,h` and that pair.
The surviving `p` is free and becomes unused after `G` is removed.
Delete `G` and `H`: the former is constant, and the latter has its
`h` pin fixed. These are two distinct deletions, even if `H` is a
child of `G` or they share descendants. No third deletion is needed.
Removing the used unpaired `w` loses one member of `I`; the only
newly unpaired coordinate is unused `p`. Thus

    Delta s<=-2,  Delta I<=-1,  Delta q=-1,  Delta a=-1,
    Delta W<=-2-5-3=-10<=-9.

This slightly weaker estimate than the public note's `-11` avoids any
need to count a successor of `H` and suffices for the theorem.

**AP2: both parents are paired.** They are in different pairs, since
PC excludes partners sharing `G`. Write the pairs as `{p,h}` and
`{r,k}`, with `G` reading `p,r`. Let `H,K` be the unique children
of `h,k`, respectively. They are gates, each distinct from `G`.
Choose a valid point of the first pair with `p` equal to a forcing
value for `G`, and fix that whole pair to this point. Such a point
always exists: either value of `p` occurs in the three-point domain.
The second member `r` becomes unused when the now constant `G` is
deleted; fix its partner `k` safely, so `r` remains a free unused
parameter. These restrictions remove three parameters and two pairs.

Delete `G,H,K`. If `H!=K`, these are three distinct gates. If `H=K`,
their common gate reads `h` and `k`: it has two distinct pins and both
parameters have it as their unique child. Therefore this common gate
cannot also read `G`. Since `G` has a child, choose one such child `J`;
then `G,H,J` are three distinct gates. The pins `h,k` are fixed, so
`H` is deletable, and the `G` pin of `J` is constant, so `J` is
deletable. This proves at least three distinct deletions in every
overlap pattern. The only surviving newly unpaired parameter is `r`,
which is unused; no old unpaired parameter becomes used. Hence

    Delta s<=-3,  Delta I<=0,  Delta q=-2,  Delta a=-1,
    Delta W<=-3-6=-9.

All these restrictions preserve the product form of the domain and
affine injectivity of `Phi`: a safe single-coordinate fixation removes
one parameter and one pair, a whole-pair fixation to a valid point
removes two parameters and one pair, and an unpaired fixation removes
one parameter. In particular their displayed changes in `a` are exact.

## 5. Coverage, termination, and the inequality

Every successful reduction satisfies

    Delta a in {0,-1},        Delta W <= 9 Delta a.

The zero-dimensional-loss reductions are D1--D4 and P0, P2, PC, PE;
the unit-loss reductions are GP, AP1, AP2. Each preserves all
representation invariants, never increases `M` or `s`, and strictly
decreases the nonnegative integer `M+s`. Each preparation is finite:
affine moves transfer existing gates into equations without introducing
new gates, and at most two role exchanges expose the selected first
AND-type gate. If a simplification becomes applicable during a
preparation, perform it as the successful reduction and restart.

Coverage is exhaustive. After simplifications and pair cleanup, a
nonempty acyclic part has a first nondegenerate nonaffine gate, whose
parents can be exposed as parameters. It cannot be the output, and
the status of its two parents is exactly one of GP, AP1, AP2.
If the acyclic part were empty, the output would be affine, which is
excluded by the guard. Output constants or literals produced by a
reduction are handled by the same domain invariant; if such an output
would conflict with that invariant, the purported pre-reduction state
was impossible. No estimate presumes an extra successor of an output.

There are at most `n+s_initial` successful iterations. Consequently the
process must leave the guard, with a final `a_T<=d+1`. If `n<=d+1`
initially, the asserted lower bound is already nonpositive and is
immediate. Otherwise sum the displayed inequalities over the process:

    s_initial + 5n >= W_0
                      >= W_T + 9(n-a_T)
                      >= 9(n-d-1),

since `W_T>=0`. Therefore

    s_initial >= 4n - 9(d+1),

as claimed.

The size in this inequality is the original number of charged gates.
Every equation was introduced one-for-one from an existing charged
gate. Role exchanges change the representation and coordinate map
without adding computational objects. Every other step is a restriction
or a deletion. The argument thus covers all sixteen tables, constants,
unary and repeated-pin gates, terminal outputs, arbitrary fan-out,
shared descendants, and unrestricted nonuniform circuits in the stated
full-B2 model.
