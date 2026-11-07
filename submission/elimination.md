# Full B2 gate lower bound for affine dispersers

## 1. Affine maps with promised pairs: the 4n-10d theorem

This self-contained argument follows Mishra--Santhanam, ECCC TR26-260 revision 1, November 3, 2026,
Sections 2--4, specialized to the pinned full B2 convention. 

**Theorem.** For every total Boolean f on n bits nonconstant on every
dimension-d affine subspace, d>=1, any full-B2 circuit computing f has
at least 4n-10d gates. No uniformity hypothesis is needed for this finite
gate statement. If d>=n its right side is nonpositive, so assume n>d.

The 16 binary Boolean tables consist of two constants, four literals of one input, two XOR-type tables u+v+c, and eight AND-type tables (u+b)(v+c)+e over F_2, with b,c,e in F_2. The latter eight have exactly one exceptional table entry (weight one or three); the remaining two nondegenerate balanced tables are XOR and its complement. Thus this classification covers every nondegenerate gate. Internal complements are absorbed into descendant truth tables. Allowing a final output literal for free strengthens the lower-bound model; no extra computational gate is inserted without charge.

### 1.1 New state and restricted correctness domain

Keep m independent binary coordinates u and an injective AFFINE map
Phi(u)=b+L*u into the original n input coordinates. Keep disjoint pairs P,
q=|P|, a=m-q, and the domain

 Z(P)={u: p*h=0 for every pair {p,h} in P}.

The represented computation is a fair affine system A*v=N*u+c (invertible
square A, at most three nonconstant terms per row), followed by acyclic
full-B2 gates with no feedback to the system. Its output is a constant
or literal of a coordinate, auxiliary value, or acyclic gate. Count one
object for each equation and acyclic gate, total s. The correctness
invariant is C(u)=C_original(Phi(u)) for ALL u IN Z(P). No agreement is
asserted outside Z(P). Initially Phi is the identity, P is empty, and
s is the original circuit's gate count.

Four coordinate operations preserve this invariant when accompanied by
the specified computation changes:

* Exchange an unpaired x with an affine value Q=x+ell(u_without_x),
  using x=w+ell in Phi. Paired coordinates do not change, so the
  bijection takes the new Z(P) onto the old Z(P); m,q,a stay fixed.
* Fix an unpaired x to a constant, or another coordinate's literal when
  needed. This is an injective affine restriction preserving all pair
  constraints. It lowers m and a by one.
* For {p,h} in P, impose p=0, h=0, or p=h+1, and remove the eliminated
  coordinate and the pair. These are exactly the three lines inside
  p*h=0. Each preserves all other pairs, lowers m and q by one and
  preserves a. Only these paired restrictions are permitted here.
* Add a pair of distinct unpaired coordinates. Phi and m stay fixed;
  Z(P) shrinks, q rises by one and a falls by one. This operation may
  justify making a gate constant ON THE NEW DOMAIN.

All coordinate maps are injective and affine; compositions therefore keep
Phi injective and affine. The paired restriction p=h+1 is legal even
though neither coordinate is fixed individually: h(h+1)=0 over F2.
Complementing an unpaired coordinate before pairing is also an affine
bijection preserving the old Z(P).

Let I count used unpaired coordinates, where usage means occurrence in a
row or as a gate input; output references are not counted. Give every
used unpaired coordinate weight 6, every unused one weight 0, and paired
coordinates no individual weight. Use

 W=s+4q+6I.

W is nonnegative. Initially W<=|C_original|+6n. Since the two possible
unpaired weights depend only on whether a coordinate is used, arbitrary
rerouting fanout cannot increase a coordinate's weight if it was already
used. Rerouting to a parameter literal only targets a parameter occurring
in the removed row or gate, so an originally unused coordinate remains
unused. This is the reason the old degree-1 versus degree-2 penalty is
absent in this new invariant.

### 1.2 Terminal outputs and computation-preserving simplifications

Fix one member of each pair to zero. The resulting affine subspace U0 is
inside Z(P), of dimension m-q=a. If the output is affine in u, it is
constant on an affine subspace of U0 of dimension at least a-1: either
it is already constant, or impose one affine equation on U0. If the
output is an AND-type expression in two affine values, write it
(l1+alpha)(l2+beta)+gamma. If a factor is nonconstant on U0, impose
that factor equal to zero, leaving dimension a-1 and constant output.
If both factors are constant, the output is already constant. These
arguments cover final literals as well. Thus while a>d, neither kind
of output can occur for a circuit computing f: its constant restriction
would give a dimension-d constant subspace after injective affine Phi.

The following deletions preserve m,q,a and lower W by at least one.
The elementary row and determinant arguments are supplied here.

1. Reroute and delete a constant or single-input-literal gate. This
   includes constant/unary tables, repeated pins, fixed inputs, and ANDs
   on affine values with dependent linear parts. A gate known constant
   only on Z(P) can likewise be replaced by that constant, because
   correctness is required only there. No previously unused unpaired
   coordinate becomes used; deleting the gate lowers s by one.
2. Delete a non-output acyclic gate with no children. Its removal cannot
   create a use of any coordinate.
3. For a row with at most two terms, invertibility ensures an auxiliary
   v occurs. Write v=U+c, with U another single value or zero. Reroute
   v to U+c and delete v and this row. Adding the row to rows containing
   v leaves one nonzero in its auxiliary column; determinant expansion
   shows the remaining auxiliary matrix is invertible. Replacing a term
   by one term cannot increase row sparsity or introduce a previously
   unused coordinate. This removes one charged object.
4. If unpaired x has exactly one computational occurrence, in a row E, Nx is nonzero. Choose an
   auxiliary Q whose coefficient of x in A^(-1)N is 1 and exchange x,Q.
   The new auxiliary matrix replaces Q's column by Nx. Multiplication
   by A^(-1) gives the identity with one column replaced by a vector
   with a 1 at its diagonal position, so it remains invertible. No row
   term changes. Old x becomes an auxiliary occurring only in E; remove
   x,E by determinant expansion on its single-occurrence column. Its
   output reference cannot be the designated output while a>d, by the
   terminal argument. The old parameter had weight 6, its replacement
   has weight at most 6, and a row is removed, so Delta W<=-1.

Affine gates whose inputs are constants, coordinates or affine-system
values can be moved into the system as one new defining row. The new
auxiliary matrix is block triangular with diagonal entry 1, and the row
has at most three terms. This exchanges one gate for one equation, hence
preserves s. It does not change whether any parameter is used.

For a general affine exchange Q=x+ell with unpaired x, both x before
and Q after occur in the system: Nx!=0, and the old Q column of A is
nonzero. Both parameter weights are therefore exactly 6. Swapping roles
changes no row or gate terms, so Delta W=0 and s,m,q,a stay fixed.
This exchange is valid regardless of old x's degree. Only Phi records
the affine substitution; no dense expression or free gate is inserted
into the computational representation.

### 1.3 Paired cleanup, before choosing a minimal AND

Apply deletions and move eligible affine gates, then inspect these cases
in order. Each paired cleanup has Delta a=0 and Delta W<=0 and strictly
lowers m+s. Rerouting may remove additional objects; counts below are
lower bounds for distinct original objects, not pin counts.

(A) Paired p has partner h. If p is unused, impose h=0 and erase h and
the pair. Surviving p becomes unpaired but unused, so Delta W<=-4.
If p has at least two distinct children, impose p=0 and erase p and
its pair. Delete two affected children. Acyclic children have a fixed
pin; row children become short and use deletion 3. If propagation removes
one early it is already a counted deletion; remaining affected rows
stay short. Only newly unpaired h can acquire a new individual weight,
at most 6. Thus Delta W<=-4-2+6=0. After (A) fails, every paired
coordinate has exactly one child.

(B) Partners p,h share their one child H. If H is a gate, evaluate its
full two-input Boolean table on (0,0),(0,1),(1,0). Two outputs agree.
The line through the two corresponding points is one of p=0, h=0,
p=h+1. Restrict to that line, making H constant, then delete H. If H
is a row it has the form v+p+h=c, since it has three terms and some
auxiliary. Impose p=h+1 and delete v=c+1. The surviving coordinate
is now unused. All rerouting is to constants. Therefore Delta W<=-4-1.
After (B) fails, partners have different children.

(C) Paired p occurs in a row Hp, and partner h has unique child Hh!=Hp.
Impose h=0, remove h and the pair, and delete Hh. Parameter p now is
unpaired and still has at most the one row occurrence Hp: Hh did not
contain p. If Hp already disappeared by propagation, it is the second
counted deletion; at most one newly unpaired coordinate contributes weight.
Otherwise delete a short Hp,
or use deletion 4 when it has three terms. This removes a second object
and leaves at most one new unpaired coordinate of weight 6. Other old
unpaired coordinates cannot become newly used. Thus Delta W<=-4-2+6=0.
After (C) fails, all paired coordinates are isolated, occur only in the
acyclic part, and have a different child from their partner.

### 1.4 Exhaustive minimal-AND step

If a>d, some acyclic gate remains: otherwise the output is affine.
Choose the first one in topological order, G. It is AND-type because
eligible affine gates were moved into the system. Its inputs are affine
values with independent linear parts, or deletion 1 would apply.
If a parent Q is auxiliary, its solved affine form depends on some
parameter x. Such x occurs in a row and is unpaired by absence of (C).
Exchange Q,x at zero W cost. Check deletions and paired cases again;
if any applies, perform it and end this iteration. Otherwise a second
exchange turns the other parent into a distinct parameter: independence
supplies a coordinate other than the first new input. At most two
exchanges suffice. We now have G=(x+alpha)(y+beta)+gamma on distinct
parameters. It is not the output by the terminal argument and has at
least one acyclic child. There are three possibilities.

(D) Both x,y are unpaired. Complement either as necessary, in Phi and
all its computational occurrences, to make G=x*y+gamma. This preserves
usage and does not change any already paired coordinate. Add pair
{x,y}; on the new Z(P), G equals gamma. Delete G and one of its
children, which has a constant input. Shared descendants are rerouted
once; G and its child are two distinct original objects. Both inputs
lose their weight 6 and the new pair contributes 4. No old unused
coordinate becomes used. Thus

 Delta m=0, Delta q=1, Delta a=-1, Delta s<=-2,
 Delta W<=-2+4-6-6=-10.

(E1) Exactly one input is paired; write x unpaired, y paired with h.
The unique child Hh of h differs from G by (B). Impose h=0 (erase
its pair), then fix unpaired x=alpha to force G constant. Delete G
and Hh. These are distinct, even if Hh is a child of G or they have
shared descendants. Coordinate y is unused afterwards, because G was
its only child. The old unpaired x loses weight 6; newly unpaired y
ends with weight zero. No old unused parameter becomes used. Hence

 Delta m=-2, Delta q=-1, Delta a=-1, Delta s<=-2,
 Delta W<=-2-4-6=-12<=-10.

(E2) Both inputs are paired. Their pairs are distinct; if they were
partners, they would share G and (B) would apply. Let hx,hy be their
respective partners. Both have unique children different from G.
Impose hx=hy=0 and erase both pairs. Fix now unpaired x=alpha,
forcing G constant. Delete G and the unique child Hx of hx. They
are distinct; hy's child may equal Hx, and is not counted again.
Initially x,y had no individual weights. Finally x is removed and y
is unused, so neither contributes an individual weight. Thus

 Delta m=-3, Delta q=-2, Delta a=-1, Delta s<=-2,
 Delta W<=-2-8=-10.

All other initially unpaired coordinates retain weight at most their old
weight in (E1)--(E2). This is because each nonconstant rerouting targets
an already occurring parent. Intermediate release of paired coordinates
need not be monotone; the displayed Delta W compares whole steps.

### 1.5 Termination and conclusion

These cases cover every state while a>d: either a terminal form gives a
constant dimension-d affine subspace, a deletion/paired cleanup applies,
or the first acyclic AND is handled by (D), (E1) or (E2). Preparations
are finite: topological affine moves, at most two exchanges, and a
restart when a strict deletion/cleanup is found. Each successful step
has Delta a in {0,-1}, Delta W<=10 Delta a, and strictly lowers m+s.
No step increases m or s. Thus the process cannot run indefinitely.
For a circuit computing the disperser, a terminal contradiction is
impossible, so it reaches a_T<=d. Summing the potential inequalities,

 |C_original|+6n >= W_0 >= W_T+10(n-a_T) >=10(n-d).

This proves |C_original|>=4n-10d. All computation counts refer to original
charged gates or their one-for-one row replacements. Affine exchanges
add no gates, and the promised-domain pair operation changes only the
validity domain for the representation. The argument permits arbitrary
sharing, every binary Boolean truth table, final literals and constant
normalization. In particular it is a full-B2 bound, not a basis-restricted
or formula bound, and it assumes no direct-sum additivity.

## 2. Counting the constant closure: 4n-8d

The following finite refinement retains the affine promised-domain state of
Section 1. It uses one additional deletion count in the two-paired-input
case. This refinement is supplied here; coefficient 4 itself is credited
to the recent Mishra--Santhanam argument in Section 1.

**Strengthened theorem.** Under the same hypotheses, every circuit
computing f has at least 4n-8d gates.

### 2.1 Four objects in case (E2)

Use the state at the beginning of Section 1.4(E2), after all preliminary
deletions and paired cleanups are unavailable. Let G read paired x,y,
whose distinct partners hx,hy have unique children Hx,Hy. Both children
are acyclic gates and differ from G. Gate G has a child, and every
non-output gate has at least one child, since dead-gate deletion was
exhausted. Impose hx=hy=0 and then x=alpha, as in (E2). The new a is
old a-1>=d. Any constant output after these restrictions would therefore
make the original f constant on a dimension-d affine subspace: the new
Z(P) contains its usual affine space of dimension a_new=a_old-1. Such a constant
output is impossible for a circuit computing the disperser.

Under these restrictions, G is constant. Gates Hx,Hy have fixed pins,
so they can be removed by literal/constant rerouting. In fact at least
FOUR DISTINCT ORIGINAL gates can be removed. The proof counts sharing
explicitly.

First suppose Hx!=Hy. We already have three distinct affected objects
G,Hx,Hy. Choose a child D of G. If D is different from Hx and Hy,
these three together with D are four affected gates. Otherwise rename
so D=Hx. Since Hx reads hx and G, and both are fixed constants, Hx
is constant. If it is the output we have the constant-output contradiction
above. Otherwise it has a child E. Acyclicity makes E different from
G,Hx. If E!=Hy, it is the fourth gate. If E=Hy, then Hy's two inputs
are Hx and hy, both constants, so Hy is constant. It cannot be the
output; hence it has a child F. The path G,Hx,Hy,F makes F distinct
from all three, and supplies the fourth gate. Each counted child has
a constant pin when removed.

Now suppose Hx=Hy=B. Since B reads both distinct coordinates hx,hy,
its two inputs are exactly those values; in particular B is not a child
of G. Both G and B become constant. Choose a child D of G, giving
three distinct gates G,B,D. Gate B cannot be a constant output, and
has some child E. If E!=D, then G,B,D,E are four distinct affected
gates: E cannot equal G because G reads x,y directly, and cannot equal
B by acyclicity. Otherwise B's child E=D reads both G and B and is
therefore constant. It cannot be the output, so a child F exists.
Acyclicity makes F distinct from G,B,D, supplying the fourth gate.
This also covers the situation when B has further children: any such
child different from D already gives the fourth object.

All selected objects were present before the step. If propagation removes
one earlier, it is already one of the counted deletions. Constant
propagation along the displayed paths cannot add computation or create
a use of an originally unused parameter. Nothing assumes disjoint
subcircuits, separate partner children, or separate descendants.

### 2.2 Reweight all cases

Keep a=m-q, but use potential

 V=s+2q+4I,

where I still counts used unpaired parameters, regardless of their degree.
V is nonnegative and initially V<=|C|+4n. The deletions, affine moves
and exchanges of Section 1 remain legal. Every deletion lowers V by
at least one: old used unpaired weight 4 is replaced by at most 4, and
no old unused coordinate becomes used. A general affine exchange swaps
two used weights of 4 and costs zero.

For completeness, the complete step estimates with these weights are:

* Paired unused cleanup (A): Delta a=0, Delta V<=-2.
* Paired degree-at-least-two cleanup (A): erase a pair, remove two
  objects, release at most one unpaired weight 4. Delta a=0 and
  Delta V<=-2-2+4=0.
* Partners with the same child (B): Delta a=0, Delta V<=-2-1=-3.
* Paired row cleanup (C): Delta a=0, Delta V<=-2-2+4=0.
* Both minimal-AND inputs unpaired (D): add one pair, remove G and
  a child, lose two unpaired weights 4. Delta a=-1 and
  Delta V<=-2+2-4-4=-8.
* Exactly one input paired (E1): remove a pair and two objects, and
  lose one unpaired weight 4; the released surviving input is unused.
  Delta a=-1 and Delta V<=-2-2-4=-8.
* Both inputs paired (E2): erase two pairs, remove at least four
  distinct objects by 2.1, and end with no new used unpaired weight
  from x,y. Delta a=-1 and Delta V<=-4-4=-8.

The whole-step weight accounting is used, including temporary release of
paired coordinates. Terminal arguments and strictly decreasing m+s are
unchanged. While a>d each successful step therefore satisfies
Delta V<=8 Delta a. A terminal constant affine subspace is impossible
for a circuit computing f, so summing until a_T<=d gives

 |C|+4n >= V_0 >= V_T+8(n-a_T) >=8(n-d),

and hence |C|>=4n-8d. The sole extra hypothesis used by the strengthened
count was precisely the already available a>d and disperser property;
no new circuit model, source hardness or uniformity assumption is added.

## 3. Count permanently inert coordinates: 4n-8d+4

**Finite refinement.** For d>=1, every full-B2 circuit computing a
function nonconstant on every dimension-d affine subspace has at least
4n-8d+4 gates. This improves the preceding finite bound by four gates;
it does not change its leading coefficient.

Use exactly the terminating reduction already proved. If n<=d, the new
right side is at most 4-4d<=0 and the claim is immediate. Otherwise each
step changes a by 0 or -1, so the terminal value is exactly a_T=d.

Call an unpaired coordinate inert if it appears in no equation, no gate
input, and no output reference. Once inert, it remains inert throughout
this reduction. Literal reroutings cannot introduce it because they use
an already occurring parent. An affine exchange uses a solved auxiliary
value independent of every inert coordinate: its column in N is zero.
Hence the inverse substitution also does not introduce an inert coordinate
into the computational representation. Pair creation only uses gate parents,
and restrictions only fix used or paired coordinates. No later operation
uses an inert coordinate as a new parameter or output value.

Count the following events by N_0,N_B,N_1,N_2 respectively:

* paired-unused cleanup (A), which leaves the unused member unpaired;
* common-child cleanup (B), which leaves its surviving coordinate unused;
* case (E1), which leaves its former paired gate input y unpaired and unused;
* case (E2), which leaves its surviving gate input y unpaired and unused.

Each event creates a distinct permanently inert coordinate. None is an
output coordinate. Before the event, a>d rules out an output literal of
a parameter by the terminal argument. During the event, a rerouting cannot
introduce that survivor as the output: its only old child is either made
constant (B,E1,E2), or it has no old child at all (A-unused). Other removed
objects did not read it. Thus these coordinates satisfy the full inertness
definition, not merely the degree-zero convention that ignores output.
Let E=N_0+N_B+N_1+N_2.

There are E distinct inert unpaired coordinates at termination, so E+q_T<=a_T=d.
Also E<=d-1: fix one coordinate of every remaining pair to zero and fix all
other parameters except the E inert ones. This is an affine E-dimensional
subspace inside Z(P), on which the represented output is constant. Its image
under injective affine Phi has the same dimension, so E>=d would contradict
the disperser property.

For a separate counting ledger, put J=s-2q. It need not be nonnegative.
The already proved deletions and preparations give these complete changes:

* ordinary deletion: Delta J<=-1, Delta a=0;
* paired-unused cleanup: Delta J<=2, Delta a=0;
* paired degree-at-least-two cleanup: Delta J<=0, Delta a=0;
* common-child cleanup: Delta J<=1, Delta a=0;
* paired-row cleanup: Delta J<=0, Delta a=0;
* case (D): Delta J<=-4, Delta a=-1;
* case (E1): Delta J<=0, Delta a=-1;
* case (E2): Delta J<=0, Delta a=-1.

Affine moves and exchanges preserve J. For (E2) the four-object count of
Section 2 is essential; for (E1) its original two-object count suffices.
Summing gives, with J_0=|C|,

    |C| >= 4(n-d)+s_T-2q_T-(2N_0+N_B+4N_1+4N_2)
         >= 4n-4d+s_T-2q_T-4E.                         (7)

If q_T=0, E<=d-1 makes (7) at least 4n-8d+4.
If q_T>=2, E<=d-q_T makes it at least
4n-8d+2q_T+s_T>=4n-8d+4.

For q_T=1, the remaining possibility E=d-1 is impossible. There are
exactly d-1 unpaired coordinates (since a_T equals the number of unpaired
coordinates plus q_T); they would all be inert. Hence the output depends
only on the two coordinates p,h in the sole pair, regardless of how many
gates or rows remain. On the three legal points (0,0),(0,1),(1,0), two
output values agree. Their affine line is one of p=0,h=0,p=h+1, and lies
inside the promised domain. Varying that line and all d-1 inert coordinates
gives a constant-output affine space of dimension d, a contradiction.
Thus for q_T=1 we have E<=d-2, and (7) is at least
4n-8d+6+s_T, which is stronger than required.

This proves 4n-8d+4 for all n and d>=1. The ledger uses the same legal
full-B2 reductions, with their sharing and terminal arguments, and tracks
inertness including output references. It neither assumes nonnegativity
of J nor identifies it with gate count. The new refinement is independent
of any source's prior review outcome.

