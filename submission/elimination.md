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


## 4. One additional gate for d>=2: 4n-8d+5

**Theorem.** If d>=2 and f is nonconstant on every dimension-d affine
subspace, every full-B2 circuit for f has at least 4n-8d+5 gates.
This is only a one-gate finite refinement; its leading coefficient is four.
For n<=d the right side is at most 5-4d<0, so assume n>d.
Run the preceding reduction, stopping immediately when a first reaches d.

### 4.1 Three deletions in (E1)

In (E1), fix h=0 and unpaired x=alpha as before. Gate G becomes constant,
and the distinct gate Hh has its h pin fixed. Choose a child D of G.
If D!=Hh, then G,Hh,D are three distinct deleted original gates.
If D=Hh, both pins of Hh are constants (h and G), so Hh is constant.
The new a is old a-1>=d. Thus Hh cannot be the output: a constant output
would make f constant on a dimension-d affine space in the new promise.
Otherwise Hh has a child E, distinct from G,Hh by acyclicity; E has a
constant pin and is deleted. Again three distinct objects are removed.
Earlier propagation counts a removed original object once. The existing
whole-step inertness and weight accounting is unchanged. In particular
case (E1) has Delta J=Delta(s-2q)<=-3+2=-1, rather than merely <=0.

With K the number of ordinary deletion steps, summing the ledger gives

 |C| >= 4(n-d)+s_T-2q_T
          -(2N_0+N_B+3N_1+4N_2)+K
      = 4n-4d+s_T-2q_T-4E+2N_0+3N_B+N_1+K.       (8)

Extra deletions in any whole step can only strengthen (8). Counts are
integers. Affine moves and exchanges do not change J. The event counts,
permanent inertness, E<=d-1, and E+q_T<=d are as proved in Section 3.

### 4.2 An elementary affine-plane pigeonhole fact

Any coloring of 18 distinct points of F_2^5 with two colors contains a
monochromatic affine plane. One color has at least nine points. Nine
points determine 36 unordered distinct-point pairs, but only 31 nonzero
pairwise sums are available. Two pairs have the same sum. They cannot
share one endpoint without being identical, so their four endpoints are
distinct and sum to zero. Those four points form an affine plane, and
all have the same color. In particular this applies to

 {x*hx=0, y*hy=0} x F_2,

which has 3*3*2=18 points on five independent coordinates. No retrieval,
probabilistic assumption or finite enumeration is needed for this fact.

### 4.3 The tight q_T=0 case must include (E1)

Suppose q_T=0 and E=d-1. The last increase of E must occur at the last
paid step. Indeed, at a state with E=d-1 and a>d, delete the inert
coordinates from consideration. If U is the number of unpaired
coordinates, there are 2^(U-E)*3^q legal active assignments and
(U-E)+q=a-E>=2, so there are at least three legal active assignments.
Two have the same Boolean output. Their two-point affine line lies in
the promised domain; varying the E inert coordinates yields a constant
affine space of dimension E+1=d, a contradiction. Thus E cannot reach
d-1 while a>d. The process stops immediately at a=d.

Ordinary deletions, affine moves, exchanges, (A) with used coordinate,
(C), and (D) create none of the four counted inert events. (A-unused)
and (B) do not lower a, and cannot be the event first reaching E=d-1
while a>d. Consequently the last paid step is (E1) or (E2).

It cannot be (E2). Immediately before such a last step, q=2,
a=d+1 and E=d-2. The number of unpaired coordinates is U=a-q=d-1,
so exactly one unpaired coordinate is outside the E tracked inert ones.
After ignoring the inert coordinates, the legal active domain is the
18-point domain on five coordinates in 4.2. Correctness and inertness
make the Boolean output a function of those active coordinates alone.
By 4.2 it is constant on an affine plane inside that domain. Varying
all d-2 inert coordinates gives a constant affine space of dimension d.
The injective affine map Phi preserves its dimension, contradicting
f's disperser property. Therefore the last step is (E1), and N_1>=1.

Equation (8), with q_T=0, E=d-1 and N_1>=1, now gives

 |C| >= 4n-8d+4+s_T+N_1 >=4n-8d+5.

If instead q_T=0 and E<=d-2, (8) gives the stronger 4n-8d+8.

### 4.4 Other terminal pair counts

For q_T=1, Section 3 proved E<=d-2. Equation (8) gives
|C|>=4n-8d+6+s_T, which suffices.

For q_T=2, if E<=d-3 the same equation gives at least 4n-8d+8.
If E=d-2, we must have s_T>=1. Otherwise the output is a constant or
literal of a parameter, with no equations or gates. All unpaired
coordinates are exactly the E inert coordinates, so an output literal
cannot refer to one of them. For a literal of a paired coordinate, set
that coordinate to zero and let its partner vary; in the other pair,
fix either member to zero and let its partner vary. This is a legal
constant-output affine plane in the active paired coordinates. A
constant output likewise has such a plane. Varying the E inert
coordinates then gives a forbidden constant dimension-d affine space.
Thus s_T>=1, and (8) gives 4n-8d+4+s_T>=4n-8d+5.

Finally q_T>=3 and E<=d-q_T give, directly from (8),
|C|>=4n-8d+2q_T+s_T>=4n-8d+6.
These cases exhaust the terminal states and prove the theorem.

Only the elementary plane fact and the refined (E1) count are new to
this section. The unrestricted full-B2 model, permitted affine maps,
all-length family, and coefficient-four asymptotic guarantee are unchanged.


## 5. Terminal geometry and an unrestricted Ramsey lemma

The following geometric facts and Lemma R are used in the additive-seven
refinement proved in Section 6 below. Its finite refutation is regenerated
and checked by the declared audit; no SAT verdict or previous review is
a premise.

### 5.1 Two elementary terminal facts

An affine function on q pairs is constant on a product of q legal lines: choose two equal-contribution points on each pair. Here is another useful pigeonhole fact. If E=d-2 and a>=d+2 at any
state, a constant dimension-d affine subspace is unavoidable. Ignore
inert coordinates, and select four independent units, where an unpaired
coordinate is one unit and a promised pair is one unit. Fix all other
active coordinates to zero. If r of the four units are pairs, the
remaining domain has respectively

 r:                    0    1    2    3    4
 points:              16   24   36   54   81
 ambient dimension:    4    5    6    7    8
 same-color points:    8   12   18   27   41
 unordered pairs:     28   66  153  351  820
 nonzero sums:        15   31   63  127  255.

In every column, equal pair sums among same-colored points give a
monochromatic plane, exactly as in 4.2. Adjoin E inert directions and
apply affine injective Phi. This proves the fact. The same argument for
three units with at least two pairs uses 18 points in dimension 5 or
27 points in dimension 6: respectively 36>31 or 91>63 pair sums. Thus
when E=d-2,a=d+1, one must have q<=1.

### 5.2 An unrestricted terminal Ramsey lemma

**Lemma R.** On independent coordinates x,y,hx,hy,p,hp,u,z let
D={x*hx=y*hy=p*hp=0}. For any AND-type G=(x+alpha)(y+beta)+gamma
and ANY Boolean H(G,hx,hy,p,hp,u,z), with no other dependence on x,y,
there is a constant affine three-flat inside D. There is no circuit-size,
formula or complexity restriction on H. The domain has 108 points.

If alpha=0, fix x=0 and hy=0, legal line choices in distinct pairs.
Then G=gamma and y is unused. The other five coordinates hx,p,hp,u,z
range over {p*hp=0}, a 24-point subset of F_2^5. At least twelve points
have one output, and their 66 unordered pairs exceed the 31 nonzero
differences. Two distinct pairs of the same difference have four distinct
endpoints and form an affine plane. Adjoin the unused y direction to
obtain the three-flat. The beta=0 case is symmetric. This proof handles
three of the four forcing orientations analytically for arbitrary H.

For alpha=beta=1 absorb gamma into H's first argument and set
G=(x+1)(y+1), with no computational gate change asserted. We supply a
small exact refutation of every putative coloring with no constant
three-flat. The following reductions explain its complete scope.

**Quotient and clauses.** Index original points by
x+2y+4hx+8hy+16p+32hp+64u+128z. The 108 legal points project to 84
distinct states (G,hx,hy,p,hp,u,z). Every H induces an arbitrary binary
color on these states. The audit assigns state IDs in first-encounter
order while scanning legal integer points increasingly, using key
(G,a>>2). It enumerates all 14391 legal three-flats. For each flat,
project its eight points and remove repeated state IDs. A coloring
without a monochromatic flat must satisfy both clauses

    OR_(v in projected flat) color_v,
    OR_(v in projected flat) NOT color_v.

There are 11349 distinct projected sets. This is an exact encoding:
no monochromatic projected set means that very original flat is not
monochromatic, even when different original points share a state.
No value of H outside the 84-state image is needed.

**Normalize a twelve-point core.** Write t=p+2hp+4u+8z and order its
legal values increasingly, omitting p=hp=1. A layer is a twelve-bit
coloring indexed by these t values. Let C0 be the layer G=0,hx=hy=0.
Fix x=1,hx=hy=0 and vary y. G stays zero and y is unused. Thus any
constant affine plane in C0 lifts to a forbidden three-flat. C0 must
avoid monochromatic planes in this twelve-point domain.

The audit enumerates its 39 affine planes by all four distinct points
whose XOR is zero, and all 4096 binary colorings. Exactly 48 avoid
monochromatic planes. This is also transparent from their fibers:
each pair of the three (p,hp) states, together with the two free bits,
is a three-cube. A coloring with no constant plane in a three-cube
must be balanced, because five same-colored points give 10 pairs but
only seven nonzero differences. Hence each free four-point fiber is
balanced. Such a fiber is a nonconstant affine function of u,z. Two
fibers cannot have the same linear part, since then their union is an
affine function on a three-cube and has a constant plane. Conversely,
two different nonzero linear parts yield a balanced quadratic on the
three-cube with no constant plane: its polar rank is two and its linear
part is nonzero on the one-dimensional radical. Alternatively the
stated finite enumeration directly checks both directions. Thus the
48 possibilities are the six orders of three distinct nonzero linear
parts and the eight constant choices.

Use these domain-preserving affine transformations: permute the three
legal pair states (any such permutation extends uniquely to an affine
bijection of F_2^2), apply any invertible linear map M to (u,z), and
add any affine function T(p,hp) to (u,z). Its values at the three legal
pair states can be any three vectors in F_2^2, since those points are
affinely independent. A global output complement is also allowed.
The audit enumerates all 6*6*4^3*2=4608 transformations. The orbit of
core code 238 is exactly the 48 valid core colorings, so C0 may be fixed
to that code. The 96 transformations fixing this code form its stabilizer.
These operations act affinely on the original input coordinates,
preserve D, and preserve the factorization through G.

**Four possibilities for the second core.** Let C1 be the layer
G=1,hx=hy=0. If two different t points have C0=C1 with the same color,
the full (x,y) plane times their two-point affine line is a forbidden
three-flat. Thus C1 may agree with C0 at most once of each color.
There are 49 possibilities: complement C0, optionally toggling one
position where C0 is zero and one where it is one. The 96-element
stabilizer acts on these possibilities. Exact orbit enumeration gives
four representatives 1809,1811,1813,3857, of orbit sizes 12,24,12,1.
The audit regenerates and checks this partition rather than assuming
that a search or an external classification found every case.

**Residual normalization.** Name the remaining layers

    A=(G=0,hx=1,hy=0), B=(0,0,1),
    D1=(1,1,0), E1=(1,0,1), F1=(1,1,1).

(The layer names D1,E1 are local to this paragraph, not gate-step types.)
In each fixed C1 case retain the core stabilizer elements also fixing C1.
In addition, the affine shear (u,z)->(u,z)+R*hx+S*hy, for arbitrary
R,S in F_2^2, preserves the domain, both cores and the G factorization.
It translates A,D1 by R, B,E1 by S and F1 by R+S. Finally swapping
x,y and hx,hy exchanges A,B and D1,E1 and preserves G.

Choose a coloring with lexicographically minimum ordered pair of integer
layer codes (A,B) in its orbit under this entire residual group. This
joint choice justifies all of the following conditions simultaneously;
it does not assume that successive greedy normalizations commute:

* A is minimum under the C1 core stabilizer composed with a uniform
  free-coordinate translation R. Its numbers of possible codes in the
  four cases are 270,364,262,49.
* B is minimum under the four translations S, giving 1072 possible codes.
* A<=B, because the layer swap is in the residual group.

The audit computes these canonical codes by complete finite group orbits
on all 4096 layer colorings. The groups are explicit affine subgroups:
the core stabilizers are subgroups, and uniform free translations are
normalized by the invertible free linear maps, so their compositions
are closed. Each of these A transformations is realized by a global
transformation in the residual group. The independent S action fixes A;
therefore the lexicographic minimum has the asserted B condition too.

**Reconstruct the normalized CNF.** For each of the four C1 codes,
append the C0 and C1 fixed colors to the projected-flat clauses. The
canonical A and B sets are encoded exactly by minimal forbidden binary
prefixes: descend the twelve bits from most to least significant, and
for each prefix with no allowed completion append the clause excluding
that prefix. Fully allowed subcubes need no clause. This is precisely
membership in the computed finite sets.

Encode A<=B by eleven fresh prefix flags, IDs85..95. The first prefix
is true. At a compared bit a,b, with current flag e, include
NOT e OR NOT a OR b; omit NOT e for the first bit. If another bit
follows, the next flag e' is forced by both clauses
NOT e OR a OR b OR e' and NOT e OR NOT a OR NOT b OR e'. These force
continuation when preceding bits are equal. They exclude the first
1>0 comparison and admit every A<=B assignment by choosing the flags
to be the actual equal-prefix indicators. Thus this auxiliary encoding
adds exactly the intended condition, up to existential flags.

Apply exact ordinary unit propagation and constant simplification until
no unit remains, rejecting contradictions if encountered. Normalize each
remaining clause by variable index and sign, remove duplicates and sort
by length then literal tuple. The resulting four CNFs have respectively
5130,5049,4773,4408 clauses on at most 95 variables. The unit counts are
24,24,25,26. These transformations preserve satisfiability; any original
counterexample has a normalized representative in one of these cases.
The audit regenerates every clause from the domain and normalizations.

**Generate and check exact resolution.** The audit regenerates all
four normalized CNFs and performs deterministic finite clause-learning
search, with100000 conflicts allowed per case. A cap, error, model or
incomplete run fails; no external solver or stored refutation is used.

Every stored clause is a regenerated root or a derived consequence.
Analysis starts from a stored conflict clause C and a reason clause D
containing a trail literal l, with -l in C. The code checks both pivot
memberships and the valid stored reason index, and computes EXACTLY
(C minus{-l}) union(D minus{l}). This resolvent follows from both parents:
falsifying its remaining literals would force both l and -l. Thus every
intermediate and learned clause is entailed by roots, by induction.
Consistency/falsification assertions add checks but are not needed for
this entailment argument.

At positive levels it resolves backwards to one current-level literal,
stores that entailed clause, backtracks and propagates. Search weights,
saved phases and fixed restarts select inferences only. At level zero,
every trail assignment has a recorded reason; reverse-trail resolution
must produce the EMPTY set, explicitly asserted before UNSAT is returned.
The other two UNSAT paths are a root empty clause or two contradictory
root units, whose single resolution is checked. Hence every UNSAT return
has a sound empty-clause derivation, irrespective of heuristic search
completeness. The wrapper requires such a return for ALL four actual
cases. Their complete normalization above proves Lemma R.



### 5.3 The read-once seven-coordinate tool Q

**Essential-coordinate fact.** On either (two pairs,three free) or
(three pairs,one free), a Boolean output with no legal constant three-flat
must use every one of the seven coordinates. If a paired coordinate is
ignored, set its partner zero, making the ignored coordinate free. The
remaining domains are respectively24 points in F_2^5 or18 points in F_2^5.
If a free coordinate is ignored, the remaining domains are36 points in
F_2^6 or27 points in F_2^6. Majority colors yield respectively12,9,18,14
points, with pair counts66>31,36>31,153>63,91>63. Equal pair sums give
a constant plane. Adjoin the ignored free direction to obtain a constant
three-flat. This proves the fact without computation.

**Lemma Q.** Every six-gate read-once full-B2 formula on seven coordinates,
two promised pairs and three free, has a legal constant three-flat.
A degenerate gate would ignore a leaf and is covered by the preceding
fact. Otherwise normalize to AND/OR/XOR gates and signed leaves by
De Morgan and XOR identities along the unique paths. Root complementation
does not affect monochromatic flats. Translate the three free coordinates
to remove their signs; retain all four paired signs. These operations
preserve affine flats and the legal domain, and cover every original
nondegenerate B2 truth table without adding gates.

Generate all10395 unordered seven-leaf trees recursively, taking the
least label on the left of each proper split. Within-pair swaps, pair
exchange and permutations of the three free labels form48 symmetries
preserving the domain. First-encounter orbit marking yields350 shapes:
each representative marks all its48 images, and any still-unmarked tree
starts a new orbit. This covers all trees; the reported count is checked,
not used as an unsupported classification assumption.

For each shape, evaluate ALL16*3^6=11664 paired-sign/gate assignments
as exact integer bit masks on the72 legal points. Across350 shapes
this checks4082400 assignments. The same affine-flat generator yields
all5049 contained three-flats. No precomputed common-flat hints are used.
For each shape start with every assignment uncovered. Choose its least
uncovered assignment bit, scan the generated flats for a monochromatic
one at that bit, and compute the FULL assignment mask constant on that
flat by intersecting its eight point masks, and likewise their complements.
The mask must contain the selected bit. Remove exactly that full mask
from the uncovered set and repeat, requiring an empty remainder. A reused
flat cannot contain an uncovered bit after its full mask was removed;
the code checks uniqueness and the finite flat-count bound. A missing
witness, error or incomplete run fails. Successful complete coverage
proves the lemma for every normalized assignment, and inverse normalizations
recover every original formula.

For both finite tools, the highest-pivot RREF generator is exhaustive.
It chooses increasing pivots; each direction has its pivot1, zeros at
other pivots and arbitrary lower nonpivot entries, and the offset is zero
at all pivots. These are unique reduced forms for affine flats. Extending
the affine span visits all its points. A contained final flat has contained
prefix spans, so pruning illegal points loses none. The canonical offset
of a legal flat is one of its points and is visited. Thus the14391 and5049
counts are checks of a complete normal-form enumeration, not the reason
for believing that enumeration complete.

## 6. Additive seven from Ramsey R and read-once Q

**Theorem.** For d>=3, a total Boolean f nonconstant on every dimension-d
affine subspace requires at least4n-8d+7 gates in the same strengthened
full-B2 model. For n<=d the right-hand side is negative, so assume n>d.
Run Sections1–4's exact reduction, stopping the first time a=d.

With E=N0+NB+N1+N2 and ordinary-deletion count K, let Lambda>=0 be the
integer sum of additional deletions beyond ALL complete J=s-2q bounds.
The exact telescoping form of equation(8), for R_ex=|C|-(4n-8d), is

 R_ex=4(d-E)-2q_T+s_T+2N0+3NB+N1+K+Lambda.                 (12)

The minima for used A/C,D,E1,E2 are respectively2,2,3,4 objects.
Affine moves/exchanges preserve J. In particular Lambda=0 forces every
whole step to attain its minimum. Normalize initially to a no-larger
strengthened circuit if necessary; original extra gates only help.
Suppose R_ex<=6.

**Terminal q_T>=3.** Since E+q_T<=d, (12) is at least2q_T. Thus q_T>=4
is impossible, while q_T=3 forces E=d-3,s_T=0 and all other terms zero.
All unpaired coordinates are inert. A constant or paired-coordinate output
literal is constant on a product of three legal pair lines: for a literal,
fix that coordinate zero and vary its partner. Adjoin E inert directions
to get a forbidden dimension-d space. No equation or gate remains to
allow any other output form.

**Terminal q_T=2.** If E<=d-3, (12) gives at least8. Otherwise E=d-2,
and R_ex<=6 gives s_T<=2. The four paired coordinates must all be used:
if one is ignored, fix its partner zero and vary it freely; two equal
outputs on the other pair give a constant line, hence a constant plane.
Adjoin E inert directions for a contradiction. But at most two objects
can depend on at most three parameters. Indeed at most two affine rows
can be oriented as two counted XOR gates: a nonzero determinant term
makes their diagonal entries1, and invertibility forbids both off-diagonal
entries being1. Thus their dependency is acyclic, with at most two parents
per row. With any remaining acyclic gates, the entire two-object output
cone still has at most three parameter inputs. Constants and signs do
not increase this count. This contradicts all four being used.

**Terminal q_T=1.** Section3 gives E<=d-2; E<=d-3 would make (12)>=10.
Thus E=d-2 and every other term is zero, including s_T,K,Lambda,N0,NB,N1.
There are N2=d-2>=1 E2 events. Immediately after the LAST E2, the
four-unit fact gives a<=d+1, so a=d or d+1.

If a=d, the reduction has stopped with q=1. Immediately before E2,
a=d+1,q=3,E=d-3. There are three pairs and one non-inert free coordinate.
Exactly four objects remain and all four are acyclic gates, because the
exact E2 count already removes four acyclic gates and leaves s_T=0.
The essential-coordinate fact requires all seven active inputs. A connected
fan-in2 cone with four gates can use at most five input coordinates,
a contradiction.

If after E2 one has a=d+1, the three-unit fact gives q<=1. If q=1,
its predecessor has a=d+2,q=3,E=d-3: three pairs and two free coordinates.
The E2 minimal AND reads paired x,y from distinct pairs. Paired cleanup
has been exhausted, so x,y occur only in G and in no affine row or output
reference. The output is arbitrary H(G,hx,hy,p,hp,u,z). Lemma R gives a
constant three-flat, which together with E inert directions is forbidden.

If instead q=0 after E2, the only remaining paid step is one D, reaching
a=d and q_T=1. No paired cleanup can occur while q=0 or after the immediate
stop. All six objects before E2 must be acyclic: the exact E2 and D
steps remove respectively four and two acyclic objects, rows never turn
into acyclic gates, and K=Lambda=s_T=0 leaves no other deletion. There
are two pairs and three non-inert free coordinates. All seven must be
used. Seven input uses and five non-output gate uses exhaust all12 pins;
thus every input and internal output is used once, giving a read-once
six-gate formula. Lemma Q supplies the forbidden constant three-flat.
This exhausts q_T=1.

**Terminal q_T=0.** If E<=d-2, (12) gives at least8. Hence E=d-1.
Section4.3 proves the final counted inert event is the final paid E1,
so N1>=1. The bound R_ex<=6 forces N0=NB=0 and N1<=2. Since d>=3,
there is a penultimate counted E-increment event, either E1 or E2.
Immediately after it, E=d-2 and a>=d+1. The four-unit and three-unit
facts force a=d+1,q<=1. The future final E1 requires a pair; no D can
occur before it without reaching a=d and stopping, and no cleanup can
remove the sole pair. Therefore q=1 after this penultimate event.

If that event is E2, its predecessor has a=d+2,q=3,E=d-3. As above,
its paired x,y are read only by G and Lemma R gives the contradiction.
If it is E1, then N1=2 and all remaining budget terms s_T,K,Lambda are
zero. Its predecessor has a=d+2,q=2,E=d-3 and six acyclic objects, all
removed in the two exact three-acyclic-object E1 steps. No affine row
could disappear in these minima or turn into a gate. All seven active
coordinates are used, so the same12-pin equality gives a read-once
six-gate formula. Lemma Q again supplies the forbidden three-flat.

All terminal cases contradict R_ex<=6. Integrality proves R_ex>=7.
Only the four normalized R refutations and the complete read-once Q
coverage are computational premises; every normalization, sharing case,
row count, inert-coordinate lift and terminal-history reduction is supplied.
