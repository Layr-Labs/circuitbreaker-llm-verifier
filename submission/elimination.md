# Part I: Exact reduction and finite Ramsey tools

# Full B2 gate lower bound for affine dispersers

## 1. Exact unrestricted reduction

Use the supplied Mishra--Santhanam-style affine promised-pair reduction.
We prove its representation, legality, deletion counts and termination;
obsolete weaker-potential corollaries are omitted. In the final theorem
d>=4 and n>d; n<=d gives a nonpositive target immediately.

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

Literal reroutings target only parents of removed objects; they never
introduce a formerly unused parameter. Track actual charged objects,
pairs and effective dimension; no degree penalty or extra free gate is
used. All exchanges preserve the computational term count.

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

The following deletions preserve m,q,a and remove at least one charged object.
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
   terminal argument. This removes one charged row and preserves m,q,a.

Affine gates whose inputs are constants, coordinates or affine-system
values can be moved into the system as one new defining row. The new
auxiliary matrix is block triangular with diagonal entry 1, and the row
has at most three terms. This exchanges one gate for one equation, hence
preserves s. It does not change whether any parameter is used.

For a general affine exchange Q=x+ell with unpaired x, both x before
and Q after occur in the system: Nx!=0, and the old Q column of A is
nonzero. Swapping roles changes no row or gate terms, so s,m,q,a stay fixed.
This exchange is valid regardless of old x's degree. Only Phi records
the affine substitution; no dense expression or free gate is inserted
into the computational representation.

### 1.3 Paired cleanup, before choosing a minimal AND

Apply deletions and move eligible affine gates, then inspect these cases
in order. Each paired cleanup has Delta a=0 and strictly
lowers m+s. Rerouting may remove additional objects; counts below are
lower bounds for distinct original objects, not pin counts.

(A) Paired p has partner h. If p is unused, impose h=0 and erase h and
the pair. Surviving p becomes unpaired but unused.
If p has at least two distinct children, impose p=0 and erase p and
its pair. Delete two affected children. Acyclic children have a fixed
pin; row children become short and use deletion 3. If propagation removes
one early it is already a counted deletion; remaining affected rows
stay short. The only newly unpaired coordinate is h. After (A) fails, every paired
coordinate has exactly one child.

(B) Partners p,h share their one child H. If H is a gate, evaluate its
full two-input Boolean table on (0,0),(0,1),(1,0). Two outputs agree.
The line through the two corresponding points is one of p=0, h=0,
p=h+1. Restrict to that line, making H constant, then delete H. If H
is a row it has the form v+p+h=c, since it has three terms and some
auxiliary. Impose p=h+1 and delete v=c+1. The surviving coordinate
is now unused. All rerouting is to constants. Exactly one object is guaranteed removed; all rerouting is to constants.
After (B) fails, partners have different children.

(C) Paired p occurs in a row Hp, and partner h has unique child Hh!=Hp.
Impose h=0, remove h and the pair, and delete Hh. Parameter p now is
unpaired and still has at most the one row occurrence Hp: Hh did not
contain p. If Hp already disappeared by propagation, it is the second
counted deletion; at most one coordinate becomes unpaired.
Otherwise delete a short Hp,
or use deletion 4 when it has three terms. This removes a second object
and leaves at most one newly unpaired coordinate. Other unused
coordinates cannot become used.
After (C) fails, all paired coordinates are isolated, occur only in the
acyclic part, and have a different child from their partner.

### 1.4 Exhaustive minimal-AND step

If a>d, some acyclic gate remains: otherwise the output is affine.
Choose the first one in topological order, G. It is AND-type because
eligible affine gates were moved into the system. Its inputs are affine
values with independent linear parts, or deletion 1 would apply.
If a parent Q is auxiliary, its solved affine form depends on some
parameter x. Such x occurs in a row and is unpaired by absence of (C).
Exchange Q,x preserving s,m,q,a. Check deletions and paired cases again;
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
once; G and its child are two distinct original objects. No old unused coordinate becomes used. Thus

 Delta m=0, Delta q=1, Delta a=-1, Delta s<=-2.

(E1) Exactly one input is paired; write x unpaired, y paired with h.
The unique child Hh of h differs from G by (B). Impose h=0 (erase
its pair), then fix unpaired x=alpha to force G constant. Delete G
and Hh. These are distinct, even if Hh is a child of G or they have
shared descendants. Coordinate y is unused afterwards, because G was
its only child. No old unused parameter becomes used. Hence

 Delta m=-2, Delta q=-1, Delta a=-1, Delta s<=-2.

(E2) Both inputs are paired. Their pairs are distinct; if they were
partners, they would share G and (B) would apply. Let hx,hy be their
respective partners. Both have unique children different from G.
Impose hx=hy=0 and erase both pairs. Fix now unpaired x=alpha,
forcing G constant. Delete G and the unique child Hx of hx. They
are distinct; hy's child may equal Hx, and is not counted again.
Finally x is removed and y is unused. Thus

 Delta m=-3, Delta q=-2, Delta a=-1, Delta s<=-2.

No nonconstant rerouting introduces an unused parameter. Cases (D),(E1)
and(E2) all lower a by exactly1; cleanup/deletion leaves a fixed.

### 1.5 Termination

While a>d, terminal forms contradict dispersion, or one of the exhaustive
cleanup/deletion/AND cases applies. Preparations use finitely many
 topological affine moves and at most two exchanges; restart only upon a
strict deletion or cleanup. Every successful step strictly lowers m+s,
and neither ever increases. Thus a disperser reaches EXACTLY a=d after
finitely many legal steps. All objects remain charged original gates or
one-for-one rows. Sharing, all B2 tables and final literals are covered;
no direct-sum or computational-basis restriction is used.

## 2. Count the constant closure

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

## 3. Permanently inert coordinates and ledger

Use exactly the terminating reduction already proved. Each step changes a by0 or-1, so the terminal value is exactly a_T=d.

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
For q_T=1, the remaining possibility E=d-1 is impossible. There are
exactly d-1 unpaired coordinates (since a_T equals the number of unpaired
coordinates plus q_T); they would all be inert. Hence the output depends
only on the two coordinates p,h in the sole pair, regardless of how many
gates or rows remain. On the three legal points (0,0),(0,1),(1,0), two
output values agree. Their affine line is one of p=0,h=0,p=h+1, and lies
inside the promised domain. Varying that line and all d-1 inert coordinates
gives a constant-output affine space of dimension d, a contradiction.
Thus q_T=1 implies E<=d-2.

## 4. Refined E1 count and terminal timing

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
whole-step inertness is unchanged. In particular
case (E1) has Delta J=Delta(s-2q)<=-3+2=-1, rather than merely <=0.

All additional object losses are nonnegative integer slack. Part II
counts ordinary deletions by K and all such excess by Lambda; affine
moves/exchanges preserve J. Inertness and E bounds remain as proved.

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

## 5. Terminal geometry and an unrestricted Ramsey lemma

The following geometric facts and Lemma R are used in the +8
refinement proved in Part II below. Its finite refutation is regenerated
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

# Appendix B: Six-variable Ramsey lemma

## B.1 Reduce every six-variable coloring to explicit finite cases

**Lemma S.** Every coloring f:F_2^6->{0,1} has a monochromatic affine
three-flat. Suppose otherwise. Each color class affinely spans F_2^6:
if one were contained in a proper affine hyperplane, the complementary
parallel hyperplane would lie entirely in the other color and contain
a three-flat.

Some three-flat has at least six points of one color. Indeed choose a
color with at least32 points. Its unordered pairs outnumber the63
nonzero differences. Two pairs with the same difference are disjoint
and form a monochromatic affine plane P. Its direction has16 cosets.
If each of the other15 cosets had at most one point of this color,
there would be at most4+15=19 such points. Some other coset therefore
has at least two. Its union with P is a three-flat with at least six.
Choose a three-flat with the greatest monochromatic count in either
color and rename that color zero. The count is6 or7, since8 would
already contradict the assumed counterexample.

If the maximum is6, an invertible affine map makes the six zeros in
this three-flat be the points0,1,2,3,4,5 and the two ones6,7, using
integer coordinate order with the first three coordinate bits varying.
An affine group on three bits is transitive on pairs of distinct points,
which proves this normalization. The zeros span all six dimensions,
so extend the first three basis directions by choosing zero points
outside the successive spans. This makes8,16,32 zero as well.
Every three-flat has at least two points of each color, by maximality.
For each of its eight points omitted, the other seven are therefore
not all zero and not all one. These are the 178560 regenerated
seven-literal clauses. The eight core units and three basis units give
178571 clauses on64 color variables. This is the first refuted case.

Now suppose the maximum is7. Make its three-flat the points0..7,
with its unique one at7. Choose a zero point outside this span and
make it8. In the coset8+F_2^3 the zero set contains no affine plane:
a zero plane there and the parallel plane in the original core avoiding
point7 would together form a zero three-flat. Any five points in F_2^3
contain a plane, since their ten pairs outnumber the seven nonzero
differences. Hence the coset has between one and four zeros. A four-set
without a plane is affinely independent. Affine permutations of F_2^3
are transitive on one-, two-, three-point sets and independent four-sets.

These permutations are realizable while preserving the original core.
For an invertible linear A on three bits, the map u->Au+b with
b=7+A7 fixes its unique one. An independent shear by coordinate bit3
adds any chosen vector on the coset8+F_2^3 and nothing on the core.
Thus the coset zero set may be exactly one of

    {0}, {0,1}, {0,1,2}, {0,1,2,4}

in its three-bit local coordinates. These are four cases, not an
assumption about an arbitrary coloring. They fix all first16 colors.
The canonical patterns include zeros0,1,2,4,8. The full-span property
allows a further zero outside this four-space to be chosen as point16.

Enumerate all65536 colorings of the coset16+F_2^4. Together with the
fixed first16 colors, require that each of the620 affine three-flats
in F_2^5 be mixed. Exact enumeration leaves respectively
3136,5760,6232,5296 patterns. It uses all five-bit three-flats, not just
coordinate flats. Truth tables are integer bit vectors and clause
satisfaction is evaluated exactly.

Use the following explicit subgroup to reduce these patterns. Enumerate
all168 invertible three-bit linear maps A, put b=7+A7, and all eight
shears R. On four bits (u,h) use

    (u,h) -> (Au+b+Rh,h).

Retain those preserving the fixed first16 colors. Their numbers in the
four cases are168,48,24,24. On the coset16+F_2^4, add an independent
translation T in F_2^4. It is the shear by coordinate bit4 in the full
space and fixes the first16 points. These maps form a group: the
invertible block affine maps are closed under composition, the
color-preserving subset is a stabilizer, and the independent shears
form its invariant translation subgroup. Every map is invertible and
preserves monochromatic affine subspaces.

The auditor generates all these permutations. For each still-unseen
pattern, it computes all its images, verifies they remain in the
allowed set, and chooses the least numerical image whose bit0 is zero.
Every allowed pattern has a zero, since an all-one coset would itself
contain a monochromatic three-flat. Translation makes bit0 zero, so
this minimum always exists. The resulting orbit counts are6,18,28,30,
with each original allowed pattern included. These computations give
82 explicitly regenerated representatives; no externally supplied
coloring-class table is assumed.

Only after this normalization choose a zero outside the first five-space
as point32. The zero class still spans six dimensions. The independent
shear by the last coordinate can make that point32 while preserving all
first32 colors. This ordering matters: no claim is made that the earlier
normalizations preserve a previously fixed point32.

For each of the82 representatives, fix its first32 colors and assign a
Boolean variable v_(a-31) to each point a=32..63. Regenerate all11160
three-flats in F_2^6. Each gives the two clauses requiring at least one
point of each color. Substitute the fixed colors, drop satisfied clauses
and fixed false literals, remove identical literal tuples, and sort by
length then tuple. Include the unit -v_1 for the chosen zero at32.
A counterexample would satisfy one of these82 residual CNFs, or the
maximum-six CNF above. This proves complete case coverage.

The affine-flat generator is the same exhaustive highest-pivot RREF
construction as Part I Section5.2, applied now to the full cubes. Increasing
pivots and offsets zero in pivot positions give each flat once; all
points are generated by extending the current affine span. The counts
620 and11160 are independently checked by the auditor.

## B.2 Direct checked refutations and the d=3 case

The audit regenerates all83 CNFs above and runs the exact checked
resolution algorithm proved in Section5.2, requiring a derived EMPTY
clause in every case. A model, error, conflict cap or incomplete search
cannot pass. These complete actual derivations prove Lemma S; no stored
CNF, RUP payload or SAT verdict is trusted.

For d=3, no eligible function exists at n>=6, since restricting to any
six-coordinate cube and applying S gives a constant three-flat. For
n<=4,4n-16<=0. At n=5 every input is essential: ignoring one leaves a
four-cube, whose majority's28 pairs exceed15 differences and force a
constant plane; adjoining the ignored direction gives a three-flat.
Five essential inputs need at least4 gates by the proved connectivity
bound. Thus C>=4n-8d+8 also holds at d=3. No sharpness claim is needed.


# Part II: Additive8 for d>=3

**Theorem.** In the pinned full-B2 model, every total Boolean f on n
bits nonconstant on every affine d-space, d>=3, satisfies C(f)>=4n-8d+8. Appendix B handles d=3;
assume d>=4 below.
Use Part I Sections1–5.2's exact promised-domain reduction. Write a=m-q,
s for charged objects, E=N0+NB+N1+N2 for distinct permanently inert
coordinates, K for ordinary deletions and Lambda for the nonnegative
integer slack in ALL complete J=s-2q step bounds. Stop at a=d.
Since n-d=N_D+N1+N2 and E=N0+NB+N1+N2, summing the complete
J bounds with refined E1 change -1 gives the EXACT excess

 R=|C|-(4n-8d)=4(d-E)-2q_T+s_T+2N0+3NB+N1+K+Lambda.       (11)

Zero slack in used A/C,D,E1,E2 means exactly2,2,3,4 distinct removed
objects, respectively, including propagation. Affine moves/exchanges
preserve J; no nonnegativity of J is assumed. For n<=d the target is
nonpositive. Otherwise suppose R<=7 and contradict every terminal case.

## 1. Geometry and object-count facts

The proved pair-sum facts and timing give E<=d-1,E+q_T<=d,
E=d-2=>a<=d+1, E=d-2,a=d+1=>q<=1, and first reaching E=d-1
is the final paid E1. Lemma R applies to three pairs plus two free
coordinates whenever the distinguished AND input coordinates have no
other occurrence; no condition on the remaining function H is needed.
Appendix B's Lemma S gives the additional fact

 E=d-3 => a<=d+2.                                      (12)

Otherwise choose a legal line in each pair and fix surplus free
coordinates to obtain a six-cube; its constant three-flat, together with
E inert directions, contradicts dispersion. A computationally unused
unpaired coordinate may be counted as an additional inert direction.

No function without a legal constant three-flat on the domains
(one pair,four free), (two pairs,three free), (three pairs,one free)
can ignore ANY displayed coordinate. If a paired coordinate is ignored,
fix its partner zero; the remaining domain sizes/dimensions are
16/4,24/5,18/5. If a free coordinate is ignored they are24/5,36/6,27/6.
The majority's pair counts exceed nonzero differences respectively by
28>15,66>31,36>31,66>31,153>63,91>63. A constant plane plus the
ignored direction is a forbidden three-flat. These are essential-input
statements on the legal domain itself.

T essential parameters require at least T-1 charged objects, including
fair affine rows. A nonzero determinant term matches each row to its
own auxiliary, making it an output with at most two remaining parent
pins. This orientation may be cyclic but adds no free operation. Every
essential parameter belongs to the undirected output component: other
components' auxiliary blocks are independently invertible. With s'
objects and at least T parameters, connectivity needs s'+T-1 edges,
while the2s' parent pins supply at most2s' edges. Hence s>=s'>=T-1.
For acyclic circuits, dead-gate removal gives the same count for all
used inputs. At most two rows can be replaced by at most two counted
XOR gates: determinant matching makes their diagonal entries1 and
invertibility excludes both off-diagonal entries being1. The resulting
one-/two-row dependency is acyclic, with at most two parents per row;
constants/complements are absorbed in full-B2 tables. Every row is paid.

## 2. Finite formula lemmas

**Q.** Every six-gate read-once B2 formula on seven coordinates, two
promised pairs and three free, has a legal constant three-flat. Normalize
to AND/OR/XOR and signed leaves; remove the three free signs by
translations and retain four paired signs. Enumerate all10395 unordered
seven-leaf trees, choosing each proper split with least label on the
left. Within-pair swaps/exchange and free S3 form48 symmetries, leaving350
representatives. Generate72 legal points and all5049 three-flats.
Exact masks cover all16*3^6=11664 assignments per representative,
4082400 total, by Section5's complete checked flat procedure. Inverse
normalizations preserve affine flats and exhaust every original formula.

**T7.** On the same domain, every circuit with<=7 gates, each paired
coordinate occurring once, and a distinguished AND-type G on literal
coordinate inputs has a constant legal three-flat. Remove dead/unary/
constant-pin gates. With<=6 gates all seven coordinates are essential
by Section1, so pin equality gives exactly6 and a read-once tree: Q.
At7 gates, seven inputs plus six non-output uses consume13 of14 pins.
Exactly one vertex B repeats; its cone is read-once and the outside is
a tree using it twice. After contracting G, a six-component cone would
leave a parent reading its same value twice; reduce that gate and use Q.

Some G positions are analytic. With two free inputs, fix a forcing value,
choosing the repeated one if needed; the other is unused, and the two
pairs/one-free residual18-point domain forces a plane. For partner
inputs choose their constant legal line; its line parameter is unused,
and the other pair/three-free24-point domain forces a plane. For inputs
from distinct pairs, alpha=0 in G=(x+alpha)(y+beta)+gamma allows x=0,
hy=0, leaving y unused and four free coordinates forcing a plane;
beta=0 is symmetric. For free x and paired y, normalize x's sign by
translation. If beta=0 and x is not repeated, y=0 makes x unused and
the remaining pair/three-free domain forces a plane.

Contract G to component0. The remaining cross case is its NOR forcing
orientation, components1,2 its partners and3,4,5 free. The mixed case
has G=x*(y+1), partner1, other pair2,3 and free4,5. A repeated vertex
is free, G, or an internal cone on these components; paired components
cannot repeat. Enumerate every subset of size1..5, every inside tree
and every outside tree on the remaining components and two cone copies6,7.
For mixed G, also enumerate repeated original x outside G, all seven-
component outside trees and both paired forcing orientations. This
exhausts single-extra-use graphs. The cross symmetry is partner exchange
and free S3 (order12); mixed swaps the other pair and two free leaves
(order4). Choose least subset in each orbit, then canonicalize trees by
its set stabilizer; this exactly quotients the whole group. Repeated-x
is canonicalized directly. The regenerated counts are4208 and12468.

Shared-cone signs require care. If B contains G, push G's output sign
along its unique inside path to B, using De Morgan/XOR and off-path
negations; both B occurrence signs toggle. Push the first B occurrence
sign through the outside path to the global output; off-path signs of
the other occurrence and leaves are still enumerated or free translations.
If B excludes G, first fix G's sign in the outside tree with B as signed
atoms, then complement the inside B tree to fix its first occurrence;
it contains no G. For input B, use the same outside path when B=G,
or translate its free coordinate. Thus only the RELATIVE second-B sign
remains beyond paired signs; no two independent cone signs are assumed.
Cross masks have3 sign bits and6 trits; mixed have4 and6. Repeated-x
mixed has5 sign bits, retaining its extra x sign because translation
already normalized G's input, and both G orientations. The checked totals
are24541056 and181643472 assignments. Generate72 points/5049 flats;
mixed uses permutation(1,3,2,4,0,5,6) of cross coordinates. Values are
cached only for identical sign layouts, operation positions and inputs.
Section5's coverage proves every assignment, establishing T7.

**T8.** Every seven-gate read-once formula on eight coordinates in four
pairs, with distinguished AND-type G on distinct pairs, has a legal
constant three-flat. A good forcing orientation sets x=0,hy=0, leaving
y unused and an18-point two-pair/one-free residual domain forcing a plane.
For the remaining NOR orientation, contract G and enumerate all seven-
component trees,64 paired signs and6 operation trits. Partner exchange,
swaps within the other pairs and their exchange give16 symmetries,880
representatives and41057280 assignments. Regenerate81 legal points and
all2106 three-flats and check full exact coverage. The same sign/tree
normalization as above, now without sharing, is complete.

## 3. Exceptional tail and predecessors

**F.** A five-gate read-once formula on six inputs, one pair and four
free, either has a legal three-flat or is equivalent to

 OR(x,u) XOR OR(y,v) XOR OR(w,z), x*y=0.                (13)

Allowed equivalences are free complements/permutations, pair exchange
and output complement. Gate-level, it has three disjoint AND/OR cherries
joined by two XOR gates, up to internal complements. To prove this,
normalize AND/OR/XOR, enumerate945 six-leaf trees modulo pair exchange
and free S4 (48 symmetries,41 representatives), four paired signs and
3^5 gate choices. On48 legal points and all1770 three-flats, exact mask
coverage leaves16 assignments on two trees with eight legal truth tables.
The audit checks each exception's three AND/OR cherries/two XOR joins
and its FULL64-point truth table in the32-element orbit of(13) under
four free complements and output complement. Each paired cherry's free
linear coefficient is1, so free-input/cherry-output complements turn it
into OR; likewise the free/free cherry. XOR joins absorb output signs.
Both tree associations are retained in the following checks. Thus this
is gate-level normalization, not substitution of a merely promise-equal
tail into a different surrounding circuit.

**A predecessor.** Restoring p paired with a free tail input h before
a tight used-A cleanup p=0 deleting two objects always gives a legal
three-flat. If both objects are gates, p has exactly two children; prune
its occurrences to obtain the tail, or inversely graft two p copies on
arbitrary tail edges. Normalize forced outputs by parent compensation.
Old-containing grafts have tables4,6,12,14 (identity at p=0); a new-only
p,p subtree needs zero-preserving tables2,4,6,8,14, giving both0 and p
on equal inputs. Identity padding12 covers a simplifiable second graft.

There are at most two extra rows, with triangular matched system.
Independent Q=p+r,R=p+s give XOR grafts. Triangular Q=p+R,R=p+r
reduces to Q=r,R=p+r: both cannot be used independently in the tail,
or at p=0 an old input repeats, violating six-input/five-gate pin equality.
At most one effective graft plus padding suffices. With Q=p+r and a
p-gate reading Q, the used result is a Boolean function of p,r, again
one graft plus padding; reading another old value gives two grafts.
No gate feeds a row. Constants/signs are absorbed. This exhausts all
zero-/one-/two-row possibilities; unused outputs only simplify them.
After restoring the pair, h cannot be translated. Retain an extra h
leaf-sign bit and all16 tail exceptions; other free complements remain
legal. Four h choices, both associations and77 two-graft structures give
616 cases/163328 assignments, all checked by exact three-flat coverage.

**C predecessor.** Reverse tight paired-C: p is read only in Hp, h only
in Hh!=Hp, and h=0 removes both with the allowed role exchange. The
surviving tail has no rows. If Hh is a gate, Hp solves Q=p+r+c, r an
old parameter or0. At h=0, Hh must reroute to a nonconstant old value,
else propagation deletes a tail gate; its output has only one tail use,
else that value repeats. Thus Hh is one h-identity-at-zero graft at
any edge, including Q, and the old Q leaf is p+r+c.

If both are rows, both cannot contain both auxiliaries (singular2x2
matrix). Independent Q=p+r,R=h+t is an h-XOR graft. For R=h+r,Q=p+R,
if R is unused use that model; otherwise its unique use is the old r leaf,
so old r->r+h and old Q->p+r+h. For Q=p+r,R=h+Q, at most one of Q,R
has a tail use after h=0 by pin equality, giving a single graft or an
ineffective h (already yielding a plane plus unused direction). These
exhaust invertible two-row cases. Absorb constants into old paired signs,
a Q leaf sign and other legal free complements; never translate p.
The audit takes four Q slots, two associations, r among five old inputs
or0, a Q sign and tables4,6,12,14 at eleven edges, plus the joint row
model for each nonzero r. Old paired signs affect tail leaves, NOT the
row formula. All96 cases/34432 assignments have checked three-flats.

**E2 predecessor.** A tight E2 deleting four gates from a nine-gate
read-once formula on ten inputs, three pairs/four free, to yield F has
a legal FOUR-flat. The removed/newly inert inputs are x,y,hx,hy.
Contract G to A and prune A,hx,hy. No old-containing branch can die,
since all six old inputs remain essential; it reroutes to the branch
or complement. Conversely graft A,hx,hy in that fixed order at arbitrary
edges of either tail association; reverse pruning proves completeness.

Write G=(x+alpha)(y+beta)+gamma, absorb gamma, and force x=alpha,
hx=hy=0, so A=0. Normalize every new-only subtree's forced value to0,
compensating its parent; its gates have tables2,4,6,8,14. Old-containing
grafts are identity when the new child is0 and have tables4,6,14.
Tail normalization is safe at gate level: conjugating an intervening
D(a,b) to1+D(1+a,b) preserves identity-at-b=0 and permutes those tables.
Old free complements/permutations are legal before E2. All four G
orientations, new tables and graft positions are enumerated. Old-pair/
cherry-partner exchange, the other free exchange and hx,hy exchange
form8 symmetries, leaving1450 structures/188616 assignments. The audit
validates360 supplied flats by offset/four basis vectors,16 distinct
legal points and uniqueness, then checks full mask coverage. They need
not exhaust all four-flats; a valid covering subset proves existence.

## 4. Exclude the exceptional five-gate state by its history

There is no state in the reduction with

    E=d-3, a=d+2, q=1, s=5,

which has five acyclic gates, no earlier N0,NB,N1, no ordinary deletion
or slack, and is followed by steps deleting exactly these five acyclic
gates. Its six active coordinates form one pair and four free. They
are all essential by Section1, hence the five-gate pin equality makes
the state read-once. If it lacks a constant three-flat it is an F tail.

Consider the last object-removing step before this state. Zero-cost
affine preparations cannot hide a later transition: the state has no
rows, and moving an affine gate into the system or exchanging roles
cannot remove a row. Any subsequent row removal would itself be a later
object-removing step. Complements are already included in the signs.
There is no earlier E1,N0 or B by the stated event counts, and no
ordinary loss or extra propagation. Thus this last step is used A,C,D
or E2. If A or C, the predecessor facts give a constant three-flat,
which with E inert directions is forbidden. If D, its preceding state
has E=d-3,a=d+3, violating (12).

If E2, its preceding state has E=d-4,a=d+3,q=3,s=9. All nine objects
are acyclic: four guaranteed deletions in E2 and the five subsequent
acyclic deletions exhaust them. Its six paired coordinates are used
uniquely by the cleanup conditions. Its four other active coordinates
are used too: an unused one would give E'=d-3 at a=d+3 and contradict
(12). Ten input uses and eight non-output gates exhaust eighteen pins,
so the state is read-once. The E2 predecessor fact supplies an active
four-flat; with d-4 inert directions this is a forbidden d-flat.
If there was no preceding removing step, the initial state has no pair
and cannot be the displayed state. This completes the history argument.

## 5. Exact finite coverage without hint packets

Every formula statement above is actually checked. The highest-pivot
RREF generator enumerates every contained affine flat once: pivots
increase, other pivot entries and offset pivot bits are zero, every
reduced free entry is visited, and extending the affine span generates
all its points. A contained final flat has contained prefixes, so pruning
illegal points loses none. The extension obligation instead uses360
explicit offset/basis witnesses, each validated as16 distinct legal
points with no duplicate flat; exhaustive four-flat enumeration is not
required for a valid covering subset.

For each regenerated topology, exact integer masks evaluate EVERY
normalized gate/sign assignment. A flat's mask is the intersection of
its point-output masks union the intersection of their complements:
EXACTLY the assignments constant there. Start with the full assignment
mask L. Apply validated common flat IDs as acceleration only. While L
is nonzero, take its least remaining bit b, extract that assignment's
point-color signature, and scan the regenerated flats (the validated
360-flat dictionary for extension). A candidate must be monochromatic
at b and its full mask must actually contain b. Remove its mask from L.
An already used flat cannot cover a remaining bit; check this and the
finite iteration bound. If no flat exists, any assertion fails, or the
runtime ends early, the obligation fails. Acceptance requires L=0 for
EVERY topology. The selected bits guide witness search only; full-mask
coverage exhausts all assignments. No per-formula hint, compressed
refutation, search-status verdict or external solver is trusted.

The sixteen F exceptions are computed/classified from all normalized
tails, not supplied as a trusted list. The two Ramsey obligations use
the separately proved exact resolution derivations. Their normalization
and the universal reductions are the supplied arguments, not inferred
from a successful finite runtime.

## 6. Exhaust all terminal states with R<=7

For q_T=0, E<=d-2 gives R>=8. Otherwise E=d-1,N1>=1 and
s_T+2N0+3NB+N1+K+Lambda<=3, so NB=0.
If N1=1,N0=0, the last E2 reaches E=d-2,a=d+1,q=1 before the
final E1; its predecessor E=d-3,a=d+2,q=3 is excluded by R,
independently of any unused slack.

If N1=N0=1, every other budget term is0. A penultimate E2 is excluded
likewise. Otherwise it is N0, whose predecessor is E=d-3,a=d+1,q=2,s=3.
Zero slack means BOTH members were unused. Such a pair cannot persist
from an earlier completed iteration, since unused-pair cleanup has
priority; affine preparations preserve occurrences. Its members became
unused at the immediately preceding removing step. D gives predecessor
E=d-3,a=d+2,q=1,s=5, the forbidden history of Section4. E2 gives
E=d-4,a=d+2,q=4,s=7, all acyclic. At least one of its two additional
active coordinates is used, or two extra inert directions contradict
the four-unit fact. Nine inputs plus six non-output gates need15 pins,
exceeding14. A/C gives E=d-3,a=d+1,q=3,s=5; all seven active inputs
are essential, requiring at least6 objects. Other step types are
excluded by their positive K,NB or N1 counts.

If N1=2, N0=0 and other terms sum to<=1. A penultimate E2 is excluded
by R. Otherwise the final two inert steps are E1, with no intervening
D or paired cleanup because the last E1 needs the remaining pair.
Their predecessor has a=d+2,E=d-3,q=2 and at most7 objects, at least
6 acyclic. Replace its possible single row by one counted XOR. All
seven inputs are essential, paired ones unique, and the selected
minimal AND has coordinate inputs: T7 applies. If N1=3, the budget
is exact and the corresponding predecessor instead has6 acyclic objects;
use Q. Larger N1 exceeds the budget. Thus q_T=0 is impossible.

For q_T=1, E<=d-2. E<=d-3 gives R>=10. Otherwise
s_T+2N0+3NB+N1+K+Lambda<=1, so N0=NB=0. If N1=0, consider the
last E2. An after-a=d predecessor has q=3, seven essential inputs and
at most5 objects, fewer than the required6. If after-a=d+1, after-q<=1.
After-q=1 gives a predecessor excluded by R. After-q=0 requires a
final D; the E2 predecessor has two pairs, three free coordinates,
at most7 objects, and at most one row beyond the6 guaranteed acyclic
deletions. Replace that row by XOR and apply T7.

If N1=1, all other terms are0,N2=d-3, and E1 is the last inert step.
If its after-a=d, its predecessor has q=2, six active inputs and only
3 gates. The four paired inputs are used; an ignored free input would,
with the constant plane forced on an18-point two-pair/one-free domain,
give a three-flat. Six inputs and two non-output uses exceed6 pins.
If after-a=d+1, after-q<=1. After-q=0 is followed by D and gives the
forbidden five-gate history. After-q=1 must lose that pair by used A/C
before D, else final q=2; unused A or B would make N0 or NB positive.
The predecessor has two pairs, three free inputs and7 objects, from
3 E1,2 A/C and2 D deletions. At most2 rows occur; replace them by
counted XORs and apply T7. These exhaust q_T=1.

For q_T=2, E<=d-2; E<=d-3 gives R>=8. At E=d-2 one has s_T>=3.
Indeed two rows give affine output; one row plus a gate, or two gates,
depend on at most3 of the4 paired inputs. An unused paired input with
partner fixed zero supplies a free direction; on the other pair choose
a constant line. Their plane with E inert directions contradicts
dispersion; zero or one object is simpler. Thus R<=7 forces s_T=3
and all other budget terms0, so every inert event is E2. At the last
one, after-a=d+1 requires after-q=1 before final D; its predecessor
is excluded by R. If after-a=d, the predecessor has q=4,s=7 and all
eight paired inputs used. There are at most3 rows and no parameter
terms in them, since paired-row cleanup is exhausted and unpaired
parameters are inert. Any invertible system of at most3 such rows has
a short row: three full rows would all be111 and singular; smaller
cases are immediate. Short-row deletion was exhausted, so there are
no rows. Seven acyclic gates, eight input uses and six non-output
uses exhaust14 pins, giving a read-once state. Apply T8.

For q_T=3, E<d-3 gives R>=10. At E=d-3, R<=7 permits at most1
object and all unpaired inputs are inert. Its constant/literal/one-gate/
one-row output depends on at most2 paired inputs. If in distinct pairs,
fix them zero along legal lines; if partners, use their constant pair
line. Complete the other pair lines to obtain a constant three-flat,
a contradiction. Finally q_T>=4 gives R>=2q_T+s_T>=8. No terminal
case has R<=7; integrality proves the theorem for d>=4.

