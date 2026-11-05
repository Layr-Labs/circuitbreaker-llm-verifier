# Gate lower bound for an affine disperser

We give the parameter-method argument of Mishra--Santhanam (TR26-260,
Sections 2-4) needed by proof.md, with explicit bookkeeping. An affine
disperser of dimension d is nonconstant on every dimension-d affine
subspace, and hence on every larger affine subspace.

**Lemma.** For such f on n bits, every full-B2 circuit has size at least
(7/2)n-7d.

## 1. Representation and potential

Start with an arbitrary circuit; it need not be a formula or minimum.
Allowing a final output literal for free only strengthens the lemma.
Absorb internal-edge negations in binary truth tables. Retained gates
stay charged until removed. Every nondegenerate two-input table is either
u+v+c or (u+b)(v+c)+e over F_2.

Keep M independent parameters u, a set P of q disjoint parameter pairs,
and an injective map to original input space

    Phi(u)=b+L*u+sum_({p,h} in P) c_(p,h)*p*h.

Coefficients are vectors over F_2; zero pair coefficients are permitted.
The computational representation is an affine system A*v=N*u+c followed
by acyclic B2 gates, with no feedback into the system. The auxiliary matrix
A is square and invertible. Each row has at most three distinct nonconstant
terms. The output is a constant or literal of a parameter, auxiliary value,
or gate. The invariant is that its value is f(Phi(u)) on every assignment.
Let s be number of equations plus acyclic gates, and a=M-q.

Degree counts distinct child rows and gates, not pins or the output
reference. An unpaired parameter has weight 0,6,7 when its degree is
0,1,>=2 respectively. Paired parameters have no individual weight. Set

    H=2s+3q+sum_(unpaired z) weight(z).

Initially M=a=n, q=0, s=|C| and H<=2|C|+7n. Always H>=0.
Fixing one member of every pair makes Phi affine and still injective on
M-q parameters. Its image therefore contains an affine space of dimension
a. Thus f(Phi(u)) is nonconstant whenever a>=d.

The following changes preserve injectivity and the displayed map form:

* Replace an unpaired x by w+ell(other parameters), with ell affine and
  w new. This is a bijection, preserving M,q,a.
* Replace an unpaired z by w+(x+b)(y+c)+e for distinct unpaired x,y,z,
  record pair {x,y}. This is a bijection; M is unchanged, q rises by 1,
  a falls by 1. No old quadratic term used z.
* Fix an unpaired x to an affine function of the others: M and a fall by 1.
  In actual reductions we use only constants or another parameter's literal.
* For pair {p,h}, fix p to a constant or to h+c, then remove p and the pair.
  Its product becomes affine. M,q each fall by 1, so a is unchanged.

The first two changes are not licenses to add free gates. Their exact
implementations appear below.

While a>d, an affine output is impossible. A constant output contradicts
the affine-space invariant. For a nonconstant affine output choose a
parameter with nonzero coefficient; if paired, first fix its partner,
preserving a. Solve for that now unpaired parameter to make the output
constant, losing only 1 in a. The remaining image still contains a
dimension-d affine space, a contradiction. Similarly, an output literal
of an AND on two unpaired parameters is impossible: fix a forcing input,
losing only 1 in a. This explicitly deals with terminal outputs.

## 2. Legal deletions and pivots

An affine gate without an AND ancestor moves into the system by adding
its output as auxiliary and its defining equation as one row. The row has
at most three terms. The new auxiliary matrix is block triangular with
last diagonal entry 1. Size s is unchanged and parameter degrees cannot
increase. Move all such gates in topological order.

Four deletions preserve M,q,a and lower H by at least 1:

(D1) A constant or single-parent-literal gate is rerouted to that value
and removed. This covers unary/constant tables, fixed or repeated pins,
and an AND on affine parents with dependent linear parts. In the last case
the two affine values are constant or equal up to a constant, so the table
has one effective input. At most the one surviving parent can gain weight,
by at most 1 (6 to 7), against a saving of 2. Constant rerouting gives no gain.

(D2) A gate without children is removed unless it is the output. This
saves 2 and cannot increase parameter weight.

(D3) A row with at most two terms contains an auxiliary, since A is
invertible. Write it v=U+c, for U zero or one other value, substitute this
in all occurrences of v and erase the row and v. This preserves sparsity.
Adding the row to all rows using v, then expanding the determinant on
its column, proves the remaining auxiliary matrix invertible. Only U
can gain weight, by at most 1; the row removal saves 2.

(D4) An unpaired degree-one parameter x used only in a row E is pivoted
with an auxiliary Q that depends on x. Such Q exists because N_x is
nonzero and A^(-1)N_x is nonzero. Write Q=x+ell(other parameters).
Exchange their roles, and substitute x=w+ell in Phi, taking Q as the new
parameter w. The new auxiliary matrix is invertible: after multiplication
by A^(-1), its determinant is the Q-coordinate of A^(-1)N_x, namely 1.
No equation terms change. Now auxiliary x appears only in E and is not
the output while a>d (that would be affine). Delete x and E; determinant
expansion on its single-occurrence column proves invertibility again.
The old parameter weight 6 is replaced by at most 7 and one row saves 2.

The same affine pivot can be used at zero potential cost when x has degree
>=2 and Q is a parent of an AND gate. Q has a nonzero column in A and is
used in the gate, so the new parameter has degree >=2. Both weights are 7,
and s,M,q,a remain fixed. The affine expression is substituted only into
Phi: no dense equation or uncharged gate is inserted into the computation.

A quadratic pivot is used only as follows. Let G=(x+b)(y+c)+e have only
one child D=G+z+h, where x,y,z are distinct unpaired parameters and z has
degree 1. Substitute z=w+G(x,y)+h in Phi and pair x,y. Then D=w; reroute
its uses to w and delete D and now-unused G. No other equation or gate
used z, so none needs an extra computation. The changes satisfy

    Delta s=-2, Delta M=0, Delta q=1, Delta a=-1,
    Delta H <= -4+3+(7-6)-weight(x)-weight(y).

For degrees of x,y at least 2 this is <=-14. Arbitrary fan-out of D
still gives new parameter weight at most 7.

## 3. Preliminary cases, in priority order

Work while a>d and restart after each reduction. Move eligible affine
gates before inspecting these cases. Each reduction has Delta a in
{0,-1,-2}, Delta H<=14*Delta a and strictly reduces M+s.

1. Perform any available D1-D4 deletion. Then Delta a=0 and Delta H<=-1.
   When none applies, gates have two distinct nonconstant parents and
   nondegenerate tables; ANDs on affine parents have independent linear
   parts; non-output gates have children; equations have three terms;
   and isolated unpaired parameters do not occur in equations.

2. Suppose an AND G on parameters has an unpaired isolated input z.
   Fix the other input w to force G constant and delete G. Parameter z
   becomes unused, though it may stay in Phi. If w is unpaired, its old
   weight is at least 6: Delta a=-1 and Delta H<=-6-6-2=-14.
   If w is paired, remove its pair and its partner becomes unpaired of
   weight at most 7: Delta a=0 and Delta H<=-3-6-2+7=-4.
   Rerouting a constant creates no occurrences.

   With cases 1-2 absent, no child can have both an isolated unpaired
   parameter and another parameter as parents. An equation is excluded
   by case 1; an XOR gate on parameters would already be in the system;
   an AND gate is case 2. This controls rerouting costs in later cases.

3. If a paired p is unused, fix its partner and erase the pair. The
   surviving p has weight zero: Delta a=0, Delta H<=-3.
   If a paired p has degree >=2, fix p, erase its pair and delete two of
   its distinct children. Gates have fixed pins; equations have at most
   two terms and are eliminated by D3. If one child disappears through
   earlier rerouting, it is already one of the two counted deletions.
   A term substitution in a short row never increases its term count.
   No originally isolated unpaired parameter co-occurs with p, by cases
   1-2. Only the now unpaired partner can gain weight, at most 7:
   Delta a=0, Delta H<=-3-4+7=0.
   When case 3 is absent all paired parameters have exactly one child.

4. If partners p,h share their one child B, impose p=h+c to make B
   constant and remove the pair. For a row v+p+h=b or an XOR any c
   works; for (p+alpha)(h+beta)+gamma choose c=alpha+beta+1.
   Delete B and reroute the constant. The surviving h is unused.
   Delta a=0, Delta H<=-3-2=-5.

5. Suppose paired p occurs in a row E_p, and its partner h has child B_h.
   These children are distinct by case 4. Fix h and delete B_h (D1 or D3).
   No previously isolated unpaired co-parent occurs there. Now p has at
   most its original row E_p as a child: B_h did not contain p, by p's
   uniqueness. If E_p survives, p has weight 6; pair removal and deletion
   have changed H by at most -3-2+6=1. Apply D4 to p,E_p for a further
   decrease of at least 1. Thus Delta a=0 and Delta H<=0. If E_p already
   vanished, count that deletion instead; p is unused and the estimate
   is better. Total s decreases at least 2.

6. Suppose an AND G on parameters has paired input p, partner h, and
   other input w. All paired parameters are isolated and absent from the
   affine system, so h's child B_h is a gate. Case 4 gives B_h!=G and
   w!=h; case 1 gives w!=p. Fix h and then fix w to force G constant.
   Delete G and B_h. The now unpaired p becomes unused.
   If w was unpaired, case 2 gives old weight 7, so Delta a=-1 and
   Delta H<=-3-7-4=-14. If w was paired, remove that distinct pair too;
   only its partner can gain weight, at most 7: Delta a=0 and
   Delta H<=-6-4+7=-3. No old isolated unpaired co-parent of B_h exists
   by cases 1-2, so no further increase has been omitted.

Restrictions reduce M and do not add gates or equations. Unused parameters
are retained in M until explicitly restricted; syntactic usage is not
substituted for affine dimension.

## 4. Prepare a minimal AND

Assume cases 1-6 are absent. The acyclic part is nonempty, or the output
would be affine. Its first gate G is AND-type since all eligible affine
gates were moved. Its parents are affine values with independent linear
parts. If a parent Q is auxiliary, it depends on some parameter x occurring
in the system. Cases 1,3,5 imply that x is unpaired and has degree >=2.
Pivot Q with x at zero potential cost, as above. If the second parent is
auxiliary, its linear independence supplies a parameter other than the
first new parent for a second pivot. Recheck preliminary cases after each
pivot, ending the iteration if one applies. At most two pivots are needed.
Otherwise G now reads two unpaired parameters x,y of degree >=2.

If a child J also reads x, replace h(x,G(x,y)) by its composed binary table
on (x,y), rerouting its inputs to x,y. This is exactly one permitted B2
gate. Only y gains occurrences and its weight already is 7. Do the same
with x,y exchanged. After each such operation, move newly eligible affine
gates and recheck cases 1-6. Each rerouting removes a child of G, and none
of these operations creates a new one; this preparation is finite.
If no preliminary reduction applies, no child of G reads x or y. G is not
the output by the terminal argument and has at least one child by case 1.
Write u=deg(x)>=v=deg(y)>=2 and h=deg(G)>=1.

## 5. Exhaustive remaining cases

7. If u=v=2 and the second child B of x,y is common, restrict (x,y) to
   a line making both G and B constant. An AND table is constant on
   three of the four inputs. A second Boolean table agrees on at least
   two of those three; the two points form an affine line. If B is a row,
   its auxiliary value is x+y+constant, so the same argument applies.
   The restriction fixes one parameter or identifies x with y+c.
   Delete G and B; the surviving parameter is unused. No weight grows
   under constant rerouting. Delta a=-1, Delta H<=-7-7-4=-18.

For later cases, if u=v=2 the other children B_x,B_y are distinct.

8. Let h_1 count G's children whose other parent is isolated and unpaired.
   Fix x to force G constant. Delete its u children and G's h children,
   a total of u+h distinct objects: no child of G reads x after preparation.
   Each row affected by x becomes short and each gate has a fixed pin.
   Any object deleted early through propagation is already a counted
   deletion; all remaining affected rows stay short. Thus the count is
   a lower bound even with coincidences and shared descendants.
   Removing x loses weight 7. A surviving old unpaired parameter gains
   weight only if it was isolated and its unique child is rerouted to it.
   Such a parameter cannot co-occur with x by cases 1-2. Only the h_1
   other parents of G's children can gain, each by at most 1. Hence

       Delta H <= -7-2(u+h)+h_1.

   If u>=3, this is <=-7-6-h<=-14. If u=2, also v=2. Since B_x!=B_y,
   and no child of G reads y, y gets no new occurrence and loses G:
   its weight drops by at least 1. Therefore in this case

       Delta H <= -12-2h+h_1.

   This is <=-14 for h>=2, or h=1,h_1=0. Delta a=-1.
   The sole uncovered configuration is u=v=2,h=h_1=1. Detect this
   before doing case 8: G's unique child D also reads an isolated
   unpaired z, distinct from x,y. Use case 9 or 10 instead.

9. If D is AND-type, write D=(G+c)(z+b)+e. Choose x0,y0 with G(x0,y0)=c
   and fix x=x0,y=y0. Such values exist because G is nondegenerate.
   Delete G,D,B_x,B_y, four distinct objects. G,D are rerouted to
   constants, and z becomes unused. B_x,B_y have no old isolated
   unpaired co-parent by cases 1-2. Thus no other weight grows and

       Delta a=-2, Delta H<=-7-7-6-8=-28.

   D is allowed to be the output. If its becoming constant leaves a>=d,
   this is an impossible state by the invariant; if a<d, termination is
   permitted. No extra successor of D is assumed.

10. Otherwise D is XOR-type. Apply the quadratic pivot of section 2.
    It deletes G,D, pairs x,y and replaces weight-6 z by a new parameter
    of weight at most 7. Thus

       Delta a=-1, Delta H<=-4+3+(7-6)-7-7=-14.

    This includes an output D, since its reference becomes the new
    parameter. All shared descendants are rerouted to that same value.

These cases are exhaustive because every nondegenerate B2 operation is
AND-type or XOR-type. The reduction to this local configuration was
proved, not assumed on the basis of finite pattern tests.

## 6. Termination and gate count

Each successful case preserves the representation, never increases M or s,
and strictly decreases their nonnegative integer sum. Preparations are
finite (topological moves, at most two pivots, finitely many reroutings).
There are at most n+|C| iterations. By exhaustive coverage while a>d,
the process stops with a_T<=d. Summing Delta H<=14*Delta a gives

    2|C|+7n >= H_0 >= H_T+14(n-a_T) >= 14(n-d),

because H_T>=0. Thus |C|>=(7/2)n-7d. If initially n<=d the right side
is nonpositive, and the conclusion is immediate.

Rows are introduced only one-for-one from charged affine gates. Affine
pivots change roles without adding computation; the quadratic pivot
eliminates an existing two-gate pattern. Therefore the potential proof
bounds the original gate count, including under arbitrary sharing.

## 7. What retuning these estimates can achieve

For a potential s+beta*q+alpha*I-sigma*ell, let I count used unpaired
parameters and ell isolated unpaired parameters; alpha>=sigma>=0,
beta>=0. If theta is the guaranteed decrease per lost unit of a, the
ordinary case estimates at (u,h,h_1)=(3,1,1) and (2,1,0) require

    theta<=alpha+4-sigma,       theta<=alpha+3+sigma.

Hence the derived coefficient theta-alpha is at most
min(4-sigma,3+sigma)<=7/2. Our H/2 has alpha=7/2,beta=3/2,sigma=1/2,
theta=7 and attains that value. This is only an optimization bound for
these local estimates, not a universal barrier to stronger circuit proofs.


## 8. An optional refinement for an AND child

This local refinement does not change the leading theorem. Assume the
representation and preliminary exclusions of Sections 1-4. Let G be an
AND-type gate on distinct unpaired parameters x,y of degrees u,v>=2,
with exactly one child D. No child of G reads x or y. Suppose D is
AND-type with other parent z, an unpaired degree-one parameter. Write

    D=(G+c)(z+b)+e.

Let j be the number of common children of x and y other than G; thus
0<=j<=min(u-1,v-1). These children may include affine rows and may share
descendants arbitrarily. Output references do not count as children.

There is always a legal reduction with

    Delta a=-2,     Delta H<=-20-2(u+v-j).                 (8.1)

Choose x0,y0 with G(x0,y0)=c, possible since G is nondegenerate, and
fix both parameters. This preserves the injective matching-quadratic map
on the remaining parameters: neither fixed parameter belongs to a pair.
G becomes c, D becomes e, and all their uses are rerouted to constants.
Delete all original children of x and y and D. The union of the two
child sets has u+v-1-j members, because G is their common child and
there are precisely j other common children. D is outside this union
by the preparation hypothesis. Hence u+v-j distinct charged objects
are deleted. A gate with a fixed input is deleted by D1, and a row
losing a parameter is short and deleted by D3. If an object has already
vanished by an earlier rerouting, count its deletion exactly once; no
claim requires its survival until a chosen position in the order.

x and y lose total weight 14. Parameter z becomes unused, losing weight
6. No initially isolated unpaired parameter occurs in any of the other
deleted objects: an affine row is excluded by D4, an affine gate on
parameters has already moved into the system, and an AND on parameters
with an isolated input is excluded by preliminary case 2. The only such
parameter at D is z, and D is rerouted to a constant. Thus no initially
isolated parameter gains weight. Initially unused parameters stay unused;
other unpaired parameters already have maximal weight 7. Pairs are
unchanged. These facts give (8.1), with no omitted rerouting charge.
Both M and a fall by two, q is unchanged, and s strictly decreases.

There is a stronger one-restriction option in a specified subcase. Let m
be the majority output value of G on its four assignments (the value
attained three times). If c=m, then there is a legal reduction with

    Delta a=-1,     Delta H<=-15-2*max(u,v).               (8.2)

For each input of an AND-type G there is a fixed value forcing its
majority output. Fix the higher-degree one, say x, to that value. Delete
its u original children and D, which are distinct. G and D are constants
and z becomes unused. The same no-isolated-coparent argument bounds all
other rerouting changes. Removing x loses 7, making z unused loses 6,
and the u+1 deletions lose 2(u+1), proving (8.2). The unfixed y may
keep its full weight 7; no decrease in its weight is needed.

These reductions allow D to be the output. If their resulting a is
at least d and their output is constant, the disperser invariant rules
out the starting configuration; otherwise they are ordinary permitted
reductions, including a final step crossing below d. They add no gate,
no auxiliary equation and no new computational operation. In particular,
(8.1) gives -32 at u=v=3,j=0 and -30 at u=v=3,j=1, compared with -28
for two units of dimension at the baseline rate. When j=2 it gives
only -28. Equation (8.2) gives at most -21 when max(u,v)=3.

Neither subcase is guaranteed to occur frequently in an arbitrary
circuit. XOR children and the other ordinary-case bottleneck remain.
Therefore these local savings alone imply no coefficient larger than
7/2 and no exponent larger than 1. They are supplemental partial results.


## 9. Removing an AND with a paired input without dimension loss

The following stronger alternative to preliminary case 6 does not change
the leading score. Assume preliminary cases 1-5 are absent, and let G be
an AND-type acyclic gate with paired parameter p as one parent and other
parent w. The other parent may be a parameter, an auxiliary value or a gate. If it is a
parameter it is distinct from p and its partner h; if it were an isolated
unpaired parameter, preliminary case 2 would apply. Parameter h has a
unique child B, distinct from G. All paired parameters have degree one
and occur in no affine row by preliminary cases 3 and 5.

While a>d, G cannot be the output: fixing paired p to the value forcing
G constant would preserve a and contradict the affine-space invariant.
It has at least one child, as dead non-output gates were removed. Let t
be its number of distinct children, and t_1 count the children whose
other parent is an initially isolated unpaired parameter.

Fix p to force G constant, remove p and its pair, and reroute/delete G
and its t original children. These t+1 objects are distinct. The children
are in the acyclic part; no gate feeds the affine system. No child of G
reads p, because p's only child was G. Every child becomes constant or
unary and is deleted by D1. If an earlier rerouting has already deleted
an original object, count its deletion once. The original affine system
is unchanged. This transformation preserves every represented value on
the remaining assignments, injectivity and the matching-quadratic map
form. M and q both fall by one, so Delta a=0.

No initially isolated unpaired parameter was the other parent w of G:
this is impossible when w is a gate, and otherwise follows from case 2.
Thus only the t_1 isolated unpaired other parents of the children of G
can gain weight. Each gains at most one, and all are distinct. Other
old unpaired weights cannot grow. Removing the pair loses three and
the t+1 charged deletions lose 2(t+1). Only newly unpaired h needs a
separate allowance.

If B is not a child of G, none of the deleted objects uses h and no
rerouting introduces h, so retain B and leave any further cleanup for the next iteration.
Then h remains isolated, with weight six. Consequently

    Delta H <= -3+6-2(t+1)+t_1 = 1-2t+t_1 <= 0,

since t_1<=t and t>=1. If B is a child of G, permit h its full possible
weight seven. B cannot contribute to t_1, because its other parent h
was paired before the change. Therefore t_1<=t-1 and

    Delta H <= -3+7-2(t+1)+t_1 = 2-2t+t_1 <= 1-t <= 0.

The argument covers shared descendants, B above or below other child
gates, arbitrary rerouting fan-out, and children that are the output.
It imposes no requirement that w be unpaired. It uses only one existing
permitted paired restriction and ordinary gate deletions, and strictly
reduces M+s. Hence an AND with a paired input can be removed at no loss
of a under these preliminary exclusions.

This is a local improvement over losing one unit of a in the earlier
unpaired-other-input branch. Its potential drop may be exactly zero,
so it does not supply a stronger leading coefficient on its own. The
ordinary unpaired XOR configurations remain outside this reduction.


## 10. A smaller finite dimension loss

**Theorem.** If f is nonconstant on every dimension-d affine subspace of
F_2^n, every full-B2 circuit for f has at least

    (7/2)n-(13/2)d                                         (10.1)

gates. This strengthens the finite loss in Section 6; the leading
coefficient remains 7/2. For n<=d the right side is nonpositive.
Assume n>d below. All maps, affine systems, semantic invariants, legal
pivots and terminal arguments are those of Sections 1-4. We change the
potential and three reductions, and recheck every case.

Use a=M-q as before and the nonnegative potential

    J=2s+2q+sum_(unpaired u) v(u),
    v(u)=0,5,6 for degree 0,1,at least 2 respectively.       (10.2)

Initially J<=2|C|+6n. We prove that while a>d a reduction strictly
lowers M+s and satisfies Delta J<=13*Delta a. Degrees still count
distinct child objects, not pins or output references. Moving an
eligible affine gate into the system leaves s unchanged and cannot
increase J. The two affine exchanges used to prepare a minimal AND
are between weight-six unpaired parameters and auxiliaries, and cost
zero. The composition reroutings of Section 4 add occurrences only to
parameters already of weight six, so also do not increase J.

### 10.1 Preliminary deletions and restrictions

D1, D3 and D4 each save two by removing one charged object and permit
at most one unit of unpaired-weight increase (five to six). Thus each
lowers J by at least one without changing a. D2 saves two. Their
matrix/injectivity justifications are unchanged; no new operation is
assumed.

For the isolated-input AND case 2, let G have an unpaired isolated
input z, and other parameter w. If w is unpaired, fix w to force G
constant. G cannot be the output while a>d: that restriction would
make the output constant and leave a-1>=d. Hence G has a child D.
Delete G and D, leaving any other degenerate objects for later cleanup.
The two are distinct; D cannot read z because z's only child was G.
Removing w loses at least five, making z unused loses five, and the
two deletions save four. At most the initially isolated unpaired other
parent of D gains one. Consequently

    Delta a=-1,  Delta J<=-5-5-4+1=-13.                    (10.3)

If D also reads fixed w it becomes constant and incurs no gain. If
w is paired, the original one-gate reduction already suffices: remove
its pair, make z unused and delete G. The newly unpaired partner gets
at most six, so Delta a=0 and Delta J<=-2-5-2+6=-3. These estimates
also cover both inputs initially isolated. Constant rerouting at G
creates no other occurrences.

An unused paired parameter is handled by fixing its partner and erasing
the pair; the unused survivor has weight zero and Delta J<=-2. For a
paired parameter of degree at least two, fix it, erase its pair and
delete two distinct original child objects. The partner gains at most
six, against pair removal two and deletions four: Delta J<=0 and
Delta a=0. As in Section 3, the absence of basic/isolated-input cases
excludes old isolated unpaired coparents of these objects. Short-row
deletions and earlier removal of a counted object are treated exactly
as before.

When partners share their unique child, impose the same affine line
that makes that child constant. The remaining parameter is unused;
pair removal and deletion give Delta J<=-2-2=-4 at Delta a=0.

For paired p in a row, fix its partner and delete that partner's distinct
child. Pair removal and this deletion cost -2-2, and p becomes unpaired
of weight at most five, so the increase is at most one. If its row
survives, apply D4 for a further decrease of at least one. If the row
already vanished, its deletion gives an even better estimate. Thus
case 5 again has Delta a=0 and Delta J<=0.

Use the dimension-preserving paired-AND operation of Section 9 instead
of the old dimension-losing case 6. Its proof uses only legal paired
restriction and specified gate deletions. If G has t>=1 children and
t_1 initially isolated unpaired other parents, and B is the partner's
unique child, the reweighted estimates are

    -2+5-2(t+1)+t_1 = 1-2t+t_1 <= 0       if B is outside,
    -2+6-2(t+1)+t_1 = 2-2t+t_1 <= 0       if B is inside.

In the second line t_1<=t-1. Retain B when it is outside the prescribed
deletion set. An output G is impossible because fixing its paired
input preserves a and makes the output constant. This handles all of
case 6 at Delta a=0, including auxiliary/gate other parents.

### 10.2 Prepare and eliminate a minimal AND

The preparation in Section 4 remains valid with weight six replacing
seven. If a preliminary case becomes applicable, perform its reduction
and restart. Otherwise the minimal AND G has distinct unpaired parents
x,y of degrees u>=v>=2, no child of G reads x or y, G is not the output,
and h=deg(G)>=1.

Case 7, with u=v=2 and a shared second child, has Delta a=-1 and
Delta J<=-6-6-4=-16<=-13 by the same affine-line restriction.
For the remaining cases with u=v=2 their other children B_x,B_y are
distinct.

For case 8 let h_1 count children of G whose other input is an initially
isolated unpaired parameter. Fix x to force G constant and delete its
u original children and G's h original children. The count u+h is
valid after preparation. The same no-isolated-coparent argument gives

    Delta J<=-6-2(u+h)+h_1.

For u>=3 this is at most -6-6-h<=-13. If u=2 then v=2 and y loses
G without a new occurrence, so its weight drops by at least one:

    Delta J<=-11-2h+h_1.

This is <=-13 if h>=2, or if h=1,h_1=0. Every reduction in this
paragraph loses exactly one unit of a. The sole remaining configuration
has u=v=2,h=h_1=1: G's unique child D reads an isolated unpaired
parameter z, distinct from x,y.

If D is XOR-type, use the existing quadratic pivot. Deleting G,D costs
four, adding the pair costs two, x,y lose six each, and the new parameter
replaces weight-five z by weight at most six. Thus

    Delta a=-1,  Delta J<=-4+2-6-6+(6-5)=-13.              (10.4)

### 10.3 The AND child with an isolated parameter

Suppose instead D is AND-type, with parents G,z as just described.
The old two-fixing reduction is replaced by the following exhaustive
choices. D cannot be the output while a>d: fixing unpaired z to a
forcing value would preserve an affine image of dimension a-1>=d and
make the output constant. Let t=deg(D)>=1.

**At least two children.** If t>=2, fix z to force D constant. Delete
D and all its t original children; G is then dead and is deleted too.
These t+2 objects are distinct. Let eta count how many of x,y have their
other child among D's children, and let t_1 count D's children with an
initially isolated unpaired other parent. The two categories are disjoint,
so eta+t_1<=t. If an x/y second child is outside this set it is retained;
that parameter loses G and drops from six to five. If its second child
is in the set allow its full possible weight six after rerouting. Other
old isolated unpaired parents gain at most one each; old high-degree
unpaired parameters cannot gain weight. Thus

    Delta a=-1,
    Delta J<=-5-2(t+2)-(2-eta)+t_1
             =-11-2t+eta+t_1<=-11-t<=-13.                 (10.5)

No new charge is hidden in optional propagation: delete the specified
objects only and retain other degenerate objects for the next iteration.

**One child without x/y.** Suppose t=1 and let E be D's child. If E
reads neither x nor y, it is distinct from B_x,B_y. Choose constants
x0,y0 making G a value that forces D constant, and fix both. Delete
G,D,B_x,B_y,E: five distinct original objects. Their ordering need not
keep an object that has already disappeared; count every deletion once.
x,y lose twelve, z becomes unused and loses five, and five deletions
save ten. No old isolated unpaired coparent occurs at B_x or B_y by
preliminary exclusions. At E at most its initially isolated unpaired
other parent can gain one. Hence

    Delta a=-2,  Delta J<=-6-6-5-10+1=-26.                (10.6)

This includes E as output. The condition that E reads neither x nor y
is what makes the fifth deletion distinct.

**One child that reads x/y and is AND-type.** Suppose E also reads x;
then E=B_x, while B_y is distinct. If E is AND-type, the chain G,D,E
can be replaced by two charged B2 gates without a restriction. There
is a value x=c at which E(D,x) is constant. At the other value x=1-c,
G(x,y) is a unary function of y, so

    H(y,z)=E(D(G(1-c,y),z),1-c)

is one B2 truth table. Compute H by one gate and then the function
that is the constant E value when x=c and H when x=1-c by one gate.
This computes the original E at all eight assignments. G and D each
have only the next chain gate as child; keep E's output identity and
replace the chain by these two gates. x now has one child instead of
two, y and z gain no new child, and no other parameter gains an
occurrence. Thus Delta a=0, Delta s=-1 and Delta J<=-2-1=-3.
Even if H is degenerate, retain it as one charged B2 gate until later
normalization. There is no free negation, uncharged gate or change of
input map. The construction with y in place of x is identical.

**One child that reads x/y and is XOR-type.** If E reads x and is
XOR-type, fix y to force G constant. Delete G and the other child B_y
(D1 or D3). D becomes constant or a literal of z. If constant, E is a
literal of x and can also be deleted. If D is a literal of z, E is
x+z+c. At this point x and z have no other computational occurrence:
x's original children were G,E, and z's only original child was D.
Use the affine reparameterization w=x+z+c, replacing z by w+x+c in
Phi. Reroute E to w and delete E, with D deleted as well. This is a
bijection on the remaining parameters; x,z were unpaired, so Phi
retains the matching-quadratic form. In either branch G,D,E,B_y are
four distinct deleted objects. y loses six. Originally x,z had total
weight eleven; afterwards the one survivor used in place of E has
weight at most six and the other is unused, or z is unused and x is
the literal survivor. This loses at least five more. No initially
isolated unpaired coparent of B_y is available under the preliminary
exclusions. Consequently

    Delta a=-1,  Delta J<=-6-5-8=-19<=-13.                 (10.7)

Swap x and y for the other orientation. AND-type and XOR-type exhaust
all nondegenerate possibilities for E; degenerate gates were removed
beforehand. Together (10.5)-(10.7) exhaust the AND-child configuration.

### 10.4 Termination and the finite theorem

Every transformation uses the earlier legal substitutions, pivots and
charged deletions, or the explicitly supplied two-gate B2 composition.
It preserves the represented function on every remaining assignment,
the invertible affine system, injectivity and the matching-quadratic
map. No reduction increases M or s, and every successful reduction
strictly decreases M+s. Preparations remain finite. Therefore the
exhaustive process stops at a_T<=d. Summing Delta J<=13*Delta a and
using J_T>=0 gives

    2|C|+6n >= J_0 >= J_T+13(n-a_T) >=13(n-d).

Rearranging proves (10.1). The original full-B2 normalization and
terminal arguments cover constants, unary/negated projections, repeated
pins, arbitrary sharing and outputs. This proof changes the finite loss
coefficient, not the leading 7/2 coefficient.

## 11. Affine maps with promised pairs: the 4n-10d theorem

This is a separate strengthened argument. It replaces, rather than silently
reuses, the matching-quadratic invariant of Sections 1--10. It supplies the
argument of Mishra--Santhanam, ECCC TR26-260 revision 1, November 3, 2026,
Sections 2--4, specialized to the pinned full B2 convention. The old
supplemental lemmas remain valid in their own stated representation.

**Theorem.** For every total Boolean f on n bits nonconstant on every
dimension-d affine subspace, d>=1, any full-B2 circuit computing f has
at least 4n-10d gates. No uniformity hypothesis is needed for this finite
gate statement. If d>=n its right side is nonpositive, so assume n>d.

### 11.1 New state and restricted correctness domain

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

### 11.2 Terminal outputs and computation-preserving simplifications

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
Their validity uses the same elementary row and determinant arguments
as Section 2, summarized here to make this new invariant explicit.

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
4. If unpaired x occurs only in one row E, Nx is nonzero. Choose an
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

### 11.3 Paired cleanup, before choosing a minimal AND

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
counted deletion; surviving p is unused. Otherwise delete a short Hp,
or use deletion 4 when it has three terms. This removes a second object
and leaves at most one new unpaired coordinate of weight 6. Other old
unpaired coordinates cannot become newly used. Thus Delta W<=-4-2+6=0.
After (C) fails, all paired coordinates are isolated, occur only in the
acyclic part, and have a different child from their partner.

### 11.4 Exhaustive minimal-AND step

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

### 11.5 Termination and conclusion

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
