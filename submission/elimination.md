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
