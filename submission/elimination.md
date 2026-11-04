# Ternary-pair gate elimination for an affine disperser

**Theorem.** Let d be a nonnegative integer. Suppose f:F_2^n -> F_2 is
nonconstant on every affine subspace of dimension strictly greater than d.
Every circuit for f in the full binary basis has at least

    4n - 9(d+1)

gates. This includes arbitrary sharing, repeated pins, constants, unary
tables, degenerate gates and a single terminal output.

The fair affine system, role exchanges and deletion machinery below are
adapted from Mishra--Santhanam, ECCC TR26-260, 2026-10-19, Sections 2-4.
The new domain of disjoint three-point pairs comes from the team's
astra-quadratic lane; the weights and the affine-disperser variant are
the team's ter-table table. This is a new proof of that table, rather than
an invocation of the earlier paper's circuit theorem.

## 1. State, semantics and potential

Start with any acyclic B2 circuit C for f. Count every retained gate,
including constant and unary gates, until its deletion. We allow a free
final constant or literal and negations on edges: an edge negation is
absorbed in the truth table of its child. Proving a lower bound in this
slightly more permissive output convention proves it in the pinned model.
The ten binary tables depending on both arguments are precisely the two
tables u+v+c and the eight tables (u+b)(v+c)+e over F_2.

A state consists of the following objects.

* There are m parameters u and an affine map Phi:F_2^m -> F_2^n that is
  injective on the entire parameter space, not just on the allowed domain.
* A set P of q disjoint unordered pairs of parameters is marked. A pair
  (p,r) has an excluded corner (e_p,e_r) in F_2^2. Its allowed set is the
  other three points. Dom is the Cartesian product of these sets and the
  unrestricted unpaired bits. Put a=m-q.
* The computation is a fair affine system A v = N u + c followed by an
  acyclic list of B2 gates. A is square and invertible. Each equation has
  at most three distinct nonconstant terms, after cancellation over F_2.
  Gates read parameters, auxiliary values v, or earlier gates; no gate
  feeds an equation. The output is a constant or a literal of an internal
  value. On Dom the output equals f(Phi(u)). Outside Dom no agreement is
  asserted. The affine system itself is uniquely solvable everywhere.

Let s be the number of equations plus acyclic gates. Degree counts
distinct child equations and gates, not pins or the output reference.
An equation is a child of every value appearing in it. Let I count the
unpaired parameters with positive degree, and set

    W = s + 5I + 3q.

Initially Phi is the identity, Dom=F_2^n, m=a=n, q=0, s=|C|, and
W_0<=|C|+5n. Always W>=0. Unused parameters are retained in Phi and m;
usage is never substituted for dimension.

For each pair, fixing p to its **safe value** 1-e_p leaves both values of
r allowed. Fix one member of each pair safely. The resulting affine
parameter subspace has dimension m-q=a, lies in Dom, and has an injective
affine image under Phi. Thus f(Phi(u)) is nonconstant on Dom whenever
a>d. This argument uses only affine dispersal.

Here and below Delta means new minus old. We run the algorithm only at

    a >= d+2.

At such a state the output cannot be affine in the parameters. On an
a-dimensional affine subspace just constructed, an affine Boolean
function is constant either everywhere or on a hyperplane of dimension
a-1>d. Its image would contradict the hypothesis on f. The same argument
excludes an output literal of an AND-type gate with affine parents: on
that subspace its factors are affine; a nonconstant factor can be set to
zero on a hyperplane, and if both factors are constant the output already
is constant. This also excludes parameter and auxiliary output literals,
including a parameter with no counted child. We will use the constant
output contradiction after restrictions even when the new a is d+1.

## 2. Legal operations and their implementation

Restrictions are made by composing Phi with an injective affine map tau
from the new full parameter space to the old one. The new Phi remains
affine and injective. Every restriction below maps its new product domain
inside the old Dom. No nonlinear reparameterisation of Phi is used.

Fixing an unpaired bit to a constant removes one parameter, preserves q,
and decreases a by one. Fixing one member of a pair safely removes that
parameter and the pair, leaves its partner unrestricted, and preserves a.
Restricting a paired three-point set to a line through two of its allowed
points removes one parameter and its pair, also preserving a. Such a line
is implemented by fixing one coordinate or identifying the coordinates
up to a constant. Each implementation substitutes only constants or
literals into the computational representation, so creates no gates or
longer equations.

We also exchange an unpaired parameter x with an auxiliary Q having solved
expression Q=x+ell(u_without_x), with ell affine. Substitute
x=w+ell in Phi, rename Q as parameter w and x as an auxiliary, and leave
the list of equations and gates otherwise unchanged. This is an invertible
affine coordinate change. Because x was an unrestricted bit, for each
fixed assignment of the other parameters it is a bijection of that bit.
Consequently its new domain is exactly the same product of marked
three-point sets and unrestricted bits, with w replacing x; no constraint
on a pair is changed, even if ell mentions paired bits.

For completeness, if Q is column j of A and N_x is the column of x,
then the new auxiliary matrix replaces column j of A by N_x. Multiplying
by A^(-1) gives an identity matrix whose column j is A^(-1)N_x. Its
determinant is that vector's j-coordinate, namely 1. The matrix stays
invertible; equation sparsity is unchanged. We do not substitute the
dense ell in the computation. The old x was used, and Q has a nonzero
column in A, so the new unpaired w is used. Other parameters have the
same children. Thus an exchange alone has Delta I=Delta s=Delta q=0.

Rerouting replaces every use, including the output, of a value that has
become a constant or a literal of another value by that constant or
literal. It never makes a previously unused parameter used: any surviving
parameter receiving an occurrence was a parent of the deleted object and
therefore was already used. This observation concerns I, which has the
same weight at every positive degree. Increasing a used parameter's
fan-out has no potential cost.

Rerouting a gate to an input creates no cycle: a gate input precedes it
topologically, and subsequent substitutions still use earlier surviving
values. Parameter and auxiliary values precede the acyclic part. Cycles
inside the affine system are permitted, with unique solvability proved
by invertibility. No nonlinear gate ever feeds that system.

## 3. Deletions and affine preparation

Each deletion below has Delta a=Delta q=0, Delta s=-1, and Delta I<=0,
so Delta W<=-1.

**D1.** Delete a constant or unary gate, a gate with a fixed or repeated
pin, or an AND-type gate whose affine parents have dependent linear
parts. Two dependent affine Boolean forms are constant or equal up to a
constant, so in the last case the output is a constant or a literal of
one of the existing parents. Reroute to it. This covers all sixteen truth
tables without retaining an uncharged extra operation.

**D2.** Delete a gate with no child unless it is the output. Its parents
can only lose occurrences. The terminal exclusions above prevent a
minimal AND or an affine value from being an exceptional output later.

**D3.** An equation with at most two terms has an auxiliary term, since
a zero row of A is impossible. Write it v=U+c, where U is zero or one
other value. Substitute this literal or constant in every occurrence of
v and delete v and the equation. Adding that row to the other rows using
v preserves determinant; afterwards v's column has one entry. Expanding
along that column proves the remaining matrix invertible. Term counts
cannot increase, and no unused parameter becomes used.

**D4.** If an unpaired isolated x has unique child equation E, its column
N_x is nonzero. Invertibility gives A^(-1)N_x!=0, so choose a dependent
auxiliary Q and exchange x with Q as in Section 2. The now auxiliary x
occurs only in E. It is not the output, by the affine-output exclusion.
Delete x and E without substitution. Expansion on its one-entry column
proves the remaining matrix invertible, and determines the removed x
uniquely from the surviving values. The old used unpaired x is replaced
by at most one used unpaired parameter Q, so Delta I<=0. This includes
the possibility that Q becomes unused after deleting E.

**M.** Move every affine gate whose parents are parameters or auxiliaries
into the affine system: replace J=b_1 t_1+b_2 t_2+c by its equation.
The auxiliary matrix becomes block triangular with an additional diagonal
1. One charged gate becomes one charged equation. Parameter usage does
not increase, so Delta W<=0 with m,s,q,a unchanged. Apply this in
topological order, then apply available deletions. It is finite because
each move removes a gate from the acyclic list.

The deletion arguments also justify the following batch counting rule.
After fixing a parameter, every original child gate has a constant pin
and every original child equation has at most two terms. Each can be
deleted. Literal substitutions never lengthen an equation or restore a
nonconstant fixed pin. If an object is removed during another scheduled
deletion, its deletion has already been counted; it is not counted twice.
Hence any specified set of distinct such objects supplies that many
deletions even with common parents, common descendants or cyclic affine
dependencies. Unscheduled simplifications can only help the estimates.

## 4. Preliminary paired rows

Exhaust D1-D4 and M before the following rows. Restart after a paired
row. All operations in these rows preserve a and do not increase W.

**P0: an unused paired p.** Fix its partner r safely and remove the pair.
The newly unpaired p remains unused. Thus Delta m=Delta q=-1,
Delta I<=0, Delta s<=0, and Delta W<=-3. Any simplification of r's
children can be left to the next iteration.

**P2: a paired p with at least two children.** Fix p safely and delete
at least two distinct children by D1 or D3. The batch counting rule
applies to equations as well as gates. Its partner r is the only newly
unpaired parameter and contributes at most one to I. Consequently

    Delta m=-1, Delta q=-1, Delta s<=-2, Delta I<=1,
    Delta a=0, Delta W<=-2-3+5=0.

With P0 and P2 absent, every paired parameter has exactly one child.

**PC: partners p,r share their child H.** If H is a gate, its two pins
are p,r. Among the three allowed assignments, two give the same Boolean
value at H; restrict the pair to the affine line consisting of those
two points. Then H is constant. If H is an equation, after D3 closure
it has exactly three terms and reads v+p+r=c for an auxiliary v. Apply
the same pigeonhole argument to p+r; v becomes constant. Delete H by
D1 or D3 and reroute the constant. The one surviving coordinate of the
line is unused because p,r had no other child. Thus

    Delta m=-1, Delta q=-1, Delta s<=-1, Delta I<=0,
    Delta a=0, Delta W<=-1-3=-4.

This choice of a line is inside the allowed three-point set; an arbitrary
old-MS diagonal restriction is not being assumed legal.

**PE: an isolated paired p occurs in an equation E.** Its partner r has
a unique different child H, by PC. Fix r safely and delete H. H may be
a gate or an equation. H does not contain p, since p's only child was E,
so its literal rerouting cannot introduce p elsewhere. The newly unpaired
p still occurs only in E, unless E has already been removed, in which
case that is a second counted deletion. Otherwise apply D4 to p,E.
This removes E and replaces p by at most one used unpaired auxiliary.
The two distinct objects H,E are deleted, and only this new unpaired
parameter can add to I. Therefore

    Delta m=-1, Delta q=-1, Delta s<=-2, Delta I<=1,
    Delta a=0, Delta W<=-2-3+5=0.

In particular, after these rows every paired parameter is isolated and
absent from the affine system. Every parameter appearing in an equation
is unpaired and used. All these statements are syntactic and do not
assume that auxiliary dependencies are acyclic.

## 5. Produce a minimal AND with parameter parents

Assume no earlier row applies and a>=d+2. If there were no acyclic gate,
the output would be affine, which is impossible. Choose the first gate
G in topological order. By M and D1 it is AND-type on two affine values
with independent linear parts. It is not the output by Section 1, and
D2 gives it a child.

**AR.** If a parent Q is auxiliary, choose an unpaired parameter x with
coefficient 1 in Q and exchange them. Such x exists: Q is nonconstant,
and every parameter occurring in the system is unpaired, by PE. A
parameter on which Q depends has a nonzero N column and so is used.
The exchange preserves the product domain, fairness and W by Section 2.
If G already has one parameter parent, choose x different from it. The
linear independence of the two parent forms guarantees such a coefficient
in Q. An invertible affine change preserves that independence. Thus at
most two exchanges make both parents parameters.

After each exchange recheck D1-D4, M and the paired rows; if a reducing
row applies, do it and end the iteration. Otherwise continue with G.
No gate operations or nonlinear ancestor relations have changed. In the
remaining state G has distinct parameter parents and at least one child,
and is not the output. There is no need to restructure children also
reading those parameters: I charges every used unpaired bit equally.

## 6. Eliminate this AND: all parameter cases

Write G=(v+b)(w+c)+e. Each paired parent, if present, is isolated. A
paired parent cannot have its partner as the other parent of G, because
that would be PC. The two pairs below, when present, are distinct.

**AP1: one paired parent p and one unpaired parent w.** Let r be p's
partner, with unique child H. H is a gate, by PE, and H!=G, by PC.
Fix r safely, then fix w to c to force G constant. These two restrictions
remove one pair and two parameters, so Delta a=-1. Delete H and G.
The newly unpaired p becomes unused; the old used unpaired w is removed.
No other previously unused parameter becomes used. Thus Delta I<=-1.

There is a third distinct deletable gate, including every overlap:
choose a child K of G. If a child K!=H exists, G's constant value makes
K deletable. It may also read w or some other scheduled value; it still
has a fixed pin. If all children of G equal H, H reads G and r. Once
both are fixed, H is constant. If H is not output, D2 closure in the
starting state gives a child J of H. Acyclicity makes J distinct from
G,H, and H's constant value makes J deletable. If H is output, after
these restrictions the output is constant on a domain with a'=a-1>d,
contradicting Section 1. Thus this last possibility cannot occur at a
loop state. Early deletion of K or J already supplies the same third
distinct deletion. H's other parent may be G, a descendant of G, or an
unrelated value; the preceding split covers all possibilities.

Consequently

    Delta m=-2, Delta q=-1, Delta s<=-3, Delta I<=-1,
    Delta a=-1, Delta W<=-3-3-5=-11 <= -9.

**AP2: both parents p,w are paired.** Let r,h be their respective
partners, with unique child gates H_r,H_h. Neither child is G, since
the only parents of G are p,w. Fix r and h safely, making p,w unpaired,
then fix w to c to kill G. Both pairs and three parameters disappear,
so Delta a=-1. The surviving p is unused and w is removed. The old
paired coordinates had no I charge, and no other unpaired usage can
appear, so Delta I<=0.

If H_r!=H_h, the three distinct gates G,H_r,H_h are deletable. If
H_r=H_h=H, that gate reads r,h. It cannot read G as a third pin. It
becomes constant after the two safe fixings. A child K of G exists and
is distinct from G,H. The fixed G pin makes K deletable. This gives
three distinct deletions also in the shared-child case. A child of G
may equal one of H_r,H_h in the first case; those three gates were
already distinct and suffice. Literal rerouting, shared descendants,
and an output among these children do not change the count. Hence

    Delta m=-3, Delta q=-2, Delta s<=-3, Delta I<=0,
    Delta a=-1, Delta W<=-3-6=-9.

**GP: both parents x,y are unpaired.** A nondegenerate AND-type truth
table has one exceptional corner and the same value at its other three
corners. Add the new pair (x,y), excluding exactly that exceptional
corner, and restrict Dom to its three allowed points. Phi does not
change at all. The pairs remain disjoint. This is a domain restriction,
not a bijective quadratic change of coordinates.

On the new Dom, G is constant. Delete G and all its k>=1 distinct child
gates. Each child now has a constant pin and is a unary table or a
constant, even if its other parent is x or y, another child of G, or a
shared descendant. Process them topologically; if one was already
deleted by propagation it is counted once. No equation reads G. The
parents x,y were each used because they read into G, and now lose their
unpaired I charges. All other unpaired parameters were already used
before any literal rerouting to them. Thus

    Delta m=0, Delta q=1, Delta s<=-(1+k), Delta I<=-2,
    Delta a=-1, Delta W<=-(1+k)-10+3=-8-k <= -9.

These three cases cover all parameter-parent types and all eight
AND-type operations; their polarities only select safe or killing
constants and the exceptional corner. Every fixing is legal: AP1 fixes
an old unpaired w arbitrarily, AP2 fixes w only after its partner has
been fixed safely, and GP pairs only previously unpaired bits.

## 7. Complete accounting, coverage and termination

For reference the upper bounds on changes used above are tabulated here.
The Delta I column bounds the total change of I, including every
rerouting and newly unpaired parameter. Negative bounds need no equality.

| Row | Delta s | Delta I | Delta q | Delta a | Delta W upper bound |
|---|---:|---:|---:|---:|---:|
| D1 | -1 | 0 | 0 | 0 | -1 |
| D2 | -1 | 0 | 0 | 0 | -1 |
| D3 | -1 | 0 | 0 | 0 | -1 |
| D4 | -1 | 0 | 0 | 0 | -1 |
| M | 0 | 0 | 0 | 0 | 0 |
| AR | 0 | 0 | 0 | 0 | 0 |
| P0 | 0 | 0 | -1 | 0 | -3 |
| P2 | -2 | 1 | -1 | 0 | 0 |
| PC | -1 | 0 | -1 | 0 | -4 |
| PE | -2 | 1 | -1 | 0 | 0 |
| AP1 | -3 | -1 | -1 | -1 | -11 |
| AP2 | -3 | 0 | -2 | -1 | -9 |
| GP, k>=1 | -(1+k) | -2 | 1 | -1 | -8-k |

For general nonnegative weights W=g*s+i*I+p*q, the corresponding
bounds are -g for each D row, 0 for M and AR, -p for P0,
-2g+i-p for P2 and PE, -g-p for PC, -3g-i-p for AP1,
-3g-2p for AP2, and -(1+k)g-2i+p for GP. At (g,i,p)=(1,5,3),
each reducing row satisfies Delta W<=9 Delta a and Delta a is 0 or -1.

While a>=d+2, do finite M preparation and the first available deletion
or paired row. If none applies, Section 5 supplies a minimal AND and at
most two affine exchanges, with the stated rechecks. Then exactly one
of AP1, AP2, GP applies. This proves universal coverage; it is not an
inference from finitely many tested local graphs.

Each reducing row preserves the state invariants, never increases m or
s, and strictly decreases the nonnegative integer m+s. In P0 even when
s is unchanged, m falls by one. GP decreases s, although m is unchanged.
Every other row either deletes an object or restricts parameters. The
preparations are finite: M removes acyclic gates and at most two AR
exchanges occur before a reducing row. There are therefore at most
n+|C| reducing iterations. Since coverage holds at every loop state,
the algorithm terminates with a_T<=d+1. If initially n<=d+1 the claimed
lower bound is already nonpositive.

Sum Delta W<=9 Delta a, including the nonincreasing preparations:

    |C|+5n >= W_0 >= W_T+9(n-a_T) >= 9(n-d-1).

Rearranging proves |C|>=4n-9(d+1). Each affine row was introduced
one-for-one from a charged gate, exchanges changed roles without extra
computation, and restrictions never added a gate. Thus this inequality
bounds the original gate count in the full B2 model, not a potential
mistaken for a gate count and not a formula or restricted-basis size.
