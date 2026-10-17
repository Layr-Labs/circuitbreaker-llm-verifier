# E. Gate elimination at the integer weights

This is a proof of the local assertion in C7, rather than an invocation of the
published circuit theorem. Its organization follows Li-Yang, ECCC TR21-023
(2021-02-20), proof of Theorem 1.1, Cases 0-8. The numerical weights are fixed
throughout, so no optimization or LP theorem is used.

## E0. Notation, accounting and choice of cases

A u-variable or u-gate has exactly u uses; u+ and u- mean at least and at most u.
Uses are actual pins, including C1's terminal convention. A descendant mentioned
in a count is an immediate consumer, unless a path is specified. A gate fed by
one constant is removed by constant/unary normalization. "Kill" an AND-type gate
means set a pin to its killing value; its table ((a+c)(b+d))+e shows explicitly
that this is always possible, independently of the other pin. XOR gates have
table a+b+c, so an affine identification a=b+c' fixes either chosen output value.

Every state below is live, with s>2d+2 at the beginning of a block. Up to three
restrictions keep s>=2d. Consequently the output is never constant, and the
topmost AND gates used in counts have real gate descendants (C2).
Apply C3-C5 first: remove dead, constant, unary, repeated-pin and useless gates.
If a packing is nonempty use C6. Thus before the following cases the packing is
empty, every gate genuinely depends on its two pins, and no gate is useless.
A gate G with one consumer H sharing another pin with G is useless in this
sense, since composition removes G without increasing the number of pins of H.

All additional normalizations have nonnegative measure decrement. Here are
frequently used sufficient accounts, in the order (gates removed, influential
bits lost, quadratic equations lost, upper bound on potential increase):

    (4,1,0,4): 20+48-4=64;
    (3,1,0,-1): 15+48+1=64;
    (2,1,1,2): 10+48+9-2=65;
    (4,1,0,2): 20+48-2=66;
    (5,1,-1,0): 25+48-9=64;
    (6,2,1,6): 30+96+9-6=129>=128;
    (4,2,2,3): 20+96+18-3=131>=128.

Losing at least two influential bits with one restriction, three with two, or
four with three also suffices: 96>=64, 144>=128, and 192=192. When such a count
is used, gate normalizations contribute nonnegative decrement; removal of
quadratic equations only helps. Distinctness of the counted influential bits
is checked below through protection, use counts, and their specified pins.

For gate counts under constant propagation use C3's constant-frontier argument.
It counts pins, not distinct named consumers. A coincident consumer with two
constant pins becomes constant and passes a constant to a further consumer,
giving the same total lower bound. Potential bounds are charged per removed
gate's affected other predecessor. A dead gate costs at most two potential
units; a constant/unary gate at most one. Sharper zero bounds below specify
why its predecessors cannot cause new troubled structure.

Cases are tried in their numbered order. Later cases assume all earlier ones
unavailable in the starting state. A case describes a whole block of at most
three restrictions; it need not rerun the ordering between restrictions in that
block. Each multi-restriction account below includes its preceding steps.

## E1. Protected bits with large use count or a nonlinear use

If a protected free bit x feeds at least two gates, fix it. Its two constant-pin
uses eliminate at least two gates. If x feeds an AND-type gate, choose its killing
value; this removes that gate and a consumer, also at least two gates. There is
no constant output. Each deletion has potential increase at most one. One
influential bit and its quadratic equation are lost, giving (2,1,1,2), hence 65.

## E2. A protected bit without a circuit use

If protected x has no use, fix its mate y. The equation for the pair becomes
affine; x is now unprotected with no use. Both x and y cease to be influential.
This gives at least 96 in one restriction. Thereafter every protected bit has
exactly one use, and that use feeds an XOR gate.

## E3. Large total degree at a primary-input AND

Suppose x feeds an AND gate G and deg(x)+deg(G)>=4. Fix x to kill G.
Besides G, all other uses of x and all uses of G receive constants. There are
at least deg(x)-1+deg(G)>=3 such pin edges. C3's frontier count gives at least
three further deletions, even if a consumer of x is also a consumer of G.
Four constant/unary deletions have potential increase at most four. One bit is
lost, giving (4,1,0,4). After this case, any free bit feeding an AND has degree
one or two.

## E4. A one-use bit in a two-primary-input AND

Let G have free inputs x,y with deg(x)=1. They are unprotected by E1. Fix y to
kill G. Both fixed y and now unused x cease to be influential, giving 96.

## E5. A two-primary-input AND

The remaining such gate G has free inputs x,y both of degree two and, by E3,
degree one itself. Both bits are unprotected. Thus G is troubled. Let B be x's
other consumer, C be y's other consumer, and D be G's consumer. Useless-gate
normalization implies B!=D and C!=D. Fixing x to kill G deletes G,B,D;
fixing y instead deletes G,C,D. If either route introduces no new troubled
gate while deleting these three, the old G disappears from the empty packing,
so dPhi<=-1. The account (3,1,0,-1) gives 64.

It remains to show that whenever both routes would introduce trouble, a suitable
alternative works. The descriptions below give the needed graph exhaustive
split, the alternative restrictions, and their accounts.

### E5.1. B=C

Both gates G and B depend on the same two bits. C6 supplies an affine line
fixing both outputs. Both bits are unprotected and have no other uses, so one
restriction loses two influential bits.

### E5.2. D feeds both B and C

Try fixing x and delete G,D,B in that order. Only D's other input can produce
new trouble: the other affected bit y becomes one-use, and a gate predecessor
of D retains its uses under bypassing. Therefore that input must be a free bit
u passed through D and then B. For u to end with exactly two uses, B must
originally have one use, to an AND gate E, while D also feeds C. The gate E
is distinct from C because E is fed by B and C is not (an extra use would
contradict the degrees just obtained). Let z be E's other free input after
bypassing; it is already its other input before the trial. Instead fix z to
kill E and y to kill G. Now B is dead and x has no use. The distinct bits
x,y,z lose influence in two restrictions. The distinctions follow because
z is an independent pin of E, while x and y already have their two stated
uses at G,B and G,C. Thus the account is 144.

We may now assume that D does not feed both B and C.

### E5.3. C feeds D and D feeds B (or the symmetric chain)

Fix x to kill G and delete G,D,B. The only potential owners are y, now of
degree one, and C, whose degree does not fall when D is bypassed to B.
Neither becomes or feeds new troubled structure. This contradicts the assumed
failure of the three-gate route. By symmetry, after E5.2-E5.3 we may label
the bits so that B and D have no edge in either direction.

### E5.4. Locating a new troubled gate

Write B=B(x,I1), D=D(G,I2). Delete G,B,D after fixing x to kill G.
C3 shows that any newly troubled gate is either a predecessor I1,I2, or is
fed by a free bit passed from one of them. Since B and D are disconnected,
a gate predecessor retains its number of uses and its pins do not change in
this normalization. It therefore cannot newly be troubled. A new troubled
gate E must instead be fed by a free bit I1 or I2.

#### E5.4a. E receives the free bit z=I2 through D

The AND gate E has one use and was a consumer of D. Its other input is either
a free bit t already present, or B passing such a bit t after x is fixed.
Because z must end with two uses, there are only two degree possibilities:
z has degree two and D degree one, or z has degree one and D degree two.

If z originally has degree two, D has one use.
If E's other input t is direct, fix t to kill E. Delete E,D,G and a consumer
of E. The bits x,y,z each go from two uses to at most one, and D,G are dead
without other potential owners. Only the other predecessor of E's consumer
can add potential. Thus at least four gates go with dPhi<=1, giving 67.
If instead E's other input is B=B(x,t), then:
* For AND-type B, fix t to kill B and y to kill G; x becomes unused.
  Three influential bits disappear in two restrictions.
* For XOR-type B, impose x=t+c, choosing c to make B kill E. Remove
  B,E,D,G and a consumer of E. The bits y,z drop to at most one use;
  only t and the other predecessor of E's consumer can add potential.
  At least five gates go, dPhi<=2, and x is removed, giving 71.
The identification is legal because x is unprotected.

Now suppose z originally has one use and D has two.
If D is AND type, fix z to kill D. Delete D, the now dead G, and two
consumers of D. The dead G has no potential cost since x,y become one-use;
the other three deletions cost at most three. This gives 65.
If z is protected, fix x to kill G, delete G,B,D, then fix z to kill E.
Both consumers of D now receive z, so this second restriction also degenerates
D's other consumer and E's consumer. At least six constant/unary deletions
are obtained. Two influential bits and z's quadratic equation are lost, with
dPhi<=6, giving 129 in two restrictions.

The only remaining possibility is unprotected z and XOR-type D. Write
G=(x+a)(y+b)+c and D=G+z+e. Impose

    z=(x+a)(y+b)+c+e+h

with h chosen to kill E. This is a legal quadratic restriction: x,y,z are
distinct unprotected bits, and z's only occurrence is in D. D becomes constant
h; remove it, then the now dead G, E, D's other consumer F, and E's consumer H.
The first three removals have zero positive potential cost:
D's other predecessors are G and the removed z; G's free predecessors x,y
become one-use; E's other predecessor is either t becoming one-use, or B
fed by now one-use x. H costs at most one. The old troubled G disappears,
giving a decrease of one. If removing F has zero positive potential cost,
at least five gates go with net dPhi<=0. The new quadratic equation costs nine,
so the account is (5,1,-1,0), exactly 64.

Here is a complete resolution of the possible obstruction at F. To introduce
trouble when F is bypassed, F's other predecessor must be a free bit u and a
consumer A must be an AND gate with other free input v. Thus the original local
edges are D->F, u->F, F->A, v->A. (A is fed by F, rather than directly by D.)
The candidate troubled gate A is fed by u,v after bypassing F. We now choose
alternatives instead of the quadratic restriction.

If E gets t through B:
* When B has one use, fix y and z to kill G,D,E. B is dead and x unused,
  losing x,y,z in two restrictions.
* Otherwise B has exactly two uses and t has one, since t must end two-use.
  If B is AND type, fix x to kill B, leaving t unused; t is unprotected by E1.
  This loses two influential bits in one restriction.
* If B is XOR type and t unprotected, fix y to kill G, then impose
  x=t+c to fix B. The now sole occurrences of x,t are in B, so x,y,t
  lose influence in two restrictions.
* If that one-use t is protected, first fix its mate, which is distinct from
  x,y (both unprotected), then use the preceding alternative. Four distinct
  influential bits disappear in three restrictions.

If t directly feeds E, the bits u,v are distinct from x,y: after the proposed
quadratic restriction x,y are at most one-use, while u,v would feed troubled A.
They are also distinct from z, whose only original use is D. In particular
F is neither B nor C.
* If u=t, fix x and bypass B,D. E and F then have the same two inputs z,u,
  both unprotected: u originally feeds E,F, hence is unprotected by E1.
  C6 gives an affine line making both outputs constant. The bits x,z,u
  disappear in two restrictions.
* If u!=t and F has one use, fix v to kill A. Then F is dead.
  Delete A,F, a consumer of A and the other consumer of v.
  The last two named consumers may coincide; C3's frontier argument still
  supplies two deletions. No potential comes from A (its other predecessor
  is a gate) or F (D is XOR and u becomes one-use); the last two cost at most
  two. This gives (4,1,0,2), hence 66.
* Otherwise F has exactly two uses and u originally one. Indeed F cannot
  feed G,D,E, so none of its other uses disappears in the earlier proposed
  deletions; ending with two uses of u forces these degrees.
  If F is AND type, fix u to kill F and t to kill E. D becomes dead and
  z unused. Three influential bits disappear in two restrictions.
* If F is XOR type and u unprotected, fix t to kill E. D has only consumer F.
  Since F=G+z+u+c, impose z=u+c' and bypass F as G plus a constant,
  deleting D,F. The bits t,z,u lose influence in two restrictions.
* If that u is protected, first fix its mate. This mate is distinct from
  x,y,z,t,v: x,y,z are unprotected, and t,v feed AND gates, hence are
  unprotected by E1. Use the preceding construction. Four bits disappear
  in three restrictions.

This exhausts the obstruction at F, so the quadratic account or an alternative
always applies.

#### E5.4b. E receives t=I1 through B, and no bit through D

E's other input must be a free bit z already feeding it. If B has at least two
uses, it is XOR type by E3. Impose x=t+c, with c making B kill E.
Rewire G's x pin to t. This replacement does not introduce trouble at G,
since E is also killed and the trial's required two-use count for t is updated
on the same node. Remove B,E and the further constant-pin consumers: B has
another consumer and E has one. Thus at least four gates go, each costing at
most one potential unit (only t can add one at B). This gives 64.
If B has one use, fix z to kill E and y to kill G. B becomes dead and x unused,
losing x,y,z in two restrictions.

All possible primary-input AND gates have now been handled.

## E6. Cleaning a protected one-use XOR input

Let protected x have its sole use at XOR gate P. Fix x to either constant and
normalize. At least P disappears, one influential bit and one quadratic
equation are lost. If a second gate disappears, it is not a dead gate:
bypassing P preserves every other predecessor's uses unless P was the selected
output, in which case its predecessor becomes the selected output and is live.
Thus the first two removals are constant/unary, dPhi<=2, giving 65.

Suppose only P disappears. Its other predecessor R can add at most one potential
unit. If it does not, the decrement is at least

    5+48+9=62.                                                   (E.1)

If R is an AND gate, it cannot newly become troubled: its pins do not change,
and its use count does not fall on bypassing P, except in the output case.
That exception would leave another use of R unable to reach the output without
a cycle containing an AND gate; dead-gate trimming excluded it. Hence positive
potential requires R to be a free bit u of original degree at most two.

If u has degree one, P has two primary inputs x,u. The following elementary
alternative will also be used later:
* If u is unprotected, fix the mate y of x, then impose x=u+c to fix P.
  Three bits x,y,u lose influence in two restrictions.
* If u is x's protected mate, impose x=u+c directly. Its quadratic equation
  becomes affine, P is constant, and both bits lose influence in one restriction.
* If u is protected with a different mate v, fix both mates y,v and then
  impose x=u+c. Four distinct bits x,u,y,v lose influence in three restrictions.
Every identification uses a permissible unprotected bit, or identifies a
protected couple as allowed in C2.

If u has degree two, P has one use: otherwise passing u through P would leave
it with at least three uses, so it could not feed a newly troubled gate.
Let A be u's other consumer and B be P's consumer. They are distinct, since
otherwise P is useless. A does not acquire new pins and its use count does not
fall, so the new troubled gate must be B, an AND gate whose other pin is a free
bit t. Fix t to kill B; P is now dead and x unused. Then fix x's mate to make
x unprotected. Three distinct bits (t,x,its mate) lose influence in two
restrictions; t cannot be either x or its mate because their only uses have
already been specified and protected bits do not feed AND gates.

Consequently, either this case yields a block with average decrement at least
64, or every constant restriction on any protected bit removes exactly one
XOR gate, introduces no new troubled gate, and has decrement at least 62.
Call the latter operation a clean peel. Subsequent cases always pair its
possibly insufficient decrement with a quantitatively stronger operation.
A clean peel leaves the packing empty and makes the bit's mate unprotected.

## E7. A two-use primary input beside a one-use affine gate

There is a topmost AND gate G with one free input x of degree two and another
input Q which is a one-use gate. Here and below topmost means there is no AND
ancestor; Q is in the fair XOR prefix. The bit x is unprotected by E1 and cannot
also feed Q, since that would make Q useless. Fix x to kill G. Remove G,
the now dead Q, a consumer of G, and x's other consumer. C3 supplies the
four-gate count under coincidences. The dead Q costs at most two potential
units, G zero, and the last two at most one each. This gives (4,1,0,4), hence 64.

## E8. An affine-gate input to a topmost AND

If no AND remains, the output is affine, contrary to C2. Choose a topmost AND
gate G. By E4-E5 its inputs are not both free bits. Label them P,Q with Q a
gate. If P too is a gate, choose the labels so deg(P)<=deg(Q). Both gate inputs
are XOR-prefix affine functions of the free bits. They are nonconstant by
normalization. All references to affine support below mean actual nonzero
coefficients after solving the fair prefix, not merely syntactic paths.

### E8.1. Q depends on an unprotected bit x

Use C5's reconstruction to fix Q to its killing value for G, solving for x.
This removes x and keeps gate count and quadratic count unchanged during the
reconstruction itself, without increasing potential.

If Q originally has at least two uses, its new constant has at least three
constant-pin uses: the old consumers and the preceding node on the reversed
path. There is no self-pin before reconstruction (C5), so the new path consumer
is distinct from each old consumer. Along with a consumer of killed G, C3
gives at least four constant/unary deletions with dPhi<=4. The decrement is 64.
If Q has one use and P is a free bit t, E1 and E7 imply t is unprotected and
one-use. Killing G leaves t unused; together with x, two influential bits go.
They are distinct, since t's sole use is G and x reaches Q.

If Q has one use and P is a gate, the ordering makes P one-use too. After
reconstruction delete G, dead P, the reversed-path predecessor of Q, and a
consumer of G. The dead P costs at most two, G zero, and the other two at most
one each. These are four distinct deletions, or the same count by the frontier
argument: P cannot be a consumer of G (a nonlinear cycle) or the reversed-path
predecessor of Q (it had the sole consumer G and so could not feed Q).
This gives dPhi<=4 and decrement 64.

We may therefore assume, for EVERY topmost AND with the specified degree
ordering, that its higher-degree gate input Q depends only on protected bits.
For equal-degree gate inputs both depend only on protected bits, since their
labels can be exchanged. This universal assumption will rule out potential
obstructions at other topmost AND gates, not just at the selected G.

### E8.2. Only protected support at the larger-degree affine input

Choose a protected bit x_j in Q's nonempty support and let x_k be its mate.

#### E8.2a. A gate input has a protected couple as its two pins

If P or Q is fed by two coupled free bits, both are one-use by E1-E2.
An identification with the appropriate sign fixes that XOR gate and dissolves
their quadratic equation. Both bits cease to be influential. This gives 96.

#### E8.2b. Q has at least two uses

Clean-peel x_k, losing at least 62 and leaving x_j unprotected.
If x_k does not directly feed Q, Q remains a gate with the same uses and
continues to depend on x_j. Its coefficient of x_j is unchanged by fixing
a distinct bit. If x_k feeds Q, that deleted gate is the unique gate removed
by the clean peel; its other predecessor Q' passes to all its old consumers.
The node Q' is a gate depending on x_j and has at least as many uses as Q:
it loses one use to Q and gains Q's at least two uses. It cannot be a
primary input, since that would be x_j and would put the coupled pair at Q,
already covered by E8.2a. Use Q' in place of Q.

Reconstruct this surviving gate to kill G, solving for x_j. At least four
constant/unary gates disappear, as in E8.1. Name B the first gate on the
reversed path, H a consumer of G, and the other constant consumers as needed.
If deleting G and B creates no positive potential, the other two credited
deletions cost at most two, giving a second-step decrement of at least 66.
Together with the clean peel this is 128 in two restrictions.

The exceptional possibilities of positive potential at G or B are not assumed
away: they are proved and resolved in E8.2e below.

#### E8.2c. Q has one use and P is a free bit t

The bit t is unprotected and one-use by E1,E3,E7. Clean-peel x_k and reconstruct
the surviving Q (or its passed predecessor Q') to kill G. The bits x_k,x_j,t
are distinct and all lose influence in two restrictions, giving at least 144.

#### E8.2d. P and Q are one-use gates

Both supports consist only of protected bits, by E8.1. We first consider whether
a bit in one support has its mate NOT directly feeding the other gate.

Suppose x_j belongs to Q's support and its mate x_k does not directly feed P
(the symmetric choice is identical). Clean-peel x_k. P is not the removed gate
and stays one-use. Reconstruct Q, or its passed predecessor, to kill G.
As in E8.1 this deletes four gates, with potential at most four. If potential
is at most two, the two-step decrement is at least 62+66=128. Positive potential
at G or reversed-path B is deferred to E8.2e. Assume neither contributes.
G's consumer contributes at most one. Therefore potential at least three
requires both predecessors of dead P each to contribute one.

A gate predecessor T cannot do so: it reaches G and hence is in its XOR prefix,
so is not AND type. If P's two predecessors were already free bits before
peeling, both are protected one-use bits (P depends on each), and neither can
create troubled structure. Thus just before removing P its two pins are free
bits t1,t2, at least one passed through the reversed-path B or the sole peeled
gate P'. The pins are distinct; repeated pins or a constant XOR would have
been normalized.

If either t1 or t2 is protected, after the two-step normalization it feeds a
new troubled gate or is the shared bit of a new packed pair. Fix it to kill
one such gate. The two bit uses and a consumer of the killed gate give at
least three constant/unary deletions, costing at most three potential units,
and its quadratic equation is lost. This third restriction gains at least 69.
The total is at least 62+64+69=195>=192 in three restrictions.
The third bit is distinct from x_k,x_j: the former was fixed, the latter became
an affine variable during reconstruction, whereas the third is still protected.

It remains to exclude the configuration with both t1,t2 unprotected.
If both were passed, one must have come through B and the other through P',
since only those gates have been bypassed into P. Relabel so t2 comes through
B and t1 through P'. Before peeling, the relevant equations are

    Q=B+t2+c1,       P=x_k+t1+B+c2.

Here B is a gate: otherwise B=x_j would be a protected bit with uses at P,Q,
contrary to E1. Since Q has no unprotected support, B depends on t2.
The distinct free bit t1 cannot cancel t2, nor can protected x_k, so P also
depends on t2, contradicting its protected-only support.

Finally suppose t1 directly fed P before peeling. If t2 was also direct or
was passed only through P', the equation of P makes its coefficient on t1
one, contradicting E8.1. Otherwise t2 came through B, possibly also P'.
Immediately after peeling, B feeds both P and Q, and t2 feeds Q.
Since only P' was removed, at least one of B->P or B->Q was a direct edge
before peeling. If B->P was direct then P=B+t1+c in that original state.
Either P depends on unprotected t1, or B does and Q therefore does as well;
the other possible extra term in Q is x_k or the distinct t2, so cannot cancel
t1. If B->Q was the direct edge, exchange the argument with P and Q.
This contradiction excludes both unprotected pins.

It follows that, after the successful cases above are excluded, every bit in
P's support has its mate directly feeding Q, and conversely for Q's support.

If one gate, say P, has two primary inputs, they are protected and one-use.
They are not coupled, by E8.2a. Fix both to make P kill G. Remove P,G, dead Q,
and a consumer of G. Potential costs are respectively zero, zero, at most two,
and at most one. Two distinct quadratic equations and two bits are lost:
the account (4,2,2,3) gives 131>=128 in two restrictions.

If both P and Q have at most one primary input, pick x_j in P's nonempty
support. Its mate x_k directly feeds Q by the conclusion above. Q's other
pin R is a gate, for otherwise the preceding case applies. It is nonconstant.
The one-use path x_k->Q->G cannot reach R without a cycle containing the AND
gate G; hence the reachability fact in C5 says R is independent of x_k.
Thus Q depends on x_k. Choose a bit x_i in R's nonempty support. It is also
in Q's support, is different from x_k, and different from x_j: x_j's sole
use is P, whose sole use is G, so it too cannot reach R. Therefore Q depends
on two protected bits x_k,x_i which are not coupled. Their two distinct mates
must directly feed P, contradicting that P has at most one primary input.
This exhausts the ordinary branch of E8.2d.

#### E8.2e. The deferred potential obstructions

For clarity retain the names of the initial state before the clean peel:
x_j is in Q's support, x_k its mate, and P' the sole XOR gate removed by fixing
x_k. After that peel, Q either remains, or its other predecessor replaces it.
In this subsection identify the surviving affine gate with Q. Reconstruct Q
to its killing value for G, solving for now unprotected x_j. Let B be the
preceding node of Q on the chosen path, now the first gate on the reversed
path; B may be the new gate representing old x_j. Delete G first, then B,
then the other constant-pin gates. We examine every way G or B can contribute
positive potential.

(1) Positive potential at G.

Since G was topmost, an AND predecessor is impossible. Its other predecessor
must be a free bit x_m which goes from degree three to two when G disappears
and then feeds a troubled AND T. Immediately after peeling it therefore has
degree three. Let its degree before peeling be h. Only P' was bypassed;
the possible values are h=1,2,3.

If h=3, it could not directly feed G originally by E3, and so must have passed
to G through P'. To become troubled at T it must also feed T after peeling.
A direct original use at T contradicts E3 again. If it was passed to T through
P', P' had at least two uses, increasing x_m's degree to at least four rather
than three. Thus h=3 is impossible.

If h=1, the only use of x_m and the only use of protected x_k were the two
pins of XOR P'. The elementary one-use alternatives of E6 apply to that pair,
giving a block of one, two or three restrictions with at least 64 per
restriction, according to whether they are coupled or separately protected.

If h=2, P' has exactly two uses and pins x_m,x_k. Here x_m is unprotected by
E1. Either x_m directly fed G before peeling, or it was passed through P'.
In the direct case, P' must already have the AND consumer T. Fix x_m and x_k to choose P' as
T's killing value, also choosing x_m's value to kill G; the XOR equation
then uniquely chooses x_k's value. Thus P',T,G become constant. With their
constant-pin descendants, at least six constant/unary deletions follow by
C3. Two bits and x_k's quadratic equation are lost, dPhi<=6, giving 129.
The killing choices are compatible because x_m first fixes G and x_k
independently controls P'.

If x_m was passed to G through P', and the original Q had degree two, then
P' is also degree two and depends on unprotected x_m. E8.1 applies at the
original G with the two affine gate inputs exchanged, giving 64 directly.
If the original Q has degree at least three, perform the clean peel (62),
then reconstruct Q. There are at least four constant-pin consumers of the
new constant, plus a consumer of G, giving at least five constant/unary
deletions with dPhi<=5. This gains at least 25+48-5=68. Together the
two restrictions gain at least 130>=128.
These cover all positive-potential possibilities at G. We may assume its
deletion contributes no positive potential for the remaining analysis.

(2) Positive potential at B when B->Q was an original edge.

After reversing the XOR equation, B's other pin must be a free bit x_l passed
from Q's other pin; a gate pin cannot become troubled by topmostness.
Originally x_l either directly fed Q or was passed there through P'.
The original affine equation is correspondingly

    Q=B+x_l+c,       or       Q=B+x_k+x_l+c.                     (E.2)

If x_l is unprotected, E8.1 says Q is independent of it, so B depends on x_l.
If B has no uses immediately before deletion, it must be the new gate Z for
old x_j: any other reversed-path gate feeds the next gate on that path.
Before reconstruction B was then the primary bit x_j, making (E.2) depend
on x_l, a contradiction.

Thus B has at least one use. Bypassing it does not lower x_l's use count;
new trouble must be an AND gate T receiving x_l through B, with other free
pin t. There are only the following original-edge possibilities, since P'
is the only gate removed before reconstruction.
* If B and t already directly fed T, T was a topmost AND with affine gate
  input B depending on unprotected x_l. Its other input is a primary bit,
  so E8.1 applies there, contradicting the universal assumption.
* If B was passed to T through P', t was direct. The XOR gate P' depends
  on x_l too, because its other input x_k is protected and distinct from
  x_l. Apply E8.1 to T with gate input P'.
* Otherwise B was direct and t was passed through P'=t+x_k+c.
  If t originally had degree one, E6's one-use alternatives apply to t,x_k.
  If t originally had degree two, P' must have had exactly one use to
  end with the two uses needed for trouble. B is a gate of degree at
  least one and depends on unprotected x_l. Thus at topmost T with
  gate inputs B,P', B can be the degree-ordered larger input (including
  equality); E8.1 applies again.
No larger original degree of t is possible, since passing it through at least
one use of P' would leave at least three uses.

If instead x_l is protected, after the reconstruction/normalization it feeds
a newly troubled gate or is the inner bit of a packed pair. Fix x_l to kill
such a gate. There are three constant/unary deletions: the killed gate, its
consumer, and the other consumer of x_l (with coincidences counted by C3).
The first has zero positive potential cost since its other free input goes
from degree two to one. The latter two cost at most two. The old troubled
gate or packed pair is destroyed, reducing potential by at least one.
Hence this third restriction gains at least 15+48+9-1=71.

For this protected alternative the first step was at least 62.
The reconstruction removes at least B,G and a consumer of G; since G has
zero potential cost, the potential increase is at most two, giving at least
15+48-2=61. Summing all three gives 62+61+71=194>=192.
If more normalization gates occur their nonnegative gains only help.
The protected x_l is distinct from fixed x_k and affine x_j, so its quadratic
equation is still a separate equation when the third restriction is applied.

(3) The B->B'->Q edge contraction.

There may have been one intermediate gate B' between B and Q before peeling,
where B'=P' was the sole peeled gate. Its equation was B'=B+x_k+c.
Immediately after the peel the edge is B->Q, and the other pin x_l directly
fed Q originally. Its original equation is Q=B+x_k+x_l+c.
If x_l is unprotected, Q's coefficient on it is zero by E8.1, while B's is
one. B cannot be dead before its deletion, for then it was the primary x_j
before reconstruction and could not depend on the distinct x_l. For a
troubled T reached through B, its edges can only be direct or pass through
the same sole P'. The three alternatives in (2) apply, using this displayed
equation for the dependency contradiction. If x_l is protected, the same
third restriction and account 62+61+71 apply. Thus contraction does not
leave an additional case or an unproved hypothesis.

(4) The B->Q'->Q contraction when Q itself was peeled.

Here the original Q had pins Q',x_k, and Q' had pins B,x_l.
After peeling, Q' takes Q's role; the original equation is again
Q=B+x_k+x_l+c. The preceding argument is identical with the surviving Q'
named Q: its support on x_j remains nonzero, original Q is independent of
every unprotected x_l, and the only removed gate is P'=old Q. A potential
owner is either a protected x_l, handled by the third restriction, or an
unprotected x_l with B's coefficient one, handled by the same direct,
through-P', or one-use-pair alternatives of (2).
This lists both possibilities for the one removed XOR vertex on the last
edge of the path; no further contractions are possible in a clean peel.

E8.2e resolves every positive-potential obstruction that was deferred in
E8.2b or E8.2d. Hence those cases provide their ordinary two-step 128 account
or one of the explicit alternative blocks.

## E9. Exhaustion and conclusion of the local assertion

Normalization and packing are handled first. E1-E2 leave protected bits
one-use at XOR gates; E3-E5 remove every possibility of an AND with two free
inputs. There must then be a topmost AND with an affine gate input, since an
all-XOR output is excluded by C2. E6 provides either a successful block or
the clean-peel property used thereafter. E7 handles the exceptional
two-use primary input beside a one-use gate. For the degree-ordered topmost
AND, E8.1 treats unprotected support; otherwise E8.2 handles coupled pins,
larger gate degree, a primary other input, and two one-use gate inputs.
The latter case ends either in an explicit restriction block or in the
contradiction that a gate with at most one primary input must have two.
E8.2e supplies all postponed potential exceptions.

Every chosen truth-table value is either an AND killing value, a solved XOR
equation, or the affine line proved in C6. Quadratic substitutions occur only
in E5.4a, with three distinct unprotected bits and no nonlinear cycle.
Reconstructions occur only in the fair XOR prefix and solve for an unprotected
bit by C5. Protected affine identifications only identify an existing couple.
All other operations are constant restrictions. They preserve C2's source
class, remove exactly one free dimension apiece, and use at most three
restrictions in a block. The accounts above prove Delta M>=64 times the number
of restrictions. This proves C7's local assertion for the full B2 circuit
model with arbitrary sharing and depth.
