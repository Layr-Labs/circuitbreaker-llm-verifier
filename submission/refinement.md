# A sharp additive constant for d>=3

**Theorem.** For every n and d>=3, if a total Boolean function f on n
inputs is nonconstant on every dimension-d affine subspace, then

    C_B2(f) >= 4n-8d+8.                                 (T)

The constant8 cannot be increased uniformly over this range. The example
xy+uv+z at n=5,d=3 has circuit size4 and no constant three-flat, as proved
in proof.md Appendix B.3. All B2 tables, fanout, sharing, constants and
charged affine rows are included. This is an additive refinement, not an
increase in the leading coefficient.

Use the complete reduction in elimination.md Sections1--5.2. Its state
has independent parameters, disjoint promised pairs, an invertible sparse
affine system and an acyclic part, with correctness on the entire promise
under injective affine Phi. Write a=m-q, s for the number of charged
objects, and E for the distinct permanently inert coordinates counted by
N0,NB,N1,N2. Initially a=n; stop immediately at a=d. No inert coordinate
is a designated output. The proved invariants include E<=d-1 and E+q<=d
at termination. Ordinary deletions are counted by K. Define Lambda as
the sum of the nonnegative integer slack in the complete J=s-2q step
bounds of elimination.md Sections3--4, including any additional deletions.
Affine moves and role exchanges have zero J change. The exact ledger is

 R=|C|-(4n-8d)
  =4(d-E)-2q_T+s_T+2N0+3NB+N1+K+Lambda.                (11)

For used paired cleanups A and C, zero slack means exactly two removed
objects. For D,E1,E2 it means exactly2,3,4 respectively. These are counts
of distinct original objects, including propagated removals. They do not
assume J nonnegative. The proof of the three-object E1 bound and all
sharing arguments are in the supplied elimination proof.

The cases n<=d give a nonpositive target and are immediate. The d=3
case follows from the complete six-variable Ramsey proof and circuit
argument in Appendix B. Hence assume n>d and d>=4 and suppose R<=7.
We will contradict every possible terminal state. No preceding +6 or
+7 theorem is assumed in this argument.

## 1. Geometry and object-count facts

We use Lemma R of elimination.md Section5.2 with its regenerated
unrestricted-coloring refutation. Its three pairs and two free coordinates
and distinguished AND-type gate force a legal constant three-flat,
regardless of the rest of the computation. The elementary facts already
proved there give

* E=d-2 implies a<=d+1;
* E=d-2,a=d+1 implies q<=1;
* the first reach of E=d-1 is the final paid E1 step.

These facts apply to intermediate states and require no tightness of the
object ledger. One may count an additional computationally unused
unpaired parameter as an extra inert direction in these arguments.
Appendix B adds a useful stronger fact:

    E=d-3 implies a<=d+2.                              (12)

Otherwise choose a legal line in each pair and fix any surplus free
parameters, leaving an affine six-cube of effective units. Its coloring
has a constant three-flat. With E inert directions this is a forbidden
d-flat. The same reasoning applies with an additional unused unpaired
direction counted among the E directions.

We also use three small essential-input observations. A function with no
legal constant three-flat on the following domains depends essentially
on every displayed coordinate:

    one pair and four free coordinates (six coordinates);
    two pairs and three free coordinates (seven coordinates);
    three pairs and one free coordinate (seven coordinates).

For an ineffective paired coordinate, fix its partner to zero and leave
that coordinate free. The remaining domains in the three cases have
respectively16,24,18 points in dimensions4,5,5. For an ineffective free
coordinate, the remaining domains have respectively24,36,27 points in
dimensions5,6,6. In each case the pairs in a same-color majority outnumber
the available nonzero differences: 28>15,66>31,36>31,66>31,153>63,
91>63. A constant plane plus the ineffective direction is a constant
three-flat. This proves all three observations.

A fair affine system and acyclic gates with T globally essential parameter
inputs have at least T-1 charged objects. Choose a nonzero determinant
term of the invertible auxiliary matrix and match each row to a distinct
auxiliary occurring in it. Regard that auxiliary as the row's output.
At most two other variable terms remain as parents. This orientation
may be cyclic; it inserts no free computation. In the undirected
component of the designated output every essential parameter occurs.
Other components have independent auxiliary systems, which are uniquely
solvable because the auxiliary matrix is block diagonal by components
and invertible. If the output component contains s' charged objects, it
has s' matched internal vertices and at least T parameter vertices.
Connectivity needs at least s'+T-1 distinct edges, while its two parent
pins per object give at most2s'. Thus s>=s'>=T-1. Constants strengthen
this inequality. For acyclic circuits the same pin argument counts all
used inputs once dead gates are removed, even without essentiality.

Finally at most two charged rows can be replaced by at most two ordinary
XOR-type gates. Match their auxiliary columns as above. The matched
2-by-2 matrix has diagonal ones. Both off-diagonal entries cannot be
one, since that would be singular. Thus the matched dependency is
acyclic; each row defines its output from at most two parameters or
earlier outputs. Constants and complements are absorbed in full-B2
truth tables. The one-row case is immediate. This counts every row.

## 2. Three finite formula lemmas

**Q.** Every six-gate read-once formula on seven independent coordinates,
with two promised pairs and three free coordinates, has a legal constant
three-flat. The complete finite proof is as follows. Normalize B2 gates
to AND/OR/XOR with signed leaves by De Morgan and XOR identities; ignore
a global output complement. Translate the three free coordinates to
remove their signs and enumerate all four paired signs. Enumerate all
10395 unordered labeled seven-leaf binary trees. The two pair flips,
pair exchange and free-coordinate S3 give a48-element affine symmetry
group. Canonicalization leaves350 representatives. Each has16*3^6=11664
operation/sign choices. Generate all72 legal points and all5049 contained
three-flats by the complete RREF procedure. The retained witness cover
has3586 plane IDs and7528 decoded bytes. Exact masks check coverage of
all4082400 representative/assignment combinations. All IDs and lengths
are bounded and the decoder checks exact end of data. These steps and
the elementary normalization supply a proof, not a random search.

**T7.** Consider a full-B2 circuit on seven coordinates in the same
two-pair/three-free domain. Suppose it has at most seven gates, every
paired coordinate occurs once, and it has a distinguished AND-type gate
G on literal coordinate inputs. Then it has a legal constant three-flat.
A circuit with at most six gates is either already reducible or, since
all seven coordinates are essential by Section1, has exactly six gates
and is read-once by the pin equality; apply Q. Remove dead, unary and
constant-pin gates first. For seven gates, seven input vertices and six
non-output gates need13 uses among14 pins. Thus exactly one vertex B
has two children and all others one. The cone of B is a read-once tree,
and the outside computation is a tree with two occurrences of its value.
An ell=6 cone after contracting G would leave only a parent reading its
same value twice; remove that degenerate gate and use Q.

Several G placements have direct proofs. If its two inputs are free,
fix a forcing value of an input, choosing the repeated input if one is
repeated. The other free input is now unused. The remaining two pairs
and one free coordinate form18 points in F_2^5, forcing a constant plane;
adjoin the unused direction. If its inputs are partners, choose a legal
pair line making G constant. Its line parameter is unused and the
other pair plus three free coordinates have24 points in F_2^5, again
forcing a plane. Paired input uniqueness justifies unusedness, including
when G itself is shared.

For cross-pair inputs x,y, if alpha=0 in G=(x+alpha)(y+beta)+gamma,
fix x=0 and hy=0. Then y is unused and the residual four free coordinates
force a plane. The beta=0 case is symmetric. For a mixed input (free x,
paired y), remove the free-input sign by translation. If beta=0 and x
is not repeated, fix y=0: x is unused and the other pair with three free
coordinates forces a plane. The remaining orientations and repeated-x
cases are exhaustively checked below. Complements of G's output are
absorbed in the subsequent tree normalization.

Contract G into component0. In the cross case the other components are
its two partners1,2 and the three free coordinates3,4,5; G is the bad
forcing NOR of its two paired coordinates. In the mixed case they are
its partner1, the other pair2,3 and free coordinates4,5; G=x*(y+1).
The single repeated vertex is either a free input or G itself, or an
internal read-once cone on a subset of these components. Enumerate every
subset of size ell=1..5, every inside tree, and every outside tree on the
remaining components and two copies6,7 of the cone. A paired component
cannot be repeated. In the mixed case there is also the possibility
that the repeated vertex is the original free input x of G, occurring
once outside G. Enumerate every seven-component outside tree for this
case and both paired-input forcing orientations of G. These exhaust the
single-extra-use graph, since any other sharing would require a second
extra pin.

The cross group is the partner exchange and free S3, of order12. The
mixed group exchanges the other pair and the two other free coordinates,
of order4. First choose the least subset in each group orbit, then use
its set stabilizer to canonicalize the two trees. This is equivalent to
quotienting by the whole group: two graphs on that chosen subset can
be related only by a map preserving the subset. The additional repeated-x
case is canonicalized directly. The counts are4208 and12468 structures.
The checker constructs them; it does not trust a supplied topology list.

Here is the complete sign normalization for the shared cone. Original
free-coordinate signs can be removed by translations; paired signs
remain enumerated. If B contains G, push the sign at G through its unique
path inside B to B's output, using De Morgan/XOR and off-path subtree
negations. This fixes the G sign and toggles both B occurrence signs.
Then push the first B occurrence sign through the outside tree to the
global output. Off-path negations affect only the other B sign and
original leaves outside B; these are independently enumerated or free
translations. If B does not contain G, first fix G's sign through the
outside tree, treating B occurrences as signed atoms. Then complement
the inside B computation if necessary to fix its first occurrence sign;
its leaves do not contain G. De Morgan/XOR keeps every gate in
AND/OR/XOR form. An input B is handled by the same outside path if it
is G, or by translating its free coordinate. Thus only the relative
second-B sign remains, in addition to the paired signs.

The cross masks use three sign bits and six operation trits. Mixed masks
use four sign bits and six trits, except repeated-x outside G, which
uses five sign bits, including its extra leaf sign and both G forcing
orientations. That extra x sign is retained because its coordinate is
already used to normalize G's input sign. Every original possibility
has a normalized representative; no independent input signs of the two
shared cones are falsely assumed. The checked totals are24541056 and
181643472 structure/assignment combinations.

Generate the72-point domain and all5049 three-flats once. For the mixed
case apply the explicit coordinate permutation (1,3,2,4,0,5,6) to the cross
domain and its flats. The audit evaluates all masks and verifies their
flat coverage. Inside-tree and bottom-gate values are cached only when
the sign layout, operation position and actual inputs agree. Caching
changes no enumeration or tested assignment. This proves T7.

**T8.** A seven-gate read-once formula on eight independent coordinates
in four promised pairs, with a distinguished AND-type gate on coordinates
from distinct pairs, has a legal constant three-flat. Good forcing
orientations are analytic: fix x=0,hy=0, leaving y unused and an18-point
two-pair/one-free residual domain forcing a plane. In the remaining NOR
orientation, contract G, enumerate every tree on the remaining seven
components, all64 paired signs and six operation trits. The partner swap,
the two other pair flips and their exchange form a16-element group.
Its880 representatives have41057280 assignments altogether. Generate
all81 legal points and2106 contained three-flats and check full mask
coverage. This is the same complete tree/sign normalization as above,
without sharing. The supplied witness data give582 exceptional flat
choices after the bounded common-flat scan. It proves T8.

All three finite lemmas use standard-library exact arithmetic. The
RREF generator enumerates every contained affine flat once, using zero
pivot-offset bits and every reduced free entry. The accompanying
common-flat mechanism and data validation are specified in Section5.

## 3. Classify the small exceptional tail and its predecessors

**F.** A five-gate read-once formula on six independent coordinates,
one promised pair and four free, either has a legal constant three-flat
or is equivalent to

    OR(x,u) XOR OR(y,v) XOR OR(w,z), with x*y=0.          (13)

Equivalence uses only free-coordinate complements/permutations, a swap
of the paired coordinates, and output complementation. The stronger
fact needed below is gate-level: there are three disjoint two-leaf
AND/OR cherries joined by two XOR gates, up to internal complements.

To prove F, normalize every gate to AND/OR/XOR with signed leaves.
Enumerate all945 six-leaf trees and quotient by the pair exchange and
free S4 (order48), giving41 representatives. For each enumerate four
paired signs and3^5 gate choices; free signs are removed by translations.
Generate all48 legal points and1770 legal three-flats. Exact mask coverage
leaves exactly16 assignments on two tree shapes, with eight legal truth
tables. The auditor decodes all16 assignments and checks that their
cherries are AND/OR and both joins XOR. It also evaluates the FULL64-point
truth tables and verifies their membership in the32-element orbit of
(13) under the four free complements and output complement. Thus the
claim is stronger than equality merely on the promise. In each paired
cherry its free-input linear coefficient is one; complementing that free
input and, when needed, the cherry output turns it into OR. The free/free
cherry is likewise normalized. The XOR joins absorb cherry output signs.
These identities give the stated gate-level equivalence. Its two tree
associations are both enumerated in the predecessor checks.

We now record three exactly checked predecessor facts. The exceptional tail is normalized as in F. In each of the following
three predecessor classes the indicated constant flat is unavoidable.

**A predecessor.** Add a parameter p paired with one of the four free
tail coordinates h. Reverse a used paired-A cleanup p=0 that deletes
exactly two objects and leaves the five-gate exceptional tail. There is
a constant three-flat in this two-pair/three-free predecessor domain.
For two acyclic deleted gates, the input p has exactly two children;
more would give extra deletions. Pruning its two occurrences leaves
the read-once tail. Conversely graft two copies of p at arbitrary edges
of that tail. Normalize each deleted gate's forced output to its old
child at p=0, compensating the parent truth table. An old-containing
graft then has truth table4,6,12 or14 (identity at p=0). A new-only p,p
subtree needs only the five nondegenerate zero-preserving tables2,4,6,8,14:
on p,p these already give both possible values0 and p. All positions are
enumerated, including degenerate identity padding12.

Rows cause no missing case. There are at most two extra rows, and their
matched auxiliary system is triangular. Independent rows Q=p+r and
R=p+s give XOR grafts. In the triangular case Q=p+R,R=p+r, substitution
gives Q=r and R=p+r. Both cannot have independent uses in the tail after
p=0, since that would repeat the same old input and violate the six-input,
five-gate pin equality. Keeping at most one requires at most one B2
graft, padded by identity. With one row Q=p+r and one p-gate using Q,
their used output is a two-input Boolean function of p,r, again at most
one graft plus padding. If that gate reads another old value, there are
two ordinary grafts. An acyclic gate cannot feed a row. Constants are
absorbed in truth tables and signs. This exhausts zero-, one- and two-row
cases; unused outputs only simplify them.

Do not translate h before restoring its pair. Instead the auditor retains
an extra h leaf-sign bit, as well as all16 normalized tail exceptions.
The other free-coordinate complements preserve the prior promise.
For each of four choices of h and both tail associations there are77
distinct two-graft structures. Their616 cases enumerate163328 assignments,
including identity padding, and every one has a checked constant
three-flat. Thus the row simplifications and all acyclic possibilities
are covered with their charged counts.

**C predecessor.** Reverse a paired-C cleanup: p is read only in its
row Hp, its partner h only in Hh!=Hp; h=0 deletes Hh and Hp, with the
permitted affine role exchange. The result is the same five-gate tail.
There is a legal constant three-flat before this cleanup.

The surviving tail has no rows. If Hh is acyclic, the single Hp solves
Q=p+r+c, where r is another parameter or zero. Hh has h as one pin;
when h=0 it must reroute to a nonconstant old value, otherwise extra
propagation would delete a tail gate. Its output has one tail use, since
two would repeat an old value after rerouting. It is therefore one h
identity-at-zero graft at an arbitrary tail edge, including the Q edge.
The old Q leaf is p+r+c. No row can read an acyclic gate output.

If Hh is a row, both extra objects are rows. Both cannot contain both
auxiliaries: their2-by-2 auxiliary matrix would be singular. Up to names
they are independent Q=p+r,R=h+t, or triangular. Independent rows give
the preceding model with an h-XOR graft on t. For R=h+r,Q=p+R, either R
has no tail use (an h-XOR graft on Q), or its unique tail use is the old
r leaf. The latter model is

    old r -> r+h,  old Q -> p+r+h.                     (14)

For Q=p+r,R=h+Q, at most one of Q,R has a tail use after h=0, by the
same read-once pin equality; this is the single-graft model or an
ineffective h, which already gives a plane plus its unused direction.
These are all invertible two-row possibilities. Row constants are
absorbed in the old paired signs, a Q leaf sign, and the permitted other
free complements; no forbidden translation of p is made.

The checker chooses each of four free Q slots, either tail association,
r from the other five parameters or zero, an extra Q leaf sign, and
all identity-at-h=0 graft tables4,6,12,14 at its eleven edges. It also
checks (14) for every nonzero r choice. Old paired signs are applied to
old tail leaves, not incorrectly to the linear row computation. These
96 packed cases contain34432 assignments and all have a checked legal
constant three-flat.

**E2 predecessor.** A nine-gate read-once formula on ten coordinates,
three pairs and four free, whose tight E2 deletes four objects to leave
the exceptional six-coordinate tail, has a legal constant FOUR-flat.
The two G inputs x,y and their partners hx,hy are the four removed or
newly inert coordinates. Contract G into atom A and prune A,hx,hy.
Because all six old coordinates remain essential, no old-containing
branch can be killed by a constant; it reroutes to that branch or its
complement. Thus the tree is obtained by grafting A,hx,hy at arbitrary
edges of one of the two tail associations, in that fixed order. Removing
these leaves in reverse order proves completeness of the graft construction.

Write G=(x+alpha)(y+beta)+gamma, absorb gamma in its uses and fix
x=alpha,hx=hy=0. Then A=0. Normalize each new-only subtree's forced
value to zero by complementing its output and compensating its parent.
Its gates are among the five zero-preserving nondegenerate tables.
Each new gate with an old-containing child must be identity at its other
child's forced value zero; after its output sign is absorbed its tables
are4,6,14. The internal signs making the old tail canonical are harmless:
conjugating an intervening graft D(a,b) to1+D(1+a,b) preserves its
identity-at-b=0 property and permutes these three tables. Original free
complements and permutations are legal before E2 too. This proves the
full gate-level normalization, not replacement by an equal output truth
table in a different surrounding circuit.

All four G forcing orientations, all permitted new tables and all graft
positions are enumerated. Swapping the old pair together with its free
cherry partners, exchanging the other two free inputs, and exchanging
hx,hy give an eight-element symmetry group. Canonicalization leaves1450
structures and188616 assignments. The auditor checks a constant four-flat
for each assignment using360 supplied flats, each independently validated
by its offset and four independent basis vectors, sixteen distinct legal
points and actual output masks. No assertion that this dictionary exhausts
all four-flats is needed. Its1195 exceptional choices supplement the
bounded common-flat scan. This proves the E2 predecessor fact.

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

## 5. Certificate execution and exact coverage

Every finite formula statement above is checked by the declared audit.
The package contains two bounded JSON certificate containers. Each blob
uses canonical ASCII base85 and bounded XZ decoding with a16-MiB memory
limit, followed by a checked magic header, counts and integer indices. All decompression is bounded, integers and array accesses are
range checked, repeated variables or invalid flats are rejected, and
exact end of data is required. The reader regenerates the topology
representatives and domains and never imports research scripts or a SAT
solver. The exceptional-tail classification supplies its own sixteen
cases rather than trusting a stored list.

To keep flat witnesses small, each new cover has a dictionary of64
valid common flat IDs. The T7 covers apply these masks in order, stopping
if coverage is already complete, then apply each exceptional supplied
mask. Every exceptional mask must cover a previously uncovered assignment
and their final union must be full. For the other covers, use this
bounded adaptive procedure: for a still-uncovered parameter-assignment bit,
scan these64 in order and take the first actual monochromatic flat
not already used. If none qualifies, consume one exceptional supplied
ID and check that it really covers this bit. Apply the flat's exact
monochromatic mask to all assignments at once. Each common ID can be
used once; every supplied exceptional ID must be consumed, used once,
and actually cover the chosen bit. The number of iterations is at most
64 plus the exceptional list length. Coverage ends only at the all-ones
assignment mask. The common list is untrusted witness data, not a
mathematical assumption or solver search. Bounded scans over explicitly
regenerated affine flats and exact bit operations check every enumerated
function. Root normalization and universal circuit reductions are the
supplied mathematical arguments above, not inferred from finite tests.

## 6. Exhaust the terminal states with R<=7

For q_T=0, if E<=d-2 equation(11) gives R>=8. Otherwise E=d-1 and
N1>=1. The remaining budget is

    s_T+2N0+3NB+N1+K+Lambda <=3.

Thus NB=0. If N1=1,N0=0, all prior inert events are E2. The last one
reaches E=d-2,a=d+1,q=1, as the only remaining paid step is E1.
Immediately before it E=d-3,a=d+2,q=3. Paired cleanup is exhausted,
so its distinguished G has no other dependence on its two inputs.
Lemma R gives a constant three-flat, contradicting dispersal, regardless
of the budget's unused ordinary losses or slack.

If N1=1,N0=1 the entire budget is tight; s_T=K=Lambda=0. If the
penultimate inert event is E2, the preceding argument applies. Otherwise
it is the sole N0. Immediately before N0, E=d-3,a=d+1,q=2,s=3: after
N0 the three-unit fact forces q=1 and the final E1 deletes three gates.
Zero slack in N0 means both members of its unused pair were already
computationally unused; fixing a used partner would have removed an
additional object. That pair cannot have been unused in an earlier
completed iteration, because unused-pair cleanup is exhausted first.
Its two members became unused at the immediately preceding removing
step; affine preparations preserve their occurrences.

If that step was D, the preceding state has E=d-3,a=d+2,q=1,s=5,
all five objects acyclic. It is exactly the forbidden history state in
Section4. If that step was E2, the prior state has E=d-4,a=d+2,q=4,s=7,
all acyclic. All eight paired coordinates are used. At least one of its
two other active coordinates is used: otherwise two extra inert directions
would give E'=d-2,a=d+2, contradicting the four-unit fact. Nine inputs
and six non-output gates need fifteen pins, more than fourteen. If the
step was A or C, the prior state has E=d-3,a=d+1,q=3,s=5. Its seven
active coordinates are all essential by Section1; charged connectivity
requires at least six objects. All possible preceding steps are excluded.

If N1=2, then N0=0 and the other budget terms sum to at most one. If
the penultimate inert event is E2 use Lemma R again. Otherwise the last
two inert events are E1 steps, with no paired cleanup or D between them:
only one pair remains for the final E1. Before these steps a=d+2,
E=d-3,q=2, with at most seven objects. At least six guaranteed deletions
are acyclic; there is at most one row. Replace that row by a counted
XOR gate. All seven coordinates are essential, paired coordinates are
unique, and the selected minimal AND is on coordinate inputs. Apply T7.
For N1=3 the budget is exact, and the last two E1 steps instead have
six objects and no rows; use Q. A penultimate E2 was already excluded.
There is no larger N1 within the budget. This exhausts q_T=0.

For q_T=1, E<=d-2. If E<=d-3, (11) gives R>=10. Otherwise

    s_T+2N0+3NB+N1+K+Lambda <=1.

Thus N0=NB=0. If N1=0, consider the last E2. Afterwards a<=d+1.
If a=d it is the last paid step: immediately before it q=3 and there
are seven active coordinates, three pairs and one free. They are all
essential, so at least six objects are needed, but the four E2 deletions
and remaining budget allow at most five. If a=d+1 then q<=1 afterwards.
For q=1 the preceding three-pair/two-free state is excluded by Lemma R,
even if a paired cleanup subsequently removes that pair. For q=0 the
last paid step is D, so the preceding E2 state has two pairs, three free
coordinates and at most seven objects. At most one row can occur beyond
the six guaranteed acyclic deletions; replace it by XOR and apply T7.

If N1=1, every other budget term is zero and there are d-3 E2 events.
The last inert event is this E1. If it ends at a=d, its preceding q=2
state has six active coordinates and only three gates. All four paired
coordinates are used; both free coordinates are essential, since an
unused free direction plus the plane forced by an18-point two-pair/one-
free domain would be a three-flat. Six input uses and two non-output
uses exceed six pins. If the E1 ends at a=d+1, its after-q is at most
one. After-q=0 is followed by D and gives the forbidden five-gate history
state immediately before E1. After-q=1 must lose this pair by one A/C
cleanup before D, else the final D would end with q=2. The state before
E1 therefore has seven active coordinates, two pairs and three free,
and seven objects: three E1 deletions, two A/C and two D. At most two
rows occur. Replace them by two counted XORs and use T7. These exhaust
all q_T=1 possibilities, including R<7.

For q_T=2, E<=d-2. If E<=d-3, (11) already gives R>=8. If E=d-2,
then s_T>=3. Indeed two affine rows give affine output. With one row
and one gate, the solved row depends on at most two parameters and the
gate on at most three. Two gates and no rows likewise depend on at most
three of the four paired coordinates. An unused paired coordinate,
with its partner fixed zero, supplies an unused direction; on the other
pair choose a constant line. The active constant plane with E inert
directions contradicts dispersal. Zero or one object is even simpler.
Thus R<=7 forces s_T=3 and every other budget term zero. All inert
events are E2.

At the last E2, if the after-a is d+1 then its after-q must be one so
that the final D reaches terminal q=2. The preceding q=3 state is excluded
by Lemma R. If the after-a is d, immediately before it q=4,s=7 and all
eight active coordinates are paired. There are at most three rows and
no parameter inputs in them: paired row cleanup is exhausted and all
unpaired parameters are inert. An invertible system of at most three
such rows has a short row. For three auxiliaries, if all rows had three
terms they would all be111 and singular; the smaller cases are immediate.
Short-row deletion was exhausted, so there are no rows. Seven acyclic
gates, eight used inputs and six non-output uses exhaust fourteen pins.
Apply T8 to this read-once state, giving a forbidden constant three-flat.

For q_T=3, E<=d-3. If E<d-3, (11) gives at least10. Otherwise R<=7
allows at most one object. All unpaired coordinates are inert. A constant,
literal, one gate or one affine row depends on at most two paired
coordinates. If those coordinates are in different pairs, set each to
zero along legal lines; if they are partners, choose their constant
pair line by pigeonhole. Complete with any line in the remaining pairs.
This gives a constant three-flat and contradicts dispersal.

Finally q_T>=4 gives R>=2q_T+s_T>=8 from (11) and E+q_T<=d.
There is no terminal state with R<=7. Integrality proves (T), and the
five-input example in Appendix B.3 proves its stated additive sharpness.

## 7. A degree-sensitive bound, reproduced from the declared source

**Theorem.** Let 4<=B<=d. If f has algebraic degree at least B on every
dimension-d affine subspace, then

    C_B2(f)>=4n-8d+4B+3.                                (15)

This argument is supplied here in full. Its material antecedent is the
declared8eef source; no degree theorem or prior verdict is assumed.
Algebraic degree means degree of the unique multilinear polynomial over
F_2. Affine substitution cannot increase it.

For n<d, or n=d and B<d, the right side is negative. For n=d=B it is3,
and the elementary degree bound below gives C>=B-1>=3. Hence assume n>d
and run the same charged reduction to terminal a=d. Its represented
output is defined on all independent parameter assignments, not only
the promise. Put r=d-E and a0=r-q, where q is the terminal pair count.
Every legal product of pair lines, with all unpaired parameters retained,
is an affine d-space. Its output ignores the E inert directions and
has degree at least B. Thus r>=B.

The terminal object count satisfies s>=B-1. Solve the charged affine
system only for this degree upper bound, granting its affine values for
free. Each non-affine B2 gate then needs one multiplication, with free
affine shifts. A Boolean computation with j multiplications and arbitrary
free affine operations has degree at most j+1: for j=0 it is affine;
replace the first quadratic multiplication output by a new binary input
z and apply induction to the remaining j-1 multiplications. The resulting
ANF has degree at most j and is multilinear in z. Substituting the
quadratic expression raises degree by at most one. Since j<=s, the
asserted lower bound on s follows. This grants affine computation only
to upper-bound degree of the already charged representation; it removes
no object from the lower-bound ledger.

Let w count globally essential unpaired parameters and v pairs with both
members globally essential. Then w<=a0. In each unilateral-essential
pair set its essential member to zero, leaving its ineffective partner
free. In each ineffective pair fix either member; in each both-essential
pair choose any legal line. This remains a full d-space, although the
output depends on at most w+v effective variables. Hence w+v>=B.
If v>0 equality is impossible. Fix all other pair lines. For each of
one both-essential pair's three lines, the restriction depends on at
most B variables and must have full degree B. Its top coefficient,
the XOR of its truth table, must be one. XOR those three coefficients:
every legal point of this pair occurs twice for every assignment of the
other B-1 variables, so their sum is zero, a contradiction. Therefore

    v>0 implies w+v>=B+1.

The represented computation has at least w+2v essential parameters.
The fully charged connectivity proof in Section1 gives s>=w+2v-1.
Dropping nonnegative event terms from (11), its offset above4n-8d is
at least2r+2a0+s. If q=0, r>=B and s>=B-1 give at least5B-1>=4B+3.
If q>0,v>0, then r>=w+v>=B+1 and
s>=2B+1-a0, so the offset is at least4B+3+a0>=4B+3.
If q>0,v=0, then w>=B, a0>=B, r>=B+1 and s>=B-1 give at least
5B+1>=4B+3. All possibilities prove (15), including cyclic fair rows,
sharing, constants and ineffective variables.
