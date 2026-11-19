# A full-B2 eventual lower bound with relaxed terminal stars

For the single language defined in `family.md`, the theorem is

    C_B2(L_N) >= (241360436914016647911933609899531/
                  50000000000000000000000000000000) N
                =4.82720873828033295823867219799062 N

at EVERY sufficiently large integer actual bit length N. The explicit
eventual remainder is zero; its starting threshold is unspecified.
The model permits all sixteen two-input Boolean truth tables, unrestricted
sharing and depth, repeated input pins and constant/unary/terminal cases.
`family.md` supplies one fixed polynomial-time evaluator, hence an NP
verifier with empty witness, at every length. No advice or certificate
bits are given to the circuit.

The argument is self-contained in the seven listed scientific files.
`graph.md` proves Gaussian comparison, the conversion from crossing
matchings to cutwidth, a finite J47 path bound and the physical multigraph
reduction. `curvature.md` pays three disjoint pair sets in one covariance
budget. `source.md` pays the actual-base ordering probabilities and two
different threshold deficits. This file proves one simultaneous ordering
and population/energy comparison. `family.md` proves its circuit and
norm-graph consumers. `certify.py` replays the finite inequalities; code
execution does not replace any of the general arguments.

Scientific provenance: the norm-one/circuit route builds on
McGuire--Schlesinger--Williams, the full B2 normalization of Li--Yang,
and the norm-graph algebra of Kollár--Rónyai--Szabó, with the exact forms
proved in `family.md`. The Gaussian matching-cut/cover method is due to
the public CircuitBreaker lineage including galarnon42, Arthur Ostapenko,
mooselumph and DarkWingMcQuack. The positive radial and random-strip
construction is recorded at Puthranps/cubic-cutwidth-retiming commits
`a901153f78d426225407c02dec82995e51192f6d` and
`6b299fde41307ff47eb4a699e140fd009e61a51d`. The one-pass weak exposure
principle also appears in public CircuitBreaker proof
`2ce5e00ab57185fbf7bfdb504b316eb8738e9fda`. All mathematical premises used
here are expressed and proved in these supplied files. Public availability,
past verdicts and numerical screens are not proof premises. The text and
standalone implementation in this release are independently authored.

## Fixed inputs and objective

For a simple cubic base on n vertices let F(B) be the full maximum
matching size of the crossing graph of B, equivalently its minimum
vertex cover number. The original keys Y are the actual finite radial
source from `graph.md`; excess E is the sum of its complete-alpha
adjacent excesses, used only for marking and covariance accounting.
The relevant exact inputs are

    U=26129513284633533/200000000000000000,
    eta=19/40000, theta=eta/2,
    K=28238380378679347/10^18,
    p=247011/10^12, u=301/10^10, w=4695/10^12,
    q_inner=283/10^9, q_annulus=273249/10^12,
    b=451347/20000000000,
    delta=4052/10^9.

Here U bounds the reference zero-threshold matching expectation per
vertex; K E is the original-base zero-threshold covariance gain.
The p,q,u,w are probabilities on the ORIGINAL actual base field, rather
than on a retimed or independent reference field.

## One safe ordering, with a disjoint neutral-group saving

Use every seven- and nine-key positive support defined in `source.md`,
including the relaxed-ending variant proved there. Seven-key lightness
uses eta; original nine-key lightness uses theta throughout. For the
relaxed variant the four terminal child edges may use eta; all other
edges retain theta. Each strict
support gives a one-unit decrease when its outside root is added to B.
Also allow a weak support: a neighbor a of a root r together with both
of a's other neighbors. If all three are in B, a is an isolated vertex
of the crossing graph after r is deleted. The original crossing edge ra
therefore adds one to the matching optimum of that graph; after flipping
r at most one can be added. This proves the root marginal is nonpositive,
with no condition on the exterior.

For every root v replace its key by the minimum of its original key and
the maxima of all eligible strict and weak supports. This is a single
radius-three function of the original keys. At any threshold every newly
added root has its support entirely in the ORIGINAL B, so no other added
root lies in that support. Submodularity preserves the nonpositive
marginal of every weak proposal and the -1 marginal of every strict
proposal as other roots are added. Thus the whole update is safe.

There is an extra gain from original zero marginals. Let m and all its
neighbors be originally outside B. Qualify a neighbor r if its two
other neighbors are inside B and the lower-key one, a, has both its
other neighbors inside B. This chooses a canonical private crossing
leaf a. Write R_m for the qualified roots, k_m=|R_m|<=3. Each r has
unique original outside neighbor m, so the moving sets R_m are disjoint
between centers. Their private leaves are distinct as well: a crossing
leaf cannot attach to a second outside root.

Remove the k_m private leaf edges by their constant perfect matchings.
Originally m is isolated. After flipping R_m its k_m roots become
leaves at m, adding one matching edge if k_m>0. Consequently

    F(B)-F(B+R_m)=(k_m-1)_+,

and every individual member had original marginal zero. A fixed disjoint
addition set's total marginal decreases as the initial B grows, by
summing its singleton marginals. Therefore these group gains persist
after earlier disjoint groups. Strict roots cannot belong to any R_m
because their original marginal is -1, not zero. Process the groups,
then all strict roots, then the remaining safe proposals to obtain

    F(B)-F(updated B)
      >=number of original strict roots+sum_m(k_m-1)_+.

Choose an unordered pair r,s among m's neighbors and one private candidate
at each root: twelve canonical names. Each successful name selects its
private candidate as the lower-key child. A center with k_m qualified
roots has exactly binomial(k_m,2) successful pair names. Weight each by
2/3, since `(2/3)binomial(k,2)<=(k-1)_+` for k=0,1,2,3. The nominal
weight per center is eight. Its eleven-edge light event is precisely
the actual-source cone of (S10); the negative guards enforce original
neutrality and disjointness from the strict flags.

Let L7,L9,LR,LW be the eligible naming weights, with LR the relaxed
nine-key count and LR>=L9. At |t|<=.03, the expected
ordering saving is therefore at least

    q L7+(p-u)L9+u LR+w LW,                                 (P1)

with q=q_inner on |t|<=.02 and q=q_annulus otherwise. The same actual
keys are used by every part of the update, and no independence of
successful events is assumed.

## Naming counts from inspected stars

An edge is high if D>eta, medium if theta<D<=eta, and low otherwise.
Let n_ab count vertices incident to a high and b medium edges, a+b<=3.
Let s_H count non-high edges with two high-marked endpoints, and s_S
low edges with two endpoints marked by high or medium edges. At one
cutoff, write N for its number of marked vertices and r for its marked-
unmarked boundary. Then

    r=sum_(d>0)(3-d)n_d-2s.

For an inspected naming, either its root is marked or the first marked
star on a branch follows an unmarked predecessor. In a named tree of
nominal weight W, any one directed predecessor position hits every
directed base edge with total weight W/3: fix that directed image and
extend the two remaining parts by locally injective choices. The number
depends only on the tree. Remote projection collisions are allowed
in this upper count; lightness supplies faithfulness afterwards.

The seven-key inspection has three nonroot complete stars and W=6.
Its loss is at most 6N+6r, giving

    L7>=6n-sum c7(a)n_ab+12s_H,
    c7=(0,18,12,6).                                 (P2)

The weak-pair inspection has four nonroot complete stars and W=8.
At the lower cutoff it gives

    LW>=8n-sum cW(a+b)n_ab+(64/3)s_S,
    cW=(0,88/3,56/3,8).                             (P3)

Negative right sides are harmless; the true counts are nonnegative.

A stronger nine-key count uses the unmarked induced graph. Let U_i be
the number of unmarked vertices with i marked neighbors, U=sum U_i,
and J the number of U_1--U_1 edges. At a U_0 root all three first
neighbors are unmarked. If j are of type U_1 and k of type U_2, their
available child counts are 1,0,2 respectively. The eligible count is
the sum of the three pairwise products of those child counts. For k=0
the values at j=0,1,2,3 are 12,8,5,3. In particular they dominate both
`9-2j-17k` and `11-3j-23k`. If k>0 their simultaneous surpluses are:

| (j,k) | eligible | first surplus | second surplus |
|---|---:|---:|---:|
| (0,1) | 4 | 12 | 16 |
| (1,1) | 2 | 12 | 17 |
| (2,1) | 1 | 13 | 19 |
| (0,2) | 0 | 25 | 35 |
| (1,2) | 0 | 27 | 38 |
| (0,3) | 0 | 42 | 58 |

The extra is at least 12k and 16k respectively. Write E01,E02 for
U_0--U_1 and U_0--U_2 edges. There are at most two U-stubs at U_1
and one at U_2, with `2U_1-E01>=2J`. Retaining those missing-stub
terms together with the table, rather than bounding them separately,
proves the simultaneous inequalities

    L9>=9U-13r+12U_2+30U_3+4J,
    L9>=11U-17r+16U_2+40U_3+6J.                    (P4)

For example the first residual is
`2(2U_1-E01)+17(U_2-E02)+12E02>=4J+12U_2`; the second uses
`3(2U_1-E01)+23(U_2-E02)+16E02>=6J+16U_2`.
Since `N+r=sum_(a+b>0)(4-a-b)n_ab-2s_S`, the main parts of (P4)
are equivalently

    9n-sum w9(a+b)n_ab+26s_S,
    11n-sum w11(a+b)n_ab+34s_S,
    w9=(0,35,22,9), w11=(0,45,28,11).

Both inequalities refer to the SAME naming count, so they are combined
only with nonnegative weights whose sum is p-u.

The relaxed count additionally satisfies

    LR>=12n-sum cR(a,b)n_ab+24s_S+12U_2+24U_3,
    cR(a,b)=(0,36,24,12)[a+b]+8(3-a-b)1[a>0].       (P4R)

At a U_0 root all four core stars are unmarked, so all twelve names
exist before terminal testing. A candidate ending vertex incident to
a high edge deletes at most four names: choose the other selected
branch and one of its two children. A high vertex of type (a,b) has
at most 3-a-b unmarked neighbors, and each intermediate neighbor
belongs to at most two U_0 roots. Thus at most
`8 sum_(a>0)(3-a-b)n_ab` names are lost from 12U_0. Using
`U_0=n-|S|-r+U_2+2U_3` and the marked-boundary identity gives (P4R).
Repeated rejection reasons only weaken this incidence upper bound.

## Unmarked excess pays a capped boundary improvement

Normalize x_e=D_e/eta and put A=eta K. In this section all finite-alpha
root debits from (S3) are collected as an arbitrarily small multiple
of n. This means source resources are chosen first to achieve that
accuracy; a fixed nonzero root debit is not said to decay with n.

At a U_1 vertex let x be its low boundary excess and y,z its two U--U
edge excesses. The first constraint (S3) gives `y+z>=x/2` up to the
uniform debit. Put `Bcap=sum_(S--U_1 edges)min(x,1/4)` and let L_UU
be the total U--U excess. Counting at most two visits per edge gives
`L_UU>=Bcap/4`. Alternatively cap each visited y,z at 1/8. Their sum
is at least min(x,1/4)/2. A twice visited edge contributes at most
`2min(x_e,1/8)<=x_e+1/8`; these are exactly the J edges. Thus
`L_UU>=Bcap/2-J/8`. Averaging these two bounds for the SAME energy gives

    L_UU>=beta Bcap-J/16, beta=3/8.                 (P5)

Choose nonnegative v9,v11 with v9+v11=p-u and set

    z_s=13v9+17v11+12u+(32/3)w.

The shared-edge terms in (P3)-(P4R) supply exactly `2z_s s_S`.
We now give one explicit allocation of `K E+2z_s s_S`. Keep the high
and medium baselines `A sum(a/2+b/4)n_ab`. Allocate half the medium
surplus above 1/2 to each endpoint of type 02 or11, at most two such
payments per edge. At eligible low endpoints of type 10,01,02, a
shared S--S low edge pays `A x/2+z_s`; a boundary S--U edge pays
`A x+beta A min(x,1/4)`. Type11 receives just A x/2.
Other low endpoint payments may be discarded.

All the extra boundary payments landing at U_1 fit (P5), except its
debit A J/16. At U_2 and U_3 there are two and three boundary edges,
each capped at beta A/4. Thus these formal payments fit the ONE amount

    K E+2z_s s_S+(3A/16)U_2+(9A/32)U_3+(A/16)J
                         +arbitrarily_small_error*n. (P6)

No high, medium or shared-low energy is recharged by the U--U allocation.
For total eligible low incidence l, the formal endpoint payment is at least

    F_beta(l)=A l/2+min(A l/2+beta A min(l,1/4),z_s). (P7)

To verify it, first sum the half-energy baseline. The remaining edge
payments are bounded below by minima of z_s and
`A x/2+beta A min(x,1/4)`. The sum of these minima is at least the
minimum of the sum and z_s, and clipped incidences satisfy
`sum min(x_i,1/4)>=min(sum x_i,1/4)`.

The exact constants below satisfy
`(7/32)A<z_s<(11/32)A`. Therefore

    F_beta(1/4)=(1+beta)A/4,
    F_beta(1/2)=A/4+z_s,
    F_beta(11/20)=11A/40+z_s.

The function is increasing and concave with slope at least A/2
everywhere. Also `F_beta(l)>=A min(l,1/20)`: on [0,1/20] its first
branch has slope (1+beta)A and the cap exceeds the amount needed;
beyond it use monotonicity.

## Paired endpoint payments, with shared endpoints included

By (S3), a high edge joining types 10 or11 has off-edge incidence at
least 8/5. Their baseline incidences total one, and each extra low or
medium-surplus unit receives at least A/2. Such a pair therefore pays
bonus 3A/10. These high pairs are vertex-disjoint. Their number is at
least half of n10+n11 minus half of the high stubs at all other types.

A medium 01--01 pair has low incidence at least 1/4 at each end and
sum at least 4/5. Monotonicity and concavity of (P7) put the smallest
sum at `(1/4,11/20)`. Above the two 01 local baselines its bonus is

    Delta=z_s+c0-beta A/4, c0=A/40.

A medium 01--02 edge also pays c0, with no hidden duplication when
two such edges share a 02 endpoint. Write y=l01-1/4, l for the 02
low incidence and s_e,s_f for its two medium surpluses above 1/2.
The paired inequality gives

    y+l+s_f>=1/20+(8/5)s_e.

With two 01 neighbors sum the two inequalities to get
`y1+y2+2l>=1/10+(3/5)(s_e+s_f)`. The available low payments are (P7) at 02, at least `A min(l,1/20)`,
and the 01 increments, at least `Ay/2` each, plus `A(s_e+s_f)/2`.
If l>=1/20, the 02 payment alone supplies A/20, enough for two edges.
Otherwise one neighbor is paid by `Ay/2+Al+A(s_e+s_f)/2`, at least
half its constraint. For two neighbors, `A(y1+y2)/2+Al` is at least
half their summed constraint. Thus the same A/40 per edge is paid
within the SAME allocation, without assuming slope A at 01.
The type01 leaves of these mixed pairs cannot belong to a 01--01 pair,
and no medium-pair endpoint belongs to the high-pair type set {10,11}.

Let P,Q count medium 01--01 and 01--02 edges, and let R count medium
stubs at other types. Then `2P+Q>=n01-R` and
`P>=(n01-2n02-R)/2`. Since Delta>=2c0, for xi=2/3 the combined bonus
is at least

    c0(n01-R)+(g/2)(n01-2n02-R),
    g=(2/3)(Delta-2c0), B=c0+g/2.

This uses `(Delta-2c0)P` once, interpolating its signed count with
P>=0. Together with the high pairs the local formal-payment coefficients
are exactly:

| (a,b) | G_ab |
|---|---|
| 00 | 0 |
| 01 | A(1/2+beta/4)+B |
| 02 | A/2-g |
| 03 | 3A/4-3B |
| 10 | 9A/10+z_s |
| 11 | 9A/10-B |
| 12 | 17A/20-2B |
| 20 | 7A/10 |
| 21 | 19A/20-B |
| 30 | 21A/20 |

This table is a lower bound for the same sum of formal endpoint payments
in (P6), rather than a second independent energy estimate.

## Exact dual and nested thresholds

At the fixed inputs choose

    v11=((31/20+beta/2)A-92p-4u-(232/3)w)/26,
    v9=p-u-v11.

Both weights are positive. The exact arithmetic also checks the stated
range of z_s, g>0 and all ten type comparisons below. For the inner
range |t|<=.02 put q=q_inner,d=0. For the annulus
`.02<=|t|<=.03` put q=q_annulus and `d=b(.02)^2`.
By (S13) the extra actual-base threshold credit is at least
`d(6n-28h7)`. Replacing the positive part by this weaker expression
permits linear accounting. Since `h7=sum a n_ab/2`, its type debit is
14da. The comparisons are

    G_ab >= q c7(a)+v9 w9(a+b)+v11 w11(a+b)
                                      +u cR(a,b)+w cW(a+b)+14da. (P8)

They bind at 00,01,02 and have strictly positive slack at the other
seven types. The simultaneous naming surpluses in (P4) also exceed
all three overlap costs of (P6):

    12v9+16v11+12u>3A/16,
    30v9+40v11+24u>9A/32,
    4v9+6v11>A/16.                                  (P9)

Sum (P8)-(P9), (P2)-(P4) and (P6). The shared-edge allowance cancels
exactly once; the favorable 12q s_H may be dropped. Equations (P1)
and (S13), telescoped after the original zero-threshold Price comparison,
give total central saving per vertex at least

    6q+9v9+11v11+12u+8w+6d.

Each exact range gains `(35/13)u` over the respective values
`165375882580092155427/41600000000000000000000000` and
`165195157204092155427/41600000000000000000000000`.
The annular floor exceeds `4052/10^9`, with strict gap greater
than `7/10^11`.
All these are direct rational comparisons reproduced by `certify.py`.

For |t|>=.03, the same-threshold cover estimate (S11) pays
`67/10^7>delta`. The combined rule is pointwise safe. Choose the
finitely many source, conditioning, cover and mesh errors strictly
inside the smallest available gap, uniformly before graph size.
The radius-three inheritance, finite-range concentration, tail and
tie treatment and (G4) in `graph.md` now imply the uniform simple-cubic
bound with

    gamma=U-delta=26128702884633533/200000000000000000.

Its physical reduction supplies one nonnegative nondecreasing h=o(D)
for all cubic multigraphs with the same gamma.

## Actual input length and final rounding

Apply the circuit/norm-graph theorem proved in `family.md` to this physical
graph hypothesis. It uses the same one total validated-prefix norm-one
language, charges every prefix and leftover bit, and gives at all large
integer lengths a uniform sublinear loss below
`c0=1+1/(2gamma)`. Its fixed evaluator supplies the empty-witness NP
verifier; constants, degenerate truth tables, shared gates, repeated pins
and terminal outputs are included in its circuit normalization.

The declared c* is the strict downward 10^-32 rounding of c0. In particular

    0<c0-c*<10^-32.

The same sublinear loss divided by actual N tends to zero along EVERY
integer length, so some fixed N0 absorbs that positive rounding gap
for all N>=N0. This gives the exact declaration with identically zero
eventual remainder, without asserting an explicit value of N0.
