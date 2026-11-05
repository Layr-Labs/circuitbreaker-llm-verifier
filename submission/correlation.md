# Affine correlation breaking by indexed alternating extraction

The starting organization is Xin Li's alternating
extraction construction, ECCC TR15-121, Sections 3--4. We reproduce the
affine conditioning and indexed row argument with different, explicitly
charged linear extractor thresholds. No correlation-breaker theorem is
used as an external premise.

All ranks and entropies below are binary. An affine source is uniform
on a nonempty affine subspace. Its entropy is the dimension of its
direction space. Joint distributions and total variation errors, rather
than assertions about a typical sample, are the objects being compared.

## 1. Fixed parameters and internal extractor contracts

Let N be the entire input length. Put

    L=floor(log_2 N), J=ceil(log_2(L+2)),
    k=floor(2^600*L^11/J), m=floor(k/(2^100*L^2)),
    h=the least power of two >=2^90*L^2,
    ell=2^110*L^2, C_L=1,
    Qlen=64h^2*ell,
    d=2^40*(h^2*ell+L),
    E=2^90*L^2, E_1=2^20*h^2*ell+E,
    epsilon=2^-E, epsilon_1=2^-E_1.

All floors, ceilings and powers are exact integer operations. The initial short linear affine
extractor Short from seeded.md has output length m and seed length s.
Its error is at most 2N^-8 and s<=2^20*log_2N eventually. Write

    M=2^s, l=log_2 h, b=ceil(s/l).

The case l=0 is an exceptional guarded shape. An index i in [M] is
represented by the s-bit binary expression of i-1. Divide it from most
significant end into b blocks of length l, padding the final block on
the right with zeros. Let a_(i,j) in [h] be one plus the integer encoded
by block j. Lexicographic order of all b blocks agrees with index order.

Every internal extractor is the complete bounded-order diagonal
construction in seeded.md. It computes its exact field width w and
entropy threshold H_call=ceil((2^20*output*w)^(33/32)), and guards
capacity, reconstruction advice, ambient entropy and supplied seed.
Larger seeds ignore their suffix, preserving strong extraction.
The calls and their exact required entropies are:

| name | ambient | output | supplied seed | error | entropy threshold |
|---|---:|---:|---:|---:|---:|
| Xrefresh | N | ell | d | epsilon_1 | H_Xr |
| Yrefresh | m | Qlen | ell | epsilon | H_Yr |
| Ynext | m | d | ell | epsilon | H_Yn |
| Xalt | N | ell | ell | epsilon | H_Xa |
| Qalt | Qlen | ell | ell | epsilon | H_Qa |

Put H_alt=max(H_Xa,H_Qa), H_X=max(H_Xr,H_Xa),
H_Y=max(H_Yr,H_Yn). These are the ACTUAL coded thresholds,
not a claim of constant-fraction extraction at entropy O(output).
In particular the larger error exponent in Xrefresh is charged
in its actual H_Xr. These values are fixed numerical metadata at N.

Let

    B_0=2^20*C_L*b*h*(d+Qlen+h*ell).

The numerical guards include

    m>=4B_0+4H_Y, k>=4h*m+4B_0+4H_X,
    h^2>=C_L, d>=64h^2*ell, m>=4d,
    Qlen>=8h^2*ell+8H_alt,
    k-h*m-B_0>=4h^2*ell+4H_alt,

and every actual entropy/seed/capacity inequality in the table. They
are deterministic conditions independent of the source direction.

In powers of L, d has order L^6, Qlen has order L^6, h has order L^2,
and b has order L/log L. Thus B_0 has order L^9/log L, as does m.
The main proof gives B_0<=2^466*L^9/log_2 L and
m>=2^497*L^9/log_2 L eventually, paying this borderline comparison.
Also h*m<=2^-9*k. The actual entropy thresholds have growth powers

    H_Xr,H_Yn,H_Yr=Theta(L^(33/4)),
    H_Xa,H_Qa=Theta(L^(33/8)).

All lie strictly below their available source budgets. Ordinary seeds
are below 2^17*E, which ell pays with a fixed factor; Xrefresh's
seed is below 2^38*h^2*ell<d/4. The final resilience order/error
charges are at most 2^85*L^2 and 2^87*L^2, below h and E.
The main proof supplies the complete numerical audit.

## 2. Affine facts and exact good-seed accounting

If X is affine and Z=L(X) is a linear observation, choose a complement
to the kernel of L inside the direction space of X. Then

    X=A+B,

where A and B are independent affine sources, A determines Z bijectively,
and B is uniform on the kernel directions. Thus H(B)=H(X|Z), which is
a constant rank independent of the observed value. Affine translations
can be allocated to either term. This decomposition is used only in
the proof and supplies no algorithmic advice.

Conditioning on a fixed value of any linear observation gives either an
empty fiber or an affine source. The rank of another linear observation
can decrease by at most the number of bits fixed. This includes
conditional ranks: H(U|Z) decreases by at most r upon fixing r further
linear bits. If the new observation is a function of Z, it does not
reduce H(X|Z).

A fixed-seed linear extractor applied to an affine source produces the
uniform distribution on an affine image. Such an image is either fully
uniform or at distance at least 1/2 from the full output cube. Therefore
a strong linear extractor of error eta, at a sufficient entropy budget,
has at most 2eta bad uniform seeds for any fixed source direction.
All other seeds give EXACT uniform output. Affine translations do not
change the good set.

We use this rank fact with side information. Suppose B is an affine
residual independent of an affine variable A, S is a function of A,
and the extractor on B has a good seed S. Its output is uniform even
given every variable determined by A. Adding an A-dependent affine
translation to that output does not change this conclusion. Conversely,
if a uniform seed is a function of B and an extractor acts linearly on
an independent affine A, its good-seed output is uniform even together
with arbitrary B-side information containing that seed. The source
DIRECTION after previous A observations must be fixed independently of
the newly sampled B seed. The order of observations below ensures this;
we never condition on a later seed-dependent A observation for free.

For a deficient affine seed of s bits and entropy at least s-r,
independent of the source residual, the bad-seed probability is at most
2^(r+1)*eta. This follows by bounding its point probabilities by 2^r
 times those of a uniform s-bit seed and applying the same direction-rank
bad-set bound. Thus Xrefresh can tolerate a defect of at most

    r_def=64h^2*ell.

Section 5.2 proves this smaller defect by revealing current smaller-row
seeds and seed-exposed alternating prefixes, NOT their entire current
Q words. Those Q words are recorded only AFTER the fresh target selected
response, before the next target Ynext direction test (Section 5.4).
The timing matters: full current Q words would give the older h*Qlen
defect, which is not used here. Since E_1=2^20*h^2*ell+E,
2^(r_def+1)*epsilon_1<=epsilon eventually, with fixed slack.
The guard d>=r_def makes the seed-entropy statement legal.

All estimates can be expressed by discarding bad-seed fibers of their
bounded total probability. Before discarding, fix the relevant seed
values; the remaining observations are linear, so each retained
component is affine. The result is a mixture of affine fibers, not a
claim that a nonlinear conditioned event itself is affine. Keeping
transcripts in the joint distribution and later projecting them away
cannot increase variation distance. This is the meaning of the error
accounting in the next sections.

## 3. Alternating extraction with linearly correlated Q sources

Given X and a Qlen-bit word Q, define S_1 to be the first ell bits of Q.
For t=1,...,h perform

    R_t=Xalt(X,S_t),
    S_(t+1)=Qalt(Q,R_t) if t<h.

Write LA(X,Q)=(R_1,...,R_h). Its evaluation is deterministic. Different
rows may use arbitrarily correlated Q words; no independence between
their seeds is assumed.

Here is the precise fact used. Suppose X is affine, Z is a linear
function of it, and Q_1,...,Q_u (u<=h) are linear functions of Z. There
is a distinguished Q_v with uniform first ell bits, and with entropy
at least Qlen minus already charged observations. Assume its remaining
entropy is at least 4h^2*ell+4H_alt, and H(X|Z)>=4h^2*ell+4H_alt.
We may also track an affine U, linear in Z, with enough conditional
entropy given the Q words for later extraction. After at most h rounds,
the distinguished R_t is within a coarse 8h*epsilon error of uniform
even jointly with all the Q words, all current S words and all earlier
R words. The next distinguished S_(t+1) is similarly uniform jointly
with the current R words and previous S words; the current Q words are
NOT included in that second uniformity assertion.

For clarity the accounting proof is given. On a current affine fiber
write X=A+B using Z as in Section 2. Current Q words and seeds are
functions of A; B is independent of all of them. At an X step, a uniform
distinguished seed has bad-direction probability at most 2epsilon.
Fix all current seed values, then fix the A components of the u responses.
These observations affect A only. For a good distinguished seed its B
response is uniform, so adding its fixed A response leaves R_v uniform
even jointly with every still-random A variable. All other responses
may be correlated with R_v; they are not included among the earlier
responses in this uniformity assertion.

Fix the B components of all u responses next. The new Q step has seeds
that are functions of B, with their A components fixed. The distinguished
seed R_v is uniform before those B values are recorded and independent
of the remaining A source. At most another 2epsilon of its seed values
fail the Qalt direction-rank test. For every other value, its Q output
is exactly uniform even together with the other B responses. Each new
S_i is a linear function of the corresponding Q_i once its R_i seed
is fixed. This restores the same affine form for the next step.

The word "uniform" here always denotes a joint statement on the
retained good-seed components. Recording arbitrary other B-side
information after the rank test does not make the output nonuniform:
for each good seed the output of the independent A source is exactly
uniform, regardless of the other B values. Nor is there an assertion
that R_v remains uniform after revealing a simultaneously equal R_i.
The next S uniformity is obtained by fresh extraction from A instead.

Per step, fixing the current u seeds and the A response components
costs at most 2h*ell of Q_v entropy. Fixing B components costs at most
h*ell of residual X entropy. The entropy of a tracked U conditional on
all Q words loses at most h*ell, since seed observations are functions
of the Q words. After h rounds these losses are at most 2h^2*ell,
h^2*ell and h^2*ell respectively. They fit the stated hypotheses and
leave H_alt for every Xalt and Qalt call. Accumulating the two rank
exceptions per step gives at most 4h*epsilon; the looser 8h*epsilon
covers initial component errors and endpoints.

## 4. Exact indexed row algorithm

For the source input X, enumerate all s-bit seeds of Short in increasing
binary order and put

    Y_i=Short(X,i-1), 1<=i<=M.

Every Y_i has length m and is a linear function of X. For each row i,
start P_(i,1) with the first d bits of Y_i. For j=1,...,b do

    T_(i,j)=Xrefresh(X,P_(i,j)),
    Q_(i,j)=Yrefresh(Y_i,T_(i,j)),
    (R_(i,j,1),...,R_(i,j,h))=LA(X,Q_(i,j)),
    O_(i,j)=R_(i,j,a_(i,j)).

If j<b set

    P_(i,j+1)=Ynext(Y_i,O_(i,j)).

Finally output W_i=O_(i,b). Every field/code primitive and truncation
was specified in seeded.md. There are polynomially many rows and a
polynomial number of operations per row, since s=O(log N), b<=s and
all remaining lengths and loop counts are polynomial in N. All calls
are on projections or linear observations of the SAME N-bit input;
there is no independent external randomness in the final evaluator.

## 5. Lexicographic separation lemma

Fix any set T of at most h row indices and a distinguished v in T for
which Y_v is exactly uniform. For this independence assertion only
indices at most v matter. Replace the working tuple by {i in T:i<=v}
and rename that prefix T, so v=max T. Every reference below to ALL
current Q words means precisely this target and its smaller rows,
whose current T seeds have been processed. Larger indices are excluded
from the working tuple and from the conditioning assertions. The final
chain rule applies this proof separately to each sorted initial prefix
of an arbitrary at-most-h good-row set. Other Y_i in T can be arbitrary affine
images and may even equal Y_v. Let Z be their concatenation, including
Y_v. Then H(X|Z)>=k-h*m. At any intermediate affine fiber decompose
X=A+B with A determined by Z and B its independent residual.

For a prefix length j, split the smaller indices i<v into those whose
first j index blocks are strictly less than v's and those whose first
j blocks equal v's. Denote the former set by L_j. Thus L_0 is empty and
L_b contains every smaller index in T.

The proof transcript is refined by recording the distinguished row's
OWN previous refresh seeds as well. At the beginning of stage j it
contains its earlier P, T and Q words, its seed-exposed alternating
prefixes through the previously SELECTED slots, and its previous selected
responses O_(v,r), r<j. The last of these are the fixed Ynext seeds.
It does not yet contain the CURRENT P_(v,j). Consequently that current
P is a linear function of the original Y_v on each transcript fiber;
it is not treated as a linear function merely because it was produced
by a nonlinear seed-dependent computation. Unused responses AFTER the
target's selected slot are never observed. Their computation in the
evaluator does not alter any selected output.

For a smaller row i at stage r one can equivalently truncate its
alternating transcript after max(a_(v,r)-1,a_(i,r)) slots. Expose all
seeds in order before their responses and refine observations by A and
B response parts when invoked below. At an output-independence assertion
only rows whose advice prefix is already strictly smaller contribute
their simultaneous selected responses. A still equal-prefix response
at the target's slot is excluded from that assertion. Before refreshing
the NEXT target seed it may be recorded as B-side information, after
its A part was fixed. These are distinct assertions, not a promise of
independence from an equal simultaneous response.

We claim the following inductively, with a coarse total error at most
2^12*b*h^2*epsilon. At the beginning of stage j, P_(v,j) is exactly
uniform on the retained affine components even jointly with that complete
PAST proof transcript (including the distinguished row's own past
observations), all completed smaller-row proof transcripts, and the
current P_(i,j) for i in L_(j-1).
At the end of that stage, O_(v,j) is uniform jointly with the selected
outputs O_(i,j) for i in L_j, and with the seed/linear observations
needed to run the following stage. Other equal-prefix current outputs
are deliberately NOT asserted independent of it. After the Ynext
refresh the beginning assertion holds for stage j+1.

The proof uses the following explicit order; it covers both an earlier
prefix difference and a first difference in the current block.

### 5.1 Completed transcripts remain affine and their ranks are charged

Once a row seed P_(i,j) is fixed, T_(i,j) is linear in X. Once that T
is fixed, Q_(i,j) is linear in the original Y_i and hence in Z. In the
look-ahead chain, fix the S seed before its R response, and that R seed
before the next S. Every fixed observation is then linear on the
current affine fiber. Decompose responses into their A and B parts
when needed; each part is a linear observation of the corresponding
independent component. This proves affine preservation, rather than
assuming that a nonlinear transcript fiber is affine without exposing
its seeds.

For a smaller row in one stage, the current P has d bits, Q has Qlen
bits, and at most h+1 responses and h seeds of ell bits are needed.
The overall A-side rank losses over all stages are bounded by a fixed
multiple of b*h*(d+Qlen+h*ell), and B-side losses by a fixed multiple
of b*h^2*ell. The distinguished row's own retained transcript adds at most a fixed
multiple of d+Qlen+h*ell per stage, including its selected response and
all prefix seed/response parts. Since h>=2, this is absorbed by doubling
the h-row bound. Separate A/B observation parts only multiply it by a
fixed factor. For a concrete conservative count, each row/stage records at most
2d bits for current/next seeds, Qlen for its postponed Q word,
2ell for T parts, h seed words and 2h response parts of ell bits,
and at most four endpoint ell words. Both the distinguished row and
all smaller rows together have at most h members. Doubling for any
repeated split or retained-prefix bookkeeping gives a total below
64b*h*(d+Qlen+h*ell), for h>=2. This includes the target's own past
transcript and the delayed full-Q observation. The factor 2^20 in
B_0=2^20*b*h*(d+Qlen+h*ell), now with C_L=1, exceeds this bound;
both A-side and B-side losses are far below B_0. In particular retaining
the old target refresh seed is charged, not a free conditioning step.

Consequently the distinguished original Y_v retains at least m-B_0
entropy whenever it is used as a source, while the residual B retains
at least k-h*m-B_0. The strong numerical guards make these exceed all
required H_Y and H_X budgets. Conditional ranks of Y_v after
additionally observing h Q words lose at most h*Qlen and fit the same
coarse budget. The refreshed Q_v itself starts with Qlen entropy and
loses only the alternating-extraction costs of Section 3; the choice
Qlen=64h^2*ell leaves ample margin.

### 5.2 Reveal seed-exposed prefixes, postponing full Q words

Assume the beginning assertion for stage j. First fix the current seeds
P_(i,j) of indices in L_(j-1). The target P_(v,j) is uniform jointly
with them by induction, so this step does not reduce its entropy.
Process those already separated rows only far enough to obtain the
look-ahead prefix through max(a_(v,j)-1,a_(i,j)). The current Yrefresh
output Q_(i,j) is COMPUTED, but its whole value is NOT recorded yet.
No next-stage P_(i,j+1) is recorded before the current target refresh.

Here is the complete observation order and charge. With P_(i,j) fixed,
T_(i,j)=Xrefresh(X,P_(i,j)) is linear in X. Fix its A and B response
parts separately (at most 2ell bits). With T fixed, Q_i is a fixed
linear function of original Y_i, hence of Z and the A component.
Record its first ell bits S_(i,1), not all Qlen bits. Before each
R_(i,t), record its ell-bit seed S_(i,t); with that seed fixed the
response is linear in X, so record its A and B parts separately.
Before each later S_(i,t+1), its seed R_(i,t) is already fixed, so S
is a linear observation of Q_i and hence A. These are affine fibers,
and A/B residuals remain independent after their separate linear
observations. All matrices are fixed by already recorded seeds.

For at most h rows and h slots, a generous upper bound for new recorded
bits is 2h*ell+3h^2*ell, plus at most h*ell for endpoint seeds, below
8h^2*ell. Selected outputs are included among these responses. Copies
or repeated observations cannot increase the bound. The full Q words,
whose h*Qlen bits would be much larger, are deliberately absent.
The coarse r_def=64h^2*ell also covers split components and all endpoint
choices. Thus on each resulting affine component the target P is a
linear function of Z of entropy at least d-r_def, while its independent
B residual has the guarded entropy. This rank statement includes every
new observation, rather than treating a computed but unrecorded Q as
conditioned information. The matrices of the B observations were fixed
by recorded seeds before their data was exposed, so the B direction is
fixed before applying the remaining target seed's deficient-seed test.

Apply Xrefresh to the target. The deficient-seed bound and E_1 give
T_(v,j) uniform up to epsilon jointly with Z and the recorded prefixes.
The unrecorded smaller Q words are linear functions of Z with matrices
fixed by their recorded T values, so this uniformity also includes them.
No full-Q observation was needed to compute or describe the prefixes.

Now fix the target seed and still equal-prefix row seeds; their A
response components are exposed separately. For a good target seed its
B response is uniform independently of the remaining A source. Before
recording that B value, apply the Yrefresh good-direction test; at most
2epsilon of seed values fail. Record T_(v,j) and the other B responses
after their A parts. Each current Q word is now a LINEAR function of
its original Y row, with a matrix fixed by its T value. The target Q
is exactly uniform of length Qlen. Other Q words may be arbitrarily
linearly correlated. None is yet conditioned as a whole merely because
it has been computed.

The already recorded smaller prefixes are retained. Current equal-prefix
seeds are not recorded before the weak-seed Xrefresh: they could equal
the entire target seed. These distinctions are required for the stated
entropy and independence assertions.

### 5.3 The first differing digit gives an earlier response slot

Apply the alternating argument of Section 3 to the target and the
currently equal-prefix smaller rows. Run through the target slot
r=a_(v,j). Fix the other rows' seeds at each step and only responses
strictly before r before obtaining the target R_(v,j,r). The earlier
separated rows have already been handled by Section 5.2.

A newly separated smaller row has a_(i,j)<a_(v,j). Its selected response
therefore occurs strictly before slot r and is among the recorded past
responses. For a still equal-prefix row the two slots are equal; its
simultaneous response is not included in a claim of target independence.
Thus the target selected response is uniform jointly with the selected
responses of all i in L_j, up to at most 8h*epsilon plus the preceding
refresh exceptions. The target is also independent of the remaining
A variables at this X step, including the original Y_v.

The first differing digit is the entire reason that index advice is
present. The argument does not require the initial Y rows to be
independent, distinct or merely weakly correlated.

### 5.4 Refreshing the target seed for the next digit

First expose ALL current Q words for rows in the current tuple,
including the target, after obtaining the target's selected Xalt response
and BEFORE Ynext's good-direction test. The target response is uniform
jointly with the remaining A variables, including those Q words.
Every Q matrix was fixed by its T seed BEFORE the selected response
was obtained, so fixing these at most h*Qlen linear A bits changes the
Y_v direction independently of the fresh target B-side seed. Charge
them to B_0's h*Qlen term. These observations are postponed here,
not omitted forever: their later inclusion preserves the full past-Q
transcript required at the next-stage beginning assertion.
The earlier alternating prefixes were exposed seed-before-response
only through their specified slots. No unused later target response
is included.

For each newly separated smaller row, its selected response was obtained
before the target response, so its Ynext seed is already fixed. The
Ynext result is then a linear observation of its original Y_i. Fix these
next seeds, charging at most h*d of A rank. The matrices of those
observations were fixed BEFORE the target response was obtained, so
the remaining Y_v direction does not depend on the target seed.
The same is true of rows separated at an earlier prefix: their current
transcripts were processed before the target refresh.

Fix A components of the simultaneous responses from still equal-prefix
rows before their B components. They consume the already charged short
A observations but impose no seed-dependent new restriction on Y_v
when B response values are subsequently recorded. The target selected
response is a uniform B-side seed independent of the remaining affine
Y_v source. At most 2epsilon of seed values fail Ynext's rank test.
For each other seed the next target P_(v,j+1) is exactly uniform even
together with arbitrary remaining B-side response information. This is
Section 2's linear good-seed statement; it is NOT the false assertion
that the target response stays uniform after another equal response
is revealed.

To pass to the next stage, retain the target's CURRENT seed P_(v,j),
its T and Q observations and its alternating prefix up to its selected
slot. The Q observation was charged before Ynext's good-direction test as
specified above; it is not a new observation appended after that test. Retain O_(v,j) itself as the Ynext seed after applying
its good-direction test. Strong linear extraction gives uniformity of
P_(v,j+1) jointly WITH that seed, so recording it does not reveal the
new output. All these fibers remain affine because their extractor seeds
were fixed before their linear data outputs. No target responses later
than its selected slot are introduced.

We have proved the beginning assertion at stage j+1, including the
current next seeds of every index now in L_j and the target's own
charged past transcript. At the FINAL stage no next seed is needed:
the assertion that O_(v,b) is uniform excludes O_(v,b) itself from the
conditioning transcript. Its desired independence from smaller selected
outputs was established before this next-seed recording operation. All actual restrictions
were seed-exposed linear observations, so the rank budgets and affine
structure used in this induction are maintained. The initial assertion
at j=1 follows from uniform Y_v and the empty set L_0.

### 5.5 Error bound and final row independence

Per stage, the deficient-seed exception is at most epsilon, the two
refresh rank exceptions at most 4epsilon, and the h alternating steps
cost at most 8h*epsilon, with harmless further endpoint rank tests.
A bound 2^12*h^2*epsilon per stage dominates these charges and the
mixture bookkeeping. Errors add in joint distributions, and all future
evaluation is deterministic, so no uncharged conditioning amplifies
them. The entropy-defect amplification that is needed was paid
explicitly in Section 5.2 using E_1.

At j=b every smaller index has a different full advice string and
belongs to L_b. Hence

    (W_v,(W_i)_(i in T,i<v))

is within 2^12*b*h^2*epsilon of a uniform ell-bit W_v independent of
the displayed previous selected outputs. Apply this statement in order
to every element of T whose initial row is uniform. If every initial
row in T is uniform, the chain rule for variation distance yields a
joint error at most 2^12*b*h^3*epsilon from h independent uniform
ell-bit rows. The looser common bound

    gamma=2^20*b^2*h^4*2^-E

will be used below. Taking first bits and taking a smaller subset cannot
increase this variation distance.

## 6. From a single affine source to a non-oblivious bit-fixing source

Short is a strong linear affine extractor, so for any affine source X
of entropy at least k there is a FIXED set G of at least
(1-4N^-8)M good row indices. Each Y_i with i in G is exactly uniform.
The good set depends only on the direction space of X, not on its
random sample. Apply Section 5 to every at-most-h subset of G. It proves
that the first bits of W_i on G are gamma-close to h-wise independent
uniform bits in the marginal sense: each such small projection is
within gamma of the uniform cube. A separate proof coupling for each
projection suffices for this assertion; no globally consistent coupling
of all M rows to an independent matrix is assumed.

All other first bits may be arbitrary functions of X and correlated
with the good bits. This is allowed by resilient.md. Their number is
at most 4M*N^-8. Since M<=N^(2^20) eventually, for all sufficiently large
N this is at most M^(1-delta), delta=1/2^20. Indeed
4N^-8<=M^(-1/2^20) follows from M<=N^(2^20), which gives M^(-1/2^20)>=N^-1, and N large.

The main parameter proof establishes s<=2^21*L, z,w<=2^21*L
for the resilient metadata and R<=2^64*L. Its required marginal
order is at most Rw<=2^85*L^2, below h>=2^90*L^2. Its exact error
charge, including the polynomial prefactor in gamma, is at most
2^87*L^2 eventually, below E=2^90*L^2. Thus both numerical guards
pass with fixed slack. This uses the actual polynomial bound on M;
it does not infer the error condition from a leading liminf or an
unquantified small-error assertion.

Feeding these M first bits into F_M from resilient.md gives both output
values with positive probability. Thus the final deterministic Boolean
function is nonconstant on EVERY affine subspace of dimension at least
k, including every affine translation. The same fixed evaluator is used
for all such sources and all sufficiently large integer lengths.

## 7. Scope and implementation obligations

The argument is universal: it does not follow from checking sampled
indices, a bounded set of affine sources or the all-zero input. It pays
for correlated rows, equal advice prefixes, repeated or equal seed
values, the distinction between A and B response observations, and all
actual internal entropy thresholds. The final family must implement the
specified ordering and include source-independent numerical guards for
all code/design and rank-budget inequalities. A finite audit can only
verify its own bounded identities and small cases; it cannot replace
this conditioning argument.
