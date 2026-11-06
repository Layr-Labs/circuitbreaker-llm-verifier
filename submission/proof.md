# Fixed logarithmic-loss family and native full-B2 theorem

This file and its five companion proofs support the exact declaration in
bound.json. Mathematical arguments are supplied below.

## 0. Component map

The component map and provenance are retained in sampler-and-disperser.md.

# A full-B2 circuit-to-affine-flat bridge

All algebra is over F_2. A circuit is a finite acyclic, single-output DAG;
each gate is any one of the sixteen two-input Boolean tables and costs one.
Wires, reuse, depth, sharing and fan-out are unrestricted. Free constants
are allowed, a stronger model for the lower bound.

Attribution: the upstream source named in the supplied proof is GordoAR,
submission 6622943f-7b2d-4f8f-b73c-979c5232484a, version
226b74d1311ff6711c3c373ab7c8e1fc4d26eb0e, proof.md, Chapters II–III;
its public antecedents are 5e9fcf7f-8338-4729-9faa-a18994661af6 and
bf313f35-d97b-4836-b1e8-30afd1157ac2. The supplied proof attributes the
affine-system framework to Mishra–Santhanam, ECCC TR26-260 (2026-10-19).
This reconstruction uses the supplied circuit-reduction.md; upstream sources
were not independently inspected. No upstream theorem is assumed below.

## Lemma

Let d >= 0 be an integer. If f:F_2^n -> F_2 is nonconstant on every affine
subspace of dimension greater than d, then

    C_fullB2(f) >= 4n - 8(d+1).                         (1)

Inputs need not all be essential.

### Representation and potential

During the proof we permit signed references to signals. Signs entering a
gate are absorbed in its truth table. A constant, unary or repeated-parent
gate is replaced by a constant or signed earlier signal. These proof
representations do not assert free native input negation: we start with the
native count and need no conversion of the terminal representation back
to a circuit. The output can be a constant or signed parameter, auxiliary
or gate output.

The six inessential tables are two constants and four unary tables. The
remaining ten are XOR, XNOR and the eight tables

    c + (x+b)(y+e),     b,e,c in F_2.                   (2)

Call the latter product gates. These exhaust the basis.

Use M parameter coordinates u, of which 2q are in disjoint pairs. Each pair
ranges over three of its four corners, and every unpaired (free) coordinate
ranges over F_2. Write D for this product domain. Maintain an affine map
Phi:F_2^M -> F_2^n injective on the entire cube, and a representation of
f(Phi(u)) for every u in D. Equality outside D is not required.

Auxiliaries v satisfy a square system

    A v = N u + c,     A invertible,                   (3)

with at most three distinct nonconstant variables in each row, after
cancellation. A row has an auxiliary, since A is invertible. Remaining
binary gates read parameters, auxiliaries or earlier gates. Gate outputs
never enter (3). Let s be rows plus gates. A parameter's children are the
distinct rows and gates containing it, not pins and not the output reference.
Let I count free parameters having a child. Set

    a=M-q,             V=s+4I+2q.                     (4)

Initially M=n, q=0, Phi is the identity, the system is empty and s=|C|.
Thus a=n and V<=|C|+4n, including circuits omitting inputs.

For each pair, fixing one endpoint opposite its forbidden coordinate
leaves the other endpoint unrestricted. Such safe fixings give a flat
inside D of dimension a; Phi maps it to an a-flat in the original inputs.
Consequently a constant represented output is impossible if a>d. If
a>=d+2, an affine output is also impossible: on a safe a-flat its constant
fiber has dimension at least a-1. A product gate with affine parents cannot
be the output either: on that flat a nonconstant affine factor has a zero
fiber of dimension a-1, and constant factors give a constant output. This covers
constants, parameter/auxiliary literals and terminal gates with affine parents.

Until a<=d+1, lose eight units of V per lost unit of a and never increase
V at fixed a.

### Algebra and cost-free preparations

All coordinate changes below compose Phi with an affine injection. Replacing
a free x by w=x+ell(other parameters), with ell affine, is a bijection of
the whole coordinate cube and of the domain's free fibers; ell may involve
paired coordinates. Fixing a free coordinate lowers a by one. A safe fixing
of a paired endpoint, or restriction of its pair to any two allowed points
parametrized as a line, lowers M,q by one and preserves a. Restricting two
free coordinates to three corners increases q by one and lowers a by one.

Move an affine gate with no product-gate ancestor into (3), processing its
affine predecessors first. Its output becomes a new auxiliary with its
one defining row. The new auxiliary matrix is block triangular with new
diagonal entry one. The row has width at most three; s and parameter use
do not change. Degenerate gates are deleted instead.

The following simplifications decrease s by one, preserve a, and do not
increase I:

* Delete a constant or literal gate, substituting its output. This includes
  constant pins, equal/complementary pins, unary tables, and a gate with
  affine parents whose linear parts are dependent. In that last case a
  parent is constant or the parents are equal up to a bit, so its result is
  a constant or parent literal.
* Delete a gate with no child unless it supplies the output.
* For a row of width at most two, solve for one auxiliary v as a constant
  or signed other variable, substitute, and delete that row and column.
  Other widths do not increase. Add this row to every other row containing
  v, then expand the determinant along v's column: the remaining auxiliary
  matrix is invertible.

A substitution cannot give an unused free parameter its first child: that
parameter already occurred in the object supplying the substituted alias.
Thus the stated use accounting also holds through iterated propagation.

We need one further elimination. An auxiliary Q with coefficient one on a
free x has the form Q=x+ell(other parameters). Exchange their roles: Q is
free and x auxiliary, without expanding any equation. If b=A^(-1)N_x,
then b_Q=1. Replacing A's Q-column with N_x gives determinant
det(A)b_Q=1, preserving invertibility and row widths. Update Phi by
x=Q+ell. Both x before and Q after have a child, so the exchange preserves
I, s, M, q and V.

If free x has exactly one child, a row E, N_x is nonzero. Some Q therefore
has b_Q=1. Exchange x,Q. Now auxiliary x occurs only in E and, during the
loop a>=d+2, cannot be the output by the output exclusion. Delete
E and x, which E determines but nothing else uses. Expanding along
x's one-entry column proves invertibility of the remaining system. The
old free x is replaced by at most one used free Q, so I does not increase.
This deletion loses one row charge and preserves a.

### Clean all pairs

Perform available affine transfers and the preceding eliminations. Then
apply the first applicable pair rule below, propagate constant/literal
substitutions and short rows, and restart. A fixed parameter deletes each
of its old gate children, since a binary table with one fixed pin is unary,
and each old row child, since that row shortens to width at most two.
Count each removed object once. All four rules preserve a.

1. If endpoint p is unused, safely fix its partner h. The new free p is
   unused; q falls by one, I does not increase, and V falls by at least two.
2. If an endpoint has at least two children, safely fix it. At least two
   distinct objects disappear, and only its partner can become newly used
   and free. Hence Delta s<=-2, Delta I<=1, Delta q=-1, so Delta V<=0.
3. Otherwise each endpoint has exactly one child. If partners share that
   child B, its gate output or sole auxiliary is a function of the pair
   (a row B has just the endpoints and that auxiliary). Two of the three allowed points
   give the same bit. Restrict the pair to their affine line. Delete B
   after substituting that constant bit. No other object used the pair,
   so the new line coordinate is unused: Delta V<=-1-2=-3.
4. Otherwise suppose p's sole child is a row E. Its partner h has a distinct
   sole child B. Safely fix h and delete B. If propagation deletes E too,
   two distinct charges have been consumed. If not, now-free p still occurs
   only in E: its only original occurrence could not spread without E's
   deletion. Apply the one-child-row elimination. At most one new used
   free parameter arises, allowing the role exchange. As in rule 2,
   Delta V<=-2+4-2=0.

When these fail, every paired endpoint has exactly one gate child; partners
have different children, and paired parameters do not occur in (3).
In particular every auxiliary depends only on free parameters.

### Prepare a product source

Continue only while a>=d+2. If all gates transfer away, the output is affine,
a contradiction. Choose a product gate G with no product ancestor; its
predecessors are affine and have been transferred. After the simplifications
its parents have nonzero independent linear parts. Turn any auxiliary
parent Q into a free parameter by a role exchange. Such a free pivot exists
because Q is nonconstant and all paired columns of N vanish. If the other
parent is free, independence supplies a pivot different from it; if paired,
Q is independent of that coordinate. With two auxiliary parents, exchange
the first and then use this same argument for the second.

After each exchange, check for simplifications and pair cleanup. If any
is possible, perform a completed action and restart. Otherwise keep this
same G: at most two exchanges make its parents distinct parameters. G is
not the output, by the product-output exclusion, and has a gate child by
dead-gate cleanup. It cannot have a row child.

### Restrictions and the sharing cases

Each restriction lowers a by one. The new safe dimension
is at least d+1, so the restricted output is nonconstant. Freeze old gates
and incidences immediately before the restriction. Propagate topologically,
recording old outputs as constants or signed aliases. Every old gate hit
by a constant pin is removed. If an old output becomes constant, all its
old children are removed too, even if it first became an alias. Charges
refer to distinct removed old objects, not to propagation events.

**Two free parents x,y.** Restrict to the three corners on which G has its
majority value, creating a pair. G and its k>=1 distinct children disappear,
even with shared descendants or children that also read x or y. Both x,y
leave I; no unused old free parameter acquires a first use. Thus

    Delta q=1, Delta M=0, Delta I<=-2, Delta s<=-(1+k),
    Delta a=-1,        Delta V<=-(1+k)-8+2<=-8.        (5)

**Paired p and free x.** Write h for p's partner and B for h's sole gate
child; B!=G. Safely fix h and fix x to kill one factor in (2) at G. G is
constant, B disappears, and newly free p loses its sole use. If G has a
child other than B, it is a third distinct deletion. Otherwise B is G's
sole child and reads the two constants h,G, so B is constant too. It cannot
be the output, since the restricted output is nonconstant. Its old child
exists by cleanup, differs from B by acyclicity and from G because G reads
only parameters, and is a third distinct deletion. Hence

    Delta M=-2, Delta q=-1, Delta I<=-1, Delta s<=-3,
    Delta a=-1,        Delta V<=-3-4-2=-9.            (6)

**Paired p and paired r.** By pair rule 3 their pairs differ. Let h,k be
their partners. Safely fix h,k and then fix
now-free p to kill G's product. G is constant and r is unused. Write B_h,B_k
for the partners' sole gate children, neither equal to G. Four distinct
old gates disappear, as follows.

If B_h!=B_k, already G,B_h,B_k disappear. A child of G outside this set is
a fourth. Otherwise choose a child B_h of G. It reads h,G and is constant.
If it has a child outside the set, that is fourth; otherwise its only
possible child is B_k, which then reads k,B_h and is also constant. None
of these constant gates is the output. Cleanup supplies a further child,
and acyclicity and G's parameter-only inputs put that child outside the
three-gate set.

If B_h=B_k=B, B reads h,k and is constant. G and B read only parameters,
so neither reads the other. Neither is the output; each has a child.
If different children can be chosen, these provide two further deletions.
Otherwise both child sets are the same singleton {J}; J reads G,B and is
constant. J is not the output, so it has a child, distinct from J by
acyclicity and from G,B by their parameter-only inputs. This is fourth.

Only p,r become free; p is then fixed and r unused. Therefore

    Delta M=-3, Delta q=-2, Delta I<=0, Delta s<=-4,
    Delta a=-1,        Delta V<=-4-4=-8.              (7)

These cases exhaust the overlaps of partner children and descendants.
Every charge belongs to one restriction and one distinct old gate.

### Termination and conclusion

Every completed action preserves representation on the current D and
strictly decreases the nonnegative integer M+s: eliminations decrease s,
pair rules decrease M, and (5)–(7) decrease s. Neither M nor s increases.
Preparatory transfers decrease the number of gates; there are finitely
many. A fixed-source preparation uses at most two role exchanges before
an action. Thus preparation cannot cycle. The geometric exclusions prevent
termination inside a>=d+2, so the process reaches a_T<=d+1. Each completed
action satisfies

    V_before - V_after >= 8(a_before-a_after).

As V_T>=0, summing along the one history gives

    |C|+4n >= V_start >= V_T+8(n-a_T) >= 8(n-d-1),

which proves (1). If n<=d+1 the right side of (1) is nonpositive and no
process is needed. Constants, literals, unused inputs and repeated pins are
covered. Each unit of s is an original gate charge, moved or consumed once.

# Exact finite assembly of the two affine breakers

The dependency order is noncircular:

1. Finite seeded-linear extraction (Lemma P, linear-and-bilinear-primitives.md
   Section 7) and the constant-advice breaker CA# (restricted-nonmalleability-
   and-breaking.md Section 8.1).
2. The affine transformation of affine-and-standard-wrappers.md Section 2.1 at
   degree 20 around CA#: NB's inner breaker F# (Section 4.2 there).
3. The CZ/IP advice condenser (Lemma U), the AGHP extractor, Ext0, Ext2 and F#: the
   standard breaker NB (affine-and-standard-wrappers.md Section 4).
4. The same transformation at the sampler's degree with C=NB and the slots
   A# of outer Section 3 (NB's seed and threshold by Lemma B9 of Section 5):
   the outer generalized affine breaker.

Every constant below is specified by integer/rational operations or a search
whose termination is proved. No search tests an extractor property or a source
distribution. Logarithms are binary. Let lg(x) be the least nonnegative integer
j with 2^j >= x, for positive rational x. All real-looking powers appearing in
error bounds are analytical bounds, not instructions to compute floating-point
numbers. All chosen parameters can be computed with exact integer arithmetic.

## 1. Inputs and their precise scope

The finite theorem in `../circuit-log-primitives-20261105T0130/seeded-linear.md`
is used with its exact K recurrence and exact seed length, not an imported
GUV/CGL asymptotic theorem. Its output may be truncated, and its seed may be
padded by ignoring a suffix. `constant-advice-breaker.md` Section 8 supplies
a one-tampering standard breaker for arbitrary correlated source and seed
pairs, retaining both complete seeds.

The deterministic bilinear primitive and the AGHP extractor with its weak-seed
lift are the sibling modules `bilinear-extraction.md` and `weak-seed-hashing.md`
in that same primitives directory. The affine composition is Sections 1--4 of
`../circuit-log-extractor-research-20261105T0106/notes/affine-transformation.md`
(included as affine-and-standard-wrappers.md Sections 1--2).
The standard breaker is Theorem NB (affine-and-standard-wrappers.md Section 4);
it replaces the earlier equal-length and general-input wrappers. Theorem D3
and Corollary D3' (linear-and-bilinear-primitives.md Section 7) extend the
seeded-linear theorem (j-way splits, fixed R, padding to U1, admissible field
dimensions); the base sampler (L) and Lemma P use them.

For the advice condenser we use `cz-conditional-entropy.md` Section 8 with
eta=1/6711000 (671eta<gamma), and its expansion input from
`finite-approximate-ring.md` Section 5: gamma=1/10000 and prime threshold
2^20000. The CZ source and prime
size guards below imply that threshold. Its canonical prime setup is the one
in `prime-fields-and-encoding.md`. From `ordinary-condenser.md` only
its scalar sum-product estimate (30) (Section 6) is used, by Lemma RG. These
are declared mathematical inputs, not proved again or strengthened here.
The ring theorem has exactly the expansion contract left conditional in the
CZ file. Substitution with the displayed gamma and eta closes that premise.

This note does NOT assume normalized contract (F) in
`recursive-sampler-schedule.md`, nor any closed-form degree envelope: the
outer breaker's seed and threshold are the slots d_O and k_O of the outer A#
(outer Section 3), Lemma B9 bounds NB's seed and threshold at the outer call,
and the outer unit Om is a closed form.

## 2. A reusable exact seeded-extractor envelope

For integer height h>=11, B=2^h, source length u>=1 and error exponent e>=1,
define

    b=e+h+4,
    rho_u=b+2 lg(4u)+2,
    S_u=least power of two >=(4h+1)rho_u,
    G_u=2B(4b+S_u+5),
    P_u=(h+7)S_u+18b+18.

Instantiate the finite recursive primitive for output z>=1 with leaf output
ceil(z/B), then truncate. Its exact threshold K(u,z) and seed length D(u,z)
are computed by its own algorithm and satisfy

    K(u,z)<=2z+G_u,       D(u,z)<=6z/B+P_u.                 (2.1)

The call is legal whenever its computed K(u,z)<=u. The inequality
128 P_u <= G_u holds if B>=64(h+25), since
P_u<=(h+25)(S_u+b+1) and G_u>=2B(S_u+b+1).

A useful envelope, valid with the definitions above, is

    G_u <= C(e+h+lg(u+1)+6),    C=64B(h+1).                (2.2)

Indeed rho_u=e+h+10+2lg u <=2(e+h+lg(u+1)+6),
S_u<2(4h+1)rho_u, and
4b+S_u+5 <= (16h+9)(e+h+lg(u+1)+6).
Both G and P are nondecreasing in u. No unspecified extractor constant
appears in these envelopes. The selected family does not call (2.1) or
(2.2): every seeded primitive of its breakers is a Lemma P primitive
(linear-and-bilinear-primitives.md Section 7).

## 3. The affine assemblies and the standard-breaker interface

### 3.1 The standard interface required here

A standard one-tampering breaker C is used through its actual seed length
d_C and threshold k_C at each call: arbitrary correlated source pair
independent of arbitrary correlated seed pair, honest seed uniform, distinct
advice, both seeds retained, error <=2^-e; its guards are explicit (k_C<=u
for CA#; k_C<=u and d2<=u for NB). No interface constant c is used. The
constant-advice breaker CA# (restricted-nonmalleability-and-breaking.md
Section 8.1) is used only at its one fixed advice length 4, inside NB's
inner breaker F#. The standard breaker NB (affine-and-standard-wrappers.md
Section 4) works at every a>=1.

### 3.2 The two affine assemblies

The assemblies are instances of affine-and-standard-wrappers.md Section 2.1
with the slots of F# (Section 4.2 there) and of the outer A# (outer Section 3);
no selector remains. If a guard of an assembly or breaker fails (k>n; for NB,
k_C>n or d2>n), the map outputs the all-zero string and no claim is made (the
totality convention).

## 4. The advice-condenser constant

NB uses Lemma U (affine-and-standard-wrappers.md Section 3.5) at
eta=1/6711000, the CZ rate allowed by cz-conditional-entropy.md Section 2 (1)
with gamma=1/10000 (671eta<1/10000), with error parameter sig=Ee+4 at each
call; its guards are NB's (I2). No other condenser constant is used.

## 5. The standard breaker NB and Lemma B9

Theorem NB (affine-and-standard-wrappers.md Section 4) is a standard breaker
in the sense of Section 3.1, with guards k_C<=n and d2<=n, at every a,z,E>=1.
In the dependency order above it uses CA# inside F#, the advice condenser of
Lemma U at the rate of Section 4, the AGHP extractor (with Lemma W2) and the
primitives Ext0, Ext2.
The outer assembly uses NB only through the outer A# (outer Section 3), and
Lemma B9 gives NB's seed and threshold at the outer call.

**Lemma B9 (NB at the outer call).** At the outer call (outer Section 3) NB is
called at source length n'=d_O<=n (the outer seed), output z=dx_O=203*Om,
advice a=pL+2+lg C<=Om and exponent e, where Om=e+2L+40 is the outer unit.
Then

    d_C < 262915*Om,     k_C < 5437485715179*Om,

and NB's guards hold (k_C<=n' is the outer ledger line "k_C <= d";
d2<=d_C<=sC<=d_O=n', sC being the outer C-seed slot).

Proof (NB's algorithm, affine-and-standard-wrappers.md Section 4.3;
Om>=2^40). Every additive constant below is absorbed by adding 10^-6 to a
coefficient, as 8lg(Om)+240<=10^-6*Om. Ee=e+4<=Om, sig=Ee+4<=Om, and NB's unit is
w_N=Ee+32+6+2lg(4n')<=Om+14<=(1+10^-6)Om (lg(4n')<=L+2). Coefficients
of Om:
Step 1: F# runs with unit w_F=Ee+2lg(4n')+32<=Om (w_F>2^20), so
m1=d_F<=80884w_F+1277283z/2^20+1<=dF*Om and
k_F<=173878w_F+35129230z/2^20+1<=kF*Om with
dF=1329254838032009/16384000000 [81131.2768] and kF=1480121339211317/8192000000 [180678.8745]
(each +10^-6).
Step 2 (Lemma P for Ext2=P[9,20] at unit w_N: c2=105/38912, pi2=45875/152,
Lam2=512000000000/322687697779 [1.5867], g2=7929856000000000/322687697779 [24574.40]):
m0<c2*m1+pi2*w_N<=M0*Om with M0=c2*dF+pi2(1+10^-6)+10^-6
=17472928791820587/33554432000000 [520.7339], and
K2<=Lam2*m1+g2*w_N<=[153303.27]*Om.
Step 3 (c0, pi0 of Ext0 as in affine-and-standard-wrappers.md Section 4.1,
eta=1/6711000): as m0<=M0*Om, w_N<=(1+10^-6)Om, Ee,sig<=Om and
12<=68*10^-6*Om, the first term of step 3's max is at most
(11(M0+pi0(1+10^-6)+3+4*10^-6)+23+24*10^-6)Om/(eta-110c0)+1<=R0*Om with
R0=(11(M0+pi0(1+10^-6)+3+4*10^-6)+23+24*10^-6)/(eta-110c0)+10^-6
=2356978245279602446984399837/10805871837000000
[218120136980.4478]; the other two terms are at most (4Om+2)/eta+1 and 64/eta+1, both
below R0*Om, so r<=R0*Om; v=10r.
Step 4 (Lemma P for Ext0=P[32,40] at unit w_N: Lam0=1844674407370955161600000000000000000000000000000000/820478914983440471652837469810666679283466230552321
[2.2483], g0=33671969068562343477256178892800000000000000000000000000000000/63113762691033882434833651523897436867958940811717
[533512305919.8]):
K00<=Lam0*v+g0*w_N<=[5437485715178.7]*Om and d0<c0*v+pi0*w_N<=[2428.90]*Om.
Steps 5-7: e_L<=(d0+Ee)/2+Ee<=d0/2+(3/2)Om; r_L<2(lg(r*m0+1)+m0+2e_L) and
lg(r*m0+1)<2lg(Om)+60 (R0*M0<2^47), so
d1=2r_L<=4(M0+[2428.90]+3)Om+10^-6*Om<=[11810.54]*Om;
d2=d0+d1+m1+K2+Ee<=[248674.99]*Om; d_C=d0+d1+d2<=[262914.43]*Om<262915*Om;
k_C=max(K00, k_F+2v+z+2Ee)<=max(K00, (kF+20R0+203+2+10^-6)Om), and
kF+20R0+203+2+10^-6<4362402920492.9<5437485715178.7<5437485715179.
Each bracket is an upper bound for the displayed expression evaluated with
the displayed values before it (exact fractions, or the earlier brackets).

## 6. Uniformity, exactness, and remaining obligations

No parameter search remains: CA#, F# and the outer A# use closed-form
slots, NB's steps are exact formulas, and constant formation is finite. Exact seed lengths and entropy thresholds
are computed using the original primitive recurrence, even though the
conservative envelopes decide a legal pad length. The evaluator never
receives entropy witnesses, good-row labels, source decompositions, coupling
coins, proof histories, or a per-length truth table.

For fixed t and outer-logarithmic m,a,E, auxiliary dimensions are O(log n).
The seeded-linear algorithms are polynomial in their dimensions and B;
the finite CZ/AGHP/bilinear setup costs 2^O(local dimension),
so their complete composition is polynomial in the actual outer n. Fixed
constants can be enormous. This does not assert polynomial bit-time jointly
in arbitrary t and n: the exponent of n may depend on t. The outer
assembly fixes t=4C-1 once, with its constants.

For ordinary seeded primitives or standard breakers with the complete
source tuple independent of the complete seed tuple, unused suffixes may be
retained by their conditional kernel given all used prefixes. The outputs
depend only on the source tuple and those prefixes. This applies to the
primitive pads of CA#, F# and the outer A# (restricted-nonmalleability-and-breaking.md
Section 8.1, affine-and-standard-wrappers.md Section 4.2, outer Section 3) and to the
zero-additive-shift independent-side use of NB's inner breaker
(affine-and-standard-wrappers.md Section 4.6, step 7). It does not assert a generic
kernel principle when arbitrary additive B-tuples correlate with suffixes.

# Tagged affine-permutation logarithmic-loss family: outer assembly

All constants below are fixed once; none is length-dependent advice.
The affine-breaker contract (F') is proved below from the outer A# (Section 3) and Lemma
B9.

This is an independently written composition of the finite seeded-linear,
scalar-growth, graph and affine-breaker constructions
credited in their companion modules. The sumset framework is due to
Chattopadhyay--Liao, https://arxiv.org/abs/2110.12652v1 , with the
advice-condenser breaker architecture following Xin Li,
https://eccc.weizmann.ac.il/report/2023/023/revision/1/ . The native full-B2
bridge is credited in minimal-circuit-bridge.md. No earlier accepted result
or model judgment is a mathematical premise.

All logarithms are binary. Write lg(x)=ceil(log2 x) for positive x, computed
by integer bit length and a power-of-two comparison. All string lengths,
indices and loop counts are integers. All arithmetic parameter algorithms
use exact integers and rational comparisons, not floating point.

## 1. The exact inputs to this assembly

(L) The base sampler's extractor is Theorem D3 of
linear-and-bilinear-primitives.md Section 7 at j=4, R=8, h'=32, leaf
10L, error exponent L, admissible dimensions and U1=2^(92+lg L), for
actual input length n and L=lg(n+1)>=2^25 (sampler-and-disperser.md Section
5.3). Section 5.3 there proves its seed D<419.275L, padded to s=B0*L with
B0=421 by ignoring a suffix, and its exact threshold K_0<kfac*L with
kfac=2^64(8/5)^32(10+(86/27)(1+2^-14)+10(1+2^-14)+2^-14),
ceil(kfac)=1455412008444481506303478516; its output
2^64*10*L>=delta_O*L is truncated to d=delta_O*L bits (delta_O of outer
Section 3). The actual guard K_0<=n is required (it follows from HL<=n). The
sampler at failure 1/100 uses entropy k_S=K_0+7<kfac*L+7. This extractor is
source-linear and strong, with a uniform polynomial evaluator.

(G') Theorem G' of sampler-and-disperser.md Section 4, with P=8970,
p=P*B0=3776370, C=4(q+1)=2516740, A=2^(pL) and Lambda_S=A^(1/2)/2^15 (Lambda
there): a deterministic graph Gamma from left pL-bit strings to right s-bit
strings, of degree exactly C, such that whenever L>=2^31, every set of at
least Lambda_S left vertices has more than 2^(1-L)*2^s distinct neighbors.
Gamma is the four rows of one AB+C step followed by a short-start Cayley
walk and an s-bit prefix; its uniform evaluator costs 2^O(pL). With the
base sampler of (L) it gives (23') there: for every k_S-source and every
test of density at most epsilon0=2^-L, fewer than Lambda_S rows are
all-bad, except on source probability 1/100.

(F') For every n passing the guards (5'') of outer Section 4, the outer A#
(outer Section 3) with C=NB, at the request (O2) and the unit Om, returns an
exact seed length d_F, source entropy requirement k_F and a one-bit
generalized affine breaker, where

  d_F=d_O<=delta_O*L,   k_F=k_O<=kappa_O*L.                                 (F')

For arbitrary correlated tuples (A0,...,At) independent of
(B0,...,Bt,Y0,...,Yt), with H_inf(A0)>=k_F, Y0 uniform on d_F bits, and
alpha0 different from every alphaj, its honest output on Xj=Aj+Bj is
2^-E-close to uniform jointly with every tampered output and EVERY seed
Y0,...,Yt. The field is F_2, so addition is XOR. Neither source copies nor
seed copies are required to be deterministic functions of their honest copy.
For the fixed t and a,E=O(log n), its complete evaluator is polynomial in
the actual n (the evaluator of affine-and-standard-wrappers.md Section 2.1
with NB, whose cost is stated in Section 4.7 there).

## D1. Algorithm and finite hypotheses

Write L=lg(n+1), s_tag=log2 A=pL and w=2, with p>=128 and L>=1024.
Build the canonical field F_(2^L), and split x into L-bit coefficients of
P_x, padding its final block with zeros. For a full pL-bit row alpha,
let beta,gamma be its first and second L-bit blocks. Define
H_alpha(x)=Prefix_2(gamma*P_x(beta)), a two-bit linear map in x.
Each (beta,gamma) has exactly A/2^(2L) row preimages. For Delta!=0,
P_Delta is nonzero of degree ceil(n/L)-1<n/L, so its root fraction is
at most n/(L*2^L)<=1/L. Off those roots, multiplication by P_Delta(beta)
permutes gamma and exactly 1/4 of gamma values have zero two-bit
prefix. Thus at most A*(1/L+1/4)<=257A/1024 rows have H_alpha(Delta)=0.
The root bound follows by factoring X-beta at distinct roots. The supplied
finite-field proof makes the canonical search halt in 2^O(L) time.
Full alpha remains in the advice label despite repeated (beta,gamma).
No independence of hash values across different rows is asserted.

Use the existing source-linear sampler Samp(x,alpha,j), j in [C], seed
length d, tolerating test density epsilon0=2^-L with fewer than Lambda_S
all-bad rows ((23') in (G')), failure probability 1/100, at source
entropy k_S. The arbitrary enumerated field index alpha replaces [A].
For t=4C-1 and a=s_tag+w+lg(C), request the padded all-seeds-retained affine
breaker F of output length one, entropy K, and error at most

    epsilon = 2^(-L)/(2^33 C A^6).

The advice for (x,alpha,j) is the concatenation of the s_tag-bit alpha,
lg(C)-bit j-1, and w-bit H_alpha(x). Define

    R_alpha(x)=XOR_(j=1)^C F(x,Samp(x,alpha,j),advice(x,alpha,j)),
    D(x)=1 iff 2 sum_alpha R_alpha(x)>A.

Assume the input source entropy satisfies

    k >= max(k_S, K+4w+4Cd+2lg(1/epsilon)).                 (D1)

All concatenations, zero padding and majority ties are specified. There
are no seed or advice inputs to D. The tags are internal computations on
the actual n-bit x.

## D2. Bad seeds with arbitrary correlated left copies

Fix an arbitrary correlated tuple (A_0,...,A_r), fixed advice differing
from the honest advice, and deterministic seed maps L_i of A_i. Assume
Hinf(A_0)>=K+ell+lg(1/epsilon), where the joint seed leakage has at most
2^ell values. Call w bad if some b_0,...,b_r,w_1,...,w_r make the honest
bit in the tuple

  F(A_0+b_0,w+L_0(A_0),advice_0),
  (F(A_i+b_i,w_i+L_i(A_i),advice_i))_(i=1)^r

more than theta away from a fresh independent uniform bit jointly with
the other outputs. Its bad-seed density is at most2epsilon/theta.

Choose one witnessing tuple for each bad w. With W uniform independent
of the full left tuple, set B_i=b_i(W), Y_0=W+L_0(A_0), and
Y_i=w_i(W)+L_i(A_i). Conditional on all L_i(A_i), the complete left and
right tuples are independent, and Y_0 is uniform. The bounded-range
inequality in the retained sampler proof loses at most ell bits of average
conditional entropy of A_0, even when the leakage involves the other left
copies. Apply its conditional lift with breaker error epsilon and
exceptional probability epsilon. The retained honest seed and leakage
determine W. Averaging the local distances over W, whose independence
leaves the original left-tuple distribution unchanged, gives the claimed
density. The arbitrary-copy quantifiers and full seed retention are those
of the supplied actual contract (F'); no copies need entropy except A_0.
All witnesses are finite proof objects, never evaluator searches.

## D3. Actual parities and pointwise collision counts

Fix independent sources X_1,X_2 of entropy at least k, and write
A_0=X_1, A_1=Lambda X_1, B_0=X_2, B_1=Lambda X_2+b, with Lambda
invertible. Suppose tau(x)=Lambda x+b has no fixed point on the sum
support. Row alpha collides when H_alpha(x_1+x_2)=H_alpha(tau(x_1+x_2)).
The tag root bound gives at most N_col=257A/1024 collisions at EVERY point
(x_1,x_2). We use actual bits throughout; no random replacement is made.

For each nonempty set Q of at most four distinct row-copy indices, expand
its at most4C calls, hence at most t=4C-1 tampered calls. If some alpha occurs in only one copy in Q, choose
such an index as honest. Otherwise every alpha in Q occurs in BOTH
copies; call such Q fully doubled and choose its least index as honest.
Enumerate the C honest neighbor choices. Condition on just the SELECTED
left tags (H_alpha(A_sigma))_((sigma,alpha) in Q), at most 2^(4w) cells.
Discard cells of mass <epsilon/2^(4w), of total probability at most epsilon.
On every remaining cell X_1, and also Lambda X_1, has honest entropy at
least k-4w-lg(1/epsilon). Enumerate just the selected right tags,
at most 2^(4w) choices. All selected advice strings are then fixed.

If Q has a singleton alpha, its honest advice label alpha,j differs from
all other selected calls' labels, regardless of collisions at other rows.
If Q is fully doubled, the selected profiles determine all its collisions.
On profiles with no selected collision the honest advice is again distinct
from all tampered advice; two copies of the same label have distinct tags.
Tampered advice strings may coincide with one another in either case.
Section D2 therefore applies in all these cases, with joint seed leakage
at most4Cd. There are at most sum_(i=1)^4 binom(2A,i)<=64A^4 sets Q,
and hence at most64C A^4*2^(8w) instances including both profiles and honest
choices. It suffices that

  k>=max(k_S,K+4w+4Cd+2lg(1/epsilon)).              (D1')

Choose theta=1/(1024A^2), epsilon=2^(-L)/(2^33*C A^6), and let B be the
union of their bad-seed sets. As 64*2^16*2*1024=2^33, Section D2 gives

  |B|/2^d<=64C A^4*2^(8w)*2epsilon/theta=2^33 C A^6 epsilon=2^(-L).    (D2')

B depends on the left source law and tau, not the right source point;
the program never constructs B. Apply the sampler to both B_0 and B_1,
each of entropy k. Except on probability2/100 both conclusions hold.
Fix such an x_2. Intersect the good row sets to obtain G, omitting
q<=2(Lambda_S-1)<2^(-14)A^(1/2) rows. Since 2^(-14)A^(-1/2)<=1/256,
g=|G|>=255A/256>=A/2. Each sigma,alpha in {0,1}xG has its own
neighbor j_(sigma,alpha) whose honest seed lies outside B. The two copies
need not choose the same j. Source-linearity gives every seed identity
with shifts B_sigma=x_2 or Lambda x_2+b and offsets Samp(B_sigma,alpha,j).
The second sampler is essential; its seeds need not be fixed translates
of the first sampler's seeds.

Put eta=2(theta+epsilon)<=3/(1024A^2). For each Q on G, using its honest
index's good j, the actual parity has absolute expectation at most eta
if Q is not fully doubled. If fully doubled, the bound is instead

  eta+Pr[at least one selected alpha collides].       (D3a)

Indeed good cells/profiles without collisions have bias at most2theta;
exceptional entropy cells have total mass at most epsilon and collision
cells cost at most their probability. These costs sum to at most
2theta+epsilon+Pr[collision]<=eta+Pr[collision].
Bad-cell bounds apply individually to each Q; simultaneous good cells
are not assumed. All probabilities here are over X_1 for the fixed x_2.
The pointwise collision bound implies sum_(alpha in G) Pr[collision at alpha]
<=N_col. No independence of rows, copies, or collision events is asserted.

## D4. A robust quartic sign certificate

Let X_alpha,Y_alpha be the actual signs in the two copies on G, with
S_0=sum X_alpha and S_1=sum Y_alpha. Put Z_alpha=X_alpha Y_alpha,
Z=sum Z_alpha and zeta=E Z. Section D3 gives

  |zeta|<=g*eta+N_col.

In S_0^3 S_1, the triples reducing to one row have coefficient 3g-2.
Every other surviving parity has a singleton row label, hence bias <=eta.
The same applies to S_0 S_1^3. In S_0^2 S_1^2, all-even tuples contribute
g^2. Fully doubled two-row sets contribute
4 sum_(alpha<beta) E Z_alpha Z_beta=2(E Z^2-g).
Every other surviving parity again has a singleton. Thus, with
|e31|,|e13|,|e22|<=g^4 eta,

  E S_0^3 S_1=(3g-2)zeta+e31,
  E S_0 S_1^3=(3g-2)zeta+e13,
  E S_0^2 S_1^2=g^2-2g+2E Z^2+e22.                 (D3')

Marginal fourth moments have only single-copy surviving parities, so
E S_sigma^4<=3g^2+g^4 eta<=4g^2. No collision error is needed there.
The elementary variance identity gives E Z^2>=(E Z)^2=zeta^2.

Define Q_+(x,y)=xy(x-y)^2 and Q_-(x,y)=-xy(x+y)^2.
Both have central term -2x^2y^2. Set r=|zeta|/g. Since
N_col/g<=(257/1024)*(256/255)=257/1020 and eta<=3/(1024A^2),
we have 0<=r<13/50. Both quartic expectations are at most

  2(3g-2)|zeta|-2g^2+4g-4zeta^2+4g^4 eta
  <=[6r-4r^2-2+4/g+12/1024]g^2<=-g^2/2.          (D4')

For the final inequality, 6r-4r^2 increases on [0,13/50],
and 6r-4r^2-2 equals -444/625 at13/50; also 4/g<1/1000.
This uses the exact matched-pair second moment rather than bounding each
collision term separately. All identities concern actual deterministic signs.

Let Y_sigma be the ACTUAL sum over ALL rows in that copy. Pointwise
|Y_sigma-S_sigma|<=delta=q<2^(-14)A^(1/2). Since g>=A/2,

  delta/sqrt(g)<2^(-14)A^(1/2)/(A/2)^(1/2)=2^(-27/2)<1/10000.

The L4 norms of S_sigma are at most2sqrt(g), and those of Y_sigma at
most3sqrt(g). The norm triangle inequality follows from Holder applied
to (|X|+|Z|)|X+Z|^3 and division by the cube of the L4 norm (zero norm
is trivial). For four factors, Holder itself follows by Cauchy-Schwarz
twice: E|abcd|<=(E a^4 E b^4 E c^4 E d^4)^(1/4).

Telescoping the four factors in any degree-four monomial shows its
expectation changes by at most4delta(3sqrt(g))^3. Each Q has sum of
absolute coefficients four. Therefore

  |E Q_+(Y)-E Q_+(S)|, |E Q_-(Y)-E Q_-(S)|
    <=432delta g^(3/2)<(432/10000)g^2<g^2/8.            (D5')

Both actual expectations are strictly negative. If the actual majority
XOR were identically zero, Y_0Y_1>=0 pointwise and Q_+(Y)>=0, impossible.
If that XOR were identically one, Y_0Y_1<=0 and Q_-(Y)>=0, also impossible.
The weak product inequalities include all ties: the evaluator outputs one
exactly for a negative sign sum, and zero for a nonnegative sum.
Thus both actual majority-XOR values occur. No sixth or eighth moments,
small-tail inference, or independence of the full sign vector is used.

All sums and collision events above are actual deterministic quantities.
For an affine S=a+V of dimension at least k, take X_1 uniform on V and
X_2 independently uniform on a+V. Their sum lies in S. Good x_2 have
probability at least98/100. Fix any such x_2; the preceding contradiction
proves both values of D(x)+D(tau(x)) occur on S. Thus the property holds
for every invertible affine tau without fixed points on S, and includes
every nonzero translation.


## 2. Fix every constant before seeing n

The base sampler (L) is fixed (j=4, R=8, h'=32, leaf 10L, B0=421). Put

  P3=2990,  P=3*P3=8970,  p=P*B0=3776370,
  q=floor((P3-1)B0/2)=629184,  C=4(q+1)=2516740 (the sampler degree),
  t=4C-1=10066959,  j=lg t=24 (2^23<t<=2^24).

These constants and the slots of Section 3 depend on no input length,
source, experiment or oracle.

For each actual n>=1 let L=lg(n+1), A=2^(pL), s=B0L and d=delta_O*L (delta_O
from Section 3). Request

  t=4C-1, m=1,
  a=pL+w+lg(C), w=2,
  E=(6p+1)L+lg(2^33*C).                                       (O2)

Here lg C=22 (2^21<C<=2^22) and lg(26t)=28 (26t=261740934 lies in
(2^27,2^28]), so A#'s exponent is e=E+lg(26t)=(6p+1)L+e0 with
e0=33+22+28=83. Then 2^-E<=2^(-L)/(2^33*C A^6)=epsilon. The advice
concatenates the pL-bit alpha, lg(C)-bit j-1, and w-bit polynomial-evaluation
tag from Section D1. It has exactly a bits and is evaluated internally.
Below, w=2 is the tag length and Om is A#'s unit.

## 3. The outer affine breaker A# and (F')

For actual n with L=lg(n+1) passing the guards (5''), put E and e as in
(O2) (e=E+lg(26t)=(6p+1)L+e0), the outer unit Om=e+2L+40<=(6p+4)L
(L>=e0+40=123), q=e, j=lg t=24, and

    dy=68Om,  dx=203Om,  r=684553404Om,  sC=262915Om,  s0=132Om,
    kC=5437485715179Om,  d_O=8085569343327Om,  k_O=6898917797719389Om.

**Lemma A#.** Run the evaluator of affine-and-standard-wrappers.md Section
2.1 at degree t=4C-1 with C=NB (output dx, advice a of (O2), exponent e, seed
padded from d_C to sC, Section 2.1's d0), L0=P[16,4] (input n, output sC, seed
s0), Lr=P[27,4] (input n, output r, seed dx), E=P[6,4] (input d_O, output dx, seed
dy), Lm=P[1,6] (input r, output dy, seed dx), all Lemma P primitives at exponent
e; output m=1<=r. Then (F') holds with d_F=d_O and k_F=k_O:

    d_F<=delta_O*L,  delta_O=8085569343327(6p+4)=183204641348636071248,
    k_F<=kappa_O*L,  kappa_O=6898917797719389(6p+4)=156317224818312605105136.

Proof. Section 2.1's contract (3) and legality are the ledger below (slot
convention of Lemma P with unit Om; valid for Om>=998272, and
Om>=(6p+1)L>=2^50), with T_R=(t+1)(s0+sC+dx)+j(t+1)(2dy+2dx) and
T_L=(t+1)sC+2j(t+1)dy expanded into their slots. NB's seed and threshold at
this call are Lemma B9's. Every primitive input is n, d_O or r, all <=n since
d_O<=delta_O*L<=HL<=n and r<=k_O; so (P) holds with w=Om
(e+h+6+2lg(4u)<=e+2L+38 for h<=28). (4) gives error
(13*2^j-8)2^-e<26t*2^-e<=2^-E. d_O=8085569343327Om<=delta_O*L and k_O<=kappa_O*L
because Om<=(6p+4)L. A# has no search: Om, every slot and every primitive are
closed forms.

    dy >= D_E(d,dx)                        [67.6498]Om + [0.033] <= 68Om
    dx >= D_r(n,r)                         [202.9174]Om + [0.001] <= 203Om
    dx >= D_m(r,dy)                        [96]Om + [0.875] <= 203Om
    s0 >= D_0(n,sC)                        [131.8996]Om + [0.001] <= 132Om
    (3a) r >= (t+1)dy + K_m + q            [684553393.2000]Om + [10066961.200] <= 684553404Om
    dy <= r                                [68]Om + [1] <= 684553404Om
    K_m <= r                               [112.2000]Om + [1.200] <= 684553404Om
    (3c) d >= (t+1)(s0+sC) + k_C + q       [8085569342300]Om + [20133921] <= 8085569343327Om
    (3d) d >= T_R + K_E + t dx + q         [2783121833184.4253]Om + [1006696004.619] <= 8085569343327Om
    s0 <= d                                [132]Om + [1] <= 8085569343327Om
    K_E <= d                               [4826.4253]Om + [5.619] <= 8085569343327Om
    k_C <= d                               [5437485715179]Om + [1] <= 8085569343327Om
    (3b) k >= T_L + K_r + t r + q          [6898917797718884.7779]Om + [503350361.359] <= 6898917797719389Om
    k >= K_0(n,sC)                         [93255182.8492]Om + [99.775] <= 6898917797719389Om

The outer evaluator computes Om and the slots of this section in closed
form, with no search, so (F') holds at every n passing (5''). Pad the breaker
to d-bit seeds by taking its first d_F bits. All other finite guards are still
checked explicitly by the evaluator.

All original d-bit seeds remain retained by the following consequence of
(F')'s universal right-side quantifiers. A simple kernel argument given only
the used seeds would not suffice when additive shifts correlate with suffixes.

Fix the admissible left tuple A. For an honest raw seed y and a deterministic
profile r=(B0,...,Bt,Y1,...,Yt), let Delta(y,r) be the TV between the output
tuple over A and an independent uniform honest output together with the
unchanged tampered-output tuple. For each y choose a maximizing profile r_y
from this finite domain. Uniform y with profile r_y is a legal (F') instance,
since the whole right side is independent of A and every retained seed is a
function of y. Retention of the honest seed gives

  2^(-d_F)*sum_y max_r Delta(y,r) <= 2^-E.

For an actual full right-side variable R independent of A jointly with the
uniform honest raw seed y, conditioning on (y,R) fixes a profile. The joint
output distance retaining (y,R) is the average of its Delta(y,r) values,
and is therefore bounded by the displayed average of maxima. Thus the
WHOLE right side, including all padded seed suffixes, may be retained with
no extra error. No independence among those suffixes is assumed. The finite
maxima are proof devices, never evaluator searches. This proves the padded
d-seed guarantee even when B-tuples correlate arbitrarily with suffixes.

## 4. One fixed all-length guard and entropy coefficient

Put

  H=max(4C*delta_O+kappa_O+12p+3, ceil(kfac)+7).                       (O3)

Here 4C*delta_O=1844313796271065383810766080, kappa_O=156317224818312605105136,
12p+3=45316443 and ceil(kfac)+7=1455412008444481506303478523 (from (L)), so
H=1844470113495883696461187659.
H is a fixed positive integer. Lemma A# gives K=k_F<=kappa_O*L. Since w=2,

  4w+2lg(1/epsilon)=8+2((6p+1)L+33+lg C)=(12p+2)L+74+2lg C<=(12p+3)L

when L>=74+2lg C. Thus, with d=delta_O*L,

  K+4w+4Cd+2lg(1/epsilon)
    <= [kappa_O+4C*delta_O+12p+3]L<=HL.

Also k_S<kfac*L+7<=(ceil(kfac)+7)L<=HL by (L) and (O3). Hence k=HL supplies (D1). The exact
primitive entropy guard K_0<=n follows from HL<=n. Since p=8970*B0>=128,
A>=2^128 and A>=(n+1)^128 as required.

The evaluator tests these literal finite guards:

  n>=1;  L>=2^31;  HL<=n.                                      (5'')

L>=2^31 implies every other lower bound on L used above: 1024,
e0+40=123, 74+2lg C=118, 58 (Lemma W'), the base guard L>=2^25 of (L), and the
walk guard (W') of sampler-and-disperser.md Section 3 (true for every L>=2^31). The evaluator also computes and checks
the supplied algorithms' exact parameters and guards. The already proved
specialization (L), (G') and (F') guarantee that these checks pass when
(5'') holds. No test of an unknown extractor/graph property, source
distribution, or existence witness occurs.

There is one fixed threshold N beyond which (5'') holds for EVERY integer n.
Indeed L tends to infinity with n, and the only length guard besides HL<=n
is L>=2^31 (Lemma W'), so it holds eventually. Finally
H*lg(n+1)<=n eventually: choose j=floor(log2 n), so 2^j<=n<2^(j+1) and
L<=j+2. The elementary exponential sequence 2^j/(j+2) is increasing with
consecutive ratio 2(j+2)/(j+3)>=4/3, hence unbounded, so it tends to infinity:
H(j+2)<=2^j for every j from some fixed j0 on. This also proves eventual
k_F<=n.

This is an all-integer-length argument, not a favorable subsequence or a
liminf used to invent a remainder. A small-input zero branch is allowed
because only an eventual lower bound is claimed.

## 5. The fixed evaluator

Precompute the constants in Section 2 by their finite algorithms. On an
input x of actual length n:

1. If n=0, output 0. Compute L and the parameters above. If any finite
   guard in (5'') or its exact supplied parameter check fails, output 0.
2. Build the exact base extractor E_n from (L), seed-padded to s bits.
   Its data argument is x and its output length is d.
3. Build the (G') graph on left pL-bit strings: the four rows of the AB+C
   step, the short-start Cayley walk and the s-bit prefix (sampler-and-
   disperser.md Section 4). Its degree is exactly C.
4. Build the canonical field of degree L and the tags in Section D1.
   For each alpha in lexicographic {0,1}^{pL} and j=1,...,C compute

     y_(alpha,j)=E_n(x,Gamma(alpha,j)),
     b_(alpha,j)=F_n(x,Prefix(y_(alpha,j),d_F),
                       alpha || encode(j-1) || H_alpha(x)),

   using the ONE breaker parameter setting (O2). Let R_alpha be the XOR
   of these C bits.
5. Output 1 iff strictly more than half of the A row bits are one;
   otherwise output zero. This is exactly Section D1's evaluator.

F_n in this pseudocode denotes the raw d_F-seed breaker; the preceding
proof uses its padded d-seed form. The two evaluations are identical.
No arbitrary seed, affine basis, irreducible polynomial advice, entropy
witness, good row, tampering, or source decomposition is supplied to this
algorithm. Constants, field bases, finite prime searches, orderings and
zero branches are fixed by the component algorithms.

The base sampler tests density epsilon0=2^-L. By (23'), for every test of
density at most epsilon0 fewer than Lambda_S rows are wholly bad, except on
source probability 1/100. Source
linearity survives graph index selection. The degree and advice lengths match Section D1, and (D1) was checked
above. Sections D2--D4 prove the affine-permutation property at k=HL,
and in particular the directional property.

## 6. Uniform polynomial time in actual input length

All one-time constant computations terminate, independent of n. Their
finite cost is part of the fixed machine constant, not length-dependent
advice. Parameter arithmetic has O(log n) bits beyond these fixed constant
lengths. There are A*C <= C*[2(n+1)]^p breaker calls, a fixed polynomial.
The base linear extractor is polynomial in n at fixed h. The graph's
prime p_m has m=pL/3 bits and is found in 2^O(pL) bit operations, so its
evaluator over all A rows is also polynomial in n. The walk fields have fixed-algorithm arithmetic.

By (F'), each breaker call is polynomial in n for the fixed t=4C-1 and
a,E=O(log n). Its local prime searches and all finite field setup costs
2^O(log n) with fixed coefficients. The tag field search and arithmetic cost 2^O(L), and all A polynomial
evaluations cost a fixed polynomial in n. Counting the majority uses an O(log A)
bit counter. Therefore the whole fixed total evaluator is polynomial in
the ACTUAL n. The polynomial degree and multiplier can be extremely large
but do not depend on n. No enumeration of extractor truth tables or
verification over all source distributions is hidden in this cost.

## 7. Native circuit bound, NP eligibility and expression

The assembly is nonconstant under every fixed-point-free affine permutation
on a source flat of dimension at least k=HL. In particular, for every
nonzero v and every affine k-flat S its translation derivative is
nonconstant. We prove a finite bridge under the WEAKER hypothesis

  for every nonzero v and affine k-flat S, some x in S has f(x)!=f(x+v).

For integers k>=1 this implies C_fullB2(f)>=4n-4k. For n<=k this is trivial.
Otherwise use the retained affine-map, three-corner-pair bridge with d=k-1,
stopping at a=M-q=k. Separation implies nonconstancy on every flat of
dimension at least k, by choosing a nonzero tangent translation of any
hypothetically constant flat.

At every state with a>=k an unused parameter is impossible: for a free
parameter keep its coordinate direction in a safe a-flat; for a paired
parameter fix its partner safely. Independence gives zero derivative
there, and injective Phi preserves its nonzero direction. Absence includes
the output reference. Pair rule1 is therefore impossible; rule3 creates an
unused free line coordinate at the same a and is impossible as well.
Mixed and paired/paired productive actions create the unused coordinate
identified in the retained ledgers at new a>=k, so only FF actions occur.

Put W=s-2q. Each FF action removes at least two objects and adds one pair,
dropping W by at least 4 as a drops by one. Basic deletions, transfers,
role exchanges and one-child-row elimination do not increase W. Each
remaining pair rule2 or4 removes at least two objects and one pair, also
not increasing W. Termination thus gives exactly n-k FF actions and

  W_T<=W_0-4(n-k)<=|C|-4(n-k).

All parameters lie in the output's undirected incidence component;
otherwise the invertible auxiliary equations split by components and
the output is independent of a parameter, already excluded. In each
component the row/auxiliary block is square: full row and column rank
of invertible A force both inequalities for that block. Each gate-output
variable belongs with its defining object. Thus the output component has
M+s_0 variable vertices and s_0 hyperedges of size at most three. A connected
hypergraph of s_0 such edges has at most 2s_0+1 vertices, since adjoining
an intersecting edge adds at most two. Hence s>=s_0>=M-1. At a=k,

  W_T=s-2q>=k-q-1>=-1.

The only integer case with W_T<0 is q=k, M=2k, s=M-1. The tight incidence
tree proof and paired read-once zero-period lemma, proved in
sampler-and-disperser.md, exclude exactly this case under our weaker
separation hypothesis. Thus W_T>=0, proving the finite bound. Its uniform
additive constant zero is optimal for this weaker hypothesis: n=k=1,
f(x)=x satisfies separation and has zero gates. This does not assert
optimality under the stronger nonconstant-derivative hypothesis.

For the fixed evaluator C_fullB2(f_n)>=4n-4H*lg(n+1) at every integer actual
n>=N, the fixed threshold of Section 4. For n>=1,
lg(n+1)<=floor(log2 n)+1<=log2 n+1, so with K=K0=4H

  C_fullB2(f_n)>=4n-K*log2(n)-K0   for every integer n>=N,           (O4)
  K=K0=7377880453983534785844750636.

K0 cannot be lowered with K=4H from this kernel: at n=2^i, lg(n+1)=log2 n+1.
Here K0 is the declared additive constant 4H; the base-extractor threshold
K_0 of (L) and the L0 threshold K0 of affine-and-standard-wrappers.md Section
2.1 are different quantities.
All constants and zero guard branches belong to one fixed total P
evaluator. Set f_0=0. The NP verifier accepts (x,empty) exactly when this
evaluator outputs one and rejects every nonempty witness, with witness
bound zero. Tags, seeds and field descriptions are internal, not circuit
inputs. The common threshold holds at every integer length. The declared
explicit remainder is exactly K*log2(n)+K0 with these 28-digit integers, and
log2 n/n tends to zero, giving (4,1).
