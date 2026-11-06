# Fixed logarithmic-loss family and native full-B2 theorem

This file and its four companion proofs establish the theorem declared in
bound.json: H=7412884525208893, K=29651538100835572. Mathematical arguments are
supplied below.

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

    C_fullB2(f) >= 4n - 8(d+1). (1)

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

    c + (x+b)(y+e), b,e,c in F_2. (2)

Call the latter product gates. These exhaust the basis.

Use M parameter coordinates u, of which 2q are in disjoint pairs. Each pair
ranges over three of its four corners, and every unpaired (free) coordinate
ranges over F_2. Write D for this product domain. Maintain an affine map
Phi:F_2^M -> F_2^n injective on the entire cube, and a representation of
f(Phi(u)) for every u in D. Equality outside D is not required.

Auxiliaries v satisfy a square system

    A v = N u + c, A invertible, (3)

with at most three distinct nonconstant variables in each row, after
cancellation. A row has an auxiliary, since A is invertible. Remaining
binary gates read parameters, auxiliaries or earlier gates. Gate outputs
never enter (3). Let s be rows plus gates. A parameter's children are the
distinct rows and gates containing it, not pins and not the output reference.
Let I count free parameters having a child. Set

    a=M-q, V=s+4I+2q. (4)

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
    Delta a=-1, Delta V<=-(1+k)-8+2<=-8. (5)

**Paired p and free x.** Write h for p's partner and B for h's sole gate
child; B!=G. Safely fix h and fix x to kill one factor in (2) at G. G is
constant, B disappears, and newly free p loses its sole use. If G has a
child other than B, it is a third distinct deletion. Otherwise B is G's
sole child and reads the two constants h,G, so B is constant too. It cannot
be the output, since the restricted output is nonconstant. Its old child
exists by cleanup, differs from B by acyclicity and from G because G reads
only parameters, and is a third distinct deletion. Hence

    Delta M=-2, Delta q=-1, Delta I<=-1, Delta s<=-3,
    Delta a=-1, Delta V<=-3-4-2=-9. (6)

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
    Delta a=-1, Delta V<=-4-4=-8. (7)

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

# Exact finite direct-inner and affine-outer assembly

The dependency order is noncircular:

1. Finite seeded-linear extraction PXT7g and exact fields/towers, together
   with the local alternating-extraction operations of restricted-nonmalleability-and-breaking.md
   local operations/conditional tools inSections1--3 andthefive-roundtheorem8.2.
2. The five-round antichain14-copy joint theorem in restricted Section8.2, instantiated as
   NB's innerF# by affine-and-standard-wrappers.md4.2.
3. The CZ/weak-AGHP advice condenser, the AGHP extractor, Ext0, Ext2 and F#: the
   standard breaker NB (affine-and-standard-wrappers.md Section 4).
4. The same transformation at the sampler's degree with C=NB and the slots
   A# of outer Section 3 (NB's seed and threshold by Lemma B9 of Section 5):
   the outer generalized affine breaker.

Constants use exact integer/rational operations or proved terminating searches,
never tests of source distributions or extractor properties. For positive
rational x, lg(x) is the least nonnegative integer j with2^j>=x.
Powers in error bounds are analytical, not floating-point computations.

## 1. Inputs and their precise scope

All selected inputs are supplied in this package. Linear-and-bilinear-primitives.md
contains seeded-linear, bilinear-extraction and weak-seed-hashing (including
the weak-seed lift); Section7 gives Theorem D3, Corollary D3' and Lemma P,
with j-way splits, fixed R, padding U1 and admissible dimensions. The base
uses7f's directly proved rational construction. Output truncation and seed
padding retain the exact K recurrence, seed lengths and all suffixes.
Restricted8.2 supplies constant-advice-breaker8.2's five-round14-copy
independent-side contract, retaining every seed, used at NB4.6 step7.
Affine-and-standard-wrappers.md Sections1--2 supply the affine transformation;
Section4 supplies Theorem NB in place of the earlier wrappers.

CZ Section8 in restricted-nonmalleability-and-breaking.md uses
eta=268435456/1975781445633. Its ring input is sampler-and-disperser.md
SP.5, gamma=256/20481, prime threshold2^3200000. The CZ source/prime guards
imply that threshold; prime-fields-and-encoding is included in restricted1--3.
The actual four-row source-to-incidence injection is sampler SP.6.2;
SP.7 supplies the full multiplicity/incidence proof and characteristic guards.
Thus the exact ring contract closes CZ's conditional premise, without a
scalar sum-product or hitting-peeling assumption.

No normalized (F) contract from recursive-sampler-schedule.md or external
degree envelope is assumed. Outer Section3's A# fixes d_O,k_O; Section5's
Lemma B9 bounds NB. Om is a closed form.

## 2. Primitive convention

Selected slots use PXT7g for NB Ext0/Ext2 and PXG7g for inner/outer maps;
explicitly named AGHP calls use the separately proved matrix extractor. Each has
its explicit source, seed and dimension guards. No older unselected asymptotic/envelope contract is used.

## 3. The affine assemblies and the standard-breaker interface

### 3.1 The standard interface required here

A standard one-tampering breaker C is used in the outer assembly through
its actual seed lengthd_C and thresholdk_C: arbitrary correlated source
pair independent of an arbitrary correlated seed pair, honestseed uniform,
distinct advice,bothcomplete seeds retained,error<=2^-e. HereC=NB has
explicitguardsk_C<=u,d2<=u and works at every advicelengtha>=1.
NB's innerF# instead has the five-round antichain14-copy tuple contract of wrappers4.2
and restricted8.2 at eight fixed labels; NB4.6 step7 uses exactly that
zero-shift contract. The outer assembly retains arbitrary additive shifts.

### 3.2 The selected outer affine assembly

The outerA# is wrappers2.1 withC=NB and the explicit outerSection3 slots.
F# is the five-round antichain multicopy map of wrappers4.2. All parameters are
closed forms;failedpublicguards return the all-zero string and no claim
is made at those exceptional calls. The total evaluator supplies all such
branches, and every guard holds at every sufficiently large actual length.

## 4. Fix the uniform-source advice condenser first

Fixgamma=256/20481,eta=gamma/(92+2^-20)=268435456/1975781445633,theta=eta/7=268435456/13830470119431.
The selected energy proof absorbs its fixed2^64 factor using the explicit
guardeta*logp>=64*2^20. Every actual NB call haseta*r>=112*2^20,as proved
in wrappers(I3),so this guard and the independent ring-prime guard hold.
The full supplied ring/energy proofs establish conditionalentropy19eta/6.
NB's repaired honest Ext0 output is exactly uniform, so it uses Lemma U
of affine-and-standard-wrappers.md Section3.5. It does not use the
general-source capacity/deficiency choice theta=eta/200.
Put sig=Ee+4 and s=floor((eta*r-17sig-8)/7). NB's finite parameter
algorithm proves all Lemma U guards, s>=m0+2e_L and s>=a+40.
Its advice row error is at most14*2^-sig<2^-Ee. Advice is the actual
a-bit label, not a silently enlarged source deficiency.

## 5. The standard breaker NB and Lemma B9

Theorem NB (affine-and-standard-wrappers.md Section 4) is a standard breaker
in the sense of Section 3.1, with guards k_C<=n and d2<=n, at every a,z,E>=1.
In the dependency order above it uses the direct multicopy F#, the advice condenser with
the constants of Section 4, the AGHP extractor and the primitives Ext0, Ext2.
The outer assembly uses NB only through the outer A# (outer Section 3), and
Lemma B9 gives NB's seed and threshold at the outer call.

**Lemma B9 (NB at the outer call).** Source length n'=d_O<=n,
outputz=120Om,advicea=s_tag+2+lgC<=Om,exponente,Om=e+2L+42.
Then d_C<15528Om,k_C<878485822Om. The outer ledger supplies k_C<=n'
and d2<=d_C<=sC<=n'. Put eps=10^-6; Om>=2^40,
8lgOm+240<=epsOm, by its k=40 endpoint and doubling.
Ee=e+4 and w_N=max(2^20,Ee+44+2lg(4n'))<=Om+12<=(1+eps)Om.
PXT heights 21,6 are below38. F#'s unit
w_F=max(2^20,Ee+2lg(4n')+32,ceil(z/128))<=Om since
n'<2^(L+1), z=120Om<=128Om. Its twelve charged rows give
m1<=dF*Om,dF=50163402688259/16384000000,k_F<=kF*Om,kF=408315228965111423/2048000000.
These bounds include eps for the final slot ceiling. Output-dependent
units and all seed suffixes remain part of the algorithm.
Ext2=PXT(vNB2): c=51013/2624832,pi=35547607831/410130000,lam=17767192194518570813999/9746610401172264448374,g=2070240070850580720873907300087/1327975667159721031090957500.
M0=c2*dF+pi2*(1+eps)+eps=3929016421673820221791/26878279680000000000,
K2bar=Lam2*dF+g2*(1+eps)+eps=77676703017926309621768070123836563758329/10878776665372434686697123840000000000.
Ext0=PXT(vNB0): c=178499/332396494848,pi=21194095677774558751799/35005298115870720000,lam=519667411934432347261058277919896373370713568062639645846889719179264/229450303876896571514656611690229888519369238617337648096710545579689,g=630807802392579255480360185075574373081703258406567294485991735041393207381934/5736257596922414287866415292255747212984230965433441202417763639492225.
Write its coefficients c0,pi0,Lam0,g0 and theta=268435456/13830470119431>7c0.
Since (38/7)Ee+a+53<=7Om, the actual NB algorithm gives
Rnew=(M0+pi0*(1+eps)+7+7eps)/(theta-7c0)+eps=142832751296224170798321244199355508760447/2946518045565944935120557812500000,
d0bar=7c0*Rnew+pi0*(1+eps)+eps=220416610059112382413132012556161675845116278861/279832077243571327106014911194726400000000000.
Thus r<=Rnew*Om,d0<=d0bar*Om. BOTH source branches satisfy
7Lam0*Rnew+g0*(1+eps)<878485822,
kF+14Rnew+130+eps<878485822; here z+2Ee<=130Om.
The AGHP call has2e_L<=d0+3Ee. Its logarithmic term satisfies
2lgOm+lg(ceil(Rnew*M0)+1)+1<Om/4: the fixed constant is <240,
and the k=40 endpoint and doubling prove the all-Om claim. Hence
d1<=[2(M0+d0bar+3)+1]Om=1195307250130881564319189254088833775739392132057087/637597387999477268811054975157184102400000000000Om,
d_C=2d0+2d1+m1+K2+Ee<=897594492524117129061773915931651543872298481870768102505997116170560467/57806029473145679053277134586964211417245445679808646348800000000000Om<15528Om.
These exact inequalities refer to the actual NB algorithm, including
both source branches, complete retained seeds and every field guard.
All base fields have degree O(log actual n); larger extensions use7g's towers.

## 6. Uniformity, exactness, and remaining obligations

F#/A# slots and NB parameters are closed formulas. The actual primitive
recurrence computes seed lengths and thresholds; envelopes choose padding.
The evaluator receives no proof witnesses, coupling coins or truth tables.

For fixed t and m,a,E=O(log n), all auxiliary dimensions are O(log n).
Closed slots and LemmaB9 give F#, Ext0 and Ext2 logarithmic enumerated
field degrees, so setup costs2^O(logactualn). CZ/AGHP/bilinear setup is
2^O(local dimension), making the composition polynomial in actual n.
Its exponent may depend on the fixed t=4C-1; no joint bound in t,n is claimed.

For ordinary seeded primitives or standard breakers with the complete
source tuple independent of the complete seed tuple, unused suffixes may be
retained by their conditional kernel given all used prefixes. The outputs
depend only on the source tuple and those prefixes. This applies to the
primitive pads of F# and the outer A# (restricted-nonmalleability-and-breaking.md
Section 8.2, affine-and-standard-wrappers.md Section 4.2, outer Section 3) and to the
zero-additive-shift independent-side use of NB's inner breaker
(affine-and-standard-wrappers.md Section 4.6, step 7). It does not assert a generic
kernel principle when arbitrary additive B-tuples correlate with suffixes.

# Tagged affine-permutation logarithmic-loss family: outer assembly

All constants below are fixed once; none is length-dependent advice.
The affine-breaker contract (F') is proved below from the outer A# (Section 3) and Lemma
B9.

Credits: sumset framework, Chattopadhyay--Liao,
https://arxiv.org/abs/2110.12652v1 ; advice-condenser architecture, Xin Li,
https://eccc.weizmann.ac.il/report/2023/023/revision/1/ ; native full-B2 bridge,
minimal-circuit-bridge.md. All selected arguments are supplied here.

All logarithms are binary. Write lg(x)=ceil(log2 x) for positive x, computed
by integer bit length and a power-of-two comparison. All string lengths,
indices and loop counts are integers. All arithmetic parameter algorithms
use exact integers and rational comparisons, not floating point.

Foractualintegern>=1,flg(n)=floor(log2n),computedasbitlength(n)-1.
ONLYtheactualouterlengthusesflg.Everylgforfielddegrees,seeds,graphconstants
remainsCEILINGlog2.Small-lengthzerobranchesremainexplicit.

## 1. The exact inputs to this assembly

Fix Q_S=293/2 before input length. Every sampler/base unit is the integer
S=floor(L/Q_S)=floor(2L/293),L=floor(log2(actual input length)).
All field/seed dimensions use integer ceilings. Even for rational Q_S,
S<=L/Q_S<S+1 and L<Q_S(S+1).

(L) PutactualL=flg(n),S=floor(actualL/Q_S)>=2^40. Base7f uses
h=19,b=14,alpha=114792754196981/549755289600,Ro=817286867215042604416161689/110012249584831199340134400,initialrate2^40,
s=2ceil((2378/7)S/2),r_F=S+8lg(s+1), and error <2^(-r_F).
Its complete seed is padded to s; its source threshold is K0+7,
ceilkfac=1085568053810387218. Output4^19ceil(alphaS)>=delta_O*actualL,
then truncate to delta_O*actualL. Sampler5.3 proves these exact claims.

(G') TheoremG' ofsampler-and-disperser.md Section4 usesM=18912/7,
p=3M=56736/7,q=1181,C=4728,s=2ceil(B0S/2),r_F=S+8lg(s+1),s_tag=3(s+2q*r_F),A=2^s_tag,Lambda_S=A^(1/2)/2^15.
Itisa deterministicgraphfroms_tag-bitrowstos-bitneighbors,degreeC.
ForS>=2^40everyrowsetofsize>=Lambda_Shasmorethan2^(1-r_F)2^sneighbors.
ItsfourAB+Cimagesfollowedbytheexactshort-startwalkcost2^O(s_tag).
Withbase(L),sampler(23')atthresholdK0+7hasfailure1/100andlessthan
Lambda_Sall-badrowsfortestdensityepsilon0=2^(-r_F).

(F') For every n passing the guards (5'') of outer Section 4, the outer A#
(outer Section 3) with C=NB, at the request (O2) and the unit Om, returns an
exact seed length d_F, source entropy requirement k_F and a one-bit
generalized affine breaker, where

  d_F=d_O<=delta_O*L, k_F=k_O<=kappa_O*L. (F')

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

Write L=flg(n),S=floor(L/Q_S),s=2ceil((2378/7)S/2),r_F=S+8lg(s+1),m_walk=s+2q*r_F,s_tag=3m_walk=log2 A,w=2. Require S>=2^40,L>=2048 and s_tag>=2L; the latter follows from the fixed constants below.
Build the canonical field F_(2^L), and split x into L-bit coefficients of
P_x, padding its final block with zeros. For a full s_tag-bit row alpha,
let beta,gamma be its first and second L-bit blocks. Define
H_alpha(x)=Prefix_2(gamma*P_x(beta)), a two-bit linear map in x.
Each (beta,gamma) has exactly A/2^(2L) row preimages. For Delta!=0,
P_Delta is nonzero of degree ceil(n/L)-1<n/L, so its root fraction is
at most n/(L*2^L)<2/L<=1/1024. Heren<2^(L+1),L>=2048;
thedegreeisalso<2^L,sotherootboundapplies. Off those roots, multiplication by P_Delta(beta)
permutes gamma and exactly 1/4 of gamma values have zero two-bit
prefix. Thus at most A*(2/L+1/4)<=257A/1024 rows have H_alpha(Delta)=0.
The root bound follows by factoring X-beta at distinct roots. The supplied
finite-field proof makes the canonical search halt in 2^O(L) time.
Full alpha remains in the advice label despite repeated (beta,gamma).
No independence of hash values across different rows is asserted.

Use the existing source-linear sampler Samp(x,alpha,j), j in [C], seed
length d, tolerating test density epsilon0=2^(-r_F) with fewer than Lambda_S
all-bad rows ((23') in (G')), failure probability 1/100, at source
entropy k_S. The arbitrary enumerated field index alpha replaces [A].
For t=4C-1 and a=s_tag+w+lg(C), request the padded all-seeds-retained affine
breaker F of output length one, entropy K, and error at most

    epsilon = 2^(-r_F)/(2^33 C A^6).

The advice for (x,alpha,j) is the concatenation of the s_tag-bit alpha,
lg(C)-bit j-1, and w-bit H_alpha(x). Define

    R_alpha(x)=XOR_(j=1)^C F(x,Samp(x,alpha,j),advice(x,alpha,j)),
    D(x)=1 iff 2 sum_alpha R_alpha(x)>A.

Assume the input source entropy satisfies

    k >= max(k_S, K+4w+4Cd+2lg(1/epsilon)). (D1)

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

  k>=max(k_S,K+4w+4Cd+2lg(1/epsilon)). (D1')

Choose theta=1/(1024A^2), epsilon=2^(-r_F)/(2^33*C A^6), and let B be the
union of their bad-seed sets. As 64*2^16*2*1024=2^33, Section D2 gives

  |B|/2^d<=64C A^4*2^(8w)*2epsilon/theta=2^33 C A^6 epsilon=2^(-r_F). (D2')

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

  eta+Pr[at least one selected alpha collides]. (D3a)

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
  E S_0^2 S_1^2=g^2-2g+2E Z^2+e22. (D3')

Marginal fourth moments have only single-copy surviving parities, so
E S_sigma^4<=3g^2+g^4 eta<=4g^2. No collision error is needed there.
The elementary variance identity gives E Z^2>=(E Z)^2=zeta^2.

Define Q_+(x,y)=xy(x-y)^2 and Q_-(x,y)=-xy(x+y)^2.
Both have central term -2x^2y^2. Set r=|zeta|/g. Since
N_col/g<=(257/1024)*(256/255)=257/1020 and eta<=3/(1024A^2),
we have 0<=r<13/50. Both quartic expectations are at most

  2(3g-2)|zeta|-2g^2+4g-4zeta^2+4g^4 eta
  <=[6r-4r^2-2+4/g+12/1024]g^2<=-g^2/2. (D4')

For the final inequality, 6r-4r^2 increases on [0,13/50],
and 6r-4r^2-2 equals -444/625 at13/50; also g>=A/2>=2^(2L-1),L>=2048 gives 4/g<=2^(3-2L)<1/1000.
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
    <=432delta g^(3/2)<(432/10000)g^2<g^2/8. (D5')

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

The base(L) is fixed7f with h=19,b=14,alpha=114792754196981/549755289600,Ro=817286867215042604416161689/110012249584831199340134400,
initialrate2^40,B0=2378/7,M=18912/7,p=56736/7,q=1181,C=4728,t=18911.
For actual n>=1,L=flg(n),S=floor(L/Q_S), set
s=2ceil(B0S/2),delta_v=s-B0S in[0,2),ellc=lg(s+1),r_F=S+8ellc,
m_walk=s+2q*r_F=MS+16qellc+delta_v,
s_tag=3m_walk=pS+48qellc+3delta_v,A=2^s_tag,d=delta_O*L. Request
 t=4C-1,m=1,a=s_tag+2+lgC,w=2,E=6s_tag+r_F+33+lgC. (O2)
Here lgC=13,lg(26t)=19,j=15,
e=E+lg(26t)=(6p+1)S+(288q+8)ellc+18delta_v+65.
Then2^-E<=2^(-r_F)/(2^33*C*A^6). All tags and seeds are internal.

(O2a) Pay the logarithmic terms jointly. On2^k<=S<2^(k+1),k>=40,
B0<511 ands<B0S+2 giveellc<=k+10; alsoL>=Q_S*S. Atk40,
2(288q+8)*50+172<2^-20*Q_S*2^40.
The left increases by2(288q+8), while the right doubles. Thus
2(288q+8)ellc+max(65+42+36,74+2lgC+72)<=2^-20*L.
This pays18delta_v<36 inOm and36delta_v<72 in native entropy.
For base7f the k40 comparisons are
8*50+64=464<1048576,10*40+468=868<1048576,
Q_S+234+9*40=1481/2<1048576=2^-20*2^40.
Their left sides grow by at most10 per k, and their right sides double.
Sampler5.3 derives these bounds on a-S,rho-S,rho0-(Q_S+1)S
from the integer definitions, including lgS<=k+1, for every S.

## 3. The outer affine breaker A# and (F')

For actual n passing(5''), Om=e+2L+42<=c_O*L,
c_O=(6p+1)/Q_S+2+2^-20=718220036099/2150629376. (O2a) pays its entire overhead.
Set q=e,j=15, and integer-multiple slots
 dy=28Om,dx=120Om,r=529567Om,sC=15528Om,s0=53Om,
 kC=878485822Om,d_O=1173153695Om,k_O=10368947086Om.

**Lemma A#.** Run wrappers2.1 at t=18911,C=NB and (O2)'s advice,
output and error, padding NB's seed to sC. Use PXG(vL0),PXG(vLr),PXG(vE)
for L0,Lr,E at exponent e. Lm is the supplied strong uniform-seed AGHP
at input r, output dy, exponent e, padded to dx. Then(F') holds with
 delta_O=ceil(1173153695c_O)=391784144017,
 kappa_O=ceil(10368947086c_O)=3462793558725.
Inputs n,d_O,r are <=n since d_O<=delta_O*L<=H*L<=n and r<=k_O.
For u<=n<2^(L+1),lg(4u)<=L+3; heights<=15
give e+height+6+2lg(4u)<=e+2L+37<=Om. Om>=2^40.
All PXG terminal rates lie in[5/3,64] and outputs sC,r,dx in unit Om
are <=256*2^height, respectively. Their full field/seed guards follow7g.
NB's unit is legal since 120<=128; B9 supplies all inner requirements.
For Lm,lg(r*dy+1)<=33+2lgOm and66+4lgOm<=10^-6*Om
by the k=40 endpoint and doubling. Its seed length is
<=(2*28+4+10^-6)Om<=dx and its threshold<=(28+2)Om.
It is strong source-linear, retaining its ENTIRE seed and suffix;
its O(log actual n)-degree canonical field is uniformly computable.

All outer slots are exact integer multiples of Om. The causal ledger is
TR=(t+1)(s0+sC+2j(dy+dx)),TL=(t+1)(sC+2jdy).
In units Om these are TR=378637152,TL=309551616.
The fourteen finite requirements, with all primitive ceilings paid, are:
 dy_E: [91472297/3276800]Om <= 28Om
 dx_Lr: [78154067393/655360000]Om <= 120Om
 dx_Lm: [60000001/1000000]Om <= 120Om
 s0_L0: [429200549/8192000]Om <= 53Om
 r_3a: [529567]Om <= 529567Om
 r_dy: [28]Om <= 529567Om
 r_Km: [30]Om <= 529567Om
 d_3c: [1173153695]Om <= 1173153695Om
 d_3d: [302308258394307521/793496200]Om <= 1173153695Om
 d_s0: [53]Om <= 1173153695Om
 d_KE: [60419513404921/793496200]Om <= 1173153695Om
 d_kC: [878485822]Om <= 1173153695Om
 k_3b: [6122804222559220666117678248032344268544820774852995336361467747715816011/590494306916503663291892073888697082992058534720787733166530000]Om <= 10368947086Om
 k_K0: [1000567152964649786888363/749663429184000]Om <= 10368947086Om
Exact primitive coefficients:
 E: c=1/32,pi=79184297/3276800,lam=9248365971/63479696,g=46546964448421/793496200
 Lr: c=1/16384,pi=56971387393/655360000,lam=3288142835770757990472143911420012578108615980093077945813653/59049430691650366329189207388869708299205853472078773316653,g=9014022244772565976201376736406068328927759934067810031092991997686011/590494306916503663291892073888697082992058534720787733166530000
 L0: c=1/8192,pi=413672549/8192000,lam=555159733986554411/47978459467776,g=171256234672895969723/148272039000
The3c row is(t+1)(s0+sC)+kC+q;3d includes TR+K_E+t*dx+q;
3b includes TL+K_r+t*r+q. Also k_O>=d_O,r>=1. All fourteen rows
are for the actual seed maps; full seed/source/padding guards and BOTH
NB entropy branches remain charged. Error(13*2^j-8)2^-e<26t2^-e<=2^-E.
Parameters are fixed before source data; there is no entropy witness,
generator search or per-length selector in the evaluator.

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

  H=max(4Cdelta_O+kappa_O+T,ceil(kfac/Q_S)+7),
  T=ceil((12p+2)/Q_S+2^-20)=664. (O3)
Rawbranch=7412884525208893,basebranch=7410020845122105,H=7412884525208893.
SinceS<=L/Q_Sand(O2a),
4w+2lg(1/epsilon)=(12p+2)S+2(288q+8)ellc+36delta_v+74+2lgC<=TL.
Hencek_F+4w+4Cd+2lg(1/epsilon)<=[kappa_O+4Cdelta_O+T]L<=HL.
AlsoK0+7<=kfacS+7<=(ceil(kfac/Q_S)+7)L<=HL,proving(D1).
Sincep>3*Q_S,S>=L/Q_S-1,L>=p+128,
pS>3L-p>=2L ands_tag>=pS. Therefore each ordered pair of L-bit tag blocks has
exactly A/2^(2L) preimages among the complete row indices. AlsoA>=2^(2L).
The collision proof and all bad-row counts in D1--D4 use this exact
coverage and A itself; they require no stronger lower power of n.
FlooractualL/bitcountsandall-lengthguardsremainexplicit.

The evaluator tests these literal finite guards:

  n>=1; L=flg(n)>=max(2048,p+128); S=floor(L/Q_S)>=2^40;
  s_tag>=2L; HL<=n. (5'')

The sampler-unit guard implies all lower bounds in7f,sampler3--5,
and all unit rounding inequalities above. The evaluator also computes
and checks each supplied algorithm's exact parameters and guards.
The proved(L),(G'),(F') show they all pass under(5''); no source
property, good-row witness, extractor oracle or existential search is
an evaluator input. S tends to infinity at EVERY integer actual length.
FinallyHL<=neventually:L=floor(log2n),n>=2^L,and2^L/L
tendstoinfinitywithconsecutiveratio>=4/3forL>=2.
Thus one fixed eventual threshold makes all guards hold for all n.

These guards hold at all sufficiently large integer lengths; smaller lengths use zero.

## 5. The fixed evaluator

Precompute the constants in Section 2 by their finite algorithms. On an
input x of actual length n:

1. If n=0, output 0. Compute L,S and the parameters above. If any finite
   guard in (5'') or its exact supplied parameter check fails, output 0.
2. Build the exact base extractor E_n from (L), seed-padded to s bits.
   Its data argument is x and its output length is d.
3. Build the (G') graph on left s_tag-bit strings: the four rows of the AB+C
   step, the short-start Cayley walk and the s-bit prefix (sampler-and-
   disperser.md Section 4). Its degree is exactly C.
4. Build the canonical field of degree L and the tags in Section D1.
   For each alpha in lexicographic {0,1}^{s_tag} and j=1,...,C compute

     y_(alpha,j)=E_n(x,Gamma(alpha,j)),
     b_(alpha,j)=F_n(x,Prefix(y_(alpha,j),d_F),
                       alpha || encode(j-1) || H_alpha(x)),

   using the ONE breaker parameter setting (O2). Let R_alpha be the XOR
   of these C bits.
5. Output 1 iff strictly more than half of the A row bits are one;
   otherwise output zero. This is exactly Section D1's evaluator.

F_n denotes the raw d_F-seed breaker, identical to its padded d-seed form.
All seeds, fields, primes, orderings and zero branches are internal fixed
component algorithms; no advice or source witness is supplied.

The base sampler tests density epsilon0=2^(-r_F). By (23'), for every test of
density at most epsilon0 fewer than Lambda_S rows are wholly bad, except on
source probability 1/100. Source
linearity survives graph index selection. The degree and advice lengths match Section D1, and (D1) was checked
above. Sections D2--D4 prove the affine-permutation property at k=HL,
and in particular the directional property.

## 6. Uniform polynomial time in actual input length

Constant computation terminates independently of n; parameter arithmetic
adds O(log n) bits. There are A*C <= C*[2(n+1)]^p breaker calls.
Base field degrees O(logn) cost2^O(logn) to enumerate; towers build the
O(n) initial relative extension without2^O(n) enumeration. The graph prime
p_m at m=s_tag/3 is found in2^O(s_tag) bit operations. All A graph rows
and fixed-algorithm walk arithmetic therefore take polynomial time.

By(F'), each breaker call is polynomial in n for fixed t=4C-1,a,E=O(log n).
Ext0 input/output/error and7g(PXT) enumerated degrees are O(log n):
setup and local prime searches cost2^O(log n), with towers for larger
extensions. Tags cost2^O(L); all A evaluations are polynomial in n.
Majority uses an O(log A)-bit counter. The complete total evaluator is
polynomial in ACTUAL n with fixed degree and multiplier, without truth-table
or source-distribution searches.

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

For the fixed evaluator the finite native theorem gives
 C_fullB2(f_n)>=4n-4H*flg(n) at EVERYactualinteger n>=N,
whereNisthefixedcommonguardthreshold.Sinceflg(n)<=log2n,
 C_fullB2(f_n)>=4n-Klog2n,K=4H=29651538100835572. (O4)
The additive remainder constant isZERO. Base reverse thresholdK_0 and
primitiveK0aredifferentobjects;neitherdenotesasubtractedconstantin(O4).
AllconstantsandzeroguardsdefineonefixedtotalPevaluator.Setf_0=0.
TheNPverifieraccepts(x,empty)exactlywhenthisevaluatoroutputs1andrejects
every nonemptywitness,withwitnessboundzero.Tags,seeds,fielddescriptions
areinternal,notcircuitinputs.ThecommonthresholdholdsatEVERYintegerlength.
ThedeclaredremainderisexactlyKlog2n,withnoadditiveloss.Log2n/ntendsto
zero,sobothjointclaimsestablishprojection(4,1).Theevaluatorreallyuses
flooractualLratherthantheformerceil-lengthconstruction.Nofavorable
subsequence,per-lengthadvice,smoothingwitnessorpriorverdictisused.

The actual stretched row length is polynomially bounded as follows. AtS>=2^40,
48q*lg(s+1)+6<=pS by the k40 comparison48q*50+6<p*2^40 and doubling.
Thus s_tag=pS+48qellc+3delta_v<=2pS<=pL sinceQ_S>=2 andL>=Q_SS. Also r_F=O(S).
Every enumeration previously bounded by2^O(pS) is now2^O(s_tag)<=2^O(pL),
a fixed polynomial in actualn. Allfielddegrees, parsinglengths and advice
labels use actuals_tag andm_walk, neverjusttheirleadingcoefficients.
