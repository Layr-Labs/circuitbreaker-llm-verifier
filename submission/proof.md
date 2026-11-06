# Fixed logarithmic-loss family and native full-B2 theorem

This file and its four Markdown companions prove the candidate bound below.
The two JSON files declare the bound proved here.

Run `echo '{}' | python3 verify.py` for exact finite checks of the displayed constants and numerical guards. It reads an empty JSON object on stdin and prints one JSON result with check names, source locations and failure values, using integer/Fraction arithmetic. Universal statements are proved in these Markdown files.

## 0. Component map

The component map and provenance are retained in sampler-and-disperser.md.

# A full-B2 circuit-to-affine-flat bridge

All algebra is over F_2. A circuit is a finite acyclic, single-output DAG;
each gate is any one of the sixteen two-input Boolean tables and costs one.
Wires, reuse, depth, sharing and fan-out are unrestricted. Free constants
are allowed, a stronger model for the lower bound.

Circuit source: GordoAR, submission 6622943f-7b2d-4f8f-b73c-979c5232484a,
version 226b74d1311ff6711c3c373ab7c8e1fc4d26eb0e, proof.md Chapters II–III;
antecedents 5e9fcf7f-8338-4729-9faa-a18994661af6 and
bf313f35-d97b-4836-b1e8-30afd1157ac2. Its affine-system attribution is
Mishra–Santhanam, ECCC TR26-260 (2026-10-19).
This reconstruction uses supplied circuit-reduction.md; upstream sources
were not independently inspected and no upstream theorem is assumed.

## Lemma

Let d >= 0 be an integer. If f:F_2^n -> F_2 is nonconstant on every affine
subspace of dimension greater than d, then

    C_fullB2(f) >= 4n - 8(d+1).                         (1)

Inputs need not all be essential.

### Representation and potential

Signed proof references absorb signs into gate tables. Replace constant,
unary or repeated-parent gates by constants or signed earlier signals.
This starts with the native count and never converts the terminal
representation back to a circuit, so it assumes no free native negation.
Output may be a constant or signed parameter, auxiliary or gate output.

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

Fixing either endpoint opposite its forbidden coordinate leaves its partner
free; these safe fixings give an a-flat in D, mapped injectively by Phi.
Thus a constant output is impossible for a>d. If a>=d+2, an affine
output is impossible: on a safe a-flat its constant fiber has dimension >=a-1.
A product of affine parents is also excluded: on that flat a nonconstant factor's zero
fiber has dimension a-1; constant factors give a constant output.
This excludes constants, parameter/auxiliary literals and terminal
gates with affine parents. Until a<=d+1, lose eight units of V per lost
unit of a, and never increase V at fixed a.

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
free x is Q=x+ell(other parameters). Exchange their roles: Q is
free and x auxiliary, without expanding any equation. If b=A^(-1)N_x,
then b_Q=1. Replacing A's Q-column with N_x gives determinant
det(A)b_Q=1, preserving invertibility and row widths. Update Phi by
x=Q+ell. Both x before and Q after have a child, so the exchange preserves
I, s, M, q and V.

If free x has exactly one child, a row E, N_x is nonzero. Some Q thus
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
   child B, its gate output or sole auxiliary is determined by the pair
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

After each exchange, simplify and clean pairs if possible, completing an
action before restarting. Otherwise retain G: at most two exchanges give
distinct parameter parents. Output exclusion and dead-gate cleanup give
G a gate child; it cannot have a row child.

### Restrictions and the sharing cases

Each restriction lowers a by one. The new safe dimension
is >= d+1, so the restricted output is nonconstant. Freeze old gates
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
old gates disappear, thus.

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

Only p,r become free; p is then fixed and r unused. Thus

    Delta M=-3, Delta q=-2, Delta I<=0, Delta s<=-4,
    Delta a=-1,        Delta V<=-4-4=-8.              (7)

These cases exhaust the overlaps of partner children and descendants.
Every charge belongs to one restriction and one distinct old gate.

### Termination and conclusion

Every completed action preserves representation on the current D and
strictly decreases the nonnegative integer M+s: eliminations decrease s,
pair rules decrease M, and (5)–(7) decrease s. Neither M nor s increases.
Preparatory transfers decrease the count of gates; there are finitely
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

# Exact finite assembly of the standard and affine breakers

The acyclic component map is in sampler-and-disperser.md Section 0.
Constants use exact integer/rational operations or proved terminating searches,
never extractor-property or source-distribution tests. Logarithms are binary:
lg(x) is the least nonnegative integer j with 2^j >= x for rational x>0.
Analytical error powers require no floating-point computation.

## 1. Inputs and their precise scope

The component map supplies all proofs. Historical paths:
`../circuit-log-primitives-20261105T0130/seeded-linear.md` and
`../circuit-log-extractor-research-20261105T0106/notes/affine-transformation.md`
(Sections 1--4, included in wrappers Sections 1--2).
Seeded-linear uses exact K/seed recurrences, permitting truncation/padding;
constant-advice-breaker.md8.2 gives seed-retaining 14-copy F# (8.1 unselected;
algorithm restricted1--7, extension8.2, wrappers4.2).
D3/D3' (primitives7) gives j-way splits, fixed R, U1 padding and admissible
fields to Lemma P; base (L) uses7f. Lemma P/PXT/PXG (7--7g) precede F#;
CZ/weak-AGHP, AGHP, Ext0,Ext2,F# precede NB (wrappers4), then wrappers2.1
with C=NB, A# (outer3), and B9 (Section5) supplies the outer breaker.
Advice uses cz-conditional-entropy.md8, eta=256000000/2109543000001;
finite-approximate-ring.md5g (sampler-and-disperser.md Section SP5g) gives gamma=256/20481 for primes>=2^3200000,
closing CZ's expansion premise. The CZ source and prime size guards below imply
that threshold. prime-fields-and-encoding.md gives setup.
ordinary-condenser.md contributes only (30'), proved by M/I in the final
sum-product Section SP6 and used in RG. No GUV/CGL asymptotic theorem or
recursive-sampler-schedule.md contract (F) is assumed: d_O,k_O and Om are
explicit A# slots, with B9 bounding NB.

## 2. Primitive convention

All selected breaker primitives are LemmaP, exact PX or directly proved PXT/PXG primitives with
explicit guards. No older unselected asymptotic/envelope contract is used.

## 3. The affine assemblies and the standard-breaker interface

### 3.1 The standard interface required here

A standard one-tampering breaker C is used through its actual seed length
d_C and threshold k_C at each call: arbitrary correlated source pair
independent of arbitrary correlated seed pair, honest seed uniform, distinct
advice, both seeds retained, error <=2^-e; its guards are explicit (k_C<=u
and d2<=u for NB). No interface constant c is used. The direct multicopy constant-advice breaker
(restricted-nonmalleability-and-breaking.md Section8.2) supplies NB's inner
F# at fixed advice length4; it needs only independent source/seed tuples. The standard breaker NB (affine-and-standard-wrappers.md
Section 4) works at every a>=1.

### 3.2 The selected assemblies

Only the outer A# uses affine-and-standard-wrappers.md Section2.1.
The inner F# uses the direct multicopy theorem, since NB's call has zero
additive shifts. Its independent-tuples interface suffices exactly for
wrappers4.6 step7. The outer A#, NB and every component retain their zero
branches when a supplied parameter guard fails; no selector remains.

## 4. Fix the uniform-source advice condenser first

Fix eta=256000000/2109543000001, theta=eta/8. The sharp finite-ring theorem supplies
gamma=256/20481; 103eta<gamma since
2109543*1000000<2109543000001. Also eta<1/8000.
(verify.py: CZ.expansion, CZ.eta)
NB's repaired honest Ext0 output is exactly uniform, so it uses Lemma U
of affine-and-standard-wrappers.md Section3.5. It does not use the
general-source capacity/deficiency choice theta=eta/200.
Put sig=Ee+4 and s=floor((eta*r-17sig-1)/8). NB's finite parameter
algorithm proves all Lemma U guards, s>=m0+2e_L and s>=a+40.
Its advice row error is at most14*2^-sig<2^-Ee. Advice is the actual
a-bit label, not a silently enlarged source deficiency.

## 5. The standard breaker NB and Lemma B9

Theorem NB (affine-and-standard-wrappers.md Section 4) is a standard breaker
in the sense of Section 3.1, with guards k_C<=n and d2<=n, at every a,z,E>=1.
In the dependency order above it uses direct multicopy F#, the advice condenser with
the constants of Section 4, the AGHP extractor and the primitives Ext0, Ext2.
Only the outer A# uses NB; B9 bounds that call.

**Lemma B9 (NB at the outer call).** NB is called at source length
n'=d_O<=n/2, output z=82Om, advice a=pS+2+lgC<=Om, exponent e,
where Om=e+2L+40>=2^40 (supplied by S>=2^40 in (5'')). Then d_C<13792Om and k_C<2874390323Om.
Its guards hold: k_C<=n' by the outer ledger, d2<=d_C<=sC<=d_O=n'.

Proof. Om>=2^40. The outer slots and (O3) give 2delta_O<=H, hence
n'<=delta_O L<=HL/2<=n/2<2^L and lg(4n')<=L+2.
Here eps=71/137438953472. On [2^k,2^(k+1)), k>=40,
8lg(Om)+240<=8(k+1)+240<=eps*2^k: equality at k=40,
and the numerator increases by8 while 2^k doubles.
(verify.py: B9.log-endpoint, B9.log-induction)
Ee=e+4<=Om and w_N=max(2^20,Ee+44+2lg(4n'))<=Om+12<=Om+14<=(1+eps)Om.
The height bound38 covers Ext0=PXT[22,11300207/1000000,1] and Ext2=PXT[5,19574499/1000000,6/5].
The capped F# unit w_F=max(2^20,Ee+2lg(4n')+32,ceil(z/128))<=Om:
the logarithmic term is <=e+2L+40 and z=82Om<=128Om.
Its charged slots give m1<=dF*Om and k_F<=kF*Om, where
 dF=1140949273153177/274877906944+eps, kF=605411727703975008315/274877906944+eps.
(verify.py: F#.seed-F through F#.input-d; certificates in wrappers4.2)
With Section7g's exact Exti coefficients define
 M0=c2*dF+pi2*(1+eps)+eps, K2bar=Lam2*dF+g2*(1+eps)+eps,
 Rnew=(M0+pi0*(1+eps)+7+7eps)/(theta-7c0)+eps,
 d0bar=7c0*Rnew+pi0*(1+eps)+eps, d1bar=2(M0+d0bar+3)+1.
Here theta=32000000/2109543000001>7c0.
(verify.py: B9.denominator)
NB's r numerator is
m0+pi0*w_N+(41/8)Ee+a+50; its last three terms are <=7Om.
Ceilings cost <=eps*Om, giving m0<=M0*Om, K2<=K2bar*Om,
r<=Rnew*Om, d0<=d0bar*Om. Rational substitution gives
 7Lam0*Rnew+g0*(1+eps)<2874390323,
 kF+14Rnew+84+eps<2874390323.                              (B9a)
(verify.py: B9a.K00, B9a.Kother)
They bound K00 and k_F+2v+z+2Ee; the first is larger.
For AGHP, 2e_L<=d0+3Ee and
lg(r*m0+1)<=2lgOm+lg(ceil(Rnew*M0)+1)+1<Om/4.
The fixed log constant is <240, and 2(k+1)+240<2^k/4 for k>=40
by endpoint and doubling. Hence d1<=d1bar*Om. Finally
(verify.py: B9.log-constant, B9.quarter-endpoint)
 2d0bar+2d1bar+dF+K2bar+1<13792.                          (B9b)
(verify.py: B9b.seed)
All base/hash degrees are O(log actual n); larger extensions use towers.
NB applies U7 only to exactly uniform repaired V, uses W2 for its weak
seed, and retains both complete original seeds.

## 6. Uniformity, exactness, and remaining obligations

F# and A# have closed slots; NB uses exact formulas. Primitive recurrences
compute actual seeds and thresholds; envelopes only determine padding.
The evaluator receives no entropy witnesses, good-row labels, decompositions,
coupling coins, proof histories or per-length truth tables.

For fixed t and m,a,E=O(log n), every local dimension is O(log n).
At selected PXT/PXG calls (linear-and-bilinear-primitives.md Section 7g),
F#, Ext0 and Ext2 have O(log n) base/hash degrees by their slots and B9.
Their canonical setup costs 2^O(log n); larger extensions use towers.
CZ/AGHP/bilinear setup costs 2^O(local dimension), hence polynomial in
actual n here. This is not joint polynomial time for arbitrary t,n:
the outer assembly fixes t=4C-1 and all constants once.

For ordinary seeded primitives or standard breakers with the complete
source tuple independent of the complete seed tuple, append ignored suffixes by their
kernel given all used prefixes: outputs depend only on source and prefixes.
This covers F#/A# primitive pads (wrappers Section 4.2, outer Section 3)
and NB's shift-free inner call (wrappers Section 4.6, step 7), not arbitrary
additive B-tuples correlated with suffixes.

# Tagged affine-permutation logarithmic-loss family: outer assembly

All constants below are fixed once; none is length-dependent advice.
The affine-breaker contract (F') is proved below from the outer A# (Section 3) and Lemma
B9.

Composition credits: Chattopadhyay--Liao's sumset framework,
https://arxiv.org/abs/2110.12652v1 ; Xin Li's advice-condenser architecture,
https://eccc.weizmann.ac.il/report/2023/023/revision/1/ .
Companion modules credit the finite seeded-linear, scalar-growth, graph,
breaker and native full-B2 constructions. No accepted result or model judgment
is a proof premise.

All logarithms are binary. Write lg(x)=ceil(log2 x) for positive x, computed
by integer bit length and a power-of-two comparison. All string lengths,
indices and loop counts are integers. All arithmetic parameter algorithms
use exact integers and rational comparisons, not floating point.

## 1. The exact inputs to this assembly

(L) Put actualL=floor(log2 n), S=floor(actualL/144)>=2^40. Base7f fixes
its displayed alpha,Ro,Rl,R0, 26 levels,
U1=2^(128+lgS), a=S+64, and error<2^-S. Sampler5.3 proves
D<=392S, padded to s=392S, K0<=kfacS, ceilkfac=12484765614829834828034.
(verify.py: base.seed, base.padding, base.kfac-ceiling)
Output 2^52ceil(alphaS)>=delta_O*actualL.
(verify.py: base.output)
All padding, output and
runtime guards are proved there. Failure1/100 needs K0+7<=HL.

(G') Theorem G' of sampler-and-disperser.md Section 4, with P=645,
p=P*B0=252840, C=4(q+1)=167780, A=2^(pS) and Lambda_S=A^(1/2)/2^15 (Lambda
there): a deterministic graph Gamma from left pS-bit strings to right s-bit
strings, of degree exactly C, such that whenever S>=2^40, every set of at
least Lambda_S left vertices has more than 2^(1-S)*2^s distinct neighbors.
Gamma is the four rows of one AB+C step followed by a short-start Cayley
walk and an s-bit prefix; its uniform evaluator costs 2^O(pS). With the
base sampler of (L) it gives (23') there: for every k_S-source and every
test of density at most epsilon0=2^-S, fewer than Lambda_S rows are
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

For n>=2 write L=floor(log2 n),S=floor(L/144),s_tag=log2 A=pS,w=2. Require S>=2^40 and s_tag>=128(L+1); the latter follows from the fixed constants below.
Build the canonical field F_(2^L), and split x into L-bit coefficients of
P_x, padding its final block with zeros. For a full pS-bit row alpha,
let beta,gamma be its first and second L-bit blocks. Define
H_alpha(x)=Prefix_2(gamma*P_x(beta)), a two-bit linear map in x.
Each (beta,gamma) has exactly A/2^(2L) row preimages. For Delta!=0,
P_Delta is nonzero of degree ceil(n/L)-1<n/L, so its root fraction is <= n/(L*2^L)<2/L, since 2^L<=n<2^(L+1). Off those roots, multiplication by P_Delta(beta)
permutes gamma and exactly 1/4 of gamma values have zero two-bit
prefix. Thus at most A*(2/L+1/4)<=257A/1024 rows (L>=2048) have H_alpha(Delta)=0.
The root bound follows by factoring X-beta at distinct roots. The supplied
finite-field proof makes the canonical search halt in 2^O(L) time.
Full alpha remains in the advice label despite repeated (beta,gamma).
No independence of hash values across different rows is asserted.

Use the existing source-linear sampler Samp(x,alpha,j), j in [C], seed
length d, tolerating test density epsilon0=2^-S with fewer than Lambda_S
all-bad rows ((23') in (G')), failure probability 1/100, at source
entropy k_S. The arbitrary enumerated field index alpha replaces [A].
For t=4C-1 and a=s_tag+w+lg(C), request the padded all-seeds-retained affine
breaker F of output length one, entropy K, and error at most

    epsilon = 2^(-S)/(2^33 C A^6).

The advice for (x,alpha,j) is the concatenation of the s_tag-bit alpha,
lg(C)-bit j-1, and w-bit H_alpha(x). Define

    R_alpha(x)=XOR_(j=1)^C F(x,Samp(x,alpha,j),advice(x,alpha,j)),
    D(x)=1 iff 2 sum_alpha R_alpha(x)>A.

Assume the input source entropy satisfies

    k >= max(k_S, K+4w+4Cd+2lg(1/epsilon)).                 (D1)

All concatenations, zero padding and majority ties are specified. No seed or advice inputs to D. The tags are internal computations on
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
all other selected calls' labels, irrespective of collisions at other rows.
If Q is fully doubled, the selected profiles determine all its collisions.
On profiles with no selected collision the honest advice is again distinct
from all tampered advice; two copies of the same label have distinct tags.
Tampered advice strings may coincide with one another in either case.
Section D2 thus applies in all these cases, with joint seed leakage
at most4Cd. There are at most sum_(i=1)^4 binom(2A,i)<=64A^4 sets Q,
and hence at most64C A^4*2^(8w) instances including both profiles and honest
choices. It suffices that

  k>=max(k_S,K+4w+4Cd+2lg(1/epsilon)).              (D1')

Choose theta=1/(1024A^2), epsilon=2^(-S)/(2^33*C A^6), and let B be the
union of their bad-seed sets. As 64*2^16*2*1024=2^33, Section D2 gives
(verify.py: D2.union)

  |B|/2^d<=64C A^4*2^(8w)*2epsilon/theta=2^33 C A^6 epsilon=2^(-S).    (D2')

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
(verify.py: D4.margin)

For the final inequality, 6r-4r^2 increases on [0,13/50],
and 6r-4r^2-2 equals -444/625 at13/50; also 4/g<1/1000.
This uses the exact matched-pair second moment rather than bounding each
collision term separately. All identities concern actual deterministic signs.

Let Y_sigma be the ACTUAL sum over ALL rows in that copy. Pointwise
|Y_sigma-S_sigma|<=delta=q<2^(-14)A^(1/2). Since g>=A/2,

  delta/sqrt(g)<2^(-14)A^(1/2)/(A/2)^(1/2)=2^(-27/2)<1/10000.
(verify.py: D4.delta-squared)

The L4 norms of S_sigma are at most2sqrt(g), and those of Y_sigma at
most3sqrt(g). The norm triangle inequality follows from Holder applied
to (|X|+|Z|)|X+Z|^3 and division by the cube of the L4 norm (zero norm
is trivial). For four factors, Holder itself follows by Cauchy-Schwarz
twice: E|abcd|<=(E a^4 E b^4 E c^4 E d^4)^(1/4).

Telescoping the four factors in any degree-four monomial shows its
expectation changes by at most4delta(3sqrt(g))^3. Each Q has sum of
absolute coefficients four. Thus

  |E Q_+(Y)-E Q_+(S)|, |E Q_-(Y)-E Q_-(S)|
    <=432delta g^(3/2)<(432/10000)g^2<g^2/8.            (D5')
(verify.py: D5.perturbation)

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

The base(L) is7f with its fixed constants, B0=392 and S=floor(L/144). Put P3=215, P=645,
p=252840, q=41944, C=167780, t=671119. For actual n>=2 set L=floor(log2 n),
S=floor(L/144), A=2^(pS), s=B0S, d=delta_O L. Request
 t=4C-1, m=1, a=pS+2+lgC, w=2,
 E=(6p+1)S+lg(2^33 C).                                   (O2)
Here lgC=18, lg(26t)=25 (26t=17449094), j=lg t=20.
(verify.py: outer.lgC, outer.lg26t, outer.j)
Thus e=E+lg(26t)=(6p+1)S+76 and
2^-E<=epsilon=2^-S/(2^33 C A^6). Advice is the pS-bit alpha,
lgC-bit index and two-bit polynomial tag, computed internally.
Actual L remains in tags, primitive units and entropy.

## 3. The outer affine breaker A# and (F')

For actual n passing(5''), put E,e as(O2), Om=e+2L+40<=c_O L,
 c_O=(6p+1)/144+2+eps=99439673353/9437184, eps=2^-20,
since S<=L/144 and 116<=eps*L.
(verify.py: outer.cO, guard.outer-intercept)
Put q=e, j=20, with integer slots
 dy=23Om, dx=82Om, r=15435805Om, sC=13792Om, s0=39Om,
 kC=2874390323Om, d_O=12156847042Om, k_O=59628486398890Om.

**Lemma A#.** Run wrappers Section2.1 at t=671119 with C=NB,
output dx, advice a of(O2), exponent e, seed padded to sC. Use
L0=PXG[13,836983/500000,7/20](input n,output sC,seed s0),
Lr=PXG[22,1149889/500000,171/200](input n,output r,seed dx),
E=PXG[5,2153593/1000000,5/32](input d_O,output dx,seed dy), and
Lm=PXT[1,7290999/1000000,1](input r,output dy,seed dx), all at exponent e.
Return m=1 bit. Then(F') holds with
 delta_O=ceil(12156847042*c_O)=128096781715697,
(verify.py: outer.delta)
 kappa_O=ceil(59628486398890*c_O)=628305775328681190.
(verify.py: outer.kappa)

Here is the exact ledger in units Om; c_V,pi_V,Lam_V,g_V are7g's
coefficients for the four ordered vectors above. Put
 TR=(t+1)(39+13792+20*(46+164))=12100964720,
 TL=(t+1)(13792+40*23)=9873517440.
(verify.py: outer.TR, outer.TL)
Rational cross multiplication gives every comparison:
 23>=c_E*82+pi_E; 82>=c_r*15435805+pi_r;
 82>=c_m*23+pi_m; 39>=c_0*13792+pi_0;
 15435805>=(t+1)*23+Lam_m*23+g_m+1;
 12156847042>=(t+1)*(39+13792)+2874390323+1;
 12156847042>=TR+Lam_E*82+g_E+t*82+1;
 59628486398890>=TL+Lam_r*15435805+g_r+t*15435805+1;
 59628486398890>=Lam_0*13792+g_0.                                  (A-ledger)
(verify.py: A.seed-E, A.seed-Lr, A.seed-Lm, A.seed-L0, A.r, A.d-NB, A.d-E, A.k-Lr, A.k-L0)
PXG's terminal-rate bounds hold for L0,Lr,E; its output guard follows
from 13792<=256*2^13, 15435805<=256*2^22 and 82<=256*2^5 in units Om.
(verify.py: PXG.L0.output, PXG.Lr.output, PXG.E.output)

Integer Om slots pay all ceilings; +1 pays q=e<=Om.
B9 supplies NB's seed/threshold. Also s0<=d_O,dy<=r,1<=r,r,d_O<=k_O.
The ledger gives Km<=r,KE,kC<=d_O,K0,Kr<=k_O.
Since k_O<=kappa_O L<=HL<=n and d_O<=delta_O L<=HL<=n,
inputs n,d_O,r are <=n and entropy guards hold.
Thus lg(4u)<=L+3; heights<=23 give
 e+h+6+2lg(4u)<=e+2L+35<=Om.
Om>=2^40 supplies7g's unit guard. Fixed logarithmic outputs give O(log n)
base/hash degrees even for n-bit inputs; canonical searches and towers
are polynomial in n. Error is
(13*2^20-8)2^-e<26t*2^-e<=2^-E,
proving(F') at every guarded length with raw seeds retained.
Pad d_F=d_O to d=delta_O L.

For suffix retention, fix an admissible left tuple A. For raw honest seed y
and right profile r=(B0,...,Bt,Y1,...,Yt), let Delta(y,r) be the output-tuple
TV over A from an independent uniform honest output with tampered tuple
unchanged. Choose a maximizing r_y. Uniform y with profile r_y is a legal
independent-right(F') instance; retaining y gives
 2^(-d_F)*sum_y max_r Delta(y,r)<=2^-E.
Any full right variable R jointly independent of A with uniform y fixes a
profile conditional on(y,R). Its retained-output TV averages Delta(y,r),
hence is <=that average of maxima. Thus ALL of R, including correlated
padded suffixes, is retained. Finite maxima are proof devices, not searches.

## 4. One fixed all-length guard and entropy coefficient

Put

  H=max(4Cdelta_O+kappa_O+T,ceil(kfac/144)+7),
  T=ceil((12p+2)/144+eps)=21071.                       (O3)
(verify.py: T)
Rawbranch=86596617920367272901,basebranch=86699761214096075202,H=86699761214096075202.
(verify.py: H.raw, H.base, H)
Allareexactfixedintegers. SinceS<=L/144 and110<=epsL,
4w+2lg(1/epsilon)=(12p+2)S+74+2lgC<=TL. Therefore
k_F+4w+4Cd+2lg(1/epsilon)<=[kappa_O+4Cdelta_O+T]L<=HL.
Also kS<=kfacS+7<=(ceil(kfac/144)+7)L<=HL,so k=HLsupplies(D1)
andK0<=n. Since p>129*144, S>=L/144-1 and L>=p+128 follow
from S>=2^40, pS>=129L-p>=128(L+1). Integer n<2^(L+1)
gives n+1<=2^(L+1), hence A>=(n+1)^128. Also L>=2048, so
all tag guards hold, including two full L-bit blocks in every row.
(verify.py: guard.tag-slope, guard.tag-length, guard.root-length)

The evaluator tests these literal finite guards:

  n>=2; S=floor(L/144)>=2^40; HL<=n.                   (5'')

The sampler-unit guard implies all lower bounds in7f,sampler3--5,
and all unit rounding inequalities above. The evaluator also computes
and checks each supplied algorithm's exact parameters and guards.
The proved(L),(G'),(F') show they all pass under(5''); no source
property, good-row witness, extractor oracle or existential search is
an evaluator input. S tends to infinity at EVERY integer actual length.
Take N=2^(144*2^40). For every integer n>=N, L>=144*2^40>=256
and S>=2^40. Also H<2^128 and 2^j>=2^128*j for every integer j>=256:
it holds at256 and persists since 2j>=j+1. Thus HL<=2^L<=n.
(verify.py: guard.H, guard.all-n-anchor, guard.all-n-doubling)
Thus this single N makes all guards hold at EVERY integer n>=N.

This proves the bound at all integer lengths n>=N; smaller lengths
may use the zero branch.

## 5. The fixed evaluator

Precompute the constants in Section 2 by their finite algorithms. On an
input x of actual length n:

1. If n<2, output 0. Compute L,S and the parameters above. If any finite
   guard in (5'') or its exact supplied parameter check fails, output 0.
2. Build the exact base extractor E_n from (L), seed-padded to s bits.
   Its data argument is x and its output length is d.
3. Build the (G') graph on left pS-bit strings: the four rows of the AB+C
   step, the short-start Cayley walk and the s-bit prefix (sampler-and-
   disperser.md Section 4). Its degree is exactly C.
4. Build the canonical field of degree L and the tags in Section D1.
   For each alpha in lexicographic {0,1}^{pS} and j=1,...,C compute

     y_(alpha,j)=E_n(x,Gamma(alpha,j)),
     b_(alpha,j)=F_n(x,Prefix(y_(alpha,j),d_F),
                       alpha || encode(j-1) || H_alpha(x)),

   using the ONE breaker parameter setting (O2). Let R_alpha be the XOR
   of these C bits.
5. Output 1 iff strictly more than half of the A row bits are one;
   otherwise output zero. This is exactly Section D1's evaluator.

F_n uses the raw d_F-seed; ignoring the padded d-seed suffix gives the same
evaluation. Component algorithms fix constants, fields, primes, orderings
and zero branches; no seed, basis, field advice, entropy witness, good row,
tampering or source decomposition is input.

The base sampler tests density epsilon0=2^-S. By (23'), for every test of
density at most epsilon0 fewer than Lambda_S rows are wholly bad, except on
source probability 1/100. Source
linearity survives graph index selection. The degree and advice lengths match Section D1, and (D1) was checked
above. Sections D2--D4 prove the affine-permutation property at k=HL,
and in particular the directional property.

## 6. Uniform polynomial time in actual input length

One-time constant computations terminate independently of n and belong to
the machine constant. Parameter arithmetic uses O(log n) bits beyond fixed
constants. There are A*C<=C*n^p breaker calls.
The exact base has O(log n) base/hash degrees, so enumeration costs
2^O(log n); its O(n) initial relative extension uses polynomial-time towers,
not 2^O(n) enumeration. The graph prime p_m has m+1 bits, m=pS/3;
its search costs 2^O(pS), polynomial over all A rows. Walk fields have
fixed-algorithm arithmetic.
By (F'), each call is polynomial in n for fixed t=4C-1 and a,E=O(log n).
Ext0 input/output/error and enumerated binary degrees are O(log n) by7g
(PXT/PXG); setup costs2^O(log n). Larger extensions use towers; local primes
and other field setup also cost 2^O(log n). Tags cost 2^O(L), all A
polynomial evaluations cost polynomial time, and majority uses an O(log A)
counter. The total evaluator is polynomial in actual n, with fixed,
possibly enormous degree/multiplier. It enumerates neither extractor
truth tables nor source distributions.

## 7. Native circuit bound, NP eligibility and expression

The assembly is nonconstant under every fixed-point-free affine permutation
on a source flat of dimension at least k=HL. Hence, for every
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
the output reference. Pair rule1 is thus impossible; rule3 creates an
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

For the fixed evaluator and every integer n>=N, the bridge gives
C_fullB2(f_n)>=4n-4HL. Since L=floor(log2 n)<=log2 n,

  C_fullB2(f_n)>=4n-K*log2(n)   for every integer n>=N,             (O4)
  K=4H=346799044856384300808, K0=0.
(verify.py: K)

Here the additive remainder is zero. Component entropy thresholds named
K_0 or K0 are separate quantities.
All constants and zero guard branches belong to one fixed total P
evaluator. Set f_0=0. The NP verifier accepts (x,empty) exactly when this
evaluator outputs one and rejects every nonempty witness, with witness
bound zero. Tags, seeds and field descriptions are internal, not circuit
inputs. The common threshold holds at every integer length. The declared
explicit remainder is exactly K*log2(n) with these fixed integers, and
log2 n/n tends to zero, giving (4,1).
