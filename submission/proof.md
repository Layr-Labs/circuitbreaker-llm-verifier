# Fixed logarithmic-loss family and native full-B2 theorem

This file and its five companion proofs support the exact declaration in
bound.json. Mathematical arguments are supplied below.

## 0. Component map

Material source: kaankolcu, submission0165533a-47ba-4e5b-a0c8-3cc7605b3da0,
immutable56cb98a8efa90be26a2804e41f37e9ec143e48f5. Its complete family,
NB, graph and native bridge are retained here. Used-prefix entropy and
bad-set lifting build on submission865cdc91-9181-45f7-a233-ad62151d36f7.
No independent rediscovery or prior-judgment premise is claimed. The new
analysis and explicit coefficients are stated in outer SectionsD1--D3,4,7.

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

1. Finite seeded-linear extraction and the constant-advice standard breaker.
2. Lemma A (Section 3) at degree 21 around it: NB's inner breaker F.
3. The CZ/IP advice condenser, the AGHP extractor and F: the standard
   breaker NB (affine-and-standard-wrappers.md Section 4).
4. Lemma A at the sampler's degree with C=NB (Section 5): the outer
   generalized affine breaker.

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
`../circuit-log-extractor-research-20261105T0106/notes/affine-transformation.md`.
The standard breaker is Theorem NB (affine-and-standard-wrappers.md Section 4);
it replaces the earlier equal-length and general-input wrappers.

For the advice condenser we use `cz-conditional-entropy.md` Section 8 with
eta=1/10000000, and its expansion input from `finite-approximate-ring.md`
Section 5: gamma=1/10000 and prime threshold 2^20000. The CZ source and prime
size guards below imply that threshold. Its canonical prime setup is the one
in `prime-fields-and-encoding.md`. The sampler uses the suppliedRGsupportgrowth proof and short-startwalk. These
are declared mathematical inputs, not proved again or strengthened here.
The ring theorem has exactly the expansion contract left conditional in the
CZ file. Substitution with the displayed gamma and eta closes that premise.

This note does NOT assume normalized contract (F) in
`recursive-sampler-schedule.md`, nor any closed-form degree envelope: Lemma A
gives the outer breaker's seed and threshold as Qw and K_c*w, Lemma B bounds
NB's two coefficients at the outer call, and outer Lemma F1 bounds w.

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
P_u<=(h+25)(S_u+b+1) and G_u>=2B(S_u+b+1). If B>=256(h+25), the same inequalities
give the stronger512 P_u<=G_u, which the asymmetric assembly below uses.

A useful envelope, valid with the definitions above, is

    G_u <= C(e+h+lg(u+1)+6),    C=64B(h+1).                (2.2)

Indeed rho_u=e+h+10+2lg u <=2(e+h+lg(u+1)+6),
S_u<2(4h+1)rho_u, and
4b+S_u+5 <= (16h+9)(e+h+lg(u+1)+6).
Both G and P are nondecreasing in u. No unspecified extractor constant
appears in these envelopes.

## 3. A polynomial-degree affine assembly from any linear standard interface

### 3.1 The standard interface required here

Suppose a standard one-tampering breaker C is supplied with an integer c>=1.
On source length u, advice length a, output z and requested exponent e, its
actual seed length d_C and source threshold k_C satisfy

    d_C,k_C <= c[z+a+e+lg(u+1)+1].                          (3.1)

It accepts arbitrary correlated source pairs independent of arbitrary
correlated seed pairs, honest seed uniform, and retains both seeds with
error <=2^-e. Its only remaining domain guard is k_C<=u; any other finite
guard must already be built into the interface. The constant-advice breaker
is used only at its one fixed advice length 4, inside NB's inner breaker.
The standard breaker NB (affine-and-standard-wrappers.md Section 4) works at
every a>=1.

### 3.2 Lemma A: a two-parameter asymmetric affine assembly

Let n,t,a,m,E>=1 and integers delta0,kappa>=1 be given. Choose a fixed
integerR>=ceil((t+1)/128)+3 (default equality), and a supplied integer
h_A>=14 satisfying2^h_A>=256(h_A+25),2^h_A>=12max(delta0,R)
(default the least such height). Put

 e=E+lg(26t),j=lg t,B=2^h_A,
 T_c=(t+1)(2+2delta0)+ceil(65j(t+1)/32),
 K_c=max(2delta0+1,T_c+(t+2)R+2),
 Q=max((t+1)(1+delta0)+kappa+1,T_c+t+4).                     (A1)

Select a dyadic widthw passing

 w>=ceil(m/R)+a+e+1, w>=lg(Qw+1), w>=G_U,U=max(n,Qw),         (A2)

whereG_U is Section2's height-h_A budget at exponente. Let C satisfy the
standard retained-seed contract of Section3.1, and suppose its call at
source d=Qw,outputw,advicea,exponente has actualseed<=delta0*w and
threshold<=kappa*w. Set

 d=Qw,k=K_c*w,s0=dx=w,dy=w/128,d0=delta0*w,r=Rw.              (A3)

All dimensions are integral: w>=G_U>=8B(e+h_A+4)>128, so passing dyadicw
is divisible by128. Pad C's seed to d0 and use the common h_A,e primitive:
L0 inputn/outputd0/seed s0; Lr inputn/outputr/seed dx;
E0 inputd/outputdx/seed dy; Lm inputr/outputdy/seed dx.
Return Prefix_m of the final row from wrappers Section2.1. Ifk<=n, this is
a degree-t generalized affine breaker of outputm,seedd,thresholdk,error
<=2^-E, retaining all seeds under that section's arbitrary-correlated-copy
left/right contract. Otherwise return0^m with no extraction claim.

Proof. The supplied height gives512P_u<=G_u. All ambient lengths are<=U:
L0,Lr readn; E0 readsd; Lm readsr=Rw<=K_c*w=k<=n. HenceP_u<=w/512 and
G_u<=w for every actual call. SinceB>=12max(delta0,R),

 D_L0,D_Lr<=w/2+w/512<w,
 D_E0<=6w/B+w/512<=w/256<=dy,
 D_Lm<=6dy/B+w/512<=w/256<w,
 K0<=(2delta0+1)w,Kr<=(2R+1)w,KE<=3w,Km<=(1+1/64)w.        (A4)

HereB>=2^14>=3072 justifies6/B<=1/512. Pad each actual seed to its stated
length. Legality:K0,Kr<=k<=n;KE<=3w<=d;Km<2w<=r sinceR>=4;
C's source threshold<=kappa*w<=d and its seed fitsd0.
Withq=e, check every wrapper contract(3). The history budget is

 T=(t+1)(2+2delta0)w+j(t+1)(4w/128+2w)<=T_c*w.

The first test givesm<=r ande<=w. Also s0<=d,dy<=r, and
r>=((t+1)/128+3)w>=(t+1)dy+Km+e. Finally
k>=max(K0,T+Kr+t*r+e), since the latter sum is at most
(T_c+(t+2)R+2)w. The two source-seed requirements are

d>=(t+1)(s0+d0)+k_C+e andd>=T+KE+t*dx+e,

by the two maxima inQ. Every primitive and C has error<=2^-e. The general
retained-seed proof therefore gives(13*2^j-8)2^-e<26t*2^-e<=2^-E.
It applies to these unequal dimensions without any change in its invariant.
The construction is a literal specialization, not an unproved asymmetric
transformation. Callers may add tests, as NB's inner F does. The selector
halts by Section3.3, and the outer bound is made explicit in LemmaF1 below.

### 3.2a Inner Lemma A_CA64: smaller constant-advice seed padding

The wrapper'sC map isEXPLICITLY C(y,r1,alpha)=Prefix_(v/4)(C_base(y,r1,alpha)),
whereC_base is the statedconstant-advicealgorithm requestedatoutputv.
ThusC has outputdx=v/4 exactly,matchingLr'sseedtype. Projectingboth
honestandtamperedC_baseoutputs preserves itsstrongordinarybreaker
guarantee andallretainedseeds byTVcontraction. The unusedbaseoutput
suffixisdiscarded,notrevealed orchargedinthehistory;onlytheactualdx-bit
C outputsenterinitialization. The conservativeC_base seedandthreshold
boundsbelowremainvalidforthisprefixmap.

UsedONLYbyNB'sinnerF. Fix t21,j5,delta1,kappa72,R17,h_A14,B_i2^14.
Atpassingdyadicwidthv use(A2),s0=dx=v/4,dy=v/128,d0=v,r=Rv,
Tc114,Kc507,Q122,d=Qv,k=Kcv. Allfractionaldimensionsareintegralsincev>128.
Use the samefourprimitives andthe SAMEconstant-advicealgorithm, but its
innerFextratests arenowv>=64(e+j0),v>=64Gca,j0=lg130=8.
HereGca isrestrictedSection8'sbudgetatexponente+j0,sourceu=Qv.

Theconstant-advicewidthselector(15) passesatwca=v/64: itcoverserrore+j0,
128wca=2v>=outputv,andU=max(Qv,2v,56v/64)=Qv. Hence its ACTUALleastwidth
is<=v/64 (noequalityisasserted). Its actualseed56wca<=7v/8<v andthreshold
36*128wca<=72v<=Qv. Its finalerror130*2^-(e+j0)<=2^-e.
All16advicevalues andthecomponentproofareunchanged.

B_i>=256(14+25),B_i>=48R,so512P_u<=G_u<=v forallambients<=max(n,Qv).
L0,Lr seedcost<=6max(1,R)v/B_i+v/512<=v/8+v/512<v/4.
E0 outputs v/4,soitsseed<=3v/(2B_i)+v/512<v/128;Lmseed<v/4.
ThresholdsK0<=3v,Kr<=(2R+1)v,KE<=3v/2,Km<=(1+1/64)v.
Underk<=n,LmlongrowRv<=k<=n;otherambientsarenorQv. ThusK0,Kr<=k,
KE<=d,Km<=r andC'snewsource/seedboundsabovefit.
ChargeALLinitialretainedobjects:
 T0=22(s0+2d0+dx)=55v,
 T=T0+5*22(4dy+2dx)=113.4375v<=114v.
Firsttestgivesm<=r,e<=v; r>=22dy+Km+e sinceR17>=3.
Theuniformseedsourceconditionsare
 22(s0+d0)+72v+e<=100.5v<122v=d,
 T+KE+21dx+e<=(114+3/2+21/4+1)v=121.75v<122v=d.
TheleftthresholdconditionisT+Kr+21r+e<=(114+23R+2)v=507v=k.
Allwrappercontracts(3)hold. The unchangedretained-sideinduction yields
error(13*2^5-8)2^-e<26*21*2^-e<=2^-E. Callswithk>n returnzeroasbefore.
Thecommonselectorstillterminates;the64basedextratests areweakerthan
itsformer128basedtests,whose finiteupperenvelopesremainvalid.

### 3.3 The selector halts, with a polynomial bound

With h=h_A, B=2^h_A and C=64B(h_A+1) as in (2.2), put

    V=e+h+m+a+lg(n+1)+lg(Q+1)+10.

Let w* be the least power of two >=16 C^2 V. Then
w*<32 C^2 V. As n,Q,w*>=1,

 lg(max(n,Qw*)+1)<=lg(n+1)+lg(Q+1)+lg(w*+1),
 lg(w*+1)<=6+2lg C+lg V<=6+2C+V.

It follows that the bracket in (2.2) is at most 2V+2C+6<=4CV,
so G_max(n,Qw*)<=4C^2V<=w*/4. The other two tests in (A2) also hold:
the first is bounded by V, and the second by 2V+2C+6<=4CV<=w*.
Thus w<=w*<32C^2V. The selector consists only of finite dyadic trials.

For parameter choices with k>n, define the affine map as the all-zero m-bit
string and make no extraction claim. This is a total computable family.

## 4. Fix the advice-condenser constants first

Set

    eta=1/10000000,     theta=eta/200,     beta=theta/2.

For the advice condenser use 10 blocks of r bits and advice capacity
b=10s, where s=floor(theta*r). Pad an actual advice of length a<=b with zeros
to b bits before every CZ call. Thus the CZ input length is u=r+b and the
honest coordinate entropy is r-2b. Its deficit is
rho=3b/(r+b)<=30theta=3eta/20<eta/2.
Also eta*r-(3-eta)b>=17eta*r/20>2, since b<=eta*r/20 andeta*r>=64;
thus r-2b>=(1-eta)(r+b+1)+1, including the rounded flat-source guard. If eta*r>=64, the prime guard, the rounded
flat-source entropy guard and the additional 2^20000 prime threshold all
hold. The nonmalleable output length is h_N=floor(eta*(r+b)), and error is
at most 91*2^-eta(r+b), by the stated CZ input.

Set delta=2^-s and B_A=19b=190s in the advice-condenser lemma. If
r>=max(ceil(64/eta),ceil(28/theta)), then s>=4, 2^s>=10, h_N<=r, and

    h_N >=192s = B_A+s+log(1/delta).

Indeed h_N>=eta*r-1, 192s<=192theta*r, and
eta*r-192theta*r=eta*r/25>=64/25>1. Also eta*r-190s>=10theta*r,
so both error terms below have at least the stated exponent. Its precise error
therefore obeys

    e_A<=14*2^-s+91*2^(-eta*r+190s)
        <=105*2^(-theta*r/2)<=2^(-theta*r/4).               (4.1)

The last step uses theta*r>=28 and 105<2^7. It supplies good-row entropy
s>=beta*r. These constants are fixed before any breaker parameter is chosen.

## 5. The standard breaker NB and Lemma B

Theorem NB (affine-and-standard-wrappers.md Section 4) is a standard breaker
with the contract of Section 3.1 at c=c_std=3470181145477985256748824, with
the single guard k_C<=n, at every a,z,E>=1. In the dependency order above it
uses the constant-advice breaker inside Lemma A at degree 21 (its inner F),
the advice condenser with the constants of Section 4, and the AGHP extractor.
The outer assembly uses NB only through Lemma A, and Lemma B gives the two
coefficients that Lemma A needs at the outer call.

**Lemma B (NB at the outer call).** OuterLemmaA usesheight21,delta512,
kappa1710000000000andtheextraGpluswidthtestinouterSection2.
Forpassingw withQw<=n,NBseed<512w,threshold<1710000000000w.

Letn'=Qw,Ee=e+4,L'=lg(n'+1),B=2^21,bplus=e+37,rho=bplus+2lg(4n)+2,
Gplus=2B(4bplus+S+5),S=pow2ceil(85rho). Actualheight21primitiveb=e+25
istwelve smaller,sotheextraGplustest impliesoriginalA2tests.
w>=Gplus>=8Bbplus,Ee,a,L'<=bplus<=w/(8B),bplus>=64.
InnerFusesA_CA64,R_i17,h_i14,Q_i122,K_i507. Attrialv=w/16,
K_iv=507w/16<32w<=Qw=n',soalllongrowambientsarelegalbeforeNB'sguard.
Firsttestceil(w/17)+Ee+15<=w/17+e+20<w/16,as1/(8B)<1/16-1/17.
Logtest<=L'<w/16. Innerb_i=Ee+28=e+32<=bplus,rho_i<=rho,
soG_i<=2^(14-21)Gplus<=w/128<w/16.
Constant-advice64Gca<=2^(11+6-21)Gplus<=w/16;
64(e+22)<64bplus<=8w/B<w/16. Everypassingwidth>w/32byfirstoutput/rowtest,
soexactv=w/16,d_F=122w/16,k_F=507w/16<32w.

NBheight40hasrho0<=rho+11<2rho,S0<644rho,
P0<30316rho<179w/B sincew>=170B*rho. OtherNBwidthtestsare<4w,
while4w<d_F<8w,soactualNBwidthW=8w,d_C<=63W=504w<512w.
Step5givesr<=k'(4W+2P0+4Ee+a+40)+1,k'=8000000000/3.
w>=8Bbplus,bplus>=64yield
 2P0+4Ee+a+40<(358+1/2+1/8+40/512)w/B<360w/B,
 2v_NB=20r<=(640k'+7200k'/B+20/w)w<1706675821942w.
Wrappers(K)overhead<=2261123072w,as(xE+xL+x0)/(8*2^21)=2261123072;
k_F+w+2Ee<34w. Thusk_C<(1706675821942+2261123072+34)w
=1708936945048w<1710000000000w. Q>=1710000000001 provesNB'sownguard.
AllgenericLF/NBtablesremainconservative fornewinnerheight14<=20,
Q122<=337<=773,Ki507<=1507665,andweaker64basedconstant-advicetests.

## 6. Uniformity, exactness, and remaining obligations

All parameter searches are now proved terminating: the dyadic selectors of
Lemma A (Section 3.3), of the constant-advice breaker and of NB's inner
breaker (Lemma F there), and finite constant formation. Exact seed lengths and entropy thresholds
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
primitive pads in Section 3 and to the zero-additive-shift independent-side
use of NB's inner breaker (affine-and-standard-wrappers.md Section 4.6, step 7). It does not assert a generic
kernel principle when arbitrary additive B-tuples correlate with suffixes.

# Tagged affine-permutation logarithmic-loss family: outer assembly

All constants below are fixed once; none is length-dependent advice.
The affine-breaker contract (F') is proved below from Lemmas A, B and F1.

This retained peer composition supplies the complete finite seeded-linear,
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

(L) Use D3_exact inprimitiveSection7a atj4,R11,h'=46,h92,Bbase2^92,
leaf17L,errorL,exactbinarydimensions,U1=2^(132+lgL),L>=2^26.
SamplerSection5.3 provesD<19809411/32768L<606L,K0<kfacL,
kfac=3691498529356500954201689234792331874121926549859336401/10376293541461622784,ceil(kfac)=355762731133810044533343678033065949.
OutputT*L=17*2^92*L=84179922671405858693140447232L is not truncated. ExactK0<=n is tested.
Failure1/100 needskS<kfacL+7,not6TL+7. This isstrong,source-linear,uniformP:all enumeratedbase/hashfielddegrees
areO(L),so setup2^O(L) ispolynomial inactualn;largerextensionfields
arebuiltbytowers,notenumeration. All breakerprimitivesretainoriginalH.

(G') Use the suppliedRGsupportgraph in samplerSections3--4 atB0606,P32990,
p=8970B0=5435820,q=2989B0/2=905667,C=4(q+1)=3622672.
AtL>=2^31 it expands everyLambda_S=A^(1/2)/2^15 leftvertices beyond
2^(1-L)*2^s rightvertices,s=B0L,A=2^(pL). Its fixeduniformalgorithm and
samplerSection5.2 givefewerthanLambda_S all-bad rows foranytestdensity2^-L,
exceptonsourceprobability1/100 atentropy kS.

(F') For every n passing the guards (5') of outer Section 4, Lemma A with
C=NB, at the request (O2) of outer Section 2, returns an exact seed length
d_F, source entropy requirement k_F and a one-bit generalized affine
breaker, where (outer Lemma F1)

  d_F=Qw_A<=TL,   k_F=K_c*w_A<=TL.                                      (F')

For arbitrary correlated tuples (A0,...,At) independent of
(B0,...,Bt,Y0,...,Yt), with H_inf(A0)>=k_F, Y0 uniform on d_F bits, and
alpha0 different from every alphaj, its honest output on Xj=Aj+Bj is
2^-E-close to uniform jointly with every tampered output and EVERY seed
Y0,...,Yt. The field is F_2, so addition is XOR. Neither source copies nor
seed copies are required to be deterministic functions of their honest copy.
For the fixed t and a,E=O(log n), its complete evaluator is polynomial in
the actual n (Lemma A's evaluator with NB, whose cost is stated in
affine-and-standard-wrappers.md Section 4.7).

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

    R_alpha(x)=XOR_(j=1)^C F_use(x,Prefix_dF(Samp(x,alpha,j)),advice(x,alpha,j)),
    D(x)=1 iff 2 sum_alpha R_alpha(x)>A.

Assume the input source entropy satisfies

    k >= max(k_S, K+4w+4C*d_F+2lg(1/epsilon)).                 (D1)

Here F_use is the raw d_F-seed map in(F'); the padded F evaluates precisely
this prefix expression. Set K=k_F=K_c*w_A, the actual raw threshold in(F').
The evaluator is fixed by the present constants and selectors, including all guards.
All concatenations, zero padding and majority ties are specified. There
are no seed or advice inputs to D. The tags are internal computations on
the actual n-bit x.

## D2. Bad seeds with arbitrary correlated left copies

All seeds in this lemma have length d_F and F_use denotes the raw(F')
map. Each seed map is the d_F-bit prefix of its corresponding sampler map.
Fix an arbitrary correlated tuple (A_0,...,A_r), fixed advice differing
from the honest advice, and deterministic seed maps L_i of A_i. Assume
Hinf(A_0)>=K+ell+lg(1/epsilon), where the joint seed leakage has at most
2^ell values. Call w bad if some b_0,...,b_r,w_1,...,w_r make the honest
bit in the tuple

  F_use(A_0+b_0,w+L_0(A_0),advice_0),
  (F_use(A_i+b_i,w_i+L_i(A_i),advice_i))_(i=1)^r

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
Section D2 therefore applies in all these cases, with joint USED-prefix seed leakage
at most4C*d_F. There are at most sum_(i=1)^4 binom(2A,i)<=64A^4 sets Q,
and hence at most64C A^4*2^(8w) instances including both profiles and honest
choices. It suffices that

  k>=max(k_S,K+4w+4C*d_F+2lg(1/epsilon)).              (D1')

Choose theta=1/(1024A^2), epsilon=2^(-L)/(2^33*C A^6), and let B_used be the
union of their bad-seed sets in the d_F-bit seed space. As 64*2^16*2*1024=2^33, Section D2 gives

  |B_used|/2^d_F<=64C A^4*2^(8w)*2epsilon/theta=2^33 C A^6 epsilon=2^(-L).    (D2')

Lift B_used to B=B_used x {0,1}^(d-d_F). Its full-d density is EXACTLY
the used-prefix density, and a full sampler seed avoids B precisely when
its used prefix avoids B_used. This lift is a proof object. No ignored
source-dependent suffix is conditioned on or assumed independent.
B depends on the left source law and tau, not the right source point;
the program never constructs B. Apply the sampler to both B_0 and B_1,
each of entropy k. Except on probability2/100 both conclusions hold.
Fix such an x_2. Intersect the good row sets to obtain G, omitting
q<=2(Lambda_S-1)<2^(-14)A^(1/2) rows. Since 2^(-14)A^(-1/2)<=1/256,
g=|G|>=255A/256>=A/2. Each sigma,alpha in {0,1}xG has its own
neighbor j_(sigma,alpha) whose honest seed lies outside B. The two copies
need not choose the same j. Source-linearity gives every seed identity
with shifts B_sigma=x_2 or Lambda x_2+b and offsets
Prefix_dF(Samp(B_sigma,alpha,j)). Prefixes preserve source-linearity.
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

Fix baseheight92 (46four-waylevels),leaf17L,B0=606,P=8970,
 T=17*2^92=84179922671405858693140447232,p=5435820,C=4(q+1)=3622672,t=4C-1=14490687,j=24,
 R=ceil((t+1)/128)+3=113212,suppliedh_A21,B=2^21,
 T_c=1026(t+1)+ceil(65*j*(t+1)/32)=15573866928,
 K_c=T_c+(t+2)R+2=1656093749998,
 Q=max(513(t+1)+1710000000000+1,T_c+t+4)=1717433722945.
Allinput-independent;B>=256(21+25),B>=12R,max(Q,K_c)=Q.
Graphdegreeis4(q+1),not2p/3.

For each actual n>=1 let L=lg(n+1), A=2^(pL), s=B0L and d=TL. Request

  t=4C-1, m=1,
  a=pL+w+lg(C), w=2,
  E=(6p+1)L+lg(2^33*C).                                       (O2)

Here lgC=22,lg(26t)=29,e=(6p+1)L+e0,e0=84.
AdviceconcatenatespL-bitalpha,lgC-bitj-1,two-bittag,exactlya bits.
The sameerrorrequest2^-E<=2^-L/(2^33*C*A^6) holds.
For this OUTER call only,addthewidthtestw_A>=Gplus_U,U=max(n,Qw_A),where
B=2^21,bplus=e+37,rho_U=bplus+2lg(4U)+2,Splus_U=pow2ceil(85rho_U),
Gplus_U=2B(4bplus+Splus_U+5). Actualouterprimitiveb=e+h_A+4=e+25 istwelve
smaller,soG_U<=Gplus_U. LemmaA permits extra tests; primitives themselves
retaintherequestederror e. Thisispartofthefixedouterselector,effecting
nochangeinwrappercontract,advice orrandomness assumptions.

Below,w_A denotes the leastdyadicouterwidth passing(A2) ANDtheGplus test;
w=2 in(O2) is onlythetaglength.

## 3. Termination and the exact feedback inequality

Put

 g1=1044p+514=5674996594,g0=4(e0+37)+170(e0+43)+5=22079.
Exactfeedback4B*max(Q,K_c)*(g1+1)<=T:                         (F1)
 81758985114152743184053043200<=84179922671405858693140447232.

**Lemma F1.** If L>=g0 and HL<=n (H from (O3) in Section 4), the least power of
two w_A passing (A2) andtheextraGplus test at (O2) satisfies w_A<4B(g1+1)L. Hence Qw_A<TL<=n,
d_F=Qw_A<=TL and k_F=K_c*w_A<=TL.

Proof. For n>=1, lg(4n)=lg n+2<=L+2. For the outerbudget useb=bplus=e+37 ratherthanactuale+25, rho_n=b+2lg(4n)+2<=e+43+2L and
S_n=pow2ceil((4h_A+1)rho_n)<2(4h_A+1)rho_n. Substituting e=(6p+1)L+e0,

  Gplus_n=2B(4b+S_n+5)<2B[4(e+37)+170(e+43+2L)+5]=2B(g1*L+g0).

Let w*=pow2ceil(2B(g1*L+g0)); then Gplus_n<w*<4B(g1*L+g0)<=4B(g1+1)L, using
L>=g0. By (F1), Qw*<TL<=HL<=n, so U=max(n,Qw*)=n and the third test of (A2)
andtheextraGplus test pass atw*. The first passes since m+a+e+1=(7p+1)L+4+lg C+e0<2B*g1*L. The
second: lg(Qw*+1)<=lg(n+1)=L<w*. So w_A<=w*, and
Qw_A, K_c*w_A<4B*max(Q,K_c)(g1+1)L<=TL<=HL<=n.

So Lemma B's hypothesis holds at every n passing (5'), Lemma A applies there
with C=NB, and (F') holds. Pad the breaker to d-bit seeds by taking its first
d_F bits. All other finite guards are still checked explicitly by the
evaluator.

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

## 4. The original guard and a sharper used-prefix entropy requirement

Keep the evaluator's ORIGINAL guard constant

 H=144709410103943034489126839574476258663392646546947510738880357178915.                                  (O3)

The output-amplifying selector in (A2), and delta0=512,kappa=1710000000000,
fix the present evaluator and its selected dimensions. They differ from
the predecessor. The sampler usesD3at46four-waylevels,R11,leaf17L,B0606,P8970,RGsupportgraph.
The literal original zero guard H is retained and exceeds(4C+1)T+12p+3.
The proof of Lemma F1 gives the sharper literal bound
w_A<4B(g1 L+g0). Since w_A is an integer, put

 D_*=K_c+4C Q=24886797895968186158,
 W_1=4B g1=47605321828401152, W_0=4B g0-1=185212076031,
 H_1=D_*W_1+12p+2=1184744023115942154450568286702883858,
 H_0=D_*W_0+74+2lg C=4609335504076190523052657779016.

Then w_A<=W_1 L+W_0. The tagged-parity conditioning costs 4w=8 bits,
and 2lg(1/epsilon)=(12p+2)L+66+2lg C. Thus the actual requirement(D1)
with USED-prefix leakage is at most

 K+4C*d_F+4w+2lg(1/epsilon)
 =D_*w_A+(12p+2)L+74+2lg C <=H_1 L+H_0.

Defineproof-onlyk(n)=H_1 L+H_0. Exactintegerssatisfy
 H_1=1184744023115942154450568286702883858>=355762731133810044533343678033065956=ceil(kfac)+7,H_0>=7.
HencekS<kfacL+7<=k(n). AlsoH_1+H_0<H,so k(n)<HL<=n atliteralguard.
Allrawsourceandseedlegalitychecks hold.
Ignored sampler suffixes never entered the left entropy conditioning.
The conditional-copy bad-seed proof and lifted set supplied(D1) directly;
no claim of suffix independence or prior-verdict premise is used.

The evaluator tests these literal finite guards:

  n>=1;L>=2^31;HL<=n.                                      (5')

This covers baseL>=2^26,exactouterwidth lowerbounds andwalkguard in
samplerSection3. All exactalgorithm/sourceguardsarecheckedandpass.
No test of anunknownextractorproperty orsource distributionoccurs.

There is one fixed threshold N beyond which (5') holds for EVERY integer n.
Indeed L tends to infinity with n, the walk guard holds at every
L>=2^31 (sampler-and-disperser.md Section3), and the fixed lower bound on
L then holds too. Finally H*lg(n+1)<=n eventually: choose j=floor(log2 n),
so 2^j<=n<2^(j+1) and L<=j+2. The elementary exponential sequence 2^j/(j+2)
is increasing (consecutive ratio 2(j+2)/(j+3)>1) and unbounded, so it tends to
infinity: H(j+2)<=2^j for every j from some fixed j0 on. This also proves
eventual k_F<=n.

This is an all-integer-length argument, not a favorable subsequence or a
liminf used to invent a remainder. A small-input zero branch is allowed
because only an eventual lower bound is claimed.

## 5. The fixed evaluator

Precompute the constants in Section 2 by their finite algorithms. On an
input x of actual length n:

1. If n=0, output 0. Compute L and the parameters above. If any finite
   guard in (5') or its exact supplied parameter check fails, output 0.
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
above. Sections D2--D4 prove the affine-permutation property at k=k(n)=H_1 L+H_0,
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
on a source flat of dimension at least k=k(n)=H_1 L+H_0. In particular, for every
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

For this SAME fixed evaluator, at every integer n beyond the original
common threshold, C_fullB2(f_n)>=4n-4k(n)=4n-4H_1*lg(n+1)-4H_0.
Since lg(n+1)<=log2(n)+1 for every integer n>=1, set

 K=4H_1=4738976092463768617802273146811535432,
 K0=4(H_1+H_0)=4738994529805784922564365357442651496.

Then C_fullB2(f_n)>=4n-K*log2(n)-K0 for EVERY integer actual n>=N. (O4)
Both coefficients are exact37-digit positive integers, each smaller than
the preceding exact-leaf construction's59-digit coefficients.
This is a strict improvement of the declared explicit remainder, for the
present fixed family, not merely an unencoded Big-O constant. No optimality
of these conservative coefficients is asserted.
The NP verifier accepts exactly the empty witness when this fixed total
polynomial-time evaluator outputs1 and rejects every nonempty witness.
Tags, fields and seeds are internal, not circuit inputs. All stated
all-integer eventual guards and uniformity arguments apply to this fixed family.
The explicit remainder is exactly K*log2(n)+K0, with projection(4,1).
