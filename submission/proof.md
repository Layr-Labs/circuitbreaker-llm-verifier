# Provenance
Uses the proofs and architecture of 64a9e40a-e9cc-4568-8982-fc9d9bfcb56e. This variant changes the walk width/count, Q and derived slots. No verdict is a premise.
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
This reconstruction uses the supplied circuit-reduction.md.
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
1. PXT/PXG7g and PXTM7h and exact fields feed restricted Sections1--3 and8.2.
2. Restricted8.3's five-round antichain12-copy corollary feeds wrappers4.2's F#.
3. CZ/E6, AGHP, Ext0,Ext2,F# feed wrappers4's NB.
4. NB with B9's complete seed/threshold feeds outer3's A# at the sampler degree.
Constants use exact integer/rational operations or proved terminating searches.
For positive rational x, lg(x) is the least nonnegative integer j with2^j>=x.
## 1. Inputs and their precise scope
All inputs are included. Linear Sections2--4,7a,7g--7h prove the arbitrary-source
PXT/PXG/PXTM vectors and exact fields; the same file supplies bilinear/weak-seed/AGHP.
Only BASE uses affine Section7f.
Sampler Part SP proves finite rings and complete incidence: gamma=18230601359788607954569119/193428131138340667952988160,
prime guard2^(2^90). Restricted supplies prime/encoding, CZ smoothing,
conditional tools and the five-round antichain8.2 theorem and its seven-label corollary8.3 retaining all13
seeds. Wrappers4.2 instantiates the zero-shift12-copy F# in NB4.6 step7.
Wrappers2,4 prove the all-seeds affine transformation and NB. Outer3 builds
A# at t319, B9 pays its full NB call. Sampler3--6 prove W_A/G_A/AS2 and the
general conditional bad-seed lemmas. Sampler0 resolves legacy
component names to these included bodies. Every component seed is retained.
## 2. Primitive convention
Selected PXG calls use pi_* from the legal-call bound in7g; their heights
are at most13. Establish threshold guards by PT1 before applying PXG*.
PXT generic coefficients remain unchanged.
Selected slots use the proved PXT/PXG7g and PXTM7h primitives;
explicitly named AGHP calls use the separately proved matrix extractor. Each has
its explicit source, seed and dimension guards.
## 3. The affine assemblies and the standard-breaker interface
### 3.1 The standard interface required here
A standard one-tampering breaker C is used in the outer assembly through
its actual seed lengthd_C and thresholdk_C: arbitrary correlated source
pair independent of an arbitrary correlated seed pair, honestseed uniform,
distinct advice,bothcomplete seeds retained,error<=2^-e. HereC=NB has
explicitguardsk_C<=u,d2<=u,(AG1) s>=a+40,(AG2) eta*r>=50(a+sig+2) of wrappers4.3
and works at every advicelengtha>=1 passing them; failed guards give 0^z.
NB's innerF# instead has the five-round antichain12-copy tuple contract of wrappers4.2
and restricted8.3 at seven fixed labels; NB4.6 step7 uses exactly that
zero-shift contract. The outer assembly retains arbitrary additive shifts.
### 3.2 The selected outer affine assembly
The outerA# is wrappers2.1 withC=NB and the explicit outerSection3 slots.
F# is the map of wrappers4.2. Failed public guards return the all-zero
string; every guard holds at every sufficiently large actual length.
## 4. Fix the uniform-source advice condenser first
Fixgamma=18230601359788607954569119/193428131138340667952988160,eta=gamma/(1081/25+7*2^-80)=91153006798943039772845595/41819161952109252411436045792,theta=eta/6=30384335599647679924281865/83638323904218504822872091584.
Every actual NB call haseta*r>=96*2^80,as proved
in wrappers(I3),so this guard and the independent ring-prime guard hold.
NB's repaired honest Ext0 output is exactly uniform, so it uses the six-block Lemma U
of affine-and-standard-wrappers.md Section3.5.
Put sig=Ee+4 and s=floor((eta*r-15sig-8)/6). NB's finite parameter
algorithm proves s>=m0+2e_L; with its explicit guards (AG1) s>=a+40 and
(AG2) eta*r>=50(a+sig+2) it proves all Lemma U and E6 guards (wrappers(I1),(I2)).
Its advice row error is at most14*2^-sig<2^-Ee.
## 5. The standard breaker NB and Lemma B9
Theorem NB (affine-and-standard-wrappers.md Section 4) is a standard breaker
in the sense of Section 3.1, with guards k_C<=n, d2<=n, (AG1) and (AG2), at every a,z,E>=1.
It uses F#, the Section4 advice condenser, AGHP, Ext0 and Ext2; the outer
assembly uses NB only through A#, and Lemma B9 gives its seed and threshold.
**Lemma B9 (NB at the ACTUAL outer call).** Source n'=d_O, output
z=80Om, advice a=s_tag+2+lg C, exponent e, Om=e+2L+42. Put
 eps_N=10^-18, alpha=433/(433+2Q_S)+eps_N.
The exact slots below have d_O=floor((448929385877401036698799/20000000000000000)Om)<2^30Om, and O2a gives
Om<=c_O L<1024S: c_O*Q_S*(1+2^-80)<1024 and S>=2^80 (5'').
Hence 4n'<2^32Om<2^42S and lg(4n')<=42+lg S. With Ee=e+4, the actual NB unit satisfies
 w_N=max(2^80,Ee+44+2lg(4n'))<=alpha*Om.
Indeed Om=E+2L+55 and E<=433S+29, L>=Q_S*S imply
 alpha*Om-E >= eps_N*(433+2Q_S)S+55alpha-29(1-alpha).
The constant is positive. The nonfloor part of w_N is at most
E+146+2lg S. The old dyadic certificate
5770*80+29028=490628<eps_N*(433+2Q_S)*2^80
still dominates 146+2lg S at every S>=2^80, hence under the new guard.
Also Om>=(406+2Q_S)S, and exact arithmetic gives
alpha*(406+2Q_S)>1. Thus alpha*Om>=S
pays the floor2^80. All constants and units remain unchanged.
The inner unit w_F=max(2^80,Ee+2lg(4n')+32,ceil(z/128))<=alpha*Om,
because its nonfloor terms are bounded by w_N, while
2^80<=S<=alpha*Om pays the new fixed floor, and
ceil(z/128)<=z/128+1<=alpha*Om follows from 80/128+eps_N<=alpha.
The selected inner heights7,7,2,2 obey their existing <=8 guard; NB
binary depths17,7 obey Ee+depth+6<=Ee+44. All generic component proofs apply.
For each integer slot ceil(A*w_F), use ceil(A*w_F)<=A*alpha*Om+1.
Thus m1<=D_F*Om and k_F<=K_F*Om, with
 D_F=alpha*d_inner+eps_N=286593853274667651243506652628649073664170408349904783083334691919/137006972430731656764180686500000000000000000000000000000000000,
 K_F=alpha*k_inner+eps_N=67166075557838086671173137770684533165414468874927623638194221709043/24910358623769392138941943000000000000000000000000000000000000.
Let (c_P,pi_P,Lam_P,g_P) be the exact PTM.1/PTM.2 for both arity vectors coefficients at vNB0,vNB2.
For Ext0 use the finite PTM.2 bound d0<=6c0*r+p0*w_N+D0,
 p0=7276956212161842440191/22879274691000000000, D0=0. Define
 M0=c2*D_F+pi2*alpha+eps_N,
 K2bar=Lam2*D_F+g2*alpha+eps_N,
 Rnew=(M0+p0*alpha+(11/2)*alpha+(53+D0)eps_N)/(theta-6c0)+eps_N,
 d0bar=6c0*Rnew+p0*alpha+(D0+1)eps_N,
 d1bar=2*(M0+d0bar+3alpha)+eps_N,
 DNB=2d0bar+2d1bar+D_F+K2bar+alpha.
All denominators are positive and p0>=16. Exact substitution gives
 ceil(M0)=84, ceil(K2bar)=6154,
 ceil(Rnew)=1350795, ceil(d0bar)=402,
 ceil(d1bar)=977, ceil(DNB)=11003.
Ee<w_N<=alpha*Om and 53+D0<=(53+D0)eps_N*Om give
(11/2)Ee+53+D0<=(11/2)*alpha*Om+(53+D0)eps_N*Om; the r formula has no advice term.
Advice guards. By (O2), e=3s_tag+r_F+42, so3a=3s_tag+27=Ee-(r_F+19) EXACTLY;
hence a<Ee/3. With sig=Ee+4, w_N>=Ee+44 and p0>=300 (6p0>=1800):
(AG1) wrappers(I1) gives s>=m0+6c0*r+p0*w_N+D0+3Ee+122/3>=3Ee+40>9a+40>=a+40,
 termwise for every S; at S=2^80, s/(a+40)>1470.
(AG2) eta*r>=6(theta-6c0)r>=6p0*w_N+33Ee+318>=1800w_N+318, while
 50(a+sig+2)<(50/3)Ee+50Ee+300<=(200/3)w_N+300<1800w_N+318, termwise for every S;
 at S=2^80, eta*r/(50(a+sig+2))>44.
Both guards hold at the only NB call (outer A#), so Theorem NB applies there.
The r ceiling costs <=eps_N*Om, so r_NB<=Rnew*Om. Both source branch coefficients are at most KC-eps_N,
 KC=18905931729686812016364677/1000000000000000000. These exact inequalities are recomputed from the displayed recurrences.
For AGHP, 2e_L<=d0+3Ee and ceil(2*(Rnew+2)*M0)+1<2^240 give
2lg(Rrow*m0+1)<=480+4lg Om<=eps_N*Om. The k80 base is
804<eps_N*2^80, with slope4<=804; doubling covers every Om>=2^80.
Therefore d1<=2m0+2d0+6Ee+eps_N*Om<=d1bar*Om and d_C<=DNB*Om.
Define U18(x)=ceil(10^18*x)/10^18 and Slot(x)=U18(x)+2eps_N.
Put SC=Slot(DNB)=11002278246233485278233/1000000000000000000, KC=Slot(max(K00bar,Kotherbar))=18905931729686812016364677/1000000000000000000,
where K00bar=6Lam0*Rnew+g0*alpha and Kotherbar=K_F+12Rnew+80+2alpha.
Since eps_N*Om>=1, floor(SC*Om)>=d_C and floor(KC*Om)>=k_C.
These are the actual padded seed and entropy slots used below.
The outer d3c contract ensures k_Cslot<=n'=d_O and d2<=d_C<=s_C<=n'.
The NB algorithm still has eta*r_NB>=6p0*w_N>=96*2^80 and (AG1),(AG2) above,
and r_NB>=6p0*w_N/eta supplies the separate ring-prime guard u>=2^90.
All base/hash degrees remain O(log actual n); large relative fields use towers.
## 6. Uniformity, exactness, and remaining obligations
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
Sections 8.2--8.3, affine-and-standard-wrappers.md Section 4.2, outer Section 3) and to the
zero-additive-shift independent-side use of NB's inner breaker
(affine-and-standard-wrappers.md Section 4.6, step 7). It does not assert a generic
kernel principle when arbitrary additive B-tuples correlate with suffixes.
# Tagged affine-permutation logarithmic-loss family: outer assembly
The affine-breaker contract (F') is proved below from the outer A# (Section 3) and Lemma
B9.
Credits: sumset framework, Chattopadhyay--Liao,
https://arxiv.org/abs/2110.12652v1 ; advice-condenser architecture, Xin Li,
https://eccc.weizmann.ac.il/report/2023/023/revision/1/ ; native full-B2 bridge,
minimal-circuit-bridge.md.
All logarithms are binary. Write lg(x)=ceil(log2 x) for positive x, computed
by integer bit length and a power-of-two comparison.
Foractualintegern>=1,flg(n)=floor(log2n),computedasbitlength(n)-1.
ONLYtheactualouterlengthusesflg.Everylgforfielddegrees,seeds,graphconstants
remainsCEILINGlog2.
## 1. The exact inputs to this assembly
Fix Q_S=3388944861463313528361373/1250000000000000000000000. With L=floor(log2(actual input length)), S=floor(L/Q_S)
is an integer; S<=L/Q_S<S+1 and L<Q_S(S+1). All other dimensions retain
their specified integer rounding.
(L_A) Put L=floor(log2 n), S=floor(L/Q_S)>=2^160, ell_w=lg(S+4), s=2ceil((beta*S-ell_w)/2),
ellc=lg(s+4) and r_F=S. Linear Section7f supplies the literal
three-condenser/kernel base E_n, source-linear at each seed and strong
for every affine source of dimension at least
 K_B=ceil(A_B*d+B_B*rho), rho=r_F+225+2lg S,
 d=floor((448929385877401036698799/20000000000000000)Om), with Om defined below and d<=3629829255L.
Its complete seed length D<=s, error<2^-r_F and public dimensions are
AF1. The full AF2--AF4 proof gives K_B<=1162477739635L<=HL.
(G_A) Sampler Sections3--5 use q_walk=19,C=80,m_walk=s+38r_F,
s_tag=3m_walk<=144S+960ellc,A=A_rows=2^s_tag and
Lambda_S=A^(1/2)/2^17. The full four-row AB+C incidence proof followed
by W_A expands every Lambda_S-row set to more than5*2^-r_F*2^s
neighbors. Direct affine averaging gives AS2: for every uniform affine
source of dimension>=K_B and fixed test of density<=2^-r_F, fewer than
Lambda_S rows are wholly bad except on probability<=2/5.
The evaluated sampler is source-linear and costs2^O(s_tag).
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
Write L=flg(n),S=floor(L/Q_S),ell_w=lg(S+4), s=2ceil((beta*S-ell_w)/2),ellc=lg(s+4),r_F=S,m_walk=s+38r_F,s_tag=3m_walk=log2 A,w=2. Require S>=2^160,L>=2048 and s_tag>=2L; the latter follows from the fixed constants below.
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
Use the source-linear AFFINE sampler Samp(x,alpha,j), j in [C], seed
length d, with AS2 at dimension k_S=K_B: for every uniform affine source
and fixed test of density<=epsilon0=2^-r_F, fewer than Lambda_S rows
are wholly bad except on probability<=2/5. Only uniform affine sampler
sources are used below. The enumerated field index alpha replaces [A].
For t=4C-1 and a=s_tag+w+lg(C), request the padded all-seeds-retained affine
breaker F of output length one, entropy K, and error at most
    epsilon = 2^(-r_F)/(2^22 C A^3).
The advice for (x,alpha,j) is the concatenation of the s_tag-bit alpha,
lg(C)-bit j-1, and w-bit H_alpha(x). Define
    R_alpha(x)=XOR_(j=1)^C F(x,Samp(x,alpha,j),advice(x,alpha,j)),
    D(x)=1 iff 2 sum_alpha R_alpha(x)>A.
Assume the input source entropy satisfies
    k >= max(k_S, K+4w+4Cd+2lg(1/epsilon)). (D1)
All concatenations, zero padding and majority ties are specified. There
are no seed or advice inputs to D. The tags are internal computations on
the actual n-bit x.
## D2. Mean error with arbitrary correlated left copies
Fix an arbitrary correlated tuple(A_0,...,A_r), fixed advice differing
from honest advice, and deterministic seed maps L_i of A_i. Assume
Hinf(A_0)>=K+ell+lg(1/epsilon), with joint seed leakage at most2^ell values.
For each seed y let f(y) be the MAXIMUM TV, over all finite choices of
b_0,...,b_r,y_1,...,y_r, of the honest bit in
 F(A_0+b_0,y+L_0(A_0),advice_0),
 (F(A_i+b_i,y_i+L_i(A_i),advice_i))_(i=1)^r
from a fresh uniform bit jointly with the other outputs. Then
 E_(uniform y) f(y)<=2epsilon. (D2mean)
Indeed choose a maximizing witness tuple for EVERY y; finite maxima exist.
Take Y uniform independent of the left tuple and let all b_i,y_i depend on
Y through those witnesses. Conditional on the complete L_i leakage, left
and right tuples are independent and honest seed Y+L_0(A_0) is uniform.
Bounded-range conditional entropy loses at most ell average bits. Its tail
at probability epsilon and the supplied conditional lift of(F') give TV
at most2epsilon, retaining all seeds and leakage. The honest seed and
leakage determine Y. Drop other retained data, then leakage: contraction
bounds TV retaining Y and tampered outputs. Since Y is independent of the
original left tuple, that TV is exactly the average of f(y). This proves
(D2mean). Markov also gives density<=2epsilon/theta for{f>theta}.
Witnesses and this function are finite proof objects, never evaluator searches.
## D3. Actual parities and pointwise collision counts
Fix independent uniform affine sources X_1,X_2 of dimension at least k, and write
A_0=X_1, A_1=Lambda X_1, B_0=X_2, B_1=Lambda X_2+b, with Lambda
invertible. Suppose tau(x)=Lambda x+b has no fixed point on the sum
support. Row alpha collides when H_alpha(x_1+x_2)=H_alpha(tau(x_1+x_2)).
The tag root bound gives at most N_col=257A/1024 collisions at EVERY point
(x_1,x_2). We use actual bits throughout; no random replacement is made.
Only parities of even size occur in the quartic and covariance identities below.
For each set Q of two or four distinct row-copy indices, expand
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
at most4Cd. Assign each Q its already chosen honest row-copy, independent of tag
profiles and of x_2. Set a_Q=1 for |Q|=4. For |Q|=2 set a_Q=A^2 if Q is a matched
pair {(0,alpha),(1,alpha)}, and a_Q=A otherwise. For each good left-tag cell,
right-tag profile with legal distinct honest advice, and honest neighbor
j, let f_i(y) be D2's maximum local distance for that conditional left law.
Omit colliding fully-doubled profiles; they will be charged as collisions.
Bad entropy cells have total mass<=epsilon per Q. The entropy hypothesis is
 k>=max(k_S,K+4w+4Cd+2lg(1/epsilon)). (D1')
Thus(D2mean) applies to EVERY included instance. Define
 Psi(y)=sum_i a_Q f_i(y).
There are C neighbor choices and at most2^(8w) profiles per Q. Since there are A matched pairs, at most2A^2 other pairs and at most(2/3)A^4 four-sets,
 sum_Q a_Q <= (2/3)A^4+3A^3 <= (2/3+3/256)A^4,
using A>=2^256. Hence E Psi<=2*2^16*(2/3+3/256) C A^4 epsilon (w=2).
Set B={y:Psi(y)>A/47} and epsilon=2^(-r_F)/(2^22 C A^3).
The exact inequality 47*2*(2/3+3/256)<64 and Markov give
 |B|/2^d<=2^22 C A^3 epsilon=2^(-r_F). (D2')
B depends on the left source law and tau, not the right source point;
the program never constructs B. Apply AS2 to BOTH B_0=X_2 and
B_1=Lambda X_2+b, uniform affine sources of the same dimension because
Lambda is invertible. Their total exceptional mass is<=4/5;
therefore good x_2 have mass at least1/5>0.
Fix such an x_2. Intersect the good row sets to obtain G, omitting
q<=2(Lambda_S-1)<2^(-16)A^(1/2) rows. Since 2^(-16)A^(-1/2)<=1/256,
g=|G|>=255A/256>=A/2. Each sigma,alpha in {0,1}xG has its own
neighbor j_(sigma,alpha) whose honest seed lies outside B. The two copies
need not choose the same j. Source-linearity gives every seed identity
with shifts B_sigma=x_2 or Lambda x_2+b and offsets Samp(B_sigma,alpha,j).
The second sampler is essential; its seeds need not be fixed translates
of the first sampler's seeds.
For EACH row-copy fix one such neighbor j_(sigma,alpha), once for all
Q assigned that honest row. Its offset seed y_row is fixed by x_2 and lies
outside B. Define b_Q as twice the sum of the relevant legal f_i(y_row)
over good left cells (dropping their probability weights), plus epsilon.
Actual parity bias is<=b_Q for a Q with a singleton; for fully doubled Q
it is<=b_Q+Pr[at least one selected collision]. Good profiles use the local
TV bound, hence parity bias<=2f_i; bad entropy cells cost<=epsilon and
collision profiles cost their probability. No individual small b_Q is claimed.
Group Q and profiles by their fixed honest row-copy. There are2g such
rows on G; each included instance occurs at most once in its row's Psi.
Consequently
 Bsum=sum_(Q on G,|Q| in {2,4}) a_Q b_Q
 <=2 sum_(row-copies on G) Psi(y_row)+32A^4 epsilon
 <=4gA/47+2^(-r_F)A/(2^17 C)
 <=gA*(4/47+2^-17)<3gA/35. (D3a)
Here g>=1,C>=1 and2^(-r_F)<=1. No independence of different offset
seeds is assumed. All expectations below are over X_1 at this fixed x_2.
The pointwise collision bound implies sum_(alpha in G) Pr[collision at alpha]
<=N_col. No independence of rows, copies, or collision events is asserted.
## D4. A robust quartic sign certificate
Let X_alpha,Y_alpha be actual signs on G and S_0=sum X_alpha,S_1=sum Y_alpha.
Put Z_alpha=X_alpha Y_alpha,Z=sum Z_alpha,zeta=E Z. Only matched pairs occur in zeta, each with weight A^2. Section D3 gives
 |zeta|<=N_col+Bsum/A^2<N_col+3g/(35A).
Since g>=255A/256 and A>=2^256,
 r=|zeta|/g<=257/1020+3/(35A)<13/50.
For Q_+(x,y)=xy(x-y)^2 and Q_-(x,y)=-xy(x+y)^2, expand the ACTUAL
ordered quartic monomials. All-even terms contribute-2g^2. A matched-pair
parity has coefficient plus or minus2(3g-2), and fully doubled two-row
parities contribute-8 sum_(alpha<beta) E Z_alpha Z_beta.
Every remaining parity has a singleton row and is bounded by b_Q.
A size-four parity has coefficient6,6,or-8 according to its two-copy
counts(3,1),(1,3),(2,2). For size two, a mixed nonmatched pair has absolute
coefficient6g-4; a same-copy pair has coefficient-4g. For a nonmatched pair its coefficient is at most6A and a_Q=A; for a four-set it is at most8 and a_Q=1. Thus every remaining absolute coefficient is at most8a_Q.
Thus the exact identities give
 E Q_+(S)=2(3g-2)zeta-2g^2+4g-4E Z^2+u_+,
 E Q_-(S)=-2(3g-2)zeta-2g^2+4g-4E Z^2+u_-, |u_+|,|u_-|<=8Bsum.
All doubled collision terms are retained exactly, not charged separately.
The variance identity E Z^2>=zeta^2 now bounds BOTH expectations by
 2(3g-2)|zeta|-2g^2+4g-4zeta^2+8Bsum
 <=[6r-4r^2-2+4/g+24A/(35g)]g^2
 <=[-444/625+1/1000+6144/8925]g^2<-g^2/50. (D4')
Here6r-4r^2 increases up to13/50, A/g<=256/255, and
4/g<=2^(3-2L)<1/1000 forL>=2048.
Marginal fourth moments have all-even part3g^2-2g. Remaining parity
coefficients are24 for size four and12g-16 for size two, each bounded by
24a_Q. Therefore
 E S_sigma^4<=3g^2+24Bsum
 <[3+72*256/(35*255)]g^2<16g^2.
Their L4 norms are consequently at most2sqrt(g). This weaker norm bound
is sufficient for the deterministic removed-row perturbation below.
Let Y_sigma be the ACTUAL sum over ALL rows in that copy. Pointwise
|Y_sigma-S_sigma|<=delta=q<2^(-16)A^(1/2). Since g>=A/2,
  delta/sqrt(g)<2^(-16)A^(1/2)/(A/2)^(1/2)=2^(-31/2)<1/40000.
The L4 norms of S_sigma are at most2sqrt(g), and those of Y_sigma at
most3sqrt(g). The norm triangle inequality follows from Holder applied
to (|X|+|Z|)|X+Z|^3 and division by the cube of the L4 norm (zero norm
is trivial). For four factors, Holder itself follows by Cauchy-Schwarz
twice: E|abcd|<=(E a^4 E b^4 E c^4 E d^4)^(1/4).
Telescoping the four factors in any degree-four monomial shows its
expectation changes by at most4delta(3sqrt(g))^3. Each Q has sum of
absolute coefficients four. Therefore
  |E Q_+(Y)-E Q_+(S)|, |E Q_-(Y)-E Q_-(S)|
    <=432delta g^(3/2)<(432/40000)g^2<g^2/80. (D5')
The expectation margins are -g^2/50 and the perturbations are less than g^2/80, so both actual expectations are strictly negative. If the actual majority
XOR were identically zero, Y_0Y_1>=0 pointwise and Q_+(Y)>=0, impossible.
If that XOR were identically one, Y_0Y_1<=0 and Q_-(Y)>=0, also impossible.
The weak product inequalities include all ties: the evaluator outputs one
exactly for a negative sign sum, and zero for a nonnegative sum.
Thus both actual majority-XOR values occur. No sixth or eighth moments,
small-tail inference, or independence of the full sign vector is used.
All sums and collision events above are actual deterministic quantities.
For an affine S=a+V of dimension at least k, take X_1 uniform on V and
X_2 independently uniform on a+V. Their sum lies in S. Good x_2 have
probability at least1/5. Fix any such x_2; the preceding contradiction
proves both values of D(x)+D(tau(x)) occur on S. Thus the property holds
for every invertible affine tau without fixed points on S, and includes
every nonzero translation.
## 2. Fix every constant before seeing n
Fix Q_S=3388944861463313528361373/1250000000000000000000000, beta=54/7-2^-100, q_walk=19, C=80, p=144,
ACTUAL outer degree t=4C-1=319 and j=lg(t+1)=lg(t)=9.
For actual n>=1 compute L=floor(log2 n), S=floor(L/Q_S),
 ell_w=lg(S+4), s=2ceil((beta*S-ell_w)/2)<=8S, 0<=delta_v=s-beta*S+ell_w<2,
 R(s)=0,
 ellc=lg(s+4), r_F=S,
 m_walk=s+2q_walk*r_F<=48S+320ellc,
 s_tag=3m_walk<=144S+960ellc, A=A_rows=2^s_tag,
 delta_O=3629829255; target d=floor((448929385877401036698799/20000000000000000)Om) after Om below.
Here A_rows (or A) is a row count, A_B=1/G_B is the BASE source
coefficient, gamma_B is the exact degree coefficient defined in linear7f,
and G_B is its rank slack. None denotes CZ's gamma or the graph Gamma.
Request the outer breaker at
 t=319, output m=1, a=s_tag+2+lg C, w=2,
 E=3s_tag+r_F+22+lg C=9s+343S+29<=433S+29. (O2)
Because lg C=7 and13*2^9-8=6648<2^13, its internal exponent is
 e=E+13<=433S+5768ellc+54,
 Om=e+2L+42<=433S+5768ellc+2L+96.
Thus 2^-E<=epsilon=2^-r_F/(2^22*C*A^3), including the ceiling in lg C.
Tags and seeds are internal. R=0<=8ellc, so all retained coarse upper calibrations remain valid. Sampler Section3 proves the walk count under the new S>=2^160 guard.
(O2a) For 2^k<=S<2^(k+1), k>=40, s+4<9S implies ellc<=k+5.
Use the common bound J=max(111+lg C+lg(26t),146+2lg C)=160.
At k40 the exact inequality is
 2*5768*45+160=519280
 <2^-20*(3388944861463313528361373/1250000000000000000000000)*2^40.
The left side grows by11536 when k grows by1; its positive intercept
makes this increase at most its value at k40, while the right side doubles.
Induction proves
 2*5768*ellc+J<2^-20*L
at EVERY allowed integer S, since L>=Q_S*S. This bound pays the logarithmic/fixed terms in 8+2E. For Om use
the separate tighter bound below. Here s<=8S, so the retained upper calibrations incur no positive rounding loss.
For Om alone start the dyadic induction at k=80. Its exact base is
 5768*85+96=490376 < 10^-15*(3388944861463313528361373/1250000000000000000000000)*2^80.
The left side increases by5768 per dyadic step, at most its positive
base; the right side doubles. Hence for EVERY integer S>=2^80,
 5768*ellc+96 < 10^-15*L.
This proves the sharper c_O below; the 8+2E allowance stays 2^-20.
The separate k80 bases 1067<10^-18*2^80 and817+Q_S<10^-18*2^80 prove
rho-S<eps_B*S and rho0-(Q_S+1)S<eps_B*S in linear7f, eps_B=10^-18.
The row-stretch base 48q_walk*45+6=41046<p*2^40 and the same induction
give 48q_walk*ellc+6<=pS. Hence
 s_tag<=2pS<=(2p/Q_S)L. (O2b)
This proves the runtime bound.
Constants precede n; each dimension uses its specified integer rounding.
## 3. The outer affine breaker A# and (F')
By(O2a), Om=e+2L+42<=c_O*L, where
 c_O=(3p+1)/Q_S+2+10^-15=548027889722926630445667607463313528361373/3388944861463313528361373000000000000000.
Construct wrappers2.1 DIRECTLY at t=319,j=9,C=NB, exponent e,output1
and advice a=s_tag+2+lg C. Use U18 and Slot from B9. All coefficients
are fixed before n. In the following definitions c_P,pi_P,Lam_P,g_P are the exact7g
coefficients, using pi_P=pi_*, and alpha is B9's unit bound:
 SC=11002278246233485278233/1000000000000000000, KC=18905931729686812016364677/1000000000000000000, S0=30949277744229673497/500000000000000000, Dslot=448929385877401036698799/20000000000000000, Kslot=18017882413063933821617/3125000000000000.
Equivalently SC=Slot(DNB), KC=Slot(max(K00bar,Kotherbar)),
 S0=Slot(c_L0*SC+pi_L0),
 TRbar=320*(S0+SC+18*(20+80)), TLbar=320*(SC+18*20),
 Dslot=Slot(max(320*(S0+SC)+KC+alpha,TRbar+Lam_E*80+g_E+319*80+alpha,
             S0,Lam_E*80+g_E,KC)),
 Kslot=Slot(max(TLbar+Lam_Lr*6423+g_Lr+319*6423+alpha,Lam_L0*SC+g_L0)).
Actual integer dimensions are
 dy=20Om, dx=80Om, r=6423Om,
 s_C=floor((11002278246233485278233/1000000000000000000)Om), s0=floor((30949277744229673497/500000000000000000)Om),
 k_Cslot=floor((18905931729686812016364677/1000000000000000000)Om), d_O=floor((448929385877401036698799/20000000000000000)Om), k_O=floor((18017882413063933821617/3125000000000000)Om).
Use PXG(vL0) for L0, PXG(vLr) for Lr, PXG(vE) for E and the supplied
uniform-seed AGHP extractor for Lm(input r,output dy, exponent e), padded
and retained to seed dx. B9 supplies NB's complete seed and both source branches.
This gives
 delta_O=ceil(448929385877401036698799/20000000000000000c_O)=3629829255,
 kappa_O=ceil(18017882413063933821617/3125000000000000c_O)=932377715.
Every source and input guard is checked separately. Since H>=delta_O
by(O3), d_O<=3629829255L<=HL<=n. Also r<=k_O<=kappa_O L<=HL<=n
(k_O>=r by k3b); wrappers2.1(3) has no k_O>=d_O requirement.
E's actual input is d_O; NB's actual source is n'=d_O; Lr,L0 read n.
For every primitive input u<=n<2^(L+1), lg(4u)<=L+3. Actual outer PXG
heights13,9,6 satisfy e+h+6+2lg(4u)<=e+2L+25<=Om.
PTM.1/PTM.2, PT1/PT2 and PXG include all their internal ceilings, and Om>=2^40.
PXG requires 16/15<=R_h<=64 and output<=256*2^h*Om: every selected
terminal rate and the outputs h,a,s (F#) and s_C,r,dx (outer) obey it.
For Lm, r*dy+1<2^33*Om^2 by 6423*20+1<2^33.
Hence lg(r*dy+1)<=33+2lg Om. The exact base390<10^-18*2^80,
slope4<=390 and doubling prove66+4lg Om<10^-18*Om for EVERY Om>=2^80.
Its exact seed2[lg(r*dy+1)+dy+2e] fits
 (2*20+4+10^-18)Om<=80Om,
and its threshold dy+2e<=(20+2)Om<=r. Every original seed and
padded suffix remains retained by the AGHP strongness proof.
The causal same-side ledgers of wrappers2.3 are
 TR=320(s0+s_C+18(dy+dx))<=12864176801721944625227/3125000000000000Om,
 TL=320(s_C+18dy)<=11362278246233485278233/3125000000000000Om;
the full transcript still retains the final (t+1)dx bits; no seed is dropped.
Let(c_P,pi_P,Lam_P,g_P) be the exact7g PXG coefficients, P in{E,Lr,L0}.
The finite row right sides are listed below. For a fractional slot A,
U18(rhs)+eps_N<=A implies floor(A*Om)>=rhs*Om because eps_N*Om>=1.
Integer-unit slots use rhs<=A directly. B9 gives q=e<=alpha*Om.
 dy_E: rhs <= 20. Here rhs=c_E*80+pi_E.
 dx_Lr: rhs <= 80. Here rhs=c_Lr*6423+pi_Lr.
 dx_Lm: rhs <= 80. Here rhs=2*20+4+10^-18.
 s0_L0: U18(rhs)+eps_N <= 30949277744229673497/500000000000000000. Here rhs=c_L0*SC+pi_L0.
 r_3a: rhs <= 6423. Here rhs=320*20+(20+2)+1.
 r_dy: rhs <= 6423. Here rhs=20.
 r_Km: rhs <= 6423. Here rhs=20+2.
 d_3c: U18(rhs)+eps_N <= 448929385877401036698799/20000000000000000. Here rhs=320*(S0+SC)+KC+alpha.
 d_3d: U18(rhs)+eps_N <= 448929385877401036698799/20000000000000000. Here rhs=TRbar+Lam_E*80+g_E+319*80+alpha.
 d_s0: U18(rhs)+eps_N <= 448929385877401036698799/20000000000000000. Here rhs=S0.
 d_KE: U18(rhs)+eps_N <= 448929385877401036698799/20000000000000000. Here rhs=Lam_E*80+g_E.
 d_kC: U18(rhs)+eps_N <= 448929385877401036698799/20000000000000000. Here rhs=KC.
 k_3b: U18(rhs)+eps_N <= 18017882413063933821617/3125000000000000. Here rhs=TLbar+Lam_Lr*6423+g_Lr+319*6423+alpha.
 k_K0: U18(rhs)+eps_N <= 18017882413063933821617/3125000000000000. Here rhs=Lam_L0*SC+g_L0.
The exact outer coefficients are defined by Section7g at vE,vLr,vL0.
For a second finite witness, verify.py rounds each coefficient X upward to
ceil(10^40 X)/10^40 and rechecks all12 inner and14 outer rows and B9.
The exact rows above and defining recurrences are authoritative; these
one-sided numerical enclosures are redundant arithmetic checks.
All floor slots include the explicit rounding buffer; all dimensions are fixed before source data, and d3c includes q=e. Together with B9 these cover every source,
seed, output, public cap, relative-field and padding guard. The proof of
wrappers2 gives error(13*2^9-8)2^-e=6648*2^-e<2^13*2^-e=2^-E,
with arbitrary within-side correlations and arbitrary additive right shifts.
All original d-bit seeds remain retained by the following consequence of
(F')'s universal right-side quantifiers.
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
WHOLE right side, including all padded seed suffixes, is retained free.
## 4. One fixed all-length guard and entropy coefficient
Put eps=2^-20, eps_B=10^-18, eps_R=2^-80 and fix
 T=ceil((6p+2)/Q_S+eps)=320,
 H_raw=4C*delta_O+kappa_O+T=1162477739635,
 H_base=ceil(A_B*delta_O+B_B*(1+eps_B)/Q_S+eps_R)=1162477739611,
 H=max(H_raw,H_base)=1162477739635. (O3)
Here H>=delta_O, which supplies the separate d_O<=n guard.
By (O2a), 4w+2lg(1/epsilon)<=8+2E
 <=866S+11536ellc+88<=T*L.
Therefore k_F+4w+4Cd+2lg(1/epsilon)
 <=[kappa_O+4C*delta_O+T]L<=H*L.
AF4 separately proves k_S=K_B<=H_base*L<=H*L. The final eps_R in
H_base pays the K_B ceiling. Thus the FULL conservative (D1) holds at
k=HL, retaining the tag cost, all4C seed leakages and BOTH entropy tails.
For integer floors S>=L/Q_S-1 and 135>3Q_S, L>=p+128 imply
 135S>=135L/Q_S-135>3L-135>=2L.
Since linear7f proves s>=7S, s_tag=3(s+38S)>=135S>2L, and every ordered pair of L-bit tag blocks has
exactly A/2^(2L) preimages among the COMPLETE row indices. The tag
collision proof uses A itself, with no stronger lower power of n.
Even s,m_walk make Lambda_S=2^(3m_walk/2-17) an integer. Also
m_walk>=45*2^160>134217728, paying the incidence characteristic guard.
The affine-only sampler is used only on uniform affine B_0 and B_1.
The total evaluator tests these literal finite guards:
 n>=1; L=floor(log2 n)>=max(2048,p+128);
 S=floor(L/Q_S)>=2^160; s_tag>=2L; HL<=n. (5'')
It also computes/tests every actual AF1 and primitive length/seed/entropy
parameter guard, returning zero on any failed public check before setup.
AF1--AF4, B9 and all fourteen outer requirements prove that those guards
pass under (5'').
Since S tends to infinity at EVERY integer actual length, all fixed
lower thresholds eventually hold. Finally n>=2^L and2^L/L tends to
infinity, its consecutive ratio being2L/(L+1)>=4/3 for L>=2.
At L=128, H*128<2^128 by exact arithmetic. The same increasing ratio
proves HL<=2^L<=n at EVERY L>=128. Thus S>=2^80 already supplies
this native guard, and gives ONE common fixed threshold N for ALL
integer n>=N. All exceptional lengths, including n=0, use the zero branch.
## 5. The fixed evaluator
Precompute the constants in Section 2 by their finite algorithms. On an
input x of actual length n:
1. If n=0, output0. Compute L=floor(log2 n), S=floor(L/Q_S).
   If L<max(2048,p+128) or S<2^160, output0 before any logarithm of S.
   Otherwise compute the remaining public parameters. If s_tag<2L,
   HL>n or any exact supplied parameter guard fails, output0.
2. Build the exact base extractor E_n from (L_A), seed-padded to s bits.
   Its data argument is x and its output length is d.
3. Build the (G_A) graph on left s_tag-bit strings: the four rows of the AB+C
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
F_n is the raw d_F-seed breaker, identical to its padded form.
The affine base sampler tests density epsilon0=2^-r_F. By AS2, for every
uniform affine source of dimension>=K_B and every fixed test of that
density, fewer than Lambda_S rows are wholly bad except on mass<=2/5.
Source-linearity survives graph index selection. The degree and advice lengths match Section D1, and (D1) was checked
above. Sections D2--D4 prove the affine-permutation property at k=HL,
and in particular the directional property.
## 6. Uniform polynomial time in actual input length
One-time constant computation terminates independently of n; parameter
arithmetic adds O(log n) bits. There are
 A*C<=C*n^ceil(2p/Q_S)=C*n^107
breaker calls. Affine base degrees b0,b=O(log n) permit2^O(log n) canonical
searches; towers build its initial O(n)-dimensional extension. All later
ambient/matrix/kernel dimensions are O(log n), so unit-vector evaluation and
Gaussian elimination are polynomial. The graph prime p_m has m+1 bits,
m=s_tag/3, and is found in2^O(s_tag); evaluating all A rows and the walk
therefore costs polynomial time in actual n.
By(F'), each breaker call is polynomial in n for fixed t=4C-1,a,E=O(log n).
Tags cost2^O(L); all A evaluations are polynomial in n.
Majority uses an O(log A)-bit counter. The total evaluator is
polynomial in ACTUAL n.
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
f(x)=x satisfies separation and has zero gates.
For the fixed evaluator, at EVERY actual integer n>=N,
 C_fullB2(f_n)>=4n-4H*flg(n)>=4n-Klog2n,K=4H=4649910958540. (O4)
Here flg(n)<=log2n, N is the common guard threshold, and the additive
remainder is ZERO. BASE K_B and primitive K0 are entropy requirements.
All constants and exceptional zero branches define one fixed total P
evaluator, with f_0=0. The NP verifier accepts (x,empty) exactly when it
outputs1, and rejects nonempty witnesses (witness bound zero). Tags, seeds
and fields are internal. The remainder is exactlyKlog2n; log2n/n tends to
zero, giving both joint projections(4,1).
The actual stretched row length obeys(O2b):
48q_walk*ellc+6<=pS and s_tag<=2pS<=(2p/Q_S)L. The exact k40
comparison41046<144*2^40 and doubling prove this at every integer S.
Thus every2^O(s_tag) enumeration is polynomial in actual n, including
all row, prime and field work.
# Logarithmic width trim eliminates the constant loss
The construction above uses beta=54/7-2^-100, ell_w=lg(S+4), s=2ceil((beta*S-ell_w)/2), r_F=S,
and the strengthened S>=2^160 guard. BASE's original field-degree
coefficients, gamma_B, all exact slots and Q_S remain fixed.
Put L=floor(log2 n), h=log2 L, a=A_B, b=B_B, Q=Q_S,
u=448929385877401036698799/20000000000000000,
v=18017882413063933821617/3125000000000000, J=320u+v.
The new width satisfies s<beta*S-ell_w+2, with ell_w=lg(S+4)>=160.
Thus E<=E1*S-9ell_w+47, E1=2887/7-9*2^-100, and
Om<=E1*S+2L-9ell_w+102. Also rho<=L/Q+2h+227.
Define A=J*(E1/Q+2)+2E1/Q, D=9(J+2), C=102J+102.
Then 320d+k_O+8+2E<=AL-D*ell_w+C.
Since ell_w>=160 and D*160-C=1338J+2778>2, this is at most AL-2.
For BASE put A0=a*u*(E1/Q+2)+b/Q, B0=2b, C0=102a*u+227b.
Dropping the negative width term gives a*d+b*rho<=A0*L+B0*h+C0.
Exact arithmetic gives Delta=A-A0>0. On 2^j<=L<2^(j+1),
h/L<=(j+1)/2^j tends to zero, with successive ratio at most3/4 for j>=1.
Therefore Delta*L>=B0*h+C0+3 at every sufficiently large L.
Consequently K_B=ceil(a*d+b*rho)<=AL-2 eventually.
Set the proof-only integer k=ceil(AL)-1. It is at least AL-1, hence
satisfies BOTH entropy requirements, and it is at most AL.
For all large L, k>=1 and k<=AL<=HL<=n. No construction parameter depends on k.
The native bridge yields C_B2(f_n)>=4n-4k>=4n-4AL>=4n-4A*log2(n).
The exact log-loss coefficient is c1=24823644797417699755355511905034617540900032714711932563336512548891785726331/5601343869123682735034745026567550992753937612800000000000000000. No constant loss is declared.
The BASE deficit and all finite guards are absorbed by one fixed eventual
threshold for every integer n beyond it, not a sparse subsequence.
