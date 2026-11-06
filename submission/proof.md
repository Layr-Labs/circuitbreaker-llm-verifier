# Fixed logarithmic-loss family and native full-B2 theorem

This file and its four included companions define ONE fixed total language
in P, with the empty-witness NP verifier below, at every actual integer
input length. For every sufficiently large integer n, the complete theorem is
 C_fullB2(f_n)>=4n-4H*floor(log2 n)>=4n-250489174635004*log2 n,
 H=62622293658751, with ZERO additive loss and joint projection(4,1).
The model has all16 binary tables, sharing, repeated pins and unrestricted
fan-out; all retained gates count once. The construction and proofs follow.
All arguments are supplied below.

bound.json declares the theorem proved here.

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

1. PXT/PXG7g and exact fields feed restricted Sections1--3 and8.2.
2. Restricted8.2's five-round antichain14-copy theorem feeds wrappers4.2's F#.
3. CZ/weak-AGHP, AGHP, Ext0,Ext2,F# feed wrappers4's NB.
4. NB with B9's complete seed/threshold feeds outer3's A# at the sampler degree.

Constants use exact integer/rational operations or proved terminating searches.
For positive rational x, lg(x) is the least nonnegative integer j with2^j>=x.

## 1. Inputs and their precise scope

All inputs are included. Linear Sections2--4,7a,7g prove the arbitrary-source
PXT/PXG vectors and exact fields; the same file supplies bilinear/weak-seed/AGHP.
Only BASE uses affine Section7f.
Sampler Part SP proves finite rings and complete incidence: gamma=256/20481,
prime guard2^3200000. Restricted supplies prime/encoding, CZ smoothing,
conditional tools and the five-round antichain8.2 theorem retaining all15
seeds. Wrappers4.2 instantiates its zero-shift14-copy F# in NB4.6 step7.
Wrappers2,4 prove the all-seeds affine transformation and NB. Outer3 builds
A# at t335, B9 pays its full NB call. Sampler3--6 prove W_A/G_A/AS2 and the
general conditional bad-seed lemmas. Sampler0 resolves legacy
component names to these included bodies. Every component seed is retained.

## 2. Primitive convention

Selected slots use the proved PXT/PXG7g primitives;
explicitly named AGHP calls use the separately proved matrix extractor. Each has
its explicit source, seed and dimension guards.

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
F# is the map of wrappers4.2. Failed public guards return the all-zero
string; every guard holds at every sufficiently large actual length.

## 4. Fix the uniform-source advice condenser first

Fixgamma=256/20481,eta=gamma/(80+7*2^-20)=268435456/1718070947847,theta=eta/7=268435456/12026496634929.
Every actual NB call haseta*r>=112*2^20,as proved
in wrappers(I3),so this guard and the independent ring-prime guard hold.
NB's repaired honest Ext0 output is exactly uniform, so it uses Lemma U
of affine-and-standard-wrappers.md Section3.5.
Put sig=Ee+4 and s=floor((eta*r-17sig-8)/7). NB's finite parameter
algorithm proves all Lemma U guards, s>=m0+2e_L and s>=a+40.
Its advice row error is at most14*2^-sig<2^-Ee.

## 5. The standard breaker NB and Lemma B9

Theorem NB (affine-and-standard-wrappers.md Section 4) is a standard breaker
in the sense of Section 3.1, with guards k_C<=n and d2<=n, at every a,z,E>=1.
It uses F#, the Section4 advice condenser, AGHP, Ext0 and Ext2; the outer
assembly uses NB only through A#, and Lemma B9 gives its seed and threshold.

**Lemma B9 (NB at the ACTUAL outer call).** Source length n'=d_O<=n,
output z=128Om, advice a=s_tag+2+lg C<=Om, exponent e and Om=e+2L+42.
The exact NB algorithm satisfies d_C<=40059Om and k_C<=627773399Om.
Use eps_N=10^-6. On2^h<=Om<2^(h+1), h>=40, lg Om<=h+1. The exact base
8*41+240=568<eps_N*2^40 and 8<=568 imply by doubling induction
8lg Om+240<eps_N*Om for EVERY integer Om>=2^40.

The NB unit is w_N=max(2^20,Ee+44+2lg(4n'))<=Om+12<=(1+eps_N)Om,
where Ee=e+4. The ACTUAL Ext0,Ext2 heights are21,8: e_PXT+h+6<=Ee+27<=Ee+44.
The inner F# unit is w_F=max(2^20,Ee+2lg(4n')+32,ceil(z/128))<=Om,
since lg(4n')<=L+3 and z=128Om<=128Om. Its selected heights7,8,2,2
satisfy (Ee+7)+h+6<=Ee+21<=Ee+32. All12 ceiling-charged requirements
in wrappers4.2 hold for the supported selective75a/44h+45s/11s theorem. Its
public cap covers z=128Om.
Including each separate outer ceiling gives
 D_F=d_inner+eps_N=9761830863471/4096000000,
 K_F=k_inner+eps_N=2119802542774625759/16384000000.
Thus m1<=D_F*Om and k_F<=K_F*Om. Every slot is ceil(A*w_F), with
ceil(A*w_F)<=A*Om+1<=(A+eps_N)Om.

For P in {Ext0,Ext2}, let(c_P,pi_P,Lam_P,g_P) be the EXACT PT1/PT2
PXT recurrences of7g using the full selected vector vNB0 or vNB2. Put
 M0=c2*D_F+pi2*(1+eps_N)+eps_N,
 K2bar=Lam2*D_F+g2*(1+eps_N)+eps_N,
 Rnew=(M0+pi0*(1+eps_N)+7+7eps_N)/(theta-7c0)+eps_N,
 d0bar=7c0*Rnew+pi0*(1+eps_N)+eps_N,
 d1bar=2*(M0+d0bar+3)+1,
 DNB=2d0bar+2d1bar+D_F+K2bar+1.
The exact rational comparisons give
 ceil(M0)=64, ceil(K2bar)=33062,
 ceil(Rnew)=35599135, ceil(d0bar)=725,
 ceil(d1bar)=1583, ceil(DNB)=40059.
The positivity theta-7c0>0 and pi0>=16 follow by exact substitution.
Since(38/7)Ee+a+53<=7Om (38/7+1+53/2^40<7), NB's own r ceiling
has r_NB<=Rnew*Om and its seed has d0<=d0bar*Om. Both source branches obey
 ceil(7Lam0*Rnew+g0*(1+eps_N))=627773399<=627773399,
 ceil(K_F+14Rnew+128+2+eps_N)=627770494<=627773399.
The second branch includes2v=14r_NB, the FULL z=128Om and2Ee<=2Om.

For the actual AGHP degree,2e_L<=d0+3Ee and
ceil(Rnew*M0)+1<2^240. Thus lg(r_NB*m0+1)<=240+2lg Om.
The same doubling induction gives240+2lg Om<Om/4 at every Om>=2^40;
its k40 bound is322<2^38 and its affine slope2 is at most its base322.
Consequently d1=2lg(r_NB*m0+1)+2m0+4e_L<=d1bar*Om.
It follows exactly that d_C=2d0+2d1+m1+K2+Ee<=DNB*Om<=40059Om.
Taking the ceiling of the MAXIMUM of the TWO source branches gives
k_C<=627773399Om. The outer source n'=d_O has n'>=k_Cslot and
d2<=d_C<=s_C<=n', by the actual outer d3c contract. Its CZ guard:
eta*r_NB>=7*pi0*w_N>=112*2^20>64*2^20; the independent prime/ring
u>=3200000 follows separately from r_NB>=7*pi0*w_N/eta.
All local base/hash degrees are O(log actual n); larger relative fields use
explicit towers.

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
Section 8.2, affine-and-standard-wrappers.md Section 4.2, outer Section 3) and to the
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


Fix Q_S=93656792881/31250000000. With L=floor(log2(actual input length)), S=floor(L/Q_S)
is an integer; S<=L/Q_S<S+1 and L<Q_S(S+1). All other dimensions retain
their specified integer ceilings.


(L_A) Put L=floor(log2 n), S=floor(L/Q_S)>=2^40, s=8S,
ellc=lg(s+255) and r_F=S+8ellc. Linear Section7f supplies the literal
three-condenser/kernel base E_n, source-linear at each seed and strong
for every affine source of dimension at least
 K_B=ceil(A_B*d+B_B*rho), rho=r_F+185+2lg S,
 d=186361825982L.
Its complete seed length D<=s, error<2^-r_F and public dimensions are
AF1. The full AF2--AF4 proof gives K_B<=62622293611440L<=HL.

(G_A) Sampler Sections3--5 use q_walk=20,C=84,m_walk=s+40r_F,
s_tag=3m_walk=144S+960ellc,A=A_rows=2^s_tag and
Lambda_S=A^(1/2)/2^15. The full four-row AB+C incidence proof followed
by W_A expands every Lambda_S-row set to more than256*2^-r_F*2^s
neighbors. Direct affine averaging gives AS2: for every uniform affine
source of dimension>=K_B and fixed test of density<=2^-r_F, fewer than
Lambda_S rows are wholly bad except on probability<=1/128<1/100.
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

Write L=flg(n),S=floor(L/Q_S),s=8S,ellc=lg(s+255),r_F=S+8ellc,m_walk=s+40r_F,s_tag=3m_walk=log2 A,w=2. Require S>=2^40,L>=2048 and s_tag>=2L; the latter follows from the fixed constants below.
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
are wholly bad except on probability<=1/128. Only uniform affine sampler
sources are used below. The enumerated field index alpha replaces [A].
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

Fix independent uniform affine sources X_1,X_2 of dimension at least k, and write
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
the program never constructs B. Apply AS2 to BOTH B_0=X_2 and
B_1=Lambda X_2+b, uniform affine sources of the same dimension because
Lambda is invertible. Their total exceptional mass is<=1/64<2/100;
thus the original conservative98/100 good-mass bound holds.
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

Fix Q_S=93656792881/31250000000, beta=8, q_walk=20, C=84, p=144,
ACTUAL outer degree t=4C-1=335 and j=lg(t+1)=lg(t)=9.
For actual n>=1 compute L=floor(log2 n), S=floor(L/Q_S),
 s=2ceil(beta*S/2)=8S, delta_v=s-beta*S=0,
 ellc=lg(s+255), r_F=S+8ellc,
 m_walk=s+2q_walk*r_F=48S+320ellc,
 s_tag=3m_walk=144S+960ellc, A=A_rows=2^s_tag,
 d=delta_O*L, delta_O=186361825982.
Here A_rows (or A) is a row count, A_B=1/G_B is the BASE source
coefficient, gamma_B is the exact degree coefficient defined in linear7f,
and G_B is its rank slack. None denotes CZ's gamma or the graph Gamma.
Request the outer breaker at
 t=335, output m=1, a=s_tag+2+lg C, w=2,
 E=6s_tag+r_F+33+lg C=865S+5768ellc+40. (O2)
Because lg C=7 and lg(26t)=14, its internal exponent is
 e=E+14=865S+5768ellc+54,
 Om=e+2L+42=865S+5768ellc+2L+96.
Thus 2^-E<=epsilon=2^-r_F/(2^33*C*A^6), including the ceiling in lg C.
Tags and seeds are internal.

(O2a) For 2^k<=S<2^(k+1), k>=40, s+255<9S implies ellc<=k+5.
Use the common bound J=max(111+lg C+lg(26t),146+2lg C)=160.
At k40 the exact inequality is
 2*(288q_walk+8)*45+160=519280
 <2^-20*(93656792881/31250000000)*2^40.
The left side grows by11536 when k grows by1; its positive intercept
makes this increase at most its value at k40, while the right side doubles.
Induction proves
 2*(288q_walk+8)*ellc+J<2^-20*L
at EVERY allowed integer S, since L>=Q_S*S. This single bound pays both
the logarithmic/fixed terms in Om and those in 8+2E. In the general even
rounding formula it also pays 18delta_v<36 in Om and36delta_v<72 in
8+2E; here both losses are zero, rather than silently discarded.
The separate exact bases 627<628 and457+Q_S<628 similarly prove
rho-S<eps_B*S and rho0-(Q_S+1)S<eps_B*S in linear7f, eps_B=157/2^38.

The row-stretch base 48q_walk*45+6=43206<p*2^40 and the same induction
give 48q_walk*ellc+6<=pS. Hence
 s_tag<=2pS<=(2p/Q_S)L. (O2b)
This proves the runtime bound.
Constants precede n; every literal dimension uses its integer ceiling.

## 3. The outer affine breaker A# and (F')

By(O2a), Om=e+2L+42<=c_O*L, where
 c_O=(6p+1)/Q_S+2+2^-20=28540732624160767793/98206265251987456.
Construct wrappers2.1 DIRECTLY at t=335,j=9,C=NB, exponent e,output1
and advice a=s_tag+2+lg C. The actual integer slots, all multiples of Om, are
 dy=23Om, dx=128Om, r=7754Om, s_C=40059Om, s0=66Om,
 k_Cslot=627773399Om, d_O=641255400Om, k_O=16241565Om.
Use PXG(vL0) for L0, PXG(vLr) for Lr, PXG(vE) for E and the supplied
uniform-seed AGHP extractor for Lm(input r,output dy, exponent e), padded
and retained to seed dx. B9 supplies NB's complete seed and both source branches.
This gives
 delta_O=ceil(641255400c_O)=186361825982,
 kappa_O=ceil(16241565c_O)=4720128221.

Every source and input guard is checked separately. Since H>=delta_O
by(O3), d_O<=186361825982L<=HL<=n. Also r<=k_O<=kappa_O L<=HL<=n
(k_O>=r by k3b); wrappers2.1(3) has no k_O>=d_O requirement.
E's actual input is d_O; NB's actual source is n'=d_O; Lr,L0 read n.
For every primitive input u<=n<2^(L+1), lg(4u)<=L+3. Actual outer PXG
heights13,10,7 satisfy e+h+6+2lg(4u)<=e+2L+25<=Om.
PT1/PT2 and PXG include all their internal ceilings, and Om>=2^40.
PXG requires 5/3<=R_h<=64 and output<=256*2^h*Om: every selected
terminal rate and the outputs h,a,s (F#) and s_C,r,dx (outer) obey it.

For Lm, r*dy+1<2^33*Om^2 by 7754*23+1<2^33.
Hence lg(r*dy+1)<=33+2lg Om. The exact base230<10^-6*2^40,
slope4<=230 and doubling prove66+4lg Om<10^-6*Om for EVERY Om>=2^40.
Its exact seed2[lg(r*dy+1)+dy+2e] fits
 (2*23+4+10^-6)Om<=128Om,
and its threshold dy+2e<=(23+2)Om<=r. Every original seed and
padded suffix remains retained by the AGHP strongness proof.

The causal same-side ledgers of wrappers2.3 are
 TR=336(s0+s_C+18(dy+dx))=14395248Om,
 TL=336(s_C+18dy)=13598928Om;
the full transcript still retains the final (t+1)dx bits; no seed is dropped.
Let(c_P,pi_P,Lam_P,g_P) be the exact7g PXG coefficients, P in{E,Lr,L0}.
All14 actual outer contracts, with each entry a coefficient of Om, are:
 dyE: c_E*128+pi_E = 1873801744661/81920000000 <= 23.
 dxLr: c_Lr*7754+pi_Lr = 524260011067/4096000000 <= 128.
 dxLm: 2*23+4+10^-6 = 50000001/1000000 <= 128.
 s0L0: c_L0*40059+pi_L0 = 5398344931107/81920000000 <= 66.
 r3a: 336*23+(23+2)+1 = 7754 <= 7754.
 rdy: 23 = 23 <= 7754.
 rKm: 23+2 = 25 <= 7754.
 d3c: 336*(66+40059)+627773399+1 = 641255400 <= 641255400.
 d3d: 14395248+Lam_E*128+g_E+335*128+1 = 3999675062474462620533515527204102346840388588870317981/20893817942201854053401252024079286059435000000 <= 641255400.
 ds0: 66 = 66 <= 641255400.
 dKE: Lam_E*128+g_E = 3698007423722437707671335361718934508486364391755317981/20893817942201854053401252024079286059435000000 <= 641255400.
 dkC: 627773399 = 627773399 <= 641255400.
 k3b: 13598928+Lam_Lr*7754+g_Lr+335*7754+1 = 223800384504610220967142956193786819127016176247870826983555436728163166669472983/13779484203925446582049200206934514186682087304576545318878148468372500000 <= 16241565.
 kK0: Lam_L0*40059+g_L0 = 17895328351185854991553481359154609775430728196297718808623222730719373669041880327972178193563/1128673682381429960532940215543555567962781041816663341970004289690991410444122020000000 <= 16241565.

The exact ordered-rate coefficients are:
 E: c=1/64, pi=1709961744661/81920000000, lam=383083890765061901644121019259633079301820629/1392921196146790270226750134938619070629, g=2962486353453518856514623004740438996226868784075317981/20893817942201854053401252024079286059435000000.
 Lr: c=1/512, pi=462228011067/4096000000, lam=16670719523435709643360183324952805096690379602654586762385146637349/5511793681570178632819680082773805674672834921830618127551259387349, g=297544887570049586061488953613818219850859251040564607595610307782060976972983/13779484203925446582049200206934514186682087304576545318878148468372500000.
 L0: c=1/4096, pi=4597164931107/81920000000, lam=8738285431844011242294546190055939749528597292084775036978496211627004554161989/43112058150551182602480527713657584719739535592691495109625832302940848374489, g=8731098516515071521977559586991945464068464250077242404681723773085171196148995755722998193563/1128673682381429960532940215543555567962781041816663341970004289690991410444122020000000.

### Finite arithmetic record (AE)

Put D=10^30. For each exact coefficient X=(c,pi,Lam,g) of Section7g,
the corresponding integer U below satisfies 0<=U-D*X<1. Thus X<=U/D
(substitute in7g and clear positive denominators).

|P|U_c|U_pi|U_Lam|U_g|
|--|--:|--:|--:|--:|
|MI|15625000000000000000000000000|19855104953381347656250000000000|161031740811401330451388976743295517|73037684530187750234843083256621119076|
|MF|7812500000000000000000000000|22631303007055664062500000000000|103418831453529232817521350859207107|93173177675197125697976376983119065279|
|MG|500000000000000000000000000000|11666844485009765625000000000000|2391681944430931694660141602800|39394028558175927226513284049993|
|MJ|500000000000000000000000000000|13190239285009765625000000000000|1920833485722495524047743386509|34285567035651210055726935960821|
|NB0|530616196700178120301327|592136594532380340570963654772561|2256729554285217511381208971049|65409998795809916859882020113887618236|
|NB2|5659213360720997139859371165|49733715511033725803411997054618|8672210208530838670164602160180|12393534922848505839074369021832717|
|L0|244140625000000000000000000|56117735975427246093750000000000|202687735327530246635785669917922|7735715515305545910904398963656310660|
|Lr|1953125000000000000000000000|112848635514404296875000000000000|3024554344110830185943283670239|21593325495107149243049659348726672|
|E|15625000000000000000000000000|20873556453381347656250000000000|275021940813865938639768072128807697|141787698239191369082266729293421446702|

Each of the12 inner and14 outer RHSs is nondecreasing in its coefficients;
substituting U/D leaves each <=its slot. In B9 substitute successively in M0,
K2bar,Rnew,d0bar,d1bar,DNB and both source branches. The denominator
theta-7U_c(NB0)/D is positive and decreasing in c0. The successive upper ceilings are
64,33062,35599135,725,1583,40059,627773399,627770494, respectively.
The exact equality rows r3a,d3c
remain integer identities. B9 and7g prove the remaining parameter guards.

All slots are integers fixed before source data; d3c includes q_A. Together with B9 these cover every source,
seed, output, public cap, relative-field and padding guard. The proof of
wrappers2 gives error(13*2^9-8)2^-e=6648*2^-e<8710*2^-e=26t*2^-e<=2^-E,
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

Put eps=2^-20, eps_B=157/2^38, eps_R=2^-40 and fix
 T=ceil((12p+2)/Q_S+eps)=578,
 H_raw=4C*delta_O+kappa_O+T=62622293658751,
 H_base=ceil(A_B*delta_O+B_B*(1+eps_B)/Q_S+eps_R)=62622293611440,
 H=max(H_raw,H_base)=62622293658751. (O3)
Here H>=delta_O, which supplies the separate d_O<=n guard.
By (O2a), 4w+2lg(1/epsilon)<=8+2E
 =1730S+11536ellc+88<=T*L.
Therefore k_F+4w+4Cd+2lg(1/epsilon)
 <=[kappa_O+4C*delta_O+T]L<=H*L.
AF4 separately proves k_S=K_B<=H_base*L<=H*L. The final eps_R in
H_base pays the K_B ceiling. Thus the FULL conservative (D1) holds at
k=HL, retaining the tag cost, all4C seed leakages and BOTH entropy tails.

For integer floors S>=L/Q_S-1 and p>3Q_S, L>=p+128 imply
 pS>=pL/Q_S-p>3L-p>=2L.
Consequently s_tag>=pS>2L, and every ordered pair of L-bit tag blocks has
exactly A/2^(2L) preimages among the COMPLETE row indices. The tag
collision proof uses A itself, with no stronger lower power of n.
Even s,m_walk make Lambda_S=2^(3m_walk/2-15) an integer. Also
m_walk>=48*2^40>134217728, paying the incidence characteristic guard.
The affine-only sampler is used only on uniform affine B_0 and B_1.

The total evaluator tests these literal finite guards:
 n>=1; L=floor(log2 n)>=max(2048,p+128);
 S=floor(L/Q_S)>=2^40; s_tag>=2L; HL<=n. (5'')
It also computes/tests every actual AF1 and primitive length/seed/entropy
parameter guard, returning zero on any failed public check before setup.
AF1--AF4, B9 and all fourteen outer requirements prove that those guards
pass under (5'').

Since S tends to infinity at EVERY integer actual length, all fixed
lower thresholds eventually hold. Finally n>=2^L and2^L/L tends to
infinity, its consecutive ratio being2L/(L+1)>=4/3 for L>=2.
At L=128, H*128<2^128 by exact arithmetic. The same increasing ratio
proves HL<=2^L<=n at EVERY L>=128. Thus S>=2^40 already supplies
this native guard, and gives ONE common fixed threshold N for ALL
integer n>=N. All exceptional lengths, including n=0, use the zero branch.

## 5. The fixed evaluator

Precompute the constants in Section 2 by their finite algorithms. On an
input x of actual length n:

1. If n=0, output0. Compute L=floor(log2 n), S=floor(L/Q_S).
   If L<max(2048,p+128) or S<2^40, output0 before any logarithm of S.
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
density, fewer than Lambda_S rows are wholly bad except on mass<=1/128.
Source-linearity survives graph index selection. The degree and advice lengths match Section D1, and (D1) was checked
above. Sections D2--D4 prove the affine-permutation property at k=HL,
and in particular the directional property.

## 6. Uniform polynomial time in actual input length

One-time constant computation terminates independently of n; parameter
arithmetic adds O(log n) bits. There are
 A*C<=C*n^ceil(2p/Q_S)=C*n^97
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
 C_fullB2(f_n)>=4n-4H*flg(n)>=4n-Klog2n,K=4H=250489174635004. (O4)
Here flg(n)<=log2n, N is the common guard threshold, and the additive
remainder is ZERO. BASE K_B and primitive K0 are entropy requirements.
All constants and exceptional zero branches define one fixed total P
evaluator, with f_0=0. The NP verifier accepts (x,empty) exactly when it
outputs1, and rejects nonempty witnesses (witness bound zero). Tags, seeds
and fields are internal. The remainder is exactlyKlog2n; log2n/n tends to
zero, giving both joint projections(4,1).

The actual stretched row length obeys(O2b):
48q_walk*ellc+6<=pS and s_tag<=2pS<=(2p/Q_S)L. The exact k40
comparison43206<144*2^40 and doubling prove this at every integer S.
Thus every2^O(s_tag) enumeration is polynomial in actual n, including
all row, prime and field work.
