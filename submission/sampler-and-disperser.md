# Finite support graph, linear sampler, and conditional breaking

Fix Q_S=293/2 before input length. Every sampler/base unit is the integer
S=floor(L/Q_S)=floor(2L/293),L=floor(log2(actual input length)).
All field/seed dimensions use integer ceilings. Even for rational Q_S,
S<=L/Q_S<S+1 and L<Q_S(S+1).

All logarithms are binary, addition is XOR unless a field is specified, and
statistical distance is TV = half the L1 distance. A k-source is a finite
random variable with every atom at most 2^(-k); k may be real. Graphs have
labelled edges and may have loops or repeated edges. Every constant is fixed
independently of the eventual outer input length N.

The support-growth Lemma S1' (Section SP.6.2), the generalized
seeded-linear extractor (Theorem D3) and the affine breaker yield the finite
graph/sampler and conditional-breaking reductions below. proof.md supplies
the tagged affine-permutation outer assembly and finite guards.

## 1. Supplied contracts and elementary probability facts

### 1.1 Finite seeded-linear contract

The supplied seeded-linear module constructs a deterministic strong extractor
E: {0,1}^n x {0,1}^D -> {0,1}^(2^h ell), linear in its first argument
at every fixed seed. Its parameters are positive integers n,h,ell,e, with

  a_E=e+h+4, U=4n, rho_E=a_E+2 ceil(log U)+2,
  R=4h, s0=least power of two >=2rho_E,
  S=least power of two >=(R+1)rho_E,
  K_h=ell+2a_E,
  K_(i-1)=ceil[(2R/(R-1))(K_i+a_E+S/2+1)] (i=h,...,1).

The required guard is K_0<=n. Put
m_h=ceil[S ceil(K_(h-1)/(S-rho_E))/2] and let W be the least power of
two >=m_h. The exact complete seed length is D=s0+hS+W. The contract is

  TV((E(X,U_D),U_D),(U_(2^h ell),U_D)) <= 2^(-e),
  D <= 6ell+(h+7)S+18a_E+12,
  K_0 <= 2*2^h(ell+4a_E+S+4), (1)

for every k-source X with k>=K_0 and an independent uniform seed. Its
canonical field/setup/evaluation algorithms have a supplied O(n^8) bit-time
bound on admissible parameters. A target shorter than 2^h ell is obtained
by taking a prefix. A longer prescribed seed is implemented by ignoring its
suffix. Neither operation harms linearity or strongness. Its construction and proof are supplied by the seeded-linear module.
The selected base sampler uses the directly proved rational exact recursion
of linear-and-bilinear-primitives.md Section7f,with its own finite seed and
threshold bounds and actual-length polynomial field setup. Its strong
source-linear contract is the one used by the sampler proof below.

### 1.2 All-seeds-retained affine-breaker contract

Let F:{0,1}^n x {0,1}^d x {0,1}^a -> {0,1}^m be one fixed
deterministic map, with a>=1. Its supplied degree-r contract at source
entropy K and error eta means the following universal statement. For
arbitrary correlated tuples

  Left=(A,A_1,...,A_r),
  Right=(B,B_1,...,B_r,Y,Y_1,...,Y_r),

assume Left and Right are independent, H_inf(A)>=K, Y is uniform on d
bits, and fixed advice strings alpha_i all differ from alpha. Set

  V=F(A+B,Y,alpha), V_i=F(A_i+B_i,Y_i,alpha_i).

Then, with U_m independent of all displayed side information,

  TV((V,V_1,...,V_r,Y,Y_1,...,Y_r),
     (U_m,V_1,...,V_r,Y,Y_1,...,Y_r)) <= eta. (4)

Only A needs entropy and only Y uniformity. No A_i need be a function
of A. All within-side correlations, repetitions, and seed equalities are allowed. Tampered advice
strings need not differ from one another. For a smaller degree, pad with
repeated permissible copies and discard their outputs/seeds. This chosen
branch directly uses the all-seeds-retained contract of the actual affine
construction; it does not need a dummy-tampering adapter for a weaker
contract that retains only tampered seeds.

### 1.3 Probability identities

For any common marginal Z, grouping the TV sum by z gives

  TV((V,Z),(V',Z))=sum_z P(Z=z) TV(V|z,V'|z). (5)

A common deterministic map or random conditional kernel cannot increase TV:
sum differences over each fiber and use the triangle inequality. In all
ideal laws below a displayed U is fresh uniform independent of the entire
retained side information.

Define average conditional min-entropy by

  Havg_inf(A|Z)=-log sum_z P(Z=z) max_x P(A=x|Z=z).

If Z has at most 2^ell possible values, then

  sum_z max_x P(A=x,Z=z) <= 2^ell max_x P(A=x),
  Havg_inf(A|Z) >= H_inf(A)-ell. (6)

If Havg_inf(A|Z)>=K+log(1/delta), with 0<delta<=1, the total probability
of z with H_inf(A|Z=z)<K is at most delta: on that set the maximum
conditional atom exceeds 2^-K, while its expectation is at most
 delta*2^-K. These facts require no unstated independence.

## 2. An explicit small-bias Cayley graph

Fix a field dimension r_F>=1 and put c=2r_F; Section3 chooses r_F=L+8lg(s+1),s=2ceil(B0L/2).

When r_F is a power of two, one admissible choice is the canonical tower field of binary dimension r_F
(Section 3 does not use it). Here is its complete
finite choice. Start with F_2 and trace-one element a_F=1. Repeatedly
adjoin w satisfying w^2+w+a_F=0 and replace a_F by a_F*w. The image of
x -> x^2+x is exactly the trace-zero subspace: its kernel has size two,
each image has trace zero, and the trace is nonzero because the current
a_F has trace one. Thus the quadratic has no root and is irreducible.
The two conjugates of w are w,w+1, so the relative trace of a_F*w is
a_F, and trace transitivity proves absolute trace one at the next level.
Trace transitivity follows by grouping the Frobenius-power sum. This
constructs every power-of-two dimension with a fixed ordered binary basis,
without an irreducible-polynomial search.

In one tower step write x=x0+x1*w,y=y0+y1*w. With
p=x0*y0, q=x1*y1, r=(x0+x1)(y0+y1), t=a_F*q, their product is
(p+t)+(r+p)*w. Four recursive products and linear extra work give O(r_F^2)
bit operations for multiplication in the final field.

Section3 uses r_F=L+8lg(s+1),s=2ceil(B0L/2) and, for every r_F (also powers of two), the canonical field of
weak-seed Section 4 of linear-and-bilinear-primitives.md (the first irreducible binary polynomial of degree
r_F in lexicographic order, found in 2^(r_F)poly(r_F) bit operations) with its
polynomial basis; schoolbook multiplication costs O(r_F^2) bit operations.
Everything below uses only field arithmetic and a fixed binary basis.

For each m>=1 and each label (u,w) in F^2, define an m-bit increment

  g(u,w)_j=<w,u^j>, j=0,...,m-1, (8)

where brackets are the ordinary coordinate dot product and u^0=1 also
when u=0. On vertices F_2^m use the 2^c labelled edges x -> x+g(u,w).
A fixed label is an involution, so the transition matrix P is symmetric
and doubly stochastic, with the uniform stationary distribution.

For any nonzero a in F_2^m, the character (-1)^<a,x> has eigenvalue
E_(u,w)(-1)^<w,p_a(u)>, where p_a(T)=sum_j a_j*T^j. The polynomial is
nonzero of degree at most m-1. Averaging over w annihilates each nonzero
coordinate vector p_a(u), leaving its fraction of zeros. A nonzero degree-j
polynomial has at most j roots, by repeated division by T-alpha. Hence

  0<=eigenvalue<=(m-1)/2^r_F. (9)

The binary characters form an orthogonal basis, so ||P-J||<=lambda=(m-1)/2^r_F,
where J projects onto constants in the uniform-normalized Euclidean inner
product. A labelled step is evaluated by successive field powers and dot
products in O(m*r_F^2) bit operations. The label count 2^c is not the
sampler degree, and the algorithm does not enumerate those labels.

## 3. An exact walk length with logarithmic stretching

Fix B0=2378/7=2(q+8)/7,q=1181,M=18912/7=B0+2q. These are fixed rational
leading coefficients; q is an integer. For integerunitL>=2^40,where
L=S=floor(flg(actualN)/Q_S), put
 s=2ceil(B0L/2),delta_v=s-B0L in[0,2),ellc=lg(s+1),h=8ellc,
 r_F=L+h,c=2r_F,m=s+2q*r_F=ML+16qellc+delta_v,m_s=s. (S6)
Every actual bitlength is computed by these integer formulas; s and m are even.
UseSection2's canonical binary field of degree r_F. An m-bit string is a
uniform s-bit start and q independent complete2r_F-bit labels.
Equation9 gives lambda=(s-1)2^(-r_F),||P-J||<=lambda. (S7)

**Lemma W'.** If |U|>2^(9m/16-8) and |N|<=2^(1-r_F)2^s,
someu inU has a walk position outsideN. This holds for every integerL>=2^40.

Proof. LetR=N in the vertexspace F_2^s,rho=|R|/2^s<=2^(1-r_F).
LetPi be its coordinate projection in the uniform-normalized innerproduct.
BecauseP is symmetric, ||Pi P Pi||<=||Pi J Pi||+||P-J||<=rho+lambda.
The probability a uniform q-stepwalk lies entirely inR is therefore
 <1_R,(Pi P Pi)^q 1_R><=rho(rho+lambda)^q.
Consequently the number of all-bad m-bit strings is at most
 2^(m+1-r_F)((s+1)2^(-r_F))^q. (S9)
This is <=2^(9m/16-8) if
 (7/16)s+9-r_F<=q(r_F/8-log2(s+1)). (S10)
Substitute r_F=L+8ellc and q=7B0/2-8. The difference of the left
and right sides is exactly(7/16)delta_v+9-8ellc+q(log2(s+1)-ellc)
<7/8+9-8ellc<0,because ellc>=2 and ellc>=log2(s+1). Thus all-bad strings are fewer
than the supportU. This is a finite inequality at every allowedL,
not an asymptotic approximation or independent-row assertion.

## 4. The one-step support graph and property (G')

Letp=3M=56736/7,m=s+2q*r_F=(18912/7)L+16qellc+delta_v,
t_tag=3m=pL+48qellc+3delta_v,A=2^t_tag,A0=2^s,
WiththemapE_manditsfourrowsin
sampler-and-disperser.md Section SP.6.2, define

    Gamma(alpha,(i,j)) = first s bits of W(row_i(alpha),j),
                         i=1,...,4, j=0,...,q, (S11)

a graph of degree exactly C=4(q+1)=4728: neighbor index j in [C] is the j-th
label (i,j') in lexicographic order (i first, then j'); no padding is used.

**Theorem G'.** If L>=2^40 satisfies (W'), then for every Q subset [A] with

    |Q| >= Lambda := A^(1/2)/2^15 = 2^(3m/2-15) (S12)

we have |Gamma(Q,[C])| > 2^(1-r_F)*2^s = 2*epsilon0*A0, epsilon0=2^(-r_F).

Proof. Lambda is an integer since s=2ceil(B0L/2) and m=s+2q*r_F are even.
m>=18912/7*2^40>134217728 and |Q|>=2^(3m/2-15),so S1' applies:
|U|>2^(9m/16-8). If N=Gamma(Q,[C]) had |N|<=2^(1-r_F)2^s, Lemma W' would
give u in U and j with first s bits of W(u,j) outside N; but u is a row of
some alpha in Q, so that prefix is Gamma(alpha,(i,j)) in N. Contradiction.

The only guards are L>=2^40 and (W'); both hold for every L>=2^40.

### 4.1 Evaluator and running time

Given L: compute s,m,c,q,m_s by (S6); find p_m by prime-fields Section 2
(2^O(m) bit operations, m=t_tag/3, polynomial in the outer n); build F_(2^r_F) by
weak-seed Section 4 of linear-and-bilinear-primitives.md (2^r_F poly(r_F)=2^L poly(L) bit
operations); for each alpha evaluate its four rows, q steps of (8) per row
(O(m_s r_F^2) bit operations each) and the s-bit prefixes. All A=2^t_tag left
vertices together cost 2^O(t_tag), a fixed polynomial in n since
2^L<=2(n+1). No field enumeration beyond that canonical search, no witness,
and no source-dependent search is used. The outer evaluator tests L>=2^40
among its guards.

## 5. Extractor-to-sampler conversion and graph amplification

### 5.1 Basic conversion with the correct entropy direction

Let E:{0,1}^n x [A0] -> {0,1}^d extract with error epsilon0 from every
source of min-entropy at least K_E, with an independent uniform index.
Strongness is not needed here. For 0<delta<1 it is an
(epsilon0,delta)-sampler at any threshold

  k_S>=K_E+ceil(log(1/delta)): (20)

for every such source X and every test set B of density <=epsilon0,

  P_X[ |{j:E(X,j) in B}|/A0>2epsilon0 ]<=delta. (21)

If the bad-input event had probability p_bad>delta, conditioning X on it
would leave entropy at least H_inf(X)+log p_bad>=K_E, because each
conditional atom is at most its old atom divided by p_bad. Every point
of that conditional source has acceptance probability >2epsilon0. But
extraction bounds the average acceptance by density(B)+epsilon0<=2epsilon0,
a contradiction. Thus sampler entropy increases by the positive overhead
in (20). No inequality with the opposite sign is being invoked. Source
linearity survives because the evaluated map is unchanged.

### 5.2 Amplification, including integer thresholds

Suppose S0 obeys (21) and a fixed map Gamma:[A] x [C] -> [A0] has

  |Gamma(Q,[C])|>2epsilon0*A0 whenever |Q|>=Lambda.

For 0<epsilon1<=epsilon0 and an integer Lambda>=1 require

  Lambda-1<=2epsilon1*A. (22)

Define Samp(x,alpha,j)=S0(x,Gamma(alpha,j)). For every k_S-source X
and every B of density at most epsilon1,

  P_X[ |{alpha: Samp(X,alpha,j) in B for every j}|>
                                      2epsilon1*A ]<=delta. (23)

Proof. Outside the exceptional delta mass of source points, the old bad
indices R_x={u:S0(x,u) in B} have size <=2epsilon0*A0. If Lambda or
more new indices had their entire neighborhoods in R_x, the expansion
hypothesis would give a neighborhood larger than R_x. Hence there are at
most Lambda-1 such new indices, proving (23). Gamma has no randomness in
this argument. Its construction and uniform evaluator, not just its
existence, are required. As it only changes the old seed index, source
linearity is preserved regardless of whether Gamma or its condenser is
linear.

Here S0 is the base sampler of Section 5.3, whose seed is padded to
s=2ceil(B0L/2) bits there. For Theorem G' of Section 4, this proof uses
the test set B only through density(B)<=epsilon0 and the graph only through
expansion of sets of size >=Lambda=A^(1/2)/2^15. It therefore proves, with
Samp(x,alpha,j)=S0(x,Gamma(alpha,j)):

    for every k_S-source X and every B of density <=epsilon0,
    P_X[ |{alpha: Samp(X,alpha,j) in B for every j}| >= Lambda ] <= delta. (23')

Outside the delta-mass, R_x={u:S0(x,u) in B} has |R_x|<=2epsilon0*A0; Lambda
all-bad alpha would form a set Q with Gamma(Q,[C]) inside R_x, contradicting
Theorem G'. No epsilon1 is needed.

### 5.3 Concrete exact base at the sampler unit

Local L denotes S=floor(flg(actualN)/Q_S)>=2^40; actualL=flg(actualN).
Use7f's h=19,b=14,alpha=114792754196981/549755289600,Ro=817286867215042604416161689/110012249584831199340134400,R0=2^40,
s_V=2ceil((2378/7)L/2),r_F=L+8lg(s_V+1),a=r_F+64,
rho=r_F+322+2lgL,U1=2^(128+lgL),error<2^(-r_F).
For eps=2^-20 and2^k<=L<2^(k+1),k>=40,
s_V+1<(2378/7)L+3<511*2^(k+1)+3<=2^(k+10),
so lg(s_V+1)<=k+10 and lgL<=k+1. Thus a-L<=8k+144.
Use the conservative envelope rho-L<=10k+468: indeed
rho=r_F+322+2lgL<=a+322+2lgL, giving 8k+144+322+2(k+1).
The equality from7f is rho=a+2lgU1+2=r_F+322+2lgL;
the 322 already includes a-r_F=64. At L=2^40 the actual overhead
is794 (the larger expression a+322+2lgL gives858); both are<=868.
Sinceceilkfac=1085568053810387218<2^84,
lgK0<=84+lgL andlgU<=actualL+3<Q_S(L+1)+3, so
rho0-(Q_S+1)L<=Q_S+234+9k. The k40 inequalities of(O2a)
and doubling prove a,rho<=(1+eps)L,rho0<=(Q_S+1+eps)L.
Section7f then gives r<=rcL,s<=scL,So<=spoL,K0<=kfacL.
Put i0=((R0+1)/R0)(Q_S+1+eps)+eps. The COMPLETE seed satisfies
 D=s0+hSo+s+r<=(i0+hspo+sc+rc)L=(2378/7)L.
Pad to the actual even s_V=2ceil(B0L/2),B0=2378/7, retaining ALL suffixes.
The displayed rational Ro makes this equality exact. Also
(R0+1)ceilkfac+i0<2^128, so initial and later padding fitU1;
R0<=rho0 follows fromR0=2^40<=L+64. The outer bound
K0<=kfacL<=H*actualL<=actualN supplies the source-length guard.
Since actualL<Q_S(L+1) andL>=2^40, L>=(1-eps)actualL/Q_S.
The exact output calibration is
4^19*alpha*(1-eps)=Q_S*391784144017.
Hence4^19ceil(alphaL)>=delta_O*actualL; truncate to that integer length.
For source failure1/100, kS=K0+7<=kfacL+7
<=(ceil(kfac/Q_S)+7)actualL<=H*actualL; 2^-7<1/100.
The total base error is3*4^19*2^-a<2^(-r_F), matching graph density.
All bounds hold at every guarded integer length. Fields and towers have
uniform polynomial running time by7f; every level uses a fresh seed.

## 6. Conditional affine breaking and the simultaneous bad-seed lemma

### 6.1 Conditional/average-entropy lift with arbitrary copies

Suppose (4) holds at degree r, entropy K, error eta. Let Z be finite,
and replace unconditional independence by the vector Markov condition

  (A,A_1,...,A_r) -- Z -- (B,B_1,...,B_r,Y,Y_1,...,Y_r).

Assume (Y,Z)=(U_d,Z), meaning Y is uniform independently of Z, and
Havg_inf(A|Z)>=K+log(1/delta). With the same fixed advice, the joint
distance in (4), now retaining Z as well, is at most eta+delta.

Proof. By (6), except on z of total mass <=delta, the conditional law of
A has min-entropy >=K. For each such z, conditional independence of the
complete source tuples, uniformity of Y, and every additive identity and
advice inequality give a legal instance of (4). Apply it with error eta.
Use the trivial distance one on exceptional z and average by (5), obtaining
eta*(1-P[bad])+P[bad]<=eta+delta. The arbitrary source-copy quantifiers
are unchanged. Choosing eta=delta=e_F/2 yields joint error e_F at entropy
K+log(2/e_F). This step uses no new evaluated algorithm.

### 6.2 Bad seeds and witness consistency

Fix a source A on n bits; fixed differing advice alpha,alpha_1,...,alpha_r;
and deterministic maps L_0,...,L_r:{0,1}^n -> {0,1}^d. Suppose their
joint leakage Z=(L_0(A),...,L_r(A)) has range size at most 2^ell and

  H_inf(A)>=K+ell+log(1/delta). (27)

For fixed b,w,w_1,...,w_r, define random outputs over A only by

  R_0(b,w)=F(A+b,w+L_0(A),alpha),
  R_i(b,w_i)=F(A+b,w_i+L_i(A),alpha_i).

For 0<theta<=1, call w bad if some b,w_1,...,w_r give

  TV((R_0(b,w),R_1(b,w_1),...,R_r(b,w_r)),
     (U_m,R_1(b,w_1),...,R_r(b,w_r)))>theta. (28)

The density of bad seeds among all d-bit strings is at most

  (eta+delta)/theta. (29)

Proof. Every witness space is finite. For each bad w select one witness
(g(w),f_1(w),...,f_r(w)) by a fixed ordering; define these functions
arbitrarily on other w. Take W uniform independent of A and put

  B=B_i=g(W), A_i=A,
  Y=W+L_0(A), Y_i=f_i(W)+L_i(A),
  X=A+B, X_i=A_i+B_i, Z=(L_0(A),...,L_r(A)).

The use of f_i(W), not W, in Y_i is essential: these are exactly the
chosen tampered unshifted seeds of (28). Since Z is a function of A,
conditional on Z the complete left tuple is independent of W, and the
complete right tuple is a function of W and fixed leakage values. Thus
the vector Markov condition holds. Given Z, Y is a fixed translate of
uniform W, so (Y,Z)=(U_d,Z). Equations (6) and (27) supply the requisite
average conditional entropy.

Apply Section 6.1. The retained pair (Y,Z) determines W=Y+L_0(A), because
L_0(A) is a coordinate of Z. Project the guarantee to the outputs and W;
this is legitimate on both real and ideal laws because W is determined
from retained side information. Its TV is at most eta+delta. Given W=w,
A retains its original unconditional law, because A and W were independent.
Consequently (5) identifies the projected TV with the average, over w,
of exactly the local distances selected in (28). Every bad w contributes
more than theta, proving (29). The witnesses are proof devices, not evaluator inputs or subroutines. Taking
ell=(r+1)d always suffices; an actual smaller leakage range may be used.

## 8. Selected outer use

Proof.md composes this graph/sampler and conditional breaker with tags,
simultaneous bad-seed accounting and its four-moment argument.

Attribution: extractor-to-sampler conversion is due to Zuckerman (1997);
the low-rate-condense/high-rate-walk architecture to Zuckerman, Theory of
Computing 3 (2007), Sections 2--3; the one AB+C support step to
Barak--Impagliazzo--Wigderson (credited in SP.8); the row-XOR/majority architecture to
Chattopadhyay--Liao, arXiv:2110.12652v1, Section 5. The affine interface
is compared to Li, ECCC TR23-023 revision 1. The finite proofs here are
independently written, not a claim of authorship of those constructions.

## 0. Atomic proof map and selected construction

Legacy component names mean the included bodies mapped below. Source URLs
and hashes record provenance, not premises.

- linear-and-bilinear-primitives.md supplies seeded-linear.md Sections1--5,8 and the selected finite scaling,
  its Section7(TheoremD3,CorollaryD3',LemmaP) and7a(exactfields),7f(block-count tree and terminal hash),7g(PXT/PXG for F#/NB/outer),weak-seed-hashing.md1--2,4
  and bilinear-extraction.md. The finite weaker linear primitive is selected,
  without a GUV/CGL premise.
- This file's SP.1--SP.5 prove the ring input; SP.7 proves M/I, SP.6.2 proves S1'.
- restricted-nonmalleability-and-breaking.md supplies prime-fields-and-encoding.md
  1--3, cz-conditional-entropy.md1--8 plus selected local conditionaltools1--3
  and the complete five-round antichain theorem8.2.
  Theselectedring/CZguardisSP.5,gamma=256/20481,eta=268435456/1975781445633,
  p>=2^3200000; CZ3 is elementary,
  so the support-growth proof's use of it is noncircular.
- affine-and-standard-wrappers.md supplies compact-affine-wrappers.md. Its2
  gives the used m<=r affine transformation,3 the advice reduction,3.4 the
  CZ/weak-AGHP reduction, and4 the standard breaker NB. Constant-advice breaking is
  used by NB's five-round antichain14-copy inner F# (restricted8.2,wrappers4.2).
- sampler-and-disperser.md retains the Cayley walk, sampler, seed-retention
  and conditional breaking, with the one-step support graph (Theorem G').
  proof.md now gives its tagged four-moment affine-permutation assembly. The
  original ordinary-majority outer alternative is omitted.
- proof.md retains minimal-circuit-bridge.md and breaker-parameter-assembly-v2.md,
  now with the outer A# and Lemma B9, whose contract (F') supplies the
  directional outer assembly.

Acyclic chain: linear extraction (Theorem D3, Lemma P) -> local core/antichain
-> F#; ring/CZ/weakAGHP -> advice condenser; condenser/AGHP/Ext0/Ext2/F#
-> NB; NB/affine transformation/PXG7g -> A# (F'). Incidence -> S1' ->
graph/walk; Theorem D3 -> sampler. Sampler/tags/four moments -> fixed
affine-permutation family; native B2 ->4n-K*log2(n), zero additive remainder.
The final two components specify its parameters and evaluator.

Aliases: condenser-integration.md Section 6.2, (N1) = wrappers Section 3.4,
(11)--(12), with the joint smoothed law there and in CZ Sections 2,8.
Preparation/audit files and circuit-reduction.md record provenance; the
selected interfaces and circuit bridge are proved in this package.

The direct circuit proof source is GordoAR, submission 6622943f-7b2d-4f8f-b73c-979c5232484a:
https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/226b74d1311ff6711c3c373ab7c8e1fc4d26eb0e/submission
Attribution record:
https://api.yukon.org/api/submissions/6622943f-7b2d-4f8f-b73c-979c5232484a/attribution
Its upstream affine-system framework is Mishra--Santhanam:
https://eccc.weizmann.ac.il/report/2026/260/revision/1/download .
The circuit proof is independently written and proved below, not adopted
through the source's status or an assumed license to republish its prose.

## Paired read-once zero-period lemma

A proper read-once full-B2 formula on 2q distinct leaves, q>=1, with any
pairing of its leaves, has a nonzero zero-period on some product of q
allowed lines in the paired three-corner domain D. A zero-period v means
F(u+v)=F(u) everywhere on that flat. Flip input bits so every pair forbids
00; this preserves properness and the assertion. A safe axis line fixes
one endpoint to 1; the other allowed line is the diagonal 01,10.

First suppose all gates are AND or XNOR (1+x+y). Both have neutral input
1. Induct on q. If a root child contains m>=1 complete pairs internally,
fix its leaves with partners outside that child to 1. Those fixings use
at most one endpoint of each crossing pair. Neutrality reduces the child
to a formula on its 2m internally paired leaves, where m<q. Induction gives
a zero-period supported only on those leaves. The other root child is
unchanged. Extend the chosen internal lines by the safe crossing lines
and arbitrary allowed lines on the other pairs. This proves this case.

Otherwise every pair crosses the root, and each child has q leaves.
If one child contains an XNOR gate, choose a leaf from each of that gate's
two child subtrees, say u,v. Fix every other leaf of the chosen root child
to 1. Leave u,v free by fixing their partners in the other root child to
1. These are consistent one-endpoint-per-pair fixings. Neutrality reduces
the chosen child to XNOR(u,v). Flipping u and v preserves it and leaves
the other child unchanged, giving the desired period. If neither child
contains XNOR, both are pure ANDs. For an AND root choose one pair's
diagonal; the whole output is identically zero. For an XNOR root choose
all pair diagonals: flipping all q free coordinates swaps the two AND
monomials, preserving their XNOR. This includes q=1 and completes induction.

Now take arbitrary proper B2 gates and suppose no such zero-period exists.
For any product gate, consider a child whose leaf set contains at most
one endpoint of every pair. Fix all those leaves to safe 1. Its output
cannot be a controlling value at that gate's pin: otherwise the entire
other child is killed. Choose a leaf z in that other child, fix its mate
to 1, and complete to a product flat leaving z free. This is consistent
with the earlier fixings. Flipping z would preserve the whole formula,
a contradiction. Thus this child's all-ones output is noncontrolling.

Complement internal signals as needed so every signal takes value 1 at
the global all-ones input, absorbing complements in parent truth tables.
Leaves stay unchanged; complementing the root preserves zero-periods.
Every affine proper gate becomes XNOR. Every product gate whose two child
leaf sets each contain at most one endpoint per pair becomes AND: both
pins have noncontrolling value 1, and its output at 11 is 1.
Choose a minimal subtree V containing a complete pair. Every proper
subtree of V has at most one endpoint per pair. Consequently all gates
within V are AND or XNOR. Fix V's leaves with mates outside V to 1;
neutrality leaves a nonempty formula on its internally paired leaves.
Apply the auxiliary induction. Its period changes only those internal
leaves, leaving the rest of the full formula unchanged. Extend the flat
by the safe crossing lines and arbitrary remaining pair lines. This
contradicts the supposition and proves the lemma. The requirement that
ALL leaves be paired is essential to this proof.

## Tight incidence tree and terminal exclusion

Use the representation and incidence graph of proof.md Section 7. Suppose
at its terminal state a=k, q=k, M=2k, s=M-1. All M parameters belong to the
output component. The component bound s_0>=M-1 forces s_0=s. There are M+s
variable vertices and s object vertices, with at most 3s incidences.
Connectedness and M=s+1 force exactly 3s incidences and a bipartite tree;
every object has three distinct variable neighbors.

The output is a gate output. Indeed an affine output on M>=2 parameters
has a nonzero kernel direction v. For any direction v the domain D has
a safe a-flat parallel to v: in each pair choose the allowed line parallel
to its nonzero projection (three corners supply all three directions),
or any allowed line if the projection is zero; keep free coordinates.
An affine output would therefore have zero derivative on such a flat,
contradicting the separation hypothesis of Section 7.

Every auxiliary variable has degree at least two. If it had degree one,
its sole incidence would be a row, since it is not the output. Delete
that row and its uniquely supported auxiliary column. Invertibility of
the remaining block follows by expansion on that column. The represented
function, parameters, Phi and D are unchanged, while s drops to M-2,
contradicting the same component bound. Every nonoutput gate variable has
degree at least two (its defining gate and a consumer); dead gates have
already been removed. The output gate variable has degree one. Parameters
have degree at least one. Their total degree is therefore at least
M+2(s-1)+1=3s, so equality forces every parameter to have degree one and
every nonoutput gate/auxiliary variable to have degree two.

Root this incidence tree at the output variable. Each object has one
rootward and two outward variables. For a row the rootward variable is
auxiliary: parameters are leaves and gate variables never occur in rows.
Solve the row for that auxiliary, obtaining an XOR or XNOR of its two
outward variables. For a gate the rootward variable is its designated
output: the native directed path through consumers to the final output
exists by dead-gate removal and must be the unique incidence path.
Thus no product gate is inverted. Orienting the tree computes a proper
full-B2 read-once formula from all M parameters, each used exactly once.
It agrees with the represented function (even off D). The paired lemma
supplies a nonzero zero-period on a safe k-flat in D. Injective affine Phi
maps it to a nonzero original direction and a k-flat, contradicting
separation. The alleged tight terminal state is impossible.

# SP. Finite rings and direct incidence support

The following proofs supply E_p(256/20481) from cz-conditional-entropy.md
for every prime p >= 2^3200000, without CZ Lemma 7.3, Bourgain's filling
lemma or an unspecified constant. No external review is claimed.

Attribution: the ring/stabilizer framework builds on the
self-improving approximate-ring construction of Terence Tao, *The
sum-product phenomenon in arbitrary rings*, Contributions to Discrete
Mathematics 4(2), 2009, pp. 59--82, especially Lemma 1.4, Proposition 3.1,
and Theorem 1.2. Official source:
https://cdm.ucalgary.ca/article/view/61994
The sumset tools below are Ruzsa's inequalities and Petridis's minimal-
growth argument. The explicit powers, normalization, proofs and application
below are independently supplied here. We work only in commutative rings
with identity. No inclusion AB+AC subset A(B+C) is used: that inclusion is
false in general. Scalar distributions and every product of sets below
are distinguished explicitly.

Notation: nA denotes an n-fold sumset, with 0A={0}; AA={ab:a,b in A};
aA is multiplication by the single element a. A unit is an invertible
ring element. All sets used in cardinality estimates are finite and
nonempty unless indicated otherwise. Products in F_p^d are coordinatewise.

## SP.1. Fully quantified additive tools

### SP.1.1 Ruzsa triangle and covering

For finite nonempty subsets X,Y,Z of an abelian group,

  |X-Z| |Y| <= |X-Y| |Y-Z|. (SP-1)

For each w in X-Z choose one representation w=x_w-z_w. The map
(w,y) -> (x_w-y,y-z_w) is injective: summing the two coordinates recovers
w and the chosen x_w then recovers y. This proves (SP-1).

If |X+Y| <= M|Y|, there is T subset X with |T|<=M and

  X subset Y-Y+T. (SP-2)

Choose a maximal collection of disjoint translates t+Y, t in X. Their
number is at most M. Every other x+Y intersects one of them, which gives
x in t+Y-Y. The same conclusion follows from |X-Y|<=M|Y| by using -Y.

### SP.1.2 Explicit iterated sumset estimate

If |U+V|<=K|U| for finite nonempty U,V and K>=1, then for all nonnegative
integers m,n with m+n>0,

  |mV-nV| <= K^(m+n) |U|. (SP-3)

Here is a proof rather than an invocation of Plunnecke's theorem. Choose
nonempty X subset U minimizing k=|X+V|/|X|. Thus k<=K and |Y+V|>=k|Y|
for every Y subset X, including Y empty. For any finite set C, order its
elements and add the translates X+c successively. For the next c, let
Y={x in X:x+c belongs to the previous union X+C_old}. Then Y+V+c is
contained in the previous union X+V+C_old. Consequently the newly added
part of X+V+c has size at most |X+V|-|Y+V|<=k(|X|-|Y|). The latter
parenthesis is exactly the newly added size of X+c. Summing proves

  |X+V+C| <= k |X+C|.

Induction gives |X+jV|<=k^j|X| for every j>=0. Apply (SP-1) to mV,-X,nV
to obtain |mV-nV|<=k^(m+n)|X|<=K^(m+n)|U|. This includes m=0 or n=0.

In particular, |A+A|<=K|A| implies

  |mA-nA| <= K^(m+n)|A|. (SP-4)

If P contains a set C of size N and |P-P|<=L N, then N<=|P|<=L N,
and (SP-3) with U=P,V=-P yields

  |mP-nP| <= L^(m+n+1)N. (SP-5)

The hypotheses of (SP-5) include C subset P. The bound |P|<=|P-P| follows
by fixing one element of P. Interchanging m and n in (SP-3) gives the stated
signs. The harmless extra power in (SP-5) avoids any unrecorded constants.

### SP.1.3 A sharper form of (SP-5)

If C subset P, |C|=N and |P-P|<=LN, then, for k=m+n>=1,

    |mP-nP| <= L^(m+n) N. (SP-5s)

Indeed put Q=|P|>=N and k0=|P-P|/Q<=LN/Q. Applying (SP-3) with
U=P,V=-P, and interchanging the signs if necessary, gives
|mP-nP|<=k0^k Q<=(LN)^k/Q^(k-1)<=L^k N. All sets are nonempty;
k0>=1 since a translate of P lies in P-P. This improves (SP-5) without
changing its hypotheses. The older bound (SP-5) remains a valid corollary.

## SP.2. A minus stabilizer and direct Ruzsa unit actions

Let C0 be a nonempty set of units,|C0|=M,K>=2,with
 |C0+C0|,|C0C0|<=KM,
 |(C0-C0) outside units|<M/(32768K^17). (SP-Mi0)
Then there areb0 in C0,B0 subsetC0 and a finiteunitalringT such that
 |B0|>=M/(2K), B0 subsetb0T, |T|<=32768K^18M. (SP-Mi1)
Multiplicativeenergy andCauchy--Schwarz give
sum_(a,b in C0)|aC0 intersectbC0|>=M^3/K. Every summand is at mostM.
Pairs with intersection at leastM/(2K) number at leastM^2/(2K).
Chooseb0 with at leastM/(2K) suchpartnersB0,includingb0.
Normalize C=b0^(-1)C0,G=b0^(-1)B0. All hypotheses are preserved;
1 in G subsetC,|G|>=M/(2K),and|gC intersectC|>=M/(2K) foreachg inG.
ForH=gC intersectC,triangle(SP-1)with middle-H gives
 |gC-C|<=|gC+H||C+H|/|H|<=2K^3M.
Covering(SP-2) in BOTH directions therefore gives
 gC subsetC-C+X_g, C subsetgC-gC+Y_g,
 |X_g|,|Y_g|<=2K^3. (SP-Mi2)
Here|gC|=M sinceg is a unit;|C-gC|=|gC-C| supplies the reversecover.

DefineS={x:|xC-C|<=8192K^15M}. Collisionsofxg1-g2 withgi inG
number at least|G|^4/(8192K^15M). For fixedd=g1-g1',at most|G|^2
representations contribute. Thus at leastM/(32768K^17) distinctd
inG-G satisfyxd inG-G. By(SP-Mi0),one is a UNITd.
Write d=g1-g2,xd=g3-g4,and putD=g1C uniong2C uniong3C uniong4C.
By(SP-Mi2),|(-C)+D|<=8K^3M. Apply the FULL two-set estimate(SP-3)
withU=-C,V=D,m=n=2 to obtain
 |2D-2D|<=(8K^3)^4M=4096K^12M.
Matched productsxdC-dC are contained ing3C-g4C-g1C+g2C,which is
contained in2D-2D. Fixed-unit dilation consequently proves
 |xC-C|=|xdC-dC|<=4096K^12M. (SP-Mi3)
No product-of-sets distributivity is used;D need not have small doubling.
RepeatingG-collisioncount at(SP-Mi3) gives at leastM/(16384K^14)
differences. Fewer thanM/(32768K^17)<=M/(32768K^14) are nonunits,
leaving more thanM/(32768K^14) unitchoices. Each pair(d,xd)withd a
unit determinesx. The pairset has size at most|C-C|^2<=K^4M^2,
soS is finite; doublecounting gives|S|<=32768K^18M.
Foranyfinite nonemptyA,the additive extension inequality is
 |A+C-C|<=K^2|A-C|. (SP-Mi4)
To proveit,choose nonemptyP subset-C minimizingk=|P-C|/|P|<=K.
ThePetridisextension argument in theproof of(SP-3),withU=-C,V=-C,
gives|A+P-C|<=k|A+P|. Triangle(SP-1)atX=A-C,Y=-P,Z=-C gives
|A+C-C|<=|A-C+P||P-C|/|P|<=k^2|A+P|<=K^2|A-C|.
HereP subset-C impliesA+P subsetA-C. No cardinality ofxC isassumed.

Translationby1 preservesS: (x+1)C-C subsetxC+C-C,so(SP-Mi4)withA=xC
and(SP-Mi3) bounditby4096K^14M<=8192K^15M. Forunitg inG,
triangle(SP-1)atX=gxC,Y=gC,Z=C gives
 |gxC-C|<=|gxC-gC||gC-C|/M
          <=4096K^12M*2K^3M/M=8192K^15M.
Onlyg,notx,isrequiredtobeaunit. HencegS subsetS.
Multiplicationbyg andtranslationby1 are ambientbijections restricting
to permutations of finiteS; their inverses preserveS. EveryG-monomialm,
includingemptyproduct1,is a unit andits multiplicationM_m preservesS
in both directions. Conjugating translationby1 givestranslationbym.
ThusS is invariant under monomialtranslations,their negatives andfinite
integerlinearcombinations. Since0 inS,the additivegroupT generatedby
G-monomials is contained inS. It contains1 andG andis closed under
products,so is a finiteunital subring. B0=b0G subsetb0T proves(SP-Mi1).
The ring claim concernsT,not arbitraryclosure ofS.

## SP.3. Coordinate-fiber bound on bad differences

Let C subset F_p^d be nonempty, |C+C|<=K|C|, and let pi_j denote the
j-th coordinate projection. Put

  Z_j=(C-C) intersect ker(pi_j).

Choose one representative c from C for each value of pi_j(C). Their
translates c+Z_j are disjoint, so

  |pi_j(C)| |Z_j| <= |C+Z_j| <= |2C-C| <= K^3|C|. (SP-17)

A vector in F_p^d is a nonunit exactly when at least one coordinate is
zero. Taking a union bound in (SP-17) proves

  |(C-C) outside the units|
      <= K^3 |C| sum_(j=1)^d 1/|pi_j(C)|. (SP-18)

This bounds the number of distinct bad differences, not just the number
of pairs yielding them. That distinction is needed in Section SP.2.

## SP.4. Four-coordinate hypotheses

Letp beprimeandA a nonempty subsetofF_p^4. The hypotheses are:

(a) |A|<p^(7/3-1/16);
(b) at least half of A consists of units;
(c) every D subset A with |D|>=p^(-1/16)|A| has both coordinate-pair
    projections pi_12(D), pi_34(D) of size at least p^(5/4);
(d) no such D is contained in either kind of plane

    {(x,y,c1*x,c2*y):x,y in F_p},
    {(x,y,c2*y,c1*x):x,y in F_p},

    for any c1,c2 in F_p.

## SP.5. Four-coordinate expansion from the minus stabilizer

Under exactly Section SP.4's hypotheses(a)--(d),forp>=2^3200000,
 max(|A+A|,|AA|)>p^(256/20481)|A|. (SP-19m)
Putgamma=256/20481 and assume failure. WriteN=|A|,C0=A intersectunits,
M=|C0|,K=2p^gamma. ThenM>=N/2,|C0+C0|,|C0C0|<=KM.
Sincep^(1/16)>2,C0 qualifies for(c),so everyindividualcoordinate
projection has at leastp^(1/4)values. Equation(SP-18) gives
 Z=|(C0-C0) outsideunits|<=4K^3M p^(-1/4).
Now1/4-20gamma=1/81924 and
 131072K^20p^(-1/4)=2^37p^(-1/81924)<1, (SP-21m)
since37*81924=3031188<3200000. ThusZ<M/(32768K^17).
Section SP.2 givesb0,B0,T with
 |B0|>=M/(2K)>=N/(8p^gamma),B0 subsetb0T,|T|<=32768K^18M.
The subsetB0 qualifies for(c)--(d),sincep^(1/16-gamma)>8.
Also32768K^18M<=2^33p^(18gamma)N<p^(1/4)N,
because18*256/20481+33/3200000<1/4. Consequently
 |T|<p^(7/3-1/16+1/4)=p^(121/48)<p^3.
An additive subgroupofF_p^4 is anF_p-vector subspace,so dimT<=2.
Unit scaling preserves pairedprojectioncardinalities. By(c),each
projection ofT has size at leastp^(5/4),hence dimensionat leasttwo.
ThusdimT=2 and both projections bijectontoF_p^2. T is thegraph of
a unitalF_p-linear ringisomorphismphi:F_p^2->F_p^2. Imagesof(1,0),(0,1)
are nonzero complementaryidempotents,so are thatpair in eitherorder.
F_p-linearity makesphi identityorcoordinateswap. Forb0=(b1,b2,b3,b4),
allbi nonzero,b0T is therefore oneof
 {(x,y,(b3/b1)x,(b4/b2)y):x,y inF_p},
 {(x,y,(b3/b2)y,(b4/b1)x):x,y inF_p}.
These are the forbiddenplanes(d),yet qualifyingB0 liesinb0T.
This contradiction proves(SP-19m)withallsubsetquantifiers.

The selectedfamily usesgamma=256/20481,delta=2^-20,
 eta=gamma/(92+delta)=268435456/1975781445633,theta=eta/7=268435456/13830470119431.
Its capeta<=1/7000 holdsbycrossmultiplication. The CZ energy proof only
needs2^64p^(92eta)<=p^gamma,which followsfrometa*logp>=64*2^20
and(92+delta)eta=gamma. This is a direct finite absorption,not an
inference from a weaker93eta theorem. Every actual advice-map call has
eta*r>=7*pi0*w_N>=112*2^20>64*2^20 sincepi0>=16,w_N>=2^20.
Thus u=r+a gives the new guard; p>2^u also meets the independent
ring thresholdp>=2^3200000. All remaining clipping,projection,conditional
19eta/6,W2,N7/U7 and runtime arguments retain their full supplied proofs.

# SP.6. Direct point-line incidence support growth

Section SP.7's full incidence proof gives the AB+C four-row support bound.
The evaluator uses fixed integer arithmetic, canonical primes, AB+C and
bit truncation; no scalar hitting/peeling lemma is used.

## SP.6.0. One support-growth step

### SP.6.1 Exact finite endpoint

Require m>=134217728; this conservative guard implies m>=4. The direct
incidence count below proves the finite endpoint |U|>2^(9m/16-8), keeping
its explicit additive constant rather than replacing it by a fixed exponent loss.

### SP.6.2 The one-step map and its support lemma

For an integerm>=1,letp_m be the least prime strictly above2^m. The supplied
prime-fields Sections1--2 inrestricted-nonmalleability-and-breaking.md prove
2^m<p_m<2^(m+1) and deterministic trialdivision bitcost
O(2^(3m/2+2)(m+1)^2). Read an m-bitstring as its binaryinteger in[0,2^m),
which injects intoF_(p_m). Define
 E_m(x1,x2,x3)=(x1*x2+x3 modp_m) mod2^m. (SP-S2)
For a3m-bitstringalpha=(alpha1,alpha2,alpha3),thefourROWS arealpha1,
alpha2,alpha3,E_m(alpha). This is the same actual deterministic row map.

**Lemma S1' (direct incidence support).** For everyintegerm>=134217728 and
Q subset{0,1}^(3m) with|Q|>=2^(3m/2-15),the unionU of itsfourrowimages
satisfies|U|>2^(9m/16-8).

Proof. SupposeN=|U|<=2^(9m/16-8). N>=1 sinceQ is nonempty. RegardU
as an injective subset ofF_p,p=p_m,and put
 W={y inF_p:the lowm bits of its canonical representative belongtoU}.
Sincep<2^(m+1),each member ofU has atmosttwo preimages,so|W|<=2N.
TakeP=U xW andL={y=a*x+c:a,c inU}. There areatmost2N^2 points and
EXACTLYN^2 distinctnonvertical lines: slope andintercept determinea,c,
includinga=0. Every triple inQ determines its line(a=alpha1,c=alpha3)
and its point(x=alpha2,y=alpha1*alpha2+alpha3 modp). The fourth-row
condition puts y inW. This map into point-line incidences is INJECTIVE:
the line determinesalpha1,alpha3 and the point's first coordinate
 determinesalpha2. Thus|Q|<=I(P,L). Zero slopes require no deletion.

Choose d=ceil(N^(4/3)) as an ANALYSIS parameter,not an evaluator operation.
SinceN<=2^(9m/16),N^(4/3)<=2^(3m/4),andm>=4 gives
 1<=d<=2^(3m/4)+1<2^m<p.
Indeed2^m-2^(3m/4)=2^(3m/4)(2^(m/4)-1)>=8>1 at m>=4.
The complete incidence theorem(SP-I) below therefore applies. Set
 v=floor(2N^2/(d+2))<=2N^(2/3).
As d-1<N^(4/3) andd<=N^(4/3)+1,it gives
 |Q|<=I(P,L)<=2N^2(v+1)+d(d-1)+N^2
 <5N^(8/3)+3N^2+N^(4/3)<=9N^(8/3)<16N^(8/3).
The middle inequality usesN>=1. Taking logarithms with the lower bound
onQ gives log2N>9m/16-57/8>9m/16-8, a contradiction.
All counts concern arbitraryQ and its actual row images; no independence
is asserted. Complete polynomial/multiplicity theorem (SP-I) is proved below.
At the actual call s=2ceil((2378/7)S/2),r_F=S+8lg(s+1),
m=s+2*1181*r_F>=(18912/7)S>=(18912/7)*2^40>134217728.

## SP.7. Complete polynomial incidence proofs

The complete multiplicity and polynomial-incidence argument below is adapted
from published1d75672c-5080-4151-89d2-cd776a09d06b,actualimmutable009e40109ba442d132016de118d39c2fdf92f206,
and Mark Lewko,arXiv2609.27023v3,Lemma2/Section4(CCBY4.0). The full proof
is supplied;neither source acceptance,the paper nor its Lean formalization
is an external premise. The direct support count in SP.6.2 is this candidate's
application. Write t_+=max(t,0),andp=infinity incharacteristiczero.

Algebra used below. Polynomial rings over fields are unique factorization
domains. One elementary justification is induction on the variables:
over a UFD, the product of primitive polynomials is primitive (reduce
modulo each prime divisor of the coefficients; the resulting polynomial
ring over a domain has no zero divisors). Clearing denominators then
transfers factorization and its uniqueness from the fraction field's
Euclidean univariate polynomial ring. This is Gauss's lemma, and also
shows that primitive irreducibles remain irreducible on passing to the
fraction field. Degrees, including degree in a chosen group of variables,
add under multiplication of nonzero polynomials, by their leading forms.
A nonzero univariate polynomial of degree j has at most j distinct roots:
divide successively by the distinct linear factors.

### SP.7.1 Lemma M: multiplicities on parallel lines

Let K be a field, 0<=d<p an integer, and F in K[x,y] nonzero and squarefree
(no repeated irreducible factors), with total degree <=d. Let lambda_i be
finitely many distinct parallel lines, with q_i on lambda_i, such that
F restricts to a nonzero polynomial on every lambda_i. For a nonzero common
direction v let mu_i be the order of u in F(q_i+uv). Then

    sum_i (mu_i-1)_+ <= d(d-1). (SP-M)

Proof. Empty collections and constant F are immediate. Interchange x,y
if the lines are vertical, then shear X=x,Y=y-a*x, where a is their common
slope. The lines become Y=b_i with distinct b_i, and q_i=(xi_i,b_i).
The invertible substitution preserves degree and squarefreeness. Factor
the transformed F as A(Y)G(X,Y), removing all factors independent of X.
Since no restriction is zero, A(b_i)!=0; multiplicities are unchanged.
If G is constant, all mu_i=0. Otherwise each irreducible factor f of G has
1<=deg_X f<=d<p, so f_X is nonzero of smaller X-degree. Write G=fQ.
Squarefreeness gives f not dividing Q; from G_X=f_X Q+f Q_X, primality
shows f cannot divide G_X. Thus G_X!=0 and gcd(G,G_X)=1.

Put h_i=(mu_i-1)_+, d0=deg G and d1=deg G_X<=d0-1. Both
G(xi_i+u,b_i) and G_X(xi_i+u,b_i) are divisible by u^h_i, since
(u^mu g)'=mu*u^(mu-1)g+u^mu g'. Let Pi_T be the space of polynomials of
total degree <=T, of dimension binom(T+2,2). For an integer T sufficiently
large, in particular T>=d0+d1, let W_T be the common kernel of the
sum h_i coefficient functionals [u^j]H(xi_i+u,b_i), 0<=j<h_i.
They are independent: for each (i,j) the polynomial

    (X-xi_i)^j product_(v!=i) (Y-b_v)/(b_i-b_v)

has just that prescribed coefficient equal to one and all others zero.
Take T also at least all their degrees. Therefore codim W_T=sum h_i.
Every multiple of G or G_X lies in this kernel, so
 J_T=G Pi_(T-d0)+G_X Pi_(T-d1) is contained in W_T.
Coprimality and degree additivity give the intersection of these summands
as G G_X Pi_(T-d0-d1). Multiplication is injective. Hence

 codim J_T=binom(T+2,2)-binom(T-d0+2,2)-binom(T-d1+2,2)
             +binom(T-d0-d1+2,2)=d0*d1.

Thus sum h_i=codim W_T<=codim J_T=d0*d1<=d(d-1). This includes d1=0.
The restriction d<p is essential: x^p-y has zero x-derivative and order p
on every horizontal line through (xi,xi^p).

### SP.7.2 Lemma I: the d-form incidence bound

Let P be m points of k^2 and L be n distinct nonvertical lines over a field
k of characteristic p. For every integer 1<=d<p put v=floor(2n/(d+2)). Then

 I(P,L)<=m(v+1)+d(d-1)+n<=m+n+2mn/(d+2)+d(d-1). (SP-I)

Proof. Empty P or L is immediate. Seek S=sum_(j=0)^v z^j S_j(x,y), with
deg S_j<=d, satisfying S(x,ax+b,a)=0 as a POLYNOMIAL IDENTITY in k[x]
for each y=ax+b in L. There are (v+1)binom(d+2,2)>n(d+1) coefficients and
at most n(d+1) homogeneous linear conditions, so a nonzero S exists.
Replace S by the product of its distinct irreducible factors in k[x,y,z].
Neither degree increases. Each old zero restriction forces a factor's
restriction to vanish, since k[x] is a domain, and that factor survives.
Thus S is squarefree, deg_(x,y) S<=d and deg_z S<=v, with all identities.

For q=(q1,q2) in P let r(q) count its incident lines and write
 R_q(u,z)=S(q1+u,q2+zu,z)=sum_j c_(q,j)(z)u^j.
A monomial x^alpha y^beta z^gamma contributes to u^j with z-degree at
most gamma+j, so deg c_(q,j)<=v+j. Every incident slope a is a root of
every c_(q,j), by the imposed identity. Distinct lines through q have
distinct slopes. Root counting gives c_(q,j)=0 when j<r(q)-v. Therefore
u^((r(q)-v)_+) divides R_q in k[z,u].

Work over K=k(z), keeping z indeterminate, and set
lambda_q={(q1+u,q2+zu):u in K}. These parallel lines are distinct:
equality of their intercepts q2-z*q1 forces equality of both coordinates.
S stays squarefree in K[x,y]. Indeed each irreducible factor involving
x or y is primitive over k[z], hence irreducible over k(z) by Gauss's
lemma. Distinct such factors f,g cannot become associates: clearing
denominators would give af=bg for nonzero a,b in k[z]; primality of f and
f not dividing b would force f dividing g. Factors in z alone become units.

Let E={q:R_q=0}, h=|E|. On the remaining lines apply (SP-M), using
mu_q>=(r(q)-v)_+. Their incidences total at most
 (m-h)(v+1)+sum_(q notin E)(r(q)-v-1)_+
 <=(m-h)(v+1)+d(d-1).
For q in E, substitute u=x-q1 in R_q=0. Division in y by the monic
linear polynomial l_q=y-q2-z(x-q1) shows l_q divides S in k[x,y,z].
These factors are irreducible (their quotient ring is k[x,z]), pairwise
non-associate, and of z-degree one. Their product divides S, so h<=v.
A given line containing s points of E contributes s<=1+binom(s,2);
each pair of distinct points lies on at most one line. Consequently
 I(E,L)<=n+binom(h,2)<=n+h(v+1).
Adding gives (SP-I). All zero restrictions above are polynomial identities,
not merely pointwise zeros on a finite field.

## SP.8. Sources and selected scope

The multiplicity and d-form incidence proofs are adapted with attribution
from Mark Lewko,*A Szemeredi-Trotter Theorem in Arbitrary Fields*,
arXiv2609.27023v3(7Oct2026),https://arxiv.org/html/2609.27023v3,
CCBY4.0,and the complete transcription in publishedYukon1d75672c,
commit009e40109ba442d132016de118d39c2fdf92f206. All needed algebra,counting,
polynomial identities and characteristic restrictions are supplied here.
No paper retrieval or formalization/verdict is a mathematical premise.
The direct Q-to-incidence injection,zero-slope handling and finite xi
calibration above are this candidate's changes.

The AB+C support/sampler architecture has its historical hitting/nearhalf
origins in b15a467b and f7447a6f/04e91ff4/7d6289bb through the credited
published precursors. The earlier exact-M witness,scalar10/6,9/4,R7/RG4
peeling and their constants are NOT invoked by the current direct proof.
The failed/unpublishedc1b72 RG4/9/4 route and ownsubmittedfinite-rg-antichain
were material research precursors,credited in the public note;the current
support theorem is the independent direct count just proved. The separate
finite ring/CZ proof remains above,with its own eta_CZ,not a scalarRG exponent.
