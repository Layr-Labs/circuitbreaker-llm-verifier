# Finite support graph, linear sampler, and conditional breaking

For this fixed candidate define the rational constant Q_S=171528522917/1000000000.
Every selected sampler/base unit is the INTEGER S=floor(L/Q_S)
=floor(1000000000L/171528522917),where L=floor(log2(actual input length)).
This constant is fixed before the input length. All field and seed
dimensions remain integers with their specified ceilings; only this
public floor divisor is rational. The floor facts S<=L/Q_S<S+1 and
L<Q_S(S+1) require no integrality of Q_S.


All logarithms are binary, addition is XOR unless a field is specified, and
statistical distance is TV = half the L1 distance. A k-source is a finite
random variable with every atom at most 2^(-k); k may be real. Graphs have
labelled edges and may have loops or repeated edges. Every constant is fixed
independently of the eventual outer input length N.

The support-growth Lemma S1' (sum-product file Section 5), the generalized
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
  K_(i-1)=ceil[(2R/(R-1))(K_i+a_E+S/2+1)]  (i=h,...,1).

The required guard is K_0<=n. Put
m_h=ceil[S ceil(K_(h-1)/(S-rho_E))/2] and let W be the least power of
two >=m_h. The exact complete seed length is D=s0+hS+W. The contract is

  TV((E(X,U_D),U_D),(U_(2^h ell),U_D)) <= 2^(-e),
  D <= 6ell+(h+7)S+18a_E+12,
  K_0 <= 2*2^h(ell+4a_E+S+4),                               (1)

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
     (U_m,V_1,...,V_r,Y,Y_1,...,Y_r)) <= eta.                 (4)

Only A needs entropy and only Y uniformity. No A_i need be a function
of A. All within-side correlations, repetitions, and seed equalities are allowed. Tampered advice
strings need not differ from one another. For a smaller degree, pad with
repeated permissible copies and discard their outputs/seeds. This chosen
branch directly uses the all-seeds-retained contract of the actual affine
construction; it does not need a dummy-tampering adapter for a weaker
contract that retains only tampered seeds.

### 1.3 Probability identities

For any common marginal Z, grouping the TV sum by z gives

  TV((V,Z),(V',Z))=sum_z P(Z=z) TV(V|z,V'|z).                 (5)

A common deterministic map or random conditional kernel cannot increase TV:
sum differences over each fiber and use the triangle inequality. In all
ideal laws below a displayed U is fresh uniform independent of the entire
retained side information.

Define average conditional min-entropy by

  Havg_inf(A|Z)=-log sum_z P(Z=z) max_x P(A=x|Z=z).

If Z has at most 2^ell possible values, then

  sum_z max_x P(A=x,Z=z) <= 2^ell max_x P(A=x),
  Havg_inf(A|Z) >= H_inf(A)-ell.                             (6)

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

  g(u,w)_j=<w,u^j>, j=0,...,m-1,                             (8)

where brackets are the ordinary coordinate dot product and u^0=1 also
when u=0. On vertices F_2^m use the 2^c labelled edges x -> x+g(u,w).
A fixed label is an involution, so the transition matrix P is symmetric
and doubly stochastic, with the uniform stationary distribution.

For any nonzero a in F_2^m, the character (-1)^<a,x> has eigenvalue
E_(u,w)(-1)^<w,p_a(u)>, where p_a(T)=sum_j a_j*T^j. The polynomial is
nonzero of degree at most m-1. Averaging over w annihilates each nonzero
coordinate vector p_a(u), leaving its fraction of zeros. A nonzero degree-j
polynomial has at most j roots, by repeated division by T-alpha. Hence

  0<=eigenvalue<=(m-1)/2^r_F.                                 (9)

The binary characters form an orthogonal basis, so ||P-J||<=lambda=(m-1)/2^r_F,
where J projects onto constants in the uniform-normalized Euclidean inner
product. A labelled step is evaluated by successive field powers and dot
products in O(m*r_F^2) bit operations. The label count 2^c is not the
sampler degree, and the algorithm does not enumerate those labels.

## 3. An exact walk length with logarithmic stretching

Fix B0=2886/7=2(q+8)/7,q=1435,M=22976/7=B0+2q. These are fixed rational
leading coefficients; q is an integer. For integerunitL>=2^40,where
L=S=floor(flg(actualN)/Q_S), put
 s=2ceil(B0L/2),delta_v=s-B0L in[0,2),ellc=lg(s+1),h=8ellc,
 r_F=L+h,c=2r_F,m=s+2q*r_F=ML+16qellc+delta_v,m_s=s.     (S6)
Every actual bitlength is computed by these integer formulas; s and m are even.
UseSection2's canonical binary field of degree r_F. An m-bit string is a
uniform s-bit start and q independent complete2r_F-bit labels.
Equation9 gives lambda=(s-1)2^(-r_F),||P-J||<=lambda.    (S7)

**Lemma W'.** If |U|>2^(9m/16-8) and |N|<=2^(1-r_F)2^s,
someu inU has a walk position outsideN. This holds for every integerL>=2^40.

Proof. LetR=N in the vertexspace F_2^s,rho=|R|/2^s<=2^(1-r_F).
LetPi be its coordinate projection in the uniform-normalized innerproduct.
BecauseP is symmetric, ||Pi P Pi||<=||Pi J Pi||+||P-J||<=rho+lambda.
The probability a uniform q-stepwalk lies entirely inR is therefore
 <1_R,(Pi P Pi)^q 1_R><=rho(rho+lambda)^q.
Consequently the number of all-bad m-bit strings is at most
 2^(m+1-r_F)((s+1)2^(-r_F))^q.                       (S9)
This is <=2^(9m/16-8) if
 (7/16)s+9-r_F<=q(r_F/8-log2(s+1)).                 (S10)
Substitute r_F=L+8ellc and q=7B0/2-8. The difference of the left
and right sides is exactly(7/16)delta_v+9-8ellc+q(log2(s+1)-ellc)
<7/8+9-8ellc<0,because ellc>=2 and ellc>=log2(s+1). Thus all-bad strings are fewer
than the supportU. This is a finite inequality at every allowedL,
not an asymptotic approximation or independent-row assertion.

## 4. The one-step support graph and property (G')

Letp=3M=68928/7,m=s+2q*r_F=(22976/7)L+16qellc+delta_v,
t_tag=3m=pL+48qellc+3delta_v,A=2^t_tag,A0=2^s,
WiththemapE_manditsfourrowsin
sum-product-and-ordinary-condensation.md Section 5.2, define

    Gamma(alpha,(i,j)) = first s bits of W(row_i(alpha),j),
                         i=1,...,4, j=0,...,q,                         (S11)

a graph of degree exactly C=4(q+1)=5744: neighbor index j in [C] is the j-th
label (i,j') in lexicographic order (i first, then j'); no padding is used.

**Theorem G'.** If L>=2^40 satisfies (W'), then for every Q subset [A] with

    |Q| >= Lambda := A^(1/2)/2^15 = 2^(3m/2-15)                         (S12)

we have |Gamma(Q,[C])| > 2^(1-r_F)*2^s = 2*epsilon0*A0, epsilon0=2^(-r_F).

Proof. Lambda is an integer since s=2ceil(B0L/2) and m=s+2q*r_F are even.
m>=22976/7*2^40>134217728 and |Q|>=2^(3m/2-15),so S1' applies:
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

  k_S>=K_E+ceil(log(1/delta)):                                (20)

for every such source X and every test set B of density <=epsilon0,

  P_X[ |{j:E(X,j) in B}|/A0>2epsilon0 ]<=delta.              (21)

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

  Lambda-1<=2epsilon1*A.                                    (22)

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
    P_X[ |{alpha: Samp(X,alpha,j) in B for every j}| >= Lambda ] <= delta.  (23')

Outside the delta-mass, R_x={u:S0(x,u) in B} has |R_x|<=2epsilon0*A0; Lambda
all-bad alpha would form a set Q with Gamma(Q,[C]) inside R_x, contradicting
Theorem G'. No epsilon1 is needed.

### 5.3 Concrete exact 21-rate base at the sampler unit

LocalL denotesS=floor(flg(actualN)/Q_S),notactualL; requireL>=2^40.
Useexactlythe21fixedrates(BV1),thresholds/algorithm(BV2),seedbound(BV3)
andentropyenvelope(BV4)oflinear-and-bilinear-primitives.md7f.
Putalpha=143727/8192,ell=ceil(alphaL),s=2ceil((2886/7)L/2),
r_F=L+8lg(s+1),a=r_F+64,U1=2^(104+lgL),initialrate4194304.
The suppliedproof giveserror2^(-r_F),strongnessintheENTIREseed,
source-linearity,andK0<=kfacL with
 kfac=581302553339164243437027173803912995895567111274880012243132339467033732659373740652365025888423942522306200680204809994521099360826639777547706587916571928949002377176466364561327657269301059849188262854091298185651004796169152476126713493642878916770730240846278880270779773/327674839471526300445919333576497538376944763349112855904437547752471792258638945810178422694610941392702778902462993256439327280659247048616860681277336987030802469914093264447439582943740869788805612832851683987264825348901651917168036621058048000000000000,ceilkfac=1774022547097874534<2^84.
Writeeps=2^-20 andsp_i=(R_i+1)(1+eps)+eps.
For2^k<=L<2^(k+1),k>=40,theexact7f/O2a endpoint-and-doubling proof
givesa+3<=(1+eps)L,rho1<=(1+eps)L,
rho0<=(Q_S+1)L+Q_S+234+9k<=(Q_S+1+eps)L.
EveryS_i<=sp_iL ands0<=[(4194305/4194304)(Q_S+1+eps)+eps]L.
Hence the COMPLETEretainedseed satisfies
 D<=[(4194305/4194304)(Q_S+1+eps)+sum_(i=1)^20 sp_i
      +(1+15/4)sp_21+2eps]L
  =(442688414866272231069137/1073741824000000000000)L<(2886/7)L<=s.                        (BV6)
Padtoactualeven s,retaining/ignoringthefullsuffix. Since s0<200L and
4194305ceilkfac+200<2^104,allinternalpaddingfitsU1. Theguard
4194304<=rho0 isautomatic. K0<=kfacL<=HactualL<=actualN from(O3).
Output4^21ceil(alphaL)>=4^21alphaL. SinceactualL<Q_S(L+1),
L>=2^40 givesL>=(1-eps)actualL/Q_S;the exact inequality
 4^21*alpha*(1-eps)>=Q_S*449853528996
provesoutput>=delta_OactualL,then truncateto thatinteger. Allfield and
hashdegrees remainintegers with the proved actual-length polynomial setup.
Failure1/100 needskS=K0+7<=kfacL+7
 <=(ceil(kfac/Q_S)+7)actualL<=HactualL.
These statements hold at EVERY actual length passingthemainliteralguard.

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

  H_inf(A)>=K+ell+log(1/delta).                              (27)

For fixed b,w,w_1,...,w_r, define random outputs over A only by

  R_0(b,w)=F(A+b,w+L_0(A),alpha),
  R_i(b,w_i)=F(A+b,w_i+L_i(A),alpha_i).

For 0<theta<=1, call w bad if some b,w_1,...,w_r give

  TV((R_0(b,w),R_1(b,w_1),...,R_r(b,w_r)),
     (U_m,R_1(b,w_1),...,R_r(b,w_r)))>theta.                 (28)

The density of bad seeds among all d-bit strings is at most

  (eta+delta)/theta.                                        (29)

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

The source-linear graph/sampler and conditional breaking above are used
by proof.md's tagged affine-permutation outer assembly. That file supplies the
modified advice, simultaneous bad-seed accounting and four-moment proof.
The unused ordinary four-moment outer evaluator and affine consequence
have been omitted; no reference to them is a mathematical dependency.

Attribution: extractor-to-sampler conversion is due to Zuckerman (1997);
the low-rate-condense/high-rate-walk architecture to Zuckerman, Theory of
Computing 3 (2007), Sections 2--3; the one AB+C support step to
Barak--Impagliazzo--Wigderson (credited in the sum-product file); the row-XOR/majority architecture to
Chattopadhyay--Liao, arXiv:2110.12652v1, Section 5. The affine interface
is compared to Li, ECCC TR23-023 revision 1. The finite proofs here are
independently written, not a claim of authorship of those constructions.

## 0. Atomic proof map and selected construction

This map describes proof.md and its companions. Legacy component names refer
to their included bodies, not external files. Source URLs, scratch prefixes
and component hashes give provenance, not mathematical premises.

- linear-and-bilinear-primitives.md supplies seeded-linear.md Sections1--5,8 and the selected finite scaling,
  its Section7(TheoremD3,CorollaryD3',LemmaP) and7a(exactfields),7f(roundedrationalleafbase),7g(exactfinitePXTforCA#/F#/NB),weak-seed-hashing.md1--2,4
  and bilinear-extraction.md. The finite weaker linear primitive is selected,
  without a GUV/CGL premise.
- sum-product-and-ordinary-condensation.md supplies the full finite-ring
  Sections1,3g,4,5,5g;itsSection6 supplies complete polynomial incidence
  lemmasM/I,and5.2 proves the direct AB+C support lemmaS1'.
- restricted-nonmalleability-and-breaking.md supplies prime-fields-and-encoding.md
  1--3, cz-conditional-entropy.md1--8 plus selected local conditionaltools1--3
  and the complete five-round antichain theorem8.2.
  Theselectedring/CZguardisTheorem5g,gamma=256/20481,eta=268435456/1975781445633,
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

The selected acyclic chain is finite linear extraction (Theorem D3, Lemma P)
-> local core and five-round antichain theorem -> inner F#; finite ring/CZ + weakAGHP ->
advice condenser; advice condenser + AGHP + Ext0, Ext2 + F# -> standard
breaker NB; finite affine transformation with NB and PXT7g primitives ->
outer breaker A# (F'); full polynomial incidence proof -> direct support Lemma S1' -> one-step graph + walk, with Theorem D3 ->
sampler; sampler + tags and four moments -> one
fixed affine-permutation family; native B2 bridge -> the eventual inequality
4n-K*log2(n), with zero additive remainder. All parameter and evaluator steps for that exact chain are
specified in the final two components.

Precise legacy equation alias: condenser-integration.md Section 6.2 and its
interface (N1) are affine-and-standard-wrappers.md Section 3.4,
equations (11)--(12), with the joint smoothed conditional-entropy law spelled
out there and in the retained CZ Sections 2,8. They are not generic ordinary
extractor assumptions. References to original preparation/audit files in the
sampler source list are provenance only; all their mathematical interfaces
are reproved in the retained sampler sections and tagged outer assembly. The broader original
circuit-reduction.md is provenance; the complete used bridge follows below.

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
