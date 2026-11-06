# Finite support graph, linear sampler, and conditional breaking

All logarithms are binary, addition is XOR unless a field is specified, and
statistical distance is TV = half the L1 distance. A k-source is a finite
random variable with every atom at most 2^(-k); k may be real. Graphs have
labelled edges and may have loops or repeated edges. Every constant is fixed
independently of the eventual outer input length N.

The support-growth Lemma S1' (final sum-product part Section SP5s), the generalized
seeded-linear extractor (Theorem D3) and the affine breaker yield the finite
graph/sampler and conditional-breaking reductions below. proof.md supplies
the tagged affine-permutation outer assembly and finite guards.

## 1. Supplied contracts and elementary probability facts

### 1.1 Finite seeded-linear contract

The supplied seeded-linear module constructs a deterministic strong extractor
E: {0,1}^n x {0,1}^D -> {0,1}^(2^h ell), linear in its first argument
at each fixed seed. Its parameters are positive integers n,h,ell,e, with

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

for every k-source X with k>=K_0 and an independent uniform seed. Seeded-linear proves this with
O(n^8) bit-time on admissible parameters. Prefix truncation below 2^h ell
and ignored seed suffixes preserve linearity/strongness.
The selected base uses linear-and-bilinear-primitives.md Section7f's
rational exact recursion, seed/threshold bounds and polynomial setup. Its strong
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

Only A needs entropy and Y uniformity. The A_i need not be functions of A;
within-side correlations, repetitions and seed equalities are unrestricted.
Tampered advice may coincide. Smaller degrees pad by permissible copies,
then discard outputs/seeds. The actual affine contract retains all seeds.

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
of z with H_inf(A|Z=z)<K is <= delta: on that set the maximum
conditional atom exceeds 2^-K, while its expectation is <=
 delta*2^-K. These facts require no unstated independence.

## 2. An explicit small-bias Cayley graph

Fix a field dimension r_F>=1 and put c=2r_F; Section 3 chooses r_F=L.

For power-of-two r_F, the canonical tower is an admissible field
(unused in Section 3). Start at F_2,a_F=1; adjoin w^2+w+a_F=0,
then replace a_F by a_F*w. The map x->x^2+x has two-element kernel,
hence half-field image; each image has trace zero, and the trace is nonzero
because a_F has trace one. The image is exactly the trace-zero subspace, so the quadratic is irreducible. Conjugates w,w+1 give
relative trace a_F to a_F*w; grouping Frobenius sums proves trace
transitivity and absolute trace one. This inductively gives every
power-of-two dimension with fixed ordered basis, without polynomial search.

In one tower step write x=x0+x1*w,y=y0+y1*w. With
p=x0*y0, q=x1*y1, r=(x0+x1)(y0+y1), t=a_F*q, their product is
(p+t)+(r+p)*w. Four recursive products and linear extra work give O(r_F^2)
bit operations for multiplication in the final field.

Section 3 uses r_F=L and, for every L (also when L is a power of two), the canonical field of
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

## 3. A walk with a short start

Fix B0=392 (the base seed coefficient of Section 5.3), P3=215, and an
integer sampler unit L>=1 (L=S=floor(floor(log2 actualN)/144) in proof.md). Write lg(x)=ceil(log2 x) for positive x. Put

    s=B0*L,   m=P3*s,   r_F=L,   c=2r_F=2L,
    q=floor((m-s)/c)=floor((P3-1)B0/2)=41944,   m_s=m-q*c.          (S6)

Since(P3-1)B0=83888 is even,2q=(P3-1)B0 and m_s=s. In particular m_s<=s+L. Use Section 2 with
field dimension r_F=L (the canonical field of weak-seed Section 4 of
linear-and-bilinear-primitives.md) and vertex space F_2^(m_s): an m-bit
string is an m_s-bit start x_0 followed by q labels (u,w) in F_(2^L)^2 of c
bits each, and x_j=x_(j-1)+g(u_j,w_j) with g from (8) at length m_s. Output
W(x,j)=x_j for j=0,...,q. Equation (9) at m=m_s gives

    0<=eigenvalue<=(m_s-1)/2^r_F,   so ||P-J||<=lambda:=(m_s-1)/2^r_F.    (S7)

Walk margin. With xi=1/427 from sampler-and-disperser.md, final sum-product part
Section SP5s,

    mu=(P3-1)xi-1/2=214/427-1/2=1/854>0.                       (S8)

**Lemma W'.** AssumeL>=1024andthewalkguard

    (P3-1)*lg((B0+2)L) <= mu*L.                                     (W')

(It holds for every L>=2^32: on [2^k,2^(k+1)),
lg((B0+2)L)<=10+k for the selected B0=392. At k32,
214*(10+32)=8988<2^32/854; left increases by214 while right doubles.) If U subset {0,1}^m has |U|>2^((1/2+xi)m) and N subset
{0,1}^s has |N|<=2^(1-L)*2^s, then some u in U has a walk position W(u,j)
whose first s bits lie outside N.

Proof. Let R subset F_2^(m_s) be the strings whose s-bit prefix lies in N;
each prefix has 2^(m_s-s) extensions, so its density is
rho=|N|/2^s<=2^(1-L). Call x in {0,1}^m all-bad if every W(x,j) lies in R. A
uniform m-bit string is exactly a uniform start with q independent uniform
labels, i.e. a q-step walk with transition matrix P from the uniform
distribution. Let Pi be the coordinate projection onto R, in the
uniform-normalized inner product. Since ||Pi J Pi||=rho and
||Pi(P-J)Pi||<=lambda by (S7), ||Pi P Pi||<=rho+lambda, and

    P[every x_j in R] = <1_R,(Pi P Pi)^q 1_R> <= ||1_R||^2 (rho+lambda)^q
                      = rho(rho+lambda)^q,

q=0 included. Hence

    #all-bad <= 2^m rho (rho+lambda)^q <= 2^(m+1-L) ((m_s+1)2^(-r_F))^q,  (S9)

using rho<=2^(1-L)=2*2^(-r_F) and (S7). It suffices that
(S9)<=2^((1/2+xi)m)<|U|: then not every point of U is all-bad. Taking logs and
substituting m=m_s+2q*r_F, (S9)<=2^((1/2+xi)m) is

    (1/2-xi)m_s+1-L <= q(2xi*r_F-ell),   ell=log(m_s+1).               (S10)

Leftside:(1/2-xi)m_s+1-L< s/2+1,becausem_s<=s+L. Right side:
m_s+1<=(B0+2)L, so (W') gives ell<=lg((B0+2)L)<=mu*L/(P3-1)<xi*L, and
D=2xi*L-ell>0. With q>=((P3-1)B0-1)/2, (P3-1)xi=1/2+mu and (P3-1)ell<=mu*L,

    q*D >= ((P3-1)B0-1)xi*L-(P3-1)B0*ell/2
        >= (1/2+mu)B0*L-xi*L-mu*B0*L/2 = s/2+(mu*B0/2-xi)L >= s/2+1,

because (mu*B0/2-xi)L=(100/427)L>=1 for L>=1024. So (S10) holds.

## 4. The one-step support graph and property (G')

Let p=3*P3*B0=P*B0 with P=645, so p=252840; put t=p*L, A=2^t, A0=2^s and
m=t/3=P3*s as in (S6). With the map E_m and the four rows of
sampler-and-disperser.md, final sum-product part Section SP5s.2, define

    Gamma(alpha,(i,j)) = first s bits of W(row_i(alpha),j),
                         i=1,...,4, j=0,...,q,                         (S11)

a graph of degree exactly C=4(q+1)=167780: neighbor index j in [C] is the j-th
label (i,j') in lexicographic order (i first, then j'); no padding is used.

**Theorem G'.** If L>=4096 satisfies (W'), then for every Q subset [A] with

    |Q| >= Lambda := A^(1/2)/2^15 = 2^(pL/2-15)                         (S12)

we have |Gamma(Q,[C])| > 2^(1-L)*2^s = 2*epsilon0*A0, epsilon0=2^-L.

Proof. Lambda is an integer since p is even and pL/2>=15. Also
m=pL/3=P3*B0*L>=215*392*4096=345210880>300000000 and |Q|>=2^(3m/2-15), so
Lemma S1' there applies: the union U of the row images has
|U|>2^((1/2+xi)m). If N=Gamma(Q,[C]) had |N|<=2^(1-L)2^s, Lemma W' would
give u in U and j with first s bits of W(u,j) outside N; but u is a row of
some alpha in Q, so that prefix is Gamma(alpha,(i,j)) in N. Contradiction.

The only guards are L>=4096 and (W'); both hold for every L>=2^32.

### 4.1 Evaluator and running time

Given L: compute s,m,c,q,m_s by (S6); find p_m by prime-fields Section 2
(2^O(m) bit operations, m=pL/3, polynomial in the outer n); build F_(2^L) by
weak-seed Section 4 of linear-and-bilinear-primitives.md (2^L poly(L) bit
operations); for each alpha evaluate its four rows, q steps of (8) per row
(O(m_s L^2) bit operations each) and the s-bit prefixes. All A=2^(pL) left
vertices together cost 2^O(pL), a fixed polynomial in n since
2^L<=2^actualL<=n (local L=S). No field enumeration beyond that canonical search, no witness,
and no source-dependent search is used. The outer evaluator tests L>=2^32
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
conditional atom is <= its old atom divided by p_bad. Every point
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
most Lambda-1 such new indices, proving (23). Gamma has no randomness here. Its construction and uniform evaluator, not just its
existence, are required. As it only changes the old seed index, source
linearity is preserved irrespective of whether Gamma or its condenser is
linear.

Here S0 is the base sampler of Section 5.3, whose seed is padded to
s=B0*L bits there. For Theorem G' of Section 4, this proof uses
the test set B only through density(B)<=epsilon0 and the graph only through
expansion of sets of size >=Lambda=A^(1/2)/2^15. It thus proves, with
Samp(x,alpha,j)=S0(x,Gamma(alpha,j)):

    for every k_S-source X and every B of density <=epsilon0,
    P_X[ |{alpha: Samp(X,alpha,j) in B for every j}| >= Lambda ] <= delta.  (23')

Outside the delta-mass, R_x={u:S0(x,u) in B} has |R_x|<=2epsilon0*A0; Lambda
all-bad alpha would form a set Q with Gamma(Q,[C]) inside R_x, contradicting
Theorem G'. No epsilon1 is needed.

### 5.3 Concrete exact base at the sampler unit

Local L denotes S=floor(floor(log2 actualN)/144), not actualL. Require L>=2^40.
Use7f's 26 four-way levels and its fixed alpha,Ro,Rl,R0, with
ell=ceil(alphaL), a=L+64, error<2^-L, U1=2^(128+lgL).
Put eps=2^-20, spo=(Ro+1)(1+eps)+eps, spl=(Rl+1)(1+eps)+eps.
On 2^k<=L<2^(k+1), k>=40, 600+2k<=2^(k-20): this holds at40
and persists by doubling. It bounds both rho1's overhead
322+2(k+1) and rho0's overhead 287+(k+1).
Hence a,rho1<=(1+eps)L, ell<=(alpha+eps)L, So<=spoL, Sl<=splL.
Section7f proves K0<=kfacL and ceilkfac=12484765614829834828034<2^75.
Since actualL<=144L+143 and U=4actualN, lgU<=actualL+3<=144(L+1)+2 and lgK0<= 75+lgL,
 rho0<=(144+1)L+287+lgL<=(145+eps)L,
 s0<=[((R0+1)/R0)(145+eps)+eps]L=:i0*L.
Section7f bounds the COMPLETE retained AGHP seed by (2alpha+4+8eps)L.
Adding the initial and all condenser seeds, with one extra eps allowance,
 D<=[((R0+1)/R0)(145+eps)+25spo+spl+2alpha+4+10eps]L
   =392L.
Pad to B0L with B0=392, retaining the suffix. The explicit parameter
formula for Ro in7f makes this seed equality exact. Rational substitution
gives (R0+1)ceilkfac+i0<2^128, so initial and internal blocks fit U1.
R0<=rho0 follows from rho0>=L+64 and R0=2^40. The outer H bounds
K0<=kfacL<=H*actualL<=actualN, supplying the initial entropy guard.
Output is 2^52ceil(alphaL)>=2^52alphaL. Since actualL<144(L+1),
L>=(1-eps)actualL/144. The exact output comparison is
 2^52alpha(1-eps)-144*128096781715697=0.
Thus output>=delta_O*actualL and may be truncated to that integer.
Every base/hash degree is O(log actualN); larger extensions use towers.
For failure 1/100, kS=K0+7<=kfacL+7<=(ceil(kfac/144)+7)actualL<=H*actualL.
All these bounds hold at every guarded integer length, by the dyadic
endpoint and doubling argument; no recurrence sampling is a premise.

## 6. Conditional affine breaking and the simultaneous bad-seed lemma

### 6.1 Conditional/average-entropy lift with arbitrary copies

Suppose (4) holds at degree r, entropy K, error eta. Let Z be finite,
and replace unconditional independence by the vector Markov condition

  (A,A_1,...,A_r) -- Z -- (B,B_1,...,B_r,Y,Y_1,...,Y_r).

Assume (Y,Z)=(U_d,Z), meaning Y is uniform independently of Z, and
Havg_inf(A|Z)>=K+log(1/delta). With the same fixed advice, the joint
distance in (4), now retaining Z as well, is <= eta+delta.

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

The density of bad seeds among all d-bit strings is <=

  (eta+delta)/theta.                                        (29)

Proof. Every witness space is finite. For each bad w select one witness
(g(w),f_1(w),...,f_r(w)) by a fixed ordering; define these functions
arbitrarily on other w. Take W uniform independent of A and put

  B=B_i=g(W), A_i=A,
  Y=W+L_0(A), Y_i=f_i(W)+L_i(A),
  X=A+B, X_i=A_i+B_i, Z=(L_0(A),...,L_r(A)).

The use of f_i(W), not W, in Y_i is essential: these are exactly the
chosen tampered unshifted seeds of (28). Since Z is determined by A,
conditional on Z the complete left tuple is independent of W, and the
complete right tuple is determined by W and fixed leakage values. Thus
the vector Markov condition holds. Given Z, Y is a fixed translate of
uniform W, so (Y,Z)=(U_d,Z). Equations (6) and (27) supply the requisite
average conditional entropy.

Apply Section 6.1. The retained pair (Y,Z) determines W=Y+L_0(A), because
L_0(A) is a coordinate of Z. Project the guarantee to the outputs and W;
this is legitimate on both real and ideal laws because W is determined
from retained side information. Its TV is <= eta+delta. Given W=w,
A retains its original unconditional law, because A and W were independent.
Thus (5) identifies the projected TV with the average, over w,
of exactly the local distances selected in (28). Every bad w contributes
more than theta, proving (29). The witnesses are proof devices, not evaluator inputs or subroutines. Taking
ell=(r+1)d always suffices; an actual smaller leakage range may be used.

## 8. Selected outer use

proof.md applies this source-linear sampler and conditional breaking with
tagged advice, simultaneous bad-seed accounting and four moments.
The omitted ordinary-majority alternative is unused.

Attribution: extractor-to-sampler conversion is due to Zuckerman (1997);
the low-rate-condense/high-rate-walk architecture to Zuckerman, Theory of
Computing 3 (2007), Sections 2--3; the one AB+C support step to
Barak--Impagliazzo--Wigderson (credited in the final sum-product part); the row-XOR/majority architecture to
Chattopadhyay--Liao, arXiv:2110.12652v1, Section 5. The affine interface
is compared to Li, ECCC TR23-023 revision 1. The finite proofs here are
independent, not a claim of authorship of those constructions.

## 0. Atomic proof map and selected construction

Legacy filenames label included components; they are not external premises.
- linear-and-bilinear-primitives.md: seeded-linear.md Sections 1--6,8; Section7 (TheoremD3,
  CorollaryD3',LemmaP),7a (exact fields),7f (rounded rational base),7g (PXT/PXG);
  weak-seed-hashing.md1--2,4; bilinear-extraction.md. The finite weaker linear primitive is selected, without a GUV/CGL premise.
- sampler-and-disperser.md, final sum-product part: finite-approximate-ring.md
  Sections SP1,SP3g,SP4,SP5,SP5g; ordinary-condenser.md scalar (30'), Section SP6
  (M and I); RG and S1' (one AB+C step).
- restricted-nonmalleability-and-breaking.md: prime-fields-and-encoding.md1--3;
  cz-conditional-entropy.md1--8; constant-advice-breaker.md1--8.
  Ring/CZ uses Theorem SP5g, gamma=256/20481, eta=256000000/2109543000001,
  p>=2^3200000. CZ3 is elementary, so S1' is noncircular.
- affine-and-standard-wrappers.md: compact-affine-wrappers.md; Section2
  affine transformation m<=r,3 advice,3.4 CZ/weak-AGHP,4 NB.
  Single-copy CA# is unselected; F# uses restricted Section8.2.
- This file: Cayley walk, sampler, seed-retention, conditional breaking, graph G', paired/tree lemmas.
- proof.md: minimal-circuit-bridge.md and breaker-parameter-assembly-v2.md;
  A#, B9 and tagged four-moment affine-permutation assembly.

Acyclic chain: finite linear extraction (Theorem D3, PXT/PXG) -> multicopy F#;
finite ring/CZ + IP -> advice; advice + AGHP + Ext0,Ext2 + F# -> NB;
affine transformation with NB/PXT/PXG -> A# (F');
scalar (30) -> RG,S1' -> graph/walk with D3 -> sampler;
sampler/tags/four moments -> fixed family -> native B2 bound
4n-K*log2(n), K=4H, K0=0. All parameter and evaluator steps are specified
in the final two components, sampler-and-disperser.md and proof.md.
Legacy condenser-integration.md6.2/(N1) means wrappers3.4/(11)--(12),
with the joint smoothed entropy law there and in CZ2,8.
Preparation/audit filenames and circuit-reduction.md are provenance only;
all used interfaces and the bridge are reproved here.

Circuit source: GordoAR, 6622943f-7b2d-4f8f-b73c-979c5232484a:
https://github.com/Layr-Labs/circuitbreaker-llm-verifier/tree/226b74d1311ff6711c3c373ab7c8e1fc4d26eb0e/submission
https://api.yukon.org/api/submissions/6622943f-7b2d-4f8f-b73c-979c5232484a/attribution
Affine-system antecedent: Mishra--Santhanam,
https://eccc.weizmann.ac.il/report/2026/260/revision/1/download .
Source status and license are not proof premises.

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
subtree of V has at most one endpoint per pair. Thus all gates
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
An affine output would thus have zero derivative on such a flat,
contradicting the separation hypothesis of Section 7.

Every auxiliary variable has degree at least two. If it had degree one,
its sole incidence would be a row, since it is not the output. Delete
that row and its uniquely supported auxiliary column. Invertibility of
the remaining block follows by expansion on that column. The represented
function, parameters, Phi and D are unchanged, while s drops to M-2,
contradicting the same component bound. Every nonoutput gate variable has
degree at least two (its defining gate and a consumer); dead gates have
already been removed. The output gate variable has degree one. Parameters
have degree at least one. Their total degree is thus at least
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


# Final part: Sum-product and ordinary condensation

# Explicit finite approximate rings and four-coordinate expansion

Finite proofs below give E_p(256/20481) for every prime p>=2^3200000,
without CZ Lemma 7.3, Bourgain filling or unspecified constants.
Condenser proofs are separate; no external review is claimed.

Attribution: the ring/stabilizer framework builds on the
self-improving approximate-ring construction of Terence Tao, *The
sum-product phenomenon in arbitrary rings*, Contributions to Discrete
Mathematics 4(2), 2009, pp. 59--82, especially Lemma 1.4, Proposition 3.1,
and Theorem 1.2. Official source:
https://cdm.ucalgary.ca/article/view/61994
Sumset tools: Ruzsa inequalities and Petridis minimal growth. Explicit
powers, normalization and proofs are supplied here. Rings are commutative
with identity. AB+AC subset A(B+C) is false and unused; scalar distributions
and set products are distinguished.

Notation: nA denotes an n-fold sumset, with 0A={0}; AA={ab:a,b in A};
aA is multiplication by the single element a. A unit is an invertible
ring element. All sets used in cardinality estimates are finite and
nonempty unless indicated otherwise. Products in F_p^d are coordinatewise.

## SP1. Fully quantified additive tools

### SP1.1 Ruzsa triangle and covering

For finite nonempty subsets X,Y,Z of an abelian group,

  |X-Z| |Y| <= |X-Y| |Y-Z|.                              (1)

For each w in X-Z choose one representation w=x_w-z_w. The map
(w,y) -> (x_w-y,y-z_w) is injective: summing the two coordinates recovers
w and the chosen x_w then recovers y. This proves (1).

If |X+Y| <= M|Y|, there is T subset X with |T|<=M and

  X subset Y-Y+T.                                        (2)

Choose a maximal collection of disjoint translates t+Y, t in X. Their
number is <= M. Every other x+Y intersects one of them, which gives
x in t+Y-Y. The same conclusion follows from |X-Y|<=M|Y| by using -Y.

### SP1.2 Explicit iterated sumset estimate

If |U+V|<=K|U| for finite nonempty U,V and K>=1, then for all nonnegative
integers m,n with m+n>0,

  |mV-nV| <= K^(m+n) |U|.                                (3)

Here is a proof rather than an invocation of Plunnecke's theorem. Choose
nonempty X subset U minimizing k=|X+V|/|X|. Thus k<=K and |Y+V|>=k|Y|
for every Y subset X, including Y empty. For any finite set C, order its
elements and add the translates X+c successively. For the next c, let
Y={x in X:x+c belongs to the previous union X+C_old}. Then Y+V+c is
contained in the previous union X+V+C_old. Thus the newly added
part of X+V+c has size at most |X+V|-|Y+V|<=k(|X|-|Y|). The latter
parenthesis is exactly the newly added size of X+c. Summing proves

  |X+V+C| <= k |X+C|.

Induction gives |X+jV|<=k^j|X| for every j>=0. Apply (1) to mV,-X,nV
to obtain |mV-nV|<=k^(m+n)|X|<=K^(m+n)|U|. This includes m=0 or n=0.

Hence, |A+A|<=K|A| implies

  |mA-nA| <= K^(m+n)|A|.                                 (4)

If P contains a set C of size N and |P-P|<=L N, then N<=|P|<=L N,
and (3) with U=P,V=-P yields

  |mP-nP| <= L^(m+n+1)N.                                 (5)

The hypotheses of (5) include C subset P. The bound |P|<=|P-P| follows
by fixing one element of P. Interchanging m and n in (3) gives the stated
signs. The harmless extra power in (5) avoids any unrecorded constants.

### SP1.3 A sharper form of (5)

If C subset P, |C|=N and |P-P|<=LN, then, for k=m+n>=1,

    |mP-nP| <= L^(m+n) N.                                (5s)

Indeed put Q=|P|>=N and k0=|P-P|/Q<=LN/Q. Applying (3) with
U=P,V=-P, and interchanging the signs if necessary, gives
|mP-nP|<=k0^k Q<=(LN)^k/Q^(k-1)<=L^k N. All sets are nonempty;
k0>=1 since a translate of P lies in P-P. This improves (5) without
changing its hypotheses. The older bound (5) remains a valid corollary.

## SP3g. A minus stabilizer and direct Ruzsa unit actions

Let C0 be a nonempty set of units,|C0|=M,K>=2,with
 |C0+C0|,|C0C0|<=KM,
 |(C0-C0) outside units|<M/(32768K^17).                 (Mi0)
Then there areb0 in C0,B0 subsetC0 and a finiteunitalringT such that
 |B0|>=M/(2K), B0 subsetb0T, |T|<=32768K^18M.           (Mi1)
Multiplicativeenergy andCauchy--Schwarz give
sum_(a,b in C0)|aC0 intersectbC0|>=M^3/K. Every summand is at mostM.
Pairs with intersection at leastM/(2K) number at leastM^2/(2K).
Chooseb0 with at leastM/(2K) suchpartnersB0,includingb0.
Normalize C=b0^(-1)C0,G=b0^(-1)B0. All hypotheses are preserved;
1 in G subsetC,|G|>=M/(2K),and|gC intersectC|>=M/(2K) foreachg inG.
ForH=gC intersectC,triangle(1)with middle-H gives
 |gC-C|<=|gC+H||C+H|/|H|<=2K^3M.
Covering(2) in BOTH directions thus gives
 gC subsetC-C+X_g, C subsetgC-gC+Y_g,
 |X_g|,|Y_g|<=2K^3.                                   (Mi2)
Here|gC|=M sinceg is a unit;|C-gC|=|gC-C| supplies the reversecover.

DefineS={x:|xC-C|<=8192K^15M}. Collisionsofxg1-g2 withgi inG
number at least|G|^4/(8192K^15M). For fixedd=g1-g1',at most|G|^2
representations contribute. Thus at leastM/(32768K^17) distinctd
inG-G satisfyxd inG-G. By(Mi0),one is a UNITd.
Write d=g1-g2,xd=g3-g4,and putD=g1C uniong2C uniong3C uniong4C.
By(Mi2),|(-C)+D|<=8K^3M. Apply the FULL two-set estimate(3)
withU=-C,V=D,m=n=2 to obtain
 |2D-2D|<=(8K^3)^4M=4096K^12M.
Matched productsxdC-dC are contained ing3C-g4C-g1C+g2C,which is
contained in2D-2D. Fixed-unit dilation thus proves
 |xC-C|=|xdC-dC|<=4096K^12M.                           (Mi3)
No product-of-sets distributivity is used;D need not have small doubling.
RepeatingG-collisioncount at(Mi3) gives at leastM/(16384K^14)
differences. Fewer thanM/(32768K^17)<=M/(32768K^14) are nonunits,
leaving more thanM/(32768K^14) unitchoices. Each pair(d,xd)withd a
unit determinesx. The pairset has size at most|C-C|^2<=K^4M^2,
soS is finite; doublecounting gives|S|<=32768K^18M.
Foranyfinite nonemptyA,the additive extension inequality is
 |A+C-C|<=K^2|A-C|.                                   (Mi4)
To proveit,choose nonemptyP subset-C minimizingk=|P-C|/|P|<=K.
ThePetridisextension argument in theproof of(3),withU=-C,V=-C,
gives|A+P-C|<=k|A+P|. Triangle(1)atX=A-C,Y=-P,Z=-C gives
|A+C-C|<=|A-C+P||P-C|/|P|<=k^2|A+P|<=K^2|A-C|.
HereP subset-C impliesA+P subsetA-C. No cardinality ofxC isassumed.

Translationby1 preservesS: (x+1)C-C subsetxC+C-C,so(Mi4)withA=xC
and(Mi3) bounditby4096K^14M<=8192K^15M. Forunitg inG,
triangle(1)atX=gxC,Y=gC,Z=C gives
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
products,so is a finiteunital subring. B0=b0G subsetb0T proves(Mi1).
The ring claim concernsT,not arbitraryclosure ofS.

## SP4. Coordinate-fiber bound on bad differences

Let C subset F_p^d be nonempty, |C+C|<=K|C|, and let pi_j denote the
j-th coordinate projection. Put

  Z_j=(C-C) intersect ker(pi_j).

Choose one representative c from C for each value of pi_j(C). Their
translates c+Z_j are disjoint, so

  |pi_j(C)| |Z_j| <= |C+Z_j| <= |2C-C| <= K^3|C|.         (17)

A vector in F_p^d is a nonunit exactly when at least one coordinate is
zero. Taking a union bound in (17) proves

  |(C-C) outside the units|
      <= K^3 |C| sum_(j=1)^d 1/|pi_j(C)|.                (18)

This bounds the count of distinct bad differences, not just the count of pairs yielding them. That distinction is needed in Section SP3g.

## SP5. Four-coordinate hypotheses

Letp beprimeandA a nonempty subsetofF_p^4. The hypotheses are:

(a) |A|<p^(7/3-1/16);
(b) at least half of A consists of units;
(c) every D subset A with |D|>=p^(-1/16)|A| has both coordinate-pair
    projections pi_12(D), pi_34(D) of size at least p^(5/4);
(d) no such D is contained in either kind of plane

    {(x,y,c1*x,c2*y):x,y in F_p},
    {(x,y,c2*y,c1*x):x,y in F_p},

    for any c1,c2 in F_p.

## SP5g. Four-coordinate expansion from the minus stabilizer

Under exactly Section SP5's hypotheses(a)--(d),forp>=2^3200000,
 max(|A+A|,|AA|)>p^(256/20481)|A|.                     (19m)
Putgamma=256/20481 and assume failure. WriteN=|A|,C0=A intersectunits,
M=|C0|,K=2p^gamma. ThenM>=N/2,|C0+C0|,|C0C0|<=KM.
Sincep^(1/16)>2,C0 qualifies for(c),so everyindividualcoordinate
projection has at leastp^(1/4)values. Equation(18) gives
 Z=|(C0-C0) outsideunits|<=4K^3M p^(-1/4).
Now1/4-20gamma=1/81924 and
 131072K^20p^(-1/4)=2^37p^(-1/81924)<1,               (21m)
since37*81924=3031188<3200000. ThusZ<M/(32768K^17).
Section SP3g givesb0,B0,T with
 |B0|>=M/(2K)>=N/(8p^gamma),B0 subsetb0T,|T|<=32768K^18M.
The subsetB0 qualifies for(c)--(d),sincep^(1/16-gamma)>8.
Also32768K^18M<=2^33p^(18gamma)N<p^(1/4)N,
because18*256/20481+33/3200000<1/4. Thus
 |T|<p^(7/3-1/16+1/4)=p^(121/48)<p^3.
An additive subgroupofF_p^4 is anF_p-vector subspace,so dimT<=2.
Unit scaling preserves pairedprojectioncardinalities. By(c),each
projection ofT has size at leastp^(5/4),hence dimensionat leasttwo.
ThusdimT=2 and both projections bijectontoF_p^2. T is thegraph of
a unitalF_p-linear ringisomorphismphi:F_p^2->F_p^2. Imagesof(1,0),(0,1)
are nonzero complementaryidempotents,so are thatpair in eitherorder.
F_p-linearity makesphi identityorcoordinateswap. Forb0=(b1,b2,b3,b4),
allbi nonzero,b0T is thus oneof
 {(x,y,(b3/b1)x,(b4/b2)y):x,y inF_p},
 {(x,y,(b3/b2)y,(b4/b1)x):x,y inF_p}.
These are the forbiddenplanes(d),yet qualifyingB0 liesinb0T.
This contradiction proves(19m)withallsubsetquantifiers.

The selected family uses eta=256000000/2109543000001, gamma=256/20481.
The CZ cap eta<=1/8000 holds; 103eta<gamma follows from
2109543*1000000<2109543000001. Theta=eta/8=32000000/2109543000001.
The separate explicit guard u>=3200000
in the advice map ensures actual p>2^u meets this theorem's threshold.
The weaker eta*u>=64 guard is not used to infer that prime threshold.
AllremainingfiniteCZ/uniform-sourceconditions retain their suppliedproofs.
The new threshold is proved directly,not inferredfrom a weaker theorem.

# Scalar sum-product and one support-growth step

Section SP6 proves the scalar sum-product estimate (30') used by Lemmas RG
and S1' in Section SP5s. No somewhere-condenser or scalar entropy-growth
theorem is assumed. All
logarithms are binary, TV is half-L1, and H_inf is min-entropy. Statements
with real entropy thresholds mean literal atom bounds, so integer ceilings are
not silently omitted.

## SP5s. One support-growth step

This section replaces the original reduction to scalar entropy growth. It
uses the two branches of (30') in Section SP6.3, with weighted cost
16*2+18*2=68; the sumset inequality (6) and the
three-edge-walk lemma (7) of restricted-nonmalleability-and-breaking.md
Sections 3.1 and 3.2, cited below as (R6) and (R7) (not this file's (6) and
(7)), used here in (F_p,+) and (F_p^*,x); and prime-fields Sections 1-2,
recalled in Section SP5s.2. U_D is uniform on D.

Constants: e=eta_RG=10/711, xi=1/427, r=(1-2xi)/(1+2xi)=425/429.
The two separate branch margins and the support margin are
 1-71e=1/711, r-34e=52105/101673, e/2-3xi=2/303597.      (RG0)
In particular 800000/711>252, 800000*52105/101673>122, and
300000000*(2/303597)>17+60/711. These are rational comparisons.

### SP5s.1 Lemma RG

**Lemma RG.** Let p be prime and M an integer with log M>=800000 and
p>M^(2/(1+2xi)). Let A1,A2,C0 be subsets of F_p^* of size exactly M and
X1,X2,X3 independent and uniform on them. Then for every W subset F_p with
|W|<=4M,

    P[X1*X2+X3 in W] <= 2M^(-e).                                      (RG)

Proof. Call a nonempty D subset F_p^* plus-good if P[U_D+U_V in W]<=M^(-e)
for every V subset F_p^* with |V|=M and every W subset F_p with |W|<=4M,
and times-good if P[U_D*U_V in W]<=M^(-e) for all such V,W. A disjoint
union of times-good sets is times-good: for fixed V,W its hitting
probability is the size-weighted average of the pieces'. Start with R=A1
and repeat: stop if R is plus-good or |R|<=M^(1-e); otherwise find a
times-good D subset R (below) and move D out of R. Each round removes a
point, so the procedure stops; it is a proof device, never evaluated.

Finding D. R is not plus-good: some V,W as above have
P[U_R+U_V in W]>alpha:=M^(-e), so the bipartite graph on R x V with edges
{(x,v):x+v in W} has more than alpha|R|M edges. (R7) gives D subset R and
C' subset V with |D|>=alpha|R|/8, |C'|>=alpha M/8, each (x,v) in D x C'
joined by at least alpha^5|R|M/2^18 walks x,v1,x1,v. A walk gives
x+v=s1-s2+s3 with s1=x+v1, s2=x1+v1, s3=x1+v in W; given (x,v) the triple
fixes v1=s1-x and x1=s2-v1, and the triple fixes the sum. So distinct sums
use disjoint sets of triples, and |D+C'|alpha^5|R|M/2^18<=|W|^3<=64M^3.

Suppose D is not times-good: some V',W' have P[U_D*U_V' in W']>alpha. The
edges {(x,v):xv in W'} have products in W' minus 0, a subset of F_p^*, and
the same counting in (F_p^*,x) (xv=s1*s2^(-1)*s3) gives D' subset D and
B' subset V' with |D'|>=alpha|D|/8, |B'|>=alpha M/8 and
|D'B'|alpha^5|D|M/2^18<=64M^3. As |R|>M^(1-e), |D|>M^(1-2e)/8 and
|D'|>M^(1-3e)/64, (R6) additively (with D' subset D) and multiplicatively
give, for K+=|D'+D'|/|D'| and Kx=|D'D'|/|D'|,

    K+<=|D+C'|^2/(|C'||D'|)<=2^54 M^(3+12e)/(|R|^2|D|)<2^57 M^(16e),
    Kx<=|D'B'|^2/(|B'||D'|)<=2^54 M^(3+12e)/|D|^3<2^63 M^(18e).

N'=|D'|>=2: (1-3e)log M>=800000*(681/711)>7, and
N'>M^(1-3e)/64. Also N'<=M and p/N'>=p/M>M^r.
Put X=K+*Kx<2^120 M^(34e) and apply the two branches of Section SP6.3.
If p>=3XN', then N'<=33X^2<33*2^240 M^(68e), so
 M^(1-71e)<64*33*2^240<2^252,
contradicting (RG0). If p<3XN', then
 M^r<p/N'<3X<3*2^120 M^(34e)<2^122 M^(34e),
contradicting the second margin in (RG0). This treats both branches
without weakening the lower bound on N' by a factor r in its exponent.
Thus D is times-good.

Distribution. Let G be the times-good union and R the final remainder.
Given X1 in G and X3=c0, P[X1X2 in W-c0]<=M^(-e) (|W-c0|=|W|); average over
c0. Given X1 in R with R plus-good, fix X2=a2: X1a2+X3=a2(X1+a2^(-1)X3),
where a2^(-1)C0 is a size-M subset of F_p^* and |a2^(-1)W|=|W|, so the bound
is again M^(-e); average over a2. If instead |R|<=M^(1-e), the R branch
has probability at most M^(-e). Empty branches contribute zero. The
G contribution is <=M^(-e); the remainder contributes <=M^(-e),
proving (RG).

### SP5s.2 The one-step map and its support lemma

For an integer m>=1 let p_m be the least prime above 2^m; prime-fields
Sections 1-2 (restricted-nonmalleability-and-breaking.md, first part) prove
2^m<p_m<2^(m+1) and find p_m by trial division in O(2^(3m/2+2)(m+1)^2) bit
operations. Read an m-bit string as a binary integer in [0,2^m), hence as an
element of F_(p_m) (injective, since 2^m<p_m). Define

    E_m(x1,x2,x3) = (x1*x2+x3 mod p_m) mod 2^m,                        (S2)

the low m bits of the canonical representative in {0,...,p_m-1}. For alpha
in {0,1}^(3m) write alpha=(alpha1,alpha2,alpha3) with m-bit blocks; its four
ROWS are alpha1, alpha2, alpha3 and E_m(alpha).

**Lemma S1' (support growth).** Let m>=300000000. For every Q subset
{0,1}^(3m) with |Q|>=2^(3m/2-15), the union U of the four row images of Q
(a subset of {0,1}^m) satisfies |U|>2^((1/2+xi)m).

Proof. Suppose N=|U|<=2^((1/2+xi)m). As Q is a subset of U^3,
N>=2^(m/2-5), so M:=|U\{0}|>=N/2>=2^(m/2-6) and log M>=m/2-6>=149999994
>=800000. Also M<=N<p_m^(1/2+xi) (as 2^m<p_m), so p_m>M^(2/(1+2xi)). Let
X1,X2,X3 be independent and uniform on U, viewed in F_(p_m). Lower bound: on
Q the fourth row lies in U, so P[E_m(X) in U]>=|Q|/N^3>=2^(-15-3xi*m).
Upper bound: E_m(X) lies in U exactly when Y=X1X2+X3 lies in
U*={y<p_m: y mod 2^m in U}, and |U*|<=2N since p_m<2^(m+1). P[some Xi=0]
<=3/N, and conditioned on no zero the Xi are independent and uniform on
U\{0}, of size M. The hitting-probability Lemma RG with |W|=|U*|<=2N<=4M gives
(the conservative upper bound below weakens its factor2):

    P[E_m(X) in U] <= 3/N+2M^(-eta_RG) <= 3M^(-eta_RG) <= 4M^(-eta_RG),

since 3/N<=3/M<=M^(-eta_RG) (M^(1-eta_RG)>=3). So 4M^(-eta_RG)>=2^(-15-3xi*m), that
is, eta_RG*log M<=17+3xi*m, and with log M>=m/2-6: (eta_RG/2-3xi)m<=17+6eta_RG. But
(eta_RG/2-3xi)m=2m/303597>=600000000/303597>18>17+6eta_RG.
Contradiction.

Only set counting and independent uniform variables on U appear; no
independence of the rows of an actual source is used. sampler-and-
disperser.md Sections 3-4 use Lemma S1' with the explicit guard
m=P3*B0*L>=300000000, automatic at the outer L>=2^40.

## SP6. Scalar sum-product from polynomial incidences

We prove every incidence estimate used here. The argument is adapted from
M. Lewko, *A Szemeredi-Trotter Theorem in Arbitrary Fields*,
arXiv:2609.27023v3, Lemma 2 and Section 4, equation (14), and Elekes's
point-line construction. No theorem from either source is a premise.
Write t_+=max(t,0); characteristic zero is denoted p=infinity.

Algebra: polynomial rings over fields are UFDs by induction in the variables.
Over a UFD, primitive products are primitive: reduce modulo each coefficient
prime and use that polynomial rings over domains have no zero divisors.
Clearing denominators transfers unique factorization from the fraction field's
Euclidean univariate ring (Gauss's lemma), and preserves primitive
irreducibility over that field. Leading forms show that total and
group-variable degrees add on nonzero products. Successive division by
distinct linear factors bounds a degree-j univariate polynomial's roots by j.

### SP6.1 Lemma M: multiplicities on parallel lines

Let K be a field, 0<=d<p an integer, and F in K[x,y] nonzero and squarefree
(no repeated irreducible factors), with total degree <=d. Let lambda_i be
finitely many distinct parallel lines, with q_i on lambda_i, such that
F restricts to a nonzero polynomial on every lambda_i. For a nonzero common
direction v let mu_i be the order of u in F(q_i+uv). Then

    sum_i (mu_i-1)_+ <= d(d-1).                              (M)

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
Take T also at least all their degrees. Thus codim W_T=sum h_i.
Every multiple of G or G_X lies in this kernel, so
 J_T=G Pi_(T-d0)+G_X Pi_(T-d1) is contained in W_T.
Coprimality and degree additivity give the intersection of these summands
as G G_X Pi_(T-d0-d1). Multiplication is injective. Hence

 codim J_T=binom(T+2,2)-binom(T-d0+2,2)-binom(T-d1+2,2)
             +binom(T-d0-d1+2,2)=d0*d1.

Thus sum h_i=codim W_T<=codim J_T=d0*d1<=d(d-1). This includes d1=0.
The restriction d<p is essential: x^p-y has zero x-derivative and order p
on every horizontal line through (xi,xi^p).

### SP6.2 Lemma I: the d-form incidence bound

Let P be m points of k^2 and L be n distinct nonvertical lines over a field
k of characteristic p. For every integer 1<=d<p put v=floor(2n/(d+2)). Then

 I(P,L)<=m(v+1)+d(d-1)+n<=m+n+2mn/(d+2)+d(d-1).            (I)

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
distinct slopes. Root counting gives c_(q,j)=0 when j<r(q)-v. Thus
u^((r(q)-v)_+) divides R_q in k[z,u].

Work over K=k(z), keeping z indeterminate, and set
lambda_q={(q1+u,q2+zu):u in K}. These parallel lines are distinct:
equality of their intercepts q2-z*q1 forces equality of both coordinates.
S stays squarefree in K[x,y]. Indeed each irreducible factor involving
x or y is primitive over k[z], hence irreducible over k(z) by Gauss's
lemma. Distinct such factors f,g cannot become associates: clearing
denominators would give af=bg for nonzero a,b in k[z]; primality of f and
f not dividing b would force f dividing g. Factors in z alone become units.

Let E={q:R_q=0}, h=|E|. On the remaining lines apply (M), using
mu_q>=(r(q)-v)_+. Their incidences total at most
 (m-h)(v+1)+sum_(q notin E)(r(q)-v-1)_+
 <=(m-h)(v+1)+d(d-1).
For q in E, substitute u=x-q1 in R_q=0. Division in y by the monic
linear polynomial l_q=y-q2-z(x-q1) shows l_q divides S in k[x,y,z].
These factors are irreducible (their quotient ring is k[x,z]), pairwise
non-associate, and of z-degree one. Their product divides S, so h<=v.
A given line containing s points of E contributes s<=1+binom(s,2);
each pair of distinct points lies on at most one line. Thus
 I(E,L)<=n+binom(h,2)<=n+h(v+1).
Adding gives (I). All zero restrictions above are polynomial identities,
not merely pointwise zeros on a finite field.

### SP6.3 Elekes theorem, with both characteristic branches

For A subset F_p^*, N=|A|>=1, let Kplus=|A+A|/N, Ktimes=|AA|/N and
X=Kplus*Ktimes>=1. Then

 p>=3XN implies N<=33X^2; p<3XN implies p/N<3X;
 min(N,p/N)<=33 Kplus^2 Ktimes^2.                         (30')

Proof. Take P=(A+A) times (AA), of size m=XN^2, and the N^2 lines
L={y=a(x-b):a,b in A}. They are distinct: their slopes determine a and,
since a!=0, their intercepts determine b. For each a,b,c in A the point
(b+c,ac) lies on its line, and distinct c give distinct points on that
line. Thus I(P,L)>=N^3. If p>=3XN, choose d=ceil(3XN)-2. Because XN>=1
and p is an integer, 1<=d<=p-2 and d+2>=3XN, while d<3XN.
Applying (I) gives
 N^3<=(X+1)N^2+(2/3)N^3+9X^2N^2,
so N<=3X+3+27X^2<=33X^2, using X>=1. Otherwise p/N<3X<=33X^2,
proving both branches and (30'). No deletion of zero is hidden: A is
assumed nonzero; Lemma RG applies it to D' subset F_p^*.

## SP8. Sources

Barak, Impagliazzo, Wigderson, *Extracting Randomness Using Few Independent
Sources*, Sections 3.1--3.4 and Appendix A:
https://www.math.ias.edu/~avi/PUBLICATIONS/MYPAPERS/BIW04/biw.pdf
(the AB+C growth step and the additive/multiplicative peeling architecture);
The finite BSG and sumset tools are proved in the companion module, with
attribution there. The hitting-probability goodness and additive-first
peeling architecture adapt published b15a467b-c75f-4152-befe-2b634e34a474
and the exact-size support route of mpjunior92's f7447a6f-a996-404e-bb7f-df63d9488837.
Section SP6 adapts Lewko's multiplicity and incidence proofs (v3, 7 Oct 2026):
https://arxiv.org/html/2609.27023v3 ; source license CC BY 4.0.
The final construction is due to G. Elekes, *On the number of sums and
products*, Acta Arith. 81 (1997), 365-367. The cost-68 branch calibration
and the explicit constants here are derived in this package. The preprint's
Lean formalization is not assumed. The earlier ring/CZ growth argument
is unchanged; eta_RG is distinct from its parameter eta_CZ.
