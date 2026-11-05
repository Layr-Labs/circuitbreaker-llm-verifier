# Proof component: compact-affine-wrappers.md

Original finite proof component. Local component filename references
resolve through proof.md Section 0.

<!-- COMPONENT compact-affine-wrappers.md SHA256 dd69a123e00a3c5c456c880219e457d27edc08e43cf469713b78856171625758 -->

# Finite affine transformation and breaker wrappers

Original proof of finite reductions, importing the supplied local primitive
proofs. Global parameter assembly is separate. An independent exact-version
audit is required before consolidated evidence claims. Originals and M=1000
are unchanged. Logs are binary, TV is half L1, + is XOR, U_b is independent
uniform, D(V|H)=TV((V,H),(U_|V|,H)), and
Havg(X|H)=-log sum_h P(h)max_x P(x|h). Only positive-probability conditional
branches are used. Lengths are integers; entropy bounds are literal atom caps.

Attribution: CL's affine transformation (Theorem 6.1),
https://arxiv.org/abs/2110.12652v1 ; Li's generalized copies, advice reduction
and equal/general-length wrappers (Definitions 3.3--3.4, Lemmas 4.2,5.1,
Theorem 5.2), https://eccc.weizmann.ac.il/report/2023/023/revision/1/download/ .
BKSSW condensation: https://eccc.weizmann.ac.il/report/2010/037/download .
CZ's restricted polynomial core:
https://eccc.weizmann.ac.il/report/2014/102/download/ . The local replacements
are credited in their modules. No third-party proof prose is copied.

## 1. Imports and elementary tools

Imports mean the actual algorithms with ALL finite guards in seeded-linear.md,
near-full-seeded.md, weak-seed-hashing.md, bilinear-extraction.md (in
/tmp/circuit-log-primitives-20261105T0130/), and the adjoining
constant-advice-breaker.md, ordinary-condenser.md, cz-conditional-entropy.md,
finite-approximate-ring.md, prime-fields-and-encoding.md. A self-contained
consolidation must include those proofs once; none is an unrestricted
asymptotic black box.

### 1.1 Probability facts

If full tuples L,R are independent given H, revealing a function of one side
and H preserves the conditional product: normalize only that factor in
P(L,R|H)=P(L|H)P(R|H). This permits alternating rectangular histories and
arbitrary correlation inside either side. A mixed XOR becomes left-sided
after its right contribution is fixed.

A T with at most 2^b values obeys
sum_t P(t)max_x P(x|t)=sum_t max_x P(x,t)<=2^b max_x P(x).
The same conditional inequality gives Havg loss at most b. Markov shows
Havg(X|H)>=K+s implies Hinf(X|H=h)>=K outside mass 2^-s.
A right-only revelation does not change conditional left laws; alternatively
charging the entire history length is always safe. These statements allow
randomized T.

TV with a common history marginal is average conditional TV. A common kernel
contracts TV by the half-L1 triangle inequality. Consequently joint laws at
distance zeta obey D_P(V|H)<=D_Q(V|H)+2zeta: one term restores H's marginal.
If a coupling changes only V and keeps H fixed, disagreement costs only once.

The common mass min(P(z),Q(z)) and residual masses give a coupling with
P(Z!=Z*)=TV(P,Q), realizable by a kernel from already sampled Z. If Z=f(A)
and auxiliaries have law K(aux|Z), this preserves the original A law and
Law(A|Z,aux)=Law(A|Z). An independent B stays independent of the enlarged
left side. Right-only revelations with auxiliaries fixed preserve that identity,
by factoring the joint density. Labels need not have small ranges.

If S=f(L), R is independent of L, and W=F(S,R), then D(W|L)=D(W|S):
the conditional W law depends on L only through S. The conditional-product
version permits a seed-retaining extractor to retain the entire left side.

### 1.2 Strong merging and robustness

Let source tuple (U,U_1,...,U_t) and seed tuple (S,S_1,...,S_t) be independent
given H. Let F be strong (K,eps), output b bits, 2^-q<=eps. For I,J subset
[t], suppose

  D(S|H,S_I)<=eta,
  Havg(U|U_J,H)>=K+|I\J|b+q.                              (1)

With V=F(U,S), V_j=F(U_j,S_j),

  D(V|V_(I union J),S,S_1,...,S_t,H)<=eta+2eps.             (2)

Proof: reveal U_J,S_I, then V_(I\J). The latter are source-only with S_I
fixed, costing |I\J|b; the product split survives. Entropy fails on at most
eps mass. Uniform-seed extraction costs eps elsewhere. On each such history,
conditional output-distance with the seed exposed is a [0,1]-valued function
of S; changing its marginal to uniform costs its conditional seed distance.
Source revelations are independent of S given H,S_I, so their averaged cost
remains eta. Append remaining seeds by their conditional kernel given H,S_I,S,
independently of the source. Outputs in J are determined by U_J and the retained
seeds. Discard surplus history to prove (2).

The same expectation argument for an all-seeds-retained breaker of error e
gives error e+zeta when its good seed marginal is zeta from uniform: preserve
the kernel of other seeds given the good one. Average source entropy K+q adds
at most eps for bad histories. Thus a strong ordinary breaker of error eps
has robust error eta+2eps. Extra seed copies can be retained when independent
of the full source pair given the retained seed pair/history, as in Section 2.

## 2. Generalized all-seeds-retained affine transformation

### 2.1 Finite contract and algorithm

Fix n,d,t,a,m>=1, eps in (0,1), q=ceil(log(1/eps)), H=ceil(log t).
Supply these deterministic maps, each of error <=eps:

  C:  {0,1}^d x {0,1}^d0 x {0,1}^a -> {0,1}^dx;
  L0: {0,1}^n x {0,1}^s0 -> {0,1}^d0, threshold K0;
  Lr: {0,1}^n x {0,1}^dx -> {0,1}^r, threshold Kr;
  E:  {0,1}^d x {0,1}^dy -> {0,1}^dx, threshold KE;
  Lm: {0,1}^r x {0,1}^dx -> {0,1}^dy, threshold Km.

Every L is strong and linear in its source for each fixed seed; E is strong.
C is a strong ordinary one-tampering breaker of threshold KC: arbitrary
correlated source pair independent of arbitrary correlated seed pair, honest
seed uniform, distinct advice, both complete seeds retained. Its first input
really has length d. Padded seeds may ignore suffixes while retaining them.
All actual primitive guards must hold. Put

  T0=(t+1)(s0+2d0+dx), T=T0+H(t+1)(4dy+2dx).

Require positive legal dimensions, k<=n, and

  s0<=d, dy<=r, m<=r,
  r>=(t+1)dy+Km+q,
  k>=max(K0,T+Kr+t*r+q),
  d>=(t+1)(s0+d0)+KC+q,
  d>=T+KE+t*dx+q.                                         (3)

The evaluator on (x,y,alpha) is

  s=Prefix_s0(y); r1=L0(x,s); s2=C(y,r1,alpha); w0=Lr(x,s2).
  For i=1,...,H:
    p=Prefix_dy(w_(i-1)); qm=E(y,p); v=Lm(w_(i-1),qm);
    qr=E(y,v); w_i=Lr(x,qr).
  Return Prefix_m(w_H).

The chosen assembly has r=(t+5)w and w>=m, so no long-output extension is used.
Let X_j=A_j+B_j, j=0,...,t. Suppose the arbitrary correlated A tuple is
independent of the arbitrary correlated (B,Y) tuple, Hinf(A_0)>=k, Y_0 is
uniform d-bit, and alpha_0 differs from every alpha_j. Other components need
no entropy, uniformity, or functional relationship to the honest components.
For F the displayed evaluator,

  D(F(X_0,Y_0,alpha_0) |
    (F(X_j,Y_j,alpha_j))_(j=1)^t,Y_0,...,Y_t)
      <=delta_H=(13*2^H-8)eps<=26t eps.                    (4)

Distance can be capped at 1. We retain the original seed without adding a
synthetic tampering. Proof objects below are never computed by the evaluator.

### 2.2 Initialization

Write L=A tuple, R=(B,Y) tuple. Fixed-seed linearity splits every L call
into A and B contributions, e.g. R1_j=R1A_j+R1B_j. Reveal all
S_j=Prefix_s0(Y_j), then R1B_j; call this H0, length <=(t+1)(s0+d0).
Strong extraction and independence of A_0 from the entire right tuple give
D(R1A_0|all S)<=eps. Given S, appending R1B is independent of A_0; translate
by R1B_0 to obtain D(R1_0|H0)<=eps. The entire R1 tuple is left-sided given
H0 and independent of the Y tuple.

Y_0 has average entropy >=d-(t+1)(s0+d0)>=KC+q. Section 1.2 applied to C
therefore gives D(S2_0|S2_j,H0,R1_0,R1_j)<=3eps for every j>=1.
Other R1 coordinates can be appended because, given H0 and the retained pair,
they are independent of the full Y tuple and hence both C outputs. Reveal
all R1, a left revelation since R1B was fixed. Then S2 is right-sided.
Apply (2) to Lr with sources A, seeds S2, I={j}, J empty. The k bound in
(3) supplies Kr+r+q after the history, giving error 5eps. Let Z0 contain
S,R1B,R1,S2 across all copies, of length <=T0. The full sides remain
independent and

  D(W0A_0|W0A_j,Z0)<=5eps.

### 2.3 Doubling invariant and full correlated-copy proof

Let u_i=min(2^i,t). After round i we maintain a rectangular Zi of length
<=T0+i(t+1)(4dy+2dx), with WiA left-sided, WiB right-sided, and for all
J subset [t] of size <=u_i,

  D(WiA_0|WiA_J,Zi)<=delta_i, delta_0=5eps.                 (5)

Assume it at i-1, denote its transcript Z, error delta, and outputs WA_j+WB_j.
For |J|=u_i choose J1,J2 of size u_(i-1) with union J; they may overlap.
Smaller sets follow by projection. Reveal all PB_j=Prefix_dy(WB_j), right-only.
Set P_j=Prefix_dy(WA_j)+PB_j. By (5), independent right retention and known
translation, D(P_0|P_J1,Z,PB)<=delta. With sources Y and seeds P, apply (2)
to E, I=J1,J empty. Its error is delta+2eps retaining all P. Reveal P and
call the resulting history H1; Qm is now right-sided.

Couple WA_0 to W* with (W*,WA_J2,Z) having independent uniform first
coordinate and P(W*!=WA_0)<=delta. A conditional maximal coupling given
WA_J2,Z is sampled from left information only, preserving the product split.
Given WA_J2,Z,PB the P coordinates in J2 are determined. Revealing the rest
costs at most (t+1-|J2|)dy, so

  Havg(W*|WA_J2,H1)>=r-(t+1-|J2|)dy.                       (6)

Apply (2) to Lm with source tuple (W*,WA_1,...,WA_t), seeds Qm,
I=J1,J=J2. Its seed error is delta+2eps. Its entropy guard follows because
(t+1-|J2|)+|J1\J2|<=t+1 and r>=(t+1)dy+Km+q. The output on W* has error
delta+4eps against the tampered Lm outputs, retaining H1,Qm. Replace only
this good output by the one on actual WA_0. History and tampered outputs
are identical, so disagreement adds delta once. Thus VA_0=Lm(WA_0,Qm_0)
has error 2delta+4eps against VA_J with H1,Qm retained.

Reveal Qm, a right revelation. Append the full right side by conditional
independence and translate by VB_j=Lm(WB_j,Qm_j), at zero error cost.
Reveal VB. Apply (2) to E(Y_j,V_j), I=J,J empty: error 2delta+6eps with
all V retained. Reveal V, left-sided, and apply (2) to Lr(A_j,Qr_j),
I=J,J empty: error 2delta+8eps. Reveal Qr, right-sided, to obtain Zi.
The new messages (PB,P,Qm,VB,V,Qr) total (t+1)(4dy+2dx) bits. Every E
call is covered by d-T>=KE+t*dx+q; every Lr call by k-T>=Kr+t*r+q.
The Lm entropy came separately from (6), not from the original A entropy.

This proves delta_i=2delta_(i-1)+8eps, hence delta_i=(13*2^i-8)eps,
including non-power-of-two t. At i=H, append the entire right tuple to the
bound against all tampered left outputs using the conditional product given
ZH. Translate by their right contributions, discard excess history, and take
prefixes. This proves (4), retaining every Y_j. No equality between different
A_j or B_j, or determinism of tampered copies, entered the argument.

### 2.4 Actual primitive use

Use seeded-linear.md's exact entropy recurrence/seed length for the L maps
and, if desired, E, with output truncation or ignored uniform seed suffixes.
Its K_0<=input-length guard and finite constants remain mandatory. For the
inner fixed-row F use constant-advice-breaker.md as C at actual source length
d and seed length d0. For advice a it has 2^a rounds and error
(8*2^a+2)eps_C, so choose eps_C to make this <=eps and check its actual
thresholds. This exponential advice dependence is used only at fixed inner
advice. The later Section 6 breaker supplies the outer growing-advice use.
The assembly is separate and cannot infer these finite tests from big-O claims.
The conditional-average lift needed later is proved in
compact-sampler-disperser.md Section 6.1, not assumed from marginal seed
uniformity.

## 3. Finite advice-condenser reduction

### 3.1 Contract

Fix ell>=2,r,s>=1, b=ell*s, n=ell*r, 1<=a<=b, 2^s>=ell, n>=b.
Supply N on ell independent (r+a)-bit inputs, output h<=r bits, such that
inputs each of entropy >=r-2b and coordinate maps fixed-point-free on their
honest supports satisfy

  TV((N(U),N(f(U))),(U_h,N(f(U))))<=eN.                     (7)

Put B=(ell+9)b, choose delta in (0,1), and require

  h-B-log(1/delta)>=s.                                    (8)

C_adv(x,alpha), x=(x1,...,xell), outputs the ell raw r-bit blocks followed
by N(x1||alpha,...,xell||alpha), zero-padded to r bits. For arbitrary
correlated X,X' with Hinf(X)>=n-b and alpha!=alpha', its matrix pair is within

  eA=(ell+2)2^-s+2^-b+delta+2^B eN                         (9)

of a convex mixture in which each component has a fixed row i satisfying
Hinf(A_i|A'_i=u)>=s for every u. Both matrices may change in the nearby law.

### 3.2 Flat, light, and dense pieces

A law capped by 1/M for integer M is a convex mixture of uniform M-point
sets: in the capped probability simplex, an extreme point cannot have two
coordinates inside their bounds (perturb them oppositely), or exactly one
(the sum constraint excludes it). Its vertices have exactly M entries 1/M.
Independently choose each f(x) according to Law(X'|X=x); this random function
realizes any conditional tampering law. By convexity it therefore suffices to
treat X uniform on S of size 2^(n-b) and fixed f.

Call (u,u') heavy in coordinate i if its pair probability is >=2^(-r-3s).
At most 2^(r+3s) such pairs exist. Partition points with a nonheavy coordinate
by their first such coordinate into T_i. Discard T_i of mass <2^-s, costing
ell*2^-s. In any surviving T_i joint pair atoms are <=2^(-r-2s); revealing
the r-bit partner leaves average conditional atoms <=2^-2s. Outside
conditional mass 2^-s the honest row therefore has entropy >=s. Split by
those partner values. Weighted failures over all surviving T_i cost <=2^-s.
If the all-heavy part has mass <2^-s, discard it and finish at total loss
(ell+2)2^-s. Otherwise its set SH has size >=2^(n-b-s).

Let Wi(u) be its possible heavy partners. Fewer than 2^(r-b-3s) coordinate
values have |Wi(u)|>2^(b+6s). Removing points with any such value discards
at most ell*2^(n-b-3s) points, a fraction <=ell*2^-2s<=2^-s of SH.
The remainder S* has size >=(1-2^-s)2^(n-b-s)>=2^(n-2b). Its coordinate
supports Si and Pset=product Si satisfy |Si|>=2^(r-2b), |S*|/|Pset|>=2^-2b.
Independently choose each gi(u) uniformly from Wi(u). Every (x,f(x)) in S*
is matched by coordinatewise g with probability >=2^(-(ell+6)b). While
more than 2^-b|S*| points remain, averaging finds a matching subset V of
size >=2^(-(ell+7)b)|S*|>=2^-B|Pset|. Remove it and repeat. This partitions
all but 2^-b|S*| points into dense pieces V, each with coordinatewise f=g.
The partition is a finite proof device, never an evaluator operation.

### 3.3 Conditioning near-independent output

If (W,T) is e-close to (U_h,T) and event E has probability p>0, there is a
law (W*,T*) within delta+e/p of (W,T)|E, with the same T|E marginal and

  Hinf(W*|T*=t)>=h-log(1/p)-log(1/delta).                  (10)

Indeed excess(w,t)=[P(w,t)-2^-h P(t)]_+ sums to at most e. Values t with
P(t|E)<delta P(t) have conditional mass <=delta. For the other t,

 P(w|t,E)<=2^-h/(p delta)+excess(w,t)/(p P(t|E)).

A probability vector's distance to atom cap c>=1/domain-size equals its
excess above c: remove that mass and distribute it below the cap, where
capacity suffices. Apply c=2^-h/(p delta), capped at 1 if necessary, average
over good t, and replace bad t with any capped law. This proves (10), even
if E is not determined by (W,T).

For each V apply (7) to independent uniforms on Si with alpha appended.
Extend gi arbitrarily off Si. The maps u||alpha -> gi(u)||alpha' have no
fixed points. If (7) requires
global fixed-point-freeness, map every input of suffix other than alpha to
0^r||alpha. Condition on V, probability p>=2^-B. Equations (8),(10) give
s conditional entropy in the final row, at cost delta+2^B eN. Pad that row
and append all other rows through their original conditional kernel given
the final-row pair. TV and the selected-row property are preserved. Replace
discarded pieces by any pair with a uniform honest selected row. The losses
are ell*2^-s for small T_i, <=2^-s for light conditional failures,
<=2^-s for heavy fibers, <=2^-b for the remainder, and the smoothing cost.
This proves (9).

### 3.4 Actual restricted ten-source input

Use ell=10, eta=1/10000000, u=r+a, rho=(a+2b)/u, and require

  rho<=eta/2, r-2b>=(1-eta)(u+1)+1,
  eta*u>=64, h=floor(eta*u)>=1, h<=r.                     (11)

The actual CZ polynomial map on the first eight injectively embedded u-bit
inputs is the one in cz-conditional-entropy.md Section 2. Select p as the
least prime above 2^u and encode its pair output W in nF=2(u+1) bits.
finite-approximate-ring.md supplies precisely E_p(1/10000) for p>=2^20000;
(11) gives p^eta>=2^64 and the stronger prime threshold, while
1000eta=1/10000. The required integer source threshold
log ceil(p^(1-eta)) is implied by (11). Distinct appended advice makes the
first eight coordinate maps fixed-point-free, the restricted theorem's
hypothesis. Its pair (W,W') is zeta<=45p^-eta from a law (W*,W'*) with
Hinf(W*|W'*)>=kappa=4eta log p on every branch.

Pad T=(U9,U10) with two zeros to nF bits, and define N=IP_nF,h(W,T) using
bilinear-extraction.md. T has entropy >=2(1-rho)u and is independent of
(W,W'). Couple to the smooth pair independently of T and condition on W'*.
Strong IP retaining T gives uniformity jointly with W'*,T, hence with the
actual deterministic function IP(W'*,T') of them. Its error is
(1/2)sqrt((2^h-1)2^(nF-kappa-2(1-rho)u)). Restore the original pair and
its tampered-output target, paying zeta twice. Since
kappa+2(1-rho)u-nF>=3eta*u-2 and h<=eta*u, (7) holds with

  eN<=91*2^(-eta*u).                                      (12)

The factor 1/2 in IP cancels the factor 2 from its square root. Only this
fixed-point-free case is claimed; no generic unchanged-input nonmalleability
is needed. The deep proofs and exact field algorithms are the named imports.
The advice amplification 2^B and guard (8) remain mandatory.

## 4. Imported ordinary and weak-seed contracts

Use ordinary-condenser.md Sections 1--5 at fixed rational 0<beta<9/10,
tau=9/10, nu=1/10. Use its ACTUAL fixed integer schedule and evaluator,
including both finite guards (12)--(13) there; below them its zero-row
fallback has no entropy guarantee. Write its computed parameters as
ellC,TC,gamma,m*, D=(ellC+1)^TC, q=ceil(r/ellC^TC). The input requirement
is kA>=beta*r and the exact mixture error is

  eC=TC*2^(-gamma*q/4)+D*2^(-q/10).                       (13)

Each mixture component has a fixed row of entropy 4q/5. To check the stronger
mixture interface: an entropy-(9/10)q selected row M_I conditioned on any
index of probability >=2^-q/10 has atoms <=2^-4q/5. All other indices total
<=D2^-q/10; replace their selected rows by uniform. This gives (13), rather
than incorrectly identifying selected-row marginal entropy with conditional
entropy. All scalar proofs, constants, row orders and guards are in that import.

Use precisely weak-seed-hashing.md's matrix map EH(Y,S), with

  m0=floor(q/64), eH=ceil(7q/32),
  BH=ceil(log(d*m0+1)), rH=BH+m0+2eH;
  m0>=1, 2rH<=q, m0+2eH<=floor(d/10).                     (14)

It uses 2rH bits of the q-bit seed and retains all q bits. For independent
S,Y with entropies >=4q/5 and >=floor(d/10), it gives

  D(EH(Y,S)|S)<=eHout=2^(q/5-eH)<=2^(-3q/160).             (15)

Only this map is selected. Its finite-field implementation is supplied locally.
A useful upper bound respects approximation: the support A of C(U_r) has
<=2^r matrices. Every mixture component gives A mass <=2^(r-4q/5), by
projecting to its good row. Thus, if eC<1,

  q<=(r+log(1/(1-eC)))/(4/5).                             (16)

Synthetic rows need not preserve actual support. Formula (16) does not replace
any direct guard (14).

## 5. Equal-length strong breaker

### 5.1 Finite statement and evaluator

Let (X,X') be independent of (Y,Y'), all d-bit, with honest entropies
kX,kY in [0,d], alpha!=alpha'. Choose h,v,r,q,m>=1, h<=d, v=10r<h.
On h-bit prefixes use bilinear-extraction.md with v output bits and require
kX-d+h,kY-d+h>=0. Its separate retained-input error is

  eI=(1/2)sqrt((2^v-1)2^(h-(kX-d+h)-(kY-d+h))).            (17)

Use Section 3 at ell=10 on the v-bit honest input, with advice capacity a,
good-row entropy kA=s and error eA from (9), checking every guard including
(8),(11). Uniform input satisfies its entropy v-b. Pad shorter advice
injectively. Use Section 4's C with kA>=beta*r, actual finite guards and
error eC. Let L=11D. Use EH and m0 from (14), checking its guards.

Supply F from Section 2 on d-bit sources and m0-bit seeds, m-bit outputs,
degree t=2L-2, row advice length >=max(1,ceil(log L)), error eF. Its inner
ordinary breaker is the old constant-advice map, never the present construction.
Choose positive sY,sW,sX,sZ and require

  KY=kY-2h-sY-m0-sW >= floor(d/10),
  KX=kX-2v-sX-m-sZ >= F's source threshold.                (18)

All F dimensions and guards (3) must hold. With fixed lexicographic rows and
injective fixed-width encode, the evaluator is

  V=IP(Prefix_h(x),Prefix_h(y)); R=C_adv(V,alpha);
  S=concatenation of C(R_i), i=1,...,11;
  W_j=EH(y,S_j); Z_j=F(x,W_j,encode(j)), j=1,...,L;
  Q=Z_1+...+Z_L.                                         (19)

The same map gives Q' on primed inputs. The claim retains the ACTUAL primed
output and the entire right pair:

  D(Q|Q',Y,Y')<=2(eI+eA+eC)+eHout+eF
                    +2^-sY+2^-sW+2^-sX+2^-sZ.             (20)

### 5.2 Source-preserving mixture couplings

Reveal T=(Prefix_h(Y),Prefix_h(Y')), a right revelation. The original left
pair is unchanged and independent of the conditional right pair; outside
mass 2^-sY, Hinf(Y|T)>=kY-2h-sY. Given T, V,V' are left functions.
The retained-right-prefix IP guarantee gives D(V|T)<=eI: adding the primed
prefix does not change the conditional left law. Couple V to Vbar uniform
given T via a kernel of (V,V',T), keeping V' and the whole original source
law unchanged. Disagreement is <=eI.

For each T apply C_adv to (Vbar,V'). Couple the matrix pair to its exact
mixture (A,A'), costing <=eA, and sample label M with selected row i=i(M).
Reveal M,A'_i. Then A_i has entropy >=kA on every branch, the conditional
right pair is unchanged, and the sides remain independent. All auxiliaries
are kernels of the original (V,V',T), not extra uncontrolled source leakage.

For every (T,M,A'_i), apply C to A_i. Couple only this honest D-row block
to the exact somewhere mixture, costing <=eC, leaving all other honest
blocks and every tampered block unchanged. Sample its label N. In the full
resulting matrices Sbar,Sbar', a selected row J has entropy >=4q/5 given
H0=(T,M,A'_i,N). Its partner Sbar'_J is fixed, a row of C(A'_i). The
full enlarged left side is independent of the right pair given H0; every
auxiliary is still generated by a kernel of (V,V',T). The three couplings
bound disagreement of the full real/synthetic row pair by eI+eA+eC. All
entropy tracking below refers to the unchanged real source law.

### 5.3 Right leakage and average weak-seed error

Reveal Wbar'_J=EH(Y',Sbar'_J), right-only since its row is fixed. The good
row law does not change. On each good T branch this m0-bit revelation leaves
Y entropy >=KY except on conditional mass 2^-sW. The two right bad sets
thus cost <=2^-sY+2^-sW. On every other (H0,Wbar'_J) branch, independence
and (14),(18) permit (15): Wbar_J=EH(Y,Sbar_J) is eHout-uniform retaining
Sbar_J. Section 1.1 lifts this to the entire enlarged left side, since the
conditional output depends on that side only through Sbar_J. The right pair
is not exposed at this stage.

Now reveal original (V,V') and the full left auxiliary record, but not X,X'.
Every synthetic row becomes fixed, making the complete (Wbar,Wbar') tuple
right-only. The original source pair and seed tuple are independent on each
branch. Conditional averaging of the retained-left guarantee bounds the
average distance of Wbar_J from uniform by eHout over earlier good branches.
There is no pointwise seed assertion and no extra seed-error tail loss.

### 5.4 Exact no-extra-leakage identity

Let G be the history now. Then

  Law(X,X'|G)=Law(X,X'|T,V,V').                            (21)

Indeed the density, before summing the right pair, factors as

 P_T(t)P_left(x,x')1_{(V,V')=v(x,x',t)}
   K(aux|v,t)P_right(y,y'|t)1_{Wbar'_J(aux,y')=w'}.

After t,v,aux,w' are fixed, summing y,y' gives a factor independent of x,x';
normalization yields (21). This, not a small label range, justifies the
auxiliaries' zero additional entropy charge. For each T the real pair V,V'
has <=2^(2v) values and X retains its original entropy kX. Hence except on
original (T,V,V') mass 2^-sX, Hinf(X|G)>=kX-2v-sX. The original marginal
is unchanged, so this remains the total mass after averaging auxiliaries and
right revelations. Do not condition the earlier weak-extraction step on this
X-good event; only charge its eventual bad branches.

Finally reveal Zbar'_J=F(X',Wbar'_J,encode(J)). Its seed and advice are fixed,
so this is a left-only m-bit revelation. On preceding good-X branches it
leaves entropy >=KX except on additional mass <=2^-sZ. Conditional products
persist. The right seed law on each pre-revelation branch is unaffected, so
its average distance bound remains eHout. Restricting to surviving branches
only decreases the weighted sum of these nonnegative distances.

### 5.5 Breaking, full right retention, and restoration

Use (X,Wbar_J,encode(J)) as the good F call. The 2L-2 tamperings are
(X,Wbar_j,encode(j)) and (X',Wbar'_j,encode(j)) for j!=J. All advice differs
from encode(J). They fit Section 2 with zero right additive components and
source copies X or X'; the entire source tuple is independent of the seed
tuple. F already retains all seeds at degree 2L-2. Its robust guarantee makes
Zbar_J uniform up to eF plus that branch's seed distance, jointly with all
other row outputs and seeds; Wbar'_J,Zbar'_J were fixed in the history.

Append (Y,Y') with no extra cost: on this branch seeds are deterministic
functions of that pair, sources are independent of it, and row outputs use
only sources and seeds. Thus Law(Y,Y'|all seeds,history) is independent of
row outputs, and is a common kernel for the real/ideal laws. Merely knowing
that side information contains seeds would not suffice without this argument.
XOR the honest outputs and separately the tampered outputs. The uniform good
output remains uniform after its known XOR shift. Average seed errors and
charge the four bad sets. The synthetic error is
eHout+eF+2^-sY+2^-sW+2^-sX+2^-sZ. Return to real evaluation at cost
<=eI+eA+eC for the joint law and the same again to restore the actual
tampered-output target. This proves (20).

## 6. Actual general-input standard breaker

Let (X,X') on n bits be independent of (Y,Y') on d bits, Hinf(X)>=d,
Y uniform, d<=n, and alpha!=alpha'. Use EXACTLY near-full-seeded.md at
source length n, entropy d, integer exponent E>=1: its five-call seeded-linear
map with h=11, e=b=E+4 and its complete computed seed length hs. Check its
finite guards d>=16(H+1), d>=1000P, with H,P defined by that module's
literal recurrence; each constituent call's legal threshold is proved there.
It outputs L0=floor(9d/10) bits with modern strong error eta=2^-E and hs<d/100.
This is the actual supplied primitive, not an optimal-extractor assumption.
Choose tau in (0,1), eL, and require

  hs+L0<=d, hs+log(1/tau)<=(1/10)L0.                      (22)

Instantiate Section 5 at length L0, honest entropies >=(9/10)L0, with every
finite guard and error (20)<=eL. The general map takes disjoint first hs and
last L0 slices Y0,B of Y (and the same primed slices), computes
A=Ext(X,Y0), A'=Ext(X',Y0'), then Q=CB_equal(A,B,alpha), likewise Q'.
Its output guarantee is

  D(Q|Q',Y,Y')<=2eta+tau+eL.                              (23)

Proof: reveal T=(Y0,Y0'). The left output pair is independent of the full
conditional right pair. Strong extraction gives average D(A|T)<=eta because
its conditional law depends only on Y0. B is uniform independent of Y0;
revealing the hs-bit Y0' therefore leaves entropy >=L0-hs-log(1/tau)
except on T-mass tau, by the entropy rule applied separately for each Y0.
On every T couple A to U_L0, preserving the kernel of A' given A. The
synthetic left pair is within conditional distance eta_T, averaging <=eta,
and remains independent of the unchanged right pair. Its honest input is
uniform; it need not still be an Ext evaluation. Apply Section 5 on good T.
Lift its retained B,B' to Y,Y' by the original conditional kernel given
T,B,B', independent of the synthetic left pair; outputs depend on the right
pair only through B,B'. Average, charge tau, and restore the real output law
and its actual tampered-output target at eta each. Drop T to obtain (23).
For eta<=eps/8,tau=eL=eps/4 this is <=3eps/4. The exact simultaneous choice
and all-length threshold are proved separately in breaker-parameter-assembly.md.

## 7. Scope

The dependency chain is acyclic: seeded primitives -> old constant-advice
breaker -> inner affine breaker; CZ+IP -> advice condenser; ordinary condenser
and weak hash plus inner affine -> equal-length breaker; near-full plus equal
length -> general breaker, available for a later outer affine call. Never use
the latter to instantiate its own inner breaker. Evaluators never compute proof
labels or partitions. Imports with 2^O(r) setup are polynomial in outer N only
after a fixed r=O(log N) bound. Global assembly, exceptional-length fallback,
outer sampler/family/circuit proof, package fit and independent exact-version
audit are outside this module; no final-candidate or external-review claim follows.
