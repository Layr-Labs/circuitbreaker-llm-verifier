# Finite affine transformation and breaker wrappers

Logs are binary, TV is half L1, + is XOR, U_b is independent
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

Imports are the supplied finite algorithms and proofs mapped in proof.md
Section 0, with all guards; no asymptotic black boxes are used.

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

The chosen assembly has r=Rw, R>=t+5 and r>=m by (A2), so truncation suffices.
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
advice. The standard breaker NB of Section 4 supplies the outer growing-advice use.
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
bilinear-extraction.md. T has entropy >=2(1-rho)u. The full pair (T,T')
is independent of (W,W') by the independent coordinate-pair inputs here.
Couple to the smooth pair independently of (T,T') and condition on W'*.
Strong IP retains T. Append T' through the common kernel Law(T'|T),
independent of the smoothed first-eight pair; this preserves TV distance.
Thus uniformity holds jointly with W'*,T,T', and IP(W'*,T') is a
deterministic function of retained variables. Its error is
(1/2)sqrt((2^h-1)2^(nF-kappa-2(1-rho)u)). Restore the original pair and
its tampered-output target, paying zeta twice. Since
kappa+2(1-rho)u-nF>=3eta*u-2 and h<=eta*u, (7) holds with

  eN<=91*2^(-eta*u).                                      (12)

The factor 1/2 in IP cancels the factor 2 from its square root. Only this
fixed-point-free case is claimed; no generic unchanged-input nonmalleability
is needed. The deep proofs and exact field algorithms are the named imports.
The advice amplification 2^B and guard (8) remain mandatory.

## 4. The standard breaker NB

Based on submission0165533a-47ba-4e5b-a0c8-3cc7605b3da0. The inner F now
uses the output-amplifying selector (A2); all generic coefficient bounds
below remain conservative upper bounds. The proof keeps left-outcome
revelations distinct from unrevealed right-side Ext2 repair randomness.

This section supplies the standard breaker of proof.md interface (3.1). It needs
a general source X and a UNIFORM seed Y only. The rate-theta advice row is used
as the SOURCE of the AGHP extractor of weak-seed-hashing Section 1, seeded by a
fresh slice of Y, so no rate-to-0.8 condensation and no ordinary condenser is
used. Inside this section F names NB's inner breaker, lg(x) is the least j>=0
with 2^j>=x, pow2ceil(x)=2^lg(x)<2x, and [1.1] cites Section 1.1.

### 4.1 Statement

Fix eta=1/10000000, theta=1/10^12 (proof.md Section 4), t0=21, a0=4, and
h0=49, the least h with 2^h>=480/theta (2^48<4.8*10^14<=2^49); B0=2^h0. Let F be
the inner breaker of Section 4.2. For integers n,z,a,E>=1, Section 4.3 computes
a seed length d_C and threshold k_C, and Section 4.4 defines a deterministic
NB: {0,1}^n x {0,1}^(d_C) x {0,1}^a -> {0,1}^z. If k_C<=n, then for every pair
(X,X') on n bits independent of every pair (Y,Y') on d_C bits, with
Hinf(X)>=k_C, Y uniform, copies arbitrarily correlated within each pair, and
alpha!=alpha',

    D(Q|Q',Y,Y') <= 15*2^(-E-4) < 2^-E,  Q=NB(X,Y,alpha), Q'=NB(X',Y',alpha').   (NB)

Moreover d_C<=k_C<=c_std[z+a+E+lg(n+1)+1] with c_std=489603919856837338526670350
(Section 4.7). So NB satisfies interface (3.1) with c=c_std at every a,z,E>=1,
with the single guard k_C<=n, retaining both complete seeds. If k_C>n the map
outputs 0^z and no claim is made (the totality convention of proof.md Section 3.3).

### 4.2 The inner breaker F

F is proof.md Lemma A (Section 3.2) at degree t=21, advice length 4, delta0=1,
kappa=36, with row multiplier R=513 and C the constant-advice breaker of restricted-nonmalleability-and-
breaking.md Section 8 at advice length 4 (M=16), run at exponent e+j0 with
j0=lg 130=8, so its error (8M+2)2^-(e+j0)<=2^-e. Then (A1) gives j=lg 21=5, R=513,
h_A=13 (2^12<12*513; 2^13>=64*38 and 2^13>=12*513), T_c=22*(2+2+30)=748,
K_c=max(3,748+23*513+2)=12549 and Q=max(514,22*2+36+1,748+25)=773. F's selector is the
least power of two w passing (A2) and the two further tests

    w>=128(e+j0),    w>=128*G^ca,                                      (F2)

where G^ca is restricted Section 8's G_u at exponent e+j0 and u=Qw. Then that
file's selector (15), at source length Qw, output w and exponent e+j0, returns
exactly w_ca=w/128. Any passing w_ca has 128w_ca>=w. And w/128 is a power of two
passing all three tests of (15): w/128>=e+j0, 128(w/128)=w, and
U=max(Qw,w,56w/128)=Qw gives G_U=G^ca<=w/128. By (16) there, its actual seed is
56w_ca=7w/16<w and its threshold is 36*128w_ca=36w<=Qw. Lemma A (remark (i)
allows the extra tests) therefore gives F: seed d_F=773w_F, threshold
k_F=12549w_F, error <=2^-E_F at requested exponent E_F, all seeds retained, for
arbitrary correlated (A_0..A_21) independent of (B_0..B_21,Y_0..Y_21).

**Lemma F.** At source length n, output z, advice 4 and exponent E_F=Ee, with
L=lg(n+1), the selected w_F exists and

    w_F <= fz*z+fE*Ee+fL*L+f0,                                          (LF)
    fz=16777216/8387115, fE=857069330628608/8387115, fL=58274116272128/8387115, f0=38930683166457856/2795705.

Proof. Each test reads R_i(w)<=w with R_i nondecreasing in w. Here e=Ee+10
(lg(26*21)=10), Lemma A's b=e+h_A+4=Ee+27, U=max(n,773w), lg(773w+1)<=10+lg w,
and lg(4U)<=12+L+lg w:

  R1=ceil(z/513)+Ee+15<=z+Ee+15 (first test of (A2));   R2=lg(773w+1)<=10+lg w;
  R3=G_U=16384(4b+S_U+5), S_U=pow2ceil(53rho_U)<106rho_U, rho_U<=Ee+2L+2lg w+53,
     so R3<16384(110Ee+212L+212lg w+5731);
  R4=128(Ee+18);
  R5=128G^ca with b'=Ee+33, rho'<=Ee+2lg w+59, S'<90rho', so
     R5<524288(94Ee+180lg w+5447).

Each R_i is at most a constant plus a constant times lg w, so large powers of
two pass every test and w_F exists. At w_F/2 some test fails, so w_F<2R_i(w_F)
for some i (trivially if w_F=1). As R2,R4<=R5 and all R_i>=0,

  w_F<2(R1+R3+R5)<2z+cE*Ee+cL*L+clw*lg w_F+c0,
  cE=2(1+16384*110+524288*94), cL=2*16384*212,
  clw=2(16384*212+524288*180)=195690496, c0=2(15+16384*5731+524288*5447).

For every integer w>=1, lg w<=w/2^40+41 (w<=2^41 gives lg w<=41; above it,
log2 w<=w/2^40+40). With q=clw/2^40<1/1000, w_F(1-q)<=2z+cE*Ee+cL*L+c0+41clw.
Dividing by 1-q gives (LF), each f being the reduced c/(1-q).

### 4.3 Parameter algorithm

Put Ee=E+4. Every operation is exact integer/rational arithmetic.

1. (d_F,k_F): F at source length n, output z, advice 4, exponent Ee (4.2).
2. P0=(h0+7)S0+18al0+12, al0=Ee+h0+4, rho0=al0+2lg(4n)+2,
   S0=pow2ceil((4h0+1)rho0): the additive seed term in bound (2) of the
   recursive primitive (linear-and-bilinear-primitives.md, first part) at input
   n, height h0, error 2^-Ee.
3. w=pow2ceil(max(d_F, 4Ee+2P0+a+64, 2^12)).
4. m1=w, r2=2w, m0=4w.
5. r=ceil((4/(3theta))(m0+2P0+4Ee+a+40)); s=floor(theta*r); b=10s; v=10r.
6. ell0=ceil(v/B0); K00 and d0 are the exact threshold K_0 (recurrence (1))
   and exact seed length s0+h0*S+W of the recursive primitive at input n,
   height h0, leaf ell0, error 2^-Ee.
7. e_L=d0+2Ee; r_L=pow2ceil(lg(r*m0+1)+m0+2e_L); d1=2r_L.
8. d2=2w+3Ee+d0+d1.
9. d_C=d0+d1+d2;  k_C=max(K00, k_F+2v+z+2Ee).

### 4.4 Evaluator

Slice y=y0y1y2 with lengths d0,d1,d2. Components:
* Ext0: the recursive primitive of step 6, output truncated to v bits; strong.
* C_adv: Section 3.1 with ell=10 blocks of r bits, capacity b, delta=2^-s, and
  the CZ/IP map of Section 3.4 at eta=1/10000000, advice zero-padded to b bits.
  It outputs 11 rows R_1..R_11 of r bits.
* LinExt: the AGHP extractor of weak-seed-hashing Section 1 from r bits to m0
  bits at error 2^-e_L, over the tower field of dimension r_L (seeded-linear
  Section 2); seed d1=2r_L bits.
* Ext2: the same extractor over the tower field of dimension r2=2w, from d2
  bits to m1=w bits at error 2^-Ee; seed 2r2=4w=m0 bits.
* F: Section 4.2; it reads the first d_F bits of its w-bit seed.

    V=Ext0(x,y0);  (R_1..R_11)=C_adv(V,alpha);
    W_j=LinExt(R_j,y1);  U_j=Ext2(y2,W_j);  Z_j=F(x,U_j,encode(j)), j=1..11;
    Q=Z_1+...+Z_11,                                                     (NBE)

with encode(j) the 4-bit binary label of j. The same map gives Q' on primed
inputs. The proof fixes a good row i inside each mixture component; the
evaluator XORs all rows with distinct labels and never needs i.

AGHP at a larger power-of-two dimension. Weak-seed Section 1 fixes
r=ceil(log(N+1))+m+2e and any presentation of F_(2^r) as an r-dimensional
binary space with the coordinate dot product. Its proof uses r only through
the bias bound (N-1)/2^r<2^(-m-2e), which only shrinks as r grows, and the
presentation only through bin(), so any F_2-linear coordinate bijection
works; the tower basis is one. So its (1) holds at every power of two at least
ceil(log(N+1))+m+2e: step 7 gives this for LinExt (N=r*m0), and (I3) gives it
for Ext2 (N=d2*m1).

### 4.5 Finite inequalities, each proved here

* (I1) s>=m0+2e_L. By (2) of the recursive primitive,
  d0<=6ell0+P0<=60r/B0+6+P0<=theta*r/8+6+P0, as B0>=480/theta. Step 5 gives
  (3/4)theta*r>=m0+2P0+4Ee+a+40, so m0+2e_L=m0+2d0+4Ee<=theta*r-a-28<s.
* (I2) Advice condenser. theta*r>=(4/3)m0>=28 and theta*r>=(16/3)Ee, so
  r>=ceil(28/theta)>=ceil(64/eta), and proof.md Section 4 applies verbatim:
  all guards of (8) and (11) hold and (4.1) gives e_A<=2^(-theta*r/4)<=2^-Ee.
  Capacity: s>=theta*r-1>=(4/3)(a+40)-1>=a, so b=10s>=a; also s>=27.
* (I3) 2w>=lg(d2*w+1)+w+2Ee. From w>=4Ee+2P0+a+64: theta*r<(4/3)(5w)+theta
  <=7w, so r<=7w/theta; d0<=7w/8+6+w/2<=2w; by (I1),
  d1<=4(lg(r*m0+1)+theta*r)<=28w+8lg w+4lg(28/theta)+4<=29w (lg(28/theta)=45,
  w>=2^12); so d2<=34w, lg(d2*w+1)<=2lg w+6<=w/2, and 2Ee<=w/2.
* (I4) d2-d0-d1-m1-Ee=m1+2Ee (definition of d2).
* (I5) k_C-2v-z-2Ee>=k_F (definition of k_C).
* (I6) (d0+Ee)-e_L=-Ee (definition of e_L).
* (I7) d_F<=w=m1 (step 3).
* (I8) B0*ell0>=v; and K00<=k_C<=Hinf(X)<=n.
* (I9) k_F<=k_C<=n, so F's guard holds; F sees 20<=t0 tamperings below.
* (I10) d_C=2d0+2d1+2w+3Ee<=65w<=20r=2v<=k_C, by (I3)'s bounds, Ee<=w/4 and
  r>=16w/(3theta).

### 4.6 Proof of (NB)

Left is (X,X') with every variable computed from it and from already revealed
right values; right is (Y,Y'). Initially they are independent. Every revelation
below is a function of one side and the history, or a kernel of such values, so
the conditional product persists [1.1]. Write Y=(Y0,Y1,Y2), Y'=(Y0',Y1',Y2').

1. Reveal T0=(Y0,Y0'), right-only. Given T0, V and V'=Ext0(X',Y0') are left
   functions. As X is independent of T0, D(V|T0)=E_(Y0) TV(Ext0(X,Y0),U_v)<=2^-Ee
   (strong primitive; K00<=Hinf(X) by (I8)). Couple V to Vbar, uniform given T0,
   by the common-mass kernel of (V,T0) [1.1]; the source law and V' are
   untouched. Cost eps0<=2^-Ee.
2. For each T0 apply Section 3.1 to (Vbar,V'): Vbar is uniform on 10r bits
   (entropy >=10r-b), V' is arbitrary, alpha!=alpha', all guards by (I2). Couple
   the matrix pair to the exact mixture (A,A') with label M by a kernel of
   (Vbar,V',T0), cost e_A<=2^-Ee; i=i(M) is the component's good row. Reveal M
   and A'_i (left kernels): H1=(T0,M,A'_i). On every branch Hinf(A_i|H1)>=s.
3. LinExt step. Y1 is independent of Y0 and Y0' has d0 bits, so
   Havg(Y1|T0)>=d1-d0 and Hinf(Y1|T0=t)>=d1-d0-Ee outside T0-mass 2^-Ee [1.1].
   Left and right are independent given T0, so Law(Y1|H1)=Law(Y1|T0) and A_i is
   independent of (Y1,Y1') given H1. On good T0, LinExt with a uniform seed has
   error 2^-e_L on A_i (entropy s>=m0+2e_L, (I1)); with seed Y1 of deficit
   d0+Ee, weak-seed (3) gives TV((W_i,Y1),(U_m0,Y1)|H1)<=2^(d0+Ee-e_L)=2^-Ee
   (I6). Reveal (Y1,Y1'), right-only: H2=(H1,Y1,Y1'). W'_i=LinExt(A'_i,Y1') is
   fixed. Given H2, W_i is a function of the left variable A_i, and its distance
   from uniform is f(Y1)=TV(LinExt(A_i|H1,Y1),U_m0), whose average over H2 is the
   H1-average above. So D(W_i|H2)<=2*2^-Ee (bad T0 charged 1). Couple W_i to
   Wbar_i, uniform given H2, by a kernel of (W_i,H2); Wbar_i is left-sided given
   H2. Cost delta_W<=2*2^-Ee.
4. Ext2 step. U'_i=Ext2(Y2',W'_i) is a right function of m1 bits; reveal it:
   H3=(H2,U'_i). Wbar_i stays uniform (right-only revelation). Leakage of Y2:
   for h=(y0,y0',mu,a',y1,y1',u'), as (M,A'_i) is independent of the right pair
   given T0,

     P(Y2=y2,H3=h) <= P(Y0=y0,Y1=y1,Y2=y2)P(M=mu,A'_i=a'|T0=(y0,y0'))
                   = 2^(-d_C) P(mu,a'|y0,y0').

   Summing the maximum over y2 across all h, the (mu,a') sum is 1 for each
   (y0,y0'), so the total is at most 2^(d0+d1)2^(d0+d1+m1)2^(-d_C)
   =2^(d0+d1+m1-d2). Hence Havg(Y2|H3)>=m1+3Ee (I4), and Hinf(Y2|H3=h)>=m1+2Ee
   outside mass 2^-Ee. Put U_i=Ext2(Y2,Wbar_i): the AGHP extractor (field by
   (I3)) with uniform seed Wbar_i independent of Y2 gives D(U_i|Wbar_i,H3)<=2^-Ee
   on good branches. Y2 is independent of every left variable given H3
   (X,X',V,V',Vbar,A,A',M,Wbar_i, coins), so Law(U_i|H3,all left)=
   Law(U_i|H3,Wbar_i). Couple U_i to Ubar_i, uniform given (H3, all left), by a
   kernel of (U_i,H3,Wbar_i); cost delta_U<=2*2^-Ee. This kernel reads only
   values present at step 4, so Ubar_i is defined before step 5, and given the
   left record it is right-sided.
5. Reveal V,V' and the left auxiliary OUTCOMES (Vbar,A,A',M,Wbar_i),
   not X,X'. No repair coin is revealed: in particular the Ext2 repair
   randomness of step4 remains on the right, unconditioned. Every W_j=LinExt(A_j,Y1) and W'_j=LinExt(A'_j,Y1')
   is now fixed, so every U_j=Ext2(Y2,W_j), U'_j=Ext2(Y2',W'_j) and Ubar_i is
   right-sided. Call the history G. No extra leakage: the joint density is
   P_left(x,x') 1{(V,V')=(Ext0(x,y0),Ext0(x',y0'))} P_right(y,y') K, where K is
   the product of the coupling-outcome kernels of steps1-3 and the revealed
   history indicators; step4's unrevealed repair kernel sums to1 in this
   marginal. Each factor reads only (v,v',t0), earlier auxiliaries and
   right values, never (x,x') directly. Fixing G's values and summing the right
   variables not in G leaves a factor free of (x,x'), so
   Law(X,X'|G)=Law(X,X'|T0,V,V'). X is independent of T0 and (V,V') has 2v bits,
   so Hinf(X|G)>=k_C-2v-Ee outside mass 2^-Ee. Step4 made Ubar_i uniform
   given H3 and all left variables; thus it remains uniform given G.
6. Reveal Z'_i=F(X',U'_i,encode(i)): U'_i is fixed in G, so this is a z-bit left
   revelation. With H5=(G,Z'_i), Hinf(X|H5)>=k_C-2v-z-2Ee>=k_F (I5) outside total
   mass 2*2^-Ee, and Ubar_i stays uniform.
7. F step. Given H5 the source tuple (X for the honest call and the other 10
   honest rows; X' for the 10 tampered rows j!=i) is independent of the seed
   tuple (Ubar_i; U_j,U'_j for j!=i). The honest seed's first d_F bits are
   uniform (I7); encode(i) differs from every other label; the 20 tamperings
   are at most t0=21 (duplicate one); the tampered row i is fixed in H5 and is
   not a tampering. F's contract (wrappers (4) via Lemma A, additive parts zero)
   makes Zbar_i=F(X,Ubar_i,encode(i)) 2^-Ee-close to uniform jointly with every
   other row output and every seed. With Qbar=Zbar_i+sum_(j!=i)Z_j and
   Q'=sum_j Z'_j (Z'_i fixed), D(Qbar|Q',seeds,H5)<=2^-Ee.
8. Retain Y,Y'. Given (H5, seed tuple), the outputs depend only on (X,X') and the
   seeds. Given H5, the left source pair is independent of the right pair
   and the unrevealed Ext2 repair randomness: its kernel is right-sided
   once Wbar_i is fixed. All seeds are functions of that right record.
   Hence Law(Y,Y'|H5,seeds,outputs)=Law(Y,Y'|H5,seeds) is a
   common kernel of the real and ideal laws, and appending (Y,Y'), including
   unused seed suffixes, costs nothing [1.1]. Averaging over H5 and charging the
   two X tails gives synthetic error <=3*2^-Ee.
9. Restore. Real and synthetic (Q,Q',Y,Y') agree unless V!=Vbar, the matrix pair
   differs from (A,A'), W_i!=Wbar_i or U_i!=Ubar_i, so
   zeta<=eps0+e_A+delta_W+delta_U<=6*2^-Ee. By [1.1] (D_P<=D_Q+2zeta),
   D(Q|Q',Y,Y')<=3*2^-Ee+12*2^-Ee=15*2^(-E-4)<2^-E.

### 4.7 Constants

Let L=lg(n+1) and s'=2(4h0+1)=394. Each line substitutes the previous ones, using
pow2ceil(x)<2x, S<2(4h+1)rho for a primitive's power-of-two S, lg(4n)<=L+2,
max<=sum for nonnegative terms, and Ee=E+4. All fractions are exact; decimals
in brackets are upper bounds for reading only.

* (P) P0<pE*Ee+pL*L+p0: pE=(h0+7)s'+18=22082, pL=2(h0+7)s'=44128,
  p0=(h0+7)s'(h0+10)+18(h0+4)+12=1302742 (rho0<=Ee+2L+h0+10).
* (w) w<2(773w_F+4Ee+2P0+a+4160)<wz*z+wE*Ee+wL*L+2a+w0 by(LF):
  wz=2*773fz=25937575936/8387115, wE=2(773fE+4+2pE)=1325029926036018608/8387115,
  wL=2(773fL+2pL)=90093264183152768/8387115, w0=2(773f0+2p0+4160)=60186850766933403416/2795705.
* (r) With k'=4/(3theta)=(4/3)10^12,
  r<=k'(4w+2P0+4Ee+a+40)+1<rz*z+rE*Ee+rL*L+ra*a+r0, where
  rz=4k'wz=83000242995200000000000/5032269, rE=k'(4wE+2pE+4)=4240096059668935801600000000000/5032269,
  rL=k'(4wL+2pL)=288299037556666009600000000000/5032269, ra=9k'=12000000000000, r0=k'(4w0+2p0+40)+1=192597928281608070467200001677423/1677423.
* (K) Bound(3), ell0<v/B0+1 and al0=Ee+h0+4 give K00<2v+xE*Ee+xL*L+x0,
  xE=2B0(4+s')=448108162923364352, xL=4B0s'=887209126591987712,
  x0=2B0(5+4(h0+4)+s'(h0+10))=26416989514248486912.
* (kC) k_C<=k_F+2v+z+2Ee+max(0,K00-2v), k_F=12549w_F,v=10r, so
  k_C<kz*z+kE*Ee+kL*L+ka*a+k0, with
  kz=12549fz+20rz+1=754547663592784694273827/2287395, kE=12549fE+20rE+2+xE=38546328840175448261824476991046/2287395,
  kL=12549fL+20rL+xL=2620902370822156012302941290496/2287395, ka=20ra=240000000000000, k0=12549f0+20r0+x0+8=1750890277247691061813215107334812/762465.
* (c) alpha=max(kz,kE,kL,ka)=kE, beta=k0+4kE. Hence
  k_C<=alpha(z+a+E+L)+beta<=(alpha+beta/5)(z+a+E+L+1), since z+a+E+L+1>=5.
  With(I10), c_std=ceil(alpha+beta/5)=489603919856837338526670350.
* (dC) The generic upper bound remains d_C<=65w<130(773w_F+4Ee+2P0+a+4160).
  Concrete outer calibration uses the exact seed layout in proof.md Lemma B.

Running time. At fixed h0 every component is polynomial in its dimensions: the
recursive primitive (its O(n^8) bound), the tower fields of dimensions r_L and
2w (O(dim^2) per product), the CZ/IP setup in 2^O(r+b) bit operations (prime
fields Section 2), and F (Lemma A's evaluator around the fixed-advice
constant-advice breaker). Every NB dimension is linear in n's logarithm and in
z,a,E with fixed coefficients, so at the outer call, where they are O(log of
the outer length), each NB call is polynomial in the outer length.

Attribution: the advice condenser and row-XOR pattern follow Li (Section 3 and
the header above); AGHP is credited in weak-seed-hashing. Using the advice row
as the source of a uniformly seeded linear extractor, and the proof above, are
written for this package.

## 5. Scope

The dependency chain is acyclic: seeded primitives -> old constant-advice
breaker -> NB's inner affine breaker F (proof.md Lemma A at degree 21); CZ+IP
-> advice condenser; advice condenser, AGHP extraction and F -> standard
breaker NB, available for the later outer affine call (proof.md Lemmas A, B).
Never use NB to instantiate its own inner breaker. Evaluators never compute
proof labels or partitions. Imports with 2^O(r) setup are polynomial in outer
N only after a fixed r=O(log N) bound. Global assembly, exceptional-length
fallback, outer sampler/family/circuit proof and package fit are outside this
module.
