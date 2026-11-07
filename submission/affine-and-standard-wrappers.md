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
are credited in their modules.

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
  Havg(U|U_J,H)>=K+|I\J|b+q. (1)

With V=F(U,S), V_j=F(U_j,S_j),

  D(V|V_(I union J),S,S_1,...,S_t,H)<=eta+2eps. (2)

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

  C: {0,1}^d x {0,1}^d0 x {0,1}^a -> {0,1}^dx;
  L0: {0,1}^n x {0,1}^s0 -> {0,1}^d0, threshold K0;
  Lr: {0,1}^n x {0,1}^dx -> {0,1}^r, threshold Kr;
  E: {0,1}^d x {0,1}^dy -> {0,1}^dx, threshold KE;
  Lm: {0,1}^r x {0,1}^dx -> {0,1}^dy, threshold Km.

Every L is strong and linear in its source for each fixed seed; E is strong.
C is a strong ordinary one-tampering breaker of threshold KC: arbitrary
correlated source pair independent of arbitrary correlated seed pair, honest
seed uniform, distinct advice, both complete seeds retained. Its first input
really has length d. Padded seeds may ignore suffixes while retaining them.
All actual primitive guards must hold. Put

  TR0=(t+1)(s0+d0+dx), TL0=(t+1)d0, T0=TR0+TL0,
  TR=TR0+H(t+1)(2dy+2dx)-(t+1)dx, TL=TL0+H(t+1)2dy, T=TR+TL+(t+1)dx.

Require positive legal dimensions, k<=n, and

  s0<=d, dy<=r, m<=r,
  r>=(t+1)dy+Km+q,
  k>=max(K0,TL+Kr+t*r+q),
  d>=(t+1)(s0+d0)+KC+q,
  d>=TR+KE+t*dx+q. (3)

The evaluator on (x,y,alpha) is

  s=Prefix_s0(y); r1=L0(x,s); s2=C(y,r1,alpha); w0=Lr(x,s2).
  For i=1,...,H:
    p=Prefix_dy(w_(i-1)); qm=E(y,p); v=Lm(w_(i-1),qm);
    qr=E(y,v); w_i=Lr(x,qr).
  Return Prefix_m(w_H).

The outer A# of proof.md outer Section3 has m<=r, so no long-output
extension is used here. Section4.2 constructs F# separately from restricted8.2.
Let X_j=A_j+B_j, j=0,...,t. Suppose the arbitrary correlated A tuple is
independent of the arbitrary correlated (B,Y) tuple, Hinf(A_0)>=k, Y_0 is
uniform d-bit, and alpha_0 differs from every alpha_j. Other components need
no entropy, uniformity, or functional relationship to the honest components.
For F the displayed evaluator,

  D(F(X_0,Y_0,alpha_0) |
    (F(X_j,Y_j,alpha_j))_(j=1)^t,Y_0,...,Y_t)
      <=delta_H=(13*2^H-8)eps<=26t eps. (4)

Distance can be capped at 1. We retain the original seed without adding a
synthetic tampering.

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

  D(WiA_0|WiA_J,Zi)<=delta_i, delta_0=5eps. (5)

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

  Havg(W*|WA_J2,H1)>=r-(t+1-|J2|)dy. (6)

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
The new messages (PB,P,Qm,VB,V,Qr) total (t+1)(4dy+2dx) bits. Every E call is covered by d-TR>=KE+t*dx+q; every Lr call by
k-TL>=Kr+t*r+q.
The Lm entropy came separately from (6), not from the original A entropy.

This proves delta_i=2delta_(i-1)+8eps, hence delta_i=(13*2^i-8)eps,
including non-power-of-two t. At i=H, append the entire right tuple to the
bound against all tampered left outputs using the conditional product given
ZH. Translate by their right contributions, discard excess history, and take
prefixes. This proves (4), retaining every Y_j.

Causal same-side ledger. On a rectangular history a message from the other
side leaves the conditional source law unchanged. An adaptive b-bit own-side
message multiplies average maximum atom by at most 2^b: for each old history
sum maximum joint atoms over its at most 2^b messages, then average. Induction
thus charges only own-side bits, even for history-dependent functions.
Initialization sends S,R1B,S2 right (TR0 bits) and R1 left (TL0 bits).
Each round sends PB,Qm,VB,Qr right ((t+1)(2dy+2dx) bits), and P,V left
(2(t+1)dy bits). Before the second E call of round i, that round's Qr has
NOT been revealed: right cost is TR0+i(t+1)(2dy+2dx)-(t+1)dx<=TR.
The first E call precedes Qm and VB as well, so costs less. If H=0 there are no E calls, so the E guard d-TR>=KE+t*dx+q is vacuous and TR's -(t+1)dx term is used only for H>=1. C is called before R1,S2 are revealed, at right cost
(t+1)(s0+d0). Every Lr call has left cost <=TL. These are exactly the
three source constraints (3), with q for the entropy tail. Lm still uses
(6). The full transcript length is T: Qr is subsequently revealed in full,
so rectangularity, (5), the error recurrence and retention of ALL original
seeds are unchanged. The outer assembly uses this ledger; F# uses
restricted8.2's ledger and CA# its separate eight-label ledger.

### 2.4 Actual primitive use

Selected linear primitives use PXT/PXG7g with their actual recurrence and
explicit integer guards. The outer assembly uses Section2.1 withC=NB.
The selected innerF# uses the five-round antichain theorem8.2 with eight
labels and14tamperings,instantiated in4.2. Its zero-shift independent-
side joint contract is exactly the one used by NB4.6 step7. The full
conditional-average lift remains sampler6.1 and all original source/seed
quantifiers and retained padding suffixes remain explicit.

## 3. Finite advice-condenser reduction

### 3.1 Contract

Fix ell>=2,r,s>=1, b=ell*s, n=ell*r, 1<=a<=b, 2^s>=ell, n>=b.
Supply N on ell independent (r+a)-bit inputs, output h<=r bits, such that
inputs each of entropy >=r-2b and coordinate maps fixed-point-free on their
honest supports satisfy

  TV((N(U),N(f(U))),(U_h,N(f(U))))<=eN. (7)

Put B=(ell+9)b, choose delta in (0,1), and require

  h-B-log(1/delta)>=s. (8)

C_adv(x,alpha), x=(x1,...,xell), outputs the ell raw r-bit blocks followed
by N(x1||alpha,...,xell||alpha), zero-padded to r bits. For arbitrary
correlated X,X' with Hinf(X)>=n-b and alpha!=alpha', its matrix pair is within

  eA=(ell+2)2^-s+2^-b+delta+2^B eN (9)

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

### 3.3 Conditioning near-independent output

If (W,T) is e-close to (U_h,T) and event E has probability p>0, there is a
law (W*,T*) within delta+e/p of (W,T)|E, with the same T|E marginal and

  Hinf(W*|T*=t)>=h-log(1/p)-log(1/delta). (10)

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

### 3.4 A balanced seven-input weak-AGHP advice map

Fix eta=268435456/1718070947847,integers r,a>=1,D>=0,u=r+a,
h=floor((8/7)eta*u)-8. For independenthonest rawr-bit U1,...,U7,
eachofentropy>=r-D,andfixed coordinatewise partners,appenddistinct
fixeda-bit advices to the first6 pairs. Assume
 (a+D)/u<=eta/2,r-D>=(1-eta)(u+1)+1,
 eta*u>=64*2^20,u>=3200000,64D<=eta*u,1<=h<=r. (11)
N_D is the following actual fixed-parameter map. After checkingguards,
choose the canonical primeabove2^u;compute E6 on the first6 embedded
inputs andencode W in nW=2(u+1)bits. Set
 eL=ceil(D/2)+ceil(eta*u)+2,R_L=lg(nW*h+1)+h+2eL,d_L=2R_L.
Build the canonical degree-R_L binary field andapply the supplied AGHP
matrixextractor to W,outputh,errorexponent eL,using firstd_Lbits of
T=U7. IgnorebutRETAINtheotherrawseedbits. Theprimedmap uses SAME
parameters withWprime,Tprime. A failedliteralguard returnszero.

Its fieldseed fits: h<=u,ceil(D/2)<=D/2+1/2,
lg(nW*h+1)<=1+2lg(u+1),andceil(eta*u)<=eta*u+1 give
 d_L<=4lg(u+1)+(44/7)eta*u+2D+12
     <=4lg(u+1)+(1415/224)eta*u+12<=7eta*u+4lg(u+1)+12.
Foru in[2^k,2^(k+1)),k>=20,4lg(u+1)+12<=4(k+1)+12
<=2^k/4096<=u/4096,by96<=256 at20 anddoubling.
Thusd_L<=(7/6000+1/4096)u<(1-eta/2)u<=r.
h>=1 followsfrometa*u>=64;h<=r alsofollowsfrometa<=1/6000
and a<=eta*u/2. Theliteralchecks make thealgorithm total.

Distinctappendedadvices makeallfirst6 mapsfixed-point-free,so restricted
E6 appliesonitsactualfirst-fourdomain. Its source threshold is supplied
by(11). The directlyproved E_p(256/20481) theorem appliesatu>=3200000;
(80+7*2^-20)eta=256/20481,eta*logp>=64*2^20. The first-six joint output(W,Wprime)
is zeta<=45p^-eta from(W*,Wprime*) with conditional honestentropy
 kappa=(19/6)eta logp>=(19/6)eta*u.
The pair(T,Tprime) is independentof this entire first-six pair,andremains
so during smoothing. OnEVERY Wprime* fiber,the exactAGHP source guard is
 h+2eL<=(22/7)eta*u+D-1<=1415/448 eta*u-1<(19/6)eta*u.
Here 19/6-1415/448=11/1344>0. LemmaW2's actual squared-collision
calculation,with FULLr-bit seeddeficitD andignoredsuffixretained,gives
error<=2^((D-1)/2)2^-eL<=2^(-ceil(eta*u)-5/2)<=(1/4)2^(-eta*u).
AverageoverWprime*. AppendTprime throughLaw(Tprime|T),independentof
thesmoothedpair;theprimedAGHPoutput isdeterminedbyWprime*,Tprime.
Restoring the real jointlaw andits ideal tamperedoutputtarget costs2zeta.
Consequently
 TV((N_D(U),N_D(Uprime)),(U_h,N_D(Uprime)))<=91*2^(-eta*u). (12)
BothCOMPLETErawseventhcoordinates are retainable. The seventhmapneednotbefixed-point-free.

This is a theoremonlyfor actualfixedappended-advice/rawr-bit inputs.
No earlier general-source template isinstantiatedoutside thatdomain;
inparticularrawseedentropy isnot hiddeninanunknownsuffix. Forselected
actualNallr,u,R_LareO(logN),sothecanonicalfield/primealgorithms cost
2^O(logN);largerrelativeextensionsusetowers.

### 3.5 Uniform-source advice lemma U

Fix ell=7, integers r,a,s,sig>=1 with sig>=4, a>=1, s>=7, r>=s+sig.
Let X be exactly uniform on ell*r bits, X' arbitrarily correlated, and
alpha!=alpha' be a-bit advice strings. Let eta=268435456/1718070947847,u=r+a, and assume
 (a+sig)/u<=eta/2, r-sig>=(1-eta)(u+1)+1,
 eta*u>=64*2^20,u>=3200000,64sig<=eta*u,1<=h=floor((8/7)eta*u)-8<=r,
 eta*u>=7s+17sig+8. (U0)
Use Section3.4's same weak-AGHP algorithm N_D with D=sig at these new input parameters.
Its proof applies with deficitD=sig and rho=(a+sig)/u: the first-six
source threshold follows from r-sig>=(1-eta)(u+1)+1. The independent
seventh raw seed has entropy>=r-sig. Its exact weak-AGHP guards and
whole-seed/partner retention were proved in3.4. Thus eN<=91*2^(-eta*u),
with outputh=floor((8/7)eta*u)-8. All prime/field/seed guards follow from(U0).
Return the ell raw blocks and N, zero-padded to r bits.
The matrix pair is within14*2^-sig of a convex mixture with a fixed row
of conditional min-entropy at least s in every component.

Proof. Realize X' by a convex mixture of fixed functions f(X), sampling
each value from its conditional law; every component keeps X uniform.
Call (z,z') heavy in coordinate i when its pair mass is at least
2^(-r-s-2sig). Partition points with a nonheavy pair by the first such
coordinate into Ti. Discard pieces of mass<2^-sig, cost at most7*2^-sig.
In a remaining Ti every joint atom is<2^(-r-s-sig); summing over at most
2^r partner values gives average conditional maximum atom<2^(-s-sig).
Thus outside conditional partner mass2^-sig the honest row has entropy
at least s. Split by partner values. Summed failures cost at most2^-sig.
If the all-heavy SH has mass<2^-sig discard it at that cost; otherwise
|SH|>=2^(7r-sig). For every honest coordinate z its original marginal
mass is exactly2^-r. Hence its heavy-partner set Wi(z) has size at most
2^(s+2sig), for EVERY z. No exceptional heavy fiber needs removal.
Let Si be the coordinate supports of SH and Pset their product. Then
|Si|>=2^(r-sig) and |SH|/|Pset|>=2^-sig. Choose gi(z) uniformly from Wi(z)
independently for all coordinates and z. Each point of SH matches g in
all coordinates with probability at least2^(-7s-14sig). While a remaining
set has more than2^-sig|SH| points, averaging supplies a matched subset V
of size greater than2^(-7s-16sig)|Pset|. Remove it and repeat, leaving
total discarded fraction at most2^-sig of SH. This finite partition is
only an analysis object. Put B=7s+16sig; each V has density>2^-B.
On independent uniforms Si, append alpha and use coordinatewise gi with
alpha'. Distinct suffixes make all coordinate maps fixed-point-free;
Section3.4 and(U0) apply with deficit(a+sig)/u. Section3.3 conditioning
with delta=2^-sig and density>2^-B gives conditional output entropy at
least h-B-sig>=s, at TV cost at most2^-sig+2^B eN.
Keep the conditioned partner marginal and append all raw rows through
their original conditional kernel given the final-row pair. This preserves
the selected-row entropy property and distance. Replace discarded pieces
by arbitrary laws with a uniform honest selected row. WhenSH isdiscarded,the old12*2^-sigma majorant remainsvalid. Otherwise
losses7+1+1+1=10 times2^-sig plus2^B eN suffice. New(U0) gives
2^B eN<=91/256*2^-sig<2^-sig,independentlyofh.
It ALSOproves the entropyguardh>=8s+17sig: h>=8eta*u/7-9,
whereas8s+17sig<=8eta*u/7-17sig/7-64/7. Theirgap is at least
17sig/7+1/7>0. Thus h-B-sig>=s as usedabove. The conservative
14*2^-sig error follows.
All source/partner changes are joint laws and no conditional independence
between honest and tampered raw blocks is assumed.

## 4. The standard breaker NB

This section supplies the standard breaker of proof.md Section 3.1. It needs
a general source X and a UNIFORM seed Y only. The rate-theta advice row is used
as the SOURCE of the AGHP extractor of weak-seed-hashing Section 1, seeded by a
fresh slice of Y, so no rate-to-0.8 condensation and no ordinary condenser is
used. Inside this section F names NB's inner breaker, lg(x) is the least j>=0
with 2^j>=x, pow2ceil(x)=2^lg(x)<2x, and [1.1] cites Section 1.1.

### 4.1 Statement

Fix eta=268435456/1718070947847,theta=eta/7,t0=14,a0=4. Use
Ext0=PXT(vNB0) andExt2=PXT(vNB2),
oflinear-and-bilinear-primitives.md Section7g. Ext0 has
 c0=7885491/14859351621632,pi0=686887390833373758731/1160886845440000000,
 7c0=7885491/2122764517376<theta=268435456/12026496634929.
Let F be the inner breaker F# of Section 4.2. For integers n,z,a,E>=1, Section
4.3 computes a seed length d_C and threshold k_C, and Section 4.4 defines a
deterministic NB: {0,1}^n x {0,1}^(d_C) x {0,1}^a -> {0,1}^z. If k_C<=n and
d2<=n, then for every pair (X,X') on n bits independent of every pair (Y,Y')
on d_C bits, with Hinf(X)>=k_C, Y uniform, copies arbitrarily correlated
within each pair, and alpha!=alpha',

    D(Q|Q',Y,Y') <= 15*2^(-E-4) < 2^-E, Q=NB(X,Y,alpha), Q'=NB(X',Y',alpha'). (NB)

So NB is a standard breaker in the sense of proof.md Section 3.1 at every
a,z,E>=1, with the guards k_C<=n and d2<=n, retaining both complete seeds. If a
guard fails the map outputs 0^z and no claim is made (the totality convention
of proof.md Section 3).

### 4.2 Five-round antichain inner breaker F#

For source lengthn,requested outputz>=1,advice among NB labels1,...,8,
and exponentE_F>=1, compute the public INTEGER
 w=max(2^20,E_F+2lg(4n)+32,ceil(z/128)),e=E_F+7,q=e,t=14,c=15.
Every actual slot is ceil(A*w),with these fixed rational coefficients:
 a=ceil((416707761/16777216)w).
 s=ceil((197434393/8388608)w).
 s0=ceil((416707207/16777216)w).
 hc=ceil((1427000461/4194304)w).
 d=ceil((39367206957/16777216)w).
 k=ceil((135983578456005/1048576)w).
Write h=hc. Use these exact Section7g PXG maps,all at exponent e:
 I=PXG(vMI)(inputn,outputh,seeds0),
 F=PXG(vMF)(inputn,outputh,seeda),
 G=PXG(vMG)(inputd,outputa,seeds),
 J=PXG(vMJ)(inputh,outputs,seeda).
Pad their seeds to the integer slots and retain suffixes. Ifk>n or the
label is outside1,...,8,return0^z. Otherwise use the ACTUAL five-round
algorithm of restricted8.2: startW=I(X,prefix_s0 Y);at each of the five
coordinates computeT=prefix_s W,A=G(Y,T),Z=J(W,A),B=G(Y,Z),chooseC=B
for codebit1 andC=A for codebit0,andsetW=F(X,C). Returnprefix_z W.
The fixed label codes are11000,10100,10010,01100,01010,00110,10001,01001.
For NB's4-bit encode(j),j1,...,8,the corresponding code is the listedjth
word.
Putd_F=d,k_F=k. Restricted8.2 proves this actual algorithm for an arbitrary
source15-tuple independent of an arbitrary seed15-tuple, honestentropy>=k,
honestseed uniform,and honestlabel different from all others. Correlations inside each tuple are arbitrary; each nonhonest label occurs at most twice. NB has exactly two copies of each other label on every conditioned branch. Error is
76*2^(-E_F-7)<2^-E_F jointly with ALL14 other outputs and ALL15 complete
original seeds. NB4.6step7 uses this zero-shift contract.

Put delta=2^-20. At integerw>=2^20,ceil(Aw)<=(A+delta)w. Every printed
RHS below equals W+sum(q_i*A_i)+delta*sum(q_i) for the indicated integer
contract;the left slot is at least its coefficient timesw.
The exact primitive coefficients and all twelve requirements are:
 I: c=1/64,pi=1599217176853/81920000000,lam=1568143971086491170292166323152135869981653/9070022988439151740886493662857981653,g=1124756448423955841592228055333615006074427855275773/15872540229768515546551363910001467892750000
 F: c=1/64,pi=1599219880213/81920000000,lam=1213960815548266455022223884715200563987507451/7022328621044905155799428172344991757451,g=6506464572399320001629692138918186988605209404561451967/91830451198279528960454060715280661443590000000
 G: c=1/2,pi=227679317709/20480000000,lam=327106124366733/121254709366733,g=104448828912292844066081/2425094187334660000000
 J: c=1/2,pi=267667674829/20480000000,lam=14667502769683/7526814769683,g=9799323330499711966857/284346335743580000000
 aJ: [32555293767181/1310720000000]w <= (416707761/16777216)w.
 aF: [520884700662653/20971520000000]w <= (416707761/16777216)w.
 sG: [61698247744877/2621440000000]w <= (197434393/8388608)w.
 s0I: [520884008602493/20971520000000]w <= (416707207/16777216)w.
 hJ: [18720580962946050035828441/55024374387671040000000]w <= (1427000461/4194304)w.
 hcap: [128]w <= (1427000461/4194304)w.
 h>=s: [197434401/8388608]w <= (1427000461/4194304)w.
 dG: [372926502978212403736162262541/158930972661164277760000000]w <= (39367206957/16777216)w.
 d>=s0: [416707223/16777216]w <= (39367206957/16777216)w.
 kF: [10146035772289684767844044859350022129429056129772894664186681/78236605846495813729380125199476234568772485120000000]w <= (135983578456005/1048576)w.
 kI: [539601204907705514303621917321251220050015262682170284187/4160891185992437739435160740823424799277056000000]w <= (135983578456005/1048576)w.
 k>=d: [39367206973/16777216]w <= (135983578456005/1048576)w.
The rows prove bothaseeds,s,s0,h>=(11+Lam_J)s+(g_J+1)w,h>=128w>=z,
h>=s,d>=15s0+(75+Lam_G)a+(g_G+1)w,d>=s0,
k>=(44+Lam_F)h+45s+(g_F+1)w,k>=K_I andk>=d.
These are exactly selective8.2(1) at2c+t=44,3c=45,5c=75 andnu<=11,includingq<=w. They implyKJ<=h,
KG<=d,KI,KF<=k;K_I>=h,so all primitive inputs are<=n whenk<=n.
Selectedheights<=8 ande=E_F+7 give
 e+h_primitive+6+2lg(4u)<=E_F+21+2lg(4n)<=w.
The PXT/PXG internal source,seed,ambient and field guards are proved in7g.
The PUBLIC output cap inw defines the map for everyz. Degrees are logarithmic;larger relative extensions use towers,as7a/7g.

### 4.3 Parameter algorithm

Put Ee=E+4 and the unit w_N=max(2^20,Ee+38+6+2lg(4n)).

1. (d_F,k_F): F# at source length n, output z, advice 4, exponent Ee (4.2);
   m1:=d_F.
2. m0,K2: the exact seed length and threshold of Ext2 at input length n,
   output m1, exponent Ee.
3. sig=Ee+4; r=ceil((m0+pi0*w_N+(38/7)Ee+a+53)/(theta-7c0));
   s=floor((eta*r-17sig-8)/7); b=a;
   v=7r.
4. d0,K00: the exact seed length and threshold of Ext0 at input n, output v,
   exponent Ee.
5. e_L=ceil((d0+Ee-1)/2)+Ee; r_L=lg(r*m0+1)+m0+2e_L; d1=2r_L.
6. d2=d0+d1+m1+K2+Ee.
7. d_C=d0+d1+d2; k_C=max(K00, k_F+2v+z+2Ee).

The guards are k_C<=n and d2<=n (Section 4.1).

### 4.4 Evaluator

Slice y=y0y1y2 with lengths d0,d1,d2. Components:
* Ext0: the strong linear primitive of step 4, output truncated to v bits.
* C_adv: LemmaU with ell=7 blocks of r bits, accuracy sig, and
  the weak-AGHP map of Section3.4 at eta=268435456/1718070947847 and uniform-source LemmaU
  of Section3.5, sig=Ee+4, actual advice length b=a.
  It outputs 8 rows R_1..R_8 of r bits.
* LinExt: the AGHP extractor of weak-seed-hashing Section 1 from r bits to m0
  bits at error 2^-e_L, over the canonical degree-r_L binary field (weak-seed
  Section 4); seed d1=2r_L bits.
* Ext2: the strong primitive of step 2 applied to y2 zero-padded to n bits,
  with seed W_j (m0 bits) and output m1=d_F bits.
* F: F# of Section 4.2; it reads all of U_j.

    V=Ext0(x,y0); (R_1..R_8)=C_adv(V,alpha);
    W_j=LinExt(R_j,y1); U_j=Ext2(y2,W_j); Z_j=F(x,U_j,encode(j)), j=1..8;
    Q=Z_1+...+Z_8, (NBE)

with encode(j) the 4-bit binary label of j. The same map gives Q' on primed
inputs. The proof fixes a good row i inside each mixture component; the
evaluator XORs all rows with distinct labels and never needs i.

AGHP at its exact integer degree. Weak-seed Section1 defines
r=lg(N+1)+m+2e over any r-dimensional binary field presentation,
with its coordinate dot product. Step5 uses exactly this degree at
N=r*m0,m=m0,e=e_L. Its root-counting and strongness proof therefore
applies directly, as does the squared-collision weak-seed LemmaW2.
Choose the first irreducible degree-r_L binary polynomial as in weak-seed
Section4; at every selected actual outer use r_L=O(logactualN), so its
2^r_L*poly(r_L) setup is polynomial in actualN. Source-linearity and full seed
retention are unchanged, and d2 is computed from this actual d1.

### 4.5 Finite inequalities, each proved here

* (I1) s>=m0+2e_L. PXT gives d0<=7c0*r+pi0*w_N and
  2e_L<=d0+3Ee. The definition of r and floor give
  s>=m0+7c0*r+pi0*w_N+3Ee+a+288/7,
  hence s>=m0+2e_L and s>=a+40. The denominator theta-7c0 is positive.
* (I2) Lemma U. Since eta*r>=7m0+7pi0*w_N+38Ee+7a+371,
  eta*r>=2(a+sig)+2 and eta*r>=64, for sig=Ee+4. Thus
  (a+sig)/(r+a)<=eta/2 and eta*(r+a+1)>=a+sig+2,
  which is exactly r-sig>=(1-eta)(r+a+1)+1. Also a<eta*r,
  eta<=1/6000, and sig<eta*r/2; (8/7)eta*(1+eta)<1,so h=floor((8/7)eta*(r+a))-8<=r
  and r>=s+sig. The floor definition gives eta*r>=7s+17sig+8;thereforeh>=8s+17sig bytheproofin(U0). Also pi0>=16 bytheprintedconstant,andw_N>=Ee+44,so eta*r>=7pi0*w_N>=112(Ee+44)>64sig;the newdeficitguard holds.
  In particular s>=a+40>=41. Also u=r+a>=7pi0*w_N/eta>3200000,as in(I3),so its explicit ring guard holds. Lemma U gives e_A<=14*2^-sig<2^-Ee.
* (I3) Ext2 is legal: K2<=d2<=n; the conservative unit w_N covers
  actual heights21 and8 of both exact PXT primitives. The prime guard
  u=r+a>=r>=7pi0*w_N/eta>=7*2^20/eta>3200000, sincepi0>=1,w_N>=2^20; this proves the independent prime guard rather than inferring it from64/eta. Also eta*r>=7*pi0*w_N>=112*2^20>64*2^20,usingpi0>=16,w_N>=2^20,so the strengthened CZ constant-absorption guard holds at EVERY actual call.
* (I4) d2-d0-d1-m1-Ee=K2 (step 6).
* (I5) k_C-2v-z-2Ee>=k_F (definition of k_C).
* (I6) ((d0+Ee-1)/2)-e_L<=-Ee by definition of e_L; Lemma W2 applies.
* (I7) d_F=m1 (step 1).
* (I8) Ext0's output is at least v (its leaf); and K00<=k_C<=Hinf(X)<=n.
* (I9) k_F<=k_C<=n, so F's guard holds; F sees 14<=t0 tamperings below.

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
2. For each T0 apply Lemma U to(Vbar,V'): Vbar is exactly uniform on7r bits,
   V' is arbitrary, alpha!=alpha', and every finite guard holds by(I2). Couple
   the matrix pair to the exact mixture (A,A') with label M by a kernel of
   (Vbar,V',T0), cost e_A<=2^-Ee; i=i(M) is the component's good row. Reveal M
   and A'_i (left kernels): H1=(T0,M,A'_i). On every branch Hinf(A_i|H1)>=s.
3. LinExt step. Y1 is independent of Y0 and Y0' has d0 bits, so
   Havg(Y1|T0)>=d1-d0 and Hinf(Y1|T0=t)>=d1-d0-Ee outside T0-mass 2^-Ee [1.1].
   Left and right are independent given T0, so Law(Y1|H1)=Law(Y1|T0) and A_i is
   independent of (Y1,Y1') given H1. On good T0, LinExt with a uniform seed has
   error 2^-e_L on A_i (entropy s>=m0+2e_L, (I1)); with seed Y1 of deficit
   d0+Ee, weak-seed LemmaW2 gives TV((W_i,Y1),(U_m0,Y1)|H1)<=2^((d0+Ee-1)/2-e_L)<=2^-Ee
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
   =2^(d0+d1+m1-d2). Hence Havg(Y2|H3)>=K2+Ee (I4), and Hinf(Y2|H3=h)>=K2
   outside mass 2^-Ee. Put U_i=Ext2(Y2,Wbar_i): Ext2 is a strong (K2,2^-Ee)
   extractor (zero-padding Y2 to n bits keeps its min-entropy), and its seed
   Wbar_i is uniform and independent of Y2 given H3, so D(U_i|Wbar_i,H3)<=2^-Ee
   on good branches. Y2 is independent of every left variable given H3
   (X,X',V,V',Vbar,A,A',M,Wbar_i, coins), so Law(U_i|H3,all left)=
   Law(U_i|H3,Wbar_i). Couple U_i to Ubar_i, uniform given (H3, all left), by a
   kernel of (U_i,H3,Wbar_i); cost delta_U<=2*2^-Ee. This kernel reads only
   values present at step 4, so Ubar_i is defined before step 5, and given the
   left record it is right-sided.
5. Reveal V,V' and the whole LEFT auxiliary record (Vbar,A,A',M,Wbar_i,
   only the left-generated coupling coins of steps1-3), not X,X'.
   The step4 right-repair coins are NOT revealed and do not belong to G. Every W_j=LinExt(A_j,Y1) and W'_j=LinExt(A'_j,Y1')
   is now fixed, so every U_j=Ext2(Y2,W_j), U'_j=Ext2(Y2',W'_j) and Ubar_i is
   right-sided. Call the history G. No extra leakage: the joint density is
   P_left(x,x') 1{(V,V')=(Ext0(x,y0),Ext0(x',y0'))} P_right(y,y') K, where K is
   the product of the coupling kernels of steps 1-4 and of the revealed values'
   indicators; each factor of K reads only (v,v',t0), earlier auxiliaries and
   right values, never (x,x') directly. Fixing G's values and summing the right
   variables not in G leaves a factor free of (x,x'), so
   The step4 right-repair kernel is summed over its unrevealed coins and
   Ubar_i when forming this marginal; its total conditional mass is1. Thus
   Law(X,X'|G)=Law(X,X'|T0,V,V'). X is independent of T0 and (V,V') has 2v bits,
   so Hinf(X|G)>=k_C-2v-Ee outside mass 2^-Ee. Ubar_i stays uniform
   given G: step4 made it uniform conditional on H3 and the ENTIRE left
   record, and G reveals only H3 and that left record, excluding repair coins.
6. Reveal Z'_i=F(X',U'_i,encode(i)): U'_i is fixed in G, so this is a z-bit left
   revelation. With H5=(G,Z'_i), Hinf(X|H5)>=k_C-2v-z-2Ee>=k_F (I5) outside total
   mass 2*2^-Ee, and Ubar_i stays uniform.
7. F step. Given H5 the source tuple (X for the honest call and the other 7
   honest rows; X' for the 7 tampered rows j!=i) is independent of the seed
   tuple (Ubar_i; U_j,U'_j for j!=i). The honest seed Ubar_i is uniform on
   d_F=m1 bits (I7); encode(i) differs from every other label; the 14 tamperings
   are exactly t0=14,with exactly two copies of each nonhonest label. This verifies the selective-transcript contract on every conditioned branch; the tampered row i is fixed in H5 and is
   not a tampering. F's direct multicopy contract (Section4.2 and restricted8.2, zero shifts)
   makes Zbar_i=F(X,Ubar_i,encode(i)) 2^-Ee-close to uniform jointly with every
   other row output and every seed. With Qbar=Zbar_i+sum_(j!=i)Z_j and
   Q'=sum_j Z'_j (Z'_i fixed), D(Qbar|Q',seeds,H5)<=2^-Ee.
8. Retain Y,Y'. Given (H5, seed tuple), the outputs depend only on (X,X') and the
   seeds, and (X,X') is independent of (Y,Y', coins), since the seeds are
   right-sided given H5. So Law(Y,Y'|H5,seeds,outputs)=Law(Y,Y'|H5,seeds) is a
   common kernel of the real and ideal laws, and appending (Y,Y'), including
   unused seed suffixes, costs nothing [1.1]. Averaging over H5 and charging the
   two X tails gives synthetic error <=3*2^-Ee.
9. Restore. Real and synthetic (Q,Q',Y,Y') agree unless V!=Vbar, the matrix pair
   differs from (A,A'), W_i!=Wbar_i or U_i!=Ubar_i, so
   zeta<=eps0+e_A+delta_W+delta_U<=6*2^-Ee. By [1.1] (D_P<=D_Q+2zeta),
   D(Q|Q',Y,Y')<=3*2^-Ee+12*2^-Ee=15*2^(-E-4)<2^-E.

### 4.7 Running time

At actual length N, closed slots and B9 bound every NB/F#/CA# dimension
and enumerated field degree by O(logN); the primitive recurrences and
CZ/weak-AGHP setup obey the same bound. Each enumeration costs2^O(logN);
larger fields use towers. The fixed number of calls is therefore polynomial
in actualN, without a joint bound for unrelated n,z,a,E.

Attribution: Li (Section 3 and header), advice condenser and row-XOR;
weak-seed-hashing, AGHP.

## 5. Scope

Acyclic chain: PXT7g/PXG -> conditionaltools1--3, multicopy8.2 -> F#
(Section4.2); CZ/weak-AGHP -> advice condenser; that condenser, AGHP,
Ext0,Ext2,F# -> NB -> outerA# (proof.md, Lemma B9). Global guards and the
family/circuit assembly are in proof.md.

Quantitative credits: dbc05fc4-aed7-487a-b708-d1e5b938fbf4 supplies the
uniform-source LemmaU, collision-based W2, same-sideTR/TLaccounting. Publishedb15supplies103eta;current(80+7*2^-20)eta/conditional19eta/6 isdirectlyreproved with its explicit finite guard. Publishedc78568c1-1aa4-4da9-8556-877fb4218f6f supplies the whole-tuple MC-LA/retention method of its eight-round predecessor. The five-round schedule,separated-set invariant,75a/44h+45s/11s entropy ledger with a selective-transcript proof and at most two copies per nonhonest label and complete proof are supplied here in restricted8.2.
