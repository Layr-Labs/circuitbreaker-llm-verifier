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
Imports are the supplied finite algorithms and proofs mapped in sampler-and-disperser.md, with all guards; no asymptotic black boxes are used.
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
restricted8.3's ledger and CA# its separate eight-label ledger.
### 2.4 Actual primitive use
Selected linear primitives use PXT/PXG7g and PXTM7h with their actual recurrence and
explicit integer guards. The outer assembly uses Section2.1 withC=NB.
The selected innerF# uses the five-round antichain corollary8.3 with seven
labels and12tamperings,instantiated in4.2. Its zero-shift independent-
side joint contract is exactly the one used by NB4.6 step7. The full
conditional-average lift remains sampler6.1 and all original source/seed
quantifiers and retained padding suffixes remain explicit.
## 3. Finite advice-condenser reduction
### 3.3 Conditioning a near atom-cap law
Let P be the law of(W,T), e-close to a law Q with Q_T=P_T and
Q(W=w|T=t)<=2^-h. For an event E of P-probability p>0 and0<delta<=1, there is a
law (W*,T*) within delta+e/p of (W,T)|E, with the same T|E marginal and
  Hinf(W*|T*=t)>=h-log(1/p)-log(1/delta). (10)
Indeed excess(w,t)=[P(w,t)-2^-h P(t)]_+ sums to at most e, by the cap on Q and Q_T=P_T. Values t with
P(t|E)<delta P(t) have conditional mass <=delta. For the other t,
 P(w|t,E)<=2^-h/(p delta)+excess(w,t)/(p P(t|E)).
A probability vector's distance to atom cap c>=1/domain-size equals its
excess above c: remove that mass and distribute it below the cap, where
capacity suffices. Apply c=2^-h/(p delta), capped at 1 if necessary, average
over good t, and replace bad t with any capped law. This proves (10), even
if E is not determined by (W,T).
### 3.4 Six-input conditional-entropy advice map
Fix eta=31516696117353382584589982250/14379938999401373506041953467633 and integers r,a,sig>=1,u=r+a.
For six independent raw r-bit inputs, each of entropy>=r-sig, append the
fixed advice alpha to each and apply restricted Section8's E6. Use its
canonical prime p above2^u, and encode its two field coordinates W in
Rrow=2(u+1) bits. The partner uses the same parameters and alpha'. Assume
 (a+sig)/u<=eta/1738, r-sig>=(1-eta/1738)(u+1)+1,
 eta*u>=64*2^80, u>=2^90. (11')
Distinct advice makes all six coordinate partners fixed-point-free. These
are precisely the supplied E6 source/domain/prime guards. Its proved joint
law is within45p^-eta of a law with conditional honest entropy
 kappa=(7/6)eta*log2 p >=(7/6)eta*u.
We may require the partner marginal to equal the real Wprime marginal:
keep each ideal conditional law and change its marginal to the real one,
using uniform F_p^2 on a formerly absent partner fiber. This costs at most
another45p^-eta by contraction. The resulting law keeps conditional
entropy>=kappa and is within90*2^(-eta*u). No hashing input or seed is used.
Rrow is an actual integer width; the six raw rows are zero-padded to it.
The E6 and field algorithms are supplied in restricted Sections1--8;
at actual outer uses u=O(logactualN), so their setup is polynomial time.
### 3.5 Six-block uniform-source advice lemma U
Fix r,a,s,sig>=1, sig>=4,s>=7,r>=s+sig,u=r+a. Let X be
exactly uniform on6r bits and Xprime arbitrary, with distinct advice.
Assume(11') and
 eta*u>=6s+15sig+8. (U0')
Return six raw blocks zero-padded to Rrow, then W from3.4.
Thus there are seven actual rows. The matrix pair is within
14*2^-sig of a mixture whose selected row i<=7 has conditional entropy
at least s given its partner on every component and partner fiber.
Realize Xprime by a mixture of fixed functions of X, keeping X uniform.
Call a coordinate pair heavy when its probability is>=2^(-r-s-2sig).
Partition points with a nonheavy coordinate by the first such coordinate.
Discard pieces smaller than2^-sig, cost<=6*2^-sig. On a remaining piece,
each joint atom is<2^(-r-s-sig). Summing over at most2^r partner values
bounds average conditional maximum atom by2^(-s-sig). Outside partner
mass2^-sig the selected raw row has entropy>=s. Padding preserves it.
If the all-heavy set SH has mass<2^-sig, discard it. Otherwise
|SH|>=2^(6r-sig). Uniform raw marginals give at most2^(s+2sig) heavy
partners for EVERY honest coordinate. Its coordinate supports Si obey
|Si|>=2^(r-sig); with Pset=product Si, |SH|/|Pset|>=2^-sig.
Choose each coordinate partner map uniformly among its heavy partners.
Each SH point matches with probability>=2^(-6s-12sig). As long as the
remaining set has size>2^-sig|SH|, averaging gives a matched subset V
of density>2^-B in Pset, B=6s+14sig. Remove it and repeat; discarded
remainder mass is<=2^-sig. This is a finite analysis partition only.
On independent uniforms Si, use3.4 with deficit sig and coordinate maps.
Apply3.3's excess/conditioning argument to the near conditional-entropy
law, whose partner marginal is real. That proof only needs the atom cap
2^-kappa, not a uniform honest marginal. At delta=2^-sig it yields entropy
>=kappa-B-sig>=s, since(U0') implies kappa>=7s+15sig.
The cost is<=2^-sig+2^B*90*2^(-eta*u)
 <=(1+90/256)*2^-sig. Keep the conditioned partner marginal and append
all raw rows through their original conditional kernel given the W pair;
this preserves selected-row entropy and TV. Replace discarded pieces by any law with a uniform selected honest
row. Total losses are at most(6+1+1+1+1+90/256)*2^-sig<14*2^-sig.
No independence between honest and tampered raw coordinates is assumed.
## 4. The standard breaker NB
This section supplies the standard breaker of proof.md Section 3.1. It needs
a general source X and a UNIFORM seed Y only. The rate-theta advice row is used
as the SOURCE of the AGHP extractor of weak-seed-hashing Section 1, seeded by a
fresh slice of Y, so no rate-to-0.8 condensation and no ordinary condenser is
used. Inside this section F names NB's inner breaker, lg(x) is the least j>=0
with 2^j>=x, pow2ceil(x)=2^lg(x)<2x, and [1.1] cites Section 1.1.
### 4.1 Statement
Fix eta=31516696117353382584589982250/14379938999401373506041953467633,theta=eta/6,t0=12,a0=4. Use
Ext0=PXTM(JS,vNB0) andExt2=PXTM(JY,vNB2),
oflinear-and-bilinear-primitives.md Sections7g/7h. Ext0 is PXTM with7 mixed-arity levels and has
 c0=46154467221/4282519852296875,pi0=43580913511330227025807/137040635273500000000,
 6c0=276926803326/4282519852296875<theta=5252782686225563764098330375/14379938999401373506041953467633.
The sharper finite seed form uses p0=43580913511330227025807/137040635273500000000, D0=0,
so d0<=c0*v+p0*w_N+D0 by PTM.2.
Let F be the inner breaker F# of Section 4.2. For integers n,z,a,E>=1, Section
4.3 computes a seed length d_C and threshold k_C, and Section 4.4 defines a
deterministic NB: {0,1}^n x {0,1}^(d_C) x {0,1}^a -> {0,1}^z. If k_C<=n,
d2<=n and the advice guards (AG1) s>=a+40, (AG2) eta*r>=1738(a+sig+2)
of Section4.3 hold, then for every pair (X,X') on n bits independent of every pair (Y,Y')
on d_C bits, with Hinf(X)>=k_C, Y uniform, copies arbitrarily correlated
within each pair, and alpha!=alpha',
    D(Q|Q',Y,Y') <= 15*2^(-E-4) < 2^-E, Q=NB(X,Y,alpha), Q'=NB(X',Y',alpha'). (NB)
So NB is a standard breaker in the sense of proof.md Section 3.1 at every
a,z,E>=1, with the guards k_C<=n, d2<=n, (AG1) and (AG2), retaining both complete seeds. If a
guard fails the map outputs 0^z and no claim is made (the totality convention
of proof.md Section 3).
### 4.2 Five-round antichain inner breaker F#
For source lengthn,requested outputz>=1,advice among NB labels1,...,7,
and exponentE_F>=1, compute the public INTEGER
 w=max(2^80,E_F+2lg(4n)+32,ceil(z/128)),e=E_F+7,q=e,t=12,c=13.
Every actual slot is ceil(A*w),with these fixed rational coefficients:
 a=ceil((5163898551600000001/200000000000000000)w).
 s=ceil((4892360580400000001/200000000000000000)w).
 s0=ceil((6454128230458398527/250000000000000000)w).
 hc=ceil((299837216629350022733/1000000000000000000)w).
 d=ceil((1059041890539034873923/500000000000000000)w).
 k=ceil((2726577018902136557601609/1000000000000000000)w).
Write h=hc. Use these exact Section7g PXG maps,all at exponent e:
 I=PXG(vMI)(inputn,outputh,seeds0),
 F=PXG(vMF)(inputn,outputh,seeda),
 G=PXG(vMG)(inputd,outputa,seeds),
 J=PXG(vMJ)(inputh,outputs,seeda).
Pad their seeds to the integer slots and retain suffixes. Ifk>n or the
label is outside1,...,7,return0^z. Otherwise use the ACTUAL five-round
algorithm of restricted8.3: startW=I(X,prefix_s0 Y);at each of the five
coordinates computeT=prefix_s W,A=G(Y,T),Z=J(W,A),B=G(Y,Z),chooseC=B
for codebit1 andC=A for codebit0,andsetW=F(X,C). Returnprefix_z W.
The seven fixed label codes are11000,10100,10010,01100,01010,00110,10001.
For NB's4-bit encode(j),j1,...,7,the corresponding code is the listedjth
word.
Putd_F=d,k_F=k. Restricted8.3 proves this actual algorithm for an arbitrary
source13-tuple independent of an arbitrary seed13-tuple, honestentropy>=k,
honestseed uniform,and honestlabel different from all others. Correlations inside each tuple are arbitrary; each nonhonest label occurs at most twice. NB has exactly two copies of each other label on every conditioned branch. Error is
76*2^(-E_F-7)<2^-E_F jointly with ALL12 other outputs and ALL13 complete
original seeds. NB4.6step7 uses this zero-shift contract.
Put delta=10^-18. For integerw>=2^80,ceil(Aw)<=(A+delta)w;
the left slot is at least its coefficient timesw.
The coefficients (c,pi,Lam,g) are the Section7g recurrences at vMI,vMF,vMG,vMJ,
using pi=pi_* from the legal-call PXG* seed bound. First the PT1 threshold
rows establish K0<=u; then the PXG* seed rows apply, without a circular premise.
With x'=x+delta for each slot coefficient x, the unrounded rows are
 aJ=c_J*s'+pi_J, aF=c_F*h'+pi_F, sG=c_G*a'+pi_G,
 s0I=c_I*h'+pi_I, hJ=(9+Lam_J)s'+g_J+1, hcap=128, [h>=s]=s',
 dG=13s0'+(65+Lam_G)a'+g_G+1, [d>=s0]=s0',
 kF=(38+Lam_F)h'+39s'+g_F+1, kI=Lam_I*h'+g_I, [k>=d]=d'.
Their displayed numerical bounds are rounded UP to grid 10^-18:
 aJ: [25819492758000000003/1000000000000000000]w <= (5163898551600000001/200000000000000000)w.
 aF: [12909746378916797053/500000000000000000]w <= (5163898551600000001/200000000000000000)w.
 sG: [24461802902000000003/1000000000000000000]w <= (4892360580400000001/200000000000000000)w.
 s0I: [12908256460916797053/500000000000000000]w <= (6454128230458398527/250000000000000000)w.
 hJ: [299837216629350022731/1000000000000000000]w <= (299837216629350022733/1000000000000000000)w.
 hcap: [128]w <= (299837216629350022733/1000000000000000000)w.
 h>=s: [12230901451000000003/500000000000000000]w <= (299837216629350022733/1000000000000000000)w.
 dG: [529520945269517436961/250000000000000000]w <= (1059041890539034873923/500000000000000000)w.
 d>=s0: [25816512921833594109/1000000000000000000]w <= (1059041890539034873923/500000000000000000)w.
 kF: [2726577016766440756266743/1000000000000000000]w <= (2726577018902136557601609/1000000000000000000)w.
 kI: [2726577018902136557601607/1000000000000000000]w <= (2726577018902136557601609/1000000000000000000)w.
 k>=d: [2118083781078069747847/1000000000000000000]w <= (2726577018902136557601609/1000000000000000000)w.
The rows prove bothaseeds,s,s0,h>=(9+Lam_J)s+(g_J+1)w,h>=128w>=z,
h>=s,d>=13s0+(65+Lam_G)a+(g_G+1)w,d>=s0,
k>=(38+Lam_F)h+39s+(g_F+1)w,k>=K_I andk>=d.
These are exactly selective8.3 at2c+t=38,3c=39,5c=65 andnu<=9,includingq<=w. They implyKJ<=h,
KG<=d,KI,KF<=k;K_I>=h,so all primitive inputs are<=n whenk<=n.
Selectedheights<=8 ande=E_F+7 give
 e+h_primitive+6+2lg(4u)<=E_F+21+2lg(4n)<=w.
The PXG/PXTM source,seed,ambient and field guards are proved in7g--7h.
The PUBLIC output cap inw defines the map for everyz. Degrees are logarithmic;larger relative extensions use towers,as7a/7g.
### 4.3 Parameter algorithm
Put Ee=E+4 and the unit w_N=max(2^80,Ee+38+6+2lg(4n)).
1. (d_F,k_F): F# at source length n, output z, advice 4, exponent Ee (4.2);
   m1:=d_F.
2. m0,K2: the exact seed length and threshold of Ext2 at input length n,
   output m1, exponent Ee.
3. sig=Ee+4; r=ceil((m0+p0*w_N+D0+(11/2)Ee+53)/(theta-6c0));
   s=floor((eta*r-15sig-8)/6); b=a;
   v=6r.
4. d0,K00: the exact seed length and threshold of Ext0 at input n, output v,
   exponent Ee.
5. Rrow=2(r+a+1); e_L=ceil((d0+Ee-1)/2)+Ee; r_L=lg(Rrow*m0+1)+m0+2e_L; d1=2r_L.
6. d2=d0+d1+m1+K2+Ee.
7. d_C=d0+d1+d2; k_C=max(K00, k_F+2v+z+2Ee).
The guards are k_C<=n, d2<=n and the advice guards
 (AG1) s>=a+40, (AG2) eta*r>=1738(a+sig+2) (Section 4.1), public integer tests
 on the computed s,r,a,sig. r has no advice term; proof.md Lemma B9 discharges both.
### 4.4 Evaluator
Slice y=y0y1y2 with lengths d0,d1,d2. Components:
* Ext0: the strong linear primitive of step 4, output truncated to v bits.
* C_adv: Section3.5 on six r-bit blocks, accuracy sig=Ee+4,
  actual advice length a. It returns seven Rrow-bit rows; its selected good row has index<=7.
* LinExt: supplied AGHP from Rrow bits to m0 bits, exponent e_L,
  canonical degree r_L and complete seed d1=2r_L bits.
* Ext2: the strong primitive of step 2 applied to y2 zero-padded to n bits,
  with seed W_j (m0 bits) and output m1=d_F bits.
* F: F# of Section 4.2; it reads all of U_j.
    V=Ext0(x,y0); (R_1..R_7)=C_adv(V,alpha);
    W_j=LinExt(R_j,y1); U_j=Ext2(y2,W_j); Z_j=F(x,U_j,encode(j)), j=1..7;
    Q=Z_1+...+Z_7, (NBE)
with encode(j) the 4-bit binary label of j. The same map gives Q' on primed
inputs. The proof fixes a good row i inside each mixture component; the
evaluator XORs all rows with distinct labels and never needs i.
AGHP at its exact integer degree. Weak-seed Section1 defines
r=lg(N+1)+m+2e over any r-dimensional binary field presentation,
with its coordinate dot product. Step5 uses exactly this degree at
N=Rrow*m0,m=m0,e=e_L. Its root-counting and strongness proof therefore
applies directly, as does the squared-collision weak-seed LemmaW2.
Choose the first irreducible degree-r_L binary polynomial as in weak-seed
Section4; at every selected actual outer use r_L=O(logactualN), so its
2^r_L*poly(r_L) setup is polynomial in actualN. Source-linearity and full seed
retention are unchanged, and d2 is computed from this actual d1.
### 4.5 Finite inequalities, each proved here
* (I1) PXT gives d0<=6c0*r+p0*w_N+D0 and2e_L<=d0+3Ee.
  The r and s definitions imply
  s>=m0+6c0*r+p0*w_N+D0+3Ee+122/3.
  Thus s>=m0+2e_L; s>=a+40 is guard (AG1). The denominator theta-6c0 is positive.
* (I2) eta*r>=6m0+6p0*w_N+6D0+33Ee+318, so eta*r>=64;
  guard (AG2) gives eta*r>=1738(a+sig+2).
  Hence(a+sig)/(r+a)<=eta/1738 and
  r-sig>=(1-eta/1738)(r+a+1)+1. Also eta<=1/453 and the s definition
  gives r>=s+sig and eta*(r+a)>=6s+15sig+8.
  Thus(U0') holds. The explicit p0>=16,w_N>=2^80 give
  eta*r>=6p0*w_N>=96*2^80>64*2^80 and u>=2^90.
  All E6 domain/prime guards hold: eta*r>=96*2^80 and eta<=1/453
  imply r>=43488*2^80>2^90. Its conditional entropy kappa is
  more than7s+15sig by(U0'), paying the conditioning loss. LemmaU
  gives error<=14*2^-sig<2^-Ee.
* (I3) Ext2 is legal: K2<=d2<=n; the conservative unit w_N covers
  actual binary depths17 and7 of both PXTM maps. The prime guard
  u=r+a>=r>=6p0*w_N/eta>=6*2^80/eta>2^90, sincep0>=1,w_N>=2^80; this proves the independent prime guard rather than inferring it from64/eta. Also eta*r>=6*p0*w_N>=96*2^80>64*2^80,usingp0>=16,w_N>=2^80,so the strengthened CZ constant-absorption guard holds at EVERY actual call.
* (I4) d2-d0-d1-m1-Ee=K2 (step 6).
* (I5) k_C-2v-z-2Ee>=k_F (definition of k_C).
* (I6) ((d0+Ee-1)/2)-e_L<=-Ee by definition of e_L; Lemma W2 applies.
* (I7) d_F=m1 (step 1).
* (I8) Ext0's output is at least v (its leaf); and K00<=k_C<=Hinf(X)<=n.
* (I9) k_F<=k_C<=n, so F's guard holds; F sees 12<=t0 tamperings below.
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
2. For each T0 apply Lemma U to(Vbar,V'): Vbar is exactly uniform on6r bits,
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
7. F step. Given H5 the source tuple (X for the honest call and the other 6
   honest rows; X' for the 6 tampered rows j!=i) is independent of the seed
   tuple (Ubar_i; U_j,U'_j for j!=i). The honest seed Ubar_i is uniform on
   d_F=m1 bits (I7); encode(i) differs from every other label; the 12 tamperings
   are exactly t0=12,with exactly two copies of each nonhonest label. This verifies the selective-transcript contract on every conditioned branch; the tampered row i is fixed in H5 and is
   not a tampering. F's direct multicopy contract (Section4.2 and restricted8.3, zero shifts)
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
CZ/E6 setup obey the same bound. Each enumeration costs2^O(logN);
larger fields use towers. The fixed number of calls is therefore polynomial
in actualN, without a joint bound for unrelated n,z,a,E.
Attribution: Li (Section 3 and header), advice condenser and row-XOR;
weak-seed-hashing, AGHP.
## 5. Scope
Acyclic chain: PXT/PXG7g and PXTM7h -> conditionaltools1--3, multicopy8.2--8.3 -> F#
(Section4.2); CZ/E6 -> advice condenser; that condenser, AGHP,
Ext0,Ext2,F# -> NB -> outerA# (proof.md, Lemma B9). Global guards and the
family/circuit assembly are in proof.md.
Quantitative credits: dbc05fc4-aed7-487a-b708-d1e5b938fbf4 supplies the
uniform-source LemmaU, collision-based W2, same-sideTR/TLaccounting. Publishedb15supplies103eta;current(37373/869+7*2^-80)eta/conditional7eta/6 isdirectlyreproved with its explicit finite guard. Publishedc78568c1-1aa4-4da9-8556-877fb4218f6f supplies the whole-tuple MC-LA/retention method of its eight-round predecessor. The five-round schedule,separated-set invariant,65a/38h+39s/9s selected entropy ledger with a selective-transcript proof and at most two copies per nonhonest label and complete proof are supplied here in restricted8.2--8.3.
