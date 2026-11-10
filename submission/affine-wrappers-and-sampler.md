# Proof layout
The wrapper and sampler parts below retain the aliases
affine-and-standard-wrappers.md and sampler-and-disperser.md in internal
references. Sampler SP.1,SP.3--SP.4 and SP.7--SP.8 are supplied in the sampler
appendices of linear-and-bilinear-primitives.md. Historical archive filenames
are unchanged;these aliases do not request missing attachments.

# Wrapper part
# Finite affine transformation and breaker wrappers
Logs are binary,TV is half L1,+ is XOR,U_b is independent
uniform,D(V|H)=TV((V,H),(U_|V|,H)),and
Havg(X|H)=-log sum_h P(h)max_x P(x|h). Only positive-probability conditional
branches are used. Lengths are integers;entropy bounds are literal atom caps.
Attribution: CL's affine transformation (Theorem 6.1),
https://arxiv.org/abs/2110.12652v1 ;Li's generalized copies,advice reduction
and equal/general-length wrappers (Definitions 3.3--3.4,Lemmas 4.2,5.1,
Theorem 5.2),https://eccc.weizmann.ac.il/report/2023/023/revision/1/download/ .
BKSSW condensation: https://eccc.weizmann.ac.il/report/2010/037/download .
CZ's restricted polynomial core:
https://eccc.weizmann.ac.il/report/2014/102/download/ . The local replacements
are credited in their modules.
## 1. Imports and elementary tools
Imports are the supplied finite algorithms and proofs mapped in sampler-and-disperser.md,with all guards;no asymptotic black boxes are used.
### 1.1 Probability facts
If full L,R are independent given H,a revelation with kernel depending
only on(L,H) multiplies only the left density factor. Thus conditional
product survives at every positive-probability history;alternate sides by
induction. Coins may be included in the side tuple. A mixed XOR becomes
left-only after its right contribution is fixed.
For T with at most2^b values,
sum_t P(t)max_x P(x|t)=sum_t max_x P(x,t)<=2^b max_x P(x).
Conditioning on H first gives Havg loss<=b,also for randomized T.
Markov on maximum conditional atom gives Hinf(X|H=h)>=K outside
mass2^-s whenever Havg(X|H)>=K+s. Right-only revelations leave left laws.
TV with a common history marginal is average conditional TV;kernels
contract by half-L1. Two joint laws at distance zeta obey
D_P(V|H)<=D_Q(V|H)+2zeta by triangle,including restoring H's marginal.
Changing only V under a coupling with H fixed costs its disagreement once.
Common mass min(P(z),Q(z)) constructs maximal coupling via a kernel from
sampled Z. If Z=f(A),sampling auxiliaries with kernelK(aux|Z,H) preserves
Law(A|Z,H,aux)=Law(A|Z,H);a right tuple independent of A given H stays so.
Right revelations with auxiliaries fixed preserve that identity by the
same product factorization. Labels need not have bounded ranges.
If S=f(L),R independent of L,W=F(S,R),then D(W|L)=D(W|S) because
its conditional law depends only on S. The conditional-product version
lets seed-retaining extraction retain the entire left tuple.
### 1.2 Strong merging and robustness
Let source tuple (U,U_1,...,U_t) and seed tuple (S,S_1,...,S_t) be independent
given H. Let F be strong (K,eps),output b bits,2^-q<=eps. For I,J subset
[t],suppose
D(S|H,S_I)<=eta,
Havg(U|U_J,H)>=K+|I\J|b+q. (1)
With V=F(U,S),V_j=F(U_j,S_j),
D(V|V_(I union J),S,S_1,...,S_t,H)<=eta+2eps. (2)
Proof: reveal U_J,S_I,then V_(I\J). The latter are source-only with S_I
fixed,costing |I\J|b;the product split survives. Entropy fails on at most
eps mass. Uniform-seed extraction costs eps elsewhere. On each such history,
conditional output-distance with the seed exposed is a [0,1]-valued function
of S;changing its marginal to uniform costs its conditional seed distance.
Source revelations are independent of S given H,S_I,so their averaged cost
remains eta. Append remaining seeds by their conditional kernel given H,S_I,S,
independently of the source. Outputs in J are determined by U_J and the retained
seeds. Discard surplus history to prove (2).
The same expectation argument for an all-seeds-retained breaker of error e
gives error e+zeta when its good seed marginal is zeta from uniform: preserve
the kernel of other seeds given the good one. Average source entropy K+q adds
at most eps for bad histories. Thus a strong ordinary breaker of error eps
has robust error eta+2eps. Extra seed copies can be retained when independent
of the full source pair given the retained seed pair/history,as in §2.
## 2. Generalized all-seeds-retained affine transformation
### 2.1 Selected affine-left contract and algorithm
Fix n,d,t,a,m>=1,eps in (0,1),q=ceil(log(1/eps)),H=ceil(log t).
Supply these deterministic maps,each of error <=eps:
C: {0,1}^d x {0,1}^d0 x {0,1}^a -> {0,1}^dx;
L0: {0,1}^n x {0,1}^s0 -> {0,1}^d0,threshold K0;
Lr: {0,1}^n x {0,1}^dx -> {0,1}^r,threshold Kr;
E: {0,1}^d x {0,1}^dy -> {0,1}^dx,threshold KE;
Lm: {0,1}^r x {0,1}^dx -> {0,1}^dy,threshold Km.
Every L is strong on uniform affine sources and source-linear at fixed seed;E is strong on arbitrary sources.
C is a strong ordinary one-tampering breaker of threshold KC: arbitrary
correlated source pair independent of arbitrary correlated seed pair,honest
seed uniform,distinct advice,both complete seeds retained. Its first input
really has length d. Padded seeds may ignore suffixes while retaining them.
All actual primitive guards must hold. Put
TR0=(t+1)(s0+d0+dx),TL0=(t+1)d0,T0=TR0+TL0,
TR=TR0+H(t+1)(2dy+2dx)-(t+1)dx,TL=TL0+H(t+1)2dy,T=TR+TL+(t+1)dx.
Require positive legal dimensions,k<=n,and
s0<=d,dy<=r,m<=r,
r>=(t+1)dy+Km,
k>=max(K0,TL+Kr+t*r),
d>=(t+1)(s0+d0)+KC+q,
d>=TR+KE+t*dx+q. (3)
The evaluator on (x,y,alpha) is
s=Prefix_s0(y);r1=L0(x,s);s2=C(y,r1,alpha);w0=Lr(x,s2).
For i=1,...,H:
p=Prefix_dy(w_(i-1));qm=E(y,p);v=Lm(w_(i-1),qm);
qr=E(y,v);w_i=Lr(x,qr).
Return Prefix_m(w_H).
The outer A# of proof.md outer Section3 has m<=r,so no long-output
extension is used here. Section4.2 constructs F# separately from restricted8.2.
The selected tuple is A_i=T_i Z+c_i for uniform affine Z,invertible T_i,
independent of arbitrary correlated(B,Y),with honest Y0 uniform and
distinct honest advice. Section2.5 proves its complete retained-right
contract. No arbitrary-left theorem is asserted by this selected wrapper.
### 2.2 Initialization
Write L=A tuple,R=(B,Y) tuple. Fixed-seed linearity splits every L call
into A and B contributions,e.g. R1_j=R1A_j+R1B_j. Reveal all
S_j=Prefix_s0(Y_j),then R1B_j;call this H0,length <=(t+1)(s0+d0).
Strong affine extraction and independence of A_0 from the entire right tuple give
D(R1A_0|all S)<=eps. Given S,appending R1B is independent of A_0;translate
by R1B_0 to obtain D(R1_0|H0)<=eps. The entire R1 tuple is left-sided given
H0 and independent of the Y tuple.
Y_0 has average entropy >=d-(t+1)(s0+d0)>=KC+q. §1.2 applied to C
thus gives D(S2_0|S2_j,H0,R1_0,R1_j)<=3eps for every j>=1.
Other R1 coordinates can be appended because,given H0 and the retained pair,
they are independent of the full Y tuple and hence both C outputs. Reveal
all R1,a left revelation since R1B was fixed. Then S2 is right-sided.
Apply (2) to Lr with sources A,seeds S2,I={j},J empty. The k bound in
(3) supplies Kr+r after the history;the exact-affine lift of2.5 gives the conservative error5eps. Let Z0 contain
S,R1B,R1,S2 across all copies,of length <=T0. The full sides remain
independent and
D(W0A_0|W0A_j,Z0)<=5eps.
### 2.3 Histories and causal ledger
After round i,set Zi by appending(PB,P,Qm,VB,V,Qr) across all copies in
that order. PB,Qm,VB,Qr are right messages,P,V are left messages once
earlier opposite contributions are fixed. The full sides remain independent
at each history by1.1. The history length is at most
T0+i(t+1)(4dy+2dx). Each actual left message is affine on its conditional
left fiber. The statistical doubling proof is2.5;the ledger below
charges the source budgets used there.
Causal same-side ledger. On a rectangular history a message from the other
side leaves the conditional source law unchanged. An adaptive b-bit own-side
message multiplies average maximum atom by at most 2^b: for each old history
sum maximum joint atoms over its at most 2^b messages,then average. Induction
thus charges only own-side bits,even for history-dependent functions.
Initialization sends S,R1B,S2 right (TR0 bits) and R1 left (TL0 bits).
Each round sends PB,Qm,VB,Qr right ((t+1)(2dy+2dx) bits),and P,V left
(2(t+1)dy bits). Before the second E call of round i,that round's Qr has
NOT been revealed: right cost is TR0+i(t+1)(2dy+2dx)-(t+1)dx<=TR.
The first E call precedes Qm and VB as well,so costs less. If H=0 there are no E calls,so the E guard d-TR>=KE+t*dx+q is vacuous and TR's -(t+1)dx term is used only for H>=1. C is called before R1,S2 are revealed,at right cost
(t+1)(s0+d0). Every Lr call has left cost <=TL. These are exactly the
right source constraints(3),with q;the affine left guard has no tail. Lm uses the exact affine source guard of2.5. The full transcript length is T: Qr is subsequently revealed in full,
so rectangularity,the invariant of2.5,and retention of ALL original
seeds are unchanged. The outer assembly uses this ledger;F# uses
restricted8.4's ledger and CA# its separate eight-label ledger.
### 2.4 Actual primitive use
Selected L0/Lr use linear7f-local (Lr retains its padded suffix);Lm uses AGHPaff. Other primitives use PXT/PXG7g and PXTM7h with explicit integer guards. The outer assembly uses Section2.1 withC=NB.
The selected innerF# uses the four-round antichain corollary8.4 with five
labels and8tamperings,instantiated in4.2. Its zero-shift independent-
side joint contract is exactly the one used by NB4.6 step7. The full
conditional-average lift remains sampler6.1 and all original source/seed
quantifiers and retained padding suffixes remain explicit.
### 2.5 Affine-left rank specialization
Let A_i=T_i Z+c_i for uniform affine Z,invertible T_i. Use evaluator2.1 with its strong affine L0,Lr,Lm and guards
r>=(t+1)dy+Km, k>=max(K0,TL+Kr+t*r). (3aff)
All other guards,including right C/E entropy tails q,remain. Then the
error is at most delta_H=(9*3^H-4)eps. At H=9 the multiplier is
177143<2^18. This is an affine-left theorem;no arbitrary-left theorem is asserted. Only this supplied affine-left contract is consumed.

An exact-affine version of1.2 needs no q when its conditional source
is a linear image of a uniform affine variable and all revealed source
coordinates are affine. Honest entropy>=K+|I\J|b on every fiber
implies entropy>=K after the source-only tampered outputs,on EVERY
nonempty fiber:each b-bit affine constraint reduces the dimension of
any linear image by at most b. The same seed-distance/extraction proof
of(2),without Markov,gives eta+eps;use the conservative eta+2eps.
Source-coordinate maps may be noninvertible. Fix retained tampered
seeds S_I first,so their source-output constraints have fixed affine
coefficients. Other seeds append by their conditional kernel as in1.2.
This lift uses only a strong affine-source extractor when every fiber
is affine,as for the selected Lm.

On an actual rectangular history all left messages R1,P,V are affine
in Z,by source-linearity with earlier right messages fixed. The remaining
left variable is uniform on an affine fiber,of codimension at most the
left bits revealed. Right messages preserve this law. Thus every Lr
call,including at most t tampered r-bit outputs,has the exact entropy
required by(3aff);initial L0 needs only k>=K0. Histories and own-side
ledgers are those of2.2--2.3.

Maintain D(WiA_0|WiA_J,Zi)<=delta_i for |J|<=min(2^i,t).
Initialization2.2 using the exact-affine Lr lift keeps the conservative
bound delta_0=5eps. At any fixed history z and set J,the rank of WA_0
conditional on WA_J is constant over nonempty WA_J-value fibers.
If it is full r,WA_0 is uniform independent of WA_J. Otherwise every
conditional law is uniform on a proper affine subspace and has distance
at least1/2 from U_r. Hence the rank-bad HISTORY event has probability
at most2delta for a tuple at distance delta. It is a function of old
history,not of a sampled value of WA_J.

For a round write J=J1 union J2 with each size<=u_(i-1),delta=delta_(i-1).
Right PB revelation preserves conditional left laws. The previous TV
invariant gives D(P_0|P_J1,Z,PB)<=delta. Ordinary right E,with its q
tail,gives Qm seed error<=delta+2eps retaining all P;call this H1.
Restrict to old histories rank-good for J2,whose complement has mass
at most2delta. Given WA_J2,Z,PB,the P coordinates in J2 are determined;
remaining P cost at most(t+1-|J2|)dy exact affine codimensions in WA_0,
initially of dimension r. With Qm fixed,the tampered Lm outputs in
J1\J2 cost at most|J1\J2|dy additional codimensions. Total cost is
at most(t+1)dy,so(3aff) supplies Km on every nonempty fiber.
Apply the exact-affine lift to the ACTUAL WA tuple,I=J1,J=J2. Its
weighted seed-distance contribution on good old histories is bounded
by the unconditional delta+2eps:all terms are nonnegative and the event
is in old retained Z. No conditional error is renormalized. Add at
most2eps extraction and2delta bad mass,giving VA distance<=3delta+4eps.
No auxiliary W* coupling is used. Append/translate right contributions,
then ordinary right E adds2eps and exact-affine Lr conservatively adds2eps.
Thus delta_i=3delta_(i-1)+8eps=(9*3^i-4)eps. All averages use the
actual history marginal and fixed-history rectangular independence.
At H append the full right tuple to the final left-output bound using
that product,translate WB,take prefixes and discard surplus history.
This retains every original seed and the entire right-side variable.
Only the selected affine Lm class,r and k changed;right entropy bounds,
transcript lengths and distinct-advice requirements remain as before.
## 3. Finite advice-condenser reduction
### 3.3 Conditioning a near atom-cap law
Let P be the law of(W,T),e-close to a law Q with Q_T=P_T and
Q(W=w|T=t)<=2^-h. For an event E of P-probability p>0 and0<delta<=1,there is a
law (W*,T*) within delta+e/p of (W,T)|E,with the same T|E marginal and
Hinf(W*|T*=t)>=h-log(1/p)-log(1/delta). (10)
Indeed excess(w,t)=[P(w,t)-2^-h P(t)]_+ sums to at most e,by the cap on Q and Q_T=P_T. Values t with
P(t|E)<delta P(t) have conditional mass <=delta. For the other t,
P(w|t,E)<=2^-h/(p delta)+excess(w,t)/(p P(t|E)).
A probability vector's distance to atom cap c>=1/domain-size equals its
excess above c: remove that mass and distribute it below the cap,where
capacity suffices. Apply c=2^-h/(p delta),capped at 1 if necessary,average
over good t,and replace bad t with any capped law. This proves (10),even
if E is not determined by (W,T).
### 3.4 Mixed-source conditional-entropy advice map
Fix the selected eta of4.1 and r,a,sig>=1,u=r+a. Use four independent raw
blocks of lengths r,r,2r-m,2r,m=floor(eta*r/2),each with entropy
at least its length minus sig.
Split the third block into r,r-m and the fourth into r,r bits;zero-pad
the third block's second part to r,then append the same advice to all six coordinates. Apply restricted Section8's
canonical prime p>2^u and mixed-source E6;encode W in Rrow=2(u+1) bits.
The partner has the same parameters and different advice alpha'. Assume
(a+sig)/u<=eta/(11321180837/12500000),r-sig>=(1-eta/(11321180837/12500000))(u+1)+1,
eta*u>=64*2^80,u>=2^90. (11')
The scalar guard also gives2r-sig>=2(r-sig)>=(2-2eta/(11321180837/12500000))(u+1)+2,
This pays Z56. Subtracting m<=eta*(u+1)/2 pays Z34
with threshold(2-(11371180837/22642361674)eta)(u+1)+1,as restricted Section8 proves. Distinct advices make the two
scalar maps fixed-point-free. The joint partner maps may be arbitrary.
The supplied theorem gives error45p^-eta and conditional honest entropy
kappa=(5/4)eta*logp>=(5/4)eta*u.
Change the ideal partner marginal to the real one,retaining each ideal
conditional law and using uniform Kp on formerly absent partner fibers.
Contraction bounds the additional cost by45p^-eta. The resulting law
keeps entropy>=kappa and is within90*2^(-eta*u),with real partner marginal.
No hashing input or extra randomness is used. The field algorithms are
supplied;at actual calls u=O(log actual N),so their setup is polynomial.
### 3.5 Four-block uniform-source advice lemma U
Fix r,a,s,sig>=1,sig>=4,s>=7,r>=s+sig,u=r+a. Put m=floor(eta*r/2),v=6r-m. Let X be uniform onv bits
and Xprime arbitrary,with distinct advice. Group their raw coordinates
as lengths ell=(r,r,2r-m,2r). Assume(11') and
eta*u>=4s+11sig+8. (U0')
Return five Rrow-bit rows in this FIXED order: the four raw groups
zero-padded,encoded W from3.4. Both joint raw rows have length<=2r and fit Rrow. The matrix pair is within14*2^-sig of a mixture with a
selected row i<=5 having entropy>=s given its partner on every component
and partner fiber. All five rows participate in the five-label NB interface.
Represent Xprime by a mixture of deterministic functions of X while
keeping X uniform. A group-pair is heavy if its probability is at least
2^(-ell_i-s-2sig). Partition points with a nonheavy group by the first
such group. Discard pieces smaller than2^-sig,costing<=4*2^-sig.
On a remaining piece each joint atom for that group is below
2^(-ell_i-s-sig). Summing over at most2^ell_i partner values bounds its
average conditional maximum atom by2^(-s-sig). Outside partner mass
2^-sig the selected raw row has entropy>=s. Padding is injective.
If the all-heavy set SH has mass<2^-sig,discard it. Otherwise
|SH|>=2^(v-sig). Uniform marginals give at most2^(s+2sig) heavy partner
values for every honest group value. Its supports Si satisfy
|Si|>=2^(ell_i-sig);for Pset=product Si,|SH|/|Pset|>=2^-sig.
Choose one partner function per group uniformly among its heavy choices.
Each SH point matches with probability>=2^(-4s-8sig). While the remaining
set exceeds2^-sig|SH|,averaging gives a matched V with density>2^-B
in Pset,B=4s+10sig. Remove it and repeat,discarding remainder mass
<=2^-sig. This partition is analysis only;it is not computed by C_adv.
On independent uniforms Si use3.4. The two joint partner functions are
allowed by its mixed-source theorem;no independence inside them is
assumed. Lemma3.3 conditions the near atom-cap law on V,with real partner
marginal and delta=2^-sig. Its honest entropy is at least
kappa-B-sig>=s,since(U0') gives kappa>=5s+(55/4)sig+10>=5s+11sig.
Its cost is at most2^-sig+2^B*90*2^(-eta*u)<=(1+90/256)*2^-sig.
Keep the conditioned partner marginal and append the raw rows using their
original conditional kernel given the W pair.
This preserves selected-row entropy and TV. On discarded pieces use any
law with a uniform selected honest row. Total losses are at most
(4+1+1+1+1+90/256)*2^-sig<14*2^-sig. Honest and tampered groups may be
arbitrarily correlated with each other;only the honest input is uniform.
## 4. The standard breaker NB
This section supplies the standard breaker of proof.md §3.1. It needs
a general source X and a UNIFORM seed Y only. The rate-theta advice row is used
as the SOURCE of the canonical extractor7g,with weak-seed bound W2,seeded by a
fresh slice of Y,so no rate-to-0.8 condensation and no ordinary condenser is
used. Inside this section F names NB's inner breaker,lg(x) is the least j>=0
with 2^j>=x,pow2ceil(x)=2^lg(x)<2x,and [1.1] cites §1.1.
### 4.1 Statement
Put DG=11321180837/12500000. Fix eta=((2-nu)/(14D_E+h_R+5/2)),theta=eta/4,t0=8,a0=4. Use
Ext0=PXTM(JS,vNB0) andExt2=PXTM(JY,vNB2),
oflinear-and-bilinear-primitives.md Sections7g/7h. Ext0 is PXTM with6 mixed-arity levels and has
c0=1/15000,pi0=10043823040395303/50000000000000,
bc0=420896407082656137773873317112941211/1052961358118344326462975118002150864000<theta=1440680823407964056583650439595673/701974238745562884308650078668100576.
Put b=6-eta/2;bc0 means b*c0. Since c0<1,v<=b*r+1 implies
d0<=bc0*r+p0*w_N+D0+c0.
The sharper finite seed form uses p0=10043823040395303/50000000000000,D0=0,
so d0<=c0*v+p0*w_N+D0 by PTM.2.
Let F be the inner breaker F# of §4.2. For integers n,z,a,E>=1,Section
4.3 computes a seed length d_C and threshold k_C,and §4.4 defines a
deterministic NB: {0,1}^n x {0,1}^(d_C) x {0,1}^a -> {0,1}^z. If k_C<=n,
d2<=n,the primitive guards (PG) below,and the advice guards (AG1) s>=a+40,(AG2) eta*(r+a)>=(11321180837/12500000)(a+sig+2)
of Section4.3 hold,then for every pair (X,X') on n bits independent of every pair (Y,Y')
on d_C bits,with Hinf(X)>=k_C,Y uniform,copies arbitrarily correlated
within each pair,and alpha!=alpha',
D(Q|Q',Y,Y') <= 15*2^(-E-4) < 2^-E,Q=NB(X,Y,alpha),Q'=NB(X',Y',alpha'). (NB)
So NB is a standard breaker in the sense of proof.md §3.1 at every
a,z,E>=1,with the guards k_C<=n,d2<=n,(PG),(AG1) and (AG2),retaining both complete seeds. If a
guard fails the map outputs 0^z and no claim is made (the totality convention
of proof.md §3).
### 4.2 Four-round antichain inner breaker F#
For source lengthn,requested outputz>=1,advice among NB labels1,...,5,
and exponentE_F>=1,compute the public INTEGER
w=max(2^80,E_F+2lg(4n)+32,ceil(z/128)),e=E_F+7,q=e,t=8,c=9. Writew_F=w.
Every actual slot is ceil(A*w),with these fixed rational coefficients:
a=ceil((4344053288370334047/250000000000000000)w).
az=ceil((885064049710693167/62500000000000000)w).
s=ceil((7039201288905882671/500000000000000000)w).
s0=ceil((4344028571587233297/250000000000000000)w).
hc=ceil((22185866417274888401/200000000000000000)w).
d=ceil((159171125524877564171/200000000000000000)w).
k=ceil((2480712103654093997382983/1000000000000000000)w).
sz=ceil((14078402577811765339/1000000000000000000)w).
Write h=hc. Use the exact7h mixed maps and7g binary maps at exponent e:
I=PXTM(JI,vMI)(inputn,outputh,seeds0),
F=PXTM(JI,vMF)(inputn,outputh,seeda),
G=PXTM(JG,vMG)(inputd,outputa,seeds),
J=PXG(vMJ)(inputh,outputs,seeda),
Fz=PXTM(JFZ,vMFZ)(inputn,outputz,seedaz),
Gz=PXG(vMGZ)(inputd,outputaz,seedsz),
Jz=PXG(vMJZ)(inputh,outputsz,seedaz).
Pad seeds to integer slots and retain suffixes. Ifk>n or z>h or100000000000z>2886344746477w or
label outside1,...,5,return0^z. Otherwise use restricted8.4: start
W=I(X,prefix_s0 Y);in coordinates1--3 computeT=prefix_s W,A=G(Y,T),
Z=J(W,A),B=G(Y,Z),chooseC=B for bit1,A for bit0,andsetW=F(X,C).
In coordinate4 computeA=Gz(Y,prefix_sz W),Z=Jz(W,A),B=Gz(Y,Z);
chooseC=B for bit1,A for bit0,andsetW=Fz(X,C). Return z-bitW.
The five label codes are1100,1010,1001,0110,0101. NB's4-bit encode(j),
j1,...,5,selects the listedjthword. Putd_F=d,k_F=k.
For z<=h and100000000000z<=2886344746477w,restricted8.4 proves this evaluator for arbitrary
source9-tuple independent of arbitrary seed9-tuple,honestentropy>=k,
honestseed uniform,and honestlabel different from all others. Within-side
correlations are arbitrary;each nonhonest label occurs at most twice.
NB has exactly two copies of each other label on each conditioned branch.
Error76*2^(-E_F-7)<2^-E_F retains ALL8 other outputs and ALL9 original seeds.
Put delta=10^-18. For integerw>=2^80,ceil(Aw)<=(A+delta)w;
each left slot is at least its coefficient timesw.
The coefficients(c,pi,Lam,g)are7g/7h's exact recurrences at each role's
ordered rates,pi=4+sum_i(R_i+1). First the source threshold rows establish
K0<=u;then seed rows apply. With x'=x+delta,the unrounded rows are
aJ=c_J*s'+pi_J,aF=c_F*h'+pi_F,sG=c_G*a'+pi_G,
s0I=c_I*h'+pi_I,hJ=(5+Lam_J)s'+g_J+1,
hJz=(5+Lam_Jz)sz'+g_Jz+1,hcap=100,[h>=s]=s',
azF=c_Fz*(2886344746477/100000000000)+pi_Fz,
azJz=c_Jz*sz'+pi_Jz,szGz=c_Gz*az'+pi_Gz,[sz<=s]=sz',[az<=a]=az',
dEarly=9s0'+(27+Lam_G)a'+g_G+1,
dLast=9s0'+27a'+(9+Lam_Gz)az'+g_Gz+1,[d>=s0]=s0',
kFearly=(22+Lam_F)h'+13s'+g_F+1,
kFlast=18h'+18s'+(8+Lam_Fz)(2886344746477/100000000000)+g_Fz+1,
kI=Lam_I*h'+g_I,[k>=d]=d'.
Exact rational substitution gives every rhs<=its selected left coefficient;
verify.py checks all rows. They imply restricted8.4 withq<=w,including
h>=5sz+KJz+q and the separate last KGz/KFz guards. No comparison
KFz<=KF,KGz<=KG orKJz<=KJ is assumed. The seed rows give both a-seeds,
s,s0,az seeds forFz/Jz and sz seed forGz. Actual100000000000z<=2886344746477w.
All thresholds are paid directly: KJ,KJz<=h,KG,KGz<=d,KI,KF,KFz<=k.
K_I>=h,k>=d andk<=n bound every input byn. SelectedlgB<=8 ande=E_F+7
givee+lgB+6+2lg(4u)<=E_F+21+2lg(4n)<=w. Every output is<=256B*w:
100000000000z<=2886344746477w,B_Fz=12;other output-slot bounds are explicit and smaller.
The terminal505,rate/arity and dimension guards are checked at each new
vector;all source/seed/ambient/field guards follow7g--7h before extraction.

### 4.3 Parameter algorithm
Put Ee=E+4 and the unit w_N=max(2^80,Ee+38+6+2lg(4n)).
1. (d_F,k_F): F# at source length n,output z,advice 4,exponent Ee (4.2);
m1:=d_F.
2. m0,K2: the exact seed length and threshold of Ext2 at input length n,
output m1,exponent Ee.
3. sig=Ee+4;r0=ceil((m0+p0*w_N+D0+(23/4)Ee+56-bc0*a)/(theta-bc0));
r=r0-a; if r<1 set d_C=1,k_C=n+1 and use the constant0^z map,
so the source guard fails without calling Ext0;
s=floor((eta*r0-11sig-8)/4);b=a;
v=6r-floor(eta*r/2).
4. d0,K00: the exact seed length and threshold of Ext0 at input n,output v,
exponent Ee.
5. Rrow=2(r+a+1);e_L=ceil((d0+Ee-1)/2)+Ee;d1=m0+4e_L+lg(4(Rrow+1)).
6. d2=d0+d1+m1+K2+Ee.
7. d_C=d0+d1+d2;k_C=max(K00,k_F+10m0+3+z+2Ee).
The guards are k_C<=n,d2<=n,the primitive guards (PG) below,and the advice guards
(AG1) s>=a+40,(AG2) eta*(r+a)>=(11321180837/12500000)(a+sig+2) (§4.1),public integer tests
on the computed s,r,a,sig. r0 includes the explicit -bc0*a term;proof.md Lemma B9 discharges both.
(PG) Public integer tests require z<=h and100000000000z<=2886344746477w_F for F#,and both Ext2's output m1<=256*BY*w_N
and Ext0's output v<=256*BM*w_N;their K0<=n guards;and each call's
h<=13,terminal j<=18,R_h<=64,kap<=505 and the unit guard of7g.
The selected fixed vectors satisfy the rate/height tests. The parameter
algorithm computes literal lengths before testing PG;it invokes no seed
envelope to decide these tests. On any failure return0^z. NB's theorem
has this expanded guard domain;B9 proves every test for its actual calls.
### 4.4 Evaluator
Slice y=y0y1y2 with lengths d0,d1,d2. Components:
* Ext0: the strong linear primitive of step 4,output truncated to v bits.
* C_adv: Section3.5 on four raw groups of lengths r,r,2r-floor(eta*r/2),2r,accuracy sig=Ee+4,
actual advice length a. It returns five Rrow-bit rows;its selected good row has index<=5.
* LinExt: canonical family7g from Rrow bits to m0 bits,exponent e_L,
with its complete d1=m0+4e_L+lg(4(Rrow+1)) seed bits.
* Ext2: the strong primitive of step 2 applied to y2 zero-padded to n bits,
with seed W_j (m0 bits) and output m1=d_F bits.
* F: F# of §4.2;it reads all of U_j.
V=Ext0(x,y0);(R_1..R_5)=C_adv(V,alpha);
W_j=LinExt(R_j,y1);U_j=Ext2(y2,W_j);Z_j=F(x,U_j,encode(j)),j=1..5;
Q=Z_1+...+Z_5,(NBE)
with encode(j) the 4-bit binary label of j. The same map gives Q' on primed
inputs. The evaluator XORs all five labeled rows;the good index i is analysis only.
Step5 uses the exact canonical7g family at input dimension Rrow,output
m0 and exponent e_L. Its uniform collision estimate and W2 apply to every
conditional source law used below. B9 bounds Rrow,m0,e_L by O(log actual N),
so its generator is polynomial in actual N. Fixed-seed source-linearity and
full seed retention hold by7g;d2 uses this actual d1.
### 4.5 Finite inequalities,each proved here
* (I1) PTM gives d0<=bc0*r+p0*w_N+D0+c0 and2e_L<=d0+3Ee.
The r0 ceiling gives theta*r0>=m0+bc0*r+p0*w_N+D0+(23/4)Ee+56.
Thus s>=m0+bc0*r+p0*w_N+D0+3Ee+42,using floor(x)>=x-1 and sig=Ee+4.
Since c0<1,this is at least
m0+d0+3Ee+41. Thus s>=m0+2e_L;(AG1) gives s>=a+40.
The denominator theta-bc0 is positive.
* (I2) Write u=r+a=r0. The r0 formula and positive raw r give eta*u>=4m0+4p0*w_N+4D0+23Ee+224,so eta*u>=64.
Guard(AG2) gives eta*u>=(11321180837/12500000)(a+sig+2).
Hence(a+sig)/u<=eta/(11321180837/12500000) and
r-sig=u-a-sig>=(1-eta/(11321180837/12500000))(u+1)+1.
Doubling pays Z56;subtracting m<=eta*r/2<=eta*(u+1)/2 pays Z34.
Also a+sig<=eta*u/DG-2 and eta/4+eta/DG<1 imply
r=u-a>=s+sig: indeed s<eta*u/4 and r>=(1-eta/DG)u+sig+2.
Thus the positive-raw-length branch passes and Lemma U's r>=s+sig holds.
The s definition gives eta*u>=4s+11sig+8,which is(U0').
Explicit p0>=16,w_N>=2^80 give eta*u>=4p0*w_N>=64*2^80.
With eta<=1/120 this implies u>=7680*2^80>2^90.
Its conditional entropy kappa>=(5/4)eta*u>=5s+(55/4)sig+10>=5s+11sig
pays Lemma U's matching and conditioning losses;its error is<=14*2^-sig<2^-Ee.
* (I3) Ext2 is legal: K2<=d2<=n;w_N covers the actual binary depths14 and6.
Independently u=r0>=4p0*w_N/eta>=64*2^80/eta>=7680*2^80>2^90,
using p0>=16,eta<=1/120. The strengthened CZ guard eta*(r+a)>=64*2^80
thus holds at every actual call,along with the independent prime guard.
* (I4) d2-d0-d1-m1-Ee=K2 (step 6).
* (I5) k_C-(10m0+3)-z-2Ee>=k_F (definition of k_C).
* (I6) ((d0+Ee-1)/2)-e_L<=-Ee by definition of e_L;Lemma W2 applies.
* (I7) d_F=m1 (step 1).
* (I8) Ext0's output is at least v (its leaf);and K00<=k_C<=Hinf(X)<=n.
* (I9) k_F<=k_C<=n and PG gives z<=h and100000000000z<=2886344746477w_F,so F's guards hold;F sees 8<=t0 tamperings below.
### 4.6 Proof of (NB)
Left is (X,X') with every variable computed from it and from already revealed
right values;right is (Y,Y'). Initially they are independent. Every revelation
below is a function of one side and the history,or a kernel of such values,so
the conditional product persists [1.1]. Write Y=(Y0,Y1,Y2),Y'=(Y0',Y1',Y2').
1. Reveal T0=(Y0,Y0'),right-only. Given T0,V and V'=Ext0(X',Y0') are left
functions. As X is independent of T0,D(V|T0)=E_(Y0) TV(Ext0(X,Y0),U_v)<=2^-Ee
(strong primitive;K00<=Hinf(X) by (I8)). Couple V to Vbar,uniform given T0,
by the common-mass kernel of (V,T0) [1.1];the source law and V' are
untouched. Cost eps0<=2^-Ee.
2. For each T0 apply Lemma U to(Vbar,V'): Vbar is exactly uniform onv bits,
V' is arbitrary,alpha!=alpha',and every finite guard holds by(I2). Couple
the matrix pair to the exact mixture (A,A') with label M by a kernel of
(Vbar,V',T0),cost e_A<=2^-Ee;i=i(M) is the component's good row. Reveal M
and A'_i (left kernels): H1=(T0,M,A'_i). On every branch Hinf(A_i|H1)>=s.
3. LinExt step. Y1 is independent of Y0 and Y0' has d0 bits,so
Havg(Y1|T0)>=d1-d0 and Hinf(Y1|T0=t)>=d1-d0-Ee outside T0-mass 2^-Ee [1.1].
Left and right are independent given T0,so Law(Y1|H1)=Law(Y1|T0) and A_i is
independent of (Y1,Y1') given H1. On good T0,LinExt with a uniform seed has
error 2^-e_L on A_i (entropy s>=m0+2e_L,(I1));with seed Y1 of deficit
d0+Ee,weak-seed LemmaW2 gives TV((W_i,Y1),(U_m0,Y1)|H1)<=2^((d0+Ee-1)/2-e_L)<=2^-Ee
(I6). Reveal (Y1,Y1'),right-only: H2=(H1,Y1,Y1'). W'_i=LinExt(A'_i,Y1') is
fixed. Given H2,W_i is a function of the left variable A_i,and its distance
from uniform is f(Y1)=TV(LinExt(A_i|H1,Y1),U_m0),whose average over H2 is the
H1-average above. So D(W_i|H2)<=2*2^-Ee (bad T0 charged 1). Couple W_i to
Wbar_i,uniform given H2,by a kernel of (W_i,H2);Wbar_i is left-sided given
H2. Cost delta_W<=2*2^-Ee.
4. Ext2 step. U'_i=Ext2(Y2',W'_i) is a right function of m1 bits;reveal it:
H3=(H2,U'_i). Wbar_i stays uniform (right-only revelation). Leakage of Y2:
for h=(y0,y0',mu,a',y1,y1',u'),as (M,A'_i) is independent of the right pair
given T0,
P(Y2=y2,H3=h) <= P(Y0=y0,Y1=y1,Y2=y2)P(M=mu,A'_i=a'|T0=(y0,y0'))
= 2^(-d_C) P(mu,a'|y0,y0').
Summing the maximum over y2 across all h,the (mu,a') sum is 1 for each
(y0,y0'),so the total is at most 2^(d0+d1)2^(d0+d1+m1)2^(-d_C)
=2^(d0+d1+m1-d2). Hence Havg(Y2|H3)>=K2+Ee (I4),and Hinf(Y2|H3=h)>=K2
outside mass 2^-Ee. Put U_i=Ext2(Y2,Wbar_i): Ext2 is a strong (K2,2^-Ee)
extractor (zero-padding Y2 to n bits keeps its min-entropy),and its seed
Wbar_i is uniform and independent of Y2 given H3,so D(U_i|Wbar_i,H3)<=2^-Ee
on good branches. Y2 is independent of every left variable given H3
(X,X',V,V',Vbar,A,A',M,Wbar_i,coins),so Law(U_i|H3,all left)=
Law(U_i|H3,Wbar_i). Couple U_i to Ubar_i,uniform given (H3,all left),by a
kernel of (U_i,H3,Wbar_i);cost delta_U<=2*2^-Ee. This kernel reads only
values present at step 4,so Ubar_i is defined before step 5,and given the
left record it is right-sided.
5. Use a finite left record instead of revealing the raw blocks.
Put Tpast=(T0,Y1,Y1'),Rfuture=(Y2,Y2'),and let Laux be the entire
left-generated record through step3,including all matrix and coupling
variables. Conditional on Tpast,the joint law factors as
P(X,X',Laux|Tpast) P(Rfuture|Tpast): every earlier added kernel reads
only left variables and Tpast. Tpast itself is independent of (X,X').
Put Lrec=(i,(W_j)_(j!=i),(W'_j)_(j=1,...,5),Wbar_i),b_G=10m0+3.
Here W_j=LinExt(A_j,Y1),W'_j=LinExt(A'_j,Y1') are synthetic seeds.
These left variables have at most5*2^(10m0)<=2^b_G values.
The unrepaired W_i stays latent: the honest synthetic seed Ubar_i
uses only Wbar_i in the right repair;restoration still pays delta_W.
Reveal Lrec and U'_i=Ext2(Y2',W'_i),and call G=(Tpast,Lrec,U'_i).
Since W'_i is fixed by Lrec,the last revelation is right-only. Thus
Law(X,X',Laux,Rfuture|G)=Law(X,X',Laux|Tpast,Lrec)
 times Law(Rfuture|Tpast,W'_i,U'_i).
Given H3 and all left variables,this same right law depends on them only
through (Tpast,W'_i,U'_i). Thus the common-mass repair of step4
can and does use a kernel reading only
(U_i,Tpast,W'_i,U'_i,Wbar_i),with U_i=Ext2(Y2,Wbar_i).
For every G this kernel couples its actual U_i marginal to U_m1.
Its conditional marginal is identical to that used in step4,so the
averaged disagreement remains delta_U<=2*2^-Ee,including bad branches.
The whole (Rfuture,repair coins,Ubar_i) is independent of (X,X',Laux)
given G,and Ubar_i is uniform given G. The repair coins are not revealed.
In particular every other U_j,U'_j is right-sided given G. No raw V,V',
mixture label M,A'_i or unbounded auxiliary record is revealed in G.
Also Law(X,X'|G)=Law(X,X'|Tpast,Lrec). Since Tpast is independent
of X and Lrec has at most2^b_G values,[1.1] gives
Havg(X|G)>=k_C-b_G,hence Hinf(X|G)>=k_C-b_G-Ee
outside mass2^-Ee. 6. Reveal Z'_i=F(X',U'_i,encode(i)): U'_i is fixed in G,so this is a z-bit left
revelation. With H5=(G,Z'_i),Hinf(X|H5)>=k_C-b_G-z-2Ee>=k_F (I5) outside total
mass 2*2^-Ee,and Ubar_i stays uniform.
7. F step. Given H5 the source tuple (X for the honest call and the other 4
honest rows;X' for the 4 tampered rows j!=i) is independent of the seed
tuple (Ubar_i;U_j,U'_j for j!=i). The honest seed Ubar_i is uniform on
d_F=m1 bits (I7);encode(i) differs from every other label;the 8 tamperings
are exactly t0=8,with exactly two copies of each nonhonest label. This verifies the selective-transcript contract on every conditioned branch;the tampered row i is fixed in H5 and is
not a tampering. F's direct multicopy contract (Section4.2 and restricted8.4,zero shifts)
makes Zbar_i=F(X,Ubar_i,encode(i)) 2^-Ee-close to uniform jointly with every
other row output and every seed. With Qbar=Zbar_i+sum_(j!=i)Z_j and
Q'=sum_j Z'_j (Z'_i fixed),D(Qbar|Q',seeds,H5)<=2^-Ee.
8. Retain Y,Y'. Given (H5,seed tuple),the outputs depend only on (X,X') and the
seeds,and (X,X') is independent of (Y,Y',coins),since the seeds are
right-sided given H5. So Law(Y,Y'|H5,seeds,outputs)=Law(Y,Y'|H5,seeds) is a
common kernel of the real and ideal laws,and appending (Y,Y'),including
unused seed suffixes,costs nothing [1.1]. Averaging over H5 and charging the
two X tails gives synthetic error <=3*2^-Ee.
9. Restore. Real and synthetic (Q,Q',Y,Y') agree unless V!=Vbar,the matrix pair
differs from (A,A'),W_i!=Wbar_i or U_i!=Ubar_i,so
zeta<=eps0+e_A+delta_W+delta_U<=6*2^-Ee. By [1.1] (D_P<=D_Q+2zeta),
D(Q|Q',Y,Y')<=3*2^-Ee+12*2^-Ee=15*2^(-E-4)<2^-E.
### 4.7 Running time
At actual length N,closed slots and B9 bound every NB/F# dimension
and enumerated field degree by O(logN);the primitive recurrences and
CZ/E6 setup obey the same bound. Each enumeration costs2^O(logN);
larger fields use towers. The fixed number of calls is thus polynomial
in actualN,without a joint bound for unrelated n,z,a,E.
Attribution: Li (§3 and header),advice condenser and row-XOR;
weak-seed second-moment accounting;canonical LinExt7g.
## 5. Scope
Acyclic chain: PXT/PXG7g and PXTM7h -> conditionaltools8.2.2,multicopy8.2--8.4 -> F#
(Section4.2);CZ/E6 -> advice condenser;that condenser,canonical LinExt/W2,
Ext0,Ext2,F# -> NB -> outerA# (proof.md,Lemma B9). Global guards and the
family/circuit assembly are in proof.md.
Credits: dbc05fc4-aed7-487a-b708-d1e5b938fbf4 (U,W2,TR/TL);
b15 (103eta precursor);c78568c1-1aa4-4da9-8556-877fb4218f6f
(whole-tuple MC-LA). Current E6 and four-round selective transcript,
entropy and seed-retention proofs are supplied above and in restricted8.2--8.4.

# Sampler part
# Finite support graph,linear sampler,and conditional breaking
Fix Q_S=42787674682353/20000000000000. With L=floor(log2(actual input length)),S=floor(L/Q_S)
is an integer;S<=L/Q_S<S+1 and L<Q_S(S+1). All other dimensions retain
their specified integer rounding.
All logarithms are binary,addition is XOR unless a field is specified,and
statistical distance is TV = half the L1 distance. A k-source is a finite
random variable with every atom at most 2^(-k);k may be real. Graphs have
labelled edges and may have loops or repeated edges. Every constant is fixed
independently of the eventual outer input length N.
The support-growth Lemma S1' (§SP.6.2),the affine-only
three-condenser/kernel base of linear Section7f and the affine breaker yield the finite
graph/sampler and conditional-breaking reductions below. proof.md supplies
the tagged affine-permutation outer assembly and finite guards.
## 1. Supplied contracts and elementary probability facts
### 1.1 Selected seeded-linear contracts
The supplied general-source maps are PXG(v) from linear7g and
PXTM(j,v) from linear7h. Each fixed ordered rate/arity vector defines
an exact integer threshold K_0(u,z,e),complete seed length D(u,z,e),
and deterministic map Ext:{0,1}^u x {0,1}^D -> {0,1}^z,source-linear
for every fixed seed. On every legal call,for every source X with
H_inf(X)>=K_0 and an independent uniform complete seed Y,
TV((Ext(X,Y),Y),(U_z,Y))<=2^-e. (1)
All ignored suffixes are retained;a shorter output takes a prefix.
The literal algorithms,threshold and seed slices include all ceilings.
Before field setup require K_0<=u. PXTM hasj_i>=2,R_i>j_i-1 and
w>=max(2^20,e+lg(product j_i)+6+2lg(4u)). Section7h proves its
source bound K_0<=Lam*z+g*w. BOTH PXG and PXTM use the clean seed
bound D<=z/B+[4+sum_i(R_i+1)]w only under the additional output-cap,
terminal-rate/arity/kap and height guards of7g. Closed slots and B9
check these hypotheses for every selected call.
Ext0 and Ext2 use their respective PXTM vectors;inner and outer PXG
maps remain binary. Canonical condenser fields are supplied in7a;leaf matrices in7g.
Enumerated base degrees are O(log actual N) at selected calls,so exhaustive base
field setup costs N^O(1);larger relative extensions use polynomial-cost
towers. Linear5/proof.md prove uniform polynomial time in actual N,
without asserting polynomial exact-degree search for arbitrary u.
The selected BASE instead uses linear7f's affine-only three-condenser/
kernel algorithm with its complete seed and integer threshold. That
affine contract alone is consumed in Sections3--5;it is not substituted
for the general-source contracts needed inside the breakers.
### 1.2 All-seeds-retained affine-breaker contract
Let F:{0,1}^n x {0,1}^d x {0,1}^a -> {0,1}^m be one fixed
deterministic map,with a>=1. Its supplied degree-r contract at source
entropy K and error eta means the following universal statement. For
arbitrary correlated tuples
Left=(A,A_1,...,A_r),
Right=(B,B_1,...,B_r,Y,Y_1,...,Y_r),
assume Left and Right are independent,H_inf(A)>=K,Y is uniform on d
bits,and fixed advice strings alpha_i all differ from alpha. Set
V=F(A+B,Y,alpha),V_i=F(A_i+B_i,Y_i,alpha_i).
Then,with U_m independent of all displayed side information,
TV((V,V_1,...,V_r,Y,Y_1,...,Y_r),
(U_m,V_1,...,V_r,Y,Y_1,...,Y_r)) <= eta. (4)
Only A needs entropy and only Y uniformity. No A_i need be a function
of A. All within-side correlations,repetitions,and seed equalities are allowed. Tampered advice
strings need not differ from one another. For a smaller degree,pad with
repeated permissible copies and discard their outputs/seeds. This
branch uses the all-seeds-retained affine contract.
### 1.3 Probability identities
For any common marginal Z,grouping the TV sum by z gives
TV((V,Z),(V',Z))=sum_z P(Z=z) TV(V|z,V'|z). (5)
A common deterministic map or random conditional kernel cannot increase TV:
sum differences over each fiber and use the triangle inequality. In all
ideal laws below a displayed U is fresh uniform independent of the entire
retained side information.
Define average conditional min-entropy by
Havg_inf(A|Z)=-log sum_z P(Z=z) max_x P(A=x|Z=z).
If Z has at most 2^ell possible values,then
sum_z max_x P(A=x,Z=z) <= 2^ell max_x P(A=x),
Havg_inf(A|Z) >= H_inf(A)-ell. (6)
If Havg_inf(A|Z)>=K+log(1/delta),with 0<delta<=1,the total probability
of z with H_inf(A|Z=z)<K is at most delta: on that set the maximum
conditional atom exceeds 2^-K,while its expectation is at most
delta*2^-K. These facts require no unstated independence.
## 2. An explicit small-bias Cayley graph
Fix a field dimension r_F>=1 and put c=2r_F;Section3 chooses r_F=S,s=2ceil((50/7-2^-79)*S/2).
When r_F is a power of two,one admissible choice is the canonical tower field of binary dimension r_F
(§3 does not use it). Here is its complete
finite choice. Start with F_2 and trace-one element a_F=1. Repeatedly
adjoin w satisfying w^2+w+a_F=0 and replace a_F by a_F*w. The image of
x -> x^2+x is exactly the trace-zero subspace: its kernel has size two,
each image has trace zero,and the trace is nonzero because the current
a_F has trace one. Thus the quadratic has no root and is irreducible.
The two conjugates of w are w,w+1,so the relative trace of a_F*w is
a_F,and trace transitivity proves absolute trace one at the next level.
Trace transitivity follows by grouping the Frobenius-power sum. This
constructs every power-of-two dimension with a fixed ordered binary basis,
without an irreducible-polynomial search.
In one tower step write x=x0+x1*w,y=y0+y1*w. With
p=x0*y0,q=x1*y1,r=(x0+x1)(y0+y1),t=a_F*q,their product is
(p+t)+(r+p)*w. Four recursive products and linear extra work give O(r_F^2)
bit operations for multiplication in the final field.
Section3 uses r_F=S,s=2ceil((50/7-2^-79)*S/2) and,for every r_F (also powers of two),the canonical field of
weak-seed §4 of linear-and-bilinear-primitives.md (the first irreducible binary polynomial of degree
r_F in lexicographic order,found in 2^(r_F)poly(r_F) bit operations) with its
polynomial basis;schoolbook multiplication costs O(r_F^2) bit operations.
Everything below uses only field arithmetic and a fixed binary basis.
For each m>=1 and each label (u,w) in F^2,define an m-bit increment
g(u,w)_j=<w,u^j>,j=0,...,m-1,(8)
where brackets are the ordinary coordinate dot product and u^0=1 also
when u=0. On vertices F_2^m use the 2^c labelled edges x -> x+g(u,w).
A fixed label is an involution,so the transition matrix P is symmetric
and doubly stochastic,with the uniform stationary distribution.
For any nonzero a in F_2^m,the character (-1)^<a,x> has eigenvalue
E_(u,w)(-1)^<w,p_a(u)>,where p_a(T)=sum_j a_j*T^j. The polynomial is
nonzero of degree at most m-1. Averaging over w annihilates each nonzero
coordinate vector p_a(u),leaving its fraction of zeros. A nonzero degree-j
polynomial has at most j roots,by repeated division by T-alpha. Hence
0<=eigenvalue<=(m-1)/2^r_F. (9)
The binary characters form an orthogonal basis,so ||P-J||<=lambda=(m-1)/2^r_F,
where J projects onto constants in the uniform-normalized Euclidean inner
product. A labelled step is evaluated by successive field powers and dot
products in O(m*r_F^2) bit operations.
## 3. The affine-base walk
Fix Q_S as fixed above,beta=50/7-2^-79,q_walk=17,C=4(q_walk+1)=72 and
p=144>3(beta+2q_walk),a retained upper budget. For actual n let L=floor(log2 n),
S=floor(L/Q_S)>=2^160,s=2ceil(beta*S/2)<=8S,ellc=lg(s+4),
R(s)=0,
r_F=S and m=s+2q_walk*r_F<=48S+320ellc.
The integers s,m are even. Use Section2's graph on F_2^s with field
degree r_F. An m-bit string specifies an s-bit start and q_walk complete,
independent 2r_F-bit labels. Its positions are W(u,j),j=0,...,q_walk.
The transition matrix is symmetric and ||P-J||<=lambda=(s-1)2^-r_F.
Lemma W_A. If |U|>2^(9m/16-8) and N subset F_2^s has
|N|<=5*2^(-r_F)2^s,some walk position of some u in U lies outside N.
Proof. Put rho_N=|N|/2^s and let Pi project coordinates onto N.
Then ||Pi P Pi||<=rho_N+lambda. A uniform walk stays in N with
probability at most rho_N*(rho_N+lambda)^q_walk. Thus the number
of all-bad m-bit strings is at most
5*2^(m-r_F)*((s+4)2^-r_F)^q_walk.
Its logarithm minus9m/16-8 equals
(7/16)s+log2(5)+8-r_F-q_walk*(r_F/8-log2(s+4))
<=-(7/16)*2^-79*S+100+20log2 S<0.
Here s<=(50/7-2^-79)S+2 and s+4<=16S.
For 2^j<=S<2^(j+1),j>=160,the negative term in magnitude is
at least (7/16)*2^(j-100),while the positive part is at most100+20(j+1).
At j=160 the first exceeds3320;it doubles at each step while the
second increases by20 and is positive. Induction proves strict negativity
at EVERY allowed S. U cannot consist only
of all-bad strings. This proves W_A at every allowed integer S.
## 4. The support graph and property (G_A)
Set s_tag=3m<=144S+960ellc,A_rows=2^s_tag,A0=2^s,
Lambda_S=2^(3m/2-17),epsilon0=2^-r_F. A_rows is the number of rows;
it is unrelated to the source coefficient A_B in linear Section7f.
Use the SAME four AB+C row maps of §SP.6.2 and define
Gamma(alpha,(i,j))=W(row_i(alpha),j),i=1,...,4,j=0,...,17.
The index order is lexicographic in (i,j),giving degree exactly C=72.
Theorem G_A. For every set Q0 of at least Lambda_S left vertices,
|Gamma(Q0,[C])|>5*epsilon0*A0.
Proof. Since m is even Lambda_S is an integer,and
m>=41*2^160>134217728. The complete direct incidence support lemma S1'
gives a union U of row images with |U|>2^(9m/16-8). If its neighborhood
had size <=5*2^(-r_F)2^s,W_A would give a row image and a walk
position outside that neighborhood,a contradiction. Canonical prime selection,field construction,AB+C and
the increments are the supplied algorithms and cost2^O(s_tag).
Proof.md outer Sections2,6 prove s_tag<=2pS<=(2p/Q_S)L,so the
complete row enumeration is polynomial in actual n.
## 5. The evaluated source-linear sampler
Let E be the actual three-condenser/kernel base of linear-and-bilinear-
primitives.md Section7f,output d,seed padded to s bits. On legal calls it
has strong error<=epsilon0=2^-r_F for every uniform affine source of
dimension>=K_B. Define Samp(x,alpha,j)=E(x,Gamma(alpha,j)),using Section4's
actual graph. This map is source-linear at every alpha,j. AF1--AF4 supply
its complete seed,padding,entropy and runtime bounds. The outer proof uses
this strong base error and G_A directly for fixed rank-good neighbors.
## 6. Conditional affine breaking
### 6.1 Conditional/average-entropy lift with arbitrary copies
Suppose(4)holds at degree r,entropy K,error eta. Put
L=(A,A_1,...,A_r),R=(B,B_1,...,B_r,Y,Y_1,...,Y_r).
Assume L--Z--R,(Y,Z)=(U_d,Z),fixed legal advice,and
Havg_inf(A|Z)>=K+log(1/delta). By(6),outside Z-mass at most delta,
Hinf(A|Z=z)>=K. Conditional independence of the full tuples and uniform Y
then meet every hypothesis of(4). Averaging by(5),joint TV retaining Z is
at most eta*(1-P[bad])+P[bad]<=eta+delta. The arbitrary source-copy
quantifiers persist. Taking eta=delta=e_F/2 needs entropy K+log(2/e_F).
This is conditional averaging,not a new evaluated algorithm.
## 8. Selected outer use
Proof.md composes this graph/sampler and conditional breaker with tags,
fixed rank-good neighbors,whole-right profile postselection and its four-moment argument.
Attribution: extractor-to-sampler conversion is due to Zuckerman (1997);
the low-rate-condense/high-rate-walk architecture to Zuckerman,Theory of
Computing 3 (2007),§§2--3;the one AB+C support step to
Barak--Impagliazzo--Wigderson (credited in Part SP);the row-XOR/majority architecture to
Chattopadhyay--Liao,arXiv:2110.12652v1,§5. The affine interface
is compared to Li,ECCC TR23-023 revision 1.
## 0. Atomic proof map and selected construction
Legacy names denote included proofs: linear-and-bilinear-primitives.md
Sections1--5,7a--h,8 cover seeded-linear,weak-seed-hashing and bilinear tools;
restricted-nonmalleability-and-breaking.md Sections1--8 cover prime fields,
encoding,local conditional tools and antichain8.2;wrappers Sections2--4 cover
affine transformation,advice,CZ/E6 and NB,with restricted8.4's F# interface.
Part SP.1--5 proves rings,SP.7 proves M/I,and SP.6.2 proves S1'.
No missing GUV/CGL theorem is a premise. The dependency order is
linear primitives -> local core -> F#;rings/CZ/E6 -> advice condenser;
condenser/primitives/F# -> NB -> affine transformation;
SP.7 -> S1' -> graph/walk -> proof.md's tagged rank-good four-moment assembly.
## Paired read-once zero-period lemma
A proper read-once full-B2 formula on 2q distinct leaves,q>=1,with any
pairing of its leaves,has a nonzero zero-period on some product of q
allowed lines in the paired three-corner domain D. A zero-period v means
F(u+v)=F(u) everywhere on that flat. Flip input bits so every pair forbids
00;this preserves properness and the assertion. A safe axis line fixes
one endpoint to 1;the other allowed line is the diagonal 01,10.
First suppose all gates are AND or XNOR (1+x+y). Both have neutral input
1. Induct on q. If a root child contains m>=1 complete pairs internally,
fix its leaves with partners outside that child to 1. Those fixings use
at most one endpoint of each crossing pair. Neutrality reduces the child
to a formula on its 2m internally paired leaves,where m<q. Induction gives
a zero-period supported only on those leaves. The other root child is
unchanged. Extend the chosen internal lines by the safe crossing lines
and arbitrary allowed lines on the other pairs. This proves this case.
Otherwise every pair crosses the root,and each child has q leaves.
If one child contains an XNOR gate,choose a leaf from each of that gate's
two child subtrees,say u,v. Fix every other leaf of the chosen root child
to 1. Leave u,v free by fixing their partners in the other root child to
1. These are consistent one-endpoint-per-pair fixings. Neutrality reduces
the chosen child to XNOR(u,v). Flipping u and v preserves it and leaves
the other child unchanged,giving the desired period. If neither child
contains XNOR,both are pure ANDs. For an AND root choose one pair's
diagonal;the whole output is identically zero. For an XNOR root choose
all pair diagonals: flipping all q free coordinates swaps the two AND
monomials,preserving their XNOR. This includes q=1 and completes induction.
Now take arbitrary proper B2 gates and suppose no such zero-period exists.
For any product gate,consider a child whose leaf set contains at most
one endpoint of every pair. Fix all those leaves to safe 1. Its output
cannot be a controlling value at that gate's pin: otherwise the entire
other child is killed. Choose a leaf z in that other child,fix its mate
to 1,and complete to a product flat leaving z free. This is consistent
with the earlier fixings. Flipping z would preserve the whole formula,
a contradiction. Thus this child's all-ones output is noncontrolling.
Complement internal signals as needed so every signal takes value 1 at
the global all-ones input,absorbing complements in parent truth tables.
Leaves stay unchanged;complementing the root preserves zero-periods.
Every affine proper gate becomes XNOR. Every product gate whose two child
leaf sets each contain at most one endpoint per pair becomes AND: both
pins have noncontrolling value 1,and its output at 11 is 1.
Choose a minimal subtree V containing a complete pair. Every proper
subtree of V has at most one endpoint per pair. Consequently all gates
within V are AND or XNOR. Fix V's leaves with mates outside V to 1;
neutrality leaves a nonempty formula on its internally paired leaves.
Apply the auxiliary induction. Its period changes only those internal
leaves,leaving the rest of the full formula unchanged. Extend the flat
by the safe crossing lines and arbitrary remaining pair lines. This
contradicts the supposition and proves the lemma. The requirement that
ALL leaves be paired is essential to this proof.
## Tight incidence tree and terminal exclusion
Use the representation and incidence graph of proof.md §7. Suppose
at its terminal state a=k,q=k,M=2k,s=M-1. All M parameters belong to the
output component. The component bound s_0>=M-1 forces s_0=s. There are M+s
variable vertices and s object vertices,with at most 3s incidences.
Connectedness and M=s+1 force exactly 3s incidences and a bipartite tree;
every object has three distinct variable neighbors.
The output is a gate output. Indeed an affine output on M>=2 parameters
has a nonzero kernel direction v. For any direction v the domain D has
a safe a-flat parallel to v: in each pair choose the allowed line parallel
to its nonzero projection (three corners supply all three directions),
or any allowed line if the projection is zero;keep free coordinates.
An affine output would thus have zero derivative on such a flat,
contradicting the separation hypothesis of §7.
Every auxiliary variable has degree at least two. If it had degree one,
its sole incidence would be a row,since it is not the output. Delete
that row and its uniquely supported auxiliary column. Invertibility of
the remaining block follows by expansion on that column. The represented
function,parameters,Phi and D are unchanged,while s drops to M-2,
contradicting the same component bound. Every nonoutput gate variable has
degree at least two (its defining gate and a consumer);dead gates have
already been removed. The output gate variable has degree one. Parameters
have degree at least one. Their total degree is thus at least
M+2(s-1)+1=3s,so equality forces every parameter to have degree one and
every nonoutput gate/auxiliary variable to have degree two.
Root this incidence tree at the output variable. Each object has one
rootward and two outward variables. For a row the rootward variable is
auxiliary: parameters are leaves and gate variables never occur in rows.
Solve the row for that auxiliary,obtaining an XOR or XNOR of its two
outward variables. For a gate the rootward variable is its designated
output: the native directed path through consumers to the final output
exists by dead-gate removal and must be the unique incidence path.
Thus no product gate is inverted. Orienting the tree computes a proper
full-B2 read-once formula from all M parameters,each used exactly once.
It agrees with the represented function (even off D). The paired lemma
supplies a nonzero zero-period on a safe k-flat in D. Injective affine Phi
maps it to a nonzero original direction and a k-flat,contradicting
separation. The alleged tight terminal state is impossible.
# SP. Finite rings and direct incidence support
The following proofs supply E_p(((2-nu)/(14+(h_R+5/2)/D_E))) from restricted-nonmalleability-and-breaking.md
for every prime p >= 2^(2^90),without CZ Lemma 7.3,Bourgain's filling
lemma or an unspecified constant.
Attribution: the ring/stabilizer framework builds on the
self-improving approximate-ring construction of Terence Tao,*The
sum-product phenomenon in arbitrary rings*,Contributions to Discrete
Mathematics 4(2),2009,pp. 59--82,especially Lemma 1.4,Proposition 3.1,
and Theorem 1.2. Official source:
https://cdm.ucalgary.ca/article/view/61994
The sumset tools below are Ruzsa's inequalities and Petridis's minimal-
growth argument. We work only in commutative rings
with identity. No inclusion AB+AC subset A(B+C) is used: that inclusion is
false in general. Scalar distributions and every product of sets below
are distinguished explicitly.
Notation: nA denotes an n-fold sumset,with 0A={0};AA={ab:a,b in A};
aA is multiplication by the single element a. A unit is an invertible
ring element. All sets used in cardinality estimates are finite and
nonempty unless indicated otherwise. Products in a product of fields are coordinatewise. The selected ring below
is Kp x Kp,using QUADRATIC-FIELD multiplication inside each coordinate;
its four F_p basis coordinates do not multiply coordinatewise.
## SP.2. Dyadic overlaps and weighted unit actions
Let C0 be M units in a commutative ring,K>=2,with
|C0-C0|,|C0/C0|<=KM. Put H=1+ceil(log2(2K)). Assume either
Z=|(C0-C0) outside units|<M/(2^20*K^14*H^2),
or the ring is a product of fields,each coordinate fiber of C0 has size
at most F_j,and W=sum_j F_j<M/(2^20*K^14*H^2). (SP-Mi0)
Then some b0 in C0,B0 subset C0 and finite unital ring T satisfy
|B0|>=M/(4KH),B0 subset b0T,|T|<=2^17*K^12*M. (SP-Mi1)
Product energy equals quotient energy: a/b=c/d iff ad=bc. Cauchy on C0/C0 gives sum_(a,b)|aC0 intersect bC0|>=M^3/K.
Normalize intersections by M. Values below1/(2K) contribute at most
M^2/(2K). Partition the rest into H ranges[mu,2mu),mu=2^-j,
j=0,...,ceil(log2(2K));the top includes1. One range contributes
at least M^2/(2KH),hence has at least M^2/(4KHmu) pairs. Choose b0 with
at least M/(4KHmu) partners B0,adding b0 if absent. Normalize
C=b0^-1 C0,G=b0^-1 B0. Put kappa=4KHmu,L=K^2/mu. Thus
1 in G subset C,|G|>=M/kappa,1/(4K)<=mu<=1,
|gC intersect C|>=mu M for every g in G.
Triangle with middle W=gC intersect C gives |gC-C|<=|gC-W||W-C|/|W|<=K^2*M/mu=LM. (SP-Mi2)
Let S={x:|xC-C|<=2048K^2L^4 M}. Cauchy counts at least
|G|^4/(2048K^2L^4 M) ordered quadruples satisfying xg1-g2=xg1'-g2'.
Put d=g1-g1',so xd=g2-g2'. For each nonunit d at most |G|^2
quadruples occur,giving bad count<=Z|G|^2 in the first case.
In a product of fields,a nonunit d has some zero coordinate j,hence
g1_j=g1prime_j. There are at most |G|F_j such ordered pairs.
For each pair choose g2prime in G;the collision equation forces the
ENTIRE element g2=x*(g1-g1prime)+g2prime uniquely. Its membership in G
can only reduce the count. Thus this coordinate contributes at most
|G|^2 F_j collisions,and the bad count is at most |G|^2 W,where
W=sum_j F_j. This also holds for nonunit x. Unit normalization and
taking subsets preserve the fiber caps.
In either case write V=Z or W. Since
4096kappa^2 K^2L^4=2^16*K^12*H^2/mu^2<=2^20*K^14*H^2 and |G|>=M/kappa,
V<|G|^2/(4096K^2L^4 M).
At least half the quadruples have unit d. Choose one and write
d=g1-g1',xd=g2-g2'. For D=union_{g in {g1,g1',g2,g2'}}gC,
|D-C|<=4LM;(SP-3) gives |2D-2D|<=(4L)^4 M.
Both xdC-dC and xdC+dC lie in this set. Unit dilation yields
|xC-C|,|xC+C|<=256L^4 M. (SP-Mi3)
Repeat Cauchy at this sharper bound. The bad count bound remains
|G|^2 V<|G|^4/(4096K^2L^4M)<=|G|^4/(512L^4M),since K>=2.
Hence each x in S supplies at least |G|^4/(512L^4M) unit-d quadruples.
Every such quadruple determines at most one x. There are at most |G|^4
quadruples in total,proving S finite and
|S|<=512L^4M<=2^17*K^12*M.
For every finite nonempty A the extension estimate is
|A+C-C|<=K^2|A+C|. (SP-Mi4)
Choose nonempty P subset C minimizing k=|P-C|/|P|<=K. The Petridis
extension in SP.1 gives |A+P-C|<=k|A+P|. Triangle with
X=A-C,Y=-P,Z=-C gives |A+C-C|<=k^2|A+P|<=K^2|A+C|.
Thus (x+1)C-C subset xC+C-C and the plus bound in(SP-Mi3),K^2<=L imply translation
by1 preserves S. For W=gC intersect C,triangle gives
|gxC-C|<=|gxC-W||W-C|/|W|<=256(K/mu)L^4M<=1024K^2L^4M,
since |gxC-W|<=|gxC-gC| and |W-C|<=KM. Thus every g in G preserves S. These bijections permute finite S;their
inverses preserve S. Conjugation gives translation by every G-monomial
and its negative. Since0 belongs to S,the additive group T generated by
these monomials lies in S. It contains1,G and is closed under products.
B0=b0G subset b0T proves (SP-Mi1).
## SP.5. Quadratic-field expansion by the minus stabilizer
Under SP.4,for p>=2^(2^90),
max(|A-A|,|C0/C0|)>p^gamma|A|. (SP-19m)
Assume failure. Set N=|A|,M=|C0|>=N/2,K=2p^gamma and
H=1+ceil(log2(2K)). Write B=log2p>=2^90. Since gamma<1/7,
H<=4+gamma*B<=B. By(c0),W=sum_j F_j<=8p^((5/2)eta0).
Using M>=N/2>=p^(2-h_R*eta0)/2 gives
2^20*K^14*H^2*W/M<=2^41*H^2*p^-nu<1. (SP-21m)
Indeed14gamma+(h_R+5/2)eta0=2-nu. On2^j<=B<2^(j+1),j>=90,
41+2log2B<43+2j,while nu*B>=2^(j-80). Base223<1024 and
doubling dominate slope2. SP.2's fiber case thus gives B0,T with
|B0|>=M/(4KH)>=N/(16H*p^gamma),
|T|<=2^17*K^12*M<=2^29*p^(12gamma)*N<p^4.
Here N<=p^2,12gamma<12/7 and29/2^90<2/7 prove the last inequality.
B0 qualifies for(c1),(d),since4+log2H<nu*B anddelta_R=gamma+nu.
Unit scaling preserves projection cardinalities.
The unital subring T is an F_p vector subspace: addition has characteristic
p,so scalar multiplication is repeated addition. Its two projections
are F_p subspaces of Kp and contain more than p elements by(c1),hence both
are all Kp. We now classify this SUBDIRECT product without a black box.
Let J={b:(0,b) in T}. It is an ideal of Kp: for every c in Kp the second
surjection supplies(a,c) in T,so(0,b)(a,c)=(0,bc) is in T.
If J has a nonzero element,J=Kp. Subtracting(0,c) from each(a,c) shows
T=Kp x Kp,contrary to |T|<p^4. Thus J=0;the first projection is
bijective and T is the graph of a unital F_p-field automorphism phi of Kp.
The canonical quadratic-field proof gives phi=identity or conjugation.
Writing b0=(b1,b2),its unit scaling is consequently a graph
(z,c*z) or(z,c*conj(z)),c=b2/phi(b1). Yet B0 qualifies for(d) and lies
in that graph. This contradiction proves(SP-19m),with all subset quantifiers.
The selected family uses gamma=((2-nu)/(14+(h_R+5/2)/D_E)),nu=2^-80,
eta=((2-nu)/(14D_E+h_R+5/2)),theta=eta/4=((2-nu)/(4*(14D_E+h_R+5/2))).
Its cap eta<=1/120 holds by cross multiplication. CZ uses unit deficit
11346180837/11321180837+nu,relative energy exponent192560074229/11321180837+4nu and hereditary cutoff
p^(2-h_R*eta). The finite guard eta*logp>=64/nu pays2^64.
Every actual NB call has eta*(r+a)>=4*p0*w_N>=64*2^80.
Since eta<=1/120,r+a>=7680*2^80>2^90,independently paying the ring guard.
The image cutoff,unit deletion,clipping and fiber partition costs are bounded explicitly in restricted Sections6--7.
The supplied conditional proof gives45p^-eta,the matched marginal90p^-eta,
and the four-block advice lemma14*2^-sig,with all actual algorithm guards.
# SP.6. Direct point-line incidence support growth
§SP.7's full incidence proof gives the AB+C four-row support bound.
The evaluator uses fixed integer arithmetic,canonical primes,AB+C and
bit truncation;no scalar hitting/peeling lemma is used.
## SP.6.0. One support-growth step
### SP.6.1 Exact finite endpoint
Require m>=134217728;this conservative guard implies m>=4. The direct
incidence count below proves the finite endpoint |U|>2^(9m/16-8),keeping
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
Q subset{0,1}^(3m) with|Q|>=2^(3m/2-17),the unionU of itsfourrowimages
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
The complete incidence theorem(SP-I) below thus applies. Set
v=floor(2N^2/(d+2))<=2N^(2/3).
As d-1<N^(4/3) andd<=N^(4/3)+1,it gives
|Q|<=I(P,L)<=2N^2(v+1)+d(d-1)+N^2
<5N^(8/3)+3N^2+N^(4/3)<=9N^(8/3)<16N^(8/3).
The middle inequality usesN>=1. Taking logarithms with the lower bound
onQ gives log2N>9m/16-63/8>9m/16-8,a contradiction.
All counts concern arbitraryQ and its actual row images;no independence
is asserted. Complete polynomial/multiplicity theorem (SP-I) is proved below.
At the selected call s=2ceil((50/7-2^-79)*S/2),r_F=S,
m=s+34S>=41S>=41*2^160>134217728.
