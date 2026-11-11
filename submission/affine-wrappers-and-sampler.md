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
### 2.1 Selected affine-left,structured-right contract and algorithm
Fix n,d,t,a,m>=1,eps in (0,1),q=ceil(log(1/eps)),H=ceil(log t)>=1.
Supply these deterministic maps,each of error <=eps:
C: {0,1}^d x {0,1}^d0 x {0,1}^a -> {0,1}^dx;
L0_shared: {0,1}^n x {0,1}^s0 -> {0,1}^d0,source dimension k;
K0 is its final-stage threshold under a rank-good prefix;
Lr: {0,1}^U1 x {0,1}^dx -> {0,1}^r,threshold Kr;
E: {0,1}^d x {0,1}^dy -> {0,1}^dx,threshold KE;
Lm: {0,1}^r x {0,1}^dx -> {0,1}^dy,threshold Km.
The two-stage kernel maps in7f-local are strong on uniform affine sources and source-linear at fixed seed;E is strong on arbitrary sources.
C is the strong two-sided-affine one-tampering breaker of threshold KC: a source pair jointly affine in one uniform affine variable,independent of a seed pair jointly affine in an independent uniform affine variable,honest
seed uniform,distinct advice,both complete seeds retained. Its first input
really has length d. Padded seeds may ignore suffixes while retaining them.
All actual primitive guards must hold. Put
TR0=(t+1)(s0+d0+dx),TL0=(t+1)d0,T0=TR0+TL0,
TR=TR0+(H-1)(t+1)(2dy+2dx)+(t+1)dy,TL=TL0+2(H-1)(t+1)dy.
Terminal history length is T0+(H-1)(t+1)(4dy+2dx)+(t+1)(3dy+dx).
Require positive legal dimensions,k<=n,and
s0<=d,dy<=r,m<=dy,
r>=(t+1)dy+Km,
k>=max(K0,TL+Kr+t*r),
d>=(t+1)(s0+d0)+KC+q+1,
d>=TR+KE+t*dx+q+1. (3)
The evaluator on (x,y,alpha) is
s=Prefix_s0(y);(xstar,r1)=L0_shared(x,s) as7f-local;
s2=C(y,r1,alpha);w0=Lr(xstar,s2).
For i=1,...,H-1:
p=Prefix_dy(w_(i-1));qm=E(y,p);v=Lm(w_(i-1),qm);
qr=E(y,v);w_i=Lr(xstar,qr).
Finally p=Prefix_dy(w_(H-1));qm=E(y,p);v=Lm(w_(H-1),qm).
Return Prefix_m(v). No final qr or Lr is computed.
The outer A# of proof.md outer Section3 has m=1<=dy,so no long-output
extension is used here. Section4.2 constructs F# separately from restricted8.6.
The selected tuple is A_i=T_i Z+c_i for uniform affine Z,invertible T_i,
independent of an original right tuple(B,Y) jointly affine in one uniform affine variable U,with honest Y0 uniform and
distinct honest advice. Section2.5 proves its complete retained-right
contract. This is an affine-left contract.
### 2.2 Initialization
Write original L=Z,R=(B,Y) tuple. For each honest source direction choose a k-subspace only in analysis. Let E0 be rank preservation by the honest shared C0 on that subspace. It depends only on the honest b_pre-bit prefix of Y0;its complement has mass beta<=eps/4. It is never tested by the evaluator. Work first under normalized E0. Original Z remains uniform affine and independent of the right tuple;H_inf(Y0|E0)>=d-1 since Pr(E0)>=3/4.
Fix only the honest C0 prefix. Its compressed honest affine image has rank>=k,and the remaining honest initial kernel seed bits are uniform independent of Z. The two-stage strong affine L0 gives error<=eps retaining the complete honest S0,uniformly over good prefixes. Now append all other S_j by their actual right conditional kernel given S0,independent of Z. Do not assert honest kernel-seed uniformity given all S_j.
Given S,each compressed coordinate splits as Astar_j=M_(S_j)T_jZ+cstar_j and Bstar_j=M_(S_j)B_j. These maps may be noninvertible and prefix-dependent. R1A is affine in Z and R1B is right-only. Reveal S then all R1B;call this H0,length <=(t+1)(s0+d0). Appending R1B by its right kernel and translating gives D(R1_0|H0)<=eps. Given H0 the entire R1 tuple is left-sided and independent of Y. The stored Xstar values are latent,never revealed or charged as messages.
The original right tuple is affine in U. On each fixed S,the rank-good event E0 is determined and every R1B map is affine in U. Hence U|H0 is uniform affine on every nonempty fiber. Honest Y0 started with rank d and loses at most(t+1)(s0+d0),so its rank is at least KC+q by(3). It is independent of the left seed tuple R1 given H0. For each j,the pair(R1_0,R1_j) is jointly affine in Z on H0. Apply the affine-section construction of restricted8.2.4 to this pair: replace its honest marginal by full uniform,extend the conditional other-coordinate kernel affinely off its support,and keep the source Y tuple independently. The repaired seed pair is uniform affine and joint distance is exactly the honest seed distance. The two-sided NB theorem applies on each H0 fiber since its honest Y rank is already KC+q. Deterministic evaluation and conditional-uniform TV transfer cost at most twice that distance. Averaging gives D(S2_0|S2_j,H0,R1_0,R1_j)<=3eps. This is a bound for the REAL law after transfer,not a permanent change of Z or its source equations.
Other R1 coordinates can be appended because,given H0 and the retained pair,
they are independent of the full Y tuple and hence both C outputs. Reveal
all R1,a left revelation since R1B was fixed. Then S2 is right-sided.
Apply the exact-affine lift to Lr with fixed-prefix sources Astar,seeds S2,I={j},J empty. The k bound in
(3) supplies Kr+r after the history;the exact-affine lift of2.5 gives the conservative error5eps. Let Z0 contain
S,R1B,R1,S2 across all copies,of length <=T0. The full sides remain
independent and
D(W0A_0|W0A_j,Z0)<=5eps.
### 2.3 Histories and causal ledger
After complete round i<=H-1,set Zi by appending(PB,P,Qm,VB,V,Qr) across all copies in
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
(2(t+1)dy bits). Before the second E of a complete round i<=H-1,the current Qr is not yet revealed. Its right cost is TR0+i(t+1)(2dy+2dx)-(t+1)dx. The terminal first E has cost TR=TR0+(H-1)(t+1)(2dy+2dx)+(t+1)dy and dominates all earlier E calls. The normalized good right law starts at entropy d-1;the extra1 in both guards(3) pays this once. C is called before R1,S2 are revealed,at right cost
(t+1)(s0+d0). Every Lr call ends by round H-1 and has left cost <=TL. The terminal P disclosure occurs later and is paid only by the Lm rank guard. These are exactly the
right source constraints(3),with q;the affine left guard has no tail. Lm uses the exact affine source guard of2.5. The terminal history is (Z_(H-1),PB,P,Qm,VB) across all copies. PB,Qm,VB are right messages,P is left. Its honest V is the target,never a revealed message. This rectangular history retains every original right seed. The outer assembly uses this ledger;F# uses
restricted8.6's ledger.
### 2.4 Actual primitive use
Selected L0/Lr/Lm use linear7f-local with the supplied outer length context;all padded suffixes are retained. Only E and LinExt use ordinary7g/7h;the NB/F# maps use7f-AK. The outer assembly uses Section2.1 withC=NB.
The selected innerF# uses the two-sided-affine antichain variant8.6 with five
labels and8tamperings,instantiated in4.2. Its zero-shift independent-
side joint contract is exactly the one used by NB4.6 step7. The full
conditional-average lift remains sampler6.1 and all original source/seed
quantifiers and retained padding suffixes remain explicit.
### 2.5 Affine-left rank specialization
Start with A_i=T_i Z+c_i,uniform affine Z,invertible T_i,independent of an original right tuple(B,Y) jointly affine in one uniform affine U. Use the shared-prefix evaluator2.1 and guards
r>=(t+1)dy+Km, k>=max(K0,TL+Kr+t*r). (3aff)
Use both right guards with q+1. On normalized E0 the terminal error is at most(9*3^H-8)eps. Restoring the original law adds beta. At H=9 its multiplier is at most177139+1/4<2^18.

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
in original Z,by source-linearity with earlier right messages and all shared prefixes fixed. Conditional on those prefixes E0 is right-only,so rectangular independence holds. The honest Astar has rank>=k;each b affine constraints on Z reduce the rank of this image by at most b,even for noninvertible prefix-dependent coordinate maps. The remaining
left variable is uniform on an affine fiber,of codimension at most the
left bits revealed. Right messages preserve this law. Thus every Lr
call on its stored Astar,including at most t tampered r-bit outputs,has the exact entropy
required by(3aff);initial L0 needs only k>=K0. Histories and own-side
ledgers are those of2.2--2.3.

For i<=H-1 maintain D(WiA_0|WiA_J,Zi)<=delta_i for |J|<=min(2^i,t).
Initialization2.2 using the exact-affine Lr lift keeps the conservative
bound delta_0=5eps. At any fixed history z and set J,the rank of WA_0
conditional on WA_J is constant over nonempty WA_J-value fibers.
If it is full r,WA_0 is uniform independent of WA_J. Otherwise every
conditional law is uniform on a proper affine subspace and has distance
at least1/2 from U_r. Hence the rank-bad HISTORY event has probability
at most2delta for a tuple at distance delta. It is a function of old
history,not of a sampled value of WA_J.

For a round choose J=J1 union J2 disjoint,|J1|<=|J2|<=u_(i-1),delta=delta_(i-1).
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
for i<=H-1 the ordinary right E adds2eps and exact-affine Lr adds2eps.
Thus delta_i=3delta_(i-1)+8eps=(9*3^i-4)eps. All averages use the
actual history marginal and fixed-history rectangular independence.
At terminal step H use the SAME partition/rank argument from Z_(H-1),but stop after VA. Its bound is3delta_(H-1)+4eps=(9*3^H-8)eps. For t271,partition the full target set into135 and136;both are covered by the previous256-copy invariant. The retained terminal history is (Z_(H-1),PB,P,Qm,VB),not the honest V. Given it,all VA coordinates are left-only and the full original right tuple is independent of them. Append that tuple by its conditional right kernel in both real/ideal laws,translate every VB and take m<=dy prefixes. This retains all tampered target outputs,every full seed and all original right information. E0 is a function of that retained right variable;mix its good-law distance with at most beta from its complement once. No history error is renormalized.
The prefix-conditioned affine invariant proves this generalized wrapper directly;no globally independent compressed-tuple hypothesis is invoked. Transcript lengths and distinct-advice requirements are unchanged.
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
(a+sig)/u<=eta/(4697/100),r-sig>=(1-eta/(4697/100))(u+1)+1,
eta*u>=64*2^80,u>=2^90. (11')
The scalar guard also gives2r-sig>=2(r-sig)>=(2-2eta/(4697/100))(u+1)+2,
This pays Z56. Subtracting m<=eta*(u+1)/2 pays Z34
with threshold(2-(5097/9394)eta)(u+1)+1,as restricted Section8 proves. Distinct advices make the two
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
## 4. The two-sided affine-source breaker NB
NB retains the advice-condenser,LinExt/W2 and row-XOR algorithms. Ext0/Ext2/I/F/Fz/G/Gz/J/Jz are affine-source AK maps;their rank proof is below. All calls receive original outer context(n_outer,L,Om). Here F=F#,lg is the ceiling binary logarithm,and[1.1]is Section1.1.
### 4.1 Statement
Put DG=4697/100,eta=(2-nu)/(14D_E+h_R+5/2),theta=eta/4,t0=8. Ext0 is AK_X of7f-AK;Ext2 is AK_Ext2 of7f-AK. Ext0 has c0=0,pi0=p0=2gamma_X+eps_N,bc0=0,D0=0. It uses d0=ceil(p0*w_N) padded seed bits,which contain its literal two kernel seeds under PG;thus d0<p0*w_N+1. It does not use an output-dependent seed term.

For positive public parameters the algorithm below computes d_C,k_C and NB:{0,1}^n x{0,1}^d_C x{0,1}^a ->{0,1}^z. Suppose its public k_C<=n,d2<=n,primitive/context/mass guards PG and advice guards AG1/AG2 hold. Let (X,X') be fixed affine images of one uniform affine Z,independent of a pair(Y,Y') jointly affine in an independent uniform affine variable,with rank(X)>=k_C and Y uniform. For distinct advice alpha,alpha',
D(NB(X,Y,alpha)|NB(X',Y',alpha'),Y,Y')<=15*2^(-E-4)<2^-E. (NB)
All actual source and seed quantifiers are part of this statement. No theorem for general-source X or arbitrary nonlinear source conditioning is asserted. Failed public guards return0^z before affected setup.
### 4.2 Four-round antichain inner breaker F#
For source lengthn,requested outputz>=1,advice among NB labels1,...,5,
and exponentE_F>=1,compute the public INTEGER
w=max(2^80,E_F+2lg(4n)+32,ceil(z/128)),e=E_F+7,q=e,t=8,c=9. Writew_F=w.
Every actual slot is ceil(A*w),with these fixed rational coefficients:
Use proof.md's ten fixed slot coefficients with this unit w.
Write h1=hc. At exponent e,use
I_h1=AK_I(inputn,outputh1,seeds0),I_s=AK_Is(inputn,outputs,seeds0s),
F_l=AK_F(inputn,outputl,seeda),
G=AK_G(inputd,outputa,seeds),
J_l=AK_J(inputl,outputs,seeda),
Fz=AK_Fz(inputn,outputz,seedaz),
Gz=AK_Gz(inputd,outputaz,seedsz),
Jz=AK_Jz(inputh2,outputsz,seedaz).
Pad seeds to integer slots and retain suffixes. Ifk>n or z>h2 or62500000000z>1771192802333w or
label outside1,...,5,return0^z. Otherwise use restricted8.6. For the selected public code b1...b4 set l_j=s for bit0,h1 for the FIRST bit1,h2 for the SECOND bit1. Start W=I_(l_1)(X,prefix_(s0 if b1=1 else s0s) Y). In early coordinate j compute T=prefix_s W,A=G(Y,T),Z=J_(max(l_j,h2))(pad_(max(l_j,h2)) W,A),B=G(Y,Z),select C=B for1,A for0,and set W=F_(l_(j+1))(X,C). In coordinate4 compute A=Gz(Y,prefix_sz W),Z=Jz(pad_h2 W,A),B=Gz(Y,Z),select C likewise,and set W=Fz(X,C). Return z-bit W. pad appends zeros. I_h1/I_s have separate gamma/seeds;F_l share gamma and padding. Thresholds are KI_l/KF_l.
The five label codes are1100,1010,0110,0101,0011. NB's4-bit encode(j),
j1,...,5,selects the listedjthword. Putd_F=d,k_F=k.
For z<=h2 and62500000000z<=1771192802333w,restricted8.6 proves this evaluator for a source9-tuple jointly affine in one uniform affine variable,independent of a seed9-tuple jointly affine in an independent uniform affine variable,honest rank>=k,
honestseed uniform,and honestlabel different from all others. Both whole tuples must have their stated common-affine representations;they are mutually independent. Each nonhonest label occurs at most twice.
NB has exactly two copies of each other label on each conditioned branch.
Error76*2^(-E_F-7)<2^-E_F retains ALL8 other outputs and ALL9 original seeds.
Put delta=10^-18. For integerw>=2^80,ceil(Aw)<=(A+delta)w;
each left slot is at least its coefficient timesw.
For I/Is/F/Fz/G/Gz/J/Jz the coefficients(c,pi,Lam,g) are the affine AK envelopes c=0,pi=2gamma+eps_N,Lam=1/G,g=B+eps_N. First the source threshold rows establish
K0<=u;then seed rows apply. With x'=x+delta,the unrounded rows are
aJ=c_J*s'+pi_J,aF=pi_F,sG=c_G*a'+pi_G,
s0I=pi_I,s0Is=pi_Is,[s0>=s0s]=s0s',
hJ1=(9+Lam_J)s'+g_J+1,hJ2=(5+Lam_J)s'+g_J+1,
hJz=(5+Lam_Jz)sz'+g_Jz+1,[h1>=h2]=h2',hcap=100,[h2>=s]=s',
azF=c_Fz*(1771192802333/62500000000)+pi_Fz,
azJz=pi_Jz,szGz=pi_Gz,[sz<=s]=sz',[az<=a]=az',
dEarly=4s0'+5s0s'+(40+Lam_G)a'+g_G+1,
dLast=4s0'+5s0s'+42a'+(11+Lam_Gz)az'+g_Gz+1,[d>=s0]=s0',
kI=Lam_I*h1'+g_I,kIs=Lam_Is*s'+g_Is,[k>=d]=d'.
Pay hJ1 on h1,hJ2/hJz/hcap on h2,and order rows on their named left slots. Additionally for each row(A,B,C,t) of restricted8.6,pay on k the value A*h1'+B*h2'+C*s'+Lam_F*l_t'+g_F+1 for t=1,2,3,where l=(h1,h2,s);for t=4 use A*h1'+B*h2'+C*s'+(8+Lam_Fz)zcap+g_Fz+1. These20rows complete the43row witness.
Exact rational substitution gives every rhs<=its selected left coefficient;
verify.py checks all rows. They imply restricted8.6 withq<=w,including
h2>=5sz+KJz+q and the separate last KGz/KFz guards. No comparison
KFz<=KF_h,KGz<=KG orKJz<=KJ is assumed. The seed rows give both a-seeds,
s,s0,az seeds forFz/Jz and sz seed forGz. Actual62500000000z<=1771192802333w.
All thresholds are paid directly: KJ,KJz<=h2<=h1,KG,KGz<=d,KI_h1,KI_s,KF_h1,KF_h2,KF_s,KFz<=k.
KI_h1>=h1,k>=d andk<=n bound every input byn. Every selected inner map uses7f-AK with the original outer context,rho=e+225+2lg Om<=w+200 eventually,and its separately paid padded seed/source/ambient guards. Right source ranks use8.6's extra affine history disclosures. B9 proves all public guards before extraction.
### 4.3 Parameter algorithm
Put Ee=E+4 and the unit w_N=max(2^80,Ee+38+6+2lg(4n)).
1. (d_F,k_F): F# at source length n,output z,advice 4,exponent Ee (4.2);
m1:=d_F;save its inner slots as s_F,sz_F (distinct from the raw advice entropy s below).
2. m0=ceil((2gamma_Ext2+eps_N)w_N) is the full padded Ext2 seed length;K2 is its literal threshold at input n,output m1,exponent Ee.
3. sig=Ee+4;r0=ceil((m0+p0*w_N+D0+(23/4)Ee+56-bc0*a)/(theta-bc0));
r=r0-a; if r<1 set d_C=1,k_C=n+1 and use the constant0^z map,
so the source guard fails without calling Ext0;
s=floor((eta*r0-11sig-8)/4);b=a;
v=6r-floor(eta*r/2).
4. d0=ceil(p0*w_N);K00 is AK_X's literal source threshold at input n,output v,exponent Ee and original outer context. Retain the full padded d0-bit seed.
5. Rrow=2(r+a+1);e_L=ceil((d0+Ee-1)/2)+Ee;d1=m0+4e_L+lg(4(Rrow+1)).
6. d2=d0+d1+m1+K2+Ee.
7. d_C=d0+d1+d2;k_C=max(K00,k_F+2v+5s_F+sz_F+z+2Ee+100).
The guards are k_C<=n,d2<=n,the primitive guards (PG) below,and the advice guards
(AG1) s>=a+40,(AG2) eta*(r+a)>=(4697/100)(a+sig+2) (§4.1),public integer tests
on the computed s,r,a,sig. r0 includes the explicit -bc0*a term;proof.md Lemma B9 discharges both.
(PG) Public integer tests require z<=h2,62500000000z<=1771192802333w_F,d2<=n and every literal source threshold below its actual input. Ext0/Ext2/I/F/Fz/G/Gz/J/Jz test7f-AK's ambient,threshold,field and padded-seed guards with original context. Also require m0+p0*w_N+(23/4)Ee>=16w_N. Parameters are computed before testing these conditions;no extraction guarantee is used to decide them. B9 proves every such public test eventually for the selected all-length family. On failure return0^z.
### 4.4 Evaluator
Slice y=y0y1y2 with lengths d0,d1,d2. Components:
* Ext0: the strong affine-source kernel map of step4,with full padded seed and output v.
* C_adv: Section3.5 on four raw groups of lengths r,r,2r-floor(eta*r/2),2r,accuracy sig=Ee+4,
actual advice length a. It returns five Rrow-bit rows;its selected good row has index<=5.
* LinExt: canonical family7g from Rrow bits to m0 bits,exponent e_L,
with its complete d1=m0+4e_L+lg(4(Rrow+1)) seed bits.
* Ext2: the strong affine AK_Ext2 primitive of step 2 applied to y2 zero-padded to n bits,
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
* (I1) AK padding gives d0<p0*w_N+1 and2e_L<=d0+3Ee.
The r0 ceiling gives theta*r0>=m0+bc0*r+p0*w_N+D0+(23/4)Ee+56.
Thus s>=m0+bc0*r+p0*w_N+D0+3Ee+42,using floor(x)>=x-1 and sig=Ee+4.
Since d0<p0*w_N+1,this is at least
m0+d0+3Ee+41. Thus s>=m0+2e_L;(AG1) gives s>=a+40.
The denominator theta-bc0 is positive.
* (I2) Write u=r+a=r0. The r0 formula and positive raw r give eta*u>=4m0+4p0*w_N+4D0+23Ee+224,so eta*u>=64.
Guard(AG2) gives eta*u>=(4697/100)(a+sig+2).
Hence(a+sig)/u<=eta/(4697/100) and
r-sig=u-a-sig>=(1-eta/(4697/100))(u+1)+1.
Doubling pays Z56;subtracting m<=eta*r/2<=eta*(u+1)/2 pays Z34.
Also a+sig<=eta*u/DG-2 and eta/4+eta/DG<1 imply
r=u-a>=s+sig: indeed s<eta*u/4 and r>=(1-eta/DG)u+sig+2.
Thus the positive-raw-length branch passes and Lemma U's r>=s+sig holds.
The s definition gives eta*u>=4s+11sig+8,which is(U0').
The explicit mass guard m0+p0*w_N+(23/4)Ee>=16w_N and w_N>=2^80 give eta*u>=4(m0+p0*w_N+(23/4)Ee)>=64*w_N>=64*2^80.
With eta<=1/120 this implies u>=7680*2^80>2^90.
Its conditional entropy kappa>=(5/4)eta*u>=5s+(55/4)sig+10>=5s+11sig
pays Lemma U's matching and conditioning losses;its error is<=14*2^-sig<2^-Ee.
* (I3) Ext2 is legal: K2<=d2<=n;Ext2's affine AK guards are paid.
Independently u=r0>=4(m0+p0*w_N+(23/4)Ee)/eta>=64*2^80/eta>=7680*2^80>2^90,
using the mass guard and eta<=1/120. The strengthened CZ guard eta*(r+a)>=64*2^80
thus holds at every actual call,along with the independent prime guard.
* (I4) d2-d0-d1-m1-Ee=K2 (step 6).
* (I5) k_C-2v-5s_F-sz_F-z-2Ee>=k_F+100 (definition of k_C).
* (I6) ((d0+Ee-1)/2)-e_L<=-Ee by definition of e_L;Lemma W2 applies.
* (I7) d_F=m1 (step 1).
* (I8) Ext0's requested AK output is v;and K00<=k_C<=Hinf(X)<=n.
* (I9) k_F<=k_C<=n and PG gives z<=h2 and62500000000z<=1771192802333w_F,so F's guards hold;F sees 8<=t0 tamperings below.
### 4.6 Proof of (NB)
Left is (X,X') with every variable computed from it and from already revealed
right values;right is (Y,Y'). Initially they are independent. Every revelation
below is a function of one side and the history,or a kernel of such values,so
the conditional product persists [1.1]. Write Y=(Y0,Y1,Y2),Y'=(Y0',Y1',Y2').
1. Reveal T0=(Y0,Y0'),right-only. Given T0,V and V'=Ext0(X',Y0') are left
functions. As X is independent of T0,D(V|T0)=E_(Y0) TV(Ext0(X,Y0),U_v)<=2^-Ee
(strong affine primitive;X is uniform affine and K00<=rank(X) by(I8)). Couple V to Vbar,uniform given T0,
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
4. Affine Ext2 and the full seed-tuple repair. Put Tpast=(T0,Y1,Y1'),Rfuture=(Y2,Y2'). Let Aux contain all generated left matrix,label and coupling variables through step3,excluding raw X,X'. Its generation depends on X,X' only through V,V' and fixed Tpast. Given Tpast the full left record and Rfuture are independent. In H2,W'_i is fixed;reveal U'_i=Ext2(Y2',W'_i),a right affine message. The original seed pair is common-affine,and Tpast is linear. Thus on every nonempty H3=(H2,U'_i) fiber the right pair is uniform affine. Honest Y starts with rank d_C;fixing the2(d0+d1)+m1 revealed bits gives rank(Y2|H3)>=d_C-2(d0+d1)-m1=K2+Ee by(I4). Extra left messages do not reduce that rank.
Wbar_i is uniform and independent of this right law given H3. Strong affine Ext2 therefore gives D(U_i|Wbar_i,H3)<=2^-Ee,where U_i=Ext2(Y2,Wbar_i). The right conditional law given the complete left record depends on it only through H3 and the fixed coefficient Wbar_i. Consequently the same averaged error remains after retaining the complete left record.
5. Now fix Ghat=(Tpast,V,V',Aux,U'_i). Aux is a kernel from(V,V',Tpast);restricted8.5 therefore leaves the raw latent left variable uniform affine with honest rank loss at most2v. The right future conditional law is uniform affine and independent of that residual left variable. All W_j,W'_j,Wbar_i coefficients are fixed. Hence the whole tuple consisting of U_i,all other U_j,U'_j,and original Y,Y' is uniform affine on this fiber. Use restricted8.2.4's affine section/kernel extension on its U_i coordinate,keeping every Ghat-fixed coordinate fixed. The honest coordinate becomes full uniform Ubar_i;the conditional kernel of all other coordinates is exactly the original one on U_i's support. The repaired WHOLE tuple is uniform affine and its joint TV cost,averaged over Ghat,is delta_U<=2^-Ee (conservatively<=2*2^-Ee). Its construction depends on raw X,X' only through fixed V,V',Aux,Tpast,so the repaired entire seed/right tuple remains independent of the residual affine left fiber. Raw Y,Y' and their equations may change OFF support in this auxiliary law;no claim that their marginal is untouched is made. A maximal coupling matches the complete tuple whenever the honest marginal matches,so this entire repair costs only delta_U in the final joint TV transfer. This repair is proof-only;the evaluator still computes the actual U_j.
6. Flatten the same-label tampered call with U'_i fixed by Ghat. In each round reveal T=Prefix_s(W) (Prefix_sz last). For a B choice additionally reveal Z=J(W,A) (Jz last);for an A choice reveal no Z. Then reveal the final z-bit output Z'_i. Every code has exactly two B choices. If its last bit is A,the record has5s_F+sz_F bits;otherwise it has4s_F+2sz_F<=5s_F+sz_F since sz_F<=s_F. At fixed U'_i,initial I_l is affine in raw latent Z. At each round T is affine;after fixing T,the deterministic G value A is constant. For B,J is affine in W at fixed A;fixing its Z makes the selected G seed C constant. For A,the selected C=A is already constant,regardless of the unused computed J value. Thus the next I_l/F_l/Fz state is affine without revealing that unused value. Induction makes every disclosed coordinate affine on the raw left latent. Each nonempty fiber loses at most5s_F+sz_F+z rank. This left-only record preserves product with the right future at fixed U'_i. By(I5),honest X retains rank k_F;all other X/X' copies are affine images of its same latent. Call this resulting history H5.

7. F step. Given the enlarged history the source tuple (X for the honest call and the other 4
honest rows;X' for the 4 tampered rows j!=i) is independent of the seed
tuple (Ubar_i;U_j,U'_j for j!=i). The honest seed Ubar_i is uniform on
d_F=m1 bits (I7);encode(i) differs from every other label;the 8 tamperings
are exactly t0=8,with exactly two copies of each nonhonest label. This verifies the two-sided-affine selective-transcript contract on every conditioned branch;the tampered row i is fixed in H5 and is
not a tampering. F's two-sided-affine multicopy contract (Section4.2 and restricted8.6,zero shifts)
makes Zbar_i=F(X,Ubar_i,encode(i)) 2^-Ee-close to uniform jointly with every
other row output and every seed. With Qbar=Zbar_i+sum_(j!=i)Z_j and
Q'=sum_j Z'_j (Z'_i fixed),D(Qbar|Q',seeds,H5)<=2^-Ee.
8. Retain Y,Y'. Given (Ghat,flattening record,seed tuple),the outputs depend only on (X,X') and the
seeds,and (X,X') is independent of (Y,Y',coins),since the seeds are
right-sided given H5. So Law(Y,Y'|Ghat,flattening record,seeds,outputs)=Law(Y,Y'|Ghat,flattening record,seeds) is a
common kernel of the real and ideal laws,and appending (Y,Y'),including
unused seed suffixes,costs nothing [1.1]. Averaging the enlarged histories gives synthetic error<=2^-Ee;the conservative bound3*2^-Ee is retained. The affine rank guard has no source-tail event.
9. Restore. Under the maximal whole-tuple repair coupling,real and synthetic (Q,Q',Y,Y') agree unless V!=Vbar,the matrix pair
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
Acyclic chain: AK7f plus conditionaltools8.2.2 and two-sided multicopy8.6 -> F#
(Section4.2);CZ/E6 -> advice condenser;that condenser,canonical LinExt/W2,
Ext0,Ext2,F# -> NB -> outerA# (proof.md,Lemma B9). Global guards and the
family/circuit assembly are in proof.md.
Credits: dbc05fc4-aed7-487a-b708-d1e5b938fbf4 (U,W2,TR/TL);
b15 (103eta precursor);c78568c1-1aa4-4da9-8556-877fb4218f6f
(whole-tuple MC-LA). Current E6 and four-round selective transcript,
entropy and seed-retention proofs are supplied above and in restricted8.2--8.6.

# Sampler part
# Finite support graph,linear sampler,and conditional breaking
Fix Q_S=185349590146229/100000000000000. With L=floor(log2(actual input length)),S=floor(L/Q_S)
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
The selected ordinary maps are outer E=PXTM[JE=(2,),R] of7h and LinExt's canonical matrix family of7g. Each has an exact threshold K0(u,z,e),full seed length D,and a deterministic fixed-seed linear map. For every source X with Hinf(X)>=K0 and an independent uniform full seed Y,
 TV((Ext(X,Y),Y),(U_z,Y))<=2^-e. (1)
E uses7h's literal one-level ceilings;require K0<=u before setup. Its unit w>=max(2^20,e+7+2lg(4u)) gives K0<=Lam*z+g*w. Its seed envelope D<=z/2+(R+5)w requires7g's terminal rate/cap guards;closed slots and B9 pay them. LinExt has K0=z+2e,D=z+4e+lg(4(u+1)),and mean squared seed-TV<=2^(-2e)/2,by7g. Ignored suffixes are retained and shorter output is truncated. No mixed-arity ordinary call remains selected.
Outer E and LinExt retain their supplied ordinary proofs. Ext0/Ext2/I/F/Fz/G/Gz/J/Jz use7f-AK only on the proved uniform affine fibers. Canonical condenser fields are supplied in7a;leaf matrices in7g.
Enumerated base degrees are O(log actual N) at selected calls,so exhaustive base
field setup costs N^O(1);larger relative extensions use polynomial-cost
towers. Linear5/proof.md prove uniform polynomial time in actual N,
without asserting polynomial exact-degree search for arbitrary u.
The selected BASE instead uses linear7f's affine-only three-condenser/
kernel algorithm with its complete seed and integer threshold. That
affine contract alone is consumed in Sections3--5;it is not substituted
for the general-source contracts needed inside the breakers.
### 1.2 Selected all-right-retained affine-breaker contract
For the selected deterministic F and degree r,let Left=(A,A_1,...,A_r) be invertible affine images of one uniform affine Z,and Right=R=(B,B_1,...,B_r,Y,Y_1,...,Y_r) an affine image of an independent uniform affine U. Require rank(A)>=K,honest Y uniform and honest advice different from every tampered advice. The supplied structured-right theorem2.5 gives,with V=F(A+B,Y,alpha),V_i=F(A_i+B_i,Y_i,alpha_i),
TV((V,V_1,...,V_r,R),(U_m,V_1,...,V_r,R))<=eta. (4)
U_m is independent of all displayed right information. Tampered advice may repeat. A smaller degree is padded by repeated permitted copies. These are the actual source classes used by proof.md D3;no arbitrary-source-copy contract is asserted.
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
## 2. Small-bias Cayley graph
See the sampler graph appendix in linear-and-bilinear-primitives.md.
## 3. The affine-base walk
Fix Q_S as fixed above,beta=48/7-2^-79,q_walk=16,C=4(q_walk+1)=68 and
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
Here s<=(48/7-2^-79)S+2 and s+4<=16S.
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
Gamma(alpha,(i,j))=W(row_i(alpha),j),i=1,...,4,j=0,...,16.
The index order is lexicographic in (i,j),giving degree exactly C=68.
Theorem G_A. For every set Q0 of at least Lambda_S left vertices,
|Gamma(Q0,[C])|>5*epsilon0*A0.
Proof. Since m is even Lambda_S is an integer,and
m>=38*2^160>134217728. The complete direct incidence support lemma S1'
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
### 6.1 Class-preserving conditional averaging
If at every fixed value Z=z outside mass delta the complete conditional tuples satisfy the selected structured-affine contract,its TV guarantee eta averages to at most eta+delta,retaining Z. Every source-class,rank,independence and honest-seed hypothesis must be checked on those conditional fibers. Average entropy alone does not establish affinity after arbitrary conditioning. The actual D3 argument instead fixes linear left leakage,retains the whole original affine right variable,and applies the selected contract before any right-tag event;it needs no unsupported nonlinear-fiber lift.
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
affine transformation,advice,CZ/E6 and NB,with restricted8.6's F# interface.
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
4897/4697+nu,relative energy exponent80649/4697+4nu and hereditary cutoff
p^(2-h_R*eta). The finite guard eta*logp>=64/nu pays2^64.
Every actual NB call has eta*(r+a)>=4*(m0+p0*w_N+(23/4)Ee)>=64*w_N>=64*2^80.
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
At the selected call s=2ceil((48/7-2^-79)*S/2),r_F=S,
m=s+32S>=38S>=38*2^160>134217728.
