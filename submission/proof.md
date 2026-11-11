# Provenance and selected structured-source construction
Ports: DX=2560817118466567590901/590295810358705651712,DY=4820629127359044686663/590295810358705651712,RATE=1336016016424583574738741/590295810358705651712.
The complete terminal half-round construction and supplied proofs are retained from published db9ad14b-9db0-43bd-bc58-a68dd49cd916 and the subsequently developed wider-cap candidate. Historical credits in all modules remain; no verdict is a premise. The present specialization changes auxiliary extractor contracts only where the actual application has a proved common affine source. It does not assert a general-source theorem for the new affine primitives.

# A full-B2 circuit-to-affine-flat bridge
All algebra is over F_2. A circuit is a finite acyclic,single-output DAG;
each gate is any one of the sixteen two-input Boolean tables and costs one.
Wires,reuse,depth,sharing and fan-out are unrestricted. Free constants
are allowed,a stronger model for the lower bound.
Attribution: the upstream source named in the supplied proof is GordoAR,
submission 6622943f-7b2d-4f8f-b73c-979c5232484a,version
226b74d1311ff6711c3c373ab7c8e1fc4d26eb0e,proof.md,Chapters II–III;
its public antecedents are 5e9fcf7f-8338-4729-9faa-a18994661af6 and
bf313f35-d97b-4836-b1e8-30afd1157ac2. The supplied proof attributes the
affine-system framework to Mishra–Santhanam,ECCC TR26-260 (2026-10-19).
This reconstruction reproduces that circuit-reduction argument below.
Fix an integer d>=0 and suppose f is nonconstant on every affine flat
of dimension>d. The following reduction supplies the operations used in Section7.
### Representation and potential
During the proof we permit signed references to signals. Signs entering a
gate are absorbed in its truth table. A constant,unary or repeated-parent
gate is replaced by a constant or signed earlier signal. These proof
representations do not assert free native input negation: we start with the
native count and need no conversion of the terminal representation back
to a circuit. The output can be a constant or signed parameter,auxiliary
or gate output.
The six inessential tables are two constants and four unary tables. The
remaining ten are XOR,XNOR and the eight tables
c + (x+b)(y+e),b,e,c in F_2. (2)
Call the latter product gates. These exhaust the basis.
Use M parameter coordinates u,of which 2q are in disjoint pairs. Each pair
ranges over three of its four corners,and every unpaired (free) coordinate
ranges over F_2. Write D for this product domain. Maintain an affine map
Phi:F_2^M -> F_2^n injective on the entire cube,and a representation of
f(Phi(u)) for every u in D. Equality outside D is not required.
Auxiliaries v satisfy a square system
A v = N u + c,A invertible,(3)
with at most three distinct nonconstant variables in each row,after
cancellation. A row has an auxiliary,since A is invertible. Remaining
binary gates read parameters,auxiliaries or earlier gates. Gate outputs
never enter (3). Let s be rows plus gates. A parameter's children are the
distinct rows and gates containing it,not pins and not the output reference.
Set a=M-q. (4)
Initially M=n,q=0,Phi is the identity,the system is empty and s=|C|.
Thus a=n,including circuits omitting inputs.
For each pair,fixing one endpoint opposite its forbidden coordinate
leaves the other endpoint unrestricted. Such safe fixings give a flat
inside D of dimension a; Phi maps it to an a-flat in the original inputs.
Consequently a constant represented output is impossible if a>d. If
a>=d+2,an affine output is also impossible: on a safe a-flat its constant
fiber has dimension at least a-1. A product gate with affine parents cannot
be the output either: on that flat a nonconstant affine factor has a zero
fiber of dimension a-1,and constant factors give a constant output. This covers
constants,parameter/auxiliary literals and terminal gates with affine parents.
### Algebra and cost-free preparations
All coordinate changes below compose Phi with an affine injection. Replacing
a free x by w=x+ell(other parameters),with ell affine,is a bijection of
the whole coordinate cube and of the domain's free fibers; ell may involve
paired coordinates. Fixing a free coordinate lowers a by one. A safe fixing
of a paired endpoint,or restriction of its pair to any two allowed points
parametrized as a line,lowers M,q by one and preserves a. Restricting two
free coordinates to three corners increases q by one and lowers a by one.
Move an affine gate with no product-gate ancestor into (3),processing its
affine predecessors first. Its output becomes a new auxiliary with its
one defining row. The new auxiliary matrix is block triangular with new
diagonal entry one. The row has width at most three; s and parameter use
do not change. Degenerate gates are deleted instead.
The following simplifications decrease s by one,preserve a:
* Delete a constant or literal gate,substituting its output. This includes
constant pins,equal/complementary pins,unary tables,and a gate with
affine parents whose linear parts are dependent. In that last case a
parent is constant or the parents are equal up to a bit,so its result is
a constant or parent literal.
* Delete a gate with no child unless it supplies the output.
* For a row of width at most two,solve for one auxiliary v as a constant
or signed other variable,substitute,and delete that row and column.
Other widths do not increase. Add this row to every other row containing
v,then expand the determinant along v's column: the remaining auxiliary
matrix is invertible.
We need one further elimination. An auxiliary Q with coefficient one on a
free x has the form Q=x+ell(other parameters). Exchange their roles: Q is
free and x auxiliary,without expanding any equation. If b=A^(-1)N_x,
then b_Q=1. Replacing A's Q-column with N_x gives determinant
det(A)b_Q=1,preserving invertibility and row widths. Update Phi by
x=Q+ell. The exchange preserves s,M,q.
If free x has exactly one child,a row E,N_x is nonzero. Some Q thus
has b_Q=1. Exchange x,Q. Now auxiliary x occurs only in E and,during the
loop a>=d+2,cannot be the output by the output exclusion. Delete
E and x,which E determines but nothing else uses. Expanding along
x's one-entry column proves invertibility of the remaining system. This deletion loses one row charge and preserves a.
### Clean all pairs
Perform available affine transfers and the preceding eliminations. Then
apply the first applicable pair rule below,propagate constant/literal
substitutions and short rows,and restart. A fixed parameter deletes each
of its old gate children,since a binary table with one fixed pin is unary,
and each old row child,since that row shortens to width at most two.
Count each removed object once. All four rules preserve a.
1. If endpoint p is unused,safely fix its partner h. The new free p is
unused and q falls by one.
2. If an endpoint has at least two children,safely fix it. At least two
distinct objects disappear;its partner becomes free and q falls by one.
3. Otherwise each endpoint has exactly one child. If partners share that
child B,its gate output or sole auxiliary is a function of the pair
(a row B has just the endpoints and that auxiliary). Two of the three
allowed points give the same bit. Restrict the pair to their affine line.
Delete B after substituting that constant bit. No other object used the
pair,so the new line coordinate is unused.
4. Otherwise suppose p's sole child is a row E. Its partner h has a distinct
sole child B. Safely fix h and delete B. If propagation deletes E too,
two distinct charges have been consumed. If not,now-free p still occurs
only in E: its only original occurrence cannot spread without E's
deletion. Apply the one-child-row elimination,consuming E's distinct
charge. Thus at least two objects disappear and q falls by one.
When these fail,every paired endpoint has exactly one gate child; partners
have different children,and paired parameters do not occur in (3).
In particular every auxiliary depends only on free parameters.
### Prepare a product source
Continue only while a>=d+2. If all gates transfer away,the output is affine,
a contradiction. Choose a product gate G with no product ancestor; its
predecessors are affine and have been transferred. After the simplifications
its parents have nonzero independent linear parts. Turn any auxiliary
parent Q into a free parameter by a role exchange. Such a free pivot exists
because Q is nonconstant and all paired columns of N vanish. If the other
parent is free,independence supplies a pivot different from it; if paired,
Q is independent of that coordinate. With two auxiliary parents,exchange
the first and then use this same argument for the second.
After each exchange,check for simplifications and pair cleanup. If any
is possible,perform a completed action and restart. Otherwise keep this
same G: at most two exchanges make its parents distinct parameters. G is
not the output,by the product-output exclusion,and has a gate child by
dead-gate cleanup. It cannot have a row child.
### Restrictions and the sharing cases
Each restriction lowers a by one. The new safe dimension
is at least d+1,so the restricted output is nonconstant. Freeze old gates
and incidences immediately before the restriction. Propagate topologically,
recording old outputs as constants or signed aliases. Every old gate hit
by a constant pin is removed. If an old output becomes constant,all its
old children are removed too,even if it first became an alias. Charges
refer to distinct removed old objects,not to propagation events.
**Two free parents x,y.** Restrict to the three corners on which G has its
majority value,creating a pair. G and its k>=1 distinct children disappear,
even with shared descendants or children that also read x or y. Thus
Delta q=1,Delta M=0,Delta s<=-(1+k),Delta a=-1.
**Paired p and free x.** Safely fix p's partner h and fix x to kill one
factor of G. Its constant output deletes G. Newly free p loses its sole
use and is absent from the output reference,so is unused. The net
changes are Delta M=-2,Delta q=-1,Delta a=-1.
**Paired p and paired r.** Their pairs differ by pair rule3. Safely fix
their partners h,k,then fix now-free p to kill G's product. The newly
free r is unused,since its sole child G is constant and paired coordinates
were absent from all affine rows. Here Delta M=-3,Delta q=-2,Delta a=-1.
These recipes exhaust the parent types;all substituted coordinates are
restricted through the affine injections specified above.
### Termination and conclusion
Every completed action preserves representation on the current D and
strictly decreases the nonnegative integer M+s: eliminations decrease s,
pair rules decrease M,and each parent-type action decreases M+s. Neither M nor s increases.
Preparatory transfers decrease the number of gates; there are finitely
many. A fixed-source preparation uses at most two role exchanges before
an action. Thus preparation cannot cycle. The geometric exclusions prevent
termination inside a>=d+2,so the process reaches a_T<=d+1.

# Exact finite affine-source and affine-outer assembly
## Selected component contracts
All logarithms are binary with lg=ceil(log2),flg=floor(log2). The supplied modules are linear-and-bilinear-primitives.md (fields,condensers,kernels,ordinary strong linear primitives and incidence), restricted-nonmalleability-and-breaking.md (ring/advice ingredients,conditional tools and the affine-source antichain breaker), and affine-wrappers-and-sampler.md (structured NB,affine outer wrapper and graph/sampler proofs). Its layout map resolves historical component aliases to these included bodies.

The outer wrapper is used on A_i=T_i Z+c_i for uniform affine Z,invertible T_i,and an independent original right tuple (B_i,Y_i) that is an affine image of one uniform affine variable U. Honest Y_0 is uniform. Both complete source-side right information and all original seeds are retained. Only this structured-right contract is asserted. In D3,all B_i and sampler offsets Y_i are affine in the original X_2,so it applies; right tag profiles are not conditioned before it is invoked.

Its C=NB uses two independent common-affine input pairs:the data pair and the seed pair,with honest seed uniform,distinct advice,complete seeds retained. Ext0/Ext2 and inner I/F/Fz/G/Gz/J/Jz use the kernel-only affine primitive AK. Only outer E remains an ordinary strong extractor;LinExt retains its ordinary matrix proof. Restricted8.6 uses h1/h2/s state widths and compresses real histories in honest-B rounds,with whole-state histories in honest-A rounds;its weighted source costs are the20rows of restricted8.6. I_s has its separate fixed gamma and s0s seed;F_s shares F_h gamma/padding;B states have hc bits at their first B coordinate,h2 at their second. In particular,nonlinear right conditioning after C is not called affine. The exact source-rank invariants and the same-label tampered-call flattening are supplied in the NB and F# proofs. The affine-section law8.2.4 and the real right flattening8.6 justify the selected maps;the independent whole-seed repair in NB4.6 is supplied.

AK receives the public original outer context (n,L,Om) and its actual source ambient N (d_O,hc,h2 or inner d_F). It uses no initial condenser: N<=U1=2^(110+lg Om). The two final condensers/kernel rows of7f use rho=e'+225+2lg Om and the public threshold K0=ceil(Az+B rho+xi/G). Its fixed-seed map is linear. Every call tests public source,ambient,field and seed guards before setup; failed public guards return zeros. Its complete rank/error/rounding proof is in7f-AK. Ordinary maps use the proved7g/7h contracts including the wider PT2 cap2^20.

## Fixed rational witness
Set Q_S=185349590576413/100000000000000,DG=937/20,eps_N=10^-18,eps_alpha=10^-40,e_B=1632/7-6*2^-79 and alpha=e_B/(e_B+Q_S)+eps_alpha. All constants are chosen once before n. They are not advice or per-input optimization. The new affine modes are:
AKX: gamma=1195163119977/500000000000.
AKI: gamma=2024732155211/1000000000000.
AKIs: gamma=200144044789/100000000000.
AKF: gamma=1025619179769/500000000000.
AKFz: gamma=2042677938593/1000000000000.
AKJ: gamma=1025619179769/500000000000.
AKJz: gamma=1025619179769/500000000000.
AKG: gamma=84371277291/31250000000.
AKGz: gamma=3293455168861/1000000000000.
AKExt2: gamma=879762190037/250000000000.

For AK,pi=2gamma+eps_N,Lam=1/G,g=B+eps_N,c=0. These reserve the fixed rho and ceiling overheads,as proved in7f-AK. The selected ordinary ordered rates are:
E: 104504625270713/100000000000000.

With unit w of wrappers4.2,the fixed inner slots are:
a=ceil((1025619179769000000500001/250000000000000000000000)w).
az=ceil((2051238359538000000500001/500000000000000000000000)w).
s=ceil((1646727584430500000500001/250000000000000000000000)w).
s0s=ceil((2001440447890000000500001/500000000000000000000000)w).
s0=ceil((2024732155211000000500001/500000000000000000000000)w).
h2=ceil((4539921934571560106723659/31250000000000000000000)w).
hc=ceil((10726571453573620214197319/62500000000000000000000)w).
d=ceil((132774371422722818411173363/500000000000000000000000)w).
k=ceil((736562066178482901395649371/200000000000000000000000)w).
sz=ceil((3293455168861000000500001/500000000000000000000000)w).


The forty-three monotone slot rows,with AK's c=0 envelopes and8.6's right costs40a/42a+11az,are displayed in wrappers4.2. Exact substitution satisfies every row. Their eps_N allowances pay the actual slot ceilings; equality in those nonstrict rows is permitted. This is an explicit feasible witness,not a convergence or optimality claim.

## B9: actual structured NB call and eventual public guards
NB's actual source ambient is n'=d_O=floor(Dslot*Om),target output z=ceil(DX*Om),advice a=s_tag+2+lg C,and requested exponent e. It sets Ee=e+4; F# uses its supplied exponent Ee and its primitive accuracy Ee+7. Let w_N=max(2^80,Ee+44+2lg(4n')) and w_F=max(2^80,Ee+2lg(4n')+32,ceil(z/128)). Source and seed lengths are the literal integer algorithms in wrappers4.2--4.4,with new d0=ceil(p0*w_N) and p0=2gamma_X+eps_N. AK ignores only a padded suffix and retains it. In particular d0 is independent of its output v,so the definition of r0 is acyclic.

Write c_P,pi_P,Lam_P,g_P for the four displayed envelopes. Let d_inner,k_inner denote the selected inner seed/source coefficients. Put D_F=alpha*d_inner+eps_N,K_F=alpha*k_inner+eps_N,
M0=c_2*D_F+pi_2*alpha+eps_N,K2bar=Lam_2*D_F+g_2*alpha+eps_N,
h_R=3+2/DG+2nu,D_E=17+8/DG+5nu,eta=(2-nu)/(14D_E+h_R+5/2),theta=eta/4,nu=2^-80,
Amin=(816/7-3*2^-79)/(1632/7-6*2^-79+Q_S),braw=6-eta/2.
Ext0 has c0=0,pi0=p0,bc0=0 and no output-dependent seed charge. Define
Rnew=(M0+p0*alpha+(23/4)*alpha+56eps_N)/theta-Amin+eps_N,
d0bar=p0*alpha+2eps_N,d1bar=M0+2d0bar+6alpha+eps_N,
DNB=2d0bar+2d1bar+D_F+K2bar+alpha,
K00bar=Lam_X*braw*Rnew+g_X*alpha+eps_N,
Kotherbar=K_F+2braw*Rnew+(5s_inner+sz_inner)alpha+DX+2alpha+100eps_N.
The new source payment is for V,V',the5s+sz selected-branch flattening summaries,the same-label z-bit output,and conservative finite allowances. The structural proof,not this arithmetic formula,justifies that payment.

SC=Slot(DNB)=318651372719292578353/500000000000000000,KC=Slot(max(K00bar,Kotherbar))=21300700210604855727513/200000000000000000,where Slot(x)=ceil(10^18*x)/10^18+2eps_N. Let f=1-2^-40 and SD=f*pi_Ext2. The exact strict density certificate is
4(p0+SD+(23/4)f)-3DG/2-10^-12=211237398838087255551/274877906944000000000000000000>0.

Here are complete eventual,all-length justifications for using the envelopes. As S tends to infinity through actual integer lengths,L/S->Q_S,Om/S->e_B+Q_S,Ee/Om->alpha_0=e_B/(e_B+Q_S),a/Om->Amin and both w_N/Om,w_F/Om->alpha_0. All logarithmic terms are O(log S)=o(S),all fixed integer ceilings differ by less than one,and d_O/Om->Dslot>0. The fixed positive reserve alpha-alpha_0=eps_alpha therefore makes w_N,w_F and every new AK rho at most alpha*Om at every sufficiently large length. Since Dslot>1,n'>=Om eventually. The definitions of w_N,w_F and accuracies Ee,Ee+7 then give rho<=w+200 for every AK role,and eps_N*w>=2^20 eventually. The explicit eps_N reserve proves every displayed padded AK seed and threshold bound. Only public dimension/field tests are evaluated;source rank is an analytic promise,not a data-dependent test. No finite grid is substituted for these limits.

Compute AK_Ext2's literal threshold before using its guarantee. Its output is m1=ceil(d_inner*w_F),and its threshold is at most K2bar*Om eventually by7f-AK's unconditional ceiling calculation and the eps_N reserve. K2bar<Dslot therefore pays K2<=n'. Its field/ambient guards are paid below. Its padded seed is exactly m0=ceil(pi_Ext2*w_N),so m0/Om<=M0 and m0/w_N>=pi_Ext2>SD. There is no mixed-arity or output-dependent seed term in this selected Ext2 call.

The actual r0 ceiling and r=r0-a then give r/Om<=Rnew eventually,with the eps_N margin paying all bounded ceilings. The d0 padding gives d0/Om<=d0bar. From e_L=ceil((d0+Ee-1)/2)+Ee and d1=m0+4e_L+lg(4(Rrow+1)),where Rrow=2(r+a+1),one obtains d1/Om<=d1bar eventually;the remaining logarithm is o(Om),paid by its eps_N margin. Thus d_C<=floor(SC*Om). Since v=6r-floor(eta*r/2)<=braw*r+1,the AK source threshold is at most K00bar*Om eventually. The new k_F+2v+5s+sz+z+2Ee+100 is at most Kotherbar*Om eventually: each of its finitely many ceilings is paid by100eps_N*Om. Therefore k_C<=floor(KC*Om). The outer source rows put KC and SC below Dslot,so all actual source and seed dimensions are legal before their guarantees are invoked.

Every AK ambient N is below2^24*Om<U1;its literal threshold is at most its charged source slot,also below2^24*Om. With gamma in(2,16],xi<2 and rho<=alpha*Om,the two condenser outputs and their positive codimension are below2^26*Om<U1 eventually. All actual fields,kernels and padding are therefore defined. The unchanged kernel-duality proof supplies their strong affine guarantee. All inner AK source thresholds are paid by the displayed rows before extraction. In particular8.6's rank ledger includes the added B-choice A values. Source ambients are n',h or d_F;each is<=n',and all field/seed guards use the same original context.

For advice/ring guards,the r0 definition still gives theta*r0>=m0+p0*w_N+(23/4)Ee+56. With d0<p0*w_N+1,the same elementary floor calculation proves s>=m0+2e_L and the full Lemma U capacity inequality. Moreover p0+SD+(23/4)f>16 and Ee/w_N->1,so the new public guard m0+p0*w_N+(23/4)Ee>=16w_N holds eventually and eta*r0 exceeds64*2^80 and r+a exceeds2^90. The above strict density certificate and m0/w_N>SD give eta*(r+a)>=DG*(a+sig+2) eventually,since(a+sig+2)/w_N->3/2. Also s/w_N tends to p0+pi_2+3>14,whereas a/w_N->1/2;hence s>=a+40 eventually. These prove AG1,AG2,positive r,raw-block capacities and the independent prime/energy guards without using a small p0 alone. DG=937/20 is used in every supplied ring/energy coefficient;all other generic hypotheses are retained.

Each statement holds eventually at every integer S,not only a subsequence: if |g(S)|/S->0,then for any fixed positive reserve its defining inequality holds for all S beyond one fixed threshold. Taking the maximum of the finitely many thresholds gives one common threshold. At smaller lengths the total evaluator tests every literal public guard and returns zero on failure before affected setup. No claimed lower bound at those exceptional lengths is needed.

# Tagged affine-permutation family and outer calibration
The supplied BASE extractor has output d_O,seed s and source threshold K_B. Its exact canonical algorithm and AF1--AF4 proof are7f. The supplied graph G_A has degree C=68,rows A=2^s_tag and the expansion/rank properties proved in the sampler part. The selected affine outer breaker has the structured-right contract stated above,error and rank guards of wrappers2.5,and the following exact slots. All conditions will hold eventually at every integer actual length.

## D1. Algorithm and finite hypotheses
Write L=flg(n),S=floor(L/Q_S),s=2ceil(beta*S/2),ellc=lg(s+4),r_F=S,m_walk=s+32r_F,s_tag=3m_walk=log2 A,w=2. Require S>=2^160,L>=2048 and s_tag>=2L; the latter follows from the fixed constants below.
Build the canonical field F_(2^L),and split x into L-bit coefficients of
P_x,padding its final block with zeros. For a full s_tag-bit row alpha,
let beta,gamma be its first and second L-bit blocks. Define
H_alpha(x)=Prefix_2(gamma*P_x(beta)),a two-bit linear map in x.
Each (beta,gamma) has exactly A/2^(2L) row preimages. For Delta!=0,
P_Delta is nonzero of degree ceil(n/L)-1<n/L,so its root fraction is
at most n/(L*2^L)<2/L<=1/1024. Heren<2^(L+1),L>=2048;
thedegreeisalso<2^L,sotherootboundapplies. Off those roots,multiplication by P_Delta(beta)
permutes gamma and exactly 1/4 of gamma values have zero two-bit
prefix. Thus at most A*(2/L+1/4)<=257A/1024 rows have H_alpha(Delta)=0.
The root bound follows by factoring X-beta at distinct roots. The supplied
finite-field proof makes the canonical search halt in 2^O(L) time.
Full alpha remains in the advice label despite repeated (beta,gamma).
No independence of hash values across different rows is asserted.
Use the source-linear sampler Samp(x,alpha,j),j in [C],of sampler Section5,
seed length d,affine base entropy requirement k_S=K_B and error2^-r_F.
Only uniform affine sources are used. The field index alpha replaces[A].
For t=4C-1 and a=s_tag+w+lg(C),request the padded all-seeds-retained affine
breaker F of output length one,affine-left threshold K,and error at most
epsilon = 1/(2^12 A^2).
The advice for (x,alpha,j) is the concatenation of the s_tag-bit alpha,
lg(C)-bit j-1,and w-bit H_alpha(x). Define
R_alpha(x)=XOR_(j=1)^C F(x,Samp(x,alpha,j),advice(x,alpha,j)),
D(x)=1 iff 2 sum_alpha R_alpha(x)>A.
For the affine application below,assume the source entropy satisfies
k >= max(k_S,K+4w+4Cd). (D1)
All concatenations,zero padding and majority ties are specified. There
are no seed or advice inputs to D. The tags are internal computations on
the actual n-bit x.
## D2. Interface
Linear left leakage and tags leave affine fibers. D3 uses(F_aff') on
each fiber and enumerates right profiles without conditioning the right law.
## D3. Fixed rank-good neighbors and actual parities
Fix independent uniform affine sources X_1,X_2 of dimension at least k.
Put A_0=X_1,A_1=Lambda X_1,B_0=X_2,B_1=Lambda X_2+b,with Lambda
invertible and tau(x)=Lambda x+b fixed-point-free on the sum support.
For each sigma define Bad_sigma as base seeds gamma for which the linear
map E_n(.,gamma) has rank<d on the direction of B_sigma. Its output on
that affine source is uniform on an affine subspace. If rank deficient,
its distance from U_d is at least1/2. The strong affine base error epsilon0
therefore gives |Bad_sigma|<=2epsilon0*2^s.
The graph G_A expands every Lambda_S-row set to more than5epsilon0*2^s
neighbors. Consequently fewer than Lambda_S rows have every neighbor in
Bad_sigma. Intersect the two good row sets to get a FIXED G,omitting
q<=2(Lambda_S-1)<2^-16*A^(1/2) rows. Thus g=|G|>=255A/256.
For every row-copy sigma,alpha in {0,1}xG choose a fixed neighbor j whose
base map has full rank on B_sigma. This choice depends on the source law,
not its sampled point. Its offset Samp(B_sigma,alpha,j) is exactly uniform
on d bits. All expectations below are over the joint X_1,X_2 law.
For a two- or four-set Q of distinct row-copy indices,expand its at most
4C calls. Choose an honest row-copy as in the singleton/fully-doubled
rule below and use its fixed rank-good neighbor;degree t=4C-1 suffices.
Condition the selected left tags,at most4w linear bits,and then the at
most4Cd complete left seed-leakage bits. Every nonempty fiber is uniform
affine and each honest A_sigma has entropy>=k-4w-4Cd>=K. The right law
X_2 is unconditioned and independent of this left fiber. Source-linearity
splits every seed into its fixed left leakage and right offset;the honest
seed remains uniform by the fixed full-rank choice.
Now enumerate the at most2^(4w) selected RIGHT tag profiles. For each fixed
profile with honest advice distinct from all tampered advice,apply(F_aff')
and its WHOLE-right retention to this conditional left law and R=X_2.
All right offsets and seeds are affine images of the original uniform affine X_2,independent of the fixed left fiber. Thus the structured-right hypothesis of the selected outer theorem holds before any right-profile event is tested. This bounds the absolute parity expectation restricted to the right-profile
event by2epsilon: that indicator is a function of retained R,and the ideal
honest bit is independent uniform. We do NOT condition the right law on
that event or assert uniform seed after such conditioning. Sum these
restricted expectations,then average the left fibers. Thus each legal
parity has bias<=b_Q=2^(4w+1)epsilon. If an alpha occurs only once in Q,
choose such an index as honest;all profiles are legal because alpha,j
identify that call. Otherwise Q is fully doubled;choose its least index.
Only profiles with a selected collision can then be illegal. Charge their
probability separately,giving bias<=b_Q+Pr[at least one selected collision].
The tag root count bounds collisions at EVERY(X_1,X_2) point by
N_col=257A/1024. No independence of collision events is asserted.
Use a_Q=1 for four-sets;for two-sets use A^2 for matched copies of one row
and A otherwise. The same counting gives sum_Q a_Q<=(2/3+3/256)A^4.
Hence Bsum=sum_(Q on G) a_Q*b_Q
<=(2/3+3/256)*2^9*epsilon*A^4
<=521A^2/6144<3gA/35,by epsilon=1/(2^12*A^2),w=2 andg>=255A/256.
The graph,rank argument and error budget require no pointwise bad-seed
set and no exceptional right-source probability. The two actual row-copy
arrays and every collision term are retained for the quartic identities.
## D4. A robust quartic sign certificate
Let X_alpha,Y_alpha be actual signs on G and S_0=sum X_alpha,S_1=sum Y_alpha.
Put Z_alpha=X_alpha Y_alpha,Z=sum Z_alpha,zeta=E Z. Only matched pairs occur in zeta,each with weight A^2. §D3 gives
|zeta|<=N_col+Bsum/A^2<N_col+3g/(35A).
Since g>=255A/256 and A>=2^256,
r=|zeta|/g<=257/1020+3/(35A)<13/50.
For Q_+(x,y)=xy(x-y)^2 and Q_-(x,y)=-xy(x+y)^2,expand the ACTUAL
ordered quartic monomials. All-even terms contribute-2g^2. A matched-pair
parity has coefficient plus or minus2(3g-2),and fully doubled two-row
parities contribute-8 sum_(alpha<beta) E Z_alpha Z_beta.
Every remaining parity has a singleton row and is bounded by b_Q.
A size-four parity has coefficient6,6,or-8 according to its two-copy
counts(3,1),(1,3),(2,2). For size two,a mixed nonmatched pair has absolute
coefficient6g-4; a same-copy pair has coefficient-4g. For a nonmatched pair its coefficient is at most6A and a_Q=A; for a four-set it is at most8 and a_Q=1. Thus every remaining absolute coefficient is at most8a_Q.
Thus the exact identities give
E Q_+(S)=2(3g-2)zeta-2g^2+4g-4E Z^2+u_+,
E Q_-(S)=-2(3g-2)zeta-2g^2+4g-4E Z^2+u_-,|u_+|,|u_-|<=8Bsum.
All doubled collision terms are retained exactly,not charged separately.
The variance identity E Z^2>=zeta^2 now bounds BOTH expectations by
2(3g-2)|zeta|-2g^2+4g-4zeta^2+8Bsum
<=[6r-4r^2-2+4/g+24A/(35g)]g^2
<=[-444/625+1/1000+6144/8925]g^2<-g^2/50. (D4')
Here6r-4r^2 increases up to13/50,A/g<=256/255,and
4/g<=2^(3-2L)<1/1000 forL>=2048.
Marginal fourth moments have all-even part3g^2-2g. Remaining parity
coefficients are24 for size four and12g-16 for size two,each bounded by
24a_Q. Thus
E S_sigma^4<=3g^2+24Bsum
<[3+72*256/(35*255)]g^2<16g^2.
Their L4 norms are consequently at most2sqrt(g). This weaker norm bound
is sufficient for the deterministic removed-row perturbation below.
Let Y_sigma be the ACTUAL sum over ALL rows in that copy. Pointwise
|Y_sigma-S_sigma|<=delta=q<2^(-16)A^(1/2). Since g>=A/2,
delta/sqrt(g)<2^(-16)A^(1/2)/(A/2)^(1/2)=2^(-31/2)<1/40000.
The L4 norms of S_sigma are at most2sqrt(g),and those of Y_sigma at
most3sqrt(g). The norm triangle inequality follows from Holder applied
to (|X|+|Z|)|X+Z|^3 and division by the cube of the L4 norm (zero norm
is trivial). For four factors,Holder itself follows by Cauchy-Schwarz
twice: E|abcd|<=(E a^4 E b^4 E c^4 E d^4)^(1/4).
Telescoping the four factors in any degree-four monomial shows its
expectation changes by at most4delta(3sqrt(g))^3. Each Q has sum of
absolute coefficients four. Thus
|E Q_+(Y)-E Q_+(S)|,|E Q_-(Y)-E Q_-(S)|
<=432delta g^(3/2)<(432/40000)g^2<g^2/80. (D5')
The expectation margins are -g^2/50 and the perturbations are less than g^2/80,so both actual expectations are strictly negative. If the actual majority
XOR were identically zero,Y_0Y_1>=0 pointwise and Q_+(Y)>=0,impossible.
If that XOR were identically one,Y_0Y_1<=0 and Q_-(Y)>=0,also impossible.
The weak product inequalities include all ties: the evaluator outputs one
exactly for a negative sign sum,and zero for a nonnegative sum.
Thus both actual majority-XOR values occur. No sixth or eighth moments,
small-tail inference,or independence of the full sign vector is used.
All sums and collision events above are actual deterministic quantities.
For an affine S=a+V of dimension at least k,take X_1 uniform on V and
X_2 independently uniform on a+V. Their sum is uniform on S,and D(x)+D(tau(x)) depends only on that sum.
The joint-law contradiction therefore proves both values occur on S. Thus the property holds
for every invertible affine tau without fixed points on S,and includes
every nonzero translation.

## 2. Fix every constant before seeing n
Fix M=2^70,beta=48/7-2^-79,q_walk=16,C=4(q_walk+1)=68,p=144,t=4C-1=271 and j=lg(t)=9. For actual n>=1 set L=flg(n),S=floor(L/Q_S),s=2ceil(beta*S/2),ellc=lg(s+4),r_F=S,R(s)=0. Test L>=max(2048,p+128),S>=2^160 before further logarithms. Set
m_walk=s+32S,s_tag=3m_walk,A=2^s_tag,E=2s_tag+12,e=E+18,Om=e+L-43=6s+192S+L-13.
These are literal integer dimensions;every later public guard is tested before affected setup. Error2^-E equals1/(2^12 A^2),and9*3^9-8+1/4<2^18 pays the wrapper's internal exponent. Tags,fields and seeds are internal computations,not inputs or advice.

Because beta*S<=s<beta*S+2,E1=6beta+192=1632/7-6*2^-79 gives Om<E1*S+L-1<=(E1/Q_S+1)L-1. In particular Om<=c_O*L eventually,where c_O=(433/Q_S)+1+10^-15 and433>E1. Also s>=6S,m_walk>=38S,and3(beta+32)<p<2p;hence s_tag<=2pS eventually. With S<=L/Q_S andceil(2p/Q_S)=156,the row count is at most n^156 at every sufficiently large actual length. Since3(beta+32)>2Q_S,the tag length is at least2L eventually. The evenness of s,m_walk makes Lambda_S=2^(3m_walk/2-17) an integer. All fixed graph/field thresholds are eventually met as S grows. The sampler's complete finite walk proof still applies for every S>=2^160.

## 3. The outer affine breaker A# and (F_aff')
Use the structured-right wrapper2.1/2.5 with C=NB,t=271,j=9 and a terminal one-bit output. Let U18(x)=ceil(10^18*x)/10^18 and Slot(x)=U18(x)+2eps_N. The initial shared affine compression remains latent. C/E read the original right input of ambient d_O;Lr reads the same stored left compressed image;Lm reads its actual r-bit left state. No final Qr or r-bit Lr state is computed.

The fixed affine modes are gamma_0=2000386497751571/1000000000000000,gamma_R=2186340882999/1000000000000,gamma_M=2186340882999/1000000000000. For each set xi=gamma/(gamma-1),G=1-(xi-1)^2,A=1/G,B=xi*gamma/G. Their normalized source/seed coefficients are L0:(0,1+2gamma_0*alpha,A_0,B_0*alpha),Lr:(0,2gamma_R*alpha,A_R,B_R*alpha),Lm:(0,2gamma_M*alpha,A_M,B_M*alpha). For E use its ordinary7h coefficients with pi_E,g_E multiplied by alpha; c_E,Lam_E unchanged.

Define S0=Slot(1+2gamma_0*alpha),S_R=Slot(2gamma_R*alpha),S_M=Slot(2gamma_M*alpha),KM=Slot(A_M*(DY+eps_N)+B_M*alpha),
TRbar=272*(S0+SC+17*(DY+DX)),TLbar=272*(SC+16*DY),
Dslot=Slot(max(272*(S0+SC)+KC+alpha,TRbar+Lam_E*DX+g_E+271*DX+alpha,S0,Lam_E*DX+g_E,KC)),
Kslot=Slot(max(TLbar+A_R*RATE+B_R*alpha+271*RATE,A_0*SC+B_0*alpha)).
The exact selected values are Dslot=5624049341930981484811/20000000000000000,Kslot=415036693217829377743421/500000000000000000,S0=4969217542531825033/1000000000000000000,S_M=2169096471234001529/500000000000000000,KM=21010562287279495589/500000000000000000.
Use literal dimensions dy=ceil(DY*Om),dx=ceil(DX*Om),r=ceil(RATE*Om),s_C=floor(SC*Om),s0=floor(S0*Om),k_Cslot=floor(KC*Om),d_O=floor(Dslot*Om),k_O=floor(Kslot*Om). All ports and source slots are independently rounded.

The exact port certificates are DY-c_E*DX-pi_E>16/M,DX-2gamma_R*alpha>16/M,DX-S_M>16/M and RATE-272DY-KM>16/M,with M=2^70. Their positive fixed margins eventually pay every bounded extra ceiling,including272dyceilings in the Lm rank guard. Every selected lambda in E,L0,Lr,Lm is below2^20. The fixed counts t=271,j=9 imply that converting normalized ports to literal ceilings changes each source-ledger right side by less than2^30: the transcript additions number at most4j(t+1)+t,and every source term increases by at most its selected lambda. Each Slot's surplus eps_N*Om before its final floor eventually exceeds2^31,paying these changes,the floor and both extra entropy bits. Thus all literal wrapper source rows hold:
r>=(t+1)dy+floor(KM*Om),
d_O>=(t+1)(s0+s_C)+k_Cslot+e+1,
d_O>=TR+K_E+t*dx+e+1,
k_O>=TL+K_R+t*r,and k_O>=K_0.
The separate d_KE row puts E's threshold below its actual d_O input before invoking its ordinary guarantee. The s0 and KC rows also give s0,k_Cslot<=d_O eventually. NB's seed/source prefix bounds are B9. Seeds padded beyond their literal primitive lengths remain fully retained.

The generic7f-local proof supplies the actual L0/R/M source and seed bounds with public context(n,L,Om). Its pre-field checks use k_O,K0<2^24Om,xi/G<2^12 and 2/N0+eps_N/4<2eps_N,N0=2^60. These exact scalar inequalities hold for the witness. The shared C0 output is below2^86Om<U1 and each later condenser output below2^26Om<U1. As Om=O(L),U1=2^(110+lg Om)<=2^L<=n eventually. Lm needs no initial condenser because RATE+eps_N<2^12 and r<=2^12Om<U1. Its separate floor(KM*Om) source slot and full padded seed are paid before extraction. The public guards are tested before any affected setup;all hold eventually by these fixed strict margins and the supplied kernel-duality argument.

For ordinary E,its actual data ambient is d_O=O(Om) and its unit requirement is below w_N. Its source threshold is paid by d_KE. Its output dx obeys the ordinary cap DX*440<256*2*233 and actual caller unit w_N>=233S,while Om<=440S eventually. Its rate,height,arity and kap<=2^20 guards are explicit. Thus its strong general-source guarantee is invoked legally;no affine claim about the post-C right law is needed.

Before C is called,H0 fixes the original right prefixes and their right linear C0 contributions. The original tuple(B,Y) is affine in an independent uniform affine U,so U|H0 is uniform affine on every nonempty fiber: E0 is determined by the fixed honest prefix. Honest Y0 loses at most(t+1)(s0+s_C) rank and hence meets structured NB's threshold by the d3c row. The R1 pair is also common-affine in the independent original left Z. Its honest marginal has the supplied averaged near-uniform guarantee. Wrapper2.2 supplies the affine-section seed-pair repair and3eps TV lift,so the two-sided NB theorem applies here. Later E calls are ordinary,so subsequent nonlinear right conditioning causes no unsupported affine invocation. The left rank proof and nine-level statistical recurrence of wrappers2.5 are unchanged. This establishes(F_aff') with error at most2^-E,retaining all original right information and every seed.

The combined raw source budget is still4C*d_O+k_O+8<=(272Dslot+Kslot)Om+8. Put c_O=(3p+1)/Q_S+1+10^-15=43485349590576413185349590576413/185349590576413000000000000000,delta_O=ceil(Dslot*c_O)=65973643,kappa_O=ceil(Kslot*c_O)=194745677. The strict smaller logarithmic coefficient is established in Section8 below,independently of these conservative coarse witnesses.

## 4. One fixed all-length guard and entropy coefficient
Put eps=2^-18,eps_B=10^-18,eps_R=2^-80 and fix
T=1,
H_raw=4C*delta_O+kappa_O+T=18139576574,
H_base=ceil(A_B*delta_O+B_B*(1+eps_B)/Q_S+eps_R)=18139576070,
H=max(H_raw,H_base)=18139576574. (O3)
Here H>=delta_O,which supplies the separate d_O<=n guard.
Here4w=8<=L=T*L. Thus k_F+4w+4Cd
<=[kappa_O+4C*delta_O+T]L<=H*L.
AF4 separately proves k_S=K_B<=H_base*L<=H*L. The final eps_R in
H_base pays the K_B ceiling. This proves(D1) eventually at k=HL,including all4C
linear seed leakages;the affine-fiber argument pays no entropy tails.
For integer floors S>=L/Q_S-1 and 123>3Q_S,L>=p+128 imply
123S>=123L/Q_S-123>3L-123>=2L.
Since linear7f proves s>=6S,s_tag=3(s+32S)>=114S>2L,and every ordered pair of L-bit tag blocks has
exactly A/2^(2L) preimages among the COMPLETE row indices. The tag
collision proof uses A itself,with no stronger lower power of n.
Even s,m_walk make Lambda_S=2^(3m_walk/2-17) an integer. Also
m_walk>=38*2^160>134217728,paying the incidence characteristic guard.
The affine-only sampler is used only on uniform affine B_0 and B_1.
The total evaluator tests these literal finite guards:
n>=1; L=floor(log2 n)>=max(2048,p+128);
S=floor(L/Q_S)>=2^160; s_tag>=2L; HL<=n. (5'')
It also computes/tests every actual AF1 and primitive length/seed/entropy
parameter guard,returning zero on any failed public check before setup.
AF1--AF4,the eventual B9 bounds and the outer requirements prove that every public guard passes at every sufficiently large integer length. Any additional finite failure returns zero before the affected setup.
Since S tends to infinity at EVERY integer actual length,all fixed
lower thresholds eventually hold. Finally n>=2^L and2^L/L tends to
infinity,its consecutive ratio being2L/(L+1)>=4/3 for L>=2.
At L=128,H*128<2^128 by exact arithmetic. The same increasing ratio
proves HL<=2^L<=n at EVERY L>=128. Thus all sufficiently large S supply
this native guard,and gives ONE common fixed threshold N for ALL
integer n>=N. All exceptional lengths,including n=0,use the zero branch.
## 5. The fixed evaluator
Precompute the constants in §2 by their finite algorithms. On an
input x of actual length n:
1. If n=0,output0. Compute L=floor(log2 n),S=floor(L/Q_S).
If L<max(2048,p+128) or S<2^160,output0 before any logarithm of S.
Otherwise compute the remaining public parameters. If s_tag<2L,
HL>n or any exact supplied parameter guard fails,output0.
2. Build the exact base extractor E_n from (L_A),seed-padded to s bits.
Its data argument is x and its output length is d.
3. Build the (G_A) graph on left s_tag-bit strings: the four rows of the AB+C
step,the short-start Cayley walk and the s-bit prefix (sampler-and-
disperser.md §4). Its degree is exactly C.
4. Build the canonical field of degree L and the tags in §D1.
For each alpha in lexicographic {0,1}^{s_tag} and j=1,...,C compute
y_(alpha,j)=E_n(x,Gamma(alpha,j)),
b_(alpha,j)=F_n(x,Prefix(y_(alpha,j),d_F),
alpha || encode(j-1) || H_alpha(x)),
using the ONE breaker parameter setting (O2). Let R_alpha be the XOR
of these C bits.
5. Output 1 iff strictly more than half of the A row bits are one;
otherwise output zero. This is exactly §D1's evaluator.
F_n is the raw d_F-seed breaker,identical to its padded form.
Source-linearity survives graph index selection. The degree and advice lengths match §D1,and (D1) was checked
above. §§D2--D4 prove the affine-permutation property at k=HL,
and in particular the directional property.
## 6. Uniform polynomial time in actual input length
All fixed rational/constant searches terminate before n;parameter arithmetic
has O(log n) bits. There are A*C<=C*n^ceil(2p/Q_S)=C*n^156 calls.
BASE and AInit binary field degrees are O(log n),so their canonical
irreducible searches,condenser evaluation and kernel row reduction cost
polynomial time in actual n. Both outer left maps use the explicitly analyzed AInit.
For every used PXT/PXG/PXTM leaf,2^29*w bounds its terminal ambient N,
output ell and accuracy a by O(log n),Remaining matrix calls have the stated logarithmic terminal bounds.
The canonical matrix algorithm is polynomial in2^N,
2^(ell+4a+lg(4(N+1))),N,ell,a,hence polynomial in n;weights have
polynomial bitlength. NB's LinExt also has Rrow,m0,e_L=O(log n).
Graph primes have O(s_tag)=O(log n) bits and searches take2^O(s_tag);
all row/walk evaluations,tags and O(log A)-bit majority counting are
polynomial. Polynomially many fixed-depth calls preserve polynomial time.
The n^156 bound counts rows,not the entire evaluator time exponent.
## 7. Native circuit bound,NP eligibility and expression
The assembly is nonconstant under every fixed-point-free affine permutation
on a source flat of dimension at least k=HL. In particular,for every
nonzero v and every affine k-flat S its translation derivative is
nonconstant. We prove a finite bridge under the WEAKER hypothesis
for every nonzero v and affine k-flat S,some x in S has f(x)!=f(x+v).
For integers k>=1 this implies C_fullB2(f)>=4n-4k. For n<=k this is trivial.
Otherwise use the retained affine-map,three-corner-pair bridge with d=k-1,
stopping at a=M-q=k. Separation implies nonconstancy on every flat of
dimension at least k,by choosing a nonzero tangent translation of any
hypothetically constant flat.
At every state with a>=k an unused parameter is impossible: for a free
parameter keep its coordinate direction in a safe a-flat; for a paired
parameter fix its partner safely. Independence gives zero derivative
there,and injective Phi preserves its nonzero direction. Absence includes
the output reference. Pair rule1 is thus impossible; rule3 creates an
unused free line coordinate at the same a and is impossible as well.
Mixed and paired/paired productive actions create the unused coordinate
identified in the retained ledgers at new a>=k,so only FF actions occur.
Put W=s-2q. Each FF action removes at least two objects and adds one pair,
dropping W by at least 4 as a drops by one. Basic deletions,transfers,
role exchanges and one-child-row elimination do not increase W. Each
remaining pair rule2 or4 removes at least two objects and one pair,also
not increasing W. Termination thus gives exactly n-k FF actions and
W_T<=W_0-4(n-k)<=|C|-4(n-k).
All parameters lie in the output's undirected incidence component;
otherwise the invertible auxiliary equations split by components and
the output is independent of a parameter,already excluded. In each
component the row/auxiliary block is square: full row and column rank
of invertible A force both inequalities for that block. Each gate-output
variable belongs with its defining object. Thus the output component has
M+s_0 variable vertices and s_0 hyperedges of size at most three. A connected
hypergraph of s_0 such edges has at most 2s_0+1 vertices,since adjoining
an intersecting edge adds at most two. Hence s>=s_0>=M-1. At a=k,
W_T=s-2q>=k-q-1>=-1.
The only integer case with W_T<0 is q=k,M=2k,s=M-1. The tight incidence
tree proof and paired read-once zero-period lemma,proved in
sampler-and-disperser.md,exclude exactly this case under our weaker
separation hypothesis. Thus W_T>=0,proving the finite bound. Its uniform
additive constant zero is optimal for this weaker hypothesis: n=k=1,
f(x)=x satisfies separation and has zero gates.
For all sufficiently large actual integers n,the bridge gives the coarse
bound4n-4H*flg(n); K=4H=72558306296. The refined remainder is proved below.
All constants and exceptional zero branches define one fixed total P
evaluator,with f_0=0. The NP verifier accepts (x,empty) exactly when it
outputs1,and rejects nonempty witnesses (witness bound zero). Tags,seeds
and fields are internal. The coarse remainder is exactlyKlog2n; log2n/n tends to
zero,giving both joint projections(4,1).
The actual stretched row length obeys(O2b):
48q_walk*ellc+6<=pS and s_tag<=2pS<=(2p/Q_S)L. The exact k40
comparison36726<144*2^40 and doubling prove this at every integer S.
Thus every2^O(s_tag) enumeration is polynomial in actual n,including
all row,prime and field work.
## 8. Exact logarithmic loss without a fixed loss
Keep the single fixed evaluator above,with the new unit Om=e+L-43.
Put Q=Q_S,E1=1632/7-6*2^-79,J=272Dslot+Kslot and
A=J*(E1/Q+1). The exact declared logarithmic coefficient is
c1=117307207108553817427868119232279499613765836146688331087664551/2991709403838695980010008563679232000000000000000000=4A. There is no constant monomial in bound.json.
Since s<beta*S+2 and e=6s+192S+30,
Om=6s+192S+L-13<E1*S+L-1<=(E1/Q+1)L-1.
Both entropy floor slots therefore satisfy
4C*d_O+k_O+8<=J*Om+8<A*L-J+8<A*L-1,
as the exact fixed J>9. This pays the raw source requirement,including
all left tags/seed leakage,without a fixed positive loss.
For BASE let a=A_B,b=B_B,u=Dslot. Its entropy requirement is bounded by
A0*L+B0*log2(L)+C0,where
A0=a*u*(E1/Q+1)+b/Q,B0=2b,C0=(M+113)a*u+227b.
This remains an upper bound with the smaller unit. The exact rational
margin A-A0 is positive. Thus there is one integer j0 such that for all
j>=j0,(A-A0)2^j>=B0(j+1)+C0+2:choose a base large enough;the
left doubles and the right has fixed increment B0. For2^j<=L<2^(j+1),
this gives the BASE requirement<=A*L-2. Hence it holds for EVERY
sufficiently large integer L,not only the dyadic sequence.
Choose proof-only k=ceil(A*L)-1. Then k>=A*L-1 pays both source
budgets,and k<=A*L. Also A<H-1,so eventually1<=k<=H*L<=n.
D3 applies whenever those two explicit source budgets hold,so its full
affine-permutation property holds at this k. The evaluator never receives k.
The native bridge gives C_B2(f_n)>=4n-4k>=4n-4A*L
>=4n-c1*log2(n) for every sufficiently large ACTUAL integer n.
All public guards and exceptional zero branches are those proved above;
the small-length branch returns zero before the remaining parameters.
The fixed evaluator uses the stated unit and primitive modes. Its exact
logarithmic coefficient is explicitly declared above;no fixed loss term
occurs. No comparison routine or checker verdict is a
mathematical premise. Since log2(n)=o(n),the joint projection is(4,1).

Material source:771814cc-1e2c-4f72-8e53-b9e1a4125e14;historical credits remain above. Selected parameters and the seed-margin proof are fixed above. Verdicts are not premises.
