# Selecting a uniform family by its proved source certificate

J_N=family.evaluate uses ALL N bits, the supplied multiplicity-weighted
SU extractor, half-accuracy linear cubic reservoir, plateau/tail budgets,
variable active gaps and minimum weighted source selection. At EVERY
N>=0, t_N<=u_N<=d_plan and C_B2(J_N)>=4N-8t_N+8>=4N-8u_N+8.
Fallback is output0 and t_N=u_N=d_plan=N+1. All definitions and
finite/eligible/all-length proofs below concern this ONE fixed program.

The public ranges are h[.1,.2],a[.4,.7],wm[1.3,2],w0[2.4,3.4],
X[.25,.6],p[.2,.35]; beta remains[.15,.53]. X is root spacing;
D=w0+X and later active gaps grow as the entropy schedule declines.
All maps/finite guards use the NEW rho^-1 SU field bound. Its complete
multiplicity/soft-weight proof is in extractors.md; output maps use
their actual smaller fields. Its cost-aware interpolation and entropy proof are supplied in full.

# Exact energy and weighted-confidence witness

Senergy is defined by the two finite certificates(PC) below. Use

 a=6249/10000,D=55933/20000,w0=25419/10000,
 wm=D/2+10^-6,X=D-w0=1019/4000,
 p=a/2-(102/10^8)*w0-10^-13,
 h=23868452849/200000000000.                              (FM1)

All range/prefilter inequalities hold exactly, including2wm>D,
a/wm<.53,a/w0>.15 and p<a/2-.00000102w0.
w0+X=55933/20000=2.79665<D_energy=139833/50000=2.79666.
Put beta=a/w,alpha=1/(1-beta)+.00001 and
f(w)=w*(1-beta)/(1+alpha*beta).

Put e=1/100000,x=w/a>1,Qpoly=x^2+e(x-1),
g(x)=x(x-1)^2/Qpoly, so f(w)=a*g(w/a). Direct differentiation gives

 g'(x)=(x-1)*[e(2x-1)(x-1)+x^3+x^2]/Qpoly^2,
 g''(x)=2*[e^2(x-1)^3+3e*x*(x-1)^2+x^3]/Qpoly^3.

All terms are positive. Since e<1,Qpoly>=x^2, dropping negative terms
in g' bounds it by1+3/x+1/x^3<=5<=11. Hence f is increasing, convex
and11-Lipschitz. Convexity gives f(mid-t)+f(mid+t)>=2f(mid);
integrating each symmetric pair makes floored midpoints a LOWER
integral bound.

The audit executes M=65536,Delta=(w0-wm)/M,scale10^15. Floored
f(midpoints) give a LOWER integral; floor(scale*f(right))+1 gives an
UPPER integral. All remainders and four original Fraction terms for
EACH sum are checked. Exact totals:

 MID=60580679463606431770,RIGHT=60581188431472866242.

The energy ratio R=(864h/109)*Delta*MID/scale is
2232320000187364700858960654909830977
/2232320000000000000000000000000000000>1.
With W=w0-wm,mu=R-1-8hW/10^12>0,
J=1+ceil(176hW^2/mu)=331581749348. Since864/109<8, the SAFE larger
11-Lipschitz conversion loss yields
R-8hW/10^12-88hW^2/J>1+mu/2>1.
This certifies original-scale left boxes without executing J terms.

For the second certificate use the same UPPER right mesh. At left w_i,
Iupper=Delta*sum_(j>=i)RIGHT_j/scale>=integral_w_i^w0 f.
Set r_lower=max(0,109/216-4h*Iupper). With Phi defined below,
Phi(r_lower)<=Phi(r(w_i)). The scale for Phi is s=18*10^12.
Write r_lower=n/A. If n=0 use q_i=0. If54n<=A use
q_i=floor(6n*s/A). Otherwise use
q_i=floor(s*sqrt((216n-A)/(108A)))-s/18.
Thus q_i/s<=Phi(r_lower) in EVERY branch; the square-root branch is
computed by integer square root. All division remainders and both integer
square bounds are checked exactly. Because Phi(r(w))/G(w)^2 increases,

 Jlow=Delta*sum_i floor(scale*q_i/[s*(D-w_i)^2])/scale
 <=integral_wm^w0 Phi(r(w))/(D-w)^2 dw.

The exact integer sum is137445282220259173291; the three branches have
2,4814,60720 cases. The audit also checks the running suffix ends at zero.
It verifies the final confidence upper bound

 2h*[17/(18X)-Jlow]
 =9389094581652324953659313755360504831957
 /30051532800000000000000000000000000000000<p.

Both certificates use actual exact finite arithmetic.

# Finite certificate theorem for the current public-grid family

cubic.md/extractors.md prove the finite primitives; elimination.md proves
the full-B2 transfer. These conditions are implemented in family.py.

Fix a parameter tuple accepted by all original numerical component guards, the tuple-row guard, and the additional guards ell*P<=N, d<=N and Amax<ceil((w0+X)Z)-ceil(w0Z). Let the actual interval tree have depth H, products V_j, arities t_j=2^A_j, block counts b_j, component thresholds K_j, active depths a_0=0<...<a_hp=H, fixed-node response maps E_j and confidence pairs(p1_u,p2_u). Define

 rho_H=0,rho_j=rho_(j+1)+binom(t_j,b_j),w_j=max(1,rho_(j+1)),m_j=w_j*P,

 T=2^max_u(V_(a_(u+1))-V_(a_u)),
 C=1+sum_(u<hp)(T_u-1), T_u=2^(V_(a_(u+1))-V_(a_u)), S=ceil(d/C), Dtot=P*sum_(j<H-1)(t_j-1)*w_j,
 D_res=P*(sum_(j<H)w_j+max_(j<H)w_j), s0=a_1, F_j=1 for0<=j<=s0,and F_(j+1)=2b_jF_j for j>=s0.

The exact geometry conditions are S-Dtot>=max_(j<H)K_j+1,
 S//F_(j+1)-Dtot>=K_j+1 for every s0<=j<H,
 (S//F_j-Dtot)//(2t_j)>=K_j+D_res+1 for every s0<=j<H,
 S//F_H-Dtot>ceil(N/2^V_H).

All actual blocks use injective prefix embeddings into their local ambient n_j. These preserve conditional affine entropy; target input length stays N. For every accepted tuple the factory clipping guards give Wmin^3>=(globalL+1)^2 and Wmax^20<=Lpad^19,while Lpad<=L_j<=globalL and Lc_j<=globalL+1. Thus (Lc_j)^(2/3)<=W_j<=(Lc_j)^(19/20) pointwise. Each actual context satisfies the checked size and rate interval conditions. Every finite primitive use below therefore lies in its stated domain; eventual estimates are unnecessary for this domain check.

For an arbitrary affine source of dimension d, the active-pair lemma
below with targets S=ceil(d/C)>c_hp has threshold at most d. It gives
consecutive active g,g' with both useful entropies>=S after fixing all
preceding coordinates at their conditional ranks.

Every subsequent prefix fixing inside g' preserves H(X_g|X_g') exactly by the affine fiber dimension identity; every rank-r response fixing decreases it by at most r. Force earlier siblings on ancestors of g' to fail at width P, using the exact challenge-response deficiency lemma while the parent's entropy is at least S-Dtot>maxK. The sibling's coordinates are fixed, so its challenge is constant. A full uniform response exists because a strong linear extractor of error1/4 has some seed with uniform affine output: a nonuniform affine output has distance at least1/2. Fix m_j=w_j*P response coordinates at that parent. This costs exactly that rank and leaves all previously forced events true pointwise.

Descend ordinarily from g'. Every possible g' has depth at least s0=a_1. At ordinary depth j>=s0,use the product F_j that starts at this first possible descent level,even when the actual descent starts later. Retain a nice source with entropy at least S/F_j-Dspent. The binary ordinary-descent lemma below uses high-child threshold parent_entropy/(2t_j). Either it finds b_j hidden blocks or it keeps a child at entropy greater than parent_entropy/(2b_j). Fix its predecessors and force them to FAIL BEFORE the next move,while the parent projection is at least the retained child and hence at least S//F_(j+1)-Dspent>K_j. Here F_(j+1)=2b_jF_j. Charge all response constraints to Dtot;division by2b_j scales old rank losses down. Every previous deterministic FAIL persists pointwise,and later prefixes inside g' preserve the outer conditional entropy. The complete winning/losing entropy calculation is supplied below.

The procedure cannot reach a leaf, by the final geometry condition. At its stopping internal node e of depth j, every hidden block has entropy at least (S//F_j-Dtot)//(2t_j)>=K_j+D_res+1. It therefore remains a block source above K_j after every possible extra rank-D_res conditioning of the full source. The conditional entropy at g given g' remains at least S-Dtot>K_g. All earlier-sibling tests on the path to e fail identically. No property is assumed about retaining a large earlier threshold after all later restrictions.

The pointwise primitive error used by the controller is obtained from actual constructor parameters. Let epsilon=2^-z from small_parameters, Jnear from near_parameters, w from reservoir_parameters, g from block_parameters and b its block count. Set etaNear=Jnear*epsilon, etaCubic=2*2^(-floor(w/2000)), etaGroup=ceil(log2 g)*(2etaNear+etaCubic), EB=(b-1)*(2epsilon+etaGroup). The smaller extractor's proof gives error epsilon at entropy K/2. Iterating it uses entropy at least K/2+floor(K/4), which the raw-source row proof preserves by its explicit rank guard (2g+ceil(log2g))*w<=floor(K/4). Its hybrid therefore gives etaNear. The cubic bound follows from the finite numerical rank and output guards in the cubic lemma; the floored exponent only weakens it. The row mixture induction gives etaGroup by summing the finite per-round error. The full-original-block history induction gives EB by summing each transfer's 2epsilon+etaGroup. These uses require every displayed finite primitive hypothesis, not an eventual-error simplification. The finite proofs use exactly their displayed numerical guards for these inequalities. Their later eventual simplifications of errors are used only in the separate analytic feasibility argument, never to certify a selected finite tuple.

Let EBmax be the maximum over the H levels, R=max2^seed_j and eta=(2R+2)^H*(EBmax+2^-P). At e the hidden row is EBmax-close to uniform on every allowed full-source fiber. Appending correlated rows preserves that approximation by coupling. Going up a level incurs at most R*(2^-P+eta_previous) for selecting the wrong child and adds eta_previous for the inherited challenge. The recurrence is bounded by the displayed eta; its starting value dominates2^-P. Resiliency at parent j decreases by at most m_j,as proved by the variable-width induction below and starts with D_res. Thus all true-path challenges retain sufficient resiliency, and path deviation probability is at most H*R*(2^-P+eta).

Put r_u=w_(a_u), the actual parent response capacity. For a particular possible active stopping index u, the conditional entropy at g gives rare-FAIL lower bounds 2^(-r_u*p1_u) and2^(-r_u*p2_u). Monotonicity gives both tests FAIL with probability at least2^(-r_u*p2_u); the event first FAIL and second PASS has probability at least2^(-r_u*p1_u)-R*(2^-p2_u+eta). Earlier active first-test FAIL events have probability at most sum_(v<u)R*(2^-p1_v+eta). Require for EVERY u that the minimum of the two lower bounds is strictly larger than this sum plus the path deviation bound. Then each output event contains a point whose selected path is correct and whose first active FAIL is g. Both output bits occur on the fixed affine subsource, hence on the original d-dimensional source. This is a finite disperser lemma; its hypotheses are exact rational/integer inequalities.

# Root plateau and tail row budgets with variable active spacing

Write D=w0+X, where X is ROOT spacing. Put y=(L-V)/L and

 g(y)=min(1,sqrt(6y)),
 E_j=min(E0,ceil(h ceil_sqrt(6max(L-V_j,Lpad)) ceil_sqrt(B))),
 E0=ceil(h ceil_sqrt(L) ceil_sqrt(B)),ell_j=2^E_j,ell=2^E0.

Each ordinary binomial cap is ell_j/(L+1). E_j is nonincreasing, so
r_j=max(1,sum_(j<k<H)binom(t_k,b_k))<=ell_j exactly: there are<L
future terms with no larger cap. Global ell*P<=N is retained. The Lpad
floor keeps tail arities large enough; eventually Lpad<L/6. In bulk,
E_j=(h*g(y)+o(1))*sqrt(LB). This profile is constant until y=1/6,
then decreases with sqrt(y), without any future-row maximum mismatch.

Store planned entropy logs W_j. Start a_0=0; from active a_u choose the
FIRST later planned j with V_j>=V_(a_u)+ceil(DZ)-W_(a_u). Repeat with
its own W_j. Last active is H; later planned levels are discarded.
At least one edge is required. Each gap is its target plus at most Amax.
Used arities are capped below the positive ROOT target
ceil(DZ)-ceil(w0Z), the smallest target, so cannot skip a crossing.
The actual retained tree/response/challenge uses these levels.

Backward confidence starts following=Qbase=2^ceil_sqrt(L). At active
parent a_u set p1_u=following,p2_u=4r_(a_u)*p1_u, then
following=4r_(a_u)*p2_u. Require r_(a_u)*p2_u<=P//4 for EVERY u.
output_events uses these full parent response_rows, including padded
rows. Finite rare-FAIL arguments above are unchanged. Confidence
integers have O(LE0+sqrt(L)) bits, even on rejected plans; dyadic
allocations occur only after this guard and have polynomial bit length.

For each fixed rational tuple in the public box require a/wm<.53,
a/w0>.15,2wm>D and p<a/2-.00000102w0. Set G(w)=D-w>0. Define

 F(y)=integral_0^y t/g(t)dt
     =2y^(3/2)/(3sqrt6) for y<=1/6,
     =y^2/2+1/216 for y>=1/6;
 Hfun(y)=integral_0^y g(t)dt
     =2sqrt6*y^(3/2)/3 for y<=1/6,
     =y-1/18 for y>=1/6.

F(1)=109/216,F(1/6)=1/54,Hfun(1)=17/18. For
r(w)=max(0,109/216-4h*integral_w^w0 f), let Y=F^-1(r). Then

 Phi(r)=Hfun(Y)=6r if0<=r<=1/54,
                 sqrt(2r-1/108)-1/18 if r>=1/54.

Phi is increasing and continuously differentiable at1/54 (both
one-sided derivatives6); its first branch is also regular at0.
Define Senergy by two FINITE exact rational certificates:
(1) an original-scale10^12 left-box sum with(864h/109)*sum>1;
(2) a lower rectangle Jlow<=integral_wm^w0 Phi(r(w))/G(w)^2 dw,
constructed by the exact right-suffix/branch/isqrt procedure in the
witness section, satisfying

 2h*[17/(18X)-Jlow]<p.                              (PC)

The energy gives r(wm)=0,r(w0)=109/216 and endpoint Hfun values0,17/18.
All premises are actual finite meshes/integers. Phi(r)/G^2 increases.
Every argument below concerns EACH FIXED tuple with fixed margins.

# Fixed-tuple ordinary tree and retained geometry

Let c_true=(1-beta)/(1+alpha*beta),alpha=1/(1-beta)+.00001,
f(w)=w*c_true(a/w),u_j=L-V_j. Clipped ambients have
L_j=max(u_j,Lpad)+O(1). Wmin^3>=(L+1)^2 and Lpad=ceil(Wmax^(20/19))
put all components in the supplied uniform entropy window. Fixed context
numerators/denominators are O(Z)=o(Lpad); all primitive guards eventually
pass uniformly over these levels. The full half-accuracy near seed fits
the linear cubic output by extractors.md. Least-prime SU estimates give

 b_j=(1/c_true+o(1))*L_j/W_j=Theta(L_j/W_j).

The binomial sandwich (t/b)^b<=binom(t,b)<=t^b gives
A_j<=E_j/b_j+log2 b_j and the feasible lower trial
floor((E_j-log2(L+1)-1)/b_j). Its minimum order
Z sqrt(B)/sqrt(L) grows faster than B. Consequently
A_j=(1+o(1))*E_j/b_j uniformly on used levels, except a terminal cap which
ends planning. Amax=O(Z sqrt(B)/sqrt(Lpad)+B)=o(Z), so the root gap cap is inactive.
All binomial guards and response sizes pass, the latter because p<a/2
leaves a fixed gap to default output exponent a. Unused terminal levels
need no response constructor.

q_j=ceil_log2(2b_j)<=B/4 eventually, since
log2(L/Wmin)=B/4-(log2 B)/4+O(1). Thus
W_j>=max(ceil(wmZ),ceil(w0Z)-jB/4). F is continuously
differentiable across y=1/6. For u_j=L-V_j>=Lpad the two-sided arity
bounds and Taylor's formula give, uniformly,

 F(u_j/L)-F(u_(j+1)/L)
 >=(h-o(1))*f(W_j/Z)*Z sqrt(B)/L^(3/2).             (PP)

Indeed F'=y/g and A_j/L=(h*g(y)+o(1))*f(W_j/Z)*
Z sqrt(B)/(y L^(3/2)); A_j/u_j→0 controls the relative remainder
also on the tail branch. Clipped tail arities are at least
c0 Z sqrt(B)/sqrt(Lpad), so there are only
O(Lpad^(3/2)/(Z sqrt(B)))+1=o(Z/B) clipped steps. If height exceeded
((w0-wm)/(1/4)+tau)Z/B, the first (w0-wm)Z/((1/4)B) raw-schedule
indices would sum PP to at least4h times the certified left integral,
strictly greater than F(1)=109/216. Removing clipped terms and mesh
endpoints loses o(1), because Z^2/(L^(3/2)*sqrt(B))→1.
This is impossible. For every fixed tau>0,

 H<((w0-wm)/(1/4)+tau)Z/B.                          (PH)

Before the entropy floor, sum q_j telescopes exactly into W_s-W_k.
Crossing overshoot is O(B); the floor cannot occur before
(w0-wm)Z/((1/4)B)+O(1) steps, so(PH) leaves at most
(tau+o(1))Z/B floor steps. Uniformly for s<=k<=H,

 log2 product_(s<=j<k)(2b_j)<=W_s-W_k+(tau/4+o(1))Z. (VF)

Each active gap is at most ceil(DZ)-W_parent+Amax, hence its maximum
is<=(D-wm+o(1))Z. Last active H lags the terminal by less than its
next target; include the unused terminal arity to get
log2(c_H+1)<=DZ-W_H+o(Z)<W_H-fixed_margin*Z,
since W_H>=wmZ and2wm>D. W_H means its stored planned value, including
its terminal value if H is the last used depth.

The ORIGINAL geometry dimension remains a separate larger finite guard.
For any Dgeo>w0+D-wm, trial d=2^ceil(Dgeo Z) has
log2 C<=(D-wm+o(1))Z. VF gives all response/hidden margins the positive
leading gap Dgeo-(D-wm)-w0-tau/4, with tau chosen small. Rank/resiliency
costs have exponent p<wm, and arity costs o(Z); the strict leaf margin
pays the last guard. Thus its exact minimum d_plan<=trial<N eventually.
The final2.79666 envelope below bounds WEIGHTED u, not this d_plan.

# Weighted active count and confidence integral

Let G(v)=D-W(v)/Z along the planned schedule. It is nondecreasing,
bounded below by X+o(1); g((L-v)/L) is nonincreasing. On an active gap,
Delta W/Z=O(B/Amin)=o(1), because its length is O(Z), each ordinary
entropy decrement is O(B), and Amin>>B. Thus each weighted left
rectangle, except an endpoint O(sqrt(LB))=o(Z), is bounded by1+o(1)
times its preceding integral rectangle. Ceiling errors are o(Z). Hence

 sum_u log2 r_(a_u)
 <=(h sqrt(LB)/Z)*(1+o(1))*integral_0^L
                 g((L-v)/L)/G(v)dv+o(Z).

As L sqrt(LB)/Z^2→1 this bounds twice the normalized sum by
2h times the limit of integral_0^1 g(y)/[D-w_N(y)]dy.
Fix delta>0 and restrict y>=delta. There q_j/B=1/4+o(1),
b_j=Theta(L/Z), both binomial arity bounds agree relatively, and the
mesh tends to zero. Summing their finite differences gives

 dF(y)/dw=4h f(w),
 F(y)=109/216-4h*integral_w^w0 f+o(1),              (PO)

uniformly on that bulk. These are derived Riemann-sum limits. The floor
cannot be reached at a fixed positive y because witnessed energy would
make F(y)<0. The actual completion w_end>wm is the zero of PO; extension
by zero to wm gives Y(w). The omitted y<delta contributes at most
O(delta^(3/2)/X); the Lpad floor adds o(1). First take N→infinity,
then delta→0. Integration by parts with Hfun'=g therefore gives

 limsup 2sum_u log2 r_(a_u)/Z
 <=2h*integral_0^1 g(y)/[D-w(y)]dy
 =2h*[17/(18X)-integral_wm^w0 Phi(r(w))/(D-w)^2 dw]
 <=2h*[17/(18X)-Jlow]<p.                            (PE)

Phi(r(w)) is absolutely continuous, including its zero and branch
transition, with endpoints0 and17/18 and G>=X>0. Thus the equality
also treats completion before wm. The actually executed LOWER mesh
bounds the subtracted integral in the needed direction.

Backward confidence obeys log2 max(r_u p2_u)<=ceil_sqrt(L)+4hp+
2sum log r_u. Its endpoint terms are o(Z), so PE gives r_u p2_u<=P//4
eventually with fixed margin. Qbase dominates every polynomial in L;
p2_u-r_u p1_u=3r_u p1_u and prior p1_(u-1)-r_u p2_u=3r_u p2_u
therefore dominate log(HR). Earlier p1 are larger. The right error
margin gives EBmax<<2^(-cP) for every fixed c. With R polynomial in N,
eta=(2R+2)^H(EBmax+2^-P)<=2^(-P+O(L^2)). The one-event lower bound
is at least .5*2^(-r_u p1_u)>=2^(-r_u p2_u). Earlier-event/path/eta
losses are negligible using r_u p2_u<=P/4. All exact output margins
are positive. Every finite guard eventually passes at EVERY sufficiently
large integer N for this FIXED tuple. Weighted source selection below
then supplies the smaller claimed envelope.

# One fixed uniform family, NP verifier and infimum

Fix W(N)=0 and V(x,w)=1 iff w is empty and family.evaluate(x)=1.
Every nonempty witness is rejected. Thus J_|x|(x)=1 iff an accepted
witness of length at most W(|x|) exists, in both directions.
The grid enumerates BASE then every six-tuple with denominators<=B,
O(B^12) tuples in fixed order. Every FIXED rational tuple eventually
occurs in every later grid. Finite guards protect all selected tuples;
no uniform threshold over shrinking grid margins is assumed.
For each fixed feasible tuple and D'>w0+X, the WEIGHTED certificate
proved below is<=2^ceil(D'Z) eventually. Therefore minimum-u selection
and FM1 give u_N<=2^ceil(D_energy Z) for every sufficiently large N.
Since Z=o(log N), this yields4N-N^epsilon for each fixed epsilon>0,
and liminf C_B2(J_N)/N>=4 for the same new fixed family.

# Exact active-pair threshold from prefix ranks

At active edge u, each actual child v_j has public inner target I_j>=2;
the edge has outer target O_u>=2. For an affine
source let p_k=H(X_[0,k+1))-H(X_[0,k)), each in{0,1}. After fixing the
prefix before a child [l,r) of a parent ending at t, its useful inner
rank is sum_(l<=k<r)p_k and its outer conditional rank is sum_(r<=k<t)p_k.
Pairs depend only on this profile. Every binary profile is realized
by the span of its selected coordinate vectors. Deleting1s cannot
create a pair, since both sums can only decrease.

For each ACTUAL active interval v define beta(v) as follows. A final
leaf has beta(v)=its actual length. At edge u, let v's ordered children
have values b_0,...,b_(T_u-1). Put o=O_u-1,A_j=min(I_j-1,b_j), and

 beta(v)=max(sum_j A_j,
   max_(j:b_j>=I_j)(sum_(i<j)A_i+b_j+min(o,sum_(i>j)b_i))).       (DP)

An empty inner maximum is omitted. Then beta(v) is EXACTLY the largest
number of pivots in v avoiding every prescribed pair in its subtree.
For the upper bound, each child avoids later pairs and has at most b_j
pivots by induction. If no child reaches its I_j, their total is<=sum A_j.
Otherwise let j be first reaching I_j. Its predecessors cost<=sum_(i<j)A_i,
it costs<=b_j, and its suffix costs<=o to avoid the current pair, as well
as<=sum_(i>j)b_i. This is one term of(DP).

For attainment, choose optimal bad child profiles inductively. The first
term uses A_j pivots in each child, deleting pivots as needed; no child
qualifies. For a j term, use A_i in each predecessor, b_j in child j,
and distribute min(o,sum_(i>j)b_i) pivots among later children, with each
at most its b_i. Subsets of their optimal profiles achieve every such
integer rank. Earlier children are below their I_i; at j and thereafter all
suffixes have at most o pivots. All child profiles avoid deeper pairs.
This attains the term. This proves equality, including empty intervals; the rank-profile identity
covers correlations/offsets. beta(root)+1 is the exact forcing dimension
(vacuous N+1 if beta(root)=N), only for this pair criterion.

For common inner/outer S, no-pair ranks have uniform bounds
B_hp=c_hp, B_u=min(c_u,B_(u+1)+(T_u-1)(S-1)),
c_u=ceil(N/2^V_(a_u)),T_u=2^(V_(a_(u+1))-V_(a_u)).
This follows from(DP): if B_next<S sum<=T_u B_next; otherwise all-low
sum<=T_u(S-1), while a first high child costs at most B_next+(T_u-1)(S-1),
including its suffix if not last. Thus beta(root)+1 is at most
c_hp+1+(S-1)*sum(T_u-1). For original S=ceil(d/C)>c_hp this is
<=1+C(S-1)<=d, supplying every original active-budget use below.

Since2^V_H<=N and arities>=2, the active tree has<=2N-1 nodes.
Prefix/suffix sums evaluate(DP) in linear many polynomial-size integer
operations, uniformly from N, without advice or source rank queries.

For ordinary descent write k=H(parent). Among t ordered children, call
one high when its conditional entropy>=k/(2t). If at least b are high,
any b of them form the needed block-wise affine source: omitting other
preceding coordinates can only increase their conditional entropy.
Otherwise the low children contribute<k/2 in total, so one of the fewer
than b high children has conditional entropy>k/(2b). The b=1 losing case
is impossible. Fix that child's predecessors, then force them to FAIL
before moving again. The parent projection retains at least this child's
entropy, so its local forcing guard holds. Charge responses to Dtot;
dividing by2b scales previously spent rank losses down. Consequently

 F_j=1 for j<=s0, F_(j+1)=2b_jF_j thereafter,
 H(current)>=S/F_j-Dspent,
 hidden entropy>=(S//F_j-Dtot)//(2t_j).

Starting products at s0 is conservative for a later actual start. Every
FAIL persists under later restrictions; the exact hidden guard preserves
the tuple after every rank-D_res conditioning. Primitive, resiliency and
output arguments remain as supplied above.

# Closed exact minimum

Write rank=Dtot,res=D_res and let Q be the maximum of

 rank+max_(j<H)K_j+1;
 F_(j+1)*(rank+K_j+1),for every s0<=j<H;
 F_j*(rank+2t_j*(K_j+res+1)),for every s0<=j<H;
 F_H*(rank+ceil(N/2^V_H)+1).

For integer S, floor(S/F)>=c iff S>=Fc. Applying this also to
the nested hidden floor gives exactly the listed maximum Q; the strict
leaf guard contributes+1. Thus Q is the least sufficient S. The original
relation S=ceil(d/C) gives S>=Q iff d>C(Q-1), so its exact stored minimum
is d(p)=1+C(Q-1). This minimizes the ORIGINAL geometry predicates
for a given plan, not its true dispersion threshold. The controller
rejects d(p)>N and retains every geometry/primitive/error/output
test, then selects least u(p) by(SEL1)-(SEL2). At the selected plan
write d_plan=d(p). The finite proofs apply to each candidate before choice.

# Lower source thresholds within a fixed NEW tuned map

Here SAME compares entropy requirements for one NEW tuned map; it does
not assert equality to a preceding family or its primitive outputs.

At depth j retain tuning K_j and instantiate extractors.md's finite
T(c),(PREFIX),(BLOCK) definitions on the actual accepted tuples.
Put Theta_j=T(response_parameters(n_j,K_j,m_j,mode_j)) and
Hsrc_j=(BLOCK) at block_parameters(n_j,K_j,mode_j), where n_j is the
selected ambient. All tuples are public and already pass their original
guards. The supplied proof gives Theta_j<=K_j,Hsrc_j<=K_j, strong linear
response extraction at Theta_j and block extraction at Hsrc_j, with the
SAME original finite errors. No map, field, parameter or bit changes.

Set rank=Dtot,res=D_res,O=rank+max_(j<H)Theta_j+1. At active edge u let
s=a_(u+1), F^u_j=1 through s and F^u_(j+1)=2b_jF^u_j thereafter.
Define I^0_u as the maximum of

 rank+max_(a_u<j<H)Theta_j+1;
 F^u_(j+1)*(rank+Theta_j+1), for s<=j<H;
 F^u_j*(rank+2t_j*(Hsrc_j+res+1)), for s<=j<H;
 F^u_H*(rank+ceil(N/2^V_H)+1).

This is a sufficient stage target by the supplied finite descent: every
uniform-response guard uses Theta_j; every resilient hidden block uses
Hsrc_j+res+1. Other occurrences of K_j are TUNING parameters and remain
unchanged. Outer O pays ancestor guards at depths<=a_u; later ancestors
use the inner source. The exact ordinary products start at s. Prefixes
preserve outer entropy and every response is charged to the SAME rank
budget. Thus I^0_u,O give the old source route with the refined primitives.

# Actual-node and position-sensitive ordinary budgets

Compute a second ordinary target l(v) at every ACTUAL ordinary node v.
At a leaf set l(v)=length(v)+1. At internal depth j, its child v_i has
actual length n_i. Put q=Hsrc_j+D_res+1 and

 g_i=l(v_0) if i=0;
 g_i=max(l(v_i),Theta_j+1)+i*m_j if i>0;
 U_i=min(n_i,g_i-1),
 bad(v)=sum_i min(U_i,q-1)
        +sum of the largest b_j-1 values among {max(0,U_i-q+1)}_i,
 l(v)=min(length(v)+1,bad(v)+1).                       (POS)

The added sum is0 for b_j=1. bad(v) is the maximum total for ranks<=U_i
with fewer than b_j high ranks>=q. A feasible parent at rank>=l(v) cannot use its capacity branch. If
fewer than b_j children are high and all ranks<=U_i, their sum<=bad(v)
contradicts l(v). Thus either the b high blocks give the resilient row,
or some child has rank>U_i and hence>=g_i by its actual length bound. Fix its predecessors.
For i>0 at most i responses cost i*m_j, during which its parent has
entropy>=g_i-i*m_j>=Theta_j+1. The child retains>=l(v_i). For i=0 there
are no predecessor challenges to force, so no parent-response entropy
guard is needed; simply retain that child. Recursion cannot reach a leaf
or make a losing move at H-1. Every response, full-source resilience rank and error is as in the finite proof.

For a node v at depth s let i_j(v) be its path's child index at depth j,
and set A(v)=sum_(j<min(s,H-1))i_j(v)*m_j. Thus A(v)<=Dtot.
For active edge u and its actual child v at depth s=a_(u+1), put

 M(v)=max({1} union {Theta_j+1:a_u<j<s and i_j(v)>0}),
 J(v)=max(3,A(v)+max(l(v),M(v))), I(v)=min(I^0_u,J(v)). (POS2)

Use this node-specific inner target in(DP), keeping outer O. The old
route applies if I^0_u is smaller. Otherwise force its ancestor siblings
at cost<=A(v). Outer O pays depths<=a_u; at intermediate depths with
positive i_j(v), the retained child pays Theta_j through M(v). With index0
there is nothing to force. The child then retains>=l(v), so(POS) applies.
Costs stay disjoint and total<=Dtot; J(v)>length(v) at leaves excludes useful leaves.

Let t_N=beta(root)+1 from(DP), fallback N+1. As Theta_j,Hsrc_j<=K_j,
every I(v),O<=Q; products starting later are smaller. The common-Q bound
gives t_N<=c_hp+1+(Q-1)(C-1)<=d_plan, since Q>=c_hp+1.
All I(v)>=3; a rank2 coordinate source avoids these pairs, so t_N>=3.
The supplied full-B2 theorem yields4N-8t_N+8 at EVERY length. (DP)
attains a coordinate source of rank t_N-1 avoiding prescribed pairs,
without a monochromaticity claim. Both trees have<=2N-1 nodes;
sorting excesses and threshold arithmetic are polynomial-time from N.
The controller below never queries an unknown source rank.

# Selecting by a weighted proved source certificate

For each candidate passing ALL original finite guards, d_plan=1+C(Q-1)
remains its original geometry dimension. Set c=ceil(N/2^V_H),rank=Dtot,
res=D_res, and use the exact actual (BLOCK) and response entropies
Hsrc_j,Theta_j already defined. Compute Y_H=c+1, then at ALL levels

 Y_j=b_jY_(j+1)+b_j(Theta_j+(t_j-1)m_j)
                         +t_j(Hsrc_j+res)+1.         (SEL1)

POS gives l(v)<=Y_j: each child cap U_i<=Y_next+Theta_j+(t_j-1)m_j,
so bad<=(b_j-1)U+t_j(Hsrc_j+res). Also Y_j>=Y_next. Put

 O=rank+max Theta+1,
 S_u=rank+max(Y_(a_(u+1)),max Theta+1)+3,
 S=rank+max(Y_(a_1),max Theta+1)+3.

POS2 gives every inner I(v)<=J(v)<=S_u at stage u: A(v)<=rank,
l(v)<=Y_childdepth and M(v)<=max Theta+1. Outer target is O<=S_u.
From DP a parent with uniform stage bound S_u has
beta_parent<=B_next+T_u*(S_u-1)+(O-1): all-low costs at most
T_u(S_u-1); otherwise predecessors cost at most that, the first high
child at most B_next, and suffix at most O-1. No monotonic relation
between B_next and S_u is needed. Induction therefore gives

 t(p)<=c+1+sum_u[T_u*(S_u-1)+(O-1)]=weighted(p).

The common-S and original geometry bounds still apply. The program uses

 u(p)=min(d_plan,c+1+(S-1)(C-1),weighted(p)).         (SEL2)

This proves t<=u<=d_plan at EVERY accepted length, even if another
expression exceeds N. source_certificate implements these exact formulas;
tree_parameters chooses least u with first ties among every fully guarded
candidate, without earlier geometry pruning. It defines ONE fixed family
J_N. Fallback sets all three=N+1 and output0. The supplied full-B2 bound
gives C_B2(J_N)>=4N-8t_N+8>=4N-8u_N+8 at EVERY length. Accepted plans
have H<=L and all constructor values<=poly(N); Y_j and added weighted
arithmetic have O(L^2) bits, preserving a fixed polynomial time bound.

For a fixed tuple the response wrapper retunes to
Kresp_j=max(2,ceil_root(m_j^q_j,p_j)), with beta_j=p_j/q_j,
and checks Kresp_j<=K_j before constructing the FULL requested map.
Thus Theta_j=T(actual response)<=Kresp_j exactly. As log m_j=(p+o(1))Z,
log Kresp_j=(p/beta_j+o(1))Z<=(p*w0/a+o(1))Z. Its strict upper
coefficient is<w0/2<wm. This accounts for the NON-polylogarithmic h
in the SU advice: it is never discarded as a poly(log N) factor.
Rank and res have exponent p<wm, and Hsrc_j<=K_j=2^W_j.
Expand SEL1 into a sum
of its nonnegative forcing terms and the terminal term, with intervening
products b_j bounded by VF's products2b_j. The hidden term at k has
logarithm<=W_k+o(Z); rank terms have exponent p and Theta has exponent p*w0/a, both<wm.
The leaf has log(c+1)<W_H-fixed_margin*Z. VF therefore gives,
uniformly in j,

 log2 Y_j<=W_j+(tau/4+o(1))Z.

The sum has only H=poly(log N) terms. Take tau arbitrarily small to
obtain log2 S_u<=W_(a_(u+1))+o(Z). Each active gap has logarithm
log2 T_u<=DZ-W_(a_u)+o(Z), and W_child<=W_parent. Hence every
T_u S_u<=2^((D+o(1))Z); the hp O and leaf terms are smaller. Thus
u_N<=weighted(p)<=2^((D+o(1))Z). In particular for every fixed
D'>D=w0+X, u_N<=2^ceil(D'Z) eventually, although d_plan may be larger.
Each fixed feasible rational tuple eventually occurs in every grid; its
own finite eventual-acceptance proof suffices. No uniform margin over
changing selected tuples, source rank query or huge-tree execution is used.

# Exact response-rank charge through the winning node

At each path depth j at most t_j-1 earlier siblings cost full m_j=w_jP
responses. Prefix coordinate ranks are separate. Thus
Dtot=sum_(j<H-1)(t_j-1)m_j. A useful child is not a leaf by its inner
leaf target. POS cannot lose at H-1, since every child target exceeds
its actual leaf width. In the original 2b route the strict final geometry
inequality gives the same contradiction before forcing there. Hence no
response at H-1 is charged or performed. Earlier FAILs persist.
D_res=sum_(j<H)m_j+max m_j; w_j<=ell gives
Dtot<=H*tmax*ell*P,D_res<=(H+1)ell*P. These are the SAME full-source
rank budgets throughout; no terminal arity is charged.

# Public node-specific response widths and complete rank accounting

Leaves have rho_H=0 rows. Each internal challenge inherits rho_(j+1)
rows and appends binom(t_j,b_j), proving rho_j's recurrence independently
of source values or child choices. A challenge tested at parent j,
including a distant active descendant, has<=rho_(j+1) rows. Set
w_j=max(1,rho_(j+1)). Pad challenges with zero rows and use ONE FULL
m_j=w_j*P-bit linear response at all child/confidence tests, comparing
row prefixes. Outputs are not padded to pretend full uniformity.

The factory builds arities/blocks/products, selects H, checks total rows
<=ell, then constructs full m_j responses only for j<H. Every actual
map passes its finite SU/field/source/seed guards. Targets satisfy
m_j<=ell*P<=N, and accepted/rejected paths retain polynomial bounds.
The local exact challenge-response lemma applies separately at width w_j:
fixing a full uniform affine response costs exactly m_j. All earlier
prefix ranks remain separately charged; both useful entropies lose at
most Dtot from these responses.

Here is the variable-width resilience induction,without importing a global-width output hypothesis. At the internal winning node e,the block row is EB-close to uniform after EVERY additional full-source affine restriction of rank<=D_res,by the hidden entropy guard. Appending arbitrary correlated rows preserves a somewhere-uniform-row coupling. Suppose the true-child challenge at u+1 has error eta_child and remaining conditioning rank r. On any extra fiber of rank<=r-m_u,conditioning on one complete parent response at one seed adds rank<=m_u. The child row is still eta_child-close to uniform there. Matching all challenge rows implies matching its uniform row,with probability<=2^-P+eta_child. Union over the at most R response seeds bounds erroneous child selection by R*(2^-P+eta_child). Inheriting its challenge adds eta_child,so eta_parent<=(2R+1)*eta_child+R*2^-P. Parent block rows and padding are arbitrary appended rows and do not spoil the coupling. Induct upward,spending at most one m_u per level. Starting with D_res=sum_(j<H)m_j+max m_j leaves at least max m_j even after every possible propagation step. Therefore one additional rare-event response conditioning,of rank<=m_g,also fits. The uniform error is bounded by eta=(2R+2)^H*(EB+2^-P),and path deviation by H*R*(2^-P+eta),as before. All restrictions are on the FULL source;no independence of original blocks or uncharged outer projection is assumed.

At a useful active parent g,fixing its useful child coordinates leaves enough outer entropy for the FULL m_g linear response. For each fixed child value some seed is uniform,hence rare-FAIL lower bounds are2^(-w_g*p1),2^(-w_g*p2). These exact bounds are used by output_events with row_widths=(w_(a_u))_u. Upper bounds R*(2^-p+eta) require only one resilient uniform challenge row after the full m_g response conditioning just charged. Both confidence maps are identical,so FAIL monotonicity applies. Earlier active and path errors are bounded exactly as in the primary finite proof. A useful child cannot be a leaf;empty leaf challenges are therefore not used as uniform rows in this argument. Every needed finite inequality remains explicitly checked.

# Uniform runtime

All six parameters/plans depend only on N. C<=1+HN and F_H,d_min have
polynomial bit length; primitive guards bound block counts and dyadic
arithmetic with ell*P<=N. The O(B^12) grid and all accepted/rejected
primitive algorithms share one fixed polynomial time exponent in N.
There is no source advice, witness, oracle or epsilon-selected program.
Auxiliary BASE numerical claims are unused.

# Exact arithmetic and component scope

extractors.md covers the full beta interval. Integer roots use exact
binary search (or the equivalent power-of-two shortcut). Arity search needs O(log L) binomial
comparisons per level, t<=2N,b=O(log N), hence polynomial-size integers.

For midpoint arithmetic write w=U/V,a=p_a/q_a,z=U*q_a-p_a*V. Clearing
positive denominators gives

 f(w)=100000*U*z²/[V*(100000*q_a²*U²+p_a*V*z)].

Choose V=lcm(denominator(wm),denominator(Delta/2),denominator(D)), initial
U=(wm+Delta/2)*V and increment Delta*V. All are exact integers and
(U+j*increment)/V is precisely the jth midpoint. Integer division after
multiplication by10^15 gives the exact lower floor; remainder bounds
and original-rational crosschecks are actually run. No tolerance or
floating integration enters the certificate. For a right endpoint T/V,
use integer division plus1 for the upper f value. For the square-root
rectangle at left endpoint T/V, multiply y_floor*V^2 by scale and
divide by(18*10^12)*(D*V-T)^2>0. All remainders and square bounds
are checked.

# Exact finite challenge-response lemma

This finite lemma permits any public row capacity: instantiate its ell
by the actual w_j and its Pmax by P at each current node.

The following test follows Shaltiel's affine challenge-response method. A matrix distribution is eta-somewhere-random if it is within total variation eta of a convex combination of distributions having at least one uniform row; different mixture components may have different good rows. Joint affinity is not required in this definition. Resiliency R of a challenge relative to an affine source X means the same eta-somewhere-random conclusion holds after fixing any possible value of any binary linear map of X of rank at most R.

At a fixed node j, use its injectively prefix-padded input at ambient
length n_j and call response_parameters(n_j,K_j,ell*Pmax,mode_j).
This wrapper retunes Kresp=max(2,ceil_root((ell*Pmax)^q,p)) for beta=p/q,
requires Kresp<=K_j, and constructs the FULL requested map at Kresp.
A returned tuple c_j supplies EXACTLY ell*Pmax output bits and d_j seed
bits; the finite optional-output theorem proves fixed-seed linearity and
strong error1/4 at entropy Theta_j=T(c_j)<=Kresp<=K_j. Arrange these bits as
ell rows of Pmax columns. All confidences use this SAME full map E_j
and its row prefixes, not a default map trimmed after field selection.
Different nodes may have different maps. Set R=max_(j<H)2^d_j, polynomial
in N by the checked seed guards. Only this common upper response count
is used when combining nodes. Injective embeddings preserve all affine
sources and conditional entropies in the four properties below.

For an input x, coordinate block s, challenge function ch on that block with at most ell rows and Pmax columns, and1<=p<=Pmax, pad ch(x_s) with zero rows to ell rows. Return FAIL exactly when some enumerated response E_j(x,y) agrees with that challenge on the first p columns of every row. Otherwise return PASS. Prefix-zero-pad the actual node input to its specified n_j before the response map. Each test is deterministic polynomial time; challenges are computed by the containing algorithm.

The following properties hold.

1. If p1<=p2, FAIL at p2 implies FAIL at p1, by restricting any matching response to fewer columns.

2. If X is affine with H(X)>=Theta_j and X_s is fixed, there is an affine subsource of EXACT deficiency ell*p on which the test always fails. Strong linear extraction has error1/4, so at least one seed has a completely uniform ell-by-Pmax response on X: a nonuniform affine output would have distance at least1/2. Fix the ell*p specified response coordinates equal to the constant challenge slice. Uniformity makes this a nonempty affine fiber of codimension ell*p, of probability2^(-ell*p). On that fiber one matching response is guaranteed. The argument does not require ch itself to be linear.

3. More generally, if the affine conditional entropy H(X|X_s)>=Theta_j, FAIL has probability at least2^(-ell*p). Fix each possible X_s value. Its fiber remains affine and has entropy at least Theta_j; property2 applies to that fiber. Average the resulting lower bounds. The choice of responding seed may depend on the fiber, which is allowed because FAIL enumerates all seeds.

4. If ch(X_s) is eta-somewhere-random with resiliency ell*p relative to X, then FAIL has probability at most R*(2^-p+eta). For any fixed response seed y, the sliced response E_j is a linear map of X with rank at most ell*p. Condition on any of its possible values v. Resiliency makes the challenge eta-close to a mixture with a uniform row. Matching the whole slice to v then has probability at most2^-p+eta, since matching that component's good row alone has probability2^-p. Average over v, and union bound over all R seeds.

All source restrictions here are affine and have explicitly charged rank deficiency. Property3 is a rare-event lower bound, not a high-probability classification of arbitrary blocks. Property4 is a conditional somewhere-randomness hypothesis, not a claim that every nonempty block automatically passes. The final disperser proof must establish that hypothesis for its recursively chosen challenges and compare its error with the rare FAIL probabilities. Replacing R by n^2 without checking the actual seed length would be unjustified; only its polynomial-in-n bound is used in the intended asymptotic parameter comparisons.


