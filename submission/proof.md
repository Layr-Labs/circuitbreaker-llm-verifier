# Selecting a uniform family by its proved source certificate

V_N=family.evaluate uses ALL N actual input bits. It uses a LINEAR-output cubic primitive, accuracy zeta=beta/2-1/1000000, grouping1-beta and
slice beta-1/2000000, then selects least source
certificate u(p), ties first, with depth-dependent row budgets and actual-width backward confidence, as(SEL1)-(SEL2) below. This differs from the
old minimum-geometry controller. On its selected plan let d_plan=d(p)
and t_N be the exact POS/DP threshold. Then t_N<=u_N<=d_plan and
C_B2(V_N)>=4N-8t_N+8>=4N-8u_N+8 for EVERY N>=0. Fallback sets
all three=N+1 and output0. Quantitative envelopes below concern u_N,
not the old geometry dimension of the newly selected plan. No accepted
large-tree execution or actual optimal circuit-size comparison is claimed.

The public ranges are h[.1,.2],a[.4,.7],wm[1.4,2],w0[2.6,3.4],
X[.25,.6],p[.2,.35]; beta remains[.15,.53]. X is root spacing;
D=w0+X and later active gaps grow as the entropy schedule declines.
All maps/finite guards use their actual new parameters.

# Exact energy and weighted-confidence witness

Svar is defined by the two finite integral certificates(VC) below.
Use the fixed rational tuple

 a=257/500,D=16791/5000,w0=14821/5000,
 wm=D/2+10^-6,X=D-w0=197/500,
 p=a/2-(102/10^8)*w0-10^-13,
 h=151515063/10^9.                                  (FM1)

All public range and prefilter inequalities hold exactly, including
2wm>D,a/wm<.53,a/w0>.15 and p<a/2-.00000102w0. Also
w0+X=16791/5000=3.3582<D_var=335821/100000=3.35821.
Put beta=a/w,alpha=1/(1-beta)+.00001 and
f(w)=w*(1-beta)/(1+4alpha*beta).

Put e=4/100000,x=w/a>1,Qpoly=x^2+3x+e(x-1),
g(x)=x(x-1)^2/Qpoly, so f(w)=a*g(w/a). Direct differentiation gives

 g'(x)=(x-1)*[e(2x-1)(x-1)+x^3+7x^2]/Qpoly^2,
 g''(x)=2*[e^2(x-1)^3+3e(x-1)^2(3x+1)+16x^3]/Qpoly^3.

All terms are positive. Since e<1,Qpoly>=x^2, dropping negative terms
in g' bounds it by1+9/x+1/x^3<=11. Hence f is increasing, convex
and11-Lipschitz. Convexity gives f(mid-t)+f(mid+t)>=2f(mid);
integrating each symmetric pair makes floored midpoints a LOWER
integral bound.

The declared audit executes M=65536 midpoint LOWER terms for f and
right-endpoint UPPER terms for f(w)/(D-w), at scale10^15. Since f is
increasing and positive and D>w0, f/(D-w) is increasing. Upward
rounding uses floor(value)+1, a valid strict upper bound even at integers.
All remainders and four original-Fraction terms for EACH sum are checked:

 sum floor(10^15*f(midpoint))=56097677780302274820,
 sum(floor(10^15*f(right)/(D-right))+1)=71278839078432219198.

For Delta=(w0-wm)/65536 the exact energy lower ratio is

 R=6h*Delta*first_sum/10^15
 =1638432443231495508108281207827851
 /1638400000000000000000000000000000>1.

For W=w0-wm,mu=R-1-6hW/10^12>0,
J=1+ceil(132hW^2/mu)=1668012. The integral-minus-left-sum loss is
at most11W^2/J and original-scale rounding loses at most W/10^12, so
R-6hW/10^12-66hW^2/J>1+mu/2>1. This supplies(VC)'s left mesh.
The executed right mesh bounds the confidence integral above by

 8h^2*Delta*second_sum/10^15
 =1051426283230247612817448222391735077569669
 /4096000000000000000000000000000000000000000<p.

Thus both certificates(VC) hold by exact finite arithmetic. No J-box
execution or floating quadrature is a premise. The later fixed-tuple
argument gives u_N<=2^ceil(D_var Z) at EVERY sufficiently large integer
length, not an executed accepted tree or a subsequence claim.

# Finite certificate theorem for the current public-grid family

The current family.py implements these conditions. cubic.md and extractors.md supply the complete finite primitive proofs; elimination.md supplies the full unrestricted gate argument. No complete asymptotically accepted tree execution is claimed.

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

# Profiled row budgets and variable active spacing

Write D=w0+X for this configuration coordinate combination; X is now
ROOT spacing, not every later gap. At ordinary level j put u_j=L-V_j,
E_j=ceil(h ceil_sqrt(max(u_j,Lpad)) ceil_sqrt(B)), ell_j=2^E_j.
The binomial cap is ell_j/(L+1). E_j is nonincreasing and<=E_0, so
r_j=max(1,sum_(j<k<H)binom(t_k,b_k))<=ell_j: there are<L future
terms with no larger cap. Global ell=2^E_0 and ell*P<=N remain guards.
The Lpad floor is essential for the clipped terminal region.

Store the planned entropy logarithms W_j. Start active a_0=0. From a_u
choose the FIRST later planned level with
V_j>=V_(a_u)+ceil(DZ)-W_(a_u). Repeat, using its own W_j. The last such
level is H; unused planned levels are discarded. Require at least one
edge. All gaps are their starting target plus at most Amax. Used arities
are capped below the ROOT target ceil(DZ)-ceil(w0Z), which is positive
eventually and the smallest target. Thus they cannot skip a crossing.
Every finite accepted tree/response/challenge uses these actual levels.

Backward confidence uses Qbase=2^ceil_sqrt(L), starting following=Qbase.
At active parent a_u set p1_u=following,p2_u=4r_(a_u)*p1_u, then
following=4r_(a_u)*p2_u. Require r_(a_u)*p2_u<=P//4 for every u.
output_events uses EXACTLY these full parent response_rows, including
padded rows. The finite rare-FAIL proof above therefore applies without
change. Confidence integers have O(L E_0+sqrt(L)) bits even for rejected
plans; dyadic denominators are formed only after their displayed guard.

Define Svar as fixed rational tuples in the public compact ranges with
a/wm<.53,a/w0>.15,2wm>D, and p<a/2-.00000102w0, admitting finite
positive-mesh certificates

 6h*lower_left_integral(f)>1,
 8h^2*upper_right_integral(f(w)/(D-w))<p.             (VC)

Both use exact rational floors: left values rounded down at10^12,
right values rounded UP at10^15. D>w0 since X>0. No numerical
quadrature or optimizer is a premise. The next sections prove eventual
acceptance and the NEW weighted source envelope for each FIXED tuple.

# Fixed-tuple ordinary tree and geometry

Put c_true=(1-beta)/(1+4alpha*beta),alpha=1/(1-beta)+.00001,
f(w)=w*c_true(a/w). Clipped ambients have L_j=max(u_j,Lpad)+O(1).
The Wmin^3>=(L+1)^2 and Lpad=ceil(Wmax^(20/19)) guards put every
component in its supplied entropy window. Fixed context denominators
and numerators are O(Z)=o(Lpad); all primitive guards eventually pass
uniformly over this tuple's levels. The full half-accuracy near seed
fits the linear cubic output by extractors.md's strict exponent gap.
The least-prime SU seed gives, uniformly,

 b_j=(1/c_true+o(1))*L_j/W_j=Theta(L_j/W_j).

The inequalities (t/b)^b<=binom(t,b)<=t^b give
A_j<=E_j/b_j+log2 b_j and the feasible lower trial
floor((E_j-log2(L+1)-1)/b_j). Its minimum has order
Z sqrt(B)/sqrt(L), which grows faster than B. Thus log b and rounding
are negligible, and A_j=(1+o(1))*E_j/b_j uniformly at used levels,
except a terminal cap which ends the plan. Also
Amax=O(Z sqrt(B)/sqrt(Lpad)+B)=o(Z); the root gap cap is inactive.
All row guards and response sizes pass; p<a/2 leaves a positive margin
to default output exponent a. No response is made for an unused terminal.

q_j=ceil_log2(2b_j)<=B/4 eventually, since
log2(L/Wmin)=B/4-(log2 B)/4+O(1). Consequently
W_j>=max(ceil(wmZ),ceil(w0Z)-jB/4). For u_j>=Lpad, Taylor's formula gives

 u_j^(3/2)-u_(j+1)^(3/2)
 >=(3h/2-o(1))*f(W_j/Z)*Z sqrt(B).                   (VP)

A_j/u_j tends uniformly to zero; its quadratic Taylor error is relative
o(1). In the clipped region A_j>=c0 Z sqrt(B)/sqrt(Lpad), so it has
O(Lpad^(3/2)/(Z sqrt(B)))+1=o(Z/B) steps. If height exceeded
((w0-wm)/(1/4)+tau)Z/B, the first (w0-wm)Z/((1/4)B) raw-schedule
indices, split into the fixed left boxes from(VC), would sum(VP) to
more than L^(3/2): its normalized coefficient is6h times that left sum.
Deleting clipped indices and endpoint rounding lose o(L^(3/2)), since
Z^2/sqrt(B)~L^(3/2). This is impossible. Hence, for every fixed tau>0,

 H<((w0-wm)/(1/4)+tau)Z/B.                          (VH)

For ANY ordinary start s, before the entropy floor, sum_(s<=j<k)q_j
is W_s-W_k. Crossing overshoot is O(B). The floor cannot be reached
before (w0-wm)Z/((1/4)B)+O(1) steps, so(VH) leaves at most
(tau+o(1))Z/B floor steps. Therefore, uniformly in s<=k<=H,

 log2 product_(s<=j<k)(2b_j)<=W_s-W_k+(tau/4+o(1))Z. (VF)

All active gaps are between the root target and
ceil(DZ)-W_(a_u)+Amax; their maximum is<=(D-wm+o(1))Z.
The last active H lags the terminal by less than its next target;
including the unused terminal arity gives
log2(c_H+1)<=D Z-W_H+o(Z)<W_H-fixed_margin*Z,
since W_H>=wmZ and2wm>D. Here W_H is the stored planned entropy,
including its terminal value if H is the last used depth.

Original geometry is retained as a FINITE guard, with its own larger
d_plan. For any Dgeo>w0+D-wm, trial d=2^ceil(Dgeo Z) has
log2 C<=(D-wm+o(1))Z. Use(VF) from s0 and choose tau sufficiently
small. Every response/hidden geometry margin has positive leading gap
Dgeo-(D-wm)-w0-tau/4. Rank and resiliency costs have logarithm
(p+o(1))Z<wmZ; one arity costs o(Z). The strict leaf margin just proved
pays the final geometry guard. Thus its exact minimum d_plan<=trial<N
eventually. This does NOT claim d_plan has the smaller final envelope.

# Weighted confidence estimate

Let F(v)=sqrt(max(L-v,Lpad)) and G(v)=D-W(v)/Z on the planned
piecewise-constant entropy schedule. F is nonincreasing; G is
nondecreasing and bounded below by X+o(1)>0. On one active gap,
Delta W/Z=O(B/Amin)=o(1), uniformly, since each entropy decrement is
O(B), each ordinary step is at least Amin>>B and gap length is O(Z).
Thus adjacent G ratios tend to1. Each active left rectangle, except
an endpoint of size O(sqrt(LB))=o(Z), is bounded by1+o(1) times the
preceding integral rectangle for F/G. The local E_j rounding contributes
o(Z). Consequently

 sum_u log2 r_(a_u)
 <=(h sqrt(B)/Z)*(1+o(1))*integral_0^L F(v)/G(v) dv+o(Z). (VR)

To identify this integral, fix any delta>0 and restrict u=L-v>=delta L.
Here q_j/B=1/4+o(1), uniformly: b_j=Theta(L/Z) with constants depending
on delta. The two-sided arity estimate above, with vanishing relative
rounding errors, gives for z=(u/L)^(3/2),w=W/Z,

 Delta z/(-Delta w)=6h f(w)+o(1),                  (VO)

until the entropy floor. The mesh size tends to zero. Summing these
uniform estimates gives the corresponding integral relation, or its
Riemann-Stieltjes form dz=6h f(w)dw while w decreases. The floor cannot
be reached with u>=delta L for arbitrarily small fixed delta: otherwise
its accumulated decrease is6h integral_wm^w0 f+o(1)>1 by(VC).
Thus the floor region, if any, contributes only to the tail as delta→0.
Changing variables in(VR), sqrt(u)du=(2/3)L^(3/2)dz and
Z^2~L^(3/2)sqrt(B), yields

 limsup 2sum_u log2 r_(a_u)/Z
 <=8h^2 integral_wm^w0 f(w)/(D-w) dw <p.             (VE)

An actual completion before wm only shortens this nonnegative integral.
The omitted u<delta L part is O(delta^(3/2)/X) in the normalized bound;
the Lpad floor adds o(1). First take N→infinity, then delta→0.
This justifies both the variable-spacing Riemann sum and its constant.

The confidence recursion has log2 max(r_u p2_u)<=
ceil_sqrt(L)+4hp+2sum log2 r_u. Its endpoint terms are o(Z), so(VE)
gives r_u p2_u<=P//4 eventually with a fixed margin. Qbase grows faster
than every polynomial in L. Thus p2_u-r_u p1_u=3r_u p1_u and
p1_(u-1)-r_u p2_u=3r_u p2_u dominate log(HR); earlier p1 are larger.
The right error margin gives EBmax<<2^(-cP) for every fixed c.
eta=(2R+2)^H(EBmax+2^-P)<=2^(-P+O(L^2)), with R polynomial in N.
The one-event bound is at least .5*2^(-r_u p1_u)>=2^(-r_u p2_u).
All earlier-event/path/eta losses are negligible relative to the latter,
using r_u p2_u<=P/4. Every exact output_events margin is positive.
All primitive, row, geometry and confidence guards therefore pass at
EVERY sufficiently large integer N for this FIXED rational tuple.

# One fixed uniform family, NP verifier and infimum

Fix W(N)=0 and V(x,w)=1 iff w is empty and family.evaluate(x)=1.
Every nonempty witness is rejected. Thus V_|x|(x)=1 iff an accepted
witness of length at most W(|x|) exists, in both directions.
The grid enumerates BASE then every six-tuple with denominators<=B,
O(B^12) tuples in fixed order. Every FIXED rational tuple eventually
occurs in every later grid. Finite guards protect all selected tuples;
no uniform threshold over shrinking grid margins is assumed.
For each fixed feasible tuple and D'>w0+X, the WEIGHTED certificate
proved below is<=2^ceil(D'Z) eventually. Therefore minimum-u selection
and FM1 give u_N<=2^ceil(D_var Z) for every sufficiently large N.
Since Z=o(log N), this yields4N-N^epsilon for each fixed epsilon>0,
and liminf C_B2(V_N)/N>=4 for the same new fixed family.

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
V_N. Fallback sets all three=N+1 and output0. The supplied full-B2 bound
gives C_B2(V_N)>=4N-8t_N+8>=4N-8u_N+8 at EVERY length. Accepted plans
have H<=L and all constructor values<=poly(N); Y_j and added weighted
arithmetic have O(L^2) bits, preserving a fixed polynomial time bound.

For a fixed feasible tuple, Theta_j=T(actual requested response) has
log2 Theta_j<=(p+o(1))Z: its m_j=r_jP has that bound, and all other
factors400hds^2 and the additive log term in T are poly(log N), from
the explicit finite constructors and fixed rational contexts. Rank and
res also have exponent p. Hsrc_j<=K_j=2^W_j. Expand SEL1 into a sum
of its nonnegative forcing terms and the terminal term, with intervening
products b_j bounded by VF's products2b_j. The hidden term at k has
logarithm<=W_k+o(Z); rank/Theta terms are smaller since p<wm.
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

The uniform bounds in extractors.md cover the full beta interval and all
five modes. For the root shortcut, n=2^e and d|e gives floor(n^(1/d))=2^(e/d);
other inputs use exact binary search. Arity search needs O(log L) binomial
comparisons per level, t<=2N,b=O(log N), hence polynomial-size integers.

For midpoint arithmetic write w=U/V,a=p_a/q_a,z=U*q_a-p_a*V. Clearing
positive denominators gives

 f(w)=100000*U*z²/[V*(100000*q_a²*U²+300000*p_a*q_a*U*V+4*p_a*V*z)].

Choose V=lcm(denominator(wm),denominator(Delta/2),denominator(D)), initial
U=(wm+Delta/2)*V and increment Delta*V. All are exact integers and
(U+j*increment)/V is precisely the jth midpoint. Integer division after
multiplication by10^15 gives the exact lower floor; remainder bounds
and original-rational crosschecks are actually run. No tolerance or
floating integration enters the certificate. For a right endpoint T/V,
multiply its f numerator by V and denominator by (D*V-T)>0,
then use integer division plus1. The audit checks those remainders
and the original Fraction expression too.

# Exact finite challenge-response lemma

This finite lemma permits any public row capacity: instantiate its ell
by the actual w_j and its Pmax by P at each current node.

The following test follows Shaltiel's affine challenge-response method. A matrix distribution is eta-somewhere-random if it is within total variation eta of a convex combination of distributions having at least one uniform row; different mixture components may have different good rows. Joint affinity is not required in this definition. Resiliency R of a challenge relative to an affine source X means the same eta-somewhere-random conclusion holds after fixing any possible value of any binary linear map of X of rank at most R.

At a fixed node j, use its injectively prefix-padded input at ambient
length n_j and call response_parameters(n_j,K_j,ell*Pmax,mode_j).
A returned tuple c_j supplies EXACTLY ell*Pmax output bits and d_j seed
bits; the finite optional-output theorem proves fixed-seed linearity and
strong error1/4 at entropy Theta_j=T(c_j)<=K_j. Arrange these bits as
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


