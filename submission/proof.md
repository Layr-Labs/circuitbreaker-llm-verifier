# Selecting a uniform family by its proved source certificate

H_N=family.evaluate uses ALL N actual input bits. It uses a LINEAR-output cubic primitive, grouping1-beta/2 and
slice beta/2-1/2000000, then selects least source
certificate u(p), ties first, as(SEL1)-(SEL2) below. This differs from the
old minimum-geometry controller. On its selected plan let d_plan=d(p)
and t_N be the exact POS/DP threshold. Then t_N<=u_N<=d_plan and
C_B2(H_N)>=4N-8t_N+8>=4N-8u_N+8 for EVERY N>=0. Fallback sets
all three=N+1 and output0. Quantitative envelopes below concern u_N,
not the old geometry dimension of the newly selected plan. No accepted
large-tree execution or actual optimal circuit-size comparison is claimed.

# Exact midpoint witness for the linear-output reservoir

Define Slin as the six-tuples in the public ranges passing the analytic
prefilter for which SOME positive integer J certifies
8h*Delta*sum_(i=0)^(J-1) floor(10^12*f(wm+i*Delta))/10^12>1,
Delta=(w0-wm)/J. Here

 beta=a/w, alpha=1/(1-beta)+1/100000,
 f(w)=w*(1-beta/2)/(1+4alpha*beta).

The grouping exponent is now1-beta/2. This is a NEW analytic set;
old Slin numerical assertions for a different grouping are not transferred.
Use the fixed rational tuple

 a=283691/500000, X=1070533/1000000,
 wm=535267/500000, w0=635907/200000,
 p=a/4-(102/10^8)*w0-10^-13,
 h=p*X/(2000001/10^6)-10^-13.                         (LM1)

All public-range and strict prefilter inequalities hold by rational
substitution; notably wm-X=10^-6 and a/wm<53/100. Also
w0+X=1062517/250000=4.250068<D_lin=425007/100000=4.25007.

Put e=4/100000,x=w/a>1,y=x-1 and
D=x^2+3x+e(x-1), g(x)=x(x-1)(2x-1)/(2D), so f(w)=a*g(w/a).
Direct differentiation gives

 g'(x)=[2y^4+(20+4e)y^3+(38+3e)y^2+24y+4]/(2D^2),
 g''(x)=[(2e^2+17e+28)y^3+(24e+84)y^2+84y+28-4e]/D^3.

Both are positive for x>1,0<e<1. Also D>=x^2; dropping negative terms
in g' written in x bounds it by1+8/x+3/x^3<=12. Thus f is increasing,
convex and12-Lipschitz. On a midpoint interval, convexity gives
f(mid-t)+f(mid+t)>=2f(mid), so midpoint floors are a LOWER integral bound.

The actual declared audit executes M=65536 midpoint terms, scale10^15:

 sum floor(10^15*f(wm+(i+1/2)*(w0-wm)/M))
      =51160959093340155045.

It checks every integer remainder and four original-Fraction terms.
Writing W=w0-wm, Delta=W/M, the resulting exact ratio is

 R=8h*Delta*sum/10^15
  =2184539449213972839229554826777979266235691
   /2184534425600000000000000000000000000000000>1.

The audit also verifies mu=R-1-8hW/10^12>0 and
J=1+ceil(192hW^2/mu)=28195070. Since f'<=12, the integral minus its
J-term left sum is at most12W^2/J; rounding left values down at10^12
loses at most W/10^12. Therefore the ORIGINAL left-box ratio is at least

 R-8hW/10^12-96hW^2/J>1+mu/2>1.                     (LM2)

This proves(LM1) belongs to Slin without executing J left boxes. No
floating-point integration or optimizer is a premise. The fixed-tuple
feasibility proof below applies at EVERY sufficiently large integer N;
eventual denominator-grid inclusion and minimum-u selection then yield
u_N<=2^ceil(D_lin Z). No accepted large-tree execution is inferred.

# Finite certificate theorem for the current public-grid family

The current family.py implements these conditions. cubic.md and extractors.md supply the complete finite primitive proofs; elimination.md supplies the full unrestricted gate argument. No complete asymptotically accepted tree execution is claimed.

Fix a parameter tuple accepted by all original numerical component guards, the tuple-row guard, and the additional guards ell*P<=N, d<=N and Amax<X. Let the actual interval tree have depth H, products V_j, arities t_j=2^A_j, block counts b_j, component thresholds K_j, active depths a_0=0<...<a_hp=H, fixed-node response maps E_j and confidence pairs(p1_u,p2_u). Define

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

For a particular possible active stopping index u, the conditional entropy at g gives rare-FAIL lower bounds 2^(-ell*p1_u) and2^(-ell*p2_u). Monotonicity gives both tests FAIL with probability at least2^(-ell*p2_u); the event first FAIL and second PASS has probability at least2^(-ell*p1_u)-R*(2^-p2_u+eta). Earlier active first-test FAIL events have probability at most sum_(v<u)R*(2^-p1_v+eta). Require for EVERY u that the minimum of the two lower bounds is strictly larger than this sum plus the path deviation bound. Then each output event contains a point whose selected path is correct and whose first active FAIL is g. Both output bits occur on the fixed affine subsource, hence on the original d-dimensional source. This is a finite disperser lemma; its hypotheses are exact rational/integer inequalities.

# Analytic feasibility and eventual acceptance of each fixed tuple

Let the six rational coordinates h,a,wm,w0,X,p lie in the compact ranges in family.py. Define Slin as EXACTLY the six-tuples in these ranges satisfying the strict prefilter conditions a/wm<.53, a/w0>.15, wm>X, (2.000001)h/X<p<a/4-.00000102w0, and by existence of a finite M-box increasing-potential certificate

 (2h/delta)*Delta*sum_(r=0)^(M-1)floor(10^12*f(wm+rDelta))/10^12>1,
 Delta=(w0-wm)/M,delta=1/4,
 f(w)=w*(1-a/(2w))/(1+4(1/(1-a/w)+.00001)(a/w)).

No numerical integration is a premise. Each FIXED feasible tuple eventually passes all primitive guards by the supplied rational-window proof;its fixed context numerator/denominator times Z are o(Lpad). Current arities use actual public block counts and binomial budgets,as proved below. No threshold uniform over shrinking grid margins is assumed.

For this eventual-feasibility calculation ONLY,fix any rational comparison D>w0+X and put d_trial=2^ceil(DZ). In the following numerical calculation S means ceil(d_trial/C). D is not a public configuration coordinate or a searched value. The exact-minimum lemma below proves the actual chosen geometry dimension is no larger than this trial whenever the trial passes.

The current planner starts W_0=ceil(w0Z). After constructing the actual public b_j at level j,it sets q_j=ceil(log2(2b_j)) and W_(j+1)=max(ceil(wmZ),W_j-q_j). Thus the entropy schedule pays its actual rounded ordinary cost rather than a fixed decline. The unused terminal arity needs no response constructor. The historical config_entropy_log and BASE level_entropy_log are arithmetic witnesses only. Uniform b_j=O(L/Wmin) implies q_j<=B/4 eventually: log2(L/Wmin)=B/4-(log2 B)/4+O(1),whose negative diverging term absorbs the constant and ceiling. Hence W_j>=max(ceil(wmZ),ceil(w0Z)-jB/4),and the potential certificate implies H<((w0-wm)/delta+tau)*L^(3/4)/B^(3/4) for each fixed tau>0. Before the floor is first reached,sum_(r<j)q_r=W_0-W_j,so log2 F_j<=W_0-W_j. At the floor-crossing step the overshoot is at most B/4; the crossing index is at least(w0-wm)Z/(delta B)+O(1). Afterwards the H bound allows at most(tau+o(1))Z/B extra steps. Thus on the floor branch log2 F_j<=(w0-wm+delta*tau+o(1))Z. Choose delta*tau<D-X-w0. Since log2 S=(D-X+o(1))Z,every log2(S/F_j)-W_j has a positive leading margin. Rank budgets have logarithm(p+o(1))Z<wmZ;one step and hidden division cost o(Z);leaf logarithm is(X+o(1))Z<wmZ. All exact geometry predicates follow.

For completeness,the current controller uses the LAST planned active level H. For the last planned level J,the first-crossing definition gives V_J-V_H<Xlog,where Xlog=ceil(X*Z). The terminal planning test gives N<2^(V_J+A_J). Therefore the largest actual H-leaf has length at most2^(Xlog+Amax)+1,and its logarithm is at most(X+o(1))*Z since Amax=o(Z). The strict current prefilter wm>X supplies a positive leaf margin. No k-based stopping threshold is used. The global13Z guard remains an explicit small-length fallback. All primitive applications use the injective clipped ambients and the exact finite source certificate; no threshold is assumed merely from actual leaf length.

The supplied primitive error estimates imply EBmax<=2^(-2^((a/4-.00000102w0-o(1))Z)), below2^-P with a strict margin. Quantized confidence has log(ell^3*p2)<=((2.000001)h/X+o(1))Z<pZ. The original diverging confidence ratio separates each earlier p1 from ell times a later p2. R is polynomial in N; H*log(2R+2)=O((log N)^2)=o(P), and the smallest p2 grows faster than every polynomial in log N. Thus every exact output-event inequality eventually holds. This proves eventual acceptance separately for every FIXED analytically feasible rational tuple.

# One fixed uniform family, NP verifier and infimum

Fix W(N)=0 and V(x,w)=1 iff w is empty and family.evaluate(x)=1.
Every nonempty witness is rejected. Thus H_|x|(x)=1 iff an accepted
witness of length at most W(|x|) exists, in both directions.
The grid enumerates BASE then every six-tuple with denominators<=B,
O(B^12) tuples in fixed order. Every FIXED rational tuple eventually
occurs in every later grid. Finite guards protect all selected tuples;
no uniform threshold over shrinking grid margins is assumed.
For each fixed feasible tuple and D>w0+X, its d(p)<=2^ceil(DZ)
eventually. As u(p)<=d(p), minimum-u selection and then taking infima
give limsup log2(u_N)/Z<=CstarLin=inf_(Slin)(w0+X)<=4.250068.
This is a bound on the new selected source certificate, not its d_plan.
The improved explicit D_lin bound below still has Z=o(log N), yielding
4N-N^epsilon for every fixed epsilon>0 and liminf C_B2(H_N)/N>=4.
The new cubic output and grouping change the selected maps/tree; this
is not an identity claim about either preceding family.

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

# Selecting by a proved source certificate

For each candidate passing ALL original finite guards, d(p)=1+C(Q-1)
is its original geometry dimension; guards still use this d(p).
Let c=ceil(N/2^V_H). The program computes the exact (BLOCK) entropies
Hsrc_j and response entropies Theta_j from their actual constructor
tuples. Set Y_H=c+1 and, backward to s0=a_1,

 Y_j=b_jY_(j+1)+b_j(Theta_j+(t_j-1)m_j)
                    +t_j(Hsrc_j+D_res)+1.             (SEL1)

POS gives l(v)<=Y_j by backward induction: child caps are bounded by
U=Y_next+Theta_j+(t_j-1)m_j; bad<=(b_j-1)U+t_j(Hsrc_j+D_res).
Since b_j>=1 and all
summands are nonnegative, Y_j>=Y_(j+1). Thus every active inner I(v),
outer O, and leaf target c+1 is at most

 S=Dtot+max(Y_s0,max_j Theta_j+1)+3.

The common-S active budget gives t(p)<=c+1+(S-1)(C-1). The original
geometry proof also gives t(p)<=d(p). Hence

 u(p)=min(d(p),c+1+(S-1)(C-1))                       (SEL2)

is a valid source certificate, even if the second expression exceeds N.
source_certificate computes exactly(SEL1)-(SEL2). The new controller
checks every fully guarded accepted candidate and chooses minimum u(p),
ties first. In particular it DOES NOT prune by a previous smaller d(p).
The selected plan determines ONE fixed evaluator H_N. Let u_N=u(p),
d_plan=d(p) and t_N its exact DP threshold. On fallback all three=N+1.
Thus C_B2(H_N)>=4N-8t_N+8>=4N-8u_N+8 at EVERY length.
Each u(p)<=d(p), so minimum-u is bounded by minimum geometry on
THESE changed candidate maps. The new Slin envelopes therefore bound u_N.

For accepted plans H<=L,C<=1+HN,b_j<=2N,Hsrc_j<=K_j<=N,m_j<=N,
so Y_j has O(L^2) bits. New arithmetic and constructors take polynomial
time over O(B^12) candidates. Every numerical rejection and finite hypothesis is checked;
the empty-witness verifier therefore has one fixed polynomial exponent.

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

At a useful active parent g,fixing its useful child coordinates leaves enough outer entropy for the FULL m_g linear response. For each fixed child value some seed is uniform,hence rare-FAIL lower bounds are2^(-w_g*p1),2^(-w_g*p2). They are at least the conservative2^(-ell*p1),2^(-ell*p2) used by unchanged output_events. Upper bounds R*(2^-p+eta) require only one resilient uniform challenge row after the full m_g response conditioning just charged. Both confidence maps are identical,so FAIL monotonicity applies. Earlier active and path errors are bounded exactly as in the primary finite proof. A useful child cannot be a leaf;empty leaf challenges are therefore not used as uniform rows in this argument. Every needed finite inequality remains explicitly checked.

# Uniform runtime

All six parameters/plans depend only on N. C<=1+HN and F_H,d_min have
polynomial bit length; primitive guards bound block counts and dyadic
arithmetic with ell*P<=N. The O(B^12) grid and all accepted/rejected
primitive algorithms share one fixed polynomial time exponent in N.
There is no source advice, witness, oracle or epsilon-selected program.
Auxiliary BASE numerical claims are unused.

# Uniform discrete feasibility details for the current six-set

For a fixed tuple,strict clamps give beta_j=aZ/W_j eventually. Put c_true(beta)=(1-beta/2)/(1+4alpha*beta),alpha=1/(1-beta)+.00001,and f(w)=w*c_true(a/w). This f is increasing. The least prime SU degree and supplied prime asymptotics give its DEFAULT full-output seed d_seed<=(1+4alpha*beta+o(1))*L_j uniformly:the field inequality has logarithm L_j/d_SU+4beta W_j+o(W_j),d_SU=(alpha+o(1))*L_j/W_j,and d_SU tends to infinity even at Lpad. Hence actual b_j=1+ceil(d_seed/log2 g_j)<=(1/c_true+o(1))*L_j/W_j. These estimates hold uniformly over levels of this FIXED tuple.

At every level log2 b=O(B),ell_log grows faster than B,and b>=L_j/O(W_j). The binomial sandwich (t/b)^b<=binom(t,b)<=t^b shows A<=ell_log/b+log2 b and ensures every integer A<=floor((ell_log-log2(L+1)-1)/b) is row-feasible. This lower trial exceeds ceil_log2 b eventually,uniformly. Its upper bound is O(ell_log*Wmax/Lpad+B)=o(Z),with square o(base*L);thus the gap cap is inactive. The remaining floorlog-v+1 cap either produces the unused terminal or loses at most one final step. Consequently at every used nonterminal level A>=(c_true-o(1))*ell_log*W_j/L_j-1. When u_j=L-V_j>=Lpad,L_j=u_j+O(1),so u_j^2-u_(j+1)^2>=2h*f(W_j/Z)*base*L-o(base*L),uniformly. In the clipped range A>=c0*ell_log*Wmin/Lpad for some fixed c0>0,so only O(Lpad^2/(base*L))+1=o(L/base) clipped steps occur. These bounds also prove Amax=o(Z).

Let h0=(w0-wm)/delta and use the fixed M-box witness for this tuple. If J>=(h0+tau)*L^(3/4)/B^(3/4),with fixed tau>0,then the first ceil(h0*Z/B) indices exist eventually. Each box contains Delta*Z/(delta*B)+O(1) raw-schedule indices. Monotonicity lower bounds their f values by the displayed integer floors/10^12. Deleting clipped indices loses o(L/base) bounded terms. Since base*Z/(L*B)->1,summing the arity potential drops gives its normalized witnessed coefficient minus o(1),strictly greater than1,contradicting initial potential L^2. The fixed box rounding errors and summed uniform o(base*L) terms are o(L^2). Hence H<=J<(h0+tau)*L^(3/4)/B^(3/4) eventually.

The planner computes b before A and chooses the largest integer A>=max(1,ceil_log2 b) with A<gap,A<=floorlog-v+1 and binom(2^A,b)<=floor(ell/(L+1)). An unused terminal A may instead be the minimum when v+minimum>floorlog. Every USED level meets the binomial cap. There are at most L used levels,so total rows<=L*floor(ell/(L+1))<ell at every length,independent of asymptotics. No response is built for the unused terminal.

The first active gap is at least ceil(XZ),all gaps are at most ceil(XZ)+Amax,and there are at most L gaps. Hence log2 C=(X+o(1))*Z. The current exact-cost partial-product proof above gives the positive geometry margins for each fixed comparison D>w0+X. In particular the pre-floor margin is D-X-w0,and the post-floor margin is at least D-X-w0-delta*tau. These arguments derive current feasibility directly,without transferring the auxiliary BASE quantitative theorem.

# Exact arithmetic and component scope

The uniform bounds in extractors.md cover the full beta interval and all
five modes. For the root shortcut, n=2^e and d|e gives floor(n^(1/d))=2^(e/d);
other inputs use exact binary search. Arity search needs O(log L) binomial
comparisons per level, t<=2N,b=O(log N), hence polynomial-size integers.

For midpoint arithmetic write w=U/V,a=p_a/q_a,z=U*q_a-p_a*V. Clearing
positive denominators gives

 f(w)=100000*U*z*(2Uq_a-p_aV)/[2V*(100000*q_a²*U²+300000*p_a*q_a*U*V+4*p_a*V*z)].

Choose V=lcm(denominator(wm),denominator(Delta/2)), initial
U=(wm+Delta/2)*V and increment Delta*V. All are exact integers and
(U+j*increment)/V is precisely the jth midpoint. Integer division after
multiplication by10^15 gives the exact lower floor; remainder bounds
and original-rational crosschecks are actually run. No tolerance or
floating integration enters the certificate.

# Exact finite challenge-response lemma

This finite lemma permits any public row capacity: instantiate its ell
by the actual w_j and its Pmax by P at each current node.

The following test follows Shaltiel's affine challenge-response method. A matrix distribution is eta-somewhere-random if it is within total variation eta of a convex combination of distributions having at least one uniform row; different mixture components may have different good rows. Joint affinity is not required in this definition. Resiliency R of a challenge relative to an affine source X means the same eta-somewhere-random conclusion holds after fixing any possible value of any binary linear map of X of rank at most R.

At a fixed node of depth j, take its injectively prefix-padded input at ambient length n_j and its local threshold K_j, within the uniform component window proved below. The supplied constant-error linear extractor has m=floor(K_j^(beta_mode_j)) output bits, d_j=O(log n_j) seed bits and error1/4. Trim to ell by Pmax to obtain E_j. All confidences at this fixed node use the SAME E_j and SAME Pmax. Different nodes at different depths may use different maps. Set R=max_{j<H}2^d_j<=N^6, a common upper bound on every response count. The following four properties apply separately at each node; only the common bound R is used when combining them. Injective linear embeddings preserve every affine source and conditional entropy to which the properties are applied.

For an input x, coordinate block s, challenge function ch on that block with at most ell rows and Pmax columns, and1<=p<=Pmax, pad ch(x_s) with zero rows to ell rows. Return FAIL exactly when some enumerated response E_j(x,y) agrees with that challenge on the first p columns of every row. Otherwise return PASS. Prefix-zero-pad the actual node input to its specified n_j before the response map. Each test is deterministic polynomial time; challenges are computed by the containing algorithm.

The following properties hold.

1. If p1<=p2, FAIL at p2 implies FAIL at p1, by restricting any matching response to fewer columns.

2. If X is affine with H(X)>=K_j and X_s is fixed, there is an affine subsource of EXACT deficiency ell*p on which the test always fails. Strong linear extraction has error1/4, so at least one seed has a completely uniform ell-by-Pmax response on X: a nonuniform affine output would have distance at least1/2. Fix the ell*p specified response coordinates equal to the constant challenge slice. Uniformity makes this a nonempty affine fiber of codimension ell*p, of probability2^(-ell*p). On that fiber one matching response is guaranteed. The argument does not require ch itself to be linear.

3. More generally, if the affine conditional entropy H(X|X_s)>=K_j, FAIL has probability at least2^(-ell*p). Fix each possible X_s value. Its fiber remains affine and has entropy at least K_j; property2 applies to that fiber. Average the resulting lower bounds. The choice of responding seed may depend on the fiber, which is allowed because FAIL enumerates all seeds.

4. If ch(X_s) is eta-somewhere-random with resiliency ell*p relative to X, then FAIL has probability at most R*(2^-p+eta). For any fixed response seed y, the sliced response E_j is a linear map of X with rank at most ell*p. Condition on any of its possible values v. Resiliency makes the challenge eta-close to a mixture with a uniform row. Matching the whole slice to v then has probability at most2^-p+eta, since matching that component's good row alone has probability2^-p. Average over v, and union bound over all R seeds.

All source restrictions here are affine and have explicitly charged rank deficiency. Property3 is a rare-event lower bound, not a high-probability classification of arbitrary blocks. Property4 is a conditional somewhere-randomness hypothesis, not a claim that every nonempty block automatically passes. The final disperser proof must establish that hypothesis for its recursively chosen challenges and compare its error with the rare FAIL probabilities. Replacing R by n^2 without checking the actual seed length would be unjustified; only its polynomial-in-n bound is used in the intended asymptotic parameter comparisons.


