# Additive8 and the same exact public-grid tree family

F_N=family.evaluate is the unchanged source85 total public-grid language
on ALL actual N bits. Its selected exact dimension is d_N; fallback is
N+1. This package proves C_B2(F_N)>=4N-8d_N+8 for EVERY N>=0, and the
retained explicit midpoint envelope eventually at EVERY integer length.
The leading projection is(4,1). Tree and primitive constructions are
inherited; the supplied generic additive8 refinement and its exact finite
obligations are transferred to this same evaluator. No prior verdict or
external solver status is a premise.

For N=0 fallback gives bound0; positive fallback gives-4N. At acceptance,
the full finite source proof below gives dispersion at d_N. The strict
leaf guard forces S>=2. There are at least two active depths, so
C=1+sum_u(2^(V_(a_(u+1))-V_(a_u))-1)>=2. Also every W_j>=1, hence
K_j=2^W_j>=2 and Q>=rank+max_j K_j+1>=3. Thus d_N=1+C*(Q-1)>=5.
Apply the complete supplied additive8 theorem for d>=4 at n=N.
No bits are fixed or omitted. The asymptotic guarantee follows from
the retained fixed-tuple feasibility and exact midpoint witness below.

All finite formula/Ramsey obligations are actually checked. Audit evidence
is scoped to those finite predicates and the rational witness, not a
complete accepted large-tree execution. Historical BASE numerical claims
and comparisons are omitted. The supplied current finite-source,
primitive, geometry and uniform-runtime proofs establish the declaration.

# A stronger exact convex-midpoint certificate for the SAME controller

No evaluator or geometry guard changes. Write W=w0-wm. The rational tuple

 a=550791690219/10^12, X=104292796651/10^11,
 wm=X+10^-12, w0=1683067911/500000000,
 p=a/4-(102/10^8)*w0-10^-13,
 h=p*X/(2000001/10^6)-10^-13                         (M1)

lies in the existing public ranges and passes every strict analytic
prefilter. These are exact rational definitions. In particular wm>X,
a/wm<53/100, a/w0>3/20, (2000001/10^6)h/X<p,
and p<a/4-(102/10^8)w0. Its sum is

 w0+X=440906378851/10^11=4.40906378851.

Here is a finite certificate of membership in the SAME analytic set S6
used in the retained proof, even though that set is stated using left
boxes with scale10^12. Put e=4/100000, x=w/a>1 and

 g(x)=x*(x-1)^2/[x*(x+3)+e*(x-1)], f(w)=a*g(w/a).

This is exactly the retained potential f, by clearing its denominator.
Direct polynomial differentiation, with D=x²+3x+e*(x-1), gives

 g'(x)=(x-1)*[e*(2x-1)*(x-1)+x³+7x²]/D²,
 g''(x)=2*[e²*(x-1)³+3e*(x-1)²*(3x+1)+16x³]/D³.

All displayed denominators and numerators are positive for x>1. Thus
f is increasing and convex on [wm,w0]. Also 0<f'(w)<=11: since e<1
and D>=x², drop the negative terms in the numerator of g' to bound it
by x*(x³+9x²+1)/x^4=1+9/x+1/x³<=11.

For each subinterval of width Delta, convexity gives
f(mid-t)+f(mid+t)>=2f(mid). Integrating for0<=t<=Delta/2
shows its integral is at least Delta*f(mid). Consequently its exact
floored midpoint sum is a certified LOWER bound, not a quadrature guess.
At M=65536, Delta=W/M, scale S=10^15, the actual declared audit computes

 sum_(j=0)^(M-1) floor(S*f(wm+(j+1/2)*Delta))
       =49109014265147760617.

The integer formula is the retained common-denominator formula, evaluated
at these exact midpoints. It crosschecks selected terms with the original
rational expression. The resulting R=8h*Delta*sum/S is exactly

 40960020509752098178821891163443417586661219884120373286607
 /40960020480000000000000000000000000000000000000000000000000

and is greater than1. The declared audit checks this exact rational
identity, both strict prefilter inequalities and

 mu=R-1-8h*W/10^12>0.                                (M2)

To connect to the EXISTING left-box definition without changing S6,
let J=1+ceil(176*h*W²/mu). On each J-mesh interval the bound f'<=11
makes the integral minus the left sum at most11*W²/J overall.
Rounding each left value down at the ORIGINAL scale10^12 loses at most
W/10^12. Since the integral is at least R/(8h), its original left-box
ratio is therefore at least

 R-8h*W/10^12-88*h*W²/J >1+mu/2>1.

Thus this explicit finite integer J certifies membership in the original
set. No J-term execution is needed: the displayed deterministic error
inequality proves its finite left-box certificate from the actual
midpoint sum. No numerical integration, floating optimizer or limit
swap is a premise. The audit verifies J's exact integer inequality.

Apply the retained fixed-tuple feasibility theorem to(M1). Each fixed
rational coordinate, however large its denominator, eventually lies in
the unchanged denominator-cap grid. For every fixed rational D>w0+X,
the SAME minimum-certificate controller eventually has d_N<=2^ceil(DZ).
In particular take

 D_new=8818127579/2000000000=4.4090637895>w0+X.

Apply the supplied gate theorem after the finite source proof. Its
explicit D_new envelope holds at every sufficiently large integer length,
for the identical evaluator. No practical large-tree run or optimality
of this rational witness is claimed.

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

Start with an arbitrary uniform affine source of entropy at least d. The exact integer active-budget lemma proved below applies with C=1+sum(T_u-1),S=ceil(d/C) and targets K_u=S+(S-1)*sum_(r>=u)(T_r-1). The root target is at most d,and the ending target S exceeds the largest leaf by the assumed last geometry inequality. The lemma fixes each preceding prefix at its charged conditional rank and gives consecutive active nodes g,g' with BOTH useful entropies at least S,and coordinates before g' fixed. Its last-child case and integer entropy inequalities are supplied below;there is no multiplicative active entropy-decay assumption.

Every subsequent prefix fixing inside g' preserves H(X_g|X_g') exactly by the affine fiber dimension identity; every rank-r response fixing decreases it by at most r. Force earlier siblings on ancestors of g' to fail at width P, using the exact challenge-response deficiency lemma while the parent's entropy is at least S-Dtot>maxK. The sibling's coordinates are fixed, so its challenge is constant. A full uniform response exists because a strong linear extractor of error1/4 has some seed with uniform affine output: a nonuniform affine output has distance at least1/2. Fix m_j=w_j*P response coordinates at that parent. This costs exactly that rank and leaves all previously forced events true pointwise.

Descend ordinarily from g'. Every possible g' has depth at least s0=a_1. At ordinary depth j>=s0,use the product F_j that starts at this first possible descent level,even when the actual descent starts later. Retain a nice source with entropy at least S/F_j-Dspent. The binary ordinary-descent lemma below uses high-child threshold parent_entropy/(2t_j). Either it finds b_j hidden blocks or it keeps a child at entropy greater than parent_entropy/(2b_j). Fix its predecessors and force them to FAIL BEFORE the next move,while the parent projection is at least the retained child and hence at least S//F_(j+1)-Dspent>K_j. Here F_(j+1)=2b_jF_j. Charge all response constraints to Dtot;division by2b_j scales old rank losses down. Every previous deterministic FAIL persists pointwise,and later prefixes inside g' preserve the outer conditional entropy. The complete winning/losing entropy calculation is supplied below.

The procedure cannot reach a leaf, by the final geometry condition. At its stopping internal node e of depth j, every hidden block has entropy at least (S//F_j-Dtot)//(2t_j)>=K_j+D_res+1. It therefore remains a block source above K_j after every possible extra rank-D_res conditioning of the full source. The conditional entropy at g given g' remains at least S-Dtot>K_g. All earlier-sibling tests on the path to e fail identically. No property is assumed about retaining a large earlier threshold after all later restrictions.

The pointwise primitive error used by the controller is obtained from actual constructor parameters. Let epsilon=2^-z from small_parameters, Jnear from near_parameters, w from reservoir_parameters, g from block_parameters and b its block count. Set etaNear=Jnear*epsilon, etaCubic=2*2^(-floor(isqrt(w)/50)), etaGroup=ceil(log2 g)*(2etaNear+etaCubic), EB=(b-1)*(2epsilon+etaGroup). The smaller extractor's proof gives error epsilon at entropy K/2. Iterating it uses entropy at least K/2+floor(K/4), which the raw-source row proof preserves by its explicit rank guard (2g+ceil(log2g))*w<=floor(K/4). Its hybrid therefore gives etaNear. The cubic bound follows from the finite numerical rank and output guards in the cubic lemma; the floored exponent only weakens it. The row mixture induction gives etaGroup by summing the finite per-round error. The full-original-block history induction gives EB by summing each transfer's 2epsilon+etaGroup. These uses require every displayed finite primitive hypothesis, not an eventual-error simplification. The finite proofs use exactly their displayed numerical guards for these inequalities. Their later eventual simplifications of errors are used only in the separate analytic feasibility argument, never to certify a selected finite tuple.

Let EBmax be the maximum over the H levels, R=max2^seed_j and eta=(2R+2)^H*(EBmax+2^-P). At e the hidden row is EBmax-close to uniform on every allowed full-source fiber. Appending correlated rows preserves that approximation by coupling. Going up a level incurs at most R*(2^-P+eta_previous) for selecting the wrong child and adds eta_previous for the inherited challenge. The recurrence is bounded by the displayed eta; its starting value dominates2^-P. Resiliency at parent j decreases by at most m_j,as proved by the variable-width induction below and starts with D_res. Thus all true-path challenges retain sufficient resiliency, and path deviation probability is at most H*R*(2^-P+eta).

For a particular possible active stopping index u, the conditional entropy at g gives rare-FAIL lower bounds 2^(-ell*p1_u) and2^(-ell*p2_u). Monotonicity gives both tests FAIL with probability at least2^(-ell*p2_u); the event first FAIL and second PASS has probability at least2^(-ell*p1_u)-R*(2^-p2_u+eta). Earlier active first-test FAIL events have probability at most sum_(v<u)R*(2^-p1_v+eta). Require for EVERY u that the minimum of the two lower bounds is strictly larger than this sum plus the path deviation bound. Then each output event contains a point whose selected path is correct and whose first active FAIL is g. Both output bits occur on the fixed affine subsource, hence on the original d-dimensional source. This is a finite disperser lemma; its hypotheses are exact rational/integer inequalities.

# Analytic feasibility and eventual acceptance of each fixed tuple

Let the six rational coordinates h,a,wm,w0,X,p lie in the compact ranges in family.py. Define S6 as EXACTLY the six-tuples in these ranges satisfying the strict prefilter conditions a/wm<.53, a/w0>.15, wm>X, (2.000001)h/X<p<a/4-.00000102w0, and by existence of a finite M-box increasing-potential certificate

 (2h/delta)*Delta*sum_(r=0)^(M-1)floor(10^12*f(wm+rDelta))/10^12>1,
 Delta=(w0-wm)/M,delta=1/4,
 f(w)=w*(1-a/w)/(1+4(1/(1-a/w)+.00001)(a/w)).

No numerical integration is a premise. Each FIXED feasible tuple eventually passes all primitive guards by the supplied rational-window proof;its fixed context numerator/denominator times Z are o(Lpad). Current arities use actual public block counts and binomial budgets,as proved below. No threshold uniform over shrinking grid margins is assumed.

For this eventual-feasibility calculation ONLY,fix any rational comparison D>w0+X and put d_trial=2^ceil(DZ). In the following numerical calculation S means ceil(d_trial/C). D is not a public configuration coordinate or a searched value. The exact-minimum lemma below proves the actual chosen geometry dimension is no larger than this trial whenever the trial passes.

The current planner starts W_0=ceil(w0Z). After constructing the actual public b_j at level j,it sets q_j=ceil(log2(2b_j)) and W_(j+1)=max(ceil(wmZ),W_j-q_j). Thus the entropy schedule pays its actual rounded ordinary cost rather than a fixed decline. The unused terminal arity needs no response constructor. The historical config_entropy_log and BASE level_entropy_log are arithmetic witnesses only. Uniform b_j=O(L/Wmin) implies q_j<=B/4 eventually: log2(L/Wmin)=B/4-(log2 B)/4+O(1),whose negative diverging term absorbs the constant and ceiling. Hence W_j>=max(ceil(wmZ),ceil(w0Z)-jB/4),and the potential certificate implies H<((w0-wm)/delta+tau)*L^(3/4)/B^(3/4) for each fixed tau>0. Before the floor is first reached,sum_(r<j)q_r=W_0-W_j,so log2 F_j<=W_0-W_j. At the floor-crossing step the overshoot is at most B/4; the crossing index is at least(w0-wm)Z/(delta B)+O(1). Afterwards the H bound allows at most(tau+o(1))Z/B extra steps. Thus on the floor branch log2 F_j<=(w0-wm+delta*tau+o(1))Z. Choose delta*tau<D-X-w0. Since log2 S=(D-X+o(1))Z,every log2(S/F_j)-W_j has a positive leading margin. Rank budgets have logarithm(p+o(1))Z<wmZ;one step and hidden division cost o(Z);leaf logarithm is(X+o(1))Z<wmZ. All exact geometry predicates follow.

For completeness,the current controller uses the LAST planned active level H. For the last planned level J,the first-crossing definition gives V_J-V_H<Xlog,where Xlog=ceil(X*Z). The terminal planning test gives N<2^(V_J+A_J). Therefore the largest actual H-leaf has length at most2^(Xlog+Amax)+1,and its logarithm is at most(X+o(1))*Z since Amax=o(Z). The strict current prefilter wm>X supplies a positive leaf margin. No k-based stopping threshold is used. The global13Z guard remains an explicit small-length fallback. All primitive applications use the injective clipped ambients and the exact finite source certificate; no threshold is assumed merely from actual leaf length.

The supplied primitive error estimates imply EBmax<=2^(-2^((a/4-.00000102w0-o(1))Z)), below2^-P with a strict margin. Quantized confidence has log(ell^3*p2)<=((2.000001)h/X+o(1))Z<pZ. The original diverging confidence ratio separates each earlier p1 from ell times a later p2. R is polynomial in N; H*log(2R+2)=O((log N)^2)=o(P), and the smallest p2 grows faster than every polynomial in log N. Thus every exact output-event inequality eventually holds. This proves eventual acceptance separately for every FIXED analytically feasible rational tuple.

# One fixed uniform family, NP verifier and infimum

Fix W(N)=0 and V(x,w)=1 iff w is empty and family.evaluate(x)=1.
Every nonempty witness is rejected. Thus F_|x|(x)=1 iff an accepted
witness of length at most W(|x|) exists, in both directions.

The controller enumerates BASE then every six-tuple with coordinate
denominators<=B in fixed order, retaining smallest accepted d, ties first.
There are O(B^12) tuples. Every FIXED rational tuple eventually occurs
in every later grid. Finite guards protect all length-selected tuples;
no uniform threshold over shrinking grid margins is assumed.
For each fixed feasible tuple and fixed rational D>w0+X, the proof above
makes the actual minimum at most2^ceil(DZ) eventually. Taking D down
to w0+X and then infima proves limsup log2(d_N)/Z<=Cstar6=inf_(S6)(w0+X).
The exact midpoint witness gives Cstar6<=4.40906378851 and the explicit
eventual D_new bound. Since Z=o(log N), every fixed epsilon>0 eventually
has8d_N<=N^epsilon; the SAME program gives4N-N^epsilon.

# Exact integer active budgets and binary ordinary descent

The current controller has six public coordinates and no searched D. Its dimension certificate uses integers throughout. Let T_u be the actual number of descendants between consecutive active depths,and define

 C=1+sum_(u<hp)(T_u-1), S=ceil(d/C).

The last geometry inequality implies S>leaf_cap=ceil(N/2^V_H). Set

 K_u=S+(S-1)*sum_(r=u)^(hp-1)(T_r-1).

Then K_hp=S exceeds every actual final leaf length. Also K_0=1+C(S-1)<=d by the defining ceiling. All these are public proof quantities; the evaluator does not query an input source or its entropies.

Suppose at active interval u the nice parent has entropy at least K_u. Conditional entropies of its T_u ordered children sum to that entropy. There is a first child with conditional entropy at least S: if all were below S,their integer entropies would sum to at most T_u(S-1),while K_u>=T_u(S-1)+1. If this is child j,fix its j predecessors at a possible value. Their total rank cost is at most j(S-1). The selected child now has entropy at least S.

If its entropy is at least K_(u+1),continue there. Otherwise it cannot be the LAST child. In that case j=T_u-1 and the fixed-prefix parent projection equals the child projection,with entropy at least K_u-(T_u-1)(S-1)=K_(u+1),contradicting the stopping condition. Consequently j<=T_u-2,and child entropy is at most K_(u+1)-1. The parent conditional entropy given it is at least

 K_u-j(S-1)-(K_(u+1)-1)
   =(T_u-1-j)(S-1)+1>=S.

Thus stopping gives both useful entropies at least S,with all coordinates before the chosen child fixed. Continuing through every active interval would preserve entropy at least K_hp=S in a leaf of length below S,impossible. Prefix cost here is its exact conditional affine entropy/rank,not the number of fixed input bits. The source is allowed to have arbitrary correlations and an affine offset. Prefix fixing inside the selected useful child subsequently preserves the outer useful conditional entropy exactly.

Ordinary descent can also use sharper factors. Call a child high if its conditional entropy is at least H(parent)/(2t). If b children are high,they give the hidden block-wise affine source at that threshold. Otherwise fewer than b are high;the low children contribute less than half of the parent entropy,so some high child has conditional entropy greater than H(parent)/(2b). The b=1 losing case is impossible. Fix its predecessors and retain that child. The parent projection after prefix fixing has at least the retained child entropy. Therefore the exact local forcing guard still applies while that entropy exceeds its LOCAL K_j. Charge every response constraint to the same Dtot. Division by2b only scales previously spent ranks down. This proves the ordinary invariant with

 F_j=1 for j<=s0,F_(j+1)=2b_jF_j for j>=s0,
 H(current node)>=S/F_j-Dspent,
 hidden-block entropy at stopping >=(S//F_j-Dtot)//(2t_j).

Using F_j starting at s0,even when descent starts later,is conservative: every actual ordinary move has depth at least s0,so no entropy loss for a nonexistent earlier move is needed. All sibling forcing remains interleaved before the next ordinary move;earlier deterministic FAIL events persist. The hidden tuple survives every additional rank-D_res conditioning by the exact geometry inequality. The finite primitive/error,resiliency and rare-event arguments above are unchanged.

# Closed exact minimum

Write rank=Dtot,res=D_res and let Q be the maximum of

 rank+max_(j<H)K_j+1;
 F_(j+1)*(rank+K_j+1),for every s0<=j<H;
 F_j*(rank+2t_j*(K_j+res+1)),for every s0<=j<H;
 F_H*(rank+ceil(N/2^V_H)+1).

Each geometry predicate is equivalent to S being at least its listed integer. For example floor(S/F)>=c iff S>=Fc. The hidden nested floor is equivalent to S>=F_j*(rank+2t_j*(K_j+res+1));the strict leaf inequality contributes+1. Thus Q is exactly the least sufficient S,and Q>leaf_cap. Since S=ceil(d/C),S>=Q iff d>C(Q-1),or equivalently

 d_min=1+C*(Q-1)

is at most d. This is the EXACT least integer satisfying these geometry predicates. It is not a claim of an optimal true dispersion threshold. The controller computes this integer,rejects it if>N,copies its cached numerical plan with this dimension,and retains the separate geometry checker and exact primitive-error/output-event tests. The minimum accepted dimension is d_N;no accepted tuple gives d_N=N+1 and output0. Consequently the finite source proof and supplied full-B2 additive8 theorem, with accepted d_N>=5 as proved above, give4N-8d_N+8 for EVERY N,including the nonpositive fallback bound.

# First possible ordinary descent and finite dominance

The useful child g' of a consecutive active pair always has depth at least s0=a_1. Forcing earlier siblings on its ancestors does not make an ordinary entropy-dividing move: its retained child contains g' and hence has entropy at least S-Dspent. The unchanged parent guard S-Dtot>maxK therefore pays every earlier forcing operation. All such response ranks remain in the SAME global Dtot; none are omitted. Likewise conditional outer entropy and full-source resilience retain their original ranks.

Only the subsequent ordinary descent can divide source entropy by2b_j. All its depths j satisfy j>=s0. Define F_j=1 up to s0 and then multiply actual2b_j. If the actual useful child starts deeper,the larger product that starts at s0 remains conservative. The parent,hidden and leaf inequalities proved above therefore need only the stated j>=s0 terms. Forcing earlier siblings at these later levels remains interleaved before the next ordinary move. Prefix conditioning preserves the useful outer conditional entropy exactly;all response ranks and the separate rare-event/output certificates remain charged unchanged.

# Exact public response-rank charge before the internal winning node

Forcing visits disjoint ancestor/descent depths on one path,at most t_j-1 siblings per depth,each costing w_j*P. Prefix entropy costs remain separate. The total is

 Dtot=P*sum_(j=0)^(H-2)(t_j-1)*w_j.

Here the final internal depth H-1 is correctly excluded. The useful child g' cannot be a leaf,since its entropy is at least S>leaf_cap. Thus initial ancestor forcing has depth at most H-2. Carry the ordinary invariant with this budget until a winning hidden-tuple node is found or depth H-1 is reached. BEFORE any response fixing there,the ordinary high/low dichotomy either gives the hidden tuple or a losing child with conditional entropy greater than parent_entropy/(2b_(H-1)). The latter exceeds S/F_H-Dspent>=S/F_H-Dtot>leaf_cap by the strict final geometry predicate. But that child is a leaf,impossible. Hence depth H-1 is necessarily winning:its appended hidden block row supplies the challenge without forcing its children. Thus no operation excluded by the budget is performed.

Use D_res=P*(sum_(j<H)w_j+max_(j<H)w_j) as defined above. The SAME Dtot bounds both useful entropy losses;H=1 has an empty sum and still cannot accept an impossible leaf source. Here Dtot<=ell*P*sum_(j<H-1)(t_j-1)<=H*tmax*ell*P and D_res<=(H+1)*ell*P;the exact Q is monotone in both budgets on SAME plans. This compares certified expressions,not language complexities. No terminal arity is charged. Its logarithm remains bounded by the previous estimate,so fixed-tuple asymptotics are unchanged. Fixtures certify conditional geometry,not primitive availability.

# Public node-specific response widths and complete rank accounting

At depth j every actual challenge has exactly rho_j rows:leaves have zero rows,and the algorithm inherits all rows of one child then appends exactly binom(t_j,b_j) block rows. This proves the recurrence independently of source values or selected children. Thus any challenge tested at parent j has at most rho_(j+1) rows,including a more distant next-active descendant. The public response capacity is w_j=max(1,rho_(j+1)). The evaluator pads the challenge to w_j rows,constructs a FULL m_j=w_j*P-bit strong linear response map,and compares row prefixes. Both confidence tests and all child tests at that parent use this SAME full map and seed space. It does not pad a smaller extractor output with constants and call that full output uniform.

Default block constructors are unchanged. The factory first builds all arities,thresholds,blocks and products,selects H,checks total rows<=ell,then computes w and response parameters only for j<H. This avoids unused tail response checks. All arrays and targets are public functions of N. Actual take=m_j must pass every finite SU/source/field/seed guard;the general finite theorem in extractors.md therefore applies even if the previous global-capacity response could not pass. Since w_j<=ell,all targets are<=ell*P<=N. Complete polynomial bounds for accepted and rejected plans remain valid.

The auxiliary exact linear challenge-response lemma is a LOCAL statement for any fixed row count. Instantiate it separately at each level with w_j. When a fixed sibling challenge is forced,the full response has m_j uniform affine coordinates at some seed,so its preimage has that rank. Disjoint ancestor/descent depths and the unchanged last-level losing-child contradiction give Dtot=sum_(j<H-1)(t_j-1)m_j. Coordinate-prefix entropy losses are still paid separately. Both useful inner and outer conditional entropy lose at most this SAME Dtot.

Here is the variable-width resilience induction,without importing a global-width output hypothesis. At the internal winning node e,the block row is EB-close to uniform after EVERY additional full-source affine restriction of rank<=D_res,by the hidden entropy guard. Appending arbitrary correlated rows preserves a somewhere-uniform-row coupling. Suppose the true-child challenge at u+1 has error eta_child and remaining conditioning rank r. On any extra fiber of rank<=r-m_u,conditioning on one complete parent response at one seed adds rank<=m_u. The child row is still eta_child-close to uniform there. Matching all challenge rows implies matching its uniform row,with probability<=2^-P+eta_child. Union over the at most R response seeds bounds erroneous child selection by R*(2^-P+eta_child). Inheriting its challenge adds eta_child,so eta_parent<=(2R+1)*eta_child+R*2^-P. Parent block rows and padding are arbitrary appended rows and do not spoil the coupling. Induct upward,spending at most one m_u per level. Starting with D_res=sum_(j<H)m_j+max m_j leaves at least max m_j even after every possible propagation step. Therefore one additional rare-event response conditioning,of rank<=m_g,also fits. The uniform error is bounded by eta=(2R+2)^H*(EB+2^-P),and path deviation by H*R*(2^-P+eta),as before. All restrictions are on the FULL source;no independence of original blocks or uncharged outer projection is assumed.

At a useful active parent g,fixing its useful child coordinates leaves enough outer entropy for the FULL m_g linear response. For each fixed child value some seed is uniform,hence rare-FAIL lower bounds are2^(-w_g*p1),2^(-w_g*p2). They are at least the conservative2^(-ell*p1),2^(-ell*p2) used by unchanged output_events. Upper bounds R*(2^-p+eta) require only one resilient uniform challenge row after the full m_g response conditioning just charged. Both confidence maps are identical,so FAIL monotonicity applies. Earlier active and path errors are bounded exactly as in the primary finite proof. A useful child cannot be a leaf;empty leaf challenges are therefore not used as uniform rows in this argument. Every needed finite inequality remains explicitly checked.

# Uniform runtime

All six parameters and cached plans are public functions of actual N. The new coefficient C<=1+H*N has polynomial bit size. Primitive guards bound the block counts;the products F_H and every rejected d_min>N have polynomially many bits. Exact dyadic error arithmetic retains ell*P<=N and its supplied bounds. The grid has O(B^12) tuples. There is one fixed polynomial time exponent in actual N for all accepted and rejected paths. No source-dependent advice,witness,oracle,epsilon-selected program,ignored target bits or target padding enters the definition. Auxiliary BASE numerical claims are not used.


# Uniform discrete feasibility details for the current six-set

For a fixed tuple,strict clamps give beta_j=aZ/W_j eventually. Put c_true(beta)=(1-beta)/(1+4alpha*beta),alpha=1/(1-beta)+.00001,and f(w)=w*c_true(a/w). This f is increasing. The least prime SU degree and supplied prime asymptotics give its DEFAULT full-output seed d_seed<=(1+4alpha*beta+o(1))*L_j uniformly:the field inequality has logarithm L_j/d_SU+4beta W_j+o(W_j),d_SU=(alpha+o(1))*L_j/W_j,and d_SU tends to infinity even at Lpad. Hence actual b_j=1+ceil(d_seed/log2 g_j)<=(1/c_true+o(1))*L_j/W_j. These estimates hold uniformly over levels of this FIXED tuple.

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

 f(w)=100000*U*z²/[V*(100000*q_a²*U²+300000*p_a*q_a*U*V+4*p_a*V*z)].

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


