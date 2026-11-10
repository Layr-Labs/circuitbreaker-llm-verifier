# Selecting a uniform family by its proved source certificate

Trim_N=family.evaluate is one total uniform P language on ALL N actual
input bits. Its default SU output is trimmed by(L_j+1)^12, with exact
alpha=1/(1-beta); full response targets invert this capacity. Public
plateau/tail row budgets and variable active gaps are retained. The
controller selects least coarse stage certificate u_N, first ties. Its
analysis certificate v_N is defined in the final section below. At EVERY
N>=0, C_B2(Trim_N)>=4N-8v_N+7>=4N-8u_N+7; on fallback
output0,v_N=u_N=N+1. The evaluator and its selection are unchanged. The generic
full-B2 proof is supplied in elimination.md. No additive8 casework,
exact DP threshold or prior function identity is asserted.

The public ranges are h[.1,.2],a[.4,.7],wm[1.3,2],w0[2.4,3.4],
X[.25,.6],p[.2,.35]; beta remains[.15,.53]. X is root spacing;
D=w0+X and later active gaps grow as the entropy schedule declines.
All maps/finite guards use the NEW rho^-1 SU field bound. Its complete
cost-aware proof is in extractors.md. The logarithmic output trim pays
its advice without an artificial positive exponent gap. Seed enumeration
uses the finite margin1/10; actual full-seed coefficient is1/(1-beta).

# Exact energy and weighted-confidence witness

Strim is defined by the two finite certificates(PC) below. Use

 a=6249/10000,D=111863/40000,w0=25419/10000,
 wm=D/2+10^-6,X=D-w0=10187/40000,
 p=a/2-(102/10^8)*w0-10^-13,
 h=119340301699/1000000000000.                              (FM1)

All range/prefilter inequalities hold exactly, including2wm>D,
a/wm<.53,a/w0>.15 and p<a/2-.00000102w0.
w0+X=111863/40000=2.796575<D_trim=139829/50000=2.79658.
Put beta=a/w,alpha=1/(1-beta) and
f(w)=w*(1-beta)^2=(w-a)^2/w.

For w>a, f'(w)=1-a^2/w^2 is positive and<=1, and
f''(w)=2a^2/w^3>0. Thus f is increasing, convex and1-Lipschitz;
we retain the SAFE looser11-Lipschitz conversion constants below.
Convexity makes floored midpoint sums LOWER integral bounds.

The audit executes M=262144,Delta=(w0-wm)/M,scale10^15. Floored
f(midpoints) give a LOWER integral; floor(scale*f(right))+1 gives an
UPPER integral. All remainders and four original Fraction terms for
EACH sum are checked. Exact totals:

 MID=242318959769412276707,RIGHT=242319468752867876256.

The energy ratio R=(864h/109)*Delta*MID/scale is
1785857496442746346634861447341648444053
/1785856000000000000000000000000000000000>1.
With W=w0-wm,mu=R-1-8hW/10^12>0,
J=1+ceil(176hW^2/mu)=32782626. Since864/109<8, the SAFE larger
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

The exact integer sum is549999513468646731662; the three branches have
2,19260,242882 cases. The audit also checks the running suffix ends at zero.
It verifies the final confidence upper bound

 2h*[17/(18X)-Jlow]
 =3754643731142310263976092824071804383363479
 /12017074176000000000000000000000000000000000<p.

Both certificates use actual exact finite arithmetic. The pinned runner mounts the effective snapshot under /snapshot/files/.
The isolated entrypoint cannot import its sibling by the bare name.
This audit loads the actual immutable family.py explicitly with importlib,
checks its full SHA256, and calls its constructors directly. Earlier
packaging explanations that only the entrypoint was mounted were
inaccurate; the earlier passing embedded-code computations remain valid.

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

The common-S pair lemma below applies with S=ceil(d/C)>c. For any
actual active start, its ordinary factors are no larger than those
starting at s0. The displayed original geometry conditions therefore
dominate every corresponding(ST) target, since actual Theta,Hsrc<=K.
The same explicit ST ancestor/descent proof then supplies an internal
winning block source of entropy>K_j+D_res after all charged responses.
The strict leaf guard prevents a losing move to a leaf. This proves the
original d_plan route independently of the sharper stage certificate.

The pointwise primitive error used by the controller is obtained from actual constructor parameters. Let epsilon=2^-z from small_parameters, Jnear from near_parameters, w from reservoir_parameters, g from block_parameters and b its block count. Set etaNear=Jnear*epsilon, etaCubic=2*2^(-floor(w/2000)), etaGroup=ceil(log2 g)*(2etaNear+etaCubic), EB=(b-1)*(2epsilon+etaGroup). The smaller extractor's proof gives error epsilon at entropy K/2. Iterating it uses entropy at least K/2+floor(K/4), which the raw-source row proof preserves by its explicit rank guard (2g+ceil(log2g))*w<=floor(K/4). Its hybrid therefore gives etaNear. The cubic bound follows from the finite numerical rank and output guards in the cubic lemma; the floored exponent only weakens it. The row mixture induction gives etaGroup by summing the finite per-round error. The full-original-block history induction gives EB by summing each transfer's 2epsilon+etaGroup. These uses require every displayed finite primitive hypothesis, not an eventual-error simplification. The finite proofs use exactly their displayed numerical guards for these inequalities. Their later eventual simplifications of errors are used only in the separate analytic feasibility argument, never to certify a selected finite tuple.

For EBmax=max_j EB_j and R=max_j 2^seed_j, the full variable-width
resilience induction below proves eta=(2R+2)^H*(EBmax+2^-P) and path
deviation<=H*R*(2^-P+eta). The final challenge-response lemma gives
rare-FAIL lower probabilities2^(-w_(a_u)*p1_u),2^(-w_(a_u)*p2_u),
upper probabilities R*(2^-p+eta), and FAIL monotonicity. For EVERY
possible active stopping index u, output_events checks that both output
events' lower bounds exceed earlier-active errors plus path deviation,
using these actual parent row widths. Thus both output bits occur. All
pointwise error/guard claims here are finite, not eventual estimates.

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
Define Strim by two FINITE exact rational certificates:
(1) an original-scale10^12 left-box sum with(864h/109)*sum>1;
(2) a lower rectangle Jlow<=integral_wm^w0 Phi(r(w))/G(w)^2 dw,
constructed by the exact right-suffix/branch/isqrt procedure in the
witness section, satisfying

 2h*[17/(18X)-Jlow]<p.                              (PC)

The energy gives r(wm)=0,r(w0)=109/216 and endpoint Hfun values0,17/18.
All premises are actual finite meshes/integers. Phi(r)/G^2 increases.
Every argument below concerns EACH FIXED tuple with fixed margins.

# Fixed-tuple ordinary tree and retained geometry

Let c_true=(1-beta)^2,alpha=1/(1-beta),
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
The final2.79658 envelope below bounds WEIGHTED u, not this d_plan.

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

# One fixed uniform family and NP verifier

Fix W(N)=0 and V(x,w)=1 iff w is empty and family.evaluate(x)=1.
Every nonempty witness is rejected. Thus Trim_|x|(x)=1 iff an accepted
witness of length at most W(|x|) exists, in both directions.
The grid enumerates BASE then every six-tuple with denominators<=B,
O(B^12) tuples in fixed order. Every FIXED rational tuple eventually
occurs in every later grid. Finite guards protect all selected tuples;
no uniform threshold over shrinking grid margins is assumed.
For each fixed feasible tuple and D'>w0+X, the STAGE certificate
proved below is<=2^ceil(D'Z) eventually. Therefore minimum-u selection
and FM1 give u_N<=2^ceil(D_trim Z) for every sufficiently large N.
Since Z=o(log N), this yields4N-N^epsilon for each fixed epsilon>0,
and liminf C_B2(Trim_N)/N>=4 for the same new fixed family.

# Original geometry predicate
The stored d_plan is independently checked by all displayed finite
geometry guards with S=ceil(d_plan/C). Only their satisfaction is needed,
not an assertion that d_plan minimizes true dispersion. In particular
d_plan>=S>=max K_j+1>=3. The source helper and strict first-tie
controller are unchanged.

# Stage targets and a weighted source certificate

Keep every actual component map. For its returned SU tuple define
T(c)=400mhds^2+ceil_log2(32m)+3. The finite proof gives strong error1/4
at this entropy, not only at ordinary tuning K. Let Theta_j=T(response
at its actual inverse-output Kresp) and Hsrc_j=(BLOCK) from extractors.md.
Both are<=K_j by explicit constructor guards. Changing the source entropy
in this proof does not change any field, seed, output bit or error.
Put rank=Dtot,res=D_res,O=rank+max_j Theta_j+1,c=ceil(N/2^V_H).

For active edge u with child depth s=a_(u+1), define F^u_j=1 for j<=s
and F^u_(j+1)=2b_jF^u_j for j>=s. Its INNER target I_u is the maximum of

 O;
 F^u_(j+1)*(rank+Theta_j+1), for s<=j<H;
 F^u_j*(rank+2t_j*(Hsrc_j+res+1)), for s<=j<H;
 F^u_H*(rank+c+1).                                  (ST)

All targets are public integers, with I_u>=O. The exact factors in(ST)
are computed by source_certificate's inner/factor loop.

Here is the sufficient active-pair argument, without an exact-threshold
claim. For an arbitrary affine source, let p_k be its prefix-rank
increments, each0 or1. After fixing the prefix before a child, its inner
conditional rank is the sum of p_k in that child. The parent's outer
conditional rank given this child is the sum in the following children,
by the affine fiber rank identity. A final active leaf has at most c
pivots. Suppose every child subtree avoids all later useful pairs and
has at most B_next pivots. If no child reaches I_u, total rank is at
most T_u*(I_u-1). Otherwise take the FIRST child reaching I_u. Its
predecessors cost at most(T_u-1)*(I_u-1), it costs at most B_next, and
its suffix must have rank<=O-1 to avoid the current pair. Thus a uniform
upper bound for a bad parent is

 B_parent<=B_next+T_u*(I_u-1)+(O-1).

No relation between B_next and I_u is assumed. Backward induction proves
that every affine source of dimension at least

 weighted=c+1+sum_u[T_u*(I_u-1)+(O-1)]               (WC)

has consecutive useful active nodes with inner rank>=I_u and outer>=O,
after the preceding coordinate prefix is fixed. This is an upper forcing
certificate, not an exact maximum, attainment or monochromaticity claim.

To use such a pair, force earlier ancestor siblings to FAIL with their
FULL m_j-bit responses. Outer O pays the parents through a_u; inner I_u
pays intermediate parents. Every forced response costs at most m_j,
with total<=rank. Prefixes inside the chosen child preserve the outer
conditional rank, and linear response fixings lose at most their rank.
Every required responding map has source entropy>Theta_j by(ST).
All previous forced events persist pointwise under later restrictions.

For completeness, ordinary descent obeys the following finite calculation.
If parent rank is k, call a child high at conditional rank>=k/(2t_j).
With b_j high children, omitting other previous children only increases
their conditional entropies, giving the desired block source. Otherwise
low children contribute<k/2 in total, so one of the fewer than b_j high
children has rank>k/(2b_j). Fix its predecessors and force them to FAIL
BEFORE descending. The parent retains at least this child's entropy;
(ST)'s next-factor term supplies its Theta_j guard. The b_j=1 losing
case is impossible. Dividing by2b_j scales old rank losses down, so at
ordinary depth j starting from s the retained rank is at least
I_u/F^u_j-Dspent. Its hidden blocks have conditional ranks at least
floor((floor(I_u/F^u_j)-rank)/(2t_j))>=Hsrc_j+res+1.
The last term of(ST) rules out reaching a leaf. In particular no losing
move at H-1 is charged. Thus the SAME rank and full-source resilience
budgets Dtot,D_res suffice. At the winning node its row is supplied by
(BLOCK), including after all rank-res allowed conditionings.

The pointwise error, resilient challenge, correct-path and output-event
proofs in the primary finite theorem apply verbatim, using these actual
Theta,Hsrc targets instead of the larger tuning entropies. Every stopping
active index is covered by the actual checked output_events inequalities.
Both output bits therefore occur on the original affine source.

The original geometry uses a sharper SPECIAL common-target bound.
When I=O=S at every stage, bad-parent rank is at most
B_next+(T_u-1)(S-1). If B_next<S, no child qualifies and
T_u*B_next<=B_next+(T_u-1)(S-1). If B_next>=S, an all-low total
T_u(S-1) satisfies the bound. Otherwise the first high child costs
B_next; if it is last its predecessors cost at most(T_u-1)(S-1),
while if not last its predecessors and suffix together cost at most
(T_u-2)(S-1)+(S-1). Thus the same bound holds. Starting at c<S,
root bad rank is at most C(S-1), C=1+sum(T_u-1). For
S=ceil(d/C), C(S-1)<d. The original strict leaf guard implies S>c.
This proves the common-S pair use in the primary geometry proof,
without importing an exact DP or attainment assertion.

Original geometry independently proves dispersion at d_plan. The program
uses the valid finite certificate

 u(p)=min(d_plan,weighted),                          (SC)

checks EVERY candidate's original geometry/primitive/error/output guards,
then selects least u(p), first ties. At any accepted length u(p)<=d_plan<=N
and u(p)>=1. On fallback u_N=N+1 and output0. The complete generic
full-B2 theorem in elimination.md gives

 C_B2(Trim_N)>=4N-8u_N+4 for EVERY integer N>=0.      (GT)

Fallback satisfies(GT) directly, so no false disperser premise is used.

# All-length envelope for the stage certificate

For each fixed feasible tuple, the response inverse is
Kresp=max(2,ceil_root((m_j*(L_j+1)^12)^q,p)) for beta_j=p/q.
Theta_j<=Kresp<=K_j exactly. Since log m_j=(p_conf+o(1))Z,
log Kresp=(p_conf/beta_j+o(1))Z<=(p_conf*w0/a+o(1))Z.
This is strictly below w0/2<wm. It accounts for the actual SU h factor;
it never calls h poly(log N). Rank/res have exponent p_conf<wm and
Hsrc_j<=K_j=2^W_j. The extra12log(L_j+1) is o(Z).

For ANY s, the entropy schedule/product proof(VF) gives, uniformly in j,
log2 F^u_j<=W_s-W_j+(tau/4+o(1))Z for every fixed tau>0. Each hidden
term in(ST) consequently has log<=W_s+(tau/4+o(1))Z; t_j costs o(Z).
The Theta term uses W_(j+1) in its next factor, which is still>=wmZ,
strictly above the global Theta exponent. The leaf term is smaller by
the proved log(c+1)<W_H-fixed_margin*Z. O also has exponent<wm.
Taking tau arbitrarily small shows log I_u<=W_s+o(Z), uniformly in u.

The active gap has log T_u<=DZ-W_(a_u)+o(Z) and W_child<=W_parent.
Therefore each T_u*I_u<=2^((D+o(1))Z). There are only O(L) stages;
O and leaf terms have smaller exponents. Hence weighted<=2^((D+o(1))Z),
even though d_plan may have the larger original geometry envelope.
For any fixed D'>D=w0+X, that fixed tuple eventually passes all guards
and gives u_N<=2^ceil(D'Z). Every fixed rational tuple eventually occurs
in each later denominator grid; no uniformity over shrinking selected
margins is assumed. The exact witness supplies the stated explicit D'.
Since Z=o(log N), (GT) yields4N-N^epsilon for every fixed epsilon>0
at every sufficiently large integer N, and liminf C_B2(Trim_N)/N>=4.

The stage factors have at most O(L^2) bits (H<=L,b_j<=2N). The helper
uses O(H^2) polynomial-size integer operations and already-checked
constructor tuples. Taking the minimum and the O(B^12) grid retain one
fixed polynomial runtime. No source ranks, pivots or unknown witnesses
are queried by the evaluator; those only establish its certificate.

# Exact response-rank charge
At depth j, at most t_j−1 earlier siblings cost full responses m_j=w_jP.
Thus rank=Dtot=sum_(j<H−1)(t_j−1)m_j. ST forbids a losing move at H−1:
it would retain rank>c in a leaf of length<=c. Its response map still
exists and is charged to res=D_res=sum_(j<H)m_j+max m_j. Prefix ranks
are separately charged; every earlier FAIL persists. Since w_j<=ell,
rank<=H*tmax*ell*P and res<=(H+1)ell*P.

# Public variable response widths and complete rank accounting

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

 f(w)=z²/(q_a²*U*V).

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
length n_j, set L_j=ceil_log2(n_j+1), and call response_parameters(n_j,K_j,ell*Pmax,mode_j).
This wrapper retunes Kresp=max(2,ceil_root((ell*Pmax*(L_j+1)^12)^q,p)) for beta=p/q,
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



# Sharper finite certificate for the UNCHANGED selected function

Keep family.py and audit.py byte-identical to the Trim construction above.
In particular source_certificate still returns the COARSE(SC) certificate
u(p), and tree_parameters still selects its least value, with first ties.
No tighter value below is substituted into that controller. Thus this
section strengthens the analysis of exactly the SAME Boolean function.

For the actual selected plan, keep the same targets I_u>=O from(ST) and
leaf capacity c. Define

 sharp(p)=c+1+sum_u[(T_u-1)*(I_u-1)],
 v_N=min(d_plan,sharp(p)).                              (SH)

On fallback define v_N=u_N=N+1. Each nonfallback target and dimension
is a positive integer, so v_N>=1. This is a public integer formula,
computable from the returned plan without any source or unknown witness.

Here is the improved UPPER counting lemma, not an attainment claim.
Suppose every child subtree has bad-source capacity at most B, meaning
its pivot rank is at most B if it avoids all later useful pairs. At the
current edge put I=I_u>=O and T=T_u. If B<I, every child has fewer than
I pivots, and their total is at most TB<=B+(T-1)*(I-1). If B>=I and
all children are low, their total is at most T*(I-1), also bounded by
B+(T-1)*(I-1). In the remaining case choose the FIRST high child.
If it is last, its predecessors cost at most(T-1)*(I-1) and it costs
at most B. If it is not last, its predecessors cost at most(T-2)*(I-1),
it costs at most B, and the suffix must cost at most O-1 to avoid the
current useful pair. Since O<=I, the same bound follows. Consequently

 B_parent<=B+(T-1)*(I-1) whenever I>=O.                 (SB)

Induct backwards starting at final active leaf capacity c. Every affine
source of dimension at least sharp(p) therefore has a useful consecutive
active pair at this plan's unchanged I_u,O. The preceding(ST) proof
forces all ancestors and ordinary descent from exactly those targets;
no later part of the challenge/error argument changes. Hence the actual
selected function disperses all affine sources of dimension at least
v_N, using either the original d_plan route or the new sharp route.
For a nonfallback plan c>=1, there is at least one active edge with
T_u>=2, and I_u>=O>=2: every actual SU entropy charge is positive and
O=rank+max Theta+1. Hence sharp(p)>=3. The original geometry has
d_plan>=S>=max K_j+1>=3. Thus v_N>=3, so the complete
d>=3 refinement in elimination.md Section6 applies. It uses the
complete unrestricted Ramsey R refutation and seven-variable read-once
Q cover, both fully validated by the declared audit. These scoped
computations are required premises of the supplied universal terminal
reduction, not samples or prior verdicts. This gives for every integer N>=0

 C_B2(Trim_N)>=4N-8v_N+7>=4N-8u_N+7.                    (SGT)

Fallback satisfies this negative bound directly, including N=0.
For each accepted plan,
coarse_weight-sharp(p)=sum_u[(I_u-1)+(O-1)]>=0. Thus v_N<=u_N for the
SAME selected plan, and(SGT) is pointwise at least the preceding(GT).
The inequality can be equal if the d_plan clamp already dominates; no
unconditional strict improvement to actual minimum circuit size is
asserted. Strict improvement occurs to the two unclamped source upper
formulas whenever the displayed positive sum is nonzero.

Finally v_N<=u_N transfers the proved EVERY-length D_trim envelope to
v_N with no new grid-selection or uniform-margin assumption. Thus the
same fixed Trim_N family satisfies the joint full declaration, including
4N-N^epsilon at every sufficiently large integer N and liminf ratio4.
The claimed asymptotic score remains[4,1]; this is a finite analytical
refinement of one fixed evaluator, not a change in the asymptotic score.

# Complete finite terminal predicates for additive7
The complete elimination.md supplies all legal affine-pair moves, inert
ledger and terminal reductions. The obsolete nonnegative baseline
potential proof is omitted; its independent J=s−2q ledger is retained
and proves every bound used here. No mathematical step relies on that
omitted potential or on a prior accepted claim.
The audit validates the ENTIRE RamseyR and read-onceQ predicates in
Sections5.2/6.1, including regenerated domains, flats, normalization
groups, CNF roots, ordinary RUP propagation and every cover assignment.
The source certificates' decoded bytes are unchanged. Their joint
length-prefixed archive is XZ/base85 encoded in the audit, bounded by
200000 decompressed bytes; both original bounded decoders then run on
exactly their original envelopes. All original clause, index, witness,
end-of-data and coverage checks remain. The numerical mesh and actual
module constructor probes also execute. No accepted tree, full field
table or raw extractor output is executed. Successful computations
establish these finite predicates; the supplied mathematical reductions
establish the universal gate/family statements.
