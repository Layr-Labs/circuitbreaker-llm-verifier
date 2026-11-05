# One fixed uniform family with log(N)^9 loglog(N)^5 / logloglog(N) loss

For the single total language evaluated by family.py, we prove that absolute B>0,N0 exist such that for EVERY actual integer input length N>=N0,

 C_B2(F_N)>=4N-B*(log2 N)^9*(log2(log2 N))^5/log2(log2(log2 N)).

Both canonical claims assert this ENTIRE declaration and all stated eventual weakenings for the SAME family. Their protected coefficient/exponent guarantee and projection remain(4,1). No larger endpoint is claimed. All required extraction arguments are supplied in seeded.md and disperser.md; the complete full-B2 gate theorem is in elimination.md. Citations identify provenance, not unexplained premises or instructions to retrieve external papers.

## Total program and actual input

The entire evaluator, scalar plans, arithmetic and compilation are in family.py. N<128 outputs0, including N=0. For larger N put L=ceil(log2 N),z=ceil(log2 L),D_int=max(2,ceil(10z/7)),C_W=31500,C_Y=86000. The actual prefix is the standard Elias gamma code of positive S: j=floor(log2 S) zeros followed by the full MSB-first binary representation of S. Its length is H(S)=2j+1. Truncation,more than L leading zeros,or S>2^L rejects. Thus1<=S<=2^L<=2N on every admitted input. Every header is handled by the same program. Let P2 round upward to a power of two. Exact scales are

 (a,B_cap,E0) is selected by the exact feasible optimization in C29,s_idx=log2 a,
 (U,V,T)=unpack_scales(S),
 h=P2(U L^2),e=P2(V L^2 z),ell=P2(8ze),
 eprime=100h*a*ell+e+2,q0=(2h*(a-1)+C_Wez^2)ell,
 d=max(q0,(2h*a+10)ell,2P2(8*(ell+eprime+L+10))),
 Y=100B_cap*h*d+C_Yez^2*P2(d),
 k=hY+100B_cap*h*a*ell+10d+C_Wez^2ell.

The index alphabet is the selected power of two2<=a<=h,with s_idx=log2a and b=ceil(d0/s_idx) and q0=(2h*(a-1)+C_W*e*z^2)ell. h is ROW independence,e is error precision; they are not equated. The bounded exact k is computed and ceil_log(k) is checked before source binding and k>N rejects. All source/seed/prefix/compiler/layout guards are exact. Failed plans return0. The initial d0<=5L cap precedes initial-seed enumeration. ell>=u ensures a whole biased group fits one completed row.

Initial A1 and all B5 field/matrix data are canonical by C1, with no input-data bits. Short E1 is data-free. A single paid internal package follows the header. First compute the four ordinary scalar plans,put s_star=max_i s_i,and replace each field degree by s_star while retaining its own nprime,H,m0,error,calls. Recheck exact advice/entropy and ell-seed guards;reject failures. C6 proves eventual success and universal sharing. The one package has t_star-a(t_star) bits,t_star=D_int*s_star,and C5/C3 decode/compile it once. Each map is separately reparameterized with its own source/output metadata and that shared field/cyclic action;bad data yield bounded ring or zero maps. No further data are read. The engine checks its cursor,removes the ENTIRE P and pads payload with P zeros to N coordinates.

Good primitive binary packages exist at every required degree by finite-field cyclicity, and C5 represents every such polynomial at large degree (the small fallback represents a reciprocal choice). They give exactly the universal A1 field/cyclic-matrix hypotheses. The one good internal tuple is only a PAID hardness restriction,universal across all four maps by C6. The evaluator neither searches for nor tests their primitivity. Arbitrary data still define bounded ring/compiler calculations or zero maps. There is no advice, witness data or external table defining the language.

## Supplied chain and exact resource closure

A1 proves strong binary-linear extraction over a field with a specified basis and a cyclic matrix. B1 supplies adaptive rank condensation;B2.1 supplies the fresh-seed composition and all dimension/error dependence. Each DyadicPrefixPlan requests FULLP2(M) and returns only firstM bits,with every stipulated seed retained. Output projection preserves strong error and fixed-seed linearity;discarded suffixes enter no later seed,message,branch or conditioning transcript.

W/T keepM=ell and useK=C_Wez^2ell. E2/E3 keepM=q0,d and useK=C_Yez^2P2(M). The smaller W/T threshold has orderL^4z^5,not the worst row thresholdL^6z^8. Current C29 provesH_WT<=1001z/100,H_E2E3<=136z/5,R<=229z/10,s_star<4267e/1000. Thus W/T advice<=(978120143/62500)Reqez^2<K/2=15750Reqez^2;row advice<=(132891448/3125)Reqez^2<K/2=43000Reqez^2. Both strict gaps pay the ACTUAL4rho_log+10,where rho_log=e+ceil_logK+13+ceil_logReq for each ownK. All first FULL candidates pass ONEcall eventually. Common-degree revalidation retains each map's ownnprime,H,K,output,error,calls;one good paid tuple works for all four. Actualseed<79ze/10 is still paid byell>=8ze. ShortE1 keeps direct B2.3 with source10d,error ep and exact actual seed.

The complete B3--B4 refreshed look-ahead proof in disperser.md usesC_W for W/T/X/Q thresholds andC_Y for E2/E3/Y thresholds. q0 leaves its exact W/T residual. Y pays100bhd and both actual FULL E2/E3 thresholds throughC_Yez^2P2(d);k separately payshY,transcript losses,10d andC_Wez^2ell. The maximum d paysseed/transcript/Q guards,andep=100haell+e+2 pays the weak-seed rank failure. For ANY at mosth good initial indices,the wholeell-bit rows are jointly within100bh^2*2^-e of independent uniform rows. There is no independence premise on matrices or discarded suffixes. B7 supplies all exact guards and all-length closure.

Initial dimension10 uses the selected precision E0. The exact selector enforces E0>=ceil_sqrt(L)+ceil(L/125) and E0<=L,initial validity,d0<=5L,b<=B_cap,and125E0-d0>=125ceil_sqrt(L). The initial requested output is M_init=P2(Y);only its first Y bits are used. Because E0>=L/125,generic B1 has condensed nprime=O(k). Thus H=O(k^(1/10)),field degreeO(L),and first-candidate adviceO(Y*k^(1/10)*L) has strict L exponent89/10<9. One call eventually suffices uniformly across selected precisions. Good canonical initial data meet the actual scalar hypotheses. The deficient-row fraction is at most2*2^-E0.

A2 and B5 supply the logarithmic-external-seed bootstrap in full:1024-wise entropy sampling,clipped probability-tree weights,repeated-position collision moments,seedwise smoothing,conditional two-block repair,direct GENERAL-source A1 on the second block and shared-multiplier weak hashing on the first. Rare transcripts and every seed are paid; a constant-bounded number of fresh-seed repetitions recovers the output. Whole G is generally nonlinear. B6 uses only its proved strong general extraction below r^-80 and shifted agreement below r^-200,never a rank promise about G.

## Whole-row robust input and one hard header

For each selected initial index i,the evaluator computes the complete W_i row and takes its first u bits as ONE biased group. B4's joint-row guarantee passes to these prefixes and biased values by data processing. Inner Bonferroni events depend on O(r^2) GROUPS and hence O(r^2) GOOD ROWS. The number of compared terms remains exp(O(r^2u)); independence h therefore needs O(r^2),while error logarithm e still needs O(r^2u). B6 proves this separation,ideal probability estimates,exceptional partitions and finite absolute constants without an external fooling theorem.

Use S_groups=B*M<=M0=2^d0 rows. Exact layout maximizes r subject to m+b_seed<=d0,with m=floor(r/125),b_seed=floor(T*ceil(log2 r)/(30*125^2)),u=2b_seed+2ceil(log2 r)+10. Every bad row gives at most one bad group; monotonicity bounds arbitrary dependence by setting bad-group bits all0/all1. Coverage loses at most2^O_S(log r). Under the cap r<=1000L and at fixed S,u=O(Tz). One fixed sufficiently large S0 meets h>=A_t*r^2 and gamma<=2^(-A_gamma*r^2u),all bias/seed/source constants,and ell>=u.

At this S0,the initial bad-row count is at most2M0*2^-E0. Coverage gives q/S_groups<=2^(-E0+O_S(z)),while S_groups^(1/125)<=2^(d0/125). The EXACT selected margin E0-d0/125>=ceil_sqrt(L) makes their product at most2^(-sqrt(L)+O_S(z)),which tends to zero. Hence the required bad-group allowance holds eventually. The fixed gamma code fits at EVERY sufficiently large integer N and every guard holds there. No fixed constant error exponent or favorable rounding is assumed for the selected initial plan.

## Uniformity, NP and gate application

Canonical initial search is polynomial on every admitted header by s0<=5L/11 and fixed dimension10. B5 searches are polynomial in r by its uniform O(log r) degrees and fixed helper dimension4; C1 proves these bounds exactly. C5 enumerates only small irreducibles of degree O(log t),uses bounded CRT and has fixed polynomial complexity in t. The one actual package t_star-a(t_star)>=t_star/2 at large degree implies t_star<=2N; small fallback degree is bounded. C3 kernel/module reconstruction is binary linear algebra and ring arithmetic of polynomial size. It enumerates no2^t elements. Bad indices/rank failures are bounded then zero.

Growing D_int=O(loglog N) is not a runtime exponent. The evaluator streams ONLY its first nprime coefficient monomials,with polynomial tuple/exponentiation cost and no full grid. All admitted finite dimensions are polynomial in N. Tower arithmetic uses explicit stacks; B5 designs use exact polynomial-time CE. The d0 cap bounds2^d0 by(2N)^5. Threshold integers have at most uBM+O(r)<=uM0+O(r) bits; u<=ell<=N gives a fixed polynomial N bound. Thus all headers/data,cutoff/rejection paths and normal evaluation use ONE fixed polynomial-time algorithm.

Let L_yes={x:evaluate(x)=1}. The NP witness bound is0. Accept exactly empty w when the same evaluator returns1. Every yes input has that witness,and every accepted allowed witness implies a yes input. Both directions concern this SAME fixed language; witnesses are not circuit inputs.

At S0,each internal s=O(e)=O(L^2z),D_int=O(z),so current P=H(S)+D_int*s_star-a(D_int*s_star)=O_S0(L^2z^2)=o(k),including every header/data bit. Fix one good shared internal package and the header. Injective zero padding preserves payload affine dimension. B1--B7 give an affine disperser at dimension k on N-P actual variables. Apply elimination.md's complete4n-8d theorem,covering all16 B2 tables,sharing,constants,repeated pins,affine exchanges and terminal cleanup. Restriction adds no gates:

 C_B2(F_N)>=4(N-P)-8k=4N-(4P+8k).

Since k=O_S0(L^9z^5/log z),the exact loss is O(L^9z^5/log z). Rounding ratios L/log2 N and z/log2 log2 N tend to1 at all integer lengths,so absolute B,N0 prove the declaration. For this SAME family any fixed extra loglog exponent absorbs B,any fixed extra log exponent absorbs the loglog factor,and N^epsilon dominates the fixed polylogarithm. These are exactly the other stated inequalities and joint liminf projection.

## Provenance and finite evidence

This revision builds on C9 explicit envelopes6a291850-3ba9-440e-83b7-fe169de7d8d2,06f789be tight source allocation,46681146 gamma coding,77e43640 common field,71e6d94e unit coding and c1f8eee2 whole-row closure. C10 source84bbd016-28f8-4a06-802f-9932903344ee supplied growing index words/divisor;C11 source6b4b7dde-5bf3-4091-b601-16bd48e0af8d supplied exact integer accounting;C12 source716f7ddf-5173-4f9c-b4e0-32f636197d64 supplied initial prefix/larger alphabet;C13 source16909354-df8a-4331-94e0-31ae44bce9c9 root-scale allocation is prior work;C14 source3c0edbbe-db2c-403f-bbeb-a99ebd4060f5 supplied packed parameters;C15 source613ac482-02fc-449b-83e3-74fa80219735 supplied exact Y reserve;C16 exact weak-seed precision is prior work;C17 exact Q/maximum seed budget is prior work;C18 reduced entropy coefficient is prior work;C19 reduced seed envelope is prior work;C20 initial precision/round cap is prior work;C21 exact alphabet optimization is prior work;C22 joint initial allocation is prior work;C23 single-call dyadic internal prefixes are prior work;C24 tight shared-field coefficient is prior work;C25 full-request entropy payment is prior work;C26 reduced full-request coefficient is prior work;C27 seed-budget dimension is prior work;C28 split source thresholds are prior work; C29 tunes dimension for the look-ahead source. C1--C4 derive from local/submitted reciprocal89c0424d,companion a50a59d7 and earlier f72b2eba; whole-row grouping/9,5 closure from c1f8eee2; growing dimension from61c0ea3e and charged direct maps from90a83107. Those endpoint verdicts are no premises. Older source/code and primary mathematical provenance are identified in the complete supporting proofs and public dependency note. Failed/unpublished sources are neither accepted results nor proofs of negation. No prior_results premise is used.

The declared audit establishes ONLY its scoped finite algebra,compiler,codec,packing and scalar predicates. Large scalar plans allocate no corresponding input/field/seed. Large canonical searches and the eventual complete family are not executed. Unused half-output code/checks and duplicate narrative were removed; no removed predicate is claimed retained. Universal correctness rests on the supplied mathematics. Advisory critique,configured-provider judgment and human promotion remain separate.

The complete published inherited material provider is3c0edbbe-db2c-403f-bbeb-a99ebd4060f5 (all source bytes verified at8bd222454102e93425ed43d6bdf9221306aef984). A content audit,not a graph-path assumption,finds57 inherited runtime definitions unchanged with PolyPlan modified and explicit own scalar helpers/prefix classes for resource candidates,parent allocation,current joint selection and hidden-suffix outputs,identical imports/module statements,byte-identical gate proof,A1,A2 and B5-B6. The earlier C17 carry audit also matched C1-C10/B1-B2;current C18-C29 explicitly updates coefficient,seed accounting,initial dimension/precision and allocation selection in C6/B2.1-B2.2,so those larger block identities are no longer asserted. The remaining scientific changes are the supplied own C15-C29 guards/resource arguments. Own C15 parent613ac482-02fc-449b-83e3-74fa80219735 is also a published direct provider;all eight source bytes were verified atf8c8605f1c647f6d0739fa4f388725e55ab0ec4b. Own C16 parentbf1f92a0-1e66-4e77-84a2-1752d9636b41 is presently unpublished and requires a later formal link after both endpoints publish. Older ancestry is credited below;formal direct-provider scope is not a claim of independence or a certified contribution graph.

Historical credit is retained in the public note and provenance paragraphs;all mathematical premises are supplied here.


# C15--C17. Exact source,precision and seed payments

The exact Y-reserve method derives from own613ac482-02fc-449b-83e3-74fa80219735,weak-seed precision frombf1f92a0-1e66-4e77-84a2-1752d9636b41,and exact Q/maximum-seed accounting from8fadd573-138b-49e7-8891-7baf47e8adc9. Their judgments are no premises. Current coefficients distinguish the actual source requirements.

CurrentY=100B_cap*h*d+C_Yez^2P2(d). Sinceb<=B_cap,d>=q0,residualY-100bhd>=C_Yez^2P2(d) meets BOTH E2/E3 actual thresholds. The full B4 proof fixes at most100hd targetY bits per round,and pays the entirehY observation from X separately. No fractional Y reserve or1000bhd guard is used. InitialY supplies its firstd bits. Every seeded output is linear after its seed is fixed;the supplied B3/B4 proof includes actual old/current/next seeds and targetA-parts. Final-round refresh and100bh^2*eta error use those exact premises.

Precisionep=100haell+e+2 gives weak-seed failure2^(100haell)*2^(1-ep)<=2^(-e-1). q0-2h(a-1)ell=C_Wez^2ell pays every T call,including equality,as each B3 transition fixes at mosth*ell seeds andh*ell A-parts and onlya-1 transitions precede its last W output. X's W threshold is paid separately by k. d=max(q0,(2ha+10)ell,Dseed),whereDseed=2P2(8*(ell+ep+L+10)) is the EXACT ShortSeedPlan seed. It also pays10d>=ell+2ep+10. Qpads to d. Coincident/shared maps are permitted by the complete proof. Every finite actual source/advice/seed predicate is rechecked;current closure proves eventual validity rather than importing a judgment.

# Scalar bootstrap used by both allocations

Fix a packed code, so h=Theta(L^2), e=Theta(L^2z), ell=Theta(L^2z^2). This calculation applies to either fixed pair of positive coefficients once a feasible scalar witness gives selected k=O(L^9z^5/logz). Initial validity and the seed/round guard imply d0>=PVseed>=L and B>=L/j. Since j<=log2h=O(z), initially B=Omega(L/z).

The exact formulas give k>=100*B*h^2*d and d>=2*h*a*ell. Hence k>=200*B*h^3*a*ell and a=O(z^4/logz). Therefore j=log2a=O(logz), improving the lower bound to B=Omega(L/logz). Also q>=C_W*e*z^2*ell=Omega(L^4z^5), so k>=100*B*h^2*q gives B=O(L/logz). Thus B=Theta(L/logz). Then d<=k/(100*B*h^2)=O(L^4z^5), while d>=q gives the matching lower bound. Finally a<=d/(2*h*ell)=O(z^3).

The row reserve divided by100*B*h*d is at most(C_Y/50)*e*z^2/(B*h)=O(z^3logz/L), tending to zero. Thus Y=Theta(L^7z^5/logz). The remaining X terms have ratio at most a*ell/(h*d)+11/(100*B*h^2) to hY, because C_W*e*z^2*ell<=q<=d. This ratio tends to zero using d>=2*h*a*ell. Consequently k=Theta(hY)=Theta(L^9z^5/logz), including the matching lower bound. No internal seeded-map validity was used: these are positive scalar formulas and the selector's initial guards. Below the reference obtains its upper bound from a direct witness, and the current allocation obtains it from a feasible initial comparison, before either invokes internal field estimates.

# Supplied C28 parent construction and independent closure

The reference is the actual C28 source94b25c1d-24ec-4e06-8a24-8e78d206f2ce: C_Wp=36500,C_Yp=73000,D_p=max(2,ceil(29z/16)),ell=P2(8ze). It uses the same finite initial selector with these resource constants, full-request map thresholds, exact X/Q and Y residuals, header, layout and canonical rules. The current parameterized split helper can evaluate its resource class by explicitly passing36500,73000. No source judgment supplies its closure.

For a direct feasible witness set a0=P2(max(2,floor(z^3/ceil_logz))), j0=log2a0, B0=ceil(337L/(100j0)), E00=ceil(L/37)+10. Its q,d have orderL^4z^5, Y has orderL^7z^5/logz, and k has orderL^9z^5/logz. The row reserve is O(L^6z^8), smaller than100B0hd. Initial condensed length is O(k); full D10 advice has strict L exponent8.9<9. Its PV seed leading bound is78L/37; direct seed contributes44*176L/(166*37). Total10346L/3071+o(L) is below337L/100 by a strict linear gap327L/307100. Error margin is at least(31/462500+o(1))L, exceedingceil_sqrtL. Thus actual capacity,5L,round and exact error guards hold at every sufficiently large integer length.

The exact positive lower bounds give the parent's selected a=O(z^3), B=Theta(L/logz), d=Theta(L^4z^5), Y=Theta(L^7z^5/logz), k=Theta(L^9z^5/logz) by the supplied bootstrap. W/T have nprime=O(L^4z^5), so5^29>2^64 gives ceilroot5 and H<=91z/10. Rows have nprime=O(L^6z^8), so10^29>2^96 gives ceilroot10 and H<=91z/5. With R<=291z/10,s_star<427e/100, the parent W/T advice11307387/625Reqez^2 is below18250Reqez^2; row advice22614774/625Reqez^2 is below36500Reqez^2. Both gaps pay the actual4rho_log+10 with rho_log=e+ceil_logK+13+ceil_logReq, since the term is O(e+logs) and Req*z^2 tends to infinity. First full candidates and common-field revalidation pass ONEcall. Actual seed<79ze/10 is paid by ell8. The conditional B3/B4 proof uses its own two coefficients and exact residuals; robust error, one good paid header, canonical/runtime/NP and full gate arguments above therefore supply its complete resource guarantee independently. Complete P_p=O(L^2z^2).

# C29. Tune dimension for the look-ahead source

C29 builds on actual C28 source94b25c1d-24ec-4e06-8a24-8e78d206f2ce. Current constants are C_W=31500, C_Y=86000 and D=max(2,ceil(10z/7)); ell=P2(8ze) is retained. The split helper now has explicit coefficient arguments, with these current defaults. Historical C28 fixtures pass36500,73000. Initial and B5 dimensions retain their fixed values. All current X/Q/W/T thresholds use C_W, while Y/E2/E3 use C_Y exactly.

Let the parent's selected alphabet be a_p=2^j. Its upper k_p=O(L^9z^5/logz) and lower k_p>=cL^9z^5/j give j=Omega(logz), and its bootstrap gives j=O(logz). Choose a'=a_p/2, j'=j-1, B'=ceil(d0_p/j'), E0'=E0_p. Then B' is in the current enumeration range, B'/B_p<=1+o(1), and B'>=B_p/6 eventually. Set r=31500/36500=63/73>1/2. The exact dyadic identity Dseed=2048haell follows as before from e+L+12<=ell and8*((100ha+1)ell+e+L+12) between512haell and1024haell. It halves when a halves. The q rank term contracts in[1/3,1/2], while its entropy term contracts by r. Thus q_p/3<=q'<=r q_p and d_p/3<=d'<=r d_p.

The larger row coefficient increases some reserve terms. Its normalized cost still tends to zero. Along the parent and this witness, B=Theta(L/logz), h=Theta(L^2), e=Theta(L^2z), and

  (C_Y e z^2 P2(d))/(100Bhd) <= (C_Y/50)*e*z^2/(B*h)
    = O(z^3*logz/L) -> 0.

Therefore Y=(1+o(1))*100Bhd for both allocations, even with current C_Y=86000. Also k=(1+o(1))*hY: the remaining B*a term has ratio at most aell/(hd)<=1/(2h^2), and10d+C_Wez^2ell<=11d has ratio O(1/(Bh^2)). Consequently Y'<=r(1+o(1))Y_p and k'<=r(1+o(1))k_p. A fixed lower comparison k'>=k_p/40 follows termwise from B'>=B_p/6,d'>=d_p/3,a'=a_p/2,C_Y'>=C_Yp and C_W'/C_Wp>1/2.

Eventually k'<=68k_p/73,Y'<Y_p, and ceil_logk drops by at most6. With A_p=ceil_logk_p+E0_p+20 and A'>=A_p-6, A_p>=100 ensures (68k_p/73)/(A_p-6)<=k_p/A_p. Rounded PV block count does not increase. PV seed, condensed length, initial fullP2Y, per-call error,H,fielddegree and direct seed do not increase. Both INITIAL full candidates pass uniformly by E0>=L/125 and the strict8.9 advice exponent. Hence actual d0'<=d0_p; the witness preserves every capacity/range/5L/round/exact-error guard. Current selected k<=r(1+o(1))k_p.

The exact positive resource lower bounds now bootstrap current a=O(z^3),B=Theta(L/logz),d=Theta(L^4z^5),Y=Theta(L^7z^5/logz),k=Theta(L^9z^5/logz) before internal validity is used. Fixed positive C_W,C_Y do not change these orders. W/T have nprime=O(L^4z^5); rows have nprime=O(L^6z^8). The integer inequalities7^10>2^28 and19^10>2^42 imply fixed positive exponential margins over these two source orders with D>=10z/7. Thus their eventual ceiling roots are at most7 and19. It follows that

  H_WT<=7D<=10z+7<=1001z/100,
  H_rows<=19D<=190z/7+19<=136z/5,
  R=16D+16<=160z/7+32<=229z/10.

Since R=16D+16 and D tends to infinity, each field ratio s_i/e tends to64/15. The common maximum therefore satisfies s_star<4267e/1000 eventually. Each map keeps its own H,K,error and output metadata when field degree is raised. The full-request advice bounds are

  W/T: (978120143/62500)Reqez^2 < K/2=15750Reqez^2,
  rows: (132891448/3125)Reqez^2 < K/2=43000Reqez^2.

The positive gaps6254857/62500 and1483552/3125 timesReqez^2 pay the ACTUAL4rho_log+10, with rho_log=e+ceil_logK+13+ceil_logReq for that own K. This term is O(e+logK+logReq), not O(logK). Req*z^2 tends to infinity and logs O(z)=o(e), so its ratio to either gap tends to zero. K/4>=Req eventually. Every first full candidate passes ONEcall and common-field checks succeed.

PV seed<5e plus directseed<(4267/1000)(D+1)e gives total<(4267/700)ze+(6767/500)e<7ze eventually, paid byell>=8ze. Full stipulated seed bits remain retained; fixed output projection preserves strong error and binary linearity, with no discarded suffix observed. The same dimension and common field keep one universal paid primitive package. Complete P=O(L^2z^2)=o(k_p), so

  limsup (4P+8k)/(4P_p+8k_p) <= 63/73 < 1

at one fixed good paid code. All exact source/initial/robust error guards, canonical searches, bad-data and cutoff paths, streamed dimension runtime, empty-witness NP and the full16-table4n-8d gate proof retain their supplied arguments. The same-family divisor theorem, eventual weakenings and native(4,1) remain. Functions,fields,layouts and true minimum circuit sizes need not match the parent. Finite evidence separates initial-neighbor and actual per-component shared-port predicates; no complete enormous positive joint search or universal theorem is executed.
