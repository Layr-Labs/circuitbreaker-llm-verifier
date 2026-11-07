# One fixed total uniform family with quarter-log loss

For the evaluator F_N=family.evaluate on exactly N input bits, set L=ceil(log2 N), B=ceil(log2(L+1)), U=ceil(L^(3/4)), C=13, Z=U*ceil(B^(1/4)), E=C*Z, construction k=2^E and analytic dispersion dimension dstar=2^ceil(44119Z/10000) for N>=2. There exists N0 such that for EVERY integer N>=N0,

 C_B2(F_N)>=4N-8dstar+4.

Here dstar=2^((44119/10000+o(1))*L^(3/4)*(log2 L)^(1/4))=N^o(1). For this SAME fixed family, for EVERY fixed real epsilon>0 there is N_epsilon such that C_B2(F_N)>=4N-N^epsilon for EVERY integer N>=N_epsilon. The program and its time exponent do not depend on epsilon.

B2 includes all16 binary Boolean truth tables, arbitrary sharing and unrestricted depth. Every actual input bit is charged. All fields, parameters and choices are deterministic, without advice, field-polynomial input, oracle, promise or certificate. The family is total uniform P, hence NP with the verifier accepting precisely empty witnesses when its evaluator returns1.

The full-B2 coefficient4 is the existing Mishra--Santhanam theorem, completely supplied in elimination.md. The construction follows Shaltiel2011. Its new constant-error component follows the Shaltiel--Umans2005 interleaved-curve extractor reconstruction, explicitly proved here with the seed available to the predictor so as to establish strong extraction. Trevisan's design reconstruction supplies the separate small-error helper; Rao's reduction supplies the affine row condenser. The cubic ingredient follows Ben-Sasson--Kopparty and Haramaty--Shpilka. Prime-degree field construction follows Shoup; the retained analytic prime-interval proof follows Newman/Zagier. These are credited existing methods. This package does not claim discovery of the coefficient4 theorem or a stronger literature record. All nonstandard premises used below are supplied, not imported from URLs or prior verdicts.

The preceding tworegime7.61 submission6909b73a-5368-4953-a4fd-2575bed23bd9 supplies the pure finite-mode component interfaces, weighted depth potential and interleaved source framework. The preceding adaptive76157.7615 submissiona8c9dfad-4eec-472d-a621-775e12be4821 supplies the retuned SU/helper, declining-entropy and interleaved-source framework. The preceding adaptive7637.763 submissioneecdce4c-621f-4691-bc9e-8ccae3c2a0b9 supplies the refined active-gap/confidence and declining-entropy source framework. The preceding adaptive7657.765 submissionc8d2e270-f7af-4d2a-994f-bd499ed91046 supplies the retuned components, declining entropy schedule and complete interleaved-source argument. The immediate adaptiveentropy7.78 submissiona8a6788c-e829-43b9-a6a7-067faacb3923 supplies declining per-depth thresholds, proportional arities and interleaved source restrictions. The preceding local837.83 submissionba28d831-631a-490a-bbb7-5724df4a2ade supplies the retuned SU/helper and local-depth framework. Earlier local847.84 submission365b18a5-45be-4ee4-bc57-905a7dec3401 and localretune7.85 submission3f27447f-decd-4010-9eae-d28ee6019ce2 supply the refined local arity, gap and component framework. Earlier localfanout9.88 submissionee6f816a-a382-473f-8cfe-a1432ae0eb41 supplies local component lengths, varying arity and the discrete potential. Earlier scaledgap9.99 submissionf8d14c29-a0b4-44f9-87f9-e0e41060f66d supplies the scaled-gap/short-leaf framework and complete components. Earlier retuned81 submission797a6131-d919-48b3-9adc-d37dfac6c93d supplies the full SU/helper/design/confidence framework. Earlier12.82 submission9eff8e17-bd9c-4430-9fbc-b57b1e940cd9,12.83 submission2cf52c7e-3a94-4daa-af5f-d9bf954c06db,12.84 submission297eabee-7fc1-4965-93ca-2f63c08d2297 and12.89 submissionfb0f619f-e7a4-427f-9f99-884fabce6e4e supply the retuned components, thresholds and separated analytic-dimension framework. Earlier13 submission3b3f37d2-7a79-49e7-9a86-94c5eb25f8a4 supplies the13-scale evaluator and complete sharper ordinary-source framework. Earlier14 submissionfc7dbc37-e987-4acc-9e1b-a668bfa0cf14 supplies the raised relative component-threshold framework. Earlier20 submission1c2133c3-eda3-43c3-87e7-c715241e46b0 supplies the finer helper/design/confidence framework. Earlier22 submission96f9752d-967c-4740-8e4f-b3c015e8b667 supplies the sharper prime-degree design estimate. Earlier30 submission4403dc7f-49e8-4060-aab3-5309781859f5 supplies the enlarged smaller-output helper and raw-reservoir balance. Preceding46 submissiond46242dd-8e96-4c95-b7fa-cec67e494e08,68 submission42c2d818-9c1e-4a22-bb28-9e1497d77801 and156 submission875e388c-de8c-4349-b538-562dd42aabfe supply the augmented row, fanout and rational SU/source frameworks. Corrected170 submissionf52757ed-faf4-4adf-873e-8fabef827d1c and original170 submissionf9ac2386-70ff-4cf7-b67b-e846204ec5fb supply the quantized grid. Corrected224af68190e-22b9-4b70-a87b-397dfe692edf is a consulted corrected reference sharing the original224 framework; the component-local correction was applied commonly after original224 rejection, with no new origin or causal route through170 inferred from AF. Original224 submission0d6b9757-9407-42b1-9677-feb0c20af62e,296 submission5abcbc66-bd77-4f01-9b66-9c0f8b22105d,336 submission07e114e5-fa5e-4b96-a079-cbce7e77612e,2048 submission6ab40f99-b8a7-41d9-aadb-2aaa864a692a,3072 replacement764612b3-6862-4e50-b958-89632a4fb0e3 and5120 submission9ff6c2f2-36e4-481b-ae6f-15da0135f777 supply the shared confidence/SU/source framework. Original156f49104cd-7581-4704-bff5-e264dc45309f and68fa6537a5-253f-476c-bc28-32c1288e810f were cancelled for the exact-definition repair. Failed/cancelled or provisional submissions are not mathematical premises. The earlier own5 submission4138000d-8f17-4740-9995-59ee9b9fc772 supplies the large-row condenser. The earlier own36 submission5c80725e-ac7d-47e8-a736-1ce338fd7edf supplies the retuned SU extractor. The earlier own120 submission2c2108c8-4bc1-4b20-8ded-58c2a4295f97 supplies the full-output row reduction. The earlier own168 submissionbb65320f-4cec-4526-8895-6a857ebaef26 supplies the finely charged source entropy analysis. The earlier own submission59e21d17-0da7-44e2-a4f6-ade8d1e69f92 supplies the fully charged active/hidden/final entropy budgets. The earlier own submission0aaaa563-b3e6-41fc-a8c6-8846ba5132c9 supplies the stronger active entropy analysis. The earlier own submission5a30b469-95fb-4a85-aeff-84aa3617a384 supplies the sharper PNT seed estimate. The earlier own submission4cd2c215-066c-41a4-814f-14273960727a supplies the explicit C=8192 slower-fanout analysis. The preceding submitted4c2e1fa2-f0cd-481e-929e-72536cf8c8b3 supplies the explicit square-root-log entropy analysis. The earlier79262846-138d-4327-8058-031bf03f5d8c supplies the generalized extractor parameter lemma and balanced tree proof. The earlier f10654cb-2d81-4a33-b929-707695ca4beb supplies the complete supplied seed-aware SU extractor and subpolynomial-entropy tree framework. The earlier own e6a47252-1b22-49cd-bcc1-1c5d3ba51ea1 and2547b8e7-92aa-466c-a684-3a91b8e64319 supply the tree/component framework, exact unequal-interval accounting, and generic4n-8d+4 transfer. The preceding threeregime submissione4cd9ed1-3891-4f3d-9057-367dd14e21fb supplies the retained public-rate component,weighted depth and interleaved-source framework. The preceding rategrid submissionca81655d-8f22-46ca-851d-7f1b4ac978e3 supplies the retained public-rate component,weighted depth and interleaved-source framework. The preceding continuousrate submission37f9dd90-d357-4850-a3a9-37a0006173c9 supplies the retained public-rate component,weighted depth and interleaved-source framework. The preceding widerate submission749b17fa-bda9-4659-b805-bb56fdbde4eb supplies the retained public-rate component,weighted depth and interleaved-source framework. The preceding finerate submission73c321bf-41d7-4a57-a6b6-0395d8cf5f83 supplies the retained public-rate component,weighted depth and interleaved-source framework. The preceding balancedrate submission735dd188-b4cb-46b8-9f08-ad2105470177 supplies the retained public-rate component,weighted depth and interleaved-source framework. The preceding refinedboxes submission7c277bc5-fff2-44e5-8d39-d5c714f84aca supplies the retained public-rate component,weighted depth and interleaved-source framework. The preceding heightscale submission1f66aaef-eb95-417a-9813-19a9231bbfd6 supplies the retained public-rate component,weighted depth and interleaved-source framework. The preceding lowrate submission8828f5de-5df2-435d-94ce-b34bd2dee433 supplies the retained public-rate component,weighted depth and interleaved-source framework. The preceding tightlow submissionee277847-cd12-4e32-91d3-f4dcf543fbda supplies the retained public-rate component,weighted depth and interleaved-source framework. The preceding efficientlow submission5f645a43-b8a9-4cb0-9541-8dc11355b983 supplies the retained public-rate component,weighted depth and interleaved-source framework. This candidate computes a public rational extraction rate directly from W_j,Z,with explicit uniform component margins and context-size bounds. Its complete new constants and finite box certificate are supplied below. No prior verdict or old quantitative theorem is transferred to this changed evaluator.

cubic.md supplies the complete field constructor, analytic prime-interval argument, cubic isolation/bias and vector extractor. extractors.md supplies local codes, both strong seeded-extractor regimes, the new Singer-cycle construction, affine row condensation and block-source extraction. This file supplies the complete tree and source analysis. family.py implements exact integer roots and total numerical fallbacks. The small component audit does not execute the asymptotically nontrivial full tree or establish any universal theorem.

The global tree logarithm is L=ceil(log2 N). At depth j the component ambient is n_j with L_j=ceil(log2 n_j), and its DISTINCT component-local logarithm is Lc_j=ceil(log2(n_j+1)). In the SU theorem in extractors.md its local L denotes Lc_j, so the degree uses Lc_j as implemented by optimal_parameters. Lc_j=L_j or L_j+1; no exact equality is claimed. The enforced seed bound uses L_j, uniformly justified by the strict coefficient margin and min L_j tending to infinity.

# A concrete tree disperser with local component lengths

At actual input length N>=2 put L=ceil(log2 N), B=ceil(log2(L+1)), Z=ceil(L^(3/4))*ceil(B^(1/4)), base=ceil(L^(1/4))*ceil(B^(3/4)), E=13Z, Wmin=ceil(527Z/500), Wmax=ceil(33743Z/10000), k=2^E, ell_log=ceil(7*ceil(sqrt L)*ceil(sqrt B)/100), ell=2^ell_log, Pmax=2^ceil(269991Z/2000000), and X=ceil(25927Z/25000). Put Lpad=ceil(Wmax^(20/19)), computed by the least integer whose19th power is at least Wmax^20. Before allocating component thresholds or2^Lpad require E<=floor(log2 N), Lpad<=floor(log2 N), Wmin^3>=(L+1)^2. Otherwise return0. These guards hold eventually and prevent disproportionate allocations on small inputs. Define output0 also at N=0,1.

Plan levels beginning V_0=0. At each planned depth j, define the deterministic local ambient length

 n_j=max(ceil(N/2^V_j),2^Lpad), L_j=ceil(log2 n_j),
 W_j=max(Wmin,Wmax-ceil(25001*j*B/100000)), K_j=2^W_j,
 beta_j=max(3/20,min(13/25,27Z/(50W_j))),
 mode_j=the reduced positive numerator/denominator pair of beta_j,
 c_(mode_j)=(9999/10000)(1-beta_j)/(1+4alpha_j*beta_j+1/100000),
 alpha_j=1/(1-beta_j)+1/100000,
 A_j=max(1,floor(7*c_(mode_j)*W_j*base*L/(100*Z*L_j))), t_j=2^A_j.

Since W_j>=Wmin=ceil(527Z/500)>27Z/(50*(13/25)),the upper clamp never applies and beta_j W_j>=27Z/50. Use this context throughout every extractor,helper,reservoir and block-source invocation at that depth. The complete rational parameter formulas are specified in extractors.md and computed by regime in family.py. Context is a pure explicit parameter,including every cache key. Before any variable integer powers,component constructors require max(pRate,qRate)<=100Lc_j and3/20<=pRate/qRate<=13/25. This guard holds for every tree context: before reduction pRate=27Z<26Wmin<=26L_j and qRate=50W_j<=50L_j,while clamped contexts(3,20),(13,25) also fit. Here W_j<=Wmax<=Lpad^(19/20)<=L_j follows from the clipping identity. Reducing by gcd never increases either integer. No source value,epsilon,advice or external field data selects the context.

If V_j+A_j<=floor(log2 N), append V_(j+1)=V_j+A_j and continue; otherwise stop at J=j. Retain this terminal A_J too. Thus every level through J has nonempty intervals, and N<2^(V_J+A_J). Let Amax=max_{0<=j<=J}A_j and tmax=2^Amax. Require10000*(X+Amax)<=10371Z, otherwise return0.

For every j<J build the supplied constant-error extractor at component input length n_j and threshold K_j. Its local degree uses Lc_j=ceil(log2(n_j+1)), DISTINCT from L_j and the global L. Require its full seed count d_j to satisfy d_j<=seedguard_(mode_j)*L_j, with the exact rational context guards in extractors.md and its output length to be at least ell*Pmax. Trim this linear map to an ell by Pmax rectangle, denoted E_j. Build the supplied affine block-source extractor at block length n_j and threshold K_j; let b_j be its prescribed block count, and require at least Pmax output bits. Every component's exact integer numerical guard must succeed; otherwise return0. These are parameter checks, not tests of extraction on particular inputs.

Let the preliminary active levels be0 and, for each positive integer multiple iX reached by the V_j, the FIRST depth j with V_j>=iX. Under the accepted gap guard Amax<X, these crossing depths are distinct. Write them a_0=0<a_1<... . Choose H to be the last preliminary active depth such that floor(N/2^V_H)^1000>=k, and retain precisely the active depths through H. If fewer than two exist, return0. Write their count minus one as H'. Require t_j>=b_j for j<H and sum_{j=0}^{H-1}binom(t_j,b_j)<=ell, otherwise return0.

At level j, node i is the actual half-open interval[floor(iN/2^V_j),floor((i+1)N/2^V_j)),0<=i<2^V_j. Its t_j children partition it EXACTLY. Bits outside the interval are never interpreted as additional target bits. For a component invocation at this node, put the interval's actual bits in the PREFIX of a vector of length n_j and append zeros. For its block-source rows, put each selected child's actual bits in a separate prefix-zero-padded vector of length n_j. These are injective linear embeddings of the respective actual block coordinates. They preserve all affine entropies and conditional entropies; no target domain is padded, ignored or restricted.

For u=0,...,H'-1 set Tconf=ceil(ell_log/2000000) and

 p_u^(1)=2^(4000002*(H'-u)*Tconf),
 p_u^(2)=2^(4000002*(H'-u)*Tconf+2000001*Tconf).

Require all p2<=Pmax, otherwise return0. The exact ratios p2/(ell*p1) and earlier p1/(ell*later p2) equal G=2^(2000001*Tconf-ell_log)>=2^Tconf. At a FIXED node all confidences use the SAME map E_j and SAME rectangle width Pmax. Different levels have different component ambient lengths and may use different E_j. The common response seed-count upper bound R=max_{j<H}2^d_j will suffice throughout the proof.

Compute challenge matrices bottom-up. Leaf matrices are empty. At an internal level-j node, test its children in order at confidence Pmax using E_j on the prefix-padded node input. Select the first PASS child, or the first child if none passes. Concatenate the selected child's challenge with the Pmax-bit rows BE_j on EVERY lexicographically ordered b_j-tuple of children, each child padded to n_j. The total number of rows is at most sum_{r=j}^{H-1}binom(t_r,b_r)<=ell and is deterministic at each depth. All tests pad the challenge matrix with zero rows to ell. Follow the resulting selected path. At active depth a_u compare the selected descendant at depth a_(u+1) using the two specified confidence tests with map E_(a_u). At the first first-confidence FAIL output the second-confidence result, encoded FAIL0/PASS1. If there is none output0. This is exactly family.evaluate.

All planned ambient lengths are <=N. The component input embeddings are linear and injective on their actual blocks. The tree has at most N nonempty nodes per depth and J<=L, hence O(N log N) nodes, a sufficient polynomial bound even with unary levels. The tuple-count guard bounds each enumeration by ell<=N eventually; otherwise ell itself is at most a fixed polynomial in N by its formula. Each response has at most(N+1)^12 seeds by its Q^6 finite-field guard, and all component evaluators have a fixed polynomial time exponent at lengths <=N. Planning and guards require at most L levels, only polynomial-size integers, and polynomial-time component constructors. Thus all successful evaluations and all fallback evaluations take time polynomial in actual N, with a fixed program and exponent independent of epsilon. The complete refined full-B2 transfer applies to dispersion at dstar=2^ceil(44119Z/10000), giving4N-8dstar+4 and every stated same-family4N-N^epsilon bound.

# Challenge-response test and its exact properties

The following test follows Shaltiel's affine challenge-response method. A matrix distribution is eta-somewhere-random if it is within total variation eta of a convex combination of distributions having at least one uniform row; different mixture components may have different good rows. Joint affinity is not required in this definition. Resiliency R of a challenge relative to an affine source X means the same eta-somewhere-random conclusion holds after fixing any possible value of any binary linear map of X of rank at most R.

At a fixed node of depth j, take its injectively prefix-padded input at ambient length n_j and its local threshold K_j, within the uniform component window proved below. The supplied constant-error linear extractor has m=floor(K_j^(beta_mode_j)) output bits, d_j=O(log n_j) seed bits and error1/4. Trim to ell by Pmax to obtain E_j. All confidences at this fixed node use the SAME E_j and SAME Pmax. Different nodes at different depths may use different maps. Set R=max_{j<H}2^d_j<=N^6, a common upper bound on every response count. The following four properties apply separately at each node; only the common bound R is used when combining them. Injective linear embeddings preserve every affine source and conditional entropy to which the properties are applied.

For an input x, coordinate block s, challenge function ch on that block with at most ell rows and Pmax columns, and1<=p<=Pmax, pad ch(x_s) with zero rows to ell rows. Return FAIL exactly when some enumerated response E_j(x,y) agrees with that challenge on the first p columns of every row. Otherwise return PASS. Prefix-zero-pad the actual node input to its specified n_j before the response map. Each test is deterministic polynomial time; challenges are computed by the containing algorithm.

The following properties hold.

1. If p1<=p2, FAIL at p2 implies FAIL at p1, by restricting any matching response to fewer columns.

2. If X is affine with H(X)>=K_j and X_s is fixed, there is an affine subsource of EXACT deficiency ell*p on which the test always fails. Strong linear extraction has error1/4, so at least one seed has a completely uniform ell-by-Pmax response on X: a nonuniform affine output would have distance at least1/2. Fix the ell*p specified response coordinates equal to the constant challenge slice. Uniformity makes this a nonempty affine fiber of codimension ell*p, of probability2^(-ell*p). On that fiber one matching response is guaranteed. The argument does not require ch itself to be linear.

3. More generally, if the affine conditional entropy H(X|X_s)>=K_j, FAIL has probability at least2^(-ell*p). Fix each possible X_s value. Its fiber remains affine and has entropy at least K_j; property2 applies to that fiber. Average the resulting lower bounds. The choice of responding seed may depend on the fiber, which is allowed because FAIL enumerates all seeds.

4. If ch(X_s) is eta-somewhere-random with resiliency ell*p relative to X, then FAIL has probability at most R*(2^-p+eta). For any fixed response seed y, the sliced response E_j is a linear map of X with rank at most ell*p. Condition on any of its possible values v. Resiliency makes the challenge eta-close to a mixture with a uniform row. Matching the whole slice to v then has probability at most2^-p+eta, since matching that component's good row alone has probability2^-p. Average over v, and union bound over all R seeds.

All source restrictions here are affine and have explicitly charged rank deficiency. Property3 is a rare-event lower bound, not a high-probability classification of arbitrary blocks. Property4 is a conditional somewhere-randomness hypothesis, not a claim that every nonempty block automatically passes. The final disperser proof must establish that hypothesis for its recursively chosen challenges and compare its error with the rare FAIL probabilities. Replacing R by n^2 without checking the actual seed length would be unjustified; only its polynomial-in-n bound is used in the intended asymptotic parameter comparisons.


# Correctness proof for the subpolynomial-entropy tree construction

This proof uses the exact tree construction above, implemented by family.py, and the complete component arguments in cubic.md, extractors.md and cubic.md. The construction and the two-stage source analysis follow Shaltiel2011; the extractor components are the separately supplied replacements.

Write H for tree depth and H' for the number of consecutive active intervals. All statements in this proof concern sufficiently large actual N, with n=N,L,k,local K_j,ell,Pmax and local arities/ambients as defined in the construction. An affine subsource means conditioning a uniform affine source on possible values of binary linear maps. Its deficiency is the rank lost, and its probability inside the original source is2^(-deficiency). Linear images stay affine. Conditioning on a rank-D map decreases any linear-image entropy, and any conditional entropy of linear images, by at most D. To verify the conditional assertion, write H(A|B)=dim image(A,B)-dim image(B); the first term drops by at most D and the second cannot increase. Entropy chain rules are therefore exact dimension identities in this setting.

## Numerical feasibility with uniformly bounded rational rates

Write S=13Z,w0=33743/10000=3.3743,wm=527/500=1.054,a=27/50=.54,delta=25001/100000. Wmin=ceil(wm Z),Wmax=ceil(w0 Z),W_j=max(Wmin,Wmax-ceil(delta jB)). Uniformly over planned depths for sufficiently large actual N,Z=(1+o(1))L^(3/4)B^(1/4),base=(1+o(1))L^(1/4)B^(3/4),ell_log=(7/100+o(1))sqrt(LB),X=(25927/25000+o(1))Z. Lpad=ceil(Wmax^(20/19))=o(L),and Wmin^3>>(L+1)^2.

The exact identity L_j=max(L-V_j,Lpad),with Lc_j=L_j or L_j+1,gives (Lc_j)^(2/3)<=W_j<=(Lc_j)^(19/20). The explicit uniform rational-rate component proof therefore applies at every level; the context-size guard was proved above. Full-seed coefficient at most5.333355 has an exact guard gap1/100000,and guard at most5.333365<6. Q,h^d<=(n_j+1)^6 hold uniformly. Response seed count is at most N^6 eventually,and at most(N+1)^12 by numerical guards at every accepted length. Every response/block output has logarithm at least beta_j W_j-O(1)>=.54Z-O(1).

For all allowed contexts .08<c<.5. The fixed height multiplier7/100 gives (7/100)c<.035<.06. Let M=Amax include terminal arity. Then

 M<=.06(Wmax/Z)base L/Lpad+1=O(L^(35/76)B^(37/76)),

so M=o(Z),M^2=o(base L),M<X,and X+M<=1.0371Z eventually. Consecutive active log gaps are at most X+M,and H'ell_log<=(1750/25927+o(1))Z. First active crossing has shortest block log L-X-M+O(1)>>S/1000. If another active depth follows H,maximality gives maxleaf<=2^(X+M)(k^(1/1000)+1)+1. Otherwise V_J-V_H<X and N<2^(V_J+A_J),giving maxleaf<=2^(X+M)+1. Both cases imply maxleaf<=2^((1.0371+.013+o(1))Z)<2^(1.051Z),strictly below hidden entropy2^(1.0542Z).

For u_j=L-V_j>=Lpad,the exact floor gives

 A_j>=(7/100)c_j W_j base L/(Z u_j)-1,
 u_j^2-u_(j+1)^2>=(7/50)c_j(W_j/Z)base L-2L-M^2.

Below clipping,each step is at least floor(.005Wmin base L/(Z Lpad)),so only O(Lpad^2/(base L))+1=o(L/base) additional steps occur. Put h0=(w0-wm)/delta=232030/25001. Here wm>a/(13/25),so the upper clamp never applies. The strict bound w0<a/(3/20) ensures that the lower clamp also never applies for sufficiently large N. Thus beta_j=aZ/W_j throughout every planned level eventually;clamping still defines the total algorithm at finite exceptions.

For wm<=w<=w0 put beta(w)=a/w,alpha(w)=1/(1-beta(w))+1/100000,

 c(w)=(9999/10000)(1-beta(w))/(1+4alpha(w)beta(w)+1/100000),
 f(w)=w*c(w).

The function f is increasing: beta decreases with w;1-beta increases,and the positive denominator1+4beta/(1-beta)+4beta/100000+1/100000 decreases. Thus c and w are positive increasing factors. No differential equation or numerical integration is assumed.

Use4000 FIXED boxes. Set Delta=(w0-wm)/4000=23203/40000000,w_r=wm+rDelta for0<=r<4000,and A_r=floor(10^12*f(w_r)). These are explicit rational/integer expressions independent of N. Exact finite arithmetic gives

 SUM_(r=0)^3999 A_r=3080208477059576,
 (7Delta/(50delta))*SUM A_r/10^12
   =62536317631561674187/62502500000000000000>1.

The declared audit recomputes every rational value,checks each floor by both integer inequalities,and verifies this exact sum and ratio. This is a finite integer certificate for a displayed constant,not empirical evidence for a universal tree claim.

Suppose J>=(23203/2500)L^(3/4)/B^(3/4). Since23203/2500>h0,the first m=ceil(h0*Z/B) planned indices then exist for all sufficiently large N. Before deleting clipped indices,each fixed box contains Delta*Z/(delta B)+O(1) indices of the raw schedule Wmax-ceil(delta jB);the possible final floor point changes only O(1) terms. On a retained index in[w_r,w_(r+1)),monotonicity gives f(W_j/Z)>=f(w_r)>=A_r/10^12. Deleting the o(L/base) clipped indices changes the normalized sum by o(1),because f is uniformly bounded on this fixed interval. Since base*Z/(LB)->1,and the arity multiplier is7/100,the remaining normalized potential drop is at least the displayed coefficient minus o(1),which exceeds1 and contradicts initial potential L^2. Count-rounding errors across4000 fixed boxes,and the terms M^2=o(base L),also total o(1). Thus H<=J<(23203/2500)L^(3/4)/B^(3/4)=9.2812L^(3/4)/B^(3/4) eventually.

At each depth,log g_j=(1-beta_j+o(1))W_j and b_j<=2+d_j/log g_j. Because (7/100)*base*L/Z=(1+o(1))*ell_log,the exact coefficient identity gives

 b_j A_j<=[c_j*seedguard_j/(1-beta_j)+o(1)]ell_log+2M
          =(.9999+o(1))ell_log+2M<.99995ell_log.

Therefore sum_j binom(t_j,b_j)<=ell eventually. Since1-beta_j>=.48,b_j=O(L/Wmin) uniformly,so log2(4b_j)<=delta B and

 log2 product_(j<H)(4b_j)<=(23203/2500)delta Z
   =580098203Z/250000000<2.3204Z.

Target dstar>=2^(4.4119Z),active log gap<=1.0371Z,and log2(4^H')=o(Z) give active stopping entropy at least2^(3.3747Z). A nice source at ordinary depth j descended from active depth jstart retains at least2^(3.3747Z-delta(j-jstart)B)-Dtot. Its log margin against W_j's declining branch is .0004Z-O(1);against the floor it is at least(3.3747-2.3204-1.054)Z-O(1)=.0003Z-O(1). One further ordinary division costs delta B=o(Z),and hidden-child division costs at most M+2=o(Z). Every interleaved forcing and hidden source therefore has PER-NODE entropy at least2^(W_j+.0001Z),after all charged ranks. No global minimum substitutes for a larger earlier threshold.

Forcing rank Dtot<=H*tmax*ell*Pmax<=2^(.137Z),and later resilience rank(H+1)ell Pmax<=2^(.137Z),are negligible relative to every K_j. Conditional entropy H(X_g|X_g') stays above2^(3.3746Z)>K_g;globally hidden entropy stays above2^(1.0542Z)>maxleaf2^(1.051Z).

For every context,eBE=beta_j/4-102/100000000. Since beta_j W_j>=aZ and W_j<=w0 Z+1,

 eBE W_j>=aZ/4-(102/100000000)(w0 Z+1)
          =.134996558214Z-.00000102>=.134996Z eventually.

Each block error is at most2^(-2^(.134996Z))<2^-Pmax. Output logarithm.54Z-O(1) exceeds log2(ell Pmax). Since(2000001/1000000)(1750/25927)<.134995<.1349955,log2(ell^3 p2)<=.134995Z<log2 Pmax eventually. Confidence grows faster than every fixed polynomial in L;R<=N^6 gives H log2(2R+2)=O(L^2)=o(Pmax). All uniform component/helper/row/rare-event/numerical guards hold for all sufficiently large actual integer lengths.

## Locate useful nodes and force earlier siblings during descent

Fix any affine source X on actual N coordinates with entropy at least dstar. A source is nice for interval s at threshold K if coordinates before s are fixed and H(X_s)>=K. Initially the root is nice at dstar. All conditional entropies below are exact affine dimensions on every possible fiber.

For a nice block partitioned into T ordered children, let k_i be each child's conditional entropy given all earlier children; their sum is at least K. Select the first child with k_i>=K/(4T), and fix its preceding children to a possible joint value. Those earlier children contributed less than K/4. If the selected child's entropy is at least K/4, continue at threshold K/4. Otherwise the remaining conditional entropy in the parent outside that child is at least K/2. At this stopping step, both the child entropy and H(parent|child) are at least K/(4T).

Apply this split on consecutive ACTIVE levels, with T<=2^(1.0371Z). Continuing to the final active leaf is impossible: the leaf would have entropy at least dstar/4^H', while its length is below2^(3Z). At a stopping step we therefore obtain consecutive active nodes g,g', an affine subsource, and fixed coordinates before g', with

 H(X_g')>=2^(3.3747Z), H(X_g|X_g')>=2^(3.3747Z).

Indeed the stopping threshold is at least dstar/(4^H'4*2^(1.0371Z))>=2^((4.4119-1.0371-o(1))Z)>2^(3.3747Z). Every subsequent prefix-fixing operation inside g' preserves H(X_g|X_g') exactly: conditioning on functions of X_g' selects its affine fibers without changing their common conditional dimension. Every additional response constraint of rank D can decrease that dimension by at most D.

We now interleave the ordinary entropy descent and the forcing of earlier siblings. This order is necessary because component entropy thresholds depend on depth. Define the common prospective rank budget Dtot=H*tmax*ell*Pmax. Initially force all earlier children on each ancestor of g'. Their actual input coordinates are fixed, so their challenge matrices are constant. At a parent of depth j the affine parent projection has entropy at least H(X_g') minus ranks spent so far, which exceeds K_j<=2^Wmax. The empty-block challenge-response property with its LOCAL map E_j therefore permits fixing exactly ell*Pmax response coordinates to make that child test FAIL at confidence Pmax. Pull this affine fiber back to the full source. Its deficiency is the same rank, and all previous FAIL events persist pointwise. This costs at most Dtot over these pairs and the later pairs below. It preserves both useful entropies up to the charged ranks.

Start the ordinary descent at s=g' of depth jstart. At depth j let F_j be the product of4b_r over ordinary steps already taken. Maintain a nice source with

 H(X_s)>=2^(3.3747Z)/F_j - Dspent,
 F_j<=2^(delta(j-jstart)B), Dspent<=Dtot.

Suppose at least b_j children have conditional entropy at least H(X_s)/(4t_j). Those children, in order, form a hidden b_j-block affine block-wise source: conditioning on fewer predecessor blocks cannot decrease conditional entropy. Stop at this node e. Otherwise fewer than b_j indices are high; the low indices together contribute at most H(X_s)/4. Some high child has conditional entropy at least H(X_s)/(4b_j). Fix its preceding children and move into it. Coordinates before the chosen child are now fixed and H(X_g|X_g') is preserved, since the new prefix restrictions are inside g'.

BEFORE the next ordinary move, force every earlier child at this parent to FAIL using E_j and confidence Pmax. At this moment the parent entropy is at least the just-retained child entropy, and the per-node estimates above give

 H(parent)>=2^(3.3747Z-delta(j-jstart+1)B)-Dspent>K_j.

The extra one-step cost delta B=o(Z) fits within the fixed .0001Z per-node reserve. Each earlier child's challenge is constant because that child's coordinates were fixed. Thus the exact-deficiency property applies while this parent still has enough entropy for its particular K_j. Charge its at most t_j*ell*Pmax ranks; previously forced FAIL events persist. The nice-child invariant survives: division by4b_j scales earlier rank losses DOWN, and adding new rank losses gives the conservative bound by total Dspent displayed above. Response constraints can reduce H(X_g|X_g') only by their charged ranks. This completes the induction for the interleaved procedure.

It cannot reach a leaf. The global retained threshold is greater than2^((3.3747-2.3204)Z)-Dtot>2^(1.0542Z), larger than every leaf length. Therefore it stops at an INTERNAL node e of depth je with a hidden b_je-block affine block-wise source. The additional division4t_je has log at most M+2=o(Z), so the stronger PER-NODE estimate gives its hidden-block threshold at least2^(W_je+.0001Z), after all forcing ranks. Let X* denote this final affine source. It satisfies

 e subset g' subset g; g,g' consecutive active nodes;
 H(X*_g|X*_g')>=2^(3.3746Z)>K_(depth(g));
 the selected tuple at e is block-wise at threshold2^(W_je+.0001Z);
 all coordinates before e are fixed;
 every earlier-sibling test on every ancestor of e FAILs identically.

Ancestors above g' were handled initially; each lower ancestor was handled at its ordinary move. There are at most H*tmax pairs in total, so the claimed budget Dtot<=2^(.137Z) holds. Later rank conditioning of at most(H+1)ell Pmax cannot exhaust the hidden tuple's per-node margin, hence its block-wise threshold remains at least K_je on every possible full-history fiber. The same reasoning preserves the stated conditional entropy at g.

We do NOT require that X* has entropy at least K_j in every earlier ancestor: a larger threshold may have been used earlier in the process. Forcing was applied while its entropy condition held, and its resulting deterministic FAIL event persists. The upward PASS bounds below use only linearity of each E_j, the uniform response count bound and challenge resiliency. The final rare FAIL lower bound is used at g, where the separate conditional-entropy bound above explicitly exceeds its own K_(depth(g)).


## Propagate resilient challenges upward

Put R=max_{j<H}2^d_j<=N^6. For an ancestor s of e at distance a from e define

 eta_a=(2R+2)^a *2^-Pmax,
 D_a=(H-a+1)*ell*Pmax.

We claim ch_s(X*_s) is eta_a-somewhere-random with resiliency D_a relative to the FULL source X*. Resiliency is the rank-conditioning definition in the challenge-response section of proof.md.

For a=0, conditioning on any rank at most(H+1)*ell*Pmax leaves the hidden tuple at e above2^(W_je+.0001Z)-(H+1)*ell*Pmax>=K_je. Node e at depth j lists BE_j on every b_j-tuple of children, using injective prefix padding to n_j. On the correct tuple one listed row is within2^-Pmax of uniform. Appending arbitrary correlated rows preserves that approximation: couple the row to its ideal uniform distribution and extend the other rows using their conditional distribution. This proves the base case for the whole challenge, whether or not its recursively selected child was useful.

For the induction step let z be the child of s toward e. Fix any rank-D_a conditioning of X*. The inductive guarantee for ch_z leaves an additional ell*Pmax of resiliency, and its error remains eta_(a-1). Challenge-response therefore tests z with FAIL probability at most

 delta=R*(2^-Pmax+eta_(a-1)).

Every earlier child of s FAILs identically, even after this extra conditioning. Thus, except for probability delta, s selects z and its challenge contains the entire challenge of z in the prescribed leading rows. Replacing the exceptional choices changes the distribution by at most delta. The child challenge's error and arbitrary appended rows add at most eta_(a-1). Their sum is at most(2R+1)*eta_(a-1)<=eta_a. This proves the induction. Variable mixture rows cause no issue: the definition needs only a uniform row in each component, and row count at a given level is fixed.

Uniformly for every a<=H,

 eta_a<=2^(-Pmax/2)

eventually, since H*log2(2R+2)=O(L^2)=o(Pmax). In particular every true-path challenge has resiliency at least ell*Pmax. Applying the PASS bound at each edge to e and union bounding shows the selected path reaches e, and hence g', with probability at least1-epsilon_path, where

 epsilon_path<=H*R*(2^-Pmax+eta_H)<=2^(-Pmax/3).

All probabilities are with respect to the single fixed affine subsource X*. Left siblings never interfere; the true child is chosen whenever its test passes.

## Both output values occur with the quantized confidence grid

Let g be at active depth a_u. Define fixed-node tests T_a=CR_(ch_g',p_u^(a))(X*_g,g') for a=1,2. When the selected path reaches g', these are exactly g's output tests. Write p1,p2 for the two confidences, Eell=log2 ell, Tconf=ceil(Eell/2000000), D=ell*p2, and G=2^(2000001*Tconf-Eell). The exact integer grid gives p2/(ell*p1)=G>=2^Tconf, tending to infinity uniformly over all active levels.

Conditional entropy H(X*_g|X*_g')>=2^(3.3746Z)>K_(depth(g)) gives Pr[T_a=FAIL]>=2^(-ell*p_a). The resilient PASS bound gives Pr[T2=FAIL]<=R*(2^-p2+eta_H). Monotonicity makes T2=FAIL imply T1=FAIL. Thus

 Pr[T1=FAIL,T2=FAIL]>=2^-D,
 Pr[T1=FAIL,T2=PASS]>=2^(-ell*p1)-R*(2^-p2+eta_H)
                       >=(1/2)*2^(-ell*p1)>=2^-D

eventually. The gap G diverges, log R=O(log n) is negligible compared with the smallest p2, and Pmax dominates every confidence exponent. The smallest p2 is at least2^(6000003*Tconf), whose logarithm is((6000003/2000000)*(7/100)+o(1))*sqrt(LB), so that p2 grows faster than every fixed polynomial in L.

For any strict active ancestor of g, the first-confidence logarithm is at least4000002*(H'-u+1)*Tconf. Subtracting log2 D=Eell+4000002*(H'-u)*Tconf+2000001*Tconf leaves2000001*Tconf-Eell=log2 G. Its first confidence is therefore at least G*D. Since G<=ell^2 eventually (indeed log2 G=(1/2000000+o(1))*Eell), we have G*D<=ell^3*p2. The numerical bounds log2(ell^3*p2)<=(26999/200000)*Z and log2 Pmax>=(269991/2000000)*Z show Pmax/(G*D) tends to infinity uniformly over all active levels.

Every strict-ancestor test has the specified resiliency. Union bounding FAIL gives

 H'*R*(2^(-G*D)+eta_H)<=2^(-sqrt(G)*D)

eventually: eta_H<=2^(-Pmax/2), the strict gap269991/2000000-26999/200000 pays for that term, and log(H'R)=O(log n) is negligible compared with G*D. If g is the root there are no ancestors and that probability is zero.

The probability of a deviation from the true path to e is at most2^(-Pmax/3). Its sum with2^(-sqrt(G)*D) is below2^-D eventually. Each output event above therefore has a point whose selected path reaches g' and whose earlier active nodes all PASS their first test. On these respective points the evaluator chooses g as the first active FAIL and outputs FAIL=0 and PASS=1. Thus both bits occur on X*, a subsource of the arbitrary affine source X of entropy at least dstar, proving dispersion for the same fixed family at every sufficiently large integer length. No uncharged node choice, extra output or input promise is used.

## Actual-length full-B2 consequence

The tree already acts on all N input bits, without outer padding. Its certified dispersion dimension is dstar=2^ceil(44119Z/10000)<=N eventually. The evaluator still uses construction k=2^(13Z), now with fixed depth-dependent component thresholds K_j=2^W_j. The complete refined full-B2 theorem in elimination.md therefore gives

 C_B2(F_N)>=4N-8dstar+4.

Since log2 dstar=(44119/10000+o(1))*L^(3/4)*(log2 L)^(1/4), for every fixed epsilon>0 we have8dstar-4<=N^epsilon at every sufficiently large integer N. Thus the SAME fixed deterministic uniform P family satisfies C_B2(F_N)>=4N-N^epsilon eventually, separately for every fixed epsilon>0. No exponent or program is selected in response to epsilon. Empty-witness NP verification and actual-length accounting hold for this same family. No explicit threshold for the full dispersion/circuit theorem is asserted.
