# One fixed total uniform family with quarter-log loss

For the evaluator F_N=family.evaluate on exactly N input bits, set L=ceil(log2 N), B=ceil(log2(L+1)), U=ceil(L^(3/4)), C=13, Z=U*ceil(B^(1/4)), E=C*Z, construction k=2^E and analytic dispersion dimension dstar=2^ceil(196Z/25) for N>=2. There exists N0 such that for EVERY integer N>=N0,

 C_B2(F_N)>=4N-8dstar+4.

Here dstar=2^((196/25+o(1))*L^(3/4)*(log2 L)^(1/4))=N^o(1). For this SAME fixed family, for EVERY fixed real epsilon>0 there is N_epsilon such that C_B2(F_N)>=4N-N^epsilon for EVERY integer N>=N_epsilon. The program and its time exponent do not depend on epsilon.

B2 includes all16 binary Boolean truth tables, arbitrary sharing and unrestricted depth. Every actual input bit is charged. All fields, parameters and choices are deterministic, without advice, field-polynomial input, oracle, promise or certificate. The family is total uniform P, hence NP with the verifier accepting precisely empty witnesses when its evaluator returns1.

The full-B2 coefficient4 is the existing Mishra--Santhanam theorem, completely supplied in elimination.md. The construction follows Shaltiel2011. Its new constant-error component follows the Shaltiel--Umans2005 interleaved-curve extractor reconstruction, explicitly proved here with the seed available to the predictor so as to establish strong extraction. Trevisan's design reconstruction supplies the separate small-error helper; Rao's reduction supplies the affine row condenser. The cubic ingredient follows Ben-Sasson--Kopparty and Haramaty--Shpilka. Prime-degree field construction follows Shoup; the retained analytic prime-interval proof follows Newman/Zagier. These are credited existing methods. This package does not claim discovery of the coefficient4 theorem or a stronger literature record. All nonstandard premises used below are supplied, not imported from URLs or prior verdicts.

The immediate preceding localretune7.85 submission3f27447f-decd-4010-9eae-d28ee6019ce2 supplies the fully retuned SU/helper and local-component source framework. This candidate refines local fanout and active spacing, with all exact source budgets rederived for analytic loss196/25. Earlier localfanout9.88 submissionee6f816a-a382-473f-8cfe-a1432ae0eb41 supplies local component lengths, varying arity and the discrete potential. Earlier scaledgap9.99 submissionf8d14c29-a0b4-44f9-87f9-e0e41060f66d supplies the scaled-gap/short-leaf framework and complete components. Earlier retuned81 submission797a6131-d919-48b3-9adc-d37dfac6c93d supplies the full SU/helper/design/confidence framework. Earlier12.82 submission9eff8e17-bd9c-4430-9fbc-b57b1e940cd9,12.83 submission2cf52c7e-3a94-4daa-af5f-d9bf954c06db,12.84 submission297eabee-7fc1-4965-93ca-2f63c08d2297 and12.89 submissionfb0f619f-e7a4-427f-9f99-884fabce6e4e supply the retuned components, thresholds and separated analytic-dimension framework. Earlier13 submission3b3f37d2-7a79-49e7-9a86-94c5eb25f8a4 supplies the13-scale evaluator and complete sharper ordinary-source framework. Earlier14 submissionfc7dbc37-e987-4acc-9e1b-a668bfa0cf14 supplies the raised relative component-threshold framework. Earlier20 submission1c2133c3-eda3-43c3-87e7-c715241e46b0 supplies the finer helper/design/confidence framework. Earlier22 submission96f9752d-967c-4740-8e4f-b3c015e8b667 supplies the sharper prime-degree design estimate. Earlier30 submission4403dc7f-49e8-4060-aab3-5309781859f5 supplies the enlarged smaller-output helper and raw-reservoir balance. Preceding46 submissiond46242dd-8e96-4c95-b7fa-cec67e494e08,68 submission42c2d818-9c1e-4a22-bb28-9e1497d77801 and156 submission875e388c-de8c-4349-b538-562dd42aabfe supply the augmented row, fanout and rational SU/source frameworks. Corrected170 submissionf52757ed-faf4-4adf-873e-8fabef827d1c and original170 submissionf9ac2386-70ff-4cf7-b67b-e846204ec5fb supply the quantized grid. Corrected224af68190e-22b9-4b70-a87b-397dfe692edf is a consulted corrected reference sharing the original224 framework; the component-local correction was applied commonly after original224 rejection, with no new origin or causal route through170 inferred from AF. Original224 submission0d6b9757-9407-42b1-9677-feb0c20af62e,296 submission5abcbc66-bd77-4f01-9b66-9c0f8b22105d,336 submission07e114e5-fa5e-4b96-a079-cbce7e77612e,2048 submission6ab40f99-b8a7-41d9-aadb-2aaa864a692a,3072 replacement764612b3-6862-4e50-b958-89632a4fb0e3 and5120 submission9ff6c2f2-36e4-481b-ae6f-15da0135f777 supply the shared confidence/SU/source framework. Original156f49104cd-7581-4704-bff5-e264dc45309f and68fa6537a5-253f-476c-bc28-32c1288e810f were cancelled for the exact-definition repair. Failed/cancelled or provisional submissions are not mathematical premises. The earlier own5 submission4138000d-8f17-4740-9995-59ee9b9fc772 supplies the large-row condenser. The earlier own36 submission5c80725e-ac7d-47e8-a736-1ce338fd7edf supplies the retuned SU extractor. The earlier own120 submission2c2108c8-4bc1-4b20-8ded-58c2a4295f97 supplies the full-output row reduction. The earlier own168 submissionbb65320f-4cec-4526-8895-6a857ebaef26 supplies the finely charged source entropy analysis. The earlier own submission59e21d17-0da7-44e2-a4f6-ade8d1e69f92 supplies the fully charged active/hidden/final entropy budgets. The earlier own submission0aaaa563-b3e6-41fc-a8c6-8846ba5132c9 supplies the stronger active entropy analysis. The earlier own submission5a30b469-95fb-4a85-aeff-84aa3617a384 supplies the sharper PNT seed estimate. The earlier own submission4cd2c215-066c-41a4-814f-14273960727a supplies the explicit C=8192 slower-fanout analysis. The preceding submitted4c2e1fa2-f0cd-481e-929e-72536cf8c8b3 supplies the explicit square-root-log entropy analysis. The earlier79262846-138d-4327-8058-031bf03f5d8c supplies the generalized extractor parameter lemma and balanced tree proof. The earlier f10654cb-2d81-4a33-b929-707695ca4beb supplies the complete supplied seed-aware SU extractor and subpolynomial-entropy tree framework. The earlier own e6a47252-1b22-49cd-bcc1-1c5d3ba51ea1 and2547b8e7-92aa-466c-a684-3a91b8e64319 supply the tree/component framework, exact unequal-interval accounting, and generic4n-8d+4 transfer. This version uses useful ACTIVE-node entropy k^(471/1300), hidden block-source entropy k^(4003/13000), and final retained entropy k^(4001/13000)>Kmin=2^(4Z), then balances the tree degree, entropy and response confidence parameters, reducing the LOSS SCALE from5*L^(3/4)*sqrt(log2 L) to(196/25)*L^(3/4)*(log2 L)^(1/4), by decoupling active spacing from challenge-row count and supplying the stronger precision/error budgets for ONE fixed family. The seed-aware extractor and all component arguments are supplied again. The earlier2d32f686-6419-451e-8758-22dd489737e0 is a comparator and research motivation: its different algebraic selector has an explicitly positive-power loss, whereas this new family has the supplied subpolynomial-loss theorem. No earlier result is admitted as a proof premise.

cubic.md supplies the complete field constructor, analytic prime-interval argument, cubic isolation/bias and vector extractor. extractors.md supplies local codes, both strong seeded-extractor regimes, the new Singer-cycle construction, affine row condensation and block-source extraction. This file supplies the complete tree and source analysis. family.py implements exact integer roots and total numerical fallbacks. The small component audit does not execute the asymptotically nontrivial full tree or establish any universal theorem.

The global tree logarithm is L=ceil(log2 N). At depth j the component ambient is n_j with L_j=ceil(log2 n_j), and its DISTINCT component-local logarithm is Lc_j=ceil(log2(n_j+1)). In the SU theorem in extractors.md its local L denotes Lc_j, so the degree uses Lc_j as implemented by optimal_parameters. Lc_j=L_j or L_j+1; no exact equality is claimed. The enforced seed bound uses L_j, uniformly justified by the strict coefficient margin and min L_j tending to infinity.

# A concrete tree disperser with local component lengths

At actual input length N>=2 put L=ceil(log2 N), B=ceil(log2(L+1)), Z=ceil(L^(3/4))*ceil(B^(1/4)), base=ceil(L^(1/4))*ceil(B^(3/4)), E=13Z, W=4Z, k=2^E, Kmin=2^W, ell_log=ceil(sqrt L)*ceil(sqrt B), ell=2^ell_log, Pmax=2^ceil(1279Z/2000), and X=ceil(391Z/125). Put Lpad=ceil(W^(20/19)), computed by the least integer whose19th power is at least W^20. Before allocating2^W or2^Lpad require E<=floor(log2 N), Lpad<=floor(log2 N), W^3>=(L+1)^2. Otherwise return0. These guards hold eventually and prevent disproportionate allocations on small inputs. Define output0 also at N=0,1.

Plan levels beginning V_0=0. At each planned depth j, define the deterministic local ambient length

 n_j=max(ceil(N/2^V_j),2^Lpad), L_j=ceil(log2 n_j),
 A_j=max(1,floor(443*base*L/(2500*L_j))), t_j=2^A_j.

If V_j+A_j<=floor(log2 N), append V_(j+1)=V_j+A_j and continue; otherwise stop at J=j. Retain this terminal A_J too. Thus every level through J has nonempty intervals, and N<2^(V_J+A_J). Let Amax=max_{0<=j<=J}A_j and tmax=2^Amax. Require1000*(X+Amax)<=3129Z, otherwise return0.

For every j<J build the supplied constant-error extractor at component input length n_j and threshold Kmin. Its local degree uses Lc_j=ceil(log2(n_j+1)), DISTINCT from L_j and the global L. Require its full seed count d_j to satisfy20000*d_j<=162223*L_j and its output length to be at least ell*Pmax. Trim this linear map to an ell by Pmax rectangle, denoted E_j. Build the supplied affine block-source extractor at block length n_j and threshold Kmin; let b_j be its prescribed block count, and require at least Pmax output bits. Every component's exact integer numerical guard must succeed; otherwise return0. These are parameter checks, not tests of extraction on particular inputs.

Let the preliminary active levels be0 and, for each positive integer multiple iX reached by the V_j, the FIRST depth j with V_j>=iX. Under the accepted gap guard Amax<X, these crossing depths are distinct. Write them a_0=0<a_1<... . Choose H to be the last preliminary active depth such that floor(N/2^V_H)^1000>=k, and retain precisely the active depths through H. If fewer than two exist, return0. Write their count minus one as H'. Require t_j>=b_j for j<H and sum_{j=0}^{H-1}binom(t_j,b_j)<=ell, otherwise return0.

At level j, node i is the actual half-open interval[floor(iN/2^V_j),floor((i+1)N/2^V_j)),0<=i<2^V_j. Its t_j children partition it EXACTLY. Bits outside the interval are never interpreted as additional target bits. For a component invocation at this node, put the interval's actual bits in the PREFIX of a vector of length n_j and append zeros. For its block-source rows, put each selected child's actual bits in a separate prefix-zero-padded vector of length n_j. These are injective linear embeddings of the respective actual block coordinates. They preserve all affine entropies and conditional entropies; no target domain is padded, ignored or restricted.

For u=0,...,H'-1 set Tconf=ceil(ell_log/2000000) and

 p_u^(1)=2^(4000002*(H'-u)*Tconf),
 p_u^(2)=2^(4000002*(H'-u)*Tconf+2000001*Tconf).

Require all p2<=Pmax, otherwise return0. The exact ratios p2/(ell*p1) and earlier p1/(ell*later p2) equal G=2^(2000001*Tconf-ell_log)>=2^Tconf. At a FIXED node all confidences use the SAME map E_j and SAME rectangle width Pmax. Different levels have different component ambient lengths and may use different E_j. The common response seed-count upper bound R=max_{j<H}2^d_j will suffice throughout the proof.

Compute challenge matrices bottom-up. Leaf matrices are empty. At an internal level-j node, test its children in order at confidence Pmax using E_j on the prefix-padded node input. Select the first PASS child, or the first child if none passes. Concatenate the selected child's challenge with the Pmax-bit rows BE_j on EVERY lexicographically ordered b_j-tuple of children, each child padded to n_j. The total number of rows is at most sum_{r=j}^{H-1}binom(t_r,b_r)<=ell and is deterministic at each depth. All tests pad the challenge matrix with zero rows to ell. Follow the resulting selected path. At active depth a_u compare the selected descendant at depth a_(u+1) using the two specified confidence tests with map E_(a_u). At the first first-confidence FAIL output the second-confidence result, encoded FAIL0/PASS1. If there is none output0. This is exactly family.evaluate.

All planned ambient lengths are <=N. The component input embeddings are linear and injective on their actual blocks. The tree has at most N nonempty nodes per depth and J<=L, hence O(N log N) nodes, a sufficient polynomial bound even with unary levels. The tuple-count guard bounds each enumeration by ell<=N eventually; otherwise ell itself is at most a fixed polynomial in N by its formula. Each response has at most N^9 seeds by its enforced full-seed bound, and all component evaluators have a fixed polynomial time exponent at lengths <=N. Planning and guards require at most L levels, only polynomial-size integers, and polynomial-time component constructors. Thus all successful evaluations and all fallback evaluations take time polynomial in actual N, with a fixed program and exponent independent of epsilon. The complete refined full-B2 transfer applies to dispersion at dstar=2^ceil(196Z/25), giving4N-8dstar+4 and every stated same-family4N-N^epsilon bound.

# Challenge-response test and its exact properties

The following test follows Shaltiel's affine challenge-response method. A matrix distribution is eta-somewhere-random if it is within total variation eta of a convex combination of distributions having at least one uniform row; different mixture components may have different good rows. Joint affinity is not required in this definition. Resiliency R of a challenge relative to an affine source X means the same eta-somewhere-random conclusion holds after fixing any possible value of any binary linear map of X of rank at most R.

At a fixed node of depth j, take its injectively prefix-padded input at ambient length n_j and threshold Kmin, within the uniform component window proved below. The supplied constant-error linear extractor has m=floor(Kmin^(16/25)) output bits, d_j=O(log n_j) seed bits and error1/4. Trim to ell by Pmax to obtain E_j. All confidences at this fixed node use the SAME E_j and SAME Pmax. Different nodes at different depths may use different maps. Set R=max_{j<H}2^d_j<=N^9, a common upper bound on every response count. The following four properties apply separately at each node; only the common bound R is used when combining them. Injective linear embeddings preserve every affine source and conditional entropy to which the properties are applied.

For an input x, coordinate block s, challenge function ch on that block with at most ell rows and Pmax columns, and1<=p<=Pmax, pad ch(x_s) with zero rows to ell rows. Return FAIL exactly when some enumerated response E_j(x,y) agrees with that challenge on the first p columns of every row. Otherwise return PASS. Prefix-zero-pad the actual node input to its specified n_j before the response map. Each test is deterministic polynomial time; challenges are computed by the containing algorithm.

The following properties hold.

1. If p1<=p2, FAIL at p2 implies FAIL at p1, by restricting any matching response to fewer columns.

2. If X is affine with H(X)>=Kmin and X_s is fixed, there is an affine subsource of EXACT deficiency ell*p on which the test always fails. Strong linear extraction has error1/4, so at least one seed has a completely uniform ell-by-Pmax response on X: a nonuniform affine output would have distance at least1/2. Fix the ell*p specified response coordinates equal to the constant challenge slice. Uniformity makes this a nonempty affine fiber of codimension ell*p, of probability2^(-ell*p). On that fiber one matching response is guaranteed. The argument does not require ch itself to be linear.

3. More generally, if the affine conditional entropy H(X|X_s)>=Kmin, FAIL has probability at least2^(-ell*p). Fix each possible X_s value. Its fiber remains affine and has entropy at least Kmin; property2 applies to that fiber. Average the resulting lower bounds. The choice of responding seed may depend on the fiber, which is allowed because FAIL enumerates all seeds.

4. If ch(X_s) is eta-somewhere-random with resiliency ell*p relative to X, then FAIL has probability at most R*(2^-p+eta). For any fixed response seed y, the sliced response E_j is a linear map of X with rank at most ell*p. Condition on any of its possible values v. Resiliency makes the challenge eta-close to a mixture with a uniform row. Matching the whole slice to v then has probability at most2^-p+eta, since matching that component's good row alone has probability2^-p. Average over v, and union bound over all R seeds.

All source restrictions here are affine and have explicitly charged rank deficiency. Property3 is a rare-event lower bound, not a high-probability classification of arbitrary blocks. Property4 is a conditional somewhere-randomness hypothesis, not a claim that every nonempty block automatically passes. The final disperser proof must establish that hypothesis for its recursively chosen challenges and compare its error with the rare FAIL probabilities. Replacing R by n^2 without checking the actual seed length would be unjustified; only its polynomial-in-n bound is used in the intended asymptotic parameter comparisons.


# Correctness proof for the subpolynomial-entropy tree construction

This proof uses the exact tree construction above, implemented by family.py, and the complete component arguments in cubic.md, extractors.md and cubic.md. The construction and the two-stage source analysis follow Shaltiel2011; the extractor components are the separately supplied replacements.

Write H for tree depth and H' for the number of consecutive active intervals. All statements in this proof concern sufficiently large actual N, with n=N,L,k,Kmin,ell,Pmax and local arities/ambients as defined in the construction. An affine subsource means conditioning a uniform affine source on possible values of binary linear maps. Its deficiency is the rank lost, and its probability inside the original source is2^(-deficiency). Linear images stay affine. Conditioning on a rank-D map decreases any linear-image entropy, and any conditional entropy of linear images, by at most D. To verify the conditional assertion, write H(A|B)=dim image(A,B)-dim image(B); the first term drops by at most D and the second cannot increase. Entropy chain rules are therefore exact dimension identities in this setting.

## Numerical feasibility with independent active spacing

Write S=log2 k=13Z, W=4Z and a=443/2500. All estimates are uniform over all sufficiently large actual integer lengths and all planned levels. Exact rounded formulas give Z=(1+o(1))*L^(3/4)*B^(1/4), base=(1+o(1))*L^(1/4)*B^(3/4), ell_log=(1+o(1))*sqrt(LB), X=(391/125+o(1))*Z. Since W=o(L), Lpad=(1+o(1))*W^(20/19)=o(L), and W^3>>(L+1)^2, the early guards hold eventually.

The local logarithm is EXACTLY L_j=max(L-V_j,Lpad): ceil(log2 ceil(N/2^v))=L-v at every integer v<=floor(log2 N). In particular Lpad<=L_j<=L, and Lc_j=L_j or L_j+1. The guards W^3>=(L+1)^2 and W^20<=Lpad^19 imply

 (Lc_j)^(2/3)<=W<=(Lc_j)^(19/20).

The supplied component theorem therefore applies uniformly at EVERY planned level. The new degree/output parameters alpha=277779/100000 and beta=16/25 satisfy beta+1/alpha=6944464/6944475<1. Their exact full-seed leading coefficient is1+4alpha*beta=633683/78125=8.1111424, below guard162223/20000=8.11115. Because min L_j>=Lpad tends to infinity and Lc_j/L_j<=1+1/Lpad, the full-seed guard holds uniformly. The larger component guard Q,h^d<=(n_j+1)^12 holds eventually by ds=(633683/78125+o(1))*Lc_j<12Lc_j and d*log h=(1+o(1))*Lc_j. The old power8 guard would not suffice and is explicitly replaced in both implementation and complete component proof. This increases only a fixed polynomial time exponent.

With Amax=M=max_j A_j including the terminal arity, M<=a*base*L/Lpad+1=O(L^(35/76)*B^(37/76)), so M=o(Z), M^2=o(base*L) and M<X eventually. Thus X+M<=3.129Z eventually. The first depth crossing each successive multiple of X gives distinct active levels; every consecutive active log-branch product is at most X+M. Also H'<=L/X+1, so

 H'*ell_log<=(125/391+o(1))*Z.

The first crossing has shortest block log at least floor(log2 N)-X-M>>S/1000, so H'>=1 eventually. If another preliminary active depth follows H, maximality gives shortest next leaf <k^(1/1000), so every current leaf is at most2^(X+M)*(k^(1/1000)+1)+1. If none follows before J, V_J-V_H<X and N<2^(V_J+A_J), so every leaf is at most2^(X+M)+1. Both cases give

 maxleaf<=2^((3.129+.013+o(1))*Z)<k^(1/4).

Every leaf also has length at least k^(1/1000). No comparison between that lower cutoff and the active branch gap is needed.

The complete quadratic depth argument is as follows. Put u_j=L-V_j. Every accepted step beginning with u_j>=Lpad has A_j>=a*base*L/u_j-1 and u_(j+1)=u_j-A_j>=0. Its nonnegative quadratic potential decreases by

 u_j^2-u_(j+1)^2=2u_j*A_j-A_j^2>=2a*base*L-2L-M^2.

This is(2a-o(1))*base*L>0 eventually. Thus there are at most L^2/((2a-o(1))*base*L) such steps. After u_j<Lpad, A_j=Aclip=floor(a*base*L/Lpad), so there are at most Lpad/Aclip+1=O(Lpad^2/(base*L))+1=o(L/base) further steps. Therefore

 H<=J<=(1250/443+o(1))*L^(3/4)/B^(3/4)<=(1411/500)*L^(3/4)/B^(3/4)

at every sufficiently large actual N, since1250/443<1411/500. This includes bottom clipping and rounding.

The prescribed b_j obeys b_j<=2+d_j/log2 g, with log2 g=(9/25+o(1))*W=(36/25+o(1))*Z. Consequently

 b_j*A_j<=((162223/20000)*(443/2500)/(36/25)+o(1))*ell_log+2M
              =(.998122069444...+o(1))*ell_log<.999*ell_log.

The exact leading fraction is71864789/72000000<4991/5000. Since log2 J=o(ell_log), sum_j binom(t_j,b_j)<=J*2^(max b_j*A_j)<=ell eventually, and t_j>=b_j. Uniformly max b_j=O(L/W), so log2(4b_j)<=(25001/100000)*B eventually. Combining with H and (1411/500)*(25001/100000)=35276411/50000000<706/1000 yields

 log2 product_{j<H}(4b_j)<=(706/1000)*Z,
 max consecutive active log-branch product<=(3129/1000)*Z.

Since dstar>=2^(7.84Z), the4^H' factor has log=o(Z), and7.84-3.129=4.711>4.710, every active stopping threshold exceeds2^(4.710Z)=k^(471/1300). Ordinary descent loses at most.706Z in product4b_j and one extra log-arity M=o(Z), leaving more than2^(4.003Z)=k^(4003/13000). Both thresholds exceed maximum leaf length. First restrictions cost at most H*tmax*ell*Pmax<=2^(.641Z), and later resilient-conditioning rank is at most(H+1)*ell*Pmax<=2^(.641Z), because log H+M+ell_log=o(Z). As.641<4<4.001<4.003, subtraction leaves at least2^(4.001Z)=k^(4001/13000)>Kmin=2^(4Z), including the later rank restrictions.

The outputs have logarithm(64/25+o(1))*Z, above log2(ell*Pmax)<.641Z. Their block-source error is at most2^(-Kmin^(7999949/50000000))<2^-Pmax, since4*(7999949/50000000)=.63999592>1279/2000=.6395. Group exponent9/25 and slice639999/1000000 sum to999999/1000000<1, paying the complete raw-row rank guard. Helper precision159999/1000000, design exponent2000001/1000000, single seed319998159999/1000000000000 and near seed319999159999/1000000000000 satisfy the strict near-seed versus cubic-output inequality319999159999/1000000000000<639999/2000000. The supplied full component proof establishes all remaining field/design/helper guards uniformly in the entropy window. All local ambients and nested row dimensions tend to infinity.

Finally (2000001/1000000)*(125/391)<3197/5000, and Tconf=(1/2000000+o(1))*ell_log, give uniformly

 log2(ell^3*p2)<=(3197/5000)*Z<log2 Pmax.

All confidence guards hold and every confidence exceeds every fixed polynomial in L. Now R=max_j2^d_j<=N^9 and H*log2(2R+2)=O(L^2)=o(Pmax). The rare-output/path/challenge estimates below therefore retain their strict margins, with these changed constants. Every numerical guard succeeds for every sufficiently large actual N.


## Locate nested useful nodes on an affine subsource

Fix any affine source X on n bits with H(X)>=dstar. Say X is nice for an interval s at threshold K if all coordinates before s are fixed and H(X_s)>=K. Initially s is the root, so this holds with K=dstar.

The following elementary split is useful. Partition a nice block s into T consecutive children. Let k_i be the entropy of child i conditioned on all earlier children. Their sum is at least K. Choose the first i with k_i>=K/(4T), and fix all earlier children to any possible joint value. Their entropy was below K/4, so the new source loses at most K/4 dimensions. It is nice for child i, whose entropy is at least K/(4T). If that entropy is at least K/4, continue with threshold K/4. If it is below K/4, the remaining conditional entropy inside s outside that child is at least3K/4-K/4=K/2. Thus both the child's entropy and H(X_s|X_child) are at least K/(4T). All these statements concern possible affine fibers, so no probabilistic source selection is assumed.

Apply this split on consecutive ACTIVE levels, with T=2^(V_(a_(u+1))-V_(a_u))<=2^(3129Z/1000). Each continuing step reduces its threshold by at most4. If it continued to the last active level, the leaf entropy would be at least dstar/4^H', greater than k^(471/1300), while its length is below k^(1/4). This is impossible. At a stopping step obtain consecutive active nodes g and g', with g' below g, and an affine subsource on which coordinates before g' are fixed and

 H(X_g')>=k^(471/1300),   H(X_g|X_g')>=k^(471/1300).

Indeed the stopping threshold is at least dstar/(4^H'*4*2^(3129Z/1000)), which is at least k^(196/325-3129/(1000*C)-o(1))>k^(471/1300) by the explicit numerical estimates.

Now descend from g' one ordinary level at a time. Suppose the current nice interval s has entropy at least K. At depth j write b=b_j and t=t_j. If at least b child conditional entropies are >=K/(4t), those b children, in order, are a hidden b-block affine block-wise source at that threshold: conditioning on fewer predecessor blocks cannot decrease conditional entropy. Otherwise, fewer than b indices contribute high conditional entropy; all low indices together contribute at most K/4. Some high index therefore has conditional entropy at least K/(4b). Fix its preceding children and continue in it with threshold K/(4b). These new constraints are inside g', because coordinates before g' were already fixed. They preserve H(X_g|X_g') exactly: conditioning on functions of X_g' selects its fibers but does not change their common affine dimension.

This descent cannot reach a leaf. After at most H steps its retained threshold is at least k^(471/1300)/product_{j<H}(4b_j)>=k^(471/1300-706/(1000*C))>k^(4003/13000), larger than leaf length. Thus at an INTERNAL node e below g', a hidden b-block block-wise source is found. Its threshold is at least k^(471/1300)/(product_{j<H}(4b_j)*4tmax)>=k^(471/1300-706/(1000*C)-o(1))>k^(4003/13000). We have an affine subsource X0 with

 e subset g' subset g;
 g,g' on consecutive active levels;
 H((X0)_g|(X0)_g')>=k^(4003/13000);
 (X0)_e hidden b_(depth(e))-block affine block-wise at threshold k^(4003/13000);
 all coordinates before e fixed.

## Force empty left siblings to fail

For every ancestor s of e and every child z of s lying before e, X0_z is fixed. There are at most H*tmax such pairs. Iteratively apply the exact-deficiency empty-block property of challenge-response to node input X_s and challenge ch_z, at confidence Pmax. It fixes ell*Pmax linear response coordinates and makes that test FAIL identically. An affine fiber in the node projection pulls back to an affine subsource of the full source with the same deficiency. At each stage H(X_s)>=H(X_e)>=k^(4003/13000)-H*tmax*ell*Pmax>Kmin, so that property applies.

The total deficiency is at most H*tmax*ell*Pmax<=k^(641/(1000C)). Since641/(1000C)<4001/13000<4003/13000, subtracting this rank from the retained entropy leaves at least k^(4001/13000) eventually. Previously forced FAIL tests persist pointwise under later restrictions, as all challenge and response functions are fixed deterministic maps. Let X* be the final affine subsource. It still has the above nested nodes, with both useful entropy thresholds at least k^(4001/13000). All left-sibling tests at every ancestor of e FAIL identically. Coordinates before e remain fixed.

## Propagate resilient challenges upward

Put R=max_{j<H}2^d_j<=N^9. For an ancestor s of e at distance a from e define

 eta_a=(2R+2)^a *2^-Pmax,
 D_a=(H-a+1)*ell*Pmax.

We claim ch_s(X*_s) is eta_a-somewhere-random with resiliency D_a relative to the FULL source X*. Resiliency is the rank-conditioning definition in the challenge-response section of proof.md.

For a=0, conditioning on any rank at most(H+1)*ell*Pmax still leaves the hidden block-wise threshold at e at least k^(4001/13000)-(H+1)*ell*Pmax>=Kmin. Node e at depth j lists BE_j on every b_j-tuple of children, using injective prefix padding to n_j. On the correct tuple one listed row is within2^-Pmax of uniform. Appending arbitrary correlated rows preserves that approximation: couple the row to its ideal uniform distribution and extend the other rows using their conditional distribution. This proves the base case for the whole challenge, whether or not its recursively selected child was useful.

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

Conditional entropy H(X*_g|X*_g')>=k^(4001/13000)>=Kmin gives Pr[T_a=FAIL]>=2^(-ell*p_a). The resilient PASS bound gives Pr[T2=FAIL]<=R*(2^-p2+eta_H). Monotonicity makes T2=FAIL imply T1=FAIL. Thus

 Pr[T1=FAIL,T2=FAIL]>=2^-D,
 Pr[T1=FAIL,T2=PASS]>=2^(-ell*p1)-R*(2^-p2+eta_H)
                       >=(1/2)*2^(-ell*p1)>=2^-D

eventually. The gap G diverges, log R=O(log n) is negligible compared with the smallest p2, and Pmax dominates every confidence exponent. The smallest p2 is at least2^(6000003*Tconf), whose logarithm is(6000003/2000000+o(1))*sqrt(LB), so that p2 grows faster than every fixed polynomial in L.

For any strict active ancestor of g, the first-confidence logarithm is at least4000002*(H'-u+1)*Tconf. Subtracting log2 D=Eell+4000002*(H'-u)*Tconf+2000001*Tconf leaves2000001*Tconf-Eell=log2 G. Its first confidence is therefore at least G*D. Since G<=ell^2 eventually (indeed log2 G=(1/2000000+o(1))*Eell), we have G*D<=ell^3*p2. The numerical bounds log2(ell^3*p2)<=(3197/5000)*Z and log2 Pmax>=(1279/2000)*Z show Pmax/(G*D) tends to infinity uniformly over all active levels.

Every strict-ancestor test has the specified resiliency. Union bounding FAIL gives

 H'*R*(2^(-G*D)+eta_H)<=2^(-sqrt(G)*D)

eventually: eta_H<=2^(-Pmax/2), the strict gap1279/2000-3197/5000 pays for that term, and log(H'R)=O(log n) is negligible compared with G*D. If g is the root there are no ancestors and that probability is zero.

The probability of a deviation from the true path to e is at most2^(-Pmax/3). Its sum with2^(-sqrt(G)*D) is below2^-D eventually. Each output event above therefore has a point whose selected path reaches g' and whose earlier active nodes all PASS their first test. On these respective points the evaluator chooses g as the first active FAIL and outputs FAIL=0 and PASS=1. Thus both bits occur on X*, a subsource of the arbitrary affine source X of entropy at least dstar, proving dispersion for the same fixed family at every sufficiently large integer length. No uncharged node choice, extra output or input promise is used.

## Actual-length full-B2 consequence

The tree already acts on all N input bits, without outer padding. Its certified dispersion dimension is dstar=2^ceil(196Z/25)<=N eventually. The evaluator still uses construction k=2^(13Z) and component Kmin=2^(4Z). The complete refined full-B2 theorem in elimination.md therefore gives

 C_B2(F_N)>=4N-8dstar+4.

Since log2 dstar=(196/25+o(1))*L^(3/4)*(log2 L)^(1/4), for every fixed epsilon>0 we have8dstar-4<=N^epsilon at every sufficiently large integer N. Thus the SAME fixed deterministic uniform P family satisfies C_B2(F_N)>=4N-N^epsilon eventually, separately for every fixed epsilon>0. No exponent or program is selected in response to epsilon. Empty-witness NP verification and actual-length accounting hold for this same family. No explicit threshold for the full dispersion/circuit theorem is asserted.
