# One fixed total uniform family with quarter-log loss

For the evaluator F_N=family.evaluate on exactly N input bits, set L=ceil(log2 N), B=ceil(log2(L+1)), U=ceil(L^(3/4)), C=20, Z=U*ceil(B^(1/4)), E=C*Z, and k=2^(2ceil(E/2)) for N>=2. There exists N0 such that for EVERY integer N>=N0,

 C_B2(F_N)>=4N-8k+4.

Here k=2^((20+o(1))*L^(3/4)*(log2 L)^(1/4))=N^o(1). For this SAME fixed family, for EVERY fixed real epsilon>0 there is N_epsilon such that C_B2(F_N)>=4N-N^epsilon for EVERY integer N>=N_epsilon. The program and its time exponent do not depend on epsilon.

B2 includes all16 binary Boolean truth tables, arbitrary sharing and unrestricted depth. Every actual input bit is charged. All fields, parameters and choices are deterministic, without advice, field-polynomial input, oracle, promise or certificate. The family is total uniform P, hence NP with the verifier accepting precisely empty witnesses when its evaluator returns1.

The full-B2 coefficient4 is the existing Mishra--Santhanam theorem, completely supplied in elimination.md. The construction follows Shaltiel2011. Its new constant-error component follows the Shaltiel--Umans2005 interleaved-curve extractor reconstruction, explicitly proved here with the seed available to the predictor so as to establish strong extraction. Trevisan's design reconstruction supplies the separate small-error helper; Rao's reduction supplies the affine row condenser. The cubic ingredient follows Ben-Sasson--Kopparty and Haramaty--Shpilka. Prime-degree field construction follows Shoup; the retained analytic prime-interval proof follows Newman/Zagier. These are credited existing methods. This package does not claim discovery of the coefficient4 theorem or a stronger literature record. All nonstandard premises used below are supplied, not imported from URLs or prior verdicts.

The immediate preceding22 submission96f9752d-967c-4740-8e4f-b3c015e8b667 supplies the sharper prime-degree design estimate. Earlier30 submission4403dc7f-49e8-4060-aab3-5309781859f5 supplies the enlarged smaller-output helper and raw-reservoir balance. Preceding46 submissiond46242dd-8e96-4c95-b7fa-cec67e494e08,68 submission42c2d818-9c1e-4a22-bb28-9e1497d77801 and156 submission875e388c-de8c-4349-b538-562dd42aabfe supply the augmented row, fanout and rational SU/source frameworks. Corrected170 submissionf52757ed-faf4-4adf-873e-8fabef827d1c and original170 submissionf9ac2386-70ff-4cf7-b67b-e846204ec5fb supply the quantized grid. Corrected224af68190e-22b9-4b70-a87b-397dfe692edf is a consulted corrected reference sharing the original224 framework; the component-local correction was applied commonly after original224 rejection, with no new origin or causal route through170 inferred from AF. Original224 submission0d6b9757-9407-42b1-9677-feb0c20af62e,296 submission5abcbc66-bd77-4f01-9b66-9c0f8b22105d,336 submission07e114e5-fa5e-4b96-a079-cbce7e77612e,2048 submission6ab40f99-b8a7-41d9-aadb-2aaa864a692a,3072 replacement764612b3-6862-4e50-b958-89632a4fb0e3 and5120 submission9ff6c2f2-36e4-481b-ae6f-15da0135f777 supply the shared confidence/SU/source framework. Original156f49104cd-7581-4704-bff5-e264dc45309f and68fa6537a5-253f-476c-bc28-32c1288e810f were cancelled for the exact-definition repair. Failed/cancelled or provisional submissions are not mathematical premises. This candidate jointly sharpens the grid, smaller-output exponent and relative-prime slack, and uses degree ceil(12971Lc/(10000log K)), SU output K^(229/1000), grouping K^(771/1000), slice K^(2289/10000), precision K^(1429/25000), fanout703/200 and fully charged confidence/source budgets. The quarter-log constant drops from22 to20. The earlier own5 submission4138000d-8f17-4740-9995-59ee9b9fc772 supplies the large-row condenser. The earlier own36 submission5c80725e-ac7d-47e8-a736-1ce338fd7edf supplies the retuned SU extractor. The earlier own120 submission2c2108c8-4bc1-4b20-8ded-58c2a4295f97 supplies the full-output row reduction. The earlier own168 submissionbb65320f-4cec-4526-8895-6a857ebaef26 supplies the finely charged source entropy analysis. The earlier own submission59e21d17-0da7-44e2-a4f6-ade8d1e69f92 supplies the fully charged active/hidden/final entropy budgets. The earlier own submission0aaaa563-b3e6-41fc-a8c6-8846ba5132c9 supplies the stronger active entropy analysis. The earlier own submission5a30b469-95fb-4a85-aeff-84aa3617a384 supplies the sharper PNT seed estimate. The earlier own submission4cd2c215-066c-41a4-814f-14273960727a supplies the explicit C=8192 slower-fanout analysis. The preceding submitted4c2e1fa2-f0cd-481e-929e-72536cf8c8b3 supplies the explicit square-root-log entropy analysis. The earlier79262846-138d-4327-8058-031bf03f5d8c supplies the generalized extractor parameter lemma and balanced tree proof. The earlier f10654cb-2d81-4a33-b929-707695ca4beb supplies the complete supplied seed-aware SU extractor and subpolynomial-entropy tree framework. The earlier own e6a47252-1b22-49cd-bcc1-1c5d3ba51ea1 and2547b8e7-92aa-466c-a684-3a91b8e64319 supply the tree/component framework, exact unequal-interval accounting, and generic4n-8d+4 transfer. This version uses useful ACTIVE-node entropy k^(4/5), hidden block-source entropy k^(3/4), and final retained entropy k^(2/3)>Kmin=k^(1/2), then balances the tree degree, entropy and response confidence parameters, reducing the LOSS SCALE from5*L^(3/4)*sqrt(log2 L) to20*L^(3/4)*(log2 L)^(1/4), by decoupling active spacing from challenge-row count and supplying the stronger precision/error budgets for ONE fixed family. The seed-aware extractor and all component arguments are supplied again. The earlier2d32f686-6419-451e-8758-22dd489737e0 is a comparator and research motivation: its different algebraic selector has an explicitly positive-power loss, whereas this new family has the supplied subpolynomial-loss theorem. No earlier result is admitted as a proof premise.

cubic.md supplies the complete field constructor, analytic prime-interval argument, cubic isolation/bias and vector extractor. extractors.md supplies local codes, both strong seeded-extractor regimes, the new Singer-cycle construction, affine row condensation and block-source extraction. This file supplies the complete tree and source analysis. family.py implements exact integer roots and total numerical fallbacks. The small component audit does not execute the asymptotically nontrivial full tree or establish any universal theorem.

The tree logarithm is L=ceil(log2 n). Define the DISTINCT component-local logarithm Lc=ceil(log2(n+1)). In the SU parameter construction and proof in extractors.md, its locally denoted L means Lc, so the exact SU degree uses Lc as implemented by optimal_parameters. At powers of two, Lc=L+1; otherwise they agree for n>=2. No exact equality between these logarithms or their rounded degrees is claimed. The global full-seed guard uses tree L, with a strict eventual margin established by the supplied comparison Lc=L or L+1. All nonstandard proofs and both original joint claims use this explicit convention.

# A concrete subpolynomial-entropy Shaltiel tree disperser

The following fixed construction is proved correct below, using the complete supplied component arguments. The approach is Shaltiel2011's challenge-response disperser, with global threshold k=2^(20*ceil(L^(3/4))*ceil((ceil(log2(L+1)))^(1/4))) and explicitly supplied extractor components instead of external black boxes.

At actual input length N>=2 set n=N, L=ceil(log2 n), B=ceil(log2(L+1)), U=ceil(L^(3/4)), Z=U*ceil(B^(1/4)), C=20, A=ceil((703/200)*ceil(L^(1/4))*ceil(B^(3/4))), t=2^A and J=max{j:t^j<=n}. All powers and roots use exact integer comparisons. No target input is padded or ignored. Define

 E=C*Z, Kmin=2^ceil(E/2), k=Kmin^2=2^E,
 q=ceil(sqrt(ceil(L/B))),
 ell_log=ceil(sqrt L)*ceil(sqrt B), ell=2^ell_log,
 Pmax=2^ceil(57Z/100).

The active spacing q is deliberately INDEPENDENT of ell_log. C is even, so E is even and k=2^E exactly. The early guard k<=n is evaluated as an exponent comparison before allocating Kmin or Pmax.

Build the constant-error linear response map at length n and entropy Kmin, trim its output to ell*Pmax bits, and view it as ell rows of Pmax columns. Let Rresp=2^d be its number of seeds. Build the affine block-source extractor at block length n, entropy Kmin, block count b specified by its row-reduction construction, and trim its output to Pmax bits.

Before constructing Kmin, require2ceil(E/2)<=floor(log2 n), otherwise return0; this prevents allocating oversized integers for infeasible small inputs and holds eventually. Use numerical feasibility checks for these constructions (the integer guards in the supplied subpolynomial-entropy extractor, enough output bits, group sizes at least2, and locally constructed field/graph/design syntax). Whenever a required numerical construction is undefined or has too few bits, return0. The following proof must show these checks succeed eventually; extraction of a particular input is not tested by a guard. Assign output0 also at N=0,1.

Consider the ordered t-ary interval tree on[n], using zero-based coordinate positions. At level j, node i is the half-open interval[floor(i*n/t^j),floor((i+1)*n/t^j)), for0<=i<t^j. Its t consecutive children partition it exactly. Lengths at a fixed level differ by at most1. Let H be the largest multiple of q in[0,J] for which floor(n/t^H)^4>=k; thus even the shortest leaf has length at least k^(1/4). Truncate the tree at depth H. If H<q, return0. Active levels are0,q,2q,...,H; put H'=H/q. Require t>=b and H*binom(t,b)<=ell, otherwise return0. For u=0,...,H'-1 define Eell=log2 ell, Tconf=ceil(Eell/20000), and the two exact confidence powers

 p_u^(1)=2^(40002*(H'-u)*Tconf),
 p_u^(2)=2^(40002*(H'-u)*Tconf+20001*Tconf).

Put G=2^(20001*Tconf-Eell). Since Eell<=20000*Tconf, G>=2^Tconf tends to infinity. The ratio p2/(ell*p1) is EXACTLY G. For consecutive active levels the earlier first confidence, divided by ell times the later second confidence, is also EXACTLY G. These two divergent gaps replace the more widely spaced earlier confidence grid.

Require all confidences <=Pmax, otherwise return0. Every test uses the same global response map and Pmax rectangle. A node's response input is its actual bits in their original[n] positions, with outside coordinates zero; this is an analysis subblock projection, not an additional target input. Its challenge function depends only on its own interval input.

Compute challenge matrices bottom-up. A leaf's matrix is empty. At an internal node s with children s1,...,st, test every child with the challenge-response procedure at confidence Pmax and node input x_s. Choose the first child for which the test returns PASS; if none passes, choose s1. Let i(s) be this deterministic choice. The challenge ch_s is the concatenation of ch_(s_i(s)) and, in lexicographic tuple order, all the rows

 BE(x_(s_i1),...,x_(s_ib)),  i1<...<ib in[t],

with each child block padded with zeros to length n before BE. At most(H-level(s))*binom(t,b) rows occur, so the global ell bound holds. The empty leaf matrix causes no ambiguity: challenge-response pads every matrix with zero rows to ell, as specified in its lemma.

The selected path follows i(s) from the root to a leaf. For each nonleaf active node s at level uq on that path, let s' be its selected descendant at level(u+1)q. Compute o_s^(a)=CR_(ch_s',p_u^(a))(x_s,s') for a=1,2. Choose the first active node on the path, starting at the root, whose o_s^(1) is FAIL. Output its o_s^(2), encoded FAIL=0 and PASS=1. If there is no such node, output0.

The tree has at most O(n) nodes. Challenge matrices, tuple lists and response tables have polynomial size for successful parameters; all component evaluators are polynomial time. The supplied optimal constant-error extractor has response seed length at most(10941/5000)ceil(log2 n) eventually, by its sharpened prime-degree estimate. The evaluator requires this numerical bound; otherwise it returns0. The block count b=O(L^(1/4)/(log L)^(1/4)); the explicit tuple-count guard bounds every enumeration by ell. Both seed enumeration and all transfers therefore have fixed polynomial running time. Here n=N exactly. With the stated total fallbacks, this is a fixed uniform P family on every actual input length. It does not choose a different program or time exponent for different requested lower bounds.

The correctness argument proves dispersion at dimension k on the actual N input bits. The supplied refined full-B2 theorem gives C(F_N)>=4N-8k+4. Since log2 k=(C+o(1))*L^(3/4)*(log2 L)^(1/4)=o(log N), the loss is eventually at most N^epsilon for each fixed epsilon>0. This is a statement about the same program for every epsilon, with no explicit full threshold claimed.

# Challenge-response test and its exact properties

The following test follows Shaltiel's affine challenge-response method. A matrix distribution is eta-somewhere-random if it is within total variation eta of a convex combination of distributions having at least one uniform row; different mixture components may have different good rows. Joint affinity is not required in this definition. Resiliency R of a challenge relative to an affine source X means the same eta-somewhere-random conclusion holds after fixing any possible value of any binary linear map of X of rank at most R.

Fix a global input length n, entropy threshold log Kmin=(C/2+o(1))*L^(3/4)*(log2 L)^(1/4), row bound ell and column count Pmax. The constant-error linear extractor with threshold Kmin supplies m=floor(Kmin^(229/1000)) bits and d=O(log n) seed bits. Assume ell*Pmax<=m and trim its output to that rectangle. Denote this fixed linear response map by E(x,y), for y in{0,1}^d. Let R=2^d, a polynomial in n. Use the SAME E and SAME Pmax for every confidence parameter below.

For an input x, coordinate block s, challenge function ch on that block with at most ell rows and Pmax columns, and1<=p<=Pmax, pad ch(x_s) with zero rows to ell rows. Return FAIL exactly when some enumerated response E(x,y) agrees with that challenge on the first p columns of every row. Otherwise return PASS. If an input block is shorter than n, pad it with zeros in its specified ambient positions before the response map. Each test is deterministic polynomial time; challenges are computed by the containing algorithm.

The following properties hold.

1. If p1<=p2, FAIL at p2 implies FAIL at p1, by restricting any matching response to fewer columns.

2. If X is affine with H(X)>=Kmin and X_s is fixed, there is an affine subsource of EXACT deficiency ell*p on which the test always fails. Strong linear extraction has error1/4, so at least one seed has a completely uniform ell-by-Pmax response on X: a nonuniform affine output would have distance at least1/2. Fix the ell*p specified response coordinates equal to the constant challenge slice. Uniformity makes this a nonempty affine fiber of codimension ell*p, of probability2^(-ell*p). On that fiber one matching response is guaranteed. The argument does not require ch itself to be linear.

3. More generally, if the affine conditional entropy H(X|X_s)>=Kmin, FAIL has probability at least2^(-ell*p). Fix each possible X_s value. Its fiber remains affine and has entropy at least Kmin; property2 applies to that fiber. Average the resulting lower bounds. The choice of responding seed may depend on the fiber, which is allowed because FAIL enumerates all seeds.

4. If ch(X_s) is eta-somewhere-random with resiliency ell*p relative to X, then FAIL has probability at most R*(2^-p+eta). For any fixed response seed y, the sliced response is a linear map of X with rank at most ell*p. Condition on any of its possible values v. Resiliency makes the challenge eta-close to a mixture with a uniform row. Matching the whole slice to v then has probability at most2^-p+eta, since matching that component's good row alone has probability2^-p. Average over v, and union bound over all R seeds.

All source restrictions here are affine and have explicitly charged rank deficiency. Property3 is a rare-event lower bound, not a high-probability classification of arbitrary blocks. Property4 is a conditional somewhere-randomness hypothesis, not a claim that every nonempty block automatically passes. The final disperser proof must establish that hypothesis for its recursively chosen challenges and compare its error with the rare FAIL probabilities. Replacing R by n^2 without checking the actual seed length would be unjustified; only its polynomial-in-n bound is used in the intended asymptotic parameter comparisons.


# Correctness proof for the subpolynomial-entropy tree construction

This proof uses the exact tree construction above, implemented by family.py, and the complete component arguments in cubic.md, extractors.md and cubic.md. The construction and the two-stage source analysis follow Shaltiel2011; the extractor components are the separately supplied replacements.

Write H for tree depth and H'=H/q for the number of active intervals. All statements in this proof concern sufficiently large actual N, with n,t,L,k,Kmin,ell,Pmax as defined in the construction. An affine subsource means conditioning a uniform affine source on possible values of binary linear maps. Its deficiency is the rank lost, and its probability inside the original source is2^(-deficiency). Linear images stay affine. Conditioning on a rank-D map decreases any linear-image entropy, and any conditional entropy of linear images, by at most D. To verify the conditional assertion, write H(A|B)=dim image(A,B)-dim image(B); the first term drops by at most D and the second cannot increase. Entropy chain rules are therefore exact dimension identities in this setting.

## Numerical feasibility with independent active spacing

Write S=log2 k=CZ, which is EXACT because C=20 is even. All following estimates hold uniformly at sufficiently large integer N. The exact rounded parameters satisfy

 Z=(1+o(1))*L^(3/4)*B^(1/4),
 A=(703/200+o(1))*L^(1/4)*B^(3/4),
 log2 ell=(1+o(1))*sqrt(L)*sqrt(B),
 q=(1+o(1))*sqrt(L)/sqrt(B).

Here B=ceil(log2(L+1)) tends to infinity. In particular S=o(L), so the early guard k<=N succeeds eventually. The total family and its polynomial time exponent remain fixed. Also at depth q the shortest block has logarithmic length L-O(Z), which exceeds S/4 eventually because Z=o(L), so H>=q eventually.

Since H<=L/A and q>=sqrt(L/B),

 H'*(log2 ell)<=L*(log2 ell)/(A*q)<=(200/703+o(1))*Z.

Also log2 ell=o(Z), Tconf=(1/20000+o(1))*log2 ell, and (20001/10000)*(200/703)=20001/35150<1139/2000. Every confidence p2, multiplied by ell^3, therefore obeys

 log2(ell^3*p2)<=(40002H'+20001)*Tconf+3*log2 ell<=(1139/2000)*Z

eventually. Thus every confidence fits Pmax=2^ceil(57Z/100). This strict POWER gap, rather than an incorrect little-o comparison of their logarithms, pays for all rare-event comparisons below.

We have A*q<=(879/250)*Z eventually. The shortest leaf has length floor(N/t^H)>=k^(1/4). Eventually H+q<=J: otherwise the largest leaf would be below t^(q+1)+1, contradicting k^(1/4)>>t^(q+1), since A(q+1)<=(879/250)*Z+o(Z)=(879/(250*C)+o(1))*S<S/4. Maximality then gives floor(N/t^(H+q))<k^(1/4), so each leaf is shorter than

 t^q*(k^(1/4)+1)+1<=k^(1/4+879/(250*C)+o(1))<k^(1/2).

The tree can have unequal interval sizes; this argument only uses their exact floor/ceiling bounds and ordered partition.

The retuned strong extractor has full seed at most(10941/5000)L eventually. For the precision-budget block source, R=floor(Kmin^(229/1000)) and g=floor(Kmin^(771/1000)); hence log2 g>=S/(12971/5000) eventually. Its exact prescribed block count satisfies

 b<=2+(10941/5000)*L/(S/(12971/5000))<=(5677/1000)*L/S.

Because A=(703/200+o(1))*L^(1/4)*B^(3/4) and log2 ell>=(sqrt L)*(sqrt B), this implies

 b*A<=(5677*703/(200000*C)+o(1))*(log2 ell)<(999/1000)*log2 ell.

With log2 H=o(log2 ell), H*binom(t,b)<=H*t^b<=ell eventually. Also b<t. This proves the challenge-row guard without coupling active spacing q to the logarithm of the row count.

For the ordinary source descent, log2(4b)<=B eventually and H<=L/A<=(200/703)*L^(3/4)/B^(3/4). Therefore

 log2((4b)^H)<=(57/200)*Z=57*S/(200*C),
 (4b)^H<=k^(57/(200*C)), t^q<=k^(879/(250*C)).

The active-level factor4^H' has logarithm O(L^(1/4)/B^(1/4))=o(S). Thus the active stopping threshold is at least k^(1-879/(250*C)-o(1))>k^(4/5), and the ordinary hidden-source threshold remains above k^(4/5-57/(200*C)-o(1))>k^(3/4). All margins are strict fixed exponent comparisons.

Rank deficiencies are NOT asserted k^o(1). Instead log2 ell+A+log2(H+1)=o(Z), and Pmax=2^ceil(57Z/100), so

 H*t*ell*Pmax<=2^((29/50)*Z)=k^(29/(50C)),
 (H+1)*ell*Pmax<=k^(29/(50C))

eventually. Since29/(50C)<1/2<2/3<3/4, the hidden source survives the first restrictions with entropy at least k^(2/3); subsequent resilient-challenge rank conditioning still leaves entropy at least Kmin=sqrt(k). Each claimed inequality follows by subtracting at most the displayed rank, not by an uncharged entropy assumption.

The response and block-source outputs have length floor(Kmin^(229/1000)). Their logarithms are(229*C/2000+o(1))*Z, exceeding log2(ell*Pmax)<=(29/50)*Z because 229*C/2000>29/50. The precision-budget block-source error is at most2^(-Kmin^(1143/20000)). Since

 Kmin^(1143/20000)=2^((1143*C/40000)*Z), Pmax=2^ceil(57Z/100),

this error is below2^-Pmax eventually. The strict inequality 1143*C/40000=1143/2000>57/100 is essential: it replaces a false assertion that every positive power of Kmin automatically dominates Pmax. The uniform optimal-extractor window and all small-error/cubic/group guards hold, as log Kmin=(C/2+o(1))*L^(3/4)*(log2 L)^(1/4) lies between L^(2/3) and L^(19/20) eventually, and all nested row lengths grow as fixed positive powers of Kmin.

Every confidence is much larger than log N, while log2 Rresp=O(L). The bound ell^3*p2<=2^((1139/2000)*Z) and Pmax=2^ceil(57Z/100) is uniform over all active levels. It shows Pmax/3 dominates every D=ell*p2 and every ell^2*D, as required to separate true output events from path/ancestor errors. All numerical parameter guards therefore pass at every sufficiently large integer length.


## Locate nested useful nodes on an affine subsource

Fix any affine source X on n bits with H(X)>=k. Say X is nice for an interval s at threshold K if all coordinates before s are fixed and H(X_s)>=K. Initially s is the root, so this holds with K=k.

The following elementary split is useful. Partition a nice block s into T consecutive children. Let k_i be the entropy of child i conditioned on all earlier children. Their sum is at least K. Choose the first i with k_i>=K/(4T), and fix all earlier children to any possible joint value. Their entropy was below K/4, so the new source loses at most K/4 dimensions. It is nice for child i, whose entropy is at least K/(4T). If that entropy is at least K/4, continue with threshold K/4. If it is below K/4, the remaining conditional entropy inside s outside that child is at least3K/4-K/4=K/2. Thus both the child's entropy and H(X_s|X_child) are at least K/(4T). All these statements concern possible affine fibers, so no probabilistic source selection is assumed.

Apply this split on consecutive ACTIVE levels, with T=t^q. Each continuing step reduces its threshold by at most4. If it continued to the last active level, the leaf entropy would be at least k/4^H', greater than k^(4/5), while its length is below k^(1/2). This is impossible. At a stopping step obtain consecutive active nodes g and g', with g' below g, and an affine subsource on which coordinates before g' are fixed and

 H(X_g')>=k^(4/5),   H(X_g|X_g')>=k^(4/5).

Indeed the stopping threshold is at least k/(4^H'*4t^q), which is at least k^(1-879/(250*C)-o(1))>k^(4/5) by the explicit numerical estimates.

Now descend from g' one ordinary level at a time. Suppose the current nice interval s has entropy at least K. If at least b child conditional entropies are >=K/(4t), those b children, in order, are a hidden b-block affine block-wise source at that threshold: conditioning on fewer predecessor blocks cannot decrease conditional entropy. Otherwise, fewer than b indices contribute high conditional entropy; all low indices together contribute at most K/4. Some high index therefore has conditional entropy at least K/(4b). Fix its preceding children and continue in it with threshold K/(4b). These new constraints are inside g', because coordinates before g' were already fixed. They preserve H(X_g|X_g') exactly: conditioning on functions of X_g' selects its fibers but does not change their common affine dimension.

This descent cannot reach a leaf. After at most H steps its retained threshold is at least k^(4/5)/(4b)^H>=k^(4/5-57/(200*C))>k^(3/4), larger than leaf length. Thus at an INTERNAL node e below g', a hidden b-block block-wise source is found. Its threshold is at least k^(4/5)/((4b)^H*4t)>=k^(4/5-57/(200*C)-o(1))>k^(3/4). We have an affine subsource X0 with

 e subset g' subset g;
 g,g' on consecutive active levels;
 H((X0)_g|(X0)_g')>=k^(3/4);
 (X0)_e hidden b-block affine block-wise at threshold k^(3/4);
 all coordinates before e fixed.

## Force empty left siblings to fail

For every ancestor s of e and every child z of s lying before e, X0_z is fixed. There are at most H*t such pairs. Iteratively apply the exact-deficiency empty-block property of challenge-response to node input X_s and challenge ch_z, at confidence Pmax. It fixes ell*Pmax linear response coordinates and makes that test FAIL identically. An affine fiber in the node projection pulls back to an affine subsource of the full source with the same deficiency. At each stage H(X_s)>=H(X_e)>=k^(3/4)-H*t*ell*Pmax>Kmin, so that property applies.

The total deficiency is at most H*t*ell*Pmax<=k^(29/(50C)). Since29/(50C)<2/3<3/4, subtracting this rank from the retained entropy leaves at least k^(2/3) eventually. Previously forced FAIL tests persist pointwise under later restrictions, as all challenge and response functions are fixed deterministic maps. Let X* be the final affine subsource. It still has the above nested nodes, with both useful entropy thresholds at least k^(2/3). All left-sibling tests at every ancestor of e FAIL identically. Coordinates before e remain fixed.

## Propagate resilient challenges upward

Put R=Rresp. For an ancestor s of e at distance a from e define

 eta_a=(2R+2)^a *2^-Pmax,
 D_a=(H-a+1)*ell*Pmax.

We claim ch_s(X*_s) is eta_a-somewhere-random with resiliency D_a relative to the FULL source X*. Resiliency is the rank-conditioning definition in the challenge-response section of proof.md.

For a=0, conditioning on any rank at most(H+1)*ell*Pmax still leaves the hidden block-wise threshold at e at least k^(2/3)-(H+1)*ell*Pmax>=Kmin. Node e lists BE on every b-tuple of children. On the correct tuple one listed row is within2^-Pmax of uniform. Appending arbitrary correlated rows preserves that approximation: couple the row to its ideal uniform distribution and extend the other rows using their conditional distribution. This proves the base case for the whole challenge, whether or not its recursively selected child was useful.

For the induction step let z be the child of s toward e. Fix any rank-D_a conditioning of X*. The inductive guarantee for ch_z leaves an additional ell*Pmax of resiliency, and its error remains eta_(a-1). Challenge-response therefore tests z with FAIL probability at most

 delta=R*(2^-Pmax+eta_(a-1)).

Every earlier child of s FAILs identically, even after this extra conditioning. Thus, except for probability delta, s selects z and its challenge contains the entire challenge of z in the prescribed leading rows. Replacing the exceptional choices changes the distribution by at most delta. The child challenge's error and arbitrary appended rows add at most eta_(a-1). Their sum is at most(2R+1)*eta_(a-1)<=eta_a. This proves the induction. Variable mixture rows cause no issue: the definition needs only a uniform row in each component, and row count at a given level is fixed.

Uniformly for every a<=H,

 eta_a<=2^(-Pmax/2)

eventually, since H*log2(2R+2)=O(L^(7/4)/(log L)^(3/4))=o(Pmax). In particular every true-path challenge has resiliency at least ell*Pmax. Applying the PASS bound at each edge to e and union bounding shows the selected path reaches e, and hence g', with probability at least1-epsilon_path, where

 epsilon_path<=H*R*(2^-Pmax+eta_H)<=2^(-Pmax/3).

All probabilities are with respect to the single fixed affine subsource X*. Left siblings never interfere; the true child is chosen whenever its test passes.

## Both output values occur with the quantized confidence grid

Let g be at active level uq. Define fixed-node tests T_a=CR_(ch_g',p_u^(a))(X*_g,g') for a=1,2. When the selected path reaches g', these are exactly g's output tests. Write p1,p2 for the two confidences, Eell=log2 ell, Tconf=ceil(Eell/20000), D=ell*p2, and G=2^(20001*Tconf-Eell). The exact integer grid gives p2/(ell*p1)=G>=2^Tconf, tending to infinity uniformly over all active levels.

Conditional entropy H(X*_g|X*_g')>=k^(2/3)>=Kmin gives Pr[T_a=FAIL]>=2^(-ell*p_a). The resilient PASS bound gives Pr[T2=FAIL]<=R*(2^-p2+eta_H). Monotonicity makes T2=FAIL imply T1=FAIL. Thus

 Pr[T1=FAIL,T2=FAIL]>=2^-D,
 Pr[T1=FAIL,T2=PASS]>=2^(-ell*p1)-R*(2^-p2+eta_H)
                       >=(1/2)*2^(-ell*p1)>=2^-D

eventually. The gap G diverges, log R=O(log n) is negligible compared with the smallest p2, and Pmax dominates every confidence exponent. The smallest p2 is at least2^(60003*Tconf), whose logarithm is(60003/20000+o(1))*sqrt(LB), so that p2 grows faster than every fixed polynomial in L.

For any strict active ancestor of g, the first-confidence logarithm is at least40002*(H'-u+1)*Tconf. Subtracting log2 D=Eell+40002*(H'-u)*Tconf+20001*Tconf leaves20001*Tconf-Eell=log2 G. Its first confidence is therefore at least G*D. Since G<=ell^2 eventually (indeed log2 G=(1/20000+o(1))*Eell), we have G*D<=ell^3*p2. The numerical bounds log2(ell^3*p2)<=(1139/2000)*Z and log2 Pmax>=(57/100)*Z show Pmax/(G*D) tends to infinity uniformly over all active levels.

Every strict-ancestor test has the specified resiliency. Union bounding FAIL gives

 H'*R*(2^(-G*D)+eta_H)<=2^(-sqrt(G)*D)

eventually: eta_H<=2^(-Pmax/2), the strict gap57/100-1139/2000 pays for that term, and log(H'R)=O(log n) is negligible compared with G*D. If g is the root there are no ancestors and that probability is zero.

The probability of a deviation from the true path to e is at most2^(-Pmax/3). Its sum with2^(-sqrt(G)*D) is below2^-D eventually. Each output event above therefore has a point whose selected path reaches g' and whose earlier active nodes all PASS their first test. On these respective points the evaluator chooses g as the first active FAIL and outputs FAIL=0 and PASS=1. Thus both bits occur on X*, a subsource of the arbitrary dimension-k source X, proving dispersion for the same fixed family at every sufficiently large integer length. No uncharged node choice, extra output or input promise is used.

## Actual-length full-B2 consequence

The tree already acts on all N input bits, without outer padding. Its global threshold k=2^(20*ceil(L^(3/4))*ceil((ceil(log2(L+1)))^(1/4)))<=N eventually. The complete refined full-B2 theorem in elimination.md therefore gives

 C_B2(F_N)>=4N-8k+4.

Since log2 k=(C+o(1))*L^(3/4)*(log2 L)^(1/4), for every fixed epsilon>0 we have8k-4<=N^epsilon at every sufficiently large integer N. Thus the SAME fixed deterministic uniform P family satisfies C_B2(F_N)>=4N-N^epsilon eventually, separately for every fixed epsilon>0. No exponent or program is selected in response to epsilon. Empty-witness NP verification and actual-length accounting hold for this same family. No explicit threshold for the full dispersion/circuit theorem is asserted.
