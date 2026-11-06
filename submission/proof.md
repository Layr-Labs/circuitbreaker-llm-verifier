# One fixed total uniform family with a smaller error order

For the evaluator F_N=family.evaluate on exactly N input bits, there exists N0 such that for EVERY integer N>=N0,

 C_B2(F_N)>=4N-8N^(1/50)-16N^(1/100)-4.

Here B2 includes all16 binary Boolean truth tables with arbitrary sharing and unrestricted circuit depth; the output is one Boolean bit. Actual input length includes every input bit. All fields and parameters are constructed deterministically, with no advice, field-polynomial inputs, oracle, promise or uncharged certificate. The evaluator is total uniform P; the NP verifier accepts precisely empty witnesses when the evaluator returns1. Empty witnesses exist exactly for members of its language. The family does not vary with a requested bound.

The leading full-B2 coefficient4 is the Mishra--Santhanam affine-disperser theorem, supplied in elimination.md. The construction route is Shaltiel2011, using supplied replacement extractor arguments following Trevisan, Rao, Ben-Sasson--Kopparty and Haramaty--Shpilka. These are credited existing methods, not a claim to discover a new4n theorem. cubic.md supplies the prime-degree constructor (Shoup specialization) and analytic prime-interval proof retained from earlier declared source material. Every nonstandard premise required by the chain is proved in the evidence; URLs and prior verdicts are not mathematical premises.

The direct preceding own submission2547b8e7-92aa-466c-a684-3a91b8e64319 supplies this tree/extractor proof framework. This version fixes the entropy exponent at1/100 and uses unequal interval lengths, removing all outer input padding. The components and correctness argument are supplied again rather than admitted by a prior verdict. The earlier own submission2d32f686-6419-451e-8758-22dd489737e0 is used as a benchmark comparator and research motivation. The present construction changes the family rather than retuning that selector. No numeric certificate, verdict or loss-optimality assertion about the previous family is a premise here; the new family is proved from the supplied evidence.

The proof files divide only for size: cubic.md proves field existence/construction and prime intervals; cubic.md proves coefficient isolation, cubic bias and vector extraction; extractors.md proves the local linear code, reconstruction, seeded extractors, affine row condenser and block-source extraction; the following sections specify the final tree and prove its all-length transfer. family.py implements these maps with exact integer roots and total parameter fallbacks.

# A concrete polynomial-entropy variant of Shaltiel's tree disperser

The following fixed construction is proved correct below, using the complete supplied component arguments. The approach is Shaltiel2011's challenge-response disperser, with global entropy threshold ceil(n^(1/100))^2 and explicitly supplied polynomial-entropy extractor components instead of its external black boxes.

At actual input length N>=2 set n=N, L=ceil(log2 n), A=ceil(L^(3/10)), t=2^A, and J=max{j:t^j<=n}. No target input is padded or ignored. All roots/ceilings are exact integer comparisons: e.g. A is the least integer with A^10>=L^3. Define

 Kmin=ceil(n^(1/100)), k=Kmin^2,
 q=ceil(sqrt L), ell=2^q,
 Pmax=2^ceil(L^(4/5)).

Build the constant-error linear response map at length n and entropy Kmin, trim its output to ell*Pmax bits, and view it as ell rows of Pmax columns. Let Rresp=2^d be its number of seeds. Build the affine block-source extractor at block length n, entropy Kmin, block count b specified by its row-reduction construction, and trim its output to Pmax bits.

Use numerical feasibility checks for these constructions (the required source parameter Kmin>=n^(1/100), enough output bits, group sizes at least2, and locally constructed field/graph/design syntax). Whenever a required numerical construction is undefined or has too few bits, return0. The following proof must show these checks succeed eventually; extraction of a particular input is not tested by a guard. Assign output0 also at N=0,1.

Consider the ordered t-ary interval tree on[n], using zero-based coordinate positions. At level j, node i is the half-open interval[floor(i*n/t^j),floor((i+1)*n/t^j)), for0<=i<t^j. Its t consecutive children partition it exactly. Lengths at a fixed level differ by at most1. Let H be the largest multiple of q in[0,J] for which floor(n/t^H)^4>=k; thus even the shortest leaf has length at least k^(1/4). Truncate the tree at depth H. If H<q, return0. Active levels are0,q,2q,...,H; put H'=H/q. Require t>=b and H*binom(t,b)<=ell, otherwise return0. For u=0,...,H'-1 define

 p_u^(1)=ell^(6(H'-u)+3),
 p_u^(2)=ell^(6(H'-u)+6).

Require all confidences <=Pmax, otherwise return0. Every test uses the same global response map and Pmax rectangle. A node's response input is its actual bits in their original[n] positions, with outside coordinates zero; this is an analysis subblock projection, not an additional target input. Its challenge function depends only on its own interval input.

Compute challenge matrices bottom-up. A leaf's matrix is empty. At an internal node s with children s1,...,st, test every child with the challenge-response procedure at confidence Pmax and node input x_s. Choose the first child for which the test returns PASS; if none passes, choose s1. Let i(s) be this deterministic choice. The challenge ch_s is the concatenation of ch_(s_i(s)) and, in lexicographic tuple order, all the rows

 BE(x_(s_i1),...,x_(s_ib)),  i1<...<ib in[t],

with each child block padded with zeros to length n before BE. At most(H-level(s))*binom(t,b) rows occur, so the global ell bound holds. The empty leaf matrix causes no ambiguity: challenge-response pads every matrix with zero rows to ell, as specified in its lemma.

The selected path follows i(s) from the root to a leaf. For each nonleaf active node s at level uq on that path, let s' be its selected descendant at level(u+1)q. Compute o_s^(a)=CR_(ch_s',p_u^(a))(x_s,s') for a=1,2. Choose the first active node on the path, starting at the root, whose o_s^(1) is FAIL. Output its o_s^(2), encoded FAIL=0 and PASS=1. If there is no such node, output0.

The tree has at most O(n) nodes. Challenge matrices, tuple lists and response tables have polynomial size for successful parameters; all component evaluators are polynomial time. The response seed length is O(log n) because Kmin is polynomial in n. The block count b is bounded independently of n. Here n=N exactly. With the stated total fallbacks, this is a fixed uniform P family on every actual input length. It does not choose a different program or time exponent for different requested lower bounds.

The correctness argument below proves dispersion at dimension k=ceil(N^(1/100))^2 on the actual N input bits. The supplied refined full-B2 theorem gives C(F_N)>=4N-8k+4. With z=N^(1/100) and k=ceil(z)^2<=(z+1)^2, this is at least4N-8N^(1/50)-16N^(1/100)-4, and hence at least4N-9N^(1/50) eventually. No explicit full threshold is claimed.

# Challenge-response test and its exact properties

The following test follows Shaltiel's affine challenge-response method. A matrix distribution is eta-somewhere-random if it is within total variation eta of a convex combination of distributions having at least one uniform row; different mixture components may have different good rows. Joint affinity is not required in this definition. Resiliency R of a challenge relative to an affine source X means the same eta-somewhere-random conclusion holds after fixing any possible value of any binary linear map of X of rank at most R.

Fix a global input length n, entropy threshold Kmin>=n^(1/100), row bound ell and column count Pmax. The constant-error linear extractor with threshold Kmin supplies m=floor(Kmin^(1/8)) bits and d=O(log n) seed bits. Assume ell*Pmax<=m and trim its output to that rectangle. Denote this fixed linear response map by E(x,y), for y in{0,1}^d. Let R=2^d, a polynomial in n. Use the SAME E and SAME Pmax for every confidence parameter below.

For an input x, coordinate block s, challenge function ch on that block with at most ell rows and Pmax columns, and1<=p<=Pmax, pad ch(x_s) with zero rows to ell rows. Return FAIL exactly when some enumerated response E(x,y) agrees with that challenge on the first p columns of every row. Otherwise return PASS. If an input block is shorter than n, pad it with zeros in its specified ambient positions before the response map. Each test is deterministic polynomial time; challenges are computed by the containing algorithm.

The following properties hold.

1. If p1<=p2, FAIL at p2 implies FAIL at p1, by restricting any matching response to fewer columns.

2. If X is affine with H(X)>=Kmin and X_s is fixed, there is an affine subsource of EXACT deficiency ell*p on which the test always fails. Strong linear extraction has error1/4, so at least one seed has a completely uniform ell-by-Pmax response on X: a nonuniform affine output would have distance at least1/2. Fix the ell*p specified response coordinates equal to the constant challenge slice. Uniformity makes this a nonempty affine fiber of codimension ell*p, of probability2^(-ell*p). On that fiber one matching response is guaranteed. The argument does not require ch itself to be linear.

3. More generally, if the affine conditional entropy H(X|X_s)>=Kmin, FAIL has probability at least2^(-ell*p). Fix each possible X_s value. Its fiber remains affine and has entropy at least Kmin; property2 applies to that fiber. Average the resulting lower bounds. The choice of responding seed may depend on the fiber, which is allowed because FAIL enumerates all seeds.

4. If ch(X_s) is eta-somewhere-random with resiliency ell*p relative to X, then FAIL has probability at most R*(2^-p+eta). For any fixed response seed y, the sliced response is a linear map of X with rank at most ell*p. Condition on any of its possible values v. Resiliency makes the challenge eta-close to a mixture with a uniform row. Matching the whole slice to v then has probability at most2^-p+eta, since matching that component's good row alone has probability2^-p. Average over v, and union bound over all R seeds.

All source restrictions here are affine and have explicitly charged rank deficiency. Property3 is a rare-event lower bound, not a high-probability classification of arbitrary blocks. Property4 is a conditional somewhere-randomness hypothesis, not a claim that every nonempty block automatically passes. The final disperser proof must establish that hypothesis for its recursively chosen challenges and compare its error with the rare FAIL probabilities. Replacing R by n^2 without checking the actual seed length would be unjustified; only its polynomial-in-n bound is used in the intended asymptotic parameter comparisons.


# Correctness proof for the polynomial-entropy tree construction

This proof uses the exact tree construction above, implemented by family.py, and the complete component arguments in cubic.md, extractors.md and cubic.md. The construction and the two-stage source analysis follow Shaltiel2011; the extractor components are the separately supplied replacements.

Write H for tree depth and H'=H/q for the number of active intervals. All statements in this proof concern sufficiently large actual N, with n,t,L,k,Kmin,ell,Pmax as defined in the construction. An affine subsource means conditioning a uniform affine source on possible values of binary linear maps. Its deficiency is the rank lost, and its probability inside the original source is2^(-deficiency). Linear images stay affine. Conditioning on a rank-D map decreases any linear-image entropy, and any conditional entropy of linear images, by at most D. To verify the conditional assertion, write H(A|B)=dim image(A,B)-dim image(B); the first term drops by at most D and the second cannot increase. Entropy chain rules are therefore exact dimension identities in this setting.

## Numerical feasibility

There is no outer padding: n=N and L=ceil(log2 n). Consequently

 log2 t=Theta(L^(3/10)), q=Theta(L^(1/2)),
 H=Theta(L^(7/10)), H'=Theta(L^(1/5)),
 k=n^(1/50)+O(n^(1/100)), Kmin=n^(1/100)+O(1).

The shortest leaf has length floor(n/t^H)>=k^(1/4). Eventually H+q<=J: otherwise the largest leaf would have length less than t^(q+1)+1, contradicting k^(1/4)>>t^(q+1). Maximality then gives floor(n/t^(H+q))<k^(1/4), so every leaf has length less than t^q*(k^(1/4)+1)+1. Since t^q=n^o(1), this is below k^(1/2) eventually. Unequal lengths cause no issue in the source splitting below, which uses only the ordered partition and its entropy chain rule.

The block-source component has b=O(1), as Kmin>=n^(1/100). Hence H*binom(t,b)<=H*t^b<=ell eventually. Its output floor(Kmin^(1/10)) exceeds Pmax, and its error is below2^-Pmax, because every fixed positive power of n eventually exceeds2^(L^(4/5)). Likewise ell*Pmax fits within floor(Kmin^(1/8)) response bits.

The logarithm of the largest confidence is at most(6H'+6)q=O(L^(7/10)), below log2 Pmax=ceil(L^(4/5)). Thus every numerical construction guard succeeds eventually. Also

 H*t*ell*Pmax=n^o(1),
 (H+1)*ell*Pmax=n^o(1),
 log2 Rresp=O(L).

Every confidence used is much larger than log n. Finally 4^H', (4b)^H and t^q are all n^o(1). These comparisons hold uniformly over the finite ranges of levels, not just at one chosen node.

## Locate nested useful nodes on an affine subsource

Fix any affine source X on n bits with H(X)>=k. Say X is nice for an interval s at threshold K if all coordinates before s are fixed and H(X_s)>=K. Initially s is the root, so this holds with K=k.

The following elementary split is useful. Partition a nice block s into T consecutive children. Let k_i be the entropy of child i conditioned on all earlier children. Their sum is at least K. Choose the first i with k_i>=K/(4T), and fix all earlier children to any possible joint value. Their entropy was below K/4, so the new source loses at most K/4 dimensions. It is nice for child i, whose entropy is at least K/(4T). If that entropy is at least K/4, continue with threshold K/4. If it is below K/4, the remaining conditional entropy inside s outside that child is at least3K/4-K/4=K/2. Thus both the child's entropy and H(X_s|X_child) are at least K/(4T). All these statements concern possible affine fibers, so no probabilistic source selection is assumed.

Apply this split on consecutive ACTIVE levels, with T=t^q. Each continuing step reduces its threshold by at most4. If it continued to the last active level, the leaf entropy would be at least k/4^H', greater than k^(4/5), while its length is below k^(1/2). This is impossible. At a stopping step obtain consecutive active nodes g and g', with g' below g, and an affine subsource on which coordinates before g' are fixed and

 H(X_g')>=k^(4/5),   H(X_g|X_g')>=k^(4/5).

Indeed the stopping threshold is at least k/(4^H'*4t^q), which exceeds k^(4/5) by the numerical estimates.

Now descend from g' one ordinary level at a time. Suppose the current nice interval s has entropy at least K. If at least b child conditional entropies are >=K/(4t), those b children, in order, are a hidden b-block affine block-wise source at that threshold: conditioning on fewer predecessor blocks cannot decrease conditional entropy. Otherwise, fewer than b indices contribute high conditional entropy; all low indices together contribute at most K/4. Some high index therefore has conditional entropy at least K/(4b). Fix its preceding children and continue in it with threshold K/(4b). These new constraints are inside g', because coordinates before g' were already fixed. They preserve H(X_g|X_g') exactly: conditioning on functions of X_g' selects its fibers but does not change their common affine dimension.

This descent cannot reach a leaf. After at most H steps its retained threshold is at least k^(4/5)/(4b)^H>k^(3/4), larger than leaf length. Thus at an INTERNAL node e below g', a hidden b-block block-wise source is found. Its threshold is at least k^(4/5)/((4b)^H*4t)>=k^(3/4). We have an affine subsource X0 with

 e subset g' subset g;
 g,g' on consecutive active levels;
 H((X0)_g|(X0)_g')>=k^(3/4);
 (X0)_e hidden b-block affine block-wise at threshold k^(3/4);
 all coordinates before e fixed.

## Force empty left siblings to fail

For every ancestor s of e and every child z of s lying before e, X0_z is fixed. There are at most H*t such pairs. Iteratively apply the exact-deficiency empty-block property of challenge-response to node input X_s and challenge ch_z, at confidence Pmax. It fixes ell*Pmax linear response coordinates and makes that test FAIL identically. An affine fiber in the node projection pulls back to an affine subsource of the full source with the same deficiency. At each stage H(X_s)>=H(X_e)>=k^(3/4)-H*t*ell*Pmax>Kmin, so that property applies.

The total deficiency is at most H*t*ell*Pmax=n^o(1). Previously forced FAIL tests persist pointwise under later restrictions, as all challenge and response functions are fixed deterministic maps. Let X* be the final affine subsource. It still has the above nested nodes, with both useful entropy thresholds at least k^(2/3). All left-sibling tests at every ancestor of e FAIL identically. Coordinates before e remain fixed.

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

eventually, since H*log2(2R+2)=O(L^(17/10))=o(Pmax). In particular every true-path challenge has resiliency at least ell*Pmax. Applying the PASS bound at each edge to e and union bounding shows the selected path reaches e, and hence g', with probability at least1-epsilon_path, where

 epsilon_path<=H*R*(2^-Pmax+eta_H)<=2^(-Pmax/3).

All probabilities are with respect to the single fixed affine subsource X*. Left siblings never interfere; the true child is chosen whenever its test passes.

## Both output values occur

Let g be at active level uq. For a=1,2 define FIXED-node tests

 T_a=CR_(ch_g',p_u^(a))(X*_g,g').

On the event that the selected path reaches g', these are exactly g's two output tests. Write p1=p_u^(1), p2=p_u^(2), so p2=ell^3*p1. Conditional entropy H(X*_g|X*_g')>=k^(2/3)>=Kmin gives

 Pr[T_a=FAIL]>=2^(-ell*p_a).

The resilient PASS bound also gives Pr[T2=FAIL]<=R*(2^-p2+eta_H). Monotonicity makes T2=FAIL imply T1=FAIL. Consequently the event(T1=FAIL,T2=FAIL) has probability at least2^(-ell*p2). The event(T1=FAIL,T2=PASS) has probability at least

 2^(-ell*p1)-R*(2^-p2+eta_H)
   >=(1/2)*2^(-ell*p1)>=2^(-ell*p2)

for sufficiently large n. Here p2/(ell*p1)=ell^2 tends to infinity; log R=O(log n) is negligible compared with p2; and Pmax dominates every confidence exponent.

Let D=ell*p2=ell^(6(H'-u)+7). For any active strict ancestor of g, its first confidence is at least ell^2*D. The probability that ANY such fixed true-path test FAILs is at most

 H'*R*(2^(-ell^2*D)+eta_H) <=2^(-ell*D)

eventually. If g is the root there are no such tests and the probability is zero. Also Pmax/3 is much larger than every possible D, since log2 D=O(L^(7/10)) whereas log2 Pmax=Theta(L^(4/5)). Therefore

 2^-D -epsilon_path-2^(-ell*D)>0.

For each desired bit v, intersect the corresponding event(T1=FAIL,T2=v) with the event that the path reaches g' and all active ancestors of g have first test PASS. The displayed inequality proves positive probability. On that intersection, g is exactly the first active node allowed to supply an output, and supplies v. Thus both Boolean outputs occur on X*. Since X* is a nonempty affine subsource of the original X, both occur on X as well.

This proves that the algorithm on the actual n=N coordinates is an affine disperser for every dimension-k affine source, at every sufficiently large integer N. The family is fixed and every chosen subspace was a proof device, not an input-dependent modification of the evaluator.

## Actual-length full-B2 consequence

The tree already acts on all N input bits, without outer padding. Since k=ceil(N^(1/100))^2<=N eventually, the refined full-B2 theorem supplied in elimination.md gives

 C_B2(F_N)>=4N-8k+4
          >=4N-8N^(1/50)-16N^(1/100)-4.

For the second inequality put z=N^(1/100) and use k<=(z+1)^2. In particular this is at least4N-9N^(1/50) whenever z>=17, since z^2>=16z+4 there. The complete dispersion threshold may be larger; no explicit threshold for the entire theorem is claimed. The deterministic total evaluator is uniform P, hence NP with empty witnesses, and all conclusions concern this same fixed family.
