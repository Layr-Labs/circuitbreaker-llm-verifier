# Locally evaluable linear binary code with a polynomial list bound

This is the standard polynomial-evaluation/trace-concatenation idea, with only a list SIZE bound required; efficient list decoding is not used by the extractor evaluator.

Let n>=1 and let0<rho<=1/4 be rational. Choose a prime integer s with q=2^s>1024(n+1)/rho^6. Use the deterministic prime-degree binary field construction to form F_q. Pack the n input bits, followed by at most s-1 zero bits, into coefficients of a polynomial p in F_q[T]. Its degree is at most D=max(1,ceil(n/s)-1). In particular D<=n+1 and D<q. Define the binary word indexed by(a,b) in F_q^2 as

 C_p(a,b)=Tr_(F_q/F2)(b*p(a)).

This encoding is binary linear and injective. Evaluation at any index takes polynomial time in n and s, by Horner evaluation and Frobenius trace. No array of q^2 entries is constructed. Prime search and field construction are polynomial in s for the chosen prime degree, assuming the supplied coarse prime-interval and constructor arguments. One may choose the least prime s exceeding the displayed binary-log threshold; then s=O(log(n+1)+log(1/rho)), eventually and globally after enlarging a fixed finite constant. Thus an index uses2s bits.

**List bound.** Any binary received word has agreement at least1/2+rho with at most q distinct encoded messages.

For a received word r(a,b), put h_a(b)=(-1)^r(a,b), and for v in F_q define c_a(v)=E_b h_a(b)*(-1)^Tr(bv). Agreement implies E_a c_a(p(a))>=2rho. Since each summand is at most1, the fraction of a for which c_a(p(a))>=rho is at least rho: if a fraction theta had this property, the average would be at most theta+(1-theta)*rho, and theta<rho would make that smaller than2rho. For each a, define V_a={v:c_a(v)>=rho}. Orthogonality of the q trace characters gives sum_v c_a(v)^2=1, so |V_a|<=rho^-2. Every good message therefore satisfies p(a) in V_a on at least rho*q field elements.

Set T=floor(rho*q/2). Choose a nonzero bivariate polynomial Q(A,Y) of weighted degree at most T, giving A weight1 and Y weightD, that vanishes on all pairs(a,v) with v in V_a. There are at most q/rho^2 homogeneous linear constraints on its coefficients. The monomials A^iY^j with 0<=i<=floor(T/2), 0<=j<=floor(T/(2D)) all have weighted degree at most T, so their count is at least T^2/(4D). Since rho*q is much larger than1, T>=rho*q/4. The number of available monomials is consequently at least rho^2*q^2/(64D)>q/rho^2, by q>1024(n+1)/rho^6 and D<=n+1. Linear algebra therefore gives such a nonzero Q.

For any good message, Q(A,p(A)) has degree at most T and vanishes on at least rho*q>T points. It is the zero polynomial. Polynomial division by the monic polynomial Y-p(A) over F_q[A] shows that Y-p(A) divides Q(A,Y). Distinct messages give distinct monic factors; viewing Q over F_q(A)[Y], the number of such factors is at most deg_Y Q<=T/D<q. If deg_YQ=0 there are no good messages. This proves the stated (looser) bound q in every case. No decoding algorithm, factor-extraction complexity or probabilistic construction is assumed.

A numerical detail to retain in later use: agreement and correlation differ by a factor2. The threshold at each row is rho, while global correlation is2rho. The list size costs at most s bits of advice. The code's evaluation complexity depends on the FIELD BIT LENGTH s, not on its exponentially larger number q of elements.


# Linear strong extraction from a design and a list-bounded code

This section proves the standard Trevisan reconstruction principle. The encoder lemma is the corresponding section of extractors.md. The following sections supply the explicit designs, parameter regimes and iteration.

Let C be an injective binary linear encoding of n-bit messages whose code positions are indexed by ell bits. Suppose at most L codewords agree with any received word on a fraction at least1/2+rho, where rho=epsilon/(2m), 0<epsilon<=1/4 and m>=1. Let S1,...,Sm be specified ell-element subsets of a d-element seed universe, ordered internally in a specified way. Put

 B=max_(1<=i<=m) sum_(j<i) 2^|Si intersect Sj|.

Define Ext(x,y)_i=C(x)[y restricted to Si]. If a source X on n bits has min-entropy at least

 k >= d+B+log2 L+log2(4m/epsilon),

then (Y,Ext(X,Y)), for independent uniform d-bit Y, is within total variation epsilon of uniform(Y,U_m). For each fixed seed the map x->Ext(x,y) is binary linear.

Proof. Suppose the total variation exceeds epsilon. A Boolean test D(y,z) distinguishes the two distributions with probability gap exceeding epsilon. Insert the m hybrids that replace suffix coordinates of Ext by independent uniform bits. One adjacent pair has signed gap exceeding epsilon/m; complement D if necessary to make it positive. Fix this coordinate i and let d_b(y,z_<i) be the conditional acceptance probability after assigning coordinate i the bit b and randomizing the later output bits. The difference between the adjacent hybrids is exactly

 (1/2)*E[(-1)^Ext_i * (d_0(Y,Ext_<i)-d_1(Y,Ext_<i))].

Consider a randomized predictor that outputs0 with probability(1+d_0-d_1)/2 on input(y,z_<i). Its prediction advantage over1/2 is precisely the displayed difference. Averaging its finite table of random choices permits fixing one deterministic predictor P with advantage alpha>epsilon/m over X,Y. No efficiency is required of P; it is fixed once and does not depend on the individual source sample.

Write a(x)=Pr_Y[P(Y,Ext(x,Y)_<i)=Ext(x,Y)_i]-1/2. It lies in[-1/2,1/2] and has expectation alpha. The set G of messages with a(x)>=alpha/2 has X-probability at least alpha/2, hence greater than epsilon/(2m). Indeed if its probability were below alpha/2, the expectation would be at most (alpha/2)*(1/2)+(1-alpha/2)*(alpha/2)<alpha. This intentionally uses a loose bound.

For each x in G, fix the seed coordinates outside Si to an assignment z for which the conditional prediction advantage remains at least alpha/2. Supply z as at most d advice bits. For each j<i, with these outside coordinates fixed, the previous code bit C(x)[Y|Sj] is a function only of the coordinates in Si intersect Sj. Supply its whole truth table as at most2^|Si intersect Sj| bits, with a fixed lexicographic table order. All tables together use at most B bits. Tables may be padded to a common fixed advice length.

Given any such advice, each t in{0,1}^ell determines the entire seed by combining t with z, and the tables determine every preceding output bit. Running the fixed predictor P therefore defines one received word r(t). For the particular advice belonging to x, this word agrees with C(x) on at least1/2+alpha/2>1/2+epsilon/(2m)=1/2+rho of positions. By the code hypothesis at most L different messages have this property for any fixed advice. There are at most2^(d+B) advice strings. Thus

 |G|<=2^(d+B)*L,
 Pr[X in G]<=2^(d+B)*L*2^-k<=epsilon/(4m),

contradicting the previously established strict lower bound epsilon/(2m). This proves strong extraction, with the seed included in the tested distribution. The seed was not assumed independent after fixing message-specific advice; independence is used only in the original expectation and hybrid argument.

Finally each output code bit is linear in x. Fixing a seed fixes its index, so all output bits jointly form a linear transformation of x. All explicit evaluators need only the locally computable encoder and the specified design; the hypothetical predictor, advice, received words and list are proof devices and are never computed by the family.


# Constant-error linear extractor for polynomial entropy, with logarithmic seed

Fix any rational a>0. For input length n and entropy parameter k>=n^a, there is a deterministic polynomial-time construction of a linear strong extractor with error1/4, output m=floor(k^(1/8)) and seed length O_a(log n), at every sufficiently large length. Finite small lengths can be assigned a specified total fallback. The proofs below invoke only the explicitly supplied code and reconstruction lemmas, field constructor and prime-interval argument.

Use epsilon=1/4 and rho=1/(8m) in the corresponding section of extractors.md, with local code index length ell=2s and list bound L=2^s. Let b=floor(log2(k)/4), t=b+1, and d=ceil(32ell^2/t). Eventually1<=t<=ell and d>=2ell. Construct m ordered ell-subsets of[d] greedily in lexicographic order, requiring all pair intersections to be at most b.

For any fixed ell-set S, a uniformly random ell-set T has intersection at least t only if it contains one of the binom(ell,t) t-subsets of S. For each particular t-subset the inclusion probability is at most(ell/(d-ell))^t<=(2ell/d)^t. The elementary factorial bound t!>=(t/3)^t gives

 Pr[|S intersect T|>=t]
   <=(3ell/t)^t*(2ell/d)^t
   <=(3/16)^t <4^-t.

(The factorial bound follows, for example, by summing log j>=integral_1^t log x dx, which gives t!>=(t/e)^t, and e<3.) At any greedy step the union of forbidden events has probability less than m*4^-t. Since t>log2(k)/4 and m<=k^(1/8), this is at most k^(-3/8)<1. A next set always exists. The intersection rule also excludes duplicates.

The reconstruction advice bound is B<=m*2^b<=k^(3/8). Prime search gives s=O(log n+log m)=O(log n), so ell=O(log n), d=O(ell^2/log k)=O_a(log n), log2L=s=O(log n), and log2(4m/epsilon)=O(log n). Therefore

 d+B+log2L+log2(4m/epsilon) < k

for all sufficiently large n. The general reconstruction lemma proves strong error at most1/4. Every output coordinate is a locally evaluated linear code bit, so the extractor is linear in its source input for fixed seed.

A deterministic constructor can enumerate all at most2^d subsets of[d], testing their cardinality and intersections, until the next accepted ell-set is found. Because a is fixed, d<=C_a log n eventually; this enumeration and the m local code evaluations have polynomial running time with a FIXED exponent. This is a proof of polynomial time, not a claim that the exponent is small. Reject parameters k<n^a with the specified fallback. For smaller n, use a fixed computable parameter guard ensuring the displayed entropy inequality and existence estimates; if they fail, return a specified zero-output routine. Such guarded lengths do not enter the eventual guarantee. Do not apply this greedy construction at merely polylogarithmic entropy and still call its seed or enumeration logarithmic/polynomial; polynomial entropy is essential here.

This component can replace an O(log n)-seed constant-error call only when its output demand is at most k^(1/8). The separate exponentially small-error transfer extractor has different design and seed requirements. It is not established by this lemma.


# Small-error linear extraction with a polynomial seed budget

Fix a rational a>0. Let a source consist of N bits, with entropy at least r/2+floor(r/4), where N^a<=r<=N. For all sufficiently large parameters, construct a linear strong extractor whose output length is M=floor(r/4), error at most2^(-r^(1/1000)/3), and seed length O(r^(1/10+4/1000)). All source hypotheses here are AFFINE when output iteration is used. Integer powers in the construction mean exact floor roots of integer powers; floating-point arithmetic is unnecessary.

First make a smaller-output extractor. Put z=floor(r^(1/1000)), epsilon=2^-z, m=floor(r^(9/10)), rho=epsilon/(2m). Use the locally evaluated code with index length ell=2s and list bound2^s from the corresponding section of extractors.md. Then ell=O(log N+z+log m)=O(r^(1/1000)), and ell=Omega(r^(1/1000)), because r is polynomially large in N.

Choose a prime degree p>=ceil(log2 ell), and construct F_Q with Q=2^p. Choose the least such prime; the supplied prime-interval proof ensures Q<=ell^2 eventually, after a harmless fixed factor if the lower endpoint is rounded. One may use Q<=4ell^2 throughout to absorb that factor. Let A=ceil(log_Q m), meaning the least integer A with Q^A>=m. Eventually A is bounded by a FIXED constant, since log Q>=log ell=Omega(log r). Also A<=ell. Encode integers0,...,m-1 as distinct coefficient vectors of polynomials P_i in F_Q[T] of degree at most A-1.

For the first ell field elements alpha in their binary integer order define the ordered set

 S_i={(alpha,P_i(alpha)):alpha among those ell elements}.

The universe is F_Q^2, of size d=Q^2=O(ell^4). Each set has EXACTLY ell elements, not Q padded elements. Distinct polynomials agree at at most A-1 field points, by the root bound. Hence the reconstruction advice B is at most m*2^(A-1)=O(r^(9/10)), with a fixed, possibly large constant. The other costs d, s and log2(4m/epsilon) are O(r^(4/1000)). Their sum with B is o(r), therefore below r/2 eventually. The reconstruction lemma proves strong error at most epsilon. For fixed seed this extractor is linear in the source input.

Produce M bits by J=ceil(M/m) successive applications with independent fresh seeds, trimming the final output block. At each stage, condition on all preceding seeds and output blocks. For an affine initial source, this is either an impossible history or a uniform affine source whose entropy dropped by no more than the number of preceding output bits: every map is linear once the previous seeds are fixed. Thus its entropy stays at least r/2 as long as the initial entropy was at least r/2+M. In the applications below the actual initial entropy is at least r-r^(1/2), which exceeds that amount eventually. More generally the iterated guarantee needs initial entropy at least r/2+M; do not erase this extra requirement.

For that initial entropy, the strong-extractor conclusion applies to every possible history. Replacing the successive output blocks by fresh independent uniform blocks, one at a time, incurs total variation at most J*epsilon by the triangle inequality and conditioning. Therefore the concatenated output is a strong extractor with error at most J*2^-z<=2^(-r^(1/1000)/3) eventually. Its total seed length is J*d=O(r^(1/10+4/1000)), and it is linear for every fixed complete seed. Every design and code bit is explicitly locally evaluated in time polynomial in N and r. Huge numbers q in the code specify index spaces only; their full tables are never enumerated.



# Affine somewhere-random extraction from the supplied components

This is a modified Rao-style row condenser. It uses the vector cubic extractor section of cubic.md, the small-error seeded extractor below, linear algebra, and total-variation inequalities. All asymptotic parameter estimates concern sufficiently large lengths; family.py supplies total numerical fallbacks. The source is uniformly distributed on an affine subspace of binary t-by-r matrices, and at least one r-bit row is uniform. Assume1<=t<=r^(1/100).

## One row-condensation step

For t>=2 put w=floor(r^(2/5)), M=floor(r/4). Take the first w bits of each row. Pair these slices consecutively, duplicating the last row if t is odd. Apply the half-entropy cubic extractor to every pair; trim its output to the seed length of the iterated small-error linear extractor on the FULL tr-bit input, with parameter r and output M. Use each resulting seed to apply that linear extractor to the full matrix. The new matrix has ceil(t/2) rows of M bits.

The seed length required is O(r^(1/10+4/1000)), whereas the cubic extractor produces Theta(r^(1/5)) bits, so trimming is permitted eventually. Construction is deterministic; the seeds are functions of the input, and their validity will follow from conditioning rather than an unjustified independence assertion.

Let Z be all tw slice bits. The affine source X splits as A(Z)+B, where B is uniform on the kernel of the slice map and is independent of Z. To obtain this decomposition, choose a basis of that kernel and extend it to a basis of the direction space of X. The remaining basis coordinates biject with the image of the slice map, hence are determined by Z; assign the affine offset to A. Since entropy(X)>=r and rank(slice)<=tw,

 H(B)>=r-tw.

A specified pair contains an originally uniform row, so its slice pair is an affine source of entropy at least w. Its cubic output seed is within beta=2*2^(-sqrt(w)/50) of uniform. It is a function of Z and thus independent of B. This statement concerns the distribution of that seed before fixing Z; it does not assert that a fixed seed value remains uniform.

For a FIXED seed u, the linear extractor applied to B has affine output. Its distance from uniform is therefore either zero or at least1/2 (a proper affine subspace has at most half the points). Strong extraction with uniform seed and error eta implies that the fraction of seeds with nonuniform output on B is at most2eta. The actual nearly uniform seed has bad probability at most2eta+beta. Here eta<=2^(-r^(1/1000)/3), since H(B)>=r-tw eventually exceeds r/2+M: t<=r^(1/100) and w<=r^(2/5) give tw<=r^(41/100)=o(r).

Condition on all slice values Z=z. Every seed is now fixed; the entire new matrix is an affine image of B plus a fixed shift. For every z at which the specified good pair supplied a good seed, its corresponding new row is uniform. Thus, except for probability at most2eta+beta, the conditional output is an affine somewhere-random source. Replacing exceptional conditional distributions by arbitrary affine somewhere-random ones proves that the unconditioned output is within2eta+beta of a convex combination of such sources. All rows, including possibly bad ones, stay jointly affine after conditioning; this is the feature needed for iteration.

## Iterate until one row remains

Repeat this procedure until the row count is1. The next row length is floor(r_j/4), while t_(j+1)=ceil(t_j/2). The number of steps is J=ceil(log2 t). If initial r is large and t<=r^(1/100), then at every step

 r_j>=r/(2*4^j)>=r/(8t^2),
 t_j<=t<=r^(1/100)<=r_j^(1/50),

where repeated integer division gives r_j=floor(r/4^j)>=r/4^j-1, and r/4^j grows without bound. The one-step proof above used only tw=o(r_j), the seed-budget exponent gap, and polynomial entropy relative to total input length. It continues to hold for the relaxed t_j<=r_j^(1/50), since t_j*w_j<=r_j^(21/50)=o(r_j). In particular r_j is still polynomially large in the length t_j*r_j of its input matrix.

Each deterministic step preserves total-variation closeness, and each component of a mixture obeys the same uniform one-step error bound. Summing errors over the J steps gives at most2^(-r^(1/2000)) eventually: the smallest r_j is at least r^(98/100)/8, so its r_j^(1/1000) dominates r^(1/2000), and J<=log2 r contributes only a logarithmic factor. The cubic-seed error has a still larger exponent. The final single row is uniform in the ideal mixture, and has at least r/(8t^2)>=r^(97/100) bits eventually. Output its first floor(r^(97/100)) bits.

Therefore, for all sufficiently large integer r and1<=t<=r^(1/100), the construction extracts floor(r^(97/100)) bits from every affine t-by-r somewhere-random source, with error at most2^(-r^(1/2000)). For t=1 simply return the requested prefix of its uniform row. Every field operation, local code evaluation, design construction and number of rounds has polynomial running time in tr; no exponentially large code table or received-word list is constructed.

The following block-source construction uses this row extractor. The complete challenge-response and tree arguments are in proof.md.


# A block-source extractor with supplied polynomial-entropy components

This component follows Shaltiel's affine block-source reduction. Fix a=1/5. Let each input block have n bits, and let integer K>=n^a. The joint source (X1,...,Xb) is affine and has conditional entropy at least K in each block after fixing any possible values of its predecessors. Conditional entropy is constant over possible affine fibers, so this is a pointwise condition, not an average entropy assumption.

Use the constant-error linear extractor on one block at entropy K and error1/4. Trim its output to R=floor(K^(1/10)) bits; its seed has d0=O(log n) bits. Enumerate all v=2^d0 seeds, obtaining a v-by-R matrix W1 of extractor values on X1. This is an affine matrix source. At least one row is uniform: for fixed seed the output is affine, nonuniform affine output has distance at least1/2, and average distance at most1/4 precludes all seeds being bad.

Let g=floor(R^(1/100)), eventually at least2. Let E be the SMALLER-output small-error linear extractor at input length n and entropy parameter K, trimmed to R output bits. Its original output floor(K^(9/10)) is longer than R, its seed length is d1=O(K^(4/1000)), and its strong error is epsilon1=2^-floor(K^(1/1000)). It needs entropy only K/2, hence applies to the independent parts below. Do not substitute the iterated near-linear-output routine here; the short seed bound is for the single small-output call.

To turn Wi with u rows into Wi+1, partition its rows into ceil(u/g) consecutive groups, padding the last group with zero rows to g if necessary. Apply the corresponding section of extractors.md to each g-by-R group, trimming its output to d1 seed bits. This is possible eventually: its output length is floor(R^(97/100))=K^(97/1000) up to floors, much larger than d1=O(K^(4/1000)). Apply E to the NEXT input block using each resulting seed, and use the R-bit results as the new matrix rows.

Choose b=1+min{j:g^j>=v}. Then b=O(log v/log g)=O(log n/log K)=O_a(1); all quantities and the minimum are computable integers. After b-1 transfers, one row remains. Output it, or any specified shorter prefix.

## Conditional correctness invariant

For each i there is a set of histories (X1,...,X_(i-1)) of probability at least1-(i-1)*tau such that, conditioned on every such history, Wi is an affine somewhere-random matrix. Here tau=2epsilon1+epsilon2, where epsilon2<=2^(-R^(1/2000)) is the group extractor's error. For i=1 the history is empty and the initial construction established the invariant exactly.

Fix a good history before step i. Given that history, Wi is a linear function of Xi plus a fixed shift: all seeds used at previous steps become fixed when the earlier input blocks are fixed. The joint pair(Xi,X_(i+1)) is affine. Split it by linear algebra as X_(i+1)=A(Xi)+B, where B is an affine source independent of Xi and has entropy H(X_(i+1)|Xi,history)>=K. This follows by choosing a complement to the kernel of the projection onto Xi; it holds even if the blocks are correlated.

One group in Wi contains a uniform row, so its extracted seed is within epsilon2 of uniform. It is a function of Xi, independent of B. Strong linear extraction on B implies that at most2epsilon1 of uniform seeds yield a nonuniform R-bit output on B. Therefore that group's actual seed is bad with probability at most tau. For every Xi at which it is good, conditioning additionally on Xi fixes all next seeds and A(Xi). The new matrix is jointly affine as a function of B, with the specified group's row uniform. This establishes the invariant at step i+1 after discarding probability at most tau of each previously good history. The invariant concerns FULL previous input-block histories, so it does not assume that an extractor output reveals all of its source block.

At the last step, outside probability (b-1)*tau, the single output row is uniform conditioned on the full previous-block history. Its unconditioned distance from uniform is therefore at most(b-1)*tau. Since b is bounded for fixed a, R=K^(1/10) up to floors, and epsilon2 has exponent K^(1/20000) up to constant factors, the final error is at most2^(-K^(1/100000)) for all sufficiently large K. Output length is floor(K^(1/10)).

## Running time and scope

Enumerating2^d0 seeds is polynomial in n because d0=O(log n) for K>=n^(1/5). Every transfer reduces the row count by a factor g up to a ceiling, uses a polynomial-time group extractor and local linear code evaluations, and there are b=O(1) transfers. The total evaluation and construction time is polynomial in n with a fixed exponent. The construction does not enumerate a codeword's exponentially large index space. Input blocks of length below n may be padded with zeros, preserving affinity and entropy.

This proves the stated block-source component under its conditional entropy hypothesis. The total tree construction, source restrictions and actual-length transfer are supplied in proof.md.