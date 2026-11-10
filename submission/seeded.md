# Supplied extraction arguments for the polylogarithmic branch

All dimensions and entropies are binary. A1 is proved for any good field/matrix tuple. Internal calls share one actual input-prefix tuple by C6; C1 constructs the initial/B5 tuples canonically. The deterministic fields in B1 require no such data.
Charged uses fix good internal prefix data for hardness while remaining total on all prefix strings. C1 instead constructs initial/B5 good data canonically; those calls consume no input-data bits.

# A1. a strong linear seeded extractor

## Exact construction and parameters

Let the source length be n bits. Fix positive integers D,h,m,s with
binomial(h+D-1,D)>=n, and let q=2^s. Let K be any field of size2^s with a specified binary basis whose first vector is1; for example K=F_2[T]/Q for a monic irreducible binary Q of degree s. Let A be a D-by-D matrix over K
whose action cycles through K^D minus zero. Such matrices exist, as
justified below. This property is a mathematical hypothesis on A; the
generic A1 evaluation does not require a search or test of A. C1 separately supplies bounded canonical searches for selected calls; C3 supplies the internal reconstruction.

Enumerate all D-variate monomials of total degree below h, first by total
degree and then lexicographically. For x in F_2^n, let P_x be the K-valued
polynomial with its first n coefficients x_1,...,x_n and all remaining
coefficients zero. This is an injective F_2-linear encoding with degree
at most h-1. The seed consists of v in K^D and z in F_2^s. Output

    E(x;v,z)_i = <bits(P_x(A^i v)),z>,   1<=i<=m.

The bracket is the binary dot product. For each seed this is linear in
x. Evaluation takes polynomial time in n,m,D,h,s and the number of
monomials, without enumerating K^D or finding a primitive element.

Take 0<epsilon<=1/4 with rho=epsilon/(8m) reciprocal to an integer.
Choose any integer r>D with r>=32. Sufficient hypotheses are

    q^(r-D) >= (2^20*r*h/rho^4)^r,
    k >= 16*m*r*h*s + 4*ceil(log2(1/rho)) + 10.             (1)

Claim to prove: every source X on n bits with min-entropy at least k
satisfies TV((U_seed,E(X,U_seed)),(U_seed,U_m))<=epsilon. This is strong,
not merely a statement about the output marginal.

## Two elementary counting tools

1. Polynomial list recovery. Let q points each have a list of at most L
field values. A degree-at-most b polynomial is accepted if its value
belongs to the list on at least a points. If a^2>=2bqL, at most 2qL/a
polynomials are accepted. Indeed, for T accepted polynomials count
coinciding accepted evaluations in qL buckets. Cauchy--Schwarz gives at
least ((Ta)^2/(qL)-Ta)/2 pairs. Distinct polynomials agree at no more
than b positions, giving at most T(T-1)b/2. Thus
T(a^2/(qL)-b)<=a and T<=2qL/a. No efficient list-recovery algorithm is
needed: the reconstruction below may enumerate all polynomials.

2. Limited-independence lower tail. If q Bernoulli variables are 2r-wise
independent, have common expectation p and mean mu=pq>=r, then

    Pr[sum Y_j < mu/2] <= (C*r/mu)^r

for an absolute C (one can take a loose C=2^12). Expand the 2r-th
central moment. Terms with an index appearing once vanish; remaining
terms have at most r distinct indices. Grouping by their multiplicities
and using |Y-p|<=1 and E|Y-p|^a<=2p for a>=2 bounds the moment by
(C*r*mu)^r when mu>=r. Markov at (mu/2)^(2r) gives the stated form after
increasing C. For an explicit check, write Z_j=Y_j-p. Since |Z_j|<=1 and
E|Z_j|^a<=p for a>=2, terms using exactly ell distinct indices contribute
at most binomial(q,ell)*ell^(2r)*p^ell. For 1<=ell<=r and mu>=r,
mu^ell/ell!<=mu^r/r!. Thus the moment is at most
r*r^(2r)*mu^r/r! <= (2e*r*mu)^r, using r!>=(r/e)^r and r<=2^r.
Markov contributes a factor 4^r, so C=32 already suffices; C=2^12
is safe. The factorial bound follows by integrating log on [1,r].

## A q-ary list-predictor contradiction

Suppose a predictor F receives a seed v and prior field symbols
P_x(A v),...,P_x(A^(i-1)v), and returns a list of at most rho^-2 field
elements. Assume its list contains P_x(A^i v) with probability greater
than rho over x drawn from X and uniform v. Then the set of source
strings on which its probability is at least rho/2 has mass greater
than rho/2. Fix such an x. We reconstruct x with short advice and
success probability at least 1/2 over coins independent of x.

At u=A^i v, the predictor can be evaluated with known u and values
P_x(A^(-i+1)u),...,P_x(A^-1 u). Its success density on K^D is at least
rho/2. Dependence on the seed is harmless: the reconstruction explicitly
passes A^-i u to F.

Choose 2r random distinct field elements t_1,...,t_(2r), and independent
uniform points y_1,...,y_(2r) in K^D. Interpolate degree-at-most 2r-1 curves
p_1,p_2 such that

    p_1(t_j)=y_j for all j,
    p_2(t_j)=A y_j for j<=r, and y_j for j>r.

Define C_(2j+1)=A^j p_1 and C_(2j+2)=A^j p_2 for j>=0. Each curve is
uniform among degree-at-most 2r-1 curves, so its evaluations at distinct
field arguments are 2r-wise independent uniform points. Consecutive
curves share r random distinct argument positions, independent of the
current curve: the last r positions for p_1,p_2 and the first r for
p_2,A p_1. A matrix power preserves these statements. If y_1 is nonzero,
the curves up to index 2q^D cover all nonzero K^D points, since A cycles
through them. The event y_1=0 has probability q^-D.

Advice gives the coefficients of P_x composed with the first 2i
curves. Each composition has univariate degree at most 2rh, so at most
6mrh field coefficients, or 6mrhs bits, suffice. The reconstruction
knows its random curves and all construction parameters as public data.

Proceed through later curves in order, starting at index 2i+1. The two
extra initial curves ensure that the preceding curve and all earlier
2a-shifted curves are known even when i=1. At a point u of current curve
C_j, the prior points A^-a u lie on C_(j-2a), so their values are already
known. Evaluate F to obtain a list at each field argument of C_j. The
moment bound shows that fewer than rho*q/4 successful predictions has
probability at most

    (2^14*r/(rho*q))^r <= 1/(64*q^D)

under (1). With at least a=rho*q/4 successful positions, the list-recovery
bound applies at b=2rh, L=rho^-2 and yields at most 8rho^-3 candidate
univariate polynomials, including the true composition.

The preceding curve's known values at its r shared argument positions
identify the correct candidate. An incorrect degree-at-most 2rh candidate
agrees at no more than 2rh arguments. The shared r-subset is independent
of the current curve, hence of its candidate list (the fixed true x and
predictor do not introduce additional randomness). Its chance to agree
at all shared positions is at most (4rh/q)^r. Union over candidates
bounds the chance of an incorrect survivor by

    8rho^-3*(4rh/q)^r <= 1/(64*q^D).

For both displayed inequalities, the new exact field condition (1) gives
q^(1-D/r)>=2^20*r*h*rho^-4 directly by taking a positive r-th root. The code chooses s=ceil(r*ceil(log2 base)/(r-D)), where base=2^20*r*h*rho^-4, using integer arithmetic. Thus s*(r-D)>=r*ceil(log2 base) and q^(r-D)>=base^r. Since r>D, this also implies q>=base>=64*r*h*rho^-4. No old eighth-power premise is needed. The first left-hand
side times q^D is at most (2^-6*rho^3/h)^r <=2^-6r<1/64.
The second left-hand side times q^D is at most
8rho^-3*(2^-18*rho^4)^r<1/64, since r>=32 and rho<=1.
Also q>=64*r*h*rho^-4 gives a^2>=2bqL, as needed for list recovery. Union over at most 2q^D steps, plus the
coverage exception y_1=0, leaves success probability above 1/2. On that
event P_x is known on every nonzero vector. These evaluations determine
a polynomial of total degree below h<q-1 uniquely: a nonzero low-degree
polynomial cannot vanish on all nonzero points. For example, restrict
to every coordinate-parallel line avoiding zero and use the root bound,
then recover the missing value on a line through zero. Thus x is recovered.

For the positive-mass set of good x, each has some advice of at most
8mrhs bits yielding recovery with probability at least 1/2. Fixing the
random coins by averaging gives a deterministic reconstruction with at
most 2^(8mrhs) outputs capturing mass greater than rho/4. A min-entropy-k
source gives these outputs mass at most 2^(8mrhs-k), which is less than
rho/4 under (1). Contradiction. Therefore no q-ary list predictor of the
specified size succeeds with probability greater than rho.

## Strong binary output

If the strong binary conclusion failed with distance greater than epsilon,
the standard hybrid argument would give an i and predictor B(v,z,prior
bits) with advantage greater than epsilon/m=8rho for the i-th bit.
For each x,v, view its predictions as a received Hadamard word indexed
by z. Averaging shows that a set of x,v of mass greater than rho has
prediction advantage at least rho over z. For each fixed v and prior
field symbols, form the list of field elements whose binary dot-product
codewords agree with this received word on at least 1/2+rho of seeds z.
By orthogonality of characters and Parseval this list has size at most
1/(4rho^2), hence at most rho^-2. For more than a rho fraction of x,v it
contains the actual next symbol. This contradicts the q-ary property.
The argument passes v and z as seed data throughout and therefore proves
strong extraction, rather than assuming a weak extractor is strong.

## Good field and matrix data exist; input-prefix use

The roots of X^(2^s)-X form a field of size 2^s; the Frobenius-orbit
argument supplies irreducible polynomials of every degree (the union of
proper subfields has fewer elements than the full field). Over K an
extension of degree D exists by the same argument. The multiplicative
group of any finite field is cyclic: in a finite abelian group its exponent
is attained as an element order, and if the exponent were below the group
size, the polynomial X^exponent-1 would have too many roots. Multiplication
by a generator in a K-basis of the extension gives the required A.

With the GENERIC full-matrix encoding, the wrapper reads Q and A as actual input bits,
with their lengths predetermined by the total length. For bad data it
still runs the same bounded ring/matrix algorithm, defining a total
polynomial-time function. A lower-bound restriction fixes good data,
whose existence was proved; all data bits must be paid for in actual
input length. That generic full-matrix encoding requires(D^2+1)s bits. C2 uses companion encoding requiring(D+1)s bits. The CURRENT family uses ONE C5 small-factor unit package shared by all four internal maps via C6, requiring D*s_common-a(D*s_common) bits; C1 computes initial/B5 tuples canonically with zero input-data bits. Every actually read bit is charged in proof.md.




# A2. A strong entropy sampler with polynomial error

This section supplies the exact sampling variant used by B5. It is a distributional statement about a deterministic, seeded bit-selection algorithm; no source probabilities, repairs or good-event flags are computed by the evaluator. Sampling/block-source conversion and backward extraction are related to Nisan--Zuckerman and Zuckerman's randomness-optimal sampling architecture (David Zuckerman, Randomness-Optimal Oblivious Sampling, Random Structures and Algorithms11(4),1997; https://www.cs.utexas.edu/~diz/pubs/sampler.pdf). No theorem from that paper is used as a premise. The arguments needed here are given below rather than assumed from those works.

## Statement and probability-tree bound

Let n be a power of two, 0<sigma<=1 a fixed rational, b=ceil(log2(64/sigma)), and t>=1 an integer. Let t<=ell<=sigma*n/(64b). Select ell positions I_1,...,I_ell in [n], uniformly and 4t-wise independently, allowing repetition. For every source X on n bits with min-entropy at least sigma*n, the joint distribution of the seed and sampled word is epsilon-close to (seed,V_seed), where EVERY conditional V_seed has min-entropy at least

 kappa=floor(sigma*ell/8),

and one may take

 epsilon = 2^(-sigma*n/8)
   +(128*t*b^2/(sigma^2*ell))^t
   +C_t*n^t*(128b/(7sigma*ell))^(2t),
 C_t=2t*(2t)^(4t)*2^(2t).                              (A2.1)

The seed marginal is unchanged and uniform. Thus this is a STRONG smoothed-entropy statement, not merely a bound on the marginal sampled word.

For the proof, use any fixed order of the n input bits. On a supported input x let p_i be its conditional bit probability given the previous bits and w_i=-log2 p_i. Then sum_i w_i=-log2 Pr[X=x]>=sigma*n pointwise. Conditional probabilities on unsupported prefixes can be completed arbitrarily; they do not affect the source distribution. Put T_i=(w_i-b)_+. Conditional on any supported prefix,

 E[2^(T_i/2)|prefix]<=1+2^(1-b)<=1+sigma/32.

Indeed an outcome of probability p<2^-b contributes 2^(-b/2)*sqrt(p), at most2^-b, and there are at most two outcomes. Other outcomes contribute their original probabilities. Iterated conditioning and Markov consequently give

 Pr[sum_i T_i>sigma*n/2]
 <=exp(sigma*n/32)*2^(-sigma*n/4)<=2^(-sigma*n/8).

The last inequality follows from ln2>=1/2, for example by integrating1/x on[1,2]. Outside this event the clipped weights sum to at least sigma*n/2.

For a fixed such x, normalize its clipped weights to z_i=min(w_i,b)/b in[0,1]. Their mean is at least sigma/(2b). The 2t-th centered moment of their sum at 2t-wise independent uniform sample positions is at most(8t*ell)^t. To check the constant, expand the moment: terms containing a singleton index vanish, and terms with j<=t distinct indices contribute at most binomial(ell,j)*j^(2t). Since ell>=t, summing these is at most t*t^(2t)*ell^t/t! <=(2e*t*ell)^t<=(8t*ell)^t. The factorial estimate follows by integrating log, as in A1. Markov at half the mean shows that the sampled clipped weight is below sigma*ell/4 with probability at most the second term of(A2.1). This argument is uniform in x, so it can be averaged over X.

## Repeated positions are explicitly paid

Let C=sum_(j<k)1[I_j=I_k] and mu=binomial(ell,2)/n. We claim

 E[(C-mu)^(2t)]<=C_t*n^t.                              (A2.2)

Its expectation agrees with independent uniform positions because every expanded term uses at most4t position variables. In the independent case expand into ordered sequences of2t centered equality indicators. A term with a vertex appearing in only one indicator vanishes: condition on all other positions and average that centered indicator. Every remaining multigraph therefore has weighted minimum degree at least two, at most2t vertices and at most t components.

For a fixed remaining graph with v vertices and c components, expand each centered indicator into equality minus1/n. For a chosen subset of retained edges, its equality probability is n^(-rank), where rank is v minus the number of components, including isolated vertices. Deleting an edge raises the number of components by at most one. Thus omitted-edge count plus retained rank is at least v-c. Every summand in this expansion has magnitude at most n^(-(v-c)), and there are at most2^(2t) summands. Overcount graph types by(2t)^(4t), label choices by ell^v<=n^v, and sum at most2t possible vertex counts. Since c<=t, this gives exactly(A2.2). It counts repeated edges and shared vertices, not only disjoint pairs.

Our ell bound gives mu<=sigma*ell/(128b). Markov applied to(A2.2) yields

 Pr[C>sigma*ell/(16b)]
 <=C_t*n^t*(128b/(7sigma*ell))^(2t).

The number of repeated draws beyond the distinct selected positions is at most C. Removing repeated clipped weights therefore loses at most bC. Outside both sampling exceptions, the sum of ORIGINAL surprisals over DISTINCT selected positions is at least

 sigma*ell/4-sigma*ell/16>=sigma*ell/8.

This step is necessary: counting the surprisal of a repeated position twice would not prove entropy of the sampled word.

## Turning the good event into conditional min-entropy

Fix a sampler seed and a consistent output word y. Its positions determine a distinct set S and the values at S. Define a modified probability-tree measure Q_y by forcing those selected bits to their specified values with probability one, and following the original conditional probabilities at unselected positions. Complete zero-probability prefixes arbitrarily. This is a probability measure. For every original supported x consistent with y,

 Pr[X=x]=Q_y(x)*product_(i in S)p_i(x_i|x_<i).

On the good event the product is at most2^(-sigma*ell/8). Summing over good supported x proves that the probability of output y AND the good event is at most2^(-sigma*ell/8)<=2^-kappa. Inconsistent duplicate outputs have probability zero. Thus, for EACH fixed seed, the good subdistribution has point masses at most2^-kappa. Fill its removed mass under that same cap to obtain a full distribution: the output alphabet has capacity2^ell*2^-kappa>=1, so sufficient capacity always remains. The total variation cost is at most the removed mass. Fill may use words inconsistent with repeated positions; this is allowed in a comparison distribution, not an operation performed by the sampler. Averaging over the unchanged seed proves(A2.1).

## Exact polynomial-seed implementation and guards

The implementation uses sigma=1/2000,b=17,t=256 and4t=1024. With n=2^w it reads1024 independent field coefficients and evaluates the degree-at-most1023 polynomial at distinct field points0,...,ell-1. Its field values are the selected positions. For any at most1024 distinct arguments, the evaluation map is surjective by the Vandermonde determinant: the product of pairwise argument differences is nonzero. Hence those positions are independent uniform field elements. If fewer arguments exist the entire smaller tuple is independent. The current sampler uses C1's canonical degree-w irreducible modulus, with existence supplied in A1 and polynomial(r) construction proved in C1. The generic sampling routine also remains bounded on arbitrary monic data.

In B5 set n=P2(r),w=ceil(log2 r),ell=floor(n/(64*2000*17)) and kappa=floor(ell/16000). The plan checks ell>=t and the following integer inequalities:

 n>=16000*152*w,
 4*(128*t*17^2*2000^2)^t*r^150<=ell^t,
 4*C_t*n^t*(128*17*2000)^(2t)*r^150<=(7ell)^(2t).

The first bounds the exponential term in(A2.1) by r^-152<=r^-150/4 for r>=2; the other two bound the polynomial terms by r^-150/4 each. Therefore sampler error is at most(3/4)r^-150. The first comparison is logarithmic: the planner never allocates2^n. All integer bit lengths are polynomial in log r, with t fixed. The guards eventually hold for EVERY integer r: ell is a fixed positive fraction of n up to flooring, r<=n<2r, and t=256>150. Small or failed plans are rejected by the outer total evaluator.


# C1. Canonical fixed-dimension field and cyclic-matrix searches

The current family retains ONE charged tuple serving the FOUR growing-dimension internal maps by C6. Initial dimension10 and B5's prefix/sampler fields and dimension4 helper are instead constructed by the following deterministic canonical searches. This changes the total family and removes their data bits; it assumes no unproved efficient finite-field or primitive-element algorithm. Search costs are explicitly bounded below.

## Least binary irreducible modulus

For degree s>=1 enumerate lower coefficient vectors0,...,2^s-1 in that order, forming monic Q. Test divisibility by EVERY monic binary polynomial of each degree1,...,floor(s/2), using ordinary leading-term XOR division. Accept the first Q with no such divisor. A reducible monic polynomial has a monic proper factor of degree at most s/2 (take a smallest irreducible factor). Conversely any enumerated divisor makes it reducible. Thus this is an exact irreducibility test, including degree1 where no divisor test is needed. A1 proved existence of an irreducible at every degree, so search terminates. At most2^s candidates and fewer than2^(s/2+1) trial factors occur per candidate. Cost is2^(3s/2+O(1))*poly(s), with fixed polynomial bit cost and bounded memory. The least choice is determined only by s; no external table or input data is consulted.

## Least primitive companion action without factoring its period

Fix D>=1 and the canonical field F_q,q=2^s. Enumerate all q^D coefficient tuples(c_0,...,c_(D-1)) in little-endian numeric order, representing monic T^D+sum c_j T^j. Its companion action on a vector v is

 v'_j=c_j*v_(D-1)+(v_(j-1) if j>0 else0).

For each tuple start at v=(1,0,...,0) and iterate this action for period T0=q^D-1 steps. Reject if zero occurs, or if the start returns before step T0. Accept exactly if the start returns at T0 with no earlier rejection. A successful orbit has T0 distinct nonzero states: if another state repeated earlier, determinism would force an earlier return to the start when following the final return. Hence it contains ALL nonzero vectors, and the companion matrix cycles them. No factorization of T0, random primality test or primitive-element oracle is used.

A successful tuple exists. By A1, a degree-D extension field exists and its multiplicative group is cyclic. A generator alpha has degree D over F_q: a smaller-degree minimal polynomial would put it in a field of q^d elements,d<D, contradicting order q^D-1. In the basis1,alpha,...,alpha^(D-1), multiplication by alpha is exactly a companion action of the enumerated form and has the required orbit. Thus the search terminates and provides precisely the A1 cyclic-matrix hypothesis. The first successful tuple is canonical. There are q^D candidates and at most q^D iterations each, each of cost polynomial(D,s). Total cost is q^(2D)*poly(D,s), with polynomial(D,s) memory and no full orbit stored.

The implementation may cache only a fixed bounded number of degree or(s,D) results during execution. This is an ordinary deterministic computation cache, not an external advice source. It constructs a Seeded object with the caller's actual scalar parameters; cached matrices do not change those parameters.

## Why these searches are uniformly polynomial in actual N

For the INITIAL call,dimension D=10 is fixed. Binding occurs only after the exact joint selector checks actuald0<=5L and every ready guard. Its seed contains at least11*s0 bits,even with multiple calls. Thus s0<=5L/11 on every admitted header. Binary-modulus search costs at most2^((15/22)L+O(1))*poly(L);companion search costs at most2^((100/11)L)*poly(L). Since2^L<2N,these are fixed polynomials inN uniformly over all admitted headers. No eventual-good-header assumption is used for this runtime bound.

B5 parameters depend only on r, not the header: w=ceil(log2 r), prefix degree1000w, sampler degree w and helper D=4. Its helper field degree s2 is at most2048w on every admitted call. For an explicit coarse bound, D1=5,000,000w^2+2000w<=2^23*w^2, so log2 M2=ceil_log(D1)<=23+2w for w>=1. H=4ceil(ell^(1/4))<=8r^(1/4),R=80, hence ceil_log(RH)<=10+w. The exact base logarithm is at most32+4*(160w+10)+4*(23+2w)+(10+w)=174+649w<=823w. Its field degree ceil(80*base_log/76)<=2*base_log<=1646w<2048w. All these are scalar inequalities; no giant field is allocated to prove them.

Consequently B5 binary-modulus searches cost2^O(w) and its fixed-dimension helper companion search costs2^(8s2)*poly(s2)=2^O(w). With r>=2,w<=2log2 r, these are fixed polynomials in r. Since r=O(log N) under the initial cap, they are polynomial in actual N. Constants/exponents are large and no practical full-family runtime is claimed.

The internal growing-dimension packages are deliberately RETAINED as actual input bits. Their field degree is Theta(L^2z), so the above exhaustive search would not have the initial polynomial bound there. No canonical-search claim is made for them. Good choices exist by A1 and are fixed only in the charged hard restriction; arbitrary admitted internal data remain bounded ring/matrix computations.

## Exact data saving

For any fixed internal package encoding at identical scalar N,S plans, canonicalization changes no seed length, extraction threshold, r/layout parameter or k. Replacing ONLY the paid initial/B5 packages removes exactly101*s0+1001*w+17*s2 construction bits. C2 separately accounts for internal-package compression. The current complete prefix consists of the H(S)-bit gamma header and the shared internal package from C6. The lower-bound proof uses this current P in4(N-P)-8k. Thus the certified finite expression improves by four times those removed bits at every common admitted scalar plan. The main O(L^9z^5) remainder and protected(4,1) projection are unchanged. This is a stronger exact length accounting, not a new asymptotic exponent claim.


# C2. Auxiliary companion encoding and primitive existence

For any binary field K of size2^s and dimension D, a paid companion package can encode a monic modulus Q in s bits followed by D field coefficients c_j. The package has(D+1)s bits and constructs the companion matrix with subdiagonal ones and last column c_j. This was the preceding encoding; C3 below now reconstructs an isomorphic field and companion action from one binary package instead. The generic companion parser remains for finite identity checks.

Good companion data exist for EVERY s,D, including growing D_int. A1 supplies a degree-D extension and a multiplicative generator alpha. Its minimal polynomial over K has degree D, since a smaller field has too few nonzero elements for its order. Multiplication by alpha in its power basis has the companion form and cycles all nonzero vectors. Any chosen binary basis of K with first vector1 supports the same A1 proof: field operations, root counts and coefficient encoding are unchanged, and binary dot products over all z are all additive characters in that basis.

For arbitrary paid coefficients the companion matrix and bounded ring evaluation still define a total fixed-seed linear map. No internal exhaustive search is performed. The full-matrix encoding has(D^2+1)s bits; companion encoding has(D+1)s bits, a saving(D^2-D)s with scalar seeds, source thresholds and outputs held fixed. Current C3 uses this existence argument and companion action, and its exact loss comparison below defines separate prefixes for each representation. No predecessor verdict is a premise.

# C3. One primitive binary package and polynomial reconstruction

## Actual data and good-package existence

For an internal A1 map keep its scalar SeedParams s,D,m,n,h,error and all call/seed/source parameters. Put t=D*s. C5 decodes t-a(t) actual input bits to a monic constant-one polynomial R of degree t, or to failure. On failure return the zero seeded map before reconstruction. Its fixed small-degree fallback uses C4; t=2 is canonical with zero bits. C3 below reconstructs from the decoded R; its field/ring arguments do not depend on which reciprocal representative was selected.

Thus monic and constant coefficients are fixed to1. Good packages exist: in a field of size2^t choose a multiplicative generator alpha. Its minimal polynomial over F2 has degree t by the same order argument as C2, is monic and has constant1. For t>2 this primitive polynomial is encoded by the C5 unit code when t>=32, or one of its reciprocal pair is encoded by the C4 small-degree fallback; for t=2 the canonical primitive polynomial is encoded. We neither search for nor test primitivity. Such a good R is used only as a PAID input restriction, universal for every required source. All other admitted R are evaluated by the bounded reconstruction below.

## Exact binary linear algebra and the derived subfield

For arbitrary encoded R, work in the t-bit quotient ring A=F2[X]/R. Multiplication is ordinary bounded polynomial reduction. The map x->x^(2^s)-x is binary linear, because repeated squaring is additive in characteristic two. Compute its t columns on1,X,...,X^(t-1) using s repeated squarings per column. XOR Gaussian column elimination retains an original-coordinate combination for each pivot; dependent columns give an independent kernel basis. The kernel basis has length t-rank, spans all solutions and is computed without enumerating ring elements. BinaryColumns.solve reduces a target by pivots and returns its unique input-coordinate combination when the selected columns are independent, or failure if outside their span. These facts follow by the invariant that every stored pivot equals its recorded combination of original columns. Dependent combinations have a new last column index, proving independence and completeness by rank-nullity.

If the fixed-point kernel does not have dimension s, return the ZERO seeded map. Otherwise choose its binary basis deterministically with FIRST vector1: start with1 and append sorted kernel-basis vectors only when independent of the selected span. The unit belongs to the kernel. This process produces exactly s independent vectors spanning it. Call this subring B. It is closed under addition and multiplication: for fixed x,y, Frobenius^(s)(x+y)=x+y and Frobenius^(s)(xy)=xy. It contains0 and1. Encoding its binary coordinates is XOR of the selected vectors; decoding is the exact independent-column solve. Products therefore always decode. The implementation returns0 defensively if an unexpected solve fails, still defining a bounded total algorithm; it cannot fail under the just-proved subring conditions.

For GOOD primitive R, A is the field F_(2^t), and B is its size2^s subfield. To justify its dimension, s divides t; X^(2^s)-X divides X^(2^t)-X, as seen by iterating the s-fold Frobenius or by polynomial division. All its2^s distinct roots lie in F_(2^t), are closed under field operations, and are exactly the fixed-point kernel. Thus the kernel dimension is exactly s, so the first guard passes. The selected basis has field identity coordinate1, as needed by monomial evaluation.

## Relative power basis and companion action

Let alpha be X mod R. Construct the t vectors alpha^j*b_i for0<=j<D,0<=i<s, ordered by j then i, where b_i is the selected B basis. Use the same binary column elimination. If their rank is not t, return the zero seeded map. Otherwise they form a binary basis of A; decoding a vector returns D consecutive s-bit B-coordinates. Decode alpha^D to coefficients c_0,...,c_(D-1) and form the companion matrix from C2.

For GOOD R, alpha has degree D over B. A smaller degree d<D would place it in a field with2^(sd) elements, contradicting its order2^t-1. Hence1,alpha,...,alpha^(D-1) are linearly independent over B, and all t module vectors form a binary basis. The second guard therefore passes. Multiplication by alpha in this relative basis is exactly the derived companion action, and has one orbit on all nonzero vectors. This is PRECISELY the universal A1 matrix hypothesis. The derived B field has2^s elements and an allowed binary basis, so all A1 arguments and seed/error conditions apply unchanged. No field representation or basis cost is uncharged advice: the computation is fully determined by the t-a(t) read bits.

For BAD R, either guard may fail, yielding a zero map, or the program still computes bounded ring/subring/matrix operations. Each fixed seed gives a binary linear map of the coefficient input: all basis choices and matrix/ring data depend only on R, not the payload, and each polynomial coefficient is a payload bit multiplied by a fixed ring value. Zero maps are also linear. No field, primitive or cyclic hypothesis is asserted for arbitrary data, and none is needed for totality.

RelativeSeeded uses the SAME coefficient/seed/monomial algorithm as Seeded, with multiplication and powering delegated to the derived B-coordinate operations. Binary character output uses the specified B basis. The generic binary-polynomial Seeded methods are retained unchanged in meaning for canonical initial/B5 and auxiliary tests. Unit coordinate1 is preserved explicitly. Encoded data are consumed once by Prefix.primitive_seeded; C6 reparameterizes that shared compilation for T,W,E2,E3 without more reads.

## Uniform polynomial runtime

The compiler uses t*s squarings of t-bit ring elements and binary elimination on t columns, plus t relative-basis vectors and another elimination. Encoding, decoding, multiplication and exponentiation have polynomial bit cost in t,s,D. The first nprime monomials are still streamed with polynomial delay. There is NO enumeration of2^t ring elements, no primitive-polynomial search and no factorization. Failed guards are bounded computations followed by zero output.

On an admitted plan, actual package length t-a(t) is counted in P and P+k<=N. C5 proves a(t)<=t/2 for t>=32, so t<=2N; the small-degree fallback is bounded by a constant. Also s,D<=t. All compiler dimensions and loops are therefore polynomial in actual N with fixed exponents, uniformly over every header and data string. Canonical initial/B5 searches retain their separate C1 polynomial bounds. Growing dimension never becomes a runtime exponent. Cached external data or credentials are not used.

## Current prefix

C3 compiles any successfully decoded R without changing scalar SeedParams, seed length, calls, source/output dimensions, k or layout. The current C5 representation has actual length t-a(t), and C5 below supplies all exact comparison prefixes. The compiler always uses only its decoded charged data. No predecessor verdict proves a premise.

# C4. A bounded reciprocal-representative package code

Let t>=3 and n=t-1. The middle coefficients of a monic constant-one binary polynomial form an n-bit word b. Taking the reciprocal polynomial reverses this word. Reciprocal primitive polynomials are primitive, because inverse field generators have the same order. A primitive degree-t polynomial cannot be self-reciprocal for t>2: if alpha^-1 were a conjugate alpha^(2^j),0<=j<t, then its order2^t-1 would divide2^j+1, which is positive and at most2^(t-1)+1<2^t-1. Contradiction. Thus every primitive reciprocal pair has one unique word b smaller numerically than its reversal.

There are(2^n-2^ceil(n/2))/2 such representatives, strictly fewer than2^(n-1). Hence n-1=t-2 bits suffice to index them. The following direct decoder is polynomial in t, not a search through exponentially many polynomials. Partition representatives by their first unequal OUTER pair j=0,...,floor(n/2)-1. The j earlier outer pairs are equal, with j free bits. At pair j the more significant bit is0 and its partner is1. The n-2j-2 inner bits are arbitrary. This block has2^(n-j-2) elements.

Given the t-2-bit index, subtract these block sizes in increasing j until it fits one. Within that block use the low n-2j-2 bits for the inner word and the remaining j bits for the equal outer pairs. Set the unequal pair to0,1 and assemble b, then R=X^t+1+X*b. If no block fits, return failure. Every canonical nonpalindromic word appears exactly once, by its unique first unequal pair and unique remaining bits. The number of decoded words agrees with the above count by summing the block sizes. All operations use O(t)-bit integers and O(t) iterations, so runtime and memory are polynomial in t. No primitive test or field search is performed by the decoder.

For t=2, the primitive polynomial X^2+X+1 is self-reciprocal; use it canonically with t-2=0 data bits. Current internal dimensions have t=D_int*s>=2. Thus good packages exist at every relevant scalar t, including this exceptional case. Failed indices give the zero seeded map. Successful indices use the C3 bounded compiler on the decoded R; primitive good representatives satisfy all its hypotheses.

C4 saves exactly ONE further actual bit per internal package versus C3's unrestricted middle-bit encoding, hence four prefix bits and16 units in the exact gate-loss expression. The reciprocal-only package has t-2 bits, so its saving versus companion encoding is s+2. All fields, seeds, calls, entropy thresholds, k and layouts are unchanged. This C4 encoding is the fixed small-degree fallback of C5. The current large-degree format and exact prefixes are supplied in C5 below. The9,5 asymptotic claim and protected(4,1) projection are unchanged. This is a finite encoding improvement only.

In the current whole-row revision, independence h and precision e are separate. Internal A1 field degrees use e; the unchanged encodings C1--C4 still have prefix O(L^2z^2). The new robust analysis closes k=O(L^9z^5). Their exact encoding comparisons are only at identical scalar plans and do not claim unchanged k between different scale revisions.


# C5. A short small-factor unit package with a polynomial decoder

## Deterministic savings and the small-factor modulus

Put t=D*s. For t<32 retain the C4 decoder with a(t)=min(2,t); degree1 may use the canonical X+1. Every actual internal scalar parameter has t>=32: its positive per-call precision gives base_log>=32+4>=36, hence s>=36 and D_int>=2 yields t>=72. This fallback also defines all small compiler fixtures. For t>=32 put

 k=floor(log2 t)-2, a(t)=floor(log2(k+1)).

Then k>=3,a(t)>=2. Form the set F_k of ALL monic irreducible binary polynomials of degrees1,...,k, and their squarefree product M. Write d=deg M. For each degree j, the j distinct roots of each monic degree-j irreducible lie in F_(2^j); roots of different irreducibles are disjoint. Thus j*I_j<=2^j, where I_j counts those irreducibles, and

 d=sum_(f in F_k)deg f<=sum_(j=1)^k 2^j<2^(k+1)<=t/2.

The splitting/root assertion follows from the Frobenius field argument in A1. All degree1 irreducibles X and X+1 are included. Define

 Phi=product_(f in F_k)(2^(deg f)-1).

No construction or large Phi integer is needed by the scalar plan: it uses only t.bit_length and a(t). The modulus/decoder are built AFTER all actual ready guards pass.

## An elementary exact density bound

The ratio Phi/2^d is product_(f in F_k)(1-2^(-deg f)). Expand the reciprocal as the finite product of convergent positive geometric series. By unique factorization of monic binary polynomials it is the sum of2^(-degree H) over all monic H whose irreducible factors belong to F_k. This includes EVERY monic polynomial of each degree0,...,k. There are exactly2^j such polynomials at degree j, so each of these k+1 degrees contributes1. Hence

 product_(f in F_k)(1-2^(-deg f))^(-1)>=k+1,
 Phi<=2^d/(k+1)<=2^(d-a(t)).                            (C5.1)

This is a positive counting argument, not an analytic Mertens theorem, asymptotic prime-distribution assertion or unsupported field algorithm. It supplies the exact fixed-length residue budget at every t>=32.

## Exact decoder and complete good-package coverage

Read t-a(t) actual bits. After M is constructed, split them into t-d low bits q and d-a(t) high bits I. Let Q=X^(t-d)+q(X), a monic quotient. Decode I in mixed radices b_f=2^(deg f)-1, in the fixed enumerated factor order: residue r_f=(I mod b_f)+1, then I=floor(I/b_f). Radix1 for degree1 causes no problem and forces residue1. If a nonzero I remains after all factors, return failure, hence a zero seeded map. Otherwise the residues are a unique tuple of nonzero classes, and every such tuple occurs at an index less than Phi. By(C5.1) all these indices fit the d-a(t) allocated bits.

For each f let M_f=M/f and u_f be the inverse of M_f modulo f. Distinct monic irreducibles are coprime, so this inverse exists. The implementation computes it as(M_f mod f)^(2^(deg f)-2) in the degree-deg f field; no factorization is needed. Define

 H=sum_(f in F_k) M_f * ((r_f*u_f) mod f).

Every summand has degree below d, and reduction modulo f leaves exactly r_f: all other M_g terms contain f. Thus H is the UNIQUE polynomial of degree below d with those residues. Uniqueness follows because the difference of two solutions is divisible by squarefree M and has degree below d. Return

 R=Q*M+H.

It is monic of degree t. Since X and X+1 are factors of M and their residues are1, R has constant1 and R(1)=1. All successful outputs are coprime to M, but no primitivity or irreducibility is claimed for arbitrary codes. The bounded C3 compiler handles them, or returns a zero map at a rank guard.

EVERY primitive binary polynomial R of degree t>=32 is represented. Such R is irreducible of degree t>k, hence coprime to M. Its ordinary division by M has a monic quotient of degree t-d and a unit remainder; the unique nonzero residue tuple has a mixed-radix index less than Phi<=2^(d-a(t)). Both pieces fit the package exactly. Good primitive packages therefore exist by A1/C3 and supply the same universal size2^s field/cyclic-matrix hypotheses for every required source. No search for a primitive polynomial occurs. The small-degree fallback retains its previously proved good-package coverage.

## Bounded algorithms on every admitted header and data string

Enumerate all2^j monic binary polynomials of each degree j<=k. A candidate is irreducible exactly when no previously found irreducible of degree at most j/2 divides it. A reducible monic has such a factor; that degree is smaller, so it has already been enumerated. This supplies the complete list with no field-size oracle. There are O(2^k) candidates and at most O(2^(k/2)) trial factors per candidate, with polynomial(k) division cost. Since2^k<=t/4, this is polynomial in t. Multiplying factors, exact polynomial division, field powering at degrees<=k, mixed-radix arithmetic and the displayed CRT formula also have fixed polynomial bit cost in t. All integers/polynomials have O(t) bits apart from polynomially many temporary arrays. A fixed bounded computation cache keyed only by k is deterministic reuse, not external advice.

The scalar package length t-a(t) is computed by bit lengths ONLY; no field/list/product of length t is allocated during scalar planning. At binding, for t>=32 we have a(t)<=t/2, so t-a(t)>=t/2. The actual prefix guard P+k_source<=N therefore gives t<=2N. Thus sieve/CRT construction and C3 compilation are uniformly polynomial in actual N, including invalid mixed-radix codes and bad R. Degrees below32 cost a fixed constant and use the C4 bounded decoder. This proves totality on all admitted data, not only on good restrictions. Zero maps and compiled maps remain fixed-seed binary linear because their entire representation depends only on the paid prefix.

## Exact saving within the unshared representation

The unshared C5 format has one package of t_i-a(t_i) bits per map. Its exact comparisons with reciprocal/companion/full encodings follow by subtracting their lengths at identical scalar plans. Since a(t_i)>=2 and tends to infinity at a fixed good header, it never uses more bits than the reciprocal format and eventually saves unboundedly many. C6 below changes the CURRENT input layout further by sharing one common package; its current exact prefix is defined there. No comparison equates k across distinct scale revisions or uses a predecessor verdict as a premise.

# C6. One paid common field and cyclic action for all internal maps

## Common scalar field with exact guards

The four maps retain common D_int and precision e. W/T threshold is2944ez^2ell; E2/E3 threshold is9096ez^2P2(M). Compute each full-output plan's ownnprime,K,H,m0,error,calls and keptM. Raise only field degree to the common maximum s_star, then recheck exact entropy/advice and ell-seed predicates. No discarded output suffix becomes a transcript or seed. One paid C5 package serves the common dimension and field.

Current C32 supplies H_WT<=301z/100,H_rows<=93z/10,R<=573z/100,s_star<5334e/1000. W/T advice459985491/312500Reqez^2<K/2=1472Reqez^2; row advice142121763/31250Reqez^2<K/2=4548Reqez^2. Each positive gap pays the exact4rho_log+10 term with ownK. Every first full candidate passes ONEcall and eventual actualseed<8ze is paid byell>=8ze. All finite actual guards remain checked; common-degree sharing does not replace the smaller per-map H by the largest row H.

## A single good tuple is universal across the four calls

Read ONE C5 package of length

 t_star-a(t_star), t_star=D*s_star.

Use C3 to compile its decoded R once. For a primitive good R the compiler produces a size2^s_star field with a specified unit-first binary basis and a D-by-D cyclic matrix A. That SAME field and A satisfy A1 for every one of the four already checked parameter tuples,including all conditional sources and repeated calls. A1 is a theorem for ANY good field/cyclic matrix satisfying its scalar conditions; it is not a probabilistic existence statement choosing a different matrix for each source. There is therefore no union-bound or independent-construction-data premise to pay.

B3--B4 use strong extraction on arbitrary conditional affine sources with fresh independent SEED variables. They already allow arbitrary correlations among nondesignated maps/streams. Sharing their fixed field/matrix constants changes no seed independence. Every map retains its own source/output metadata and full-seed error guarantee. The source/error/rank transcript proof is therefore unchanged after the explicit scalar rechecks above.

The implementation compiles one template using s_star,D (other template SeedParams are unused by compilation). Each bound internal map is a fresh RelativeSeeded object with its OWN nprime,H,m0,error parameters but the SAME immutable-by-use field,ring and matrix objects. Reparameterization checks the dimension and field degree. If the shared package fails a decoder/compiler guard,each map receives a ZeroSeeded object with its OWN output parameters. No internal map reads additional prefix data. Generic unshared encodings remain only as auxiliary component paths.

## Exact cursor, runtime and saving

After the H(S)-bit gamma header,the input cursor consumes the single common package once. Initial A1,short E1,B5 and all four internal wrappers then consume zero further data bits. The engine checks the actual cursor and removes the ENTIRE prefix. Current

 P=H(S)+t_star-a(t_star).

As before,a(t_star)<=t_star/2 at actual large degree,so the ready prefix guard implies t_star<=2N. The one sieve/CRT/kernel/module compilation is polynomial in N. Reparameterizing four metadata objects is polynomial and does not copy/compute an exponential table. Bad tuples still give bounded ring or zero-map evaluation; good tuples exist at every required degree and are paid restrictions. All remaining canonical/search/seed/row runtime bounds are unchanged. No source-dependent data are constructed.

For an EXACT comparison retain each individual pre-common package length b_i=D*s_i-a(D*s_i). The common package length is exactly b_j for any index j attaining max s_i: its format depends only on D and that field degree. At any common admitted scalar plan define

 P_individual=H(S)+sum_i b_i.

Then P_individual-P=sum_(i!=j)b_i>0. Thus the certified loss expression4P_individual+8k decreases by FOUR times the sum of the other three package lengths. k,initial d0 and layout are identical; actual internal seed slices may grow but are revalidated. This is a comparison of proven resource expressions,not an assertion about equality of two Boolean functions or their true circuit complexities.

At any fixed good scale,e dominates every other logarithm in the field formula. Each individual s_i/e equals4R/(R-D)+o(1),which tends to16/3 because R=4D+32 eventually andD tends to infinity (currentD=max(2,ceil(10z/7))). The maximum has the same limit. Since a(t)=O(loglog t) and L=o(D*e),current P divided by P_individual tends to1/4. In particular current P=O_S(L^2z^2)=o(k),and the whole-row9,5 remainder and protected(4,1) projection stay unchanged. The exact prefix reduction and single compilation are new relative to71e6d94e; no end-to-end practical runtime claim is made.

# C7. Self-delimiting paid scale header

The CURRENT program reads the standard Elias gamma code of positive S. If j=floor(log2 S), the code is j zeros followed by the full (j+1)-bit MSB-first binary S, and H(S)=2j+1. The first1 determines j and hence the exact end of the header; appending arbitrary package/payload bits never changes a completed parse. The parser rejects more than L leading zeros, missing first1, a truncated binary suffix or S>2^L. It reads at most2L+1 bits and uses O(L)-bit arithmetic. Inputs below128 return0 before parsing. Every scale in1..2^L has exactly one admitted code whenever enough input bits exist. All other prefixes simply define value0 in the same fixed total language.

PolyPlan computes H(S) from the parsed positive scale. All extraction/scalar/layout parameters depend on the SAME N,S as before; only the actual header length and hence prefix/payload cursor change. The engine consumes the entire H(S)-bit code and then the one C6 package. Its actual prefix is P_gamma=H(S)+t_star-a(t_star). Every readiness cap uses this exact P_gamma. There is no scan of the paid package for a delimiter and no external scale advice.

For the one fixed sufficiently large good integer S0 from B7, H(S0) is a fixed constant. Once2^L>=S0 and N>=H(S0), its code fits at EVERY integer length. Fix those actual header bits and one good shared package. The N-P_gamma free coordinates, zero padding, affine disperser and finite gate argument apply unchanged. Thus C_B2(F_N)>=4N-4P_gamma-8k for the CURRENT gamma-coded family. This is one fixed language for all lengths, not a promise or subsequence.

At identical admitted scalar plans the preceding fixed-width representation paid P_fixed=L+t_star-a(t_star). Consequently the difference between its certified loss expression and the current one is exactly4*(L-H(S)); the comparison is strictly beneficial for every fixed S once L>H(S), and tends to infinity at S0. It may have the other sign at arbitrary growing scales, so no pointwise improvement over all headers is claimed. C6's four-versus-one package comparison remains valid with H(S) in both expressions and still has ratio1/4 at fixed S0. Current P remains O_S(L^2z^2)=o(k); the9,5 asymptotic theorem and(4,1) projection are unchanged. Earlier representation comparisons are resource comparisons at the specified same header convention, not equations for actual circuit complexity across distinct languages. This gamma-header change and its full charged accounting are new relative to77e43640-f90b-45ac-a18c-0233345a54c6.

# C11. Retained exact source-budget method

C11 source6b4b7dde-5bf3-4091-b601-16bd48e0af8d following84bbd016 supplied exact integer source-budget accounting;the current C17 d maximum is a later change. Its historical single-threshold X budget had form k=hY+100B_cap*h*a*ell+10d+C_old*e*z^2*ell;the current exact budget uses C_W at that position. C12 retains that method,changes the alphabet and removes Y's outer rounding safely via an explicitly dyadic initial output and row prefix. Its complete current proof follows. No prior judgment is a premise.

# C12. Dyadic initial output prefix and a larger alphabet

Let (a,B_cap,E0) be the exact joint selector in proof.md C32,s_idx=log2a. Reject a>h. All current resource envelopes are B7.1,with raw integer Y=4B_cap*h*d+C_Yez^2P2(d) and exact X budget k. Compute k as a bounded scalar integer,check ceil_log(k)<=floor(log2 N),and store that exact k. Every source/seed/prefix guard remains. The conservative finite cap passes eventually for fixed-scale polylogarithmic k.

Set M_init=P2(Y). Request M_init bits from the canonical initial dimension10 AffineSeedPlan. The engine EXPLICITLY takes only its first Y bits before any index-round source operation. E2/E3 receive this Y-bit prefix with their correct ambient metadata. No extra M_init-Y bits are revealed in the B3/B4 source/rank accounting. A full M_init-bit uniform row has an exactly uniform Y-bit prefix;bad-seed probability for the full row therefore bounds the bad-prefix probability. Fixed-seed full rank implies prefix rank. Thus initial deficient-row fraction<=2*2^-E0 remains. On arbitrary finite headers the actual plan guards can reject a request too large for k;this defines value0 in the same total language.

BecauseY<=M_init<2Y andk~hY withh tending to infinity,k/4>=M_init eventually. Joint selection derivesE0>=L/125 and enforcesE0<=L;generic B1 has nprime=O(k),so the D10 full-candidate advice has strict L exponent89/10<9. ONE initial call eventually suffices uniformly across selected precisions. All finite plans already check exact validity,5L seed cap and the explicit robust error margin. The Y-bit truncation is deterministic polynomial work.

Internal arbitrary requested outputs now request one full dyadicP2(M) and explicitly return only firstM bits. B2.2/C6 prove separate W/T and row source/ceiling/advice bounds with the same shared field andactualseed<8ze paid byell>=8ze. All X/Y/Q residuals and B3/B4 row-error bounds apply to returned M lengths. Current all-length source/initial/robust/prefix/runtime closure is supplied in proof.md C32;no selectedd=q equality is needed.

The earlier C12 comparison concerned its original h/e convention. C14 below supplies its historical matched-bias comparison;the current comparator is in proof.md C32;no unchanged-h/e claim is made across the later allocations.

Current k remains Theta_S(L^9z^5/log z),so4(N-P)-8k yields the same full logarithmic-divisor theorem and every same-family weakening at ALL sufficiently large actual integer lengths. The native projection stays(4,1). C12's new contributions are the full-output/prefix logic,larger alphabet and strict leading-constant comparison;all required extraction,robustness and gate arguments remain supplied without an accepted-result premise.

# C13. Retained resource-decoupling motivation

C13 source16909354-df8a-4331-94e0-31ae44bce9c9 used square/fourth-root row/error scales from one header. C14 replaces that coupling by independent paid components. The current good-parameter proof and comparison are below;the root-scale formula is used only for the explicit predecessor comparator,not as a definition of current inputs.

# C14. One paid code for independent row,error and bias scales

The ordinary positive gamma header code remains S,with S<=2^L. Write S-1 in binary and send its bit at position3j to bitj of U-1,position3j+1 to bitj of V-1,and position3j+2 to bitj of T-1. Thus all three decoded parameters are positive. The inverse interleaves the three finite binary words and adds1. Every bit has exactly one destination and inverse,proving a bijection between positive integer codes and positive integer triples. Zero high bits cause no ambiguity. Each component is at most S,so every scalar bit length remains polynomial in actual N. Decoding scans O(log S) bits and no invalid negative word is shifted indefinitely. The code is range-checked before decoding.

Current h=P2(U L^2),e=P2(V L^2z). The layout is called with scale T^2,so its integer square root is EXACTLY T and its seed bias is floor(T*ceil(log2r)/(30*125^2)). The header length is still H(S)=2floor(log2S)+1;no separate parameter bits are omitted. All other C12 exact source/seed/initial-prefix formulas are retained. Any fixed triple has one fixed finite code,which fits S<=2^L and the actual header at every sufficiently large integer length. The code is a normal paid restriction,not advice or a promise. Arbitrary codes,malformed inputs,bad packages and failed guards remain part of the one total language.

For any fixed positive triple,all L/z resource orders are unchanged. The alphabet guard passes eventually since a=O(z^3)=o(U L^2). q0=Theta(L^4z^5),eprime=O(L^4z^5),d=Theta(L^4z^5) by the C32 joint-selection proof,Y~4B_cap*h*d,k~hY. Internal largest K=O_{U,V,T}(L^6z^8),field/advice/seed bounds,common-degree revalidation,single full dyadic internal calls and dyadic initial M_init with explicit Y prefix retain their supplied proofs. Initial advice is o(k),logM_init=o(L),and the exact joint guard givesd0<=5L with the required robust margin. Canonical data/search caps and fixed polynomial runtime are retained. Components<=2N and T^2 has O(log N) bits;normal-plan dimension/source/prefix caps bound all actual work. Large bias causing an invalid layout is rejected by bounded scalar checks before evaluator binding.

Let U0 be a sufficiently large power of two for the row constants,and K0 a sufficiently large power of two for the precision constants in B7,also takingK0>=U0. Choose T sufficiently large for B5/bias. Then the fixed triple(U0,K0*T,T) meets every robust hypothesis,all explicit guards eventually,and the strict bad-row/coverage margin. Choose one such triple and its interleaved S0. The complete gate/disperser proof yields4(N-P)-8k for actual length at EVERY sufficiently large integer N,with P=H(S0)+shared_data=O_{S0}(L^2z^2)=o(k). This proves the same full logarithmic-divisor theorem and same-family weakenings unconditionally.

The C14 packed-code allocation is retained,while proof.md C32 supplies the CURRENT comparison and explicitly rechecks the weakened Y reserve in the complete B4 proof. No historical matched-bias equality is asserted as a current comparator.

# Current resource accounting

The complete C15/C16/C17 source-reserve,weak-seed precision and exact Q/seed-budget arguments are supplied in proof.md. The full initial-prefix and packed-header arguments above remain part of the current proof.
