# Supplied extraction arguments for the polylogarithmic branch

All dimensions and entropies are binary. Construction data for A1 are actual
input-prefix bits. The deterministic fields in B1 require no such data.
Every use below fixes one tuple of good prefix data for the hardness proof,
while the evaluator remains total on all prefix strings.

# A1. a strong linear seeded extractor

## Exact construction and parameters

Let the source length be n bits. Fix positive integers D,h,m,s with
binomial(h+D-1,D)>=n, and let q=2^s. Let Q be any monic irreducible binary
polynomial of degree s and K=F_2[T]/Q. Let A be a D-by-D matrix over K
whose action cycles through K^D minus zero. Such matrices exist, as
justified below. This property is a mathematical hypothesis on A; the
finite evaluator does not search for or test A.

Enumerate all D-variate monomials of total degree below h, first by total
degree and then lexicographically. For x in F_2^n, let P_x be the K-valued
polynomial with its first n coefficients x_1,...,x_n and all remaining
coefficients zero. This is an injective F_2-linear encoding with degree
at most h-1. The seed consists of v in K^D and z in F_2^s. Output

    E(x;v,z)_i = <bits(P_x(A^i v)),z>,   1<=i<=m.

The bracket is the binary dot product. For each seed this is linear in
x. Evaluation takes polynomial time in n,m,D,h,s and the number of
monomials, without enumerating K^D or finding a primitive element.

Take 0<epsilon<=1/4 with rho=epsilon/(8m) reciprocal to an integer. Let
r=16D+16. Sufficient deliberately loose hypotheses are

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

The circuit-family wrapper reads Q and A as actual input bits,
with their lengths predetermined by the total length. For bad data it
still runs the same bounded ring/matrix algorithm, defining a total
polynomial-time function. A lower-bound restriction fixes good data,
whose existence was proved; all data bits must be paid for in actual
input length. The matrix and field data require (D^2+1)s bits. Every one of those bits is charged in proof.md.




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

The implementation uses sigma=1/2000,b=17,t=256 and4t=1024. With n=2^w it reads1024 independent field coefficients and evaluates the degree-at-most1023 polynomial at distinct field points0,...,ell-1. Its field values are the selected positions. For any at most1024 distinct arguments, the evaluation map is surjective by the Vandermonde determinant: the product of pairwise argument differences is nonzero. Hence those positions are independent uniform field elements. If fewer arguments exist the entire smaller tuple is independent. A good monic irreducible modulus of degree w is actual charged construction data; existence is supplied in A1. Even bad monic data give bounded w-bit ring values and a total polynomial-time sampling algorithm.

In B5 set n=P2(r),w=ceil(log2 r),ell=floor(n/(64*2000*17)) and kappa=floor(ell/16000). The plan checks ell>=t and the following integer inequalities:

 n>=16000*152*w,
 4*(128*t*17^2*2000^2)^t*r^150<=ell^t,
 4*C_t*n^t*(128*17*2000)^(2t)*r^150<=(7ell)^(2t).

The first bounds the exponential term in(A2.1) by r^-152<=r^-150/4 for r>=2; the other two bound the polynomial terms by r^-150/4 each. Therefore sampler error is at most(3/4)r^-150. The first comparison is logarithmic: the planner never allocates2^n. All integer bit lengths are polynomial in log r, with t fixed. The guards eventually hold for EVERY integer r: ell is a fixed positive fraction of n up to flooring, r<=n<2r, and t=256>150. Small or failed plans are rejected by the outer total evaluator.
