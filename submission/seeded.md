# Linear seeded components, with all material arguments

This manuscript supplies a folded affine condenser and a
basic algebraic short-seed extractor, followed by the linear code-bit
extractor needed inside alternating extraction. All fields have
characteristic two. Logarithms of integer sizes are binary.

The folded evaluation is the folded-Wronskian rank-condenser idea of
Forbes--Shpilka and Guruswami--Kopparty (see Forbes--Guruswami, ECCC
TR14-162, Construction 2.7). The algebraic orbit construction adapts
Shaltiel--Umans, JACM 52(2), 2005, Sections 4.2--4.5. The following
arguments are supplied independently with the actual conservative
parameters. In particular, the optimized small-error construction and
algorithmic list decoder from that literature are not premises here.
The reconstruction below has no runtime requirement; only the extractor
itself must run in uniform polynomial time.

## 1. Deterministic finite fields and their multiplicative generators

For any finite base field F of size q and integer s>=1, an irreducible
monic polynomial of degree s can be found by exhaustive search. An exact
test for a candidate P is

    X^(q^s)=X modulo P,
    gcd(P,X^(q^(s/p))-X)=1 for each prime divisor p of s.

Indeed, X^(q^s)-X has no repeated root (its derivative is -1), and its
roots in a splitting field form the field of q^s elements. Their
Frobenius orbits give exactly the monic irreducibles with degree dividing
s. The first test therefore makes P a product of distinct such factors;
every proper divisor t of s divides s/p for some prime p|s, so the gcd
tests exclude all proper factors. Conversely an irreducible of degree s
passes both tests.

Existence follows without advice. If I_t counts degree-t irreducibles,
then q^s=sum_(t|s)t I_t and t I_t<=q^t. For s>=2 the contribution of
proper divisors is at most sum_(t=1)^floor(s/2)q^t<q^s. Degree one is
immediate. Thus I_s>0. The exhaustive search has at most q^s candidates;
Frobenius and gcd arithmetic take polynomially many bit operations in
s and log q per candidate. Bounds on q^s below establish the required
runtime in the ORIGINAL input length, even when the encoded source has
a smaller length.

The multiplicative group of a finite field is cyclic. To see this, let e
be its exponent. In a finite abelian group there is an element of order
e: combine elements realizing the maximal power of each prime in e.
Every group element is a root of X^e-1, so the polynomial root bound
makes the group size at most e. Lagrange's theorem gives the reverse
inequality. Thus its exponent equals its size and a generator exists.
Trial division factors q^s-1 in at most O(sqrt(q^s)) iterations; enumerate
field elements and test g^((q^s-1)/p)!=1 for each prime divisor p to
find the first generator. This is deterministic and polynomial in q^s.
All searches and encodings use fixed increasing coefficient orders.

## 2. Folded evaluation preserves affine entropy

Let the original bit length be N and let 1<=k<=N. Put

    q_1=2^ceil(log_2(N^16*k)), w_1=log_2 q_1.

For small exceptional arguments use any stated guard in the main family;
here N>=2. Construct F_(q_1) by Section 1, and choose a generator gamma.
Pack all N input bits injectively into coefficients of a polynomial f
of degree less than M=ceil(N/w_1); the final coefficient is zero padded.
For a seed a in F_(q_1), let

    Cond_a(f)=(f(a),f(gamma*a),...,f(gamma^(k-1)*a)).

The output has n'=k*w_1 bits. Each fixed-seed map is F_2-linear.

Take any affine input source of dimension at least k, and choose a
k-dimensional subspace V of its direction space. Regard its elements
as polynomials over F_(q_1), and let W be their span over that larger
field. Its dimension r is at most k. Row reduction gives a basis
f_1,...,f_r for W with distinct leading degrees d_1,...,d_r, all below M.
The determinant

    D(a)=det[f_i(gamma^j*a)]_(0<=j<r,1<=i<=r)

is a nonzero polynomial. Its leading coefficient is the product of the
basis leading coefficients times the Vandermonde determinant on the
r distinct elements gamma^(d_i). They are distinct because gamma has
order q_1-1>M. Its degree is sum d_i<=r(M-1)<=kN. Consequently at most
kN seeds make D(a)=0.

For every other seed, the first r output field coordinates are injective
on W, hence on V. The entire F_2-linear output has rank at least k on
the source direction space. An affine shift changes only its translation,
not its entropy. Thus for all but a fraction

    kN/q_1 <= N^-15

of seeds, the condensed source is affine with entropy at least k. It is
not asserted to be uniform on every output coordinate. This distinction
is important: a rank condenser alone does not provide the uniform seed
needed in alternating extraction.

Computing k evaluations uses polynomial time in N. Field searches cost
polynomial time in q_1, and q_1<=2N^17. Output length n'=O(k log N).
No field description or generator is an input coordinate or advice.

## 3. A strong algebraic extractor on the condensed source

In this section N remains the ORIGINAL input length. Assume

    L=floor(log_2 N), J=ceil(log_2(L+2)),
    k=floor(2^10000*L^25/J), k<=N,
    n'=k*w_1, m=floor(k/(2^300*L^2)),
    h=floor(k^(1/32)), d=64,
    q=2^ceil(log_2(N^128)), w=log_2 q.

All statements are for sufficiently large integers N. Integer-root
rounding changes these powers by factors tending to one. In particular

    binom(h+d-1,d)*w >= n'

holds eventually, since its left side has order k^2 log N whereas n'
is O(k log N). Pack the n' input bits, with a zero-padded final field
coefficient, into the first degree-lexicographic monomials of total
degree less than h in d variables over F_q. Denote the resulting
polynomial by P_x. This coefficient map is F_2-linear and injective.

Construct the degree-d extension of F_q and a multiplicative generator.
In its fixed power basis let A be the matrix of multiplication by that
generator. It is invertible, has order q^d-1, and its powers carry every
nonzero vector through all nonzero vectors of F_q^d.

For uniform seed v in F_q^d and uniform j in F_2^w, output m bits

    Ext(x;v,j)_i = <P_x(A^i v),j>, 1<=i<=m,

where the brackets are the ordinary binary coordinate inner product in
the chosen basis of F_q. Its seed length is (d+1)w=65w. Every fixed-seed
map is F_2-linear in x. Computing the polynomial at m points, constructing
the fields, and finding A all take N^O(1) time: q^d<=2^64 N^8192 and d
is the fixed integer 64. This runtime is measured in N, not n'.

We prove that this is a strong extractor of error at most N^-8 for every
source on n' bits of min-entropy at least k. This general statement is
used only on the affine condensed sources. The proof does not assume
that their larger-field span is an affine space over F_q.

### 3.1 Two elementary counting estimates

First, if a collection of distinct univariate polynomials of degree at
most D each agrees with at least B pairs in a set of P distinct pairs
(t,z), then its size is at most 2P/B when B^2>=2DP. Select exactly B
agreement pairs per polynomial. Two selected sets intersect in at most
D pairs, by the polynomial root bound. If there are L polynomials and
r_u is the selected incidence at pair u, Cauchy--Schwarz gives

    (LB)^2 <= P*sum_u r_u^2 <= P*(LB+L(L-1)D).

Rearranging gives L(B^2-PD)<=P(B-D)<=PB, and hence the claim. This
is a list SIZE bound only. Exhaustive enumeration suffices in the proof;
no polynomial-time decoding algorithm is being assumed.

Second, fix a positive integer r. For 2r-wise independent Bernoulli
variables with common mean p, let Z be their sum and mu its mean. In
expanding E(Z-mu)^(2r), every term containing an index just once vanishes.
The remaining terms have at most r distinct indices. For each power
at least two, the absolute centered moment of a Bernoulli variable is
at most 2p. Counting the possible index patterns therefore gives, for
mu>=1,

    E(Z-mu)^(2r) <= (2r)^(2r)*sum_(s=1)^r(2mu)^s
                 <= r*(2r)^(2r)*(2mu)^r.

Markov's inequality bounds Pr[Z<mu/2] by K_r*mu^-r, where K_r is a
fixed finite constant depending only on r. This coarse bound is ample
for the parameters below; no asymptotically sharp concentration result
is needed.

### 3.2 Reconstructing a fixed input from a successful q-ary predictor

Put rho=N^-8/(2m), and fix an index i* in [m]. Consider an arbitrary
predictor which, given v and the first i*-1 field outputs, supplies a
list of at most rho^-2 candidate next field values. Lists may be rounded
down to an integer size. The predictor is allowed to depend on v; this
is what will give STRONG extraction. Fix an input x for which it predicts
correctly on at least a rho/2 fraction of v.

Set r=128 and choose 2r distinct positions t_1,...,t_(2r) uniformly in
F_q and independent uniform vectors y_1,...,y_(2r) in F_q^d. Let p_1
and p_2 be the unique vector-valued curves of degree at most 2r-1 with

    p_1(t_i)=y_i for all i;
    p_2(t_i)=A*y_i for i<=r, and p_2(t_i)=y_i for i>r.

Define P_(2s+1)=A^s p_1 and P_(2s+2)=A^s p_2. Each curve has 2r-wise
independent uniform values at distinct field positions. Adjacent curves
agree on r uniformly random distinct positions whose set is independent
of the CURRENT curve. For example, conditional on all t_i the 2r values
of either curve are independent uniform vectors because A is invertible;
thus its polynomial coefficients are uniform and independent of the
chosen t_i. For an odd-current curve the shared positions are the first
r anchors, and for an even-current curve they are the last r anchors.

The polynomial P_x composed with any curve has degree at most

    D=(2r-1)(h-1).

An advice string gives the coefficients of these univariate polynomials
for the first 2i* curves, INCLUDING when i*=1. Its length is at most

    a_adv=2m*(D+1)*w.

This initialization supplies known values on a preceding curve even for
the first predicted output. In particular, the case i*=1 receives two
initial curve compositions; it is not a zero-advice reconstruction of
multiple possible constant polynomials. The uniformly charged 2m-curve
budget already covers all these compositions for every i*<=m.

Starting at curve number 2i*+1, reconstruct successive curves. To predict P_x(u)
for u on a new curve, apply the fixed predictor at seed A^(-i*)u with
preceding values P_x(A^(-i*+1)u),...,P_x(A^-1u). These are known from
previous curves because A^-j P_i=P_(i-2j). At the first reconstruction
step the smallest required preceding index is
(2i*+1)-2(i*-1)=3, within the supplied initial range when i*>1;
when i*=1 there are no preceding predictor inputs, while curve 2 supplies
the known shared anchors. Every later step has both its predictor inputs
and the immediately preceding anchor curve already reconstructed. This remains true when the
predictor depends on its seed. No value of x itself is encoded in the
predictor beyond the separately counted advice.

The indicator of a correct list prediction on a uniform point has mean
at least rho/2. On the current curve its q indicators are 2r-wise
independent. Section 3.1 shows that, except with probability at most
K_r*(rho*q/2)^-r, at least rho*q/4 positions give a correct candidate.
There are at most q*rho^-2 distinct predicted pairs. With
B=floor(rho*q/4)>=rho*q/8, the list-size estimate applies eventually,
because q>>D*rho^-4. It leaves at most 16rho^-3 univariate polynomials
of degree at most D that could be the true composition.

The candidate list depends only on the current curve and the fixed
predictor and x. A wrong candidate agrees with the true composition at
at most D positions. The r independent-of-current-curve shared anchor
positions identify the correct one, except with probability at most

    16rho^-3*(D/(q-r))^r <=16rho^-3*(2D/q)^r.

These probabilities need not be conditioned on success of earlier
reconstruction steps: bound each failure event using the TRUE preceding
values, and take a union bound over all steps. On their common complement
the actual reconstruction agrees with those values inductively.

Both per-step failure probabilities are eventually smaller than
1/(16q^d). Indeed k<=N gives h<=N^(1/32), m<=N, rho>=N^-9/2,
and N^128<=q<2N^128. With fixed r=128, the first error is at most a
fixed constant times N^(-119r); the second is at most a fixed constant
times N^(27-(4095/32)*r). Both exponents are far below -128d=-8192.
Thus a union bound over 2q^d steps costs at most 1/8 eventually.

With probability 1-q^-d, y_1 is nonzero. On that event, the orbit of
that point through the curves visits every nonzero vector in F_q^d.
Combining this with the reconstruction success event has probability
at least 1-1/8-q^-d>1/2 eventually. The nonzero evaluations uniquely
determine P_x. For completeness, a nonzero d-variate polynomial of total
degree at most h-1 has at most (h-1)q^(d-1) zeros: induction on d,
expanding in the last variable, proves this familiar root estimate.
Since h<<q, it cannot vanish at all q^d-1 nonzero vectors. Thus no value
at zero is needed. Finally the injective coefficient map determines x.
All reconstruction searches may be exhaustive; they are counting tools,
not subroutines of the extractor evaluator.

### 3.3 Entropy count and conversion to binary output

Fix the predictor and fix reconstruction randomness for which many good
inputs are reconstructed. The advice has at most 2^a_adv possible values.
The preceding success probability at least 1/2 implies that a source of
min-entropy k gives mass at most 2^(a_adv-k+1) to the inputs predictable
on a rho/2 fraction of seeds. All other inputs have prediction success
less than rho/2. Consequently the predictor's total success probability
is at most

    2^(a_adv-k+1)+rho/2 <=rho

when k>=a_adv+log_2(4/rho). This inequality holds eventually: a_adv is
O(m*k^(1/32)*log N). Dividing by k gives order
L^(-7/32)*J^(-1/32), up to fixed positive factors; this tends to zero.
Here k=Theta(L^25/J) and m=Theta(L^23/J).
This proves the q-ary unpredictability statement, including predictors
that see the seed v.

For a field value z, the binary word (<z,j>)_(j in F_2^w) is its Hadamard
encoding. Distinct words are orthogonal in the sign representation.
For any proposed binary word, Bessel's inequality says that at most
rho^-2 Hadamard words agree with it on a fraction at least 1/2+rho.
This follows directly by summing the squared correlations with the
orthonormal characters; each such correlation is at least 2rho.

Suppose the binary extractor output together with its seed (v,j) has
variation distance greater than 2rho*m from a uniform m-bit output with
the same seed. The telescoping hybrid over output bits gives some index
and a binary predictor, seeing (v,j) and the previous output bits, with
average prediction advantage greater than 2rho. For at least a rho
fraction of pairs (x,v), its prediction word as j varies agrees with
the correct Hadamard word on more than a 1/2+rho fraction of positions.
Otherwise the average advantage would be at most rho+rho/2<2rho.
Using the preceding field outputs, the predictor can form this prediction
word for all j and list its at most rho^-2 close Hadamard codewords.
It is therefore a q-ary predictor seeing v with success greater than
rho, a contradiction. Thus the strong binary error is at most

    2rho*m=N^-8.

All code operations are linear in the input for fixed (v,j), as required.

## 4. The combined short linear affine extractor

Apply Cond_a to the N-bit source, followed by the Section 3 extractor
on its n'=k*w_1 output bits, using independent seed (v,j). Call their
composition Short. For all but N^-15 of a, the condensed affine source
has min-entropy at least k. Strong extraction on each of those sources
and averaging over a prove strong error at most N^-15+N^-8<=2N^-8.
Its output length is m=floor(k/(2^300*L^2)); its total seed length is

    d_short=w_1+65w,

which is at most 16384 log_2 N at sufficiently large integer N. Both fixed
seed maps are linear. Therefore for an affine source every fixed-seed
output is uniform on an affine image. If that image has rank less than
m, its distance from uniform is at least 1/2. It follows that a fraction
at most 4N^-8 of seeds are bad, and EVERY other seed gives an EXACTLY
uniform m-bit output. This exact-rank fact will supply the good rows.
The good set depends on the source direction space, never its sample.

The algorithms use only polynomially bounded fields and a fixed degree
64 extension over them. No field search is over an advice-selected
polynomial; all choices are deterministic searches whose existence was
proved in Section 1. The seed is an internal enumerated choice in the
final construction and is not an additional target-circuit input.

## 5. Internal linear code-bit extractors

The remaining alternating-extraction calls need longer seeds and very
small error. They use the complete linear Reed--Solomon/Hadamard
code-bit argument reproduced next. Their entropy constants will be
retuned to a small fixed output fraction; no near-full output theorem
is used as an unexplained primitive.

### 5.1 A binary linear code with a small combinatorial list

Given n data bits, group them into w-bit field elements (zero-pad the last),
interpret them as coefficients of a polynomial P of degree below h=ceil(n/w)
over F_(2^w), and index code coordinates by pairs (a,b) of field elements.
Output the binary dot product b . P(a), in a fixed binary basis. This is
linear in the original n bits. Its codeword length is 2^(2w), and a code bit
can be computed in poly(n,w) bit operations if the field is supplied.

For distinct data words, P-P' is nonzero and has at most h-1 roots. Averaging
the sign of b.(P(a)-P'(a)) first over b shows that the normalized inner product
of their sign codewords is between 0 and (h-1)/2^w. Choose 2^w>=h/gamma^2.

Any fixed sign word has correlation at least 2gamma with at most
1/(3gamma^2) codewords. Indeed, for M such codewords v_i and target v,
with normalized Euclidean norms 1, Cauchy-Schwarz gives

    (2gamma M)^2 <= ||sum_i v_i||^2 <= M+gamma^2 M(M-1).

Thus M <= (1-gamma^2)/(3gamma^2). No list-decoding algorithm is needed:
only this counting bound enters the extractor proof. Index length ell=2w.

### 5.2 Reconstruction estimate for an arbitrary design

Let S_1,...,S_m be ell-element subsets of [D] with pairwise intersection
at most a. For uniform seed y in F_2^D, output

    E(x,y)_i = Code(x)[ y restricted to S_i ].

This is linear in x for each seed. A sufficient condition for it to be a
strong (k,epsilon) extractor is the following conservative inequality:

    k > D + m*2^a + log_2(m) + 4 log_2(8m/epsilon) + 10,

where the code is set up with gamma=epsilon/(8m).

Proof: if a distinguisher separates
(Y,E(X,Y)) from (Y,U_m) by more than epsilon, hybridizing the output bits
gives a fixed index i and a predictor for its true code bit, from the seed
and the preceding true bits, with average advantage at least epsilon/m.
The hybrid gap equals one half of E[true_sign*Delta], so a prediction
with expected sign Delta has correctness advantage equal to that gap.
Formally let Delta be the difference between the distinguisher's two
candidate-bit values, averaged over the remaining independent uniform
output bits, oriented so the hybrid difference is positive. Use a randomized
prediction whose expected sign is Delta. Its correlation with the true sign
is the required hybrid correlation. Assign its random choices once for every
possible seed/prefix argument and fix the resulting finite random table by
averaging. Some single fixed deterministic predictor retains the average
correlation. This predictor may depend on the source distribution, but is
then common to all x; neither it nor its fixed random table is charged as
per-x advice. Taking sign(Delta) directly would not justify the correlation.

At least epsilon/(2m) of source mass has prediction advantage at least
epsilon/(2m) after averaging the seed (using advantage bounded by 1).
For each such x, choose values of seed bits outside S_i that retain this
advantage. Store those values using at most D bits. Each preceding code
bit, after the outside seed is fixed, depends only on S_j intersect S_i;
its entire truth table costs at most 2^a advice bits. All preceding bits
therefore require at most m*2^a bits. The predictor with this advice gives
a sign word on the ell-bit code index with correlation at least
epsilon/m, hence at least 2gamma.

For each fixed advice the code list bound leaves at most 1/(3gamma^2)
possible x. Thus the set of such x has cardinality at most
2^(D+m*2^a)/(3gamma^2), and its probability under any min-entropy-k source
is at most this cardinality times 2^-k. The displayed conservative k
inequality makes this less than epsilon/(2m), a contradiction. The stated
inequality includes more logarithmic slack than the count requires.

### 5.3 Polynomial-graph design and the output fraction

Take an integer a>=1. Let Q be a power of two with ell<=Q<2ell.
Use a fixed ell-point subset T of F_Q and let S_i be the graphs over T
of distinct degree-at-most-a polynomials. There are Q^(a+1) indexed graphs,
which are distinct when a<ell, as eventually holds for the fixed A=4096
used below. Their universe has D=ell*Q<2ell^2 elements, and intersections have size
at most a by the polynomial root bound. Enumerate their coefficients
lexicographically. As long as m<=Q^(a+1), the design is explicit.

For a desired input entropy k polynomially large in the ambient length n,
choose epsilon=2^(-ceil(k^e)) for a sufficiently small positive constant e.
Then w=O(log n+k^e), ell=2w and D=O((log n+k^e)^2)=o(k) if e<1/2.
Choose fixed a large enough that Q^(a+1)>=k, and output
m=floor(k/2^(a+4)). The advice term m*2^a<=k/16 and all other terms are
o(k). This gives a constant-fraction output, a seed of length O(k^(2e)),
and error 2^(-ceil(k^e)) at every sufficiently large length.

The constant fraction can be extremely small because a may be large.
Rao's condensation must therefore be retuned to this fraction; it cannot
retain the published nearly-full-length output parameters.


### 5.4 Large code fields do not require exhaustive search

The code field degree can exceed O(log N), so Section 1 exhaustive
search is NOT used for those fields. Round the code degree w upward to
a power of two and construct its field by an explicit quadratic tower.
Start with F_2 and trace-one element c=1. Given a field K with a chosen
absolute-trace-one element c, adjoin t with t^2+t+c=0. The map z->z^2+z
on K has kernel {0,1}; its image has size |K|/2, is contained in the
trace-zero kernel, and equals that kernel because trace(c)=1. Thus the
polynomial is irreducible. In the new field, t*c has relative trace c
and absolute trace one, so it supplies the next extension's element.

Multiplication of (a+b*t)(d+e*t) is

    (ad+be*c)+(ae+bd+be)*t.

Three products ad,be,(a+b)(d+e), and one product be*c in the base field
suffice. Hence T(2w)<=4T(w)+poly(w), giving polynomial bit cost in w.
The representation and trace-one elements are completely deterministic.
No factorization search at an exponentially large field size is hidden.
For the small polynomial-graph field, degree is O(log ell), so Section 1
exhaustive search remains polynomial in N.

### 5.5 Exact internal contract and padding seeds

For ambient length n, requested output m and error 2^-E, set

    w_target=bitlength(n)+2bitlength(8m)+2E+10,
    w=the least power of two >=w_target, ell_code=2w,
    Q=the least power of two >=ell_code,
    D=ell_code*Q, a=4096.

Use the tower field of degree w for the linear code and the
polynomial-graph field of size Q for the design. The stated w ensures
2^w>=ceil(n/w)/gamma^2 with gamma=2^-E/(8m). Guard
m<=Q^(a+1) and a<ell_code. Define the sufficient entropy threshold

    H_call=2D+2^(a+4)*m+8*(E+ceil(log_2(m+2))+10).

It strictly exceeds the Section 5.2 reconstruction right side. Thus
this exact evaluator is a strong linear (H_call,2^-E) extractor.
If a caller supplies s>=D seed bits, it uses the first D and ignores the
rest. A uniform s-bit seed has a uniform independent D-bit prefix;
appending unused seed coordinates to the joint distribution preserves
strong extraction. All fixed-seed output maps remain linear.

In the main construction write L=floor(log_2 N). Low-error calls use
E=2^200*L^2 and supplied seed ell=2^500*L^4. Since every requested
output has bit length O(log L), w_target<=3E eventually. Rounding
therefore gives w<=6E, ell_code<=12E, Q<=24E and
D<=288E^2<2^409*L^4. Thus D is less than the supplied ell by a fixed
factor exceeding 2^91. The design has capacity at least a fixed
positive multiple of L^(2*(a+1)), far exceeding the largest requested
output, which has order L^20. All low-error entropy thresholds are
at most C_L times the requested output, C_L=2^4110: the dominant
term is 2^4100 times that output, D is below ell/2^91, and the other
terms are lower order than the smallest output ell.

For the exceptionally low-error first refresh,
E_1=2^20*h^3*ell+E and d_SR=2^64*(h^3*ell+L)^2. Eventually its
w_target<=4E_1, so its seed D<=512E_1^2<=2^51*h^6*ell^2<d_SR/8192.
Its larger H_call is charged against C_L*d_SR, not C_L*ell.
The full main parameter proof includes the remaining exact guards.

All powers, roots, ceilings and design order are fixed deterministic
operations. The asymptotic assertions here require their actual coded
seed/capacity inequalities to be checked in the family metadata. Failure
at an exceptional small length returns the documented zero branch;
the main proof must establish their eventual success through every
sufficiently large integer N.

