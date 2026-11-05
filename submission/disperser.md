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

    q >= (2^20*r*h/rho^4)^8,
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

For both displayed inequalities, note r=16D+16 and therefore
q^(1-D/r)>=q^(1/8)>=2^20*r*h*rho^-4. The first left-hand
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



# B1. Deterministic linear rank condensation

## Fields without per-length advice

For every power-of-two integer s construct F_(2^s) by an Artin--Schreier
tower. At degree 1 the field is F_2. Given the degree-h field, adjoin T
with T^2+T+a=0, where a is its top basis vector, the product of all
previous tower generators (and a=1 when h=1). Its absolute trace is one.
Indeed Tr_(2h)/h(u+vT)=v; applying this identity successively proves the
assertion for the top basis vector. The image of x->x^2+x is exactly the
trace-zero hyperplane: its kernel has size two, its image has codimension
one, and trace annihilates it. Thus the quadratic is irreducible at each
step. The embedded degree-h field occupies the low h coordinates.
Multiplication is

 (u+vT)(r+wT)=(ur+a*vw)+(uw+vr+vw)T.

The implementation evaluates its five smaller-field products by an explicit
postorder stack. There is no unbounded interpreter recursion. The recurrence
uses polynomial time in the field degree; shifts and XORs have polynomial
bit cost.

Let s,D be powers of two with D>=2. In the degree-sD tower let theta be
the last quadratic generator. It lies outside the degree-sD/2 subfield.
Every proper intermediate field over F_(2^s) is contained in that subfield,
since its degree divides the power of two D. Hence theta has relative
degree D. Put q=2^s and

 E(Y)=product_(j=0)^(D-1) (Y+theta^(q^j)).

Frobenius permutes its roots, so its coefficients belong to the embedded
F_q. It is the irreducible minimal polynomial of theta. Computing the
product and the q-th powers requires polynomial time in sD and s. The
routine constructs E from its integer dimensions; no input prefix of
length comparable with the ambient input is hidden here.

## The polynomial graph and its expansion

Fix a field F_q, an irreducible E of degree D>=2, and h=2^a. For a
polynomial f of degree below D, define f_i=f^(h^i) mod E. For a seed y
in F_q output the m symbols f_0(y),...,f_(m-1)(y). Each fixed-seed map
is F_2-linear. View the seed and output together as a graph neighbor
Gamma(f,y)=(y,f_0(y),...,f_(m-1)(y)). Let

 A=q-(D-1)(h-1)m.

For every integer J<=h^m and every set S of J polynomials, its union of
neighbors has size at least A*J. Here and below a vacuous nonpositive A
causes no difficulty, although our parameter choices make A positive.

To prove this, suppose a set T of fewer than A*J neighbors contains all
neighbors of at least J polynomials. Interpolate a nonzero polynomial

 Q(Y,Z_1,...,Z_m)=sum_(j=0)^(J-1) p_j(Y) M_j(Z_1,...,Z_m),

where deg p_j<A and M_j has the base-h digits of j as exponents. There
are A*J unknown coefficients and fewer constraints. Among solutions choose
one of minimal Y-degree. Not all p_j are divisible by E: otherwise divide
Q by E; E has no F_q root, so the quotient still vanishes on T and has
smaller Y-degree. For a polynomial f whose neighbors all lie in T,
Q(Y,f_0(Y),...,f_(m-1)(Y)) has degree at most
(A-1)+(D-1)(h-1)m<q and vanishes at all q field points. It is identically
zero. Reducing modulo E, f is a root over the extension field of

 Q*(Z)=sum_(j=0)^(J-1) (p_j mod E)*Z^j.

This is nonzero and has degree at most J-1. It cannot have J distinct
roots. This contradiction proves expansion. Interpolation is used only in
the proof, not as an algorithm that enumerates source points.

## Exact parameters and rank preservation

For n input bits, an integer entropy threshold K>=1, and error log e,
let

 a=ceil(log2 n)+ceil(log2 K)+e+10,
 s=the least power of two >=4a,
 D=the least power of two >=max(2,ceil(n/s)),
 m=ceil(K/a), h=2^a.

Pack the input bits injectively into D field coefficients, padding zero.
Use the deterministic tower and E above. The seed has s<8a bits and the
output m*s<=8K+8a bits. Since D<=2n+2 and m<=K+1, the chosen s gives
(D-1)(h-1)m<=2^-e*q, directly by taking binary logarithms; the additive
10 leaves ample room at the smallest positive integers. Thus A>=(1-2^-e)q.

For an affine input of dimension exactly K, each seed output is an affine
image of rank r_y<=K. Its graph-neighbor count is sum_y 2^r_y. Expansion
therefore gives E_y[2^(r_y-K)]>=1-2^-e. Whenever r_y<K this summand is at
most 1/2, so the fraction of rank-deficient seeds is at most 2^(1-e).
On every other seed the output is an affine source of dimension K. For
larger input dimension decompose it into parallel K-dimensional flats and
average; alternatively the image rank bound follows from any such flat.
This is the precise affine rank-condensation statement used below.

# B2. Two supplied strong linear extraction routines

An affine distribution has a vector-space complement description: if Z is
an affine function of X, then X=A+B, with A an affine function of Z and
B an independent uniform linear source. Its dimension is H(X|Z).
Conditioning on a fixed linear observation of t bits reduces any remaining
affine rank by at most t. An observation determined by Z does not reduce
H(X|Z). These statements follow by choosing a basis of the projection
kernel and extending it to a basis of the direction space.

For the initial and short-output E1 maps, consider an affine extractor plan for source length n, entropy K, output
M<=K/2 and precision 2^-e. First use B1 with precision e+10. On every
good first seed its output is affine of entropy at least K. For its ambient
length n'=O(K+log n+e), use A1 with the plan's stated dimension D and
per-call error 2^-(e+ceil(log2 K)+10). Choose the largest power-of-two
output m0<=min(M,K/4) satisfying A1's condition with threshold floor(K/2).
The executable plan halves m0 until it succeeds. Failed finite parameter
cases are rejected by the outer ready flag.

Repeat A1 on the SAME condensed input with fresh independent seeds until
M bits have been produced, truncating the last call. Earlier fixed seeds
and outputs leave an affine source of entropy at least K-M>=K/2. A
sequential coupling loses at most the sum of call errors; there are at
most K calls. Adding the B1 bad-seed probability gives total strong error
at most 2^-e. Each fixed combined seed is a linear map of the original
input. Padding a longer stipulated seed with ignored bits preserves this.

We prove a budget retaining the dimension explicitly. Put
H=D*ceil((n')^(1/D)), R=16D+16. A1's field degree satisfies
s=O(e+log K+log H+log R). Substituting a sufficiently small fixed
multiple of K/[R*H*(e+log K+log H+log R)] for m0 in A1's condition
shows that the largest admissible power of two is at least a fixed
multiple of that value, unless M is smaller and one call suffices.
The additive logarithmic term is o(K) in all our uses. Thus the combined
seed length is at most

 C(log n+e+log K)
 +C(D+1)*R*H*(e+log K+log H+log R)^2.                    (B2.1)

Here C is independent of the displayed dimensions. This proves eventual
success of the halving loop as well as its seed cost. Only actually
computed scalar dimensions enter prefix size or executable seed slicing.

Initial extraction keeps D=10000 FIXED. It has source k, output Y and
precision 4L+10, where B7 specifies k,Y and L=ceil(log2 N). For every
FIXED scale S, one call eventually suffices: its advice has order
Y*k^(.0001)*L=o(k): its L exponent is 14.0015<15.
Its combined seed length d0 is at most 10^9 L eventually. Directly from
A1, the field-degree leading term is at most 128L, multiplied by
D+1=10001, plus O(L) condensation bits. Log Y and log k are o(L).
The finite plan checks this cap before enumerating any initial seeds.

## B2.2 Direct half-output weak-design extraction

The other four internal maps use HalfSeedPlan, not the repeated A1 routine. Let n,K,M,E be positive integers with K<=n, M<=K/2. Put B=E+ceil(log2 M)+10 and require K>=32B. First apply B1 with precision E+10. Write nprime for its exact output length. Set

 s=P2(8*(E+ceil(log2 nprime)+ceil(log2 M)+10)),
 a=2s, T=100a^2.

Use the deterministic degree-s Artin--Schreier tower field from B1; s is a power of two. No field modulus, matrix or row table is input data for this map. Construct M rows of size a in [T] using the exact conditional-expectation weak-design algorithm proved in B5. Its proof depends only on a and the number of rows, so it applies unchanged here and gives sum_(j<i)2^|A_i intersection A_j|<=1.5i. The seed is the B1 seed followed by T uniform bits; unused stipulated bits are ignored. The exact seed length is B1.s+T and data length is zero.

Pack the condensed nprime bits injectively as coefficients of a polynomial P_X over F_(2^s). A row seed consists of independent y,z in that field. Its output bit is the parity of the tower-basis coordinates of z*P_X(y). This is a nonzero binary linear functional (it takes value1 on the field element1). Every fixed-seed map is binary linear in the original source. A nonzero polynomial difference has degree below ceil(nprime/s), hence at most nprime roots. At all other y, multiplication by its value makes z uniform, so the parity difference is uniform.

For any source of min-entropy H, independent copies coincide with probability at most2^-H. The collision calculation just given, followed by Cauchy--Schwarz averaged over the ENTIRE seed(y,z), bounds strong one-bit distance by

 (1/2)*sqrt(2*(2^-H+nprime/2^s)).

This is at most2^-B when H>=2B+2 and s>=2B+ceil(log2 nprime)+2, with room in the constants. The chosen s has this property. The argument works for arbitrary min-entropy sources, not only affine distributions; the conditional sources below therefore need not remain affine.

For a condensed source of min-entropy K, apply the direct table hybrid of B5 to the M output bits, retaining every seed bit in the output. At step i, fix seed coordinates outside A_i and the later independent uniform hybrid bits. List all previous outputs for every assignment on their overlaps with A_i. This table depends only on X and the fixed outside seed, not on the current uniform row seed. Its length is at most1.5i<=1.5M<=3K/4. Table values of probability below2^(-length-B) have total mass at most2^-B. On every other value, conditional min-entropy is at least K-length-B>=K/4-B>=7B, using K>=32B. Thus the one-bit bound applies to every test conditioned on this table, with error at most2^-B. Each hybrid step loses at most2*2^-B, including exceptional tables. Summing M steps loses at most2M*2^-B<=2^(-E-9).

For every good B1 seed the condensed affine source has entropy at least K. The B1 bad-seed fraction is at most2^(1-(E+10))=2^(-E-9). Averaging good-seed bounds and adding bad mass yields total strong error at most2^(-E-8)<2^-E. This proves the stated HalfSeedPlan extractor, INCLUDING its whole seed. For larger affine sources the rank statement follows from any K-dimensional flat and averaging parallel cosets.

All construction is deterministic and polynomial in the displayed finite parameters: the B5 row algorithm uses at most M*T coordinate decisions with at most M prior rows and a+1 binomial terms per decision, with polynomial bit lengths. Field operations use the postorder tower algorithm of B1. Horner evaluation uses ceil(nprime/s) field coefficients per output bit. Selected row coordinates are extracted by iterating set bits, not scanning T coordinates after construction. This is a runtime improvement only; it preserves exactly the same row set as the stated conditional-expectation rule.

For the four internal maps, E=h dominates log n, log K and log nprime eventually at each fixed scale. Their exact seed length is O(E^2), with an absolute finite constant, and K>=32B eventually. This replaces the former O(z^3*E^2) half-output bound. The short E1 map STILL uses A1 with D=z and one call: its advice is O(ell*z^2*Eprime), and seed is O(log N+z*Eprime). The initial map STILL uses fixed D=10000 and its O(log N) seed cap. Neither is replaced by the half-output routine.

At fixed good construction data, an affine output either is exactly uniform or has distance at least1/2 from uniform. Strong error eta bounds bad-seed fraction by2eta. An independent seed of min-entropy a-t can hit a bad set with probability at most2^(t+1)*eta. These rank observations used by B3--B4 are unchanged because the new maps are fixed-seed linear.

# B3. Parallel affine look-ahead extraction

Let at most h streams Q_i be affine functions of Z, each padded to a common
ambient bound. Let a<=h be the look-ahead time limit, separate from the
number of streams. One Q_v is exactly uniform on q0=(2h*a+10)ell bits.
Let H(X|Z) be at least h*a*ell+10ell. The extractors W on X and T
on Q_i have output and seed length ell, source entropy threshold 10ell,
and strong error eta. Set

 S_i,1=first ell bits of Q_i,
 R_i,t=W(X,S_i,t),
 S_i,t+1=T(Q_i,R_i,t).

For t<=a, the designated R_v,t is uniform conditional on all Q_i, all
current and earlier S values, and all earlier R values, except on total
history mass at most 4a*eta. More precisely, after the preceding messages
and the current seeds are fixed, split X=A+B by Z. Outside that error
mass the map from B to R_v,t has full output rank. It therefore remains
uniform under any further linear fixing determined by Z. This stronger
rank form is what is used in B4.

Here is a direct proof, including the conditioning order. Initially the
designated S_v,1 is uniform. At a step fix all current S_i,t. These are
linear observations of Z once earlier R values are fixed, so the B kernel
is unchanged. Its bad-seed set for W is fixed, and the marginal uniformity
of S_v,t bounds the bad history mass by 2eta. On a good seed the B-part
of R_v,t is exactly uniform and independent of A. This holds even though
the other current seeds may be arbitrarily correlated with it.

To advance, fix the A-parts of ALL current R_i,t, including the designated
one. Each is a linear observation of Z because its seed is already fixed.
The designated R_v,t is now a constant translate of a uniform B-part,
independent of A. Apply T to Q_v with that seed. Fix R_v,t and then the
other B-parts. The output S_v,t+1 is linear in A and exactly uniform except
for another bad-seed mass 2eta; the other B-fixings cannot affect it.
All next seeds are now linear functions of Z. This proves the induction.

Each completed step fixes at most h*ell bits from the seeds and h*ell
A-parts, hence reduces Q_v's rank by at most 2h*ell. The B dimension falls
by at most h*ell. Before the last needed transition, q0-2h*(a-1)ell>=10ell and
the remaining B dimension exceeds 10ell. The initial uniform seed and the
sequential rank conclusions give the claimed cumulative error. All
conditioned components are affine: seeds are fixed before a seeded map's
output is conditioned on, and every subsequent observation is then linear.
A good-history flag is determined by a fixed transcript value and does
not impose a further nonlinear condition inside that affine component.

# B4. Many almost independent affine rows

The following explicit variant of index-based alternating extraction is
supplied in full because its conditioning order matters. It follows the
architecture of Li 2015, but refreshes the look-ahead source on EVERY
index round, including the last. No theorem from that paper is used as
an unexplained extractor premise.

Let k,h,ell,d be the parameters of B7 and let Y_i be the initial B2
extractor output on seed i, of length Y. There are M0=2^d0 seeds.
The index alphabet is a=2 and b=d0. Write i as its d0-bit binary word,
and add one to each bit to obtain t_i,1,...,t_i,b in {1,2}. This is
injective and preserves numerical order. The independence order h is
not the index alphabet; it still bounds the number of parallel streams.

For each i start D_i,1=first d bits of Y_i. On every round j compute

 R0_i,j=E1(X,D_i,j),
 Q_i,j=E2(Y_i,R0_i,j),
 (R_i,j,1,...,R_i,j,a)=look-ahead(X,Q_i,j),
 V_i,j=R_i,j,t_i,j.

For j<b put D_i,j+1=E3(Y_i,V_i,j). The row W_i is V_i,b. It is enough
for the evaluator to compute each look-ahead stream only up to its selected
digit. The full prefix messages used in the proof are its identical
mathematical continuation and do not change the selected output.

The source/seed/output specifications are:

* E1: ambient N, entropy 10d, output ell, seed d, error eta'=2^-e'.
* E2: ambient Y, entropy (4h*a+20)ell, output q0, seed ell,
  error eta=2^-e.
* E3: ambient Y, entropy 2d, output d, seed ell, error eta.
* W: ambient N, entropy 10ell, output ell, seed ell, error eta.
* T: ambient d, entropy 10ell, output ell, seed ell, error eta.

Short Q_i are padded at the END before T is evaluated. Construction B2
provides all five maps with these seed lengths, as verified in B7.
Assume the following conservative bounds, also checked eventually there:

 Y>=1000*b*h*d,
 k-h*Y-100*b*h*a*ell>=10d,
 d>=(2h*a+10)ell,
 e'>=100h*a*ell+2e.                                      (B4.1)

For every affine source of entropy at least k, all but at most
2M0/N^2 initial rows Y_i are exactly uniform, eventually. For ANY set of
at most h of those good indices, its W_i rows are jointly within
100*b*h^2*eta of independent uniform ell-bit rows.

## Proof by ordered prefixes

Fix such a set I, and fix a designated v in it. Let Z concatenate its
original Y_i; these are linear functions of X, and H(X|Z)>=k-h*Y.
Let E_j be predecessors i<v whose first j index digits are strictly
smaller than v's, and H_j the predecessors with the same j-digit prefix.
At j=0, E_0 is empty. At j=b, H_b is empty because the full index words
are distinct. A same-prefix predecessor's next digit is never larger
than v's digit.

We maintain an adaptive transcript of fixed seed and linear-output values.
At the start of round j, all previous chosen R values are fixed, so all
current D_i,j are linear functions of Z. The seeds D_i,j for i in E_(j-1)
are fixed as part of the transcript. Outside accumulated bad history mass,
D_v,j is exactly uniform conditional on that transcript. Moreover

 H(Y_v)>=Y-100(j-1)h*d,
 H(X|Z)>=k-h*Y-100(j-1)h*a*ell.                     (B4.2)

Initially this follows from uniformity of Y_v. The additional transcript
values described below make the induction stronger, not weaker: whenever
an A-part is fixed it is a linear function of Z, and whenever a seeded
output is fixed its seed has already been fixed.

First process the already-diverged rows E_(j-1). Their current seeds are
fixed. Fix R0_i,j, after which Q_i,j is a LINEAR function of Z, but do NOT
fix the whole Q_i,j. Instead fix its first seed S_i,1, then R_i,1, then
S_i,2, then R_i,2 if needed. In general each next S is linear in Z once
the preceding R is fixed, and each R is linear in X once its seed is fixed.
This conditions only the actual seeds and outputs up to at most a steps.
The full intermediate row Q_i,j need not be revealed.

For at most h already-diverged streams, this costs no more than
(2a+1)h*ell<=100h*a*ell bits from D_v,j and O(h*a*ell) from H(X|Z).
Thus D_v,j retains entropy at least d-100h*a*ell. Every chosen R of
these predecessors is fixed, which is precisely what makes their next
seed linear in Z. No same-prefix seed of length d has yet been fixed.

Split X=A+B by Z in the resulting affine component. The target seed is a
function of A and independent of B. Its weak-seed probability of choosing
an E1 map of deficient rank on B is at most
2^(100h*a*ell+1)*eta'<=eta, by (B4.1). Fix D_v,j first. On a good value
the B-part of R0_v,j is uniform independently of A. Now fix the current
seeds of H_(j-1); these are A-only observations and cannot spoil that
rank. Fix the A-contributions to R0_i,j for EVERY remaining stream,
including v itself. This last inclusion ensures the target R0 is purely
a translate of a B-variable, rather than a seed varying with Y_v.

The uniform R0_v,j is now independent of all Y_i. E2 therefore produces
a uniform Q_v,j except for bad-seed mass 2eta. Fix R0_v,j and then the
other B-contributions to R0_i,j. All Q_i,j are now linear functions of Z;
Q_v,j remains uniform. Its dimension is q0, and the source B entropy is
more than sufficient for B3 by (B4.1)--(B4.2).

Apply B3 to the target and the same-prefix streams up to the target's
current digit t. Outside mass at most 4a*eta, its final B-part is uniform
and stays so under further Z-linear fixings. A predecessor with digit
less than t has its chosen R already fixed in this prefix. Its next seed
E3(Y_i,chosen R) is therefore linear in Z. Fix these next seeds and those
of the already-diverged predecessors. Also fix the Q list and the A-parts
of the current R messages for all streams, including the target. These
are all A-only fixings after the current look-ahead seeds are fixed;
therefore the target final B-rank survives them. The remaining Y_v entropy
is still at least 2d by the first bound in (B4.1).

For j<b apply E3 to Y_v using the uniform target chosen R as an independent
seed. Fix this seed, and then fix the other current B-parts. Outside another
bad mass 2eta the resulting D_v,j+1 is linear in A and exactly uniform.
The other B-fixings leave this uniformity intact. Do NOT fix the next
seeds of still-same-prefix predecessors: those may be correlated with it.
They are fixed only on the next round, AFTER the target seed is fixed.
This gives precisely the induction assertion conditional on E_j next seeds.

All original Y_v rank losses in one round are bounded by 100h*d: the
seed observations and next-seed observations cost O(hd), while Q lists
and look-ahead A-parts cost O(h^2*a*ell)<=O(hd). Only O(h*a*ell) output bits
per round can reduce H(X|Z); observations determined by Z do not. These
bounds include the target's own previous seeds, Q values and A-parts.
Thus (B4.2) is preserved with the displayed conservative constants. The
failure mass per round is at most 10h*eta.

On the LAST round the same E1/E2 refresh is used. Every predecessor is
already diverged or has digit less than t, since H_b is empty. Hence all
its FINAL chosen W_i values are fixed before the target R_v,b,t. The
target row is uniform conditional on them outside the accumulated error;
no final E3 call is needed. Applying this conclusion to the ordered indices
and coupling one row at a time proves joint error at most
100b*h^2*eta. This proof permits arbitrary correlations among original
Y_i and among all nondesignated streams.

# B5. A supplied weak-design extractor with a squared-log seed

We supply the exact recycling argument needed here, following the
Nisan--Wigderson/Trevisan weak-design architecture. No extractor theorem
from that literature is assumed. Write w=ceil(log2 r), delta=10^-12,
m=floor(delta*r), s=1000w, c=ceil(s/4), a=2s. A single good actual-input
modulus defines F_(2^s). Pack the r-bit v injectively into coefficients
of P_v over that field, degree below ceil(r/s). For field seeds y,z put

 C(v;y,z)=first c coordinates of z*P_v(y),
 E0(v;y,z)=parity of the binary coordinates of z*P_v(y).

Every fixed-seed map is binary linear. For a nonzero difference v, P_v
has at most r roots. At every other y, multiplication by P_v(y) makes z
uniform. Thus every c-bit offset has probability at most

 2^-c+r/2^s.                                             (B5.1)

For ANY K-min-entropy source, independent copies coincide with probability
at most2^-K. The same collision calculation and Cauchy--Schwarz give
strong distance of C from uniform, INCLUDING its field seeds, at most

 (1/2)*sqrt(2^c*(2^-K+r/2^s)).                           (B5.2)

For the one-bit map E0 the identical bound has2 in place of2^c. Therefore
E0 has strong error at most r^-400 whenever K>=1000w, eventually.
For K=2m the C-prefix has strong error at most r^-300; (B5.1) is at most
r^-200. These estimates use s=1000w,c=250w and c=o(m).

## Constructing the weak design by exact conditional expectations

Let T=100a^2. Construct M=m-c rows A_i, each an a-subset of [T]. Given
rows with j<i, construct an a-subset A_i satisfying

 sum_(j<i)2^|A_i intersection A_j| <=(3/2)i.              (B5.3)

The code compares integers2*sum<=3i, and permits repeated rows. Here is
an elementary existence proof for every step. For a uniform a-subset A
and any fixed a-subset B, expansion of2^|A intersection B| gives

 E[2^|A intersection B|]
 =sum_(l=0)^a binom(a,l)*(a)_l/(T)_l
 <=(1+a/(T-a+1))^a
 <=exp(a^2/(T-a+1))<1.02.

Summing over i prior rows gives an average below1.02i<1.5i for i>0,
so a suitable row exists; i=0 is immediate. This gives the initial potential for the deterministic conditional-
expectation construction below.

To construct a row, maintain selected coordinates P, undecided coordinates
R, and k=a-|P| slots still to fill. Initially P is empty, R=[T], k=a.
Use the average over UNIFORM k-subsets of R as potential. If t=|R|,
then for a previous row B, h=|P intersection B| and l=|R intersection B|,
expanding2^overlap into subset indicators gives its exact expected weight

 2^h * sum_(j=0)^min(k,l) binom(l,j)*binom(t-j,k-j) / binom(t,k).

Sum these numerators over all previous rows, using the common positive
denominator binom(t,k). At the next coordinate, compute the include and
exclude potentials with this formula and choose the SMALLER by integer
cross-multiplication, choosing include on a tie. Their weighted average,
with weights k/t and1-k/t, is the current potential. Thus the chosen
potential never increases. If k=0 exclude all remaining coordinates;
if k=t include them all. The result has size exactly a and its final
integer score is at most its initial average below1.02i, hence satisfies
(B5.3). For i=0 all potentials are zero and the same algorithm applies.
No search over row subsets, nonconstructive choice, or unpaid advice is
used. The executable assertion checks the same final weight inequality.

There are at most T decisions per row, i previous rows and at most a+1
terms per expected-weight sum. Binomial integers have O(a log(T+1)) bits,
and row masks have T bits. With at most M<=r rows and T=100a^2, the
entire construction uses polynomial(r,w) time and memory. Rows are computed
once in Design.bind and reused. Evaluating a one-bit map also takes
polynomial(r,w) operations. This improves the preceding exhaustive
construction, whose2^O(w log(w+1)) time was only polynomial in actual N.

## Strong extraction by a direct hybrid and table argument

For U uniform on T bits, let NW(X,U)_i=E0(X,U restricted to A_i), in the
increasing coordinate order of A_i. We prove strong error at most r^-90
for ANY source with min-entropy at least7m/4, eventually.

Compare hybrids with the first i outputs genuine and all later outputs
uniform, keeping the ENTIRE U seed in the output. To bound one hybrid
difference, fix U outside A_i and fix the later independent uniform bits.
For each j<i, list E0(X,U|A_j) for ALL assignments on A_j intersection
A_i, keeping its coordinates outside A_i at the fixed values. This is
a table depending only on X and the fixed outside seed. Its total length
ell_i is at most the sum in(B5.3), hence at most1.5i<=1.5m. Given the
table, every previous output is a deterministic lookup from U|A_i.
The table is independent of U|A_i, which is still uniform.

For each fixed outside seed, table values of probability below
2^(-ell_i-100w) have total mass at most2^-100w. On every other table
value, the conditional source has min-entropy at least

 7m/4-ell_i-100w >=m/4-100w>=1000w

for sufficiently large r. The strong one-bit guarantee for E0 bounds
any test of the current bit, the seed and the now table-determined earlier
outputs by r^-400. Adding the exceptional table mass, the one-hybrid
difference is at most r^-400+2^-100w, uniformly in the fixed outside
seed and later bits. Averaging and summing at most M<=r hybrids gives
at most r^-90 eventually. This proof works on general min-entropy sources;
there is no assertion that conditioning on a table preserves uniformity.

## Appending to the hash prefix and fitting the seed

Use independent seed portions for C and NW and output

 G(v;seed)=C(v;first2s bits) concatenated with NW(v;next T bits).

For each fixed hash seed, hash-prefix values of probability below
2^(-c-100w) have mass at most2^-100w. On every other value a2m-source
retains min-entropy at least2m-c-100w>=7m/4 eventually. Thus the strong
NW guarantee applies after conditioning on that prefix. Its seed is fresh
and independent. Combining this error, exceptional mass and(B5.2), the
WHOLE G is strong on2m-sources with error epsilon<=r^-80. Its hash prefix
alone gives pair-shift agreement eta<=r^-200, for ANY pair and offset.
Fixed-seed linearity is retained by every output bit.

The entire seed has length2s+T<=C_design*w^2, with a fixed explicit finite
constant C_design. It fits

 b=floor(floor(sqrt(S_header))*w^2/(30*10^24))

for one sufficiently large fixed header S_header. Unused bits are ignored.
The only design construction data are the s actual field-modulus bits;
there is no design A1 package or hidden weak-design table.

Provenance: [Trevisan](https://theory.stanford.edu/~trevisan/pubs/extractor-full.pdf)
and the [weak-design approach of Raz--Reingold--Vadhan](https://people.seas.harvard.edu/~salil/research/trev-abs.html).
The collision proof, greedy design, entropy conditioning, exact algorithm
and runtime actually used above are all supplied here.

# B6. A robust function for almost independent bits

We supply a specialized bounded-independence argument below. No general
AC0 fooling theorem or non-oblivious bit-fixing extractor is assumed.
This proof uses the partition architecture of Chattopadhyay--Zuckerman,
with B5's much smaller shifted agreement and an explicit two-level
Bonferroni comparison. Only constant separation from zero and one is
needed for dispersal.

## Exact layout and bias parameter

Let the available bit string have length M0=2^d0. Choose the largest integer
r with

 m=floor(delta*r), b=floor(floor(sqrt(S_header))*ceil(log2 r)^2/(30*10^24)),
 u=2b+2ceil(log2 r)+10,
 u*2^(m+b)<=M0.                                         (B6.1)

S_header is the actual fixed input scale for this hardness argument.
The implementation finds r by integer binary search, using
ceil(log2 u)+m+b<=d0, an equivalent comparison. Put R=2^r, B=2^b,
M=2^m, S=BM. Use the first u*S input bits as S consecutive u-bit groups.
For EVERY fixed header S_header, b=Theta(log(r)^2)=o(r), and
u=Theta(log(r)^2). Hence the maximal layout has r=Theta(log M0).
At an adjacent integer r, m increases by at most one. The ceil-logarithm
w increases by at most one, and b increases by at most C_S*w, for a
fixed constant C_S. The adjacent ratio of u is at most two eventually.
Maximality therefore implies

 M0*2^(-C'_S*log r) < u*B*M <= M0.                   (B6.1a)

A constant-fraction coverage claim is NOT used. Instead the uncovered
factor is 2^poly(log r), subpolynomial in M0. We have

 b/m -> 0, r=Theta(log M0), u=Theta(log(r)^2).

For the sufficiently large fixed header of B7, B5's seed fits b.
All computation and the displayed comparisons use exact integers.

Choose an INTEGER J by bisection, as the largest 0<=J<2^u satisfying

 R*(2^(uB)-(2^u-J)^B)^M <= 2^(uBM).                      (B6.2)

Define g on a group to be one iff its positively weighted binary integer
is at least J. This is monotone in its bits. On uniform group bits let
p1=J/2^u, p2=(1-p1)^B, p3=(1-p2)^M, and alpha=Rp3.
Then

 alpha -> 1,    Mp2~ln R,    Bp1~ln M.                  (B6.3)

Here is a check that bisection hides no real constant or accuracy premise.
The continuous root with Rp3=1 has
p2=1-R^(-1/M)~ln R/M and p1~(ln M-ln ln R)/B.
The derivative of ln p3 with respect to p1 near that root is O(B ln R).
The grid step 2^-u is O(1/(B^2*r^2)), so its product with that derivative
tends to zero. The root lies between consecutive grid points chosen in
(B6.2), proving alpha->1 from below. These estimates also give 0<J<2^u
and the other two comparisons eventually. The implemented comparisons
are exact integer exponentiations on at most uBM+O(r) bits. Since
uBM<=M0, they require polynomial time in M0.

For each v in F_2^r the B seed outputs z_j=G(v,j) define a partition of
the S groups into M blocks

 P_v,w={(j,z_j+w): j in F_2^b},  w in F_2^m.

The group indexed (j,a) has array index j*M+a; addition is XOR. Each block
has size B and these blocks partition the groups. Blocks from DIFFERENT
partitions meet in at most eta*B groups by B5(ii), for ANY two shifts.
Define the monotone function

 H=AND_v OR_w AND_(group in P_v,w) g(group bits).         (B6.4)

Its loops have at most R*M*B iterations, with each group map computable in polynomial time in r. G includes
B5's one-time polynomial(r) weak-design construction and polynomial(r)
work per value. Since r=O(log M0), the total is polynomial time in M0. Its polynomial exponent is very
large; no modest-length efficiency assertion is made.

## Two elementary probability tools

Bonferroni follows by expanding product(1-I_i) and truncating the
alternating sum. For any finite event family, truncations of consecutive
parities bound the probability that none occurs; their gap is the
corresponding factorial moment sum S_t=sum_(|T|=t)Pr[all events in T].
If the events are independent with mean count lambda, S_t<=lambda^t/t!.
We also use t!>=(t/e)^t, obtained for example by integrating log x.

We record and prove the needed multiplicative form of Janson's inequality.
For conjunction events E_i in independent Bernoulli coordinates, let
p_i=Pr E_i<=tau<1, and let Delta=sum_(i<j, supports intersect) Pr[E_i E_j].
Then

 product_i(1-p_i) <= Pr[no E_i]
 <= product_i(1-p_i)*exp(Delta/(1-tau)).                 (B6.5)

Harris association is elementary here: induct on the last independent
coordinate; conditional covariances of increasing functions are nonnegative,
and their conditional expectations are increasing in that last coordinate.
It gives positive association of decreasing events too, proving the lower
bound. For the upper bound expose E_i in order. Write A for nonoccurrence
of earlier NONneighbors (their supports are disjoint from E_i), and B'
for nonoccurrence of earlier neighbors. A is independent of E_i. Since
E_i E_j is increasing and A decreasing, Harris gives
Pr[E_i E_j|A]<=Pr[E_i E_j]. Hence

 Pr[E_i|A B'] >= Pr[E_i B'|A]
              >= p_i-sum_(earlier neighbors j) Pr[E_i E_j].

Multiply the resulting upper bounds on conditional nonoccurrence
probabilities, and use 1+x<=exp x. This gives (B6.5). This proof also
applies to conjunctions of independent group variables with the same bias.

## Arbitrary bad coordinates

Assume at most q of the original M0 bits are bad, and the remaining bits
have the following almost independence: every set of at most t of them
is within gamma of uniform. The bad bits may have arbitrary joint dependence on the good bits,
including additional conditional randomness. It is enough that, eventually,

 q<=S^(1-delta), t>=A_t*r^2*u, gamma<=2^(-A_gamma*r^2*u).     (B6.6)

Here A_t,A_gamma are fixed finite positive constants, depending only on
delta and the fixed outer truncation order 100. Their existence is PROVED
by the finite degree and coefficient estimates below, not assumed as an
external extractor premise. B7 explains how the fixed total language
contains an actual input header large enough to meet all such constants.

Call a group bad if it contains any bad bit. There are at most q bad
groups. Define H0,H1 by setting ALL bits of bad groups to zero or one,
respectively. Monotonicity gives H0<=H<=H1 pointwise. The good groups
retain their asserted almost independence. On ideal independent good
bits, each good group has failure probability p1.

Call a partition exceptional if one of its blocks has more than
3epsilon*B bad groups, where epsilon is the strong error of B5. There
are fewer than 2^(2m)*M=2^(3m)<=R^(3delta) exceptional partitions.
Indeed, otherwise some shift w has at least 2^(2m) such v. Take their
uniform source of min-entropy 2m and apply the strong extractor G, keeping
the seed in the output. Its probability of landing in the bad-group set
exceeds 3epsilon, whereas the uniform probability is q/S<=epsilon
for large r; this contradicts strong error epsilon. The latter inequality
uses q/S<=S^-delta, which is exponentially small in r, whereas
epsilon=r^-80.

We will ignore exceptional partitions temporarily. Let F_v be the event
that partition v's OR-of-ANDs is zero. For any partition, H0 has at least
M-q wholly unaffected disjoint blocks. Under the independent model

 Pr[F_v for H0] <= (1-p2)^(M-q)=(1+o(1))/R.              (B6.7)

The last identity uses q/M->0 exponentially: from (B6.1),(B6.6),
log2(q/M)<= (1-delta)b-delta*m=-Theta(delta^2*r).
The same bound, with an additional negligible term, will be proved below
under the almost-independent good-bit distribution. Since H1's event
F_v is contained in H0's, discarding all exceptional partitions changes
either output probability by o(1), using their R^(3delta) count.

## Ideal independent probabilities on the remaining partitions

Fix any CONSTANT c<=100 distinct nonexceptional partitions and either
bad-group setting. For H0 discard blocks made impossible by a zero group;
at most q per partition are discarded. Other block events have probability
p2. For H1 each block loses at most 3epsilon*B random groups, and its
all-one probability is p2*(1+o(1)), uniformly. The product of individual
nonoccurrence probabilities over the c partitions is

 D=(p3)^c*(1+o(1)).

For H0 its logarithmic discrepancy is O(c*q*p2)=o(1). For H1 it is
O(c*epsilon*r^2)=o(1), using (B6.3). A block meets at most B blocks of
each other partition. The number of dependent pairs is at most c^2*M*B.
Their joint probabilities are at most
p2^(2-6epsilon-eta)=O(p2^2), since (epsilon+eta)r->0. Thus

 Delta=O(c^2*B*r^2/M)=o(1).

The ratio B/M tends to zero exponentially in r. Applying the supplied
inequality (B6.5) proves, uniformly in these partitions and the bad set,

 Pr[intersection of these c F_v]=(p3)^c*(1+o(1)).        (B6.8)

There are R-o(R) retained partitions. Their outer factorial moment sums
of each fixed order c<=100 tend to 1/c!, by alpha->1. The two outer
Bonferroni truncations at consecutive orders 99 and 100 consequently
place the retained function's independent-model expectation strictly
between .36 and .38 for sufficiently large r. These bounds follow also
by the elementary alternating series for exp(-1); the truncation gap
is at most 1/99!, far below .001.

## Comparison under limited independence

We give the details because pairwise independence of the good bits would
not suffice. Throughout, c is fixed, at most 100. Constants below may
depend on c and delta, but never on r, M0 or the chosen bad set.

Consider the no-block-event probability for these c partitions. For j
DISTINCT block events, j<=C*r, their union contains at least

 jB(1-3epsilon)-j^2*eta*B/2

random good groups; impossible H0 events simply have probability zero.
Using independent group bits, the joint block probability is at most
2*p2^j for large r. Therefore its factorial moment of order j is at most
2*(c*M*p2)^j/j!. Its mean parameter cMp2 is O(r). Choose an EVEN
J_inner=C*r, rounded up, with C sufficiently large that the moment at
that order is at most exp(-H*r), where H is any prescribed fixed constant.
The factorial estimate above permits an arbitrarily large H by increasing C.

To compare a joint block event under the two distributions, view it as
NO group failures on a union of at most jB groups. In the independent
model the failure count has mean at most jB*p1=O(r^2). Truncate its
Bonferroni expansion at an EVEN K_inner=C'*r^2. Choosing C' sufficiently
large makes the independent factorial remainder at most exp(-H'*r^2),
for any prescribed H'. The same remainder under the other distribution
can be compared termwise. Each term involves at most K_inner good groups,
hence at most K_inner*u=O(r^2*u) good underlying bits. Its expectation
differs by at most gamma by (B6.6). The total number of terms is at most
exp(O(r^2*u)), since log(jB)=O(u) and K_inner=O(r^2). Consequently EACH joint block probability differs by at most

 2exp(-H'*r^2)+gamma*exp(O(r^2*u)).

There are at most (cM)^j=exp(O(r^2)) block tuples of order j<=J_inner.
Choose H' large enough before using the estimate. The differences of all
block factorial moments through J_inner are then exp(-H''r^2) plus
an error gamma*exp(O(r^2*u)). Apply the block-level Bonferroni bounds to
both distributions; their gap is controlled by the J_inner moment just
bounded. The probability of the c-fold F_v intersection differs by at
most exp(-H*r)+o(exp(-H*r)).

All required bit degrees are at most C_deg*r^2*u, and the coefficient
counts are at most exp(C_count*r^2*u), for FIXED finite constants after
J_inner and K_inner above have been chosen. These constants are uniform
in the selected bad set and c<=100. Set A_t>=C_deg and choose
A_gamma*ln(2)>C_count+1. Then gamma*exp(C_count*r^2*u)<=exp(-r^2*u),
which is smaller than exp(-H*r^2) for every fixed H eventually. This
proves the existence and adequacy of the constants in (B6.6). They do
not depend on N, an input header, or the particular source. Choose H>1000. After weighting by
at most R^c=2^(cr) choices, every outer factorial-moment discrepancy
through order 100 tends to zero. Thus the same .36--.38 expectation
bounds hold under the actual good-bit distribution for the retained
H0 and H1.

The one-partition comparison for (B6.7) is identical with its impossible
H0 blocks discarded. All remaining blocks are disjoint and wholly good,
so their factorial moments have the same O(r) mean bound. It yields
Pr[F_v for H0]<2/R for large r, plus an error exp(-1000r). Union over
R^(3delta) exceptional partitions is o(1). Restoring those partitions
therefore puts BOTH E[H0] and E[H1] between .35 and .4 eventually.
Since H0<=H<=H1, both values of H have positive probability, regardless
of how the bad bits depend on the good ones. This proves the robust
non-oblivious bit-fixing disperser needed here.

# B7. Polylogarithmic design seed and a fixed scale header

## Fixed algorithm on every header

Let N be the ACTUAL input length, L=ceil(log2 N), z=ceil(log2 L).
The first L input bits encode S-1, so 1<=S<=2^L<=2N. The algorithm
evaluates every S without knowing or searching for a good scale. Let
P2(t) be the least power of two at least positive integer t. The exact power-of-two parameters are

 k=P2(S^13*L^15*z^16), h=e=P2(S*L^2*z^2),
 ell=P2(S^3*L^4*z^4), e'=P2(S^5*L^6*z^6),
 d=P2(S^9*L^10*z^12), Y=P2(S^11*L^13*z^14),
 a=2, b=d0.                                              (B7.1)

All operations are exact integer arithmetic. Initial A1 dimension is 10000 FIXED, E1 has A1 dimension z, and the other four internal maps use B2.2, and B5 uses its supplied weak-design recycling
map rather than a design A1 package. The raw log2 k is checked before expansion. Rejected large-length
plans return zero; only the finite cutoff uses the legacy evaluator.
Bad field data on an admitted plan still give bounded ring/matrix evaluation.

Eligibility is proved for the ONE fixed algorithm on every header. The
lower-bound argument separately proves that one fixed integer S0 makes a
hard actual-input restriction at every sufficiently large length. Header
and data bits are charged; S0 is not an external input, advice or witness.

## Uniform estimates for each fixed scale

The four half-output maps use B2.2 with seed cost at most C_seed*e^2 plus lower terms. For the single A1-based internal map E1, D=z, H=O(z),R=O(z), and field degree O(eprime). Its short one-call advice is at most C_advice*ell*z^2*eprime plus lower terms, and its seed cost is O(z*eprime+log N)=o(d). The implied constants are independent of fixed S once N is large enough, since log S/z->0.

Initial output Y is much smaller than k. Its single-call advice has order
Y*k^(.0001)*L=o(k): its L exponent is 14.0015<15. Its field degree is at most
(128+o(1))L and d0<=10^9 L eventually for every fixed S. This cap is
checked before seed enumeration. Since the alphabet is BINARY, b=d0,
so b<=B0*L with B0=2*10^9; no divisor log h is claimed.

Power-of-two rounding costs at most two. The key comparisons are

 e^2 <=4S^2L^4z^4,       ell>=S^3L^4z^4,
 ell*z^2*e' <=4S^8L^10z^12,  d>=S^9L^10z^12,
 h*a*ell <=8S^4L^6z^6,     e'>=S^5L^6z^6,
 bhd <=4B0S^10L^13z^14,     Y>=S^11L^13z^14,
 hY <=4S^12L^15z^16,        k>=S^13L^15z^16.                (B7.2)

Each required constant ratio gains a power of S. A sufficiently large
fixed S meets all seed, advice, weak-seed and row-rank requirements. For
example S>=3200 leaves weak-seed room, S>=8000B0 leaves Y>=1000bhd
room, and S>=16 gives hY<=k/4. The other source loss bh*a*ell has
smaller L power than k, as does10d. Also q0=O(h*a*ell)=o(d), and2e=o(e').
The code retains the actual checks with h*a, not a mistaken h^2 alphabet.

## One fixed hard scale

The layout still has r<(d0+1)*10^12<=c_r L, c_r=2*10^21.
For every fixed S, eventually ceil(log2 r)<=2z and
u<=C0*sqrt(S)*z^2, for an absolute finite C0. B6 requires

 t>=A_t*r^2*u, gamma<=2^(-A_gamma*r^2*u).

Thus h=e>=S*L^2*z^2 meets both requirements, including the polynomial
prefactor in100bh^2*2^-e, once sqrt(S) exceeds the corresponding fixed
constants times C0*c_r^2. B5's seed also fits once sqrt(S) exceeds its
fixed C_design bound. Choose ONE integer S0 exceeding those finitely
many bounds and the seed/advice/coefficient bounds in (B7.2).
Executable checks h>=r^2*u,e>=100r^2*u are resource checks, not assertions
that every admitted header certifies dispersion; this S0 meets the full
proved constants eventually.

B6.1a gives S_groups>=M0*2^(-O_S0(log r))/u. Since r=O(L),
the initial bad-coordinate bound q<=2M0/N^2 implies

 q/S_groups<=N^(-2+o(1)), S_groups^delta<=2^(.001L).

Their product tends to zero. Thus q<=S_groups^(1-delta) eventually.
Every B6 hypothesis holds. This ONE S0 fits the header once2^L>=S0,
and every ready guard holds eventually, for ALL sufficiently large
actual integer lengths.

This is a lower-bound restriction, not a nonconstructive definition of
the language. The evaluator already computes all header values. It never
receives S0 externally. The same algorithm and language work at every
length and have no good-header promise.

## Actual prefix and uniform computation

The first L bits are the header, followed by A1 packages initial and E1, then B5's hash modulus. The other four internal maps use deterministic tower fields and have zero data bits. An A1 package has exactly (D^2+1)s bits at its stated dimension. After reading the header the engine computes
all lengths, starts its field cursor AFTER that header, and removes the
ENTIRE header-plus-data prefix before padding the payload to N coordinates.

For fixed S0, the largest internal precision is e'=Theta(S0^5 L^6z^6).
Including every matrix, modulus and the header, the exact computed prefix is

 P(N,S0)=L+O(S0^5 L^6z^8)=o(k(N,S0)).

The large condenser modulus is still computed deterministically. One good
tuple works simultaneously for all affine sources at the fixed scale and
length. Zero padding injects every k-flat on N-P actual payload bits into
the N-coordinate source. The payload restriction is an affine disperser.
All its actual input bits are charged and P+k<=N eventually.

Uniform eligibility also covers all other headers. S<=2N, and k>N plans
reject early. Every remaining scalar is polynomially bounded in N. The
component/source checks and total-prefix guard bound admitted fields and
matrices. The d0 cap bounds all seed enumeration; r=O(log N) bounds
robust-function loops by fixed polynomials. Internal D=z introduces no
N^D loop: only the first n' coefficient monomials are streamed with
polynomial delay. B5's weak design is explicitly computed with exact
conditional expectations in polynomial(r) time, not supplied as advice. Field degrees and all bit operations are polynomially
bounded. Rejected large plans give zero. Small finite cutoff lengths use
the fixed halting legacy routine. Thus eligibility does not depend on
finding S0, or on any per-length advice or unproved efficient construction.

## Provenance

A1 supplies the Shaltiel--Umans reconstruction, B1 the GUV graph/binary
linearity and deterministic tower. B4 follows Li's index-based affine
architecture with a supplied binary alphabet, reduced prephase transcript,
correct A/B conditioning and every-round refresh. B6 follows the
Chattopadhyay--Zuckerman partition architecture with the supplied shifted-
agreement generator and moment comparison. No modern extractor theorem
or prior acceptance is a premise.

Primary provenance: [GUV](https://people.eecs.berkeley.edu/~venkatg/pubs/papers/PV-condenser.pdf),
[Li2015](https://eccc.weizmann.ac.il/report/2015/121/download), and
[Chattopadhyay--Zuckerman](https://annals.math.princeton.edu/wp-content/uploads/annals-v189-n3-p01-s.pdf).
The finite full-B2 gate proof is separately supplied in elimination.md.
