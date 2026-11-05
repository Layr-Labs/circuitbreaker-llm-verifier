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

# B2. Streamed strong linear extraction

An affine distribution has a vector-space complement description: if Z is
an affine function of X, then X=A+B, with A an affine function of Z and
B an independent uniform linear source. Its dimension is H(X|Z).
Conditioning on a fixed linear observation of t bits reduces any remaining
affine rank by at most t. An observation determined by Z does not reduce
H(X|Z). These statements follow by choosing a basis of the projection
kernel and extending it to a basis of the direction space.

Consider the following extractor plan for source length n, entropy K,
output M<=K/2 and error precision 2^-e. First apply B1 with precision
e+10. On every good first seed its output is affine of entropy at least
K. For the resulting ambient length n'=O(K+log n+e), use A1 with D=10000
and per-call error 2^-(e+ceil(log2 K)+10). Choose the largest power-of-two
output m0<=min(M,K/4) satisfying A1's explicit entropy condition with
threshold floor(K/2). The executable plan halves m0 until it succeeds;
failed finite parameter cases are rejected by the outer ready flag.

Repeat this A1 extractor on the SAME condensed input with independent
fresh seeds until M bits have been produced, truncating the last call.
After fixing earlier seeds and outputs the source remains affine, of
entropy at least K-M>=K/2. A sequential coupling thus loses at most the
sum of the call errors. There are at most K calls. Adding the B1 bad-seed
probability gives total strong error at most 2^-e. Every fixed combined
seed produces a linear map of the original input. Padding a prescribed
seed length with ignored bits preserves both statements.

For completeness the seed budget is not assumed optimal without proof.
Write H=10000*ceil((n')^(1/10000)) and R=160016. A1's field degree is
O(e+log K+log H). The largest admissible m0 is at least a fixed constant
multiple of K/[R*H*(e+log K+log H)], unless the desired M is smaller,
in which case one call suffices. This follows by substituting that value
in A1's displayed inequalities; the additive logarithmic term is o(K)
in our uses. Thus the total seed length is bounded by

 C(log n+e+log K)+C(n')^(1/10000)*(e+log K+log n')^2,       (B2.1)

for a fixed absolute C. All uses below have K>>log n+e and so n'=O(K).
This also proves that the halving loop succeeds eventually.

There is one initial use with output Y and precision 2^-(4L+10),
where L=ceil(log2 N), n=N, and k,Y,T are specified in B7. It needs just one
A1 call eventually: its entropy advice is
T^(51+.0056+1+o(1))=o(T^56)=o(k). Its field degree is O(L), and its
combined seed length d0 is at most 10^9 L eventually. The generous
constant follows directly from A1: the field-degree leading term is at
most 128L, multiplied by D+1=10001, plus O(L) condensation bits. The
finite plan checks this cap before enumerating seeds.

At fixed good construction data, an affine output either is exactly
uniform or has distance at least 1/2 from uniform. A strong error bound
eta therefore implies a bad-seed fraction at most 2eta. If a seed of length
a has min-entropy at least a-t and is independent of the affine source,
its probability of hitting this bad set is at most 2^(t+1)*eta. We use
this elementary rank observation, rather than treating a nearly uniform
output as a new affine source without justification.

# B3. Parallel affine look-ahead extraction

Let at most h streams Q_i be affine functions of Z, each padded to a common
ambient bound. One designated Q_v is exactly uniform on q0=(2h^2+10)ell
bits. Let H(X|Z) be at least h^2*ell+10ell. The extractors W on X and T
on Q_i have output and seed length ell, source entropy threshold 10ell,
and strong error eta. Set

 S_i,1=first ell bits of Q_i,
 R_i,t=W(X,S_i,t),
 S_i,t+1=T(Q_i,R_i,t).

For t<=h, the designated R_v,t is uniform conditional on all Q_i, all
current and earlier S values, and all earlier R values, except on total
history mass at most 4h*eta. More precisely, after the preceding messages
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
by at most h*ell. Before the last needed transition their remaining
ranks exceed 10ell by the hypotheses. The initial uniform seed and the
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
Let b=ceil(d0/log2 h). Write i as a d0-bit word and append zero bits at
the RIGHT to make b base-h digits. Write its digits plus one as
 t_i,1,...,t_i,b in [h]. This padded representation is injective and
preserves numerical order.

For each i start D_i,1=first d bits of Y_i. On every round j compute

 R0_i,j=E1(X,D_i,j),
 Q_i,j=E2(Y_i,R0_i,j),
 (R_i,j,1,...,R_i,j,h)=look-ahead(X,Q_i,j),
 V_i,j=R_i,j,t_i,j.

For j<b put D_i,j+1=E3(Y_i,V_i,j). The row W_i is V_i,b. It is enough
for the evaluator to compute each look-ahead stream only up to its selected
digit. The full prefix messages used in the proof are its identical
mathematical continuation and do not change the selected output.

The source/seed/output specifications are:

* E1: ambient N, entropy 10d, output ell, seed d, error eta'=2^-e'.
* E2: ambient Y, entropy (4h^2+20)ell, output q0, seed ell,
  error eta=2^-e.
* E3: ambient Y, entropy 2d, output d, seed ell, error eta.
* W: ambient N, entropy 10ell, output ell, seed ell, error eta.
* T: ambient d, entropy 10ell, output ell, seed ell, error eta.

Short Q_i are padded at the END before T is evaluated. Construction B2
provides all five maps with these seed lengths, as verified in B7.
Assume the following conservative bounds, also checked eventually there:

 Y>=1000*b*h*d,
 k-h*Y-100*b*h^2*ell>=10d,
 d>=(2h^2+10)ell,
 e'>=100h^3*ell+2e.                                      (B4.1)

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
 H(X|Z)>=k-h*Y-100(j-1)h^2*ell.                     (B4.2)

Initially this follows from uniformity of Y_v. The additional transcript
values described below make the induction stronger, not weaker: whenever
an A-part is fixed it is a linear function of Z, and whenever a seeded
output is fixed its seed has already been fixed.

First process the already-diverged rows E_(j-1). Their current seeds are
fixed, so fix their R0, then Q, then all h look-ahead R messages in order.
All these fixings are affine; the S messages are determined. They cost at
most 100h^3*ell bits from D_v,j, and at most O(h^2*ell) from H(X|Z).
In particular D_v,j now has entropy at least d-100h^3*ell. No same-prefix
seed of length d has been fixed at this point. This distinction is essential.

Split X=A+B by Z in the resulting affine component. The target seed is a
function of A and independent of B. Its weak-seed probability of choosing
an E1 map of deficient rank on B is at most
2^(100h^3*ell+1)*eta'<=eta, by (B4.1). Fix D_v,j first. On a good value
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
current digit t. Outside mass at most 4h*eta, its final B-part is uniform
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
and look-ahead A-parts cost O(h^3 ell)=o(hd). Only O(h^2 ell) output bits
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

# B5. A strong extractor with very small shifted agreement

This section constructs the exact partition generator needed in B6, rather
than invoking Trevisan extraction or a resilient-function black box.
Let delta=10^-12. For a large integer r put m=floor(delta*r), and later
pad a seed to b~delta^2*r/30 bits. We construct G:{0,1}^r x {0,1}^b ->
{0,1}^m with the following two properties:

 (i) G is strong on every source of min-entropy 2m, with error
     epsilon<=2^(-r^.08).
 (ii) For every distinct v,w and every offset a in F_2^m,
      Pr_z[G(v,z)+G(w,z)=a]<=eta<=2^(-r^.09).

All bounds are eventual in r. The mathematical hypotheses on the A1 field
and traversing-matrix data are supplied by the actual charged prefix.

Set s=100ceil(r^.1), c=ceil(s/4), and take a field F_(2^s). Encode v
injectively as a polynomial P_v over that field by packing successive s-bit
coefficients. For a pair of independent uniform field seeds y,z put

 C(v;y,z)=first c binary coordinates of z*P_v(y).

For v!=0, the polynomial P_v is nonzero of degree below ceil(r/s)<=r.
At at most r values of y it vanishes; at all other y multiplication by
P_v(y) makes z uniform in the field. Consequently for every c-bit a,

 Pr[C(v;y,z)=a]<=2^-c+r/2^s.                              (B5.1)

This proves (ii) for any map whose output starts with C, since C is linear
in v for each seed. It also makes C a strong extractor on a 2m-source.
Indeed the average conditional collision probability for C(X;y,z) is
at most 2^-c+2^-2m+r/2^s. Cauchy--Schwarz gives strong distance at most

 (1/2)*sqrt(2^c*(2^-2m+r/2^s)).                           (B5.2)

This collision calculation explicitly includes the seed, not just the
output marginal.

Append m-c output bits obtained by fresh independent A1 seeds on the SAME
r-bit input. Use D=100, per-call output the largest admissible power of two
at most r^.875, and per-call error 2^-(ceil(r^.1)+ceil(log2 r)+10).
The explicit entropy threshold is floor(m/2). Its monomial degree is
O(r^.011) (the input ambient is r), and its field degree is O(r^.1).
Thus the advice exponent is at most .875+.011+.1=.986<1; it fits m/2
for every fixed positive delta eventually. There are O(r^.125) calls,
and their total seed length is O(r^.225). Including the 2s hash bits,
this is o(r), hence it fits b~delta^2*r/30 eventually. Any unused seed
bits are ignored. All output bits remain linear in the input for fixed
whole seed.

Here A1 must work for general min-entropy sources, not only affine ones;
its proof does. The conditioning issue in repetition is handled explicitly.
For fixed earlier seeds, an earlier output prefix of at most m bits has
probability below 2^(-m-r^.1) on a set of total mass at most 2^-r^.1.
The same argument after fixing the c hash bits uses a threshold with c
additional bits. Outside those exceptional fixings, conditioning on the
hash and earlier outputs leaves min-entropy at least
2m-c-(m-c)-O(r^.1)>=m-O(r^.1)>=m/2. Fresh next seeds are independent.
Thus a sequential coupling costs at most a polynomial factor times
2^-r^.1, in addition to (B5.2). This is at most 2^-r^.08 eventually.
It proves (i) for the appended output jointly with the hash and all seeds.
The sharper shifted-agreement estimate (B5.1) is at most 2^-r^.09
because c~25r^.1 and s~100r^.1. The constants .08 and .09 leave room
for every integer rounding and polynomial factor.

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

 m=floor(delta*r), b=floor(delta^2*r/30),
 u=2b+2ceil(log2 r)+10,
 u*2^(m+b)<=M0.                                         (B6.1)

The implementation finds r by integer binary search, using
ceil(log2 u)+m+b<=d0, an equivalent comparison. Put R=2^r, B=2^b,
M=2^m, S=BM. Use the first u*S input bits as S consecutive u-bit groups.
One step of r changes 2^(m+b) by a factor at most four, and changes u
by a factor tending to one. Thus uS>M0/5 eventually, while uS<=M0;
in particular S>=M0/poly(log M0). We have

 b/m -> delta/30,   r=Theta(log M0),   u=Theta(r).

The seed length required in B5 is less than b eventually. All constants
in these asymptotic statements depend only on the fixed delta.

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

Its loops have at most R*M*B iterations, with each group map and G
computable in polynomial time in r. Since r=O(log M0) with a fixed
constant, H is polynomial time in M0. Its polynomial exponent is very
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
is within gamma of uniform. The bad bits may be arbitrary functions of
the good bits. It is enough that, eventually,

 q<=S^(1-delta),  t>=r^(7/2),  gamma<=2^(-r^(7/2)).             (B6.6)

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
uses q/S<=S^-delta and epsilon decaying only stretched exponentially in r.

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
hence at most K_inner*u=O(r^3) good underlying bits. Its expectation
differs by at most gamma by (B6.6). The total number of terms is at most
exp(O(r^3)). Consequently EACH joint block probability differs by at most

 2exp(-H'*r^2)+gamma*exp(O(r^3)).

There are at most (cM)^j=exp(O(r^2)) block tuples of order j<=J_inner.
Choose H' large enough before using the estimate. The differences of all
block factorial moments through J_inner are then exp(-H''r^2) plus
an error gamma*exp(O(r^3)). Apply the block-level Bonferroni bounds to
both distributions; their gap is controlled by the J_inner moment just
bounded. The probability of the c-fold F_v intersection differs by at
most exp(-H*r)+o(exp(-H*r)).

All required bit degrees are O(r^3)<r^(7/2), and all coefficient counts
are exp(O(r^3)), so the two explicit conditions on t,gamma in (B6.6)
justify every comparison, eventually. Choose H>1000. After weighting by
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

# B7. Exact global parameters and all-length coverage

Let N be the ACTUAL input length, L=ceil(log2 N), z=ceil(log2 L), T=2^z.
The exact integer scales are

 k=T^56, h=T^4, ell=T^9, d=T^45,
 Y=T^51, e=T^4, e'=T^22.                                (B7.1)

Y is the initial row width, chosen with the explicit source-rank and
row-rank budgets of B4.
All powers and comparisons are integer operations. Initial extraction and
the five maps are exactly B2 and B4 with their stated dimensions. The
algorithm compares 56z with floor(log2 N) before allocating k. All field
degrees and halving choices use A1's explicit formulas.

The initial B2 extraction has output Y=T^51 and source entropy k=T^56.
Eventually it uses one A1 call: its advice has scale
Y*k^(1/10000)*L=T^(52.0056+o(1))=o(k). Its source still has entropy at
least k/2 for that call. The field-degree leading bound in B2 is unchanged:
error precision is 4L+10, while log Y=O(log L)=o(L). Thus d0<=10^9 L
eventually, and the guard applies before seed enumeration.

For the remaining maps the condensed ambient is O(K), so (B2.1) yields:

* Q-side and X-side: K=10T^9, precision T^4, seed cost
  O(log N+T^8.0009+lower terms)=o(ell).
* E2: K=Theta(T^17), precision T^4, seed cost
  T^(8.0017+o(1))=o(ell).
* E3: K=2T^45, precision T^4, seed cost
  T^(8.0045+o(1))=o(ell).
* E1: K=10T^45, precision T^22, seed cost
  T^(44.0045+o(1))=o(d).

Hence every stipulated seed string can be padded with ignored bits. The
internal field seeds are evaluated without enumeration. We have
b=ceil(d0/log2 h)=T^(1+o(1)). The ACTUAL conditions of (B4.1) hold:

 bhd=T^(50+o(1))=o(Y),
 hY=T^55=o(k), bh^2ell=T^(18+o(1))=o(k), d=T^45=o(k),
 e'=T^22 >>100h^3ell=100T^21 and 2e=2T^4,
 d=T^45 >>q0=Theta(T^17).

The checks are Y>=1000bhd and k-hY-100bh^2ell>=10d, plus the stated source,
seed, output and precision checks. They prove every entropy budget in
B3--B4 eventually.

For the robust layout r=Theta(d0)=O(L), with fixed possibly enormous
constants depending on delta. The evaluator checks

 h>=ceil(sqrt(r^7)),  e>=100ceil(sqrt(r^7)).

Both hold eventually because h=e=T^4 and r=O(T). The row joint error
100bh^2*2^-e is at most 2^(-r^(7/2)) eventually. This controls every
O(r^3)-degree observation and exp(O(r^3)) coefficient sum in B6.
There are at most q=2M0/N^2 bad first-bit rows. Since uS>M0/5,
u=O(L), and d0<=10^9 L,

 q/S<=10u/N^2,  S^delta<=M0^delta<=2^(.001L).

Their product tends to zero, giving q<=S^(1-delta) eventually. All B6
conditions hold, so the final output takes both values on every affine
input of entropy k.

## Prefix accounting

Only A1's moduli and traversing matrices and B5's hash modulus are read
from the ACTUAL input prefix. B1's large extension modulus is computed
deterministically. Prefix order is the six A1 packages (initial,Q-side,
X-side,E1,E2,E3), then the hash modulus and its A1 package. Each package
has exactly (D^2+1)s bits, with D=10000 in B2 and D=100 in B5. Repeated
calls share their construction data and use independent seed substrings.

The largest field-degree scale is E1's T^22. Other precisions are T^4,
L=Theta(T), and design O(r^.1+log r). Therefore the ENTIRE prefix has
its exact computed integer length with asymptotic

 P=T^(22+o(1))=k^(11/28+o(1))<k^.4.

In particular P=o(k), and P+k<=N eventually. Good tuples exist by A1's
field/traversing-matrix argument and B5's field existence. One tuple works
simultaneously for all affine sources at a given length. It is not chosen
as uncharged advice for each source. The payload is zero-padded injectively,
so every k-flat on the N-P genuine payload coordinates remains a k-flat
in the N-coordinate ambient source. The restricted payload function is
an affine disperser.

All guards pass as N ranges over EVERY sufficiently large integer:
L<=T<2L introduces only fixed factors, while every comparison above has
a strict exponent gap. Failed finite lengths use the fixed legacy evaluator.
Invalid data at admitted lengths still give bounded ring/matrix evaluation,
not a promised-input language. Uniform total polynomial time is detailed
in proof.md.

## Provenance and distinction of results

A1 supplies the Shaltiel--Umans reconstruction argument. B1 supplies the
Guruswami--Umans--Vadhan graph argument with binary linearity and a
supplied deterministic tower representation. B4 follows Li's affine
alternating-extraction architecture with the conditioning order and
last-round refresh proved above. B6 follows the Chattopadhyay--Zuckerman
partition architecture with the supplied shifted-agreement map and
specialized moment comparison. These supplied components and proofs,
not earlier acceptance or a modern extractor black box, are the premises.

Primary provenance: [GUV](https://people.eecs.berkeley.edu/~venkatg/pubs/papers/PV-condenser.pdf),
[Li2015](https://eccc.weizmann.ac.il/report/2015/121/download), and
[Chattopadhyay--Zuckerman](https://annals.math.princeton.edu/wp-content/uploads/annals-v189-n3-p01-s.pdf).
The full-B2 finite gate proof is separately supplied in elimination.md.
