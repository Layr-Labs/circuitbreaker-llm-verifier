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



