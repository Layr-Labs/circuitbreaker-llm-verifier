# Supplied algebraic and extraction arguments

Sections 1--9 retain the earlier symmetric-polynomial construction and its
image-polynomial preliminaries. The main family's large-length branch uses
A1--A6 below, together with those proved preliminaries. Small lengths use
the legacy evaluator. All construction data are actual prefix bits.

# Algebraic disperser argument for the payload

This supplies the BK construction argument used in proof.md, specialized to
F_2 and prime field degree. We use k for field degree, r for block count,
and t for the symmetric-polynomial order. All subspace dimensions are over
F_2. This file establishes the needed nonstandard algebra rather than
assuming a cited construction is a certified black box.

**Lemma.** Let K be the field of size 2^k with k prime. Suppose 1<=t<r and
k*t>=r(r+1). For any nonzero linear pi:K->F_2 put w_i=r+2-i and

    f(a_1,...,a_r)=pi(sum_(|I|=t) product_(i in I) a_i^(2^(w_i)-1)).

Then f is nonconstant on every affine subspace of K^r of dimension
strictly greater than kr/t+kt+r(r+1).

This is the binary case of BK's Theorem 7.3. We prove it below, including
the coefficient-isolation and block-selection statements it needs.

## 1. Finite-field and polynomial preliminaries

Let q=2^k. A polynomial function on K^r has a unique representative with
each variable degree <q, obtained by reducing modulo X_i^q-X_i. Uniqueness
follows by induction on the number of variables using the root bound
for a univariate polynomial. Thus a nonzero positive-degree coefficient
in this representative proves that the function is nonconstant.

A linearized polynomial is sum_(j=0)^(k-1) a_j X^(2^j). These are exactly
the F_2-linear maps K->K. Each is linear by Frobenius. The coefficient
space and the space of linear maps both have F_2-dimension k^2; the map
from coefficients to functions is injective by the root bound, so is
surjective as well.

Tr(z)=sum_(h=0)^(k-1) z^(2^h) takes values in F_2. It is not zero as a
function: its nonzero polynomial has degree 2^(k-1)<q. For mu!=0 the
functional z->Tr(mu*z) is nonzero, since multiplication by mu is a
bijection. The resulting linear map from mu in K to the dual of K is
injective between two k-dimensional spaces, hence bijective. Therefore
pi(z)=Tr(mu*z) for a unique nonzero mu.

Let wt(b) be the number of 1s in the binary expansion of an exponent b
in the reduced range 0,...,q-1. Multiplication and reduction cannot increase
the total number of power-of-two factors: replace two equal powers by
one doubled power, with cyclic wrap at k because X^(2^k)=X. Each such
carry reduces the number of factors. Consequently an exponent formed
from h power-of-two factors has wt at most h. If h<k and the final weight
is h, there were no carries: the chosen residues modulo k were all
distinct. Taking a Frobenius power only cyclically shifts residues.
This also proves the subadditivity of binary degree used below, including
after reduction modulo X^q-X.

## 2. Image polynomials and their zero patterns

For a linear subspace V of dimension d define P_V(X)=product_(v in V)(X-v).
This is a monic linearized polynomial of degree 2^d. For an induction proof,
start with P_{0}(X)=X; adjoining a vector b to U gives

    P_(U+<b>)(X)=P_U(X)*P_U(X+b)
                =P_U(X)^2+P_U(b)*P_U(X),

which is linearized. Its kernel on K is V, so its image V' has dimension
k-d. Let Q=P_{V'}. Then Q(P_V(X)) is monic of degree q and vanishes on K,
hence is X^q-X. Also P_V(Q(X))=X^q-X: composing their difference with
P_V and using P_V(X^q-X)=P_V(X)^q-P_V(X) proves the difference is zero
as a polynomial (a nonzero polynomial stays nonzero after nonconstant
polynomial substitution). It follows that Q(K) is contained in V; its
kernel dimension k-d shows its image dimension is d, so Q(K)=V.
Thus Q is a monic image polynomial of degree 2^(k-d), with 2^(k-d)
distinct roots. Write Q=sum a_j X^(2^j), with coefficients extended
periodically modulo k. We use this for d>=2.

There is no cyclic run of d zero coefficients. Otherwise a Frobenius
power of Q, reduced modulo X^q-X, moves that run to the final d
positions. The resulting nonzero polynomial would have degree at most
2^(k-d-1) and still have all 2^(k-d) roots of Q, a contradiction.
In particular a_0!=0 and a_(k-d)!=0; the latter is also monicity.
The interval k-d+1,...,k-1 is a zero run of length d-1.

For prime k and d>1, this is the **only** cyclic zero run of length d-1.
Indeed another such run can be shifted to the final positions by a
Frobenius power 2^h with 0<h<k. The shifted polynomial has degree at most
2^(k-d), has the same roots, and so is cQ for nonzero c (subtract leading
terms and use the root bound). Hence every z in V satisfies z^(2^h)=c*z.
Choose nonzero z0 in V; then every z/z0 is fixed by Frobenius^h and by
Frobenius^k. Since gcd(h,k)=1, it is fixed by Frobenius itself (use Bezout
on powers of this invertible automorphism), so z/z0 is in F_2.
That would force dim(V)<=1, a contradiction. Thus every zero gap other
than the final one has length at most d-2. This is the required
prime-degree subfield exclusion, with its proof.

## 3. A nonzero coefficient for independent subspaces

Let V_1,...,V_t be linear subspaces of K of dimensions d_i>1. Let
1<=e_i<d_i-1, h_i=e_i+1, delta_i=d_i-h_i>0. Suppose

    sum_i delta_i > k-max_i d_i+1.

For image polynomials Q_i of V_i we prove that the reduced polynomial

    Tr(mu * product_i Q_i(X_i)^(2^(h_i)-1))

has a nonzero coefficient at a monomial product_i X_i^(b_i) with
wt(b_i)=h_i. Notice h_i<d_i<=k, so h_i<k.

Write Q_i=sum_j a_(i,j) X^(2^j). Define its k-by-k matrix T_i by

    T_i[u,v]=a_(i,v-u)^(2^u),

with indices modulo k. This is the coefficient of X^(2^v) in Q_i(X)^(2^u).
For A_i={0,...,e_i}, choose a column set B_i of size h_i and let
b_i=sum_(v in B_i)2^v. Expanding the h_i linearized factors and using the
no-carry observation, the coefficient of this monomial is

    sum_(rho=0)^(k-1) mu^(2^rho)
        * product_i perm(T_i[A_i,B_i-rho])^(2^rho).             (1)

A full-weight monomial can only use each column once, so the sum over
bijections is exactly the permanent. Constants or carries have smaller
binary degree and cannot contribute. In characteristic two permanent
equals determinant, although the following argument only needs unique
triangular contributions.

We construct B_i so the rho=0 summand in (1) is nonzero and every other
summand has a zero factor. Reorder these independent sources so d_1 is
maximum. Let J be the first index with sum_(i<=J)delta_i>=k-d_1+1.
It exists by the assumption. For i<=J put ell_i=sum_(j<i)delta_j;
for i>J put ell_i=0. For i<=J, ell_i<=k-d_1<=k-d_i and there is no
wraparound before or at the first crossing. In fact the crossing sum
is at most k-d_1+(d_J-e_J-1)<=k-2.

For each source abbreviate d=d_i,e=e_i,ell=ell_i and a_j=a_(i,j).
If a_ell!=0, take the ordered columns

    B=(ell, k-d+1, ..., k-d+e).

If a_ell=0, take adjacent nonzero coefficients a_l,a_j with
l<ell<j<=k-d. Put z=j-l-1, v=min(z,e), and take

    B=(j, l+1,...,l+v, k-d+v+1,...,k-d+e).

The final interval is empty if v=e. The zero gap is not the final gap,
so z<=d-2. These columns are distinct. In the second case, if z<e then
j<=k-d separates the first columns from the final interval; if z>=e
there is no final interval. In both cases the minor T[A,B] is lower
triangular with nonzero diagonal. For its first block the diagonal uses
a_j followed by copies of a_l under Frobenius; entries above it are in
the chosen zero gap. The remaining diagonal uses a_(k-d), and entries
above it are in the final zero gap. In the first case the diagonal uses
a_ell followed by a_(k-d). Thus each permanent at rho=0 is nonzero.

A column b of B becomes an all-zero column of T[A,B-rho] whenever

    rho in [b+1,b+delta_i] modulo k.                          (2)

Indeed the coefficient indices b-rho-u, for 0<=u<=e_i, all lie in
the final zero interval [k-d_i+1,k-1] modulo k.

The intervals (2) cover every rho=1,...,k-1 across the sources:
source 1, whose ell_1=0, covers [1,delta_1] and [k-d_1+2,k-1].
For each 1<i<=J they cover [ell_i+1,ell_i+delta_i], as follows.
If a_ell!=0, ell is a chosen column. If a_ell=0 and z<e, the chosen
columns include the entire interval [l+1,j], hence include ell. If z>=e,
the columns l+1,...,l+e give the contiguous covered interval
[l+2,l+d-1]. It reaches j, since j=l+z+1<=l+d-1. Adding column j extends
coverage through j+delta_i. This contains [ell+1,ell+delta_i].
There is no omitted gap since l<ell<j. Together the first J sources
therefore cover [1,sum_(i<=J)delta_i], meeting the tail interval from
source 1. Additional sources cannot spoil the zero-column property.

For each nonzero rho, at least one minor in (1) has a zero column.
Thus only the rho=0 term survives. It is a product of nonzero field
values, and hence is nonzero. This proves the independent-subspace
coefficient assertion with all its hypotheses.

## 4. Decompose an arbitrary affine subspace into blocks

Let A be an affine subspace of K^r. Gaussian elimination on an F_2-basis,
respecting the order of blocks, gives linear subspaces V_i of K and
linear maps sigma_ij:V_j->K (j<i), with a fixed gamma in A, such that

    A={(gamma_i+z_i+sum_(j<i)sigma_ij(z_j))_(i=1)^r : z_i in V_i},
    dim(A)=sum_i dim(V_i).                                  (3)

For completeness, filter the direction space of A by successively
requiring its first 0,1,...,r blocks to vanish. Split each successive
quotient by choosing an F_2-basis. Its leading block images form V_i;
later coordinates of those basis vectors give the sigma_ji. Their direct
sum is the original direction space, proving (3) and the dimension sum.
This is ordinary linear Gaussian elimination, not a circuit direct sum.

Let Q_i map K onto V_i as above (Q_i=0 if V_i={0}). The compositions
sigma_ij composed with Q_j are linear maps on K, so have linearized
polynomial representations Q_ij. Set

    R_i(X)=gamma_i+Q_i(X_i)+sum_(j<i)Q_ij(X_j).

The map R:K^r->K^r has image A. Write d_i=dim(V_i) and D=dim(A).
We will exhibit a positive-degree coefficient of the reduced polynomial
f(R(X)).

## 5. Select t high-dimension blocks

Assume D>kr/t+kt+r(r+1). Put S={i:d_i>r+1}. Because each d_i<=k,
D>kt+r(r+1) implies |S|>=t. Also

    sum_(i in S)d_i >= D-r(r+1) > kr/t+kt.

Choose t indices I0 in S with largest d_i (break ties by index). Their
sum is at least t/r times the sum over S, since |S|<=r. Therefore

    sum_(i in I0)(d_i-(r+1))
        > k + k*t^2/r - t*(r+1) >= k,                      (4)

where the last inequality is k*t>=r(r+1).
For each selected i, e_i=r+1-i is at least 1, and e_i+1=w_i<=r+1<d_i.
Thus e_i<d_i-1 and d_i>1. The independent-subspace assertion applies,
since sum_(i in I0)(d_i-w_i)>k>=k-max_(i in I0)d_i+1.
It supplies a monomial

    M=product_(i in I0) X_i^(b_i),    wt(b_i)=w_i,

with nonzero coefficient in
Tr(mu*product_(i in I0)Q_i(X_i)^(2^(w_i)-1)).

## 6. Triangular mixing does not alter that coefficient

The same monomial has exactly the same coefficient when Q_i(X_i) is
replaced by R_i(X) in the I0 summand. Expand each exponent as its w_i
Frobenius factors, and expand the trace. Among selected indices, consider
the largest i first. No lower-index factor mentions X_i. Achieving binary
degree w_i in X_i requires all w_i factors from R_i to contribute an
X_i term, hence a term from Q_i(X_i), not a constant or an earlier
variable. Carries cannot increase degree. Having used those factors,
repeat for the next-largest selected index. Descending induction forces
the diagonal contribution in every block, proving equality of coefficients
for each trace summand and hence their sum. This does not assume the
block coordinates of A are independent.

## 7. No other symmetric summand cancels it

The remaining point is isolation of M from summands indexed by J!=I0
with |J|=t. In a trace expansion of the J summand, index j supplies
w_j=r+2-j linearized factors, each mentioning only variables of index
<=j or a constant. Suppose a product choice could reduce to M. Let
n_ji count factors from j choosing variable i in I0. Necessarily

    sum_(i in I0)n_ji <= w_j,
    sum_(j in J)n_ji >= w_i for every i in I0,
    n_ji>0 implies i<=j.                                   (5)

The second inequality is binary-degree subadditivity, also valid after
reduction. Since w_i strictly decreases with i, i<=j implies w_i>=w_j.
A short weighted-capacity argument now replaces any implicit uniqueness
assumption:

    t <= sum_(i in I0,j in J) n_ji/w_i
       <= sum_(i in I0,j in J) n_ji/w_j <= t.                 (6)

The first inequality uses the demands in (5); the last uses the supplies.
If a positive n_ji has i<j, the middle inequality is strict, impossible.
Thus every used factor has i=j. Each target index has positive demand,
so I0 is contained in J. Their equal cardinalities imply J=I0.
Consequently no other symmetric summand contributes to M.

The coefficient of M in f(R(X)) is therefore the nonzero coefficient
from section 3. M has positive degree, so the reduced polynomial is
nonconstant. Since R(K^r)=A, f is nonconstant on A. This proves the lemma.

## 8. Exact application and source scope

For proof.md, k is a prime supplied by the evaluator's length rule,
K=F_2[Z]/(Q) for an accepted prefix, and pi is its constant-coordinate
functional. The guard supplies r>t>=1 and k*t>=r(r+1). The zero-padded
payload embedding is injective, so every payload affine space of dimension
d_N=floor(kr/t)+kt+r(r+1)+1 maps to an affine space satisfying the lemma.
No field-encoding loss besides the counted prefix k is suppressed.

The source route is BK Lemmas 3.7,3.10,3.13; Theorem 6.1; Lemma 7.2;
Theorem 7.3 and Claim 7.4. We specialized the zero-pattern argument to
prime k and used (6) to give a direct proof of the final combinatorial
uniqueness. Only elementary finite-field algebra, root counting, linear
algebra and binary polynomial expansion remain outside this file. The
finite-field existence and polynomial-time evaluation arguments are
supplied in proof.md. No efficient search for a hard truth table or for
an irreducible prefix is assumed.


## 9. A sharper threshold when t squared is at least r

For the same polynomial, assume in addition t^2>=r. Then f is
nonconstant on every affine space of dimension strictly greater than

    k*t + r(r+3)/2.                                       (9.1)

Here all r blocks remain counted, and no family or field construction
changes. We replace only the block-selection step in Section 5.

Use the block decomposition in Section 4. Put w_i=r+2-i as before,
S=sum_i w_i=r(r+3)/2, delta_i=d_i-w_i, and D=sum_i d_i.
If D>kt+S, then sum_i delta_i=D-S>kt. Since each delta_i<k,
there must be at least t positive delta_i: otherwise their entire
positive sum is at most (t-1)k, contradicting the displayed inequality.
Choose I0 to be the t indices of largest positive delta_i, breaking
ties by index. If there are h positive indices, then h<=r and averaging
gives

    sum_(i in I0) delta_i
        >= (t/h)*sum_(delta_i>0) delta_i
        > (t/r)*k*t >= k.                                 (9.2)

For every selected i, d_i>w_i=e_i+1, where e_i=r+1-i>=1.
Thus d_i>1 and 1<=e_i<d_i-1, exactly the individual hypotheses
of Section 3. Its sum hypothesis also holds because the sum in (9.2)
is greater than k, hence greater than k-max_(i in I0)d_i+1.
Apply the independent-subspace coefficient assertion to these t blocks.
Sections 6 and 7 apply to any such selected set and preserve and isolate
its nonzero positive-degree monomial. Consequently f is nonconstant
on the original affine space, proving (9.1).

This argument uses the excess d_i-w_i directly. It avoids discarding
all blocks below a common dimension threshold and paying for that
discard in addition to the later averaging loss. It is a strengthening
of the supplied construction analysis, not a change to the polynomial
or to the input encoding.


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



# A2. a cubic trace disperser

Let p be a sufficiently large odd integer, K=F_(2^p), mu!=0, and
F(x)=Tr(mu*x^7). Claim: F is nonconstant on every affine subspace of
dimension at least d=ceil(9p/20)+10. It suffices to consider a space of
exactly dimension d. Put l=p-d.

The image-polynomial argument in submission/disperser.md Sections 1--2
supplies a monic linearized polynomial Q(X)=sum_(j=0)^l a_j X^(2^j)
whose image is the direction space, with a_0,a_l nonzero. Extend the
coefficient indices modulo p. The final interval [l+1,p-1] consists of
d-1 zero coefficients, and every OTHER cyclic gap has at most d-2 zeros. Here is the extension
of the earlier prime-degree stabilizer proof to the odd-degree case.
A second length-(d-1) gap would imply Q(X)^(2^h)=cQ(X) as functions,
for some 0<h<p and c!=0, by the identical root-count argument in Section 2.
All ratios z/z0 for nonzero z0 in the image space are then fixed by
Frobenius^h and Frobenius^p. Bezout on powers of this invertible
field automorphism makes them fixed by Frobenius^g, g=gcd(h,p).
The polynomial X^(2^g)-X has at most 2^g roots, so the direction space
has dimension at most g. Since p is odd and g is a proper divisor,
g<=p/3<d, a contradiction. No prime distribution theorem is used.

The affine image is Q(K)+c. In the reduced polynomial Tr(mu*(Q(X)+c)^7),
terms involving c have binary degree at most two. Thus a nonzero reduced
coefficient of binary weight three proves nonconstancy on that affine
space. For B a three-element set of residues, set beta=sum_(b in B)2^b.
The coefficient of X^beta is

 sum_(r=0)^(p-1) mu^(2^r) * perm(T[{0,1,2},B-r])^(2^r),
 T[u,v]=a_(v-u)^(2^u).

This follows by expanding the three linearized factors Q,Q^2,Q^4;
a full-weight term has distinct Frobenius exponents, hence is indexed
by a bijection of columns and contributes the permanent. Carries lower
binary weight, and other Frobenius trace terms shift all columns. The
argument is exactly the general coefficient identity already proved in
submission/disperser.md Section 3. In characteristic two permanent is
determinant. We choose B making the r=0 determinant nonzero and all
other determinants zero.

## A middle nonzero coefficient

If there is a nonzero a_j with

 p-2d+4 <=j<=d-4,

choose B={0,j+1,l+2}. The unshifted indicator matrix is triangular,
with diagonal a_0,a_j^2,a_l^4 nonzero. Its zero entries use the final
zero interval and j>=1,j<=l-1, true for these parameters.

For 1<=r<=d-3 the column at 0 is all zero: its coefficient indices
are p-r,p-r-1,p-r-2, inside [l+1,p-1]. For p-d+3<=r<=p-1,
the column at l+2 is all zero for the same reason. The remaining
r satisfy d-2<=r<=l+2. In that interval the indices of the middle
column are j+1-r,j-r,j-1-r. They are negative, at most -1 since
j<=d-4, and at least -d+1 since j>=p-2d+4. Modulo p they therefore
belong to [l+1,p-1]. Hence its column is zero. Only r=0 contributes.

## Two separated clusters

Otherwise every nonzero coefficient lies in one of

 [0,A], [B0,l], where A=p-2d+3, B0=d-3.

Both clusters are nonempty, since a_0,a_l!=0. Let alpha be the last
nonzero index in the first cluster, and j2 the first in the second;
put beta=l-j2. Then 0<=alpha,beta<=p-2d+3=A. All nonzero column
positions for the three-row matrix lie in the two intervals

 I=[0,alpha+2], J=[j2,l+2],

because a column v is nonzero only if one of a_v,a_(v-1),a_(v-2)
is nonzero. The intervals are disjoint and do not wrap. Both cyclic
gaps between them exceed A+2 for sufficiently large p: the internal
gap is at least j2-alpha-2 >=3d-p-8, while the wrap gap is
p-(l+2)=d-2. Both exceed p-2d+5 because d>=.45p+10. Therefore two
column positions separated by at most A+2 can both be nonzero only
when they lie in the SAME interval I or J. This remains true after
any cyclic translation.

If alpha>beta, choose B={0,l+1,alpha+2}. The first and last positions
are the endpoints of I, separated by its full width alpha+2. A cyclic
translation can place both in I only with translation zero, and cannot
place both in the shorter interval J. Thus every nonzero shift has
one of those columns zero. The unshifted matrix is triangular, with
nonzero diagonal a_0,a_l^2,a_alpha^4.

If beta>alpha, choose B={j2,alpha+1,l+2}. The first and last positions
are endpoints of J, and the same width argument applies. The unshifted
matrix is triangular, with diagonal a_j2,a_alpha^2,a_l^4. The indicated
zero entries use the gap between clusters, whose length is greater than
two. Thus the unshifted determinant is nonzero.

Finally let alpha=beta and choose B={0,j2+1,alpha+2}. The two endpoint
positions 0,alpha+2 can both be nonzero only for zero translation,
or for translation +j2 which takes I to J. Consequently only r=0
and r=p-j2 can possibly have nonzero determinant. The zero-shift
matrix is upper triangular with diagonal a_0,a_j2^2,a_alpha^4.
For r=p-j2 its columns at 0 and alpha+2 have shapes

 (a_j2,0,0) and (0,0,a_l^4),

respectively, so its determinant is a_j2*a_(2j2)^2*a_l^4, up to a
sign irrelevant in characteristic two. We show a_(2j2)=0, with the
index interpreted modulo p. If 2j2<p, then 2j2>=2d-6>l and hence
is in the final zero interval. If 2j2>=p, use alpha=l-j2. The zero
gap between alpha and j2 has length j2-alpha-1=2j2-l-1>=p-l-1=d-1.
This would be a second cyclic gap of at least d-1 zeros, forbidden by
the odd-degree image-polynomial property. Hence this second case
cannot occur. The exceptional shifted determinant is zero, completing
the proof of the cubic disperser claim.

## Instantiating a half-rate base extractor

For slice length w, set p=2w+1. Encode a pair of w-bit slices as
the first 2w field coordinates, padding one coordinate by zero. If one
slice is uniform, the affine source has entropy at least w=(p-1)/2,
which is >=.49p for sufficiently large w. The cubic-bias component
then implies every nonzero trace projection of x^7 has bias less than
2^-floor(sqrt(p)/1000): a hypothetical larger bias would give a constant
space of dimension >=.49p-246(floor(sqrt(p)/1000)+1)^2, greater than
ceil(.45p)+10. Taking floor(p^(1/8)) independent trace projections gives
TV error <=(1/2)*2^(p^(1/8)/2)*2^-floor(sqrt(p)/1000).

The polynomial basis coefficients 1,T,... yield independent trace
functionals by nondegeneracy of the trace pairing. The needed modulus
is an actual input-prefix field datum in the eventual family. The evaluator computes p=2w+1 directly. It need not find
an irreducible modulus: the hard restriction fixes a good one, whose
existence follows from elementary finite-field theory. All prefix bits
must count in total input length. The wrapper in proof.md supplies the exact data layout and charges its size.


# A3. biased cubics have a large constant restriction

**Claim.** If a cubic f on a k-dimensional affine space has bias at least
2^-L, L>=1, then it is constant on an affine subspace of dimension at least
k-246(L+1)^2. Coordinate changes identify the source with F2^k.

For a quadratic q write R(q) for half the rank of its alternating polar
form. Symplectic elimination gives q=sum_{j=1}^R a_j b_j + ell in independent
linear coordinates, plus a constant. Its bias is either zero or 2^-R:
when the linear part has a component outside the paired coordinates it
cancels; otherwise summing each pair gives magnitude 1/2. In particular,
nonzero bias at least 2^-u implies R<=u. Fixing one coordinate in every
pair, and then the remaining linear form, makes q constant at codimension
at most R+1. R is subadditive and unchanged by linear terms.

## Dense sets and bounded-rank directions

A set A of density mu in F2^k contains, in its 2j-fold sum, a linear
subspace V of codimension <=2 log2(1/mu)+1, with
j<=ceil(2 log2(1/mu))+2. Here a sum uses exactly 2j elements, repetitions
allowed. To see this, iteratively restrict A to an affine hyperplane
whenever a nonzero normalized Fourier coefficient has magnitude at least
3mu/4. A suitable hyperplane has density at least 7mu/4 (in particular
at least 3mu/2). Continue until density exceeds 1/2 or all nontrivial
coefficients are smaller than 3mu/4. The number of restrictions is at
most 2 log2(1/mu)+1. In the first case A+A fills the direction space,
and repetitions pad to 2j. In the second, translate the current affine
space to a vector space, and let B be its intersection with A, of density
nu>=mu. The normalized 2j-fold convolution at every point is at least
nu^(2j) - (3nu/4)^(2j-2)*nu, by Parseval. This is positive for the stated
j since (3/4)^(2j-2)<mu<=nu. Translating back adds an even number of the
same offset, hence gives the direction space V inside 2j A. This proves
the assertion without a separate additive-combinatorics premise.

Put delta=bias(f)>=2^-L. The identity
E_y E_x (-1)^(f(x+y)+f(x))=delta^2 implies that at least delta^2/2 of
directions have derivative bias at least delta^2/2. Call these directions A;
their quadratic derivative rank is at most 2L+1. The map
T:y -> quadratic homogeneous part of f(x+y)+f(x) is F2-linear; its polar
forms add exactly. Thus every direction in 2j A has derivative rank at
most 2j(2L+1). With mu>=2^(-2L-1), the preceding paragraph gives V of
codimension <=4L+3, and every derivative in V has rank at most
(8L+12)(2L+1)<=60(L+1)^2. Denote this last upper bound by R.

## A space of bounded-rank alternating forms

For any linear space M of quadratic homogeneous parts, all of rank <=R,
there is a subspace W of codimension <=4R on which all their polar forms
vanish. Here and below rank means number of symplectic pairs.

Choose g in M of maximal rank r<=R, and let its 2r independent symplectic
coordinates be a_1,...,a_(2r). Set these coordinates to zero, obtaining W1.
For any h in M let its restriction to W1 have s symplectic pairs, in
coordinates b_1,...,b_(2s). The a's and b's are independent. Completing
squares in the alternating form, equivalently replacing each b by b plus
a linear combination of a's, splits h into s independent pairs plus a
form h0 on the remaining coordinates, including the a's. Hence
R(h0)=R(h)-s<=r-s. The same split applies to g+h, since g uses only a's.
Consequently R(g+h)=s+R(g+h0)>=s+r-R(h0)>=2s. Maximality gives s<=r/2.
Apply the argument to the restricted space on W1, whose maximal rank has
halved. The total imposed codimension is bounded by
2r+2(r/2)+2(r/4)+...<=4r<=4R. Linear terms play no role in this argument.

Apply this lemma to M=T(V), obtaining W. Put U=V intersect W, so its
codimension is at most 4L+3+240(L+1)^2 <=244(L+1)^2. For every y in U,
the derivative of f along y, restricted to any affine coset of U, has
no quadratic part: translating the source changes a quadratic derivative
only by linear terms. It follows, by differentiating a multilinear
polynomial in a basis of U, that f restricted to each coset of U has
degree at most two. The signed biases on these cosets average to the
signed bias of f. At least one coset therefore has absolute bias >=delta.
Its quadratic restriction has rank <=L and is constant after at most
L+1 further affine constraints. Total codimension is at most
244(L+1)^2+L+1 <=246(L+1)^2, proving the claim.

## Consequence using A2

Suppose a cubic F on n bits is nonconstant on every affine space of
dimension greater than .45n+11. On any source of dimension k>=n/2,
if its bias were >=2^-L, the claim would produce a constant space of
dimension at least n/2-246(L+1)^2. For example, taking
L=floor(sqrt(n)/1000), at sufficiently large n this exceeds .45n+11.
Thus bias is smaller than 2^-floor(sqrt(n)/1000).

If F_j(x)=Tr(mu_j*x^7), with linearly independent nonzero trace
coefficients over F_(2^n) and n odd, then every nonzero binary linear
combination has this property by A2. Taking ell=floor(n^(1/8)) coordinates,
Parseval/Cauchy--Schwarz gives total variation from uniform at most
(1/2)*2^(ell/2)*2^-floor(sqrt(n)/1000), hence stretched-exponential error.
The exact odd-degree-padding base extractor is specified at the end of A2.


# A4. affine somewhere-random condensation

## A near-half-entropy linear strong seeded extractor

Let an affine input have ambient length n<=r^1.1 and entropy >=3r/4.
Set D=100000, h=D*ceil(n^(1/D)), m0 the largest power of two <=r^.999,
and epsilon0=2^-floor(r^.00001). Use the low-degree seeded extractor component
with output m0, error epsilon0 and minimal power of two q meeting its bound.
Here rho=epsilon0/(8m0) is reciprocal to an integer. Its field degree
s=O(log r+r^.00001), while h=O(r^.000011). Thus its entropy requirement
16m0*(16D+16)*h*s+O(log(1/rho)) is o(r), in particular <r/4 eventually.
Its seed length is O(r^.00001), with constants depending on D only.

Concatenate independent calls on the SAME input, using fresh independent
seeds each time, until M=floor(r/2) output bits have been produced (truncate
the last call). There are O(r^.001) calls. Given all prior seeds and outputs,
the remaining input is affine of dimension >=3r/4-M>=r/4. Thus each call
is strong with error epsilon0 on every conditional input distribution.
A hybrid replacing successive outputs by independent uniforms gives total
error <=O(r^.001)*epsilon0. The output for each fixed concatenated seed is
linear. The total seed length is O(r^.00101), hence <=r^.002 eventually.
The product of TV error and two bounds the fraction of bad seeds for ANY
fixed affine source: a fixed-seed affine output either is exactly uniform
or has TV distance at least 1/2 from uniform.

One may apply this extractor at entropy >=r-r^.3 as well: it exceeds 3r/4
for large r. Good field and matrix data are mathematical construction data;
actual evaluator prefixes must encode and pay for them. The many repeated
calls share these data, using independent seeds only.

## One condensation round

Input: an affine t-by-r binary matrix with at least one uniform row and
t<=r^.1. Let w=floor(r^.1). Slice the first w bits from each row, and pair
rows, duplicating the last one when t is odd. For each pair, concatenate
the two w-bit slices and pad by zeros to p bits, where p=2w+1 is odd. Apply the cubic trace extractor from the accompanying supplied
coefficient and bias arguments. One pair includes a uniform w-bit slice,
so its entropy is at least w; w/p=(p-1)/(2p)>=.49 eventually. The output
of that pair is therefore close to uniform. Keep enough output
bits to serve as a seed for the near-half extractor above, padding the
seeded construction's seed space with ignored bits if necessary. This
fits because the base extractor supplies p^.125=r^(.0125+o(1)) bits,
whereas the seed length is <=r^.002. Base error is stretched exponential
in w, for example <=2^-r^.04 eventually.

Apply the near-half extractor to the ENTIRE t-by-r input, separately with
each of the ceil(t/2) computed seeds. Output ceil(t/2) rows of length
floor(r/2). Its ambient input length rt<=r^1.1, as required.

Proof: let Z be all the sliced bits, at most tw<=r^.2. Split the original
affine input as A+B, where A is a deterministic affine function of Z
and B is a uniform linear source independent of Z, of entropy at least
r-tw>=r-r^.2. This is obtained by choosing a vector-space complement to
the kernel of the slicing map. For each fixed Z all computed seeds are
fixed; every output row is linear in B plus a fixed offset. Hence the
full output conditioned on Z is affine. B is fixed as a distribution,
independent of which value of Z is chosen.

For a uniformly random seed the near-half extractor's output on B is
exactly uniform except on a fraction at most twice its strong error.
The seed from the good pair has TV distance at most beta from uniform,
so the probability over Z of choosing a bad seed is at most beta plus
twice the strong error. Conditional on any Z outside this bad event,
one output row is uniform. Consequently the output is within that error
of a convex combination of affine somewhere-random matrices. This
argument retains joint affine structure after conditioning; no independence
between the two paired slices, or between different seeds, is presumed.

## Iteration

Iterate ceil(log2 t) rounds, always pairing rows and halving row length.
Write r_j=floor(r/2^j), t_j=ceil(t/2^j). For nonterminal rounds and large
r, t_j<=2 r_j^.1; replace the one-round t<=r^.1 condition by t<=2r^.1,
which changes only a factor two in the slice-entropy budget and in the
ambient-length bound. To avoid that ambient-length factor, pad the basic
seeded extractor's allowed input n to ceil(2r^1.1), which leaves h's
asymptotic exponent and all entropy inequalities unchanged.

There are O(log r) rounds. The final row length is at least r/(4t)>=r^.9/4.
At every stage the row length is at least this value. The repeated
seeded extractor error, of order polynomial(r_j)*2^(-r_j^.00001),
dominates the base-extractor error. After r_j>=r^.9/4 and summing
O(log r) rounds, it is bounded by 2^(-r^.000001). Total variation
errors add under deterministic transformations, and convex mixtures can
be analyzed componentwise. The final single row is therefore within
2^-r^.000001 of uniform, for sufficiently large r. Its length is at least
r^.9/4; one may uniformly truncate to floor(r^.89). Runtime is polynomial
in rt and in the lengths of the supplied field/matrix data. This gives
a somewhere-random extractor on at most r^.1 affine rows, with output
r^.89 and stretched-exponential error.



# A5. block-wise affine extraction

Let N be an ambient length bound, and K an entropy threshold satisfying
log K >= (log N)^.7, for sufficiently large N. All logs are binary.
Blocks have length at most N and are zero-padded to N for seeded calls.
We require an extractor on b-block affine sources: conditional on every
possible fixing of earlier blocks, the next block has entropy at least K.

Choose D=ceil(100 log N/log K), h=D*ceil(N^(1/D)), r the largest power
of two <=sqrt(K), and R=16D+16. Then h=K^(.01+o(1)), D=K^o(1).
The low-degree extractor with output r and constant error 1/8 has
s=O(log K+log D), seed length d=(D+1)s<=C log N for an absolute constant
C. Its advice/entropy requirement is K^(.51+o(1))<K. Denote this
linear strong seeded extractor by E0. Its 2^d rows on the first block
form an affine somewhere-random r-bit-row matrix: at least one row is
exactly uniform, because strong error <1/2 and nonuniform affine images
have distance at least 1/2 from uniform.

Use a second low-degree extractor E1 with the same D,h,r but error
2^-floor(K^.001). Its field degree is O(K^.001+log K), its entropy
requirement is K^(.511+o(1))<K, and its seed length is K^(.001+o(1)),
hence <=K^.002 eventually. Its fraction of bad seeds on every fixed
affine source of entropy at least K is at most twice its error.

Let T be the largest power of two <=r^.1. Group the current matrix's
rows into consecutive chunks of at most T rows, padding the final chunk
by repeating its last row. The affine somewhere-random extractor component
maps each chunk to floor(r^.89) bits, with error <=2^-r^.000001. This
supplies more than K^.002 seed bits, since r^.89=K^(.445+o(1)).
Apply E1 to the next input block with the prefix of each chunk output
as its seed, obtaining a new matrix of ceil(v/T) rows of length r.
Repeat with successive input blocks until there is one row. The number
of needed blocks is

  b=1+ceil(d/log2 T)=O(log N/log K).

Here row count starts at 2^d, T is a power of two, and repeated ceiling
agrees with ceil(2^d/T^j), so that this exact b suffices. Extra blocks,
if supplied, can be ignored.

## Why correlated affine blocks suffice

If (A,B) is affine, then B=L(A)+B0 for an affine map L and a linear
source B0 independent of A, of dimension H(B|A). This follows by
choosing a complement to the kernel of the projection onto A. For
fixed A=a, any fixed-seed linear extraction of B is a translate of the
same extraction of B0. Thus every seed good for B0 gives a uniform
output for every possible fixing a.

Assume a current row matrix W is an affine somewhere-random source,
and the next block B has H(B|W)>=K. Let B=L(W)+B0 as above. One chunk
of W contains a uniform row, so its chunk extractor produces a seed
within beta=2^-r^.000001 of uniform. A uniformly random seed fails on B0
with probability at most 2epsilon1. Consequently, for a set G of W
values of probability at least 1-beta-2epsilon1, that chunk's seed is
good for B0. Conditional on any w in G, all new seeds are fixed,
all new rows are linear images of the same conditional affine B,
and one row is exactly uniform. Hence the conditional output matrix
is affine somewhere-random. There is no assumption that different
rows or chunks are independent.

For the full b-block source, the induction conditions on actual earlier
input blocks as needed. Conditional on any successful history, W is
an affine linear image of the most recently unconditioned input block,
and H(B|W,earlier history)>=K because conditioning on all preceding
input blocks would still leave entropy K. The joint pair (W,B) is
affine. The previous paragraph therefore produces a new mixture of
affine somewhere-random matrices with failure probability at most
beta+2epsilon1. Iterating and summing errors gives final distance
<=b*(beta+2epsilon1). This is at most 2^-K^.0000001 eventually, because
r^.000001=K^(.0000005+o(1)), epsilon1=2^-K^.001, and b=K^o(1).

The single final row has length r>=sqrt(K)/2 and may be truncated.
For the global choices below, put S=ceil((L^3*z)^(1/4)), z=ceil(log2 L).
Then log2 k=1000000000000*S, log2 K=500000000000*S, and
log2 pmax=1000*S. The domain log K>=(log N)^.7 holds eventually.
Both sqrt(K) and K^.0000001 exceed pmax: their binary logarithms
are respectively 250000000000*S and 50000*S, larger than 1000*S.
Thus BE supplies pmax bits with error <=2^-pmax.


# A6. subpolynomial-entropy affine dispersal

## Parameters and algorithm

Let L=ceil(log2 N), z=ceil(log2 L), S=ceil((L^3*z)^(1/4)), and
set the exact integer scales

 k=2^(1000000000000*S),  t=2^ceil((L*z^3)^(1/4)),
 ell=2^ceil(sqrt(L*z)),  pmax=2^(1000*S),
 q=floor(sqrt(floor(L/z))).

Fractional powers inside ceilings/floors are computed by integer root
comparisons. Use a complete t-ary tree with ambient padded length
n=t^ceil(L/log2 t), between N and Nt. Truncate at the first depth h
whose node length is below floor(k^.4), and put hA=floor(h/q).
The following asymptotics hold as L tends to infinity:

 S=(1+o(1))*L^.75*(log2 L)^.25,
 h=(1+o(1))*L^.75*(log2 L)^(-.75),
 hA=(1+o(1))*L^.25*(log2 L)^(-.25).

The final active node has length at most k^.4*t^q. Its binary logarithm
divided by log2 k tends to .4+10^(-12)<.5, so that node is shorter
than sqrt(k) eventually. Only nodes having a next active level suggest
an output.

Use BE at entropy K=2^floor(log2 k/2), ambient n, output pmax and error
<=2^-pmax. Its block count b=1+ceil(d0/log2 T), with d0<=10000*log2 n
and log2 T=25000000000*S. Therefore

 limsup b/[L^.25*(log2 L)^(-.25)] <= 4*10^(-7).

At each node, try all increasing b-tuples of children, padding blocks
to length n. There are at most t^b candidates. The displayed bound
implies limsup log2(t^b)/log2 ell <=4*10^(-7)<1; multiplying by h
still leaves h*t^b<=ell eventually. Every challenge has at most ell
rows of pmax bits; zero-pad missing rows.

Let E be the constant-error linear strong extractor at entropy K on n
bits, output ell*pmax, obtained by truncating E0. This fits because
log2(ell*pmax)=1000*S+o(S), while log2(K^.5)=250000000000*S.
Its seed length is <=C log n for an absolute constant C; use V=n^C
as a uniform upper bound on the number of responses.

For a node s and descendant s', define CR(ch,p;x_s,s') to fail if some
seed y makes the first p bits of EVERY padded challenge row ch(x_s')
coincide with the corresponding row of E(x_s,y). Here x_s is padded to
length n with trailing zeros; response output is interpreted as ell rows
of pmax bits. Otherwise CR passes. The test uses at most V responses.
Larger confidence p gives a smaller fail event.

Compute challenges bottom-up. At each internal node s, add as rows all
BE outputs on its b-tuples of children. If its children themselves have
challenges, select the first child s' passing CR at confidence pmax,
and include all that child's challenge rows as well; if none pass choose
the first child arbitrarily. At the bottom internal level there are no
inherited rows. Define challenges at leaves as all-zero matrices, so
leaf tests are total even though the proof will stop above them.

The selected first-passing children define a root-to-leaf path. For an
active node s on level u*q with a next active node s' on that path, put

  p_u^(a)=ell^(6(hA-u)+3a), a=1,2,
  o_s^(a)=CR(ch_s',p_u^(a);x_s,s').

These confidence values are <=ell^(12hA+6), whose binary logarithm
divided by log2 pmax tends to 12/1000<1. Thus they are all <=pmax
eventually, and even the smallest exceeds (C+10)log n. Choose the first active node on the path whose o^(1) fails,
and output its o^(2). If none exists output zero. Pass and fail are
encoded as distinct Boolean values.

## Three elementary facts about affine sources

Every affine conditioning on linear observations decreases dimension by
at most the rank of those observations. Conditional distributions remain
affine and their dimensions do not depend on the observed values.
For affine blocks, conditional entropy equals the dimension of a kernel.
These facts follow from elementary linear algebra.

(1) Empty-block failure. If an affine source Y has entropy >=K and its
block Y_s' is fixed, some seed of E maps Y to a uniform response matrix:
otherwise every fixed-seed affine image has distance >=1/2, contradicting
strong error <=1/8. Equating its ell*p sliced bits with the fixed
challenge defines an affine subsource of deficiency ell*p on which CR
fails. More generally if H(Y|Y_s')>=K, conditioning on each block value
shows CR fails with probability at least 2^(-ell*p).

(2) Resilient passing. Say a challenge is eta-somewhere-random with
resiliency R relative to Y when, after every possible linear observation
of rank <=R on Y, at least one of its rows has distance <=eta from
uniform. The good row may depend on the observed value. If R>=ell*p,
then CR fails with probability at most V*(2^-p+eta): fix each response's
sliced value (a linear observation of rank <=ell*p); the challenge has
one nearly uniform p-bit row, whose chance to equal that fixed response
row is at most 2^-p+eta; then union over seeds. Resiliency relative to Y
implies resiliency relative to any projection of Y.

(3) A hidden block source. If b children in increasing order have
conditional entropy at least K given preceding children, their tuple
is a b-block affine source. A node's added BE rows then have a nearly
uniform row, even after an additional rank R observation, provided those
conditional entropies were >=K+R. Ignoring other preceding children can
only increase conditional entropy.

## Structural subsource and two good active nodes

Start with any affine source X of entropy >=k on the padded n coordinates.
We prove that it has a subsource X0 and nodes e subset g' subset g, with
g,g' consecutive active nodes, such that the coordinates left of e are
fixed, e is a hidden b-block source with threshold k^.75, and
H(X0_g|X0_g')>=k^.75.

First maintain a node s at an active level, with all earlier coordinates
fixed and entropy at least k/4^j in s. Partition into T=t^q next-active
children. For threshold a=k/4^j, let i be the first child whose entropy
conditional on earlier children is at least a/(4T). The total entropy
in earlier children is <=a/4. Fix them. If the chosen child's remaining
entropy is >=a/4, descend to it and increment j. Otherwise its entropy
is <a/4 while the parent retains >=3a/4; so conditional parent entropy
outside this child is >=a/2, and the pair is strongly structured, with
both child and conditional parent entropy >=a/(4T). The descent cannot
reach the last active level: k/4^hA >k^.9, whereas that level's node
length is <sqrt(k). Thus it stops at a pair g,g', with both entropies
>=k/(4^(hA+1)*t^q)>=k^.8, and all coordinates left of g' fixed.

Now maintain a node s below g', with threshold a initially k^.8 and
with coordinates left of s fixed. In its t-child partition, if at least
b conditional child entropies are >=a/(4t), it is a hidden block source
and we stop. Otherwise the <b large conditional entropies carry at least
3a/4 total entropy, and some child has conditional entropy >=a/(4b).
Fix all preceding children and descend to that child, replacing a by
a/(4b). These fixings lie inside g', so they do not change the kernel
dimension H(X_g|X_g'). There are at most h descents, and
k^.8/(4b)^h >= k^.77 eventually. Indeed,
[h*log2(4b)]/[log2 k] <=1/(4000000000000)+o(1), because
log2 b<=(1/4)log2 L-(1/4)log2 log2 L+O(1) (or is smaller, which only helps).
Thus h*log2(4b)<=.02 log2 k eventually, leaving at least k^.78.
Leaves have length <k^.4, so this process stops at an internal node e.
Its hidden-source threshold is at least k^.8/((4b)^h*4t)>=k^.75.
The conditional parent entropy is inherited and >=k^.8>=k^.75.

## Making earlier empty children fail

For every proper ancestor s of e and each child s' strictly left of e,
the block s' is fixed. Apply fact (1) to the source projected on s, with
confidence pmax, and lift the resulting linear restriction to X0.
There are at most ht such restrictions, each of rank ell*pmax. The
entire lost dimension is <=ht*ell*pmax=2^((1000+o(1))*S)=o(k^.75).
At every intermediate step the source in s has entropy >=k^.75-o(k^.75)
>=K, so each use of fact (1) is justified. The final subsource Y still
has hidden block entropy >=k^(2/3) at e and conditional entropy
H(Y_g|Y_g')>=k^(2/3). All empty earlier-child tests now fail identically,
and they keep failing on every further subsource.

## Propagating resilient challenges and reaching g'

If s is an ancestor of e at distance a, its challenge on Y_s has error
eta_a and resiliency (h-a+1)*ell*pmax relative to Y, where
eta_0=2^-pmax and eta_a<=(V+1)^a*2^a*2^-pmax.
At e, every rank <=(h+1)ell*pmax observation leaves hidden conditional
entropies >=k^(2/3)-(h+1)ell*pmax>=K; fact (3) proves the base case.
For the induction step, condition on any observation of rank
(h-a+1)ell*pmax. Its correct child challenge retains resiliency ell*pmax,
so fact (2) says the correct-child test fails with probability at most
V*(2^-pmax+eta_(a-1)). All earlier children fail identically. Thus except
on that event the parent inherits the correct child's rows, giving
eta_a <= eta_(a-1)+V*(2^-pmax+eta_(a-1)). This is bounded by the stated
recurrence. Hence all these errors are <=2^(-pmax/2), because
h log(2(V+1))=O(L^1.75)=o(pmax).

Apply fact (2) at every consecutive edge on the path from root to g'.
Each correct child passes with probability at least 1-2^(-pmax/3).
Earlier children always fail. A union bound shows that the selected
path reaches g' except with probability <=2^(-pmax/4).

## Both outputs survive selection of the first active node

Let g have active index u, and abbreviate p1=p_u^(1), p2=p_u^(2).
Analyze r_a=CR(ch_g',p_a;Y_g,g') independently of the selected path.
Conditional entropy H(Y_g|Y_g')>=K implies
Pr[r_a fails]>=2^(-ell*p_a), by fact (1). Resilient passing implies
Pr[r_a fails]<=V*(2^-p_a+2^(-pmax/2)).
Since p2=ell^3*p1, the event r1 fails AND r2 passes has probability at
least 2^(-ell*p1)-V*(2^-p2+2^(-pmax/2)). The event r1 fails AND r2 fails
is exactly r2 fails, by monotonicity, with probability >=2^(-ell*p2).
Both probabilities are therefore >=2^(-ell^(6(hA-u)+7)) for large N.

At each active ancestor with index v<u, its next active node on the path
to g is an ancestor of e. Fact (2) shows its r1 test fails with probability
at most V*(2^-p_v^(1)+2^(-pmax/2)). The nearest such ancestor has
p_(u-1)^(1)=ell^(6(hA-u)+9), so union over them gives failure probability
<=2^(-ell^(6(hA-u)+8)). If there are no ancestors this error is zero.
Combining with the path-reaching error, each of the two joint events
at g retains positive probability after requiring that the selected
path reach g' and that all active ancestors pass at their first tests.
Indeed pmax=2^(1000*S) dominates all exponents
ell^(6hA+8)=2^((6+o(1))*S). On these events the algorithm selects g and
outputs the corresponding distinct r2 value. Thus it is nonconstant on
Y, hence on X. This proves affine dispersal at entropy k=n^o(1). The exact integer
parameters, guards, and charged prefix layout are specified in proof.md
and engine.py.
