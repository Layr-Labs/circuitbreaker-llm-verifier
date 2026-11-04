# A complete certificate-prefix proof with an explicit remainder

## Part A. Family, eligibility and final theorem

## Claim and provenance

For the single total Boolean language defined below, with **N denoting the
entire input length**, we prove

    C_B2(F_N) >= (7/2)N - R(N),       R(N)/N -> 0.

More precisely, for every integer N>=100000,

    C_B2(F_N) >= (7/2)N - 86*N^(4/5).

The joint pair is (7/2,1). This applies existing mathematics; it does not
claim discovery of a new circuit lower bound. The gate argument is the
parameter method of Mishra and Santhanam, ECCC TR26-260, October 19, 2026,
Theorem 4.1, pp. 11-28. The algebraic construction is based on Ben-Sasson
and Kopparty (BK), *Affine Dispersers from Subspace Polynomials*, author
version March 12, 2010, Sections 3, 6 and 7. The necessary arguments are
supplied in our notation in Parts C and B of this proof.

Our input-prefix wrapper removes the deterministic irreducible-polynomial
construction prerequisite of the repository's BK/Shoup seed. It explicitly
pays for all prefix bits. The **whole** function is not an affine disperser:
only a fixed-prefix restriction is asserted to be one. No conjecture,
restricted circuit basis, direct-sum assumption or exponent gain is used.
We assert no exponent above 1, endpoint, or optimality over all families.

### A.1. One evaluator on every bit string

The first bit of each field block is its constant coefficient. On x of length N:

1. If N=0, output 0.
2. Let B be the smallest positive integer with B^5 > 32*N^3. Let k be the
   smallest prime at least B. Use deterministic trial division to find k.
3. If k>=N, output 0. Otherwise put m=N-k, r=ceil(m/k), t=ceil(sqrt(r)).
   If t>=r or k*t<r*(r+1), output 0.
4. Read the first k bits b_0,...,b_(k-1) as

       Q(Z)=Z^k + sum_(j=0)^(k-1) b_j Z^j    in F_2[Z].

   Check Z^(2^k)=Z modulo Q and gcd(Q,Z^2+Z)=1. If either fails, output 0.
5. Use K=F_2[Z]/(Q). Append rk-m zero bits to the remaining m payload bits
   and interpret the r blocks of k bits as a_1,...,a_r in K. This is an
   injective F_2-linear map of the payload.
6. For i=1,...,r set e_i=r+1-i and y_i=a_i^(2^(e_i+1)-1). Output the
   constant coefficient of the elementary symmetric polynomial of order t
   in y_1,...,y_r, in its representative of degree less than k.

For explicit evaluation, initialize E_0=1 and E_1=...=E_t=0. For each i,
update E_j <- E_j+y_i*E_(j-1), in decreasing order j=min(i,t),...,1.
Splitting subsets by whether they contain i proves this division-free
recurrence. Binary powering takes O(e_i+1) multiplications per y_i.

This specifies every choice and every input length. Q is part of x, not
advice or a value supplied by an unspecified construction algorithm.

### A.2. Exact prefix test and existence of hard prefixes

In an algebraic closure of F_2, Z^(2^k)-Z has 2^k distinct roots because
its derivative is 1. The roots are closed under addition, multiplication
and nonzero inversion, hence form a field of size 2^k. Frobenius z->z^2
permutes them and its k-th iterate is the identity. A Frobenius orbit of
length d gives an irreducible degree-d polynomial over F_2: its product
of linear factors is Frobenius-fixed, while a proper factor over F_2 would
require a proper invariant subset of the orbit. Conversely, an irreducible
factor's roots form one orbit. Thus factor degrees divide k.

Since k is prime, the orbits have size 1 or k. Only 0 and 1 are fixed,
and 2^k>2, so an irreducible degree-k polynomial exists for every such k.
No algorithm selecting it is needed in the uniform evaluator.

The first prefix test says Q divides the squarefree polynomial Z^(2^k)-Z;
the second excludes the factors Z and Z+1. Remaining irreducible factors
all have degree k. Since Q has degree k, it must itself be irreducible.
Conversely, every irreducible degree-k Q passes both tests. The test is
therefore exact. Compute its first part by k squarings modulo Q and its
second by polynomial Euclid, both in polynomial bit time.

### A.3. Uniform P and NP with empty witnesses

Bertrand's prime theorem gives k<2B, and so

    2*N^(3/5) < k < 4*N^(3/5)+2.

The prime scan and trial division are polynomial in N. B is computed by
integer binary search and fifth-power comparisons. All subsequent integer,
polynomial and field operations have polynomial bit cost. Field multiplication
can use schoolbook O(k^2) operations; the powers take O(r^2) multiplications
in total and the symmetric recurrence O(rt). A coarse O((N+1)^6) bound
suffices for the whole algorithm, including rejected inputs.

Let A be this fixed evaluator, L={x:A(x)=1}, W(N)=0, and let V(x,w)
accept exactly when w is empty and A(x)=1. If x is in L the empty witness
is accepted; if a permitted witness is accepted then A(x)=1 and x is in L.
This proves both directions with one uniform polynomial-time verifier.
Witnesses are not inputs to the target circuits.

### A.4. Explicit length threshold and guard bounds

Fix any integer N>=100000, and write s=N^(1/5), so s>=10.
By the definition of B and Bertrand's prime theorem,

    2*s^3 < k < 4*s^3+2 <= 5*s^3.

Thus k/N <= 5/s^2 <= 1/20, and m=N-k >= N/2 >0. The block count obeys

    s^2/10 <= m/k <= r < s^2/2+1 <= s^2,
    r+1 < s^2/2+2 <= (3/5)*s^2.

In particular r>=10. Since r is an integer at least 3,
ceil(sqrt(r))<=r-1, so 1<=t<r. Also

    t >= sqrt(r) >= s/sqrt(10),
    k*t > (2/sqrt(10))*s^4 > (3/5)*s^4,
    r*(r+1) <= (9/25)*s^4 < (3/5)*s^4.

The strict comparison 2/sqrt(10)>3/5 follows by squaring positive
quantities: 2/5>9/25. All guards therefore pass at every integer
N>=100000. No statement about a sparse sequence of lengths is used.

### A.5. Hard restriction and its precise loss

For each sufficiently large N fix any irreducible Q of degree k and fix
the first k bits of x to its coefficients. The remaining function G_Q
has m inputs. Put

    d_N = floor(kr/t) + kt + r(r+1) + 1.

Part B proves that, for prime field degree k, r>t>=1 and
k*t>=r(r+1), the function

    (a_i) -> pi(sum_(|I|=t) product_(i in I) a_i^(2^(r+2-i)-1))

is nonconstant on every affine subspace of K^r of dimension greater than
kr/t+kt+r(r+1), for any nonzero F_2-linear pi. Our pi is the constant
coefficient, with pi(1)=1. The payload embedding is injective and linear,
so G_Q is an affine disperser of dimension d_N on its m actual input bits.
Zero padding preserves dimension and does not add circuit inputs.

The estimates in Section A.4 give k=O(N^(3/5)) and d_N=O(N^(4/5)).
The following explicit bounds avoid an unspecified constant in the final loss.
Since kr<m+k=N, and t>=s/sqrt(10),

    kr/t < sqrt(10)*s^4 < 4*s^4.

Moreover t<=sqrt(r)+1<=s+1 and k<=5*s^3, so

    k*t <= (11/2)*s^4.

Finally r*(r+1)<s^4 and 1<=s^4. Therefore

    d_N = floor(kr/t)+k*t+r*(r+1)+1 <= 12*s^4.

These estimates are deliberately conservative. They hold at every
integer N>=100000. No computation is needed to verify them.

Fixing prefix inputs of any circuit for F_N yields a circuit for G_Q with
no additional gates. Part C proves the full-B2 bound

    |C| >= (7/2)m-7d

for a dimension-d affine disperser. It includes arbitrary sharing,
repeated pins, constants, degenerate gates and terminal-output cases.
Apply it to the restricted circuit (or its no-larger normalized form):

    C_B2(F_N) >= (7/2)(N-k)-7d_N
               = (7/2)N-((7/2)k+7d_N).

For N>=100000 put R(N)=(7/2)k+7d_N. As s>=10,

    (7/2)*k <= (35/2)*s^3 <= (7/4)*s^4 < 2*s^4,
    R(N) <= 2*s^4+84*s^4 = 86*N^(4/5).

This proves the displayed explicit theorem. In particular R(N)/N<=86*N^(-1/5)
tends to zero. Define R arbitrarily large on the finitely many smaller lengths;
this does not change the limit. The choice of Q in this lower-bound
restriction is not advice to the evaluator: every circuit must work for
all prefixes, while the uniform algorithm simply reads and tests its input.

For each epsilon>0, take N_epsilon=max(100000,ceil((86/epsilon)^5)).
Then R(N)/N<=epsilon for every integer N>=N_epsilon. This gives exactly
C_B2(F_N)>=(7/2-epsilon)N for all those integers. Both coordinates in
bound.json belong to this one theorem and one language.

### A.6. Evidence and limits

Part B contains the coefficient separation and block-mixing proof.
Part C contains the complete gate argument. This is a proof-only submission:
it declares no computational experiment and no LP certificate. No claimed
review, checker success or prior acceptance is a mathematical premise.

The family and the proofs in Parts B and C derive from the seven-file
candidate retained at `tests/fixtures/candidate-3_5` in the required source
revision. Those fixture bytes also identify the public submission
`f4f2c464-3045-4c32-a6e4-32478c05f343`. The present submission consolidates
its proof and adds the explicit length and error estimates in Part A.
It declares that source as a contribution dependency. It does not reuse
any prior judgment. The scored pair remains (7/2,1); the explicit remainder
is not a larger coefficient or a new accepted-expression class.

Sources, for attribution and version identification:

* Mishra--Santhanam, ECCC TR26-260, 2026-10-19, Theorem 4.1 and Sections 2-4,
  https://eccc.weizmann.ac.il/report/2026/260/download/
  SHA256 38e9ff3cea4ac0399723756f00c2425036c438887ab5ae4abb10d09fa971a3bc.
* Ben-Sasson--Kopparty, author preprint 2010-03-12, Lemmas 3.7, 3.10, 3.13,
  Theorem 6.1, Lemma 7.2, Theorem 7.3 and Claim 7.4,
  https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
  SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf.
* Li--Yang, ECCC TR21-023, 2021-02-20, Sections 1.1, 2.1, 2.4 and 3.3,
  for the protected full-B2 conventions. Their LP and lower-bound case
  analysis are not invoked by this submission.

## Part B. Algebraic disperser argument for the payload

This supplies the BK construction argument used in Part A, specialized to
F_2 and prime field degree. We use k for field degree, r for block count,
and t for the symmetric-polynomial order. All subspace dimensions are over
F_2. This part establishes the needed nonstandard algebra rather than
assuming a cited construction is a certified black box.

**Lemma.** Let K be the field of size 2^k with k prime. Suppose 1<=t<r and
k*t>=r(r+1). For any nonzero linear pi:K->F_2 put w_i=r+2-i and

    f(a_1,...,a_r)=pi(sum_(|I|=t) product_(i in I) a_i^(2^(w_i)-1)).

Then f is nonconstant on every affine subspace of K^r of dimension
strictly greater than kr/t+kt+r(r+1).

This is the binary case of BK's Theorem 7.3. We prove it below, including
the coefficient-isolation and block-selection statements it needs.

### B.1. Finite-field and polynomial preliminaries

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

### B.2. Image polynomials and their zero patterns

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

### B.3. A nonzero coefficient for independent subspaces

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

### B.4. Decompose an arbitrary affine subspace into blocks

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

### B.5. Select t high-dimension blocks

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

### B.6. Triangular mixing does not alter that coefficient

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

### B.7. No other symmetric summand cancels it

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
from Section B.3. M has positive degree, so the reduced polynomial is
nonconstant. Since R(K^r)=A, f is nonconstant on A. This proves the lemma.

### B.8. Exact application and source scope

For Part A, k is a prime supplied by the evaluator's length rule,
K=F_2[Z]/(Q) for an accepted prefix, and pi is its constant-coordinate
functional. The guard supplies r>t>=1 and k*t>=r(r+1). The zero-padded
payload embedding is injective, so every payload affine space of dimension
d_N=floor(kr/t)+kt+r(r+1)+1 maps to an affine space satisfying the lemma.
No field-encoding loss besides the counted prefix k is suppressed.

The source route is BK Lemmas 3.7,3.10,3.13; Theorem 6.1; Lemma 7.2;
Theorem 7.3 and Claim 7.4. We specialized the zero-pattern argument to
prime k and used (6) to give a direct proof of the final combinatorial
uniqueness. Only elementary finite-field algebra, root counting, linear
algebra and binary polynomial expansion remain outside this part. The
finite-field existence and polynomial-time evaluation arguments are
supplied in Part A. No efficient search for a hard truth table or for
an irreducible prefix is assumed.

## Part C. Gate lower bound for an affine disperser

We give the parameter-method argument of Mishra--Santhanam (TR26-260,
Sections 2-4) needed by Part A, with explicit bookkeeping. An affine
disperser of dimension d is nonconstant on every dimension-d affine
subspace, and hence on every larger affine subspace.

**Lemma.** For such f on n bits, every full-B2 circuit has size at least
(7/2)n-7d.

### C.1. Representation and potential

Start with an arbitrary circuit; it need not be a formula or minimum.
Allowing a final output literal for free only strengthens the lemma.
Absorb internal-edge negations in binary truth tables. Retained gates
stay charged until removed. Every nondegenerate two-input table is either
u+v+c or (u+b)(v+c)+e over F_2.

Keep M independent parameters u, a set P of q disjoint parameter pairs,
and an injective map to original input space

    Phi(u)=b+L*u+sum_({p,h} in P) c_(p,h)*p*h.

Coefficients are vectors over F_2; zero pair coefficients are permitted.
The computational representation is an affine system A*v=N*u+c followed
by acyclic B2 gates, with no feedback into the system. The auxiliary matrix
A is square and invertible. Each row has at most three distinct nonconstant
terms. The output is a constant or literal of a parameter, auxiliary value,
or gate. The invariant is that its value is f(Phi(u)) on every assignment.
Let s be number of equations plus acyclic gates, and a=M-q.

Degree counts distinct child rows and gates, not pins or the output
reference. An unpaired parameter has weight 0,6,7 when its degree is
0,1,>=2 respectively. Paired parameters have no individual weight. Set

    H=2s+3q+sum_(unpaired z) weight(z).

Initially M=a=n, q=0, s=|C| and H<=2|C|+7n. Always H>=0.
Fixing one member of every pair makes Phi affine and still injective on
M-q parameters. Its image therefore contains an affine space of dimension
a. Thus f(Phi(u)) is nonconstant whenever a>=d.

The following changes preserve injectivity and the displayed map form:

* Replace an unpaired x by w+ell(other parameters), with ell affine and
  w new. This is a bijection, preserving M,q,a.
* Replace an unpaired z by w+(x+b)(y+c)+e for distinct unpaired x,y,z,
  record pair {x,y}. This is a bijection; M is unchanged, q rises by 1,
  a falls by 1. No old quadratic term used z.
* Fix an unpaired x to an affine function of the others: M and a fall by 1.
  In actual reductions we use only constants or another parameter's literal.
* For pair {p,h}, fix p to a constant or to h+c, then remove p and the pair.
  Its product becomes affine. M,q each fall by 1, so a is unchanged.

The first two changes are not licenses to add free gates. Their exact
implementations appear below.

While a>d, an affine output is impossible. A constant output contradicts
the affine-space invariant. For a nonconstant affine output choose a
parameter with nonzero coefficient; if paired, first fix its partner,
preserving a. Solve for that now unpaired parameter to make the output
constant, losing only 1 in a. The remaining image still contains a
dimension-d affine space, a contradiction. Similarly, an output literal
of an AND on two unpaired parameters is impossible: fix a forcing input,
losing only 1 in a. This explicitly deals with terminal outputs.

### C.2. Legal deletions and pivots

An affine gate without an AND ancestor moves into the system by adding
its output as auxiliary and its defining equation as one row. The row has
at most three terms. The new auxiliary matrix is block triangular with
last diagonal entry 1. Size s is unchanged and parameter degrees cannot
increase. Move all such gates in topological order.

Four deletions preserve M,q,a and lower H by at least 1:

(D1) A constant or single-parent-literal gate is rerouted to that value
and removed. This covers unary/constant tables, fixed or repeated pins,
and an AND on affine parents with dependent linear parts. In the last case
the two affine values are constant or equal up to a constant, so the table
has one effective input. At most the one surviving parent can gain weight,
by at most 1 (6 to 7), against a saving of 2. Constant rerouting gives no gain.

(D2) A gate without children is removed unless it is the output. This
saves 2 and cannot increase parameter weight.

(D3) A row with at most two terms contains an auxiliary, since A is
invertible. Write it v=U+c, for U zero or one other value, substitute this
in all occurrences of v and erase the row and v. This preserves sparsity.
Adding the row to all rows using v, then expanding the determinant on
its column, proves the remaining auxiliary matrix invertible. Only U
can gain weight, by at most 1; the row removal saves 2.

(D4) An unpaired degree-one parameter x used only in a row E is pivoted
with an auxiliary Q that depends on x. Such Q exists because N_x is
nonzero and A^(-1)N_x is nonzero. Write Q=x+ell(other parameters).
Exchange their roles, and substitute x=w+ell in Phi, taking Q as the new
parameter w. The new auxiliary matrix is invertible: after multiplication
by A^(-1), its determinant is the Q-coordinate of A^(-1)N_x, namely 1.
No equation terms change. Now auxiliary x appears only in E and is not
the output while a>d (that would be affine). Delete x and E; determinant
expansion on its single-occurrence column proves invertibility again.
The old parameter weight 6 is replaced by at most 7 and one row saves 2.

The same affine pivot can be used at zero potential cost when x has degree
>=2 and Q is a parent of an AND gate. Q has a nonzero column in A and is
used in the gate, so the new parameter has degree >=2. Both weights are 7,
and s,M,q,a remain fixed. The affine expression is substituted only into
Phi: no dense equation or uncharged gate is inserted into the computation.

A quadratic pivot is used only as follows. Let G=(x+b)(y+c)+e have only
one child D=G+z+h, where x,y,z are distinct unpaired parameters and z has
degree 1. Substitute z=w+G(x,y)+h in Phi and pair x,y. Then D=w; reroute
its uses to w and delete D and now-unused G. No other equation or gate
used z, so none needs an extra computation. The changes satisfy

    Delta s=-2, Delta M=0, Delta q=1, Delta a=-1,
    Delta H <= -4+3+(7-6)-weight(x)-weight(y).

For degrees of x,y at least 2 this is <=-14. Arbitrary fan-out of D
still gives new parameter weight at most 7.

### C.3. Preliminary cases, in priority order

Work while a>d and restart after each reduction. Move eligible affine
gates before inspecting these cases. Each reduction has Delta a in
{0,-1,-2}, Delta H<=14*Delta a and strictly reduces M+s.

1. Perform any available D1-D4 deletion. Then Delta a=0 and Delta H<=-1.
   When none applies, gates have two distinct nonconstant parents and
   nondegenerate tables; ANDs on affine parents have independent linear
   parts; non-output gates have children; equations have three terms;
   and isolated unpaired parameters do not occur in equations.

2. Suppose an AND G on parameters has an unpaired isolated input z.
   Fix the other input w to force G constant and delete G. Parameter z
   becomes unused, though it may stay in Phi. If w is unpaired, its old
   weight is at least 6: Delta a=-1 and Delta H<=-6-6-2=-14.
   If w is paired, remove its pair and its partner becomes unpaired of
   weight at most 7: Delta a=0 and Delta H<=-3-6-2+7=-4.
   Rerouting a constant creates no occurrences.

   With cases 1-2 absent, no child can have both an isolated unpaired
   parameter and another parameter as parents. An equation is excluded
   by case 1; an XOR gate on parameters would already be in the system;
   an AND gate is case 2. This controls rerouting costs in later cases.

3. If a paired p is unused, fix its partner and erase the pair. The
   surviving p has weight zero: Delta a=0, Delta H<=-3.
   If a paired p has degree >=2, fix p, erase its pair and delete two of
   its distinct children. Gates have fixed pins; equations have at most
   two terms and are eliminated by D3. If one child disappears through
   earlier rerouting, it is already one of the two counted deletions.
   A term substitution in a short row never increases its term count.
   No originally isolated unpaired parameter co-occurs with p, by cases
   1-2. Only the now unpaired partner can gain weight, at most 7:
   Delta a=0, Delta H<=-3-4+7=0.
   When case 3 is absent all paired parameters have exactly one child.

4. If partners p,h share their one child B, impose p=h+c to make B
   constant and remove the pair. For a row v+p+h=b or an XOR any c
   works; for (p+alpha)(h+beta)+gamma choose c=alpha+beta+1.
   Delete B and reroute the constant. The surviving h is unused.
   Delta a=0, Delta H<=-3-2=-5.

5. Suppose paired p occurs in a row E_p, and its partner h has child B_h.
   These children are distinct by case 4. Fix h and delete B_h (D1 or D3).
   No previously isolated unpaired co-parent occurs there. Now p has at
   most its original row E_p as a child: B_h did not contain p, by p's
   uniqueness. If E_p survives, p has weight 6; pair removal and deletion
   have changed H by at most -3-2+6=1. Apply D4 to p,E_p for a further
   decrease of at least 1. Thus Delta a=0 and Delta H<=0. If E_p already
   vanished, count that deletion instead; p is unused and the estimate
   is better. Total s decreases at least 2.

6. Suppose an AND G on parameters has paired input p, partner h, and
   other input w. All paired parameters are isolated and absent from the
   affine system, so h's child B_h is a gate. Case 4 gives B_h!=G and
   w!=h; case 1 gives w!=p. Fix h and then fix w to force G constant.
   Delete G and B_h. The now unpaired p becomes unused.
   If w was unpaired, case 2 gives old weight 7, so Delta a=-1 and
   Delta H<=-3-7-4=-14. If w was paired, remove that distinct pair too;
   only its partner can gain weight, at most 7: Delta a=0 and
   Delta H<=-6-4+7=-3. No old isolated unpaired co-parent of B_h exists
   by cases 1-2, so no further increase has been omitted.

Restrictions reduce M and do not add gates or equations. Unused parameters
are retained in M until explicitly restricted; syntactic usage is not
substituted for affine dimension.

### C.4. Prepare a minimal AND

Assume cases 1-6 are absent. The acyclic part is nonempty, or the output
would be affine. Its first gate G is AND-type since all eligible affine
gates were moved. Its parents are affine values with independent linear
parts. If a parent Q is auxiliary, it depends on some parameter x occurring
in the system. Cases 1,3,5 imply that x is unpaired and has degree >=2.
Pivot Q with x at zero potential cost, as above. If the second parent is
auxiliary, its linear independence supplies a parameter other than the
first new parent for a second pivot. Recheck preliminary cases after each
pivot, ending the iteration if one applies. At most two pivots are needed.
Otherwise G now reads two unpaired parameters x,y of degree >=2.

If a child J also reads x, replace h(x,G(x,y)) by its composed binary table
on (x,y), rerouting its inputs to x,y. This is exactly one permitted B2
gate. Only y gains occurrences and its weight already is 7. Do the same
with x,y exchanged. After each such operation, move newly eligible affine
gates and recheck cases 1-6. Each rerouting removes a child of G, and none
of these operations creates a new one; this preparation is finite.
If no preliminary reduction applies, no child of G reads x or y. G is not
the output by the terminal argument and has at least one child by case 1.
Write u=deg(x)>=v=deg(y)>=2 and h=deg(G)>=1.

### C.5. Exhaustive remaining cases

7. If u=v=2 and the second child B of x,y is common, restrict (x,y) to
   a line making both G and B constant. An AND table is constant on
   three of the four inputs. A second Boolean table agrees on at least
   two of those three; the two points form an affine line. If B is a row,
   its auxiliary value is x+y+constant, so the same argument applies.
   The restriction fixes one parameter or identifies x with y+c.
   Delete G and B; the surviving parameter is unused. No weight grows
   under constant rerouting. Delta a=-1, Delta H<=-7-7-4=-18.

For later cases, if u=v=2 the other children B_x,B_y are distinct.

8. Let h_1 count G's children whose other parent is isolated and unpaired.
   Fix x to force G constant. Delete its u children and G's h children,
   a total of u+h distinct objects: no child of G reads x after preparation.
   Each row affected by x becomes short and each gate has a fixed pin.
   Any object deleted early through propagation is already a counted
   deletion; all remaining affected rows stay short. Thus the count is
   a lower bound even with coincidences and shared descendants.
   Removing x loses weight 7. A surviving old unpaired parameter gains
   weight only if it was isolated and its unique child is rerouted to it.
   Such a parameter cannot co-occur with x by cases 1-2. Only the h_1
   other parents of G's children can gain, each by at most 1. Hence

       Delta H <= -7-2(u+h)+h_1.

   If u>=3, this is <=-7-6-h<=-14. If u=2, also v=2. Since B_x!=B_y,
   and no child of G reads y, y gets no new occurrence and loses G:
   its weight drops by at least 1. Therefore in this case

       Delta H <= -12-2h+h_1.

   This is <=-14 for h>=2, or h=1,h_1=0. Delta a=-1.
   The sole uncovered configuration is u=v=2,h=h_1=1. Detect this
   before doing case 8: G's unique child D also reads an isolated
   unpaired z, distinct from x,y. Use case 9 or 10 instead.

9. If D is AND-type, write D=(G+c)(z+b)+e. Choose x0,y0 with G(x0,y0)=c
   and fix x=x0,y=y0. Such values exist because G is nondegenerate.
   Delete G,D,B_x,B_y, four distinct objects. G,D are rerouted to
   constants, and z becomes unused. B_x,B_y have no old isolated
   unpaired co-parent by cases 1-2. Thus no other weight grows and

       Delta a=-2, Delta H<=-7-7-6-8=-28.

   D is allowed to be the output. If its becoming constant leaves a>=d,
   this is an impossible state by the invariant; if a<d, termination is
   permitted. No extra successor of D is assumed.

10. Otherwise D is XOR-type. Apply the quadratic pivot of Section C.2.
    It deletes G,D, pairs x,y and replaces weight-6 z by a new parameter
    of weight at most 7. Thus

       Delta a=-1, Delta H<=-4+3+(7-6)-7-7=-14.

    This includes an output D, since its reference becomes the new
    parameter. All shared descendants are rerouted to that same value.

These cases are exhaustive because every nondegenerate B2 operation is
AND-type or XOR-type. The reduction to this local configuration was
proved, not assumed on the basis of finite pattern tests.

### C.6. Termination and gate count

Each successful case preserves the representation, never increases M or s,
and strictly decreases their nonnegative integer sum. Preparations are
finite (topological moves, at most two pivots, finitely many reroutings).
There are at most n+|C| iterations. By exhaustive coverage while a>d,
the process stops with a_T<=d. Summing Delta H<=14*Delta a gives

    2|C|+7n >= H_0 >= H_T+14(n-a_T) >= 14(n-d),

because H_T>=0. Thus |C|>=(7/2)n-7d. If initially n<=d the right side
is nonpositive, and the conclusion is immediate.

Rows are introduced only one-for-one from charged affine gates. Affine
pivots change roles without adding computation; the quadratic pivot
eliminates an existing two-gate pattern. Therefore the potential proof
bounds the original gate count, including under arbitrary sharing.

### C.7. What retuning these estimates can achieve

For a potential s+beta*q+alpha*I-sigma*ell, let I count used unpaired
parameters and ell isolated unpaired parameters; alpha>=sigma>=0,
beta>=0. If theta is the guaranteed decrease per lost unit of a, the
ordinary case estimates at (u,h,h_1)=(3,1,1) and (2,1,0) require

    theta<=alpha+4-sigma,       theta<=alpha+3+sigma.

Hence the derived coefficient theta-alpha is at most
min(4-sigma,3+sigma)<=7/2. Our H/2 has alpha=7/2,beta=3/2,sigma=1/2,
theta=7 and attains that value. This is only an optimization bound for
these local estimates, not a universal barrier to stronger circuit proofs.
