# Complete supplied proof and attribution

This three-file submission reproduces the coefficient-4 certificate-prefix theorem
from peer submission 52d8d43a-9613-45dc-9691-a406b7ae2c94,
whose sources are submissions 5e9fcf7f and f4f2c464. It does not claim a new
lower bound or reuse any judgment. The complete argument is supplied here.
Part A defines the family and proves eligibility and all-length accounting;
Part B supplies the algebraic disperser proof; Part C supplies the circuit
elimination proof. Section numbers are local to their part. Mentions in
Parts B/C of the source filenames identify their original provenance; all
mathematical content is contained in this proof.md.

# Part A: Family and all-length theorem

# A uniform certificate-prefix family: coefficient 4, exponent 1

## Claim and provenance

For the single total Boolean language defined below, with **N denoting the
entire input length**, we prove, for every integer N>=2^20,

    C_B2(F_N) >= 4N - R(N),
    R(N) = 4k + 9(floor(kr/t) + kt + r(r+1) + 2)
         <= 182 N^(4/5).

Here k,r,t are the deterministic length parameters below. In particular
liminf over all integer N of C_B2(F_N)/N >= 4. For any real epsilon>0 an
explicit sufficient threshold is

    N_epsilon = max(2^20, ceil((182/epsilon)^5)).

The family, evaluator and algebraic disperser proof are reused from
mooselumph's f4f2c464 package, with attribution; they are reproduced here
rather than treated as an external acceptance premise. Its family and
polynomial-time evaluator are unchanged. The BK algebraic argument is in
Part B of this proof. The gate argument in Part C of this proof (Sections 1-7) and
Sections 1-5 of this file are reproduced from public submission 5e9fcf7f
(ArthurOstapenko). That submission credits Mishra--Santhanam, ECCC
TR26-260, 2026-10-19, Sections 2-4, for the state, fair affine system,
affine exchanges and deletions. Its conclusion is 4n-9(d+1) for any
function nonconstant on affine spaces of dimension >d. This package is an
independent audit and repackaging of that theorem. It claims no new
lower bound.

The input-prefix wrapper pays for all k prefix bits. The whole family
need not be an affine disperser: a fixed-prefix restriction is one. The
lower bound uses arbitrary nonuniform circuits in the full binary basis,
with all sixteen truth tables, unrestricted fan-out and depth, repeated
pins, and exact computation on every input. Its normalization and terminal
cases are proved in Part C of this proof without adding uncharged gates.

## 1. One evaluator on every bit string

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

## 2. Exact prefix test and existence of hard prefixes

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

## 3. Uniform P and NP with empty witnesses

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

## 4. Explicit guards at every integer N>=2^20

Let X=N^(2/5). The integer B is the first positive integer strictly greater
than 2N^(3/5), so B<=2N^(3/5)+1. Bertrand's postulate gives k<2B. Hence
for every N>=2^20,

    2N^(3/5) < k < 4N^(3/5)+2 <= 5N^(3/5) < N/2.

The last inequality follows from X>=256>10. Thus m=N-k>N/2 and the
payload is nonempty. The ceiling defining r gives

    N^(2/5)/10 <= m/k <= r < m/k+1 < N/k+1
                                   < X/2+1 <= X.

In particular r>=3. Therefore t=ceil(sqrt(r)) satisfies 1<=t<r, and

    N^(1/5)/sqrt(10) <= sqrt(r) <= t <= sqrt(r)+1 <= 2N^(1/5).

The remaining guard follows without an asymptotic estimate:

    sqrt(r)(r+1) <= N^(1/5)(N^(2/5)+1) <= 2N^(3/5) < k,
    r(r+1) <= k sqrt(r) <= kt.

Consequently all length guards pass for **every** integer N>=2^20.
Small lengths, including N=0 and all failed guards, were already defined
by the evaluator. These bounds do not assume that N is a field block
multiple or that N lies on any selected subsequence.

## 5. Hard restriction and explicit remainder

Fix any integer N>=2^20. Section 2 supplies a monic irreducible Q of degree
k; fix the first k input bits to its coefficients. The remaining function
G_Q has exactly m=N-k actual inputs. Put

    d_N = floor(kr/t) + kt + r(r+1) + 1.

Part B of this proof proves, for prime k, r>t>=1 and kt>=r(r+1), that

    (a_i) -> pi(sum_(|J|=t) product_(i in J) a_i^(2^(r+2-i)-1))

is nonconstant on every affine subspace of K^r of dimension strictly
greater than kr/t+kt+r(r+1). Our pi is the constant-coefficient functional
and has pi(1)=1. The zero-padding payload embedding is injective and
linear, so G_Q is nonconstant on every payload affine subspace of dimension
d_N and every larger dimension. We use the slightly weaker hypothesis
of nonconstancy at dimensions strictly greater than d_N in the gate
lemma. The padding is internal to the evaluator and supplies no extra
circuit inputs.

Fixing prefix inputs of an arbitrary circuit for F_N gives a circuit for
G_Q with no additional gates. Constant nodes cost zero; any gates that
become unary, constant, or have repeated pins remain counted until their
legal removal. The full model treatment in Part C of this proof therefore
applies and gives

    C_B2(F_N) >= 4m - 9(d_N+1)
               = 4N - [4k + 9(floor(kr/t)+kt+r(r+1)+2)]
               = 4N - R(N).

The choice of Q in this lower-bound restriction is not advice to the
uniform evaluator: every circuit must handle every prefix, and the
algorithm simply reads and tests its input prefix.

Here is a numerical upper bound on the entire remainder. The strict
ceiling inequality gives kr<m+k=N. Using Section 4 and sqrt(10)<4,

    kr/t <= sqrt(10) N^(4/5) <= 4N^(4/5),
    kt <= 10N^(4/5),
    r(r+1) <= 2N^(4/5).

Hence d_N<=16N^(4/5)+1 and

    R(N) <= 20N^(3/5) + 144N^(4/5) + 18
          <= 182N^(4/5).

Every k-bit prefix cost is present in the term 4k; no linear input loss
is hidden. The ratio R(N)/N is at most 182N^(-1/5), which tends to zero.
For any epsilon>0 and every integer

    N >= max(2^20, ceil((182/epsilon)^5)) = N_epsilon,

we have R(N)<=epsilon N and C_B2(F_N)>=(4-epsilon)N. For a globally
defined remainder on smaller lengths one may set R(N)=4N there; the
nonnegativity of circuit size then gives the same inequality, while the
asymptotic ratio remains unchanged. Both finite coordinates in bound.json
belong to this one language and theorem.

## 6. Evidence and source scope

Part C of this proof proves the state invariants, all deletion and paired
rows, affine preparation, the three minimal-AND cases, every overlap used
in their accounting, coverage, termination and the final gate inequality.
Part B of this proof reproduces the coefficient-separation and block-mixing
proof from the source package with only an attribution addition. Sections
1-3 above reproduce its family, evaluator, prefix test and eligibility
argument. Section 4 and the explicit remainder estimates are re-derived.
There are no declared experiments or supplementary Python evidence files.
No claim of earlier acceptance or independent audit is a premise.

Sources, for attribution and version identification:

* mooselumph, CircuitBreaker package f4f2c464, `proof.md` Sections 1-3 and
  the payload construction in Section 5, and Part B of this proof in its entirety:
  source of the reused certificate-prefix family, evaluator and reproduced
  BK argument. The copied material is checked as mathematics here; its
  reported judge status does not establish the gate theorem.
* ArthurOstapenko, CircuitBreaker submission 5e9fcf7f: source of the
  ternary-pair domain, safe restrictions, potential W=s+5I+3q, the gate
  argument of Part C of this proof Sections 1-7, and Sections 1-5 of this
  file. Here it is checked as mathematics, and its reported judge status
  is not a premise.
* Mishra--Santhanam, ECCC TR26-260, 2026-10-19, Sections 2-4, for the fair
  semi-circuit representation, affine role exchange and D1-D4 machinery;
  https://eccc.weizmann.ac.il/report/2026/260/download/
  SHA256 38e9ff3cea4ac0399723756f00c2425036c438887ab5ae4abb10d09fa971a3bc.
* Ben-Sasson--Kopparty, author preprint 2010-03-12, Lemmas 3.7, 3.10, 3.13,
  Theorem 6.1, Lemma 7.2, Theorem 7.3 and Claim 7.4, for the reproduced
  algebraic construction;
  https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
  SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf.
* Li--Yang, ECCC TR21-023, 2021-02-20, Sections 1.1, 2.1, 2.4 and 3.3,
  for the pinned full-B2 conventions. Its lower-bound case analysis and
  LP are not invoked.


# Part B: Algebraic disperser

# Algebraic disperser argument for the payload

Reproduced from mooselumph's f4f2c464 CircuitBreaker package, with this
attribution addition only. Its construction is due to Ben-Sasson--Kopparty
in the versions and sections identified below and in proof.md. No reported
review outcome is used as a mathematical premise.

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

# Part C: Full-B2 circuit lower bound

# Ternary-pair gate elimination for an affine disperser

**Theorem.** Let d be a nonnegative integer. Suppose f:F_2^n -> F_2 is
nonconstant on every affine subspace of dimension strictly greater than d.
Every circuit for f in the full binary basis has at least

    4n - 9(d+1)

gates. This includes arbitrary sharing, repeated pins, constants, unary
tables, degenerate gates and a single terminal output.

Provenance. Sections 1-7 are reproduced unchanged from public CircuitBreaker
submission 5e9fcf7f (ArthurOstapenko), which credits Mishra--Santhanam,
ECCC TR26-260, 2026-10-19, Sections 2-4, for the fair affine system, role
exchanges and deletion machinery. Every step needed is proved in this
file; no external circuit theorem is invoked. Section 8 is an independent
audit added in this package.

## 1. State, semantics and potential

Start with any acyclic B2 circuit C for f. Count every retained gate,
including constant and unary gates, until its deletion. We allow a free
final constant or literal and negations on edges: an edge negation is
absorbed in the truth table of its child. Proving a lower bound in this
slightly more permissive output convention proves it in the pinned model.
The ten binary tables depending on both arguments are precisely the two
tables u+v+c and the eight tables (u+b)(v+c)+e over F_2.

A state consists of the following objects.

* There are m parameters u and an affine map Phi:F_2^m -> F_2^n that is
  injective on the entire parameter space, not just on the allowed domain.
* A set P of q disjoint unordered pairs of parameters is marked. A pair
  (p,r) has an excluded corner (e_p,e_r) in F_2^2. Its allowed set is the
  other three points. Dom is the Cartesian product of these sets and the
  unrestricted unpaired bits. Put a=m-q.
* The computation is a fair affine system A v = N u + c followed by an
  acyclic list of B2 gates. A is square and invertible. Each equation has
  at most three distinct nonconstant terms, after cancellation over F_2.
  Gates read parameters, auxiliary values v, or earlier gates; no gate
  feeds an equation. The output is a constant or a literal of an internal
  value. On Dom the output equals f(Phi(u)). Outside Dom no agreement is
  asserted. The affine system itself is uniquely solvable everywhere.

Let s be the number of equations plus acyclic gates. Degree counts
distinct child equations and gates, not pins or the output reference.
An equation is a child of every value appearing in it. Let I count the
unpaired parameters with positive degree, and set

    W = s + 5I + 3q.

Initially Phi is the identity, Dom=F_2^n, m=a=n, q=0, s=|C|, and
W_0<=|C|+5n. Always W>=0. Unused parameters are retained in Phi and m;
usage is never substituted for dimension.

For each pair, fixing p to its **safe value** 1-e_p leaves both values of
r allowed. Fix one member of each pair safely. The resulting affine
parameter subspace has dimension m-q=a, lies in Dom, and has an injective
affine image under Phi. Thus f(Phi(u)) is nonconstant on Dom whenever
a>d. This argument uses only affine dispersal.

Here and below Delta means new minus old. We run the algorithm only at

    a >= d+2.

At such a state the output cannot be affine in the parameters. On an
a-dimensional affine subspace just constructed, an affine Boolean
function is constant either everywhere or on a hyperplane of dimension
a-1>d. Its image would contradict the hypothesis on f. The same argument
excludes an output literal of an AND-type gate with affine parents: on
that subspace its factors are affine; a nonconstant factor can be set to
zero on a hyperplane, and if both factors are constant the output already
is constant. This also excludes parameter and auxiliary output literals,
including a parameter with no counted child. We will use the constant
output contradiction after restrictions even when the new a is d+1.

## 2. Legal operations and their implementation

Restrictions are made by composing Phi with an injective affine map tau
from the new full parameter space to the old one. The new Phi remains
affine and injective. Every restriction below maps its new product domain
inside the old Dom. No nonlinear reparameterisation of Phi is used.

Fixing an unpaired bit to a constant removes one parameter, preserves q,
and decreases a by one. Fixing one member of a pair safely removes that
parameter and the pair, leaves its partner unrestricted, and preserves a.
Restricting a paired three-point set to a line through two of its allowed
points removes one parameter and its pair, also preserving a. Such a line
is implemented by fixing one coordinate or identifying the coordinates
up to a constant. Each implementation substitutes only constants or
literals into the computational representation, so creates no gates or
longer equations.

We also exchange an unpaired parameter x with an auxiliary Q having solved
expression Q=x+ell(u_without_x), with ell affine. Substitute
x=w+ell in Phi, rename Q as parameter w and x as an auxiliary, and leave
the list of equations and gates otherwise unchanged. This is an invertible
affine coordinate change. Because x was an unrestricted bit, for each
fixed assignment of the other parameters it is a bijection of that bit.
Consequently its new domain is exactly the same product of marked
three-point sets and unrestricted bits, with w replacing x; no constraint
on a pair is changed, even if ell mentions paired bits.

For completeness, if Q is column j of A and N_x is the column of x,
then the new auxiliary matrix replaces column j of A by N_x. Multiplying
by A^(-1) gives an identity matrix whose column j is A^(-1)N_x. Its
determinant is that vector's j-coordinate, namely 1. The matrix stays
invertible; equation sparsity is unchanged. We do not substitute the
dense ell in the computation. The old x was used, and Q has a nonzero
column in A, so the new unpaired w is used. Other parameters have the
same children. Thus an exchange alone has Delta I=Delta s=Delta q=0.

Rerouting replaces every use, including the output, of a value that has
become a constant or a literal of another value by that constant or
literal. It never makes a previously unused parameter used: any surviving
parameter receiving an occurrence was a parent of the deleted object and
therefore was already used. This observation concerns I, which has the
same weight at every positive degree. Increasing a used parameter's
fan-out has no potential cost.

Rerouting a gate to an input creates no cycle: a gate input precedes it
topologically, and subsequent substitutions still use earlier surviving
values. Parameter and auxiliary values precede the acyclic part. Cycles
inside the affine system are permitted, with unique solvability proved
by invertibility. No nonlinear gate ever feeds that system.

## 3. Deletions and affine preparation

Each deletion below has Delta a=Delta q=0, Delta s=-1, and Delta I<=0,
so Delta W<=-1.

**D1.** Delete a constant or unary gate, a gate with a fixed or repeated
pin, or an AND-type gate whose affine parents have dependent linear
parts. Two dependent affine Boolean forms are constant or equal up to a
constant, so in the last case the output is a constant or a literal of
one of the existing parents. Reroute to it. This covers all sixteen truth
tables without retaining an uncharged extra operation.

**D2.** Delete a gate with no child unless it is the output. Its parents
can only lose occurrences. The terminal exclusions above prevent a
minimal AND or an affine value from being an exceptional output later.

**D3.** An equation with at most two terms has an auxiliary term, since
a zero row of A is impossible. Write it v=U+c, where U is zero or one
other value. Substitute this literal or constant in every occurrence of
v and delete v and the equation. Adding that row to the other rows using
v preserves determinant; afterwards v's column has one entry. Expanding
along that column proves the remaining matrix invertible. Term counts
cannot increase, and no unused parameter becomes used.

**D4.** If an unpaired isolated x has unique child equation E, its column
N_x is nonzero. Invertibility gives A^(-1)N_x!=0, so choose a dependent
auxiliary Q and exchange x with Q as in Section 2. The now auxiliary x
occurs only in E. It is not the output, by the affine-output exclusion.
Delete x and E without substitution. Expansion on its one-entry column
proves the remaining matrix invertible, and determines the removed x
uniquely from the surviving values. The old used unpaired x is replaced
by at most one used unpaired parameter Q, so Delta I<=0. This includes
the possibility that Q becomes unused after deleting E.

**M.** Move every affine gate whose parents are parameters or auxiliaries
into the affine system: replace J=b_1 t_1+b_2 t_2+c by its equation.
The auxiliary matrix becomes block triangular with an additional diagonal
1. One charged gate becomes one charged equation. Parameter usage does
not increase, so Delta W<=0 with m,s,q,a unchanged. Apply this in
topological order, then apply available deletions. It is finite because
each move removes a gate from the acyclic list.

The deletion arguments also justify the following batch counting rule.
After fixing a parameter, every original child gate has a constant pin
and every original child equation has at most two terms. Each can be
deleted. Literal substitutions never lengthen an equation or restore a
nonconstant fixed pin. If an object is removed during another scheduled
deletion, its deletion has already been counted; it is not counted twice.
Hence any specified set of distinct such objects supplies that many
deletions even with common parents, common descendants or cyclic affine
dependencies. Unscheduled simplifications can only help the estimates.

## 4. Preliminary paired rows

Exhaust D1-D4 and M before the following rows. Restart after a paired
row. All operations in these rows preserve a and do not increase W.

**P0: an unused paired p.** Fix its partner r safely and remove the pair.
The newly unpaired p remains unused. Thus Delta m=Delta q=-1,
Delta I<=0, Delta s<=0, and Delta W<=-3. Any simplification of r's
children can be left to the next iteration.

**P2: a paired p with at least two children.** Fix p safely and delete
at least two distinct children by D1 or D3. The batch counting rule
applies to equations as well as gates. Its partner r is the only newly
unpaired parameter and contributes at most one to I. Consequently

    Delta m=-1, Delta q=-1, Delta s<=-2, Delta I<=1,
    Delta a=0, Delta W<=-2-3+5=0.

With P0 and P2 absent, every paired parameter has exactly one child.

**PC: partners p,r share their child H.** If H is a gate, its two pins
are p,r. Among the three allowed assignments, two give the same Boolean
value at H; restrict the pair to the affine line consisting of those
two points. Then H is constant. If H is an equation, after D3 closure
it has exactly three terms and reads v+p+r=c for an auxiliary v. Apply
the same pigeonhole argument to p+r; v becomes constant. Delete H by
D1 or D3 and reroute the constant. The one surviving coordinate of the
line is unused because p,r had no other child. Thus

    Delta m=-1, Delta q=-1, Delta s<=-1, Delta I<=0,
    Delta a=0, Delta W<=-1-3=-4.

This choice of a line is inside the allowed three-point set; an arbitrary
old-MS diagonal restriction is not being assumed legal.

**PE: an isolated paired p occurs in an equation E.** Its partner r has
a unique different child H, by PC. Fix r safely and delete H. H may be
a gate or an equation. H does not contain p, since p's only child was E,
so its literal rerouting cannot introduce p elsewhere. The newly unpaired
p still occurs only in E, unless E has already been removed, in which
case that is a second counted deletion. Otherwise apply D4 to p,E.
This removes E and replaces p by at most one used unpaired auxiliary.
The two distinct objects H,E are deleted, and only this new unpaired
parameter can add to I. Therefore

    Delta m=-1, Delta q=-1, Delta s<=-2, Delta I<=1,
    Delta a=0, Delta W<=-2-3+5=0.

In particular, after these rows every paired parameter is isolated and
absent from the affine system. Every parameter appearing in an equation
is unpaired and used. All these statements are syntactic and do not
assume that auxiliary dependencies are acyclic.

## 5. Produce a minimal AND with parameter parents

Assume no earlier row applies and a>=d+2. If there were no acyclic gate,
the output would be affine, which is impossible. Choose the first gate
G in topological order. By M and D1 it is AND-type on two affine values
with independent linear parts. It is not the output by Section 1, and
D2 gives it a child.

**AR.** If a parent Q is auxiliary, choose an unpaired parameter x with
coefficient 1 in Q and exchange them. Such x exists: Q is nonconstant,
and every parameter occurring in the system is unpaired, by PE. A
parameter on which Q depends has a nonzero N column and so is used.
The exchange preserves the product domain, fairness and W by Section 2.
If G already has one parameter parent, choose x different from it. The
linear independence of the two parent forms guarantees such a coefficient
in Q. An invertible affine change preserves that independence. Thus at
most two exchanges make both parents parameters.

After each exchange recheck D1-D4, M and the paired rows; if a reducing
row applies, do it and end the iteration. Otherwise continue with G.
No gate operations or nonlinear ancestor relations have changed. In the
remaining state G has distinct parameter parents and at least one child,
and is not the output. There is no need to restructure children also
reading those parameters: I charges every used unpaired bit equally.

## 6. Eliminate this AND: all parameter cases

Write G=(v+b)(w+c)+e. Each paired parent, if present, is isolated. A
paired parent cannot have its partner as the other parent of G, because
that would be PC. The two pairs below, when present, are distinct.

**AP1: one paired parent p and one unpaired parent w.** Let r be p's
partner, with unique child H. H is a gate, by PE, and H!=G, by PC.
Fix r safely, then fix w to c to force G constant. These two restrictions
remove one pair and two parameters, so Delta a=-1. Delete H and G.
The newly unpaired p becomes unused; the old used unpaired w is removed.
No other previously unused parameter becomes used. Thus Delta I<=-1.

There is a third distinct deletable gate, including every overlap:
choose a child K of G. If a child K!=H exists, G's constant value makes
K deletable. It may also read w or some other scheduled value; it still
has a fixed pin. If all children of G equal H, H reads G and r. Once
both are fixed, H is constant. If H is not output, D2 closure in the
starting state gives a child J of H. Acyclicity makes J distinct from
G,H, and H's constant value makes J deletable. If H is output, after
these restrictions the output is constant on a domain with a'=a-1>d,
contradicting Section 1. Thus this last possibility cannot occur at a
loop state. Early deletion of K or J already supplies the same third
distinct deletion. H's other parent may be G, a descendant of G, or an
unrelated value; the preceding split covers all possibilities.

Consequently

    Delta m=-2, Delta q=-1, Delta s<=-3, Delta I<=-1,
    Delta a=-1, Delta W<=-3-3-5=-11 <= -9.

**AP2: both parents p,w are paired.** Let r,h be their respective
partners, with unique child gates H_r,H_h. Neither child is G, since
the only parents of G are p,w. Fix r and h safely, making p,w unpaired,
then fix w to c to kill G. Both pairs and three parameters disappear,
so Delta a=-1. The surviving p is unused and w is removed. The old
paired coordinates had no I charge, and no other unpaired usage can
appear, so Delta I<=0.

If H_r!=H_h, the three distinct gates G,H_r,H_h are deletable. If
H_r=H_h=H, that gate reads r,h. It cannot read G as a third pin. It
becomes constant after the two safe fixings. A child K of G exists and
is distinct from G,H. The fixed G pin makes K deletable. This gives
three distinct deletions also in the shared-child case. A child of G
may equal one of H_r,H_h in the first case; those three gates were
already distinct and suffice. Literal rerouting, shared descendants,
and an output among these children do not change the count. Hence

    Delta m=-3, Delta q=-2, Delta s<=-3, Delta I<=0,
    Delta a=-1, Delta W<=-3-6=-9.

**GP: both parents x,y are unpaired.** A nondegenerate AND-type truth
table has one exceptional corner and the same value at its other three
corners. Add the new pair (x,y), excluding exactly that exceptional
corner, and restrict Dom to its three allowed points. Phi does not
change at all. The pairs remain disjoint. This is a domain restriction,
not a bijective quadratic change of coordinates.

On the new Dom, G is constant. Delete G and all its k>=1 distinct child
gates. Each child now has a constant pin and is a unary table or a
constant, even if its other parent is x or y, another child of G, or a
shared descendant. Process them topologically; if one was already
deleted by propagation it is counted once. No equation reads G. The
parents x,y were each used because they read into G, and now lose their
unpaired I charges. All other unpaired parameters were already used
before any literal rerouting to them. Thus

    Delta m=0, Delta q=1, Delta s<=-(1+k), Delta I<=-2,
    Delta a=-1, Delta W<=-(1+k)-10+3=-8-k <= -9.

These three cases cover all parameter-parent types and all eight
AND-type operations; their polarities only select safe or killing
constants and the exceptional corner. Every fixing is legal: AP1 fixes
an old unpaired w arbitrarily, AP2 fixes w only after its partner has
been fixed safely, and GP pairs only previously unpaired bits.

## 7. Complete accounting, coverage and termination

For reference the upper bounds on changes used above are tabulated here.
The Delta I column bounds the total change of I, including every
rerouting and newly unpaired parameter. Negative bounds need no equality.

| Row | Delta s | Delta I | Delta q | Delta a | Delta W upper bound |
|---|---:|---:|---:|---:|---:|
| D1 | -1 | 0 | 0 | 0 | -1 |
| D2 | -1 | 0 | 0 | 0 | -1 |
| D3 | -1 | 0 | 0 | 0 | -1 |
| D4 | -1 | 0 | 0 | 0 | -1 |
| M | 0 | 0 | 0 | 0 | 0 |
| AR | 0 | 0 | 0 | 0 | 0 |
| P0 | 0 | 0 | -1 | 0 | -3 |
| P2 | -2 | 1 | -1 | 0 | 0 |
| PC | -1 | 0 | -1 | 0 | -4 |
| PE | -2 | 1 | -1 | 0 | 0 |
| AP1 | -3 | -1 | -1 | -1 | -11 |
| AP2 | -3 | 0 | -2 | -1 | -9 |
| GP, k>=1 | -(1+k) | -2 | 1 | -1 | -8-k |

For general nonnegative weights W=g*s+i*I+p*q, the corresponding
bounds are -g for each D row, 0 for M and AR, -p for P0,
-2g+i-p for P2 and PE, -g-p for PC, -3g-i-p for AP1,
-3g-2p for AP2, and -(1+k)g-2i+p for GP. At (g,i,p)=(1,5,3),
each reducing row satisfies Delta W<=9 Delta a and Delta a is 0 or -1.

While a>=d+2, do finite M preparation and the first available deletion
or paired row. If none applies, Section 5 supplies a minimal AND and at
most two affine exchanges, with the stated rechecks. Then exactly one
of AP1, AP2, GP applies. This proves universal coverage; it is not an
inference from finitely many tested local graphs.

Each reducing row preserves the state invariants, never increases m or
s, and strictly decreases the nonnegative integer m+s. In P0 even when
s is unchanged, m falls by one. GP decreases s, although m is unchanged.
Every other row either deletes an object or restricts parameters. The
preparations are finite: M removes acyclic gates and at most two AR
exchanges occur before a reducing row. There are therefore at most
n+|C| reducing iterations. Since coverage holds at every loop state,
the algorithm terminates with a_T<=d+1. If initially n<=d+1 the claimed
lower bound is already nonpositive.

Sum Delta W<=9 Delta a, including the nonincreasing preparations:

    |C|+5n >= W_0 >= W_T+9(n-a_T) >= 9(n-d-1).

Rearranging proves |C|>=4n-9(d+1). Each affine row was introduced
one-for-one from a charged gate, exchanges changed roles without extra
computation, and restrictions never added a gate. Thus this inequality
bounds the original gate count in the full B2 model, not a potential
mistaken for a gate count and not a formula or restricted-basis size.

## 8. Audit clarifications (added; no new premise)

These remarks restate why the steps most likely to be questioned are
sound. They add no row and change no constant.

* **Unused coordinates.** GP, PC and P0 can leave an unpaired parameter
  with no child. It stays in m and a and is never fixed. This is not a
  hidden loss. The circuit equals f(Phi(u)) on Dom and does not read
  such a coordinate, and Dom is a product with that bit unrestricted.
  So f(Phi(u)) is independent of it on the a-dimensional safe subspace of
  Section 1. If more than d coordinates were unused, fixing every other
  safe-subspace coordinate would give a constant restriction of dimension
  greater than d. Hence they never exceed d, and they are counted inside
  a_T<=d+1 at termination. No row makes an unused parameter used.
* **Three-point domains.** Only the a-dimensional affine subspace obtained
  by fixing one member of each pair to its safe value is used for
  dispersal. The three-point product set is never treated as having
  dimension log2(3) per pair.
* **Exchanges with paired bits in ell.** In an AR or D4 exchange, x is
  unpaired and unrestricted. For each assignment of the other parameters
  the substitution x=w+ell is a bijection of that bit, so Dom is
  unchanged as a product set.
* **Per-row dimension change.** Every row has Delta a in {0,-1}. No
  single row needs -18.
* **Accounting summary.** With W=s+5I+3q, the rows with Delta a=0 have
  Delta W<=0 and the rows with Delta a=-1 have Delta W<=-9. The tight
  rows are P2/PE (0), AP2 (-9) and GP with k=1 (-9). Summation gives
  |C|+5n>=9(n-d-1), that is |C|>=4n-9(d+1).
