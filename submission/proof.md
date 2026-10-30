# CircuitBreaker: a complete coefficient-four proof

This submission reuses the complete construction and proof from public Yukon
submission `33f76827-c1af-454c-bd5b-304495c3db1a`, commit
`8d2c1d9d8bd13915d64b84a4853b0f951d799908`, with explicit attribution.
Its family and algebraic argument originate in submission
`f4f2c464-3045-4c32-a6e4-32478c05f343`; its three-point domain method
originates in the note of `5e9fcf7f-8338-4729-9faa-a18994661af6`.
All three are declared dependencies. This package consolidates the arguments
into one proof, adds the explicit remainder in Part D, and spells out the
converse in the prefix test. It does not claim a new leading coefficient.
No previous judgment is used. Every mathematical premise needed for the
conclusion is proved here or is an explicitly identified elementary theorem.

Parts A, B, C and D form one proof. Numbered section references within a part
refer to that part; cross-part references identify the part explicitly.


## Part A. Fixed family and all-length theorem

This file gives the complete family, eligibility, parameter and all-length
assembly. Its circuit conclusion uses exactly the class lemma stated in
Section 5 and proved in Part C of this same file. A missing or unproved
class lemma is a gap in the main theorem; neither a prior review nor the
public coefficient-four idea note establishes that lemma.

## 1. Provenance and the precise theorem

The certificate-prefix family and the reproduced BK disperser argument are
reused from mooselumph's public submission
`f4f2c464-3045-4c32-a6e4-32478c05f343`, candidate commit
`8db0d3171ce7d19a48503ae34b5945331da55535`. They also appear in the
repository's `tests/fixtures/candidate-3_5` fixture. This reuse concerns
mathematical text and construction, not any previous validity judgment.
The public idea note of submission
`5e9fcf7f-8338-4729-9faa-a18994661af6` suggested the coefficient-four
elimination target. It supplies an idea rather than an
accepted theorem premise; the necessary argument is reproduced and must
stand on its own in Part C.

For the single language below, with N the entire actual input bit length,
the resulting theorem is

    C_B2(F_N) >= 4N - 52*N^(4/5) for every integer N>=2^20

for every sufficiently large integer N. Thus
liminf_(N->infinity) C_B2(F_N)/N >= 4. The coefficient and exponent refer
to this same family. The whole family is not asserted to be an affine
disperser: the property is established after a legitimate prefix restriction.

## 2. One exact evaluator at every length

Every integer choice below uses ordinary nonnegative integer arithmetic.
On a binary input x of length N, the deterministic evaluator A does this:

1. If N=0, return 0.
2. Find the least positive integer B with B^5>32*N^3. Find the least prime
   integer k>=B by testing successive integers using deterministic trial
   division. Here primality means k>=2 with no divisor from 2 through
   floor(sqrt(k)).
3. If k>=N, return 0. Otherwise set m=N-k, r=ceil(m/k), and let t be the
   least nonnegative integer with t^2>=r. If t>=r or k*t<r*(r+1), return 0.
4. Read the first k input bits b_0,...,b_(k-1) and form the monic polynomial

       Q(Z)=Z^k + sum_(j=0)^(k-1) b_j Z^j  over F_2.

   Check both Z^(2^k)=Z modulo Q and gcd(Q,Z^2+Z)=1. If either fails,
   return 0. Compute the first test by k repeated squarings modulo Q,
   starting at Z modulo Q; no integer of bit length 2^k is constructed.
5. Work in K=F_2[Z]/(Q). Append exactly rk-m zero bits to the payload,
   divide into r consecutive blocks of k bits, and regard block i as
   a_i=sum_(j=0)^(k-1) a_(i,j) Z^j. Thus the first bit in every block
   is its constant coefficient. The prefix is not part of any field block.
6. Put w_i=r+2-i and y_i=a_i^(2^(w_i)-1), for i=1,...,r. Return the
   constant coefficient of

       S_t(y_1,...,y_r)=sum_(I subset {1,...,r}, |I|=t) product_(i in I) y_i.

The exponents descend strictly with i. They are exactly those in
Part B. The output functional pi is the constant coefficient of
the unique representative of degree below k; it is F_2-linear and pi(1)=1.
In particular it is nonzero. We do not replace this functional with an
inconsistent trace or reverse the exponent order.

Step 6 has a polynomial-time implementation. Initialize E_0=1 and all
E_j=0 for 1<=j<=t. For i=1,...,r, first compute y_i by binary powering;
then update E_j <- E_j+y_i*E_(j-1) for j=min(i,t),...,1 in descending
order. Partitioning subsets according to whether they contain i proves
that after processing i, E_j is the order-j elementary symmetric
polynomial in the first i values. Hence E_t is exactly the specified sum.
No division is used. Each exponent 2^(w_i)-1 has w_i one-bits, so binary
powering uses O(w_i) quotient-ring multiplications.

All failed guards return the specified constant, and the valid branch is
also specified on every input. There is no hidden choice of a modulus,
truth table, advice string, witness or hard prefix in this evaluator.

## 3. Exact irreducibility test, existence, and uniform eligibility

We justify the special prime-degree test, including its converse. In an
algebraic closure of F_2 the polynomial H(Z)=Z^(2^k)-Z has derivative 1,
so it has exactly 2^k distinct roots. Its roots form a field: Frobenius
shows closure under addition and multiplication, and if a is a nonzero
root, (a^(-1))^(2^k)=a^(-1). Call this field K_0. The automorphism
sigma(a)=a^2 has sigma^k=id on K_0. An orbit has length dividing k, since
its period divides any exponent returning that element to itself.

For an orbit O, the monic polynomial product_(a in O)(Z-a) is fixed by
squaring all its coefficients. These coefficients are in K_0 and satisfy
c^2=c, so belong to F_2. This polynomial is irreducible over F_2: the roots
of any nonconstant F_2 factor would be a nonempty subset of O invariant
under sigma, and an orbit has no proper such subset. Conversely any monic
irreducible factor of H must consist of one orbit, by this construction
and irreducibility. Therefore every irreducible factor of H has degree
1 or k, because k is prime. Its degree-one roots are exactly 0 and 1.
There are 2^k>2 roots, so at least one orbit has length k. Thus an
irreducible monic degree-k Q exists for every prime k used by A.

The first test says Q divides the squarefree polynomial H. Consequently Q
is a product of distinct irreducible factors of degrees 1 or k. The gcd
test excludes exactly the two degree-one factors Z and Z+1. Since deg Q=k,
a passing Q is one irreducible degree-k factor. Conversely, for an irreducible degree-k Q, the quotient F_2[Z]/(Q)
is a field of 2^k elements. Lagrange's theorem for its multiplicative
group implies z^(2^k)=z for every element z, including zero. Applying
this to the residue of Z shows Q divides H. It has neither degree-one
factor (k>=2), so it passes. This proves the test in both directions. In particular every
accepted quotient is a field, and a passing prefix exists at every length
whose length guards pass. The existence proof makes no algorithmic search
for that prefix part of A.

All length comparisons are exact. Binary search finds B by fifth-power
comparisons with 32*N^3. Integer division computes r, and binary search
with square comparisons computes t. For N>=1, B>=3. Bertrand's theorem
in the form "for every integer B>1 there is a prime p with B<p<2B"
therefore applies, giving

    2*N^(3/5) < B <= k < 2B <= 4*N^(3/5)+2.                (1)

(The last inequality can be made strict, but is not needed.) Thus the
prime scan and its trial divisions take polynomial time in N. The
standard theorem's hypotheses are checked here; no prime-gap estimate
stronger than Bertrand is needed.

In the branch k<N, the degree-k quotient representatives use k bits.
Schoolbook multiplication followed by reduction modulo the monic Q takes
O(k^2) bit operations. Polynomial Euclid and k modular squarings take
polynomial time even on rejected prefixes. The powers in Step 6 use
O(r^2) multiplications in total, and the symmetric recurrence O(rt).
Here k,r,t<=N and all scalar integer operands have polynomial bit length.
Together with the prime search, a coarse fixed O((N+1)^6) bit bound
suffices for A on every string. The early guard covers the other branch.

Let L={x:A(x)=1}, and F_N be its characteristic function on {0,1}^N.
Take the one fixed witness bound W(N)=0 and verifier

    V(x,w)=1 iff w is empty and A(x)=1.

If x is in L, the empty witness is permitted and accepted. If a permitted
witness is accepted then A(x)=1, so x is in L. Both directions hold with
one polynomial-time verifier. In fact L is in P. Circuit inputs are only
the N bits of x; witness bits are never circuit inputs.

## 4. Every sufficiently large integer length is hard

Inequality (1) holds for every positive integer N. It gives
k=Theta(N^(3/5)) and k/N->0. Hence k<N eventually, and m=N-k=Theta(N).
Then the exact ceiling definitions give

    r=Theta(N^(2/5)),     t=Theta(N^(1/5)),
    r < (1/2)N^(2/5)+1.

In particular r->infinity, so t<r eventually. To check the other guard,
use t>=sqrt(r) and the displayed upper bound for r:

    r(r+1)/t <= sqrt(r)(r+1)
      <= sqrt((1/2)N^(2/5)+1)*((1/2)N^(2/5)+2).

The last quantity divided by N^(3/5) tends to 1/(2*sqrt(2)), whereas
k/N^(3/5)>2. It follows that k*t>=r(r+1) for every sufficiently large
integer N. These are uniform inequalities on all integers, not estimates
along a prime or sparse subsequence. Finite guarded lengths require no
special circuit claim.

For each such N fix any accepted irreducible prefix Q; Section 3 proves
that one exists. With these k original input bits fixed, the remaining
function G_Q has exactly m=N-k Boolean inputs. Put

    d_N = floor(kr/t) + kt + r(r+1) + 1.                  (2)

The proved lemma in Part B applies because k is prime, 1<=t<r,
and kt>=r(r+1). It says the field function in Step 6 is nonconstant on
every affine subspace of K^r with dimension greater than
kr/t+kt+r(r+1), for the exact nonzero functional pi already specified.

Appending the fixed zero bits defines an injective F_2-linear map from
F_2^m into K^r. The image of each affine input subspace has exactly its
original dimension. By (2), G_Q is therefore nonconstant on every affine
subspace of dimension d_N (and on all larger ones). This is the dimension-d_N
affine-disperser property used below. The rk-m<k padding bits are fixed
zeros and do not contribute to m or N. They lose no affine dimension;
we apply the gate lemma at m, never at the artificially enlarged rk.

Finally, kr<=m+k=N gives kr=O(N). From the uniform estimates for k,r,t,

    floor(kr/t)=O(N^(4/5)),   kt=O(N^(4/5)),
    r(r+1)=O(N^(4/5)),       d_N=O(N^(4/5)).             (3)

In particular d_N/m->0, so 1<=d_N<=m eventually. No vacuous
high-dimensional disperser condition is being used to obtain a positive
leading coefficient.

## 5. Full-B2 model, gate lemma, and the actual-length bound

The sole class lower bound needed from Part C is:

    For 1<=d<=m, if g:{0,1}^m->{0,1} is nonconstant on every
    dimension-d affine subspace, then every exact deterministic
    full-B2 circuit C for g has |C| >= 4m-9(d+1).           (4)

Here full B2 means the pinned general model with all sixteen binary truth
tables, acyclic circuits, one output, arbitrary sharing and fan-out, two
pins per gate including repeated pins, free input/wire/output selection,
and the gate counting and constant-substitution normalization in rules.md.
Unary, constant and negated-projection retained gates are covered by (4),
as are terminal outputs. The complete proof of (4), rather than any
citation to an unproved four-coefficient theorem, is required evidence in
Part C. No formula, monotone or direct-sum transfer is assumed.

Take any B2 circuit for F_N. Fix its k prefix inputs to the coefficients
of the chosen Q. It now computes G_Q exactly on all m payload inputs.
Substitution adds no gates: it only attaches fixed constants to former
input pins, and the no-larger normalized circuit is within the model of
(4). Constants introduced here are dealt with by (4)'s normalization;
no free computational basis is being appended. Therefore (2)--(4) give

    C_B2(F_N) >= 4m-9(d_N+1)
              = 4N - [4k+9(d_N+1)].                     (5)

All k certificate bits are thus charged at their actual input length.
Restriction does not need an efficient rule selecting Q: an exact circuit
for F_N must work on every prefix. The evaluator itself reads and tests
its input prefix and was fixed before any asymptotic quantifier.

By (1) and (3), R(N)=4k+9(d_N+1)=O(N^(4/5))=o(N) eventually.
Extend R arbitrarily over the finitely many failed-guard lengths. For any
real epsilon>0 choose N_epsilon large enough that all the preceding
conditions hold and R(N)/N<=epsilon for every integer N>=N_epsilon.
Then (5) yields

    C_B2(F_N) >= (4-epsilon)N for every integer N>=N_epsilon.

This proves precisely the claimed liminf and the coefficient-four,
exponent-one polynomial_log guarantee with both logarithmic exponents zero.
In the equivalent relative-epsilon version rendered by `expressions.py`,
replace epsilon here by 4*epsilon. There is one fixed family throughout.
The explicit threshold and remainder are proved in Part D below. No finite
experiment or stronger exponent is asserted.

## 6. Source scope and supplied evidence

Part B reproduces the needed algebraic argument specialized to
binary fields of prime degree, including coefficient isolation and the
triangular mixing argument. Its source attribution is Ben-Sasson and
Kopparty, *Affine Dispersers from Subspace Polynomials*, author preprint
2010-03-12, Lemmas 3.7, 3.10, 3.13, Theorem 6.1, Lemma 7.2, Theorem 7.3
and Claim 7.4:
https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
(SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf,
as recorded in the reused package). The reproduced argument, not retrieval
of that link, supplies its mathematical premises.

The predecessor's parameter-method background is Mishra--Santhanam,
ECCC TR26-260, 2026-10-19, Theorem 4.1 and Sections 2--4. Its earlier
7/2 conclusion does not imply (4). The proof of the improved gate lemma
is supplied in Part C. Li--Yang, ECCC TR21-023,
2021-02-20, Sections 1.1, 2.1, 2.4 (footnotes 3--4) and 3.3 identify the
model pinned by rules.md; no Li--Yang LP bound is used here.

All supporting arguments are in this one file. The package is proof-only: `experiments` and `checks` are empty. There is
no declared computation whose output is a mathematical premise. The two
main claims are identical canonical renderings of `bound.json`, and their
identical sorted evidence lists contain all supplied proof and bound files.

## Part B. Algebraic disperser argument for the payload

This text reuses, with attribution, the algebraic proof in public submission
`f4f2c464-3045-4c32-a6e4-32478c05f343` by mooselumph, candidate commit
`8db0d3171ce7d19a48503ae34b5945331da55535`. The repository fixture
`tests/fixtures/candidate-3_5/disperser.md` has the same supplied text.
The argument is reproduced here as mathematical evidence; no previous
acceptance, score, or review is used as a premise. The coefficient-four
change is entirely in the gate lemma in Part C, not this algebraic lemma.

This supplies the BK construction argument used in proof.md, specialized to
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

Here and below, exponent reduction fixes exponent zero; a positive exponent
is reduced to the member of 1,...,q-1 congruent to it modulo q-1.
This is the reduction modulo X^q-X and preserves the value also at zero.

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
bijections is exactly the permanent. Each full-weight contribution uses distinct residues in each variable;
a carry strictly reduces that variable's weight and cannot contribute. In characteristic two permanent
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
algebra and binary polynomial expansion remain outside this part. The
finite-field existence and polynomial-time evaluation arguments are
supplied in proof.md. No efficient search for a hard truth table or for
an irreducible prefix is assumed.

## Part C. A coefficient-four gate bound from three-point pair domains

## Statement and credit

**Theorem.** Let `d` be a nonnegative integer, and let
`f : F_2^n -> F_2` be nonconstant on every affine subspace of dimension at
least `d`. Every deterministic acyclic, one-output circuit for `f` over
all sixteen binary Boolean truth tables, with arbitrary sharing, has at
least

    4n - 9(d+1)

gates. Input variables, wires, and output selection have zero cost; each
retained gate, including a trivial gate, has cost one. This theorem is a
class theorem; eligibility and an all-length family instantiating it must
be supplied separately.

The three-point-pair-domain idea and the potential `s+5I+3q` are taken
from the public candidate note titled "Coefficient 4: three-point pair
domains give 4n - 9(d+1) for every affine disperser", associated with
submission `5e9fcf7f-8338-4729-9faa-a18994661af6`. That note credits GPT 6
Astra, GPT 6.1 Sol, and Claude Opus 5.5 for the idea, accounting, and
writing. The proof below reconstructs and proves the claim independently;
neither its submission status nor its reported audits are premises.
It uses the fair affine-system/parameter framework of Mishra--Santhanam,
ECCC TR26-260 (19 October 2026), Sections 2--4. Everything needed from
that framework is proved below. Three-point fibers were also used in
weighted gate elimination by Golovnev--Kulikov, ECCC TR15-170 (ITCS 2016);
no theorem from that work is a premise here.

All arithmetic and affine maps below are over `F_2`.

## 1. Representation and the domain invariant

We work in a relaxation of the circuit model, starting from the given
circuit without increasing its size. There are `M` independent Boolean
parameters, a collection `P` of `q` disjoint pairs of these parameters,
and a domain `D` in parameter space. For every pair `{p,h}` exactly one
corner `(alpha,beta)` is excluded; the allowed values of that pair are
the other three points of `F_2^2`. Unpaired parameters range freely.
The domain is the product of these pair domains and the free bits.

There is an affine map

    Phi : F_2^M -> F_2^n

which is injective on the **entire** parameter cube, not just on `D`.
There are auxiliary variables `v` specified by an affine system

    A v = N u + c,

where `A` is square and invertible. Each equation has at most three
distinct nonconstant terms after cancellation. An acyclic collection of
binary Boolean gates follows the system. Gate inputs can be parameters,
auxiliaries, or earlier gates. No gate is an input to an equation.
The output is a constant or a literal of a parameter, auxiliary, or
gate. Permitting the final output literal for free is a relaxation of
the stated model and therefore strengthens the lower-bound theorem.

The representation is required to compute `f(Phi(u))` for `u in D`.
Its behavior outside `D` is unrestricted. Auxiliaries nevertheless have
unique affine values on the entire cube because `A` is invertible.
The representation's size `s` is its number of equations plus retained
gates. Initially `M=n`, `Phi` is the identity, `q=0`, there are no
equations, and `s` is the original circuit size.

A value's degree is the number of distinct equations in which it occurs
as a nonzero term, plus the number of distinct gates having it as a
parent; repeated pins count once. Equations have no distinguished
output term for this definition. The output reference is not counted.
Let `I` be the number of **used unpaired
parameters**, namely unpaired parameters of positive degree. Put

    a = M-q,                 W = s + 5I + 3q.

Always `W>=0`, and initially `W<=s_initial+5n`. We retain unused
parameters in `M` unless a specified restriction removes them; syntactic
nonuse does not by itself justify a dimension decrease or increase.

For a pair excluding `(alpha,beta)`, fixing `p=1+alpha` leaves `h`
free, and fixing `h=1+beta` leaves `p` free. We call these the safe
values. Fixing one member of every pair safely gives an affine subspace
of the parameter cube, contained in `D`, of dimension `M-q=a`.
Its image under `Phi` is an affine subspace of exactly that dimension.
In particular, a constant output is impossible whenever `a>=d`.

We perform reductions only while

    a >= d+2.                                      (guard)

Two useful terminal exclusions follow directly, without requiring that
the following affine slices themselves be implemented as reductions.
On a safe affine subspace of dimension `a`, any affine output is
constant on a subspace of dimension at least `a-1`. Thus a constant or
affine output, including a literal of a parameter or auxiliary, is
impossible under the guard. An AND-type operation on two affine values
also has a constant affine slice of codimension at most one: write it
`(L+b)(R+c)+e`. If `L` is nonconstant, set `L=b`; if `L` is constant,
the operation is either constant or affine in `R`, to which the previous
argument applies. Consequently a gate of this form with affine parents
cannot be the output under the guard. Each resulting slice has dimension
at least `a-1>=d`, and `Phi` preserves its dimension.

## 2. Boolean tables and legal simplifications

The sixteen two-input tables consist of two constants, four literals
of one input, two tables `U+V+e`, and eight tables

    (U+b)(V+c)+e.

These last eight are called AND-type. Their majority value is `e`,
and their sole exceptional input is `(1+b,1+c)`. Constant nodes introduced
by substitution are free. Negations are absorbed into descendant truth
tables or into affine equation constants; a terminal negation is covered
by the relaxed output-literal convention. No retained negation or
projection gate is silently omitted from the size count.

Here are four simplifications. Each deletes at least one gate or
equation, keeps `M,q,a` unchanged, and does not increase `I`.

**D1: a degenerate gate.** A gate which is a constant or a literal of
one actual parent is rerouted to that value and deleted. This includes
constant or unary tables, fixed pins, and repeated pins. It also applies
when both parents are affine values with dependent linear parts: they
are constant or equal up to a constant, so composing the two-input table
gives a constant or a literal of a nonconstant actual parent. A rerouted
parameter was already used in the deleted gate, so it cannot newly
enter `I`.

**D2: a dead gate.** Delete a gate of degree zero unless it supplies the
output. This cannot increase parameter usage.

**D3: a short equation.** An equation with at most two distinct terms
contains an auxiliary, since otherwise its row in `A` would be zero.
Solve it as `v=U+c` with `U` either one other value or zero. Substitute
in all other equations, gates and the output, and delete this equation
and `v`. Substitution replaces a single term by at most one term, so
equation lengths do not increase. Ordinary row elimination followed by
expansion along the `v` column proves that the remaining auxiliary
matrix is invertible. If `U` is a parameter, it was already used in the
deleted equation, so `I` does not increase.

**D4: an isolated unpaired parameter in an equation.** Suppose an
unpaired parameter `x` has just one child, an equation `E`. Some
auxiliary `Q` has coefficient one on `x` in its affine solution, because
the column `N_x` is nonzero and `A^{-1}N_x` is nonzero. Exchange the
roles of `x` and `Q`, leaving the equations' terms unchanged. Taking
`w=Q` as the new unpaired parameter, write

    Q=x+ell(other parameters),   x=w+ell(other parameters),

and make this affine substitution in `Phi`. The coordinate change is
a bijection on the full cube. It preserves the product domain: all
paired coordinates stay fixed by the change, and for every choice of
them the old free coordinate `x` is exchanged bijectively with the new
free coordinate `w`. Hence `Phi` stays affine and injective.

The new auxiliary matrix replaces the column `A_Q` by `N_x`. Its
determinant is nonzero because the `Q` coordinate of `A^{-1}N_x` is
one. The old parameter `x`, now auxiliary, occurs only in `E`; it is
not the output, since that would be affine under the guard. Delete it
and `E`. Expanding the determinant on its single-occurrence column
proves invertibility of the remaining matrix. No gate or equation has
been added. Before deletion a used unpaired parameter was exchanged
for another used unpaired parameter: `Q` occurs in some equation
because the original column `A_Q` is nonzero. Thus `I` cannot increase.

All four simplifications give `Delta W<=-1` and `Delta a=0`.
The deletion arguments count original distinct objects. In particular,
if propagation already deletes a second designated object, it counts
as deleted; if it survives, a short equation stays short under D3 and
a fixed-pin gate stays degenerate until it is deleted. These observations
will justify simultaneous deletion counts below.

An affine gate whose parents are parameters or auxiliaries can instead
be moved into the system: add its defining equation and make its output
a new auxiliary. The new matrix is block triangular with a final
diagonal entry one; the equation has at most three terms. This removes
one gate and adds one equation, so `s` does not change, and `I` cannot
increase. Repeat in topological order whenever possible. This is a
transfer of one charged object to another, not an additional free
computation.

The role exchange used in D4 is also legal without the final deletion,
whenever an unpaired `x` has coefficient one in the solution for an
auxiliary `Q`. It preserves `s,M,q,a,I`, for the same reasons. We call
this an affine role exchange. The possibly dense affine expression is
substituted only in `Phi`; the computational equations retain their
original terms. No dense affine gate is inserted for free.

## 3. Clean up the paired parameters

At the start of each iteration, move eligible affine gates and perform
any available D1--D4 simplification. A simplification is a successful
iteration, after which we restart. Otherwise try the following pair
cases in order, restarting after any successful case.

**P0: an unused member.** If paired `p` is unused, fix its partner `h`
safely, and remove `h` and this pair. The surviving free `p` is unused.
Thus `Delta M=Delta q=-1`, `Delta a=0`, `Delta I<=0`, and
`Delta W<=-3`.

**P2: a member of degree at least two.** Fix such a member `p` safely,
remove it and the pair, and delete at least two distinct children of
`p`. A gate child has a fixed pin and is deleted by D1; an equation
child becomes short and is deleted by D3. The other member becomes an
unpaired parameter, so it is the only possible new member of `I`.
Consequently

    Delta s<=-2,  Delta I<=1,  Delta q=-1,  Delta a=0,
    Delta W<=-2+5-3=0.

When neither P0 nor P2 applies, every paired member has exactly one
child.

**PC: partners have the same child.** Suppose partners `p,h` have the
same unique child `B`. If `B` is a gate, its Boolean table has the same
value on two of the three allowed points, by the pigeonhole principle.
Those two points form an affine line contained in the allowed domain.
Restrict the pair to that line. The line is expressible by fixing one
coordinate or identifying the coordinates up to a constant; remove one
parameter and the pair. Now `B` is constant and can be deleted. If `B`
is an equation, its other term is an auxiliary and it has form
`v+p+h=c`. Among the three allowed points two have the same `p+h`;
on their line `v` is constant, and D3 deletes the equation. The
surviving parameter is unused. Therefore

    Delta s<=-1,  Delta I<=0,  Delta q=-1,  Delta a=0,
    Delta W<=-4.

The restricted `Phi` is still affine and injective. This includes the
case that the deleted gate supplies the output: if the resulting
constant output conflicts with the domain invariant, the supposed
configuration simply cannot occur under the guard.

**PE: a paired member occurs in an equation.** Let `p` have sole child
`E`, an equation, and let its partner `h` have sole child `B`. By PC,
`E` and `B` are distinct. Fix `h` safely and remove it and the pair.
Delete `B` using D1 or D3. It contained no `p`, since `p` had only
child `E`; none of these substitutions creates new occurrences of `p`.
If `E` survives, `p` is an unpaired parameter isolated in `E`, so D4
deletes `E`. If `E` has already disappeared, the second distinct deletion has
already occurred. Propagation can make `p` used elsewhere, but no other
previously unused unpaired parameter becomes used. The only possible increase in
`I` is one, for `p` or its replacement in D4. Thus

    Delta s<=-2,  Delta I<=1,  Delta q=-1,  Delta a=0,
    Delta W<=0.

If none of these pair cases applies, every paired parameter has exactly
one child, that child is a gate, and partners have distinct children.
In particular no paired parameter occurs in any equation. All auxiliary
values therefore depend only on unpaired parameters.

## 4. Expose a first AND-type gate

Assume none of the preceding reductions applies under the guard. The
acyclic gate part is nonempty, since otherwise the output would be
affine. Its first gate `G` has only parameter or auxiliary parents.
It is neither degenerate nor affine, by D1 and the topological moves,
and hence is AND-type. Its parents' affine linear parts are independent,
or D1 would apply.

If one parent `Q` is auxiliary, its affine form depends on an unpaired
parameter `x`: it is nonconstant, and paired parameters occur in no
equation. Use an affine role exchange to make `Q` an unpaired parameter.
If the other parent is still auxiliary, linear independence guarantees
that its form involves some unpaired parameter other than the first
new parent; exchange using that parameter. Thus at most two exchanges
make both parents of `G` parameters. Paired coordinates and all their
occurrences are unchanged. Newly exposed degeneracies or eligible affine
gates can be simplified before continuing; otherwise these exchanges
have zero cost and the same `G` has two distinct parameter parents.

The terminal exclusion in Section 1 shows that `G` is not the output.
D2 consequently gives at least one gate child of `G`. Gates cannot
have equation children. No preparation which reroutes children of `G`
away from its input parameters is needed.

Exactly three cases remain, according to the pair status of its parents.

**GP: both parents are unpaired.** Write

    G=(x+b)(y+c)+e.

Restrict `(x,y)` to the three-point majority fiber, excluding
`(1+b,1+c)`, and record this as a new pair. No parameter or coordinate
map is changed. On the new domain `G=e`, so delete `G` and each of
its original children by constant substitution and D1. There is at
least one child, and it is distinct from `G`. The previously used
unpaired parameters `x,y` cease to contribute to `I`. Ordinary
rerouting cannot create another used parameter, since any parameter
rerouted from a deleted gate was already used there. Hence

    Delta s<=-2,  Delta I<=-2,  Delta q=1,  Delta a=-1,
    Delta W<=-2-10+3=-9.

The new computation need only be correct on the new domain. In
particular `G` is not asserted constant on the full parameter cube.
The coordinate map remains affine and injective on that full cube.

**AP1: exactly one parent is paired.** Let the parents be paired `p`
and unpaired `w`, and let `h` be the partner of `p`. Let `H` be the
unique child of `h`; by the pair cleanup `H` is a gate distinct from
`G`. Fix `w` to the value forcing the AND-type `G` to its majority
constant. Also fix `h` safely in its pair. Remove `w,h` and that pair.
The surviving `p` is free and becomes unused after `G` is removed.
Delete `G` and `H`: the former is constant, and the latter has its
`h` pin fixed. These are two distinct deletions, even if `H` is a
child of `G` or they share descendants. No third deletion is needed.
Removing the used unpaired `w` loses one member of `I`; the only
newly unpaired coordinate is unused `p`. Thus

    Delta s<=-2,  Delta I<=-1,  Delta q=-1,  Delta a=-1,
    Delta W<=-2-5-3=-10<=-9.

This slightly weaker estimate than the public note's `-11` avoids any
need to count a successor of `H` and suffices for the theorem.

**AP2: both parents are paired.** They are in different pairs, since
PC excludes partners sharing `G`. Write the pairs as `{p,h}` and
`{r,k}`, with `G` reading `p,r`. Let `H,K` be the unique children
of `h,k`, respectively. They are gates, each distinct from `G`.
Choose a valid point of the first pair with `p` equal to a forcing
value for `G`, and fix that whole pair to this point. Such a point
always exists: either value of `p` occurs in the three-point domain.
The second member `r` becomes unused when the now constant `G` is
deleted; fix its partner `k` safely, so `r` remains a free unused
parameter. These restrictions remove three parameters and two pairs.

Delete `G,H,K`. If `H!=K`, these are three distinct gates. If `H=K`,
their common gate reads `h` and `k`: it has two distinct pins and both
parameters have it as their unique child. Therefore this common gate
cannot also read `G`. Since `G` has a child, choose one such child `J`;
then `G,H,J` are three distinct gates. The pins `h,k` are fixed, so
`H` is deletable, and the `G` pin of `J` is constant, so `J` is
deletable. This proves at least three distinct deletions in every
overlap pattern. The only surviving newly unpaired parameter is `r`,
which is unused; no old unpaired parameter becomes used. Hence

    Delta s<=-3,  Delta I<=0,  Delta q=-2,  Delta a=-1,
    Delta W<=-3-6=-9.

All these restrictions preserve the product form of the domain and
affine injectivity of `Phi`: a safe single-coordinate fixation removes
one parameter and one pair, a whole-pair fixation to a valid point
removes two parameters and one pair, and an unpaired fixation removes
one parameter. In particular their displayed changes in `a` are exact.

## 5. Coverage, termination, and the inequality

Every successful reduction satisfies

    Delta a in {0,-1},        Delta W <= 9 Delta a.

The zero-dimensional-loss reductions are D1--D4 and P0, P2, PC, PE;
the unit-loss reductions are GP, AP1, AP2. Each preserves all
representation invariants, never increases `M` or `s`, and strictly
decreases the nonnegative integer `M+s`. Each preparation is finite:
affine moves transfer existing gates into equations without introducing
new gates, and at most two role exchanges expose the selected first
AND-type gate. If a simplification becomes applicable during a
preparation, perform it as the successful reduction and restart.

Coverage is exhaustive. After simplifications and pair cleanup, a
nonempty acyclic part has a first nondegenerate nonaffine gate, whose
parents can be exposed as parameters. It cannot be the output, and
the status of its two parents is exactly one of GP, AP1, AP2.
If the acyclic part were empty, the output would be affine, which is
excluded by the guard. Output constants or literals produced by a
reduction are handled by the same domain invariant; if such an output
would conflict with that invariant, the purported pre-reduction state
was impossible. No estimate presumes an extra successor of an output.

There are at most `n+s_initial` successful iterations. Consequently the
process must leave the guard, with a final `a_T<=d+1`. If `n<=d+1`
initially, the asserted lower bound is already nonpositive and is
immediate. Otherwise sum the displayed inequalities over the process:

    s_initial + 5n >= W_0
                      >= W_T + 9(n-a_T)
                      >= 9(n-d-1),

since `W_T>=0`. Therefore

    s_initial >= 4n - 9(d+1),

as claimed.

The size in this inequality is the original number of charged gates.
Every equation was introduced one-for-one from an existing charged
gate. Role exchanges change the representation and coordinate map
without adding computational objects. Every other step is a restriction
or a deletion. The argument thus covers all sixteen tables, constants,
unary and repeated-pin gates, terminal outputs, arbitrary fan-out,
shared descendants, and unrestricted nonuniform circuits in the stated
full-B2 model.

## Part D. Explicit threshold, remainder and quantifiers

We now replace Part A's asymptotic estimates by a finite bound. This part
uses only inequalities on integers and real fifth roots, not numerical
experiments. Let N be any integer at least 2^20 and put z=N^(1/5), so z>=16.
The definitions in Part A and Bertrand's theorem give

    2*z^3 < k < 4*z^3+2 <= (65/16)*z^3.

Indeed 2<=z^3/16 for z>=16. In particular

    k/N < 65/4096 < 1/32,       m=N-k > (31/32)*N,
    r=ceil(m/k) >= m/k > (31/130)*z^2,
    r < N/k+1 < (1/2)*z^2+1 <= (129/256)*z^2.

The last ceiling bound is valid whether or not k divides m. The lower
bound implies r>61, so t=ceil(sqrt(r))<r. Moreover

    (12/25)*z < sqrt(r) <= t <= sqrt(r)+1 <= (13/16)*z.

For the lower estimate use (12/25)^2 < 31/130. For the upper estimate use
129/256 < (3/4)^2 and 1<=z/16. Consequently

    k*t > (24/25)*z^4,
    r*(r+1) <= (129/256)*(130/256)*z^4
             = (16770/65536)*z^4 < (24/25)*z^4.

Thus all length guards pass at every integer N>=2^20. The prefix test
still rejects invalid moduli, as it must; a passing modulus exists by
Part A. The prime is selected from N, not from the input values.

Let d=d_N=floor(kr/t)+kt+r(r+1)+1. Since kr<=m+k=N and the estimates
above hold, we have

    floor(kr/t) <= (25/12)*z^4,
    kt <= (845/256)*z^4,
    r(r+1) <= (16770/65536)*z^4,
    d+1 <= [25/12+845/256+16770/65536+2/65536]*z^4.

The last line uses z^4>=16^4=65536 to bound the two constant units.
For completeness the bracket is less than 6. Hence d+1<6*z^4 while
m>(31/32)*z^5 >= (31/2)*z^4, proving 1<=d<m at these lengths.
The affine-disperser hypothesis is therefore nonvacuous throughout.

Part C applied to the m-bit valid-prefix restriction gives exactly

    C_B2(F_N) >= 4N - R(N),     R(N)=4k+9(d+1).

Using 4k <= (65/4)*z^3 <= (65/64)*z^4 gives

    R(N)/z^4 <= 65/64 + 9*(25/12+845/256+16770/65536+2/65536)
               = 848297/16384 < 52.

Therefore the main finite theorem is

    C_B2(F_N) >= 4N - 52*N^(4/5)    for every integer N>=2^20.

At N>=2^20 this right side is at least (3/4)*N, since N^(1/5)>=16.
For every real epsilon>0, take

    N_epsilon=max(2^20, ceil((52/epsilon)^5)).

Then 52*N^(-1/5)<=epsilon for all integers N>=N_epsilon, and hence
C_B2(F_N)>=(4-epsilon)*N at all these lengths. In relative form, for
eta>0 choose N_eta=max(2^20,ceil((13/eta)^5)) to obtain
C_B2(F_N)>=(1-eta)*4N. Thus the full guarantee g(N)=4N and the finite
projection (4/1,1/1) describe exactly the same family and theorem.

The class lemma in Part C is an intermediate lemma, not an independent
score coordinate. Part B is a lemma about the payload restriction, not
about the whole prefix language. The joint claims require all four parts.
No superlinear or endpoint theorem, optimality claim, formal verification,
or inherited acceptance is asserted.
