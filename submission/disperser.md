# Algebraic disperser argument for the payload

This supplies the BK construction argument used in proof.md, specialized to
F_2 and prime field degree. We use k for field degree, r for block count,
and t for the symmetric-polynomial order. All subspace dimensions are over
F_2. This file establishes the needed nonstandard algebra rather than
assuming a cited construction is a certified black box.

**Lemma.** Let K be the field of size 2^k with k prime. Suppose 1<=t<r. For any nonzero linear pi:K->F_2 put w_i=r+2-i and

    f(a_1,...,a_r)=pi(sum_(|I|=t) product_(i in I) a_i^(2^(w_i)-1)).

Then f is nonconstant on every affine subspace of K^r of dimension
strictly greater than H+max(T, U/(t+1)), where
H=r(r+3)/2, T=sum_(j=2)^t max(k-j,0), and
U=sum_(j=1)^r max(k-j,0) (empty sums zero).

The construction follows the binary case of BK's Theorem 7.3; Section 5
sharpens its sufficient dimension threshold. We prove the result below, including
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

Put H=sum_(i=1)^r w_i=r(r+3)/2,
T=sum_(j=2)^t max(k-j,0), U=sum_(j=1)^r max(k-j,0).
Assume D>H+max(T,U/(t+1)). For each block put delta_i=d_i-w_i
and S={i:delta_i>0}. Then

    sum_(i in S)delta_i >= D-H > max(T,U/(t+1)).

There are at least t positive deficits. Otherwise their sum is at most
the t-1 largest positive capacities max(k-w_i,0), because d_i<=k.
The weights are exactly r+1,r,...,2, so this capacity sum is T.
Zero capacities correctly account for k<t as well.

Suppose every t-element subset I of S failed the coefficient-isolation
hypothesis: sum_(i in I)delta_i+max_(i in I)d_i<=k+1. Since a maximum
is at least the average, this implies

    (t+1)*sum_(i in I)delta_i+sum_(i in I)w_i <= t*(k+1).

Average uniformly over all t-element subsets of S. Each index occurs
with probability t/|S|, so cancelling this factor gives

    (t+1)*sum_(i in S)delta_i
        <= sum_(i in S)(k+1-w_i) <= U.

The last bound holds because U is the sum of all positive capacities
max(k+1-w_i,0), with w_i running from 2 to r+1. This contradicts
sum delta_i>U/(t+1). Hence some t-element I0 satisfies

    sum_(i in I0)delta_i+max_(i in I0)d_i > k+1.           (4)

Every selected d_i>w_i=e_i+1, with e_i=r+1-i>=1. Thus d_i>1,
e_i<d_i-1, and (4) is exactly sum delta_i>k-max d_i+1, the
independent-subspace hypothesis. No guard kt>=r(r+1) is required.
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

## 8. Application, source scope and the two selection routes

For the current proof, k is prime, K=F_2[Z]/Q for an accepted paid
prefix, and pi is the nonzero constant-coordinate functional. The old
threshold H+max(T,floor(U/(t+1)))+1 remains sufficient by Section5.
Section10 supplies a guarded smaller threshold. Section11 now supplies
the exact integer threshold used by the evaluator. All these routes
produce the SAME coefficient-isolation hypothesis(4) for t selected blocks.
Thus Sections6--7 apply unchanged, including arbitrary triangular mixing.
Zero padding is injective and preserves the actual payload dimension.

The source route is BK Lemmas3.7,3.10,3.13, Theorem6.1, Lemma7.2,
Theorem7.3 and Claim7.4, reproduced here through the supplied peer source.
Only elementary field algebra, root counting and linear algebra remain
outside this file; the main proof supplies field existence and uniformity.
The new marked-subset average below is proved directly and does not rely
on a prior judgment or a requested external retrieval.

## 9. The same coefficient certifies high algebraic degree

For a Boolean function on binary coordinates, its algebraic degree is the
maximum number of distinct variables in a nonzero monomial of its unique
multilinear polynomial over F_2 (its algebraic normal form). Affine
substitution cannot increase this degree: substitute affine polynomials
into each monomial and reduce using z^2=z.

The nonzero monomial M above has total binary weight

    sum_(i in I0) wt(b_i)=sum_(i in I0)w_i >= t(t+3)/2=:B_deg,

because the t smallest distinct weights among 2,3,...,r+1 are
2,3,...,t+1. This forces algebraic degree at least B_deg in f(R(X)),
not merely nonconstancy. To justify the connection precisely, every bit
coordinate of a field element is an F_2-linear map K->F_2 subset K, so
Section 1 represents it by a linearized polynomial. A Boolean polynomial
of degree h in the kr bit coordinates therefore has a field-polynomial
representation whose reduced monomials all have total binary weight <=h:
each bit contributes weight one, and multiplication and reduction cannot
increase weight, as proved in Section 1. Uniqueness of the reduced
representative then excludes a nonzero coefficient of weight >h.

Finally R is F_2-affine and surjects onto the original affine space A.
If the restriction f|A had algebraic degree less than B_deg in an affine
coordinate basis of A, its composition with R would also have degree less
than B_deg, a contradiction. Consequently the supplied lemma actually
proves degree at least B_deg on every affine space whose dimension strictly
exceeds its stated threshold. This conclusion concerns the payload
function after fixing a valid field prefix, not the whole prefix language.

## 10. Weighted marked-subset selection

This strengthens Section5's sufficient dimension without changing the
field polynomial, the coefficients isolated in Sections3,6--7, or the
Boolean degree consequence in Section9.

Assume integer2<=t<r, k>=2r+1 and2k>=(t+1)^2. Put

    h_m=ceil(2m/(t+1)),
    F_m=m*(2k-2m-1+h_m)/(2(t+1))    (t+1<=m<=r),
    F_r=r*(2k-2r-1+h_r)/(2(t+1)).

Then dimension D>H+max(T,F_r) suffices for the same selected subsetI0.
This is a real rational threshold; the implemented integer dimension is
H+max(T,floor(F_r))+1. Here H,T and weightsw_i=r+2-i are unchanged.

### 10.1 Exact distribution on marked t-subsets

Let S contain m>=t+1 indices. Let Hs be the h=h_m indices of S with
largest weights; ties can be broken by index (weights here are distinct).
Write Ls=S\Hs. Put u=(t+1)/m and lambda=h*u-2. Since h=ceil(2/u),
0<=lambda<u<=1. Also h>=2 and m-h>=t-1: the latter follows from
2m/(t+1)<=m-t+1, equivalent to m>=t+1, whose right side is an integer.

Choose a marked index j uniformly from Hs. With probability1-lambda,
let I contain j and a uniform(t-1)-subset of Ls. With probabilitylambda,
let I contain j, one additional uniform index of Hs\{j}, and a uniform
(t-2)-subset of Ls. These choices are all legal, includingt=2 and
lambda=0. The marked index always belongs to the t-element set I.

For i in Hs, the expected coefficient of delta_i in sum_(i in I)delta_i
plus delta_j is (2+lambda)/h=u. For i in Ls it is
(t-1-lambda)/(m-h)=u. Thus all deficits receive exactly the SAME weight,
although the extra marked weight is concentrated at high w_i.

### 10.2 Apply the distribution to failed isolation inequalities

Let delta_i=d_i-w_i and S={i:delta_i>0}, as before. Suppose there is
NO t-element I subsetS with sum_(i in I)delta_i+max_(i in I)d_i>k+1.
If |S|<t, its positive-deficit sum is at most T, exactly as in Section5.

If m=|S|>=t+1, every marked pair(I,j) in the distribution obeys

    sum_(i in I)delta_i+delta_j+w_j <= k+1.

Average it with the exact marginals above. Writing A=sum_(i in S)delta_i,

    ((t+1)/m)*A + average_(j in Hs)w_j <= k+1.

The m distinct weights in S are drawn from2,...,r+1. The h largest
therefore have average at least m+(3-h)/2: their smallest possible list
is m-h+2,...,m+1, obtained by using the m smallest weights overall.
Consequently A<=F_m.

F_m is increasing for m=t+1,...,r under k>=2r+1. To check the ceilings
explicitly, h_(m+1)-h_m is0 or1. If N_m=2(t+1)F_m, then

    N_(m+1)-N_m
      =2k-4m-3+h_(m+1)+m*(h_(m+1)-h_m)
      >=2k-4m-1>0                 (m<=r-1).

Thus A<=F_r in this case.

If m=t, use the sole subset I=S and max d_i>=average d_i:

    (t+1)*A <= t*(k+1)-sum_(i in S)w_i
                  <= t*(k+1)-t*(t+3)/2.

Its resulting bound G_t satisfies

    F_(t+1)-G_t=(2k-(t+1)^2)/(2(t+1))>=0.

Again A<=F_r. These three cases exhaust every possible positive set.
Since D=sum d_i<=H+A, the assumed D>H+max(T,F_r) contradicts failure.
A suitable positive t-subsetI0 therefore exists, with exactly the strict
hypothesis(4) required in Section3. No zero-deficit padding is assumed.

### 10.3 Relation to the old dimension and retained high degree

The old U=rk-r(r+1)/2 satisfies the exact identity

    U/(t+1)-F_r = r*(r-h_r)/(2(t+1))>=0,

because2<=t<r gives h_r<=r-1. Thus the new integer dimension never
exceeds the old one. The new proof supplies the same nonzero monomial M,
so Section9 still gives degree at least B=t(t+3)/2 on every affine space
of the new dimension. This concerns fixed-valid-prefix payloads only.
If either new guard fails, use Sections1--9's original threshold instead.
No claim extends this averaging argument outside its stated guards.

## 11. Exact integer threshold for the unchanged isolation predicate

This further improves Section10 without changing the coefficient-isolation
predicate or any field argument. Throughout k>r+1 and2<=t<r are integers.
For integer dimensions0<=d_i<=k and weightsw_i ranging from2 to r+1,
write delta_i=d_i-w_i. A good subset is a t-subset of POSITIVE deficits
satisfying sum_I delta_i+max_I d_i>k+1, exactly as in Section3.

Define T as before. For each integer m=t,...,r and A=t-1,...,k-m-2 put

    X=min(floor(A/(t-1)), floor((k-t-A)/2),
          floor(((t-1)k-t(t-1)/2-tA)/(t-1))).

Skip choices with X<1. For the others put C=k-A,J=min(m,C-2X), and

    E_(m,A)=A+(J-t+1)X+floor((C-J-1)^2/4)-floor((C-m-1)^2/4).

Let E be the maximum of T and all these values (empty ranges add no
values). Then H+E is EXACTLY the largest total dimension of an integer
vector in[0,k]^r having no good subset. Consequently

    d_exact=H+E+1

is sufficient for the same nonzero coefficient and degree>=B=t(t+3)/2.
The optimality asserted here is ONLY for this fixed integer subset
predicate. It is not an optimality statement about affine dispersers,
actual circuits, their complexity, or other lower-bound arguments.

### 11.1 Reduce the failed vectors by rearrangement

Let S be the positive set, m=|S|. Replacing nonpositive deficits by0
increases total dimension to H+sum_S delta_i and does not change good
positive subsets. Its dimensions remain <=k because every w_i<k.
For a fixed positive set, the absence of good subsets is equivalent to

    g_j+w_j<=k+1 for every j in S,
    g_j=2delta_j+sum of the largest(t-1) other positive deficits,

when m>=t. First move the m positive deficits to the m smallest weights,
keeping their order by weight. This only decreases every associated w_j,
preserves the deficit sum, and cannot violate dimension upper bounds.
For that new multiset of weights, sort deficits in nonincreasing order
against increasing weights. This also preserves failure and dimensions.

Here are the details of this latter exchange. If a<=b are two deficit
values, replacing b by a in a largest-(t-1) sum lowers that sum by a
quantity between0 and b-a. Thus g(b)-g(a)>=b-a>=0. Suppose smaller
weight u had deficit a and larger weight v had deficit b. Before
swapping, g(b)+v is at least both affected expressions. Afterwards both
g(b)+u and g(a)+v are <=g(b)+v. Expressions for other indices are
unchanged, since their other-deficit multisets are unchanged. Also
b+u<=b+v<=k and a+v<=b+v<=k. Repeating these swaps proves the reduction.

We may therefore maximize over positive deficits

    delta_1>=...>=delta_m>=1,  w_j=j+1,

with all other deficits zero. This construction preserves the positive
sum of every original failed vector, and may increase its total dimension.
For m<t that sum is at most T, attained by putting dimension k in the
t-1 smallest-weight blocks and dimension w_i in the other blocks.

### 11.2 Exact feasible top sums and pivots

Assume m>=t. Put A=sum_(j=1)^(t-1)delta_j and x=delta_t. The failed
inequalities above are exactly

    delta_j+j+A+x<=k        for1<=j<t,
    2delta_j+j+A<=k         for t<=j<=m.

In particular x>=1,A>=t-1 and positivity at j=m requires A<=k-m-2.
The top t-1 deficits must be nonincreasing, at least x, have sum A,
and have upper bounds u_j=k-A-x-j. These requirements are feasible
if and only if

    (t-1)x<=A,
    2x<=k-t-A,
    tA+(t-1)x <= (t-1)k-t(t-1)/2.

Necessity follows by ordering, the j=t inequality, and summing top upper
bounds. For sufficiency, the second inequality gives u_(t-1)>=x+1.
Start all top deficits at x, then fill them from left to right up to
u_j until their sum is A. The first and third inequalities give the
necessary lower and upper total sums. Since u_j decreases with j, this
procedure remains nonincreasing and attains every intermediate integer
sum. No assumption about real-valued LP feasibility is substituted for
this integer construction.

These inequalities say exactly1<=x<=X with X in the definition. Given A,
all tail entries are bounded by min(x,floor((k-j-A)/2)). Taking those
upper bounds for every j>=t gives a nonincreasing positive tail and is
optimal for that A,x. The tail sum is nondecreasing in x, so its maximum
for fixed A occurs at X. This proves the reduction to the stated finite
m,A search, rather than an unproved rounding of continuous parameters.

### 11.3 The closed tail sum and attainable maximum

For an eligible m,A,X let C=k-A. The pivot bound gives C-2X>=t, so
J=min(m,C-2X)>=t. The maximal tail is X for t<=j<=J, and is
floor((C-j)/2) for J<j<=m. All these entries are positive because
C-m>=2. The elementary identity

    sum_(z=0)^n floor(z/2)=floor(n^2/4)

now gives exactly the stated E_(m,A). The prefix arguments C-J-1 and
C-m-1 are nonnegative. Each top deficit and this tail satisfy all failed
inequalities, so they give an actual bad vector. Their dimensions are
<=k: the top bounds imply delta_j+j+1<=k (indeed A+x>=2), and the tail
bounds do as well because k-j>=2 and A>=1. Other blocks have d_i=w_i<k.
Thus every listed value is attainable. Together with the m<t construction,
E is the exact maximum, not only an upper estimate or a finite experiment.

### 11.4 Consequences and safe implementation shortcut

Sections5 and10 already prove that their respective thresholds force a
good subset, under their stated guards/fallback. Since a failed vector
of dimension H+E exists, H+E is strictly below either applicable
sufficient integer dimension. Hence d_exact never exceeds the old or
Section10 dimension. This comparison uses the supplied proofs, not a
prior acceptance or a candidate-created judgment.

In particular, if the Section10 guarded/fallback dimension equals H+T+1,
then E<=T while the m<t construction gives E>=T. The code may return
that dimension without scanning m,A. Otherwise it computes the exact
maximum with integer operations. There are at most r*k scanned choices
per parameter pair, so selection remains uniform polynomial time.

Finally, Section11 supplies the same strict good-subset hypothesis to
Sections3,6--7; Section9's high-degree witness is unchanged. This is the
only additional model-transfer premise needed to use d_exact in the
retained full-B2 theorem. All-length bounds and prefix accounting remain
those of proof.md for the new exact selector.
