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
strictly greater than max(A,B), where

    A=(t-1)k+(r-t+1)(r+1),   B=kr/t+r(r+1).

The construction and coefficient-isolation route are the binary case of
BK's Theorem 7.3. The threshold above sharpens the block selection in the
supplied predecessor proof; no priority over the original literature is
asserted. We prove the strengthened statement below.

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

Assume D>max(A,B), with A and B as in the lemma. Put
S={i:d_i>r+1} and s=|S|. The guard k*t>=r(r+1) and t<r
imply k>r+1. If s<t, then

    D <= s*k+(r-s)(r+1)
      <= (t-1)k+(r-t+1)(r+1)=A,

since the middle bound increases with s. This contradicts D>A, so s>=t.
The low-dimension blocks contribute at most (r-s)(r+1). Thus D>B gives

    sum_(i in S)d_i >= D-(r-s)(r+1)
                     > kr/t+s(r+1).

Choose t indices I0 in S with largest d_i, breaking ties by index. Their
sum is at least t/s times the sum over S: for a sorted list, the mean
of its largest t members is at least its full mean. Consequently

    sum_(i in I0)d_i > kr/s+t(r+1),
    sum_(i in I0)(d_i-(r+1)) > kr/s >= k,                    (4)

because s<=r. Strictness comes from D>B even when s=r. Averaging over
the actual s, rather than replacing it by r before subtracting the low
blocks, is the sharpened step.

For each selected i, e_i=r+1-i is at least 1, and e_i+1=w_i<=r+1<d_i.
Thus e_i<d_i-1 and d_i>1. Equation (4) implies

    sum_(i in I0)(d_i-w_i)>k>=k-max_(i in I0)d_i+1.

The independent-subspace assertion therefore applies with every hypothesis.
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

For the fixed family in proof.md, k is a prime supplied by the evaluator's length rule,
K=F_2[Z]/(Q) for an accepted prefix, and pi is its constant-coordinate
functional. The guard supplies r>t>=1 and k*t>=r(r+1). The zero-padded
payload embedding is injective, so every payload affine space of dimension
d_select=1+max((t-1)k+(r-t+1)(r+1),floor(kr/t)+r(r+1)) maps to an affine space satisfying the lemma.
No field-encoding loss besides the counted prefix k is suppressed.

The source route is BK Lemmas 3.7,3.10,3.13; Theorem 6.1; Lemma 7.2;
Theorem 7.3 and Claim 7.4. We specialized the zero-pattern argument to
prime k and used (6) to give a direct proof of the final combinatorial
uniqueness. Only elementary finite-field algebra, root counting, linear
algebra and binary polynomial expansion remain outside this file. The
finite-field existence and polynomial-time evaluation arguments are
supplied in proof.md. No efficient search for a hard truth table or for
an irreducible prefix is assumed.

## 9. Comparison with the supplied predecessor threshold

Write d_old=floor(kr/t)+kt+r(r+1)+1 and

d_new=1+max((t-1)k+(r-t+1)(r+1),floor(kr/t)+r(r+1)).

Both are strict-threshold integer dimensions. Direct subtraction gives

 d_old-d_new=min(kt, floor(kr/t)+k+(t-1)(r+1))>0.

Thus the same field, symmetric polynomial, payload injection and input
family receive a strictly smaller sufficient disperser dimension at every
accepted parameter tuple. No evaluator changes or comparisons between
different functions are involved. Under the chosen all-length parameters,
both dimensions remain O(N^(4/5)); their leading asymptotic grammar
expression is unchanged. Only the finite threshold and resulting gate
remainder are strengthened. This proof uses exact inequalities and does
not rely on a finite dimension sweep or a source judgment.

## 10. Use each block's individual Frobenius weight

There is a stronger threshold for the very same symmetric polynomial.
Set W=sum_(i=1)^r w_i=r(r+3)/2 and H=(t-1)(t+2)/2. We prove
nonconstancy on every affine subspace of dimension strictly greater than

    W+max(kr/t, (t-1)k-H).                                  (7)

For the triangular decomposition (3), put delta_i=max(d_i-w_i,0).
Since d_i<=w_i+delta_i, their total T=sum_i delta_i satisfies T>=D-W.
Under (7), T>kr/t and T>(t-1)k-H. The guard implies k>r+1, as
proved in Section 5. Thus all capacities k-w_i are positive. If at
most t-1 excesses were positive, their sum would be at most the sum
of the t-1 largest capacities k-w_i. The weights are precisely the
integers2,3,...,r+1, so that sum is

    (t-1)k-(2+3+...+t)=(t-1)k-H.

This contradicts the strict lower bound on T. Consequently at least t
excesses are positive. Choose t indices with the largest delta_i,
breaking ties by index. The top-t mean is at least the full r-term
mean, giving

    sum_selected delta_i >= (t/r)T > k.

For every selected index, d_i>w_i=e_i+1>=2, with e_i=r+1-i>=1.
Thus all individual hypotheses of Section 3 hold, and

    sum_selected(d_i-w_i)>k>=k-max_selected d_i+1.

Apply that section to obtain a nonzero full-weight coefficient. Sections
6 and 7 require only the selected individual weights and the triangular
index structure, not a common d_i>r+1. Their descending-factor argument
and weighted-capacity isolation therefore apply unchanged. This proves
(7) for every field and tuple satisfying the original lemma hypotheses.

The sufficient integer dimension is

 d_weight=1+W+max(floor(kr/t),(t-1)k-H).

It strictly improves Section 5's actual-s dimension d_select. To see this,
let A=(t-1)k+(r-t+1)(r+1) and B0=floor(kr/t)+r(r+1). Exact identities are

 W+(t-1)k-H = A-(r-t)(r-t+1)/2,
 W+floor(kr/t) = B0-r(r-1)/2.

Because 1<=t<r, both decreases are positive integers, and the first is
no larger than the second. Subtracting these respective losses before
taking a maximum gives

 d_select-d_weight >= (r-t)(r-t+1)/2 >= 1.                   (8)

This is an all-parameter theorem, not an inference from samples. Applied
to the existing fixed payload injection, d_weight is a sufficient
disperser dimension; all field-prefix, circuit and all-length arguments
remain unchanged. The earlier thresholds remain supplied supporting
lemmas. No literature priority is asserted for this sharpening.

## 11. Boundary-dimension blocks and the largest positive image

We can sharpen the same construction further by allowing selected d_i=w_i.
First record the needed coefficient argument. For a dimension-d subspace V,
with 2<=d<k and image polynomial Q, the first d Frobenius coefficient rows
of Q(X)^(2^u),u=0,...,d-1, have rank d over K. If a nonzero linear
combination of these rows vanished, the nonzero linearized polynomial
L(Y)=sum_(u=0)^(d-1)c_u Y^(2^u), of degree at most2^(d-1), would vanish
on Q(K)=V, a set of2^d points. The root bound forbids this. Hence some
d-column minor is nonzero. Its determinant equals the coefficient of the
corresponding full-weight monomial in Q(X)^(2^d-1), by the same no-carry
expansion as Section3.

Suppose now that selected blocks have d_i>=w_i, at least one has positive
excess delta_i=d_i-w_i, and the positive excesses have sum S satisfying

    S>k-max_(positive i)d_i+1.                               (9)

Apply Section3's explicit rho-isolation construction just to the positive
blocks. It selects minors nonzero at rho=0 and makes a factor zero for
EVERY rho!=0. For a zero-excess block d_i=w_i, select any nonzero minor
from the preceding rank argument. Multiplying these factors into equation
(1) preserves all zero terms and leaves the rho=0 term nonzero. Thus the
same coefficient conclusion holds with these additional boundary blocks.
Sections6–7 then isolate it in the full symmetric polynomial exactly as
before. No stronger common dimension hypothesis is needed.

Set W=r(r+3)/2,H=(t-1)(t+2)/2, and

 A0=(t-1)k-H-(r-t+1),  B0=r(k-1)/(t+1).

We claim dispersion above W+max(A0,B0). For arbitrary dimensions d_i,
put delta_i=max(d_i-w_i,0), T=sum_i delta_i. Then T>=D-W>B0>0.
There must be at least t blocks with d_i>=w_i. Otherwise, if their
number is s<=t-1, every other block has d_i<=w_i-1 by integrality.
Using k>r+1, the capacities k-w_i+1 are all positive, so

 D<=W-r+sum_(s blocks)(k-w_i+1)
  <=W-r+(t-1)(k+1)-H=W+A0,

contradicting the assumed dimension. The middle maximum takes the t-1
smallest weights2,...,t; the claim also covers t=1 with an empty sum.

Choose t blocks from {i:d_i>=w_i} with the largest delta_i, breaking
ties by index. Their positive excesses sum to S>=(t/r)T, so

    S>t(k-1)/(t+1)>0.

Among these positive-excess blocks some has delta_i>=S/t, since their
number is at most t. Its dimension is at least S/t+2. Writing d_max for
the largest POSITIVE-excess image dimension therefore gives

 k-d_max+1<=k-S/t-1<S.

The last inequality is exactly S>t(k-1)/(t+1). This proves (9), using
only a positive-excess block for d_max. The extended coefficient argument
above applies to all selected blocks, proving the claimed dispersion.

The sufficient integer dimension is

 d_boundary=1+W+max(A0,floor(r(k-1)/(t+1))).

Let d_weight be Section10's threshold. Then

 d_weight-d_boundary >= min(r-t+1,
             floor(kr/t)-floor(r(k-1)/(t+1))) >=2.           (10)

For the last inequality, r-t+1>=2. The difference of the two unfloored
terms is r(k+t)/(t(t+1)); since r>=t+1 and k>=r+2, it is at least
(2t+3)/t>2. Taking floors decreases this difference by less than one,
so its integer difference is at least2. Subtracting the two positive
branch losses before taking a maximum proves the first inequality.
All claims hold for every guarded tuple, not just asymptotic choices.

## 12. Distinct positive weights sharpen the remaining branch

Assume in addition k>=t(t+1)/2. This holds at every sufficiently large
length of the fixed family because k=Theta(N^(3/5)) and t^2=Theta(N^(2/5)).
Keep W,H,A0 and the eligible-block argument of Section11, but replace
its rational B0 by

    Bstar=r*k/(t+1)-r/2.

Suppose D>W+max(A0,Bstar). At least t blocks have d_i>=w_i, by the
unchanged capacity proof. Select the t largest eligible excesses and
let S be their positive sum. Then

 S>=(t/r)T>t*k/(t+1)-t/2>=t(t-1)/2.

Let u,1<=u<=t, count their positive excesses. The corresponding weights
are distinct integers at least2, so their sum is at least u(u+3)/2.
Averaging the POSITIVE image dimensions gives

 d_max>=S/u+(u+3)/2.

For 1<=j<t, the difference between adjacent values of the right-hand
expression is S/(j(j+1))-1/2>=0. Thus its minimum over u<=t occurs at t,
and d_max>=S/t+(t+3)/2. Consequently

 k-d_max+1<=k-S/t-(t+1)/2<S,

where the last inequality is exactly S>t*k/(t+1)-t/2. This is the
positive-group isolation premise (9). Attach boundary blocks as in
Section11, and apply Sections6–7. This proves dispersion above the
new threshold without any finite computation.

The sufficient integer dimension is

 d_sharp=1+W+max(A0,floor(r*(2k-t-1)/(2*(t+1)))).

Its second branch is no larger than Section11's because their rational
difference is r(t-1)/(2(t+1))>=0. Hence d_sharp<=d_boundary. The first
branch is unchanged, so no strict improvement is asserted at every length.
For r=t^2 with t>=3, however, the new rational second branch exceeds A0 by

    k/(t+1)+t^2-t/2>1.

Both maxima therefore use their second branch, and taking floors yields

 d_boundary-d_sharp>=floor(t^2*(t-1)/(2*(t+1)))>=2.           (11)

There are infinitely many such lengths in this ONE fixed family. Extend
r(N)=ceil(N/k(N))-1 to all N>=1. The selected k(N) is nondecreasing, so
r(N+1)<=r(N)+1. Also r(N)->infinity by the proved Bertrand bounds.
For every sufficiently large target R, the first N with r(N)>=R has
r(N)=R: an upward step cannot skip an integer. Thus every sufficiently
large square R=t^2 is attained. These lengths tend to infinity, so all
original guards and the new eventual k condition hold there. This proves
strict improvements on infinitely many actual integer lengths, while
the main circuit theorem still holds at EVERY sufficiently large length.
No sparse-length argument is used to replace the main liminf quantifier.


## 13. Exact integer balancing of the selected excesses

Under the same additional k>=t(t+1)/2, set

 L=floor(t*k/(t+1)-t/2), q=floor(L/t), a=L-t*q,
 B_int=r*q+a, d_int=1+W+max(A0,B_int).                       (12)

Here q is a scalar quotient, unrelated to gate-proof pair counts; 0<=a<t.
We prove dispersion at dimension d_int for the identical field function.
The eligible-block argument of Section11 is unchanged: D>W+A0 ensures
at least t blocks with d_i>=w_i. Order all r nonnegative integer excesses
in decreasing order, and select the first t eligible blocks. They have
the same excess sum S as the first t of all r: any ineligible excess is
zero and eligible zero-excess blocks fill the selection if necessary.

If S<=L, the t-th largest excess is at most q=floor(L/t). Every
remaining excess is at most q. Thus the full sum T satisfies

 T<=S+(r-t)q<=L+(r-t)q=r*q+a=B_int.

But T>=D-W>B_int, a contradiction. Hence S>L>=0 and there are u
positive selected blocks, 1<=u<=t. Let d_max be their largest dimension.
If the coefficient-isolation condition (9) failed, then S+d_max<=k+1.
Their distinct weights sum to at least u(u+3)/2, so

 S+u(u+3)/2 <= sum_positive d_i <=u*d_max<=u(k+1-S).

Consequently S<=u*k/(u+1)-u/2. The function of the integer u on the
right is nondecreasing for 1<=u<=t: its next difference is
k/((u+1)(u+2))-1/2>=0 by k>=t(t+1)/2. Therefore
S<=t*k/(t+1)-t/2, and integrality forces S<=L, again a contradiction.
Condition (9) holds. Append selected boundary blocks as in Section11
and apply Sections6–7 to obtain nonconstancy. This proves the new
threshold using only integer inequalities, without a new experiment.

The old second branch equals floor(r*(L+theta)/t) for some 0<=theta<1.
Since r>t and a>=0,

 B_int=r*q+a <=r*q+floor(r*a/t)=floor(r*L/t)
             <=floor(r*(2k-t-1)/(2*(t+1))).                  (13)

Thus d_int<=d_sharp at every tuple in scope, and all previous thresholds
remain valid supporting lemmas. The integer balancing bound itself is
exact for arbitrary nonnegative integer excesses: a entries q+1 and the
remaining r-a entries q have top-t sum L and total B_int. This exactness
is only for that scalar relaxation; no optimal disperser threshold or
actual optimal circuit complexity is claimed.

The old strict square-block comparison also remains true with d_int in
place of d_sharp, by transitivity. No new infinite strict-improvement
claim about d_int versus d_sharp is needed. For an explicit guarded
old32-rule example N=565701,k=5657,r=100,t=10 one has L=5137,q=513,a=7,
B_int=51307 and d_int=56458, versus d_sharp=56528. This improves the
same-function circuit certificate by560 gates at the same additive term.
At the old32-rule tuple N=1048576 the first branch dominates and d_int=d_sharp=98362.
Finally B_int differs from the old second branch by at most O(r), so
it does not alter d=k*sqrt(r)+r^2/2+O(k+r) or any all-length asymptotic
comparison in proof.md. These statements concern proved certificates.


## 14. Joint choice of excess sum and largest image dimension

Keep the hypotheses of Section13. The selected t blocks need not be
those with the largest excesses: we can optimize their sum and their
largest positive image dimension together. Define

 B_joint=r*(2k-r-1)/(2*(t+1)),
 d_joint=1+W+max(A0,min(B_int,floor(B_joint))).                (14)

We prove dispersion above W+max(A0,B_joint). Combining this with
Section13 proves the integer dimension (14), which is never larger
than d_int. The min combines two separately proved sufficient thresholds;
there is no assertion that both inequalities hold at the smaller one.

Let v>=1 be the total number of positive excesses, with total T.
As before D>W+A0 supplies t eligible blocks; D>W+B_joint implies
T>B_joint>0. Suppose no eligible t-subset satisfies isolation (9).
If v>=t, take any t-subset I of the positive indices and put
S_I=sum_(i in I)delta_i. Failure of (9) implies

 S_I+max_(i in I)d_i<=k+1,
 (t+1)S_I+sum_(i in I)w_i<=t(k+1).

Average this last inequality over ALL t-subsets of the v positive
indices. Each index occurs in the same fraction t/v of the subsets.
Multiplying the resulting inequality by v/t yields

 (t+1)T+sum_positive w_i<=v(k+1).

Since these v weights are distinct integers at least2,

 T<=v*(k-(v+1)/2)/(t+1) <=r*(k-(r+1)/2)/(t+1)=B_joint.

The second inequality holds because h(j)=j*(k-(j+1)/2) has
h(j+1)-h(j)=k-j-1>0 for j<r, by k>r+1. This contradicts T>B_joint.
If v<t, select all v positive blocks and fill with eligible boundary
blocks. Failure then gives

 T<=v*k/(v+1)-v/2<=t*k/(t+1)-t/2
   =h(t)/(t+1)<=h(r)/(t+1)=B_joint.

The first comparison uses the Section13 monotonicity under
k>=t(t+1)/2; the last uses the same increasing h. Again contradiction.
Thus some eligible t-subset satisfies (9). Its positive blocks provide
rho-isolation, its boundary blocks nonzero minors, and Sections6–7
isolate the coefficient in the complete symmetric polynomial. This
proves the stated threshold without an experiment or a new family.

For this fixed all-length family the second branch is eventually
unnecessary. Indeed r<=t^2 and k=Theta(r^(3/2)),t=Theta(sqrt(r)), so

 B_joint-A0 <= k/(t+1)-r*(r+1)/(2*(t+1))+H+r-t+1.

The positive terms on the right are O(r), while the negative term is
-Theta(r^(3/2)). Therefore B_joint<A0 at every sufficiently large
integer length and

 d_joint=1+W+A0
        =(t-1)k+(r*(r+1)-t*(t-1))/2+1.                    (15)

All these are inequalities of proved thresholds, not complexity ratios.
There is strict improvement over d_int at infinitely many ACTUAL lengths.
At r=t^2 let theta=t*k/(t+1)-t/2-L in [0,1), and L=tq+a.
The preceding real branch Bstar satisfies

 Bstar-B_int=(r/t-1)*a+r*theta/t <r-t+1<r,
 Bstar-A0=k/(t+1)+t^2-t/2.

Hence B_int-A0>k/(t+1)-t/2>=0 under the existing k hypothesis.
Since both are integers, once (15) holds, d_int-d_joint>=1, and in
fact it is at least floor(k/(t+1)-t/2)+1. Section12 proves that these
square r occur at infinitely many actual lengths; all eventual conditions
hold there. The main theorem continues to cover ALL sufficiently large
integer lengths, including those at which no strict improvement occurs.

For the old32-rule guarded tuple N565701,k5657,r100,t10, floor(B_joint)=50968
and A0=50768, so d_joint=56119, versus d_int56458 and d_sharp56528.
The additive-seven same-family circuit expression rises from1794176
to1796888, a2712-gate improvement over integer balancing. At the old32-rule tuple N1048576
the first branch still dominates and d_joint=98362. The reduction to
(15) is eventual, not asserted at the first example. Asymptotically
(15) retains k*sqrt(r)+r^2/2+O(k+r), so every previously proved
all-length leading expression and certificate-ratio limit remains valid.
