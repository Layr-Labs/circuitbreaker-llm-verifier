# A uniform certificate-prefix family: coefficient 4, exponent 1

## Claim and provenance

For the single total Boolean language defined below, with **N denoting the
entire input length**, we prove

    C_B2(F_N) >= 4N - R(N),       R(N)/N -> 0.

The joint pair is (4,1). The leading coefficient is inherited mathematics.
refinement.md supplies the sharp finite bound4m-8d+8
for d>=3, without increasing that leading coefficient. The leading gate argument is the
parameter method of Mishra and Santhanam, ECCC TR26-260, revision 1, November 3, 2026,
Sections 2–4. The algebraic construction is based on Ben-Sasson
and Kopparty (BK), *Affine Dispersers from Subspace Polynomials*, author
version March 12, 2010, Sections 3, 6 and 7. The necessary arguments are
supplied in our notation in `elimination.md` and Appendix A of `proof.md`.

The paid input prefix supplies the field polynomial without an efficient
construction premise. Only a fixed-prefix restriction is a disperser.
All claims use the full B2 model and actual input length N.

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

## 4. All sufficiently large integer lengths pass the guards

We have k=Theta(N^(3/5)), m=N-k, r=Theta(N^(2/5)), t=Theta(N^(1/5)).
Consequently k<N and t<r eventually. Also r<N/k+1<(1/2)N^(2/5)+1 and

    r(r+1)/t <= sqrt(r)(r+1)
              <= sqrt((1/2)N^(2/5)+1)*((1/2)N^(2/5)+2).

Dividing the last expression by N^(3/5) gives limit 1/(2*sqrt(2)),
whereas k/N^(3/5)>2. Thus k*t>=r(r+1) eventually. These estimates hold
for every sufficiently large integer N, not just selected lengths.


## 5. A stronger certificate for this same fixed language

Fix any irreducible degree-k prefix, whose existence was proved above.
The remaining payload polynomial has nonzero constant-coordinate
functional pi(1)=1. Use the sharper coefficient-selection lemma in
Appendix A and its degree conclusion. For the unchanged evaluator's (k,m,r,t), define

    H=r(r+3)/2,
    T=sum_(j=2)^t max(k-j,0), U=sum_(j=1)^r max(k-j,0),
    D_N=H+max(T,floor(U/(t+1)))+1,
    B_N=t(t+3)/2.

Every dimension-D_N payload affine space has degree at least B_N.
Zero padding is injective and preserves its dimension. Eventually
4<=B_N<=D_N<m, since B_N=O(N^(2/5)) and D_N=O(N^(4/5)). Fixing a
valid paid prefix, applying refinement.md Section7 and recovering k
prefix gates by Section5.1 proves

    C_B2(F_N)>=4N-3k-8D_N+4B_N+3.                       (16)

The evaluator is unchanged. Sharper dimension and degree arguments are
credited to the declared8eef source; its alternative family is unused.

For an explicit same-function comparison, put

    d_old=floor(kr/t)+kt+r(r+1)+1,
    L_old=4k+8d_old-8, L_new=3k+8D_N-4B_N-3.

L_old is exactly the preceding +8 certificate's complete loss, including
the paid prefix. We prove L_new/L_old->1/2 over all actual integer lengths,
without assuming a prime-gap asymptotic. The already proved Bertrand
bounds give k=Theta(N^(3/5)), r=Theta(N^(2/5)), t=ceil(sqrt(r)). Hence
kt and r² are Theta(N^(4/5)), t²-r=O(sqrt(r)), and k>r+1 eventually.
Then T=kt+O(k+t²) and U/(t+1)=kt+O(k+r²/t). Floors add only O(1).
Consequently

    D_N=r²/2+kt+O(N^(3/5)),
    d_old=r²+2kt+O(N^(3/5)), B_N=O(N^(2/5)).

Thus 2L_new-L_old=O(N^(3/5)), while L_old=Theta(N^(4/5)). The ratio
limit proves a strictly better certificate at all sufficiently large
lengths, comparing one family's proved losses, not optimal circuit sizes.
Both losses are o(N), preserving joint leading(4,1).
Relative to the immediate degree-sensitive certificate, the new complete
loss is exactly k smaller, an unbounded Theta(N^(3/5)) improvement.

Put R(N)=L_new eventually, extending it arbitrarily on the finitely
many guarded lengths. It is o(N). For every epsilon>0 take N_epsilon
large enough that R(N)<=4*epsilon*N for every N>=N_epsilon. Then
C_B2(F_N)>=(1-epsilon)*4N at all those actual lengths. This is the
full guarantee and its native projection jointly for the one language.

### 5.1 Recover the paid prefix gates

Let P be the payload restriction at any valid prefix b. Eventually
P has degree at least B_N>=4 and is nonconstant. Every irreducible Q
of degree k>1 has Q(1)=1. Flipping any single coefficient bit toggles
Q(1) to0, making Q divisible by Z+1 and therefore rejected on every
payload. After fixing any subset of prefix coordinates to b, every
remaining prefix bit is still essential: set the other free prefix
bits to b and choose a payload with P=1; flipping that bit changes
output to0. This applies in any fixed order of the k coordinates.

For completeness, if a full-B2 circuit depends on x and its restriction
x=c is not a constant or literal, that restriction has a circuit with
at least one fewer gate. Normalize constant/unary internal gates by
substitution, absorbing complements into descendant tables. If a
complement reaches the output, complement its driving gate's table;
a literal output is the only exception. Choose a gate using x, which
exists by dependence. Fixing x makes this gate constant or unary in
its other pin (also with repeated pins). Remove it by the same rules.
A terminal complement is absorbed into its surviving driving gate;
if no driving gate survived, the restricted function would be a
constant or literal, excluded here. All tables remain among the16
allowed B2 tables, and sharing costs no extra gate: each affected
consumer table is changed in place. Additional simplifications only
reduce size further.

At each successive prefix fixing the remaining function has P as a
further restriction, hence cannot be constant or a literal. Apply this
gate-removal observation k times to any circuit for F_N. It proves

    C_B2(F_N)>=k+C_B2(P).

This is an explicit elimination argument, not direct-sum additivity.
Compared with simple restriction, (16) gains k gates. The prefix is still part of N;
its cost in the complete loss is3k rather than4k.

## 6. Sources and evidence

All material arguments are reproduced in the supplied proofs. Sources:
Mishra--Santhanam, ECCC TR26-260 revision1,2026-11-03, Sections2–4;
BK author preprint2010-03-12, Lemmas3.7/3.10/3.13, Theorem6.1,
Lemma7.2, Theorem7.3/Claim7.4, SHA256
922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf;
Li--Yang, ECCC TR21-023,2021-02-20, Sections1.1/2.1/2.4/3.3 for
the protected model, not their LP. Declared source8eef supplies the
sharper selection, degree bridge and degree bound reproduced here;
its alternative family is not adopted. Prior verdicts are not premises.

# Appendix A. Algebraic disperser argument for the payload

Here k is prime field degree, r the block count, t the symmetric order;
all dimensions are over F_2. The needed BK argument is reproduced.

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

## 8. Exact application and source scope

For proof.md, k is a prime supplied by the evaluator's length rule,
K=F_2[Z]/(Q) for an accepted prefix, and pi is its constant-coordinate
functional. The unchanged evaluator guards supply r>t>=2; no parameter is reselected. The zero-padded
payload embedding is injective, so every payload affine space of dimension
D_N=H+max(T,floor(U/(t+1)))+1 maps to an affine space
satisfying the lemma, with T and H as above.
No field-encoding loss besides the counted prefix k is suppressed.

The source route is BK Lemmas 3.7,3.10,3.13; Theorem 6.1; Lemma 7.2;
Theorem 7.3 and Claim 7.4. We specialized the zero-pattern argument to
prime k and used (6) to give a direct proof of the final combinatorial
uniqueness. Only elementary finite-field algebra, root counting, linear
algebra and binary polynomial expansion remain outside this file. The
finite-field existence and polynomial-time evaluation arguments are
supplied in proof.md. No efficient search for a hard truth table or for
an irreducible prefix is assumed.

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


# Appendix B. Six-variable Ramsey lemma and sharpness at d=3

For d=3 the finite bound is

    C_B2(f) >= 4n-8d+8 = 4n-16.                         (10)

This appendix proves (10) and its sharpness at d=3. The general +8
theorem is proved independently in refinement.md.
The ingredient is the classical six-variable Ramsey statement, with a
complete proof by finite reduction and independently checked refutations
supplied here. No cited classification or previous verdict is a premise.

## B.1 Reduce every six-variable coloring to explicit finite cases

**Lemma S.** Every coloring f:F_2^6->{0,1} has a monochromatic affine
three-flat. Suppose otherwise. Each color class affinely spans F_2^6:
if one were contained in a proper affine hyperplane, the complementary
parallel hyperplane would lie entirely in the other color and contain
a three-flat.

Some three-flat has at least six points of one color. Indeed choose a
color with at least32 points. Its unordered pairs outnumber the63
nonzero differences. Two pairs with the same difference are disjoint
and form a monochromatic affine plane P. Its direction has16 cosets.
If each of the other15 cosets had at most one point of this color,
there would be at most4+15=19 such points. Some other coset therefore
has at least two. Its union with P is a three-flat with at least six.
Choose a three-flat with the greatest monochromatic count in either
color and rename that color zero. The count is6 or7, since8 would
already contradict the assumed counterexample.

If the maximum is6, an invertible affine map makes the six zeros in
this three-flat be the points0,1,2,3,4,5 and the two ones6,7, using
integer coordinate order with the first three coordinate bits varying.
An affine group on three bits is transitive on pairs of distinct points,
which proves this normalization. The zeros span all six dimensions,
so extend the first three basis directions by choosing zero points
outside the successive spans. This makes8,16,32 zero as well.
Every three-flat has at least two points of each color, by maximality.
For each of its eight points omitted, the other seven are therefore
not all zero and not all one. These are the 178560 regenerated
seven-literal clauses. The eight core units and three basis units give
178571 clauses on64 color variables. This is the first refuted case.

Now suppose the maximum is7. Make its three-flat the points0..7,
with its unique one at7. Choose a zero point outside this span and
make it8. In the coset8+F_2^3 the zero set contains no affine plane:
a zero plane there and the parallel plane in the original core avoiding
point7 would together form a zero three-flat. Any five points in F_2^3
contain a plane, since their ten pairs outnumber the seven nonzero
differences. Hence the coset has between one and four zeros. A four-set
without a plane is affinely independent. Affine permutations of F_2^3
are transitive on one-, two-, three-point sets and independent four-sets.

These permutations are realizable while preserving the original core.
For an invertible linear A on three bits, the map u->Au+b with
b=7+A7 fixes its unique one. An independent shear by coordinate bit3
adds any chosen vector on the coset8+F_2^3 and nothing on the core.
Thus the coset zero set may be exactly one of

    {0}, {0,1}, {0,1,2}, {0,1,2,4}

in its three-bit local coordinates. These are four cases, not an
assumption about an arbitrary coloring. They fix all first16 colors.
The canonical patterns include zeros0,1,2,4,8. The full-span property
allows a further zero outside this four-space to be chosen as point16.

Enumerate all65536 colorings of the coset16+F_2^4. Together with the
fixed first16 colors, require that each of the620 affine three-flats
in F_2^5 be mixed. Exact enumeration leaves respectively
3136,5760,6232,5296 patterns. It uses all five-bit three-flats, not just
coordinate flats. Truth tables are integer bit vectors and clause
satisfaction is evaluated exactly.

Use the following explicit subgroup to reduce these patterns. Enumerate
all168 invertible three-bit linear maps A, put b=7+A7, and all eight
shears R. On four bits (u,h) use

    (u,h) -> (Au+b+Rh,h).

Retain those preserving the fixed first16 colors. Their numbers in the
four cases are168,48,24,24. On the coset16+F_2^4, add an independent
translation T in F_2^4. It is the shear by coordinate bit4 in the full
space and fixes the first16 points. These maps form a group: the
invertible block affine maps are closed under composition, the
color-preserving subset is a stabilizer, and the independent shears
form its invariant translation subgroup. Every map is invertible and
preserves monochromatic affine subspaces.

The auditor generates all these permutations. For each still-unseen
pattern, it computes all its images, verifies they remain in the
allowed set, and chooses the least numerical image whose bit0 is zero.
Every allowed pattern has a zero, since an all-one coset would itself
contain a monochromatic three-flat. Translation makes bit0 zero, so
this minimum always exists. The resulting orbit counts are6,18,28,30,
with each original allowed pattern included. These computations give
82 explicitly regenerated representatives; no externally supplied
coloring-class table is assumed.

Only after this normalization choose a zero outside the first five-space
as point32. The zero class still spans six dimensions. The independent
shear by the last coordinate can make that point32 while preserving all
first32 colors. This ordering matters: no claim is made that the earlier
normalizations preserve a previously fixed point32.

For each of the82 representatives, fix its first32 colors and assign a
Boolean variable v_(a-31) to each point a=32..63. Regenerate all11160
three-flats in F_2^6. Each gives the two clauses requiring at least one
point of each color. Substitute the fixed colors, drop satisfied clauses
and fixed false literals, remove identical literal tuples, and sort by
length then tuple. Include the unit -v_1 for the chosen zero at32.
A counterexample would satisfy one of these82 residual CNFs, or the
maximum-six CNF above. This proves complete case coverage.

The affine-flat generator is the same exhaustive highest-pivot RREF
construction as elimination.md Section5.2, applied now to the full cubes. Increasing
pivots and offsets zero in pivot positions give each flat once; all
points are generated by extending the current affine span. The counts
620 and11160 are independently checked by the auditor.

## B.2 Check every finite refutation without a solver

the ramsey6 blob of refinement_certificates.txt is a UTF-8 JSON envelope containing base85/XZ
encoded untrusted witnesses. The decoded payload is29062 bytes. Its
header R6UP and one-byte count83 are followed by the cases in the exact
order above: the maximum-six case, then each four-case orbit enumeration
in order of the first still-unseen pattern. For each case unsigned
varints give the number of selected root clauses, their strictly positive
index differences in the regenerated list, and the number of derived
clauses. Each derived clause has a varint length followed by one-byte
literals encoded as2*(abs(v)-1)+(v<0). The final clause of each case is
empty. The decoder bounds lengths, counts, decompression and variable
ranges, checks distinct variables and exact end of data, and rejects an
early empty clause. Root selection does not introduce extra assumptions:
unsatisfiability of a subset implies unsatisfiability of the full CNF.

There are7544 selected root clauses and5185 derived clauses across the
83 cases. Each derived clause must pass exactly the reverse-unit-
propagation check proved sound in elimination.md Section5.2: negate its literals and
apply unit propagation to the previously checked database. Only a
contradiction licenses adding the clause. The empty clause must pass
this same test. Thus induction proves every accepted clause from the
regenerated roots, and the final contradiction proves the CNF unsatisfiable.
The runtime uses standard-library Python, no SAT solver, solver verdict,
proof hints, native executable, network or classification oracle. The
complete declared audit regenerates the groups, patterns, root clauses
and these refutations. This establishes Lemma S.

The six-variable statement is classical. For motivation see the normality
statement in Courtois et al., Theorem5.2,
https://discovery.ucl.ac.uk/id/eprint/10134213/1/Courtois_gen_constr_loops.pdf .
That paper's classification-based proof is not used here. All finite
case coverage and checked evidence needed here are supplied above.

## B.3 Deduce the finite circuit bound and prove sharpness

If n>=6, restricting any six input coordinates and fixing the others
would give a coloring of F_2^6. Lemma S supplies a constant three-flat,
contrary to the hypothesis. Thus every eligible function at d=3 has
n<=5. For n<=4 the right side4n-16 is nonpositive. For n=5 all five
inputs must be essential. Otherwise on the other four inputs one color
has at least eight points; its28 pairs outnumber the15 nonzero
differences, giving a constant plane. Adjoining the unused input gives
a forbidden constant three-flat. A binary circuit with five essential
inputs requires at least four gates: in the connected output component,
T essential inputs and s gates require at least T+s-1 distinct edges,
while the two input pins per gate give at most2s. Hence s>=T-1.
This covers arbitrary B2 tables, sharing, constants and terminal literals,
and proves (10) in every case.

The constant8 is sharp for d=3. Take

    h(x,y,u,v,z)=xy+uv+z.

Two ANDs and two XORs compute h. Its five inputs are essential, so the
connectivity argument gives the matching lower bound four. Suppose h
were constant on an affine three-flat with direction U. Its polar form
is B(a,b)=a_x*b_y+a_y*b_x+a_u*b_v+a_v*b_u. The four-point second
finite difference of a constant function is zero, so B vanishes on U.
Projection of U onto the first four coordinates is totally isotropic for
the nondegenerate displayed alternating form. For such a subspace W,
W is contained in its orthogonal complement; nondegeneracy gives
2*dim(W)<=4. Therefore U of dimension3 intersects the one-dimensional
kernel of this projection, the z-direction, nontrivially. But h toggles
along that direction, contradicting constancy. Thus h has no constant
three-flat and has circuit size four. Substituting n=5,d=3 in any
4n-8d+K bound gives4>=K-4, so K<=8.

This example proves sharpness8 at the fixed dimension3 and, together
with refinement.md, uniformly over all d>=3. No sharpness of any leading
circuit coefficient is asserted.
