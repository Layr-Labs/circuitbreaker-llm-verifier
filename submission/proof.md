# Self-contained 31/10 circuit lower bound

## Main theorem and proof map

Let G_N be the single total Boolean family defined in F4 below, with actual
input length N (field-description prefix included). Its yes-language
L={z:G_|z|(z)=1} is in P and hence in NP. In the full deterministic acyclic
B2 model with one output, all sixteen two-pin truth tables, arbitrary sharing,
depth and nonuniformity,

    liminf_(N->infinity) C_B2(G_N)/N >= 31/10.

This is a literature reproduction, not a new coefficient. All mathematical
premises needed for this theorem are supplied in this evidence package:
F1-F5 define the family and prove eligibility at every length; D1-D9 prove
its specialized affine-disperser property; C1-C7 prove the circuit reduction
and potential induction, whose complete local case proof is the companion
elimination.md, E0-E9. F6 combines them and proves the all-integer-length
asymptotic bound while charging every input bit. The proof is read together
with that explicitly named appendix; no paper's specialized theorem is taken
as an unproved premise. The elementary facts used are polynomial division and
root counting, finite-dimensional linear algebra, and finite sums/products;
their needed applications are derived below.

The scored conclusion is exactly the displayed liminf statement. Bounds inside
the proof are supporting lemmas for this family, not additional scored claims.
In particular the optional LP arithmetic does not establish the circuit theorem.
Sources, exact versions, and changes from the rejected submission are documented
in source-obligations.md.

The field description is genuine input; no deterministic field-construction
algorithm, advice string or truth-table selection is required for eligibility.
For the lower bound only, a valid prefix is fixed as a restriction. Its
sublinear length is included in the loss. A weaker claim is not being substituted
for the declared coefficient: the remainder divided by actual N tends to zero.

---


---

# F. One uniform family with a field description among its input bits

This module supplies an explicit family without invoking a deterministic
irreducible-polynomial construction theorem. It uses the disperser theorem proved
in module D and the circuit theorem proved in module C. The field description is
part of the actual input, and its complete length is charged in the final bound.

## F1. Elementary field facts and existence

A polynomial over F_2 has unique factorization, by the Euclidean algorithm.
An irreducible polynomial of degree a defines a field with 2^a elements by taking
the quotient F_2[Z]/(q). Multiplication by a nonzero u permutes its nonzero
elements. Taking their product and cancelling it gives u^(2^a-1)=1; hence every
element satisfies u^(2^a)=u. The map u -> u^2 is an automorphism: it is injective
in a field, hence bijective on a finite field.

For completeness, a field of size 2^a exists for every positive integer a. Form
a splitting field of X^(2^a)-X by adjoining roots of its irreducible factors one
at a time. Each adjoining step is a quotient by an irreducible polynomial and
there are finitely many steps. The polynomial has derivative 1, so has exactly
2^a distinct roots there. Its roots are closed under addition and multiplication,
because the 2^a-th power map respects both operations; every nonzero root has
its inverse among the roots. They therefore form a field K of size 2^a.

The fixed elements of the square map in any characteristic-two field are exactly
0 and 1, since X^2-X=X(X+1). If an element u in K has an orbit of length b under
squaring, then

    product_(j=0)^(b-1) (X-u^(2^j))

has coefficients fixed by squaring and therefore in F_2. It is the minimal
polynomial of u: a polynomial over F_2 vanishing at u vanishes at its entire
orbit, so cannot have degree less than b. This proves both irreducibility of this
orbit polynomial and that the degree of u equals its orbit length. In particular,
an element of degree b is a root of X^(2^a)-X exactly when b divides a.

When a is prime and a>1, K contains more than two elements. An element outside
F_2 has an orbit length dividing a and different from 1, hence equal to a.
Its minimal polynomial is an irreducible monic polynomial of degree a over F_2.
This proves the existence needed below. It does not choose a polynomial as advice
for an evaluator: the evaluator receives polynomial coefficients in its input.

For an odd prime a, define Tr(u)=sum_(j=0)^(a-1) u^(2^j). Its square equals itself
by cyclically shifting the summands. Thus Tr takes values in F_2; it is F_2-linear
and Tr(1)=1. The product pairing Tr(uv) is nondegenerate, since if u is nonzero,
choosing v=u^(-1) gives Tr(uv)=1.

## F2. A polynomial-time validity test for the supplied field

Let a be prime and q a monic polynomial of degree a. The following test is exact:

    q(0)=1, q(1)=1, and Z^(2^a) = Z modulo q.

The last congruence is computed by a iterations of modular squaring starting
from Z. If it holds, q divides X^(2^a)-X, a squarefree polynomial. Consequently q
is squarefree and every irreducible factor of q has degree dividing a, by the
orbit argument in F1. The first two tests exclude degree-one factors. As a is
prime, each remaining factor has degree a; since q itself has degree a it has
exactly one such factor and is irreducible. Conversely, irreducible q of prime
degree a has no root in F_2 and the quotient field has size 2^a, so passes all
three tests. The test is therefore equivalent to irreducibility.

Long multiplication, reduction modulo q, and the Euclidean algorithm have
polynomial bit cost in a. The displayed test uses a modular squarings and two
bit sums, so is deterministic and polynomial-time, without a search for q.

## F3. An elementary prime-growth estimate

Only the following weak estimate is needed: for all sufficiently large integers
m there is a prime between m and 4m. Here is a proof to avoid a prime-distribution
premise. Write P(s) for the product of all primes at most the integer s. We have
P(s)<=4^s. Induct on s. The even case s=2v>2 follows from P(2v)=P(2v-1). In the
odd case s=2v+1, every prime in (v+1,2v+1] divides binom(2v+1,v). The two central
binomial coefficients are equal, hence each is less than 2^(2v), so

    P(2v+1) <= P(v+1)*binom(2v+1,v) < 4^(v+1)*2^(2v)=4^(2v+1).

The cases s=1,2 start the induction.

Suppose there is no prime in (m,4m]. The integer B=binom(4m,2m) satisfies
B>=2^(4m)/(4m+1), since it is the largest of the 4m+1 binomial coefficients.
For every prime p, its exponent in B is

    sum_(j>=1) (floor(4m/p^j)-2*floor(2m/p^j));

each summand is zero or one, so p to that exponent is at most 4m. The contribution
from primes p<=2*sqrt(m) is at most (4m)^(2*sqrt(m)). For larger primes no square
power contributes. By the assumed missing interval all of these primes are at
most m, and their contribution is at most P(m)<=4^m. Thus

    4^m <= (4m+1)*(4m)^(2*sqrt(m)).

Taking logarithms and dividing by m contradicts the limit
log(m)/sqrt(m) -> 0. Therefore the asserted interval contains a prime for every
sufficiently large m. This estimate is weaker than Bertrand's bound but suffices.

## F4. Definition at every actual input length N

For N=0, set G_N(empty)=0. For each N>=1, perform the following fixed algorithm.
Search among the integers a<=floor(N/2) for the least prime satisfying

    a^5 > 32*N^3.

Primality is decided by trial division. If there is no such prime, output 0 for
all inputs of that length. Otherwise put k=a, n=N-k, r=ceil(n/k), and
 t=ceil(sqrt(r)), using exact integer arithmetic. If t>=r or
k*t<r*(r+1), again output 0 for that length. These are explicit small-length
conventions and are eventually inactive, as F5 proves.

Split the N-bit input into k prefix bits followed by n data bits x. Interpret
those k bits as the coefficients q_0,...,q_(k-1) of

    q(Z)=Z^k+sum_(j=0)^(k-1) q_j*Z^j.

If q fails F2, output 0. If it passes, work in K=F_2[Z]/(q) using its polynomial
basis. Append r*k-n zero bits to x, partition the resulting string into r blocks
of k bits, and read each block as a field element x_i, with its first bit the
constant coefficient. Let

    e_i=r+1-i,  alpha_i=1+2+...+2^(e_i),
    H(x_1,...,x_r)=sum_(I subset {1,...,r}, |I|=t) product_(i in I) x_i^(alpha_i).

Define G_N(q,x)=Tr(H(x_1,...,x_r)). The trace is the k-term sum in F1. This is
one Boolean output, not a vector of field coordinates. The exponents are in the
descending order used in the proof of module D.

The preceding steps define one total language at every integer length. The
polynomial q is ordinary input, not a witness, oracle, or chosen language
parameter. Invalid encodings have the explicitly defined value zero. For any
valid q the remaining n-bit restriction is the stated field polynomial.

## F5. Uniform complexity and parameters

The prime-growth estimate of F3, applied to ceil(2*N^(3/5)), implies that the
searched prime eventually exists and k=O(N^(3/5)); in particular k<=N/2 for all
sufficiently large N. Its defining strict inequality gives k>2*N^(3/5), so

    n=N-O(N^(3/5)), r=O(N^(2/5)), t=O(N^(1/5)), and r tends to infinity.

Consequently t<r eventually. Also

    r*(r+1)/t <= r^(3/2)+sqrt(r)
                    <= (1/(2*sqrt(2)))*N^(3/5)+o(N^(3/5)) < k,

using r<=N/k+1<(1/2)*N^(2/5)+1. Thus the second exceptional-length test is also
eventually inactive. Every sufficiently large integer length is covered; no
subsequence is selected.

All searches are explicitly bounded by N. Modular arithmetic costs polynomial
in k, the largest exponent alpha_i has at most r+1 bits, and binary powering
uses O(r) field multiplications per block. The elementary symmetric polynomial
is computed in O(r*t) field operations by descending updates

    E_0=1, E_j=0 (j>0); E_j <- E_j + y_i*E_(j-1), j=min(i,t),...,1.

Splitting subsets according to whether they contain i proves this recurrence.
It uses no division in characteristic two. The trace uses k squarings and adds.
Thus one deterministic algorithm computes G_N in time polynomial in N, at all
lengths including invalid encodings and the exceptional lengths. Define an NP
verifier to reject nonempty witnesses and otherwise run this evaluator. With
witness bound p(N)=0, membership gives an accepted empty witness, and acceptance
implies that the evaluator returned one. This establishes both directions of
NP membership and in fact P membership.

## F6. Hard restriction, all bits charged, and final coefficient

For every sufficiently large N, fix any irreducible q of degree k; F1 proves at
least one exists. This is a restriction used in a lower-bound argument, not an
algorithmic choice in the language. The full circuit must be correct on every
prefix, including that q. The zero-padding map from F_2^n to K^r is injective and
preserves affine dimensions. Module D proves that the restricted n-bit function
is nonconstant on every affine subspace of dimension greater than

    B=k*r/t+k*t+r*(r+1).

Set d=floor(B)+1. The parameter bounds above imply d=O(N^(4/5))=o(n), since
k*r/t<=k*sqrt(r), k*t<=k*(sqrt(r)+1), and r^2+r=O(N^(4/5)).

Specializing k prefix inputs of a B2 circuit to constants cannot increase its
number of gates. Gates fed by a constant become legal trivial or unary truth
tables; repeated pins and output degeneracies are handled in C1. No field gate
is introduced in this Boolean circuit transformation. Therefore module C's
n-bit bound implies

    C_B2(G_N) >= (31/10)*n-(261/10)*d-128/5
               = (31/10)*N - [(31/10)*k+(261/10)*d+128/5].

The bracket is O(N^(4/5))=o(N). In epsilon form, for each epsilon>0 choose a
threshold after all the proved parameter conditions hold and after that bracket
divided by N is at most epsilon. Every integer N above that threshold then
satisfies C_B2(G_N)>=(31/10-epsilon)*N. Neither the family, its evaluator nor its
encoding depends on epsilon. The field prefix, all data bits and all small-length
conventions have been accounted for. This reduction claims only the finite
coefficient 31/10 and does not assert a superlinear endpoint.

---

# D. The affine-disperser theorem, with supporting arguments

The presentation specializes the algebraic construction to characteristic two.
All statements used below are proved here or in the elementary field facts of F.
The source of the construction is Ben-Sasson and Kopparty's 2010-03-12 preprint,
Sections 3, 6 and 7; the citation is attribution, not an unproved premise.

## D1. Polynomial representations of linear maps and spaces

Let K have 2^k elements. A polynomial of degree less than 2^k represents the zero
function on K only if it is zero: a nonzero polynomial has at most its degree
many roots, as follows by repeatedly factoring X-a at a root. By induction on
variables, a polynomial of degree less than 2^k in each variable is zero on K^r
only if every coefficient is zero. Division by X_j^(2^k)-X_j in each variable
therefore gives a unique representation of every represented function. We call
this the reduced polynomial. A positive-degree coefficient in its reduced
polynomial proves that a function is nonconstant.

A linearized polynomial is sum_(j=0)^(k-1) a_j*X^(2^j). It defines an F_2-linear
map K->K. Conversely every such linear map has exactly one such representation:
the coefficient space has F_2-dimension k^2, as does the space of linear maps, and
evaluation is injective by the preceding root bound. Compositions and linear
combinations of linear maps thus have reduced linearized representations.

For a linear space U<=K of dimension d, put P_U(X)=product_(u in U)(X-u).
This is monic and linearized, of degree 2^d, with nonzero coefficient of X.
To prove this, induct on a basis: for U=W+F_2*v with v not in W,

    P_U(X)=P_W(X)*P_W(X+v)=P_W(X)^2+P_W(v)*P_W(X).

The induction gives linearization and monicity; P_W(v) is nonzero and preserves
the nonzero linear coefficient. The kernel is exactly U. The image W'=P_U(K)
has dimension k-d. Let Q_U=P_(W'). Then

    Q_U(P_U(X))=X^(2^k)-X.

Both sides are monic of degree 2^k and have exactly all elements of K as roots.
Also P_U(Q_U(X))=X^(2^k)-X. Indeed, composing the difference with P_U gives
zero, because P_U(X^(2^k)-X)=P_U(X)^(2^k)-P_U(X); a nonzero polynomial composed
with the nonconstant P_U cannot be zero, by the product rule for degrees.
It follows that Q_U(K)<=U. Since Q_U has kernel W' of dimension k-d, its image
has dimension d and equals U. For d=0 use its reduced zero map. For d>=1
its degree is less than 2^k, so write

    Q_U(X)=sum_(j=0)^(k-1) a_j*X^(2^j).

Then a_(k-d)=1, a_j=0 for k-d<j<k, and a_0 is nonzero. The polynomial Q_U has
exactly 2^(k-d) distinct roots. No algorithm for finding Q_U is needed: it is
used to analyze an arbitrary subspace, not to evaluate the Boolean family.

## D2. The cyclic zero-block lemma

Assume d>=1. Indices of a_j are read modulo k. No d consecutive coefficients of Q_U are all
zero. If they were, raising Q_U to an appropriate power of two and reducing
modulo X^(2^k)-X would rotate those zeros to the highest d positions. The result
would be a nonzero polynomial of degree at most 2^(k-d-1), still vanishing on
the 2^(k-d) roots of Q_U, contradicting the root bound. The rotation formula is

    (Q_U(X))^(2^s) reduced = sum_j a_(j-s)^(2^s)*X^(2^j).

Suppose additionally that k is prime and d>1. A block of d-1 zeros can occur
only in the displayed terminal block k-d+1,...,k-1. Otherwise rotate that block
to those positions by a nonzero shift s in {1,...,k-1}. The rotated polynomial
Q' has degree at most 2^(k-d) and the same 2^(k-d) roots as Q_U. Comparing leading
coefficients and subtracting shows Q'=c*Q_U for a nonzero c in K; a smaller-degree
nonzero difference could not have so many roots. Therefore every y in U satisfies

    y^(2^s)=c*y.

For a fixed nonzero b in U, each y/b is fixed by the 2^s-th power map. It is also
fixed by the 2^k-th power map. Since gcd(s,k)=1 and the square map is invertible,
Bezout's identity for s and k shows y/b is fixed by the square map itself. Thus
y/b belongs to F_2, so U<=F_2*b and dim U<=1, a contradiction. This proves the
stronger zero-block assertion. In particular, any zero block between consecutive
nonzero coefficients outside the terminal block has length at most d-2.

## D3. Binary weight and a coefficient formula

For a nonnegative integer b<2^k, let wt(b) be its number of binary ones. In a
reduced monomial, applying a power of two rotates binary positions modulo k and
does not increase weight. Multiplication does not increase the sum of weights.
A direct way to check the last statement is to represent an exponent as a
multiset of powers 2^j, wrap indices modulo k, and replace two copies at j by
one copy at j+1. Every collision decreases the number of copies. This process
terminates with distinct positions, including the all-ones exponent 2^k-1 when
appropriate, and gives the reduced positive exponent. Thus its weight is no
larger than the number of original factors.

Let Q_i(X)=sum_j a_(i,j)*X^(2^j), and define its k-by-k matrix

    M_i[u,v]=a_(i,v-u)^(2^u),   u,v in {0,...,k-1}.

This is the coefficient of X^(2^v) in Q_i(X)^(2^u), reduced. Let A_i and B_i be
sets of equal size h_i with 1<=h_i<k. Put alpha_i=sum_(u in A_i)2^u and
beta_i=sum_(v in B_i)2^v. In

    R=Tr(mu*product_i Q_i(X_i)^(alpha_i)), reduced,

the coefficient of product_i X_i^(beta_i) equals

    sum_(s=0)^(k-1) mu^(2^s)*product_i Perm(M_i[A_i,B_i-s])^(2^s).       (D.1)

Here Perm is the sum over bijections between rows and columns of the product of
the selected entries; indices in B_i-s are modulo k. To prove (D.1), expand the
trace and then expand each product of h_i linearized factors. A target exponent
of weight h_i<k can arise only when all h_i selected binary positions are
distinct: any collision would lower weight. Consequently each contribution is
one bijection of those positions to the h_i factors, exactly the permanent.
The variables X_i are separate, so their coefficients multiply. Summing over
trace shifts proves the formula, including all reduction effects.

If at shift zero all these minors are triangular with nonzero diagonals, their
permanents are nonzero. If at every other shift at least one minor has a zero
column, all other summands in (D.1) vanish. The coefficient is then nonzero for
every nonzero mu. This observation will be used with explicit column sets.

## D4. A nonzero coefficient for independent subspaces

Let U_i<=K have dimensions d_i, with k prime, and choose integers

    1<=e_i<d_i-1,  h_i=e_i+1,  w_i=d_i-h_i>0.

Let Q_i be the image polynomial from D1. Suppose

    sum_i w_i >= k-max_i d_i+1.                                      (D.2)

Then Tr(mu*product_i Q_i(X_i)^(1+2+...+2^(e_i))) has, in its reduced polynomial,
a nonzero coefficient product_i X_i^(beta_i) with wt(beta_i)=e_i+1.

Here is the complete choice of columns proving this assertion. Reorder the
independent factors so that d_1=max_i d_i. Let j_* be the first index for which
sum_(i<=j_*) w_i >= k-d_1+1. For i<=j_*, put L_i=sum_(j<i)w_j, so

    0<=L_i<=k-d_1<=k-d_i.

For i>j_*, choose B_i={0,k-d_i+1,...,k-d_i+e_i}. These give triangular minors
at shift zero and need not be used to cover other shifts.

For i<=j_* there are two choices. If a_(i,L_i) is nonzero, choose

    B_i={L_i,k-d_i+1,...,k-d_i+e_i}.

If it is zero, let b<L_i<c be the nearest nonzero coefficient indices on either
side. They exist because a_(i,0) and a_(i,k-d_i) are nonzero. Their intervening
zero block is not the terminal block. Thus q=c-b-1<=d_i-2 by D2. Put
s_i=min(q,e_i) and choose the ordered columns

    B_i=(c, b+1,...,b+s_i, k-d_i+s_i+1,...,k-d_i+e_i).

These columns are distinct: c<=k-d_i, the middle block lies before c, and the
last block lies after k-d_i. Use row set A_i={0,...,e_i}.
In the first choice, entries above the diagonal in any column after the first
are terminal zero coefficients; the diagonal is a_(i,L_i) followed by Frobenius powers
of the nonzero a_(i,k-d_i). In the second choice the first s_i+1 diagonal
entries are a_(i,c), then Frobenius powers of a_(i,b), all nonzero. Above-diagonal entries
in that block belong to the gap b+1,...,c-1. The remaining columns have terminal
zeros above the diagonal and Frobenius powers of a_(i,k-d_i) on the diagonal. Thus all minors are
triangular with nonzero diagonals.

It remains to eliminate every shift v in {1,...,k-1}. A column indexed by u in
B_i becomes all zero in the minor M_i[A_i,B_i-v] whenever

    v in [u+1,u+w_i] modulo k,                                     (D.3)

because u-v, u-v-1,...,u-v-e_i then lie in the terminal zero block.
For i=1 the standard columns therefore cover

    [1,w_1] and [k-d_1+2,k-1].

For any i<=j_* we claim the union of intervals (D.3) covers
[L_i+1,L_i+w_i]. If a_(i,L_i) is nonzero this follows from u=L_i.
If it is zero and q<e_i, the middle block together with c includes L_i, so the
same argument applies. If q>=e_i, the columns b+1,...,b+e_i cover the interval
[b+2,b+d_i-1]. This contains c, since q<=d_i-2. Adding column c extends the
covered interval through c+w_i, with no gap. Because b<L_i<c, the union covers
[L_i+1,L_i+w_i] as claimed. All these intervals are ordinary integer intervals
before wrapping: L_(j_*+1)<=k-2, since L_(j_*)<=k-d_1 and w_(j_*)<=d_1-2.

The consecutive covered segments therefore reach at least k-d_1+1; the first
factor's terminal segment covers k-d_1+2 through k-1. Every nonzero shift is
covered. Formula (D.1) now has exactly one surviving, nonzero summand. This
proves the independent-subspace assertion with all coefficient choices specified.

## D5. Triangular decomposition of an arbitrary affine space

Let A be an F_2-affine subspace of K^r, and choose gamma in A. Successively
project A-gamma on the first block and choose a linear section of that projection;
then apply the same operation to its kernel, which has first block zero. Continue
through the blocks. Linear sections exist by extending a basis of each kernel to
a basis of its domain. This constructs spaces U_i<=K and linear maps sigma_ij,
j<i, with

    A={ (gamma_i+y_i+sum_(j<i) sigma_ij(y_j))_(i=1)^r : y_i in U_i },
    dim A=sum_i dim U_i.

The parameterization is one-to-one: at the first differing y_i its i-th block
differs, after all earlier blocks have been fixed. Thus the dimension identity
also follows directly. This is block Gaussian elimination, not an assumption of
independence of the original coordinates.

Let Q_i map K onto U_i as in D1. The maps sigma_ij composed with Q_j have reduced
linearized polynomials Q_ij by D1. Define

    R_i(X_1,...,X_i)=gamma_i+Q_i(X_i)+sum_(j<i)Q_ij(X_j).

Their joint image is exactly A. Constants and zero-dimensional U_i cause no
problem: Q_i is then the zero map in this parameterization.

## D6. Preservation of the highest binary weights

Choose I={i_1<...<i_t} and positive integers e_i. Let alpha_i=1+2+...+2^(e_i).
In the reduced polynomial

    Tr(mu*product_(i in I) R_i(X_1,...,X_i)^(alpha_i)),

the coefficient of a monomial product_(i in I) X_i^(beta_i) with
wt(beta_i)=e_i+1 equals its coefficient with every R_i replaced by Q_i(X_i).
Expand a fixed trace summand into its e_i+1 Frobenius-linear factors per i.
Only the factors belonging to the largest index i_t can contain X_(i_t).
To reach its required weight e_(i_t)+1, all those factors must choose their
Q_(i_t) part, with no constant or earlier-variable part. D3 shows that reduction
cannot restore weight lost by any other choice. These chosen factors involve no
earlier variable. Repeat for i_(t-1), and then downward. Exactly the independent
product remains. The argument applies separately to every trace shift, and proves
the stated coefficient identity, including e_i=1.

## D7. A resource-assignment lemma excluding competing terms

Fix r and let p_i=r+2-i. Suppose I,J are t-element subsets of {1,...,r}.
For each j in J imagine p_j objects, each assigned to an integer lambda<=j.
If every i in I receives at least p_i objects, then I=J.
It suffices to allow lambda in {1,...,j}; a constant term in a polynomial expansion
can be assigned arbitrarily to 1 and only increases available counts.

Counting all objects gives sum_(i in I)p_i<=sum_(j in J)p_j, hence
S_1=sum_(i in I)i-sum_(j in J)j>=0. For any positive integer h, assignment and
lambda<=j give

    sum_(j in J) p_j*j^h >= sum_(i in I)p_i*i^h.

Substitute p_i=r+2-i and rearrange to get S_(h+1)>=(r+2)*S_h. By induction all
S_h are nonnegative. This forces max I>=max J: otherwise for large h the largest
power on the J side exceeds the entire I sum. Assignment also forces max I<=max J,
since a larger target could receive no object. Let their equal maximum be m.
Only source m can feed target m, and its capacity equals the demand p_m, so all
its objects are assigned to m. Remove that source and target and apply the same
argument to the remaining equal-size sets. Induction proves I=J.

## D8. The arbitrary-affine-source theorem and its proof

Let k be prime and let 1<=t<r satisfy k>=r*(r+1)/t. These hypotheses imply k>r+1, so all e_i+1<k as needed in D3. Put

    e_i=r+1-i, alpha_i=1+2+...+2^(e_i),
    f(z_1,...,z_r)=Tr(mu*sum_(|I|=t) product_(i in I) z_i^(alpha_i)), mu!=0.

We claim f is nonconstant on every affine A<=K^r with

    dim A > B=k*r/t+k*t+r*(r+1).                                   (D.4)

Decompose A by D5, and write d_i=dim U_i. Let S={i:d_i>r+1}.
If |S|<t, then sum_i d_i<=k*t+r*(r+1), contradicting (D.4). Also

    sum_(i in S) d_i >= dim A-r*(r+1) > k*r/t+k*t.

Take I_0 to be indices of the t largest d_i among S. Since |S|<=r,

    sum_(i in I_0)(d_i-(r+1))
       > k+k*t^2/r-t*(r+1) >= k.

For i in I_0 we have 1<=e_i<d_i-1, and replacing r+1 by e_i+1 only increases
that sum. The independent-subspace result D4 thus applies. It gives a monomial
M=product_(i in I_0) X_i^(beta_i), wt(beta_i)=e_i+1, with nonzero coefficient in
Tr(mu*product_(i in I_0)Q_i(X_i)^(alpha_i)). By D6 its coefficient is the same
and nonzero in the corresponding term with R_i in place of Q_i.

No other t-subset J contributes to M. To see this, expand its trace summands
and then each R_j^(alpha_j) into p_j=e_j+1=r+2-j Frobenius-linear factors. A
chosen monomial from any such factor is a constant or a single-variable binary
power with variable index at most j. To reach M, variable i requires at least
wt(beta_i)=p_i such factors, since reduction cannot increase binary weight.
Assign each factor to its chosen variable, and assign constants to 1. The
resource-assignment lemma D7 forces J=I_0. Hence the nonzero coefficient just
found cannot cancel against another term of the symmetric sum.

The reduced polynomial of f composed with the surjection onto A therefore has
a positive-degree monomial with nonzero coefficient. By D1 it is a nonconstant
function. Its image equals f(A), proving (D.4). This proves the disperser theorem
rather than invoking its published statement as a premise.

## D9. Application to the actual data bits

In F, k is odd prime, mu=1, and the tested parameters satisfy the hypotheses of
D8. Zero padding and the polynomial basis give an injective F_2-linear map from
the n data bits to K^r. An affine subspace of the n-bit space keeps its dimension
under that injection, so D8 applies directly to its image. With d=floor(B)+1,
the n-bit restriction is nonconstant on every affine subspace of dimension d.
The polynomial description prefix is handled by the circuit restriction argument
in F6; it is never ignored when counting the full input length N.

---

# C. Circuit lower bound: model, invariants, and induction

The local elimination proof is in module E. This module proves its algebraic,
combinatorial and accounting prerequisites, including the initial potential
bound. The source of the elimination strategy is Li-Yang, ECCC TR21-023,
2021-02-20, Sections 2-4. No circuit theorem from that source is assumed here.

## C1. The complete gate model and normalization conventions

A gate has two input pins and an arbitrary Boolean truth table. Every such table
has algebraic normal form c+a*u+b*v+d*u*v over F_2. The sixteen tables split into:
two constants, four projections or negated projections, two functions u+v+c,
and eight functions ((u+a)*(v+b))+c. Call the last two classes XOR and AND type.
All retained gates cost one, including constants and unary gates; wires, inputs,
output selection and constant nodes introduced by restrictions cost zero.

The single output is a distinguished terminal wire. For use counts, include its
one use of its source node. The terminal is not a gate and does not compute an
extra function. Thus a gate with no uses is genuinely dead, whereas a selected
output gate is not deleted as dead. Repeated pins are counted as two uses until
normalization replaces the gate by its equivalent unary or constant table.

A unary gate can be bypassed by modifying each affected descendant truth table.
Complementing a pin does not change arity or add a gate. If the selected output is
a negated projection of another gate, complement that gate's value and all its
other occurrences, including its defining truth table, then select it as output.
This is a change of a variable representing a gate value, not a free inverter.
If the predecessor is a primary input, a costed negated projection must be
retained. The live disperser states used below cannot be functions of one free
input or constants, as C2 proves, so that terminal exception is never needed in
the credited elimination steps.

We use semicircuits during the proof: an XOR-only prefix may have directed
cycles, followed by an acyclic part containing the nonlinear gates. If its gate
values are v, its equations are A*v=b(x) over F_2. Require A invertible, so every
assignment of the free inputs determines exactly one value of every prefix gate.
Call this condition fairness. All nonlinear descendants are then evaluated in
acyclic order. Original acyclic B2 circuits are a special case. Gate counts and
truth tables are unchanged by this enlargement of the proof's state space.

## C2. Allowed restrictions and survival of the function

Partition the original n data variables into free variables F, affine variables
L, and quadratic variables Q. Each variable in Q is an AND-type function of a
pair of distinct free variables, with different quadratic equations using
disjoint pairs. Each variable in L is an affine function of the free and quadratic
variables. Only free variables have uses in the semicircuit. A free variable in
a quadratic pair is protected; its mate is its couple. Other free variables are
unprotected. This is a read-once quadratic source. It is parameterized
one-to-one by its s=|F| free bits.

Fix one free bit in each of the q quadratic pairs. All remaining equations become
affine, so the source contains an affine subspace of dimension s-q>=ceil(s/2).
An n-bit affine disperser for dimension d therefore remains nonconstant whenever
s>=2d. If its value on the source depended on only one free bit u, fix u as well,
choosing u itself as the fixed member of its quadratic pair if it is protected.
If unprotected, q<=floor((s-1)/2). In either case the resulting constant affine
subspace has dimension at least d when s>=2d. Thus a live source of that size
cannot have a constant or one-input output.

If the source's value is an affine function of the free bits, first fix one bit
per quadratic pair and then impose the one affine equation that fixes this
function's value. The remaining affine subspace has dimension at least
s-q-1>=ceil(s/2)-1. Hence it cannot be affine when s>=2d+1. The same argument
applies to an AND-type function of two affine forms. On the affine subspace
obtained after fixing the quadratic pairs, if a form is nonconstant, fixing it
to its killing value is one consistent affine equation. If both forms are
constant the output is already constant; if one is a nonkilling constant, use
the other nonconstant form. Thus at most one consistent affine equation is
needed. Consequently a topmost nonlinear
gate with only affine ancestors cannot itself be the output in the live states
where s>2d+2. It has a genuine gate descendant, not just the terminal wire.

A constant restriction on a free bit removes one free dimension. If protected,
its quadratic equation becomes affine and its couple becomes unprotected.
An affine restriction on an unprotected bit, solved as that bit equal to an
affine function of the others, also removes one dimension and preserves the
source form. An affine identification of a protected couple is allowed: its
quadratic expression becomes affine in their one remaining bit. Finally, an
unprotected bit may be made an AND-type function of two other unprotected bits,
provided its circuit occurrences can be rewired to that function without a
nonlinear cycle. It becomes a quadratic variable and the two right-hand bits
become a new protected couple. Disjointness is preserved. All four operations
parameterize a subset of the previous source with exactly one fewer free bit.
Their validity is explicit in the local cases in E; witnesses play no role.

## C3. Packing and integer measure

A troubled gate is an AND-type gate with exactly one use and with two distinct
free inputs, each having exactly two uses. A troubled pair consists of two such
gates sharing a free input. A packing is a set of pairs with no gate used twice.
Let T be the number of troubled gates and P the packing; define

    Phi=T-|P|,       M=5*g+48*i+9*q+Phi.

Here g counts gates, q counts quadratic equations, and i counts free inputs that
have a use or are protected. All terms are nonnegative. A decrease is written
Delta M=M_before-M_after; an increase of Phi is written dPhi=Phi_after-Phi_before.
Thus deleting a gates, losing b influential inputs and losing c quadratic
equations gives Delta M=5*a+48*b+9*c-dPhi. A new quadratic equation has c=-1.

Retain every old packed pair whose two gates remain troubled and adjacent.
If one ceases to be troubled, removing its pair causes no positive potential
change: two old gates and one pair contributed one, and at most one unpaired
old gate remains. A pair cannot separate while both gates remain troubled under
the local normalizations: its shared free input had exactly those two gate uses,
so changing either shared occurrence also changes or deletes one of the pair's
gates. Newly created adjacent troubled gates can be paired.

The following normalization bounds hold, with any additional losses of i or q
only helping:

- Deleting a dead gate frees at most its two predecessor uses. Each affected
  predecessor can itself become one troubled gate, or, if a free input, can
  newly feed at most two troubled gates that can be packed together. Thus
  dPhi<=2 and Delta M>=3.
- A gate trivialized by one constant pin can be replaced by its constant value.
  Only its other predecessor loses a use; hence dPhi<=1 and Delta M>=4.
- A gate made unary by a constant pin is bypassed. Only its other predecessor
  can become or newly feed troubled gates. The same bound dPhi<=1 holds.
  In a cyclic prefix this is Gaussian elimination of the gate variable using
  its equation with unit coefficient on its left side. A self-occurrence with
  the other pin constant would give a zero matrix row and contradict fairness.
- If a one-use gate G and its consumer H share another input u, substitute the
  two-input truth table of G into H. H remains a two-input Boolean table and
  G is removed. The shared input u loses one occurrence; the other input of G
  keeps its occurrence, now at H. New trouble can only be H itself or gates
  newly fed by u; if two appear they share u and can be packed, so dPhi<=1.
  In an XOR prefix the equation defining G has coefficient one on G (there
  is no self-pin after C5's constant cleanup). Eliminating that equation and
  substituting into H is a Schur-complement pivot with determinant unchanged.
  It therefore preserves fairness, including when other XOR paths form cycles.
  In the acyclic part it is direct truth-table composition. The result is still
  a two-pin table; any repeated or trivial pin is normalized in turn.
- If a gate's two pins coincide, its table becomes constant or unary in that
  one predecessor. Replacing the duplicate pin by a constant and using the
  preceding rules affects only that predecessor, so dPhi<=1. Constants and
  unary output cases follow C1 and the live-state exclusion in C2.

Every such normalization strictly decreases g, hence terminates, and has
nonnegative Delta M. None assumes a tree: each shared predecessor is a single
node and its actual use count is updated in place.

### Coincident consumers during constant propagation

Here is the edge count used for coincidences in the local cases. In an affected
acyclic cone whose vertices were live before rewiring, let K be the removed gates
whose final values are constant, and U the other removed gates. Let h be the
initial constant-pin edges entering this cone. Every original outgoing edge
from K ends at a removed gate: otherwise a retained gate would have a constant
pin and would still require normalization. No K vertex is the output, because
the function is nonconstant. Each K vertex originally had at least one outgoing
edge, since the cone was trimmed to the output. Counting constant-pin edges
with multiplicity gives

    h + number of edges from K <= 2*|K|+|U|,

as a nonconstant removed gate has at most one constant predecessor. Hence
h<=|K|+|U|. This accounts for coincident consumers: two constant pins at one
consumer make it constant, and its onward edge pays for the merged occurrence.
Acyclicity ensures that the classification by final constant values is
well-defined even if the gates were removed in a different order.

This count is used for coincidences among nonlinear gates and their descendants,
which lie in the acyclic part. It does not assert that the XOR prefix is acyclic.
Prefix counts in E8 use the distinct immediate consumers of the reconstructed
constant, including the reversed-path predecessor (C5). Additional credited
consumers come from the killed nonlinear G; any coincidence with an already
counted consumer lies in G's acyclic descendant cone and is handled by the
displayed edge count. An unused newly added reverse-path gate is never needed
to supply an extra descendant count. Thus sharing cannot turn the named local
degree accounts into a formula assumption.

Each constant/unary deletion costs at most one potential unit. When its other
predecessor is a constant, a disappearing gate, or a free input of at most two
uses which drops to at most one, its sharper positive increase bound is zero.

## C4. Initial potential without an external gate-count lemma

Start with a gate-minimal acyclic circuit for the n-bit disperser. Trim all gates
not reaching its output and normalize repeated pins and redundant unary tables;
this cannot increase g. Choose a maximum packing of its troubled gates.
Form a multigraph whose vertices are troubled gates and whose edges are their
shared free input variables. Each vertex has degree at most two: its two inputs
each have at most one other use. Parallel edges are allowed, but loops are not,
since the pins have been normalized to distinct inputs. Components are paths,
cycles, isolated vertices, or the two-edge compact cycle.

A path with L vertices uses L+1 distinct input variables and has a matching of
size floor(L/2). Its potential ceil(L/2) is at most (L+1)/2. A cycle with L
vertices uses L distinct variables and has potential ceil(L/2). Thus an even
cycle is charged by half its inputs; an odd cycle has an excess of one half.
Different components use disjoint variables. Let b be the number of odd cycles.
Then Phi<=n/2+b/2.

In each odd cycle choose one variable to remain free. Its two incident troubled
gates have distinct other input variables, since an odd cycle has at least three
vertices. Fix each of those two neighbors to the respective killing value for
its AND-type gate, and fix all other variables of that cycle. All its gate
outputs are constant. Do this independently in every odd cycle, and fix every
input outside them. Each retained free bit has exactly its two incident gate
uses, both killed, so no free bit can affect the circuit output. Acyclic
evaluation therefore makes the output constant on an affine subspace of
dimension b. The disperser condition forces b<d. Consequently

    Phi<=n/2+(d-1)/2<=n/2+5*d/2,
    M_initial<=5*g+48*n+n/2+5*d/2.                                (C.1)

This proves the needed initial bound, with slack, including compact pairs and
arbitrary use sharing. The maximum packing exists by matching each path or
cycle; its computation is not an assumption about the Boolean family's evaluator.

## C5. Constant gates, dead cyclic components, and reconstruction

All XOR-prefix values are affine functions of the free bits, obtained by solving
A*v=b(x). A coefficient of x in a gate value can be nonzero only if there is a
path from x to that gate. Indeed, nodes unreachable from x form a closed
subsystem; ordering them before the reachable nodes makes A block triangular,
so their block is invertible and their values cannot change when x changes.
The same argument applies to a new right-hand forcing term in one gate equation.

If the affine value of a prefix gate I contains an unprotected free bit x, choose
a simple path x=I_0 -> I_1 -> ... -> I_l=I. Each defining equation on it has form
I_j=I_(j-1)+T_j+c_j. Reverse these equations to
I_(j-1)=I_j+T_j+c_j, fix I=b, and make the old input x a gate Z. Since the affine
coefficient of x in I was one, the added condition I=b solves uniquely for x in
terms of the other free bits. The reversed system is equivalent to that original
system under this condition, and has a unique solution. It is therefore fair.
Exactly one old gate I is removed and one gate Z added, so g is unchanged.
The unprotected input x becomes affine, i decreases by at least one, and q does
not increase. Side inputs along the path are retained once each; internal path
nodes remain XOR type. At an empty packing the reconstruction itself creates
no troubled gate: a primary input is replaced by a gate, and changes of XOR
orientations do not turn a gate predecessor into a free input. A free input
with two existing troubled-gate uses could not also occur on this XOR path.
This is the reconstruction used in E.

Semantic constant gates in a cyclic prefix need some care: deleting a guessed
constant equation alone need not preserve rank. Here is a rank-preserving
elimination. Suppose gate I has the constant value c. Select a row J with
(A^(-1))_(I,J)=1. Such a row exists because A is invertible. There is a directed
path from J to I, by the forcing-term reachability observation above. Delete row
J, fix I=c, and orient the equations on a simple J-to-I path backwards, leaving
all other equations unchanged. The resulting square matrix on the remaining
gate variables is invertible: its determinant is the nonzero cofactor obtained
by deleting row J and column I from A. The old full solution satisfies it, so
its unique solution is exactly the old solution with I omitted. Hence it is a
fair semicircuit with one fewer gate. All reoriented gates remain XOR type.
Only the two predecessor uses from the deleted row J can lower free-input use
counts. The same potential accounting as C3 gives dPhi<=2, so Delta M>=3.
This supplies a valid constant elimination even when the prefix has cycles.

In the acyclic part a semantic constant gate can simply be replaced by its value,
with the same dPhi<=2 bound. A gate whose output is a unary function of one pin
can also be removed there by changing descendant truth tables, with dPhi<=2.
In the XOR prefix a nonconstant unary dependence on one pin means its other pin
is a semantic constant, which is removed first by the rank-preserving operation.
These operations strictly reduce g, so can be repeated to reach a state in which
every gate is nonconstant and genuinely uses both of its pins on the free-input
space. In particular no XOR-prefix gate then has a self-pin: its equation would
force the other pin to be constant. Repeated pins are excluded for the same
reason. Therefore a reconstruction starts with distinct pins and no self-pins,
and its new constant control has one additional distinct consumer on the reversed
path. These observations justify the consumer counts in E.

Dead gates in the acyclic part can be deleted in reverse order. A dead XOR
component can be removed as a block: the live prefix block remains invertible
because no dead node feeds a live node and A is block triangular. Removing m
gates changes at most 2m free-input predecessor occurrences, so dPhi<=2m and
Delta M>=3m. Thus all gates can be assumed to reach the selected output.
These are free simplifications with nonnegative measure change, not restrictions
of the source or new premises about the function.

## C6. Adjacent troubled gates and binary truth-table choices

For any AND-type truth table a(u,v), one output value occurs on three of the
four pairs. For any other binary truth table b(u,v), two of those three pairs
have the same b value. These two points form an affine line in F_2^2. Therefore
one of u=0, u=1, v=0, v=1, u=v, u=v+1 makes both a and b constant. This proves
all binary truth-table choices used below without an unstated finite enumeration.

For a noncompact troubled pair, write its inputs as (x,y) and (y,z), with x,y,z
distinct. If y is protected, fix it to kill the first gate. Both gates and a
descendant of the killed gate are removed by constant/unary normalization;
coincidences are handled by C3. This loses one influential input and one quadratic
equation and removes at least three gates, with dPhi<=3. Hence Delta M>=69>=64.
If y is unprotected, fix x and z to kill their respective gates. Then x,y,z all
cease to be influential: x and z are restricted, and y has lost its only uses
and is unprotected. Two restrictions give at least 3*48=144>=2*64, in addition
to nonnegative normalization gains.

For a compact pair with the same two inputs, a protected input is handled as
above. If both are unprotected, use the affine line from the truth-table argument
to make both gates constant. Both variables cease to be influential after one
restriction, giving at least 96>=64. Thus any nonempty packing gives an admissible
one- or two-restriction step. No packed component is assumed to be a formula.

## C7. Induction and conversion to gate count

Module E proves the following local assertion at the fixed integer weights above:
if the source has s>2d+2 free bits, then either a simplification strictly decreases
g without increasing M, or there is a block of t=1,2,3 admissible restrictions
with Delta M>=64*t, followed by nonnegative normalizations. C2 guarantees that
all credited operations retain a nonconstant function, and that a topmost AND
gate with affine ancestors has an actual descendant in such a state.

Begin with all n data bits free and the maximum packing used in C4. Whenever a
nonempty packing occurs use C6. Otherwise simplify as in C3 and C5 and apply the
first available case in E. A free simplification decreases the finite gate
count; a restriction decreases s. Thus there cannot be an infinite sequence
with neither progress measure decreasing. Continue until at least n-2d-2
restrictions have been applied. A final block has at most three restrictions,
so its source still has at least 2d free bits and remains nonconstant. Since
M is nonnegative at the end, summing the measure decrements gives

    M_initial>=64*(n-2d-2).

If n-2d-2<=0 the asserted measure inequality is immediate. Combining it with
(C.1) and rearranging yields

    5*g >= (31/2)*n-(261/2)*d-128,
    g >= (31/10)*n-(261/10)*d-128/5.                              (C.2)

This is a bound on the original acyclic B2 circuit: initial normalization could
only reduce its gate count, and enlarged fair semicircuits were used only as
proof states. All gates, constants, unary effects and output selection have
been accounted for. The potential is converted explicitly to a gate count; it
is never declared to be one. Applying (C.2) to the data restriction in F6 gives
the claimed full-input asymptotic coefficient.
