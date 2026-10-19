# A self-contained 7/2 lower bound for one total P family

## Exact theorem, scope and proof map

For the single total Boolean family G_N defined in F4 at every actual input
length N, including its k-bit field-description prefix and every exceptional
length, let L={z:G_|z|(z)=1}. Then L is in P and hence NP, and

    liminf_(N->infinity) C_B2(G_N)/N >= 7/2.

Here B2 is exactly the benchmark model: deterministic acyclic circuits,
one Boolean output, all sixteen two-pin Boolean truth tables, unrestricted
depth and sharing, possibly nonuniform, exact on every input. Gates are
counted once; inputs, wires and output selection are not gates. Retained
negated projections are charged, and restriction-introduced constants are
handled explicitly. The inputs of these circuits are the actual N bits,
not NP witnesses.

This is an attributed reconstruction and instantiation of the 7/2 advance of
Om Swostik Mishra and Rahul Santhanam, ECCC TR26-260 (2026-10-19). It does not
claim an original coefficient or their stronger logarithmic remainder.
The construction and its elementary supporting proofs are reused from this
participant's earlier field-prefix proposal. The earlier 31/10 circuit
case analysis is replaced in full by module C below.

All mathematical premises needed here are proved in this file. F1-F5 give
one fixed evaluator, both directions of NP membership, field validity and
all-length parameter coverage. D1-D9 prove the specialized affine-disperser
property from subspace-polynomial coefficients. C1-C8 prove the full-model
gate-elimination argument with an injective parameter map, sparse invertible
affine equations and ten complete cases. F6 applies it to a hard restriction,
charges every prefix bit and proves the eventual epsilon quantifiers.

The only scored assertion is the exact theorem above, encoded as reduced
rational strings in claims[0]. The finite diagnostic program certifies only
its explicitly bounded checks. It and the cited papers are not substitutes
for the proofs in this file. Source versions, hypotheses and contribution
provenance are separately recorded in source-obligations.md.

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
 t=ceil(sqrt(r)), using exact integer arithmetic. The defining inequality excludes k=2 for N>=1,
so every accepted field degree is odd. If t>=r or
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
NP membership and in fact P membership. For a conservative explicit time bound,
all parameters are at most N, there are at most O(N^2) field operations of
O(N^2) elementary bit cost each, and the bounded primality scan and integer
bookkeeping are polynomial as well; O(N^6) bit operations suffices for the
whole fixed evaluator. No enormous integer exponent is expanded in unary.

## F6. The hard restriction, finite bound, and every asymptotic quantifier

At every sufficiently large N, choose any irreducible q of the searched odd
prime degree k; F1 proves existence. This is only a restriction of a circuit
that must compute G_N on every prefix, not advice to the language evaluator.
The injective zero-padding map from the n=N-k data bits to K^r preserves affine
dimensions. Module D proves that this restriction is nonconstant on every
affine space of dimension greater than

    B=k*r/t+k*t+r*(r+1).

Let d=floor(B)+1. D9 therefore gives a disperser for dimension d. If d>n
the finite circuit bound below is nonpositive; eventually d=o(n), so this
issue disappears. Module C applies directly to this very n-bit restriction
and gives C_B2(f_N)>=(7/2)n-7d.

Fixing the k input bits of an N-bit Boolean circuit to q cannot increase its
gate count: substitute constants and retain or normalize the same binary
truth tables. No field gate or uncharged negation is added. Thus

    C_B2(G_N) >= (7/2)n-7d
                = (7/2)N - [(7/2)k+7d].                         (F.1)

For clarity the asymptotic loss has an explicit eventual envelope. F3 supplies
a prime in (m,4m] for every sufficiently large integer m. Taking
m=ceil(2N^(3/5)) shows that the least searched prime satisfies, eventually,

    2N^(3/5)<k<=4ceil(2N^(3/5))<=12N^(3/5).

As proved in F5, all length checks are inactive after one fixed threshold N0.
Also r<=N/k+1<=1.5N^(2/5) and t<=sqrt(r)+1<=3N^(1/5).
Using r/t<=sqrt(r), the definition of d gives

    d<=k*sqrt(r)+k*t+r^2+r+1
      <=(15+36+9/4+3/2+1)N^(4/5)
      <64N^(4/5).

These bounds hold at every integer N>=N0 after enlarging N0 if necessary;
N0 is independent of epsilon and the circuit. Consequently

    (7/2)k+7d <=42N^(3/5)+448N^(4/5)<=512N^(4/5).

In particular, (F.1) implies
C_B2(G_N)>=(7/2)N-512N^(4/5) for every N>=N0. For any epsilon>0, choose

    N_epsilon=max(N0, ceil(max(1,512/epsilon)^5)).

Then for every integer N>=N_epsilon,
C_B2(G_N)>=(7/2-epsilon)N. The same language, input encoding and evaluator
work for all these choices. This proves the declared liminf, with all actual
input bits charged. It is a finite linear coefficient, not a superlinear result.

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

Assume d>=1. Indices of a_j are read modulo k. If d=k, then Q_U=X,
so both zero-block assertions below are immediate. For d<k, no d consecutive
coefficients of Q_U are all zero. If they were, raising Q_U to an appropriate power of two and reducing
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

# C. A self-contained 7/2 circuit argument

This module reconstructs the argument of Mishra and Santhanam, ECCC TR26-260,
2026-10-19, Sections 2-4 and Theorem 4.1. All transformations and estimates needed
here are proved below. The circuit argument is independent of how the affine
disperser is evaluated. In particular it uses neither the optional Li-Yang LP
nor the more recent logarithmic-dimension disperser cited in that paper.

## C1. The original model and the larger proof state

Start with any deterministic acyclic, one-output circuit with two pins per gate
and any of the sixteen binary Boolean truth tables. Shared gates are counted
once. Repeated pins, constant and unary tables, arbitrary fan-out and arbitrary
depth are allowed. Its only inputs are the n data bits. The circuit need not
be produced uniformly.

Over F_2 every binary table has the unique expression c+a*u+b*v+d*u*v.
There are two constant, four unary, two XOR-type and eight AND-type tables.
An XOR-type table is u+v+c. An AND-type table is (u+alpha)(v+beta)+gamma.
These identities follow either by evaluating the four input pairs or by the
unique multilinear interpolation on F_2^2. Negating a wire changes constants
in its consumers' truth tables, without creating an extra input or gate.

For the proof we permit the final output to be a constant or a literal of an
existing value. This is a relaxation of the benchmark model: a retained
negated-projection gate in the original circuit still counts one, while
forgetting that cost can only make a lower bound harder. No transfer in the
opposite direction is being assumed. Thus proving the bound in these relaxed
intermediate states proves it for the original costed circuit.

A state has m independent Boolean parameters u, a set P of q disjoint unordered
pairs of those parameters, and an injective map into the original data space,

    Phi(u)=b+L*u+sum_({p,r} in P) c_(p,r)*p*r in F_2^n.             (C.1)

The vectors b,c_(p,r) have length n. A parameter in a pair is paired; all
others are unpaired. Write a=m-q. Parameters may remain in Phi even if the
circuit no longer reads them; losing all circuit uses does not delete a
parameter or reduce m.

The circuit part of a state consists of:
(i) r equations A*v=N*u+c, with r auxiliary variables v and invertible
r-by-r matrix A over F_2, each equation containing at most three distinct
nonconstant terms; and
(ii) an acyclic sequence of two-pin B2 gates reading parameters, auxiliary
values, constants or earlier gates. No gate feeds the equation system.
The output is a constant or a literal of a parameter, auxiliary value or gate.
Set s=r+(number of gates). The state invariant is

    C(u)=f(Phi(u)) for every u.                                   (C.2)

Invertibility is exactly what makes the auxiliary values unique at every
parameter assignment. Equations have no distinguished output variable; this
will allow exchanges of a parameter and an auxiliary value. The initial state
has Phi the identity, no pairs or equations, m=a=n and s=|C_original|.

The degree of a value is its number of distinct child gates and equations,
counting an equation when that value has a nonzero coefficient in it. Two pins
of one gate give one child. The designated output is not a child. Its possible
exceptions are dealt with in C3, rather than deleting an output as unused.
For unpaired parameters define w(0)=0, w(1)=6, w(j)=7 for j>=2. Put

    W=2s+3q+sum_(x unpaired) w(deg(x)).                            (C.3)

W is nonnegative. Below Delta means value_after minus value_before.

## C2. Changes of the parameters and their representation

The following four operations preserve the form (C.1).

AR: replace an unpaired x by w+phi(u_without_x), where w is fresh and phi
is affine in the other parameters. This is a bijective change, with inverse
w=x+phi. It leaves m,q,a unchanged.

QR: for distinct unpaired x,y,z, replace z by
w+(x+alpha)(y+beta)+gamma, with w fresh, and add the pair {x,y}.
This is also bijective. It leaves m unchanged, adds one to q, and reduces a
by one. Because z,x,y were unpaired, the substitution changes only a linear
term of Phi and adds just the one allowed disjoint quadratic pair.

UF: replace an unpaired x by an affine function of the remaining parameters.
This is an injective parameterization of a restriction, reducing m and a
by one and preserving q. The algorithm uses only constants and x=y+c.

PF: for a pair {p,r}, either fix p to a constant or impose p=r+c; then remove
p and its pair. This reduces m and q by one, leaving a unchanged. The old
quadratic term becomes affine: p*r=c*r in the first case and
p*r=(r+c)*r=(1+c)*r in the second. The other pairs stay disjoint.

In each case the new map is Phi composed with the indicated injective
parameter map. Thus it stays injective. In AR and QR its image stays the
same; in UF and PF its image is restricted. Gates/equations will be changed
only in the ways justified next so that (C.2) also remains true.

An affine gate whose pins are parameters, constants or auxiliary values can
be moved into the system: replace J=b1*t1+b2*t2+c by its equation, making
J auxiliary. The new auxiliary matrix is block triangular with diagonal
blocks A and 1, hence invertible. There are at most three terms, s is unchanged,
and no parameter obtains an extra child. Repeat in topological order, moving
every affine gate with no nonlinear ancestor. This never increases W.

Here is the implementable AR exchange. Suppose an auxiliary value Q has solved
form Q=x+phi(u_without_x), with x unpaired. Leave every equation and gate
unchanged, but regard Q as fresh parameter w and x as auxiliary; make the AR
change only in Phi. If Q is auxiliary coordinate j and N_x the old parameter
column of x, put t=A^(-1)N_x. Its j-th coordinate is 1. The new auxiliary
matrix is A with column j replaced by N_x, namely A*B, where B is the identity
with column j replaced by t. Expansion along that column gives det B=t_j=1.
Fairness and the three-term bound follow. Only x's parameter weight is replaced
by Q's. In particular if both have degree at least two, Delta W=0.
Every auxiliary column is nonzero, so Q occurs in at least one equation.

Fixing a parameter to a constant or a literal of another parameter does not
change A and cannot increase any row's number of terms. It therefore preserves
fairness and the row bound without adding gates.

## C3. Source survival and output exceptions

Fix one member of each pair in (C.1) to a constant. The remaining map is affine
and injective on m-q parameters; its image is an affine subspace of dimension a.
Thus if f is a disperser for dimension d, f composed with Phi is nonconstant
whenever a>=d.

When a>d its value cannot be affine in the parameters. If a nonconstant affine
output contains an unpaired x, one UF equation makes that output constant,
leaving affine dimension at least a-1>=d, a contradiction. If every occurring
parameter is paired, first fix the partner of one occurring parameter by PF;
a does not change and the chosen parameter becomes unpaired, so the same
argument applies. A constant output was already excluded by survival.
This use of a general affine UF is a statement about the source alone and
does not require an implementation in the circuit.

Similarly the output cannot be an AND-type gate on two unpaired parameters:
fixing one input to its killing value makes it constant at cost one in a,
again leaving dimension at least d. Consequently an output parameter or
auxiliary literal, and an output minimal AND on unpaired parameters, cannot
be mistakenly deleted as an ordinary unused value during a live iteration.

## C4. Free deletions, including all rank and weight accounts

Rerouting means substituting a constant or literal of an existing value for
every use of a value that has become equal to it. For an acyclic gate there
is no effect on the affine system. For a system equation v=U+c, with U either
an auxiliary value, a parameter or zero, add that row to every other row
containing v, substitute in gate inputs and the output, then delete its row
and auxiliary column v. The row operations preserve determinant. Column v
now has its only 1 in the deleted row, so cofactor expansion shows that the
remaining auxiliary matrix is invertible. Replacing a term by a constant
or one literal cannot increase row size.

We shall use the following deletions, each with Delta a=Delta m=Delta q=0
and Delta s=-1.

D1. A gate that is a constant or a literal of one pin is rerouted and removed.
This includes constant or unary tables, a constant pin, repeated pins, and an
AND on affine inputs with dependent linear parts. In the latter situation
both inputs are affine functions of a common single bit (or constants), so
their binary table is constant or a literal of that bit, hence of a nonconstant
pin if one exists. It needs no new gate. At most one initially isolated
unpaired parameter can gain children: the one to which the gate is rerouted.
Its weight rises by at most one, from 6 to 7. A used parameter cannot become
newly used when it previously had degree zero. Thus Delta W<=-1. Rerouting to
a constant gives the sharper Delta W<=-2.

D2. A gate with no children that is not the selected output may be removed.
No parameter gains children, so Delta W<=-2.

D3. An equation with at most two nonconstant terms contains an auxiliary
variable, since otherwise A has a zero row. Solve it as v=U+c and use the
rank-preserving deletion above. At most one initially isolated unpaired
parameter, the literal U when it is a parameter, can increase its weight
from 6 to 7. Thus Delta W<=-1, or <=-2 if none increases.

D4. Let an unpaired degree-one parameter x have its sole occurrence in row E.
The vector N_x has just one nonzero entry, so A^(-1)N_x is nonzero. Choose an
auxiliary coordinate Q with coefficient 1 and exchange x,Q by AR. The old x
is now an auxiliary occurring only in row E and nowhere in the circuit output.
Delete x and E. Cofactor expansion along its singleton auxiliary column again
proves fairness; all other equations and gates are unchanged. Weight 6 is
replaced by at most 7, so Delta W<=-2-6+7=-1. Injectivity and the image of
Phi were preserved by the AR exchange; no source dimension was lost.

All four operations preserve s's interpretation as equations plus gates,
rather than equating a potential with a gate count. Their lower-bound credits
are valid with sharing: consumers are rerouted in place, never duplicated.

For the case counts below, process only the specifically named gates/rows.
A named gate that has already disappeared by rerouting is already counted as
a deletion; it is never counted a second time. An undeleted named equation
still contains an auxiliary variable by invertibility. Constant propagation
through other gates may be deferred to the next D1-D4 iteration. Thus each
distinct named gate or row disappears and contributes a distinct unit to
the decrease of s, even when consumers or descendant subgraphs coincide.
The explicit distinctness conditions are recorded in the individual cases.

## C5. Six preliminary cases

At the start of each iteration move eligible affine gates into the system.
Apply the first available case in the order below. In every case, Delta a is
0,-1 or -2 and Delta W<=14 Delta a; also m+s strictly decreases.

Case 1. A D1-D4 deletion is possible. Apply it. This decreases s, preserves a,
and decreases W. If none applies, every remaining gate has two distinct
nonconstant pins and is XOR- or AND-type; equations have exactly three terms;
every nonoutput gate has a child; and no unpaired degree-one parameter occurs
in the system. An AND on affine inputs has independent linear parts.

Case 2. An AND G on parameters z,w has z unpaired of degree one.
Fix w to kill G, delete G and reroute to its constant. Then z becomes unused
and its weight falls from 6 to zero. No other existing unpaired parameter gains
a child. If w is unpaired, its weight was at least 6 and UF gives
Delta a=-1, Delta W<=-6-6-2=-14. If w is paired with r, PF removes the pair,
r's new weight is at most 7, Delta a=0 and Delta W<=-3-6-2+7=-4.

If Cases 1-2 fail, no row or gate has both an isolated unpaired parameter and
a different parameter as parents. For rows this follows from D4. An affine
gate on two parameters would have been moved into the system; an AND on such
parents is Case 2. An isolated unpaired parameter cannot be the output by C3.
We will repeatedly use this no-isolated-parent fact to bound weight increases.

Case 3. A paired p with partner r has degree zero or at least two.
If p is unused, fix r by PF. Its pair disappears and p becomes unpaired but
still unused. No deletion is needed: Delta a=0, Delta W<=-3, Delta m=-1.
If p has at least two distinct children, fix p by PF and delete two of those
children. Each is a gate with a constant pin or a row of at most two terms.
The no-isolated-parent fact prevents any initially unpaired weight from
increasing. Only newly unpaired r can acquire weight, at most 7. Consequently
Delta a=0, Delta s<=-2 and Delta W<=-3-4+7=0.
Thus after this case is excluded every paired parameter has degree exactly one.

Case 4. Two partners p,r have the same unique child H.
Impose p=r+c by PF to make H constant. If H is a row it is v+p+r=c0,
so any c works; the same holds for XOR. For AND
(p+alpha)(r+beta)+gamma choose c=alpha+beta+1, making its product zero.
Delete H and reroute to the constant. Surviving r is unused, so contributes
zero after becoming unpaired. Delta a=0, Delta W<=-3-2=-5.

Case 5. Paired p occurs in the system, with partner r.
Their unique children H_p,H_r are distinct by Case 4. Fix r by PF and delete
H_r. It contains no initially isolated unpaired parameter, by the
no-isolated-parent fact. The new unpaired p has weight at most 6. Removing
the pair and one child changes W by at most -3-2+6=1.
If p is still isolated in row H_p, use D4, gaining at least one more unit.
If it has already become unused because its row was deleted while processing
H_r, both named rows/gates have already disappeared and p's weight is zero,
which is stronger. Overall Delta m=-1, Delta s<=-2, Delta a=0, Delta W<=0.

Case 6. An AND G on parameters p,w has p paired with r.
Both p,r have unique children, G and H_r respectively. Case 4 excludes
H_r=G and w=r; distinct pins exclude w=p. Case 5 implies H_r is a gate.
Case 2 implies that w, if unpaired, has degree at least two.
Fix r by PF; then fix w so G is killed, using UF or PF as appropriate.
Delete the distinct gates H_r and G. The surviving p becomes unused.
No initially unpaired weight increases: G is rerouted to a constant, and
the other input of H_r cannot be an initially isolated unpaired parameter.
If w was unpaired, Delta m=-2, Delta q=-1, Delta a=-1, and
Delta W<=-7-3-4=-14.
If w was paired with r', disjointness gives two different removed pairs;
only r' can acquire positive new weight, at most 7. Now Delta m=-2,
Delta q=-2, Delta a=0, and Delta W<=-6-4+7=-3.

The distinctness and weight arguments in Cases 3,5,6 use direct children,
not an assumption that all downstream gates are different. Additional
normalizations may be postponed and then charged by D1-D4.

## C6. Preparing the first nonlinear gate

Suppose a>d and no preliminary case applies. If the acyclic part were empty,
the output would be affine, contradicting C3. After moving all eligible
affine gates, its first gate G is AND-type, with parameters or auxiliary
values as pins. Case 1 ensures their two linear parts are independent.

If a pin Q is auxiliary, choose a parameter x with coefficient 1 in Q's
solved affine expression. It occurs in N, hence in some system row.
Cases 1 and 5 make x unpaired of degree at least two. The auxiliary Q occurs
in a row and in G, so has degree at least two. Exchange x,Q by AR. All
of m,s,q,a,W stay unchanged. Move eligible affine gates and check the
preliminary cases again, ending this iteration if any reduction applies.

If the other pin is auxiliary, it depends on some parameter other than the
first new pin: otherwise the two linear parts would be dependent. Make a
second exchange. If no preliminary reduction intervenes, G now has two
distinct unpaired parameters x,y with degrees at least two. Independent
linear parts stay independent under the invertible affine coordinate change.

For any child J of G that also reads x, replace its operation h(x,G(x,y))
by that two-variable truth table on x,y. This removes J as a child of G,
retains it as a child of x, and may add a child to y. The latter already
has weight 7, so W does not increase. Do the symmetric step for y as needed.
After each step move eligible affine gates and recheck preliminary cases.
If a deletion/restriction occurs, end the current iteration with its established
bounds. Otherwise each rerouting removes one child of this fixed G, and none
of the other preparatory operations creates one; this finite process stops.

We may therefore suppose G has inputs x,y, both unpaired with degrees at
least two, and none of G's children reads x or y. By C3 G is not the selected
output; by Case 1 it has k>=1 distinct children. Swap the names x,y if necessary
so deg(x)>=deg(y). We now use exactly one of Cases 7-10.

## C7. The four final cases

Case 7. Both degrees are two and the second child B of x,y is the same.
If B is an equation, it is v+x+y=c. Otherwise it is a two-pin table h(x,y).
An AND-type table G has one value on three of the four points of F_2^2.
Two of these three points have equal h values (use h=x+y for the row case).
Those two points form an affine line. Its equation is a fixing of x, a fixing
of y, or x=y+c, all permitted UF restrictions. On it G and B are constant
(or B's auxiliary v is constant). Delete the two distinct objects and reroute
them to constants. The removed parameter and the surviving now-unused one
both lose weight 7. Hence Delta a=-1, Delta s<=-2, Delta W<=-14-4=-18.
If Case 7 is absent and both degrees are two, write their different second
children as B_x and B_y.

Case 8. A general killing restriction outside the bottleneck described below.
Let h=deg(x)>=deg(y)>=2, and let k_1<=k count G's children whose other pin is
an isolated unpaired parameter. Fix x to kill G. Delete all h children of x
and all k children of G. These h+k original objects are distinct: no child of
G reads x. Rows and gates touching x become deletable; after G is replaced
by its constant, each child of G has a constant pin and is deletable.

Parameter x loses weight 7. An initially unused parameter cannot gain a child.
An initially degree-at-least-two unpaired parameter already has maximum weight.
An initially isolated unpaired parameter can increase its weight by one only
when its unique child is deleted and rerouted to a literal of it. By the
no-isolated-parent fact that child is not a child of x. It is therefore among
the k_1 indicated children of G. Even with sharing there are at most k_1
such weight increases, each at most one.

If h=2, then deg(y)=2 as well. It has no occurrence in B_x because Case 7
is excluded, and no occurrence in a child of G by preparation. Rerouting G
to a constant removes one of its two uses and creates none, so y's weight
drops by at least one. Thus, with Delta a=-1,

    Delta W <= -7-2h-2k+k_1                 if h>=3,
    Delta W <= -12-2k+k_1                  if h=2.                (C.4)

For h>=3 the first expression is <=-7-6-k<=-14. For h=2,k>=2,
the second is <=-12-k<=-14. For h=2,k=1,k_1=0 it is exactly <=-14.
The only remaining configuration has both degrees two, a unique child D of G,
and D's other pin an unpaired isolated z. This is the bottleneck.
Here z differs from x,y, and B_x,B_y are distinct and different from D.

Case 9. The bottleneck with D AND-type.
Write D=(G+delta)(z+eta)+theta. Since G is a nonconstant AND table,
choose x0,y0 with G(x0,y0)=delta. Fix x=x0,y=y0 by two UF restrictions.
G and D become constants. Delete G,D,B_x,B_y, four distinct original objects.
The second children are distinct from each other by exclusion of Case 7,
and distinct from D by preparation; each is deletable after its parameter
is fixed. Parameter z loses its only child. No initially unpaired parameter
increases its weight: G,D are rerouted to constants, and B_x,B_y cannot
have another isolated unpaired parameter as parent. The weight losses
of x,y,z are 7,7,6. Thus Delta a=-2 and
Delta W<=-7-7-6-8=-28=14 Delta a.

Case 10. The bottleneck with D XOR-type, D=G+z+c.
Apply QR to z=w+G(x,y)+c and pair x,y. Its change is bijective on parameters;
since z occurred only at D, it affects no system equation or other gate.
D becomes w and is rerouted to w, then removed. G is now unused and not the
output, so is removed too. No original unpaired parameter obtains a new
child: only the fresh w can acquire D's former children. The old weights
7,7,6 for x,y,z disappear; the new pair contributes 3 and w at most 7.
Consequently Delta m=0, Delta q=1, Delta a=-1, Delta s=-2, and

    Delta W<=-4-7-7-6+3+7=-14=14 Delta a.                         (C.5)

This is where the new proof improves on a restriction-only treatment of the
bottleneck. It changes the coordinates bijectively while charging the new pair.
It does not add a free AND gate or pretend that unrestricted quadratic
substitutions preserve the disperser.

## C8. Coverage, termination and the finite circuit bound

The preliminary cases are tested in order; all uses of their failure above
are therefore justified. If they fail, C6 reaches a prepared minimal AND,
unless an already proved preliminary reduction intervenes. Cases 7-10
cover its degrees, child multiplicity, and all ten genuine binary gate
tables at its bottleneck. Constant, unary and repeated-pin cases were
already D1. An affine-only output or a terminal minimal AND is excluded
by C3. Thus whenever a>d an iteration exists satisfying

    Delta a in {0,-1,-2}, Delta W<=14 Delta a, and Delta(m+s)<0.   (C.6)

Preparatory moves and exchanges preserve m,s,a and do not increase W.
Each moving step is topological, at most two exchanges are needed for the
chosen G, and every child-rerouting removes a child of that fixed gate.
Hence preparation terminates. The subsequent preliminary or final reduction
strictly decreases m+s: D1-D4 decrease s, QR decreases s, and every other
case reduces m without adding any gates/equations. There cannot be an
infinite sequence of iterations. All parameters and equations remain finite.

If 1<=d<n, run while a>d, starting at a=n. The process must reach a_T<=d,
since otherwise another iteration would exist and decrease the nonnegative
integer m+s. The last iteration may reduce a by two; this is harmless.
Summing (C.6) and using W_T>=0 gives

    2|C_original|+7n >= W_0
                        >= W_T+14(n-a_T)
                        >= 14(n-d).

Therefore every original full-B2 circuit satisfies

    |C_original| >= (7/2)n-7d.                                   (C.7)

If d>=n this lower bound is nonpositive and follows from gate-count
nonnegativity. No step used uniformity of the circuits. Equations and
reparameterizations were only proof states with explicitly charged size;
the source language is still evaluated by the fixed algorithm in F.
