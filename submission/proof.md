# A complete 4n - 65 n^(4/5) circuit lower bound for one total language in P

## Part I. The family, uniformity, model and complete reduction

### I.1. Statement and provenance

We define a single total Boolean family F_n below, at the actual input bit
length n, and prove

    C_B2(F_n) >= 4n - 65 n^(4/5)    for every integer n >= 16807.

Here B2 is the full basis of all sixteen binary truth tables; circuits are
deterministic, acyclic, single-output, and may have arbitrary depth, fan-out,
shared descendants and repeated pins. Each retained gate is counted once.
The circuits may be nonuniform. Thus the declared eventual guarantee holds,
and liminf over all integer n of C_B2(F_n)/n is at least 4. The coefficient
and exponent are the coordinates (4,1) of this one joint theorem.

The evaluator for F_n is one fixed deterministic polynomial-time algorithm.
In particular its yes-language is in P and hence in NP, with a verifier that
requires an empty witness. All lengths below 1024, including zero, have
F_n=0; no hard subsequence or length-dependent advice is used.

This is an application and reconstruction of published methods, not a claim
of an original circuit lower-bound record. Part II proves the needed
Ben-Sasson--Kopparty (2010) affine-disperser construction. Part III proves the
Mishra--Santhanam (ECCC TR26-260, revision 1, 3 November 2026) finite circuit
bound 4m-10d. We combine them using an input prefix that specifies the field
modulus. The prefix removes the need to assume an efficient algorithm for
constructing an irreducible polynomial, and its full bit cost is charged.
The complete arguments, rather than only their theorem citations, appear
in this document. No finite search, LP certificate, prior Yukon judgment,
or logarithmic-dimension disperser is used as a premise.

### I.2. A fixed total evaluator

Given x in {0,1}^n, perform the following algorithm.

1. If n<1024, output 0.
2. Let k be the least prime strictly larger than 2 n^(3/5). This is an
   exact integer prescription: scan integers k in increasing order,
   use trial division to test primality, and test k^5>32 n^3.
   Put m=n-k, r=ceil(m/k), and t=ceil(sqrt(r)).
3. Interpret the first k bits b_0,...,b_(k-1) as the lower coefficients
   of the monic polynomial

       q(T)=T^k + sum_{j=0}^{k-1} b_j T^j in F_2[T].

   Its leading coefficient is implicit; the prefix has exactly k bits.
   Test both

       T^(2^k) = T modulo q(T),    gcd(q(T), T^2-T)=1.

   If either test fails, output 0. The first test is implemented by k
   repeated squarings modulo q, never by writing a degree-2^k polynomial.
4. If both tests pass, form K=F_2[T]/(q), using the polynomial basis
   1,T,...,T^(k-1). Append rk-m zero bits to the remaining m data bits
   and divide them into r blocks of k bits. Interpret a block in K with
   its first bit as the constant coefficient. Denote the resulting field
   elements by z_1,...,z_r. This map phi:F_2^m -> K^r is linear and injective.
5. For i=1,...,r, set a_i=2^(r+2-i)-1 and y_i=z_i^(a_i). Compute

       E_t = sum_{I subset {1,...,r}, |I|=t} product_{i in I} y_i.

   Use the following polynomial-time recurrence, not enumeration of
   subsets: initially E_0=1 and E_j=0 for 1<=j<=t; for each new y_i,
   update E_j <- E_j+y_i E_(j-1) for j=min(i,t),...,1 in descending order.
6. Output Tr(E_t)=sum_{j=0}^{k-1} E_t^(2^j), regarded as a bit in F_2.

The parameters are valid whenever n>=1024: Part II.6 proves that k is an
odd prime with 0<k<n, r>=3, 1<=t<r and k>=r(r+1)/t. Thus there is no
undefined block size or symmetric polynomial in this branch. Invalid
polynomial prefixes have the already specified output 0.

The algorithm fixes the family before any circuit, accuracy parameter or
asymptotic quantifier is chosen. It makes no choice of a field modulus:
that modulus is supplied in the actual input and is tested deterministically.

### I.3. Why the field test is correct and a valid prefix exists

We give the required elementary finite-field facts to make both the
evaluator and the lower-bound restriction precise.

Let h be any monic irreducible binary polynomial of degree j. Its quotient
field E=F_2[T]/(h) has exactly 2^j elements. Lagrange's theorem in its
finite multiplicative group shows z^(2^j)=z for every z in E. Frobenius
sigma:z->z^2 is an automorphism of E. For alpha=T mod h, its orbit length
under sigma is exactly j. Indeed the length l divides j. The polynomial
product_{i=0}^{l-1}(X-sigma^i(alpha)) is fixed coefficientwise by sigma;
the fixed elements of sigma are precisely 0 and 1, because they are the
roots of X^2-X. This polynomial is in F_2[X], vanishes at alpha, and so
is divisible by its degree-j minimal polynomial h. Therefore l>=j,
whence l=j. It follows that

    h divides T^(2^k)-T  if and only if  j divides k.           (I1)

The polynomial T^(2^k)-T is squarefree: its formal derivative is -1=1.
As k is prime, all of its irreducible factors have degree 1 or k by (I1).
Consequently, a monic degree-k q satisfying the first test is squarefree
and has only factors of those degrees. The gcd test removes its possible
degree-one factors, since the only monic irreducible linear polynomials
over F_2 are T and T+1. What remains must be a single degree-k irreducible
factor. Conversely every irreducible degree-k polynomial passes both
tests. This proves their equivalence to irreducibility for this algorithm.

Existence requires no search procedure. Factor T^(2^k)-T into monic
irreducibles in F_2[T], whose Euclidean division and unique factorization
are elementary. Its only possible linear factors are T and T+1, each
occurring once. Its degree is 2^k>2. By (I1) at least one remaining
factor has degree k. Thus there is an irreducible degree-k polynomial
and hence at least one valid k-bit prefix at every relevant input length.

For any valid prefix, K really is a field of size 2^k. Frobenius to the
k-th power is the identity, so the trace in the evaluator is fixed by
Frobenius and belongs to F_2. Also Tr(1)=1 because k is odd. The output
therefore is a well-defined bit and the trace is a nonzero linear map.

### I.4. Uniform polynomial time and NP membership

Bertrand's postulate gives k<=4 n^(3/5), as proved with the exact
rounding in Part II.6. The search by integer arithmetic and trial division
is thus polynomial in the input length n. The field-validity tests use
polynomial gcd and k repeated modular squarings on polynomials of degree
less than k; elementary binary polynomial arithmetic implements them in
polynomial time, even for prefixes that ultimately fail the test.

For a valid prefix, addition in K is bitwise XOR and multiplication is
ordinary polynomial multiplication followed by division by q. Each costs
polynomial time in k. The exponent a_i has r+2-i<=r+1 binary digits,
so square-and-multiply evaluates each y_i in O(r) field multiplications.
The symmetric-polynomial recurrence takes O(rt) field operations. Its
correctness follows by partitioning j-subsets of {1,...,i} into those
containing i and those not containing i; descending updates preserve
the previous row of the recurrence. Trace evaluation takes k squarings
and additions. Since k,r,t<=n, these are all bounded by one fixed
polynomial in n. Field elements are k-bit strings throughout, and there
is no exponential-size object or truth-table lookup in the evaluator.

Define the total language L={x:F_|x|(x)=1}. Let A be exactly the algorithm
above and define V(x,w)=1 precisely when w is empty and A(x)=1. Set the
polynomial witness bound W(n)=0. If x is in L, the empty witness is
accepted. If any allowed witness is accepted, it is empty and A(x)=1,
so x is in L. This proves both directions of the NP specification with
one uniform polynomial-time verifier. Its witnesses are not circuit inputs.

### I.5. All-length lower bound and actual bit accounting

Fix any integer n>=16807=7^5 and the length parameters from I.2. Choose
one existing irreducible q of degree k using I.3. This choice is only for
analyzing a circuit restriction, not advice supplied to the evaluator.
On the slice whose first k input bits specify this q, the data-only
function h:F_2^m -> F_2 is exactly B composed with phi from Part II.

Put

    D=kr/t+kt+r(r+1),    d=floor(D)+1.

The exact algebraic theorem in Part II proves that h is nonconstant on
every dimension-d affine subspace. Part II.6 gives

    k<=4 n^(3/5),       d<=6 n^(4/5)+1.

These bounds imply 1<=d<=m for n>=7^5. To see the upper bound, write
a=n^(1/5)>=7. Then m>=a^5-4a^3, and

    a^5-4a^3-(6a^4+1)=a^3(a^2-6a-4)-1>0.

Here a^2-6a-4>=3 for a>=7. Thus Theorem C in Part III applies and
every circuit for this slice has at least 4m-10d gates.

Take any full B2 circuit for F_n, of size S. Fix its k prefix input
bits to the coefficients of q. The resulting circuit computes h.
It has no more than S gates in the relaxed convention with free
constant sources and complemented wires used in Theorem C. Substitution
does not duplicate any gate or subgraph. Theorem C proves its lower
bound even in that more permissive model, so

    S >= 4m-10d
      >= 4n - 60 n^(4/5) - 16 n^(3/5) - 10
      >= 4n - 65 n^(4/5).                                  (I2)

For the last inequality, a>=4 gives 16a^3<=4a^4 and 10<=a^4, so
16a^3+10<=5a^4. Every prefix bit is accounted for through m=n-k;
the field blocks and zero padding introduce no extra circuit inputs.
Neither a direct-sum theorem nor an assumption about disjoint circuitry
is used. The same argument works for every integer n>=16807.

Since S was arbitrary, (I2) proves the claimed eventual lower bound.
For every epsilon>0, choosing an integer threshold at least 16807
and at least (65/epsilon)^5 gives C_B2(F_n)>=(4-epsilon)n at all
larger integer n. This is the declared joint finite projection (4,1).
It makes no assertion of superlinear growth or an endpoint score.

### I.6. What remains to prove in this document

Part II below supplies the entire affine-disperser argument invoked in
I.5, including the coefficient lemmas and the parameter bounds. Part III
supplies the entire finite circuit argument, including its normalizations,
reparameterizations, exceptional output cases, counting and termination.
These parts are proofs, not hypotheses or computational certificates.
Bertrand's postulate and elementary field, polynomial and linear algebra
are the only standard background facts used without a separate proof.

## Part II. A complete binary Ben-Sasson--Kopparty disperser lemma

This part supplies the algebraic premise used in Part I. The construction
and the coefficient argument are from Eli Ben-Sasson and Swastik Kopparty,
*Affine Dispersers from Subspace Polynomials*, author preprint dated
12 March 2010, Sections 3, 6 and 7, especially Theorems 6.1 and 7.3:
<https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf>.
The argument below specializes to characteristic two and prime extension
degree, proves all the needed subspace-polynomial statements, and uses a
short weighted-counting proof for the final triangular substitution step.
It does not require the reader to retrieve that source. No new disperser
theorem is claimed.

### 1. The exact algebraic statement

Let k be an odd prime, let K be any field with 2^k elements, and let r,t be
integers satisfying

    1 <= t < r,                 k >= r(r+1)/t.

For 1<=i<=r put

    e_i = r+1-i,    w_i = e_i+1 = r+2-i,    a_i = 2^(e_i+1)-1.

Write Tr(z)=sum_{s=0}^{k-1} z^(2^s), which is an F_2-linear map K -> F_2.
Define

    B(z_1,...,z_r)
      = Tr( sum_{I subset {1,...,r}, |I|=t} product_{i in I} z_i^(a_i) ).

**Disperser theorem.** On every affine F_2-subspace A of K^r with

    dim_F2 A > D := kr/t + kt + r(r+1),

the function B takes both values 0 and 1. In particular, for every injective
F_2-linear map phi:F_2^m -> K^r, B composed with phi is nonconstant on every
affine subspace of F_2^m whose dimension is strictly larger than D.

An affine subspace means a translate of a linear subspace. Its dimension
is always over F_2. The injective-map assertion follows from the first
assertion because an injective linear map preserves affine dimension.

### 2. Polynomial functions and binary degree

Put q=2^k. A polynomial in s variables over K has a unique representative
whose degree in each variable is less than q, modulo the relations
X_j^q=X_j. These relations preserve its values on K^s. Uniqueness follows
by induction on s: a nonzero univariate polynomial of degree below q has
fewer than q roots; apply this to one variable at a time. Consequently a
reduced polynomial having a nonzero coefficient on a nonconstant monomial
represents a nonconstant function.

For a nonnegative integer b let wt(b) be its binary Hamming weight. Call
wt(b) the binary degree of X^b. Reducing a positive exponent modulo
X^q=X gives the unique exponent in {1,...,q-1} congruent to the original
exponent modulo q-1. Binary degree is subadditive under multiplication of monomials and
does not increase under reduction or Frobenius powers. Here is a direct justification
covering wraparound as well as ordinary carries. Express an exponent as
a sum of powers 2^j, fold indices j modulo k, and whenever two copies
occur at one index replace them by one copy at the next index modulo k.
Every such carry reduces the number of summands by one. The process
terminates with a nonempty subset of the k indices for a positive
exponent. Its value lies in {1,...,q-1}, is congruent to the original
exponent modulo q-1, and has at most as many binary digits. For exponent
zero the claim is immediate. Frobenius powers just rotate the indices.

In particular, a product of L monomials each of binary degree one has
reduced binary degree at most L. If L<k and the final binary degree is L,
there were no carries: the chosen indices must be distinct.

Every F_2-linear map K -> K is uniquely represented by a linearized
polynomial

    Q(X) = sum_{j=0}^{k-1} a_j X^(2^j),       a_j in K.

Indeed these polynomials define linear maps, their coefficient space and
End_F2(K) both have F_2-dimension k^2, and evaluation is injective by the
univariate root bound, since every nonzero such polynomial has degree
at most 2^(k-1)<q. An affine linear map adds a constant term.

Tr has values in F_2 because Tr(z)^2=Tr(z), and is linear. Since k is odd,
Tr(1)=1, so it is nonzero.

### 3. Image polynomials of linear subspaces

For a linear subspace L of K of dimension d define

    P_L(X) = product_{u in L} (X-u).

This polynomial is monic, linearized, and of degree 2^d. To prove the
linearized assertion, induct on d. For L={0}, it is X. If
L=L_0 direct-sum F_2 v, then

    P_L(X) = P_L0(X) P_L0(X+v)
           = P_L0(X)^2 + P_L0(v) P_L0(X),

which is linearized by the induction hypothesis. All roots are distinct,
so its coefficient of X, which is its formal derivative, is nonzero.

The image L'=P_L(K) has dimension k-d. Define Q_L=P_L'. Then

    Q_L(P_L(X)) = X^q-X.

Both sides are monic of degree q and have all q elements of K as roots.
The reverse composition identity also holds. To see this, apply P_L
to the displayed identity. Its coefficients are in K and its monomials
have powers of two, so P_L(X^q-X)=P_L(X)^q-P_L(X). Thus the polynomial
P_L(Q_L(Y))-(Y^q-Y), after substitution Y=P_L(X), is zero. Composition
with a nonconstant polynomial is injective, by degrees, so

    P_L(Q_L(X)) = X^q-X.

It follows that Q_L(K) is contained in L. Q_L has kernel L', of dimension
k-d, so its image has dimension d and is exactly L. We have therefore
constructed a monic linearized image polynomial

    Q_L(X) = sum_{j=0}^{k-1} a_j X^(2^j)

with a_0!=0, a_(k-d)=1, and a_j=0 for j>k-d. Its kernel has
2^(k-d) distinct elements. The d=0 case also exists as a linear map
(the zero map); image-polynomial statements below are used only for d>1.

#### The zero-gap property

For d>1 the only cyclic block of d-1 consecutive zero coefficients of
Q_L is the displayed terminal block k-d+1,...,k-1. In particular, every
zero run between two successive nonzero coefficients in the interval
[0,k-d] has length at most d-2.

First, there cannot be d consecutive cyclic zero coefficients. A cyclic
rotation of such a block to the last d positions is obtained by taking
a Frobenius power of Q_L and reducing modulo X^q-X. It would give a
nonzero polynomial of degree at most 2^(k-d-1) with the same
2^(k-d) roots as Q_L, contradicting the root bound.

Suppose instead a different block of d-1 consecutive zero coefficients
exists. Rotate it to the last d-1 positions using Frobenius power 2^s,
where 1<=s<k. The reduced polynomial Q' is nonzero, has degree at most
2^(k-d), and vanishes at every root of Q_L. Hence Q'=c Q_L for some
nonzero c in K. On the image L this says

    u^(2^s) = c u       for every u in L.

Choose a nonzero v in L. Every u/v is fixed by the s-th power of the
Frobenius automorphism z -> z^2. Frobenius has order k: it has k-th
power equal to the identity, and for 0<j<k the polynomial X^(2^j)-X
cannot vanish on all q elements. Since k is prime, s and k are coprime;
an element fixed by the s-th power is therefore fixed by Frobenius
itself (Bezout's identity, applied to powers of this automorphism).
Its fixed field is F_2, the roots of X^2-X. Thus L is contained in
v F_2 and has dimension at most one, a contradiction. This proves the
zero-gap property, including cyclic blocks.

### 4. A product of independent subspace polynomials

We need the following special case of the independent-source argument.

**Product lemma.** Let L_1,...,L_s be linear subspaces of K, of dimensions
d_i>1, and let Q_i be their image polynomials. Let integers e_i satisfy
1<=e_i<d_i-1. If

    sum_i (d_i-e_i-1) >= k,

then the reduced polynomial

    R(X_1,...,X_s)
      = Tr( product_i Q_i(X_i)^(2^(e_i+1)-1) )

has a nonzero monomial product_i X_i^(beta_i) with
wt(beta_i)=e_i+1 for every i.

Write Q_i(X)=sum_j a_(i,j) X^(2^j), with indices j taken modulo k.
Define a k by k matrix

    M_i[u,v] = a_(i,v-u)^(2^u).

This is the coefficient of X^(2^v) in the reduced polynomial
Q_i(X)^(2^u). For A_i={0,...,e_i}, and any set B_i of e_i+1 distinct
indices, put beta_i=sum_{b in B_i}2^b. Expanding the trace gives the
following coefficient of product_i X_i^(beta_i):

    sum_{h=0}^{k-1}
       ( product_i perm(M_i[A_i,B_i-h]) )^(2^h).              (1)

Here perm is the sum over bijections between the indicated row and
column sets, without signs. To verify (1), expand each power as the
product of e_i+1 Frobenius powers of Q_i. A term with final binary
degree e_i+1 must select distinct bit indices, by Section 2, so its
choices are precisely those bijections. The h-th trace summand rotates
all indices by h and raises coefficients to 2^h, giving (1).
The inequality e_i+1<d_i<=k ensures that this argument uses fewer
than k factors for each variable.

We construct B_i so the h=0 summand is nonzero and every other summand
has a zero factor. Reorder the pairs (Q_i,e_i) so d_1=max_i d_i.
Put L_i=d_i-e_i-1>=1. Let s_0 be the first index with

    sum_{i=1}^{s_0} L_i >= k-d_1+1.

It exists by the hypothesis. For i<=s_0 put ell_i=sum_{j<i}L_j.
Then ell_1=0 and

    0 <= ell_i <= k-d_1 <= k-d_i.

Moreover the prefix ending at s_0 is at most k-2, because its
preceding prefix is at most k-d_1 and L_(s_0)<=d_1-2.
For i>s_0 set ell_i=0 instead; such additional factors will not be
needed to kill any trace summand.

For each i use the following ordered list of columns; writing a_j
for a_(i,j) temporarily makes the notation shorter.

* If a_(ell_i)!=0, take

      B_i = (ell_i, k-d_i+1, ..., k-d_i+e_i).

* If a_(ell_i)=0, let b<ell_i<j be the adjacent nonzero coefficient
  indices surrounding ell_i. They exist in [0,k-d_i] because both
  endpoints of that interval have nonzero coefficients. Put
  z=j-b-1 and v=min(z,e_i), and take

      B_i = (j, b+1, ..., b+v,
             k-d_i+v+1, ..., k-d_i+e_i).

  The zero-gap property gives z<=d_i-2.

The entries of each list are distinct. In the first case its first
entry is in [0,k-d_i] and the others are larger. In the second case
j and b+1,...,b+v are distinct and at most k-d_i, while the remaining
entries are larger.

Order rows as 0,...,e_i and columns in the listed order. The matrix
M_i[A_i,B_i] is lower triangular with nonzero diagonal. In the first
case the first diagonal entry comes from a_(ell_i), and all subsequent
diagonal entries come from a_(k-d_i)=1. Entries above the diagonal
come from the terminal zero block. In the second case the first
diagonal entry comes from a_j, the next v from a_b, and the rest
from a_(k-d_i). Above-diagonal entries in the initial columns use
indices in the zero interval [b+1,b+v]; those in the final columns
use the terminal zero block. Frobenius powers preserve zero and
nonzero status. Thus every permanent at h=0 in (1) is nonzero.

For a column b_0 in B_i, the shifted column b_0-h of M_i[A_i,B_i-h]
is zero whenever

    h in [b_0+1, b_0+L_i] modulo k.                         (2)

Indeed its entries have coefficient indices b_0-h-u for 0<=u<=e_i;
for h in that interval all of these belong to the terminal zero
block modulo k.

For i=1, a_0!=0 and ell_1=0. Its columns in (2) kill the entire union

    [1,L_1] union [k-d_1+2,k-1].                            (3)

For 1<i<=s_0 its columns kill [ell_i+1,ell_i+L_i]. This is immediate
if a_(ell_i)!=0. In the second case, if z<e_i, B_i contains ell_i,
which again proves it by (2). If z>=e_i, the columns b+1,...,b+e_i
kill the contiguous interval [b+2,b+d_i-1]. It contains j because
z<=d_i-2. Adding the interval [j+1,j+L_i] killed by column j gives
[b+2,j+L_i], which contains [ell_i+1,ell_i+L_i] as b<ell_i<j.

The consecutive prefix intervals cover [1,sum_{i<=s_0}L_i], reaching
at least k-d_1+1. Together with (3) they cover every integer h from
1 to k-1. Thus for every such h at least one matrix in (1) has an
all-zero column and zero permanent. Only its h=0 summand remains;
it is nonzero. This proves the product lemma.

### 5. Triangular parametrization of an arbitrary affine source

Let A=gamma+U be an affine subspace of K^r. Successive coordinate
projections yield linear subspaces Y_i of K and linear maps
sigma_(ij):Y_j -> K, j<i, such that

    A = { (gamma_i+y_i+sum_{j<i}sigma_(ij)(y_j))_(i=1)^r
          : y_i in Y_i },
    dim A = sum_i dim Y_i.                                  (4)

For completeness, project U onto its first K-coordinate and choose a
linear section of that projection. Subtract its image to leave the
kernel, whose first coordinate is zero. Project that kernel onto the
second coordinate and choose another section, and continue. The
successive sections form a direct sum, and a section chosen at step j
has zero coordinates before j. Their coordinate functions give (4).

Let d_i=dim Y_i and Q_i:K -> Y_i be an image polynomial if d_i>0,
and the zero map if d_i=0. The linear maps sigma_(ij) composed with
Q_j are represented by linearized polynomials Q_(ij), by Section 2.
Define affine-linearized polynomials

    H_i(X_1,...,X_r)
      = gamma_i+Q_i(X_i)+sum_{j<i}Q_(ij)(X_j).

Their joint image on K^r is A.

Assume dim A>D. Let S={i:d_i>r+1}. Since every d_i<=k,

    sum_{i in S}d_i >= dim A-r(r+1) > kr/t+kt.

In particular |S|>=t. Choose I, of size t, consisting of t largest
dimensions within S. Since |S|<=r,

    sum_{i in I} (d_i-(r+1))
       >= (t/r) sum_{i in S} d_i-t(r+1)
        > k+kt^2/r-t(r+1)
       >= k,                                                (5)

where the last step uses k>=r(r+1)/t. For i in I we have
1<=e_i<d_i-1, and e_i+1<=r+1. The product lemma therefore applies
to the Q_i, i in I. It supplies a monomial

    M = product_{i in I} X_i^(beta_i),   wt(beta_i)=w_i,

with nonzero coefficient in Tr(product_{i in I}Q_i(X_i)^(a_i)).

We claim its coefficient remains exactly the same in the reduced
polynomial B(H_1,...,H_r). This will prove nonconstancy on A.
Expand a summand indexed by J, |J|=t, and a trace summand indexed
by h. Its product is

    product_{j in J} product_{u=0}^{e_j} H_j^(2^(u+h)),

reduced modulo X_i^q=X_i. Each of its w_j factors contributed by j
is an affine-linearized polynomial involving only variables X_i
with i<=j. Any selection of one monomial from each factor that
contributes to M must supply at least w_i factors to variable X_i
for each i in I, by binary-degree nonincrease. A positive exponent
does not reduce to zero, so a term involving a variable outside I
cannot contribute to M.

Charge a selected factor involving X_i by 1/w_i and a selected
constant by zero. The charge is at least

    sum_{i in I} w_i/w_i = t.

On the other hand w_i=r+2-i strictly decreases with i. A factor
from donor j can involve only i<=j and hence has charge
1/w_i<=1/w_j. The charge is at most

    sum_{j in J} w_j/w_j = t.

Equality forces every selected factor from donor j to use precisely
X_j: a constant has charge zero, and any i<j has strictly smaller
charge than 1/w_j. It also forces J=I. Thus the only contributions
to the coefficient of M come from J=I and from the diagonal
Q_i(X_i) terms in H_i; they are exactly the contributions in the
product lemma. Their sum is nonzero.

The reduced polynomial B(H_1,...,H_r) is therefore nonconstant by
Section 2. Its image is B(A), so B takes both binary values on A.
This completes the exact disperser theorem.

### 6. Parameters used by the submission

Let N>=1024 be the actual input length and write a=N^(1/5)>=4.
Let k be the least prime strictly larger than 2a^3, let m=N-k,
let r=ceil(m/k), and let t=ceil(sqrt(r)). Bertrand's postulate,
applied to floor(2a^3), gives

    2a^3 < k <= 4a^3.

In particular m>0, and

    a^2/4-1 <= r=ceil(N/k)-1 < N/k < a^2/2.

Thus r>=3 and 1<=t<r. Also

    r(r+1)/t <= sqrt(r)(r+1)
               < a(a^2/2+1) <= (3/2)a^3 < k.

The hypotheses of Section 1 hold. Since kr<N and k<=4a^3,

    kr/t <= k sqrt(r) < sqrt(kN) <= 2a^4,
    kt <= k(sqrt(r)+1) < 2a^4+4a^3,
    r(r+1) < a^4/4+a^2/2.

Consequently

    D=kr/t+kt+r(r+1)
       < (17/4)a^4+4a^3+a^2/2 < 6a^4 = 6N^(4/5).

The last inequality follows from a>=4: 4/a<=1 and
1/(2a^2)<=1/32, so the coefficient is at most
17/4+1+1/32=169/32<6. Define the exact integer

    d = floor(kr/t+kt+r(r+1))+1.

Then d>D and d<=6N^(4/5)+1. Padding m data bits with zeros into
r blocks of k bits gives an injective linear map into K^r. For
every irreducible monic degree-k modulus supplied in the input,
the resulting data-only function is a dimension-d affine disperser.
This assertion is independent of which irreducible modulus is used.
The validity test, total evaluator and input-length accounting for
the full input-modulus language are given in Part I.

## Part III. The finite circuit lower bound for an affine disperser

### Attribution and precise statement

The theorem proved in this part is the circuit lemma of Om Swostik Mishra and
Rahul Santhanam, *An Explicit 4n-O(log(n)) Size Lower Bound against General
Boolean Circuits*, ECCC TR26-260, revision 1, accepted 3 November 2026,
Theorem 4.1 and Sections 2--4:
https://eccc.weizmann.ac.il/report/2026/260/revision/1/download/.
The argument below reconstructs the parameter-pair method and its complete
case analysis. No new circuit-lower-bound theorem is claimed. In particular,
this submission does not use that paper's separate invocation of a
logarithmic-dimension explicit affine disperser. The proof in this part needs
no computational certificate and takes no unproved circuit lemma as a premise.

**Theorem C.** Let `h:F_2^M -> F_2`, where `M>=1`, and let `d` be an integer
with `1<=d<=M`. Suppose that `h` is nonconstant on every affine subspace of
`F_2^M` of dimension `d`. Every acyclic, single-output Boolean circuit for `h`
using all sixteen binary truth tables has at least

    4M - 10d

gates. The assertion also holds when constant sources and negated wires or
a negated output are allowed for free. Fan-out, depth and sharing are
unrestricted, and repeated input pins are allowed. Each gate is counted once.

Nonconstancy on all dimension-`d` subspaces also implies nonconstancy on every
larger-dimensional affine subspace: a constant larger subspace would contain
a constant dimension-`d` subspace. If `4M-10d<=0`, the claimed inequality
already follows from nonnegativity of size; the proof below covers the other
cases without any asymptotic assumptions on `M` or `d`.

### 1. The input circuit and the accounting representation

Fix an input circuit `C_0` of size `S` computing `h`. We may regard each wire
as carrying a value or its complement. A complement on an input pin is
absorbed into the receiving binary truth table. A complemented designated
output will be handled directly in the terminal argument, without assuming
the existence of an uncounted output gate. Constant, unary and degenerate
binary gates in the original circuit still contribute to `S`; deleting them
later is an accounted reduction. Thus the ordinary full `B_2` model embeds
in the more permissive representation with initial size exactly `S`.

All arithmetic below is in `F_2`. A binary Boolean function is either affine,
or has the form

    (z+alpha)(t+beta)+gamma,    alpha,beta,gamma in F_2.          (C1)

Indeed its algebraic normal form is `azt+bz+ct+e`; `a=0` is the affine case,
and for `a=1` completing the product gives (C1). An affine function may be
constant or unary. When the two sources of a gate coincide, or one source is
constant, the gate's value is a constant or a literal of the other source.
These observations classify all sixteen truth tables, including repeated
pins, without excluding any gate type from the input model.

A state of the proof consists of the following data.

* There are `m` independent parameters `u=(u_1,...,u_m)` and an injective
  affine map `Phi:F_2^m -> F_2^M` describing the original input assignment.
* `P` is a collection of `q` disjoint unordered pairs of parameters. Define

      Z(P)={u : pr=0 for every {p,r} in P},     a=m-q.           (C2)

  A parameter in a pair is called paired; other parameters are unpaired.
  Even a parameter no longer read by the state remains a coordinate of `Phi`
  unless an explicitly stated restriction removes it.
* An affine system `Av=Nu+c` determines `r` auxiliary values
  `v=(v_1,...,v_r)`. The matrix `A` is invertible. Every equation, after
  cancellations, contains at most three nonconstant symbols from `u,v`.
* A finite acyclic list of binary gates reads constants, parameters,
  auxiliary values and earlier gates. No gate feeds back into the affine
  system. The designated output is a constant or a literal of one of those
  values. Put `s=r+(number of these gates)`.

The essential invariant is

    C(u)=C_0(Phi(u))   for every u in Z(P).                    (C3)

Here `C(u)` means: solve the invertible affine system, then evaluate the
acyclic list. The state is an algebraic record of the proof, not a claim
that arbitrary affine maps are free gates in the original circuit. Its
equations originate in counted gates or accounted changes of variables;
`s` never increases. Initially `m=M`, `Phi` is the identity, `P` and the
affine system are empty, and the list is precisely `C_0`, so `s=S`.

A child of a symbol is a gate that reads it or an affine equation in which
it has coefficient one. Count distinct children, not pins. A gate's own
definition is not its child. In the affine system there are no designated
output symbols for individual equations: an equation is an unordered linear
constraint. The degree of a parameter is its number of children. The
designated output alone is not a child; whenever it is a literal of a
parameter or auxiliary value the terminal rule below applies before any
operation that could discard that value.

Call a parameter used if its degree is positive and isolated if its degree
is one. Let `I` denote the number of used unpaired parameters. Our potential
is

    W=s+4q+6I >= 0.                                         (C4)

In particular `W_0<=S+6M`. We will construct steps with

    Delta a in {0,-1},    Delta W <= 10 Delta a,              (C5)

each strictly decreasing the nonnegative integer `m+s`, unless a
dimension-at-least-`d` constant affine restriction of `C_0` has been found.
Changes are always new value minus old value.

### 2. Restrictions and changes of coordinates

The following operations preserve (C3) and the injectivity and affine form
of `Phi`.

**Unpaired fixing.** Set an unpaired parameter to a constant. This removes
one parameter, leaves `P` unchanged, and lowers `a` by one.

**Paired fixing.** For `{p,r} in P`, impose one of

    p=0,       r=0,       p=r+1.                              (C6)

Remove the substituted parameter and the pair. Each displayed line is
contained in `pr=0`. The other parameter becomes unpaired. Both `m` and `q`
decrease by one, so `a` is unchanged.

**Creating a pair.** For distinct unpaired parameters `x,y`, insert
`{x,y}` in `P`. This shrinks the allowed set to its subset where `xy=0`,
does not change `m`, and decreases `a` by one.

**Changing an unpaired coordinate.** An assignment

    x=w+ell(u excluding x),                                 (C7)

where `x` is unpaired and `ell` is affine, is a bijection between old and new
parameter spaces, with all coordinates other than `x` unchanged. The fresh
`w` is unpaired. Compose `Phi` with this bijection. It preserves `m,q,a` and
maps the new allowed set bijectively to the old one because paired
coordinates are unchanged. We use such a change only as implemented below,
so no long affine expression is expanded into a sparse equation.

Each fixing in (C6) or of an unpaired parameter also composes `Phi` with
an injective affine map from the remaining coordinates, whose image lies
in the old allowed set. Creating a pair uses the identity map. Direct
substitution in the state, or an operation proved equivalent to it below,
therefore preserves (C3). Fixing constants or replacing one parameter by
another plus a constant does not increase the number of symbols in any
equation, and leaves the auxiliary matrix `A` unchanged.

Complementing an unpaired coordinate is the special case `x=w+1` of (C7).
In equations this only changes constant terms; in gates the complement is
absorbed into their truth tables. It preserves all degrees and `W`.

### 3. Elementary state operations

#### 3.1 Moving an affine gate into the system

If an affine gate `J=b_1t_1+b_2t_2+c_0` has only constants, parameters or
auxiliary values as sources, replace it by the equation

    J+b_1t_1+b_2t_2=c_0

and make `J` an auxiliary value. After collecting terms this equation has
at most three symbols. Its enlarged auxiliary matrix has block form
`[[A,0],[*,1]]`, so is invertible. For every assignment the new auxiliary
value equals the old gate value, and all uses of `J` now refer to that
value. One gate becomes one equation, so `s,m,q,a` are unchanged. No unused
parameter becomes used: the equation mentions only the old gate's sources.
Thus `W` cannot increase. Moving gates in topological order transfers all
affine gates having no nonaffine ancestor into the system after finitely
many such moves.

#### 3.2 Exchanging a parameter and an auxiliary value

Suppose an auxiliary value `Q` depends on an unpaired parameter `x` in
the solved system. Its affine expression has the form

    Q=x+ell(u excluding x).

Use (C7), with the new parameter `w` equal to the old value of `Q`, and
exchange the roles of `x` and `Q`: `Q` becomes the fresh parameter and `x`
becomes an auxiliary value. Leave all equations and gates unchanged apart
from these role names. This gives exactly the same solutions in changed
coordinates and preserves (C3).

Here is a direct verification of invertibility. Let `j` be the old column
of `Q` in `A`, let `N_x` be the column of `x`, and put `t=A^{-1}N_x`.
The assumed coefficient of `x` in `Q` means `t_j=1`. The new auxiliary
matrix is `A B`, where `B` is the identity matrix with column `j` replaced
by `t`. Expansion in that column gives `det B=t_j=1`; the new matrix is
invertible. Row supports do not change, so sparsity is preserved.

No other parameter changes its children. The old `x` was used, since
`N_x!=0`, and the new parameter `Q` is used, since column `j` of the old
invertible `A` was nonzero. Both have weight six. Consequently

    Delta m=Delta q=Delta a=Delta s=Delta W=0.                (C8)

#### 3.3 Deleting a unary or constant gate

If a gate is identically a constant or a literal of one of its sources,
replace each use and the designated output, if necessary, by that constant
or literal and remove the gate. This applies to every gate with a constant
source, repeated sources up to complement, or a unary truth table. It also
applies to a nonaffine gate whose two source values are affine functions
with dependent linear parts: over `F_2` either one is constant, or the two
differ by a constant, so any function of them is unary in one source.
Finally, if a gate is constant on the current allowed set `Z(P)`, it may
be replaced by that constant, because (C3) only asserts equality there.
The later case analysis uses this last option only with an explicit reason
for constancy.

Each removal decreases `s` by one and leaves `m,q,a` fixed. A source to
which uses are redirected was already in that gate. Consequently a
previously unused parameter cannot become used, even when shared fan-out
causes the degree of an already used parameter to increase. Thus

    Delta W <= -1,    Delta a=0.                             (C9)

Complemented wires require no new gates: change each recipient's truth
table, or the output literal. If the output becomes an affine value or
a nonaffine gate with affine sources, the terminal rule in Section 4
applies before any further deletion.

#### 3.4 Deleting an irrelevant gate

A gate with no children that is not the designated output can be removed.
It has the same changes (C9). Repeated application discards irrelevant
parts of the acyclic list without making assumptions about how descendants
are shared.

#### 3.5 Deleting a short affine equation

An equation with at most two nonconstant symbols has an auxiliary symbol:
otherwise its row of `A` would be zero. It can therefore be written as

    v=U+c_0,

where `v` is auxiliary and `U` is zero or another single symbol. Replace
every occurrence of `v` by `U+c_0`, and delete `v` and this equation. This
does not expand a term into several terms, so all rows remain of support
at most three. On the auxiliary matrix, elimination clears the column of
`v` by row additions; it then has just the pivot one in the deleted row.
Deleting that row and column leaves an invertible matrix, as follows by
determinant expansion along the cleared column. The remaining assignments
have unique extensions and all remaining output values agree.

If the replacement is a parameter, it was already present in the removed
equation. Thus no unused parameter becomes used. The changes again obey
(C9). This argument includes a constant auxiliary value (`U=0`), and
includes cancellations when two occurrences become equal.

#### 3.6 Removing an equation read by an isolated unpaired parameter

Suppose an unpaired parameter `x` has exactly one child, an affine equation
`E`, and the terminal rule does not apply. Then its column `N_x` has exactly
one nonzero entry. Invertibility gives `A^{-1}N_x!=0`; choose any coordinate
`Q` with coefficient one and perform the exchange of Section 3.2.

The old symbol `x` is now auxiliary and occurs only in `E`: it has no
acyclic-gate use by the isolation assumption and is not the designated
output by the terminal check. Delete this auxiliary value and `E`.
Its column has just one nonzero entry, so determinant expansion shows
the remaining auxiliary matrix is invertible. The deleted equation uniquely
determines the discarded auxiliary value from the remaining data, without
affecting any retained value. The new parameter is `Q`; `m` and `q` have
not changed. Before the exchange `x` had weight six; afterwards `Q` has
weight at most six. No other parameter gains a child. The removal therefore
has (C9).

This deletion is why it is useful to count a sparse affine system as
equations rather than impose a directed output on each equation. It neither
discards a free parameter nor adds an uncounted circuit operation.

### 4. The terminal restriction

Assume `a>d`. Choose one parameter from each pair and set it to zero.
The resulting coordinate subspace `U_0` lies in `Z(P)` and has dimension
`m-q=a`. If the state's output is affine in the parameters, it is constant
on an affine subspace of `U_0` of dimension at least `a-1`: either it was
already constant, or take a nonempty affine level hyperplane.

The same conclusion holds when the output is a literal of a nonaffine
gate with affine source values. By (C1) that output has form

    (ell_1+alpha)(ell_2+beta)+gamma'

on parameters. If both factors are constant on `U_0`, the output is
constant there. If some factor is nonconstant, restrict it to zero, losing
exactly one dimension and making the output constant. Because `a>d` and
the dimensions are integers, `a-1>=d`.

The image of this subspace under the injective affine map `Phi` has the
same dimension. Equation (C3) makes `C_0` constant on that image. This
contradicts the assumption that `C_0` computes a dimension-`d` affine
disperser. Thus these terminal forms cannot arise while `a>d` for a correct
input circuit. They include a constant output, an input or auxiliary
literal as output, a final negation after all other gates disappear, and
a minimal nonaffine gate that is itself the output. No successor gate is
assumed to exist in those cases.

### 5. Three preliminary reductions for pairs

Before these reductions, use the terminal test and the operations of
Section 3 to delete any eligible gates or equations and move affine gates
with affine sources into the system. A deletion is already a complete
progress step of the form (C5), with `Delta a=0`. Pure moves of affine gates
do not themselves count as a progress step, but their number is bounded by
the finite current gate list. When no deletion is available, use the first
applicable reduction below.

#### P1. A paired parameter has degree different from one

Let `{p,r}` be its pair.

If `deg(p)=0`, fix `r=0` by paired fixing. Do not yet simplify gates.
The new unpaired parameter `p` is unused and has weight zero. No other
parameter becomes used merely by replacing `r` by a constant. Thus

    Delta m=-1, Delta q=-1, Delta a=0,
    Delta s=0, Delta W<=-4.

If `deg(p)>=2`, fix `p=0` and remove two distinct children of `p`. A child
gate has a constant source and is removable by Section 3.3. A child
equation, which previously had at most three symbols, now has at most two
and is removable by Section 3.5. Removing the first cannot increase the
number of symbols in the other; if some additional simplification deletes
the second first, its deletion has already been counted. The invertibility
arguments prevent a remaining equation from imposing a nontrivial relation
only on free parameters.

Only the surviving partner `r` can acquire a new individual weight, and
that weight is at most six. Therefore

    Delta m=-1, Delta q=-1, Delta a=0,
    Delta s<=-2, Delta W<=-2-4+6=0.

In both subcases `m+s` strictly decreases. Once P1 is inapplicable, every
paired parameter has exactly one child.

#### P2. Two members of a pair have the same child

Let `{p,r}` have common child `H`. By failure of P1, it is their only child.
If `H` is a gate, it reads only `p,r` (up to complements). Among the three
allowed assignments `(0,0),(0,1),(1,0)` two give the same Boolean value.
The affine line through that pair of assignments is one of (C6). Restrict
to that line and delete `H` as a constant gate.

If `H` is an equation, it has three symbols, since short equations were
already removed. Invertibility requires an auxiliary symbol, so it is
`v+p+r=c_0`. Impose `p=r+1`. It becomes `v=c_0+1`; remove the equation
and reroute `v` to this constant.

In either case the remaining parameter of the pair is unused: its only
child was `H` and the replacement was constant. No originally unpaired
parameter becomes used. Hence

    Delta m=-1, Delta q=-1, Delta a=0,
    Delta s=-1, Delta W<=-1-4=-5.

After P1 and P2 have failed, the two members of each pair are isolated
and have different children.

#### P3. A paired parameter is mentioned by an affine equation

Let `p` have the unique equation child `E_p`, and let its partner `r`
have unique child `H_r`. By failure of P2 these are different. Fix `r=0`
and remove `H_r`, using Section 3.3 or 3.5 as appropriate.

The deletion cannot introduce a use of `p`, because `p` did not occur
in `H_r`: the replacement is a constant or one of the other symbols of
`H_r`. Thus `p` is now unpaired and still has at most the one child
`E_p`. If its equation has become short, remove it by Section 3.5.
Otherwise remove it by Section 3.6. In the latter operation `p` is
replaced as a parameter by an auxiliary value. In either event the
additional weight created by dissolving the pair is at most six, and no
previously unused unpaired parameter becomes used.

There are two distinct removed children. Therefore

    Delta m=-1, Delta q=-1, Delta a=0,
    Delta s<=-2, Delta W<=-2-4+6=0.

It is possible that a preceding simplification already removes `E_p`;
then both removals are already accounted for and the estimate is only
stronger. If `p` itself becomes unused, its weight is zero. These
possibilities do not require another parameter fixing.

When none of P1--P3 applies, every paired parameter is isolated, occurs
only in the acyclic gate list, and has a different child from its partner.
In particular, every parameter mentioned in the affine system is unpaired.

### 6. Exposing a minimal nonaffine gate

Suppose `a>d`, no terminal form has occurred, all deletions are exhausted,
and none of P1--P3 applies. The gate list is nonempty; otherwise the output
would be an affine value. Choose its first gate `G` in topological order.
It cannot be affine, because all affine gates with affine sources have
been transferred to the system. It is therefore of form (C1), and its
sources are parameters or auxiliary values. Their linear parts are
independent: if dependent, the gate would be unary or constant in one
source and Section 3.3 would apply.

We can make these two sources into distinct parameters using at most two
exchanges from Section 3.2. If a source `Q` is auxiliary, its nonzero linear
part depends on a parameter `x`. Such `x` is mentioned in the affine
system (a zero column `N_x` could give no dependence), so it is unpaired
by failure of P3. Exchange `Q` with `x`.

After the exchange, the newly exposed source is a used unpaired parameter.
If the other source is still auxiliary, its linear part is independent of
the first source's linear part, since an invertible affine change preserves
independence. Consequently it depends on a parameter other than the first
exposed source. That parameter is again mentioned in the system and is
unpaired. Exchange this second source with such a parameter. This argument
also covers the possibility that the first exchange turned an originally
parameter-valued second source into an auxiliary value. After at most two
exchanges both inputs are parameters, and neither exchange changes `W`.

If an exchange makes a deletion, a terminal form or one of P1--P3
available, carry it out and finish the progress step there. This avoids
any loop of coordinate exchanges. Otherwise we reach a gate

    G=(x+alpha)(y+beta)+gamma                              (C10)

with distinct parameter sources. The terminal rule rules out `G` as the
designated output. As irrelevant gates have been removed, `G` has at
least one child. Every such child is in the acyclic gate list, so is
different from `G`; no gate feeds the affine system.

### 7. Removing the exposed gate

#### A. Both sources are unpaired

Complement unpaired `x` and `y`, if necessary, so that (C10) becomes
`G=xy+gamma`. These are the harmless coordinate changes specified in
Section 2. Insert `{x,y}` as a new pair. On the new allowed set,
`G=gamma`. Remove `G` and at least one of its children; the latter has
a constant source after `G` is replaced. These are two distinct gates.

Both `x,y` had positive degree, so each had weight six. Pairing removes
both individual weights and adds four through `q`. No other unused
parameter becomes used when a constant gate and a unary child are
removed. Thus

    Delta m=0, Delta q=1, Delta a=-1,
    Delta s<=-2, Delta W<=-2+4-12=-10.                      (C11)

#### B. Exactly one source is paired

Write `x` for the unpaired source, `y` for the paired one, and `r` for
the partner of `y`. Failure of P1--P3 says that `y` has only child `G`,
while `r` has one gate child `H_r` different from `G`.

First fix `r=0` by paired fixing. Next fix the now-unpaired `x` to
`alpha`, making `G=gamma`. Remove `G`, and remove `H_r` because it has
the constant source `r`. The gate `H_r` may be a descendant of `G`,
and its other source may coincide with another named value; this causes
no double count, since `H_r!=G`, and it remains removable in all cases.
The remaining `y` is unused after its sole child `G` disappears.

The initially used unpaired `x` loses weight six; `y` has weight zero
both before the whole step (paired) and after it (unused). No other
unpaired parameter becomes used. Hence

    Delta m=-2, Delta q=-1, Delta a=-1,
    Delta s<=-2, Delta W<=-2-4-6=-12<=-10.                  (C12)

#### C. Both sources are paired

The sources cannot be partners of each other, since they have the same
child `G` and P2 was inapplicable. Write their distinct pairs as
`{x,r_x}` and `{y,r_y}`. All four parameters are distinct. The unique
child `H_x` of `r_x` is a gate different from `G`.

Fix both `r_x` and `r_y` to zero by paired fixing, and then fix the newly
unpaired `x` to `alpha`. Remove the resulting constant gate `G` and
remove `H_x`, which has a constant source. It is unnecessary to claim
any additional removal for the child of `r_y`: that child may coincide
with `H_x`. Only the distinct gates `G,H_x` are counted. The remaining
`y` becomes unused because `G` was its only child.

No newly unpaired parameter has positive weight at the end: `x` was
fixed and `y` is unused. No initially unused unpaired parameter becomes
used. It follows that

    Delta m=-3, Delta q=-2, Delta a=-1,
    Delta s<=-2, Delta W<=-2-8=-10.                        (C13)

These cases exhaust all possibilities for the two sources. Their
estimates compare the endpoints of the entire grouped step; temporary
individual weights after dissolving a pair are not ignored or counted
twice.

### 8. Termination and the lower bound

Every progress step above has `Delta a` zero or minus one, obeys (C5),
and strictly decreases `m+s`. No operation increases `s` or `m`.
Deleting a gate or equation lowers `s`; creating a pair in case A lowers
`s` by at least two; P1--P3 and cases B,C lower `m`. Moving an affine
gate into the system changes neither size, and can be performed only
finitely often before the current gate list is exhausted. Between progress
steps at most two exchanges are needed to expose a chosen gate, or an
earlier applicable progress step is taken. Thus the construction cannot
continue forever without either reducing `a` to at most `d` or producing
the terminal affine subspace.

For a direct contradiction assume `S<4M-10d`. Then

    W_0 <= S+6M < 10(M-d).

After any number of nonterminal progress steps, summing (C5) yields

    0 <= W <= W_0+10(a-M) < 10(a-d).                       (C14)

Therefore `a>d` always. The nonnegative integer `m+s` cannot decrease
indefinitely, so a terminal case must occur. Section 4 supplies an affine
subspace of dimension at least `d` on which `C_0` is constant, contradicting
that it computes `h`. This proves `S>=4M-10d` and hence Theorem C.

### 9. Scope of the proof

The proof is a universal argument on a finite circuit, not an enumeration
of selected configurations. Arbitrary shared descendants are retained in
the state. A deletion counts an existing gate or equation once, and the
case analysis explicitly identifies the distinct objects it charges.
Degree counts use distinct children, so repeated pins do not generate
spurious credit. Normalization includes every affine, constant and unary
truth table, affine dependencies between sources, constant sources after
substitution, and final-output cases. The enlarged bookkeeping representation
is used only with the invariant (C3), the exact initial size `S`, and
nonincreasing size under the operations proved above. It supplies no
unproved direct-sum assertion and never identifies a potential value with
the original gate count.

The theorem applies separately at any given integer input length and has
no uniformity requirement on the circuit. To obtain an eligible family
one must additionally supply a uniform language construction and prove
its affine-disperser property. Those obligations are proved in Parts I and II of this document.
