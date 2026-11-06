# Deterministic prime-degree fields over F2

This is a prime-degree specialization of the deterministic construction
route of Victor Shoup, *New Algorithms for Finding Irreducible Polynomials
over Finite Fields*, Mathematics of Computation54(189),1990,435--447,
https://www.shoup.net/papers/detirred.pdf, Sections2--3. The supplied
algorithm and argument below prove the particular statement used here;
no appeal to an external factoring or irreducible-construction oracle is
needed. All operations used by the constructor are specified below.

## 1. Extract an irreducible factor using binary linear algebra

Polynomials are represented by binary integers, with bit j the coefficient
of X^j. Euclidean division, multiplication and gcd are exact finite bit
operations. On a squarefree monic polynomial f of degree D, form the
D-dimensional quotient algebra A=F2[X]/f. The map b->b^2 is F2-linear.
Its matrix has column j the coefficients of X^(2j) modf. Gaussian
elimination computes the kernel V of that matrix minus the identity,
with pivots and free columns in increasing order and the usual free-bit
basis. The factor algorithm returns f if dimV=1. Otherwise it selects
the first basis vector b other than0 or1, computes g=gcd(f,b), divides
h=f/g, and repeats on the numerically smaller of g,h. All choices are
fixed and every repeat strictly reduces degree.

For correctness, write the mathematical factorization f=product f_i
into distinct monic irreducibles. Existence follows by repeated division
in the Euclidean domain; this does not assume a factoring algorithm.
The Chinese remainder map A->product F2[X]/f_i is an isomorphism:
Euclid gives inverses of the pairwise coprime factors and hence explicit
idempotents for each component. Each component is a field, whose only
solutions of z^2=z are0 and1. Thus V corresponds exactly to one binary
choice per irreducible factor, and its dimension equals their number.
Dimension1 means f is irreducible. Otherwise a nonconstant b in V has
some component0 and some1. Therefore its gcd g includes exactly the
zero components, while h includes the others. Both have positive degree
less than D, are squarefree, and divide f. Repetition terminates at an
irreducible divisor. A nonconstant basis vector exists whenever dimV>1,
since the constants form a one-dimensional subspace. The construction
requires no enumeration of quotient-field elements.

Each matrix has D rows of D bits. Its entries and Gaussian elimination
are computable in polynomial time. There are at most D descents. A
coarse O((D+1)^5) bit-operation bound covers this elementary version;
no fast factorization algorithm is required.

## 2. A cyclotomic field and a nonresidue

Let k>=3 be odd prime. Compute m, the least positive integer with
2^m=1 modk, by repeated multiplication modulo k. Fermat gives
m|k-1, so2<=m<=k-1 and gcd(m,k)=1. Set e=v_k(2^m-1), using repeated
exact integer division. In particular1<=e<=m. Start with the polynomial
Phi_k=1+X+...+X^(k-1), and extract an irreducible factor f_1 by Section1.
For i=2,...,e form f_(i-1)(X^k) and extract an irreducible factor f_i.
Every f_i has degree m, as shown next. At the end put f=f_e,
K=F2[A]/f and a=A modf.

We use the elementary cyclicity of the multiplicative group of a finite
field. For completeness, for every prime-power ell^u dividing its order
Q-1, some element has order divisible by ell^u. Otherwise every nonzero
element would be a root of X^((Q-1)/ell)-1, contradicting the root bound.
Taking the appropriate powers isolates such an ell-primary element.
The product of these elements for the different primes has order Q-1:
commuting elements of coprime orders have product order equal to their
product, as follows by raising to each order. Lagrange's theorem supplies
the divisibilities used here. Thus a finite field contains all d roots
of unity when d divides Q-1, and primitive ones exist.

The roots of Phi_k are precisely the primitive k-th roots of unity. Their
Frobenius orbits have size ord_k(2)=m, so every irreducible factor of Phi_k
has degree m. This polynomial is squarefree since X^k-1 has derivative
X^(k-1) in characteristic2 and does not vanish at its nonzero roots.
Inductively a root beta of f_(i-1)(X^k) has beta^k primitive of order
k^(i-1), so beta has order exactly k^i. For i<=e, k^i divides2^m-1, so
all these roots are in the field of2^m elements. Their orbit length is
ord_(k^i)(2), which divides m and is divisible by ord_k(2)=m, hence is m.
The composition is squarefree: its roots are nonzero, k is odd, and the
original irreducible polynomial has simple roots. Thus every factor
extracted by Section1 has the required degree and primitive-root order.
In particular a has order k^e in K. It is not a k-th power in K: a k-th
root would have order k^(e+1), impossible since that does not divide
|K*|=2^m-1. These statements apply also to e=1, with no lifting loop.

All factored polynomials have degree at most mk<k^2. There are e<=k
factor extractions. Neither an integer factorization of2^m-1 nor an
exponentially large polynomial is used; only division by the given k
and polynomials of degree less than k^2 occur.

## 3. A prime extension and descent by trace

The polynomial Y^k-a is irreducible over K. To prove this directly, take
a root beta in an algebraic closure. K contains the k-th roots of unity.
Put q=2^m. Since beta^k=a and a^q=a, beta^q/beta is a k-th root of unity.
It cannot be1, as that would put beta in K (the roots of Z^q-Z are exactly
K), contradicting the nonresidue property. As k is prime it has order k.
Frobenius^m fixes that ratio. The conjugates beta,beta^q,... therefore
form an orbit of length exactly k. Their minimal polynomial over K has
degree k and divides Y^k-a, proving irreducibility.

We may consequently calculate in E=K[Y]/(Y^k-a), representing elements
as k coefficients in K. Its dimension over F2 is mk, which is less than
k^2. Let beta be the residue of Y. Define the explicit element

    gamma=sum_(i=0)^(m-1) beta^(2^(ki)).

The automorphism sigma:z->z^(2^k) permutes this sum cyclically, because
sigma^m is the identity on E. Thus gamma^(2^k)=gamma. Also gamma is not
in K. Indeed write2^(ki)=k u_i+r_i with1<=r_i<k. The residues r_i are
distinct as i runs0,...,m-1: Fermat gives2^k=2 modk, so they are exactly
the first m powers of2 modulo k. Therefore in the K-basis1,beta,...,
beta^(k-1),

    gamma=sum_i a^(u_i) beta^(r_i)

has distinct nonzero positive-degree terms, with no constant term.
It cannot lie in K. In particular gamma is not in F2. Since its binary
Frobenius orbit has size dividing prime k, that size is exactly k.
Hence its minimal polynomial over F2 is

    Q_k(Z)=product_(j=0)^(k-1)(Z+gamma^(2^j)).

The product is computable using represented E arithmetic. Frobenius
permutes the factors, so each coefficient is fixed by squaring and is
therefore0 or1. It is monic, degree k and irreducible, because gamma has
orbit length k. This defines one completely deterministic Q_k for every
odd prime k, with no field-polynomial input, advice, random choice,
number-theory conjecture or external algorithm.

The represented arithmetic never enumerates E. A K element has m bits,
and an E element has k such coefficients. Reduce Y^k to a in products.
For squaring, characteristic2 kills all cross terms and a term of degree
2i>=k needs just one such reduction. Compute successive sigma powers by
k squarings, repeated m times, rather than writing enormous exponents.
Compute the final product with ordinary polynomial multiplication over E.
These exact procedures use polynomially many bit operations. Combining
Sections1--3, the generous bound O((k+1)^12) covers the implementation;
its largest factor matrix has fewer than k^2 rows and columns. Only the
selected degree's Q_k need be constructed by the language evaluator.


## 4. Retained compressed-prefix reference class

Let k>=5 be odd prime, K a field of2^k elements by the elementary finite-field facts in proof.md Section2, and
Tr(z)=sum_(i=0)^(k-1)z^(2^i). Tr takes values in F_2 and Tr(1)=1,
so its kernel has M=2^(k-1) elements. There exists nonzero alpha with
Tr(alpha)=Tr(alpha^-1)=0. Suppose otherwise. The monic polynomial

    T(X)=Tr(X)/X=sum_(i=0)^(k-1)X^(2^i-1)

has precisely the M-1 distinct nonzero trace-zero elements as roots.
All would then have inverse trace1. Thus T divides the polynomial

    U(X)=X^M(Tr(1/X)+1)=X^M+sum_(i=0)^(k-1)X^(M-2^i).

U is monic degree M, while T has degree M-1, so the quotient is X+c.
The coefficient of X^(M-1) forces c=1. But U has coefficient1 at
X^(M-2), while (X+1)T has coefficient0 there. Its exponents are powers
of2 or powers of2 minus1; M-2 is neither when M>=8. This contradiction
proves the claim using polynomial division and the elementary root bound.

Neither alpha nor alpha^-1 is1 since Tr(1)=1. Replace alpha by alpha^-1
if necessary to arrange also Tr(1/(alpha+1))=0. This is always possible:

    Tr(1/(alpha^-1+1))=Tr(alpha/(alpha+1))
                       =1+Tr(1/(alpha+1)).

The first two zero traces are preserved under inversion. Since k is
prime and alpha is not in F_2, its Frobenius orbit has size k. Its monic
irreducible polynomial Q has constant1, coefficient of Z^(k-1) equal
Tr(alpha)=0, and coefficient of Z equal Tr(alpha^-1)=0. These follow by
expanding the product over the conjugates; the product of conjugates is
nonzero in F_2, hence1. Also Q(1)=1 since Q has no root1. The logarithmic
derivative identity, obtained by differentiating that product, gives

    Q'(1)/Q(1)=sum_(i=0)^(k-1)1/(1+alpha^(2^i))
              =Tr(1/(1+alpha))=0.

Consequently Q+1 is divisible by Z^2 (its constant and linear coefficients
vanish) and by (Z+1)^2 (its value and derivative at1 vanish). These
coprime factors give

    Q=1+Z^2(Z+1)^2 R=1+(Z^4+Z^2)R.

R is monic degree k-4. Its coefficient at degree k-5 is0 because Q's
coefficient at degree k-1 is0. All remaining coefficients are exactly
k-5 free bits, as in the decoder. Hence every eligible k has at least
one encoded irreducible polynomial. Multiplication by Z^4+Z^2 is
injective in F_2[Z], so different bit prefixes encode different Q. The
prefix test from proof.md Section2 remains exact for every decoded monic Q.
This proof constructs no field polynomial by advice or an external
algorithm: it proves existence for restricting circuits; evaluation just
decodes its actual input and checks it in polynomial time.

