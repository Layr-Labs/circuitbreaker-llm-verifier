# Uniform field construction and coarse prime intervals

In an algebraic closure of F2, the polynomial Z^(2^u)-Z has exactly2^u
distinct roots, since its derivative is1. Those roots are closed under
addition, multiplication and nonzero inversion, and form a field.
Frobenius z->z^2 permutes them. A Frobenius orbit gives a monic irreducible
polynomial over F2 by taking the product of its linear factors: the
coefficients are fixed by Frobenius, and a proper binary factor would
require a proper invariant subset of the orbit. Conversely the roots of
an irreducible polynomial form one orbit. These elementary facts apply
to every degree u and underlie the supplied constructor.


## Appendix A. Deterministic prime-degree fields over F2

This is a prime-degree specialization of the deterministic construction
route of Victor Shoup, *New Algorithms for Finding Irreducible Polynomials
over Finite Fields*, Mathematics of Computation54(189),1990,435--447,
https://www.shoup.net/papers/detirred.pdf, Sections2--3. The supplied
algorithm and argument below prove the particular statement used here;
no appeal to an external factoring or irreducible-construction oracle is
needed. All operations used by the constructor are specified below.

### A.1. Extract an irreducible factor using binary linear algebra

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

### A.2. A cyclotomic field and a nonresidue

Let k>=3 be odd prime. Compute m, the least positive integer with
2^m=1 modk, by repeated multiplication modulo k. Fermat gives
m|k-1, so2<=m<=k-1 and gcd(m,k)=1. Set e=v_k(2^m-1), using repeated
exact integer division. In particular1<=e<=m. Start with the polynomial
Phi_k=1+X+...+X^(k-1), and extract an irreducible factor f_1 by Appendix A.1.
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
extracted by Appendix A.1 has the required degree and primitive-root order.
In particular a has order k^e in K. It is not a k-th power in K: a k-th
root would have order k^(e+1), impossible since that does not divide
|K*|=2^m-1. These statements apply also to e=1, with no lifting loop.

All factored polynomials have degree at most mk<k^2. There are e<=k
factor extractions. Neither an integer factorization of2^m-1 nor an
exponentially large polynomial is used; only division by the given k
and polynomials of degree less than k^2 occur.

### A.3. A prime extension and descent by trace

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
Appendix A.1--A.3, the generous bound O((k+1)^12) covers the implementation;
its largest factor matrix has fewer than k^2 rows and columns. Only the
selected degree's Q_k need be constructed by the language evaluator.



## Coarse prime intervals from a familiar standard theorem

Use the classical prime number theorem in its standard form
pi(x)~x/log x, for real x tending to infinity (NIST DLMF27.2(iii),
https://dlmf.nist.gov/27.2 , Eq27.2.3). This familiar standard result is
stated explicitly rather than re-proved. Its ONLY needed consequence
is proved here: for every fixed real r>1,

 pi(rx)-pi(x)=(r-1+o(1))*x/log x>0

for EVERY sufficiently large real x. Hence a prime lies in(x,rx], and
the least prime>=x has ratio to x tending to1. For each fixed relative
error this holds at EVERY threshold above one constant, so it also
applies uniformly to public context-dependent thresholds whose minimum
tends to infinity. No effective threshold, prime oracle or acceptance
verdict is assumed. All searches below use actual trial division, and
finite exceptions enlarge a fixed runtime constant or take stated guarded
fallbacks. This supplies every relative-prime estimate used by the
constructors and the current fixed-tuple proof.

# Field preliminaries and the trace-cubic extractor

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



# Prime-degree trace-cubic disperser: explicit isolation proof

Construction and isolation method originate in Ben-Sasson--Kopparty, Affine Dispersers from Subspace Polynomials, Section4.1. The interval argument below is stated with integer endpoints and proves the needed coefficient isolation directly, including the equal-island case where a permanent vanishes without an all-zero column.

Let K=F_(2^k), k an odd prime. For every nonzero mu in K, f_mu(x)=Tr(mu*x^7) is nonconstant on every affine binary subspace of dimension at least d0=ceil(2k/5)+10, whenever d0<=k. If d0>k the assertion is vacuous. The function has Boolean degree at most3, since x^7=x*x^2*x^4 and Frobenius and trace are binary linear.

It suffices to consider dimension exactly d=d0. Parameterize an affine space as Q(K)+gamma, where Q=sum_i a_i X^(2^i) is the monic image polynomial of its direction space, of degree2^(k-d). The supplied subspace-polynomial construction proves a0 and a_(k-d) nonzero, all a_(k-d+1),...,a_(k-1) zero, and that this final zero interval of length d-1 is the only cyclic zero interval of that length. The last property follows from primality: any second such interval would make a nonidentity Frobenius power of Q proportional to Q, forcing all ratios in its image into F2, contrary to d>1. No cyclic zero interval has length d, by the root bound after rotating coefficients. These field preliminaries are proved in the image-polynomial section above.

For three distinct columns B in Z/kZ, the coefficient of X^(sum_{b in B}2^b) in the reduced polynomial Tr(mu*(Q(X)+gamma)^7) equals

 sum_(rho=0)^(k-1) mu^(2^rho)
    * perm((a_(b-u-rho)^(2^u))_(u=0,1,2;b in B))^(2^rho).

Indeed expand the three Frobenius factors, then the trace. A weight-three reduced exponent uses three distinct residues: a repeated residue causes a carry and lowers binary weight. Constants also lower weight. Thus the contributing assignments are precisely the permanent, and reduction or gamma cannot introduce further contributions of weight three. Permanent equals determinant in characteristic two. We will choose B whose rho=0 permanent is nonzero and all other permanents are zero. Then the displayed coefficient is nonzero for every mu!=0. Uniqueness of the reduced field polynomial establishes nonconstancy.

All coefficient indices below are cyclic modulo k. Put L=k-2d+4 and U=d-4.

First suppose some nonzero a_j has L<=j<=U (with0<=j<=k-d as necessarily holds). Take B={0,j+1,k-d+2}. These are distinct columns. At rho0, the first column is supported only in row0, the last only in row2, and the middle has nonzero row1. Its permanent is nonzero. For1<=rho<=d-3 the first column lies wholly in the final zero interval. For k-d+3<=rho<=k-1 the last column does. For every remaining rho, d-2<=rho<=k-d+2, each of j+1-rho,j-rho,j-1-rho belongs to[-d+1,-1]: the bounds j<=d-4 and j>=k-2d+4 give both endpoints. Thus the middle column is zero. Empty rho intervals need no argument. This settles this case.

Otherwise there are no nonzero coefficients in[L,U]. Since a0 and a_(k-d) are nonzero, there are two outer support islands. Let a be the last nonzero index before L, and b the first nonzero index after U. They exist; otherwise an endpoint would belong to the prohibited interval. Every nonzero coefficient is in [0,a] or[b,k-d], and both endpoints of each island are nonzero. Write their lengths as A=a+1 and C=k-d-b+1, and the intervening zero gaps as G=b-a-1 and H=d-1. Then

 A,C<=k-2d+4,
 G>=3d-k-7,
 G-max(A,C)>=5d-2k-11>=39,
 H>G.

The last strict inequality uses uniqueness of the longest zero gap. In particular both zero gaps are longer than either island by more than4. Internal zeros in either island only help the following support arguments.

If A>C take B={0,k-d+1,a+2}. At rho0, rows0,1,2 have nonzero pivots a0,a_(k-d),a_a, respectively, after this column ordering; every above-pivot entry is zero, so the permanent is nonzero. Consider its first and third columns at arbitrary rho. A nonzero coefficient u in the first and v in the third would have cyclic difference v-u in[A-1,A+3]. Both zero gaps exceed max(A,C)+4, so those two support points must lie in one island, without wraparound. The shorter island cannot contain such a difference. In the longer island the sole possible difference is A-1, attained only by u=0,v=a. The respective rows must be0 and2, forcing rho=0. Thus for every nonzero rho, one of these columns is wholly zero.

If C>A use B={b,a+1,k-d+2}. The same argument now uses its first and third columns with differences in[C-1,C+3]. Only the endpoints b,k-d of the longer second island can occur, again in rows0 and2 at rho0. At rho0 the middle row1 pivot is a_a; its other required upper entries are zero. This gives the same isolation.

Finally let A=C and take B={0,b+1,a+2}. Its rho0 permanent has pivots a0,a_b,a_a in rows0,1,2 and is nonzero. The first/third-column support argument leaves only two possible shifts: rho0, matching endpoints of the first island, and rho=-b modk, matching endpoints of the second. At the second shift the first column is supported only in row0 and the third only in row2. A nonzero permanent would therefore require the middle row1 coefficient a_(2b). Equal island lengths give

 k=2A+G+H, b=A+G, k-d=b+A-1.

Since H>G, we have k-d<2b<k. Hence a_(2b)=0, and that remaining permanent vanishes as well. This case needs only a zero permanent, not an all-zero column.

In all cases exactly one trace summand survives at the chosen weight-three exponent. Its coefficient is nonzero. This proves the asserted dispersion at d0. Together with the corresponding section of cubic.md, it yields a quantitative affine extractor whenever the source dimension exceeds d0 by a linear fraction of k. The following vector-output argument and cubic.md supply the output and uniform field-degree details.


# Linear bias-to-flat bound for cubics

We give a self-contained binary specialization of Amichai Lampert,
“Slice rank and analytic rank for trilinear forms,” Discrete Analysis
2025:25, arXiv:2404.19704v2 (8 October2025). The weaker constants below
are derived here; no external rank theorem is a premise.

A slice of a trilinear T:U×V×W->F2 is a product of a linear form in
one argument and a bilinear form in the other two. Its slice rank is
the least number of slices summing to T. Restricting an argument to a
codimension-c subspace increases the rank by at most c when lifting:
choose coordinates extending that subspace and collect the terms with
each removed coordinate. This also proves every lifting assertion below.

**Dense bounded-rank lemma.** If A(u)=T(u,.,.) satisfies
Pr_u[rank A(u)>t]<1/(4t), for integer t>=1, then T has slice rank<=4t.

Write E={u:rank A(u)<=t}, with complement density e<1/(4t).
If no element of E has rank t, induction at t-1 applies. At t=1, if
all A(u) on E vanish, the linear map A vanishes everywhere: a proper
kernel has density<=1/2, whereas E has density>1/2.
Otherwise choose u0 in E of rank t and write
A(u0)=sum_(i=1)^t alpha_i(v) beta_i(w) with independent forms.
Let V0=intersection ker alpha_i and W0=intersection ker beta_i.
For u in E intersection(E+u0), both B=A(u) and B+A(u0) have rank<=t.
If B|V0×W0 has rank s, choose an invertible s-block there and eliminate
its other row/column entries. These operations add kernel vectors for
A(u0), so they leave A(u0) unchanged. Thus B splits into that s-block
and B0, while A(u0) remains rank t in the complementary block. Hence

 rank B0<=t-s,
 t>=rank(B+A(u0))>=s+t-rank B0>=2s.

The restricted tensor therefore has rank<=floor(t/2) at all
u in E intersection(E+u0), except a set of density<2e. For t>=2,
2e<1/(2t)<=1/[4floor(t/2)], so induction bounds its slice rank by
4floor(t/2)<=2t. Lifting V0,W0 adds at most2t slices, proving4t.
For t=1 the restricted linear map is zero on a set of density>1/2,
hence everywhere; lifting costs only2 slices. This completes induction.

**Moment lemma.** Put Z={(u,v):T(u,v,.)=0}, density mu>0, r=-log2 mu.
Then T has slice rank<=9r+49.

Let R_u=rank T(u,.,.) and R'_v=rank T(.,v,.). Under the UNIFORM law
on Z, the u marginal is 2^(-R_u)*|V|/|Z|. Consequently

 E_Z 2^R_u=|U||V|/|Z|=2^r=E_Z 2^R'_v.

Conditioning this sum on v gives some v0 with

 2^R'_(v0)+E_(u in U')2^R_u<=2^(r+1),
 U'=ker T(.,v0,.).

Thus codim U'<=r+1 and Pr_(u in U')[R_u>t]<2^(r+1-t).
Choose t=2ceil(r)+10<=2r+12. The latter probability is less than
2^(-r-9)<1/(8r+48)<=1/(4t); the elementary inequality
2^(r+9)>8r+48 holds for every r>=0. The dense lemma gives rank<=4t
on U'. Lifting U' adds<=r+1, giving9r+49 as asserted.

For a Boolean quadratic, let R be its alternating polar rank.
Symplectic Gaussian elimination expresses it as sum_(i<=R/2)x_i y_i
plus an affine function. Nonzero absolute bias delta forces
R<=2log2(1/delta), by summing signs over the independent pairs.
It also gives a constant affine subspace of codimension<=R+1:
fix the at most R+1 linear forms on which this expression depends.
These statements include affine/constant quadratics.

**Cubic theorem.** A Boolean cubic f with absolute bias delta>0 is
constant on an affine space of codimension<=128ceil(log2(2/delta)).

Write b(g)=E_x(-1)^g(x) and D_a g(x)=g(x+a)+g(x). Averaging a
squared signed bias gives b(g)^2=E_a b(D_a g). Two applications of
Jensen therefore yield

 delta^4<=E_(a,b)b(D_b D_a f),
 delta^8<=E_(a,b)b(D_b D_a f)^2
          =E_(a,b,c)(-1)^T(a,b,c),
 T(a,b,c)=D_c D_b D_a f(0).

For a multilinear Boolean representative of degree<=3, T is trilinear
and independent of the base point. Averaging in c shows its last
expectation is exactly Pr_(a,b)[T(a,b,.)=0]. Its analytic parameter r
is therefore<=8log2(1/delta). The moment lemma decomposes T into at
most9r+49 slices. Identify its three argument spaces with the original
space and intersect the kernels of ALL their unary forms, obtaining
L of codimension<=9r+49 and T|L³=0. Thus f is quadratic on EVERY coset
of L: a remaining cubic monomial in a basis of L would have a nonzero
third derivative in three corresponding basis directions. Translation
does not change that cubic part. This uses no division by3! in F2.

Some coset has absolute bias>=delta, since the global signed bias is
the average of its coset biases. The quadratic fact gives a constant
subspace there with total codimension at most

 9r+49+2log2(1/delta)+1<=74log2(1/delta)+50
                          <=128ceil(log2(2/delta)).

If this exceeds ambient dimension, a point suffices. This proves all
cases without inferring dimensions from finite tests.

If a cubic is nonconstant on every dimension-d0 affine subspace, its
restriction to any dimension-D affine source, D>d0, consequently has

 absolute bias <4*2^(-(D-d0)/128).                    (LIN-BIAS)

Otherwise ceil(log2(2/delta))<=(D-d0)/128 would supply a constant d0-flat.
The zero-bias case is automatic.


# A vector trace-cubic extractor at entropy one half

This section combines the supplied prime-cubic isolation and cubic-bias arguments with the deterministic field constructor and prime-interval proof in cubic.md. All thresholds below are eventual; no explicit universal small-w threshold is claimed.

At input length2w, choose the least prime k>=2w. The supplied analytic prime theorem implies k<11w/5 for every sufficiently large integer w. Trial division finds this k in polynomial time: the eventual linear bound, together with the finitely many smaller w, gives one global constant upper bound on k/w. Construct the deterministic irreducible polynomial of degree k by the prime-degree constructor. Interpret the2w input bits as the first2w coordinates of K=F2[Z]/Q_k, filling all other coordinates with zero. Compute x^7 in K. Output its first m=floor(w/10000) binary coordinates (or fewer if necessary for the finitely many small lengths). Call this function B_w. It is one total deterministic polynomial-time family; choosing any fixed total convention for small w does not affect the eventual statement.

The code returns a tuple only when k is prime, k>=2w,
m=floor(w/10000)>=1, 5k<=11w and10(w-d0)>=w. For EVERY returned
tuple the following argument is finite; prime asymptotics only prove
that these guards eventually pass. Indeed the eventual prime bound
implies w-d0>=3w/25-11>=w/10.

Let the input be uniform on any affine subspace of binary dimension D>=w. Its embedding in K has the same dimension. Put d0=ceil(2k/5)+10. Eventually

 D-d0 >= w-(2/5)*(11w/5)-11
        =3w/25-11 >= w/10.

Each nonzero binary linear combination of output coordinates is a nonzero linear functional pi of x^7. Trace pairing expresses pi(z)=Tr(mu*z) for a nonzero mu. The prime-cubic lemma supplies dispersion at d0 for this same cubic, and the cubic-bias lemma gives the bound

 |E[(-1)^(pi(x^7))]| < 4*2^(-(D-d0)/128)
                       <=4*2^(-w/1280).

For any distribution P on m bits, normalized Fourier/Parseval and Cauchy--Schwarz give

 TV(P,uniform) <=(1/2)*sqrt(sum_(a!=0) |E[(-1)^(a.P)]|^2).

To see normalization explicitly, write p(z)-2^-m=2^-m sum_(a!=0) bias_a*(-1)^(a.z), square and sum over z, then multiply by2^m and apply Cauchy--Schwarz to the L1 distance. Thus the output error is at most

 2*2^(m/2-w/1280) <=2*2^(-w/2000),

because m/2<=w/20000 and 1/1280-1/20000>1/2000. This is an exponentially small error in w, with output length Theta(w). The input affine source need not have independent halves. In particular it applies to concatenating two w-bit affine rows whenever either row is uniform, since joint entropy is then at least w.

This vector extractor supplies the seed inside the affine row condenser proved in extractors.md. Its use requires the conditional-affinity and entropy arguments supplied there.