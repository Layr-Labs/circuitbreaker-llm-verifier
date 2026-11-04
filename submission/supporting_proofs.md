# CircuitBreaker round 7: compact complete proof, part 2 of 2

This file is mandatory supporting mathematical evidence for proof.md.
The fixed total evaluator, actual-length theorem, native gate model, and
disjoint gate proofs are in Chapters I--IV of that file. Chapters V--VIII
below discharge every extraction and composition premise used there.

# Chapter V. Atomic linear strong seeded extractors

This note proves two uniform, input-linear strong seeded extractor constructions from elementary finite-field, coding, set-system, and counting arguments. It supplies the short-seed response component and a longer-seed, tiny-error component; it does not itself establish an affine disperser or a circuit lower bound.

The construction uses the coding-and-design reconstruction method associated with Trevisan, but every statement needed here is proved below. Primary provenance: Luca Trevisan, *Extractors and Pseudorandom Generators*, STOC 1999 / JACM 2001, https://lucatrevisan.github.io/pubs/extractor-full.pdf . No theorem is imported from that paper.

All logarithms are base two unless explicitly written `ln`. For distributions P,Q on a finite set, total variation distance is `TV(P,Q)=1/2 sum_z |P(z)-Q(z)|`, equivalently the largest difference in probabilities of a subset. A distribution X on n-bit strings has minentropy at least K when every atom has probability at most `2^(-K)`.

## 1. Main statements

### Theorem A: polynomial output and a logarithmic seed

Fix a rational constant `0<kappa<=1`. There is one deterministic polynomial-time algorithm constructing and evaluating maps

`Ext_n:{0,1}^n x {0,1}^{D_n} -> {0,1}^{m_n}`

with the following properties, for every integer `n>=2^(9/kappa)`:

* `m_n=floor(n^(kappa/3))`;
* `D_n <= (55300/kappa) ceil(log n)`;
* for every seed y, `x -> Ext_n(x,y)` is F2-linear;
* for every source X of minentropy at least `n^kappa`, independent of a uniform seed Y,
  `TV((Y,Ext_n(X,Y)),(Y,U_(m_n))) <= 1/8`.

The polynomial running-time exponent is an absolute constant once kappa is fixed, though it is large. The construction covers every integer n above the stated threshold. It can be defined arbitrarily, for example as zero, at smaller lengths while remaining total polynomial time.

The number of seeds is at most `n^(110600/kappa)` for n>=2. In particular, at `kappa=3/8`, the output is `floor(n^(1/8))` and the number of seeds is at most `n^300000`.

### Theorem B: tiny error and a constant fraction of the entropy

Fix constants `0<beta<kappa<=1`, with rational exponents for uniform integer parameter selection. Let `R=ceil(n^beta)`, let `d=ceil(1/beta)`, and let K be a given integer satisfying `n^kappa<=K<=n`. Put

`m=floor(K/(4*2^d))`.

For every sufficiently large integer n there is a deterministic polynomial-time, input-linear strong extractor with output length m, seed length at most

`2304 n^(2 beta)`,

and error at most `2^(-R)`. In particular its output length is at least

`K/(8*2^d)`

once `K>=8*2^d`.

The same conclusion applies when K is any prescribed polynomial lower bound with exponent kappa up to a positive constant factor: only the sufficiently-large threshold changes. More generally, Section 6 gives exact finite inequalities from which any desired parameter regime can be checked directly.

Theorem B is not claimed to use a logarithmic seed. Its seed is small enough for constructions that extract a seed of length `n^(2 beta)` from a separate component, and beta can be chosen as any sufficiently small positive fixed constant.

## 2. Deterministic field constructions needed by the code and the design

### 2.1 Fields whose Boolean dimension is a power of two

For every `s=2^j`, a field of size `2^s`, with an explicit F2-basis and deterministic arithmetic, can be constructed in time polynomial in s.

Start with F2. Suppose a field K of size `2^s` with a basis has been constructed. Every element z satisfies `z^(2^s)=z`: for nonzero z, multiplication by z permutes the nonzero field elements, so multiplying all of them and cancelling their nonzero product gives `z^(2^s-1)=1`.

Define

`Tr_s(z)=z+z^2+...+z^(2^(s-1))`.

This is F2-linear and takes values in F2, since its square equals itself. It is not the zero map: as a polynomial it is nonzero and has degree `2^(s-1)<2^s`, so the elementary root bound prevents it from vanishing at every field element. Consequently at least one basis vector has trace one. Find the first such basis vector c by computing the traces of all basis vectors.

The map `z -> z^2+z` is F2-linear with kernel `{0,1}`, so its image has `2^(s-1)` elements. Its image is contained in the kernel of `Tr_s`, by telescoping. The latter kernel has the same size because the trace is a nonzero linear functional. Thus the image is exactly the trace-zero subspace. In particular `Z^2+Z+c` has no root in K and is irreducible of degree two. The quotient

`K[Z]/(Z^2+Z+c)`

is a field of size `2^(2s)`, with basis given by the old basis and its multiples by the residue class of Z.

This constructs the next tower stage. To justify the algorithmic claim, multiplication of pairs of old field elements requires at most five multiplications in the preceding field and a linear number of bit additions. The recurrence `T(2s)<=5T(s)+O(s)` is polynomial in s. A trace requires s squarings and additions; scanning s basis vectors is polynomial. Summing over the doubling stages stays polynomial in the final dimension. This proof also supplies a deterministic, fixed choice of every extension coefficient.

The elementary root bound used here follows by repeatedly dividing by `X-a`: a nonzero polynomial of degree e has at most e distinct roots over a field. A quotient by an irreducible polynomial is a field because Euclid's algorithm gives a multiplicative inverse modulo the polynomial for every nonzero residue class.

### 2.2 Small fields of an arbitrary Boolean dimension

When `q=2^h` is polynomially bounded in the external input length, a field of size q can be constructed deterministically in polynomial time by enumerating monic degree-h binary polynomials and selecting the first irreducible one. To test a candidate, try division by every monic polynomial of degrees 1 through `floor(h/2)`. A reducible degree-h polynomial has a nonconstant factor of at most that degree. This uses at most `2^h` candidates and fewer than `2^(h/2+1)` divisors per candidate, hence `q^(3/2)` times a polynomial in h bit operations.

Here is an elementary existence proof so the search is not an unproved premise. Take a finite splitting field for `X^(2^h)-X`; it can be obtained by repeatedly adjoining a root of an irreducible factor. The polynomial has derivative one and hence `2^h` distinct roots. Its roots form a field under addition, multiplication and inversion, by the Frobenius identities, so this field has exactly `2^h` elements. Frobenius has its h-th power equal to the identity. If every element had orbit of size less than h, each orbit size would be a proper divisor of h and at most h/2. An element of orbit size e is a root of `X^(2^e)-X`. The union of all such root sets has at most

`sum_(e=1)^floor(h/2) 2^e < 2^(floor(h/2)+1) <= 2^h`

elements when h>=2, a contradiction. Therefore an orbit of size h exists. The product of the linear factors associated with that orbit is Frobenius-fixed, has coefficients in F2, and is an irreducible polynomial of degree h. The case h=1 is immediate.

## 3. A binary linear code with a proved list-size bound

Fix integers n,m>=1 and R>=1. Write `epsilon=2^(-R)`, `L=ceil(log n)` and `M=ceil(log m)`, with ceilings zero when the argument is one. Put

`B=3+L+2M+2R`,

and let s be the least power of two not smaller than B. Construct the field K of size `q=2^s` by Section 2.1. Then

`q >= 8 n m^2 / epsilon^2`,  and  `s<2B`.

Split the n input bits into `u=ceil(n/s)` coefficient blocks, padding the last block with zeros. Interpret these as elements `c_0,...,c_(u-1)` of K and write

`P_x(T)=sum_(j=0)^(u-1) c_j T^j`.

Define a binary codeword indexed by pairs `(a,b) in K x K`:

`C(x)_(a,b)=Tr_s(b P_x(a))`.

The codeword length is `q^2=2^(2s)`, so its coordinate index is exactly `a0=2s` bits. This length need not be materialized. One coordinate can be evaluated in time polynomial in n and s by Horner evaluation and trace. Each coordinate is F2-linear in x.

For distinct x,x', the difference polynomial is nonzero, with degree at most u-1. At a point a where that polynomial is nonzero, the function of b

`b -> Tr_s(b(P_x(a)+P_(x')(a)))`

is a nonzero linear functional and therefore balanced. At a root it is zero. Thus the normalized correlation of the two sign codewords is at most

`eta=(u-1)/q <= n/q <= epsilon^2/(8m^2)`.

Let an arbitrary binary received word have agreement at least

`1/2+gamma`,  where `gamma=epsilon/(2m)`,

with J distinct codewords. Normalize their sign vectors and the received sign vector to Euclidean norm one, writing the code vectors as v_1,...,v_J and the received vector as w. Each inner product `<w,v_i>` is at least `2 gamma`. Cauchy-Schwarz and the pairwise correlation bound give

`4 gamma^2 J^2 <= ||sum_i v_i||^2 <= J+eta J(J-1)`.

Hence

`J <= 1/(4 gamma^2-eta) <= (8/7)m^2/epsilon^2 < 2m^2/epsilon^2`.  (List)

This is a list-size statement only. No efficient list-decoding algorithm is needed for the extractor proof. All code evaluation remains deterministic polynomial time in n and s.

## 4. The reconstruction theorem, proved directly

Suppose sets `S_1,...,S_m` in a universe `[D]` each have exactly `a0=2s` elements. Give each set a fixed ordering. For a seed `y in {0,1}^D`, define

`Ext(x,y)_i=C(x)_(y|S_i)`.

For each i define the overlap-advice quantity

`A_i=sum_(j<i) 2^(|S_i intersect S_j|)`.

### Reconstruction theorem

If X has minentropy at least K and

`K >= max_i A_i + 1+3 log m+3R`,                       (Rec)

then, for uniform independent Y,

`TV((Y,Ext(X,Y)),(Y,U_m)) <= epsilon=2^(-R)`.

Proof. Suppose the distance exceeds epsilon. A Boolean test distinguishes the two distributions with probability difference greater than epsilon. Replace the m output bits one at a time by uniform bits. By telescoping, for some index i the test distinguishes the hybrid with the first i actual bits from the hybrid with the first i-1 actual bits by more than `epsilon/m`; complement the test if necessary so this difference is positive.

This gives a deterministic predictor for the i-th code bit from y and the previous i-1 code bits, succeeding with probability greater than `1/2+epsilon/m`. Here is an explicit derivation. Choose a uniform candidate bit b and uniform suffix bits; run the distinguishing test on that candidate. Predict b if the test accepts, and its complement otherwise. For the true next bit c, conditional success is

`1/2 + [Test(c)-Test(1-c)]/2`,

averaged over the suffix. Its bias is exactly the hybrid probability difference. Fix the candidate bit and suffix randomness to a choice attaining at least the average success. The resulting predictor P is a deterministic Boolean function of y and the true previous bits. It can be computationally arbitrary; the proof uses it only as a fixed mathematical function.

Write y as a0 bits v on S_i and bits z outside S_i. For any x,z, each earlier output bit j is a Boolean function of only the bits of v indexed by `S_i intersect S_j`; all its other seed bits are fixed by z. Record its complete truth table. The concatenation of these tables is advice `sigma(x,z)` of precisely A_i bits. Given fixed z and this advice, the predictor defines a binary word `W_(z,sigma)` indexed by v in `{0,1}^(a0)`.

Let `p(x,z)` be its agreement with C(x), averaged over uniform v. The predictor guarantee says

`E_(X,z) p(X,z) > 1/2+epsilon/m`.

Call `(x,z)` good when `p(x,z)>=1/2+epsilon/(2m)`. Since agreement is at most one, the probability of a good pair is greater than `epsilon/m`. Indeed, with `gamma=epsilon/(2m)` and good-pair probability p, the expectation is at most `1/2+gamma+(1/2-gamma)p`; exceeding `1/2+2gamma` forces `p>2gamma=epsilon/m`.

Fix one z for which the probability over X of goodness is greater than `epsilon/m`. There are at most `2^(A_i)` advice strings. For each one, (List) bounds the number of codewords agreeing enough with its received word by `2m^2/epsilon^2`. Each x has probability at most `2^(-K)`. Therefore

`Pr_X[good] <= 2^(A_i) * (2m^2/epsilon^2) * 2^(-K) <= epsilon/m`

by (Rec), a contradiction. This proves the theorem.

This proof counts one fixed predictor, index and outside-seed string. They are not encoded separately for each x. Only the earlier-bit truth tables depend on x, and all their bits are included in A_i. It does not assume X is flat or uniform.

For every fixed y the extractor is F2-linear in x because every selected code coordinate is linear. Thus the theorem establishes a strong **linear** seeded extractor, not merely an arbitrary seeded extractor.

## 5. A greedy logarithmic-universe set system and proof of Theorem A

Let `0<kappa<=1`, `n>=2^(9/kappa)`, `m=floor(n^(kappa/3))`, `R=3`, and construct the code above. Put `L=ceil(log n)` and

`d=ceil(kappa L/3)`,  `D=ceil(128 a0^2/d)`.

We have `1<=d<=a0` and `D>=2a0`. We construct m sets of size a0 with pairwise intersection strictly less than d by greedily choosing the lexicographically first eligible set from all a0-subsets of `[D]`.

To prove existence at each stage, let S be one already selected set and T a uniformly random a0-subset. A union bound over its d-subsets gives

`Pr[|S intersect T|>=d] <= binom(a0,d) * (2a0/D)^d <= (6a0^2/(Dd))^d < 16^(-d)`.

For the first inequality, a specified d-subset is contained in T with probability `(a0)_d/(D)_d`, bounded by `(2a0/D)^d` because `D>=2a0`. For the second, `binom(a0,d)<=(3a0/d)^d`: the inequality `ln(d!)>=integral_1^d ln t dt` yields `d!>=(d/e)^d`, and `e<3` follows directly from its power series.

The union of bad events from fewer than m prior sets has probability less than `m16^(-d)`. Since

`m<=n^(kappa/3)<=2^(kappa L/3)<=2^d`,

this is at most `2^(-3d)<1`. An eligible set therefore always exists. Intersections are at most d-1, so in particular

`A_i <= m2^d < 3n^(2kappa/3)`.

The last inequality uses `L<log n+1` and `d<kappa L/3+1`, giving `2^d<2^(1+kappa/3)n^(kappa/3)<3n^(kappa/3)`.

The reconstruction overhead is

`1+3log m+3R <= 10+kappa log n`.

Set `z=n^(kappa/3)>=8`. Then `n^kappa=z^3`, while the entire right side of (Rec) is less than

`3z^2+10+3log z <= 3z^2+10+3z <=4z^2<z^3`.

The middle inequality holds for z>=8. Thus minentropy `n^kappa` satisfies (Rec), proving error at most `1/8`.

For the seed estimate, write `M=ceil(log m)`. Since `L>=9` and kappa<=1,

`B=9+L+2M <=11+(1+2kappa/3)L <3L`.

The code field dimension obeys `s<2B<6L`, so `a0=2s<12L`. Also `d>=kappa L/3`. Consequently

`D<=55296 L/kappa+1 <=55300 L/kappa`.

Enumerating all subsets of `[D]` takes at most `2^D` steps. For fixed kappa this is polynomial in n, and greedy checking and code evaluation add only polynomial factors. More explicitly `L<=2log n` for n>=2, so

`2^D <= n^(110600/kappa)`.

All integer parameters can be computed by integer comparisons and binary search; for rational kappa, a statement such as `m=floor(n^(kappa/3))` is an exact integer-root condition. A small-length fallback makes one total algorithm. This completes Theorem A with no nonuniform choices.

## 6. Algebraic graph sets, exact finite inequalities, and Theorem B

For arbitrary n,K,m,R first construct the code, with index length a0. Let q0 be the least power of two at least a0, so `a0<=q0<2a0`. Construct the field F of size q0 by Section 2.2. Select the first a0 field elements in their fixed binary ordering as domain points.

Choose an integer d with `0<=d<a0` and `q0^(d+1)>=m`. Enumerate the first m distinct coefficient vectors of polynomials p over F of degree at most d. For each p use the set

`S_p={(u,p(u)): u is one of the a0 selected domain points}`

in the universe `F x F`, of size `D=q0^2<4a0^2`. Each set has size a0. Two distinct polynomials have at most d common graph points by the root bound, so

`A_i <=(m-1)2^d <=m2^d`.

Consequently the exact sufficient conditions for this strong linear extractor are

`d<a0`,
`q0^(d+1)>=m`,
`K >=m2^d+1+3log m+3R`.                              (Graph-Rec)

Under these conditions its output length is m, its seed length is less than `4a0^2`, and its error is at most `2^(-R)`. Its evaluation and construction take deterministic polynomial time in n, m and R when d is a fixed constant: the code field uses the polynomial-time tower construction, q0 is polynomial in `log n+log m+R`, and only the first m coefficient vectors and their a0 graph values are generated. The exponentially long binary code is never listed.

For Theorem B set `R=ceil(n^beta)`, `d=ceil(1/beta)`, and `m=floor(K/(4*2^d))`. For all sufficiently large n, `m>=1`, `d<a0`, and `m<=K<=n`. The code dimension is at least `2R`, so `a0>=4R>=4n^beta`. Therefore

`q0^(d+1)>=n^(beta(d+1))>=n>=m`.

The advice term is at most K/4. Since beta<kappa and `K>=n^kappa`,

`1+3log m+3R <=1+3log n+3ceil(n^beta) <=K/2`

for all sufficiently large n. This checks (Graph-Rec) with room to spare. These eventual inequalities are elementary: for positive alpha, `log n=o(n^alpha)` follows by writing `n=2^v` and using the exponential power series to dominate v; and `n^beta/n^kappa=n^(-(kappa-beta))` tends to zero. They hold at every sufficiently large integer n, not just on a subsequence. If an explicit threshold is desired, it can instead be defined as the first integer beyond which the displayed monotone domination bounds hold, or be replaced by any sufficiently large elementary power-of-two bound.

For the seed length, eventually `R>=5+3ceil(log n)`. Then

`B=3+ceil(log n)+2ceil(log m)+2R <=5+3ceil(log n)+2R <=3R`,

and hence `a0=2s<4B<=12R`. Therefore

`D<4a0^2<=576R^2<=2304n^(2beta)`.

Finally, once `K/(4*2^d)>=2`, the floor loses at most a factor two, so `m>=K/(8*2^d)`. This proves Theorem B. All constants depend only on the fixed exponents, and all construction choices are deterministic.

## 7. Exact surjectivity on affine sources from the short-seed construction

Let Ext be any input-linear strong seeded extractor with error strictly less than one half on a class containing a uniform affine source X on A. For a fixed seed y, the image of A under the linear map `x -> Ext(x,y)` is an affine subspace, and the output distribution is uniform on that image. If its rank is r_y<m, its total variation distance from the m-bit uniform distribution is

`1-2^(r_y-m)>=1/2`.

The strong total variation distance, including the independent uniform seed, is exactly the average of these conditional distances. Thus if every seed had deficient rank the distance would be at least one half, contradicting the extractor guarantee. At least one seed therefore restricts to a surjective affine map on A.

For Theorem A all seeds can be enumerated in polynomial time, giving a polynomial list of linear maps such that every affine source of the required dimension makes at least one list member exactly uniform on all m output bits. No efficient search for the successful seed on a given unknown source is asserted or required.



# Chapter VI. Cubic bias gives a constant flat of squared-logarithmic codimension

## Theorem

Let f:F_2^n -> F_2 have Boolean algebraic degree at most three and absolute bias beta>0. Then f is constant on an affine flat of codimension at most

    64 * (1+log_2(1/beta))^2.

All ingredients are proved below. This theorem supplies a stronger elementary bias upgrade than the inverse-polynomial upgrade in the original BK argument. It does not, by itself, construct a sublinear-threshold affine disperser or improve the benchmark's full-B2 coefficient.

Write b=log_2(1/beta). Complement f if necessary so that its signed bias is beta; this does not change its degree, derivative forms, or constant-flat dimensions.

## 1. Fourier and entropy preliminaries

For a real function g on F_2^n, put

    hat_g(xi)=E_x g(x)*(-1)^(xi dot x).

The characters are orthonormal, because the average of a nontrivial character is zero. Expanding sums therefore gives Fourier inversion and Parseval:

    g(x)=sum_xi hat_g(xi)*(-1)^(xi dot x),
    sum_xi hat_g(xi)^2=E_x g(x)^2.

For the normalized convolution

    (g*h)(z)=E_x g(x)*h(z+x),

direct substitution gives hat_(g*h)(xi)=hat_g(xi)*hat_h(xi).

We use Shannon entropy in bits. The chain rule follows by factoring joint probabilities into conditional probabilities and taking logarithms. Subadditivity follows from the chain rule and concavity of the binary entropy function

    H_2(p)=-p log_2 p-(1-p)log_2(1-p).

Indeed, its second derivative on (0,1) is negative, so averaging conditional probabilities cannot decrease entropy. Also, a distribution whose maximum point mass is at most 2^(-q) has entropy at least q, by averaging the pointwise inequality -log_2 p>=q. These elementary facts are the only entropy properties used below.

## 2. Dense sets contain a low-codimension subspace after few sums

**Lemma.** If A is a nonempty subset of F_2^n of density alpha and lambda=log_2(1/alpha), then there is a linear subspace V of codimension at most 8*lambda and an even positive integer m with

    lambda+3 <= m < lambda+5

such that every v in V is a sum of m members of A, with repetitions permitted.

**Proof.** Let p=1_A/alpha, so E p=1 and E p^2=1/alpha. Let

    Lambda={xi!=0 : |hat_p(xi)|>=1/2}.

Choose a linearly independent spanning subset xi_1,...,xi_d of Lambda. If X is uniform on A, the d-bit vector

    Y=(xi_1 dot X,...,xi_d dot X)

has maximum point mass at most 2^(n-d)/|A|=2^(-d)/alpha. Thus H(Y)>=d-lambda. Every coordinate of Y has bias of magnitude at least 1/2, and so its entropy is at most H_2(1/4). Moreover,

    H_2(1/4)=2-(3/4)log_2 3 <= 7/8,

because 3>=sqrt(8). Subadditivity now gives d-lambda<=7d/8, hence d<=8lambda.

Let V be the annihilator of the span of Lambda. Choose m to be the least even integer at least lambda+3; then m<lambda+5. For v in V, Fourier inversion gives

    p^(*m)(v)=sum_xi hat_p(xi)^m*(-1)^(xi dot v).

The zero-frequency term is 1. Every term from Lambda is nonnegative, because its character is 1 on V and m is even. The absolute value of the remaining sum is at most

    sum_(xi notin Lambda union {0}) |hat_p(xi)|^m
       <= 2^(-(m-2))*sum_xi hat_p(xi)^2
       = 2^(-(m-2))/alpha
       <= 1/2.

Consequently p^(*m)(v)>=1/2>0. Since p is nonnegative and supported on A, positive convolution value means exactly that v is a sum of m elements of A. QED.

## 3. Alternating forms of bounded rank have a common isotropic subspace

An alternating bilinear form B on a vector space E satisfies B(x,x)=0 for every x. In characteristic two it is symmetric. Its radical is

    rad(B)={x in E : B(x,y)=0 for every y in E},

and its rank is dim(E)-dim(rad(B)).

Every alternating form has a decomposition into nondegenerate two-dimensional summands and its radical. To see this, if B!=0 choose u,v with B(u,v)=1. The span of u,v is nondegenerate, and

    x -> x+B(x,v)u+B(x,u)v

projects onto its orthogonal complement. Repeat. In particular, the rank is even, and a form of rank d has a nondegenerate subspace of dimension d.

**Rank inequality.** For alternating forms P,Q on E,

    rank(P)+rank(P+Q)
       >= rank(Q)+2*rank(P restricted to rad(Q)).

**Proof.** Put K=rad(Q) and d=rank(P|K). Choose a d-dimensional subspace K0 of K on which P is nondegenerate. The decomposition

    E=K0 direct_sum K0^(perp_P)

exists because P|K0 is nondegenerate. In this decomposition P has block form D direct_sum P0, where D has rank d. Since K0 is in rad(Q), Q has block form 0 direct_sum Q0 and rank(Q0)=rank(Q). Thus

    rank(P)+rank(P+Q)
       =2d+rank(P0)+rank(P0+Q0)
       >=2d+rank(Q0),

where the final step is the elementary matrix inequality rank(A+B)<=rank(A)+rank(B). QED.

**Common-isotropic lemma.** Let L be an F_2-linear space of alternating forms on E, and suppose every member has rank at most R. Then there is W<=E of codimension at most 2R such that every form in L vanishes on W x W.

**Proof.** If every form is zero, take W=E. Otherwise choose Q in L of maximum rank r<=R and restrict every form to E1=rad(Q), losing r dimensions. For any P in L, both P and P+Q have rank at most r. The rank inequality therefore gives rank(P|E1)<=r/2. The restricted forms still form a linear space. Repeat on that space. At each nontrivial stage the maximum rank is at most half that of the preceding stage. The process terminates with all restricted forms zero and total codimension at most

    r+r/2+r/4+... <= 2r <= 2R.

QED.

## 4. Derivative polar forms of a cubic are linear in the direction

For a direction a define

    q_a(x)=f(x+a)+f(x).

Because f has degree at most three, q_a has degree at most two. Let B_a be its polar form:

    B_a(u,v)=q_a(u+v)+q_a(u)+q_a(v)+q_a(0).

It is an alternating bilinear form. The map a -> B_a is linear. This can be checked on algebraic-normal-form monomials: monomials of degree at most two contribute zero to B_a; for a cubic monomial x_i*x_j*x_k, the quadratic coefficients in q_a are a_i*x_j*x_k+a_j*x_i*x_k+a_k*x_i*x_j, which are linear in a. The associated trilinear form is the third derivative of f. In particular,

    B_(a1+...+am)=B_(a1)+...+B_(am).

## 5. Bias gives many low-rank derivative forms

Let h(x)=(-1)^f(x). For a uniform direction a put

    Z_a=E_x h(x+a)*h(x)=E_x(-1)^q_a(x).

Averaging over a,x gives E_a Z_a=(E h)^2=beta^2. Since Z_a<=1, the set

    A={a : Z_a>=beta^2/2}

has density alpha at least beta^2/(2-beta^2), and therefore at least beta^2/2. Hence

    lambda=log_2(1/alpha) <= 2b+1.

If q is a quadratic with nonzero absolute bias zeta, its polar rank is 2*log_2(1/zeta). Here is the full reason. The alternating decomposition from Section 3 puts its quadratic part into the form sum_(i=1)^r u_i*v_i. The remaining terms are affine, because their polar form is zero. Translate each pair to absorb its linear terms. Any nonzero remaining linear term on the radical would force zero bias. Thus it is absent when zeta!=0, and independent averaging over the r pairs gives zeta=2^(-r), while polar rank is 2r.

For a in A, the derivative has nonzero bias at least beta^2/2. Consequently

    rank(B_a) <= 2*log_2(2/beta^2)=4b+2.

Apply the dense-sumset lemma to A. It gives a subspace V of codimension at most 8lambda<=16b+8, and an even m<lambda+5<=2b+6, with V subset of mA. Linearity of the derivative forms and subadditivity of rank imply, for every v in V,

    rank(B_v) <= m*(4b+2)
       <= (2b+6)*(4b+2)
       =8b^2+28b+12=:R.

## 6. A large subspace on whose cosets f is quadratic

Restrict the linear space of forms {B_v:v in V} to the ambient vector space V itself. Their ranks are still at most R. By the common-isotropic lemma, there is W<=V with codimension within V at most 2R such that

    B_v(w1,w2)=0  for every v in V and w1,w2 in W.

In particular, since W<=V, the third derivative of f vanishes on W x W x W. It follows that the restriction of f to every affine coset of W has degree at most two. To justify this directly, parameterize a coset by any basis of W. Each coefficient of a degree-three monomial in the resulting Boolean polynomial equals the third derivative in the corresponding three basis directions, so every such coefficient is zero. Affine translation does not change the degree-three coefficients.

The ambient codimension of W is at most

    (16b+8)+2*(8b^2+28b+12)
       =16b^2+72b+32.

The signed bias beta of f is the average of its signed biases on the cosets of W. Thus at least one coset has signed bias at least beta. On that coset f is quadratic, and the quadratic normal-form argument above shows that fixing at most b+1 further affine coordinates makes it constant: fix one coordinate in each of the at most b quadratic pairs, then fix the value of the remaining affine function if it is nonconstant.

The total codimension is at most

    16b^2+73b+33 <= 64(b+1)^2.

This proves the theorem.

## 7. A fully proved error consequence

Suppose a cubic f has no constant affine flat of dimension greater than d0. Restrict it to any affine flat A of dimension D=d0+Delta, Delta>0, and let beta be its absolute bias there. If beta=0 there is nothing to prove. Otherwise the theorem implies

    Delta <= 64*(1+log_2(1/beta))^2.

Consequently

    beta <= min(1, 2^(1-sqrt(Delta)/8)).

No additional inverse theorem is used for this consequence.



# Chapter VII. Atomic cubic trace disperser and rate-half extractor

The cubic trace construction is due to Ben-Sasson and Kopparty, *Affine Dispersers from Subspace Polynomials*, Theorem 2.6 and Section 4.1. This file gives the complete binary coefficient proof, using interval notation with integer endpoints. The extractor upgrade uses the separately proved squared-logarithm cubic inverse theorem in Chapter VI.

## 1. Cubic disperser theorem

Let ell>=201 be odd, let K be any field of size q=2^ell, and let mu be any nonzero element of K. Then

    f_mu(x)=Tr_(K/F2)(mu*x^7)

is nonconstant on every F_2-affine flat in K of dimension greater than 2*ell/5+10.

All field identifications with Boolean coordinate spaces are F_2-linear. The Boolean algebraic degree of f_mu is at most three, because x^7=x*x^2*x^4 and both squaring maps are F_2-linear.

## 2. Linearized polynomials and image maps

A polynomial of the form sum_j a_j X^(2^j) is F_2-linear on K. A nonzero polynomial of degree D over a field has at most D roots: divide by X-a at each distinct root and induct. It follows that every polynomial function on K has a unique representative of degree less than q after reduction modulo X^q-X. In particular, a nonzero coefficient of a positive-degree reduced monomial certifies that the function is nonconstant.

For an F_2-linear subspace U<=K of dimension d define

    P_U(X)=product_(u in U)(X+u).

It is monic of degree 2^d and linearized. To prove the latter, induct on d. If U=U0 direct_sum F_2*b then

    P_U(X)=P_(U0)(X)*P_(U0)(X+b)
          =P_(U0)(X)^2+P_(U0)(b)*P_(U0)(X).

The roots are distinct, so the linear coefficient is nonzero. Put W=P_U(K), which has dimension ell-d, and let Q=P_W. The composite Q(P_U(X)) is monic of degree q and vanishes on K, so it equals X^q+X. Also P_U(Q(X))=X^q+X. Indeed, composing the difference of these last two polynomials with P_U gives zero: linearization and coefficients in K imply

    P_U(X^q+X)=P_U(X)^q+P_U(X).

A nonzero polynomial cannot become the zero polynomial after composition with a nonconstant polynomial, because degrees multiply. Thus the difference itself is zero.

Consequently P_U(Q(x))=0 for x in K, so Q(K)<=U. The kernel of Q is W, of dimension ell-d, hence its image has dimension d and equals U. Write

    Q(X)=sum_(j=0)^b a_j X^(2^j),    b=ell-d.

Its endpoints a_0,a_b are nonzero. Extend the coefficient sequence cyclically modulo ell, setting a_j=0 for b<j<ell. The number of roots of Q in K is 2^b.

For an affine flat gamma+U, the polynomial gamma+Q parametrizes it surjectively.

## 3. The unique long coefficient gap

Assume d>ell/3. The cyclic coefficient sequence of Q has no run of d consecutive zero coefficients, and its only run of d-1 consecutive zeros is the designated tail

    b+1,...,ell-1.

For the first statement, a Frobenius rotation of d consecutive zeros into the highest positions would give a nonzero polynomial of degree at most 2^(ell-d-1) with the same 2^(ell-d) roots as Q, contradicting the root bound.

For the second statement, suppose another run of d-1 zeros exists. A nontrivial Frobenius rotation moves it to the designated tail. The resulting reduced polynomial Q'=Q^(2^k) mod (X^q-X), with 0<k<ell, has degree at most 2^b and the same 2^b roots as Q. Comparing leading coefficients and subtracting gives Q'=cQ for a nonzero c. Thus every y in U satisfies y^(2^k)=c*y.

Choose beta in U nonzero. Then every z=y/beta satisfies z^(2^k)=z. Put g=gcd(k,ell). Frobenius is an invertible automorphism of K whose ell-th power is the identity. An integer Bezout combination of k and ell therefore gives z^(2^g)=z. This equation has at most 2^g roots, so dim_F2(U)<=g. Since ell is odd and g is a proper divisor, g<=ell/3, contradicting d>ell/3.

In particular, every zero gap strictly between two nonzero coefficients a_0,...,a_b has length at most d-2. The designated tail has length d-1 and is strictly longer.

## 4. The coefficient to be isolated

Let A=gamma+U be an affine flat of dimension greater than 2*ell/5+10. We may restrict to an affine subflat of dimension

    d=ceil(2*ell/5)+10.

This is possible under the hypothesis, including when 2*ell/5 is an integer. For ell>=201 we have d<ell/2, d>ell/3, and b=ell-d>d.

We will find a three-element subset B of {0,...,ell-1} such that the reduced polynomial

    R(X)=Tr(mu*(gamma+Q(X))^7) mod (X^q-X)

has a nonzero coefficient at X^(sum_(v in B)2^v).

A positive monomial obtained from a sum of h powers of two has reduced binary weight at most h. To see this, reduce positions modulo ell, then merge any repeated position into one unit at the next cyclic position. Each merge reduces the number of units and preserves the monomial function. The process terminates because the number of units decreases. A positive exponent never becomes exponent zero. Therefore any term of (gamma+Q)*(gamma^2+Q^2)*(gamma^4+Q^4) that uses a constant factor cannot contribute to a reduced exponent of binary weight three.

For a three-element column set B, define the matrix

    M(B)_(u,v)=a_(v-u)^(2^u),    u in {0,1,2}, v in B,

with coefficient indices taken modulo ell. Expansion gives the target coefficient exactly as

    sum_(tau=0)^(ell-1) mu^(2^tau)*Perm(M(B-tau))^(2^tau).       (1)

Indeed, to attain binary weight three the three chosen power positions must be distinct. Their assignments to the three Frobenius factors are precisely the permanent terms. Taking the trace rotates positions by tau and applies the corresponding Frobenius to the coefficient.

It is enough to choose B so that Perm(M(B))!=0 and Perm(M(B-tau))=0 for every nonzero tau. Multiplication by nonzero mu then leaves the tau=0 term nonzero, without any possible cancellation.

Let S={j in [0,b]:a_j!=0}. It contains 0 and b.

## 5. Case I: a nonzero coefficient in the middle interval

Set

    L=ell-2d+4,    U=d-4.

For ell>=201 and our chosen d, 0<L<=U<b. If there is j in S intersect [L,U], take

    B={0,j+1,b+2}.

At tau=0, column 0 has only its row-0 entry nonzero, and column b+2 has only its row-2 entry nonzero. The remaining row-1 entry of column j+1 is a_j^2!=0. Thus the permanent is the single nonzero product a_0*a_j^2*a_b^4.

For 1<=tau<=d-3, all three coefficient indices in column 0 of M(B-tau) lie in the designated zero tail, so that column is zero. For b+3<=tau<=ell-1, the same is true of column b+2.

The remaining shifts satisfy d-2<=tau<=b+2. In their middle column the coefficient indices are j+1-tau,j-tau,j-1-tau. The choices L<=j<=U give

    -d+1 <= j-1-(b+2),
    j+1-(d-2) <= -1.

Thus all three indices lie in [-d+1,-1], which is the designated tail modulo ell. The middle column is zero. This proves that every nonzero shift has zero permanent.

## 6. Case II: two short support intervals

Assume S has no element in [L,U]. Let

    a=max(S intersect [0,L-1]),
    c=min(S intersect [U+1,b]).

These exist because 0,b are in S. All nonzero coefficients lie in the two intervals [0,a] and [c,b]. Put

    A0=a+1,    C0=b-c+1,
    G=c-a-1,  T=ell-b-1=d-1.

Their endpoints are nonzero, the intervening gap of length G is zero, and the designated tail of length T is zero. Section 3 gives G<T. Also

    A0,C0 <= ell-2d+4,
    G >= 3d-ell-7.

Therefore

    G-max(A0,C0) >= 5d-2ell-11 >= 39.

The support of a nonzero column, viewed by its column index z, is contained in

    P=[0,A0+1] union [c,b+2]       modulo ell,

because its entries use coefficient indices z,z-1,z-2. The two intervals in P are separated in either cyclic direction by more than max(A0,C0)+1 positions: the smallest forward distances between the intervals are G-1 and T-1, both larger than max(A0,C0)+1.

Consequently, if two columns have index separation at most max(A0,C0)+1 and both are nonzero, their indices must lie in the same interval of P. This is only a necessary condition for nonzero columns, which is all we use; internal coefficient zeros cannot create additional possibilities.

### Case II.a: A0>C0

Choose

    B={0,b+1,a+2}={0,b+1,A0+1}.

At shift zero, column 0 forces row 0, column A0+1 forces row 2 with entry a_a^4, and column b+1 supplies the nonzero row-1 entry a_b^2. The permanent is nonzero.

For a general shift put z=-tau modulo ell. If the columns indexed by z and z+(A0+1) are both nonzero, they must lie in the same interval of P. The second interval has width C0+1<A0+1, so it cannot contain both. The first interval has width exactly A0+1, so the only possibility is z=0. Thus every nonzero shift has a zero column.

### Case II.b: A0<C0

Choose

    B={c,a+1,b+2}.

At shift zero, column c forces row 0 with entry a_c, column b+2 forces row 2 with entry a_b^4, and column a+1 supplies the nonzero row-1 entry a_a^2. The permanent is nonzero.

The first and last columns have separation C0+1. By the same interval argument, both can be nonzero only when their indices are exactly c and b+2, which forces tau=0. Every nonzero shift has a zero column.

### Case II.c: A0=C0

Choose

    B={0,c+1,a+2}={0,c+1,A0+1}.

At shift zero, the first and last columns force rows 0 and 2, with nonzero entries a_0 and a_a^4. The middle column has nonzero row-1 entry a_c^2, so the permanent is nonzero.

The first and last columns have separation A0+1. Since both support intervals have that exact width, both columns can be nonzero only when z=-tau is 0 or c. The case z=0 is the desired shift. In the exceptional case z=c, the shifted column set is

    {c,2c+1,b+2},

because b=c+A0-1. The columns c and b+2 again force rows 0 and 2. The only possible permanent term uses the row-1 middle entry a_(2c)^2. It is zero: from

    ell=2A0+G+T,    c=A0+G,    b=2A0+G-1,

and G<T, we obtain

    b<2c<ell.

Thus 2c lies in the designated zero tail. The exceptional permanent is zero, even though the middle column need not be entirely zero. All remaining nonzero shifts have a zero first or last column.

This finishes Case II and the coefficient isolation.

## 7. Completion of the cubic theorem

In every case equation (1) has exactly one nonzero summand, the one for tau=0. The reduced polynomial R therefore has a nonzero coefficient on a positive-degree monomial and is nonconstant on K. Since gamma+Q maps K onto the selected affine subflat of A, f_mu is nonconstant there and hence on A.

The proof works for every nonzero mu; it does not select mu depending on the flat. No character-sum estimate, unproved extractor theorem, or computation is a premise.

## 8. Quantitative rate-1/2 extractor

Combine Section 1 with the fully proved cubic inverse theorem from Chapter VI: a cubic of bias beta is constant on a flat of codimension at most 64(1+log_2(1/beta))^2. Therefore, if an affine source has dimension D and

    Delta = D-floor(2ell/5+10)>0,

every nonzero mu satisfies

    |E Tr-character(mu*x^7)| <= min(1,2^(1-sqrt(Delta)/8)).

Here the notation on the left means |E_x (-1)^(Tr(mu*x^7))|.

For odd ell>=2^16 and D>=(ell-1)/2, one has Delta>=ell/20. Choose any F_2-independent field elements mu_1,...,mu_M with

    M=floor(sqrt(ell)/64),

and output

    Q(x)=(Tr(mu_1*x^7),...,Tr(mu_M*x^7)).

Then Q is an affine extractor for dimension at least (ell-1)/2, with

    M >= sqrt(ell)/128,
    TV(Q(X),Uniform_M) <= 2^(-sqrt(ell)/64).

For the error bound, every nonzero output character has the bias estimate above. Fourier inversion, Parseval, and Cauchy--Schwarz give

    TV <= 2^(M/2-sqrt(Delta)/8)
       <= 2^(sqrt(ell)/128-sqrt(ell)/(8sqrt(20)))
       <= 2^(-sqrt(ell)/64).

The last inequality follows already from sqrt(20)<5 and 1/40-1/128>1/64. The output bound uses floor(z)>=z/2 for z>=2; our ell lower bound is more than sufficient.

This gives an elementary substitute for a rate-1/2 affine extractor with polynomial output and stretched-exponential error, after choosing an odd field length and padding a two-row slice by at most one zero coordinate. It does not yet give the entire Shaltiel construction or a benchmark family lower bound.

## 9. Evaluator scope

Given a monic irreducible modulus of odd degree ell and independent multipliers, evaluation is deterministic polynomial time: compute x^7 by squarings and two field multiplications, then compute the required traces. The multipliers may be the first M vectors of the polynomial basis.

For an all-length benchmark family, one must still specify how valid moduli are supplied or constructed and count all prefix bits if they are supplied as restrictions. That uniform-family accounting is not silently omitted here. The extractor theorem itself is field-isomorphism invariant and holds for every valid modulus.


# Chapter VIII. Direct recursive composition for a fixed reciprocal exponent

Fix an integer M>=2 once and for all. Use exactly Chapter I's parameter definitions and deterministic conventions; the submitted language fixes M=100. Write kappa=1/M, P=Pconf, q=q_tree, R=Rseed, and epsilon=2^(-ceil(n^beta)). The symbols ell (challenge row count) and h=2w+1 (field degree) are distinct. All logarithms are base two. An affine source is uniform on an affine flat; its entropy H is the flat's dimension. We prove that, for every sufficiently large integer n and every valid degree-h modulus, the specified payload evaluator is nonconstant on every affine flat of dimension at least K=ceil(n^kappa). Every asymptotic comparison below has fixed constants depending only on M. No parameter grows with n except those explicitly defined from n in Chapter I.

## 1. Atomic components and affine facts

Chapter V's component A at entropy exponent kappa/2 is input-linear for each seed, outputs L0=floor(n^(kappa/6)) bits, and has strong statistical error at most 1/8 on every source of entropy at least K0=ceil(n^(kappa/2)). Its R seeds can be enumerated in polynomial time, and

    R<=n^(221200M)<=n^Aexp.

On an affine source, a fixed-seed output is uniform on an affine image. If this image is proper, its distance from uniform is at least 1/2. Since strong error averages the fixed-seed distances, at least 3/4 of the seeds give exactly uniform output. In particular one such seed exists; truncating its output preserves uniformity.

Chapter V's component B, with overlap d=1000M and error epsilon, has seed length at most 2304*n^(2beta), is input-linear for each seed, and supports every output length z<=K'/(4*2^d) when K'>=n^(kappa/9), for sufficiently large n. The finite Bpars tests in Chapter I check the exact reconstruction conditions for every actual invocation. On an affine source of dimension at least K', at most a 2*epsilon fraction of uniform seeds fail to give exactly uniform output, by the same affine-image argument.

Chapters VI and VII supply the following cubic map Q for every valid odd-degree field representation of degree h=2w+1. On any affine source of dimension at least w it outputs floor(sqrt(h)/64) bits within distance 2^(-sqrt(h)/64) of uniform. A pair of w-bit affine slices, one uniform, has joint entropy at least w; append one zero bit and apply Q. Since w=ceil(n^(6beta)), its output length is Omega(n^(3beta)), longer than every B seed, and its error is at most epsilon, eventually. The same valid modulus serves every cubic invocation. The guards in Chapter I ensure the finite field-size, output-length, and error conditions whenever evaluation proceeds.

We will use the following elementary affine facts. A consistent conditioning on at most c linear equations loses at most c dimensions. This also bounds the loss of any projected entropy and of any conditional entropy of projected blocks: intersect the kernel of the earlier projection with the new equations. For a jointly affine pair (U,V), a linear section of the projection gives V=L(U)+Y, with L affine and Y an independent linear source of dimension H(V|U). Likewise, for a linear observation Z of X,

    X=a(Z)+Y,

where Y is a fixed linear source independent of Z, of dimension at least H(X)-length(Z). Conditioning further on a function of U alone preserves the affine conditional entropy H(V|U), because all nonempty U-fibers have the same direction kernel.

Deterministic maps contract statistical distance. If outside an event of probability eta a distribution is a mixture of members of a class, replacing its bad conditional distributions by arbitrary members makes it eta-close to such a mixture. We use this fact explicitly when nonlinear extraction prevents unconditional affinity.

## 2. Merger with deterministic lengths

Pad each input chunk to exactly G=floor(n^gamma) rows with zero rows. Suppose its affine matrix source has row length L0 and at least one uniform row. At step j there are T_j=ceil(G/2^j) rows of length L_j. Slice the first w bits of every row and call the entire collection Z. Pair the slices in order, using a zero partner for a final unmatched slice. Apply Q to each pair and retain the seed length of

    B_(n, floor(L_j/2), L_(j+1)),
    L_(j+1)=floor(L_j/(16*2^d)).

For each pair, use its seed to apply this B map to the entire current matrix, concatenated and padded with zeros to n bits; this produces one output row. This is exactly the halving operation specified in Chapter I.

Choose a pair containing a uniform current row. Its Q seed is epsilon-close to uniform. Fixing all slices, including every other pair's slices, gives X=a(Z)+Y with the same independent linear source Y on every fiber, and

    dim Y>=L_j-T_j*w>=L_j/2.

At most 2*epsilon of uniform seeds are bad for B on Y, so the chosen pair's actual seed is bad with probability at most 3*epsilon. Outside this bad slice event, all seeds are fixed and every output row is affine in the remaining source. The chosen pair's row is exactly uniform, since it is a constant plus B(Y,seed). Thus the unconditional output is within 3*epsilon of a mixture of affine matrices having a uniform row. It need not itself be affine.

Iterate this argument component by component in the mixture. Distance contraction adds at most 3*epsilon per step. After b=ceil(log G) steps, every good mixture component consists of a single uniform row, so the output is within 3b*epsilon of uniform.

Here are all the length checks, without a circular lower-bound assumption. Put Ddiv=16*2^d=2^(C-1). Repeated division by this integer gives the exact identity

    L_j=floor(L0/Ddiv^j).

Since b<=gamma*log n+1 and L0>=n^(kappa/6)/2 eventually, uniformly over j<=b,

    L_j>=2^(-C)*n^(kappa/6-(C-1)*gamma)-1
       >=2^(-C-1)*n^(kappa/6-C*gamma)
       =2^(-C-1)*n^(497*kappa/3000)
       >=n^(kappa/8)

eventually. The last two eventual comparisons use fixed positive exponent gaps. Further,

    G*L0<=n^(gamma+kappa/6)<n,
    T_j*w<=2*n^(gamma+6beta)=o(n^(kappa/8)),

because gamma+kappa/6<1 and gamma+6beta<kappa/8. Thus matrix padding is legal and every slice cost is at most floor(L_j/2) for large n. The B entropy parameter is at least n^(kappa/8)/3>=n^(kappa/9). Its output is at most floor(L_j/2)/(4*2^d) once L_j>=4. Every component hypothesis holds. The final output length exceeds 2304*n^(2beta), because kappa/8>2beta, and

    eta_M:=3b*epsilon<=2^(-n^beta/2)

eventually. These estimates apply to every padded chunk, including chunks originally containing only one row.

## 3. Extraction from B affine blocks

Consider B jointly affine raw blocks X_1,...,X_B, each of length at most n and each with conditional entropy at least K0 given its predecessors. Pad shorter inputs to n bits. Enumerate A's seeds on X_1; the resulting R rows of length L0 are jointly affine and have a uniform row. For each next block X_i, partition the current rows into consecutive chunks, pad each chunk to G, merge it, and truncate its output to the seed length of B_(n,K0,L0). Apply this B map to X_i with each resulting seed. The new rows all have length L0; their number is the previous number divided by G, rounded up.

The regeneration output is permitted because L0<=K0/(4*2^d) eventually, and the merger output supplies its seed. Since G>=n^(gamma/2) eventually and (B-1)*gamma/2>Aexp, after B-1 rounds the row count is

    ceil(R/G^(B-1))=1.

For completeness, the affinity induction conditions on exact raw prefixes, not on nonlinear row values. After processing block i, except on previous raw prefixes (X_1,...,X_(i-1)) of probability delta_i, the current matrix conditioned on that prefix is affine and has a uniform row. Initially delta_1=0. Fix a good prefix through i-2. One chunk of the current matrix has a uniform row; its merger seed, as a function of X_(i-1), is eta_M-close to uniform. In the jointly affine pair (X_(i-1),X_i) on this prefix, write

    X_i=L(X_(i-1))+Y,

where Y is independent and has dimension at least K0. The seed is bad for B on Y with probability at most eta_M+2*epsilon. For each exact X_(i-1) outside that event, every seed is fixed; the new matrix is affine in X_i and the designated row is exactly uniform. Hence

    delta_i<=delta_(i-1)+eta_M+2*epsilon.

At the last round, on every good prefix the single row is uniform. Its unconditional error is at most

    eta_BE:=(B-1)*(eta_M+2*epsilon)<=2^(-n^beta/3)

eventually. Denote this deterministic block extractor by BE. Its output length is L0. No independence of the raw blocks was assumed.

## 4. Tree, challenges, and the complete evaluator

Use the balanced t-ary coordinate tree of height H from Chapter I, with active depths 0,q,...,Jq. The exact parameter identities are

    tau*q=kappa/20,   tau*H=1-kappa/20.

Thus t^H<=2^H*n^(1-kappa/20)<=n eventually, every interval is nonempty, and each leaf has at most ceil(n/t^H)<=n^(kappa/20)+1 coordinates. The tree has O(H*t^H) nodes.

For each A seed j, truncate A_n(x,j) to ell*P bits and arrange them in an ell-by-P response matrix L_j(x). On a node s, pad its input x_s to n bits before applying L_j. For a descendant c of s and p<=P, set

    CR_p(x_s,c)=FAIL

exactly when the first p columns of the padded challenge ch_c equal the first p columns of some L_j(x_s); otherwise set PASS. FAIL is 0 and PASS is 1.

Compute challenges from the leaves upward. Leaves have empty challenges. At an internal node, choose the first child passing CR_P, or child 1 if none passes. Start the challenge with that child's challenge. Append the BE outputs on every increasing B-tuple of children, truncated to P bits, and pad this appended part to F=t^B rows. A depth-z node's challenge therefore has exactly (H-z)*F rows; pad it to ell rows in comparisons. Equal child row counts ensure that every inherited row has a fixed position.

Starting at the root, follow these child choices. For an active node s at depth u*q, let c be the next active node on this computed path and form the two tests CR_(p_u^(1))(x_s,c) and CR_(p_u^(2))(x_s,c). At the first active node whose first test fails, output its second test. Output zero if none fails. No test is made at the final active level. This is the evaluator in Chapter I; all proof-dependent choices below are absent from it.

We record the remaining parameter margins once. With

    e=2*tau*(B+1)=kappa^2/2000000000,
    pi=beta/4=kappa/4000,

we have ell=Ht^B<=n^e eventually. The explicit definitions give

    e+pi<kappa/6,
    tau+e+pi<beta,
    (120M+1)*e<pi.

Consequently ell*P<=L0 eventually, so responses and truncated BE rows fit. Also p_0^(2)=ell^(120M)<=P and

    ell*p_0^(2)=ell^(120M+1)=o(P).

The total length H*t*ell*P+(H+1)*ell*P is o(n^beta). The smallest active confidence is ell^9, which exceeds any fixed multiple of log n because ell>=H*n^(tau*B) grows as a fixed positive power. Finally eta_BE<=2^(-P) eventually. These are elementary fixed-power comparisons, valid at every sufficiently large integer n. The guards in Chapter I check finite component feasibility; the source and probability margins are mathematical consequences, not runtime oracles.

## 5. Two response bounds

Let Y be an affine source on a node s. If H(Y|Y_c)>=K0, then

    Pr[CR_p(Y,c)=FAIL]>=2^(-ell*p).

Indeed condition on any exact Y_c. Its challenge becomes fixed, and the conditional source on s has dimension at least K0. One enumerated A seed is surjective on ell*P bits, hence on the first p columns, and matches the fixed challenge with probability 2^(-ell*p). The event FAIL is the union over all seeds. The chosen seed may depend on Y_c; averaging preserves the lower bound.

Conversely, suppose that after every consistent conditioning on at most ell*p global linear equations, ch_c has a fixed row eta-close to uniform on P bits. Then

    Pr[CR_p(Y,c)=FAIL]<=R*(2^(-p)+eta).

For each response seed, condition on its first p columns, which imposes at most ell*p linear equations. The good challenge row matches its prescribed p bits with probability at most 2^(-p)+eta; the entire matrices cannot match more often. Average and sum over R seeds. A good row may depend on the conditional distribution, but is one fixed row within that distribution.

## 6. Nested structural selection

Let X be any affine source of dimension at least K. Call it nice at an interval s if all coordinates strictly before s are fixed and H(X_s) exceeds the stated threshold.

First find consecutive active nodes g,g'. Start at the root with threshold K. At an active node with threshold k', partition it into its T=t^q next-active descendants and let k_i be their successive conditional entropies. Choose the first i with k_i>=k'/(4T); such an i exists because their sum is H(X_s)>=k'. Fix the preceding blocks to a possible value. Their total conditional entropy is less than k'/4, so the parent loses at most that many dimensions and the chosen child's entropy is at least k'/(4T). If its entropy is at least k'/4, continue on it with threshold k'/4. Otherwise the parent's conditional entropy outside this child is at least k'/2. Stop, with g the parent and g' this child.

Continued descent to the leaves would retain entropy at least K/4^J, larger than leaf length, so stopping occurs. At the stopping pair, both child entropy and H(X_g|X_g') are at least

    K/(4^(J+1)*t^q)>=c_M*n^(19*kappa/20)>=n^(9*kappa/10)

eventually, for a positive constant c_M. In particular g' is internal, and all coordinates before it are fixed.

Next descend inside g' to expose B affine blocks. Start with threshold n^(9*kappa/10). At a current node of entropy at least k', inspect the t children's successive conditional entropies. If at least B of them are at least k'/(4t), stop at e and select an increasing B-tuple of such children. Their conditional entropy given only earlier selected children is at least the conditional entropy given all earlier children, so they form the required block source.

If fewer than B children meet that threshold, all small children together account for less than k'/4. Some remaining child has conditional entropy at least k'/(4B). Fix all preceding children and continue on it with threshold k'/(4B). These fixings remain inside g', hence preserve H(X_g|X_g') exactly by the common-kernel fact. They maintain niceness. This descent also stops before a leaf: after H steps the retained entropy would still be at least n^(9*kappa/10)/(4B)^H, larger than leaf length. At the stopping node each selected block has conditional entropy at least

    n^(9*kappa/10)/((4B)^H*4t)>=n^(3*kappa/4)

eventually; the inequality follows from 9kappa/10-tau>3kappa/4. Let X0 be the resulting nonempty affine subsource. All coordinates before e are fixed. Its selected B children have conditional entropy at least T0=n^(3*kappa/4), and H((X0)_g|(X0)_g')>=T0. All restrictions used here are nested.

## 7. Force preceding empty children to fail

For every proper ancestor s of e, process every child c preceding the child on the path to e. All c's input coordinates were fixed, so its padded challenge is a constant ell-by-P matrix. On the current source the projection to s retains entropy at least T0 minus previous restriction costs, which is at least K0. Choose a response surjective on ell*P bits and restrict it to that constant challenge. This is a consistent affine restriction of codimension at most ell*P and makes CR_P(x_s,c)=FAIL identically. Earlier failures persist under later restrictions.

There are at most H*t operations, costing at most D0=H*t*ell*P. Denote the final source by X'. Each selected block's conditional entropy and the outside entropy at g lose at most D0. After any further (H+1)*ell*P global linear equations they still exceed K0, because T0=n^(3*kappa/4), K0=ceil(n^(kappa/2)), and the total equation loss is o(n^beta). This is one nested affine history. Its restrictions are for proving dispersal; they are not additional circuit-deletion charges.

## 8. Resilient challenges and navigation

For an ancestor s of e lying a edges above e, including e at a=0, we claim that after every consistent conditioning of X' on at most (H-a+1)*ell*P global linear equations, ch_s has a fixed row within distance

    eta_a=(2*(R+1)^a-1)*2^(-P)

of uniform. At a=0, the chosen B children retain entropy K0 under all allowed equations, so their BE row has error at most eta_BE<=2^(-P). Its appended position is fixed because every child challenge has the same length.

For a>=1, let c be the child of s leading to e. One extra response conditioning costs at most ell*P, within c's induction budget. The upper response bound gives

    Pr[CR_P(x_s,c)=FAIL]<=R*(2^(-P)+eta_(a-1)).

Every preceding child fails identically, even after the conditioning. Thus s chooses c except on this event. When it does, ch_s begins with exactly ch_c. A good row of ch_c has a fixed position in that prefix because all child challenges have equal lengths. Coupling the two row values changes their distributions by at most the exceptional probability. The resulting error is at most

    eta_(a-1)+R*(2^(-P)+eta_(a-1))
       =(2*(R+1)^a-1)*2^(-P)=eta_a.

This proves the claim with its full conditioning budget. Since H is fixed and R is polynomial in n,

    eta_a<=(R+1)^(a+1)*2^(-P)<=2^(-P/2)

eventually, uniformly for a<=H. Apply the upper response bound to each intended edge on the fixed root-to-g' path and take a union bound. The earlier-child forced failures ensure that a passing intended child is selected. Therefore

    Pr[computed path does not visit g']<=2^(-P/4)

for large n. This uses the resilience statement on the actual affine source X'; it does not assert unconditional affinity of nonlinear challenges.

## 9. Both output values occur

Suppose g has active depth u*q. Define fixed-path tests r_a(x)=CR_(p_u^(a))(x_g,g'), and abbreviate p1=p_u^(1), p2=p_u^(2). They equal g's actual tests whenever the computed path visits g'. The outside entropy at g is at least K0, so the lower response bound gives

    Pr[r_a=FAIL]>=2^(-ell*p_a).

The resilient challenge at g' and the upper response bound give

    Pr[r_2=FAIL]<=R*(2^(-p2)+2^(-P/2)).

Because p2=ell^3*p1, ell*p2=o(P), and p1 dominates log R, this is at most (1/2)*2^(-ell*p1) for large n. Consequently

    Pr[r_1=FAIL,r_2=PASS]>=(1/2)*2^(-ell*p1).

FAIL at confidence p2 implies FAIL at p1, so also

    Pr[r_1=FAIL,r_2=FAIL]>=2^(-ell*p2).

Both displayed desired events have probability at least 2^(-ell*p2) eventually. For every earlier active level v<u, use its fixed intended next-active descendant on the path to g. That descendant is an ancestor of e, so its first test fails with probability at most

    R*(2^(-p_v^(1))+2^(-P/2)).

If u>0 the smallest such earlier confidence is p_(u-1)^(1)=ell^6*p1, whereas ell*p2=ell^4*p1. Thus the union over at most J earlier active tests is o(2^(-ell*p2)). If u=0 this bad event is empty. The navigation error 2^(-P/4) is also o(2^(-ell*p2)), since ell*p2<=ell*p_0^(2)=o(P). These comparisons absorb all polynomial factors R and the fixed factors H,J: ell^9 dominates log n and the indicated exponent differences tend to infinity as positive powers of n.

Subtract the two bad probabilities from either desired event. Positive probability remains that the computed path reaches g', every earlier first test passes, g's first test fails, and its second test has the chosen Boolean value. The evaluator then selects g and returns that value. Both outputs occur on X', and therefore on the original affine source X.

This proves ordinary affine dispersal at threshold ceil(n^(1/M)) for every sufficiently large integer n and every valid degree-h modulus. The evaluator is deterministic and polynomial time for each fixed M: the tree depth and block count are constants, its numbers of nodes, tuples, responses, and field-coordinate operations are bounded by fixed powers of n, and every component is uniform. Chapter I gives its total small-length convention and counted all-length wrapper. Applying Chapters II and III exactly as in Chapter I yields, for the permanently fixed M=100 language,

    C_B2(F_N)>=4N-8*N^(1/100)-6*N^(3/50000)-23

at every sufficiently large integer N. More generally each separately fixed M>=2 gives 4N-9*N^(1/M) eventually. This does not allow M to vary with N, and proves neither an n^o(1) error for one fixed language nor a leading coefficient above 4.

Provenance: the merger and affine residual-source method follow Anup Rao,
*Extractors for Low-Weight Affine Sources*, author manuscript,
https://homes.cs.washington.edu/~anuprao/pubs/lowweight.pdf . Recursive
challenge-response follows Ronen Shaltiel, *Dispersers for affine sources
with sub-polynomial entropy*, ECCC TR11-127 (2011),
https://eccc.weizmann.ac.il/report/2011/127/download/ . Every merger,
conditioning, structural-selection, and recursive statement used here is
proved above for the displayed parameters. These citations provide credit;
no external theorem is an unproved premise, and the stronger published
subpolynomial-entropy statement is not claimed here.
