# A polynomial-entropy affine disperser, with all component proofs

This file proves the construction used by family.py at ambient length N.
The components and their parameters are supplied in full. Fix STRENGTH=1024,
write sigma=STRENGTH and delta=1/sigma throughout. These are constants,
not input-dependent parameters. Mathematical
provenance: BK for the cubic/image polynomials; Trevisan and NW for seeded
code-bit extraction; Rao for affine somewhere-random condensation; Shaltiel
TR11-127 for challenge-response. All refinements below are proved directly.

## Finite-field and image-polynomial background

All vector dimensions are over F_2. Let K have size 2^n. A polynomial
function on K has a unique representative of degree below 2^n, by division
by X^(2^n)-X and the root bound. Linearized polynomials sum_j a_j X^(2^j)
are exactly the F_2-linear maps K->K: their coefficient space and the
linear-map space both have dimension n^2, and the root bound gives
injectivity of the map from polynomials to functions.

Tr(z)=sum_(j=0)^(n-1)z^(2^j) is a nonzero map to F_2. It is nonzero
by the root bound, since its nonzero polynomial has degree 2^(n-1)<2^n.
For nonzero mu, z->Tr(mu*z) is nonzero by bijectivity of multiplication.
Thus mu->(z->Tr(mu*z)) is an injective map between two n-dimensional
spaces, and every nonzero functional pi has that form with mu nonzero.

For a d-dimensional linear subspace V, let P_V(X)=product_(v in V)(X-v).
It is monic and linearized. Induct from P_{0}=X: adjoining b gives
P_(U+<b>)=P_U(X)^2+P_U(b)P_U(X). Its kernel is V and its image V' has
dimension n-d. Put Q=P_(V'). Then Q(P_V(X))=X^(2^n)-X: it is monic,
has that degree and vanishes on K. Also P_V(Q(X))=X^(2^n)-X. To see
this, compose the proposed difference with P_V and use
P_V(X^(2^n)-X)=P_V(X)^(2^n)-P_V(X); substitution by a nonconstant
polynomial is injective. Therefore Q(K) is contained in V, and dimension
counting shows Q(K)=V. It has degree 2^(n-d), is monic, and has distinct
roots, so its coefficient a_0 is nonzero. Its coefficient a_(n-d)=1 and
all later coefficients vanish.

There is no cyclic run of d zero coefficients of Q: rotating such a run
to the final positions by a Frobenius power and reducing modulo
X^(2^n)-X would give a nonzero polynomial of degree at most
2^(n-d-1) with all 2^(n-d) roots of Q, contradicting the root bound.

For a reduced exponent its binary weight is the number of one bits.
Products of h power-of-two factors have reduced weight at most h:
repeated residues carry, replacing two factors by one doubled factor,
with cyclic wrap at n. Every carry decreases the number of factors.
For h<n, weight h therefore requires h distinct residues with no carries.
Frobenius rotates these residues. These observations justify the
weight-three coefficient calculation in the cubic proof below.

## Statement

For a degree-at-most-three polynomial f on F_2^n, let
delta=|E_x (-1)^f(x)|>0 and L=ceil(log_2(2/delta)). Then f is constant on
an affine subspace of dimension at least n-512 L^4.

The bound may be vacuous. It is only its logarithmic dependence on bias that
matters. All ranks and dimensions below are over F_2.

## Quadratic preliminaries

For a quadratic q define its alternating polar form
B_q(u,v)=q(u+v)+q(u)+q(v)+q(0). Repeatedly choosing a pair on which B_q=1,
and taking their common orthogonal complement, gives a symplectic basis:
q is a sum of r disjoint products plus an affine linear function, and B_q
has rank 2r. Summing the sign over the radical either gives zero or leaves
absolute bias 2^-r. Fixing one coordinate in each product leaves an affine
linear function; fixing at most one additional coordinate makes it constant.
Consequently nonzero bias at least eta implies rank <=2 log_2(1/eta)
and a constant affine subspace of codimension at most log_2(1/eta)+1.

## A many-summands density lemma

If A is a subset of F_2^n of density mu>0, there is a linear subspace V of
codimension at most ceil(log_(3/2)(1/mu)) contained in 2k A, where
k=ceil(log_2(1/mu)/2)+2. Here 2k A means sums of exactly 2k elements,
allowing repetitions.

Start with the ambient space and A. If a nontrivial normalized Fourier
coefficient has absolute value at least mu/2, one of its two hyperplane
fibers has relative density at least (3/2)mu. Restrict to that fiber,
translate it to a linear subspace, and repeat. At most the stated number
of restrictions is possible. Call the resulting density nu>=mu and
translated subset A' of the resulting linear space V. All nontrivial
Fourier coefficients of its indicator have absolute value <nu/2.

The normalized 2k-fold convolution at any x in V has constant Fourier
term nu^(2k). By Parseval, the sum of the absolute values of its remaining
terms is at most

    (nu/2)^(2k-2) sum_gamma hat(1_A')(gamma)^2
    <= nu^(2k-1)/2^(2k-2) < nu^(2k).

The last inequality follows from 2^(2k-2)>1/mu>=1/nu. Thus the convolution
is positive at every x. Translating back causes no shift in an even sum,
so V is contained in 2k A as asserted. Each density-increment threshold
uses the CURRENT relative density. The integer k is chosen from the
original density, and remains sufficient because the final density is
at least the original one.

## A bounded-rank alternating matrix space lemma

Let M be a linear space of alternating forms on U, every member having
rank at most 2r. There is W<=U of codimension at most 2r^2 such that every
form in M vanishes on W x W.

Greedily start V_0=U. If some B_i in M is nonzero on V_(i-1), choose
u_i,v_i in V_(i-1) with B_i(u_i,v_i)=1, and set
V_i=V_(i-1) intersect Rad_U(B_i). Each step costs at most 2r dimensions.
There cannot be t=r+1 steps. To see this, restrict the formal form
sum_(i=1)^t a_i B_i to the ordered list u_1,v_1,...,u_t,v_t. Its Pfaffian
is a homogeneous polynomial of degree t. Its coefficient of a_1...a_t
is 1: B_1 vanishes on every later chosen vector, so a term using a_1
must match the first pair to itself; removing that pair gives the same
argument for B_2 and so on. This is the unique such contribution.

Reducing this polynomial modulo a_i^2-a_i cannot cancel that coefficient:
every monomial has total degree t, and a monomial that reduces to the
product of all t variables must already be exactly that product. A
nonzero multilinear polynomial is nonzero at some Boolean point (prove
by induction on its last variable). At that point the Pfaffian is 1,
so the restricted alternating form is nonsingular and has rank 2t.
This contradicts the rank bound for M. The identity determinant=Pfaffian^2
follows directly from the determinant expansion: permutations with a
cycle of length greater than two cancel with their reversed cycle in
characteristic two, leaving the matchings. A nonzero determinant on the
chosen list also proves the list is independent, so no independence
assumption was hidden in the construction.

Therefore the process stops after at most r steps and gives W as claimed.

## Applying the lemmas to a cubic

Write D_h f(x)=f(x+h)+f(x), and b_h=E_x (-1)^D_h f(x). The identity
E_h b_h=delta^2 follows by averaging over the pair (x,x+h). Hence the
set A={h:b_h>=delta^2/2} has density at least delta^2/2: otherwise even
bounding b_h by 1 on A and by delta^2/2 elsewhere makes its average too
small. Every D_h f for h in A is quadratic with polar rank at most
2 log_2(2/delta^2)<=4L.

The map h -> B_(D_h f) is linear: in the unique multilinear polynomial
representation its entries are the coefficients of the cubic terms
contracted against h. Equivalently it is the trilinear third derivative
T_f(h,u,v). Apply the density lemma to A. Its density mu is at least
2^(1-2L), so the resulting subspace V has codimension at most 4L and
every h in V is a sum of 2k elements of A with k<=L+2. Rank subadditivity
therefore bounds every B_(D_h f) by (2L+4)*4L<=24L^2.

Restrict these forms to V. They are a linear matrix space of rank at
most 24L^2. Apply the matrix lemma with r=12L^2 to obtain W<=V of
codimension at most 288L^4 in V, on which every such form vanishes.
In particular T_f vanishes on W^3. Thus f has degree at most two on
every coset of W; translation does not alter its highest cubic part.

The signed biases over the cosets average to the signed global bias,
so some coset has absolute bias at least delta. The quadratic argument
then gives a constant affine restriction within it of additional
codimension at most L+1. Total codimension is at most
288L^4+5L+1<=294L^4<=512L^4.

## Cubic extractor corollary

If a degree-three f on F_2^m is nonconstant on every affine space of
dimension greater than .45m+21, then on every affine source of dimension
at least m/2 its bias is less than 2^(-m^(1/8)) for all sufficiently
large m. Otherwise L<=m^(1/8)+2 and the lemma supplies a constant space
of dimension at least m/2-512(m^(1/8)+2)^4>.45m+21, a contradiction.

If the same disperser statement holds for every nonzero linear output
test of a polynomial map, projecting that map to floor(m^(1/16)) output
coordinates yields statistical distance at most
2^(m^(1/16)/2) * 2^(-m^(1/8)) / 2, by Parseval and Cauchy-Schwarz.
This is stretched exponentially small. This output is sublinear; the
published Rao parameters cannot be used verbatim.


## Cubic disperser by relaxed coefficient isolation

Let K=F_(2^n), n odd and sufficiently large, and pi:K->F_2 nonzero linear.
Then f(x)=pi(x^7) is nonconstant on every affine subspace of dimension
greater than 9n/20+21.

## Algebra already available in the preliminary coefficient expansion

Write pi(x)=Tr(mu*x) with mu nonzero. For an affine space A of dimension d,
write A=c+V and construct its image polynomial Q with Q(K)=V as in the
preliminaries above. It is linearized of degree 2^(n-d), is monic,
and its kernel has dimension n-d. Write q=n-d and
Q(X)=sum_j a_j X^(2^j). We have a_0 and a_q nonzero and a_j=0 for
q+1<=j<=n-1. No cyclic zero run has length d.

There is no other zero run of length d-1 when d>n/3. Otherwise a nonzero
Frobenius rotation of Q has the same roots and degree at most 2^q and
therefore equals a nonzero scalar times Q. Its image V then has
z^(2^h)=c_0*z for some 0<h<n and every z in V. Dividing by any nonzero
z_0 in V makes every ratio fixed by Frobenius^h and Frobenius^n, hence by
Frobenius^g for g=gcd(h,n). The latter implication follows from Bezout
on powers of the invertible Frobenius automorphism. The polynomial
X^(2^g)-X has at most 2^g roots, so dim(V)<=g. Since n is odd, every
proper divisor g of n is at most n/3, a contradiction.

For a three-element column set B, the coefficient of X^(sum_(b in B)2^b)
in Tr(mu*(c+Q(X))^7), reduced modulo X^(2^n)-X, is

    sum_(r=0)^(n-1) mu^(2^r)
      * perm( ( a_(b-r-u)^(2^u) )_(u=0,1,2; b in B) )^(2^r).

All indices are modulo n. This is the same expansion as the independent
source coefficient formula in the preliminary coefficient expansion. Three distinct
binary exponents cannot be produced by terms using c, or by a carry,
because those terms have binary weight at most two. To prove nonzero
it suffices to find B whose r=0 pattern is triangular with nonzero
diagonal and whose other permanents vanish.

## Pattern isolation with slack

By restriction take d=ceil(9n/20)+20. It is enough to consider sufficiently
large n, for example n>=10000. Then d>n/3, q=n-d, and all the inequalities
below have ample integer slack. All coefficient intervals are ordinary
integer intervals within [0,n-1] unless explicitly called cyclic.

Case 1: some a_j is nonzero with n-2d+4<=j<=d-4.
Take B={0,j+1,q+2}. At r=0 the first column has only its top entry
nonzero, the last column has only its bottom entry nonzero, and the
middle column's center entry is nonzero. Thus its permanent is nonzero.
For shifts 1<=r<=d-3 the first column is zero; for q+3<=r<=n-1 the
last column is zero. For the remaining d-2<=r<=q+2 all indices of the
middle column lie between -d+1 and -1, by the bounds on j, and hence
are in the terminal zero run modulo n. Every other permanent is zero.

Case 2: that interval contains no nonzero coefficient. Let A be the last
nonzero coefficient index before it and B_0 the first nonzero index after
it. Put alpha=A+1 and beta=q-B_0+1. The support of the coefficients is
contained in two clumps [0,alpha-1] and [B_0,q], whose endpoints are
nonzero and whose lengths alpha,beta are at most n-2d+4. The two gaps
between these clumps are entirely zero. Both gap lengths exceed
max(alpha,beta)+5: the internal gap is at least 3d-n-7 and the terminal
gap has length d-1. These bounds follow from d>=9n/20+20.

By reflection of coefficient indices assume alpha>=beta. This loses
nothing: replacing a_j by a_(q-j) converts a pattern indexed by column b,
shift r and row u into the original pattern at column q+2-b, shift -r
and row 2-u. Thus zero columns and triangular nonzero permanents transfer
back under reflection, regardless of actual Frobenius powers of entries.

Choose first and last columns 0 and alpha+1. At shift r, their coefficient
windows are [x-2,x] and [x+alpha-1,x+alpha+1], where x=-r modulo n.
Because both zero gaps are longer than alpha+5, these windows cannot
hit different clumps. Within a clump [c,c+L-1], both windows can hit it
only if c<=x<=c+L-alpha. For the first clump L=alpha, this forces x=0.
For the second L=beta, it is impossible if alpha>beta, and forces
x=B_0 if alpha=beta. Holes inside clumps can only eliminate possibilities.

If alpha>beta, choose the remaining middle column q+1. At shift zero
its center entry a_q is nonzero, the first column only has a_0 at its
top, and the last column only has a_(alpha-1) at its bottom. Thus the
permanent is nonzero; at every nonzero shift one of the first/last
columns is zero by the preceding paragraph.

If alpha=beta, choose the middle column B_0+1. Its center entry is
a_(B_0), and the r=0 permanent is again triangular and nonzero.
The only possible other shift is r=-B_0. There the first column has
only its top entry a_(B_0) nonzero, and the last column has only its
bottom entry a_q nonzero. Its permanent therefore uses the middle
entry a_(2B_0 mod n). This entry is zero. Indeed the internal zero gap
has length B_0-alpha, strictly smaller than d-1 by the unique-long-gap
property. Using B_0=q-alpha+1 gives

    2alpha >= n-2d+3, hence 2B_0<=n-1.

Also B_0>=d-3, so 2B_0>=2d-6>q for these parameters. Consequently
q+1<=2B_0<=n-1, which is in the terminal zero run. This last permanent
vanishes even though it need not have a zero column.

All cases give a nonzero weight-three coefficient and prove nonconstancy.
The construction is polynomial-time once the field is supplied, since
it uses a fixed power 7 and a linear output functional.

## Consequence with the separate bias lemma

For a source of dimension at least n/2, bias >=2^(-n^(1/8)) would yield,
by the cubic restriction lemma, a constant space of dimension at least
n/2-512(n^(1/8)+2)^4>9n/20+21 eventually. This contradicts the disperser
statement. All nonzero output tests of a projected x^7 have this property,
so outputting floor(n^(1/16)) fixed field coordinates gives a half-rate
affine extractor with stretched exponential error by the XOR/Parseval bound.


## A binary linear code with a small combinatorial list

Given n data bits, group them into w-bit field elements (zero-pad the last),
interpret them as coefficients of a polynomial P of degree below h=ceil(n/w)
over F_(2^w), and index code coordinates by pairs (a,b) of field elements.
Output the binary dot product b . P(a), in a fixed polynomial basis. This is
linear in the original n bits. Its codeword length is 2^(2w), and a code bit
can be computed in poly(n,w) bit operations if the field is supplied.

For distinct data words, P-P' is nonzero and has at most h-1 roots. Averaging
the sign of b.(P(a)-P'(a)) first over b shows that the normalized inner product
of their sign codewords is between 0 and (h-1)/2^w. Choose 2^w>=h/gamma^2.

Any fixed sign word has correlation at least 2gamma with at most
1/(3gamma^2) codewords. Indeed, for M such codewords v_i and target v,
with normalized Euclidean norms 1, Cauchy-Schwarz gives

    (2gamma M)^2 <= ||sum_i v_i||^2 <= M+gamma^2 M(M-1).

Thus M <= (1-gamma^2)/(3gamma^2). No list-decoding algorithm is needed:
only this counting bound enters the extractor proof. Index length ell=2w.

## Reconstruction estimate for an arbitrary design

Let S_1,...,S_m be ell-element subsets of [D] with pairwise intersection
at most a. For uniform seed y in F_2^D, output

    E(x,y)_i = Code(x)[ y restricted to S_i ].

This is linear in x for each seed. A sufficient condition for it to be a
strong (k,epsilon) extractor is the following conservative inequality:

    k > D + m*2^a + log_2(m) + 4 log_2(8m/epsilon) + 10,

where the code is set up with gamma=epsilon/(8m).

Proof: if a distinguisher separates
(Y,E(X,Y)) from (Y,U_m) by more than epsilon, hybridizing the output bits
gives a fixed index i and a predictor for its true code bit, from the seed
and the preceding true bits, with average advantage at least epsilon/m.
The hybrid gap equals one half of E[true_sign*Delta], so a prediction
with expected sign Delta has correctness advantage equal to that gap.
Formally let Delta be the difference between the distinguisher's two
candidate-bit values, averaged over the remaining independent uniform
output bits, oriented so the hybrid difference is positive. Use a randomized
prediction whose expected sign is Delta. Its correlation with the true sign
is the required hybrid correlation. Assign its random choices once for every
possible seed/prefix argument and fix the resulting finite random table by
averaging. Some single fixed deterministic predictor retains the average
correlation. This predictor may depend on the source distribution, but is
then common to all x; neither it nor its fixed random table is charged as
per-x advice. Taking sign(Delta) directly would not justify the correlation.

At least epsilon/(2m) of source mass has prediction advantage at least
epsilon/(2m) after averaging the seed (using advantage bounded by 1).
For each such x, choose values of seed bits outside S_i that retain this
advantage. Store those values using at most D bits. Each preceding code
bit, after the outside seed is fixed, depends only on S_j intersect S_i;
its entire truth table costs at most 2^a advice bits. All preceding bits
therefore require at most m*2^a bits. The predictor with this advice gives
a sign word on the ell-bit code index with correlation at least
epsilon/m, hence at least 2gamma.

For each fixed advice the code list bound leaves at most 1/(3gamma^2)
possible x. Thus the set of such x has cardinality at most
2^(D+m*2^a)/(3gamma^2), and its probability under any min-entropy-k source
is at most this cardinality times 2^-k. The displayed conservative k
inequality makes this less than epsilon/(2m), a contradiction. The stated
inequality includes more logarithmic slack than the count requires.

## Long seed, stretched exponentially small error

Take an integer a>=1. Let Q be a power of two with ell<=Q<2ell.
Use a fixed ell-point subset T of F_Q and let S_i be the graphs over T
of distinct degree-at-most-a polynomials. There are Q^(a+1) indexed graphs,
which are distinct when a<ell, as eventually holds for the fixed A=4096
used below. Their universe has D=ell*Q<2ell^2 elements, and intersections have size
at most a by the polynomial root bound. Enumerate their coefficients
lexicographically. As long as m<=Q^(a+1), the design is explicit.

For a desired input entropy k polynomially large in the ambient length n,
choose epsilon=2^(-ceil(k^e)) for a sufficiently small positive constant e.
Then w=O(log n+k^e), ell=2w and D=O((log n+k^e)^2)=o(k) if e<1/2.
Choose fixed a large enough that Q^(a+1)>=k, and output
m=floor(k/2^(a+4)). The advice term m*2^a<=k/16 and all other terms are
o(k). This gives a constant-fraction output, a seed of length O(k^(2e)),
and error 2^(-ceil(k^e)) at every sufficiently large length.

The constant fraction can be extremely small because a may be large.
Rao's condensation must therefore be retuned to this fraction; it cannot
retain the published nearly-full-length output parameters.

## Short seed, constant error, polynomial output

Fix a positive rational kappa and assume k>=n^kappa. Choose m=floor(k^(1/8)),
epsilon=1/4, and a=floor(log_2(k)/4). Here w=O(log n), ell=2w<=C log_2 n
for an absolute C, and m*2^a<=k^(3/8). Construct a design as follows.

Choose a sufficiently large constant power-of-two alphabet size Q,
depending on kappa. Use graph sets of Q-ary words of length ell, and greedily
select words with pairwise agreement at most a. For any previously selected
word, the fraction of all words agreeing in more than a coordinates is at
most 2^ell Q^(-(a+1)). For sufficiently large n, a>=(kappa/8)log_2 n. Taking
log_2 Q >= 8(C+3)/kappa makes this fraction at most n^-3. Since m<=n,
the forbidden fractions sum to less than one. Greedy selection succeeds.

Enumerating all Q^ell words takes n^O_kappa(1) time. The seed universe has
D=Q*ell=O_kappa(log n) coordinates. The reconstruction inequality holds
eventually because D+log(m)=O_kappa(log n)=o(k), and m*2^a=o(k).
This supplies precisely the O(log n)-seed, polynomial-output, linear,
constant-error ingredient needed for enumerating challenge responses.

The fixed exponent of the polynomial runtime is large; no practical speed
claim is made. Both parameters and runtime must be stated honestly in the
ultimate family, with total small-length guards as needed.

## Fields and uniformity

For degree O(log N), brute-force irreducible search is polynomial in N;
Frobenius/gcd testing is polynomial in the degree. In particular the
polynomial-graph design uses degree z=ceil(log2 ell), retaining Q<2ell.
Code-field degrees, which may be polynomial in N, are rounded upward to
powers of two, changing code-index length by at most a constant factor.
They admit a fully explicit
quadratic tower: start at F_2, and extend by t^2+t+c for an element c of
absolute trace one. This is irreducible since the image of z->z^2+z is
exactly the trace-zero hyperplane. In the extension, t*c has absolute
trace one because its relative trace is c. Iterating constructs the tower
in time polynomial in its final degree. For example three base-field
products ac, bd, (a+b)(c+d), plus multiplication of bd by the fixed
trace-one modulus coefficient, give a four-call recurrence and O(w^2)
base-bit cost up to polynomial overhead. A faster recurrence is unnecessary.

The cubic disperser needs odd degrees, supplied by the validated prefix.


## Affine somewhere-random and block extractors

All parameters below depend on a single ambient length N. All data vectors
passed to seeded extractors are padded to length N. Constants are fixed:

    e=1/1024, A=4096, c=2^-4106, gamma=2^-20,
    K=floor(N^(delta/2)), R=floor(K^(1/8)),
    G=2^floor(floor(log2 R)/(2*2^20)).

Only sufficiently large N matters for the extractor assertions; family.py
explicitly guards undefined or insufficient shapes.

## Affine independent residual

For an affine source X and a linear observation T(X), choose a complement
to the observation kernel in the underlying linear space. Then
X=L(T(X))+B, where L can include a fixed shift, B is uniform on a linear
space and independent of T(X), and H(B)=H(X)-H(T(X)). Likewise if
(X_1,X_2) is an affine block source with H(X_2|X_1)>=k, write
X_2=L(X_1)+B with B independent of X_1 and H(B)>=k. This follows by
extending a basis for the kernel to a basis of the source space.

A linear image of an affine source is uniform on an affine space. If its
image is a proper subspace of the full output space, its distance from
uniform is at least 1/2. Consequently a strong linear seeded extractor
with error epsilon has at most 2epsilon fraction of seeds whose image
on any fixed entropy-eligible affine source is not EXACTLY uniform.

## Retuned somewhere-random extractor

Input an affine G by R matrix with at least one exactly uniform row.
Let L=log2G, r_0=R and r_(i+1)=floor(r_i/2^4106). At step i the row
count is t_i=G/2^i. Use a long seeded linear extractor E_i on N data
bits, entropy floor(r_i/2), output r_(i+1) and error
epsilon_i=2^-ceil(floor(r_i/2)^e). The polynomial-graph design order is
A. Its seed length d_i is O(log(N)^2+r_i^(2e)); since r_i is a fixed
positive power of N, this is O(r_i^(2e)) eventually. The reconstruction
inequality holds: output*2^A <=r_i/1024 and all remaining terms are o(r_i).

Take w_i=floor(r_i^(1/10)) and m_i=2w_i-1. Concatenate the first w_i
bits of each adjacent pair of rows and drop the last bit. This is an
m_i-bit word, interpreted in a supplied degree-m_i field. Apply x->x^7
and take its first d_i coordinate bits as a seed. Use that seed on the
whole original t_i by r_i matrix through E_i. This produces the pair's
new row. Do it for every pair.

For large N, all r_i>=R^.99. Indeed L<=gamma log2R/2, and floor division
retains at least r_i/2^4107 while r_i is large, giving
r_L>=R^(1-4107gamma/2)>R^.99. Thus t_i<=r_i^gamma eventually, and
t_i*w_i=o(r_i). Also d_i=O(r_i^(2e))<m_i^(1/16) eventually.

Condition all t_i slices at once. The independent residual B has entropy
at least r_i-t_i*w_i>=r_i/2. The pair containing a uniform row has an
affine slice of entropy at least w_i-1=m_i/2-1/2. The cubic extractor
therefore makes its seed beta_i-close to uniform with
beta_i<=2^(-r_i^(1/100)) eventually. The half-rate-minus-one claim follows
from the same constant-space contradiction as half rate: the missing one
bit is negligible compared to the linear gap to .45m_i.

For this pair, the set of bad seeds for E_i on B has uniform probability
at most 2epsilon_i. Its actual seed hits this set with probability at most
2epsilon_i+beta_i. For every other slice fixing, all output rows are affine
linear functions of B, with the distinguished row exactly uniform. The
output is therefore within 2epsilon_i+beta_i of a mixture of affine
somewhere-random sources. A bad component may be replaced by any valid
affine somewhere-random source, which proves the distance statement.

Apply this argument to each mixture component and repeat for L rounds.
Deterministic processing cannot increase distance. The total error is at
most sum_i(2epsilon_i+beta_i)<=2^(-R^(e/2)) eventually. The last row has
length at least R^.99, so return its first floor(sqrt R) bits. This gives
a polynomial-time affine somewhere-random extractor for exactly G rows.
Fewer rows can be padded with zero rows before calling it.

## Extracting from a fixed number of affine blocks

Put alpha=1024*sigma, C=2^(alpha+4), and b=2^(alpha+log2(sigma)+40). A b-block affine source consists of vectors
X_1,...,X_b in F_2^N with H(X_j|X_1,...,X_(j-1))>=K for every j.
Use the short-seed linear strong extractor E' with input N, entropy K,
output R and error 1/4. Its design has intersection
floor(log2K/4), alphabet 2^alpha, and code index length at most 8log2N
eventually. Its seed length D' satisfies D'<=C log2N eventually.
The actual fixed alphabet satisfies the greedy existence inequality:
a=floor(log2K/4)>=delta*log2N/16 eventually, so alpha*a>=64log2N,
whereas ell_code+log2R<9log2N. Thus R*2^ell_code<2^(alpha*(a+1)).
The seed constant C bounds 8*2^alpha with slack.

Enumerate all seeds to form the v=2^D' by R matrix M of E'(X_1,y).
This matrix is affine, and at least one row is exactly uniform: otherwise
the strong average error would be at least 1/2 instead of at most 1/4.

Partition M into groups of G rows, pad the last group with zeros, and
apply the retuned somewhere-random extractor to each group. Its output
has floor(sqrt R) bits. Truncate each seed to the length needed by a long
linear seeded extractor E on entropy K, input N, output R and error
epsilon=2^-ceil(K^e), again with design order A. This seed length is
O(K^(2e))<sqrt R eventually. Apply E to the next fresh block with every
resulting seed to produce a new matrix with ceil(v/G) rows of length R.
Repeat until only one row remains, then return that row.

The total number of stages is at most ceil(D'/log2G). Since
log2G >=(gamma/4)log2R eventually and log2R >=(1/16)log2K,
log2K>=(delta/3)log2N eventually. More sharply, using its limiting
coefficient delta/2, the stage count is at most 128*sigma*C/gamma+O(1)
eventually with harmless slack in the preceding inequalities. A looser
256*sigma*C/gamma+O(1) bound follows directly from log2K>=delta*log2N/3
and is also sufficient. Since b/C=sigma*2^36, this number plus one is
below b for large N. Each needed input block
exists; unused final blocks are ignored.

For correctness condition the blocks already consumed, leaving the current
matrix affine with a uniform row, except on an error mass tracked below.
Only the group containing that uniform row matters. Its seed is within
2^(-R^(e/2)) of uniform. Write the next block as a function of all prior
blocks plus an independent affine residual B of entropy at least K.
Except with probability 2^(-R^(e/2))+2epsilon, the useful group's seed is
good for E on B. After fixing the current consumed block, the next matrix
is an affine image of the fresh block, and its useful row is exactly
uniform. Future blocks retain their conditional entropy K. This establishes
the induction on mixtures by coupling away the exceptional mass at each
stage. It does not assume that the unconditioned nonlinear-seed matrix is
affine, and it does not union-bound over all groups.

With at most b stages, the final output is within

    b*(2^(-R^(e/2))+2epsilon) <=2^(-K^((1/8)*e/4))

of uniform for sufficiently large N. This is the desired block extractor
BE. Its output can be truncated to
pmax=floor(N^(delta/2^24)), and its error is then at most 2^-pmax eventually:
the exponent of K^((1/8)*e/4) in N is delta/2^16, strictly larger.

## Field certificates and complexity accounting

All somewhere-random calls use the same R,G and row lengths r_i. The
odd field degrees m_i=2floor(r_i^(1/10))-1 are therefore a deterministic
list of O(log N) degrees. The largest is O(R^(1/10))=O(N^(delta/160)).
A prefix supplying one monic irreducible polynomial per distinct degree
has total length O(N^(delta/160)log N)=o(N^delta). Validate each exactly by
the general Frobenius/gcd test. No field polynomial is advice.

Seeded code fields use the explicit quadratic towers. Polynomial-graph
design fields have degree O(log N), so exhaustive irreducible search is
polynomial in N. All full-seed enumerations use at most N^C seeds;
all other loops have polynomially bounded lengths with fixed exponents.
The very large constants are part of one fixed uniform polynomial-time
algorithm, and no practical runtime claim is made.


## The complete challenge-response construction

Fix ambient length N and use the preceding block extractor with

    alpha=1024*sigma, b=2^(alpha+log2(sigma)+40),
    E_tree=alpha+3log2(sigma)+90, epsilon_tree=2^-E_tree,
    H=2^(E_tree+2), q=2^(alpha+2log2(sigma)+80),
    J=2^(14+log2(sigma)),
    k=ceil(N^delta), K=floor(N^(delta/2)), t=floor(N^epsilon_tree),
    ell=H*t^b+1, pmax=floor(N^(delta/2^24)).

Then epsilon_tree*b=delta^2/2^50, epsilon_tree*q=delta/1024,
and H/q=4096*sigma<J. These identities are exact fixed-constant identities.

The data is an N-bit vector. All parameter statements in this argument
are for sufficiently large integer N. Every displayed constant is fixed,
so factors depending only on these constants do not affect exponent gaps.

## Balanced tree and challenges

Split each interval into t consecutive intervals of lengths differing by
at most one. Stop at intervals shorter than floor(N^(delta/4)). For large N
every internal interval has length greater than 2t. Tree depth is below H:
at depth z interval lengths are at most N/t^z+2, and log t is eventually
at least (epsilon_tree/2)log N. The number of nodes is polynomial in N.
Levels divisible by q are active. Their number is at most 4096*sigma+1<J.

For an interval s let x_s mean the N-bit vector with coordinates outside
s zeroed. Let E' be the short-seed extractor on entropy K, padded input
N, output R from the preceding note. Write S=2^D' for its number of
seeds. We have S<=N^C eventually for the fixed C=2^(alpha+4). Only the first
ell*pmax output bits are used, as an ell by pmax response matrix; this
truncation is allowed since ell*pmax<R eventually. It stays linear and
strong with error at most 1/4.

For any challenge ch on a subinterval s' of s, and confidence p<=pmax,
CR_p(x_s,s') is FAIL if the first p bits of every challenge row agree with
the first p bits of every corresponding response row for at least one
seed of E'. Missing challenge rows are zero padded to ell. Otherwise it
is PASS. Here each response is E'(x_s,seed) truncated as just described.

Compute challenge matrices bottom up. A leaf has no rows. At an internal
node s with children s_1,...,s_t, first form one BE(x_(s_i1),...,x_(s_ib))
row for every ordered-index b-subset i1<...<ib. Then find the first child
for which CR_pmax(x_s,s_i) is PASS, choosing the first child if none
passes, and append that child's challenge rows AFTER the node's own BE
rows. This ordering is deliberate: the node's own rows have fixed indices
even when child subtrees have different final depths. At most H*t^b<ell
rows result. Computation is polynomial in N because b,H,C are fixed.

The same first-PASS rule gives a path from the root to a leaf. For every
pair of consecutive active nodes on this path, at levels uq and (u+1)q,
compute CR at confidences

    p_u^(a)=ell^(6*(J-u)+3a), a=1,2.

Choose the first such active node whose confidence-1 result is FAIL and
output its confidence-2 result, interpreting FAIL as 0 and PASS as 1.
If there is no such node output 0. The maximal confidence times ell is
at most ell^(6J+7), which is o(pmax). Indeed ell=N^(delta^2/2^50+o(1)), and
(6J+7)*delta^2/2^50 <delta/2^24. Thus all calls are well defined eventually.

## Two elementary CR facts

If an affine source X on interval s has conditional entropy
H(X_s|X_(s'))>=K, then CR_p has FAIL probability at least 2^(-ell*p).
For each value of X_(s'), its challenge is fixed and the conditional
remaining source is affine of entropy at least K. Some seed of the linear
strong extractor is exactly uniform on its ell*pmax output bits, since
otherwise its strong error would be at least 1/2. Equating the ell*p
observed response bits with the fixed challenge cuts out an affine
subsource of codimension ell*p. Average over X_(s').

In particular, if X_(s') is constant and H(X_s)>=K, there is an affine
subsource of deficiency ell*p on which CR_p always FAILS. This assertion
requires no extraction property of the challenge function.

Say a challenge matrix is eta-somewhere-uniform with resilience r if,
after every possible linear observation of rank at most r, at least one
fixed row is eta-close to uniform. Which row may depend on the observation
and its value. If this holds with r>=ell*p, then CR_p has FAIL probability
at most S*(2^-p+eta). Fix each candidate response slice, a linear
observation of rank at most ell*p; the challenge's near-uniform row matches
that fixed row with probability at most 2^-p+eta. Union-bound over seeds.

These assertions apply with respect to a larger ambient affine source as
well: observing a linear function of one interval is an observation of the
whole vector. FAIL at a larger confidence implies FAIL at any smaller
confidence by truncation.

## Finding a useful configuration of intervals

Take any affine source X of entropy at least k. We construct an affine
subsource and intervals e subset g' subset g, where g,g' are consecutive
active nodes, with three properties:

1. H(X_g|X_(g'))>=k^(3/4).
2. The children of e contain b ordered blocks forming an affine block-wise
   source with threshold at least k^(3/4).
3. Every coordinate strictly before e is fixed.

For the first step suppose a node s has entropy at least a and all earlier
coordinates fixed. Divide it into T=t^q next-active intervals; these full
descendants exist during this search because the current entropy remains
a fixed constant times N^delta, far above the leaf cutoff even after one
q-level jump. Let a_i be the conditional entropy contributed by interval i
after all earlier ones. Choose the first i with a_i>=a/(4T), and fix the
earlier intervals. This costs at most a/4 entropy. The chosen interval has
entropy at least a/(4T). If its entropy is at least a/4, continue within
it with a replaced by a/4. Otherwise the entropy remaining outside the
chosen interval is at least a/2, giving a strongly structured pair.

The heavy-interval case cannot continue through all active levels: entropy
would remain at least k/4^J, whereas the eventual intervals are shorter
than the leaf cutoff N^(delta/4). On success let the pair be g,g'. Both
H(X_(g')) and H(X_g|X_(g')) are at least
k/(4^(J+1)*t^q)>=k^(4/5) eventually, since t^q<=N^(delta/1024).

Starting at g', search down ordinary levels. At a node of entropy a,
if b child conditional entropies are at least a/(4t), those children
are the desired hidden block source. Otherwise fewer than b children
contribute that much, and the total contribution of all smaller children
is at most a/4. Some child therefore contributes at least a/(4b).
Fix its earlier coordinates and continue in it. These restrictions lie
inside g', so H(X_g|X_(g')) is unchanged: conditional affine entropy is
the dimension of a fixed projection kernel. After at most H such moves,
the preserved entropy is at least k^(4/5)/(4b)^H, which still greatly
exceeds the leaf size. The search must therefore find the hidden block
source at some internal e. Its threshold is at least
k^(4/5)/((4b)^H*4t)>=k^(3/4) eventually. All coordinates before e are fixed.
For the selected b children, conditioning only on earlier SELECTED children
costs no more entropy than conditioning on the entire earlier prefix, so
they have the required block-wise entropy threshold.

For every edge on the root-to-e path, force the CR_pmax tests on all
children left of that path to FAIL, using the empty-interval fact. There
are at most H*t such tests, each costing ell*pmax dimensions. Total loss
is at most H*t*ell*pmax=o(k^(3/4)). At every step the parent still has
entropy at least K, so the refinement is legal. The final affine subsource
X* retains conditional entropy H(X*_g|X*_(g'))>=k^(2/3), hidden-block
threshold at e at least k^(2/3), and all those left-child failures hold
identically. Denote this fixed final subsource by X* for the rest of the
argument. No probabilistic condition was silently substituted for an
affine restriction.

## Resilient challenge propagation

At e, after any linear observation of rank at most (H+1)*ell*pmax, the
chosen b children retain conditional entropy at least
k^(2/3)-(H+1)*ell*pmax>=K. Their BE row is therefore 2^-pmax-close to
uniform. Its position among the node's own rows is fixed, so the whole
challenge has that somewhere-uniform property and resilience.

Move upward a edges from e. Inductively spend ell*pmax resilience per
edge. For any observation of the remaining allowed rank, the true child's
challenge remains resilient to the further response slice. Its test passes
except with probability S*(2^-pmax+eta). Every left child fails identically,
even after the observation. The parent therefore selects the true child
except with that probability. Its useful row occupies a fixed offset after
the parent's own fixed list of rows. The new row's error is at most
(S+1)*eta+S*2^-pmax. Thus after a edges error is bounded by
(2S+2)^(a+1)*2^-pmax and resilience remains at least ell*pmax.

Since a<=H and log S=O(log N) with fixed constants, all these errors are
at most 2^(-pmax/2) eventually. Applying the CR pass estimate on every
edge to g' and union-bounding shows that the selected path visits g'
except with probability at most H*2^(-pmax/3). The fixed left failures
prevent any earlier choice. They persist through all source observations.

For any confidence p among the active-node confidences, the same bound
gives FAIL probability at most 2^(-p/2) eventually: log S is negligible
compared to the smallest such p, and pmax is asymptotically much larger
than the largest confidence. These statements also hold for active jumps,
which simply test a deeper descendant instead of one child.

## Both output values occur

Let g have active index u. Write r_a for the fixed test
CR_(p_u^(a))(X*_g,g'), whether or not the algorithm's path reaches g'.
The conditional entropy H(X*_g|X*_(g'))>=K gives
Pr[r_a=FAIL]>=2^(-ell*p_u^(a)). The upper pass estimate gives
Pr[r_2=FAIL]<=2^(-p_u^(2)/2).

For output 0, r_2=FAIL implies r_1=FAIL by monotonicity, with probability
at least 2^(-ell*p_u^(2)). For output 1, subtraction gives

    Pr[r_1=FAIL, r_2=PASS]
    >=2^(-ell*p_u^(1))-2^(-p_u^(2)/2)
    >=2^(-ell*p_u^(2))

eventually, since p_u^(2)=ell^3*p_u^(1) and ell tends to infinity.

An earlier active node at index u'<u fails its confidence-1 test with
probability at most 2^(-p_(u')^(1)/2). The nearest possible earlier one
has p_(u-1)^(1)=ell^3*p_u^(2); all others have even larger confidence.
Their total failure probability is at most
J*2^(-ell^3*p_u^(2)/2), negligible compared to
2^(-ell*p_u^(2)). The path-misses-g' probability H*2^(-pmax/3) is also
negligible compared to it, because pmax/(ell*p_u^(2)) tends to infinity.

Subtract these two exceptional probabilities from either positive r-event.
A positive probability remains that the path visits g', that no earlier
active node overrides g, and that g selects and emits the requested value.
Thus both output bits occur on X*. They occur on X as well. The construction
is nonconstant on every affine source of entropy at least ceil(N^delta).

## Runtime and finite guards

All full-seed loops have S=N^O(1) iterations. The b-subset list has at most
t^b=N^O(1) entries; the tree has polynomially many nodes; BE and every
field/code operation have polynomial cost. All exponents and enormous
constants are fixed. The evaluator returns zero
when e.g. t<b, the tree's cutoff is too small, a response is longer than
R, a requested seed is too long, or any listed shape is undefined. Those
guards are eventually false by the strict exponent inequalities above.
The actual-length prefix charge and final gate theorem are assembled in proof.md.
