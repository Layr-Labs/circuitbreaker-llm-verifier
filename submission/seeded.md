# Linear extraction with a bounded-order diagonal orbit

This manuscript supplies all field, condenser and extractor premises.
The folded rank argument is the folded-Wronskian idea (Forbes--Guruswami,
ECCC TR14-162, Construction 2.7). The reconstruction organization adapts
Shaltiel--Umans, JACM52(2),2005,Sections4.2--4.5, with an important change:
a diagonal orbit supplies a Vandermonde information set, so neither a
primitive element of the full multiplicative group nor a transitive
orbit through the entire field extension is needed. The arguments and
changed algorithms are proved here; no optimized extractor theorem or
efficient list decoder is an external premise.

All entropies and dimensions are binary. A strong extractor of error
eta means that its output together with an independent uniform seed is
within total variation eta of independent uniform output and that same
seed. Every fixed-seed map below is F_2-linear in its data. The internal
long calls use precisely the same elementary construction as Short.

## 1. Explicit fields of large degree, without exhaustive field search

For every power-of-two degree w, build F_(2^w) by quadratic towers.
Start with F_2 and absolute-trace-one element c=1. Given field K and
trace-one c, adjoin t with t^2+t+c=0. The map z->z^2+z has kernel
{0,1} and image equal to the trace-zero subspace: it is contained there,
has |K|/2 elements, and trace is a nonzero binary linear functional
because c has trace one. Thus the polynomial is irreducible.
In the new field, t*c has relative trace c and absolute trace one,
providing the next tower's distinguished element.

Multiplication is

    (a+bt)(d+et)=(ad+be*c)+(ae+bd+be)t.

The base-field products ad,be,(a+b)(d+e),be*c give all terms. Thus
T(2w)<=4T(w)+poly(w) and bit arithmetic is polynomial in w. The code
builds the tower and performs arithmetic using explicit stacks; it
relies on no unbounded recursion limit. Its basis and all choices are
fixed, uniform and deterministic. These fields may be far larger than
any fixed polynomial in the original input length; their representation
and arithmetic costs are polynomial in w, not in 2^w.

Small fields in the finite audit or graph design may use exhaustive
irreducible search instead. Its exact tests are X^(2^d)=X modulo P
and gcd(P,X^(2^(d/p))-X)=1 for each prime p|d. Frobenius orbits and the
factorization of X^(2^d)-X prove these tests, and
2^d=sum_(t|d)t I_t together with the proper-divisor bound
sum_(t<=d/2)2^t<2^d proves existence of a degree-d irreducible.
The graph field degree is O(log M)=O(log N), with a fixed constant,
so its exhaustive search is polynomial in actual N.

## 2. A sufficiently large order is easy to find here

For an integer U>=1 and q>U^2+2, a field F_q contains an explicitly
findable element gamma of multiplicative order greater than U.
Indeed at most sum_(j=1)^U j=U(U+1)/2 nonzero elements have order
at most U: each is a root of some X^j-1, which has at most j roots.
Try the first U^2+1 nonzero elements in the fixed binary ordering,
and for each test gamma,gamma^2,...,gamma^U for equality with one.
At least one succeeds. The algorithm uses O(U^3) field operations,
without factoring q-1 or enumerating F_q. The slightly stronger coded
field-size guard guarantees that this many distinct nonzero candidates
fit in the field. It will always hold at the parameters used below.

Fix ambient data length A, entropy threshold H<=A, requested output m,
and error parameter E>=1. Let

    T=E+bitlength(A)+bitlength(m+2)+10,
    w=the least power of two >=128T, q=2^w,
    d=128, h=floor(H^(1/64)), U=max(A,h^128).

Then h^128<=H^2<=A^2, so U<=A^2. Since w>=128T, q exceeds U^2+2
with ample slack for A>=2; the finite exceptional cases are guarded.
Finding gamma of order greater than U takes a fixed polynomial in A
and w. Even when E is polylogarithmic and q is superpolynomial in A,
there is no hidden exhaustive search of q elements.

## 3. Folded evaluation preserves affine entropy

Pack the A input bits injectively into coefficients of a polynomial f
over F_q, of degree below a=ceil(A/w). Pad its last coefficient by zeros.
For uniform seed b in F_q, output H field evaluations

    Cond_b(f)=(f(b),f(gamma*b),...,f(gamma^(H-1)*b)).

Each fixed-seed map is binary linear, with H*w output bits. On any
affine source of binary dimension at least H choose an H-dimensional
direction subspace V and take its span W over F_q. Its dimension r<=H.
Row reduction gives a basis with distinct leading degrees t_1,...,t_r<a.
The determinant of their first r folded evaluations is a nonzero
polynomial in b: its leading coefficient is the product of the leading
coefficients times the Vandermonde determinant on gamma^(t_i).
These elements are distinct because gamma has order greater than A.
Its degree is at most H*A. For all but H*A/q seeds, the folded map is
injective on W, hence on V, so its affine output has entropy at least H.
Affine offsets change no rank. It is not claimed already fully uniform.

The bad fraction is at most A^2/q<=2^-E/2 at all sufficiently large
parameters, directly from w>=128T. Evaluation and field-order search
are polynomial in A,w. The condenser supplies no advice and adds no
input bits to the final target function.

## 4. A diagonal information set for the compact source

Assume the exact capacity guard

    binom(h+127,128)>=H.

Encode the H*w compact input bits injectively as the first H
coefficients of degree-lexicographically ordered monomials of total
degree less than h in d=128 variables over F_q. Call the polynomial
P_x. Each coordinate exponent lies between 0 and h-1. Set

    D0=diag(gamma,gamma^h,...,gamma^(h^127)).

This is an invertible linear map. For uniform seed v in F_q^128 and
uniform binary vector j of length w, output m bits

    Ext(x;v,j)_i=<P_x(D0^i v),j>, 1<=i<=m.

The brackets are binary inner product in the fixed tower basis. The
map is linear in x for each seed. Its computation uses polynomially
many field operations in A,w,m; every h power has bit length polynomial
in log A, and the dimension is the fixed integer 128.

The key reconstruction property is elementary. Write each monomial
alpha=(alpha_0,...,alpha_127) as the base-h integer
c(alpha)=sum_j alpha_j h^j. Distinct monomials have distinct codes
between 0 and h^128-1. The roots gamma^(c(alpha)) are therefore distinct.
For a vector v with every coordinate nonzero,

    P_x(D0^t v)=sum_alpha a_alpha*v^alpha*(gamma^(c(alpha)))^t.

The H-by-H matrix at t=0,...,H-1 is a Vandermonde matrix multiplied by
nonzero column scales v^alpha. It is invertible. Thus these H orbit
values determine all H coefficients, and hence x. A uniform v has all
coordinates nonzero with probability at least 1-128/q. Transitivity on
F_q^128 is unnecessary. The evaluator never computes this interpolation;
it is an information-theoretic reconstruction argument.

## 5. Full strong-extraction proof

Put eta=2^-E and rho=eta/(4m). We show the algebraic extractor has
strong error at most eta/2 for any compact source of min-entropy H,
under the exact reconstruction-advice guard below. Its proof also works
when a predictor sees the seed, which is essential for strong extraction.

Two elementary estimates will be used. First, if distinct univariate
polynomials of degree <=D each agree with >=B pairs in a set of P
pairs, their number is at most 2P/B when B^2>=2DP. Choose B agreement
pairs each; two selected sets intersect in <=D pairs. Cauchy--Schwarz
on incidence counts gives
(LB)^2<=P*(LB+L(L-1)D), yielding the assertion. Enumeration is allowed;
no efficient list-decoding algorithm is claimed.

Second, for 2r-wise independent Bernoulli variables of mean p, with
sum Z and mean mu>=1, expanding E(Z-mu)^(2r) eliminates terms with
singleton indices and leaves at most r distinct indices. Each repeated
centered Bernoulli moment has absolute value at most 2p. Hence

    E(Z-mu)^(2r)<=r*(2r)^(2r)*(2mu)^r.

Markov bounds Pr[Z<mu/2] by K_r*mu^-r, with
K_r<=r*(2r)^(2r)*2^(3r). Fix r=128.

### 5.1 Reconstruct a fixed predictable input

Fix an output index i* and a q-ary list predictor which, given seed v
and the preceding i*-1 field values, produces at most rho^-2 candidates
for P_x(D0^i* v). Fix an x correct on at least a rho/2 fraction of v.
Choose 2r distinct uniform field positions t_1,...,t_(2r), and independent
uniform vectors y_1,...,y_(2r) in F_q^128. Let p1,p2 be the unique vector
curves of degree <=2r-1 with

    p1(t_i)=y_i for all i,
    p2(t_i)=D0*y_i for i<=r, and p2(t_i)=y_i for i>r.

Define curves P_(2s+1)=D0^s p1, P_(2s+2)=D0^s p2. Their values at
any 2r distinct positions are independent uniform vectors. Adjacent
curves agree on r of the random anchor positions. The relevant anchor
set is independent of the CURRENT curve: conditional on all positions,
that curve's 2r values are independent uniform because D0 is invertible.
Thus its polynomial coefficient vector is uniform independently of
those positions. This statement remains true for every adjacent pair.

Advice supplies the composition P_x with the first 2i* curves,
INCLUDING two curves when i*=1. Its bit length is at most

    a_adv=2m*(255*(h-1)+1)*w.

This provides an anchor curve even at the first reconstructed output.
Each composition has univariate degree <=D=255*(h-1). For a new curve
point u apply the predictor at D0^(-i*)u. Its previous field inputs are
P_x(D0^(-i*+1)u),...,P_x(D0^-1 u), available from earlier curves since
D0^-j P_t=P_(t-2j). At the first new curve, index 2i*+1, the smallest
required predecessor is three when i*>1; when i*=1 there are no prefix
inputs and the supplied second curve provides anchors.

A current curve has, except with probability K_r*(rho*q/2)^-r, at
least rho*q/4 correct candidate positions. There are at most q*rho^-2
predicted pairs. The list estimate leaves at most 16rho^-3 possible
univariate compositions, since q>>D*rho^-4. A wrong composition agrees
with the true one at at most D positions. The r shared anchors select
the right one except with probability at most

    16rho^-3*(2D/q)^r.

The list depends only on the current curve, fixed predictor and x, not
on the independent-of-current anchor set. Bound these events with TRUE
predecessor values and union bound; no conditioning on earlier successful
reconstructions is needed. On their common complement the actual process
agrees with the true predecessors inductively.

Here m<=A and H<=A whenever the advice guard holds, h<=A^(1/64),
and rho>=2^(-T-2). Since q>=2^(128T), both error bounds are eventually
smaller than 1/(64A^2). For example, their binary logarithms are at most
a fixed constant minus 16000T, whereas log_2(A^2)<=2T. These estimates
include all rounding constants and fixed r, and q>>D*rho^-4 as well.
Reconstructing through at most 2H+2m<=4A^2 curves therefore succeeds
with probability at least 15/16 eventually. With additional probability
loss at most 128/q, y_1 has no zero coordinate. On that event the odd
curves at t_1 reveal P_x(D0^s y_1) for s=0,...,H-1, which by Section 4
determine x. Overall reconstruction succeeds with probability >1/2.

### 5.2 Entropy count and binary conversion

Fix the predictor. Averaging reconstruction randomness over predictable
x shows that their total mass under a min-entropy-H source is at most
2^(a_adv-H+1): for fixed randomness and advice a deterministic decoder
returns at most one x, and there are at most 2^a_adv advice strings.
All other x have prediction success below rho/2. Thus total q-ary
prediction success is at most

    2^(a_adv-H+1)+rho/2<=rho

under the exact coded guard

    H>=a_adv+E+ceil(log_2(16m)).

Distinct field values have orthogonal Hadamard sign words
(<z,j>)_(j in F_2^w). Bessel's inequality implies that at most rho^-2
such words agree with any proposed prediction word on >=1/2+rho of
its positions. If strong binary distance exceeded 2rho*m, the output
hybrid would yield an index and binary predictor with average advantage
>2rho, seeing (v,j) and the previous output bits. On more than a rho
fraction of (x,v), its word over j then agrees with the true Hadamard
word on >1/2+rho positions: otherwise average advantage is at most
rho/2+rho<2rho. Given the preceding FIELD values it can form that word
for all j and list its <=rho^-2 candidates. This is a q-ary predictor
seeing v with success >rho, a contradiction. The strong binary error
is therefore at most 2rho*m=eta/2.

## 6. The complete linear extractor and its long-call budget

Compose the folded condenser (independent seed b) with the algebraic
extractor (independent seeds v,j), using the same field and gamma.
Bad condenser seeds have fraction <=eta/2; all others give compact
affine sources with entropy at least H. The Section 5 strong error
on each gives total <=eta. The combined seed has exactly 130w bits.
All fixed-seed maps are binary linear, so on an affine source a rank-
deficient output has distance >=1/2 from uniform. At most 2eta of
uniform seeds are bad; every other image is EXACTLY uniform. Their
good set is fixed by the direction space, not the sample or offset.

For a long call of ambient length A, output m and error 2^-E, use the
same w as above and set the sufficient entropy threshold to

    H_call=ceil((2^20*m*w)^(33/32)).

The evaluator computes the rational power exactly by an integer root
of an integer power. It guards H_call<=A, coefficient capacity, the
exact advice inequality and supplied seed length. As m*w grows,
H_call^(63/64)/(m*w) tends to infinity, because
(33/32)*(63/64)=2079/2048>1. Thus its actual advice
O(m*H_call^(1/64)*w) is o(H_call), and the guard passes eventually.
The polynomial encoding has capacity Theta(H_call^2), exceeding its
H_call field coefficients. Exact guards cover exceptional shapes.

A larger supplied seed is handled by ignoring its suffix; the code
masks the inner-product seed as well as the field/vector seeds. This
preserves strong extraction and fixed-seed linearity. The total seed
is O(E+log A+log m), rather than quadratic in that quantity. No field-
size enumeration or transitive-orbit search is hidden in this bound.
The main proof pays the resulting entropy thresholds separately; it
does not assert constant-fraction output at entropy O(m).

For the initial Short call use H=k and E=16ceil(log_2(A+2)), with its
separately proved exact capacity/advice guards. Its error is <=A^-16,
hence <=2A^-8. Its seed length is <=2^20 log_2 A eventually. All calls
have fixed polynomial cost in their original ambient length and E;
the main family has E polynomially bounded in actual N. The finite
audit checks bounded order searches and small diagonal Vandermonde
matrices, not any universal extractor statement.
