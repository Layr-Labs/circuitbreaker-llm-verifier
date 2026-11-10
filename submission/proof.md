# A uniform certificate-prefix family: coefficient 4, exponent 1

## Claim and provenance

For the single total Boolean language defined below, with **N denoting the
entire input length**, we prove

    C_B2(F_N) >= 4N - R(N),       R(N)/N -> 0.

The joint pair is (4,1). The leading coefficient is inherited mathematics.
Section6 of elimination.md supplies an additive-seven theorem using one finite Ramsey tool,
without increasing that leading coefficient. The leading gate argument is the
parameter method of Mishra and Santhanam, ECCC TR26-260, revision 1, November 3, 2026,
Sections 2–4. The algebraic construction is based on Ben-Sasson
and Kopparty (BK), *Affine Dispersers from Subspace Polynomials*, author
version March 12, 2010, Sections 3, 6 and 7. The necessary arguments are
supplied in our notation in `elimination.md` and `disperser.md`.

Our input-prefix wrapper removes the deterministic irreducible-polynomial
construction prerequisite of the repository's BK/Shoup seed. It explicitly
pays for all prefix bits. The **whole** function is not an affine disperser:
only a fixed-prefix restriction is asserted to be one. No conjecture,
restricted circuit basis, direct-sum assumption or exponent gain is used.
We assert no exponent above 1, endpoint, or optimality over all families.

## 1. One evaluator on every bit string

The first bit of each field block is its constant coefficient. On x of length N:

1. If N=0, output 0.
2. Let B be the smallest positive integer with
   B^5>((27+18sqrt(2))/64)*N^3. Equivalently require
   g=64B^5-27N^3>0 and g^2>648N^6, using exact integers. Let k be the
   smallest prime at least B. Use deterministic trial division to find k.
3. If k>=N, output 0. Otherwise put m=N-k, r=ceil(m/k), t=ceil(sqrt(r)).
   If t>=r or k<=r+1, output 0.
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

Put c=(27+18sqrt(2))/64 and b=c^(1/5)>0. For B>=2 (in particular
at every sufficiently large length), Bertrand's prime theorem gives
k<2B. For B=1, k=2 still satisfies the displayed bound. Thus

    b*N^(3/5) < k < 2*b*N^(3/5)+2.

The prime scan and trial division are polynomial in N. B is computed by
integer binary search and the exact positive-gap square comparison.
The target real quantity c*N^3 is strictly increasing, so B(N)
and k(N) are nondecreasing, although the implementation uses no real arithmetic. All subsequent integer,
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
Consequently k<N and t<r eventually. Also k/r=Theta(N^(1/5)) tends
to infinity, so k>r+1 at every sufficiently large length. Similarly
k/t^2 tends to infinity and supplies k>=t(t+1)/2 eventually. No product
inequality k*t>=r(r+1) is assumed: Section15 of disperser.md verifies
that the complete algebraic proof needs only the weaker k>r+1.
Every one of these conditions holds uniformly over all sufficiently
large integer N, rather than only a selected subsequence.

## 5. Hard restriction and its precise loss

For each sufficiently large N fix any irreducible Q of degree k and fix
the first k bits of x to its coefficients. The remaining function G_Q
has m inputs. Put

    W = r(r+3)/2, H = (t-1)(t+2)/2,
    A0 = (t-1)k-H-(r-t+1),
    L = floor(t*k/(t+1)-t/2), q=floor(L/t), a=L-t*q,
    B_int=r*q+a, B_joint=r*(2k-r-1)/(2*(t+1)),
    d_N = 1+W+max(A0,min(B_int,floor(B_joint))).

`disperser.md` proves that, for prime field degree k, r>t>=1, k>r+1, and
k>=t(t+1)/2, the function

    (a_i) -> pi(sum_(|I|=t) product_(i in I) a_i^(2^(r+2-i)-1))

is nonconstant on every affine subspace of K^r of dimension at least
d_N, by its Sections13–14 refinements, for any nonzero
F_2-linear pi. Our pi is the constant
coefficient, with pi(1)=1. The payload embedding is injective and linear,
so G_Q is an affine disperser of dimension d_N on its m actual input bits.
Zero padding preserves dimension and does not add circuit inputs.

Since kr=O(N), r=Theta(N^(2/5)) and t=Theta(N^(1/5)),

    d_N=O(N^(4/5)),          k=O(N^(3/5)).

We use this conservative sublinear disperser dimension, not later
logarithmic-dimension constructions. Either suffices for the leading pair.

Fixing an accepted prefix yields G_Q, nonconstant for every sufficiently
large N because d_N=o(N) and d_N<=m. Section6 of `elimination.md` proves,
for d>=3, in the strengthened model with free signed wires and a free
output literal,

    C(G_Q)>=4m-8d_N+7.

All16 binary truth tables, constants under substitution, repeated pins,
unary effects, arbitrary sharing and terminal-output cases are covered.
Here d_N>=3 eventually. The prefix-removal lemma below gives the stronger transfer

    C_B2(F_N)>=C_free(G_Q)+k
             >=4(N-k)-8d_N+7+k=4N-3k-8d_N+7.

Its remainder R(N)=3k+8d_N-7 is O(N^(4/5))=o(N). Extend R on the
finitely many guarded lengths. The exact finite tools support only their scoped predicates, with full
universal normalization and gate transfer supplied in elimination.md.

### Prefix-removal lemma, including terminal signs

Let C be any fan-in2 circuit depending essentially on k specified inputs,
and fix those inputs to constants so that its remaining function is
nonconstant. In the strengthened model allowing free wire signs, constants
under restriction and an output literal, at least k distinct charged gates
can be removed. Consequently |C|>=C_free(C restricted)+k. Applying this
to a normalized no-larger strengthened circuit also proves a lower bound
for the original full-B2 model.

Trim to ancestors of the output. The k essential inputs are present. In
the ORIGINAL graph, retain each fixed input as its own vertex even when
two assigned values coincide. Let S contain these roots, original constant
nodes, and gates evaluated as constants by ordinary constant propagation.
Let T be the number of charged gates in S. Consider the undirected
components of the induced S graph that contain a fixed input. Initially
there are k such input components. Adding a constant gate merges at most
two existing components, reducing their number by at most one. Original
constant nodes alone merge none. Hence there are at least k-T components
containing fixed inputs (with a vacuous bound if T>=k).

Each such component has a path in the original trimmed graph to the
output. The output is not in S because the restricted function is
nonconstant. Take the first gate B outside S along such a path. It has a
constant parent and is not constant, so its restricted binary table is a
literal of its other parent. It can therefore be deleted by literal
rerouting. It has at most one parent in S: two constant parents would
make B constant. Thus one B cannot be the chosen exit of two distinct S
components. There are at least k-T distinct such boundary gates. All T
constant gates and these boundary gates can be removed simultaneously,
with constants replaced by zero-cost substitution nodes and literal signs
absorbed into descendants; the terminal output literal is explicitly free.
This removes at least T+(k-T)=k gates. Later simplification can only remove
more. This argument counts original gates once and permits fan-out,
coincident pins, original constant nodes and arbitrary shared descendants.

For this family, EVERY prefix bit is essential at all sufficiently large N.
Fix any accepted irreducible Q of prime degree k>1 and a payload with
G_Q=1, which exists by nonconstancy. Since Q has no root1, Q(1)=1.
Flipping any one lower coefficient gives Q'(1)=0, so Q' has factor Z+1
and fails the exact gcd prefix test. The same payload then outputs0.
This proves essentiality of each of the k specified bits. Fixing the
original accepted Q leaves nonconstant G_Q, so every premise of the lemma
holds. The proof does not assume that all payload bits are essential.

For each epsilon>0, take N_epsilon large enough that the guards pass and
R(N)/N<=epsilon for every N>=N_epsilon. This gives exactly
C_B2(F_N)>=(4-epsilon)N for all those integers. Both coordinates in
bound.json belong to this one theorem and one language.

## 6. Evidence and limits

`disperser.md` contains the coefficient separation and block-mixing proof.
`elimination.md` contains the complete additive-seven gate argument.
The declared audit.py experiment regenerates and refutes exactly four
normalized Ramsey CNFs and all4082400 normalized read-once formula
assignments required by its two finite tools; universal
normalization and transfer are proved in the same document. No claimed
review, checker success or prior acceptance is a mathematical premise.

Sources, for attribution and version identification:

* Mishra--Santhanam, ECCC TR26-260 revision 1, November 3, 2026, Sections 2–4.
  https://eccc.weizmann.ac.il/report/2026/260/ . The gate argument needed here is reproduced in elimination.md.
* Ben-Sasson--Kopparty, author preprint 2010-03-12, Lemmas 3.7, 3.10, 3.13,
  Theorem 6.1, Lemma 7.2, Theorem 7.3 and Claim 7.4,
  https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
  SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf.
* Li--Yang, ECCC TR21-023, 2021-02-20, Sections 1.1, 2.1, 2.4 and 3.3,
  for the protected full-B2 conventions. Their LP and lower-bound case
  analysis are not invoked by this submission.


## Exact finite scope and retained source

The evaluator framework, original algebraic proof and analytical gate foundation
come from635e6e27-f9ec-4693-8ccc-c99bd5bdf90f, commit
66e35592d7e5f899fdee0c4f8786fb4ba8ed6c4d. Sections1–4 of elimination.md
retain that source's complete additive-five argument. Section5 reproduces
the geometry and single unrestricted Ramsey R tool from94a62e03, commit
f972a0f5ca675db0d9978d82741c505983f24497, including its complete finite
normalization and exact-resolution soundness argument. audit.py extracts
its R84 generators, exact resolution kernel and read-once Q enumeration, with no family import.
All four actual normalized CNFs must derive the empty clause. A cap,
model, error or incomplete run fails. No prior verdict or stored opaque
SAT assertion is a premise.

Section6 supplies the compact additive-seven transfer. Its only finite
tools are the unrestricted Ramsey R and read-once formula Q lemmas.
The proof gives complete terminal pair-count and last-event cases,
including the conversion from exact acyclic object/pin counts to read-once
formulas. The RREF generator and full assignment-mask coverage are
justified. No six-variable Ramsey theorem,83-case computation,
near-read-once audit or additive-eight proof is needed or retained.

Together with prefix removal at the current field threshold, this gives
4N-3k-8d_N+7, two gates above the corresponding additive-five
certificate4N-3k-8d_N+5 at identical new-family parameters.
The native leading expression and projection remain(4,1). No leading
coefficient record or literature priority is claimed.

## Weaker algebraic guard and the current fixed family

The CURRENT total language uses c=(27+18sqrt(2))/64 and the monotone
length rule B^5>c*N^3. Its evaluator decides this comparison exactly as
specified in Section1, and explicitly rejects k<=r+1. Section15 of
disperser.md verifies the precise hypothesis strengthening: the old
product guard was used only to deduce k>r+1, and the entire supplied
coefficient and mixing argument applies with this direct condition.
The additional k>=t(t+1)/2 holds eventually by Section4.

This is ONE new fixed language throughout the joint declaration. The
earlier submitted32 rule and the separately prepared11/4-margin rule
define different languages; no proof swaps among them by coefficient
or length. Prime-prefix verification, field representation and payload
polynomial are unchanged. All inherited field, algebraic and gate
arguments are supplied and their hypotheses checked for current parameters.
Old32 numerical examples in disperser.md are labeled as supporting tuples.
They are not outputs of the current parameter routine.

Joint selection gives the current d_N from Section5 and eventually

 d_N=(t-1)k+(r(r+1)-t(t-1))/2+1.

The monotone k(N) and unbounded r(N)=ceil(N/k(N))-1 still give upward
steps at most one, so every sufficiently large square r is attained.
All supporting strict square-block threshold comparisons therefore also
hold at infinitely many actual current-family lengths. The MAIN circuit
theorem applies to every sufficiently large integer length.

## Balanced endpoint field sizing and a smaller uniform envelope

Let a_N=k/N^(3/5), so a_N lies in [b,2b+o(1)] by Section3. With
r=N/k+O(1), t=sqrt(r)+O(1), the joint threshold uniformly gives

 d_N=k*sqrt(r)+r^2/2+O(k+r),
 d_N/N^(4/5)=f(a_N)+o(1), f(a)=sqrt(a)+1/(2*a^2).

The derivative of f has the sign of a^(5/2)-2, so it has one minimum
and its maximum on any compact interval is at an endpoint. Our explicit
constant satisfies

 b^(5/2)=3*(sqrt(2)+1)/8,
 f(2b)-f(b)=sqrt(b)*((sqrt(2)-1)-3/(8*b^(5/2)))=0.

Define the fixed algebraic number M=f(b)=f(2b). Continuity handles the
o(1) endpoint extension uniformly. Therefore

 limsup d_N/N^(4/5)<=M,
 limsup (3k+8d_N-7)/N^(4/5)<=8M.

The exact finite theorem proves, for every epsilon>0 and every sufficiently
large integer N,

 C_B2(F_N)>=4N-(8M+epsilon)*N^(4/5).

This is a uniform certificate for the single CURRENT language. The scored
leading pair remains(4,1). To compare the previously prepared margin-rule
envelope, note b<1 because27+18sqrt(2)<54<64. Also (2b)^(5/2)
=4sqrt(2)*3*(sqrt(2)+1)/8=(3/2)*(2+sqrt(2))>2. Thus 2b lies strictly
between the minimum point of f and2, where f is increasing. Hence

 M=f(2b)<f(2)=sqrt(2)+1/8,
 8M<8sqrt(2)+1.

Numerically8M is about12.173705, versus12.313709 from the prepared
margin rule and65/4=16.25 from the earlier32-rule joint family. These
are comparisons of PROVED uniform remainder envelopes for different
languages, not pointwise orderings of their unknown optimal complexities.
No optimality over all constructions or literature priority is asserted.
Only Bertrand's theorem is used for primes; no asymptotic prime number
theorem or conjectural prime-gap hypothesis is invoked.

The finite R84/Q7 predicates and their universal gate transfer are retained
unchanged. They are rerun for final bytes; no new finite computation is a
premise of the guard relaxation or the balanced endpoint calculation.
