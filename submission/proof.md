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
   g=64B^5-27N^3>0 and g^2>648N^6, using exact integers.
   Let p be the smallest prime at least B.
3. Inspect every prime k with p<=k<=2p and k<N. For each put
   m=N-k, r=ceil(m/k), t=ceil(sqrt(r)), and compute d=d_joint
   from the exact integer formula in Section5. A candidate is eligible
   iff t<r, k>r+1, 2k>=t(t+1), and 3<=d<=m. If none is eligible,
   output0. Otherwise select the lexicographically least pair
   (3k+8d-7,k), and use its (k,m,r,t) for every following step.
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

Put c=(27+18sqrt(2))/64 and b=c^(1/5)>0. For B>=2, Bertrand
gives p<2B. For B=1, p=2 satisfies the displayed bound too. Hence

 b*N^(3/5)<p<2*b*N^(3/5)+2,
 p<=k<=2p<4*b*N^(3/5)+4 whenever an eligible degree is selected.

The prime scan for p and exhaustive inspection of all integers through
2p are polynomial in N; trial division is polynomial in these integer
values. B uses integer binary search and exact positive-gap square
comparisons. No floating-point selection or search for a hard truth table
occurs. Each candidate's ceiling square root, threshold and tie-breaking
key use integer arithmetic. Field multiplication uses schoolbook O(k^2)
operations, the powers O(r^2) multiplications, and symmetric recurrence
O(rt). A coarse O((N+1)^6) total bit bound still suffices, including
rejected inputs and the exhaustive polynomial-size prime interval.

Let A be this fixed evaluator, L={x:A(x)=1}, W(N)=0, and let V(x,w)
accept exactly when w is empty and A(x)=1. If x is in L the empty witness
is accepted; if a permitted witness is accepted then A(x)=1 and x is in L.
This proves both directions with one uniform polynomial-time verifier.
Witnesses are not inputs to the target circuits.

## 4. All sufficiently large integer lengths pass the guards

First consider the baseline candidate k=p. The bounds in Section3
give m=N-p=Theta(N), r=Theta(N^(2/5)), t=Theta(N^(1/5)). Thus
p<N, t<r, p>r+1 and 2p>=t(t+1) eventually. The supplied disperser
threshold d_joint is O(N^(4/5)), hence 3<=d_joint<=m eventually.
Thus the eligible set is nonempty at EVERY sufficiently large integer N.

For the selected k, p<=k<=2p gives the same uniform Theta estimates:
k=Theta(N^(3/5)), m=Theta(N), r=Theta(N^(2/5)), t=Theta(N^(1/5)).
All stated eligibility hypotheses hold by the explicit filtering, and
Section15 of disperser.md proves those are sufficient. The joint branch
is eventually below A0 uniformly over this bounded k range, so its
closed eventual dimension and all sublinear estimates apply. No monotonicity
of the selected degree k(N), or sparse-length restriction, is needed.

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
* Erdős, PNAS35(7), July1949, pp374–384, DOI10.1073/pnas.35.7.374,
  Sections1 (equation3) and3, classical prime number theorem. The exact
  standard premise and every corollary used are stated in the final section.
  https://www.renyi.hu/~p_erdos/1949-02.pdf
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

## One fixed family chooses its best explicit finite certificate

The CURRENT language uses c=(27+18sqrt(2))/64 and B^5>c*N^3 to
find the baseline prime p. It then exhaustively inspects every prime
in [p,2p] and uses the eligible degree minimizing (3k+8d_joint-7,k).
All eligibility tests and the objective are the explicit integer formulas
above. This is ONE deterministic uniform evaluator fixed in advance,
not a choice of language depending on an accuracy request, a truth table,
an unproved circuit complexity, or advice. All input and prefix costs
are assessed at the selected k and actual length N.

The complete proof for that selected degree is supplied independently:
field validity, payload nonconstancy, essential-prefix removal and the
unrestricted additive-seven gate theorem yield4N-3k-8d_joint+7. The
minimum is taken over explicit CERTIFICATES, never actual circuit size.
For the baseline p define its analogous scalar certificate remainder R_p.
By Section4 it is eligible eventually, and therefore R_N<=R_p for all
sufficiently large N. This orders the explicit certified remainders of
two parameter-selection rules. It does not order unknown optimal circuit
complexities for their different resulting languages pointwise.

Selected k(N) need not be nondecreasing. Accordingly no infinite-actual-
square-length claim is asserted for this optimizer. Disperser.md states
such first-crossing conclusions only conditionally for monotone wrappers.
The current main theorem and uniform envelope hold at EVERY sufficiently
large integer length, without any such condition.

## Uniform envelope inherited from the eligible baseline prime

For the baseline p let a_N=p/N^(3/5), lying in [b,2b+o(1)]. Its
r_p=N/p+O(1), t_p=sqrt(r_p)+O(1) give

 d_p=p*sqrt(r_p)+r_p^2/2+O(p+r_p),
 d_p/N^(4/5)=f(a_N)+o(1), f(a)=sqrt(a)+1/(2a^2).

The derivative of f has the sign of a^(5/2)-2; hence its maximum on
a compact interval is at an endpoint. Our constant satisfies

 b^(5/2)=3*(sqrt(2)+1)/8,
 f(2b)-f(b)=sqrt(b)*((sqrt(2)-1)-3/(8*b^(5/2)))=0.

Put M=f(b)=f(2b). Uniform continuity handles the o(1) endpoint
extension. Thus limsup R_p/N^(4/5)<=8M. Since R_N<=R_p eventually,

 limsup R_N/N^(4/5)<=8M,
 C_B2(F_N)>=4N-(8M+epsilon)*N^(4/5)

for every epsilon>0 at every sufficiently large integer N. This is a
certificate for the SINGLE current language, with leading pair(4,1).
No argument incorrectly bounds the SELECTED a by2b: selected k may
be as large as2p, so its direct ratio bound is4b. The scalar minimum
comparison is the reason the sharper baseline envelope still applies.

For comparison with preceding uniform envelopes, b<1 since
27+18sqrt(2)<54<64, while (2b)^(5/2)=(3/2)*(2+sqrt(2))>2.
Thus2b lies strictly between the minimum point of f and2, where f
increases. Hence M=f(2b)<f(2)=sqrt(2)+1/8. Numerically8M is about
12.173705, versus12.313709 from the earlier margin-rule envelope and
16.25 from the old32-rule joint family. These are comparisons of proved
uniform certificates for distinct languages, not leading-score gains or
pointwise circuit-complexity orderings. No global optimality or literature
priority is asserted. This preceding envelope uses only Bertrand; the
stronger final envelope below additionally uses the classical prime number theorem.

The finite R84/Q7 predicates and universal gate transfer are unchanged.
A fresh final intake, declared Docker run and independent complete critique
bind all evidence to this optimizer. No finite scan proves a universal
analytic statement; the exhaustive prime scan IS the deterministic family
definition and is polynomial in actual input length.


## The classical prime number theorem sharpens the SAME family

We now use one additional established standard number-theory theorem:

 PNT: pi(x)/(x/log(x)) ->1 as real x->infinity,

where pi(x) counts primes not exceeding x. This is the classical prime
number theorem, not a conjectural error term or unproved gap hypothesis.
A primary proof is Erdős, *On a New Method in Elementary Number Theory
Which Leads to An Elementary Proof of the Prime Number Theorem*,
PNAS35(7), July1949,374–384, DOI10.1073/pnas.35.7.374, Section1
Eq(3) and Section3. The standard theorem is the external premise; all
corollaries and its application to the current fixed evaluator follow here.
No model or runtime result is substituted for that mathematical premise.

For each FIXED delta>0, PNT gives

 pi((1+delta)x)-pi(x)=(delta+o(1))*x/log(x)>0

for every sufficiently large real x. This follows by dividing both PNT
asymptotics by x/log(x) and using log((1+delta)x)/log(x)->1.
Consequently the least prime P(x)>=x has P(x)/x->1: it is at least1,
and for each delta it is at most1+delta eventually. This does not assert
any quantitative rate or efficient subpolynomial prime search.

For our baseline prime p=P(B), B/N^(3/5)->b, so p/N^(3/5)->b.
Let alpha=4^(1/5)=2^(2/5), the global minimizer of f(a). The power
identities in the preceding section show b<alpha<2b. Set
q_N=P(alpha*N^(3/5)). Then q_N/N^(3/5)->alpha, and the strict endpoint
gaps imply p<=q_N<=2p at every sufficiently large integer N. The
uniform eligibility argument of Section4 applies to q_N too: all guards,
2q_N>=t_q(t_q+1) and3<=d_q<=N-q_N eventually hold. Therefore q_N
is one of the ACTUAL inspected eligible candidates, not advice supplied
to an unspecified evaluator.

The exact objective minimum gives R_N<=R_(q_N), and the uniform dimension
expansion yields

 limsup R_N/N^(4/5)<=8*f(alpha)=10*2^(1/5).

For the matching scalar lower bound, the SELECTED a_N=k/N^(3/5) stays
in a fixed compact positive interval by p<=k<=2p and the earlier Bertrand
bounds. Its uniformly valid dimension expansion gives

 R_N/N^(4/5)=8*f(a_N)+o(1)>=8*f(alpha)+o(1),

because f has global minimum at alpha. Thus the PROVED CERTIFICATE
REMAINDER, for the same exact optimizer and unchanged payload, satisfies

 R_N/N^(4/5) ->10*2^(1/5).

The constant is about11.486984, strictly smaller than8M because b and2b
are distinct from the unique minimizer alpha. Combining with the finite
gate theorem, for every epsilon>0 and every sufficiently large integer N,

 C_B2(F_N)>=4N-(10*2^(1/5)+epsilon)*N^(4/5).

This is a lower-bound certificate for ONE fixed total uniform NP family
at all actual lengths. It is NOT an asymptotic equality for C_B2(F_N),
not a leading coefficient/exponent gain, and not an optimality theorem
among all constructions. The code is byte-identical to the preceding
optimizer; only its proved analytic guarantee is strengthened. No sparse
subsequence replaces the all-integer-length quantifier.
