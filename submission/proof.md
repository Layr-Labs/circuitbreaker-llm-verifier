# A uniform certificate-prefix family: coefficient 4, exponent 1

## Claim and provenance

For the single total Boolean language defined below, with **N denoting the
entire input length**, we prove

    C_B2(F_N) >= 4N - R(N),       R(N)/N -> 0.

The joint pair is (4,1). The leading coefficient is inherited mathematics.
Section 4 of elimination.md supplies a finite one-gate refinement,
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
2. Let B be the smallest positive integer with B^5 > 32*N^3. Let k be the
   smallest prime at least B. Use deterministic trial division to find k.
3. If k>=N, output 0. Otherwise put m=N-k, r=ceil(m/k), t=ceil(sqrt(r)).
   If t>=r or k*t<r*(r+1), output 0.
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

Bertrand's prime theorem gives k<2B, and so

    2*N^(3/5) < k < 4*N^(3/5)+2.

The prime scan and trial division are polynomial in N. B is computed by
integer binary search and fifth-power comparisons. All subsequent integer,
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
Consequently k<N and t<r eventually. Also r<N/k+1<(1/2)N^(2/5)+1 and

    r(r+1)/t <= sqrt(r)(r+1)
              <= sqrt((1/2)N^(2/5)+1)*((1/2)N^(2/5)+2).

Dividing the last expression by N^(3/5) gives limit 1/(2*sqrt(2)),
whereas k/N^(3/5)>2. Thus k*t>=r(r+1) eventually. These estimates hold
for every sufficiently large integer N, not just selected lengths.

## 5. Hard restriction and its precise loss

For each sufficiently large N fix any irreducible Q of degree k and fix
the first k bits of x to its coefficients. The remaining function G_Q
has m inputs. Put

    d_N = floor(kr/t) + kt + r(r+1) + 1.

`disperser.md` proves that, for prime field degree k, r>t>=1 and
k*t>=r(r+1), the function

    (a_i) -> pi(sum_(|I|=t) product_(i in I) a_i^(2^(r+2-i)-1))

is nonconstant on every affine subspace of K^r of dimension greater than
kr/t+kt+r(r+1), for any nonzero F_2-linear pi. Our pi is the constant
coefficient, with pi(1)=1. The payload embedding is injective and linear,
so G_Q is an affine disperser of dimension d_N on its m actual input bits.
Zero padding preserves dimension and does not add circuit inputs.

Since kr=O(N), r=Theta(N^(2/5)) and t=Theta(N^(1/5)),

    d_N=O(N^(4/5)),          k=O(N^(3/5)).

We use this conservative sublinear disperser dimension, not later
logarithmic-dimension constructions. Either suffices for the leading pair.

Fixing prefix inputs of any circuit for F_N yields a circuit for G_Q with
no additional gates. Sections 1–3 of `elimination.md` prove the supplied full-B2 bound
4m-8d+4. Sections 4–5 refine this, for d>=2, to

    |C| >= 4m-8d+6.

Section 5 includes the exact terminal-cover certificate verified by the
declared audit; its universal reduction and normalization are supplied there.

for a dimension-d affine disperser. It includes arbitrary sharing,
repeated pins, constants, degenerate gates and terminal-output cases.
Here d_N>=2 by its displayed formula and r>t>=1. Apply it to the
restricted circuit (or its no-larger normalized form):

    C_B2(F_N) >= 4(N-k)-8d_N+6
               = 4N-(4k+8d_N-6).

Set R(N)=4k+8d_N-6 eventually and extend it on the finitely many guarded
lengths. Then R(N)=O(N^(4/5))=o(N). The choice of Q in this lower-bound
restriction is not advice to the evaluator: every circuit must work for
all prefixes, while the uniform algorithm simply reads and tests its input.

For each epsilon>0, take N_epsilon large enough that the guards pass and
R(N)/N<=epsilon for every N>=N_epsilon. This gives exactly
C_B2(F_N)>=(4-epsilon)N for all those integers. Both coordinates in
bound.json belong to this one theorem and one language.

## 6. Evidence and limits

`disperser.md` contains the coefficient separation and block-mixing proof.
`elimination.md` contains the complete gate argument. `audit.py` supplies
supplementary finite checks of Boolean identities, rational estimates and
the prefix test; these do not replace either universal proof. No claimed
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

## Package provenance

The prefix construction, algebraic proof and finite audit are reused from the repository candidate-3_5 fixture, originally submission db5f9c6e-607f-424a-89f8-5a49d29d07df. The complete gate proof is Sections 1–3 of submission 85e46e2, source commit 69361e0f153025db560934b9159c53d69cd895f1. It establishes 4m-8d+4, which suffices with the sublinear BK dimension. No prior judgment is a premise. The older potential arithmetic audit is supplementary arithmetic only and does not certify the new gate argument.


## Additional finite diagnostics in this package

The historical diagnostics below originated in submission
2e89196b-709d-4eaa-97e3-1fbb061b5644, source commit
7557a2e956d84db6e2e09efbdc6e35d0e32fdeaf, and are carried here through
the immediate aab75a1e source identified below. The same family and leading
guarantee are retained; the present finite theorem uses the complete gate source.
The additional audit exhausts all invertible binary auxiliary matrices of
orders one through three, every nonzero parameter column, every legal pivot
coordinate, every right-hand side and both values of the new parameter.
It independently solves both systems and checks agreement under the role
exchange. It also checks all sixteen Boolean tables on the three-point pair
domain, all eight AND normalizations, forcing restrictions and seven complete
potential inequalities. These diagnostics can detect transcription or local
algebra errors; they do not prove sparsity, exhaustiveness, deletion counts,
termination or any asymptotic statement. Those obligations are discharged by
the written argument, for arbitrary dimensions and unrestricted sharing.
The historical quadratic-pivot and 7/2-potential tests retained in audit.py
are explicitly supplementary and are not the coefficient-four proof.

## Consolidated supplied finite theorem

This package directly reuses the simple seven-file family and audit from
aab75a1e-b0fb-4f16-bce8-1c3c1279bf73, immutable commit
824745a80accf0e71c45ef6a92b9f923a250613e. That package extracted only
Section 1 of the gate source and used 4m-10d. Here elimination.md retains
all three gate sections from 85e46e26-d8d2-42c7-8a30-5c09ac427208,
immutable commit 69361e0f153025db560934b9159c53d69cd895f1. Its four-object
constant-closure count, potential s+2q+4I, and permanent-inert-coordinate
ledger establish the stronger finite bound 4m-8d+4. These are supplied
source results, not new discoveries or accepted mathematical premises.
Their complete proofs are included and must stand on their own.

The construction and asymptotic leading expression remain those of the
existing simple certificate-prefix family. Applying the supplied finite
theorem gives 4N-4k-8d_N+4, with the same sublinear k and d_N. This
consolidation improves the original checkout's 31/10 guarantee, and the
direct seven-file predecessor's finite remainder, but does not exceed
previous provisional coefficient-four expressions or the complete gate
source. Supplementary potential and terminal-ledger arithmetic in audit.py
certifies only its finite predicates, not the universal closure or inertness
argument. No source's review verdict or score is used as a premise.

## New one-gate finite refinement

Section 4 of elimination.md strengthens the complete supplied generic
theorem by one gate when d>=2: 4m-8d+5. It proves three deletions in
(E1), excludes the tight final(E2) state by a nine-point pair-sum
pigeonhole argument in F_2^5, and checks the remaining terminal pair
counts. Since d_N>=2 eventually, the main finite theorem for this same
family is therefore 4N-4k-8d_N+5. The earlier supplied +4 theorem above
is retained as a supporting source result; the new theorem supersedes
its finite remainder. No improvement over the leading coefficient four
or the existing expression score is claimed.

## Exact finite-cover refinement

Section 5 of elimination.md strengthens the finite theorem again to
4m-8d+6 for d>=2. It reduces a hypothetical +5 equality to a seven-gate
read-once tree on eight active coordinates with three promised pairs,
and uses the analytic orientation split and the exact cover embedded in
audit.py. The declared audit verifies all 19945440 remaining normalized
shape/assignment combinations using exact bitwise masks; normalization and full-B2
coverage, arbitrary sharing before the terminal reduction, inertness,
small dimensions and all quantifiers are proved in the Markdown. This
updates the same family's eventual bound to 4N-4k-8d_N+6. The preceding
+4 and +5 statements remain supporting weaker results. The leading
guarantee and expression score are still (4,1).

Section 5.5 also proves that the uniform additive constant six is sharp
for an inequality valid at every finite n and integer d>=2: xy XOR z
on three inputs has circuit size two and no constant affine plane.
This separate unscored sharpness claim concerns only that universal
additive constant, not the asymptotic coefficient or other theorem forms.

## Analytic orientation split and smaller certificate

The new first part of elimination.md Section 5.2 proves a stronger local
statement: when either paired input of G has forcing value zero, ANY
function of (G,hx,hy,p,hp,u,z), with no other dependence on x,y, is
constant on a legal affine three-flat. The proof fixes two legal paired
lines, obtains one unused direction, and uses 66>31 pair differences
among twelve same-colored points in a 24-point domain in F_2^5. There
is no circuit-size or read-once hypothesis for this orientation lemma.

The remaining read-once case has both forcing values one. Only four
paired leaf signs remain free, reducing its exact enumeration from
64*3^6 to 16*3^6 assignments per representative. The 16480 retained
plane witnesses are embedded as untrusted data in audit.py, which
checks their actual coverage. No separate certificate file is needed.
This strengthens the local supporting statement and reduces computation
without changing the evaluator, finite main bound, additive sharpness,
all-length quantifiers or leading expression (4,1). The immediate
material predecessor is ecba81b9-ec6a-4743-8179-10b5f3c9c197; its full
seven-gate reduction and prior finite +6 result are retained. No prior
review supplies mathematical authority.
