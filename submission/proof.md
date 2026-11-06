# Eventual coefficient 4 with a quantified N^(4/5) remainder

## Result and provenance

We prove that the single total language F defined below has fixed real
K_R>0 and integer N_0 such that, for every integer N>=N_0 at its actual
input length,

    C_B2(F_N) >= 4N - K_R*N^(4/5).

More precisely, the proof obtains C_B2(F_N)>=4N-R(N) with
R(N)=O(N^(4/5)) and then selects K_R and N_0 from that displayed remainder.
Its complete finite projection is (4,1), while the eventual Big-O declaration
is stronger than a bare 4N-o(N) liminf statement. This is an application of
published mathematics, not a claim to have discovered a new circuit lower
bound. The gate-elimination argument is Mishra--Santhanam, ECCC TR26-260,
Revision 1, November 3, 2026, Sections 2--4 and Theorem 4.1. Their newer
4n-10d theorem is deliberately distinguished from the original report's
3.5n-7d statement; only Revision 1 is used for the coefficient 4.

The hard function is a certificate-prefix wrapper around the descending-power
Ben-Sasson--Kopparty (BK) symmetric-polynomial construction. The complete
algebraic disperser argument is supplied in `disperser.md`, and the complete
gate argument in `elimination.md`. The prefix wrapper makes the field
polynomial part of the input rather than nonuniform advice and charges all of
its bits. The whole function is not asserted to be an affine disperser; only
an explicitly fixed valid-prefix restriction is.

## 1. Definition at every input length

Bits are indexed in their written order, and the first bit of each field block
is its constant coefficient. On an input x of length N, the fixed evaluator A
performs the following steps.

1. If N=0, return 0. Otherwise compute, by integer binary search, the least
   positive integer B with B^5>32N^3. Let k be the least prime at least B,
   found by incrementing and trial division.
2. If k>=N, return 0. Otherwise set m=N-k, r=ceil(m/k), and
   t=ceil(sqrt(r)). If t>=r or kt<r(r+1), return 0.
3. Read the first k bits b_0,...,b_(k-1) and form the monic polynomial

       Q(Z)=Z^k+sum_(j=0)^(k-1)b_j Z^j over F_2.

   Test both Z^(2^k)=Z (mod Q) and gcd(Q,Z^2+Z)=1. Return 0 if either
   fails. These computations use repeated squaring and polynomial Euclid.
4. Let K=F_2[Z]/(Q), represented in the basis 1,Z,...,Z^(k-1). Append
   rk-m zero bits to the remaining m input bits, split the result into r
   blocks of k bits, and let a_1,...,a_r in K be their images. This is an
   injective F_2-linear embedding of the payload.
5. For 1<=i<=r put w_i=r+2-i and

       y_i=a_i^(2^(w_i)-1).

   Return the constant coefficient of the elementary symmetric polynomial
   of order t in y_1,...,y_r, taking the representative of degree below k.

Evaluate the symmetric polynomial without division in characteristic two.
Initialize E_0=1 and E_1=...=E_t=0. After reading y_i, update

    E_j <- E_j+y_i E_(j-1)

in decreasing order j=min(i,t),...,1. Splitting t-subsets according to
whether they contain i proves the recurrence; the answer is the constant
coefficient of E_t. Binary exponentiation computes each y_i. Consequently
`family.py` is only an executable restatement of this definition, not an extra
source of choices.

Every length and every bit string is covered, including invalid prefixes and
small or otherwise guarded lengths. Q is read from x, never supplied as
per-length advice.

## 2. The prefix test is exact and a valid prefix exists

Work in an algebraic closure of F_2. The polynomial Z^(2^k)-Z is squarefree,
since its derivative is 1. Its roots are closed under addition,
multiplication, and nonzero inversion, so they form a field of size 2^k.
Frobenius z->z^2 permutes this field. The product of the linear factors in a
Frobenius orbit is fixed by Frobenius and hence lies in F_2[Z]. A proper
factor over F_2 would select a proper invariant subset of the orbit, so this
orbit polynomial is irreducible. Conversely, the roots of an irreducible
factor form one orbit.

Every irreducible factor degree therefore divides the orbit length k. For
prime k, the only possibilities are 1 and k. The fixed points are exactly 0
and 1, while there are 2^k>2 roots. Hence an irreducible degree-k polynomial
exists for every prime k.

The first test says that Q divides the squarefree polynomial Z^(2^k)-Z.
The second excludes its factors Z and Z+1. All remaining irreducible factors
have degree k. As Q itself has degree k, these tests hold exactly when Q is
irreducible. In particular, some input prefix passes at each length for which
the other guards pass; no algorithm for finding that prefix is needed by the
uniform evaluator or the lower-bound argument.

## 3. One polynomial-time evaluator and empty-witness NP verifier

By definition of B, B^5>32N^3, hence B>2N^(3/5). The previous integer has
fifth power at most 32N^3, giving B<=2N^(3/5)+1. Bertrand's postulate gives
k<2B. Therefore

    2N^(3/5)<k<4N^(3/5)+2.

The prime search and trial division are polynomial in the input length. So
are binary search, polynomial Euclid, repeated field squaring, schoolbook
field multiplication, and the O(rt) symmetric-polynomial recurrence. The
exponents have bit length O(r), and a coarse polynomial bound covering all
accepted and rejected inputs is immediate; the exact exponent is irrelevant.

Let L={x:A(x)=1}. Define the witness bound W(N)=0 and the fixed verifier
V(x,w)=1 exactly when w is empty and A(x)=1. If x is in L, its empty witness
is accepted. Conversely, every accepted permitted witness certifies A(x)=1,
so x is in L. This proves both directions and L in P, hence NP. Witness bits
are not circuit inputs, and circuits remain unrestricted nonuniform full-B2
circuits.

## 4. All sufficiently large lengths pass the guards

The estimates for k imply

    m=N-k=Theta(N),   r=Theta(N^(2/5)),   t=Theta(N^(1/5)).

Consequently k<N and t<r eventually. For the last guard, use
r<N/k+1<(1/2)N^(2/5)+1 and t>=sqrt(r). Then

    r(r+1)/t <= sqrt(r)(r+1)
              <= sqrt((1/2)N^(2/5)+1)*((1/2)N^(2/5)+2).

Dividing this expression by N^(3/5) gives limit 1/(2*sqrt(2)), whereas
k/N^(3/5)>2. Thus kt>=r(r+1) eventually. These limits range over every
sufficiently large integer N, not over a sparse subsequence.

## 5. Hard prefix restriction and exact input-length accounting

For each sufficiently large N, fix any irreducible degree-k polynomial Q and
fix the first k circuit inputs to its low coefficients. Such Q exists by
Section 2. The restricted function G_Q has m=N-k actual input bits. Define

    d_N=floor(kr/t)+kt+r(r+1)+1.

For prime k, 1<=t<r, and kt>=r(r+1), `disperser.md` proves that

    (a_i) -> pi(sum_(|I|=t) product_(i in I) a_i^(2^(r+2-i)-1))

is nonconstant on every affine subspace of K^r whose dimension is strictly
greater than kr/t+kt+r(r+1), for every nonzero F_2-linear functional pi.
Our pi is the constant coefficient, which is nonzero because pi(1)=1. The
zero-padded payload embedding is injective and linear, so G_Q is an affine
disperser for dimension d_N on its m input bits. Padding is part of this
embedding and adds no circuit inputs.

The growth estimates above give

    kr=O(N),   kt=O(N^(4/5)),   r(r+1)=O(N^(4/5)),

and also kr/t=O(N^(4/5)); hence d_N=O(N^(4/5)). We use this explicit BK
dimension rather than assuming any later logarithmic-dimension construction.

Fixing inputs of any full-B2 circuit for F_N yields a full-B2 circuit for G_Q
with no additional gates. `elimination.md` proves the size bound

    |C| >= 4m-10d

for every affine disperser of dimension d on its m input bits, including
arbitrary sharing, repeated pins, constants, degenerate tables, and terminal
outputs. Applying it to the restricted circuit, or its no-larger normalized
form, yields

    C_B2(F_N) >= 4(N-k)-10d_N = 4N-R(N),
    R(N)=4k+10d_N=O(N^(4/5)).

The lower bound fixes Q only to exhibit a restriction that every circuit for
F_N must compute. It is not advice to A, which reads and tests every prefix.

The displayed relation R(N)=O(N^(4/5)) is a proved property of the explicit
length parameters, not an inference from a liminf or from the numerical pair.
Consequently there exist a fixed real K_R>0 and an integer N_0, at least the
eventual guard threshold, such that

    0 <= R(N) <= K_R*N^(4/5) for every integer N>=N_0.

Both K_R and N_0 depend only on this fixed theorem. Therefore

    C_B2(F_N) >= 4N-K_R*N^(4/5) for every integer N>=N_0,

which is exactly the complete v2 eventual_lower_bound guarantee with leading
4N and Big-O remainder exponent 4/5.

For any epsilon>0, choose N_epsilon beyond the eventual guard threshold and
large enough that R(N)/N<=epsilon for all N>=N_epsilon; this is possible
because N^(4/5)/N tends to zero. Then C_B2(F_N)>=(4-epsilon)N for every such
integer N. Together with the displayed eventual bound above, this establishes
both the complete declaration and its joint projection (4,1) for this one
language and one theorem.

## 6. Evidence boundary and source identity

`elimination.md` reproduces the required affine-system invariant, all four
parameter operations, all deletions, all five cases, the terminal-output
argument, termination, and the final potential summation. `disperser.md`
reproduces the coefficient-isolation, block-decomposition, and block-mixing
arguments. `audit.py` checks only bounded Boolean, rational, prefix, and
executable-family predicates; it is not a substitute for either universal
proof. No checker result, earlier judgment, leaderboard row, or retrieved
review is a mathematical premise.

Sources used for the stated theorem:

* Om Swostik Mishra and Rahul Santhanam, *An Explicit 4n-O(log(n)) Size Lower
  Bound against General Boolean Circuits*, ECCC TR26-260, **Revision 1**,
  2026-11-03, Sections 2--4, Lemma 2.2, Lemma 4.2, and Theorem 4.1.
  https://eccc.weizmann.ac.il/report/2026/260/revision/1/download/
  SHA256 a123943ed05c9569385d3b9d98c695630a110559f33a4c990ec44e2a22b4d921.
  The original 2026-10-19 PDF and its 3.5 result are not used for this bound.
* Eli Ben-Sasson and Swastik Kopparty, *Affine Dispersers from Subspace
  Polynomials*, author version 2010-03-12, Lemmas 3.7, 3.10, 3.13, Theorem
  6.1, Lemma 7.2, Theorem 7.3, and Claim 7.4. The needed arguments are
  specialized and supplied in `disperser.md`.
  https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
  SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf.
* Jiatu Li and Tianqi Yang, ECCC TR21-023, 2021-02-20, Sections 1.1, 2.1,
  2.4 (including footnotes 3--4), and 3.3, only for the protected full-B2
  normalization conventions. Their 3.1 lower-bound case analysis is not used.
  https://eccc.weizmann.ac.il/report/2021/023/download/
  SHA256 4530cbaf38f0f9a09704e94f7d1f57109895b71b92d0752a285e07d1109a4062.
