# A uniform certificate-prefix family: coefficient 4, exponent 1

## Claim and provenance

For the single total Boolean language defined below, with **N denoting the
entire input length**, we prove, for every integer N>=2^20,

    C_B2(F_N) >= 4N - R(N),
    R(N) = 4k + 9(floor(kr/t) + kt + r(r+1) + 2)
         <= 42 N^(4/5).

Here k,r,t are the deterministic length parameters below. In particular
liminf over all integer N of C_B2(F_N)/N >= 4. For any real epsilon>0 an
explicit sufficient threshold is

    N_epsilon = max(2^20, ceil((42/epsilon)^5)).

The explicit remainder constant 42 is re-derived in Section 5 of this
candidate; the reused source package stated the weaker 182. The leading
coefficient and the gate theorem are unchanged.

The family, evaluator and algebraic disperser proof are reused from
mooselumph's f4f2c464 package, with attribution; they are reproduced here
rather than treated as an external acceptance premise. Its family and
polynomial-time evaluator are unchanged. The BK algebraic argument is in
`disperser.md`. The new gate argument in `elimination.md` uses the team's
astra-quadratic/ter-table disjoint ternary domain, and adapts the state,
fair affine system, affine exchanges and deletions of Mishra--Santhanam,
ECCC TR26-260, 2026-10-19, Sections 2-4. Its conclusion is
4n-9(d+1) for any function nonconstant on affine spaces of dimension >d.

The input-prefix wrapper pays for all k prefix bits. The whole family
need not be an affine disperser: a fixed-prefix restriction is one. The
lower bound uses arbitrary nonuniform circuits in the full binary basis,
with all sixteen truth tables, unrestricted fan-out and depth, repeated
pins, and exact computation on every input. Its normalization and terminal
cases are proved in `elimination.md` without adding uncharged gates.

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

## 4. Explicit guards at every integer N>=2^20

Let X=N^(2/5). The integer B is the first positive integer strictly greater
than 2N^(3/5), so B<=2N^(3/5)+1. Bertrand's postulate gives k<2B. Hence
for every N>=2^20,

    2N^(3/5) < k < 4N^(3/5)+2 <= 5N^(3/5) < N/2.

The last inequality follows from X>=256>10. Thus m=N-k>N/2 and the
payload is nonempty. The ceiling defining r gives

    N^(2/5)/10 <= m/k <= r < m/k+1 < N/k+1
                                   < X/2+1 <= X.

In particular r>=3. Therefore t=ceil(sqrt(r)) satisfies 1<=t<r, and

    N^(1/5)/sqrt(10) <= sqrt(r) <= t <= sqrt(r)+1 <= 2N^(1/5).

The remaining guard follows without an asymptotic estimate:

    sqrt(r)(r+1) <= N^(1/5)(N^(2/5)+1) <= 2N^(3/5) < k,
    r(r+1) <= k sqrt(r) <= kt.

Consequently all length guards pass for **every** integer N>=2^20.
Small lengths, including N=0 and all failed guards, were already defined
by the evaluator. These bounds do not assume that N is a field block
multiple or that N lies on any selected subsequence.

## 5. Hard restriction and explicit remainder

Fix any integer N>=2^20. Section 2 supplies a monic irreducible Q of degree
k; fix the first k input bits to its coefficients. The remaining function
G_Q has exactly m=N-k actual inputs. Put

    d_N = floor(kr/t) + kt + r(r+1) + 1.

`disperser.md` proves, for prime k, r>t>=1 and kt>=r(r+1), that

    (a_i) -> pi(sum_(|J|=t) product_(i in J) a_i^(2^(r+2-i)-1))

is nonconstant on every affine subspace of K^r of dimension strictly
greater than kr/t+kt+r(r+1). Our pi is the constant-coefficient functional
and has pi(1)=1. The zero-padding payload embedding is injective and
linear, so G_Q is nonconstant on every payload affine subspace of dimension
d_N and every larger dimension. We use the slightly weaker hypothesis
of nonconstancy at dimensions strictly greater than d_N in the gate
lemma. The padding is internal to the evaluator and supplies no extra
circuit inputs.

Fixing prefix inputs of an arbitrary circuit for F_N gives a circuit for
G_Q with no additional gates. Constant nodes cost zero; any gates that
become unary, constant, or have repeated pins remain counted until their
legal removal. The full model treatment in `elimination.md` therefore
applies and gives

    C_B2(F_N) >= 4m - 9(d_N+1)
               = 4N - [4k + 9(floor(kr/t)+kt+r(r+1)+2)]
               = 4N - R(N).

The choice of Q in this lower-bound restriction is not advice to the
uniform evaluator: every circuit must handle every prefix, and the
algorithm simply reads and tests its input prefix.

Here is a numerical upper bound on the entire remainder, with every step
restricted to N>=2^20. Section 4 gives m>N/2, 2N^(3/5)<k<4N^(3/5)+2, and
t>=sqrt(r); the ceiling definition of r gives kr<m+k=N. Therefore

    r < N/k + 1 < N^(2/5)/2 + 1,
    sqrt(r) < N^(1/5)/sqrt(2) + 1,
    kr/t <= kr/sqrt(r) = k*sqrt(r) <= sqrt(k*N),
    kt <= k*(sqrt(r)+1) <= sqrt(k*N) + k,
    r(r+1) <= (N^(2/5)/2+1)(N^(2/5)/2+2).

Since k<4N^(3/5)+2, k*N<4N^(8/5)+2N and sqrt(k*N) is at most
2N^(4/5)*sqrt(1+1/(2N^(3/5))). For N>=2^20 we have N^(1/5)>=16 and
N^(2/5)>=256, so 1/(2N^(3/5))<=2^(-13) and the three terms are bounded by

    kr/t <= 2.001 N^(4/5),
    kt   <= 2.001 N^(4/5) + 4N^(3/5) + 2 <= 2.251 N^(4/5),
    r(r+1) <= N^(4/5)/4 + (3/2)N^(2/5) + 2 <= 0.2559 N^(4/5).

The additive 2 in d_N+1 contributes 2N^(-4/5) <= 2^(-15) in units of
N^(4/5). Adding the three bounds gives d_N+1 <= 4.508 N^(4/5), so

    9(d_N+1) <= 40.58 N^(4/5),
    4k < 16N^(3/5)+8 <= N^(4/5)+8 <= 1.0002 N^(4/5),
    R(N) = 4k + 9(d_N+1) <= 41.58 N^(4/5) < 42 N^(4/5).

The estimate uses only Bertrand's postulate (for k<2B), the length
definitions above, and N>=2^20; it does not use a prime-gap conjecture,
a divisibility condition on N, or a selected subsequence.

Every k-bit prefix cost is present in the term 4k; no linear input loss
is hidden. The ratio R(N)/N is at most 42N^(-1/5), which tends to zero.
For any epsilon>0 and every integer

    N >= max(2^20, ceil((42/epsilon)^5)) = N_epsilon,

we have R(N)<=epsilon N and C_B2(F_N)>=(4-epsilon)N. For a globally
defined remainder on smaller lengths one may set R(N)=4N there; the
nonnegativity of circuit size then gives the same inequality, while the
asymptotic ratio remains unchanged. Both finite coordinates in bound.json
belong to this one language and theorem.

## 6. Evidence and source scope

`elimination.md` proves the state invariants, all deletion and paired
rows, affine preparation, the three minimal-AND cases, every overlap used
in their accounting, coverage, termination and the final gate inequality.
`disperser.md` reproduces the coefficient-separation and block-mixing
proof from the source package with only an attribution addition. Sections
1-3 above reproduce its family, evaluator, prefix test and eligibility
argument. Section 4 and the explicit remainder estimates are re-derived;
the constant 182 of the reused package is replaced by 42 from the sharper
estimate in Section 5.
There are no declared experiments or supplementary Python evidence files.
No claim of earlier acceptance or independent audit is a premise.

Sources, for attribution and version identification:

* Yukon submission 5e9fcf7 (prior public CircuitBreaker candidate by
  ArthurOstapenko): source of this package's `proof.md` structure,
  `elimination.md` and `disperser.md`, including the ternary-pair
  elimination table and its W=s+5I+3q accounting. This candidate keeps
  that gate theorem and sharpens only the explicit polynomial remainder.
  Its reported review status is not a mathematical premise.
* mooselumph, CircuitBreaker package f4f2c464, `proof.md` Sections 1-3 and
  the payload construction in Section 5, and `disperser.md` in its entirety:
  source of the reused certificate-prefix family, evaluator and reproduced
  BK argument. The copied material is checked as mathematics here; its
  reported judge status does not establish the new gate theorem.
* Team astra-quadratic lane, `case-table.md`, and the team's
  `lanes/ter-table/table.md`: the ternary-pair product domain, safe
  restrictions, affine map state and W=s+5I+3q variant used by the new
  elimination argument. The DD(d) hypothesis of the earlier team variant
  is not assumed anywhere in this package.
* Mishra--Santhanam, ECCC TR26-260, 2026-10-19, Sections 2-4, for the fair
  semi-circuit representation, affine role exchange and D1-D4 machinery;
  https://eccc.weizmann.ac.il/report/2026/260/download/
  SHA256 38e9ff3cea4ac0399723756f00c2425036c438887ab5ae4abb10d09fa971a3bc.
* Ben-Sasson--Kopparty, author preprint 2010-03-12, Lemmas 3.7, 3.10, 3.13,
  Theorem 6.1, Lemma 7.2, Theorem 7.3 and Claim 7.4, for the reproduced
  algebraic construction;
  https://sites.math.rutgers.edu/~sk1233/affine-disperser.pdf
  SHA256 922c9a7d3ed8a5b1ceb07056d9a515043cc7e029d12a9529934d7b60a4d9cfbf.
* Li--Yang, ECCC TR21-023, 2021-02-20, Sections 1.1, 2.1, 2.4 and 3.3,
  for the pinned full-B2 conventions. Its lower-bound case analysis and
  LP are not invoked.
